const std = @import("std");

pub const TransportAxis = enum {
    column,
    row,
    vertical,
};

/// Site-file NCNG semantics for horizontal exchange.
pub const HorizontalExchange = enum {
    enabled,
    standalone,
};

pub const SaltEquilibriumMode = enum {
    static_concentrations,
    dynamic_transport,
};

/// Free ions and inorganic carbon species in exact REDIST order.
pub const PrimaryIonFlux = struct {
    aluminum_mol_per_step: f64 = 0,
    iron_mol_per_step: f64 = 0,
    hydrogen_mol_per_step: f64 = 0,
    calcium_mol_per_step: f64 = 0,
    magnesium_mol_per_step: f64 = 0,
    sodium_mol_per_step: f64 = 0,
    potassium_mol_per_step: f64 = 0,
    hydroxide_mol_per_step: f64 = 0,
    sulfate_mol_per_step: f64 = 0,
    chloride_mol_per_step: f64 = 0,
    carbonate_mol_per_step: f64 = 0,
    bicarbonate_mol_per_step: f64 = 0,
};

/// Dissolved metal hydrolysis and ion-pair species in source order.
pub const MetalComplexFlux = struct {
    aluminum_monohydroxide_mol_per_step: f64 = 0,
    aluminum_dihydroxide_mol_per_step: f64 = 0,
    aluminum_trihydroxide_mol_per_step: f64 = 0,
    aluminum_tetrahydroxide_mol_per_step: f64 = 0,
    aluminum_sulfate_mol_per_step: f64 = 0,
    iron_monohydroxide_mol_per_step: f64 = 0,
    iron_dihydroxide_mol_per_step: f64 = 0,
    iron_trihydroxide_mol_per_step: f64 = 0,
    iron_tetrahydroxide_mol_per_step: f64 = 0,
    iron_sulfate_mol_per_step: f64 = 0,
    calcium_hydroxide_mol_per_step: f64 = 0,
    calcium_carbonate_mol_per_step: f64 = 0,
    calcium_bicarbonate_mol_per_step: f64 = 0,
    calcium_sulfate_mol_per_step: f64 = 0,
    magnesium_hydroxide_mol_per_step: f64 = 0,
    magnesium_carbonate_mol_per_step: f64 = 0,
    magnesium_bicarbonate_mol_per_step: f64 = 0,
    magnesium_sulfate_mol_per_step: f64 = 0,
    sodium_carbonate_mol_per_step: f64 = 0,
    sodium_sulfate_mol_per_step: f64 = 0,
    potassium_sulfate_mol_per_step: f64 = 0,
};

pub const PhosphorusComplexFlux = struct {
    phosphate_mol_per_step: f64 = 0,
    phosphoric_acid_mol_per_step: f64 = 0,
    iron_hydrogen_phosphate_mol_per_step: f64 = 0,
    iron_dihydrogen_phosphate_mol_per_step: f64 = 0,
    calcium_phosphate_mol_per_step: f64 = 0,
    calcium_hydrogen_phosphate_mol_per_step: f64 = 0,
    calcium_dihydrogen_phosphate_mol_per_step: f64 = 0,
    magnesium_hydrogen_phosphate_mol_per_step: f64 = 0,
};

/// The micropore branch has one additional soluble H-silicate pool.
pub const MicroporeSaltFlux = struct {
    primary_ions: PrimaryIonFlux = .{},
    metal_complexes: MetalComplexFlux = .{},
    hydrogen_silicate_mol_per_step: f64 = 0,
    nonband_phosphorus: PhosphorusComplexFlux = .{},
    band_phosphorus: PhosphorusComplexFlux = .{},
};

pub const MacroporeSaltFlux = struct {
    primary_ions: PrimaryIonFlux = .{},
    metal_complexes: MetalComplexFlux = .{},
    nonband_phosphorus: PhosphorusComplexFlux = .{},
    band_phosphorus: PhosphorusComplexFlux = .{},
};

pub const SaltFlux = struct {
    micropore: MicroporeSaltFlux = .{},
    macropore: MacroporeSaltFlux = .{},
};

pub const Inputs = struct {
    transport_axis: TransportAxis,
    horizontal_exchange: HorizontalExchange,
    salt_equilibrium_mode: SaltEquilibriumMode,
    current_layer_thickness_m: f64,
    layer_activity_threshold_m: f64,
    current_cell_flux: SaltFlux,
    positive_neighbor_flux: SaltFlux,
};

pub const State = struct {
    net_flux: SaltFlux,
};

/// Aggregates dynamic salt and phosphorus transfer between adjacent cells.
///
/// Traceability: REDIST.F lines 3566--3765 under enclosing gates 3350 and
/// 3363. The source's `ISALTG != 0` branch is the explicit dynamic mode.
/// Micropores retain 50 source-order molar pools: 12 primary ions, 21 metal
/// complexes, H-silicate, eight non-band phosphorus complexes, and eight band
/// complexes. Macropores retain the same order without H-silicate (49 pools).
/// Net transfer is current-cell face flux minus the already-selected positive
/// neighbor `N6` flux. Candidate state commits atomically after finite
/// source-ordered evaluation.
pub fn aggregate(inputs: Inputs, state: *State) !void {
    if (inputs.transport_axis != .vertical and
        inputs.horizontal_exchange == .standalone)
    {
        return;
    }
    if (inputs.salt_equilibrium_mode == .static_concentrations) return;
    try validateInputs(inputs, state.*);
    if (inputs.current_layer_thickness_m <= inputs.layer_activity_threshold_m)
        return;

    var candidate = state.net_flux;
    try addMicroporeDifference(
        &candidate.micropore,
        inputs.current_cell_flux.micropore,
        inputs.positive_neighbor_flux.micropore,
    );
    try addMacroporeDifference(
        &candidate.macropore,
        inputs.current_cell_flux.macropore,
        inputs.positive_neighbor_flux.macropore,
    );
    state.net_flux = candidate;
}

fn addMicroporeDifference(
    candidate: *MicroporeSaltFlux,
    current: MicroporeSaltFlux,
    positive: MicroporeSaltFlux,
) !void {
    try addDifference(
        &candidate.primary_ions,
        current.primary_ions,
        positive.primary_ions,
    );
    try addDifference(
        &candidate.metal_complexes,
        current.metal_complexes,
        positive.metal_complexes,
    );
    candidate.hydrogen_silicate_mol_per_step = try checkedDifferenceIncrement(
        candidate.hydrogen_silicate_mol_per_step,
        current.hydrogen_silicate_mol_per_step,
        positive.hydrogen_silicate_mol_per_step,
    );
    try addDifference(
        &candidate.nonband_phosphorus,
        current.nonband_phosphorus,
        positive.nonband_phosphorus,
    );
    try addDifference(
        &candidate.band_phosphorus,
        current.band_phosphorus,
        positive.band_phosphorus,
    );
}

fn addMacroporeDifference(
    candidate: *MacroporeSaltFlux,
    current: MacroporeSaltFlux,
    positive: MacroporeSaltFlux,
) !void {
    try addDifference(
        &candidate.primary_ions,
        current.primary_ions,
        positive.primary_ions,
    );
    try addDifference(
        &candidate.metal_complexes,
        current.metal_complexes,
        positive.metal_complexes,
    );
    try addDifference(
        &candidate.nonband_phosphorus,
        current.nonband_phosphorus,
        positive.nonband_phosphorus,
    );
    try addDifference(
        &candidate.band_phosphorus,
        current.band_phosphorus,
        positive.band_phosphorus,
    );
}

fn addDifference(candidate: anytype, current: anytype, positive: anytype) !void {
    inline for (@typeInfo(@TypeOf(candidate.*)).@"struct".fields) |field|
        @field(candidate.*, field.name) = try checkedDifferenceIncrement(
            @field(candidate.*, field.name),
            @field(current, field.name),
            @field(positive, field.name),
        );
}

fn checkedDifferenceIncrement(
    initial: f64,
    current: f64,
    positive: f64,
) !f64 {
    const with_current = initial + current;
    if (!std.math.isFinite(with_current))
        return error.NonFiniteAdjacentSaltResult;
    const result = with_current - positive;
    if (!std.math.isFinite(result))
        return error.NonFiniteAdjacentSaltResult;
    return result;
}

fn validateInputs(inputs: Inputs, state: State) !void {
    if (!std.math.isFinite(inputs.current_layer_thickness_m) or
        !std.math.isFinite(inputs.layer_activity_threshold_m))
    {
        return error.NonFiniteAdjacentSaltInput;
    }
    if (inputs.current_layer_thickness_m < 0 or
        inputs.layer_activity_threshold_m < 0)
    {
        return error.InvalidAdjacentSaltThickness;
    }
    try validateSalt(inputs.current_cell_flux);
    try validateSalt(inputs.positive_neighbor_flux);
    try validateSalt(state.net_flux);
}

fn validateSalt(value: SaltFlux) !void {
    try validateMicropore(value.micropore);
    try validateMacropore(value.macropore);
}

fn validateMicropore(value: MicroporeSaltFlux) !void {
    try validateLeaf(value.primary_ions);
    try validateLeaf(value.metal_complexes);
    if (!std.math.isFinite(value.hydrogen_silicate_mol_per_step))
        return error.NonFiniteAdjacentSaltInput;
    try validateLeaf(value.nonband_phosphorus);
    try validateLeaf(value.band_phosphorus);
}

fn validateMacropore(value: MacroporeSaltFlux) !void {
    try validateLeaf(value.primary_ions);
    try validateLeaf(value.metal_complexes);
    try validateLeaf(value.nonband_phosphorus);
    try validateLeaf(value.band_phosphorus);
}

fn validateLeaf(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name)))
            return error.NonFiniteAdjacentSaltInput;
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| {
        if (@typeInfo(field.type) == .@"struct")
            @field(result, field.name) = filled(field.type, value)
        else
            @field(result, field.name) = value;
    }
    return result;
}

fn expectFilled(actual: anytype, expected: f64) !void {
    inline for (@typeInfo(@TypeOf(actual)).@"struct".fields) |field| {
        if (@typeInfo(field.type) == .@"struct")
            try expectFilled(@field(actual, field.name), expected)
        else
            try std.testing.expectEqual(expected, @field(actual, field.name));
    }
}

fn expectSumZero(first: anytype, second: @TypeOf(first)) !void {
    inline for (@typeInfo(@TypeOf(first)).@"struct".fields) |field| {
        if (@typeInfo(field.type) == .@"struct")
            try expectSumZero(
                @field(first, field.name),
                @field(second, field.name),
            )
        else
            try std.testing.expectEqual(
                @as(f64, 0),
                @field(first, field.name) + @field(second, field.name),
            );
    }
}

fn baseInputs(current: f64, positive: f64) Inputs {
    return .{
        .transport_axis = .column,
        .horizontal_exchange = .enabled,
        .salt_equilibrium_mode = .dynamic_transport,
        .current_layer_thickness_m = 0.2,
        .layer_activity_threshold_m = 0.1,
        .current_cell_flux = filled(SaltFlux, current),
        .positive_neighbor_flux = filled(SaltFlux, positive),
    };
}

test "pool layout preserves 50 micropore and 49 macropore inventories" {
    const primary_count = @typeInfo(PrimaryIonFlux).@"struct".fields.len;
    const complex_count = @typeInfo(MetalComplexFlux).@"struct".fields.len;
    const phosphorus_count =
        @typeInfo(PhosphorusComplexFlux).@"struct".fields.len;
    try std.testing.expectEqual(@as(usize, 12), primary_count);
    try std.testing.expectEqual(@as(usize, 21), complex_count);
    try std.testing.expectEqual(@as(usize, 8), phosphorus_count);
    try std.testing.expectEqual(
        @as(usize, 50),
        primary_count + complex_count + 1 + 2 * phosphorus_count,
    );
    try std.testing.expectEqual(
        @as(usize, 49),
        primary_count + complex_count + 2 * phosphorus_count,
    );
}

test "all dynamic salt pools follow current minus positive-neighbor flux" {
    var state = State{ .net_flux = filled(SaltFlux, 100) };
    try aggregate(baseInputs(5, 2), &state);
    try expectFilled(state.net_flux, 103);
}

test "reversed adjacent salt faces are exactly antisymmetric" {
    var first = State{ .net_flux = .{} };
    var second = State{ .net_flux = .{} };
    try aggregate(baseInputs(7, 3), &first);
    try aggregate(baseInputs(3, 7), &second);
    try expectSumZero(first.net_flux, second.net_flux);
}

test "equal salt faces produce exact zero increment" {
    var state = State{ .net_flux = .{} };
    try aggregate(baseInputs(9, 9), &state);
    try expectFilled(state.net_flux, 0);
}

test "static standalone inactive and vertical modes preserve source gates" {
    var state = State{ .net_flux = filled(SaltFlux, 9) };
    var inputs = baseInputs(100, 0);
    inputs.salt_equilibrium_mode = .static_concentrations;
    inputs.current_layer_thickness_m = std.math.nan(f64);
    inputs.current_cell_flux = filled(SaltFlux, std.math.nan(f64));
    try aggregate(inputs, &state);
    try expectFilled(state.net_flux, 9);

    inputs.salt_equilibrium_mode = .dynamic_transport;
    inputs.horizontal_exchange = .standalone;
    try aggregate(inputs, &state);
    try expectFilled(state.net_flux, 9);

    inputs.horizontal_exchange = .enabled;
    inputs.current_layer_thickness_m = 0.1;
    inputs.layer_activity_threshold_m = 0.1;
    inputs.current_cell_flux = filled(SaltFlux, 100);
    try aggregate(inputs, &state);
    try expectFilled(state.net_flux, 9);

    inputs.transport_axis = .vertical;
    inputs.horizontal_exchange = .standalone;
    inputs.current_layer_thickness_m = 0.2;
    try aggregate(inputs, &state);
    try expectFilled(state.net_flux, 109);
}

test "invalid and overflow failures preserve salt state atomically" {
    var state = State{ .net_flux = filled(SaltFlux, 5) };
    var inputs = baseInputs(3, 1);
    inputs.current_cell_flux.macropore.band_phosphorus
        .magnesium_hydrogen_phosphate_mol_per_step = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteAdjacentSaltInput,
        aggregate(inputs, &state),
    );
    try expectFilled(state.net_flux, 5);

    inputs.current_cell_flux = filled(SaltFlux, std.math.floatMax(f64));
    inputs.positive_neighbor_flux =
        filled(SaltFlux, -std.math.floatMax(f64));
    state.net_flux = filled(SaltFlux, std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteAdjacentSaltResult,
        aggregate(inputs, &state),
    );
    try expectFilled(state.net_flux, std.math.floatMax(f64));
}
