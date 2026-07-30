const std = @import("std");

pub const SaltEquilibriumMode = enum {
    static,
    dynamic,
};

pub const ExternalBoundaryWindow = struct {
    first_column: usize,
    last_column_inclusive: usize,
    first_row: usize,
    last_row_inclusive: usize,
};

pub const Face = enum(u8) {
    east,
    west,
    south,
    north,

    pub const count: usize = @typeInfo(Face).@"enum".fields.len;
};

pub const MultiplicityOneFlux = struct {
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
    po4_mol_p_per_step: f64 = 0,
};

pub const MultiplicityTwoFlux = struct {
    bicarbonate_mol_per_step: f64 = 0,
    aluminum_hydroxide_1_mol_per_step: f64 = 0,
    aluminum_sulfate_mol_per_step: f64 = 0,
    iron_hydroxide_1_mol_per_step: f64 = 0,
    iron_sulfate_mol_per_step: f64 = 0,
    calcium_hydroxide_mol_per_step: f64 = 0,
    calcium_carbonate_mol_per_step: f64 = 0,
    calcium_sulfate_mol_per_step: f64 = 0,
    magnesium_hydroxide_mol_per_step: f64 = 0,
    magnesium_carbonate_mol_per_step: f64 = 0,
    magnesium_sulfate_mol_per_step: f64 = 0,
    sodium_carbonate_mol_per_step: f64 = 0,
    sodium_sulfate_mol_per_step: f64 = 0,
    potassium_sulfate_mol_per_step: f64 = 0,
    calcium_po4_mol_p_per_step: f64 = 0,
};

pub const MultiplicityThreeFlux = struct {
    aluminum_hydroxide_2_mol_per_step: f64 = 0,
    iron_hydroxide_2_mol_per_step: f64 = 0,
    calcium_bicarbonate_mol_per_step: f64 = 0,
    magnesium_bicarbonate_mol_per_step: f64 = 0,
    iron_hpo4_mol_p_per_step: f64 = 0,
    calcium_hpo4_mol_p_per_step: f64 = 0,
    magnesium_hpo4_mol_p_per_step: f64 = 0,
};

pub const MultiplicityFourFlux = struct {
    aluminum_hydroxide_3_mol_per_step: f64 = 0,
    iron_hydroxide_3_mol_per_step: f64 = 0,
    phosphoric_acid_mol_p_per_step: f64 = 0,
    iron_h2po4_mol_p_per_step: f64 = 0,
    calcium_h2po4_mol_p_per_step: f64 = 0,
};

pub const MultiplicityFiveFlux = struct {
    aluminum_hydroxide_4_mol_per_step: f64 = 0,
    iron_hydroxide_4_mol_per_step: f64 = 0,
};

pub const BoundaryFlux = struct {
    water_m3_per_step: f64 = 0,
    nitrate_nitrogen_g_n_per_step: f64 = 0,
    /// Source `XQSH0P`: snow-carried PO4 included by REDIST in this runoff
    /// phosphate mass expression and again in the later snow block.
    snow_po4_mol_p_per_step: f64 = 0,
    multiplicity_one: MultiplicityOneFlux = .{},
    multiplicity_two: MultiplicityTwoFlux = .{},
    multiplicity_three: MultiplicityThreeFlux = .{},
    multiplicity_four: MultiplicityFourFlux = .{},
    multiplicity_five: MultiplicityFiveFlux = .{},
};

pub const ConductivityCoefficients = struct {
    hydrogen_dS_m2_per_mol: f64,
    hydroxide_dS_m2_per_mol: f64,
    aluminum_dS_m2_per_equivalent: f64,
    iron_dS_m2_per_equivalent: f64,
    calcium_dS_m2_per_equivalent: f64,
    magnesium_dS_m2_per_equivalent: f64,
    sodium_dS_m2_per_mol: f64,
    potassium_dS_m2_per_mol: f64,
    carbonate_dS_m2_per_equivalent: f64,
    bicarbonate_dS_m2_per_mol: f64,
    sulfate_dS_m2_per_equivalent: f64,
    chloride_dS_m2_per_mol: f64,
    nitrate_dS_m2_per_mol: f64,
};

pub const Inputs = struct {
    grid_column_count: usize,
    grid_row_count: usize,
    external_boundary_window: ExternalBoundaryWindow,
    salt_equilibrium_mode: SaltEquilibriumMode,
    flux_by_face: [Face.count][]const BoundaryFlux,
    negligible_water_m3_by_cell: []const f64,
    phosphorus_g_p_per_mol_p: f64,
    nitrogen_g_n_per_mol_n: f64,
    conductivity: ConductivityCoefficients,
};

pub const State = struct {
    cumulative_phosphorus_output_g_p: f64,
    cumulative_ion_component_output_mol: f64,
    ion_component_output_mol_by_cell: []f64,
    last_runoff_electrical_conductivity_dS_per_m_by_cell: []f64,
};

/// Accounts for dynamic runoff salts and calculates runoff conductivity.
///
/// Traceability: REDIST.F lines 785--846 (`PSS`, `SS1`--`SS4`, `SSR`,
/// `TPOU`, `TIONOU`, `UIONOU`, and `ECNDQ`). The exact 1/2/3/4/5 component
/// multiplicities and source addition order are retained. The operation is
/// gated by dynamic salt equilibrium and the enclosing strict water-volume
/// threshold. Runtime coefficients replace the source conductivity literals.
pub fn apply(inputs: Inputs, state: *State) !void {
    const cell_count = try validateDimensions(inputs, state.*);
    try validateInputsAndState(inputs, state.*, cell_count);
    if (inputs.salt_equilibrium_mode == .static) return;
    try preflightUpdates(inputs, state.*);

    const window = inputs.external_boundary_window;
    for (window.first_column..window.last_column_inclusive + 1) |column| {
        for (window.first_row..window.last_row_inclusive + 1) |row| {
            const cell = row * inputs.grid_column_count + column;
            for (std.meta.tags(Face)) |face| {
                if (!isExternalFace(face, column, row, window)) continue;
                const flux = inputs.flux_by_face[@intFromEnum(face)][cell];
                if (@abs(flux.water_m3_per_step) <=
                    inputs.negligible_water_m3_by_cell[cell]) continue;
                const result = calculate(flux, inwardSign(face), inputs) catch unreachable;
                state.cumulative_phosphorus_output_g_p -= result.phosphorus_g_p;
                state.cumulative_ion_component_output_mol -=
                    result.ion_components_mol;
                state.ion_component_output_mol_by_cell[cell] -=
                    result.ion_components_mol;
                state.last_runoff_electrical_conductivity_dS_per_m_by_cell[cell] =
                    result.electrical_conductivity_dS_per_m;
            }
        }
    }
}

const Result = struct {
    phosphorus_g_p: f64,
    ion_components_mol: f64,
    electrical_conductivity_dS_per_m: f64,
};

fn calculate(flux: BoundaryFlux, direction: f64, inputs: Inputs) !Result {
    var phosphorus_mol_p = try sourceSum(&.{
        flux.multiplicity_one.po4_mol_p_per_step,
        flux.multiplicity_two.calcium_po4_mol_p_per_step,
        flux.multiplicity_three.iron_hpo4_mol_p_per_step,
        flux.multiplicity_three.calcium_hpo4_mol_p_per_step,
        flux.multiplicity_three.magnesium_hpo4_mol_p_per_step,
        flux.multiplicity_four.phosphoric_acid_mol_p_per_step,
        flux.multiplicity_four.iron_h2po4_mol_p_per_step,
        flux.multiplicity_four.calcium_h2po4_mol_p_per_step,
    });
    var phosphorus_g_p = try checkedProduct(direction, inputs.phosphorus_g_p_per_mol_p);
    phosphorus_g_p = try checkedProduct(phosphorus_g_p, phosphorus_mol_p);
    phosphorus_mol_p = flux.snow_po4_mol_p_per_step;
    phosphorus_g_p = try checkedSum(
        phosphorus_g_p,
        try checkedProduct(
            try checkedProduct(direction, inputs.phosphorus_g_p_per_mol_p),
            phosphorus_mol_p,
        ),
    );

    const group_one = try checkedProduct(
        direction,
        try sourceStructSum(flux.multiplicity_one),
    );
    const group_two = try checkedProduct(
        try checkedProduct(direction, 2),
        try sourceStructSum(flux.multiplicity_two),
    );
    const group_three = try checkedProduct(
        try checkedProduct(direction, 3),
        try sourceStructSum(flux.multiplicity_three),
    );
    var group_four = try checkedProduct(
        try checkedProduct(direction, 4),
        try sourceStructSum(flux.multiplicity_four),
    );
    group_four = try checkedSum(
        group_four,
        try checkedProduct(
            try checkedProduct(direction, 5),
            try sourceStructSum(flux.multiplicity_five),
        ),
    );
    var ion_components_mol = try checkedSum(group_one, group_two);
    ion_components_mol = try checkedSum(ion_components_mol, group_three);
    ion_components_mol = try checkedSum(ion_components_mol, group_four);

    return .{
        .phosphorus_g_p = phosphorus_g_p,
        .ion_components_mol = ion_components_mol,
        .electrical_conductivity_dS_per_m = try electricalConductivity(flux, inputs),
    };
}

fn electricalConductivity(flux: BoundaryFlux, inputs: Inputs) !f64 {
    const water = flux.water_m3_per_step;
    const one = flux.multiplicity_one;
    const two = flux.multiplicity_two;
    const coefficient = inputs.conductivity;
    var total: f64 = 0;
    total = try checkedSum(total, coefficient.hydrogen_dS_m2_per_mol * @max(0, one.hydrogen_mol_per_step / water));
    total = try checkedSum(total, coefficient.hydroxide_dS_m2_per_mol * @max(0, one.hydroxide_mol_per_step / water));
    total = try checkedSum(total, coefficient.aluminum_dS_m2_per_equivalent * @max(0, one.aluminum_mol_per_step * 3 / water));
    total = try checkedSum(total, coefficient.iron_dS_m2_per_equivalent * @max(0, one.iron_mol_per_step * 3 / water));
    total = try checkedSum(total, coefficient.calcium_dS_m2_per_equivalent * @max(0, one.calcium_mol_per_step * 2 / water));
    total = try checkedSum(total, coefficient.magnesium_dS_m2_per_equivalent * @max(0, one.magnesium_mol_per_step * 2 / water));
    total = try checkedSum(total, coefficient.sodium_dS_m2_per_mol * @max(0, one.sodium_mol_per_step / water));
    total = try checkedSum(total, coefficient.potassium_dS_m2_per_mol * @max(0, one.potassium_mol_per_step / water));
    total = try checkedSum(total, coefficient.carbonate_dS_m2_per_equivalent * @max(0, one.carbonate_mol_per_step * 2 / water));
    total = try checkedSum(total, coefficient.bicarbonate_dS_m2_per_mol * @max(0, two.bicarbonate_mol_per_step / water));
    total = try checkedSum(total, coefficient.sulfate_dS_m2_per_equivalent * @max(0, one.sulfate_mol_per_step * 2 / water));
    total = try checkedSum(total, coefficient.chloride_dS_m2_per_mol * @max(0, one.chloride_mol_per_step / water));
    total = try checkedSum(
        total,
        coefficient.nitrate_dS_m2_per_mol *
            @max(0, flux.nitrate_nitrogen_g_n_per_step /
                (water * inputs.nitrogen_g_n_per_mol_n)),
    );
    return total;
}

fn preflightUpdates(inputs: Inputs, state: State) !void {
    var cumulative_phosphorus_output_g_p =
        state.cumulative_phosphorus_output_g_p;
    var cumulative_ion_component_output_mol =
        state.cumulative_ion_component_output_mol;
    const window = inputs.external_boundary_window;
    for (window.first_column..window.last_column_inclusive + 1) |column| {
        for (window.first_row..window.last_row_inclusive + 1) |row| {
            const cell = row * inputs.grid_column_count + column;
            var cell_ion_output_mol = state.ion_component_output_mol_by_cell[cell];
            for (std.meta.tags(Face)) |face| {
                if (!isExternalFace(face, column, row, window)) continue;
                const flux = inputs.flux_by_face[@intFromEnum(face)][cell];
                if (@abs(flux.water_m3_per_step) <=
                    inputs.negligible_water_m3_by_cell[cell]) continue;
                const result = try calculate(flux, inwardSign(face), inputs);
                cumulative_phosphorus_output_g_p = try checkedSum(
                    cumulative_phosphorus_output_g_p,
                    -result.phosphorus_g_p,
                );
                cumulative_ion_component_output_mol = try checkedSum(
                    cumulative_ion_component_output_mol,
                    -result.ion_components_mol,
                );
                cell_ion_output_mol =
                    try checkedSum(cell_ion_output_mol, -result.ion_components_mol);
            }
        }
    }
}

fn validateDimensions(inputs: Inputs, state: State) !usize {
    if (inputs.grid_column_count == 0 or inputs.grid_row_count == 0)
        return error.InvalidRunoffSaltDimensions;
    const cell_count = std.math.mul(
        usize,
        inputs.grid_column_count,
        inputs.grid_row_count,
    ) catch return error.InvalidRunoffSaltDimensions;
    inline for (inputs.flux_by_face) |values|
        if (values.len != cell_count) return error.InvalidRunoffSaltDimensions;
    if (inputs.negligible_water_m3_by_cell.len != cell_count or
        state.ion_component_output_mol_by_cell.len != cell_count or
        state.last_runoff_electrical_conductivity_dS_per_m_by_cell.len != cell_count)
        return error.InvalidRunoffSaltDimensions;
    const window = inputs.external_boundary_window;
    if (window.first_column > window.last_column_inclusive or
        window.first_row > window.last_row_inclusive or
        window.last_column_inclusive >= inputs.grid_column_count or
        window.last_row_inclusive >= inputs.grid_row_count)
        return error.InvalidRunoffSaltBoundaryWindow;
    return cell_count;
}

fn validateInputsAndState(inputs: Inputs, state: State, cell_count: usize) !void {
    if (!positiveFinite(inputs.phosphorus_g_p_per_mol_p) or
        !positiveFinite(inputs.nitrogen_g_n_per_mol_n))
        return error.InvalidRunoffSaltInput;
    try validatePositiveStruct(inputs.conductivity, error.InvalidRunoffSaltInput);
    inline for (inputs.flux_by_face) |values| {
        for (values) |flux| {
            inline for (.{
                flux.water_m3_per_step,
                flux.nitrate_nitrogen_g_n_per_step,
                flux.snow_po4_mol_p_per_step,
            }) |value| if (!std.math.isFinite(value))
                return error.InvalidRunoffSaltInput;
            try validateFiniteStruct(flux.multiplicity_one, error.InvalidRunoffSaltInput);
            try validateFiniteStruct(flux.multiplicity_two, error.InvalidRunoffSaltInput);
            try validateFiniteStruct(flux.multiplicity_three, error.InvalidRunoffSaltInput);
            try validateFiniteStruct(flux.multiplicity_four, error.InvalidRunoffSaltInput);
            try validateFiniteStruct(flux.multiplicity_five, error.InvalidRunoffSaltInput);
        }
    }
    if (!std.math.isFinite(state.cumulative_phosphorus_output_g_p) or
        !std.math.isFinite(state.cumulative_ion_component_output_mol))
        return error.InvalidRunoffSaltState;
    for (0..cell_count) |cell| {
        if (!nonnegativeFinite(inputs.negligible_water_m3_by_cell[cell]))
            return error.InvalidRunoffSaltInput;
        if (!std.math.isFinite(state.ion_component_output_mol_by_cell[cell]) or
            !nonnegativeFinite(
                state.last_runoff_electrical_conductivity_dS_per_m_by_cell[cell],
            ))
            return error.InvalidRunoffSaltState;
    }
}

fn isExternalFace(
    face: Face,
    column: usize,
    row: usize,
    window: ExternalBoundaryWindow,
) bool {
    return switch (face) {
        .east => column == window.last_column_inclusive,
        .west => column == window.first_column,
        .south => row == window.last_row_inclusive,
        .north => row == window.first_row,
    };
}

fn inwardSign(face: Face) f64 {
    return switch (face) {
        .east, .south => -1,
        .west, .north => 1,
    };
}

fn sourceStructSum(value: anytype) !f64 {
    var total: f64 = 0;
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        total = try checkedSum(total, @field(value, field.name));
    return total;
}

fn sourceSum(values: []const f64) !f64 {
    var total: f64 = 0;
    for (values) |value| total = try checkedSum(total, value);
    return total;
}

fn checkedProduct(left: f64, right: f64) !f64 {
    const result = left * right;
    if (!std.math.isFinite(result)) return error.NonFiniteRunoffSaltResult;
    return result;
}

fn checkedSum(current: f64, increment: f64) !f64 {
    const result = current + increment;
    if (!std.math.isFinite(result)) return error.NonFiniteRunoffSaltResult;
    return result;
}

fn validateFiniteStruct(value: anytype, comptime failure: anyerror) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name))) return failure;
}

fn validatePositiveStruct(value: anytype, comptime failure: anyerror) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!positiveFinite(@field(value, field.name))) return failure;
}

fn positiveFinite(value: f64) bool {
    return std.math.isFinite(value) and value > 0;
}

fn nonnegativeFinite(value: f64) bool {
    return std.math.isFinite(value) and value >= 0;
}

fn allOneFlux() BoundaryFlux {
    return .{
        .water_m3_per_step = 1,
        .nitrate_nitrogen_g_n_per_step = 14,
        .snow_po4_mol_p_per_step = 1,
        .multiplicity_one = filled(MultiplicityOneFlux, 1),
        .multiplicity_two = filled(MultiplicityTwoFlux, 1),
        .multiplicity_three = filled(MultiplicityThreeFlux, 1),
        .multiplicity_four = filled(MultiplicityFourFlux, 1),
        .multiplicity_five = filled(MultiplicityFiveFlux, 1),
    };
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn unitConductivity() ConductivityCoefficients {
    return filled(ConductivityCoefficients, 1);
}

fn sourceConductivity() ConductivityCoefficients {
    return .{
        .hydrogen_dS_m2_per_mol = 0.337,
        .hydroxide_dS_m2_per_mol = 0.192,
        .aluminum_dS_m2_per_equivalent = 0.056,
        .iron_dS_m2_per_equivalent = 0.051,
        .calcium_dS_m2_per_equivalent = 0.060,
        .magnesium_dS_m2_per_equivalent = 0.053,
        .sodium_dS_m2_per_mol = 0.050,
        .potassium_dS_m2_per_mol = 0.070,
        .carbonate_dS_m2_per_equivalent = 0.072,
        .bicarbonate_dS_m2_per_mol = 0.044,
        .sulfate_dS_m2_per_equivalent = 0.080,
        .chloride_dS_m2_per_mol = 0.076,
        .nitrate_dS_m2_per_mol = 0.071,
    };
}

test "REDIST dynamic runoff salt groups preserve source multiplicities" {
    const east = [_]BoundaryFlux{allOneFlux()};
    const zero = [_]BoundaryFlux{.{}};
    var ion_output = [_]f64{0};
    var conductivity = [_]f64{0};
    var state: State = .{
        .cumulative_phosphorus_output_g_p = 0,
        .cumulative_ion_component_output_mol = 0,
        .ion_component_output_mol_by_cell = &ion_output,
        .last_runoff_electrical_conductivity_dS_per_m_by_cell = &conductivity,
    };
    try apply(.{
        .grid_column_count = 1,
        .grid_row_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .salt_equilibrium_mode = .dynamic,
        .flux_by_face = .{ &east, &zero, &zero, &zero },
        .negligible_water_m3_by_cell = &.{0},
        .phosphorus_g_p_per_mol_p = 31,
        .nitrogen_g_n_per_mol_n = 14,
        .conductivity = sourceConductivity(),
    }, &state);
    try std.testing.expectEqual(@as(f64, 279), state.cumulative_phosphorus_output_g_p);
    try std.testing.expectEqual(@as(f64, 93), state.cumulative_ion_component_output_mol);
    try std.testing.expectEqual(@as(f64, 93), ion_output[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 1.691), conductivity[0], 1.0e-15);
}

test "west face reverses ledgers while signed carrier ratios retain conductivity" {
    var west_flux = allOneFlux();
    west_flux.water_m3_per_step = -1;
    west_flux.nitrate_nitrogen_g_n_per_step = -14;
    west_flux.snow_po4_mol_p_per_step = -1;
    west_flux.multiplicity_one = filled(MultiplicityOneFlux, -1);
    west_flux.multiplicity_two = filled(MultiplicityTwoFlux, -1);
    west_flux.multiplicity_three = filled(MultiplicityThreeFlux, -1);
    west_flux.multiplicity_four = filled(MultiplicityFourFlux, -1);
    west_flux.multiplicity_five = filled(MultiplicityFiveFlux, -1);
    const west = [_]BoundaryFlux{west_flux};
    const zero = [_]BoundaryFlux{.{}};
    var ion_output = [_]f64{100};
    var conductivity = [_]f64{0};
    var state: State = .{
        .cumulative_phosphorus_output_g_p = 300,
        .cumulative_ion_component_output_mol = 100,
        .ion_component_output_mol_by_cell = &ion_output,
        .last_runoff_electrical_conductivity_dS_per_m_by_cell = &conductivity,
    };
    try apply(.{
        .grid_column_count = 1,
        .grid_row_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .salt_equilibrium_mode = .dynamic,
        .flux_by_face = .{ &zero, &west, &zero, &zero },
        .negligible_water_m3_by_cell = &.{0},
        .phosphorus_g_p_per_mol_p = 31,
        .nitrogen_g_n_per_mol_n = 14,
        .conductivity = unitConductivity(),
    }, &state);
    try std.testing.expectEqual(@as(f64, 579), state.cumulative_phosphorus_output_g_p);
    try std.testing.expectEqual(@as(f64, 193), state.cumulative_ion_component_output_mol);
    try std.testing.expectEqual(@as(f64, 193), ion_output[0]);
    try std.testing.expectEqual(@as(f64, 21), conductivity[0]);
}

test "static mode and equal water threshold leave salt state unchanged" {
    const east = [_]BoundaryFlux{allOneFlux()};
    const zero = [_]BoundaryFlux{.{}};
    var ion_output = [_]f64{2};
    var conductivity = [_]f64{3};
    var state: State = .{
        .cumulative_phosphorus_output_g_p = 4,
        .cumulative_ion_component_output_mol = 5,
        .ion_component_output_mol_by_cell = &ion_output,
        .last_runoff_electrical_conductivity_dS_per_m_by_cell = &conductivity,
    };
    const base: Inputs = .{
        .grid_column_count = 1,
        .grid_row_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .salt_equilibrium_mode = .static,
        .flux_by_face = .{ &east, &zero, &zero, &zero },
        .negligible_water_m3_by_cell = &.{0},
        .phosphorus_g_p_per_mol_p = 31,
        .nitrogen_g_n_per_mol_n = 14,
        .conductivity = unitConductivity(),
    };
    try apply(base, &state);
    var threshold = base;
    threshold.salt_equilibrium_mode = .dynamic;
    threshold.negligible_water_m3_by_cell = &.{1};
    try apply(threshold, &state);
    try std.testing.expectEqual(@as(f64, 4), state.cumulative_phosphorus_output_g_p);
    try std.testing.expectEqual(@as(f64, 5), state.cumulative_ion_component_output_mol);
    try std.testing.expectEqual(@as(f64, 2), ion_output[0]);
    try std.testing.expectEqual(@as(f64, 3), conductivity[0]);
}

test "runtime grid indexing conserves the boundary ion ledger" {
    const east = [_]BoundaryFlux{ .{}, allOneFlux(), .{}, allOneFlux() };
    const zero = [_]BoundaryFlux{ .{}, .{}, .{}, .{} };
    var ion_output = [_]f64{ 0, 0, 0, 0 };
    var conductivity = [_]f64{ 0, 0, 0, 0 };
    var state: State = .{
        .cumulative_phosphorus_output_g_p = 0,
        .cumulative_ion_component_output_mol = 0,
        .ion_component_output_mol_by_cell = &ion_output,
        .last_runoff_electrical_conductivity_dS_per_m_by_cell = &conductivity,
    };
    try apply(.{
        .grid_column_count = 2,
        .grid_row_count = 2,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 1,
            .first_row = 0,
            .last_row_inclusive = 1,
        },
        .salt_equilibrium_mode = .dynamic,
        .flux_by_face = .{ &east, &zero, &zero, &zero },
        .negligible_water_m3_by_cell = &.{ 0, 0, 0, 0 },
        .phosphorus_g_p_per_mol_p = 31,
        .nitrogen_g_n_per_mol_n = 14,
        .conductivity = unitConductivity(),
    }, &state);

    try std.testing.expectEqual(@as(f64, 558), state.cumulative_phosphorus_output_g_p);
    try std.testing.expectEqual(@as(f64, 186), state.cumulative_ion_component_output_mol);
    try std.testing.expectEqualSlices(f64, &.{ 0, 93, 0, 93 }, &ion_output);
    try std.testing.expectEqualSlices(f64, &.{ 0, 21, 0, 21 }, &conductivity);
}

test "late invalid salt flux leaves every ledger unchanged" {
    var east = [_]BoundaryFlux{ allOneFlux(), allOneFlux() };
    east[1].multiplicity_five.iron_hydroxide_4_mol_per_step =
        std.math.nan(f64);
    const zero = [_]BoundaryFlux{ .{}, .{} };
    var ion_output = [_]f64{ 1, 2 };
    var conductivity = [_]f64{ 3, 4 };
    var state: State = .{
        .cumulative_phosphorus_output_g_p = 5,
        .cumulative_ion_component_output_mol = 6,
        .ion_component_output_mol_by_cell = &ion_output,
        .last_runoff_electrical_conductivity_dS_per_m_by_cell = &conductivity,
    };
    try std.testing.expectError(error.InvalidRunoffSaltInput, apply(.{
        .grid_column_count = 2,
        .grid_row_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 1,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .salt_equilibrium_mode = .dynamic,
        .flux_by_face = .{ &east, &zero, &zero, &zero },
        .negligible_water_m3_by_cell = &.{ 0, 0 },
        .phosphorus_g_p_per_mol_p = 31,
        .nitrogen_g_n_per_mol_n = 14,
        .conductivity = unitConductivity(),
    }, &state));
    try std.testing.expectEqual(@as(f64, 5), state.cumulative_phosphorus_output_g_p);
    try std.testing.expectEqual(@as(f64, 6), state.cumulative_ion_component_output_mol);
    try std.testing.expectEqualSlices(f64, &.{ 1, 2 }, &ion_output);
    try std.testing.expectEqualSlices(f64, &.{ 3, 4 }, &conductivity);
}

test "invalid salt runtime window and coefficients fail explicitly" {
    const zero = [_]BoundaryFlux{.{}};
    var ion_output = [_]f64{0};
    var conductivity = [_]f64{0};
    var state: State = .{
        .cumulative_phosphorus_output_g_p = 0,
        .cumulative_ion_component_output_mol = 0,
        .ion_component_output_mol_by_cell = &ion_output,
        .last_runoff_electrical_conductivity_dS_per_m_by_cell = &conductivity,
    };
    var coefficients = unitConductivity();
    coefficients.nitrate_dS_m2_per_mol = 0;
    const inputs: Inputs = .{
        .grid_column_count = 1,
        .grid_row_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .salt_equilibrium_mode = .dynamic,
        .flux_by_face = .{ &zero, &zero, &zero, &zero },
        .negligible_water_m3_by_cell = &.{0},
        .phosphorus_g_p_per_mol_p = 31,
        .nitrogen_g_n_per_mol_n = 14,
        .conductivity = coefficients,
    };
    try std.testing.expectError(error.InvalidRunoffSaltInput, apply(inputs, &state));

    var invalid_window = inputs;
    invalid_window.conductivity = unitConductivity();
    invalid_window.external_boundary_window.first_column = 1;
    try std.testing.expectError(
        error.InvalidRunoffSaltBoundaryWindow,
        apply(invalid_window, &state),
    );
}
