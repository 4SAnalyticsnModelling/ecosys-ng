const std = @import("std");

pub const DisturbanceMode = enum {
    no_profile_effects,
    freeze_thaw,
    freeze_thaw_and_erosion,
    freeze_thaw_and_organic_matter,
    freeze_thaw_erosion_and_organic_matter,
};

pub const TransportAxis = enum {
    east_west,
    north_south,
    vertical,
};

/// Signed fertilizer-pool erosion flux, mol per model step.
pub const FertilizerFlux = struct {
    nonband_ammonium_mol_per_step: f64 = 0,
    nonband_ammonia_mol_per_step: f64 = 0,
    nonband_urea_mol_per_step: f64 = 0,
    nonband_nitrate_mol_per_step: f64 = 0,
    band_ammonium_mol_per_step: f64 = 0,
    band_ammonia_mol_per_step: f64 = 0,
    band_urea_mol_per_step: f64 = 0,
    band_nitrate_mol_per_step: f64 = 0,
};

pub const BoundaryFlux = struct {
    total_sediment_Mg_per_step: f64,
    fertilizer: FertilizerFlux,
};

pub const Inputs = struct {
    disturbance_mode: DisturbanceMode,
    transport_axis: TransportAxis,
    sediment_activity_threshold_Mg_per_step: f64,
    local_flux_by_boundary_side: []const BoundaryFlux,
    positive_neighbor_total_sediment_Mg_per_step_by_boundary_side: []const f64,
};

pub const State = struct {
    net_incoming: FertilizerFlux,
};

/// Accounts local fertilizer-pool erosion on one horizontal axis.
///
/// Traceability: REDIST.F lines 2922--2929, under the erosion/axis/side gate
/// at lines 2888--2895. Runtime boundary-side slices replace fixed `NN=1,2`.
/// An active side adds its eight local non-band/band molar inventories in
/// source order. Candidate state commits only after every addition is finite.
pub fn account(inputs: Inputs, state: *State) !void {
    if (!erosionEnabled(inputs.disturbance_mode) or
        inputs.transport_axis == .vertical)
    {
        return;
    }
    try validateInputs(inputs, state.*);

    var candidate = state.net_incoming;
    for (
        inputs.local_flux_by_boundary_side,
        inputs.positive_neighbor_total_sediment_Mg_per_step_by_boundary_side,
    ) |local, positive_neighbor_total| {
        if (@abs(local.total_sediment_Mg_per_step) <=
            inputs.sediment_activity_threshold_Mg_per_step and
            @abs(positive_neighbor_total) <=
                inputs.sediment_activity_threshold_Mg_per_step)
        {
            continue;
        }
        try addFlux(&candidate, local.fertilizer);
    }
    state.net_incoming = candidate;
}

fn erosionEnabled(mode: DisturbanceMode) bool {
    return mode == .freeze_thaw_and_erosion or
        mode == .freeze_thaw_erosion_and_organic_matter;
}

fn validateInputs(inputs: Inputs, state: State) !void {
    if (inputs.local_flux_by_boundary_side.len == 0)
        return error.InvalidFertilizerErosionDimensions;
    if (inputs.positive_neighbor_total_sediment_Mg_per_step_by_boundary_side.len !=
        inputs.local_flux_by_boundary_side.len)
    {
        return error.FertilizerErosionDimensionMismatch;
    }
    if (!std.math.isFinite(inputs.sediment_activity_threshold_Mg_per_step))
        return error.NonFiniteFertilizerErosionInput;
    if (inputs.sediment_activity_threshold_Mg_per_step < 0)
        return error.InvalidErosionActivityThreshold;
    try validateFlux(state.net_incoming);
    for (inputs.local_flux_by_boundary_side) |flux| {
        if (!std.math.isFinite(flux.total_sediment_Mg_per_step))
            return error.NonFiniteFertilizerErosionInput;
        try validateFlux(flux.fertilizer);
    }
    for (inputs.positive_neighbor_total_sediment_Mg_per_step_by_boundary_side) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteFertilizerErosionInput;
}

fn validateFlux(flux: FertilizerFlux) !void {
    inline for (@typeInfo(FertilizerFlux).@"struct".fields) |field|
        if (!std.math.isFinite(@field(flux, field.name)))
            return error.NonFiniteFertilizerErosionInput;
}

fn addFlux(candidate: *FertilizerFlux, contribution: FertilizerFlux) !void {
    inline for (@typeInfo(FertilizerFlux).@"struct".fields) |field| {
        const result = @field(candidate.*, field.name) +
            @field(contribution, field.name);
        if (!std.math.isFinite(result))
            return error.NonFiniteFertilizerErosionResult;
        @field(candidate.*, field.name) = result;
    }
}

fn filledFlux(value: f64) FertilizerFlux {
    var result: FertilizerFlux = undefined;
    inline for (@typeInfo(FertilizerFlux).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn expectFlux(actual: FertilizerFlux, expected: f64) !void {
    inline for (@typeInfo(FertilizerFlux).@"struct".fields) |field|
        try std.testing.expectEqual(expected, @field(actual, field.name));
}

test "active erosion sides add all nonband and band fertilizer pools" {
    const local = [_]BoundaryFlux{
        .{ .total_sediment_Mg_per_step = 2, .fertilizer = filledFlux(3) },
        .{ .total_sediment_Mg_per_step = -4, .fertilizer = filledFlux(5) },
    };
    const positive = [_]f64{ 0, 0 };
    var state = State{ .net_incoming = filledFlux(100) };

    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &positive,
    }, &state);
    try expectFlux(state.net_incoming, 108);
}

test "strict threshold and positive-neighbor activation match source gate" {
    const local = [_]BoundaryFlux{
        .{ .total_sediment_Mg_per_step = 1, .fertilizer = filledFlux(2) },
        .{ .total_sediment_Mg_per_step = 1, .fertilizer = filledFlux(3) },
    };
    const positive = [_]f64{ 2, 1 };
    var state = State{ .net_incoming = .{} };

    try account(.{
        .disturbance_mode = .freeze_thaw_erosion_and_organic_matter,
        .transport_axis = .north_south,
        .sediment_activity_threshold_Mg_per_step = 1,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &positive,
    }, &state);
    try expectFlux(state.net_incoming, 2);
}

test "runtime accepted-side sum conserves every fertilizer inventory exactly" {
    const local = [_]BoundaryFlux{
        .{ .total_sediment_Mg_per_step = 2, .fertilizer = filledFlux(2) },
        .{ .total_sediment_Mg_per_step = -3, .fertilizer = filledFlux(-3) },
        .{ .total_sediment_Mg_per_step = 4, .fertilizer = filledFlux(4) },
        .{ .total_sediment_Mg_per_step = 5, .fertilizer = filledFlux(5) },
    };
    const positive = [_]f64{ 0, 0, 0, 0 };
    var state = State{ .net_incoming = .{} };

    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &positive,
    }, &state);
    try expectFlux(state.net_incoming, 8);
}

test "disabled erosion and vertical axes bypass unused inputs" {
    var state = State{ .net_incoming = filledFlux(9) };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_organic_matter,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
        .local_flux_by_boundary_side = &.{},
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &.{},
    }, &state);
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .vertical,
        .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
        .local_flux_by_boundary_side = &.{},
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &.{},
    }, &state);
    try expectFlux(state.net_incoming, 9);
}

test "invalid input overflow and dimensions preserve state atomically" {
    var local = [_]BoundaryFlux{
        .{ .total_sediment_Mg_per_step = 2, .fertilizer = filledFlux(2) },
        .{ .total_sediment_Mg_per_step = 3, .fertilizer = filledFlux(3) },
    };
    local[1].fertilizer.band_nitrate_mol_per_step = std.math.nan(f64);
    const positive = [_]f64{ 0, 0 };
    var state = State{ .net_incoming = filledFlux(4) };
    const inputs = Inputs{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &positive,
    };

    try std.testing.expectError(
        error.NonFiniteFertilizerErosionInput,
        account(inputs, &state),
    );
    try expectFlux(state.net_incoming, 4);

    local[1].fertilizer = filledFlux(std.math.floatMax(f64));
    state.net_incoming = filledFlux(std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteFertilizerErosionResult,
        account(inputs, &state),
    );
    try expectFlux(state.net_incoming, std.math.floatMax(f64));

    const short_positive = [_]f64{0};
    try std.testing.expectError(
        error.FertilizerErosionDimensionMismatch,
        account(.{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .east_west,
            .sediment_activity_threshold_Mg_per_step = 1,
            .local_flux_by_boundary_side = &local,
            .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &short_positive,
        }, &state),
    );
    try expectFlux(state.net_incoming, std.math.floatMax(f64));
}
