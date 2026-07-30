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

/// Signed mineral flux entering from one erosion boundary.
pub const MineralFlux = struct {
    total_sediment_Mg_per_step: f64 = 0,
    sand_Mg_per_step: f64 = 0,
    silt_Mg_per_step: f64 = 0,
    clay_Mg_per_step: f64 = 0,
    cation_exchange_capacity_mol_per_step: f64 = 0,
    anion_exchange_capacity_mol_per_step: f64 = 0,
};

pub const Inputs = struct {
    disturbance_mode: DisturbanceMode,
    transport_axis: TransportAxis,
    /// Cell-specific `ZEROS` gate, Mg per model step.
    sediment_activity_threshold_Mg_per_step: f64,
    /// Current-cell `X*ER(N,NN,N2,N1)` [boundary_side].
    local_flux_by_boundary_side: []const MineralFlux,
    /// Positive-neighbor `XSEDER(N,NN,N5,N4)` [boundary_side].
    positive_neighbor_total_sediment_Mg_per_step_by_boundary_side: []const f64,
};

pub const State = struct {
    net_incoming: MineralFlux,
};

/// Accounts local mineral erosion fluxes on one horizontal axis.
///
/// Traceability: REDIST.F lines 2888--2909 (`TSEDER`--`TAECER`). Runtime
/// boundary-side slices replace the source's fixed `NN=1,2`. A side is active
/// when either local or positive-neighbor total sediment magnitude is strictly
/// greater than the cell threshold; only the local six inventories are added
/// here. Vertical axes and disturbance modes without erosion bypass the block.
/// Candidate state commits only after all source-ordered additions are finite.
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
        try addFlux(&candidate, local);
    }
    state.net_incoming = candidate;
}

fn erosionEnabled(mode: DisturbanceMode) bool {
    return mode == .freeze_thaw_and_erosion or
        mode == .freeze_thaw_erosion_and_organic_matter;
}

fn validateInputs(inputs: Inputs, state: State) !void {
    if (inputs.local_flux_by_boundary_side.len == 0)
        return error.InvalidErosionIncomingDimensions;
    if (inputs.positive_neighbor_total_sediment_Mg_per_step_by_boundary_side.len !=
        inputs.local_flux_by_boundary_side.len)
    {
        return error.ErosionIncomingDimensionMismatch;
    }
    if (!std.math.isFinite(inputs.sediment_activity_threshold_Mg_per_step))
        return error.NonFiniteErosionIncomingInput;
    if (inputs.sediment_activity_threshold_Mg_per_step < 0)
        return error.InvalidErosionActivityThreshold;
    try validateFlux(state.net_incoming);
    for (inputs.local_flux_by_boundary_side) |flux| try validateFlux(flux);
    for (inputs.positive_neighbor_total_sediment_Mg_per_step_by_boundary_side) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteErosionIncomingInput;
}

fn validateFlux(flux: MineralFlux) !void {
    inline for (@typeInfo(MineralFlux).@"struct".fields) |field|
        if (!std.math.isFinite(@field(flux, field.name)))
            return error.NonFiniteErosionIncomingInput;
}

fn addFlux(candidate: *MineralFlux, contribution: MineralFlux) !void {
    inline for (@typeInfo(MineralFlux).@"struct".fields) |field| {
        const result = @field(candidate.*, field.name) +
            @field(contribution, field.name);
        if (!std.math.isFinite(result))
            return error.NonFiniteErosionIncomingResult;
        @field(candidate.*, field.name) = result;
    }
}

fn filledFlux(value: f64) MineralFlux {
    var result: MineralFlux = undefined;
    inline for (@typeInfo(MineralFlux).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn expectFlux(actual: MineralFlux, expected: f64) !void {
    inline for (@typeInfo(MineralFlux).@"struct".fields) |field|
        try std.testing.expectEqual(expected, @field(actual, field.name));
}

test "horizontal erosion adds source-ordered local mineral inventories" {
    const local = [_]MineralFlux{
        .{
            .total_sediment_Mg_per_step = 2,
            .sand_Mg_per_step = 3,
            .silt_Mg_per_step = 4,
            .clay_Mg_per_step = 5,
            .cation_exchange_capacity_mol_per_step = 6,
            .anion_exchange_capacity_mol_per_step = 7,
        },
        .{
            .total_sediment_Mg_per_step = 8,
            .sand_Mg_per_step = 9,
            .silt_Mg_per_step = 10,
            .clay_Mg_per_step = 11,
            .cation_exchange_capacity_mol_per_step = 12,
            .anion_exchange_capacity_mol_per_step = 13,
        },
    };
    const positive_neighbor_totals = [_]f64{ 0, 0 };
    var state = State{ .net_incoming = filledFlux(100) };

    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &positive_neighbor_totals,
    }, &state);

    try std.testing.expectEqual(
        @as(f64, 110),
        state.net_incoming.total_sediment_Mg_per_step,
    );
    try std.testing.expectEqual(@as(f64, 112), state.net_incoming.sand_Mg_per_step);
    try std.testing.expectEqual(@as(f64, 114), state.net_incoming.silt_Mg_per_step);
    try std.testing.expectEqual(@as(f64, 116), state.net_incoming.clay_Mg_per_step);
    try std.testing.expectEqual(
        @as(f64, 118),
        state.net_incoming.cation_exchange_capacity_mol_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 120),
        state.net_incoming.anion_exchange_capacity_mol_per_step,
    );
}

test "strict activity gate may be triggered by positive neighbor sediment" {
    const local = [_]MineralFlux{
        .{
            .total_sediment_Mg_per_step = 1,
            .sand_Mg_per_step = 2,
            .silt_Mg_per_step = 3,
            .clay_Mg_per_step = 4,
            .cation_exchange_capacity_mol_per_step = 5,
            .anion_exchange_capacity_mol_per_step = 6,
        },
        .{
            .total_sediment_Mg_per_step = 1,
            .sand_Mg_per_step = 7,
            .silt_Mg_per_step = 7,
            .clay_Mg_per_step = 7,
            .cation_exchange_capacity_mol_per_step = 7,
            .anion_exchange_capacity_mol_per_step = 7,
        },
    };
    const positive_neighbor_totals = [_]f64{ 2, 1 };
    var state = State{ .net_incoming = .{} };

    try account(.{
        .disturbance_mode = .freeze_thaw_erosion_and_organic_matter,
        .transport_axis = .north_south,
        .sediment_activity_threshold_Mg_per_step = 1,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &positive_neighbor_totals,
    }, &state);

    try std.testing.expectEqual(
        @as(f64, 1),
        state.net_incoming.total_sediment_Mg_per_step,
    );
    try std.testing.expectEqual(@as(f64, 2), state.net_incoming.sand_Mg_per_step);
    try std.testing.expectEqual(
        @as(f64, 6),
        state.net_incoming.anion_exchange_capacity_mol_per_step,
    );
}

test "accepted side sum exactly conserves every mineral inventory" {
    const local = [_]MineralFlux{
        filledFlux(2),
        filledFlux(-3),
        filledFlux(4),
        filledFlux(5),
    };
    const positive_neighbor_totals = [_]f64{ 0, 0, 0, 0 };
    var state = State{ .net_incoming = .{} };

    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &positive_neighbor_totals,
    }, &state);

    try expectFlux(state.net_incoming, 8);
}

test "non-erosion modes and vertical transport bypass unused inputs" {
    var state = State{ .net_incoming = filledFlux(9) };
    const empty_flux = [_]MineralFlux{};
    const empty_neighbor = [_]f64{};

    try account(.{
        .disturbance_mode = .freeze_thaw,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
        .local_flux_by_boundary_side = &empty_flux,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &empty_neighbor,
    }, &state);
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .vertical,
        .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
        .local_flux_by_boundary_side = &empty_flux,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &empty_neighbor,
    }, &state);
    try expectFlux(state.net_incoming, 9);
}

test "nonfinite input and late overflow preserve state atomically" {
    var local = [_]MineralFlux{ filledFlux(2), filledFlux(3) };
    local[1].anion_exchange_capacity_mol_per_step = std.math.nan(f64);
    const positive_neighbor_totals = [_]f64{ 0, 0 };
    var state = State{ .net_incoming = filledFlux(4) };
    const base_inputs = Inputs{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &positive_neighbor_totals,
    };

    try std.testing.expectError(
        error.NonFiniteErosionIncomingInput,
        account(base_inputs, &state),
    );
    try expectFlux(state.net_incoming, 4);

    local[1] = filledFlux(std.math.floatMax(f64));
    state.net_incoming = filledFlux(std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteErosionIncomingResult,
        account(base_inputs, &state),
    );
    try expectFlux(state.net_incoming, std.math.floatMax(f64));
}

test "dimension and threshold failures precede mutation" {
    const local = [_]MineralFlux{filledFlux(2)};
    const no_neighbors = [_]f64{};
    var state = State{ .net_incoming = filledFlux(5) };

    try std.testing.expectError(
        error.ErosionIncomingDimensionMismatch,
        account(.{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .north_south,
            .sediment_activity_threshold_Mg_per_step = 1,
            .local_flux_by_boundary_side = &local,
            .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &no_neighbors,
        }, &state),
    );
    try expectFlux(state.net_incoming, 5);

    const neighbor = [_]f64{0};
    try std.testing.expectError(
        error.InvalidErosionActivityThreshold,
        account(.{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .north_south,
            .sediment_activity_threshold_Mg_per_step = -1,
            .local_flux_by_boundary_side = &local,
            .positive_neighbor_total_sediment_Mg_per_step_by_boundary_side = &neighbor,
        }, &state),
    );
    try expectFlux(state.net_incoming, 5);
}
