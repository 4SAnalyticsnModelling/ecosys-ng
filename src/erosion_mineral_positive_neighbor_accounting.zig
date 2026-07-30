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

pub const NeighborConnection = enum {
    blocked,
    connected,
};

pub const MineralFlux = struct {
    total_sediment_Mg_per_step: f64 = 0,
    sand_Mg_per_step: f64 = 0,
    silt_Mg_per_step: f64 = 0,
    clay_Mg_per_step: f64 = 0,
    cation_exchange_capacity_mol_per_step: f64 = 0,
    anion_exchange_capacity_mol_per_step: f64 = 0,
};

pub const NeighborSide = struct {
    local_total_sediment_Mg_per_step: f64,
    positive_neighbor_flux: MineralFlux,
    connection: NeighborConnection,
};

pub const Inputs = struct {
    disturbance_mode: DisturbanceMode,
    transport_axis: TransportAxis,
    sediment_activity_threshold_Mg_per_step: f64,
    neighbor_by_boundary_side: []const NeighborSide,
};

pub const State = struct {
    net_erosion: MineralFlux,
};

/// Subtracts connected positive-neighbor mineral erosion on one horizontal axis.
///
/// Traceability: REDIST.F lines 3051--3057, within gates 2888 and 2894--2895.
/// Runtime sides replace fixed `NN=1,2`. A side first requires local or
/// positive-neighbor sediment magnitude strictly above the threshold, then the
/// `IFLBH == 0` connection. The six positive-neighbor inventories are
/// subtracted in source order and committed only after finite evaluation.
pub fn account(inputs: Inputs, state: *State) !void {
    if (!erosionEnabled(inputs.disturbance_mode) or
        inputs.transport_axis == .vertical)
    {
        return;
    }
    try validateInputs(inputs, state.*);
    var candidate = state.net_erosion;
    for (inputs.neighbor_by_boundary_side) |side| {
        if (@abs(side.local_total_sediment_Mg_per_step) <=
            inputs.sediment_activity_threshold_Mg_per_step and
            @abs(side.positive_neighbor_flux.total_sediment_Mg_per_step) <=
                inputs.sediment_activity_threshold_Mg_per_step)
        {
            continue;
        }
        if (side.connection == .connected)
            try subtractFlux(&candidate, side.positive_neighbor_flux);
    }
    state.net_erosion = candidate;
}

fn erosionEnabled(mode: DisturbanceMode) bool {
    return mode == .freeze_thaw_and_erosion or
        mode == .freeze_thaw_erosion_and_organic_matter;
}

fn validateInputs(inputs: Inputs, state: State) !void {
    if (inputs.neighbor_by_boundary_side.len == 0)
        return error.InvalidMineralNeighborErosionDimensions;
    if (!std.math.isFinite(inputs.sediment_activity_threshold_Mg_per_step))
        return error.NonFiniteMineralNeighborErosionInput;
    if (inputs.sediment_activity_threshold_Mg_per_step < 0)
        return error.InvalidErosionActivityThreshold;
    try validateFlux(state.net_erosion);
    for (inputs.neighbor_by_boundary_side) |side| {
        if (!std.math.isFinite(side.local_total_sediment_Mg_per_step))
            return error.NonFiniteMineralNeighborErosionInput;
        try validateFlux(side.positive_neighbor_flux);
    }
}

fn validateFlux(flux: MineralFlux) !void {
    inline for (@typeInfo(MineralFlux).@"struct".fields) |field|
        if (!std.math.isFinite(@field(flux, field.name)))
            return error.NonFiniteMineralNeighborErosionInput;
}

fn subtractFlux(candidate: *MineralFlux, contribution: MineralFlux) !void {
    inline for (@typeInfo(MineralFlux).@"struct".fields) |field| {
        const result = @field(candidate.*, field.name) -
            @field(contribution, field.name);
        if (!std.math.isFinite(result))
            return error.NonFiniteMineralNeighborErosionResult;
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

test "connected positive neighbor subtracts six mineral inventories" {
    const sides = [_]NeighborSide{
        .{
            .local_total_sediment_Mg_per_step = 2,
            .positive_neighbor_flux = filledFlux(3),
            .connection = .connected,
        },
        .{
            .local_total_sediment_Mg_per_step = 4,
            .positive_neighbor_flux = filledFlux(5),
            .connection = .connected,
        },
    };
    var state = State{ .net_erosion = filledFlux(100) };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    }, &state);
    try expectFlux(state.net_erosion, 92);
}

test "connection and strict activity gates are independent" {
    const sides = [_]NeighborSide{
        .{
            .local_total_sediment_Mg_per_step = 2,
            .positive_neighbor_flux = filledFlux(7),
            .connection = .blocked,
        },
        .{
            .local_total_sediment_Mg_per_step = 1,
            .positive_neighbor_flux = .{
                .total_sediment_Mg_per_step = 1,
                .sand_Mg_per_step = 9,
                .silt_Mg_per_step = 9,
                .clay_Mg_per_step = 9,
                .cation_exchange_capacity_mol_per_step = 9,
                .anion_exchange_capacity_mol_per_step = 9,
            },
            .connection = .connected,
        },
        .{
            .local_total_sediment_Mg_per_step = 2,
            .positive_neighbor_flux = .{
                .total_sediment_Mg_per_step = 1,
                .sand_Mg_per_step = 4,
                .silt_Mg_per_step = 4,
                .clay_Mg_per_step = 4,
                .cation_exchange_capacity_mol_per_step = 4,
                .anion_exchange_capacity_mol_per_step = 4,
            },
            .connection = .connected,
        },
    };
    var state = State{ .net_erosion = filledFlux(10) };
    try account(.{
        .disturbance_mode = .freeze_thaw_erosion_and_organic_matter,
        .transport_axis = .north_south,
        .sediment_activity_threshold_Mg_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    }, &state);
    try std.testing.expectEqual(
        @as(f64, 9),
        state.net_erosion.total_sediment_Mg_per_step,
    );
    try std.testing.expectEqual(@as(f64, 6), state.net_erosion.sand_Mg_per_step);
}

test "shared face local addition and neighbor subtraction conserve exactly" {
    const shared_face_flux = filledFlux(7);
    const receiving_cell = shared_face_flux;
    var source_cell = State{ .net_erosion = .{} };
    const sides = [_]NeighborSide{.{
        .local_total_sediment_Mg_per_step = 0,
        .positive_neighbor_flux = shared_face_flux,
        .connection = .connected,
    }};
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    }, &source_cell);
    inline for (@typeInfo(MineralFlux).@"struct".fields) |field|
        try std.testing.expectEqual(
            @as(f64, 0),
            @field(receiving_cell, field.name) +
                @field(source_cell.net_erosion, field.name),
        );
}

test "disabled erosion and vertical axes bypass unused neighbors" {
    var state = State{ .net_erosion = filledFlux(9) };
    try account(.{
        .disturbance_mode = .freeze_thaw,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
        .neighbor_by_boundary_side = &.{},
    }, &state);
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .vertical,
        .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
        .neighbor_by_boundary_side = &.{},
    }, &state);
    try expectFlux(state.net_erosion, 9);
}

test "invalid input and subtraction overflow preserve state atomically" {
    var sides = [_]NeighborSide{.{
        .local_total_sediment_Mg_per_step = 2,
        .positive_neighbor_flux = filledFlux(3),
        .connection = .connected,
    }};
    sides[0].positive_neighbor_flux.anion_exchange_capacity_mol_per_step =
        std.math.nan(f64);
    var state = State{ .net_erosion = filledFlux(5) };
    const inputs = Inputs{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    };
    try std.testing.expectError(
        error.NonFiniteMineralNeighborErosionInput,
        account(inputs, &state),
    );
    try expectFlux(state.net_erosion, 5);

    sides[0].positive_neighbor_flux = filledFlux(-std.math.floatMax(f64));
    state.net_erosion = filledFlux(std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteMineralNeighborErosionResult,
        account(inputs, &state),
    );
    try expectFlux(state.net_erosion, std.math.floatMax(f64));
}
