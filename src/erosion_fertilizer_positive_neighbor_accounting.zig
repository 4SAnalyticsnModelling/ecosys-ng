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

pub const NeighborSide = struct {
    local_total_sediment_Mg_per_step: f64,
    positive_neighbor_total_sediment_Mg_per_step: f64,
    positive_neighbor_fertilizer: FertilizerFlux,
    connection: NeighborConnection,
};

pub const Inputs = struct {
    disturbance_mode: DisturbanceMode,
    transport_axis: TransportAxis,
    sediment_activity_threshold_Mg_per_step: f64,
    neighbor_by_boundary_side: []const NeighborSide,
};

pub const State = struct {
    net_erosion: FertilizerFlux,
};

/// Subtracts connected positive-neighbor fertilizer-pool erosion.
///
/// Traceability: REDIST.F lines 3070--3077 under gates 2888, 2894--2895, and
/// 3051. Runtime sides replace fixed `NN=1,2`. After the strict sediment gate
/// and connection gate, eight positive-neighbor non-band/band molar pools are
/// subtracted in source order. Candidate state commits atomically.
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
            @abs(side.positive_neighbor_total_sediment_Mg_per_step) <=
                inputs.sediment_activity_threshold_Mg_per_step)
        {
            continue;
        }
        if (side.connection == .connected)
            try subtractFlux(&candidate, side.positive_neighbor_fertilizer);
    }
    state.net_erosion = candidate;
}

fn erosionEnabled(mode: DisturbanceMode) bool {
    return mode == .freeze_thaw_and_erosion or
        mode == .freeze_thaw_erosion_and_organic_matter;
}

fn validateInputs(inputs: Inputs, state: State) !void {
    if (inputs.neighbor_by_boundary_side.len == 0)
        return error.InvalidFertilizerNeighborErosionDimensions;
    if (!std.math.isFinite(inputs.sediment_activity_threshold_Mg_per_step))
        return error.NonFiniteFertilizerNeighborErosionInput;
    if (inputs.sediment_activity_threshold_Mg_per_step < 0)
        return error.InvalidErosionActivityThreshold;
    try validateFlux(state.net_erosion);
    for (inputs.neighbor_by_boundary_side) |side| {
        if (!std.math.isFinite(side.local_total_sediment_Mg_per_step) or
            !std.math.isFinite(side.positive_neighbor_total_sediment_Mg_per_step))
        {
            return error.NonFiniteFertilizerNeighborErosionInput;
        }
        try validateFlux(side.positive_neighbor_fertilizer);
    }
}

fn validateFlux(flux: FertilizerFlux) !void {
    inline for (@typeInfo(FertilizerFlux).@"struct".fields) |field|
        if (!std.math.isFinite(@field(flux, field.name)))
            return error.NonFiniteFertilizerNeighborErosionInput;
}

fn subtractFlux(candidate: *FertilizerFlux, contribution: FertilizerFlux) !void {
    inline for (@typeInfo(FertilizerFlux).@"struct".fields) |field| {
        const result = @field(candidate.*, field.name) -
            @field(contribution, field.name);
        if (!std.math.isFinite(result))
            return error.NonFiniteFertilizerNeighborErosionResult;
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

test "connected active sides subtract all eight fertilizer pools" {
    const sides = [_]NeighborSide{
        .{
            .local_total_sediment_Mg_per_step = 2,
            .positive_neighbor_total_sediment_Mg_per_step = 3,
            .positive_neighbor_fertilizer = filledFlux(3),
            .connection = .connected,
        },
        .{
            .local_total_sediment_Mg_per_step = 4,
            .positive_neighbor_total_sediment_Mg_per_step = 5,
            .positive_neighbor_fertilizer = filledFlux(5),
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

test "strict activity and connection gates remain independent" {
    const sides = [_]NeighborSide{
        .{
            .local_total_sediment_Mg_per_step = 2,
            .positive_neighbor_total_sediment_Mg_per_step = 7,
            .positive_neighbor_fertilizer = filledFlux(7),
            .connection = .blocked,
        },
        .{
            .local_total_sediment_Mg_per_step = 1,
            .positive_neighbor_total_sediment_Mg_per_step = 1,
            .positive_neighbor_fertilizer = filledFlux(9),
            .connection = .connected,
        },
        .{
            .local_total_sediment_Mg_per_step = 2,
            .positive_neighbor_total_sediment_Mg_per_step = 1,
            .positive_neighbor_fertilizer = filledFlux(4),
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
    try expectFlux(state.net_erosion, 6);
}

test "shared face fertilizer transfer conserves every pool exactly" {
    const shared_flux = filledFlux(7);
    const receiving_cell = shared_flux;
    var source_cell = State{ .net_erosion = .{} };
    const sides = [_]NeighborSide{.{
        .local_total_sediment_Mg_per_step = 0,
        .positive_neighbor_total_sediment_Mg_per_step = 7,
        .positive_neighbor_fertilizer = shared_flux,
        .connection = .connected,
    }};
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    }, &source_cell);
    inline for (@typeInfo(FertilizerFlux).@"struct".fields) |field|
        try std.testing.expectEqual(
            @as(f64, 0),
            @field(receiving_cell, field.name) +
                @field(source_cell.net_erosion, field.name),
        );
}

test "disabled and vertical modes bypass unused neighbor input" {
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
        .positive_neighbor_total_sediment_Mg_per_step = 3,
        .positive_neighbor_fertilizer = filledFlux(3),
        .connection = .connected,
    }};
    sides[0].positive_neighbor_fertilizer.band_nitrate_mol_per_step =
        std.math.nan(f64);
    var state = State{ .net_erosion = filledFlux(5) };
    const inputs = Inputs{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_Mg_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    };
    try std.testing.expectError(
        error.NonFiniteFertilizerNeighborErosionInput,
        account(inputs, &state),
    );
    try expectFlux(state.net_erosion, 5);

    sides[0].positive_neighbor_fertilizer =
        filledFlux(-std.math.floatMax(f64));
    state.net_erosion = filledFlux(std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteFertilizerNeighborErosionResult,
        account(inputs, &state),
    );
    try expectFlux(state.net_erosion, std.math.floatMax(f64));
}
