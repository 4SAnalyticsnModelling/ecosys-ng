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

pub const BoundarySide = enum {
    first,
    second,
};

pub const MineralFlux = struct {
    total_sediment_megagrams_per_step: f64 = 0,
    sand_megagrams_per_step: f64 = 0,
    silt_megagrams_per_step: f64 = 0,
    clay_megagrams_per_step: f64 = 0,
    cation_exchange_capacity_mol_per_step: f64 = 0,
    anion_exchange_capacity_mol_per_step: f64 = 0,
};

pub const Inputs = struct {
    disturbance_mode: DisturbanceMode,
    transport_axis: TransportAxis,
    boundary_side: BoundarySide,
    sediment_activity_threshold_megagrams_per_step: f64,
    /// Null when the geometry-derived opposite-neighbor coordinate is absent.
    opposite_neighbor_first_side_flux: ?MineralFlux,
};

pub const State = struct {
    net_erosion: MineralFlux,
};

/// Subtracts mineral erosion from the first side's opposite neighbor.
///
/// Traceability: REDIST.F lines 3192--3206. The optional flux represents valid
/// runtime `N4B,N5B` geometry; `.first` represents source `NN == 1`. Unlike the
/// preceding positive-neighbor branch, this branch has no `IFLBH` connection
/// gate. Its single sediment magnitude must be strictly above the cell
/// threshold. All six inventories commit only after finite evaluation.
pub fn account(inputs: Inputs, state: *State) !void {
    if (!erosionEnabled(inputs.disturbance_mode) or
        inputs.transport_axis == .vertical or
        inputs.boundary_side != .first)
    {
        return;
    }
    const flux = inputs.opposite_neighbor_first_side_flux orelse return;
    try validateInputs(inputs.sediment_activity_threshold_megagrams_per_step, flux, state.*);
    if (@abs(flux.total_sediment_megagrams_per_step) <=
        inputs.sediment_activity_threshold_megagrams_per_step)
    {
        return;
    }

    var candidate = state.net_erosion;
    inline for (@typeInfo(MineralFlux).@"struct".fields) |field| {
        const result = @field(candidate, field.name) - @field(flux, field.name);
        if (!std.math.isFinite(result))
            return error.NonFiniteOppositeNeighborMineralErosionResult;
        @field(candidate, field.name) = result;
    }
    state.net_erosion = candidate;
}

fn erosionEnabled(mode: DisturbanceMode) bool {
    return mode == .freeze_thaw_and_erosion or
        mode == .freeze_thaw_erosion_and_organic_matter;
}

fn validateInputs(threshold: f64, flux: MineralFlux, state: State) !void {
    if (!std.math.isFinite(threshold))
        return error.NonFiniteOppositeNeighborMineralErosionInput;
    if (threshold < 0) return error.InvalidErosionActivityThreshold;
    try validateFlux(flux);
    try validateFlux(state.net_erosion);
}

fn validateFlux(flux: MineralFlux) !void {
    inline for (@typeInfo(MineralFlux).@"struct".fields) |field|
        if (!std.math.isFinite(@field(flux, field.name)))
            return error.NonFiniteOppositeNeighborMineralErosionInput;
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

test "active first-side opposite neighbor subtracts six mineral inventories" {
    var state = State{ .net_erosion = filledFlux(10) };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .opposite_neighbor_first_side_flux = filledFlux(3),
    }, &state);
    try expectFlux(state.net_erosion, 7);
}

test "sediment activity threshold is strict" {
    var state = State{ .net_erosion = filledFlux(10) };
    var flux = filledFlux(4);
    flux.total_sediment_megagrams_per_step = -1;
    try account(.{
        .disturbance_mode = .freeze_thaw_erosion_and_organic_matter,
        .transport_axis = .north_south,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .opposite_neighbor_first_side_flux = flux,
    }, &state);
    try expectFlux(state.net_erosion, 10);
    flux.total_sediment_megagrams_per_step = -1.0001;
    try account(.{
        .disturbance_mode = .freeze_thaw_erosion_and_organic_matter,
        .transport_axis = .north_south,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .opposite_neighbor_first_side_flux = flux,
    }, &state);
    try std.testing.expectEqual(
        @as(f64, 11.0001),
        state.net_erosion.total_sediment_megagrams_per_step,
    );
    try std.testing.expectEqual(@as(f64, 6), state.net_erosion.sand_megagrams_per_step);
}

test "shared face opposite-neighbor transfer conserves exactly" {
    const shared = filledFlux(7);
    const receiving_cell = shared;
    var source_cell = State{ .net_erosion = .{} };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .opposite_neighbor_first_side_flux = shared,
    }, &source_cell);
    inline for (@typeInfo(MineralFlux).@"struct".fields) |field|
        try std.testing.expectEqual(
            @as(f64, 0),
            @field(receiving_cell, field.name) +
                @field(source_cell.net_erosion, field.name),
        );
}

test "geometry side process and axis gates bypass unused input" {
    var state = State{ .net_erosion = filledFlux(9) };
    const bypass = [_]Inputs{
        .{
            .disturbance_mode = .freeze_thaw,
            .transport_axis = .east_west,
            .boundary_side = .first,
            .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
            .opposite_neighbor_first_side_flux = null,
        },
        .{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .vertical,
            .boundary_side = .first,
            .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
            .opposite_neighbor_first_side_flux = null,
        },
        .{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .east_west,
            .boundary_side = .second,
            .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
            .opposite_neighbor_first_side_flux = null,
        },
        .{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .east_west,
            .boundary_side = .first,
            .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
            .opposite_neighbor_first_side_flux = null,
        },
    };
    for (bypass) |inputs| try account(inputs, &state);
    try expectFlux(state.net_erosion, 9);
}

test "invalid input and subtraction overflow preserve state atomically" {
    var flux = filledFlux(3);
    flux.anion_exchange_capacity_mol_per_step = std.math.nan(f64);
    var state = State{ .net_erosion = filledFlux(5) };
    var inputs = Inputs{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .opposite_neighbor_first_side_flux = flux,
    };
    try std.testing.expectError(
        error.NonFiniteOppositeNeighborMineralErosionInput,
        account(inputs, &state),
    );
    try expectFlux(state.net_erosion, 5);

    inputs.opposite_neighbor_first_side_flux =
        filledFlux(-std.math.floatMax(f64));
    state.net_erosion = filledFlux(std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteOppositeNeighborMineralErosionResult,
        account(inputs, &state),
    );
    try expectFlux(state.net_erosion, std.math.floatMax(f64));
}
