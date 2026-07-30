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

pub const OppositeNeighborFlux = struct {
    total_sediment_Mg_per_step: f64,
    fertilizer: FertilizerFlux,
};

pub const Inputs = struct {
    disturbance_mode: DisturbanceMode,
    transport_axis: TransportAxis,
    boundary_side: BoundarySide,
    sediment_activity_threshold_Mg_per_step: f64,
    /// Null when the geometry-derived opposite-neighbor coordinate is absent.
    opposite_neighbor_first_side_flux: ?OppositeNeighborFlux,
};

pub const State = struct {
    net_erosion: FertilizerFlux,
};

/// Subtracts fertilizer erosion from the first side's opposite neighbor.
///
/// Traceability: REDIST.F lines 3218--3225 under gates 3192--3193. Optional
/// geometry replaces `N4B,N5B`; `.first` represents `NN == 1`. The opposite
/// sediment magnitude alone activates all eight non-band/band molar pools.
/// Candidate state commits only after every source-ordered subtraction is
/// finite.
pub fn account(inputs: Inputs, state: *State) !void {
    if (!erosionEnabled(inputs.disturbance_mode) or
        inputs.transport_axis == .vertical or
        inputs.boundary_side != .first)
    {
        return;
    }
    const flux = inputs.opposite_neighbor_first_side_flux orelse return;
    try validateInputs(
        inputs.sediment_activity_threshold_Mg_per_step,
        flux,
        state.*,
    );
    if (@abs(flux.total_sediment_Mg_per_step) <=
        inputs.sediment_activity_threshold_Mg_per_step)
    {
        return;
    }

    var candidate = state.net_erosion;
    inline for (@typeInfo(FertilizerFlux).@"struct".fields) |field| {
        const result = @field(candidate, field.name) -
            @field(flux.fertilizer, field.name);
        if (!std.math.isFinite(result))
            return error.NonFiniteOppositeNeighborFertilizerErosionResult;
        @field(candidate, field.name) = result;
    }
    state.net_erosion = candidate;
}

fn erosionEnabled(mode: DisturbanceMode) bool {
    return mode == .freeze_thaw_and_erosion or
        mode == .freeze_thaw_erosion_and_organic_matter;
}

fn validateInputs(threshold: f64, flux: OppositeNeighborFlux, state: State) !void {
    if (!std.math.isFinite(threshold) or
        !std.math.isFinite(flux.total_sediment_Mg_per_step))
    {
        return error.NonFiniteOppositeNeighborFertilizerErosionInput;
    }
    if (threshold < 0) return error.InvalidErosionActivityThreshold;
    try validateFertilizer(flux.fertilizer);
    try validateFertilizer(state.net_erosion);
}

fn validateFertilizer(flux: FertilizerFlux) !void {
    inline for (@typeInfo(FertilizerFlux).@"struct".fields) |field|
        if (!std.math.isFinite(@field(flux, field.name)))
            return error.NonFiniteOppositeNeighborFertilizerErosionInput;
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

test "active opposite neighbor subtracts all eight fertilizer pools" {
    var state = State{ .net_erosion = filledFlux(10) };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .boundary_side = .first,
        .sediment_activity_threshold_Mg_per_step = 1,
        .opposite_neighbor_first_side_flux = .{
            .total_sediment_Mg_per_step = -2,
            .fertilizer = filledFlux(3),
        },
    }, &state);
    try expectFlux(state.net_erosion, 7);
}

test "opposite sediment threshold is strict and independent of pool values" {
    var state = State{ .net_erosion = filledFlux(10) };
    const flux = OppositeNeighborFlux{
        .total_sediment_Mg_per_step = 1,
        .fertilizer = filledFlux(100),
    };
    try account(.{
        .disturbance_mode = .freeze_thaw_erosion_and_organic_matter,
        .transport_axis = .north_south,
        .boundary_side = .first,
        .sediment_activity_threshold_Mg_per_step = 1,
        .opposite_neighbor_first_side_flux = flux,
    }, &state);
    try expectFlux(state.net_erosion, 10);
}

test "shared face fertilizer transfer conserves all pools exactly" {
    const shared = filledFlux(7);
    const receiving_cell = shared;
    var source_cell = State{ .net_erosion = .{} };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .boundary_side = .first,
        .sediment_activity_threshold_Mg_per_step = 1,
        .opposite_neighbor_first_side_flux = .{
            .total_sediment_Mg_per_step = 7,
            .fertilizer = shared,
        },
    }, &source_cell);
    inline for (@typeInfo(FertilizerFlux).@"struct".fields) |field|
        try std.testing.expectEqual(
            @as(f64, 0),
            @field(receiving_cell, field.name) +
                @field(source_cell.net_erosion, field.name),
        );
}

test "outer and geometry gates bypass unused fertilizer input" {
    var state = State{ .net_erosion = filledFlux(9) };
    const bypass = [_]Inputs{
        .{
            .disturbance_mode = .freeze_thaw,
            .transport_axis = .east_west,
            .boundary_side = .first,
            .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
            .opposite_neighbor_first_side_flux = null,
        },
        .{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .vertical,
            .boundary_side = .first,
            .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
            .opposite_neighbor_first_side_flux = null,
        },
        .{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .east_west,
            .boundary_side = .second,
            .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
            .opposite_neighbor_first_side_flux = null,
        },
        .{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .east_west,
            .boundary_side = .first,
            .sediment_activity_threshold_Mg_per_step = std.math.nan(f64),
            .opposite_neighbor_first_side_flux = null,
        },
    };
    for (bypass) |inputs| try account(inputs, &state);
    try expectFlux(state.net_erosion, 9);
}

test "invalid and overflow failures preserve fertilizer state atomically" {
    var flux = filledFlux(3);
    flux.band_nitrate_mol_per_step = std.math.nan(f64);
    var state = State{ .net_erosion = filledFlux(5) };
    var inputs = Inputs{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .boundary_side = .first,
        .sediment_activity_threshold_Mg_per_step = 1,
        .opposite_neighbor_first_side_flux = .{
            .total_sediment_Mg_per_step = 2,
            .fertilizer = flux,
        },
    };
    try std.testing.expectError(
        error.NonFiniteOppositeNeighborFertilizerErosionInput,
        account(inputs, &state),
    );
    try expectFlux(state.net_erosion, 5);

    inputs.opposite_neighbor_first_side_flux = .{
        .total_sediment_Mg_per_step = 2,
        .fertilizer = filledFlux(-std.math.floatMax(f64)),
    };
    state.net_erosion = filledFlux(std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteOppositeNeighborFertilizerErosionResult,
        account(inputs, &state),
    );
    try expectFlux(state.net_erosion, std.math.floatMax(f64));
}
