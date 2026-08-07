const std = @import("std");

pub const runoff_species_count = 42;

pub const BoundaryAxis = enum {
    east_west,
    north_south,
};

pub const Inputs = struct {
    axis: BoundaryAxis,
    source_runoff_water_m3_per_step: f64,
    backward_boundary_water_m3_per_step: f64,
    minimum_runoff_water_m3_per_step: f64,
    source_runoff_solute_flux_amount_per_step: []const f64,
};

/// Exact source-order backward-side (`NN=2`) translation of TRNSFRS.F lines
/// 4229--4366 for a valid internal neighbor returned by the boundary locator.
/// Unlike the forward branch, inactive source runoff leaves backward storage
/// untouched because the outer legacy `ELSE` only zeros side 2. Neutral amount
/// naming retains the documented free-phosphate unit ambiguity.
pub fn partitionAndAccumulate(
    inputs: Inputs,
    boundary_flux_amount_per_step: []f64,
    accumulated_boundary_flux_amount_per_step: []f64,
) !void {
    _ = inputs.axis;
    if (inputs.source_runoff_solute_flux_amount_per_step.len != runoff_species_count or
        boundary_flux_amount_per_step.len != runoff_species_count or
        accumulated_boundary_flux_amount_per_step.len != runoff_species_count)
        return error.BackwardRunoffSolutePartitionDimensionMismatch;
    inline for (.{ inputs.source_runoff_water_m3_per_step, inputs.backward_boundary_water_m3_per_step, inputs.minimum_runoff_water_m3_per_step }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteBackwardRunoffSolutePartitionInput;
    if (inputs.source_runoff_water_m3_per_step < 0 or inputs.backward_boundary_water_m3_per_step < 0 or
        inputs.minimum_runoff_water_m3_per_step < 0)
        return error.InvalidBackwardRunoffSolutePartitionInput;
    for (inputs.source_runoff_solute_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteBackwardRunoffSolutePartitionInput;
    for (accumulated_boundary_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteBackwardRunoffSolutePartitionInput;

    if (!(inputs.source_runoff_water_m3_per_step > inputs.minimum_runoff_water_m3_per_step)) return;
    const directional_fraction = inputs.backward_boundary_water_m3_per_step /
        inputs.source_runoff_water_m3_per_step;
    for (0..runoff_species_count) |species| {
        const flux = inputs.source_runoff_solute_flux_amount_per_step[species] * directional_fraction;
        if (!std.math.isFinite(flux) or
            !std.math.isFinite(accumulated_boundary_flux_amount_per_step[species] + flux))
            return error.NonFiniteBackwardRunoffSolutePartitionResult;
    }
    for (0..runoff_species_count) |species| {
        const flux = inputs.source_runoff_solute_flux_amount_per_step[species] * directional_fraction;
        boundary_flux_amount_per_step[species] = flux;
        accumulated_boundary_flux_amount_per_step[species] += flux;
    }
}

fn fixture(source_flux: []const f64) Inputs {
    return .{ .axis = .north_south, .source_runoff_water_m3_per_step = 5, .backward_boundary_water_m3_per_step = 2, .minimum_runoff_water_m3_per_step = 0.01, .source_runoff_solute_flux_amount_per_step = source_flux };
}

test "TRNSFRS backward side partitions and accumulates all 42 fields" {
    var source: [runoff_species_count]f64 = undefined;
    for (&source, 0..) |*value, species| value.* = @floatFromInt(species + 1);
    var boundary = [_]f64{9} ** runoff_species_count;
    var accumulated = [_]f64{10} ** runoff_species_count;
    try partitionAndAccumulate(fixture(&source), &boundary, &accumulated);
    try std.testing.expectEqual(@as(f64, 0.4), boundary[0]);
    try std.testing.expectEqual(@as(f64, 10.4), accumulated[0]);
    try std.testing.expectEqual(@as(f64, 16.8), boundary[runoff_species_count - 1]);
    try std.testing.expectEqual(@as(f64, 26.8), accumulated[runoff_species_count - 1]);
}

test "inactive runoff preserves backward instantaneous and accumulated storage" {
    const source = [_]f64{8} ** runoff_species_count;
    var boundary = [_]f64{9} ** runoff_species_count;
    var accumulated = [_]f64{10} ** runoff_species_count;
    var inputs = fixture(&source);
    inputs.source_runoff_water_m3_per_step = 0;
    try partitionAndAccumulate(inputs, &boundary, &accumulated);
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** runoff_species_count), &boundary);
    try std.testing.expectEqualSlices(f64, &([_]f64{10} ** runoff_species_count), &accumulated);
}

test "late invalid accumulator leaves both backward arrays atomic" {
    const source = [_]f64{8} ** runoff_species_count;
    var boundary = [_]f64{9} ** runoff_species_count;
    var accumulated = [_]f64{10} ** runoff_species_count;
    accumulated[runoff_species_count - 1] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteBackwardRunoffSolutePartitionInput, partitionAndAccumulate(fixture(&source), &boundary, &accumulated));
    try std.testing.expectEqual(@as(f64, 9), boundary[0]);
    try std.testing.expectEqual(@as(f64, 10), accumulated[0]);
}

test "backward partition rejects a short source topology" {
    const source = [_]f64{8} ** (runoff_species_count - 1);
    var boundary = [_]f64{9} ** runoff_species_count;
    var accumulated = [_]f64{10} ** runoff_species_count;
    try std.testing.expectError(error.BackwardRunoffSolutePartitionDimensionMismatch, partitionAndAccumulate(fixture(&source), &boundary, &accumulated));
}
