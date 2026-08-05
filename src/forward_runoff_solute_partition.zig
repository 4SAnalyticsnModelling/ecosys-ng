const std = @import("std");

pub const runoff_species_count = 42;

pub const BoundaryAxis = enum {
    east_west,
    north_south,
};

pub const Inputs = struct {
    /// Retained explicitly for the caller's runtime directional axis (`N`).
    axis: BoundaryAxis,
    source_runoff_water_m3_per_step: f64,
    forward_boundary_water_m3_per_step: f64,
    minimum_runoff_water_m3_per_step: f64,
    source_runoff_solute_flux_amount_per_step: []const f64,
};

/// Exact source-order forward-side (`NN=1`) translation of TRNSFRS.F lines
/// 4070--4206. Instantaneous flux is overwritten and its active contribution
/// is added to the REDIST accumulator. Most fields are mol step-1; neutral
/// amount naming retains the documented free-phosphate unit ambiguity.
pub fn partitionAndAccumulate(
    inputs: Inputs,
    boundary_flux_amount_per_step: []f64,
    accumulated_boundary_flux_amount_per_step: []f64,
) !void {
    _ = inputs.axis;
    if (inputs.source_runoff_solute_flux_amount_per_step.len != runoff_species_count or
        boundary_flux_amount_per_step.len != runoff_species_count or
        accumulated_boundary_flux_amount_per_step.len != runoff_species_count)
        return error.ForwardRunoffSolutePartitionDimensionMismatch;
    inline for (.{ inputs.source_runoff_water_m3_per_step, inputs.forward_boundary_water_m3_per_step, inputs.minimum_runoff_water_m3_per_step }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteForwardRunoffSolutePartitionInput;
    if (inputs.source_runoff_water_m3_per_step < 0 or inputs.forward_boundary_water_m3_per_step < 0 or
        inputs.minimum_runoff_water_m3_per_step < 0)
        return error.InvalidForwardRunoffSolutePartitionInput;
    for (inputs.source_runoff_solute_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteForwardRunoffSolutePartitionInput;
    for (accumulated_boundary_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteForwardRunoffSolutePartitionInput;

    if (!(inputs.source_runoff_water_m3_per_step > inputs.minimum_runoff_water_m3_per_step)) {
        @memset(boundary_flux_amount_per_step, 0);
        return;
    }
    const directional_fraction = inputs.forward_boundary_water_m3_per_step /
        inputs.source_runoff_water_m3_per_step;
    for (0..runoff_species_count) |species| {
        const flux = inputs.source_runoff_solute_flux_amount_per_step[species] * directional_fraction;
        if (!std.math.isFinite(flux) or
            !std.math.isFinite(accumulated_boundary_flux_amount_per_step[species] + flux))
            return error.NonFiniteForwardRunoffSolutePartitionResult;
    }
    for (0..runoff_species_count) |species| {
        const flux = inputs.source_runoff_solute_flux_amount_per_step[species] * directional_fraction;
        boundary_flux_amount_per_step[species] = flux;
        accumulated_boundary_flux_amount_per_step[species] += flux;
    }
}

fn fixture(source_flux: []const f64) Inputs {
    return .{ .axis = .east_west, .source_runoff_water_m3_per_step = 4, .forward_boundary_water_m3_per_step = 1, .minimum_runoff_water_m3_per_step = 0.01, .source_runoff_solute_flux_amount_per_step = source_flux };
}

test "TRNSFRS forward side partitions and accumulates all 42 fields" {
    var source: [runoff_species_count]f64 = undefined;
    for (&source, 0..) |*value, species| value.* = @floatFromInt(species + 1);
    var boundary = [_]f64{9} ** runoff_species_count;
    var accumulated = [_]f64{10} ** runoff_species_count;
    try partitionAndAccumulate(fixture(&source), &boundary, &accumulated);
    try std.testing.expectEqual(@as(f64, 0.25), boundary[0]);
    try std.testing.expectEqual(@as(f64, 10.25), accumulated[0]);
    try std.testing.expectEqual(@as(f64, 10.5), boundary[runoff_species_count - 1]);
    try std.testing.expectEqual(@as(f64, 20.5), accumulated[runoff_species_count - 1]);
}

test "inactive runoff zeros instantaneous flux without changing accumulation" {
    const source = [_]f64{8} ** runoff_species_count;
    var boundary = [_]f64{9} ** runoff_species_count;
    var accumulated = [_]f64{10} ** runoff_species_count;
    var inputs = fixture(&source);
    inputs.source_runoff_water_m3_per_step = 0;
    try partitionAndAccumulate(inputs, &boundary, &accumulated);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** runoff_species_count), &boundary);
    try std.testing.expectEqualSlices(f64, &([_]f64{10} ** runoff_species_count), &accumulated);
}

test "late invalid accumulator leaves instantaneous and accumulated flux atomic" {
    const source = [_]f64{8} ** runoff_species_count;
    var boundary = [_]f64{9} ** runoff_species_count;
    var accumulated = [_]f64{10} ** runoff_species_count;
    accumulated[runoff_species_count - 1] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteForwardRunoffSolutePartitionInput, partitionAndAccumulate(fixture(&source), &boundary, &accumulated));
    try std.testing.expectEqual(@as(f64, 9), boundary[0]);
    try std.testing.expectEqual(@as(f64, 10), accumulated[0]);
}

test "forward partition rejects phantom runoff species" {
    const source = [_]f64{8} ** (runoff_species_count + 1);
    var boundary = [_]f64{9} ** runoff_species_count;
    var accumulated = [_]f64{10} ** runoff_species_count;
    try std.testing.expectError(error.ForwardRunoffSolutePartitionDimensionMismatch, partitionAndAccumulate(fixture(&source), &boundary, &accumulated));
}
