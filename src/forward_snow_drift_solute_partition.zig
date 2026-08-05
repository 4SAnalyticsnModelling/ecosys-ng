const std = @import("std");

/// Compact snow-drift topology: 33 salt/complex plus eight phosphate fields.
pub const snow_drift_species_count = 41;

pub const BoundaryAxis = enum { east_west, north_south };

pub const Inputs = struct {
    axis: BoundaryAxis,
    source_snow_drift_water_m3_per_step: f64,
    forward_boundary_snow_water_m3_per_step: f64,
    minimum_snow_drift_water_m3_per_step: f64,
    source_snow_drift_solute_flux_amount_per_step: []const f64,
};

/// Exact forward-side (`NN=1`) translation of TRNSFRS.F lines 4481--4614.
/// Instantaneous flux is overwritten and then accumulated for REDIST, in
/// amount step-1. H4SiO4 has no snow-drift slot; the neutral amount name is
/// required because legacy free-phosphate entries are not uniformly molar.
pub fn partitionAndAccumulate(inputs: Inputs, boundary_flux_amount_per_step: []f64, accumulated_flux_amount_per_step: []f64) !void {
    _ = inputs.axis;
    if (inputs.source_snow_drift_solute_flux_amount_per_step.len != snow_drift_species_count or
        boundary_flux_amount_per_step.len != snow_drift_species_count or
        accumulated_flux_amount_per_step.len != snow_drift_species_count)
        return error.ForwardSnowDriftSolutePartitionDimensionMismatch;
    inline for (.{ inputs.source_snow_drift_water_m3_per_step, inputs.forward_boundary_snow_water_m3_per_step, inputs.minimum_snow_drift_water_m3_per_step }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteForwardSnowDriftSolutePartitionInput;
    if (inputs.source_snow_drift_water_m3_per_step < 0 or
        inputs.forward_boundary_snow_water_m3_per_step < 0 or
        inputs.minimum_snow_drift_water_m3_per_step < 0)
        return error.InvalidForwardSnowDriftSolutePartitionInput;
    for (inputs.source_snow_drift_solute_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteForwardSnowDriftSolutePartitionInput;
    for (accumulated_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteForwardSnowDriftSolutePartitionInput;

    if (!(inputs.source_snow_drift_water_m3_per_step > inputs.minimum_snow_drift_water_m3_per_step)) {
        @memset(boundary_flux_amount_per_step, 0);
        return;
    }
    const fraction = inputs.forward_boundary_snow_water_m3_per_step /
        inputs.source_snow_drift_water_m3_per_step;
    for (0..snow_drift_species_count) |species| {
        const flux = inputs.source_snow_drift_solute_flux_amount_per_step[species] * fraction;
        if (!std.math.isFinite(flux) or !std.math.isFinite(accumulated_flux_amount_per_step[species] + flux))
            return error.NonFiniteForwardSnowDriftSolutePartitionResult;
    }
    for (0..snow_drift_species_count) |species| {
        const flux = inputs.source_snow_drift_solute_flux_amount_per_step[species] * fraction;
        boundary_flux_amount_per_step[species] = flux;
        accumulated_flux_amount_per_step[species] += flux;
    }
}

fn fixture(source: []const f64) Inputs {
    return .{ .axis = .east_west, .source_snow_drift_water_m3_per_step = 4, .forward_boundary_snow_water_m3_per_step = 1, .minimum_snow_drift_water_m3_per_step = 0.01, .source_snow_drift_solute_flux_amount_per_step = source };
}

test "TRNSFRS forward snow drift partitions and accumulates all 41 fields" {
    var source: [snow_drift_species_count]f64 = undefined;
    for (&source, 0..) |*value, species| value.* = @floatFromInt(species + 1);
    var boundary = [_]f64{9} ** snow_drift_species_count;
    var accumulated = [_]f64{10} ** snow_drift_species_count;
    try partitionAndAccumulate(fixture(&source), &boundary, &accumulated);
    try std.testing.expectEqual(@as(f64, 0.25), boundary[0]);
    try std.testing.expectEqual(@as(f64, 10.25), accumulated[0]);
    try std.testing.expectEqual(@as(f64, 10.25), boundary[40]);
    try std.testing.expectEqual(@as(f64, 20.25), accumulated[40]);
}

test "inactive snow drift zeros instantaneous forward flux only" {
    const source = [_]f64{8} ** snow_drift_species_count;
    var boundary = [_]f64{9} ** snow_drift_species_count;
    var accumulated = [_]f64{10} ** snow_drift_species_count;
    var inputs = fixture(&source);
    inputs.source_snow_drift_water_m3_per_step = 0;
    try partitionAndAccumulate(inputs, &boundary, &accumulated);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** snow_drift_species_count), &boundary);
    try std.testing.expectEqual(@as(f64, 10), accumulated[0]);
}

test "late invalid snow-drift accumulator is atomic" {
    const source = [_]f64{8} ** snow_drift_species_count;
    var boundary = [_]f64{9} ** snow_drift_species_count;
    var accumulated = [_]f64{10} ** snow_drift_species_count;
    accumulated[40] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteForwardSnowDriftSolutePartitionInput, partitionAndAccumulate(fixture(&source), &boundary, &accumulated));
    try std.testing.expectEqual(@as(f64, 9), boundary[0]);
    try std.testing.expectEqual(@as(f64, 10), accumulated[0]);
}

test "forward snow partition rejects phantom H4SiO4 slot" {
    const source = [_]f64{8} ** (snow_drift_species_count + 1);
    var boundary = [_]f64{9} ** snow_drift_species_count;
    var accumulated = [_]f64{10} ** snow_drift_species_count;
    try std.testing.expectError(error.ForwardSnowDriftSolutePartitionDimensionMismatch, partitionAndAccumulate(fixture(&source), &boundary, &accumulated));
}
