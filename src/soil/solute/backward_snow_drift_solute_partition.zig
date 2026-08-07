const std = @import("std");

pub const snow_drift_species_count = 41;
pub const BoundaryAxis = enum { east_west, north_south };

pub const Inputs = struct {
    axis: BoundaryAxis,
    source_snow_drift_water_m3_per_step: f64,
    backward_boundary_snow_water_m3_per_step: f64,
    minimum_snow_drift_water_m3_per_step: f64,
    source_snow_drift_solute_flux_amount_per_step: []const f64,
};

/// Exact valid-neighbor backward-side (`NN=2`) translation of TRNSFRS.F
/// lines 4618--4752. H4SiO4 is absent from the compact 41-field snow topology.
/// Inactive outer snow drift leaves backward storage untouched.
pub fn partitionAndAccumulate(inputs: Inputs, boundary_flux_amount_per_step: []f64, accumulated_flux_amount_per_step: []f64) !void {
    _ = inputs.axis;
    if (inputs.source_snow_drift_solute_flux_amount_per_step.len != snow_drift_species_count or
        boundary_flux_amount_per_step.len != snow_drift_species_count or
        accumulated_flux_amount_per_step.len != snow_drift_species_count)
        return error.BackwardSnowDriftSolutePartitionDimensionMismatch;
    inline for (.{ inputs.source_snow_drift_water_m3_per_step, inputs.backward_boundary_snow_water_m3_per_step, inputs.minimum_snow_drift_water_m3_per_step }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteBackwardSnowDriftSolutePartitionInput;
    if (inputs.source_snow_drift_water_m3_per_step < 0 or
        inputs.backward_boundary_snow_water_m3_per_step < 0 or
        inputs.minimum_snow_drift_water_m3_per_step < 0)
        return error.InvalidBackwardSnowDriftSolutePartitionInput;
    for (inputs.source_snow_drift_solute_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteBackwardSnowDriftSolutePartitionInput;
    for (accumulated_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteBackwardSnowDriftSolutePartitionInput;

    if (!(inputs.source_snow_drift_water_m3_per_step > inputs.minimum_snow_drift_water_m3_per_step)) return;
    const fraction = inputs.backward_boundary_snow_water_m3_per_step /
        inputs.source_snow_drift_water_m3_per_step;
    for (0..snow_drift_species_count) |species| {
        const flux = inputs.source_snow_drift_solute_flux_amount_per_step[species] * fraction;
        if (!std.math.isFinite(flux) or !std.math.isFinite(accumulated_flux_amount_per_step[species] + flux))
            return error.NonFiniteBackwardSnowDriftSolutePartitionResult;
    }
    for (0..snow_drift_species_count) |species| {
        const flux = inputs.source_snow_drift_solute_flux_amount_per_step[species] * fraction;
        boundary_flux_amount_per_step[species] = flux;
        accumulated_flux_amount_per_step[species] += flux;
    }
}

fn fixture(source: []const f64) Inputs {
    return .{ .axis = .north_south, .source_snow_drift_water_m3_per_step = 5, .backward_boundary_snow_water_m3_per_step = 2, .minimum_snow_drift_water_m3_per_step = 0.01, .source_snow_drift_solute_flux_amount_per_step = source };
}

test "TRNSFRS backward snow drift partitions and accumulates all 41 fields" {
    var source: [snow_drift_species_count]f64 = undefined;
    for (&source, 0..) |*value, species| value.* = @floatFromInt(species + 1);
    var boundary = [_]f64{9} ** snow_drift_species_count;
    var accumulated = [_]f64{10} ** snow_drift_species_count;
    try partitionAndAccumulate(fixture(&source), &boundary, &accumulated);
    try std.testing.expectEqual(@as(f64, 0.4), boundary[0]);
    try std.testing.expectEqual(@as(f64, 10.4), accumulated[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 16.4), boundary[40], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 26.4), accumulated[40], 1e-14);
}

test "inactive snow drift preserves backward instantaneous and accumulation" {
    const source = [_]f64{8} ** snow_drift_species_count;
    var boundary = [_]f64{9} ** snow_drift_species_count;
    var accumulated = [_]f64{10} ** snow_drift_species_count;
    var inputs = fixture(&source);
    inputs.source_snow_drift_water_m3_per_step = 0;
    try partitionAndAccumulate(inputs, &boundary, &accumulated);
    try std.testing.expectEqual(@as(f64, 9), boundary[0]);
    try std.testing.expectEqual(@as(f64, 10), accumulated[0]);
}

test "late invalid backward snow accumulation is atomic" {
    const source = [_]f64{8} ** snow_drift_species_count;
    var boundary = [_]f64{9} ** snow_drift_species_count;
    var accumulated = [_]f64{10} ** snow_drift_species_count;
    accumulated[40] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteBackwardSnowDriftSolutePartitionInput, partitionAndAccumulate(fixture(&source), &boundary, &accumulated));
    try std.testing.expectEqual(@as(f64, 9), boundary[0]);
    try std.testing.expectEqual(@as(f64, 10), accumulated[0]);
}

test "backward snow partition rejects phantom H4SiO4 slot" {
    const source = [_]f64{8} ** (snow_drift_species_count + 1);
    var boundary = [_]f64{9} ** snow_drift_species_count;
    var accumulated = [_]f64{10} ** snow_drift_species_count;
    try std.testing.expectError(error.BackwardSnowDriftSolutePartitionDimensionMismatch, partitionAndAccumulate(fixture(&source), &boundary, &accumulated));
}
