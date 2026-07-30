const std = @import("std");

pub const Dimensions = struct {
    root_axis_count: usize,
    soil_layer_count: usize,
    first_soil_layer: usize,
    end_soil_layer: usize,
};

/// UPTAKE.F 1339--1342. Adds one M-substep into runtime-sized root water
/// uptake totals in source root-axis-major, soil-layer-minor order.
pub fn addSubstep(
    dimensions: Dimensions,
    substep_uptake_m3_per_step_by_root_layer: []const f64,
    total_uptake_m3_per_step_by_root_layer: []f64,
) !void {
    const count = std.math.mul(
        usize,
        dimensions.root_axis_count,
        dimensions.soil_layer_count,
    ) catch return error.RootWaterUptakeAggregationDimensionMismatch;
    if (substep_uptake_m3_per_step_by_root_layer.len != count or
        total_uptake_m3_per_step_by_root_layer.len != count or
        dimensions.first_soil_layer > dimensions.end_soil_layer or
        dimensions.end_soil_layer > dimensions.soil_layer_count)
        return error.RootWaterUptakeAggregationDimensionMismatch;
    for (0..dimensions.root_axis_count) |root_axis| {
        for (dimensions.first_soil_layer..dimensions.end_soil_layer) |layer| {
            const index = root_axis * dimensions.soil_layer_count + layer;
            const substep = substep_uptake_m3_per_step_by_root_layer[index];
            const total = total_uptake_m3_per_step_by_root_layer[index];
            if (!std.math.isFinite(substep) or !std.math.isFinite(total))
                return error.InvalidRootWaterUptakeAggregationInput;
            const updated = total + substep;
            if (!std.math.isFinite(updated))
                return error.NonFiniteRootWaterUptakeAggregationResult;
            total_uptake_m3_per_step_by_root_layer[index] = updated;
        }
    }
}

test "UPTAKE root uptake aggregation preserves runtime source traversal" {
    const substep = [_]f64{ 99, -1, -2, 99, 99, 3, 4, 99 };
    var total = [_]f64{ 10, 10, 10, 10, 20, 20, 20, 20 };
    try addSubstep(
        .{
            .root_axis_count = 2,
            .soil_layer_count = 4,
            .first_soil_layer = 1,
            .end_soil_layer = 3,
        },
        &substep,
        &total,
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 10, 9, 8, 10, 20, 23, 24, 20 },
        &total,
    );
}

test "empty active layer interval leaves totals unchanged" {
    const substep = [_]f64{ 1, 2 };
    var total = [_]f64{ 3, 4 };
    try addSubstep(
        .{
            .root_axis_count = 1,
            .soil_layer_count = 2,
            .first_soil_layer = 1,
            .end_soil_layer = 1,
        },
        &substep,
        &total,
    );
    try std.testing.expectEqualSlices(f64, &.{ 3, 4 }, &total);
}

test "runtime dimension mismatch fails explicitly" {
    const substep = [_]f64{ 1, 2 };
    var total = [_]f64{3};
    try std.testing.expectError(
        error.RootWaterUptakeAggregationDimensionMismatch,
        addSubstep(
            .{
                .root_axis_count = 1,
                .soil_layer_count = 2,
                .first_soil_layer = 0,
                .end_soil_layer = 2,
            },
            &substep,
            &total,
        ),
    );
}
