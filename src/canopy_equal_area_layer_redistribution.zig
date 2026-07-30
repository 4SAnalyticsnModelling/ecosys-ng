const std = @import("std");

pub const Inputs = struct {
    canopy_height_m: f64,
    combined_leaf_area_m2: f64,
    combined_stalk_area_m2: f64,
    combined_standing_dead_area_m2: f64,
    leaf_area_by_layer_m2: []const f64,
    stalk_area_by_layer_m2: []const f64,
    standing_dead_area_by_layer_m2: []const f64,
    negligible_area_m2: f64,
};

/// HOUR1 control/comment lines 1775--1786 and executable lines 1787--1827.
/// Caller-provided scratch has the same runtime `layer_count + 1` extent.
pub fn apply(
    inputs: Inputs,
    layer_boundary_height_m: []f64,
    next_layer_boundary_height_m: []f64,
) !void {
    try validate(inputs, layer_boundary_height_m, next_layer_boundary_height_m);
    const layer_count = inputs.leaf_area_by_layer_m2.len;
    layer_boundary_height_m[layer_count] = inputs.canopy_height_m + 0.01;
    next_layer_boundary_height_m[layer_count] =
        layer_boundary_height_m[layer_count];
    next_layer_boundary_height_m[0] = 0.0;
    const target_area_m2 =
        (inputs.combined_leaf_area_m2 + inputs.combined_stalk_area_m2 +
            inputs.combined_standing_dead_area_m2) /
        @as(f64, @floatFromInt(layer_count));
    if (target_area_m2 > inputs.negligible_area_m2) {
        var upper = layer_count;
        while (upper > 1) : (upper -= 1) {
            const layer = upper - 1;
            const area_m2 = inputs.leaf_area_by_layer_m2[layer] +
                inputs.stalk_area_by_layer_m2[layer] +
                inputs.standing_dead_area_by_layer_m2[layer];
            if (area_m2 > 1.01 * target_area_m2) {
                const thickness_m = layer_boundary_height_m[layer + 1] -
                    layer_boundary_height_m[layer];
                next_layer_boundary_height_m[layer] =
                    layer_boundary_height_m[layer] +
                    0.5 * @min(1.0, (area_m2 - target_area_m2) / area_m2) *
                        thickness_m;
            } else if (area_m2 < 0.99 * target_area_m2) {
                const lower_area_m2 =
                    inputs.leaf_area_by_layer_m2[layer - 1] +
                    inputs.stalk_area_by_layer_m2[layer - 1] +
                    inputs.standing_dead_area_by_layer_m2[layer - 1];
                const lower_thickness_m = layer_boundary_height_m[layer] -
                    layer_boundary_height_m[layer - 1];
                next_layer_boundary_height_m[layer] =
                    if (lower_area_m2 > inputs.negligible_area_m2)
                        @min(
                            inputs.canopy_height_m,
                            layer_boundary_height_m[layer] -
                                0.5 * @min(
                                    1.0,
                                    (target_area_m2 - area_m2) / lower_area_m2,
                                ) * lower_thickness_m,
                        )
                    else
                        next_layer_boundary_height_m[layer + 1];
            } else {
                next_layer_boundary_height_m[layer] =
                    layer_boundary_height_m[layer];
            }
        }
        upper = layer_count;
        while (upper > 1) : (upper -= 1)
            layer_boundary_height_m[upper - 1] =
                next_layer_boundary_height_m[upper - 1];
    } else {
        var upper = layer_count;
        while (upper > 1) : (upper -= 1)
            layer_boundary_height_m[upper - 1] = 0.0;
    }
}

fn validate(inputs: Inputs, boundaries: []const f64, scratch: []const f64) !void {
    const layer_count = inputs.leaf_area_by_layer_m2.len;
    const boundary_count = try std.math.add(usize, layer_count, 1);
    if (layer_count < 2 or
        inputs.stalk_area_by_layer_m2.len != layer_count or
        inputs.standing_dead_area_by_layer_m2.len != layer_count or
        boundaries.len != boundary_count or scratch.len != boundary_count)
        return error.CanopyLayerRedistributionDimensionMismatch;
    inline for (.{
        inputs.canopy_height_m,
        inputs.combined_leaf_area_m2,
        inputs.combined_stalk_area_m2,
        inputs.combined_standing_dead_area_m2,
        inputs.negligible_area_m2,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopyLayerRedistributionInput;
    inline for (.{
        inputs.leaf_area_by_layer_m2,
        inputs.stalk_area_by_layer_m2,
        inputs.standing_dead_area_by_layer_m2,
        boundaries,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopyLayerRedistributionInput;
}

test "high top-layer area moves its lower boundary upward" {
    var boundaries = [_]f64{ 0, 0.5, 1.0, 1.5 };
    var scratch: [4]f64 = undefined;
    try apply(.{
        .canopy_height_m = 1.5,
        .combined_leaf_area_m2 = 9,
        .combined_stalk_area_m2 = 0,
        .combined_standing_dead_area_m2 = 0,
        .leaf_area_by_layer_m2 = &.{ 1, 1, 7 },
        .stalk_area_by_layer_m2 = &.{ 0, 0, 0 },
        .standing_dead_area_by_layer_m2 = &.{ 0, 0, 0 },
        .negligible_area_m2 = 1e-12,
    }, &boundaries, &scratch);
    try std.testing.expect(boundaries[2] > 1.0);
    try std.testing.expectEqual(@as(f64, 1.51), boundaries[3]);
}

test "negligible canopy zeros internal boundaries only" {
    var boundaries = [_]f64{ 0, 0.5, 1.0, 1.5 };
    var scratch: [4]f64 = undefined;
    try apply(.{
        .canopy_height_m = 1.5,
        .combined_leaf_area_m2 = 0,
        .combined_stalk_area_m2 = 0,
        .combined_standing_dead_area_m2 = 0,
        .leaf_area_by_layer_m2 = &.{ 0, 0, 0 },
        .stalk_area_by_layer_m2 = &.{ 0, 0, 0 },
        .standing_dead_area_by_layer_m2 = &.{ 0, 0, 0 },
        .negligible_area_m2 = 1e-12,
    }, &boundaries, &scratch);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0, 1.51 }, &boundaries);
}
