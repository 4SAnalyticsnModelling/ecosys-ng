const std = @import("std");

pub const State = struct {
    community_height_m: []f64,
    /// Cell-major `[cell][layer_boundary]`, including source boundary zero.
    layer_boundary_height_m: []f64,
    /// Cell-major `[cell][canopy_layer]`.
    living_leaf_area_m2: []f64,
    living_stalk_area_m2: []f64,
    standing_dead_area_m2: []f64,
    living_leaf_carbon_g_c: []f64,
};

/// Exact source-order translation of legacy `STARTS` lines 665--673.
pub fn initialize(
    state: State,
    cell_count: usize,
    canopy_layer_count: usize,
) !void {
    if (cell_count == 0 or canopy_layer_count == 0)
        return error.InvalidCommunityCanopyDimensions;
    const layer_value_count = std.math.mul(
        usize,
        cell_count,
        canopy_layer_count,
    ) catch return error.DimensionOverflow;
    const boundary_count_per_cell = std.math.add(
        usize,
        canopy_layer_count,
        1,
    ) catch return error.DimensionOverflow;
    const boundary_value_count = std.math.mul(
        usize,
        cell_count,
        boundary_count_per_cell,
    ) catch return error.DimensionOverflow;
    if (state.community_height_m.len != cell_count or
        state.layer_boundary_height_m.len != boundary_value_count or
        state.living_leaf_area_m2.len != layer_value_count or
        state.living_stalk_area_m2.len != layer_value_count or
        state.standing_dead_area_m2.len != layer_value_count or
        state.living_leaf_carbon_g_c.len != layer_value_count)
    {
        return error.CommunityCanopyDimensionMismatch;
    }

    for (0..cell_count) |cell| {
        state.community_height_m[cell] = 0.0;
        const boundary_base = cell * boundary_count_per_cell;
        state.layer_boundary_height_m[boundary_base] = 0.0;
        for (0..canopy_layer_count) |layer| {
            state.layer_boundary_height_m[boundary_base + layer + 1] = 0.0;
            const index = cell * canopy_layer_count + layer;
            state.living_leaf_area_m2[index] = 0.0;
            state.living_stalk_area_m2[index] = 0.0;
            state.standing_dead_area_m2[index] = 0.0;
            state.living_leaf_carbon_g_c[index] = 0.0;
        }
    }
}

test "STARTS resets community canopy in exact boundary and layer order" {
    const cells = 2;
    const layers = 3;
    var community_height = [_]f64{7.0} ** cells;
    var boundary_height = [_]f64{7.0} ** (cells * (layers + 1));
    var leaf_area = [_]f64{7.0} ** (cells * layers);
    var stalk_area = [_]f64{7.0} ** (cells * layers);
    var dead_area = [_]f64{7.0} ** (cells * layers);
    var leaf_carbon = [_]f64{7.0} ** (cells * layers);

    try initialize(.{
        .community_height_m = &community_height,
        .layer_boundary_height_m = &boundary_height,
        .living_leaf_area_m2 = &leaf_area,
        .living_stalk_area_m2 = &stalk_area,
        .standing_dead_area_m2 = &dead_area,
        .living_leaf_carbon_g_c = &leaf_carbon,
    }, cells, layers);

    inline for (.{
        community_height,
        boundary_height,
        leaf_area,
        stalk_area,
        dead_area,
        leaf_carbon,
    }) |values| {
        for (values) |value| try std.testing.expectEqual(@as(f64, 0.0), value);
    }
}

test "dimension mismatch is detected before any canopy reset" {
    var community_height = [_]f64{ 1, 2 };
    var short_boundaries = [_]f64{ 3, 4, 5 };
    var layers = [_]f64{ 6, 7, 8, 9 };
    const before = community_height;
    try std.testing.expectError(
        error.CommunityCanopyDimensionMismatch,
        initialize(.{
            .community_height_m = &community_height,
            .layer_boundary_height_m = &short_boundaries,
            .living_leaf_area_m2 = &layers,
            .living_stalk_area_m2 = &layers,
            .standing_dead_area_m2 = &layers,
            .living_leaf_carbon_g_c = &layers,
        }, 2, 2),
    );
    try std.testing.expectEqualSlices(f64, &before, &community_height);
}
