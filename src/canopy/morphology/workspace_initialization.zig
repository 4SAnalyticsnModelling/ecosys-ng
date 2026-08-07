const std = @import("std");

pub const Dimensions = struct {
    branch_count: usize,
    canopy_layer_count: usize,
    node_position_count: usize,
    surface_component_count: usize,
};

pub const Workspace = struct {
    standing_stalk_area_m2: []f64,
    branch_surface_area_m2: []f64,
    standing_dead_surface_area_m2: []f64,
    leaf_area_m2: []f64,
    node_height_m: []f64,
    maximum_node_height_m: []f64,
    sheath_height_m: []f64,
    living_leaf_carbon_g_c: []f64,
    senescing_leaf_carbon_g_c: []f64,
    living_leaf_nitrogen_g_n: []f64,
    living_leaf_phosphorus_g_p: []f64,
    living_sheath_carbon_g_c: []f64,
    senescing_sheath_carbon_g_c: []f64,
    living_sheath_nitrogen_g_n: []f64,
    living_sheath_phosphorus_g_p: []f64,
    node_carbon_g_c: []f64,
    node_nitrogen_g_n: []f64,
    node_phosphorus_g_p: []f64,
    layer_leaf_area_m2: []f64,
    layer_leaf_carbon_g_c: []f64,
    layer_leaf_nitrogen_g_n: []f64,
    layer_leaf_phosphorus_g_p: []f64,
    node_carbon_pool_3_g_c: []f64,
    node_carbon_dioxide_g_c: []f64,
    node_bicarbonate_g_c: []f64,
    node_carbon_pool_4_g_c: []f64,
    node_surface_area_m2: []f64,
};

pub const InitializationError = error{
    ExtentOverflow,
    ExtentMismatch,
    MissingBaseNodePosition,
};

/// Translates `startq.f` lines 546--586.
///
/// Flattened order is branch-major. Within each branch, node position precedes
/// canopy layer where both are present; surface component is the innermost
/// dimension for surface-area arrays.
pub fn initialize(
    dimensions: Dimensions,
    workspace: Workspace,
) InitializationError!void {
    if (dimensions.node_position_count == 0) return error.MissingBaseNodePosition;
    const branch_layer = try product(dimensions.branch_count, dimensions.canopy_layer_count);
    const branch_layer_surface = try product(branch_layer, dimensions.surface_component_count);
    const layer_surface = try product(
        dimensions.canopy_layer_count,
        dimensions.surface_component_count,
    );
    const branch_node = try product(dimensions.branch_count, dimensions.node_position_count);
    const branch_nonbase_node = try product(
        dimensions.branch_count,
        dimensions.node_position_count - 1,
    );
    const branch_node_layer = try product(branch_node, dimensions.canopy_layer_count);
    const branch_node_layer_surface = try product(
        try product(branch_nonbase_node, dimensions.canopy_layer_count),
        dimensions.surface_component_count,
    );
    try requireExtent(workspace.standing_stalk_area_m2, branch_layer);
    try requireExtent(workspace.branch_surface_area_m2, branch_layer_surface);
    try requireExtent(workspace.standing_dead_surface_area_m2, layer_surface);
    inline for (.{
        workspace.leaf_area_m2,
        workspace.node_height_m,
        workspace.maximum_node_height_m,
        workspace.sheath_height_m,
        workspace.living_leaf_carbon_g_c,
        workspace.senescing_leaf_carbon_g_c,
        workspace.living_leaf_nitrogen_g_n,
        workspace.living_leaf_phosphorus_g_p,
        workspace.living_sheath_carbon_g_c,
        workspace.senescing_sheath_carbon_g_c,
        workspace.living_sheath_nitrogen_g_n,
        workspace.living_sheath_phosphorus_g_p,
        workspace.node_carbon_g_c,
        workspace.node_nitrogen_g_n,
        workspace.node_phosphorus_g_p,
    }) |values| try requireExtent(values, branch_node);
    inline for (.{
        workspace.node_carbon_pool_3_g_c,
        workspace.node_carbon_dioxide_g_c,
        workspace.node_bicarbonate_g_c,
        workspace.node_carbon_pool_4_g_c,
    }) |values| try requireExtent(values, branch_nonbase_node);
    inline for (.{
        workspace.layer_leaf_area_m2,
        workspace.layer_leaf_carbon_g_c,
        workspace.layer_leaf_nitrogen_g_n,
        workspace.layer_leaf_phosphorus_g_p,
    }) |values| try requireExtent(values, branch_node_layer);
    try requireExtent(workspace.node_surface_area_m2, branch_node_layer_surface);

    for (0..dimensions.branch_count) |branch| {
        for (0..dimensions.canopy_layer_count) |layer| {
            workspace.standing_stalk_area_m2[branch * dimensions.canopy_layer_count + layer] = 0.0;
            for (0..dimensions.surface_component_count) |component| {
                const branch_surface_index =
                    (branch * dimensions.canopy_layer_count + layer) *
                    dimensions.surface_component_count + component;
                workspace.branch_surface_area_m2[branch_surface_index] = 0.0;
                if (branch == 0) {
                    workspace.standing_dead_surface_area_m2[
                        layer * dimensions.surface_component_count + component
                    ] = 0.0;
                }
            }
        }
        for (0..dimensions.node_position_count) |node| {
            const node_index = branch * dimensions.node_position_count + node;
            inline for (.{
                workspace.leaf_area_m2,
                workspace.node_height_m,
                workspace.maximum_node_height_m,
                workspace.sheath_height_m,
                workspace.living_leaf_carbon_g_c,
                workspace.senescing_leaf_carbon_g_c,
                workspace.living_leaf_nitrogen_g_n,
                workspace.living_leaf_phosphorus_g_p,
                workspace.living_sheath_carbon_g_c,
                workspace.senescing_sheath_carbon_g_c,
                workspace.living_sheath_nitrogen_g_n,
                workspace.living_sheath_phosphorus_g_p,
                workspace.node_carbon_g_c,
                workspace.node_nitrogen_g_n,
                workspace.node_phosphorus_g_p,
            }) |values| values[node_index] = 0.0;
            for (0..dimensions.canopy_layer_count) |layer| {
                const layer_node_index =
                    node_index * dimensions.canopy_layer_count + layer;
                workspace.layer_leaf_area_m2[layer_node_index] = 0.0;
                workspace.layer_leaf_carbon_g_c[layer_node_index] = 0.0;
                workspace.layer_leaf_nitrogen_g_n[layer_node_index] = 0.0;
                workspace.layer_leaf_phosphorus_g_p[layer_node_index] = 0.0;
            }
            if (node != 0) {
                const nonbase_node_index =
                    branch * (dimensions.node_position_count - 1) + node - 1;
                workspace.node_carbon_pool_3_g_c[nonbase_node_index] = 0.0;
                workspace.node_carbon_dioxide_g_c[nonbase_node_index] = 0.0;
                workspace.node_bicarbonate_g_c[nonbase_node_index] = 0.0;
                workspace.node_carbon_pool_4_g_c[nonbase_node_index] = 0.0;
                for (0..dimensions.canopy_layer_count) |layer| {
                    for (0..dimensions.surface_component_count) |component| {
                        const surface_index =
                            (nonbase_node_index * dimensions.canopy_layer_count + layer) *
                            dimensions.surface_component_count + component;
                        workspace.node_surface_area_m2[surface_index] = 0.0;
                    }
                }
            }
        }
    }
}

fn product(left: usize, right: usize) InitializationError!usize {
    return std.math.mul(usize, left, right) catch error.ExtentOverflow;
}

fn requireExtent(values: []f64, expected: usize) InitializationError!void {
    if (values.len != expected) return error.ExtentMismatch;
}

test "runtime morphology workspace resets source-owned extents" {
    const dimensions = Dimensions{
        .branch_count = 2,
        .canopy_layer_count = 2,
        .node_position_count = 3,
        .surface_component_count = 2,
    };
    var storage: [186]f64 = undefined;
    @memset(&storage, 7.0);
    var cursor: usize = 0;
    const workspace = makeWorkspace(dimensions, &storage, &cursor);
    try std.testing.expectEqual(storage.len, cursor);

    try initialize(dimensions, workspace);

    try std.testing.expectEqualSlices(
        f64,
        &[_]f64{0.0} ** storage.len,
        &storage,
    );
}

test "extent mismatch fails before mutation" {
    const dimensions = Dimensions{
        .branch_count = 0,
        .canopy_layer_count = 0,
        .node_position_count = 1,
        .surface_component_count = 0,
    };
    var storage: [1]f64 = .{42.0};
    var cursor: usize = 0;
    var workspace = makeWorkspace(dimensions, storage[0..0], &cursor);
    workspace.leaf_area_m2 = &storage;
    try std.testing.expectError(error.ExtentMismatch, initialize(dimensions, workspace));
    try std.testing.expectEqual(@as(f64, 42.0), storage[0]);
}

fn makeWorkspace(
    dimensions: Dimensions,
    storage: []f64,
    cursor: *usize,
) Workspace {
    const branch_layer = dimensions.branch_count * dimensions.canopy_layer_count;
    const branch_layer_surface = branch_layer * dimensions.surface_component_count;
    const layer_surface = dimensions.canopy_layer_count * dimensions.surface_component_count;
    const branch_node = dimensions.branch_count * dimensions.node_position_count;
    const branch_nonbase_node =
        dimensions.branch_count * (dimensions.node_position_count - 1);
    const branch_node_layer = branch_node * dimensions.canopy_layer_count;
    const branch_node_layer_surface =
        branch_nonbase_node * dimensions.canopy_layer_count *
        dimensions.surface_component_count;
    var result: Workspace = undefined;
    inline for (std.meta.fields(Workspace)) |field| {
        const length =
            if (std.mem.eql(u8, field.name, "standing_stalk_area_m2"))
                branch_layer
            else if (std.mem.eql(u8, field.name, "branch_surface_area_m2"))
                branch_layer_surface
            else if (std.mem.eql(u8, field.name, "standing_dead_surface_area_m2"))
                layer_surface
            else if (std.mem.startsWith(u8, field.name, "layer_leaf_"))
                branch_node_layer
            else if (std.mem.eql(u8, field.name, "node_surface_area_m2"))
                branch_node_layer_surface
            else if (std.mem.startsWith(u8, field.name, "node_carbon_pool_") or
            std.mem.eql(u8, field.name, "node_carbon_dioxide_g_c") or
            std.mem.eql(u8, field.name, "node_bicarbonate_g_c"))
                branch_nonbase_node
            else
                branch_node;
        @field(result, field.name) = storage[cursor.* .. cursor.* + length];
        cursor.* += length;
    }
    return result;
}
