const std = @import("std");

pub const State = struct {
    leaf_carbon_g_c_by_layer_branch: []f64,
};

pub const Inputs = struct {
    branch_count: usize,
    canopy_layer_count: usize,
    leaf_node_pool_count: usize,
    leaf_carbon_g_c_by_layer_node_branch: []const f64,
};

fn aggregateIndex(inputs: Inputs, layer: usize, branch: usize) usize {
    return layer * inputs.branch_count + branch;
}

fn nodeIndex(inputs: Inputs, layer: usize, node: usize, branch: usize) usize {
    return (layer * inputs.leaf_node_pool_count + node) * inputs.branch_count + branch;
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.branch_count == 0 or inputs.canopy_layer_count == 0 or inputs.leaf_node_pool_count == 0) return error.CanopyLeafAggregationDimensionMismatch;
    const layer_nodes = std.math.mul(usize, inputs.canopy_layer_count, inputs.leaf_node_pool_count) catch return error.CanopyLeafAggregationDimensionOverflow;
    const count = std.math.mul(usize, layer_nodes, inputs.branch_count) catch return error.CanopyLeafAggregationDimensionOverflow;
    if (inputs.leaf_carbon_g_c_by_layer_node_branch.len != count) return error.CanopyLeafAggregationDimensionMismatch;
    for (inputs.leaf_carbon_g_c_by_layer_node_branch) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidCanopyLeafAggregationInput;
}

fn validateState(state: State, inputs: Inputs) !void {
    const count = std.math.mul(usize, inputs.canopy_layer_count, inputs.branch_count) catch return error.CanopyLeafAggregationDimensionOverflow;
    if (state.leaf_carbon_g_c_by_layer_branch.len != count) return error.CanopyLeafAggregationDimensionMismatch;
    for (state.leaf_carbon_g_c_by_layer_branch) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidCanopyLeafAggregationState;
}

/// Exact GROSUB 8791--8802 reset and aggregation of leaf-node carbon into
/// canopy layer x branch totals. Runtime node pools replace K=0..25. Traversal
/// remains branch, layer, node and carbon mass is g C.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateInputs(inputs);
    try validateState(state, inputs);
    try validateState(workspace, inputs);
    @memcpy(workspace.leaf_carbon_g_c_by_layer_branch, state.leaf_carbon_g_c_by_layer_branch);
    for (0..inputs.branch_count) |branch| for (0..inputs.canopy_layer_count) |layer| for (0..inputs.leaf_node_pool_count) |_| {
        workspace.leaf_carbon_g_c_by_layer_branch[aggregateIndex(inputs, layer, branch)] = 0.0;
    };
    for (0..inputs.branch_count) |branch| for (0..inputs.canopy_layer_count) |layer| for (0..inputs.leaf_node_pool_count) |node| {
        const target = aggregateIndex(inputs, layer, branch);
        workspace.leaf_carbon_g_c_by_layer_branch[target] += inputs.leaf_carbon_g_c_by_layer_node_branch[nodeIndex(inputs, layer, node, branch)];
    };
    try validateState(workspace, inputs);
    @memcpy(state.leaf_carbon_g_c_by_layer_branch, workspace.leaf_carbon_g_c_by_layer_branch);
}

test "GROSUB canopy leaf aggregation preserves branch layer node traversal" {
    var totals = [_]f64{ 9, 9, 9, 9 };
    var work = [_]f64{0} ** 4;
    const inputs: Inputs = .{ .branch_count = 2, .canopy_layer_count = 2, .leaf_node_pool_count = 3, .leaf_carbon_g_c_by_layer_node_branch = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } };
    try apply(.{ .leaf_carbon_g_c_by_layer_branch = &totals }, .{ .leaf_carbon_g_c_by_layer_branch = &work }, inputs);
    try std.testing.expectEqualSlices(f64, &.{ 9, 12, 27, 30 }, &totals);
}
