const std = @import("std");

pub const State = struct {
    total_root_carbon_g_c_by_compartment_layer: []f64,
};

pub const Inputs = struct {
    association_compartment_count: usize,
    soil_layer_count: usize,
    shared_planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    root_mobile_carbon_g_c_by_compartment_layer: []const f64,
};

fn validateInputs(inputs: Inputs) !void {
    if (inputs.association_compartment_count == 0 or inputs.soil_layer_count == 0 or inputs.shared_planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count) return error.RootCarbonAugmentationDimensionMismatch;
    const count = std.math.mul(usize, inputs.association_compartment_count, inputs.soil_layer_count) catch return error.RootCarbonAugmentationDimensionOverflow;
    if (inputs.root_mobile_carbon_g_c_by_compartment_layer.len != count) return error.RootCarbonAugmentationDimensionMismatch;
    for (inputs.root_mobile_carbon_g_c_by_compartment_layer) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidRootCarbonAugmentationInput;
}

fn validateState(state: State, inputs: Inputs) !void {
    const count = std.math.mul(usize, inputs.association_compartment_count, inputs.soil_layer_count) catch return error.RootCarbonAugmentationDimensionOverflow;
    if (state.total_root_carbon_g_c_by_compartment_layer.len != count) return error.RootCarbonAugmentationDimensionMismatch;
    for (state.total_root_carbon_g_c_by_compartment_layer) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidRootCarbonAugmentationState;
}

/// Exact GROSUB 8508--8511 augmentation of structural root/mycorrhizal carbon
/// by non-structural carbon. Layout is [association compartment][soil layer],
/// traversal is compartment then shared NU..NI layer, and masses are g C.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateInputs(inputs);
    try validateState(state, inputs);
    try validateState(workspace, inputs);
    @memcpy(workspace.total_root_carbon_g_c_by_compartment_layer, state.total_root_carbon_g_c_by_compartment_layer);
    for (0..inputs.association_compartment_count) |compartment| {
        for (inputs.shared_planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
            const index = compartment * inputs.soil_layer_count + layer;
            workspace.total_root_carbon_g_c_by_compartment_layer[index] += inputs.root_mobile_carbon_g_c_by_compartment_layer[index];
        }
    }
    try validateState(workspace, inputs);
    @memcpy(state.total_root_carbon_g_c_by_compartment_layer, workspace.total_root_carbon_g_c_by_compartment_layer);
}

test "GROSUB root carbon augmentation uses every runtime compartment and rooted layer" {
    var totals = [_]f64{ 1, 2, 3, 4, 5, 6 };
    var workspace_values = [_]f64{0} ** 6;
    const inputs: Inputs = .{ .association_compartment_count = 2, .soil_layer_count = 3, .shared_planting_layer_index = 1, .deepest_rooted_layer_index = 2, .root_mobile_carbon_g_c_by_compartment_layer = &.{ 10, 20, 30, 40, 50, 60 } };
    try apply(.{ .total_root_carbon_g_c_by_compartment_layer = &totals }, .{ .total_root_carbon_g_c_by_compartment_layer = &workspace_values }, inputs);
    try std.testing.expectEqualSlices(f64, &.{ 1, 22, 33, 4, 55, 66 }, &totals);
}
