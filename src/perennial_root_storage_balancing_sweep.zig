const std = @import("std");

pub const State = struct {
    mobile_carbon_g_c_by_domain_layer: []f64,
    seasonal_storage_carbon_g_c: []f64,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    plant_deepest_rooted_layer_index: usize,
    is_perennial: bool,
    active_root_carbon_g_c_by_domain_layer: []const f64,
    plant_total_root_carbon_g_c: f64,
    storage_deficit_threshold_g_c_per_g_root_c: f64,
    exchange_fraction_per_h: f64,
    biological_timestep_h: f64,
    presence_threshold_g_c: f64,
};

fn validate(state: State, inputs: Inputs) !void {
    const values = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.RootStorageBalancingDimensionOverflow;
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count or inputs.plant_deepest_rooted_layer_index >= inputs.soil_layer_count or state.mobile_carbon_g_c_by_domain_layer.len != values or state.seasonal_storage_carbon_g_c.len != 1 or inputs.active_root_carbon_g_c_by_domain_layer.len != values) return error.RootStorageBalancingDimensionMismatch;
    inline for (.{ inputs.plant_total_root_carbon_g_c, inputs.storage_deficit_threshold_g_c_per_g_root_c, inputs.exchange_fraction_per_h, inputs.biological_timestep_h, inputs.presence_threshold_g_c, state.seasonal_storage_carbon_g_c[0] }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootStorageBalancingInput;
    for (inputs.active_root_carbon_g_c_by_domain_layer) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootStorageBalancingInput;
    for (state.mobile_carbon_g_c_by_domain_layer) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootStorageBalancingInput;
}

/// Exact GROSUB 7379--7400 perennial root-to-seasonal-storage balancing.
/// Traversal is biological domain N then rooted layer L. The plant-wide
/// seasonal C pool is advanced after each layer, preserving source order.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validate(state, inputs);
    try validate(workspace, inputs);
    @memcpy(workspace.mobile_carbon_g_c_by_domain_layer, state.mobile_carbon_g_c_by_domain_layer);
    workspace.seasonal_storage_carbon_g_c[0] = state.seasonal_storage_carbon_g_c[0];
    if (inputs.is_perennial) for (0..inputs.biological_domain_count) |domain| for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
        if (layer > inputs.plant_deepest_rooted_layer_index) continue;
        const index = domain * inputs.soil_layer_count + layer;
        const layer_root_g_c = inputs.active_root_carbon_g_c_by_domain_layer[index];
        if (layer_root_g_c <= inputs.presence_threshold_g_c or inputs.plant_total_root_carbon_g_c <= inputs.presence_threshold_g_c or workspace.seasonal_storage_carbon_g_c[0] >= inputs.storage_deficit_threshold_g_c_per_g_root_c * inputs.plant_total_root_carbon_g_c) continue;
        const layer_fraction = layer_root_g_c / inputs.plant_total_root_carbon_g_c;
        const represented_total_root_g_c = inputs.plant_total_root_carbon_g_c * layer_fraction;
        const combined_g_c = layer_root_g_c + represented_total_root_g_c;
        const mobile_g_c = @max(0, workspace.mobile_carbon_g_c_by_domain_layer[index]);
        const represented_storage_g_c = @max(0, workspace.seasonal_storage_carbon_g_c[0] * layer_fraction);
        const gradient_g_c = (represented_storage_g_c * layer_root_g_c - mobile_g_c * represented_total_root_g_c) / combined_g_c;
        const signed_root_change_g_c = @min(0, inputs.exchange_fraction_per_h * gradient_g_c * inputs.biological_timestep_h);
        workspace.mobile_carbon_g_c_by_domain_layer[index] += signed_root_change_g_c;
        workspace.seasonal_storage_carbon_g_c[0] -= signed_root_change_g_c;
        if (!std.math.isFinite(workspace.mobile_carbon_g_c_by_domain_layer[index]) or workspace.mobile_carbon_g_c_by_domain_layer[index] < 0 or !std.math.isFinite(workspace.seasonal_storage_carbon_g_c[0])) return error.RootStorageBalancingWouldOverdraw;
    };
    @memcpy(state.mobile_carbon_g_c_by_domain_layer, workspace.mobile_carbon_g_c_by_domain_layer);
    state.seasonal_storage_carbon_g_c[0] = workspace.seasonal_storage_carbon_g_c[0];
}

test "GROSUB perennial root storage sweep conserves C in N L source order" {
    var mobile = [_]f64{ 4, 2, 3, 1 };
    var storage = [_]f64{5};
    var wm = [_]f64{0} ** 4;
    var ws = [_]f64{0};
    const before = storage[0] + mobile[0] + mobile[1] + mobile[2] + mobile[3];
    try apply(.{ .mobile_carbon_g_c_by_domain_layer = &mobile, .seasonal_storage_carbon_g_c = &storage }, .{ .mobile_carbon_g_c_by_domain_layer = &wm, .seasonal_storage_carbon_g_c = &ws }, .{ .biological_domain_count = 2, .soil_layer_count = 2, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .plant_deepest_rooted_layer_index = 1, .is_perennial = true, .active_root_carbon_g_c_by_domain_layer = &.{ 20, 10, 5, 5 }, .plant_total_root_carbon_g_c = 100, .storage_deficit_threshold_g_c_per_g_root_c = 0.1, .exchange_fraction_per_h = 0.5, .biological_timestep_h = 0.25, .presence_threshold_g_c = 1e-12 });
    try std.testing.expectApproxEqAbs(before, storage[0] + mobile[0] + mobile[1] + mobile[2] + mobile[3], 1e-12);
    try std.testing.expect(storage[0] > 5);
}

test "GROSUB root storage balancing rolls back late donor overdraw" {
    var mobile = [_]f64{ 1, 0.01 };
    var storage = [_]f64{0};
    var wm = [_]f64{0} ** 2;
    var ws = [_]f64{0};
    const before = mobile;
    try std.testing.expectError(error.RootStorageBalancingWouldOverdraw, apply(.{ .mobile_carbon_g_c_by_domain_layer = &mobile, .seasonal_storage_carbon_g_c = &storage }, .{ .mobile_carbon_g_c_by_domain_layer = &wm, .seasonal_storage_carbon_g_c = &ws }, .{ .biological_domain_count = 1, .soil_layer_count = 2, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .plant_deepest_rooted_layer_index = 1, .is_perennial = true, .active_root_carbon_g_c_by_domain_layer = &.{ 10, 10 }, .plant_total_root_carbon_g_c = 20, .storage_deficit_threshold_g_c_per_g_root_c = 1, .exchange_fraction_per_h = 100, .biological_timestep_h = 1, .presence_threshold_g_c = 0 }));
    try std.testing.expectEqualDeep(before, mobile);
    try std.testing.expectEqual(@as(f64, 0), storage[0]);
}
