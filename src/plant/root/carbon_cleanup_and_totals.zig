const std = @import("std");

pub const State = struct {
    primary_carbon_g_c_by_domain_layer_axis: []f64,
    secondary_carbon_g_c_by_domain_layer_axis: []f64,
    primary_length_m_by_domain_layer_axis: []f64,
    secondary_length_m_by_domain_layer_axis: []f64,
    mobile_carbon_g_c_by_domain_layer: []f64,
    deepest_rooted_layer_index_by_axis: []usize,
    axis_closed_by_domain_axis: []bool,
    plant_deepest_rooted_layer_index: []usize,
    primary_carbon_total_g_c_by_domain_layer: []f64,
    primary_length_total_m_by_domain_layer: []f64,
    secondary_carbon_total_g_c_by_domain_layer: []f64,
    secondary_length_total_m_by_domain_layer: []f64,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    profile_bottom_layer_index: usize,
    host_domain_index: usize,
    layer_active: []const bool,
    primary_tip_branch_active_by_layer_axis: []const bool,
};

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateState(model_state: State, axis_values: usize, domain_layers: usize, domain_axes: usize, axes: usize) !void {
    inline for (.{ model_state.primary_carbon_g_c_by_domain_layer_axis, model_state.secondary_carbon_g_c_by_domain_layer_axis }) |values| {
        if (values.len != axis_values) return error.RootCarbonTotalsDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootCarbonTotalsState;
    }
    inline for (.{ model_state.primary_length_m_by_domain_layer_axis, model_state.secondary_length_m_by_domain_layer_axis }) |values| {
        if (values.len != axis_values) return error.RootCarbonTotalsDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootCarbonTotalsState;
    }
    inline for (.{ model_state.mobile_carbon_g_c_by_domain_layer, model_state.primary_carbon_total_g_c_by_domain_layer, model_state.primary_length_total_m_by_domain_layer, model_state.secondary_carbon_total_g_c_by_domain_layer, model_state.secondary_length_total_m_by_domain_layer }) |values| {
        if (values.len != domain_layers) return error.RootCarbonTotalsDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootCarbonTotalsState;
    }
    for (model_state.mobile_carbon_g_c_by_domain_layer) |value| if (value < 0) return error.InvalidRootCarbonTotalsState;
    if (model_state.deepest_rooted_layer_index_by_axis.len != axes or model_state.axis_closed_by_domain_axis.len != domain_axes or model_state.plant_deepest_rooted_layer_index.len != 1) return error.RootCarbonTotalsDimensionMismatch;
}

fn axisLayerIndex(inputs: Inputs, domain: usize, layer: usize, axis: usize) usize {
    return (domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
}

/// Exact GROSUB 7301--7334 cleanup and root-total/state publication. The
/// host primary-tip branch contributes RTLGZ/WTRTZ at 7322--7323 and the
/// enclosing active-axis branch contributes the same values again at
/// 7328--7329, preserving the source's branch-sensitive double accumulation.
pub fn apply(model_state: State, workspace: State, inputs: Inputs) !void {
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.RootCarbonTotalsDimensionOverflow;
    const axis_values = std.math.mul(usize, domain_layers, inputs.root_axis_count) catch return error.RootCarbonTotalsDimensionOverflow;
    const domain_axes = std.math.mul(usize, inputs.biological_domain_count, inputs.root_axis_count) catch return error.RootCarbonTotalsDimensionOverflow;
    try validateState(model_state, axis_values, domain_layers, domain_axes, inputs.root_axis_count);
    try validateState(workspace, axis_values, domain_layers, domain_axes, inputs.root_axis_count);
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or inputs.host_domain_index >= inputs.biological_domain_count or inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count or inputs.profile_bottom_layer_index < inputs.planting_layer_index or inputs.profile_bottom_layer_index >= inputs.soil_layer_count or inputs.layer_active.len != inputs.soil_layer_count or inputs.primary_tip_branch_active_by_layer_axis.len != inputs.soil_layer_count * inputs.root_axis_count) return error.RootCarbonTotalsDimensionMismatch;
    for (model_state.deepest_rooted_layer_index_by_axis) |layer| if (layer >= inputs.soil_layer_count) return error.InvalidRootCarbonTotalsState;
    copyState(workspace, model_state);
    @memset(workspace.primary_carbon_total_g_c_by_domain_layer, 0);
    @memset(workspace.primary_length_total_m_by_domain_layer, 0);
    @memset(workspace.secondary_carbon_total_g_c_by_domain_layer, 0);
    @memset(workspace.secondary_length_total_m_by_domain_layer, 0);

    for (0..inputs.biological_domain_count) |domain| for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
        if (!inputs.layer_active[layer]) continue;
        const domain_layer = domain * inputs.soil_layer_count + layer;
        for (0..inputs.root_axis_count) |axis| {
            const value_index = axisLayerIndex(inputs, domain, layer, axis);
            const domain_axis = domain * inputs.root_axis_count + axis;
            const active_axis = layer <= workspace.deepest_rooted_layer_index_by_axis[axis] and !workspace.axis_closed_by_domain_axis[domain_axis];
            if (active_axis) {
                const host_tip = domain == inputs.host_domain_index and inputs.primary_tip_branch_active_by_layer_axis[layer * inputs.root_axis_count + axis];
                if (host_tip) {
                    if (workspace.primary_carbon_g_c_by_domain_layer_axis[value_index] < 0) {
                        workspace.mobile_carbon_g_c_by_domain_layer[domain_layer] += workspace.primary_carbon_g_c_by_domain_layer_axis[value_index];
                        workspace.primary_carbon_g_c_by_domain_layer_axis[value_index] = 0;
                    }
                    if (workspace.secondary_carbon_g_c_by_domain_layer_axis[value_index] < 0) {
                        workspace.mobile_carbon_g_c_by_domain_layer[domain_layer] += workspace.secondary_carbon_g_c_by_domain_layer_axis[value_index];
                        workspace.secondary_carbon_g_c_by_domain_layer_axis[value_index] = 0;
                    }
                    workspace.primary_length_total_m_by_domain_layer[domain_layer] += workspace.primary_length_m_by_domain_layer_axis[value_index];
                    workspace.primary_carbon_total_g_c_by_domain_layer[domain_layer] += workspace.primary_carbon_g_c_by_domain_layer_axis[value_index];
                    workspace.deepest_rooted_layer_index_by_axis[axis] = @min(workspace.deepest_rooted_layer_index_by_axis[axis], inputs.profile_bottom_layer_index);
                    if (layer == workspace.deepest_rooted_layer_index_by_axis[axis]) workspace.axis_closed_by_domain_axis[domain_axis] = true;
                }
                workspace.primary_length_total_m_by_domain_layer[domain_layer] += workspace.primary_length_m_by_domain_layer_axis[value_index];
                workspace.primary_carbon_total_g_c_by_domain_layer[domain_layer] += workspace.primary_carbon_g_c_by_domain_layer_axis[value_index];
            }
            workspace.plant_deepest_rooted_layer_index[0] = @max(workspace.plant_deepest_rooted_layer_index[0], workspace.deepest_rooted_layer_index_by_axis[axis]);
            workspace.secondary_carbon_total_g_c_by_domain_layer[domain_layer] += workspace.secondary_carbon_g_c_by_domain_layer_axis[value_index];
            workspace.secondary_length_total_m_by_domain_layer[domain_layer] += workspace.secondary_length_m_by_domain_layer_axis[value_index];
        }
    };
    for (workspace.mobile_carbon_g_c_by_domain_layer) |value| if (!std.math.isFinite(value) or value < 0) return error.RootCarbonCleanupWouldOverdrawMobilePool;
    try validateState(workspace, axis_values, domain_layers, domain_axes, inputs.root_axis_count);
    copyState(model_state, workspace);
}

fn makeState(values: *[9][2]f64, deepest: *[1]usize, closed: *[1]bool, plant_deepest: *[1]usize) State {
    return .{ .primary_carbon_g_c_by_domain_layer_axis = &values[0], .secondary_carbon_g_c_by_domain_layer_axis = &values[1], .primary_length_m_by_domain_layer_axis = &values[2], .secondary_length_m_by_domain_layer_axis = &values[3], .mobile_carbon_g_c_by_domain_layer = &values[4], .deepest_rooted_layer_index_by_axis = deepest, .axis_closed_by_domain_axis = closed, .plant_deepest_rooted_layer_index = plant_deepest, .primary_carbon_total_g_c_by_domain_layer = &values[5], .primary_length_total_m_by_domain_layer = &values[6], .secondary_carbon_total_g_c_by_domain_layer = &values[7], .secondary_length_total_m_by_domain_layer = &values[8] };
}

test "GROSUB negative root carbon cleanup conserves carbon and preserves double primary total" {
    var values: [9][2]f64 = std.mem.zeroes([9][2]f64);
    values[0] = .{ 1, -0.25 };
    values[1] = .{ 2, -0.5 };
    values[2] = .{ 3, 4 };
    values[3] = .{ 5, 6 };
    values[4] = .{ 2, 2 };
    var deepest = [_]usize{1};
    var closed = [_]bool{false};
    var plant_deepest = [_]usize{0};
    var work: [9][2]f64 = std.mem.zeroes([9][2]f64);
    var wd = [_]usize{0};
    var wc = [_]bool{false};
    var wp = [_]usize{0};
    const carbon_before = values[0][1] + values[1][1] + values[4][1];
    try apply(makeState(&values, &deepest, &closed, &plant_deepest), makeState(&work, &wd, &wc, &wp), .{ .biological_domain_count = 1, .soil_layer_count = 2, .root_axis_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .profile_bottom_layer_index = 1, .host_domain_index = 0, .layer_active = &.{ true, true }, .primary_tip_branch_active_by_layer_axis = &.{ false, true } });
    try std.testing.expectApproxEqAbs(carbon_before, values[0][1] + values[1][1] + values[4][1], 1e-12);
    try std.testing.expectEqual(@as(f64, 0), values[0][1]);
    try std.testing.expectEqual(@as(f64, 0), values[1][1]);
    try std.testing.expectEqual(@as(f64, 0), values[5][1]);
    try std.testing.expectApproxEqAbs(8, values[6][1], 1e-12);
    try std.testing.expectEqual(@as(f64, 0), values[7][1]);
    try std.testing.expectApproxEqAbs(6, values[8][1], 1e-12);
    try std.testing.expect(closed[0]);
    try std.testing.expectEqual(@as(usize, 1), plant_deepest[0]);
}

test "GROSUB negative cleanup rolls back mobile carbon overdraw" {
    var values: [9][2]f64 = std.mem.zeroes([9][2]f64);
    values[0][1] = -2;
    values[4][1] = 1;
    values[2][1] = 1;
    const before = values;
    var deepest = [_]usize{1};
    var closed = [_]bool{false};
    var plant = [_]usize{0};
    var work: [9][2]f64 = std.mem.zeroes([9][2]f64);
    var wd = [_]usize{0};
    var wc = [_]bool{false};
    var wp = [_]usize{0};
    try std.testing.expectError(error.RootCarbonCleanupWouldOverdrawMobilePool, apply(makeState(&values, &deepest, &closed, &plant), makeState(&work, &wd, &wc, &wp), .{ .biological_domain_count = 1, .soil_layer_count = 2, .root_axis_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .profile_bottom_layer_index = 1, .host_domain_index = 0, .layer_active = &.{ true, true }, .primary_tip_branch_active_by_layer_axis = &.{ false, true } }));
    try std.testing.expectEqualDeep(before, values);
}
