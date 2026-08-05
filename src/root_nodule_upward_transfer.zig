const std = @import("std");

pub const State = struct {
    structural_carbon_g_c_by_layer: []f64,
    structural_nitrogen_g_n_by_layer: []f64,
    structural_phosphorus_g_p_by_layer: []f64,
    mobile_carbon_g_c_by_layer: []f64,
    mobile_nitrogen_g_n_by_layer: []f64,
    mobile_phosphorus_g_p_by_layer: []f64,
    deepest_rooted_layer_index_by_axis: []usize,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    nitrogen_fixation_type: i32,
    trigger_active_by_domain_layer_axis: []const bool,
    initial_deepest_rooted_layer_index_by_axis: []const usize,
    retracted_depth_m_by_domain_axis: []const f64,
    seeding_depth_m: f64,
    layer_bottom_depth_m: []const f64,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    axis_primary_sink_strength_m_by_domain_layer_axis: []const f64,
    axis_secondary_sink_strength_m_by_domain_layer_axis: []const f64,
    total_root_sink_strength_m_by_domain_layer: []const f64,
    sink_presence_threshold_m: f64,
};

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateState(model_state: State, layers: usize, axes: usize) !void {
    inline for (.{ model_state.structural_carbon_g_c_by_layer, model_state.structural_nitrogen_g_n_by_layer, model_state.structural_phosphorus_g_p_by_layer, model_state.mobile_carbon_g_c_by_layer, model_state.mobile_nitrogen_g_n_by_layer, model_state.mobile_phosphorus_g_p_by_layer }) |values| {
        if (values.len != layers) return error.RootNoduleTransferDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootNoduleTransferState;
    }
    if (model_state.deepest_rooted_layer_index_by_axis.len != axes) return error.RootNoduleTransferDimensionMismatch;
    for (model_state.deepest_rooted_layer_index_by_axis) |layer| if (layer >= layers) return error.InvalidRootNoduleTransferState;
}

fn validateInputs(inputs: Inputs, domain_layers: usize, axis_layers: usize, domain_axes: usize) !void {
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count or
        inputs.trigger_active_by_domain_layer_axis.len != axis_layers or inputs.initial_deepest_rooted_layer_index_by_axis.len != inputs.root_axis_count or inputs.retracted_depth_m_by_domain_axis.len != domain_axes or inputs.layer_bottom_depth_m.len != inputs.soil_layer_count or inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.axis_primary_sink_strength_m_by_domain_layer_axis.len != axis_layers or inputs.axis_secondary_sink_strength_m_by_domain_layer_axis.len != axis_layers or inputs.total_root_sink_strength_m_by_domain_layer.len != domain_layers)
        return error.RootNoduleTransferDimensionMismatch;
    inline for (.{ inputs.seeding_depth_m, inputs.minimum_active_layer_thickness_m, inputs.sink_presence_threshold_m }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootNoduleTransferInput;
    for (inputs.initial_deepest_rooted_layer_index_by_axis) |layer| if (layer < inputs.planting_layer_index or layer >= inputs.soil_layer_count) return error.InvalidRootNoduleTransferInput;
    inline for (.{ inputs.retracted_depth_m_by_domain_axis, inputs.layer_bottom_depth_m, inputs.layer_thickness_m, inputs.axis_primary_sink_strength_m_by_domain_layer_axis, inputs.axis_secondary_sink_strength_m_by_domain_layer_axis, inputs.total_root_sink_strength_m_by_domain_layer }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootNoduleTransferInput;
}

fn layerIndex(inputs: Inputs, domain: usize, layer: usize) usize {
    return domain * inputs.soil_layer_count + layer;
}
fn axisLayerIndex(inputs: Inputs, domain: usize, layer: usize, axis: usize) usize {
    return (domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
}
fn transfer(current: *f64, upper: *f64, fraction: f64) void {
    const amount = fraction * current.*;
    current.* -= amount;
    upper.* += amount;
}

fn transferLayer(model_state: State, inputs: Inputs, domain: usize, layer: usize, axis: usize) !void {
    const axis_layer = axisLayerIndex(inputs, domain, layer, axis);
    const domain_layer = layerIndex(inputs, domain, layer);
    const total_sink = inputs.total_root_sink_strength_m_by_domain_layer[domain_layer];
    const fraction = if (total_sink > inputs.sink_presence_threshold_m)
        (inputs.axis_primary_sink_strength_m_by_domain_layer_axis[axis_layer] + inputs.axis_secondary_sink_strength_m_by_domain_layer_axis[axis_layer]) / total_sink
    else
        1;
    if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidRootNoduleTransferFraction;
    if (inputs.nitrogen_fixation_type >= 1 and inputs.nitrogen_fixation_type <= 3) {
        transfer(&model_state.structural_carbon_g_c_by_layer[layer], &model_state.structural_carbon_g_c_by_layer[layer - 1], fraction);
        transfer(&model_state.structural_nitrogen_g_n_by_layer[layer], &model_state.structural_nitrogen_g_n_by_layer[layer - 1], fraction);
        transfer(&model_state.structural_phosphorus_g_p_by_layer[layer], &model_state.structural_phosphorus_g_p_by_layer[layer - 1], fraction);
        transfer(&model_state.mobile_carbon_g_c_by_layer[layer], &model_state.mobile_carbon_g_c_by_layer[layer - 1], fraction);
        transfer(&model_state.mobile_nitrogen_g_n_by_layer[layer], &model_state.mobile_nitrogen_g_n_by_layer[layer - 1], fraction);
        transfer(&model_state.mobile_phosphorus_g_p_by_layer[layer], &model_state.mobile_phosphorus_g_p_by_layer[layer - 1], fraction);
    }
    model_state.deepest_rooted_layer_index_by_axis[axis] = @max(inputs.planting_layer_index, layer - 1);
}

/// Exact GROSUB 7262--7293 nodule C/N/P upward transfer and NINR update.
/// Traversal is controlling domain (N), shared-NI layer (L), runtime root axis
/// (NR), then descending layer (LL). Nodule pools are plant-wide as in Fortran.
pub fn apply(model_state: State, workspace: State, inputs: Inputs) !void {
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.RootNoduleTransferDimensionOverflow;
    const axis_layers = std.math.mul(usize, domain_layers, inputs.root_axis_count) catch return error.RootNoduleTransferDimensionOverflow;
    const domain_axes = std.math.mul(usize, inputs.biological_domain_count, inputs.root_axis_count) catch return error.RootNoduleTransferDimensionOverflow;
    try validateState(model_state, inputs.soil_layer_count, inputs.root_axis_count);
    try validateState(workspace, inputs.soil_layer_count, inputs.root_axis_count);
    try validateInputs(inputs, domain_layers, axis_layers, domain_axes);
    copyState(workspace, model_state);
    for (0..inputs.biological_domain_count) |domain| for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |trigger_layer| for (0..inputs.root_axis_count) |axis| {
        const trigger = axisLayerIndex(inputs, domain, trigger_layer, axis);
        // NINR is shared mutable plant state. A preceding biological domain may
        // already have withdrawn this axis, so read the staged value here.
        if (!inputs.trigger_active_by_domain_layer_axis[trigger] or trigger_layer != workspace.deepest_rooted_layer_index_by_axis[axis]) continue;
        var layer = trigger_layer;
        while (layer > inputs.planting_layer_index) : (layer -= 1) {
            const depth_m = inputs.retracted_depth_m_by_domain_axis[domain * inputs.root_axis_count + axis];
            if (!(inputs.layer_thickness_m[layer - 1] > inputs.minimum_active_layer_thickness_m and (depth_m < inputs.layer_bottom_depth_m[layer - 1] or depth_m < inputs.seeding_depth_m))) break;
            try transferLayer(workspace, inputs, domain, layer, axis);
        }
    };
    try validateState(workspace, inputs.soil_layer_count, inputs.root_axis_count);
    copyState(model_state, workspace);
}

fn makeState(values: *[6][3]f64, deepest: *[1]usize) State {
    return .{ .structural_carbon_g_c_by_layer = &values[0], .structural_nitrogen_g_n_by_layer = &values[1], .structural_phosphorus_g_p_by_layer = &values[2], .mobile_carbon_g_c_by_layer = &values[3], .mobile_nitrogen_g_n_by_layer = &values[4], .mobile_phosphorus_g_p_by_layer = &values[5], .deepest_rooted_layer_index_by_axis = deepest };
}

test "GROSUB legume nodule upward transfer conserves C N P and publishes NINR" {
    var values = [6][3]f64{ .{ 1, 2, 4 }, .{ 0.1, 0.2, 0.4 }, .{ 0.01, 0.02, 0.04 }, .{ 3, 4, 8 }, .{ 0.3, 0.4, 0.8 }, .{ 0.03, 0.04, 0.08 } };
    var deepest = [_]usize{2};
    var work: [6][3]f64 = std.mem.zeroes([6][3]f64);
    var work_deepest = [_]usize{0};
    var before: [6]f64 = undefined;
    for (0..6) |element| before[element] = values[element][0] + values[element][1] + values[element][2];
    try apply(makeState(&values, &deepest), makeState(&work, &work_deepest), .{ .biological_domain_count = 1, .soil_layer_count = 3, .root_axis_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 2, .nitrogen_fixation_type = 2, .trigger_active_by_domain_layer_axis = &.{ false, false, true }, .initial_deepest_rooted_layer_index_by_axis = &.{2}, .retracted_depth_m_by_domain_axis = &.{0.05}, .seeding_depth_m = 0.1, .layer_bottom_depth_m = &.{ 0.2, 0.4, 0.6 }, .layer_thickness_m = &.{ 0.2, 0.2, 0.2 }, .minimum_active_layer_thickness_m = 0.01, .axis_primary_sink_strength_m_by_domain_layer_axis = &.{ 0, 1, 1 }, .axis_secondary_sink_strength_m_by_domain_layer_axis = &.{ 0, 1, 1 }, .total_root_sink_strength_m_by_domain_layer = &.{ 0, 4, 4 }, .sink_presence_threshold_m = 0 });
    for (0..6) |element| try std.testing.expectApproxEqAbs(before[element], values[element][0] + values[element][1] + values[element][2], 1e-12);
    try std.testing.expectEqual(@as(usize, 0), deepest[0]);
}

test "GROSUB nonlegume skips nodule mass but still updates deepest axis" {
    var values = [6][3]f64{ .{ 1, 2, 4 }, .{ 1, 2, 4 }, .{ 1, 2, 4 }, .{ 1, 2, 4 }, .{ 1, 2, 4 }, .{ 1, 2, 4 } };
    const before = values;
    var deepest = [_]usize{2};
    var work: [6][3]f64 = std.mem.zeroes([6][3]f64);
    var wd = [_]usize{0};
    try apply(makeState(&values, &deepest), makeState(&work, &wd), .{ .biological_domain_count = 1, .soil_layer_count = 3, .root_axis_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 2, .nitrogen_fixation_type = 0, .trigger_active_by_domain_layer_axis = &.{ false, false, true }, .initial_deepest_rooted_layer_index_by_axis = &.{2}, .retracted_depth_m_by_domain_axis = &.{0}, .seeding_depth_m = 0.1, .layer_bottom_depth_m = &.{ 0.2, 0.4, 0.6 }, .layer_thickness_m = &.{ 0.2, 0.2, 0.2 }, .minimum_active_layer_thickness_m = 0, .axis_primary_sink_strength_m_by_domain_layer_axis = &.{ 0, 0, 0 }, .axis_secondary_sink_strength_m_by_domain_layer_axis = &.{ 0, 0, 0 }, .total_root_sink_strength_m_by_domain_layer = &.{ 0, 0, 0 }, .sink_presence_threshold_m = 0 });
    try std.testing.expectEqualDeep(before, values);
    try std.testing.expectEqual(@as(usize, 0), deepest[0]);
}

test "GROSUB nodule transfer is atomic on invalid late fraction" {
    var values = [6][3]f64{ .{ 1, 2, 4 }, .{ 1, 2, 4 }, .{ 1, 2, 4 }, .{ 1, 2, 4 }, .{ 1, 2, 4 }, .{ 1, 2, 4 } };
    const before = values;
    var deepest = [_]usize{2};
    var work: [6][3]f64 = std.mem.zeroes([6][3]f64);
    var wd = [_]usize{0};
    try std.testing.expectError(error.InvalidRootNoduleTransferFraction, apply(makeState(&values, &deepest), makeState(&work, &wd), .{ .biological_domain_count = 1, .soil_layer_count = 3, .root_axis_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 2, .nitrogen_fixation_type = 1, .trigger_active_by_domain_layer_axis = &.{ false, false, true }, .initial_deepest_rooted_layer_index_by_axis = &.{2}, .retracted_depth_m_by_domain_axis = &.{0}, .seeding_depth_m = 0.1, .layer_bottom_depth_m = &.{ 0.2, 0.4, 0.6 }, .layer_thickness_m = &.{ 0.2, 0.2, 0.2 }, .minimum_active_layer_thickness_m = 0, .axis_primary_sink_strength_m_by_domain_layer_axis = &.{ 0, 1, 1 }, .axis_secondary_sink_strength_m_by_domain_layer_axis = &.{ 0, 0, 1 }, .total_root_sink_strength_m_by_domain_layer = &.{ 0, 0.5, 2 }, .sink_presence_threshold_m = 0 }));
    try std.testing.expectEqualDeep(before, values);
}
