const std = @import("std");

pub const State = struct {
    primary_carbon_g_c: []f64,
    primary_nitrogen_g_n: []f64,
    primary_phosphorus_g_p: []f64,
    secondary_carbon_g_c: []f64,
    secondary_nitrogen_g_n: []f64,
    secondary_phosphorus_g_p: []f64,
    primary_length_m: []f64,
    mobile_carbon_g_c: []f64,
    mobile_nitrogen_g_n: []f64,
    mobile_phosphorus_g_p: []f64,
    protein_carbon_g_c: []f64,
    total_root_carbon_g_c: []f64,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    trigger_active_by_domain_layer_axis: []const bool,
    deepest_rooted_layer_index_by_axis: []const usize,
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

fn validateState(state: State, axis_layers: usize, domain_layers: usize) !void {
    inline for (.{ state.primary_carbon_g_c, state.primary_nitrogen_g_n, state.primary_phosphorus_g_p, state.secondary_carbon_g_c, state.secondary_nitrogen_g_n, state.secondary_phosphorus_g_p, state.primary_length_m }) |values| {
        if (values.len != axis_layers) return error.PrimaryRootUpwardTransferDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootUpwardTransferState;
    }
    inline for (.{ state.mobile_carbon_g_c, state.mobile_nitrogen_g_n, state.mobile_phosphorus_g_p, state.protein_carbon_g_c, state.total_root_carbon_g_c }) |values| {
        if (values.len != domain_layers) return error.PrimaryRootUpwardTransferDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootUpwardTransferState;
    }
}

fn validateInputs(inputs: Inputs, axis_layers: usize, domain_layers: usize, domain_axes: usize) !void {
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count or
        inputs.trigger_active_by_domain_layer_axis.len != axis_layers or inputs.deepest_rooted_layer_index_by_axis.len != inputs.root_axis_count or inputs.retracted_depth_m_by_domain_axis.len != domain_axes or
        inputs.layer_bottom_depth_m.len != inputs.soil_layer_count or inputs.layer_thickness_m.len != inputs.soil_layer_count or inputs.axis_primary_sink_strength_m_by_domain_layer_axis.len != axis_layers or inputs.axis_secondary_sink_strength_m_by_domain_layer_axis.len != axis_layers or inputs.total_root_sink_strength_m_by_domain_layer.len != domain_layers)
        return error.PrimaryRootUpwardTransferDimensionMismatch;
    inline for (.{ inputs.seeding_depth_m, inputs.minimum_active_layer_thickness_m, inputs.sink_presence_threshold_m }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootUpwardTransferInput;
    for (inputs.deepest_rooted_layer_index_by_axis) |layer| if (layer < inputs.planting_layer_index or layer >= inputs.soil_layer_count) return error.InvalidPrimaryRootUpwardTransferInput;
    inline for (.{ inputs.retracted_depth_m_by_domain_axis, inputs.layer_bottom_depth_m, inputs.layer_thickness_m, inputs.axis_primary_sink_strength_m_by_domain_layer_axis, inputs.axis_secondary_sink_strength_m_by_domain_layer_axis, inputs.total_root_sink_strength_m_by_domain_layer }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootUpwardTransferInput;
}

fn layerIndex(inputs: Inputs, domain: usize, layer: usize) usize {
    return domain * inputs.soil_layer_count + layer;
}
fn axisLayerIndex(inputs: Inputs, domain: usize, layer: usize, axis: usize) usize {
    return (domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
}

fn transferMass(current: *f64, upper: *f64, fraction: f64) void {
    const amount = fraction * current.*;
    current.* -= amount;
    upper.* += amount;
}

fn transferLayer(state: State, inputs: Inputs, controlling_domain: usize, trigger_layer: usize, layer: usize, axis: usize) !void {
    const controller_axis_layer = axisLayerIndex(inputs, controlling_domain, layer, axis);
    const controller_domain_layer = layerIndex(inputs, controlling_domain, layer);
    const total_sink = inputs.total_root_sink_strength_m_by_domain_layer[controller_domain_layer];
    const fraction = if (total_sink > inputs.sink_presence_threshold_m)
        (inputs.axis_primary_sink_strength_m_by_domain_layer_axis[controller_axis_layer] + inputs.axis_secondary_sink_strength_m_by_domain_layer_axis[controller_axis_layer]) / total_sink
    else
        1;
    if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidPrimaryRootUpwardTransferFraction;

    // GROSUB NN loop: the controlling domain/axis fraction moves every
    // biological domain before the source layer is cleared.
    for (0..inputs.biological_domain_count) |domain| {
        const current_axis = axisLayerIndex(inputs, domain, layer, axis);
        const upper_axis = axisLayerIndex(inputs, domain, layer - 1, axis);
        state.primary_carbon_g_c[upper_axis] += state.primary_carbon_g_c[current_axis];
        state.primary_nitrogen_g_n[upper_axis] += state.primary_nitrogen_g_n[current_axis];
        state.primary_phosphorus_g_p[upper_axis] += state.primary_phosphorus_g_p[current_axis];
        state.secondary_carbon_g_c[upper_axis] += state.secondary_carbon_g_c[current_axis];
        state.secondary_nitrogen_g_n[upper_axis] += state.secondary_nitrogen_g_n[current_axis];
        state.secondary_phosphorus_g_p[upper_axis] += state.secondary_phosphorus_g_p[current_axis];
        state.primary_length_m[upper_axis] += state.primary_length_m[current_axis];
        state.primary_carbon_g_c[current_axis] = 0;
        state.primary_nitrogen_g_n[current_axis] = 0;
        state.primary_phosphorus_g_p[current_axis] = 0;
        state.secondary_carbon_g_c[current_axis] = 0;
        state.secondary_nitrogen_g_n[current_axis] = 0;
        state.secondary_phosphorus_g_p[current_axis] = 0;
        state.primary_length_m[current_axis] = 0;

        const current = layerIndex(inputs, domain, layer);
        const upper = layerIndex(inputs, domain, layer - 1);
        transferMass(&state.mobile_carbon_g_c[current], &state.mobile_carbon_g_c[upper], fraction);
        transferMass(&state.mobile_nitrogen_g_n[current], &state.mobile_nitrogen_g_n[upper], fraction);
        transferMass(&state.mobile_phosphorus_g_p[current], &state.mobile_phosphorus_g_p[upper], fraction);
        // Source line 7163 intentionally reads WSRTL at original L while
        // subtracting at descending LL.
        const protein_source = layerIndex(inputs, domain, trigger_layer);
        const protein_transfer = fraction * state.protein_carbon_g_c[protein_source];
        if (protein_transfer > state.protein_carbon_g_c[current]) return error.PrimaryRootUpwardTransferWouldOverdrawState;
        state.protein_carbon_g_c[current] -= protein_transfer;
        state.protein_carbon_g_c[upper] += protein_transfer;
        transferMass(&state.total_root_carbon_g_c[current], &state.total_root_carbon_g_c[upper], fraction);
    }
}

/// Exact GROSUB 7127--7174 primary/secondary structural and mobile C/N/P
/// withdrawal. Traversal is controlling biological domain (N), rooted layer
/// (L), root axis (NR), descending withdrawal layer (LL), then all domains
/// (NN). NI is one shared plant bound. Masses are g C/N/P and length is m.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.PrimaryRootUpwardTransferDimensionOverflow;
    const axis_layers = std.math.mul(usize, domain_layers, inputs.root_axis_count) catch return error.PrimaryRootUpwardTransferDimensionOverflow;
    const domain_axes = std.math.mul(usize, inputs.biological_domain_count, inputs.root_axis_count) catch return error.PrimaryRootUpwardTransferDimensionOverflow;
    try validateState(state, axis_layers, domain_layers);
    try validateState(workspace, axis_layers, domain_layers);
    try validateInputs(inputs, axis_layers, domain_layers, domain_axes);
    copyState(workspace, state);
    for (0..inputs.biological_domain_count) |domain| for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |trigger_layer| for (0..inputs.root_axis_count) |axis| {
        const trigger = axisLayerIndex(inputs, domain, trigger_layer, axis);
        if (!inputs.trigger_active_by_domain_layer_axis[trigger] or trigger_layer != inputs.deepest_rooted_layer_index_by_axis[axis]) continue;
        var layer = trigger_layer;
        while (layer > inputs.planting_layer_index) : (layer -= 1) {
            const depth = inputs.retracted_depth_m_by_domain_axis[domain * inputs.root_axis_count + axis];
            if (!(inputs.layer_thickness_m[layer - 1] > inputs.minimum_active_layer_thickness_m and (depth < inputs.layer_bottom_depth_m[layer - 1] or depth < inputs.seeding_depth_m))) break;
            try transferLayer(workspace, inputs, domain, trigger_layer, layer, axis);
        }
    };
    try validateState(workspace, axis_layers, domain_layers);
    copyState(state, workspace);
}

fn makeState(v: *[12][4]f64) State {
    return .{ .primary_carbon_g_c = &v[0], .primary_nitrogen_g_n = &v[1], .primary_phosphorus_g_p = &v[2], .secondary_carbon_g_c = &v[3], .secondary_nitrogen_g_n = &v[4], .secondary_phosphorus_g_p = &v[5], .primary_length_m = &v[6], .mobile_carbon_g_c = &v[7], .mobile_nitrogen_g_n = &v[8], .mobile_phosphorus_g_p = &v[9], .protein_carbon_g_c = &v[10], .total_root_carbon_g_c = &v[11] };
}

test "GROSUB upward withdrawal conserves structural mobile C N P and length" {
    var v: [12][4]f64 = std.mem.zeroes([12][4]f64);
    // Two domains, two layers, one axis: deep indexes 1 and 3.
    inline for (0..7) |field| {
        v[field][1] = 2;
        v[field][3] = 3;
    }
    v[7] = .{ 1, 4, 2, 6 };
    v[8] = .{ 0.1, 0.4, 0.2, 0.6 };
    v[9] = .{ 0.01, 0.04, 0.02, 0.06 };
    v[10] = .{ 1, 2, 1, 2 };
    v[11] = .{ 1, 4, 2, 6 };
    var w: [12][4]f64 = std.mem.zeroes([12][4]f64);
    const structural_before = [_]f64{5} ** 7;
    try apply(makeState(&v), makeState(&w), .{ .biological_domain_count = 2, .soil_layer_count = 2, .root_axis_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .trigger_active_by_domain_layer_axis = &.{ false, true, false, false }, .deepest_rooted_layer_index_by_axis = &.{1}, .retracted_depth_m_by_domain_axis = &.{ 0.1, 0.1 }, .seeding_depth_m = 0.15, .layer_bottom_depth_m = &.{ 0.2, 0.4 }, .layer_thickness_m = &.{ 0.2, 0.2 }, .minimum_active_layer_thickness_m = 0.01, .axis_primary_sink_strength_m_by_domain_layer_axis = &.{ 0, 1, 0, 0 }, .axis_secondary_sink_strength_m_by_domain_layer_axis = &.{ 0, 1, 0, 0 }, .total_root_sink_strength_m_by_domain_layer = &.{ 0, 4, 0, 0 }, .sink_presence_threshold_m = 1e-12 });
    for (0..7) |field| try std.testing.expectApproxEqAbs(structural_before[field], v[field][0] + v[field][1] + v[field][2] + v[field][3], 1e-12);
    try std.testing.expectEqual(@as(f64, 0), v[0][1]);
    try std.testing.expectEqual(@as(f64, 0), v[0][3]);
    try std.testing.expectApproxEqAbs(5, v[6][0] + v[6][2], 1e-12);
    try std.testing.expectApproxEqAbs(13, v[7][0] + v[7][1] + v[7][2] + v[7][3], 1e-12);
    try std.testing.expectApproxEqAbs(1.3, v[8][0] + v[8][1] + v[8][2] + v[8][3], 1e-12);
    try std.testing.expectApproxEqAbs(0.13, v[9][0] + v[9][1] + v[9][2] + v[9][3], 1e-12);
}

test "GROSUB upward withdrawal rolls back invalid late sink fraction" {
    var v: [12][4]f64 = std.mem.zeroes([12][4]f64);
    var w: [12][4]f64 = std.mem.zeroes([12][4]f64);
    v[10] = .{ 0, 10, 0, 0.1 };
    const before = v;
    try std.testing.expectError(error.InvalidPrimaryRootUpwardTransferFraction, apply(makeState(&v), makeState(&w), .{ .biological_domain_count = 2, .soil_layer_count = 2, .root_axis_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .trigger_active_by_domain_layer_axis = &.{ false, true, false, false }, .deepest_rooted_layer_index_by_axis = &.{1}, .retracted_depth_m_by_domain_axis = &.{ 0, 0 }, .seeding_depth_m = 0.1, .layer_bottom_depth_m = &.{ 0.2, 0.4 }, .layer_thickness_m = &.{ 0.2, 0.2 }, .minimum_active_layer_thickness_m = 0, .axis_primary_sink_strength_m_by_domain_layer_axis = &.{ 0, 1, 0, 0 }, .axis_secondary_sink_strength_m_by_domain_layer_axis = &.{ 0, 1, 0, 0 }, .total_root_sink_strength_m_by_domain_layer = &.{ 0, 1, 0, 0 }, .sink_presence_threshold_m = 0 }));
    try std.testing.expectEqualDeep(before, v);
}
