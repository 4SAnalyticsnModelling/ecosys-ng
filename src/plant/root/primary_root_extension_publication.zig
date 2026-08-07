const std = @import("std");

pub const State = struct {
    axis_total_carbon_g_c: []f64,
    axis_total_nitrogen_g_n: []f64,
    axis_total_phosphorus_g_p: []f64,
    axis_depth_m: []f64,
    layer_carbon_g_c: []f64,
    layer_nitrogen_g_n: []f64,
    layer_phosphorus_g_p: []f64,
    layer_length_m: []f64,
    layer_protein_carbon_g_c: []f64,
    layer_total_root_carbon_g_c: []f64,
    mobile_carbon_g_c: []f64,
    mobile_nitrogen_g_n: []f64,
    mobile_phosphorus_g_p: []f64,
    total_water_potential_megapascal: []f64,
    osmotic_water_potential_megapascal: []f64,
    turgor_water_potential_megapascal: []f64,
    primary_radius_m: []f64,
    deepest_rooted_layer_index_by_axis: []usize,
    /// Source transient `RTLGZ`, retained for downstream morphology totals.
    primary_root_length_accumulator_m: []f64,
    /// Source transient `WTRTZ`, retained for downstream morphology totals.
    primary_root_carbon_accumulator_g_c: []f64,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    minimum_rooted_layer_index: usize,
    profile_bottom_layer_index: usize,
    /// Source `L1` selected earlier for each `[domain][layer]`.
    next_lower_layer_index_by_domain_layer: []const usize,
    tip_active_by_domain_layer_axis: []const bool,
    gross_growth_g_c_by_domain_layer_axis: []const f64,
    net_growth_g_c_by_domain_layer_axis: []const f64,
    net_growth_g_n_by_domain_layer_axis: []const f64,
    net_growth_g_p_by_domain_layer_axis: []const f64,
    specific_length_m_per_g_c_by_domain: []const f64,
    population_count: f64,
    extension_water_response_by_domain_layer: []const f64,
    seeding_depth_m: f64,
    layer_bottom_depth_m: []const f64,
    layer_thickness_m: []const f64,
    extension_presence_threshold_m: f64,
    root_carbon_presence_threshold_g_c: f64,
    protein_carbon_per_nitrogen_g_c_per_g_n: f64,
    protein_carbon_per_phosphorus_g_c_per_g_p: f64,
    axis_sink_fraction_by_domain_layer_axis: []const f64,
};

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateState(state: State, domain_axes: usize, axis_layers: usize, domain_layers: usize, axes: usize) !void {
    inline for (.{ state.axis_total_carbon_g_c, state.axis_total_nitrogen_g_n, state.axis_total_phosphorus_g_p, state.axis_depth_m }) |values| if (values.len != domain_axes) return error.PrimaryRootPublicationDimensionMismatch;
    inline for (.{ state.layer_carbon_g_c, state.layer_nitrogen_g_n, state.layer_phosphorus_g_p, state.layer_length_m }) |values| if (values.len != axis_layers) return error.PrimaryRootPublicationDimensionMismatch;
    inline for (.{ state.layer_protein_carbon_g_c, state.layer_total_root_carbon_g_c, state.mobile_carbon_g_c, state.mobile_nitrogen_g_n, state.mobile_phosphorus_g_p, state.total_water_potential_megapascal, state.osmotic_water_potential_megapascal, state.turgor_water_potential_megapascal, state.primary_radius_m }) |values| if (values.len != domain_layers) return error.PrimaryRootPublicationDimensionMismatch;
    if (state.deepest_rooted_layer_index_by_axis.len != axes) return error.PrimaryRootPublicationDimensionMismatch;
    if (state.primary_root_length_accumulator_m.len != 1 or state.primary_root_carbon_accumulator_g_c.len != 1) return error.PrimaryRootPublicationDimensionMismatch;
    inline for (@typeInfo(State).@"struct".fields) |field| switch (field.type) {
        []f64 => for (@field(state, field.name)) |value| if (!std.math.isFinite(value)) return error.NonFinitePrimaryRootPublicationState,
        else => {},
    };
}

fn validateInputs(inputs: Inputs, axis_layers: usize, domain_layers: usize) !void {
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count or inputs.minimum_rooted_layer_index >= inputs.soil_layer_count or inputs.profile_bottom_layer_index >= inputs.soil_layer_count or
        inputs.next_lower_layer_index_by_domain_layer.len != domain_layers or
        inputs.tip_active_by_domain_layer_axis.len != axis_layers or inputs.gross_growth_g_c_by_domain_layer_axis.len != axis_layers or inputs.net_growth_g_c_by_domain_layer_axis.len != axis_layers or inputs.net_growth_g_n_by_domain_layer_axis.len != axis_layers or inputs.net_growth_g_p_by_domain_layer_axis.len != axis_layers or inputs.axis_sink_fraction_by_domain_layer_axis.len != axis_layers or
        inputs.specific_length_m_per_g_c_by_domain.len != inputs.biological_domain_count or inputs.extension_water_response_by_domain_layer.len != domain_layers or inputs.layer_bottom_depth_m.len != inputs.soil_layer_count or inputs.layer_thickness_m.len != inputs.soil_layer_count)
        return error.PrimaryRootPublicationDimensionMismatch;
    inline for (.{ inputs.population_count, inputs.seeding_depth_m, inputs.extension_presence_threshold_m, inputs.root_carbon_presence_threshold_g_c, inputs.protein_carbon_per_nitrogen_g_c_per_g_n, inputs.protein_carbon_per_phosphorus_g_c_per_g_p }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootPublicationInput;
    if (inputs.population_count == 0) return error.InvalidPrimaryRootPublicationInput;
    inline for (.{ inputs.gross_growth_g_c_by_domain_layer_axis, inputs.net_growth_g_c_by_domain_layer_axis, inputs.net_growth_g_n_by_domain_layer_axis, inputs.net_growth_g_p_by_domain_layer_axis, inputs.axis_sink_fraction_by_domain_layer_axis }) |values| for (values) |value| if (!std.math.isFinite(value)) return error.InvalidPrimaryRootPublicationInput;
    for (inputs.gross_growth_g_c_by_domain_layer_axis) |value| if (value < 0) return error.InvalidPrimaryRootPublicationInput;
    for (inputs.axis_sink_fraction_by_domain_layer_axis) |value| if (value < 0 or value > 1) return error.InvalidPrimaryRootPublicationInput;
    inline for (.{ inputs.specific_length_m_per_g_c_by_domain, inputs.extension_water_response_by_domain_layer, inputs.layer_bottom_depth_m, inputs.layer_thickness_m }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootPublicationInput;
    for (inputs.extension_water_response_by_domain_layer) |value| if (value > 1) return error.InvalidPrimaryRootPublicationInput;
    for (1..inputs.biological_domain_count) |domain| for (0..inputs.soil_layer_count) |layer| for (0..inputs.root_axis_count) |axis|
        if (inputs.tip_active_by_domain_layer_axis[(domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis])
            return error.PrimaryRootTipOutsideHostDomain;
    for (0..inputs.biological_domain_count) |domain| for (0..inputs.soil_layer_count) |layer| {
        const lower = inputs.next_lower_layer_index_by_domain_layer[domain * inputs.soil_layer_count + layer];
        if (lower >= inputs.soil_layer_count or (layer < inputs.profile_bottom_layer_index and lower <= layer))
            return error.InvalidPrimaryRootPublicationLowerLayer;
    };
}

fn axisIndex(inputs: Inputs, domain: usize, axis: usize) usize {
    return domain * inputs.root_axis_count + axis;
}
fn layerIndex(inputs: Inputs, domain: usize, layer: usize) usize {
    return domain * inputs.soil_layer_count + layer;
}
fn axisLayerIndex(inputs: Inputs, domain: usize, layer: usize, axis: usize) usize {
    return (domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
}

fn applyAxis(state: State, inputs: Inputs, domain: usize, layer: usize, axis: usize) !void {
    const source = axisLayerIndex(inputs, domain, layer, axis);
    if (!inputs.tip_active_by_domain_layer_axis[source]) return;
    const aggregate = axisIndex(inputs, domain, axis);
    const domain_layer = layerIndex(inputs, domain, layer);
    const net_c = inputs.net_growth_g_c_by_domain_layer_axis[source];
    const net_n = inputs.net_growth_g_n_by_domain_layer_axis[source];
    const net_p = inputs.net_growth_g_p_by_domain_layer_axis[source];
    const lower_layer = inputs.next_lower_layer_index_by_domain_layer[domain_layer];
    var extension_m = inputs.gross_growth_g_c_by_domain_layer_axis[source] * inputs.specific_length_m_per_g_c_by_domain[domain] / inputs.population_count * inputs.extension_water_response_by_domain_layer[domain_layer];
    if (net_c < 0 and state.axis_total_carbon_g_c[aggregate] > inputs.root_carbon_presence_threshold_g_c) extension_m += net_c * (state.axis_depth_m[aggregate] - inputs.seeding_depth_m) / state.axis_total_carbon_g_c[aggregate];
    if (layer < inputs.profile_bottom_layer_index) extension_m = @min(inputs.layer_thickness_m[lower_layer], extension_m);
    if (!std.math.isFinite(extension_m)) return error.NonFinitePrimaryRootPublication;
    var crosses = false;
    if (extension_m > inputs.extension_presence_threshold_m and layer < inputs.profile_bottom_layer_index) {
        const fraction_current = std.math.clamp((inputs.layer_bottom_depth_m[layer] - state.axis_depth_m[aggregate]) / extension_m, 0, 1);
        crosses = fraction_current < 1;
    }
    const target_layer = if (crosses) lower_layer else layer;
    const target = axisLayerIndex(inputs, domain, target_layer, axis);
    const target_domain_layer = layerIndex(inputs, domain, target_layer);
    state.axis_total_carbon_g_c[aggregate] += net_c;
    state.axis_total_nitrogen_g_n[aggregate] += net_n;
    state.axis_total_phosphorus_g_p[aggregate] += net_p;
    state.axis_depth_m[aggregate] += extension_m;
    state.layer_carbon_g_c[target] += net_c;
    state.layer_nitrogen_g_n[target] += net_n;
    state.layer_phosphorus_g_p[target] += net_p;
    state.layer_length_m[target] += extension_m;
    state.layer_protein_carbon_g_c[target_domain_layer] += @min(inputs.protein_carbon_per_nitrogen_g_c_per_g_n * state.layer_nitrogen_g_n[target], inputs.protein_carbon_per_phosphorus_g_c_per_g_p * state.layer_phosphorus_g_p[target]);
    if (crosses) {
        state.layer_total_root_carbon_g_c[target_domain_layer] += state.layer_carbon_g_c[target];
        state.primary_radius_m[target_domain_layer] = state.primary_radius_m[domain_layer];
        const fraction = inputs.axis_sink_fraction_by_domain_layer_axis[source];
        const carbon_transfer = fraction * state.mobile_carbon_g_c[domain_layer];
        const nitrogen_transfer = fraction * state.mobile_nitrogen_g_n[domain_layer];
        const phosphorus_transfer = fraction * state.mobile_phosphorus_g_p[domain_layer];
        state.mobile_carbon_g_c[domain_layer] -= carbon_transfer;
        state.mobile_nitrogen_g_n[domain_layer] -= nitrogen_transfer;
        state.mobile_phosphorus_g_p[domain_layer] -= phosphorus_transfer;
        state.mobile_carbon_g_c[target_domain_layer] += carbon_transfer;
        state.mobile_nitrogen_g_n[target_domain_layer] += nitrogen_transfer;
        state.mobile_phosphorus_g_p[target_domain_layer] += phosphorus_transfer;
        state.total_water_potential_megapascal[target_domain_layer] = state.total_water_potential_megapascal[domain_layer];
        state.osmotic_water_potential_megapascal[target_domain_layer] = state.osmotic_water_potential_megapascal[domain_layer];
        state.turgor_water_potential_megapascal[target_domain_layer] = state.turgor_water_potential_megapascal[domain_layer];
        state.primary_root_length_accumulator_m[0] += state.layer_length_m[target];
        state.primary_root_carbon_accumulator_g_c[0] += state.layer_carbon_g_c[target];
        if (!std.math.isFinite(state.primary_root_length_accumulator_m[0]) or !std.math.isFinite(state.primary_root_carbon_accumulator_g_c[0])) return error.NonFinitePrimaryRootPublication;
        state.deepest_rooted_layer_index_by_axis[axis] = @max(inputs.minimum_rooted_layer_index, layer + 1);
    }
    inline for (.{ state.axis_total_carbon_g_c[aggregate], state.axis_total_nitrogen_g_n[aggregate], state.axis_total_phosphorus_g_p[aggregate], state.axis_depth_m[aggregate], state.layer_carbon_g_c[target], state.layer_nitrogen_g_n[target], state.layer_phosphorus_g_p[target], state.layer_length_m[target] }) |value| if (!std.math.isFinite(value) or value < 0) return error.PrimaryRootPublicationWouldOverdrawState;
}

/// Exact GROSUB 6978--7080 primary-root extension and state publication.
/// Traversal is biological domain (N), shared plant rooted layers (L), then
/// runtime root axes (NR). C/N/P masses are g, lengths/depths are m, and water
/// potentials are MPa. Publication is atomic through caller-owned workspace.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.PrimaryRootPublicationDimensionOverflow;
    const axis_layers = std.math.mul(usize, domain_layers, inputs.root_axis_count) catch return error.PrimaryRootPublicationDimensionOverflow;
    const domain_axes = std.math.mul(usize, inputs.biological_domain_count, inputs.root_axis_count) catch return error.PrimaryRootPublicationDimensionOverflow;
    try validateState(state, domain_axes, axis_layers, domain_layers, inputs.root_axis_count);
    try validateState(workspace, domain_axes, axis_layers, domain_layers, inputs.root_axis_count);
    try validateInputs(inputs, axis_layers, domain_layers);
    copyState(workspace, state);
    for (0..inputs.biological_domain_count) |domain| for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| for (0..inputs.root_axis_count) |axis| try applyAxis(workspace, inputs, domain, layer, axis);
    copyState(state, workspace);
}

fn makeState(comptime layers: usize, v: *[19][layers]f64, deepest: *[1]usize) State {
    return .{ .axis_total_carbon_g_c = v[0][0..1], .axis_total_nitrogen_g_n = v[1][0..1], .axis_total_phosphorus_g_p = v[2][0..1], .axis_depth_m = v[3][0..1], .layer_carbon_g_c = &v[4], .layer_nitrogen_g_n = &v[5], .layer_phosphorus_g_p = &v[6], .layer_length_m = &v[7], .layer_protein_carbon_g_c = &v[8], .layer_total_root_carbon_g_c = &v[9], .mobile_carbon_g_c = &v[10], .mobile_nitrogen_g_n = &v[11], .mobile_phosphorus_g_p = &v[12], .total_water_potential_megapascal = &v[13], .osmotic_water_potential_megapascal = &v[14], .turgor_water_potential_megapascal = &v[15], .primary_radius_m = &v[16], .deepest_rooted_layer_index_by_axis = deepest, .primary_root_length_accumulator_m = v[17][0..1], .primary_root_carbon_accumulator_g_c = v[18][0..1] };
}

test "GROSUB primary extension crossing conserves C N P and mobile pools" {
    var v: [19][2]f64 = std.mem.zeroes([19][2]f64);
    v[0][0] = 10;
    v[1][0] = 1;
    v[2][0] = 0.1;
    v[3][0] = 0.19;
    v[4][0] = 5;
    v[5][0] = 0.5;
    v[6][0] = 0.05;
    v[7][0] = 1;
    v[10] = .{ 4, 1 };
    v[11] = .{ 0.4, 0.1 };
    v[12] = .{ 0.04, 0.01 };
    v[13][0] = -0.2;
    v[14][0] = -0.4;
    v[15][0] = 0.2;
    v[16][0] = 0.001;
    var deepest = [_]usize{0};
    var w: [19][2]f64 = std.mem.zeroes([19][2]f64);
    var wd = [_]usize{0};
    const before_mobile = v[10][0] + v[10][1];
    try apply(makeState(2, &v, &deepest), makeState(2, &w, &wd), .{ .biological_domain_count = 1, .soil_layer_count = 2, .root_axis_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 0, .minimum_rooted_layer_index = 0, .profile_bottom_layer_index = 1, .next_lower_layer_index_by_domain_layer = &.{ 1, 1 }, .tip_active_by_domain_layer_axis = &.{ true, false }, .gross_growth_g_c_by_domain_layer_axis = &.{ 1, 0 }, .net_growth_g_c_by_domain_layer_axis = &.{ 1, 0 }, .net_growth_g_n_by_domain_layer_axis = &.{ 0.1, 0 }, .net_growth_g_p_by_domain_layer_axis = &.{ 0.01, 0 }, .specific_length_m_per_g_c_by_domain = &.{1}, .population_count = 1, .extension_water_response_by_domain_layer = &.{ 1, 1 }, .seeding_depth_m = 0, .layer_bottom_depth_m = &.{ 0.2, 0.4 }, .layer_thickness_m = &.{ 0.2, 0.2 }, .extension_presence_threshold_m = 1e-12, .root_carbon_presence_threshold_g_c = 1e-12, .protein_carbon_per_nitrogen_g_c_per_g_n = 2, .protein_carbon_per_phosphorus_g_c_per_g_p = 20, .axis_sink_fraction_by_domain_layer_axis = &.{ 0.25, 0 } });
    try std.testing.expectEqual(@as(usize, 1), deepest[0]);
    try std.testing.expectApproxEqAbs(before_mobile, v[10][0] + v[10][1], 1e-12);
    try std.testing.expectApproxEqAbs(0.5, v[11][0] + v[11][1], 1e-12);
    try std.testing.expectApproxEqAbs(0.05, v[12][0] + v[12][1], 1e-12);
    try std.testing.expectApproxEqAbs(11, v[0][0], 1e-12);
    try std.testing.expectApproxEqAbs(1.1, v[1][0], 1e-12);
    try std.testing.expectApproxEqAbs(0.11, v[2][0], 1e-12);
    try std.testing.expectApproxEqAbs(1, v[4][1], 1e-12);
    try std.testing.expectApproxEqAbs(0.2, v[17][0], 1e-12);
    try std.testing.expectApproxEqAbs(1, v[18][0], 1e-12);
}

test "GROSUB primary extension publication is atomic on late invalid growth" {
    var v: [19][2]f64 = std.mem.zeroes([19][2]f64);
    var deepest = [_]usize{0};
    v[0][0] = 1;
    v[1][0] = 1;
    v[2][0] = 1;
    var w: [19][2]f64 = std.mem.zeroes([19][2]f64);
    var wd = [_]usize{0};
    const before = v;
    try std.testing.expectError(error.InvalidPrimaryRootPublicationInput, apply(makeState(2, &v, &deepest), makeState(2, &w, &wd), .{ .biological_domain_count = 1, .soil_layer_count = 2, .root_axis_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .minimum_rooted_layer_index = 0, .profile_bottom_layer_index = 1, .next_lower_layer_index_by_domain_layer = &.{ 1, 1 }, .tip_active_by_domain_layer_axis = &.{ true, true }, .gross_growth_g_c_by_domain_layer_axis = &.{ 0, std.math.nan(f64) }, .net_growth_g_c_by_domain_layer_axis = &.{ 0, 0 }, .net_growth_g_n_by_domain_layer_axis = &.{ 0, 0 }, .net_growth_g_p_by_domain_layer_axis = &.{ 0, 0 }, .specific_length_m_per_g_c_by_domain = &.{1}, .population_count = 1, .extension_water_response_by_domain_layer = &.{ 1, 1 }, .seeding_depth_m = 0, .layer_bottom_depth_m = &.{ 0.2, 0.4 }, .layer_thickness_m = &.{ 0.2, 0.2 }, .extension_presence_threshold_m = 0, .root_carbon_presence_threshold_g_c = 0, .protein_carbon_per_nitrogen_g_c_per_g_n = 1, .protein_carbon_per_phosphorus_g_c_per_g_p = 1, .axis_sink_fraction_by_domain_layer_axis = &.{ 0, 0 } }));
    try std.testing.expectEqualDeep(before, v);
}

test "GROSUB primary extension uses selected nonadjacent lower layer but source NINR L plus one" {
    var v: [19][3]f64 = std.mem.zeroes([19][3]f64);
    v[0][0] = 10;
    v[1][0] = 1;
    v[2][0] = 0.1;
    v[3][0] = 0.19;
    v[10][0] = 4;
    v[11][0] = 0.4;
    v[12][0] = 0.04;
    v[16][0] = 0.001;
    var deepest = [_]usize{0};
    var w: [19][3]f64 = std.mem.zeroes([19][3]f64);
    var wd = [_]usize{0};
    try apply(makeState(3, &v, &deepest), makeState(3, &w, &wd), .{
        .biological_domain_count = 1,
        .soil_layer_count = 3,
        .root_axis_count = 1,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 0,
        .minimum_rooted_layer_index = 0,
        .profile_bottom_layer_index = 2,
        .next_lower_layer_index_by_domain_layer = &.{ 2, 2, 2 },
        .tip_active_by_domain_layer_axis = &.{ true, false, false },
        .gross_growth_g_c_by_domain_layer_axis = &.{ 0.3, 0, 0 },
        .net_growth_g_c_by_domain_layer_axis = &.{ 1, 0, 0 },
        .net_growth_g_n_by_domain_layer_axis = &.{ 0.1, 0, 0 },
        .net_growth_g_p_by_domain_layer_axis = &.{ 0.01, 0, 0 },
        .specific_length_m_per_g_c_by_domain = &.{1},
        .population_count = 1,
        .extension_water_response_by_domain_layer = &.{ 1, 1, 1 },
        .seeding_depth_m = 0,
        .layer_bottom_depth_m = &.{ 0.2, 0.2005, 0.5 },
        .layer_thickness_m = &.{ 0.2, 0.0005, 0.3 },
        .extension_presence_threshold_m = 1e-12,
        .root_carbon_presence_threshold_g_c = 1e-12,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 2,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 20,
        .axis_sink_fraction_by_domain_layer_axis = &.{ 0.25, 0, 0 },
    });
    try std.testing.expectEqual(@as(f64, 0), v[4][1]);
    try std.testing.expectEqual(@as(f64, 1), v[4][2]);
    try std.testing.expectEqual(@as(usize, 1), deepest[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), v[17][0], 1e-15);
    try std.testing.expectEqual(@as(f64, 1), v[18][0]);
}
