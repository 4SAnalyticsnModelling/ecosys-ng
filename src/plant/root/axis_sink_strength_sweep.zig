const std = @import("std");
const metabolism = @import("plant_root_metabolism.zig");

pub const State = struct {
    sink_strength_by_domain_layer_axis: []metabolism.RootAxisSinkStrength,
    total_sink_strength_m_by_domain: []f64,
    layer_sink_strength_m_by_domain_layer: []f64,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_index: usize,
    /// Plant-wide inclusive `NI(NZ)` mapped to a zero-based layer.
    deepest_rooted_layer_index: usize,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    /// Flattened [domain][layer][axis] in exact GROSUB N,L,NR order.
    geometry_by_domain_layer_axis: []const metabolism.SourceOrderRootAxisSinkInputs,
    parameters: metabolism.SecondaryRootParameters,
};

/// Exact grosub.f lines 5865--5953 runtime root-axis sink-strength sweep.
/// This owns traversal and accumulation while reusing the verified source-
/// order scalar geometry kernel. Units of RTSK1, RTSK2, RTNT, and RLNT are m.
pub fn apply(state: State, inputs: Inputs) !void {
    const domain_layer_count = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch
        return error.RootAxisSinkDimensionOverflow;
    const axis_value_count = std.math.mul(usize, domain_layer_count, inputs.root_axis_count) catch
        return error.RootAxisSinkDimensionOverflow;
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or
        inputs.root_axis_count == 0 or inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.geometry_by_domain_layer_axis.len != axis_value_count or
        state.sink_strength_by_domain_layer_axis.len != axis_value_count or
        state.total_sink_strength_m_by_domain.len != inputs.biological_domain_count or
        state.layer_sink_strength_m_by_domain_layer.len != domain_layer_count)
        return error.RootAxisSinkDimensionMismatch;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count)
        return error.InvalidRootedLayerRange;
    if (!std.math.isFinite(inputs.minimum_active_layer_thickness_m) or
        inputs.minimum_active_layer_thickness_m < 0)
        return error.InvalidRootAxisSinkSweepInput;
    for (inputs.layer_thickness_m) |thickness_m|
        if (!std.math.isFinite(thickness_m) or thickness_m < 0)
            return error.InvalidRootLayerThickness;
    try inputs.parameters.validate();

    for (state.total_sink_strength_m_by_domain) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidRootAxisSinkAccumulator;
    for (state.layer_sink_strength_m_by_domain_layer) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidRootAxisSinkAccumulator;

    // Full preflight with the exact N outer, L middle, NR inner additions.
    for (0..inputs.biological_domain_count) |domain| {
        var next_domain_total = state.total_sink_strength_m_by_domain[domain];
        const deepest_layer = inputs.deepest_rooted_layer_index;
        for (inputs.planting_layer_index..deepest_layer + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const domain_layer = domain * inputs.soil_layer_count + layer;
            var next_layer_total = state.layer_sink_strength_m_by_domain_layer[domain_layer];
            for (0..inputs.root_axis_count) |axis| {
                const index = domain_layer * inputs.root_axis_count + axis;
                const geometry = inputs.geometry_by_domain_layer_axis[index];
                if (geometry.primary_biological_domain != (domain == 0))
                    return error.RootAxisBiologicalDomainMismatch;
                const strength = try metabolism.sourceOrderRootAxisSinkStrength(inputs.parameters, geometry);
                next_domain_total += strength.primary_m;
                next_layer_total += strength.primary_m;
                next_domain_total += strength.secondary_m;
                next_layer_total += strength.secondary_m;
                inline for (.{ next_domain_total, next_layer_total }) |value|
                    if (!std.math.isFinite(value)) return error.NonFiniteRootAxisSinkAccumulator;
            }
        }
    }

    for (0..inputs.biological_domain_count) |domain| {
        const deepest_layer = inputs.deepest_rooted_layer_index;
        for (inputs.planting_layer_index..deepest_layer + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const domain_layer = domain * inputs.soil_layer_count + layer;
            for (0..inputs.root_axis_count) |axis| {
                const index = domain_layer * inputs.root_axis_count + axis;
                const strength = try metabolism.sourceOrderRootAxisSinkStrength(inputs.parameters, inputs.geometry_by_domain_layer_axis[index]);
                state.sink_strength_by_domain_layer_axis[index] = strength;
                state.total_sink_strength_m_by_domain[domain] += strength.primary_m;
                state.layer_sink_strength_m_by_domain_layer[domain_layer] += strength.primary_m;
                state.total_sink_strength_m_by_domain[domain] += strength.secondary_m;
                state.layer_sink_strength_m_by_domain_layer[domain_layer] += strength.secondary_m;
            }
        }
    }
}

fn testGeometry(primary_domain: bool, layer_top_m: f64, root_depth_m: f64, secondary_count: f64) metabolism.SourceOrderRootAxisSinkInputs {
    return .{
        .root_profile_type = 1,
        .primary_axis_count_multiplier = 2,
        .primary_root_radius_m = 0.01,
        .primary_root_depth_from_surface_m = root_depth_m,
        .layer_top_depth_m = layer_top_m,
        .layer_thickness_m = 0.2,
        .secondary_root_origin_offset_m = 0,
        .seeding_depth_m = 0.05,
        .hypocotyledon_height_m = 0,
        .canopy_height_m = 1,
        .secondary_axis_count = secondary_count,
        .secondary_root_radius_m = 0.001,
        .average_secondary_root_length_m = 0.1,
        .negligible_sink_m = 1e-12,
        .primary_biological_domain = primary_domain,
    };
}

test "GROSUB N L NR sweep accumulates axis values into layer and domain totals" {
    const domains = 2;
    const layers = 2;
    const axes = 2;
    const geometry_values = [_]metabolism.SourceOrderRootAxisSinkInputs{
        testGeometry(true, 0, 0.1, 2),    testGeometry(true, 0, 0.1, 3),
        testGeometry(true, 0.2, 0.3, 2),  testGeometry(true, 0.2, 0.3, 3),
        testGeometry(false, 0, 0.1, 2),   testGeometry(false, 0, 0.1, 3),
        testGeometry(false, 0.2, 0.3, 2), testGeometry(false, 0.2, 0.3, 3),
    };
    var strengths: [domains * layers * axes]metabolism.RootAxisSinkStrength = undefined;
    var domain_totals = [_]f64{0} ** domains;
    var layer_totals = [_]f64{0} ** (domains * layers);
    try apply(.{ .sink_strength_by_domain_layer_axis = &strengths, .total_sink_strength_m_by_domain = &domain_totals, .layer_sink_strength_m_by_domain_layer = &layer_totals }, .{ .biological_domain_count = domains, .soil_layer_count = layers, .root_axis_count = axes, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .layer_thickness_m = &.{ 0.2, 0.2 }, .minimum_active_layer_thickness_m = 0.001, .geometry_by_domain_layer_axis = &geometry_values, .parameters = metabolism.compatibilitySecondaryRootParameters() });
    for (0..domains) |domain| {
        var expected: f64 = 0;
        for (layer_totals[domain * layers .. (domain + 1) * layers]) |value| expected += value;
        try std.testing.expectApproxEqAbs(expected, domain_totals[domain], 1e-15);
    }
    try std.testing.expect(strengths[0].primary_m > 0);
    try std.testing.expectEqual(@as(f64, 0), strengths[4].primary_m);
}

test "GROSUB NI is one plant-wide bound shared by biological domains" {
    const domains = 2;
    const layers = 3;
    const axes = 2;
    var geometry_values = [_]metabolism.SourceOrderRootAxisSinkInputs{testGeometry(true, 0, 0.5, 2)} ** (domains * layers * axes);
    for (0..domains) |domain| for (0..layers) |layer| for (0..axes) |axis| {
        const index = (domain * layers + layer) * axes + axis;
        geometry_values[index].primary_biological_domain = domain == 0;
        geometry_values[index].layer_top_depth_m = @as(f64, @floatFromInt(layer)) * 0.2;
    };
    var strengths = [_]metabolism.RootAxisSinkStrength{.{ .primary_m = -1, .secondary_m = -1 }} ** (domains * layers * axes);
    var domain_totals = [_]f64{0} ** domains;
    var layer_totals = [_]f64{0} ** (domains * layers);

    try apply(
        .{
            .sink_strength_by_domain_layer_axis = &strengths,
            .total_sink_strength_m_by_domain = &domain_totals,
            .layer_sink_strength_m_by_domain_layer = &layer_totals,
        },
        .{
            .biological_domain_count = domains,
            .soil_layer_count = layers,
            .root_axis_count = axes,
            .planting_layer_index = 0,
            .deepest_rooted_layer_index = 2,
            .layer_thickness_m = &.{ 0.2, 0.2, 0.2 },
            .minimum_active_layer_thickness_m = 0.001,
            .geometry_by_domain_layer_axis = &geometry_values,
            .parameters = metabolism.compatibilitySecondaryRootParameters(),
        },
    );

    try std.testing.expect(layer_totals[1] > 0);
    try std.testing.expect(layer_totals[2] > 0);
    try std.testing.expect(strengths[(0 * layers + 1) * axes].secondary_m > 0);
    try std.testing.expect(layer_totals[3] > 0);
    try std.testing.expect(layer_totals[4] > 0);
    try std.testing.expect(layer_totals[5] > 0);
}

test "GROSUB invalid plant-wide NI fails before mutation" {
    var strengths = [_]metabolism.RootAxisSinkStrength{.{ .primary_m = 7, .secondary_m = 8 }} ** 4;
    var domain_totals = [_]f64{ 3, 4 };
    var layer_totals = [_]f64{ 5, 6 };
    const before_strengths = strengths;
    const before_domains = domain_totals;
    const before_layers = layer_totals;
    const geometry_values = [_]metabolism.SourceOrderRootAxisSinkInputs{
        testGeometry(true, 0, 0.1, 2),
        testGeometry(true, 0, 0.1, 2),
        testGeometry(false, 0, 0.1, 2),
        testGeometry(false, 0, 0.1, 2),
    };
    try std.testing.expectError(error.InvalidRootedLayerRange, apply(
        .{ .sink_strength_by_domain_layer_axis = &strengths, .total_sink_strength_m_by_domain = &domain_totals, .layer_sink_strength_m_by_domain_layer = &layer_totals },
        .{ .biological_domain_count = 2, .soil_layer_count = 1, .root_axis_count = 2, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .layer_thickness_m = &.{0.2}, .minimum_active_layer_thickness_m = 0.001, .geometry_by_domain_layer_axis = &geometry_values, .parameters = metabolism.compatibilitySecondaryRootParameters() },
    ));
    try std.testing.expectEqualDeep(before_strengths, strengths);
    try std.testing.expectEqualDeep(before_domains, domain_totals);
    try std.testing.expectEqualDeep(before_layers, layer_totals);
}

test "GROSUB strict thickness gate leaves inactive axis values unread and unchanged" {
    var strengths = [_]metabolism.RootAxisSinkStrength{.{ .primary_m = std.math.nan(f64), .secondary_m = std.math.nan(f64) }};
    var domain_totals = [_]f64{0};
    var layer_totals = [_]f64{0};
    var invalid = testGeometry(true, 0, 0.1, 2);
    invalid.primary_root_radius_m = std.math.nan(f64);
    try apply(.{ .sink_strength_by_domain_layer_axis = &strengths, .total_sink_strength_m_by_domain = &domain_totals, .layer_sink_strength_m_by_domain_layer = &layer_totals }, .{ .biological_domain_count = 1, .soil_layer_count = 1, .root_axis_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 0, .layer_thickness_m = &.{0.001}, .minimum_active_layer_thickness_m = 0.001, .geometry_by_domain_layer_axis = &.{invalid}, .parameters = metabolism.compatibilitySecondaryRootParameters() });
    try std.testing.expect(std.math.isNan(strengths[0].primary_m));
}

test "GROSUB runtime root-axis sweep is atomic on invalid late axis" {
    const domains = 3;
    const layers = 7;
    const axes = 11;
    var geometry_values = [_]metabolism.SourceOrderRootAxisSinkInputs{testGeometry(true, 0, 0.5, 2)} ** (domains * layers * axes);
    for (0..domains) |domain| for (0..layers) |layer| for (0..axes) |axis| {
        const index = (domain * layers + layer) * axes + axis;
        geometry_values[index].primary_biological_domain = domain == 0;
        geometry_values[index].layer_top_depth_m = @as(f64, @floatFromInt(layer)) * 0.2;
    };
    var strengths = [_]metabolism.RootAxisSinkStrength{.{ .primary_m = 0, .secondary_m = 0 }} ** (domains * layers * axes);
    var domain_totals = [_]f64{0} ** domains;
    var layer_totals = [_]f64{0} ** (domains * layers);
    const thickness = [_]f64{0.2} ** layers;
    const state = State{ .sink_strength_by_domain_layer_axis = &strengths, .total_sink_strength_m_by_domain = &domain_totals, .layer_sink_strength_m_by_domain_layer = &layer_totals };
    const inputs = Inputs{ .biological_domain_count = domains, .soil_layer_count = layers, .root_axis_count = axes, .planting_layer_index = 0, .deepest_rooted_layer_index = 6, .layer_thickness_m = &thickness, .minimum_active_layer_thickness_m = 0.001, .geometry_by_domain_layer_axis = &geometry_values, .parameters = metabolism.compatibilitySecondaryRootParameters() };
    try apply(state, inputs);
    const before_strengths = strengths;
    const before_domains = domain_totals;
    geometry_values[geometry_values.len - 1].average_secondary_root_length_m = 0;
    try std.testing.expectError(error.InvalidRootAxisSinkInput, apply(state, inputs));
    try std.testing.expectEqualDeep(before_strengths, strengths);
    try std.testing.expectEqualDeep(before_domains, domain_totals);
}

test "GROSUB sink totals retain exact N outer L middle NR inner addition order" {
    const domains = 2;
    const layers = 2;
    const axes = 3;
    var geometry = [_]metabolism.SourceOrderRootAxisSinkInputs{testGeometry(true, 0, 0.3, 1)} ** (domains * layers * axes);
    for (0..domains) |domain| for (0..layers) |layer| for (0..axes) |axis| {
        const index = (domain * layers + layer) * axes + axis;
        geometry[index].primary_biological_domain = domain == 0;
        geometry[index].layer_top_depth_m = @as(f64, @floatFromInt(layer)) * 0.2;
        geometry[index].secondary_axis_count = switch (index % 4) {
            0 => 1.0e16,
            1 => 1,
            2 => 1.0e-8,
            else => 3,
        };
    };
    var strengths: [domains * layers * axes]metabolism.RootAxisSinkStrength = undefined;
    var domain_totals = [_]f64{0} ** domains;
    var layer_totals = [_]f64{0} ** (domains * layers);
    const parameters = metabolism.compatibilitySecondaryRootParameters();
    try apply(.{ .sink_strength_by_domain_layer_axis = &strengths, .total_sink_strength_m_by_domain = &domain_totals, .layer_sink_strength_m_by_domain_layer = &layer_totals }, .{ .biological_domain_count = domains, .soil_layer_count = layers, .root_axis_count = axes, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .layer_thickness_m = &.{ 0.2, 0.2 }, .minimum_active_layer_thickness_m = 0, .geometry_by_domain_layer_axis = &geometry, .parameters = parameters });

    var expected_domain = [_]f64{0} ** domains;
    var expected_layer = [_]f64{0} ** (domains * layers);
    for (0..domains) |domain| for (0..layers) |layer| for (0..axes) |axis| {
        const index = (domain * layers + layer) * axes + axis;
        const value = try metabolism.sourceOrderRootAxisSinkStrength(parameters, geometry[index]);
        expected_domain[domain] += value.primary_m;
        expected_layer[domain * layers + layer] += value.primary_m;
        expected_domain[domain] += value.secondary_m;
        expected_layer[domain * layers + layer] += value.secondary_m;
    };
    try std.testing.expectEqualDeep(expected_domain, domain_totals);
    try std.testing.expectEqualDeep(expected_layer, layer_totals);
}
