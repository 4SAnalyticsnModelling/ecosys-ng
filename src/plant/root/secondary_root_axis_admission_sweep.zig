const std = @import("std");

pub const State = struct {
    /// Flattened `[domain][layer][axis]`; inactive entries remain untouched.
    active_by_domain_layer_axis: []bool,
    sink_fraction_by_domain_layer_axis: []f64,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    /// Source `NINR(NR,NZ)`, flattened `[domain][axis]`.
    deepest_secondary_root_layer_by_domain_axis: []const usize,
    /// Source `NRX(N,NR) != 0`, flattened `[domain][axis]`.
    axis_inactive_by_domain_axis: []const bool,
    /// Source `RTSK2(N,L,NR)`, m.
    secondary_sink_strength_m_by_domain_layer_axis: []const f64,
    /// Source `RLNT(N,L)`, m; includes primary and secondary sink strength.
    total_sink_strength_m_by_domain_layer: []const f64,
    negligible_sink_strength_m_by_domain: []const f64,
};

/// Exact grosub.f lines 6038--6055 N,L,NR secondary-axis admission sweep.
/// Sink strengths are m and FRTN is dimensionless.
pub fn apply(state: State, inputs: Inputs) !void {
    const domain_axis_count = std.math.mul(usize, inputs.biological_domain_count, inputs.root_axis_count) catch
        return error.SecondaryRootAdmissionDimensionOverflow;
    const domain_layer_count = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch
        return error.SecondaryRootAdmissionDimensionOverflow;
    const value_count = std.math.mul(usize, domain_layer_count, inputs.root_axis_count) catch
        return error.SecondaryRootAdmissionDimensionOverflow;
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or
        inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.deepest_secondary_root_layer_by_domain_axis.len != domain_axis_count or
        inputs.axis_inactive_by_domain_axis.len != domain_axis_count or
        inputs.secondary_sink_strength_m_by_domain_layer_axis.len != value_count or
        inputs.total_sink_strength_m_by_domain_layer.len != domain_layer_count or
        inputs.negligible_sink_strength_m_by_domain.len != inputs.biological_domain_count or
        state.active_by_domain_layer_axis.len != value_count or
        state.sink_fraction_by_domain_layer_axis.len != value_count)
        return error.SecondaryRootAdmissionDimensionMismatch;
    if (!std.math.isFinite(inputs.minimum_active_layer_thickness_m) or inputs.minimum_active_layer_thickness_m < 0)
        return error.InvalidSecondaryRootAdmissionThreshold;
    for (inputs.layer_thickness_m) |thickness_m|
        if (!std.math.isFinite(thickness_m) or thickness_m < 0)
            return error.InvalidSecondaryRootAdmissionLayerThickness;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count)
        return error.InvalidSecondaryRootAdmissionLayerRange;
    for (inputs.deepest_secondary_root_layer_by_domain_axis) |deepest_layer|
        if (deepest_layer >= inputs.soil_layer_count)
            return error.InvalidSecondaryRootAxisDepth;
    for (inputs.negligible_sink_strength_m_by_domain) |threshold|
        if (!std.math.isFinite(threshold) or threshold < 0)
            return error.InvalidSecondaryRootAdmissionThreshold;

    // Full source-order preflight before any state mutation.
    for (0..inputs.biological_domain_count) |domain| {
        const deepest_layer = inputs.deepest_rooted_layer_index;
        for (inputs.planting_layer_index..deepest_layer + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const domain_layer = domain * inputs.soil_layer_count + layer;
            const total_m = inputs.total_sink_strength_m_by_domain_layer[domain_layer];
            if (!std.math.isFinite(total_m) or total_m < 0)
                return error.InvalidSecondaryRootSinkStrength;
            for (0..inputs.root_axis_count) |axis| {
                const domain_axis = domain * inputs.root_axis_count + axis;
                if (layer > inputs.deepest_secondary_root_layer_by_domain_axis[domain_axis] or
                    inputs.axis_inactive_by_domain_axis[domain_axis]) continue;
                const index = domain_layer * inputs.root_axis_count + axis;
                const strength_m = inputs.secondary_sink_strength_m_by_domain_layer_axis[index];
                if (!std.math.isFinite(strength_m) or strength_m < 0 or strength_m > total_m)
                    return error.InvalidSecondaryRootSinkStrength;
                const fraction = if (total_m > inputs.negligible_sink_strength_m_by_domain[domain])
                    strength_m / total_m
                else
                    1.0;
                if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
                    return error.InvalidSecondaryRootSinkFraction;
            }
        }
    }

    for (0..inputs.biological_domain_count) |domain| {
        const deepest_layer = inputs.deepest_rooted_layer_index;
        for (inputs.planting_layer_index..deepest_layer + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const domain_layer = domain * inputs.soil_layer_count + layer;
            const total_m = inputs.total_sink_strength_m_by_domain_layer[domain_layer];
            for (0..inputs.root_axis_count) |axis| {
                const domain_axis = domain * inputs.root_axis_count + axis;
                if (layer > inputs.deepest_secondary_root_layer_by_domain_axis[domain_axis] or
                    inputs.axis_inactive_by_domain_axis[domain_axis]) continue;
                const index = domain_layer * inputs.root_axis_count + axis;
                state.active_by_domain_layer_axis[index] = true;
                state.sink_fraction_by_domain_layer_axis[index] =
                    if (total_m > inputs.negligible_sink_strength_m_by_domain[domain])
                        inputs.secondary_sink_strength_m_by_domain_layer_axis[index] / total_m
                    else
                        1.0;
            }
        }
    }
}

test "GROSUB N L NR admission preserves NINR NRX and FRTN" {
    const domains = 2;
    const layers = 3;
    const axes = 3;
    var active = [_]bool{false} ** (domains * layers * axes);
    var fractions = [_]f64{-1} ** (domains * layers * axes);
    const strengths = [_]f64{
        1, 2, 3, 4, 5,  6,  7,  8,  9,
        2, 4, 6, 8, 10, 12, 14, 16, 18,
    };
    try apply(.{ .active_by_domain_layer_axis = &active, .sink_fraction_by_domain_layer_axis = &fractions }, .{
        .biological_domain_count = domains,
        .soil_layer_count = layers,
        .root_axis_count = axes,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 2,
        .layer_thickness_m = &.{ 0.2, 0.2, 0.2 },
        .minimum_active_layer_thickness_m = 0.001,
        .deepest_secondary_root_layer_by_domain_axis = &.{ 0, 2, 2, 2, 1, 2 },
        .axis_inactive_by_domain_axis = &.{ false, false, true, false, false, false },
        .secondary_sink_strength_m_by_domain_layer_axis = &strengths,
        .total_sink_strength_m_by_domain_layer = &.{ 10, 20, 30, 20, 40, 60 },
        .negligible_sink_strength_m_by_domain = &.{ 1e-12, 1e-12 },
    });

    try std.testing.expectApproxEqAbs(@as(f64, 0.1), fractions[0], 1e-15);
    try std.testing.expect(!active[3]); // domain 0, layer 1, axis 0 exceeds NINR.
    try std.testing.expectEqual(@as(f64, -1), fractions[5]); // NRX suppresses axis 2.
    try std.testing.expectApproxEqAbs(@as(f64, 14.0 / 60.0), fractions[15], 1e-15);
    try std.testing.expectApproxEqAbs(strengths[15], fractions[15] * 60.0, 1e-14);
}

test "GROSUB FRTN zero-total fallback remains one" {
    var active = [_]bool{false};
    var fractions = [_]f64{-1};
    try apply(.{ .active_by_domain_layer_axis = &active, .sink_fraction_by_domain_layer_axis = &fractions }, .{
        .biological_domain_count = 1,
        .soil_layer_count = 1,
        .root_axis_count = 1,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 0,
        .layer_thickness_m = &.{0.2},
        .minimum_active_layer_thickness_m = 0.001,
        .deepest_secondary_root_layer_by_domain_axis = &.{0},
        .axis_inactive_by_domain_axis = &.{false},
        .secondary_sink_strength_m_by_domain_layer_axis = &.{0},
        .total_sink_strength_m_by_domain_layer = &.{0},
        .negligible_sink_strength_m_by_domain = &.{1e-12},
    });
    try std.testing.expect(active[0]);
    try std.testing.expectEqual(@as(f64, 1), fractions[0]);
}

test "GROSUB secondary admission fails atomically on invalid late axis" {
    var active = [_]bool{true} ** 4;
    var fractions = [_]f64{0.25} ** 4;
    const before_active = active;
    const before_fractions = fractions;
    try std.testing.expectError(error.InvalidSecondaryRootSinkStrength, apply(
        .{ .active_by_domain_layer_axis = &active, .sink_fraction_by_domain_layer_axis = &fractions },
        .{
            .biological_domain_count = 2,
            .soil_layer_count = 1,
            .root_axis_count = 2,
            .planting_layer_index = 0,
            .deepest_rooted_layer_index = 0,
            .layer_thickness_m = &.{0.2},
            .minimum_active_layer_thickness_m = 0.001,
            .deepest_secondary_root_layer_by_domain_axis = &.{ 0, 0, 0, 0 },
            .axis_inactive_by_domain_axis = &.{ false, false, false, false },
            .secondary_sink_strength_m_by_domain_layer_axis = &.{ 1, 2, 3, std.math.nan(f64) },
            .total_sink_strength_m_by_domain_layer = &.{ 4, 8 },
            .negligible_sink_strength_m_by_domain = &.{ 1e-12, 1e-12 },
        },
    ));
    try std.testing.expectEqualDeep(before_active, active);
    try std.testing.expectEqualDeep(before_fractions, fractions);
}
