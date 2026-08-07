const std = @import("std");

pub const State = struct {
    /// Source `ICHK1(N,NR)`, flattened `[domain][axis]`.
    tip_processed_by_domain_axis: []bool,
    /// Source `RTN1(N,L)`, flattened `[domain][layer]`.
    primary_axis_count_by_domain_layer: []f64,
    /// True only where the source enters line 6446 for this axis.
    tip_active_by_domain_layer_axis: []bool,
    primary_sink_fraction_by_domain_layer_axis: []f64,
};

pub const Workspace = State;

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_index: usize,
    /// Source `NI(NZ,NY,NX)`: one plant-wide inclusive rooted bound.
    deepest_rooted_layer_index: usize,
    bulk_density_megagrams_per_m3_by_layer: []const f64,
    layer_top_depth_m: []const f64,
    layer_bottom_depth_m: []const f64,
    soil_profile_surface_offset_m: f64,
    profile_bottom_layer_index: usize,
    root_depth_from_surface_m_by_domain_axis: []const f64,
    primary_axis_count_multiplier: f64,
    primary_sink_strength_m_by_domain_layer_axis: []const f64,
    total_sink_strength_m_by_domain_layer: []const f64,
    /// Source `ZEROP(NZ,NY,NX)`: one plant-wide sink threshold.
    negligible_sink_strength_m: f64,
};

fn copyState(target: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        @memcpy(@field(target, field.name), @field(source, field.name));
}

fn validateShape(state: State, domain_axes: usize, domain_layers: usize, values: usize) !void {
    if (state.tip_processed_by_domain_axis.len != domain_axes or
        state.primary_axis_count_by_domain_layer.len != domain_layers or
        state.tip_active_by_domain_layer_axis.len != values or
        state.primary_sink_fraction_by_domain_layer_axis.len != values)
        return error.PrimaryRootTipDimensionMismatch;
}

/// Exact grosub.f lines 6438--6466 N,L,NR primary-tip transaction.
/// Depths and sink strengths are m; bulk density is Mg m-3; FRTN is
/// dimensionless. Domain zero is the zero-based mapping of source `N.EQ.1`.
pub fn apply(state: State, workspace: Workspace, inputs: Inputs) !void {
    const domain_axes = std.math.mul(usize, inputs.biological_domain_count, inputs.root_axis_count) catch return error.PrimaryRootTipDimensionOverflow;
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.PrimaryRootTipDimensionOverflow;
    const values = std.math.mul(usize, domain_layers, inputs.root_axis_count) catch return error.PrimaryRootTipDimensionOverflow;
    try validateShape(state, domain_axes, domain_layers, values);
    try validateShape(workspace, domain_axes, domain_layers, values);
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or
        inputs.bulk_density_megagrams_per_m3_by_layer.len != inputs.soil_layer_count or
        inputs.layer_top_depth_m.len != inputs.soil_layer_count or inputs.layer_bottom_depth_m.len != inputs.soil_layer_count or
        inputs.root_depth_from_surface_m_by_domain_axis.len != domain_axes or
        inputs.primary_sink_strength_m_by_domain_layer_axis.len != values or
        inputs.total_sink_strength_m_by_domain_layer.len != domain_layers)
        return error.PrimaryRootTipDimensionMismatch;
    if (inputs.profile_bottom_layer_index >= inputs.soil_layer_count or
        inputs.planting_layer_index > inputs.deepest_rooted_layer_index or
        inputs.deepest_rooted_layer_index >= inputs.soil_layer_count)
        return error.InvalidPrimaryRootTipLayerRange;
    inline for (.{ inputs.soil_profile_surface_offset_m, inputs.primary_axis_count_multiplier, inputs.negligible_sink_strength_m }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootTipInput;
    inline for (.{ inputs.bulk_density_megagrams_per_m3_by_layer, inputs.layer_top_depth_m, inputs.layer_bottom_depth_m, inputs.root_depth_from_surface_m_by_domain_axis, inputs.total_sink_strength_m_by_domain_layer }) |slice| for (slice) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootTipInput;
    for (0..inputs.soil_layer_count) |layer|
        if (inputs.layer_bottom_depth_m[layer] < inputs.layer_top_depth_m[layer]) return error.InvalidPrimaryRootTipLayerGeometry;
    for (state.primary_axis_count_by_domain_layer) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootTipState;

    copyState(workspace, state);
    for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
        const domain_layer = layer;
        for (0..inputs.root_axis_count) |axis| {
            const domain_axis = axis;
            const index = domain_layer * inputs.root_axis_count + axis;
            const adjusted_depth_m = if (inputs.bulk_density_megagrams_per_m3_by_layer[layer] > 0)
                inputs.root_depth_from_surface_m_by_domain_axis[domain_axis] - inputs.soil_profile_surface_offset_m
            else
                inputs.root_depth_from_surface_m_by_domain_axis[domain_axis];
            if (!std.math.isFinite(adjusted_depth_m)) return error.InvalidPrimaryRootTipInput;
            if (adjusted_depth_m > inputs.layer_top_depth_m[layer] and !workspace.tip_processed_by_domain_axis[domain_axis]) {
                workspace.primary_axis_count_by_domain_layer[domain_layer] += inputs.primary_axis_count_multiplier;
                if (!std.math.isFinite(workspace.primary_axis_count_by_domain_layer[domain_layer])) return error.NonFinitePrimaryRootAxisCount;
                if (adjusted_depth_m <= inputs.layer_bottom_depth_m[layer] or layer == inputs.profile_bottom_layer_index) {
                    const strength_m = inputs.primary_sink_strength_m_by_domain_layer_axis[index];
                    const total_m = inputs.total_sink_strength_m_by_domain_layer[domain_layer];
                    if (!std.math.isFinite(strength_m) or strength_m < 0 or strength_m > total_m) return error.InvalidPrimaryRootSinkStrength;
                    workspace.tip_processed_by_domain_axis[domain_axis] = true;
                    workspace.tip_active_by_domain_layer_axis[index] = true;
                    workspace.primary_sink_fraction_by_domain_layer_axis[index] = if (total_m > inputs.negligible_sink_strength_m) strength_m / total_m else 1;
                }
            }
        }
    }
    copyState(state, workspace);
}

test "GROSUB primary tips traverse layers once and publish source FRTN" {
    var processed = [_]bool{ false, false, true, true };
    var counts = [_]f64{ 0, 0, 7, 7, 7, 7 };
    var active = [_]bool{false} ** 12;
    var fractions = [_]f64{-1} ** 12;
    var w_processed = [_]bool{false} ** 4;
    var w_counts = [_]f64{0} ** 6;
    var w_active = [_]bool{false} ** 12;
    var w_fractions = [_]f64{0} ** 12;
    const state = State{ .tip_processed_by_domain_axis = &processed, .primary_axis_count_by_domain_layer = &counts, .tip_active_by_domain_layer_axis = &active, .primary_sink_fraction_by_domain_layer_axis = &fractions };
    const workspace = Workspace{ .tip_processed_by_domain_axis = &w_processed, .primary_axis_count_by_domain_layer = &w_counts, .tip_active_by_domain_layer_axis = &w_active, .primary_sink_fraction_by_domain_layer_axis = &w_fractions };
    try apply(state, workspace, .{
        .biological_domain_count = 2,
        .soil_layer_count = 3,
        .root_axis_count = 2,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 2,
        .bulk_density_megagrams_per_m3_by_layer = &.{ 1.2, 1.3, 1.4 },
        .layer_top_depth_m = &.{ 0, 0.2, 0.4 },
        .layer_bottom_depth_m = &.{ 0.2, 0.4, 0.6 },
        .soil_profile_surface_offset_m = 0,
        .profile_bottom_layer_index = 2,
        .root_depth_from_surface_m_by_domain_axis = &.{ 0.15, 0.45, 99, 99 },
        .primary_axis_count_multiplier = 3,
        .primary_sink_strength_m_by_domain_layer_axis = &.{ 2, 2, 4, 4, 6, 6, 0, 0, 0, 0, 0, 0 },
        .total_sink_strength_m_by_domain_layer = &.{ 10, 20, 30, 1, 1, 1 },
        .negligible_sink_strength_m = 1e-12,
    });
    try std.testing.expect(active[0]);
    try std.testing.expect(active[5]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), fractions[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), fractions[5], 1e-15);
    try std.testing.expectEqual(@as(f64, 6), counts[0]);
    try std.testing.expectEqual(@as(f64, 3), counts[1]);
    try std.testing.expectEqual(@as(f64, 10), counts[2]);
    try std.testing.expectEqual(@as(f64, 7), counts[3]);
}

test "GROSUB invalid late primary tip rolls back processed and RTN1 state" {
    var processed = [_]bool{ false, false };
    var counts = [_]f64{ 1, 2 };
    var active = [_]bool{ false, false, false, false };
    var fractions = [_]f64{ -1, -1, -1, -1 };
    var w_processed = [_]bool{false} ** 2;
    var w_counts = [_]f64{0} ** 2;
    var w_active = [_]bool{false} ** 4;
    var w_fractions = [_]f64{0} ** 4;
    const before_processed = processed;
    const before_counts = counts;
    try std.testing.expectError(error.InvalidPrimaryRootSinkStrength, apply(
        .{ .tip_processed_by_domain_axis = &processed, .primary_axis_count_by_domain_layer = &counts, .tip_active_by_domain_layer_axis = &active, .primary_sink_fraction_by_domain_layer_axis = &fractions },
        .{ .tip_processed_by_domain_axis = &w_processed, .primary_axis_count_by_domain_layer = &w_counts, .tip_active_by_domain_layer_axis = &w_active, .primary_sink_fraction_by_domain_layer_axis = &w_fractions },
        .{ .biological_domain_count = 1, .soil_layer_count = 2, .root_axis_count = 2, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .bulk_density_megagrams_per_m3_by_layer = &.{ 1, 1 }, .layer_top_depth_m = &.{ 0, 0.2 }, .layer_bottom_depth_m = &.{ 0.2, 0.4 }, .soil_profile_surface_offset_m = 0, .profile_bottom_layer_index = 1, .root_depth_from_surface_m_by_domain_axis = &.{ 0.1, 0.3 }, .primary_axis_count_multiplier = 1, .primary_sink_strength_m_by_domain_layer_axis = &.{ 1, 1, 1, std.math.nan(f64) }, .total_sink_strength_m_by_domain_layer = &.{ 2, 2 }, .negligible_sink_strength_m = 1e-12 },
    ));
    try std.testing.expectEqualDeep(before_processed, processed);
    try std.testing.expectEqualDeep(before_counts, counts);
}
