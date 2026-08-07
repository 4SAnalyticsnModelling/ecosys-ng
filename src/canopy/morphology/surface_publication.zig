const std = @import("std");

pub const State = struct {
    node_leaf_area_m2: []const f64,
    node_layer_leaf_area_m2: []const f64,
    layer_stalk_area_m2: []const f64,
    node_layer_inclination_leaf_surface_m2: []f64,
    layer_inclination_stalk_surface_m2: []f64,
};

pub const Inputs = struct {
    solar_angle_sine: f64,
    current_branch: usize,
    main_branch: usize,
    first_surface_node: usize,
    end_surface_node: usize,
    node_count: usize,
    layer_count: usize,
    inclination_count: usize,
    azimuth_count: usize,
    leaf_fraction_by_inclination: []const f64,
    side_branch_angle_sine: f64,
};

/// grosub.f lines 4005--4049. At daylight, clears and republishes leaf surfaces
/// in node/layer/inclination order, then clears stalk surfaces and publishes one
/// inclination class. Main stalk is vertical; side branches use their angle.
/// Runtime azimuth count replaces the source uniform factor 0.25.
pub fn publish(state: State, inputs: Inputs) !bool {
    if (!std.math.isFinite(inputs.solar_angle_sine)) return error.NonFiniteCanopySurfaceInput;
    if (inputs.solar_angle_sine <= 0) return false;
    try validate(state, inputs);

    @memset(state.node_layer_inclination_leaf_surface_m2, 0);
    const inverse_azimuth_count = 1.0 / @as(f64, @floatFromInt(inputs.azimuth_count));
    for (inputs.first_surface_node..inputs.end_surface_node) |node| {
        if (state.node_leaf_area_m2[node] <= 0) continue;
        for (0..inputs.layer_count) |layer| {
            const layer_area_m2 = state.node_layer_leaf_area_m2[node * inputs.layer_count + layer];
            for (0..inputs.inclination_count) |inclination| {
                const output = (node * inputs.layer_count + layer) * inputs.inclination_count + inclination;
                state.node_layer_inclination_leaf_surface_m2[output] = @max(0, inputs.leaf_fraction_by_inclination[inclination] * inverse_azimuth_count * layer_area_m2);
            }
        }
    }

    @memset(state.layer_inclination_stalk_surface_m2, 0);
    const stalk_inclination = if (inputs.current_branch == inputs.main_branch)
        inputs.inclination_count - 1
    else blk: {
        const width_radians = (0.5 * std.math.pi) / @as(f64, @floatFromInt(inputs.inclination_count));
        break :blk @min(inputs.inclination_count - 1, @as(usize, @intFromFloat(@floor(std.math.asin(inputs.side_branch_angle_sine) / width_radians))));
    };
    for (0..inputs.layer_count) |layer| {
        const output = layer * inputs.inclination_count + stalk_inclination;
        state.layer_inclination_stalk_surface_m2[output] = inverse_azimuth_count * state.layer_stalk_area_m2[layer];
    }
    return true;
}

fn validate(state: State, inputs: Inputs) !void {
    if (inputs.node_count == 0 or inputs.layer_count == 0 or inputs.inclination_count == 0 or inputs.azimuth_count == 0 or inputs.first_surface_node > inputs.end_surface_node or inputs.end_surface_node > inputs.node_count)
        return error.InvalidCanopySurfaceDimensions;
    if (!std.math.isFinite(inputs.side_branch_angle_sine) or inputs.side_branch_angle_sine < 0 or inputs.side_branch_angle_sine > 1)
        return error.InvalidCanopySurfaceInput;
    if (state.node_leaf_area_m2.len != inputs.node_count or state.layer_stalk_area_m2.len != inputs.layer_count or inputs.leaf_fraction_by_inclination.len != inputs.inclination_count)
        return error.CanopySurfaceDimensionMismatch;
    const node_layers = std.math.mul(usize, inputs.node_count, inputs.layer_count) catch return error.CanopySurfaceDimensionOverflow;
    const node_layer_inclinations = std.math.mul(usize, node_layers, inputs.inclination_count) catch return error.CanopySurfaceDimensionOverflow;
    const layer_inclinations = std.math.mul(usize, inputs.layer_count, inputs.inclination_count) catch return error.CanopySurfaceDimensionOverflow;
    if (state.node_layer_leaf_area_m2.len != node_layers or state.node_layer_inclination_leaf_surface_m2.len != node_layer_inclinations or state.layer_inclination_stalk_surface_m2.len != layer_inclinations)
        return error.CanopySurfaceDimensionMismatch;
    var fraction_total: f64 = 0;
    for (inputs.leaf_fraction_by_inclination) |fraction| {
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidCanopyInclinationFraction;
        fraction_total += fraction;
    }
    if (!std.math.isFinite(fraction_total) or @abs(fraction_total - 1) > 1e-12) return error.InvalidCanopyInclinationFraction;
    inline for (.{ state.node_leaf_area_m2, state.node_layer_leaf_area_m2, state.layer_stalk_area_m2 }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidCanopySurfaceState;
}

test "daylight publishes leaf then vertical main-stalk surfaces" {
    var leaf_surface: [12]f64 = @splat(9);
    var stalk_surface: [6]f64 = @splat(9);
    try std.testing.expect(try publish(.{
        .node_leaf_area_m2 = &.{ 0, 3 },
        .node_layer_leaf_area_m2 = &.{ 0, 0, 1, 2 },
        .layer_stalk_area_m2 = &.{ 4, 8 },
        .node_layer_inclination_leaf_surface_m2 = &leaf_surface,
        .layer_inclination_stalk_surface_m2 = &stalk_surface,
    }, .{ .solar_angle_sine = 0.5, .current_branch = 4, .main_branch = 4, .first_surface_node = 1, .end_surface_node = 2, .node_count = 2, .layer_count = 2, .inclination_count = 3, .azimuth_count = 4, .leaf_fraction_by_inclination = &.{ 0.2, 0.3, 0.5 }, .side_branch_angle_sine = 0 }));
    try std.testing.expectEqual(@as(f64, 0), leaf_surface[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), leaf_surface[6], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), leaf_surface[11], 1e-15);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 1, 0, 0, 2 }, &stalk_surface);
}

test "side branch angle selects runtime inclination class" {
    var leaf_surface: [5]f64 = undefined;
    var stalk_surface: [5]f64 = undefined;
    _ = try publish(.{ .node_leaf_area_m2 = &.{0}, .node_layer_leaf_area_m2 = &.{0}, .layer_stalk_area_m2 = &.{6}, .node_layer_inclination_leaf_surface_m2 = &leaf_surface, .layer_inclination_stalk_surface_m2 = &stalk_surface }, .{ .solar_angle_sine = 1, .current_branch = 2, .main_branch = 0, .first_surface_node = 0, .end_surface_node = 1, .node_count = 1, .layer_count = 1, .inclination_count = 5, .azimuth_count = 3, .leaf_fraction_by_inclination = &.{ 0.2, 0.2, 0.2, 0.2, 0.2 }, .side_branch_angle_sine = std.math.sin(std.math.pi / 4.0) });
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 2, 0, 0 }, &stalk_surface);
}

test "night is strict no-op before topology reads" {
    var leaf_surface = [_]f64{std.math.nan(f64)};
    var stalk_surface = [_]f64{std.math.nan(f64)};
    try std.testing.expect(!try publish(.{ .node_leaf_area_m2 = &.{}, .node_layer_leaf_area_m2 = &.{}, .layer_stalk_area_m2 = &.{}, .node_layer_inclination_leaf_surface_m2 = &leaf_surface, .layer_inclination_stalk_surface_m2 = &stalk_surface }, .{ .solar_angle_sine = 0, .current_branch = 99, .main_branch = 0, .first_surface_node = 9, .end_surface_node = 0, .node_count = 0, .layer_count = 0, .inclination_count = 0, .azimuth_count = 0, .leaf_fraction_by_inclination = &.{}, .side_branch_angle_sine = std.math.nan(f64) }));
    try std.testing.expect(std.math.isNan(leaf_surface[0]));
}

test "runtime node thirty and late invalid state are handled atomically" {
    var leaf_area: [31]f64 = @splat(0);
    var layer_area: [31]f64 = @splat(0);
    var leaf_surface: [62]f64 = @splat(7);
    var stalk_surface: [2]f64 = @splat(7);
    leaf_area[30] = 1;
    layer_area[30] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidCanopySurfaceState, publish(.{ .node_leaf_area_m2 = &leaf_area, .node_layer_leaf_area_m2 = &layer_area, .layer_stalk_area_m2 = &.{1}, .node_layer_inclination_leaf_surface_m2 = &leaf_surface, .layer_inclination_stalk_surface_m2 = &stalk_surface }, .{ .solar_angle_sine = 1, .current_branch = 0, .main_branch = 0, .first_surface_node = 1, .end_surface_node = 31, .node_count = 31, .layer_count = 1, .inclination_count = 2, .azimuth_count = 4, .leaf_fraction_by_inclination = &.{ 0.5, 0.5 }, .side_branch_angle_sine = 0 }));
    try std.testing.expectEqual(@as(f64, 7), leaf_surface[61]);
    try std.testing.expectEqual(@as(f64, 7), stalk_surface[1]);
}
