const std = @import("std");

pub const WorkingRadiation = struct {
    shortwave_megajoules_m2_h: f64,
    par_umol_m2_s: f64,
};

pub const Boundaries = struct {
    diffuse_transmittance: []f64,
    forward_scattered_shortwave_megajoules_m2_h: []f64,
    forward_scattered_par_umol_m2_s: []f64,
    backscattered_shortwave_megajoules_m2_h: []const f64,
    backscattered_par_umol_m2_s: []const f64,
};

/// HOUR1 lines 1681--1685. Initializes bottom boundary index zero in exact
/// source assignment order.
pub fn initialize(working: *WorkingRadiation, boundaries: Boundaries) !void {
    try validateBoundaryDimensions(boundaries);
    working.shortwave_megajoules_m2_h = 0.0;
    working.par_umol_m2_s = 0.0;
    boundaries.diffuse_transmittance[0] = 1.0;
    boundaries.forward_scattered_shortwave_megajoules_m2_h[0] = 0.0;
    boundaries.forward_scattered_par_umol_m2_s[0] = 0.0;
}

/// HOUR1 lines 1686--1692 for one zero-based layer in ascending traversal.
/// Returns false for a snow- or water-submerged layer.
pub fn admitLayer(
    layer: usize,
    layer_bottom_height_m: []const f64,
    snow_depth_m: f64,
    surface_water_ice_depth_m: f64,
    depth_tolerance_m: f64,
    working: *WorkingRadiation,
    boundaries: Boundaries,
) !bool {
    try validateLayer(
        layer,
        layer_bottom_height_m,
        snow_depth_m,
        surface_water_ice_depth_m,
        depth_tolerance_m,
        working.*,
        boundaries,
    );
    if (layer_bottom_height_m[layer] < snow_depth_m - depth_tolerance_m or
        layer_bottom_height_m[layer] <
            surface_water_ice_depth_m - depth_tolerance_m)
        return false;
    const below_boundary = layer;
    const current_boundary = layer + 1;
    working.shortwave_megajoules_m2_h =
        working.shortwave_megajoules_m2_h *
        boundaries.diffuse_transmittance[below_boundary] +
        boundaries.forward_scattered_shortwave_megajoules_m2_h[below_boundary] +
        boundaries.backscattered_shortwave_megajoules_m2_h[below_boundary];
    working.par_umol_m2_s =
        working.par_umol_m2_s *
        boundaries.diffuse_transmittance[below_boundary] +
        boundaries.forward_scattered_par_umol_m2_s[below_boundary] +
        boundaries.backscattered_par_umol_m2_s[below_boundary];
    boundaries.forward_scattered_shortwave_megajoules_m2_h[current_boundary] = 0.0;
    boundaries.forward_scattered_par_umol_m2_s[current_boundary] = 0.0;
    return true;
}

fn validateBoundaryDimensions(boundaries: Boundaries) !void {
    const boundary_count = boundaries.diffuse_transmittance.len;
    if (boundary_count < 2 or
        boundaries.forward_scattered_shortwave_megajoules_m2_h.len != boundary_count or
        boundaries.forward_scattered_par_umol_m2_s.len != boundary_count or
        boundaries.backscattered_shortwave_megajoules_m2_h.len != boundary_count or
        boundaries.backscattered_par_umol_m2_s.len != boundary_count)
        return error.CanopyUpwardBoundaryDimensionMismatch;
}

fn validateLayer(
    layer: usize,
    heights: []const f64,
    snow_depth_m: f64,
    water_ice_depth_m: f64,
    tolerance_m: f64,
    working: WorkingRadiation,
    boundaries: Boundaries,
) !void {
    try validateBoundaryDimensions(boundaries);
    if (heights.len + 1 != boundaries.diffuse_transmittance.len or
        layer >= heights.len)
        return error.CanopyUpwardLayerOutOfRange;
    inline for (.{
        snow_depth_m,
        water_ice_depth_m,
        tolerance_m,
        working.shortwave_megajoules_m2_h,
        working.par_umol_m2_s,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopyUpwardTraversalInput;
    for (heights) |value| if (!std.math.isFinite(value))
        return error.InvalidCanopyUpwardTraversalInput;
    inline for (.{
        boundaries.diffuse_transmittance,
        boundaries.forward_scattered_shortwave_megajoules_m2_h,
        boundaries.forward_scattered_par_umol_m2_s,
        boundaries.backscattered_shortwave_megajoules_m2_h,
        boundaries.backscattered_par_umol_m2_s,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopyUpwardTraversalInput;
}

test "initialization then exposed layer advances upward scattered radiation" {
    var transmission = [_]f64{ 9, 0.8, 0.7 };
    var forward_sw = [_]f64{ 9, 5, 6 };
    var forward_par = [_]f64{ 8, 7, 6 };
    const back_sw = [_]f64{ 2, 3, 4 };
    const back_par = [_]f64{ 20, 30, 40 };
    var working: WorkingRadiation = .{
        .shortwave_megajoules_m2_h = 99,
        .par_umol_m2_s = 98,
    };
    const boundaries: Boundaries = .{
        .diffuse_transmittance = &transmission,
        .forward_scattered_shortwave_megajoules_m2_h = &forward_sw,
        .forward_scattered_par_umol_m2_s = &forward_par,
        .backscattered_shortwave_megajoules_m2_h = &back_sw,
        .backscattered_par_umol_m2_s = &back_par,
    };
    try initialize(&working, boundaries);
    try std.testing.expectEqual(@as(f64, 0), working.shortwave_megajoules_m2_h);
    try std.testing.expectEqual(@as(f64, 1), transmission[0]);
    try std.testing.expect(try admitLayer(
        0,
        &.{ 1, 2 },
        0,
        0,
        0,
        &working,
        boundaries,
    ));
    try std.testing.expectEqual(@as(f64, 2), working.shortwave_megajoules_m2_h);
    try std.testing.expectEqual(@as(f64, 20), working.par_umol_m2_s);
    try std.testing.expectEqual(@as(f64, 0), forward_sw[1]);
    try std.testing.expectEqual(@as(f64, 0), forward_par[1]);
}

test "submerged layer leaves traversal state unchanged" {
    var transmission = [_]f64{ 1, 0.8 };
    var forward_sw = [_]f64{ 1, 2 };
    var forward_par = [_]f64{ 3, 4 };
    const back_sw = [_]f64{ 5, 6 };
    const back_par = [_]f64{ 7, 8 };
    var working: WorkingRadiation = .{
        .shortwave_megajoules_m2_h = 9,
        .par_umol_m2_s = 10,
    };
    try std.testing.expect(!(try admitLayer(0, &.{0}, 1, 0, 0.1, &working, .{
        .diffuse_transmittance = &transmission,
        .forward_scattered_shortwave_megajoules_m2_h = &forward_sw,
        .forward_scattered_par_umol_m2_s = &forward_par,
        .backscattered_shortwave_megajoules_m2_h = &back_sw,
        .backscattered_par_umol_m2_s = &back_par,
    })));
    try std.testing.expectEqual(@as(f64, 9), working.shortwave_megajoules_m2_h);
    try std.testing.expectEqual(@as(f64, 2), forward_sw[1]);
}
