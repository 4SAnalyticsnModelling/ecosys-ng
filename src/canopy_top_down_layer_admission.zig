const std = @import("std");

pub const WorkingDiffuseRadiation = struct {
    shortwave_mj_per_m2_h: f64,
    par_umol_per_m2_s: f64,
};

pub const LayerBoundaries = struct {
    diffuse_transmittance: []const f64,
    forward_scattered_shortwave_mj_per_m2_h: []f64,
    forward_scattered_par_umol_per_m2_s: []f64,
    backscattered_shortwave_mj_per_m2_h: []f64,
    backscattered_par_umol_per_m2_s: []f64,
};

pub const InterceptionResets = struct {
    diffuse_interception_fraction: f64 = 0,
    direct_species_interception_fraction: f64 = 0,
    diffuse_species_interception_fraction: f64 = 0,
};

/// HOUR1 lines 1191--1202 for one caller-selected layer in the source
/// descending traversal. Boundary arrays have `layer_count + 1` entries.
pub fn admitLayer(
    layer: usize,
    layer_bottom_height_m: []const f64,
    snow_depth_m: f64,
    surface_water_ice_depth_m: f64,
    depth_tolerance_m: f64,
    working: *WorkingDiffuseRadiation,
    boundaries: LayerBoundaries,
) !?InterceptionResets {
    try validate(
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
        return null;

    working.shortwave_mj_per_m2_h =
        working.shortwave_mj_per_m2_h *
        boundaries.diffuse_transmittance[layer + 1] +
        boundaries.forward_scattered_shortwave_mj_per_m2_h[layer + 1];
    working.par_umol_per_m2_s =
        working.par_umol_per_m2_s *
        boundaries.diffuse_transmittance[layer + 1] +
        boundaries.forward_scattered_par_umol_per_m2_s[layer + 1];
    boundaries.forward_scattered_shortwave_mj_per_m2_h[layer] = 0.0;
    boundaries.forward_scattered_par_umol_per_m2_s[layer] = 0.0;
    boundaries.backscattered_shortwave_mj_per_m2_h[layer] = 0.0;
    boundaries.backscattered_par_umol_per_m2_s[layer] = 0.0;
    return .{};
}

fn validate(
    layer: usize,
    heights: []const f64,
    snow_depth_m: f64,
    water_ice_depth_m: f64,
    tolerance_m: f64,
    working: WorkingDiffuseRadiation,
    boundaries: LayerBoundaries,
) !void {
    if (heights.len == 0 or layer >= heights.len)
        return error.CanopyLayerOutOfRange;
    const boundary_count = try std.math.add(usize, heights.len, 1);
    if (boundaries.diffuse_transmittance.len != boundary_count or
        boundaries.forward_scattered_shortwave_mj_per_m2_h.len != boundary_count or
        boundaries.forward_scattered_par_umol_per_m2_s.len != boundary_count or
        boundaries.backscattered_shortwave_mj_per_m2_h.len != boundary_count or
        boundaries.backscattered_par_umol_per_m2_s.len != boundary_count)
        return error.CanopyLayerBoundaryDimensionMismatch;
    inline for (.{
        snow_depth_m,
        water_ice_depth_m,
        tolerance_m,
        working.shortwave_mj_per_m2_h,
        working.par_umol_per_m2_s,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopyLayerAdmissionInput;
    for (heights) |value| if (!std.math.isFinite(value))
        return error.InvalidCanopyLayerAdmissionInput;
    inline for (.{
        boundaries.diffuse_transmittance,
        boundaries.forward_scattered_shortwave_mj_per_m2_h,
        boundaries.forward_scattered_par_umol_per_m2_s,
        boundaries.backscattered_shortwave_mj_per_m2_h,
        boundaries.backscattered_par_umol_per_m2_s,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopyLayerAdmissionInput;
}

test "exposed layer advances diffuse radiation then resets accumulators" {
    var forward_sw = [_]f64{ 9, 1, 2 };
    var forward_par = [_]f64{ 8, 10, 20 };
    var back_sw = [_]f64{ 7, 6, 5 };
    var back_par = [_]f64{ 4, 3, 2 };
    var working: WorkingDiffuseRadiation = .{
        .shortwave_mj_per_m2_h = 4,
        .par_umol_per_m2_s = 100,
    };
    const resets = (try admitLayer(
        1,
        &.{ 0, 2 },
        1,
        0.5,
        0.1,
        &working,
        .{
            .diffuse_transmittance = &.{ 0.1, 0.2, 0.5 },
            .forward_scattered_shortwave_mj_per_m2_h = &forward_sw,
            .forward_scattered_par_umol_per_m2_s = &forward_par,
            .backscattered_shortwave_mj_per_m2_h = &back_sw,
            .backscattered_par_umol_per_m2_s = &back_par,
        },
    )).?;
    try std.testing.expectEqual(@as(f64, 4), working.shortwave_mj_per_m2_h);
    try std.testing.expectEqual(@as(f64, 70), working.par_umol_per_m2_s);
    try std.testing.expectEqual(@as(f64, 0), forward_sw[1]);
    try std.testing.expectEqual(@as(f64, 0), forward_par[1]);
    try std.testing.expectEqual(@as(f64, 0), back_sw[1]);
    try std.testing.expectEqual(@as(f64, 0), back_par[1]);
    try std.testing.expectEqual(@as(f64, 0), resets.diffuse_interception_fraction);
}

test "submerged layer leaves working and boundary state unchanged" {
    var forward_sw = [_]f64{ 9, 1 };
    var forward_par = [_]f64{ 8, 2 };
    var back_sw = [_]f64{ 7, 3 };
    var back_par = [_]f64{ 6, 4 };
    var working: WorkingDiffuseRadiation = .{
        .shortwave_mj_per_m2_h = 5,
        .par_umol_per_m2_s = 50,
    };
    const result = try admitLayer(0, &.{0}, 1, 0, 0.1, &working, .{
        .diffuse_transmittance = &.{ 0.2, 0.5 },
        .forward_scattered_shortwave_mj_per_m2_h = &forward_sw,
        .forward_scattered_par_umol_per_m2_s = &forward_par,
        .backscattered_shortwave_mj_per_m2_h = &back_sw,
        .backscattered_par_umol_per_m2_s = &back_par,
    });
    try std.testing.expect(result == null);
    try std.testing.expectEqual(@as(f64, 5), working.shortwave_mj_per_m2_h);
    try std.testing.expectEqual(@as(f64, 9), forward_sw[0]);
}
