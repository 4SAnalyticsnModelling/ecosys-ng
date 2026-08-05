const std = @import("std");

pub const Inputs = struct {
    ground_direct_incidence_fraction: f64,
    direct_shortwave_megajoules_per_m2_h: f64,
    ground_sky_incidence_fraction: []const f64,
    diffuse_shortwave_megajoules_per_m2_h: f64,
    cell_area_m2: f64,
};

pub const SpeciesDiagnostics = struct {
    live_shortwave_megajoules_h: []f64,
    standing_dead_shortwave_megajoules_h: []f64,
    live_par_umol_s: []f64,
    standing_dead_par_umol_s: []f64,
};

/// HOUR1 control lines 1739--1751 and executable lines 1752--1762. Applies
/// the no-canopy branch and resets runtime-species canopy diagnostics.
pub fn apply(inputs: Inputs, diagnostics: SpeciesDiagnostics) !f64 {
    const species_count = try validate(inputs, diagnostics);
    var ground_shortwave_megajoules_per_m2_h =
        @abs(inputs.ground_direct_incidence_fraction) *
        inputs.direct_shortwave_megajoules_per_m2_h;
    for (inputs.ground_sky_incidence_fraction) |incidence|
        ground_shortwave_megajoules_per_m2_h +=
            @abs(incidence) * inputs.diffuse_shortwave_megajoules_per_m2_h;
    const ground_shortwave_megajoules_h =
        ground_shortwave_megajoules_per_m2_h * inputs.cell_area_m2;
    for (0..species_count) |species| {
        diagnostics.live_shortwave_megajoules_h[species] = 0.0;
        diagnostics.standing_dead_shortwave_megajoules_h[species] = 0.0;
        diagnostics.live_par_umol_s[species] = 0.0;
        diagnostics.standing_dead_par_umol_s[species] = 0.0;
    }
    return ground_shortwave_megajoules_h;
}

fn validate(inputs: Inputs, diagnostics: SpeciesDiagnostics) !usize {
    inline for (.{
        inputs.direct_shortwave_megajoules_per_m2_h,
        inputs.diffuse_shortwave_megajoules_per_m2_h,
        inputs.cell_area_m2,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidNoCanopyRadiationInput;
    if (inputs.cell_area_m2 == 0 or
        !std.math.isFinite(inputs.ground_direct_incidence_fraction) or
        @abs(inputs.ground_direct_incidence_fraction) > 1)
        return error.InvalidNoCanopyRadiationInput;
    if (inputs.ground_sky_incidence_fraction.len == 0)
        return error.ZeroNoCanopySkyZoneExtent;
    for (inputs.ground_sky_incidence_fraction) |value|
        if (!std.math.isFinite(value) or @abs(value) > 1)
            return error.InvalidNoCanopyRadiationInput;
    const species_count = diagnostics.live_shortwave_megajoules_h.len;
    if (species_count == 0 or
        diagnostics.standing_dead_shortwave_megajoules_h.len != species_count or
        diagnostics.live_par_umol_s.len != species_count or
        diagnostics.standing_dead_par_umol_s.len != species_count)
        return error.NoCanopyRadiationDimensionMismatch;
    return species_count;
}

test "no canopy ground radiation accumulates sky zones then resets species" {
    var live_sw = [_]f64{ 1, 2, 3 };
    var dead_sw = [_]f64{ 4, 5, 6 };
    var live_par = [_]f64{ 7, 8, 9 };
    var dead_par = [_]f64{ 10, 11, 12 };
    const ground_shortwave_megajoules_h = try apply(.{
        .ground_direct_incidence_fraction = -0.5,
        .direct_shortwave_megajoules_per_m2_h = 2,
        .ground_sky_incidence_fraction = &.{ -0.25, 0.75 },
        .diffuse_shortwave_megajoules_per_m2_h = 1,
        .cell_area_m2 = 10,
    }, .{
        .live_shortwave_megajoules_h = &live_sw,
        .standing_dead_shortwave_megajoules_h = &dead_sw,
        .live_par_umol_s = &live_par,
        .standing_dead_par_umol_s = &dead_par,
    });
    try std.testing.expectEqual(@as(f64, 20), ground_shortwave_megajoules_h);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, &live_sw);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, &dead_sw);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, &live_par);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, &dead_par);
}

test "dimension mismatch leaves species diagnostics unchanged" {
    var live_sw = [_]f64{42};
    var empty = [_]f64{};
    try std.testing.expectError(
        error.NoCanopyRadiationDimensionMismatch,
        apply(.{
            .ground_direct_incidence_fraction = 0,
            .direct_shortwave_megajoules_per_m2_h = 0,
            .ground_sky_incidence_fraction = &.{0},
            .diffuse_shortwave_megajoules_per_m2_h = 0,
            .cell_area_m2 = 1,
        }, .{
            .live_shortwave_megajoules_h = &live_sw,
            .standing_dead_shortwave_megajoules_h = &empty,
            .live_par_umol_s = &live_sw,
            .standing_dead_par_umol_s = &live_sw,
        }),
    );
    try std.testing.expectEqual(@as(f64, 42), live_sw[0]);
}
