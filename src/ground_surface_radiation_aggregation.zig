const std = @import("std");

pub const Inputs = struct {
    direct_shortwave_megajoules_per_m2_h: f64,
    direct_par_umol_per_m2_s: f64,
    bottom_direct_transmittance: f64,
    working_diffuse_shortwave_megajoules_per_m2_h: f64,
    working_diffuse_par_umol_per_m2_s: f64,
    bottom_diffuse_transmittance: f64,
    bottom_forward_scattered_shortwave_megajoules_per_m2_h: f64,
    bottom_forward_scattered_par_umol_per_m2_s: f64,
    ground_direct_incidence_fraction: f64,
    ground_sky_incidence_fraction: []const f64,
    cell_area_m2: f64,
};

pub const Result = struct {
    transmitted_direct_shortwave_megajoules_per_m2_h: f64,
    transmitted_diffuse_shortwave_megajoules_per_m2_h: f64,
    transmitted_direct_par_umol_per_m2_s: f64,
    transmitted_diffuse_par_umol_per_m2_s: f64,
    ground_shortwave_megajoules_per_m2_h: f64,
    ground_par_umol_per_m2_s: f64,
    ground_shortwave_megajoules_h: f64,
};

/// HOUR1 lines 1608--1618. Runtime sky-zone accumulation preserves the
/// direct initialization followed by source-order diffuse additions.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    const transmitted_direct_shortwave_megajoules_per_m2_h =
        inputs.direct_shortwave_megajoules_per_m2_h * inputs.bottom_direct_transmittance;
    const transmitted_diffuse_shortwave_megajoules_per_m2_h =
        inputs.working_diffuse_shortwave_megajoules_per_m2_h *
        inputs.bottom_diffuse_transmittance +
        inputs.bottom_forward_scattered_shortwave_megajoules_per_m2_h;
    const transmitted_direct_par_umol_per_m2_s =
        inputs.direct_par_umol_per_m2_s * inputs.bottom_direct_transmittance;
    const transmitted_diffuse_par_umol_per_m2_s =
        inputs.working_diffuse_par_umol_per_m2_s *
        inputs.bottom_diffuse_transmittance +
        inputs.bottom_forward_scattered_par_umol_per_m2_s;
    var ground_shortwave_megajoules_per_m2_h =
        @abs(inputs.ground_direct_incidence_fraction) *
        transmitted_direct_shortwave_megajoules_per_m2_h;
    var ground_par_umol_per_m2_s =
        @abs(inputs.ground_direct_incidence_fraction) *
        transmitted_direct_par_umol_per_m2_s;
    for (inputs.ground_sky_incidence_fraction) |incidence| {
        ground_shortwave_megajoules_per_m2_h +=
            @abs(incidence) * transmitted_diffuse_shortwave_megajoules_per_m2_h;
        ground_par_umol_per_m2_s +=
            @abs(incidence) * transmitted_diffuse_par_umol_per_m2_s;
    }
    const ground_shortwave_megajoules_h =
        ground_shortwave_megajoules_per_m2_h * inputs.cell_area_m2;
    return .{
        .transmitted_direct_shortwave_megajoules_per_m2_h = transmitted_direct_shortwave_megajoules_per_m2_h,
        .transmitted_diffuse_shortwave_megajoules_per_m2_h = transmitted_diffuse_shortwave_megajoules_per_m2_h,
        .transmitted_direct_par_umol_per_m2_s = transmitted_direct_par_umol_per_m2_s,
        .transmitted_diffuse_par_umol_per_m2_s = transmitted_diffuse_par_umol_per_m2_s,
        .ground_shortwave_megajoules_per_m2_h = ground_shortwave_megajoules_per_m2_h,
        .ground_par_umol_per_m2_s = ground_par_umol_per_m2_s,
        .ground_shortwave_megajoules_h = ground_shortwave_megajoules_h,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (.{
        inputs.direct_shortwave_megajoules_per_m2_h,
        inputs.direct_par_umol_per_m2_s,
        inputs.bottom_direct_transmittance,
        inputs.working_diffuse_shortwave_megajoules_per_m2_h,
        inputs.working_diffuse_par_umol_per_m2_s,
        inputs.bottom_diffuse_transmittance,
        inputs.bottom_forward_scattered_shortwave_megajoules_per_m2_h,
        inputs.bottom_forward_scattered_par_umol_per_m2_s,
        inputs.cell_area_m2,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidGroundRadiationAggregationInput;
    if (inputs.cell_area_m2 == 0 or
        inputs.bottom_direct_transmittance > 1 or
        inputs.bottom_diffuse_transmittance > 1 or
        !std.math.isFinite(inputs.ground_direct_incidence_fraction) or
        @abs(inputs.ground_direct_incidence_fraction) > 1)
        return error.InvalidGroundRadiationAggregationInput;
    if (inputs.ground_sky_incidence_fraction.len == 0)
        return error.ZeroGroundSkyZoneExtent;
    for (inputs.ground_sky_incidence_fraction) |value|
        if (!std.math.isFinite(value) or @abs(value) > 1)
            return error.InvalidGroundRadiationAggregationInput;
}

test "direct and runtime sky-zone radiation aggregate in source order" {
    const result = try compute(.{
        .direct_shortwave_megajoules_per_m2_h = 2,
        .direct_par_umol_per_m2_s = 100,
        .bottom_direct_transmittance = 0.5,
        .working_diffuse_shortwave_megajoules_per_m2_h = 1,
        .working_diffuse_par_umol_per_m2_s = 40,
        .bottom_diffuse_transmittance = 0.8,
        .bottom_forward_scattered_shortwave_megajoules_per_m2_h = 0.2,
        .bottom_forward_scattered_par_umol_per_m2_s = 8,
        .ground_direct_incidence_fraction = -0.5,
        .ground_sky_incidence_fraction = &.{ -0.25, 0.75 },
        .cell_area_m2 = 10,
    });
    try std.testing.expectEqual(@as(f64, 1), result.transmitted_direct_shortwave_megajoules_per_m2_h);
    try std.testing.expectEqual(@as(f64, 1), result.transmitted_diffuse_shortwave_megajoules_per_m2_h);
    try std.testing.expectEqual(@as(f64, 1.5), result.ground_shortwave_megajoules_per_m2_h);
    try std.testing.expectEqual(@as(f64, 65), result.ground_par_umol_per_m2_s);
    try std.testing.expectEqual(@as(f64, 15), result.ground_shortwave_megajoules_h);
}

test "empty sky-zone extent fails explicitly" {
    try std.testing.expectError(error.ZeroGroundSkyZoneExtent, compute(.{
        .direct_shortwave_megajoules_per_m2_h = 0,
        .direct_par_umol_per_m2_s = 0,
        .bottom_direct_transmittance = 1,
        .working_diffuse_shortwave_megajoules_per_m2_h = 0,
        .working_diffuse_par_umol_per_m2_s = 0,
        .bottom_diffuse_transmittance = 1,
        .bottom_forward_scattered_shortwave_megajoules_per_m2_h = 0,
        .bottom_forward_scattered_par_umol_per_m2_s = 0,
        .ground_direct_incidence_fraction = 0,
        .ground_sky_incidence_fraction = &.{},
        .cell_area_m2 = 1,
    }));
}
