const std = @import("std");

pub const Inputs = struct {
    base_horizontal_shortwave_mj_m2_h: f64,
    horizontal_shortwave_change_mj_m2_h: f64,
    day_fraction: f64,
    hour_fraction: f64,
    extraterrestrial_shortwave_mj_m2_h: f64,
    solar_zenith_sine: f64,
    diffuse_sky_projection: f64,
};

pub const Result = struct {
    direct_shortwave_mj_m2_h: f64,
    diffuse_shortwave_mj_m2_h: f64,
    direct_par_umol_m2_s: f64,
    diffuse_par_umol_m2_s: f64,
    horizontal_shortwave_mj_m2_h: f64,
    horizontal_par_umol_m2_s: f64,
};

const maximum_direct_shortwave_mj_m2_h = 4.167;
const direct_visible_fraction = 0.42;
const diffuse_visible_fraction = 0.58;
const par_umol_m2_s_per_mj_m2_h = 1269.4;

/// HOUR1 lines 994--1022. Includes the source RADM interpolation assembly
/// and preserves the SSIN daylight/night branch.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    const horizontal_shortwave_mj_m2_h =
        inputs.base_horizontal_shortwave_mj_m2_h +
        inputs.day_fraction * inputs.hour_fraction *
            inputs.horizontal_shortwave_change_mj_m2_h;
    if (!std.math.isFinite(horizontal_shortwave_mj_m2_h) or
        horizontal_shortwave_mj_m2_h < 0)
        return error.InvalidInterpolatedCanopyRadiation;
    if (inputs.solar_zenith_sine <= 0)
        return zeroResult();

    const diffuse_horizontal_mj_m2_h = @max(
        0,
        @min(
            0.5 * (inputs.extraterrestrial_shortwave_mj_m2_h -
                horizontal_shortwave_mj_m2_h),
            horizontal_shortwave_mj_m2_h,
        ),
    );
    const direct_shortwave_mj_m2_h = @min(
        maximum_direct_shortwave_mj_m2_h,
        (horizontal_shortwave_mj_m2_h - diffuse_horizontal_mj_m2_h) /
            inputs.solar_zenith_sine,
    );
    const diffuse_shortwave_mj_m2_h =
        diffuse_horizontal_mj_m2_h / inputs.diffuse_sky_projection;
    const direct_par_umol_m2_s =
        direct_shortwave_mj_m2_h * direct_visible_fraction *
        par_umol_m2_s_per_mj_m2_h;
    const diffuse_par_umol_m2_s =
        diffuse_shortwave_mj_m2_h * diffuse_visible_fraction *
        par_umol_m2_s_per_mj_m2_h;
    return .{
        .direct_shortwave_mj_m2_h = direct_shortwave_mj_m2_h,
        .diffuse_shortwave_mj_m2_h = diffuse_shortwave_mj_m2_h,
        .direct_par_umol_m2_s = direct_par_umol_m2_s,
        .diffuse_par_umol_m2_s = diffuse_par_umol_m2_s,
        .horizontal_shortwave_mj_m2_h = direct_shortwave_mj_m2_h * inputs.solar_zenith_sine +
            diffuse_shortwave_mj_m2_h * inputs.diffuse_sky_projection,
        .horizontal_par_umol_m2_s = direct_par_umol_m2_s * inputs.solar_zenith_sine +
            diffuse_par_umol_m2_s * inputs.diffuse_sky_projection,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteCanopyRadiationInput;
    if (inputs.base_horizontal_shortwave_mj_m2_h < 0 or
        inputs.extraterrestrial_shortwave_mj_m2_h < 0 or
        inputs.day_fraction < 0 or inputs.day_fraction > 1 or
        inputs.hour_fraction < 0 or inputs.hour_fraction > 1 or
        inputs.solar_zenith_sine < 0 or inputs.solar_zenith_sine > 1 or
        inputs.diffuse_sky_projection <= 0)
        return error.InvalidCanopyRadiationInput;
}

fn zeroResult() Result {
    return .{
        .direct_shortwave_mj_m2_h = 0,
        .diffuse_shortwave_mj_m2_h = 0,
        .direct_par_umol_m2_s = 0,
        .diffuse_par_umol_m2_s = 0,
        .horizontal_shortwave_mj_m2_h = 0,
        .horizontal_par_umol_m2_s = 0,
    };
}

test "hourly assembly precedes direct diffuse partition" {
    const result = try compute(.{
        .base_horizontal_shortwave_mj_m2_h = 1,
        .horizontal_shortwave_change_mj_m2_h = 2,
        .day_fraction = 0.5,
        .hour_fraction = 0.5,
        .extraterrestrial_shortwave_mj_m2_h = 3,
        .solar_zenith_sine = 0.75,
        .diffuse_sky_projection = 2,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), result.horizontal_shortwave_mj_m2_h, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), result.direct_shortwave_mj_m2_h, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.375), result.diffuse_shortwave_mj_m2_h, 1e-15);
}

test "night zeros outputs even when assembled shortwave is nonzero" {
    const result = try compute(.{
        .base_horizontal_shortwave_mj_m2_h = 1,
        .horizontal_shortwave_change_mj_m2_h = 0,
        .day_fraction = 0,
        .hour_fraction = 0,
        .extraterrestrial_shortwave_mj_m2_h = 1,
        .solar_zenith_sine = 0,
        .diffuse_sky_projection = 2,
    });
    inline for (@typeInfo(Result).@"struct".fields) |field|
        try std.testing.expectEqual(@as(f64, 0), @field(result, field.name));
}

test "invalid interpolated shortwave fails explicitly" {
    try std.testing.expectError(error.InvalidInterpolatedCanopyRadiation, compute(.{
        .base_horizontal_shortwave_mj_m2_h = 0,
        .horizontal_shortwave_change_mj_m2_h = -2,
        .day_fraction = 1,
        .hour_fraction = 1,
        .extraterrestrial_shortwave_mj_m2_h = 1,
        .solar_zenith_sine = 0.5,
        .diffuse_sky_projection = 2,
    }));
}
