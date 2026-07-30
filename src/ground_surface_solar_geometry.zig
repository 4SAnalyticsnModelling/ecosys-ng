const std = @import("std");

pub const Inputs = struct {
    solar_zenith_sine: f64,
    negligible_solar_sine: f64,
    solar_noon_hour: f64,
    source_hour: f64,
    ground_azimuth_rad: f64,
    ground_slope_cosine: f64,
    ground_slope_sine: f64,
    total_canopy_area_m2: f64,
};

pub const DaylightGeometry = struct {
    solar_azimuth_rad: f64,
    ground_incidence_fraction: f64,
    canopy_solar_elevation_rad: ?f64,
};

/// HOUR1 lines 1042--1049. A null result represents the enclosing source
/// night branch, where these daylight-only locals are not evaluated.
pub fn compute(inputs: Inputs) !?DaylightGeometry {
    try validate(inputs);
    if (inputs.solar_zenith_sine <= inputs.negligible_solar_sine)
        return null;
    const solar_azimuth_rad =
        0.2618 * (inputs.solar_noon_hour - inputs.source_hour) + 4.7124;
    const solar_zenith_cosine =
        @sqrt(1.0 - inputs.solar_zenith_sine * inputs.solar_zenith_sine);
    const relative_azimuth_cosine =
        @cos(inputs.ground_azimuth_rad - solar_azimuth_rad);
    const ground_incidence_fraction = @max(
        0,
        @min(
            1,
            inputs.ground_slope_cosine * inputs.solar_zenith_sine +
                inputs.ground_slope_sine * solar_zenith_cosine *
                    relative_azimuth_cosine,
        ),
    );
    return .{
        .solar_azimuth_rad = solar_azimuth_rad,
        .ground_incidence_fraction = ground_incidence_fraction,
        .canopy_solar_elevation_rad = if (inputs.total_canopy_area_m2 > 0)
            std.math.asin(inputs.solar_zenith_sine)
        else
            null,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteGroundSolarGeometryInput;
    if (inputs.solar_zenith_sine < 0 or inputs.solar_zenith_sine > 1 or
        inputs.negligible_solar_sine < 0 or
        inputs.ground_slope_cosine < -1 or inputs.ground_slope_cosine > 1 or
        inputs.ground_slope_sine < -1 or inputs.ground_slope_sine > 1 or
        inputs.total_canopy_area_m2 < 0)
        return error.InvalidGroundSolarGeometryInput;
}

test "flat ground incidence equals solar zenith sine" {
    const result = (try compute(.{
        .solar_zenith_sine = 0.8,
        .negligible_solar_sine = 1.0e-12,
        .solar_noon_hour = 12,
        .source_hour = 10,
        .ground_azimuth_rad = 0,
        .ground_slope_cosine = 1,
        .ground_slope_sine = 0,
        .total_canopy_area_m2 = 2,
    })).?;
    try std.testing.expectApproxEqAbs(@as(f64, 5.236), result.solar_azimuth_rad, 1e-15);
    try std.testing.expectEqual(@as(f64, 0.8), result.ground_incidence_fraction);
    try std.testing.expectApproxEqAbs(
        std.math.asin(@as(f64, 0.8)),
        result.canopy_solar_elevation_rad.?,
        1e-15,
    );
}

test "daylight without canopy leaves canopy angle unset" {
    const result = (try compute(.{
        .solar_zenith_sine = 0.5,
        .negligible_solar_sine = 0,
        .solar_noon_hour = 12,
        .source_hour = 12,
        .ground_azimuth_rad = 0,
        .ground_slope_cosine = 1,
        .ground_slope_sine = 0,
        .total_canopy_area_m2 = 0,
    })).?;
    try std.testing.expect(result.canopy_solar_elevation_rad == null);
}

test "night skips all daylight geometry" {
    try std.testing.expect((try compute(.{
        .solar_zenith_sine = 0,
        .negligible_solar_sine = 1.0e-12,
        .solar_noon_hour = 12,
        .source_hour = 1,
        .ground_azimuth_rad = 0,
        .ground_slope_cosine = 1,
        .ground_slope_sine = 0,
        .total_canopy_area_m2 = 1,
    })) == null);
}
