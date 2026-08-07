const std = @import("std");

pub const ScatteringSide = enum(u8) {
    forward = 1,
    backward = 2,
};

pub const Inputs = struct {
    solar_azimuth_rad: f64,
    solar_elevation_rad: f64,
    solar_zenith_sine: f64,
    solar_zenith_cosine: f64,
    leaf_inclination_cosine: []const f64,
    leaf_inclination_sine: []const f64,
    azimuth_class_count: usize,
};

pub const Outputs = struct {
    /// Azimuth-major, then inclination, matching source M -> N traversal.
    direct_incidence_fraction: []f64,
    horizontal_relative_incidence: []f64,
    scattering_side: []ScatteringSide,
};

const source_pi = 3.1416;
const source_half_pi = 1.5708;

/// `hour1.f` lines 1076--1097. Runtime angular extents retain the source
/// azimuth-major then inclination traversal; four azimuth classes reproduce
/// the original `(M-0.5)*3.1416/4` centers exactly.
pub fn compute(inputs: Inputs, outputs: Outputs) !void {
    try validate(inputs, outputs);
    const azimuth_width_rad =
        source_pi / @as(f64, @floatFromInt(inputs.azimuth_class_count));
    for (0..inputs.azimuth_class_count) |azimuth| {
        const leaf_azimuth_rad = inputs.solar_azimuth_rad +
            (@as(f64, @floatFromInt(azimuth + 1)) - 0.5) *
                azimuth_width_rad;
        const relative_azimuth_cosine =
            @cos(leaf_azimuth_rad - inputs.solar_azimuth_rad);
        for (inputs.leaf_inclination_cosine, inputs.leaf_inclination_sine, 0..) |
            inclination_cosine,
            inclination_sine,
            inclination,
        | {
            const index =
                azimuth * inputs.leaf_inclination_cosine.len + inclination;
            const signed_incidence = inclination_cosine *
                inputs.solar_zenith_sine +
                inclination_sine * inputs.solar_zenith_cosine *
                    relative_azimuth_cosine;
            if (signed_incidence < -1 or signed_incidence > 1)
                return error.FoliarIncidenceOutsideAcosDomain;
            const incidence = @abs(signed_incidence);
            outputs.direct_incidence_fraction[index] = incidence;
            outputs.horizontal_relative_incidence[index] =
                incidence / inputs.solar_zenith_sine;
            const signed_surface_angle = if (inclination_cosine >
                inputs.solar_zenith_sine)
                std.math.acos(signed_incidence)
            else
                -std.math.acos(signed_incidence);
            const scattering_angle = if (signed_surface_angle > -source_half_pi)
                inputs.solar_elevation_rad + 2.0 * signed_surface_angle
            else
                inputs.solar_elevation_rad -
                    2.0 * (source_pi + signed_surface_angle);
            outputs.scattering_side[index] =
                if (scattering_angle > 0 and scattering_angle < source_pi)
                    .forward
                else
                    .backward;
        }
    }
}

fn validate(inputs: Inputs, outputs: Outputs) !void {
    const inclination_count = inputs.leaf_inclination_cosine.len;
    if (inclination_count == 0 or inputs.azimuth_class_count == 0)
        return error.ZeroFoliarAngleExtent;
    const output_count = try std.math.mul(
        usize,
        inclination_count,
        inputs.azimuth_class_count,
    );
    if (inputs.leaf_inclination_sine.len != inclination_count or
        outputs.direct_incidence_fraction.len != output_count or
        outputs.horizontal_relative_incidence.len != output_count or
        outputs.scattering_side.len != output_count)
        return error.FoliarAngleDimensionMismatch;
    inline for (.{
        inputs.solar_azimuth_rad,
        inputs.solar_elevation_rad,
        inputs.solar_zenith_sine,
        inputs.solar_zenith_cosine,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteFoliarAngleInput;
    if (inputs.solar_zenith_sine <= 0 or inputs.solar_zenith_sine > 1 or
        inputs.solar_zenith_cosine < 0 or inputs.solar_zenith_cosine > 1)
        return error.InvalidFoliarAngleInput;
    inline for (.{
        inputs.leaf_inclination_cosine,
        inputs.leaf_inclination_sine,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < -1 or value > 1)
            return error.InvalidFoliarAngleInput;
}

test "four by four foliar angles preserve source traversal and equations" {
    var incidence: [16]f64 = undefined;
    var relative: [16]f64 = undefined;
    var side: [16]ScatteringSide = undefined;
    try compute(.{
        .solar_azimuth_rad = 4.7124,
        .solar_elevation_rad = std.math.asin(@as(f64, 0.8)),
        .solar_zenith_sine = 0.8,
        .solar_zenith_cosine = 0.6,
        .leaf_inclination_cosine = &.{ 1, 0.8, 0.6, 0 },
        .leaf_inclination_sine = &.{ 0, 0.6, 0.8, 1 },
        .azimuth_class_count = 4,
    }, .{
        .direct_incidence_fraction = &incidence,
        .horizontal_relative_incidence = &relative,
        .scattering_side = &side,
    });
    try std.testing.expectEqual(@as(f64, 0.8), incidence[0]);
    try std.testing.expectEqual(@as(f64, 1), relative[0]);
    try std.testing.expect(side[0] == .forward);
}

test "runtime angular dimensions are not fixed to four" {
    var incidence: [6]f64 = undefined;
    var relative: [6]f64 = undefined;
    var side: [6]ScatteringSide = undefined;
    try compute(.{
        .solar_azimuth_rad = 0,
        .solar_elevation_rad = 0.5,
        .solar_zenith_sine = 0.5,
        .solar_zenith_cosine = @sqrt(0.75),
        .leaf_inclination_cosine = &.{ 1, 0 },
        .leaf_inclination_sine = &.{ 0, 1 },
        .azimuth_class_count = 3,
    }, .{
        .direct_incidence_fraction = &incidence,
        .horizontal_relative_incidence = &relative,
        .scattering_side = &side,
    });
}
