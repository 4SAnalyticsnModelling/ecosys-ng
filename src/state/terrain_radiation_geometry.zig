const std = @import("std");

pub const Inputs = struct {
    horizontal_cell_width_m: []const f64,
    vertical_cell_width_m: []const f64,
    aspect_degrees: []const f64,
    slope_degrees: []const f64,
    sky_azimuth_rad: []const f64,
    sky_elevation_sine: []const f64,
    sky_elevation_cosine: []const f64,
    calculation_floor: f64,
    geometry_floor: f64,
    degrees_per_radian: f64,
    minimum_slope_m_per_m: f64,
};

pub const State = struct {
    area_scaled_calculation_floor_m2: []f64,
    area_scaled_geometry_floor_m2: []f64,
    ground_azimuth_rad: []f64,
    absolute_ground_azimuth_sine: []f64,
    absolute_ground_azimuth_cosine: []f64,
    ground_slope_sine_m_per_m: []f64,
    ground_slope_cosine: []f64,
    /// Cell-major `[cell][sky_azimuth]`, dimensionless.
    incident_sky_projection_fraction: []f64,
};

/// Exact supplemental translation of legacy `STARTS` lines 269--274 and
/// 315--321. The intervening directional slope/routing and elevation block is
/// already owned by `terrain_hydrology.zig`.
pub fn derive(state: State, inputs: Inputs) !void {
    const cell_count = inputs.horizontal_cell_width_m.len;
    const sky_count = inputs.sky_elevation_sine.len;
    if (cell_count == 0 or sky_count == 0 or
        inputs.vertical_cell_width_m.len != cell_count or
        inputs.aspect_degrees.len != cell_count or
        inputs.slope_degrees.len != cell_count or
        inputs.sky_azimuth_rad.len != sky_count or
        inputs.sky_elevation_cosine.len != sky_count or
        state.area_scaled_calculation_floor_m2.len != cell_count or
        state.area_scaled_geometry_floor_m2.len != cell_count or
        state.ground_azimuth_rad.len != cell_count or
        state.absolute_ground_azimuth_sine.len != cell_count or
        state.absolute_ground_azimuth_cosine.len != cell_count or
        state.ground_slope_sine_m_per_m.len != cell_count or
        state.ground_slope_cosine.len != cell_count)
    {
        return error.TerrainRadiationDimensionMismatch;
    }
    const projection_count = std.math.mul(
        usize,
        cell_count,
        sky_count,
    ) catch return error.DimensionOverflow;
    if (state.incident_sky_projection_fraction.len != projection_count)
        return error.TerrainRadiationDimensionMismatch;

    inline for (.{
        inputs.calculation_floor,
        inputs.geometry_floor,
        inputs.degrees_per_radian,
        inputs.minimum_slope_m_per_m,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteTerrainRadiationInput;
    }
    if (inputs.calculation_floor <= 0 or inputs.geometry_floor <= 0 or
        inputs.degrees_per_radian <= 0 or
        inputs.minimum_slope_m_per_m < 0 or
        inputs.minimum_slope_m_per_m >= 1)
        return error.InvalidTerrainRadiationInput;
    inline for (.{
        inputs.horizontal_cell_width_m,
        inputs.vertical_cell_width_m,
        inputs.aspect_degrees,
        inputs.slope_degrees,
        inputs.sky_azimuth_rad,
        inputs.sky_elevation_sine,
        inputs.sky_elevation_cosine,
    }) |values| {
        for (values) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteTerrainRadiationInput;
        }
    }
    for (0..cell_count) |cell| {
        if (inputs.horizontal_cell_width_m[cell] <= 0 or
            inputs.vertical_cell_width_m[cell] <= 0 or
            inputs.aspect_degrees[cell] < 0 or
            inputs.aspect_degrees[cell] > 360 or
            inputs.slope_degrees[cell] < 0 or
            inputs.slope_degrees[cell] > 90)
            return error.InvalidTerrainRadiationInput;
    }
    for (0..sky_count) |sky| {
        if (@abs(inputs.sky_elevation_sine[sky]) > 1 or
            @abs(inputs.sky_elevation_cosine[sky]) > 1)
            return error.InvalidTerrainRadiationInput;
    }

    for (0..cell_count) |cell| {
        const area_m2 = inputs.horizontal_cell_width_m[cell] *
            inputs.vertical_cell_width_m[cell];
        state.area_scaled_calculation_floor_m2[cell] =
            inputs.calculation_floor * area_m2;
        state.area_scaled_geometry_floor_m2[cell] =
            inputs.geometry_floor * area_m2;
        state.ground_azimuth_rad[cell] =
            inputs.aspect_degrees[cell] / inputs.degrees_per_radian;
        state.absolute_ground_azimuth_sine[cell] =
            @abs(@sin(state.ground_azimuth_rad[cell]));
        state.absolute_ground_azimuth_cosine[cell] =
            @abs(@cos(state.ground_azimuth_rad[cell]));
        state.ground_slope_sine_m_per_m[cell] = @max(
            inputs.minimum_slope_m_per_m,
            @sin(inputs.slope_degrees[cell] / inputs.degrees_per_radian),
        );
        state.ground_slope_cosine[cell] = @sqrt(
            1.0 - state.ground_slope_sine_m_per_m[cell] *
                state.ground_slope_sine_m_per_m[cell],
        );
        for (0..sky_count) |sky| {
            const relative_azimuth_cosine =
                @cos(state.ground_azimuth_rad[cell] -
                    inputs.sky_azimuth_rad[sky]);
            state.incident_sky_projection_fraction[cell * sky_count + sky] =
                @max(0.0, @min(
                    1.0,
                    state.ground_slope_cosine[cell] *
                        inputs.sky_elevation_sine[sky] +
                        state.ground_slope_sine_m_per_m[cell] *
                            inputs.sky_elevation_cosine[sky] *
                            relative_azimuth_cosine,
                ));
        }
    }
}

test "STARTS terrain radiation geometry preserves source calculations" {
    var scalar_outputs = [_]f64{0.0} ** 7;
    var projections = [_]f64{0.0} ** 4;
    try derive(.{
        .area_scaled_calculation_floor_m2 = scalar_outputs[0..1],
        .area_scaled_geometry_floor_m2 = scalar_outputs[1..2],
        .ground_azimuth_rad = scalar_outputs[2..3],
        .absolute_ground_azimuth_sine = scalar_outputs[3..4],
        .absolute_ground_azimuth_cosine = scalar_outputs[4..5],
        .ground_slope_sine_m_per_m = scalar_outputs[5..6],
        .ground_slope_cosine = scalar_outputs[6..7],
        .incident_sky_projection_fraction = &projections,
    }, .{
        .horizontal_cell_width_m = &.{10.0},
        .vertical_cell_width_m = &.{20.0},
        .aspect_degrees = &.{45.0},
        .slope_degrees = &.{30.0},
        .sky_azimuth_rad = &.{ 3.1416 / 4.0, 3.0 * 3.1416 / 4.0, 5.0 * 3.1416 / 4.0, 7.0 * 3.1416 / 4.0 },
        .sky_elevation_sine = &.{ @sin(3.1416 / 4.0), @sin(3.1416 / 4.0), @sin(3.1416 / 4.0), @sin(3.1416 / 4.0) },
        .sky_elevation_cosine = &.{ @cos(3.1416 / 4.0), @cos(3.1416 / 4.0), @cos(3.1416 / 4.0), @cos(3.1416 / 4.0) },
        .calculation_floor = 1.0e-15,
        .geometry_floor = 1.0e-6,
        .degrees_per_radian = 57.2958,
        .minimum_slope_m_per_m = 1.745e-3,
    });

    try std.testing.expectEqual(@as(f64, 2.0e-13), scalar_outputs[0]);
    try std.testing.expectApproxEqAbs(
        @as(f64, 2.0e-4),
        scalar_outputs[1],
        1.0e-19,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 45.0 / 57.2958),
        scalar_outputs[2],
        1.0e-15,
    );
    for (projections) |projection| {
        try std.testing.expect(projection >= 0 and projection <= 1);
    }
}

test "runtime cells and sky classes determine projection extent" {
    var cell_fields = [_]f64{0.0} ** (7 * 2);
    var projections = [_]f64{0.0} ** 6;
    try derive(.{
        .area_scaled_calculation_floor_m2 = cell_fields[0..2],
        .area_scaled_geometry_floor_m2 = cell_fields[2..4],
        .ground_azimuth_rad = cell_fields[4..6],
        .absolute_ground_azimuth_sine = cell_fields[6..8],
        .absolute_ground_azimuth_cosine = cell_fields[8..10],
        .ground_slope_sine_m_per_m = cell_fields[10..12],
        .ground_slope_cosine = cell_fields[12..14],
        .incident_sky_projection_fraction = &projections,
    }, .{
        .horizontal_cell_width_m = &.{ 10, 20 },
        .vertical_cell_width_m = &.{ 30, 40 },
        .aspect_degrees = &.{ 0, 270 },
        .slope_degrees = &.{ 0, 5 },
        .sky_azimuth_rad = &.{ 0.5, 1.5, 2.5 },
        .sky_elevation_sine = &.{ 0.5, 0.7, 0.9 },
        .sky_elevation_cosine = &.{ 0.8660254, 0.7141428, 0.4358899 },
        .calculation_floor = 1.0e-15,
        .geometry_floor = 1.0e-6,
        .degrees_per_radian = 57.2958,
        .minimum_slope_m_per_m = 1.745e-3,
    });
    for (projections) |projection|
        try std.testing.expect(std.math.isFinite(projection));
}
