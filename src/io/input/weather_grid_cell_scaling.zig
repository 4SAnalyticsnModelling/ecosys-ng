const std = @import("std");

pub const Intensive = struct {
    rainfall_m_per_h: f64,
    snowfall_water_equivalent_m_per_h: f64,
    surface_irrigation_m_per_h: f64,
    subsurface_irrigation_m_per_h: f64,
    sky_longwave_megajoules_per_m2_h: f64,
};

pub const Extensive = struct {
    rainfall_m3_per_h: f64,
    snowfall_water_equivalent_m3_per_h: f64,
    surface_irrigation_m3_per_h: f64,
    subsurface_irrigation_m3_per_h: f64,
    rain_plus_surface_irrigation_m3_per_h: f64,
    rain_plus_snow_m3_per_h: f64,
    sky_longwave_megajoules_per_h: f64,
};

/// Exact WTHR conversion from per-area forcing to one runtime grid cell.
pub fn scale(cell_area_m2: f64, forcing: Intensive) !Extensive {
    if (!std.math.isFinite(cell_area_m2) or cell_area_m2 <= 0)
        return error.InvalidWeatherCellArea;
    inline for (@typeInfo(Intensive).@"struct".fields) |field| {
        const value = @field(forcing, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteWeatherCellForcing;
        if (value < 0) return error.InvalidWeatherCellForcing;
    }
    const result: Extensive = .{
        .rainfall_m3_per_h = forcing.rainfall_m_per_h * cell_area_m2,
        .snowfall_water_equivalent_m3_per_h = forcing.snowfall_water_equivalent_m_per_h * cell_area_m2,
        .surface_irrigation_m3_per_h = forcing.surface_irrigation_m_per_h * cell_area_m2,
        .subsurface_irrigation_m3_per_h = forcing.subsurface_irrigation_m_per_h * cell_area_m2,
        .rain_plus_surface_irrigation_m3_per_h = (forcing.rainfall_m_per_h +
            forcing.surface_irrigation_m_per_h) * cell_area_m2,
        .rain_plus_snow_m3_per_h = (forcing.rainfall_m_per_h +
            forcing.snowfall_water_equivalent_m_per_h) * cell_area_m2,
        .sky_longwave_megajoules_per_h = forcing.sky_longwave_megajoules_per_m2_h * cell_area_m2,
    };
    inline for (@typeInfo(Extensive).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.WeatherCellScalingOverflow;
    return result;
}

/// Runtime-cell batch form. All inputs and candidates validate before any
/// destination element is replaced.
pub fn scaleCells(
    cell_area_m2: []const f64,
    forcing_by_cell: []const Intensive,
    extensive_by_cell: []Extensive,
) !void {
    const count = cell_area_m2.len;
    if (count == 0 or forcing_by_cell.len != count or
        extensive_by_cell.len != count)
        return error.WeatherCellScalingDimensionMismatch;
    // Two-pass evaluation preserves atomic output without allocating a
    // second whole-domain candidate. The pure scalar calculation is repeated
    // only at this low-cost forcing boundary.
    for (cell_area_m2, forcing_by_cell) |area, forcing|
        _ = try scale(area, forcing);
    for (cell_area_m2, forcing_by_cell, extensive_by_cell) |area, forcing, *out|
        out.* = scale(area, forcing) catch unreachable;
}

test "WTHR cell scaling preserves all extensive carrier identities" {
    const result = try scale(25, .{
        .rainfall_m_per_h = 0.001,
        .snowfall_water_equivalent_m_per_h = 0.002,
        .surface_irrigation_m_per_h = 0.003,
        .subsurface_irrigation_m_per_h = 0.004,
        .sky_longwave_megajoules_per_m2_h = 0.5,
    });
    try std.testing.expectEqual(@as(f64, 0.025), result.rainfall_m3_per_h);
    try std.testing.expectEqual(
        @as(f64, 0.05),
        result.snowfall_water_equivalent_m3_per_h,
    );
    try std.testing.expectEqual(
        @as(f64, 0.075),
        result.surface_irrigation_m3_per_h,
    );
    try std.testing.expectEqual(
        @as(f64, 0.1),
        result.subsurface_irrigation_m3_per_h,
    );
    try std.testing.expectEqual(
        @as(f64, 0.1),
        result.rain_plus_surface_irrigation_m3_per_h,
    );
    try std.testing.expectEqual(
        @as(f64, 0.075),
        result.rain_plus_snow_m3_per_h,
    );
    try std.testing.expectEqual(@as(f64, 12.5), result.sky_longwave_megajoules_per_h);
}

test "runtime-cell batch accepts heterogeneous areas" {
    var output = [_]Extensive{undefined} ** 2;
    try scaleCells(
        &.{ 10, 100 },
        &.{
            .{
                .rainfall_m_per_h = 0.001,
                .snowfall_water_equivalent_m_per_h = 0,
                .surface_irrigation_m_per_h = 0,
                .subsurface_irrigation_m_per_h = 0,
                .sky_longwave_megajoules_per_m2_h = 1,
            },
            .{
                .rainfall_m_per_h = 0.001,
                .snowfall_water_equivalent_m_per_h = 0,
                .surface_irrigation_m_per_h = 0,
                .subsurface_irrigation_m_per_h = 0,
                .sky_longwave_megajoules_per_m2_h = 1,
            },
        },
        &output,
    );
    try std.testing.expectEqual(@as(f64, 0.01), output[0].rainfall_m3_per_h);
    try std.testing.expectEqual(@as(f64, 0.1), output[1].rainfall_m3_per_h);
    try std.testing.expectEqual(@as(f64, 10), output[0].sky_longwave_megajoules_per_h);
    try std.testing.expectEqual(@as(f64, 100), output[1].sky_longwave_megajoules_per_h);
}

test "invalid late cell leaves batch output unchanged" {
    var output = [_]Extensive{
        std.mem.zeroes(Extensive),
        std.mem.zeroes(Extensive),
    };
    output[0].rainfall_m3_per_h = 99;
    const before = output;
    try std.testing.expectError(
        error.NonFiniteWeatherCellForcing,
        scaleCells(
            &.{ 10, 20 },
            &.{
                .{
                    .rainfall_m_per_h = 0.001,
                    .snowfall_water_equivalent_m_per_h = 0,
                    .surface_irrigation_m_per_h = 0,
                    .subsurface_irrigation_m_per_h = 0,
                    .sky_longwave_megajoules_per_m2_h = 1,
                },
                .{
                    .rainfall_m_per_h = 0,
                    .snowfall_water_equivalent_m_per_h = 0,
                    .surface_irrigation_m_per_h = 0,
                    .subsurface_irrigation_m_per_h = 0,
                    .sky_longwave_megajoules_per_m2_h = std.math.nan(f64),
                },
            },
            &output,
        ),
    );
    try std.testing.expectEqualDeep(before, output);
}

test "extensive overflow fails scalar scaling" {
    try std.testing.expectError(
        error.WeatherCellScalingOverflow,
        scale(std.math.floatMax(f64), .{
            .rainfall_m_per_h = 2,
            .snowfall_water_equivalent_m_per_h = 0,
            .surface_irrigation_m_per_h = 0,
            .subsurface_irrigation_m_per_h = 0,
            .sky_longwave_megajoules_per_m2_h = 0,
        }),
    );
}
