const std = @import("std");

pub const Totals = struct {
    shortwave_radiation_mj_per_m2_day: f64,
    maximum_air_temperature_c: f64,
    minimum_air_temperature_c: f64,
    maximum_vapor_pressure_kpa: f64,
    minimum_vapor_pressure_kpa: f64,
    wind_travel_m_per_day: f64,
    total_water_input_mm_per_day: f64,
};

pub const Hour = struct {
    current_solar_angle_sine: f64,
    seasonal_diffuse_radiation_sine: f64,
    direct_shortwave_mj_per_m2_h: f64,
    diffuse_shortwave_mj_per_m2_h: f64,
    air_temperature_c: f64,
    ambient_vapor_pressure_kpa: f64,
    wind_travel_m_per_h: f64,
    rainfall_m_per_h: f64,
    snowfall_water_equivalent_m_per_h: f64,
    surface_irrigation_m_per_h: f64,
    subsurface_irrigation_m_per_h: f64,
};

/// Exact WTHR daily forcing totals from wthr.f:490-513.
pub fn accumulate(totals: *Totals, hour: Hour) !void {
    inline for (@typeInfo(Totals).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals.*, field.name)))
            return error.NonFiniteDailyAtmosphericTotal;
    inline for (@typeInfo(Hour).@"struct".fields) |field|
        if (!std.math.isFinite(@field(hour, field.name)))
            return error.NonFiniteDailyAtmosphericForcing;
    if (hour.current_solar_angle_sine < 0 or
        hour.current_solar_angle_sine > 1 or
        hour.seasonal_diffuse_radiation_sine < 0 or
        hour.seasonal_diffuse_radiation_sine > 1 or
        hour.direct_shortwave_mj_per_m2_h < 0 or
        hour.diffuse_shortwave_mj_per_m2_h < 0 or
        hour.ambient_vapor_pressure_kpa < 0 or
        hour.wind_travel_m_per_h < 0 or
        hour.rainfall_m_per_h < 0 or
        hour.snowfall_water_equivalent_m_per_h < 0 or
        hour.surface_irrigation_m_per_h < 0 or
        hour.subsurface_irrigation_m_per_h < 0)
        return error.InvalidDailyAtmosphericForcing;

    var next = totals.*;
    if (hour.current_solar_angle_sine > 0) {
        next.shortwave_radiation_mj_per_m2_day +=
            hour.direct_shortwave_mj_per_m2_h *
            hour.current_solar_angle_sine +
            hour.diffuse_shortwave_mj_per_m2_h *
                hour.seasonal_diffuse_radiation_sine;
    }
    next.maximum_air_temperature_c =
        @max(next.maximum_air_temperature_c, hour.air_temperature_c);
    next.minimum_air_temperature_c =
        @min(next.minimum_air_temperature_c, hour.air_temperature_c);
    next.maximum_vapor_pressure_kpa =
        @max(next.maximum_vapor_pressure_kpa, hour.ambient_vapor_pressure_kpa);
    next.minimum_vapor_pressure_kpa =
        @min(next.minimum_vapor_pressure_kpa, hour.ambient_vapor_pressure_kpa);
    next.wind_travel_m_per_day += hour.wind_travel_m_per_h;
    next.total_water_input_mm_per_day +=
        (hour.rainfall_m_per_h +
            hour.snowfall_water_equivalent_m_per_h +
            hour.surface_irrigation_m_per_h +
            hour.subsurface_irrigation_m_per_h) * 1000;
    inline for (@typeInfo(Totals).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.DailyAtmosphericTotalOverflow;
    totals.* = next;
}

pub fn reset() Totals {
    return .{
        .shortwave_radiation_mj_per_m2_day = 0,
        .maximum_air_temperature_c = -100,
        .minimum_air_temperature_c = 100,
        .maximum_vapor_pressure_kpa = 0,
        .minimum_vapor_pressure_kpa = 100,
        .wind_travel_m_per_day = 0,
        .total_water_input_mm_per_day = 0,
    };
}

fn exampleHour() Hour {
    return .{
        .current_solar_angle_sine = 0.5,
        .seasonal_diffuse_radiation_sine = 0.8,
        .direct_shortwave_mj_per_m2_h = 2,
        .diffuse_shortwave_mj_per_m2_h = 1,
        .air_temperature_c = 15,
        .ambient_vapor_pressure_kpa = 1.2,
        .wind_travel_m_per_h = 3600,
        .rainfall_m_per_h = 0.001,
        .snowfall_water_equivalent_m_per_h = 0.002,
        .surface_irrigation_m_per_h = 0.003,
        .subsurface_irrigation_m_per_h = 0.004,
    };
}

test "hour updates exact WTHR totals and extrema" {
    var totals = reset();
    try accumulate(&totals, exampleHour());
    try std.testing.expectEqual(
        @as(f64, 1.8),
        totals.shortwave_radiation_mj_per_m2_day,
    );
    try std.testing.expectEqual(@as(f64, 15), totals.maximum_air_temperature_c);
    try std.testing.expectEqual(@as(f64, 15), totals.minimum_air_temperature_c);
    try std.testing.expectEqual(@as(f64, 1.2), totals.maximum_vapor_pressure_kpa);
    try std.testing.expectEqual(@as(f64, 3600), totals.wind_travel_m_per_day);
    try std.testing.expectEqual(@as(f64, 10), totals.total_water_input_mm_per_day);
}

test "zero solar angle suppresses both direct and diffuse radiation total" {
    var totals = reset();
    var hour = exampleHour();
    hour.current_solar_angle_sine = 0;
    try accumulate(&totals, hour);
    try std.testing.expectEqual(
        @as(f64, 0),
        totals.shortwave_radiation_mj_per_m2_day,
    );
}

test "multiple hours retain correct extrema and accumulate travel" {
    var totals = reset();
    var hour = exampleHour();
    try accumulate(&totals, hour);
    hour.air_temperature_c = -5;
    hour.ambient_vapor_pressure_kpa = 0.4;
    try accumulate(&totals, hour);
    try std.testing.expectEqual(@as(f64, 15), totals.maximum_air_temperature_c);
    try std.testing.expectEqual(@as(f64, -5), totals.minimum_air_temperature_c);
    try std.testing.expectEqual(@as(f64, 1.2), totals.maximum_vapor_pressure_kpa);
    try std.testing.expectEqual(@as(f64, 0.4), totals.minimum_vapor_pressure_kpa);
    try std.testing.expectEqual(@as(f64, 7200), totals.wind_travel_m_per_day);
}

test "invalid late forcing leaves every daily total unchanged" {
    var totals = reset();
    const before = totals;
    var hour = exampleHour();
    hour.subsurface_irrigation_m_per_h = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteDailyAtmosphericForcing,
        accumulate(&totals, hour),
    );
    try std.testing.expectEqualDeep(before, totals);
}
