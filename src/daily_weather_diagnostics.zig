const std = @import("std");

pub const State = struct {
    solar_radiation_mj_per_m2: f64,
    maximum_air_temperature_c: f64,
    minimum_air_temperature_c: f64,
    maximum_vapor_pressure_kpa: f64,
    minimum_vapor_pressure_kpa: f64,
    vapor_concentration_m3_per_m3: f64,
    wind_travel_m: f64,
    total_water_input_mm: f64,
    net_radiation_change_mj_per_m2: f64,
    air_temperature_change_k: f64,
    vapor_concentration_change_m3_per_m3: f64,
};

pub const Hour = struct {
    source_hour: u8,
    first_scene_day: bool,
    solar_angle_sine: f64,
    seasonal_daylight_sine: f64,
    direct_shortwave_mj_per_m2: f64,
    diffuse_shortwave_mj_per_m2: f64,
    air_temperature_c: f64,
    air_temperature_k: f64,
    previous_air_temperature_k: f64,
    vapor_pressure_kpa: f64,
    wind_speed_m_per_h: f64,
    rainfall_m: f64,
    snowfall_water_equivalent_m: f64,
    surface_irrigation_m: f64,
    subsurface_irrigation_m: f64,
    net_radiation_mj_per_m2: f64,
    previous_net_radiation_mj_per_m2: f64,
    previous_vapor_concentration_m3_per_m3: f64,
};

/// Atomic per-cell WTHR daily diagnostic update.
pub fn accumulate(state: *State, hour: Hour) !void {
    if (hour.source_hour < 1 or hour.source_hour > 24)
        return error.InvalidDailyWeatherSourceHour;
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (!std.math.isFinite(@field(state.*, field.name)))
            return error.NonFiniteDailyWeatherState;
    inline for (@typeInfo(Hour).@"struct".fields) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(hour, field.name)))
            return error.NonFiniteDailyWeatherInput;
    }
    if (hour.solar_angle_sine < -1 or hour.solar_angle_sine > 1 or
        hour.seasonal_daylight_sine < -1 or
        hour.seasonal_daylight_sine > 1 or
        hour.direct_shortwave_mj_per_m2 < 0 or
        hour.diffuse_shortwave_mj_per_m2 < 0 or
        hour.air_temperature_k <= 0 or hour.vapor_pressure_kpa < 0 or
        hour.wind_speed_m_per_h < 0 or hour.rainfall_m < 0 or
        hour.snowfall_water_equivalent_m < 0 or
        hour.surface_irrigation_m < 0 or hour.subsurface_irrigation_m < 0)
        return error.InvalidDailyWeatherInput;

    var next = state.*;
    if (hour.solar_angle_sine > 0) {
        next.solar_radiation_mj_per_m2 +=
            hour.direct_shortwave_mj_per_m2 * hour.solar_angle_sine +
            hour.diffuse_shortwave_mj_per_m2 *
                hour.seasonal_daylight_sine;
    }
    next.maximum_air_temperature_c =
        @max(next.maximum_air_temperature_c, hour.air_temperature_c);
    next.minimum_air_temperature_c =
        @min(next.minimum_air_temperature_c, hour.air_temperature_c);
    next.maximum_vapor_pressure_kpa =
        @max(next.maximum_vapor_pressure_kpa, hour.vapor_pressure_kpa);
    next.minimum_vapor_pressure_kpa =
        @min(next.minimum_vapor_pressure_kpa, hour.vapor_pressure_kpa);
    next.vapor_concentration_m3_per_m3 =
        hour.vapor_pressure_kpa * 2.173e-3 / hour.air_temperature_k;
    next.wind_travel_m += hour.wind_speed_m_per_h;
    next.total_water_input_mm += 1000 *
        (hour.rainfall_m + hour.snowfall_water_equivalent_m +
            hour.surface_irrigation_m + hour.subsurface_irrigation_m);
    if (hour.first_scene_day and hour.source_hour == 1) {
        next.net_radiation_change_mj_per_m2 = 0;
        next.air_temperature_change_k = 0;
        next.vapor_concentration_change_m3_per_m3 = 0;
    } else {
        next.net_radiation_change_mj_per_m2 =
            hour.net_radiation_mj_per_m2 -
            hour.previous_net_radiation_mj_per_m2;
        next.air_temperature_change_k =
            hour.air_temperature_k - hour.previous_air_temperature_k;
        next.vapor_concentration_change_m3_per_m3 =
            next.vapor_concentration_m3_per_m3 -
            hour.previous_vapor_concentration_m3_per_m3;
    }
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.DailyWeatherDiagnosticOverflow;
    state.* = next;
}

fn initialState() State {
    return .{
        .solar_radiation_mj_per_m2 = 0,
        .maximum_air_temperature_c = -100,
        .minimum_air_temperature_c = 100,
        .maximum_vapor_pressure_kpa = 0,
        .minimum_vapor_pressure_kpa = 100,
        .vapor_concentration_m3_per_m3 = 0,
        .wind_travel_m = 0,
        .total_water_input_mm = 0,
        .net_radiation_change_mj_per_m2 = 0,
        .air_temperature_change_k = 0,
        .vapor_concentration_change_m3_per_m3 = 0,
    };
}

fn exampleHour() Hour {
    return .{
        .source_hour = 12,
        .first_scene_day = false,
        .solar_angle_sine = 0.5,
        .seasonal_daylight_sine = 0.8,
        .direct_shortwave_mj_per_m2 = 2,
        .diffuse_shortwave_mj_per_m2 = 1,
        .air_temperature_c = 20,
        .air_temperature_k = 293.15,
        .previous_air_temperature_k = 292.15,
        .vapor_pressure_kpa = 1.2,
        .wind_speed_m_per_h = 100,
        .rainfall_m = 0.001,
        .snowfall_water_equivalent_m = 0.002,
        .surface_irrigation_m = 0.003,
        .subsurface_irrigation_m = 0.004,
        .net_radiation_mj_per_m2 = 3,
        .previous_net_radiation_mj_per_m2 = 2.5,
        .previous_vapor_concentration_m3_per_m3 = 1e-6,
    };
}

test "WTHR daily diagnostics include all four water carriers" {
    var state = initialState();
    try accumulate(&state, exampleHour());
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.8),
        state.solar_radiation_mj_per_m2,
        1e-15,
    );
    try std.testing.expectEqual(@as(f64, 20), state.maximum_air_temperature_c);
    try std.testing.expectEqual(@as(f64, 20), state.minimum_air_temperature_c);
    try std.testing.expectEqual(@as(f64, 10), state.total_water_input_mm);
    try std.testing.expectEqual(@as(f64, 100), state.wind_travel_m);
    try std.testing.expectEqual(@as(f64, 0.5), state.net_radiation_change_mj_per_m2);
    try std.testing.expectEqual(@as(f64, 1), state.air_temperature_change_k);
}

test "night hour excludes radiation but retains weather totals" {
    var state = initialState();
    var hour = exampleHour();
    hour.solar_angle_sine = 0;
    try accumulate(&state, hour);
    try std.testing.expectEqual(@as(f64, 0), state.solar_radiation_mj_per_m2);
    try std.testing.expectEqual(@as(f64, 10), state.total_water_input_mm);
}

test "first scene source hour clears all hourly changes" {
    var state = initialState();
    var hour = exampleHour();
    hour.source_hour = 1;
    hour.first_scene_day = true;
    try accumulate(&state, hour);
    try std.testing.expectEqual(@as(f64, 0), state.net_radiation_change_mj_per_m2);
    try std.testing.expectEqual(@as(f64, 0), state.air_temperature_change_k);
    try std.testing.expectEqual(
        @as(f64, 0),
        state.vapor_concentration_change_m3_per_m3,
    );
}

test "invalid late forcing rolls back all daily diagnostics" {
    var state = initialState();
    const before = state;
    var hour = exampleHour();
    hour.subsurface_irrigation_m = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteDailyWeatherInput,
        accumulate(&state, hour),
    );
    try std.testing.expectEqualDeep(before, state);
}
