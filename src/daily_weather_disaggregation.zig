const std = @import("std");
const execution_calendar_date = @import("execution_calendar_date.zig");
const DailyForcing = @import("weather.zig").DailyForcing;
const HourlyForcing = @import("weather.zig").HourlyForcing;

const twilight_sine = 0.06976;
const PrecipitationPhase = @import("weather.zig").PrecipitationPhase;

pub const Result = struct {
    forcing: HourlyForcing,
    rainfall_m: f64,
    snowfall_m: f64,
    daylength_h: f64,
};

/// Exact daily-weather curve equations from DAY.F and WTHR.F. The caller
/// supplies adjacent days; at a file boundary it should repeat the current
/// day, matching the reference boundary behavior.
pub fn disaggregateHour(previous: DailyForcing, current: DailyForcing, next: DailyForcing, day_of_year: u16, execution_year: u16, hour: u8, latitude_degrees_north: f64, solar_noon_hour: f64, altitude_m: f64, phytotron: bool) !Result {
    return disaggregateHourWithPhase(previous, current, next, day_of_year, execution_year, hour, latitude_degrees_north, solar_noon_hour, altitude_m, phytotron, .{});
}

pub fn disaggregateHourWithPhase(previous: DailyForcing, current: DailyForcing, next: DailyForcing, day_of_year: u16, execution_year: u16, hour: u8, latitude_degrees_north: f64, solar_noon_hour: f64, altitude_m: f64, phytotron: bool, phase: PrecipitationPhase) !Result {
    try phase.validate();
    _ = execution_calendar_date.fromDayOfYear(
        day_of_year,
        execution_year,
    ) catch return error.InvalidDailyWeatherTime;
    if (hour < 1 or hour > 24) return error.InvalidDailyWeatherTime;
    if (!std.math.isFinite(latitude_degrees_north) or latitude_degrees_north < -90 or latitude_degrees_north > 90 or !std.math.isFinite(solar_noon_hour) or !std.math.isFinite(altitude_m)) return error.InvalidDailyWeatherGeometry;
    const daylength_h = if (phytotron) 24.0 else calculateDaylength(day_of_year, latitude_degrees_north);
    const maximum_radiation = if (phytotron)
        current.shortwave_radiation_mj_per_m2_per_day
    else if (daylength_h > 0)
        current.shortwave_radiation_mj_per_m2_per_day / (daylength_h * 0.658)
    else
        0.0;
    const hour_f: f64 = @floatFromInt(hour);
    const shortwave = if (phytotron)
        maximum_radiation / 24.0
    else if (daylength_h > 0)
        @max(0.0, maximum_radiation * @sin((hour_f - (solar_noon_hour - daylength_h / 2.0)) * 3.1416 / daylength_h))
    else
        0.0;

    const temperature_curve = curveParameters(
        previous.maximum_air_temperature_c,
        current.maximum_air_temperature_c,
        current.minimum_air_temperature_c,
        next.minimum_air_temperature_c,
    );
    const vapor_curve = curveParameters(
        previous.mean_vapor_pressure_kpa,
        current.mean_vapor_pressure_kpa,
        current.saturation_vapor_pressure_at_minimum_kpa,
        next.saturation_vapor_pressure_at_minimum_kpa,
    );
    const temperature_c = curveValue(temperature_curve, hour_f, solar_noon_hour, daylength_h);
    const raw_vapor = curveValue(vapor_curve, hour_f, solar_noon_hour, daylength_h);
    const saturation = 0.61 * @exp(5360.0 * (0.003661 - 1.0 / (273.15 + temperature_c))) * @exp(-altitude_m / 7272.0);
    const precipitation = if (hour >= 13 and hour <= 16) current.precipitation_m_per_day / 4.0 else 0.0;
    var rainfall: f64 = 0;
    var snowfall: f64 = 0;
    if (temperature_c > phase.snowfall_temperature_threshold_c) rainfall = precipitation else {
        snowfall = precipitation;
        if (snowfall < phase.minimum_snowfall_water_equivalent_m) snowfall = 0;
    }
    const forcing: HourlyForcing = .{
        .air_temperature_c = temperature_c,
        .vapor_pressure_kpa = @min(saturation, raw_vapor),
        .precipitation_m = rainfall + snowfall,
        .rainfall_m = rainfall,
        .snowfall_water_equivalent_m = snowfall,
        .shortwave_radiation_mj_per_m2 = shortwave,
        .wind_speed_m_per_h = @max(3600.0, current.wind_speed_m_per_h),
        .longwave_radiation_mj_per_m2 = null,
    };
    try validate(forcing);
    return .{ .forcing = forcing, .rainfall_m = rainfall, .snowfall_m = snowfall, .daylength_h = daylength_h };
}

const Curve = struct { average_before: f64, average_current: f64, average_after: f64, amplitude_before: f64, amplitude_current: f64, amplitude_after: f64 };

fn curveParameters(previous_max_or_mean: f64, current_max_or_mean: f64, current_minimum: f64, next_minimum: f64) Curve {
    const average_before = 0.5 * (previous_max_or_mean + current_minimum);
    const average_current = 0.5 * (current_max_or_mean + current_minimum);
    const average_after = 0.5 * (current_max_or_mean + next_minimum);
    return .{
        .average_before = average_before,
        .average_current = average_current,
        .average_after = average_after,
        .amplitude_before = average_before - current_minimum,
        .amplitude_current = average_current - current_minimum,
        .amplitude_after = average_after - next_minimum,
    };
}

fn curveValue(curve: Curve, hour: f64, solar_noon: f64, daylength: f64) f64 {
    const sunrise = solar_noon - daylength / 2.0;
    const denominator = solar_noon + 9.0 - daylength / 2.0;
    if (hour < sunrise) return curve.average_before + curve.amplitude_before * @sin((hour + solar_noon - 3.0) * 3.1416 / denominator + 1.5708);
    if (hour > solar_noon + 3.0) return curve.average_after + curve.amplitude_after * @sin((hour - solar_noon - 3.0) * 3.1416 / denominator + 1.5708);
    return curve.average_current + curve.amplitude_current * @sin((hour - sunrise) * 3.1416 / (3.0 + daylength / 2.0) - 1.5708);
}

pub fn calculateDaylength(day_of_year: u16, latitude_degrees: f64) f64 {
    const effective_day = if (day_of_year == 366) 365.5 else @as(f64, @floatFromInt(day_of_year));
    const declination_degrees = @sin((effective_day + 100.0) * 0.9863 * 1.7453e-2) * -23.47;
    const latitude_radians = latitude_degrees * 1.7453e-2;
    const declination_radians = declination_degrees * 1.7453e-2;
    const azimuth_component = @sin(latitude_radians) * @sin(declination_radians);
    const declination_component = @cos(latitude_radians) * @cos(declination_radians);
    const ratio = azimuth_component / declination_component;
    if (ratio >= 1.0 - twilight_sine) return 24.0;
    if (ratio <= -1.0 + twilight_sine) return 0.0;
    return 12.0 * (1.0 + 2.0 / 3.1416 * std.math.asin(twilight_sine + ratio));
}

fn validate(forcing: HourlyForcing) !void {
    inline for (@typeInfo(HourlyForcing).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(forcing, field.name))) return error.NonFiniteDisaggregatedWeather;
    if (forcing.air_temperature_c < -273.15 or forcing.vapor_pressure_kpa < 0 or forcing.precipitation_m < 0 or forcing.shortwave_radiation_mj_per_m2 < 0 or forcing.wind_speed_m_per_h < 0) return error.InvalidDisaggregatedWeather;
}

test "daily precipitation is assigned to reference hours and conserved when rain" {
    const day: DailyForcing = .{
        .maximum_air_temperature_c = 20,
        .minimum_air_temperature_c = 10,
        .mean_vapor_pressure_kpa = 1,
        .saturation_vapor_pressure_at_minimum_kpa = 1.2,
        .precipitation_m_per_day = 0.02,
        .shortwave_radiation_mj_per_m2_per_day = 12,
        .wind_speed_m_per_h = 100,
    };
    var total: f64 = 0;
    for (1..25) |hour| total += (try disaggregateHour(day, day, day, 180, 2000, @intCast(hour), 53.69, 12, 645, false)).rainfall_m;
    try std.testing.expectApproxEqAbs(@as(f64, 0.02), total, 1.0e-12);
}

test "polar day and night remain finite" {
    const day: DailyForcing = .{ .maximum_air_temperature_c = -5, .minimum_air_temperature_c = -15, .mean_vapor_pressure_kpa = 0.2, .saturation_vapor_pressure_at_minimum_kpa = 0.1, .precipitation_m_per_day = 0, .shortwave_radiation_mj_per_m2_per_day = 5, .wind_speed_m_per_h = 0 };
    const summer = try disaggregateHour(day, day, day, 180, 2000, 12, 81.8, 12, 290, false);
    const winter = try disaggregateHour(day, day, day, 1, 2000, 12, 81.8, 12, 290, false);
    try std.testing.expectEqual(@as(f64, 24), summer.daylength_h);
    try std.testing.expectEqual(@as(f64, 0), winter.daylength_h);
}

test "daily disaggregation preserves DAY modulo-four chronology" {
    const day: DailyForcing = .{ .maximum_air_temperature_c = 5, .minimum_air_temperature_c = -5, .mean_vapor_pressure_kpa = 0.2, .saturation_vapor_pressure_at_minimum_kpa = 0.1, .precipitation_m_per_day = 0, .shortwave_radiation_mj_per_m2_per_day = 5, .wind_speed_m_per_h = 3600 };
    _ = try disaggregateHour(day, day, day, 366, 1900, 12, 53.5, 12, 645, false);
    try std.testing.expectError(
        error.InvalidDailyWeatherTime,
        disaggregateHour(day, day, day, 366, 1901, 12, 53.5, 12, 645, false),
    );
}
