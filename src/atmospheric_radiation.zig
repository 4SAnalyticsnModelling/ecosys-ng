const std = @import("std");
const HourlyForcing = @import("weather.zig").HourlyForcing;
const Timestamp = @import("weather.zig").Timestamp;
const execution_calendar_date = @import("execution_calendar_date.zig");

pub const Result = struct {
    forcing: HourlyForcing,
    solar_angle_sine: f64,
    next_hour_solar_angle_sine: f64,
    solar_azimuth_radians: f64,
    extraterrestrial_shortwave_mj_per_m2: f64,
    cloudiness_fraction: f64,
    sky_emissivity: f64,
};

pub fn prepare(forcing: HourlyForcing, timestamp: Timestamp, latitude_degrees_north: f64, solar_noon_hour: f64, phytotron: bool) !Result {
    const day = try validateRadiationDay(timestamp);
    if (!std.math.isFinite(latitude_degrees_north) or latitude_degrees_north < -90 or latitude_degrees_north > 90 or !std.math.isFinite(solar_noon_hour)) return error.InvalidRadiationGeometry;
    const air_temperature_k = forcing.air_temperature_c + 273.15;
    if (!std.math.isFinite(air_temperature_k) or air_temperature_k <= 0 or !std.math.isFinite(forcing.vapor_pressure_kpa) or forcing.vapor_pressure_kpa < 0) return error.InvalidRadiationForcing;
    var adjusted = forcing;
    var solar_sine: f64 = 0;
    var next_sine: f64 = 0;
    var extraterrestrial: f64 = 0;
    var solar_azimuth_radians: f64 = 0;
    var cloudiness: f64 = 0;
    var emissivity: f64 = 0;
    if (!phytotron) {
        const geometry = solarGeometry(day, latitude_degrees_north);
        const hour: f64 = @as(f64, @floatFromInt(timestamp.hour)) + @as(f64, @floatFromInt(timestamp.minute)) / 60.0;
        solar_azimuth_radians = (std.math.pi / 12.0) * (solar_noon_hour - hour) + 3.0 * std.math.pi / 2.0;
        solar_sine = @max(0.0, geometry.azimuth_component + geometry.declination_component * @cos(0.2618 * (solar_noon_hour - (hour - 0.5))));
        next_sine = @max(0.0, geometry.azimuth_component + geometry.declination_component * @cos(0.2618 * (solar_noon_hour - (hour + 0.5))));
        if (adjusted.shortwave_radiation_mj_per_m2 <= 0) solar_sine = 0;
        extraterrestrial = 4.896 * @max(0.0, solar_sine);
        adjusted.shortwave_radiation_mj_per_m2 = @min(extraterrestrial, adjusted.shortwave_radiation_mj_per_m2);
        cloudiness = if (extraterrestrial > 0) std.math.clamp(2.33 - 3.33 * adjusted.shortwave_radiation_mj_per_m2 / extraterrestrial, 0.2, 1.0) else 0.2;
        emissivity = 0.625 * @max(1.0, std.math.pow(f64, 1000.0 * forcing.vapor_pressure_kpa / air_temperature_k, 0.131));
        emissivity *= 1.0 + 0.242 * std.math.pow(f64, cloudiness, 0.583);
    } else {
        solar_sine = if (adjusted.shortwave_radiation_mj_per_m2 > 0) 1 else 0;
        next_sine = 1;
        cloudiness = 0;
        emissivity = 0.97;
    }
    adjusted.longwave_radiation_mj_per_m2 = if (forcing.longwave_radiation_mj_per_m2) |observed|
        if (observed > 0) observed else emissivity * 2.04e-10 * std.math.pow(f64, air_temperature_k, 4)
    else
        emissivity * 2.04e-10 * std.math.pow(f64, air_temperature_k, 4);
    try validate(adjusted);
    return .{
        .forcing = adjusted,
        .solar_angle_sine = solar_sine,
        .next_hour_solar_angle_sine = next_sine,
        .solar_azimuth_radians = solar_azimuth_radians,
        .extraterrestrial_shortwave_mj_per_m2 = extraterrestrial,
        .cloudiness_fraction = cloudiness,
        .sky_emissivity = emissivity,
    };
}

fn validateRadiationDay(timestamp: Timestamp) !u16 {
    const day = timestamp.day_of_year orelse
        return error.RadiationRequiresDayOfYear;
    const maximum_day: u16 = if (timestamp.year) |year| maximum: {
        if (year == 0) return error.InvalidRadiationGeometry;
        break :maximum if (execution_calendar_date.isLeapYear(year)) 366 else 365;
    } else 366;
    if (day == 0 or day > maximum_day) return error.InvalidRadiationGeometry;
    return day;
}

fn solarGeometry(day_of_year: u16, latitude_degrees: f64) struct { azimuth_component: f64, declination_component: f64 } {
    const effective_day = if (day_of_year == 366) 365.5 else @as(f64, @floatFromInt(day_of_year));
    const declination_degrees = @sin((effective_day + 100.0) * 0.9863 * 1.7453e-2) * -23.47;
    const latitude = latitude_degrees * 1.7453e-2;
    const declination = declination_degrees * 1.7453e-2;
    return .{ .azimuth_component = @sin(latitude) * @sin(declination), .declination_component = @cos(latitude) * @cos(declination) };
}

fn validate(forcing: HourlyForcing) !void {
    inline for (@typeInfo(HourlyForcing).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(forcing, field.name))) return error.NonFiniteAtmosphericRadiation;
    if (forcing.shortwave_radiation_mj_per_m2 < 0 or forcing.longwave_radiation_mj_per_m2.? < 0) return error.InvalidAtmosphericRadiation;
}

test "outdoor radiation is capped and longwave is finite" {
    const forcing: HourlyForcing = .{ .air_temperature_c = 20, .vapor_pressure_kpa = 1, .precipitation_m = 0, .shortwave_radiation_mj_per_m2 = 100, .wind_speed_m_per_h = 3600, .longwave_radiation_mj_per_m2 = null };
    const result = try prepare(forcing, .{ .year = 2024, .day_of_year = 180, .month = null, .day_of_month = null, .hour = 12, .minute = 0 }, 53.69, 12, false);
    try std.testing.expect(result.forcing.shortwave_radiation_mj_per_m2 <= result.extraterrestrial_shortwave_mj_per_m2);
    try std.testing.expect(result.forcing.longwave_radiation_mj_per_m2.? > 0);
}

test "observed longwave is retained" {
    const forcing: HourlyForcing = .{ .air_temperature_c = 5, .vapor_pressure_kpa = 0.5, .precipitation_m = 0, .shortwave_radiation_mj_per_m2 = 0, .wind_speed_m_per_h = 3600, .longwave_radiation_mj_per_m2 = 0.8 };
    const result = try prepare(forcing, .{ .year = null, .day_of_year = 1, .month = null, .day_of_month = null, .hour = 1, .minute = 0 }, 81.8, 12, false);
    try std.testing.expectEqual(@as(f64, 0.8), result.forcing.longwave_radiation_mj_per_m2.?);
}

test "radiation dates preserve DAY modulo-four chronology" {
    try std.testing.expectEqual(
        @as(u16, 366),
        try validateRadiationDay(.{ .year = 1900, .day_of_year = 366, .month = null, .day_of_month = null, .hour = 1, .minute = 0 }),
    );
    try std.testing.expectError(
        error.InvalidRadiationGeometry,
        validateRadiationDay(.{ .year = 1901, .day_of_year = 366, .month = null, .day_of_month = null, .hour = 1, .minute = 0 }),
    );
    try std.testing.expectError(
        error.InvalidRadiationGeometry,
        validateRadiationDay(.{ .year = 0, .day_of_year = 1, .month = null, .day_of_month = null, .hour = 1, .minute = 0 }),
    );
    try std.testing.expectEqual(
        @as(u16, 366),
        try validateRadiationDay(.{ .year = null, .day_of_year = 366, .month = null, .day_of_month = null, .hour = 1, .minute = 0 }),
    );
}
