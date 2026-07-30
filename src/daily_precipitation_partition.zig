const std = @import("std");

pub const Inputs = struct {
    source_hour: u8,
    air_temperature_c: f64,
    daily_precipitation_water_equivalent_m: f64,
    daily_wind_travel_m_per_h: f64,
    snowfall_temperature_threshold_c: f64,
    first_precipitation_hour: u8 = 13,
    last_precipitation_hour: u8 = 16,
    minimum_snowfall_water_equivalent_m: f64 = 0.1e-3,
    minimum_wind_travel_m_per_h: f64 = 3600,
};

pub const Result = struct {
    wind_travel_m_per_h: f64,
    rainfall_m: f64,
    snowfall_water_equivalent_m: f64,
};

/// Exact daily-forcing wind and precipitation block from wthr.f:140-159.
pub fn partition(inputs: Inputs) !Result {
    if (inputs.source_hour < 1 or inputs.source_hour > 24 or
        inputs.first_precipitation_hour < 1 or
        inputs.last_precipitation_hour > 24 or
        inputs.first_precipitation_hour > inputs.last_precipitation_hour)
        return error.InvalidDailyPrecipitationHour;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == u8) continue;
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteDailyPrecipitationInput;
    }
    if (inputs.daily_precipitation_water_equivalent_m < 0 or
        inputs.daily_wind_travel_m_per_h < 0 or
        inputs.minimum_snowfall_water_equivalent_m < 0 or
        inputs.minimum_wind_travel_m_per_h < 0)
        return error.InvalidDailyPrecipitationInput;

    var rainfall_m: f64 = 0;
    var snowfall_m: f64 = 0;
    if (inputs.source_hour >= inputs.first_precipitation_hour and
        inputs.source_hour <= inputs.last_precipitation_hour)
    {
        const pulse_hour_count =
            @as(u16, inputs.last_precipitation_hour) -
            @as(u16, inputs.first_precipitation_hour) + 1;
        const hourly_precipitation =
            inputs.daily_precipitation_water_equivalent_m /
            @as(f64, @floatFromInt(pulse_hour_count));
        if (inputs.air_temperature_c >
            inputs.snowfall_temperature_threshold_c)
        {
            rainfall_m = hourly_precipitation;
        } else {
            snowfall_m = hourly_precipitation;
            if (snowfall_m <
                inputs.minimum_snowfall_water_equivalent_m)
                snowfall_m = 0;
        }
    }
    return .{
        .wind_travel_m_per_h = @max(
            inputs.minimum_wind_travel_m_per_h,
            inputs.daily_wind_travel_m_per_h,
        ),
        .rainfall_m = rainfall_m,
        .snowfall_water_equivalent_m = snowfall_m,
    };
}

fn exampleInputs() Inputs {
    return .{
        .source_hour = 13,
        .air_temperature_c = 5,
        .daily_precipitation_water_equivalent_m = 0.02,
        .daily_wind_travel_m_per_h = 100,
        .snowfall_temperature_threshold_c = 0,
    };
}

test "source pulse hours thirteen through sixteen conserve daily rain" {
    var inputs = exampleInputs();
    var total: f64 = 0;
    for (1..25) |hour| {
        inputs.source_hour = @intCast(hour);
        total += (try partition(inputs)).rainfall_m;
    }
    try std.testing.expectEqual(
        inputs.daily_precipitation_water_equivalent_m,
        total,
    );
}

test "hours outside inclusive pulse window receive no precipitation" {
    var inputs = exampleInputs();
    inputs.source_hour = 12;
    var result = try partition(inputs);
    try std.testing.expectEqual(@as(f64, 0), result.rainfall_m);
    inputs.source_hour = 17;
    result = try partition(inputs);
    try std.testing.expectEqual(@as(f64, 0), result.rainfall_m);
}

test "temperature equality follows snowfall branch and wind uses floor" {
    var inputs = exampleInputs();
    inputs.air_temperature_c = inputs.snowfall_temperature_threshold_c;
    const result = try partition(inputs);
    try std.testing.expectEqual(@as(f64, 0), result.rainfall_m);
    try std.testing.expectEqual(
        @as(f64, 0.005),
        result.snowfall_water_equivalent_m,
    );
    try std.testing.expectEqual(
        @as(f64, 3600),
        result.wind_travel_m_per_h,
    );
}

test "daily snowfall pulse below strict cutoff is discarded" {
    var inputs = exampleInputs();
    inputs.air_temperature_c = -1;
    inputs.daily_precipitation_water_equivalent_m = 0.000396;
    var result = try partition(inputs);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.snowfall_water_equivalent_m,
    );
    inputs.daily_precipitation_water_equivalent_m = 0.0004;
    result = try partition(inputs);
    try std.testing.expectEqual(
        @as(f64, 0.0001),
        result.snowfall_water_equivalent_m,
    );
}

test "nonfinite late forcing fails before returning partition" {
    var inputs = exampleInputs();
    inputs.minimum_snowfall_water_equivalent_m = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteDailyPrecipitationInput,
        partition(inputs),
    );
}
