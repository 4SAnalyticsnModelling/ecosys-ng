const std = @import("std");

pub const Daily = struct {
    maximum_air_temperature_c: f64,
    minimum_air_temperature_c: f64,
    first_vapor_pressure_kpa: f64,
    second_vapor_pressure_kpa: f64,
    shortwave_radiation_mj_per_m2_day: f64,
};

pub const Inputs = struct {
    previous: Daily,
    current: Daily,
    next: Daily,
    current_day_is_scene_begin: bool,
    current_day_is_year_end: bool,
    current_daylength_h: f64,
    radiation_input_type: i8,
    radiation_peak_shape_factor: f64 = 0.658,
};

pub const Parameters = struct {
    maximum_hourly_radiation_mj_per_m2_h: f64,
    previous_to_current_temperature_average_c: f64,
    current_temperature_average_c: f64,
    current_to_next_temperature_average_c: f64,
    previous_temperature_amplitude_c: f64,
    current_temperature_amplitude_c: f64,
    next_temperature_amplitude_c: f64,
    previous_to_current_vapor_average_kpa: f64,
    current_vapor_average_kpa: f64,
    current_to_next_vapor_average_kpa: f64,
    previous_vapor_amplitude_kpa: f64,
    current_vapor_amplitude_kpa: f64,
    next_vapor_amplitude_kpa: f64,
};

/// Exact DAY three-record interpolation preparation. Only previous/current/
/// next records are required, preserving bounded out-of-core weather memory.
pub fn derive(inputs: Inputs) !Parameters {
    inline for (@typeInfo(Daily).@"struct".fields) |field| {
        inline for (.{ inputs.previous, inputs.current, inputs.next }) |day| {
            const value = @field(day, field.name);
            if (!std.math.isFinite(value))
                return error.NonFiniteDailyWeatherInterpolationInput;
            if (comptime std.mem.indexOf(u8, field.name, "radiation") != null)
                if (value < 0) return error.InvalidDailyWeatherInterpolationInput;
        }
    }
    if (!std.math.isFinite(inputs.current_daylength_h) or
        inputs.current_daylength_h < 0 or inputs.current_daylength_h > 24 or
        !std.math.isFinite(inputs.radiation_peak_shape_factor) or
        inputs.radiation_peak_shape_factor <= 0)
        return error.InvalidDailyWeatherInterpolationInput;

    const previous = if (inputs.current_day_is_scene_begin)
        inputs.current
    else
        inputs.previous;
    const next = if (inputs.current_day_is_year_end)
        inputs.current
    else
        inputs.next;
    const radiation_peak =
        if (inputs.radiation_input_type >= -1)
            if (inputs.current_daylength_h > 0)
                inputs.current.shortwave_radiation_mj_per_m2_day /
                    (inputs.current_daylength_h *
                        inputs.radiation_peak_shape_factor)
            else
                0
        else
            inputs.current.shortwave_radiation_mj_per_m2_day;
    const temperature_average_1 =
        0.5 * (previous.maximum_air_temperature_c +
            inputs.current.minimum_air_temperature_c);
    const temperature_average_2 =
        0.5 * (inputs.current.maximum_air_temperature_c +
            inputs.current.minimum_air_temperature_c);
    const temperature_average_3 =
        0.5 * (inputs.current.maximum_air_temperature_c +
            next.minimum_air_temperature_c);
    const vapor_average_1 =
        0.5 * (previous.first_vapor_pressure_kpa +
            inputs.current.second_vapor_pressure_kpa);
    const vapor_average_2 =
        0.5 * (inputs.current.first_vapor_pressure_kpa +
            inputs.current.second_vapor_pressure_kpa);
    const vapor_average_3 =
        0.5 * (inputs.current.first_vapor_pressure_kpa +
            next.second_vapor_pressure_kpa);
    const result: Parameters = .{
        .maximum_hourly_radiation_mj_per_m2_h = radiation_peak,
        .previous_to_current_temperature_average_c = temperature_average_1,
        .current_temperature_average_c = temperature_average_2,
        .current_to_next_temperature_average_c = temperature_average_3,
        .previous_temperature_amplitude_c = temperature_average_1 -
            inputs.current.minimum_air_temperature_c,
        .current_temperature_amplitude_c = temperature_average_2 -
            inputs.current.minimum_air_temperature_c,
        .next_temperature_amplitude_c = temperature_average_3 - next.minimum_air_temperature_c,
        .previous_to_current_vapor_average_kpa = vapor_average_1,
        .current_vapor_average_kpa = vapor_average_2,
        .current_to_next_vapor_average_kpa = vapor_average_3,
        .previous_vapor_amplitude_kpa = vapor_average_1 - inputs.current.second_vapor_pressure_kpa,
        .current_vapor_amplitude_kpa = vapor_average_2 - inputs.current.second_vapor_pressure_kpa,
        .next_vapor_amplitude_kpa = vapor_average_3 - next.second_vapor_pressure_kpa,
    };
    inline for (@typeInfo(Parameters).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.DailyWeatherInterpolationOverflow;
    return result;
}

fn exampleDay(maximum: f64, minimum: f64, first_vapor: f64, second_vapor: f64, radiation: f64) Daily {
    return .{
        .maximum_air_temperature_c = maximum,
        .minimum_air_temperature_c = minimum,
        .first_vapor_pressure_kpa = first_vapor,
        .second_vapor_pressure_kpa = second_vapor,
        .shortwave_radiation_mj_per_m2_day = radiation,
    };
}

test "DAY derives exact three-record temperature vapor and radiation terms" {
    const result = try derive(.{
        .previous = exampleDay(10, 0, 1, 0.5, 5),
        .current = exampleDay(20, 4, 2, 1, 12),
        .next = exampleDay(30, 8, 3, 1.5, 15),
        .current_day_is_scene_begin = false,
        .current_day_is_year_end = false,
        .current_daylength_h = 12,
        .radiation_input_type = 1,
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, 12) / (12 * 0.658),
        result.maximum_hourly_radiation_mj_per_m2_h,
        1e-15,
    );
    try std.testing.expectEqual(
        @as(f64, 7),
        result.previous_to_current_temperature_average_c,
    );
    try std.testing.expectEqual(@as(f64, 12), result.current_temperature_average_c);
    try std.testing.expectEqual(
        @as(f64, 14),
        result.current_to_next_temperature_average_c,
    );
    try std.testing.expectEqual(@as(f64, 3), result.previous_temperature_amplitude_c);
    try std.testing.expectEqual(@as(f64, 8), result.current_temperature_amplitude_c);
    try std.testing.expectEqual(@as(f64, 6), result.next_temperature_amplitude_c);
    try std.testing.expectEqual(
        @as(f64, 1),
        result.previous_to_current_vapor_average_kpa,
    );
}

test "scene begin and year end reuse current record exactly" {
    const result = try derive(.{
        .previous = exampleDay(-99, -99, 99, 99, 0),
        .current = exampleDay(20, 4, 2, 1, 12),
        .next = exampleDay(99, 99, 99, 99, 0),
        .current_day_is_scene_begin = true,
        .current_day_is_year_end = true,
        .current_daylength_h = 12,
        .radiation_input_type = 1,
    });
    try std.testing.expectEqual(
        result.current_temperature_average_c,
        result.previous_to_current_temperature_average_c,
    );
    try std.testing.expectEqual(
        result.current_temperature_average_c,
        result.current_to_next_temperature_average_c,
    );
}

test "zero daylength and direct hourly radiation type preserve source branches" {
    var inputs: Inputs = .{
        .previous = exampleDay(10, 0, 1, 0.5, 5),
        .current = exampleDay(20, 4, 2, 1, 12),
        .next = exampleDay(30, 8, 3, 1.5, 15),
        .current_day_is_scene_begin = false,
        .current_day_is_year_end = false,
        .current_daylength_h = 0,
        .radiation_input_type = 1,
    };
    var result = try derive(inputs);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.maximum_hourly_radiation_mj_per_m2_h,
    );
    inputs.radiation_input_type = -2;
    result = try derive(inputs);
    try std.testing.expectEqual(
        @as(f64, 12),
        result.maximum_hourly_radiation_mj_per_m2_h,
    );
}

test "invalid next record fails before returning interpolation parameters" {
    try std.testing.expectError(
        error.NonFiniteDailyWeatherInterpolationInput,
        derive(.{
            .previous = exampleDay(10, 0, 1, 0.5, 5),
            .current = exampleDay(20, 4, 2, 1, 12),
            .next = exampleDay(30, 8, 3, std.math.nan(f64), 15),
            .current_day_is_scene_begin = false,
            .current_day_is_year_end = false,
            .current_daylength_h = 12,
            .radiation_input_type = 1,
        }),
    );
}
