const std = @import("std");

pub const State = struct {
    air_temperature_c: f64,
    air_temperature_k: f64,
    vapor_concentration_m3_per_m3: f64,
    previous_air_temperature_k: f64,
    previous_vapor_concentration_m3_per_m3: f64,
};

pub const Inputs = struct {
    source_hour: u8,
    solar_noon_source_hour: f64,
    first_scene_day: bool,
    maximum_temperature_change_c: f64,
    minimum_temperature_change_c: f64,
    initial_mean_annual_air_temperature_k: f64,
};

pub const Result = struct {
    daily_average_change_c: f64,
    daily_amplitude_change_c: f64,
    diurnal_amplitude_factor: f64,
    applied: bool,
};

/// Atomic WTHR diurnal climate temperature adjustment, including the exact
/// previous-hour snapshot needed by subsequent daily-change diagnostics.
pub fn apply(state: *State, inputs: Inputs) !Result {
    if (inputs.source_hour < 1 or inputs.source_hour > 24)
        return error.InvalidWeatherTemperatureSourceHour;
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (!std.math.isFinite(@field(state.*, field.name)))
            return error.NonFiniteWeatherTemperatureState;
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteWeatherTemperatureInput;
    if (state.air_temperature_k <= 0 or
        inputs.initial_mean_annual_air_temperature_k <= 0 or
        inputs.solar_noon_source_hour < 1 or
        inputs.solar_noon_source_hour >= 25)
        return error.InvalidWeatherTemperatureInput;
    if (inputs.maximum_temperature_change_c == 0 and
        inputs.minimum_temperature_change_c == 0)
        return .{
            .daily_average_change_c = 0,
            .daily_amplitude_change_c = 0,
            .diurnal_amplitude_factor = 0,
            .applied = false,
        };

    const average_change_c = 0.5 *
        (inputs.maximum_temperature_change_c +
            inputs.minimum_temperature_change_c);
    const amplitude_change_c = 0.5 *
        (inputs.maximum_temperature_change_c -
            inputs.minimum_temperature_change_c);
    const source_hour_f64: f64 = @floatFromInt(inputs.source_hour);
    const diurnal_factor = @sin(
        0.2618 *
            (source_hour_f64 -
                (inputs.solar_noon_source_hour + 3)) +
            1.5708,
    );
    const next_temperature_c = state.air_temperature_c +
        average_change_c + amplitude_change_c * diurnal_factor;
    const next_temperature_k = next_temperature_c + 273.15;
    if (!std.math.isFinite(next_temperature_c) or
        !std.math.isFinite(next_temperature_k) or next_temperature_k <= 0)
        return error.WeatherTemperatureChangeOverflow;

    state.previous_air_temperature_k =
        if (inputs.first_scene_day and inputs.source_hour == 1)
            inputs.initial_mean_annual_air_temperature_k
        else
            state.air_temperature_k;
    state.previous_vapor_concentration_m3_per_m3 =
        if (inputs.first_scene_day and inputs.source_hour == 1)
            0
        else
            state.vapor_concentration_m3_per_m3;
    state.air_temperature_c = next_temperature_c;
    state.air_temperature_k = next_temperature_k;
    return .{
        .daily_average_change_c = average_change_c,
        .daily_amplitude_change_c = amplitude_change_c,
        .diurnal_amplitude_factor = diurnal_factor,
        .applied = true,
    };
}

fn exampleState() State {
    return .{
        .air_temperature_c = 20,
        .air_temperature_k = 293.15,
        .vapor_concentration_m3_per_m3 = 0.001,
        .previous_air_temperature_k = 290,
        .previous_vapor_concentration_m3_per_m3 = 0.0009,
    };
}

test "maximum and minimum changes form exact average and amplitude" {
    var state = exampleState();
    const result = try apply(&state, .{
        .source_hour = 15,
        .solar_noon_source_hour = 12,
        .first_scene_day = false,
        .maximum_temperature_change_c = 4,
        .minimum_temperature_change_c = 2,
        .initial_mean_annual_air_temperature_k = 280,
    });
    try std.testing.expect(result.applied);
    try std.testing.expectEqual(@as(f64, 3), result.daily_average_change_c);
    try std.testing.expectEqual(@as(f64, 1), result.daily_amplitude_change_c);
    try std.testing.expectApproxEqAbs(
        @as(f64, 1),
        result.diurnal_amplitude_factor,
        1e-8,
    );
    try std.testing.expectApproxEqAbs(@as(f64, 24), state.air_temperature_c, 1e-8);
    try std.testing.expectEqual(@as(f64, 293.15), state.previous_air_temperature_k);
    try std.testing.expectEqual(
        @as(f64, 0.001),
        state.previous_vapor_concentration_m3_per_m3,
    );
}

test "first scene source hour uses annual temperature and zero vapor snapshot" {
    var state = exampleState();
    _ = try apply(&state, .{
        .source_hour = 1,
        .solar_noon_source_hour = 12,
        .first_scene_day = true,
        .maximum_temperature_change_c = 1,
        .minimum_temperature_change_c = 1,
        .initial_mean_annual_air_temperature_k = 285,
    });
    try std.testing.expectEqual(@as(f64, 285), state.previous_air_temperature_k);
    try std.testing.expectEqual(
        @as(f64, 0),
        state.previous_vapor_concentration_m3_per_m3,
    );
}

test "zero temperature changes leave all state untouched" {
    var state = exampleState();
    const before = state;
    const result = try apply(&state, .{
        .source_hour = 12,
        .solar_noon_source_hour = 12,
        .first_scene_day = false,
        .maximum_temperature_change_c = 0,
        .minimum_temperature_change_c = 0,
        .initial_mean_annual_air_temperature_k = 285,
    });
    try std.testing.expect(!result.applied);
    try std.testing.expectEqualDeep(before, state);
}

test "nonphysical resulting temperature rolls back state" {
    var state = exampleState();
    const before = state;
    try std.testing.expectError(
        error.WeatherTemperatureChangeOverflow,
        apply(&state, .{
            .source_hour = 15,
            .solar_noon_source_hour = 12,
            .first_scene_day = false,
            .maximum_temperature_change_c = -1000,
            .minimum_temperature_change_c = -1000,
            .initial_mean_annual_air_temperature_k = 285,
        }),
    );
    try std.testing.expectEqualDeep(before, state);
}
