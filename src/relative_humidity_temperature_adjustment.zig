const std = @import("std");

pub const State = struct {
    vapor_pressure_kpa: f64,
    saturated_vapor_pressure_kpa: f64,
};

pub const Inputs = struct {
    air_temperature_k: f64,
    altitude_m: f64,
    target_humidity_fraction: f64,
    saturation_prefactor_kpa: f64 = 0.61,
    clausius_clapeyron_coefficient_k: f64 = 5360,
    reference_inverse_temperature_per_k: f64 = 3.661e-3,
    altitude_scale_m: f64 = 7272,
};

/// Exact WTHR adjustment that preserves relative humidity after a temperature
/// change only when the requested humidity multiplier is unity.
pub fn apply(state: *State, inputs: Inputs) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const value = @field(state.*, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteHumidityAdjustmentState;
        if (value < 0) return error.InvalidHumidityAdjustmentState;
    }
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteHumidityAdjustmentInput;
    }
    if (inputs.air_temperature_k <= 0 or inputs.altitude_scale_m <= 0 or
        inputs.saturation_prefactor_kpa <= 0 or
        inputs.clausius_clapeyron_coefficient_k <= 0 or
        inputs.reference_inverse_temperature_per_k <= 0 or
        inputs.target_humidity_fraction < 0)
        return error.InvalidHumidityAdjustmentInput;
    if (inputs.target_humidity_fraction != 1) return;
    if (state.saturated_vapor_pressure_kpa <= 0)
        return error.NonPositivePreviousSaturationPressure;

    const next_saturation_kpa =
        inputs.saturation_prefactor_kpa *
        @exp(inputs.clausius_clapeyron_coefficient_k *
            (inputs.reference_inverse_temperature_per_k -
                1 / inputs.air_temperature_k)) *
        @exp(-inputs.altitude_m / inputs.altitude_scale_m);
    const next_vapor_pressure_kpa =
        state.vapor_pressure_kpa * next_saturation_kpa /
        state.saturated_vapor_pressure_kpa;
    if (!std.math.isFinite(next_saturation_kpa) or
        !std.math.isFinite(next_vapor_pressure_kpa) or
        next_saturation_kpa < 0 or next_vapor_pressure_kpa < 0)
        return error.HumidityTemperatureAdjustmentOverflow;
    state.* = .{
        .vapor_pressure_kpa = next_vapor_pressure_kpa,
        .saturated_vapor_pressure_kpa = next_saturation_kpa,
    };
}

test "unity humidity target preserves relative humidity after warming" {
    var state: State = .{
        .vapor_pressure_kpa = 1,
        .saturated_vapor_pressure_kpa = 2,
    };
    try apply(&state, .{
        .air_temperature_k = 303.15,
        .altitude_m = 0,
        .target_humidity_fraction = 1,
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.5),
        state.vapor_pressure_kpa / state.saturated_vapor_pressure_kpa,
        1e-15,
    );
    try std.testing.expect(state.saturated_vapor_pressure_kpa > 2);
}

test "altitude exponentially lowers recomputed saturation pressure" {
    var sea_level: State = .{
        .vapor_pressure_kpa = 1,
        .saturated_vapor_pressure_kpa = 2,
    };
    var high: State = sea_level;
    try apply(&sea_level, .{
        .air_temperature_k = 293.15,
        .altitude_m = 0,
        .target_humidity_fraction = 1,
    });
    try apply(&high, .{
        .air_temperature_k = 293.15,
        .altitude_m = 7272,
        .target_humidity_fraction = 1,
    });
    try std.testing.expectApproxEqRel(
        sea_level.saturated_vapor_pressure_kpa / @exp(@as(f64, 1)),
        high.saturated_vapor_pressure_kpa,
        1e-15,
    );
}

test "nonunity humidity target leaves state for later scaling block" {
    var state: State = .{
        .vapor_pressure_kpa = 1,
        .saturated_vapor_pressure_kpa = 2,
    };
    const before = state;
    try apply(&state, .{
        .air_temperature_k = 303.15,
        .altitude_m = 100,
        .target_humidity_fraction = 1.2,
    });
    try std.testing.expectEqualDeep(before, state);
}

test "invalid previous saturation rolls back humidity state" {
    var state: State = .{
        .vapor_pressure_kpa = 1,
        .saturated_vapor_pressure_kpa = 0,
    };
    const before = state;
    try std.testing.expectError(
        error.NonPositivePreviousSaturationPressure,
        apply(&state, .{
            .air_temperature_k = 300,
            .altitude_m = 0,
            .target_humidity_fraction = 1,
        }),
    );
    try std.testing.expectEqualDeep(before, state);
}
