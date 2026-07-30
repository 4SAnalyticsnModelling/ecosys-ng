const std = @import("std");

/// Inputs for the HOUR1 subhourly atmospheric interpolation at one grid cell.
/// Temperatures are K, vapor pressures are kPa, and fractions are dimensionless.
pub const Inputs = struct {
    is_start_day: bool,
    source_hour: u8,
    current_air_temperature_k: f64,
    current_vapor_pressure_kpa: f64,
    base_air_temperature_k: f64,
    base_vapor_pressure_kpa: f64,
    day_fraction: f64,
    hour_fraction: f64,
    air_temperature_change_k: f64,
    vapor_pressure_change_kpa: f64,
};

pub const AtmosphericState = struct {
    air_temperature_k: f64,
    vapor_pressure_kpa: f64,
};

pub const InterpolationError = error{
    NonFiniteInput,
    InvalidAirTemperature,
    InvalidVaporPressure,
};

/// Translates HOUR1 lines 2395-2401 (TKAM and VPAM assignment).
pub fn interpolate(inputs: Inputs) InterpolationError!AtmosphericState {
    inline for (std.meta.fields(Inputs)) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) {
            return error.NonFiniteInput;
        }
    }

    const state = if (inputs.is_start_day and inputs.source_hour == 1)
        AtmosphericState{
            .air_temperature_k = inputs.current_air_temperature_k,
            .vapor_pressure_kpa = inputs.current_vapor_pressure_kpa,
        }
    else
        AtmosphericState{
            .air_temperature_k = inputs.base_air_temperature_k +
                inputs.day_fraction * inputs.hour_fraction * inputs.air_temperature_change_k,
            .vapor_pressure_kpa = inputs.base_vapor_pressure_kpa +
                inputs.day_fraction * inputs.hour_fraction * inputs.vapor_pressure_change_kpa,
        };

    if (!std.math.isFinite(state.air_temperature_k) or state.air_temperature_k <= 0.0) {
        return error.InvalidAirTemperature;
    }
    if (!std.math.isFinite(state.vapor_pressure_kpa) or state.vapor_pressure_kpa < 0.0) {
        return error.InvalidVaporPressure;
    }
    return state;
}

test "start day first source hour uses current atmospheric state" {
    const state = try interpolate(.{
        .is_start_day = true,
        .source_hour = 1,
        .current_air_temperature_k = 289.0,
        .current_vapor_pressure_kpa = 1.7,
        .base_air_temperature_k = 250.0,
        .base_vapor_pressure_kpa = 0.2,
        .day_fraction = 0.75,
        .hour_fraction = 0.5,
        .air_temperature_change_k = 20.0,
        .vapor_pressure_change_kpa = 2.0,
    });

    try std.testing.expectEqual(@as(f64, 289.0), state.air_temperature_k);
    try std.testing.expectEqual(@as(f64, 1.7), state.vapor_pressure_kpa);
}

test "other hours preserve legacy interpolation operation order" {
    const state = try interpolate(.{
        .is_start_day = false,
        .source_hour = 1,
        .current_air_temperature_k = 289.0,
        .current_vapor_pressure_kpa = 1.7,
        .base_air_temperature_k = 280.0,
        .base_vapor_pressure_kpa = 1.0,
        .day_fraction = 0.5,
        .hour_fraction = 0.25,
        .air_temperature_change_k = 8.0,
        .vapor_pressure_change_kpa = 0.8,
    });

    try std.testing.expectEqual(@as(f64, 281.0), state.air_temperature_k);
    try std.testing.expectApproxEqAbs(@as(f64, 1.1), state.vapor_pressure_kpa, 1.0e-15);
}

test "non-finite atmospheric input fails explicitly" {
    try std.testing.expectError(error.NonFiniteInput, interpolate(.{
        .is_start_day = false,
        .source_hour = 2,
        .current_air_temperature_k = 289.0,
        .current_vapor_pressure_kpa = 1.7,
        .base_air_temperature_k = std.math.nan(f64),
        .base_vapor_pressure_kpa = 1.0,
        .day_fraction = 0.5,
        .hour_fraction = 0.25,
        .air_temperature_change_k = 8.0,
        .vapor_pressure_change_kpa = 0.8,
    }));
}
