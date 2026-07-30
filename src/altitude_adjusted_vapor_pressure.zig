const std = @import("std");

pub const Inputs = struct {
    unconstrained_ambient_vapor_pressure_kpa: f64,
    air_temperature_k: f64,
    altitude_m: f64,
    saturation_pressure_at_reference_kpa: f64 = 0.61,
    clausius_clapeyron_temperature_k: f64 = 5360,
    reference_inverse_temperature_per_k: f64 = 3.661e-3,
    altitude_scale_m: f64 = 7272,
};

pub const Result = struct {
    saturated_vapor_pressure_kpa: f64,
    ambient_vapor_pressure_kpa: f64,
};

/// Exact WTHR VPS and VPK cap from wthr.f:122-139 and 187-200.
pub fn constrain(inputs: Inputs) !Result {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteAltitudeVaporPressureInput;
    }
    if (inputs.unconstrained_ambient_vapor_pressure_kpa < 0 or
        inputs.air_temperature_k <= 0 or
        inputs.saturation_pressure_at_reference_kpa <= 0 or
        inputs.clausius_clapeyron_temperature_k <= 0 or
        inputs.reference_inverse_temperature_per_k <= 0 or
        inputs.altitude_scale_m <= 0)
        return error.InvalidAltitudeVaporPressureInput;
    const saturation =
        inputs.saturation_pressure_at_reference_kpa *
        @exp(
            inputs.clausius_clapeyron_temperature_k *
                (inputs.reference_inverse_temperature_per_k -
                    1 / inputs.air_temperature_k),
        ) *
        @exp(-inputs.altitude_m / inputs.altitude_scale_m);
    if (!std.math.isFinite(saturation))
        return error.AltitudeVaporPressureOverflow;
    return .{
        .saturated_vapor_pressure_kpa = saturation,
        .ambient_vapor_pressure_kpa = @min(
            inputs.unconstrained_ambient_vapor_pressure_kpa,
            saturation,
        ),
    };
}

test "ambient pressure is capped at altitude-adjusted saturation" {
    const result = try constrain(.{
        .unconstrained_ambient_vapor_pressure_kpa = 20,
        .air_temperature_k = 283.15,
        .altitude_m = 727.2,
    });
    try std.testing.expectEqual(
        result.saturated_vapor_pressure_kpa,
        result.ambient_vapor_pressure_kpa,
    );
}

test "ambient pressure below saturation is preserved" {
    const result = try constrain(.{
        .unconstrained_ambient_vapor_pressure_kpa = 0.2,
        .air_temperature_k = 293.15,
        .altitude_m = 0,
    });
    try std.testing.expectEqual(
        @as(f64, 0.2),
        result.ambient_vapor_pressure_kpa,
    );
}

test "higher altitude lowers saturation at fixed temperature" {
    const low = try constrain(.{
        .unconstrained_ambient_vapor_pressure_kpa = 0,
        .air_temperature_k = 280,
        .altitude_m = 0,
    });
    const high = try constrain(.{
        .unconstrained_ambient_vapor_pressure_kpa = 0,
        .air_temperature_k = 280,
        .altitude_m = 2000,
    });
    try std.testing.expect(
        high.saturated_vapor_pressure_kpa <
            low.saturated_vapor_pressure_kpa,
    );
}

test "runtime scientific coefficients are honored" {
    const result = try constrain(.{
        .unconstrained_ambient_vapor_pressure_kpa = 10,
        .air_temperature_k = 300,
        .altitude_m = 100,
        .saturation_pressure_at_reference_kpa = 0.7,
        .clausius_clapeyron_temperature_k = 5000,
        .reference_inverse_temperature_per_k = 0.0035,
        .altitude_scale_m = 7000,
    });
    const expected =
        0.7 * @exp(5000 * (0.0035 - 1.0 / 300.0)) *
        @exp(-100.0 / 7000.0);
    try std.testing.expectApproxEqAbs(
        expected,
        result.saturated_vapor_pressure_kpa,
        2e-15,
    );
}

test "nonfinite late coefficient fails before returning pressure" {
    try std.testing.expectError(
        error.NonFiniteAltitudeVaporPressureInput,
        constrain(.{
            .unconstrained_ambient_vapor_pressure_kpa = 1,
            .air_temperature_k = 280,
            .altitude_m = 0,
            .altitude_scale_m = std.math.nan(f64),
        }),
    );
}
