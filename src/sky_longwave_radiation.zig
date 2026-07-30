const std = @import("std");

pub const Inputs = struct {
    observed_longwave_mj_per_m2_h: f64,
    sky_emissivity_fraction: f64,
    air_temperature_k: f64,
    stefan_boltzmann_mj_per_m2_h_k4: f64 = 2.04e-10,
};

pub const Result = struct {
    longwave_mj_per_m2_h: f64,
    used_observation: bool,
};

/// Exact WTHR longwave selection from wthr.f:265-280.
///
/// The source uses measured longwave only when it is strictly greater than
/// zero. A zero observation is therefore a missing-value signal and invokes
/// the emissivity-temperature calculation.
pub fn select(inputs: Inputs) !Result {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSkyLongwaveInput;
    }
    if (inputs.observed_longwave_mj_per_m2_h < 0 or
        inputs.sky_emissivity_fraction < 0 or
        inputs.sky_emissivity_fraction > 1 or
        inputs.air_temperature_k <= 0 or
        inputs.stefan_boltzmann_mj_per_m2_h_k4 <= 0)
        return error.InvalidSkyLongwaveInput;

    if (inputs.observed_longwave_mj_per_m2_h > 0) {
        return .{
            .longwave_mj_per_m2_h = inputs.observed_longwave_mj_per_m2_h,
            .used_observation = true,
        };
    }
    const calculated =
        inputs.sky_emissivity_fraction *
        inputs.stefan_boltzmann_mj_per_m2_h_k4 *
        std.math.pow(f64, inputs.air_temperature_k, 4);
    if (!std.math.isFinite(calculated))
        return error.SkyLongwaveCalculationOverflow;
    return .{
        .longwave_mj_per_m2_h = calculated,
        .used_observation = false,
    };
}

test "strictly positive observed longwave bypasses calculation" {
    const result = try select(.{
        .observed_longwave_mj_per_m2_h = 0.8,
        .sky_emissivity_fraction = 0.7,
        .air_temperature_k = 290,
    });
    try std.testing.expect(result.used_observation);
    try std.testing.expectEqual(
        @as(f64, 0.8),
        result.longwave_mj_per_m2_h,
    );
}

test "zero observation invokes exact emissivity temperature equation" {
    const result = try select(.{
        .observed_longwave_mj_per_m2_h = 0,
        .sky_emissivity_fraction = 0.85,
        .air_temperature_k = 280,
    });
    const expected =
        0.85 * 2.04e-10 * std.math.pow(f64, @as(f64, 280), 4);
    try std.testing.expect(!result.used_observation);
    try std.testing.expectEqual(expected, result.longwave_mj_per_m2_h);
}

test "runtime radiation coefficient controls calculated result" {
    const result = try select(.{
        .observed_longwave_mj_per_m2_h = 0,
        .sky_emissivity_fraction = 0.97,
        .air_temperature_k = 300,
        .stefan_boltzmann_mj_per_m2_h_k4 = 2.1e-10,
    });
    try std.testing.expectEqual(
        0.97 * 2.1e-10 * std.math.pow(f64, @as(f64, 300), 4),
        result.longwave_mj_per_m2_h,
    );
}

test "invalid measured and derived inputs fail immediately" {
    try std.testing.expectError(
        error.InvalidSkyLongwaveInput,
        select(.{
            .observed_longwave_mj_per_m2_h = -0.1,
            .sky_emissivity_fraction = 0.8,
            .air_temperature_k = 280,
        }),
    );
    try std.testing.expectError(
        error.NonFiniteSkyLongwaveInput,
        select(.{
            .observed_longwave_mj_per_m2_h = 0,
            .sky_emissivity_fraction = std.math.nan(f64),
            .air_temperature_k = 280,
        }),
    );
}

test "derived overflow cannot silently enter energy balance" {
    try std.testing.expectError(
        error.SkyLongwaveCalculationOverflow,
        select(.{
            .observed_longwave_mj_per_m2_h = 0,
            .sky_emissivity_fraction = 1,
            .air_temperature_k = std.math.floatMax(f64),
        }),
    );
}
