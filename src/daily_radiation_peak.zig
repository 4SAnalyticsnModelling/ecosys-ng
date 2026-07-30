const std = @import("std");

pub const Inputs = struct {
    climate_type_code: i32,
    daily_shortwave_radiation_mj_per_m2_day: f64,
    daylength_h: f64,
    hourly_curve_shape_factor: f64 = 0.658,
};

/// Exact DAY RMAX preparation from day.f:238-256.
pub fn derive(inputs: Inputs) !f64 {
    inline for (.{
        inputs.daily_shortwave_radiation_mj_per_m2_day,
        inputs.daylength_h,
        inputs.hourly_curve_shape_factor,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteDailyRadiationPeakInput;
    if (inputs.daily_shortwave_radiation_mj_per_m2_day < 0 or
        inputs.daylength_h < 0 or inputs.daylength_h > 24 or
        inputs.hourly_curve_shape_factor <= 0)
        return error.InvalidDailyRadiationPeakInput;

    const result =
        if (inputs.climate_type_code >= -1)
            if (inputs.daylength_h > 0)
                inputs.daily_shortwave_radiation_mj_per_m2_day /
                    (inputs.daylength_h *
                        inputs.hourly_curve_shape_factor)
            else
                0
        else
            inputs.daily_shortwave_radiation_mj_per_m2_day;
    if (!std.math.isFinite(result))
        return error.DailyRadiationPeakOverflow;
    return result;
}

test "outdoor daily radiation is normalized by daylength and shape" {
    try std.testing.expectApproxEqAbs(
        @as(f64, 12) / (12 * 0.658),
        try derive(.{
            .climate_type_code = 0,
            .daily_shortwave_radiation_mj_per_m2_day = 12,
            .daylength_h = 12,
        }),
        1e-15,
    );
}

test "open-top chamber boundary follows outdoor normalization" {
    const result = try derive(.{
        .climate_type_code = -1,
        .daily_shortwave_radiation_mj_per_m2_day = 10,
        .daylength_h = 10,
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, 1) / 0.658,
        result,
        1e-15,
    );
}

test "phytotron code keeps source radiation value for later division" {
    const result = try derive(.{
        .climate_type_code = -2,
        .daily_shortwave_radiation_mj_per_m2_day = 24,
        .daylength_h = 0,
    });
    try std.testing.expectEqual(@as(f64, 24), result);
}

test "zero outdoor daylength produces zero radiation peak" {
    const result = try derive(.{
        .climate_type_code = 0,
        .daily_shortwave_radiation_mj_per_m2_day = 5,
        .daylength_h = 0,
    });
    try std.testing.expectEqual(@as(f64, 0), result);
}

test "nonfinite or invalid shape input fails before division" {
    try std.testing.expectError(
        error.NonFiniteDailyRadiationPeakInput,
        derive(.{
            .climate_type_code = 0,
            .daily_shortwave_radiation_mj_per_m2_day = 5,
            .daylength_h = std.math.nan(f64),
        }),
    );
    try std.testing.expectError(
        error.InvalidDailyRadiationPeakInput,
        derive(.{
            .climate_type_code = 0,
            .daily_shortwave_radiation_mj_per_m2_day = 5,
            .daylength_h = 12,
            .hourly_curve_shape_factor = 0,
        }),
    );
}
