const std = @import("std");

pub const Inputs = struct {
    phytotron: bool,
    source_hour: u8,
    solar_noon_hour: f64,
    daylength_h: f64,
    maximum_hourly_shortwave_megajoules_per_m2_h: f64,
    pi_radians: f64 = 3.1416,
};

/// Exact WTHR daily-weather shortwave curve from wthr.f:85-105.
pub fn radiationAtHour(inputs: Inputs) !f64 {
    if (inputs.source_hour < 1 or inputs.source_hour > 24)
        return error.InvalidDailyShortwaveSourceHour;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == bool or field.type == u8) continue;
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteDailyShortwaveInput;
    }
    if (inputs.solar_noon_hour < 1 or inputs.solar_noon_hour >= 25 or
        inputs.daylength_h < 0 or inputs.daylength_h > 24 or
        inputs.maximum_hourly_shortwave_megajoules_per_m2_h < 0 or
        inputs.pi_radians <= 0)
        return error.InvalidDailyShortwaveInput;

    const result = if (inputs.phytotron)
        inputs.maximum_hourly_shortwave_megajoules_per_m2_h / 24.0
    else if (inputs.daylength_h > 0)
        @max(
            0.0,
            inputs.maximum_hourly_shortwave_megajoules_per_m2_h *
                @sin(
                    (@as(f64, @floatFromInt(inputs.source_hour)) -
                        (inputs.solar_noon_hour - inputs.daylength_h / 2.0)) *
                        inputs.pi_radians / inputs.daylength_h,
                ),
        )
    else
        0.0;
    if (!std.math.isFinite(result))
        return error.DailyShortwaveCalculationOverflow;
    return result;
}

fn outdoors() Inputs {
    return .{
        .phytotron = false,
        .source_hour = 12,
        .solar_noon_hour = 12,
        .daylength_h = 12,
        .maximum_hourly_shortwave_megajoules_per_m2_h = 2,
    };
}

test "solar noon reaches exact daily-weather radiation maximum" {
    const result = try radiationAtHour(outdoors());
    try std.testing.expectApproxEqAbs(@as(f64, 2), result, 2e-11);
}

test "pre-sunrise sine values clamp to zero" {
    var inputs = outdoors();
    inputs.source_hour = 1;
    const result = try radiationAtHour(inputs);
    try std.testing.expectEqual(@as(f64, 0), result);
}

test "zero daylength yields zero outdoor radiation" {
    var inputs = outdoors();
    inputs.daylength_h = 0;
    const result = try radiationAtHour(inputs);
    try std.testing.expectEqual(@as(f64, 0), result);
}

test "phytotron distributes maximum radiation uniformly over 24 hours" {
    var inputs = outdoors();
    inputs.phytotron = true;
    inputs.maximum_hourly_shortwave_megajoules_per_m2_h = 24;
    inputs.daylength_h = 0;
    inline for (.{ 1, 12, 24 }) |hour| {
        inputs.source_hour = hour;
        try std.testing.expectEqual(
            @as(f64, 1),
            try radiationAtHour(inputs),
        );
    }
}

test "nonfinite late forcing fails before producing radiation" {
    var inputs = outdoors();
    inputs.maximum_hourly_shortwave_megajoules_per_m2_h = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteDailyShortwaveInput,
        radiationAtHour(inputs),
    );
}
