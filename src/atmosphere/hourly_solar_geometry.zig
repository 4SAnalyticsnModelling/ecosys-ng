const std = @import("std");

pub const Inputs = struct {
    outdoors: bool,
    source_hour: u8,
    solar_noon_hour: f64,
    latitude_declination_sine_product: f64,
    latitude_declination_cosine_product: f64,
    incoming_shortwave_megajoules_per_m2_h: f64,
    solar_hour_angle_radians_per_h: f64 = 0.2618,
    solar_constant_megajoules_per_m2_h: f64 = 4.896,
};

pub const Result = struct {
    current_solar_angle_sine: f64,
    next_hour_solar_angle_sine: f64,
    extraterrestrial_shortwave_megajoules_per_m2_h: f64,
    capped_incoming_shortwave_megajoules_per_m2_h: f64,
};

/// Exact WTHR hourly solar geometry from wthr.f:220-237 and 252-263.
pub fn calculate(inputs: Inputs) !Result {
    if (inputs.source_hour < 1 or inputs.source_hour > 24)
        return error.InvalidSolarGeometrySourceHour;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == bool or field.type == u8) continue;
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteHourlySolarGeometryInput;
    }
    if (inputs.solar_noon_hour < 1 or inputs.solar_noon_hour >= 25 or
        inputs.latitude_declination_sine_product < -1 or
        inputs.latitude_declination_sine_product > 1 or
        inputs.latitude_declination_cosine_product < -1 or
        inputs.latitude_declination_cosine_product > 1 or
        inputs.incoming_shortwave_megajoules_per_m2_h < 0 or
        inputs.solar_hour_angle_radians_per_h <= 0 or
        inputs.solar_constant_megajoules_per_m2_h <= 0)
        return error.InvalidHourlySolarGeometryInput;

    if (!inputs.outdoors) {
        return .{
            .current_solar_angle_sine = if (inputs.incoming_shortwave_megajoules_per_m2_h <= 0) 0 else 1,
            .next_hour_solar_angle_sine = 1,
            .extraterrestrial_shortwave_megajoules_per_m2_h = 0,
            .capped_incoming_shortwave_megajoules_per_m2_h = inputs.incoming_shortwave_megajoules_per_m2_h,
        };
    }

    const hour: f64 = @floatFromInt(inputs.source_hour);
    var current_sine = @max(
        0.0,
        inputs.latitude_declination_sine_product +
            inputs.latitude_declination_cosine_product *
                @cos(
                    inputs.solar_hour_angle_radians_per_h *
                        (inputs.solar_noon_hour - (hour - 0.5)),
                ),
    );
    const next_sine = @max(
        0.0,
        inputs.latitude_declination_sine_product +
            inputs.latitude_declination_cosine_product *
                @cos(
                    inputs.solar_hour_angle_radians_per_h *
                        (inputs.solar_noon_hour - (hour + 0.5)),
                ),
    );
    if (inputs.incoming_shortwave_megajoules_per_m2_h <= 0)
        current_sine = 0;
    const extraterrestrial =
        inputs.solar_constant_megajoules_per_m2_h * @max(0.0, current_sine);
    const capped = @min(
        extraterrestrial,
        inputs.incoming_shortwave_megajoules_per_m2_h,
    );
    if (!std.math.isFinite(extraterrestrial))
        return error.HourlySolarGeometryOverflow;
    return .{
        .current_solar_angle_sine = current_sine,
        .next_hour_solar_angle_sine = next_sine,
        .extraterrestrial_shortwave_megajoules_per_m2_h = extraterrestrial,
        .capped_incoming_shortwave_megajoules_per_m2_h = capped,
    };
}

fn outdoorInputs() Inputs {
    return .{
        .outdoors = true,
        .source_hour = 12,
        .solar_noon_hour = 12,
        .latitude_declination_sine_product = 0.1,
        .latitude_declination_cosine_product = 0.8,
        .incoming_shortwave_megajoules_per_m2_h = 10,
    };
}

test "outdoor incoming shortwave is capped by extraterrestrial radiation" {
    const result = try calculate(outdoorInputs());
    try std.testing.expect(result.current_solar_angle_sine > 0);
    try std.testing.expectEqual(
        result.extraterrestrial_shortwave_megajoules_per_m2_h,
        result.capped_incoming_shortwave_megajoules_per_m2_h,
    );
}

test "zero incoming radiation clears current but not next solar sine" {
    var inputs = outdoorInputs();
    inputs.incoming_shortwave_megajoules_per_m2_h = 0;
    const result = try calculate(inputs);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.current_solar_angle_sine,
    );
    try std.testing.expect(result.next_hour_solar_angle_sine > 0);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.extraterrestrial_shortwave_megajoules_per_m2_h,
    );
}

test "night geometry clamps both outdoor solar sines to zero" {
    var inputs = outdoorInputs();
    inputs.source_hour = 1;
    inputs.latitude_declination_sine_product = -0.5;
    inputs.latitude_declination_cosine_product = 0.1;
    const result = try calculate(inputs);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.current_solar_angle_sine,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.next_hour_solar_angle_sine,
    );
}

test "phytotron preserves incoming radiation and source angle flags" {
    var inputs = outdoorInputs();
    inputs.outdoors = false;
    const lit = try calculate(inputs);
    try std.testing.expectEqual(@as(f64, 1), lit.current_solar_angle_sine);
    try std.testing.expectEqual(@as(f64, 1), lit.next_hour_solar_angle_sine);
    try std.testing.expectEqual(
        inputs.incoming_shortwave_megajoules_per_m2_h,
        lit.capped_incoming_shortwave_megajoules_per_m2_h,
    );

    inputs.incoming_shortwave_megajoules_per_m2_h = 0;
    const dark = try calculate(inputs);
    try std.testing.expectEqual(@as(f64, 0), dark.current_solar_angle_sine);
    try std.testing.expectEqual(@as(f64, 1), dark.next_hour_solar_angle_sine);
}

test "invalid late geometry carrier fails immediately" {
    var inputs = outdoorInputs();
    inputs.latitude_declination_cosine_product = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteHourlySolarGeometryInput,
        calculate(inputs),
    );
}
