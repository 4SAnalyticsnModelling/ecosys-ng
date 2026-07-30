const std = @import("std");

pub const Cause = packed struct {
    scheduled_at_solar_noon: bool = false,
    canopy_temperature: bool = false,
    upper_material_temperature: bool = false,
    deeper_soil_temperature: bool = false,
};

pub const Inputs = struct {
    source_hour: u8,
    solar_noon_source_hour: f64,
    scheduled_fire_today: bool,
    canopy_air_temperature_k: f64,
    material_temperature_k: []const f64,
    /// Inclusive NU-equivalent index. Indices through this value use the
    /// surface-fire threshold; deeper indices use the lower-soil threshold.
    upper_material_last_index: usize,
    surface_combustion_threshold_k: f64,
    deeper_soil_combustion_threshold_k: f64,
};

pub const Result = struct {
    active: bool,
    cause: Cause,
};

/// Evaluates WTHR fire activation and commits the cell flag only after every
/// runtime layer and derived source-hour value has passed validation.
pub fn evaluateAndCommit(
    fire_active_this_hour: *bool,
    inputs: Inputs,
) !Result {
    const result = try evaluate(inputs);
    fire_active_this_hour.* = result.active;
    return result;
}

pub fn evaluate(inputs: Inputs) !Result {
    if (inputs.source_hour < 1 or inputs.source_hour > 24)
        return error.InvalidFireSourceHour;
    inline for (.{
        inputs.solar_noon_source_hour,
        inputs.canopy_air_temperature_k,
        inputs.surface_combustion_threshold_k,
        inputs.deeper_soil_combustion_threshold_k,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteFireActivationInput;
    if (inputs.solar_noon_source_hour < 1 or
        inputs.solar_noon_source_hour >= 25 or
        inputs.canopy_air_temperature_k <= 0 or
        inputs.surface_combustion_threshold_k <= 0 or
        inputs.deeper_soil_combustion_threshold_k <= 0)
        return error.InvalidFireActivationInput;
    if (inputs.material_temperature_k.len == 0 or
        inputs.upper_material_last_index >= inputs.material_temperature_k.len)
        return error.FireMaterialDimensionMismatch;

    var result: Result = .{
        .active = false,
        .cause = .{},
    };
    const solar_noon_integer_hour: u8 =
        @intFromFloat(@trunc(inputs.solar_noon_source_hour));
    result.cause.scheduled_at_solar_noon =
        inputs.scheduled_fire_today and
        inputs.source_hour == solar_noon_integer_hour;
    result.cause.canopy_temperature =
        inputs.canopy_air_temperature_k >
        inputs.surface_combustion_threshold_k;

    for (inputs.material_temperature_k, 0..) |temperature_k, material| {
        if (!std.math.isFinite(temperature_k) or temperature_k <= 0)
            return error.InvalidFireMaterialTemperature;
        if (material <= inputs.upper_material_last_index) {
            result.cause.upper_material_temperature =
                result.cause.upper_material_temperature or
                temperature_k > inputs.surface_combustion_threshold_k;
        } else {
            result.cause.deeper_soil_temperature =
                result.cause.deeper_soil_temperature or
                temperature_k > inputs.deeper_soil_combustion_threshold_k;
        }
    }
    result.active =
        result.cause.scheduled_at_solar_noon or
        result.cause.canopy_temperature or
        result.cause.upper_material_temperature or
        result.cause.deeper_soil_temperature;
    return result;
}

fn baseInputs(temperatures_k: []const f64) Inputs {
    return .{
        .source_hour = 12,
        .solar_noon_source_hour = 12.9,
        .scheduled_fire_today = false,
        .canopy_air_temperature_k = 300,
        .material_temperature_k = temperatures_k,
        .upper_material_last_index = 1,
        .surface_combustion_threshold_k = 500,
        .deeper_soil_combustion_threshold_k = 700,
    };
}

test "operation 22 activates only at truncated solar-noon source hour" {
    const temperatures = [_]f64{ 300, 310, 320, 330 };
    var inputs = baseInputs(&temperatures);
    inputs.scheduled_fire_today = true;
    const noon = try evaluate(inputs);
    try std.testing.expect(noon.active);
    try std.testing.expect(noon.cause.scheduled_at_solar_noon);

    inputs.source_hour = 13;
    const after_noon = try evaluate(inputs);
    try std.testing.expect(!after_noon.active);
}

test "runtime upper and deeper temperatures use distinct strict thresholds" {
    const temperatures = [_]f64{ 500, 499, 700, 701, 300, 250 };
    var inputs = baseInputs(&temperatures);
    inputs.upper_material_last_index = 1;
    const result = try evaluate(inputs);
    try std.testing.expect(result.active);
    try std.testing.expect(!result.cause.upper_material_temperature);
    try std.testing.expect(result.cause.deeper_soil_temperature);

    inputs.canopy_air_temperature_k = 501;
    const canopy = try evaluate(inputs);
    try std.testing.expect(canopy.cause.canopy_temperature);
}

test "invalid late runtime layer leaves prior fire flag unchanged" {
    const temperatures =
        [_]f64{ 300, 310, 320, std.math.nan(f64), 330 };
    var active = true;
    try std.testing.expectError(
        error.InvalidFireMaterialTemperature,
        evaluateAndCommit(&active, baseInputs(&temperatures)),
    );
    try std.testing.expect(active);
}

test "thermal activation is independent of scheduled disturbance" {
    const temperatures = [_]f64{ 300, 501, 320 };
    var inputs = baseInputs(&temperatures);
    inputs.source_hour = 4;
    inputs.scheduled_fire_today = false;
    const result = try evaluate(inputs);
    try std.testing.expect(result.active);
    try std.testing.expect(result.cause.upper_material_temperature);
    try std.testing.expect(!result.cause.scheduled_at_solar_noon);
}
