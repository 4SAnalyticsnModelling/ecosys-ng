const std = @import("std");

pub const Inputs = struct {
    execution_day: u32,
    first_execution_day: u32,
    source_hour: u8,
    current_horizontal_shortwave_megajoules_per_m2_h: f64,
    previous_horizontal_shortwave_megajoules_per_m2_h: f64,
    current_air_temperature_k: f64,
    previous_air_temperature_k: f64,
    current_atmospheric_vapor_concentration_m3_per_m3: f64,
    previous_atmospheric_vapor_concentration_m3_per_m3: f64,
};

pub const Change = struct {
    horizontal_shortwave_change_megajoules_per_m2_h: f64,
    air_temperature_change_k: f64,
    atmospheric_vapor_concentration_change_m3_per_m3: f64,
};

/// Exact WTHR hourly forcing changes from wthr.f:515-526.
///
/// DRADN is the change in incoming horizontal shortwave RADN. It is not net
/// radiation. The first source hour of the first execution day clears all
/// changes irrespective of the supplied previous snapshot.
pub fn calculate(inputs: Inputs) !Change {
    if (inputs.execution_day == 0 or inputs.first_execution_day == 0 or
        inputs.source_hour < 1 or inputs.source_hour > 24)
        return error.InvalidHourlyForcingChangeTime;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == u32 or field.type == u8) continue;
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteHourlyForcingChangeInput;
    }
    if (inputs.current_horizontal_shortwave_megajoules_per_m2_h < 0 or
        inputs.previous_horizontal_shortwave_megajoules_per_m2_h < 0 or
        inputs.current_air_temperature_k <= 0 or
        inputs.previous_air_temperature_k <= 0 or
        inputs.current_atmospheric_vapor_concentration_m3_per_m3 < 0 or
        inputs.previous_atmospheric_vapor_concentration_m3_per_m3 < 0)
        return error.InvalidHourlyForcingChangeInput;

    if (inputs.execution_day == inputs.first_execution_day and
        inputs.source_hour == 1)
        return std.mem.zeroes(Change);

    const result: Change = .{
        .horizontal_shortwave_change_megajoules_per_m2_h = inputs.current_horizontal_shortwave_megajoules_per_m2_h -
            inputs.previous_horizontal_shortwave_megajoules_per_m2_h,
        .air_temperature_change_k = inputs.current_air_temperature_k -
            inputs.previous_air_temperature_k,
        .atmospheric_vapor_concentration_change_m3_per_m3 = inputs.current_atmospheric_vapor_concentration_m3_per_m3 -
            inputs.previous_atmospheric_vapor_concentration_m3_per_m3,
    };
    inline for (@typeInfo(Change).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.HourlyForcingChangeOverflow;
    return result;
}

fn exampleInputs() Inputs {
    return .{
        .execution_day = 100,
        .first_execution_day = 100,
        .source_hour = 2,
        .current_horizontal_shortwave_megajoules_per_m2_h = 0.8,
        .previous_horizontal_shortwave_megajoules_per_m2_h = 0.3,
        .current_air_temperature_k = 290,
        .previous_air_temperature_k = 288,
        .current_atmospheric_vapor_concentration_m3_per_m3 = 0.0012,
        .previous_atmospheric_vapor_concentration_m3_per_m3 = 0.001,
    };
}

test "ordinary hour calculates exact signed forcing differences" {
    const result = try calculate(exampleInputs());
    try std.testing.expectEqual(
        @as(f64, 0.5),
        result.horizontal_shortwave_change_megajoules_per_m2_h,
    );
    try std.testing.expectEqual(@as(f64, 2), result.air_temperature_change_k);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.0002),
        result.atmospheric_vapor_concentration_change_m3_per_m3,
        1e-18,
    );
}

test "changes remain signed when hourly forcing decreases" {
    var inputs = exampleInputs();
    inputs.current_horizontal_shortwave_megajoules_per_m2_h = 0.1;
    inputs.current_air_temperature_k = 285;
    inputs.current_atmospheric_vapor_concentration_m3_per_m3 = 0.0008;
    const result = try calculate(inputs);
    try std.testing.expect(result.horizontal_shortwave_change_megajoules_per_m2_h < 0);
    try std.testing.expect(result.air_temperature_change_k < 0);
    try std.testing.expect(
        result.atmospheric_vapor_concentration_change_m3_per_m3 < 0,
    );
}

test "first execution hour clears all changes" {
    var inputs = exampleInputs();
    inputs.source_hour = 1;
    const result = try calculate(inputs);
    try std.testing.expectEqualDeep(std.mem.zeroes(Change), result);
}

test "hour one on later execution day uses differences" {
    var inputs = exampleInputs();
    inputs.execution_day = 101;
    inputs.source_hour = 1;
    const result = try calculate(inputs);
    try std.testing.expectEqual(
        @as(f64, 0.5),
        result.horizontal_shortwave_change_megajoules_per_m2_h,
    );
}

test "nonfinite late previous carrier fails immediately" {
    var inputs = exampleInputs();
    inputs.previous_atmospheric_vapor_concentration_m3_per_m3 =
        std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteHourlyForcingChangeInput,
        calculate(inputs),
    );
}
