const std = @import("std");

pub const Segment = enum {
    previous_day_night,
    current_day,
    next_day_night,
};

pub const Inputs = struct {
    source_hour: u8,
    solar_noon_hour: f64,
    daylength_h: f64,
    previous_average: f64,
    current_average: f64,
    next_average: f64,
    previous_amplitude: f64,
    current_amplitude: f64,
    next_amplitude: f64,
    pi_radians: f64 = 3.1416,
    half_pi_radians: f64 = 1.5708,
};

pub const Result = struct {
    value: f64,
    segment: Segment,
};

/// Exact WTHR three-segment temperature/vapor curve from wthr.f:107-137.
pub fn evaluate(inputs: Inputs) !Result {
    if (inputs.source_hour < 1 or inputs.source_hour > 24)
        return error.InvalidDiurnalCurveSourceHour;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == u8) continue;
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteDiurnalCurveInput;
    }
    if (inputs.solar_noon_hour < 1 or inputs.solar_noon_hour >= 25 or
        inputs.daylength_h < 0 or inputs.daylength_h > 24 or
        inputs.pi_radians <= 0 or inputs.half_pi_radians <= 0)
        return error.InvalidDiurnalCurveInput;

    const hour: f64 = @floatFromInt(inputs.source_hour);
    const sunrise_hour = inputs.solar_noon_hour - inputs.daylength_h / 2;
    const night_denominator =
        inputs.solar_noon_hour + 9 - inputs.daylength_h / 2;
    if (night_denominator == 0)
        return error.SingularDiurnalCurveGeometry;

    const result: Result = if (hour < sunrise_hour)
        .{
            .value = inputs.previous_average +
                inputs.previous_amplitude *
                    @sin(
                        (hour + inputs.solar_noon_hour - 3) *
                            inputs.pi_radians / night_denominator +
                            inputs.half_pi_radians,
                    ),
            .segment = .previous_day_night,
        }
    else if (hour > inputs.solar_noon_hour + 3)
        .{
            .value = inputs.next_average +
                inputs.next_amplitude *
                    @sin(
                        (hour - inputs.solar_noon_hour - 3) *
                            inputs.pi_radians / night_denominator +
                            inputs.half_pi_radians,
                    ),
            .segment = .next_day_night,
        }
    else
        .{
            .value = inputs.current_average +
                inputs.current_amplitude *
                    @sin(
                        (hour - sunrise_hour) * inputs.pi_radians /
                            (3 + inputs.daylength_h / 2) -
                            inputs.half_pi_radians,
                    ),
            .segment = .current_day,
        };
    if (!std.math.isFinite(result.value))
        return error.DiurnalCurveCalculationOverflow;
    return result;
}

fn exampleInputs() Inputs {
    return .{
        .source_hour = 5,
        .solar_noon_hour = 12,
        .daylength_h = 12,
        .previous_average = 4,
        .current_average = 10,
        .next_average = 6,
        .previous_amplitude = 2,
        .current_amplitude = 8,
        .next_amplitude = 3,
    };
}

test "hour before sunrise uses previous-day night curve" {
    const result = try evaluate(exampleInputs());
    try std.testing.expectEqual(
        Segment.previous_day_night,
        result.segment,
    );
}

test "sunrise equality uses current-day curve" {
    var inputs = exampleInputs();
    inputs.source_hour = 6;
    const result = try evaluate(inputs);
    try std.testing.expectEqual(Segment.current_day, result.segment);
    try std.testing.expectApproxEqAbs(
        inputs.current_average - inputs.current_amplitude,
        result.value,
        1e-10,
    );
}

test "three hours after noon equality remains current-day curve" {
    var inputs = exampleInputs();
    inputs.source_hour = 15;
    const at_boundary = try evaluate(inputs);
    try std.testing.expectEqual(Segment.current_day, at_boundary.segment);

    inputs.source_hour = 16;
    const after_boundary = try evaluate(inputs);
    try std.testing.expectEqual(
        Segment.next_day_night,
        after_boundary.segment,
    );
}

test "same curve accepts vapor-pressure averages and amplitudes" {
    var inputs = exampleInputs();
    inputs.source_hour = 12;
    inputs.previous_average = 0.8;
    inputs.current_average = 1.2;
    inputs.next_average = 1;
    inputs.previous_amplitude = 0.1;
    inputs.current_amplitude = 0.3;
    inputs.next_amplitude = 0.2;
    const result = try evaluate(inputs);
    try std.testing.expect(std.math.isFinite(result.value));
    try std.testing.expectEqual(Segment.current_day, result.segment);
}

test "singular and nonfinite geometry fail immediately" {
    var inputs = exampleInputs();
    inputs.solar_noon_hour = 3;
    inputs.daylength_h = 24;
    try std.testing.expectError(
        error.SingularDiurnalCurveGeometry,
        evaluate(inputs),
    );
    inputs = exampleInputs();
    inputs.next_amplitude = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteDiurnalCurveInput,
        evaluate(inputs),
    );
}
