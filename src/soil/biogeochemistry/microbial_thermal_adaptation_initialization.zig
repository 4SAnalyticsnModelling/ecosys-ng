const std = @import("std");

pub const Parameters = struct {
    branch_temperature_c: f64,
    minimum_effective_temperature_c: f64,
    maximum_effective_temperature_c: f64,
    cold_response_c_per_c: f64,
    warm_response_c_per_c: f64,
};

/// Exact source-order translation of legacy `STARTS` lines 510--514.
pub fn derive(
    thermal_adaptation_offset_c: []f64,
    mean_annual_soil_temperature_c: []const f64,
    parameters: Parameters,
) !void {
    const cell_count = mean_annual_soil_temperature_c.len;
    if (cell_count == 0 or thermal_adaptation_offset_c.len != cell_count)
        return error.MicrobialThermalAdaptationDimensionMismatch;
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteMicrobialThermalAdaptationParameter;
    }
    if (parameters.minimum_effective_temperature_c >
        parameters.branch_temperature_c or
        parameters.branch_temperature_c >
            parameters.maximum_effective_temperature_c or
        parameters.cold_response_c_per_c < 0 or
        parameters.warm_response_c_per_c < 0)
    {
        return error.InvalidMicrobialThermalAdaptationParameter;
    }
    for (mean_annual_soil_temperature_c) |temperature_c| {
        if (!std.math.isFinite(temperature_c))
            return error.NonFiniteMeanAnnualSoilTemperature;
        const candidate =
            if (temperature_c <= parameters.branch_temperature_c)
                parameters.cold_response_c_per_c *
                    (parameters.branch_temperature_c -
                        @max(
                            parameters.minimum_effective_temperature_c,
                            temperature_c,
                        ))
            else
                parameters.warm_response_c_per_c *
                    (parameters.branch_temperature_c -
                        @min(
                            parameters.maximum_effective_temperature_c,
                            temperature_c,
                        ));
        if (!std.math.isFinite(candidate))
            return error.MicrobialThermalAdaptationOverflow;
    }

    for (mean_annual_soil_temperature_c, thermal_adaptation_offset_c) |
        temperature_c,
        *offset_c,
    | {
        if (temperature_c <= parameters.branch_temperature_c) {
            offset_c.* = parameters.cold_response_c_per_c *
                (parameters.branch_temperature_c -
                    @max(
                        parameters.minimum_effective_temperature_c,
                        temperature_c,
                    ));
        } else {
            offset_c.* = parameters.warm_response_c_per_c *
                (parameters.branch_temperature_c -
                    @min(
                        parameters.maximum_effective_temperature_c,
                        temperature_c,
                    ));
        }
    }
}

fn sourceParameters() Parameters {
    return .{
        .branch_temperature_c = 15.0,
        .minimum_effective_temperature_c = 0.0,
        .maximum_effective_temperature_c = 30.0,
        .cold_response_c_per_c = 0.333,
        .warm_response_c_per_c = 0.167,
    };
}

test "STARTS thermal adaptation reproduces both source branches and clamps" {
    var offsets = [_]f64{0.0} ** 7;
    try derive(
        &offsets,
        &.{ -10.0, 0.0, 10.0, 15.0, 20.0, 30.0, 40.0 },
        sourceParameters(),
    );

    try std.testing.expectApproxEqAbs(@as(f64, 4.995), offsets[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 4.995), offsets[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.665), offsets[2], 1e-15);
    try std.testing.expectEqual(@as(f64, 0.0), offsets[3]);
    try std.testing.expectApproxEqAbs(@as(f64, -0.835), offsets[4], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -2.505), offsets[5], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -2.505), offsets[6], 1e-15);
}

test "runtime parameters control thermal adaptation without fixed constants" {
    var offsets = [_]f64{ 9.0, 9.0, 9.0 };
    try derive(&offsets, &.{ 5.0, 12.0, 25.0 }, .{
        .branch_temperature_c = 12.0,
        .minimum_effective_temperature_c = 2.0,
        .maximum_effective_temperature_c = 22.0,
        .cold_response_c_per_c = 0.4,
        .warm_response_c_per_c = 0.2,
    });
    for (offsets, [_]f64{ 2.8, 0.0, -2.0 }) |actual, expected| {
        try std.testing.expectApproxEqAbs(expected, actual, 1e-15);
    }
}

test "late nonfinite temperature preserves every offset" {
    var offsets = [_]f64{ 7.0, 8.0 };
    const before = offsets;
    try std.testing.expectError(
        error.NonFiniteMeanAnnualSoilTemperature,
        derive(
            &offsets,
            &.{ 10.0, std.math.nan(f64) },
            sourceParameters(),
        ),
    );
    try std.testing.expectEqualSlices(f64, &before, &offsets);
}
