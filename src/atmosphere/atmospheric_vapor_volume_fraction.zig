const std = @import("std");

pub const Inputs = struct {
    ambient_vapor_pressure_kpa: f64,
    air_temperature_k: f64,
    vapor_pressure_to_volume_temperature_coefficient_k: f64 = 2.173e-3,
};

/// Exact WTHR conversion VPA=VPK*2.173E-03/TKA from wthr.f:509.
pub fn calculate(inputs: Inputs) !f64 {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteAtmosphericVaporCarrierInput;
    }
    if (inputs.ambient_vapor_pressure_kpa < 0 or
        inputs.air_temperature_k <= 0 or
        inputs.vapor_pressure_to_volume_temperature_coefficient_k <= 0)
        return error.InvalidAtmosphericVaporCarrierInput;
    const volume_fraction =
        inputs.ambient_vapor_pressure_kpa *
        inputs.vapor_pressure_to_volume_temperature_coefficient_k /
        inputs.air_temperature_k;
    if (!std.math.isFinite(volume_fraction))
        return error.AtmosphericVaporCarrierOverflow;
    return volume_fraction;
}

test "WTHR pressure temperature conversion is exact" {
    const result = try calculate(.{
        .ambient_vapor_pressure_kpa = 1.2,
        .air_temperature_k = 300,
    });
    try std.testing.expectEqual(
        @as(f64, 1.2) * 2.173e-3 / 300,
        result,
    );
}

test "zero ambient vapor pressure produces zero volume fraction" {
    try std.testing.expectEqual(
        @as(f64, 0),
        try calculate(.{
            .ambient_vapor_pressure_kpa = 0,
            .air_temperature_k = 273.15,
        }),
    );
}

test "volume fraction decreases at fixed pressure as temperature rises" {
    const cold = try calculate(.{
        .ambient_vapor_pressure_kpa = 1,
        .air_temperature_k = 270,
    });
    const warm = try calculate(.{
        .ambient_vapor_pressure_kpa = 1,
        .air_temperature_k = 300,
    });
    try std.testing.expect(cold > warm);
}

test "runtime conversion coefficient is honored" {
    const result = try calculate(.{
        .ambient_vapor_pressure_kpa = 1,
        .air_temperature_k = 250,
        .vapor_pressure_to_volume_temperature_coefficient_k = 0.002,
    });
    try std.testing.expectEqual(@as(f64, 8e-6), result);
}

test "invalid or nonfinite forcing cannot enter interpolation carrier" {
    try std.testing.expectError(
        error.InvalidAtmosphericVaporCarrierInput,
        calculate(.{
            .ambient_vapor_pressure_kpa = 1,
            .air_temperature_k = 0,
        }),
    );
    try std.testing.expectError(
        error.NonFiniteAtmosphericVaporCarrierInput,
        calculate(.{
            .ambient_vapor_pressure_kpa = std.math.nan(f64),
            .air_temperature_k = 280,
        }),
    );
}
