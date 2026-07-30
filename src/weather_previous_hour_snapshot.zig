const std = @import("std");

pub const CurrentForcing = struct {
    shortwave_radiation_mj_per_m2_h: f64,
    air_temperature_k: f64,
    atmospheric_vapor_concentration_m3_per_m3: f64,
};

pub const Snapshot = struct {
    previous_shortwave_radiation_mj_per_m2_h: f64,
    previous_air_temperature_k: f64,
    previous_atmospheric_vapor_concentration_m3_per_m3: f64,
};

pub const Inputs = struct {
    execution_day: u32,
    first_execution_day: u32,
    source_hour: u8,
    annual_mean_air_temperature_k: f64,
    current: CurrentForcing,
};

/// Exact WTHR previous-hour carrier update from wthr.f:69-83 and 167-181.
///
/// Both daily and hourly weather paths execute the same source branch. At
/// hour one of the first execution day, radiation and vapor are reset to zero
/// and temperature is seeded from annual mean air temperature. Every other
/// hour snapshots the current forcing before the next forcing is calculated.
pub fn capture(inputs: Inputs) !Snapshot {
    if (inputs.execution_day == 0 or inputs.first_execution_day == 0 or
        inputs.source_hour < 1 or inputs.source_hour > 24)
        return error.InvalidWeatherSnapshotTime;
    inline for (.{
        inputs.annual_mean_air_temperature_k,
        inputs.current.shortwave_radiation_mj_per_m2_h,
        inputs.current.air_temperature_k,
        inputs.current.atmospheric_vapor_concentration_m3_per_m3,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteWeatherSnapshotInput;
    }
    if (inputs.annual_mean_air_temperature_k <= 0 or
        inputs.current.air_temperature_k <= 0 or
        inputs.current.shortwave_radiation_mj_per_m2_h < 0 or
        inputs.current.atmospheric_vapor_concentration_m3_per_m3 < 0)
        return error.InvalidWeatherSnapshotInput;

    if (inputs.execution_day == inputs.first_execution_day and
        inputs.source_hour == 1)
    {
        return .{
            .previous_shortwave_radiation_mj_per_m2_h = 0,
            .previous_air_temperature_k = inputs.annual_mean_air_temperature_k,
            .previous_atmospheric_vapor_concentration_m3_per_m3 = 0,
        };
    }
    return .{
        .previous_shortwave_radiation_mj_per_m2_h = inputs.current.shortwave_radiation_mj_per_m2_h,
        .previous_air_temperature_k = inputs.current.air_temperature_k,
        .previous_atmospheric_vapor_concentration_m3_per_m3 = inputs.current.atmospheric_vapor_concentration_m3_per_m3,
    };
}

fn exampleInputs() Inputs {
    return .{
        .execution_day = 100,
        .first_execution_day = 100,
        .source_hour = 1,
        .annual_mean_air_temperature_k = 281,
        .current = .{
            .shortwave_radiation_mj_per_m2_h = 0.5,
            .air_temperature_k = 290,
            .atmospheric_vapor_concentration_m3_per_m3 = 0.001,
        },
    };
}

test "first execution hour uses annual temperature and zero carriers" {
    const snapshot = try capture(exampleInputs());
    try std.testing.expectEqual(
        @as(f64, 0),
        snapshot.previous_shortwave_radiation_mj_per_m2_h,
    );
    try std.testing.expectEqual(
        @as(f64, 281),
        snapshot.previous_air_temperature_k,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        snapshot.previous_atmospheric_vapor_concentration_m3_per_m3,
    );
}

test "later hours snapshot every current forcing carrier" {
    var inputs = exampleInputs();
    inputs.source_hour = 2;
    const snapshot = try capture(inputs);
    try std.testing.expectEqual(
        inputs.current.shortwave_radiation_mj_per_m2_h,
        snapshot.previous_shortwave_radiation_mj_per_m2_h,
    );
    try std.testing.expectEqual(
        inputs.current.air_temperature_k,
        snapshot.previous_air_temperature_k,
    );
    try std.testing.expectEqual(
        inputs.current.atmospheric_vapor_concentration_m3_per_m3,
        snapshot.previous_atmospheric_vapor_concentration_m3_per_m3,
    );
}

test "hour one after first execution day snapshots current forcing" {
    var inputs = exampleInputs();
    inputs.execution_day = 101;
    const snapshot = try capture(inputs);
    try std.testing.expectEqual(
        inputs.current.air_temperature_k,
        snapshot.previous_air_temperature_k,
    );
}

test "invalid late carrier fails before producing interpolation state" {
    var inputs = exampleInputs();
    inputs.source_hour = 2;
    inputs.current.atmospheric_vapor_concentration_m3_per_m3 =
        std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteWeatherSnapshotInput,
        capture(inputs),
    );
}

test "source hours are strictly one through twenty four" {
    var inputs = exampleInputs();
    inputs.source_hour = 0;
    try std.testing.expectError(
        error.InvalidWeatherSnapshotTime,
        capture(inputs),
    );
    inputs.source_hour = 25;
    try std.testing.expectError(
        error.InvalidWeatherSnapshotTime,
        capture(inputs),
    );
}
