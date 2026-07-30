const std = @import("std");

pub const State = struct {
    temperature_k: *f64,
    temperature_c: *f64,
    water_vapor_volume_m3: *f64,
};

/// Exact translation of legacy `STARTS` lines 1222--1224.
pub fn initialize(
    state: State,
    atmospheric_temperature_k: f64,
    atmospheric_temperature_c: f64,
) !void {
    if (!std.math.isFinite(atmospheric_temperature_k) or
        !std.math.isFinite(atmospheric_temperature_c))
        return error.NonFiniteLayerThermalInput;
    if (atmospheric_temperature_k <= 0.0)
        return error.InvalidLayerThermalInput;

    state.temperature_k.* = atmospheric_temperature_k;
    state.temperature_c.* = atmospheric_temperature_c;
    state.water_vapor_volume_m3.* = 0.0;
}

test "STARTS copies atmospheric layer temperatures and zeros vapor" {
    var temperature_k: f64 = 0.0;
    var temperature_c: f64 = 0.0;
    var vapor_m3: f64 = 9.0;
    try initialize(.{
        .temperature_k = &temperature_k,
        .temperature_c = &temperature_c,
        .water_vapor_volume_m3 = &vapor_m3,
    }, 283.15, 10.0);
    try std.testing.expectEqual(@as(f64, 283.15), temperature_k);
    try std.testing.expectEqual(@as(f64, 10.0), temperature_c);
    try std.testing.expectEqual(@as(f64, 0.0), vapor_m3);
}

test "non-finite temperature fails before mutation" {
    var temperature_k: f64 = 1.0;
    var temperature_c: f64 = 2.0;
    var vapor_m3: f64 = 3.0;
    try std.testing.expectError(
        error.NonFiniteLayerThermalInput,
        initialize(.{
            .temperature_k = &temperature_k,
            .temperature_c = &temperature_c,
            .water_vapor_volume_m3 = &vapor_m3,
        }, std.math.nan(f64), 10.0),
    );
    try std.testing.expectEqual(@as(f64, 1.0), temperature_k);
    try std.testing.expectEqual(@as(f64, 2.0), temperature_c);
    try std.testing.expectEqual(@as(f64, 3.0), vapor_m3);
}
