const std = @import("std");

pub const State = struct {
    solar_radiation_megajoules_per_m2: f64,
    maximum_air_temperature_c: f64,
    minimum_air_temperature_c: f64,
    maximum_vapor_pressure_kpa: f64,
    minimum_vapor_pressure_kpa: f64,
    wind_travel_m: f64,
    total_water_input_mm: f64,
    maximum_soil_temperature_c_by_layer: []f64,
    minimum_soil_temperature_c_by_layer: []f64,
};

/// Exact DAY reset for daily weather diagnostics and layer extrema.
pub fn reset(state: *State) !void {
    if (state.maximum_soil_temperature_c_by_layer.len == 0 or
        state.maximum_soil_temperature_c_by_layer.len !=
            state.minimum_soil_temperature_c_by_layer.len)
        return error.DailyWeatherResetLayerDimensionMismatch;
    state.solar_radiation_megajoules_per_m2 = 0;
    state.maximum_air_temperature_c = -100;
    state.minimum_air_temperature_c = 100;
    state.maximum_vapor_pressure_kpa = 0;
    state.minimum_vapor_pressure_kpa = 100;
    state.wind_travel_m = 0;
    state.total_water_input_mm = 0;
    @memset(state.maximum_soil_temperature_c_by_layer, -9999);
    @memset(state.minimum_soil_temperature_c_by_layer, 9999);
}

test "DAY resets weather diagnostics to exact source sentinels" {
    var maximum = [_]f64{ 1, 2, 3 };
    var minimum = [_]f64{ -1, -2, -3 };
    var state: State = .{
        .solar_radiation_megajoules_per_m2 = 10,
        .maximum_air_temperature_c = 30,
        .minimum_air_temperature_c = 10,
        .maximum_vapor_pressure_kpa = 2,
        .minimum_vapor_pressure_kpa = 1,
        .wind_travel_m = 500,
        .total_water_input_mm = 20,
        .maximum_soil_temperature_c_by_layer = &maximum,
        .minimum_soil_temperature_c_by_layer = &minimum,
    };
    try reset(&state);
    try std.testing.expectEqual(@as(f64, 0), state.solar_radiation_megajoules_per_m2);
    try std.testing.expectEqual(@as(f64, -100), state.maximum_air_temperature_c);
    try std.testing.expectEqual(@as(f64, 100), state.minimum_air_temperature_c);
    try std.testing.expectEqual(@as(f64, 0), state.maximum_vapor_pressure_kpa);
    try std.testing.expectEqual(@as(f64, 100), state.minimum_vapor_pressure_kpa);
    try std.testing.expectEqual(@as(f64, 0), state.wind_travel_m);
    try std.testing.expectEqual(@as(f64, 0), state.total_water_input_mm);
    for (maximum) |value|
        try std.testing.expectEqual(@as(f64, -9999), value);
    for (minimum) |value|
        try std.testing.expectEqual(@as(f64, 9999), value);
}

test "runtime layer count has no source JZ ceiling" {
    const allocator = std.testing.allocator;
    const count = 257;
    const maximum = try allocator.alloc(f64, count);
    defer allocator.free(maximum);
    const minimum = try allocator.alloc(f64, count);
    defer allocator.free(minimum);
    @memset(maximum, 0);
    @memset(minimum, 0);
    var state: State = .{
        .solar_radiation_megajoules_per_m2 = 0,
        .maximum_air_temperature_c = 0,
        .minimum_air_temperature_c = 0,
        .maximum_vapor_pressure_kpa = 0,
        .minimum_vapor_pressure_kpa = 0,
        .wind_travel_m = 0,
        .total_water_input_mm = 0,
        .maximum_soil_temperature_c_by_layer = maximum,
        .minimum_soil_temperature_c_by_layer = minimum,
    };
    try reset(&state);
    try std.testing.expectEqual(@as(f64, -9999), maximum[count - 1]);
    try std.testing.expectEqual(@as(f64, 9999), minimum[count - 1]);
}

test "dimension failure leaves scalar and layer state unchanged" {
    var maximum = [_]f64{ 1, 2 };
    var minimum = [_]f64{3};
    const maximum_before = maximum;
    const minimum_before = minimum;
    var state: State = .{
        .solar_radiation_megajoules_per_m2 = 10,
        .maximum_air_temperature_c = 30,
        .minimum_air_temperature_c = 10,
        .maximum_vapor_pressure_kpa = 2,
        .minimum_vapor_pressure_kpa = 1,
        .wind_travel_m = 500,
        .total_water_input_mm = 20,
        .maximum_soil_temperature_c_by_layer = &maximum,
        .minimum_soil_temperature_c_by_layer = &minimum,
    };
    try std.testing.expectError(
        error.DailyWeatherResetLayerDimensionMismatch,
        reset(&state),
    );
    try std.testing.expectEqual(@as(f64, 10), state.solar_radiation_megajoules_per_m2);
    try std.testing.expectEqualSlices(f64, &maximum_before, &maximum);
    try std.testing.expectEqualSlices(f64, &minimum_before, &minimum);
}
