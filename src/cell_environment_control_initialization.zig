const std = @import("std");

pub const State = struct {
    minimum_surface_elevation_m: *f64,
    maximum_soil_surface_depth_m: *f64,
    maximum_temperature_offset_c: []f64,
    minimum_temperature_offset_c: []f64,
    radiation_multiplier: []f64,
    wind_speed_multiplier: []f64,
    vapor_pressure_multiplier: []f64,
    precipitation_multiplier: []f64,
    irrigation_multiplier: []f64,
    carbon_dioxide_multiplier: []f64,
    precipitation_ammonium_multiplier: []f64,
    precipitation_nitrate_multiplier: []f64,
    urea_formulation_index: []u8,
    ammonium_band_active: []u8,
    nitrate_band_active: []u8,
    phosphate_band_active: []u8,
    soil_state_active: []u8,
    active_plant_flag_count: []usize,
};

/// Exact source-order translation of legacy `STARTS` lines 368--389.
///
/// Climate fields use cell-major `[cell][modifier_period]` storage. The source
/// uses twelve monthly periods; the translated owner accepts the compulsory
/// runtime period count.
pub fn initialize(
    state: State,
    cell_count: usize,
    modifier_period_count: usize,
) !void {
    if (cell_count == 0 or modifier_period_count == 0)
        return error.InvalidEnvironmentControlDimensions;
    const modifier_count = std.math.mul(
        usize,
        cell_count,
        modifier_period_count,
    ) catch return error.DimensionOverflow;
    inline for (.{
        state.maximum_temperature_offset_c,
        state.minimum_temperature_offset_c,
        state.radiation_multiplier,
        state.wind_speed_multiplier,
        state.vapor_pressure_multiplier,
        state.precipitation_multiplier,
        state.irrigation_multiplier,
        state.carbon_dioxide_multiplier,
        state.precipitation_ammonium_multiplier,
        state.precipitation_nitrate_multiplier,
    }) |values| {
        if (values.len != modifier_count)
            return error.EnvironmentControlDimensionMismatch;
    }
    inline for (.{
        state.urea_formulation_index,
        state.ammonium_band_active,
        state.nitrate_band_active,
        state.phosphate_band_active,
        state.soil_state_active,
    }) |values| {
        if (values.len != cell_count)
            return error.EnvironmentControlDimensionMismatch;
    }
    if (state.active_plant_flag_count.len != cell_count)
        return error.EnvironmentControlDimensionMismatch;

    state.minimum_surface_elevation_m.* = 0.0;
    state.maximum_soil_surface_depth_m.* = 0.0;
    for (0..cell_count) |cell| {
        for (0..modifier_period_count) |period| {
            const index = cell * modifier_period_count + period;
            state.maximum_temperature_offset_c[index] = 0.0;
            state.minimum_temperature_offset_c[index] = 0.0;
            state.radiation_multiplier[index] = 1.0;
            state.wind_speed_multiplier[index] = 1.0;
            state.vapor_pressure_multiplier[index] = 1.0;
            state.precipitation_multiplier[index] = 1.0;
            state.irrigation_multiplier[index] = 1.0;
            state.carbon_dioxide_multiplier[index] = 1.0;
            state.precipitation_ammonium_multiplier[index] = 1.0;
            state.precipitation_nitrate_multiplier[index] = 1.0;
        }
        state.urea_formulation_index[cell] = 0;
        state.ammonium_band_active[cell] = 0;
        state.nitrate_band_active[cell] = 0;
        state.phosphate_band_active[cell] = 0;
        state.soil_state_active[cell] = 1;
        state.active_plant_flag_count[cell] = 0;
    }
}

test "STARTS initializes source monthly controls and cell flags" {
    const cells = 2;
    const periods = 12;
    var minimum_elevation = std.math.nan(f64);
    var maximum_depth = std.math.nan(f64);
    var maximum_temperature = [_]f64{7.0} ** (cells * periods);
    var minimum_temperature = [_]f64{7.0} ** (cells * periods);
    var radiation = [_]f64{7.0} ** (cells * periods);
    var wind = [_]f64{7.0} ** (cells * periods);
    var vapor = [_]f64{7.0} ** (cells * periods);
    var precipitation = [_]f64{7.0} ** (cells * periods);
    var irrigation = [_]f64{7.0} ** (cells * periods);
    var carbon_dioxide = [_]f64{7.0} ** (cells * periods);
    var ammonium = [_]f64{7.0} ** (cells * periods);
    var nitrate = [_]f64{7.0} ** (cells * periods);
    var urea = [_]u8{9} ** cells;
    var ammonium_band = [_]u8{9} ** cells;
    var nitrate_band = [_]u8{9} ** cells;
    var phosphate_band = [_]u8{9} ** cells;
    var soil_active = [_]u8{9} ** cells;
    var plant_flags = [_]usize{9} ** cells;

    try initialize(.{
        .minimum_surface_elevation_m = &minimum_elevation,
        .maximum_soil_surface_depth_m = &maximum_depth,
        .maximum_temperature_offset_c = &maximum_temperature,
        .minimum_temperature_offset_c = &minimum_temperature,
        .radiation_multiplier = &radiation,
        .wind_speed_multiplier = &wind,
        .vapor_pressure_multiplier = &vapor,
        .precipitation_multiplier = &precipitation,
        .irrigation_multiplier = &irrigation,
        .carbon_dioxide_multiplier = &carbon_dioxide,
        .precipitation_ammonium_multiplier = &ammonium,
        .precipitation_nitrate_multiplier = &nitrate,
        .urea_formulation_index = &urea,
        .ammonium_band_active = &ammonium_band,
        .nitrate_band_active = &nitrate_band,
        .phosphate_band_active = &phosphate_band,
        .soil_state_active = &soil_active,
        .active_plant_flag_count = &plant_flags,
    }, cells, periods);

    try std.testing.expectEqual(@as(f64, 0.0), minimum_elevation);
    try std.testing.expectEqual(@as(f64, 0.0), maximum_depth);
    for (maximum_temperature) |value|
        try std.testing.expectEqual(@as(f64, 0.0), value);
    for (minimum_temperature) |value|
        try std.testing.expectEqual(@as(f64, 0.0), value);
    inline for (.{
        radiation,
        wind,
        vapor,
        precipitation,
        irrigation,
        carbon_dioxide,
        ammonium,
        nitrate,
    }) |values| {
        for (values) |value| try std.testing.expectEqual(@as(f64, 1.0), value);
    }
    try std.testing.expectEqualSlices(u8, &.{ 1, 1 }, &soil_active);
    try std.testing.expectEqualSlices(usize, &.{ 0, 0 }, &plant_flags);
}

test "dimension mismatch preserves all state including leading scalars" {
    var minimum_elevation: f64 = 8.0;
    var maximum_depth: f64 = 9.0;
    var short = [_]f64{ 1.0, 2.0, 3.0 };
    var cell_float = [_]f64{ 4.0, 5.0 };
    var flags = [_]u8{ 6, 7 };
    var plant_flags = [_]usize{ 8, 9 };

    try std.testing.expectError(
        error.EnvironmentControlDimensionMismatch,
        initialize(.{
            .minimum_surface_elevation_m = &minimum_elevation,
            .maximum_soil_surface_depth_m = &maximum_depth,
            .maximum_temperature_offset_c = &short,
            .minimum_temperature_offset_c = &cell_float,
            .radiation_multiplier = &cell_float,
            .wind_speed_multiplier = &cell_float,
            .vapor_pressure_multiplier = &cell_float,
            .precipitation_multiplier = &cell_float,
            .irrigation_multiplier = &cell_float,
            .carbon_dioxide_multiplier = &cell_float,
            .precipitation_ammonium_multiplier = &cell_float,
            .precipitation_nitrate_multiplier = &cell_float,
            .urea_formulation_index = &flags,
            .ammonium_band_active = &flags,
            .nitrate_band_active = &flags,
            .phosphate_band_active = &flags,
            .soil_state_active = &flags,
            .active_plant_flag_count = &plant_flags,
        }, 2, 2),
    );

    try std.testing.expectEqual(@as(f64, 8.0), minimum_elevation);
    try std.testing.expectEqual(@as(f64, 9.0), maximum_depth);
    try std.testing.expectEqualSlices(f64, &.{ 1, 2, 3 }, &short);
}
