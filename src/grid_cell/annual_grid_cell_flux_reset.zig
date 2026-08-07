const std = @import("std");
const execution_calendar_date = @import("../driver/execution_calendar_date.zig");

/// One cell's annual diagnostic accumulators from DAY. All fields are
/// extensive or area-normalized only at output; unit suffixes identify their
/// elemental carrier.
pub const State = struct {
    net_carbon_g_c: f64,
    net_methane_carbon_g_c: f64,
    net_oxygen_g_o: f64,
    organic_fertilizer_carbon_g_c: f64,
    plant_carbon_sink_g_c: f64,
    crop_carbon_output_g_c: f64,
    dissolved_organic_carbon_runoff_g_c: f64,
    dissolved_organic_carbon_drainage_g_c: f64,
    dissolved_inorganic_carbon_runoff_g_c: f64,
    dissolved_inorganic_carbon_drainage_g_c: f64,
    biome_carbon_input_g_c: f64,
    precipitation_water_m3: f64,
    evaporation_water_m3: f64,
    runoff_water_m3: f64,
    sediment_output_g: f64,
    boundary_water_output_m3: f64,
    drainage_water_m3: f64,
    ion_output_mol: f64,
    fertilizer_nitrogen_g_n: f64,
    plant_nitrogen_sink_g_n: f64,
    dissolved_organic_nitrogen_runoff_g_n: f64,
    dissolved_organic_nitrogen_drainage_g_n: f64,
    dissolved_inorganic_nitrogen_runoff_g_n: f64,
    dissolved_inorganic_nitrogen_drainage_g_n: f64,
    fertilizer_phosphorus_g_p: f64,
    plant_phosphorus_sink_g_p: f64,
    dissolved_organic_phosphorus_runoff_g_p: f64,
    dissolved_organic_phosphorus_drainage_g_p: f64,
    dissolved_inorganic_phosphorus_runoff_g_p: f64,
    dissolved_inorganic_phosphorus_drainage_g_p: f64,
    carbon_dioxide_exchange_g_c: f64,
    methane_exchange_g_c: f64,
    oxygen_exchange_g_o: f64,
    nitrous_oxide_exchange_g_n: f64,
    dinitrogen_exchange_g_n: f64,
    ammonia_exchange_g_n: f64,
    soil_dinitrogen_fixation_g_n: f64,
    hydrogen_exchange_g_h: f64,
    ecosystem_respiration_g_c: f64,
    gross_primary_productivity_g_c: f64,
    net_primary_productivity_g_c: f64,
    autotrophic_respiration_g_c: f64,
    canopy_carbon_g_c: f64,
    canopy_co2_exchange_g_c: f64,
    harvest_carbon_g_c: f64,
    harvest_nitrogen_g_n: f64,
    harvest_phosphorus_g_p: f64,
    root_ammonium_uptake_g_n: f64,
    root_phosphate_uptake_g_p: f64,
    fire_carbon_dioxide_g_c: f64,
    fire_methane_g_c: f64,
    fire_nitrogen_oxide_g_n: f64,
    fire_phosphorus_oxide_g_p: f64,
    residue_fire_carbon_dioxide_g_c: f64,
    residue_fire_methane_g_c: f64,
    residue_fire_nitrogen_oxide_g_n: f64,
    residue_fire_phosphorus_oxide_g_p: f64,
    soil_fire_carbon_loss_g_c: f64,
    soil_fire_nitrogen_loss_g_n: f64,
    soil_fire_phosphorus_loss_g_p: f64,
    erosion_sediment_g: f64,
};

/// Exact annual reset condition is day one for either latitude hemisphere.
/// Erosion sediment resets only for erosion modes 1 or 3.
pub fn resetIfYearStart(
    state: *State,
    execution_day_of_year: u16,
    execution_year: u16,
    latitude_degrees_north: f64,
    erosion_mode: u8,
) !bool {
    _ = execution_calendar_date.fromDayOfYear(
        execution_day_of_year,
        execution_year,
    ) catch return error.InvalidAnnualCellResetDate;
    if (!std.math.isFinite(latitude_degrees_north) or
        latitude_degrees_north < -90 or latitude_degrees_north > 90)
        return error.InvalidAnnualCellResetLatitude;
    if (erosion_mode > 3) return error.InvalidAnnualCellResetErosionMode;
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (!std.math.isFinite(@field(state.*, field.name)))
            return error.NonFiniteAnnualCellAccumulator;
    if (execution_day_of_year != 1) return false;

    const retained_sediment =
        if (erosion_mode == 1 or erosion_mode == 3)
            @as(f64, 0)
        else
            state.erosion_sediment_g;
    inline for (@typeInfo(State).@"struct".fields) |field|
        @field(state.*, field.name) = 0;
    state.erosion_sediment_g = retained_sediment;
    return true;
}

fn filled(value: f64) State {
    var state: State = undefined;
    inline for (@typeInfo(State).@"struct".fields) |field|
        @field(state, field.name) = value;
    return state;
}

test "DAY year start clears every typed cell accumulator" {
    var state = filled(7);
    try std.testing.expect(try resetIfYearStart(&state, 1, 2000, 53.5, 1));
    inline for (@typeInfo(State).@"struct".fields) |field|
        try std.testing.expectEqual(
            @as(f64, 0),
            @field(state, field.name),
        );
}

test "annual reset is identical in northern and southern hemispheres" {
    var north = filled(1);
    var south = filled(1);
    _ = try resetIfYearStart(&north, 1, 2000, 80, 3);
    _ = try resetIfYearStart(&south, 1, 2000, -80, 3);
    try std.testing.expectEqualDeep(north, south);
}

test "erosion sediment follows exact mode one or three condition" {
    inline for (.{ @as(u8, 0), @as(u8, 2) }) |mode| {
        var state = filled(4);
        _ = try resetIfYearStart(&state, 1, 2000, 0, mode);
        try std.testing.expectEqual(@as(f64, 4), state.erosion_sediment_g);
        try std.testing.expectEqual(@as(f64, 0), state.net_carbon_g_c);
    }
    inline for (.{ @as(u8, 1), @as(u8, 3) }) |mode| {
        var state = filled(4);
        _ = try resetIfYearStart(&state, 1, 2000, 0, mode);
        try std.testing.expectEqual(@as(f64, 0), state.erosion_sediment_g);
    }
}

test "ordinary day and invalid late accumulator preserve state" {
    var state = filled(3);
    const before = state;
    try std.testing.expect(!try resetIfYearStart(&state, 2, 2000, 0, 1));
    try std.testing.expectEqualDeep(before, state);

    state.fire_phosphorus_oxide_g_p = std.math.nan(f64);
    const carbon_before = state.net_carbon_g_c;
    try std.testing.expectError(
        error.NonFiniteAnnualCellAccumulator,
        resetIfYearStart(&state, 1, 2000, 0, 1),
    );
    try std.testing.expectEqual(carbon_before, state.net_carbon_g_c);
    try std.testing.expect(std.math.isNan(state.fire_phosphorus_oxide_g_p));
}

test "annual reset dates preserve DAY modulo-four chronology" {
    var state = filled(1);
    try std.testing.expect(!try resetIfYearStart(&state, 366, 1900, 0, 0));
    try std.testing.expectError(
        error.InvalidAnnualCellResetDate,
        resetIfYearStart(&state, 366, 1901, 0, 0),
    );
    try std.testing.expectError(
        error.InvalidAnnualCellResetDate,
        resetIfYearStart(&state, 1, 0, 0, 0),
    );
}
