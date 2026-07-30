const std = @import("std");

/// Per-grid-cell values reset by HOUR1 before hourly process accumulation.
/// Mass totals use legacy grid-cell units (g element), water uses m3, canopy
/// temperatures use K, vapor pressure uses kPa, and hourly exchanges retain
/// their legacy per-timestep units.
pub const GridCellAccumulators = struct {
    soil_co2_carbon_g: f64, // UCO2S
    total_organic_matter_g: f64, // TOMT
    total_organic_nitrogen_g: f64, // TONT
    total_organic_phosphorus_g: f64, // TOPT
    water_volume_m3: f64, // UVOLW
    litter_carbon_g: f64, // URSDC
    humus_carbon_g: f64, // UORGC
    litter_nitrogen_g: f64, // URSDN
    humus_nitrogen_g: f64, // UORGN
    litter_phosphorus_g: f64, // URSDP
    humus_phosphorus_g: f64, // UORGP
    ammonium_nitrogen_g: f64, // UNH4
    nitrate_nitrogen_g: f64, // UNO3
    phosphate_phosphorus_g: f64, // UPO4
    exchangeable_phosphorus_g: f64, // UPX4
    precipitated_phosphorus_g: f64, // UPP4
    ionic_charge_equivalents: f64, // UION
    canopy_surface_temperature_k: f64, // TKCT
    canopy_air_temperature_k: f64, // TKQT
    canopy_vapor_pressure_kpa: f64, // VPQT
    canopy_evapotranspiration_m3_timestep: f64, // TEVAPP
    canopy_co2_exchange_g_c_timestep: f64, // XCNET
    canopy_ch4_exchange_g_c_timestep: f64, // XHNET
    canopy_o2_exchange_g_o_timestep: f64, // XONET
    canopy_combustion_g_c_timestep: f64, // RCGCK
};

/// Translates HOUR1 lines 2408-2432 in source assignment order.
pub fn reset(accumulators: *GridCellAccumulators) void {
    accumulators.soil_co2_carbon_g = 0.0;
    accumulators.total_organic_matter_g = 0.0;
    accumulators.total_organic_nitrogen_g = 0.0;
    accumulators.total_organic_phosphorus_g = 0.0;
    accumulators.water_volume_m3 = 0.0;
    accumulators.litter_carbon_g = 0.0;
    accumulators.humus_carbon_g = 0.0;
    accumulators.litter_nitrogen_g = 0.0;
    accumulators.humus_nitrogen_g = 0.0;
    accumulators.litter_phosphorus_g = 0.0;
    accumulators.humus_phosphorus_g = 0.0;
    accumulators.ammonium_nitrogen_g = 0.0;
    accumulators.nitrate_nitrogen_g = 0.0;
    accumulators.phosphate_phosphorus_g = 0.0;
    accumulators.exchangeable_phosphorus_g = 0.0;
    accumulators.precipitated_phosphorus_g = 0.0;
    accumulators.ionic_charge_equivalents = 0.0;
    accumulators.canopy_surface_temperature_k = 0.0;
    accumulators.canopy_air_temperature_k = 0.0;
    accumulators.canopy_vapor_pressure_kpa = 0.0;
    accumulators.canopy_evapotranspiration_m3_timestep = 0.0;
    accumulators.canopy_co2_exchange_g_c_timestep = 0.0;
    accumulators.canopy_ch4_exchange_g_c_timestep = 0.0;
    accumulators.canopy_o2_exchange_g_o_timestep = 0.0;
    accumulators.canopy_combustion_g_c_timestep = 0.0;
}

test "reset clears every hourly grid-cell accumulator" {
    var accumulators: GridCellAccumulators = undefined;
    inline for (std.meta.fields(GridCellAccumulators), 0..) |field, index| {
        @field(accumulators, field.name) = @floatFromInt(index + 1);
    }

    reset(&accumulators);

    inline for (std.meta.fields(GridCellAccumulators)) |field| {
        try std.testing.expectEqual(@as(f64, 0.0), @field(accumulators, field.name));
    }
}
