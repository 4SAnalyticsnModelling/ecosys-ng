const std = @import("std");

pub const SaltMode = enum { static_equilibrium, dynamic_transport };

pub const WaterHeatFluxes = struct {
    soil_water_m3_timestep: f64, // FLSW
    soil_vapor_m3_timestep: f64, // FLSV
    horizontal_water_m3_timestep: f64, // FLSWH
    soil_heat_j_timestep: f64, // HFLSW
    residue_water_m3_timestep: f64, // FLSWR
    residue_vapor_m3_timestep: f64, // FLSVR
    residue_heat_j_timestep: f64, // HFLSWR
    snow_water_m3_timestep: f64, // XFLWS
    liquid_water_m3_timestep: f64, // XFLWW
    vapor_water_m3_timestep: f64, // XFLWV
    ice_water_m3_timestep: f64, // XFLWI
    water_vapor_flux_m3_timestep: f64, // XWFLVW
    soil_vapor_flux_m3_timestep: f64, // XWFLVS
    vapor_heat_surface_j_timestep: f64, // XHFLV0
    vapor_heat_water_j_timestep: f64, // XHFLWW
    freezing_soil_water_m3_timestep: f64, // XWFLFS
    freezing_ice_water_m3_timestep: f64, // XWFLFI
    freezing_heat_j_timestep: f64, // XHFLF0
};

pub const GasNutrientFluxes = struct {
    carbon_dioxide_g_timestep: f64, // XCOBLS
    methane_g_timestep: f64, // XCHBLS
    oxygen_g_timestep: f64, // XOXBLS
    nitrogen_g_timestep: f64, // XNGBLS
    nitrous_oxide_g_timestep: f64, // XN2BLS
    ammonium_g_n_timestep: f64, // XN4BLW
    ammonia_g_n_timestep: f64, // XN3BLW
    nitrate_g_n_timestep: f64, // XNOBLW
    phosphate_h2po4_g_p_timestep: f64, // XH1PBS
    phosphate_hpo4_g_p_timestep: f64, // XH2PBS
};

/// Source order of dynamic salt species in HOUR1 lines 3110-3150.
pub const SaltSpecies = enum {
    aluminum,
    iron,
    hydrogen,
    calcium,
    magnesium,
    sodium,
    potassium,
    hydroxide,
    sulfate,
    chloride,
    carbonate,
    bicarbonate,
    aluminum_hydroxide_1,
    aluminum_hydroxide_2,
    aluminum_hydroxide_3,
    aluminum_hydroxide_4,
    aluminum_sulfate,
    iron_hydroxide_1,
    iron_hydroxide_2,
    iron_hydroxide_3,
    iron_hydroxide_4,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
    hydrogen_phosphate_0,
    hydrogen_phosphate_3,
    iron_phosphate_1,
    iron_phosphate_2,
    calcium_phosphate_0,
    calcium_phosphate_1,
    calcium_phosphate_2,
    magnesium_phosphate_1,
};

pub const LayerBoundaryFluxes = struct {
    water_heat: WaterHeatFluxes,
    gas_nutrients: GasNutrientFluxes,
    salt_species_mass_timestep: []f64,
};

pub const ResetError = error{SaltSpeciesCountMismatch};

/// Translates HOUR1 lines 3080-3152 over runtime-allocated soil layers.
pub fn reset(
    salt_mode: SaltMode,
    layers: []LayerBoundaryFluxes,
) ResetError!void {
    if (salt_mode == .dynamic_transport) {
        const salt_count = std.meta.fields(SaltSpecies).len;
        for (layers) |layer| {
            if (layer.salt_species_mass_timestep.len != salt_count) {
                return error.SaltSpeciesCountMismatch;
            }
        }
    }
    for (layers) |*layer| {
        inline for (std.meta.fields(WaterHeatFluxes)) |field| {
            @field(layer.water_heat, field.name) = 0.0;
        }
        inline for (std.meta.fields(GasNutrientFluxes)) |field| {
            @field(layer.gas_nutrients, field.name) = 0.0;
        }
        if (salt_mode == .dynamic_transport) {
            for (layer.salt_species_mass_timestep) |*flux| flux.* = 0.0;
        }
    }
}

test "runtime layers reset base fluxes and dynamic salt species" {
    const allocator = std.testing.allocator;
    const layers = try allocator.alloc(LayerBoundaryFluxes, 6);
    defer allocator.free(layers);
    for (layers) |*layer| {
        layer.salt_species_mass_timestep =
            try allocator.alloc(f64, std.meta.fields(SaltSpecies).len);
        @memset(layer.salt_species_mass_timestep, 4.0);
        inline for (std.meta.fields(WaterHeatFluxes)) |field| {
            @field(layer.water_heat, field.name) = 2.0;
        }
        inline for (std.meta.fields(GasNutrientFluxes)) |field| {
            @field(layer.gas_nutrients, field.name) = 3.0;
        }
    }
    defer for (layers) |layer| allocator.free(layer.salt_species_mass_timestep);

    try reset(.dynamic_transport, layers);

    for (layers) |layer| {
        inline for (std.meta.fields(WaterHeatFluxes)) |field| {
            try std.testing.expectEqual(@as(f64, 0.0), @field(layer.water_heat, field.name));
        }
        inline for (std.meta.fields(GasNutrientFluxes)) |field| {
            try std.testing.expectEqual(@as(f64, 0.0), @field(layer.gas_nutrients, field.name));
        }
        for (layer.salt_species_mass_timestep) |flux| {
            try std.testing.expectEqual(@as(f64, 0.0), flux);
        }
    }
}

test "static salt mode resets base fluxes but preserves salt state" {
    var salt = [_]f64{7.0};
    var layers = [_]LayerBoundaryFluxes{.{
        .water_heat = std.mem.zeroes(WaterHeatFluxes),
        .gas_nutrients = std.mem.zeroes(GasNutrientFluxes),
        .salt_species_mass_timestep = &salt,
    }};
    layers[0].water_heat.soil_water_m3_timestep = 5.0;
    try reset(.static_equilibrium, &layers);
    try std.testing.expectEqual(@as(f64, 0.0), layers[0].water_heat.soil_water_m3_timestep);
    try std.testing.expectEqual(@as(f64, 7.0), salt[0]);
}
