const std = @import("std");

pub const LayerPosition = enum { surface, subsurface };
pub const SaltMode = enum { static_equilibrium, dynamic_transport };

pub const GasProcessFluxes = struct {
    inhibition_factor: f64, // FINH
    soil_carbon_dioxide_g_timestep: f64, // TCO2S
    plant_carbon_dioxide_g_timestep: f64, // TCO2P
    carbon_dioxide_flux_g_timestep: f64, // TCOFLA
    methane_flux_g_timestep: f64, // TCHFLA, assigned at 3392 and 3402
    litter_carbon_dioxide_g_timestep: f64, // TLCO2P
    oxygen_flux_g_timestep: f64, // TOXFLA
    nitrous_oxide_flux_g_timestep: f64, // TN2FLA
    ammonia_flux_g_timestep: f64, // TNHFLA
    hydrogen_flux_g_timestep: f64, // THGFLA
    plant_oxygen_g_timestep: f64, // TLOXYP
    plant_methane_g_timestep: f64, // TLCH4P
    plant_nitrous_oxide_g_timestep: f64, // TLN2OP
    plant_ammonia_g_timestep: f64, // TLNH3P
    plant_hydrogen_g_timestep: f64, // TLH2GP
};

pub const PlantUptakeFluxes = struct {
    oxygen_g_timestep: f64, // TUPOXP
    soil_oxygen_g_timestep: f64, // TUPOXS
    methane_g_timestep: f64, // TUPCHS
    nitrogen_g_timestep: f64, // TUPN2S
    ammonia_non_band_g_n_timestep: f64, // TUPN3S
    ammonia_band_g_n_timestep: f64, // TUPN3B
    hydrogen_g_timestep: f64, // TUPHGS
    ammonium_non_band_g_n_timestep: f64, // TUPNH4
    nitrate_non_band_g_n_timestep: f64, // TUPNO3
    hydrogen_phosphate_non_band_g_p_timestep: f64, // TUPH2P
    dihydrogen_phosphate_non_band_g_p_timestep: f64, // TUPH1P
    ammonium_band_g_n_timestep: f64, // TUPNHB
    nitrate_band_g_n_timestep: f64, // TUPNOB
    hydrogen_phosphate_band_g_p_timestep: f64, // TUPH2B
    dihydrogen_phosphate_band_g_p_timestep: f64, // TUPH1B
    nitrogen_fixation_g_n_timestep: f64, // TUPNF
};

pub const SaltUptakeFluxes = struct {
    aluminum_g_timestep: f64, // TUPZAL
    iron_g_timestep: f64, // TUPZFE
    calcium_g_timestep: f64, // TUPZCA
    magnesium_g_timestep: f64, // TUPZMG
    sodium_g_timestep: f64, // TUPZNA
    potassium_g_timestep: f64, // TUPZKA
    sulfate_g_timestep: f64, // TUPZSO
    chloride_g_timestep: f64, // TUPZCL
    hydrogen_ion_transfer_g_timestep: f64, // TRHYSI
};

pub const LayerFluxes = struct {
    gas_processes: GasProcessFluxes,
    plant_uptake: PlantUptakeFluxes,
    salt_uptake: SaltUptakeFluxes,
};

/// Translates HOUR1 lines 3387-3430 in source order.
pub fn reset(
    layer_position: LayerPosition,
    salt_mode: SaltMode,
    fluxes: *LayerFluxes,
) void {
    if (layer_position == .surface) return;

    fluxes.gas_processes.inhibition_factor = 0.0;
    fluxes.gas_processes.soil_carbon_dioxide_g_timestep = 0.0;
    fluxes.gas_processes.plant_carbon_dioxide_g_timestep = 0.0;
    fluxes.gas_processes.carbon_dioxide_flux_g_timestep = 0.0;
    fluxes.gas_processes.methane_flux_g_timestep = 0.0;
    fluxes.gas_processes.litter_carbon_dioxide_g_timestep = 0.0;
    fluxes.plant_uptake.oxygen_g_timestep = 0.0;
    fluxes.plant_uptake.soil_oxygen_g_timestep = 0.0;
    fluxes.plant_uptake.methane_g_timestep = 0.0;
    fluxes.plant_uptake.nitrogen_g_timestep = 0.0;
    fluxes.plant_uptake.ammonia_non_band_g_n_timestep = 0.0;
    fluxes.plant_uptake.ammonia_band_g_n_timestep = 0.0;
    fluxes.plant_uptake.hydrogen_g_timestep = 0.0;
    fluxes.gas_processes.oxygen_flux_g_timestep = 0.0;
    fluxes.gas_processes.methane_flux_g_timestep = 0.0;
    fluxes.gas_processes.nitrous_oxide_flux_g_timestep = 0.0;
    fluxes.gas_processes.ammonia_flux_g_timestep = 0.0;
    fluxes.gas_processes.hydrogen_flux_g_timestep = 0.0;
    fluxes.gas_processes.plant_oxygen_g_timestep = 0.0;
    fluxes.gas_processes.plant_methane_g_timestep = 0.0;
    fluxes.gas_processes.plant_nitrous_oxide_g_timestep = 0.0;
    fluxes.gas_processes.plant_ammonia_g_timestep = 0.0;
    fluxes.gas_processes.plant_hydrogen_g_timestep = 0.0;
    inline for (std.meta.fields(PlantUptakeFluxes)[7..]) |field| {
        @field(fluxes.plant_uptake, field.name) = 0.0;
    }
    if (salt_mode == .dynamic_transport) {
        inline for (std.meta.fields(SaltUptakeFluxes)) |field| {
            @field(fluxes.salt_uptake, field.name) = 0.0;
        }
    }
}

test "subsurface dynamic reset clears all gas uptake and salt fluxes" {
    var fluxes: LayerFluxes = undefined;
    inline for (std.meta.fields(LayerFluxes)) |group| {
        inline for (std.meta.fields(group.type)) |field| {
            @field(@field(fluxes, group.name), field.name) = 3.0;
        }
    }
    reset(.subsurface, .dynamic_transport, &fluxes);
    inline for (std.meta.fields(LayerFluxes)) |group| {
        inline for (std.meta.fields(group.type)) |field| {
            try std.testing.expectEqual(
                @as(f64, 0.0),
                @field(@field(fluxes, group.name), field.name),
            );
        }
    }
}

test "surface guard preserves state and static mode preserves salt uptake" {
    var fluxes: LayerFluxes = undefined;
    fluxes.gas_processes.inhibition_factor = 5.0;
    reset(.surface, .dynamic_transport, &fluxes);
    try std.testing.expectEqual(@as(f64, 5.0), fluxes.gas_processes.inhibition_factor);

    fluxes.plant_uptake.oxygen_g_timestep = 6.0;
    fluxes.salt_uptake.aluminum_g_timestep = 7.0;
    reset(.subsurface, .static_equilibrium, &fluxes);
    try std.testing.expectEqual(@as(f64, 0.0), fluxes.plant_uptake.oxygen_g_timestep);
    try std.testing.expectEqual(@as(f64, 7.0), fluxes.salt_uptake.aluminum_g_timestep);
}
