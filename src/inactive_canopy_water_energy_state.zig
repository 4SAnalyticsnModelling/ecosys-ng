const std = @import("std");

pub const Inputs = struct {
    intercepted_water_volume_m3: f64,
    foliar_water_retention_m3_per_h: f64,
    water_timestep_h_per_step: f64,
    canopy_height_m: f64,
    snow_depth_m: f64,
    depth_tolerance_m: f64,
    ambient_air_temperature_k: f64,
    snow_air_temperature_k: f64,
    atmospheric_vapor_concentration_m3_per_m3: f64,
    emissivity: f64,
    stefan_boltzmann_megajoules_per_m2_h_k4: f64,
    canopy_radiation_fraction: f64,
    horizontal_cell_area_m2: f64,
    energy_timestep_h_per_step: f64,
    root_reference_soil_total_water_potential_mpa: f64,
    minimum_dry_matter_fraction_g_c_per_g: f64,
    canopy_nonstructural_carbon_g_per_g_c: f64,
    canopy_nonstructural_nitrogen_g_per_g_c: f64,
    canopy_nonstructural_phosphorus_g_per_g_c: f64,
    osmotic_potential_at_zero_total_mpa: f64,
    osmotic_temperature_k: f64,
    canopy_salt_concentration_mol_per_g_c: f64,
    turgor_response_shape_per_mpa: f64,
    minimum_stomatal_resistance_h_per_m: f64,
    cuticular_resistance_h_per_m: f64,
    biome_boundary_resistance_h_per_m: f64,
    active_canopy_carbon_g: f64,
    stalk_volume_m3_per_g_c: f64,
};

pub const Result = struct {
    intercepted_water_volume_m3: f64,
    net_radiation_flux_megajoules_per_step: f64,
    latent_heat_flux_megajoules_per_step: f64,
    sensible_heat_flux_megajoules_per_step: f64,
    storage_heat_flux_megajoules_per_step: f64,
    convective_water_heat_flux_megajoules_per_step: f64,
    intercepted_evaporation_m3_per_step: f64,
    transpiration_m3_per_step: f64,
    ground_vapor_flux_m3_per_step: f64,
    ground_sensible_heat_flux_megajoules_per_step: f64,
    canopy_air_temperature_k: f64,
    canopy_surface_temperature_k: f64,
    canopy_surface_temperature_c: f64,
    canopy_air_vapor_concentration_m3_per_m3: f64,
    longwave_emission_megajoules_per_step: f64,
    canopy_total_water_potential_mpa: f64,
    canopy_osmotic_water_potential_mpa: f64,
    canopy_turgor_water_potential_mpa: f64,
    stomatal_resistance_h_per_m: f64,
    aerodynamic_resistance_h_per_m: f64,
    dry_canopy_heat_capacity_megajoules_per_k: f64,
};

/// UPTAKE.F 1480--1516. Initializes the scientifically inactive/small-canopy
/// state; this is not the legacy nonlinear non-convergence fallback.
pub fn calculate(inputs: Inputs) !Result {
    try validate(inputs);
    const intercepted_water =
        inputs.intercepted_water_volume_m3 +
        inputs.foliar_water_retention_m3_per_h *
            inputs.water_timestep_h_per_step;
    const temperature =
        if (inputs.canopy_height_m >=
        inputs.snow_depth_m - inputs.depth_tolerance_m)
            inputs.ambient_air_temperature_k
        else
            inputs.snow_air_temperature_k;
    const temperature_c = temperature - 273.15;
    const emitted_longwave_coefficient =
        inputs.emissivity *
        inputs.stefan_boltzmann_megajoules_per_m2_h_k4 *
        inputs.canopy_radiation_fraction *
        inputs.horizontal_cell_area_m2 *
        inputs.energy_timestep_h_per_step;
    const longwave_emission =
        emitted_longwave_coefficient *
        temperature * temperature * temperature * temperature;
    const total_potential =
        inputs.root_reference_soil_total_water_potential_mpa;
    const absolute_potential = @abs(total_potential);
    const dry_matter_fraction =
        inputs.minimum_dry_matter_fraction_g_c_per_g +
        0.10 * absolute_potential /
            (0.05 * absolute_potential + 2.0);
    const nonstructural_solute_concentration =
        inputs.canopy_nonstructural_carbon_g_per_g_c +
        inputs.canopy_nonstructural_nitrogen_g_per_g_c +
        inputs.canopy_nonstructural_phosphorus_g_per_g_c;
    const solute_molar_mass =
        36.0 + 840.0 * @max(0, nonstructural_solute_concentration);
    const osmotic_potential =
        dry_matter_fraction /
        inputs.minimum_dry_matter_fraction_g_c_per_g *
        inputs.osmotic_potential_at_zero_total_mpa -
        8.3143 * inputs.osmotic_temperature_k *
            dry_matter_fraction *
            (nonstructural_solute_concentration / solute_molar_mass +
                inputs.canopy_salt_concentration_mol_per_g_c);
    const turgor_potential = @max(0, total_potential - osmotic_potential);
    const turgor_response =
        @exp(inputs.turgor_response_shape_per_mpa * turgor_potential);
    const stomatal_resistance =
        inputs.minimum_stomatal_resistance_h_per_m +
        (inputs.cuticular_resistance_h_per_m -
            inputs.minimum_stomatal_resistance_h_per_m) *
            turgor_response;
    const dry_heat_capacity =
        2.496 * inputs.active_canopy_carbon_g *
        inputs.stalk_volume_m3_per_g_c;
    const result = Result{
        .intercepted_water_volume_m3 = intercepted_water,
        .net_radiation_flux_megajoules_per_step = 0,
        .latent_heat_flux_megajoules_per_step = 0,
        .sensible_heat_flux_megajoules_per_step = 0,
        .storage_heat_flux_megajoules_per_step = 0,
        .convective_water_heat_flux_megajoules_per_step = 0,
        .intercepted_evaporation_m3_per_step = 0,
        .transpiration_m3_per_step = 0,
        .ground_vapor_flux_m3_per_step = 0,
        .ground_sensible_heat_flux_megajoules_per_step = 0,
        .canopy_air_temperature_k = temperature,
        .canopy_surface_temperature_k = temperature,
        .canopy_surface_temperature_c = temperature_c,
        .canopy_air_vapor_concentration_m3_per_m3 = inputs.atmospheric_vapor_concentration_m3_per_m3,
        .longwave_emission_megajoules_per_step = longwave_emission,
        .canopy_total_water_potential_mpa = total_potential,
        .canopy_osmotic_water_potential_mpa = osmotic_potential,
        .canopy_turgor_water_potential_mpa = turgor_potential,
        .stomatal_resistance_h_per_m = stomatal_resistance,
        .aerodynamic_resistance_h_per_m = inputs.biome_boundary_resistance_h_per_m,
        .dry_canopy_heat_capacity_megajoules_per_k = dry_heat_capacity,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteInactiveCanopyStateResult;
    return result;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidInactiveCanopyStateInput;
    if (inputs.intercepted_water_volume_m3 < 0 or
        inputs.water_timestep_h_per_step < 0 or
        inputs.depth_tolerance_m < 0 or
        inputs.ambient_air_temperature_k <= 0 or
        inputs.snow_air_temperature_k <= 0 or
        inputs.atmospheric_vapor_concentration_m3_per_m3 < 0 or
        inputs.emissivity < 0 or
        inputs.stefan_boltzmann_megajoules_per_m2_h_k4 < 0 or
        inputs.canopy_radiation_fraction < 0 or
        inputs.horizontal_cell_area_m2 <= 0 or
        inputs.energy_timestep_h_per_step < 0 or
        inputs.minimum_dry_matter_fraction_g_c_per_g <= 0 or
        inputs.osmotic_temperature_k <= 0 or
        inputs.canopy_salt_concentration_mol_per_g_c < 0 or
        inputs.minimum_stomatal_resistance_h_per_m < 0 or
        inputs.cuticular_resistance_h_per_m < 0 or
        inputs.biome_boundary_resistance_h_per_m < 0 or
        inputs.active_canopy_carbon_g < 0 or
        inputs.stalk_volume_m3_per_g_c < 0)
        return error.InvalidInactiveCanopyStateInput;
}

fn sourceInputs() Inputs {
    return .{
        .intercepted_water_volume_m3 = 0.1,
        .foliar_water_retention_m3_per_h = 0.02,
        .water_timestep_h_per_step = 0.5,
        .canopy_height_m = 1,
        .snow_depth_m = 0.2,
        .depth_tolerance_m = 1e-12,
        .ambient_air_temperature_k = 300,
        .snow_air_temperature_k = 270,
        .atmospheric_vapor_concentration_m3_per_m3 = 0.01,
        .emissivity = 0.95,
        .stefan_boltzmann_megajoules_per_m2_h_k4 = 2.04e-10,
        .canopy_radiation_fraction = 0.2,
        .horizontal_cell_area_m2 = 10,
        .energy_timestep_h_per_step = 0.5,
        .root_reference_soil_total_water_potential_mpa = -1,
        .minimum_dry_matter_fraction_g_c_per_g = 0.16,
        .canopy_nonstructural_carbon_g_per_g_c = 0.1,
        .canopy_nonstructural_nitrogen_g_per_g_c = 0.01,
        .canopy_nonstructural_phosphorus_g_per_g_c = 0.001,
        .osmotic_potential_at_zero_total_mpa = -1.5,
        .osmotic_temperature_k = 299,
        .canopy_salt_concentration_mol_per_g_c = 0.001,
        .turgor_response_shape_per_mpa = -2,
        .minimum_stomatal_resistance_h_per_m = 0.01,
        .cuticular_resistance_h_per_m = 1.5,
        .biome_boundary_resistance_h_per_m = 2,
        .active_canopy_carbon_g = 100,
        .stalk_volume_m3_per_g_c = 4e-6,
    };
}

test "inactive canopy above snow preserves ambient source branch" {
    const inputs = sourceInputs();
    const result = try calculate(inputs);
    try std.testing.expectApproxEqAbs(@as(f64, 0.11), result.intercepted_water_volume_m3, 1e-16);
    try std.testing.expectEqual(@as(f64, 300), result.canopy_air_temperature_k);
    try std.testing.expectEqual(@as(f64, 300), result.canopy_surface_temperature_k);
    try std.testing.expectApproxEqAbs(@as(f64, 26.85), result.canopy_surface_temperature_c, 1e-13);
    try std.testing.expectEqual(@as(f64, -1), result.canopy_total_water_potential_mpa);
    try std.testing.expectEqual(@as(f64, 2), result.aerodynamic_resistance_h_per_m);
}

test "inactive canopy below snow preserves snow temperature branch" {
    var inputs = sourceInputs();
    inputs.canopy_height_m = 0.1;
    const result = try calculate(inputs);
    try std.testing.expectEqual(@as(f64, 270), result.canopy_air_temperature_k);
    try std.testing.expectEqual(@as(f64, 270), result.canopy_surface_temperature_k);
    try std.testing.expect(result.longwave_emission_megajoules_per_step > 0);
}

test "invalid runtime dry matter divisor fails explicitly" {
    var inputs = sourceInputs();
    inputs.minimum_dry_matter_fraction_g_c_per_g = 0;
    try std.testing.expectError(
        error.InvalidInactiveCanopyStateInput,
        calculate(inputs),
    );
}
