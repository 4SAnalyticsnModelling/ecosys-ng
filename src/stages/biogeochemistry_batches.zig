//! Batched and per-tile soil and surface biogeochemistry drivers.
//!
//! Extracted verbatim from `ecosys_ng.zig` so the entry point holds only
//! `main`. Declaration bodies are unchanged.

const std = @import("std");
const ecosys = @import("ecosys_ng");
const surface_litter_convergence = @import("surface_litter_convergence.zig");
const tile_kernels = @import("tile_kernels.zig");
noinline fn runSoilBiogeochemistryBatch(context: anytype) !void {
    {
        const nitrogen_parameters = context.soil_nitrogen_parameters.*;
        var reset_nitrogen_fluxes: ecosys.soil_nitrogen_flux_workspace.ResetContext = .{ .state = context.soil_nitrogen_flux_workspace };
        try tile_kernels.runScienceCellLayers(context, &reset_nitrogen_fluxes, ecosys.soil_nitrogen_flux_workspace.resetTile);
        const nutrient_zones = context.runscript.plant_nutrient_initialization;
        const zone_fractions: ecosys.solute_charge_classification.ZoneFractions = .{
            .ammonium_non_band = 1 - nutrient_zones.initial_ammonium_band_fraction,
            .ammonium_band = nutrient_zones.initial_ammonium_band_fraction,
            .nitrate_non_band = 1 - nutrient_zones.initial_nitrate_band_fraction,
            .nitrate_band = nutrient_zones.initial_nitrate_band_fraction,
            .phosphate_non_band = 1 - nutrient_zones.initial_phosphate_band_fraction,
            .phosphate_band = nutrient_zones.initial_phosphate_band_fraction,
        };
        var nitrifier_environment_context: ecosys.soil_nitrifier_environment_step.ApplyContext = .{ .result = context.soil_nitrifier_environment, .microbial_state = context.soil_microbial, .model_grid = context.grid, .gas_state = context.gas_transport, .parameters = nitrogen_parameters };
        try tile_kernels.runScienceCellLayers(context, &nitrifier_environment_context, ecosys.soil_nitrifier_environment_step.applyTile);
        var nitrification_context: ecosys.soil_nitrification_step.ApplyContext = .{ .result = context.soil_nitrogen_flux_workspace, .reactive_nitrogen = context.soil_reactive_nitrogen, .chemistry_state = context.soil_chemistry, .model_grid = context.grid, .zone_fractions = zone_fractions, .roles = context.soil_nitrifier_environment.roles, .temperature_water_activity = context.soil_nitrifier_environment.temperature_water_activity, .nitrogen_phosphorus_activity = context.soil_nitrifier_environment.nitrogen_phosphorus_activity, .aqueous_co2_activity = context.soil_nitrifier_environment.aqueous_co2_activity, .active_biomass_g_c = context.soil_nitrifier_environment.active_biomass_g_c, .microbial_active_fraction = context.soil_nitrifier_environment.microbial_active_fraction, .parameters = nitrogen_parameters.nitrification, .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol, .timestep_h = 1, .negligible_demand_g_n = context.config.absolute_tolerance };
        try tile_kernels.runScienceCellLayers(context, &nitrification_context, ecosys.soil_nitrification_step.applyTile);
        var heterotrophic_respiration_context: ecosys.soil_heterotrophic_respiration_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .reactive_nitrogen = context.soil_reactive_nitrogen,
            .organic_state = context.soil_organic,
            .microbial_state = context.soil_microbial,
            .model_grid = context.grid,
            .retention_curve = context.soil_solver_properties.retention_curve,
            .matrix_bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3,
            .matric_plus_osmotic_potential_megapascal = context.soil_hourly_workspace.matric_plus_osmotic_potential_megapascal,
            // NITRO-BIND-BACKLOG (`DISC-NITRO-001`), source `nitro.f` 556--557,
            // 918--928 and 1000--1004. Read-only inputs for the anaerobic
            // `ECHZ` product-energy feedback on the two anaerobic
            // growth-respiration requirements.
            // `soil_heterotrophic_respiration_step.applyTile` remains the SOLE
            // writer of the anaerobic publish set: no field gains a second
            // writer here, only two scalar supply ceilings change value, and
            // only when the runtime parameter set carries a
            // `soil_anaerobic_growth_energy` record. Absent that record the
            // step reproduces the legacy constants bit for bit, which its own
            // test pins. `gas_transport` is already threaded into the
            // neighbouring steps of this same phase, and source `GH2X` is
            // computed once per layer before the loops that consume it, so
            // there is no ordering inversion and no phase change.
            .gas_state = context.gas_transport,
            .aqueous_water_volume_m3 = context.grid.matrix_liquid_water_m3,
            .parameters = nitrogen_parameters,
            .timestep_h = 1,
            .negligible_amount = context.config.absolute_tolerance,
        };
        try tile_kernels.runScienceCellLayers(context, &heterotrophic_respiration_context, ecosys.soil_heterotrophic_respiration_step.applyTile);
        const surface_parameters = context.surface_gas_parameters.*;
        var soil_oxygen_context: ecosys.soil_oxygen_step.ApplyContext = .{
            .staging = context.soil_oxygen_staging,
            .result = context.soil_microbial_oxygen,
            .flux = context.soil_nitrogen_flux_workspace,
            .reactive = context.soil_reactive_nitrogen,
            .environment = context.soil_nitrifier_environment,
            .grid = context.grid,
            .matrix_bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3,
            .porosity_fraction = context.soil_solver_properties.porosity_fraction,
            .field_capacity_fraction = context.soil_field_capacity_fraction,
            .gas = context.gas_transport,
            .root_gas_parameters = context.root_gas_parameters,
            .oxygen_parameters = nitrogen_parameters.oxygen_uptake,
            .atmospheric_oxygen_g_o_per_m3 = surface_parameters.atmospheric_concentration_g_per_m3[@intFromEnum(ecosys.gas_transport.Species.oxygen)],
            .absolute_tolerance_g_o = context.config.absolute_tolerance,
            .relative_tolerance = context.config.relative_tolerance,
            .derivative_floor = 1e-14,
            .picard_relaxation = context.config.picard_relaxation,
            .gas_max_iterations = context.iteration_limits.gas_max_iterations,
        };
        try tile_kernels.runScienceCellLayers(context, &soil_oxygen_context, ecosys.soil_oxygen_step.applyTile);
        var microbial_mixing_activity_context: ecosys.soil_microbial_layer_mixing.PrepareContext = .{
            .result = context.soil_microbial_layer_mixing,
            .respiration_fluxes = context.soil_nitrogen_flux_workspace,
            .oxygen_allocation = context.soil_microbial_oxygen,
        };
        try tile_kernels.runScienceCellLayers(context, &microbial_mixing_activity_context, ecosys.soil_microbial_layer_mixing.prepareActivityTile);
        var heterotrophic_denitrification_context: ecosys.soil_heterotrophic_denitrification_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .reactive_nitrogen = context.soil_reactive_nitrogen,
            .chemistry_state = context.soil_chemistry,
            .organic_state = context.soil_organic,
            .microbial_state = context.soil_microbial,
            .gas_state = context.gas_transport,
            .model_grid = context.grid,
            .oxygen_state = context.soil_microbial_oxygen,
            .zone_fractions = zone_fractions,
            .nitrogen_parameters = nitrogen_parameters,
            .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
            .timestep_h = 1,
            .negligible_amount = context.config.absolute_tolerance,
        };
        try tile_kernels.runScienceCellLayers(context, &heterotrophic_denitrification_context, ecosys.soil_heterotrophic_denitrification_step.applyTile);
        var autotrophic_denitrification_context: ecosys.soil_autotrophic_denitrification_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .reactive_nitrogen = context.soil_reactive_nitrogen,
            .chemistry_state = context.soil_chemistry,
            .model_grid = context.grid,
            .oxygen_state = context.soil_microbial_oxygen,
            .environment = context.soil_nitrifier_environment,
            .retention_curve = context.soil_solver_properties.retention_curve,
            .matrix_bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3,
            .zone_fractions = zone_fractions,
            .parameters = nitrogen_parameters,
            .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
            .timestep_h = 1,
            .negligible_amount = context.config.absolute_tolerance,
        };
        try tile_kernels.runScienceCellLayers(context, &autotrophic_denitrification_context, ecosys.soil_autotrophic_denitrification_step.applyTile);
        var soil_microbial_nitrogen_exchange_context: ecosys.soil_microbial_nitrogen_exchange_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .reactive_nitrogen = context.soil_reactive_nitrogen,
            .microbial_state = context.soil_microbial,
            .chemistry_state = context.soil_chemistry,
            .model_grid = context.grid,
            .matric_plus_osmotic_potential_megapascal = context.soil_hourly_workspace.matric_plus_osmotic_potential_megapascal,
            .zone_fractions = zone_fractions,
            .parameters = nitrogen_parameters,
            .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
            .timestep_h = 1,
            .negligible_amount = context.config.absolute_tolerance,
        };
        try tile_kernels.runScienceCellLayers(context, &soil_microbial_nitrogen_exchange_context, ecosys.soil_microbial_nitrogen_exchange_step.applyTile);
        var soil_microbial_phosphorus_exchange_context: ecosys.soil_microbial_phosphorus_exchange_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .phosphorus_history = context.soil_microbial_phosphorus,
            .microbial_state = context.soil_microbial,
            .chemistry_state = context.soil_chemistry,
            .model_grid = context.grid,
            .matric_plus_osmotic_potential_megapascal = context.soil_hourly_workspace.matric_plus_osmotic_potential_megapascal,
            .zone_fractions = zone_fractions,
            .parameters = nitrogen_parameters,
            .timestep_h = 1,
            .negligible_amount = context.config.absolute_tolerance,
        };
        try tile_kernels.runScienceCellLayers(context, &soil_microbial_phosphorus_exchange_context, ecosys.soil_microbial_phosphorus_exchange_step.applyTile);
        var soil_microbial_maintenance_context: ecosys.soil_microbial_maintenance_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .microbial_state = context.soil_microbial,
            .organic_state = context.soil_organic,
            .chemistry_state = context.soil_chemistry,
            .model_grid = context.grid,
            .oxygen_state = context.soil_microbial_oxygen,
            .matric_plus_osmotic_potential_megapascal = context.soil_hourly_workspace.matric_plus_osmotic_potential_megapascal,
            .parameters = nitrogen_parameters,
            .timestep_h = 1,
        };
        try tile_kernels.runScienceCellLayers(context, &soil_microbial_maintenance_context, ecosys.soil_microbial_maintenance_step.applyTile);
        var autotrophic_carbon_context: ecosys.soil_autotrophic_carbon_step.ApplyContext = .{ .result = context.soil_autotrophic_carbon, .flux_workspace = context.soil_nitrogen_flux_workspace, .oxygen_state = context.soil_microbial_oxygen, .roles = context.soil_nitrifier_environment.roles, .parameters = nitrogen_parameters };
        try tile_kernels.runScienceCellLayers(context, &autotrophic_carbon_context, ecosys.soil_autotrophic_carbon_step.applyTile);
        var soil_nitrogen_fixation_context: ecosys.soil_nonsymbiotic_nitrogen_fixation_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .microbial_state = context.soil_microbial,
            .gas_state = context.gas_transport,
            .water_volume_m3 = context.grid.matrix_liquid_water_m3,
            .parameters = nitrogen_parameters,
            .timestep_h = 1,
        };
        try tile_kernels.runScienceCellLayers(context, &soil_nitrogen_fixation_context, ecosys.soil_nonsymbiotic_nitrogen_fixation_step.applyTile);
        var diagnostic_soil_fixed_n_g_n: f64 = 0;
        for (context.soil_nitrogen_flux_workspace.fixed_dinitrogen_g_n) |value| diagnostic_soil_fixed_n_g_n += value;
        std.log.debug("soil nonsymbiotic fixation: fixed_g_n={e}", .{diagnostic_soil_fixed_n_g_n});
        var soil_microbial_substrate_uptake_context: ecosys.soil_microbial_substrate_uptake_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .organic_state = context.soil_organic,
            .microbial_state = context.soil_microbial,
            .oxygen_state = context.soil_microbial_oxygen,
            .parameters = nitrogen_parameters,
            .negligible_amount = context.config.absolute_tolerance,
        };
        try tile_kernels.runScienceCellLayers(context, &soil_microbial_substrate_uptake_context, ecosys.soil_microbial_substrate_uptake_step.applyTile);
        var soil_respiration_products_context: ecosys.soil_respiration_products_step.ApplyContext = .{ .result = context.soil_respiration_products, .microbial_state = context.soil_microbial, .respiration_fluxes = context.soil_nitrogen_flux_workspace };
        try tile_kernels.runScienceCellLayers(context, &soil_respiration_products_context, ecosys.soil_respiration_products_step.applyTile);
        if (nitrogen_parameters.methane) |methane_parameters| {
            var prepare_methane_context: ecosys.soil_methane_step.PrepareContext = .{
                .result = context.soil_methane,
                .microbial_state = context.soil_microbial,
                .respiration_products = context.soil_respiration_products,
                .gas_state = context.gas_transport,
                .water_volume_m3 = context.grid.matrix_liquid_water_m3,
                .soil_temperature_k = context.grid.soil_temperature_k,
                .matric_plus_osmotic_potential_megapascal = context.soil_hourly_workspace.matric_plus_osmotic_potential_megapascal,
                .autotrophic_substrate_index = nitrogen_parameters.nitrifier_indices.autotrophic_substrate_index,
                .hydrogenotroph_population_index = methane_parameters.hydrogenotroph_population_index,
                .methanotroph_population_index = methane_parameters.methanotroph_population_index,
                .labile_biomass_fraction = nitrogen_parameters.nitrifier_environment.labile_biomass_fraction,
                .target_nitrogen_per_carbon_g_n_per_g_c = nitrogen_parameters.heterotrophic_respiration.target_nitrogen_per_carbon_g_n_per_g_c,
                .target_phosphorus_per_carbon_g_p_per_g_c = nitrogen_parameters.heterotrophic_respiration.target_phosphorus_per_carbon_g_p_per_g_c,
                .aqueous_co2_half_saturation_g_c_per_m3 = nitrogen_parameters.nitrifier_environment.aqueous_co2_half_saturation_g_c_per_m3,
                .hydrogen_product_inhibition_g_h_per_m3 = methane_parameters.hydrogen_product_inhibition_g_h_per_m3,
                .gas_constant_kj_per_mol_k = 8.3143e-3,
                .minimum_hydrogen_concentration_g_h_per_m3 = 1.0e-3,
                .hydrogen_feedback_stoichiometric_exponent = 4,
                .thermal_adaptation_offset_k = nitrogen_parameters.microbial_thermal_adaptation_offset_k,
            };
            try tile_kernels.runScienceCellLayers(context, &prepare_methane_context, ecosys.soil_methane_step.prepareTile);
            var methane_context: ecosys.soil_methane_step.ApplyContext = .{
                .result = context.soil_methane,
                .gas_state = context.gas_transport,
                .water_volume_m3 = context.grid.matrix_liquid_water_m3,
                .fermentation_hydrogen_production_g_h = context.soil_methane.fermentation_hydrogen_production_g_h,
                .acetotrophic_methane_production_g_c = context.soil_methane.acetotrophic_methane_production_g_c,
                .hydrogenotroph_active_biomass_g_c = context.soil_methane.hydrogenotroph_active_biomass_g_c,
                .methanotroph_active_biomass_g_c = context.soil_methane.methanotroph_active_biomass_g_c,
                .temperature_water_response = context.soil_methane.temperature_water_response,
                .nutrient_limitation_fraction = context.soil_methane.nutrient_limitation_fraction,
                .aqueous_co2_limitation_fraction = context.soil_methane.aqueous_co2_limitation_fraction,
                .hydrogen_feedback_energy_kj_per_mol = context.soil_methane.hydrogen_feedback_energy_kj_per_mol,
                .methanotroph_specific_oxidation_per_h = methane_parameters.methanotroph_specific_oxidation_per_h,
                .parameters = .{
                    .hydrogenotrophic = .{ .hydrogen_half_saturation_g_h_per_m3 = methane_parameters.hydrogen_half_saturation_g_h_per_m3, .specific_co2_reduction_g_c_per_g_c_h = methane_parameters.hydrogenotrophic_specific_co2_reduction_g_c_per_g_c_h, .reference_energy_yield_kj_per_g_c = methane_parameters.hydrogenotrophic_reference_energy_yield_kj_per_g_c, .growth_energy_requirement_kj_per_g_c = methane_parameters.methanogen_growth_energy_requirement_kj_per_g_c, .minimum_growth_respiration_fraction = methane_parameters.minimum_growth_respiration_fraction, .hydrogen_supply_conversion_g_c_per_g_h = methane_parameters.hydrogen_supply_conversion_g_c_per_g_h, .fermentation_hydrogen_to_pool_fraction = methane_parameters.fermentation_hydrogen_to_pool_fraction },
                    .methane_half_saturation_g_c_per_m3 = methane_parameters.methane_half_saturation_g_c_per_m3,
                    .methane_solubility_water_to_air = methane_parameters.methane_solubility_water_to_air,
                    .gas_exchange_rate_per_step = methane_parameters.gas_exchange_rate_per_step,
                    .biomass_conversion_efficiency_g_c_per_g_c = methane_parameters.methanotroph_biomass_conversion_efficiency_g_c_per_g_c,
                    .methanotroph_growth_respiration_g_c_per_g_c = methane_parameters.methanotroph_growth_respiration_g_c_per_g_c,
                },
                .timestep_h = 1,
                .solver_options = .{ .absolute_tolerance_g_c = context.config.absolute_tolerance, .relative_tolerance = context.config.relative_tolerance, .derivative_floor = 1e-14, .picard_relaxation = context.config.picard_relaxation, .gas_max_iterations = context.iteration_limits.gas_max_iterations },
            };
            try tile_kernels.runScienceCellLayers(context, &methane_context, ecosys.soil_methane_step.applyTile);
        }
        {
            const organic_runtime_parameters = context.organic_parameters;
            var soil_microbial_assimilation_context: ecosys.soil_microbial_assimilation_step.ApplyContext = .{
                .result = context.soil_nitrogen_flux_workspace,
                .microbial_state = context.soil_microbial,
                .model_grid = context.grid,
                .matric_plus_osmotic_potential_megapascal = context.soil_hourly_workspace.matric_plus_osmotic_potential_megapascal,
                .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.microbial_nitrogen_to_carbon,
                .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.microbial_phosphorus_to_carbon,
                .parameters = nitrogen_parameters,
                .timestep_h = 1,
            };
            try tile_kernels.runScienceCellLayers(context, &soil_microbial_assimilation_context, ecosys.soil_microbial_assimilation_step.applyTile);
            var soil_microbial_turnover_context: ecosys.soil_microbial_turnover_step.ApplyContext = .{
                .result = context.soil_microbial_turnover,
                .microbial_state = context.soil_microbial,
                .organic_state = context.soil_organic,
                .maintenance_fluxes = context.soil_nitrogen_flux_workspace,
                .model_grid = context.grid,
                .clay_mass_fraction = context.soil_solver_properties.clay_mass_fraction,
                .matric_plus_osmotic_potential_megapascal = context.soil_hourly_workspace.matric_plus_osmotic_potential_megapascal,
                .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.microbial_nitrogen_to_carbon,
                .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.microbial_phosphorus_to_carbon,
                .substrate_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.substrate_nitrogen_to_carbon,
                .substrate_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.substrate_phosphorus_to_carbon,
                .parameters = nitrogen_parameters,
                .timestep_h = 1,
                .negligible_carbon_g_c = context.config.absolute_tolerance,
            };
            try tile_kernels.runScienceCellLayers(context, &soil_microbial_turnover_context, ecosys.soil_microbial_turnover_step.applyTile);
        }
        {
            const organic_runtime_parameters = context.organic_parameters.*;
            var soil_organic_priming_context: ecosys.soil_organic_priming_step.ApplyContext = .{
                .result = context.soil_organic_priming,
                .organic_state = context.soil_organic,
                .microbial_state = context.soil_microbial,
                .respiration_fluxes = context.soil_nitrogen_flux_workspace,
                .soil_temperature_k = context.grid.soil_temperature_k,
                .matric_plus_osmotic_potential_megapascal = context.soil_hourly_workspace.matric_plus_osmotic_potential_megapascal,
                .thermal_adaptation_offset_k = nitrogen_parameters.microbial_thermal_adaptation_offset_k,
                .dissolved_priming_rate_per_h = organic_runtime_parameters.soil_dissolved_priming_rate_per_h,
                .microbial_priming_rate_per_h = organic_runtime_parameters.soil_microbial_priming_rate_per_h,
                .timestep_h = 1,
                .negligible_carbon_g_c = context.config.absolute_tolerance,
            };
            try tile_kernels.runScienceCellLayers(context, &soil_organic_priming_context, ecosys.soil_organic_priming_step.applyTile);
            var soil_organic_decomposition_context: ecosys.soil_organic_decomposition_step.ApplyContext = .{
                .result = context.soil_organic_decomposition,
                .organic_state = context.soil_organic,
                .microbial_state = context.soil_microbial,
                .respiration_fluxes = context.soil_nitrogen_flux_workspace,
                .priming = context.soil_organic_priming,
                .soil_temperature_k = context.grid.soil_temperature_k,
                .matrix_bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3,
                .bulk_density_megagrams_per_m3 = context.soil_solver_properties.bulk_density_megagrams_per_m3,
                .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.microbial_nitrogen_to_carbon,
                .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.microbial_phosphorus_to_carbon,
                .substrate_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.substrate_nitrogen_to_carbon,
                .substrate_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.substrate_phosphorus_to_carbon,
                .thermal_adaptation_offset_k = nitrogen_parameters.microbial_thermal_adaptation_offset_k,
                .parameters = organic_runtime_parameters.soil_organic_decomposition,
                .timestep_h = 1,
                .negligible_carbon_g_c = context.config.absolute_tolerance,
            };
            try tile_kernels.runScienceCellLayers(context, &soil_organic_decomposition_context, ecosys.soil_organic_decomposition_step.applyTile);
            var soil_organic_sorption_context: ecosys.soil_organic_sorption_step.ApplyContext = .{
                .result = context.soil_organic_sorption,
                .organic_state = context.soil_organic,
                .microbial_state = context.soil_microbial,
                .respiration_fluxes = context.soil_nitrogen_flux_workspace,
                .water_volume_m3 = context.grid.matrix_liquid_water_m3,
                .matrix_bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3,
                .bulk_density_megagrams_per_m3 = context.soil_solver_properties.bulk_density_megagrams_per_m3,
                .anion_exchange_capacity_mol_per_megagram = context.soil_solver_properties.anion_exchange_capacity_mol_per_megagram,
                .sorption_rate_per_h = organic_runtime_parameters.soil_organic_sorption_rate_per_h,
                .adsorption_coefficient = organic_runtime_parameters.soil_organic_adsorption_coefficient,
                .timestep_h = 1,
                .negligible_amount_g = context.config.absolute_tolerance,
            };
            try tile_kernels.runScienceCellLayers(context, &soil_organic_sorption_context, ecosys.soil_organic_sorption_step.applyTile);
        }
        var soil_litter_colonization_context: ecosys.soil_litter_colonization_step.ApplyContext = .{
            .result = context.soil_litter_colonization,
            .organic_state = context.soil_organic,
            .microbial_state = context.soil_microbial,
            .respiration_fluxes = context.soil_nitrogen_flux_workspace,
            .decomposition = context.soil_organic_decomposition,
            .colonization_per_g_respired_carbon = .{
                nitrogen_parameters.microbial_turnover.woody_colonization_per_g_respired_carbon,
                nitrogen_parameters.microbial_turnover.fine_litter_colonization_per_g_respired_carbon,
                nitrogen_parameters.microbial_turnover.manure_colonization_per_g_respired_carbon,
                nitrogen_parameters.microbial_turnover.particulate_colonization_per_g_respired_carbon,
                nitrogen_parameters.microbial_turnover.humus_colonization_per_g_respired_carbon,
            },
            .negligible_carbon_g_c = context.config.absolute_tolerance,
        };
        try tile_kernels.runScienceCellLayers(context, &soil_litter_colonization_context, ecosys.soil_litter_colonization_step.applyTile);
        var chemodenitrification_context: ecosys.soil_chemodenitrification_step.ApplyContext = .{ .result = context.soil_nitrogen_flux_workspace, .reactive_nitrogen = context.soil_reactive_nitrogen, .chemistry_state = context.soil_chemistry, .model_grid = context.grid, .zone_fractions = zone_fractions, .parameters = nitrogen_parameters, .timestep_h = 1 };
        try tile_kernels.runScienceCellLayers(context, &chemodenitrification_context, ecosys.soil_chemodenitrification_step.applyTile);
        var gas_aggregation_context: ecosys.soil_biogeochemical_gas_aggregation.ProcessContext = .{
            .result = context.soil_biogeochemical_gas_fluxes,
            .autotrophic_carbon = context.soil_autotrophic_carbon,
            .respiration_products = context.soil_respiration_products,
            .methane = if (nitrogen_parameters.methane != null) context.soil_methane else null,
            .oxygen = context.soil_microbial_oxygen,
            .nitrogen_fluxes = context.soil_nitrogen_flux_workspace,
            .redox_satisfaction_fraction = context.soil_redox_satisfaction_fraction,
        };
        try tile_kernels.runScienceCellLayers(context, &gas_aggregation_context, ecosys.soil_biogeochemical_gas_aggregation.aggregateProcessTile);
        var nitrogen_commit_context: ecosys.soil_nitrogen_commit.ApplyContext = .{ .reactive_nitrogen = context.soil_reactive_nitrogen, .phosphorus_history = context.soil_microbial_phosphorus, .chemistry_state = context.soil_chemistry, .gas_state = context.gas_transport, .organic_state = context.soil_organic, .microbial_state = context.soil_microbial, .flux_workspace = context.soil_nitrogen_flux_workspace, .microbial_turnover = context.soil_microbial_turnover, .litter_colonization = context.soil_litter_colonization, .organic_sorption = context.soil_organic_sorption, .organic_decomposition = context.soil_organic_decomposition, .organic_priming = context.soil_organic_priming, .respiration_products = context.soil_respiration_products, .methane = if (nitrogen_parameters.methane != null) context.soil_methane else null, .methane_autotrophic_substrate_index = nitrogen_parameters.nitrifier_indices.autotrophic_substrate_index, .methanotroph_population_index = if (nitrogen_parameters.methane) |methane_parameters| methane_parameters.methanotroph_population_index else 0, .humus_partition_by_cell = context.topsoil_humus_partition, .soil_layer_capacity = context.grid.soil_layer_capacity, .water_volume_m3 = context.grid.matrix_liquid_water_m3, .zone_fractions = zone_fractions, .oxygen_satisfaction_fraction = context.soil_microbial_oxygen.demand_satisfaction_fraction, .redox_satisfaction_fraction = context.soil_redox_satisfaction_fraction, .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol, .phosphorus_molar_mass_g_per_mol = nitrogen_parameters.microbial_mineral_exchange.phosphorus_molar_mass_g_per_mol, .tolerance_g_n = context.config.absolute_tolerance, .dynamic_salts = context.runscript.dynamic_plant_salts, .timestep_h = 1, .hourly_signed_heterotrophic_respiration_g_c = context.daily_heterotrophic_respiration.soil_hourly_signed_g_c, .hourly_carbon_dioxide_production_g_c = context.daily_heterotrophic_respiration.soil_hourly_carbon_dioxide_production_g_c };
        try tile_kernels.runScienceCellLayers(context, &nitrogen_commit_context, ecosys.soil_nitrogen_commit.applyTile);
        var autotrophic_carbon_commit_context: ecosys.soil_autotrophic_carbon_step.CommitContext = .{ .state = context.soil_autotrophic_carbon, .gas_state = context.gas_transport, .microbial_state = context.soil_microbial, .autotrophic_substrate_index = nitrogen_parameters.nitrifier_indices.autotrophic_substrate_index, .tolerance_g_c = context.config.absolute_tolerance };
        try tile_kernels.runScienceCellLayers(context, &autotrophic_carbon_commit_context, ecosys.soil_autotrophic_carbon_step.commitTile);
        var microbial_layer_mixing_context: ecosys.soil_microbial_layer_mixing.ApplyContext = .{
            .microbial_state = context.soil_microbial,
            .active_layer_count = context.grid.active_soil_layer_count,
            .layer_volume_m3 = context.soil_solver_properties.layer_volume_m3,
            .dry_bulk_density_megagrams_per_m3 = context.soil_solver_properties.bulk_density_megagrams_per_m3,
            .layer_thickness_m = context.soil_solver_properties.layer_thickness_m,
            .total_organic_carbon_g_per_megagram = context.soil_solver_properties.total_organic_carbon_g_per_megagram,
            .substrate_unlimited_oxygen_limited_activity_g_c = context.soil_microbial_layer_mixing.substrate_unlimited_oxygen_limited_activity_g_c,
            .parameters = .{ .mixing_rate_per_h = nitrogen_parameters.microbial_layer_mixing_rate_per_h, .timestep_h = 1, .minimum_mixing_layer_thickness_m = context.runscript.soil_geometry_parameters.minimum_layer_thickness_m },
        };
        try tile_kernels.runScienceCells(context, &microbial_layer_mixing_context, ecosys.soil_microbial_layer_mixing.applyTile);
    }
}

pub noinline fn runSoilBiogeochemistryBySerialTile(context: anytype) !void {
    const plan = context.tile_plan;
    for (plan.tiles, 0..) |_, tile_index| {
        const owned_cells = try plan.ownedCells(tile_index);
        context.active_tile_cells.* = owned_cells;
        defer context.active_tile_cells.* = null;
        try runSoilBiogeochemistryBatch(context);
    }
}

noinline fn runSurfaceBiogeochemistryBatch(context: anytype, surface_parameters: anytype) !void {
    var litter_water_environment_context: ecosys.surface_litter_water_environment.ApplyContext = .{
        .result = context.surface_litter_water_environment,
        .litter_geometry = context.surface_litter_geometry,
        .litter_chemistry = context.surface_litter_chemistry,
        .litter_water_m3 = context.surface_precipitation.litter_water_m3,
        .litter_temperature_k = context.grid.surface_temperature_k,
        .field_capacity_potential_megapascal = context.surface_field_capacity_potential_megapascal,
        .wilting_point_potential_megapascal = context.surface_wilting_point_potential_megapascal,
        .mean_annual_temperature_c = context.mean_annual_temperature_c_by_cell,
        .parameters = surface_parameters.litter_water_environment,
    };
    try tile_kernels.runScienceCells(context, &litter_water_environment_context, ecosys.surface_litter_water_environment.applyTile);
    var microbial_environment_context: ecosys.surface_microbial_environment_step.ApplyContext = .{
        .result = context.surface_microbial_environment,
        .litter_water_m3 = context.surface_precipitation.litter_water_m3,
        .reference_litter_volume_m3 = context.surface_litter_geometry.expanded_total_volume_m3,
        .surface_porosity_m3_per_m3 = context.surface_litter_geometry.porosity_m3_per_m3,
        .field_capacity_m3_per_m3 = context.surface_litter_geometry.field_capacity_m3_per_m3,
        .inactive_water_threshold_m3_per_m3 = context.surface_litter_water_environment.inactive_water_threshold_m3_per_m3,
        .negligible_reference_volume_m3 = context.config.absolute_tolerance,
        .litter_temperature_k = context.grid.surface_temperature_k,
        .thermal_adaptation_offset_k = context.surface_litter_water_environment.thermal_adaptation_offset_k,
    };
    try tile_kernels.runScienceCells(context, &microbial_environment_context, ecosys.surface_microbial_environment_step.applyTile);
    var microbial_respiration_context: ecosys.surface_microbial_respiration_step.ApplyContext = .{
        .result = context.surface_microbial_respiration,
        .surface_organic = context.surface_organic,
        .biologically_active_water_m3 = context.surface_microbial_environment.biologically_active_water_m3,
        .growth_temperature_response = context.surface_microbial_environment.growth_temperature_response,
        .matric_plus_osmotic_potential_megapascal = context.surface_litter_water_environment.matric_plus_osmotic_water_potential_megapascal,
        .timestep_h = 1,
        .parameters = surface_parameters.microbial_respiration,
    };
    try tile_kernels.runScienceCells(context, &microbial_respiration_context, ecosys.surface_microbial_respiration_step.applyTile);
    var microbial_oxygen_context: ecosys.surface_microbial_oxygen_driver.ApplyContext = .{
        .result = context.surface_microbial_oxygen,
        .litter_gas = context.litter_gas_transport,
        .surface_organic = context.surface_organic,
        .respiration = context.surface_microbial_respiration,
        .geometry = context.surface_litter_geometry,
        .water_environment = context.surface_litter_water_environment,
        .litter_water_m3 = context.surface_precipitation.litter_water_m3,
        .litter_water_input_m3_per_h = context.surface_precipitation.water_to_litter_m3_per_h,
        .litter_ice_m3 = context.surface_litter_ice_m3,
        .parameters = surface_parameters,
        .timestep_h = 1,
        .solver_options = .{
            .absolute_tolerance_g_o = context.config.absolute_tolerance,
            .relative_tolerance = context.config.relative_tolerance,
            .derivative_floor = 1e-14,
            .picard_relaxation = context.config.picard_relaxation,
            .gas_max_iterations = context.iteration_limits.gas_max_iterations,
        },
    };
    try tile_kernels.runScienceCells(context, &microbial_oxygen_context, ecosys.surface_microbial_oxygen_driver.applyTile);
    var microbial_maintenance_context: ecosys.surface_microbial_maintenance_step.ApplyContext = .{
        .result = context.surface_microbial_maintenance,
        .surface_organic = context.surface_organic,
        .litter_chemistry = context.surface_litter_chemistry,
        .respiration = context.surface_microbial_respiration,
        .oxygen = context.surface_microbial_oxygen,
        .environment = context.surface_microbial_environment,
        .timestep_h = 1,
        .parameters = surface_parameters.microbial_respiration,
    };
    try tile_kernels.runScienceCells(context, &microbial_maintenance_context, ecosys.surface_microbial_maintenance_step.applyTile);
    {
        const organic_runtime_parameters = context.organic_parameters;
        var nitrogen_fixation_context: ecosys.surface_nonsymbiotic_nitrogen_fixation_step.ApplyContext = .{
            .result = context.surface_nonsymbiotic_nitrogen_fixation,
            .surface_organic = context.surface_organic,
            .maintenance = context.surface_microbial_maintenance,
            .litter_gas = context.litter_gas_transport,
            .litter_water_m3 = context.surface_precipitation.litter_water_m3,
            .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.surfaceMicrobialNitrogenToCarbon(),
            .parameters = surface_parameters.microbial_respiration,
            .timestep_h = 1,
        };
        try tile_kernels.runScienceCells(context, &nitrogen_fixation_context, ecosys.surface_nonsymbiotic_nitrogen_fixation_step.applyTile);
        var diagnostic_surface_fixed_n_g_n: f64 = 0;
        for (context.surface_nonsymbiotic_nitrogen_fixation.fixed_nitrogen_g_n) |value| diagnostic_surface_fixed_n_g_n += value;
        std.log.debug("surface nonsymbiotic fixation: fixed_g_n={e}", .{diagnostic_surface_fixed_n_g_n});
        var assimilation_context: ecosys.surface_microbial_assimilation_step.ApplyContext = .{
            .result = context.surface_microbial_assimilation,
            .surface_organic = context.surface_organic,
            .growth_temperature_response = context.surface_microbial_environment.growth_temperature_response,
            .matric_plus_osmotic_potential_megapascal = context.surface_litter_water_environment.matric_plus_osmotic_water_potential_megapascal,
            .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.surfaceMicrobialNitrogenToCarbon(),
            .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.surfaceMicrobialPhosphorusToCarbon(),
            .labile_structural_fraction = surface_parameters.microbial_respiration.labile_biomass_fraction,
            .transfer_rate_per_h = surface_parameters.microbial_respiration.nonstructural_to_structural_rate_per_h,
            .timestep_h = 1,
        };
        try tile_kernels.runScienceCells(context, &assimilation_context, ecosys.surface_microbial_assimilation_step.applyTile);
        var mineral_exchange_context: ecosys.surface_microbial_mineral_exchange_step.ApplyContext = .{
            .result = context.surface_microbial_mineral_exchange,
            .surface_organic = context.surface_organic,
            .litter_chemistry = context.surface_litter_chemistry,
            .litter_water_m3 = context.surface_precipitation.litter_water_m3,
            .growth_temperature_response = context.surface_microbial_environment.growth_temperature_response,
            .matric_plus_osmotic_potential_megapascal = context.surface_litter_water_environment.matric_plus_osmotic_water_potential_megapascal,
            .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.surfaceMicrobialNitrogenToCarbon(),
            .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.surfaceMicrobialPhosphorusToCarbon(),
            .labile_biomass_fraction = surface_parameters.microbial_respiration.labile_biomass_fraction,
            .microbial_surface_area_m2_per_g_c = 4 * std.math.pi * surface_parameters.microbial_radius_m * surface_parameters.microbial_radius_m * surface_parameters.microbial_count_per_g_c,
            .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
            .parameters = surface_parameters.mineral_exchange,
            .timestep_h = 1,
        };
        try tile_kernels.runScienceCells(context, &mineral_exchange_context, ecosys.surface_microbial_mineral_exchange_step.applyTile);
        const nutrient_zones = context.runscript.plant_nutrient_initialization;
        var topsoil_exchange_context: ecosys.surface_topsoil_mineral_exchange_step.ApplyContext = .{
            .result = context.surface_topsoil_mineral_exchange,
            .model_grid = context.grid,
            .topsoil_chemistry = context.soil_chemistry,
            .zone_fractions = .{
                .ammonium_non_band = 1 - nutrient_zones.initial_ammonium_band_fraction,
                .ammonium_band = nutrient_zones.initial_ammonium_band_fraction,
                .nitrate_non_band = 1 - nutrient_zones.initial_nitrate_band_fraction,
                .nitrate_band = nutrient_zones.initial_nitrate_band_fraction,
                .phosphate_non_band = 1 - nutrient_zones.initial_phosphate_band_fraction,
                .phosphate_band = nutrient_zones.initial_phosphate_band_fraction,
            },
            .surface_organic = context.surface_organic,
            .litter_exchange = context.surface_microbial_mineral_exchange,
            .growth_temperature_response = context.surface_microbial_environment.growth_temperature_response,
            .matric_plus_osmotic_potential_megapascal = context.surface_litter_water_environment.matric_plus_osmotic_water_potential_megapascal,
            .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.surfaceMicrobialNitrogenToCarbon(),
            .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.surfaceMicrobialPhosphorusToCarbon(),
            .labile_biomass_fraction = surface_parameters.microbial_respiration.labile_biomass_fraction,
            .microbial_surface_area_m2_per_g_c = 4 * std.math.pi * surface_parameters.microbial_radius_m * surface_parameters.microbial_radius_m * surface_parameters.microbial_count_per_g_c,
            .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
            .parameters = surface_parameters.mineral_exchange,
            .timestep_h = 1,
        };
        try tile_kernels.runScienceCells(context, &topsoil_exchange_context, ecosys.surface_topsoil_mineral_exchange_step.applyTile);
        var turnover_context: ecosys.surface_microbial_turnover_step.ApplyContext = .{
            .result = context.surface_microbial_turnover,
            .surface_organic = context.surface_organic,
            .maintenance = context.surface_microbial_maintenance,
            .growth_temperature_response = context.surface_microbial_environment.growth_temperature_response,
            .matric_plus_osmotic_potential_megapascal = context.surface_litter_water_environment.matric_plus_osmotic_water_potential_megapascal,
            .humification_fraction = context.surface_humification_fraction,
            .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.surfaceMicrobialNitrogenToCarbon(),
            .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.surfaceMicrobialPhosphorusToCarbon(),
            .labile_biomass_fraction = surface_parameters.microbial_respiration.labile_biomass_fraction,
            .decomposition_density_half_saturation_g_c_per_g_c = surface_parameters.microbial_respiration.decomposition_density_half_saturation_g_c_per_g_c,
            .humus_nitrogen_per_carbon_g_n_per_g_c = organic_runtime_parameters.substrate_nitrogen_to_carbon[4],
            .humus_phosphorus_per_carbon_g_p_per_g_c = organic_runtime_parameters.substrate_phosphorus_to_carbon[4],
            .negligible_carbon_g_c = context.config.absolute_tolerance,
            .timestep_h = 1,
            .parameters = surface_parameters.microbial_turnover,
        };
        try tile_kernels.runScienceCells(context, &turnover_context, ecosys.surface_microbial_turnover_step.applyTile);
        var priming_context: ecosys.surface_organic_priming_step.ApplyContext = .{
            .result = context.surface_organic_priming,
            .surface_organic = context.surface_organic,
            .respiration = context.surface_microbial_respiration,
            .substrate_uptake = context.surface_microbial_substrate_uptake,
            .environment = context.surface_microbial_environment,
            .matric_plus_osmotic_potential_megapascal = context.surface_litter_water_environment.matric_plus_osmotic_water_potential_megapascal,
            .dissolved_priming_rate_per_h = surface_parameters.microbial_turnover.dissolved_priming_rate_per_h,
            .microbial_priming_rate_per_h = surface_parameters.microbial_turnover.microbial_priming_rate_per_h,
            .timestep_h = 1,
            .negligible_carbon_g_c = context.config.absolute_tolerance,
        };
        try tile_kernels.runScienceCells(context, &priming_context, ecosys.surface_organic_priming_step.applyTile);
        var organic_decomposition_context: ecosys.surface_organic_decomposition_step.ApplyContext = .{
            .result = context.surface_organic_decomposition,
            .surface_organic = context.surface_organic,
            .respiration = context.surface_microbial_respiration,
            .priming = context.surface_organic_priming,
            .biologically_active_water_m3 = context.surface_microbial_environment.biologically_active_water_m3,
            .litter_bulk_volume_m3 = context.surface_litter_geometry.expanded_total_volume_m3,
            .growth_temperature_response = context.surface_microbial_environment.growth_temperature_response,
            .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.surfaceMicrobialNitrogenToCarbon(),
            .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.surfaceMicrobialPhosphorusToCarbon(),
            .particulate_nitrogen_per_carbon_g_n_per_g_c = organic_runtime_parameters.substrate_nitrogen_to_carbon[3],
            .particulate_phosphorus_per_carbon_g_p_per_g_c = organic_runtime_parameters.substrate_phosphorus_to_carbon[3],
            .labile_biomass_fraction = surface_parameters.microbial_respiration.labile_biomass_fraction,
            .timestep_h = 1,
            .negligible_carbon_g_c = context.config.absolute_tolerance,
            .parameters = surface_parameters.organic_decomposition,
        };
        try tile_kernels.runScienceCells(context, &organic_decomposition_context, ecosys.surface_organic_decomposition_step.applyTile);
        var organic_sorption_context: ecosys.surface_organic_sorption_step.ApplyContext = .{
            .result = context.surface_organic_sorption,
            .surface_organic = context.surface_organic,
            .litter_water_m3 = context.surface_precipitation.litter_water_m3,
            .litter_dry_mass_megagrams = context.surface_litter_geometry.dry_mass_megagrams,
            .timestep_h = 1,
            .negligible_mass_g = context.config.absolute_tolerance,
            .parameters = surface_parameters.organic_sorption,
        };
        try tile_kernels.runScienceCells(context, &organic_sorption_context, ecosys.surface_organic_sorption_step.applyTile);
        var litter_colonization_context: ecosys.surface_litter_colonization_step.ApplyContext = .{
            .result = context.surface_litter_colonization,
            .surface_organic = context.surface_organic,
            .decomposition = context.surface_organic_decomposition,
            .parameters = surface_parameters.litter_colonization,
            .negligible_carbon_g_c = context.config.absolute_tolerance,
        };
        try tile_kernels.runScienceCells(context, &litter_colonization_context, ecosys.surface_litter_colonization_step.applyTile);
    }
    var denitrification_context: ecosys.surface_denitrification_step.ApplyContext = .{
        .result = context.surface_denitrification,
        .surface_organic = context.surface_organic,
        .litter_chemistry = context.surface_litter_chemistry,
        .litter_gas = context.litter_gas_transport,
        .litter_water_m3 = context.surface_precipitation.litter_water_m3,
        .biologically_active_water_m3 = context.surface_microbial_environment.biologically_active_water_m3,
        .growth_temperature_response = context.surface_microbial_environment.growth_temperature_response,
        .respiration = context.surface_microbial_respiration,
        .oxygen = context.surface_microbial_oxygen,
        .microbial_parameters = surface_parameters.microbial_respiration,
        .parameters = surface_parameters.denitrification,
        .chemodenitrification_parameters = surface_parameters.chemodenitrification,
        .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
        .timestep_h = 1,
    };
    try tile_kernels.runScienceCells(context, &denitrification_context, ecosys.surface_denitrification_step.applyTile);
    var denitrification_map = surface_litter_convergence.SurfaceDenitrificationRespirationMap{
        .destination_g_c = context.surface_microbial_substrate_uptake.denitrification_respiration_g_c,
        .source_g_c = context.surface_denitrification.respiration_g_c,
    };
    try tile_kernels.runScienceCells(context, &denitrification_map, surface_litter_convergence.mapSurfaceDenitrificationRespiration);
    var substrate_uptake_context: ecosys.surface_microbial_substrate_uptake_step.ApplyContext = .{
        .result = context.surface_microbial_substrate_uptake,
        .surface_organic = context.surface_organic,
        .respiration = context.surface_microbial_respiration,
        .oxygen = context.surface_microbial_oxygen,
        .maintenance = context.surface_microbial_maintenance,
        .nitrogen_fixation = context.surface_nonsymbiotic_nitrogen_fixation,
        .parameters = surface_parameters.microbial_respiration,
        .denitrification_growth_respiration_requirement_g_c_per_g_c = surface_parameters.microbial_respiration.denitrification_growth_respiration_requirement_g_c_per_g_c,
    };
    try tile_kernels.runScienceCells(context, &substrate_uptake_context, ecosys.surface_microbial_substrate_uptake_step.applyTile);
    var metabolism_commit_context: ecosys.surface_metabolism_commit.ApplyContext = .{
        .surface_organic = context.surface_organic,
        .litter_chemistry = context.surface_litter_chemistry,
        .litter_gas = context.litter_gas_transport,
        .litter_water_m3 = context.surface_precipitation.litter_water_m3,
        .respiration = context.surface_microbial_respiration,
        .oxygen = context.surface_microbial_oxygen,
        .nitrogen_fixation = context.surface_nonsymbiotic_nitrogen_fixation,
        .substrate_uptake = context.surface_microbial_substrate_uptake,
        .denitrification = context.surface_denitrification,
        .assimilation = context.surface_microbial_assimilation,
        .mineral_exchange = context.surface_microbial_mineral_exchange,
        .topsoil_exchange = context.surface_topsoil_mineral_exchange,
        .turnover = context.surface_microbial_turnover,
        .priming = context.surface_organic_priming,
        .organic_decomposition = context.surface_organic_decomposition,
        .organic_sorption = context.surface_organic_sorption,
        .litter_colonization = context.surface_litter_colonization,
        .topsoil_organic = context.soil_organic,
        .topsoil_humus_partition = context.topsoil_humus_partition,
        .topsoil_chemistry = context.soil_chemistry,
        .model_grid = context.grid,
        .zone_fractions = .{
            .ammonium_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
            .ammonium_band = context.runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
            .nitrate_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
            .nitrate_band = context.runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
            .phosphate_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
            .phosphate_band = context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
        },
        .microbial_parameters = surface_parameters.microbial_respiration,
        .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
        .phosphorus_molar_mass_g_per_mol = surface_parameters.mineral_exchange.phosphorus_molar_mass_g_per_mol,
        .tolerance = context.config.absolute_tolerance,
        .hourly_signed_heterotrophic_respiration_g_c = context.daily_heterotrophic_respiration.surface_hourly_signed_g_c,
        .hourly_carbon_dioxide_production_g_c = context.daily_heterotrophic_respiration.surface_hourly_carbon_dioxide_production_g_c,
    };
    try tile_kernels.runScienceCells(context, &metabolism_commit_context, ecosys.surface_metabolism_commit.applyTile);
    // NITRO 4168--4275 evaluates the L=0 -> NU microbial transfer only after
    // both surface and first-soil-layer biological states are accepted.
    {
        const nitrogen_parameters = context.soil_nitrogen_parameters.*;
        var microbial_mixing_context: ecosys.surface_topsoil_microbial_mixing.ApplyContext = .{
            .surface_organic = context.surface_organic,
            .soil_microbial = context.soil_microbial,
            .active_soil_layer_count = context.grid.active_soil_layer_count,
            .surface_activity_g_c_per_step = context.surface_microbial_oxygen.oxygen_limited_activity_g_c_per_step,
            .topsoil_activity_g_c_per_step = context.soil_microbial_layer_mixing.substrate_unlimited_oxygen_limited_activity_g_c,
            .surface_volume_m3 = context.surface_litter_geometry.expanded_total_volume_m3,
            .topsoil_volume_m3 = context.soil_solver_properties.layer_volume_m3,
            .surface_area_m2 = context.canopy_cell_area_m2,
            .topsoil_thickness_m = context.soil_solver_properties.layer_thickness_m,
            .surface_dry_mass_megagrams = context.surface_litter_geometry.dry_mass_megagrams,
            .topsoil_bulk_density_megagrams_per_m3 = context.soil_solver_properties.bulk_density_megagrams_per_m3,
            .topsoil_organic_carbon_g_per_megagram = context.soil_solver_properties.total_organic_carbon_g_per_megagram,
            .parameters = .{
                .mixing_rate_per_h = nitrogen_parameters.microbial_layer_mixing_rate_per_h,
                .timestep_h = 1,
                .minimum_layer_thickness_m = context.runscript.soil_geometry_parameters.minimum_layer_thickness_m,
            },
        };
        try tile_kernels.runScienceCells(context, &microbial_mixing_context, ecosys.surface_topsoil_microbial_mixing.applyTile);
    }
    {
        const reaction_parameters = context.chemistry_reaction_parameters.*;
        const organic_runtime_parameters = context.organic_parameters;
        var fertilizer_context: ecosys.surface_litter_fertilizer_step.ApplyContext = .{
            .state = context.surface_litter_fertilizer,
            .chemistry = context.surface_litter_chemistry,
            .surface_organic = context.surface_organic,
            .litter_water_m3 = context.surface_precipitation.litter_water_m3,
            .dry_litter_volume_m3 = context.surface_litter_geometry.dry_litter_volume_m3,
            .biologically_active_water_m3 = context.surface_microbial_environment.biologically_active_water_m3,
            .active_biomass_respiration_g_c_per_step = context.surface_microbial_oxygen.oxygen_limited_activity_g_c_per_step,
            .microbial_temperature_factor = context.surface_microbial_environment.growth_temperature_response,
            .litter_dry_mass_megagrams_per_g_c = organic_runtime_parameters.surface_litter_dry_mass_megagrams_per_g_c,
            .step_duration_h = 1,
            .parameters = reaction_parameters.surface_fertilizer,
            .diagnostics = context.surface_litter_fertilizer_diagnostics,
        };
        try tile_kernels.runScienceCells(context, &fertilizer_context, ecosys.surface_litter_fertilizer_step.applyTile);
    }
}

pub noinline fn runSurfaceBiogeochemistryBySerialTile(context: anytype, surface_parameters: anytype) !void {
    const plan = context.tile_plan;
    for (plan.tiles, 0..) |_, tile_index| {
        const owned_cells = try plan.ownedCells(tile_index);
        context.active_tile_cells.* = owned_cells;
        defer context.active_tile_cells.* = null;
        try runSurfaceBiogeochemistryBatch(context, surface_parameters);
    }
}
