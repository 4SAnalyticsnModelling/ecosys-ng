//! Hourly soil chemistry convergence loop (NITRO/SOLUTE equivalent).
//!
//! Extracted verbatim from `ecosys_ng.zig` so the entry point holds only
//! `main`. Declaration bodies are unchanged.

const std = @import("std");
const ecosys = @import("ecosys_ng");
pub fn convergeHourlySoilChemistry(
    context: anytype,
    fertilizer_band_hour: ecosys.fertilizer_band_phase_coordinator.HourToken,
    failure_report: ?ecosys.solute_failure_reporter.Request,
) !void {
    const nutrient_zones = context.runscript.plant_nutrient_initialization;
    if (context.executed_weather_hours.* < 24) {
        var diagnostic_root_n_g: f64 = 0;
        if (@hasField(@TypeOf(context), "root_soil_ammonia_exchange_publication_state")) {
            const publication = context.root_soil_ammonia_exchange_publication_state;
            for (publication.non_band_exchange_g_n_per_h_by_layer, publication.band_exchange_g_n_per_h_by_layer) |non_band, band|
                diagnostic_root_n_g += @max(0, non_band) + @max(0, band);
        }
        if (@hasField(@TypeOf(context), "root_nutrient_uptake_publication_state")) {
            const publication = context.root_nutrient_uptake_publication_state;
            diagnostic_root_n_g += diagnosticPositiveRootNitrogenUptake_g(publication);
        }
        std.log.debug("soil chemistry root nitrogen deduction: g_n={e}", .{diagnostic_root_n_g});
    }
    try ecosys.mineral_fertilizer_inventory.publishWetted(
        context.mineral_fertilizer_inventory,
        context.soil_chemistry,
        context.surface_litter_chemistry,
        context.grid.matrix_liquid_water_m3,
        context.surface_precipitation.litter_water_m3,
        .{
            .ammonium_non_band = 1 - nutrient_zones.initial_ammonium_band_fraction,
            .ammonium_band = nutrient_zones.initial_ammonium_band_fraction,
            .nitrate_non_band = 1 - nutrient_zones.initial_nitrate_band_fraction,
            .nitrate_band = nutrient_zones.initial_nitrate_band_fraction,
            .phosphate_non_band = 1 - nutrient_zones.initial_phosphate_band_fraction,
            .phosphate_band = nutrient_zones.initial_phosphate_band_fraction,
        },
        context.config.absolute_tolerance,
    );
    {
        const chemistry_parameters = context.chemistry_reaction_parameters.*;
        const nitrogen_parameters = context.soil_nitrogen_parameters.*;
        const nutrient_offset_k =
            nitrogen_parameters.microbial_thermal_adaptation_offset_k;
        var staged_soil_fertilizer = try context.allocator.dupe(
            @TypeOf(context.soil_fertilizer_inventory.soil[0]),
            context.soil_fertilizer_inventory.soil,
        );
        defer context.allocator.free(staged_soil_fertilizer);
        var staged_aqueous = try context.allocator.dupe(
            @TypeOf(context.soil_chemistry.aqueous[0]),
            context.soil_chemistry.aqueous,
        );
        defer context.allocator.free(staged_aqueous);
        const staged_gaseous_mass_g = try context.allocator.dupe(
            f64,
            context.gas_transport.gaseous_mass_g,
        );
        defer context.allocator.free(staged_gaseous_mass_g);
        const staged_urease_inhibition = try context.allocator.dupe(
            f64,
            context.soil_fertilizer_inventory
                .current_urease_inhibition_fraction,
        );
        defer context.allocator.free(staged_urease_inhibition);

        for (0..context.grid.cell_count) |cell| {
            const first_layer = cell * context.grid.soil_layer_capacity;
            for (0..context.grid.active_soil_layer_count[cell]) |layer_within_cell| {
                const layer = first_layer + layer_within_cell;
                const water_volume_m3 = context.grid.matrix_liquid_water_m3[layer];
                const soil_mass_megagrams = context.soil_solver_properties.matrix_bulk_volume_m3[layer] * context.soil_solver_properties.bulk_density_megagrams_per_m3[layer];
                const bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3[layer];
                if (!std.math.isFinite(water_volume_m3) or water_volume_m3 < 0 or !std.math.isFinite(soil_mass_megagrams) or soil_mass_megagrams <= 0) return error.InvalidHourlySoilChemistryGeometry;
                if (water_volume_m3 <= context.config.absolute_tolerance) continue;
                if (!std.math.isFinite(bulk_volume_m3) or bulk_volume_m3 <= 0) return error.InvalidHourlySoilChemistryGeometry;
                const water_content_m3_per_m3 = water_volume_m3 / bulk_volume_m3;
                const biologically_active_water_m3 = context.soil_nitrogen_flux_workspace.layer_biologically_active_water_m3[layer];
                if (!std.math.isFinite(biologically_active_water_m3) or biologically_active_water_m3 < 0)
                    return error.InvalidHourlySoilChemistryGeometry;
                const temperature_response = try ecosys.soil_microbial_metabolism.growthTemperatureResponse(
                    context.grid.soil_temperature_k[layer],
                    nutrient_offset_k,
                );
                // SOLUTE legacy mapping:
                // TOQCK <- total maintenance respiration from NITRO workspace
                // VOLQ <- biologically active soil water from NITRO workspace
                // TFNQ <- microbial-temperature growth factor from soil temperature state
                const toqck_g_c_per_step = context.soil_nitrogen_flux_workspace.total_maintenance_respiration_g_c[layer];
                const volq_m3 = biologically_active_water_m3;
                const tfnq = temperature_response;
                const fertilizer_parameters = try chemistry_parameters.surface_fertilizer.forFormulation(
                    context.soil_fertilizer_inventory.formulation[layer],
                    1,
                );
                const hydrolysis = try ecosys.soil_fertilizer_dissolution.ureaHydrolysis(
                    .{
                        .broadcast_urea_mol_n = context.soil_fertilizer_inventory.soil[layer].broadcast_urea_mol_n,
                        .banded_urea_mol_n = context.soil_fertilizer_inventory.soil[layer].banded_urea_mol_n,
                        .soil_mass_megagrams = soil_mass_megagrams,
                        .water_volume_m3 = water_volume_m3,
                        .biologically_active_water_volume_m3 = volq_m3,
                        .total_microbial_respiration_activity_g_c_per_step = toqck_g_c_per_step,
                        .temperature_response = tfnq,
                        .initial_inhibitor_activity = context.soil_fertilizer_inventory.initial_urease_inhibition_fraction[layer],
                        .current_inhibitor_activity = context.soil_fertilizer_inventory.current_urease_inhibition_fraction[layer],
                        .timestep_h = 1,
                    },
                    .{
                        .minimum_half_saturation_mol_n_per_megagram = fertilizer_parameters.minimum_urea_half_saturation_mol_n_per_megagram,
                        .microbial_activity_inhibition_g_c_per_m3_h = fertilizer_parameters.microbial_activity_inhibition_g_c_per_m3_per_h,
                        .specific_hydrolysis_mol_n_per_g_c_h = fertilizer_parameters.specific_urea_hydrolysis_mol_n_per_g_c,
                        .inhibitor_decline_rate_per_h = fertilizer_parameters.urease_inhibition_decline_fraction_per_step,
                        .negligible_amount = context.config.absolute_tolerance,
                    },
                );
                const dissolved = try ecosys.soil_fertilizer_dissolution.dissolution(
                    staged_soil_fertilizer[layer],
                    hydrolysis,
                    .{
                        .ammonium_non_band = 1 - nutrient_zones.initial_ammonium_band_fraction,
                        .ammonium_band = nutrient_zones.initial_ammonium_band_fraction,
                        .nitrate_non_band = 1 - nutrient_zones.initial_nitrate_band_fraction,
                        .nitrate_band = nutrient_zones.initial_nitrate_band_fraction,
                    },
                    .{
                        .ammonium_per_h = fertilizer_parameters.ammonium_dissolution_fraction_per_step,
                        .ammonia_per_h = fertilizer_parameters.ammonia_dissolution_fraction_per_step,
                        .nitrate_per_h = fertilizer_parameters.nitrate_dissolution_fraction_per_step,
                    },
                    water_content_m3_per_m3,
                    1,
                );
                const ammonia_gas_index = try ecosys.gas_transport.massIndex(
                    layer,
                    .ammonia,
                    context.gas_transport.cell_count,
                );
                try ecosys.soil_fertilizer_dissolution.commitToRecipients(
                    &staged_soil_fertilizer[layer],
                    &staged_aqueous[layer],
                    &staged_gaseous_mass_g[ammonia_gas_index],
                    dissolved,
                    .{
                        .ammonium_non_band = 1 - nutrient_zones.initial_ammonium_band_fraction,
                        .ammonium_band = nutrient_zones.initial_ammonium_band_fraction,
                        .nitrate_non_band = 1 - nutrient_zones.initial_nitrate_band_fraction,
                        .nitrate_band = nutrient_zones.initial_nitrate_band_fraction,
                    },
                    water_volume_m3,
                    context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                );
                staged_urease_inhibition[layer] =
                    hydrolysis.next_inhibitor_activity;
            }
        }
        @memcpy(context.soil_fertilizer_inventory.soil, staged_soil_fertilizer);
        @memcpy(context.soil_chemistry.aqueous, staged_aqueous);
        @memcpy(context.gas_transport.gaseous_mass_g, staged_gaseous_mass_g);
        @memcpy(
            context.soil_fertilizer_inventory
                .current_urease_inhibition_fraction,
            staged_urease_inhibition,
        );

        for (0..context.grid.cell_count) |cell| {
            const first_layer = cell * context.grid.soil_layer_capacity;
            for (0..context.grid.active_soil_layer_count[cell]) |layer_within_cell| {
                const layer = first_layer + layer_within_cell;
                const water_volume_m3 = context.grid.matrix_liquid_water_m3[layer];
                const soil_mass_megagrams = context.soil_solver_properties.matrix_bulk_volume_m3[layer] * context.soil_solver_properties.bulk_density_megagrams_per_m3[layer];
                if (!std.math.isFinite(water_volume_m3) or water_volume_m3 < 0 or !std.math.isFinite(soil_mass_megagrams) or soil_mass_megagrams <= 0) return error.InvalidHourlySoilChemistryGeometry;
                if (water_volume_m3 <= context.config.absolute_tolerance) continue;
                // SOLUTE lines 376-379 (TUPN3S/TUPN3B): deduct previous-hour
                // root NH3 uptake from initial aqueous NH3 concentrations before
                // SOLUTE equilibration. Positive exchange = soil losing N to root.
                if (@hasField(@TypeOf(context), "root_soil_ammonia_exchange_publication_state")) {
                    const pub_state = context.root_soil_ammonia_exchange_publication_state;
                    const g_n_nonband = pub_state.non_band_exchange_g_n_per_h_by_layer[layer];
                    const g_n_band = pub_state.band_exchange_g_n_per_h_by_layer[layer];
                    if (g_n_nonband > 0) {
                        context.soil_chemistry.aqueous[layer].ammonia_non_band = @max(
                            0.0,
                            context.soil_chemistry.aqueous[layer].ammonia_non_band -
                                g_n_nonband / (14.0 * water_volume_m3),
                        );
                    }
                    if (g_n_band > 0) {
                        context.soil_chemistry.aqueous[layer].ammonia_band = @max(
                            0.0,
                            context.soil_chemistry.aqueous[layer].ammonia_band -
                                g_n_band / (14.0 * water_volume_m3),
                        );
                    }
                }
                // SOLUTE lines 374-375 (TUPNH4/TUPNHB): deduct previous-hour
                // root NH4 uptake from initial aqueous NH4 concentrations before
                // SOLUTE equilibration. Index 0 = non-band NH4, index 4 = band NH4.
                if (@hasField(@TypeOf(context), "root_nutrient_uptake_publication_state")) {
                    const pub_state = context.root_nutrient_uptake_publication_state;
                    const g_n_nonband = pub_state.uptake_g_element_per_h_by_nutrient_and_layer[0][layer];
                    const g_n_band = pub_state.uptake_g_element_per_h_by_nutrient_and_layer[4][layer];
                    if (g_n_nonband > 0) {
                        context.soil_chemistry.aqueous[layer].ammonium_non_band = @max(
                            0.0,
                            context.soil_chemistry.aqueous[layer].ammonium_non_band -
                                g_n_nonband / (14.0 * water_volume_m3),
                        );
                    }
                    if (g_n_band > 0) {
                        context.soil_chemistry.aqueous[layer].ammonium_band = @max(
                            0.0,
                            context.soil_chemistry.aqueous[layer].ammonium_band -
                                g_n_band / (14.0 * water_volume_m3),
                        );
                    }
                }
                var parameters = context.soil_chemistry_layer_parameters[layer];
                const prepared_zones = try ecosys.soil_fertilizer_dissolution
                    .prepareLayerZones(.{
                    .water_volume_m3 = water_volume_m3,
                    .soil_mass_megagrams = soil_mass_megagrams,
                    .soil_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3[layer],
                    .fractions = .{
                        .ammonium_non_band = parameters.fractions.ammonium_non_band,
                        .ammonium_band = parameters.fractions.ammonium_band,
                        .nitrate_non_band = parameters.fractions.nitrate_non_band,
                        .nitrate_band = parameters.fractions.nitrate_band,
                        .phosphate_non_band = parameters.fractions.phosphate_non_band,
                        .phosphate_band = parameters.fractions.phosphate_band,
                    },
                    .positive_soil_mass_threshold_megagrams = context.config.absolute_tolerance,
                });
                const shared_ratio = try ecosys.soil_fertilizer_dissolution
                    .normalizationBasisPerWaterVolume(
                    prepared_zones.whole_layer_normalization_basis,
                    water_volume_m3,
                );
                parameters.cation_exchange_water_ratios = .{
                    .shared_megagrams_per_m3 = shared_ratio,
                    .ammonium_non_band_megagrams_per_m3 = try ecosys.soil_fertilizer_dissolution.normalizationBasisPerWaterVolume(
                        prepared_zones.ammonium_non_band_normalization_basis,
                        prepared_zones.ammonium_non_band_water_m3,
                    ),
                    .ammonium_band_megagrams_per_m3 = try ecosys.soil_fertilizer_dissolution.normalizationBasisPerWaterVolume(
                        prepared_zones.ammonium_band_normalization_basis,
                        prepared_zones.ammonium_band_water_m3,
                    ),
                };
                parameters.non_band_phosphate_soil_mass_per_water_volume_megagrams_per_m3 = try ecosys.soil_fertilizer_dissolution.normalizationBasisPerWaterVolume(
                    prepared_zones.phosphate_non_band_normalization_basis,
                    prepared_zones.phosphate_non_band_water_m3,
                );
                parameters.band_phosphate_soil_mass_per_water_volume_megagrams_per_m3 = try ecosys.soil_fertilizer_dissolution.normalizationBasisPerWaterVolume(
                    prepared_zones.phosphate_band_normalization_basis,
                    prepared_zones.phosphate_band_water_m3,
                );
                parameters.total_carboxyl_sites_mol_per_megagram = context.chemistry_reaction_parameters.*.surface_litter.carboxyl_sites_mol_per_megagram_c *
                    1.0e-6 * context.soil_solver_properties.total_organic_carbon_g_per_megagram[layer];
                // SOLUTE line 2934: CCO21 = CCO2S/12 = (CO2S/VOLW)/12 mol/m³ water.
                // Without dissolved CO2, bicarbonate/carbonate are zero and the
                // charge balance diverges to pH ~2.7 from the very first hour.
                {
                    const co2_idx = try ecosys.gas_transport.massIndex(layer, .carbon_dioxide, context.gas_transport.cell_count);
                    context.soil_chemistry.aqueous[layer].carbon_dioxide = @max(0, context.gas_transport.dissolved_mass_g[co2_idx] / (12.0 * water_volume_m3));
                }
                const solver_options: ecosys.solute_reaction_solver.Options = .{
                    .absolute_tolerance = context.config.absolute_tolerance,
                    .relative_tolerance = context.config.relative_tolerance,
                    .picard_relaxation = context.config.picard_relaxation,
                    .max_iterations = context.iteration_limits.solute_reaction_max_iterations,
                };
                _ = ecosys.solute_reaction_solver.solveCellWithWorkspace(
                    context.soil_chemistry_solver_workspace,
                    context.soil_chemistry,
                    layer,
                    parameters,
                    solver_options,
                ) catch |err| switch (err) {
                    // Legacy STARTE/SOLUTE retain their state after exactly MRXN
                    // cycles without requiring convergence. For hourly runs over
                    // frozen soil layers (concentrated chemistry, small water
                    // volume), non-convergence is expected and the transactional
                    // rollback already restores the pre-solve cell state. Log and
                    // continue; the layer re-attempts equilibration next hour as
                    // conditions change.
                    error.SoluteReactionSolverStagnated,
                    error.SoluteReactionSolverDidNotConverge,
                    => std.log.warn(
                        "SOLUTE hourly non-convergence retained: cell={d} layer={d} water_m3={e} err={s}",
                        .{ cell, layer_within_cell, water_volume_m3, @errorName(err) },
                    ),
                    else => {
                        if (failure_report) |report| {
                            var contextual_report = report;
                            contextual_report.context.global_cell_id =
                                @intCast(cell);
                            contextual_report.context.soil_layer_id =
                                @intCast(layer_within_cell);
                            contextual_report.context.packed_cell_index =
                                @intCast(layer);
                            return ecosys.solute_failure_reporter
                                .reportPreservingSolverError(
                                context.allocator,
                                contextual_report,
                                context.soil_chemistry,
                                layer,
                                parameters,
                                solver_options,
                                err,
                            );
                        }
                        return err;
                    },
                };
            }
        }
        // After solute solver equilibrates CO2(aq) ↔ HCO3⁻, sync gas_state
        // dissolved CO2 with the post-equilibration chemistry CO2(aq).  Without
        // this update, the inventory counts pre-equilibration dissolved CO2 from
        // gas_state alongside post-equilibration micropore HCO3⁻ from
        // exportChemistry, creating a systematic carbon balance drift equal to
        // the equilibrium shift in each layer.
        for (0..context.grid.cell_count) |sync_cell| {
            const sync_first = sync_cell * context.grid.soil_layer_capacity;
            for (0..context.grid.active_soil_layer_count[sync_cell]) |sync_local_layer| {
                const sync_layer = sync_first + sync_local_layer;
                const sync_water_m3 = context.grid.matrix_liquid_water_m3[sync_layer];
                if (sync_water_m3 <= context.config.absolute_tolerance) continue;
                const sync_co2_idx = try ecosys.gas_transport.massIndex(sync_layer, .carbon_dioxide, context.gas_transport.cell_count);
                context.gas_transport.dissolved_mass_g[sync_co2_idx] = context.soil_chemistry.aqueous[sync_layer].carbon_dioxide * 12.0 * sync_water_m3;
            }
        }
        for (0..context.grid.cell_count) |cell| {
            const active_layer_count = context.grid.active_soil_layer_count[cell];
            if (active_layer_count == 0) continue;
            if (active_layer_count > context.grid.soil_layer_capacity)
                return error.InvalidHourlySoilChemistryGeometry;
            const first_layer = cell * context.grid.soil_layer_capacity;
            var layer_properties = try context.allocator.alloc(
                ecosys.fertilizer_band_nitrate_phosphate.LayerProperties,
                active_layer_count,
            );
            defer context.allocator.free(layer_properties);
            const nitrate_geometry = try context.fertilizer_band.geometry(
                cell,
                .nitrate,
            );
            const phosphate_geometry = try context.fertilizer_band.geometry(
                cell,
                .phosphate,
            );
            const nitrate_row_width_m = nitrate_geometry.row_spacing_m;
            const phosphate_row_width_m = phosphate_geometry.row_spacing_m;
            const nitrate_application: ecosys.fertilizer_band_nitrate_phosphate.BandApplication = if (nitrate_row_width_m > 0) .banded else .unbanded;
            const phosphate_application: ecosys.fertilizer_band_nitrate_phosphate.BandApplication = if (phosphate_row_width_m > 0) .banded else .unbanded;
            var nonband_fractional_change = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(nonband_fractional_change);
            var layer_top_depth_m: f64 = 0;
            for (0..active_layer_count) |local_layer| {
                const global_layer = first_layer + local_layer;
                const layer_thickness = context.soil_solver_properties.layer_thickness_m[global_layer];
                if (!std.math.isFinite(layer_thickness) or layer_thickness < 0)
                    return error.InvalidHourlySoilChemistryGeometry;
                const water_volume_m3 = context.grid.matrix_liquid_water_m3[global_layer];
                const pore_capacity_m3 = context.grid.matrix_pore_capacity_m3[global_layer];
                if (!std.math.isFinite(water_volume_m3) or water_volume_m3 < 0 or !std.math.isFinite(pore_capacity_m3) or pore_capacity_m3 < 0)
                    return error.InvalidHourlySoilChemistryGeometry;
                const water_fraction = if (pore_capacity_m3 > 0)
                    std.math.clamp(water_volume_m3 / pore_capacity_m3, 0, 1)
                else
                    0;
                const top_depth_m = layer_top_depth_m;
                layer_top_depth_m += layer_thickness;
                layer_properties[local_layer] = .{
                    .top_depth_m = top_depth_m,
                    .bottom_depth_m = layer_top_depth_m,
                    .thickness_m = layer_thickness,
                    .tortuosity = context.runscript.root_nutrient_parameters.liquid_tortuosity_coefficient * water_fraction * water_fraction,
                    .nitrate_diffusivity_m2_h = try context.runscript.root_nutrient_parameters.diffusivityM2PerH(
                        1,
                        context.grid.soil_temperature_k[global_layer],
                    ),
                    .phosphate_diffusivity_m2_h = try context.runscript.root_nutrient_parameters.diffusivityM2PerH(
                        2,
                        context.grid.soil_temperature_k[global_layer],
                    ),
                };
            }
            const nitrate_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(nitrate_nonband_pools);
            const nitrate_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(nitrate_band_pools);
            const nitrite_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(nitrite_nonband_pools);
            const nitrite_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(nitrite_band_pools);
            const fertilizer_nitrate_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(fertilizer_nitrate_nonband_pools);
            const fertilizer_nitrate_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(fertilizer_nitrate_band_pools);
            const hydrogen_phosphate_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(hydrogen_phosphate_nonband_pools);
            const hydrogen_phosphate_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(hydrogen_phosphate_band_pools);
            const dihydrogen_phosphate_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(dihydrogen_phosphate_nonband_pools);
            const dihydrogen_phosphate_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(dihydrogen_phosphate_band_pools);
            const adsorbed_oh0_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(adsorbed_oh0_nonband_pools);
            const adsorbed_oh0_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(adsorbed_oh0_band_pools);
            const adsorbed_oh1_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(adsorbed_oh1_nonband_pools);
            const adsorbed_oh1_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(adsorbed_oh1_band_pools);
            const adsorbed_oh2_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(adsorbed_oh2_nonband_pools);
            const adsorbed_oh2_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(adsorbed_oh2_band_pools);
            const adsorbed_hpo4_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(adsorbed_hpo4_nonband_pools);
            const adsorbed_hpo4_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(adsorbed_hpo4_band_pools);
            const adsorbed_h2po4_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(adsorbed_h2po4_nonband_pools);
            const adsorbed_h2po4_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(adsorbed_h2po4_band_pools);
            const aluminum_phosphate_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(aluminum_phosphate_nonband_pools);
            const aluminum_phosphate_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(aluminum_phosphate_band_pools);
            const iron_phosphate_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(iron_phosphate_nonband_pools);
            const iron_phosphate_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(iron_phosphate_band_pools);
            const dicalcium_phosphate_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(dicalcium_phosphate_nonband_pools);
            const dicalcium_phosphate_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(dicalcium_phosphate_band_pools);
            const hydroxyapatite_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(hydroxyapatite_nonband_pools);
            const hydroxyapatite_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(hydroxyapatite_band_pools);
            const monocalcium_phosphate_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(monocalcium_phosphate_nonband_pools);
            const monocalcium_phosphate_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(monocalcium_phosphate_band_pools);
            const phosphate_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(phosphate_nonband_pools);
            const phosphate_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(phosphate_band_pools);
            const phosphoric_acid_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(phosphoric_acid_nonband_pools);
            const phosphoric_acid_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(phosphoric_acid_band_pools);
            const iron_hpo4_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(iron_hpo4_nonband_pools);
            const iron_hpo4_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(iron_hpo4_band_pools);
            const iron_h2po4_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(iron_h2po4_nonband_pools);
            const iron_h2po4_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(iron_h2po4_band_pools);
            const calcium_hpo4_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(calcium_hpo4_nonband_pools);
            const calcium_hpo4_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(calcium_hpo4_band_pools);
            const calcium_h2po4_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(calcium_h2po4_nonband_pools);
            const calcium_h2po4_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(calcium_h2po4_band_pools);
            const calcium_phosphate_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(calcium_phosphate_nonband_pools);
            const calcium_phosphate_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(calcium_phosphate_band_pools);
            const magnesium_hpo4_nonband_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(magnesium_hpo4_nonband_pools);
            const magnesium_hpo4_band_pools = try context.allocator.alloc(f64, active_layer_count);
            defer context.allocator.free(magnesium_hpo4_band_pools);

            for (0..active_layer_count) |local_layer| {
                const global_layer = first_layer + local_layer;
                const layer = context.soil_chemistry.aqueous[global_layer];
                const non_band = context.soil_chemistry.non_band_phosphate[global_layer];
                const band = context.soil_chemistry.band_phosphate[global_layer];
                const water_volume_m3 = context.grid.matrix_liquid_water_m3[global_layer];
                if (!std.math.isFinite(water_volume_m3) or water_volume_m3 < 0)
                    return error.InvalidHourlySoilChemistryGeometry;
                const soil_mass_megagrams = context.soil_solver_properties.matrix_bulk_volume_m3[global_layer] * context.soil_solver_properties.bulk_density_megagrams_per_m3[global_layer];
                if (!std.math.isFinite(soil_mass_megagrams) or soil_mass_megagrams <= 0)
                    return error.InvalidHourlySoilChemistryGeometry;
                const nitrate_nonband_g_n = layer.nitrate_non_band * water_volume_m3 * 14.0;
                const nitrate_band_g_n = layer.nitrate_band * water_volume_m3 * 14.0;
                const nitrite_nonband_g_n = context.soil_reactive_nitrogen.non_band_nitrite_g_n[global_layer];
                const nitrite_band_g_n = context.soil_reactive_nitrogen.band_nitrite_g_n[global_layer];
                if (!std.math.isFinite(nitrate_nonband_g_n) or !std.math.isFinite(nitrate_band_g_n) or !std.math.isFinite(nitrite_nonband_g_n) or !std.math.isFinite(nitrite_band_g_n))
                    return error.InvalidHourlySoilChemistryGeometry;
                const fertilizer_nitrate_nonband_g_n = context.soil_fertilizer_inventory.soil[global_layer].broadcast_nitrate_mol_n * 14.0;
                const fertilizer_nitrate_band_g_n = context.soil_fertilizer_inventory.soil[global_layer].banded_nitrate_mol_n * 14.0;
                if (!std.math.isFinite(fertilizer_nitrate_nonband_g_n) or !std.math.isFinite(fertilizer_nitrate_band_g_n))
                    return error.InvalidHourlySoilChemistryGeometry;
                const hydrogen_phosphate_nonband_mol = non_band.dissolved_hpo4_mol_p_per_m3 * water_volume_m3;
                const hydrogen_phosphate_band_mol = band.dissolved_hpo4_mol_p_per_m3 * water_volume_m3;
                const dihydrogen_phosphate_nonband_mol = non_band.dissolved_h2po4_mol_p_per_m3 * water_volume_m3;
                const dihydrogen_phosphate_band_mol = band.dissolved_h2po4_mol_p_per_m3 * water_volume_m3;
                if (!(std.math.isFinite(hydrogen_phosphate_nonband_mol) and std.math.isFinite(hydrogen_phosphate_band_mol) and std.math.isFinite(dihydrogen_phosphate_nonband_mol) and std.math.isFinite(dihydrogen_phosphate_band_mol)))
                    return error.InvalidHourlySoilChemistryGeometry;
                const adsorbed_oh0_nonband_mol = non_band.deprotonated_site_mol_per_megagram * soil_mass_megagrams;
                const adsorbed_oh0_band_mol = band.deprotonated_site_mol_per_megagram * soil_mass_megagrams;
                const adsorbed_oh1_nonband_mol = non_band.hydroxyl_site_mol_per_megagram * soil_mass_megagrams;
                const adsorbed_oh1_band_mol = band.hydroxyl_site_mol_per_megagram * soil_mass_megagrams;
                const adsorbed_oh2_nonband_mol = non_band.protonated_site_mol_per_megagram * soil_mass_megagrams;
                const adsorbed_oh2_band_mol = band.protonated_site_mol_per_megagram * soil_mass_megagrams;
                const adsorbed_hpo4_nonband_mol = non_band.adsorbed_hpo4_mol_p_per_megagram * soil_mass_megagrams;
                const adsorbed_hpo4_band_mol = band.adsorbed_hpo4_mol_p_per_megagram * soil_mass_megagrams;
                const adsorbed_h2po4_nonband_mol = non_band.adsorbed_h2po4_mol_p_per_megagram * soil_mass_megagrams;
                const adsorbed_h2po4_band_mol = band.adsorbed_h2po4_mol_p_per_megagram * soil_mass_megagrams;
                const aluminum_phosphate_nonband_mol = non_band.aluminum_phosphate_solid_mol_per_m3 * water_volume_m3;
                const aluminum_phosphate_band_mol = band.aluminum_phosphate_solid_mol_per_m3 * water_volume_m3;
                const iron_phosphate_nonband_mol = non_band.iron_phosphate_solid_mol_per_m3 * water_volume_m3;
                const iron_phosphate_band_mol = band.iron_phosphate_solid_mol_per_m3 * water_volume_m3;
                const dicalcium_phosphate_nonband_mol = non_band.dicalcium_phosphate_solid_mol_per_m3 * water_volume_m3;
                const dicalcium_phosphate_band_mol = band.dicalcium_phosphate_solid_mol_per_m3 * water_volume_m3;
                const hydroxyapatite_nonband_mol = non_band.hydroxyapatite_solid_mol_per_m3 * water_volume_m3;
                const hydroxyapatite_band_mol = band.hydroxyapatite_solid_mol_per_m3 * water_volume_m3;
                const monocalcium_phosphate_nonband_mol = non_band.monocalcium_phosphate_solid_mol_per_m3 * water_volume_m3;
                const monocalcium_phosphate_band_mol = band.monocalcium_phosphate_solid_mol_per_m3 * water_volume_m3;
                if (!std.math.isFinite(adsorbed_oh0_nonband_mol) or !std.math.isFinite(adsorbed_oh0_band_mol) or !std.math.isFinite(adsorbed_oh1_nonband_mol) or !std.math.isFinite(adsorbed_oh1_band_mol) or !std.math.isFinite(adsorbed_oh2_nonband_mol) or !std.math.isFinite(adsorbed_oh2_band_mol) or !std.math.isFinite(adsorbed_hpo4_nonband_mol) or !std.math.isFinite(adsorbed_hpo4_band_mol) or !std.math.isFinite(adsorbed_h2po4_nonband_mol) or !std.math.isFinite(adsorbed_h2po4_band_mol) or !std.math.isFinite(aluminum_phosphate_nonband_mol) or !std.math.isFinite(aluminum_phosphate_band_mol) or !std.math.isFinite(iron_phosphate_nonband_mol) or !std.math.isFinite(iron_phosphate_band_mol) or !std.math.isFinite(dicalcium_phosphate_nonband_mol) or !std.math.isFinite(dicalcium_phosphate_band_mol) or !std.math.isFinite(hydroxyapatite_nonband_mol) or !std.math.isFinite(hydroxyapatite_band_mol) or !std.math.isFinite(monocalcium_phosphate_nonband_mol) or !std.math.isFinite(monocalcium_phosphate_band_mol))
                    return error.InvalidHourlySoilChemistryGeometry;
                const phosphate_nonband_mol = non_band.dissolved_po4_mol_p_per_m3 * water_volume_m3;
                const phosphate_band_mol = band.dissolved_po4_mol_p_per_m3 * water_volume_m3;
                const phosphoric_acid_nonband_mol = non_band.dissolved_h3po4_mol_p_per_m3 * water_volume_m3;
                const phosphoric_acid_band_mol = band.dissolved_h3po4_mol_p_per_m3 * water_volume_m3;
                const iron_hpo4_nonband_mol = non_band.iron_hpo4_pair_mol_per_m3 * water_volume_m3;
                const iron_hpo4_band_mol = band.iron_hpo4_pair_mol_per_m3 * water_volume_m3;
                const iron_h2po4_nonband_mol = non_band.iron_h2po4_pair_mol_per_m3 * water_volume_m3;
                const iron_h2po4_band_mol = band.iron_h2po4_pair_mol_per_m3 * water_volume_m3;
                const calcium_hpo4_nonband_mol = non_band.calcium_po4_pair_mol_per_m3 * water_volume_m3;
                const calcium_hpo4_band_mol = band.calcium_po4_pair_mol_per_m3 * water_volume_m3;
                const calcium_h2po4_nonband_mol = non_band.calcium_h2po4_pair_mol_per_m3 * water_volume_m3;
                const calcium_h2po4_band_mol = band.calcium_h2po4_pair_mol_per_m3 * water_volume_m3;
                const calcium_phosphate_nonband_mol = non_band.calcium_po4_pair_mol_per_m3 * water_volume_m3;
                const calcium_phosphate_band_mol = band.calcium_po4_pair_mol_per_m3 * water_volume_m3;
                const magnesium_hpo4_nonband_mol = non_band.magnesium_hpo4_pair_mol_per_m3 * water_volume_m3;
                const magnesium_hpo4_band_mol = band.magnesium_hpo4_pair_mol_per_m3 * water_volume_m3;
                if (!std.math.isFinite(phosphate_nonband_mol) or !std.math.isFinite(phosphate_band_mol) or !std.math.isFinite(phosphoric_acid_nonband_mol) or !std.math.isFinite(phosphoric_acid_band_mol) or !std.math.isFinite(iron_hpo4_nonband_mol) or !std.math.isFinite(iron_hpo4_band_mol) or !std.math.isFinite(iron_h2po4_nonband_mol) or !std.math.isFinite(iron_h2po4_band_mol) or !std.math.isFinite(calcium_hpo4_nonband_mol) or !std.math.isFinite(calcium_hpo4_band_mol) or !std.math.isFinite(calcium_h2po4_nonband_mol) or !std.math.isFinite(calcium_h2po4_band_mol) or !std.math.isFinite(calcium_phosphate_nonband_mol) or !std.math.isFinite(calcium_phosphate_band_mol) or !std.math.isFinite(magnesium_hpo4_nonband_mol) or !std.math.isFinite(magnesium_hpo4_band_mol))
                    return error.InvalidHourlySoilChemistryGeometry;
                nitrate_nonband_pools[local_layer] = nitrate_nonband_g_n;
                nitrate_band_pools[local_layer] = nitrate_band_g_n;
                nitrite_nonband_pools[local_layer] = nitrite_nonband_g_n;
                nitrite_band_pools[local_layer] = nitrite_band_g_n;
                fertilizer_nitrate_nonband_pools[local_layer] = fertilizer_nitrate_nonband_g_n;
                fertilizer_nitrate_band_pools[local_layer] = fertilizer_nitrate_band_g_n;
                hydrogen_phosphate_nonband_pools[local_layer] = hydrogen_phosphate_nonband_mol;
                hydrogen_phosphate_band_pools[local_layer] = hydrogen_phosphate_band_mol;
                dihydrogen_phosphate_nonband_pools[local_layer] = dihydrogen_phosphate_nonband_mol;
                dihydrogen_phosphate_band_pools[local_layer] = dihydrogen_phosphate_band_mol;
                adsorbed_oh0_nonband_pools[local_layer] = adsorbed_oh0_nonband_mol;
                adsorbed_oh0_band_pools[local_layer] = adsorbed_oh0_band_mol;
                adsorbed_oh1_nonband_pools[local_layer] = adsorbed_oh1_nonband_mol;
                adsorbed_oh1_band_pools[local_layer] = adsorbed_oh1_band_mol;
                adsorbed_oh2_nonband_pools[local_layer] = adsorbed_oh2_nonband_mol;
                adsorbed_oh2_band_pools[local_layer] = adsorbed_oh2_band_mol;
                adsorbed_hpo4_nonband_pools[local_layer] = adsorbed_hpo4_nonband_mol;
                adsorbed_hpo4_band_pools[local_layer] = adsorbed_hpo4_band_mol;
                adsorbed_h2po4_nonband_pools[local_layer] = adsorbed_h2po4_nonband_mol;
                adsorbed_h2po4_band_pools[local_layer] = adsorbed_h2po4_band_mol;
                aluminum_phosphate_nonband_pools[local_layer] = aluminum_phosphate_nonband_mol;
                aluminum_phosphate_band_pools[local_layer] = aluminum_phosphate_band_mol;
                iron_phosphate_nonband_pools[local_layer] = iron_phosphate_nonband_mol;
                iron_phosphate_band_pools[local_layer] = iron_phosphate_band_mol;
                dicalcium_phosphate_nonband_pools[local_layer] = dicalcium_phosphate_nonband_mol;
                dicalcium_phosphate_band_pools[local_layer] = dicalcium_phosphate_band_mol;
                hydroxyapatite_nonband_pools[local_layer] = hydroxyapatite_nonband_mol;
                hydroxyapatite_band_pools[local_layer] = hydroxyapatite_band_mol;
                monocalcium_phosphate_nonband_pools[local_layer] = monocalcium_phosphate_nonband_mol;
                monocalcium_phosphate_band_pools[local_layer] = monocalcium_phosphate_band_mol;
                phosphate_nonband_pools[local_layer] = phosphate_nonband_mol;
                phosphate_band_pools[local_layer] = phosphate_band_mol;
                phosphoric_acid_nonband_pools[local_layer] = phosphoric_acid_nonband_mol;
                phosphoric_acid_band_pools[local_layer] = phosphoric_acid_band_mol;
                iron_hpo4_nonband_pools[local_layer] = iron_hpo4_nonband_mol;
                iron_hpo4_band_pools[local_layer] = iron_hpo4_band_mol;
                iron_h2po4_nonband_pools[local_layer] = iron_h2po4_nonband_mol;
                iron_h2po4_band_pools[local_layer] = iron_h2po4_band_mol;
                calcium_hpo4_nonband_pools[local_layer] = calcium_hpo4_nonband_mol;
                calcium_hpo4_band_pools[local_layer] = calcium_hpo4_band_mol;
                calcium_h2po4_nonband_pools[local_layer] = calcium_h2po4_nonband_mol;
                calcium_h2po4_band_pools[local_layer] = calcium_h2po4_band_mol;
                calcium_phosphate_nonband_pools[local_layer] = calcium_phosphate_nonband_mol;
                calcium_phosphate_band_pools[local_layer] = calcium_phosphate_band_mol;
                magnesium_hpo4_nonband_pools[local_layer] = magnesium_hpo4_nonband_mol;
                magnesium_hpo4_band_pools[local_layer] = magnesium_hpo4_band_mol;
            }

            var nitrate_geometry_for_update = ecosys.fertilizer_band_nitrate_phosphate.BandGeometry{
                .total_depth_m = nitrate_geometry.lower_edge_depth_m,
                .penetration_front_depth_m = nitrate_geometry.upper_edge_depth_m,
                .layer_depth_m = nitrate_geometry.band_depth_m[0..active_layer_count],
                .layer_width_m = nitrate_geometry.band_width_m[0..active_layer_count],
                .nonband_volume_fraction = nitrate_geometry.non_band_volume_fraction[0..active_layer_count],
                .band_volume_fraction = nitrate_geometry.band_volume_fraction[0..active_layer_count],
                .nonband_fractional_change_per_timestep = nonband_fractional_change[0..active_layer_count],
            };
            var phosphate_geometry_for_update = ecosys.fertilizer_band_nitrate_phosphate.BandGeometry{
                .total_depth_m = phosphate_geometry.lower_edge_depth_m,
                .penetration_front_depth_m = phosphate_geometry.upper_edge_depth_m,
                .layer_depth_m = phosphate_geometry.band_depth_m[0..active_layer_count],
                .layer_width_m = phosphate_geometry.band_width_m[0..active_layer_count],
                .nonband_volume_fraction = phosphate_geometry.non_band_volume_fraction[0..active_layer_count],
                .band_volume_fraction = phosphate_geometry.band_volume_fraction[0..active_layer_count],
                .nonband_fractional_change_per_timestep = nonband_fractional_change[0..active_layer_count],
            };
            var nitrate_pools = ecosys.fertilizer_band_nitrate_phosphate.NitratePools{
                .nitrate_nonband_g_n = nitrate_nonband_pools,
                .nitrate_band_g_n = nitrate_band_pools,
                .nitrite_nonband_g_n = nitrite_nonband_pools,
                .nitrite_band_g_n = nitrite_band_pools,
                .fertilizer_nitrate_nonband_g_n = fertilizer_nitrate_nonband_pools,
                .fertilizer_nitrate_band_g_n = fertilizer_nitrate_band_pools,
            };
            var phosphate_pools = ecosys.fertilizer_band_nitrate_phosphate.PhosphatePools{
                .hydrogen_phosphate_nonband_mol = hydrogen_phosphate_nonband_pools,
                .hydrogen_phosphate_band_mol = hydrogen_phosphate_band_pools,
                .dihydrogen_phosphate_nonband_mol = dihydrogen_phosphate_nonband_pools,
                .dihydrogen_phosphate_band_mol = dihydrogen_phosphate_band_pools,
                .adsorbed_oh0_nonband_mol = adsorbed_oh0_nonband_pools,
                .adsorbed_oh0_band_mol = adsorbed_oh0_band_pools,
                .adsorbed_oh1_nonband_mol = adsorbed_oh1_nonband_pools,
                .adsorbed_oh1_band_mol = adsorbed_oh1_band_pools,
                .adsorbed_oh2_nonband_mol = adsorbed_oh2_nonband_pools,
                .adsorbed_oh2_band_mol = adsorbed_oh2_band_pools,
                .adsorbed_hpo4_nonband_mol = adsorbed_hpo4_nonband_pools,
                .adsorbed_hpo4_band_mol = adsorbed_hpo4_band_pools,
                .adsorbed_h2po4_nonband_mol = adsorbed_h2po4_nonband_pools,
                .adsorbed_h2po4_band_mol = adsorbed_h2po4_band_pools,
                .aluminum_phosphate_nonband_mol = aluminum_phosphate_nonband_pools,
                .aluminum_phosphate_band_mol = aluminum_phosphate_band_pools,
                .iron_phosphate_nonband_mol = iron_phosphate_nonband_pools,
                .iron_phosphate_band_mol = iron_phosphate_band_pools,
                .dicalcium_phosphate_nonband_mol = dicalcium_phosphate_nonband_pools,
                .dicalcium_phosphate_band_mol = dicalcium_phosphate_band_pools,
                .hydroxyapatite_nonband_mol = hydroxyapatite_nonband_pools,
                .hydroxyapatite_band_mol = hydroxyapatite_band_pools,
                .monocalcium_phosphate_nonband_mol = monocalcium_phosphate_nonband_pools,
                .monocalcium_phosphate_band_mol = monocalcium_phosphate_band_pools,
                .phosphate_nonband_mol = phosphate_nonband_pools,
                .phosphate_band_mol = phosphate_band_pools,
                .phosphoric_acid_nonband_mol = phosphoric_acid_nonband_pools,
                .phosphoric_acid_band_mol = phosphoric_acid_band_pools,
                .iron_hpo4_nonband_mol = iron_hpo4_nonband_pools,
                .iron_hpo4_band_mol = iron_hpo4_band_pools,
                .iron_h2po4_nonband_mol = iron_h2po4_nonband_pools,
                .iron_h2po4_band_mol = iron_h2po4_band_pools,
                .calcium_hpo4_nonband_mol = calcium_hpo4_nonband_pools,
                .calcium_hpo4_band_mol = calcium_hpo4_band_pools,
                .calcium_h2po4_nonband_mol = calcium_h2po4_nonband_pools,
                .calcium_h2po4_band_mol = calcium_h2po4_band_pools,
                .calcium_phosphate_nonband_mol = calcium_phosphate_nonband_pools,
                .calcium_phosphate_band_mol = calcium_phosphate_band_pools,
                .magnesium_hpo4_nonband_mol = magnesium_hpo4_nonband_pools,
                .magnesium_hpo4_band_mol = magnesium_hpo4_band_pools,
            };
            for (0..active_layer_count) |layer_within_cell| {
                try ecosys.fertilizer_band_nitrate_phosphate.updateLayer(
                    .{
                        .first_active_layer_index = 0,
                        .layer_index = layer_within_cell,
                        .solute_timestep_h = 1,
                        .depth_threshold_m = 0,
                        .minimum_layer_thickness_m = context.runscript.soil_geometry_parameters.minimum_layer_thickness_m,
                        .absent_band_fraction_threshold = 1.0e-12,
                        .maximum_band_volume_fraction = 0.9999,
                        .nitrate_application = nitrate_application,
                        .nitrate_row_width_m = nitrate_row_width_m,
                        .phosphate_application = phosphate_application,
                        .phosphate_row_width_m = phosphate_row_width_m,
                        .salinity_chemistry = if (context.runscript.dynamic_plant_salts)
                            .enabled
                        else
                            .disabled,
                        .layers = layer_properties,
                    },
                    &nitrate_geometry_for_update,
                    &nitrate_pools,
                    &phosphate_geometry_for_update,
                    &phosphate_pools,
                );
            }
        }
    }
    try ecosys.soil_aqueous_transport_bridge.validateCarrierVolumes(
        context.micropore_solute_state,
        context.grid.matrix_liquid_water_m3,
        context.config.absolute_tolerance,
    );
    try ecosys.soil_aqueous_transport_bridge.exportChemistry(context.soil_chemistry, context.micropore_solute_state);
    try ecosys.fertilizer_band_production.consumeUndissolved(
        context.allocator,
        context.fertilizer_band,
        fertilizer_band_hour,
        context.soil_fertilizer_inventory,
        context.mineral_fertilizer_inventory,
        context.soil_chemistry,
        context.soil_reactive_nitrogen,
        context.runscript.dynamic_plant_salts,
    );
}

fn diagnosticPositiveRootNitrogenUptake_g(publication: anytype) f64 {
    var total_g_n: f64 = 0;
    inline for ([_]usize{ 0, 1, 4, 5 }) |pool| {
        for (publication.uptake_g_element_per_h_by_nutrient_and_layer[pool]) |amount_g_n| {
            total_g_n += @max(0, amount_g_n);
        }
    }
    return total_g_n;
}
