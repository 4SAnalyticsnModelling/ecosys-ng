const std = @import("std");
const ecosys = @import("ecosys_ng");

/// Dispatches cell-local science across either the resident domain or the
/// non-contiguous Morton-owned cells of one loaded tile. Complete vertical
/// columns remain indivisible worker units.
fn runScienceCellLayers(
    context: anytype,
    kernel_context: anytype,
    comptime kernel: anytype,
) !void {
    const cells = context.active_tile_cells.* orelse
        return error.MissingActiveScienceTile;
    try context.executor.runIndexedCellLayers(
        cells,
        context.grid.soil_layer_capacity,
        kernel_context,
        kernel,
    );
}

fn runScienceCells(
    context: anytype,
    kernel_context: anytype,
    comptime kernel: anytype,
) !void {
    const cells = context.active_tile_cells.* orelse
        return error.MissingActiveScienceTile;
    try context.executor.runIndexedCells(cells, kernel_context, kernel);
}

/// Preserves an hourly stage boundary while executing that stage across
/// serial Morton tiles. No owned-cell list is allocated here: TilePlan builds
/// and validates the lists once during runtime initialization.
fn runKernelAcrossSerialTiles(
    context: anytype,
    kernel_context: anytype,
    comptime kernel: anytype,
) !void {
    const plan = context.tile_plan;
    for (plan.tiles, 0..) |_, tile_index| {
        try context.executor.runIndexedCells(
            try plan.ownedCells(tile_index),
            kernel_context,
            kernel,
        );
    }
}

/// Initialization and refresh kernels obey the same execution boundary as the
/// hourly science: tiles advance serially, while owned cells within the active
/// tile may run concurrently.
fn runKernelAcrossSerialTilePlan(
    executor: ecosys.compute.CpuExecutor,
    plan: *const ecosys.spatial_grid.TilePlan,
    kernel_context: anytype,
    comptime kernel: anytype,
) !void {
    for (plan.tiles, 0..) |_, tile_index| {
        try executor.runIndexedCells(
            try plan.ownedCells(tile_index),
            kernel_context,
            kernel,
        );
    }
}

const SurfaceDenitrificationRespirationMap = struct {
    destination_g_c: []f64,
    source_g_c: []const f64,
};

fn mapSurfaceDenitrificationRespiration(
    context: *SurfaceDenitrificationRespirationMap,
    cells: ecosys.compute.CellRange,
) !void {
    for (cells.first..cells.end) |cell| {
        for (0..ecosys.surface_microbial_respiration_step.litter_complex_count) |complex| {
            const compact =
                cell * ecosys.surface_microbial_respiration_step.litter_complex_count +
                complex;
            const unit =
                cell * ecosys.surface_microbial_respiration_step.unit_count_per_cell +
                complex *
                    ecosys.surface_microbial_respiration_step.source_population_count +
                ecosys.surface_denitrification_step.denitrifier_population;
            context.destination_g_c[unit] = context.source_g_c[compact];
        }
    }
}

fn convergeHourlySoilChemistry(
    context: anytype,
    fertilizer_band_hour: ecosys.fertilizer_band_phase_coordinator.HourToken,
    failure_report: ?ecosys.solute_failure_reporter.Request,
) !void {
    const nutrient_zones = context.runscript.plant_nutrient_initialization;
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
    if (context.chemistry_reaction_parameters.* != null) {
        const chemistry_parameters = context.chemistry_reaction_parameters.*.?;
        const nitrogen_parameters =
            context.soil_nitrogen_parameters.* orelse
            return error.MissingSoilNitrogenParameters;
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
                const soil_mass_Mg = context.soil_solver_properties.matrix_bulk_volume_m3[layer] * context.soil_solver_properties.bulk_density_megagrams_per_m3[layer];
                const bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3[layer];
                if (!std.math.isFinite(water_volume_m3) or water_volume_m3 < 0 or !std.math.isFinite(soil_mass_Mg) or soil_mass_Mg <= 0) return error.InvalidHourlySoilChemistryGeometry;
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
                        .soil_mass_Mg = soil_mass_Mg,
                        .water_volume_m3 = water_volume_m3,
                        .biologically_active_water_volume_m3 = volq_m3,
                        .total_microbial_respiration_activity_g_c_per_step = toqck_g_c_per_step,
                        .temperature_response = tfnq,
                        .initial_inhibitor_activity = context.soil_fertilizer_inventory.initial_urease_inhibition_fraction[layer],
                        .current_inhibitor_activity = context.soil_fertilizer_inventory.current_urease_inhibition_fraction[layer],
                        .timestep_h = 1,
                    },
                    .{
                        .minimum_half_saturation_mol_n_per_Mg = fertilizer_parameters.minimum_urea_half_saturation_mol_n_per_Mg,
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
                const soil_mass_Mg = context.soil_solver_properties.matrix_bulk_volume_m3[layer] * context.soil_solver_properties.bulk_density_megagrams_per_m3[layer];
                if (!std.math.isFinite(water_volume_m3) or water_volume_m3 < 0 or !std.math.isFinite(soil_mass_Mg) or soil_mass_Mg <= 0) return error.InvalidHourlySoilChemistryGeometry;
                if (water_volume_m3 <= context.config.absolute_tolerance) continue;
                var parameters = context.soil_chemistry_layer_parameters[layer];
                const shared_ratio = soil_mass_Mg / water_volume_m3;
                parameters.cation_exchange_water_ratios = .{
                    .shared_Mg_per_m3 = shared_ratio,
                    .ammonium_non_band_Mg_per_m3 = if (parameters.fractions.ammonium_non_band > 0) shared_ratio / parameters.fractions.ammonium_non_band else shared_ratio,
                    .ammonium_band_Mg_per_m3 = if (parameters.fractions.ammonium_band > 0) shared_ratio / parameters.fractions.ammonium_band else shared_ratio,
                };
                parameters.non_band_phosphate_soil_mass_per_water_volume_Mg_per_m3 = if (parameters.fractions.phosphate_non_band > 0) shared_ratio / parameters.fractions.phosphate_non_band else shared_ratio;
                parameters.band_phosphate_soil_mass_per_water_volume_Mg_per_m3 = if (parameters.fractions.phosphate_band > 0) shared_ratio / parameters.fractions.phosphate_band else shared_ratio;
                parameters.total_carboxyl_sites_mol_per_Mg = context.chemistry_reaction_parameters.*.?.surface_litter.carboxyl_sites_mol_per_Mg_c *
                    1.0e-6 * context.soil_solver_properties.total_organic_carbon_g_per_megagram[layer];
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
                ) catch |err| {
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
                };
            }
        }
    }
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

fn applyPlantStorageRemobilization(context: anytype, calendar: ecosys.plant_development.Calendar) !void {
    if (context.plant_roots.* == null or context.detailed_canopy.* == null or context.plant_growth_stages.* == null or context.branch_development.* == null or context.plant_phenology.* == null or context.plant_dormancy.* == null or context.plant_storage_remobilization_workspace.* == null) return;
    const roots = &context.plant_roots.*.?;
    const canopy = &context.detailed_canopy.*.?;
    const stages = &context.plant_growth_stages.*.?;
    const development = &context.branch_development.*.?;
    const dormancy = &context.plant_dormancy.*.?;
    const phenology = &context.plant_phenology.*.?;
    const workspace = &context.plant_storage_remobilization_workspace.*.?;
    for (0..stages.plant_count) |plant| {
        const range = try stages.branchRange(plant);
        if (range.first == range.end or !context.plant_phenology.*.?.active[plant]) continue;
        const branch = (try stages.mainLivingBranch(plant)) orelse continue;
        const plant_parameters = context.shoot_growth_plant_parameters[plant];
        const dormancy_parameters = context.development_dormancy_parameters[plant];
        const dormancy_state = dormancy.branches[branch];
        const leafoff_start_fraction = if (plant_parameters.leaf_phenology_type == 0)
            dormancy_parameters.evergreen_leafoff_remobilization_start_fraction
        else
            dormancy_parameters.deciduous_leafoff_remobilization_start_fraction;
        if (!try ecosys.plant_storage_remobilization.activationEnabled(.{
            .annual_growth_habit = plant_parameters.growth_habit == 0,
            .lifecycle_initialized = phenology.lifecycle_initialized[plant],
            .current_day_of_year = calendar.day_of_year,
            .current_year = calendar.current_year,
            .planting_day_of_year = context.development_planting_day_of_year[plant],
            .planting_year = context.development_planting_year[plant],
            .accumulated_leafout_h = dormancy_state.accumulated_leafout_h,
            .required_leafout_h = dormancy_parameters.required_leafout_h,
            .accumulated_leafoff_h = dormancy_state.accumulated_leafoff_h,
            .required_leafoff_h = dormancy_parameters.required_leafoff_h,
            .leafoff_remobilization_start_fraction = leafoff_start_fraction,
        })) continue;
        const water = try ecosys.canopy_photosynthesis.canopyWaterGrowthResponse(
            plant_parameters.root_profile_type == 0,
            canopy.plant_canopy_turgor_potential_mpa[plant],
            context.runscript.phenology_parameters.minimum_turgor_potential_mpa,
            context.plants.canopy_water_potential_mpa[plant],
            plant_parameters.stomatal_turgor_shape_per_mpa,
        );
        const temperature_factor = canopy.plant_uptake_growth_temperature_response[plant];
        const remobilization_increment_h = try ecosys.plant_storage_remobilization.remobilizationTimeIncrementH(
            temperature_factor,
            water.growth_fraction,
            1,
        );
        if (!development.leafout_initialization_enabled[branch]) {
            development.remobilization_progress_h[branch] = 0;
            development.leafout_initialization_enabled[branch] = true;
        }
        development.remobilization_progress_h[branch] += remobilization_increment_h;

        var root_mobile_carbon_g_c: f64 = 0;
        var root_mobile_nitrogen_g_n: f64 = 0;
        var root_mobile_phosphorus_g_p: f64 = 0;
        for (0..roots.soil_layer_count) |layer| {
            const root = try roots.layerIndex(plant, 0, layer);
            root_mobile_carbon_g_c += roots.mobile_carbon_g[root];
            root_mobile_nitrogen_g_n += roots.mobile_nitrogen_g[root];
            root_mobile_phosphorus_g_p += roots.mobile_phosphorus_g[root];
        }
        const slices = try workspace.refreshPlant(roots, plant);
        const growth_habit: u8 = if (context.plant_topology_controls.growth_habit_code[plant] == 0) 0 else 1;
        const transfers = try ecosys.plant_storage_remobilization.calculate(context.runscript.storage_remobilization_parameters, .{
            .growth_habit = growth_habit,
            .aboveground_turnover_type = context.canopy_layer_controls.biomass_turnover_type[plant],
            .accumulated_remobilization_h = development.remobilization_progress_h[branch],
            .remobilization_time_increment_h = remobilization_increment_h,
            .biological_timestep_h = 1,
            .storage_carbon_g_c = canopy.plant_seed_storage_carbon_g[plant],
            .storage_nitrogen_g_n = canopy.plant_seed_storage_nitrogen_g[plant],
            .storage_phosphorus_g_p = canopy.plant_seed_storage_phosphorus_g[plant],
            .shoot_mobile_carbon_g_c = canopy.branch_mobile_carbon_g[branch],
            .shoot_mobile_nitrogen_g_n = canopy.branch_mobile_nitrogen_g[branch],
            .shoot_mobile_phosphorus_g_p = canopy.branch_mobile_phosphorus_g[branch],
            .root_mobile_carbon_g_c = root_mobile_carbon_g_c,
            .root_mobile_nitrogen_g_n = root_mobile_nitrogen_g_n,
            .root_mobile_phosphorus_g_p = root_mobile_phosphorus_g_p,
            .continue_annual_remobilization_after_duration = plant_parameters.leaf_phenology_type < 2,
        });
        try ecosys.plant_storage_remobilization.commit(canopy, roots, plant, branch, slices.root_layer_indices, slices.structural_carbon_fractions, slices.mobile_carbon_fractions, transfers);
        const remobilization_duration_h = context.runscript.storage_remobilization_parameters.remobilization_duration_h[growth_habit];
        for (range.first..range.end) |lateral_branch| {
            if (lateral_branch == branch or stages.branches[lateral_branch].dead or development.remobilization_progress_h[lateral_branch] > remobilization_duration_h) continue;
            development.remobilization_progress_h[lateral_branch] += remobilization_increment_h;
            try ecosys.shoot_growth_runtime.redistributeMainBranchMobileDuringRemobilization(
                canopy,
                branch,
                lateral_branch,
                temperature_factor,
                context.runscript.branch_mobile_exchange_parameters,
                1,
            );
        }
    }
}

fn applyStorageExhaustionMortality(context: anytype, calendar: ecosys.plant_development.Calendar) !void {
    if (context.plant_roots.* == null or context.detailed_canopy.* == null or context.plant_growth_stages.* == null or context.branch_development.* == null or context.plant_phenology.* == null or context.canopy_precipitation_retention.* == null) return;
    const canopy = &context.detailed_canopy.*.?;
    const phenology = &context.plant_phenology.*.?;
    for (0..phenology.active.len) |plant| {
        if (!phenology.active[plant]) continue;
        const cell = plant / context.config.plant_populations;
        const perennial = context.plant_topology_controls.growth_habit_code[plant] != 0;
        const storage_exhausted = try ecosys.plant_mortality.sourceOrderStorageExhausted(
            canopy.plant_seed_storage_carbon_g[plant],
            canopy.plant_seed_storage_nitrogen_g[plant],
            canopy.plant_seed_storage_phosphorus_g[plant],
            perennial,
            context.runscript.plant_pool_parameters.plant_root_structural_presence_g_per_plant,
            canopy.plant_population_count[plant],
        );
        const harvest = if (perennial and storage_exhausted)
            context.plant_harvest orelse return error.StorageExhaustionRequiresRootLitterContext
        else
            context.plant_harvest;
        const died = try ecosys.plant_mortality.applyStorageExhaustion(
            canopy,
            &context.plant_growth_stages.*.?,
            &context.branch_development.*.?,
            phenology,
            &context.plant_roots.*.?,
            &context.canopy_precipitation_retention.*.?,
            context.plants,
            plant,
            perennial,
            context.canopy_cell_area_m2[cell],
            context.runscript.plant_pool_parameters.plant_root_structural_presence_g_per_plant,
        );
        if (died) {
            try ecosys.plant_harvest_runtime.applyWholePlantMortalityResidue(harvest.?, plant);
            try ecosys.plant_harvest_runtime.releaseDeadRootsToLitter(harvest.?, plant);
            const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(harvest.?, plant);
            if (exported.carbon_g != 0 or exported.nitrogen_g != 0 or exported.phosphorus_g != 0)
                return error.UnexpectedPlantMortalityExport;
            if (calendar.current_year < 0 or calendar.current_year > std.math.maxInt(u16))
                return error.PerennialReplantYearOutOfRange;
            const current_year: u16 = @intCast(calendar.current_year);
            try ecosys.plant_mortality.scheduleNextDayReplant(
                phenology,
                plant,
                calendar.day_of_year,
                current_year,
                ecosys.climate_change.daysInYear(current_year),
            );
        }
        if (died) std.log.info(
            "perennial storage exhausted: plant={d} cell={d} storage_C_g={e} storage_N_g={e} storage_P_g={e}",
            .{ plant, cell, canopy.plant_seed_storage_carbon_g[plant], canopy.plant_seed_storage_nitrogen_g[plant], canopy.plant_seed_storage_phosphorus_g[plant] },
        );
    }
}

fn applyNaturalBranchMortality(context: anytype, calendar: ecosys.plant_development.Calendar) !void {
    if (context.detailed_canopy.* == null or context.plant_growth_stages.* == null or context.branch_development.* == null or
        context.plant_phenology.* == null or context.plant_roots.* == null or context.canopy_precipitation_retention.* == null)
        return;
    const harvest = context.plant_harvest orelse return error.NaturalBranchDeathRequiresLitterContext;
    const canopy = &context.detailed_canopy.*.?;
    const growth = &context.plant_growth_stages.*.?;
    const development = &context.branch_development.*.?;
    const phenology = &context.plant_phenology.*.?;
    for (0..phenology.active.len) |plant| {
        if (!phenology.active[plant]) continue;
        const parameters = context.shoot_growth_plant_parameters[plant];
        const branches = try growth.branchRange(plant);
        var dead_count: usize = 0;
        var encountered_dead = false;
        for (branches.first..branches.end) |branch| {
            const dead = growth.branches[branch].dead or development.dead[branch];
            if (!dead) continue;
            dead_count += 1;
            encountered_dead = true;
            try ecosys.plant_harvest_runtime.applyNaturalDeadBranchResidue(
                harvest,
                plant,
                branch,
                parameters.growth_habit == 0 and parameters.leaf_phenology_type != 0,
            );
            try ecosys.shoot_growth_runtime.resetNaturalDeadBranch(growth, development, &context.plant_dormancy.*.?, branch);
        }
        if (dead_count == branches.end - branches.first) {
            const winter_annual = parameters.growth_habit == 0 and parameters.leaf_phenology_type != 0;
            if (!winter_annual)
                try ecosys.plant_harvest_runtime.applyWholePlantMortalityResidue(harvest, plant);
            const cell = plant / context.config.plant_populations;
            try ecosys.plant_mortality.applyAllBranchesDead(
                canopy,
                phenology,
                &context.plant_roots.*.?,
                &context.canopy_precipitation_retention.*.?,
                context.plants,
                plant,
                winter_annual,
                context.canopy_cell_area_m2[cell],
            );
            try ecosys.plant_harvest_runtime.releaseDeadRootsToLitter(harvest, plant);
            if (parameters.growth_habit != 0) {
                if (calendar.current_year < 0 or calendar.current_year > std.math.maxInt(u16))
                    return error.PerennialReplantYearOutOfRange;
                const year: u16 = @intCast(calendar.current_year);
                try ecosys.plant_mortality.scheduleNextDayReplant(phenology, plant, calendar.day_of_year, year, ecosys.climate_change.daysInYear(year));
            }
        }
        if (encountered_dead) {
            const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(harvest, plant);
            if (exported.carbon_g != 0 or exported.nitrogen_g != 0 or exported.phosphorus_g != 0)
                return error.UnexpectedNaturalBranchDeathExport;
        }
    }
}

fn previousRootNutrientDemand(roots: anytype, pool: ecosys.plant_root_nutrient_uptake.NutrientPool, root: usize) f64 {
    return switch (pool) {
        .ammonium_nonband => roots.previous_ammonium_demand_nonband_g_n_per_h[root],
        .ammonium_band => roots.previous_ammonium_demand_band_g_n_per_h[root],
        .nitrate_nonband => roots.previous_nitrate_demand_nonband_g_n_per_h[root],
        .nitrate_band => roots.previous_nitrate_demand_band_g_n_per_h[root],
        .phosphate_h2_nonband => roots.previous_phosphate_h2_demand_nonband_g_p_per_h[root],
        .phosphate_h2_band => roots.previous_phosphate_h2_demand_band_g_p_per_h[root],
        .phosphate_h_nonband => roots.previous_phosphate_h_demand_nonband_g_p_per_h[root],
        .phosphate_h_band => roots.previous_phosphate_h_demand_band_g_p_per_h[root],
    };
}

const DailyPlantElement = enum { nitrogen, phosphorus };

fn calculateDailyPlantElementPools(
    canopy: *const ecosys.canopy_photosynthesis.State,
    roots: *const ecosys.plant_root_system.State,
    plant: usize,
    branches: ecosys.canopy_photosynthesis.Range,
    element: DailyPlantElement,
    biological_domain_count: usize,
    root_workspace_g_by_layer: []f64,
) !ecosys.plant_daily_pool_aggregation.ElementPools {
    if (root_workspace_g_by_layer.len != roots.soil_layer_count or biological_domain_count < 1 or biological_domain_count > ecosys.plant_root_system.biological_domain_count)
        return error.DailyPlantElementRootDimensionMismatch;
    @memset(root_workspace_g_by_layer, 0);
    var root_symbiont_g: f64 = 0;
    for (0..roots.soil_layer_count) |layer| for (0..biological_domain_count) |domain| {
        const root = try roots.layerIndex(plant, domain, layer);
        root_workspace_g_by_layer[layer] += switch (element) {
            .nitrogen => roots.mobile_nitrogen_g[root],
            .phosphorus => roots.mobile_phosphorus_g[root],
        };
        root_symbiont_g += switch (element) {
            .nitrogen => roots.symbiont_structural_nitrogen_g_n[root] + roots.symbiont_mobile_nitrogen_g_n[root],
            .phosphorus => roots.symbiont_structural_phosphorus_g_p[root] + roots.symbiont_mobile_phosphorus_g_p[root],
        };
        for (0..roots.active_root_axis_count[plant]) |axis| {
            const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
            root_workspace_g_by_layer[layer] += switch (element) {
                .nitrogen => roots.axis_primary_nitrogen_g[axis_layer] + roots.axis_secondary_nitrogen_g[axis_layer],
                .phosphorus => roots.axis_primary_phosphorus_g[axis_layer] + roots.axis_secondary_phosphorus_g[axis_layer],
            };
        }
    };
    var canopy_symbiont_g: f64 = 0;
    for (branches.first..branches.end) |branch| {
        canopy_symbiont_g += switch (element) {
            .nitrogen => canopy.branch_symbiont_structural_nitrogen_g[branch] + canopy.branch_symbiont_mobile_nitrogen_g[branch],
            .phosphorus => canopy.branch_symbiont_structural_phosphorus_g[branch] + canopy.branch_symbiont_mobile_phosphorus_g[branch],
        };
    }
    return ecosys.plant_daily_pool_aggregation.calculateElement(.{
        .branch_leaf_g = switch (element) {
            .nitrogen => canopy.branch_leaf_nitrogen_g[branches.first..branches.end],
            .phosphorus => canopy.branch_leaf_phosphorus_g[branches.first..branches.end],
        },
        .branch_sheath_g = switch (element) {
            .nitrogen => canopy.branch_sheath_nitrogen_g[branches.first..branches.end],
            .phosphorus => canopy.branch_sheath_phosphorus_g[branches.first..branches.end],
        },
        .branch_stalk_g = switch (element) {
            .nitrogen => canopy.branch_stalk_nitrogen_g[branches.first..branches.end],
            .phosphorus => canopy.branch_stalk_phosphorus_g[branches.first..branches.end],
        },
        .branch_reserve_g = switch (element) {
            .nitrogen => canopy.branch_reserve_nitrogen_g[branches.first..branches.end],
            .phosphorus => canopy.branch_reserve_phosphorus_g[branches.first..branches.end],
        },
        .branch_husk_g = switch (element) {
            .nitrogen => canopy.branch_husk_nitrogen_g[branches.first..branches.end],
            .phosphorus => canopy.branch_husk_phosphorus_g[branches.first..branches.end],
        },
        .branch_ear_g = switch (element) {
            .nitrogen => canopy.branch_ear_nitrogen_g[branches.first..branches.end],
            .phosphorus => canopy.branch_ear_phosphorus_g[branches.first..branches.end],
        },
        .branch_grain_g = switch (element) {
            .nitrogen => canopy.branch_grain_nitrogen_g[branches.first..branches.end],
            .phosphorus => canopy.branch_grain_phosphorus_g[branches.first..branches.end],
        },
        .branch_mobile_g = switch (element) {
            .nitrogen => canopy.branch_mobile_nitrogen_g[branches.first..branches.end],
            .phosphorus => canopy.branch_mobile_phosphorus_g[branches.first..branches.end],
        },
        .root_g_by_layer = root_workspace_g_by_layer,
        .canopy_symbiont_g = canopy_symbiont_g,
        .root_symbiont_g = root_symbiont_g,
        .standing_dead_g = switch (element) {
            .nitrogen => canopy.plant_standing_dead_nitrogen_g[plant],
            .phosphorus => canopy.plant_standing_dead_phosphorus_g[plant],
        },
        .seed_storage_g = switch (element) {
            .nitrogen => canopy.plant_seed_storage_nitrogen_g[plant],
            .phosphorus => canopy.plant_seed_storage_phosphorus_g[plant],
        },
    });
}

fn applyRootMetabolism(context: anytype) !void {
    if (context.plant_roots.* == null or context.detailed_canopy.* == null or context.plant_growth_stages.* == null or context.branch_development.* == null or context.plant_dormancy.* == null or context.plant_root_metabolism_workspace.* == null or context.plant_storage_remobilization_workspace.* == null or context.plant_water_workspace.* == null or context.plant_litter_partition.* == null) return;
    const roots = &context.plant_roots.*.?;
    try context.root_litter_carbon_ledger.validateDimensions(
        roots.plant_count,
        ecosys.plant_root_system.biological_domain_count,
        roots.soil_layer_count,
    );
    @memset(context.root_litter_products_by_plant, std.mem.zeroes(ecosys.plant_root_metabolism.RootLitter));
    context.root_litter_carbon_ledger.resetHourly();
    const canopy = &context.detailed_canopy.*.?;
    const growth_stages = &context.plant_growth_stages.*.?;
    const branch_development = &context.branch_development.*.?;
    const parameters = context.runscript.root_metabolism_parameters;
    for (0..roots.plant_count) |plant| {
        if (!context.plant_phenology.*.?.active[plant]) continue;
        var oxygen_uptake_g_o: f64 = 0;
        var oxygen_demand_g_o: f64 = 0;
        for (0..context.root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| for (0..roots.soil_layer_count) |layer| {
            const root = try roots.layerIndex(plant, domain, layer);
            oxygen_uptake_g_o += roots.oxygen_uptake_g_o_per_h[root];
            oxygen_demand_g_o += roots.oxygen_demand_g_o_per_h[root];
        };
        const oxygen_satisfaction = try ecosys.plant_root_porosity.sourceOxygenSatisfaction(
            oxygen_uptake_g_o,
            oxygen_demand_g_o,
            context.config.absolute_tolerance,
        );
        roots.current_porosity_fraction[plant] = try ecosys.plant_root_porosity.adapt(
            roots.current_porosity_fraction[plant],
            roots.initial_porosity_fraction[plant],
            oxygen_satisfaction,
            1,
            context.runscript.root_porosity_parameters,
        );
        context.plant_water_workspace.*.?.root_porosity_fraction[plant] = roots.current_porosity_fraction[plant];
    }
    for (0..context.grid.cell_count) |cell| {
        var workspace = &context.plant_root_metabolism_workspace.*.?.per_cell[cell];
        for (0..context.config.plant_populations) |species| {
            const plant = cell * context.config.plant_populations + species;
            if (!context.plant_phenology.*.?.active[plant] or roots.roots_dead[plant]) continue;
            const traits = context.root_metabolism_plant_parameters[plant];
            try workspace.beginPlantHour(roots.active_root_axis_count[plant]);
            const branch_range = try growth_stages.branchRange(plant);
            if (branch_range.first == branch_range.end) continue;
            const main_branch = (try growth_stages.mainLivingBranch(plant)) orelse continue;
            const root_respiration_active = ecosys.plant_root_metabolism.rootRespirationActive(
                traits.growth_habit != 0,
                branch_development.stage_day[main_branch * 10 + 9] != 0,
            );
            try ecosys.plant_dormancy.advanceRemobilization(
                &context.plant_dormancy.*.?.branches[main_branch],
                .{
                    .timestep_h = 1,
                    .canopy_temperature_c = context.plants.canopy_temperature_k[plant] - 273.15,
                    .canopy_total_water_potential_mpa = context.plants.canopy_water_potential_mpa[plant],
                },
                context.development_dormancy_parameters[plant],
                ecosys.plant_growth_stages.growthHabitFromReadq(traits.growth_habit),
                try ecosys.plant_growth_stages.phenologyTypeFromReadq(traits.leaf_phenology_type),
                growth_stages.branches[main_branch].seed_number_set_end_day != 0,
            );
            var stalk_carbon_g_c: f64 = 0;
            var sapwood_carbon_g_c: f64 = 0;
            for (branch_range.first..branch_range.end) |branch| {
                stalk_carbon_g_c += canopy.branch_stalk_carbon_g[branch];
                sapwood_carbon_g_c += canopy.branch_sapwood_carbon_g[branch];
            }
            const wood = try ecosys.plant_root_metabolism.rootWoodComposition(
                context.canopy_layer_controls.biomass_turnover_type[plant] != 0,
                traits.root_profile_type > 1,
                stalk_carbon_g_c,
                sapwood_carbon_g_c,
                traits.stalk_nitrogen_to_carbon_g_n_per_g_c,
                traits.root_nitrogen_to_carbon_g_n_per_g_c,
                traits.stalk_phosphorus_to_carbon_g_p_per_g_c,
                traits.root_phosphorus_to_carbon_g_p_per_g_c,
                context.runscript.plant_pool_parameters.plant_root_structural_presence_g_per_plant,
                parameters.nonwoody_root_fraction_exponent,
            );
            const termination_feedback = try ecosys.plant_root_metabolism.annualTerminationFeedback(
                traits.growth_habit,
                branch_development.hours_without_grain_fill[main_branch],
                parameters.annual_termination_hours_without_grain_fill,
            );
            const coarse_litter = try context.plant_litter_partition.*.?.get(plant, .coarse_wood);
            const fine_litter = try context.plant_litter_partition.*.?.get(plant, .fine_root);
            const litter_kinetics: ecosys.plant_root_metabolism.RootLitterFractions = .{
                .woody_carbon = coarse_litter.carbon,
                .woody_nitrogen = coarse_litter.nitrogen,
                .woody_phosphorus = coarse_litter.phosphorus,
                .nonwoody_carbon = fine_litter.carbon,
                .nonwoody_nitrogen = fine_litter.nitrogen,
                .nonwoody_phosphorus = fine_litter.phosphorus,
            };
            var layer_top_m: f64 = 0;
            for (0..context.grid.active_soil_layer_count[cell]) |layer| {
                const soil = try context.grid.layerIndex(cell, layer);
                const layer_bottom_m = layer_top_m + context.soil_solver_properties.layer_thickness_m[soil];
                for (0..traits.biologicalDomainCount()) |domain| {
                    const root = try roots.layerIndex(plant, domain, layer);
                    const active_axis_count = roots.active_root_axis_count[plant];
                    if (active_axis_count == 0) continue;
                    try workspace.resetAxes(active_axis_count);
                    var litterfall: ecosys.plant_root_litterfall.LayerInput = .{};
                    var root_carbon_g_c: f64 = 0;
                    for (0..active_axis_count) |axis| {
                        const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
                        root_carbon_g_c += roots.axis_primary_carbon_g[axis_layer] + roots.axis_secondary_carbon_g[axis_layer];
                    }
                    var upper_host_active_root_carbon_g_c: f64 = 0;
                    if (domain == 0 and layer > 0) for (0..active_axis_count) |axis| {
                        const upper_axis_layer = try roots.layerAxisIndex(plant, 0, layer - 1, axis);
                        upper_host_active_root_carbon_g_c +=
                            roots.axis_primary_carbon_g[upper_axis_layer] +
                            roots.axis_secondary_carbon_g[upper_axis_layer];
                    };
                    const mobile_c = roots.mobile_carbon_g[root];
                    const mobile_n = roots.mobile_nitrogen_g[root];
                    const mobile_p = roots.mobile_phosphorus_g[root];
                    const recycling = try ecosys.plant_root_metabolism.secondaryRootRecyclingFractions(
                        context.development_emerged[plant],
                        if (root_carbon_g_c > 0) mobile_c / root_carbon_g_c else 1,
                        if (root_carbon_g_c > 0) mobile_n / root_carbon_g_c else 1,
                        if (root_carbon_g_c > 0) mobile_p / root_carbon_g_c else 1,
                        parameters,
                    );
                    const hydrogen_mol_per_m3 = context.soil_chemistry.aqueous[soil].hydrogen;
                    const soil_ph = if (hydrogen_mol_per_m3 > 0) -@log10(hydrogen_mol_per_m3 / 1.0e3) else 7;
                    const growth_temperature = try ecosys.plant_root_nutrient_uptake.rootGrowthTemperatureResponse(context.grid.soil_temperature_k[soil], canopy.plant_thermal_adaptation_offset_c[plant], context.runscript.canopy_stress_parameters.growth_temperature);
                    const primary_environment = try ecosys.plant_root_metabolism.rootEnvironmentResponses(
                        parameters,
                        context.grid.soil_temperature_k[soil],
                        canopy.plant_thermal_adaptation_offset_c[plant],
                        soil_ph,
                        roots.total_water_potential_mpa[root],
                        roots.turgor_water_potential_mpa[root],
                        context.runscript.phenology_parameters.minimum_turgor_potential_mpa,
                        0,
                        traits.primary_root_radius_m,
                        traits.root_profile_type == 0,
                    );
                    const secondary_environment = try ecosys.plant_root_metabolism.rootEnvironmentResponses(
                        parameters,
                        context.grid.soil_temperature_k[soil],
                        canopy.plant_thermal_adaptation_offset_c[plant],
                        soil_ph,
                        roots.total_water_potential_mpa[root],
                        roots.turgor_water_potential_mpa[root],
                        context.runscript.phenology_parameters.minimum_turgor_potential_mpa,
                        0,
                        traits.secondary_root_radius_m,
                        traits.root_profile_type == 0,
                    );
                    for (0..active_axis_count) |axis| {
                        const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
                        const axis_index = try roots.axisIndex(plant, domain, axis);
                        const primary_depth_m = roots.axis_depth_m[axis_index];
                        const tip_in_layer = primary_depth_m > layer_top_m and (primary_depth_m <= layer_bottom_m or layer + 1 == context.grid.active_soil_layer_count[cell]);
                        const secondary_reaches_layer = primary_depth_m > layer_top_m;
                        workspace.primary_active[axis] = domain == 0 and tip_in_layer and !try workspace.primaryWasProcessed(domain, axis);
                        workspace.secondary_active[axis] = secondary_reaches_layer and
                            (roots.axis_secondary_count[axis_layer] > 0 or roots.axis_secondary_carbon_g[axis_layer] > context.config.absolute_tolerance);
                        workspace.sink_strengths[axis] = try ecosys.plant_root_metabolism.sourceOrderRootAxisSinkStrength(parameters, .{
                            .root_profile_type = traits.root_profile_type,
                            .primary_axis_count_multiplier = @max(1, roots.axis_primary_count[axis_layer]),
                            .primary_root_radius_m = traits.primary_root_radius_m,
                            .primary_root_depth_from_surface_m = primary_depth_m,
                            .layer_top_depth_m = layer_top_m,
                            .layer_thickness_m = context.soil_solver_properties.layer_thickness_m[soil],
                            .secondary_root_origin_offset_m = 0,
                            .seeding_depth_m = context.plant_water_workspace.*.?.seeding_depth_m[plant],
                            .hypocotyledon_height_m = canopy.plant_hypocotyledon_height_m[plant],
                            .canopy_height_m = context.development_canopy_height_m[plant],
                            .secondary_axis_count = roots.axis_secondary_count[axis_layer],
                            .secondary_root_radius_m = traits.secondary_root_radius_m,
                            .average_secondary_root_length_m = roots.average_secondary_length_m[root],
                            .negligible_sink_m = context.config.absolute_tolerance,
                            .primary_biological_domain = domain == 0,
                        });
                    }
                    roots.sink_strength_m[root] = try ecosys.plant_root_metabolism.normalizeRootAxisSinkFractions(
                        workspace.sink_strengths[0..active_axis_count],
                        workspace.primary_sink_fractions[0..active_axis_count],
                        workspace.secondary_sink_fractions[0..active_axis_count],
                        context.config.absolute_tolerance,
                    );
                    for (0..active_axis_count) |axis| {
                        const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
                        inline for (.{ true, false }) |primary| {
                            const active = if (primary) workspace.primary_active[axis] else workspace.secondary_active[axis];
                            if (active) {
                                const axis_c = if (primary) roots.axis_primary_carbon_g[axis_layer] else roots.axis_secondary_carbon_g[axis_layer];
                                const axis_n = if (primary) roots.axis_primary_nitrogen_g[axis_layer] else roots.axis_secondary_nitrogen_g[axis_layer];
                                const axis_p = if (primary) roots.axis_primary_phosphorus_g[axis_layer] else roots.axis_secondary_phosphorus_g[axis_layer];
                                const environment = if (primary) primary_environment else secondary_environment;
                                const water_responses = try ecosys.plant_root_metabolism.sourceRootRespirationWaterResponses(
                                    traits.root_profile_type,
                                    traits.leaf_phenology_type,
                                    environment.growth_water,
                                    environment.maintenance_water,
                                );
                                const shared: ecosys.plant_root_metabolism.SecondaryRootInputs = .{
                                    .mobile_carbon_g_c = mobile_c,
                                    .nonstructural_nitrogen_g_n = mobile_n,
                                    .nonstructural_phosphorus_g_p = mobile_p,
                                    .root_carbon_g_c = root_carbon_g_c,
                                    .root_nitrogen_g_n = axis_n,
                                    .root_nitrogen_to_carbon_ratio_g_n_per_g_c = wood.growth_nitrogen_to_carbon_g_n_per_g_c,
                                    .root_phosphorus_to_carbon_ratio_g_p_per_g_c = wood.growth_phosphorus_to_carbon_g_p_per_g_c,
                                    .root_growth_yield_g_c_per_g_c = traits.root_growth_yield_g_c_per_g_c,
                                    .active_root_fraction = if (primary) workspace.primary_sink_fractions[axis] else workspace.secondary_sink_fractions[axis],
                                    .biological_timestep_h = 1,
                                    .substrate_temperature_response = if (root_respiration_active) growth_temperature else 0,
                                    .maintenance_temperature_response = if (root_respiration_active) environment.maintenance_temperature else 0,
                                    .acidity_response = environment.acidity,
                                    .substrate_feedback = termination_feedback,
                                    .oxygen_limitation = std.math.clamp(roots.oxygen_process_constraint_fraction[root], 0, 1),
                                    .substrate_water_response = water_responses.substrate,
                                    .maintenance_water_response = water_responses.maintenance,
                                };
                                const metabolism = if (primary)
                                    try ecosys.plant_root_metabolism.primaryRootMetabolism(parameters, .{ .shared = shared, .primary_tip_at_or_below_profile_bottom = layer + 1 == context.grid.active_soil_layer_count[cell] })
                                else
                                    try ecosys.plant_root_metabolism.secondaryRootMetabolism(parameters, shared);
                                const senescence_inputs: ecosys.plant_root_metabolism.SecondaryRootSenescenceInputs = .{
                                    .oxygen_unlimited_substrate_minus_maintenance_g_c_per_h = metabolism.substrate_respiration_oxygen_unlimited_g_c_per_h - metabolism.maintenance_respiration_g_c_per_h,
                                    .actual_substrate_minus_maintenance_g_c_per_h = metabolism.substrate_respiration_actual_g_c_per_h - metabolism.maintenance_respiration_g_c_per_h,
                                    .root_carbon_g_c = axis_c,
                                    .root_nitrogen_g_n = axis_n,
                                    .root_phosphorus_g_p = axis_p,
                                    .oxygen_limitation = shared.oxygen_limitation,
                                    .phenological_remobilization_enabled = context.plant_dormancy.*.?.branches[main_branch].phenological_remobilization_enabled,
                                    .root_remobilization_enabled = context.plant_dormancy.*.?.branches[main_branch].shoot_remobilization_enabled,
                                    .storage_exchange_fraction_per_h = parameters.storage_exchange_fraction_per_h,
                                    .remobilization_elapsed_h = context.plant_dormancy.*.?.branches[main_branch].remobilization_elapsed_h,
                                    .full_senescence_h = parameters.full_senescence_duration_h,
                                    .biological_timestep_h = 1,
                                    .structural_presence_threshold_g_c = context.runscript.plant_pool_parameters.plant_root_structural_presence_g_per_plant,
                                };
                                const senescence = if (primary)
                                    try ecosys.plant_root_metabolism.primaryRootSenescence(senescence_inputs, recycling)
                                else
                                    try ecosys.plant_root_metabolism.secondaryRootSenescence(senescence_inputs, recycling);
                                try litterfall.add(try ecosys.plant_root_metabolism.secondaryRootLitter(
                                    senescence,
                                    axis_c,
                                    axis_n,
                                    axis_p,
                                    wood.carbon_fraction,
                                    wood.nitrogen_fraction,
                                    wood.phosphorus_fraction,
                                    litter_kinetics,
                                ));
                                if (primary) {
                                    workspace.primary_metabolism[axis] = metabolism;
                                    workspace.primary_senescence[axis] = senescence;
                                } else {
                                    workspace.secondary_metabolism[axis] = metabolism;
                                    workspace.secondary_senescence[axis] = senescence;
                                }
                            }
                        }
                    }
                    try ecosys.plant_root_litterfall.validatePublication(context.soil_organic, soil, litterfall);
                    try ecosys.plant_root_metabolism.commitStagedLayerAxes(roots, plant, domain, layer, workspace, active_axis_count, .{
                        .primary_specific_length_m_per_g_c = traits.primary_specific_length_m_per_g_c,
                        .secondary_specific_length_m_per_g_c = traits.secondary_specific_length_m_per_g_c,
                        .plant_population_count = context.plant_water_workspace.*.?.plant_population_count[plant],
                        .seeding_depth_m = context.plant_water_workspace.*.?.seeding_depth_m[plant],
                        .current_layer_bottom_depth_m = layer_bottom_m,
                        .next_layer_thickness_m = if (layer + 1 < context.grid.active_soil_layer_count[cell]) context.soil_solver_properties.layer_thickness_m[try context.grid.layerIndex(cell, layer + 1)] else 0,
                        .extension_presence_threshold_m = context.runscript.plant_pool_parameters.plant_root_structural_presence_g_per_plant *
                            context.plant_water_workspace.*.?.plant_population_count[plant],
                        .root_extension_water_response = @min(primary_environment.extension_water, secondary_environment.extension_water),
                        .nonwoody_carbon_fraction = wood.carbon_fraction[1],
                        .nonwoody_nitrogen_fraction = wood.nitrogen_fraction[1],
                        .nonwoody_phosphorus_fraction = wood.phosphorus_fraction[1],
                        .protein_carbon_per_nitrogen_g_c_per_g_n = parameters.root_protein_carbon_per_nitrogen_g_c_per_g_n,
                        .protein_carbon_per_phosphorus_g_c_per_g_p = parameters.root_protein_carbon_per_phosphorus_g_c_per_g_p,
                    });
                    const domain_litterfall = litterfall.litter;
                    var current_mycorrhizal_litterfall = std.mem.zeroes(ecosys.plant_root_metabolism.RootLitter);
                    var upper_mycorrhizal_litterfall: ecosys.plant_root_litterfall.LayerInput = .{};
                    if (domain == 0) {
                        const concurrent_loss = try ecosys.plant_root_metabolism.commitMycorrhizalLossWithSecondaryRoots(
                            roots,
                            plant,
                            layer,
                            workspace.*,
                            active_axis_count,
                            // GROSUB FSNCP uses WTRTL(1): total active host
                            // primary plus secondary root C, not WTRT2 alone.
                            .{ root_carbon_g_c, upper_host_active_root_carbon_g_c },
                            context.runscript.plant_pool_parameters.plant_root_structural_presence_g_per_plant,
                            .{ wood.carbon_fraction, wood.nitrogen_fraction, wood.phosphorus_fraction },
                            litter_kinetics,
                        );
                        current_mycorrhizal_litterfall = concurrent_loss.current;
                        try litterfall.add(concurrent_loss.current);
                        try upper_mycorrhizal_litterfall.add(concurrent_loss.upper);
                        try ecosys.plant_root_litterfall.validatePublication(context.soil_organic, soil, litterfall);
                        if (layer > 0)
                            try ecosys.plant_root_litterfall.validatePublication(context.soil_organic, try context.grid.layerIndex(cell, layer - 1), upper_mycorrhizal_litterfall);
                    }
                    var plant_litterfall: ecosys.plant_root_litterfall.LayerInput = .{ .litter = context.root_litter_products_by_plant[plant] };
                    try plant_litterfall.add(litterfall.litter);
                    try plant_litterfall.add(upper_mycorrhizal_litterfall.litter);
                    for (0..active_axis_count) |axis| if (workspace.primary_active[axis]) try workspace.markPrimaryProcessed(domain, axis);
                    ecosys.plant_root_litterfall.publishValidated(context.soil_organic, soil, litterfall);
                    if (domain == 0 and layer > 0)
                        ecosys.plant_root_litterfall.publishValidated(context.soil_organic, try context.grid.layerIndex(cell, layer - 1), upper_mycorrhizal_litterfall);
                    context.root_litter_carbon_ledger.addValidated(plant, domain, layer, domain_litterfall);
                    if (domain == 0) {
                        // Concurrent loss is mycorrhizal (legacy N=2), even
                        // though its host trigger is evaluated in domain zero.
                        context.root_litter_carbon_ledger.addValidated(plant, 1, layer, current_mycorrhizal_litterfall);
                        if (layer > 0)
                            context.root_litter_carbon_ledger.addValidated(plant, 1, layer - 1, upper_mycorrhizal_litterfall.litter);
                    }
                    context.root_litter_products_by_plant[plant] = plant_litterfall.litter;
                    for (0..active_axis_count) |axis| {
                        if (!workspace.primary_active[axis]) continue;
                        const metabolism = workspace.primary_metabolism[axis];
                        const senescence = workspace.primary_senescence[axis];
                        const respiration = try ecosys.plant_root_metabolism.assemble(.{
                            .maintenance_demand_g_c = metabolism.maintenance_respiration_g_c_per_h,
                            .substrate_respiration_actual_g_c = metabolism.substrate_respiration_actual_g_c_per_h,
                            .substrate_respiration_oxygen_unlimited_g_c = metabolism.substrate_respiration_oxygen_unlimited_g_c_per_h,
                            .growth_respiration_actual_g_c = metabolism.growth_respiration_actual_g_c_per_h,
                            .growth_respiration_oxygen_unlimited_g_c = metabolism.growth_respiration_oxygen_unlimited_g_c_per_h,
                            .senescence_respiration_actual_g_c = senescence.respiration_actual_g_c_per_h,
                            .senescence_respiration_oxygen_unlimited_g_c = senescence.respiration_oxygen_unlimited_g_c_per_h,
                            .nitrogen_assimilation_respiration_actual_g_c = metabolism.nitrogen_assimilation_respiration_actual_g_c_per_h,
                            .nitrogen_assimilation_respiration_oxygen_unlimited_g_c = metabolism.nitrogen_assimilation_respiration_oxygen_unlimited_g_c_per_h,
                        });
                        const traversed_layer_count = layer + 1;
                        for (0..traversed_layer_count) |traversed_layer| {
                            workspace.primary_respiration_root_layer_indices[traversed_layer] = try roots.layerIndex(plant, domain, traversed_layer);
                            workspace.primary_length_m_by_layer[traversed_layer] = roots.axis_primary_length_m[try roots.layerAxisIndex(plant, domain, traversed_layer, axis)];
                        }
                        try ecosys.plant_root_metabolism.allocatePrimaryRootRespiration(
                            roots,
                            workspace.primary_respiration_root_layer_indices[0..traversed_layer_count],
                            workspace.primary_length_m_by_layer[0..traversed_layer_count],
                            workspace.primary_respiration_allocation_fractions[0..traversed_layer_count],
                            roots.axis_depth_m[try roots.axisIndex(plant, domain, axis)],
                            context.plant_water_workspace.*.?.seeding_depth_m[plant],
                            layer > roots.planting_layer_by_plant[plant],
                            respiration,
                        );
                    }
                    if (domain == 0 and layer > roots.planting_layer_by_plant[plant]) {
                        for (0..active_axis_count) |axis| {
                            if (!workspace.primary_active[axis]) continue;
                            const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
                            const axis_index = try roots.axisIndex(plant, domain, axis);
                            if (roots.axis_depth_m[axis_index] >= layer_top_m or
                                roots.axis_secondary_carbon_g[axis_layer] > context.config.absolute_tolerance) continue;
                            const withdrawal_fraction = std.math.clamp(
                                workspace.primary_sink_fractions[axis] + workspace.secondary_sink_fractions[axis],
                                0,
                                1,
                            );
                            try ecosys.plant_root_disturbance.withdrawRootAxisLayer(
                                roots,
                                plant,
                                domain,
                                axis,
                                layer,
                                layer - 1,
                                withdrawal_fraction,
                            );
                        }
                    }
                }
                const shoot_traits = context.shoot_growth_plant_parameters[plant];
                if (shoot_traits.nitrogen_fixation_type >= 1 and shoot_traits.nitrogen_fixation_type <= 3) {
                    const root = try roots.layerIndex(plant, 0, layer);
                    var host_structural_carbon_g_c: f64 = 0;
                    for (0..roots.active_root_axis_count[plant]) |axis| {
                        const axis_layer = try roots.layerAxisIndex(plant, 0, layer, axis);
                        host_structural_carbon_g_c += roots.axis_primary_carbon_g[axis_layer] + roots.axis_secondary_carbon_g[axis_layer];
                    }
                    const host_presence_threshold_g_c =
                        context.runscript.plant_pool_parameters.plant_root_structural_presence_g_per_plant *
                        canopy.plant_population_count[plant];
                    if (host_structural_carbon_g_c > host_presence_threshold_g_c) {
                        const hydrogen_mol_per_m3 = context.soil_chemistry.aqueous[soil].hydrogen;
                        const soil_ph = if (hydrogen_mol_per_m3 > 0) -@log10(hydrogen_mol_per_m3 / 1.0e3) else 7;
                        const environment = try ecosys.plant_root_metabolism.rootEnvironmentResponses(
                            parameters,
                            context.grid.soil_temperature_k[soil],
                            canopy.plant_thermal_adaptation_offset_c[plant],
                            soil_ph,
                            roots.total_water_potential_mpa[root],
                            roots.turgor_water_potential_mpa[root],
                            context.runscript.phenology_parameters.minimum_turgor_potential_mpa,
                            0,
                            traits.primary_root_radius_m,
                            traits.root_profile_type == 0,
                        );
                        const fine_root_litter = try context.plant_litter_partition.*.?.get(plant, .fine_root);
                        const physiological_maturity_reached = branch_development.stage_day[main_branch * 10 + 9] != 0;
                        const symbiosis = try ecosys.plant_root_symbiotic_fixation.calculate(.{
                            .fixation_type = shoot_traits.nitrogen_fixation_type,
                            .first_subhour = true,
                            .restoring_checkpoint = context.restoring_checkpoint,
                            .structural = .{
                                .carbon_g_c = roots.symbiont_structural_carbon_g_c[root],
                                .nitrogen_g_n = roots.symbiont_structural_nitrogen_g_n[root],
                                .phosphorus_g_p = roots.symbiont_structural_phosphorus_g_p[root],
                            },
                            .mobile = .{
                                .carbon_g_c = roots.symbiont_mobile_carbon_g_c[root],
                                .nitrogen_g_n = roots.symbiont_mobile_nitrogen_g_n[root],
                                .phosphorus_g_p = roots.symbiont_mobile_phosphorus_g_p[root],
                            },
                            .host_mobile = .{
                                .carbon_g_c = roots.mobile_carbon_g[root],
                                .nitrogen_g_n = roots.mobile_nitrogen_g[root],
                                .phosphorus_g_p = roots.mobile_phosphorus_g[root],
                            },
                            .host_structural_carbon_g_c = host_structural_carbon_g_c,
                            .host_presence_threshold_g_c = host_presence_threshold_g_c,
                            .cell_area_m2 = context.canopy_cell_area_m2[cell],
                            .temperature_response = try ecosys.plant_root_nutrient_uptake.rootGrowthTemperatureResponse(context.grid.soil_temperature_k[soil], canopy.plant_thermal_adaptation_offset_c[plant], context.runscript.canopy_stress_parameters.growth_temperature),
                            .growth_water_response = environment.growth_water,
                            .maintenance_temperature_response = environment.maintenance_temperature * environment.acidity,
                            .maintenance_water_response = environment.maintenance_water,
                            .oxygen_constraint_fraction = roots.oxygen_process_constraint_fraction[root],
                            .host_exchange_enabled = roots.mobile_carbon_g[root] >
                                context.runscript.plant_pool_parameters.branch_structural_presence_g_per_plant *
                                    canopy.plant_population_count[plant] and
                                (shoot_traits.growth_habit != 0 or !physiological_maturity_reached),
                            .timestep_h = 1,
                        }, context.runscript.symbiotic_fixation_parameters, shoot_traits.symbiont_nitrogen_to_carbon_g_n_per_g_c, shoot_traits.symbiont_phosphorus_to_carbon_g_p_per_g_c, shoot_traits.symbiont_growth_yield_g_c_per_g_c, fine_root_litter);
                        var symbiotic_litterfall: ecosys.plant_root_litterfall.LayerInput = .{};
                        try symbiotic_litterfall.add(symbiosis.litterfall);
                        try ecosys.plant_root_litterfall.validatePublication(context.soil_organic, soil, symbiotic_litterfall);
                        var plant_litterfall: ecosys.plant_root_litterfall.LayerInput = .{ .litter = context.root_litter_products_by_plant[plant] };
                        try plant_litterfall.add(symbiotic_litterfall.litter);
                        inline for (.{
                            roots.actual_respiration_g_c_per_h[root] + symbiosis.respiration_actual_g_c,
                            roots.respiration_unlimited_by_oxygen_g_c_per_h[root] + symbiosis.respiration_oxygen_unlimited_g_c,
                            roots.fixation_uptake_g_n_per_h[plant] + symbiosis.fixed_nitrogen_g_n,
                            roots.fixation_uptake_g_n_per_h_by_layer[root] + symbiosis.fixed_nitrogen_g_n,
                        }) |value| if (!std.math.isFinite(value) or value < 0) return error.NonFiniteRootSymbioticPublication;
                        roots.symbiont_structural_carbon_g_c[root] = symbiosis.structural.carbon_g_c;
                        roots.symbiont_structural_nitrogen_g_n[root] = symbiosis.structural.nitrogen_g_n;
                        roots.symbiont_structural_phosphorus_g_p[root] = symbiosis.structural.phosphorus_g_p;
                        roots.symbiont_mobile_carbon_g_c[root] = symbiosis.mobile.carbon_g_c;
                        roots.symbiont_mobile_nitrogen_g_n[root] = symbiosis.mobile.nitrogen_g_n;
                        roots.symbiont_mobile_phosphorus_g_p[root] = symbiosis.mobile.phosphorus_g_p;
                        roots.mobile_carbon_g[root] = symbiosis.host_mobile.carbon_g_c;
                        roots.mobile_nitrogen_g[root] = symbiosis.host_mobile.nitrogen_g_n;
                        roots.mobile_phosphorus_g[root] = symbiosis.host_mobile.phosphorus_g_p;
                        roots.symbiotic_respiration_actual_g_c_per_h[root] += symbiosis.respiration_actual_g_c;
                        roots.symbiotic_respiration_oxygen_unlimited_g_c_per_h[root] += symbiosis.respiration_oxygen_unlimited_g_c;
                        roots.actual_respiration_g_c_per_h[root] += symbiosis.respiration_actual_g_c;
                        roots.respiration_unlimited_by_oxygen_g_c_per_h[root] += symbiosis.respiration_oxygen_unlimited_g_c;
                        roots.fixation_uptake_g_n_per_h[plant] += symbiosis.fixed_nitrogen_g_n;
                        roots.fixation_uptake_g_n_per_h_by_layer[root] += symbiosis.fixed_nitrogen_g_n;
                        ecosys.plant_root_litterfall.publishValidated(context.soil_organic, soil, symbiotic_litterfall);
                        context.root_litter_carbon_ledger.addValidated(plant, 1, layer, symbiotic_litterfall.litter);
                        context.root_litter_products_by_plant[plant] = plant_litterfall.litter;
                    }
                }
                layer_top_m = layer_bottom_m;
            }
            if (traits.growth_habit != 0) {
                const domain_count = context.root_metabolism_plant_parameters[plant].biologicalDomainCount();
                const active_layer_count = context.grid.active_soil_layer_count[cell];
                const presence_threshold_g_c =
                    context.runscript.plant_pool_parameters.plant_root_structural_presence_g_per_plant *
                    context.plant_water_workspace.*.?.plant_population_count[plant];
                var plant_total_root_carbon_g_c: f64 = 0;
                for (0..domain_count) |domain| for (0..active_layer_count) |layer| {
                    for (0..roots.active_root_axis_count[plant]) |axis| {
                        const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
                        plant_total_root_carbon_g_c += roots.axis_primary_carbon_g[axis_layer] +
                            roots.axis_secondary_carbon_g[axis_layer];
                    }
                };

                // First pass proves that the complete source-ordered sequence
                // is admissible. The second pass commits the same deterministic
                // sequence, so a late failure cannot leave partial transfers.
                var validated_storage_g_c = canopy.plant_seed_storage_carbon_g[plant];
                for (0..domain_count) |domain| for (0..active_layer_count) |layer| {
                    var layer_active_root_carbon_g_c: f64 = 0;
                    for (0..roots.active_root_axis_count[plant]) |axis| {
                        const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
                        layer_active_root_carbon_g_c += roots.axis_primary_carbon_g[axis_layer] +
                            roots.axis_secondary_carbon_g[axis_layer];
                    }
                    const root = try roots.layerIndex(plant, domain, layer);
                    const result = try ecosys.plant_storage_remobilization.replenishDepletedSeasonalStorage(.{
                        .growth_habit = .perennial,
                        .layer_is_rooted = true,
                        .layer_active_root_carbon_g_c = layer_active_root_carbon_g_c,
                        .plant_total_root_carbon_g_c = plant_total_root_carbon_g_c,
                        .layer_mobile_carbon_g_c = roots.mobile_carbon_g[root],
                        .seasonal_storage_carbon_g_c = validated_storage_g_c,
                        .storage_deficit_threshold_g_c_per_g_root_c = context.runscript.storage_remobilization_parameters.depleted_storage_threshold_g_c_per_g_root_c,
                        .exchange_fraction_per_h = parameters.storage_exchange_fraction_per_h,
                        .biological_timestep_h = 1,
                        .presence_threshold_g_c = presence_threshold_g_c,
                    });
                    validated_storage_g_c = result.next_seasonal_storage_carbon_g_c;
                };
                for (0..domain_count) |domain| for (0..active_layer_count) |layer| {
                    var layer_active_root_carbon_g_c: f64 = 0;
                    for (0..roots.active_root_axis_count[plant]) |axis| {
                        const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
                        layer_active_root_carbon_g_c += roots.axis_primary_carbon_g[axis_layer] +
                            roots.axis_secondary_carbon_g[axis_layer];
                    }
                    const root = try roots.layerIndex(plant, domain, layer);
                    const result = try ecosys.plant_storage_remobilization.replenishDepletedSeasonalStorage(.{
                        .growth_habit = .perennial,
                        .layer_is_rooted = true,
                        .layer_active_root_carbon_g_c = layer_active_root_carbon_g_c,
                        .plant_total_root_carbon_g_c = plant_total_root_carbon_g_c,
                        .layer_mobile_carbon_g_c = roots.mobile_carbon_g[root],
                        .seasonal_storage_carbon_g_c = canopy.plant_seed_storage_carbon_g[plant],
                        .storage_deficit_threshold_g_c_per_g_root_c = context.runscript.storage_remobilization_parameters.depleted_storage_threshold_g_c_per_g_root_c,
                        .exchange_fraction_per_h = parameters.storage_exchange_fraction_per_h,
                        .biological_timestep_h = 1,
                        .presence_threshold_g_c = presence_threshold_g_c,
                    });
                    roots.mobile_carbon_g[root] = result.next_layer_mobile_carbon_g_c;
                    canopy.plant_seed_storage_carbon_g[plant] = result.next_seasonal_storage_carbon_g_c;
                };
            }
        }
    }
}

fn applyRootNutrientUptake(context: anytype) !void {
    if (context.plant_roots.* == null or context.plant_water_workspace.* == null or context.plant_root_nutrient_workspace.* == null or context.plant_root_exudation_workspace.* == null) return;
    if (context.runscript.dynamic_plant_salts and context.plant_root_salt_workspace.* == null) return error.MissingRootSaltWorkspace;
    const roots = &context.plant_roots.*.?;
    const water = &context.plant_water_workspace.*.?;
    const grid_workspace = &context.plant_root_nutrient_workspace.*.?;
    const fractions = ecosys.solute_charge_classification.ZoneFractions{
        .ammonium_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
        .ammonium_band = context.runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
        .nitrate_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
        .nitrate_band = context.runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
        .phosphate_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
        .phosphate_band = context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
    };
    const zone_fraction_by_pool = [ecosys.plant_root_nutrient_uptake.nutrient_pool_count]f64{
        fractions.ammonium_non_band,
        fractions.ammonium_band,
        fractions.nitrate_non_band,
        fractions.nitrate_band,
        fractions.phosphate_non_band,
        fractions.phosphate_band,
        fractions.phosphate_non_band,
        fractions.phosphate_band,
    };
    try context.plant_available_nutrients.refreshMineralPools(context.soil_chemistry, context.grid.matrix_liquid_water_m3, fractions);
    for (0..context.grid.cell_count) |cell| {
        var workspace = &grid_workspace.per_cell[cell];
        const salt_workspace = if (context.runscript.dynamic_plant_salts) &context.plant_root_salt_workspace.*.?.per_cell[cell] else null;
        var exudation_workspace = &context.plant_root_exudation_workspace.*.?.per_cell[cell];
        const plant_first = cell * context.config.plant_populations;
        for (0..context.grid.active_soil_layer_count[cell]) |layer| {
            const soil = try context.grid.layerIndex(cell, layer);
            const soil_pools = try context.plant_available_nutrients.mineralPools(soil);
            var total_previous = [_]f64{0} ** ecosys.plant_root_nutrient_uptake.nutrient_pool_count;
            for (0..context.config.plant_populations) |species| {
                const plant = plant_first + species;
                for (0..context.root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| {
                    const root = try roots.layerIndex(plant, domain, layer);
                    for (0..ecosys.plant_root_nutrient_uptake.nutrient_pool_count) |pool_index| {
                        const pool: ecosys.plant_root_nutrient_uptake.NutrientPool = @enumFromInt(pool_index);
                        total_previous[pool_index] += previousRootNutrientDemand(roots, pool, root);
                    }
                }
            }
            total_previous[@intFromEnum(ecosys.plant_root_nutrient_uptake.NutrientPool.ammonium_nonband)] += context.soil_reactive_nitrogen.previous_total_non_band_ammonium_demand_g_n[soil];
            total_previous[@intFromEnum(ecosys.plant_root_nutrient_uptake.NutrientPool.ammonium_band)] += context.soil_reactive_nitrogen.previous_total_band_ammonium_demand_g_n[soil];
            total_previous[@intFromEnum(ecosys.plant_root_nutrient_uptake.NutrientPool.nitrate_nonband)] += context.soil_reactive_nitrogen.previous_total_non_band_nitrate_demand_g_n[soil];
            total_previous[@intFromEnum(ecosys.plant_root_nutrient_uptake.NutrientPool.nitrate_band)] += context.soil_reactive_nitrogen.previous_total_band_nitrate_demand_g_n[soil];
            total_previous[@intFromEnum(ecosys.plant_root_nutrient_uptake.NutrientPool.phosphate_h2_nonband)] += context.soil_microbial_phosphorus.previous_total_non_band_h2po4_demand_g_p[soil];
            total_previous[@intFromEnum(ecosys.plant_root_nutrient_uptake.NutrientPool.phosphate_h2_band)] += context.soil_microbial_phosphorus.previous_total_band_h2po4_demand_g_p[soil];
            total_previous[@intFromEnum(ecosys.plant_root_nutrient_uptake.NutrientPool.phosphate_h_nonband)] += context.soil_microbial_phosphorus.previous_total_non_band_hpo4_demand_g_p[soil];
            total_previous[@intFromEnum(ecosys.plant_root_nutrient_uptake.NutrientPool.phosphate_h_band)] += context.soil_microbial_phosphorus.previous_total_band_hpo4_demand_g_p[soil];

            var competitor_count: usize = 0;
            var exudation_competitor_count: usize = 0;
            for (0..context.config.plant_populations) |species| {
                const plant = plant_first + species;
                if (!context.plant_phenology.*.?.active[plant]) continue;
                for (0..context.root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| {
                    const root = try roots.layerIndex(plant, domain, layer);
                    if (roots.root_surface_area_m2_per_plant[root] <= context.config.absolute_tolerance) continue;
                    if (roots.aqueous_volume_m3[root] > context.config.absolute_tolerance) {
                        exudation_workspace.competitors[exudation_competitor_count] = .{ .plant = plant, .domain = domain, .layer = layer };
                        exudation_competitor_count += 1;
                    }
                    const activity = try ecosys.plant_root_nutrient_uptake.activityFractions(
                        roots.protein_carbon_g[root],
                        roots.total_carbon_g[root],
                        roots.maximum_protein_carbon_g_per_g_c[root],
                        roots.mobile_carbon_g[root],
                        roots.respiration_unlimited_by_carbon_g_c_per_h[root],
                        roots.mobile_carbon_concentration_g_per_g[root],
                        roots.mobile_nitrogen_concentration_g_per_g[root],
                        roots.mobile_phosphorus_concentration_g_per_g[root],
                        context.root_nutrient_feedback_enabled[plant],
                        context.runscript.root_nutrient_parameters,
                    );
                    if (activity.protein <= 0 or activity.carbon <= 0) continue;
                    var competition: [ecosys.plant_root_nutrient_uptake.nutrient_pool_count]f64 = undefined;
                    for (0..ecosys.plant_root_nutrient_uptake.nutrient_pool_count) |pool_index| {
                        const pool: ecosys.plant_root_nutrient_uptake.NutrientPool = @enumFromInt(pool_index);
                        const biome = water.root_biome_fraction[root];
                        competition[pool_index] = if (total_previous[pool_index] > context.config.absolute_tolerance)
                            @max(context.runscript.root_nutrient_parameters.minimum_population_uptake_fraction_multiplier * biome, previousRootNutrientDemand(roots, pool, root) / total_previous[pool_index])
                        else
                            biome;
                    }
                    const temperature_response = try ecosys.plant_root_nutrient_uptake.rootGrowthTemperatureResponse(context.grid.soil_temperature_k[soil], context.detailed_canopy.*.?.plant_thermal_adaptation_offset_c[plant], context.runscript.canopy_stress_parameters.growth_temperature);
                    const water_filled_fraction = if (context.soil_solver_properties.matrix_bulk_volume_m3[soil] > 0)
                        std.math.clamp(context.grid.matrix_liquid_water_m3[soil] / context.soil_solver_properties.matrix_bulk_volume_m3[soil], 0, 1)
                    else
                        0;
                    const diffusivity = [3]f64{
                        try context.runscript.root_nutrient_parameters.diffusivityM2PerH(0, context.grid.soil_temperature_k[soil]),
                        try context.runscript.root_nutrient_parameters.diffusivityM2PerH(1, context.grid.soil_temperature_k[soil]),
                        try context.runscript.root_nutrient_parameters.diffusivityM2PerH(2, context.grid.soil_temperature_k[soil]),
                    };
                    const inputs = try workspace.inputs(competitor_count);
                    try ecosys.plant_root_nutrient_uptake.buildLayerInputs(.{
                        .traits_by_element = context.root_nutrient_traits[plant],
                        .soil_pool_g_element = soil_pools,
                        .total_soil_water_volume_m3 = context.grid.matrix_liquid_water_m3[soil],
                        .zone_fraction_by_pool = zone_fraction_by_pool,
                        .aqueous_diffusivity_m2_per_h_by_element = diffusivity,
                        .liquid_tortuosity = context.runscript.root_nutrient_parameters.liquid_tortuosity_coefficient * water_filled_fraction * water_filled_fraction,
                        .timestep_h = 1,
                        .soil_path_length_m = water.soil_path_length_m[root],
                        .root_cylinder_radius_m = water.root_cylinder_radius_m[root],
                        .root_surface_area_per_radius_m = water.root_surface_area_per_radius_m[root],
                        .root_surface_area_m2_per_plant = roots.root_surface_area_m2_per_plant[root],
                        .root_activity_fraction = std.math.clamp(activity.protein * temperature_response, 0, 1),
                        .nutrient_activity_fraction_by_element = .{ activity.nitrogen, activity.nitrogen, activity.phosphorus },
                        .oxygen_limitation_fraction = std.math.clamp(roots.oxygen_process_constraint_fraction[root], 0, 1),
                        .root_water_uptake_m3_per_plant_step = @max(0, -roots.water_uptake_m3_per_h[root]) / water.plant_population_count[plant],
                        .plant_population_count = water.plant_population_count[plant],
                        .population_competition_fraction_by_pool = competition,
                        .carbon_uptake_limitation_fraction = activity.carbon,
                    }, inputs);
                    workspace.competitors[competitor_count] = .{ .plant = plant, .domain = domain, .layer = layer, .input_by_pool = inputs };
                    competitor_count += 1;
                }
            }
            if (exudation_competitor_count > 0) {
                var substrate_fractions: [ecosys.plant_root_exudation.substrate_count]f64 = undefined;
                try ecosys.soil_heterotrophic_respiration_step.substrateComplexFractions(
                    context.soil_organic,
                    soil,
                    &substrate_fractions,
                    context.config.absolute_tolerance,
                );
                try exudation_workspace.stageLayer(
                    roots,
                    context.soil_organic,
                    soil,
                    context.soil_nitrogen_flux_workspace.layer_biologically_active_water_m3[soil],
                    &substrate_fractions,
                    context.runscript.root_exudation_parameters,
                    context.config.absolute_tolerance,
                    exudation_competitor_count,
                );
            }
            var salt_competitor_count: usize = 0;
            if (exudation_competitor_count > 0 and salt_workspace != null and context.grid.matrix_liquid_water_m3[soil] > context.config.absolute_tolerance) {
                const salt_water_filled_fraction = if (context.soil_solver_properties.matrix_bulk_volume_m3[soil] > 0)
                    std.math.clamp(context.grid.matrix_liquid_water_m3[soil] / context.soil_solver_properties.matrix_bulk_volume_m3[soil], 0, 1)
                else
                    0;
                var soil_salt_content_mol = [ecosys.plant_root_salt_exchange.species_count]f64{
                    context.soil_chemistry.aqueous[soil].aluminum * context.grid.matrix_liquid_water_m3[soil],
                    context.soil_chemistry.aqueous[soil].iron * context.grid.matrix_liquid_water_m3[soil],
                    context.soil_chemistry.aqueous[soil].calcium * context.grid.matrix_liquid_water_m3[soil],
                    context.soil_chemistry.aqueous[soil].magnesium * context.grid.matrix_liquid_water_m3[soil],
                    context.soil_chemistry.aqueous[soil].sodium * context.grid.matrix_liquid_water_m3[soil],
                    context.soil_chemistry.aqueous[soil].potassium * context.grid.matrix_liquid_water_m3[soil],
                    context.soil_chemistry.aqueous[soil].sulfate * context.grid.matrix_liquid_water_m3[soil],
                    context.soil_chemistry.aqueous[soil].chloride * context.grid.matrix_liquid_water_m3[soil],
                };
                for (0..exudation_competitor_count) |competitor_index| {
                    const competitor = exudation_workspace.competitors[competitor_index];
                    const root = try roots.layerIndex(competitor.plant, competitor.domain, competitor.layer);
                    if (roots.aqueous_volume_m3[root] <= context.config.absolute_tolerance) continue;
                    const radial_log = @log((water.soil_path_length_m[root] + water.root_cylinder_radius_m[root]) / water.root_cylinder_radius_m[root]);
                    if (!std.math.isFinite(radial_log) or radial_log <= 0) return error.InvalidRootSaltDiffusionGeometry;
                    salt_workspace.?.competitors[salt_competitor_count] = .{
                        .plant = competitor.plant,
                        .domain = competitor.domain,
                        .layer = competitor.layer,
                        .soil_inventory_fraction = water.root_biome_fraction[root],
                        .root_water_volume_m3 = roots.aqueous_volume_m3[root],
                        .water_advection_m3_per_step = @max(0, -roots.water_uptake_m3_per_h[root]) / water.plant_population_count[competitor.plant],
                        .diffusive_geometry_m = salt_water_filled_fraction * water.root_surface_area_per_radius_m[root] / radial_log,
                        .plant_population_count = water.plant_population_count[competitor.plant],
                    };
                    salt_competitor_count += 1;
                }
                if (salt_competitor_count > 0) _ = try salt_workspace.?.advance(
                    roots,
                    &soil_salt_content_mol,
                    context.grid.matrix_liquid_water_m3[soil],
                    context.grid.soil_temperature_k[soil],
                    context.runscript.root_salt_parameters,
                    salt_competitor_count,
                    .{
                        .absolute_tolerance_mol = context.config.absolute_tolerance,
                        .relative_tolerance = context.config.relative_tolerance,
                        .picard_relaxation = context.config.picard_relaxation,
                        .max_iterations = context.iteration_limits.water_heat_solute_max_iterations,
                    },
                );
                const inverse_water = 1 / context.grid.matrix_liquid_water_m3[soil];
                context.soil_chemistry.aqueous[soil].aluminum = soil_salt_content_mol[0] * inverse_water;
                context.soil_chemistry.aqueous[soil].iron = soil_salt_content_mol[1] * inverse_water;
                context.soil_chemistry.aqueous[soil].calcium = soil_salt_content_mol[2] * inverse_water;
                context.soil_chemistry.aqueous[soil].magnesium = soil_salt_content_mol[3] * inverse_water;
                context.soil_chemistry.aqueous[soil].sodium = soil_salt_content_mol[4] * inverse_water;
                context.soil_chemistry.aqueous[soil].potassium = soil_salt_content_mol[5] * inverse_water;
                context.soil_chemistry.aqueous[soil].sulfate = soil_salt_content_mol[6] * inverse_water;
                context.soil_chemistry.aqueous[soil].chloride = soil_salt_content_mol[7] * inverse_water;
            }
            if (competitor_count > 0) try workspace.advanceAssimilating(
                roots,
                soil_pools,
                competitor_count,
                context.runscript.root_metabolism_parameters.nutrient_uptake_respiration_g_c_per_g_element,
            );
            if (exudation_competitor_count > 0) {
                for (0..ecosys.plant_root_exudation.substrate_count) |substrate| {
                    var total: ecosys.plant_root_exudation.Result = .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
                    for (0..exudation_competitor_count) |competitor_index| {
                        const result = exudation_workspace.staged_results[competitor_index * ecosys.plant_root_exudation.substrate_count + substrate];
                        total.carbon_g_c += result.carbon_g_c;
                        total.nitrogen_g_n += result.nitrogen_g_n;
                        total.phosphorus_g_p += result.phosphorus_g_p;
                    }
                    const available_index = try context.plant_available_nutrients.organicIndex(soil, substrate);
                    inline for (.{
                        context.plant_available_nutrients.organic_carbon_change_g_c_per_h[available_index] - total.carbon_g_c,
                        context.plant_available_nutrients.organic_nitrogen_change_g_n_per_h[available_index] - total.nitrogen_g_n,
                        context.plant_available_nutrients.organic_phosphorus_change_g_p_per_h[available_index] - total.phosphorus_g_p,
                    }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootExudationMirrorCommit;
                }
                try exudation_workspace.commitLayer(
                    roots,
                    context.soil_organic,
                    soil,
                    exudation_competitor_count,
                );
                for (0..ecosys.plant_root_exudation.substrate_count) |substrate| {
                    var total: ecosys.plant_root_exudation.Result = .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
                    for (0..exudation_competitor_count) |competitor_index| {
                        const result = exudation_workspace.staged_results[competitor_index * ecosys.plant_root_exudation.substrate_count + substrate];
                        total.carbon_g_c += result.carbon_g_c;
                        total.nitrogen_g_n += result.nitrogen_g_n;
                        total.phosphorus_g_p += result.phosphorus_g_p;
                    }
                    const available_index = try context.plant_available_nutrients.organicIndex(soil, substrate);
                    const authoritative = context.soil_organic.dissolved[soil * ecosys.plant_root_exudation.substrate_count + substrate];
                    context.plant_available_nutrients.organic_carbon_g_c[available_index] = authoritative.carbon_g_c;
                    context.plant_available_nutrients.organic_nitrogen_g_n[available_index] = authoritative.nitrogen_g_n;
                    context.plant_available_nutrients.organic_phosphorus_g_p[available_index] = authoritative.phosphorus_g_p;
                    context.plant_available_nutrients.organic_carbon_change_g_c_per_h[available_index] -= total.carbon_g_c;
                    context.plant_available_nutrients.organic_nitrogen_change_g_n_per_h[available_index] -= total.nitrogen_g_n;
                    context.plant_available_nutrients.organic_phosphorus_change_g_p_per_h[available_index] -= total.phosphorus_g_p;
                }
            }
            if ((competitor_count > 0 or salt_competitor_count > 0) and context.runscript.dynamic_plant_salts) {
                var hydrogen_charge_mol: f64 = 0;
                for (0..competitor_count) |competitor_index| {
                    const result_base = competitor_index * ecosys.plant_root_nutrient_uptake.nutrient_pool_count;
                    hydrogen_charge_mol += try ecosys.plant_root_ion_balance.hydrogenChargeMol(
                        workspace.staged_results[result_base..][0..ecosys.plant_root_nutrient_uptake.nutrient_pool_count],
                        context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                        context.runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                        .{},
                    );
                }
                const zero_nutrient_results = [_]ecosys.plant_root_nutrient_uptake.Result{.{
                    .demand_g_element = 0,
                    .uptake_g_element = 0,
                    .oxygen_unlimited_uptake_g_element = 0,
                    .carbon_unlimited_uptake_g_element = 0,
                    .available_g_element = 0,
                }} ** ecosys.plant_root_nutrient_uptake.nutrient_pool_count;
                for (0..salt_competitor_count) |competitor_index| {
                    const competitor = salt_workspace.?.competitors[competitor_index];
                    var salts: ecosys.plant_root_ion_balance.SaltUptakeMol = .{};
                    inline for (@typeInfo(ecosys.plant_root_system.SaltSpecies).@"enum".fields) |field| {
                        const species: ecosys.plant_root_system.SaltSpecies = @enumFromInt(field.value);
                        @field(salts, field.name) = roots.salt_uptake_mol_per_h[try roots.saltIndex(competitor.plant, competitor.domain, competitor.layer, species)];
                    }
                    hydrogen_charge_mol += try ecosys.plant_root_ion_balance.hydrogenChargeMol(
                        &zero_nutrient_results,
                        context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                        context.runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                        salts,
                    );
                }
                _ = try ecosys.plant_root_ion_balance.commitHydrogenCharge(
                    &context.soil_chemistry.aqueous[soil].hydrogen,
                    context.grid.matrix_liquid_water_m3[soil],
                    hydrogen_charge_mol,
                    1,
                );
            }
        }
    }
    try context.plant_available_nutrients.publishMineralPools(context.soil_chemistry, context.grid.matrix_liquid_water_m3, fractions);
}

fn ensureCanopyGrowthNodeTopology(canopy: *ecosys.canopy_photosynthesis.State, growth_stages: *const ecosys.plant_growth_stages.State, samples_per_node: usize) !void {
    if (samples_per_node == 0 or growth_stages.branches.len != canopy.branch_node_offsets.len - 1) return error.CanopyGrowthTopologyDimensionMismatch;
    for (growth_stages.branches, 0..) |stage, branch| {
        const desired_node_count = try std.math.add(usize, stage.newest_growing_leaf_ordinal, 1);
        while (true) {
            const nodes = try canopy.nodeRange(branch);
            if (nodes.end - nodes.first >= desired_node_count) break;
            _ = try canopy.appendNode(branch, samples_per_node);
        }
    }
}

fn convergeSurfaceLitterChemistry(context: anytype) !void {
    if (context.chemistry_reaction_parameters.*) |reaction_parameters| if (context.organic_parameters.*) |*organic_runtime_parameters| {
        context.surface_litter_chemistry_diagnostics.reset();
        var litter_chemistry_context: ecosys.surface_litter_chemistry_step.ApplyContext = .{
            .state = context.surface_litter_chemistry,
            .surface_organic = context.surface_organic,
            .litter_water_m3 = context.surface_precipitation.litter_water_m3,
            .chemistry_parameters = reaction_parameters,
            .cation_selectivity_by_cell = context.surface_litter_cation_selectivity,
            .litter_dry_mass_Mg_per_g_c = organic_runtime_parameters.surface_litter_dry_mass_Mg_per_g_c,
            .dynamic_salts = context.runscript.dynamic_plant_salts,
            .solver_options = .{
                .absolute_tolerance = context.config.absolute_tolerance,
                .relative_tolerance = context.config.relative_tolerance,
                .picard_relaxation = context.config.picard_relaxation,
                // Surface litter is seeded from topsoil and must close the
                // Dynamic salts use STARTE's MRXN=1000 equilibrium ceiling.
                // The fixed-pH ISALTG=0 source branch is outside that loop;
                // its already-hourly kinetic coefficients are applied once.
                .max_iterations = context.iteration_limits.initial_solute_reaction_max_iterations,
            },
            .diagnostics = context.surface_litter_chemistry_diagnostics,
        };
        try runKernelAcrossSerialTiles(context, &litter_chemistry_context, ecosys.surface_litter_chemistry_step.applyTile);
    };
}

noinline fn runSoilBiogeochemistryBatch(context: anytype) !void {
    if (context.soil_nitrogen_parameters.*) |nitrogen_parameters| {
        var reset_nitrogen_fluxes: ecosys.soil_nitrogen_flux_workspace.ResetContext = .{ .state = context.soil_nitrogen_flux_workspace };
        try runScienceCellLayers(context, &reset_nitrogen_fluxes, ecosys.soil_nitrogen_flux_workspace.resetTile);
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
        try runScienceCellLayers(context, &nitrifier_environment_context, ecosys.soil_nitrifier_environment_step.applyTile);
        var nitrification_context: ecosys.soil_nitrification_step.ApplyContext = .{ .result = context.soil_nitrogen_flux_workspace, .reactive_nitrogen = context.soil_reactive_nitrogen, .chemistry_state = context.soil_chemistry, .model_grid = context.grid, .zone_fractions = zone_fractions, .roles = context.soil_nitrifier_environment.roles, .temperature_water_activity = context.soil_nitrifier_environment.temperature_water_activity, .nitrogen_phosphorus_activity = context.soil_nitrifier_environment.nitrogen_phosphorus_activity, .aqueous_co2_activity = context.soil_nitrifier_environment.aqueous_co2_activity, .active_biomass_g_c = context.soil_nitrifier_environment.active_biomass_g_c, .microbial_active_fraction = context.soil_nitrifier_environment.microbial_active_fraction, .parameters = nitrogen_parameters.nitrification, .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol, .timestep_h = 1, .negligible_demand_g_n = context.config.absolute_tolerance };
        try runScienceCellLayers(context, &nitrification_context, ecosys.soil_nitrification_step.applyTile);
        var heterotrophic_respiration_context: ecosys.soil_heterotrophic_respiration_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .reactive_nitrogen = context.soil_reactive_nitrogen,
            .organic_state = context.soil_organic,
            .microbial_state = context.soil_microbial,
            .model_grid = context.grid,
            .retention_curve = context.soil_solver_properties.retention_curve,
            .matrix_bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3,
            .matric_plus_osmotic_potential_mpa = context.soil_hourly_workspace.matric_plus_osmotic_potential_mpa,
            .parameters = nitrogen_parameters,
            .timestep_h = 1,
            .negligible_amount = context.config.absolute_tolerance,
        };
        try runScienceCellLayers(context, &heterotrophic_respiration_context, ecosys.soil_heterotrophic_respiration_step.applyTile);
        const surface_parameters = context.surface_gas_parameters.* orelse return error.SoilNitrogenRequiresAtmosphericGasParameters;
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
        try runScienceCellLayers(context, &soil_oxygen_context, ecosys.soil_oxygen_step.applyTile);
        var microbial_mixing_activity_context: ecosys.soil_microbial_layer_mixing.PrepareContext = .{
            .result = context.soil_microbial_layer_mixing,
            .respiration_fluxes = context.soil_nitrogen_flux_workspace,
            .oxygen_allocation = context.soil_microbial_oxygen,
        };
        try runScienceCellLayers(context, &microbial_mixing_activity_context, ecosys.soil_microbial_layer_mixing.prepareActivityTile);
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
        try runScienceCellLayers(context, &heterotrophic_denitrification_context, ecosys.soil_heterotrophic_denitrification_step.applyTile);
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
        try runScienceCellLayers(context, &autotrophic_denitrification_context, ecosys.soil_autotrophic_denitrification_step.applyTile);
        var soil_microbial_nitrogen_exchange_context: ecosys.soil_microbial_nitrogen_exchange_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .reactive_nitrogen = context.soil_reactive_nitrogen,
            .microbial_state = context.soil_microbial,
            .chemistry_state = context.soil_chemistry,
            .model_grid = context.grid,
            .matric_plus_osmotic_potential_mpa = context.soil_hourly_workspace.matric_plus_osmotic_potential_mpa,
            .zone_fractions = zone_fractions,
            .parameters = nitrogen_parameters,
            .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
            .timestep_h = 1,
            .negligible_amount = context.config.absolute_tolerance,
        };
        try runScienceCellLayers(context, &soil_microbial_nitrogen_exchange_context, ecosys.soil_microbial_nitrogen_exchange_step.applyTile);
        var soil_microbial_phosphorus_exchange_context: ecosys.soil_microbial_phosphorus_exchange_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .phosphorus_history = context.soil_microbial_phosphorus,
            .microbial_state = context.soil_microbial,
            .chemistry_state = context.soil_chemistry,
            .model_grid = context.grid,
            .matric_plus_osmotic_potential_mpa = context.soil_hourly_workspace.matric_plus_osmotic_potential_mpa,
            .zone_fractions = zone_fractions,
            .parameters = nitrogen_parameters,
            .timestep_h = 1,
            .negligible_amount = context.config.absolute_tolerance,
        };
        try runScienceCellLayers(context, &soil_microbial_phosphorus_exchange_context, ecosys.soil_microbial_phosphorus_exchange_step.applyTile);
        var soil_microbial_maintenance_context: ecosys.soil_microbial_maintenance_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .microbial_state = context.soil_microbial,
            .organic_state = context.soil_organic,
            .chemistry_state = context.soil_chemistry,
            .model_grid = context.grid,
            .oxygen_state = context.soil_microbial_oxygen,
            .matric_plus_osmotic_potential_mpa = context.soil_hourly_workspace.matric_plus_osmotic_potential_mpa,
            .parameters = nitrogen_parameters,
            .timestep_h = 1,
        };
        try runScienceCellLayers(context, &soil_microbial_maintenance_context, ecosys.soil_microbial_maintenance_step.applyTile);
        var autotrophic_carbon_context: ecosys.soil_autotrophic_carbon_step.ApplyContext = .{ .result = context.soil_autotrophic_carbon, .flux_workspace = context.soil_nitrogen_flux_workspace, .oxygen_state = context.soil_microbial_oxygen, .roles = context.soil_nitrifier_environment.roles, .parameters = nitrogen_parameters };
        try runScienceCellLayers(context, &autotrophic_carbon_context, ecosys.soil_autotrophic_carbon_step.applyTile);
        var soil_nitrogen_fixation_context: ecosys.soil_nonsymbiotic_nitrogen_fixation_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .microbial_state = context.soil_microbial,
            .gas_state = context.gas_transport,
            .water_volume_m3 = context.grid.matrix_liquid_water_m3,
            .parameters = nitrogen_parameters,
            .timestep_h = 1,
        };
        try runScienceCellLayers(context, &soil_nitrogen_fixation_context, ecosys.soil_nonsymbiotic_nitrogen_fixation_step.applyTile);
        var soil_microbial_substrate_uptake_context: ecosys.soil_microbial_substrate_uptake_step.ApplyContext = .{
            .result = context.soil_nitrogen_flux_workspace,
            .organic_state = context.soil_organic,
            .microbial_state = context.soil_microbial,
            .oxygen_state = context.soil_microbial_oxygen,
            .parameters = nitrogen_parameters,
            .negligible_amount = context.config.absolute_tolerance,
        };
        try runScienceCellLayers(context, &soil_microbial_substrate_uptake_context, ecosys.soil_microbial_substrate_uptake_step.applyTile);
        var soil_respiration_products_context: ecosys.soil_respiration_products_step.ApplyContext = .{ .result = context.soil_respiration_products, .microbial_state = context.soil_microbial, .respiration_fluxes = context.soil_nitrogen_flux_workspace };
        try runScienceCellLayers(context, &soil_respiration_products_context, ecosys.soil_respiration_products_step.applyTile);
        if (nitrogen_parameters.methane) |methane_parameters| {
            var prepare_methane_context: ecosys.soil_methane_step.PrepareContext = .{
                .result = context.soil_methane,
                .microbial_state = context.soil_microbial,
                .respiration_products = context.soil_respiration_products,
                .gas_state = context.gas_transport,
                .water_volume_m3 = context.grid.matrix_liquid_water_m3,
                .soil_temperature_k = context.grid.soil_temperature_k,
                .matric_plus_osmotic_potential_mpa = context.soil_hourly_workspace.matric_plus_osmotic_potential_mpa,
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
            try runScienceCellLayers(context, &prepare_methane_context, ecosys.soil_methane_step.prepareTile);
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
            try runScienceCellLayers(context, &methane_context, ecosys.soil_methane_step.applyTile);
        }
        if (context.organic_parameters.*) |*organic_runtime_parameters| {
            var soil_microbial_assimilation_context: ecosys.soil_microbial_assimilation_step.ApplyContext = .{
                .result = context.soil_nitrogen_flux_workspace,
                .microbial_state = context.soil_microbial,
                .model_grid = context.grid,
                .matric_plus_osmotic_potential_mpa = context.soil_hourly_workspace.matric_plus_osmotic_potential_mpa,
                .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.microbial_nitrogen_to_carbon,
                .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.microbial_phosphorus_to_carbon,
                .parameters = nitrogen_parameters,
                .timestep_h = 1,
            };
            try runScienceCellLayers(context, &soil_microbial_assimilation_context, ecosys.soil_microbial_assimilation_step.applyTile);
            var soil_microbial_turnover_context: ecosys.soil_microbial_turnover_step.ApplyContext = .{
                .result = context.soil_microbial_turnover,
                .microbial_state = context.soil_microbial,
                .organic_state = context.soil_organic,
                .maintenance_fluxes = context.soil_nitrogen_flux_workspace,
                .model_grid = context.grid,
                .clay_mass_fraction = context.soil_solver_properties.clay_mass_fraction,
                .matric_plus_osmotic_potential_mpa = context.soil_hourly_workspace.matric_plus_osmotic_potential_mpa,
                .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.microbial_nitrogen_to_carbon,
                .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.microbial_phosphorus_to_carbon,
                .substrate_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.substrate_nitrogen_to_carbon,
                .substrate_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.substrate_phosphorus_to_carbon,
                .parameters = nitrogen_parameters,
                .timestep_h = 1,
                .negligible_carbon_g_c = context.config.absolute_tolerance,
            };
            try runScienceCellLayers(context, &soil_microbial_turnover_context, ecosys.soil_microbial_turnover_step.applyTile);
        }
        if (context.organic_parameters.*) |organic_runtime_parameters| {
            var soil_organic_priming_context: ecosys.soil_organic_priming_step.ApplyContext = .{
                .result = context.soil_organic_priming,
                .organic_state = context.soil_organic,
                .microbial_state = context.soil_microbial,
                .respiration_fluxes = context.soil_nitrogen_flux_workspace,
                .soil_temperature_k = context.grid.soil_temperature_k,
                .matric_plus_osmotic_potential_mpa = context.soil_hourly_workspace.matric_plus_osmotic_potential_mpa,
                .thermal_adaptation_offset_k = nitrogen_parameters.microbial_thermal_adaptation_offset_k,
                .dissolved_priming_rate_per_h = organic_runtime_parameters.soil_dissolved_priming_rate_per_h,
                .microbial_priming_rate_per_h = organic_runtime_parameters.soil_microbial_priming_rate_per_h,
                .timestep_h = 1,
                .negligible_carbon_g_c = context.config.absolute_tolerance,
            };
            try runScienceCellLayers(context, &soil_organic_priming_context, ecosys.soil_organic_priming_step.applyTile);
            var soil_organic_decomposition_context: ecosys.soil_organic_decomposition_step.ApplyContext = .{
                .result = context.soil_organic_decomposition,
                .organic_state = context.soil_organic,
                .microbial_state = context.soil_microbial,
                .respiration_fluxes = context.soil_nitrogen_flux_workspace,
                .priming = context.soil_organic_priming,
                .soil_temperature_k = context.grid.soil_temperature_k,
                .matrix_bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3,
                .bulk_density_Mg_per_m3 = context.soil_solver_properties.bulk_density_megagrams_per_m3,
                .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.microbial_nitrogen_to_carbon,
                .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.microbial_phosphorus_to_carbon,
                .substrate_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.substrate_nitrogen_to_carbon,
                .substrate_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.substrate_phosphorus_to_carbon,
                .thermal_adaptation_offset_k = nitrogen_parameters.microbial_thermal_adaptation_offset_k,
                .parameters = organic_runtime_parameters.soil_organic_decomposition,
                .timestep_h = 1,
                .negligible_carbon_g_c = context.config.absolute_tolerance,
            };
            try runScienceCellLayers(context, &soil_organic_decomposition_context, ecosys.soil_organic_decomposition_step.applyTile);
            var soil_organic_sorption_context: ecosys.soil_organic_sorption_step.ApplyContext = .{
                .result = context.soil_organic_sorption,
                .organic_state = context.soil_organic,
                .microbial_state = context.soil_microbial,
                .respiration_fluxes = context.soil_nitrogen_flux_workspace,
                .water_volume_m3 = context.grid.matrix_liquid_water_m3,
                .matrix_bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3,
                .bulk_density_Mg_per_m3 = context.soil_solver_properties.bulk_density_megagrams_per_m3,
                .anion_exchange_capacity_mol_per_Mg = context.soil_solver_properties.anion_exchange_capacity_mol_per_Mg,
                .sorption_rate_per_h = organic_runtime_parameters.soil_organic_sorption_rate_per_h,
                .adsorption_coefficient = organic_runtime_parameters.soil_organic_adsorption_coefficient,
                .timestep_h = 1,
                .negligible_amount_g = context.config.absolute_tolerance,
            };
            try runScienceCellLayers(context, &soil_organic_sorption_context, ecosys.soil_organic_sorption_step.applyTile);
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
        try runScienceCellLayers(context, &soil_litter_colonization_context, ecosys.soil_litter_colonization_step.applyTile);
        var chemodenitrification_context: ecosys.soil_chemodenitrification_step.ApplyContext = .{ .result = context.soil_nitrogen_flux_workspace, .reactive_nitrogen = context.soil_reactive_nitrogen, .chemistry_state = context.soil_chemistry, .model_grid = context.grid, .zone_fractions = zone_fractions, .parameters = nitrogen_parameters, .timestep_h = 1 };
        try runScienceCellLayers(context, &chemodenitrification_context, ecosys.soil_chemodenitrification_step.applyTile);
        var gas_aggregation_context: ecosys.soil_biogeochemical_gas_aggregation.ProcessContext = .{
            .result = context.soil_biogeochemical_gas_fluxes,
            .autotrophic_carbon = context.soil_autotrophic_carbon,
            .respiration_products = context.soil_respiration_products,
            .methane = if (nitrogen_parameters.methane != null) context.soil_methane else null,
            .oxygen = context.soil_microbial_oxygen,
            .nitrogen_fluxes = context.soil_nitrogen_flux_workspace,
            .redox_satisfaction_fraction = context.soil_redox_satisfaction_fraction,
        };
        try runScienceCellLayers(context, &gas_aggregation_context, ecosys.soil_biogeochemical_gas_aggregation.aggregateProcessTile);
        var nitrogen_commit_context: ecosys.soil_nitrogen_commit.ApplyContext = .{ .reactive_nitrogen = context.soil_reactive_nitrogen, .phosphorus_history = context.soil_microbial_phosphorus, .chemistry_state = context.soil_chemistry, .gas_state = context.gas_transport, .organic_state = context.soil_organic, .microbial_state = context.soil_microbial, .flux_workspace = context.soil_nitrogen_flux_workspace, .microbial_turnover = context.soil_microbial_turnover, .litter_colonization = context.soil_litter_colonization, .organic_sorption = context.soil_organic_sorption, .organic_decomposition = context.soil_organic_decomposition, .organic_priming = context.soil_organic_priming, .respiration_products = context.soil_respiration_products, .methane = if (nitrogen_parameters.methane != null) context.soil_methane else null, .methane_autotrophic_substrate_index = nitrogen_parameters.nitrifier_indices.autotrophic_substrate_index, .methanotroph_population_index = if (nitrogen_parameters.methane) |methane_parameters| methane_parameters.methanotroph_population_index else 0, .humus_partition_by_cell = context.topsoil_humus_partition, .soil_layer_capacity = context.grid.soil_layer_capacity, .water_volume_m3 = context.grid.matrix_liquid_water_m3, .zone_fractions = zone_fractions, .oxygen_satisfaction_fraction = context.soil_microbial_oxygen.demand_satisfaction_fraction, .redox_satisfaction_fraction = context.soil_redox_satisfaction_fraction, .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol, .phosphorus_molar_mass_g_per_mol = nitrogen_parameters.microbial_mineral_exchange.phosphorus_molar_mass_g_per_mol, .tolerance_g_n = context.config.absolute_tolerance, .dynamic_salts = context.runscript.dynamic_plant_salts, .timestep_h = 1, .hourly_signed_heterotrophic_respiration_g_c = context.daily_heterotrophic_respiration.soil_hourly_signed_g_c, .hourly_carbon_dioxide_production_g_c = context.daily_heterotrophic_respiration.soil_hourly_carbon_dioxide_production_g_c };
        try runScienceCellLayers(context, &nitrogen_commit_context, ecosys.soil_nitrogen_commit.applyTile);
        var autotrophic_carbon_commit_context: ecosys.soil_autotrophic_carbon_step.CommitContext = .{ .state = context.soil_autotrophic_carbon, .gas_state = context.gas_transport, .microbial_state = context.soil_microbial, .autotrophic_substrate_index = nitrogen_parameters.nitrifier_indices.autotrophic_substrate_index, .tolerance_g_c = context.config.absolute_tolerance };
        try runScienceCellLayers(context, &autotrophic_carbon_commit_context, ecosys.soil_autotrophic_carbon_step.commitTile);
        var microbial_layer_mixing_context: ecosys.soil_microbial_layer_mixing.ApplyContext = .{
            .microbial_state = context.soil_microbial,
            .active_layer_count = context.grid.active_soil_layer_count,
            .layer_volume_m3 = context.soil_solver_properties.layer_volume_m3,
            .dry_bulk_density_Mg_per_m3 = context.soil_solver_properties.bulk_density_megagrams_per_m3,
            .layer_thickness_m = context.soil_solver_properties.layer_thickness_m,
            .total_organic_carbon_g_per_megagram = context.soil_solver_properties.total_organic_carbon_g_per_megagram,
            .substrate_unlimited_oxygen_limited_activity_g_c = context.soil_microbial_layer_mixing.substrate_unlimited_oxygen_limited_activity_g_c,
            .parameters = .{ .mixing_rate_per_h = nitrogen_parameters.microbial_layer_mixing_rate_per_h, .timestep_h = 1, .minimum_mixing_layer_thickness_m = context.runscript.soil_geometry_parameters.minimum_layer_thickness_m },
        };
        try runScienceCells(context, &microbial_layer_mixing_context, ecosys.soil_microbial_layer_mixing.applyTile);
    }
}

noinline fn runSoilBiogeochemistryBySerialTile(context: anytype) !void {
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
        .field_capacity_potential_mpa = context.surface_field_capacity_potential_mpa,
        .wilting_point_potential_mpa = context.surface_wilting_point_potential_mpa,
        .mean_annual_temperature_c = context.mean_annual_temperature_c_by_cell,
        .parameters = surface_parameters.litter_water_environment,
    };
    try runScienceCells(context, &litter_water_environment_context, ecosys.surface_litter_water_environment.applyTile);
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
    try runScienceCells(context, &microbial_environment_context, ecosys.surface_microbial_environment_step.applyTile);
    var microbial_respiration_context: ecosys.surface_microbial_respiration_step.ApplyContext = .{
        .result = context.surface_microbial_respiration,
        .surface_organic = context.surface_organic,
        .biologically_active_water_m3 = context.surface_microbial_environment.biologically_active_water_m3,
        .growth_temperature_response = context.surface_microbial_environment.growth_temperature_response,
        .matric_plus_osmotic_potential_mpa = context.surface_litter_water_environment.matric_plus_osmotic_water_potential_mpa,
        .timestep_h = 1,
        .parameters = surface_parameters.microbial_respiration,
    };
    try runScienceCells(context, &microbial_respiration_context, ecosys.surface_microbial_respiration_step.applyTile);
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
    try runScienceCells(context, &microbial_oxygen_context, ecosys.surface_microbial_oxygen_driver.applyTile);
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
    try runScienceCells(context, &microbial_maintenance_context, ecosys.surface_microbial_maintenance_step.applyTile);
    if (context.organic_parameters.*) |*organic_runtime_parameters| {
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
        try runScienceCells(context, &nitrogen_fixation_context, ecosys.surface_nonsymbiotic_nitrogen_fixation_step.applyTile);
        var assimilation_context: ecosys.surface_microbial_assimilation_step.ApplyContext = .{
            .result = context.surface_microbial_assimilation,
            .surface_organic = context.surface_organic,
            .growth_temperature_response = context.surface_microbial_environment.growth_temperature_response,
            .matric_plus_osmotic_potential_mpa = context.surface_litter_water_environment.matric_plus_osmotic_water_potential_mpa,
            .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.surfaceMicrobialNitrogenToCarbon(),
            .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.surfaceMicrobialPhosphorusToCarbon(),
            .labile_structural_fraction = surface_parameters.microbial_respiration.labile_biomass_fraction,
            .transfer_rate_per_h = surface_parameters.microbial_respiration.nonstructural_to_structural_rate_per_h,
            .timestep_h = 1,
        };
        try runScienceCells(context, &assimilation_context, ecosys.surface_microbial_assimilation_step.applyTile);
        var mineral_exchange_context: ecosys.surface_microbial_mineral_exchange_step.ApplyContext = .{
            .result = context.surface_microbial_mineral_exchange,
            .surface_organic = context.surface_organic,
            .litter_chemistry = context.surface_litter_chemistry,
            .litter_water_m3 = context.surface_precipitation.litter_water_m3,
            .growth_temperature_response = context.surface_microbial_environment.growth_temperature_response,
            .matric_plus_osmotic_potential_mpa = context.surface_litter_water_environment.matric_plus_osmotic_water_potential_mpa,
            .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.surfaceMicrobialNitrogenToCarbon(),
            .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.surfaceMicrobialPhosphorusToCarbon(),
            .labile_biomass_fraction = surface_parameters.microbial_respiration.labile_biomass_fraction,
            .microbial_surface_area_m2_per_g_c = 4 * std.math.pi * surface_parameters.microbial_radius_m * surface_parameters.microbial_radius_m * surface_parameters.microbial_count_per_g_c,
            .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
            .parameters = surface_parameters.mineral_exchange,
            .timestep_h = 1,
        };
        try runScienceCells(context, &mineral_exchange_context, ecosys.surface_microbial_mineral_exchange_step.applyTile);
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
            .matric_plus_osmotic_potential_mpa = context.surface_litter_water_environment.matric_plus_osmotic_water_potential_mpa,
            .microbial_nitrogen_to_carbon_g_n_per_g_c = organic_runtime_parameters.surfaceMicrobialNitrogenToCarbon(),
            .microbial_phosphorus_to_carbon_g_p_per_g_c = organic_runtime_parameters.surfaceMicrobialPhosphorusToCarbon(),
            .labile_biomass_fraction = surface_parameters.microbial_respiration.labile_biomass_fraction,
            .microbial_surface_area_m2_per_g_c = 4 * std.math.pi * surface_parameters.microbial_radius_m * surface_parameters.microbial_radius_m * surface_parameters.microbial_count_per_g_c,
            .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
            .parameters = surface_parameters.mineral_exchange,
            .timestep_h = 1,
        };
        try runScienceCells(context, &topsoil_exchange_context, ecosys.surface_topsoil_mineral_exchange_step.applyTile);
        var turnover_context: ecosys.surface_microbial_turnover_step.ApplyContext = .{
            .result = context.surface_microbial_turnover,
            .surface_organic = context.surface_organic,
            .maintenance = context.surface_microbial_maintenance,
            .growth_temperature_response = context.surface_microbial_environment.growth_temperature_response,
            .matric_plus_osmotic_potential_mpa = context.surface_litter_water_environment.matric_plus_osmotic_water_potential_mpa,
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
        try runScienceCells(context, &turnover_context, ecosys.surface_microbial_turnover_step.applyTile);
        var priming_context: ecosys.surface_organic_priming_step.ApplyContext = .{
            .result = context.surface_organic_priming,
            .surface_organic = context.surface_organic,
            .respiration = context.surface_microbial_respiration,
            .substrate_uptake = context.surface_microbial_substrate_uptake,
            .environment = context.surface_microbial_environment,
            .matric_plus_osmotic_potential_mpa = context.surface_litter_water_environment.matric_plus_osmotic_water_potential_mpa,
            .dissolved_priming_rate_per_h = surface_parameters.microbial_turnover.dissolved_priming_rate_per_h,
            .microbial_priming_rate_per_h = surface_parameters.microbial_turnover.microbial_priming_rate_per_h,
            .timestep_h = 1,
            .negligible_carbon_g_c = context.config.absolute_tolerance,
        };
        try runScienceCells(context, &priming_context, ecosys.surface_organic_priming_step.applyTile);
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
        try runScienceCells(context, &organic_decomposition_context, ecosys.surface_organic_decomposition_step.applyTile);
        var organic_sorption_context: ecosys.surface_organic_sorption_step.ApplyContext = .{
            .result = context.surface_organic_sorption,
            .surface_organic = context.surface_organic,
            .litter_water_m3 = context.surface_precipitation.litter_water_m3,
            .litter_dry_mass_Mg = context.surface_litter_geometry.dry_mass_Mg,
            .timestep_h = 1,
            .negligible_mass_g = context.config.absolute_tolerance,
            .parameters = surface_parameters.organic_sorption,
        };
        try runScienceCells(context, &organic_sorption_context, ecosys.surface_organic_sorption_step.applyTile);
        var litter_colonization_context: ecosys.surface_litter_colonization_step.ApplyContext = .{
            .result = context.surface_litter_colonization,
            .surface_organic = context.surface_organic,
            .decomposition = context.surface_organic_decomposition,
            .parameters = surface_parameters.litter_colonization,
            .negligible_carbon_g_c = context.config.absolute_tolerance,
        };
        try runScienceCells(context, &litter_colonization_context, ecosys.surface_litter_colonization_step.applyTile);
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
    try runScienceCells(context, &denitrification_context, ecosys.surface_denitrification_step.applyTile);
    var denitrification_map = SurfaceDenitrificationRespirationMap{
        .destination_g_c = context.surface_microbial_substrate_uptake.denitrification_respiration_g_c,
        .source_g_c = context.surface_denitrification.respiration_g_c,
    };
    try runScienceCells(context, &denitrification_map, mapSurfaceDenitrificationRespiration);
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
    try runScienceCells(context, &substrate_uptake_context, ecosys.surface_microbial_substrate_uptake_step.applyTile);
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
    try runScienceCells(context, &metabolism_commit_context, ecosys.surface_metabolism_commit.applyTile);
    // NITRO 4168--4275 evaluates the L=0 -> NU microbial transfer only after
    // both surface and first-soil-layer biological states are accepted.
    if (context.soil_nitrogen_parameters.*) |nitrogen_parameters| {
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
            .surface_dry_mass_Mg = context.surface_litter_geometry.dry_mass_Mg,
            .topsoil_bulk_density_Mg_per_m3 = context.soil_solver_properties.bulk_density_megagrams_per_m3,
            .topsoil_organic_carbon_g_per_megagram = context.soil_solver_properties.total_organic_carbon_g_per_megagram,
            .parameters = .{
                .mixing_rate_per_h = nitrogen_parameters.microbial_layer_mixing_rate_per_h,
                .timestep_h = 1,
                .minimum_layer_thickness_m = context.runscript.soil_geometry_parameters.minimum_layer_thickness_m,
            },
        };
        try runScienceCells(context, &microbial_mixing_context, ecosys.surface_topsoil_microbial_mixing.applyTile);
    }
    if (context.chemistry_reaction_parameters.*) |reaction_parameters| if (context.organic_parameters.*) |*organic_runtime_parameters| {
        var fertilizer_context: ecosys.surface_litter_fertilizer_step.ApplyContext = .{
            .state = context.surface_litter_fertilizer,
            .chemistry = context.surface_litter_chemistry,
            .surface_organic = context.surface_organic,
            .litter_water_m3 = context.surface_precipitation.litter_water_m3,
            .dry_litter_volume_m3 = context.surface_litter_geometry.dry_litter_volume_m3,
            .biologically_active_water_m3 = context.surface_microbial_environment.biologically_active_water_m3,
            .active_biomass_respiration_g_c_per_step = context.surface_microbial_oxygen.oxygen_limited_activity_g_c_per_step,
            .microbial_temperature_factor = context.surface_microbial_environment.growth_temperature_response,
            .litter_dry_mass_Mg_per_g_c = organic_runtime_parameters.surface_litter_dry_mass_Mg_per_g_c,
            .step_duration_h = 1,
            .parameters = reaction_parameters.surface_fertilizer,
            .diagnostics = context.surface_litter_fertilizer_diagnostics,
        };
        try runScienceCells(context, &fertilizer_context, ecosys.surface_litter_fertilizer_step.applyTile);
    };
}

noinline fn runSurfaceBiogeochemistryBySerialTile(context: anytype, surface_parameters: anytype) !void {
    const plan = context.tile_plan;
    for (plan.tiles, 0..) |_, tile_index| {
        const owned_cells = try plan.ownedCells(tile_index);
        context.active_tile_cells.* = owned_cells;
        defer context.active_tile_cells.* = null;
        try runSurfaceBiogeochemistryBatch(context, surface_parameters);
    }
}

fn commitHourlySoilSoluteContributionGenerations(
    context: anytype,
) !void {
    const generation = std.math.cast(
        u64,
        context.executed_weather_hours.*,
    ) orelse return error.LateralContributionGenerationExceedsU64;
    const workspace = context.lateral_contribution_workspace.*;
    const buffer_byte_count: usize = 64 * 1024;
    const aqueous_species_count =
        context.micropore_solute_state.species_count;
    const micropore_store = try workspace.store(
        context.allocator,
        workspace.micropore_solute,
        buffer_byte_count,
        generation,
    );
    try ecosys.hourly_lateral_contribution_io.publishTransportGeneration(
        context.allocator,
        micropore_store,
        context.tile_plan.*,
        context.grid.soil_layer_capacity,
        context.soil_transport_faces.micropore_faces,
        aqueous_species_count,
        context.micropore_solute_face_flux_mol,
    );
    try ecosys.hourly_lateral_contribution_io.verifyGenerationConservation(
        context.allocator,
        micropore_store,
        context.tile_plan.*,
        context.grid.soil_layer_capacity,
        aqueous_species_count,
    );
    try ecosys.hourly_lateral_contribution_io.commitTransportGeneration(
        context.allocator,
        micropore_store,
        context.tile_plan.*,
        context.grid.soil_layer_capacity,
        aqueous_species_count,
        context.micropore_solute_state.amount_mol,
    );

    const macropore_store = try workspace.store(
        context.allocator,
        workspace.macropore_solute,
        buffer_byte_count,
        generation,
    );
    try ecosys.hourly_lateral_contribution_io.publishTransportGeneration(
        context.allocator,
        macropore_store,
        context.tile_plan.*,
        context.grid.soil_layer_capacity,
        context.soil_transport_faces.macropore_faces,
        aqueous_species_count,
        context.macropore_solute_face_flux_mol,
    );
    try ecosys.hourly_lateral_contribution_io.verifyGenerationConservation(
        context.allocator,
        macropore_store,
        context.tile_plan.*,
        context.grid.soil_layer_capacity,
        aqueous_species_count,
    );
    try ecosys.hourly_lateral_contribution_io.commitTransportGeneration(
        context.allocator,
        macropore_store,
        context.tile_plan.*,
        context.grid.soil_layer_capacity,
        aqueous_species_count,
        context.macropore_solute_state.amount_mol,
    );
}

fn commitHourlyWaterHeatStateGeneration(
    context: anytype,
    accepted: *const ecosys.soil_water_heat_step.DeferredMappedResult,
) !void {
    const generation = std.math.cast(
        u64,
        context.executed_weather_hours.*,
    ) orelse return error.LateralContributionGenerationExceedsU64;
    const carrier_count =
        ecosys.soil_water_heat_step.deferred_grid_carrier_count;
    const component_count = try std.math.mul(
        usize,
        context.grid.layer_count,
        carrier_count,
    );
    if (accepted.grid_delta_by_layer_carrier.len != component_count)
        return error.DeferredSoilStateDeltaDimensionMismatch;
    const packed_state = try context.allocator.alloc(f64, component_count);
    defer context.allocator.free(packed_state);
    const state_fields: [carrier_count][]f64 = .{
        context.grid.matrix_liquid_water_m3,
        context.grid.macropore_liquid_water_m3,
        context.grid.liquid_water_m3,
        context.grid.matrix_air_volume_m3,
        context.grid.macropore_air_volume_m3,
        context.grid.air_volume_m3,
        context.grid.water_vapor_volume_m3,
        context.grid.matrix_ice_water_m3,
        context.grid.macropore_ice_water_m3,
        context.grid.ice_water_m3,
        context.grid.soil_temperature_k,
        context.grid.matric_potential_mpa,
    };
    for (0..context.grid.layer_count) |layer| {
        for (0..carrier_count) |carrier|
            packed_state[layer * carrier_count + carrier] =
                state_fields[carrier][layer];
    }
    const workspace = context.lateral_contribution_workspace.*;
    const store = try workspace.store(
        context.allocator,
        workspace.water_heat_vapor,
        64 * 1024,
        generation,
    );
    try ecosys.hourly_lateral_contribution_io
        .publishLayerCellDeltaGeneration(
        context.allocator,
        store,
        context.tile_plan.*,
        context.grid.soil_layer_capacity,
        carrier_count,
        accepted.grid_delta_by_layer_carrier,
    );
    try ecosys.hourly_lateral_contribution_io.commitFiniteStateGeneration(
        context.allocator,
        store,
        context.tile_plan.*,
        context.grid.soil_layer_capacity,
        carrier_count,
        packed_state,
    );
    for (0..context.grid.layer_count) |layer| {
        for (0..carrier_count) |carrier| {
            const value = packed_state[layer * carrier_count + carrier];
            if (!std.math.isFinite(value) or
                (carrier < 10 and value < 0) or
                (carrier == 10 and value <= 0))
                return error.InvalidCommittedSoilWaterHeatState;
        }
    }
    for (0..context.grid.layer_count) |layer| {
        for (0..carrier_count) |carrier| {
            const value = packed_state[layer * carrier_count + carrier];
            state_fields[carrier][layer] = value;
        }
    }
    @memcpy(
        context.transport_hydrology.micropore_water_volume_m3,
        context.grid.matrix_liquid_water_m3,
    );
    @memcpy(
        context.transport_hydrology.macropore_water_volume_m3,
        context.grid.macropore_liquid_water_m3,
    );
    @memcpy(
        context.transport_hydrology.matrix_air_volume_m3,
        context.grid.matrix_air_volume_m3,
    );
    @memcpy(
        context.transport_hydrology.macropore_air_volume_m3,
        context.grid.macropore_air_volume_m3,
    );
    @memcpy(
        context.transport_hydrology.air_volume_m3,
        context.grid.air_volume_m3,
    );
    @memcpy(
        context.transport_hydrology.water_vapor_volume_m3,
        context.grid.water_vapor_volume_m3,
    );
}

fn commitHourlyGasContributionGeneration(
    context: anytype,
    accepted: *const ecosys.gas_transport.State,
) !void {
    const generation = std.math.cast(
        u64,
        context.executed_weather_hours.*,
    ) orelse return error.LateralContributionGenerationExceedsU64;
    const species_count = ecosys.gas_transport.species_count;
    const phase_count: usize = 4;
    const carrier_count = phase_count * species_count;
    const component_count = try std.math.mul(
        usize,
        context.grid.layer_count,
        carrier_count,
    );
    const packed_amount = try context.allocator.alloc(f64, component_count);
    defer context.allocator.free(packed_amount);
    const delta = try context.allocator.alloc(f64, component_count);
    defer context.allocator.free(delta);
    const current = context.gas_transport;
    const current_phases: [phase_count][]f64 = .{
        current.gaseous_mass_g,
        current.dissolved_mass_g,
        current.macropore_dissolved_mass_g,
        current.band_dissolved_mass_g,
    };
    const accepted_phases: [phase_count][]const f64 = .{
        accepted.gaseous_mass_g,
        accepted.dissolved_mass_g,
        accepted.macropore_dissolved_mass_g,
        accepted.band_dissolved_mass_g,
    };
    for (0..context.grid.layer_count) |layer| {
        for (0..phase_count) |phase| {
            for (0..species_count) |species| {
                const phase_component = layer * species_count + species;
                const packed_component =
                    layer * carrier_count + phase * species_count + species;
                packed_amount[packed_component] =
                    current_phases[phase][phase_component];
                delta[packed_component] =
                    accepted_phases[phase][phase_component] -
                    current_phases[phase][phase_component];
                if (!std.math.isFinite(delta[packed_component]))
                    return error.NonFiniteAcceptedGasDelta;
            }
        }
    }
    const workspace = context.lateral_contribution_workspace.*;
    const store = try workspace.store(
        context.allocator,
        workspace.gas,
        64 * 1024,
        generation,
    );
    try ecosys.hourly_lateral_contribution_io
        .publishLayerCellDeltaGeneration(
        context.allocator,
        store,
        context.tile_plan.*,
        context.grid.soil_layer_capacity,
        carrier_count,
        delta,
    );
    try ecosys.hourly_lateral_contribution_io.commitTransportGeneration(
        context.allocator,
        store,
        context.tile_plan.*,
        context.grid.soil_layer_capacity,
        carrier_count,
        packed_amount,
    );
    for (0..context.grid.layer_count) |layer| {
        for (0..phase_count) |phase| {
            for (0..species_count) |species| {
                const phase_component = layer * species_count + species;
                const packed_component =
                    layer * carrier_count + phase * species_count + species;
                if (!std.math.isFinite(packed_amount[packed_component]) or
                    packed_amount[packed_component] < 0)
                    return error.InvalidCommittedGasInventory;
                current_phases[phase][phase_component] =
                    packed_amount[packed_component];
            }
        }
    }
}

fn executeHourlyScience(
    context: anytype,
    hour_of_day: u8,
    radiation_by_cell: []const ecosys.atmospheric_radiation.Result,
    forcing_by_cell: []const ecosys.weather.HourlyForcing,
    weather_header_by_cell: []const ecosys.weather.Header,
    plant_calendar_by_cell: []const ecosys.plant_development.Calendar,
    fertilizer_band_hour: ecosys.fertilizer_band_phase_coordinator.HourToken,
    gas_failure_report: ?ecosys.soil_gas_transport_step.FailureReportRequest,
    solute_failure_report: ?ecosys.solute_failure_reporter.Request,
) !void {
    if (radiation_by_cell.len != context.grid.cell_count or
        forcing_by_cell.len != context.grid.cell_count or
        weather_header_by_cell.len != context.grid.cell_count or
        plant_calendar_by_cell.len != context.grid.cell_count)
        return error.HourlyEnvironmentDimensionMismatch;
    const plant_calendar = plant_calendar_by_cell[0];
    if (hour_of_day > 23) return error.InvalidHourlyScienceHour;
    context.plant_available_nutrients.resetHourlyChanges();
    try ecosys.soil_organic_carbon_change.captureHourStart(context.soil_organic, context.soil_organic_carbon_at_hour_start_g_c);
    // Derive restart-sensitive surface state from the checkpointed snow owner.
    for (0..context.grid.cell_count) |cell| context.snow_depth_m[cell] = context.snow_transport.cumulative_depth_m[(cell + 1) * context.snow_transport.layer_capacity - 1];
    // Chemistry is concentration-based while TRNSFRS is amount-based. Export
    // before water moves so dilution does not create or destroy solute mass.
    try ecosys.soil_aqueous_transport_bridge.exportChemistry(context.soil_chemistry, context.micropore_solute_state);
    // Legacy TRNSFR copies extensive ZNH4S/ZNH4B inventories before water
    // transport. Capture mineral N at the same boundary; reconstructing
    // amount from concentration after water moves manufactures or destroys N.
    try context.mineral_nitrogen_transport.captureHourStartMatrix(
        context.soil_chemistry,
        context.soil_reactive_nitrogen,
        context.grid.matrix_liquid_water_m3,
        context.mineral_nitrogen_zone_fractions,
        context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
    );
    var forcing_context = ecosys.atmospheric_forcing.MappedApplyContext{ .state = context.atmosphere, .forcing_by_cell = forcing_by_cell };
    try runKernelAcrossSerialTiles(context, &forcing_context, ecosys.atmospheric_forcing.applyMappedTile);
    for (0..context.grid.cell_count) |cell| {
        const irrigation_m = context.irrigation_water_depth_m[cell];
        if (!std.math.isFinite(irrigation_m) or irrigation_m < 0) return error.InvalidHourlyIrrigationDepth;
        context.atmosphere.rainfall_m[cell] += irrigation_m;
        context.atmosphere.precipitation_m[cell] += irrigation_m;
    }
    @memset(context.surface_total_canopy_area_m2, 0);
    @memset(context.surface_canopy_height_m, 0);
    if (context.canopy_layer_distribution.*) |*layers| for (0..context.grid.cell_count) |cell| {
        const first = cell * layers.layer_count;
        for (0..layers.layer_count) |layer| context.surface_total_canopy_area_m2[cell] += layers.cell_leaf_area_m2[first + layer] + layers.cell_stalk_area_m2[first + layer] + layers.cell_standing_dead_area_m2[first + layer];
        const boundaries = try layers.cellBoundaries(cell);
        context.surface_canopy_height_m[cell] = @max(0, boundaries[boundaries.len - 1]);
    };
    // HOUR1 ARLFS/AREA refresh: canopy structure must follow the live
    // layer-distributed leaf inventory each hour, not its initialization value.
    if (context.canopy_layer_distribution.*) |*layers| if (context.detailed_canopy.*) |*canopy| if (context.canopy_structure.*) |*structure| {
        try layers.publishLeafAreaIndex(canopy, context.plants, context.canopy_cell_area_m2);
        var structure_context: ecosys.canopy_structure.ApplyContext = .{ .structure = structure, .plants = context.plants };
        try runKernelAcrossSerialTiles(context, &structure_context, ecosys.canopy_structure.applyLeafAreaTile);
    };
    for (weather_header_by_cell, 0..) |header, cell| context.hourly_weather_reference_height_m[cell] = header.aerodynamic_roughness_m;
    var surface_aerodynamic_context: ecosys.surface_aerodynamics.ApplyContext = .{ .state = context.surface_aerodynamics, .cell_area_m2 = context.canopy_cell_area_m2, .total_canopy_area_m2 = context.surface_total_canopy_area_m2, .canopy_height_m = context.surface_canopy_height_m, .snow_depth_m = context.snow_depth_m, .weather_reference_height_m = context.hourly_weather_reference_height_m, .wind_speed_m_per_h = context.atmosphere.wind_speed_m_per_h, .parameters = context.runscript.surface_aerodynamic_parameters };
    try runKernelAcrossSerialTiles(context, &surface_aerodynamic_context, ecosys.surface_aerodynamics.applyTile);
    try context.ground_air.refreshGeometry(context.canopy_cell_area_m2, context.surface_aerodynamics.wind_reference_height_m, context.runscript.ground_air_parameters);
    for (0..context.grid.cell_count) |cell| context.ground_air_vapor_pressure_kpa[cell] = try ecosys.ground_air_exchange.vaporPressureKpa(context.ground_air.vapor_volume_fraction[cell], context.ground_air.temperature_k[cell], context.runscript.ground_air_parameters);
    for (radiation_by_cell, forcing_by_cell, 0..) |cell_radiation, cell_forcing, cell| {
        context.hourly_extraterrestrial_shortwave_mj_per_m2[cell] =
            cell_radiation.extraterrestrial_shortwave_mj_per_m2;
        context.hourly_solar_angle_sine[cell] = cell_radiation.solar_angle_sine;
        context.hourly_solar_azimuth_radians[cell] = cell_radiation.solar_azimuth_radians;
        context.hourly_adjusted_shortwave_mj_per_m2[cell] =
            cell_forcing.shortwave_radiation_mj_per_m2;
    }
    var canopy_radiation_context: ecosys.canopy_radiation.MappedApplyContext = .{ .state = context.canopy_radiation, .horizontal_shortwave_mj_per_m2 = context.hourly_adjusted_shortwave_mj_per_m2, .extraterrestrial_horizontal_shortwave_mj_per_m2 = context.hourly_extraterrestrial_shortwave_mj_per_m2, .solar_angle_sine = context.hourly_solar_angle_sine };
    try runKernelAcrossSerialTiles(context, &canopy_radiation_context, ecosys.canopy_radiation.applyMappedTile);
    if (context.canopy_optics.*) |*optics| {
        var optics_context: ecosys.canopy_optics.ApplyContext = .{ .state = optics, .radiation = context.canopy_radiation };
        try runKernelAcrossSerialTiles(context, &optics_context, ecosys.canopy_optics.applyLeafAbsorptionTile);
    }
    try context.canopy_geometry.directSolarIncidenceMapped(
        context.hourly_solar_angle_sine,
        context.direct_incidence_fraction,
        context.direct_incidence_per_horizontal_area,
        context.direct_scattering_direction,
    );
    if (context.canopy_interception.*) |*interception| {
        if (context.canopy_layer_distribution.* != null) try ecosys.canopy_interception.refreshLayerTransmission(interception, &context.canopy_layer_distribution.*.?, &context.canopy_structure.*.?, context.canopy_geometry, context.direct_incidence_per_horizontal_area, context.canopy_cell_area_m2);
        var interception_context: ecosys.canopy_interception.ApplyContext = .{ .result = interception, .structure = &context.canopy_structure.*.?, .optics = &context.canopy_optics.*.?, .geometry = context.canopy_geometry, .direct_incidence_fraction = context.direct_incidence_fraction, .direct_incidence_per_horizontal_area = context.direct_incidence_per_horizontal_area };
        try runKernelAcrossSerialTiles(context, &interception_context, ecosys.canopy_interception.applySingleLayerTile);
        if (context.canopy_layer_distribution.* != null) try ecosys.canopy_interception.refreshAtmosphericLayerAbsorption(interception, &context.canopy_layer_distribution.*.?, &context.canopy_structure.*.?, &context.canopy_optics.*.?, context.canopy_geometry, context.canopy_radiation, context.direct_incidence_fraction, context.direct_incidence_per_horizontal_area, context.direct_scattering_direction, context.canopy_cell_area_m2, context.runscript.woody_optics_parameters);
    }
    var terrain_radiation_context: ecosys.terrain_radiation.MappedDirectSolarContext = .{ .state = context.terrain_radiation, .solar_angle_sine = context.hourly_solar_angle_sine, .solar_azimuth_radians = context.hourly_solar_azimuth_radians };
    try runKernelAcrossSerialTiles(context, &terrain_radiation_context, ecosys.terrain_radiation.applyMappedDirectSolarTile);
    var ground_context: ecosys.ground_radiation.ApplyContext = .{ .result = context.ground_radiation, .radiation = context.canopy_radiation, .interception = if (context.canopy_interception.*) |*interception| interception else null, .terrain = context.terrain_radiation };
    try runKernelAcrossSerialTiles(context, &ground_context, ecosys.ground_radiation.applyTile);
    if (context.canopy_interception.* != null and context.canopy_layer_distribution.* != null) try ecosys.canopy_interception.applyGroundReflectedUpwardSweep(&context.canopy_interception.*.?, &context.canopy_layer_distribution.*.?, &context.canopy_structure.*.?, &context.canopy_optics.*.?, context.canopy_geometry, context.ground_radiation.reflected_shortwave_mj_per_m2, context.ground_radiation.reflected_par_micromol_per_m2_per_s, context.runscript.woody_optics_parameters);
    if (context.canopy_precipitation_retention.*) |*retention| try ecosys.canopy_precipitation_retention.refreshFromModel(retention, &context.canopy_layer_distribution.*.?, &context.detailed_canopy.*.?, &context.canopy_interception.*.?, context.atmosphere.rainfall_m, context.canopy_cell_area_m2, context.canopy_layer_controls.root_profile_type, context.hourly_solar_angle_sine, context.ground_radiation.incident_shortwave_mj_per_m2, context.runscript.canopy_retention_parameters);
    try ecosys.surface_precipitation.prepareFromModel(context.surface_precipitation, context.atmosphere, context.grid, if (context.canopy_precipitation_retention.*) |*retention| retention else null, context.canopy_cell_area_m2, context.snow_depth_m, context.runscript.snow_full_cover_depth_m, context.iteration_limits.water_heat_solute_max_iterations);
    try context.snow_transport.commitAtmosphericWater(context.surface_precipitation.snow_to_snow_m3_per_h, context.surface_precipitation.rain_to_snow_m3_per_h, context.surface_precipitation.heat_to_snow_mj_per_h, context.atmosphere.air_temperature_k, context.runscript.initial_snow_density_Mg_per_m3);
    _ = try ecosys.snow_heat_conduction.solve(context.allocator, context.snow_transport, context.runscript.snow_heat_conduction_parameters, .{
        .timestep_h = 1,
        .full_snow_cover_depth_m = context.runscript.snow_full_cover_depth_m,
        .absolute_tolerance_k = context.config.absolute_tolerance,
        .relative_tolerance = context.config.relative_tolerance,
        .picard_relaxation = context.config.picard_relaxation,
        .max_iterations = context.iteration_limits.snowpack_max_iterations,
    });
    _ = try ecosys.snow_vapor_diffusion.solve(context.allocator, context.snow_transport, context.runscript.snow_vapor_diffusion_parameters, .{
        .timestep_h = 1,
        .full_snow_cover_depth_m = context.runscript.snow_full_cover_depth_m,
        .absolute_tolerance_m3 = context.config.absolute_tolerance,
        .relative_tolerance = context.config.relative_tolerance,
        .max_iterations = context.iteration_limits.snowpack_max_iterations,
    });
    const snow_vapor_equilibrium_report = try ecosys.snow_vapor_equilibrium.solve(context.allocator, context.snow_transport, context.runscript.snow_vapor_parameters, .{
        .absolute_tolerance_m3 = context.config.absolute_tolerance,
        .relative_tolerance = context.config.relative_tolerance,
        .max_iterations = context.iteration_limits.snowpack_max_iterations,
    });
    const snow_phase_change_report = try ecosys.snow_phase_change.solve(context.allocator, context.snow_transport, .{
        .ice_density_Mg_per_m3 = context.runscript.snow_ice_density_Mg_per_m3,
        .latent_heat_of_fusion_mj_per_m3 = context.runscript.snow_latent_heat_of_fusion_mj_per_m3,
        .damping_divisor = context.runscript.snow_phase_damping_divisor,
        .absolute_temperature_tolerance_k = context.config.absolute_tolerance,
        .relative_tolerance = context.config.relative_tolerance,
        .picard_relaxation = context.config.picard_relaxation,
        .max_iterations = context.iteration_limits.snowpack_max_iterations,
    });
    try ecosys.snow_melt_water_routing.calculate(.{
        .cell_count = context.grid.cell_count,
        .layer_capacity = context.snow_transport.layer_capacity,
        .active = context.snow_transport.active,
        .liquid_water_volume_m3 = context.snow_transport.liquid_water_volume_m3,
        .solid_snow_volume_m3 = context.snow_transport.solid_snow_water_equivalent_m3,
        .air_filled_volume_m3 = context.snow_transport.air_filled_volume_m3,
        .litter_cover_fraction = context.surface_precipitation.litter_cover_fraction,
        .micropore_fraction = context.surface_precipitation.matrix_fraction,
        .macropore_fraction = context.surface_precipitation.macropore_fraction,
        .topsoil_micropore_air_capacity_m3 = context.surface_precipitation.matrix_air_capacity_m3,
        .topsoil_macropore_air_capacity_m3 = context.surface_precipitation.macropore_air_capacity_m3,
        .other_micropore_water_input_m3 = context.surface_precipitation.water_to_matrix_m3_per_h,
        .other_macropore_water_input_m3 = context.surface_precipitation.water_to_macropore_m3_per_h,
        .step_fraction = 1,
    }, .{
        .downward_water_flux_m3 = context.transport_hydrology.snow_downward_water_flux_m3_per_step,
        .litter_water_flux_m3 = context.transport_hydrology.snow_to_litter_water_flux_m3_per_step,
        .soil_micropore_water_flux_m3 = context.transport_hydrology.snow_to_soil_micropore_flux_m3_per_step,
        .soil_macropore_water_flux_m3 = context.transport_hydrology.snow_to_soil_macropore_flux_m3_per_step,
    });
    @memcpy(context.transport_hydrology.snow_liquid_water_volume_m3, context.snow_transport.liquid_water_volume_m3);
    @memset(context.snow_atmospheric_input_g, 0);
    for (0..context.grid.cell_count) |cell| {
        const weather_header = weather_header_by_cell[cell];
        var rain_gas_concentration = [_]f64{0} ** 5;
        if (context.surface_gas_parameters.*) |surface_parameters| {
            const solubility = try ecosys.gas_transport.surfaceSolubilityWaterToAir(context.atmosphere.air_temperature_k[cell], surface_parameters.solubility);
            for (0..5) |species| rain_gas_concentration[species] = context.current_atmospheric_gas_concentration_g_per_m3.*[species] * solubility[species];
        }
        const nutrients = if (context.chemistry_reaction_parameters.*) |reaction_parameters|
            try ecosys.precipitation_nutrient_speciation.calculate(.{ .ph = weather_header.precipitation_ph, .ammonium_g_n_per_m3 = weather_header.precipitation_ammonium_g_per_m3, .nitrate_g_n_per_m3 = weather_header.precipitation_nitrate_g_per_m3, .phosphate_g_p_per_m3 = weather_header.precipitation_phosphate_g_per_m3 }, reaction_parameters.aqueous_constants, reaction_parameters.phosphate_constants)
        else
            [5]f64{ weather_header.precipitation_ammonium_g_per_m3 / 14, 0, weather_header.precipitation_nitrate_g_per_m3 / 14, 0, weather_header.precipitation_phosphate_g_per_m3 / 31 };
        const area_m2 = context.canopy_cell_area_m2[cell];
        const irrigation_depth_m = context.irrigation_water_depth_m[cell];
        const irrigation_volume_m3 = irrigation_depth_m * area_m2;
        const total_liquid_volume_m3 = context.surface_precipitation.rainfall_m3_per_h[cell];
        const liquid_to_snow_m3 = context.surface_precipitation.rain_to_snow_m3_per_h[cell];
        const irrigation_to_snow_m3 = if (total_liquid_volume_m3 > 0) liquid_to_snow_m3 * irrigation_volume_m3 / total_liquid_volume_m3 else 0;
        const rain_to_snow_m3 = liquid_to_snow_m3 - irrigation_to_snow_m3 + context.surface_precipitation.snow_to_snow_m3_per_h[cell];
        var irrigation_nutrients = [_]f64{0} ** 5;
        var irrigation_ions_g_per_m3 = [_]f64{0} ** 8;
        if (irrigation_depth_m > 0) {
            const first = cell * ecosys.irrigation_management_dispatch.dissolved_species_count;
            const concentration = context.irrigation_dissolved_mass_g_per_m2[first .. first + ecosys.irrigation_management_dispatch.dissolved_species_count];
            const hydrogen_mol_per_m3 = context.irrigation_hydrogen_mol_per_m2[cell] / irrigation_depth_m;
            const irrigation_ph = -std.math.log10(@max(hydrogen_mol_per_m3 / 1000.0, 1.0e-14));
            irrigation_nutrients = if (context.chemistry_reaction_parameters.*) |reaction_parameters|
                try ecosys.precipitation_nutrient_speciation.calculate(.{
                    .ph = irrigation_ph,
                    .ammonium_g_n_per_m3 = concentration[0] / irrigation_depth_m,
                    .nitrate_g_n_per_m3 = concentration[1] / irrigation_depth_m,
                    .phosphate_g_p_per_m3 = concentration[2] / irrigation_depth_m,
                }, reaction_parameters.aqueous_constants, reaction_parameters.phosphate_constants)
            else
                [5]f64{ concentration[0] / irrigation_depth_m / 14, 0, concentration[1] / irrigation_depth_m / 14, 0, concentration[2] / irrigation_depth_m / 31 };
            for (&irrigation_ions_g_per_m3, 0..) |*ion, index| ion.* = concentration[3 + index] / irrigation_depth_m;
        }
        const input = try ecosys.snow_solute_transport.atmosphericInputG(rain_to_snow_m3, irrigation_to_snow_m3, rain_gas_concentration, [_]f64{0} ** 5, nutrients, irrigation_nutrients, [_]f64{0} ** 8, irrigation_ions_g_per_m3);
        @memcpy(context.snow_atmospheric_input_g[cell * ecosys.snow_solute_transport.species_count ..][0..ecosys.snow_solute_transport.species_count], &input);
        context.direct_surface_solute_input[cell] = .{};
        const top_soil_zone_fractions =
            try context.fertilizer_band.zoneFractions(cell, 0);
        const direct_liquid_m3 = context.surface_precipitation.water_to_litter_m3_per_h[cell] + context.surface_precipitation.water_to_matrix_m3_per_h[cell] + context.surface_precipitation.water_to_macropore_m3_per_h[cell];
        if (direct_liquid_m3 > 0) {
            const total_weather_rain_m3 = @max(0, total_liquid_volume_m3 - irrigation_volume_m3);
            const snow_fraction = if (total_liquid_volume_m3 > 0) liquid_to_snow_m3 / total_liquid_volume_m3 else 0;
            const direct_rain_m3 = total_weather_rain_m3 * (1 - snow_fraction);
            const direct_irrigation_m3 = irrigation_volume_m3 * (1 - snow_fraction);
            const direct_input = try ecosys.snow_solute_transport.atmosphericInputG(direct_rain_m3, direct_irrigation_m3, rain_gas_concentration, [_]f64{0} ** 5, nutrients, irrigation_nutrients, [_]f64{0} ** 8, irrigation_ions_g_per_m3);
            const litter_fraction = std.math.clamp(context.surface_precipitation.water_to_litter_m3_per_h[cell] / direct_liquid_m3, 0, 1);
            const soil_fraction = 1 - litter_fraction;
            for (direct_input, 0..) |amount_g, species| {
                context.direct_surface_solute_input[cell].litter_g[species] = amount_g * litter_fraction;
                const band_fraction = switch (@as(ecosys.snow_solute_transport.Species, @enumFromInt(species))) {
                    .ammonium_nitrogen, .ammonia_nitrogen => top_soil_zone_fractions.ammonium_band,
                    .nitrate_nitrogen => top_soil_zone_fractions.nitrate_band,
                    .hydrogen_phosphate_phosphorus, .dihydrogen_phosphate_phosphorus => top_soil_zone_fractions.phosphate_band,
                    else => 0,
                };
                context.direct_surface_solute_input[cell].soil_nonband_g[species] = amount_g * soil_fraction * (1 - band_fraction);
                context.direct_surface_solute_input[cell].soil_band_g[species] = amount_g * soil_fraction * band_fraction;
            }
        }
        context.snow_surface_partitions[cell] = .{
            .litter_cover_fraction = context.surface_precipitation.litter_cover_fraction[cell],
            .bare_soil_fraction = 1 - context.surface_precipitation.litter_cover_fraction[cell],
            .nonband_ammonium_fraction = top_soil_zone_fractions.ammonium_non_band,
            .band_ammonium_fraction = top_soil_zone_fractions.ammonium_band,
            .nonband_nitrate_fraction = top_soil_zone_fractions.nitrate_non_band,
            .band_nitrate_fraction = top_soil_zone_fractions.nitrate_band,
            .nonband_phosphate_fraction = top_soil_zone_fractions.phosphate_non_band,
            .band_phosphate_fraction = top_soil_zone_fractions.phosphate_band,
        };
    }
    try context.atmospheric_solute_input_ledger.accumulateAcceptedHour(
        context.snow_atmospheric_input_g,
        context.direct_surface_solute_input,
    );
    _ = try ecosys.snow_transport_solver.solve(context.allocator, context.snow_transport, .{
        .atmospheric_top_input_g = context.snow_atmospheric_input_g,
        .transport_water_volume_m3 = context.transport_hydrology.snow_liquid_water_volume_m3,
        .water_flux_to_lower_m3 = context.transport_hydrology.snow_downward_water_flux_m3_per_step,
        .litter_water_flux_m3 = context.transport_hydrology.snow_to_litter_water_flux_m3_per_step,
        .soil_micropore_water_flux_m3 = context.transport_hydrology.snow_to_soil_micropore_flux_m3_per_step,
        .soil_macropore_water_flux_m3 = context.transport_hydrology.snow_to_soil_macropore_flux_m3_per_step,
        .surface_partitions = context.snow_surface_partitions,
    }, .{ .absolute_tolerance_g = context.config.absolute_tolerance, .relative_tolerance = context.config.relative_tolerance, .picard_relaxation = context.config.picard_relaxation, .max_iterations = context.iteration_limits.snowpack_max_iterations }, context.snow_surface_discharge);
    try context.snow_transport.commitMeltWater(context.transport_hydrology.snow_downward_water_flux_m3_per_step, context.transport_hydrology.snow_to_litter_water_flux_m3_per_step, context.transport_hydrology.snow_to_soil_micropore_flux_m3_per_step, context.transport_hydrology.snow_to_soil_macropore_flux_m3_per_step);
    try ecosys.snow_compaction.apply(context.allocator, context.snow_transport, .{
        .snowfall_water_equivalent_m3 = context.surface_precipitation.snow_to_snow_m3_per_h,
        .atmospheric_temperature_k = context.atmosphere.air_temperature_k,
        .timestep_h = 1,
        .initial_snow_density_Mg_per_m3 = context.runscript.initial_snow_density_Mg_per_m3,
        .ice_density_Mg_per_m3 = context.runscript.snow_ice_density_Mg_per_m3,
    }, context.runscript.snow_compaction_parameters);
    _ = try ecosys.snow_relayering.apply(context.allocator, context.snow_transport, context.config.absolute_tolerance);
    for (0..context.grid.cell_count) |cell| {
        context.snow_depth_m[cell] = context.snow_transport.cumulative_depth_m[(cell + 1) * context.snow_transport.layer_capacity - 1];
        var solid_snow_water_equivalent_m3: f64 = 0;
        for (0..context.snow_transport.layer_capacity) |layer| solid_snow_water_equivalent_m3 += context.snow_transport.solid_snow_water_equivalent_m3[cell * context.snow_transport.layer_capacity + layer];
        context.surface_precipitation.solid_snow_water_equivalent_m3[cell] = solid_snow_water_equivalent_m3;
    }
    for (0..context.grid.cell_count) |cell| {
        context.surface_precipitation.water_to_litter_m3_per_h[cell] += context.transport_hydrology.snow_to_litter_water_flux_m3_per_step[cell];
        context.surface_precipitation.water_to_matrix_m3_per_h[cell] += context.transport_hydrology.snow_to_soil_micropore_flux_m3_per_step[cell];
        context.surface_precipitation.water_to_macropore_m3_per_h[cell] += context.transport_hydrology.snow_to_soil_macropore_flux_m3_per_step[cell];
    }
    if (context.canopy_exposure.*) |*exposure| {
        var exposure_context: ecosys.canopy_exposure.ApplyContext = .{ .result = exposure, .structure = &context.canopy_structure.*.?, .interception = &context.canopy_interception.*.?, .ground_radiation = context.ground_radiation, .solar_angle_sine_by_cell = context.hourly_solar_angle_sine };
        try runKernelAcrossSerialTiles(context, &exposure_context, ecosys.canopy_exposure.applyTile);
    }
    var surface_energy_context: ecosys.surface_energy.ApplyContext = .{ .result = context.surface_energy, .grid = context.grid, .atmosphere = context.atmosphere, .ground_radiation = context.ground_radiation, .snow_depth_m = context.snow_depth_m, .exposure = if (context.canopy_exposure.*) |*exposure| exposure else null, .settings = context.surface_energy_settings };
    try runKernelAcrossSerialTiles(context, &surface_energy_context, ecosys.surface_energy.applyTile);
    for (0..context.grid.cell_count) |cell| {
        const area_m2 = context.canopy_cell_area_m2[cell];
        if (!std.math.isFinite(area_m2) or area_m2 <= 0) return error.InvalidSurfaceFireCellArea;
        context.surface_combustion_heat_mj_per_m2[cell] = context.delayed_surface_combustion_heat_mj[cell] / area_m2;
        if (!std.math.isFinite(context.surface_combustion_heat_mj_per_m2[cell]) or context.surface_combustion_heat_mj_per_m2[cell] < 0) return error.NonFiniteSurfaceCombustionHeatSource;
    }
    var surface_temperature_context: ecosys.surface_temperature_solver.ApplyContext = .{ .result = context.surface_temperature, .grid = context.grid, .atmosphere = context.atmosphere, .air_temperature_k = context.ground_air.temperature_k, .air_vapor_pressure_kpa = context.ground_air_vapor_pressure_kpa, .ground_radiation = context.ground_radiation, .surface_energy = context.surface_energy, .soil_thermal = context.soil_thermal, .exposure = if (context.canopy_exposure.*) |*exposure| exposure else null, .external_heat_mj_per_m2 = context.surface_combustion_heat_mj_per_m2, .surface_phase = .{ .liquid_water_m3 = context.surface_precipitation.litter_water_m3, .ice_water_equivalent_m3 = context.surface_litter_ice_m3, .retention_capacity_m3 = context.surface_precipitation.litter_water_capacity_m3, .horizontal_area_m2 = context.canopy_cell_area_m2, .residual_water_content_m3_per_m3 = context.runscript.soil_process_parameters.surface_residue_residual_water_content_m3_per_m3, .van_genuchten_alpha_per_m = context.runscript.soil_process_parameters.surface_residue_van_genuchten_alpha_per_m, .van_genuchten_n = context.runscript.soil_process_parameters.surface_residue_van_genuchten_n, .gravitational_water_potential_mpa_per_m = context.runscript.soil_process_parameters.gravitational_water_potential_mpa_per_m, .latent_heat_of_fusion_mj_per_m3 = context.runscript.soil_phase_heat_parameters.freeze_thaw.latent_heat_of_fusion_mj_per_m3, .pure_water_melting_temperature_k = context.runscript.soil_phase_heat_parameters.freeze_thaw.pure_water_freezing_temperature_k }, .settings = .{ .timestep_hours = 1.0, .sensible_heat_conductance_mj_per_m2_h_k = context.runscript.surface_sensible_heat_conductance_mj_per_m2_h_k, .latent_heat_conductance_mj_per_m2_h_kpa = context.runscript.surface_latent_heat_conductance_mj_per_m2_h_kpa, .liquid_water_heat_capacity_mj_per_m3_k = context.runscript.canopy_surface_exchange_parameters.liquid_water_heat_capacity_mj_per_m3_k, .latent_heat_of_vaporization_mj_per_m3 = context.runscript.ground_air_parameters.liquid_water_latent_heat_mj_per_m3, .surface_vapor_activity_fraction = context.runscript.surface_vapor_activity_fraction, .minimum_temperature_k = context.runscript.minimum_surface_temperature_k, .maximum_temperature_k = context.runscript.maximum_surface_temperature_k, .solver_options = context.nonlinear_solver_options } };
    try runKernelAcrossSerialTiles(context, &surface_temperature_context, ecosys.surface_temperature_solver.applyTile);
    @memset(context.delayed_surface_combustion_heat_mj, 0);
    @memset(context.surface_combustion_heat_mj_per_m2, 0);
    if (context.canopy_airflow.*) |*airflow| if (context.canopy_surface_exchange.*) |*exchange| if (context.canopy_surface_input_workspace.*) |*surface_workspace| if (context.detailed_canopy.*) |*canopy| if (context.canopy_precipitation_retention.*) |*retention| if (context.canopy_exposure.*) |*exposure| {
        try ecosys.plant_development.refreshCanopyHeight(canopy, context.development_canopy_height_m);
        for (0..context.grid.cell_count) |cell| {
            context.canopy_atmospheric_vapor_diffusivity_m2_per_h[cell] = context.runscript.soil_process_parameters.reference_water_vapor_diffusivity_m2_per_h * std.math.pow(f64, context.atmosphere.air_temperature_k[cell] / context.runscript.soil_process_parameters.vapor_diffusivity_reference_temperature_k, context.runscript.soil_process_parameters.vapor_diffusivity_temperature_exponent);
            context.canopy_surface_roughness_height_m[cell] = if (context.snow_depth_m[cell] > 0) context.runscript.surface_aerodynamic_parameters.snow_roughness_height_m else context.runscript.surface_aerodynamic_parameters.soil_roughness_height_m;
        }
        for (0..context.config.plant_populations * context.grid.cell_count) |plant| context.canopy_available_intercepted_water_m3[plant] = @max(0, retention.living_surface_water_m3[plant] + retention.living_retention_m3_per_h[plant]);
        var airflow_context: ecosys.canopy_airflow.ApplyContext = .{
            .state = airflow,
            .cell_area_m2 = context.canopy_cell_area_m2,
            .total_canopy_area_m2 = context.surface_total_canopy_area_m2,
            .biome_canopy_height_m = context.surface_canopy_height_m,
            .surface_roughness_height_m = context.canopy_surface_roughness_height_m,
            .species_canopy_height_m = context.development_canopy_height_m,
            .standing_dead_height_m = canopy.plant_standing_dead_height_m,
            .atmospheric_vapor_diffusivity_m2_per_h = context.canopy_atmospheric_vapor_diffusivity_m2_per_h,
            .atmospheric_temperature_k = context.atmosphere.air_temperature_k,
            .ground_air_temperature_k = context.ground_air.temperature_k,
            .species_canopy_air_temperature_k = canopy.plant_canopy_aerodynamic_temperature_k,
            .standing_dead_air_temperature_k = canopy.plant_standing_dead_aerodynamic_temperature_k,
            .bulk_richardson_coefficient_k = context.surface_aerodynamics.bulk_richardson_coefficient_k,
            .parameters = .{
                .minimum_richardson_number = context.runscript.canopy_surface_exchange_parameters.minimum_richardson_number,
                .maximum_richardson_number = context.runscript.canopy_surface_exchange_parameters.maximum_richardson_number,
                .richardson_resistance_multiplier = context.runscript.canopy_surface_exchange_parameters.richardson_resistance_multiplier,
                .minimum_canopy_resistance_h_per_m = context.runscript.canopy_surface_exchange_parameters.minimum_boundary_resistance_h_per_m,
                .maximum_canopy_resistance_h_per_m = context.runscript.canopy_surface_exchange_parameters.maximum_boundary_resistance_h_per_m,
                .canopy_drag_length_m = context.runscript.surface_gas_resistance_parameters.canopy_drag_length_m,
                .volumetric_air_heat_capacity_mj_per_m3_k = context.runscript.ground_air_parameters.volumetric_air_heat_capacity_mj_per_m3_k,
            },
        };
        try runKernelAcrossSerialTiles(context, &airflow_context, ecosys.canopy_airflow.applyTile);
        try surface_workspace.refresh(
            canopy.plant_canopy_aerodynamic_temperature_k,
            canopy.plant_canopy_aerodynamic_vapor_pressure_kpa,
            canopy.plant_minimum_water_vapor_resistance_h_per_m,
            canopy.plant_cuticular_water_vapor_resistance_h_per_m,
            context.plant_reproduction_controls.stomatal_turgor_shape,
            canopy.plant_canopy_turgor_potential_mpa,
            context.runscript.canopy_sensible_surface_resistance_h_per_m,
            context.runscript.canopy_latent_surface_resistance_h_per_m,
            context.runscript.canopy_surface_exchange_parameters,
        );
        var exchange_context: ecosys.canopy_surface_exchange.ApplyContext = .{
            .state = exchange,
            .inputs = .{
                .atmospheric_temperature_k_by_cell = context.atmosphere.air_temperature_k,
                .canopy_air_temperature_k = canopy.plant_canopy_aerodynamic_temperature_k,
                .canopy_air_vapor_fraction = surface_workspace.canopy_air_vapor_fraction,
                .bulk_richardson_coefficient_k_by_cell = context.surface_aerodynamics.bulk_richardson_coefficient_k,
                .biome_isothermal_boundary_resistance_h_per_m_by_cell = context.surface_aerodynamics.isothermal_aerodynamic_resistance_h_per_m,
                .latent_boundary_numerator_m2_per_h_by_cell = airflow.latent_boundary_numerator_m2_per_h,
                .sensible_boundary_numerator_mj_per_m_h_k_by_cell = airflow.sensible_boundary_numerator_mj_per_m_h_k,
                .canopy_surface_temperature_k = context.plants.canopy_temperature_k,
                .aerodynamic_resistance_below_biome_h_per_m_by_cell = airflow.resistance_below_biome_h_per_m,
                .aerodynamic_resistance_below_species_h_per_m = airflow.resistance_below_species_h_per_m,
                .species_canopy_radiation_fraction = retention.living_radiation_fraction,
                .sensible_surface_resistance_h_per_m = surface_workspace.sensible_surface_resistance_h_per_m,
                .latent_surface_resistance_h_per_m = surface_workspace.latent_surface_resistance_h_per_m,
                .stomatal_resistance_h_per_m = surface_workspace.stomatal_resistance_h_per_m,
                .canopy_total_water_potential_mpa = context.plants.canopy_water_potential_mpa,
                .intercepted_water_volume_m3 = context.canopy_available_intercepted_water_m3,
            },
            .parameters = context.runscript.canopy_surface_exchange_parameters,
        };
        try runKernelAcrossSerialTiles(context, &exchange_context, ecosys.canopy_surface_exchange.applyTile);
        @memcpy(canopy.plant_transpiration_m3_per_h, exchange.transpiration_m3_per_h);
        if (context.standing_dead_surface_exchange.*) |*dead_exchange| {
            for (0..context.grid.cell_count) |cell| for (0..context.config.plant_populations) |species| {
                const plant = cell * context.config.plant_populations + species;
                context.standing_dead_evaporation_m3_per_h[plant] = 0;
                dead_exchange.intercepted_water_change_m3_per_h[plant] = 0;
                dead_exchange.net_radiation_mj_per_h[plant] = 0;
                dead_exchange.sensible_heat_flux_mj_per_h[plant] = 0;
                dead_exchange.latent_heat_flux_mj_per_h[plant] = 0;
                dead_exchange.vapor_sensible_heat_flux_mj_per_h[plant] = 0;
                dead_exchange.storage_heat_flux_mj_per_h[plant] = 0;
                const standing_dead_present =
                    retention.standing_dead_radiation_fraction[plant] > 1.0e-12 and
                    canopy.plant_standing_dead_height_m[plant] > 0;
                if (!standing_dead_present) {
                    const ambient_vapor_fraction =
                        try ecosys.ground_air_exchange.vaporVolumeFraction(
                            context.atmosphere.vapor_pressure_kpa[cell],
                            context.atmosphere.air_temperature_k[cell],
                            context.runscript.ground_air_parameters,
                        );
                    const fallback = try ecosys.absent_standing_dead_canopy_state.apply(
                        false,
                        std.mem.zeroes(ecosys.absent_standing_dead_canopy_state.State),
                        .{
                            .intercepted_water_m3 = 0,
                            .intercepted_water_rate_m3_per_h = 0,
                            .timestep_h = 0,
                            .ambient_air_temperature_k = context.atmosphere.air_temperature_k[cell],
                            .ambient_vapor_volume_fraction = ambient_vapor_fraction,
                            .canopy_height_m = canopy.plant_standing_dead_height_m[plant],
                            .snow_surface_depth_m = context.snow_depth_m[cell],
                            .depth_tolerance_m = 1.0e-12,
                            .topsoil_temperature_k = context.grid.soil_temperature_k[cell * context.grid.soil_layer_capacity],
                        },
                    );
                    canopy.plant_standing_dead_aerodynamic_temperature_k[plant] =
                        fallback.canopy_air_temperature_k;
                    canopy.plant_standing_dead_aerodynamic_vapor_pressure_kpa[plant] =
                        context.atmosphere.vapor_pressure_kpa[cell];
                    canopy.plant_standing_dead_surface_temperature_k[plant] =
                        fallback.canopy_surface_temperature_k;
                    continue;
                }
                const air_temperature_k = canopy.plant_standing_dead_aerodynamic_temperature_k[plant];
                const exchange_inputs: ecosys.standing_dead_surface_exchange.Inputs = .{
                    .atmospheric_temperature_k = context.atmosphere.air_temperature_k[cell],
                    .standing_dead_air_temperature_k = air_temperature_k,
                    .standing_dead_surface_temperature_k = canopy.plant_standing_dead_surface_temperature_k[plant],
                    .standing_dead_air_vapor_fraction = canopy.plant_standing_dead_aerodynamic_vapor_pressure_kpa[plant] * context.runscript.canopy_surface_exchange_parameters.saturation_vapor_prefactor_k / air_temperature_k,
                    .bulk_richardson_coefficient_k = context.surface_aerodynamics.bulk_richardson_coefficient_k[cell],
                    .biome_isothermal_boundary_resistance_h_per_m = context.surface_aerodynamics.isothermal_aerodynamic_resistance_h_per_m[cell],
                    .aerodynamic_resistance_below_biome_h_per_m = airflow.resistance_below_biome_h_per_m[cell],
                    .aerodynamic_resistance_below_standing_dead_h_per_m = airflow.resistance_below_standing_dead_h_per_m[plant],
                    .standing_dead_radiation_fraction = retention.standing_dead_radiation_fraction[plant],
                    .latent_boundary_numerator_m2_per_h = airflow.latent_boundary_numerator_m2_per_h[cell],
                    .sensible_boundary_numerator_mj_per_m_h_k = airflow.sensible_boundary_numerator_mj_per_m_h_k[cell],
                    .sensible_surface_resistance_h_per_m = context.runscript.canopy_sensible_surface_resistance_h_per_m,
                    .latent_surface_resistance_h_per_m = context.runscript.canopy_latent_surface_resistance_h_per_m,
                    .intercepted_water_volume_m3 = @max(0, retention.standing_dead_surface_water_m3[plant] + retention.standing_dead_retention_m3_per_h[plant]),
                };
                const dead_parameters: ecosys.standing_dead_surface_exchange.Parameters = .{
                    .minimum_richardson_number = context.runscript.canopy_surface_exchange_parameters.minimum_richardson_number,
                    .maximum_richardson_number = context.runscript.canopy_surface_exchange_parameters.maximum_richardson_number,
                    .richardson_resistance_multiplier = context.runscript.canopy_surface_exchange_parameters.richardson_resistance_multiplier,
                    .minimum_boundary_resistance_h_per_m = context.runscript.canopy_surface_exchange_parameters.minimum_boundary_resistance_h_per_m,
                    .maximum_boundary_resistance_h_per_m = context.runscript.canopy_surface_exchange_parameters.maximum_boundary_resistance_h_per_m,
                    .saturation_vapor_prefactor_k = context.runscript.canopy_surface_exchange_parameters.saturation_vapor_prefactor_k,
                    .saturation_relative_humidity = context.runscript.canopy_surface_exchange_parameters.saturation_relative_humidity,
                    .saturation_temperature_coefficient_k = context.runscript.canopy_surface_exchange_parameters.saturation_temperature_k,
                    .saturation_reference_inverse_temperature_per_k = context.runscript.canopy_surface_exchange_parameters.saturation_reference_inverse_temperature_per_k,
                    .latent_heat_of_vaporization_mj_per_m3 = context.runscript.canopy_surface_exchange_parameters.latent_heat_of_vaporization_mj_per_m3,
                    .liquid_water_heat_capacity_mj_per_m3_k = context.runscript.canopy_surface_exchange_parameters.liquid_water_heat_capacity_mj_per_m3_k,
                };
                const active_dry_volume_m3 = @min(
                    context.runscript.standing_dead_sapwood_thickness_m * retention.standing_dead_surface_area_m2[plant],
                    canopy.plant_standing_dead_carbon_g[plant] * context.runscript.stalk_volume_m3_per_g_c,
                );
                const dry_heat_capacity_mj_per_k = context.runscript.standing_dead_dry_volume_heat_capacity_mj_per_m3_k * active_dry_volume_m3;
                const activation_threshold_mj_per_k = context.runscript.standing_dead_activation_heat_capacity_mj_per_m2_k * context.canopy_cell_area_m2[cell];
                if (dry_heat_capacity_mj_per_k <= activation_threshold_mj_per_k) continue;
                const wet_heat_capacity_mj_per_k = dry_heat_capacity_mj_per_k +
                    dead_parameters.liquid_water_heat_capacity_mj_per_m3_k * retention.standing_dead_surface_water_m3[plant];
                const solved_surface = try ecosys.standing_dead_surface_exchange.solveSurfaceTemperature(.{
                    .exchange_inputs = exchange_inputs,
                    .absorbed_shortwave_mj_per_h = retention.standing_dead_absorbed_shortwave_mj_per_m2[plant] * context.canopy_cell_area_m2[cell],
                    .downward_longwave_mj_per_h = context.atmosphere.longwave_radiation_mj_per_m2[cell] * retention.standing_dead_radiation_fraction[plant] * context.canopy_cell_area_m2[cell],
                    .lateral_longwave_mj_per_h = 0,
                    .ground_surface_temperature_k = context.grid.surface_temperature_k[cell],
                    .emission_coefficient_mj_per_h_k4 = context.runscript.standing_dead_emissivity * 2.04e-10 * retention.standing_dead_radiation_fraction[plant] * context.canopy_cell_area_m2[cell],
                    .dry_and_existing_water_heat_capacity_mj_per_k = wet_heat_capacity_mj_per_k,
                    .retained_precipitation_water_m3_per_h = retention.standing_dead_retention_m3_per_h[plant],
                    .retained_precipitation_heat_mj_per_h = 0,
                    .minimum_effective_heat_capacity_mj_per_k = context.runscript.standing_dead_effective_heat_capacity_floor_mj_per_m2_k * context.canopy_cell_area_m2[cell],
                }, dead_parameters, .{
                    .minimum_temperature_k = context.runscript.minimum_surface_temperature_k,
                    .maximum_temperature_k = context.runscript.maximum_surface_temperature_k,
                    .solver_options = .{
                        .absolute_tolerance = context.config.absolute_tolerance,
                        .relative_tolerance = context.config.relative_tolerance,
                        .max_iterations = try context.iteration_limits.standingDeadEnergyMaxIterations(),
                        .picard_relaxation = context.config.picard_relaxation,
                    },
                });
                canopy.plant_standing_dead_surface_temperature_k[plant] = solved_surface.temperature_k;
                const result = solved_surface.exchange;
                dead_exchange.intercepted_water_change_m3_per_h[plant] = result.intercepted_water_change_m3_per_h;
                dead_exchange.net_radiation_mj_per_h[plant] = solved_surface.net_radiation_mj_per_h;
                dead_exchange.sensible_heat_flux_mj_per_h[plant] = result.sensible_heat_flux_mj_per_h;
                dead_exchange.latent_heat_flux_mj_per_h[plant] = result.latent_heat_flux_mj_per_h;
                dead_exchange.vapor_sensible_heat_flux_mj_per_h[plant] = result.vapor_sensible_heat_flux_mj_per_h;
                dead_exchange.storage_heat_flux_mj_per_h[plant] = solved_surface.storage_heat_flux_mj_per_h;
                context.standing_dead_evaporation_m3_per_h[plant] = result.intercepted_water_change_m3_per_h;
                if (context.standing_dead_air_exchange.*) |*dead_air| {
                    var total_canopy_exposure: f64 = 0;
                    for (0..context.config.plant_populations) |population| {
                        const population_index = cell * context.config.plant_populations + population;
                        total_canopy_exposure += retention.living_radiation_fraction[population_index] + retention.standing_dead_radiation_fraction[population_index];
                    }
                    const dead_share = if (total_canopy_exposure > 1.0e-12) retention.standing_dead_radiation_fraction[plant] / total_canopy_exposure else 0;
                    if (dead_share > 1.0e-12) {
                        const air_column_height_m = @max(5.0, context.surface_aerodynamics.wind_reference_height_m[cell]);
                        const cell_air_volume_m3 = air_column_height_m * context.canopy_cell_area_m2[cell];
                        const atmospheric_sensible_conductance = airflow.sensible_boundary_numerator_mj_per_m_h_k[cell] * retention.standing_dead_radiation_fraction[plant] / result.total_aerodynamic_resistance_h_per_m * dead_share;
                        const atmospheric_vapor_conductance = @min(
                            airflow.latent_boundary_numerator_m2_per_h[cell] * retention.standing_dead_radiation_fraction[plant] / result.total_aerodynamic_resistance_h_per_m,
                            cell_air_volume_m3,
                        ) * dead_share;
                        const ground_sensible_conductance = airflow.sensible_boundary_numerator_mj_per_m_h_k[cell] / airflow.resistance_below_standing_dead_h_per_m[plant] * total_canopy_exposure * dead_share;
                        const ground_vapor_conductance = airflow.latent_boundary_numerator_m2_per_h[cell] / airflow.resistance_below_standing_dead_h_per_m[plant] * total_canopy_exposure * dead_share;
                        const dead_air_result = try ecosys.canopy_air_exchange.solveInto(dead_air, cell, species, .{
                            .initial_temperature_k = air_temperature_k,
                            .initial_vapor_fraction = canopy.plant_standing_dead_aerodynamic_vapor_pressure_kpa[plant] * context.runscript.canopy_surface_exchange_parameters.saturation_vapor_prefactor_k / air_temperature_k,
                            .atmospheric_temperature_k = context.atmosphere.air_temperature_k[cell],
                            .atmospheric_vapor_fraction = try ecosys.ground_air_exchange.vaporVolumeFraction(context.atmosphere.vapor_pressure_kpa[cell], context.atmosphere.air_temperature_k[cell], context.runscript.ground_air_parameters),
                            .ground_air_temperature_k = context.ground_air.temperature_k[cell],
                            .ground_air_vapor_fraction = context.ground_air.vapor_volume_fraction[cell],
                            .heat_capacity_mj_per_k = cell_air_volume_m3 * context.runscript.ground_air_parameters.volumetric_air_heat_capacity_mj_per_m3_k * dead_share,
                            .air_volume_m3 = cell_air_volume_m3 * dead_share,
                            .atmospheric_sensible_conductance_mj_per_h_k = atmospheric_sensible_conductance,
                            .atmospheric_vapor_conductance_m3_per_h = atmospheric_vapor_conductance,
                            .ground_sensible_conductance_mj_per_h_k = ground_sensible_conductance,
                            .ground_vapor_conductance_m3_per_h = ground_vapor_conductance,
                            .canopy_surface_sensible_heat_flux_mj_per_h = result.sensible_heat_flux_mj_per_h,
                            .canopy_surface_vapor_flux_m3_per_h = result.intercepted_water_change_m3_per_h,
                            .lateral_sensible_heat_flux_mj_per_h = context.delayed_standing_dead_combustion_heat_mj[plant],
                            .lateral_vapor_flux_m3_per_h = 0,
                        }, .{
                            .saturation_vapor_prefactor_k = context.runscript.canopy_surface_exchange_parameters.saturation_vapor_prefactor_k,
                            .saturation_relative_humidity = context.runscript.canopy_surface_exchange_parameters.saturation_relative_humidity,
                            .saturation_temperature_coefficient_k = context.runscript.canopy_surface_exchange_parameters.saturation_temperature_k,
                            .saturation_reference_inverse_temperature_per_k = context.runscript.canopy_surface_exchange_parameters.saturation_reference_inverse_temperature_per_k,
                        }, .{
                            .max_iterations = context.iteration_limits.canopy_energy_water_max_iterations,
                            .picard_relaxation = context.config.picard_relaxation,
                        });
                        canopy.plant_standing_dead_aerodynamic_temperature_k[plant] = dead_air_result.temperature_k;
                        canopy.plant_standing_dead_aerodynamic_vapor_pressure_kpa[plant] = dead_air_result.vapor_fraction * dead_air_result.temperature_k / context.runscript.canopy_surface_exchange_parameters.saturation_vapor_prefactor_k;
                    }
                }
            };
        }
        if (context.canopy_air_exchange.*) |*canopy_air| {
            for (0..context.grid.cell_count) |cell| {
                const air_column_height_m = @max(5.0, context.surface_aerodynamics.wind_reference_height_m[cell]);
                const cell_air_volume_m3 = air_column_height_m * context.canopy_cell_area_m2[cell];
                const cell_air_heat_capacity_mj_per_k = cell_air_volume_m3 * context.runscript.ground_air_parameters.volumetric_air_heat_capacity_mj_per_m3_k;
                var total_canopy_radiation_fraction: f64 = 0;
                for (0..context.config.plant_populations) |population| {
                    const population_index = cell * context.config.plant_populations + population;
                    total_canopy_radiation_fraction += retention.living_radiation_fraction[population_index] + retention.standing_dead_radiation_fraction[population_index];
                }
                for (0..context.config.plant_populations) |species| {
                    const plant = cell * context.config.plant_populations + species;
                    const exposure_fraction = retention.living_radiation_fraction[plant];
                    const canopy_share = if (total_canopy_radiation_fraction > 1.0e-12) exposure_fraction / total_canopy_radiation_fraction else 0;
                    if (exposure_fraction <= 1.0e-12 or canopy_share <= 1.0e-12) continue;
                    const total_resistance_h_per_m = exchange.total_aerodynamic_resistance_h_per_m[plant];
                    const below_species_resistance_h_per_m = airflow.resistance_below_species_h_per_m[plant];
                    const atmospheric_sensible_conductance = airflow.sensible_boundary_numerator_mj_per_m_h_k[cell] * exposure_fraction / total_resistance_h_per_m * canopy_share;
                    const atmospheric_vapor_conductance = @min(
                        airflow.latent_boundary_numerator_m2_per_h[cell] * exposure_fraction / total_resistance_h_per_m,
                        cell_air_volume_m3,
                    ) * canopy_share;
                    const ground_sensible_conductance = airflow.sensible_boundary_numerator_mj_per_m_h_k[cell] / below_species_resistance_h_per_m * exposure.canopy_exposure_fraction[cell] * canopy_share;
                    const ground_vapor_conductance = airflow.latent_boundary_numerator_m2_per_h[cell] / below_species_resistance_h_per_m * exposure.canopy_exposure_fraction[cell] * canopy_share;
                    const result = try ecosys.canopy_air_exchange.solveInto(canopy_air, cell, species, .{
                        .initial_temperature_k = canopy.plant_canopy_aerodynamic_temperature_k[plant],
                        .initial_vapor_fraction = surface_workspace.canopy_air_vapor_fraction[plant],
                        .atmospheric_temperature_k = context.atmosphere.air_temperature_k[cell],
                        .atmospheric_vapor_fraction = try ecosys.ground_air_exchange.vaporVolumeFraction(context.atmosphere.vapor_pressure_kpa[cell], context.atmosphere.air_temperature_k[cell], context.runscript.ground_air_parameters),
                        .ground_air_temperature_k = context.ground_air.temperature_k[cell],
                        .ground_air_vapor_fraction = context.ground_air.vapor_volume_fraction[cell],
                        .heat_capacity_mj_per_k = cell_air_heat_capacity_mj_per_k * canopy_share,
                        .air_volume_m3 = cell_air_volume_m3 * canopy_share,
                        .atmospheric_sensible_conductance_mj_per_h_k = atmospheric_sensible_conductance,
                        .atmospheric_vapor_conductance_m3_per_h = atmospheric_vapor_conductance,
                        .ground_sensible_conductance_mj_per_h_k = ground_sensible_conductance,
                        .ground_vapor_conductance_m3_per_h = ground_vapor_conductance,
                        .canopy_surface_sensible_heat_flux_mj_per_h = exchange.sensible_heat_flux_mj_per_h[plant],
                        .canopy_surface_vapor_flux_m3_per_h = exchange.intercepted_water_change_m3_per_h[plant] + exchange.transpiration_m3_per_h[plant],
                        .lateral_sensible_heat_flux_mj_per_h = context.delayed_live_canopy_combustion_heat_mj[plant],
                        .lateral_vapor_flux_m3_per_h = 0,
                    }, .{
                        .saturation_vapor_prefactor_k = context.runscript.canopy_surface_exchange_parameters.saturation_vapor_prefactor_k,
                        .saturation_relative_humidity = context.runscript.canopy_surface_exchange_parameters.saturation_relative_humidity,
                        .saturation_temperature_coefficient_k = context.runscript.canopy_surface_exchange_parameters.saturation_temperature_k,
                        .saturation_reference_inverse_temperature_per_k = context.runscript.canopy_surface_exchange_parameters.saturation_reference_inverse_temperature_per_k,
                    }, .{
                        .max_iterations = context.iteration_limits.canopy_energy_water_max_iterations,
                        .picard_relaxation = context.config.picard_relaxation,
                    });
                    canopy.plant_canopy_aerodynamic_temperature_k[plant] = result.temperature_k;
                    canopy.plant_canopy_aerodynamic_vapor_pressure_kpa[plant] = result.vapor_fraction * result.temperature_k / context.runscript.canopy_surface_exchange_parameters.saturation_vapor_prefactor_k;
                }
            }
        }
    };
    for (0..context.grid.cell_count) |cell| {
        context.atmospheric_vapor_fraction[cell] = try ecosys.ground_air_exchange.vaporVolumeFraction(context.atmosphere.vapor_pressure_kpa[cell], context.atmosphere.air_temperature_k[cell], context.runscript.ground_air_parameters);
        context.ground_air_canopy_resistance_h_per_m[cell] = if (context.canopy_airflow.*) |*airflow| airflow.neutral_resistance_below_biome_h_per_m[cell] else 0;
    }
    @memset(context.delayed_live_canopy_combustion_heat_mj, 0);
    @memset(context.delayed_standing_dead_combustion_heat_mj, 0);
    @memset(context.ground_air_sensible_source_mj_per_h, 0);
    @memset(context.ground_air_vapor_source_m3_per_h, 0);
    if (context.canopy_air_exchange.*) |*canopy_air| if (context.canopy_airflow.*) |*airflow| if (context.canopy_precipitation_retention.*) |*retention| {
        for (0..context.grid.cell_count) |cell| {
            var total_canopy_exposure: f64 = 0;
            for (0..context.config.plant_populations) |population| {
                const population_index = cell * context.config.plant_populations + population;
                total_canopy_exposure += retention.living_radiation_fraction[population_index] + retention.standing_dead_radiation_fraction[population_index];
            }
            for (0..context.config.plant_populations) |species| {
                const plant = cell * context.config.plant_populations + species;
                const canopy_share = if (total_canopy_exposure > 1.0e-12) retention.living_radiation_fraction[plant] / total_canopy_exposure else 0;
                if (canopy_share <= 1.0e-12) continue;
                const below_species_resistance_h_per_m = airflow.resistance_below_species_h_per_m[plant];
                if (below_species_resistance_h_per_m <= 0) return error.InvalidCanopyGroundAirResistance;
                const sensible_conductance = airflow.sensible_boundary_numerator_mj_per_m_h_k[cell] / below_species_resistance_h_per_m * total_canopy_exposure * canopy_share;
                const vapor_conductance = airflow.latent_boundary_numerator_m2_per_h[cell] / below_species_resistance_h_per_m * total_canopy_exposure * canopy_share;
                context.ground_air_sensible_source_mj_per_h[cell] += sensible_conductance * (canopy_air.temperature_k[plant] - context.ground_air.temperature_k[cell]);
                context.ground_air_vapor_source_m3_per_h[cell] += vapor_conductance * (canopy_air.vapor_fraction[plant] - context.ground_air.vapor_volume_fraction[cell]);
            }
        }
    };
    if (context.standing_dead_air_exchange.*) |*dead_air| if (context.canopy_airflow.*) |*airflow| if (context.canopy_precipitation_retention.*) |*retention| {
        for (0..context.grid.cell_count) |cell| {
            var total_canopy_exposure: f64 = 0;
            for (0..context.config.plant_populations) |population| {
                const population_index = cell * context.config.plant_populations + population;
                total_canopy_exposure += retention.living_radiation_fraction[population_index] + retention.standing_dead_radiation_fraction[population_index];
            }
            for (0..context.config.plant_populations) |species| {
                const plant = cell * context.config.plant_populations + species;
                const dead_share = if (total_canopy_exposure > 1.0e-12) retention.standing_dead_radiation_fraction[plant] / total_canopy_exposure else 0;
                if (dead_share <= 1.0e-12) continue;
                const resistance_h_per_m = airflow.resistance_below_standing_dead_h_per_m[plant];
                if (resistance_h_per_m <= 0) return error.InvalidStandingDeadGroundAirResistance;
                const sensible_conductance = airflow.sensible_boundary_numerator_mj_per_m_h_k[cell] / resistance_h_per_m * total_canopy_exposure * dead_share;
                const vapor_conductance = airflow.latent_boundary_numerator_m2_per_h[cell] / resistance_h_per_m * total_canopy_exposure * dead_share;
                context.ground_air_sensible_source_mj_per_h[cell] += sensible_conductance * (dead_air.temperature_k[plant] - context.ground_air.temperature_k[cell]);
                context.ground_air_vapor_source_m3_per_h[cell] += vapor_conductance * (dead_air.vapor_fraction[plant] - context.ground_air.vapor_volume_fraction[cell]);
            }
        }
    };
    for (0..context.grid.cell_count) |cell| {
        context.ground_air_surface_sensible_conductance_mj_per_h_k[cell] = context.runscript.surface_sensible_heat_conductance_mj_per_m2_h_k * context.canopy_cell_area_m2[cell];
        context.ground_air_surface_vapor_conductance_m3_per_h[cell] = context.runscript.surface_latent_heat_conductance_mj_per_m2_h_kpa * context.canopy_cell_area_m2[cell] / context.runscript.ground_air_parameters.liquid_water_latent_heat_mj_per_m3 * context.ground_air.temperature_k[cell] / context.runscript.ground_air_parameters.saturation_vapor_prefactor_k;
        const surface_temperature_k = context.grid.surface_temperature_k[cell];
        context.ground_air_surface_vapor_fraction[cell] = context.runscript.ground_air_parameters.saturation_vapor_prefactor_k / surface_temperature_k * context.runscript.ground_air_parameters.saturation_relative_humidity * context.runscript.surface_vapor_activity_fraction * @exp(context.runscript.ground_air_parameters.saturation_temperature_k * (context.runscript.ground_air_parameters.saturation_reference_inverse_temperature_per_k - 1 / surface_temperature_k));
    }
    try ecosys.ground_air_exchange.solve(context.ground_air, .{
        .atmospheric_temperature_k = context.atmosphere.air_temperature_k,
        .atmospheric_vapor_volume_fraction = context.atmospheric_vapor_fraction,
        .cell_area_m2 = context.canopy_cell_area_m2,
        .bulk_richardson_coefficient_k = context.surface_aerodynamics.bulk_richardson_coefficient_k,
        .neutral_atmospheric_resistance_h_per_m = context.surface_aerodynamics.isothermal_aerodynamic_resistance_h_per_m,
        .canopy_resistance_h_per_m = context.ground_air_canopy_resistance_h_per_m,
        .non_atmospheric_sensible_heat_mj_per_h = context.ground_air_sensible_source_mj_per_h,
        .non_atmospheric_vapor_flux_m3_per_h = context.ground_air_vapor_source_m3_per_h,
        .non_atmospheric_sensible_conductance_mj_per_h_k = context.ground_air_surface_sensible_conductance_mj_per_h_k,
        .non_atmospheric_sensible_source_temperature_k = context.grid.surface_temperature_k,
        .non_atmospheric_vapor_conductance_m3_per_h = context.ground_air_surface_vapor_conductance_m3_per_h,
        .non_atmospheric_vapor_source_fraction = context.ground_air_surface_vapor_fraction,
    }, context.runscript.ground_air_parameters, .{
        .absolute_tolerance = context.config.absolute_tolerance,
        .relative_tolerance = context.config.relative_tolerance,
        .max_iterations = context.iteration_limits.water_heat_solute_max_iterations,
        .picard_relaxation = context.config.picard_relaxation,
    });
    try ecosys.surface_precipitation.commitSoilIngress(context.surface_precipitation, context.grid, context.transport_hydrology, 1);
    try runKernelAcrossSerialTiles(context, context.soil_thermal_context, ecosys.soil_thermal.updateTile);
    try context.soil_hourly_workspace.refresh(context.grid, context.soil_solver_properties, context.soil_thermal, context.terrain_hydrology, context.runscript.soil_process_parameters);
    var fertilizer_band_workspace = try ecosys.fertilizer_band_production.Workspace.init(context.allocator, context.grid.soil_layer_capacity);
    defer fertilizer_band_workspace.deinit();
    try ecosys.fertilizer_band_production.prepareHour(
        context.fertilizer_band,
        fertilizer_band_hour,
        context.soil_solver_properties.layer_thickness_m,
        context.grid.matrix_liquid_water_m3,
        context.grid.matrix_pore_capacity_m3,
        context.grid.soil_temperature_k,
        context.grid.active_soil_layer_count,
        context.runscript.soil_geometry_parameters.minimum_layer_thickness_m,
        context.runscript.root_nutrient_parameters,
        &fertilizer_band_workspace,
    );
    try context.soil_hourly_workspace.fillMacroporeFaceConductance(context.soil_transport_faces, context.soil_face_geometry, context.soil_face_geometry.macropore_hydraulic_conductance_m_per_h_mpa);
    try context.soil_hourly_workspace.bindSurfaceHeatFlux(context.grid, context.surface_temperature);
    try ecosys.surface_precipitation.bindSoilHeatIngress(context.surface_precipitation, context.grid, context.soil_hourly_workspace.cell_heat_source_mj, 1);
    _ = try ecosys.subsurface_irrigation_heat.addToLayerHeatSources(
        context.soil_hourly_workspace.cell_heat_source_mj,
        context.subsurface_irrigation_water_m3,
        context.atmosphere.air_temperature_k,
        context.grid.soil_layer_capacity,
        context.runscript.canopy_surface_exchange_parameters.liquid_water_heat_capacity_mj_per_m3_k,
    );
    for (context.delayed_subsurface_combustion_heat_mj, context.soil_hourly_workspace.cell_heat_source_mj) |*delayed_heat, *heat_source| {
        heat_source.* += delayed_heat.*;
        delayed_heat.* = 0;
        if (!std.math.isFinite(heat_source.*)) return error.NonFiniteSubsurfaceCombustionHeatSource;
    }
    const water_heat_solute_max_iterations = try ecosys.iteration_control.waterHeatSoluteCeilingForCurrentState(context.iteration_limits.water_heat_solute_max_iterations, context.soil_hourly_workspace.heat_capacity_mj_per_k, context.soil_hourly_workspace.horizontal_face_area_m2, context.soil_hourly_workspace.is_top_soil_layer);
    var accepted_soil_water_heat = try ecosys.soil_water_heat_step.advanceMappedDeferred(context.allocator, context.grid, context.transport_hydrology, context.soil_transport_faces, context.soil_face_geometry, context.soil_solver_properties, context.soil_hourly_workspace, context.soil_thermal, context.soil_heat_solver_workspace, context.runscript.soil_phase_heat_parameters, .{ .max_iterations = water_heat_solute_max_iterations, .picard_relaxation = context.config.picard_relaxation, .vapor_pore_tortuosity = context.runscript.soil_process_parameters.vapor_pore_tortuosity, .osmotic_reflection_coefficient = context.runscript.soil_process_parameters.osmotic_reflection_coefficient, .absolute_tolerance = context.config.absolute_tolerance, .relative_tolerance = context.config.relative_tolerance, .boundary_topology = context.soil_boundary_topology, .geothermal_enabled_by_cell = context.geothermal_enabled_by_cell, .mean_annual_temperature_k_by_cell = context.mean_annual_temperature_k_by_cell, .geothermal_minimum_source_depth_m = context.runscript.geothermal_controls.minimum_source_depth_m, .geothermal_source_depth_below_profile_m = context.runscript.geothermal_controls.source_depth_below_profile_m, .geothermal_conductivity_m_mj_per_h_k = context.runscript.geothermal_controls.conductivity_m_mj_per_h_k, .geothermal_flux_mj_per_m2_h = context.runscript.geothermal_controls.geothermal_flux_mj_per_m2_h, .water_table_air_fraction_threshold = context.runscript.water_table_air_fraction_threshold, .active_layer_ice_fraction_threshold = context.runscript.active_layer_ice_fraction_threshold, .dense_newton_max_components = context.config.tile_cells, .matrix_external_water_source_m3_per_step = context.subsurface_irrigation_water_m3 });
    defer accepted_soil_water_heat.deinit();
    try commitHourlyWaterHeatStateGeneration(
        context,
        &accepted_soil_water_heat,
    );
    try context.landscape_boundary_ledger.accumulateAccepted(.{
        .heat_input_mj = accepted_soil_water_heat.solver.heat.boundary_heat_input_mj,
        .heat_output_mj = accepted_soil_water_heat.solver.heat.boundary_heat_output_mj,
    });
    const subsurface_irrigation_chemistry_parameters: ecosys.subsurface_irrigation_chemistry.Parameters = .{
        .molar_mass_g_per_mol = if (context.runscript.chemistry_primary_initialization) |parameters|
            .{
                .nitrogen = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                .phosphorus = context.runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                .aluminum = parameters.molar_mass_g_per_mol.aluminum,
                .iron = parameters.molar_mass_g_per_mol.iron,
                .calcium = parameters.molar_mass_g_per_mol.calcium,
                .magnesium = parameters.molar_mass_g_per_mol.magnesium,
                .sodium = parameters.molar_mass_g_per_mol.sodium,
                .potassium = parameters.molar_mass_g_per_mol.potassium,
                .sulfur = parameters.molar_mass_g_per_mol.sulfur,
                .chloride = parameters.molar_mass_g_per_mol.chloride,
            }
        else
            .{
                .nitrogen = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                .phosphorus = context.runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                .aluminum = 27,
                .iron = 55.8,
                .calcium = 40,
                .magnesium = 24.3,
                .sodium = 23,
                .potassium = 39.1,
                .sulfur = 32,
                .chloride = 35.5,
            },
        .equilibrium = if (context.chemistry_reaction_parameters.*) |parameters|
            .{
                .aqueous = parameters.aqueous_constants,
                .phosphate = parameters.phosphate_constants,
            }
        else
            null,
        .ammonium_band_fraction = context.runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
        .nitrate_band_fraction = context.runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
        .phosphate_band_fraction = context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
    };
    try ecosys.subsurface_irrigation_chemistry.addTransportedIons(
        context.irrigation_loads,
        context.micropore_solute_state,
        subsurface_irrigation_chemistry_parameters,
    );
    const phosphate_band_fraction = context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction;
    try context.soil_solute_face_parameters.refresh(context.grid, context.soil_transport_faces, context.soil_face_geometry, context.soil_solver_properties.matrix_bulk_volume_m3, 1 - phosphate_band_fraction, phosphate_band_fraction, 1, .{});
    const soil_solute_inputs: ecosys.transport_step.SoilSoluteInputs = .{
        .micropore_diffusive_conductance_m3_per_step = context.soil_solute_face_parameters.micropore_conductance_m3_per_step,
        .macropore_diffusive_conductance_m3_per_step = context.soil_solute_face_parameters.macropore_conductance_m3_per_step,
        .micropore_mobility_fraction = context.soil_solute_face_parameters.micropore_mobility_fraction,
        .macropore_mobility_fraction = context.soil_solute_face_parameters.macropore_mobility_fraction,
        .layer_volume_m3 = context.soil_solver_properties.layer_volume_m3,
        .maximum_convective_fraction = 1,
        .pore_exchange_step_fraction = 1,
        .boundary_mobility_fraction = context.soil_solute_face_parameters.boundary_mobility_fraction,
        .recharge_concentration_mol_per_m3 = context.soil_recharge_concentration_mol_per_m3,
        .boundary_net_flux_mol_by_cell = context.soil_solute_boundary_net_flux_mol,
        .micropore_face_flux_mol_by_component = context.micropore_solute_face_flux_mol,
        .macropore_face_flux_mol_by_component = context.macropore_solute_face_flux_mol,
        .solver_options = .{ .absolute_tolerance_mol = context.config.absolute_tolerance, .relative_tolerance = context.config.relative_tolerance, .picard_relaxation = context.config.picard_relaxation, .max_iterations = water_heat_solute_max_iterations },
    };
    _ = try ecosys.transport_step.solveSoilFaceFluxesDeferred(
        context.allocator,
        context.grid,
        context.soil_transport_faces,
        context.micropore_solute_state,
        context.macropore_solute_state,
        soil_solute_inputs,
    );
    try commitHourlySoilSoluteContributionGenerations(context);
    try ecosys.transport_step.advanceSoilLocalSoluteProcesses(
        context.allocator,
        context.grid,
        context.transport_hydrology,
        context.micropore_solute_state,
        context.macropore_solute_state,
        soil_solute_inputs,
    );
    try ecosys.soil_aqueous_transport_bridge.importChemistry(context.micropore_solute_state, context.soil_chemistry);
    try ecosys.subsurface_irrigation_chemistry.addPhosphate(
        context.irrigation_loads,
        context.soil_chemistry,
        context.grid.matrix_liquid_water_m3,
        subsurface_irrigation_chemistry_parameters,
    );
    try context.soil_organic_face_parameters.refresh(
        context.grid,
        context.soil_transport_faces,
        context.soil_face_geometry,
        context.soil_solver_properties.matrix_bulk_volume_m3,
        1,
        .{},
    );
    _ = try ecosys.soil_organic_transport.advance(
        context.allocator,
        context.soil_organic_transport,
        context.soil_organic,
        context.soil_transport_faces,
        .{
            .micropore_conductance_m3_per_step = context.soil_organic_face_parameters.micropore_conductance_m3_per_step,
            .macropore_conductance_m3_per_step = context.soil_organic_face_parameters.macropore_conductance_m3_per_step,
            .matrix_water_m3 = context.grid.matrix_liquid_water_m3,
            .macropore_water_m3 = context.grid.macropore_liquid_water_m3,
            .layer_bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3,
            .micropore_external_water_flux_m3_per_step = context.transport_hydrology.micropore_external_water_flux_m3_per_step,
            .macropore_external_water_flux_m3_per_step = context.transport_hydrology.macropore_external_water_flux_m3_per_step,
            .recharge_concentration_g_per_m3 = context.soil_organic_recharge_concentration_g_per_m3,
        },
        .{
            .absolute_tolerance_g = context.config.absolute_tolerance,
            .relative_tolerance = context.config.relative_tolerance,
            .picard_relaxation = context.config.picard_relaxation,
            .max_iterations = water_heat_solute_max_iterations,
        },
    );
    try ecosys.subsurface_irrigation_chemistry.addMineralNitrogen(
        context.irrigation_loads,
        context.mineral_nitrogen_transport,
        subsurface_irrigation_chemistry_parameters,
    );
    try context.mineral_nitrogen_face_parameters.refresh(
        context.grid,
        context.soil_transport_faces,
        context.soil_face_geometry,
        context.soil_solver_properties.matrix_bulk_volume_m3,
        context.mineral_nitrogen_zone_fractions,
        context.runscript.root_nutrient_parameters,
        1,
    );
    _ = try ecosys.mineral_nitrogen_transport.advance(context.allocator, context.mineral_nitrogen_transport, .{
        .matrix_water_volume_m3 = context.grid.matrix_liquid_water_m3,
        .macropore_water_volume_m3 = context.grid.macropore_liquid_water_m3,
        .layer_volume_m3 = context.soil_solver_properties.layer_volume_m3,
        .matrix_faces = context.soil_transport_faces.micropore_faces,
        .macropore_faces = context.soil_transport_faces.macropore_faces,
        .matrix_conductance_m3_per_step = context.mineral_nitrogen_face_parameters.matrix_conductance_m3_per_step,
        .macropore_conductance_m3_per_step = context.mineral_nitrogen_face_parameters.macropore_conductance_m3_per_step,
        .mobility_fraction = context.mineral_nitrogen_face_parameters.mobility_fraction,
        .matrix_external_water_flux_m3_per_step = context.transport_hydrology.micropore_external_water_flux_m3_per_step,
        .macropore_external_water_flux_m3_per_step = context.transport_hydrology.macropore_external_water_flux_m3_per_step,
        .maximum_convective_fraction = 1,
        .pore_exchange_step_fraction = 1,
        .nitrogen_molar_mass_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
        .solver_options = .{
            .absolute_tolerance_mol = context.config.absolute_tolerance,
            .relative_tolerance = context.config.relative_tolerance,
            .picard_relaxation = context.config.picard_relaxation,
            .max_iterations = water_heat_solute_max_iterations,
        },
    });
    try context.mineral_nitrogen_transport.publishMatrix(
        context.soil_chemistry,
        context.soil_reactive_nitrogen,
        context.mineral_nitrogen_zone_fractions,
        context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
    );
    try context.soil_dissolved_gas_face_parameters.refresh(
        context.grid,
        context.soil_transport_faces,
        context.soil_face_geometry,
        context.soil_solver_properties.matrix_bulk_volume_m3,
        1,
        .{},
    );
    _ = try ecosys.soil_dissolved_gas_transport.advance(
        context.allocator,
        context.soil_dissolved_gas_transport,
        context.gas_transport,
        .{
            .faces = context.soil_transport_faces,
            .micropore_conductance_m3_per_step = context.soil_dissolved_gas_face_parameters.micropore_conductance_m3_per_step,
            .macropore_conductance_m3_per_step = context.soil_dissolved_gas_face_parameters.macropore_conductance_m3_per_step,
            .micropore_water_m3 = context.grid.matrix_liquid_water_m3,
            .macropore_water_m3 = context.grid.macropore_liquid_water_m3,
            .layer_bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3,
            .micropore_external_water_flux_m3_per_step = context.transport_hydrology.micropore_external_water_flux_m3_per_step,
            .macropore_external_water_flux_m3_per_step = context.transport_hydrology.macropore_external_water_flux_m3_per_step,
            .recharge_concentration_g_per_m3 = context.soil_dissolved_gas_recharge_concentration_g_per_m3,
        },
        .{
            .absolute_tolerance = context.config.absolute_tolerance,
            .relative_tolerance = context.config.relative_tolerance,
            .picard_relaxation = context.config.picard_relaxation,
            .max_iterations = context.iteration_limits.gas_max_iterations,
        },
    );
    if (context.surface_gas_parameters.*) |surface_parameters| {
        for (0..context.grid.cell_count) |cell| {
            const snow_first = cell * context.snow_transport.layer_capacity;
            const snow_end = snow_first + context.snow_transport.layer_capacity;
            for (snow_first..snow_end) |snow_layer| {
                const snow_temperature_k = if (context.snow_transport.temperature_k[snow_layer] > 0) context.snow_transport.temperature_k[snow_layer] else context.atmosphere.air_temperature_k[cell];
                context.snow_layer_gas_diffusivity_m2_per_h[snow_layer] =
                    context.runscript.snow_vapor_diffusion_parameters.reference_vapor_diffusivity_m2_per_h *
                    std.math.pow(f64, snow_temperature_k / context.runscript.snow_vapor_diffusion_parameters.reference_temperature_k, context.runscript.snow_vapor_diffusion_parameters.temperature_exponent);
            }
            const area_m2 = context.canopy_cell_area_m2[cell];
            const litter_volume_m3 = context.surface_litter_geometry.expanded_total_volume_m3[cell];
            const litter_air_m3 = context.surface_litter_geometry.air_volume_m3[cell];
            const litter_porosity = context.surface_litter_geometry.porosity_m3_per_m3[cell];
            const atmospheric_diffusivity_m2_per_h =
                context.runscript.soil_process_parameters.reference_water_vapor_diffusivity_m2_per_h *
                std.math.pow(f64, context.atmosphere.air_temperature_k[cell] / context.runscript.soil_process_parameters.vapor_diffusivity_reference_temperature_k, context.runscript.soil_process_parameters.vapor_diffusivity_temperature_exponent);
            const litter_porous_resistance_h_per_m = if (litter_volume_m3 > 0 and litter_air_m3 > 0 and litter_porosity > 0) blk: {
                const air_fraction = std.math.clamp(litter_air_m3 / litter_volume_m3, 0, litter_porosity);
                const transport_factor = @max(
                    context.runscript.surface_gas_resistance_parameters.minimum_air_fraction,
                    context.runscript.soil_gas_transport_parameters.penman_tortuosity * air_fraction * air_fraction / litter_porosity,
                );
                break :blk (litter_volume_m3 / area_m2) / atmospheric_diffusivity_m2_per_h / transport_factor;
            } else 0;
            const litter_fraction = std.math.clamp(context.surface_precipitation.litter_cover_fraction[cell], 0, 1);
            const resistance = try ecosys.surface_gas_boundary_conductance.calculate(.{
                .cell_area_m2 = area_m2,
                .air_temperature_k = context.atmosphere.air_temperature_k[cell],
                .ground_air_temperature_k = context.ground_air.temperature_k[cell],
                .surface_temperature_k = context.grid.surface_temperature_k[cell],
                .bulk_richardson_coefficient_k = context.surface_aerodynamics.bulk_richardson_coefficient_k[cell],
                .isothermal_atmospheric_resistance_h_per_m = context.surface_aerodynamics.isothermal_aerodynamic_resistance_h_per_m[cell],
                .total_canopy_area_m2 = context.surface_total_canopy_area_m2[cell],
                .canopy_height_m = context.surface_canopy_height_m[cell],
                .atmospheric_vapor_diffusivity_m2_per_h = atmospheric_diffusivity_m2_per_h,
                .isothermal_ground_surface_resistance_h_per_m = @max(context.runscript.surface_gas_resistance_parameters.minimum_aerodynamic_resistance_h_per_m, context.runscript.canopy_sensible_surface_resistance_h_per_m),
                .bare_surface_fraction = 1 - litter_fraction,
                .litter_surface_fraction = litter_fraction,
                .litter_porous_resistance_h_per_m = litter_porous_resistance_h_per_m,
                .snow_layer_thickness_m = context.snow_transport.layer_thickness_m[snow_first..snow_end],
                .snow_layer_total_volume_m3 = context.snow_transport.total_layer_volume_m3[snow_first..snow_end],
                .snow_layer_air_volume_m3 = context.snow_transport.air_filled_volume_m3[snow_first..snow_end],
                .snow_layer_vapor_diffusivity_m2_per_h = context.snow_layer_gas_diffusivity_m2_per_h[snow_first..snow_end],
            }, context.runscript.surface_gas_resistance_parameters);
            context.litter_atmospheric_gas_conductance_m3_per_h[cell] = resistance.atmospheric_litter_gas_conductance_m3_per_h;
            context.soil_atmospheric_gas_conductance_m3_per_h[cell] = resistance.atmospheric_gas_conductance_m3_per_h;
        }
        var accepted_gas_state = try context.gas_transport.clone(
            context.allocator,
        );
        defer accepted_gas_state.deinit();
        _ = try context.soil_gas_transport.advance(.{
            .grid = context.grid,
            .hydrology = context.transport_hydrology,
            .soil_faces = context.soil_transport_faces,
            .geometry = context.soil_face_geometry,
            .matrix_bulk_volume_m3 = context.soil_solver_properties.matrix_bulk_volume_m3,
            .total_porosity_fraction = context.soil_solver_properties.porosity_fraction,
            .field_capacity_fraction = context.soil_field_capacity_fraction,
            .gas_state = &accepted_gas_state,
            .solubility_parameters = surface_parameters.solubility,
            .exchange_parameters = surface_parameters.exchange,
            .ammonium_band_fraction = context.runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
            .surface_boundary_inputs = .{
                .atmospheric_conductance_m3_per_step = context.soil_atmospheric_gas_conductance_m3_per_h,
                .cell_area_m2 = context.canopy_cell_area_m2,
                .top_layer_thickness_m = context.soil_solver_properties.layer_thickness_m,
                .atmospheric_concentration_g_per_m3 = context.current_atmospheric_gas_concentration_g_per_m3.*,
            },
            .subsurface_boundary_inputs = .{
                .topology = context.soil_boundary_topology,
                .layer_thickness_m = context.soil_solver_properties.layer_thickness_m,
                .external_concentration_g_per_m3 = context.current_atmospheric_gas_concentration_g_per_m3.*,
            },
            .parameters = context.runscript.soil_gas_transport_parameters,
            .solver_options = .{ .absolute_tolerance_g = context.config.absolute_tolerance, .relative_tolerance = context.config.relative_tolerance, .picard_relaxation = context.config.picard_relaxation, .transport_iteration_fraction = 1, .max_iterations = context.iteration_limits.gas_max_iterations },
            .failure_report = gas_failure_report,
        });
        try commitHourlyGasContributionGeneration(
            context,
            &accepted_gas_state,
        );
    }
    try ecosys.surface_precipitation.commitRuntimeIngress(context.surface_precipitation, 1);
    // WATSUB:137 refreshes ALTG from immutable site altitude and the current
    // REDIST surface boundary before any surface hydraulic-head comparison.
    try context.terrain_hydrology.refreshCurrentSurfaceElevations(context.soil_geometry);
    for (0..context.grid.cell_count) |cell| {
        const area_m2 = context.canopy_cell_area_m2[cell];
        const recovered_energy_j_per_m2 = try ecosys.surface_precipitation.recoverRainfallImpactEnergy(
            context.surface_precipitation.cumulative_rainfall_impact_energy_j[cell],
            1,
            context.runscript.rainfall_impact_parameters.conductivity_recovery_fraction_per_h,
        );
        context.surface_precipitation.cumulative_rainfall_impact_energy_j[cell] = recovered_energy_j_per_m2;
        const top_layer = cell * context.grid.soil_layer_capacity;
        const sand_fraction = context.soil_solver_properties.sand_mass_fraction[top_layer];
        const clay_fraction = context.soil_solver_properties.clay_mass_fraction[top_layer];
        const silt_fraction = @max(0, 1 - sand_fraction - clay_fraction);
        const total_rainfall_mm_per_h = context.surface_precipitation.rainfall_m3_per_h[cell] / area_m2 * 1.0e3;
        if (total_rainfall_mm_per_h <= 0) {
            context.surface_precipitation.rainfall_impact_energy_j[cell] = 0;
            const multiplier = try ecosys.surface_precipitation.rainfallConductivityMultiplier(
                recovered_energy_j_per_m2,
                silt_fraction,
                clay_fraction,
                context.runscript.rainfall_impact_parameters.conductivity_damage_per_j_per_megagram_per_megagram,
            );
            context.surface_precipitation.saturated_hydraulic_conductivity_multiplier[cell] = multiplier;
            context.soil_solver_properties.rainfall_conductivity_multiplier[top_layer] = multiplier;
            continue;
        }
        const throughfall_m3_per_h = if (context.canopy_precipitation_retention.*) |*retention|
            retention.cell_throughfall_m3_per_h[cell]
        else
            context.surface_precipitation.rainfall_m3_per_h[cell];
        const canopy_exposure = if (context.canopy_exposure.*) |*exposure|
            std.math.clamp(exposure.canopy_exposure_fraction[cell], 0, 1)
        else
            0;
        const canopy_throughfall_mm_per_h = throughfall_m3_per_h * canopy_exposure / area_m2 * 1.0e3;
        const direct_precipitation_mm_per_h = throughfall_m3_per_h * (1 - canopy_exposure) / area_m2 * 1.0e3;
        const impact = try ecosys.surface_precipitation.rainfallImpact(
            recovered_energy_j_per_m2,
            .{
                .direct_precipitation_mm_per_h = direct_precipitation_mm_per_h,
                .throughfall_mm_per_h = canopy_throughfall_mm_per_h,
                .total_precipitation_mm_per_h = total_rainfall_mm_per_h,
                .canopy_height_m = context.surface_canopy_height_m[cell],
                .excess_surface_storage_m3 = @max(0, context.surface_precipitation.litter_water_m3[cell] + context.surface_litter_ice_m3[cell] - context.surface_precipitation.litter_water_capacity_m3[cell]),
                .ground_surface_retention_m3 = context.runscript.surface_runoff_parameters.ground_surface_retention_m3_per_m2 * area_m2,
                .surface_area_m2 = area_m2,
                .bare_soil_fraction = 1 - std.math.clamp(context.surface_precipitation.litter_cover_fraction[cell], 0, 1),
                .time_fraction = 1,
                .surface_silt_megagrams_per_megagram = silt_fraction,
                .surface_clay_megagrams_per_megagram = clay_fraction,
            },
            context.runscript.rainfall_impact_parameters,
        );
        context.surface_precipitation.rainfall_impact_energy_j[cell] = impact.incremental_energy_j;
        context.surface_precipitation.cumulative_rainfall_impact_energy_j[cell] = impact.cumulative_energy_j;
        context.surface_precipitation.saturated_hydraulic_conductivity_multiplier[cell] = impact.saturated_conductivity_multiplier;
        context.soil_solver_properties.rainfall_conductivity_multiplier[top_layer] = impact.saturated_conductivity_multiplier;
    }
    try ecosys.surface_runoff.route(
        context.surface_runoff,
        context.config.grid_columns,
        context.config.grid_rows,
        context.terrain_hydrology,
        context.canopy_cell_area_m2,
        context.surface_precipitation.litter_water_m3,
        context.surface_litter_ice_m3,
        context.surface_precipitation.litter_water_capacity_m3,
        context.lateral_connection_mode_by_cell,
        .{
            .north = context.surface_runoff_boundary_fraction_by_direction[0],
            .east = context.surface_runoff_boundary_fraction_by_direction[1],
            .south = context.surface_runoff_boundary_fraction_by_direction[2],
            .west = context.surface_runoff_boundary_fraction_by_direction[3],
        },
        context.runscript.surface_runoff_parameters,
    );
    try ecosys.surface_mineral_transport.advance(
        context.allocator,
        context.surface_litter_chemistry,
        context.surface_denitrification.nitrite_g_n,
        context.config.grid_columns,
        context.config.grid_rows,
        context.surface_precipitation.litter_water_m3,
        context.surface_runoff.water_change_m3,
        .{
            .east_m3 = context.surface_runoff.east_runoff_m3_per_step,
            .west_m3 = context.surface_runoff.west_runoff_m3_per_step,
            .south_m3 = context.surface_runoff.south_runoff_m3_per_step,
            .north_m3 = context.surface_runoff.north_runoff_m3_per_step,
        },
        1,
        context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
        context.runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
        .{
            .inorganic_nitrogen_export_g_n_by_cell = context.surface_inorganic_nitrogen_export_g_n_per_h,
            .inorganic_phosphorus_export_g_p_by_cell = context.surface_inorganic_phosphorus_export_g_p_per_h,
        },
    );
    try ecosys.surface_organic_transport.advance(
        context.allocator,
        context.surface_organic,
        context.config.grid_columns,
        context.config.grid_rows,
        context.surface_precipitation.litter_water_m3,
        context.surface_runoff.water_change_m3,
        .{
            .east_m3 = context.surface_runoff.east_runoff_m3_per_step,
            .west_m3 = context.surface_runoff.west_runoff_m3_per_step,
            .south_m3 = context.surface_runoff.south_runoff_m3_per_step,
            .north_m3 = context.surface_runoff.north_runoff_m3_per_step,
        },
        1,
        .{
            .dissolved_organic_carbon_export_g_c_by_cell = context.surface_organic_carbon_export_g_c_per_h,
            .dissolved_organic_nitrogen_export_g_n_by_cell = context.surface_organic_nitrogen_export_g_n_per_h,
            .dissolved_organic_phosphorus_export_g_p_by_cell = context.surface_organic_phosphorus_export_g_p_per_h,
        },
    );
    try ecosys.surface_dissolved_gas_transport.advance(
        context.allocator,
        context.litter_gas_transport,
        context.config.grid_columns,
        context.config.grid_rows,
        context.surface_precipitation.litter_water_m3,
        context.surface_runoff.water_change_m3,
        .{
            .east_m3 = context.surface_runoff.east_runoff_m3_per_step,
            .west_m3 = context.surface_runoff.west_runoff_m3_per_step,
            .south_m3 = context.surface_runoff.south_runoff_m3_per_step,
            .north_m3 = context.surface_runoff.north_runoff_m3_per_step,
        },
        1,
        context.surface_inorganic_carbon_export_g_c_per_h,
    );
    for (0..context.grid.cell_count) |cell| {
        const top_layer = cell * context.grid.soil_layer_capacity;
        const sand_fraction = context.soil_solver_properties.sand_mass_fraction[top_layer];
        const clay_fraction = context.soil_solver_properties.clay_mass_fraction[top_layer];
        const silt_fraction = @max(0, 1 - sand_fraction - clay_fraction);
        const humus_fraction = std.math.clamp(1.82e-6 * context.soil_solver_properties.total_organic_carbon_g_per_megagram[top_layer], 0, 1);
        var root_length_density_m_per_m3: f64 = 0;
        if (context.plant_roots.*) |*roots| if (context.plant_water_workspace.*) |*workspace| {
            for (0..context.config.plant_populations) |species| {
                const plant = cell * context.config.plant_populations + species;
                for (0..context.root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| {
                    root_length_density_m_per_m3 += roots.root_length_density_m_per_m3[try roots.layerIndex(plant, domain, 0)] * workspace.plant_population_count[plant] / context.canopy_cell_area_m2[cell];
                }
            }
        };
        const erosion_properties = try ecosys.soil_erosion.deriveSurfaceProperties(
            .{ .sand_mass_fraction = sand_fraction, .silt_mass_fraction = silt_fraction, .clay_mass_fraction = clay_fraction, .humus_mass_fraction = humus_fraction, .residue_mass_fraction = 0, .root_length_density_m_per_m3 = root_length_density_m_per_m3, .surface_temperature_c = context.grid.surface_temperature_k[cell] - 273.15 },
            .{
                .reference_water_viscosity_Mg_per_m_s = context.runscript.soil_process_parameters.reference_water_viscosity_megagrams_per_m_s,
                .viscosity_temperature_intercept = context.runscript.soil_process_parameters.water_viscosity_temperature_intercept,
                .viscosity_temperature_coefficient_per_c = context.runscript.soil_process_parameters.water_viscosity_temperature_coefficient_per_c,
            },
        );
        const soil_mass_Mg = context.soil_solver_properties.matrix_bulk_volume_m3[top_layer] * context.soil_solver_properties.bulk_density_megagrams_per_m3[top_layer];
        if (!context.surface_erosion.surface_soil_mass_initialized[cell]) {
            context.surface_erosion.surface_soil_mass_Mg[cell] = soil_mass_Mg;
            context.surface_erosion.surface_soil_mass_initialized[cell] = true;
        }
        context.surface_soil_mass_at_erosion_start_Mg[cell] = context.surface_erosion.surface_soil_mass_Mg[cell];
        const local_solve =
            try ecosys.soil_erosion.calculateConvergedHourlyLocalStep(.{
                .erosion_enabled = context.site_by_cell[cell].erosionEnabled(),
                .surface_soil_bulk_density_Mg_per_m3 = context.soil_solver_properties.bulk_density_megagrams_per_m3[top_layer],
                .surface_soil_mass_Mg = context.surface_erosion.surface_soil_mass_Mg[cell],
                .surface_soil_water_m3 = context.grid.matrix_liquid_water_m3[top_layer],
                .surface_soil_pore_volume_m3 = context.grid.matrix_pore_capacity_m3[top_layer],
                .excess_surface_water_m3 = context.surface_runoff.excess_surface_water_m3[cell],
                .excess_surface_ice_m3 = context.surface_runoff.excess_surface_ice_m3[cell],
                .surface_ponding_capacity_m3 = @max(context.config.absolute_tolerance, context.runscript.surface_runoff_parameters.ground_surface_retention_m3_per_m2 * context.canopy_cell_area_m2[cell]),
                .sediment_in_surface_water_Mg = context.surface_erosion.surface_sediment_Mg[cell],
                .rainfall_kinetic_energy_j = context.surface_precipitation.rainfall_impact_energy_j[cell],
                .soil_rainfall_detachability_g_per_j = erosion_properties.rainfall_detachability_g_per_j,
                .soil_runoff_detachability = erosion_properties.runoff_detachability,
                .sediment_settling_velocity_m_per_h = erosion_properties.settling_velocity_m_per_h,
                .grid_cell_area_m2 = context.canopy_cell_area_m2[cell],
                .soil_matrix_fraction = std.math.clamp(context.soil_solver_properties.matrix_bulk_volume_m3[top_layer] / context.soil_solver_properties.layer_volume_m3[top_layer], 0, 1),
                .snow_free_fraction = 1 - std.math.clamp(context.surface_precipitation.snow_cover_fraction[cell], 0, 1),
                .runoff_velocity_m_per_s = context.surface_runoff.runoff_velocity_m_per_s[cell],
                .slope_sine = context.terrain_hydrology.slope_m_per_m[cell],
                .surface_particle_density_Mg_per_m3 = erosion_properties.particle_density_Mg_per_m3,
                .transport_capacity_coefficient = erosion_properties.transport_capacity_coefficient,
                .transport_capacity_exponent = erosion_properties.transport_capacity_exponent,
                .maximum_erodible_soil_fraction_per_step = 1,
                .water_transport_timestep_h = 1,
                .negligible_volume_m3 = context.runscript.surface_runoff_parameters.negligible_water_m3,
                .negligible_mass_Mg = context.runscript.surface_runoff_parameters.negligible_water_m3,
            }, .{
                .absolute_tolerance_Mg = context.config.absolute_tolerance,
                .relative_tolerance = context.config.relative_tolerance,
                .picard_relaxation = context.config.picard_relaxation,
                .max_iterations = context.iteration_limits.erosion_max_iterations,
            });
        const local = local_solve.local;
        context.surface_erosion.local_detachment_Mg[cell] = local.net_detachment_Mg;
        context.surface_erosion.transportable_sediment_Mg[cell] = if (context.site_by_cell[cell].erosionEnabled())
            try ecosys.soil_erosion.calculateDownslopeTransport_Mg(
                context.surface_erosion.surface_sediment_Mg[cell],
                local.net_detachment_Mg,
                context.surface_runoff.excess_surface_water_m3[cell],
                context.surface_runoff.total_runoff_m3_per_step[cell],
                context.soil_solver_properties.bulk_density_megagrams_per_m3[top_layer],
                context.runscript.surface_runoff_parameters.negligible_water_m3,
            )
        else
            0;
        const column = cell % context.config.grid_columns;
        const row = cell / context.config.grid_columns;
        context.surface_erosion.east_boundary_open[cell] = column + 1 == context.config.grid_columns and context.site_by_cell[cell].surface_runoff_boundary_fraction[1] > 0;
        context.surface_erosion.west_boundary_open[cell] = column == 0 and context.site_by_cell[cell].surface_runoff_boundary_fraction[3] > 0;
        context.surface_erosion.south_boundary_open[cell] = row + 1 == context.config.grid_rows and context.site_by_cell[cell].surface_runoff_boundary_fraction[2] > 0;
        context.surface_erosion.north_boundary_open[cell] = row == 0 and context.site_by_cell[cell].surface_runoff_boundary_fraction[0] > 0;
    }
    try ecosys.sediment_routing.route(
        &context.surface_erosion.routing,
        context.surface_erosion.transportable_sediment_Mg,
        context.surface_runoff.total_runoff_m3_per_step,
        .{ .east_m3 = context.surface_runoff.east_runoff_m3_per_step, .west_m3 = context.surface_runoff.west_runoff_m3_per_step, .south_m3 = context.surface_runoff.south_runoff_m3_per_step, .north_m3 = context.surface_runoff.north_runoff_m3_per_step },
        .{
            .east_open = context.surface_erosion.east_boundary_open,
            .west_open = context.surface_erosion.west_boundary_open,
            .south_open = context.surface_erosion.south_boundary_open,
            .north_open = context.surface_erosion.north_boundary_open,
            .lateral_connection_mode_by_cell = context.lateral_connection_mode_by_cell,
        },
        context.runscript.surface_runoff_parameters.negligible_water_m3,
    );
    try ecosys.sediment_routing.commitSurfaceSediment(context.surface_erosion.surface_sediment_Mg, context.surface_erosion.local_detachment_Mg, context.surface_erosion.routing.sediment_change_Mg);
    try ecosys.soil_erosion_organic_bridge.route(
        context.config.grid_columns,
        context.config.grid_rows,
        context.grid.soil_layer_capacity,
        context.surface_erosion.surface_soil_mass_Mg,
        context.soil_organic,
        .{
            .east_Mg = context.surface_erosion.routing.east_flux_Mg,
            .west_Mg = context.surface_erosion.routing.west_flux_Mg,
            .south_Mg = context.surface_erosion.routing.south_flux_Mg,
            .north_Mg = context.surface_erosion.routing.north_flux_Mg,
        },
        context.eroded_organic_workspace,
    );
    try ecosys.soil_erosion_organic_bridge.refreshSurfaceOrganicCarbonGPerMg(
        context.soil_organic,
        context.grid.soil_layer_capacity,
        context.surface_erosion.surface_soil_mass_Mg,
        context.soil_solver_properties.total_organic_carbon_g_per_megagram,
    );
    try ecosys.soil_erosion_fertilizer_bridge.route(
        context.config.grid_columns,
        context.config.grid_rows,
        context.surface_erosion.surface_soil_mass_Mg,
        context.soil_fertilizer_inventory,
        .{
            .east_Mg = context.surface_erosion.routing.east_flux_Mg,
            .west_Mg = context.surface_erosion.routing.west_flux_Mg,
            .south_Mg = context.surface_erosion.routing.south_flux_Mg,
            .north_Mg = context.surface_erosion.routing.north_flux_Mg,
        },
        context.eroded_fertilizer_workspace,
    );
    try ecosys.soil_erosion_chemistry_bridge.route(
        context.config.grid_columns,
        context.config.grid_rows,
        context.grid.soil_layer_capacity,
        context.surface_erosion.surface_soil_mass_Mg,
        context.grid.matrix_liquid_water_m3,
        context.soil_chemistry,
        .{
            .east_Mg = context.surface_erosion.routing.east_flux_Mg,
            .west_Mg = context.surface_erosion.routing.west_flux_Mg,
            .south_Mg = context.surface_erosion.routing.south_flux_Mg,
            .north_Mg = context.surface_erosion.routing.north_flux_Mg,
        },
        context.eroded_chemistry_workspace,
    );
    try ecosys.soil_erosion_mineral_bridge.route(
        context.config.grid_columns,
        context.config.grid_rows,
        context.grid.soil_layer_capacity,
        context.surface_erosion.surface_soil_mass_Mg,
        context.soil_solver_properties,
        .{
            .east_Mg = context.surface_erosion.routing.east_flux_Mg,
            .west_Mg = context.surface_erosion.routing.west_flux_Mg,
            .south_Mg = context.surface_erosion.routing.south_flux_Mg,
            .north_Mg = context.surface_erosion.routing.north_flux_Mg,
        },
        context.eroded_mineral_state,
    );
    // Mineral routing is the authoritative commit of transported surface
    // soil mass. Rebind the organic concentration to that accepted carrier
    // immediately; the extensive organic pools were routed above.
    try ecosys.soil_erosion_organic_bridge.refreshSurfaceOrganicCarbonGPerMg(
        context.soil_organic,
        context.grid.soil_layer_capacity,
        context.surface_erosion.surface_soil_mass_Mg,
        context.soil_solver_properties.total_organic_carbon_g_per_megagram,
    );
    const eroded_organic_export =
        try ecosys.soil_erosion_organic_bridge.exportedElements(
            context.soil_organic,
            context.eroded_organic_workspace,
        );
    const eroded_fertilizer_export =
        try ecosys.soil_erosion_fertilizer_bridge.exported(
            context.eroded_fertilizer_workspace,
            context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
        );
    const eroded_chemistry_export =
        try ecosys.soil_erosion_chemistry_bridge.exported(
            context.eroded_chemistry_workspace,
            12,
            context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
            context.runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
        );
    try context.landscape_boundary_ledger.accumulateAccepted(.{
        .carbon_output_g_c = eroded_organic_export.carbon_g_c +
            eroded_chemistry_export.inorganic_carbon_g_c,
        .nitrogen_output_g_n = eroded_organic_export.nitrogen_g_n +
            eroded_fertilizer_export.nitrogen_g_n +
            eroded_chemistry_export.nitrogen_g_n,
        .phosphorus_output_g_p = eroded_organic_export.phosphorus_g_p +
            eroded_chemistry_export.phosphorus_g_p,
        .ion_output_mol = eroded_fertilizer_export.ion_mol +
            eroded_chemistry_export.ion_mol,
    });
    try ecosys.soil_sediment_change.publishAcceptedNetSedimentMg(context.surface_soil_mass_at_erosion_start_Mg, context.surface_erosion.surface_soil_mass_Mg, context.net_sediment_Mg_per_h);
    @memcpy(context.transport_hydrology.runoff_total_m3_per_step, context.surface_runoff.total_runoff_m3_per_step);
    @memcpy(context.transport_hydrology.runoff_east_m3_per_step, context.surface_runoff.east_runoff_m3_per_step);
    @memcpy(context.transport_hydrology.runoff_west_m3_per_step, context.surface_runoff.west_runoff_m3_per_step);
    @memcpy(context.transport_hydrology.runoff_south_m3_per_step, context.surface_runoff.south_runoff_m3_per_step);
    @memcpy(context.transport_hydrology.runoff_north_m3_per_step, context.surface_runoff.north_runoff_m3_per_step);
    for (context.snow_surface_discharge, context.direct_surface_solute_input) |*snow_discharge, direct_input| {
        for (&snow_discharge.litter_g, direct_input.litter_g) |*destination, amount| destination.* += amount;
        for (&snow_discharge.soil_nonband_g, direct_input.soil_nonband_g) |*destination, amount| destination.* += amount;
        for (&snow_discharge.soil_band_g, direct_input.soil_band_g) |*destination, amount| destination.* += amount;
    }
    const ion_molar_masses: ecosys.snow_surface_discharge.IonMolarMassesGPerMol = if (context.runscript.chemistry_primary_initialization) |parameters|
        .{ .aluminum = parameters.molar_mass_g_per_mol.aluminum, .iron = parameters.molar_mass_g_per_mol.iron, .calcium = parameters.molar_mass_g_per_mol.calcium, .magnesium = parameters.molar_mass_g_per_mol.magnesium, .sodium = parameters.molar_mass_g_per_mol.sodium, .potassium = parameters.molar_mass_g_per_mol.potassium, .sulfur = parameters.molar_mass_g_per_mol.sulfur, .chloride = parameters.molar_mass_g_per_mol.chloride }
    else
        .{ .aluminum = 27, .iron = 55.8, .calcium = 40, .magnesium = 24.3, .sodium = 23, .potassium = 39.1, .sulfur = 32, .chloride = 35.5 };
    try ecosys.snow_surface_discharge.commit(context.allocator, .{ .discharge = context.snow_surface_discharge, .litter_water_volume_m3 = context.surface_precipitation.litter_water_m3, .topsoil_water_volume_m3 = context.grid.matrix_liquid_water_m3[0..context.grid.layer_count], .soil_layer_capacity = context.grid.soil_layer_capacity, .ion_molar_mass_g_per_mol = ion_molar_masses }, context.litter_gas_transport, context.gas_transport, context.surface_litter_chemistry, context.soil_chemistry);
    try runSoilBiogeochemistryBySerialTile(context);
    if (context.surface_gas_parameters.*) |surface_parameters| {
        for (0..context.grid.cell_count) |cell| {
            context.surface_charcoal_carbon_g_c[cell] = try context.surface_organic.charcoalCarbon_g_c(cell);
        }
        var litter_geometry_context: ecosys.surface_litter_geometry_step.ApplyContext = .{
            .result = context.surface_litter_geometry,
            .surface_organic = context.surface_organic,
            .water_m3 = context.surface_precipitation.litter_water_m3,
            .ice_m3 = context.surface_litter_ice_m3,
            .charcoal_carbon_g_c = context.surface_charcoal_carbon_g_c,
            .parameters = surface_parameters.litter_geometry,
        };
        try runKernelAcrossSerialTiles(context, &litter_geometry_context, ecosys.surface_litter_geometry_step.applyTile);
        @memcpy(context.surface_precipitation.litter_water_capacity_m3, context.surface_litter_geometry.water_retention_capacity_m3);
        for (0..context.grid.cell_count) |cell| {
            context.litter_gas_transport.air_volume_m3[cell] = context.surface_litter_geometry.air_volume_m3[cell];
            context.litter_gas_transport.temperature_k[cell] = context.grid.surface_temperature_k[cell];
        }
        _ = try context.surface_litter_gas_transport.advance(
            context.litter_gas_transport,
            context.surface_litter_geometry,
            context.surface_precipitation.litter_water_m3,
            context.surface_litter_ice_m3,
            context.canopy_cell_area_m2,
            context.litter_atmospheric_gas_conductance_m3_per_h,
            context.current_atmospheric_gas_concentration_g_per_m3.*,
            surface_parameters.solubility,
            surface_parameters.exchange,
            .{
                .reference_temperature_k = context.runscript.soil_gas_transport_parameters.reference_temperature_k,
                .temperature_exponent = context.runscript.soil_gas_transport_parameters.temperature_exponent,
                .free_air_diffusivity_m2_per_h = context.runscript.soil_gas_transport_parameters.free_air_diffusivity_m2_per_h,
                .penman_tortuosity = context.runscript.soil_gas_transport_parameters.penman_tortuosity,
                .minimum_air_fraction = context.runscript.surface_gas_resistance_parameters.minimum_air_fraction,
            },
            .{
                .absolute_tolerance_g = context.config.absolute_tolerance,
                .relative_tolerance = context.config.relative_tolerance,
                .picard_relaxation = context.config.picard_relaxation,
                .max_iterations = context.iteration_limits.gas_max_iterations,
            },
        );
        @memset(context.surface_microbial_substrate_uptake.denitrification_respiration_g_c, 0);
        if (context.chemistry_reaction_parameters.* != null and context.organic_parameters.* != null)
            context.surface_litter_fertilizer_diagnostics.reset();
        try runSurfaceBiogeochemistryBySerialTile(context, surface_parameters);
    }
    // Legacy SOLUTE follows the soil and litter biological source/sink commits.
    // Converge locally once; do not repeat a full sub-hourly model cycle.
    try convergeHourlySoilChemistry(
        context,
        fertilizer_band_hour,
        solute_failure_report,
    );
    try convergeSurfaceLitterChemistry(context);
    var pond_transition_context: ecosys.surface_pond_transition_step.ApplyContext = .{
        .result = context.surface_pond_transition,
        .surface_liquid_water_m3 = context.surface_precipitation.litter_water_m3,
        .surface_ice_m3 = context.surface_litter_ice_m3,
        .surface_ponding_capacity_m3 = context.surface_precipitation.litter_water_capacity_m3,
        .surface_litter_volume_m3 = context.surface_litter_geometry.dry_litter_volume_m3,
        .surface_litter_water_capacity_m3 = context.surface_litter_geometry.water_retention_capacity_m3,
        .horizontal_area_m2 = context.canopy_cell_area_m2,
        .minimum_heat_capacity_mj_per_k = context.surface_pond_minimum_heat_capacity_mj_per_k,
        .liquid_water_heat_capacity_mj_per_m3_k = context.runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_mj_per_m3_k,
    };
    try runKernelAcrossSerialTiles(context, &pond_transition_context, ecosys.surface_pond_transition_step.applyTile);
    // REDIST lines 333-614: settle represented pond particulates after the
    // biological/chemical commits and before runoff-domain redistribution.
    try ecosys.surface_pond_particulate_settling.apply(.{
        .surface_organic = context.surface_organic,
        .soil_organic = context.soil_organic,
        .surface_gas = context.litter_gas_transport,
        .soil_gas = context.gas_transport,
        .surface_nitrogen_fertilizer = context.surface_litter_fertilizer,
        .soil_nitrogen_fertilizer = context.soil_fertilizer_inventory,
        .mineral_fertilizer = context.mineral_fertilizer_inventory,
    }, .{
        .surface_sediment_Mg = context.surface_erosion.surface_sediment_Mg,
        .surface_soil_mass_Mg = context.surface_erosion.surface_soil_mass_Mg,
        .settled_sediment_Mg = context.surface_erosion.pond_settled_sediment_Mg,
    }, context.surface_pond_transition, 1);
    try ecosys.surface_pond_domain_transaction.apply(context.surface_pond_domain_workspace, .{
        .inventories = .{
            .surface_organic = context.surface_organic,
            .soil_organic = context.soil_organic,
            .surface_gas = context.litter_gas_transport,
            .soil_gas = context.gas_transport,
            .surface_nitrogen_fertilizer = context.surface_litter_fertilizer,
            .soil_nitrogen_fertilizer = context.soil_fertilizer_inventory,
            .mineral_fertilizer = context.mineral_fertilizer_inventory,
        },
        .surface_chemistry = context.surface_litter_chemistry,
        .soil_chemistry = context.soil_chemistry,
        .water_heat = .{
            .surface_liquid_water_m3 = context.surface_precipitation.litter_water_m3,
            .surface_ice_m3 = context.surface_litter_ice_m3,
            .surface_temperature_k = context.grid.surface_temperature_k,
            .surface_geometry = context.surface_litter_geometry,
            .grid = context.grid,
            .soil_thermal = context.soil_thermal,
        },
        .soil_geometry = context.soil_geometry,
        .soil_properties = context.soil_solver_properties,
        .soil_faces = context.soil_transport_faces,
        .soil_face_geometry = context.soil_face_geometry,
    }, .{
        .transitions = context.surface_pond_transition,
        .dynamic_salts = context.runscript.dynamic_plant_salts,
        .water_heat_parameters = .{
            .dry_organic_heat_capacity_mj_per_g_c_k = context.runscript.surface_pond_dry_organic_heat_capacity_mj_per_g_c_k,
            .liquid_water_heat_capacity_mj_per_m3_k = context.runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_mj_per_m3_k,
            .ice_heat_capacity_mj_per_m3_k = context.runscript.soil_phase_heat_parameters.ice_heat_capacity_mj_per_m3_k,
            .minimum_heat_capacity_mj_per_k = 0,
        },
        .minimum_heat_capacity_mj_per_k = context.surface_pond_minimum_heat_capacity_mj_per_k,
        .minimum_soil_layer_thickness_m = context.runscript.soil_geometry_parameters.minimum_layer_thickness_m,
        .horizontal_cell_width_m = context.horizontal_cell_width_m,
        .vertical_cell_width_m = context.vertical_cell_width_m,
    });
    if (context.canopy_precipitation_retention.*) |*retention| {
        if (context.canopy_surface_exchange.*) |exchange|
            try ecosys.canopy_precipitation_retention.commitSurfaceWater(retention, exchange.intercepted_water_change_m3_per_h, context.standing_dead_evaporation_m3_per_h, 1)
        else
            try ecosys.canopy_precipitation_retention.commitRetention(retention, 1);
    }
    for (0..context.grid.cell_count) |cell| context.transport_hydrology.snow_surface_carrier_volume_m3[cell] = context.surface_precipitation.solid_snow_water_equivalent_m3[cell];
    // Gas inventories share the runtime layer topology but their capacities
    // follow the newly converged water/heat state each hour.
    for (0..context.grid.layer_count) |layer_cell| {
        context.gas_transport.temperature_k[layer_cell] = context.grid.soil_temperature_k[layer_cell];
        context.gas_transport.air_volume_m3[layer_cell] = context.grid.air_volume_m3[layer_cell];
    }
    if (context.canopy_energy.*) |*energy| {
        var canopy_energy_context: ecosys.canopy_energy.ApplyContext = .{ .result = energy, .plants = context.plants, .atmosphere = context.atmosphere, .exposure = &context.canopy_exposure.*.?, .interception = &context.canopy_interception.*.?, .canopy_longwave_emissivity = context.runscript.canopy_longwave_emissivity };
        try runKernelAcrossSerialTiles(context, &canopy_energy_context, ecosys.canopy_energy.applyTile);
    }
    if (context.detailed_canopy.*) |*canopy| {
        const angular_sample_count = try std.math.mul(usize, context.canopy_geometry.leaf_inclination_sine.len, context.canopy_geometry.leaf_azimuth_radians.len);
        const samples_per_node = try std.math.mul(usize, context.runscript.canopy_layer_count, angular_sample_count);
        try ensureCanopyGrowthNodeTopology(canopy, &context.plant_growth_stages.*.?, samples_per_node);
        var biochemistry_context: ecosys.canopy_biochemistry.ApplyContext = .{
            .canopy = canopy,
            .parameters_by_plant = context.canopy_biochemistry_parameters,
            .c4_carbon_parameters = context.runscript.c4_carbon_parameters,
            .canopy_temperature_k_by_plant = context.plants.canopy_temperature_k,
            .atmospheric_co2_umol_per_mol = context.current_atmospheric_co2_umol_per_mol.*,
            .dormancy = &context.plant_dormancy.*.?,
            .branch_development = &context.branch_development.*.?,
            .growth_stages = &context.plant_growth_stages.*.?,
            .dormancy_parameters_by_plant = context.development_dormancy_parameters,
            .stress_parameters = context.runscript.canopy_stress_parameters,
            .annual_termination_hours_without_grain_fill = context.runscript.root_metabolism_parameters.annual_termination_hours_without_grain_fill,
            .presence_threshold_g_per_plant = context.runscript.plant_pool_parameters.branch_structural_presence_g_per_plant,
            .timestep_h = 1,
        };
        try runKernelAcrossSerialTiles(context, &biochemistry_context, ecosys.canopy_biochemistry.applyTile);
    }
    if (context.plant_water_balance.*) |*balance| {
        const workspace = &context.plant_water_workspace.*.?;
        try ecosys.plant_water_balance.refreshCanopyWorkspace(workspace, &context.detailed_canopy.*.?, context.plants, context.config.plant_populations);
        try ecosys.plant_water_balance.refreshRootWorkspace(workspace, &context.plant_roots.*.?, &context.detailed_canopy.*.?, context.plants, context.grid, context.soil_solver_properties, context.config.plant_populations, context.root_biological_domain_count_by_plant, context.runscript.plant_geometry_parameters.root_volume_numerator_m3_per_g_c, context.runscript.plant_geometry_parameters.root_dry_matter_fraction, context.runscript.plant_geometry_parameters.root_pi, context.runscript.root_morphology_parameters);
        try workspace.refreshActive(context.config.soil_layers);
        var water_context: ecosys.plant_water_balance.ApplyContext = .{
            .result = balance,
            .grid = context.grid,
            .plants = context.plants,
            .soil_total_water_potential_mpa = context.soil_hourly_workspace.root_referenced_total_water_potential_mpa,
            .active = workspace.active,
            .root_conductance_m_per_h_mpa = workspace.root_conductance_m_per_h_mpa,
            .maximum_uptake_m = workspace.maximum_uptake_m,
            .maximum_release_m = workspace.maximum_release_m,
            .canopy_water_capacitance_m_per_m2_mpa = workspace.canopy_water_capacitance_m_per_m2_mpa,
            .transpiration_loss_m = workspace.transpiration_loss_m,
            .settings = .{ .minimum_canopy_water_potential_mpa = -100.0, .maximum_canopy_water_potential_mpa = 0.0, .solver_options = context.nonlinear_solver_options },
        };
        try runKernelAcrossSerialTiles(context, &water_context, ecosys.plant_water_balance.applyTile);
        try ecosys.plant_water_balance.commitRootHydraulics(balance.*, &context.plant_roots.*.?, context.grid, context.soil_hourly_workspace.root_referenced_total_water_potential_mpa, context.plants, workspace.active, workspace.cell_area_m2, workspace.soil_resistance_mpa_h_per_m, workspace.root_resistance_mpa_h_per_m, workspace.leaf_osmotic_potential_at_zero_total_mpa, context.root_biological_domain_count_by_plant);
        try ecosys.plant_water_balance.updateDailyMinimumCanopyWaterPotential(&context.detailed_canopy.*.?, context.plants);
        try ecosys.plant_root_gas_exchange.refreshOxygenDemand(&context.plant_roots.*.?, context.root_gas_parameters);
    }
    try applyPlantStorageRemobilization(context, plant_calendar);
    try applyRootMetabolism(context);
    if (context.plant_roots.*) |*roots| if (context.plant_water_workspace.*) |*water| if (context.plant_phenology.*) |phenology| {
        try ecosys.plant_root_gas_transport.advance(
            roots,
            water,
            context.grid,
            context.soil_solver_properties,
            context.gas_transport,
            context.root_gas_parameters,
            context.current_atmospheric_gas_concentration_g_per_m3.*,
            context.root_biological_domain_count_by_plant,
            phenology.active,
            .{
                .liquid_tortuosity_coefficient = context.runscript.root_nutrient_parameters.liquid_tortuosity_coefficient,
                .absolute_tolerance_g = context.config.absolute_tolerance,
                .relative_tolerance = context.config.relative_tolerance,
                .maximum_iterations = context.iteration_limits.gas_max_iterations,
            },
        );
    };
    try applyStorageExhaustionMortality(context, plant_calendar);
    try applyRootNutrientUptake(context);
    if (context.plant_growth_stages.*) |*growth_stages| {
        const roots = if (context.plant_roots.*) |*value| value else return error.MissingPlantRootsForPoolAggregation;
        const canopy = if (context.detailed_canopy.*) |*value| value else return error.MissingCanopyForPoolAggregation;
        var pool_context: ecosys.plant_pool_aggregation.ApplyContext = .{ .canopy = canopy, .roots = roots, .growth_stages = growth_stages, .active_by_plant = context.plant_phenology.*.?.active, .biological_domain_count_by_plant = context.root_biological_domain_count_by_plant, .dynamic_salts = context.runscript.dynamic_plant_salts, .parameters = context.runscript.plant_pool_parameters };
        try runKernelAcrossSerialTiles(context, &pool_context, ecosys.plant_pool_aggregation.applyTile);
    }
    if (context.plant_growth_stages.* != null and context.plant_roots.* != null and context.detailed_canopy.* != null and context.branch_development.* != null and context.plant_phenology.* != null) {
        const topology_states: ecosys.plant_topology.RuntimeStates = .{ .canopy = &context.detailed_canopy.*.?, .growth_stages = &context.plant_growth_stages.*.?, .dormancy = &context.plant_dormancy.*.?, .branch_development = &context.branch_development.*.? };
        const plant_execution_year = std.math.cast(u16, plant_calendar.current_year) orelse
            return error.InvalidPlantTopologyDate;
        const new_shoot_branches = try ecosys.plant_topology.advanceShootBranches(.{ .states = topology_states, .roots = &context.plant_roots.*.?, .controls = context.plant_topology_controls, .active_by_plant = context.plant_phenology.*.?.active, .emerged_by_plant = context.development_emerged, .day_of_year = plant_calendar.day_of_year, .execution_year = plant_execution_year, .minimum_root_turgor_potential_mpa = context.runscript.phenology_parameters.minimum_turgor_potential_mpa });
        if (new_shoot_branches > 0) {
            var replacement_layers = try ecosys.canopy_layer_distribution.State.init(context.allocator, context.grid.cell_count, context.config.plant_populations, context.runscript.canopy_layer_count, context.canopy_geometry.leaf_inclination_sine.len, context.canopy_geometry.leaf_azimuth_radians.len, &context.detailed_canopy.*.?);
            const replacement_carbon_exchange = ecosys.canopy_carbon_exchange.State.init(context.allocator, context.detailed_canopy.*.?.branch_node_offsets.len - 1) catch |err| {
                replacement_layers.deinit();
                return err;
            };
            var previous_layers = context.canopy_layer_distribution.*.?;
            context.canopy_layer_distribution.*.? = replacement_layers;
            previous_layers.deinit();
            var previous_carbon_exchange = context.canopy_carbon_exchange.*.?;
            context.canopy_carbon_exchange.*.? = replacement_carbon_exchange;
            previous_carbon_exchange.deinit();
        }
        _ = try ecosys.plant_topology.advanceRootAxes(.{
            .roots = &context.plant_roots.*.?,
            .canopy = &context.detailed_canopy.*.?,
            .growth_stages = &context.plant_growth_stages.*.?,
            .branch_development = &context.branch_development.*.?,
            .active_by_plant = context.plant_phenology.*.?.active,
            .root_branching_carbon_fraction = context.plant_topology_controls.root_branching_carbon_fraction,
            .minimum_root_turgor_potential_mpa = context.runscript.phenology_parameters.minimum_turgor_potential_mpa,
        });
    }
    if (context.plant_growth_stages.*) |*growth_stages| {
        const roots = if (context.plant_roots.*) |*value| value else return error.MissingPlantRootsForEmergence;
        const canopy = if (context.detailed_canopy.*) |*value| value else return error.MissingCanopyForEmergence;
        const water_workspace = if (context.plant_water_workspace.*) |*value| value else return error.MissingPlantWaterWorkspaceForEmergence;
        try ecosys.plant_development.refreshHypocotyledonHeight(canopy, water_workspace.seeding_depth_m);
        try ecosys.plant_development.refreshEmergence(canopy, roots, growth_stages, water_workspace.seeding_depth_m, plant_calendar.day_of_year, context.runscript.phenology_parameters.emergence_area_threshold_m2_per_plant, context.runscript.phenology_parameters.emergence_root_depth_margin_m, context.development_emerged);
    }
    if (context.plant_phenology.*) |*phenology| {
        const roots = if (context.plant_roots.*) |*value| value else return error.MissingPlantRootsForPhenology;
        const canopy = if (context.detailed_canopy.*) |*value| value else return error.MissingCanopyForPhenology;
        try ecosys.plant_root_gas_exchange.fillPlantOxygenUptakeToDemandFraction(roots, context.root_biological_domain_count_by_plant, context.phenology_root_oxygen_fraction);
        var phenology_context: ecosys.plant_phenology.AdvanceContext = .{
            .phenology = phenology,
            .plants = context.plants,
            .thermal_acclimation_offset_k = canopy.plant_thermal_adaptation_offset_c,
            .canopy_turgor_potential_mpa = canopy.plant_canopy_turgor_potential_mpa,
            .root_oxygen_uptake_to_demand_fraction = context.phenology_root_oxygen_fraction,
            .emerged = context.development_emerged,
            .timestep_hours = 1,
            .parameters = context.runscript.phenology_parameters,
        };
        try runKernelAcrossSerialTiles(context, &phenology_context, ecosys.plant_phenology.advanceTile);
    }
    if (context.plant_growth_stages.*) |*growth_stages| {
        const roots = if (context.plant_roots.*) |*value| value else return error.MissingPlantRootsForDevelopment;
        const canopy = if (context.detailed_canopy.*) |*value| value else return error.MissingCanopyForDevelopment;
        const dormancy_state = if (context.plant_dormancy.*) |*value| value else return error.MissingDormancyStateForDevelopment;
        try ecosys.plant_development.refreshCanopyHeight(canopy, context.development_canopy_height_m);
        try ecosys.plant_development.refreshSoilWaterPotentials(
            context.grid,
            context.soil_hourly_workspace.landscape_total_water_potential_mpa,
            context.surface_litter_water_environment.matric_plus_osmotic_water_potential_mpa,
            context.terrain_hydrology.relative_surface_elevation_m,
            context.soil_hourly_workspace.gravitational_water_potential_mpa_per_m,
            roots,
            context.development_surface_water_potential_mpa,
            context.development_seed_layer_water_potential_mpa,
        );
        var development_context: ecosys.plant_development.AdvanceContext = .{
            .growth_stages = growth_stages,
            .dormancy_state = dormancy_state,
            .phenology_state = &context.plant_phenology.*.?,
            .branch_development = &context.branch_development.*.?,
            .species_parameters_by_plant = context.development_species_parameters,
            .dormancy_parameters_by_plant = context.development_dormancy_parameters,
            .planting_day_of_year_by_plant = context.development_planting_day_of_year,
            .planting_year_by_plant = context.development_planting_year,
            .canopy_height_m_by_plant = context.development_canopy_height_m,
            .snow_depth_m_by_cell = context.snow_depth_m,
            .canopy_temperature_k_by_plant = context.plants.canopy_temperature_k,
            .canopy_turgor_potential_mpa_by_plant = canopy.plant_canopy_turgor_potential_mpa,
            .canopy_total_water_potential_mpa_by_plant = context.plants.canopy_water_potential_mpa,
            .surface_soil_water_potential_mpa_by_cell = context.development_surface_water_potential_mpa,
            .seed_layer_soil_water_potential_mpa_by_plant = context.development_seed_layer_water_potential_mpa,
            .emerged_by_plant = context.development_emerged,
            .calendar_by_cell = plant_calendar_by_cell,
            .timestep_h = 1,
        };
        try runKernelAcrossSerialTiles(context, &development_context, ecosys.plant_development.advanceTile);
        var reproduction_context: ecosys.plant_reproduction.ApplyContext = .{ .canopy = canopy, .plants = context.plants, .growth_stages = growth_stages, .controls = context.plant_reproduction_controls, .active_by_plant = context.plant_phenology.*.?.active, .minimum_turgor_potential_mpa = context.runscript.phenology_parameters.minimum_turgor_potential_mpa, .seed_set_parameters = context.runscript.seed_set_parameters, .structural_presence_threshold_g_per_plant = context.runscript.plant_pool_parameters.branch_structural_presence_g_per_plant, .timestep_h = 1 };
        try runKernelAcrossSerialTiles(context, &reproduction_context, ecosys.plant_reproduction.applyTile);
    }
    if (context.canopy_layer_distribution.*) |*layers| {
        const canopy = if (context.detailed_canopy.*) |*value| value else return error.MissingCanopyForLayerDistribution;
        const growth_stages = if (context.plant_growth_stages.*) |*value| value else return error.MissingGrowthStagesForLayerDistribution;
        const structure = if (context.canopy_structure.*) |*value| value else return error.MissingCanopyStructureForLayerDistribution;
        try layers.refresh(canopy, growth_stages, context.canopy_layer_controls, context.development_emerged, context.canopy_geometry.leaf_inclination_sine, structure.leaf_inclination_fraction, context.hourly_solar_angle_sine, 1.0e-12);
        try layers.publishNodeSamples(canopy);
        if (context.canopy_surface_input_workspace.*) |*surface_workspace| {
            var carboxylation_context: ecosys.canopy_carboxylation.ApplyContext = .{
                .canopy = canopy,
                .layers = layers,
                .interception = &context.canopy_interception.*.?,
                .optics = &context.canopy_optics.*.?,
                .geometry = context.canopy_geometry,
                .parameters_by_plant = context.canopy_biochemistry_parameters,
                .c4_carbon_parameters = context.runscript.c4_carbon_parameters,
                .direct_incidence_fraction = context.direct_incidence_fraction,
                .canopy_temperature_k_by_plant = context.plants.canopy_temperature_k,
                .stomatal_resistance_h_per_m_by_plant = surface_workspace.stomatal_resistance_h_per_m,
                .minimum_stomatal_resistance_h_per_m_by_plant = canopy.plant_minimum_water_vapor_resistance_h_per_m,
                .shallow_root_profile_by_plant = context.canopy_layer_controls.root_profile_type,
                .canopy_turgor_potential_mpa_by_plant = canopy.plant_canopy_turgor_potential_mpa,
                .minimum_turgor_potential_mpa = context.runscript.phenology_parameters.minimum_turgor_potential_mpa,
                .atmospheric_co2_umol_per_mol = context.current_atmospheric_co2_umol_per_mol.*,
                .picard_relaxation = context.config.picard_relaxation,
                .timestep_h = 1,
                .carbon_exchange = if (context.canopy_carbon_exchange.*) |*ledger| ledger else null,
            };
            try runKernelAcrossSerialTiles(context, &carboxylation_context, ecosys.canopy_carboxylation.applyTile);
        }
    }
    if (context.detailed_canopy.*) |*canopy| {
        const shoot_execution_year = std.math.cast(u16, plant_calendar.current_year) orelse
            return error.InvalidShootGrowthDate;
        const litter_partition = if (context.plant_litter_partition.*) |*value| value else return error.MissingPlantLitterPartitionForShootGrowth;
        const canopy_surface_exchange_for_growth = if (context.canopy_surface_exchange.*) |*value| value else return error.MissingCanopySurfaceExchangeForShootGrowth;
        const canopy_surface_workspace_for_growth = if (context.canopy_surface_input_workspace.*) |*value| value else return error.MissingCanopySurfaceWorkspaceForShootGrowth;
        const canopy_retention_for_growth = if (context.canopy_precipitation_retention.*) |*value| value else return error.MissingCanopyRetentionForShootGrowth;
        const surface_gas_for_growth = context.surface_gas_parameters.* orelse return error.MissingSurfaceGasParametersForShootGrowth;
        const shoot_solar_noon_hour_by_cell = try context.allocator.alloc(u8, context.grid.cell_count);
        defer context.allocator.free(shoot_solar_noon_hour_by_cell);
        for (weather_header_by_cell, shoot_solar_noon_hour_by_cell) |header, *solar_noon_hour| {
            if (!std.math.isFinite(header.solar_noon_hour) or header.solar_noon_hour < 0 or header.solar_noon_hour > 23)
                return error.InvalidShootSolarNoonHour;
            solar_noon_hour.* = @intFromFloat(@floor(header.solar_noon_hour));
        }
        @memset(context.shoot_senescence_products_by_plant, .{});
        var shoot_growth_context: ecosys.shoot_growth_runtime.ApplyContext = .{
            .canopy = canopy,
            .growth_stages = &context.plant_growth_stages.*.?,
            .dormancy = &context.plant_dormancy.*.?,
            .development = &context.branch_development.*.?,
            .plant_parameters = context.shoot_growth_plant_parameters,
            .active_by_plant = context.plant_phenology.*.?.active,
            .emerged_by_plant = context.development_emerged,
            .roots = &context.plant_roots.*.?,
            .soil_temperature_k = context.grid.soil_temperature_k,
            .soil_layer_capacity = context.grid.soil_layer_capacity,
            .root_growth_temperature_parameters = context.runscript.canopy_stress_parameters.growth_temperature,
            .canopy_temperature_k_by_plant = context.plants.canopy_temperature_k,
            .canopy_total_water_potential_mpa_by_plant = context.plants.canopy_water_potential_mpa,
            .total_aerodynamic_resistance_h_per_m_by_plant = canopy_surface_exchange_for_growth.total_aerodynamic_resistance_h_per_m,
            .stomatal_resistance_h_per_m_by_plant = canopy_surface_workspace_for_growth.stomatal_resistance_h_per_m,
            .plant_radiation_fraction = canopy_retention_for_growth.living_radiation_fraction,
            .atmospheric_ammonia_g_n_per_m3 = surface_gas_for_growth.atmospheric_concentration_g_per_m3[@intFromEnum(ecosys.gas_transport.Species.ammonia)],
            .ammonia_solubility_at_25_c = surface_gas_for_growth.solubility.reference_water_to_air[@intFromEnum(ecosys.gas_transport.Species.ammonia)],
            .canopy_ammonia_parameters = context.runscript.canopy_ammonia_exchange_parameters,
            .partition_parameters = context.runscript.organ_partition_parameters,
            .metabolism_parameters = context.runscript.shoot_metabolism_parameters,
            .phenology_parameters = context.runscript.phenology_parameters,
            .node_growth_parameters = context.runscript.shoot_node_growth_parameters,
            .branch_mobile_exchange_parameters = context.runscript.branch_mobile_exchange_parameters,
            .storage_remobilization_duration_h_by_growth_habit = context.runscript.storage_remobilization_parameters.remobilization_duration_h,
            .cell_area_m2 = context.canopy_cell_area_m2,
            .stalk_volume_m3_per_g_c = context.runscript.stalk_volume_m3_per_g_c,
            .structural_presence_threshold_g_per_plant = context.runscript.plant_pool_parameters.branch_structural_presence_g_per_plant,
            .grain_fill_detection_threshold_g_per_plant = context.runscript.plant_pool_parameters.grain_fill_detection_g_c_per_plant,
            .day_of_year = plant_calendar.day_of_year,
            .hour_of_day = hour_of_day,
            .solar_noon_hour_by_cell = shoot_solar_noon_hour_by_cell,
            .execution_year = shoot_execution_year,
            .timestep_h = 1,
            .dormancy_parameters_by_plant = context.development_dormancy_parameters,
            .litter_partition = litter_partition,
            .senescence_recycling = .{
                .minimum_carbon_fraction = context.runscript.root_metabolism_parameters.minimum_carbon_recycling_fraction,
                .responsive_carbon_fraction = context.runscript.root_metabolism_parameters.responsive_carbon_recycling_fraction,
                .maximum_nitrogen_fraction = context.runscript.root_metabolism_parameters.maximum_nitrogen_recycling_fraction,
                .maximum_phosphorus_fraction = context.runscript.root_metabolism_parameters.maximum_phosphorus_recycling_fraction,
            },
            .senescence_products_by_plant = context.shoot_senescence_products_by_plant,
            .senescence_demand_tolerance_g_c = context.config.absolute_tolerance,
            .leaf_area_presence_tolerance_m2 = context.config.absolute_tolerance,
            .symbiosis_parameters = context.runscript.symbiotic_fixation_parameters,
            .carbon_exchange = if (context.canopy_carbon_exchange.*) |*ledger| ledger else null,
        };
        try runKernelAcrossSerialTiles(context, &shoot_growth_context, ecosys.shoot_growth_runtime.applyTile);
        try applyNaturalBranchMortality(context, plant_calendar);
        for (0..canopy.cell_count) |cell| {
            var cell_products: ecosys.canopy_photosynthesis.SenescenceProducts = .{};
            const first_plant = cell * canopy.species_count;
            for (context.shoot_senescence_products_by_plant[first_plant .. first_plant + canopy.species_count]) |products|
                ecosys.canopy_photosynthesis.addSenescenceProducts(&cell_products, products);
            try ecosys.shoot_litter_bridge.commitCell(context.surface_organic, cell, cell_products);
        }
        const roots = if (context.plant_roots.*) |*value| value else return error.MissingPlantRootsForShootRootExchange;
        var root_mycorrhizal_exchange_context: ecosys.plant_root_mycorrhizal_exchange.ApplyContext = .{
            .roots = roots,
            .cell_count = context.grid.cell_count,
            .species_count = canopy.species_count,
            .active_soil_layer_count_by_cell = context.grid.active_soil_layer_count,
            .active_by_plant = context.plant_phenology.*.?.active,
            .plant_parameters = context.root_metabolism_plant_parameters,
            .parameters = context.runscript.root_mycorrhizal_exchange_parameters,
            .timestep_h = 1,
            .dynamic_salts = context.runscript.dynamic_plant_salts,
        };
        try runKernelAcrossSerialTiles(context, &root_mycorrhizal_exchange_context, ecosys.plant_root_mycorrhizal_exchange.applyTile);
        var shoot_root_exchange_context: ecosys.plant_shoot_root_exchange.ApplyContext = .{
            .canopy = canopy,
            .roots = roots,
            .growth_stages = &context.plant_growth_stages.*.?,
            .active_by_plant = context.plant_phenology.*.?.active,
            .biomass_turnover_type_by_plant = context.canopy_layer_controls.biomass_turnover_type,
            .plant_parameters = context.root_metabolism_plant_parameters,
            .parameters = context.runscript.shoot_root_exchange_parameters,
            .root_nonwoody_fraction_exponent = context.runscript.root_metabolism_parameters.nonwoody_root_fraction_exponent,
            .structural_presence_threshold_g_per_plant = context.runscript.plant_pool_parameters.branch_structural_presence_g_per_plant,
            .timestep_h = 1,
            .dynamic_salts = context.runscript.dynamic_plant_salts,
        };
        try runKernelAcrossSerialTiles(context, &shoot_root_exchange_context, ecosys.plant_shoot_root_exchange.applyTile);
    }
    // REDIST DORGC is the source-signed accepted loss from all soil organic C owners over
    // the complete hour, including microbial, residue, dissolved, adsorbed,
    // acetate, and structural pools.
    try ecosys.soil_organic_carbon_change.publishAcceptedHourlyChange(context.soil_organic, context.soil_organic_carbon_at_hour_start_g_c, context.soil_organic_carbon_change_g_c_per_h);
    // REDIST HEATIN precipitation term. Publish only after every hourly
    // process above has accepted its state. Rainfall already includes the
    // runtime irrigation addition; using the pre-routing atmospheric depths
    // also counts canopy-retained water without double-counting snow, litter,
    // or soil ingress.
    try context.landscape_boundary_ledger.accumulateAcceptedPrecipitationHeat(
        context.atmosphere.rainfall_m,
        context.atmosphere.snowfall_water_equivalent_m,
        context.canopy_cell_area_m2,
        context.atmosphere.air_temperature_k,
    );
    // REDIST XHFLF0/XHFLV0: the snow inventory reconstructs sensible
    // capacity*temperature only, so accepted phase latent energy changes are
    // boundary adjustments rather than additional snow storage.
    try context.landscape_boundary_ledger.accumulateAcceptedSignedHeat(
        snow_phase_change_report.sensible_energy_change_mj +
            snow_vapor_equilibrium_report.sensible_energy_change_mj,
    );
    const subsurface_irrigation_input =
        try ecosys.subsurface_irrigation_heat.calculate(
            context.subsurface_irrigation_water_m3,
            context.atmosphere.air_temperature_k,
            context.grid.soil_layer_capacity,
            context.runscript.canopy_surface_exchange_parameters.liquid_water_heat_capacity_mj_per_m3_k,
        );
    try context.landscape_boundary_ledger.accumulateAccepted(.{
        .rain_m3 = subsurface_irrigation_input.water_input_m3,
        .heat_input_mj = subsurface_irrigation_input.heat_input_mj,
    });
    const subsurface_irrigation_solutes =
        try ecosys.subsurface_irrigation_chemistry.boundaryInput(
            context.irrigation_loads,
            subsurface_irrigation_chemistry_parameters,
        );
    try context.landscape_boundary_ledger.accumulateAccepted(.{
        .nitrogen_input_g_n = subsurface_irrigation_solutes.nitrogen_g_n,
        .phosphorus_input_g_p = subsurface_irrigation_solutes.phosphorus_g_p,
        .ion_input_mol = subsurface_irrigation_solutes.ion_mol,
    });
}

fn sameCalendarDay(left: ecosys.weather.Timestamp, right: ecosys.weather.Timestamp) bool {
    return left.year == right.year and left.day_of_year == right.day_of_year and left.month == right.month and left.day_of_month == right.day_of_month;
}

fn reconstructLandscapeMassBalance(context: anytype) !ecosys.mass_balance_audit.Totals {
    return ecosys.landscape_mass_balance_runtime.reconstruct(.{
        .grid = context.grid,
        .plants = context.plants,
        .snow = context.snow_transport,
        .soil_thermal = context.soil_thermal,
        .soil_properties = context.soil_solver_properties,
        .soil_gas = context.gas_transport,
        .soil_organic = context.soil_organic,
        .surface_organic = context.surface_organic,
        .mineral_nitrogen = context.mineral_nitrogen_transport,
        .soil_chemistry = context.soil_chemistry,
        .nitrogen_fertilizer = context.soil_fertilizer_inventory,
        .mineral_fertilizer = context.mineral_fertilizer_inventory,
        .micropore_solutes = context.micropore_solute_state,
        .macropore_solutes = context.macropore_solute_state,
        .surface_chemistry = context.surface_litter_chemistry,
        .surface_fertilizer = context.surface_litter_fertilizer,
        .surface = context.surface_precipitation,
        .surface_ice_water_equivalent_m3 = context.surface_litter_ice_m3,
        .surface_gas = context.litter_gas_transport,
        .surface_litter_dry_mass_Mg = context.surface_litter_geometry.dry_mass_Mg,
        .canopy_retention = if (context.canopy_precipitation_retention.*) |*value| value else null,
        .cell_area_m2 = context.canopy_cell_area_m2,
        .soil_mass_Mg_scratch = context.landscape_soil_mass_Mg_scratch,
        .parameters = .{
            .snow_ice_density_Mg_per_m3 = context.runscript.snow_ice_density_Mg_per_m3,
            .carbon_g_per_mol = 12,
            .nitrogen_g_per_mol = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
            .phosphorus_g_per_mol = context.runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
            .phosphate_zone_fractions = .{
                .ammonium_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                .ammonium_band = context.runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                .nitrate_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                .nitrate_band = context.runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                .phosphate_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
                .phosphate_band = context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
            },
            .surface_physical = .{
                .dry_organic_heat_capacity_mj_per_g_c_k = context.runscript.surface_pond_dry_organic_heat_capacity_mj_per_g_c_k,
                .liquid_water_heat_capacity_mj_per_m3_k = context.runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_mj_per_m3_k,
                .ice_heat_capacity_mj_per_m3_k = context.runscript.soil_phase_heat_parameters.ice_heat_capacity_mj_per_m3_k,
                .water_molar_mass_g_per_mol = context.runscript.soil_gas_transport_parameters.water_molar_mass_g_per_mol,
                .liquid_water_density_g_per_m3 = context.runscript.soil_gas_transport_parameters.water_density_g_per_m3,
            },
        },
    }, context.landscape_boundary_ledger);
}

fn sameWeatherTimestamp(left: ecosys.weather.Timestamp, right: ecosys.weather.Timestamp) bool {
    return left.year == right.year and
        left.day_of_year == right.day_of_year and
        left.month == right.month and
        left.day_of_month == right.day_of_month and
        left.hour == right.hour and
        left.minute == right.minute;
}

fn dayOfYearFromTimestamp(timestamp: ecosys.weather.Timestamp) !u16 {
    if (timestamp.day_of_year) |day| return day;
    const date = try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp);
    return try (ecosys.plant_management.PackedDate{ .day = date.day, .month = date.month, .year = date.year }).dayOfYear(date.year);
}

fn soilOutputCatalog(allocator: std.mem.Allocator, editor_index: usize, layers: usize) !ecosys.soil_output_catalog.Catalog {
    return switch (editor_index) {
        0 => ecosys.soil_output_catalog.carbon(allocator, layers, layers, layers),
        1 => ecosys.soil_output_catalog.water(allocator, layers),
        2 => ecosys.soil_output_catalog.nitrogen(allocator, layers, layers),
        3 => ecosys.soil_output_catalog.phosphorus(allocator),
        4 => ecosys.soil_output_catalog.heat(allocator, layers),
        5 => ecosys.soil_output_catalog.dailyCarbon(allocator, layers),
        6 => ecosys.soil_output_catalog.dailyWater(allocator, layers, layers, layers),
        7 => ecosys.soil_output_catalog.dailyNitrogen(allocator, layers),
        8 => ecosys.soil_output_catalog.dailyPhosphorus(allocator, layers),
        9 => ecosys.soil_output_catalog.dailyHeat(allocator, layers, layers),
        else => error.OutputEditorIndexOutOfBounds,
    };
}

fn plantOutputCatalog(allocator: std.mem.Allocator, editor_index: usize, layers: usize) !ecosys.plant_output_catalog.Catalog {
    return switch (editor_index) {
        0 => ecosys.plant_output_catalog.carbon(allocator),
        1 => ecosys.plant_output_catalog.water(allocator, layers),
        2 => ecosys.plant_output_catalog.nitrogen(allocator, layers),
        3 => ecosys.plant_output_catalog.phosphorus(allocator, layers),
        4 => ecosys.plant_output_catalog.heat(allocator),
        5 => ecosys.plant_output_catalog.dailyCarbon(allocator, layers),
        6 => ecosys.plant_output_catalog.dailyWater(allocator),
        7 => ecosys.plant_output_catalog.dailyNitrogen(allocator),
        8 => ecosys.plant_output_catalog.dailyPhosphorus(allocator),
        9 => ecosys.plant_output_catalog.dailyDevelopment(allocator),
        else => error.OutputEditorIndexOutOfBounds,
    };
}

pub fn main(init: std.process.Init) !void {
    var debug_allocator: std.heap.DebugAllocator(.{}) = .init;
    defer _ = debug_allocator.deinit();
    const allocator = debug_allocator.allocator();

    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 2) {
        std.log.err("usage: ecosys_ng <runscript>", .{});
        return error.MissingRunscriptPath;
    }
    const source = try std.Io.Dir.cwd().readFileAlloc(init.io, args[1], allocator, .limited(16 * 1024 * 1024));
    defer allocator.free(source);
    var runscript = try ecosys.runscript.parse(allocator, source);
    defer runscript.deinit();
    if (runscript.geospatial_bounds == null) {
        std.log.err("missing compulsory geospatial_grid record: expected minimum/maximum latitude, minimum/maximum longitude, and latitude/longitude intervals", .{});
        return error.MissingGeospatialGridRecord;
    }
    if (runscript.uses_four_value_species_default) {
        std.log.warn("four-value domain header has no plant-species count; compatibility mode selected 5 species. Add a fifth value to make this explicit", .{});
    }
    if (runscript.uses_compatibility_runtime_controls) {
        std.log.warn("runscript has no runtime control record; compatibility controls selected. Add 'runtime,workers,tile_cells,relative_tolerance,absolute_tolerance,max_nonlinear_iterations,picard_relaxation' after the domain header", .{});
    } else {
        std.log.warn("runtime max_nonlinear_iterations is retained for input compatibility but process solvers use option-derived NPH/NPG or legacy NPR/NPS/NPRS iteration ceilings", .{});
    }
    if (runscript.uses_compatibility_soil_solver_parameters) {
        std.log.warn("runscript has no soil_solver record; historical HOUR1 coefficients selected at runtime. Add an explicit soil_solver record to control retention and hydraulic-conductivity science", .{});
    }
    if (runscript.uses_compatibility_soil_process_parameters) {
        std.log.warn("runscript has no soil_process record; historical gravity, vapor-diffusion, tortuosity, and osmotic-reflection coefficients selected at runtime", .{});
    }
    if (runscript.uses_compatibility_soil_gas_transport_parameters) {
        std.log.warn("runscript has no soil_gas_transport record; historical seven-gas free-air diffusivities, temperature response, Penman tortuosity, and gas-volume constants selected at runtime", .{});
    }
    if (runscript.uses_compatibility_surface_runoff) {
        std.log.warn("runscript has no surface_runoff record; HOUR1/WATSUB-compatible retention, roughness, hydraulic-volume, time-conversion, and negligible-water controls selected at runtime", .{});
    }
    if (runscript.uses_compatibility_rainfall_impact) {
        std.log.warn("runscript has no rainfall_impact record; HOUR1 direct-rain, throughfall, ponding attenuation, and conductivity-damage coefficients selected at runtime", .{});
    }
    if (runscript.uses_compatibility_soil_phase_heat_parameters) {
        std.log.warn("runscript has no soil_phase_heat record; historical vapor-equilibrium, freeze-thaw, and heat-convection coefficients selected at runtime", .{});
    }
    if (runscript.uses_compatibility_geothermal_controls) {
        std.log.warn("runscript has no geothermal record; STARTS-compatible source depth, conductivity, and 57 mW m-2 geothermal flux selected at runtime", .{});
    }
    if (runscript.uses_compatibility_subsurface_state_controls) {
        std.log.warn("runscript has no subsurface_state record; HOUR1 THETPW=0.001 and ZERO2=1e-6 selected at runtime for water-table and active-layer diagnostics", .{});
    }
    if (runscript.uses_compatibility_soil_geometry_parameters) {
        std.log.warn("runscript has no soil_geometry record; REDIST FORGC=110000 g C Mg-1, 1.82e-6 m3 g-1 SOC volume, DENSJ=0.083, and 1e-9 m minimum layer thickness selected at runtime", .{});
    }
    if (runscript.uses_compatibility_surface_energy) {
        std.log.warn("runscript has no surface-energy record; compatibility controls selected. Add 'surface_energy,soil_emissivity,snow_emissivity,canopy_emissivity,snow_full_cover_depth_m,sensible_conductance_mj_per_m2_h_k,latent_conductance_mj_per_m2_h_kpa,vapor_activity,minimum_temperature_k,maximum_temperature_k' before the site filename", .{});
    }
    if (runscript.uses_compatibility_canopy_ammonia_exchange) {
        std.log.warn("runscript has no canopy_ammonia_exchange record; historical UPTAKE FDMP and RNH3B coefficients selected at runtime", .{});
    }
    if (runscript.uses_compatibility_plant_structure) {
        std.log.warn("runscript has no plant_structure record; source-compatible 10 root axes selected at runtime. Add 'plant_structure,root_axes_per_plant' before the site filename to remove that compatibility default", .{});
    }
    if (runscript.uses_compatibility_canopy_layers) {
        std.log.warn("runscript has no canopy_layers record; source-compatible 10 canopy layers selected at runtime. Add 'canopy_layers,count' before the site filename", .{});
    }
    if (runscript.uses_compatibility_canopy_discretization) {
        std.log.warn("runscript has no canopy_discretization record; source-compatible 4 leaf-inclination, 4 leaf-azimuth, and 4 diffuse-sky classes selected at runtime. Add 'canopy_discretization,inclination_count,leaf_azimuth_count,diffuse_sky_count'", .{});
    }
    if (runscript.uses_compatibility_canopy_geometry) {
        std.log.warn("runscript has no canopy_geometry record; STARTQ stalk volume:carbon ratio 4.0e-6 m3 g-C-1 selected at runtime. Add 'canopy_geometry,stalk_volume_m3_per_g_c' before the site filename", .{});
    }
    if (runscript.uses_compatibility_standing_dead_partition) {
        std.log.warn("runscript has no standing_dead_partition record; STARTQ coarse-wood carbon fractions and CNOPC/CPOPC weights selected as runtime compatibility data", .{});
    }
    if (runscript.uses_compatibility_plant_heat_water) {
        std.log.warn("runscript has no plant_initial_heat_water record; STARTQ temperature, vapor-pressure, humidity, and initial water-potential coefficients selected as runtime compatibility data", .{});
    }
    if (runscript.uses_compatibility_plant_geometry) {
        std.log.warn("runscript has no plant_initial_geometry record; STARTQ seed and root geometry coefficients selected as runtime compatibility data", .{});
    }
    if (runscript.uses_compatibility_phenology_initialization) {
        std.log.warn("runscript has no plant_initial_phenology record; READQ perennial scaling and STARTQ concurrent-node thresholds selected as runtime compatibility data", .{});
    }
    if (runscript.uses_compatibility_root_initialization) {
        std.log.warn("runscript has no root_initialization record; STARTQ protein, mycorrhizal radius, water-potential, active-length, and water-fraction values selected as runtime compatibility data", .{});
    }
    if (runscript.uses_compatibility_root_morphology) {
        std.log.warn("runscript has no root_morphology record; GROSUB minimum secondary-root length and elastic modulus selected as runtime compatibility data", .{});
    }
    if (runscript.uses_compatibility_standing_dead_energy) {
        std.log.warn("runscript has no standing_dead_energy record; UPTAKE FARS, dry heat capacity, emissivity, and heat-capacity thresholds selected at runtime. Add 'standing_dead_energy,sapwood_thickness_m,dry_volume_heat_capacity_mj_per_m3_k,emissivity,activation_heat_capacity_mj_per_m2_k,effective_heat_capacity_floor_mj_per_m2_k'", .{});
    }
    if (runscript.uses_compatibility_woody_optics) {
        std.log.warn("runscript has no woody_optics record; HOUR1 stalk and standing-dead SW/PAR albedos of 0.1 selected at runtime. Add 'woody_optics,stalk_sw_albedo,stalk_par_albedo,dead_sw_albedo,dead_par_albedo' before the site filename", .{});
    }
    if (runscript.uses_compatibility_canopy_retention) {
        std.log.warn("runscript has no canopy_retention record; historical HOUR1 XVOLWC capacities, low-sun extinction, and solar-share threshold selected at runtime", .{});
    }
    if (runscript.uses_compatibility_shoot_control_parameters) {
        std.log.warn("runscript has no shoot_controls record; historical STARTQ time conversion, cuticular CO2 ratio, and C3/C4 intercellular O2 values selected at runtime", .{});
    }
    if (runscript.uses_source_c4_carbon_parameters) {
        std.log.warn("runscript has no c4_carbon record; historical GROSUB bundle-sheath transfer, decarboxylation, and leakage coefficients selected at runtime", .{});
    }
    if (runscript.uses_compatibility_thermal_acclimation_parameters) {
        std.log.warn("runscript has no thermal_controls record; historical STARTQ acclimation, leafout/leafoff, and C3/C4 seed-set heat coefficients selected at runtime", .{});
    }
    if (runscript.uses_compatibility_canopy_stress_parameters) {
        std.log.warn("runscript has no canopy_stress record; historical UPTAKE CHILL/HEAT limits selected at runtime", .{});
    }
    if (runscript.uses_compatibility_phenology_parameters) {
        std.log.warn("runscript has no phenology_controls record; historical HFUNC Arrhenius, turgor, and oxygen-stress coefficients selected at runtime", .{});
    }
    if (runscript.uses_compatibility_plant_pool_parameters) {
        std.log.warn("runscript has no plant_pool_controls record; historical HFUNC ZERO/ZERO2/CNKIC/CPKIC values and static plant salts selected at runtime", .{});
    }
    if (runscript.uses_compatibility_seed_set_parameters) {
        std.log.warn("runscript has no seed_set_controls record; historical GROSUB SETC/SETN/SETP values selected at runtime", .{});
    }
    if (runscript.uses_compatibility_root_gas_parameters) {
        std.log.warn("runscript has no root_gas record; historical HOUR1/WATSUB/UPTAKE O2 coefficients selected at runtime. Add 'root_gas,reference_temperature_k,oxygen_diffusivity_m2_per_h,temperature_exponent,oxygen_solubility,activity_coefficient,temperature_intercept,temperature_coefficient_per_c,pure_water_solute_mol_per_m3,oxygen_to_carbon_ratio,minimum_water_film_m,water_film_scale_m,water_film_log_intercept,water_film_potential_exponent' before the site filename", .{});
    }
    if (runscript.uses_compatibility_root_nutrient_parameters) {
        std.log.warn("runscript has no root_nutrients record; historical HOUR1/UPTAKE diffusivity, tortuosity, C/N/P feedback, minimum population uptake fraction, and phosphorus molar mass selected at runtime. Add 'root_nutrients,reference_temperature_k,ammonium_diffusivity_m2_per_h,nitrate_diffusivity_m2_per_h,phosphate_diffusivity_m2_per_h,temperature_exponent,liquid_tortuosity_coefficient,n_feedback_n,n_feedback_p,p_feedback_p,p_feedback_n,minimum_population_fraction_multiplier,phosphorus_molar_mass_g_per_mol' before the plant_nutrients or site record", .{});
    }
    if (runscript.dynamic_plant_salts and runscript.uses_compatibility_root_salt_parameters) {
        std.log.warn("runscript enables dynamic plant salts but has no root_salts record; historical UPTAKE diffusivities and root-concentration inhibition constants selected at runtime. Add the case-insensitive root_salts record before root_metabolism, plant_nutrients, or site", .{});
    }
    if (runscript.uses_compatibility_root_exudation_parameters) {
        std.log.warn("runscript has no root_exudation record; historical UPTAKE FEXU, mobile N/P fractions, and maximum root-C concentration selected at runtime. Add 'root_exudation,maximum_root_carbon_g_c_per_m3,nitrogen_exchange_fraction,phosphorus_exchange_fraction,exchange_rate_per_h' before root_metabolism, plant_nutrients, or site", .{});
    }
    if (runscript.uses_compatibility_root_mycorrhizal_exchange_parameters) {
        std.log.warn("runscript has no root_mycorrhizal_exchange record; source GROSUB partner-water floor and exchange rate selected at runtime", .{});
    }
    if (runscript.uses_compatibility_root_porosity_parameters) {
        std.log.warn("runscript has no root_porosity record; historical GROSUB PORT maximum, oxygen-stress induction, and relaxation rates selected at runtime. Add 'root_porosity,maximum_fraction,oxygen_stress_induction_per_h,relaxation_per_h' before root_metabolism, plant_nutrients, or site", .{});
    }
    if (runscript.uses_compatibility_root_metabolism_parameters) {
        std.log.warn("runscript has no root_metabolism record; historical GROSUB root respiration, recycling, protein, and nutrient-uptake respiration constants selected at runtime. Supply the complete case-insensitive root_metabolism record before storage_remobilization, plant_nutrients, or site", .{});
    }
    if (runscript.uses_compatibility_organ_partition_parameters) {
        std.log.warn("runscript has no organ_partition record; historical GROSUB stage partition, turnover response, and reserve-redirection coefficients selected at runtime", .{});
    }
    if (runscript.uses_compatibility_shoot_metabolism_parameters) {
        std.log.warn("runscript has no shoot_metabolism record; historical GROSUB VMXC, CCKM, RMPLT, ZPLFM, CNKI, CPKI, and nitrogen-assimilation respiration coefficients selected at runtime", .{});
    }
    if (runscript.uses_source_shoot_node_growth_parameters) {
        std.log.warn("runscript has no shoot_node_growth record; historical GROSUB CNWS, CPWS, SLA2, SSL2, SNL2, SLAX, SSLX, SNLX, leaf nutrient exchange, branch reserve exchange, ZPGRM, SETN, SETP, and FLG4X values selected at runtime", .{});
    }
    if (runscript.uses_source_branch_mobile_exchange_parameters) {
        std.log.warn("runscript has no branch_mobile_exchange record; historical GROSUB interbranch mobile C/N/P exchange fractions selected as runtime compatibility data", .{});
    }
    if (runscript.uses_source_symbiotic_fixation_parameters) {
        std.log.warn("runscript has no symbiotic_fixation record; historical GROSUB WTNDI, VMXO, RMPLT, EN2F, CZKM, CPKM, ZCKI, ZPKI, SPNDL, recycling, FXRN, and CNDLI values selected at runtime", .{});
    }
    if (runscript.uses_source_plant_fire_combustion_parameters) {
        std.log.warn("runscript has no plant_fire_combustion record; historical GROSUB TCMBX, TFNCX, Arrhenius, SPCMB, and charcoal values selected at runtime", .{});
    }
    if (runscript.uses_source_soil_fire_combustion_parameters) {
        std.log.warn("runscript has no soil_fire_combustion record; historical NITRO soil-pool and charcoal combustion values selected at runtime", .{});
    }
    if (runscript.uses_compatibility_shoot_root_exchange_parameters) {
        std.log.warn("runscript has no shoot_root_exchange record; source GROSUB partner floor, annual minimum rate, leaf-partition exponent, and salt exchange rate selected at runtime", .{});
    }
    if (runscript.uses_compatibility_storage_remobilization_parameters) {
        std.log.warn("runscript has no storage_remobilization record; historical GROSUB germination/leafout storage durations, oxidation, shoot/root partition, and nutrient-equilibration rates selected at runtime", .{});
    }
    if (runscript.uses_compatibility_plant_nutrient_initialization) {
        std.log.warn("runscript has no plant_nutrients record; STARTS-compatible zero band fractions and all-H2PO4 soluble phosphate selected at runtime. Add 'plant_nutrients,ammonium_band_fraction,nitrate_band_fraction,phosphate_band_fraction,h2po4_fraction' before the site filename", .{});
    }
    if (runscript.uses_compatibility_microbial_dimensions) {
        std.log.warn("runscript has no microbial_dimensions record; NITRO-compatible 6 substrates and 7 microbial populations selected at runtime. Add 'microbial_dimensions,substrate_count,population_count' before the site filename", .{});
    }
    if (runscript.organic_initialization_file == null) {
        std.log.warn("runscript has no organic_initialization_file record; source STARTS/NITRO organic allocation and decomposition parameters selected as replaceable runtime defaults", .{});
    }
    if (runscript.surface_gas_parameter_file == null) {
        std.log.warn("runscript has no surface_gas_parameter_file record; source HOUR1/NITRO atmosphere, litter gas, and surface microbial coefficients selected as replaceable runtime defaults", .{});
    }
    if (runscript.soil_nitrogen_parameter_file == null) {
        std.log.warn("runscript has no soil_nitrogen_parameter_file record; source NITRO coefficients selected as replaceable runtime defaults", .{});
    }
    if (runscript.chemistry_initialization == null) {
        std.log.warn("runscript has no chemistry_initialization record; source STARTE water and phosphate dissociation constants selected as replaceable runtime defaults", .{});
    }
    if (runscript.chemistry_initialization != null and runscript.chemistry_primary_initialization == null) {
        std.log.warn("runscript has no chemistry_units record; only STARTE H/OH/phosphate protonation is seeded. Add runtime molar masses, concentration multipliers, and minimum extract concentrations to seed all primary ions and solids", .{});
    }
    if (runscript.chemistry_primary_initialization != null and runscript.chemistry_reaction_file == null) {
        std.log.warn("runscript has no chemistry_reaction_file record; primary chemistry is seeded but the coupled STARTE Newton-Raphson/Picard equilibrium solve cannot run", .{});
    }
    if (runscript.uses_compatibility_fertilizer_units) {
        std.log.warn("runscript has no fertilizer_units or chemistry_units record; legacy fertilizer conversion selected 14 g N mol-1 at runtime. Add 'fertilizer_units,nitrogen_molar_mass_g_per_mol' before the site filename", .{});
    }

    const runscript_directory = std.fs.path.dirname(args[1]) orelse ".";
    var organic_parameters: ?ecosys.soil_organic_parameters.OwnedParameters = null;
    defer if (organic_parameters) |*parameters| parameters.deinit();
    if (runscript.organic_initialization_file) |parameter_name| {
        const parameter_path = try resolveInputPath(allocator, init.io, runscript_directory, parameter_name);
        defer allocator.free(parameter_path);
        const parameter_source = try std.Io.Dir.cwd().readFileAlloc(init.io, parameter_path, allocator, .limited(16 * 1024 * 1024));
        defer allocator.free(parameter_source);
        organic_parameters = try ecosys.soil_organic_parameters.parse(allocator, parameter_source);
    } else {
        organic_parameters = try ecosys.soil_organic_parameters.sourceParameters(allocator);
    }
    var surface_gas_parameters: ?ecosys.surface_gas_parameters.Parameters = null;
    if (runscript.surface_gas_parameter_file) |parameter_name| {
        const parameter_path = try resolveInputPath(allocator, init.io, runscript_directory, parameter_name);
        defer allocator.free(parameter_path);
        const parameter_source = try std.Io.Dir.cwd().readFileAlloc(init.io, parameter_path, allocator, .limited(1024 * 1024));
        defer allocator.free(parameter_source);
        surface_gas_parameters = try ecosys.surface_gas_parameters.parse(parameter_source);
    } else surface_gas_parameters = try ecosys.surface_gas_parameters.sourceParameters();
    if (surface_gas_parameters != null) std.log.info("loaded runtime surface gas and microbial oxygen parameters", .{});
    var soil_nitrogen_parameters: ?ecosys.soil_nitrogen_parameters.Parameters = null;
    if (runscript.soil_nitrogen_parameter_file) |parameter_name| {
        const parameter_path = try resolveInputPath(allocator, init.io, runscript_directory, parameter_name);
        defer allocator.free(parameter_path);
        const parameter_source = try std.Io.Dir.cwd().readFileAlloc(init.io, parameter_path, allocator, .limited(1024 * 1024));
        defer allocator.free(parameter_source);
        soil_nitrogen_parameters = try ecosys.soil_nitrogen_parameters.parse(parameter_source);
    } else soil_nitrogen_parameters = try ecosys.soil_nitrogen_parameters.sourceParameters();
    if (soil_nitrogen_parameters != null) std.log.info("loaded runtime layer-resolved soil nitrogen parameters", .{});
    var chemistry_reaction_parameters: ?ecosys.soil_chemistry_parameters.Parameters = null;
    if (runscript.chemistry_reaction_file) |parameter_name| {
        const parameter_path = try resolveInputPath(allocator, init.io, runscript_directory, parameter_name);
        defer allocator.free(parameter_path);
        const parameter_source = try std.Io.Dir.cwd().readFileAlloc(init.io, parameter_path, allocator, .limited(4 * 1024 * 1024));
        defer allocator.free(parameter_source);
        chemistry_reaction_parameters = try ecosys.soil_chemistry_parameters.parse(parameter_source);
    }
    const grid_input_path = try resolveInputPath(
        allocator,
        init.io,
        runscript_directory,
        runscript.grid_input_file,
    );
    defer allocator.free(grid_input_path);
    const grid_input_source = try std.Io.Dir.cwd().readFileAlloc(
        init.io,
        grid_input_path,
        allocator,
        .limited(16 * 1024 * 1024),
    );
    defer allocator.free(grid_input_source);
    var grid_input_files = try ecosys.grid_input_files.CellFiles.parse(
        allocator,
        grid_input_source,
        try runscript.domain.columns(),
        try runscript.domain.rows(),
    );
    defer grid_input_files.deinit();

    var site_catalog = ecosys.site_catalog.Catalog.init(allocator);
    defer site_catalog.deinit();
    var topography_catalog = ecosys.topography_catalog.Catalog.init(allocator);
    defer topography_catalog.deinit();
    for (grid_input_files.site_file_by_cell) |site_name| {
        if (site_catalog.find(site_name) != null) continue;
        const site_path = try resolveInputPath(allocator, init.io, runscript_directory, site_name);
        defer allocator.free(site_path);
        const site_source = try std.Io.Dir.cwd().readFileAlloc(
            init.io,
            site_path,
            allocator,
            .limited(16 * 1024 * 1024),
        );
        defer allocator.free(site_source);
        _ = site_catalog.appendFromSource(
            site_name,
            site_source,
            1,
            1,
        ) catch |err| {
            std.log.err("site input validation failed: file='{s}' error={s}", .{ site_name, @errorName(err) });
            return err;
        };
    }
    for (grid_input_files.topography_file_by_cell) |topography_name| {
        if (topography_catalog.find(topography_name) != null) continue;
        const topography_path = try resolveInputPath(allocator, init.io, runscript_directory, topography_name);
        defer allocator.free(topography_path);
        const topography_source = try std.Io.Dir.cwd().readFileAlloc(
            init.io,
            topography_path,
            allocator,
            .limited(64 * 1024 * 1024),
        );
        defer allocator.free(topography_source);
        _ = topography_catalog.appendFromSource(topography_name, topography_source) catch |err| {
            std.log.err("topography input validation failed: file='{s}' error={s}", .{ topography_name, @errorName(err) });
            return err;
        };
    }
    var grid_environment = try ecosys.grid_environment.Assignments.init(
        allocator,
        runscript.domain,
        grid_input_files,
        site_catalog,
        topography_catalog,
    );
    defer grid_environment.deinit();
    const lateral_connection_mode_by_cell = try allocator.alloc(
        u8,
        grid_environment.site_catalog_index_by_cell.len,
    );
    defer allocator.free(lateral_connection_mode_by_cell);
    for (lateral_connection_mode_by_cell, 0..) |*mode, cell| {
        mode.* = @intFromEnum(site_catalog.entries.items[
            grid_environment.site_catalog_index_by_cell[cell]
        ].site.lateral_connection_mode);
    }
    var geospatial_grid: ?ecosys.spatial_grid.RegularGrid = if (runscript.geospatial_bounds) |bounds|
        try ecosys.spatial_grid.RegularGrid.init(allocator, bounds)
    else
        null;
    defer if (geospatial_grid) |*grid| grid.deinit();
    var tile_plan: ?ecosys.spatial_grid.TilePlan = null;
    defer if (tile_plan) |*plan| plan.deinit();
    if (geospatial_grid) |grid| {
        const cell_count = try grid.cellCount();
        if (cell_count != grid_environment.horizontal_cell_width_m.len)
            return error.GeospatialGridEnvironmentDimensionMismatch;
        for (0..cell_count) |cell| {
            const selected_site = site_catalog.entries.items[
                grid_environment.site_catalog_index_by_cell[cell]
            ].site;
            try grid.validateSiteCoordinate(
                cell,
                selected_site.latitude_degrees_north,
                selected_site.longitude_degrees_east,
            );
            // Geospatial dimensions are authoritative. Site widths remain
            // available only to exterior-boundary fallbacks.
            grid_environment.horizontal_cell_width_m[cell] =
                grid.east_west_cell_width_m[cell];
            grid_environment.vertical_cell_width_m[cell] =
                grid.north_south_cell_width_m[cell];
        }
        try ecosys.spatial_grid.validateNeighborExchangeControls(
            grid.row_count,
            grid.column_count,
            lateral_connection_mode_by_cell,
        );
        tile_plan = try ecosys.spatial_grid.TilePlan.init(
            allocator,
            grid.row_count,
            grid.column_count,
            runscript.tile_row_count,
            runscript.tile_column_count,
            runscript.lateral_flow_halo_cell_count,
        );
    }
    if (tile_plan == null) {
        tile_plan = try ecosys.spatial_grid.TilePlan.init(
            allocator,
            try runscript.domain.rows(),
            try runscript.domain.columns(),
            runscript.tile_row_count,
            runscript.tile_column_count,
            runscript.lateral_flow_halo_cell_count,
        );
    }
    const site = site_catalog.entries.items[grid_environment.site_catalog_index_by_cell[0]].site;
    // Atmospheric composition is environmental state owned by the site
    // input, never a process-parameter default. Supplying the compulsory
    // surface-gas coefficient file must not disable this conversion.
    surface_gas_parameters.?.atmospheric_concentration_g_per_m3 = try ecosys.surface_gas_parameters.atmosphericConcentrationsFromUmolPerMol(.{
        site.atmospheric_co2_umol_mol,
        site.atmospheric_methane_umol_mol,
        site.atmospheric_oxygen_umol_mol,
        site.atmospheric_nitrogen_umol_mol,
        site.atmospheric_nitrous_oxide_umol_mol,
        site.atmospheric_ammonia_umol_mol,
        0,
    }, site.mean_annual_air_temperature_c + 273.15);

    var topography = try grid_environment.buildDomainTopography(
        runscript.domain,
        grid_input_files,
        topography_catalog,
    );
    defer topography.deinit();
    const unit_by_cell = try topography.buildCellUnitMap(runscript.domain, allocator);
    defer allocator.free(unit_by_cell);
    var terrain_hydrology_state = try ecosys.terrain_hydrology.State.initMapped(
        allocator,
        topography,
        unit_by_cell,
        grid_environment.horizontal_cell_width_m,
        grid_environment.vertical_cell_width_m,
        try runscript.domain.columns(),
        try runscript.domain.rows(),
    );
    defer terrain_hydrology_state.deinit();
    var surface_runoff_state = try ecosys.surface_runoff.State.init(allocator, unit_by_cell.len);
    defer surface_runoff_state.deinit();
    var surface_pond_transition_state = try ecosys.surface_pond_transition_step.State.init(allocator, unit_by_cell.len);
    defer surface_pond_transition_state.deinit();
    const surface_inorganic_nitrogen_export_g_n_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_inorganic_nitrogen_export_g_n_per_h);
    @memset(surface_inorganic_nitrogen_export_g_n_per_h, 0);
    const surface_inorganic_phosphorus_export_g_p_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_inorganic_phosphorus_export_g_p_per_h);
    @memset(surface_inorganic_phosphorus_export_g_p_per_h, 0);
    const surface_organic_carbon_export_g_c_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_organic_carbon_export_g_c_per_h);
    @memset(surface_organic_carbon_export_g_c_per_h, 0);
    const surface_inorganic_carbon_export_g_c_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_inorganic_carbon_export_g_c_per_h);
    @memset(surface_inorganic_carbon_export_g_c_per_h, 0);
    const surface_organic_nitrogen_export_g_n_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_organic_nitrogen_export_g_n_per_h);
    @memset(surface_organic_nitrogen_export_g_n_per_h, 0);
    const surface_organic_phosphorus_export_g_p_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_organic_phosphorus_export_g_p_per_h);
    @memset(surface_organic_phosphorus_export_g_p_per_h, 0);
    var surface_erosion_state = try ecosys.soil_erosion.RuntimeState.init(allocator, try runscript.domain.columns(), try runscript.domain.rows());
    defer surface_erosion_state.deinit();
    const surface_soil_mass_at_erosion_start_Mg = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(surface_soil_mass_at_erosion_start_Mg);
    const net_sediment_Mg_per_h = try allocator.alloc(f64, unit_by_cell.len);
    defer allocator.free(net_sediment_Mg_per_h);
    @memset(surface_soil_mass_at_erosion_start_Mg, 0);
    @memset(net_sediment_Mg_per_h, 0);
    var eroded_mineral_state = try ecosys.soil_erosion_mineral_bridge.State.init(allocator, unit_by_cell.len);
    defer eroded_mineral_state.deinit();
    var soil_catalog = ecosys.soil_catalog.Catalog.init(allocator);
    defer soil_catalog.deinit();
    const catalog_index_by_cell = try allocator.alloc(usize, unit_by_cell.len);
    defer allocator.free(catalog_index_by_cell);
    for (unit_by_cell, 0..) |unit_index, cell_index| {
        const profile_name = topography.units[unit_index].soil_profile_file;
        const catalog_index = soil_catalog.find(profile_name) orelse load: {
            const profile_path = try resolveInputPath(allocator, init.io, runscript_directory, profile_name);
            defer allocator.free(profile_path);
            const profile_source = try std.Io.Dir.cwd().readFileAlloc(init.io, profile_path, allocator, .limited(64 * 1024 * 1024));
            defer allocator.free(profile_source);
            break :load try soil_catalog.appendFromSource(profile_name, profile_source, runscript.soil_solver_parameters.retention, runscript.soil_solver_parameters.profile_derivation);
        };
        catalog_index_by_cell[cell_index] = catalog_index;
    }

    const scene_options = try allocator.alloc(ecosys.options.SceneOptions, runscript.scenes.len);
    defer allocator.free(scene_options);
    const weather_forcing_geometry_by_cell = try allocator.alloc(
        ecosys.weather_grid.ForcingGeometry,
        grid_environment.site_catalog_index_by_cell.len,
    );
    defer allocator.free(weather_forcing_geometry_by_cell);
    for (weather_forcing_geometry_by_cell, 0..) |*geometry, cell| {
        const selected_site = site_catalog.entries.items[
            grid_environment.site_catalog_index_by_cell[cell]
        ].site;
        geometry.* = .{
            .altitude_m = selected_site.elevation_m,
            .latitude_degrees_north = selected_site.latitude_degrees_north,
            .longitude_degrees_east = selected_site.longitude_degrees_east,
            .phytotron = selected_site.ecosystem_type == -2,
        };
    }
    const scene_weather_assignments = try allocator.alloc(
        ecosys.weather_grid.Assignments,
        runscript.scenes.len,
    );
    defer allocator.free(scene_weather_assignments);
    var initialized_scene_weather_assignments: usize = 0;
    defer for (scene_weather_assignments[0..initialized_scene_weather_assignments]) |*assignments|
        assignments.deinit();

    var total_weather_records: usize = 0;
    for (runscript.scenes, 0..) |scene, scene_index| {
        const options_path = try resolveInputPath(allocator, init.io, runscript_directory, scene.options);
        defer allocator.free(options_path);
        const options_source = try std.Io.Dir.cwd().readFileAlloc(init.io, options_path, allocator, .limited(16 * 1024 * 1024));
        defer allocator.free(options_source);
        scene_options[scene_index] = ecosys.options.parse(options_source) catch |err| {
            std.log.err("scene {d} options validation failed for '{s}'", .{ scene_index + 1, scene.options });
            return err;
        };
        _ = try ecosys.iteration_control.Limits.fromSceneOptions(scene_options[scene_index]);
        const weather_grid_path = try resolveInputPath(
            allocator,
            init.io,
            runscript_directory,
            scene.weather_grid_file,
        );
        defer allocator.free(weather_grid_path);
        const weather_grid_source = try std.Io.Dir.cwd().readFileAlloc(
            init.io,
            weather_grid_path,
            allocator,
            .limited(16 * 1024 * 1024),
        );
        defer allocator.free(weather_grid_source);
        var weather_files = ecosys.grid_input_files.WeatherFiles.parse(
            allocator,
            weather_grid_source,
            try runscript.domain.columns(),
            try runscript.domain.rows(),
        ) catch |err| {
            std.log.err("scene {d} weather-grid validation failed for '{s}': {s}", .{
                scene_index + 1,
                scene.weather_grid_file,
                @errorName(err),
            });
            return err;
        };
        defer weather_files.deinit();
        scene_weather_assignments[scene_index] = try ecosys.weather_grid.Assignments.init(
            allocator,
            weather_files,
            weather_forcing_geometry_by_cell,
        );
        initialized_scene_weather_assignments += 1;
        for (scene_weather_assignments[scene_index].unique_file_names) |weather_name| {
            const weather_path = try resolveInputPath(allocator, init.io, runscript_directory, weather_name);
            defer allocator.free(weather_path);
            var weather_summary = ecosys.weather.scanFile(allocator, init.io, weather_path, 64 * 1024) catch |err| {
                std.log.err("scene {d} weather validation failed for '{s}'", .{ scene_index + 1, weather_name });
                return err;
            };
            defer weather_summary.deinit();
            total_weather_records = try std.math.add(
                usize,
                total_weather_records,
                weather_summary.observation_count,
            );
        }
    }
    const iteration_limits = try ecosys.iteration_control.Limits.fromSceneOptions(scene_options[0]);
    const timeline = try ecosys.simulation_timeline.summarize(scene_options, runscript.scenarios, runscript.execution_repeat_count);
    var plant_assignments: ?ecosys.plant_assignment.Assignments = null;
    defer if (plant_assignments) |*assignments| assignments.deinit();
    var plant_unit_by_cell: ?[]usize = null;
    defer if (plant_unit_by_cell) |map| allocator.free(map);
    var plant_catalog = ecosys.plant_catalog.Catalog.init(allocator);
    defer plant_catalog.deinit();
    var plant_management_catalog = ecosys.plant_management.Catalog.init(allocator);
    defer plant_management_catalog.deinit();
    var land_management_assignments: ?ecosys.land_management.Assignments = null;
    defer if (land_management_assignments) |*assignments| assignments.deinit();
    var land_management_unit_by_cell: ?[]usize = null;
    defer if (land_management_unit_by_cell) |map| allocator.free(map);
    var disturbance_catalog = ecosys.disturbance_schedule.Catalog.init(allocator);
    defer disturbance_catalog.deinit();
    const disturbance_schedule_maps = try allocator.alloc(?ecosys.disturbance_management_dispatch.ScheduleMap, runscript.scenes.len);
    @memset(disturbance_schedule_maps, null);
    defer {
        for (disturbance_schedule_maps) |*maybe_map| if (maybe_map.*) |*map| map.deinit();
        allocator.free(disturbance_schedule_maps);
    }
    var fertilizer_catalog = ecosys.fertilizer_schedule.Catalog.init(allocator);
    defer fertilizer_catalog.deinit();
    const fertilizer_schedule_maps = try allocator.alloc(?ecosys.fertilizer_management_dispatch.ScheduleMap, runscript.scenes.len);
    @memset(fertilizer_schedule_maps, null);
    defer {
        for (fertilizer_schedule_maps) |*maybe_map| if (maybe_map.*) |*map| map.deinit();
        allocator.free(fertilizer_schedule_maps);
    }
    var irrigation_catalog = ecosys.irrigation_schedule.Catalog.init(allocator);
    defer irrigation_catalog.deinit();
    const irrigation_schedule_maps = try allocator.alloc(?ecosys.irrigation_management_dispatch.ScheduleMap, runscript.scenes.len);
    @memset(irrigation_schedule_maps, null);
    defer {
        for (irrigation_schedule_maps) |*maybe_map| if (maybe_map.*) |*map| map.deinit();
        allocator.free(irrigation_schedule_maps);
    }
    var output_selection_catalog = ecosys.output_selection.Catalog.init(allocator);
    defer output_selection_catalog.deinit();
    for (runscript.scenes, 0..) |scene, scene_index| for (scene.output_editors) |editor_name| {
        if (ecosys.delimited_input.isNo(editor_name) or output_selection_catalog.find(editor_name) != null) continue;
        const path = try resolveInputPath(allocator, init.io, runscript_directory, editor_name);
        defer allocator.free(path);
        const editor_source = try std.Io.Dir.cwd().readFileAlloc(init.io, path, allocator, .limited(16 * 1024 * 1024));
        defer allocator.free(editor_source);
        _ = output_selection_catalog.appendFromSource(editor_name, editor_source) catch |err| {
            std.log.err("scene {d} output-selection validation failed for '{s}'", .{ scene_index + 1, editor_name });
            return err;
        };
    };
    if (!ecosys.delimited_input.isNo(runscript.scenes[0].land_management)) {
        const management_path = try resolveInputPath(allocator, init.io, runscript_directory, runscript.scenes[0].land_management);
        defer allocator.free(management_path);
        const management_source = try std.Io.Dir.cwd().readFileAlloc(init.io, management_path, allocator, .limited(16 * 1024 * 1024));
        defer allocator.free(management_source);
        land_management_assignments = try ecosys.land_management.parse(allocator, management_source);
        land_management_unit_by_cell = try land_management_assignments.?.buildCellUnitMap(allocator, runscript.domain);
        for (land_management_assignments.?.units) |unit| {
            if (!ecosys.delimited_input.isNo(unit.tillage_file) and disturbance_catalog.find(unit.tillage_file) == null) {
                const path = try resolveInputPath(allocator, init.io, runscript_directory, unit.tillage_file);
                defer allocator.free(path);
                const schedule_source = try std.Io.Dir.cwd().readFileAlloc(init.io, path, allocator, .limited(64 * 1024 * 1024));
                defer allocator.free(schedule_source);
                _ = try disturbance_catalog.appendFromSource(unit.tillage_file, schedule_source);
            }
            if (!ecosys.delimited_input.isNo(unit.fertilizer_file) and fertilizer_catalog.find(unit.fertilizer_file) == null) {
                const path = try resolveInputPath(allocator, init.io, runscript_directory, unit.fertilizer_file);
                defer allocator.free(path);
                const schedule_source = try std.Io.Dir.cwd().readFileAlloc(init.io, path, allocator, .limited(64 * 1024 * 1024));
                defer allocator.free(schedule_source);
                _ = try fertilizer_catalog.appendFromSource(unit.fertilizer_file, schedule_source);
            }
            if (!ecosys.delimited_input.isNo(unit.irrigation_file) and irrigation_catalog.find(unit.irrigation_file) == null) {
                const path = try resolveInputPath(allocator, init.io, runscript_directory, unit.irrigation_file);
                defer allocator.free(path);
                const schedule_source = try std.Io.Dir.cwd().readFileAlloc(init.io, path, allocator, .limited(64 * 1024 * 1024));
                defer allocator.free(schedule_source);
                _ = try irrigation_catalog.appendFromSource(unit.irrigation_file, schedule_source);
            }
        }
        fertilizer_schedule_maps[0] = try ecosys.fertilizer_management_dispatch.ScheduleMap.init(allocator, land_management_assignments.?, land_management_unit_by_cell.?, fertilizer_catalog);
        disturbance_schedule_maps[0] = try ecosys.disturbance_management_dispatch.ScheduleMap.init(allocator, land_management_assignments.?, land_management_unit_by_cell.?, disturbance_catalog);
        irrigation_schedule_maps[0] = try ecosys.irrigation_management_dispatch.ScheduleMap.init(allocator, land_management_assignments.?, land_management_unit_by_cell.?, irrigation_catalog);
    }
    if (!ecosys.delimited_input.isNo(runscript.scenes[0].plant_management)) {
        const plant_assignment_path = try resolveInputPath(allocator, init.io, runscript_directory, runscript.scenes[0].plant_management);
        defer allocator.free(plant_assignment_path);
        const plant_assignment_source = try std.Io.Dir.cwd().readFileAlloc(init.io, plant_assignment_path, allocator, .limited(16 * 1024 * 1024));
        defer allocator.free(plant_assignment_source);
        plant_assignments = try ecosys.plant_assignment.parse(allocator, plant_assignment_source, site.ecosystem_type);
        plant_unit_by_cell = try plant_assignments.?.buildCellUnitMap(allocator, runscript.domain, runscript.plant_species_count);
        for (plant_assignments.?.units) |unit| for (unit.species) |species| {
            if (plant_catalog.find(species.species_file) == null) {
                const traits_path = try resolveInputPath(allocator, init.io, runscript_directory, species.species_file);
                defer allocator.free(traits_path);
                const traits_source = try std.Io.Dir.cwd().readFileAlloc(init.io, traits_path, allocator, .limited(16 * 1024 * 1024));
                defer allocator.free(traits_source);
                _ = try plant_catalog.appendFromSource(species.species_file, traits_source);
            }
            if (plant_management_catalog.find(species.management_file) == null and !ecosys.delimited_input.isNo(species.management_file)) {
                const management_path = try resolveInputPath(allocator, init.io, runscript_directory, species.management_file);
                defer allocator.free(management_path);
                const management_source = try std.Io.Dir.cwd().readFileAlloc(init.io, management_path, allocator, .limited(16 * 1024 * 1024));
                defer allocator.free(management_source);
                _ = try plant_management_catalog.appendFromSource(species.management_file, management_source);
            }
        };
    }

    // Validate later scene assignments and populate the shared input catalogs.
    // The first scene maps remain alive for state initialization; subsequent
    // maps are checked and released immediately.
    for (runscript.scenes[1..], 1..) |scene, scene_index| {
        if (!ecosys.delimited_input.isNo(scene.land_management)) {
            const path = try resolveInputPath(allocator, init.io, runscript_directory, scene.land_management);
            defer allocator.free(path);
            const assignment_source = try std.Io.Dir.cwd().readFileAlloc(init.io, path, allocator, .limited(16 * 1024 * 1024));
            defer allocator.free(assignment_source);
            var assignments = ecosys.land_management.parse(allocator, assignment_source) catch |err| {
                std.log.err("scene {d} land-management validation failed for '{s}'", .{ scene_index + 1, scene.land_management });
                return err;
            };
            defer assignments.deinit();
            const map = try assignments.buildCellUnitMap(allocator, runscript.domain);
            defer allocator.free(map);
            for (assignments.units) |unit| {
                if (!ecosys.delimited_input.isNo(unit.tillage_file) and disturbance_catalog.find(unit.tillage_file) == null) {
                    const schedule_path = try resolveInputPath(allocator, init.io, runscript_directory, unit.tillage_file);
                    defer allocator.free(schedule_path);
                    const schedule_source = try std.Io.Dir.cwd().readFileAlloc(init.io, schedule_path, allocator, .limited(64 * 1024 * 1024));
                    defer allocator.free(schedule_source);
                    _ = try disturbance_catalog.appendFromSource(unit.tillage_file, schedule_source);
                }
                if (!ecosys.delimited_input.isNo(unit.fertilizer_file) and fertilizer_catalog.find(unit.fertilizer_file) == null) {
                    const schedule_path = try resolveInputPath(allocator, init.io, runscript_directory, unit.fertilizer_file);
                    defer allocator.free(schedule_path);
                    const schedule_source = try std.Io.Dir.cwd().readFileAlloc(init.io, schedule_path, allocator, .limited(64 * 1024 * 1024));
                    defer allocator.free(schedule_source);
                    _ = try fertilizer_catalog.appendFromSource(unit.fertilizer_file, schedule_source);
                }
                if (!ecosys.delimited_input.isNo(unit.irrigation_file) and irrigation_catalog.find(unit.irrigation_file) == null) {
                    const schedule_path = try resolveInputPath(allocator, init.io, runscript_directory, unit.irrigation_file);
                    defer allocator.free(schedule_path);
                    const schedule_source = try std.Io.Dir.cwd().readFileAlloc(init.io, schedule_path, allocator, .limited(64 * 1024 * 1024));
                    defer allocator.free(schedule_source);
                    _ = try irrigation_catalog.appendFromSource(unit.irrigation_file, schedule_source);
                }
            }
            fertilizer_schedule_maps[scene_index] = try ecosys.fertilizer_management_dispatch.ScheduleMap.init(allocator, assignments, map, fertilizer_catalog);
            disturbance_schedule_maps[scene_index] = try ecosys.disturbance_management_dispatch.ScheduleMap.init(allocator, assignments, map, disturbance_catalog);
            irrigation_schedule_maps[scene_index] = try ecosys.irrigation_management_dispatch.ScheduleMap.init(allocator, assignments, map, irrigation_catalog);
        }
        if (!ecosys.delimited_input.isNo(scene.plant_management)) {
            const path = try resolveInputPath(allocator, init.io, runscript_directory, scene.plant_management);
            defer allocator.free(path);
            const assignment_source = try std.Io.Dir.cwd().readFileAlloc(init.io, path, allocator, .limited(16 * 1024 * 1024));
            defer allocator.free(assignment_source);
            var assignments = ecosys.plant_assignment.parse(allocator, assignment_source, site.ecosystem_type) catch |err| {
                std.log.err("scene {d} plant-management validation failed for '{s}'", .{ scene_index + 1, scene.plant_management });
                return err;
            };
            defer assignments.deinit();
            const map = try assignments.buildCellUnitMap(allocator, runscript.domain, runscript.plant_species_count);
            defer allocator.free(map);
            for (assignments.units) |unit| for (unit.species) |species| {
                if (plant_catalog.find(species.species_file) == null) {
                    const traits_path = try resolveInputPath(allocator, init.io, runscript_directory, species.species_file);
                    defer allocator.free(traits_path);
                    const traits_source = try std.Io.Dir.cwd().readFileAlloc(init.io, traits_path, allocator, .limited(16 * 1024 * 1024));
                    defer allocator.free(traits_source);
                    _ = try plant_catalog.appendFromSource(species.species_file, traits_source);
                }
                if (!ecosys.delimited_input.isNo(species.management_file) and plant_management_catalog.find(species.management_file) == null) {
                    const management_path = try resolveInputPath(allocator, init.io, runscript_directory, species.management_file);
                    defer allocator.free(management_path);
                    const management_source = try std.Io.Dir.cwd().readFileAlloc(init.io, management_path, allocator, .limited(16 * 1024 * 1024));
                    defer allocator.free(management_source);
                    _ = try plant_management_catalog.appendFromSource(species.management_file, management_source);
                }
            };
        }
    }

    const config = try ecosys.config.SimulationConfig.init(
        .{
            .grid_columns = try runscript.domain.columns(),
            .grid_rows = try runscript.domain.rows(),
            .soil_layers = try soil_catalog.maximumLayerCount(),
            .plant_populations = runscript.plant_species_count,
        },
        .{ .worker_threads = runscript.worker_count, .tile_cells = runscript.tile_cell_count },
        .{
            .relative_tolerance = runscript.relative_tolerance,
            .absolute_tolerance = runscript.absolute_tolerance,
            .max_nonlinear_iterations = iteration_limits.water_heat_solute_max_iterations,
            .picard_relaxation = runscript.picard_relaxation,
        },
    );
    const resolved_output_editors = try allocator.alloc(?ecosys.output_editor_layout.Resolved, try std.math.mul(usize, runscript.scenes.len, 10));
    @memset(resolved_output_editors, null);
    defer {
        for (resolved_output_editors) |*maybe_editor| if (maybe_editor.*) |*editor| editor.deinit();
        allocator.free(resolved_output_editors);
    }
    for (runscript.scenes, 0..) |scene, scene_index| for (scene.output_editors, 0..) |editor_name, editor_index| {
        if (ecosys.delimited_input.isNo(editor_name)) continue;
        const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
        var soil_catalog_for_editor = try soilOutputCatalog(allocator, editor_index, config.soil_layers);
        defer soil_catalog_for_editor.deinit();
        var plant_catalog_for_editor = try plantOutputCatalog(allocator, editor_index, config.soil_layers);
        defer plant_catalog_for_editor.deinit();
        resolved_output_editors[scene_index * 10 + editor_index] = try ecosys.output_editor_layout.resolve(
            allocator,
            output_selection_catalog.entries.items[selection_index].selection.enabled_variables,
            soil_catalog_for_editor.variables.len,
            plant_catalog_for_editor.variables.len,
        );
    };
    var state = try ecosys.grid.GridState.init(allocator, config);
    defer state.deinit();
    var any_visualization_enabled = false;
    for (scene_options) |options|
        any_visualization_enabled =
            any_visualization_enabled or options.visualization_enabled;
    var visualization_streams = try ecosys.visualization_output.Streams.init(
        allocator,
        init.io,
        std.Io.Dir.cwd(),
        any_visualization_enabled,
        state.cell_count,
        config.plant_populations,
        config.soil_layers,
        scene_options[0].visualization_start_year,
        scene_options[0].visualization_end_year,
        64 * 1024,
    );
    defer visualization_streams.deinit();
    const visualization_species_respired_carbon_g =
        try allocator.alloc(f64, config.plant_populations);
    defer allocator.free(visualization_species_respired_carbon_g);
    const visualization_leaf_carbon_g =
        try allocator.alloc(f64, config.plant_populations);
    defer allocator.free(visualization_leaf_carbon_g);
    const visualization_stalk_carbon_g =
        try allocator.alloc(f64, config.plant_populations);
    defer allocator.free(visualization_stalk_carbon_g);
    const visualization_root_carbon_g =
        try allocator.alloc(f64, config.plant_populations);
    defer allocator.free(visualization_root_carbon_g);
    const visualization_root_profile_carbon_g =
        try allocator.alloc(f64, config.soil_layers);
    defer allocator.free(visualization_root_profile_carbon_g);
    const visualization_organic_carbon_g =
        try allocator.alloc(f64, config.soil_layers + 1);
    defer allocator.free(visualization_organic_carbon_g);
    var daily_soil_gas_flux = try ecosys.soil_daily_gas_flux.State.init(allocator, state.cell_count);
    defer daily_soil_gas_flux.deinit();
    var daily_canopy_gas_exchange = try ecosys.daily_canopy_gas_exchange.State.init(allocator, state.cell_count);
    defer daily_canopy_gas_exchange.deinit();
    var daily_heterotrophic_respiration = try ecosys.soil_daily_heterotrophic_respiration.State.init(allocator, state.cell_count, state.layer_count);
    defer daily_heterotrophic_respiration.deinit();
    var daily_carbon_export = try ecosys.soil_daily_carbon_export.State.init(allocator, state.cell_count);
    defer daily_carbon_export.deinit();
    var daily_phosphorus_export = try ecosys.soil_daily_phosphorus_export.State.init(allocator, state.cell_count);
    defer daily_phosphorus_export.deinit();
    var daily_nitrogen_export = try ecosys.soil_daily_nitrogen_export.State.init(allocator, state.cell_count);
    defer daily_nitrogen_export.deinit();
    var daily_water_ledger = try ecosys.soil_daily_water_ledger.State.init(allocator, state.cell_count);
    defer daily_water_ledger.deinit();
    var daily_heat_ledger = try ecosys.soil_daily_heat_ledger.State.init(allocator, state.cell_count, config.soil_layers);
    defer daily_heat_ledger.deinit();
    var landscape_mass_balance_state: ecosys.landscape_mass_balance_checkpoint.State = .{
        .boundary_ledger = .{},
        .monitor = null,
    };
    const landscape_soil_mass_Mg_scratch = try allocator.alloc(
        f64,
        state.layer_count,
    );
    defer allocator.free(landscape_soil_mass_Mg_scratch);
    @memset(landscape_soil_mass_Mg_scratch, 0);
    var daily_ecosystem_carbon = try ecosys.soil_daily_ecosystem_carbon.State.init(allocator, state.cell_count);
    defer daily_ecosystem_carbon.deinit();
    const daily_soil_fire_charcoal_production_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_fire_charcoal_production_g_c);
    @memset(daily_soil_fire_charcoal_production_g_c, 0);
    const daily_soil_fire_carbon_dioxide_emission_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_fire_carbon_dioxide_emission_g_c);
    @memset(daily_soil_fire_carbon_dioxide_emission_g_c, 0);
    const daily_soil_fire_methane_emission_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_fire_methane_emission_g_c);
    @memset(daily_soil_fire_methane_emission_g_c, 0);
    const daily_soil_fire_phosphorus_flux_g_p = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_fire_phosphorus_flux_g_p);
    @memset(daily_soil_fire_phosphorus_flux_g_p, 0);
    const daily_soil_fire_nitrogen_flux_g_n = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_fire_nitrogen_flux_g_n);
    @memset(daily_soil_fire_nitrogen_flux_g_n, 0);
    const daily_soil_combusted_nitrogen_g_n = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_combusted_nitrogen_g_n);
    @memset(daily_soil_combusted_nitrogen_g_n, 0);
    const daily_soil_combusted_phosphorus_g_p = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_soil_combusted_phosphorus_g_p);
    @memset(daily_soil_combusted_phosphorus_g_p, 0);
    const daily_microbial_phosphate_mineralization_g_p = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_microbial_phosphate_mineralization_g_p);
    @memset(daily_microbial_phosphate_mineralization_g_p, 0);
    const daily_microbial_nitrogen_mineralization_g_n = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_microbial_nitrogen_mineralization_g_n);
    @memset(daily_microbial_nitrogen_mineralization_g_n, 0);
    const daily_organic_fertilizer_carbon_input_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_organic_fertilizer_carbon_input_g_c);
    @memset(daily_organic_fertilizer_carbon_input_g_c, 0);
    const daily_biome_organic_carbon_input_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_biome_organic_carbon_input_g_c);
    @memset(daily_biome_organic_carbon_input_g_c, 0);
    const daily_organic_fertilizer_phosphorus_input_g_p = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_organic_fertilizer_phosphorus_input_g_p);
    @memset(daily_organic_fertilizer_phosphorus_input_g_p, 0);
    const daily_organic_fertilizer_nitrogen_input_g_n = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(daily_organic_fertilizer_nitrogen_input_g_n);
    @memset(daily_organic_fertilizer_nitrogen_input_g_n, 0);
    var plant_state = try ecosys.grid.PlantState.init(allocator, config);
    defer plant_state.deinit();
    var hourly_soil_carbon_catalog = try ecosys.soil_output_catalog.carbon(allocator, config.soil_layers, config.soil_layers, config.soil_layers);
    defer hourly_soil_carbon_catalog.deinit();
    var hourly_soil_carbon_enabled = false;
    for (runscript.scenes) |scene| hourly_soil_carbon_enabled = hourly_soil_carbon_enabled or !ecosys.delimited_input.isNo(scene.output_editors[0]);
    var hourly_soil_carbon_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_soil_carbon_enabled, state.cell_count, hourly_soil_carbon_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_soil_carbon_bank.deinit();
    var hourly_plant_carbon_catalog = try ecosys.plant_output_catalog.carbon(allocator);
    defer hourly_plant_carbon_catalog.deinit();
    var hourly_plant_carbon_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[0])) continue;
        const resolved = resolved_output_editors[scene_index * 10] orelse continue;
        hourly_plant_carbon_enabled = hourly_plant_carbon_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    const output_plant_count = try std.math.mul(usize, state.cell_count, config.plant_populations);
    var hourly_plant_carbon_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_plant_carbon_enabled, output_plant_count, hourly_plant_carbon_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_plant_carbon_bank.deinit();
    var hourly_soil_water_catalog = try ecosys.soil_output_catalog.water(allocator, config.soil_layers);
    defer hourly_soil_water_catalog.deinit();
    var hourly_soil_water_enabled = false;
    for (runscript.scenes) |scene| hourly_soil_water_enabled = hourly_soil_water_enabled or !ecosys.delimited_input.isNo(scene.output_editors[1]);
    var hourly_soil_water_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_soil_water_enabled, state.cell_count, hourly_soil_water_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_soil_water_bank.deinit();
    var hourly_plant_water_catalog = try ecosys.plant_output_catalog.water(allocator, config.soil_layers);
    defer hourly_plant_water_catalog.deinit();
    var hourly_plant_water_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[1])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 1] orelse continue;
        hourly_plant_water_enabled = hourly_plant_water_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var hourly_plant_water_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_plant_water_enabled, output_plant_count, hourly_plant_water_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_plant_water_bank.deinit();
    var daily_plant_water_catalog = try ecosys.plant_output_catalog.dailyWater(allocator);
    defer daily_plant_water_catalog.deinit();
    var daily_plant_water_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[6])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 6] orelse continue;
        daily_plant_water_enabled = daily_plant_water_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var daily_plant_water_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_plant_water_enabled, output_plant_count, daily_plant_water_catalog.variables.len, 64 * 1024, .tab);
    defer daily_plant_water_bank.deinit();
    var daily_plant_carbon_catalog = try ecosys.plant_output_catalog.dailyCarbon(allocator, config.soil_layers);
    defer daily_plant_carbon_catalog.deinit();
    var daily_plant_carbon_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[5])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 5] orelse continue;
        daily_plant_carbon_enabled = daily_plant_carbon_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var daily_plant_carbon_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_plant_carbon_enabled, output_plant_count, daily_plant_carbon_catalog.variables.len, 64 * 1024, .tab);
    defer daily_plant_carbon_bank.deinit();
    var daily_soil_carbon_catalog = try ecosys.soil_output_catalog.dailyCarbon(allocator, config.soil_layers);
    defer daily_soil_carbon_catalog.deinit();
    var daily_soil_carbon_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[5])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 5] orelse continue;
        daily_soil_carbon_enabled = daily_soil_carbon_enabled or std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null;
    }
    var daily_soil_carbon_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_soil_carbon_enabled, state.cell_count, daily_soil_carbon_catalog.variables.len, 64 * 1024, .tab);
    defer daily_soil_carbon_bank.deinit();
    var daily_soil_water_catalog = try ecosys.soil_output_catalog.dailyWater(allocator, config.soil_layers, config.soil_layers, @min(config.soil_layers, 10));
    defer daily_soil_water_catalog.deinit();
    var daily_soil_water_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[6])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 6] orelse continue;
        daily_soil_water_enabled = daily_soil_water_enabled or std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null;
    }
    var daily_soil_water_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_soil_water_enabled, state.cell_count, daily_soil_water_catalog.variables.len, 64 * 1024, .tab);
    defer daily_soil_water_bank.deinit();
    const daily_water_profile_workspace = try allocator.alloc(f64, 3 * config.soil_layers);
    defer allocator.free(daily_water_profile_workspace);
    var daily_soil_phosphorus_catalog = try ecosys.soil_output_catalog.dailyPhosphorus(allocator, config.soil_layers);
    defer daily_soil_phosphorus_catalog.deinit();
    var daily_soil_phosphorus_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[8])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 8] orelse continue;
        daily_soil_phosphorus_enabled = daily_soil_phosphorus_enabled or std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null;
    }
    var daily_soil_phosphorus_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_soil_phosphorus_enabled, state.cell_count, daily_soil_phosphorus_catalog.variables.len, 64 * 1024, .tab);
    defer daily_soil_phosphorus_bank.deinit();
    var daily_soil_nitrogen_catalog = try ecosys.soil_output_catalog.dailyNitrogen(allocator, config.soil_layers);
    defer daily_soil_nitrogen_catalog.deinit();
    var daily_soil_nitrogen_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[7])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 7] orelse continue;
        daily_soil_nitrogen_enabled = daily_soil_nitrogen_enabled or std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null;
    }
    var daily_soil_nitrogen_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_soil_nitrogen_enabled, state.cell_count, daily_soil_nitrogen_catalog.variables.len, 64 * 1024, .tab);
    defer daily_soil_nitrogen_bank.deinit();
    var daily_soil_heat_catalog = try ecosys.soil_output_catalog.dailyHeat(allocator, config.soil_layers, config.soil_layers);
    defer daily_soil_heat_catalog.deinit();
    var daily_soil_heat_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[9])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 9] orelse continue;
        daily_soil_heat_enabled = daily_soil_heat_enabled or std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null;
    }
    var daily_soil_heat_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_soil_heat_enabled, state.cell_count, daily_soil_heat_catalog.variables.len, 64 * 1024, .tab);
    defer daily_soil_heat_bank.deinit();
    const daily_soil_heat_workspace = try allocator.alloc(f64, 3 * config.soil_layers);
    defer allocator.free(daily_soil_heat_workspace);
    const daily_soil_organic_carbon_by_layer_g_c = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(daily_soil_organic_carbon_by_layer_g_c);
    @memset(daily_soil_organic_carbon_by_layer_g_c, 0);
    const daily_soil_carbonate_mol_per_m3 = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(daily_soil_carbonate_mol_per_m3);
    const daily_soil_bicarbonate_mol_per_m3 = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(daily_soil_bicarbonate_mol_per_m3);
    const daily_soil_carbonate_complexes_mol_per_m3 = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(daily_soil_carbonate_complexes_mol_per_m3);
    const daily_soil_calcite_mol_per_m3 = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(daily_soil_calcite_mol_per_m3);
    var daily_plant_nitrogen_catalog = try ecosys.plant_output_catalog.dailyNitrogen(allocator);
    defer daily_plant_nitrogen_catalog.deinit();
    var daily_plant_nitrogen_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[7])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 7] orelse continue;
        daily_plant_nitrogen_enabled = daily_plant_nitrogen_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var daily_plant_nitrogen_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_plant_nitrogen_enabled, output_plant_count, daily_plant_nitrogen_catalog.variables.len, 64 * 1024, .tab);
    defer daily_plant_nitrogen_bank.deinit();
    var daily_plant_phosphorus_catalog = try ecosys.plant_output_catalog.dailyPhosphorus(allocator);
    defer daily_plant_phosphorus_catalog.deinit();
    var daily_plant_phosphorus_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[8])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 8] orelse continue;
        daily_plant_phosphorus_enabled = daily_plant_phosphorus_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var daily_plant_phosphorus_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_plant_phosphorus_enabled, output_plant_count, daily_plant_phosphorus_catalog.variables.len, 64 * 1024, .tab);
    defer daily_plant_phosphorus_bank.deinit();
    const daily_element_root_workspace = try allocator.alloc(f64, try std.math.mul(usize, output_plant_count, config.soil_layers));
    defer allocator.free(daily_element_root_workspace);
    @memset(daily_element_root_workspace, 0);
    const zero_root_layer_flux_g_per_h = try allocator.alloc(f64, config.soil_layers);
    defer allocator.free(zero_root_layer_flux_g_per_h);
    @memset(zero_root_layer_flux_g_per_h, 0);
    var daily_plant_development_catalog = try ecosys.plant_output_catalog.dailyDevelopment(allocator);
    defer daily_plant_development_catalog.deinit();
    var daily_plant_development_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[9])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 9] orelse continue;
        daily_plant_development_enabled = daily_plant_development_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var daily_plant_development_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), daily_plant_development_enabled, output_plant_count, daily_plant_development_catalog.variables.len, 64 * 1024, .tab);
    defer daily_plant_development_bank.deinit();
    var hourly_soil_nitrogen_catalog = try ecosys.soil_output_catalog.nitrogen(allocator, config.soil_layers, config.soil_layers);
    defer hourly_soil_nitrogen_catalog.deinit();
    var hourly_soil_nitrogen_enabled = false;
    for (runscript.scenes) |scene| hourly_soil_nitrogen_enabled = hourly_soil_nitrogen_enabled or !ecosys.delimited_input.isNo(scene.output_editors[2]);
    var hourly_soil_nitrogen_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_soil_nitrogen_enabled, state.cell_count, hourly_soil_nitrogen_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_soil_nitrogen_bank.deinit();
    var hourly_plant_nitrogen_catalog = try ecosys.plant_output_catalog.nitrogen(allocator, config.soil_layers);
    defer hourly_plant_nitrogen_catalog.deinit();
    var hourly_plant_nitrogen_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[2])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 2] orelse continue;
        hourly_plant_nitrogen_enabled = hourly_plant_nitrogen_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var hourly_plant_nitrogen_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_plant_nitrogen_enabled, output_plant_count, hourly_plant_nitrogen_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_plant_nitrogen_bank.deinit();
    var hourly_soil_phosphorus_catalog = try ecosys.soil_output_catalog.phosphorus(allocator);
    defer hourly_soil_phosphorus_catalog.deinit();
    var hourly_soil_phosphorus_enabled = false;
    for (runscript.scenes) |scene| hourly_soil_phosphorus_enabled = hourly_soil_phosphorus_enabled or !ecosys.delimited_input.isNo(scene.output_editors[3]);
    var hourly_soil_phosphorus_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_soil_phosphorus_enabled, state.cell_count, hourly_soil_phosphorus_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_soil_phosphorus_bank.deinit();
    var hourly_plant_phosphorus_catalog = try ecosys.plant_output_catalog.phosphorus(allocator, config.soil_layers);
    defer hourly_plant_phosphorus_catalog.deinit();
    var hourly_plant_phosphorus_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[3])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 3] orelse continue;
        hourly_plant_phosphorus_enabled = hourly_plant_phosphorus_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var hourly_plant_phosphorus_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_plant_phosphorus_enabled, output_plant_count, hourly_plant_phosphorus_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_plant_phosphorus_bank.deinit();
    var hourly_soil_heat_catalog = try ecosys.soil_output_catalog.heat(allocator, config.soil_layers);
    defer hourly_soil_heat_catalog.deinit();
    var hourly_soil_heat_enabled = false;
    for (runscript.scenes) |scene| hourly_soil_heat_enabled = hourly_soil_heat_enabled or !ecosys.delimited_input.isNo(scene.output_editors[4]);
    var hourly_soil_heat_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_soil_heat_enabled, state.cell_count, hourly_soil_heat_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_soil_heat_bank.deinit();
    var hourly_plant_heat_catalog = try ecosys.plant_output_catalog.heat(allocator);
    defer hourly_plant_heat_catalog.deinit();
    var hourly_plant_heat_enabled = false;
    for (runscript.scenes, 0..) |scene, scene_index| {
        if (ecosys.delimited_input.isNo(scene.output_editors[4])) continue;
        const resolved = resolved_output_editors[scene_index * 10 + 4] orelse continue;
        hourly_plant_heat_enabled = hourly_plant_heat_enabled or std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null;
    }
    var hourly_plant_heat_bank = try ecosys.output_stream_bank.Bank.init(allocator, init.io, std.Io.Dir.cwd(), hourly_plant_heat_enabled, output_plant_count, hourly_plant_heat_catalog.variables.len, 64 * 1024, .tab);
    defer hourly_plant_heat_bank.deinit();
    try state.validateFinite();
    try plant_state.validateFinite();
    const executor = try ecosys.compute.CpuExecutor.init(allocator, config.worker_threads, config.tile_cells);
    var soil_thermal_state = try ecosys.soil_thermal.State.initMapped(
        allocator,
        state,
        soil_catalog.entries.items,
        catalog_index_by_cell,
        grid_environment.horizontal_cell_width_m,
        grid_environment.vertical_cell_width_m,
    );
    defer soil_thermal_state.deinit();
    var initialization_context = ecosys.model_initialization.MappedHydrologyContext{
        .grid = &state,
        .catalog_entries = soil_catalog.entries.items,
        .catalog_index_by_cell = catalog_index_by_cell,
        .horizontal_cell_width_m = grid_environment.horizontal_cell_width_m,
        .vertical_cell_width_m = grid_environment.vertical_cell_width_m,
    };
    try runKernelAcrossSerialTilePlan(
        executor,
        &tile_plan.?,
        &initialization_context,
        ecosys.model_initialization.initializeMappedHydrologyTile,
    );
    const site_by_cell = try allocator.alloc(ecosys.site.Site, state.cell_count);
    defer allocator.free(site_by_cell);
    for (site_by_cell, grid_environment.site_catalog_index_by_cell) |*selected_site, site_index|
        selected_site.* = site_catalog.entries.items[site_index].site;
    const initial_surface_elevation_m = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(initial_surface_elevation_m);
    for (site_by_cell, initial_surface_elevation_m) |selected_site, *elevation_m|
        elevation_m.* = selected_site.elevation_m;
    try terrain_hydrology_state.bindInitialSurfaceElevations(initial_surface_elevation_m);
    const surface_runoff_boundary_storage = try allocator.alloc(
        f64,
        try std.math.mul(usize, state.cell_count, 4),
    );
    defer allocator.free(surface_runoff_boundary_storage);
    const mean_annual_temperature_c_by_cell = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(mean_annual_temperature_c_by_cell);
    const mean_annual_temperature_k_by_cell = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(mean_annual_temperature_k_by_cell);
    const geothermal_enabled_by_cell = try allocator.alloc(bool, state.cell_count);
    defer allocator.free(geothermal_enabled_by_cell);
    const surface_runoff_boundary_fraction_by_direction = [4][]f64{
        surface_runoff_boundary_storage[0 * state.cell_count ..][0..state.cell_count],
        surface_runoff_boundary_storage[1 * state.cell_count ..][0..state.cell_count],
        surface_runoff_boundary_storage[2 * state.cell_count ..][0..state.cell_count],
        surface_runoff_boundary_storage[3 * state.cell_count ..][0..state.cell_count],
    };
    for (site_by_cell, 0..) |selected_site, cell| {
        mean_annual_temperature_c_by_cell[cell] =
            selected_site.mean_annual_air_temperature_c;
        mean_annual_temperature_k_by_cell[cell] =
            selected_site.mean_annual_air_temperature_c + 273.15;
        geothermal_enabled_by_cell[cell] =
            runscript.geothermal_controls.enabled and selected_site.ecosystem_type != -2;
        for (0..4) |direction| {
            surface_runoff_boundary_fraction_by_direction[direction][cell] =
                selected_site.surface_runoff_boundary_fraction[direction];
        }
    }
    var soil_boundary_topology_state = try ecosys.soil_boundary_topology.State.initMapped(
        allocator,
        &state,
        &terrain_hydrology_state,
        config.grid_columns,
        config.grid_rows,
        grid_environment.horizontal_cell_width_m,
        grid_environment.vertical_cell_width_m,
        site_by_cell,
    );
    defer soil_boundary_topology_state.deinit();
    var soil_solver_property_state = try ecosys.soil_solver_properties.State.initMapped(
        allocator,
        &state,
        soil_catalog.entries.items,
        catalog_index_by_cell,
        grid_environment.horizontal_cell_width_m,
        grid_environment.vertical_cell_width_m,
        runscript.soil_solver_parameters,
    );
    defer soil_solver_property_state.deinit();
    var soil_geometry_state = try ecosys.soil_layer_geometry.State.init(allocator, state.cell_count, config.soil_layers);
    defer soil_geometry_state.deinit();
    var surface_pond_domain_workspace = try ecosys.surface_pond_domain_transaction.Workspace.init(allocator, state.cell_count, config.soil_layers);
    defer surface_pond_domain_workspace.deinit();
    for (0..state.cell_count) |cell| {
        const first_layer: usize = 0;
        const active_count = state.active_soil_layer_count[cell];
        const layer_base = cell * config.soil_layers;
        try ecosys.soil_layer_geometry.initializeCell(&soil_geometry_state, cell, first_layer, soil_solver_property_state.layer_thickness_m[layer_base .. layer_base + active_count], 0, runscript.soil_geometry_parameters.minimum_layer_thickness_m);
    }
    const fertilizer_band_layer_upper_depth_m =
        try allocator.alloc(f64, state.layer_count);
    defer allocator.free(fertilizer_band_layer_upper_depth_m);
    const fertilizer_band_layer_lower_depth_m =
        try allocator.alloc(f64, state.layer_count);
    defer allocator.free(fertilizer_band_layer_lower_depth_m);
    for (0..state.cell_count) |cell| {
        const boundary_base = cell * (config.soil_layers + 1);
        const layer_base = cell * config.soil_layers;
        const surface_depth_m =
            soil_geometry_state.boundary_depth_m[boundary_base];
        for (0..config.soil_layers) |layer| {
            fertilizer_band_layer_upper_depth_m[layer_base + layer] =
                soil_geometry_state.boundary_depth_m[
                    boundary_base + layer
                ] - surface_depth_m;
            fertilizer_band_layer_lower_depth_m[layer_base + layer] =
                soil_geometry_state.boundary_depth_m[
                    boundary_base + layer + 1
                ] - surface_depth_m;
        }
    }
    var fertilizer_band_state =
        try ecosys.fertilizer_band_state.State.initFromPlantNutrientParameters(
            allocator,
            state.cell_count,
            config.soil_layers,
            soil_geometry_state.active_layer_count,
            fertilizer_band_layer_upper_depth_m,
            fertilizer_band_layer_lower_depth_m,
            soil_geometry_state.layer_thickness_m,
            runscript.plant_nutrient_initialization,
        );
    defer fertilizer_band_state.deinit();
    var soil_microbial_state = try ecosys.soil_microbial_state.State.init(allocator, state.cell_count, config.soil_layers, runscript.microbial_substrate_count, runscript.microbial_population_count);
    defer soil_microbial_state.deinit();
    const soil_nitrogen_process_units_per_layer = try std.math.mul(usize, runscript.microbial_substrate_count, runscript.microbial_population_count);
    var soil_reactive_nitrogen_state = try ecosys.soil_reactive_nitrogen_state.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_reactive_nitrogen_state.deinit();
    var soil_microbial_phosphorus_state = try ecosys.soil_microbial_phosphorus_state.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_microbial_phosphorus_state.deinit();
    var soil_microbial_turnover_state = try ecosys.soil_microbial_turnover_step.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_microbial_turnover_state.deinit();
    var soil_litter_colonization_state = try ecosys.soil_litter_colonization_step.State.init(allocator, state.layer_count);
    defer soil_litter_colonization_state.deinit();
    var soil_organic_sorption_state = try ecosys.soil_organic_sorption_step.State.init(allocator, state.layer_count);
    defer soil_organic_sorption_state.deinit();
    var soil_organic_decomposition_state = try ecosys.soil_organic_decomposition_step.State.init(allocator, state.layer_count);
    defer soil_organic_decomposition_state.deinit();
    var soil_organic_priming_state = try ecosys.soil_organic_priming_step.State.init(allocator, state.layer_count, runscript.microbial_population_count);
    defer soil_organic_priming_state.deinit();
    var soil_respiration_products_state = try ecosys.soil_respiration_products_step.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_respiration_products_state.deinit();
    var soil_autotrophic_carbon_state = try ecosys.soil_autotrophic_carbon_step.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_autotrophic_carbon_state.deinit();
    var soil_microbial_layer_mixing_state = try ecosys.soil_microbial_layer_mixing.State.init(allocator, state.layer_count);
    defer soil_microbial_layer_mixing_state.deinit();
    var soil_biogeochemical_gas_flux_state = try ecosys.soil_biogeochemical_gas_aggregation.State.init(allocator, state.layer_count);
    defer soil_biogeochemical_gas_flux_state.deinit();
    var soil_methane_state = try ecosys.soil_methane_step.State.init(allocator, state.layer_count);
    defer soil_methane_state.deinit();
    var organic_matter_fire_exchange_state = try ecosys.organic_matter_fire_exchange.State.init(allocator, state.layer_count, ecosys.soil_organic_initialization.microbial_substrate_count);
    defer organic_matter_fire_exchange_state.deinit();
    var surface_fire_exchange_state = try ecosys.organic_matter_fire_exchange.State.init(allocator, state.cell_count, ecosys.soil_organic_initialization.microbial_substrate_count);
    defer surface_fire_exchange_state.deinit();
    var soil_microbial_oxygen_state = try ecosys.soil_oxygen_allocation.State.init(allocator, state.cell_count, config.soil_layers, soil_nitrogen_process_units_per_layer);
    defer soil_microbial_oxygen_state.deinit();
    var soil_oxygen_staging_state = try ecosys.soil_oxygen_step.State.init(allocator, state.cell_count, config.soil_layers, soil_nitrogen_process_units_per_layer);
    defer soil_oxygen_staging_state.deinit();
    var soil_nitrogen_flux_workspace = try ecosys.soil_nitrogen_flux_workspace.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_nitrogen_flux_workspace.deinit();
    var soil_nitrifier_environment_state = try ecosys.soil_nitrifier_environment_step.State.init(allocator, state.layer_count, soil_nitrogen_process_units_per_layer);
    defer soil_nitrifier_environment_state.deinit();
    const soil_field_capacity_fraction = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(soil_field_capacity_fraction);
    for (soil_solver_property_state.retention_curve, soil_field_capacity_fraction) |curve, *field_capacity| field_capacity.* = curve.curve.field_capacity_fraction;
    const soil_redox_satisfaction_fraction = try allocator.alloc(f64, try std.math.mul(usize, state.layer_count, soil_nitrogen_process_units_per_layer));
    defer allocator.free(soil_redox_satisfaction_fraction);
    @memset(soil_redox_satisfaction_fraction, 1);
    // Runtime transport domains. Capacities derive only from parsed grid and
    // profile inputs; there is no PARAMETERS.H-equivalent extent.
    const aqueous_species_count = ecosys.solute_transport_species.AqueousSpecies.count;
    var micropore_solute_state = try ecosys.solute_transport.State.init(allocator, state.layer_count, aqueous_species_count);
    defer micropore_solute_state.deinit();
    var macropore_solute_state = try ecosys.solute_transport.State.init(allocator, state.layer_count, aqueous_species_count);
    defer macropore_solute_state.deinit();
    const soil_recharge_concentration_mol_per_m3 = try allocator.alloc(f64, try std.math.mul(usize, state.layer_count, aqueous_species_count));
    defer allocator.free(soil_recharge_concentration_mol_per_m3);
    const soil_solute_boundary_net_flux_mol = try allocator.alloc(f64, try std.math.mul(usize, state.layer_count, aqueous_species_count));
    defer allocator.free(soil_solute_boundary_net_flux_mol);
    @memset(soil_solute_boundary_net_flux_mol, 0);
    var gas_transport_state = try ecosys.gas_transport.State.init(allocator, state.layer_count);
    defer gas_transport_state.deinit();
    var soil_gas_transport_state = try ecosys.soil_gas_transport_step.State.init(allocator, state.layer_count);
    defer soil_gas_transport_state.deinit();
    var surface_litter_gas_transport_state = try ecosys.surface_litter_gas_transport_step.State.init(allocator, state.cell_count);
    defer surface_litter_gas_transport_state.deinit();
    var litter_gas_transport_state = try ecosys.gas_transport.State.init(allocator, state.cell_count);
    defer litter_gas_transport_state.deinit();
    var plant_available_nutrient_state = try ecosys.soil_plant_available_nutrients.State.init(allocator, state.layer_count);
    defer plant_available_nutrient_state.deinit();
    try plant_available_nutrient_state.initializeMapped(&state, &soil_solver_property_state, soil_catalog.entries.items, catalog_index_by_cell, runscript.plant_nutrient_initialization);
    var soil_organic_state = try ecosys.soil_organic_initialization.State.init(allocator, state.layer_count);
    defer soil_organic_state.deinit();
    const soil_organic_carbon_at_hour_start_g_c = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(soil_organic_carbon_at_hour_start_g_c);
    const soil_organic_carbon_change_g_c_per_h = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(soil_organic_carbon_change_g_c_per_h);
    @memset(soil_organic_carbon_at_hour_start_g_c, 0);
    @memset(soil_organic_carbon_change_g_c_per_h, 0);
    var surface_organic_state = try ecosys.soil_organic_initialization.State.init(allocator, state.cell_count);
    defer surface_organic_state.deinit();
    const shoot_senescence_products_by_plant = try allocator.alloc(ecosys.canopy_photosynthesis.SenescenceProducts, output_plant_count);
    defer allocator.free(shoot_senescence_products_by_plant);
    @memset(shoot_senescence_products_by_plant, .{});
    const root_litter_products_by_plant = try allocator.alloc(ecosys.plant_root_metabolism.RootLitter, output_plant_count);
    defer allocator.free(root_litter_products_by_plant);
    @memset(root_litter_products_by_plant, std.mem.zeroes(ecosys.plant_root_metabolism.RootLitter));
    var root_litter_carbon_ledger = try ecosys.plant_root_litter_ledger.State.init(
        allocator,
        output_plant_count,
        ecosys.plant_root_system.biological_domain_count,
        config.soil_layers,
    );
    defer root_litter_carbon_ledger.deinit();
    var plant_salt_harvest_workspace = try ecosys.plant_salt_harvest_adapter.Workspace.init(
        allocator,
        config.soil_layers,
    );
    defer plant_salt_harvest_workspace.deinit();
    const cumulative_harvest_salt_mol_by_plant = try allocator.alloc(
        f64,
        output_plant_count * ecosys.plant_salt_harvest.salt_count,
    );
    defer allocator.free(cumulative_harvest_salt_mol_by_plant);
    @memset(cumulative_harvest_salt_mol_by_plant, 0);
    const litter_salt_mol_by_plant_exchange_layer = try allocator.alloc(
        f64,
        output_plant_count * (config.soil_layers + 1) * ecosys.plant_salt_harvest.salt_count,
    );
    defer allocator.free(litter_salt_mol_by_plant_exchange_layer);
    @memset(litter_salt_mol_by_plant_exchange_layer, 0);
    var plant_litter_salt_ingress_state = try ecosys.plant_litter_salt_ingress.State.init(
        allocator,
        state.cell_count,
        config.soil_layers,
    );
    defer plant_litter_salt_ingress_state.deinit();
    var litter_salt_publication_state = try ecosys.litter_salt_publication.State.init(
        allocator,
        state.cell_count,
        config.plant_populations,
        config.soil_layers,
    );
    defer litter_salt_publication_state.deinit();
    const root_litter_carbon_g_c_by_layer = try allocator.alloc(f64, config.soil_layers);
    defer allocator.free(root_litter_carbon_g_c_by_layer);
    @memset(root_litter_carbon_g_c_by_layer, 0);
    const shoot_harvest_litter_carbon_g_c_by_plant = try allocator.alloc(f64, output_plant_count);
    defer allocator.free(shoot_harvest_litter_carbon_g_c_by_plant);
    @memset(shoot_harvest_litter_carbon_g_c_by_plant, 0);
    const hourly_manure_products_by_plant = try allocator.alloc(
        ecosys.grazing_manure.Products,
        output_plant_count,
    );
    defer allocator.free(hourly_manure_products_by_plant);
    @memset(hourly_manure_products_by_plant, .{});
    const shoot_harvest_litter_nitrogen_g_n_by_plant = try allocator.alloc(f64, output_plant_count);
    defer allocator.free(shoot_harvest_litter_nitrogen_g_n_by_plant);
    @memset(shoot_harvest_litter_nitrogen_g_n_by_plant, 0);
    const shoot_harvest_litter_phosphorus_g_p_by_plant = try allocator.alloc(f64, output_plant_count);
    defer allocator.free(shoot_harvest_litter_phosphorus_g_p_by_plant);
    @memset(shoot_harvest_litter_phosphorus_g_p_by_plant, 0);
    const harvest_carbon_at_hour_start_g_c_by_plant = try allocator.alloc(f64, output_plant_count);
    defer allocator.free(harvest_carbon_at_hour_start_g_c_by_plant);
    @memset(harvest_carbon_at_hour_start_g_c_by_plant, 0);
    const eroded_organic_component_count = try ecosys.soil_erosion_organic_bridge.componentCount(&soil_organic_state);
    var eroded_organic_workspace = try ecosys.eroded_constituents.PackedWorkspace.init(allocator, state.cell_count, eroded_organic_component_count);
    defer eroded_organic_workspace.deinit();
    var surface_litter_geometry_state = try ecosys.surface_litter_geometry_step.State.init(allocator, state.cell_count);
    defer surface_litter_geometry_state.deinit();
    var surface_litter_water_environment_state = try ecosys.surface_litter_water_environment.State.init(allocator, state.cell_count);
    defer surface_litter_water_environment_state.deinit();
    var surface_microbial_environment_state = try ecosys.surface_microbial_environment_step.State.init(allocator, state.cell_count);
    defer surface_microbial_environment_state.deinit();
    const surface_field_capacity_potential_mpa = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_field_capacity_potential_mpa);
    const surface_wilting_point_potential_mpa = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_wilting_point_potential_mpa);
    for (0..state.cell_count) |cell| {
        const profile = soil_catalog.entries.items[catalog_index_by_cell[cell]].profile;
        surface_field_capacity_potential_mpa[cell] = profile.field_capacity_potential_mpa;
        surface_wilting_point_potential_mpa[cell] = profile.wilting_point_potential_mpa;
    }
    const surface_litter_ice_m3 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_litter_ice_m3);
    @memset(surface_litter_ice_m3, 0);
    const surface_charcoal_carbon_g_c = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_charcoal_carbon_g_c);
    @memset(surface_charcoal_carbon_g_c, 0);
    var surface_litter_chemistry_state = try ecosys.surface_litter_chemistry.State.init(allocator, state.cell_count);
    defer surface_litter_chemistry_state.deinit();
    var surface_litter_fertilizer_state = try ecosys.surface_litter_fertilizer.State.init(allocator, state.cell_count);
    defer surface_litter_fertilizer_state.deinit();
    var surface_microbial_respiration_state = try ecosys.surface_microbial_respiration_step.State.init(allocator, state.cell_count);
    defer surface_microbial_respiration_state.deinit();
    var surface_microbial_oxygen_state = try ecosys.surface_microbial_oxygen_driver.State.init(allocator, state.cell_count);
    defer surface_microbial_oxygen_state.deinit();
    var surface_microbial_maintenance_state = try ecosys.surface_microbial_maintenance_step.State.init(allocator, state.cell_count);
    defer surface_microbial_maintenance_state.deinit();
    var surface_nonsymbiotic_nitrogen_fixation_state = try ecosys.surface_nonsymbiotic_nitrogen_fixation_step.State.init(allocator, state.cell_count);
    defer surface_nonsymbiotic_nitrogen_fixation_state.deinit();
    var surface_microbial_substrate_uptake_state = try ecosys.surface_microbial_substrate_uptake_step.State.init(allocator, state.cell_count);
    defer surface_microbial_substrate_uptake_state.deinit();
    var surface_denitrification_state = try ecosys.surface_denitrification_step.State.init(allocator, state.cell_count);
    defer surface_denitrification_state.deinit();
    var surface_microbial_assimilation_state = try ecosys.surface_microbial_assimilation_step.State.init(allocator, state.cell_count);
    defer surface_microbial_assimilation_state.deinit();
    var surface_microbial_mineral_exchange_state = try ecosys.surface_microbial_mineral_exchange_step.State.init(allocator, state.cell_count);
    defer surface_microbial_mineral_exchange_state.deinit();
    var surface_topsoil_mineral_exchange_state = try ecosys.surface_topsoil_mineral_exchange_step.State.init(allocator, state.cell_count);
    defer surface_topsoil_mineral_exchange_state.deinit();
    var surface_microbial_turnover_state = try ecosys.surface_microbial_turnover_step.State.init(allocator, state.cell_count);
    defer surface_microbial_turnover_state.deinit();
    var surface_organic_priming_state = try ecosys.surface_organic_priming_step.State.init(allocator, state.cell_count);
    defer surface_organic_priming_state.deinit();
    var surface_organic_decomposition_state = try ecosys.surface_organic_decomposition_step.State.init(allocator, state.cell_count);
    defer surface_organic_decomposition_state.deinit();
    var surface_organic_sorption_state = try ecosys.surface_organic_sorption_step.State.init(allocator, state.cell_count);
    defer surface_organic_sorption_state.deinit();
    var surface_litter_colonization_state = try ecosys.surface_litter_colonization_step.State.init(allocator, state.cell_count);
    defer surface_litter_colonization_state.deinit();
    const surface_humification_fraction = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_humification_fraction);
    const topsoil_humus_partition = try allocator.alloc([2]f64, state.cell_count);
    defer allocator.free(topsoil_humus_partition);
    for (0..state.cell_count) |cell| {
        const topsoil_material = soil_catalog.entries.items[catalog_index_by_cell[cell]].material;
        surface_humification_fraction[cell] = 0.150 + 0.300 * @min(0.333, topsoil_material.clay_mass_fraction[0]);
        // STARTS FCY=0.20 when no initialized humus inventory is available.
        // NITRO converts that humus fraction to its microbial-residue routing
        // fraction through 5*FC/(4*FC+1).
        const source_less_resistant_humus_fraction: f64 = 0.20;
        const source_microbial_litter_partition =
            5 * source_less_resistant_humus_fraction /
            (4 * source_less_resistant_humus_fraction + 1);
        topsoil_humus_partition[cell] = .{
            source_microbial_litter_partition,
            1 - source_microbial_litter_partition,
        };
    }
    var surface_litter_fertilizer_diagnostics = try ecosys.surface_litter_fertilizer_step.Diagnostics.init(allocator, state.cell_count);
    defer surface_litter_fertilizer_diagnostics.deinit();
    var soil_fertilizer_inventory = try ecosys.fertilizer_nitrogen_inventory.State.init(allocator, state.cell_count, config.soil_layers);
    defer soil_fertilizer_inventory.deinit();
    var mineral_fertilizer_inventory = try ecosys.mineral_fertilizer_inventory.State.init(allocator, state.cell_count, config.soil_layers);
    defer mineral_fertilizer_inventory.deinit();
    var eroded_fertilizer_workspace = try ecosys.eroded_constituents.PackedWorkspace.init(allocator, state.cell_count, ecosys.soil_erosion_fertilizer_bridge.component_count);
    defer eroded_fertilizer_workspace.deinit();
    var surface_litter_chemistry_diagnostics = try ecosys.surface_litter_chemistry_step.Diagnostics.init(allocator, state.cell_count);
    defer surface_litter_chemistry_diagnostics.deinit();
    const surface_litter_cation_selectivity = try allocator.alloc(ecosys.solute_cation_exchange.Selectivity, state.cell_count);
    defer allocator.free(surface_litter_cation_selectivity);
    for (surface_litter_cation_selectivity, 0..) |*selectivity, cell| {
        const profile = soil_catalog.entries.items[catalog_index_by_cell[cell]].profile;
        selectivity.* = .{
            .calcium_ammonium = profile.property(.gapon_calcium_ammonium)[0],
            .calcium_hydrogen = profile.property(.gapon_calcium_hydrogen)[0],
            .calcium_aluminum_and_iron = profile.property(.gapon_calcium_aluminum)[0],
            .calcium_magnesium = profile.property(.gapon_calcium_magnesium)[0],
            .calcium_sodium = profile.property(.gapon_calcium_sodium)[0],
            .calcium_potassium = profile.property(.gapon_calcium_potassium)[0],
        };
    }
    if (organic_parameters) |*parameters| {
        const depth_partition_factor = try allocator.alloc(f64, state.layer_count);
        defer allocator.free(depth_partition_factor);
        try ecosys.soil_organic_initialization.deriveMappedDepthPartition(
            allocator,
            &state,
            &soil_solver_property_state,
            soil_catalog.entries.items,
            catalog_index_by_cell,
            grid_environment.horizontal_cell_width_m,
            grid_environment.vertical_cell_width_m,
            terrain_hydrology_state.relative_surface_elevation_m,
            terrain_hydrology_state.minimum_surface_elevation_m,
            grid_environment.initial_water_table_depth_m,
            grid_environment.natural_water_table_surface_slope,
            parameters.humusDepthParameters(),
            depth_partition_factor,
        );
        try soil_organic_state.initializeMapped(&state, &soil_solver_property_state, soil_catalog.entries.items, catalog_index_by_cell, parameters.mapped(depth_partition_factor));
        try plant_available_nutrient_state.bindInitializedDissolvedOrganic(&soil_organic_state);
        try surface_organic_state.initializeMappedSurface(
            soil_catalog.entries.items,
            catalog_index_by_cell,
            grid_environment.horizontal_cell_width_m,
            grid_environment.vertical_cell_width_m,
            parameters.surface(),
        );
        for (0..state.cell_count) |cell| {
            const top = try state.layerIndex(cell, 0);
            const humus_base = (top * ecosys.soil_organic_initialization.substrate_count + 4) * ecosys.soil_organic_initialization.structural_fraction_count;
            const less_resistant_carbon_g_c = soil_organic_state.structural[humus_base].carbon_g_c;
            const more_resistant_carbon_g_c = soil_organic_state.structural[humus_base + 1].carbon_g_c;
            const total_humus_carbon_g_c = less_resistant_carbon_g_c + more_resistant_carbon_g_c;
            const less_resistant_fraction = if (total_humus_carbon_g_c > 0) less_resistant_carbon_g_c / total_humus_carbon_g_c else parameters.less_resistant_humus_fraction_at_surface;
            const microbial_litter_to_less_resistant = 5 * less_resistant_fraction / (4 * less_resistant_fraction + 1);
            topsoil_humus_partition[cell] = .{ microbial_litter_to_less_resistant, 1 - microbial_litter_to_less_resistant };
        }
    }
    for (0..state.layer_count) |layer| for (0..@min(runscript.microbial_substrate_count, ecosys.soil_organic_initialization.microbial_substrate_count)) |substrate| for (0..@min(runscript.microbial_population_count, ecosys.soil_organic_initialization.microbial_population_count)) |population| {
        const runtime_index = try soil_microbial_state.populationIndex(layer / config.soil_layers, layer % config.soil_layers, substrate, population);
        const organic_index = ((layer * ecosys.soil_organic_initialization.microbial_substrate_count + substrate) * ecosys.soil_organic_initialization.microbial_population_count + population) * ecosys.soil_organic_initialization.kinetic_fraction_count;
        const initialized_labile = soil_organic_state.microbial[organic_index];
        const initialized_resistant = soil_organic_state.microbial[organic_index + 1];
        const initialized_nonstructural = soil_organic_state.microbial[organic_index + 2];
        soil_microbial_state.structural[runtime_index * 2] = .{ .carbon_g_c = initialized_labile.carbon_g_c, .nitrogen_g_n = initialized_labile.nitrogen_g_n, .phosphorus_g_p = initialized_labile.phosphorus_g_p };
        soil_microbial_state.structural[runtime_index * 2 + 1] = .{ .carbon_g_c = initialized_resistant.carbon_g_c, .nitrogen_g_n = initialized_resistant.nitrogen_g_n, .phosphorus_g_p = initialized_resistant.phosphorus_g_p };
        soil_microbial_state.nonstructural[runtime_index] = .{ .carbon_g_c = initialized_nonstructural.carbon_g_c, .nitrogen_g_n = initialized_nonstructural.nitrogen_g_n, .phosphorus_g_p = initialized_nonstructural.phosphorus_g_p };
    };
    var initial_chemistry_state = try ecosys.solute_chemistry_state.State.init(allocator, state.layer_count);
    defer initial_chemistry_state.deinit();
    const soil_chemistry_layer_parameters = try allocator.alloc(ecosys.solute_chemistry_state.ReactionParameters, state.layer_count);
    defer allocator.free(soil_chemistry_layer_parameters);
    var soil_chemistry_solver_workspace = try ecosys.solute_reaction_solver.Workspace.init(allocator);
    defer soil_chemistry_solver_workspace.deinit();
    var eroded_chemistry_workspace = try ecosys.eroded_constituents.PackedWorkspace.init(allocator, state.cell_count, ecosys.soil_erosion_chemistry_bridge.component_count);
    defer eroded_chemistry_workspace.deinit();
    const chemistry_parameters = runscript.chemistry_initialization orelse
        ecosys.soil_chemistry_initialization.sourceParameters();
    {
        for (0..state.cell_count) |cell| {
            const soil_entry = soil_catalog.entries.items[catalog_index_by_cell[cell]];
            const profile = soil_entry.profile;
            for (0..profile.total_layer_count) |layer| {
                const layer_cell = try state.layerIndex(cell, layer);
                if (runscript.chemistry_primary_initialization) |primary_parameters| {
                    const profile_chemistry_inputs = try ecosys.soil_chemistry_initialization.resolveProfileEquilibriumSentinels(.{
                        .soil_ph = profile.ph[layer],
                        .ammonium_g_n_per_megagram = soil_entry.material.initial_ammonium_g_per_megagram[layer],
                        .nitrate_g_n_per_megagram = profile.property(.nitrate_g_per_megagram)[layer],
                        .phosphate_g_p_per_megagram = profile.property(.phosphate_g_per_megagram)[layer],
                        .aluminum_g_per_megagram = profile.property(.aluminum_g_per_megagram)[layer],
                        .iron_g_per_megagram = profile.property(.iron_g_per_megagram)[layer],
                        .calcium_g_per_megagram = soil_entry.material.initial_calcium_g_per_megagram[layer],
                        .magnesium_g_per_megagram = profile.property(.magnesium_g_per_megagram)[layer],
                        .sodium_g_per_megagram = profile.property(.sodium_g_per_megagram)[layer],
                        .potassium_g_per_megagram = profile.property(.potassium_g_per_megagram)[layer],
                        .sulfate_sulfur_g_s_per_megagram = profile.property(.sulfate_sulfur_g_per_megagram)[layer],
                        .chloride_g_per_megagram = profile.property(.chloride_g_per_megagram)[layer],
                        .aluminum_phosphate_g_p_per_megagram = profile.property(.aluminum_phosphate_p_g_per_megagram)[layer],
                        .iron_phosphate_g_p_per_megagram = profile.property(.iron_phosphate_p_g_per_megagram)[layer],
                        .dicalcium_phosphate_g_p_per_megagram = profile.property(.calcium_hydrogen_phosphate_p_g_per_megagram)[layer],
                        .apatite_g_p_per_megagram = profile.property(.apatite_phosphorus_g_per_megagram)[layer],
                        .aluminum_hydroxide_g_al_per_megagram = profile.property(.aluminum_hydroxide_al_g_per_megagram)[layer],
                        .iron_hydroxide_g_fe_per_megagram = profile.property(.iron_hydroxide_fe_g_per_megagram)[layer],
                        .calcium_carbonate_g_ca_per_megagram = profile.property(.calcium_carbonate_ca_g_per_megagram)[layer],
                        .calcium_sulfate_g_ca_per_megagram = profile.property(.calcium_sulfate_ca_g_per_megagram)[layer],
                    }, primary_parameters);
                    try ecosys.soil_chemistry_initialization.seedProfilePrimaryState(&initial_chemistry_state, layer_cell, profile_chemistry_inputs, soil_solver_property_state.matrix_bulk_volume_m3[layer_cell] * soil_solver_property_state.bulk_density_megagrams_per_m3[layer_cell], state.matrix_liquid_water_m3[layer_cell], primary_parameters);
                } else {
                    try ecosys.soil_chemistry_initialization.seedProfilePhosphate(&initial_chemistry_state, layer_cell, profile.ph[layer], profile.property(.phosphate_g_per_megagram)[layer], chemistry_parameters);
                }
            }
        }
        if (chemistry_reaction_parameters) |reaction_parameters| {
            const zone_volume_storage = try allocator.alloc(f64, 6 * state.layer_count);
            defer allocator.free(zone_volume_storage);
            @memset(zone_volume_storage, 0);
            const ammonium_non_band_volume_m3 = zone_volume_storage[0 * state.layer_count ..][0..state.layer_count];
            const ammonium_band_volume_m3 = zone_volume_storage[1 * state.layer_count ..][0..state.layer_count];
            const nitrate_non_band_volume_m3 = zone_volume_storage[2 * state.layer_count ..][0..state.layer_count];
            const nitrate_band_volume_m3 = zone_volume_storage[3 * state.layer_count ..][0..state.layer_count];
            const phosphate_non_band_volume_m3 = zone_volume_storage[4 * state.layer_count ..][0..state.layer_count];
            const phosphate_band_volume_m3 = zone_volume_storage[5 * state.layer_count ..][0..state.layer_count];
            const nutrient_zones = runscript.plant_nutrient_initialization;
            const fractions = ecosys.solute_charge_classification.ZoneFractions{
                .ammonium_non_band = 1 - nutrient_zones.initial_ammonium_band_fraction,
                .ammonium_band = nutrient_zones.initial_ammonium_band_fraction,
                .nitrate_non_band = 1 - nutrient_zones.initial_nitrate_band_fraction,
                .nitrate_band = nutrient_zones.initial_nitrate_band_fraction,
                .phosphate_non_band = 1 - nutrient_zones.initial_phosphate_band_fraction,
                .phosphate_band = nutrient_zones.initial_phosphate_band_fraction,
            };
            for (0..state.cell_count) |cell| {
                const profile = soil_catalog.entries.items[catalog_index_by_cell[cell]].profile;
                for (0..profile.total_layer_count) |layer| {
                    const layer_cell = try state.layerIndex(cell, layer);
                    const water_volume_m3 = state.matrix_liquid_water_m3[layer_cell];
                    const soil_mass_Mg = soil_solver_property_state.matrix_bulk_volume_m3[layer_cell] * soil_solver_property_state.bulk_density_megagrams_per_m3[layer_cell];
                    if (!std.math.isFinite(water_volume_m3) or water_volume_m3 <= 0 or !std.math.isFinite(soil_mass_Mg) or soil_mass_Mg <= 0) {
                        std.log.err("STARTE requires positive runtime soil and matrix-water volumes: cell={d} layer={d} soil_mass_Mg={e} water_volume_m3={e}", .{ cell, layer, soil_mass_Mg, water_volume_m3 });
                        return error.InvalidStarteLayerGeometry;
                    }
                    ammonium_non_band_volume_m3[layer_cell] = water_volume_m3 * fractions.ammonium_non_band;
                    ammonium_band_volume_m3[layer_cell] = water_volume_m3 * fractions.ammonium_band;
                    nitrate_non_band_volume_m3[layer_cell] = water_volume_m3 * fractions.nitrate_non_band;
                    nitrate_band_volume_m3[layer_cell] = water_volume_m3 * fractions.nitrate_band;
                    phosphate_non_band_volume_m3[layer_cell] = water_volume_m3 * fractions.phosphate_non_band;
                    phosphate_band_volume_m3[layer_cell] = water_volume_m3 * fractions.phosphate_band;
                    const shared_soil_mass_per_water_volume_Mg_per_m3 = soil_mass_Mg / water_volume_m3;
                    const phosphate_non_band_ratio = if (fractions.phosphate_non_band > 0) shared_soil_mass_per_water_volume_Mg_per_m3 / fractions.phosphate_non_band else shared_soil_mass_per_water_volume_Mg_per_m3;
                    const phosphate_band_ratio = if (fractions.phosphate_band > 0) shared_soil_mass_per_water_volume_Mg_per_m3 / fractions.phosphate_band else shared_soil_mass_per_water_volume_Mg_per_m3;
                    const ammonium_non_band_ratio = if (fractions.ammonium_non_band > 0) shared_soil_mass_per_water_volume_Mg_per_m3 / fractions.ammonium_non_band else shared_soil_mass_per_water_volume_Mg_per_m3;
                    const ammonium_band_ratio = if (fractions.ammonium_band > 0) shared_soil_mass_per_water_volume_Mg_per_m3 / fractions.ammonium_band else shared_soil_mass_per_water_volume_Mg_per_m3;
                    const selectivity = ecosys.solute_cation_exchange.Selectivity{
                        .calcium_ammonium = profile.property(.gapon_calcium_ammonium)[layer],
                        .calcium_hydrogen = profile.property(.gapon_calcium_hydrogen)[layer],
                        .calcium_aluminum_and_iron = profile.property(.gapon_calcium_aluminum)[layer],
                        .calcium_magnesium = profile.property(.gapon_calcium_magnesium)[layer],
                        .calcium_sodium = profile.property(.gapon_calcium_sodium)[layer],
                        .calcium_potassium = profile.property(.gapon_calcium_potassium)[layer],
                    };
                    const cation_exchange_capacity_mol_charge_per_Mg = soil_solver_property_state.cation_exchange_capacity_mol_per_Mg[layer_cell];
                    const anion_exchange_capacity_mol_charge_per_Mg = 10 * profile.anion_exchange_capacity_cmol_kg[layer];
                    const total_carboxyl_sites_mol_per_Mg = reaction_parameters.surface_litter.carboxyl_sites_mol_per_Mg_c *
                        1.0e-6 * soil_solver_property_state.total_organic_carbon_g_per_megagram[layer_cell];
                    try ecosys.soil_chemistry_initialization.seedProfilePhosphateSurfaceSites(&initial_chemistry_state, layer_cell, profile.property(.phosphate_g_per_megagram)[layer] / 31.0, anion_exchange_capacity_mol_charge_per_Mg, reaction_parameters.phosphate_surface, reaction_parameters.phosphate_constants.h2po4);
                    try ecosys.soil_chemistry_initialization.seedProfileCationExchange(&initial_chemistry_state, layer_cell, cation_exchange_capacity_mol_charge_per_Mg, selectivity, fractions);
                    initial_chemistry_state.water_mol_per_m3[layer_cell] = reaction_parameters.water_concentration_mol_per_m3;
                    const layer_parameters = reaction_parameters.forLayer(fractions, phosphate_non_band_ratio, phosphate_band_ratio, cation_exchange_capacity_mol_charge_per_Mg, total_carboxyl_sites_mol_per_Mg, .{ .shared_Mg_per_m3 = shared_soil_mass_per_water_volume_Mg_per_m3, .ammonium_non_band_Mg_per_m3 = ammonium_non_band_ratio, .ammonium_band_Mg_per_m3 = ammonium_band_ratio }, selectivity);
                    soil_chemistry_layer_parameters[layer_cell] = layer_parameters;
                    const initialization_solver_options: ecosys.solute_reaction_solver.Options = .{
                        .absolute_tolerance = runscript.absolute_tolerance,
                        .relative_tolerance = runscript.relative_tolerance,
                        .picard_relaxation = runscript.picard_relaxation,
                        .max_iterations = iteration_limits.initial_solute_reaction_max_iterations,
                    };
                    _ = ecosys.solute_reaction_solver.solveCell(
                        allocator,
                        &initial_chemistry_state,
                        layer_cell,
                        layer_parameters,
                        initialization_solver_options,
                    ) catch |err| return ecosys.solute_failure_reporter
                        .reportPreservingSolverError(
                        allocator,
                        .{
                            .io = init.io,
                            .directory = std.Io.Dir.cwd(),
                            .file_path = "ecosys-ng-solute-initialization-failure.bin",
                            .context = .{
                                .global_cell_id = @intCast(cell),
                                .soil_layer_id = @intCast(layer),
                                .packed_cell_index = @intCast(layer_cell),
                            },
                        },
                        &initial_chemistry_state,
                        layer_cell,
                        layer_parameters,
                        initialization_solver_options,
                        err,
                    );
                }
            }
            try plant_available_nutrient_state.bindEquilibratedMineralPools(&initial_chemistry_state, .{
                .ammonium_non_band_m3 = ammonium_non_band_volume_m3,
                .ammonium_band_m3 = ammonium_band_volume_m3,
                .nitrate_non_band_m3 = nitrate_non_band_volume_m3,
                .nitrate_band_m3 = nitrate_band_volume_m3,
                .phosphate_non_band_m3 = phosphate_non_band_volume_m3,
                .phosphate_band_m3 = phosphate_band_volume_m3,
            });
        }
        for (0..state.cell_count) |cell| {
            const profile = soil_catalog.entries.items[catalog_index_by_cell[cell]].profile;
            try ecosys.soil_chemistry_initialization.seedSurfaceLitterFromTopsoil(
                &surface_litter_chemistry_state,
                cell,
                &initial_chemistry_state,
                try state.layerIndex(cell, 0),
                profile.surface_litter_ph,
                chemistry_parameters.water_activity_product_mol2_per_m6,
                runscript.dynamic_plant_salts,
            );
        }
    }
    for (0..state.layer_count) |layer_cell| {
        micropore_solute_state.water_volume_m3[layer_cell] = state.matrix_liquid_water_m3[layer_cell];
        macropore_solute_state.water_volume_m3[layer_cell] = state.macropore_liquid_water_m3[layer_cell];
        gas_transport_state.temperature_k[layer_cell] = state.soil_temperature_k[layer_cell];
        gas_transport_state.air_volume_m3[layer_cell] = state.matrix_air_volume_m3[layer_cell];
        if (surface_gas_parameters) |parameters| {
            const solubility = try ecosys.gas_transport.surfaceSolubilityWaterToAir(state.soil_temperature_k[layer_cell], parameters.solubility);
            for (0..ecosys.gas_transport.species_count) |species| {
                const gas_index = layer_cell * ecosys.gas_transport.species_count + species;
                const atmospheric_concentration_g_per_m3 = parameters.atmospheric_concentration_g_per_m3[species];
                gas_transport_state.gaseous_mass_g[gas_index] = atmospheric_concentration_g_per_m3 * state.matrix_air_volume_m3[layer_cell];
                gas_transport_state.dissolved_mass_g[gas_index] = if (species == @intFromEnum(ecosys.gas_transport.Species.ammonia)) 0 else atmospheric_concentration_g_per_m3 * solubility[species] * state.matrix_liquid_water_m3[layer_cell];
            }
        }
    }
    try ecosys.soil_aqueous_transport_bridge.exportConcentrations(&initial_chemistry_state, soil_recharge_concentration_mol_per_m3);
    var snow_transport_state = try ecosys.snow_solute_transport.State.init(allocator, state.cell_count, runscript.snow_layer_bottom_depth_m.len);
    defer snow_transport_state.deinit();
    const snow_atmospheric_input_g = try allocator.alloc(f64, state.cell_count * ecosys.snow_solute_transport.species_count);
    defer allocator.free(snow_atmospheric_input_g);
    @memset(snow_atmospheric_input_g, 0);
    var atmospheric_solute_input_ledger_state =
        try ecosys.atmospheric_solute_input_ledger.State.init(
            allocator,
            state.cell_count,
        );
    defer atmospheric_solute_input_ledger_state.deinit();
    const snow_surface_partitions = try allocator.alloc(ecosys.snow_solute_transport.SurfacePartition, state.cell_count);
    defer allocator.free(snow_surface_partitions);
    const snow_surface_discharge = try allocator.alloc(ecosys.snow_solute_transport.SurfaceDischarge, state.cell_count);
    defer allocator.free(snow_surface_discharge);
    const direct_surface_solute_input = try allocator.alloc(ecosys.snow_solute_transport.SurfaceDischarge, state.cell_count);
    defer allocator.free(direct_surface_solute_input);
    @memset(direct_surface_solute_input, .{});
    var transport_hydrology_state = try ecosys.transport_hydrology.State.init(allocator, config.grid_columns, config.grid_rows, config.soil_layers, snow_transport_state.layer_capacity);
    defer transport_hydrology_state.deinit();
    try transport_hydrology_state.syncStorage(&state, &snow_transport_state);
    var soil_transport_faces = try ecosys.transport_hydrology.buildSoilFacesMapped(allocator, &transport_hydrology_state, &state, lateral_connection_mode_by_cell);
    defer soil_transport_faces.deinit();
    var soil_face_geometry_state = try ecosys.soil_face_geometry.State.initMapped(
        allocator,
        &state,
        &soil_transport_faces,
        soil_solver_property_state.layer_thickness_m,
        grid_environment.horizontal_cell_width_m,
        grid_environment.vertical_cell_width_m,
    );
    defer soil_face_geometry_state.deinit();
    var soil_solute_face_parameter_state = try ecosys.soil_solute_face_parameters.State.init(allocator, soil_transport_faces.micropore_faces.len);
    defer soil_solute_face_parameter_state.deinit();
    const soil_solute_face_component_count = try std.math.mul(
        usize,
        soil_transport_faces.micropore_faces.len,
        aqueous_species_count,
    );
    const micropore_solute_face_flux_mol = try allocator.alloc(
        f64,
        soil_solute_face_component_count,
    );
    defer allocator.free(micropore_solute_face_flux_mol);
    const macropore_solute_face_flux_mol = try allocator.alloc(
        f64,
        soil_solute_face_component_count,
    );
    defer allocator.free(macropore_solute_face_flux_mol);
    @memset(micropore_solute_face_flux_mol, 0);
    @memset(macropore_solute_face_flux_mol, 0);
    var soil_organic_face_parameter_state = try ecosys.soil_organic_face_parameters.State.init(allocator, soil_transport_faces.micropore_faces.len);
    defer soil_organic_face_parameter_state.deinit();
    var soil_dissolved_gas_face_parameter_state = try ecosys.soil_dissolved_gas_face_parameters.State.init(allocator, soil_transport_faces.micropore_faces.len);
    defer soil_dissolved_gas_face_parameter_state.deinit();
    var soil_dissolved_gas_transport_state = try ecosys.soil_dissolved_gas_transport.State.init(allocator, state.layer_count);
    defer soil_dissolved_gas_transport_state.deinit();
    const soil_dissolved_gas_recharge_concentration_g_per_m3 = try allocator.alloc(f64, state.layer_count * ecosys.gas_transport.species_count);
    defer allocator.free(soil_dissolved_gas_recharge_concentration_g_per_m3);
    @memset(soil_dissolved_gas_recharge_concentration_g_per_m3, 0);
    var soil_organic_transport_state = try ecosys.soil_organic_transport.State.init(allocator, state.layer_count);
    defer soil_organic_transport_state.deinit();
    try soil_organic_transport_state.initializeFromProfile(&soil_organic_state);
    const soil_organic_recharge_concentration_g_per_m3 = try allocator.alloc(f64, state.layer_count * ecosys.soil_organic_transport.component_count);
    defer allocator.free(soil_organic_recharge_concentration_g_per_m3);
    @memset(soil_organic_recharge_concentration_g_per_m3, 0);
    var mineral_nitrogen_transport_state = try ecosys.mineral_nitrogen_transport.State.init(allocator, state.layer_count);
    defer mineral_nitrogen_transport_state.deinit();
    var mineral_nitrogen_face_parameters = try ecosys.mineral_nitrogen_transport.FaceParameters.init(allocator, soil_transport_faces.micropore_faces.len);
    defer mineral_nitrogen_face_parameters.deinit();
    const mineral_nitrogen_zone_fractions: ecosys.mineral_nitrogen_transport.ZoneFractions = .{
        .ammonium_non_band = 1 - runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
        .ammonium_band = runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
        .nitrate_non_band = 1 - runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
        .nitrate_band = runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
    };
    try mineral_nitrogen_transport_state.initializeMatrix(
        &initial_chemistry_state,
        &soil_reactive_nitrogen_state,
        state.matrix_liquid_water_m3,
        mineral_nitrogen_zone_fractions,
        runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
    );
    var surface_transport_state = try ecosys.surface_solute_routing.State.init(allocator, config.grid_columns, config.grid_rows, aqueous_species_count);
    defer surface_transport_state.deinit();
    var redistribution_ledger = try ecosys.transport_redistribution.Ledger.init(
        allocator,
        micropore_solute_state.amount_mol.len,
        gas_transport_state.gaseous_mass_g.len,
        snow_transport_state.amount_g.len,
        surface_transport_state.amount_mol.len,
        micropore_solute_state.species_count,
        ecosys.gas_transport.species_count,
        ecosys.snow_solute_transport.species_count,
        aqueous_species_count,
    );
    defer redistribution_ledger.deinit();
    var soil_thermal_context: ecosys.soil_thermal.UpdateContext = .{ .thermal = &soil_thermal_state, .grid = &state, .liquid_water_heat_capacity_mj_per_m3_k = runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_mj_per_m3_k, .ice_heat_capacity_mj_per_m3_k = runscript.soil_phase_heat_parameters.ice_heat_capacity_mj_per_m3_k };
    try runKernelAcrossSerialTilePlan(
        executor,
        &tile_plan.?,
        &soil_thermal_context,
        ecosys.soil_thermal.updateTile,
    );
    try soil_thermal_state.validateFinite();
    var soil_hourly_workspace = try ecosys.soil_hourly_workspace.State.init(allocator, state.layer_count);
    defer soil_hourly_workspace.deinit();
    try soil_hourly_workspace.refresh(&state, &soil_solver_property_state, &soil_thermal_state, &terrain_hydrology_state, runscript.soil_process_parameters);
    try soil_hourly_workspace.fillMacroporeFaceConductance(&soil_transport_faces, &soil_face_geometry_state, soil_face_geometry_state.macropore_hydraulic_conductance_m_per_h_mpa);
    var soil_heat_solver_workspace = try ecosys.soil_heat_solver.Workspace.init(
        allocator,
        state.layer_count,
        soil_transport_faces.micropore_faces.len,
        config.tile_cells,
    );
    defer soil_heat_solver_workspace.deinit();
    var atmospheric_state = try ecosys.atmospheric_forcing.State.init(allocator, state.cell_count);
    defer atmospheric_state.deinit();
    const hourly_forcing_by_cell = try allocator.alloc(ecosys.weather.HourlyForcing, state.cell_count);
    defer allocator.free(hourly_forcing_by_cell);
    const hourly_radiation_by_cell = try allocator.alloc(ecosys.atmospheric_radiation.Result, state.cell_count);
    defer allocator.free(hourly_radiation_by_cell);
    const hourly_weather_header_by_cell = try allocator.alloc(ecosys.weather.Header, state.cell_count);
    defer allocator.free(hourly_weather_header_by_cell);
    const plant_calendar_by_cell = try allocator.alloc(
        ecosys.plant_development.Calendar,
        state.cell_count,
    );
    defer allocator.free(plant_calendar_by_cell);
    const hourly_weather_reference_height_m = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(hourly_weather_reference_height_m);
    const hourly_adjusted_shortwave_mj_per_m2 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(hourly_adjusted_shortwave_mj_per_m2);
    const hourly_extraterrestrial_shortwave_mj_per_m2 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(hourly_extraterrestrial_shortwave_mj_per_m2);
    const hourly_solar_angle_sine = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(hourly_solar_angle_sine);
    const hourly_solar_azimuth_radians = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(hourly_solar_azimuth_radians);
    const irrigation_water_depth_m = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(irrigation_water_depth_m);
    @memset(irrigation_water_depth_m, 0);
    const irrigation_dissolved_mass_g_per_m2 = try allocator.alloc(f64, state.cell_count * ecosys.irrigation_management_dispatch.dissolved_species_count);
    defer allocator.free(irrigation_dissolved_mass_g_per_m2);
    @memset(irrigation_dissolved_mass_g_per_m2, 0);
    const irrigation_hydrogen_mol_per_m2 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(irrigation_hydrogen_mol_per_m2);
    @memset(irrigation_hydrogen_mol_per_m2, 0);
    var irrigation_loads = try ecosys.irrigation_layer_routing.Loads.init(
        allocator,
        state.cell_count,
        state.soil_layer_capacity,
    );
    defer irrigation_loads.deinit();
    const automatic_irrigation_depth_m_by_cell_hour = try allocator.alloc(f64, state.cell_count * 24);
    defer allocator.free(automatic_irrigation_depth_m_by_cell_hour);
    @memset(automatic_irrigation_depth_m_by_cell_hour, 0);
    const minimum_canopy_water_potential_mpa_by_cell = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(minimum_canopy_water_potential_mpa_by_cell);
    @memset(minimum_canopy_water_potential_mpa_by_cell, 0);
    const soil_wilting_point_fraction = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(soil_wilting_point_fraction);
    for (soil_solver_property_state.retention_curve, soil_wilting_point_fraction) |curve, *wilting| wilting.* = curve.curve.wilting_point_fraction;
    var canopy_radiation_state = try ecosys.canopy_radiation.State.init(allocator, state.cell_count);
    defer canopy_radiation_state.deinit();
    var canopy_optics_state: ?ecosys.canopy_optics.State = null;
    defer if (canopy_optics_state) |*optics_state| optics_state.deinit();
    if (plant_assignments) |assignments| canopy_optics_state = try ecosys.canopy_optics.State.initMapped(
        allocator,
        state.cell_count,
        config.plant_populations,
        assignments,
        plant_unit_by_cell.?,
        plant_catalog,
    );
    const canopy_cell_area_m2 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(canopy_cell_area_m2);
    const surface_pond_minimum_heat_capacity_mj_per_k = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_pond_minimum_heat_capacity_mj_per_k);
    var total_grid_area_m2: f64 = 0;
    for (canopy_cell_area_m2, 0..) |*area_m2, cell| {
        area_m2.* = grid_environment.horizontal_cell_width_m[cell] *
            grid_environment.vertical_cell_width_m[cell];
        if (!std.math.isFinite(area_m2.*) or area_m2.* <= 0) return error.InvalidCanopyCellArea;
        surface_pond_minimum_heat_capacity_mj_per_k[cell] = runscript.surface_pond_activation_heat_capacity_mj_per_m2_k * area_m2.*;
        total_grid_area_m2 += area_m2.*;
    }
    if (!std.math.isFinite(total_grid_area_m2) or total_grid_area_m2 <= 0) return error.InvalidOutputGridArea;
    const initial_ground_air_temperature_k = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(initial_ground_air_temperature_k);
    @memset(initial_ground_air_temperature_k, site.mean_annual_air_temperature_c + 273.15);
    const initial_ground_air_vapor_fraction = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(initial_ground_air_vapor_fraction);
    const initial_ground_air_temperature = site.mean_annual_air_temperature_c + 273.15;
    const initial_ground_air_saturation = runscript.ground_air_parameters.saturation_vapor_prefactor_k / initial_ground_air_temperature * runscript.ground_air_parameters.saturation_relative_humidity * @exp(runscript.ground_air_parameters.saturation_temperature_k * (runscript.ground_air_parameters.saturation_reference_inverse_temperature_per_k - 1 / initial_ground_air_temperature));
    @memset(initial_ground_air_vapor_fraction, initial_ground_air_saturation);
    var surface_aerodynamic_state = try ecosys.surface_aerodynamics.State.init(allocator, state.cell_count, runscript.surface_aerodynamic_parameters.soil_roughness_height_m);
    defer surface_aerodynamic_state.deinit();
    var ground_air_state = try ecosys.ground_air_exchange.State.init(allocator, initial_ground_air_temperature_k, initial_ground_air_vapor_fraction, canopy_cell_area_m2, runscript.ground_air_parameters.minimum_air_column_height_m, runscript.ground_air_parameters);
    defer ground_air_state.deinit();
    const surface_total_canopy_area_m2 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_total_canopy_area_m2);
    @memset(surface_total_canopy_area_m2, 0);
    const surface_canopy_height_m = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_canopy_height_m);
    @memset(surface_canopy_height_m, 0);
    const ground_air_vapor_pressure_kpa = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_vapor_pressure_kpa);
    @memset(ground_air_vapor_pressure_kpa, 0);
    const atmospheric_vapor_fraction = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(atmospheric_vapor_fraction);
    @memset(atmospheric_vapor_fraction, 0);
    const ground_air_canopy_resistance_h_per_m = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_canopy_resistance_h_per_m);
    @memset(ground_air_canopy_resistance_h_per_m, 0);
    const ground_air_sensible_source_mj_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_sensible_source_mj_per_h);
    @memset(ground_air_sensible_source_mj_per_h, 0);
    const runtime_plant_count = try std.math.mul(usize, state.cell_count, config.plant_populations);
    const grazing_average_shoot_carbon_g_c = try allocator.alloc(f64, config.plant_populations);
    defer allocator.free(grazing_average_shoot_carbon_g_c);
    @memset(grazing_average_shoot_carbon_g_c, 0);
    const grazing_active_cell_count = try allocator.alloc(usize, config.plant_populations);
    defer allocator.free(grazing_active_cell_count);
    @memset(grazing_active_cell_count, 0);
    const grazing_total_shoot_carbon_g_c_by_cell = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(grazing_total_shoot_carbon_g_c_by_cell);
    const grazing_total_grazer_biomass_by_cell = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(grazing_total_grazer_biomass_by_cell);
    var plant_daily_flux_ledger = try ecosys.plant_daily_flux_ledger.State.init(allocator, runtime_plant_count);
    defer plant_daily_flux_ledger.deinit();
    const delayed_live_canopy_combustion_heat_mj = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(delayed_live_canopy_combustion_heat_mj);
    @memset(delayed_live_canopy_combustion_heat_mj, 0);
    const delayed_standing_dead_combustion_heat_mj = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(delayed_standing_dead_combustion_heat_mj);
    @memset(delayed_standing_dead_combustion_heat_mj, 0);
    const delayed_subsurface_combustion_heat_mj = try allocator.alloc(f64, state.layer_count);
    defer allocator.free(delayed_subsurface_combustion_heat_mj);
    @memset(delayed_subsurface_combustion_heat_mj, 0);
    const delayed_surface_combustion_heat_mj = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(delayed_surface_combustion_heat_mj);
    @memset(delayed_surface_combustion_heat_mj, 0);
    const surface_combustion_heat_mj_per_m2 = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(surface_combustion_heat_mj_per_m2);
    @memset(surface_combustion_heat_mj_per_m2, 0);
    const ground_air_vapor_source_m3_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_vapor_source_m3_per_h);
    @memset(ground_air_vapor_source_m3_per_h, 0);
    const ground_air_surface_sensible_conductance_mj_per_h_k = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_surface_sensible_conductance_mj_per_h_k);
    @memset(ground_air_surface_sensible_conductance_mj_per_h_k, 0);
    const ground_air_surface_vapor_conductance_m3_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_surface_vapor_conductance_m3_per_h);
    @memset(ground_air_surface_vapor_conductance_m3_per_h, 0);
    const ground_air_surface_vapor_fraction = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(ground_air_surface_vapor_fraction);
    @memset(ground_air_surface_vapor_fraction, 0);
    var canopy_geometry = try ecosys.canopy_geometry.Geometry.init(allocator, runscript.canopy_discretization);
    defer canopy_geometry.deinit();
    var canopy_structure_state: ?ecosys.canopy_structure.State = null;
    defer if (canopy_structure_state) |*structure_state| structure_state.deinit();
    if (plant_assignments) |assignments| {
        canopy_structure_state = try ecosys.canopy_structure.State.initMapped(
            allocator,
            state.cell_count,
            config.plant_populations,
            canopy_geometry.leaf_inclination_sine.len,
            assignments,
            plant_unit_by_cell.?,
            plant_catalog,
        );
        var structure_context: ecosys.canopy_structure.ApplyContext = .{ .structure = &canopy_structure_state.?, .plants = &plant_state };
        try runKernelAcrossSerialTilePlan(
            executor,
            &tile_plan.?,
            &structure_context,
            ecosys.canopy_structure.applyLeafAreaTile,
        );
    }
    var detailed_canopy_state: ?ecosys.canopy_photosynthesis.State = null;
    defer if (detailed_canopy_state) |*canopy_state| canopy_state.deinit();
    var canopy_carbon_exchange_state: ?ecosys.canopy_carbon_exchange.State = null;
    defer if (canopy_carbon_exchange_state) |*ledger| ledger.deinit();
    var canopy_layer_distribution_state: ?ecosys.canopy_layer_distribution.State = null;
    defer if (canopy_layer_distribution_state) |*layer_state| layer_state.deinit();
    var canopy_precipitation_retention_state: ?ecosys.canopy_precipitation_retention.State = null;
    defer if (canopy_precipitation_retention_state) |*retention_state| retention_state.deinit();
    var surface_precipitation_state = try ecosys.surface_precipitation.RuntimeState.init(allocator, state.cell_count);
    defer surface_precipitation_state.deinit();
    if (surface_gas_parameters) |parameters| {
        for (0..state.cell_count) |cell| surface_precipitation_state.litter_water_m3[cell] = parameters.initial_litter_water_m3_per_g_c * try surface_organic_state.totalCarbon_g_c(cell);
        try surface_litter_chemistry_state.bindMineralReferenceWater(
            surface_precipitation_state.litter_water_m3,
        );
        var geometry_context: ecosys.surface_litter_geometry_step.ApplyContext = .{ .result = &surface_litter_geometry_state, .surface_organic = &surface_organic_state, .water_m3 = surface_precipitation_state.litter_water_m3, .ice_m3 = surface_litter_ice_m3, .charcoal_carbon_g_c = surface_charcoal_carbon_g_c, .parameters = parameters.litter_geometry };
        try runKernelAcrossSerialTilePlan(
            executor,
            &tile_plan.?,
            &geometry_context,
            ecosys.surface_litter_geometry_step.applyTile,
        );
        @memcpy(surface_precipitation_state.litter_water_capacity_m3, surface_litter_geometry_state.water_retention_capacity_m3);
        for (0..state.cell_count) |cell| try ecosys.gas_transport.initializeSurfaceCell(&litter_gas_transport_state, cell, surface_litter_geometry_state.air_volume_m3[cell], surface_precipitation_state.litter_water_m3[cell], state.surface_temperature_k[cell], parameters.atmospheric_concentration_g_per_m3, parameters.solubility);
    } else if (organic_parameters) |*parameters| {
        for (0..state.cell_count) |cell| {
            const litter_carbon_g_c = try surface_organic_state.totalCarbon_g_c(cell);
            surface_precipitation_state.litter_water_capacity_m3[cell] = parameters.surface_litter_water_capacity_m3_per_g_c * litter_carbon_g_c;
            if (!std.math.isFinite(surface_precipitation_state.litter_water_capacity_m3[cell])) return error.NonFiniteSurfaceLitterWaterCapacity;
        }
    }
    var plant_root_state: ?ecosys.plant_root_system.State = null;
    defer if (plant_root_state) |*root_state| root_state.deinit();
    var plant_root_nutrient_workspace: ?ecosys.plant_root_nutrient_uptake.GridWorkspace = null;
    defer if (plant_root_nutrient_workspace) |*workspace| workspace.deinit();
    var plant_root_salt_workspace: ?ecosys.plant_root_salt_exchange.GridWorkspace = null;
    defer if (plant_root_salt_workspace) |*workspace| workspace.deinit();
    var plant_root_exudation_workspace: ?ecosys.plant_root_exudation.GridWorkspace = null;
    defer if (plant_root_exudation_workspace) |*workspace| workspace.deinit();
    var plant_root_metabolism_workspace: ?ecosys.plant_root_metabolism.GridWorkspace = null;
    defer if (plant_root_metabolism_workspace) |*workspace| workspace.deinit();
    var plant_storage_remobilization_workspace: ?ecosys.plant_storage_remobilization.Workspace = null;
    defer if (plant_storage_remobilization_workspace) |*workspace| workspace.deinit();
    var plant_litter_partition_state: ?ecosys.plant_litter_partition.State = null;
    defer if (plant_litter_partition_state) |*partition_state| partition_state.deinit();
    var plant_phenology_state: ?ecosys.plant_phenology.State = null;
    defer if (plant_phenology_state) |*phenology_state| phenology_state.deinit();
    const phenology_root_oxygen_fraction = try allocator.alloc(f64, try std.math.mul(usize, state.cell_count, config.plant_populations));
    defer allocator.free(phenology_root_oxygen_fraction);
    @memset(phenology_root_oxygen_fraction, 0);
    var branch_development_state: ?ecosys.plant_phenology.BranchDevelopmentState = null;
    defer if (branch_development_state) |*development_state| development_state.deinit();
    var plant_growth_stage_state: ?ecosys.plant_growth_stages.State = null;
    defer if (plant_growth_stage_state) |*growth_stage_state| growth_stage_state.deinit();
    var plant_dormancy_state: ?ecosys.plant_dormancy.RuntimeState = null;
    defer if (plant_dormancy_state) |*dormancy_state| dormancy_state.deinit();
    const root_nutrient_traits = try allocator.alloc([3]ecosys.plant_traits.NutrientUptake, runtime_plant_count);
    defer allocator.free(root_nutrient_traits);
    const canopy_biochemistry_parameters = try allocator.alloc(ecosys.canopy_biochemistry.Parameters, runtime_plant_count);
    defer allocator.free(canopy_biochemistry_parameters);
    const shoot_growth_plant_parameters = try allocator.alloc(ecosys.shoot_growth_runtime.PlantParameters, runtime_plant_count);
    defer allocator.free(shoot_growth_plant_parameters);
    @memset(shoot_growth_plant_parameters, ecosys.shoot_growth_runtime.inactivePlantParameters());
    @memset(canopy_biochemistry_parameters, .{
        .pathway = .c3,
        .growth_habit = 0,
        .phenology_type = 0,
        .aboveground_turnover_type = 0,
        .rubisco_carboxylation_umol_per_g_protein_s = 1,
        .rubisco_oxygenation_umol_per_g_protein_s = 0,
        .pep_carboxylation_umol_per_g_protein_s = 0,
        .rubisco_co2_half_saturation_umol_per_l = 1,
        .rubisco_o2_half_saturation_umol_per_l = 1,
        .pep_co2_half_saturation_umol_per_l = 1,
        .rubisco_leaf_protein_fraction = 0,
        .pep_leaf_protein_fraction = 0,
        .chlorophyll_electron_transport_umol_per_g_protein_s = 0,
        .c3_chlorophyll_leaf_protein_fraction = 0,
        .c4_chlorophyll_leaf_protein_fraction = 0,
        .intercellular_to_atmospheric_co2_ratio = 0.7,
    });
    const inactive_nutrient_trait: ecosys.plant_traits.NutrientUptake = .{ .maximum_rate_g_per_m2_h = 0, .half_saturation_umol_per_l = 0, .minimum_concentration_umol_per_l = 0 };
    @memset(root_nutrient_traits, .{inactive_nutrient_trait} ** 3);
    const root_nutrient_feedback_enabled = try allocator.alloc(bool, runtime_plant_count);
    defer allocator.free(root_nutrient_feedback_enabled);
    @memset(root_nutrient_feedback_enabled, false);
    const root_metabolism_plant_parameters = try allocator.alloc(ecosys.plant_root_metabolism.RuntimePlantParameters, runtime_plant_count);
    defer allocator.free(root_metabolism_plant_parameters);
    const root_biological_domain_count_by_plant = try allocator.alloc(u8, runtime_plant_count);
    defer allocator.free(root_biological_domain_count_by_plant);
    @memset(root_biological_domain_count_by_plant, 1);
    @memset(root_metabolism_plant_parameters, .{
        .root_profile_type = 0,
        .mycorrhizal_type = 0,
        .growth_habit = 0,
        .leaf_phenology_type = 0,
        .root_growth_yield_g_c_per_g_c = 1,
        .root_nitrogen_to_carbon_g_n_per_g_c = 1,
        .root_phosphorus_to_carbon_g_p_per_g_c = 1,
        .stalk_nitrogen_to_carbon_g_n_per_g_c = 1,
        .stalk_phosphorus_to_carbon_g_p_per_g_c = 1,
        .primary_root_radius_m = 1,
        .secondary_root_radius_m = 1,
        .primary_specific_length_m_per_g_c = 1,
        .secondary_specific_length_m_per_g_c = 1,
        .secondary_root_branching_per_m = 1,
        .shoot_root_equilibration_fraction_per_h = 0,
    });
    const checkpoint_species_names = try allocator.alloc([]const u8, runtime_plant_count);
    defer allocator.free(checkpoint_species_names);
    @memset(checkpoint_species_names, "inactive");
    const checkpoint_metadata_cells = try allocator.alloc(ecosys.plant_checkpoint_metadata.CellView, state.cell_count);
    defer allocator.free(checkpoint_metadata_cells);
    @memset(checkpoint_metadata_cells, .{ .species_names = &.{}, .species_alive = &.{} });
    var plant_topology_controls = try ecosys.plant_topology.Controls.init(allocator, runtime_plant_count);
    defer plant_topology_controls.deinit();
    var plant_reproduction_controls = try ecosys.plant_reproduction.Controls.init(allocator, runtime_plant_count);
    defer plant_reproduction_controls.deinit();
    var canopy_layer_controls = try ecosys.canopy_layer_distribution.Controls.init(allocator, runtime_plant_count);
    defer canopy_layer_controls.deinit();
    const development_species_parameters = try allocator.alloc(ecosys.plant_growth_stages.SpeciesParameters, runtime_plant_count);
    defer allocator.free(development_species_parameters);
    @memset(development_species_parameters, .{ .growth_habit = .annual, .phenology_type = .evergreen, .photoperiod_type = .insensitive, .maturity_group_node_count = 1, .branch_floral_node_requirement = 1, .critical_photoperiod_h = 12, .photoperiod_sensitivity_h = 0, .determinate = true, .vegetative_stage_duration = runscript.phenology_parameters.vegetative_stage_duration, .reproductive_stage_duration = runscript.phenology_parameters.reproductive_stage_duration });
    const development_dormancy_parameters = try allocator.alloc(ecosys.plant_dormancy.Parameters, runtime_plant_count);
    defer allocator.free(development_dormancy_parameters);
    @memset(development_dormancy_parameters, .{ .required_leafout_h = 0, .required_leafoff_h = 0, .leafout_temperature_threshold_c = 0, .leafoff_temperature_threshold_c = 0, .chilling_temperature_c = 0, .drought_leafout_total_water_potential_mpa = runscript.phenology_parameters.drought_leafout_total_water_potential_mpa, .combined_leafout_turgor_potential_mpa = runscript.phenology_parameters.minimum_turgor_potential_mpa, .leafoff_total_water_potential_mpa = runscript.phenology_parameters.vascular_leafoff_total_water_potential_mpa, .maximum_photoperiod_counter_h = runscript.phenology_parameters.maximum_photoperiod_counter_h, .evergreen_leafoff_remobilization_start_fraction = runscript.root_metabolism_parameters.evergreen_leafoff_remobilization_start_fraction, .deciduous_leafoff_remobilization_start_fraction = runscript.root_metabolism_parameters.deciduous_leafoff_remobilization_start_fraction, .full_senescence_duration_h = runscript.root_metabolism_parameters.full_senescence_duration_h });
    const development_planting_dates = try allocator.alloc(ecosys.plant_management.PackedDate, runtime_plant_count);
    defer allocator.free(development_planting_dates);
    @memset(development_planting_dates, .{ .day = 1, .month = 1, .year = 9999 });
    const development_planting_day_of_year = try allocator.alloc(u16, runtime_plant_count);
    defer allocator.free(development_planting_day_of_year);
    @memset(development_planting_day_of_year, 1);
    const development_planting_year = try allocator.alloc(i32, runtime_plant_count);
    defer allocator.free(development_planting_year);
    @memset(development_planting_year, 0);
    const development_canopy_height_m = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(development_canopy_height_m);
    @memset(development_canopy_height_m, 0);
    const development_surface_water_potential_mpa = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(development_surface_water_potential_mpa);
    @memset(development_surface_water_potential_mpa, 0);
    const development_seed_layer_water_potential_mpa = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(development_seed_layer_water_potential_mpa);
    @memset(development_seed_layer_water_potential_mpa, 0);
    const newly_activated_plants = try allocator.alloc(bool, runtime_plant_count);
    defer allocator.free(newly_activated_plants);
    @memset(newly_activated_plants, false);
    const reseed_population_per_m2_by_plant = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(reseed_population_per_m2_by_plant);
    @memset(reseed_population_per_m2_by_plant, 0);
    const replant_layer_bottom_depth_m = try allocator.alloc(f64, config.soil_layers);
    defer allocator.free(replant_layer_bottom_depth_m);
    var management_schedule_map: ?ecosys.plant_management_dispatch.ScheduleMap = null;
    defer if (management_schedule_map) |*schedule_map| schedule_map.deinit();
    var harvest_science_by_plant: ?[]ecosys.plant_harvest_runtime.ScienceParameters = null;
    defer if (harvest_science_by_plant) |science| allocator.free(science);
    var harvest_products_by_plant: ?[]ecosys.plant_harvest_runtime.ProductLedger = null;
    defer if (harvest_products_by_plant) |products| allocator.free(products);
    var harvest_context: ?ecosys.plant_harvest_runtime.Context = null;
    if (plant_assignments) |assignments| {
        const plant_count = try std.math.mul(usize, state.cell_count, config.plant_populations);
        plant_root_state = try ecosys.plant_root_system.State.init(allocator, plant_count, config.soil_layers, runscript.root_axes_per_plant);
        plant_root_nutrient_workspace = try ecosys.plant_root_nutrient_uptake.GridWorkspace.init(
            allocator,
            state.cell_count,
            try std.math.mul(usize, config.plant_populations, ecosys.plant_root_system.biological_domain_count),
        );
        plant_root_salt_workspace = try ecosys.plant_root_salt_exchange.GridWorkspace.init(
            allocator,
            state.cell_count,
            try std.math.mul(usize, config.plant_populations, ecosys.plant_root_system.biological_domain_count),
        );
        plant_root_exudation_workspace = try ecosys.plant_root_exudation.GridWorkspace.init(
            allocator,
            state.cell_count,
            try std.math.mul(usize, config.plant_populations, ecosys.plant_root_system.biological_domain_count),
        );
        plant_root_metabolism_workspace = try ecosys.plant_root_metabolism.GridWorkspace.init(allocator, state.cell_count, runscript.root_axes_per_plant, config.soil_layers);
        plant_storage_remobilization_workspace = try ecosys.plant_storage_remobilization.Workspace.init(allocator, plant_count, config.soil_layers);
        plant_litter_partition_state = try ecosys.plant_litter_partition.State.init(allocator, plant_count);
        plant_phenology_state = try ecosys.plant_phenology.State.init(allocator, state.cell_count, config.plant_populations);
        const branch_count_by_plant = try allocator.alloc(usize, plant_count);
        defer allocator.free(branch_count_by_plant);
        @memset(branch_count_by_plant, 0);
        for (plant_unit_by_cell.?, 0..) |unit_index, cell| {
            const active_species_count = assignments.units[unit_index].species.len;
            for (0..active_species_count) |species| branch_count_by_plant[cell * config.plant_populations + species] = 1;
        }
        var branch_count: usize = 0;
        for (branch_count_by_plant) |count| branch_count = try std.math.add(usize, branch_count, count);
        if (branch_count > 0) {
            plant_growth_stage_state = try ecosys.plant_growth_stages.State.init(allocator, branch_count_by_plant);
            plant_dormancy_state = try ecosys.plant_dormancy.RuntimeState.init(allocator, branch_count);
            const node_count_by_branch = try allocator.alloc(usize, branch_count);
            defer allocator.free(node_count_by_branch);
            @memset(node_count_by_branch, 1);
            const sample_count_by_node = try allocator.alloc(usize, branch_count);
            defer allocator.free(sample_count_by_node);
            const angular_sample_count = try std.math.mul(usize, canopy_geometry.leaf_inclination_sine.len, canopy_geometry.leaf_azimuth_radians.len);
            const radiation_sample_count_per_node = try std.math.mul(usize, runscript.canopy_layer_count, angular_sample_count);
            @memset(sample_count_by_node, radiation_sample_count_per_node);
            detailed_canopy_state = try ecosys.canopy_photosynthesis.State.init(allocator, state.cell_count, config.plant_populations, branch_count_by_plant, node_count_by_branch, sample_count_by_node);
            canopy_carbon_exchange_state = try ecosys.canopy_carbon_exchange.State.init(allocator, branch_count);
            canopy_layer_distribution_state = try ecosys.canopy_layer_distribution.State.init(allocator, state.cell_count, config.plant_populations, runscript.canopy_layer_count, canopy_geometry.leaf_inclination_sine.len, canopy_geometry.leaf_azimuth_radians.len, &detailed_canopy_state.?);
            canopy_precipitation_retention_state = try ecosys.canopy_precipitation_retention.State.init(allocator, state.cell_count, config.plant_populations);
            for (plant_unit_by_cell.?, 0..) |unit_index, cell| {
                const cell_area_m2 = grid_environment.horizontal_cell_width_m[cell] *
                    grid_environment.vertical_cell_width_m[cell];
                const profile = soil_catalog.entries.items[catalog_index_by_cell[cell]].profile;
                const layer_bottom_depth_m = try allocator.alloc(f64, profile.total_layer_count);
                defer allocator.free(layer_bottom_depth_m);
                var cumulative_depth_m: f64 = 0;
                for (layer_bottom_depth_m, 0..) |*bottom, layer| {
                    cumulative_depth_m += soil_solver_property_state.layer_thickness_m[cell * config.soil_layers + layer];
                    bottom.* = cumulative_depth_m;
                }
                for (assignments.units[unit_index].species, 0..) |assignment, species| {
                    if (ecosys.delimited_input.isNo(assignment.management_file)) continue;
                    const trait_index = plant_catalog.find(assignment.species_file) orelse return error.MissingPlantTraitProfile;
                    const management_index = plant_management_catalog.find(assignment.management_file) orelse return error.MissingPlantManagementSchedule;
                    const traits = plant_catalog.entries.items[trait_index].traits;
                    const planting = plant_management_catalog.entries.items[management_index].schedule.planting;
                    const planting_layer = try ecosys.plant_initialization.plantingLayer(planting.seed_depth_m, layer_bottom_depth_m);
                    const storage = try ecosys.plant_initialization.seedStorage(
                        traits.morphology.seed_mass_at_planting_g,
                        planting.population_per_m2,
                        cell_area_m2,
                        traits.organ_nitrogen_to_carbon_ratio.grain,
                        traits.organ_phosphorus_to_carbon_ratio.grain,
                    );
                    const plant = try detailed_canopy_state.?.plantIndex(cell, species);
                    reseed_population_per_m2_by_plant[plant] = planting.population_per_m2;
                    root_nutrient_traits[plant] = .{ traits.ammonium_uptake, traits.nitrate_uptake, traits.phosphate_uptake };
                    canopy_biochemistry_parameters[plant] = .{
                        .pathway = if (traits.functional_type.photosynthesis_pathway == 3) .c3 else .c4,
                        .growth_habit = traits.functional_type.growth_habit,
                        .phenology_type = traits.functional_type.leaf_phenology_type,
                        .aboveground_turnover_type = traits.functional_type.aboveground_turnover_type,
                        .rubisco_carboxylation_umol_per_g_protein_s = traits.photosynthesis.rubisco_carboxylase_umol_c_per_g_enzyme_s,
                        .rubisco_oxygenation_umol_per_g_protein_s = traits.photosynthesis.rubisco_oxygenase_umol_o_per_g_enzyme_s,
                        .pep_carboxylation_umol_per_g_protein_s = traits.photosynthesis.pep_carboxylase_umol_per_g_enzyme_s,
                        .rubisco_co2_half_saturation_umol_per_l = traits.photosynthesis.co2_half_saturation_umol_per_l,
                        .rubisco_o2_half_saturation_umol_per_l = traits.photosynthesis.oxygen_half_saturation_umol_per_l,
                        .pep_co2_half_saturation_umol_per_l = traits.photosynthesis.pep_co2_half_saturation_umol_per_l,
                        .rubisco_leaf_protein_fraction = traits.photosynthesis.rubisco_leaf_protein_fraction,
                        .pep_leaf_protein_fraction = traits.photosynthesis.pep_carboxylase_leaf_protein_fraction,
                        .chlorophyll_electron_transport_umol_per_g_protein_s = traits.photosynthesis.chlorophyll_electron_activity_umol_per_g_enzyme_s,
                        .c3_chlorophyll_leaf_protein_fraction = traits.photosynthesis.mesophyll_chlorophyll_leaf_protein_fraction,
                        .c4_chlorophyll_leaf_protein_fraction = traits.photosynthesis.c4_mesophyll_chlorophyll_leaf_protein_fraction,
                        .intercellular_to_atmospheric_co2_ratio = traits.photosynthesis.intercellular_to_atmospheric_co2_ratio,
                    };
                    try canopy_biochemistry_parameters[plant].validate();
                    shoot_growth_plant_parameters[plant] = ecosys.shoot_growth_runtime.parametersFromTraits(traits);
                    root_nutrient_feedback_enabled[plant] = traits.functional_type.growth_habit != 0;
                    const root_geometry = try ecosys.plant_initialization.rootGeometry(traits.roots.root_porosity_fraction, traits.roots.primary_root_radius_m, traits.roots.secondary_root_radius_m, runscript.plant_geometry_parameters);
                    root_metabolism_plant_parameters[plant] = .{
                        .root_profile_type = traits.functional_type.root_profile_type,
                        .mycorrhizal_type = traits.functional_type.mycorrhizal_type,
                        .growth_habit = if (traits.functional_type.growth_habit == 0) 0 else 1,
                        .leaf_phenology_type = traits.functional_type.leaf_phenology_type,
                        .root_growth_yield_g_c_per_g_c = traits.organ_growth_yield_g_c_per_g_c.root,
                        .root_nitrogen_to_carbon_g_n_per_g_c = traits.organ_nitrogen_to_carbon_ratio.root,
                        .root_phosphorus_to_carbon_g_p_per_g_c = traits.organ_phosphorus_to_carbon_ratio.root,
                        .stalk_nitrogen_to_carbon_g_n_per_g_c = traits.organ_nitrogen_to_carbon_ratio.stalk,
                        .stalk_phosphorus_to_carbon_g_p_per_g_c = traits.organ_phosphorus_to_carbon_ratio.stalk,
                        .primary_root_radius_m = traits.roots.primary_root_radius_m,
                        .secondary_root_radius_m = traits.roots.secondary_root_radius_m,
                        .primary_specific_length_m_per_g_c = root_geometry.primary_specific_length_m_per_g_c,
                        .secondary_specific_length_m_per_g_c = root_geometry.secondary_specific_length_m_per_g_c,
                        .secondary_root_branching_per_m = traits.roots.secondary_root_branching_per_m,
                        .shoot_root_equilibration_fraction_per_h = traits.roots.shoot_root_carbon_equilibration_fraction_per_h,
                    };
                    try root_metabolism_plant_parameters[plant].validate();
                    root_biological_domain_count_by_plant[plant] = @intCast(root_metabolism_plant_parameters[plant].biologicalDomainCount());
                    const adjusted_phenology = try ecosys.plant_initialization.adjustedPhenology(traits, runscript.phenology_initialization_parameters);
                    try plant_topology_controls.setPlant(plant, traits.phenology.branching_nonstructural_carbon_fraction, traits.roots.branching_nonstructural_carbon_fraction, traits.functional_type.growth_habit, traits.functional_type.leaf_phenology_type, adjusted_phenology.floral_initiation_node_count_after_seed);
                    try plant_reproduction_controls.setPlant(plant, traits.morphology.maximum_seed_sites, traits.morphology.maximum_seeds_per_site, traits.morphology.maximum_seed_mass_g, traits.phenology.chilling_temperature_c, traits.water_relations.stomatal_turgor_shape, traits.functional_type.root_profile_type == 0);
                    try canopy_layer_controls.setPlant(plant, traits.phenology.leaf_length_to_width_ratio, runscript.stalk_volume_m3_per_g_c, traits.functional_type.aboveground_turnover_type, traits.functional_type.root_profile_type, traits.functional_type.growth_habit == 0, @sin(traits.morphology.stem_angle_degrees * std.math.pi / 180.0));
                    const thermal_acclimation = try ecosys.plant_initialization.thermalAcclimation(
                        assignment.species_file,
                        traits.functional_type.photosynthesis_pathway,
                        traits.functional_type.thermal_adaptation_zone,
                        runscript.thermal_acclimation_parameters,
                    );
                    development_species_parameters[plant] = .{
                        .growth_habit = ecosys.plant_growth_stages.growthHabitFromReadq(traits.functional_type.growth_habit),
                        .phenology_type = try ecosys.plant_growth_stages.phenologyTypeFromReadq(traits.functional_type.leaf_phenology_type),
                        .photoperiod_type = try ecosys.plant_growth_stages.photoperiodTypeFromReadq(traits.functional_type.photoperiod_type),
                        .maturity_group_node_count = adjusted_phenology.floral_initiation_node_count_after_seed,
                        .branch_floral_node_requirement = adjusted_phenology.floral_initiation_node_count_after_seed,
                        .critical_photoperiod_h = traits.phenology.critical_photoperiod_h,
                        .photoperiod_sensitivity_h = traits.phenology.floral_induction_photoperiod_difference_h,
                        .determinate = traits.functional_type.determinacy_type == 0,
                        .vegetative_stage_duration = runscript.phenology_parameters.vegetative_stage_duration,
                        .reproductive_stage_duration = runscript.phenology_parameters.reproductive_stage_duration,
                    };
                    development_dormancy_parameters[plant] = .{
                        .required_leafout_h = traits.phenology.spring_leafout_requirement_h,
                        .required_leafoff_h = traits.phenology.autumn_leafoff_requirement_h,
                        .leafout_temperature_threshold_c = thermal_acclimation.leafout_threshold_c,
                        .leafoff_temperature_threshold_c = thermal_acclimation.leafoff_threshold_c,
                        .chilling_temperature_c = traits.phenology.chilling_temperature_c,
                        .drought_leafout_total_water_potential_mpa = runscript.phenology_parameters.drought_leafout_total_water_potential_mpa,
                        .combined_leafout_turgor_potential_mpa = runscript.phenology_parameters.minimum_turgor_potential_mpa,
                        .leafoff_total_water_potential_mpa = if (traits.functional_type.root_profile_type == 0) runscript.phenology_parameters.nonvascular_leafoff_total_water_potential_mpa else runscript.phenology_parameters.vascular_leafoff_total_water_potential_mpa,
                        .maximum_photoperiod_counter_h = runscript.phenology_parameters.maximum_photoperiod_counter_h,
                        .evergreen_leafoff_remobilization_start_fraction = runscript.root_metabolism_parameters.evergreen_leafoff_remobilization_start_fraction,
                        .deciduous_leafoff_remobilization_start_fraction = runscript.root_metabolism_parameters.deciduous_leafoff_remobilization_start_fraction,
                        .full_senescence_duration_h = runscript.root_metabolism_parameters.full_senescence_duration_h,
                    };
                    development_planting_dates[plant] = planting.date;
                    plant_phenology_state.?.active[plant] = false;
                    plant_phenology_state.?.annual_growth_habit[plant] = traits.functional_type.growth_habit == 0;
                    plant_phenology_state.?.node_initiation_rate_at_25c_per_h[plant] = adjusted_phenology.node_initiation_per_h;
                    plant_phenology_state.?.leaf_appearance_rate_at_25c_per_h[plant] = adjusted_phenology.leaf_appearance_per_h;
                    plant_phenology_state.?.initiated_node_count[plant] = adjusted_phenology.seed_initial_node_count;
                    plant_phenology_state.?.appeared_leaf_count[plant] = adjusted_phenology.seed_initial_node_count;
                    const growth_stage_range = try plant_growth_stage_state.?.branchRange(plant);
                    if (growth_stage_range.first < growth_stage_range.end) {
                        plant_growth_stage_state.?.branches[growth_stage_range.first].initiated_node_count = adjusted_phenology.seed_initial_node_count;
                        plant_growth_stage_state.?.branches[growth_stage_range.first].appeared_leaf_count = adjusted_phenology.seed_initial_node_count;
                    }
                    try plant_root_state.?.initializePlant(plant, traits, planting_layer, planting.seed_depth_m, runscript.root_initialization_parameters);
                    const seed_geometry = try ecosys.plant_initialization.seedGeometry(traits.morphology.seed_mass_at_planting_g, runscript.plant_geometry_parameters);
                    try plant_root_state.?.setSeedGeometry(plant, seed_geometry.volume_m3, seed_geometry.length_m, seed_geometry.surface_area_m2);
                    try plant_litter_partition_state.?.initializePlant(plant, traits, runscript.standing_dead_partition_parameters);
                    try ecosys.plant_initialization.initializeShootControls(
                        &detailed_canopy_state.?,
                        plant,
                        planting.population_per_m2,
                        cell_area_m2,
                        traits.water_relations.cuticular_resistance_s_per_m,
                        traits.functional_type.photosynthesis_pathway,
                        runscript.shoot_control_parameters,
                    );
                    try ecosys.plant_initialization.initializeThermalAcclimation(
                        &detailed_canopy_state.?,
                        plant,
                        thermal_acclimation,
                    );
                    try ecosys.plant_initialization.initializeSeedStorage(&detailed_canopy_state.?, plant, storage);
                    const standing_dead = try ecosys.plant_initialization.standingDeadStorage(
                        traits.morphology.standing_dead_carbon_g_per_m2,
                        cell_area_m2,
                        traits.organ_nitrogen_to_carbon_ratio.stalk,
                        traits.organ_phosphorus_to_carbon_ratio.stalk,
                        runscript.standing_dead_partition_parameters,
                    );
                    try ecosys.plant_initialization.initializeStandingDeadStorage(&detailed_canopy_state.?, plant, standing_dead);
                    try ecosys.plant_initialization.initializePlantHeatAndWater(
                        &plant_state,
                        &detailed_canopy_state.?,
                        plant,
                        site.mean_annual_air_temperature_c,
                        traits.water_relations.osmotic_potential_mpa,
                        runscript.plant_heat_water_parameters,
                    );
                }
            }
            branch_development_state = try ecosys.plant_phenology.BranchDevelopmentState.init(allocator, branch_count);
            for (plant_unit_by_cell.?, 0..) |unit_index, cell| {
                for (assignments.units[unit_index].species, 0..) |assignment, species| {
                    if (ecosys.delimited_input.isNo(assignment.management_file)) continue;
                    const trait_index = plant_catalog.find(assignment.species_file) orelse return error.MissingPlantTraitProfile;
                    const traits = plant_catalog.entries.items[trait_index].traits;
                    const plant = try detailed_canopy_state.?.plantIndex(cell, species);
                    const branch_range = try detailed_canopy_state.?.branchRange(plant);
                    const adjusted = try ecosys.plant_initialization.adjustedPhenology(traits, runscript.phenology_initialization_parameters);
                    const concurrent = try ecosys.plant_initialization.concurrentNodeSettings(traits.functional_type.aboveground_turnover_type, adjusted.leaf_appearance_per_h, adjusted.floral_initiation_node_count_after_seed, runscript.phenology_initialization_parameters);
                    try branch_development_state.?.initializeRange(
                        branch_range.first,
                        branch_range.end,
                        adjusted.floral_initiation_node_count_after_seed,
                        adjusted.seed_initial_node_count,
                        traits.functional_type.growth_habit == 0 and traits.functional_type.leaf_phenology_type == 1,
                        concurrent.perennial_node_scaling,
                        concurrent.maximum_concurrently_growing_nodes,
                    );
                }
            }
            harvest_science_by_plant = try allocator.alloc(ecosys.plant_harvest_runtime.ScienceParameters, plant_count);
            @memset(harvest_science_by_plant.?, .{
                .carbon_woody_fraction = .{ 0, 1 },
                .leaf_nitrogen_woody_fraction = .{ 0, 1 },
                .sheath_nitrogen_woody_fraction = .{ 0, 1 },
                .leaf_phosphorus_woody_fraction = .{ 0, 1 },
                .sheath_phosphorus_woody_fraction = .{ 0, 1 },
            });
            for (harvest_science_by_plant.?, shoot_growth_plant_parameters) |*harvest_science, plant_parameters|
                harvest_science.nitrogen_fixation_type = plant_parameters.nitrogen_fixation_type;
            harvest_products_by_plant = try allocator.alloc(ecosys.plant_harvest_runtime.ProductLedger, plant_count);
            @memset(harvest_products_by_plant.?, .{});
            management_schedule_map = try ecosys.plant_management_dispatch.ScheduleMap.initMapped(allocator, state.cell_count, config.plant_populations, assignments, plant_unit_by_cell.?, plant_management_catalog);
            harvest_context = .{
                .canopy_state = &detailed_canopy_state.?,
                .canopy_structure_state = &canopy_structure_state.?,
                .canopy_layer_state = &canopy_layer_distribution_state.?,
                .branch_development = &branch_development_state.?,
                .science_by_plant = harvest_science_by_plant.?,
                .products_by_plant = harvest_products_by_plant.?,
                .leaf_area_presence_tolerance_m2 = config.absolute_tolerance,
                .plant_structural_presence_threshold_g_per_plant = runscript.plant_pool_parameters.branch_structural_presence_g_per_plant,
                .plant_tissue_presence_threshold_g_per_plant = runscript.plant_pool_parameters.plant_root_structural_presence_g_per_plant,
                .canopy_biochemistry_parameters_by_plant = canopy_biochemistry_parameters,
                .plant_phenology = &plant_phenology_state.?,
                .growth_stages = &plant_growth_stage_state.?,
                .emerged_by_plant = plant_phenology_state.?.emerged,
                .root_state = &plant_root_state.?,
                .root_litter_partition = &plant_litter_partition_state.?,
                .root_litter_carbon_ledger = &root_litter_carbon_ledger,
                .shoot_litter_carbon_g_c_by_plant = shoot_harvest_litter_carbon_g_c_by_plant,
                .shoot_litter_nitrogen_g_n_by_plant = shoot_harvest_litter_nitrogen_g_n_by_plant,
                .shoot_litter_phosphorus_g_p_by_plant = shoot_harvest_litter_phosphorus_g_p_by_plant,
                .soil_organic_state = &soil_organic_state,
                .surface_organic_state = &surface_organic_state,
                .surface_nutrient_state = &surface_fire_exchange_state,
                .daily_manure_carbon_input_g_c = daily_organic_fertilizer_carbon_input_g_c,
                .daily_manure_nitrogen_input_g_n = daily_organic_fertilizer_nitrogen_input_g_n,
                .daily_manure_phosphorus_input_g_p = daily_organic_fertilizer_phosphorus_input_g_p,
                .hourly_manure_products_by_plant = hourly_manure_products_by_plant,
                .grid = &state,
                .reseed_population_per_m2_by_plant = reseed_population_per_m2_by_plant,
                .cell_area_m2_by_cell = canopy_cell_area_m2,
            };
            for (plant_unit_by_cell.?, 0..) |unit_index, cell| {
                const active_species = assignments.units[unit_index].species;
                const first = cell * config.plant_populations;
                for (active_species, 0..) |assignment, species| checkpoint_species_names[first + species] = assignment.species_file;
                checkpoint_metadata_cells[cell] = .{
                    .species_names = checkpoint_species_names[first .. first + active_species.len],
                    .species_alive = plant_phenology_state.?.active[first .. first + active_species.len],
                };
            }
        }
    }
    var terrain_radiation_state = try ecosys.terrain_radiation.State.init(allocator, topography, unit_by_cell, canopy_geometry);
    defer terrain_radiation_state.deinit();
    const direct_incidence_count = try std.math.mul(usize, canopy_geometry.leaf_inclination_sine.len, canopy_geometry.leaf_azimuth_radians.len);
    const direct_incidence_buffer_count =
        try std.math.mul(usize, state.cell_count, direct_incidence_count);
    const direct_incidence_fraction =
        try allocator.alloc(f64, direct_incidence_buffer_count);
    defer allocator.free(direct_incidence_fraction);
    @memset(direct_incidence_fraction, 0);
    const direct_incidence_per_horizontal_area =
        try allocator.alloc(f64, direct_incidence_buffer_count);
    defer allocator.free(direct_incidence_per_horizontal_area);
    @memset(direct_incidence_per_horizontal_area, 0);
    const direct_scattering_direction =
        try allocator.alloc(
            ecosys.canopy_geometry.ScatteringDirection,
            direct_incidence_buffer_count,
        );
    defer allocator.free(direct_scattering_direction);
    @memset(direct_scattering_direction, .forward);
    var canopy_interception_state: ?ecosys.canopy_interception.State = null;
    defer if (canopy_interception_state) |*interception_state| interception_state.deinit();
    if (canopy_structure_state != null and canopy_optics_state != null) canopy_interception_state = try ecosys.canopy_interception.State.init(allocator, state.cell_count, config.plant_populations, runscript.canopy_layer_count);
    var ground_radiation_state = try ecosys.ground_radiation.State.initMapped(allocator, topography, unit_by_cell, soil_catalog, catalog_index_by_cell);
    defer ground_radiation_state.deinit();
    const initial_snow_temperature_k = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(initial_snow_temperature_k);
    @memset(initial_snow_temperature_k, site.mean_annual_air_temperature_c + 273.15);
    try snow_transport_state.initializePhysicalState(ground_radiation_state.initial_snow_depth_m, canopy_cell_area_m2, initial_snow_temperature_k, runscript.snow_layer_bottom_depth_m, runscript.initial_snow_density_Mg_per_m3);
    const snow_depth_m = try allocator.dupe(f64, ground_radiation_state.initial_snow_depth_m);
    defer allocator.free(snow_depth_m);
    @memset(surface_precipitation_state.solid_snow_water_equivalent_m3, 0);
    for (0..state.cell_count) |cell| {
        for (0..snow_transport_state.layer_capacity) |layer| surface_precipitation_state.solid_snow_water_equivalent_m3[cell] += snow_transport_state.solid_snow_water_equivalent_m3[cell * snow_transport_state.layer_capacity + layer];
    }
    try transport_hydrology_state.syncStorage(&state, &snow_transport_state);
    var canopy_exposure_state: ?ecosys.canopy_exposure.State = null;
    defer if (canopy_exposure_state) |*exposure_state| exposure_state.deinit();
    if (canopy_interception_state != null) canopy_exposure_state = try ecosys.canopy_exposure.State.init(allocator, state.cell_count, config.plant_populations);
    var canopy_airflow_state: ?ecosys.canopy_airflow.State = null;
    defer if (canopy_airflow_state) |*airflow_state| airflow_state.deinit();
    var canopy_air_exchange_state: ?ecosys.canopy_air_exchange.State = null;
    defer if (canopy_air_exchange_state) |*air_state| air_state.deinit();
    var standing_dead_air_exchange_state: ?ecosys.canopy_air_exchange.State = null;
    defer if (standing_dead_air_exchange_state) |*air_state| air_state.deinit();
    var canopy_surface_exchange_state: ?ecosys.canopy_surface_exchange.State = null;
    defer if (canopy_surface_exchange_state) |*exchange_state| exchange_state.deinit();
    var standing_dead_surface_exchange_state: ?ecosys.standing_dead_surface_exchange.State = null;
    defer if (standing_dead_surface_exchange_state) |*exchange_state| exchange_state.deinit();
    var canopy_surface_input_workspace: ?ecosys.canopy_surface_exchange.SurfaceInputWorkspace = null;
    defer if (canopy_surface_input_workspace) |*workspace| workspace.deinit();
    if (detailed_canopy_state != null and canopy_exposure_state != null and canopy_precipitation_retention_state != null) {
        canopy_airflow_state = try ecosys.canopy_airflow.State.init(allocator, state.cell_count, config.plant_populations);
        canopy_air_exchange_state = try ecosys.canopy_air_exchange.State.init(allocator, state.cell_count, config.plant_populations);
        standing_dead_air_exchange_state = try ecosys.canopy_air_exchange.State.init(allocator, state.cell_count, config.plant_populations);
        canopy_surface_exchange_state = try ecosys.canopy_surface_exchange.State.init(allocator, state.cell_count, config.plant_populations);
        standing_dead_surface_exchange_state = try ecosys.standing_dead_surface_exchange.State.init(allocator, state.cell_count, config.plant_populations);
        canopy_surface_input_workspace = try ecosys.canopy_surface_exchange.SurfaceInputWorkspace.init(allocator, runtime_plant_count);
    }
    const canopy_atmospheric_vapor_diffusivity_m2_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(canopy_atmospheric_vapor_diffusivity_m2_per_h);
    @memset(canopy_atmospheric_vapor_diffusivity_m2_per_h, 0);
    const soil_atmospheric_gas_conductance_m3_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(soil_atmospheric_gas_conductance_m3_per_h);
    @memset(soil_atmospheric_gas_conductance_m3_per_h, 0);
    const litter_atmospheric_gas_conductance_m3_per_h = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(litter_atmospheric_gas_conductance_m3_per_h);
    @memset(litter_atmospheric_gas_conductance_m3_per_h, 0);
    const snow_layer_gas_diffusivity_m2_per_h = try allocator.alloc(f64, snow_transport_state.active.len);
    defer allocator.free(snow_layer_gas_diffusivity_m2_per_h);
    @memset(snow_layer_gas_diffusivity_m2_per_h, 0);
    const canopy_surface_roughness_height_m = try allocator.alloc(f64, state.cell_count);
    defer allocator.free(canopy_surface_roughness_height_m);
    @memset(canopy_surface_roughness_height_m, runscript.surface_aerodynamic_parameters.soil_roughness_height_m);
    const canopy_available_intercepted_water_m3 = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(canopy_available_intercepted_water_m3);
    @memset(canopy_available_intercepted_water_m3, 0);
    const standing_dead_evaporation_m3_per_h = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(standing_dead_evaporation_m3_per_h);
    @memset(standing_dead_evaporation_m3_per_h, 0);
    var surface_energy_state = try ecosys.surface_energy.State.init(allocator, state.cell_count);
    defer surface_energy_state.deinit();
    var ecosystem_energy_ledger_state = try ecosys.ecosystem_energy_ledger.State.init(allocator, state.cell_count);
    defer ecosystem_energy_ledger_state.deinit();
    var plant_energy_publication_state = try ecosys.plant_energy_publication.State.init(
        allocator,
        runtime_plant_count,
    );
    defer plant_energy_publication_state.deinit();
    var plant_water_publication_state = try ecosys.plant_water_publication.State.init(
        allocator,
        state.cell_count,
        config.plant_populations,
    );
    defer plant_water_publication_state.deinit();
    var canopy_water_energy_publication_state =
        try ecosys.canopy_water_energy_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer canopy_water_energy_publication_state.deinit();
    var living_canopy_layer_publication_state =
        try ecosys.living_canopy_layer_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            runscript.canopy_layer_count,
        );
    defer living_canopy_layer_publication_state.deinit();
    var standing_dead_area_publication_state =
        try ecosys.standing_dead_area_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            runscript.canopy_layer_count,
        );
    defer standing_dead_area_publication_state.deinit();
    var root_uptake_ledger_state = try ecosys.root_uptake_ledger.State.init(allocator, state.layer_count);
    defer root_uptake_ledger_state.deinit();
    var root_water_uptake_publication_state =
        try ecosys.root_water_uptake_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_water_uptake_publication_state.deinit();
    var root_gas_content_publication_state =
        try ecosys.root_gas_content_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_gas_content_publication_state.deinit();
    var root_atmosphere_gas_publication_state =
        try ecosys.root_atmosphere_gas_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_atmosphere_gas_publication_state.deinit();
    var root_internal_gas_publication_state =
        try ecosys.root_internal_gas_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_internal_gas_publication_state.deinit();
    var root_soil_gas_publication_state =
        try ecosys.root_soil_gas_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_soil_gas_publication_state.deinit();
    var root_nutrient_uptake_publication_state =
        try ecosys.root_nutrient_uptake_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_nutrient_uptake_publication_state.deinit();
    var root_salt_uptake_publication_state =
        try ecosys.root_salt_uptake_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_salt_uptake_publication_state.deinit();
    var root_exudate_publication_state =
        try ecosys.root_exudate_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_exudate_publication_state.deinit();
    var root_competition_demand_publication_state =
        try ecosys.root_competition_demand_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_competition_demand_publication_state.deinit();
    var root_nitrogen_fixation_publication_state =
        try ecosys.root_nitrogen_fixation_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            state.soil_layer_capacity,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_nitrogen_fixation_publication_state.deinit();
    var canopy_carbon_publication_state =
        try ecosys.canopy_carbon_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer canopy_carbon_publication_state.deinit();
    const canopy_net_fixation_g_c_per_h_by_plant =
        try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(canopy_net_fixation_g_c_per_h_by_plant);
    @memset(canopy_net_fixation_g_c_per_h_by_plant, 0);
    var root_soil_element_exchange_publication_state =
        try ecosys.root_soil_element_exchange_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer root_soil_element_exchange_publication_state.deinit();
    const root_soil_element_exchange_workspace =
        try allocator.alloc(f64, 3 * runtime_plant_count);
    defer allocator.free(root_soil_element_exchange_workspace);
    @memset(root_soil_element_exchange_workspace, 0);
    const root_soil_carbon_exchange_g_c_per_h_by_plant =
        root_soil_element_exchange_workspace[0..runtime_plant_count];
    const root_soil_nitrogen_exchange_g_n_per_h_by_plant =
        root_soil_element_exchange_workspace[runtime_plant_count .. 2 * runtime_plant_count];
    const root_soil_phosphorus_exchange_g_p_per_h_by_plant =
        root_soil_element_exchange_workspace[2 * runtime_plant_count .. 3 * runtime_plant_count];
    var root_gas_withdrawal_publication_state =
        try ecosys.root_gas_withdrawal_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer root_gas_withdrawal_publication_state.deinit();
    var canopy_ammonia_publication_state =
        try ecosys.canopy_ammonia_publication.State.init(
            allocator,
            runtime_plant_count,
        );
    defer canopy_ammonia_publication_state.deinit();
    var cell_litter_standing_dead_publication_state =
        try ecosys.cell_litter_standing_dead_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer cell_litter_standing_dead_publication_state.deinit();
    var manure_deposition_publication_state =
        try ecosys.manure_deposition_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer manure_deposition_publication_state.deinit();
    var canopy_fire_publication_state =
        try ecosys.canopy_fire_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
        );
    defer canopy_fire_publication_state.deinit();
    var plant_combustion_publication_state =
        try ecosys.plant_combustion_publication.State.init(
            allocator,
            runtime_plant_count,
        );
    defer plant_combustion_publication_state.deinit();
    var root_combustion_salt_publication_state =
        try ecosys.root_combustion_salt_publication.State.init(
            allocator,
            state.cell_count,
            config.plant_populations,
            config.soil_layers,
            ecosys.plant_root_system.biological_domain_count,
        );
    defer root_combustion_salt_publication_state.deinit();
    const aboveground_litter_publication_workspace =
        try allocator.alloc(f64, 3 * runtime_plant_count);
    defer allocator.free(aboveground_litter_publication_workspace);
    @memset(aboveground_litter_publication_workspace, 0);
    const aboveground_litter_carbon_g_c_per_h_by_plant =
        aboveground_litter_publication_workspace[0..runtime_plant_count];
    const aboveground_litter_nitrogen_g_n_per_h_by_plant =
        aboveground_litter_publication_workspace[runtime_plant_count .. 2 * runtime_plant_count];
    const aboveground_litter_phosphorus_g_p_per_h_by_plant =
        aboveground_litter_publication_workspace[2 * runtime_plant_count .. 3 * runtime_plant_count];
    const canopy_net_radiation_mj = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(canopy_net_radiation_mj);
    const canopy_storage_heat_mj = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(canopy_storage_heat_mj);
    const zero_plant_energy_mj = try allocator.alloc(f64, runtime_plant_count);
    defer allocator.free(zero_plant_energy_mj);
    @memset(canopy_net_radiation_mj, 0);
    @memset(canopy_storage_heat_mj, 0);
    @memset(zero_plant_energy_mj, 0);
    const fire_active_this_hour = try allocator.alloc(bool, state.cell_count);
    defer allocator.free(fire_active_this_hour);
    @memset(fire_active_this_hour, false);
    var surface_temperature_solver_state = try ecosys.surface_temperature_solver.State.init(allocator, state.cell_count);
    defer surface_temperature_solver_state.deinit();
    var canopy_energy_state: ?ecosys.canopy_energy.State = null;
    defer if (canopy_energy_state) |*energy_state| energy_state.deinit();
    if (canopy_exposure_state != null) canopy_energy_state = try ecosys.canopy_energy.State.init(allocator, state.cell_count, config.plant_populations);
    var plant_water_balance_state: ?ecosys.plant_water_balance.State = null;
    defer if (plant_water_balance_state) |*balance| balance.deinit();
    var plant_water_workspace: ?ecosys.plant_water_balance.Workspace = null;
    defer if (plant_water_workspace) |*workspace| workspace.deinit();
    if (plant_root_state != null and detailed_canopy_state != null) {
        plant_water_balance_state = try ecosys.plant_water_balance.State.init(allocator, state.cell_count, config.plant_populations, config.soil_layers);
        plant_water_workspace = try ecosys.plant_water_balance.Workspace.init(allocator, state.cell_count, config.plant_populations, config.soil_layers);
        for (0..state.cell_count) |cell|
            plant_water_workspace.?.cell_area_m2[cell] =
                grid_environment.horizontal_cell_width_m[cell] *
                grid_environment.vertical_cell_width_m[cell];
        if (plant_assignments) |assignments| for (plant_unit_by_cell.?, 0..) |unit_index, cell| for (assignments.units[unit_index].species, 0..) |assignment, species| {
            const trait_index = plant_catalog.find(assignment.species_file) orelse return error.MissingPlantTraitProfile;
            const management_index = plant_management_catalog.find(assignment.management_file) orelse continue;
            const plant = cell * config.plant_populations + species;
            const traits = plant_catalog.entries.items[trait_index].traits;
            const planting = plant_management_catalog.entries.items[management_index].schedule.planting;
            plant_water_workspace.?.leaf_osmotic_potential_at_zero_total_mpa[plant] = traits.water_relations.osmotic_potential_mpa;
            plant_water_workspace.?.plant_population_count[plant] = planting.population_per_m2 * plant_water_workspace.?.cell_area_m2[cell];
            plant_water_workspace.?.root_porosity_fraction[plant] = traits.roots.root_porosity_fraction;
            plant_water_workspace.?.root_radial_resistivity_mpa_h_per_m3[plant] = traits.roots.radial_resistivity_mpa_h_per_m3;
            plant_water_workspace.?.root_axial_resistivity_mpa_h_per_m2[plant] = traits.roots.axial_resistivity_mpa_h_per_m2;
            plant_water_workspace.?.primary_root_radius_m[plant] = traits.roots.primary_root_radius_m;
            plant_water_workspace.?.secondary_root_radius_m[plant] = traits.roots.secondary_root_radius_m;
            plant_water_workspace.?.seeding_depth_m[plant] = planting.seed_depth_m;
            plant_water_workspace.?.woody_root_fraction[plant] = 1;
            plant_water_workspace.?.dynamically_woody[plant] = traits.functional_type.root_profile_type != 0 and traits.functional_type.growth_habit > 1;
            plant_water_workspace.?.vascular_growth_habit[plant] = traits.functional_type.growth_habit != 0;
        };
    }
    if (harvest_context) |*context| {
        context.root_woody_fraction_by_plant = if (plant_water_workspace) |*workspace| workspace.woody_root_fraction else null;
        context.carbon_exchange_state = if (canopy_carbon_exchange_state) |*exchange| exchange else null;
    }
    const surface_energy_settings: ecosys.surface_energy.Settings = .{
        .soil_longwave_emissivity = runscript.soil_longwave_emissivity,
        .snow_longwave_emissivity = runscript.snow_longwave_emissivity,
        .snow_full_cover_depth_m = runscript.snow_full_cover_depth_m,
    };
    const nonlinear_solver_options: ecosys.numerics.SolverOptions = .{
        .absolute_tolerance = config.absolute_tolerance,
        .relative_tolerance = config.relative_tolerance,
        .picard_relaxation = config.picard_relaxation,
        .max_iterations = config.max_nonlinear_iterations,
    };
    var current_atmospheric_co2_umol_per_mol = site.atmospheric_co2_umol_mol;
    const current_atmospheric_gas_concentration_g_per_m3 = try allocator.create([ecosys.gas_transport.species_count]f64);
    defer allocator.destroy(current_atmospheric_gas_concentration_g_per_m3);
    current_atmospheric_gas_concentration_g_per_m3.* = if (surface_gas_parameters) |parameters| parameters.atmospheric_concentration_g_per_m3 else [_]f64{0} ** ecosys.gas_transport.species_count;
    var active_tile_cells: ?[]const usize = null;
    const runtime_tile_plan: *const ecosys.spatial_grid.TilePlan =
        &tile_plan.?;
    var executed_weather_hours: usize = 0;
    const lateral_contribution_path = try std.fmt.allocPrint(
        allocator,
        "{s}.ecosys-ng-tile-io",
        // Workspace accepts one safe directory name relative to its parent;
        // input paths may be absolute or contain separators.
        .{std.fs.path.basename(args[1])},
    );
    defer allocator.free(lateral_contribution_path);
    var lateral_contribution_workspace =
        try ecosys.hourly_lateral_contribution_io.Workspace.init(
            init.io,
            std.Io.Dir.cwd(),
            lateral_contribution_path,
        );
    defer lateral_contribution_workspace.deinit();
    var hourly_science_context = .{
        .allocator = allocator,
        .executor = executor,
        .active_tile_cells = &active_tile_cells,
        .tile_plan = runtime_tile_plan,
        .executed_weather_hours = &executed_weather_hours,
        .restoring_checkpoint = scene_options[0].resume_from_checkpoint,
        .lateral_contribution_workspace = &lateral_contribution_workspace,
        .landscape_boundary_ledger = &landscape_mass_balance_state.boundary_ledger,
        .landscape_soil_mass_Mg_scratch = landscape_soil_mass_Mg_scratch,
        .grid = &state,
        .lateral_connection_mode_by_cell = lateral_connection_mode_by_cell,
        .site_by_cell = site_by_cell,
        .mean_annual_temperature_c_by_cell = mean_annual_temperature_c_by_cell,
        .mean_annual_temperature_k_by_cell = mean_annual_temperature_k_by_cell,
        .geothermal_enabled_by_cell = geothermal_enabled_by_cell,
        .surface_runoff_boundary_fraction_by_direction = surface_runoff_boundary_fraction_by_direction,
        .atmosphere = &atmospheric_state,
        .hourly_adjusted_shortwave_mj_per_m2 = hourly_adjusted_shortwave_mj_per_m2,
        .hourly_weather_reference_height_m = hourly_weather_reference_height_m,
        .hourly_extraterrestrial_shortwave_mj_per_m2 = hourly_extraterrestrial_shortwave_mj_per_m2,
        .hourly_solar_angle_sine = hourly_solar_angle_sine,
        .hourly_solar_azimuth_radians = hourly_solar_azimuth_radians,
        .irrigation_water_depth_m = irrigation_water_depth_m,
        .irrigation_dissolved_mass_g_per_m2 = irrigation_dissolved_mass_g_per_m2,
        .irrigation_hydrogen_mol_per_m2 = irrigation_hydrogen_mol_per_m2,
        .subsurface_irrigation_water_m3 = irrigation_loads.subsurface_water_m3,
        .irrigation_loads = &irrigation_loads,
        .canopy_radiation = &canopy_radiation_state,
        .canopy_optics = &canopy_optics_state,
        .canopy_geometry = &canopy_geometry,
        .canopy_cell_area_m2 = canopy_cell_area_m2,
        .surface_aerodynamics = &surface_aerodynamic_state,
        .ground_air = &ground_air_state,
        .surface_total_canopy_area_m2 = surface_total_canopy_area_m2,
        .surface_canopy_height_m = surface_canopy_height_m,
        .ground_air_vapor_pressure_kpa = ground_air_vapor_pressure_kpa,
        .atmospheric_vapor_fraction = atmospheric_vapor_fraction,
        .ground_air_canopy_resistance_h_per_m = ground_air_canopy_resistance_h_per_m,
        .ground_air_sensible_source_mj_per_h = ground_air_sensible_source_mj_per_h,
        .delayed_live_canopy_combustion_heat_mj = delayed_live_canopy_combustion_heat_mj,
        .delayed_standing_dead_combustion_heat_mj = delayed_standing_dead_combustion_heat_mj,
        .delayed_subsurface_combustion_heat_mj = delayed_subsurface_combustion_heat_mj,
        .delayed_surface_combustion_heat_mj = delayed_surface_combustion_heat_mj,
        .surface_combustion_heat_mj_per_m2 = surface_combustion_heat_mj_per_m2,
        .ground_air_vapor_source_m3_per_h = ground_air_vapor_source_m3_per_h,
        .ground_air_surface_sensible_conductance_mj_per_h_k = ground_air_surface_sensible_conductance_mj_per_h_k,
        .ground_air_surface_vapor_conductance_m3_per_h = ground_air_surface_vapor_conductance_m3_per_h,
        .ground_air_surface_vapor_fraction = ground_air_surface_vapor_fraction,
        .direct_incidence_fraction = direct_incidence_fraction,
        .direct_incidence_per_horizontal_area = direct_incidence_per_horizontal_area,
        .direct_scattering_direction = direct_scattering_direction,
        .canopy_interception = &canopy_interception_state,
        .canopy_structure = &canopy_structure_state,
        .terrain_radiation = &terrain_radiation_state,
        .ground_radiation = &ground_radiation_state,
        .snow_depth_m = snow_depth_m,
        .canopy_exposure = &canopy_exposure_state,
        .surface_energy = &surface_energy_state,
        .surface_energy_settings = surface_energy_settings,
        .surface_temperature = &surface_temperature_solver_state,
        .nonlinear_solver_options = nonlinear_solver_options,
        .soil_thermal = &soil_thermal_state,
        .soil_thermal_context = &soil_thermal_context,
        .soil_hourly_workspace = &soil_hourly_workspace,
        .soil_heat_solver_workspace = &soil_heat_solver_workspace,
        .soil_solver_properties = &soil_solver_property_state,
        .soil_geometry = &soil_geometry_state,
        .fertilizer_band = &fertilizer_band_state,
        .surface_pond_domain_workspace = &surface_pond_domain_workspace,
        .terrain_hydrology = &terrain_hydrology_state,
        .surface_runoff = &surface_runoff_state,
        .surface_inorganic_nitrogen_export_g_n_per_h = surface_inorganic_nitrogen_export_g_n_per_h,
        .surface_inorganic_phosphorus_export_g_p_per_h = surface_inorganic_phosphorus_export_g_p_per_h,
        .surface_organic_carbon_export_g_c_per_h = surface_organic_carbon_export_g_c_per_h,
        .surface_inorganic_carbon_export_g_c_per_h = surface_inorganic_carbon_export_g_c_per_h,
        .surface_organic_nitrogen_export_g_n_per_h = surface_organic_nitrogen_export_g_n_per_h,
        .surface_organic_phosphorus_export_g_p_per_h = surface_organic_phosphorus_export_g_p_per_h,
        .surface_erosion = &surface_erosion_state,
        .surface_soil_mass_at_erosion_start_Mg = surface_soil_mass_at_erosion_start_Mg,
        .net_sediment_Mg_per_h = net_sediment_Mg_per_h,
        .eroded_mineral_state = &eroded_mineral_state,
        .soil_transport_faces = &soil_transport_faces,
        .soil_face_geometry = &soil_face_geometry_state,
        .soil_solute_face_parameters = &soil_solute_face_parameter_state,
        .micropore_solute_face_flux_mol = micropore_solute_face_flux_mol,
        .macropore_solute_face_flux_mol = macropore_solute_face_flux_mol,
        .soil_organic_face_parameters = &soil_organic_face_parameter_state,
        .soil_organic_transport = &soil_organic_transport_state,
        .soil_organic_recharge_concentration_g_per_m3 = soil_organic_recharge_concentration_g_per_m3,
        .soil_dissolved_gas_face_parameters = &soil_dissolved_gas_face_parameter_state,
        .soil_dissolved_gas_transport = &soil_dissolved_gas_transport_state,
        .soil_dissolved_gas_recharge_concentration_g_per_m3 = soil_dissolved_gas_recharge_concentration_g_per_m3,
        .mineral_nitrogen_transport = &mineral_nitrogen_transport_state,
        .mineral_nitrogen_face_parameters = &mineral_nitrogen_face_parameters,
        .mineral_nitrogen_zone_fractions = mineral_nitrogen_zone_fractions,
        .soil_boundary_topology = &soil_boundary_topology_state,
        .transport_hydrology = &transport_hydrology_state,
        .snow_transport = &snow_transport_state,
        .snow_atmospheric_input_g = snow_atmospheric_input_g,
        .atmospheric_solute_input_ledger = &atmospheric_solute_input_ledger_state,
        .snow_surface_partitions = snow_surface_partitions,
        .snow_surface_discharge = snow_surface_discharge,
        .direct_surface_solute_input = direct_surface_solute_input,
        .micropore_solute_state = &micropore_solute_state,
        .macropore_solute_state = &macropore_solute_state,
        .soil_recharge_concentration_mol_per_m3 = soil_recharge_concentration_mol_per_m3,
        .soil_solute_boundary_net_flux_mol = soil_solute_boundary_net_flux_mol,
        .gas_transport = &gas_transport_state,
        .soil_gas_transport = &soil_gas_transport_state,
        .surface_litter_gas_transport = &surface_litter_gas_transport_state,
        .soil_reactive_nitrogen = &soil_reactive_nitrogen_state,
        .soil_microbial_phosphorus = &soil_microbial_phosphorus_state,
        .soil_microbial_turnover = &soil_microbial_turnover_state,
        .soil_litter_colonization = &soil_litter_colonization_state,
        .soil_organic_sorption = &soil_organic_sorption_state,
        .soil_organic_decomposition = &soil_organic_decomposition_state,
        .soil_organic_priming = &soil_organic_priming_state,
        .soil_respiration_products = &soil_respiration_products_state,
        .daily_heterotrophic_respiration = &daily_heterotrophic_respiration,
        .soil_autotrophic_carbon = &soil_autotrophic_carbon_state,
        .soil_microbial_layer_mixing = &soil_microbial_layer_mixing_state,
        .soil_biogeochemical_gas_fluxes = &soil_biogeochemical_gas_flux_state,
        .soil_methane = &soil_methane_state,
        .soil_microbial_oxygen = &soil_microbial_oxygen_state,
        .soil_oxygen_staging = &soil_oxygen_staging_state,
        .soil_field_capacity_fraction = soil_field_capacity_fraction,
        .soil_nitrogen_flux_workspace = &soil_nitrogen_flux_workspace,
        .soil_nitrifier_environment = &soil_nitrifier_environment_state,
        .soil_microbial = &soil_microbial_state,
        .soil_redox_satisfaction_fraction = soil_redox_satisfaction_fraction,
        .soil_nitrogen_parameters = &soil_nitrogen_parameters,
        .plant_available_nutrients = &plant_available_nutrient_state,
        .iteration_limits = iteration_limits,
        .config = config,
        .runscript = &runscript,
        .site = &site,
        .horizontal_cell_width_m = grid_environment.horizontal_cell_width_m,
        .vertical_cell_width_m = grid_environment.vertical_cell_width_m,
        .current_atmospheric_co2_umol_per_mol = &current_atmospheric_co2_umol_per_mol,
        .current_atmospheric_gas_concentration_g_per_m3 = current_atmospheric_gas_concentration_g_per_m3,
        .canopy_energy = &canopy_energy_state,
        .canopy_airflow = &canopy_airflow_state,
        .canopy_air_exchange = &canopy_air_exchange_state,
        .standing_dead_air_exchange = &standing_dead_air_exchange_state,
        .canopy_surface_exchange = &canopy_surface_exchange_state,
        .standing_dead_surface_exchange = &standing_dead_surface_exchange_state,
        .canopy_surface_input_workspace = &canopy_surface_input_workspace,
        .canopy_atmospheric_vapor_diffusivity_m2_per_h = canopy_atmospheric_vapor_diffusivity_m2_per_h,
        .soil_atmospheric_gas_conductance_m3_per_h = soil_atmospheric_gas_conductance_m3_per_h,
        .litter_atmospheric_gas_conductance_m3_per_h = litter_atmospheric_gas_conductance_m3_per_h,
        .snow_layer_gas_diffusivity_m2_per_h = snow_layer_gas_diffusivity_m2_per_h,
        .canopy_surface_roughness_height_m = canopy_surface_roughness_height_m,
        .canopy_available_intercepted_water_m3 = canopy_available_intercepted_water_m3,
        .standing_dead_evaporation_m3_per_h = standing_dead_evaporation_m3_per_h,
        .plants = &plant_state,
        .plant_water_balance = &plant_water_balance_state,
        .plant_water_workspace = &plant_water_workspace,
        .plant_roots = &plant_root_state,
        .plant_root_nutrient_workspace = &plant_root_nutrient_workspace,
        .plant_root_salt_workspace = &plant_root_salt_workspace,
        .plant_root_exudation_workspace = &plant_root_exudation_workspace,
        .plant_root_metabolism_workspace = &plant_root_metabolism_workspace,
        .plant_litter_partition = &plant_litter_partition_state,
        .plant_storage_remobilization_workspace = &plant_storage_remobilization_workspace,
        .plant_harvest = if (harvest_context) |*value| value else null,
        .root_nutrient_traits = root_nutrient_traits,
        .root_nutrient_feedback_enabled = root_nutrient_feedback_enabled,
        .root_metabolism_plant_parameters = root_metabolism_plant_parameters,
        .root_biological_domain_count_by_plant = root_biological_domain_count_by_plant,
        .plant_phenology = &plant_phenology_state,
        .plant_topology_controls = &plant_topology_controls,
        .plant_reproduction_controls = &plant_reproduction_controls,
        .canopy_layer_controls = &canopy_layer_controls,
        .branch_development = &branch_development_state,
        .plant_growth_stages = &plant_growth_stage_state,
        .plant_dormancy = &plant_dormancy_state,
        .phenology_root_oxygen_fraction = phenology_root_oxygen_fraction,
        .development_species_parameters = development_species_parameters,
        .development_dormancy_parameters = development_dormancy_parameters,
        .development_planting_day_of_year = development_planting_day_of_year,
        .development_planting_year = development_planting_year,
        .development_canopy_height_m = development_canopy_height_m,
        .development_surface_water_potential_mpa = development_surface_water_potential_mpa,
        .development_seed_layer_water_potential_mpa = development_seed_layer_water_potential_mpa,
        .development_emerged = plant_phenology_state.?.emerged,
        .root_gas_parameters = runscript.root_gas_parameters,
        .detailed_canopy = &detailed_canopy_state,
        .canopy_carbon_exchange = &canopy_carbon_exchange_state,
        .canopy_biochemistry_parameters = canopy_biochemistry_parameters,
        .shoot_growth_plant_parameters = shoot_growth_plant_parameters,
        .canopy_layer_distribution = &canopy_layer_distribution_state,
        .canopy_precipitation_retention = &canopy_precipitation_retention_state,
        .surface_precipitation = &surface_precipitation_state,
        .surface_pond_transition = &surface_pond_transition_state,
        .surface_pond_minimum_heat_capacity_mj_per_k = surface_pond_minimum_heat_capacity_mj_per_k,
        .surface_litter_geometry = &surface_litter_geometry_state,
        .surface_litter_water_environment = &surface_litter_water_environment_state,
        .surface_microbial_environment = &surface_microbial_environment_state,
        .surface_microbial_respiration = &surface_microbial_respiration_state,
        .surface_microbial_oxygen = &surface_microbial_oxygen_state,
        .surface_microbial_maintenance = &surface_microbial_maintenance_state,
        .surface_nonsymbiotic_nitrogen_fixation = &surface_nonsymbiotic_nitrogen_fixation_state,
        .surface_microbial_substrate_uptake = &surface_microbial_substrate_uptake_state,
        .surface_denitrification = &surface_denitrification_state,
        .surface_microbial_assimilation = &surface_microbial_assimilation_state,
        .surface_microbial_mineral_exchange = &surface_microbial_mineral_exchange_state,
        .surface_topsoil_mineral_exchange = &surface_topsoil_mineral_exchange_state,
        .surface_microbial_turnover = &surface_microbial_turnover_state,
        .surface_organic_priming = &surface_organic_priming_state,
        .surface_organic_decomposition = &surface_organic_decomposition_state,
        .surface_organic_sorption = &surface_organic_sorption_state,
        .surface_litter_colonization = &surface_litter_colonization_state,
        .surface_humification_fraction = surface_humification_fraction,
        .soil_chemistry = &initial_chemistry_state,
        .soil_chemistry_layer_parameters = soil_chemistry_layer_parameters,
        .soil_chemistry_solver_workspace = &soil_chemistry_solver_workspace,
        .eroded_chemistry_workspace = &eroded_chemistry_workspace,
        .soil_organic = &soil_organic_state,
        .soil_organic_carbon_at_hour_start_g_c = soil_organic_carbon_at_hour_start_g_c,
        .soil_organic_carbon_change_g_c_per_h = soil_organic_carbon_change_g_c_per_h,
        .topsoil_humus_partition = topsoil_humus_partition,
        .surface_litter_fertilizer = &surface_litter_fertilizer_state,
        .soil_fertilizer_inventory = &soil_fertilizer_inventory,
        .mineral_fertilizer_inventory = &mineral_fertilizer_inventory,
        .eroded_fertilizer_workspace = &eroded_fertilizer_workspace,
        .surface_litter_fertilizer_diagnostics = &surface_litter_fertilizer_diagnostics,
        .surface_field_capacity_potential_mpa = surface_field_capacity_potential_mpa,
        .surface_wilting_point_potential_mpa = surface_wilting_point_potential_mpa,
        .surface_litter_ice_m3 = surface_litter_ice_m3,
        .surface_charcoal_carbon_g_c = surface_charcoal_carbon_g_c,
        .litter_gas_transport = &litter_gas_transport_state,
        .surface_gas_parameters = &surface_gas_parameters,
        .surface_litter_chemistry = &surface_litter_chemistry_state,
        .surface_litter_chemistry_diagnostics = &surface_litter_chemistry_diagnostics,
        .surface_litter_cation_selectivity = surface_litter_cation_selectivity,
        .surface_organic = &surface_organic_state,
        .shoot_senescence_products_by_plant = shoot_senescence_products_by_plant,
        .root_litter_products_by_plant = root_litter_products_by_plant,
        .root_litter_carbon_ledger = &root_litter_carbon_ledger,
        .eroded_organic_workspace = &eroded_organic_workspace,
        .organic_parameters = &organic_parameters,
        .chemistry_reaction_parameters = &chemistry_reaction_parameters,
    };
    var climate_state: ecosys.climate_change.State = .{};
    var executed_scene_passes: usize = 0;
    var resume_position: ?ecosys.checkpoint_manifest.SimulationInstant = null;
    var first_pass_preview = try ecosys.simulation_timeline.PassIterator.init(runscript.scenarios, runscript.execution_repeat_count);
    const first_requested_pass = first_pass_preview.next() orelse return error.EmptySimulationTimeline;
    hourly_science_context.restoring_checkpoint =
        scene_options[first_requested_pass.scene_index].resume_from_checkpoint;
    // Legacy multi-scene runs mark later scenes as continuation scenes. They
    // already inherit live state in this process; only a resume flag on the
    // first executed scene requests an external bundle restore.
    if (scene_options[first_requested_pass.scene_index].resume_from_checkpoint) {
        if (plant_root_state == null or detailed_canopy_state == null or canopy_precipitation_retention_state == null or plant_phenology_state == null or plant_growth_stage_state == null or plant_dormancy_state == null or branch_development_state == null) return error.CheckpointResumeRequiresCompletePlantState;
        var maximum_species_name_bytes: usize = 1;
        for (checkpoint_species_names) |name| maximum_species_name_bytes = @max(maximum_species_name_bytes, name.len);
        const canopy = &detailed_canopy_state.?;
        const maximum_branches = canopy.plant_branch_offsets[canopy.plant_branch_offsets.len - 1];
        const maximum_nodes = canopy.branch_node_offsets[canopy.branch_node_offsets.len - 1];
        const maximum_samples = canopy.node_sample_offsets[canopy.node_sample_offsets.len - 1];
        var restored = try ecosys.checkpoint_bundle_reader.read(
            allocator,
            init.io,
            std.Io.Dir.cwd(),
            config,
            .{
                .manifest = .{ .maximum_columns = config.grid_columns, .maximum_rows = config.grid_rows, .maximum_soil_layers = config.soil_layers, .maximum_snow_layers = snow_transport_state.layer_capacity, .maximum_plant_species_per_cell = config.plant_populations, .maximum_root_axes_per_plant = runscript.root_axes_per_plant },
                .plant_metadata = .{ .maximum_cells = state.cell_count, .maximum_species_per_cell = config.plant_populations, .maximum_species_name_bytes = maximum_species_name_bytes },
                .plant_development = .{ .maximum_cells = state.cell_count, .maximum_species = config.plant_populations, .maximum_branches = plant_growth_stage_state.?.branches.len },
                .plant_roots = .{ .maximum_plants = plant_root_state.?.plant_count, .maximum_soil_layers = config.soil_layers, .maximum_root_axes = runscript.root_axes_per_plant },
                .plant_canopy = .{ .maximum_cells = state.cell_count, .maximum_species = config.plant_populations, .maximum_branches = maximum_branches, .maximum_nodes = maximum_nodes, .maximum_samples = maximum_samples, .maximum_layers = runscript.canopy_layer_count, .maximum_inclinations = canopy_geometry.leaf_inclination_sine.len, .maximum_azimuths = canopy_geometry.leaf_azimuth_radians.len },
                .soil_biogeochemistry = .{ .maximum_cells = state.cell_count, .maximum_layers = config.soil_layers, .maximum_substrates = runscript.microbial_substrate_count, .maximum_populations = runscript.microbial_population_count },
                .soil_organic = .{ .maximum_profile_layers = state.layer_count, .maximum_surface_cells = state.cell_count },
                .transport = .{ .maximum_transport_cells = state.layer_count, .maximum_solute_species = micropore_solute_state.species_count, .maximum_snow_cells = state.cell_count, .maximum_snow_layers = snow_transport_state.layer_capacity },
                .soil_geometry = .{ .maximum_columns = config.grid_columns, .maximum_rows = config.grid_rows, .maximum_soil_layers = config.soil_layers, .maximum_snow_layers = snow_transport_state.layer_capacity, .maximum_plants = runtime_plant_count },
            },
            .{ .manifest_buffer_bytes = 16 * 1024, .section_read_buffer_bytes = 64 * 1024, .section_verify_buffer_bytes = 64 * 1024 },
        );
        const restored_generation = restored.manifest.generation;
        resume_position = restored.manifest.instant;
        try ecosys.checkpoint_bundle_reader.swapIntoLive(&restored, .{
            .grid = &state,
            .plants = &plant_state,
            .plant_development = .{ .phenology = &plant_phenology_state.?, .growth = &plant_growth_stage_state.?, .dormancy = &plant_dormancy_state.?, .branch_development = &branch_development_state.? },
            .plant_roots = &plant_root_state.?,
            .plant_canopy = .{ .canopy = &detailed_canopy_state.?, .retention = &canopy_precipitation_retention_state.?, .layer_distribution = &canopy_layer_distribution_state.? },
            .soil_biogeochemistry = .{ .microbial = &soil_microbial_state, .chemistry = &initial_chemistry_state, .available_nutrients = &plant_available_nutrient_state, .fertilizer = &soil_fertilizer_inventory, .mineral_fertilizer = &mineral_fertilizer_inventory, .fertilizer_band = &fertilizer_band_state, .reactive_nitrogen = &soil_reactive_nitrogen_state, .microbial_phosphorus = &soil_microbial_phosphorus_state },
            .soil_organic_matter = .{ .profile = &soil_organic_state, .surface = &surface_organic_state, .litter_chemistry = &surface_litter_chemistry_state, .litter_fertilizer = &surface_litter_fertilizer_state, .surface_respiration = &surface_microbial_respiration_state, .surface_denitrification = &surface_denitrification_state, .surface_fire_exchange = &surface_fire_exchange_state, .litter_salt_ingress = &plant_litter_salt_ingress_state },
            .transport = .{ .micropore = &micropore_solute_state, .macropore = &macropore_solute_state, .mineral_nitrogen = &mineral_nitrogen_transport_state, .organic = &soil_organic_transport_state, .gas = &gas_transport_state, .litter_gas = &litter_gas_transport_state, .snow = &snow_transport_state, .surface = &surface_transport_state },
            .soil_geometry_and_hydrology = .{ .geometry = &soil_geometry_state, .hydrology = &transport_hydrology_state, .surface = &surface_precipitation_state, .erosion = &surface_erosion_state, .climate = &climate_state, .eroded_minerals = &eroded_mineral_state, .runtime = .{ .soil_properties = &soil_solver_property_state, .soil_thermal = &soil_thermal_state }, .surface_boundary = .{ .ground_air = &ground_air_state, .surface_aerodynamics = &surface_aerodynamic_state }, .surface_litter_geometry = &surface_litter_geometry_state, .surface_litter_ice_m3 = surface_litter_ice_m3, .delayed_live_canopy_combustion_heat_mj = delayed_live_canopy_combustion_heat_mj, .delayed_standing_dead_combustion_heat_mj = delayed_standing_dead_combustion_heat_mj, .delayed_subsurface_combustion_heat_mj = delayed_subsurface_combustion_heat_mj, .delayed_surface_combustion_heat_mj = delayed_surface_combustion_heat_mj },
            .landscape_mass_balance = &landscape_mass_balance_state,
        });
        // The owner swap replaces these slices. Rebind consumers before the
        // restored bundle releases the former live allocation.
        hourly_science_context.development_emerged = plant_phenology_state.?.emerged;
        hourly_science_context.restoring_checkpoint = true;
        harvest_context.?.emerged_by_plant = plant_phenology_state.?.emerged;
        if (plant_assignments) |assignments| for (plant_unit_by_cell.?, 0..) |unit_index, cell| {
            const active_species = assignments.units[unit_index].species.len;
            const first = cell * config.plant_populations;
            checkpoint_metadata_cells[cell].species_alive = plant_phenology_state.?.active[first .. first + active_species];
        };
        restored.deinit();
        // The checkpointed thermal owner is authoritative. Its accepted
        // layer thickness and volume can intentionally differ from the soil
        // geometry carrier after freeze/thaw and remapping; rebuilding either
        // value here changes the next surface heat-storage equation.
        // The same rule applies to the checkpointed runtime soil-property
        // owner. Erosion texture, organic-carbon concentration, exchange
        // capacities, and rainfall conductivity damage are accepted mutable
        // state, not restart-time derivations.
        try ecosys.transport_hydrology.refreshSoilFacesFromHydrology(
            &soil_transport_faces,
            &transport_hydrology_state,
        );
        try soil_face_geometry_state.refreshMapped(
            &state,
            &soil_transport_faces,
            soil_solver_property_state.layer_thickness_m,
            grid_environment.horizontal_cell_width_m,
            grid_environment.vertical_cell_width_m,
        );
        executed_weather_hours = std.math.cast(usize, restored_generation) orelse return error.CheckpointGenerationExceedsRuntimeCounter;
    }
    var pass_iterator = try ecosys.simulation_timeline.PassIterator.init(runscript.scenarios, runscript.execution_repeat_count);
    while (pass_iterator.next()) |pass| {
        {
            const scene = runscript.scenes[pass.scene_index];
            const active_options = scene_options[pass.scene_index];
            try visualization_streams.configureScene(
                active_options.visualization_enabled,
                active_options.visualization_start_year,
                active_options.visualization_end_year,
            );
            if (resume_position) |position| {
                const before_target = pass.execution_iteration < position.execution_iteration or
                    (pass.execution_iteration == position.execution_iteration and pass.scenario_index < position.scenario_index) or
                    (pass.execution_iteration == position.execution_iteration and pass.scenario_index == position.scenario_index and pass.scenario_iteration < position.scenario_iteration) or
                    (pass.execution_iteration == position.execution_iteration and pass.scenario_index == position.scenario_index and pass.scenario_iteration == position.scenario_iteration and pass.scene_index < position.scene_index);
                if (before_target) continue;
                if (pass.execution_iteration != position.execution_iteration or pass.scenario_index != position.scenario_index or pass.scenario_iteration != position.scenario_iteration or pass.scene_index != position.scene_index) return error.CheckpointTimelinePositionNotInRunscript;
            }
            const active_iteration_limits = try ecosys.iteration_control.Limits.fromSceneOptions(active_options);
            hourly_science_context.iteration_limits = active_iteration_limits;
            hourly_science_context.nonlinear_solver_options.max_iterations = active_iteration_limits.water_heat_solute_max_iterations;

            const weather_assignment = scene_weather_assignments[pass.scene_index];
            const weather_streams = try allocator.alloc(
                *ecosys.weather.WeatherStream,
                weather_assignment.unique_file_names.len,
            );
            defer allocator.free(weather_streams);
            var opened_weather_stream_count: usize = 0;
            defer for (weather_streams[0..opened_weather_stream_count]) |stream| stream.deinit();
            const hourly_weather_streams = try allocator.alloc(
                ecosys.hourly_weather_stream.Stream,
                weather_assignment.unique_file_names.len,
            );
            defer allocator.free(hourly_weather_streams);
            for (
                weather_assignment.unique_file_names,
                weather_assignment.forcing_geometry_by_stream,
                0..,
            ) |weather_name, geometry, stream_index| {
                const weather_path = try resolveInputPath(
                    allocator,
                    init.io,
                    runscript_directory,
                    weather_name,
                );
                defer allocator.free(weather_path);
                weather_streams[stream_index] = ecosys.weather.WeatherStream.init(
                    allocator,
                    init.io,
                    weather_path,
                    64 * 1024,
                ) catch |err| {
                    std.log.err("weather stream open failed: execution={d} scenario={d} scenario_repeat={d} scene={d} file='{s}' error={s}", .{ pass.execution_iteration + 1, pass.scenario_index + 1, pass.scenario_iteration + 1, pass.scene_index + 1, weather_name, @errorName(err) });
                    return err;
                };
                hourly_science_context.restoring_checkpoint = false;
                opened_weather_stream_count += 1;
                hourly_weather_streams[stream_index] = try ecosys.hourly_weather_stream.Stream.init(
                    weather_streams[stream_index],
                    geometry.altitude_m,
                    geometry.latitude_degrees_north,
                    geometry.phytotron,
                    .{
                        .snowfall_temperature_threshold_c = active_options.snowfall_temperature_threshold_c,
                        .minimum_snowfall_water_equivalent_m = active_options.minimum_snowfall_water_equivalent_m,
                    },
                );
            }
            const weather_stream = weather_streams[0];
            const fertilizer_solar_noon_hour_by_cell = try allocator.alloc(u8, state.cell_count);
            defer allocator.free(fertilizer_solar_noon_hour_by_cell);
            for (fertilizer_solar_noon_hour_by_cell, 0..) |*solar_noon, cell| {
                const stream_index = weather_assignment.stream_index_by_cell[cell];
                if (stream_index >= weather_streams.len)
                    return error.WeatherGridStreamIndexOutOfRange;
                const value = weather_streams[stream_index].header.solar_noon_hour;
                if (!std.math.isFinite(value) or value < 0 or value >= 25)
                    return error.InvalidHourlyFertilizerSchedule;
                solar_noon.* = @intFromFloat(@trunc(value));
            }
            const weather_hour_by_stream = try allocator.alloc(
                ecosys.hourly_weather_stream.Observation,
                hourly_weather_streams.len,
            );
            defer allocator.free(weather_hour_by_stream);
            var previous_weather_timestamp: ?ecosys.weather.Timestamp = null;
            var resolved_development_year: ?u16 = null;
            var scene_weather_hours: usize = 0;
            if (resume_position) |position| {
                const completed_scene_hours = std.math.cast(usize, position.completed_scene_hours) orelse return error.CheckpointSceneHourExceedsRuntimeCounter;
                for (0..completed_scene_hours) |_| {
                    for (hourly_weather_streams, 0..) |*stream, stream_index| {
                        const skipped = try stream.next() orelse return error.CheckpointSceneHourExceedsWeatherStream;
                        if (stream_index == 0) {
                            previous_weather_timestamp = skipped.timestamp;
                        } else if (!sameWeatherTimestamp(previous_weather_timestamp.?, skipped.timestamp)) {
                            return error.WeatherGridTimestampMismatch;
                        }
                    }
                }
                const checkpoint_timestamp = previous_weather_timestamp orelse return error.CheckpointHasNoCompletedWeatherHour;
                const checkpoint_year = checkpoint_timestamp.year orelse active_options.start_date.year;
                const checkpoint_day = checkpoint_timestamp.day_of_year orelse return error.CheckpointWeatherRequiresDayOfYear;
                const checkpoint_hour = if (checkpoint_timestamp.hour == 24) 23 else checkpoint_timestamp.hour;
                if (checkpoint_year != position.year or checkpoint_day != position.day_of_year or checkpoint_hour != position.hour) return error.CheckpointWeatherTimestampMismatch;
                scene_weather_hours = completed_scene_hours;
                resume_position = null;
            }
            // EXEC IBEGIN/ISTART reset: establish a new seven-domain
            // reference only at a fresh scene boundary. A resumed scene keeps
            // the checkpointed monitor and cumulative boundary history.
            if (scene_weather_hours == 0) {
                const totals = try reconstructLandscapeMassBalance(
                    hourly_science_context,
                );
                if (landscape_mass_balance_state.monitor) |*monitor|
                    try monitor.reset(totals)
                else
                    landscape_mass_balance_state.monitor =
                        try ecosys.mass_balance_audit.Monitor.init(
                            totals,
                            config.absolute_tolerance,
                        );
            } else if (landscape_mass_balance_state.monitor == null) {
                return error.CheckpointResumeMissingMassBalanceMonitor;
            }
            const scene_expected_hours = try std.math.mul(usize, try ecosys.simulation_timeline.inclusiveDays(active_options.start_date, active_options.end_date), 24);
            while (scene_weather_hours < scene_expected_hours) {
                @memset(fire_active_this_hour, false);
                for (hourly_weather_streams, 0..) |*stream, stream_index| {
                    weather_hour_by_stream[stream_index] = try stream.next() orelse {
                        std.log.err("weather stream ended before scene: scene={d} expected_hours={d} received_hours={d} file='{s}'", .{
                            pass.scene_index + 1,
                            scene_expected_hours,
                            scene_weather_hours,
                            weather_assignment.unique_file_names[stream_index],
                        });
                        return error.WeatherStreamEndedBeforeScene;
                    };
                    if (stream_index > 0 and !sameWeatherTimestamp(
                        weather_hour_by_stream[0].timestamp,
                        weather_hour_by_stream[stream_index].timestamp,
                    )) {
                        std.log.err("weather-grid timestamp mismatch: scene={d} hour_index={d} first_file='{s}' mismatched_file='{s}'", .{
                            pass.scene_index + 1,
                            scene_weather_hours,
                            weather_assignment.unique_file_names[0],
                            weather_assignment.unique_file_names[stream_index],
                        });
                        return error.WeatherGridTimestampMismatch;
                    }
                }
                const weather_hour = weather_hour_by_stream[0];
                var timestamp = weather_hour.timestamp;
                if (timestamp.year == null) timestamp.year = active_options.start_date.year;
                // HOUR1 clears hourly plant ledgers before GROSUB management.
                // Management dispatch below may publish disturbance, harvest,
                // or combustion fluxes that must survive into this hour's
                // HCNET and output aggregation.
                if (plant_root_state) |*roots| roots.resetHourlyFluxes();
                if (canopy_carbon_exchange_state) |*ledger| ledger.reset();
                const current_year = timestamp.year.?;
                const current_day_of_year = try dayOfYearFromTimestamp(timestamp);
                if (resolved_development_year == null or resolved_development_year.? != current_year) {
                    for (development_planting_dates, 0..) |date, plant| {
                        const resolved = try date.resolve(current_year);
                        development_planting_day_of_year[plant] = try date.dayOfYear(current_year);
                        development_planting_year[plant] = resolved.year;
                    }
                    resolved_development_year = current_year;
                }
                const begins_day = if (previous_weather_timestamp) |previous| !sameCalendarDay(previous, timestamp) else true;
                if (begins_day) {
                    daily_soil_gas_flux.reset();
                    daily_canopy_gas_exchange.reset();
                    plant_water_publication_state.resetDaily();
                }
                if (begins_day) daily_heterotrophic_respiration.resetDaily();
                if (begins_day) {
                    daily_carbon_export.resetDaily();
                    daily_phosphorus_export.resetDaily();
                    daily_nitrogen_export.resetDaily();
                    daily_water_ledger.resetDaily();
                    daily_heat_ledger.reset();
                    atmospheric_solute_input_ledger_state.resetDaily();
                    @memset(daily_soil_fire_charcoal_production_g_c, 0);
                    @memset(daily_soil_fire_carbon_dioxide_emission_g_c, 0);
                    @memset(daily_soil_fire_methane_emission_g_c, 0);
                    @memset(daily_soil_fire_phosphorus_flux_g_p, 0);
                    @memset(daily_soil_fire_nitrogen_flux_g_n, 0);
                    @memset(daily_soil_combusted_nitrogen_g_n, 0);
                    @memset(daily_soil_combusted_phosphorus_g_p, 0);
                    @memset(daily_microbial_phosphate_mineralization_g_p, 0);
                    @memset(daily_microbial_nitrogen_mineralization_g_n, 0);
                    @memset(daily_organic_fertilizer_carbon_input_g_c, 0);
                    @memset(daily_organic_fertilizer_phosphorus_input_g_p, 0);
                    @memset(daily_organic_fertilizer_nitrogen_input_g_n, 0);
                    soil_fertilizer_inventory.resetDaily();
                    mineral_fertilizer_inventory.resetDaily();
                    @memset(daily_biome_organic_carbon_input_g_c, 0);
                }
                // HFUNC/GROSUB evaluates IFLGC/IFLGI at NFZ=1, the first
                // source substep of every hour. The converged hourly model
                // therefore evaluates activation once per hour, while the
                // active/lifecycle flags keep reconstruction single-shot.
                if (management_schedule_map != null) {
                    const management_date = try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp);
                    if (plant_phenology_state != null and plant_growth_stage_state != null) {
                        const activation_count = try ecosys.plant_management_dispatch.refreshPlantActivity(management_schedule_map.?, plant_management_catalog, management_date, &plant_phenology_state.?, &plant_growth_stage_state.?, plant_phenology_state.?.emerged, plant_phenology_state.?.lifecycle_initialized, newly_activated_plants);
                        if (activation_count > 0) {
                            const assignments = plant_assignments orelse return error.PlantActivationRequiresAssignments;
                            var topology_changed = false;
                            for (newly_activated_plants, 0..) |activated, plant| {
                                if (!activated) continue;
                                const cell = plant / config.plant_populations;
                                const species = plant % config.plant_populations;
                                const unit_index = plant_unit_by_cell.?[cell];
                                if (species >= assignments.units[unit_index].species.len) return error.PlantActivationSpeciesOutOfBounds;
                                const assignment = assignments.units[unit_index].species[species];
                                const traits = plant_catalog.entries.items[plant_catalog.find(assignment.species_file) orelse return error.MissingPlantTraitProfile].traits;
                                const planting = plant_management_catalog.entries.items[plant_management_catalog.find(assignment.management_file) orelse return error.MissingPlantManagementSchedule].schedule.planting;
                                var cumulative_depth_m: f64 = 0;
                                for (replant_layer_bottom_depth_m, 0..) |*bottom, layer| {
                                    cumulative_depth_m += soil_solver_property_state.layer_thickness_m[cell * config.soil_layers + layer];
                                    bottom.* = cumulative_depth_m;
                                }
                                const planting_layer = try ecosys.plant_initialization.plantingLayer(planting.seed_depth_m, replant_layer_bottom_depth_m);
                                const adjusted = try ecosys.plant_initialization.adjustedPhenology(traits, runscript.phenology_initialization_parameters);
                                const concurrent = try ecosys.plant_initialization.concurrentNodeSettings(traits.functional_type.aboveground_turnover_type, adjusted.leaf_appearance_per_h, adjusted.floral_initiation_node_count_after_seed, runscript.phenology_initialization_parameters);
                                const storage = try ecosys.plant_initialization.seedStorage(traits.morphology.seed_mass_at_planting_g, planting.population_per_m2, canopy_cell_area_m2[cell], traits.organ_nitrogen_to_carbon_ratio.grain, traits.organ_phosphorus_to_carbon_ratio.grain);
                                const thermal = try ecosys.plant_initialization.thermalAcclimation(assignment.species_file, traits.functional_type.photosynthesis_pathway, traits.functional_type.thermal_adaptation_zone, runscript.thermal_acclimation_parameters);
                                const preserve_reseed_inventories = plant_phenology_state.?.reseed_pending[plant];
                                const preserve_death_standing_dead = plant_phenology_state.?.death_replant_pending[plant];
                                const preserve_persistent_inventories = preserve_reseed_inventories or preserve_death_standing_dead;
                                const persistent_reseed_inventories = if (preserve_persistent_inventories)
                                    try ecosys.canopy_photosynthesis.capturePersistentReseedInventories(&detailed_canopy_state.?, plant)
                                else
                                    undefined;
                                const carried_standing_dead_surface_water_m3 = canopy_precipitation_retention_state.?.standing_dead_surface_water_m3[plant];

                                const old_branch_range = try detailed_canopy_state.?.branchRange(plant);
                                const old_node_range = try detailed_canopy_state.?.nodeRange(old_branch_range.first);
                                topology_changed = topology_changed or old_branch_range.end - old_branch_range.first != 1 or old_node_range.end - old_node_range.first != 1;
                                try ecosys.plant_topology.compactPlantToInitialTopology(.{
                                    .canopy = &detailed_canopy_state.?,
                                    .growth_stages = &plant_growth_stage_state.?,
                                    .dormancy = &plant_dormancy_state.?,
                                    .branch_development = &branch_development_state.?,
                                }, plant);
                                try detailed_canopy_state.?.clearPlantForReconstruction(plant);
                                try plant_root_state.?.reconstructPlant(plant, traits, planting_layer, planting.seed_depth_m, runscript.root_initialization_parameters);
                                const seed_geometry = try ecosys.plant_initialization.seedGeometry(traits.morphology.seed_mass_at_planting_g, runscript.plant_geometry_parameters);
                                try plant_root_state.?.setSeedGeometry(plant, seed_geometry.volume_m3, seed_geometry.length_m, seed_geometry.surface_area_m2);
                                try plant_litter_partition_state.?.initializePlant(plant, traits, runscript.standing_dead_partition_parameters);
                                const branch_range = try plant_growth_stage_state.?.clearPlantForReconstruction(plant);
                                try plant_dormancy_state.?.clearRangeForReconstruction(branch_range.first, branch_range.end);
                                try branch_development_state.?.clearRangeForReconstruction(branch_range.first, branch_range.end);
                                const initial_branch_end = @min(branch_range.end, branch_range.first + 1);
                                try branch_development_state.?.initializeRange(branch_range.first, initial_branch_end, adjusted.floral_initiation_node_count_after_seed, adjusted.seed_initial_node_count, traits.functional_type.growth_habit == 0 and traits.functional_type.leaf_phenology_type == 1, concurrent.perennial_node_scaling, concurrent.maximum_concurrently_growing_nodes);
                                for (plant_growth_stage_state.?.branches[initial_branch_end..branch_range.end]) |*branch| branch.dead = true;
                                plant_phenology_state.?.floral_initiated[plant] = false;
                                plant_phenology_state.?.node_initiation_rate_at_25c_per_h[plant] = adjusted.node_initiation_per_h;
                                plant_phenology_state.?.leaf_appearance_rate_at_25c_per_h[plant] = adjusted.leaf_appearance_per_h;
                                plant_phenology_state.?.initiated_node_count[plant] = adjusted.seed_initial_node_count;
                                plant_phenology_state.?.appeared_leaf_count[plant] = adjusted.seed_initial_node_count;
                                plant_phenology_state.?.node_initiation_per_timestep[plant] = 0;
                                plant_phenology_state.?.leaf_appearance_per_timestep[plant] = 0;
                                if (branch_range.first < branch_range.end) {
                                    plant_growth_stage_state.?.branches[branch_range.first].initiated_node_count = adjusted.seed_initial_node_count;
                                    plant_growth_stage_state.?.branches[branch_range.first].appeared_leaf_count = adjusted.seed_initial_node_count;
                                }
                                try ecosys.plant_initialization.initializeShootControls(&detailed_canopy_state.?, plant, planting.population_per_m2, canopy_cell_area_m2[cell], traits.water_relations.cuticular_resistance_s_per_m, traits.functional_type.photosynthesis_pathway, runscript.shoot_control_parameters);
                                try ecosys.plant_initialization.initializeThermalAcclimation(&detailed_canopy_state.?, plant, thermal);
                                if (!preserve_reseed_inventories) {
                                    try ecosys.plant_initialization.initializeSeedStorage(&detailed_canopy_state.?, plant, storage);
                                }
                                const initialized_seed_storage = ecosys.canopy_photosynthesis.ElementalMass{
                                    .carbon_g = detailed_canopy_state.?.plant_seed_storage_carbon_g[plant],
                                    .nitrogen_g = detailed_canopy_state.?.plant_seed_storage_nitrogen_g[plant],
                                    .phosphorus_g = detailed_canopy_state.?.plant_seed_storage_phosphorus_g[plant],
                                };
                                const standing_dead = try ecosys.plant_initialization.standingDeadStorage(traits.morphology.standing_dead_carbon_g_per_m2, canopy_cell_area_m2[cell], traits.organ_nitrogen_to_carbon_ratio.stalk, traits.organ_phosphorus_to_carbon_ratio.stalk, runscript.standing_dead_partition_parameters);
                                if (!preserve_reseed_inventories) {
                                    try ecosys.plant_initialization.initializeStandingDeadStorage(&detailed_canopy_state.?, plant, standing_dead);
                                }
                                try ecosys.plant_initialization.initializePlantHeatAndWater(&plant_state, &detailed_canopy_state.?, plant, site.mean_annual_air_temperature_c, traits.water_relations.osmotic_potential_mpa, runscript.plant_heat_water_parameters);
                                if (preserve_persistent_inventories) {
                                    try ecosys.canopy_photosynthesis.restorePersistentReseedInventories(&detailed_canopy_state.?, plant, persistent_reseed_inventories);
                                    if (preserve_death_standing_dead and !preserve_reseed_inventories) {
                                        detailed_canopy_state.?.plant_seed_storage_carbon_g[plant] = initialized_seed_storage.carbon_g;
                                        detailed_canopy_state.?.plant_seed_storage_nitrogen_g[plant] = initialized_seed_storage.nitrogen_g;
                                        detailed_canopy_state.?.plant_seed_storage_phosphorus_g[plant] = initialized_seed_storage.phosphorus_g;
                                    }
                                }
                                plant_phenology_state.?.reseed_pending[plant] = false;
                                plant_phenology_state.?.death_replant_pending[plant] = false;
                                plant_state.canopy_water_storage_m_per_m2[plant] = 0;
                                plant_state.leaf_area_index_m2_m2[plant] = 0;
                                plant_state.shoot_carbon_g_m2[plant] = 0;
                                @memset(plant_state.root_carbon_g_m2[plant * config.soil_layers .. (plant + 1) * config.soil_layers], 0);
                                const retention_state = &canopy_precipitation_retention_state.?;
                                inline for (@typeInfo(ecosys.canopy_precipitation_retention.State).@"struct".fields) |field| {
                                    if (field.type == []f64 and !std.mem.startsWith(u8, field.name, "cell_")) @field(retention_state, field.name)[plant] = 0;
                                }
                                if (preserve_reseed_inventories)
                                    retention_state.standing_dead_surface_water_m3[plant] = carried_standing_dead_surface_water_m3;
                            }
                            if (topology_changed) {
                                var replacement_layers = try ecosys.canopy_layer_distribution.State.init(allocator, state.cell_count, config.plant_populations, runscript.canopy_layer_count, canopy_geometry.leaf_inclination_sine.len, canopy_geometry.leaf_azimuth_radians.len, &detailed_canopy_state.?);
                                const replacement_carbon_exchange = ecosys.canopy_carbon_exchange.State.init(allocator, detailed_canopy_state.?.branch_node_offsets.len - 1) catch |err| {
                                    replacement_layers.deinit();
                                    return err;
                                };
                                var previous_carbon_exchange = canopy_carbon_exchange_state.?;
                                canopy_carbon_exchange_state.? = replacement_carbon_exchange;
                                previous_carbon_exchange.deinit();
                                var previous_layers = canopy_layer_distribution_state.?;
                                canopy_layer_distribution_state.? = replacement_layers;
                                previous_layers.deinit();
                            }
                            if (harvest_context != null) for (0..runtime_plant_count) |plant| {
                                const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                                try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                            };
                        }
                    }
                }
                if (fertilizer_schedule_maps[pass.scene_index] != null) {
                    const management_date = try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp);
                    var fertilizer_context: ecosys.fertilizer_management_dispatch.NitrogenApplyContext = .{
                        .soil = &soil_fertilizer_inventory,
                        .surface = &surface_litter_fertilizer_state,
                        .cell_area_m2 = canopy_cell_area_m2,
                        .active_soil_layer_count = state.active_soil_layer_count,
                        .soil_layer_thickness_m = soil_solver_property_state.layer_thickness_m,
                        .nitrogen_molar_mass_g_per_mol = runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                        .surface_organic = &surface_organic_state,
                        .source_hour_one_through_twenty_four = timestamp.hour,
                        .solar_noon_hour_by_cell = fertilizer_solar_noon_hour_by_cell,
                    };
                    _ = try ecosys.fertilizer_management_dispatch.dispatchDate(fertilizer_schedule_maps[pass.scene_index].?, fertilizer_catalog, management_date, &fertilizer_context, ecosys.fertilizer_management_dispatch.applyNitrogen);
                    var mineral_fertilizer_context: ecosys.fertilizer_management_dispatch.MineralApplyContext = .{
                        .inventory = &mineral_fertilizer_inventory,
                        .surface_organic = &surface_organic_state,
                        .cell_area_m2 = canopy_cell_area_m2,
                        .active_soil_layer_count = state.active_soil_layer_count,
                        .soil_layer_thickness_m = soil_solver_property_state.layer_thickness_m,
                        .source_hour_one_through_twenty_four = timestamp.hour,
                        .solar_noon_hour_by_cell = fertilizer_solar_noon_hour_by_cell,
                    };
                    _ = try ecosys.fertilizer_management_dispatch.dispatchDate(fertilizer_schedule_maps[pass.scene_index].?, fertilizer_catalog, management_date, &mineral_fertilizer_context, ecosys.fertilizer_management_dispatch.applyMinerals);
                    var organic_fertilizer_context: ecosys.fertilizer_management_dispatch.OrganicApplyContext = .{
                        .soil = &soil_organic_state,
                        .surface = &surface_organic_state,
                        .parameters = if (organic_parameters) |*parameters| parameters else return error.MissingOrganicMatterParameters,
                        .cell_area_m2 = canopy_cell_area_m2,
                        .active_soil_layer_count = state.active_soil_layer_count,
                        .soil_layer_capacity = config.soil_layers,
                        .soil_layer_thickness_m = soil_solver_property_state.layer_thickness_m,
                        .daily_organic_carbon_input_g_c = daily_organic_fertilizer_carbon_input_g_c,
                        .daily_biome_carbon_input_g_c = daily_biome_organic_carbon_input_g_c,
                        .daily_organic_phosphorus_input_g_p = daily_organic_fertilizer_phosphorus_input_g_p,
                        .daily_organic_nitrogen_input_g_n = daily_organic_fertilizer_nitrogen_input_g_n,
                        .source_hour_one_through_twenty_four = timestamp.hour,
                        .solar_noon_hour_by_cell = fertilizer_solar_noon_hour_by_cell,
                    };
                    _ = try ecosys.fertilizer_management_dispatch.dispatchDate(fertilizer_schedule_maps[pass.scene_index].?, fertilizer_catalog, management_date, &organic_fertilizer_context, ecosys.fertilizer_management_dispatch.applyOrganic);
                }
                // GROSUB applies tillage at integer solar noon, not at the
                // daily management initialization boundary.
                if (@as(f64, @floatFromInt(timestamp.hour)) == @floor(weather_stream.header.solar_noon_hour) and disturbance_schedule_maps[pass.scene_index] != null) {
                    const management_date = try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp);
                    var disturbance_context: ecosys.disturbance_management_dispatch.ApplyContext = .{
                        .roots = &plant_root_state.?,
                        .litter_partition = &plant_litter_partition_state.?,
                        .soil_organic = &soil_organic_state,
                        .grid = &state,
                        .species_count = config.plant_populations,
                        .biological_domain_count_by_plant = root_biological_domain_count_by_plant,
                        .root_nonwoody_fraction_by_plant = plant_water_workspace.?.woody_root_fraction,
                        .biomass_turnover_type_by_plant = canopy_layer_controls.biomass_turnover_type,
                        .root_profile_type_by_plant = canopy_layer_controls.root_profile_type,
                        .growth_habit_by_plant = plant_topology_controls.growth_habit_code,
                        .leaf_phenology_type_by_plant = plant_topology_controls.leaf_phenology_code,
                        .planting_day_of_year_by_plant = development_planting_day_of_year,
                        .planting_year_by_plant = development_planting_year,
                        .current_day_of_year = current_day_of_year,
                        .current_year = @intCast(current_year),
                        .plant_harvest = if (harvest_context) |*value| value else null,
                        .root_litter_carbon_ledger = &root_litter_carbon_ledger,
                        .surface_energy = &surface_energy_state,
                        .fire_active_this_hour = fire_active_this_hour,
                    };
                    _ = try ecosys.disturbance_management_dispatch.dispatchDatePhase(
                        disturbance_schedule_maps[pass.scene_index].?,
                        disturbance_catalog,
                        management_date,
                        &disturbance_context,
                        .pre_science,
                    );
                }
                if (begins_day and active_options.climate_change_mode == 2) {
                    try climate_state.advanceDay(active_options, ecosys.climate_change.daysInYear(timestamp.year.?));
                }
                for (0..state.cell_count) |cell| {
                    const stream_index = weather_assignment.stream_index_by_cell[cell];
                    if (stream_index >= weather_hour_by_stream.len)
                        return error.WeatherGridStreamIndexOutOfRange;
                    const selected_site = site_by_cell[cell];
                    const cell_weather_stream = weather_streams[stream_index];
                    hourly_weather_header_by_cell[cell] = cell_weather_stream.header;
                    var cell_radiation = try ecosys.atmospheric_radiation.prepare(
                        weather_hour_by_stream[stream_index].forcing,
                        timestamp,
                        selected_site.latitude_degrees_north,
                        cell_weather_stream.header.solar_noon_hour,
                        selected_site.ecosystem_type == -2,
                    );
                    hourly_forcing_by_cell[cell] = try climate_state.apply(
                        cell_radiation.forcing,
                        active_options,
                        timestamp,
                        cell_weather_stream.header.solar_noon_hour,
                        selected_site.elevation_m,
                    );
                    cell_radiation.forcing = hourly_forcing_by_cell[cell];
                    hourly_radiation_by_cell[cell] = cell_radiation;
                }
                current_atmospheric_co2_umol_per_mol = site.atmospheric_co2_umol_mol * try climate_state.atmosphericCo2Multiplier(active_options, current_day_of_year);
                if (!std.math.isFinite(current_atmospheric_co2_umol_per_mol) or current_atmospheric_co2_umol_per_mol <= 0) return error.InvalidAtmosphericCarbonDioxide;
                if (surface_gas_parameters) |parameters| {
                    current_atmospheric_gas_concentration_g_per_m3.* = parameters.atmospheric_concentration_g_per_m3;
                    current_atmospheric_gas_concentration_g_per_m3[@intFromEnum(ecosys.gas_transport.Species.carbon_dioxide)] *= current_atmospheric_co2_umol_per_mol / site.atmospheric_co2_umol_mol;
                }
                const management_date = try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp);
                if (begins_day and irrigation_schedule_maps[pass.scene_index] != null) {
                    @memset(minimum_canopy_water_potential_mpa_by_cell, 0);
                    try ecosys.plant_water_balance.publishFirstPlantDailyMinimumByCell(&detailed_canopy_state.?, minimum_canopy_water_potential_mpa_by_cell);
                    _ = try ecosys.irrigation_management_dispatch.planAutomatedDay(
                        irrigation_schedule_maps[pass.scene_index].?,
                        irrigation_catalog,
                        .{
                            .date = management_date,
                            .cell_area_m2 = canopy_cell_area_m2,
                            .active_layer_count = state.active_soil_layer_count,
                            .soil_layer_capacity = state.soil_layer_capacity,
                            .layer_thickness_m = soil_solver_property_state.layer_thickness_m,
                            .porosity_fraction = soil_solver_property_state.porosity_fraction,
                            .field_capacity_fraction = soil_field_capacity_fraction,
                            .wilting_point_fraction = soil_wilting_point_fraction,
                            .liquid_water_m3 = state.liquid_water_m3,
                            .ice_water_m3 = state.ice_water_m3,
                            .minimum_canopy_water_potential_mpa = minimum_canopy_water_potential_mpa_by_cell,
                        },
                        automatic_irrigation_depth_m_by_cell_hour,
                    );
                }
                if (begins_day and detailed_canopy_state != null) ecosys.plant_water_balance.resetDailyMinimumCanopyWaterPotential(&detailed_canopy_state.?);
                @memset(irrigation_water_depth_m, 0);
                @memset(irrigation_dissolved_mass_g_per_m2, 0);
                @memset(irrigation_hydrogen_mol_per_m2, 0);
                irrigation_loads.reset();
                const irrigation_source_hour: u8 = if (timestamp.hour == 0) 1 else timestamp.hour;
                if (irrigation_schedule_maps[pass.scene_index]) |irrigation_map| {
                    _ = try ecosys.irrigation_management_dispatch.accumulateScheduledHourRouted(
                        irrigation_map,
                        irrigation_catalog,
                        management_date,
                        irrigation_source_hour,
                        .{
                            .cell_area_m2 = canopy_cell_area_m2,
                            .active_layer_count = state.active_soil_layer_count,
                            .layer_thickness_m = soil_solver_property_state.layer_thickness_m,
                            .loads = &irrigation_loads,
                        },
                    );
                    _ = try ecosys.irrigation_management_dispatch.accumulatePlannedAutomatedHourRouted(
                        irrigation_map,
                        irrigation_catalog,
                        irrigation_source_hour,
                        automatic_irrigation_depth_m_by_cell_hour,
                        .{
                            .cell_area_m2 = canopy_cell_area_m2,
                            .active_layer_count = state.active_soil_layer_count,
                            .layer_thickness_m = soil_solver_property_state.layer_thickness_m,
                            .loads = &irrigation_loads,
                        },
                    );
                }
                for (0..state.cell_count) |cell| {
                    const area_m2 = canopy_cell_area_m2[cell];
                    if (!std.math.isFinite(area_m2) or area_m2 <= 0)
                        return error.InvalidIrrigationCellArea;
                    irrigation_water_depth_m[cell] =
                        irrigation_loads.surface_water_m3[cell] / area_m2;
                    irrigation_hydrogen_mol_per_m2[cell] =
                        irrigation_loads.surface_hydrogen_mol[cell] / area_m2;
                    const first = cell *
                        ecosys.irrigation_management_dispatch.dissolved_species_count;
                    for (0..ecosys.irrigation_management_dispatch.dissolved_species_count) |species|
                        irrigation_dissolved_mass_g_per_m2[first + species] =
                            irrigation_loads.surface_dissolved_mass_g[first + species] /
                            area_m2;
                }
                const previous_day_of_year: u16 = if (current_day_of_year > 1) current_day_of_year - 1 else ecosys.climate_change.daysInYear(current_year -| 1);
                for (site_by_cell, plant_calendar_by_cell) |cell_site, *calendar| {
                    const phytotron = cell_site.ecosystem_type == -2;
                    const current_daylength_h = if (phytotron) 24 else ecosys.daily_weather_disaggregation.calculateDaylength(current_day_of_year, cell_site.latitude_degrees_north);
                    const previous_daylength_h = if (phytotron) 24 else ecosys.daily_weather_disaggregation.calculateDaylength(previous_day_of_year, cell_site.latitude_degrees_north);
                    const solstice_day: u16 = if (cell_site.latitude_degrees_north >= 0) 173 else 355;
                    const maximum_daylength_h = if (phytotron) 24 else ecosys.daily_weather_disaggregation.calculateDaylength(solstice_day, cell_site.latitude_degrees_north);
                    calendar.* = .{
                        .day_of_year = current_day_of_year,
                        .current_year = current_year,
                        .current_daylength_h = current_daylength_h,
                        .previous_daylength_h = previous_daylength_h,
                        .maximum_seasonal_daylength_h = maximum_daylength_h,
                        .latitude_degrees_north = cell_site.latitude_degrees_north,
                    };
                }
                organic_matter_fire_exchange_state.resetHourly();
                surface_fire_exchange_state.resetHourly();
                @memset(shoot_harvest_litter_carbon_g_c_by_plant, 0);
                @memset(shoot_harvest_litter_nitrogen_g_n_by_plant, 0);
                @memset(shoot_harvest_litter_phosphorus_g_p_by_plant, 0);
                @memset(hourly_manure_products_by_plant, .{});
                @memcpy(
                    harvest_carbon_at_hour_start_g_c_by_plant,
                    plant_daily_flux_ledger.harvested_carbon_g,
                );
                const gas_failure_file_name = try std.fmt.allocPrint(
                    allocator,
                    "ecosys-ng-gas-failure-ex{d}-scenario{d}-repeat{d}-scene{d}-year{d}-day{d}-hour{d}-cells0-{d}.bin",
                    .{
                        pass.execution_iteration + 1,
                        pass.scenario_index + 1,
                        pass.scenario_iteration + 1,
                        pass.scene_index + 1,
                        current_year,
                        current_day_of_year,
                        timestamp.hour,
                        state.cell_count - 1,
                    },
                );
                defer allocator.free(gas_failure_file_name);
                const gas_failure_path = try std.fs.path.join(
                    allocator,
                    &.{ runscript_directory, gas_failure_file_name },
                );
                defer allocator.free(gas_failure_path);
                const solute_failure_file_name = try std.fmt.allocPrint(
                    allocator,
                    "ecosys-ng-solute-failure-ex{d}-scenario{d}-repeat{d}-scene{d}-year{d}-day{d}-hour{d}.bin",
                    .{
                        pass.execution_iteration + 1,
                        pass.scenario_index + 1,
                        pass.scenario_iteration + 1,
                        pass.scene_index + 1,
                        current_year,
                        current_day_of_year,
                        timestamp.hour,
                    },
                );
                defer allocator.free(solute_failure_file_name);
                const solute_failure_path = try std.fs.path.join(
                    allocator,
                    &.{ runscript_directory, solute_failure_file_name },
                );
                defer allocator.free(solute_failure_path);
                executeHourlyScience(
                    hourly_science_context,
                    timestamp.hour,
                    hourly_radiation_by_cell,
                    hourly_forcing_by_cell,
                    hourly_weather_header_by_cell,
                    plant_calendar_by_cell,
                    .{ .value = executed_weather_hours + 1 },
                    .{
                        .io = init.io,
                        .directory = std.Io.Dir.cwd(),
                        .file_path = gas_failure_path,
                    },
                    .{
                        .io = init.io,
                        .directory = std.Io.Dir.cwd(),
                        .file_path = solute_failure_path,
                        .context = .{
                            .execution_id = @intCast(
                                pass.execution_iteration + 1,
                            ),
                            .scenario_id = @intCast(
                                pass.scenario_index + 1,
                            ),
                            .repeat_id = @intCast(
                                pass.scenario_iteration + 1,
                            ),
                            .scene_id = @intCast(pass.scene_index + 1),
                            .scene_hour = @intCast(scene_weather_hours + 1),
                            .year = @intCast(current_year),
                            .day_of_year = @intCast(current_day_of_year),
                            .hour = @intCast(timestamp.hour),
                        },
                    },
                ) catch |err| {
                    std.log.err("hourly science failed: execution={d} scenario={d} scenario_repeat={d} scene={d} scene_hour={d} total_hour={d} year={?} day_of_year={?} month={?} day={?} hour={d} error={s}", .{ pass.execution_iteration + 1, pass.scenario_index + 1, pass.scenario_iteration + 1, pass.scene_index + 1, scene_weather_hours + 1, executed_weather_hours + 1, timestamp.year, timestamp.day_of_year, timestamp.month, timestamp.day_of_month, timestamp.hour, @errorName(err) });
                    return err;
                };
                if (harvest_context != null and plant_phenology_state != null) {
                    var transitioned = false;
                    for (plant_phenology_state.?.leafout_transition_this_step, 0..) |leafout, plant| {
                        if (!leafout or !plant_phenology_state.?.active[plant] or plant_topology_controls.growth_habit_code[plant] == 0) continue;
                        const branches = try plant_growth_stage_state.?.branchRange(plant);
                        const fully_deciduous = canopy_layer_controls.biomass_turnover_type[plant] == 0;
                        for (branches.first..branches.end) |branch| {
                            branch_development_state.?.remobilization_progress_h[branch] = 0;
                            branch_development_state.?.leafout_initialization_enabled[branch] = false;
                            try ecosys.plant_growth_stages.resetBranchForSeasonalLeafout(
                                &plant_growth_stage_state.?.branches[branch],
                                current_day_of_year,
                                fully_deciduous,
                                branch_development_state.?.initial_reproductive_stage[branch],
                            );
                        }
                        if (fully_deciduous and branches.first < branches.end) {
                            plant_phenology_state.?.initiated_node_count[plant] = branch_development_state.?.initial_reproductive_stage[branches.first];
                            plant_phenology_state.?.appeared_leaf_count[plant] = 0;
                        }
                        try ecosys.plant_harvest_runtime.applyStartOfSeasonResidue(
                            &harvest_context.?,
                            plant,
                            canopy_layer_controls.biomass_turnover_type[plant],
                            canopy_layer_controls.root_profile_type[plant],
                        );
                        transitioned = true;
                    }
                    if (transitioned) for (0..runtime_plant_count) |plant| {
                        const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                        try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                    };
                }
                // GROSUB creates an immediate grain/JHVST=2 event when the
                // main branch of a self-seeding deciduous annual enters
                // end-of-season reproductive turnover.
                if (harvest_context != null and plant_dormancy_state != null) {
                    const automatic_self_seed_count = try ecosys.plant_harvest_runtime.applyAutomaticSelfSeedingHarvests(
                        &harvest_context.?,
                        &plant_dormancy_state.?,
                        plant_topology_controls.growth_habit_code,
                        plant_topology_controls.leaf_phenology_code,
                    );
                    if (automatic_self_seed_count > 0) for (0..runtime_plant_count) |plant| {
                        const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                        try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                    };
                }
                // Tillage is a late-GROSUB transaction. Applying it here keeps
                // current-hour growth and FWODR composition, and prevents the
                // removed roots from regrowing before the next hour.
                if (@as(f64, @floatFromInt(timestamp.hour)) == @floor(weather_stream.header.solar_noon_hour) and disturbance_schedule_maps[pass.scene_index] != null) {
                    var disturbance_context: ecosys.disturbance_management_dispatch.ApplyContext = .{
                        .roots = &plant_root_state.?,
                        .litter_partition = &plant_litter_partition_state.?,
                        .soil_organic = &soil_organic_state,
                        .grid = &state,
                        .species_count = config.plant_populations,
                        .biological_domain_count_by_plant = root_biological_domain_count_by_plant,
                        .root_nonwoody_fraction_by_plant = plant_water_workspace.?.woody_root_fraction,
                        .biomass_turnover_type_by_plant = canopy_layer_controls.biomass_turnover_type,
                        .root_profile_type_by_plant = canopy_layer_controls.root_profile_type,
                        .growth_habit_by_plant = plant_topology_controls.growth_habit_code,
                        .leaf_phenology_type_by_plant = plant_topology_controls.leaf_phenology_code,
                        .planting_day_of_year_by_plant = development_planting_day_of_year,
                        .planting_year_by_plant = development_planting_year,
                        .current_day_of_year = current_day_of_year,
                        .current_year = @intCast(current_year),
                        .plant_harvest = if (harvest_context) |*value| value else null,
                        .surface_energy = &surface_energy_state,
                        .fire_active_this_hour = fire_active_this_hour,
                    };
                    _ = try ecosys.disturbance_management_dispatch.dispatchDatePhase(
                        disturbance_schedule_maps[pass.scene_index].?,
                        disturbance_catalog,
                        try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp),
                        &disturbance_context,
                        .post_science,
                    );
                    if (harvest_context != null) for (0..runtime_plant_count) |plant| {
                        const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                        try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                    };
                }
                if (harvest_context != null and management_schedule_map != null and plant_phenology_state != null) {
                    @memset(grazing_average_shoot_carbon_g_c, 0);
                    @memset(grazing_active_cell_count, 0);
                    const canopy = &detailed_canopy_state.?;
                    for (0..runtime_plant_count) |plant| {
                        if (!plant_phenology_state.?.active[plant]) continue;
                        const species = plant % config.plant_populations;
                        grazing_average_shoot_carbon_g_c[species] += canopy.plant_total_shoot_carbon_g[plant];
                        grazing_active_cell_count[species] += 1;
                    }
                    for (grazing_average_shoot_carbon_g_c, grazing_active_cell_count) |*average, count| {
                        if (count > 0) average.* /= @floatFromInt(count);
                    }
                    try ecosys.plant_management_dispatch.accumulateGrazingRebalanceTotals(
                        management_schedule_map.?,
                        plant_management_catalog,
                        management_date,
                        canopy.plant_total_shoot_carbon_g,
                        plant_phenology_state.?.active,
                        grazing_total_shoot_carbon_g_c_by_cell,
                        grazing_total_grazer_biomass_by_cell,
                    );
                    const GrazingDispatch = struct {
                        harvest: *ecosys.plant_harvest_runtime.Context,
                        average_shoot_carbon_g_c: []const f64,
                        cell_area_m2: []const f64,
                        active: []const bool,
                        turnover_type: []const u8,
                        root_profile_type: []const u8,
                        total_shoot_carbon_g_c_by_cell: []const f64,
                        total_grazer_biomass_by_cell: []const f64,
                        species_count: usize,
                        monthly_self_thinning_hour: bool,

                        fn apply(self: *@This(), plant: usize, event: ecosys.plant_management.HarvestEvent) !void {
                            if (!self.active[plant]) return;
                            if (self.monthly_self_thinning_hour and self.turnover_type[plant] != 0 and self.root_profile_type[plant] > 1) return;
                            const species = plant % self.species_count;
                            const cell = plant / self.species_count;
                            var effective_event = event;
                            if (event.kind == .animal_grazing and self.total_shoot_carbon_g_c_by_cell[cell] > 0) {
                                effective_event.cutting_height_m_or_lai_fraction =
                                    self.total_grazer_biomass_by_cell[cell] *
                                    self.harvest.canopy_state.plant_total_shoot_carbon_g[plant] /
                                    self.total_shoot_carbon_g_c_by_cell[cell];
                            }
                            _ = try ecosys.plant_harvest_runtime.applyGrazingEvent(
                                self.harvest,
                                plant,
                                effective_event,
                                self.average_shoot_carbon_g_c[species],
                                self.cell_area_m2[cell],
                            );
                        }
                    };
                    var grazing_dispatch: GrazingDispatch = .{
                        .harvest = &harvest_context.?,
                        .average_shoot_carbon_g_c = grazing_average_shoot_carbon_g_c,
                        .cell_area_m2 = canopy_cell_area_m2,
                        .active = plant_phenology_state.?.active,
                        .turnover_type = canopy_layer_controls.biomass_turnover_type,
                        .root_profile_type = canopy_layer_controls.root_profile_type,
                        .total_shoot_carbon_g_c_by_cell = grazing_total_shoot_carbon_g_c_by_cell,
                        .total_grazer_biomass_by_cell = grazing_total_grazer_biomass_by_cell,
                        .species_count = config.plant_populations,
                        .monthly_self_thinning_hour = current_day_of_year % 30 == 0 and
                            @as(f64, @floatFromInt(timestamp.hour)) == @floor(weather_stream.header.solar_noon_hour),
                    };
                    _ = try ecosys.plant_management_dispatch.dispatchDateGrazing(
                        management_schedule_map.?,
                        plant_management_catalog,
                        management_date,
                        &grazing_dispatch,
                        GrazingDispatch.apply,
                    );
                    for (0..runtime_plant_count) |plant| {
                        const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                        try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                    }
                }
                // GROSUB monthly forest self-thinning is evaluated at local
                // solar noon after the current-hour stem diameter refresh.
                if (harvest_context != null and current_day_of_year % 30 == 0 and
                    @as(f64, @floatFromInt(timestamp.hour)) == @floor(weather_stream.header.solar_noon_hour))
                {
                    for (0..runtime_plant_count) |plant| {
                        if (canopy_layer_controls.biomass_turnover_type[plant] == 0 or
                            canopy_layer_controls.root_profile_type[plant] <= 1 or
                            plant_phenology_state == null or !plant_phenology_state.?.active[plant]) continue;
                        if (management_schedule_map != null and
                            !try ecosys.plant_management_dispatch.allowsForestSelfThinning(
                                management_schedule_map.?,
                                plant_management_catalog,
                                try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp),
                                plant,
                            )) continue;
                        const thinned_fraction = try ecosys.plant_harvest_runtime.applyForestSelfThinning(&harvest_context.?, plant);
                        if (thinned_fraction > 0) {
                            const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                            try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                        }
                    }
                }
                if (harvest_context != null and management_schedule_map != null and
                    @as(f64, @floatFromInt(timestamp.hour)) == @floor(weather_stream.header.solar_noon_hour))
                {
                    const harvest_date = try ecosys.plant_management_dispatch.dateFromTimestamp(timestamp);
                    _ = try ecosys.plant_management_dispatch.dispatchDateNonGrazing(
                        management_schedule_map.?,
                        plant_management_catalog,
                        harvest_date,
                        &harvest_context.?,
                        ecosys.plant_harvest_runtime.applyEvent,
                    );
                    for (harvest_products_by_plant.?, 0..) |_, plant| {
                        const exported = try ecosys.plant_harvest_runtime.publishPlantProducts(&harvest_context.?, plant);
                        try plant_daily_flux_ledger.accumulateHarvest(plant, exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g);
                    }
                }
                try daily_soil_gas_flux.accumulateHour(
                    if (plant_root_state) |*roots| roots else null,
                    config.plant_populations,
                    config.soil_layers,
                    soil_gas_transport_state.atmospheric_flux_g_per_h,
                    surface_litter_gas_transport_state.atmospheric_flux_g_per_h,
                );
                for (0..state.cell_count) |cell| {
                    var evaporation_m3: f64 = @max(
                        0,
                        ground_air_surface_vapor_conductance_m3_per_h[cell] *
                            (ground_air_surface_vapor_fraction[cell] - ground_air_state.vapor_volume_fraction[cell]),
                    );
                    for (0..config.plant_populations) |species| {
                        const plant = cell * config.plant_populations + species;
                        if (canopy_surface_exchange_state) |exchange| {
                            evaporation_m3 += @max(0, -exchange.transpiration_m3_per_h[plant]);
                            evaporation_m3 += @max(0, -exchange.intercepted_water_change_m3_per_h[plant]);
                        }
                        evaporation_m3 += @max(0, -standing_dead_evaporation_m3_per_h[plant]);
                    }
                    var external_water_outflow_m3: f64 = 0;
                    for (0..config.soil_layers) |local_layer| {
                        const layer = try state.layerIndex(cell, local_layer);
                        external_water_outflow_m3 += @max(0, transport_hydrology_state.micropore_external_water_flux_m3_per_step[layer]);
                        external_water_outflow_m3 += @max(0, transport_hydrology_state.macropore_external_water_flux_m3_per_step[layer]);
                    }
                    const top_layer = try state.layerIndex(cell, 0);
                    const bulk_density_Mg_per_m3 = soil_solver_property_state.bulk_density_megagrams_per_m3[top_layer];
                    if (!std.math.isFinite(bulk_density_Mg_per_m3) or bulk_density_Mg_per_m3 <= 0) return error.InvalidOutputSoilBulkDensity;
                    try daily_water_ledger.accumulateCell(cell, .{
                        .rainfall_m3 = atmospheric_state.precipitation_m[cell] * canopy_cell_area_m2[cell],
                        .evaporation_m3 = evaporation_m3,
                        .runoff_m3 = @max(0, -surface_runoff_state.exported_water_m3[cell]),
                        .water_outflow_m3 = external_water_outflow_m3,
                        .lateral_water_outflow_m3 = transport_hydrology_state.artificial_drainage_outflow_m3_per_step[cell],
                        .sediment_outflow_m3 = surface_erosion_state.routing.sediment_export_Mg[cell] / bulk_density_Mg_per_m3,
                    });
                    var ionic_outflow_mol: f64 = 0;
                    const active_layers = state.active_soil_layer_count[cell];
                    const first_layer = try state.layerIndex(cell, 0);
                    for (0..active_layers) |local_layer| {
                        const base = (first_layer + local_layer) * aqueous_species_count;
                        for (soil_solute_boundary_net_flux_mol[base .. base + aqueous_species_count]) |flux_mol| ionic_outflow_mol += @max(0, -flux_mol);
                    }
                    try daily_heat_ledger.accumulateHour(
                        cell,
                        active_layers,
                        atmospheric_state.shortwave_radiation_mj_per_m2[cell],
                        atmospheric_state.air_temperature_k[cell],
                        atmospheric_state.vapor_pressure_kpa[cell],
                        atmospheric_state.wind_speed_m_per_h[cell],
                        atmospheric_state.precipitation_m[cell],
                        state.soil_temperature_k[first_layer .. first_layer + active_layers],
                        state.surface_temperature_k[cell],
                        ionic_outflow_mol,
                    );
                }
                try daily_heterotrophic_respiration.accumulateHour(&state);
                try daily_carbon_export.accumulateOrganicDrainageHour(&state, &soil_organic_transport_state);
                try daily_carbon_export.accumulateOrganicRunoffHour(surface_organic_carbon_export_g_c_per_h);
                try daily_carbon_export.accumulateInorganicRunoffHour(surface_inorganic_carbon_export_g_c_per_h);
                try daily_carbon_export.accumulateDissolvedInorganicDrainageHour(&state, &soil_dissolved_gas_transport_state);
                try daily_carbon_export.accumulateGaseousInorganicDrainageHour(&state, &soil_gas_transport_state);
                try daily_phosphorus_export.accumulateOrganicDrainageHour(&state, &soil_organic_transport_state);
                try daily_phosphorus_export.accumulateInorganicDrainageHour(
                    &state,
                    soil_solute_boundary_net_flux_mol,
                    runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                );
                try daily_phosphorus_export.accumulateRunoffHour(
                    surface_organic_phosphorus_export_g_p_per_h,
                    surface_inorganic_phosphorus_export_g_p_per_h,
                );
                try daily_nitrogen_export.accumulateOrganicDrainageHour(&state, &soil_organic_transport_state);
                try daily_nitrogen_export.accumulateInorganicDrainageHour(&state, &mineral_nitrogen_transport_state);
                try daily_nitrogen_export.accumulateRunoffHour(
                    surface_organic_nitrogen_export_g_n_per_h,
                    surface_inorganic_nitrogen_export_g_n_per_h,
                );
                const soil_process_units = soil_nitrogen_flux_workspace.process_unit_count_per_layer;
                const surface_process_units = surface_microbial_mineral_exchange_state.h2po4_exchange_g_p.len / state.cell_count;
                for (0..state.cell_count) |cell| {
                    var immobilization_g_p: f64 = 0;
                    var immobilization_g_n: f64 = 0;
                    for (0..state.active_soil_layer_count[cell]) |local_layer| {
                        const layer = try state.layerIndex(cell, local_layer);
                        const first_unit = layer * soil_process_units;
                        for (first_unit..first_unit + soil_process_units) |unit| immobilization_g_p +=
                            soil_nitrogen_flux_workspace.non_band_microbial_h2po4_exchange_g_p[unit] +
                            soil_nitrogen_flux_workspace.band_microbial_h2po4_exchange_g_p[unit] +
                            soil_nitrogen_flux_workspace.non_band_microbial_hpo4_exchange_g_p[unit] +
                            soil_nitrogen_flux_workspace.band_microbial_hpo4_exchange_g_p[unit];
                        for (first_unit..first_unit + soil_process_units) |unit| immobilization_g_n +=
                            soil_nitrogen_flux_workspace.non_band_microbial_ammonium_exchange_g_n[unit] +
                            soil_nitrogen_flux_workspace.band_microbial_ammonium_exchange_g_n[unit] +
                            soil_nitrogen_flux_workspace.non_band_microbial_nitrate_exchange_g_n[unit] +
                            soil_nitrogen_flux_workspace.band_microbial_nitrate_exchange_g_n[unit];
                    }
                    const surface_first = cell * surface_process_units;
                    for (surface_first..surface_first + surface_process_units) |unit| immobilization_g_p +=
                        surface_microbial_mineral_exchange_state.h2po4_exchange_g_p[unit] +
                        surface_microbial_mineral_exchange_state.hpo4_exchange_g_p[unit];
                    for (surface_first..surface_first + surface_process_units) |unit| immobilization_g_n +=
                        surface_microbial_mineral_exchange_state.ammonium_exchange_g_n[unit] +
                        surface_microbial_mineral_exchange_state.nitrate_exchange_g_n[unit];
                    const next = daily_microbial_phosphate_mineralization_g_p[cell] - immobilization_g_p;
                    if (!std.math.isFinite(next)) return error.NonFiniteDailyMicrobialPhosphateMineralization;
                    daily_microbial_phosphate_mineralization_g_p[cell] = next;
                    const nitrogen_next = daily_microbial_nitrogen_mineralization_g_n[cell] - immobilization_g_n;
                    if (!std.math.isFinite(nitrogen_next)) return error.NonFiniteDailyMicrobialNitrogenMineralization;
                    daily_microbial_nitrogen_mineralization_g_n[cell] = nitrogen_next;
                }
                // Late GROSUB salt removal follows all natural turnover,
                // management harvest, grazing, and tillage publications.
                // This is the single serial pre-EXTRACT ownership boundary.
                if (runscript.dynamic_plant_salts) {
                    const canopy = if (detailed_canopy_state) |*value| value else return error.MissingCanopyForPlantSaltHarvest;
                    const roots = if (plant_root_state) |*value| value else return error.MissingPlantRootsForPlantSaltHarvest;
                    for (0..runtime_plant_count) |plant| {
                        const active_domain_count = root_metabolism_plant_parameters[plant].biologicalDomainCount();
                        try root_litter_carbon_ledger.carbonByPlantLayer(
                            plant,
                            active_domain_count,
                            root_litter_carbon_g_c_by_layer,
                        );
                        var natural_shoot_litter_carbon_g_c: f64 = 0;
                        const natural_litter = shoot_senescence_products_by_plant[plant];
                        for (natural_litter.woody_carbon_g, natural_litter.nonwoody_carbon_g) |woody, nonwoody| {
                            natural_shoot_litter_carbon_g_c += woody + nonwoody;
                            if (!std.math.isFinite(natural_shoot_litter_carbon_g_c))
                                return error.NonFinitePlantSaltLitterCarbon;
                        }
                        try ecosys.plant_salt_harvest_adapter.applyPlant(
                            &plant_salt_harvest_workspace,
                            canopy,
                            roots,
                            .{
                                .cumulative_harvest_salt_mol_by_plant = cumulative_harvest_salt_mol_by_plant,
                                .litter_salt_mol_by_plant_exchange_layer = litter_salt_mol_by_plant_exchange_layer,
                            },
                            plant,
                            .{
                                .shoot_carbon_g_c = canopy.plant_total_shoot_carbon_g[plant],
                                .current_harvest_carbon_g_c = plant_daily_flux_ledger.harvested_carbon_g[plant],
                                .previous_harvest_carbon_g_c = harvest_carbon_at_hour_start_g_c_by_plant[plant],
                                .shoot_litterfall_carbon_g_c = natural_shoot_litter_carbon_g_c +
                                    shoot_harvest_litter_carbon_g_c_by_plant[plant],
                                .root_litterfall_carbon_g_c_by_layer = root_litter_carbon_g_c_by_layer,
                                .active_root_domain_count = active_domain_count,
                                .minimum_plant_mass_g_c = config.absolute_tolerance,
                            },
                        );
                    }
                    try ecosys.litter_salt_publication.refresh(
                        &litter_salt_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .mol_by_plant_exchange_layer_salt = litter_salt_mol_by_plant_exchange_layer,
                        },
                    );
                    try ecosys.plant_litter_salt_ingress.apply(
                        &plant_litter_salt_ingress_state,
                        .{
                            .species_count = config.plant_populations,
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .litter_salt_mol_by_plant_exchange_layer = litter_salt_mol_by_plant_exchange_layer,
                            .surface_water_m3 = surface_precipitation_state.litter_water_m3,
                            .soil_water_m3 = state.matrix_liquid_water_m3,
                            .minimum_water_m3 = config.absolute_tolerance,
                        },
                        &surface_litter_chemistry_state,
                        &initial_chemistry_state,
                    );
                }
                try ecosys.manure_deposition_publication.refresh(
                    &manure_deposition_publication_state,
                    hourly_manure_products_by_plant,
                );
                // EXTRACT accumulates these source-signed hourly exchanges
                // for OUTPD before DAY carries balances and clears them.
                @memset(root_soil_element_exchange_workspace, 0);
                if (detailed_canopy_state) |*canopy| {
                    const carbon_exchange = if (canopy_carbon_exchange_state) |*value| value else return error.PlantDailyLedgerRequiresCarbonExchangeState;
                    try ecosys.canopy_ammonia_publication.refresh(
                        &canopy_ammonia_publication_state,
                        .{
                            .plant_branch_offsets = canopy.plant_branch_offsets,
                            .branch_exchange_g_n_per_h = canopy.branch_canopy_ammonia_exchange_g_n_per_h,
                            .preceding_cumulative_exchange_g_n_by_plant = plant_daily_flux_ledger.ammonia_exchange_g_n,
                        },
                    );
                    for (0..runtime_plant_count) |plant| {
                        const carbon_flux = try carbon_exchange.fluxForPlant(canopy, plant);
                        const living_water_source_m3 = if (canopy_surface_exchange_state) |exchange|
                            exchange.transpiration_m3_per_h[plant] + exchange.intercepted_water_change_m3_per_h[plant]
                        else
                            0;
                        var signed_total_respiration_carbon_g = carbon_flux.signed_aboveground_respiration_g_c_per_h;
                        if (plant_root_state) |roots| {
                            for (0..root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| for (0..config.soil_layers) |layer| {
                                const root = try roots.layerIndex(plant, domain, layer);
                                signed_total_respiration_carbon_g -= roots.actual_respiration_g_c_per_h[root] +
                                    roots.symbiotic_respiration_actual_g_c_per_h[root];
                            };
                        }
                        try plant_daily_flux_ledger.accumulateHourlyExchange(plant, .{
                            .net_canopy_carbon_g = carbon_flux.net_co2_g_c_per_h,
                            .gross_primary_productivity_g = canopy.plant_gross_primary_productivity_g_c_per_h[plant],
                            .signed_total_respiration_carbon_g = signed_total_respiration_carbon_g,
                            .signed_aboveground_respiration_carbon_g = carbon_flux.signed_aboveground_respiration_g_c_per_h,
                            .canopy_and_standing_dead_water_source_m3 = living_water_source_m3 + standing_dead_evaporation_m3_per_h[plant],
                            .canopy_ammonia_exchange_g_n = 0,
                        });
                        if (plant_root_state) |*roots| {
                            var root_soil_carbon_exchange_g: f64 = 0;
                            var root_soil_nitrogen_exchange_g: f64 = 0;
                            var root_soil_phosphorus_exchange_g: f64 = 0;
                            for (0..root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| for (0..config.soil_layers) |layer| {
                                const root = try roots.layerIndex(plant, domain, layer);
                                for (0..ecosys.plant_root_system.organic_substrate_count) |fraction| {
                                    const substrate =
                                        root * ecosys.plant_root_system.organic_substrate_count +
                                        fraction;
                                    root_soil_carbon_exchange_g +=
                                        roots.exudate_carbon_exchange_g_c_per_h[substrate];
                                    root_soil_nitrogen_exchange_g +=
                                        roots.exudate_nitrogen_exchange_g_n_per_h[substrate];
                                    root_soil_phosphorus_exchange_g +=
                                        roots.exudate_phosphorus_exchange_g_p_per_h[substrate];
                                }
                            };
                            root_soil_nitrogen_exchange_g +=
                                roots.ammonium_uptake_g_n_per_h[plant] +
                                roots.nitrate_uptake_g_n_per_h[plant];
                            root_soil_phosphorus_exchange_g += roots.phosphate_uptake_g_p_per_h[plant];
                            root_soil_carbon_exchange_g_c_per_h_by_plant[plant] =
                                root_soil_carbon_exchange_g;
                            root_soil_nitrogen_exchange_g_n_per_h_by_plant[plant] =
                                root_soil_nitrogen_exchange_g;
                            root_soil_phosphorus_exchange_g_p_per_h_by_plant[plant] =
                                root_soil_phosphorus_exchange_g;
                            try plant_daily_flux_ledger.accumulateHourlyRootSoilExchange(
                                plant,
                                root_soil_carbon_exchange_g,
                                root_soil_nitrogen_exchange_g,
                                root_soil_phosphorus_exchange_g,
                                roots.fixation_uptake_g_n_per_h[plant],
                            );
                        }
                        const litterfall = try ecosys.plant_litterfall_publication.totals(
                            shoot_senescence_products_by_plant[plant],
                            root_litter_products_by_plant[plant],
                            .{
                                .carbon_g_c = shoot_harvest_litter_carbon_g_c_by_plant[plant],
                                .nitrogen_g_n = shoot_harvest_litter_nitrogen_g_n_by_plant[plant],
                                .phosphorus_g_p = shoot_harvest_litter_phosphorus_g_p_by_plant[plant],
                            },
                        );
                        aboveground_litter_carbon_g_c_per_h_by_plant[plant] =
                            litterfall.aboveground.carbon_g_c;
                        aboveground_litter_nitrogen_g_n_per_h_by_plant[plant] =
                            litterfall.aboveground.nitrogen_g_n;
                        aboveground_litter_phosphorus_g_p_per_h_by_plant[plant] =
                            litterfall.aboveground.phosphorus_g_p;
                        try plant_daily_flux_ledger.accumulateLitterPublication(
                            plant,
                            litterfall.aboveground.carbon_g_c,
                            litterfall.aboveground.nitrogen_g_n,
                            litterfall.aboveground.phosphorus_g_p,
                            litterfall.belowground.carbon_g_c,
                            litterfall.belowground.nitrogen_g_n,
                            litterfall.belowground.phosphorus_g_p,
                        );
                    }
                    try ecosys.cell_litter_standing_dead_publication.refresh(
                        &cell_litter_standing_dead_publication_state,
                        .{
                            .aboveground_litter_carbon_g_c_per_h_by_plant = aboveground_litter_carbon_g_c_per_h_by_plant,
                            .aboveground_litter_nitrogen_g_n_per_h_by_plant = aboveground_litter_nitrogen_g_n_per_h_by_plant,
                            .aboveground_litter_phosphorus_g_p_per_h_by_plant = aboveground_litter_phosphorus_g_p_per_h_by_plant,
                            .standing_dead_carbon_g_c_by_plant = canopy.plant_standing_dead_carbon_g,
                        },
                    );
                    const phenology = if (plant_phenology_state) |*value|
                        value
                    else
                        return error.RootSoilElementPublicationRequiresPhenologyState;
                    try ecosys.root_soil_element_exchange_publication.refresh(
                        &root_soil_element_exchange_publication_state,
                        .{
                            .active_by_plant = phenology.active,
                            .exchange_g_element_per_h_by_element_and_plant = .{
                                root_soil_carbon_exchange_g_c_per_h_by_plant,
                                root_soil_nitrogen_exchange_g_n_per_h_by_plant,
                                root_soil_phosphorus_exchange_g_p_per_h_by_plant,
                            },
                        },
                    );
                    @memcpy(
                        plant_daily_flux_ledger.ammonia_exchange_g_n,
                        canopy_ammonia_publication_state
                            .cumulative_exchange_g_n_by_plant,
                    );
                } else {
                    inline for (
                        root_soil_element_exchange_publication_state
                            .change_g_element_per_h_by_element_and_cell,
                    ) |values| @memset(values, 0);
                    @memset(
                        canopy_ammonia_publication_state
                            .current_exchange_g_n_per_h_by_plant,
                        0,
                    );
                    @memcpy(
                        canopy_ammonia_publication_state
                            .cumulative_exchange_g_n_by_plant,
                        plant_daily_flux_ledger.ammonia_exchange_g_n,
                    );
                    inline for (@typeInfo(
                        ecosys.cell_litter_standing_dead_publication.State,
                    ).@"struct".fields) |field|
                        if (field.type == []f64)
                            @memset(
                                @field(
                                    cell_litter_standing_dead_publication_state,
                                    field.name,
                                ),
                                0,
                            );
                }
                @memset(canopy_net_radiation_mj, 0);
                @memset(canopy_storage_heat_mj, 0);
                if (canopy_energy_state) |energy| if (canopy_surface_exchange_state) |exchange| {
                    for (0..runtime_plant_count) |plant| {
                        const cell = plant / config.plant_populations;
                        canopy_net_radiation_mj[plant] = energy.net_radiation_mj_per_m2[plant] * canopy_cell_area_m2[cell];
                        // HFLXC is the converged canopy energy residual sum;
                        // EXTRACT later subtracts its convective VFLXC term.
                        canopy_storage_heat_mj[plant] = canopy_net_radiation_mj[plant] +
                            exchange.latent_heat_flux_mj_per_h[plant] +
                            exchange.sensible_heat_flux_mj_per_h[plant] +
                            exchange.vapor_sensible_heat_flux_mj_per_h[plant];
                    }
                };
                const living_exchange = if (canopy_surface_exchange_state) |*exchange| exchange else null;
                const dead_exchange = if (standing_dead_surface_exchange_state) |*exchange| exchange else null;
                try ecosys.plant_energy_publication.refresh(
                    &plant_energy_publication_state,
                    .{
                        .living_net_radiation_mj = canopy_net_radiation_mj,
                        .living_latent_heat_mj = if (living_exchange) |exchange| exchange.latent_heat_flux_mj_per_h else zero_plant_energy_mj,
                        .living_sensible_heat_mj = if (living_exchange) |exchange| exchange.sensible_heat_flux_mj_per_h else zero_plant_energy_mj,
                        .living_storage_heat_mj = canopy_storage_heat_mj,
                        .living_convective_water_heat_mj = if (living_exchange) |exchange| exchange.vapor_sensible_heat_flux_mj_per_h else zero_plant_energy_mj,
                        .standing_dead_net_radiation_mj = if (dead_exchange) |exchange| exchange.net_radiation_mj_per_h else zero_plant_energy_mj,
                        .standing_dead_latent_heat_mj = if (dead_exchange) |exchange| exchange.latent_heat_flux_mj_per_h else zero_plant_energy_mj,
                        .standing_dead_sensible_heat_mj = if (dead_exchange) |exchange| exchange.sensible_heat_flux_mj_per_h else zero_plant_energy_mj,
                        .standing_dead_storage_heat_mj = if (dead_exchange) |exchange| exchange.storage_heat_flux_mj_per_h else zero_plant_energy_mj,
                        .standing_dead_convective_water_heat_mj = if (dead_exchange) |exchange| exchange.vapor_sensible_heat_flux_mj_per_h else zero_plant_energy_mj,
                    },
                );
                if (canopy_precipitation_retention_state) |*retention| {
                    try ecosys.plant_water_publication.refresh(
                        &plant_water_publication_state,
                        .{
                            .cell_area_m2 = canopy_cell_area_m2,
                            .internal_water_depth_m_per_m2_by_plant = plant_state.canopy_water_storage_m_per_m2,
                            .living_surface_water_m3_by_plant = retention.living_surface_water_m3,
                            .standing_dead_surface_water_m3_by_plant = retention.standing_dead_surface_water_m3,
                            .transpiration_m3_per_h_by_plant = if (living_exchange) |exchange|
                                exchange.transpiration_m3_per_h
                            else
                                zero_plant_energy_mj,
                            .living_evaporation_m3_per_h_by_plant = if (living_exchange) |exchange|
                                exchange.intercepted_water_change_m3_per_h
                            else
                                zero_plant_energy_mj,
                            .standing_dead_evaporation_m3_per_h_by_plant = if (dead_exchange) |exchange|
                                exchange.intercepted_water_change_m3_per_h
                            else
                                zero_plant_energy_mj,
                        },
                    );
                }
                try ecosys.ecosystem_energy_ledger.refresh(&ecosystem_energy_ledger_state, .{
                    .cell_area_m2 = canopy_cell_area_m2,
                    .ground_net_radiation_mj_per_m2 = surface_energy_state.net_radiation_mj_per_m2,
                    .ground_latent_heat_mj_per_m2 = surface_temperature_solver_state.latent_heat_flux_mj_per_m2,
                    .ground_sensible_heat_mj_per_m2 = surface_temperature_solver_state.sensible_heat_flux_mj_per_m2,
                    .ground_storage_heat_mj_per_m2 = surface_temperature_solver_state.storage_heat_flux_mj_per_m2,
                    .species_count = config.plant_populations,
                    .canopy_net_radiation_mj = canopy_net_radiation_mj,
                    .canopy_latent_heat_mj = if (living_exchange) |exchange| exchange.latent_heat_flux_mj_per_h else zero_plant_energy_mj,
                    .canopy_sensible_heat_mj = if (living_exchange) |exchange| exchange.sensible_heat_flux_mj_per_h else zero_plant_energy_mj,
                    .canopy_storage_heat_mj = canopy_storage_heat_mj,
                    .canopy_convective_water_heat_mj = if (living_exchange) |exchange| exchange.vapor_sensible_heat_flux_mj_per_h else zero_plant_energy_mj,
                    .standing_dead_net_radiation_mj = if (dead_exchange) |exchange| exchange.net_radiation_mj_per_h else zero_plant_energy_mj,
                    .standing_dead_latent_heat_mj = if (dead_exchange) |exchange| exchange.latent_heat_flux_mj_per_h else zero_plant_energy_mj,
                    .standing_dead_sensible_heat_mj = if (dead_exchange) |exchange| exchange.sensible_heat_flux_mj_per_h else zero_plant_energy_mj,
                    .standing_dead_storage_heat_mj = if (dead_exchange) |exchange| exchange.storage_heat_flux_mj_per_h else zero_plant_energy_mj,
                    .standing_dead_convective_water_heat_mj = if (dead_exchange) |exchange| exchange.vapor_sensible_heat_flux_mj_per_h else zero_plant_energy_mj,
                });
                if (canopy_precipitation_retention_state) |*retention| {
                    try ecosys.canopy_water_energy_publication.refresh(
                        &canopy_water_energy_publication_state,
                        .{
                            .air_temperature_k_by_cell = atmospheric_state.air_temperature_k,
                            .canopy_temperature_k_by_plant = plant_state.canopy_temperature_k,
                            .living_surface_water_m3_by_plant = retention.living_surface_water_m3,
                            .standing_dead_surface_water_m3_by_plant = retention.standing_dead_surface_water_m3,
                            .living_retention_m3_per_h_by_plant = retention.living_retention_m3_per_h,
                            .standing_dead_retention_m3_per_h_by_plant = retention.standing_dead_retention_m3_per_h,
                        },
                        retention.previous_water_energy_mj,
                    );
                    @memcpy(
                        ecosystem_energy_ledger_state.canopy_water_energy_mj,
                        canopy_water_energy_publication_state
                            .water_energy_mj_by_cell,
                    );
                    @memcpy(
                        ecosystem_energy_ledger_state.canopy_water_energy_change_mj_per_h,
                        canopy_water_energy_publication_state
                            .water_energy_change_mj_per_h_by_cell,
                    );
                } else {
                    @memset(
                        canopy_water_energy_publication_state
                            .water_energy_mj_by_plant,
                        0,
                    );
                    @memset(
                        canopy_water_energy_publication_state
                            .water_energy_change_mj_per_h_by_plant,
                        0,
                    );
                    @memset(
                        canopy_water_energy_publication_state
                            .water_energy_mj_by_cell,
                        0,
                    );
                    @memset(
                        canopy_water_energy_publication_state
                            .water_energy_change_mj_per_h_by_cell,
                        0,
                    );
                    @memset(ecosystem_energy_ledger_state.canopy_water_energy_mj, 0);
                    @memset(ecosystem_energy_ledger_state.canopy_water_energy_change_mj_per_h, 0);
                }
                if (canopy_layer_distribution_state) |*layers| {
                    const canopy = if (detailed_canopy_state) |*value|
                        value
                    else
                        return error.LivingCanopyPublicationRequiresCanopyState;
                    const phenology = if (plant_phenology_state) |*value|
                        value
                    else
                        return error.LivingCanopyPublicationRequiresPhenologyState;
                    try ecosys.living_canopy_layer_publication.refresh(
                        &living_canopy_layer_publication_state,
                        .{
                            .active_by_plant = phenology.active,
                            .plant_branch_offsets = canopy.plant_branch_offsets,
                            .branch_node_offsets = canopy.branch_node_offsets,
                            .node_leaf_area_m2_by_layer = layers.node_leaf_area_m2,
                            .node_leaf_carbon_g_c_by_layer = layers.node_leaf_carbon_g,
                            .branch_stalk_area_m2_by_layer = layers.branch_stalk_area_m2,
                        },
                    );
                    try ecosys.standing_dead_area_publication.refresh(
                        &standing_dead_area_publication_state,
                        layers.plant_standing_dead_area_m2,
                    );
                } else {
                    @memset(
                        living_canopy_layer_publication_state
                            .leaf_area_m2_by_cell_layer,
                        0,
                    );
                    @memset(
                        living_canopy_layer_publication_state
                            .leaf_carbon_g_c_by_cell_layer,
                        0,
                    );
                    @memset(
                        living_canopy_layer_publication_state
                            .stalk_area_m2_by_cell_layer,
                        0,
                    );
                    @memset(
                        living_canopy_layer_publication_state
                            .leaf_area_m2_by_cell,
                        0,
                    );
                    @memset(
                        living_canopy_layer_publication_state
                            .stalk_area_m2_by_cell,
                        0,
                    );
                    @memset(
                        standing_dead_area_publication_state
                            .area_m2_by_cell_layer,
                        0,
                    );
                    @memset(
                        standing_dead_area_publication_state.area_m2_by_cell,
                        0,
                    );
                }
                try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedSurfaceAndCanopyHeat(
                    surface_energy_state.net_radiation_mj_per_m2,
                    surface_temperature_solver_state.sensible_heat_flux_mj_per_m2,
                    surface_temperature_solver_state.latent_heat_flux_mj_per_m2,
                    surface_temperature_solver_state.vapor_sensible_heat_flux_mj_per_m2,
                    canopy_cell_area_m2,
                    ecosystem_energy_ledger_state.canopy_water_energy_change_mj_per_h,
                );
                if (plant_root_state) |*roots| {
                    const water_workspace = if (plant_water_workspace) |*workspace| workspace else return error.RootUptakeLedgerRequiresPlantWaterWorkspace;
                    try ecosys.root_uptake_ledger.refresh(
                        &root_uptake_ledger_state,
                        &state,
                        roots,
                        config.plant_populations,
                        root_biological_domain_count_by_plant,
                        water_workspace.plant_population_count,
                        canopy_cell_area_m2,
                        runscript.dynamic_plant_salts,
                    );
                    const phenology = if (plant_phenology_state) |*value|
                        value
                    else
                        return error.RootWaterPublicationRequiresPhenologyState;
                    try ecosys.root_water_uptake_publication.refresh(
                        &root_water_uptake_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .plant_population_count = water_workspace.plant_population_count,
                            .cell_area_m2 = canopy_cell_area_m2,
                            .soil_temperature_k_by_layer = state.soil_temperature_k,
                            .root_length_density_m_per_m3 = roots.root_length_density_m_per_m3,
                            .water_uptake_m3_per_h = roots.water_uptake_m3_per_h,
                        },
                    );
                    @memcpy(
                        root_uptake_ledger_state
                            .total_root_length_density_m_per_m3,
                        root_water_uptake_publication_state
                            .root_length_density_m_per_m3,
                    );
                    @memcpy(
                        root_uptake_ledger_state.total_water_uptake_m3_per_h,
                        root_water_uptake_publication_state
                            .water_uptake_m3_per_h,
                    );
                    @memcpy(
                        root_uptake_ledger_state.convective_water_heat_mj_per_h,
                        root_water_uptake_publication_state
                            .convective_water_heat_mj_per_h,
                    );
                    try ecosys.root_gas_content_publication.refresh(
                        &root_gas_content_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .gaseous_g_by_gas = .{
                                roots.gaseous_carbon_dioxide_g_c,
                                roots.gaseous_oxygen_g_o,
                                roots.gaseous_methane_g_c,
                                roots.gaseous_nitrous_oxide_g_n,
                                roots.gaseous_ammonia_g_n,
                                roots.gaseous_hydrogen_g_h,
                            },
                            .aqueous_g_by_gas = .{
                                roots.aqueous_carbon_dioxide_g_c,
                                roots.aqueous_oxygen_g_o,
                                roots.aqueous_methane_g_c,
                                roots.aqueous_nitrous_oxide_g_n,
                                roots.aqueous_ammonia_g_n,
                                roots.aqueous_hydrogen_g_h,
                            },
                        },
                    );
                    inline for (
                        root_uptake_ledger_state.total_root_gas_content_g,
                        root_gas_content_publication_state
                            .total_g_by_gas_and_layer,
                    ) |destination, source_values| @memcpy(destination, source_values);
                    try ecosys.root_atmosphere_gas_publication.refresh(
                        &root_atmosphere_gas_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .exchange_g_per_h_by_root_and_transport_gas = roots.atmosphere_to_root_gas_exchange_g_per_h,
                        },
                    );
                    inline for (
                        root_uptake_ledger_state
                            .total_atmosphere_to_root_gas_exchange_g_per_h,
                        root_atmosphere_gas_publication_state
                            .exchange_g_per_h_by_gas_and_layer,
                    ) |destination, source_values| @memcpy(destination, source_values);
                    try ecosys.root_internal_gas_publication.refresh(
                        &root_internal_gas_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .aqueous_carbon_dioxide_reaction_g_c_per_h_by_root = roots.aqueous_carbon_dioxide_reaction_g_c_per_h,
                            .oxygen_uptake_from_root_pool_g_o_per_h_by_root = roots.oxygen_uptake_from_root_pool_g_o_per_h,
                        },
                    );
                    @memcpy(
                        root_uptake_ledger_state
                            .root_pool_oxygen_uptake_g_o_per_h,
                        root_internal_gas_publication_state
                            .oxygen_uptake_g_o_per_h_by_layer,
                    );
                    try ecosys.root_soil_gas_publication.refresh(
                        &root_soil_gas_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .soil_to_root_exchange_g_per_h_by_root_and_transport_gas = roots.soil_to_root_gas_exchange_g_per_h,
                            .oxygen_uptake_from_soil_g_o_per_h_by_root = roots.oxygen_uptake_from_soil_g_o_per_h,
                        },
                    );
                    for (
                        root_uptake_ledger_state
                            .total_soil_to_root_gas_exchange_g_per_h[0..4],
                        root_soil_gas_publication_state
                            .exchange_g_per_h_by_gas_and_layer[0..4],
                    ) |destination, source_values| @memcpy(destination, source_values);
                    @memcpy(
                        root_uptake_ledger_state
                            .total_soil_to_root_gas_exchange_g_per_h[5],
                        root_soil_gas_publication_state
                            .exchange_g_per_h_by_gas_and_layer[4],
                    );
                    @memcpy(
                        root_uptake_ledger_state.soil_oxygen_uptake_g_o_per_h,
                        root_soil_gas_publication_state
                            .exchange_g_per_h_by_gas_and_layer[1],
                    );
                    try ecosys.root_nutrient_uptake_publication.refresh(
                        &root_nutrient_uptake_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .uptake_g_element_per_h_by_nutrient_and_root = .{
                                roots.ammonium_uptake_nonband_g_n_per_h,
                                roots.nitrate_uptake_nonband_g_n_per_h,
                                roots.phosphate_h2_uptake_nonband_g_p_per_h,
                                roots.phosphate_h_uptake_nonband_g_p_per_h,
                                roots.ammonium_uptake_band_g_n_per_h,
                                roots.nitrate_uptake_band_g_n_per_h,
                                roots.phosphate_h2_uptake_band_g_p_per_h,
                                roots.phosphate_h_uptake_band_g_p_per_h,
                            },
                        },
                    );
                    for (
                        root_uptake_ledger_state.total_nutrient_uptake_by_layer,
                        root_nutrient_uptake_publication_state
                            .uptake_g_element_per_h_by_nutrient_and_layer,
                    ) |destination, source_values| @memcpy(
                        destination,
                        source_values,
                    );
                    try ecosys.root_salt_uptake_publication.refresh(
                        &root_salt_uptake_publication_state,
                        .{
                            .dynamic_salts_enabled = runscript.dynamic_plant_salts,
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .uptake_mol_per_h_by_root_and_salt = roots.salt_uptake_mol_per_h,
                        },
                    );
                    @memcpy(
                        root_uptake_ledger_state.total_salt_uptake_mol_per_h,
                        root_salt_uptake_publication_state
                            .uptake_mol_per_h_by_layer_and_salt,
                    );
                    try ecosys.root_exudate_publication.refresh(
                        &root_exudate_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .exchange_g_element_per_h_by_element_root_and_fraction = .{
                                roots.exudate_carbon_exchange_g_c_per_h,
                                roots.exudate_nitrogen_exchange_g_n_per_h,
                                roots.exudate_phosphorus_exchange_g_p_per_h,
                            },
                        },
                    );
                    @memcpy(
                        root_uptake_ledger_state
                            .total_exudate_carbon_change_g_c_per_h,
                        root_exudate_publication_state
                            .change_g_element_per_h_by_element_layer_and_fraction[0],
                    );
                    @memcpy(
                        root_uptake_ledger_state
                            .total_exudate_nitrogen_change_g_n_per_h,
                        root_exudate_publication_state
                            .change_g_element_per_h_by_element_layer_and_fraction[1],
                    );
                    @memcpy(
                        root_uptake_ledger_state
                            .total_exudate_phosphorus_change_g_p_per_h,
                        root_exudate_publication_state
                            .change_g_element_per_h_by_element_layer_and_fraction[2],
                    );
                    try ecosys.root_competition_demand_publication.refresh(
                        &root_competition_demand_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .demand_g_element_per_h_by_kind_and_root = .{
                                roots.oxygen_demand_g_o_per_h,
                                roots.ammonium_demand_nonband_g_n_per_h,
                                roots.nitrate_demand_nonband_g_n_per_h,
                                roots.phosphate_h2_demand_nonband_g_p_per_h,
                                roots.phosphate_h_demand_nonband_g_p_per_h,
                                roots.ammonium_demand_band_g_n_per_h,
                                roots.nitrate_demand_band_g_n_per_h,
                                roots.phosphate_h2_demand_band_g_p_per_h,
                                roots.phosphate_h_demand_band_g_p_per_h,
                            },
                        },
                    );
                    for (
                        root_uptake_ledger_state
                            .total_competition_demand_by_layer,
                        root_competition_demand_publication_state
                            .demand_g_element_per_h_by_kind_and_layer,
                    ) |destination, source_values| @memcpy(
                        destination,
                        source_values,
                    );
                    try ecosys.root_nitrogen_fixation_publication.refresh(
                        &root_nitrogen_fixation_publication_state,
                        .{
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .active_by_plant = phenology.active,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .fixation_g_n_per_h_by_root = roots.fixation_uptake_g_n_per_h_by_layer,
                        },
                    );
                    @memcpy(
                        root_uptake_ledger_state
                            .total_nitrogen_fixation_g_n_per_h,
                        root_nitrogen_fixation_publication_state
                            .fixation_g_n_per_h_by_layer,
                    );
                } else {
                    @memset(
                        root_water_uptake_publication_state
                            .root_length_density_m_per_m3,
                        0,
                    );
                    @memset(
                        root_water_uptake_publication_state
                            .water_uptake_m3_per_h,
                        0,
                    );
                    @memset(
                        root_water_uptake_publication_state
                            .convective_water_heat_mj_per_h,
                        0,
                    );
                    inline for (
                        root_gas_content_publication_state
                            .total_g_by_gas_and_layer,
                    ) |values| @memset(values, 0);
                    inline for (
                        root_atmosphere_gas_publication_state
                            .exchange_g_per_h_by_gas_and_layer,
                    ) |values| @memset(values, 0);
                    @memset(
                        root_internal_gas_publication_state
                            .carbon_dioxide_source_g_c_per_h_by_layer,
                        0,
                    );
                    inline for (
                        root_soil_gas_publication_state
                            .exchange_g_per_h_by_gas_and_layer,
                    ) |values| @memset(values, 0);
                    inline for (
                        root_nutrient_uptake_publication_state
                            .uptake_g_element_per_h_by_nutrient_and_layer,
                    ) |values| @memset(values, 0);
                    @memset(
                        root_salt_uptake_publication_state
                            .uptake_mol_per_h_by_layer_and_salt,
                        0,
                    );
                    inline for (
                        root_exudate_publication_state
                            .change_g_element_per_h_by_element_layer_and_fraction,
                    ) |values| @memset(values, 0);
                    inline for (
                        root_competition_demand_publication_state
                            .demand_g_element_per_h_by_kind_and_layer,
                    ) |values| @memset(values, 0);
                    @memset(
                        root_nitrogen_fixation_publication_state
                            .fixation_g_n_per_h_by_layer,
                        0,
                    );
                    @memset(
                        root_internal_gas_publication_state
                            .oxygen_uptake_g_o_per_h_by_layer,
                        0,
                    );
                    @memset(root_uptake_ledger_state.total_root_length_density_m_per_m3, 0);
                    @memset(root_uptake_ledger_state.total_water_uptake_m3_per_h, 0);
                    @memset(root_uptake_ledger_state.convective_water_heat_mj_per_h, 0);
                    inline for (root_uptake_ledger_state.total_root_gas_content_g) |values| @memset(values, 0);
                    inline for (root_uptake_ledger_state.total_soil_to_root_gas_exchange_g_per_h) |values| @memset(values, 0);
                    inline for (root_uptake_ledger_state.total_aqueous_to_gaseous_root_exchange_g_per_h) |values| @memset(values, 0);
                    inline for (root_uptake_ledger_state.total_atmosphere_to_root_gas_exchange_g_per_h) |values| @memset(values, 0);
                    @memset(root_uptake_ledger_state.total_salt_uptake_mol_per_h, 0);
                    @memset(root_uptake_ledger_state.total_exudate_carbon_change_g_c_per_h, 0);
                    @memset(root_uptake_ledger_state.total_exudate_nitrogen_change_g_n_per_h, 0);
                    @memset(root_uptake_ledger_state.total_exudate_phosphorus_change_g_p_per_h, 0);
                    inline for (root_uptake_ledger_state.total_competition_demand_by_layer) |values| @memset(values, 0);
                    inline for (root_uptake_ledger_state.total_nutrient_uptake_by_layer) |values| @memset(values, 0);
                    @memset(root_uptake_ledger_state.soil_oxygen_uptake_g_o_per_h, 0);
                    @memset(root_uptake_ledger_state.root_pool_oxygen_uptake_g_o_per_h, 0);
                    @memset(root_uptake_ledger_state.total_nitrogen_fixation_g_n_per_h, 0);
                }
                var soil_fire = runscript.soil_fire_combustion_parameters;
                soil_fire.negligible_carbon_g_c = config.absolute_tolerance;
                const fire_product_parameters: ecosys.plant_soil_exchange.SubsurfaceFireParameters = .{
                    .oxygen_g_per_g_carbon = soil_fire.oxygen_g_per_g_combusted_carbon,
                    .maximum_aerobic_charcoal_fraction = soil_fire.maximum_aerobic_charcoal_fraction,
                    .maximum_anaerobic_charcoal_fraction = soil_fire.maximum_anaerobic_charcoal_fraction,
                    .oxygen_half_saturation_g_o_per_m3 = soil_fire.oxygen_half_saturation_g_o_per_m3,
                    .methane_half_saturation_g_c_per_m3 = soil_fire.methane_half_saturation_g_c_per_m3,
                    .aerobic_combustion_energy_mj_per_g_carbon = soil_fire.aerobic_combustion_energy_mj_per_g_c,
                    .anaerobic_combustion_energy_mj_per_g_carbon = soil_fire.anaerobic_combustion_energy_mj_per_g_c,
                    .methane_combustion_energy_mj_per_g_carbon = soil_fire.methane_combustion_energy_mj_per_g_c,
                };
                for (0..state.cell_count) |cell| {
                    _ = ecosys.soil_combustion.burnSurfaceOrganicStateCell(
                        &surface_organic_state,
                        &surface_fire_exchange_state,
                        cell,
                        fire_active_this_hour[cell],
                        state.surface_temperature_k[cell],
                        canopy_cell_area_m2[cell],
                        1,
                        soil_fire,
                    ) catch |err| {
                        std.log.err("surface organic matter combustion failed: scene={d} scene_hour={d} cell={d} error={s}", .{ pass.scene_index + 1, scene_weather_hours + 1, cell, @errorName(err) });
                        return err;
                    };
                    surface_fire_exchange_state.finalizeSurfaceCell(
                        cell,
                        &litter_gas_transport_state,
                        &surface_litter_chemistry_state,
                        &surface_organic_state,
                        surface_precipitation_state.litter_water_m3,
                        runscript.dynamic_plant_salts,
                        delayed_surface_combustion_heat_mj,
                        config.absolute_tolerance,
                        runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                        runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                        fire_product_parameters,
                    ) catch |err| {
                        std.log.err("surface combustion products failed: scene={d} scene_hour={d} cell={d} error={s}", .{ pass.scene_index + 1, scene_weather_hours + 1, cell, @errorName(err) });
                        return err;
                    };
                }
                // GROSUB manure and fire mineral products enter the surface
                // solution above. Re-converge only this local SOLUTE kernel;
                // never repeat the full hourly or legacy sub-hourly cycle.
                try convergeSurfaceLitterChemistry(hourly_science_context);
                for (0..state.cell_count) |cell| {
                    if (!fire_active_this_hour[cell]) continue;
                    for (0..state.active_soil_layer_count[cell]) |local_layer| {
                        const layer = try state.layerIndex(cell, local_layer);
                        _ = ecosys.soil_combustion.burnOrganicStateLayer(
                            &soil_organic_state,
                            &organic_matter_fire_exchange_state,
                            layer,
                            true,
                            state.soil_temperature_k[layer],
                            canopy_cell_area_m2[cell],
                            1,
                            soil_fire,
                        ) catch |err| {
                            std.log.err("soil organic matter combustion failed: scene={d} scene_hour={d} cell={d} layer={d} error={s}", .{ pass.scene_index + 1, scene_weather_hours + 1, cell, local_layer, @errorName(err) });
                            return err;
                        };
                    }
                }
                if (plant_root_state) |*roots| {
                    ecosys.disturbance_management_dispatch.applyRootFireCombustion(
                        roots,
                        if (detailed_canopy_state) |*canopy| canopy else null,
                        &state,
                        config.plant_populations,
                        root_biological_domain_count_by_plant,
                        canopy_cell_area_m2,
                        fire_active_this_hour,
                        1,
                        runscript.dynamic_plant_salts,
                        runscript.plant_fire_combustion_parameters,
                        &organic_matter_fire_exchange_state,
                    ) catch |err| {
                        std.log.err("root nodule combustion failed: scene={d} scene_hour={d} error={s}", .{ pass.scene_index + 1, scene_weather_hours + 1, @errorName(err) });
                        return err;
                    };
                    try ecosys.root_combustion_salt_publication.refresh(
                        &root_combustion_salt_publication_state,
                        .{
                            .dynamic_salts_enabled = runscript.dynamic_plant_salts,
                            .fire_active_by_cell = fire_active_this_hour,
                            .active_soil_layer_count_by_cell = state.active_soil_layer_count,
                            .root_domain_count_by_plant = root_biological_domain_count_by_plant,
                            .released_mol_per_h_by_root_and_salt = roots.combustion_salt_loss_mol_per_h,
                        },
                    );
                } else {
                    @memset(
                        root_combustion_salt_publication_state
                            .released_mol_per_h_by_layer_and_salt,
                        0,
                    );
                    @memset(
                        root_combustion_salt_publication_state
                            .total_released_mol_per_h_by_cell,
                        0,
                    );
                }
                // Soil organic fire products exist independently of whether
                // the run has plant roots; roots merely add to this ledger.
                for (0..state.cell_count) |cell| {
                    if (!fire_active_this_hour[cell]) continue;
                    for (0..state.active_soil_layer_count[cell]) |local_layer| {
                        const layer = try state.layerIndex(cell, local_layer);
                        organic_matter_fire_exchange_state.finalizeLayer(
                            layer,
                            &gas_transport_state,
                            &plant_available_nutrient_state,
                            &soil_organic_state,
                            &micropore_solute_state,
                            runscript.dynamic_plant_salts,
                            delayed_subsurface_combustion_heat_mj,
                            config.absolute_tolerance,
                            fire_product_parameters,
                        ) catch |err| {
                            std.log.err("subsurface combustion products failed: scene={d} scene_hour={d} cell={d} layer={d} error={s}", .{ pass.scene_index + 1, scene_weather_hours + 1, cell, local_layer, @errorName(err) });
                            return err;
                        };
                    }
                    var emitted_carbon_dioxide_g_c =
                        surface_fire_exchange_state.carbon_dioxide_emission_g_c[cell];
                    var emitted_methane_g_c =
                        surface_fire_exchange_state.methane_emission_g_c[cell];
                    var produced_charcoal_g_c =
                        surface_fire_exchange_state.charcoal_production_g_c[cell];
                    var emitted_phosphorus_g_p = surface_fire_exchange_state.gaseous_phosphorus_emission_g_p[cell];
                    var combusted_phosphorus_g_p = surface_fire_exchange_state.combusted_phosphorus_g_p[cell];
                    var emitted_nitrogen_g_n = surface_fire_exchange_state.gaseous_nitrogen_emission_g_n[cell];
                    var combusted_nitrogen_g_n = surface_fire_exchange_state.combusted_nitrogen_g_n[cell];
                    for (0..state.active_soil_layer_count[cell]) |local_layer| {
                        const layer = try state.layerIndex(cell, local_layer);
                        emitted_carbon_dioxide_g_c +=
                            organic_matter_fire_exchange_state.carbon_dioxide_emission_g_c[layer];
                        emitted_methane_g_c +=
                            organic_matter_fire_exchange_state.methane_emission_g_c[layer];
                        produced_charcoal_g_c +=
                            organic_matter_fire_exchange_state.charcoal_production_g_c[layer];
                        emitted_phosphorus_g_p += organic_matter_fire_exchange_state.gaseous_phosphorus_emission_g_p[layer];
                        combusted_phosphorus_g_p += organic_matter_fire_exchange_state.combusted_phosphorus_g_p[layer];
                        emitted_nitrogen_g_n += organic_matter_fire_exchange_state.gaseous_nitrogen_emission_g_n[layer];
                        combusted_nitrogen_g_n += organic_matter_fire_exchange_state.combusted_nitrogen_g_n[layer];
                    }
                    if (!std.math.isFinite(emitted_carbon_dioxide_g_c) or emitted_carbon_dioxide_g_c < 0 or
                        !std.math.isFinite(emitted_methane_g_c) or emitted_methane_g_c < 0)
                        return error.InvalidSoilFireCarbonEmission;
                    if (!std.math.isFinite(produced_charcoal_g_c) or produced_charcoal_g_c < 0) return error.InvalidSoilFireCharcoalProduction;
                    if (!std.math.isFinite(emitted_phosphorus_g_p) or emitted_phosphorus_g_p < 0) return error.InvalidSoilFirePhosphorusEmission;
                    if (!std.math.isFinite(combusted_phosphorus_g_p) or combusted_phosphorus_g_p < 0) return error.InvalidSoilFirePhosphorusCombustion;
                    if (!std.math.isFinite(emitted_nitrogen_g_n) or emitted_nitrogen_g_n < 0) return error.InvalidSoilFireNitrogenEmission;
                    if (!std.math.isFinite(combusted_nitrogen_g_n) or combusted_nitrogen_g_n < 0) return error.InvalidSoilFireNitrogenCombustion;
                    // REDIST fire ledgers retain their source signs:
                    // VCO2G/VCH4G are negative emissions, while VCOXFS is
                    // positive charcoal transferred back to soil organic C.
                    daily_soil_fire_carbon_dioxide_emission_g_c[cell] -= emitted_carbon_dioxide_g_c;
                    daily_soil_fire_methane_emission_g_c[cell] -= emitted_methane_g_c;
                    daily_soil_fire_charcoal_production_g_c[cell] += produced_charcoal_g_c;
                    daily_soil_fire_phosphorus_flux_g_p[cell] -= emitted_phosphorus_g_p;
                    daily_soil_combusted_phosphorus_g_p[cell] -= combusted_phosphorus_g_p;
                    daily_soil_fire_nitrogen_flux_g_n[cell] -= emitted_nitrogen_g_n;
                    daily_soil_combusted_nitrogen_g_n[cell] -= combusted_nitrogen_g_n;
                }
                if (plant_root_state) |*roots| {
                    const any_fire_active = std.mem.indexOfScalar(bool, fire_active_this_hour, true) != null;
                    if (any_fire_active) if (detailed_canopy_state) |*canopy| {
                        const canopy_fire_gas_parameters = surface_gas_parameters orelse return error.ShootFireRequiresSurfaceGasParameters;
                        ecosys.plant_shoot_fire.apply(
                            canopy,
                            &canopy_layer_controls,
                            roots,
                            canopy_cell_area_m2,
                            fire_active_this_hour,
                            1,
                            runscript.dynamic_plant_salts,
                            runscript.plant_fire_combustion_parameters,
                            ground_air_state.temperature_k,
                            ground_air_state.air_volume_m3,
                            site.atmospheric_oxygen_umol_mol,
                            site.atmospheric_methane_umol_mol,
                            canopy_fire_gas_parameters.atmospheric_concentration_g_per_m3[@intFromEnum(ecosys.gas_transport.Species.oxygen)],
                            delayed_live_canopy_combustion_heat_mj,
                            delayed_standing_dead_combustion_heat_mj,
                            &surface_organic_state,
                        ) catch |err| {
                            std.log.err("shoot combustion failed: scene={d} scene_hour={d} error={s}", .{ pass.scene_index + 1, scene_weather_hours + 1, @errorName(err) });
                            return err;
                        };
                    };
                    if (detailed_canopy_state) |*canopy| {
                        try ecosys.plant_combustion_publication.refresh(
                            &plant_combustion_publication_state,
                            .{
                                .plant_species_per_cell = config.plant_populations,
                                .fire_active_by_cell = fire_active_this_hour,
                                .signed_carbon_loss_g_c_per_h_by_plant = canopy.plant_combustion_carbon_loss_g_c_per_h,
                                .charcoal_return_g_c_per_h_by_plant = canopy.plant_fire_charcoal_production_g_c_per_h,
                                .signed_nitrogen_loss_g_n_per_h_by_plant = canopy.plant_combustion_nitrogen_loss_g_n_per_h,
                                .signed_phosphorus_loss_g_p_per_h_by_plant = canopy.plant_combustion_phosphorus_loss_g_p_per_h,
                            },
                        );
                        for (0..runtime_plant_count) |plant| {
                            const cell = plant / config.plant_populations;
                            if (!fire_active_this_hour[cell]) continue;
                            try plant_daily_flux_ledger.accumulateOxidation(
                                plant,
                                plant_combustion_publication_state.signed_carbon_balance_g_c_per_h_by_plant[plant],
                                plant_combustion_publication_state.signed_nitrogen_loss_g_n_per_h_by_plant[plant],
                                plant_combustion_publication_state.signed_phosphorus_loss_g_p_per_h_by_plant[plant],
                            );
                        }
                    }
                }
                if (detailed_canopy_state) |*canopy| {
                    try ecosys.canopy_fire_publication.refresh(
                        &canopy_fire_publication_state,
                        .{
                            .fire_active_by_cell = fire_active_this_hour,
                            .carbon_dioxide_emission_g_c_per_h_by_plant = canopy.plant_fire_carbon_dioxide_emission_g_c_per_h,
                            .methane_emission_g_c_per_h_by_plant = canopy.plant_fire_methane_emission_g_c_per_h,
                            .oxygen_consumption_g_o_per_h_by_plant = canopy.plant_fire_oxygen_consumption_g_o_per_h,
                            .charcoal_production_g_c_per_h_by_plant = canopy.plant_fire_charcoal_production_g_c_per_h,
                            .heat_release_mj_per_h_by_plant = canopy.plant_fire_heat_release_mj_per_h,
                        },
                    );
                } else {
                    inline for (@typeInfo(
                        ecosys.canopy_fire_publication.State,
                    ).@"struct".fields) |field|
                        if (field.type == []f64)
                            @memset(
                                @field(canopy_fire_publication_state, field.name),
                                0,
                            );
                }
                if (plant_root_state) |*roots| {
                    try ecosys.root_gas_withdrawal_publication.refresh(
                        &root_gas_withdrawal_publication_state,
                        .{
                            .loss_g_element_per_h_by_gas_and_plant = .{
                                roots.withdrawal_carbon_dioxide_loss_g_c_per_h,
                                roots.withdrawal_oxygen_loss_g_o_per_h,
                                roots.withdrawal_methane_loss_g_c_per_h,
                                roots.withdrawal_nitrous_oxide_loss_g_n_per_h,
                                roots.withdrawal_ammonia_loss_g_n_per_h,
                                roots.withdrawal_hydrogen_loss_g_h_per_h,
                            },
                        },
                    );
                } else {
                    inline for (
                        root_gas_withdrawal_publication_state
                            .loss_g_element_per_h_by_gas_and_cell,
                    ) |values| @memset(values, 0);
                }
                if (detailed_canopy_state) |*canopy| {
                    const exchange = if (canopy_carbon_exchange_state) |*value|
                        value
                    else
                        return error.CanopyCarbonPublicationRequiresExchangeState;
                    const phenology = if (plant_phenology_state) |*value|
                        value
                    else
                        return error.CanopyCarbonPublicationRequiresPhenologyState;
                    for (0..runtime_plant_count) |plant|
                        canopy_net_fixation_g_c_per_h_by_plant[plant] =
                            (try exchange.fluxForPlant(canopy, plant))
                                .net_co2_g_c_per_h;
                    try ecosys.canopy_carbon_publication.refresh(
                        &canopy_carbon_publication_state,
                        .{
                            .active_by_plant = phenology.active,
                            .net_fixation_g_c_per_h_by_plant = canopy_net_fixation_g_c_per_h_by_plant,
                            .oxygen_g_o_per_g_c = runscript.root_gas_parameters
                                .oxygen_to_carbon_respiration_ratio_g_o_per_g_c,
                        },
                    );
                } else {
                    @memset(
                        canopy_carbon_publication_state
                            .hourly_net_fixation_g_c_per_h_by_cell,
                        0,
                    );
                    @memset(
                        canopy_carbon_publication_state
                            .canopy_carbon_exchange_g_c_per_h_by_cell,
                        0,
                    );
                    @memset(
                        canopy_carbon_publication_state
                            .canopy_oxygen_exchange_g_o_per_h_by_cell,
                        0,
                    );
                }
                for (0..state.cell_count) |cell| {
                    const canopy_net_fixation_g_c = if (canopy_carbon_exchange_state) |ledger|
                        try ledger.netFixationForCell(&detailed_canopy_state.?, cell)
                    else
                        0;
                    try daily_canopy_gas_exchange.accumulateHour(
                        cell,
                        canopy_net_fixation_g_c,
                        canopy_fire_publication_state
                            .carbon_dioxide_emission_g_c_per_h_by_cell[cell],
                        canopy_fire_publication_state
                            .methane_emission_g_c_per_h_by_cell[cell],
                        canopy_fire_publication_state
                            .oxygen_consumption_g_o_per_h_by_cell[cell],
                        runscript.root_gas_parameters.oxygen_to_carbon_respiration_ratio_g_o_per_g_c,
                    );
                }
                const completed_hour_count = try std.math.add(usize, scene_weather_hours, 1);
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 0;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlySoilCarbonEditor;
                        const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                        const selection = output_selection_catalog.entries.items[selection_index].selection;
                        for (0..state.cell_count) |cell| {
                            const root_gas_withdrawal = if (plant_root_state) |*roots|
                                try ecosys.plant_root_disturbance.rootGasWithdrawalForCell(roots, cell, config.plant_populations)
                            else
                                ecosys.plant_root_disturbance.CellRootGasWithdrawal{};
                            const values = try hourly_soil_carbon_bank.row(cell);
                            const carbon_dioxide_concentration = values[4 .. 4 + config.soil_layers];
                            const methane_first = 5 + config.soil_layers;
                            const methane_concentration = values[methane_first .. methane_first + config.soil_layers];
                            const oxygen_first = methane_first + config.soil_layers;
                            const oxygen_concentration = values[oxygen_first .. oxygen_first + config.soil_layers];
                            for (0..config.soil_layers) |local_layer| {
                                const layer = try state.layerIndex(cell, local_layer);
                                const water_m3 = state.liquid_water_m3[layer];
                                inline for (.{ ecosys.gas_transport.Species.carbon_dioxide, ecosys.gas_transport.Species.methane, ecosys.gas_transport.Species.oxygen }) |species| {
                                    const component = layer * ecosys.gas_transport.species_count + @intFromEnum(species);
                                    const concentration = try ecosys.soil_biogeochemistry_output.dissolvedGasConcentration(gas_transport_state.dissolved_mass_g[component], water_m3);
                                    switch (species) {
                                        .carbon_dioxide => carbon_dioxide_concentration[local_layer] = concentration,
                                        .methane => methane_concentration[local_layer] = concentration,
                                        .oxygen => oxygen_concentration[local_layer] = concentration,
                                        else => unreachable,
                                    }
                                }
                            }
                            var ground_flux_g: [ecosys.gas_transport.species_count]f64 = @splat(0);
                            inline for (.{ ecosys.gas_transport.Species.carbon_dioxide, ecosys.gas_transport.Species.methane, ecosys.gas_transport.Species.oxygen }) |species|
                                ground_flux_g[@intFromEnum(species)] = try ecosys.soil_hourly_output_binding.gasBoundaryExchangeG(
                                    soil_gas_transport_state.atmospheric_flux_g_per_h,
                                    surface_litter_gas_transport_state.atmospheric_flux_g_per_h,
                                    try state.layerIndex(cell, 0),
                                    config.soil_layers,
                                    cell,
                                    species,
                                );
                            const litter_oxygen = try ecosys.soil_biogeochemistry_output.dissolvedGasConcentration(
                                litter_gas_transport_state.dissolved_mass_g[cell * ecosys.gas_transport.species_count + @intFromEnum(ecosys.gas_transport.Species.oxygen)],
                                surface_litter_gas_transport_state.water_volume_m3[cell],
                            );
                            const canopy_net_fixation_g_c_per_h = if (canopy_carbon_exchange_state) |ledger|
                                try ledger.netFixationForCell(&detailed_canopy_state.?, cell)
                            else
                                0;
                            const root_atmospheric_co2_g_c_per_h = try ecosys.root_uptake_ledger.atmosphereToRootForCellGPerH(&root_uptake_ledger_state, &state, cell, 0);
                            const root_atmospheric_oxygen_g_o_per_h = try ecosys.root_uptake_ledger.atmosphereToRootForCellGPerH(&root_uptake_ledger_state, &state, cell, 1);
                            const root_atmospheric_methane_g_c_per_h = try ecosys.root_uptake_ledger.atmosphereToRootForCellGPerH(&root_uptake_ledger_state, &state, cell, 2);
                            try ecosys.soil_biogeochemistry_output.calculateCarbonInto(.{
                                .carbon_dioxide_emission_g_c_per_h = ground_flux_g[@intFromEnum(ecosys.gas_transport.Species.carbon_dioxide)] + root_gas_withdrawal.carbon_dioxide_g_c_per_h - root_atmospheric_co2_g_c_per_h,
                                .net_carbon_exchange_g_c_per_h = canopy_net_fixation_g_c_per_h + ground_flux_g[@intFromEnum(ecosys.gas_transport.Species.carbon_dioxide)] + root_gas_withdrawal.carbon_dioxide_g_c_per_h - root_atmospheric_co2_g_c_per_h,
                                .methane_emission_g_c_per_h = ground_flux_g[@intFromEnum(ecosys.gas_transport.Species.methane)] + root_gas_withdrawal.methane_g_c_per_h - root_atmospheric_methane_g_c_per_h,
                                .oxygen_exchange_g_o2_per_h = ground_flux_g[@intFromEnum(ecosys.gas_transport.Species.oxygen)] + root_gas_withdrawal.oxygen_g_o_per_h - root_atmospheric_oxygen_g_o_per_h,
                                .local_surface_area_m2 = canopy_cell_area_m2[cell],
                                .carbon_dioxide_concentration_by_layer = carbon_dioxide_concentration,
                                .canopy_air_carbon_dioxide_umol_per_mol = current_atmospheric_co2_umol_per_mol,
                                .methane_concentration_by_layer = methane_concentration,
                                .oxygen_concentration_by_layer = oxygen_concentration,
                                .litter_oxygen_concentration = litter_oxygen,
                            }, values);
                            const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                            const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                            const file_name = try ecosys.output_record.buildCellFileName(allocator, grid_column, grid_row, current_year, editor_name);
                            defer allocator.free(file_name);
                            _ = try hourly_soil_carbon_bank.streams[cell].write(file_name, hourly_soil_carbon_catalog.variables, selection, resolved.soil_enabled, .{
                                .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                .grid_column = grid_column,
                                .grid_row = grid_row,
                                .values = values,
                            });
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 0;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlyPlantCarbonEditor;
                        if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                            const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                            const selection = output_selection_catalog.entries.items[selection_index].selection;
                            const canopy = if (detailed_canopy_state) |*value| value else return error.HourlyPlantCarbonOutputRequiresCanopyState;
                            const exchange = if (canopy_carbon_exchange_state) |*value| value else return error.HourlyPlantCarbonOutputRequiresCarbonExchangeState;
                            const surface = if (canopy_surface_exchange_state) |*value| value else return error.HourlyPlantCarbonOutputRequiresSurfaceExchangeState;
                            const surface_workspace = if (canopy_surface_input_workspace) |*value| value else return error.HourlyPlantCarbonOutputRequiresSurfaceWorkspace;
                            const layers = if (canopy_layer_distribution_state) |*value| value else return error.HourlyPlantCarbonOutputRequiresLayerDistribution;
                            for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                const plant = try canopy.plantIndex(cell, species);
                                const flux = try exchange.fluxForPlant(canopy, plant);
                                const output = try ecosys.plant_hourly_output.carbon(
                                    flux.net_co2_g_c_per_h,
                                    canopy.plant_gross_primary_productivity_g_c_per_h[plant],
                                    flux.signed_aboveground_respiration_g_c_per_h,
                                    canopy.plant_mobile_carbon_concentration_g_per_g[plant],
                                    surface_workspace.stomatal_resistance_h_per_m[plant],
                                    surface.boundary_layer_resistance_h_per_m[plant],
                                    try layers.plantLeafAreaM2(canopy, plant),
                                    canopy_cell_area_m2[cell],
                                );
                                const values = try hourly_plant_carbon_bank.row(plant);
                                values[0..7].* = output.values();
                                const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                const file_name = try ecosys.output_record.buildPlantFileName(allocator, grid_column, grid_row, species + 1, current_year, editor_name);
                                defer allocator.free(file_name);
                                _ = try hourly_plant_carbon_bank.streams[plant].write(file_name, hourly_plant_carbon_catalog.variables, selection, resolved.plant_enabled, .{
                                    .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                    .grid_column = grid_column,
                                    .grid_row = grid_row,
                                    .values = values,
                                });
                            };
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 1;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlySoilWaterEditor;
                        const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                        const selection = output_selection_catalog.entries.items[selection_index].selection;
                        for (0..state.cell_count) |cell| {
                            const values = try hourly_soil_water_bank.row(cell);
                            const liquid_fraction = values[6 .. 6 + config.soil_layers];
                            const surface_liquid_index = 6 + config.soil_layers;
                            const ice_fraction = values[surface_liquid_index + 1 ..][0..config.soil_layers];
                            var evapotranspiration_m3: f64 = 0;
                            var root_water_uptake_m3: f64 = 0;
                            for (0..config.plant_populations) |species| {
                                const plant = cell * config.plant_populations + species;
                                if (canopy_surface_exchange_state) |exchange| {
                                    evapotranspiration_m3 += @max(0, -exchange.transpiration_m3_per_h[plant]);
                                    evapotranspiration_m3 += @max(0, -exchange.intercepted_water_change_m3_per_h[plant]);
                                }
                                evapotranspiration_m3 += @max(0, -standing_dead_evaporation_m3_per_h[plant]);
                                if (plant_root_state) |roots| for (0..root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| for (0..config.soil_layers) |local_layer| {
                                    root_water_uptake_m3 += @max(0, -roots.water_uptake_m3_per_h[try roots.layerIndex(plant, domain, local_layer)]);
                                };
                            }
                            var external_water_outflow_m3: f64 = 0;
                            for (0..config.soil_layers) |local_layer| {
                                const layer = try state.layerIndex(cell, local_layer);
                                const layer_volume_m3 = soil_solver_property_state.layer_volume_m3[layer];
                                if (!std.math.isFinite(layer_volume_m3) or layer_volume_m3 <= 0) return error.InvalidOutputSoilLayerVolume;
                                liquid_fraction[local_layer] = state.liquid_water_m3[layer] / layer_volume_m3;
                                ice_fraction[local_layer] = state.ice_water_m3[layer] / layer_volume_m3;
                                external_water_outflow_m3 += @max(0, transport_hydrology_state.micropore_external_water_flux_m3_per_step[layer]);
                                external_water_outflow_m3 += @max(0, transport_hydrology_state.macropore_external_water_flux_m3_per_step[layer]);
                            }
                            const litter_capacity_m3 = surface_precipitation_state.litter_water_capacity_m3[cell];
                            const surface_liquid_fraction = if (litter_capacity_m3 > 0) surface_precipitation_state.litter_water_m3[cell] / litter_capacity_m3 else 0;
                            const surface_ice_fraction = if (litter_capacity_m3 > 0) surface_litter_ice_m3[cell] / litter_capacity_m3 else 0;
                            const top_layer = try state.layerIndex(cell, 0);
                            const bulk_density_Mg_per_m3 = soil_solver_property_state.bulk_density_megagrams_per_m3[top_layer];
                            if (!std.math.isFinite(bulk_density_Mg_per_m3) or bulk_density_Mg_per_m3 <= 0) return error.InvalidOutputSoilBulkDensity;
                            try ecosys.soil_water_output.calculateInto(.{
                                .evapotranspiration_m3 = evapotranspiration_m3,
                                .runoff_m3 = -surface_runoff_state.exported_water_m3[cell],
                                .sediment_discharge_water_m3 = surface_erosion_state.routing.sediment_export_Mg[cell] / bulk_density_Mg_per_m3,
                                .root_water_uptake_m3 = root_water_uptake_m3,
                                .external_water_outflow_m3 = external_water_outflow_m3,
                                .surface_snow_volume_m3 = surface_precipitation_state.solid_snow_water_equivalent_m3[cell],
                                .surface_ice_volume_m3 = surface_litter_ice_m3[cell],
                                .surface_liquid_water_m3 = surface_precipitation_state.litter_water_m3[cell],
                                // Dall'Amico soil/surface ice is already
                                // stored as water-equivalent volume.
                                .ice_density_megagrams_per_m3 = 1,
                                .local_surface_area_m2 = canopy_cell_area_m2[cell],
                                .total_grid_area_m2 = total_grid_area_m2,
                                .volumetric_liquid_water_fraction_by_layer = liquid_fraction,
                                .surface_volumetric_liquid_water_fraction = surface_liquid_fraction,
                                .volumetric_ice_fraction_by_layer = ice_fraction,
                                .surface_volumetric_ice_fraction = surface_ice_fraction,
                                .active_layer_boundary_depth_m = -soil_boundary_topology_state.active_layer_depth_m[cell],
                                .water_table_boundary_depth_m = -soil_boundary_topology_state.internal_water_table_depth_m[cell],
                                .soil_surface_reference_depth_m = 0,
                            }, values);
                            const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                            const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                            const file_name = try ecosys.output_record.buildCellFileName(allocator, grid_column, grid_row, current_year, editor_name);
                            defer allocator.free(file_name);
                            _ = try hourly_soil_water_bank.streams[cell].write(file_name, hourly_soil_water_catalog.variables, selection, resolved.soil_enabled, .{
                                .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                .grid_column = grid_column,
                                .grid_row = grid_row,
                                .values = values,
                            });
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 1;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlyPlantWaterEditor;
                        if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                            const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                            const selection = output_selection_catalog.entries.items[selection_index].selection;
                            const canopy = if (detailed_canopy_state) |*value| value else return error.HourlyPlantWaterOutputRequiresCanopyState;
                            const roots = if (plant_root_state) |*value| value else return error.HourlyPlantWaterOutputRequiresRootState;
                            const surface = if (canopy_surface_exchange_state) |*value| value else return error.HourlyPlantWaterOutputRequiresSurfaceExchangeState;
                            const surface_workspace = if (canopy_surface_input_workspace) |*value| value else return error.HourlyPlantWaterOutputRequiresSurfaceWorkspace;
                            for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                const plant = try canopy.plantIndex(cell, species);
                                const values = try hourly_plant_water_bank.row(plant);
                                const primary_root_first = try roots.layerIndex(plant, 0, 0);
                                const primary_root_potential = roots.total_water_potential_mpa[primary_root_first .. primary_root_first + config.soil_layers];
                                const outward_water_flux_m3_per_h = -(surface.transpiration_m3_per_h[plant] +
                                    surface.intercepted_water_change_m3_per_h[plant] +
                                    standing_dead_evaporation_m3_per_h[plant]);
                                try ecosys.plant_hourly_output.calculateWaterInto(
                                    canopy.plant_canopy_osmotic_potential_mpa[plant] + canopy.plant_canopy_turgor_potential_mpa[plant],
                                    canopy.plant_canopy_turgor_potential_mpa[plant],
                                    surface_workspace.stomatal_resistance_h_per_m[plant],
                                    surface.boundary_layer_resistance_h_per_m[plant],
                                    outward_water_flux_m3_per_h,
                                    phenology_root_oxygen_fraction[plant],
                                    primary_root_potential,
                                    canopy_cell_area_m2[cell],
                                    values,
                                );
                                const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                const file_name = try ecosys.output_record.buildPlantFileName(allocator, grid_column, grid_row, species + 1, current_year, editor_name);
                                defer allocator.free(file_name);
                                _ = try hourly_plant_water_bank.streams[plant].write(file_name, hourly_plant_water_catalog.variables, selection, resolved.plant_enabled, .{
                                    .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                    .grid_column = grid_column,
                                    .grid_row = grid_row,
                                    .values = values,
                                });
                            };
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 2;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlySoilNitrogenEditor;
                        const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                        const selection = output_selection_catalog.entries.items[selection_index].selection;
                        for (0..state.cell_count) |cell| {
                            const root_gas_withdrawal = if (plant_root_state) |*roots|
                                try ecosys.plant_root_disturbance.rootGasWithdrawalForCell(roots, cell, config.plant_populations)
                            else
                                ecosys.plant_root_disturbance.CellRootGasWithdrawal{};
                            const values = try hourly_soil_nitrogen_bank.row(cell);
                            const nitrous_oxide_concentration = values[5 .. 5 + config.soil_layers];
                            const litter_nitrous_oxide_index = 5 + config.soil_layers;
                            const ammonia_concentration = values[litter_nitrous_oxide_index + 1 ..][0..config.soil_layers];
                            for (0..config.soil_layers) |local_layer| {
                                const layer = try state.layerIndex(cell, local_layer);
                                const water_m3 = state.liquid_water_m3[layer];
                                const n2o_component = layer * ecosys.gas_transport.species_count + @intFromEnum(ecosys.gas_transport.Species.nitrous_oxide);
                                nitrous_oxide_concentration[local_layer] = try ecosys.soil_biogeochemistry_output.dissolvedGasConcentration(gas_transport_state.dissolved_mass_g[n2o_component], water_m3);
                            }
                            const first_soil_layer = try state.layerIndex(cell, 0);
                            try ecosys.soil_hourly_output_binding.writeMineralAmmoniaNitrogenProfile(
                                initial_chemistry_state.aqueous,
                                first_soil_layer,
                                config.soil_layers,
                                runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                                ammonia_concentration,
                            );
                            const nitrous_oxide_exchange_g_n_per_h = try ecosys.soil_hourly_output_binding.gasBoundaryExchangeG(
                                soil_gas_transport_state.atmospheric_flux_g_per_h,
                                surface_litter_gas_transport_state.atmospheric_flux_g_per_h,
                                first_soil_layer,
                                config.soil_layers,
                                cell,
                                .nitrous_oxide,
                            );
                            const dinitrogen_exchange_g_n_per_h = try ecosys.soil_hourly_output_binding.gasBoundaryExchangeG(
                                soil_gas_transport_state.atmospheric_flux_g_per_h,
                                surface_litter_gas_transport_state.atmospheric_flux_g_per_h,
                                first_soil_layer,
                                config.soil_layers,
                                cell,
                                .nitrogen,
                            );
                            const ammonia_exchange_g_n_per_h = try ecosys.soil_hourly_output_binding.gasBoundaryExchangeG(
                                soil_gas_transport_state.atmospheric_flux_g_per_h,
                                surface_litter_gas_transport_state.atmospheric_flux_g_per_h,
                                first_soil_layer,
                                config.soil_layers,
                                cell,
                                .ammonia,
                            );
                            const litter_base = cell * ecosys.gas_transport.species_count;
                            const litter_nitrous_oxide = try ecosys.soil_biogeochemistry_output.dissolvedGasConcentration(
                                litter_gas_transport_state.dissolved_mass_g[litter_base + @intFromEnum(ecosys.gas_transport.Species.nitrous_oxide)],
                                surface_litter_gas_transport_state.water_volume_m3[cell],
                            );
                            const litter_ammonia_mol_n_per_m3 = surface_litter_chemistry_state.cells[cell].ammonia_mol_per_m3;
                            if (!std.math.isFinite(litter_ammonia_mol_n_per_m3) or litter_ammonia_mol_n_per_m3 < 0) return error.InvalidOutputLitterAmmoniaConcentration;
                            var dissolved_inorganic_nitrogen_drainage_g_n_per_h: f64 = 0;
                            for (0..config.soil_layers) |local_layer| {
                                dissolved_inorganic_nitrogen_drainage_g_n_per_h += mineral_nitrogen_transport_state.boundary_export_g_n_per_step[try state.layerIndex(cell, local_layer)];
                            }
                            const root_atmospheric_nitrous_oxide_g_n_per_h = try ecosys.root_uptake_ledger.atmosphereToRootForCellGPerH(&root_uptake_ledger_state, &state, cell, 3);
                            const root_atmospheric_ammonia_g_n_per_h = try ecosys.root_uptake_ledger.atmosphereToRootForCellGPerH(&root_uptake_ledger_state, &state, cell, 4);
                            try ecosys.soil_biogeochemistry_output.calculateNitrogenInto(.{
                                .nitrous_oxide_emission_g_n_per_h = nitrous_oxide_exchange_g_n_per_h + root_gas_withdrawal.nitrous_oxide_g_n_per_h - root_atmospheric_nitrous_oxide_g_n_per_h,
                                .dinitrogen_emission_g_n_per_h = dinitrogen_exchange_g_n_per_h,
                                .ammonia_emission_g_n_per_h = ammonia_exchange_g_n_per_h + root_gas_withdrawal.ammonia_g_n_per_h - root_atmospheric_ammonia_g_n_per_h,
                                .dissolved_inorganic_nitrogen_runoff_g_n_per_h = surface_inorganic_nitrogen_export_g_n_per_h[cell],
                                .dissolved_inorganic_nitrogen_drainage_g_n_per_h = dissolved_inorganic_nitrogen_drainage_g_n_per_h,
                                .local_surface_area_m2 = canopy_cell_area_m2[cell],
                                .total_grid_area_m2 = total_grid_area_m2,
                                .nitrous_oxide_concentration_by_layer = nitrous_oxide_concentration,
                                .litter_nitrous_oxide_concentration = litter_nitrous_oxide,
                                .ammonia_concentration_by_layer = ammonia_concentration,
                                .litter_ammonia_concentration = litter_ammonia_mol_n_per_m3 * 14,
                            }, values);
                            const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                            const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                            const file_name = try ecosys.output_record.buildCellFileName(allocator, grid_column, grid_row, current_year, editor_name);
                            defer allocator.free(file_name);
                            _ = try hourly_soil_nitrogen_bank.streams[cell].write(file_name, hourly_soil_nitrogen_catalog.variables, selection, resolved.soil_enabled, .{
                                .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                .grid_column = grid_column,
                                .grid_row = grid_row,
                                .values = values,
                            });
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 2;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlyPlantNitrogenEditor;
                        if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                            const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                            const selection = output_selection_catalog.entries.items[selection_index].selection;
                            const canopy = if (detailed_canopy_state) |*value| value else return error.HourlyPlantNitrogenOutputRequiresCanopyState;
                            const roots = if (plant_root_state) |*value| value else return error.HourlyPlantNitrogenOutputRequiresRootState;
                            for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                const plant = try canopy.plantIndex(cell, species);
                                const root_first = try roots.layerIndex(plant, 0, 0);
                                const mycorrhiza_first = try roots.layerIndex(plant, 1, 0);
                                const has_mycorrhiza = root_biological_domain_count_by_plant[plant] == 2;
                                const first_soil_layer = try state.layerIndex(cell, 0);
                                const values = try hourly_plant_nitrogen_bank.row(plant);
                                try ecosys.plant_hourly_output.calculateNitrogenInto(
                                    roots.ammonium_uptake_g_n_per_h[plant],
                                    roots.nitrate_uptake_g_n_per_h[plant],
                                    roots.fixation_uptake_g_n_per_h[plant],
                                    canopy.plant_mobile_nitrogen_concentration_g_per_g[plant],
                                    canopy_ammonia_publication_state
                                        .current_exchange_g_n_per_h_by_plant[plant],
                                    canopy_cell_area_m2[cell],
                                    .{
                                        .ammonium_root_non_band_g_n_per_h = roots.ammonium_uptake_nonband_g_n_per_h[root_first .. root_first + config.soil_layers],
                                        .ammonium_mycorrhiza_non_band_g_n_per_h = if (has_mycorrhiza) roots.ammonium_uptake_nonband_g_n_per_h[mycorrhiza_first .. mycorrhiza_first + config.soil_layers] else zero_root_layer_flux_g_per_h,
                                        .ammonium_root_band_g_n_per_h = roots.ammonium_uptake_band_g_n_per_h[root_first .. root_first + config.soil_layers],
                                        .ammonium_mycorrhiza_band_g_n_per_h = if (has_mycorrhiza) roots.ammonium_uptake_band_g_n_per_h[mycorrhiza_first .. mycorrhiza_first + config.soil_layers] else zero_root_layer_flux_g_per_h,
                                        .nitrate_root_non_band_g_n_per_h = roots.nitrate_uptake_nonband_g_n_per_h[root_first .. root_first + config.soil_layers],
                                        .nitrate_mycorrhiza_non_band_g_n_per_h = if (has_mycorrhiza) roots.nitrate_uptake_nonband_g_n_per_h[mycorrhiza_first .. mycorrhiza_first + config.soil_layers] else zero_root_layer_flux_g_per_h,
                                        .nitrate_root_band_g_n_per_h = roots.nitrate_uptake_band_g_n_per_h[root_first .. root_first + config.soil_layers],
                                        .nitrate_mycorrhiza_band_g_n_per_h = if (has_mycorrhiza) roots.nitrate_uptake_band_g_n_per_h[mycorrhiza_first .. mycorrhiza_first + config.soil_layers] else zero_root_layer_flux_g_per_h,
                                        .layer_area_m2 = soil_hourly_workspace.horizontal_face_area_m2[first_soil_layer .. first_soil_layer + config.soil_layers],
                                    },
                                    values,
                                );
                                const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                const file_name = try ecosys.output_record.buildPlantFileName(allocator, grid_column, grid_row, species + 1, current_year, editor_name);
                                defer allocator.free(file_name);
                                _ = try hourly_plant_nitrogen_bank.streams[plant].write(file_name, hourly_plant_nitrogen_catalog.variables, selection, resolved.plant_enabled, .{
                                    .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                    .grid_column = grid_column,
                                    .grid_row = grid_row,
                                    .values = values,
                                });
                            };
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 4;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlySoilHeatEditor;
                        const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                        const selection = output_selection_catalog.entries.items[selection_index].selection;
                        for (0..state.cell_count) |cell| {
                            const values = try hourly_soil_heat_bank.row(cell);
                            const soil_temperature_c = values[13 .. 13 + config.soil_layers];
                            for (soil_temperature_c, 0..) |*temperature_c, local_layer| temperature_c.* = state.soil_temperature_k[try state.layerIndex(cell, local_layer)] - 273.15;
                            const litter_air_volume_m3 = litter_gas_transport_state.air_volume_m3[cell];
                            const litter_vapor_pressure_kpa = if (litter_air_volume_m3 > 0)
                                litter_gas_transport_state.water_vapor_mol[cell] * 8.31446261815324e-3 * litter_gas_transport_state.temperature_k[cell] / litter_air_volume_m3
                            else
                                0;
                            try ecosys.soil_heat_output.calculateInto(.{
                                .incoming_shortwave_radiation_mj_per_m2_h = atmospheric_state.shortwave_radiation_mj_per_m2[cell],
                                .air_temperature_c = atmospheric_state.air_temperature_k[cell] - 273.15,
                                .atmospheric_vapor_pressure_kpa = atmospheric_state.vapor_pressure_kpa[cell],
                                .wind_travel_m_per_h = atmospheric_state.wind_speed_m_per_h[cell],
                                .rainfall_m3 = atmospheric_state.rainfall_m[cell] * canopy_cell_area_m2[cell],
                                .irrigation_m3 = irrigation_water_depth_m[cell] * canopy_cell_area_m2[cell],
                                .local_surface_area_m2 = canopy_cell_area_m2[cell],
                                .ground_surface_net_radiation_mj = ecosystem_energy_ledger_state.ground_surface_net_radiation_mj[cell],
                                .ground_surface_latent_heat_mj = ecosystem_energy_ledger_state.ground_surface_latent_heat_mj[cell],
                                .ground_surface_sensible_heat_mj = ecosystem_energy_ledger_state.ground_surface_sensible_heat_mj[cell],
                                .ground_surface_storage_heat_mj = ecosystem_energy_ledger_state.ground_surface_storage_heat_mj[cell],
                                .ecosystem_net_radiation_mj = ecosystem_energy_ledger_state.ecosystem_net_radiation_mj[cell],
                                .ecosystem_latent_heat_mj = ecosystem_energy_ledger_state.ecosystem_latent_heat_mj[cell],
                                .ecosystem_sensible_heat_mj = ecosystem_energy_ledger_state.ecosystem_sensible_heat_mj[cell],
                                .ecosystem_storage_heat_mj = ecosystem_energy_ledger_state.ecosystem_storage_heat_mj[cell],
                                .soil_temperature_c_by_layer = soil_temperature_c,
                                .surface_soil_temperature_c = state.soil_temperature_k[try state.layerIndex(cell, 0)] - 273.15,
                                .surface_water_temperature_c = state.surface_temperature_k[cell] - 273.15,
                                .litter_temperature_c = litter_gas_transport_state.temperature_k[cell] - 273.15,
                                .litter_water_vapor_partial_pressure_kpa = litter_vapor_pressure_kpa,
                                .litter_absolute_temperature_k = litter_gas_transport_state.temperature_k[cell],
                            }, values);
                            const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                            const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                            const file_name = try ecosys.output_record.buildCellFileName(allocator, grid_column, grid_row, current_year, editor_name);
                            defer allocator.free(file_name);
                            _ = try hourly_soil_heat_bank.streams[cell].write(file_name, hourly_soil_heat_catalog.variables, selection, resolved.soil_enabled, .{
                                .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                .grid_column = grid_column,
                                .grid_row = grid_row,
                                .values = values,
                            });
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 4;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlyPlantHeatEditor;
                        if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                            const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                            const selection = output_selection_catalog.entries.items[selection_index].selection;
                            const canopy = if (detailed_canopy_state) |*value| value else return error.HourlyPlantHeatOutputRequiresCanopyState;
                            const living = if (canopy_surface_exchange_state) |*value| value else return error.HourlyPlantHeatOutputRequiresSurfaceExchangeState;
                            const dead = if (standing_dead_surface_exchange_state) |*value| value else return error.HourlyPlantHeatOutputRequiresStandingDeadExchangeState;
                            for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                const plant = try canopy.plantIndex(cell, species);
                                const temperature_function = canopy.plant_uptake_growth_temperature_response[plant];
                                const output = try ecosys.plant_hourly_output.heat(
                                    canopy_net_radiation_mj[plant] + dead.net_radiation_mj_per_h[plant],
                                    living.latent_heat_flux_mj_per_h[plant] + dead.latent_heat_flux_mj_per_h[plant],
                                    living.sensible_heat_flux_mj_per_h[plant] + dead.sensible_heat_flux_mj_per_h[plant],
                                    -canopy_storage_heat_mj[plant] + living.vapor_sensible_heat_flux_mj_per_h[plant] - dead.storage_heat_flux_mj_per_h[plant] + dead.vapor_sensible_heat_flux_mj_per_h[plant],
                                    plant_state.canopy_temperature_k[plant] - 273.15,
                                    temperature_function,
                                    canopy.plant_standing_dead_surface_temperature_k[plant] - 273.15,
                                    canopy_cell_area_m2[cell],
                                );
                                const values = try hourly_plant_heat_bank.row(plant);
                                values[0..7].* = output.values();
                                const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                const file_name = try ecosys.output_record.buildPlantFileName(allocator, grid_column, grid_row, species + 1, current_year, editor_name);
                                defer allocator.free(file_name);
                                _ = try hourly_plant_heat_bank.streams[plant].write(file_name, hourly_plant_heat_catalog.variables, selection, resolved.plant_enabled, .{
                                    .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                    .grid_column = grid_column,
                                    .grid_row = grid_row,
                                    .values = values,
                                });
                            };
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 3;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlySoilPhosphorusEditor;
                        const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                        const selection = output_selection_catalog.entries.items[selection_index].selection;
                        for (0..state.cell_count) |cell| {
                            const values = try hourly_soil_phosphorus_bank.row(cell);
                            var dissolved_inorganic_phosphorus_drainage_g_p_per_h: f64 = 0;
                            for (0..config.soil_layers) |local_layer| {
                                const layer = try state.layerIndex(cell, local_layer);
                                const base = layer * aqueous_species_count;
                                inline for (@typeInfo(ecosys.solute_transport_species.AqueousSpecies).@"enum".fields) |field| {
                                    const species: ecosys.solute_transport_species.AqueousSpecies = @enumFromInt(field.value);
                                    if (ecosys.solute_transport_species.diffusivityClass(species) == .phosphate) {
                                        dissolved_inorganic_phosphorus_drainage_g_p_per_h += @max(0, -soil_solute_boundary_net_flux_mol[base + field.value]) * runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol;
                                    }
                                }
                            }
                            const phosphorus = try ecosys.soil_biogeochemistry_output.calculatePhosphorus(
                                surface_inorganic_phosphorus_export_g_p_per_h[cell],
                                dissolved_inorganic_phosphorus_drainage_g_p_per_h,
                                total_grid_area_m2,
                            );
                            values[0..2].* = ecosys.soil_biogeochemistry_output.phosphorusValues(phosphorus);
                            const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                            const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                            const file_name = try ecosys.output_record.buildCellFileName(allocator, grid_column, grid_row, current_year, editor_name);
                            defer allocator.free(file_name);
                            _ = try hourly_soil_phosphorus_bank.streams[cell].write(file_name, hourly_soil_phosphorus_catalog.variables, selection, resolved.soil_enabled, .{
                                .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                .grid_column = grid_column,
                                .grid_row = grid_row,
                                .values = values,
                            });
                        }
                    }
                }
                if (@mod(completed_hour_count, active_options.hourly_output_interval_hours) == 0) {
                    const editor_index: usize = 3;
                    const editor_name = scene.output_editors[editor_index];
                    if (!ecosys.delimited_input.isNo(editor_name)) {
                        const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedHourlyPlantPhosphorusEditor;
                        if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                            const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                            const selection = output_selection_catalog.entries.items[selection_index].selection;
                            const canopy = if (detailed_canopy_state) |*value| value else return error.HourlyPlantPhosphorusOutputRequiresCanopyState;
                            const roots = if (plant_root_state) |*value| value else return error.HourlyPlantPhosphorusOutputRequiresRootState;
                            for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                const plant = try canopy.plantIndex(cell, species);
                                const root_first = try roots.layerIndex(plant, 0, 0);
                                const mycorrhiza_first = try roots.layerIndex(plant, 1, 0);
                                const has_mycorrhiza = root_biological_domain_count_by_plant[plant] == 2;
                                const first_soil_layer = try state.layerIndex(cell, 0);
                                const root_non_band = roots.phosphate_h2_uptake_nonband_g_p_per_h[root_first .. root_first + config.soil_layers];
                                const mycorrhiza_non_band = if (has_mycorrhiza) roots.phosphate_h2_uptake_nonband_g_p_per_h[mycorrhiza_first .. mycorrhiza_first + config.soil_layers] else zero_root_layer_flux_g_per_h;
                                const root_band = roots.phosphate_h2_uptake_band_g_p_per_h[root_first .. root_first + config.soil_layers];
                                const mycorrhiza_band = if (has_mycorrhiza) roots.phosphate_h2_uptake_band_g_p_per_h[mycorrhiza_first .. mycorrhiza_first + config.soil_layers] else zero_root_layer_flux_g_per_h;
                                var h2_phosphate_uptake_g_p_per_h: f64 = 0;
                                for (0..config.soil_layers) |layer| h2_phosphate_uptake_g_p_per_h += root_non_band[layer] + mycorrhiza_non_band[layer] + root_band[layer] + mycorrhiza_band[layer];
                                const values = try hourly_plant_phosphorus_bank.row(plant);
                                try ecosys.plant_hourly_output.calculatePhosphorusInto(
                                    h2_phosphate_uptake_g_p_per_h,
                                    canopy.plant_mobile_phosphorus_concentration_g_per_g[plant],
                                    root_non_band,
                                    mycorrhiza_non_band,
                                    root_band,
                                    mycorrhiza_band,
                                    soil_hourly_workspace.horizontal_face_area_m2[first_soil_layer .. first_soil_layer + config.soil_layers],
                                    canopy_cell_area_m2[cell],
                                    values,
                                );
                                const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                const file_name = try ecosys.output_record.buildPlantFileName(allocator, grid_column, grid_row, species + 1, current_year, editor_name);
                                defer allocator.free(file_name);
                                _ = try hourly_plant_phosphorus_bank.streams[plant].write(file_name, hourly_plant_phosphorus_catalog.variables, selection, resolved.plant_enabled, .{
                                    .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = if (timestamp.hour == 24) 23 else timestamp.hour },
                                    .grid_column = grid_column,
                                    .grid_row = grid_row,
                                    .values = values,
                                });
                            };
                        }
                    }
                }
                if (active_options.visualizationIncludesYear(current_year)) {
                    for (0..state.cell_count) |cell| {
                        const first_plant = cell * config.plant_populations;
                        const plant_end =
                            first_plant + config.plant_populations;
                        var ecosystem_gross_primary_productivity_g: f64 = 0;
                        var ecosystem_signed_autotrophic_respiration_g: f64 = 0;
                        var harvested_carbon_g: f64 = 0;
                        for (0..config.plant_populations) |species| {
                            const plant = first_plant + species;
                            ecosystem_gross_primary_productivity_g +=
                                plant_daily_flux_ledger.gross_primary_productivity_g[plant];
                            ecosystem_signed_autotrophic_respiration_g +=
                                plant_daily_flux_ledger.signed_total_respiration_carbon_g[plant];
                            harvested_carbon_g +=
                                plant_daily_flux_ledger.harvested_carbon_g[plant];
                            visualization_species_respired_carbon_g[species] =
                                -plant_daily_flux_ledger.signed_total_respiration_carbon_g[plant];
                        }
                        const ecosystem_net_primary_productivity_g =
                            ecosystem_gross_primary_productivity_g +
                            ecosystem_signed_autotrophic_respiration_g;
                        const ecosystem_carbon = try ecosys.soil_daily_ecosystem_carbon.calculate(.{
                            .net_primary_productivity_g_c = ecosystem_net_primary_productivity_g,
                            .signed_heterotrophic_respiration_g_c = daily_heterotrophic_respiration.daily_signed_g_c[cell],
                            .dissolved_organic_carbon_runoff_g_c = daily_carbon_export.dissolved_organic_carbon_runoff_g[cell],
                            .dissolved_inorganic_carbon_runoff_g_c = daily_carbon_export.dissolved_inorganic_carbon_runoff_g[cell],
                            .dissolved_organic_carbon_drainage_g_c = daily_carbon_export.dissolved_organic_carbon_drainage_g[cell],
                            .dissolved_inorganic_carbon_drainage_g_c = daily_carbon_export.dissolved_inorganic_carbon_drainage_g[cell],
                            .harvested_carbon_g_c = harvested_carbon_g,
                            .organic_fertilizer_carbon_input_g_c = daily_biome_organic_carbon_input_g_c[cell],
                        });
                        const first_layer = try state.layerIndex(cell, 0);
                        const layer_end = first_layer + config.soil_layers;
                        const grid_column =
                            runscript.domain.west_column +
                            cell % config.grid_columns;
                        const grid_row =
                            runscript.domain.north_row +
                            cell / config.grid_columns;
                        _ = try visualization_streams.calculateAndWriteHourly(
                            cell,
                            grid_column,
                            grid_row,
                            .{
                                .year = current_year,
                                .day_of_year = current_day_of_year,
                                .month = management_date.month,
                                .day = management_date.day,
                                .hour = if (timestamp.hour == 24)
                                    23
                                else
                                    timestamp.hour,
                            },
                            .{
                                .cell_area_m2 = canopy_cell_area_m2[cell],
                                .species_cumulative_gross_carbon_g = plant_daily_flux_ledger.gross_primary_productivity_g[first_plant..plant_end],
                                .species_cumulative_respired_carbon_g = visualization_species_respired_carbon_g,
                                .ecosystem_gross_primary_productivity_g = ecosystem_gross_primary_productivity_g,
                                .ecosystem_net_primary_productivity_g = ecosystem_net_primary_productivity_g,
                                .ecosystem_autotrophic_respiration_g = -ecosystem_signed_autotrophic_respiration_g,
                                .ecosystem_heterotrophic_respiration_g = -daily_heterotrophic_respiration.daily_signed_g_c[cell],
                                .net_biome_productivity_g = ecosystem_carbon.net_biome_productivity_g_c,
                                .carbon_dioxide_flux_g = (try daily_soil_gas_flux.get(cell, .carbon_dioxide)) +
                                    daily_canopy_gas_exchange.net_carbon_dioxide_uptake_g_c[cell],
                                .methane_flux_g = (try daily_soil_gas_flux.get(cell, .methane)) +
                                    daily_canopy_gas_exchange.net_methane_uptake_g_c[cell],
                                .net_radiation_mj_h = ecosystem_energy_ledger_state.ecosystem_net_radiation_mj[cell],
                                .latent_heat_flux_mj_h = ecosystem_energy_ledger_state.ecosystem_latent_heat_mj[cell],
                                .sensible_heat_flux_mj_h = ecosystem_energy_ledger_state.ecosystem_sensible_heat_mj[cell],
                                .liquid_water_m3_by_layer = state.matrix_liquid_water_m3[first_layer..layer_end],
                                .air_volume_m3_by_layer = state.matrix_air_volume_m3[first_layer..layer_end],
                                .macropore_water_m3_by_layer = state.macropore_liquid_water_m3[first_layer..layer_end],
                                .total_volume_m3_by_layer = soil_solver_property_state.layer_volume_m3[first_layer..layer_end],
                                .soil_temperature_k_by_layer = state.soil_temperature_k[first_layer..layer_end],
                            },
                        );
                    }
                }
                if (timestamp.hour == 24) {
                    try daily_ecosystem_carbon.finalizeDay(
                        &state,
                        config.plant_populations,
                        plant_daily_flux_ledger.gross_primary_productivity_g,
                        plant_daily_flux_ledger.signed_total_respiration_carbon_g,
                        plant_daily_flux_ledger.harvested_carbon_g,
                        daily_heterotrophic_respiration.daily_signed_g_c,
                        daily_biome_organic_carbon_input_g_c,
                        &daily_carbon_export,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedWater(
                        daily_water_ledger.rainfall_m3,
                        daily_water_ledger.runoff_m3,
                        daily_water_ledger.evaporation_m3,
                        daily_water_ledger.water_outflow_m3,
                        daily_water_ledger.lateral_water_outflow_m3,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedFertilizer(
                        soil_fertilizer_inventory.daily_nitrogen_input_g_n,
                        daily_organic_fertilizer_nitrogen_input_g_n,
                        mineral_fertilizer_inventory.daily_phosphorus_input_g_p,
                        daily_organic_fertilizer_phosphorus_input_g_p,
                        daily_organic_fertilizer_carbon_input_g_c,
                        daily_biome_organic_carbon_input_g_c,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedExports(
                        .{
                            daily_carbon_export.dissolved_organic_carbon_runoff_g,
                            daily_carbon_export.dissolved_inorganic_carbon_runoff_g,
                            daily_carbon_export.dissolved_organic_carbon_drainage_g,
                            daily_carbon_export.dissolved_inorganic_carbon_drainage_g,
                        },
                        .{
                            daily_nitrogen_export.dissolved_organic_nitrogen_runoff_g_n,
                            daily_nitrogen_export.dissolved_inorganic_nitrogen_runoff_g_n,
                            daily_nitrogen_export.dissolved_organic_nitrogen_drainage_g_n,
                            daily_nitrogen_export.dissolved_inorganic_nitrogen_drainage_g_n,
                        },
                        .{
                            daily_phosphorus_export.dissolved_organic_phosphorus_runoff_g_p,
                            daily_phosphorus_export.dissolved_inorganic_phosphorus_runoff_g_p,
                            daily_phosphorus_export.dissolved_organic_phosphorus_drainage_g_p,
                            daily_phosphorus_export.dissolved_inorganic_phosphorus_drainage_g_p,
                        },
                        daily_heat_ledger.ionic_outflow_mol,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedPlantLitter(
                        plant_daily_flux_ledger.carbon_sink_g,
                        plant_daily_flux_ledger.nitrogen_sink_g,
                        plant_daily_flux_ledger.phosphorus_sink_g,
                    );
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedAtmosphericGas(
                        &daily_soil_gas_flux,
                        &daily_canopy_gas_exchange,
                    );
                    const atmospheric_ion_molar_masses: ecosys.snow_surface_discharge.IonMolarMassesGPerMol =
                        if (runscript.chemistry_primary_initialization) |parameters|
                            .{
                                .aluminum = parameters.molar_mass_g_per_mol.aluminum,
                                .iron = parameters.molar_mass_g_per_mol.iron,
                                .calcium = parameters.molar_mass_g_per_mol.calcium,
                                .magnesium = parameters.molar_mass_g_per_mol.magnesium,
                                .sodium = parameters.molar_mass_g_per_mol.sodium,
                                .potassium = parameters.molar_mass_g_per_mol.potassium,
                                .sulfur = parameters.molar_mass_g_per_mol.sulfur,
                                .chloride = parameters.molar_mass_g_per_mol.chloride,
                            }
                        else
                            .{
                                .aluminum = 27,
                                .iron = 55.8,
                                .calcium = 40,
                                .magnesium = 24.3,
                                .sodium = 23,
                                .potassium = 39.1,
                                .sulfur = 32,
                                .chloride = 35.5,
                            };
                    try landscape_mass_balance_state.boundary_ledger.accumulateAcceptedAtmosphericSolutes(
                        &atmospheric_solute_input_ledger_state,
                        atmospheric_ion_molar_masses,
                    );
                    const completed_scene_day = (scene_weather_hours + 1) / 24;
                    if (try ecosys.mass_balance_audit.shouldAudit(
                        completed_scene_day,
                        active_options.daily_output_interval_days,
                    )) {
                        const totals = try reconstructLandscapeMassBalance(
                            hourly_science_context,
                        );
                        const monitor = landscape_mass_balance_state.monitor orelse
                            return error.MassBalanceMonitorNotInitialized;
                        _ = try monitor.check(
                            completed_scene_day,
                            current_year,
                            totals,
                        );
                    }
                    if (@mod(completed_scene_day, active_options.daily_output_interval_days) == 0) {
                        const carbon_editor_index: usize = 5;
                        const carbon_editor_name = scene.output_editors[carbon_editor_index];
                        if (!ecosys.delimited_input.isNo(carbon_editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + carbon_editor_index] orelse return error.UnresolvedDailyPlantCarbonEditor;
                            if (std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(carbon_editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                for (0..state.cell_count) |cell| {
                                    const values = try daily_soil_carbon_bank.row(cell);
                                    @memset(values, std.math.nan(f64));
                                    const first_layer = try state.layerIndex(cell, 0);
                                    const active_layers = state.active_soil_layer_count[cell];
                                    const layer_carbon = daily_soil_organic_carbon_by_layer_g_c[first_layer .. first_layer + config.soil_layers];
                                    @memset(layer_carbon, 0);
                                    const pools = try ecosys.soil_daily_carbon_pools.calculate(
                                        &soil_organic_state,
                                        &surface_organic_state,
                                        first_layer,
                                        active_layers,
                                        cell,
                                        layer_carbon[0..active_layers],
                                    );
                                    const plant_flux = try ecosys.soil_daily_plant_carbon_flux.calculate(&plant_daily_flux_ledger, cell, config.plant_populations);
                                    var gas_inventory_g: [ecosys.gas_transport.species_count]f64 = @splat(0);
                                    for (0..active_layers) |local_layer| {
                                        const layer = first_layer + local_layer;
                                        for (0..ecosys.gas_transport.species_count) |species| {
                                            const component = layer * ecosys.gas_transport.species_count + species;
                                            gas_inventory_g[species] +=
                                                gas_transport_state.gaseous_mass_g[component] +
                                                gas_transport_state.dissolved_mass_g[component] +
                                                gas_transport_state.macropore_dissolved_mass_g[component] +
                                                gas_transport_state.band_dissolved_mass_g[component];
                                        }
                                    }
                                    for (0..ecosys.gas_transport.species_count) |species| {
                                        const component = cell * ecosys.gas_transport.species_count + species;
                                        gas_inventory_g[species] +=
                                            litter_gas_transport_state.gaseous_mass_g[component] +
                                            litter_gas_transport_state.dissolved_mass_g[component] +
                                            litter_gas_transport_state.macropore_dissolved_mass_g[component] +
                                            litter_gas_transport_state.band_dissolved_mass_g[component];
                                    }
                                    const area_inverse = 1.0 / canopy_cell_area_m2[cell];
                                    values[0] = pools.residue_carbon_g_c * area_inverse;
                                    values[1] = pools.humus_carbon_g_c * area_inverse;
                                    values[2] = daily_organic_fertilizer_carbon_input_g_c[cell] * area_inverse;
                                    values[3] = plant_flux.plant_litterfall_g_c * area_inverse;
                                    values[4] = (try daily_soil_gas_flux.getSoilLitterBoundary(cell, .carbon_dioxide)) * area_inverse;
                                    values[5] = (try daily_soil_gas_flux.getSoilLitterBoundary(cell, .oxygen)) * area_inverse;
                                    values[6] = daily_heterotrophic_respiration.daily_carbon_dioxide_production_g_c[cell] * area_inverse;
                                    values[7] = pools.microbial_carbon_g_c * area_inverse;
                                    values[8] = pools.surface_noncharcoal_organic_carbon_g_c * area_inverse;
                                    values[9] = (try daily_soil_gas_flux.getSoilLitterBoundary(cell, .methane)) * area_inverse;
                                    values[10] = daily_carbon_export.dissolved_organic_carbon_runoff_g[cell] / total_grid_area_m2;
                                    values[11] = daily_carbon_export.dissolved_organic_carbon_drainage_g[cell] / total_grid_area_m2;
                                    values[12] = daily_carbon_export.dissolved_inorganic_carbon_runoff_g[cell] / total_grid_area_m2;
                                    values[13] = daily_carbon_export.dissolved_inorganic_carbon_drainage_g[cell] / total_grid_area_m2;
                                    values[14] = current_atmospheric_co2_umol_per_mol;
                                    values[15] = daily_ecosystem_carbon.net_biome_productivity_g_c[cell] * area_inverse;
                                    values[16] = daily_soil_fire_carbon_dioxide_emission_g_c[cell] * area_inverse;
                                    for (0..config.soil_layers) |layer| values[17 + layer] = layer_carbon[layer] * area_inverse;
                                    const after_layers = 17 + config.soil_layers;
                                    values[after_layers] = daily_soil_fire_charcoal_production_g_c[cell] * area_inverse;
                                    values[after_layers + 1] =
                                        ((try daily_soil_gas_flux.get(cell, .carbon_dioxide)) +
                                            daily_canopy_gas_exchange.net_carbon_dioxide_uptake_g_c[cell]) * area_inverse;
                                    values[after_layers + 2] =
                                        ((try daily_soil_gas_flux.get(cell, .methane)) +
                                            daily_canopy_gas_exchange.net_methane_uptake_g_c[cell]) * area_inverse;
                                    values[after_layers + 3] =
                                        ((try daily_soil_gas_flux.get(cell, .oxygen)) +
                                            daily_canopy_gas_exchange.net_oxygen_uptake_g_o[cell]) * area_inverse;
                                    // OUTSD choices 36..40 are source-reserved zeros.
                                    @memset(values[after_layers + 4 .. after_layers + 9], 0);
                                    const tail = after_layers + 9;
                                    // OUTSD choice 41 is UH2GG: accepted daily
                                    // ecosystem-atmosphere H2 exchange, not
                                    // the end-of-day H2 storage inventory.
                                    values[tail] = try ecosys.soil_daily_output.dailyHydrogenFluxGPerM2(
                                        &daily_soil_gas_flux,
                                        cell,
                                        canopy_cell_area_m2[cell],
                                    );
                                    values[tail + 1] = plant_flux.harvested_carbon_g_c * area_inverse;
                                    // With no active plant branches, OUTSD ARLFC and WTSTGT
                                    // are both authoritative zero-valued ecosystem pools.
                                    values[tail + 2] = 0;
                                    values[tail + 9] = 0;
                                    if (detailed_canopy_state) |canopy| {
                                        var total_leaf_area_m2: f64 = 0;
                                        for (0..config.plant_populations) |species| {
                                            const plant = cell * config.plant_populations + species;
                                            const branches = try canopy.branchRange(plant);
                                            for (canopy.branch_leaf_area_m2[branches.first..branches.end]) |leaf_area_m2| total_leaf_area_m2 += leaf_area_m2;
                                        }
                                        values[tail + 2] = total_leaf_area_m2 * area_inverse;
                                        values[tail + 9] =
                                            cell_litter_standing_dead_publication_state
                                                .standing_dead_carbon_g_c_by_cell[cell] *
                                            area_inverse;
                                    }
                                    values[tail + 3] = plant_flux.gross_primary_productivity_g_c * area_inverse;
                                    values[tail + 4] = plant_flux.autotrophic_respiration_g_c * area_inverse;
                                    values[tail + 5] = plant_flux.net_primary_productivity_g_c * area_inverse;
                                    values[tail + 6] = daily_heterotrophic_respiration.daily_signed_g_c[cell] * area_inverse;
                                    values[tail + 7] = daily_soil_fire_methane_emission_g_c[cell] * area_inverse;
                                    const active_carbonate = daily_soil_carbonate_mol_per_m3[first_layer .. first_layer + active_layers];
                                    const active_bicarbonate = daily_soil_bicarbonate_mol_per_m3[first_layer .. first_layer + active_layers];
                                    const active_carbonate_complexes = daily_soil_carbonate_complexes_mol_per_m3[first_layer .. first_layer + active_layers];
                                    const active_calcite = daily_soil_calcite_mol_per_m3[first_layer .. first_layer + active_layers];
                                    for (0..active_layers) |local_layer| {
                                        const layer = first_layer + local_layer;
                                        const aqueous = initial_chemistry_state.aqueous[layer];
                                        active_carbonate[local_layer] = aqueous.carbonate;
                                        active_bicarbonate[local_layer] = aqueous.bicarbonate;
                                        active_carbonate_complexes[local_layer] =
                                            aqueous.calcium_carbonate + aqueous.calcium_bicarbonate +
                                            aqueous.magnesium_carbonate + aqueous.magnesium_bicarbonate +
                                            aqueous.sodium_carbonate;
                                        active_calcite[local_layer] = initial_chemistry_state.geochemistry_solids[layer].calcite_solid_mol_per_m3;
                                    }
                                    const litter_aqueous = surface_litter_chemistry_state.cells[cell];
                                    values[tail + 8] = try ecosys.soil_inorganic_carbon_storage.calculate(.{
                                        .carbon_dioxide_and_methane_g_c = gas_inventory_g[@intFromEnum(ecosys.gas_transport.Species.carbon_dioxide)] +
                                            gas_inventory_g[@intFromEnum(ecosys.gas_transport.Species.methane)],
                                        .soil_water_volume_m3 = state.matrix_liquid_water_m3[first_layer .. first_layer + active_layers],
                                        .soil_bulk_volume_m3 = soil_solver_property_state.matrix_bulk_volume_m3[first_layer .. first_layer + active_layers],
                                        .carbonate_mol_per_m3 = active_carbonate,
                                        .bicarbonate_mol_per_m3 = active_bicarbonate,
                                        .dissolved_carbonate_complexes_mol_per_m3 = active_carbonate_complexes,
                                        .calcite_mol_per_m3 = active_calcite,
                                        .litter_water_volume_m3 = surface_precipitation_state.litter_water_m3[cell],
                                        .litter_carbonate_mol_per_m3 = litter_aqueous.carbonate_mol_per_m3,
                                        .litter_bicarbonate_mol_per_m3 = litter_aqueous.bicarbonate_mol_per_m3,
                                        .litter_dissolved_carbonate_complexes_mol_per_m3 = 0,
                                        .litter_calcite_mol_per_m3 = litter_aqueous.salt_minerals.calcite_mol_per_m3,
                                        .litter_bulk_volume_m3 = surface_litter_geometry_state.expanded_total_volume_m3[cell],
                                    }) * area_inverse;
                                    for (resolved.soil_enabled, values, 0..) |enabled, value, choice| if (enabled and !std.math.isFinite(value)) {
                                        std.log.err("daily soil carbon choice is not yet backed by an authoritative translated owner: choice={d} name={s}", .{ choice + 1, daily_soil_carbon_catalog.variables[choice].name });
                                        return error.UntranslatedDailySoilCarbonChoice;
                                    };
                                    const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                    const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                    const file_name = try ecosys.output_record.buildCellFileName(allocator, grid_column, grid_row, current_year, carbon_editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_soil_carbon_bank.streams[cell].write(file_name, daily_soil_carbon_catalog.variables, selection, resolved.soil_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = 23 },
                                        .grid_column = grid_column,
                                        .grid_row = grid_row,
                                        .values = values,
                                    });
                                }
                            }
                            if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(carbon_editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                const canopy = if (detailed_canopy_state) |*value| value else return error.DailyPlantCarbonOutputRequiresCanopyState;
                                const roots = if (plant_root_state) |*value| value else return error.DailyPlantCarbonOutputRequiresRootState;
                                for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                    const plant = try canopy.plantIndex(cell, species);
                                    const branches = try canopy.branchRange(plant);
                                    const values = try daily_plant_carbon_bank.row(plant);
                                    const root_profile = values[20 .. 20 + config.soil_layers];
                                    @memset(root_profile, 0);
                                    var root_symbiont_carbon_g: f64 = 0;
                                    for (0..config.soil_layers) |layer| {
                                        for (0..root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| {
                                            const root = try roots.layerIndex(plant, domain, layer);
                                            root_profile[layer] += roots.mobile_carbon_g[root];
                                            root_symbiont_carbon_g +=
                                                roots.symbiont_structural_carbon_g_c[root] +
                                                roots.symbiont_mobile_carbon_g_c[root];
                                            for (0..roots.active_root_axis_count[plant]) |axis| {
                                                const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
                                                root_profile[layer] += roots.axis_primary_carbon_g[axis_layer] + roots.axis_secondary_carbon_g[axis_layer];
                                            }
                                        }
                                    }
                                    var leaf_intermediate_carbon_g: f64 = 0;
                                    var canopy_symbiont_carbon_g: f64 = 0;
                                    var projected_leaf_area_m2: f64 = 0;
                                    for (branches.first..branches.end) |branch| {
                                        canopy_symbiont_carbon_g +=
                                            canopy.branch_symbiont_structural_carbon_g[branch] +
                                            canopy.branch_symbiont_mobile_carbon_g[branch];
                                        projected_leaf_area_m2 += canopy.branch_leaf_area_m2[branch];
                                        const nodes = try canopy.nodeRange(branch);
                                        for (nodes.first..nodes.end) |node| {
                                            leaf_intermediate_carbon_g +=
                                                canopy.node_c3_nonstructural_carbon_g[node] +
                                                canopy.node_c4_mesophyll_nonstructural_carbon_g[node] +
                                                canopy.node_bundle_sheath_co2_carbon_g[node] +
                                                canopy.node_bundle_sheath_bicarbonate_carbon_g[node];
                                        }
                                    }
                                    const primary_root_first = try roots.layerIndex(plant, 0, 0);
                                    const pools = try ecosys.plant_daily_pool_aggregation.calculateCarbonInto(.{
                                        .branch_leaf_carbon_g = canopy.branch_leaf_carbon_g[branches.first..branches.end],
                                        .branch_sheath_carbon_g = canopy.branch_sheath_carbon_g[branches.first..branches.end],
                                        .branch_stalk_carbon_g = canopy.branch_stalk_carbon_g[branches.first..branches.end],
                                        .branch_reserve_carbon_g = canopy.branch_reserve_carbon_g[branches.first..branches.end],
                                        .branch_husk_carbon_g = canopy.branch_husk_carbon_g[branches.first..branches.end],
                                        .branch_ear_carbon_g = canopy.branch_ear_carbon_g[branches.first..branches.end],
                                        .branch_grain_carbon_g = canopy.branch_grain_carbon_g[branches.first..branches.end],
                                        .branch_mobile_carbon_g = canopy.branch_mobile_carbon_g[branches.first..branches.end],
                                        .branch_seed_count = canopy.branch_seed_count[branches.first..branches.end],
                                        .c4_intermediate_carbon_g = leaf_intermediate_carbon_g,
                                        .canopy_symbiont_carbon_g = canopy_symbiont_carbon_g,
                                        .root_carbon_g_by_layer = root_profile,
                                        .primary_root_length_density_m_per_m3_by_layer = roots.root_length_density_m_per_m3[primary_root_first .. primary_root_first + config.soil_layers],
                                        .root_symbiont_carbon_g = root_symbiont_carbon_g,
                                        .standing_dead_carbon_g = canopy.plant_standing_dead_carbon_g[plant],
                                        .seed_storage_carbon_g = canopy.plant_seed_storage_carbon_g[plant],
                                        .projected_leaf_area_m2 = projected_leaf_area_m2,
                                        .plant_population_count = canopy.plant_population_count[plant],
                                    }, root_profile);
                                    const net_primary_productivity_g =
                                        plant_daily_flux_ledger.net_carbon_change_g[plant] +
                                        plant_daily_flux_ledger.signed_total_respiration_carbon_g[plant];
                                    const balance_carbon_g = try ecosys.plant_daily_output.carbonBalance(.{
                                        .shoot_carbon_g = pools.shoot_carbon_g,
                                        .root_carbon_g = pools.root_carbon_g,
                                        .nodule_carbon_g = pools.nodule_carbon_g,
                                        .storage_carbon_g = pools.storage_carbon_g,
                                        .standing_dead_carbon_g = pools.vegetative_residue_carbon_g,
                                        .cumulative_carbon_sink_g = plant_daily_flux_ledger.carbon_sink_g[plant],
                                        .cumulative_root_soil_carbon_exchange_g = plant_daily_flux_ledger.root_soil_carbon_exchange_g[plant],
                                        .cumulative_carbon_balance_g = plant_daily_flux_ledger.cumulative_carbon_balance_g[plant],
                                        .cumulative_harvested_carbon_g = plant_daily_flux_ledger.cumulative_harvested_carbon_g[plant],
                                        .harvested_carbon_g = plant_daily_flux_ledger.harvested_carbon_g[plant],
                                        .carbon_oxidation_g = plant_daily_flux_ledger.carbon_oxidation_g[plant],
                                        .cumulative_net_primary_productivity_g = net_primary_productivity_g,
                                    });
                                    try ecosys.plant_daily_output.calculateCarbonInto(.{
                                        .cell_area_m2 = canopy_cell_area_m2[cell],
                                        .shoot_carbon_g = pools.shoot_carbon_g,
                                        .leaf_carbon_g = pools.leaf_carbon_g,
                                        .sheath_carbon_g = pools.sheath_carbon_g,
                                        .stalk_carbon_g = pools.stalk_carbon_g,
                                        .reserve_carbon_g = pools.reserve_carbon_g,
                                        .husk_carbon_g = pools.husk_carbon_g,
                                        .ear_carbon_g = pools.ear_carbon_g,
                                        .grain_carbon_g = pools.grain_carbon_g,
                                        .root_carbon_g = pools.root_carbon_g,
                                        .nodule_carbon_g = pools.nodule_carbon_g,
                                        .vegetative_residue_carbon_g = pools.vegetative_residue_carbon_g,
                                        .grain_number = pools.grain_number,
                                        .projected_leaf_area_m2 = pools.projected_leaf_area_m2,
                                        .daily_net_carbon_change_g = plant_daily_flux_ledger.net_carbon_change_g[plant],
                                        .cumulative_carbon_uptake_g = plant_daily_flux_ledger.root_soil_carbon_exchange_g[plant],
                                        .cumulative_carbon_sink_g = plant_daily_flux_ledger.carbon_sink_g[plant],
                                        .initial_cumulative_carbon_sink_g = plant_daily_flux_ledger.initial_carbon_sink_g[plant],
                                        .signed_total_respiration_carbon_g = plant_daily_flux_ledger.signed_total_respiration_carbon_g[plant],
                                        .signed_aboveground_respiration_carbon_g = plant_daily_flux_ledger.signed_aboveground_respiration_carbon_g[plant],
                                        .carbon_pollination_factor = canopy.plant_mobile_carbon_concentration_g_per_g[plant],
                                        .harvested_carbon_g = plant_daily_flux_ledger.harvested_carbon_g[plant],
                                        .root_length_density_m_per_m3_by_layer = root_profile,
                                        .plant_population = canopy.plant_population_count[plant],
                                        .balance_carbon_g = balance_carbon_g,
                                        .storage_carbon_g = pools.storage_carbon_g,
                                        .carbon_oxidation_flux_g = plant_daily_flux_ledger.carbon_oxidation_g[plant],
                                        .net_primary_productivity_g = net_primary_productivity_g,
                                        .canopy_height_m = development_canopy_height_m[plant],
                                    }, values);
                                    const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                    const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                    const file_name = try ecosys.output_record.buildPlantFileName(allocator, grid_column, grid_row, species + 1, current_year, carbon_editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_plant_carbon_bank.streams[plant].write(file_name, daily_plant_carbon_catalog.variables, selection, resolved.plant_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = 23 },
                                        .grid_column = grid_column,
                                        .grid_row = grid_row,
                                        .values = values,
                                    });
                                };
                            }
                        }
                        const editor_index: usize = 6;
                        const editor_name = scene.output_editors[editor_index];
                        if (!ecosys.delimited_input.isNo(editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedDailySoilWaterEditor;
                            if (std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                const liquid_fraction = daily_water_profile_workspace[0..config.soil_layers];
                                const ice_fraction = daily_water_profile_workspace[config.soil_layers .. 2 * config.soil_layers];
                                const potential_count = @min(config.soil_layers, 10);
                                const matric_plus_osmotic_potential = daily_water_profile_workspace[2 * config.soil_layers ..][0..potential_count];
                                for (0..state.cell_count) |cell| {
                                    var soil_water_storage_m3: f64 = 0;
                                    for (0..config.soil_layers) |local_layer| {
                                        const layer = try state.layerIndex(cell, local_layer);
                                        const layer_volume_m3 = soil_solver_property_state.layer_volume_m3[layer];
                                        if (!std.math.isFinite(layer_volume_m3) or layer_volume_m3 <= 0) return error.InvalidOutputSoilLayerVolume;
                                        liquid_fraction[local_layer] = state.liquid_water_m3[layer] / layer_volume_m3;
                                        ice_fraction[local_layer] = state.ice_water_m3[layer] / layer_volume_m3;
                                        soil_water_storage_m3 += state.liquid_water_m3[layer] +
                                            state.water_vapor_volume_m3[layer] +
                                            state.ice_water_m3[layer];
                                        if (local_layer < potential_count) matric_plus_osmotic_potential[local_layer] = soil_hourly_workspace.matric_plus_osmotic_potential_mpa[layer];
                                    }
                                    const snow_base = cell * snow_transport_state.layer_capacity;
                                    for (0..snow_transport_state.layer_capacity) |snow_layer| {
                                        const snow = snow_base + snow_layer;
                                        soil_water_storage_m3 += snow_transport_state.solid_snow_water_equivalent_m3[snow] +
                                            snow_transport_state.liquid_water_volume_m3[snow] +
                                            snow_transport_state.vapor_water_equivalent_m3[snow] +
                                            snow_transport_state.ice_volume_m3[snow] * runscript.snow_ice_density_Mg_per_m3;
                                    }
                                    soil_water_storage_m3 +=
                                        surface_precipitation_state.litter_water_m3[cell] +
                                        surface_litter_ice_m3[cell];
                                    for (0..config.plant_populations) |species| {
                                        const plant = cell * config.plant_populations + species;
                                        soil_water_storage_m3 += plant_state.canopy_water_storage_m_per_m2[plant] * canopy_cell_area_m2[cell];
                                        if (canopy_precipitation_retention_state) |retention| {
                                            soil_water_storage_m3 += retention.living_surface_water_m3[plant] + retention.standing_dead_surface_water_m3[plant];
                                        }
                                    }
                                    const litter_capacity_m3 = surface_precipitation_state.litter_water_capacity_m3[cell];
                                    const surface_litter_thickness_m =
                                        surface_litter_geometry_state.expanded_total_volume_m3[cell] /
                                        canopy_cell_area_m2[cell];
                                    if (!std.math.isFinite(surface_litter_thickness_m) or
                                        surface_litter_thickness_m < 0)
                                        return error.InvalidOutputSurfaceLitterThickness;
                                    const active_boundary_base = cell * (soil_geometry_state.layer_capacity + 1);
                                    const active_surface_depth_m = soil_geometry_state.boundary_depth_m[active_boundary_base + soil_geometry_state.first_active_layer[cell]];
                                    const values = try daily_soil_water_bank.row(cell);
                                    try ecosys.soil_daily_output.calculateWaterInto(.{
                                        .cell_area_m2 = canopy_cell_area_m2[cell],
                                        .grid_area_m2 = total_grid_area_m2,
                                        .rainfall_m3 = daily_water_ledger.rainfall_m3[cell],
                                        .evaporation_m3 = daily_water_ledger.evaporation_m3[cell],
                                        .runoff_m3 = daily_water_ledger.runoff_m3[cell],
                                        .soil_water_storage_m3 = soil_water_storage_m3,
                                        .outflow_m3 = daily_water_ledger.water_outflow_m3[cell],
                                        .snow_depth_m = snow_depth_m[cell],
                                        .liquid_water_fraction_by_layer = liquid_fraction,
                                        .surface_excess_liquid_water_depth_m = @max(
                                            0,
                                            (surface_precipitation_state.litter_water_m3[cell] -
                                                litter_capacity_m3) /
                                                canopy_cell_area_m2[cell],
                                        ),
                                        .ice_fraction_by_layer = ice_fraction,
                                        .surface_excess_ice_water_depth_m = @max(
                                            0,
                                            (surface_litter_ice_m3[cell] -
                                                litter_capacity_m3) /
                                                canopy_cell_area_m2[cell],
                                        ),
                                        .matric_plus_osmotic_potential_mpa_by_layer = matric_plus_osmotic_potential,
                                        .surface_matric_potential_mpa = surface_litter_water_environment_state.matric_water_potential_mpa[cell],
                                        .lateral_water_outflow_m3 = daily_water_ledger.lateral_water_outflow_m3[cell],
                                        .sediment_outflow_m3 = daily_water_ledger.sediment_outflow_m3[cell],
                                        .mineral_soil_surface_depth_m = active_surface_depth_m,
                                        .surface_litter_thickness_m = surface_litter_thickness_m,
                                        .active_layer_bottom_depth_m = soil_boundary_topology_state.active_layer_depth_m[cell],
                                        .water_table_depth_m = soil_boundary_topology_state.internal_water_table_depth_m[cell],
                                    }, values);
                                    const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                    const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                    const file_name = try ecosys.output_record.buildCellFileName(allocator, grid_column, grid_row, current_year, editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_soil_water_bank.streams[cell].write(file_name, daily_soil_water_catalog.variables, selection, resolved.soil_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = 23 },
                                        .grid_column = grid_column,
                                        .grid_row = grid_row,
                                        .values = values,
                                    });
                                }
                            }
                        }
                        if (!ecosys.delimited_input.isNo(editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + editor_index] orelse return error.UnresolvedDailyPlantWaterEditor;
                            if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                const growth = if (plant_growth_stage_state) |*value| value else return error.DailyPlantWaterOutputRequiresGrowthStageState;
                                const dormancy = if (plant_dormancy_state) |*value| value else return error.DailyPlantWaterOutputRequiresDormancyState;
                                for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                    const plant = cell * config.plant_populations + species;
                                    const output = try ecosys.plant_daily_output.water(
                                        plant_daily_flux_ledger.transpiration_source_m3[plant],
                                        try ecosys.plant_development.coldOrWaterStressHours(growth, dormancy, plant),
                                        phenology_root_oxygen_fraction[plant],
                                        plant_state.canopy_water_storage_m_per_m2[plant] * canopy_cell_area_m2[cell],
                                        canopy_cell_area_m2[cell],
                                    );
                                    const values = try daily_plant_water_bank.row(plant);
                                    values[0] = output.transpiration_mm;
                                    values[1] = output.cold_or_water_stress_h;
                                    values[2] = output.oxygen_stress_factor;
                                    values[3] = output.plant_water_storage_mm;
                                    const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                    const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                    const file_name = try ecosys.output_record.buildPlantFileName(allocator, grid_column, grid_row, species + 1, current_year, editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_plant_water_bank.streams[plant].write(file_name, daily_plant_water_catalog.variables, selection, resolved.plant_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = 23 },
                                        .grid_column = grid_column,
                                        .grid_row = grid_row,
                                        .values = values,
                                    });
                                };
                            }
                        }
                        const soil_nitrogen_editor_index: usize = 7;
                        const soil_nitrogen_editor_name = scene.output_editors[soil_nitrogen_editor_index];
                        if (!ecosys.delimited_input.isNo(soil_nitrogen_editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + soil_nitrogen_editor_index] orelse return error.UnresolvedDailySoilNitrogenEditor;
                            if (std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(soil_nitrogen_editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                for (0..state.cell_count) |cell| {
                                    const values = try daily_soil_nitrogen_bank.row(cell);
                                    @memset(values, std.math.nan(f64));
                                    const first_layer = try state.layerIndex(cell, 0);
                                    const active_layers = state.active_soil_layer_count[cell];
                                    const pools = try ecosys.soil_daily_nitrogen_pools.calculate(
                                        &soil_organic_state,
                                        &surface_organic_state,
                                        first_layer,
                                        active_layers,
                                        cell,
                                    );
                                    const area_inverse = 1.0 / canopy_cell_area_m2[cell];
                                    values[0] = pools.residue_nitrogen_g_n * area_inverse;
                                    values[1] = pools.humus_nitrogen_g_n * area_inverse;
                                    values[2] = (soil_fertilizer_inventory.daily_nitrogen_input_g_n[cell] + daily_organic_fertilizer_nitrogen_input_g_n[cell]) * area_inverse;
                                    var plant_nitrogen_sink_g_n: f64 = 0;
                                    var harvested_nitrogen_g_n: f64 = 0;
                                    for (0..config.plant_populations) |species| {
                                        const plant = cell * config.plant_populations + species;
                                        plant_nitrogen_sink_g_n += plant_daily_flux_ledger.nitrogen_sink_g[plant];
                                        harvested_nitrogen_g_n += plant_daily_flux_ledger.harvested_nitrogen_g[plant];
                                    }
                                    values[3] = plant_nitrogen_sink_g_n * area_inverse;
                                    values[6] = daily_nitrogen_export.dissolved_organic_nitrogen_runoff_g_n[cell] / total_grid_area_m2;
                                    values[7] = daily_nitrogen_export.dissolved_organic_nitrogen_drainage_g_n[cell] / total_grid_area_m2;
                                    values[8] = daily_nitrogen_export.dissolved_inorganic_nitrogen_runoff_g_n[cell] / total_grid_area_m2;
                                    values[9] = daily_nitrogen_export.dissolved_inorganic_nitrogen_drainage_g_n[cell] / total_grid_area_m2;
                                    var dissolved_dinitrogen_g_n: f64 = 0;
                                    for (0..active_layers) |local_layer| {
                                        const layer = first_layer + local_layer;
                                        for (0..ecosys.gas_transport.species_count) |species| {
                                            const component = layer * ecosys.gas_transport.species_count + species;
                                            if (species == @intFromEnum(ecosys.gas_transport.Species.nitrogen)) {
                                                dissolved_dinitrogen_g_n +=
                                                    gas_transport_state.dissolved_mass_g[component] +
                                                    gas_transport_state.macropore_dissolved_mass_g[component] +
                                                    gas_transport_state.band_dissolved_mass_g[component];
                                            }
                                        }
                                    }
                                    for (0..ecosys.gas_transport.species_count) |species| {
                                        const component = cell * ecosys.gas_transport.species_count + species;
                                        if (species == @intFromEnum(ecosys.gas_transport.Species.nitrogen)) {
                                            dissolved_dinitrogen_g_n +=
                                                litter_gas_transport_state.dissolved_mass_g[component] +
                                                litter_gas_transport_state.macropore_dissolved_mass_g[component] +
                                                litter_gas_transport_state.band_dissolved_mass_g[component];
                                        }
                                    }
                                    values[10] = (try daily_soil_gas_flux.getSoilLitterBoundary(cell, .nitrous_oxide)) * area_inverse;
                                    values[11] = (try daily_soil_gas_flux.getSoilLitterBoundary(cell, .ammonia)) * area_inverse;
                                    values[12] = dissolved_dinitrogen_g_n * area_inverse;
                                    values[13] = (pools.residue_nitrogen_g_n + pools.humus_nitrogen_g_n) * area_inverse;
                                    const ammonium_profile = values[14 .. 14 + config.soil_layers];
                                    const oxidized_profile = values[14 + config.soil_layers .. 14 + 2 * config.soil_layers];
                                    const mineral = try ecosys.soil_daily_mineral_nitrogen.calculateInto(
                                        &initial_chemistry_state,
                                        &soil_reactive_nitrogen_state,
                                        first_layer,
                                        active_layers,
                                        state.matrix_liquid_water_m3,
                                        soil_solver_property_state.matrix_bulk_volume_m3,
                                        soil_solver_property_state.bulk_density_megagrams_per_m3,
                                        .{
                                            .ammonium_non_band = 1 - runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                                            .ammonium_band = runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                                            .nitrate_non_band = 1 - runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                                            .nitrate_band = runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                                        },
                                        runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
                                        config.absolute_tolerance,
                                        ammonium_profile,
                                        oxidized_profile,
                                    );
                                    values[4] = mineral.ammonium_nitrogen_g_n * area_inverse;
                                    values[5] = mineral.nitrate_plus_nitrite_nitrogen_g_n * area_inverse;
                                    const tail = 14 + 2 * config.soil_layers;
                                    const litter = surface_litter_chemistry_state.cells[cell];
                                    const litter_bulk_volume_m3 = surface_litter_geometry_state.expanded_total_volume_m3[cell];
                                    const litter_ammonium_g_n = runscript.fertilizer_nitrogen_molar_mass_g_per_mol *
                                        (litter.ammonium_mol_per_m3 * surface_precipitation_state.litter_water_m3[cell] +
                                            litter.exchange.ammonium_mol_per_Mg * surface_litter_geometry_state.dry_mass_Mg[cell]);
                                    values[tail] = if (litter_bulk_volume_m3 > config.absolute_tolerance) litter_ammonium_g_n / litter_bulk_volume_m3 else 0;
                                    values[tail + 1] = daily_soil_combusted_nitrogen_g_n[cell] * area_inverse;
                                    values[tail + 2] = harvested_nitrogen_g_n * area_inverse;
                                    values[tail + 3] = daily_microbial_nitrogen_mineralization_g_n[cell] * area_inverse;
                                    values[tail + 4] = daily_soil_fire_nitrogen_flux_g_n[cell] * area_inverse;
                                    values[tail + 5] = (try daily_soil_gas_flux.getSoilLitterBoundary(cell, .nitrogen)) * area_inverse;
                                    for (resolved.soil_enabled, values, 0..) |enabled, value, choice| if (enabled and !std.math.isFinite(value)) {
                                        std.log.err("daily soil nitrogen choice has non-finite authoritative state: choice={d} name={s}", .{ choice + 1, daily_soil_nitrogen_catalog.variables[choice].name });
                                        return error.NonFiniteDailySoilNitrogenChoice;
                                    };
                                    const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                    const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                    const file_name = try ecosys.output_record.buildCellFileName(allocator, grid_column, grid_row, current_year, soil_nitrogen_editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_soil_nitrogen_bank.streams[cell].write(file_name, daily_soil_nitrogen_catalog.variables, selection, resolved.soil_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = 23 },
                                        .grid_column = grid_column,
                                        .grid_row = grid_row,
                                        .values = values,
                                    });
                                }
                            }
                        }
                        const soil_phosphorus_editor_index: usize = 8;
                        const soil_phosphorus_editor_name = scene.output_editors[soil_phosphorus_editor_index];
                        if (!ecosys.delimited_input.isNo(soil_phosphorus_editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + soil_phosphorus_editor_index] orelse return error.UnresolvedDailySoilPhosphorusEditor;
                            if (std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(soil_phosphorus_editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                for (0..state.cell_count) |cell| {
                                    const values = try daily_soil_phosphorus_bank.row(cell);
                                    @memset(values, std.math.nan(f64));
                                    const area_inverse = 1.0 / canopy_cell_area_m2[cell];
                                    values[2] = (mineral_fertilizer_inventory.daily_phosphorus_input_g_p[cell] + daily_organic_fertilizer_phosphorus_input_g_p[cell]) * area_inverse;
                                    var plant_litterfall_phosphorus_g_p: f64 = 0;
                                    for (0..config.plant_populations) |species| plant_litterfall_phosphorus_g_p += plant_daily_flux_ledger.phosphorus_sink_g[cell * config.plant_populations + species];
                                    values[3] = plant_litterfall_phosphorus_g_p * area_inverse;
                                    const first_layer = try state.layerIndex(cell, 0);
                                    const active_layers = state.active_soil_layer_count[cell];
                                    const organic_pools = try ecosys.soil_daily_phosphorus_pools.calculate(
                                        &soil_organic_state,
                                        &surface_organic_state,
                                        first_layer,
                                        active_layers,
                                        cell,
                                    );
                                    values[0] = organic_pools.residue_phosphorus_g_p * area_inverse;
                                    values[1] = organic_pools.humus_phosphorus_g_p * area_inverse;
                                    values[5] = daily_phosphorus_export.dissolved_organic_phosphorus_runoff_g_p[cell] / total_grid_area_m2;
                                    values[6] = daily_phosphorus_export.dissolved_organic_phosphorus_drainage_g_p[cell] / total_grid_area_m2;
                                    values[7] = daily_phosphorus_export.dissolved_inorganic_phosphorus_runoff_g_p[cell] / total_grid_area_m2;
                                    values[8] = daily_phosphorus_export.dissolved_inorganic_phosphorus_drainage_g_p[cell] / total_grid_area_m2;
                                    const phosphate_band_fraction = runscript.plant_nutrient_initialization.initial_phosphate_band_fraction;
                                    const phosphate_non_band_fraction = 1.0 - phosphate_band_fraction;
                                    values[4] = try ecosys.soil_phosphate_inventory.exchangeablePhosphorus_g_p(
                                        &state,
                                        &initial_chemistry_state,
                                        &surface_litter_chemistry_state,
                                        soil_solver_property_state.matrix_bulk_volume_m3,
                                        soil_solver_property_state.bulk_density_megagrams_per_m3,
                                        surface_litter_geometry_state.dry_mass_Mg,
                                        cell,
                                        phosphate_non_band_fraction,
                                        phosphate_band_fraction,
                                        runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                                    ) * area_inverse;
                                    values[9] = try ecosys.soil_phosphate_inventory.precipitatedPhosphorus_g_p(
                                        &state,
                                        &initial_chemistry_state,
                                        &surface_litter_chemistry_state,
                                        state.matrix_liquid_water_m3,
                                        surface_precipitation_state.litter_water_m3,
                                        cell,
                                        phosphate_non_band_fraction,
                                        phosphate_band_fraction,
                                        runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                                    ) * area_inverse;
                                    values[10] = organic_pools.microbial_phosphorus_g_p * area_inverse;
                                    values[11] = try ecosys.soil_daily_fire_phosphorus.sourceSignedFlux_g_p(
                                        cell,
                                        daily_soil_fire_phosphorus_flux_g_p,
                                        plant_daily_flux_ledger.phosphorus_oxidation_g,
                                        config.plant_populations,
                                    ) * area_inverse;
                                    for (0..config.soil_layers) |local_layer| {
                                        const aqueous_index = 12 + local_layer;
                                        const sorbed_index = 12 + config.soil_layers + local_layer;
                                        if (local_layer >= active_layers) {
                                            values[aqueous_index] = 0;
                                            values[sorbed_index] = 0;
                                            continue;
                                        }
                                        const layer = first_layer + local_layer;
                                        const non_band = initial_chemistry_state.non_band_phosphate[layer];
                                        const band = initial_chemistry_state.band_phosphate[layer];
                                        values[aqueous_index] = 31.0 * (phosphate_non_band_fraction * (non_band.dissolved_hpo4_mol_p_per_m3 + non_band.dissolved_h2po4_mol_p_per_m3) + phosphate_band_fraction * (band.dissolved_hpo4_mol_p_per_m3 + band.dissolved_h2po4_mol_p_per_m3));
                                        values[sorbed_index] = try ecosys.soil_daily_output.sorbedPhosphorusConcentrationGPerM3(
                                            phosphate_non_band_fraction *
                                                (non_band.adsorbed_hpo4_mol_p_per_Mg +
                                                    non_band.adsorbed_h2po4_mol_p_per_Mg) +
                                                phosphate_band_fraction *
                                                    (band.adsorbed_hpo4_mol_p_per_Mg +
                                                        band.adsorbed_h2po4_mol_p_per_Mg),
                                            soil_solver_property_state.bulk_density_megagrams_per_m3[layer],
                                            runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                                        );
                                    }
                                    const after_layers = 12 + 2 * config.soil_layers;
                                    const litter = surface_litter_chemistry_state.cells[cell];
                                    values[after_layers] = 31.0 * (litter.hpo4_mol_p_per_m3 + litter.h2po4_mol_p_per_m3);
                                    values[after_layers + 1] = if (surface_litter_geometry_state.expanded_total_volume_m3[cell] >
                                        config.absolute_tolerance)
                                        try ecosys.soil_daily_output.sorbedPhosphorusConcentrationGPerM3(
                                            litter.phosphate_surface.adsorbed_hpo4_mol_p_per_Mg +
                                                litter.phosphate_surface.adsorbed_h2po4_mol_p_per_Mg,
                                            surface_litter_geometry_state.dry_mass_Mg[cell] /
                                                surface_litter_geometry_state.expanded_total_volume_m3[cell],
                                            runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                                        )
                                    else
                                        0;
                                    const final_fixed = after_layers + 2;
                                    values[final_fixed] = daily_soil_combusted_phosphorus_g_p[cell] * area_inverse;
                                    values[final_fixed + 1] = try ecosys.soil_phosphate_inventory.solublePhosphorus_g_p(
                                        &state,
                                        &initial_chemistry_state,
                                        &surface_litter_chemistry_state,
                                        state.matrix_liquid_water_m3,
                                        surface_precipitation_state.litter_water_m3,
                                        cell,
                                        phosphate_non_band_fraction,
                                        phosphate_band_fraction,
                                        runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
                                    ) * area_inverse;
                                    var harvested_phosphorus_g_p: f64 = 0;
                                    for (0..config.plant_populations) |species| harvested_phosphorus_g_p += plant_daily_flux_ledger.harvested_phosphorus_g[cell * config.plant_populations + species];
                                    values[final_fixed + 2] = harvested_phosphorus_g_p * area_inverse;
                                    values[final_fixed + 3] = daily_microbial_phosphate_mineralization_g_p[cell] * area_inverse;
                                    // OUTSD phosphorus choices 49 and 50 are source-reserved.
                                    values[values.len - 2] = 0;
                                    values[values.len - 1] = 0;
                                    for (resolved.soil_enabled, values, 0..) |enabled, value, choice| if (enabled and !std.math.isFinite(value)) {
                                        std.log.err("daily soil phosphorus choice is not yet backed by an authoritative translated owner: choice={d} name={s}", .{ choice + 1, daily_soil_phosphorus_catalog.variables[choice].name });
                                        return error.UntranslatedDailySoilPhosphorusChoice;
                                    };
                                    const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                    const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                    const file_name = try ecosys.output_record.buildCellFileName(allocator, grid_column, grid_row, current_year, soil_phosphorus_editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_soil_phosphorus_bank.streams[cell].write(file_name, daily_soil_phosphorus_catalog.variables, selection, resolved.soil_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = 23 },
                                        .grid_column = grid_column,
                                        .grid_row = grid_row,
                                        .values = values,
                                    });
                                }
                            }
                        }
                        const soil_heat_editor_index: usize = 9;
                        const soil_heat_editor_name = scene.output_editors[soil_heat_editor_index];
                        if (!ecosys.delimited_input.isNo(soil_heat_editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + soil_heat_editor_index] orelse return error.UnresolvedDailySoilHeatEditor;
                            if (std.mem.indexOfScalar(bool, resolved.soil_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(soil_heat_editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                const maximum_temperature_c = daily_soil_heat_workspace[0..config.soil_layers];
                                const minimum_temperature_c = daily_soil_heat_workspace[config.soil_layers .. 2 * config.soil_layers];
                                const electrical_conductivity_dS_per_m = daily_soil_heat_workspace[2 * config.soil_layers .. 3 * config.soil_layers];
                                const charge_fractions = ecosys.solute_charge_classification.ZoneFractions{
                                    .ammonium_non_band = 1 - runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                                    .ammonium_band = runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                                    .nitrate_non_band = 1 - runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                                    .nitrate_band = runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                                    .phosphate_non_band = 1 - runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
                                    .phosphate_band = runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
                                };
                                for (0..state.cell_count) |cell| {
                                    try daily_heat_ledger.validateCell(cell);
                                    const first_layer = try state.layerIndex(cell, 0);
                                    const active_layers = state.active_soil_layer_count[cell];
                                    const ledger_first = cell * config.soil_layers;
                                    @memset(maximum_temperature_c, 0);
                                    @memset(minimum_temperature_c, 0);
                                    @memset(electrical_conductivity_dS_per_m, 0);
                                    @memcpy(maximum_temperature_c[0..active_layers], daily_heat_ledger.maximum_soil_temperature_c[ledger_first .. ledger_first + active_layers]);
                                    @memcpy(minimum_temperature_c[0..active_layers], daily_heat_ledger.minimum_soil_temperature_c[ledger_first .. ledger_first + active_layers]);
                                    for (0..active_layers) |local_layer| {
                                        const coefficients = try initial_chemistry_state.activityCoefficients(first_layer + local_layer, charge_fractions);
                                        electrical_conductivity_dS_per_m[local_layer] = coefficients.electrical_conductivity_dS_per_m;
                                    }
                                    const values = try daily_soil_heat_bank.row(cell);
                                    try ecosys.soil_daily_output.calculateHeatInto(.{
                                        .total_radiation_mj_m2 = daily_heat_ledger.total_radiation_mj_per_m2[cell],
                                        .maximum_air_temperature_c = daily_heat_ledger.maximum_air_temperature_c[cell],
                                        .minimum_air_temperature_c = daily_heat_ledger.minimum_air_temperature_c[cell],
                                        .maximum_atmospheric_vapor_pressure_kpa = daily_heat_ledger.maximum_vapor_pressure_kpa[cell],
                                        .minimum_atmospheric_vapor_pressure_kpa = daily_heat_ledger.minimum_vapor_pressure_kpa[cell],
                                        .cumulative_wind_m = daily_heat_ledger.cumulative_wind_distance_m[cell],
                                        .total_precipitation_mm = daily_heat_ledger.total_precipitation_mm[cell],
                                        .maximum_soil_temperature_c_by_layer = maximum_temperature_c,
                                        .minimum_soil_temperature_c_by_layer = minimum_temperature_c,
                                        .surface_maximum_soil_temperature_c = daily_heat_ledger.surface_maximum_temperature_c[cell],
                                        .surface_minimum_soil_temperature_c = daily_heat_ledger.surface_minimum_temperature_c[cell],
                                        .electrical_conductivity_by_layer = electrical_conductivity_dS_per_m,
                                        .ionic_outflow_mol = daily_heat_ledger.ionic_outflow_mol[cell],
                                        .grid_area_m2 = total_grid_area_m2,
                                    }, values);
                                    const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                    const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                    const file_name = try ecosys.output_record.buildCellFileName(allocator, grid_column, grid_row, current_year, soil_heat_editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_soil_heat_bank.streams[cell].write(file_name, daily_soil_heat_catalog.variables, selection, resolved.soil_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = 23 },
                                        .grid_column = grid_column,
                                        .grid_row = grid_row,
                                        .values = values,
                                    });
                                }
                            }
                        }
                        inline for (.{ DailyPlantElement.nitrogen, DailyPlantElement.phosphorus }) |element| {
                            const nutrient_editor_index: usize = switch (element) {
                                .nitrogen => 7,
                                .phosphorus => 8,
                            };
                            const nutrient_editor_name = scene.output_editors[nutrient_editor_index];
                            if (!ecosys.delimited_input.isNo(nutrient_editor_name)) {
                                const resolved = resolved_output_editors[pass.scene_index * 10 + nutrient_editor_index] orelse return error.UnresolvedDailyPlantNutrientEditor;
                                if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                                    const selection_index = output_selection_catalog.find(nutrient_editor_name) orelse return error.MissingOutputSelection;
                                    const selection = output_selection_catalog.entries.items[selection_index].selection;
                                    const canopy = if (detailed_canopy_state) |*value| value else return error.DailyPlantNutrientOutputRequiresCanopyState;
                                    const roots = if (plant_root_state) |*value| value else return error.DailyPlantNutrientOutputRequiresRootState;
                                    for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                        const plant = try canopy.plantIndex(cell, species);
                                        const branches = try canopy.branchRange(plant);
                                        const root_workspace = daily_element_root_workspace[plant * config.soil_layers .. (plant + 1) * config.soil_layers];
                                        const pools = try calculateDailyPlantElementPools(
                                            canopy,
                                            roots,
                                            plant,
                                            branches,
                                            element,
                                            root_metabolism_plant_parameters[plant].biologicalDomainCount(),
                                            root_workspace,
                                        );
                                        var leaf_carbon_g: f64 = 0;
                                        var leaf_nitrogen_g: f64 = 0;
                                        var leaf_phosphorus_g: f64 = 0;
                                        for (branches.first..branches.end) |branch| {
                                            leaf_carbon_g += canopy.branch_leaf_carbon_g[branch];
                                            leaf_nitrogen_g += canopy.branch_leaf_nitrogen_g[branch];
                                            leaf_phosphorus_g += canopy.branch_leaf_phosphorus_g[branch];
                                        }
                                        const balance_inputs: ecosys.plant_daily_output.NutrientBalanceInputs = switch (element) {
                                            .nitrogen => .{
                                                .shoot_g = pools.shoot_g,
                                                .root_g = pools.root_g,
                                                .nodule_g = pools.nodule_g,
                                                .storage_g = pools.storage_g,
                                                .standing_dead_g = pools.vegetative_residue_g,
                                                .cumulative_sink_g = plant_daily_flux_ledger.nitrogen_sink_g[plant],
                                                .cumulative_root_soil_exchange_g = plant_daily_flux_ledger.root_soil_nitrogen_exchange_g[plant],
                                                .cumulative_balance_g = plant_daily_flux_ledger.cumulative_nitrogen_balance_g[plant],
                                                .cumulative_harvested_g = plant_daily_flux_ledger.cumulative_harvested_nitrogen_g[plant],
                                                .harvested_g = plant_daily_flux_ledger.harvested_nitrogen_g[plant],
                                                .oxidation_g = plant_daily_flux_ledger.nitrogen_oxidation_g[plant],
                                                .atmospheric_exchange_g = plant_daily_flux_ledger.ammonia_exchange_g_n[plant],
                                                .biological_fixation_g = plant_daily_flux_ledger.symbiotic_nitrogen_fixation_g[plant],
                                            },
                                            .phosphorus => .{
                                                .shoot_g = pools.shoot_g,
                                                .root_g = pools.root_g,
                                                .nodule_g = pools.nodule_g,
                                                .storage_g = pools.storage_g,
                                                .standing_dead_g = pools.vegetative_residue_g,
                                                .cumulative_sink_g = plant_daily_flux_ledger.phosphorus_sink_g[plant],
                                                .cumulative_root_soil_exchange_g = plant_daily_flux_ledger.root_soil_phosphorus_exchange_g[plant],
                                                .cumulative_balance_g = plant_daily_flux_ledger.cumulative_phosphorus_balance_g[plant],
                                                .cumulative_harvested_g = plant_daily_flux_ledger.cumulative_harvested_phosphorus_g[plant],
                                                .harvested_g = plant_daily_flux_ledger.harvested_phosphorus_g[plant],
                                                .oxidation_g = plant_daily_flux_ledger.phosphorus_oxidation_g[plant],
                                            },
                                        };
                                        const nutrient_inputs: ecosys.plant_daily_output.NutrientInputs = .{
                                            .cell_area_m2 = canopy_cell_area_m2[cell],
                                            .shoot_g = pools.shoot_g,
                                            .leaf_g = pools.leaf_g,
                                            .sheath_g = pools.sheath_g,
                                            .stalk_g = pools.stalk_g,
                                            .reserve_g = pools.reserve_g,
                                            .husk_g = pools.husk_g,
                                            .ear_g = pools.ear_g,
                                            .grain_g = pools.grain_g,
                                            .root_g = pools.root_g,
                                            .nodule_g = pools.nodule_g,
                                            .vegetative_residue_g = pools.vegetative_residue_g,
                                            .cumulative_uptake_g = switch (element) {
                                                .nitrogen => plant_daily_flux_ledger.root_soil_nitrogen_exchange_g[plant],
                                                .phosphorus => plant_daily_flux_ledger.root_soil_phosphorus_exchange_g[plant],
                                            },
                                            .cumulative_sink_g = switch (element) {
                                                .nitrogen => plant_daily_flux_ledger.nitrogen_sink_g[plant],
                                                .phosphorus => plant_daily_flux_ledger.phosphorus_sink_g[plant],
                                            },
                                            .nitrogen_fixation_g = plant_daily_flux_ledger.symbiotic_nitrogen_fixation_g[plant],
                                            .aboveground_litter_sink_g = switch (element) {
                                                .nitrogen => plant_daily_flux_ledger.initial_nitrogen_sink_g[plant],
                                                .phosphorus => plant_daily_flux_ledger.initial_phosphorus_sink_g[plant],
                                            },
                                            .pollination_factor = switch (element) {
                                                .nitrogen => canopy.plant_mobile_nitrogen_concentration_g_per_g[plant],
                                                .phosphorus => canopy.plant_mobile_phosphorus_concentration_g_per_g[plant],
                                            },
                                            .leaf_carbon_g = leaf_carbon_g,
                                            .leaf_nitrogen_g = leaf_nitrogen_g,
                                            .leaf_phosphorus_g = leaf_phosphorus_g,
                                            .nonstructural_carbon_g = canopy.plant_mobile_carbon_g[plant],
                                            .nonstructural_nitrogen_g = canopy.plant_mobile_nitrogen_g[plant],
                                            .nonstructural_phosphorus_g = canopy.plant_mobile_phosphorus_g[plant],
                                            .minimum_leaf_carbon_g = runscript.plant_pool_parameters.branch_structural_presence_g_per_plant * canopy.plant_population_count[plant],
                                            .other_pollination_factor = canopy.plant_mobile_phosphorus_concentration_g_per_g[plant],
                                            .gaseous_loss_g = plant_daily_flux_ledger.ammonia_exchange_g_n[plant],
                                            .harvested_g = switch (element) {
                                                .nitrogen => plant_daily_flux_ledger.harvested_nitrogen_g[plant],
                                                .phosphorus => plant_daily_flux_ledger.harvested_phosphorus_g[plant],
                                            },
                                            .balance_g = switch (element) {
                                                .nitrogen => try ecosys.plant_daily_output.nitrogenBalance(balance_inputs),
                                                .phosphorus => try ecosys.plant_daily_output.phosphorusBalance(balance_inputs),
                                            },
                                            .storage_g = pools.storage_g,
                                            .oxidation_flux_g = switch (element) {
                                                .nitrogen => plant_daily_flux_ledger.nitrogen_oxidation_g[plant],
                                                .phosphorus => plant_daily_flux_ledger.phosphorus_oxidation_g[plant],
                                            },
                                        };
                                        const values = switch (element) {
                                            .nitrogen => try daily_plant_nitrogen_bank.row(plant),
                                            .phosphorus => try daily_plant_phosphorus_bank.row(plant),
                                        };
                                        switch (element) {
                                            .nitrogen => @memcpy(values, &(try ecosys.plant_daily_output.nitrogen(nutrient_inputs)).values),
                                            .phosphorus => @memcpy(values, &(try ecosys.plant_daily_output.phosphorus(nutrient_inputs)).values),
                                        }
                                        const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                        const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                        const file_name = try ecosys.output_record.buildPlantFileName(allocator, grid_column, grid_row, species + 1, current_year, nutrient_editor_name);
                                        defer allocator.free(file_name);
                                        _ = switch (element) {
                                            .nitrogen => try daily_plant_nitrogen_bank.streams[plant].write(file_name, daily_plant_nitrogen_catalog.variables, selection, resolved.plant_enabled, .{
                                                .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = 23 },
                                                .grid_column = grid_column,
                                                .grid_row = grid_row,
                                                .values = values,
                                            }),
                                            .phosphorus => try daily_plant_phosphorus_bank.streams[plant].write(file_name, daily_plant_phosphorus_catalog.variables, selection, resolved.plant_enabled, .{
                                                .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = 23 },
                                                .grid_column = grid_column,
                                                .grid_row = grid_row,
                                                .values = values,
                                            }),
                                        };
                                    };
                                }
                            }
                        }
                        const development_editor_index: usize = 9;
                        const development_editor_name = scene.output_editors[development_editor_index];
                        if (!ecosys.delimited_input.isNo(development_editor_name)) {
                            const resolved = resolved_output_editors[pass.scene_index * 10 + development_editor_index] orelse return error.UnresolvedDailyPlantDevelopmentEditor;
                            if (std.mem.indexOfScalar(bool, resolved.plant_enabled, true) != null) {
                                const selection_index = output_selection_catalog.find(development_editor_name) orelse return error.MissingOutputSelection;
                                const selection = output_selection_catalog.entries.items[selection_index].selection;
                                const canopy = if (detailed_canopy_state) |*value| value else return error.DailyPlantDevelopmentOutputRequiresCanopyState;
                                const growth = if (plant_growth_stage_state) |*value| value else return error.DailyPlantDevelopmentOutputRequiresGrowthStageState;
                                const development_state = if (branch_development_state) |*value| value else return error.DailyPlantDevelopmentOutputRequiresBranchDevelopmentState;
                                const phenology = if (plant_phenology_state) |*value| value else return error.DailyPlantDevelopmentOutputRequiresPhenologyState;
                                for (0..state.cell_count) |cell| for (0..config.plant_populations) |species| {
                                    const plant = cell * config.plant_populations + species;
                                    const branches = try growth.branchRange(plant);
                                    const has_branch = branches.first < branches.end;
                                    const main_branch = (try growth.mainLivingBranch(plant)) orelse branches.first;
                                    if (has_branch and (main_branch >= canopy.branch_leaf_carbon_g.len or main_branch >= development_state.branch_count)) return error.DailyPlantDevelopmentBranchDimensionMismatch;
                                    const stage: ecosys.plant_growth_stages.BranchState = if (has_branch) growth.branches[main_branch] else .{};
                                    const milestones = [_]u32{
                                        stage.emergence_day,
                                        stage.floral_initiation_day,
                                        stage.stem_elongation_start_day,
                                        stage.stem_elongation_midpoint_day,
                                        stage.stem_elongation_end_day,
                                        stage.anthesis_day,
                                        stage.grain_fill_start_day,
                                        stage.seed_number_set_end_day,
                                        stage.seed_size_set_end_day,
                                        if (has_branch) development_state.stage_day[main_branch * 10 + 9] else 0,
                                    };
                                    var leaf_carbon_g: f64 = 0;
                                    var leaf_nitrogen_g: f64 = 0;
                                    var leaf_phosphorus_g: f64 = 0;
                                    for (branches.first..branches.end) |branch| {
                                        if (growth.branches[branch].dead) continue;
                                        leaf_carbon_g += canopy.branch_leaf_carbon_g[branch];
                                        leaf_nitrogen_g += canopy.branch_leaf_nitrogen_g[branch];
                                        leaf_phosphorus_g += canopy.branch_leaf_phosphorus_g[branch];
                                    }
                                    const output = try ecosys.plant_daily_output.development(.{
                                        .alive = has_branch and phenology.active[plant],
                                        .milestone_day_by_stage = &milestones,
                                        .branch_count = branches.end - branches.first,
                                        .main_branch_stage = stage.accumulated_vegetative_stage,
                                        .development_feedback = canopy.plant_nitrogen_phosphorus_fixation_constraint_fraction[plant],
                                        .leaf_structural_carbon_g = leaf_carbon_g,
                                        .leaf_structural_nitrogen_g = leaf_nitrogen_g,
                                        .leaf_structural_phosphorus_g = leaf_phosphorus_g,
                                        .plant_nonstructural_carbon_g = canopy.plant_mobile_carbon_g[plant],
                                        .plant_nonstructural_nitrogen_g = canopy.plant_mobile_nitrogen_g[plant],
                                        .plant_nonstructural_phosphorus_g = canopy.plant_mobile_phosphorus_g[plant],
                                        .minimum_leaf_carbon_g = runscript.plant_pool_parameters.branch_structural_presence_g_per_plant * canopy.plant_population_count[plant],
                                        .minimum_daily_canopy_water_potential_mpa = canopy.plant_minimum_daily_canopy_water_potential_mpa[plant],
                                        .oxygen_stress_factor = phenology_root_oxygen_fraction[plant],
                                        .temperature_function = canopy.plant_uptake_growth_temperature_response[plant],
                                    });
                                    const values = try daily_plant_development_bank.row(plant);
                                    values[0] = @floatFromInt(@intFromEnum(output.phase));
                                    values[1] = @floatFromInt(output.branch_count);
                                    values[2] = output.main_branch_stage;
                                    values[3] = output.development_feedback;
                                    values[4] = output.leaf_nitrogen_to_carbon_ratio;
                                    values[5] = output.leaf_phosphorus_to_carbon_ratio;
                                    values[6] = output.minimum_daily_canopy_water_potential_mpa;
                                    values[7] = output.oxygen_stress_factor;
                                    values[8] = output.temperature_function;
                                    const grid_column = runscript.domain.west_column + cell % config.grid_columns;
                                    const grid_row = runscript.domain.north_row + cell / config.grid_columns;
                                    const file_name = try ecosys.output_record.buildPlantFileName(allocator, grid_column, grid_row, species + 1, current_year, development_editor_name);
                                    defer allocator.free(file_name);
                                    _ = try daily_plant_development_bank.streams[plant].write(file_name, daily_plant_development_catalog.variables, selection, resolved.plant_enabled, .{
                                        .timestamp = .{ .year = current_year, .day_of_year = current_day_of_year, .month = management_date.month, .day = management_date.day, .hour = 23 },
                                        .grid_column = grid_column,
                                        .grid_row = grid_row,
                                        .values = values,
                                    });
                                };
                            }
                        }
                    }
                    if (active_options.visualizationIncludesYear(current_year)) {
                        const canopy = if (detailed_canopy_state) |*value|
                            value
                        else
                            return error.VisualizationDailyOutputRequiresCanopyState;
                        const roots = if (plant_root_state) |*value|
                            value
                        else
                            return error.VisualizationDailyOutputRequiresRootState;
                        for (0..state.cell_count) |cell| {
                            var total_leaf_area_m2: f64 = 0;
                            for (0..config.plant_populations) |species| {
                                const plant =
                                    cell * config.plant_populations + species;
                                const branches = try canopy.branchRange(plant);
                                var leaf_carbon_g: f64 = 0;
                                var stalk_carbon_g: f64 = 0;
                                for (branches.first..branches.end) |branch| {
                                    leaf_carbon_g +=
                                        canopy.branch_leaf_carbon_g[branch];
                                    stalk_carbon_g +=
                                        canopy.branch_stalk_carbon_g[branch];
                                    total_leaf_area_m2 +=
                                        canopy.branch_leaf_area_m2[branch];
                                }
                                @memset(
                                    visualization_root_profile_carbon_g,
                                    0,
                                );
                                var root_carbon_g: f64 = 0;
                                for (0..config.soil_layers) |local_layer| {
                                    for (0..root_metabolism_plant_parameters[plant].biologicalDomainCount()) |domain| {
                                        const root = try roots.layerIndex(
                                            plant,
                                            domain,
                                            local_layer,
                                        );
                                        visualization_root_profile_carbon_g[local_layer] +=
                                            roots.mobile_carbon_g[root];
                                        for (0..roots.active_root_axis_count[plant]) |axis| {
                                            const axis_layer =
                                                try roots.layerAxisIndex(
                                                    plant,
                                                    domain,
                                                    local_layer,
                                                    axis,
                                                );
                                            visualization_root_profile_carbon_g[local_layer] +=
                                                roots.axis_primary_carbon_g[axis_layer] +
                                                roots.axis_secondary_carbon_g[axis_layer];
                                        }
                                    }
                                    root_carbon_g +=
                                        visualization_root_profile_carbon_g[local_layer];
                                }
                                visualization_leaf_carbon_g[species] =
                                    leaf_carbon_g;
                                visualization_stalk_carbon_g[species] =
                                    stalk_carbon_g;
                                visualization_root_carbon_g[species] =
                                    root_carbon_g;
                            }

                            const first_layer = try state.layerIndex(cell, 0);
                            const active_layers =
                                state.active_soil_layer_count[cell];
                            const carbon_pools =
                                try ecosys.soil_daily_carbon_pools.calculate(
                                    &soil_organic_state,
                                    &surface_organic_state,
                                    first_layer,
                                    active_layers,
                                    cell,
                                    visualization_organic_carbon_g[1 .. active_layers + 1],
                                );
                            visualization_organic_carbon_g[0] =
                                carbon_pools.surface_noncharcoal_organic_carbon_g_c;
                            @memset(
                                visualization_organic_carbon_g[active_layers + 1 ..],
                                0,
                            );
                            const active_boundary_base =
                                cell * (soil_geometry_state.layer_capacity + 1);
                            const active_surface_depth_m =
                                soil_geometry_state.boundary_depth_m[
                                    active_boundary_base +
                                        soil_geometry_state.first_active_layer[cell]
                                ];
                            const first_plant =
                                cell * config.plant_populations;
                            const plant_end =
                                first_plant + config.plant_populations;
                            const grid_column =
                                runscript.domain.west_column +
                                cell % config.grid_columns;
                            const grid_row =
                                runscript.domain.north_row +
                                cell / config.grid_columns;
                            _ = try visualization_streams.calculateAndWriteDaily(
                                cell,
                                grid_column,
                                grid_row,
                                .{
                                    .year = current_year,
                                    .day_of_year = current_day_of_year,
                                    .month = management_date.month,
                                    .day = management_date.day,
                                    .hour = 23,
                                },
                                .{
                                    .cell_area_m2 = canopy_cell_area_m2[cell],
                                    .cumulative_transpiration_source_m3_by_species = plant_daily_flux_ledger.transpiration_source_m3[first_plant..plant_end],
                                    .cumulative_evapotranspiration_m3 = daily_water_ledger.evaporation_m3[cell],
                                    .cumulative_runoff_m3 = daily_water_ledger.runoff_m3[cell],
                                    .cumulative_outflow_m3 = daily_water_ledger.water_outflow_m3[cell] +
                                        daily_water_ledger.lateral_water_outflow_m3[cell],
                                    .water_table_bottom_depth_m = soil_boundary_topology_state.internal_water_table_depth_m[cell],
                                    .active_surface_depth_m = active_surface_depth_m,
                                    .snow_depth_m = snow_depth_m[cell],
                                    .total_leaf_area_m2 = total_leaf_area_m2,
                                    .leaf_carbon_g_by_species = visualization_leaf_carbon_g,
                                    .stalk_carbon_g_by_species = visualization_stalk_carbon_g,
                                    .root_carbon_g_by_species = visualization_root_carbon_g,
                                    .organic_carbon_g_by_layer_including_surface = visualization_organic_carbon_g,
                                },
                            );
                        }
                    }
                    // Legacy SOIL writes OUTPD after hour 24; DAY carries the
                    // balances and resets these fields before the next hour 1.
                    try plant_daily_flux_ledger.closeDayAfterOutput();
                }
                previous_weather_timestamp = timestamp;
                scene_weather_hours = try std.math.add(usize, scene_weather_hours, 1);
                executed_weather_hours = try std.math.add(usize, executed_weather_hours, 1);
                const scene_final_hour = scene_weather_hours == scene_expected_hours;
                if (try ecosys.checkpoint_schedule.shouldPublish(active_options.checkpoint_output_enabled, active_options.checkpoint_interval_days, current_day_of_year, timestamp.hour, scene_final_hour)) {
                    if (plant_root_state == null or detailed_canopy_state == null or canopy_precipitation_retention_state == null or plant_phenology_state == null or plant_growth_stage_state == null or plant_dormancy_state == null or branch_development_state == null) return error.CheckpointRequestedWithoutCompletePlantState;
                    try ecosys.checkpoint_bundle_writer.publish(
                        allocator,
                        init.io,
                        std.Io.Dir.cwd(),
                        @intCast(executed_weather_hours),
                        .{
                            .year = current_year,
                            .day_of_year = current_day_of_year,
                            .hour = if (timestamp.hour == 24) 23 else timestamp.hour,
                            .execution_iteration = @intCast(pass.execution_iteration),
                            .scenario_index = @intCast(pass.scenario_index),
                            .scenario_iteration = @intCast(pass.scenario_iteration),
                            .scene_index = @intCast(pass.scene_index),
                            .completed_scene_hours = @intCast(scene_weather_hours),
                        },
                        .{ .columns = config.grid_columns, .rows = config.grid_rows, .soil_layers = config.soil_layers, .snow_layers = snow_transport_state.layer_capacity, .plant_species_per_cell = config.plant_populations, .root_axes_per_plant = runscript.root_axes_per_plant },
                        .{
                            .grid = &state,
                            .plants = &plant_state,
                            .plant_metadata_cells = checkpoint_metadata_cells,
                            .plant_development = .{ .phenology = &plant_phenology_state.?, .growth = &plant_growth_stage_state.?, .dormancy = &plant_dormancy_state.?, .branch_development = &branch_development_state.? },
                            .plant_roots = &plant_root_state.?,
                            .plant_canopy = .{ .canopy = &detailed_canopy_state.?, .retention = &canopy_precipitation_retention_state.?, .layer_distribution = &canopy_layer_distribution_state.? },
                            .soil_biogeochemistry = .{ .microbial = &soil_microbial_state, .chemistry = &initial_chemistry_state, .available_nutrients = &plant_available_nutrient_state, .fertilizer = &soil_fertilizer_inventory, .mineral_fertilizer = &mineral_fertilizer_inventory, .fertilizer_band = &fertilizer_band_state, .reactive_nitrogen = &soil_reactive_nitrogen_state, .microbial_phosphorus = &soil_microbial_phosphorus_state },
                            .soil_organic_matter = .{ .profile = &soil_organic_state, .surface = &surface_organic_state, .litter_chemistry = &surface_litter_chemistry_state, .litter_fertilizer = &surface_litter_fertilizer_state, .surface_respiration = &surface_microbial_respiration_state, .surface_denitrification = &surface_denitrification_state, .surface_fire_exchange = &surface_fire_exchange_state, .litter_salt_ingress = &plant_litter_salt_ingress_state },
                            .transport = .{ .micropore = &micropore_solute_state, .macropore = &macropore_solute_state, .mineral_nitrogen = &mineral_nitrogen_transport_state, .organic = &soil_organic_transport_state, .gas = &gas_transport_state, .litter_gas = &litter_gas_transport_state, .snow = &snow_transport_state, .surface = &surface_transport_state },
                            .soil_geometry_and_hydrology = .{ .geometry = &soil_geometry_state, .hydrology = &transport_hydrology_state, .surface = &surface_precipitation_state, .erosion = &surface_erosion_state, .climate = &climate_state, .eroded_minerals = &eroded_mineral_state, .runtime = .{ .soil_properties = &soil_solver_property_state, .soil_thermal = &soil_thermal_state }, .surface_boundary = .{ .ground_air = &ground_air_state, .surface_aerodynamics = &surface_aerodynamic_state }, .surface_litter_geometry = &surface_litter_geometry_state, .surface_litter_ice_m3 = surface_litter_ice_m3, .delayed_live_canopy_combustion_heat_mj = delayed_live_canopy_combustion_heat_mj, .delayed_standing_dead_combustion_heat_mj = delayed_standing_dead_combustion_heat_mj, .delayed_subsurface_combustion_heat_mj = delayed_subsurface_combustion_heat_mj, .delayed_surface_combustion_heat_mj = delayed_surface_combustion_heat_mj },
                            .landscape_mass_balance = landscape_mass_balance_state,
                        },
                        .{ .section_write_buffer_bytes = 64 * 1024, .section_verify_buffer_bytes = 64 * 1024, .manifest_write_buffer_bytes = 16 * 1024 },
                    );
                }
            }
            if (scene_weather_hours == 0) return error.EmptyWeatherStream;
            executed_scene_passes = try std.math.add(usize, executed_scene_passes, 1);
        }
    }
    try hourly_soil_carbon_bank.finish();
    try hourly_plant_carbon_bank.finish();
    try hourly_soil_water_bank.finish();
    try hourly_plant_water_bank.finish();
    try hourly_soil_nitrogen_bank.finish();
    try hourly_plant_nitrogen_bank.finish();
    try hourly_soil_phosphorus_bank.finish();
    try hourly_plant_phosphorus_bank.finish();
    try hourly_soil_heat_bank.finish();
    try daily_soil_carbon_bank.finish();
    try daily_soil_water_bank.finish();
    try daily_soil_phosphorus_bank.finish();
    try daily_soil_nitrogen_bank.finish();
    try daily_soil_heat_bank.finish();
    try hourly_plant_heat_bank.finish();
    try daily_plant_carbon_bank.finish();
    try daily_plant_water_bank.finish();
    try daily_plant_nitrogen_bank.finish();
    try daily_plant_phosphorus_bank.finish();
    try daily_plant_development_bank.finish();
    try visualization_streams.finish();
    if (executed_weather_hours == 0) return error.EmptyWeatherStream;
    try atmospheric_state.validateFinite();
    if (canopy_optics_state) |optics_state| try optics_state.validateFinite();
    if (canopy_structure_state) |structure_state| try structure_state.validateFinite();
    if (canopy_interception_state) |interception_state| try interception_state.validateFinite();
    try ground_radiation_state.validateFinite();
    if (canopy_exposure_state) |exposure_state| try exposure_state.validateFinite();
    try surface_energy_state.validateFinite();
    try surface_temperature_solver_state.validateFinite();
    const surface_solver_diagnostics = try surface_temperature_solver_state.validateConvergence();
    std.log.info("surface Newton-Raphson/Picard diagnostics: cells={d}, iterations={d}, newton_raphson_steps={d}, picard_steps={d}, cells_using_picard={d}, maximum_absolute_residual_mj_per_m2={e}", .{ state.cell_count, surface_solver_diagnostics.total_iterations, surface_solver_diagnostics.total_newton_raphson_steps, surface_solver_diagnostics.total_picard_steps, surface_solver_diagnostics.cells_using_picard, surface_solver_diagnostics.maximum_absolute_residual_mj_per_m2 });
    if (canopy_energy_state) |energy_state| try energy_state.validateFinite();
    if (detailed_canopy_state) |canopy_state| try canopy_state.validateFinite();
    if (plant_root_state) |root_state| try root_state.validateFinite();
    if (branch_development_state) |development_state| try development_state.validateFinite();
    try terrain_radiation_state.validateFinite();
    std.log.info("ecosys-ng inputs validated: cells={d}, soil_layer_capacity={d}, soil_profiles={d}, plant_species_capacity={d}, plant_trait_profiles={d}, plant_management_schedules={d}, plant_assignment_units={d}, land_management_units={d}, disturbance_schedules={d}, fertilizer_schedules={d}, irrigation_schedules={d}, output_selections={d}, landscape_units={d}, scenarios={d}, scenes={d}, planned_days={d}, planned_hours={d}, total_weather_records={d}, executed_scene_passes={d}, executed_weather_hours={d}, latitude_deg_n={d:.3}, first_scene_water_heat_solute_max_iterations={d}, first_scene_gas_max_iterations={d}", .{
        state.cell_count,
        config.soil_layers,
        soil_catalog.entries.items.len,
        config.plant_populations,
        plant_catalog.entries.items.len,
        plant_management_catalog.entries.items.len,
        if (plant_assignments) |assignments| assignments.units.len else 0,
        if (land_management_assignments) |assignments| assignments.units.len else 0,
        disturbance_catalog.entries.items.len,
        fertilizer_catalog.entries.items.len,
        irrigation_catalog.entries.items.len,
        output_selection_catalog.entries.items.len,
        topography.units.len,
        runscript.scenarios.len,
        runscript.scenes.len,
        timeline.execution_weighted_days,
        timeline.execution_weighted_hours,
        total_weather_records,
        executed_scene_passes,
        executed_weather_hours,
        site.latitude_degrees_north,
        iteration_limits.water_heat_solute_max_iterations,
        iteration_limits.gas_max_iterations,
    });
}

fn resolveInputPath(allocator: std.mem.Allocator, io: std.Io, runscript_directory: []const u8, name: []const u8) ![]u8 {
    const direct = try std.fs.path.join(allocator, &.{ runscript_directory, name });
    const direct_file = std.Io.Dir.cwd().openFile(io, direct, .{}) catch |err| switch (err) {
        error.FileNotFound => null,
        else => return err,
    };
    if (direct_file) |file| {
        file.close(io);
        return direct;
    }
    defer allocator.free(direct);
    const parent = std.fs.path.dirname(runscript_directory) orelse return error.InputFileNotFound;
    const inherited = try std.fs.path.join(allocator, &.{ parent, name });
    errdefer allocator.free(inherited);
    const inherited_file = try std.Io.Dir.cwd().openFile(io, inherited, .{});
    inherited_file.close(io);
    return inherited;
}
