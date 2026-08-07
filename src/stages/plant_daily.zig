//! Daily plant carbon/nutrient pool accounting, mortality and litterfall.
//!
//! Extracted verbatim from `ecosys_ng.zig` so the entry point holds only
//! `main`. Declaration bodies are unchanged.

const std = @import("std");
const ecosys = @import("ecosys_ng");
pub fn applyPlantStorageRemobilization(context: anytype, calendar: ecosys.plant_development.Calendar) !void {
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
            canopy.plant_canopy_turgor_potential_megapascal[plant],
            context.runscript.phenology_parameters.minimum_turgor_potential_megapascal,
            context.plants.canopy_water_potential_megapascal[plant],
            plant_parameters.stomatal_turgor_shape_per_megapascal,
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

pub fn applyStorageExhaustionMortality(context: anytype, calendar: ecosys.plant_development.Calendar) !void {
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

pub fn applyNaturalBranchMortality(context: anytype, calendar: ecosys.plant_development.Calendar) !void {
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

/// GROSUB lines 12598–12626. Exponential decay of standing dead C/N/P into
/// surface litter at each hourly step. The woody/fine-residue partition is
/// derived from per-plant branch sapwood/total-stalk ratios (FWOOD in
/// Fortran). All three elemental partitions use the same fraction
/// (FWOODN=FWOODP=FWOOD). Output is routed into the surface organic pools via
/// the existing shoot-litter bridge (position-0=woody, position-1=fine residue).
pub fn applyStandingDeadLitterfall(context: anytype) !void {
    if (context.detailed_canopy.* == null) return;
    const canopy = &context.detailed_canopy.*.?;
    const plant_count = canopy.cell_count * canopy.species_count;
    if (plant_count == 0) return;

    const fraction_count = ecosys.standing_dead_litterfall.kinetic_fraction_count;
    const position_count = ecosys.standing_dead_litterfall.litter_position_count;

    var lf_state = try ecosys.standing_dead_litterfall.State.init(context.allocator, plant_count);
    defer lf_state.deinit();

    const partition_count = plant_count * position_count;
    const partitions = try context.allocator.alloc(f64, partition_count);
    defer context.allocator.free(partitions);
    for (0..plant_count) |plant| {
        const b_first = canopy.plant_branch_offsets[plant];
        const b_end = canopy.plant_branch_offsets[plant + 1];
        var total_stalk: f64 = 0;
        var total_sapwood: f64 = 0;
        for (b_first..b_end) |branch| {
            total_stalk += canopy.branch_stalk_carbon_g[branch];
            total_sapwood += canopy.branch_sapwood_carbon_g[branch];
        }
        const sapwood_frac = if (total_stalk > 0)
            std.math.clamp(total_sapwood / total_stalk, 0, 1)
        else
            1.0;
        partitions[plant * position_count + 0] = 1.0 - sapwood_frac;
        partitions[plant * position_count + 1] = sapwood_frac;
    }

    try ecosys.standing_dead_litterfall.apply(
        &lf_state,
        .{
            .biomass_turnover_type_by_plant = context.canopy_layer_controls.biomass_turnover_type,
            .root_profile_type_by_plant = context.canopy_layer_controls.root_profile_type,
            .canopy_growth_temperature_response_by_plant = canopy.plant_uptake_growth_temperature_response,
            .timestep_h = 1,
            .carbon_partition_by_plant_and_position = partitions,
            .nitrogen_partition_by_plant_and_position = partitions,
            .phosphorus_partition_by_plant_and_position = partitions,
        },
        canopy.plant_standing_dead_carbon_by_kinetic_g,
        canopy.plant_standing_dead_nitrogen_by_kinetic_g,
        canopy.plant_standing_dead_phosphorus_by_kinetic_g,
    );

    // Keep per-plant totals in sync with the updated per-fraction pools.
    for (0..plant_count) |plant| {
        var c: f64 = 0;
        var n: f64 = 0;
        var p: f64 = 0;
        const first = plant * fraction_count;
        for (first..first + fraction_count) |i| {
            c += canopy.plant_standing_dead_carbon_by_kinetic_g[i];
            n += canopy.plant_standing_dead_nitrogen_by_kinetic_g[i];
            p += canopy.plant_standing_dead_phosphorus_by_kinetic_g[i];
        }
        canopy.plant_standing_dead_carbon_g[plant] = c;
        canopy.plant_standing_dead_nitrogen_g[plant] = n;
        canopy.plant_standing_dead_phosphorus_g[plant] = p;
    }

    // Commit litter to surface organic (position-0 → woody, position-1 → fine residue).
    for (0..canopy.cell_count) |cell| {
        var products: ecosys.canopy_photosynthesis.SenescenceProducts = .{};
        for (0..canopy.species_count) |species| {
            const plant = cell * canopy.species_count + species;
            for (0..fraction_count) |fraction| {
                const base = (plant * fraction_count + fraction) * position_count;
                products.woody_carbon_g[fraction] += lf_state.carbon_litterfall_g_c_by_plant_fraction_position[base + 0];
                products.woody_nitrogen_g[fraction] += lf_state.nitrogen_litterfall_g_n_by_plant_fraction_position[base + 0];
                products.woody_phosphorus_g[fraction] += lf_state.phosphorus_litterfall_g_p_by_plant_fraction_position[base + 0];
                products.nonwoody_carbon_g[fraction] += lf_state.carbon_litterfall_g_c_by_plant_fraction_position[base + 1];
                products.nonwoody_nitrogen_g[fraction] += lf_state.nitrogen_litterfall_g_n_by_plant_fraction_position[base + 1];
                products.nonwoody_phosphorus_g[fraction] += lf_state.phosphorus_litterfall_g_p_by_plant_fraction_position[base + 1];
            }
        }
        try ecosys.shoot_litter_bridge.commitCell(context.surface_organic, cell, products);
    }
}

pub const DailyPlantElement = enum { nitrogen, phosphorus };

pub fn calculateDailyPlantElementPools(
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

fn calculateDailyPlantCarbonPools(
    canopy: *const ecosys.canopy_photosynthesis.State,
    roots: *const ecosys.plant_root_system.State,
    plant: usize,
    branches: ecosys.canopy_photosynthesis.Range,
    biological_domain_count: usize,
    root_carbon_g_by_layer: []f64,
    root_length_density_m_per_m3_by_layer: []f64,
) !ecosys.plant_daily_pool_aggregation.CarbonPools {
    if (root_carbon_g_by_layer.len != roots.soil_layer_count or
        root_length_density_m_per_m3_by_layer.len != roots.soil_layer_count)
    {
        return error.DailyPlantElementRootDimensionMismatch;
    }
    if (biological_domain_count < 1 or biological_domain_count > ecosys.plant_root_system.biological_domain_count)
        return error.DailyPlantElementRootDimensionMismatch;
    @memset(root_carbon_g_by_layer, 0);
    @memset(root_length_density_m_per_m3_by_layer, 0);
    for (0..roots.soil_layer_count) |layer| {
        for (0..biological_domain_count) |domain| {
            const root = try roots.layerIndex(plant, domain, layer);
            root_carbon_g_by_layer[layer] += roots.mobile_carbon_g[root];
            for (0..roots.active_root_axis_count[plant]) |axis| {
                const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
                root_carbon_g_by_layer[layer] +=
                    roots.axis_primary_carbon_g[axis_layer] +
                    roots.axis_secondary_carbon_g[axis_layer];
            }
        }
    }

    // Preserve source-order aggregation semantics used by plant_daily_output and
    // avoid reusing historical transformed outputs as inputs.
    const primary_root_first = try roots.layerIndex(plant, 0, 0);
    for (0..roots.soil_layer_count) |local_layer| {
        root_length_density_m_per_m3_by_layer[local_layer] =
            roots.root_length_density_m_per_m3[primary_root_first + local_layer];
    }

    var canopy_symbiont_carbon_g: f64 = 0;
    for (branches.first..branches.end) |branch| {
        canopy_symbiont_carbon_g +=
            canopy.branch_symbiont_structural_carbon_g[branch] +
            canopy.branch_symbiont_mobile_carbon_g[branch];
    }
    var root_symbiont_carbon_g: f64 = 0;
    for (0..roots.soil_layer_count) |layer| {
        for (0..biological_domain_count) |domain| {
            const root = try roots.layerIndex(plant, domain, layer);
            root_symbiont_carbon_g +=
                roots.symbiont_structural_carbon_g_c[root] +
                roots.symbiont_mobile_carbon_g_c[root];
        }
    }

    var leaf_intermediate_carbon_g: f64 = 0;
    for (branches.first..branches.end) |branch| {
        const nodes = try canopy.nodeRange(branch);
        for (nodes.first..nodes.end) |node| {
            leaf_intermediate_carbon_g +=
                canopy.node_c3_nonstructural_carbon_g[node] +
                canopy.node_c4_mesophyll_nonstructural_carbon_g[node] +
                canopy.node_bundle_sheath_co2_carbon_g[node] +
                canopy.node_bundle_sheath_bicarbonate_carbon_g[node];
        }
    }

    var projected_leaf_area_m2: f64 = 0;
    for (branches.first..branches.end) |branch| projected_leaf_area_m2 += canopy.branch_leaf_area_m2[branch];

    return ecosys.plant_daily_pool_aggregation.calculateCarbonInto(.{
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
        .root_carbon_g_by_layer = root_carbon_g_by_layer,
        .primary_root_length_density_m_per_m3_by_layer = root_length_density_m_per_m3_by_layer,
        .root_symbiont_carbon_g = root_symbiont_carbon_g,
        .standing_dead_carbon_g = canopy.plant_standing_dead_carbon_g[plant],
        .seed_storage_carbon_g = canopy.plant_seed_storage_carbon_g[plant],
        .projected_leaf_area_m2 = projected_leaf_area_m2,
        .plant_population_count = canopy.plant_population_count[plant],
    }, root_carbon_g_by_layer);
}

pub fn assemblePlantBalanceInputs(
    canopy: *const ecosys.canopy_photosynthesis.State,
    roots: *const ecosys.plant_root_system.State,
    root_metabolism_plant_parameters: []const ecosys.plant_root_metabolism.RuntimePlantParameters,
    plant_daily_flux_ledger: *const ecosys.plant_daily_flux_ledger.State,
    active_by_plant: []const bool,
    carbon_inputs: []ecosys.plant_daily_output.CarbonBalanceInputs,
    nitrogen_inputs: []ecosys.plant_daily_output.NutrientBalanceInputs,
    phosphorus_inputs: []ecosys.plant_daily_output.NutrientBalanceInputs,
    carbon_root_by_layer: []f64,
    carbon_root_length_density_by_layer: []f64,
    nutrient_root_by_layer: []f64,
) !void {
    if (carbon_inputs.len != active_by_plant.len or
        nitrogen_inputs.len != active_by_plant.len or
        phosphorus_inputs.len != active_by_plant.len)
        return error.InvalidPlantBalancePublicationDimensions;

    @memset(std.mem.sliceAsBytes(carbon_inputs), 0);
    @memset(std.mem.sliceAsBytes(nitrogen_inputs), 0);
    @memset(std.mem.sliceAsBytes(phosphorus_inputs), 0);

    for (0..active_by_plant.len) |plant| {
        if (!active_by_plant[plant]) continue;
        const branches = try canopy.branchRange(plant);
        const carbon_pools = try calculateDailyPlantCarbonPools(
            canopy,
            roots,
            plant,
            branches,
            root_metabolism_plant_parameters[plant].biologicalDomainCount(),
            carbon_root_by_layer,
            carbon_root_length_density_by_layer,
        );
        const net_primary_productivity_g =
            plant_daily_flux_ledger.net_carbon_change_g[plant] +
            plant_daily_flux_ledger.signed_total_respiration_carbon_g[plant];
        const carbon_balance_inputs: ecosys.plant_daily_output.CarbonBalanceInputs = .{
            .shoot_carbon_g = carbon_pools.shoot_carbon_g,
            .root_carbon_g = carbon_pools.root_carbon_g,
            .nodule_carbon_g = carbon_pools.nodule_carbon_g,
            .storage_carbon_g = carbon_pools.storage_carbon_g,
            .standing_dead_carbon_g = carbon_pools.vegetative_residue_carbon_g,
            .cumulative_carbon_sink_g = plant_daily_flux_ledger.carbon_sink_g[plant],
            .cumulative_root_soil_carbon_exchange_g = plant_daily_flux_ledger.root_soil_carbon_exchange_g[plant],
            .cumulative_carbon_balance_g = plant_daily_flux_ledger.cumulative_carbon_balance_g[plant],
            .cumulative_harvested_carbon_g = plant_daily_flux_ledger.cumulative_harvested_carbon_g[plant],
            .harvested_carbon_g = plant_daily_flux_ledger.harvested_carbon_g[plant],
            .carbon_oxidation_g = plant_daily_flux_ledger.carbon_oxidation_g[plant],
            .cumulative_net_primary_productivity_g = net_primary_productivity_g,
        };
        carbon_inputs[plant] = carbon_balance_inputs;

        const element_domains = root_metabolism_plant_parameters[plant].biologicalDomainCount();
        const nitrogen_pools = try calculateDailyPlantElementPools(
            canopy,
            roots,
            plant,
            branches,
            .nitrogen,
            element_domains,
            nutrient_root_by_layer,
        );
        const phosphorus_pools = try calculateDailyPlantElementPools(
            canopy,
            roots,
            plant,
            branches,
            .phosphorus,
            element_domains,
            nutrient_root_by_layer,
        );

        nitrogen_inputs[plant] = .{
            .shoot_g = nitrogen_pools.shoot_g,
            .root_g = nitrogen_pools.root_g,
            .nodule_g = nitrogen_pools.nodule_g,
            .storage_g = nitrogen_pools.storage_g,
            .standing_dead_g = nitrogen_pools.vegetative_residue_g,
            .cumulative_sink_g = plant_daily_flux_ledger.nitrogen_sink_g[plant],
            .cumulative_root_soil_exchange_g = plant_daily_flux_ledger.root_soil_nitrogen_exchange_g[plant],
            .cumulative_balance_g = plant_daily_flux_ledger.cumulative_nitrogen_balance_g[plant],
            .cumulative_harvested_g = plant_daily_flux_ledger.cumulative_harvested_nitrogen_g[plant],
            .harvested_g = plant_daily_flux_ledger.harvested_nitrogen_g[plant],
            .oxidation_g = plant_daily_flux_ledger.nitrogen_oxidation_g[plant],
            .atmospheric_exchange_g = plant_daily_flux_ledger.ammonia_exchange_g_n[plant],
            .biological_fixation_g = plant_daily_flux_ledger.symbiotic_nitrogen_fixation_g[plant],
        };
        phosphorus_inputs[plant] = .{
            .shoot_g = phosphorus_pools.shoot_g,
            .root_g = phosphorus_pools.root_g,
            .nodule_g = phosphorus_pools.nodule_g,
            .storage_g = phosphorus_pools.storage_g,
            .standing_dead_g = phosphorus_pools.vegetative_residue_g,
            .cumulative_sink_g = plant_daily_flux_ledger.phosphorus_sink_g[plant],
            .cumulative_root_soil_exchange_g = plant_daily_flux_ledger.root_soil_phosphorus_exchange_g[plant],
            .cumulative_balance_g = plant_daily_flux_ledger.cumulative_phosphorus_balance_g[plant],
            .cumulative_harvested_g = plant_daily_flux_ledger.cumulative_harvested_phosphorus_g[plant],
            .harvested_g = plant_daily_flux_ledger.harvested_phosphorus_g[plant],
            .oxidation_g = plant_daily_flux_ledger.phosphorus_oxidation_g[plant],
        };
    }
}
