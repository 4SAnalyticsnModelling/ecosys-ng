//! Root metabolism and root nutrient uptake for the hourly step.
//!
//! Extracted verbatim from `ecosys_ng.zig` so the entry point holds only
//! `main`. Declaration bodies are unchanged.

const std = @import("std");
const ecosys = @import("ecosys_ng");
pub fn applyRootMetabolism(context: anytype) !void {
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
        const first_domain = try roots.domainIndex(plant, 0);
        const domain_count = context.root_metabolism_plant_parameters[plant].biologicalDomainCount();
        const domain_offsets = [2]usize{ 0, domain_count };
        const oxygen_satisfaction_by_plant = [1]f64{oxygen_satisfaction};
        try ecosys.plant_root_porosity_sweep.apply(
            .{ .current_porosity_fraction_by_domain = roots.current_porosity_fraction_by_domain[first_domain..][0..domain_count] },
            .{
                .plant_count = 1,
                .biological_domain_offsets_by_plant = &domain_offsets,
                .initial_porosity_fraction_by_domain = roots.initial_porosity_fraction_by_domain[first_domain..][0..domain_count],
                .oxygen_satisfaction_fraction_by_plant = &oxygen_satisfaction_by_plant,
                .biological_timestep_h = 1,
                .parameters = context.runscript.root_porosity_parameters,
            },
        );
        context.plant_water_workspace.*.?.root_porosity_fraction[plant] =
            roots.current_porosity_fraction_by_domain[first_domain];
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
                    .canopy_total_water_potential_megapascal = context.plants.canopy_water_potential_megapascal[plant],
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
            // BIND-GROSUB-506, source `grosub.f` 506--512, from
            // `docs/binding_requests/A2_grosub_batch1.md` item A2-B1-03.
            //
            //     IF(PP.GT.ZERO)THEN
            //       WTRTA = AMAX1(0.999992087*WTRTA*XNFH, WTRT/PP)
            //     ELSE
            //       WTRTA = 0.0
            //     ENDIF
            //     XRTN1 = AMAX1(1.0, WTRTA**0.667)*PP
            //
            // THIS FIXES A LIVE SILENT DEFECT IN BOUND CODE, it does not merely
            // bind an idle kernel. `XRTN1` and `RTN1` are DIFFERENT quantities
            // in the source: `XRTN1` is the per-plant axis-count multiplier,
            // while `RTN1` is a per-layer-per-axis accumulator that GROSUB
            // BUILDS FROM `XRTN1` at 6445 and decrements by it at 7243. The
            // consumers of `XRTN1` are 5892, 5927 and 6424. Production was
            // passing `@max(1, roots.axis_primary_count[axis_layer])`, i.e.
            // `RTN1`, into the `XRTN1` slot of
            // `sourceOrderRootAxisSinkStrength`.
            //
            // Worse, `roots.axis_primary_count` has NO PRODUCTION WRITER. A1
            // re-grepped all of `src/`: every assignment
            // (`disturbance_management_dispatch.zig:707`/`:801`,
            // `plant_root_disturbance.zig:1100`, `plant_root_system.zig:1009`,
            // `water_balance.zig:849`) is inside a `test` block, and
            // `reconstructPlant` memsets it to zero. So `@max(1, 0)` pinned the
            // multiplier to the constant `1` for every plant, layer and axis
            // for the whole run, and primary root sink strength lost its entire
            // dependence on root mass per plant and on population. The failure
            // was SILENT because `1` is finite and positive, so no boundary
            // validation could fire.
            //
            // Publish set: `roots.retained_root_carbon_g_c_per_plant[plant]`
            // (`WTRTA`) and the hourly local `primary_axis_count_multiplier`.
            // The `WTRTA` field is NEW in this change and has exactly one
            // writer, this line, so no owner is replaced and nothing is double
            // mutated. It is registered in `plant_root_checkpoint` version 7
            // in the same change, because `WTRTA` is a recurrence on its own
            // previous value and an unserialized recurrence would make a
            // resumed run diverge from a continuous one, which is the Wave 2
            // `RESTART-EQUIVALENCE` obligation A8's `INIT-004` warns about.
            //
            // Ordering, and the one honest caveat. In the source, `WTRT` is
            // accumulated at 13004 (mobile `CPOOLR` plus primary `WTRT1` plus
            // secondary `WTRT2`), i.e. at the END of the plant pass, so line
            // 507 reads the PREVIOUS hour's total. Production has no persisted
            // `WTRT`, so this reassembles it here from the same three pools
            // before this hour's layer loop below mutates any of them. That is
            // the closest faithful analogue available without adding a second
            // persisted total, and A1 states the residual difference rather
            // than hiding it: any earlier step in this hour that moved root
            // carbon is already reflected here, where the source would still be
            // carrying the prior hour's sum.
            const primary_axis_count_multiplier = blk: {
                const population = context.plant_water_workspace.*.?.plant_population_count[plant];
                var total_root_carbon_g_c: f64 = 0;
                for (0..traits.biologicalDomainCount()) |domain| {
                    for (0..context.grid.active_soil_layer_count[cell]) |layer| {
                        const domain_root = try roots.layerIndex(plant, domain, layer);
                        total_root_carbon_g_c += roots.mobile_carbon_g[domain_root];
                        for (0..roots.active_root_axis_count[plant]) |axis| {
                            const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
                            total_root_carbon_g_c += roots.axis_primary_carbon_g[axis_layer] +
                                roots.axis_secondary_carbon_g[axis_layer];
                        }
                    }
                }
                var scaling_state: ecosys.primary_root_axis_scaling.State = .{
                    .retained_root_carbon_g_c_per_plant = roots.retained_root_carbon_g_c_per_plant[plant],
                };
                const scaled = try ecosys.primary_root_axis_scaling.advance(&scaling_state, .{
                    .total_root_carbon_g_c = total_root_carbon_g_c,
                    .plant_population = population,
                    .biological_timestep_h = 1,
                    // `grosub.f` 500--501 names `0.999992087` as the rate of
                    // decline in primary root number, and 512 fixes the
                    // exponent at `0.667`. Both are source literals passed
                    // explicitly rather than promoted to runtime controls, so
                    // this stays a direct translation and no example parameter
                    // file changes.
                    .hourly_retention_fraction = 0.999992087,
                    .axis_scaling_exponent = 0.667,
                });
                // `advance` validates and only then writes its own state, so a
                // failure above leaves `WTRTA` at its previous value and this
                // commit never runs.
                roots.retained_root_carbon_g_c_per_plant[plant] =
                    scaling_state.retained_root_carbon_g_c_per_plant;
                break :blk scaled.primary_root_axis_count_multiplier;
            };
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
                        roots.total_water_potential_megapascal[root],
                        roots.turgor_water_potential_megapascal[root],
                        context.runscript.phenology_parameters.minimum_turgor_potential_megapascal,
                        0,
                        traits.primary_root_radius_m,
                        traits.root_profile_type == 0,
                    );
                    const secondary_environment = try ecosys.plant_root_metabolism.rootEnvironmentResponses(
                        parameters,
                        context.grid.soil_temperature_k[soil],
                        canopy.plant_thermal_adaptation_offset_c[plant],
                        soil_ph,
                        roots.total_water_potential_megapascal[root],
                        roots.turgor_water_potential_megapascal[root],
                        context.runscript.phenology_parameters.minimum_turgor_potential_megapascal,
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
                            // BIND-GROSUB-506: the source's `XRTN1`, computed
                            // once per plant above. This slot previously read
                            // `roots.axis_primary_count`, which is the source's
                            // per-layer `RTN1` and has no production writer, so
                            // it was identically `1`.
                            .primary_axis_count_multiplier = primary_axis_count_multiplier,
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
                            roots.total_water_potential_megapascal[root],
                            roots.turgor_water_potential_megapascal[root],
                            context.runscript.phenology_parameters.minimum_turgor_potential_megapascal,
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

pub fn applyRootNutrientUptake(context: anytype) !void {
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
        const transaction_salt = if (context.runscript.dynamic_plant_salts) try context.plant_root_salt_workspace.*.?.transactionSaltBuffer(cell) else null;
        const transaction_salt_selected = if (context.runscript.dynamic_plant_salts) try context.plant_root_salt_workspace.*.?.transactionSaltSelection(cell) else null;
        var exudation_workspace = &context.plant_root_exudation_workspace.*.?.per_cell[cell];
        const transaction_organic = try context.plant_root_exudation_workspace.*.?.transactionOrganicBuffer(cell);
        const transaction_organic_selected = try context.plant_root_exudation_workspace.*.?.transactionOrganicSelection(cell);
        const plant_first = cell * context.config.plant_populations;
        const active_layer_count = context.grid.active_soil_layer_count[cell];
        const first_soil = try context.grid.layerIndex(cell, 0);
        const layer_thickness_m = context.soil_solver_properties.layer_thickness_m[first_soil..][0..active_layer_count];
        const admission_buffer = try grid_workspace.admissionBuffer(cell);
        var admitted_count: usize = 0;
        for (0..context.config.plant_populations) |species| {
            const plant = plant_first + species;
            var iterator = try ecosys.rooted_layer_admission.Iterator.init(.{
                .plant = plant,
                .active_biological_domain_count = context.root_metabolism_plant_parameters[plant].biologicalDomainCount(),
                .first_active_soil_layer = 0,
                .deepest_rooted_soil_layer = roots.current_deepest_rooted_layer_by_plant[plant],
                .layer_thickness_m = layer_thickness_m,
                .minimum_active_layer_thickness_m = context.runscript.soil_geometry_parameters.minimum_layer_thickness_m,
            });
            while (iterator.next()) |coordinate| {
                if (admitted_count == admission_buffer.len) return error.RootNutrientAdmissionCapacityExceeded;
                admission_buffer[admitted_count] = coordinate;
                admitted_count += 1;
            }
        }
        const admitted = admission_buffer[0..admitted_count];
        if (transaction_salt_selected) |selected| @memset(selected[0..admitted_count], false);
        @memset(transaction_organic_selected[0..admitted_count], false);
        const transaction_nutrient = try grid_workspace.transactionNutrientBuffer(cell);
        const transaction_nutrient_selected = try grid_workspace.transactionNutrientSelection(cell);
        @memset(transaction_nutrient_selected[0..admitted_count], false);
        for (0..context.grid.active_soil_layer_count[cell]) |layer| {
            const soil = try context.grid.layerIndex(cell, layer);
            const soil_pools = try context.plant_available_nutrients.mineralPools(soil);
            var total_previous = [_]f64{0} ** ecosys.plant_root_nutrient_uptake.nutrient_pool_count;
            for (admitted) |coordinate| {
                if (coordinate.soil_layer != layer) continue;
                const root = try roots.layerIndex(coordinate.plant, coordinate.biological_domain, layer);
                for (0..ecosys.plant_root_nutrient_uptake.nutrient_pool_count) |pool_index| {
                    const pool: ecosys.plant_root_nutrient_uptake.NutrientPool = @enumFromInt(pool_index);
                    total_previous[pool_index] += previousRootNutrientDemand(roots, pool, root);
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
            for (admitted, 0..) |coordinate, admission_index| {
                if (coordinate.soil_layer != layer) continue;
                const plant = coordinate.plant;
                const domain = coordinate.biological_domain;
                if (!context.plant_phenology.*.?.active[plant]) continue;
                const root = try roots.layerIndex(plant, domain, layer);
                if (roots.root_surface_area_m2_per_plant[root] <= context.config.absolute_tolerance) continue;
                if (roots.aqueous_volume_m3[root] > context.config.absolute_tolerance) {
                    exudation_workspace.competitors[exudation_competitor_count] = .{ .plant = plant, .domain = domain, .layer = layer };
                    exudation_workspace.admission_index_by_competitor[exudation_competitor_count] = admission_index;
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
                workspace.admission_index_by_competitor[competitor_count] = admission_index;
                competitor_count += 1;
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
                for (0..exudation_competitor_count) |competitor_index| {
                    const result_base = competitor_index * ecosys.plant_root_exudation.substrate_count;
                    const admission_index = exudation_workspace.admission_index_by_competitor[competitor_index];
                    transaction_organic[admission_index] = try ecosys.plant_root_exudation.mapTransactionResult(exudation_workspace.staged_results[result_base..][0..ecosys.plant_root_exudation.substrate_count]);
                    transaction_organic_selected[admission_index] = true;
                }
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
                    var salt_admission_index: ?usize = null;
                    for (admitted, 0..) |coordinate, admission_index| {
                        if (coordinate.plant == competitor.plant and coordinate.biological_domain == competitor.domain and coordinate.soil_layer == competitor.layer) {
                            salt_admission_index = admission_index;
                            break;
                        }
                    }
                    salt_workspace.?.admission_index_by_competitor[salt_competitor_count] = salt_admission_index orelse return error.MissingRootSaltAdmissionCoordinate;
                    salt_competitor_count += 1;
                }
                if (salt_competitor_count > 0) {
                    _ = try salt_workspace.?.stage(
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
                    for (0..salt_competitor_count) |competitor_index| {
                        const result_base = competitor_index * ecosys.plant_root_salt_exchange.species_count;
                        const admission_index = salt_workspace.?.admission_index_by_competitor[competitor_index];
                        transaction_salt.?[admission_index] = try ecosys.plant_root_salt_exchange.mapTransactionResult(salt_workspace.?.staged_exchange_mol[result_base..][0..ecosys.plant_root_salt_exchange.species_count]);
                        transaction_salt_selected.?[admission_index] = true;
                    }
                    try salt_workspace.?.commitStaged(roots, &soil_salt_content_mol, salt_competitor_count);
                }
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
            if (competitor_count > 0) {
                try workspace.stage(roots, soil_pools, competitor_count);
                for (0..competitor_count) |competitor_index| {
                    const result_base = competitor_index * ecosys.plant_root_nutrient_uptake.nutrient_pool_count;
                    const admission_index = workspace.admission_index_by_competitor[competitor_index];
                    transaction_nutrient[admission_index] = try ecosys.plant_root_nutrient_uptake.mapTransactionResult(
                        workspace.staged_results[result_base..][0..ecosys.plant_root_nutrient_uptake.nutrient_pool_count],
                        context.runscript.root_metabolism_parameters.nutrient_uptake_respiration_g_c_per_g_element,
                    );
                    transaction_nutrient_selected[admission_index] = true;
                }
                try workspace.commitStagedAssimilating(
                    roots,
                    soil_pools,
                    competitor_count,
                    context.runscript.root_metabolism_parameters.nutrient_uptake_respiration_g_c_per_g_element,
                );
            }
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
