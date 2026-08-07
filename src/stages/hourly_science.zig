//! The hourly science step and the atomic state-generation commits.
//!
//! Extracted verbatim from `ecosys_ng.zig` so the entry point holds only
//! `main`. Declaration bodies are unchanged.

const std = @import("std");
const ecosys = @import("ecosys_ng");
const biogeochemistry_batches = @import("biogeochemistry_batches.zig");
const diagnostics = @import("diagnostics.zig");
const plant_daily = @import("plant_daily.zig");
const root_processes = @import("root_processes.zig");
const soil_chemistry_convergence = @import("soil_chemistry_convergence.zig");
const surface_litter_convergence = @import("surface_litter_convergence.zig");
const tile_kernels = @import("tile_kernels.zig");
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
        context.grid.matric_potential_megapascal,
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

fn diagnosticGasNitrogen_g(state: *const ecosys.gas_transport.State) !f64 {
    var total_g_n: f64 = 0;
    inline for (.{ ecosys.gas_transport.Species.nitrogen, .nitrous_oxide, .ammonia }) |species| {
        for (0..state.cell_count) |layer| {
            const index = try ecosys.gas_transport.massIndex(layer, species, state.cell_count);
            total_g_n += state.gaseous_mass_g[index] + state.dissolved_mass_g[index] +
                state.macropore_dissolved_mass_g[index] + state.band_dissolved_mass_g[index];
        }
    }
    if (!std.math.isFinite(total_g_n)) return error.NonFiniteDiagnosticGasNitrogen;
    return total_g_n;
}

fn snowpackInternalNonSaltFluxFromSpeciesAmounts(
    amounts: []const f64,
) ecosys.snowpack_internal_solute_aggregation.SoluteFlux {
    return .{
        .carbon_dioxide_g_c_per_step = amounts[@intFromEnum(ecosys.snow_solute_transport.Species.carbon_dioxide_carbon)],
        .methane_g_c_per_step = amounts[@intFromEnum(ecosys.snow_solute_transport.Species.methane_carbon)],
        .oxygen_g_o_per_step = amounts[@intFromEnum(ecosys.snow_solute_transport.Species.oxygen)],
        .dinitrogen_g_n_per_step = amounts[@intFromEnum(ecosys.snow_solute_transport.Species.dinitrogen_nitrogen)],
        .nitrous_oxide_g_n_per_step = amounts[@intFromEnum(ecosys.snow_solute_transport.Species.nitrous_oxide_nitrogen)],
        .ammonium_g_n_per_step = amounts[@intFromEnum(ecosys.snow_solute_transport.Species.ammonium_nitrogen)],
        .ammonia_g_n_per_step = amounts[@intFromEnum(ecosys.snow_solute_transport.Species.ammonia_nitrogen)],
        .nitrate_g_n_per_step = amounts[@intFromEnum(ecosys.snow_solute_transport.Species.nitrate_nitrogen)],
        .hydrogen_phosphate_g_p_per_step = amounts[@intFromEnum(ecosys.snow_solute_transport.Species.hydrogen_phosphate_phosphorus)],
        .dihydrogen_phosphate_g_p_per_step = amounts[@intFromEnum(ecosys.snow_solute_transport.Species.dihydrogen_phosphate_phosphorus)],
    };
}

pub fn executeHourlyScience(
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
    const diagnostic_first_hour = context.executed_weather_hours.* < 24;
    var diagnostic_previous_n_g = if (diagnostic_first_hour)
        try diagnostics.diagnosticStoredNitrogen_g(context)
    else
        0;
    var diagnostic_previous_p_g = if (diagnostic_first_hour)
        try diagnostics.diagnosticStoredPhosphorus_g(context)
    else
        0;
    var diagnostic_previous_p_owners = if (diagnostic_first_hour)
        try diagnostics.diagnosticPhosphorusOwners_g(context)
    else
        @as([3]f64, @splat(0));
    var diagnostic_previous_heat_megajoules = if (diagnostic_first_hour)
        (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules
    else
        0;
    const diagnostic_relayer_phosphate_before = if (diagnostic_first_hour)
        diagnostics.diagnosticRelayerPhosphateOwners_g(context)
    else
        @as([3]f64, @splat(0));
    context.plant_available_nutrients.resetHourlyChanges();
    try ecosys.soil_organic_carbon_change.captureHourStart(context.soil_organic, context.soil_organic_carbon_at_hour_start_g_c);
    // Derive restart-sensitive surface state from the checkpointed snow owner.
    for (0..context.grid.cell_count) |cell| context.snow_depth_m[cell] = context.snow_transport.cumulative_depth_m[(cell + 1) * context.snow_transport.layer_capacity - 1];
    // Chemistry is concentration-based while TRNSFRS is amount-based. Export
    // before water moves so dilution does not create or destroy solute mass.
    try ecosys.soil_aqueous_transport_bridge.validateCarrierVolumes(
        context.micropore_solute_state,
        context.grid.matrix_liquid_water_m3,
        context.config.absolute_tolerance,
    );
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
    if (diagnostic_first_hour) {
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug("nitrogen stage: mineral_capture delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
        diagnostic_previous_n_g = current_n_g;
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: mineral_capture delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        diagnostic_previous_p_g = current_p_g;
    }
    var diagnostic_mineral_before_mol: f64 = 0;
    const diagnostic_transport_before = if (diagnostic_first_hour) try diagnostics.reconstructLandscapeMassBalance(context) else undefined;

    const diagnostic_transport_ammonium_before = if (diagnostic_first_hour) try diagnostics.diagnosticAmmoniumOwners_g_n(context) else undefined;
    if (diagnostic_first_hour) {
        for (context.mineral_nitrogen_transport.matrix.amount_mol) |amount| diagnostic_mineral_before_mol += amount;
        for (context.mineral_nitrogen_transport.macropore.amount_mol) |amount| diagnostic_mineral_before_mol += amount;
    }
    var forcing_context = ecosys.atmospheric_forcing.MappedApplyContext{ .state = context.atmosphere, .forcing_by_cell = forcing_by_cell };
    try tile_kernels.runKernelAcrossSerialTiles(context, &forcing_context, ecosys.atmospheric_forcing.applyMappedTile);
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
        try tile_kernels.runKernelAcrossSerialTiles(context, &structure_context, ecosys.canopy_structure.applyLeafAreaTile);
    };
    for (weather_header_by_cell, 0..) |header, cell| context.hourly_weather_reference_height_m[cell] = header.aerodynamic_roughness_m;
    var surface_aerodynamic_context: ecosys.surface_aerodynamics.ApplyContext = .{ .state = context.surface_aerodynamics, .cell_area_m2 = context.canopy_cell_area_m2, .total_canopy_area_m2 = context.surface_total_canopy_area_m2, .canopy_height_m = context.surface_canopy_height_m, .snow_depth_m = context.snow_depth_m, .weather_reference_height_m = context.hourly_weather_reference_height_m, .wind_speed_m_per_h = context.atmosphere.wind_speed_m_per_h, .parameters = context.runscript.surface_aerodynamic_parameters };
    try tile_kernels.runKernelAcrossSerialTiles(context, &surface_aerodynamic_context, ecosys.surface_aerodynamics.applyTile);
    try context.ground_air.refreshGeometry(context.canopy_cell_area_m2, context.surface_aerodynamics.wind_reference_height_m, context.runscript.ground_air_parameters);
    for (0..context.grid.cell_count) |cell| context.ground_air_vapor_pressure_kpa[cell] = try ecosys.ground_air_exchange.vaporPressureKpa(context.ground_air.vapor_volume_fraction[cell], context.ground_air.temperature_k[cell], context.runscript.ground_air_parameters);
    for (radiation_by_cell, forcing_by_cell, 0..) |cell_radiation, cell_forcing, cell| {
        context.hourly_extraterrestrial_shortwave_megajoules_per_m2[cell] =
            cell_radiation.extraterrestrial_shortwave_megajoules_per_m2;
        context.hourly_solar_angle_sine[cell] = cell_radiation.solar_angle_sine;
        context.hourly_solar_azimuth_radians[cell] = cell_radiation.solar_azimuth_radians;
        context.hourly_adjusted_shortwave_megajoules_per_m2[cell] =
            cell_forcing.shortwave_radiation_megajoules_per_m2;
    }
    var canopy_radiation_context: ecosys.canopy_radiation.MappedApplyContext = .{ .state = context.canopy_radiation, .horizontal_shortwave_megajoules_per_m2 = context.hourly_adjusted_shortwave_megajoules_per_m2, .extraterrestrial_horizontal_shortwave_megajoules_per_m2 = context.hourly_extraterrestrial_shortwave_megajoules_per_m2, .solar_angle_sine = context.hourly_solar_angle_sine };
    try tile_kernels.runKernelAcrossSerialTiles(context, &canopy_radiation_context, ecosys.canopy_radiation.applyMappedTile);
    if (context.canopy_optics.*) |*optics| {
        var optics_context: ecosys.canopy_optics.ApplyContext = .{ .state = optics, .radiation = context.canopy_radiation };
        try tile_kernels.runKernelAcrossSerialTiles(context, &optics_context, ecosys.canopy_optics.applyLeafAbsorptionTile);
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
        try tile_kernels.runKernelAcrossSerialTiles(context, &interception_context, ecosys.canopy_interception.applySingleLayerTile);
        if (context.canopy_layer_distribution.* != null) try ecosys.canopy_interception.refreshAtmosphericLayerAbsorption(interception, &context.canopy_layer_distribution.*.?, &context.canopy_structure.*.?, &context.canopy_optics.*.?, context.canopy_geometry, context.canopy_radiation, context.direct_incidence_fraction, context.direct_incidence_per_horizontal_area, context.direct_scattering_direction, context.canopy_cell_area_m2, context.runscript.woody_optics_parameters);
    }
    var terrain_radiation_context: ecosys.terrain_radiation.MappedDirectSolarContext = .{ .state = context.terrain_radiation, .solar_angle_sine = context.hourly_solar_angle_sine, .solar_azimuth_radians = context.hourly_solar_azimuth_radians };
    try tile_kernels.runKernelAcrossSerialTiles(context, &terrain_radiation_context, ecosys.terrain_radiation.applyMappedDirectSolarTile);
    var ground_context: ecosys.ground_radiation.ApplyContext = .{ .result = context.ground_radiation, .radiation = context.canopy_radiation, .interception = if (context.canopy_interception.*) |*interception| interception else null, .terrain = context.terrain_radiation };
    try tile_kernels.runKernelAcrossSerialTiles(context, &ground_context, ecosys.ground_radiation.applyTile);
    if (context.canopy_interception.* != null and context.canopy_layer_distribution.* != null) try ecosys.canopy_interception.applyGroundReflectedUpwardSweep(&context.canopy_interception.*.?, &context.canopy_layer_distribution.*.?, &context.canopy_structure.*.?, &context.canopy_optics.*.?, context.canopy_geometry, context.ground_radiation.reflected_shortwave_megajoules_per_m2, context.ground_radiation.reflected_par_micromol_per_m2_per_s, context.runscript.woody_optics_parameters);
    if (context.canopy_precipitation_retention.*) |*retention| try ecosys.canopy_precipitation_retention.refreshFromModel(retention, &context.canopy_layer_distribution.*.?, &context.detailed_canopy.*.?, &context.canopy_interception.*.?, context.atmosphere.rainfall_m, context.canopy_cell_area_m2, context.canopy_layer_controls.root_profile_type, context.hourly_solar_angle_sine, context.ground_radiation.incident_shortwave_megajoules_per_m2, context.runscript.canopy_retention_parameters);
    try ecosys.surface_precipitation.prepareFromModel(context.surface_precipitation, context.atmosphere, context.grid, if (context.canopy_precipitation_retention.*) |*retention| retention else null, context.canopy_cell_area_m2, context.snow_depth_m, context.runscript.snow_full_cover_depth_m, context.iteration_limits.water_heat_solute_max_iterations);
    try ecosys.surface_litter_chemistry_carrier_rebase.rebaseFromAcceptedLiquidWaterChange(
        context.surface_litter_chemistry,
        context.surface_precipitation.litter_water_m3,
        context.surface_precipitation.water_to_litter_m3_per_h,
    );
    try context.snow_transport.commitAtmosphericWater(context.surface_precipitation.snow_to_snow_m3_per_h, context.surface_precipitation.rain_to_snow_m3_per_h, context.surface_precipitation.heat_to_snow_megajoules_per_h, context.atmosphere.air_temperature_k, context.runscript.initial_snow_density_megagrams_per_m3);
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
        .ice_density_megagrams_per_m3 = context.runscript.snow_ice_density_megagrams_per_m3,
        .latent_heat_of_fusion_megajoules_per_m3 = context.runscript.snow_latent_heat_of_fusion_megajoules_per_m3,
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
        {
            const surface_parameters = context.surface_gas_parameters.*;
            const solubility = try ecosys.gas_transport.surfaceSolubilityWaterToAir(context.atmosphere.air_temperature_k[cell], surface_parameters.solubility);
            for (0..5) |species| rain_gas_concentration[species] = context.current_atmospheric_gas_concentration_g_per_m3.*[species] * solubility[species];
        }
        const reaction_parameters = context.chemistry_reaction_parameters.*;
        const nutrients = try ecosys.precipitation_nutrient_speciation.calculate(.{ .ph = weather_header.precipitation_ph, .ammonium_g_n_per_m3 = weather_header.precipitation_ammonium_g_per_m3, .nitrate_g_n_per_m3 = weather_header.precipitation_nitrate_g_per_m3, .phosphate_g_p_per_m3 = weather_header.precipitation_phosphate_g_per_m3 }, reaction_parameters.aqueous_constants, reaction_parameters.phosphate_constants);
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
            irrigation_nutrients = try ecosys.precipitation_nutrient_speciation.calculate(.{
                .ph = irrigation_ph,
                .ammonium_g_n_per_m3 = concentration[0] / irrigation_depth_m,
                .nitrate_g_n_per_m3 = concentration[1] / irrigation_depth_m,
                .phosphate_g_p_per_m3 = concentration[2] / irrigation_depth_m,
            }, reaction_parameters.aqueous_constants, reaction_parameters.phosphate_constants);
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
            const accepted_direct = try ecosys.atmospheric_solute_input_ledger.partitionAcceptedDirectLiquid(
                direct_liquid_m3,
                total_weather_rain_m3 * (1 - snow_fraction),
                irrigation_volume_m3 * (1 - snow_fraction),
                context.config.absolute_tolerance,
            );
            const direct_rain_m3 = accepted_direct.rain_m3;
            const direct_irrigation_m3 = accepted_direct.irrigation_m3;
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
    const diagnostic_snow_n_before_g = if (diagnostic_first_hour)
        diagnostics.diagnosticSnowNitrogen_g(try ecosys.landscape_mass_inventory.aggregateSnow(
            context.snow_transport,
            context.runscript.snow_ice_density_megagrams_per_m3,
        ))
    else
        0;
    const diagnostic_snow_n_input_g = if (diagnostic_first_hour)
        diagnostics.diagnosticSnowSpeciesNitrogen_g(context.snow_atmospheric_input_g)
    else
        0;
    _ = try ecosys.snow_transport_solver.solve(context.allocator, context.snow_transport, .{
        .atmospheric_top_input_g = context.snow_atmospheric_input_g,
        .transport_water_volume_m3 = context.transport_hydrology.snow_liquid_water_volume_m3,
        .water_flux_to_lower_m3 = context.transport_hydrology.snow_downward_water_flux_m3_per_step,
        .litter_water_flux_m3 = context.transport_hydrology.snow_to_litter_water_flux_m3_per_step,
        .soil_micropore_water_flux_m3 = context.transport_hydrology.snow_to_soil_micropore_flux_m3_per_step,
        .soil_macropore_water_flux_m3 = context.transport_hydrology.snow_to_soil_macropore_flux_m3_per_step,
        .surface_partitions = context.snow_surface_partitions,
    }, .{ .absolute_tolerance_g = context.config.absolute_tolerance, .relative_tolerance = context.config.relative_tolerance, .picard_relaxation = context.config.picard_relaxation, .max_iterations = context.iteration_limits.snowpack_max_iterations }, context.snow_surface_discharge);
    if (diagnostic_first_hour) {
        const snow_n_after_g = diagnostics.diagnosticSnowNitrogen_g(
            try ecosys.landscape_mass_inventory.aggregateSnow(
                context.snow_transport,
                context.runscript.snow_ice_density_megagrams_per_m3,
            ),
        );
        const discharge_n_g = diagnostics.diagnosticSnowDischargeNitrogen_g(context.snow_surface_discharge);
        std.log.debug(
            "snow nitrogen transaction: before={e} input={e} after={e} discharge={e} residual={e}",
            .{
                diagnostic_snow_n_before_g,
                diagnostic_snow_n_input_g,
                snow_n_after_g,
                discharge_n_g,
                snow_n_after_g + discharge_n_g - diagnostic_snow_n_before_g - diagnostic_snow_n_input_g,
            },
        );
    }
    try solveSnowSurfaceEnergyAndSoilTransport(
        context,
        hour_of_day,
        weather_header_by_cell,
        plant_calendar_by_cell,
        fertilizer_band_hour,
        gas_failure_report,
        solute_failure_report,
        diagnostic_first_hour,
        diagnostic_mineral_before_mol,
        diagnostic_relayer_phosphate_before,
        diagnostic_transport_ammonium_before,
        diagnostic_transport_before,
        plant_calendar,
        snow_phase_change_report,
        snow_vapor_equilibrium_report,
        &diagnostic_previous_heat_megajoules,
        &diagnostic_previous_n_g,
        &diagnostic_previous_p_g,
        &diagnostic_previous_p_owners,
    );
}

fn finalizeVegetationAndLedgers(
    context: anytype,
    hour_of_day: u8,
    weather_header_by_cell: []const ecosys.weather.Header,
    plant_calendar_by_cell: []const ecosys.plant_development.Calendar,
    diagnostic_first_hour: anytype,
    diagnostic_previous_p_g: anytype,
    plant_calendar: anytype,
    snow_phase_change_report: anytype,
    snow_vapor_equilibrium_report: anytype,
    subsurface_irrigation_chemistry_parameters: anytype,
) !void {
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
        try tile_kernels.runKernelAcrossSerialTiles(context, &canopy_energy_context, ecosys.canopy_energy.applyTile);
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
        try tile_kernels.runKernelAcrossSerialTiles(context, &biochemistry_context, ecosys.canopy_biochemistry.applyTile);
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
            .soil_total_water_potential_megapascal = context.soil_hourly_workspace.root_referenced_total_water_potential_megapascal,
            .active = workspace.active,
            .root_conductance_m_per_h_megapascal = workspace.root_conductance_m_per_h_megapascal,
            .maximum_uptake_m = workspace.maximum_uptake_m,
            .maximum_release_m = workspace.maximum_release_m,
            .canopy_water_capacitance_m_per_m2_megapascal = workspace.canopy_water_capacitance_m_per_m2_megapascal,
            .transpiration_loss_m = workspace.transpiration_loss_m,
            .settings = .{ .minimum_canopy_water_potential_megapascal = -100.0, .maximum_canopy_water_potential_megapascal = 0.0, .solver_options = context.nonlinear_solver_options },
        };
        try tile_kernels.runKernelAcrossSerialTiles(context, &water_context, ecosys.plant_water_balance.applyTile);
        try ecosys.plant_water_balance.commitRootHydraulics(balance.*, &context.plant_roots.*.?, context.grid, context.soil_hourly_workspace.root_referenced_total_water_potential_megapascal, context.plants, workspace.active, workspace.cell_area_m2, workspace.soil_resistance_mpa_h_per_m, workspace.root_resistance_mpa_h_per_m, workspace.leaf_osmotic_potential_at_zero_total_megapascal, context.root_biological_domain_count_by_plant);
        try ecosys.plant_root_water_storage_commit.commit(
            &context.plant_roots.*.?,
            context.grid,
            context.soil_chemistry,
            context.root_biological_domain_count_by_plant,
        );
        try ecosys.plant_water_balance.updateDailyMinimumCanopyWaterPotential(&context.detailed_canopy.*.?, context.plants);
        try ecosys.plant_root_gas_exchange.refreshOxygenDemand(&context.plant_roots.*.?, context.root_gas_parameters);
    }
    try plant_daily.applyPlantStorageRemobilization(context, plant_calendar);
    try root_processes.applyRootMetabolism(context);
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
    try plant_daily.applyStorageExhaustionMortality(context, plant_calendar);
    try root_processes.applyRootNutrientUptake(context);
    if (context.plant_growth_stages.*) |*growth_stages| {
        const roots = if (context.plant_roots.*) |*value| value else return error.MissingPlantRootsForPoolAggregation;
        const canopy = if (context.detailed_canopy.*) |*value| value else return error.MissingCanopyForPoolAggregation;
        var pool_context: ecosys.plant_pool_aggregation.ApplyContext = .{ .canopy = canopy, .roots = roots, .growth_stages = growth_stages, .active_by_plant = context.plant_phenology.*.?.active, .biological_domain_count_by_plant = context.root_biological_domain_count_by_plant, .dynamic_salts = context.runscript.dynamic_plant_salts, .parameters = context.runscript.plant_pool_parameters };
        try tile_kernels.runKernelAcrossSerialTiles(context, &pool_context, ecosys.plant_pool_aggregation.applyTile);
    }
    if (context.plant_growth_stages.* != null and context.plant_roots.* != null and context.detailed_canopy.* != null and context.branch_development.* != null and context.plant_phenology.* != null) {
        const topology_states: ecosys.plant_topology.RuntimeStates = .{ .canopy = &context.detailed_canopy.*.?, .growth_stages = &context.plant_growth_stages.*.?, .dormancy = &context.plant_dormancy.*.?, .branch_development = &context.branch_development.*.? };
        const plant_execution_year = std.math.cast(u16, plant_calendar.current_year) orelse
            return error.InvalidPlantTopologyDate;
        const new_shoot_branches = try ecosys.plant_topology.advanceShootBranches(.{ .states = topology_states, .roots = &context.plant_roots.*.?, .controls = context.plant_topology_controls, .active_by_plant = context.plant_phenology.*.?.active, .emerged_by_plant = context.development_emerged, .day_of_year = plant_calendar.day_of_year, .execution_year = plant_execution_year, .minimum_root_turgor_potential_megapascal = context.runscript.phenology_parameters.minimum_turgor_potential_megapascal });
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
            .minimum_root_turgor_potential_megapascal = context.runscript.phenology_parameters.minimum_turgor_potential_megapascal,
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
            .canopy_turgor_potential_megapascal = canopy.plant_canopy_turgor_potential_megapascal,
            .root_oxygen_uptake_to_demand_fraction = context.phenology_root_oxygen_fraction,
            .emerged = context.development_emerged,
            .timestep_hours = 1,
            .parameters = context.runscript.phenology_parameters,
        };
        try tile_kernels.runKernelAcrossSerialTiles(context, &phenology_context, ecosys.plant_phenology.advanceTile);
    }
    if (context.plant_growth_stages.*) |*growth_stages| {
        const roots = if (context.plant_roots.*) |*value| value else return error.MissingPlantRootsForDevelopment;
        const canopy = if (context.detailed_canopy.*) |*value| value else return error.MissingCanopyForDevelopment;
        const dormancy_state = if (context.plant_dormancy.*) |*value| value else return error.MissingDormancyStateForDevelopment;
        try ecosys.plant_development.refreshCanopyHeight(canopy, context.development_canopy_height_m);
        try ecosys.plant_development.refreshSoilWaterPotentials(
            context.grid,
            context.soil_hourly_workspace.landscape_total_water_potential_megapascal,
            context.surface_litter_water_environment.matric_plus_osmotic_water_potential_megapascal,
            context.terrain_hydrology.relative_surface_elevation_m,
            context.soil_hourly_workspace.gravitational_water_potential_mpa_per_m,
            roots,
            context.development_surface_water_potential_megapascal,
            context.development_seed_layer_water_potential_megapascal,
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
            .canopy_turgor_potential_mpa_by_plant = canopy.plant_canopy_turgor_potential_megapascal,
            .canopy_total_water_potential_mpa_by_plant = context.plants.canopy_water_potential_megapascal,
            .surface_soil_water_potential_mpa_by_cell = context.development_surface_water_potential_megapascal,
            .seed_layer_soil_water_potential_mpa_by_plant = context.development_seed_layer_water_potential_megapascal,
            .emerged_by_plant = context.development_emerged,
            .calendar_by_cell = plant_calendar_by_cell,
            .timestep_h = 1,
        };
        try tile_kernels.runKernelAcrossSerialTiles(context, &development_context, ecosys.plant_development.advanceTile);
        var reproduction_context: ecosys.plant_reproduction.ApplyContext = .{ .canopy = canopy, .plants = context.plants, .growth_stages = growth_stages, .controls = context.plant_reproduction_controls, .active_by_plant = context.plant_phenology.*.?.active, .minimum_turgor_potential_megapascal = context.runscript.phenology_parameters.minimum_turgor_potential_megapascal, .seed_set_parameters = context.runscript.seed_set_parameters, .structural_presence_threshold_g_per_plant = context.runscript.plant_pool_parameters.branch_structural_presence_g_per_plant, .timestep_h = 1 };
        try tile_kernels.runKernelAcrossSerialTiles(context, &reproduction_context, ecosys.plant_reproduction.applyTile);
    }
    // GROSUB 4182--4378: spring phenology and residue cleanup precede the
    // 4393--4409 seasonal flag reset executed by shoot_growth_runtime.
    if (context.plant_harvest) |harvest| if (context.plant_phenology.*) |*phenology_state| {
        for (phenology_state.leafout_transition_this_step, 0..) |leafout, plant| {
            if (!leafout or !phenology_state.active[plant] or context.plant_topology_controls.growth_habit_code[plant] == 0) continue;
            const branches = try context.plant_growth_stages.*.?.branchRange(plant);
            const fully_deciduous = context.canopy_layer_controls.biomass_turnover_type[plant] == 0;
            for (branches.first..branches.end) |branch| {
                context.branch_development.*.?.remobilization_progress_h[branch] = 0;
                context.branch_development.*.?.leafout_initialization_enabled[branch] = false;
                try ecosys.plant_growth_stages.resetBranchForSeasonalLeafout(
                    &context.plant_growth_stages.*.?.branches[branch],
                    plant_calendar.day_of_year,
                    fully_deciduous,
                    context.branch_development.*.?.initial_reproductive_stage[branch],
                );
            }
            if (fully_deciduous and branches.first < branches.end) {
                phenology_state.initiated_node_count[plant] = context.branch_development.*.?.initial_reproductive_stage[branches.first];
                phenology_state.appeared_leaf_count[plant] = 0;
            }
            try ecosys.plant_harvest_runtime.applyStartOfSeasonResidue(
                harvest,
                plant,
                context.canopy_layer_controls.biomass_turnover_type[plant],
                context.canopy_layer_controls.root_profile_type[plant],
            );
        }
    };
    if (context.canopy_layer_distribution.*) |*layers| {
        const canopy = if (context.detailed_canopy.*) |*value| value else return error.MissingCanopyForLayerDistribution;
        const growth_stages = if (context.plant_growth_stages.*) |*value| value else return error.MissingGrowthStagesForLayerDistribution;
        const structure = if (context.canopy_structure.*) |*value| value else return error.MissingCanopyStructureForLayerDistribution;
        const water_workspace = if (context.plant_water_workspace.*) |*value| value else return error.MissingPlantWaterWorkspaceForLayerDistribution;
        try layers.refresh(canopy, growth_stages, context.canopy_layer_controls, water_workspace.seeding_depth_m, context.canopy_geometry.leaf_inclination_sine, structure.leaf_inclination_fraction, context.hourly_solar_angle_sine, 1.0e-12);
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
                .canopy_turgor_potential_mpa_by_plant = canopy.plant_canopy_turgor_potential_megapascal,
                .leaf_carbon_presence_threshold_g_c_per_plant = context.runscript.plant_pool_parameters.branch_structural_presence_g_per_plant,
                .minimum_turgor_potential_megapascal = context.runscript.phenology_parameters.minimum_turgor_potential_megapascal,
                .atmospheric_co2_umol_per_mol = context.current_atmospheric_co2_umol_per_mol.*,
                .picard_relaxation = context.config.picard_relaxation,
                .timestep_h = 1,
                .carbon_exchange = if (context.canopy_carbon_exchange.*) |*ledger| ledger else null,
            };
            try tile_kernels.runKernelAcrossSerialTiles(context, &carboxylation_context, ecosys.canopy_carboxylation.applyTile);
        }
    }
    if (context.detailed_canopy.*) |*canopy| {
        const shoot_execution_year = std.math.cast(u16, plant_calendar.current_year) orelse
            return error.InvalidShootGrowthDate;
        const litter_partition = if (context.plant_litter_partition.*) |*value| value else return error.MissingPlantLitterPartitionForShootGrowth;
        const canopy_surface_exchange_for_growth = if (context.canopy_surface_exchange.*) |*value| value else return error.MissingCanopySurfaceExchangeForShootGrowth;
        const canopy_surface_workspace_for_growth = if (context.canopy_surface_input_workspace.*) |*value| value else return error.MissingCanopySurfaceWorkspaceForShootGrowth;
        const canopy_retention_for_growth = if (context.canopy_precipitation_retention.*) |*value| value else return error.MissingCanopyRetentionForShootGrowth;
        const surface_gas_for_growth = context.surface_gas_parameters.*;
        const shoot_solar_noon_hour_by_cell = try context.allocator.alloc(u8, context.grid.cell_count);
        defer context.allocator.free(shoot_solar_noon_hour_by_cell);
        for (weather_header_by_cell, shoot_solar_noon_hour_by_cell) |header, *solar_noon_hour| {
            if (!std.math.isFinite(header.solar_noon_hour) or header.solar_noon_hour < 0 or header.solar_noon_hour > 23)
                return error.InvalidShootSolarNoonHour;
            solar_noon_hour.* = @intFromFloat(@floor(header.solar_noon_hour));
        }
        @memset(context.shoot_senescence_products_by_plant, .{});
        @memset(context.seasonal_turnover_event_by_plant, false);
        var shoot_growth_context: ecosys.shoot_growth_runtime.ApplyContext = .{
            .canopy = canopy,
            .growth_stages = &context.plant_growth_stages.*.?,
            .dormancy = &context.plant_dormancy.*.?,
            .development = &context.branch_development.*.?,
            .plant_parameters = context.shoot_growth_plant_parameters,
            .active_by_plant = context.plant_phenology.*.?.active,
            .emerged_by_plant = context.development_emerged,
            .sowing_depth_m_by_plant = context.plant_water_workspace.*.?.seeding_depth_m,
            .roots = &context.plant_roots.*.?,
            .soil_temperature_k = context.grid.soil_temperature_k,
            .soil_layer_capacity = context.grid.soil_layer_capacity,
            .root_growth_temperature_parameters = context.runscript.canopy_stress_parameters.growth_temperature,
            .canopy_temperature_k_by_plant = context.plants.canopy_temperature_k,
            .canopy_total_water_potential_mpa_by_plant = context.plants.canopy_water_potential_megapascal,
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
            .seasonal_litterfall_rate_per_h = context.runscript.seasonal_turnover_parameters.litterfall_rate_per_h,
            .seasonal_litterfall_delay_threshold_h = context.runscript.seasonal_turnover_parameters.litterfall_delay_threshold_h,
            .dormancy_parameters_by_plant = context.development_dormancy_parameters,
            .litter_partition = litter_partition,
            .senescence_recycling = .{
                .minimum_carbon_fraction = context.runscript.root_metabolism_parameters.minimum_carbon_recycling_fraction,
                .responsive_carbon_fraction = context.runscript.root_metabolism_parameters.responsive_carbon_recycling_fraction,
                .maximum_nitrogen_fraction = context.runscript.root_metabolism_parameters.maximum_nitrogen_recycling_fraction,
                .maximum_phosphorus_fraction = context.runscript.root_metabolism_parameters.maximum_phosphorus_recycling_fraction,
            },
            .senescence_products_by_plant = context.shoot_senescence_products_by_plant,
            .seasonal_turnover_event_by_plant = context.seasonal_turnover_event_by_plant,
            .senescence_demand_tolerance_g_c = context.config.absolute_tolerance,
            .leaf_area_presence_tolerance_m2 = context.config.absolute_tolerance,
            .symbiosis_parameters = context.runscript.symbiotic_fixation_parameters,
            .carbon_exchange = if (context.canopy_carbon_exchange.*) |*ledger| ledger else null,
        };
        try tile_kernels.runKernelAcrossSerialTiles(context, &shoot_growth_context, ecosys.shoot_growth_runtime.applyTile);
        try plant_daily.applyNaturalBranchMortality(context, plant_calendar);
        for (0..canopy.cell_count) |cell| {
            var cell_products: ecosys.canopy_photosynthesis.SenescenceProducts = .{};
            const first_plant = cell * canopy.species_count;
            for (context.shoot_senescence_products_by_plant[first_plant .. first_plant + canopy.species_count]) |products|
                ecosys.canopy_photosynthesis.addSenescenceProducts(&cell_products, products);
            try ecosys.shoot_litter_bridge.commitCell(context.surface_organic, cell, cell_products);
        }
        try plant_daily.applyStandingDeadLitterfall(context);
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
        try tile_kernels.runKernelAcrossSerialTiles(context, &root_mycorrhizal_exchange_context, ecosys.plant_root_mycorrhizal_exchange.applyTile);
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
        try tile_kernels.runKernelAcrossSerialTiles(context, &shoot_root_exchange_context, ecosys.plant_shoot_root_exchange.applyTile);
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
        snow_phase_change_report.sensible_energy_change_megajoules +
            snow_vapor_equilibrium_report.sensible_energy_change_megajoules,
    );
    const subsurface_irrigation_input =
        try ecosys.subsurface_irrigation_heat.calculate(
            context.subsurface_irrigation_water_m3,
            context.atmosphere.air_temperature_k,
            context.grid.soil_layer_capacity,
            context.runscript.canopy_surface_exchange_parameters.liquid_water_heat_capacity_megajoules_per_m3_k,
        );
    try context.landscape_boundary_ledger.accumulateAccepted(.{
        .rain_m3 = subsurface_irrigation_input.water_input_m3,
        .heat_input_megajoules = subsurface_irrigation_input.heat_input_megajoules,
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
    if (diagnostic_first_hour) {
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: vegetation_and_final_ledgers delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
    }
}

fn routeSedimentAndErosion(
    context: anytype,
    hour_of_day: u8,
    weather_header_by_cell: []const ecosys.weather.Header,
    plant_calendar_by_cell: []const ecosys.plant_development.Calendar,
    fertilizer_band_hour: ecosys.fertilizer_band_phase_coordinator.HourToken,
    solute_failure_report: ?ecosys.solute_failure_reporter.Request,
    diagnostic_first_hour: anytype,
    plant_calendar: anytype,
    snow_phase_change_report: anytype,
    snow_vapor_equilibrium_report: anytype,
    subsurface_irrigation_chemistry_parameters: anytype,
    diagnostic_previous_heat_megajoules_ptr: anytype,
    diagnostic_previous_n_g_ptr: anytype,
    diagnostic_previous_p_g_ptr: anytype,
) !void {
    var diagnostic_previous_heat_megajoules = diagnostic_previous_heat_megajoules_ptr.*;
    defer diagnostic_previous_heat_megajoules_ptr.* = diagnostic_previous_heat_megajoules;
    var diagnostic_previous_n_g = diagnostic_previous_n_g_ptr.*;
    defer diagnostic_previous_n_g_ptr.* = diagnostic_previous_n_g;
    var diagnostic_previous_p_g = diagnostic_previous_p_g_ptr.*;
    defer diagnostic_previous_p_g_ptr.* = diagnostic_previous_p_g;
    try ecosys.sediment_routing.route(
        &context.surface_erosion.routing,
        context.surface_erosion.transportable_sediment_megagrams,
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
    try ecosys.sediment_routing.commitSurfaceSediment(context.surface_erosion.surface_sediment_megagrams, context.surface_erosion.local_detachment_megagrams, context.surface_erosion.routing.sediment_change_megagrams);
    try ecosys.soil_erosion_organic_bridge.route(
        context.config.lon_count,
        context.config.lat_count,
        context.grid.soil_layer_capacity,
        context.surface_erosion.surface_soil_mass_megagrams,
        context.soil_organic,
        .{
            .east_megagrams = context.surface_erosion.routing.east_flux_megagrams,
            .west_megagrams = context.surface_erosion.routing.west_flux_megagrams,
            .south_megagrams = context.surface_erosion.routing.south_flux_megagrams,
            .north_megagrams = context.surface_erosion.routing.north_flux_megagrams,
        },
        context.eroded_organic_workspace,
    );
    try ecosys.soil_erosion_organic_bridge.refreshSurfaceOrganicCarbonGPerMg(
        context.soil_organic,
        context.grid.soil_layer_capacity,
        context.surface_erosion.surface_soil_mass_megagrams,
        context.soil_solver_properties.total_organic_carbon_g_per_megagram,
    );
    try ecosys.soil_erosion_fertilizer_bridge.route(
        context.config.lon_count,
        context.config.lat_count,
        context.surface_erosion.surface_soil_mass_megagrams,
        context.soil_fertilizer_inventory,
        .{
            .east_megagrams = context.surface_erosion.routing.east_flux_megagrams,
            .west_megagrams = context.surface_erosion.routing.west_flux_megagrams,
            .south_megagrams = context.surface_erosion.routing.south_flux_megagrams,
            .north_megagrams = context.surface_erosion.routing.north_flux_megagrams,
        },
        context.eroded_fertilizer_workspace,
    );
    try ecosys.soil_erosion_chemistry_bridge.route(
        context.config.lon_count,
        context.config.lat_count,
        context.grid.soil_layer_capacity,
        context.surface_erosion.surface_soil_mass_megagrams,
        context.grid.matrix_liquid_water_m3,
        context.soil_chemistry,
        .{
            .east_megagrams = context.surface_erosion.routing.east_flux_megagrams,
            .west_megagrams = context.surface_erosion.routing.west_flux_megagrams,
            .south_megagrams = context.surface_erosion.routing.south_flux_megagrams,
            .north_megagrams = context.surface_erosion.routing.north_flux_megagrams,
        },
        context.eroded_chemistry_workspace,
    );
    try ecosys.soil_erosion_mineral_bridge.route(
        context.config.lon_count,
        context.config.lat_count,
        context.grid.soil_layer_capacity,
        context.surface_erosion.surface_soil_mass_megagrams,
        context.soil_solver_properties,
        .{
            .east_megagrams = context.surface_erosion.routing.east_flux_megagrams,
            .west_megagrams = context.surface_erosion.routing.west_flux_megagrams,
            .south_megagrams = context.surface_erosion.routing.south_flux_megagrams,
            .north_megagrams = context.surface_erosion.routing.north_flux_megagrams,
        },
        context.eroded_mineral_state,
    );
    // Mineral routing is the authoritative commit of transported surface
    // soil mass. Rebind the organic concentration to that accepted carrier
    // immediately; the extensive organic pools were routed above.
    try ecosys.soil_erosion_organic_bridge.refreshSurfaceOrganicCarbonGPerMg(
        context.soil_organic,
        context.grid.soil_layer_capacity,
        context.surface_erosion.surface_soil_mass_megagrams,
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
    try ecosys.soil_sediment_change.publishAcceptedNetSedimentMg(context.surface_soil_mass_at_erosion_start_megagrams, context.surface_erosion.surface_soil_mass_megagrams, context.net_sediment_megagrams_per_h);
    @memcpy(context.transport_hydrology.runoff_total_m3_per_step, context.surface_runoff.total_runoff_m3_per_step);
    @memcpy(context.transport_hydrology.runoff_east_m3_per_step, context.surface_runoff.east_runoff_m3_per_step);
    @memcpy(context.transport_hydrology.runoff_west_m3_per_step, context.surface_runoff.west_runoff_m3_per_step);
    @memcpy(context.transport_hydrology.runoff_south_m3_per_step, context.surface_runoff.south_runoff_m3_per_step);
    @memcpy(context.transport_hydrology.runoff_north_m3_per_step, context.surface_runoff.north_runoff_m3_per_step);
    if (diagnostic_first_hour) {
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug("nitrogen stage: erosion delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
        diagnostic_previous_n_g = current_n_g;
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: erosion delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        diagnostic_previous_p_g = current_p_g;
    }
    for (context.snow_surface_discharge, context.direct_surface_solute_input) |*snow_discharge, direct_input| {
        for (&snow_discharge.litter_g, direct_input.litter_g) |*destination, amount| destination.* += amount;
        for (&snow_discharge.soil_nonband_g, direct_input.soil_nonband_g) |*destination, amount| destination.* += amount;
        for (&snow_discharge.soil_band_g, direct_input.soil_band_g) |*destination, amount| destination.* += amount;
    }
    const ion_parameters = context.runscript.chemistry_primary_initialization.molar_mass_g_per_mol;
    const ion_molar_masses: ecosys.snow_surface_discharge.IonMolarMassesGPerMol = .{ .aluminum = ion_parameters.aluminum, .iron = ion_parameters.iron, .calcium = ion_parameters.calcium, .magnesium = ion_parameters.magnesium, .sodium = ion_parameters.sodium, .potassium = ion_parameters.potassium, .sulfur = ion_parameters.sulfur, .chloride = ion_parameters.chloride };
    try ecosys.snow_surface_discharge.commit(context.allocator, .{ .discharge = context.snow_surface_discharge, .litter_water_volume_m3 = context.surface_precipitation.litter_water_m3, .topsoil_water_volume_m3 = context.grid.matrix_liquid_water_m3[0..context.grid.layer_count], .soil_layer_capacity = context.grid.soil_layer_capacity, .ion_molar_mass_g_per_mol = ion_molar_masses }, context.litter_gas_transport, context.gas_transport, context.surface_litter_chemistry, context.soil_chemistry);
    if (diagnostic_first_hour) {
        const heat_megajoules = (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules;
        std.log.debug("heat stage: snow_discharge hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, heat_megajoules - diagnostic_previous_heat_megajoules });
        diagnostic_previous_heat_megajoules = heat_megajoules;
    }
    if (diagnostic_first_hour) {
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug("nitrogen stage: snow_discharge delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
        diagnostic_previous_n_g = current_n_g;
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: snow_discharge delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        diagnostic_previous_p_g = current_p_g;
    }
    // Earlier REDIST/transport steps mutate the fixed organic microbial
    // mirror. Refresh the runtime NITRO owner before it computes metabolism;
    // otherwise the post-NITRO mirror publication overwrites transported CNP.
    try ecosys.soil_microbial_inventory_bridge.publishFromOrganic(
        context.soil_organic,
        context.soil_microbial,
    );
    const diagnostic_biogeochemistry_n_before_g = if (diagnostic_first_hour)
        try diagnostics.diagnosticStoredNitrogen_g(context)
    else
        0;
    try biogeochemistry_batches.runSoilBiogeochemistryBySerialTile(context);
    if (diagnostic_first_hour) {
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug("nitrogen stage: soil_biogeochemistry delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
        diagnostic_previous_n_g = current_n_g;
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: soil_biogeochemistry delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        diagnostic_previous_p_g = current_p_g;
    }
    {
        const surface_parameters = context.surface_gas_parameters.*;
        const old_litter_dry_mass_megagrams = try context.allocator.dupe(
            f64,
            context.surface_litter_geometry.dry_mass_megagrams,
        );
        defer context.allocator.free(old_litter_dry_mass_megagrams);
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
        try tile_kernels.runKernelAcrossSerialTiles(context, &litter_geometry_context, ecosys.surface_litter_geometry_step.applyTile);
        try ecosys.surface_litter_chemistry_carrier_rebase.rebaseFromAcceptedDryMassChange(
            context.surface_litter_chemistry,
            old_litter_dry_mass_megagrams,
            context.surface_litter_geometry.dry_mass_megagrams,
        );
        @memcpy(context.surface_precipitation.litter_water_capacity_m3, context.surface_litter_geometry.water_retention_capacity_m3);
        for (0..context.grid.cell_count) |cell| {
            context.litter_gas_transport.air_volume_m3[cell] = context.surface_litter_geometry.air_volume_m3[cell];
            context.litter_gas_transport.temperature_k[cell] = context.grid.surface_temperature_k[cell];
        }
        if (diagnostic_first_hour) {
            const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
            std.log.debug("phosphorus stage: surface_litter_geometry delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
            diagnostic_previous_p_g = current_p_g;
        }
        const diagnostic_litter_gas_n_before_g = if (diagnostic_first_hour)
            try diagnosticGasNitrogen_g(context.litter_gas_transport)
        else
            0;
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
        if (diagnostic_first_hour) {
            var diagnostic_litter_boundary_n_g: f64 = 0;
            for (0..context.grid.cell_count) |cell| {
                const base = cell * ecosys.gas_transport.species_count;
                inline for ([_]ecosys.gas_transport.Species{ .nitrogen, .nitrous_oxide, .ammonia }) |species|
                    diagnostic_litter_boundary_n_g += context.surface_litter_gas_transport.atmospheric_flux_g_per_h[base + @intFromEnum(species)];
            }
            const diagnostic_litter_after_g = try diagnosticGasNitrogen_g(context.litter_gas_transport);
            std.log.debug("litter gas nitrogen transaction: delta_g={e} boundary_g={e} residual_g={e}", .{ diagnostic_litter_after_g - diagnostic_litter_gas_n_before_g, diagnostic_litter_boundary_n_g, diagnostic_litter_after_g - diagnostic_litter_gas_n_before_g - diagnostic_litter_boundary_n_g });
        }
        if (diagnostic_first_hour) {
            const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
            std.log.debug("nitrogen stage: litter_gas_transport delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
            diagnostic_previous_n_g = current_n_g;
            const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
            std.log.debug("phosphorus stage: litter_gas_transport delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
            diagnostic_previous_p_g = current_p_g;
            const heat_megajoules = (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules;
            std.log.debug("heat stage: litter_gas_transport hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, heat_megajoules - diagnostic_previous_heat_megajoules });
            diagnostic_previous_heat_megajoules = heat_megajoules;
        }
        @memset(context.surface_microbial_substrate_uptake.denitrification_respiration_g_c, 0);
        context.surface_litter_fertilizer_diagnostics.reset();
        try biogeochemistry_batches.runSurfaceBiogeochemistryBySerialTile(context, surface_parameters);
        if (diagnostic_first_hour) {
            const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
            std.log.debug("nitrogen stage: surface_biogeochemistry delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
            diagnostic_previous_n_g = current_n_g;
            const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
            std.log.debug("phosphorus stage: surface_biogeochemistry delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
            diagnostic_previous_p_g = current_p_g;
        }
    }
    try ecosys.soil_microbial_inventory_bridge.publishToOrganic(
        context.soil_microbial,
        context.soil_organic,
    );
    if (diagnostic_first_hour) {
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug(
            "nitrogen closed stage: biogeochemistry_and_mirror delta_g={e}",
            .{current_n_g - diagnostic_biogeochemistry_n_before_g},
        );
    }
    try context.mineral_nitrogen_transport.refreshMatrixFromReactionState(
        context.soil_chemistry,
        context.soil_reactive_nitrogen,
        context.grid.matrix_liquid_water_m3,
        context.mineral_nitrogen_zone_fractions,
        context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
    );
    if (diagnostic_first_hour) {
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug("nitrogen stage: pre_chemistry_matrix_refresh delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
        diagnostic_previous_n_g = current_n_g;
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: pre_chemistry_matrix_refresh delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        diagnostic_previous_p_g = current_p_g;
    }
    // Legacy SOLUTE follows the soil and litter biological source/sink commits.
    // Converge locally once; do not repeat a full sub-hourly model cycle.
    const diagnostic_chemistry_n_before_g = if (diagnostic_first_hour)
        try diagnostics.diagnosticStoredNitrogen_g(context)
    else
        0;
    try soil_chemistry_convergence.convergeHourlySoilChemistry(
        context,
        fertilizer_band_hour,
        solute_failure_report,
    );
    if (diagnostic_first_hour) {
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: soil_chemistry_only delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        diagnostic_previous_p_g = current_p_g;
    }
    try surface_litter_convergence.convergeSurfaceLitterChemistry(context);
    try context.mineral_nitrogen_transport.refreshMatrixFromReactionState(
        context.soil_chemistry,
        context.soil_reactive_nitrogen,
        context.grid.matrix_liquid_water_m3,
        context.mineral_nitrogen_zone_fractions,
        context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
    );
    if (diagnostic_first_hour) {
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug(
            "nitrogen closed stage: chemistry delta_g={e}",
            .{current_n_g - diagnostic_chemistry_n_before_g},
        );
    }
    if (diagnostic_first_hour) {
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug("nitrogen stage: chemistry delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
        diagnostic_previous_n_g = current_n_g;
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: chemistry delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        diagnostic_previous_p_g = current_p_g;
    }
    const diagnostic_pond_before = if (diagnostic_first_hour) try diagnostics.reconstructLandscapeMassBalance(context) else undefined;
    const diagnostic_ammonium_before = if (diagnostic_first_hour) try diagnostics.diagnosticAmmoniumOwners_g_n(context) else undefined;
    var pond_transition_context: ecosys.surface_pond_transition_step.ApplyContext = .{
        .result = context.surface_pond_transition,
        .surface_liquid_water_m3 = context.surface_precipitation.litter_water_m3,
        .surface_ice_m3 = context.surface_litter_ice_m3,
        .surface_ponding_capacity_m3 = context.surface_precipitation.litter_water_capacity_m3,
        .surface_litter_volume_m3 = context.surface_litter_geometry.dry_litter_volume_m3,
        .surface_litter_water_capacity_m3 = context.surface_litter_geometry.water_retention_capacity_m3,
        .horizontal_area_m2 = context.canopy_cell_area_m2,
        .minimum_heat_capacity_megajoules_per_k = context.surface_pond_minimum_heat_capacity_megajoules_per_k,
        .liquid_water_heat_capacity_megajoules_per_m3_k = context.runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_megajoules_per_m3_k,
    };
    try tile_kernels.runKernelAcrossSerialTiles(context, &pond_transition_context, ecosys.surface_pond_transition_step.applyTile);
    const diagnostic_heat_before_settling_megajoules = if (diagnostic_first_hour)
        (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules
    else
        0;
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
        .surface_sediment_megagrams = context.surface_erosion.surface_sediment_megagrams,
        .surface_soil_mass_megagrams = context.surface_erosion.surface_soil_mass_megagrams,
        .settled_sediment_megagrams = context.surface_erosion.pond_settled_sediment_megagrams,
    }, context.surface_pond_transition, 1);
    if (diagnostic_first_hour) {
        const components = try diagnostics.reconstructLandscapeMassBalance(context);
        std.log.debug("heat stage: pre_pond_processes hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, diagnostic_heat_before_settling_megajoules - diagnostic_previous_heat_megajoules });
        std.log.debug("heat stage: particulate_settling hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, components.heat_storage_megajoules - diagnostic_heat_before_settling_megajoules });
        diagnostic_previous_heat_megajoules = components.heat_storage_megajoules;
        std.log.debug("settling nitrogen components: residue={e} organic={e} n2={e} nh4={e} no3={e}", .{ components.residue_nitrogen_g - diagnostic_pond_before.residue_nitrogen_g, components.organic_nitrogen_g - diagnostic_pond_before.organic_nitrogen_g, components.dinitrogen_nitrogen_g - diagnostic_pond_before.dinitrogen_nitrogen_g, components.ammonium_nitrogen_g - diagnostic_pond_before.ammonium_nitrogen_g, components.nitrate_nitrogen_g - diagnostic_pond_before.nitrate_nitrogen_g });
    }
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
        .micropore_solutes = context.micropore_solute_state,
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
            .dry_organic_heat_capacity_megajoules_per_g_c_k = context.runscript.surface_pond_dry_organic_heat_capacity_megajoules_per_g_c_k,
            .liquid_water_heat_capacity_megajoules_per_m3_k = context.runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_megajoules_per_m3_k,
            .ice_heat_capacity_megajoules_per_m3_k = context.runscript.soil_phase_heat_parameters.ice_heat_capacity_megajoules_per_m3_k,
            .minimum_heat_capacity_megajoules_per_k = 0,
        },
        .minimum_heat_capacity_megajoules_per_k = context.surface_pond_minimum_heat_capacity_megajoules_per_k,
        .minimum_soil_layer_thickness_m = context.runscript.soil_geometry_parameters.minimum_layer_thickness_m,
        .horizontal_cell_width_m = context.horizontal_cell_width_m,
        .vertical_cell_width_m = context.vertical_cell_width_m,
        .ammonium_non_band_water_fraction = context.mineral_nitrogen_zone_fractions.ammonium_non_band,
        .nitrate_non_band_water_fraction = context.mineral_nitrogen_zone_fractions.nitrate_non_band,
    });
    if (diagnostic_first_hour) {
        const heat_megajoules = (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules;
        std.log.debug("heat stage: pond_domain_transaction hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, heat_megajoules - diagnostic_previous_heat_megajoules });
        diagnostic_previous_heat_megajoules = heat_megajoules;
    }
    if (diagnostic_first_hour) {
        const components = try diagnostics.reconstructLandscapeMassBalance(context);
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug("nitrogen stage: pond_transaction delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: pond_transaction delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        diagnostic_previous_p_g = current_p_g;
        std.log.debug("pond nitrogen components: residue={e} organic={e} n2={e} nh4={e} no3={e}", .{ components.residue_nitrogen_g - diagnostic_pond_before.residue_nitrogen_g, components.organic_nitrogen_g - diagnostic_pond_before.organic_nitrogen_g, components.dinitrogen_nitrogen_g - diagnostic_pond_before.dinitrogen_nitrogen_g, components.ammonium_nitrogen_g - diagnostic_pond_before.ammonium_nitrogen_g, components.nitrate_nitrogen_g - diagnostic_pond_before.nitrate_nitrogen_g });
        const owners = try diagnostics.diagnosticAmmoniumOwners_g_n(context);
        std.log.debug("pond ammonium owners delta: surface_aq={e} surface_exchange={e} surface_fertilizer={e} soil_aq={e} soil_exchange={e} soil_fertilizer={e}", .{ owners[0] - diagnostic_ammonium_before[0], owners[1] - diagnostic_ammonium_before[1], owners[2] - diagnostic_ammonium_before[2], owners[3] - diagnostic_ammonium_before[3], owners[4] - diagnostic_ammonium_before[4], owners[5] - diagnostic_ammonium_before[5] });
        diagnostic_previous_n_g = current_n_g;
    }
    const diagnostic_pond_after = if (diagnostic_first_hour) try diagnostics.reconstructLandscapeMassBalance(context) else undefined;
    // Snow discharge, nonlinear chemistry, and pond-domain transfers update
    // the concentration owner after hourly mineral-N transport has finished.
    // Rebase its extensive matrix mirror now so end-of-hour audits and restart
    // persistence observe the accepted state rather than the hour-start copy.
    try context.mineral_nitrogen_transport.refreshMatrixFromReactionState(
        context.soil_chemistry,
        context.soil_reactive_nitrogen,
        context.grid.matrix_liquid_water_m3,
        context.mineral_nitrogen_zone_fractions,
        context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
    );
    // Pond and settling steps transfer surface_organic.microbial into soil_organic.microbial
    // but cannot reach soil_microbial directly.  Propagate the update now so that
    // soil_microbial and soil_organic.microbial remain in sync before the next d-step.
    try ecosys.soil_microbial_inventory_bridge.publishFromOrganic(
        context.soil_organic,
        context.soil_microbial,
    );
    if (diagnostic_first_hour) {
        const components = try diagnostics.reconstructLandscapeMassBalance(context);
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug("nitrogen stage: post_pond_microbial_refresh delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: post_pond_microbial_refresh delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        std.log.debug("phosphorus end science: stored_g={e}", .{current_p_g});
        diagnostic_previous_p_g = current_p_g;
        std.log.debug("post-pond refresh nitrogen components: residue={e} organic={e} n2={e} nh4={e} no3={e}", .{ components.residue_nitrogen_g - diagnostic_pond_after.residue_nitrogen_g, components.organic_nitrogen_g - diagnostic_pond_after.organic_nitrogen_g, components.dinitrogen_nitrogen_g - diagnostic_pond_after.dinitrogen_nitrogen_g, components.ammonium_nitrogen_g - diagnostic_pond_after.ammonium_nitrogen_g, components.nitrate_nitrogen_g - diagnostic_pond_after.nitrate_nitrogen_g });
        diagnostic_previous_n_g = current_n_g;
    }
    try finalizeVegetationAndLedgers(
        context,
        hour_of_day,
        weather_header_by_cell,
        plant_calendar_by_cell,
        diagnostic_first_hour,
        diagnostic_previous_p_g,
        plant_calendar,
        snow_phase_change_report,
        snow_vapor_equilibrium_report,
        subsurface_irrigation_chemistry_parameters,
    );
}

fn transportDissolvedGasAndSurfaceWater(
    context: anytype,
    hour_of_day: u8,
    weather_header_by_cell: []const ecosys.weather.Header,
    plant_calendar_by_cell: []const ecosys.plant_development.Calendar,
    fertilizer_band_hour: ecosys.fertilizer_band_phase_coordinator.HourToken,
    gas_failure_report: ?ecosys.soil_gas_transport_step.FailureReportRequest,
    solute_failure_report: ?ecosys.solute_failure_reporter.Request,
    diagnostic_first_hour: anytype,
    plant_calendar: anytype,
    snow_phase_change_report: anytype,
    snow_vapor_equilibrium_report: anytype,
    subsurface_irrigation_chemistry_parameters: anytype,
    diagnostic_previous_heat_megajoules_ptr: anytype,
    diagnostic_previous_n_g_ptr: anytype,
    diagnostic_previous_p_g_ptr: anytype,
) !void {
    var diagnostic_previous_heat_megajoules = diagnostic_previous_heat_megajoules_ptr.*;
    defer diagnostic_previous_heat_megajoules_ptr.* = diagnostic_previous_heat_megajoules;
    var diagnostic_previous_n_g = diagnostic_previous_n_g_ptr.*;
    defer diagnostic_previous_n_g_ptr.* = diagnostic_previous_n_g;
    var diagnostic_previous_p_g = diagnostic_previous_p_g_ptr.*;
    defer diagnostic_previous_p_g_ptr.* = diagnostic_previous_p_g;
    try context.soil_dissolved_gas_face_parameters.refresh(
        context.grid,
        context.soil_transport_faces,
        context.soil_face_geometry,
        context.soil_solver_properties.matrix_bulk_volume_m3,
        1,
        .{},
    );
    // Non-convergence rolls back dissolved-gas state to the pre-step snapshot
    // (see aqueous_extensive_transport.advance defer block). Log and continue;
    // the transport re-attempts next hour as water-balance conditions change.
    const diagnostic_dissolved_transport_n_before_g = if (diagnostic_first_hour)
        try diagnosticGasNitrogen_g(context.gas_transport)
    else
        0;
    _ = ecosys.soil_dissolved_gas_transport.advance(
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
    ) catch |err| switch (err) {
        error.AqueousExtensiveTransportDidNotConverge,
        error.NonFiniteAqueousExtensiveTransport,
        => std.log.warn("dissolved-gas transport non-convergence retained: err={s}", .{@errorName(err)}),
        else => return err,
    };
    if (diagnostic_first_hour) {
        var diagnostic_boundary_n_g: f64 = 0;
        for (0..context.grid.layer_count) |layer| {
            const base = layer * ecosys.gas_transport.species_count;
            inline for ([_]ecosys.gas_transport.Species{ .nitrogen, .nitrous_oxide, .ammonia }) |species|
                diagnostic_boundary_n_g += context.soil_dissolved_gas_transport.boundary_net_flux_g[base + @intFromEnum(species)];
        }
        const diagnostic_after_g = try diagnosticGasNitrogen_g(context.gas_transport);
        std.log.debug("dissolved gas nitrogen transaction: delta_g={e} boundary_g={e} residual_g={e}", .{ diagnostic_after_g - diagnostic_dissolved_transport_n_before_g, diagnostic_boundary_n_g, diagnostic_after_g - diagnostic_dissolved_transport_n_before_g - diagnostic_boundary_n_g });
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug("nitrogen stage: dissolved_gas_transport delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
        diagnostic_previous_n_g = current_n_g;
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: dissolved_gas_transport delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        diagnostic_previous_p_g = current_p_g;
        const heat_megajoules = (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules;
        std.log.debug("heat stage: dissolved_gas_transport hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, heat_megajoules - diagnostic_previous_heat_megajoules });
        diagnostic_previous_heat_megajoules = heat_megajoules;
    }
    {
        const surface_parameters = context.surface_gas_parameters.*;
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
            // `hour1.f` 4688--4694 derives THREE distinct vapor diffusivities
            // from three distinct temperatures, and they are not
            // interchangeable: `WGSGA` from canopy air `TKQ` (4688--4689),
            // `WGSGR` from litter `TKS(0)` (4690--4691), and `WGSGW` from
            // snowpack `TKW` (4693--4694). `watsub.f` then consumes `WGSGR`
            // specifically as the litter diffusivity, at 656, 1847, and 2975.
            // The litter porous resistance below was using the air-temperature
            // diffusivity, which is the `WGSGA` role, so it silently applied the
            // canopy-air temperature to a litter transport path. The two differ
            // by 13.2 % at a 20 K air-litter offset and the error reverses sign
            // when litter is warmer than air, so it does not average out over a
            // diurnal cycle. `grid.surface_temperature_k` is production's
            // established `TKS(0)` analogue, as at lines 2751, 2766, and 5236.
            const litter_diffusivity_m2_per_h =
                context.runscript.soil_process_parameters.reference_water_vapor_diffusivity_m2_per_h *
                std.math.pow(f64, context.grid.surface_temperature_k[cell] / context.runscript.soil_process_parameters.vapor_diffusivity_reference_temperature_k, context.runscript.soil_process_parameters.vapor_diffusivity_temperature_exponent);
            const litter_porous_resistance_h_per_m = if (litter_volume_m3 > 0 and litter_air_m3 > 0 and litter_porosity > 0) blk: {
                const air_fraction = std.math.clamp(litter_air_m3 / litter_volume_m3, 0, litter_porosity);
                const transport_factor = @max(
                    context.runscript.surface_gas_resistance_parameters.minimum_air_fraction,
                    context.runscript.soil_gas_transport_parameters.penman_tortuosity * air_fraction * air_fraction / litter_porosity,
                );
                break :blk (litter_volume_m3 / area_m2) / litter_diffusivity_m2_per_h / transport_factor;
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
        if (diagnostic_first_hour) {
            const heat_megajoules = (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules;
            std.log.debug("heat stage: before_coupled_gas hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, heat_megajoules - diagnostic_previous_heat_megajoules });
            diagnostic_previous_heat_megajoules = heat_megajoules;
        }
        coupled_gas: {
            const diagnostic_n_before_g = if (context.executed_weather_hours.* < 24)
                try diagnosticGasNitrogen_g(context.gas_transport)
            else
                0;
            var accepted_gas_state = try context.gas_transport.clone(
                context.allocator,
            );
            defer accepted_gas_state.deinit();
            _ = context.soil_gas_transport.advance(.{
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
            }) catch |err| switch (err) {
                // Non-convergence discards the accepted_gas_state clone, leaving the
                // pre-step gas state in context.gas_transport. The commit is skipped
                // via break; the solver retries next hour as water-balance conditions change.
                error.CoupledGasSolverDidNotConverge => {
                    std.log.warn("coupled-gas solver non-convergence retained: err={s}", .{@errorName(err)});
                    break :coupled_gas;
                },
                else => return err,
            };
            const diagnostic_n_accepted_g = if (context.executed_weather_hours.* < 24)
                try diagnosticGasNitrogen_g(&accepted_gas_state)
            else
                0;
            try commitHourlyGasContributionGeneration(
                context,
                &accepted_gas_state,
            );
            if (context.executed_weather_hours.* < 24) {
                var diagnostic_boundary_n_g: f64 = 0;
                inline for ([_]ecosys.gas_transport.Species{ .nitrogen, .nitrous_oxide, .ammonia }) |species| {
                    for (0..context.grid.layer_count) |layer| {
                        diagnostic_boundary_n_g += context.soil_gas_transport.atmospheric_flux_g_per_h[
                            layer * ecosys.gas_transport.species_count + @intFromEnum(species)
                        ];
                        diagnostic_boundary_n_g += context.soil_gas_transport.subsurface_flux_g_per_h[
                            layer * ecosys.gas_transport.species_count + @intFromEnum(species)
                        ];
                    }
                }
                std.log.debug(
                    "gas nitrogen transaction: hour={d} before={e} accepted={e} committed={e} delta={e} boundary={e} residual={e}",
                    .{
                        context.executed_weather_hours.* + 1,
                        diagnostic_n_before_g,
                        diagnostic_n_accepted_g,
                        try diagnosticGasNitrogen_g(context.gas_transport),
                        diagnostic_n_accepted_g - diagnostic_n_before_g,
                        diagnostic_boundary_n_g,
                        diagnostic_n_accepted_g - diagnostic_n_before_g - diagnostic_boundary_n_g,
                    },
                );
            }
        }
    }
    if (diagnostic_first_hour) {
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug("nitrogen stage: coupled_gas_transport delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
        diagnostic_previous_n_g = current_n_g;
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: coupled_gas_transport delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        diagnostic_previous_p_g = current_p_g;
        const heat_megajoules = (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules;
        std.log.debug("heat stage: through_coupled_gas hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, heat_megajoules - diagnostic_previous_heat_megajoules });
        diagnostic_previous_heat_megajoules = heat_megajoules;
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
        context.config.lon_count,
        context.config.lat_count,
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
        context.config.lon_count,
        context.config.lat_count,
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
        context.config.lon_count,
        context.config.lat_count,
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
        context.config.lon_count,
        context.config.lat_count,
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
    if (diagnostic_first_hour) {
        const heat_megajoules = (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules;
        std.log.debug("heat stage: surface_runoff_and_dissolved_gas hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, heat_megajoules - diagnostic_previous_heat_megajoules });
        diagnostic_previous_heat_megajoules = heat_megajoules;
    }
    if (diagnostic_first_hour) {
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug("nitrogen stage: surface_runoff_transport delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
        diagnostic_previous_n_g = current_n_g;
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: surface_runoff_transport delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        diagnostic_previous_p_g = current_p_g;
    }
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
                .reference_water_viscosity_megagrams_per_m_s = context.runscript.soil_process_parameters.reference_water_viscosity_megagrams_per_m_s,
                .viscosity_temperature_intercept = context.runscript.soil_process_parameters.water_viscosity_temperature_intercept,
                .viscosity_temperature_coefficient_per_c = context.runscript.soil_process_parameters.water_viscosity_temperature_coefficient_per_c,
            },
        );
        const soil_mass_megagrams = context.soil_solver_properties.matrix_bulk_volume_m3[top_layer] * context.soil_solver_properties.bulk_density_megagrams_per_m3[top_layer];
        if (!context.surface_erosion.surface_soil_mass_initialized[cell]) {
            context.surface_erosion.surface_soil_mass_megagrams[cell] = soil_mass_megagrams;
            context.surface_erosion.surface_soil_mass_initialized[cell] = true;
        }
        context.surface_soil_mass_at_erosion_start_megagrams[cell] = context.surface_erosion.surface_soil_mass_megagrams[cell];
        const local_solve =
            try ecosys.soil_erosion.calculateConvergedHourlyLocalStep(.{
                .erosion_enabled = context.site_by_cell[cell].erosionEnabled(),
                .surface_soil_bulk_density_megagrams_per_m3 = context.soil_solver_properties.bulk_density_megagrams_per_m3[top_layer],
                .surface_soil_mass_megagrams = context.surface_erosion.surface_soil_mass_megagrams[cell],
                .surface_soil_water_m3 = context.grid.matrix_liquid_water_m3[top_layer],
                .surface_soil_pore_volume_m3 = context.grid.matrix_pore_capacity_m3[top_layer],
                .excess_surface_water_m3 = context.surface_runoff.excess_surface_water_m3[cell],
                .excess_surface_ice_m3 = context.surface_runoff.excess_surface_ice_m3[cell],
                .surface_ponding_capacity_m3 = @max(context.config.absolute_tolerance, context.runscript.surface_runoff_parameters.ground_surface_retention_m3_per_m2 * context.canopy_cell_area_m2[cell]),
                .sediment_in_surface_water_megagrams = context.surface_erosion.surface_sediment_megagrams[cell],
                .rainfall_kinetic_energy_j = context.surface_precipitation.rainfall_impact_energy_j[cell],
                .soil_rainfall_detachability_g_per_j = erosion_properties.rainfall_detachability_g_per_j,
                .soil_runoff_detachability = erosion_properties.runoff_detachability,
                .sediment_settling_velocity_m_per_h = erosion_properties.settling_velocity_m_per_h,
                .grid_cell_area_m2 = context.canopy_cell_area_m2[cell],
                .soil_matrix_fraction = std.math.clamp(context.soil_solver_properties.matrix_bulk_volume_m3[top_layer] / context.soil_solver_properties.layer_volume_m3[top_layer], 0, 1),
                .snow_free_fraction = 1 - std.math.clamp(context.surface_precipitation.snow_cover_fraction[cell], 0, 1),
                .runoff_velocity_m_per_s = context.surface_runoff.runoff_velocity_m_per_s[cell],
                .slope_sine = context.terrain_hydrology.slope_m_per_m[cell],
                .surface_particle_density_megagrams_per_m3 = erosion_properties.particle_density_megagrams_per_m3,
                .transport_capacity_coefficient = erosion_properties.transport_capacity_coefficient,
                .transport_capacity_exponent = erosion_properties.transport_capacity_exponent,
                .maximum_erodible_soil_fraction_per_step = 1,
                .water_transport_timestep_h = 1,
                .negligible_volume_m3 = context.runscript.surface_runoff_parameters.negligible_water_m3,
                .negligible_mass_megagrams = context.runscript.surface_runoff_parameters.negligible_water_m3,
            }, .{
                .absolute_tolerance_megagrams = context.config.absolute_tolerance,
                .relative_tolerance = context.config.relative_tolerance,
                .picard_relaxation = context.config.picard_relaxation,
                .max_iterations = context.iteration_limits.erosion_max_iterations,
            });
        const local = local_solve.local;
        context.surface_erosion.local_detachment_megagrams[cell] = local.net_detachment_megagrams;
        context.surface_erosion.transportable_sediment_megagrams[cell] = if (context.site_by_cell[cell].erosionEnabled())
            try ecosys.soil_erosion.calculateDownslopeTransport_megagrams(
                context.surface_erosion.surface_sediment_megagrams[cell],
                local.net_detachment_megagrams,
                context.surface_runoff.excess_surface_water_m3[cell],
                context.surface_runoff.total_runoff_m3_per_step[cell],
                context.soil_solver_properties.bulk_density_megagrams_per_m3[top_layer],
                context.runscript.surface_runoff_parameters.negligible_water_m3,
            )
        else
            0;
        const column = cell % context.config.lon_count;
        const row = cell / context.config.lon_count;
        context.surface_erosion.east_boundary_open[cell] = column + 1 == context.config.lon_count and context.site_by_cell[cell].surface_runoff_boundary_fraction[1] > 0;
        context.surface_erosion.west_boundary_open[cell] = column == 0 and context.site_by_cell[cell].surface_runoff_boundary_fraction[3] > 0;
        context.surface_erosion.south_boundary_open[cell] = row + 1 == context.config.lat_count and context.site_by_cell[cell].surface_runoff_boundary_fraction[2] > 0;
        context.surface_erosion.north_boundary_open[cell] = row == 0 and context.site_by_cell[cell].surface_runoff_boundary_fraction[0] > 0;
    }
    try routeSedimentAndErosion(
        context,
        hour_of_day,
        weather_header_by_cell,
        plant_calendar_by_cell,
        fertilizer_band_hour,
        solute_failure_report,
        diagnostic_first_hour,
        plant_calendar,
        snow_phase_change_report,
        snow_vapor_equilibrium_report,
        subsurface_irrigation_chemistry_parameters,
        &diagnostic_previous_heat_megajoules,
        &diagnostic_previous_n_g,
        &diagnostic_previous_p_g,
    );
}

fn solveSnowSurfaceEnergyAndSoilTransport(
    context: anytype,
    hour_of_day: u8,
    weather_header_by_cell: []const ecosys.weather.Header,
    plant_calendar_by_cell: []const ecosys.plant_development.Calendar,
    fertilizer_band_hour: ecosys.fertilizer_band_phase_coordinator.HourToken,
    gas_failure_report: ?ecosys.soil_gas_transport_step.FailureReportRequest,
    solute_failure_report: ?ecosys.solute_failure_reporter.Request,
    diagnostic_first_hour: anytype,
    diagnostic_mineral_before_mol: anytype,
    diagnostic_relayer_phosphate_before: anytype,
    diagnostic_transport_ammonium_before: anytype,
    diagnostic_transport_before: anytype,
    plant_calendar: anytype,
    snow_phase_change_report: anytype,
    snow_vapor_equilibrium_report: anytype,
    diagnostic_previous_heat_megajoules_ptr: anytype,
    diagnostic_previous_n_g_ptr: anytype,
    diagnostic_previous_p_g_ptr: anytype,
    diagnostic_previous_p_owners_ptr: anytype,
) !void {
    var diagnostic_previous_heat_megajoules = diagnostic_previous_heat_megajoules_ptr.*;
    defer diagnostic_previous_heat_megajoules_ptr.* = diagnostic_previous_heat_megajoules;
    var diagnostic_previous_n_g = diagnostic_previous_n_g_ptr.*;
    defer diagnostic_previous_n_g_ptr.* = diagnostic_previous_n_g;
    var diagnostic_previous_p_g = diagnostic_previous_p_g_ptr.*;
    defer diagnostic_previous_p_g_ptr.* = diagnostic_previous_p_g;
    var diagnostic_previous_p_owners = diagnostic_previous_p_owners_ptr.*;
    defer diagnostic_previous_p_owners_ptr.* = diagnostic_previous_p_owners;
    try context.snow_transport.commitMeltWater(context.transport_hydrology.snow_downward_water_flux_m3_per_step, context.transport_hydrology.snow_to_litter_water_flux_m3_per_step, context.transport_hydrology.snow_to_soil_micropore_flux_m3_per_step, context.transport_hydrology.snow_to_soil_macropore_flux_m3_per_step);
    var aggregate_fluxes = try ecosys.snow_solute_transport.calculateFluxes(
        context.allocator,
        context.snow_transport,
        context.transport_hydrology.snow_downward_water_flux_m3_per_step,
        context.transport_hydrology.snow_to_litter_water_flux_m3_per_step,
        context.transport_hydrology.snow_to_soil_micropore_flux_m3_per_step,
        context.transport_hydrology.snow_to_soil_macropore_flux_m3_per_step,
        context.snow_surface_partitions,
    );
    defer aggregate_fluxes.deinit();
    const layer_capacity = context.snow_transport.layer_capacity;
    const snow_layer_count = try std.math.mul(usize, context.grid.cell_count, layer_capacity);
    const salt_value_count = try std.math.mul(
        usize,
        snow_layer_count,
        ecosys.snowpack_internal_salt_aggregation.salt_species_count,
    );
    const salt_layer_value_count = try std.math.mul(
        usize,
        layer_capacity,
        ecosys.snowpack_internal_salt_aggregation.salt_species_count,
    );
    @memset(context.snowpack_internal_solute_flux_by_layer, .{});
    @memset(context.snowpack_internal_solute_flux_workspace, .{});
    @memset(context.snowpack_internal_salt_flux_mol_by_layer_species, 0);
    @memset(context.snowpack_internal_salt_flux_workspace_by_layer_species, 0);
    for (0..context.grid.cell_count) |cell| {
        const cell_layer_start = cell * layer_capacity;
        for (0..context.snow_transport.layer_capacity) |layer| {
            const face_layer_start = cell_layer_start + layer;
            if (layer > 0) {
                const source_layer_start =
                    (cell_layer_start + layer - 1) * ecosys.snow_solute_transport.species_count;
                const source = aggregate_fluxes.downward_g[source_layer_start .. source_layer_start + ecosys.snow_solute_transport.species_count];
                context.snowpack_internal_solute_flux_by_layer[face_layer_start] =
                    snowpackInternalNonSaltFluxFromSpeciesAmounts(source);
            }
        }
        var solute_state: ecosys.snowpack_internal_solute_aggregation.State = .{
            .net_flux_by_layer = context.snowpack_internal_solute_flux_by_layer[cell_layer_start .. cell_layer_start + layer_capacity],
        };
        const solute_workspace: ecosys.snowpack_internal_solute_aggregation.Workspace = .{
            .net_flux_by_layer = context.snowpack_internal_solute_flux_workspace[cell_layer_start .. cell_layer_start + layer_capacity],
        };
        try ecosys.snowpack_internal_solute_aggregation.aggregate(.{
            .heat_capacity_megajoules_per_k_by_layer = context.snow_transport.heat_capacity_megajoules_per_k[cell_layer_start .. cell_layer_start + layer_capacity],
            .minimum_heat_capacity_megajoules_per_k = 1,
            .upper_face_flux_by_layer = context.snowpack_internal_solute_flux_by_layer[cell_layer_start .. cell_layer_start + layer_capacity],
        }, &solute_state, solute_workspace);
        const cell_salt_start = cell_layer_start * ecosys.snowpack_internal_salt_aggregation.salt_species_count;
        const cell_salt_end = @min(
            salt_value_count,
            cell_salt_start + salt_layer_value_count,
        );
        var salt_state: ecosys.snowpack_internal_salt_aggregation.State = .{
            .net_mol_per_step_by_layer_species = context.snowpack_internal_salt_flux_mol_by_layer_species[cell_salt_start..cell_salt_end],
        };
        const salt_workspace: ecosys.snowpack_internal_salt_aggregation.Workspace = .{
            .net_mol_per_step_by_layer_species = context.snowpack_internal_salt_flux_workspace_by_layer_species[cell_salt_start..cell_salt_end],
        };
        try ecosys.snowpack_internal_salt_aggregation.aggregate(.static, .{
            .heat_capacity_megajoules_per_k_by_layer = context.snow_transport.heat_capacity_megajoules_per_k[cell_layer_start .. cell_layer_start + layer_capacity],
            .minimum_heat_capacity_megajoules_per_k = 1,
            .upper_face_mol_per_step_by_layer_species = context.snowpack_internal_salt_flux_mol_by_layer_species[cell_salt_start..cell_salt_end],
        }, &salt_state, salt_workspace);
    }
    try ecosys.snow_compaction.apply(context.allocator, context.snow_transport, .{
        .snowfall_water_equivalent_m3 = context.surface_precipitation.snow_to_snow_m3_per_h,
        .atmospheric_temperature_k = context.atmosphere.air_temperature_k,
        .timestep_h = 1,
        .initial_snow_density_megagrams_per_m3 = context.runscript.initial_snow_density_megagrams_per_m3,
        .ice_density_megagrams_per_m3 = context.runscript.snow_ice_density_megagrams_per_m3,
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
        // HOUR1-002. `hour1.f` 4713--4779 sums `ARLSS` over leaf, stalk
        // (`ARSTK`), and standing dead (`ARSTD`), and `DO 145`/`DO 155` form
        // `FRADT` from the per-plant `FRADP`/`FRADQ` shares. Passing the
        // retention owner's fractions retires `canopy_exposure`'s own
        // leaf-area-only derivation, which omitted stalk and standing-dead area
        // and hard-coded the `0.65` extinction and `0.05` solar-angle gate that
        // are runtime controls. This RETIRES a producer rather than adding one:
        // `canopy_exposure.State` has exactly one writer before and after, and
        // what changes is which upstream quantity that writer reads.
        //
        // Ordering verified rather than assumed, since a producer/consumer
        // inversion is what blocks BIND-REDIST-ROGOX:
        // `canopy_precipitation_retention.refreshFromModel` publishes
        // `living_radiation_fraction`/`standing_dead_radiation_fraction` at
        // line 3509 of this file, in the same hourly pass and 265 lines before
        // this call site, so the values read here are current for this hour.
        var exposure_context: ecosys.canopy_exposure.ApplyContext = .{ .result = exposure, .structure = &context.canopy_structure.*.?, .interception = &context.canopy_interception.*.?, .ground_radiation = context.ground_radiation, .solar_angle_sine_by_cell = context.hourly_solar_angle_sine, .radiation_fractions = if (context.canopy_precipitation_retention.*) |*retention| .{
            .living_radiation_fraction = retention.living_radiation_fraction,
            .standing_dead_radiation_fraction = retention.standing_dead_radiation_fraction,
            .species_count = retention.species_count,
        } else null };
        try tile_kernels.runKernelAcrossSerialTiles(context, &exposure_context, ecosys.canopy_exposure.applyTile);
    }
    var surface_energy_context: ecosys.surface_energy.ApplyContext = .{ .result = context.surface_energy, .grid = context.grid, .atmosphere = context.atmosphere, .ground_radiation = context.ground_radiation, .snow_depth_m = context.snow_depth_m, .exposure = if (context.canopy_exposure.*) |*exposure| exposure else null, .settings = context.surface_energy_settings };
    try tile_kernels.runKernelAcrossSerialTiles(context, &surface_energy_context, ecosys.surface_energy.applyTile);
    for (0..context.grid.cell_count) |cell| {
        const area_m2 = context.canopy_cell_area_m2[cell];
        if (!std.math.isFinite(area_m2) or area_m2 <= 0) return error.InvalidSurfaceFireCellArea;
        context.surface_combustion_heat_megajoules_per_m2[cell] = context.delayed_surface_combustion_heat_megajoules[cell] / area_m2;
        if (!std.math.isFinite(context.surface_combustion_heat_megajoules_per_m2[cell]) or context.surface_combustion_heat_megajoules_per_m2[cell] < 0) return error.NonFiniteSurfaceCombustionHeatSource;
    }
    for (0..context.grid.cell_count) |cell| {
        const vapor_water_equivalent_m3 = context.litter_gas_transport.water_vapor_mol[cell] * context.runscript.soil_gas_transport_parameters.water_molar_mass_g_per_mol / context.runscript.soil_gas_transport_parameters.water_density_g_per_m3;
        context.surface_heat_capacity_megajoules_per_k[cell] = context.runscript.surface_pond_dry_organic_heat_capacity_megajoules_per_g_c_k * try context.surface_organic.totalCarbon_g_c(cell) + context.runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_megajoules_per_m3_k * (context.surface_precipitation.litter_water_m3[cell] + vapor_water_equivalent_m3) + context.runscript.soil_phase_heat_parameters.ice_heat_capacity_megajoules_per_m3_k * context.surface_litter_ice_m3[cell];
        if (!std.math.isFinite(context.surface_heat_capacity_megajoules_per_k[cell]) or context.surface_heat_capacity_megajoules_per_k[cell] <= 0) return error.InvalidSurfaceHeatCapacity;
    }
    var surface_temperature_context: ecosys.surface_temperature_solver.ApplyContext = .{ .result = context.surface_temperature, .grid = context.grid, .atmosphere = context.atmosphere, .air_temperature_k = context.ground_air.temperature_k, .air_vapor_pressure_kpa = context.ground_air_vapor_pressure_kpa, .ground_radiation = context.ground_radiation, .surface_energy = context.surface_energy, .soil_thermal = context.soil_thermal, .surface_heat_capacity_megajoules_per_k = context.surface_heat_capacity_megajoules_per_k, .exposure = if (context.canopy_exposure.*) |*exposure| exposure else null, .external_heat_megajoules_per_m2 = context.surface_combustion_heat_megajoules_per_m2, .surface_phase = .{ .liquid_water_m3 = context.surface_precipitation.litter_water_m3, .ice_water_equivalent_m3 = context.surface_litter_ice_m3, .retention_capacity_m3 = context.surface_precipitation.litter_water_capacity_m3, .horizontal_area_m2 = context.canopy_cell_area_m2, .residual_water_content_m3_per_m3 = context.runscript.soil_process_parameters.surface_residue_residual_water_content_m3_per_m3, .van_genuchten_alpha_per_m = context.runscript.soil_process_parameters.surface_residue_van_genuchten_alpha_per_m, .van_genuchten_n = context.runscript.soil_process_parameters.surface_residue_van_genuchten_n, .gravitational_water_potential_mpa_per_m = context.runscript.soil_process_parameters.gravitational_water_potential_mpa_per_m, .latent_heat_of_fusion_megajoules_per_m3 = context.runscript.soil_phase_heat_parameters.freeze_thaw.latent_heat_of_fusion_megajoules_per_m3, .pure_water_melting_temperature_k = context.runscript.soil_phase_heat_parameters.freeze_thaw.pure_water_freezing_temperature_k }, .settings = .{ .timestep_hours = 1.0, .sensible_heat_conductance_megajoules_per_m2_h_k = context.runscript.surface_sensible_heat_conductance_megajoules_per_m2_h_k, .latent_heat_conductance_megajoules_per_m2_h_kpa = context.runscript.surface_latent_heat_conductance_megajoules_per_m2_h_kpa, .liquid_water_heat_capacity_megajoules_per_m3_k = context.runscript.canopy_surface_exchange_parameters.liquid_water_heat_capacity_megajoules_per_m3_k, .latent_heat_of_vaporization_megajoules_per_m3 = context.runscript.ground_air_parameters.liquid_water_latent_heat_megajoules_per_m3, .surface_vapor_activity_fraction = context.runscript.surface_vapor_activity_fraction, .minimum_temperature_k = context.runscript.minimum_surface_temperature_k, .maximum_temperature_k = context.runscript.maximum_surface_temperature_k, .solver_options = context.nonlinear_solver_options } };
    // The surface temperature Newton solver needs extra iterations at the
    // 273.15 K melting kink; the NPH outer-loop ceiling (water_heat_solute_max_iterations)
    // is too tight for the inner Newton convergence near phase transitions.
    surface_temperature_context.settings.solver_options.max_iterations =
        @max(surface_temperature_context.settings.solver_options.max_iterations, 50);
    // The production path retains the preceding surface state when a tile
    // solve stagnates. Clear every cell first so failed and not-yet-visited
    // cells cannot publish a phase change left by the preceding hour.
    context.surface_temperature.resetPhaseChangeDiagnostics();
    tile_kernels.runKernelAcrossSerialTiles(context, &surface_temperature_context, ecosys.surface_temperature_solver.applyTile) catch |err| switch (err) {
        // Stagnation near phase-transition temperatures is caused by the water-table
        // overflow in cell 0; the previous hour's surface_temperature values are
        // retained unchanged. Surface energy balance retries next hour.
        error.NewtonPicardStagnated => std.log.warn(
            "surface temperature Newton stagnation retained: err={s}",
            .{@errorName(err)},
        ),
        else => return err,
    };
    try ecosys.surface_litter_chemistry_carrier_rebase.rebaseFromAcceptedLiquidWaterChange(
        context.surface_litter_chemistry,
        context.surface_precipitation.litter_water_m3,
        context.surface_temperature.liquid_water_change_m3,
    );
    if (diagnostic_first_hour) {
        const heat_megajoules = (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules;
        std.log.debug("heat stage: surface_solve hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, heat_megajoules - diagnostic_previous_heat_megajoules });
        // EXEC-002: hour 2 loses ~33x the neighbouring hours here while water
        // closes to 1e-9 m3, so the carriers move correctly and only their heat
        // content is wrong. Report the surface temperature and the litter/pond
        // carriers that `surface_solve` owns, to see which one jumps.
        var diagnostic_surface_temperature_sum_k: f64 = 0;
        var diagnostic_litter_water_m3: f64 = 0;
        var diagnostic_litter_ice_m3: f64 = 0;
        var diagnostic_surface_temperature_min_k: f64 = std.math.floatMax(f64);
        var diagnostic_surface_temperature_max_k: f64 = 0;
        for (0..context.grid.cell_count) |cell| {
            diagnostic_surface_temperature_sum_k += context.grid.surface_temperature_k[cell];
            diagnostic_litter_water_m3 += context.surface_precipitation.litter_water_m3[cell];
            diagnostic_litter_ice_m3 += context.surface_litter_ice_m3[cell];
            diagnostic_surface_temperature_min_k = @min(diagnostic_surface_temperature_min_k, context.grid.surface_temperature_k[cell]);
            diagnostic_surface_temperature_max_k = @max(diagnostic_surface_temperature_max_k, context.grid.surface_temperature_k[cell]);
        }
        std.log.debug("surface carriers: hour={d} mean_surface_temperature_k={e} min_k={e} max_k={e} litter_water_m3={e} litter_ice_m3={e}", .{ context.executed_weather_hours.* + 1, diagnostic_surface_temperature_sum_k / @as(f64, @floatFromInt(context.grid.cell_count)), diagnostic_surface_temperature_min_k, diagnostic_surface_temperature_max_k, diagnostic_litter_water_m3, diagnostic_litter_ice_m3 });
        // EXEC-004: the litter solute pools are water-normalized, so a carrier of
        // exactly zero is unrepresentable and currently a hard error. Report the
        // minimum per-cell carrier so it is visible how near production comes to
        // that boundary even with evaporation disabled.
        var diagnostic_min_litter_water_m3: f64 = std.math.floatMax(f64);
        for (0..context.grid.cell_count) |cell|
            diagnostic_min_litter_water_m3 = @min(diagnostic_min_litter_water_m3, context.surface_precipitation.litter_water_m3[cell]);
        std.log.debug("litter carrier floor: hour={d} min_litter_water_m3={e}", .{ context.executed_weather_hours.* + 1, diagnostic_min_litter_water_m3 });
        // EXEC-004: the surface maintenance step needs a strictly positive litter
        // hydrogen activity to form pH. Report its range so a zero can be
        // distinguished from a merely small value.
        var diagnostic_min_hydrogen_mol_per_m3: f64 = std.math.floatMax(f64);
        var diagnostic_max_hydrogen_mol_per_m3: f64 = 0;
        for (context.surface_litter_chemistry.cells) |litter_cell| {
            diagnostic_min_hydrogen_mol_per_m3 = @min(diagnostic_min_hydrogen_mol_per_m3, litter_cell.hydrogen_mol_per_m3);
            diagnostic_max_hydrogen_mol_per_m3 = @max(diagnostic_max_hydrogen_mol_per_m3, litter_cell.hydrogen_mol_per_m3);
        }
        std.log.debug("litter hydrogen activity: hour={d} min_mol_per_m3={e} max_mol_per_m3={e}", .{ context.executed_weather_hours.* + 1, diagnostic_min_hydrogen_mol_per_m3, diagnostic_max_hydrogen_mol_per_m3 });
        // EXEC-004: whether the litter carrier rebase has taken its hold-while-dry
        // path. Zero in the shipped configuration; nonzero once surface evaporation
        // is active. Reported so reachability is measured rather than assumed.
        std.log.debug("litter dry branch executions: hour={d} count={d}", .{ context.executed_weather_hours.* + 1, ecosys.surface_litter_chemistry_carrier_rebase.dry_branch_executions });
        diagnostic_previous_heat_megajoules = heat_megajoules;
    }
    @memset(context.delayed_surface_combustion_heat_megajoules, 0);
    @memset(context.surface_combustion_heat_megajoules_per_m2, 0);
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
                .volumetric_air_heat_capacity_megajoules_per_m3_k = context.runscript.ground_air_parameters.volumetric_air_heat_capacity_megajoules_per_m3_k,
            },
        };
        try tile_kernels.runKernelAcrossSerialTiles(context, &airflow_context, ecosys.canopy_airflow.applyTile);
        try surface_workspace.refresh(
            canopy.plant_canopy_aerodynamic_temperature_k,
            canopy.plant_canopy_aerodynamic_vapor_pressure_kpa,
            canopy.plant_minimum_water_vapor_resistance_h_per_m,
            canopy.plant_cuticular_water_vapor_resistance_h_per_m,
            context.plant_reproduction_controls.stomatal_turgor_shape,
            canopy.plant_canopy_turgor_potential_megapascal,
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
                .sensible_boundary_numerator_megajoules_per_m_h_k_by_cell = airflow.sensible_boundary_numerator_megajoules_per_m_h_k,
                .canopy_surface_temperature_k = context.plants.canopy_temperature_k,
                .aerodynamic_resistance_below_biome_h_per_m_by_cell = airflow.resistance_below_biome_h_per_m,
                .aerodynamic_resistance_below_species_h_per_m = airflow.resistance_below_species_h_per_m,
                .species_canopy_radiation_fraction = retention.living_radiation_fraction,
                .sensible_surface_resistance_h_per_m = surface_workspace.sensible_surface_resistance_h_per_m,
                .latent_surface_resistance_h_per_m = surface_workspace.latent_surface_resistance_h_per_m,
                .stomatal_resistance_h_per_m = surface_workspace.stomatal_resistance_h_per_m,
                .canopy_total_water_potential_megapascal = context.plants.canopy_water_potential_megapascal,
                .intercepted_water_volume_m3 = context.canopy_available_intercepted_water_m3,
            },
            .parameters = context.runscript.canopy_surface_exchange_parameters,
        };
        try tile_kernels.runKernelAcrossSerialTiles(context, &exchange_context, ecosys.canopy_surface_exchange.applyTile);
        @memcpy(canopy.plant_transpiration_m3_per_h, exchange.transpiration_m3_per_h);
        if (context.standing_dead_surface_exchange.*) |*dead_exchange| {
            for (0..context.grid.cell_count) |cell| for (0..context.config.plant_populations) |species| {
                const plant = cell * context.config.plant_populations + species;
                context.standing_dead_evaporation_m3_per_h[plant] = 0;
                dead_exchange.intercepted_water_change_m3_per_h[plant] = 0;
                dead_exchange.net_radiation_megajoules_per_h[plant] = 0;
                dead_exchange.sensible_heat_flux_megajoules_per_h[plant] = 0;
                dead_exchange.latent_heat_flux_megajoules_per_h[plant] = 0;
                dead_exchange.vapor_sensible_heat_flux_megajoules_per_h[plant] = 0;
                dead_exchange.storage_heat_flux_megajoules_per_h[plant] = 0;
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
                    .sensible_boundary_numerator_megajoules_per_m_h_k = airflow.sensible_boundary_numerator_megajoules_per_m_h_k[cell],
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
                    .latent_heat_of_vaporization_megajoules_per_m3 = context.runscript.canopy_surface_exchange_parameters.latent_heat_of_vaporization_megajoules_per_m3,
                    .liquid_water_heat_capacity_megajoules_per_m3_k = context.runscript.canopy_surface_exchange_parameters.liquid_water_heat_capacity_megajoules_per_m3_k,
                };
                const active_dry_volume_m3 = @min(
                    context.runscript.standing_dead_sapwood_thickness_m * retention.standing_dead_surface_area_m2[plant],
                    canopy.plant_standing_dead_carbon_g[plant] * context.runscript.stalk_volume_m3_per_g_c,
                );
                const dry_heat_capacity_megajoules_per_k = context.runscript.standing_dead_dry_volume_heat_capacity_megajoules_per_m3_k * active_dry_volume_m3;
                const activation_threshold_megajoules_per_k = context.runscript.standing_dead_activation_heat_capacity_megajoules_per_m2_k * context.canopy_cell_area_m2[cell];
                if (dry_heat_capacity_megajoules_per_k <= activation_threshold_megajoules_per_k) continue;
                const wet_heat_capacity_megajoules_per_k = dry_heat_capacity_megajoules_per_k +
                    dead_parameters.liquid_water_heat_capacity_megajoules_per_m3_k * retention.standing_dead_surface_water_m3[plant];
                const solved_surface = try ecosys.standing_dead_surface_exchange.solveSurfaceTemperature(.{
                    .exchange_inputs = exchange_inputs,
                    .absorbed_shortwave_megajoules_per_h = retention.standing_dead_absorbed_shortwave_megajoules_per_m2[plant] * context.canopy_cell_area_m2[cell],
                    .downward_longwave_megajoules_per_h = context.atmosphere.longwave_radiation_megajoules_per_m2[cell] * retention.standing_dead_radiation_fraction[plant] * context.canopy_cell_area_m2[cell],
                    .lateral_longwave_megajoules_per_h = 0,
                    .ground_surface_temperature_k = context.grid.surface_temperature_k[cell],
                    .emission_coefficient_megajoules_per_h_k4 = context.runscript.standing_dead_emissivity * 2.04e-10 * retention.standing_dead_radiation_fraction[plant] * context.canopy_cell_area_m2[cell],
                    .dry_and_existing_water_heat_capacity_megajoules_per_k = wet_heat_capacity_megajoules_per_k,
                    .retained_precipitation_water_m3_per_h = retention.standing_dead_retention_m3_per_h[plant],
                    .retained_precipitation_heat_megajoules_per_h = 0,
                    .minimum_effective_heat_capacity_megajoules_per_k = context.runscript.standing_dead_effective_heat_capacity_floor_megajoules_per_m2_k * context.canopy_cell_area_m2[cell],
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
                dead_exchange.net_radiation_megajoules_per_h[plant] = solved_surface.net_radiation_megajoules_per_h;
                dead_exchange.sensible_heat_flux_megajoules_per_h[plant] = result.sensible_heat_flux_megajoules_per_h;
                dead_exchange.latent_heat_flux_megajoules_per_h[plant] = result.latent_heat_flux_megajoules_per_h;
                dead_exchange.vapor_sensible_heat_flux_megajoules_per_h[plant] = result.vapor_sensible_heat_flux_megajoules_per_h;
                dead_exchange.storage_heat_flux_megajoules_per_h[plant] = solved_surface.storage_heat_flux_megajoules_per_h;
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
                        const atmospheric_sensible_conductance = airflow.sensible_boundary_numerator_megajoules_per_m_h_k[cell] * retention.standing_dead_radiation_fraction[plant] / result.total_aerodynamic_resistance_h_per_m * dead_share;
                        const atmospheric_vapor_conductance = @min(
                            airflow.latent_boundary_numerator_m2_per_h[cell] * retention.standing_dead_radiation_fraction[plant] / result.total_aerodynamic_resistance_h_per_m,
                            cell_air_volume_m3,
                        ) * dead_share;
                        const ground_sensible_conductance = airflow.sensible_boundary_numerator_megajoules_per_m_h_k[cell] / airflow.resistance_below_standing_dead_h_per_m[plant] * total_canopy_exposure * dead_share;
                        const ground_vapor_conductance = airflow.latent_boundary_numerator_m2_per_h[cell] / airflow.resistance_below_standing_dead_h_per_m[plant] * total_canopy_exposure * dead_share;
                        const dead_air_result = try ecosys.canopy_air_exchange.solveInto(dead_air, cell, species, .{
                            .initial_temperature_k = air_temperature_k,
                            .initial_vapor_fraction = canopy.plant_standing_dead_aerodynamic_vapor_pressure_kpa[plant] * context.runscript.canopy_surface_exchange_parameters.saturation_vapor_prefactor_k / air_temperature_k,
                            .atmospheric_temperature_k = context.atmosphere.air_temperature_k[cell],
                            .atmospheric_vapor_fraction = try ecosys.ground_air_exchange.vaporVolumeFraction(context.atmosphere.vapor_pressure_kpa[cell], context.atmosphere.air_temperature_k[cell], context.runscript.ground_air_parameters),
                            .ground_air_temperature_k = context.ground_air.temperature_k[cell],
                            .ground_air_vapor_fraction = context.ground_air.vapor_volume_fraction[cell],
                            .heat_capacity_megajoules_per_k = cell_air_volume_m3 * context.runscript.ground_air_parameters.volumetric_air_heat_capacity_megajoules_per_m3_k * dead_share,
                            .air_volume_m3 = cell_air_volume_m3 * dead_share,
                            .atmospheric_sensible_conductance_megajoules_per_h_k = atmospheric_sensible_conductance,
                            .atmospheric_vapor_conductance_m3_per_h = atmospheric_vapor_conductance,
                            .ground_sensible_conductance_megajoules_per_h_k = ground_sensible_conductance,
                            .ground_vapor_conductance_m3_per_h = ground_vapor_conductance,
                            .canopy_surface_sensible_heat_flux_megajoules_per_h = result.sensible_heat_flux_megajoules_per_h,
                            .canopy_surface_vapor_flux_m3_per_h = result.intercepted_water_change_m3_per_h,
                            .lateral_sensible_heat_flux_megajoules_per_h = context.delayed_standing_dead_combustion_heat_megajoules[plant],
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
                const cell_air_heat_capacity_megajoules_per_k = cell_air_volume_m3 * context.runscript.ground_air_parameters.volumetric_air_heat_capacity_megajoules_per_m3_k;
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
                    const atmospheric_sensible_conductance = airflow.sensible_boundary_numerator_megajoules_per_m_h_k[cell] * exposure_fraction / total_resistance_h_per_m * canopy_share;
                    const atmospheric_vapor_conductance = @min(
                        airflow.latent_boundary_numerator_m2_per_h[cell] * exposure_fraction / total_resistance_h_per_m,
                        cell_air_volume_m3,
                    ) * canopy_share;
                    const ground_sensible_conductance = airflow.sensible_boundary_numerator_megajoules_per_m_h_k[cell] / below_species_resistance_h_per_m * exposure.canopy_exposure_fraction[cell] * canopy_share;
                    const ground_vapor_conductance = airflow.latent_boundary_numerator_m2_per_h[cell] / below_species_resistance_h_per_m * exposure.canopy_exposure_fraction[cell] * canopy_share;
                    const result = try ecosys.canopy_air_exchange.solveInto(canopy_air, cell, species, .{
                        .initial_temperature_k = canopy.plant_canopy_aerodynamic_temperature_k[plant],
                        .initial_vapor_fraction = surface_workspace.canopy_air_vapor_fraction[plant],
                        .atmospheric_temperature_k = context.atmosphere.air_temperature_k[cell],
                        .atmospheric_vapor_fraction = try ecosys.ground_air_exchange.vaporVolumeFraction(context.atmosphere.vapor_pressure_kpa[cell], context.atmosphere.air_temperature_k[cell], context.runscript.ground_air_parameters),
                        .ground_air_temperature_k = context.ground_air.temperature_k[cell],
                        .ground_air_vapor_fraction = context.ground_air.vapor_volume_fraction[cell],
                        .heat_capacity_megajoules_per_k = cell_air_heat_capacity_megajoules_per_k * canopy_share,
                        .air_volume_m3 = cell_air_volume_m3 * canopy_share,
                        .atmospheric_sensible_conductance_megajoules_per_h_k = atmospheric_sensible_conductance,
                        .atmospheric_vapor_conductance_m3_per_h = atmospheric_vapor_conductance,
                        .ground_sensible_conductance_megajoules_per_h_k = ground_sensible_conductance,
                        .ground_vapor_conductance_m3_per_h = ground_vapor_conductance,
                        .canopy_surface_sensible_heat_flux_megajoules_per_h = exchange.sensible_heat_flux_megajoules_per_h[plant],
                        .canopy_surface_vapor_flux_m3_per_h = exchange.intercepted_water_change_m3_per_h[plant] + exchange.transpiration_m3_per_h[plant],
                        .lateral_sensible_heat_flux_megajoules_per_h = context.delayed_live_canopy_combustion_heat_megajoules[plant],
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
    try solveSoilHeatWaterAndSoluteTransport(
        context,
        hour_of_day,
        weather_header_by_cell,
        plant_calendar_by_cell,
        fertilizer_band_hour,
        gas_failure_report,
        solute_failure_report,
        diagnostic_first_hour,
        diagnostic_mineral_before_mol,
        diagnostic_relayer_phosphate_before,
        diagnostic_transport_ammonium_before,
        diagnostic_transport_before,
        plant_calendar,
        snow_phase_change_report,
        snow_vapor_equilibrium_report,
        &diagnostic_previous_heat_megajoules,
        &diagnostic_previous_n_g,
        &diagnostic_previous_p_g,
        &diagnostic_previous_p_owners,
    );
}

fn solveSoilHeatWaterAndSoluteTransport(
    context: anytype,
    hour_of_day: u8,
    weather_header_by_cell: []const ecosys.weather.Header,
    plant_calendar_by_cell: []const ecosys.plant_development.Calendar,
    fertilizer_band_hour: ecosys.fertilizer_band_phase_coordinator.HourToken,
    gas_failure_report: ?ecosys.soil_gas_transport_step.FailureReportRequest,
    solute_failure_report: ?ecosys.solute_failure_reporter.Request,
    diagnostic_first_hour: anytype,
    diagnostic_mineral_before_mol: anytype,
    diagnostic_relayer_phosphate_before: anytype,
    diagnostic_transport_ammonium_before: anytype,
    diagnostic_transport_before: anytype,
    plant_calendar: anytype,
    snow_phase_change_report: anytype,
    snow_vapor_equilibrium_report: anytype,
    diagnostic_previous_heat_megajoules_ptr: anytype,
    diagnostic_previous_n_g_ptr: anytype,
    diagnostic_previous_p_g_ptr: anytype,
    diagnostic_previous_p_owners_ptr: anytype,
) !void {
    var diagnostic_previous_heat_megajoules = diagnostic_previous_heat_megajoules_ptr.*;
    defer diagnostic_previous_heat_megajoules_ptr.* = diagnostic_previous_heat_megajoules;
    var diagnostic_previous_n_g = diagnostic_previous_n_g_ptr.*;
    defer diagnostic_previous_n_g_ptr.* = diagnostic_previous_n_g;
    var diagnostic_previous_p_g = diagnostic_previous_p_g_ptr.*;
    defer diagnostic_previous_p_g_ptr.* = diagnostic_previous_p_g;
    var diagnostic_previous_p_owners = diagnostic_previous_p_owners_ptr.*;
    defer diagnostic_previous_p_owners_ptr.* = diagnostic_previous_p_owners;
    for (0..context.grid.cell_count) |cell| {
        context.atmospheric_vapor_fraction[cell] = try ecosys.ground_air_exchange.vaporVolumeFraction(context.atmosphere.vapor_pressure_kpa[cell], context.atmosphere.air_temperature_k[cell], context.runscript.ground_air_parameters);
        context.ground_air_canopy_resistance_h_per_m[cell] = if (context.canopy_airflow.*) |*airflow| airflow.neutral_resistance_below_biome_h_per_m[cell] else 0;
    }
    @memset(context.delayed_live_canopy_combustion_heat_megajoules, 0);
    @memset(context.delayed_standing_dead_combustion_heat_megajoules, 0);
    @memset(context.ground_air_sensible_source_megajoules_per_h, 0);
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
                const sensible_conductance = airflow.sensible_boundary_numerator_megajoules_per_m_h_k[cell] / below_species_resistance_h_per_m * total_canopy_exposure * canopy_share;
                const vapor_conductance = airflow.latent_boundary_numerator_m2_per_h[cell] / below_species_resistance_h_per_m * total_canopy_exposure * canopy_share;
                context.ground_air_sensible_source_megajoules_per_h[cell] += sensible_conductance * (canopy_air.temperature_k[plant] - context.ground_air.temperature_k[cell]);
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
                const sensible_conductance = airflow.sensible_boundary_numerator_megajoules_per_m_h_k[cell] / resistance_h_per_m * total_canopy_exposure * dead_share;
                const vapor_conductance = airflow.latent_boundary_numerator_m2_per_h[cell] / resistance_h_per_m * total_canopy_exposure * dead_share;
                context.ground_air_sensible_source_megajoules_per_h[cell] += sensible_conductance * (dead_air.temperature_k[plant] - context.ground_air.temperature_k[cell]);
                context.ground_air_vapor_source_m3_per_h[cell] += vapor_conductance * (dead_air.vapor_fraction[plant] - context.ground_air.vapor_volume_fraction[cell]);
            }
        }
    };
    for (0..context.grid.cell_count) |cell| {
        context.ground_air_surface_sensible_conductance_megajoules_per_h_k[cell] = context.runscript.surface_sensible_heat_conductance_megajoules_per_m2_h_k * context.canopy_cell_area_m2[cell];
        context.ground_air_surface_vapor_conductance_m3_per_h[cell] = context.runscript.surface_latent_heat_conductance_megajoules_per_m2_h_kpa * context.canopy_cell_area_m2[cell] / context.runscript.ground_air_parameters.liquid_water_latent_heat_megajoules_per_m3 * context.ground_air.temperature_k[cell] / context.runscript.ground_air_parameters.saturation_vapor_prefactor_k;
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
        .non_atmospheric_sensible_heat_megajoules_per_h = context.ground_air_sensible_source_megajoules_per_h,
        .non_atmospheric_vapor_flux_m3_per_h = context.ground_air_vapor_source_m3_per_h,
        .non_atmospheric_sensible_conductance_megajoules_per_h_k = context.ground_air_surface_sensible_conductance_megajoules_per_h_k,
        .non_atmospheric_sensible_source_temperature_k = context.grid.surface_temperature_k,
        .non_atmospheric_vapor_conductance_m3_per_h = context.ground_air_surface_vapor_conductance_m3_per_h,
        .non_atmospheric_vapor_source_fraction = context.ground_air_surface_vapor_fraction,
    }, context.runscript.ground_air_parameters, .{
        .absolute_tolerance = context.config.absolute_tolerance,
        .relative_tolerance = context.config.relative_tolerance,
        .max_iterations = context.iteration_limits.water_heat_solute_max_iterations,
        .picard_relaxation = context.config.picard_relaxation,
    });
    try ecosys.ground_surface_vapor_water_commit.commit(.{
        .surface_vapor_conductance_m3_per_h = context.ground_air_surface_vapor_conductance_m3_per_h,
        .surface_vapor_fraction = context.ground_air_surface_vapor_fraction,
        .accepted_ground_air_vapor_fraction = context.ground_air.vapor_volume_fraction,
        .litter_liquid_water_m3 = context.surface_precipitation.litter_water_m3,
        .soil_matrix_liquid_water_m3 = context.grid.matrix_liquid_water_m3,
        .active_soil_layer_count = context.grid.active_soil_layer_count,
        .soil_layer_capacity = context.grid.soil_layer_capacity,
        .evaporation_m3_per_h = context.ground_surface_evaporation_m3_per_h,
        .condensation_m3_per_h = context.ground_surface_condensation_m3_per_h,
        .litter_liquid_water_change_m3 = context.ground_surface_litter_water_change_m3,
        .topsoil_liquid_water_change_m3 = context.ground_surface_topsoil_water_change_m3,
    });
    try ecosys.surface_litter_chemistry_carrier_rebase.rebaseFromAcceptedLiquidWaterChange(
        context.surface_litter_chemistry,
        context.surface_precipitation.litter_water_m3,
        context.ground_surface_litter_water_change_m3,
    );
    for (0..context.grid.cell_count) |cell| {
        const topsoil = try context.grid.layerIndex(cell, 0);
        const new_water_m3 = context.grid.matrix_liquid_water_m3[topsoil];
        try ecosys.soil_chemistry_water_carrier_rebase.rebaseLayer(
            context.soil_chemistry,
            topsoil,
            new_water_m3 - context.ground_surface_topsoil_water_change_m3[cell],
            new_water_m3,
        );
    }
    try ecosys.surface_precipitation.commitSoilIngress(context.surface_precipitation, context.grid, context.transport_hydrology, 1);
    try tile_kernels.runKernelAcrossSerialTiles(context, context.soil_thermal_context, ecosys.soil_thermal.updateTile);
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
    try context.soil_hourly_workspace.fillMacroporeFaceConductance(context.soil_transport_faces, context.soil_face_geometry, context.soil_face_geometry.macropore_hydraulic_conductance_m_per_h_megapascal);
    try context.soil_hourly_workspace.bindSurfaceHeatFlux(context.grid, context.surface_temperature);
    try ecosys.surface_precipitation.bindSoilHeatIngress(context.surface_precipitation, context.grid, context.soil_hourly_workspace.cell_heat_source_megajoules, 1);
    _ = try ecosys.subsurface_irrigation_heat.addToLayerHeatSources(
        context.soil_hourly_workspace.cell_heat_source_megajoules,
        context.subsurface_irrigation_water_m3,
        context.atmosphere.air_temperature_k,
        context.grid.soil_layer_capacity,
        context.runscript.canopy_surface_exchange_parameters.liquid_water_heat_capacity_megajoules_per_m3_k,
    );
    for (context.delayed_subsurface_combustion_heat_megajoules, context.soil_hourly_workspace.cell_heat_source_megajoules) |*delayed_heat, *heat_source| {
        heat_source.* += delayed_heat.*;
        delayed_heat.* = 0;
        if (!std.math.isFinite(heat_source.*)) return error.NonFiniteSubsurfaceCombustionHeatSource;
    }
    const water_heat_solute_max_iterations = try ecosys.iteration_control.waterHeatSoluteCeilingForCurrentState(context.iteration_limits.water_heat_solute_max_iterations, context.soil_hourly_workspace.heat_capacity_megajoules_per_k, context.soil_hourly_workspace.horizontal_face_area_m2, context.soil_hourly_workspace.is_top_soil_layer);
    var accepted_soil_water_heat = try ecosys.soil_water_heat_step.advanceMappedDeferred(context.allocator, context.grid, context.transport_hydrology, context.soil_transport_faces, context.soil_face_geometry, context.soil_solver_properties, context.soil_hourly_workspace, context.soil_thermal, context.soil_heat_solver_workspace, context.runscript.soil_phase_heat_parameters, .{ .max_iterations = water_heat_solute_max_iterations, .picard_relaxation = context.config.picard_relaxation, .vapor_pore_tortuosity = context.runscript.soil_process_parameters.vapor_pore_tortuosity, .osmotic_reflection_coefficient = context.runscript.soil_process_parameters.osmotic_reflection_coefficient, .absolute_tolerance = context.config.absolute_tolerance, .relative_tolerance = context.config.relative_tolerance, .boundary_topology = context.soil_boundary_topology, .geothermal_enabled_by_cell = context.geothermal_enabled_by_cell, .mean_annual_temperature_k_by_cell = context.mean_annual_temperature_k_by_cell, .geothermal_minimum_source_depth_m = context.runscript.geothermal_controls.minimum_source_depth_m, .geothermal_source_depth_below_profile_m = context.runscript.geothermal_controls.source_depth_below_profile_m, .geothermal_conductivity_m_megajoules_per_h_k = context.runscript.geothermal_controls.conductivity_m_megajoules_per_h_k, .geothermal_flux_megajoules_per_m2_h = context.runscript.geothermal_controls.geothermal_flux_megajoules_per_m2_h, .water_table_air_fraction_threshold = context.runscript.water_table_air_fraction_threshold, .active_layer_ice_fraction_threshold = context.runscript.active_layer_ice_fraction_threshold, .dense_newton_max_components = context.config.tile_cells, .matrix_external_water_source_m3_per_step = context.subsurface_irrigation_water_m3 });
    defer accepted_soil_water_heat.deinit();
    try commitHourlyWaterHeatStateGeneration(
        context,
        &accepted_soil_water_heat,
    );
    try context.landscape_boundary_ledger.accumulateAccepted(.{
        .heat_input_megajoules = accepted_soil_water_heat.solver.heat.boundary_heat_input_megajoules,
        .heat_output_megajoules = accepted_soil_water_heat.solver.heat.boundary_heat_output_megajoules,
    });
    if (diagnostic_first_hour) {
        const heat_megajoules = (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules;
        std.log.debug("heat stage: soil_water_heat_commit hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, heat_megajoules - diagnostic_previous_heat_megajoules });
        // EXEC-002: hours 5--16 inject a smoothly decaying POSITIVE error with
        // this stage dominant, which means the storage this stage moves and the
        // boundary heat it books disagree systematically. Report both, plus the
        // solver's own convergence state, so the discrepancy can be attributed.
        std.log.debug("soil heat boundary detail: hour={d} stage_delta_megajoules={e} boundary_input_megajoules={e} boundary_output_megajoules={e} net_boundary_megajoules={e} iterations={d} max_scaled_residual={e}", .{
            context.executed_weather_hours.* + 1,
            heat_megajoules - diagnostic_previous_heat_megajoules,
            accepted_soil_water_heat.solver.heat.boundary_heat_input_megajoules,
            accepted_soil_water_heat.solver.heat.boundary_heat_output_megajoules,
            accepted_soil_water_heat.solver.heat.boundary_heat_input_megajoules -
                accepted_soil_water_heat.solver.heat.boundary_heat_output_megajoules,
            accepted_soil_water_heat.solver.heat.iterations,
            accepted_soil_water_heat.solver.heat.maximum_scaled_residual,
        });
        // EXEC-002: the surface solver's conductive flux is an INTERNAL
        // surface<->soil transfer, so it is correctly absent from the boundary
        // ledger, but only if the soil actually receives what the surface gives.
        // Report both halves of the pairing.
        var diagnostic_surface_conduction_megajoules: f64 = 0;
        var diagnostic_soil_heat_source_megajoules: f64 = 0;
        for (0..context.grid.cell_count) |cell| {
            const top = cell * context.grid.soil_layer_capacity;
            // HEAT-001: `bindSurfaceHeatFlux` multiplies by
            // `horizontal_face_area_m2[top]`, so the instrument must use the
            // same area. Using `canopy_cell_area_m2` compared two different
            // areas and reported the difference as a mismatch.
            diagnostic_surface_conduction_megajoules +=
                context.surface_temperature.conductive_heat_flux_megajoules_per_m2[cell] *
                context.soil_hourly_workspace.horizontal_face_area_m2[top];
            // HEAT-001: read the published conduction slot, not
            // `cell_heat_source_megajoules[top]`, which by this point also
            // carries precipitation ingress, subsurface irrigation and
            // delayed combustion heat and so fabricates a mismatch.
            diagnostic_soil_heat_source_megajoules +=
                context.soil_hourly_workspace.published_surface_conduction_heat_megajoules[top];
        }
        std.log.debug("surface soil conduction pairing: hour={d} surface_gives_megajoules={e} soil_receives_megajoules={e} mismatch_megajoules={e}", .{ context.executed_weather_hours.* + 1, -diagnostic_surface_conduction_megajoules, diagnostic_soil_heat_source_megajoules, -diagnostic_surface_conduction_megajoules - diagnostic_soil_heat_source_megajoules });
        diagnostic_previous_heat_megajoules = heat_megajoules;
    }
    if (context.executed_weather_hours.* < 24) std.log.debug("current soil heat boundary: input_megajoules={e} output_megajoules={e}", .{ accepted_soil_water_heat.solver.heat.boundary_heat_input_megajoules, accepted_soil_water_heat.solver.heat.boundary_heat_output_megajoules });
    if (context.executed_weather_hours.* < 24) {
        if (accepted_soil_water_heat.solver.energy) |energy| std.log.debug(
            "soil transaction energy: richards_megajoules={e} vapor_megajoules={e} phase_sensible_megajoules={e} phase_absolute_megajoules={e} spatial_heat_megajoules={e} vapor_latent_megajoules={e} freeze_thaw_latent_megajoules={e}",
            .{
                energy.richards_enthalpy_change_megajoules,
                energy.vapor_transport_enthalpy_change_megajoules,
                energy.phase_enthalpy_change_megajoules,
                energy.phase_absolute_enthalpy_change_megajoules,
                energy.spatial_heat_enthalpy_change_megajoules,
                energy.phase_latent_heat_megajoules,
                energy.heat_solver_freeze_thaw_latent_megajoules,
            },
        );
    }
    // DO 225 / DO 245: redistribute soil pool contents across layer boundaries
    // that moved due to freeze-thaw volume change (DVOLI) this hour.
    {
        const carrier_count = ecosys.soil_water_heat_step.deferred_grid_carrier_count;
        const ws = context.soil_profile_relayering_workspace;
        // SOLUTE stores persistent minerals per unit liquid-water carrier.
        // Richards/freeze-thaw changes that carrier before chemistry runs, so
        // rebase concentrations to preserve extensive solid amounts.
        for (0..context.grid.layer_count) |l| {
            ws.ice_volume_delta_m3[l] =
                accepted_soil_water_heat.grid_delta_by_layer_carrier[l * carrier_count];
        }
        try ecosys.solute_solid_carrier_rebase.rebaseFromAcceptedChange(
            context.soil_chemistry,
            context.grid.matrix_liquid_water_m3,
            ws.ice_volume_delta_m3,
        );
        for (0..context.grid.layer_count) |l| {
            ws.ice_volume_delta_m3[l] =
                accepted_soil_water_heat.grid_delta_by_layer_carrier[l * carrier_count + 7];
        }
        for (0..context.grid.layer_count) |l| {
            const lv = context.soil_solver_properties.layer_volume_m3[l];
            ws.matrix_zone_fraction[l] = if (lv > 0)
                @max(1e-6, @min(1.0, context.soil_solver_properties.matrix_bulk_volume_m3[l] / lv))
            else
                1e-6;
        }
        try ecosys.soil_profile_relayering.assembleFreezeThawChanges(
            ws,
            context.soil_geometry,
            ws.ice_volume_delta_m3,
            ws.matrix_zone_fraction,
            context.canopy_cell_area_m2,
            context.runscript.soil_geometry_parameters.ice_to_water_specific_volume_difference,
            context.config.absolute_tolerance,
        );
        // The cached `total_heat_capacity_megajoules_per_m3_k` table was last
        // built at the top of the hour (the `soil_thermal.updateTile` pass
        // above the workspace refresh). Everything between here and there --
        // the water-and-heat commit, surface ingress, the phase solver -- moves
        // liquid, ice, and vapor between layers without rebuilding it, so by
        // this point it is stale.
        //
        // `applyLayerRedistribution` is the one consumer that cannot tolerate
        // that. It forms each layer's "before" energy from the *cached* table
        // (`heat_layer_remap.zig:112`) but rebuilds the "after"
        // capacity from the *live* carriers (:159), then divides to recover
        // temperature (:204). When the two disagree by a factor r, the layer's
        // temperature is restated as r*T even at zero transferred fraction.
        // Measured on `_split1` day 5, r ranged over 0.851--1.177, which drove
        // the top two layers +42.7 K and +45.2 K in a single hour with every
        // mass carrier bit-identical, and left layer 0 at 303 K while holding
        // 1.5e5 m3 of ice. The landscape census prices that ice at its full
        // frozen enthalpy and reads the inconsistency as created heat: this is
        // the whole of the `profile_relayering` stage delta.
        //
        // Refreshing here rather than reconstructing the capacity inside the
        // remap is deliberate. A local reconstruction was tried before and
        // stagnated the phase solver, because it left every other consumer of
        // the table reading the stale value; refreshing makes all of them agree.
        try tile_kernels.runKernelAcrossSerialTiles(context, context.soil_thermal_context, ecosys.soil_thermal.updateTile);
        try ecosys.soil_profile_relayering.applyLayerRedistribution(.{
            .grid = context.grid,
            .soil_thermal = context.soil_thermal,
            .gas_transport = context.gas_transport,
            .soil_organic = context.soil_organic,
            .soil_chemistry = context.soil_chemistry,
            .reactive_nitrogen = context.soil_reactive_nitrogen,
            .soil_fertilizer_inventory = context.soil_fertilizer_inventory,
            .mineral_fertilizer_inventory = context.mineral_fertilizer_inventory,
            .soil_properties = context.soil_solver_properties,
            .plant_roots = if (context.plant_roots.*) |*roots| roots else null,
            .soil_geometry = context.soil_geometry,
            .soil_face_geometry = context.soil_face_geometry,
            .soil_transport_faces = context.soil_transport_faces,
            .water_heat_parameters = .{
                .liquid_water_heat_capacity_megajoules_per_m3_k = context.runscript.soil_phase_heat_parameters.liquid_water_heat_capacity_megajoules_per_m3_k,
                .ice_heat_capacity_megajoules_per_m3_k = context.runscript.soil_phase_heat_parameters.ice_heat_capacity_megajoules_per_m3_k,
                .minimum_heat_capacity_megajoules_per_k = 0,
            },
            .dynamic_salts = context.runscript.dynamic_plant_salts,
            .nutrient_zone_fractions = .{
                .ammonium_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                .ammonium_band = context.runscript.plant_nutrient_initialization.initial_ammonium_band_fraction,
                .nitrate_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                .nitrate_band = context.runscript.plant_nutrient_initialization.initial_nitrate_band_fraction,
                .phosphate_non_band = 1 - context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
                .phosphate_band = context.runscript.plant_nutrient_initialization.initial_phosphate_band_fraction,
            },
            .plant_populations = context.config.plant_populations,
            .minimum_layer_thickness_m = context.runscript.soil_geometry_parameters.minimum_layer_thickness_m,
            .horizontal_cell_width_m = context.horizontal_cell_width_m,
            .vertical_cell_width_m = context.vertical_cell_width_m,
            // EXEC-002: `rebase_thermal_volume_to_geometry` keeps the
            // soil-thermal layer volume equal to the committed geometry and does
            // reduce relayering's spurious heat gain (3.455e8 to 2.942e8 MJ over
            // the first day), but measured end to end it makes the day-1 audit
            // deviation worse, from 7.075e-2 to 2.006e-1 per m2, because it also
            // shifts the water/heat commit and surface solve trajectory. It
            // stays off until that interaction is understood.
            .rebase_thermal_volume_to_geometry = false,
        }, ws.changes());
    }
    if (diagnostic_first_hour) {
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        std.log.debug("phosphorus stage: after_relayer delta_g={e}", .{current_p_g - diagnostic_previous_p_g});
        const owners = try diagnostics.diagnosticPhosphorusOwners_g(context);
        std.log.debug("phosphorus relayer owners: residue_g={e} organic_g={e} phosphate_g={e}", .{ owners[0] - diagnostic_previous_p_owners[0], owners[1] - diagnostic_previous_p_owners[1], owners[2] - diagnostic_previous_p_owners[2] });
        const phosphate_owners = diagnostics.diagnosticRelayerPhosphateOwners_g(context);
        std.log.debug("phosphorus relayer chemistry: dissolved_g={e} adsorbed_g={e} precipitate_g={e}", .{ phosphate_owners[0] - diagnostic_relayer_phosphate_before[0], phosphate_owners[1] - diagnostic_relayer_phosphate_before[1], phosphate_owners[2] - diagnostic_relayer_phosphate_before[2] });
        diagnostic_previous_p_owners = owners;
        diagnostic_previous_p_g = current_p_g;
        const heat_megajoules = (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules;
        std.log.debug("heat stage: profile_relayering hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, heat_megajoules - diagnostic_previous_heat_megajoules });
        diagnostic_previous_heat_megajoules = heat_megajoules;
    }
    const subsurface_irrigation_chemistry_parameters: ecosys.subsurface_irrigation_chemistry.Parameters = .{
        .molar_mass_g_per_mol = .{
            .nitrogen = context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol,
            .phosphorus = context.runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol,
            .aluminum = context.runscript.chemistry_primary_initialization.molar_mass_g_per_mol.aluminum,
            .iron = context.runscript.chemistry_primary_initialization.molar_mass_g_per_mol.iron,
            .calcium = context.runscript.chemistry_primary_initialization.molar_mass_g_per_mol.calcium,
            .magnesium = context.runscript.chemistry_primary_initialization.molar_mass_g_per_mol.magnesium,
            .sodium = context.runscript.chemistry_primary_initialization.molar_mass_g_per_mol.sodium,
            .potassium = context.runscript.chemistry_primary_initialization.molar_mass_g_per_mol.potassium,
            .sulfur = context.runscript.chemistry_primary_initialization.molar_mass_g_per_mol.sulfur,
            .chloride = context.runscript.chemistry_primary_initialization.molar_mass_g_per_mol.chloride,
        },
        .equilibrium = .{
            .aqueous = context.chemistry_reaction_parameters.*.aqueous_constants,
            .phosphate = context.chemistry_reaction_parameters.*.phosphate_constants,
        },
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
    if (diagnostic_first_hour) {
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        var boundary_p_g: f64 = 0;
        for (context.soil_solute_boundary_net_flux_mol, 0..) |flux_mol, component| {
            const species: ecosys.solute_transport_species.AqueousSpecies = @enumFromInt(component % ecosys.solute_transport_species.AqueousSpecies.count);
            if (ecosys.solute_transport_species.diffusivityClass(species) == .phosphate)
                boundary_p_g += flux_mol * context.runscript.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol;
        }
        std.log.debug("phosphorus stage: aqueous_transport delta_g={e} boundary_net_input_g={e} residual_g={e}", .{ current_p_g - diagnostic_previous_p_g, boundary_p_g, current_p_g - diagnostic_previous_p_g - boundary_p_g });
        diagnostic_previous_p_g = current_p_g;
    }
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
            .max_iterations = context.iteration_limits.organic_transport_max_iterations,
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
    if (diagnostic_first_hour) {
        const current_n_g = try diagnostics.diagnosticStoredNitrogen_g(context);
        std.log.debug("nitrogen stage: through_mineral_transport delta_g={e}", .{current_n_g - diagnostic_previous_n_g});
        const current_p_g = try diagnostics.diagnosticStoredPhosphorus_g(context);
        const diagnostic_boundary_p_g: f64 = 0;
        std.log.debug("phosphorus stage: through_mineral_transport delta_g={e} boundary_g={e} residual_g={e}", .{ current_p_g - diagnostic_previous_p_g, diagnostic_boundary_p_g, current_p_g - diagnostic_previous_p_g - diagnostic_boundary_p_g });
        diagnostic_previous_p_g = current_p_g;
        const heat_megajoules = (try diagnostics.reconstructLandscapeMassBalance(context)).heat_storage_megajoules;
        std.log.debug("heat stage: mineral_transport hour={d} delta_megajoules={e}", .{ context.executed_weather_hours.* + 1, heat_megajoules - diagnostic_previous_heat_megajoules });
        diagnostic_previous_heat_megajoules = heat_megajoules;
        var diagnostic_mineral_after_mol: f64 = 0;
        for (context.mineral_nitrogen_transport.matrix.amount_mol) |amount| diagnostic_mineral_after_mol += amount;
        for (context.mineral_nitrogen_transport.macropore.amount_mol) |amount| diagnostic_mineral_after_mol += amount;
        var diagnostic_mineral_export_g_n: f64 = 0;
        for (context.mineral_nitrogen_transport.boundary_export_g_n_per_step) |amount| diagnostic_mineral_export_g_n += amount;
        std.log.debug("mineral nitrogen transaction: inventory_delta_g_n={e} boundary_export_g_n={e} residual_g_n={e}", .{ (diagnostic_mineral_after_mol - diagnostic_mineral_before_mol) * context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol, diagnostic_mineral_export_g_n, (diagnostic_mineral_after_mol - diagnostic_mineral_before_mol) * context.runscript.fertilizer_nitrogen_molar_mass_g_per_mol + diagnostic_mineral_export_g_n });
        const components_n = try diagnostics.reconstructLandscapeMassBalance(context);
        std.log.debug("transport nitrogen components: residue={e} organic={e} n2={e} nh4={e} no3={e}", .{ components_n.residue_nitrogen_g - diagnostic_transport_before.residue_nitrogen_g, components_n.organic_nitrogen_g - diagnostic_transport_before.organic_nitrogen_g, components_n.dinitrogen_nitrogen_g - diagnostic_transport_before.dinitrogen_nitrogen_g, components_n.ammonium_nitrogen_g - diagnostic_transport_before.ammonium_nitrogen_g, components_n.nitrate_nitrogen_g - diagnostic_transport_before.nitrate_nitrogen_g });
        std.log.debug("transport phosphorus components: residue={e} organic={e} phosphate={e}", .{ components_n.residue_phosphorus_g - diagnostic_transport_before.residue_phosphorus_g, components_n.organic_phosphorus_g - diagnostic_transport_before.organic_phosphorus_g, components_n.phosphate_phosphorus_g - diagnostic_transport_before.phosphate_phosphorus_g });
        const owners = try diagnostics.diagnosticAmmoniumOwners_g_n(context);
        std.log.debug("transport ammonium owners delta: surface_aq={e} surface_exchange={e} surface_fertilizer={e} soil_aq={e} soil_exchange={e} soil_fertilizer={e}", .{ owners[0] - diagnostic_transport_ammonium_before[0], owners[1] - diagnostic_transport_ammonium_before[1], owners[2] - diagnostic_transport_ammonium_before[2], owners[3] - diagnostic_transport_ammonium_before[3], owners[4] - diagnostic_transport_ammonium_before[4], owners[5] - diagnostic_transport_ammonium_before[5] });
        diagnostic_previous_n_g = current_n_g;
    }
    try transportDissolvedGasAndSurfaceWater(
        context,
        hour_of_day,
        weather_header_by_cell,
        plant_calendar_by_cell,
        fertilizer_band_hour,
        gas_failure_report,
        solute_failure_report,
        diagnostic_first_hour,
        plant_calendar,
        snow_phase_change_report,
        snow_vapor_equilibrium_report,
        subsurface_irrigation_chemistry_parameters,
        &diagnostic_previous_heat_megajoules,
        &diagnostic_previous_n_g,
        &diagnostic_previous_p_g,
    );
}
