const std = @import("std");
const grid_module = @import("grid.zig");
const hydrology_module = @import("transport_hydrology.zig");
const solute = @import("solute_transport.zig");
const water_solver = @import("soil_water_solver.zig");
const vapor_solver = @import("soil_vapor_solver.zig");
const phase_solver = @import("soil_phase_solver.zig");
const heat_solver = @import("soil_heat_solver.zig");
const face_geometry_module = @import("soil_face_geometry.zig");
const solver_properties_module = @import("soil_solver_properties.zig");
const workspace_module = @import("soil_hourly_workspace.zig");
const thermal_module = @import("soil_thermal.zig");
const science_module = @import("soil_process_science.zig");
const boundary_topology_module = @import("soil_boundary_topology.zig");

pub const Inputs = struct {
    water_geometry: water_solver.FaceGeometry,
    water_properties: water_solver.Properties,
    water_options: water_solver.Options,
    vapor_geometry: vapor_solver.FaceGeometry,
    vapor_properties: vapor_solver.Properties,
    vapor_options: vapor_solver.Options,
    phase_properties: phase_solver.Properties,
    phase_options: phase_solver.Options,
    heat_geometry: heat_solver.FaceGeometry,
    heat_properties: heat_solver.Properties,
    heat_options: heat_solver.Options,
    heat_workspace: ?*heat_solver.Workspace = null,
};

pub const Result = struct { water: water_solver.Result, vapor: vapor_solver.Result, phase: phase_solver.Result, heat: heat_solver.Result };

pub const deferred_grid_carrier_count: usize = 12;

pub const DeferredMappedResult = struct {
    allocator: std.mem.Allocator,
    solver: Result,
    grid_delta_by_layer_carrier: []f64,

    pub fn deinit(self: *DeferredMappedResult) void {
        self.allocator.free(self.grid_delta_by_layer_carrier);
        self.* = undefined;
    }
};

/// Converges the complete coupled water/vapor/freeze-thaw/heat transaction,
/// captures its accepted grid-state delta, then restores authoritative state.
/// Flux and boundary ledgers remain published for downstream solute transport;
/// grid state is committed later from serial Morton tile files.
pub fn advanceMappedDeferred(
    allocator: std.mem.Allocator,
    grid: *grid_module.GridState,
    hydrology: *hydrology_module.State,
    faces: *hydrology_module.SoilFaces,
    geometry: *const face_geometry_module.State,
    properties: *const solver_properties_module.State,
    workspace: *const workspace_module.State,
    thermal: *const thermal_module.State,
    heat_workspace: *heat_solver.Workspace,
    science: science_module.RuntimeParameters,
    options: MappedOptions,
) !DeferredMappedResult {
    var base = try Snapshot.capture(allocator, grid, hydrology, faces);
    defer base.deinit();
    errdefer base.restore(grid, hydrology, faces);
    const solver = try advanceMapped(
        allocator,
        grid,
        hydrology,
        faces,
        geometry,
        properties,
        workspace,
        thermal,
        heat_workspace,
        science,
        options,
    );
    var accepted = try Snapshot.capture(allocator, grid, hydrology, faces);
    defer accepted.deinit();
    const delta = try allocator.alloc(
        f64,
        try std.math.mul(
            usize,
            grid.layer_count,
            deferred_grid_carrier_count,
        ),
    );
    errdefer allocator.free(delta);
    for (0..grid.layer_count) |layer| {
        for (0..deferred_grid_carrier_count) |carrier| {
            const value =
                accepted.gridStateSlice(carrier)[layer] -
                base.gridStateSlice(carrier)[layer];
            if (!std.math.isFinite(value))
                return error.NonFiniteDeferredSoilStateDelta;
            delta[layer * deferred_grid_carrier_count + carrier] = value;
        }
    }
    base.restore(grid, hydrology, faces);
    accepted.publishLedgers(hydrology, faces);
    return .{
        .allocator = allocator,
        .solver = solver,
        .grid_delta_by_layer_carrier = delta,
    };
}

pub const MappedOptions = struct {
    max_iterations: u16,
    picard_relaxation: f64,
    vapor_pore_tortuosity: f64,
    osmotic_reflection_coefficient: f64,
    absolute_tolerance: f64,
    relative_tolerance: f64,
    boundary_topology: ?*boundary_topology_module.State = null,
    geothermal_enabled_by_cell: ?[]const bool = null,
    mean_annual_temperature_k_by_cell: []const f64 = &.{},
    geothermal_minimum_source_depth_m: f64 = 10,
    geothermal_source_depth_below_profile_m: f64 = 1,
    geothermal_conductivity_m_mj_per_h_k: f64 = 8.1e-3,
    geothermal_flux_mj_per_m2_h: f64 = 2.052e-4,
    water_table_air_fraction_threshold: f64 = 1.0e-3,
    active_layer_ice_fraction_threshold: f64 = 1.0e-6,
    dense_newton_max_components: usize,
    matrix_external_water_source_m3_per_step: []const f64 = &.{},
};

/// Builds zero-copy solver views over mapped runtime state, then performs one
/// atomic whole-hour WATSUB transaction. NPH is solely a convergence ceiling.
pub fn advanceMapped(
    allocator: std.mem.Allocator,
    grid: *grid_module.GridState,
    hydrology: *hydrology_module.State,
    faces: *hydrology_module.SoilFaces,
    geometry: *const face_geometry_module.State,
    properties: *const solver_properties_module.State,
    workspace: *const workspace_module.State,
    thermal: *const thermal_module.State,
    heat_workspace: *heat_solver.Workspace,
    science: science_module.RuntimeParameters,
    options: MappedOptions,
) !Result {
    if (options.max_iterations == 0 or !std.math.isFinite(options.picard_relaxation) or options.picard_relaxation <= 0 or options.picard_relaxation > 1 or !std.math.isFinite(options.vapor_pore_tortuosity) or options.vapor_pore_tortuosity < 0 or !std.math.isFinite(options.osmotic_reflection_coefficient) or !std.math.isFinite(options.absolute_tolerance) or options.absolute_tolerance <= 0 or !std.math.isFinite(options.relative_tolerance) or options.relative_tolerance <= 0) return error.InvalidMappedSoilStepOptions;
    try science_module.validate(science);
    if (options.boundary_topology) |topology| try topology.refreshInternalWaterTable(grid, properties.matrix_bulk_volume_m3, properties.porosity_fraction, properties.matrix_air_entry_water_fraction, properties.layer_thickness_m, properties.layer_midpoint_depth_m, properties.layer_bottom_depth_m, options.water_table_air_fraction_threshold, options.active_layer_ice_fraction_threshold);
    return advance(allocator, grid, hydrology, faces, .{
        .water_geometry = .{ .source_path_length_m = geometry.source_path_length_m, .destination_path_length_m = geometry.destination_path_length_m, .face_area_m2 = geometry.face_area_m2, .wetting_depth_factor = geometry.wetting_depth_factor, .macropore_hydraulic_conductance_m_per_h_mpa = geometry.macropore_hydraulic_conductance_m_per_h_mpa },
        .water_properties = .{ .matrix_bulk_volume_m3 = properties.matrix_bulk_volume_m3, .matrix_air_entry_water_fraction = properties.matrix_air_entry_water_fraction, .retention_curve = properties.retention_curve, .mualem_van_genuchten_parameters = properties.mualem_van_genuchten_parameters, .macropore_mualem_van_genuchten_parameters = workspace.macropore_mualem_van_genuchten_parameters, .macropore_spacing_m = workspace.macropore_spacing_m, .macropore_radius_m = workspace.macropore_radius_m, .dual_domain_exchange_enabled = workspace.macropore_exchange_enabled, .dual_domain_geometry_factor = workspace.dual_domain_geometry_factor, .dual_domain_scaling_coefficient = workspace.dual_domain_scaling_coefficient, .frozen_hydraulic_impedance_exponent = workspace.frozen_hydraulic_impedance_exponent, .gravitational_water_potential_mpa_per_m = workspace.gravitational_water_potential_mpa_per_m, .gravitational_potential_mpa = workspace.gravitational_potential_mpa, .osmotic_potential_mpa = workspace.osmotic_potential_mpa, .matrix_hydraulic_conductivity_m2_per_h_mpa = properties.matrix_hydraulic_conductivity_m2_per_h_mpa, .rainfall_conductivity_multiplier = properties.rainfall_conductivity_multiplier, .matrix_external_source_m3_per_step = options.matrix_external_water_source_m3_per_step, .hydraulic_conductivity_class_count = properties.hydraulic_conductivity_class_count, .vertical_thickness_m = properties.layer_thickness_m, .osmotic_potential_multiplier = options.osmotic_reflection_coefficient, .boundary_topology = options.boundary_topology, .boundary_face_area_m2 = workspace.horizontal_face_area_m2, .boundary_macropore_hydraulic_conductivity_m2_per_h_mpa = workspace.macropore_hydraulic_conductivity_m2_per_h_mpa, .boundary_layer_volume_m3 = properties.layer_volume_m3, .boundary_layer_midpoint_depth_m = properties.layer_midpoint_depth_m, .boundary_layer_bottom_depth_m = properties.layer_bottom_depth_m },
        .water_options = .{ .max_iterations = options.max_iterations, .absolute_tolerance_m3 = options.absolute_tolerance, .relative_tolerance = options.relative_tolerance, .picard_relaxation = options.picard_relaxation, .maximum_newton_fraction = 8.0 },
        .vapor_geometry = .{ .source_path_length_m = geometry.source_path_length_m, .destination_path_length_m = geometry.destination_path_length_m, .face_area_m2 = geometry.face_area_m2 },
        .vapor_properties = .{ .vapor_diffusivity_m2_per_h = workspace.vapor_diffusivity_m2_per_h, .air_fraction = workspace.air_fraction, .porosity_fraction = properties.porosity_fraction, .tortuosity = options.vapor_pore_tortuosity },
        .vapor_options = .{ .max_iterations = options.max_iterations, .absolute_tolerance_m3 = options.absolute_tolerance, .relative_tolerance = options.relative_tolerance, .picard_relaxation = options.picard_relaxation },
        .phase_properties = .{ .matrix_bulk_volume_m3 = properties.matrix_bulk_volume_m3, .retention_curve = properties.retention_curve, .mualem_van_genuchten_parameters = properties.mualem_van_genuchten_parameters, .macropore_mualem_van_genuchten_parameters = workspace.macropore_mualem_van_genuchten_parameters, .osmotic_potential_mpa = workspace.osmotic_potential_mpa, .saturation_water_potential_mpa = properties.saturation_water_potential_mpa, .heat_capacity_mj_per_k = workspace.heat_capacity_mj_per_k, .saturated_lateral_matrix_conductivity_m2_per_h_mpa = workspace.saturated_lateral_matrix_conductivity_m2_per_h_mpa, .face_area_m2 = workspace.horizontal_face_area_m2, .macropore_spacing_m = workspace.macropore_spacing_m, .macropore_radius_m = workspace.macropore_radius_m, .pore_exchange_enabled = &.{}, .vapor = science.vapor_equilibrium, .freeze_thaw = science.freeze_thaw, .gravitational_water_potential_mpa_per_m = workspace.gravitational_water_potential_mpa_per_m, .liquid_water_heat_capacity_mj_per_m3_k = science.liquid_water_heat_capacity_mj_per_m3_k, .ice_heat_capacity_mj_per_m3_k = science.ice_heat_capacity_mj_per_m3_k },
        .phase_options = .{ .max_iterations = options.max_iterations, .absolute_tolerance_m3 = options.absolute_tolerance, .relative_tolerance = options.relative_tolerance, .picard_relaxation = options.picard_relaxation },
        .heat_geometry = .{ .source_path_length_m = geometry.source_path_length_m, .destination_path_length_m = geometry.destination_path_length_m, .face_area_m2 = geometry.face_area_m2 },
        .heat_properties = .{ .heat_capacity_mj_per_k = workspace.heat_capacity_mj_per_k, .minimum_heat_capacity_mj_per_k = workspace.minimum_heat_capacity_mj_per_k, .bulk_density_megagrams_per_m3 = properties.bulk_density_megagrams_per_m3, .liquid_water_fraction = workspace.liquid_water_fraction, .ice_fraction = workspace.ice_fraction, .air_fraction = workspace.air_fraction, .fraction_of_pore_volume_air_filled = workspace.fraction_of_pore_volume_air_filled, .solid_conductivity_numerator_m_mj_per_h_k = thermal.solid_thermal_conductivity_numerator_m_mj_per_h_k, .solid_conductivity_denominator = thermal.solid_thermal_conductivity_denominator, .is_top_soil_layer = workspace.is_top_soil_layer, .top_snow_heat_capacity_mj_per_k = workspace.top_snow_heat_capacity_mj_per_k, .maximum_negligible_snow_heat_capacity_mj_per_k = workspace.maximum_negligible_snow_heat_capacity_mj_per_k, .snow_storage_heat_flux_mj = workspace.snow_storage_heat_flux_mj, .cell_heat_source_mj = workspace.cell_heat_source_mj, .liquid_water_heat_capacity_mj_per_m3_k = science.liquid_water_heat_capacity_mj_per_m3_k, .turbulence = science.heat_turbulence, .geothermal_boundary = if (options.geothermal_enabled_by_cell) |enabled_by_cell| .{ .topology = options.boundary_topology orelse return error.MissingGeothermalBoundaryTopology, .layer_bottom_depth_m = properties.layer_bottom_depth_m, .lower_face_area_m2 = workspace.horizontal_face_area_m2, .enabled_by_cell = enabled_by_cell, .mean_annual_temperature_k_by_cell = options.mean_annual_temperature_k_by_cell, .minimum_source_depth_m = options.geothermal_minimum_source_depth_m, .source_depth_below_profile_m = options.geothermal_source_depth_below_profile_m, .conductivity_m_mj_per_h_k = options.geothermal_conductivity_m_mj_per_h_k, .geothermal_flux_mj_per_m2_h = options.geothermal_flux_mj_per_m2_h } else null, .enthalpy_coupling = .{ .matrix_liquid_water_m3 = grid.matrix_liquid_water_m3, .matrix_ice_water_equivalent_m3 = grid.matrix_ice_water_m3, .porous_medium_volume_m3 = properties.matrix_bulk_volume_m3, .matrix_pore_capacity_m3 = grid.matrix_pore_capacity_m3, .mualem_van_genuchten = properties.mualem_van_genuchten_parameters, .gravitational_water_potential_mpa_per_m = workspace.gravitational_water_potential_mpa_per_m, .pure_water_melting_temperature_k = science.freeze_thaw.pure_water_freezing_temperature_k, .ice_water_equivalent_heat_capacity_mj_per_m3_k = science.ice_heat_capacity_mj_per_m3_k, .latent_heat_of_fusion_mj_per_m3 = science.freeze_thaw.latent_heat_of_fusion_mj_per_m3, .solver_options = .{ .max_iterations = options.max_iterations, .absolute_enthalpy_tolerance_mj = options.absolute_tolerance, .relative_enthalpy_tolerance = options.relative_tolerance }, .macropore_liquid_water_m3 = grid.macropore_liquid_water_m3, .macropore_ice_water_equivalent_m3 = grid.macropore_ice_water_m3, .macropore_porous_medium_volume_m3 = grid.macropore_pore_capacity_m3, .macropore_mualem_van_genuchten = workspace.macropore_mualem_van_genuchten_parameters } },
        .heat_options = .{ .max_iterations = options.max_iterations, .absolute_tolerance_k = options.absolute_tolerance, .relative_tolerance = options.relative_tolerance, .picard_relaxation = options.picard_relaxation, .maximum_newton_fraction = 8.0, .dense_newton_max_components = options.dense_newton_max_components },
        .heat_workspace = heat_workspace,
    });
}

/// One atomic WATSUB transport transaction. Each process converges directly;
/// there is no surrounding sub-hour full-model loop.
pub fn advance(allocator: std.mem.Allocator, grid: *grid_module.GridState, hydrology: *hydrology_module.State, faces: *hydrology_module.SoilFaces, inputs: Inputs) !Result {
    var snapshot = try Snapshot.capture(allocator, grid, hydrology, faces);
    defer snapshot.deinit();
    errdefer snapshot.restore(grid, hydrology, faces);
    const water = try water_solver.solveAndBindTransportFaces(allocator, grid, hydrology, faces, inputs.water_geometry, inputs.water_properties, inputs.water_options);
    try deriveExternalWaterFluxes(grid, hydrology, faces, snapshot.grid_matrix_water, snapshot.grid_macro_water);
    const vapor = try vapor_solver.solveAndBindTransportFaces(allocator, grid, hydrology, faces, inputs.vapor_geometry, inputs.vapor_properties, inputs.vapor_options);
    const latent_heat_mj = try allocator.alloc(f64, grid.layer_count);
    defer allocator.free(latent_heat_mj);
    const pore_exchange_m3 = try allocator.alloc(f64, grid.layer_count);
    defer allocator.free(pore_exchange_m3);
    const phase_result = try phase_solver.solve(allocator, grid, inputs.phase_properties, .{ .latent_heat_mj = latent_heat_mj, .macropore_to_matrix_water_m3 = pore_exchange_m3 }, inputs.phase_options);
    @memcpy(hydrology.micropore_water_volume_m3, grid.matrix_liquid_water_m3);
    @memcpy(hydrology.macropore_water_volume_m3, grid.macropore_liquid_water_m3);
    @memcpy(hydrology.matrix_air_volume_m3, grid.matrix_air_volume_m3);
    @memcpy(hydrology.macropore_air_volume_m3, grid.macropore_air_volume_m3);
    @memcpy(hydrology.air_volume_m3, grid.air_volume_m3);
    @memcpy(hydrology.water_vapor_volume_m3, grid.water_vapor_volume_m3);
    // The phase solver's sixth coordinate is already the WATSUB ENGY1/VHCP1
    // endpoint temperature and therefore already contains condensation and
    // fusion enthalpy. Passing `latent_heat_mj` again as a heat source would
    // apply both latent terms twice. The spatial heat solve starts from that
    // phase-consistent temperature and adds only non-phase sources,
    // conduction, and convective carrier heat.
    const heat_result = if (inputs.heat_workspace) |heat_workspace|
        try heat_solver.solveAndBindTransportFacesWithWorkspace(
            heat_workspace,
            grid,
            hydrology,
            faces,
            inputs.heat_geometry,
            inputs.heat_properties,
            inputs.heat_options,
        )
    else
        try heat_solver.solveAndBindTransportFaces(
            allocator,
            grid,
            hydrology,
            faces,
            inputs.heat_geometry,
            inputs.heat_properties,
            inputs.heat_options,
        );
    // The mapped heat residual performs the final Dall'Amico enthalpy
    // repartition at its accepted temperature. Publish that atomic phase
    // state instead of retaining the pre-conduction phase snapshot.
    @memcpy(hydrology.micropore_water_volume_m3, grid.matrix_liquid_water_m3);
    @memcpy(hydrology.macropore_water_volume_m3, grid.macropore_liquid_water_m3);
    @memcpy(hydrology.matrix_air_volume_m3, grid.matrix_air_volume_m3);
    @memcpy(hydrology.macropore_air_volume_m3, grid.macropore_air_volume_m3);
    @memcpy(hydrology.air_volume_m3, grid.air_volume_m3);
    try publishAcceptedIceVolumeChanges(hydrology, snapshot.grid_matrix_ice, snapshot.grid_macro_ice, grid.matrix_ice_water_m3, grid.macropore_ice_water_m3);
    return .{ .water = water, .vapor = vapor, .phase = phase_result, .heat = heat_result };
}

fn publishAcceptedIceVolumeChanges(hydrology: *hydrology_module.State, initial_matrix_m3: []const f64, initial_macropore_m3: []const f64, final_matrix_m3: []const f64, final_macropore_m3: []const f64) !void {
    const count = hydrology.matrix_ice_volume_change_m3_per_step.len;
    if (initial_matrix_m3.len != count or initial_macropore_m3.len != count or final_matrix_m3.len != count or final_macropore_m3.len != count or hydrology.macropore_ice_volume_change_m3_per_step.len != count or hydrology.total_ice_volume_change_m3_per_step.len != count) return error.AcceptedIceVolumeChangeDimensionMismatch;
    for (0..count) |layer| {
        const matrix_change_m3 = final_matrix_m3[layer] - initial_matrix_m3[layer];
        const macropore_change_m3 = final_macropore_m3[layer] - initial_macropore_m3[layer];
        const total_change_m3 = matrix_change_m3 + macropore_change_m3;
        if (!std.math.isFinite(matrix_change_m3) or !std.math.isFinite(macropore_change_m3) or !std.math.isFinite(total_change_m3)) return error.NonFiniteAcceptedIceVolumeChange;
    }
    for (0..count) |layer| {
        hydrology.matrix_ice_volume_change_m3_per_step[layer] = final_matrix_m3[layer] - initial_matrix_m3[layer];
        hydrology.macropore_ice_volume_change_m3_per_step[layer] = final_macropore_m3[layer] - initial_macropore_m3[layer];
        hydrology.total_ice_volume_change_m3_per_step[layer] = hydrology.matrix_ice_volume_change_m3_per_step[layer] + hydrology.macropore_ice_volume_change_m3_per_step[layer];
    }
}

const Snapshot = struct {
    allocator: std.mem.Allocator,
    grid_matrix_water: []f64,
    grid_macro_water: []f64,
    grid_total_water: []f64,
    grid_matrix_air: []f64,
    grid_macro_air: []f64,
    grid_total_air: []f64,
    grid_vapor: []f64,
    grid_matrix_ice: []f64,
    grid_macro_ice: []f64,
    grid_total_ice: []f64,
    grid_temperature: []f64,
    grid_matric_potential: []f64,
    hydrology_matrix_water: []f64,
    hydrology_macro_water: []f64,
    hydrology_matrix_air: []f64,
    hydrology_macro_air: []f64,
    hydrology_total_air: []f64,
    hydrology_vapor: []f64,
    hydrology_matrix_flux: []f64,
    hydrology_macro_flux: []f64,
    hydrology_vapor_flux: []f64,
    hydrology_heat_flux: []f64,
    hydrology_matrix_external_flux: []f64,
    hydrology_macro_external_flux: []f64,
    hydrology_matrix_ice_change: []f64,
    hydrology_macro_ice_change: []f64,
    hydrology_total_ice_change: []f64,
    face_matrix_flux: []f64,
    face_macro_flux: []f64,
    face_vapor_flux: []f64,
    face_heat_flux: []f64,
    micropore_faces: []solute.Face,
    macropore_faces: []solute.Face,

    fn capture(allocator: std.mem.Allocator, grid: *const grid_module.GridState, hydrology: *const hydrology_module.State, faces: *const hydrology_module.SoilFaces) !Snapshot {
        var result: Snapshot = undefined;
        result.allocator = allocator;
        var allocated: usize = 0;
        errdefer result.freeAllocated(allocated);
        inline for (@typeInfo(Snapshot).@"struct".fields) |field| {
            if (field.type == []f64) {
                const source: []const f64 = snapshotF64Source(field.name, grid, hydrology, faces);
                @field(result, field.name) = try allocator.dupe(f64, source);
                allocated += 1;
            } else if (field.type == []solute.Face) {
                const source: []const solute.Face = if (comptime std.mem.eql(u8, field.name, "micropore_faces")) faces.micropore_faces else faces.macropore_faces;
                @field(result, field.name) = try allocator.dupe(solute.Face, source);
                allocated += 1;
            }
        }
        return result;
    }

    fn restore(self: *const Snapshot, grid: *grid_module.GridState, hydrology: *hydrology_module.State, faces: *hydrology_module.SoilFaces) void {
        @memcpy(grid.matrix_liquid_water_m3, self.grid_matrix_water);
        @memcpy(grid.macropore_liquid_water_m3, self.grid_macro_water);
        @memcpy(grid.liquid_water_m3, self.grid_total_water);
        @memcpy(grid.matrix_air_volume_m3, self.grid_matrix_air);
        @memcpy(grid.macropore_air_volume_m3, self.grid_macro_air);
        @memcpy(grid.air_volume_m3, self.grid_total_air);
        @memcpy(grid.water_vapor_volume_m3, self.grid_vapor);
        @memcpy(grid.matrix_ice_water_m3, self.grid_matrix_ice);
        @memcpy(grid.macropore_ice_water_m3, self.grid_macro_ice);
        @memcpy(grid.ice_water_m3, self.grid_total_ice);
        @memcpy(grid.soil_temperature_k, self.grid_temperature);
        @memcpy(grid.matric_potential_mpa, self.grid_matric_potential);
        @memcpy(hydrology.micropore_water_volume_m3, self.hydrology_matrix_water);
        @memcpy(hydrology.macropore_water_volume_m3, self.hydrology_macro_water);
        @memcpy(hydrology.matrix_air_volume_m3, self.hydrology_matrix_air);
        @memcpy(hydrology.macropore_air_volume_m3, self.hydrology_macro_air);
        @memcpy(hydrology.air_volume_m3, self.hydrology_total_air);
        @memcpy(hydrology.water_vapor_volume_m3, self.hydrology_vapor);
        @memcpy(hydrology.micropore_face_flux_m3_per_step, self.hydrology_matrix_flux);
        @memcpy(hydrology.macropore_face_flux_m3_per_step, self.hydrology_macro_flux);
        @memcpy(hydrology.vapor_face_flux_m3_per_step, self.hydrology_vapor_flux);
        @memcpy(hydrology.heat_face_flux_mj_per_step, self.hydrology_heat_flux);
        @memcpy(hydrology.micropore_external_water_flux_m3_per_step, self.hydrology_matrix_external_flux);
        @memcpy(hydrology.macropore_external_water_flux_m3_per_step, self.hydrology_macro_external_flux);
        @memcpy(hydrology.matrix_ice_volume_change_m3_per_step, self.hydrology_matrix_ice_change);
        @memcpy(hydrology.macropore_ice_volume_change_m3_per_step, self.hydrology_macro_ice_change);
        @memcpy(hydrology.total_ice_volume_change_m3_per_step, self.hydrology_total_ice_change);
        @memcpy(faces.micropore_water_flux_m3_per_step, self.face_matrix_flux);
        @memcpy(faces.macropore_water_flux_m3_per_step, self.face_macro_flux);
        @memcpy(faces.vapor_flux_m3_per_step, self.face_vapor_flux);
        @memcpy(faces.heat_flux_mj_per_step, self.face_heat_flux);
        @memcpy(faces.micropore_faces, self.micropore_faces);
        @memcpy(faces.macropore_faces, self.macropore_faces);
    }

    fn publishLedgers(
        self: *const Snapshot,
        hydrology: *hydrology_module.State,
        faces: *hydrology_module.SoilFaces,
    ) void {
        @memcpy(hydrology.micropore_face_flux_m3_per_step, self.hydrology_matrix_flux);
        @memcpy(hydrology.macropore_face_flux_m3_per_step, self.hydrology_macro_flux);
        @memcpy(hydrology.vapor_face_flux_m3_per_step, self.hydrology_vapor_flux);
        @memcpy(hydrology.heat_face_flux_mj_per_step, self.hydrology_heat_flux);
        @memcpy(hydrology.micropore_external_water_flux_m3_per_step, self.hydrology_matrix_external_flux);
        @memcpy(hydrology.macropore_external_water_flux_m3_per_step, self.hydrology_macro_external_flux);
        @memcpy(hydrology.matrix_ice_volume_change_m3_per_step, self.hydrology_matrix_ice_change);
        @memcpy(hydrology.macropore_ice_volume_change_m3_per_step, self.hydrology_macro_ice_change);
        @memcpy(hydrology.total_ice_volume_change_m3_per_step, self.hydrology_total_ice_change);
        @memcpy(faces.micropore_water_flux_m3_per_step, self.face_matrix_flux);
        @memcpy(faces.macropore_water_flux_m3_per_step, self.face_macro_flux);
        @memcpy(faces.vapor_flux_m3_per_step, self.face_vapor_flux);
        @memcpy(faces.heat_flux_mj_per_step, self.face_heat_flux);
        @memcpy(faces.micropore_faces, self.micropore_faces);
        @memcpy(faces.macropore_faces, self.macropore_faces);
    }

    fn gridStateSlice(
        self: *const Snapshot,
        carrier: usize,
    ) []const f64 {
        return switch (carrier) {
            0 => self.grid_matrix_water,
            1 => self.grid_macro_water,
            2 => self.grid_total_water,
            3 => self.grid_matrix_air,
            4 => self.grid_macro_air,
            5 => self.grid_total_air,
            6 => self.grid_vapor,
            7 => self.grid_matrix_ice,
            8 => self.grid_macro_ice,
            9 => self.grid_total_ice,
            10 => self.grid_temperature,
            11 => self.grid_matric_potential,
            else => unreachable,
        };
    }

    fn deinit(self: *Snapshot) void {
        inline for (@typeInfo(Snapshot).@"struct".fields) |field| if (field.type == []f64 or field.type == []solute.Face) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }

    fn freeAllocated(self: *Snapshot, count: usize) void {
        var visited: usize = 0;
        inline for (@typeInfo(Snapshot).@"struct".fields) |field| if (field.type == []f64 or field.type == []solute.Face) {
            if (visited < count) self.allocator.free(@field(self, field.name));
            visited += 1;
        };
    }
};

fn snapshotF64Source(comptime name: []const u8, grid: *const grid_module.GridState, hydrology: *const hydrology_module.State, faces: *const hydrology_module.SoilFaces) []const f64 {
    if (comptime std.mem.eql(u8, name, "grid_matrix_water")) return grid.matrix_liquid_water_m3;
    if (comptime std.mem.eql(u8, name, "grid_macro_water")) return grid.macropore_liquid_water_m3;
    if (comptime std.mem.eql(u8, name, "grid_total_water")) return grid.liquid_water_m3;
    if (comptime std.mem.eql(u8, name, "grid_matrix_air")) return grid.matrix_air_volume_m3;
    if (comptime std.mem.eql(u8, name, "grid_macro_air")) return grid.macropore_air_volume_m3;
    if (comptime std.mem.eql(u8, name, "grid_total_air")) return grid.air_volume_m3;
    if (comptime std.mem.eql(u8, name, "grid_vapor")) return grid.water_vapor_volume_m3;
    if (comptime std.mem.eql(u8, name, "grid_matrix_ice")) return grid.matrix_ice_water_m3;
    if (comptime std.mem.eql(u8, name, "grid_macro_ice")) return grid.macropore_ice_water_m3;
    if (comptime std.mem.eql(u8, name, "grid_total_ice")) return grid.ice_water_m3;
    if (comptime std.mem.eql(u8, name, "grid_temperature")) return grid.soil_temperature_k;
    if (comptime std.mem.eql(u8, name, "grid_matric_potential")) return grid.matric_potential_mpa;
    if (comptime std.mem.eql(u8, name, "hydrology_matrix_water")) return hydrology.micropore_water_volume_m3;
    if (comptime std.mem.eql(u8, name, "hydrology_macro_water")) return hydrology.macropore_water_volume_m3;
    if (comptime std.mem.eql(u8, name, "hydrology_matrix_air")) return hydrology.matrix_air_volume_m3;
    if (comptime std.mem.eql(u8, name, "hydrology_macro_air")) return hydrology.macropore_air_volume_m3;
    if (comptime std.mem.eql(u8, name, "hydrology_total_air")) return hydrology.air_volume_m3;
    if (comptime std.mem.eql(u8, name, "hydrology_vapor")) return hydrology.water_vapor_volume_m3;
    if (comptime std.mem.eql(u8, name, "hydrology_matrix_flux")) return hydrology.micropore_face_flux_m3_per_step;
    if (comptime std.mem.eql(u8, name, "hydrology_macro_flux")) return hydrology.macropore_face_flux_m3_per_step;
    if (comptime std.mem.eql(u8, name, "hydrology_vapor_flux")) return hydrology.vapor_face_flux_m3_per_step;
    if (comptime std.mem.eql(u8, name, "hydrology_heat_flux")) return hydrology.heat_face_flux_mj_per_step;
    if (comptime std.mem.eql(u8, name, "hydrology_matrix_external_flux")) return hydrology.micropore_external_water_flux_m3_per_step;
    if (comptime std.mem.eql(u8, name, "hydrology_macro_external_flux")) return hydrology.macropore_external_water_flux_m3_per_step;
    if (comptime std.mem.eql(u8, name, "hydrology_matrix_ice_change")) return hydrology.matrix_ice_volume_change_m3_per_step;
    if (comptime std.mem.eql(u8, name, "hydrology_macro_ice_change")) return hydrology.macropore_ice_volume_change_m3_per_step;
    if (comptime std.mem.eql(u8, name, "hydrology_total_ice_change")) return hydrology.total_ice_volume_change_m3_per_step;
    if (comptime std.mem.eql(u8, name, "face_matrix_flux")) return faces.micropore_water_flux_m3_per_step;
    if (comptime std.mem.eql(u8, name, "face_macro_flux")) return faces.macropore_water_flux_m3_per_step;
    if (comptime std.mem.eql(u8, name, "face_vapor_flux")) return faces.vapor_flux_m3_per_step;
    if (comptime std.mem.eql(u8, name, "face_heat_flux")) return faces.heat_flux_mj_per_step;
    unreachable;
}

/// Recovers the external WATSUB contribution from its conservative water
/// balance before freeze/thaw and matrix/macropore phase exchange alter the
/// storages. Positive output is discharge, matching TRNSFRS boundary signs.
fn deriveExternalWaterFluxes(grid: *const grid_module.GridState, hydrology: *hydrology_module.State, faces: *const hydrology_module.SoilFaces, initial_matrix_water_m3: []const f64, initial_macropore_water_m3: []const f64) !void {
    if (initial_matrix_water_m3.len != grid.layer_count or initial_macropore_water_m3.len != grid.layer_count or hydrology.micropore_external_water_flux_m3_per_step.len != grid.layer_count or hydrology.macropore_external_water_flux_m3_per_step.len != grid.layer_count) return error.ExternalWaterFluxDimensionMismatch;
    @memset(hydrology.micropore_external_water_flux_m3_per_step, 0);
    @memset(hydrology.macropore_external_water_flux_m3_per_step, 0);
    for (faces.micropore_faces, faces.macropore_faces) |micro_face, macro_face| {
        const micro_flux = micro_face.water_flux_m3_per_step;
        const macro_flux = macro_face.water_flux_m3_per_step;
        hydrology.micropore_external_water_flux_m3_per_step[micro_face.first_cell] -= micro_flux;
        hydrology.micropore_external_water_flux_m3_per_step[micro_face.second_cell] += micro_flux;
        hydrology.macropore_external_water_flux_m3_per_step[macro_face.first_cell] -= macro_flux;
        hydrology.macropore_external_water_flux_m3_per_step[macro_face.second_cell] += macro_flux;
    }
    for (0..grid.layer_count) |layer| {
        const internal_matrix_change = hydrology.micropore_external_water_flux_m3_per_step[layer];
        const internal_macro_change = hydrology.macropore_external_water_flux_m3_per_step[layer];
        hydrology.micropore_external_water_flux_m3_per_step[layer] = internal_matrix_change - (grid.matrix_liquid_water_m3[layer] - initial_matrix_water_m3[layer]);
        hydrology.macropore_external_water_flux_m3_per_step[layer] = internal_macro_change - (grid.macropore_liquid_water_m3[layer] - initial_macropore_water_m3[layer]);
        if (!std.math.isFinite(hydrology.micropore_external_water_flux_m3_per_step[layer]) or !std.math.isFinite(hydrology.macropore_external_water_flux_m3_per_step[layer])) return error.NonFiniteExternalWaterFlux;
    }
}

test "external water flux is separated from conservative internal movement" {
    const cfg = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 2, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    @memset(grid.active_soil_layer_count, 1);
    var snow = try @import("snow_solute_transport.zig").State.init(std.testing.allocator, 2, 1);
    defer snow.deinit();
    var hydrology = try hydrology_module.State.init(std.testing.allocator, 2, 1, 1, 1);
    defer hydrology.deinit();
    var faces = try hydrology_module.buildSoilFaces(std.testing.allocator, &hydrology, &grid);
    defer faces.deinit();
    faces.micropore_faces[0].water_flux_m3_per_step = 0.1;
    faces.macropore_faces[0].water_flux_m3_per_step = 0.02;
    grid.matrix_liquid_water_m3[0] = 0.85;
    grid.matrix_liquid_water_m3[1] = 1.1;
    grid.macropore_liquid_water_m3[0] = 0.98;
    grid.macropore_liquid_water_m3[1] = 1.01;
    try deriveExternalWaterFluxes(&grid, &hydrology, &faces, &.{ 1, 1 }, &.{ 1, 1 });
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), hydrology.micropore_external_water_flux_m3_per_step[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0), hydrology.micropore_external_water_flux_m3_per_step[1], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0), hydrology.macropore_external_water_flux_m3_per_step[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), hydrology.macropore_external_water_flux_m3_per_step[1], 1e-14);
}

test "accepted WATSUB ice ledger publishes REDIST DVOLI atomically" {
    var hydrology = try hydrology_module.State.init(std.testing.allocator, 1, 1, 2, 1);
    defer hydrology.deinit();
    try publishAcceptedIceVolumeChanges(&hydrology, &.{ 0.1, 0.3 }, &.{ 0.2, 0.1 }, &.{ 0.25, 0.2 }, &.{ 0.15, 0.4 });
    try std.testing.expectApproxEqAbs(@as(f64, 0.15), hydrology.matrix_ice_volume_change_m3_per_step[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, -0.1), hydrology.matrix_ice_volume_change_m3_per_step[1], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, -0.05), hydrology.macropore_ice_volume_change_m3_per_step[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), hydrology.macropore_ice_volume_change_m3_per_step[1], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), hydrology.total_ice_volume_change_m3_per_step[0], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), hydrology.total_ice_volume_change_m3_per_step[1], 1e-14);
    hydrology.total_ice_volume_change_m3_per_step[0] = 9;
    try std.testing.expectError(error.NonFiniteAcceptedIceVolumeChange, publishAcceptedIceVolumeChanges(&hydrology, &.{ 0, 0 }, &.{ 0, 0 }, &.{ std.math.nan(f64), 0 }, &.{ 0, 0 }));
    try std.testing.expectEqual(@as(f64, 9), hydrology.total_ice_volume_change_m3_per_step[0]);
}

test "late WATSUB failure rolls back water vapor heat and shared faces" {
    const retention = @import("soil_water_retention.zig");
    const cfg = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 2, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.matrix_liquid_water_m3[0] = 0.35;
    grid.matrix_liquid_water_m3[1] = 0.15;
    @memset(grid.matrix_pore_capacity_m3, 0.5);
    grid.matrix_air_volume_m3[0] = 0.15;
    grid.matrix_air_volume_m3[1] = 0.35;
    @memcpy(grid.air_volume_m3, grid.matrix_air_volume_m3);
    grid.water_vapor_volume_m3[0] = 0.01;
    grid.soil_temperature_k[0] = 300;
    grid.soil_temperature_k[1] = 280;
    var snow = try @import("snow_solute_transport.zig").State.init(std.testing.allocator, 2, 1);
    defer snow.deinit();
    var hydrology = try hydrology_module.State.init(std.testing.allocator, 2, 1, 1, 1);
    defer hydrology.deinit();
    try hydrology.syncStorage(&grid, &snow);
    var faces = try hydrology_module.buildSoilFaces(std.testing.allocator, &hydrology, &grid);
    defer faces.deinit();
    const curve: retention.ResolvedCurve = .{ .porosity_fraction = 0.5, .curve = .{ .field_capacity_fraction = 0.3, .wilting_point_fraction = 0.1, .saturation_water_potential_mpa = -0.0005, .field_capacity_water_potential_mpa = -0.01, .wilting_point_water_potential_mpa = -1.5, .minimum_water_potential_mpa = -1.5e12, .saturation_to_field_shape = 0.5, .below_wilting_shape = 0.5 } };
    const curves = [_]retention.ResolvedCurve{ curve, curve };
    const one_face = [_]f64{1};
    const zero_face = [_]f64{0};
    const one_cell = [_]f64{ 1, 1 };
    const phase_heat_capacity = [_]f64{ 2, 2 };
    const half_cell = [_]f64{ 0.5, 0.5 };
    const zero_cell = [_]f64{ 0, 0 };
    const conductivity = [_]f64{0.002} ** 6;
    const thickness = [_]f64{ 0.1, 0.1 };
    const bools = [_]bool{ true, false };
    const bad_capacity = [_]f64{ 0, 0 };
    const saturation_potential = [_]f64{ -0.0005, -0.0005 };
    const pore_exchange_disabled = [_]bool{ false, false };
    const matrix_parameters = [_]retention.MualemVanGenuchtenParameters{
        .{ .residual_water_content_m3_per_m3 = 0.05, .saturated_water_content_m3_per_m3 = 0.5, .alpha_per_m = 3.6, .n = 1.56, .saturated_hydraulic_conductivity_m_per_h = 0.002 },
        .{ .residual_water_content_m3_per_m3 = 0.05, .saturated_water_content_m3_per_m3 = 0.5, .alpha_per_m = 3.6, .n = 1.56, .saturated_hydraulic_conductivity_m_per_h = 0.002 },
    };
    const macropore_parameters = [_]retention.MualemVanGenuchtenParameters{
        .{ .residual_water_content_m3_per_m3 = 0, .saturated_water_content_m3_per_m3 = 1, .alpha_per_m = 15, .n = 2.68, .saturated_hydraulic_conductivity_m_per_h = 0.1 },
        .{ .residual_water_content_m3_per_m3 = 0, .saturated_water_content_m3_per_m3 = 1, .alpha_per_m = 15, .n = 2.68, .saturated_hydraulic_conductivity_m_per_h = 0.1 },
    };
    const macropore_spacing_m = [_]f64{ 0.2, 0.2 };
    const macropore_radius_m = [_]f64{ 0.001, 0.001 };
    const before_water = grid.matrix_liquid_water_m3[0];
    try std.testing.expectError(error.InvalidSoilHeatCapacity, advance(std.testing.allocator, &grid, &hydrology, &faces, .{
        .water_geometry = .{ .source_path_length_m = &one_face, .destination_path_length_m = &one_face, .face_area_m2 = &one_face, .wetting_depth_factor = &zero_face, .macropore_hydraulic_conductance_m_per_h_mpa = &zero_face },
        .water_properties = .{ .matrix_bulk_volume_m3 = &one_cell, .matrix_air_entry_water_fraction = &half_cell, .retention_curve = &curves, .mualem_van_genuchten_parameters = &matrix_parameters, .macropore_mualem_van_genuchten_parameters = &macropore_parameters, .macropore_spacing_m = &macropore_spacing_m, .macropore_radius_m = &macropore_radius_m, .dual_domain_exchange_enabled = &pore_exchange_disabled, .gravitational_potential_mpa = &zero_cell, .osmotic_potential_mpa = &zero_cell, .matrix_hydraulic_conductivity_m2_per_h_mpa = &conductivity, .vertical_thickness_m = &thickness, .osmotic_potential_multiplier = 1 },
        .water_options = .{ .max_iterations = 20 },
        .vapor_geometry = .{ .source_path_length_m = &one_face, .destination_path_length_m = &one_face, .face_area_m2 = &one_face },
        .vapor_properties = .{ .vapor_diffusivity_m2_per_h = &one_cell, .air_fraction = &half_cell, .porosity_fraction = &half_cell, .tortuosity = 1 },
        .vapor_options = .{ .max_iterations = 20 },
        .phase_properties = .{ .matrix_bulk_volume_m3 = &one_cell, .retention_curve = &curves, .mualem_van_genuchten_parameters = &matrix_parameters, .macropore_mualem_van_genuchten_parameters = &macropore_parameters, .osmotic_potential_mpa = &zero_cell, .saturation_water_potential_mpa = &saturation_potential, .heat_capacity_mj_per_k = &phase_heat_capacity, .saturated_lateral_matrix_conductivity_m2_per_h_mpa = &one_cell, .face_area_m2 = &one_cell, .macropore_spacing_m = &one_cell, .macropore_radius_m = &one_cell, .pore_exchange_enabled = &pore_exchange_disabled, .vapor = .{ .vapor_density_temperature_coefficient = 2.173e-3, .molecular_weight_ratio = 0.61, .clausius_clapeyron_coefficient_k = 5360, .reference_inverse_temperature_per_k = 3.661e-3, .water_molar_mass_g_per_mol = 18, .gas_constant_j_per_mol_k = 8.3143, .latent_heat_of_vaporization_mj_per_m3 = 2450 }, .freeze_thaw = .{ .freezing_potential_numerator_k_mpa = 9.0959e4, .latent_heat_of_fusion_mj_per_m3 = 333, .ice_density_megagrams_per_m3 = 0.917, .heat_capacity_temperature_feedback_per_k = 6.2913e-3, .pure_water_freezing_temperature_k = 273.15 }, .liquid_water_heat_capacity_mj_per_m3_k = 4.19, .ice_heat_capacity_mj_per_m3_k = 1.9274 },
        .phase_options = .{ .max_iterations = 20 },
        .heat_geometry = .{ .source_path_length_m = &one_face, .destination_path_length_m = &one_face, .face_area_m2 = &one_face },
        .heat_properties = .{ .heat_capacity_mj_per_k = &bad_capacity, .minimum_heat_capacity_mj_per_k = &zero_cell, .bulk_density_megagrams_per_m3 = &one_cell, .liquid_water_fraction = &half_cell, .ice_fraction = &zero_cell, .air_fraction = &half_cell, .fraction_of_pore_volume_air_filled = &half_cell, .solid_conductivity_numerator_m_mj_per_h_k = &one_cell, .solid_conductivity_denominator = &one_cell, .is_top_soil_layer = &bools, .top_snow_heat_capacity_mj_per_k = &zero_cell, .maximum_negligible_snow_heat_capacity_mj_per_k = &zero_cell, .snow_storage_heat_flux_mj = &zero_cell, .cell_heat_source_mj = &zero_cell, .liquid_water_heat_capacity_mj_per_m3_k = 4.19, .turbulence = .{ .water_fraction_threshold = 1, .air_fraction_threshold = 1, .water_rayleigh_coefficient = 0, .air_rayleigh_coefficient = 0, .water_nusselt_denominator = 1, .air_nusselt_denominator = 1 } },
        .heat_options = .{ .max_iterations = 20 },
    }));
    try std.testing.expectEqual(before_water, grid.matrix_liquid_water_m3[0]);
    try std.testing.expectEqual(@as(f64, 0.01), grid.water_vapor_volume_m3[0]);
    try std.testing.expectEqual(@as(f64, 300), grid.soil_temperature_k[0]);
    try std.testing.expectEqual(@as(f64, 0), faces.micropore_water_flux_m3_per_step[0]);
    try std.testing.expectEqual(@as(f64, 0), faces.vapor_flux_m3_per_step[0]);
}
