const std = @import("std");
const delimited_input = @import("delimited_input.zig");
const soil_solver_properties = @import("soil_solver_properties.zig");
const soil_hourly_workspace = @import("soil_hourly_workspace.zig");
const soil_process_science = @import("soil_process_science.zig");
const soil_gas_transport_step = @import("soil_gas_transport_step.zig");
const plant_root_gas_exchange = @import("plant_root_gas_exchange.zig");
const plant_root_nutrient_uptake = @import("plant_root_nutrient_uptake.zig");
const plant_root_salt_exchange = @import("plant_root_salt_exchange.zig");
const plant_root_exudation = @import("plant_root_exudation.zig");
const plant_root_mycorrhizal_exchange = @import("plant_root_mycorrhizal_exchange.zig");
const plant_root_porosity = @import("plant_root_porosity.zig");
const plant_root_metabolism = @import("plant_root_metabolism.zig");
const plant_organ_partition = @import("plant_organ_partition.zig");
const shoot_growth_metabolism = @import("shoot_growth_metabolism.zig");
const shoot_growth_runtime = @import("shoot_growth_runtime.zig");
const plant_symbiotic_fixation = @import("plant_symbiotic_fixation.zig");
const plant_root_disturbance = @import("plant_root_disturbance.zig");
const soil_combustion = @import("soil_combustion.zig");
const plant_shoot_root_exchange = @import("plant_shoot_root_exchange.zig");
const plant_storage_remobilization = @import("plant_storage_remobilization.zig");
const plant_initialization = @import("plant_initialization.zig");
const plant_root_system = @import("plant_root_system.zig");
const plant_phenology = @import("plant_phenology.zig");
const plant_pool_aggregation = @import("plant_pool_aggregation.zig");
const canopy_photosynthesis = @import("canopy_photosynthesis.zig");
const canopy_biochemistry = @import("canopy_biochemistry.zig");
const canopy_optics = @import("canopy_optics.zig");
const canopy_precipitation_retention = @import("canopy_precipitation_retention.zig");
const soil_plant_available_nutrients = @import("soil_plant_available_nutrients.zig");
const soil_chemistry_initialization = @import("soil_chemistry_initialization.zig");
const snow_compaction = @import("snow_compaction.zig");
const snow_heat_conduction = @import("snow_heat_conduction.zig");
const snow_vapor_equilibrium = @import("snow_vapor_equilibrium.zig");
const snow_vapor_diffusion = @import("snow_vapor_diffusion.zig");
const surface_gas_boundary_conductance = @import("surface_gas_boundary_conductance.zig");
const surface_aerodynamics = @import("surface_aerodynamics.zig");
const surface_runoff = @import("surface_runoff.zig");
const surface_precipitation = @import("surface_precipitation.zig");
const ground_air_exchange = @import("ground_air_exchange.zig");
const canopy_surface_exchange = @import("canopy_surface_exchange.zig");
const plant_soil_exchange = @import("plant_soil_exchange.zig");
const spatial_grid = @import("spatial_grid.zig");

pub const Domain = struct {
    west_column: usize,
    north_row: usize,
    east_column: usize,
    south_row: usize,

    pub fn columns(self: Domain) !usize {
        if (self.east_column < self.west_column) return error.InvalidDomain;
        return try std.math.add(usize, self.east_column - self.west_column, 1);
    }

    pub fn rows(self: Domain) !usize {
        if (self.south_row < self.north_row) return error.InvalidDomain;
        return try std.math.add(usize, self.south_row - self.north_row, 1);
    }
};

pub const SceneFiles = struct {
    weather_grid_file: []const u8,
    options: []const u8,
    land_management: []const u8,
    plant_management: []const u8,
    output_editors: [10][]const u8,
};

pub const Scenario = struct {
    first_scene_index: usize,
    scene_count: usize,
    repeat_count: usize,
};

pub const GeothermalControls = struct {
    enabled: bool,
    minimum_source_depth_m: f64,
    source_depth_below_profile_m: f64,
    conductivity_m_mj_per_h_k: f64,
    geothermal_flux_mj_per_m2_h: f64,
};

pub const SoilGeometryParameters = struct {
    organic_carbon_specific_volume_m3_per_g: f64,
    organic_horizon_threshold_g_c_per_Mg: f64,
    ice_to_water_specific_volume_difference: f64,
    minimum_layer_thickness_m: f64,
};

pub const RunScript = struct {
    allocator: std.mem.Allocator,
    domain: Domain,
    geospatial_bounds: ?spatial_grid.BoundsDegrees,
    tile_row_count: usize,
    tile_column_count: usize,
    lateral_flow_halo_cell_count: usize,
    grid_input_file: []const u8,
    plant_species_count: usize,
    uses_four_value_species_default: bool,
    worker_count: usize,
    tile_cell_count: usize,
    relative_tolerance: f64,
    absolute_tolerance: f64,
    max_nonlinear_iterations: u16,
    picard_relaxation: f64,
    uses_compatibility_runtime_controls: bool,
    soil_solver_parameters: soil_solver_properties.RuntimeParameters,
    uses_compatibility_soil_solver_parameters: bool,
    soil_process_parameters: soil_hourly_workspace.RuntimeParameters,
    uses_compatibility_soil_process_parameters: bool,
    soil_gas_transport_parameters: soil_gas_transport_step.RuntimeParameters,
    uses_compatibility_soil_gas_transport_parameters: bool,
    soil_phase_heat_parameters: soil_process_science.RuntimeParameters,
    uses_compatibility_soil_phase_heat_parameters: bool,
    geothermal_controls: GeothermalControls,
    uses_compatibility_geothermal_controls: bool,
    water_table_air_fraction_threshold: f64,
    active_layer_ice_fraction_threshold: f64,
    uses_compatibility_subsurface_state_controls: bool,
    soil_geometry_parameters: SoilGeometryParameters,
    uses_compatibility_soil_geometry_parameters: bool,
    surface_pond_dry_organic_heat_capacity_mj_per_g_c_k: f64,
    surface_pond_activation_heat_capacity_mj_per_m2_k: f64,
    uses_compatibility_surface_pond_energy: bool,
    soil_longwave_emissivity: f64,
    snow_longwave_emissivity: f64,
    canopy_longwave_emissivity: f64,
    snow_full_cover_depth_m: f64,
    surface_sensible_heat_conductance_mj_per_m2_h_k: f64,
    surface_latent_heat_conductance_mj_per_m2_h_kpa: f64,
    surface_vapor_activity_fraction: f64,
    minimum_surface_temperature_k: f64,
    maximum_surface_temperature_k: f64,
    uses_compatibility_surface_energy: bool,
    snow_layer_bottom_depth_m: []f64,
    initial_snow_density_Mg_per_m3: f64,
    snow_ice_density_Mg_per_m3: f64,
    snow_latent_heat_of_fusion_mj_per_m3: f64,
    snow_phase_damping_divisor: f64,
    uses_compatibility_snow_layers: bool,
    snow_compaction_parameters: snow_compaction.Parameters,
    uses_compatibility_snow_compaction: bool,
    snow_heat_conduction_parameters: snow_heat_conduction.Parameters,
    uses_compatibility_snow_heat_conduction: bool,
    snow_vapor_parameters: snow_vapor_equilibrium.Parameters,
    uses_compatibility_snow_vapor: bool,
    snow_vapor_diffusion_parameters: snow_vapor_diffusion.Parameters,
    uses_compatibility_snow_vapor_diffusion: bool,
    surface_gas_resistance_parameters: surface_gas_boundary_conductance.Parameters,
    uses_compatibility_surface_gas_resistance: bool,
    surface_runoff_parameters: surface_runoff.Parameters,
    uses_compatibility_surface_runoff: bool,
    rainfall_impact_parameters: surface_precipitation.RainfallImpactParameters,
    uses_compatibility_rainfall_impact: bool,
    surface_aerodynamic_parameters: surface_aerodynamics.Parameters,
    uses_compatibility_surface_aerodynamics: bool,
    ground_air_parameters: ground_air_exchange.Parameters,
    uses_compatibility_ground_air: bool,
    canopy_surface_exchange_parameters: canopy_surface_exchange.Parameters,
    canopy_sensible_surface_resistance_h_per_m: f64,
    canopy_latent_surface_resistance_h_per_m: f64,
    uses_compatibility_canopy_surface_exchange: bool,
    canopy_ammonia_exchange_parameters: plant_soil_exchange.CanopyAmmoniaExchangeParameters,
    uses_compatibility_canopy_ammonia_exchange: bool,
    root_axes_per_plant: usize,
    uses_compatibility_plant_structure: bool,
    canopy_layer_count: usize,
    uses_compatibility_canopy_layers: bool,
    canopy_discretization: @import("canopy_geometry.zig").Discretization,
    uses_compatibility_canopy_discretization: bool,
    stalk_volume_m3_per_g_c: f64,
    uses_compatibility_canopy_geometry: bool,
    standing_dead_partition_parameters: plant_initialization.StandingDeadPartitionParameters,
    uses_compatibility_standing_dead_partition: bool,
    plant_heat_water_parameters: plant_initialization.PlantHeatWaterParameters,
    uses_compatibility_plant_heat_water: bool,
    plant_geometry_parameters: plant_initialization.PlantGeometryParameters,
    uses_compatibility_plant_geometry: bool,
    phenology_initialization_parameters: plant_initialization.PhenologyInitializationParameters,
    uses_compatibility_phenology_initialization: bool,
    root_initialization_parameters: plant_root_system.InitializationParameters,
    uses_compatibility_root_initialization: bool,
    root_morphology_parameters: plant_root_system.MorphologyParameters,
    uses_compatibility_root_morphology: bool,
    standing_dead_sapwood_thickness_m: f64,
    standing_dead_dry_volume_heat_capacity_mj_per_m3_k: f64,
    standing_dead_emissivity: f64,
    standing_dead_activation_heat_capacity_mj_per_m2_k: f64,
    standing_dead_effective_heat_capacity_floor_mj_per_m2_k: f64,
    uses_compatibility_standing_dead_energy: bool,
    woody_optics_parameters: canopy_optics.WoodyOpticsParameters,
    uses_compatibility_woody_optics: bool,
    canopy_retention_parameters: canopy_precipitation_retention.Parameters,
    uses_compatibility_canopy_retention: bool,
    shoot_control_parameters: plant_initialization.ShootControlParameters,
    uses_compatibility_shoot_control_parameters: bool,
    c4_carbon_parameters: canopy_photosynthesis.C4CarbonParameters,
    uses_source_c4_carbon_parameters: bool,
    thermal_acclimation_parameters: plant_initialization.ThermalAcclimationParameters,
    uses_compatibility_thermal_acclimation_parameters: bool,
    canopy_stress_parameters: canopy_biochemistry.StressParameters,
    uses_compatibility_canopy_stress_parameters: bool,
    phenology_parameters: plant_phenology.Parameters,
    uses_compatibility_phenology_parameters: bool,
    plant_pool_parameters: plant_pool_aggregation.RuntimeParameters,
    dynamic_plant_salts: bool,
    uses_compatibility_plant_pool_parameters: bool,
    seed_set_parameters: canopy_photosynthesis.SeedSetParameters,
    uses_compatibility_seed_set_parameters: bool,
    root_gas_parameters: plant_root_gas_exchange.RuntimeParameters,
    uses_compatibility_root_gas_parameters: bool,
    root_nutrient_parameters: plant_root_nutrient_uptake.RuntimeParameters,
    uses_compatibility_root_nutrient_parameters: bool,
    root_salt_parameters: plant_root_salt_exchange.Parameters,
    uses_compatibility_root_salt_parameters: bool,
    root_mycorrhizal_exchange_parameters: plant_root_mycorrhizal_exchange.Parameters,
    uses_compatibility_root_mycorrhizal_exchange_parameters: bool,
    root_exudation_parameters: plant_root_exudation.Parameters,
    uses_compatibility_root_exudation_parameters: bool,
    root_porosity_parameters: plant_root_porosity.Parameters,
    uses_compatibility_root_porosity_parameters: bool,
    root_metabolism_parameters: plant_root_metabolism.SecondaryRootParameters,
    uses_compatibility_root_metabolism_parameters: bool,
    organ_partition_parameters: plant_organ_partition.Parameters,
    uses_compatibility_organ_partition_parameters: bool,
    shoot_metabolism_parameters: shoot_growth_metabolism.Parameters,
    uses_compatibility_shoot_metabolism_parameters: bool,
    shoot_node_growth_parameters: shoot_growth_runtime.NodeGrowthParameters,
    uses_source_shoot_node_growth_parameters: bool,
    branch_mobile_exchange_parameters: shoot_growth_runtime.BranchMobileExchangeParameters,
    uses_source_branch_mobile_exchange_parameters: bool,
    symbiotic_fixation_parameters: plant_symbiotic_fixation.RuntimeParameters,
    uses_source_symbiotic_fixation_parameters: bool,
    plant_fire_combustion_parameters: plant_root_disturbance.CombustionParameters,
    uses_source_plant_fire_combustion_parameters: bool,
    soil_fire_combustion_parameters: soil_combustion.Parameters,
    uses_source_soil_fire_combustion_parameters: bool,
    shoot_root_exchange_parameters: plant_shoot_root_exchange.Parameters,
    uses_compatibility_shoot_root_exchange_parameters: bool,
    storage_remobilization_parameters: plant_storage_remobilization.Parameters,
    uses_compatibility_storage_remobilization_parameters: bool,
    plant_nutrient_initialization: soil_plant_available_nutrients.InitializationParameters,
    uses_compatibility_plant_nutrient_initialization: bool,
    microbial_substrate_count: usize,
    microbial_population_count: usize,
    uses_compatibility_microbial_dimensions: bool,
    organic_initialization_file: ?[]const u8,
    surface_gas_parameter_file: ?[]const u8,
    soil_nitrogen_parameter_file: ?[]const u8,
    chemistry_initialization: ?soil_chemistry_initialization.ProfileSolubleParameters,
    chemistry_primary_initialization: ?soil_chemistry_initialization.PrimaryInitializationParameters,
    chemistry_reaction_file: ?[]const u8,
    fertilizer_nitrogen_molar_mass_g_per_mol: f64,
    uses_compatibility_fertilizer_units: bool,
    execution_repeat_count: usize,
    scenarios: []Scenario,
    scenes: []SceneFiles,

    pub fn deinit(self: *RunScript) void {
        self.allocator.free(self.snow_layer_bottom_depth_m);
        if (self.organic_initialization_file) |name| self.allocator.free(name);
        if (self.surface_gas_parameter_file) |name| self.allocator.free(name);
        if (self.soil_nitrogen_parameter_file) |name| self.allocator.free(name);
        if (self.chemistry_reaction_file) |name| self.allocator.free(name);
        self.allocator.free(self.grid_input_file);
        for (self.scenes) |scene| {
            self.allocator.free(scene.weather_grid_file);
            self.allocator.free(scene.options);
            self.allocator.free(scene.land_management);
            self.allocator.free(scene.plant_management);
            for (scene.output_editors) |name| self.allocator.free(name);
        }
        self.allocator.free(self.scenes);
        self.allocator.free(self.scenarios);
        self.* = undefined;
    }

    pub fn gridCellCount(self: RunScript) !usize {
        return try std.math.mul(usize, try self.domain.columns(), try self.domain.rows());
    }
};

pub fn parse(allocator: std.mem.Allocator, source: []const u8) !RunScript {
    const body = extractInputBody(source);
    try rejectEmptyExplicitFields(body);
    const first_line_end = std.mem.indexOfScalar(u8, body, '\n') orelse return error.MissingRunscriptBody;
    try validatePlantNutrientsRecordArity(body[first_line_end + 1 ..]);
    var header_tokens = delimited_input.recordTokens(body[0..first_line_end]);
    var header_values: [5]usize = undefined;
    var header_count: usize = 0;
    while (header_tokens.next()) |token| {
        if (header_count == header_values.len) return error.TooManyDomainHeaderValues;
        header_values[header_count] = std.fmt.parseUnsigned(usize, token, 10) catch return error.InvalidRunscriptInteger;
        header_count += 1;
    }
    var domain: Domain = switch (header_count) {
        // Modern compact form: horizontal_cells, vertical_cells, plant_species.
        3 => .{ .west_column = 1, .north_row = 1, .east_column = header_values[0], .south_row = header_values[1] },
        // Four-value bounds form retained for the supplied examples.
        4, 5 => .{ .west_column = header_values[0], .north_row = header_values[1], .east_column = header_values[2], .south_row = header_values[3] },
        else => return error.InvalidDomainHeader,
    };
    const plant_species_count = if (header_count == 3) header_values[2] else if (header_count == 5) header_values[4] else 5;
    const uses_four_value_species_default = header_count == 4;
    if (plant_species_count == 0) return error.NoPlantSpecies;
    _ = try domain.columns();
    _ = try domain.rows();

    var tokens = delimited_input.tokens(body[first_line_end + 1 ..]);

    const first_after_header = try next(&tokens, "geospatial-grid, runtime, or grid-inputs record");
    const has_geospatial_grid = std.ascii.eqlIgnoreCase(first_after_header, "geospatial_grid");
    const geospatial_bounds: ?spatial_grid.BoundsDegrees = if (has_geospatial_grid) .{
        .minimum_latitude_degrees_north = try nextFloat(&tokens, "minimum latitude in degrees north"),
        .maximum_latitude_degrees_north = try nextFloat(&tokens, "maximum latitude in degrees north"),
        .minimum_longitude_degrees_east = try nextFloat(&tokens, "minimum longitude in degrees east"),
        .maximum_longitude_degrees_east = try nextFloat(&tokens, "maximum longitude in degrees east"),
        .latitude_interval_degrees = try nextFloat(&tokens, "latitude interval in degrees"),
        .longitude_interval_degrees = try nextFloat(&tokens, "longitude interval in degrees"),
    } else null;
    var tile_row_count: usize = 1;
    var tile_column_count: usize = 1;
    var lateral_flow_halo_cell_count: usize = 2;
    const first_after_geospatial = if (has_geospatial_grid)
        try next(&tokens, "tile-layout record")
    else
        first_after_header;
    const has_tile_layout = std.ascii.eqlIgnoreCase(first_after_geospatial, "tile_layout");
    if (has_geospatial_grid and !has_tile_layout) return error.MissingTileLayoutRecord;
    if (has_tile_layout) {
        tile_row_count = try nextUnsigned(&tokens, "tile row count");
        tile_column_count = try nextUnsigned(&tokens, "tile column count");
        lateral_flow_halo_cell_count = try nextUnsigned(&tokens, "lateral-flow halo cell count");
        if (tile_row_count == 0 or tile_column_count == 0 or lateral_flow_halo_cell_count != 2)
            return error.InvalidTileLayout;
    }
    if (geospatial_bounds) |bounds| {
        var inferred_grid = try spatial_grid.RegularGrid.init(allocator, bounds);
        defer inferred_grid.deinit();
        domain = .{
            .west_column = 1,
            .north_row = 1,
            .east_column = inferred_grid.column_count,
            .south_row = inferred_grid.row_count,
        };
    }
    const first_after_domain = if (has_tile_layout)
        try next(&tokens, "runtime or grid-inputs record")
    else
        first_after_geospatial;
    const has_runtime_record = std.ascii.eqlIgnoreCase(first_after_domain, "runtime");
    const worker_count = if (has_runtime_record) try nextUnsigned(&tokens, "worker count") else 1;
    const tile_cell_count = if (has_runtime_record) try nextUnsigned(&tokens, "tile cell count") else 256;
    const relative_tolerance = if (has_runtime_record) try nextFloat(&tokens, "relative tolerance") else 1.0e-8;
    const absolute_tolerance = if (has_runtime_record) try nextFloat(&tokens, "absolute tolerance") else 1.0e-11;
    const max_nonlinear_iterations = if (has_runtime_record) try nextUnsignedType(u16, &tokens, "maximum nonlinear iterations") else 40;
    const picard_relaxation = if (has_runtime_record) try nextFloat(&tokens, "Picard relaxation") else 0.5;
    if (worker_count == 0 or tile_cell_count == 0 or max_nonlinear_iterations == 0 or !std.math.isFinite(picard_relaxation) or picard_relaxation <= 0 or picard_relaxation > 1) return error.InvalidRuntimeControls;
    if (!std.math.isFinite(relative_tolerance) or relative_tolerance <= 0.0 or
        !std.math.isFinite(absolute_tolerance) or absolute_tolerance <= 0.0) return error.InvalidRuntimeControls;
    const first_after_runtime = if (has_runtime_record) try next(&tokens, "soil-solver, surface-energy, or site record") else first_after_domain;
    const has_soil_solver_record = std.ascii.eqlIgnoreCase(first_after_runtime, "soil_solver");
    if (has_geospatial_grid and !has_soil_solver_record)
        return error.MissingSoilSolverRecord;
    const soil_solver_parameters = if (has_soil_solver_record) try parseSoilSolverParameters(&tokens) else soil_solver_properties.compatibilityParameters();
    const first_after_soil_solver = if (has_soil_solver_record) try next(&tokens, "soil-process, surface-energy, or site record") else first_after_runtime;
    const has_soil_process_record = std.ascii.eqlIgnoreCase(first_after_soil_solver, "soil_process");
    if (has_geospatial_grid and !has_soil_process_record)
        return error.MissingSoilProcessRecord;
    var soil_process_parameters: soil_hourly_workspace.RuntimeParameters = if (has_soil_process_record) .{
        .gravitational_water_potential_mpa_per_m = try nextFloat(&tokens, "gravitational water potential in MPa m-1"),
        .reference_water_vapor_diffusivity_m2_per_h = try nextFloat(&tokens, "reference water-vapor diffusivity in m2 h-1"),
        .vapor_diffusivity_reference_temperature_k = try nextFloat(&tokens, "vapor diffusivity reference temperature in K"),
        .vapor_diffusivity_temperature_exponent = try nextFloat(&tokens, "vapor diffusivity temperature exponent"),
        .vapor_pore_tortuosity = try nextFloat(&tokens, "vapor pore tortuosity"),
        .osmotic_reflection_coefficient = try nextFloat(&tokens, "osmotic reflection coefficient"),
        .macropore_radius_m = try nextFloat(&tokens, "macropore radius in m"),
        .macropore_residual_saturation = 0,
        .macropore_van_genuchten_alpha_per_m = 15,
        .macropore_van_genuchten_n = 2.68,
        .macropore_pore_connectivity = 0.5,
        .dual_domain_geometry_factor = 3,
        .dual_domain_scaling_coefficient = 0.4,
        .frozen_hydraulic_impedance_exponent = 7,
        .surface_residue_residual_water_content_m3_per_m3 = 0,
        .surface_residue_van_genuchten_alpha_per_m = 400,
        .surface_residue_van_genuchten_n = 2.5,
        .reference_water_viscosity_megagrams_per_m_s = try nextFloat(&tokens, "reference water viscosity in Mg m-1 s"),
        .water_viscosity_temperature_intercept = try nextFloat(&tokens, "water viscosity temperature intercept"),
        .water_viscosity_temperature_coefficient_per_c = try nextFloat(&tokens, "water viscosity temperature coefficient per C"),
    } else soil_hourly_workspace.compatibilityParameters();
    const first_after_soil_process = if (has_soil_process_record) try next(&tokens, "soil-gas-transport, soil-phase-heat, surface-energy, or site record") else first_after_soil_solver;
    const has_macropore_van_genuchten_record =
        std.ascii.eqlIgnoreCase(first_after_soil_process, "macropore_van_genuchten");
    if (has_geospatial_grid and !has_macropore_van_genuchten_record)
        return error.MissingMacroporeVanGenuchtenRecord;
    if (has_macropore_van_genuchten_record) {
        soil_process_parameters.macropore_residual_saturation =
            try nextFloat(&tokens, "macropore residual saturation");
        soil_process_parameters.macropore_van_genuchten_alpha_per_m =
            try nextFloat(&tokens, "macropore van Genuchten alpha in m-1");
        soil_process_parameters.macropore_van_genuchten_n =
            try nextFloat(&tokens, "macropore van Genuchten n");
        soil_process_parameters.macropore_pore_connectivity =
            try nextFloat(&tokens, "macropore Mualem pore connectivity");
    }
    const first_after_macropore_van_genuchten = if (has_macropore_van_genuchten_record)
        try next(&tokens, "soil-gas-transport, soil-phase-heat, surface-energy, or site record")
    else
        first_after_soil_process;
    const has_dual_domain_exchange_record =
        std.ascii.eqlIgnoreCase(first_after_macropore_van_genuchten, "dual_domain_exchange");
    if (has_geospatial_grid and !has_dual_domain_exchange_record)
        return error.MissingDualDomainExchangeRecord;
    if (has_dual_domain_exchange_record) {
        soil_process_parameters.dual_domain_geometry_factor =
            try nextFloat(&tokens, "dual-domain geometry factor");
        soil_process_parameters.dual_domain_scaling_coefficient =
            try nextFloat(&tokens, "dual-domain water-transfer scaling coefficient");
    }
    const first_after_dual_domain_exchange = if (has_dual_domain_exchange_record)
        try next(&tokens, "frozen-hydraulic-impedance, soil-gas-transport, soil-phase-heat, surface-energy, or site record")
    else
        first_after_macropore_van_genuchten;
    const has_frozen_hydraulic_impedance_record =
        std.ascii.eqlIgnoreCase(first_after_dual_domain_exchange, "frozen_hydraulic_impedance");
    if (has_geospatial_grid and !has_frozen_hydraulic_impedance_record)
        return error.MissingFrozenHydraulicImpedanceRecord;
    if (has_frozen_hydraulic_impedance_record)
        soil_process_parameters.frozen_hydraulic_impedance_exponent =
            try nextFloat(&tokens, "frozen hydraulic impedance exponent");
    try soil_hourly_workspace.validateRuntimeParameters(soil_process_parameters);
    const first_after_frozen_hydraulic_impedance =
        if (has_frozen_hydraulic_impedance_record)
            try next(&tokens, "surface-residue-freeze-thaw, soil-gas-transport, soil-phase-heat, surface-energy, or site record")
        else
            first_after_dual_domain_exchange;
    const has_surface_residue_freeze_thaw_record =
        std.ascii.eqlIgnoreCase(first_after_frozen_hydraulic_impedance, "surface_residue_freeze_thaw");
    if (has_geospatial_grid and !has_surface_residue_freeze_thaw_record)
        return error.MissingSurfaceResidueFreezeThawRecord;
    if (has_surface_residue_freeze_thaw_record) {
        soil_process_parameters.surface_residue_residual_water_content_m3_per_m3 =
            try nextFloat(&tokens, "surface residue residual water content in m3 m-3");
        soil_process_parameters.surface_residue_van_genuchten_alpha_per_m =
            try nextFloat(&tokens, "surface residue van Genuchten alpha in m-1");
        soil_process_parameters.surface_residue_van_genuchten_n =
            try nextFloat(&tokens, "surface residue van Genuchten n");
    }
    try soil_hourly_workspace.validateRuntimeParameters(soil_process_parameters);
    const first_after_surface_residue_freeze_thaw =
        if (has_surface_residue_freeze_thaw_record)
            try next(&tokens, "soil-gas-transport, soil-phase-heat, surface-energy, or site record")
        else
            first_after_frozen_hydraulic_impedance;
    const has_soil_gas_transport_record = std.ascii.eqlIgnoreCase(first_after_surface_residue_freeze_thaw, "soil_gas_transport");
    if (has_geospatial_grid and !has_soil_gas_transport_record)
        return error.MissingSoilGasTransportRecord;
    var soil_gas_transport_parameters: soil_gas_transport_step.RuntimeParameters = .{};
    if (has_soil_gas_transport_record) {
        soil_gas_transport_parameters.reference_temperature_k = try nextFloat(&tokens, "gas diffusivity reference temperature in K");
        soil_gas_transport_parameters.temperature_exponent = try nextFloat(&tokens, "gas diffusivity temperature exponent");
        for (&soil_gas_transport_parameters.free_air_diffusivity_m2_per_h) |*value| value.* = try nextFloat(&tokens, "free-air gas diffusivity in m2 h-1");
        soil_gas_transport_parameters.penman_tortuosity = try nextFloat(&tokens, "Penman gas tortuosity");
        soil_gas_transport_parameters.minimum_air_filled_porosity_m3_per_m3 = try nextFloat(&tokens, "minimum gas-filled porosity");
        soil_gas_transport_parameters.water_density_g_per_m3 = try nextFloat(&tokens, "water density in g m-3");
        soil_gas_transport_parameters.water_molar_mass_g_per_mol = try nextFloat(&tokens, "water molar mass in g mol-1");
    }
    const first_after_soil_gas_transport = if (has_soil_gas_transport_record) try next(&tokens, "soil-phase-heat, surface-energy, or site record") else first_after_surface_residue_freeze_thaw;
    const has_soil_phase_heat_record = std.ascii.eqlIgnoreCase(first_after_soil_gas_transport, "soil_phase_heat");
    if (has_geospatial_grid and !has_soil_phase_heat_record)
        return error.MissingSoilPhaseHeatRecord;
    const soil_phase_heat_parameters = if (has_soil_phase_heat_record) try parseSoilPhaseHeatParameters(&tokens) else soil_process_science.compatibilityParameters();
    const first_after_soil_phase_heat = if (has_soil_phase_heat_record) try next(&tokens, "geothermal, surface-energy record or site filename") else first_after_soil_gas_transport;
    const has_geothermal_record = std.ascii.eqlIgnoreCase(first_after_soil_phase_heat, "geothermal");
    if (has_geospatial_grid and !has_geothermal_record)
        return error.MissingGeothermalRecord;
    const geothermal_controls: GeothermalControls = if (has_geothermal_record) .{
        .enabled = try nextBool(&tokens, "geothermal boundary enabled"),
        .minimum_source_depth_m = try nextFloat(&tokens, "minimum geothermal source depth in m"),
        .source_depth_below_profile_m = try nextFloat(&tokens, "geothermal source depth below profile in m"),
        .conductivity_m_mj_per_h_k = try nextFloat(&tokens, "below-profile thermal conductivity in m MJ h-1 K-1"),
        .geothermal_flux_mj_per_m2_h = try nextFloat(&tokens, "geothermal flux in MJ m-2 h-1"),
    } else .{ .enabled = true, .minimum_source_depth_m = 10, .source_depth_below_profile_m = 1, .conductivity_m_mj_per_h_k = 8.1e-3, .geothermal_flux_mj_per_m2_h = 2.052e-4 };
    if (!std.math.isFinite(geothermal_controls.minimum_source_depth_m) or geothermal_controls.minimum_source_depth_m <= 0 or !std.math.isFinite(geothermal_controls.source_depth_below_profile_m) or geothermal_controls.source_depth_below_profile_m <= 0 or !std.math.isFinite(geothermal_controls.conductivity_m_mj_per_h_k) or geothermal_controls.conductivity_m_mj_per_h_k <= 0 or !std.math.isFinite(geothermal_controls.geothermal_flux_mj_per_m2_h)) return error.InvalidGeothermalControls;
    const first_after_geothermal = if (has_geothermal_record) try next(&tokens, "subsurface-state, surface-energy record or site filename") else first_after_soil_phase_heat;
    const has_subsurface_state_record = std.ascii.eqlIgnoreCase(first_after_geothermal, "subsurface_state");
    if (has_geospatial_grid and !has_subsurface_state_record)
        return error.MissingSubsurfaceStateRecord;
    const water_table_air_fraction_threshold = if (has_subsurface_state_record) try nextFloat(&tokens, "water-table air-filled porosity threshold") else 1.0e-3;
    const active_layer_ice_fraction_threshold = if (has_subsurface_state_record) try nextFloat(&tokens, "active-layer ice fraction threshold") else 1.0e-6;
    if (!std.math.isFinite(water_table_air_fraction_threshold) or water_table_air_fraction_threshold < 0 or water_table_air_fraction_threshold > 1 or !std.math.isFinite(active_layer_ice_fraction_threshold) or active_layer_ice_fraction_threshold < 0 or active_layer_ice_fraction_threshold > 1) return error.InvalidSubsurfaceStateControls;
    const first_after_subsurface_state = if (has_subsurface_state_record) try next(&tokens, "surface-energy record or site filename") else first_after_geothermal;
    const has_soil_geometry_record = std.ascii.eqlIgnoreCase(first_after_subsurface_state, "soil_geometry");
    if (has_geospatial_grid and !has_soil_geometry_record)
        return error.MissingSoilGeometryRecord;
    const soil_geometry_parameters: SoilGeometryParameters = if (has_soil_geometry_record) .{
        .organic_carbon_specific_volume_m3_per_g = try nextFloat(&tokens, "organic carbon specific volume in m3 g-1"),
        .organic_horizon_threshold_g_c_per_Mg = try nextFloat(&tokens, "organic horizon threshold in g C Mg-1"),
        .ice_to_water_specific_volume_difference = try nextFloat(&tokens, "ice to water specific volume difference"),
        .minimum_layer_thickness_m = try nextFloat(&tokens, "minimum soil layer thickness in m"),
    } else .{
        .organic_carbon_specific_volume_m3_per_g = 1.82e-6,
        .organic_horizon_threshold_g_c_per_Mg = 110_000,
        .ice_to_water_specific_volume_difference = 0.083,
        .minimum_layer_thickness_m = 1.0e-9,
    };
    inline for (@typeInfo(SoilGeometryParameters).@"struct".fields) |field| {
        const value = @field(soil_geometry_parameters, field.name);
        if (!std.math.isFinite(value) or value <= 0) return error.InvalidSoilGeometryParameters;
    }
    if (soil_geometry_parameters.ice_to_water_specific_volume_difference >= 1) return error.InvalidSoilGeometryParameters;
    const first_after_soil_geometry = if (has_soil_geometry_record) try next(&tokens, "surface-pond-energy, surface-energy record or site filename") else first_after_subsurface_state;
    const has_surface_pond_energy_record = std.ascii.eqlIgnoreCase(first_after_soil_geometry, "surface_pond_energy");
    if (has_geospatial_grid and !has_surface_pond_energy_record)
        return error.MissingSurfacePondEnergyRecord;
    const surface_pond_dry_organic_heat_capacity_mj_per_g_c_k = if (has_surface_pond_energy_record) try nextFloat(&tokens, "surface pond dry-organic heat capacity in MJ g-C-1 K-1") else 2.496e-6;
    const surface_pond_activation_heat_capacity_mj_per_m2_k = if (has_surface_pond_energy_record) try nextFloat(&tokens, "surface pond activation heat capacity in MJ m-2 K-1") else 0.838e-3;
    if (!std.math.isFinite(surface_pond_dry_organic_heat_capacity_mj_per_g_c_k) or surface_pond_dry_organic_heat_capacity_mj_per_g_c_k <= 0 or !std.math.isFinite(surface_pond_activation_heat_capacity_mj_per_m2_k) or surface_pond_activation_heat_capacity_mj_per_m2_k <= 0) return error.InvalidSurfacePondEnergyControls;
    const first_after_surface_pond_energy = if (has_surface_pond_energy_record) try next(&tokens, "surface-energy record or site filename") else first_after_soil_geometry;
    const has_surface_energy_record = std.ascii.eqlIgnoreCase(first_after_surface_pond_energy, "surface_energy");
    if (has_geospatial_grid and !has_surface_energy_record)
        return error.MissingSurfaceEnergyRecord;
    const soil_longwave_emissivity = if (has_surface_energy_record) try nextFloat(&tokens, "soil longwave emissivity") else 0.97;
    const snow_longwave_emissivity = if (has_surface_energy_record) try nextFloat(&tokens, "snow longwave emissivity") else 0.97;
    const canopy_longwave_emissivity = if (has_surface_energy_record) try nextFloat(&tokens, "canopy longwave emissivity") else 0.97;
    const snow_full_cover_depth_m = if (has_surface_energy_record) try nextFloat(&tokens, "snow full-cover depth in m") else 0.07;
    const surface_sensible_heat_conductance_mj_per_m2_h_k = if (has_surface_energy_record) try nextFloat(&tokens, "surface sensible heat conductance in MJ m-2 h-1 K-1") else 0.43;
    const surface_latent_heat_conductance_mj_per_m2_h_kpa = if (has_surface_energy_record) try nextFloat(&tokens, "surface latent heat conductance in MJ m-2 h-1 kPa-1") else 0.0;
    const surface_vapor_activity_fraction = if (has_surface_energy_record) try nextFloat(&tokens, "surface vapor activity fraction") else 1.0;
    const minimum_surface_temperature_k = if (has_surface_energy_record) try nextFloat(&tokens, "minimum surface temperature in K") else 173.15;
    const maximum_surface_temperature_k = if (has_surface_energy_record) try nextFloat(&tokens, "maximum surface temperature in K") else 373.15;
    if (!std.math.isFinite(soil_longwave_emissivity) or soil_longwave_emissivity < 0 or soil_longwave_emissivity > 1 or
        !std.math.isFinite(snow_longwave_emissivity) or snow_longwave_emissivity < 0 or snow_longwave_emissivity > 1 or
        !std.math.isFinite(canopy_longwave_emissivity) or canopy_longwave_emissivity < 0 or canopy_longwave_emissivity > 1 or
        !std.math.isFinite(snow_full_cover_depth_m) or snow_full_cover_depth_m <= 0 or
        !std.math.isFinite(surface_sensible_heat_conductance_mj_per_m2_h_k) or surface_sensible_heat_conductance_mj_per_m2_h_k <= 0 or
        !std.math.isFinite(surface_latent_heat_conductance_mj_per_m2_h_kpa) or surface_latent_heat_conductance_mj_per_m2_h_kpa < 0 or
        !std.math.isFinite(surface_vapor_activity_fraction) or surface_vapor_activity_fraction < 0 or surface_vapor_activity_fraction > 1 or
        !std.math.isFinite(minimum_surface_temperature_k) or !std.math.isFinite(maximum_surface_temperature_k) or minimum_surface_temperature_k <= 0 or minimum_surface_temperature_k >= maximum_surface_temperature_k) return error.InvalidSurfaceEnergyControls;
    const first_after_surface_energy = if (has_surface_energy_record) try next(&tokens, "snow-layers, plant-structure record or site filename") else first_after_surface_pond_energy;
    const has_snow_layers_record = std.ascii.eqlIgnoreCase(first_after_surface_energy, "snow_layers");
    const snow_layer_count = if (has_snow_layers_record) try nextUnsigned(&tokens, "snow layer count") else 5;
    if (snow_layer_count == 0) return error.NoSnowLayers;
    const snow_layer_bottom_depth_m = try allocator.alloc(f64, snow_layer_count);
    errdefer allocator.free(snow_layer_bottom_depth_m);
    if (has_snow_layers_record) {
        for (snow_layer_bottom_depth_m) |*bottom| bottom.* = try nextFloat(&tokens, "snow layer bottom depth in m");
    } else {
        const source_bottoms = [_]f64{ 0.05, 0.125, 0.25, 0.50, 1.00 };
        @memcpy(snow_layer_bottom_depth_m, &source_bottoms);
    }
    const initial_snow_density_Mg_per_m3 = if (has_snow_layers_record) try nextFloat(&tokens, "initial snow density in Mg m-3") else 0.05;
    const snow_ice_density_Mg_per_m3 = if (has_snow_layers_record) try nextFloat(&tokens, "snow ice density in Mg m-3") else 0.92;
    const snow_latent_heat_of_fusion_mj_per_m3 = if (has_snow_layers_record) try nextFloat(&tokens, "snow latent heat of fusion in MJ m-3") else 333;
    const snow_phase_damping_divisor = if (has_snow_layers_record) try nextFloat(&tokens, "snow phase damping divisor") else 2.7185;
    var previous_snow_bottom: f64 = 0;
    for (snow_layer_bottom_depth_m) |bottom| {
        if (!std.math.isFinite(bottom) or bottom <= previous_snow_bottom) return error.InvalidSnowLayerBoundaries;
        previous_snow_bottom = bottom;
    }
    if (!std.math.isFinite(initial_snow_density_Mg_per_m3) or initial_snow_density_Mg_per_m3 <= 0 or !std.math.isFinite(snow_ice_density_Mg_per_m3) or snow_ice_density_Mg_per_m3 <= 0 or !std.math.isFinite(snow_latent_heat_of_fusion_mj_per_m3) or snow_latent_heat_of_fusion_mj_per_m3 <= 0 or !std.math.isFinite(snow_phase_damping_divisor) or snow_phase_damping_divisor <= 0) return error.InvalidSnowPhaseParameters;
    const first_after_snow_layers = if (has_snow_layers_record) try next(&tokens, "snow-compaction, plant-structure record or site filename") else first_after_surface_energy;
    const has_snow_compaction_record = std.ascii.eqlIgnoreCase(first_after_snow_layers, "snow_compaction");
    const snow_compaction_parameters: snow_compaction.Parameters = if (has_snow_compaction_record) .{
        .maximum_temperature_metamorphism_density_Mg_per_m3 = try nextFloat(&tokens, "snow metamorphism density threshold in Mg m-3"),
        .temperature_metamorphism_rate_per_h = try nextFloat(&tokens, "snow temperature metamorphism rate per h"),
        .temperature_metamorphism_exponent_per_c = try nextFloat(&tokens, "snow temperature metamorphism exponent per C"),
        .viscosity_scale_Mg_h_per_m3 = try nextFloat(&tokens, "snow viscosity scale in Mg h m-3"),
        .viscosity_temperature_exponent_per_c = try nextFloat(&tokens, "snow viscosity temperature exponent per C"),
        .viscosity_density_exponent_m3_per_Mg = try nextFloat(&tokens, "snow viscosity density exponent in m3 Mg-1"),
        .minimum_snowfall_temperature_c = try nextFloat(&tokens, "minimum snowfall density temperature in C"),
        .maximum_snowfall_temperature_c = try nextFloat(&tokens, "maximum snowfall density temperature in C"),
        .snowfall_density_temperature_coefficient_Mg_per_m3_c_pow_1_5 = try nextFloat(&tokens, "snowfall density temperature coefficient"),
    } else .{ .maximum_temperature_metamorphism_density_Mg_per_m3 = 0.25, .temperature_metamorphism_rate_per_h = 1e-5, .temperature_metamorphism_exponent_per_c = 0.04, .viscosity_scale_Mg_h_per_m3 = 0.25, .viscosity_temperature_exponent_per_c = -0.08, .viscosity_density_exponent_m3_per_Mg = 23, .minimum_snowfall_temperature_c = -15, .maximum_snowfall_temperature_c = 2, .snowfall_density_temperature_coefficient_Mg_per_m3_c_pow_1_5 = 1.7e-3 };
    const first_after_snow_compaction = if (has_snow_compaction_record) try next(&tokens, "snow-thermal, plant-structure record or site filename") else first_after_snow_layers;
    const has_snow_thermal_record = std.ascii.eqlIgnoreCase(first_after_snow_compaction, "snow_thermal");
    const snow_heat_conduction_parameters: snow_heat_conduction.Parameters = .{
        .conductivity_scale_m_mj_per_h_k = if (has_snow_thermal_record) try nextFloat(&tokens, "snow conductivity scale in m MJ h-1 K-1") else 0.0036,
        .conductivity_density_exponent_m3_per_Mg = if (has_snow_thermal_record) try nextFloat(&tokens, "snow conductivity density exponent in m3 Mg-1") else 2.650,
        .conductivity_log10_intercept = if (has_snow_thermal_record) try nextFloat(&tokens, "snow conductivity log10 intercept") else -1.652,
        .maximum_effective_density_Mg_per_m3 = if (has_snow_thermal_record) try nextFloat(&tokens, "maximum effective snow density in Mg m-3") else 0.6,
        .ice_density_Mg_per_m3 = snow_ice_density_Mg_per_m3,
    };
    const first_after_snow_thermal = if (has_snow_thermal_record) try next(&tokens, "snow-vapor, plant-structure record or site filename") else first_after_snow_compaction;
    const has_snow_vapor_record = std.ascii.eqlIgnoreCase(first_after_snow_thermal, "snow_vapor");
    const snow_vapor_parameters: snow_vapor_equilibrium.Parameters = .{
        .vapor_volume_prefactor_k = if (has_snow_vapor_record) try nextFloat(&tokens, "snow vapor volume prefactor in K") else 2.173e-3,
        .equilibrium_relative_humidity = if (has_snow_vapor_record) try nextFloat(&tokens, "snow equilibrium relative humidity") else 0.61,
        .clausius_clapeyron_temperature_k = if (has_snow_vapor_record) try nextFloat(&tokens, "snow Clausius-Clapeyron temperature in K") else 5360,
        .reference_inverse_temperature_per_k = if (has_snow_vapor_record) try nextFloat(&tokens, "snow reference inverse temperature per K") else 3.661e-3,
        .liquid_evaporation_latent_heat_mj_per_m3 = if (has_snow_vapor_record) try nextFloat(&tokens, "liquid evaporation latent heat in MJ m-3") else 2465,
        .snow_sublimation_latent_heat_mj_per_m3 = if (has_snow_vapor_record) try nextFloat(&tokens, "snow sublimation latent heat in MJ m-3") else 2834,
    };
    const first_after_snow_vapor = if (has_snow_vapor_record) try next(&tokens, "snow-vapor-transport, plant-structure record or site filename") else first_after_snow_thermal;
    const has_snow_vapor_transport_record = std.ascii.eqlIgnoreCase(first_after_snow_vapor, "snow_vapor_transport");
    const snow_vapor_diffusion_parameters: snow_vapor_diffusion.Parameters = .{
        .reference_vapor_diffusivity_m2_per_h = if (has_snow_vapor_transport_record) try nextFloat(&tokens, "snow reference vapor diffusivity in m2 h-1") else soil_process_parameters.reference_water_vapor_diffusivity_m2_per_h,
        .reference_temperature_k = if (has_snow_vapor_transport_record) try nextFloat(&tokens, "snow vapor diffusivity reference temperature in K") else soil_process_parameters.vapor_diffusivity_reference_temperature_k,
        .temperature_exponent = if (has_snow_vapor_transport_record) try nextFloat(&tokens, "snow vapor diffusivity temperature exponent") else soil_process_parameters.vapor_diffusivity_temperature_exponent,
        .minimum_air_fraction = if (has_snow_vapor_transport_record) try nextFloat(&tokens, "minimum snow air fraction") else 0,
        .vapor_sensible_heat_capacity_mj_per_m3_k = if (has_snow_vapor_transport_record) try nextFloat(&tokens, "snow vapor sensible heat capacity in MJ m-3 K-1") else 4.19,
    };
    const first_after_snow_vapor_transport = if (has_snow_vapor_transport_record) try next(&tokens, "surface-gas-resistance, plant-structure record or site filename") else first_after_snow_vapor;
    const has_surface_gas_resistance_record = std.ascii.eqlIgnoreCase(first_after_snow_vapor_transport, "surface_gas_resistance");
    const surface_gas_resistance_parameters: surface_gas_boundary_conductance.Parameters = .{
        .minimum_richardson_number = if (has_surface_gas_resistance_record) try nextFloat(&tokens, "minimum Richardson number") else -0.10,
        .maximum_richardson_number = if (has_surface_gas_resistance_record) try nextFloat(&tokens, "maximum Richardson number") else 0.05,
        .richardson_resistance_multiplier = if (has_surface_gas_resistance_record) try nextFloat(&tokens, "Richardson resistance multiplier") else 10,
        .minimum_aerodynamic_resistance_h_per_m = if (has_surface_gas_resistance_record) try nextFloat(&tokens, "minimum aerodynamic resistance in h m-1") else 1e-6,
        .maximum_aerodynamic_resistance_h_per_m = if (has_surface_gas_resistance_record) try nextFloat(&tokens, "maximum aerodynamic resistance in h m-1") else 1e6,
        .canopy_drag_length_m = if (has_surface_gas_resistance_record) try nextFloat(&tokens, "canopy drag length in m") else 2e-4,
        .minimum_air_fraction = if (has_surface_gas_resistance_record) try nextFloat(&tokens, "minimum surface gas air fraction") else 1e-6,
    };
    const first_after_surface_gas_resistance = if (has_surface_gas_resistance_record) try next(&tokens, "surface-runoff, surface-aerodynamics, plant-structure record or site filename") else first_after_snow_vapor_transport;
    const has_surface_runoff_record = std.ascii.eqlIgnoreCase(first_after_surface_gas_resistance, "surface_runoff");
    if (has_geospatial_grid and !has_surface_runoff_record)
        return error.MissingSurfaceRunoffRecord;
    const surface_runoff_parameters: surface_runoff.Parameters = .{
        .ground_surface_retention_m3_per_m2 = if (has_surface_runoff_record) try nextFloat(&tokens, "ground surface retention in m3 m-2") else 0.005,
        .runoff_roughness_h_per_m_one_third = if (has_surface_runoff_record) try nextFloat(&tokens, "runoff roughness coefficient") else 0.011,
        .maximum_hydraulic_volume_m3 = if (has_surface_runoff_record) try nextFloat(&tokens, "maximum hydraulic runoff volume in m3") else 1.0e-3,
        .manning_time_conversion_s_per_h = if (has_surface_runoff_record) try nextFloat(&tokens, "Manning seconds per hour") else 3.6e3,
        .negligible_water_m3 = if (has_surface_runoff_record) try nextFloat(&tokens, "negligible surface water in m3") else 1.0e-12,
    };
    const first_after_surface_runoff = if (has_surface_runoff_record) try next(&tokens, "rainfall-impact, surface-aerodynamics, plant-structure record or site filename") else first_after_surface_gas_resistance;
    const has_rainfall_impact_record = std.ascii.eqlIgnoreCase(first_after_surface_runoff, "rainfall_impact");
    if (has_geospatial_grid and !has_rainfall_impact_record)
        return error.MissingRainfallImpactRecord;
    const rainfall_impact_parameters: surface_precipitation.RainfallImpactParameters = .{
        .direct_energy_intercept_j_per_mm = if (has_rainfall_impact_record) try nextFloat(&tokens, "direct rainfall energy intercept in J mm-1") else 8.95,
        .direct_energy_log_coefficient_j_per_mm = if (has_rainfall_impact_record) try nextFloat(&tokens, "direct rainfall logarithmic energy coefficient in J mm-1") else 8.44,
        .throughfall_energy_height_coefficient_j_per_mm_sqrt_m = if (has_rainfall_impact_record) try nextFloat(&tokens, "throughfall energy height coefficient") else 15.8,
        .throughfall_energy_intercept_j_per_mm = if (has_rainfall_impact_record) try nextFloat(&tokens, "throughfall energy intercept in J mm-1") else 5.87,
        .maximum_canopy_height_m = if (has_rainfall_impact_record) try nextFloat(&tokens, "maximum canopy height in rainfall energy calculation") else 2.5,
        .ponding_attenuation_per_mm = if (has_rainfall_impact_record) try nextFloat(&tokens, "rainfall energy ponding attenuation per mm") else 2,
        .conductivity_damage_per_j_per_megagram_per_megagram = if (has_rainfall_impact_record) try nextFloat(&tokens, "conductivity damage per J per Mg Mg-1") else 1.0e-3,
        .conductivity_recovery_fraction_per_h = if (has_rainfall_impact_record) try nextFloat(&tokens, "rainfall conductivity recovery fraction per h") else 5.0e-4,
    };
    const first_after_rainfall_impact = if (has_rainfall_impact_record) try next(&tokens, "surface-aerodynamics, plant-structure record or site filename") else first_after_surface_runoff;
    const has_surface_aerodynamics_record = std.ascii.eqlIgnoreCase(first_after_rainfall_impact, "surface_aerodynamics");
    const surface_aerodynamic_parameters: surface_aerodynamics.Parameters = .{
        .canopy_area_attenuation = if (has_surface_aerodynamics_record) try nextFloat(&tokens, "canopy area attenuation") else 0.5,
        .minimum_reference_height_above_displacement_m = if (has_surface_aerodynamics_record) try nextFloat(&tokens, "minimum reference height above displacement in m") else 2,
        .snow_roughness_height_m = if (has_surface_aerodynamics_record) try nextFloat(&tokens, "snow roughness height in m") else 0.005,
        .soil_roughness_height_m = if (has_surface_aerodynamics_record) try nextFloat(&tokens, "soil roughness height in m") else 0.025,
        .richardson_coefficient_k_m = if (has_surface_aerodynamics_record) try nextFloat(&tokens, "bulk Richardson coefficient in K m") else 1.27e8,
        .neutral_resistance_denominator = if (has_surface_aerodynamics_record) try nextFloat(&tokens, "neutral aerodynamic resistance denominator") else 0.168,
        .minimum_wind_speed_m_per_h = if (has_surface_aerodynamics_record) try nextFloat(&tokens, "minimum wind speed in m h-1") else 1e-12,
        .minimum_aerodynamic_resistance_h_per_m = if (has_surface_aerodynamics_record) try nextFloat(&tokens, "minimum neutral aerodynamic resistance in h m-1") else 0.00139,
        .weather_height_includes_displacement = if (has_surface_aerodynamics_record) try nextBool(&tokens, "weather height includes displacement") else false,
    };
    const first_after_surface_aerodynamics = if (has_surface_aerodynamics_record) try next(&tokens, "ground-air, plant-structure record or site filename") else first_after_rainfall_impact;
    const has_ground_air_record = std.ascii.eqlIgnoreCase(first_after_surface_aerodynamics, "ground_air");
    const ground_air_parameters: ground_air_exchange.Parameters = .{
        .minimum_richardson_number = surface_gas_resistance_parameters.minimum_richardson_number,
        .maximum_richardson_number = surface_gas_resistance_parameters.maximum_richardson_number,
        .richardson_resistance_multiplier = surface_gas_resistance_parameters.richardson_resistance_multiplier,
        .minimum_aerodynamic_resistance_h_per_m = surface_gas_resistance_parameters.minimum_aerodynamic_resistance_h_per_m,
        .maximum_aerodynamic_resistance_h_per_m = surface_gas_resistance_parameters.maximum_aerodynamic_resistance_h_per_m,
        .volumetric_air_heat_capacity_mj_per_m3_k = if (has_ground_air_record) try nextFloat(&tokens, "ground-air volumetric heat capacity in MJ m-3 K-1") else 1.25e-3,
        .minimum_air_column_height_m = if (has_ground_air_record) try nextFloat(&tokens, "minimum ground-air column height in m") else 5,
        .sensible_heat_conductivity_mj_per_m_h_k = if (has_ground_air_record) try nextFloat(&tokens, "ground-air sensible heat conductivity in MJ m-1 h-1 K-1") else 1.25e-3,
        .liquid_water_latent_heat_mj_per_m3 = if (has_ground_air_record) try nextFloat(&tokens, "liquid-water latent heat in MJ m-3") else 2465,
        .saturation_vapor_prefactor_k = if (has_ground_air_record) try nextFloat(&tokens, "ground-air saturation vapor prefactor in K") else 2.173e-3,
        .saturation_relative_humidity = if (has_ground_air_record) try nextFloat(&tokens, "ground-air saturation relative humidity") else 0.61,
        .saturation_temperature_k = if (has_ground_air_record) try nextFloat(&tokens, "ground-air saturation temperature coefficient in K") else 5360,
        .saturation_reference_inverse_temperature_per_k = if (has_ground_air_record) try nextFloat(&tokens, "ground-air saturation reference inverse temperature per K") else 3.661e-3,
    };
    const first_after_ground_air = if (has_ground_air_record) try next(&tokens, "canopy-surface-exchange, plant-structure record or site filename") else first_after_surface_aerodynamics;
    const has_canopy_surface_exchange_record = std.ascii.eqlIgnoreCase(first_after_ground_air, "canopy_surface_exchange");
    const canopy_sensible_surface_resistance_h_per_m = if (has_canopy_surface_exchange_record) try nextFloat(&tokens, "canopy sensible surface resistance in h m-1") else 0.00139;
    const canopy_latent_surface_resistance_h_per_m = if (has_canopy_surface_exchange_record) try nextFloat(&tokens, "canopy latent surface resistance in h m-1") else 0.0278;
    const canopy_minimum_boundary_resistance_h_per_m = if (has_canopy_surface_exchange_record) try nextFloat(&tokens, "minimum canopy boundary resistance in h m-1") else 0.00139;
    const canopy_maximum_boundary_resistance_h_per_m = if (has_canopy_surface_exchange_record) try nextFloat(&tokens, "maximum canopy boundary resistance in h m-1") else 0.0139;
    const canopy_surface_exchange_parameters: canopy_surface_exchange.Parameters = .{
        .minimum_richardson_number = ground_air_parameters.minimum_richardson_number,
        .maximum_richardson_number = ground_air_parameters.maximum_richardson_number,
        .richardson_resistance_multiplier = ground_air_parameters.richardson_resistance_multiplier,
        .minimum_boundary_resistance_h_per_m = canopy_minimum_boundary_resistance_h_per_m,
        .maximum_boundary_resistance_h_per_m = canopy_maximum_boundary_resistance_h_per_m,
        .saturation_vapor_prefactor_k = ground_air_parameters.saturation_vapor_prefactor_k,
        .saturation_relative_humidity = ground_air_parameters.saturation_relative_humidity,
        .saturation_temperature_k = ground_air_parameters.saturation_temperature_k,
        .saturation_reference_inverse_temperature_per_k = ground_air_parameters.saturation_reference_inverse_temperature_per_k,
        .water_potential_vapor_coefficient_mol_per_m3 = if (has_canopy_surface_exchange_record) try nextFloat(&tokens, "canopy water-potential vapor coefficient in mol m-3") else 18,
        .universal_gas_constant_j_per_mol_k = if (has_canopy_surface_exchange_record) try nextFloat(&tokens, "universal gas constant in J mol-1 K-1") else 8.3143,
        .latent_heat_of_vaporization_mj_per_m3 = ground_air_parameters.liquid_water_latent_heat_mj_per_m3,
        .liquid_water_heat_capacity_mj_per_m3_k = if (has_canopy_surface_exchange_record) try nextFloat(&tokens, "liquid-water heat capacity in MJ m-3 K-1") else 4.19,
    };
    if (!std.math.isFinite(canopy_sensible_surface_resistance_h_per_m) or canopy_sensible_surface_resistance_h_per_m < 0 or !std.math.isFinite(canopy_latent_surface_resistance_h_per_m) or canopy_latent_surface_resistance_h_per_m < 0) return error.InvalidCanopySurfaceResistance;
    const first_after_canopy_surface_exchange = if (has_canopy_surface_exchange_record) try next(&tokens, "canopy-ammonia-exchange, plant-structure record or site filename") else first_after_ground_air;
    const has_canopy_ammonia_exchange_record = std.ascii.eqlIgnoreCase(first_after_canopy_surface_exchange, "canopy_ammonia_exchange");
    if (has_geospatial_grid and !has_canopy_ammonia_exchange_record)
        return error.MissingCanopyAmmoniaExchangeRecord;
    const canopy_ammonia_exchange_parameters: plant_soil_exchange.CanopyAmmoniaExchangeParameters = if (has_canopy_ammonia_exchange_record) .{
        .minimum_canopy_dry_matter_fraction = try nextFloat(&tokens, "minimum canopy dry-matter fraction"),
        .water_potential_dry_matter_increment = try nextFloat(&tokens, "water-potential dry-matter increment"),
        .water_potential_denominator_per_mpa = try nextFloat(&tokens, "water-potential denominator coefficient per MPa"),
        .water_potential_denominator_offset = try nextFloat(&tokens, "water-potential denominator offset"),
        .canopy_air_volume_ratio_m3_per_g_c = try nextFloat(&tokens, "canopy air-volume ratio in m3 g-C-1"),
        .maximum_mobile_nitrogen_transfer_fraction_per_step = try nextFloat(&tokens, "maximum mobile-N transfer fraction per step"),
        .solubility_log_intercept = try nextFloat(&tokens, "ammonia solubility log intercept"),
        .solubility_temperature_coefficient_per_c = try nextFloat(&tokens, "ammonia solubility temperature coefficient per degC"),
    } else plant_soil_exchange.compatibilityCanopyAmmoniaExchangeParameters();
    try canopy_ammonia_exchange_parameters.validate();
    const first_after_canopy_ammonia_exchange = if (has_canopy_ammonia_exchange_record) try next(&tokens, "plant-structure record or site filename") else first_after_canopy_surface_exchange;
    const has_plant_structure_record = std.ascii.eqlIgnoreCase(first_after_canopy_ammonia_exchange, "plant_structure");
    if (has_geospatial_grid and !has_plant_structure_record)
        return error.MissingPlantStructureRecord;
    const root_axes_per_plant = if (has_plant_structure_record) try nextUnsigned(&tokens, "root axes per plant") else 10;
    if (root_axes_per_plant == 0) return error.NoRootAxes;
    const first_after_plant_structure = if (has_plant_structure_record) try next(&tokens, "canopy-layers, shoot-controls, root-gas, or site record") else first_after_canopy_ammonia_exchange;
    const has_canopy_layers_record = std.ascii.eqlIgnoreCase(first_after_plant_structure, "canopy_layers");
    if (has_geospatial_grid and !has_canopy_layers_record)
        return error.MissingCanopyLayersRecord;
    const canopy_layer_count = if (has_canopy_layers_record) try nextUnsigned(&tokens, "canopy layer count") else 10;
    if (canopy_layer_count == 0) return error.NoCanopyLayers;
    const first_after_canopy_layers = if (has_canopy_layers_record) try next(&tokens, "canopy-geometry, shoot-controls, root-gas, or site record") else first_after_plant_structure;
    const has_canopy_geometry_record = std.ascii.eqlIgnoreCase(first_after_canopy_layers, "canopy_geometry");
    if (has_geospatial_grid and !has_canopy_geometry_record)
        return error.MissingCanopyGeometryRecord;
    const stalk_volume_m3_per_g_c = if (has_canopy_geometry_record) try nextFloat(&tokens, "stalk volume per carbon mass in m3 g-C-1") else 4.0e-6;
    if (!std.math.isFinite(stalk_volume_m3_per_g_c) or stalk_volume_m3_per_g_c <= 0) return error.InvalidCanopyGeometryControls;
    const first_after_canopy_geometry = if (has_canopy_geometry_record) try next(&tokens, "shoot-controls, root-gas, or site record") else first_after_canopy_layers;
    const has_canopy_discretization_record = std.ascii.eqlIgnoreCase(first_after_canopy_geometry, "canopy_discretization");
    if (has_geospatial_grid and !has_canopy_discretization_record)
        return error.MissingCanopyDiscretizationRecord;
    const canopy_discretization: @import("canopy_geometry.zig").Discretization = if (has_canopy_discretization_record) .{
        .leaf_inclination_class_count = try nextUnsigned(&tokens, "leaf inclination class count"),
        .leaf_azimuth_class_count = try nextUnsigned(&tokens, "leaf azimuth class count"),
        .diffuse_sky_sector_count = try nextUnsigned(&tokens, "diffuse sky sector count"),
    } else .{};
    try @import("canopy_geometry.zig").validateDiscretization(canopy_discretization);
    const first_after_canopy_discretization = if (has_canopy_discretization_record) try next(&tokens, "standing-dead-energy, shoot-controls, root-gas, or site record") else first_after_canopy_geometry;
    const has_standing_dead_partition_record = std.ascii.eqlIgnoreCase(first_after_canopy_discretization, "standing_dead_partition");
    if (has_geospatial_grid and !has_standing_dead_partition_record)
        return error.MissingStandingDeadPartitionRecord;
    const standing_dead_partition_parameters: plant_initialization.StandingDeadPartitionParameters = if (has_standing_dead_partition_record) .{
        .carbon_fraction = .{ try nextFloat(&tokens, "standing-dead carbon fraction 1"), try nextFloat(&tokens, "standing-dead carbon fraction 2"), try nextFloat(&tokens, "standing-dead carbon fraction 3"), try nextFloat(&tokens, "standing-dead carbon fraction 4") },
        .nitrogen_weight = .{ try nextFloat(&tokens, "standing-dead nitrogen weight 1"), try nextFloat(&tokens, "standing-dead nitrogen weight 2"), try nextFloat(&tokens, "standing-dead nitrogen weight 3"), try nextFloat(&tokens, "standing-dead nitrogen weight 4") },
        .phosphorus_weight = .{ try nextFloat(&tokens, "standing-dead phosphorus weight 1"), try nextFloat(&tokens, "standing-dead phosphorus weight 2"), try nextFloat(&tokens, "standing-dead phosphorus weight 3"), try nextFloat(&tokens, "standing-dead phosphorus weight 4") },
    } else plant_initialization.compatibilityStandingDeadPartitionParameters();
    try standing_dead_partition_parameters.validate();
    const first_after_standing_dead_partition = if (has_standing_dead_partition_record) try next(&tokens, "standing-dead-energy, shoot-controls, root-gas, or site record") else first_after_canopy_discretization;
    const has_plant_heat_water_record = std.ascii.eqlIgnoreCase(first_after_standing_dead_partition, "plant_initial_heat_water");
    if (has_geospatial_grid and !has_plant_heat_water_record)
        return error.MissingPlantInitialHeatWaterRecord;
    const plant_heat_water_parameters: plant_initialization.PlantHeatWaterParameters = if (has_plant_heat_water_record) .{
        .kelvin_offset_k = try nextFloat(&tokens, "plant Kelvin offset in K"),
        .vapor_pressure_numerator_kpa_k = try nextFloat(&tokens, "initial vapor-pressure numerator in kPa K"),
        .initial_relative_humidity = try nextFloat(&tokens, "initial plant relative humidity"),
        .vapor_pressure_exponent_temperature_k = try nextFloat(&tokens, "initial vapor-pressure exponent temperature in K"),
        .vapor_pressure_reference_inverse_temperature_per_k = try nextFloat(&tokens, "initial vapor-pressure reference inverse temperature per K"),
        .initial_total_water_potential_mpa = try nextFloat(&tokens, "initial plant total water potential in MPa"),
    } else plant_initialization.compatibilityPlantHeatWaterParameters();
    try plant_heat_water_parameters.validate();
    const first_after_plant_heat_water = if (has_plant_heat_water_record) try next(&tokens, "standing-dead-energy, shoot-controls, root-gas, or site record") else first_after_standing_dead_partition;
    const has_plant_geometry_record = std.ascii.eqlIgnoreCase(first_after_plant_heat_water, "plant_initial_geometry");
    if (has_geospatial_grid and !has_plant_geometry_record)
        return error.MissingPlantInitialGeometryRecord;
    const plant_geometry_parameters: plant_initialization.PlantGeometryParameters = if (has_plant_geometry_record) .{
        .seed_volume_m3_per_g_c = try nextFloat(&tokens, "seed volume in m3 per g C"),
        .seed_length_multiplier = try nextFloat(&tokens, "seed length multiplier"),
        .seed_shape_volume_factor = try nextFloat(&tokens, "seed shape volume factor"),
        .seed_pi = try nextFloat(&tokens, "seed geometry pi"),
        .seed_length_exponent = try nextFloat(&tokens, "seed length exponent"),
        .seed_surface_area_multiplier = try nextFloat(&tokens, "seed surface-area multiplier"),
        .root_volume_numerator_m3_per_g_c = try nextFloat(&tokens, "root volume numerator in m3 per g C"),
        .root_dry_matter_fraction = try nextFloat(&tokens, "root dry-matter fraction"),
        .root_porosity_floor = try nextFloat(&tokens, "root porosity floor"),
        .root_pi = try nextFloat(&tokens, "root geometry pi"),
    } else plant_initialization.compatibilityPlantGeometryParameters();
    try plant_geometry_parameters.validate();
    const first_after_plant_geometry = if (has_plant_geometry_record) try next(&tokens, "standing-dead-energy, shoot-controls, root-gas, or site record") else first_after_plant_heat_water;
    const has_phenology_initialization_record = std.ascii.eqlIgnoreCase(first_after_plant_geometry, "plant_initial_phenology");
    if (has_geospatial_grid and !has_phenology_initialization_record)
        return error.MissingPlantInitialPhenologyRecord;
    const phenology_initialization_parameters: plant_initialization.PhenologyInitializationParameters = if (has_phenology_initialization_record) .{
        .perennial_input_scale = try nextFloat(&tokens, "perennial phenology input scale"),
        .minimum_perennial_node_scaling = try nextFloat(&tokens, "minimum perennial node scaling"),
        .perennial_maximum_concurrently_growing_nodes = try nextUnsigned(&tokens, "perennial maximum concurrently growing nodes"),
        .early_maturity_group_maximum = try nextFloat(&tokens, "early maturity-group maximum"),
        .intermediate_maturity_group_maximum = try nextFloat(&tokens, "intermediate maturity-group maximum"),
        .early_maximum_concurrently_growing_nodes = try nextUnsigned(&tokens, "early maximum concurrently growing nodes"),
        .intermediate_maximum_concurrently_growing_nodes = try nextUnsigned(&tokens, "intermediate maximum concurrently growing nodes"),
        .late_maximum_concurrently_growing_nodes = try nextUnsigned(&tokens, "late maximum concurrently growing nodes"),
    } else plant_initialization.compatibilityPhenologyInitializationParameters();
    try phenology_initialization_parameters.validate();
    const first_after_phenology_initialization = if (has_phenology_initialization_record) try next(&tokens, "standing-dead-energy, shoot-controls, root-gas, or site record") else first_after_plant_geometry;
    const has_root_initialization_record = std.ascii.eqlIgnoreCase(first_after_phenology_initialization, "root_initialization");
    if (has_geospatial_grid and !has_root_initialization_record)
        return error.MissingRootInitializationRecord;
    const root_initialization_parameters: plant_root_system.InitializationParameters = if (has_root_initialization_record) .{
        .root_nitrogen_to_maximum_protein_multiplier = try nextFloat(&tokens, "root N:C to maximum protein multiplier"),
        .root_phosphorus_to_maximum_protein_multiplier = try nextFloat(&tokens, "root P:C to maximum protein multiplier"),
        .mycorrhizal_radius_m = try nextFloat(&tokens, "initial mycorrhizal radius in m"),
        .initial_total_water_potential_mpa = try nextFloat(&tokens, "initial root total water potential in MPa"),
        .osmotic_water_potential_decrement_mpa = try nextFloat(&tokens, "initial root osmotic-potential decrement in MPa"),
        .initial_active_length_m = try nextFloat(&tokens, "initial active root length in m"),
        .initial_water_fraction = try nextFloat(&tokens, "initial root water fraction"),
    } else plant_root_system.compatibilityInitializationParameters();
    try root_initialization_parameters.validate();
    const first_after_root_initialization = if (has_root_initialization_record) try next(&tokens, "root-morphology, standing-dead-energy, shoot-controls, root-gas, or site record") else first_after_phenology_initialization;
    const has_root_morphology_record = std.ascii.eqlIgnoreCase(first_after_root_initialization, "root_morphology");
    if (has_geospatial_grid and !has_root_morphology_record)
        return error.MissingRootMorphologyRecord;
    const root_morphology_parameters: plant_root_system.MorphologyParameters = if (has_root_morphology_record) .{
        .minimum_average_secondary_length_m = try nextFloat(&tokens, "minimum average secondary root length in m"),
        .root_elastic_modulus_mpa = try nextFloat(&tokens, "root elastic modulus in MPa"),
    } else plant_root_system.compatibilityMorphologyParameters();
    try root_morphology_parameters.validate();
    const first_after_root_morphology = if (has_root_morphology_record) try next(&tokens, "standing-dead-energy, shoot-controls, root-gas, or site record") else first_after_root_initialization;
    const has_standing_dead_energy_record = std.ascii.eqlIgnoreCase(first_after_root_morphology, "standing_dead_energy");
    if (has_geospatial_grid and !has_standing_dead_energy_record)
        return error.MissingStandingDeadEnergyRecord;
    const standing_dead_sapwood_thickness_m = if (has_standing_dead_energy_record) try nextFloat(&tokens, "standing-dead sapwood thickness in m") else 0.0025;
    const standing_dead_dry_volume_heat_capacity_mj_per_m3_k = if (has_standing_dead_energy_record) try nextFloat(&tokens, "standing-dead dry volume heat capacity in MJ m-3 K-1") else 2.496;
    const standing_dead_emissivity = if (has_standing_dead_energy_record) try nextFloat(&tokens, "standing-dead longwave emissivity") else 0.97;
    const standing_dead_activation_heat_capacity_mj_per_m2_k = if (has_standing_dead_energy_record) try nextFloat(&tokens, "standing-dead activation heat capacity in MJ m-2 K-1") else 0.838e-3;
    const standing_dead_effective_heat_capacity_floor_mj_per_m2_k = if (has_standing_dead_energy_record) try nextFloat(&tokens, "standing-dead effective heat capacity floor in MJ m-2 K-1") else 0.838e-4;
    if (!std.math.isFinite(standing_dead_sapwood_thickness_m) or standing_dead_sapwood_thickness_m <= 0 or !std.math.isFinite(standing_dead_dry_volume_heat_capacity_mj_per_m3_k) or standing_dead_dry_volume_heat_capacity_mj_per_m3_k <= 0 or !std.math.isFinite(standing_dead_emissivity) or standing_dead_emissivity < 0 or standing_dead_emissivity > 1 or !std.math.isFinite(standing_dead_activation_heat_capacity_mj_per_m2_k) or standing_dead_activation_heat_capacity_mj_per_m2_k <= 0 or !std.math.isFinite(standing_dead_effective_heat_capacity_floor_mj_per_m2_k) or standing_dead_effective_heat_capacity_floor_mj_per_m2_k <= 0) return error.InvalidStandingDeadEnergyControls;
    const first_after_standing_dead_energy = if (has_standing_dead_energy_record) try next(&tokens, "woody-optics, shoot-controls, root-gas, or site record") else first_after_root_morphology;
    const has_woody_optics_record = std.ascii.eqlIgnoreCase(first_after_standing_dead_energy, "woody_optics");
    if (has_geospatial_grid and !has_woody_optics_record)
        return error.MissingWoodyOpticsRecord;
    const woody_optics_parameters: canopy_optics.WoodyOpticsParameters = if (has_woody_optics_record) .{
        .stalk_shortwave_albedo = try nextFloat(&tokens, "stalk shortwave albedo"),
        .stalk_par_albedo = try nextFloat(&tokens, "stalk PAR albedo"),
        .standing_dead_shortwave_albedo = try nextFloat(&tokens, "standing-dead shortwave albedo"),
        .standing_dead_par_albedo = try nextFloat(&tokens, "standing-dead PAR albedo"),
    } else canopy_optics.compatibilityWoodyOpticsParameters();
    try woody_optics_parameters.validate();
    const first_after_woody_optics = if (has_woody_optics_record) try next(&tokens, "shoot-controls, root-gas, or site record") else first_after_standing_dead_energy;
    const has_canopy_retention_record = std.ascii.eqlIgnoreCase(first_after_woody_optics, "canopy_retention");
    if (has_geospatial_grid and !has_canopy_retention_record)
        return error.MissingCanopyRetentionRecord;
    const canopy_retention_parameters: canopy_precipitation_retention.Parameters = if (has_canopy_retention_record) .{
        .surface_water_capacity_m3_per_m2_by_root_profile = .{ try nextFloat(&tokens, "nonvascular canopy water capacity in m3 m-2"), try nextFloat(&tokens, "shallow-root canopy water capacity in m3 m-2"), try nextFloat(&tokens, "intermediate-root canopy water capacity in m3 m-2"), try nextFloat(&tokens, "deep-root canopy water capacity in m3 m-2") },
        .low_sun_extinction_per_area_index = try nextFloat(&tokens, "low-sun canopy extinction coefficient"),
        .minimum_solar_angle_sine_for_radiation_shares = try nextFloat(&tokens, "minimum solar-angle sine for radiation shares"),
    } else canopy_precipitation_retention.compatibilityParameters();
    try canopy_retention_parameters.validate();
    const first_after_canopy_retention = if (has_canopy_retention_record) try next(&tokens, "shoot-controls, root-gas, or site record") else first_after_woody_optics;
    const has_shoot_controls_record = std.ascii.eqlIgnoreCase(first_after_canopy_retention, "shoot_controls");
    if (has_geospatial_grid and !has_shoot_controls_record)
        return error.MissingShootControlsRecord;
    const shoot_control_parameters: plant_initialization.ShootControlParameters = if (has_shoot_controls_record) .{
        .seconds_per_hour = try nextFloat(&tokens, "seconds per hour"),
        .co2_to_water_cuticular_resistance_ratio = try nextFloat(&tokens, "CO2-to-water cuticular resistance ratio"),
        .c3_intercellular_oxygen_umol_per_mol = try nextFloat(&tokens, "C3 intercellular oxygen in umol mol-1"),
        .c4_intercellular_oxygen_umol_per_mol = try nextFloat(&tokens, "C4 intercellular oxygen in umol mol-1"),
    } else plant_initialization.compatibilityShootControlParameters();
    try shoot_control_parameters.validate();
    const first_after_shoot_controls = if (has_shoot_controls_record) try next(&tokens, "thermal-controls, root-gas, or site record") else first_after_canopy_retention;
    const has_c4_carbon_record = std.ascii.eqlIgnoreCase(first_after_shoot_controls, "c4_carbon");
    if (has_geospatial_grid and !has_c4_carbon_record)
        return error.MissingC4CarbonRecord;
    const c4_carbon_parameters: canopy_photosynthesis.C4CarbonParameters = if (has_c4_carbon_record) .{
        .bundle_sheath_water_g_per_g_c = try nextFloat(&tokens, "bundle-sheath water in g per g C"),
        .mesophyll_water_g_per_g_c = try nextFloat(&tokens, "mesophyll water in g per g C"),
        .co2_concentration_umol_per_l_per_g_c_per_g_leaf_c = try nextFloat(&tokens, "bundle-sheath CO2 concentration conversion"),
        .decarboxylation_fraction_per_h = try nextFloat(&tokens, "bundle-sheath decarboxylation fraction per hour"),
        .co2_decarboxylation_inhibition_umol_per_l = try nextFloat(&tokens, "bundle-sheath decarboxylation inhibition in umol L-1"),
        .decarboxylated_co2_fraction = try nextFloat(&tokens, "decarboxylated CO2 fraction"),
        .leakage_g_c_per_umol_per_l_g_leaf_c_h = try nextFloat(&tokens, "bundle-sheath leakage coefficient"),
        .mesophyll_feedback_half_saturation_umol_per_l = try nextFloat(&tokens, "C4 mesophyll feedback half saturation in umol L-1"),
        .co2_compensation_umol_per_l = try nextFloat(&tokens, "C4 CO2 compensation in umol L-1"),
        .electron_requirement_umol_e_per_umol_co2 = try nextFloat(&tokens, "C4 electron requirement in umol electron per umol CO2"),
    } else canopy_photosynthesis.sourceC4CarbonParameters();
    try c4_carbon_parameters.validate();
    const first_after_c4_carbon = if (has_c4_carbon_record) try next(&tokens, "thermal-controls, root-gas, or site record") else first_after_shoot_controls;
    const has_thermal_controls_record = std.ascii.eqlIgnoreCase(first_after_c4_carbon, "thermal_controls");
    if (has_geospatial_grid and !has_thermal_controls_record)
        return error.MissingThermalControlsRecord;
    const thermal_acclimation_parameters: plant_initialization.ThermalAcclimationParameters = if (has_thermal_controls_record) .{
        .adaptation_zone_pivot = try nextFloat(&tokens, "adaptation-zone pivot"),
        .cold_zone_offset_per_zone_c = try nextFloat(&tokens, "cold-zone offset per zone in C"),
        .warm_zone_offset_per_zone_c = try nextFloat(&tokens, "warm-zone offset per zone in C"),
        .base_leafout_threshold_c = try nextFloat(&tokens, "base leafout threshold in C"),
        .base_leafoff_threshold_c = try nextFloat(&tokens, "base leafoff threshold in C"),
        .maximum_leafoff_threshold_c = try nextFloat(&tokens, "maximum leafoff threshold in C"),
        .soybean_c3_seed_set_base_c = try nextFloat(&tokens, "soybean C3 seed-set base in C"),
        .other_c3_seed_set_base_c = try nextFloat(&tokens, "other C3 seed-set base in C"),
        .c4_seed_set_base_c = try nextFloat(&tokens, "C4 seed-set base in C"),
        .seed_set_adaptation_increment_c_per_zone = try nextFloat(&tokens, "seed-set adaptation increment in C per zone"),
        .soybean_loss_fraction_per_c_h = try nextFloat(&tokens, "soybean seed loss fraction per C h"),
        .other_c3_loss_fraction_per_c_h = try nextFloat(&tokens, "other C3 seed loss fraction per C h"),
        .maize_loss_fraction_per_c_h = try nextFloat(&tokens, "maize seed loss fraction per C h"),
        .other_c4_loss_fraction_per_c_h = try nextFloat(&tokens, "other C4 seed loss fraction per C h"),
    } else plant_initialization.compatibilityThermalAcclimationParameters();
    try thermal_acclimation_parameters.validate();
    const first_after_thermal_controls = if (has_thermal_controls_record) try next(&tokens, "canopy-stress, phenology-controls, root-gas record, or site filename") else first_after_c4_carbon;
    const has_canopy_stress_record = std.ascii.eqlIgnoreCase(first_after_thermal_controls, "canopy_stress");
    if (has_geospatial_grid and !has_canopy_stress_record)
        return error.MissingCanopyStressRecord;
    const canopy_stress_parameters: canopy_biochemistry.StressParameters = if (has_canopy_stress_record) .{
        .maximum_chilling_h = try nextFloat(&tokens, "maximum canopy chilling in h"),
        .heat_accumulation_threshold_c = try nextFloat(&tokens, "canopy heat accumulation threshold in C"),
        .heat_recovery_per_h = try nextFloat(&tokens, "canopy heat recovery per h"),
        .growth_temperature = .{
            .gas_constant_j_per_mol_k = try nextFloat(&tokens, "UPTAKE growth gas constant in J mol-1 K-1"),
            .temperature_scale_k = try nextFloat(&tokens, "UPTAKE growth temperature scale in K"),
            .arrhenius_log_prefactor = try nextFloat(&tokens, "UPTAKE growth Arrhenius log prefactor"),
            .activation_energy_j_per_mol = try nextFloat(&tokens, "UPTAKE growth activation energy in J mol-1"),
            .low_temperature_inactivation_j_per_mol = try nextFloat(&tokens, "UPTAKE growth low-temperature inactivation in J mol-1"),
            .high_temperature_inactivation_j_per_mol = try nextFloat(&tokens, "UPTAKE growth high-temperature inactivation in J mol-1"),
        },
    } else canopy_biochemistry.compatibilityStressParameters();
    try canopy_stress_parameters.validate();
    const first_after_canopy_stress = if (has_canopy_stress_record) try next(&tokens, "phenology-controls, root-gas record, or site filename") else first_after_thermal_controls;
    const has_phenology_controls_record = std.ascii.eqlIgnoreCase(first_after_canopy_stress, "phenology_controls");
    if (has_geospatial_grid and !has_phenology_controls_record)
        return error.MissingPhenologyControlsRecord;
    const phenology_parameters: plant_phenology.Parameters = if (has_phenology_controls_record) .{
        .gas_constant_j_per_mol_k = try nextFloat(&tokens, "phenology gas constant in J mol-1 K-1"),
        .temperature_scale_k = try nextFloat(&tokens, "phenology temperature scale in K"),
        .arrhenius_log_prefactor = try nextFloat(&tokens, "phenology Arrhenius log prefactor"),
        .activation_energy_j_per_mol = try nextFloat(&tokens, "phenology activation energy in J mol-1"),
        .low_temperature_inactivation_j_per_mol = try nextFloat(&tokens, "phenology low-temperature inactivation energy in J mol-1"),
        .high_temperature_inactivation_j_per_mol = try nextFloat(&tokens, "phenology high-temperature inactivation energy in J mol-1"),
        .minimum_turgor_potential_mpa = try nextFloat(&tokens, "minimum canopy turgor potential in MPa"),
        .oxygen_stress_exponent = try nextFloat(&tokens, "phenology oxygen-stress exponent"),
        .vegetative_stage_duration = try nextFloat(&tokens, "normalized vegetative-stage duration"),
        .reproductive_stage_duration = try nextFloat(&tokens, "normalized reproductive-stage duration"),
        .drought_leafout_total_water_potential_mpa = try nextFloat(&tokens, "drought leafout total-water-potential threshold in MPa"),
        .nonvascular_leafoff_total_water_potential_mpa = try nextFloat(&tokens, "nonvascular leafoff total-water-potential threshold in MPa"),
        .vascular_leafoff_total_water_potential_mpa = try nextFloat(&tokens, "vascular leafoff total-water-potential threshold in MPa"),
        .maximum_photoperiod_counter_h = try nextFloat(&tokens, "maximum photoperiod counter in h"),
        .emergence_area_threshold_m2_per_plant = try nextFloat(&tokens, "emergence area threshold in m2 per plant"),
        .emergence_root_depth_margin_m = try nextFloat(&tokens, "emergence root-depth margin in m"),
    } else plant_phenology.compatibilityParameters();
    try phenology_parameters.validate();
    const first_after_phenology_controls = if (has_phenology_controls_record) try next(&tokens, "plant-pool controls, root-gas record, or site filename") else first_after_canopy_stress;
    const has_plant_pool_controls_record = std.ascii.eqlIgnoreCase(first_after_phenology_controls, "plant_pool_controls");
    if (has_geospatial_grid and !has_plant_pool_controls_record)
        return error.MissingPlantPoolControlsRecord;
    const dynamic_plant_salts = if (has_plant_pool_controls_record) try delimited_input.parseYesNo(try next(&tokens, "dynamic plant salts control")) else false;
    const plant_pool_parameters: plant_pool_aggregation.RuntimeParameters = if (has_plant_pool_controls_record) .{
        .branch_structural_presence_g_per_plant = try nextFloat(&tokens, "branch structural-presence threshold in g per plant"),
        .grain_fill_detection_g_c_per_plant = try nextFloat(&tokens, "grain-fill detection threshold in g C per plant"),
        .plant_root_structural_presence_g_per_plant = try nextFloat(&tokens, "plant/root structural-presence threshold in g per plant"),
        .feedback_carbon_concentration_minimum_g_per_g = try nextFloat(&tokens, "feedback carbon-concentration minimum in g g-1"),
        .nitrogen_inhibition_g_n_per_g_c = try nextFloat(&tokens, "nitrogen inhibition in g N g-1 C"),
        .phosphorus_inhibition_g_p_per_g_c = try nextFloat(&tokens, "phosphorus inhibition in g P g-1 C"),
    } else plant_pool_aggregation.compatibilityParameters();
    try plant_pool_parameters.validate();
    const first_after_plant_pool_controls = if (has_plant_pool_controls_record) try next(&tokens, "seed-set, root-gas record, or site filename") else first_after_phenology_controls;
    const has_seed_set_controls_record = std.ascii.eqlIgnoreCase(first_after_plant_pool_controls, "seed_set_controls");
    if (has_geospatial_grid and !has_seed_set_controls_record)
        return error.MissingSeedSetControlsRecord;
    const seed_set_parameters: canopy_photosynthesis.SeedSetParameters = if (has_seed_set_controls_record) .{
        .carbon_half_saturation_g_per_g = try nextFloat(&tokens, "seed-set carbon half saturation in g g-1"),
        .nitrogen_half_saturation_g_per_g = try nextFloat(&tokens, "seed-set nitrogen half saturation in g g-1"),
        .phosphorus_half_saturation_g_per_g = try nextFloat(&tokens, "seed-set phosphorus half saturation in g g-1"),
    } else canopy_photosynthesis.compatibilitySeedSetParameters();
    try seed_set_parameters.validate();
    const first_after_seed_set_controls = if (has_seed_set_controls_record) try next(&tokens, "root-gas record or site filename") else first_after_plant_pool_controls;
    const has_root_gas_record = std.ascii.eqlIgnoreCase(first_after_seed_set_controls, "root_gas");
    if (has_geospatial_grid and !has_root_gas_record)
        return error.MissingRootGasRecord;
    var root_gas_parameters: plant_root_gas_exchange.RuntimeParameters = if (has_root_gas_record) blk: {
        var parameters = plant_root_gas_exchange.compatibilityParameters();
        parameters.reference_temperature_k = try nextFloat(&tokens, "root gas reference temperature in K");
        parameters.oxygen_aqueous_diffusivity_m2_per_h_at_reference = try nextFloat(&tokens, "oxygen aqueous diffusivity in m2 h-1 at reference temperature");
        parameters.aqueous_diffusivity_temperature_exponent = try nextFloat(&tokens, "aqueous gas diffusivity temperature exponent");
        parameters.oxygen_solubility_at_25c = try nextFloat(&tokens, "oxygen water-to-air mass solubility at 25 C");
        parameters.oxygen_solubility_activity_coefficient = try nextFloat(&tokens, "oxygen solubility activity coefficient");
        parameters.oxygen_solubility_temperature_intercept = try nextFloat(&tokens, "oxygen solubility temperature intercept");
        parameters.oxygen_solubility_temperature_coefficient_per_c = try nextFloat(&tokens, "oxygen solubility temperature coefficient per C");
        parameters.pure_water_solute_concentration_mol_per_m3 = try nextFloat(&tokens, "pure-water solute concentration in mol m-3");
        parameters.oxygen_to_carbon_respiration_ratio_g_o_per_g_c = try nextFloat(&tokens, "respiration oxygen-to-carbon ratio in g O g-1 C");
        parameters.minimum_soil_water_film_m = try nextFloat(&tokens, "minimum soil water film in m");
        parameters.soil_water_film_scale_m = try nextFloat(&tokens, "soil water film scale in m");
        parameters.soil_water_film_log_intercept = try nextFloat(&tokens, "soil water film log intercept");
        parameters.soil_water_film_potential_exponent = try nextFloat(&tokens, "soil water film potential exponent");
        const oxygen = @intFromEnum(plant_root_gas_exchange.TransportedGas.oxygen);
        parameters.aqueous_diffusivity_m2_per_h_at_reference[oxygen] = parameters.oxygen_aqueous_diffusivity_m2_per_h_at_reference;
        parameters.water_to_air_mass_solubility_at_25c[oxygen] = parameters.oxygen_solubility_at_25c;
        parameters.solubility_activity_coefficient[oxygen] = parameters.oxygen_solubility_activity_coefficient;
        parameters.solubility_temperature_intercept[oxygen] = parameters.oxygen_solubility_temperature_intercept;
        parameters.solubility_temperature_coefficient_per_c[oxygen] = parameters.oxygen_solubility_temperature_coefficient_per_c;
        break :blk parameters;
    } else plant_root_gas_exchange.compatibilityParameters();
    const first_after_root_gas_base = if (has_root_gas_record) try next(&tokens, "root-gas-transport, root-nutrients, plant-nutrients, or site filename") else first_after_seed_set_controls;
    const has_root_gas_transport_record = std.ascii.eqlIgnoreCase(first_after_root_gas_base, "root_gas_transport");
    if (has_geospatial_grid and !has_root_gas_transport_record)
        return error.MissingRootGasTransportRecord;
    if (has_root_gas_transport_record) {
        root_gas_parameters.gaseous_diffusivity_temperature_exponent = try nextFloat(&tokens, "root gaseous diffusivity temperature exponent");
        for (&root_gas_parameters.gaseous_diffusivity_m2_per_h_at_reference) |*value| value.* = try nextFloat(&tokens, "root gaseous diffusivity in m2 h-1 at reference temperature");
        for (&root_gas_parameters.aqueous_diffusivity_m2_per_h_at_reference) |*value| value.* = try nextFloat(&tokens, "root aqueous diffusivity in m2 h-1 at reference temperature");
        for (&root_gas_parameters.water_to_air_mass_solubility_at_25c) |*value| value.* = try nextFloat(&tokens, "root gas water-to-air mass solubility at 25 C");
        for (&root_gas_parameters.solubility_activity_coefficient) |*value| value.* = try nextFloat(&tokens, "root gas solubility activity coefficient");
        for (&root_gas_parameters.solubility_temperature_intercept) |*value| value.* = try nextFloat(&tokens, "root gas solubility temperature intercept");
        for (&root_gas_parameters.solubility_temperature_coefficient_per_c) |*value| value.* = try nextFloat(&tokens, "root gas solubility temperature coefficient per C");
        const oxygen = @intFromEnum(plant_root_gas_exchange.TransportedGas.oxygen);
        root_gas_parameters.oxygen_aqueous_diffusivity_m2_per_h_at_reference = root_gas_parameters.aqueous_diffusivity_m2_per_h_at_reference[oxygen];
        root_gas_parameters.oxygen_solubility_at_25c = root_gas_parameters.water_to_air_mass_solubility_at_25c[oxygen];
        root_gas_parameters.oxygen_solubility_activity_coefficient = root_gas_parameters.solubility_activity_coefficient[oxygen];
        root_gas_parameters.oxygen_solubility_temperature_intercept = root_gas_parameters.solubility_temperature_intercept[oxygen];
        root_gas_parameters.oxygen_solubility_temperature_coefficient_per_c = root_gas_parameters.solubility_temperature_coefficient_per_c[oxygen];
    }
    try plant_root_gas_exchange.validateRuntimeParameters(root_gas_parameters);
    const first_after_root_gas = if (has_root_gas_transport_record) try next(&tokens, "root-nutrients, plant-nutrients, or site filename") else first_after_root_gas_base;
    const has_root_nutrients_record = std.ascii.eqlIgnoreCase(first_after_root_gas, "root_nutrients");
    if (has_geospatial_grid and !has_root_nutrients_record)
        return error.MissingRootNutrientsRecord;
    const root_nutrient_parameters: plant_root_nutrient_uptake.RuntimeParameters = if (has_root_nutrients_record) .{
        .reference_temperature_k = try nextFloat(&tokens, "root nutrient reference temperature in K"),
        .aqueous_diffusivity_m2_per_h_at_reference = .{
            try nextFloat(&tokens, "ammonium aqueous diffusivity in m2 h-1"),
            try nextFloat(&tokens, "nitrate aqueous diffusivity in m2 h-1"),
            try nextFloat(&tokens, "phosphate aqueous diffusivity in m2 h-1"),
        },
        .aqueous_diffusivity_temperature_exponent = try nextFloat(&tokens, "root nutrient diffusivity temperature exponent"),
        .liquid_tortuosity_coefficient = try nextFloat(&tokens, "root nutrient liquid tortuosity coefficient"),
        .nitrogen_inhibition_by_nitrogen_g_n_per_g_c = try nextFloat(&tokens, "root nitrogen-uptake nitrogen inhibition constant in g N g-1 C"),
        .nitrogen_inhibition_by_phosphorus_g_p_per_g_c = try nextFloat(&tokens, "root nitrogen-uptake phosphorus inhibition constant in g P g-1 C"),
        .phosphorus_inhibition_by_phosphorus_g_p_per_g_c = try nextFloat(&tokens, "root phosphorus-uptake phosphorus inhibition constant in g P g-1 C"),
        .phosphorus_inhibition_by_nitrogen_g_n_per_g_c = try nextFloat(&tokens, "root phosphorus-uptake nitrogen inhibition constant in g N g-1 C"),
        .minimum_population_uptake_fraction_multiplier = try nextFloat(&tokens, "minimum root population uptake fraction multiplier"),
        .phosphorus_molar_mass_g_per_mol = try nextFloat(&tokens, "phosphorus molar mass in g P mol-1"),
    } else plant_root_nutrient_uptake.compatibilityRuntimeParameters();
    try root_nutrient_parameters.validate();
    const first_after_root_nutrients = if (has_root_nutrients_record) try next(&tokens, "root-salts, root-metabolism, plant-nutrients, or site filename") else first_after_root_gas;
    const has_root_salts_record = std.ascii.eqlIgnoreCase(first_after_root_nutrients, "root_salts");
    const root_salt_parameters: plant_root_salt_exchange.Parameters = if (has_root_salts_record) .{
        .reference_temperature_k = try nextFloat(&tokens, "root salt reference temperature in K"),
        .aqueous_diffusivity_m2_per_h_at_reference = .{
            try nextFloat(&tokens, "aluminum root salt diffusivity in m2 h-1"),
            try nextFloat(&tokens, "iron root salt diffusivity in m2 h-1"),
            try nextFloat(&tokens, "calcium root salt diffusivity in m2 h-1"),
            try nextFloat(&tokens, "magnesium root salt diffusivity in m2 h-1"),
            try nextFloat(&tokens, "sodium root salt diffusivity in m2 h-1"),
            try nextFloat(&tokens, "potassium root salt diffusivity in m2 h-1"),
            try nextFloat(&tokens, "sulfate root salt diffusivity in m2 h-1"),
            try nextFloat(&tokens, "chloride root salt diffusivity in m2 h-1"),
        },
        .aqueous_diffusivity_temperature_exponent = try nextFloat(&tokens, "root salt diffusivity temperature exponent"),
        .root_concentration_inhibition_mol_per_m3 = .{
            try nextFloat(&tokens, "aluminum root concentration inhibition in mol m-3"),
            try nextFloat(&tokens, "iron root concentration inhibition in mol m-3"),
            try nextFloat(&tokens, "calcium root concentration inhibition in mol m-3"),
            try nextFloat(&tokens, "magnesium root concentration inhibition in mol m-3"),
            try nextFloat(&tokens, "sodium root concentration inhibition in mol m-3"),
            try nextFloat(&tokens, "potassium root concentration inhibition in mol m-3"),
            try nextFloat(&tokens, "sulfate root concentration inhibition in mol m-3"),
            try nextFloat(&tokens, "chloride root concentration inhibition in mol m-3"),
        },
    } else plant_root_salt_exchange.compatibilityParameters();
    try root_salt_parameters.validate();
    const first_after_root_salts = if (has_root_salts_record) try next(&tokens, "root-mycorrhizal-exchange, root-exudation, root-metabolism, plant-nutrients, or site filename") else first_after_root_nutrients;
    const has_root_mycorrhizal_exchange_record = std.ascii.eqlIgnoreCase(first_after_root_salts, "root_mycorrhizal_exchange");
    if (has_geospatial_grid and !has_root_mycorrhizal_exchange_record)
        return error.MissingRootMycorrhizalExchangeRecord;
    const root_mycorrhizal_exchange_parameters: plant_root_mycorrhizal_exchange.Parameters = if (has_root_mycorrhizal_exchange_record) .{
        .minimum_partner_water_volume_ratio = try nextFloat(&tokens, "minimum mycorrhizal partner water-volume ratio"),
        .exchange_fraction_per_h = try nextFloat(&tokens, "root-mycorrhizal exchange fraction per h"),
    } else plant_root_mycorrhizal_exchange.compatibilityParameters();
    try root_mycorrhizal_exchange_parameters.validate();
    const first_after_root_mycorrhizal_exchange = if (has_root_mycorrhizal_exchange_record) try next(&tokens, "root-exudation, root-metabolism, plant-nutrients, or site filename") else first_after_root_salts;
    const has_root_exudation_record = std.ascii.eqlIgnoreCase(first_after_root_mycorrhizal_exchange, "root_exudation");
    if (has_geospatial_grid and !has_root_exudation_record)
        return error.MissingRootExudationRecord;
    const root_exudation_parameters: plant_root_exudation.Parameters = if (has_root_exudation_record) .{
        .maximum_root_carbon_concentration_g_c_per_m3 = try nextFloat(&tokens, "maximum root mobile carbon concentration in g C m-3"),
        .root_mobile_nitrogen_exchange_fraction = try nextFloat(&tokens, "root mobile nitrogen exchange fraction"),
        .root_mobile_phosphorus_exchange_fraction = try nextFloat(&tokens, "root mobile phosphorus exchange fraction"),
        .exchange_rate_per_h = try nextFloat(&tokens, "root exudation exchange rate per h"),
    } else plant_root_exudation.compatibilityParameters();
    try root_exudation_parameters.validate();
    const first_after_root_exudation = if (has_root_exudation_record) try next(&tokens, "root-porosity, root-metabolism, plant-nutrients, or site filename") else first_after_root_mycorrhizal_exchange;
    const has_root_porosity_record = std.ascii.eqlIgnoreCase(first_after_root_exudation, "root_porosity");
    if (has_geospatial_grid and !has_root_porosity_record)
        return error.MissingRootPorosityRecord;
    const root_porosity_parameters: plant_root_porosity.Parameters = if (has_root_porosity_record) .{
        .maximum_porosity_fraction = try nextFloat(&tokens, "maximum root porosity fraction"),
        .oxygen_stress_induction_fraction_per_h = try nextFloat(&tokens, "root porosity oxygen-stress induction fraction per h"),
        .relaxation_fraction_per_h = try nextFloat(&tokens, "root porosity relaxation fraction per h"),
    } else plant_root_porosity.compatibilityParameters();
    try root_porosity_parameters.validate();
    const first_after_root_porosity = if (has_root_porosity_record) try next(&tokens, "root-metabolism, plant-nutrients, or site filename") else first_after_root_exudation;
    const has_root_metabolism_record = std.ascii.eqlIgnoreCase(first_after_root_porosity, "root_metabolism");
    if (has_geospatial_grid and !has_root_metabolism_record)
        return error.MissingRootMetabolismRecord;
    const root_metabolism_parameters: plant_root_metabolism.SecondaryRootParameters = if (has_root_metabolism_record) .{
        .maximum_substrate_respiration_fraction_per_h = try nextFloat(&tokens, "maximum root substrate respiration fraction per h"),
        .substrate_respiration_half_saturation_g_c_per_g_c = try nextFloat(&tokens, "root substrate respiration half saturation in g C g-1 C"),
        .nitrogen_feedback_half_saturation_g_n_per_g_c = try nextFloat(&tokens, "root nitrogen feedback half saturation in g N g-1 C"),
        .phosphorus_feedback_half_saturation_g_p_per_g_c = try nextFloat(&tokens, "root phosphorus feedback half saturation in g P g-1 C"),
        .maintenance_respiration_g_c_per_g_n_h = try nextFloat(&tokens, "root maintenance respiration in g C g-1 N h-1"),
        .nitrogen_assimilation_respiration_g_c_per_g_n = try nextFloat(&tokens, "nitrogen assimilation respiration in g C g-1 N"),
        .minimum_carbon_recycling_fraction = try nextFloat(&tokens, "minimum root carbon recycling fraction"),
        .responsive_carbon_recycling_fraction = try nextFloat(&tokens, "responsive root carbon recycling fraction"),
        .maximum_nitrogen_recycling_fraction = try nextFloat(&tokens, "maximum root nitrogen recycling fraction"),
        .maximum_phosphorus_recycling_fraction = try nextFloat(&tokens, "maximum root phosphorus recycling fraction"),
        .storage_exchange_fraction_per_h = try nextFloat(&tokens, "root storage exchange fraction per h"),
        .nonwoody_root_fraction_exponent = try nextFloat(&tokens, "nonwoody root fraction exponent"),
        .maintenance_gas_constant_j_per_mol_k = try nextFloat(&tokens, "root maintenance gas constant in J mol-1 K-1"),
        .maintenance_enthalpy_j_per_mol_k = try nextFloat(&tokens, "root maintenance enthalpy coefficient in J mol-1 K-1"),
        .maintenance_activation_energy_j_per_mol = try nextFloat(&tokens, "root maintenance activation energy in J mol-1"),
        .maintenance_low_temperature_inactivation_energy_j_per_mol = try nextFloat(&tokens, "root maintenance low-temperature inactivation energy in J mol-1"),
        .maintenance_normalization_log_intercept = try nextFloat(&tokens, "root maintenance normalization log intercept"),
        .maximum_maintenance_temperature_response = try nextFloat(&tokens, "maximum root maintenance temperature response"),
        .shallow_root_water_response_per_mpa = try nextFloat(&tokens, "shallow-root water response per MPa"),
        .deep_root_water_response_per_mpa = try nextFloat(&tokens, "deep-root water response per MPa"),
        .maintenance_water_response_exponent = try nextFloat(&tokens, "root maintenance water response exponent"),
        .root_penetration_reference_radius_m = try nextFloat(&tokens, "root penetration reference radius in m"),
        .acidity_half_effect_hydrogen_activity_mol_per_m3 = try nextFloat(&tokens, "root acidity half-effect hydrogen activity in mol m-3"),
        .maximum_acidity_enhancement = try nextFloat(&tokens, "maximum root acidity enhancement"),
        .shallow_primary_root_sink_multiplier = try nextFloat(&tokens, "shallow primary-root sink multiplier"),
        .intermediate_primary_root_sink_multiplier = try nextFloat(&tokens, "intermediate primary-root sink multiplier"),
        .deep_primary_root_sink_multiplier = try nextFloat(&tokens, "deep primary-root sink multiplier"),
        .deeper_primary_root_sink_multiplier = try nextFloat(&tokens, "deeper primary-root sink multiplier"),
        .annual_termination_hours_without_grain_fill = try nextFloat(&tokens, "annual termination hours without grain fill"),
        .root_protein_carbon_per_nitrogen_g_c_per_g_n = try nextFloat(&tokens, "root protein carbon per nitrogen in g C g-1 N"),
        .root_protein_carbon_per_phosphorus_g_c_per_g_p = try nextFloat(&tokens, "root protein carbon per phosphorus in g C g-1 P"),
        .nutrient_uptake_respiration_g_c_per_g_element = try nextFloat(&tokens, "root nutrient uptake respiration in g C g-1 element"),
        .evergreen_leafoff_remobilization_start_fraction = try nextFloat(&tokens, "evergreen leafoff remobilization start fraction"),
        .deciduous_leafoff_remobilization_start_fraction = try nextFloat(&tokens, "deciduous leafoff remobilization start fraction"),
        .full_senescence_duration_h = try nextFloat(&tokens, "full root senescence duration in h"),
    } else plant_root_metabolism.compatibilitySecondaryRootParameters();
    try root_metabolism_parameters.validate();
    const first_after_root_metabolism = if (has_root_metabolism_record) try next(&tokens, "organ-partition, storage-remobilization, plant-nutrients, or site filename") else first_after_root_porosity;
    const has_organ_partition_record = std.ascii.eqlIgnoreCase(first_after_root_metabolism, "organ_partition");
    if (has_geospatial_grid and !has_organ_partition_record)
        return error.MissingOrganPartitionRecord;
    const organ_partition_parameters: plant_organ_partition.Parameters = if (has_organ_partition_record) .{
        .initial_leaf_fraction = try nextFloat(&tokens, "initial leaf growth partition fraction"),
        .initial_sheath_fraction = try nextFloat(&tokens, "initial sheath growth partition fraction"),
        .minimum_leaf_fraction_by_determinacy = .{ try nextFloat(&tokens, "minimum determinate leaf partition fraction"), try nextFloat(&tokens, "minimum indeterminate leaf partition fraction") },
        .minimum_sheath_fraction_by_determinacy = .{ try nextFloat(&tokens, "minimum determinate sheath partition fraction"), try nextFloat(&tokens, "minimum indeterminate sheath partition fraction") },
        .leaf_reduction_by_turnover = .{ try nextFloat(&tokens, "turnover-0 leaf partition reduction"), try nextFloat(&tokens, "turnover-1 leaf partition reduction"), try nextFloat(&tokens, "turnover-2 leaf partition reduction"), try nextFloat(&tokens, "turnover-3 leaf partition reduction"), try nextFloat(&tokens, "turnover-4 leaf partition reduction"), try nextFloat(&tokens, "turnover-5 leaf partition reduction") },
        .sheath_reduction_by_turnover = .{ try nextFloat(&tokens, "turnover-0 sheath partition reduction"), try nextFloat(&tokens, "turnover-1 sheath partition reduction"), try nextFloat(&tokens, "turnover-2 sheath partition reduction"), try nextFloat(&tokens, "turnover-3 sheath partition reduction"), try nextFloat(&tokens, "turnover-4 sheath partition reduction"), try nextFloat(&tokens, "turnover-5 sheath partition reduction") },
        .low_reserve_carbon_per_sapwood_g_c_per_g_c = try nextFloat(&tokens, "low reserve carbon per sapwood in g C g-1 C"),
        .low_reserve_redirect_fraction = try nextFloat(&tokens, "low reserve growth redirect fraction"),
    } else plant_organ_partition.compatibilityParameters();
    try organ_partition_parameters.validate();
    const first_after_organ_partition = if (has_organ_partition_record) try next(&tokens, "shoot-metabolism, shoot-root-exchange, storage-remobilization, plant-nutrients, or site filename") else first_after_root_metabolism;
    const has_shoot_metabolism_record = std.ascii.eqlIgnoreCase(first_after_organ_partition, "shoot_metabolism");
    if (has_geospatial_grid and !has_shoot_metabolism_record)
        return error.MissingShootMetabolismRecord;
    const shoot_metabolism_parameters: shoot_growth_metabolism.Parameters = if (has_shoot_metabolism_record) .{
        .maximum_mobile_carbon_oxidation_per_h = try nextFloat(&tokens, "maximum shoot mobile-carbon oxidation fraction per h"),
        .mobile_carbon_respiration_half_saturation_g_c_per_g_c = try nextFloat(&tokens, "shoot mobile-carbon respiration half saturation in g C g-1 C"),
        .maintenance_respiration_g_c_per_g_n_h = try nextFloat(&tokens, "shoot maintenance respiration in g C g-1 N h-1"),
        .minimum_leaf_nutrient_fraction = try nextFloat(&tokens, "minimum leaf nutrient fraction"),
        .nitrogen_assimilation_respiration_g_c_per_g_n = try nextFloat(&tokens, "shoot nitrogen-assimilation respiration in g C g-1 N"),
        .fixation_respiration_credit_g_c_per_g_fixed_c = try nextFloat(&tokens, "shoot fixation respiration credit in g C g-1 fixed C"),
        .mobile_nitrogen_inhibition_g_n_per_g_c = try nextFloat(&tokens, "mobile nitrogen inhibition in g N g-1 C"),
        .mobile_phosphorus_inhibition_g_p_per_g_c = try nextFloat(&tokens, "mobile phosphorus inhibition in g P g-1 C"),
    } else shoot_growth_metabolism.compatibilityParameters();
    try shoot_metabolism_parameters.validate();
    const first_after_shoot_metabolism = if (has_shoot_metabolism_record) try next(&tokens, "shoot-root-exchange, storage-remobilization, plant-nutrients, or site filename") else first_after_organ_partition;
    const has_shoot_node_growth_record = std.ascii.eqlIgnoreCase(first_after_shoot_metabolism, "shoot_node_growth");
    if (has_geospatial_grid and !has_shoot_node_growth_record)
        return error.MissingShootNodeGrowthRecord;
    const shoot_node_growth_parameters: shoot_growth_runtime.NodeGrowthParameters = if (has_shoot_node_growth_record) .{
        .protein_per_nitrogen_g_protein_per_g_n = try nextFloat(&tokens, "shoot protein per nitrogen in g protein g-1 N"),
        .protein_per_phosphorus_g_protein_per_g_p = try nextFloat(&tokens, "shoot protein per phosphorus in g protein g-1 P"),
        .leaf_mass_exponent = try nextFloat(&tokens, "specific leaf-area mass exponent"),
        .sheath_mass_exponent = try nextFloat(&tokens, "specific sheath-length mass exponent"),
        .internode_mass_exponent = try nextFloat(&tokens, "specific internode-length mass exponent"),
        .minimum_leaf_carbon_g_c_per_m2_cell = try nextFloat(&tokens, "minimum leaf carbon in g C m-2 cell"),
        .minimum_sheath_carbon_g_c_per_m2_cell = try nextFloat(&tokens, "minimum sheath carbon in g C m-2 cell"),
        .minimum_internode_carbon_g_c_per_m2_cell = try nextFloat(&tokens, "minimum internode carbon in g C m-2 cell"),
        .leaf_nutrient_exchange_fraction = try nextFloat(&tokens, "leaf nutrient exchange fraction"),
        .branch_reserve_carbon_exchange_fraction_per_h = try nextFloat(&tokens, "branch reserve carbon exchange fraction per h"),
        .branch_reserve_nutrient_exchange_fraction_per_h = try nextFloat(&tokens, "branch reserve nutrient exchange fraction per h"),
        .minimum_grain_nutrient_fraction = try nextFloat(&tokens, "minimum grain nutrient fraction"),
        .reserve_nitrogen_half_saturation_g_n_per_g_c = try nextFloat(&tokens, "reserve nitrogen half saturation in g N g-1 C"),
        .reserve_phosphorus_half_saturation_g_p_per_g_c = try nextFloat(&tokens, "reserve phosphorus half saturation in g P g-1 C"),
        .physiological_maturity_no_fill_h = try nextFloat(&tokens, "hours without grain fill to physiological maturity"),
        .annual_leafoff_delay_h_by_phenology = .{
            try nextFloat(&tokens, "evergreen annual leafoff delay in h"),
            try nextFloat(&tokens, "cold-deciduous annual leafoff delay in h"),
            try nextFloat(&tokens, "drought-deciduous annual leafoff delay in h"),
            try nextFloat(&tokens, "combined-deciduous annual leafoff delay in h"),
        },
        .leaf_storage_exchange_fraction_per_h_by_turnover = .{
            try nextFloat(&tokens, "deciduous-herbaceous leaf storage exchange fraction per h"),
            try nextFloat(&tokens, "deciduous-woody leaf storage exchange fraction per h"),
            try nextFloat(&tokens, "needleleaf-evergreen leaf storage exchange fraction per h"),
            try nextFloat(&tokens, "broadleaf-evergreen leaf storage exchange fraction per h"),
            try nextFloat(&tokens, "semi-deciduous leaf storage exchange fraction per h"),
            try nextFloat(&tokens, "semi-evergreen leaf storage exchange fraction per h"),
        },
    } else shoot_growth_runtime.sourceNodeGrowthParameters();
    try shoot_node_growth_parameters.validate();
    const first_after_shoot_node_growth = if (has_shoot_node_growth_record) try next(&tokens, "branch-mobile-exchange, shoot-root-exchange, storage-remobilization, plant-nutrients, or site filename") else first_after_shoot_metabolism;
    const has_branch_mobile_exchange_record = std.ascii.eqlIgnoreCase(first_after_shoot_node_growth, "branch_mobile_exchange");
    if (has_geospatial_grid and !has_branch_mobile_exchange_record)
        return error.MissingBranchMobileExchangeRecord;
    const branch_mobile_exchange_parameters: shoot_growth_runtime.BranchMobileExchangeParameters = if (has_branch_mobile_exchange_record) .{
        .carbon_exchange_fraction_per_h = try nextFloat(&tokens, "interbranch mobile carbon exchange fraction per h"),
        .nutrient_exchange_fraction_per_h = try nextFloat(&tokens, "interbranch mobile nutrient exchange fraction per h"),
        .remobilization_redistribution_fraction_per_h = try nextFloat(&tokens, "main-to-lateral remobilization redistribution fraction per h"),
    } else shoot_growth_runtime.sourceBranchMobileExchangeParameters();
    try branch_mobile_exchange_parameters.validate();
    const first_after_branch_mobile_exchange = if (has_branch_mobile_exchange_record) try next(&tokens, "symbiotic-fixation, shoot-root-exchange, storage-remobilization, plant-nutrients, or site filename") else first_after_shoot_node_growth;
    const has_symbiotic_fixation_record = std.ascii.eqlIgnoreCase(first_after_branch_mobile_exchange, "symbiotic_fixation");
    if (has_geospatial_grid and !has_symbiotic_fixation_record)
        return error.MissingSymbioticFixationRecord;
    const symbiotic_fixation_parameters: plant_symbiotic_fixation.RuntimeParameters = if (has_symbiotic_fixation_record) .{
        .initial_bacterial_carbon_g_c_per_m2 = try nextFloat(&tokens, "initial symbiotic bacterial carbon in g C m-2"),
        .specific_respiration_per_h = try nextFloat(&tokens, "symbiotic specific respiration per h"),
        .specific_maintenance_g_c_per_g_n_h = try nextFloat(&tokens, "symbiotic maintenance in g C g-1 N h-1"),
        .nitrogen_fixation_yield_g_n_per_g_c = try nextFloat(&tokens, "symbiotic nitrogen fixation yield in g N g-1 C"),
        .nonstructural_nitrogen_half_saturation_g_n_per_g_c = try nextFloat(&tokens, "symbiotic mobile nitrogen half saturation in g N g-1 C"),
        .nonstructural_phosphorus_half_saturation_g_p_per_g_c = try nextFloat(&tokens, "symbiotic mobile phosphorus half saturation in g P g-1 C"),
        .excess_nitrogen_inhibition_g_n_per_g_c = try nextFloat(&tokens, "symbiotic excess nitrogen inhibition in g N g-1 C"),
        .excess_nitrogen_to_phosphorus_inhibition_g_n_per_g_p = try nextFloat(&tokens, "symbiotic nitrogen phosphorus inhibition in g N g-1 P"),
        .decomposition_rate_per_h = try nextFloat(&tokens, "symbiotic decomposition rate per h"),
        .minimum_carbon_recycling_fraction = try nextFloat(&tokens, "symbiotic minimum carbon recycling fraction"),
        .carbon_recycling_range_fraction = try nextFloat(&tokens, "symbiotic responsive carbon recycling fraction"),
        .maximum_nitrogen_recycling_fraction = try nextFloat(&tokens, "symbiotic maximum nitrogen recycling fraction"),
        .maximum_phosphorus_recycling_fraction = try nextFloat(&tokens, "symbiotic maximum phosphorus recycling fraction"),
        .host_exchange_fraction_per_h_by_fixation_type = .{
            try nextFloat(&tokens, "fixation type 1 host exchange fraction per h"), try nextFloat(&tokens, "fixation type 2 host exchange fraction per h"), try nextFloat(&tokens, "fixation type 3 host exchange fraction per h"),
            try nextFloat(&tokens, "fixation type 4 host exchange fraction per h"), try nextFloat(&tokens, "fixation type 5 host exchange fraction per h"), try nextFloat(&tokens, "fixation type 6 host exchange fraction per h"),
        },
        .decomposition_control_ratio_by_fixation_type = .{
            try nextFloat(&tokens, "fixation type 1 decomposition control ratio"), try nextFloat(&tokens, "fixation type 2 decomposition control ratio"), try nextFloat(&tokens, "fixation type 3 decomposition control ratio"),
            try nextFloat(&tokens, "fixation type 4 decomposition control ratio"), try nextFloat(&tokens, "fixation type 5 decomposition control ratio"), try nextFloat(&tokens, "fixation type 6 decomposition control ratio"),
        },
    } else plant_symbiotic_fixation.sourceRuntimeParameters();
    try symbiotic_fixation_parameters.validate();
    const first_after_symbiotic_fixation = if (has_symbiotic_fixation_record) try next(&tokens, "shoot-root-exchange, storage-remobilization, plant-nutrients, or site filename") else first_after_branch_mobile_exchange;
    const has_plant_fire_combustion_record = std.ascii.eqlIgnoreCase(first_after_symbiotic_fixation, "plant_fire_combustion");
    if (has_geospatial_grid and !has_plant_fire_combustion_record)
        return error.MissingPlantFireCombustionRecord;
    const plant_fire_combustion_parameters: plant_root_disturbance.CombustionParameters = if (has_plant_fire_combustion_record) .{
        .minimum_combustion_temperature_k = try nextFloat(&tokens, "minimum root nodule combustion temperature in K"),
        .maximum_arrhenius_response = try nextFloat(&tokens, "maximum root nodule combustion Arrhenius response"),
        .gas_constant_j_per_mol_k = try nextFloat(&tokens, "root nodule combustion gas constant in J mol-1 K-1"),
        .activation_energy_j_per_mol = try nextFloat(&tokens, "root nodule combustion activation energy in J mol-1"),
        .arrhenius_intercept = try nextFloat(&tokens, "root nodule combustion Arrhenius intercept"),
        .mobile_and_leaf_specific_combustion_g_c_per_m2_h = try nextFloat(&tokens, "mobile and leaf specific combustion in g C m-2 h-1"),
        .nonwoody_structural_specific_combustion_g_c_per_m2_h = try nextFloat(&tokens, "nonwoody structural specific combustion in g C m-2 h-1"),
        .root_structural_specific_combustion_g_c_per_m2_h = try nextFloat(&tokens, "root structural specific combustion in g C m-2 h-1"),
        .woody_structural_specific_combustion_g_c_per_m2_h = try nextFloat(&tokens, "woody structural specific combustion in g C m-2 h-1"),
        .standing_dead_specific_combustion_g_c_per_m2_h = try nextFloat(&tokens, "standing-dead specific combustion in g C m-2 h-1"),
        .charcoal_activation_energy_j_per_mol = try nextFloat(&tokens, "charcoal combustion activation energy in J mol-1"),
        .charcoal_arrhenius_intercept = try nextFloat(&tokens, "charcoal combustion Arrhenius intercept"),
        .charcoal_specific_combustion_g_c_per_m2_h = try nextFloat(&tokens, "charcoal specific combustion in g C m-2 h-1"),
        .oxygen_g_per_g_combusted_carbon = try nextFloat(&tokens, "oxygen consumption in g O per g combusted C"),
        .maximum_aerobic_charcoal_fraction = try nextFloat(&tokens, "maximum aerobic charcoal fraction"),
        .maximum_anaerobic_charcoal_fraction = try nextFloat(&tokens, "maximum anaerobic charcoal fraction"),
        .oxygen_half_saturation_umol_per_mol = try nextFloat(&tokens, "fire oxygen half saturation in umol mol-1"),
        .methane_half_saturation_umol_per_mol = try nextFloat(&tokens, "fire methane half saturation in umol mol-1"),
        .aerobic_combustion_energy_mj_per_g_c = try nextFloat(&tokens, "aerobic combustion energy in MJ per g C"),
        .anaerobic_combustion_energy_mj_per_g_c = try nextFloat(&tokens, "anaerobic combustion energy in MJ per g C"),
        .methane_combustion_energy_mj_per_g_c = try nextFloat(&tokens, "methane combustion energy in MJ per g C"),
    } else plant_root_disturbance.sourceCombustionParameters();
    try plant_fire_combustion_parameters.validate();
    const first_after_plant_fire_combustion = if (has_plant_fire_combustion_record) try next(&tokens, "shoot-root-exchange, storage-remobilization, plant-nutrients, or site filename") else first_after_symbiotic_fixation;
    const has_soil_fire_combustion_record = std.ascii.eqlIgnoreCase(first_after_plant_fire_combustion, "soil_fire_combustion");
    if (has_geospatial_grid and !has_soil_fire_combustion_record)
        return error.MissingSoilFireCombustionRecord;
    const soil_fire_combustion_parameters: soil_combustion.Parameters = if (has_soil_fire_combustion_record) .{
        .minimum_combustion_temperature_k = try nextFloat(&tokens, "minimum soil combustion temperature in K"),
        .maximum_arrhenius_response = try nextFloat(&tokens, "maximum soil combustion Arrhenius response"),
        .gas_constant_j_per_mol_k = try nextFloat(&tokens, "soil combustion gas constant in J mol-1 K-1"),
        .activation_energy_j_per_mol = try nextFloat(&tokens, "soil combustion activation energy in J mol-1"),
        .arrhenius_intercept = try nextFloat(&tokens, "soil combustion Arrhenius intercept"),
        .charcoal_activation_energy_j_per_mol = try nextFloat(&tokens, "soil charcoal activation energy in J mol-1"),
        .charcoal_arrhenius_intercept = try nextFloat(&tokens, "soil charcoal Arrhenius intercept"),
        .specific_combustion_by_substrate_g_c_per_m2_h = .{
            try nextFloat(&tokens, "woody residue soil combustion in g C m-2 h-1"),
            try nextFloat(&tokens, "nonwoody residue soil combustion in g C m-2 h-1"),
            try nextFloat(&tokens, "manure soil combustion in g C m-2 h-1"),
            try nextFloat(&tokens, "particulate organic matter combustion in g C m-2 h-1"),
            try nextFloat(&tokens, "humus combustion in g C m-2 h-1"),
            try nextFloat(&tokens, "autotrophic biomass combustion in g C m-2 h-1"),
        },
        .charcoal_specific_combustion_g_c_per_m2_h = try nextFloat(&tokens, "soil charcoal combustion in g C m-2 h-1"),
        .oxygen_g_per_g_combusted_carbon = try nextFloat(&tokens, "subsurface fire oxygen in g O per g C"),
        .maximum_aerobic_charcoal_fraction = try nextFloat(&tokens, "subsurface maximum aerobic charcoal fraction"),
        .maximum_anaerobic_charcoal_fraction = try nextFloat(&tokens, "subsurface maximum anaerobic charcoal fraction"),
        .oxygen_half_saturation_g_o_per_m3 = try nextFloat(&tokens, "subsurface fire oxygen half saturation in g O m-3"),
        .methane_half_saturation_g_c_per_m3 = try nextFloat(&tokens, "subsurface fire methane half saturation in g C m-3"),
        .aerobic_combustion_energy_mj_per_g_c = try nextFloat(&tokens, "subsurface aerobic combustion energy in MJ per g C"),
        .anaerobic_combustion_energy_mj_per_g_c = try nextFloat(&tokens, "subsurface anaerobic combustion energy in MJ per g C"),
        .methane_combustion_energy_mj_per_g_c = try nextFloat(&tokens, "subsurface methane combustion energy in MJ per g C"),
    } else .{};
    try soil_fire_combustion_parameters.validate();
    const first_after_soil_fire_combustion = if (has_soil_fire_combustion_record) try next(&tokens, "shoot-root-exchange, storage-remobilization, plant-nutrients, or site filename") else first_after_plant_fire_combustion;
    const has_shoot_root_exchange_record = std.ascii.eqlIgnoreCase(first_after_soil_fire_combustion, "shoot_root_exchange");
    if (has_geospatial_grid and !has_shoot_root_exchange_record)
        return error.MissingShootRootExchangeRecord;
    const shoot_root_exchange_parameters: plant_shoot_root_exchange.Parameters = if (has_shoot_root_exchange_record) .{
        .minimum_partner_structural_ratio = try nextFloat(&tokens, "minimum shoot-root partner structural ratio"),
        .minimum_annual_carbon_exchange_fraction_per_h = try nextFloat(&tokens, "minimum annual shoot-root carbon exchange fraction per h"),
        .annual_leaf_partition_exponent = try nextFloat(&tokens, "annual leaf partition exchange exponent"),
        .salt_exchange_fraction_per_h = try nextFloat(&tokens, "shoot-root salt exchange fraction per h"),
    } else plant_shoot_root_exchange.compatibilityParameters();
    try shoot_root_exchange_parameters.validate();
    const first_after_shoot_root_exchange = if (has_shoot_root_exchange_record) try next(&tokens, "storage-remobilization, plant-nutrients, or site filename") else first_after_soil_fire_combustion;
    const has_storage_remobilization_record = std.ascii.eqlIgnoreCase(first_after_shoot_root_exchange, "storage_remobilization");
    if (has_geospatial_grid and !has_storage_remobilization_record)
        return error.MissingStorageRemobilizationRecord;
    const storage_remobilization_parameters: plant_storage_remobilization.Parameters = if (has_storage_remobilization_record) .{
        .remobilization_duration_h = .{ try nextFloat(&tokens, "annual storage remobilization duration in h"), try nextFloat(&tokens, "perennial storage remobilization duration in h") },
        .storage_carbon_oxidation_fraction_per_h = .{ try nextFloat(&tokens, "annual storage carbon oxidation fraction per h"), try nextFloat(&tokens, "perennial storage carbon oxidation fraction per h") },
        .shoot_carbon_partition_fraction = .{ try nextFloat(&tokens, "annual shoot storage partition fraction"), try nextFloat(&tokens, "perennial shoot storage partition fraction") },
        .root_carbon_partition_fraction = .{ try nextFloat(&tokens, "annual root storage partition fraction"), try nextFloat(&tokens, "perennial root storage partition fraction") },
        .perennial_nutrient_equilibration_fraction_per_h = .{
            try nextFloat(&tokens, "turnover-0 storage nutrient equilibration fraction per h"),
            try nextFloat(&tokens, "turnover-1 storage nutrient equilibration fraction per h"),
            try nextFloat(&tokens, "turnover-2 storage nutrient equilibration fraction per h"),
            try nextFloat(&tokens, "turnover-3 storage nutrient equilibration fraction per h"),
            try nextFloat(&tokens, "turnover-4 storage nutrient equilibration fraction per h"),
            try nextFloat(&tokens, "turnover-5 storage nutrient equilibration fraction per h"),
        },
        .minimum_mobile_nitrogen_per_carbon_g_n_per_g_c = try nextFloat(&tokens, "minimum mobile nitrogen per carbon in g N g-1 C"),
        .maximum_mobile_nitrogen_per_carbon_g_n_per_g_c = try nextFloat(&tokens, "maximum mobile nitrogen per carbon in g N g-1 C"),
        .minimum_mobile_phosphorus_per_carbon_g_p_per_g_c = try nextFloat(&tokens, "minimum mobile phosphorus per carbon in g P g-1 C"),
        .maximum_mobile_phosphorus_per_carbon_g_p_per_g_c = try nextFloat(&tokens, "maximum mobile phosphorus per carbon in g P g-1 C"),
        .depleted_storage_threshold_g_c_per_g_root_c = try nextFloat(&tokens, "depleted seasonal storage threshold in g C g-1 root C"),
    } else plant_storage_remobilization.compatibilityParameters();
    try storage_remobilization_parameters.validate();
    const first_after_storage_remobilization = if (has_storage_remobilization_record) try next(&tokens, "plant-nutrients record or site filename") else first_after_shoot_root_exchange;
    const has_plant_nutrients_record = std.ascii.eqlIgnoreCase(first_after_storage_remobilization, "plant_nutrients");
    if (has_geospatial_grid and !has_plant_nutrients_record)
        return error.MissingPlantNutrientsRecord;
    const plant_nutrient_initialization: soil_plant_available_nutrients.InitializationParameters = if (has_plant_nutrients_record) .{
        .initial_ammonium_band_fraction = try nextFloat(&tokens, "initial ammonium band fraction"),
        .initial_nitrate_band_fraction = try nextFloat(&tokens, "initial nitrate band fraction"),
        .initial_phosphate_band_fraction = try nextFloat(&tokens, "initial phosphate band fraction"),
        .initial_h2po4_fraction = try nextFloat(&tokens, "initial H2PO4 fraction of soluble phosphate"),
        .initial_ammonium_band_row_spacing_m = try nextFloat(&tokens, "initial ammonium band row spacing in m"),
        .initial_nitrate_band_row_spacing_m = try nextFloat(&tokens, "initial nitrate band row spacing in m"),
        .initial_phosphate_band_row_spacing_m = try nextFloat(&tokens, "initial phosphate band row spacing in m"),
    } else soil_plant_available_nutrients.compatibilityInitializationParameters();
    try plant_nutrient_initialization.validate();
    const first_after_plant_nutrients = if (has_plant_nutrients_record) try next(&tokens, "microbial-dimensions, organic-initialization record, or site filename") else first_after_storage_remobilization;
    const has_microbial_dimensions_record = std.ascii.eqlIgnoreCase(first_after_plant_nutrients, "microbial_dimensions");
    if (has_geospatial_grid and !has_microbial_dimensions_record)
        return error.MissingMicrobialDimensionsRecord;
    const microbial_substrate_count = if (has_microbial_dimensions_record) try nextUnsigned(&tokens, "microbial substrate count") else 6;
    const microbial_population_count = if (has_microbial_dimensions_record) try nextUnsigned(&tokens, "microbial population count") else 7;
    if (microbial_substrate_count == 0 or microbial_population_count == 0) return error.InvalidMicrobialDimensions;
    const first_after_microbial_dimensions = if (has_microbial_dimensions_record) try next(&tokens, "organic-initialization record or site filename") else first_after_plant_nutrients;
    const has_organic_initialization_record = std.ascii.eqlIgnoreCase(first_after_microbial_dimensions, "organic_initialization_file");
    if (has_geospatial_grid and !has_organic_initialization_record)
        return error.MissingOrganicInitializationFileRecord;
    const organic_initialization_file = if (has_organic_initialization_record) try allocator.dupe(u8, try next(&tokens, "organic initialization parameter filename")) else null;
    errdefer if (organic_initialization_file) |name| allocator.free(name);
    const first_after_organic_initialization = if (has_organic_initialization_record) try next(&tokens, "surface-gas or chemistry-initialization record or site filename") else first_after_microbial_dimensions;
    const has_surface_gas_parameter_record = std.ascii.eqlIgnoreCase(first_after_organic_initialization, "surface_gas_parameter_file");
    if (has_geospatial_grid and !has_surface_gas_parameter_record)
        return error.MissingSurfaceGasParameterFileRecord;
    const surface_gas_parameter_file = if (has_surface_gas_parameter_record) try allocator.dupe(u8, try next(&tokens, "surface gas parameter filename")) else null;
    errdefer if (surface_gas_parameter_file) |name| allocator.free(name);
    const first_after_surface_gas = if (has_surface_gas_parameter_record) try next(&tokens, "soil-nitrogen, chemistry-initialization record, or site filename") else first_after_organic_initialization;
    const has_soil_nitrogen_parameter_record = std.ascii.eqlIgnoreCase(first_after_surface_gas, "soil_nitrogen_parameter_file");
    if (has_geospatial_grid and !has_soil_nitrogen_parameter_record)
        return error.MissingSoilNitrogenParameterFileRecord;
    const soil_nitrogen_parameter_file = if (has_soil_nitrogen_parameter_record) try allocator.dupe(u8, try next(&tokens, "soil nitrogen parameter filename")) else null;
    errdefer if (soil_nitrogen_parameter_file) |name| allocator.free(name);
    const first_after_soil_nitrogen = if (has_soil_nitrogen_parameter_record) try next(&tokens, "chemistry-initialization record or site filename") else first_after_surface_gas;
    const has_chemistry_initialization_record = std.ascii.eqlIgnoreCase(first_after_soil_nitrogen, "chemistry_initialization");
    if (has_geospatial_grid and !has_chemistry_initialization_record)
        return error.MissingChemistryInitializationRecord;
    const chemistry_initialization: ?soil_chemistry_initialization.ProfileSolubleParameters = if (has_chemistry_initialization_record) .{
        .saturated_paste_phosphate_multiplier = try nextFloat(&tokens, "saturated-paste phosphate multiplier"),
        .water_activity_product_mol2_per_m6 = try nextFloat(&tokens, "water activity product in mol2 m-6"),
        .gibbsite_solubility_product_mol4_per_m12 = try nextFloat(&tokens, "gibbsite solubility product in mol4 m-12"),
        .ferric_hydroxide_solubility_product_mol4_per_m12 = try nextFloat(&tokens, "ferric hydroxide solubility product in mol4 m-12"),
        .phosphate_dissociation = .{
            .h3po4_to_h2po4_mol_per_m3 = try nextFloat(&tokens, "H3PO4 dissociation constant in mol m-3"),
            .h2po4_to_hpo4_mol_per_m3 = try nextFloat(&tokens, "H2PO4 dissociation constant in mol m-3"),
            .hpo4_to_po4_mol_per_m3 = try nextFloat(&tokens, "HPO4 dissociation constant in mol m-3"),
        },
    } else null;
    if (chemistry_initialization) |parameters| {
        inline for (.{ parameters.saturated_paste_phosphate_multiplier, parameters.water_activity_product_mol2_per_m6, parameters.gibbsite_solubility_product_mol4_per_m12, parameters.ferric_hydroxide_solubility_product_mol4_per_m12, parameters.phosphate_dissociation.h3po4_to_h2po4_mol_per_m3, parameters.phosphate_dissociation.h2po4_to_hpo4_mol_per_m3, parameters.phosphate_dissociation.hpo4_to_po4_mol_per_m3 }) |value| if (!std.math.isFinite(value) or value <= 0) return error.InvalidChemistryInitializationParameters;
        if (parameters.saturated_paste_phosphate_multiplier > 1) return error.InvalidChemistryInitializationParameters;
    }
    const first_after_chemistry_initialization = if (has_chemistry_initialization_record) try next(&tokens, "chemistry-units record or site filename") else first_after_soil_nitrogen;
    const has_chemistry_units_record = std.ascii.eqlIgnoreCase(first_after_chemistry_initialization, "chemistry_units");
    if (has_geospatial_grid and !has_chemistry_units_record)
        return error.MissingChemistryUnitsRecord;
    if (has_chemistry_units_record and chemistry_initialization == null) return error.ChemistryUnitsRequireChemistryInitialization;
    const chemistry_primary_initialization: ?soil_chemistry_initialization.PrimaryInitializationParameters = if (has_chemistry_units_record) .{
        .soluble = chemistry_initialization.?,
        .molar_mass_g_per_mol = .{
            .nitrogen = try nextFloat(&tokens, "nitrogen molar mass in g mol-1"),
            .phosphorus = try nextFloat(&tokens, "phosphorus molar mass in g mol-1"),
            .aluminum = try nextFloat(&tokens, "aluminum molar mass in g mol-1"),
            .iron = try nextFloat(&tokens, "iron molar mass in g mol-1"),
            .calcium = try nextFloat(&tokens, "calcium molar mass in g mol-1"),
            .magnesium = try nextFloat(&tokens, "magnesium molar mass in g mol-1"),
            .sodium = try nextFloat(&tokens, "sodium molar mass in g mol-1"),
            .potassium = try nextFloat(&tokens, "potassium molar mass in g mol-1"),
            .sulfur = try nextFloat(&tokens, "sulfur molar mass in g mol-1"),
            .chloride = try nextFloat(&tokens, "chloride molar mass in g mol-1"),
        },
        .minimum_ammonium_g_n_per_megagram = try nextFloat(&tokens, "minimum ammonium in g N Mg-1"),
        .minimum_calcium_g_per_megagram = try nextFloat(&tokens, "minimum calcium in g Mg-1"),
        .soil_ammonium_extract_multiplier = try nextFloat(&tokens, "soil ammonium extract multiplier"),
        .extract_mol_per_megagram_to_mol_per_m3 = try nextFloat(&tokens, "extract mol Mg-1 to mol m-3 multiplier"),
    } else null;
    if (chemistry_primary_initialization) |parameters| {
        inline for (@typeInfo(soil_chemistry_initialization.ElementMolarMassesGPerMol).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters.molar_mass_g_per_mol, field.name)) or @field(parameters.molar_mass_g_per_mol, field.name) <= 0) return error.InvalidChemistryInitializationParameters;
        inline for (.{ parameters.minimum_ammonium_g_n_per_megagram, parameters.minimum_calcium_g_per_megagram, parameters.soil_ammonium_extract_multiplier, parameters.extract_mol_per_megagram_to_mol_per_m3 }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidChemistryInitializationParameters;
        if (parameters.soil_ammonium_extract_multiplier > 1 or parameters.extract_mol_per_megagram_to_mol_per_m3 <= 0) return error.InvalidChemistryInitializationParameters;
    }
    const first_after_chemistry_units = if (has_chemistry_units_record) try next(&tokens, "chemistry-reaction-file record or site filename") else first_after_chemistry_initialization;
    const has_chemistry_reaction_file_record = std.ascii.eqlIgnoreCase(first_after_chemistry_units, "chemistry_reaction_file");
    if (has_geospatial_grid and !has_chemistry_reaction_file_record)
        return error.MissingChemistryReactionFileRecord;
    if (has_chemistry_reaction_file_record and chemistry_primary_initialization == null) return error.ChemistryReactionFileRequiresChemistryUnits;
    const chemistry_reaction_file = if (has_chemistry_reaction_file_record) try allocator.dupe(u8, try next(&tokens, "chemistry reaction parameter filename")) else null;
    errdefer if (chemistry_reaction_file) |name| allocator.free(name);
    const first_after_chemistry_reaction = if (has_chemistry_reaction_file_record) try next(&tokens, "fertilizer-units record or site filename") else first_after_chemistry_units;
    const has_fertilizer_units_record = std.ascii.eqlIgnoreCase(first_after_chemistry_reaction, "fertilizer_units");
    if (has_geospatial_grid and !has_fertilizer_units_record)
        return error.MissingFertilizerUnitsRecord;
    const fertilizer_nitrogen_molar_mass_g_per_mol = if (has_fertilizer_units_record)
        try nextFloat(&tokens, "fertilizer nitrogen molar mass in g mol-1")
    else if (chemistry_primary_initialization) |parameters|
        parameters.molar_mass_g_per_mol.nitrogen
    else
        14.0;
    if (!std.math.isFinite(fertilizer_nitrogen_molar_mass_g_per_mol) or fertilizer_nitrogen_molar_mass_g_per_mol <= 0) return error.InvalidFertilizerUnits;
    const grid_inputs_record = if (has_fertilizer_units_record)
        try next(&tokens, "grid_inputs record")
    else
        first_after_chemistry_reaction;
    if (!std.ascii.eqlIgnoreCase(grid_inputs_record, "grid_inputs"))
        return error.MissingGridInputsRecord;
    const grid_input_file = try allocator.dupe(
        u8,
        try next(&tokens, "grid-input mapping filename"),
    );
    errdefer allocator.free(grid_input_file);

    var scenarios: std.ArrayList(Scenario) = .empty;
    defer scenarios.deinit(allocator);
    var scenes: std.ArrayList(SceneFiles) = .empty;
    defer {
        for (scenes.items) |scene| freeScene(allocator, scene);
        scenes.deinit(allocator);
    }

    var execution_repeat_count: ?usize = null;
    while (true) {
        const scenario_count = try nextUnsigned(&tokens, "scenario count or 0 0 terminator");
        const group_repeat_count = try nextUnsigned(&tokens, "execution repeat count");
        if (scenario_count == 0 and group_repeat_count == 0) break;
        if (scenario_count == 0 or group_repeat_count == 0) return error.InvalidScenarioGroup;
        if (execution_repeat_count == null) execution_repeat_count = group_repeat_count;

        var scenario_index: usize = 0;
        while (scenario_index < scenario_count) : (scenario_index += 1) {
            const scene_count = try nextUnsigned(&tokens, "scene count");
            const repeat_count = try nextUnsigned(&tokens, "scenario repeat count");
            if (scene_count == 0 or repeat_count == 0) return error.EmptyScenario;
            const first_scene_index = scenes.items.len;
            var scene_index: usize = 0;
            while (scene_index < scene_count) : (scene_index += 1) {
                try scenes.append(allocator, try parseScene(allocator, &tokens));
            }
            try scenarios.append(allocator, .{
                .first_scene_index = first_scene_index,
                .scene_count = scene_count,
                .repeat_count = repeat_count,
            });
        }
    }
    if (tokens.next()) |trailing| {
        // Some supplied command-file variants retain the interactive Fortran
        // STOP record before the shell heredoc marker.
        if (!std.ascii.eqlIgnoreCase(trailing, "stop")) return error.TrailingRunscriptData;
        if (tokens.next() != null) return error.TrailingRunscriptData;
    }
    if (scenarios.items.len == 0) return error.NoScenarios;

    const owned_scenarios = try scenarios.toOwnedSlice(allocator);
    errdefer allocator.free(owned_scenarios);
    const owned_scenes = try scenes.toOwnedSlice(allocator);
    return .{
        .allocator = allocator,
        .domain = domain,
        .geospatial_bounds = geospatial_bounds,
        .tile_row_count = tile_row_count,
        .tile_column_count = tile_column_count,
        .lateral_flow_halo_cell_count = lateral_flow_halo_cell_count,
        .grid_input_file = grid_input_file,
        .plant_species_count = plant_species_count,
        .uses_four_value_species_default = uses_four_value_species_default,
        .worker_count = worker_count,
        .tile_cell_count = tile_cell_count,
        .relative_tolerance = relative_tolerance,
        .absolute_tolerance = absolute_tolerance,
        .max_nonlinear_iterations = max_nonlinear_iterations,
        .picard_relaxation = picard_relaxation,
        .uses_compatibility_runtime_controls = !has_runtime_record,
        .soil_solver_parameters = soil_solver_parameters,
        .uses_compatibility_soil_solver_parameters = !has_soil_solver_record,
        .soil_process_parameters = soil_process_parameters,
        .uses_compatibility_soil_process_parameters = !has_soil_process_record,
        .soil_gas_transport_parameters = soil_gas_transport_parameters,
        .uses_compatibility_soil_gas_transport_parameters = !has_soil_gas_transport_record,
        .soil_phase_heat_parameters = soil_phase_heat_parameters,
        .uses_compatibility_soil_phase_heat_parameters = !has_soil_phase_heat_record,
        .geothermal_controls = geothermal_controls,
        .uses_compatibility_geothermal_controls = !has_geothermal_record,
        .water_table_air_fraction_threshold = water_table_air_fraction_threshold,
        .active_layer_ice_fraction_threshold = active_layer_ice_fraction_threshold,
        .uses_compatibility_subsurface_state_controls = !has_subsurface_state_record,
        .soil_geometry_parameters = soil_geometry_parameters,
        .uses_compatibility_soil_geometry_parameters = !has_soil_geometry_record,
        .surface_pond_dry_organic_heat_capacity_mj_per_g_c_k = surface_pond_dry_organic_heat_capacity_mj_per_g_c_k,
        .surface_pond_activation_heat_capacity_mj_per_m2_k = surface_pond_activation_heat_capacity_mj_per_m2_k,
        .uses_compatibility_surface_pond_energy = !has_surface_pond_energy_record,
        .soil_longwave_emissivity = soil_longwave_emissivity,
        .snow_longwave_emissivity = snow_longwave_emissivity,
        .canopy_longwave_emissivity = canopy_longwave_emissivity,
        .snow_full_cover_depth_m = snow_full_cover_depth_m,
        .surface_sensible_heat_conductance_mj_per_m2_h_k = surface_sensible_heat_conductance_mj_per_m2_h_k,
        .surface_latent_heat_conductance_mj_per_m2_h_kpa = surface_latent_heat_conductance_mj_per_m2_h_kpa,
        .surface_vapor_activity_fraction = surface_vapor_activity_fraction,
        .minimum_surface_temperature_k = minimum_surface_temperature_k,
        .maximum_surface_temperature_k = maximum_surface_temperature_k,
        .uses_compatibility_surface_energy = !has_surface_energy_record,
        .snow_layer_bottom_depth_m = snow_layer_bottom_depth_m,
        .initial_snow_density_Mg_per_m3 = initial_snow_density_Mg_per_m3,
        .snow_ice_density_Mg_per_m3 = snow_ice_density_Mg_per_m3,
        .snow_latent_heat_of_fusion_mj_per_m3 = snow_latent_heat_of_fusion_mj_per_m3,
        .snow_phase_damping_divisor = snow_phase_damping_divisor,
        .uses_compatibility_snow_layers = !has_snow_layers_record,
        .snow_compaction_parameters = snow_compaction_parameters,
        .uses_compatibility_snow_compaction = !has_snow_compaction_record,
        .snow_heat_conduction_parameters = snow_heat_conduction_parameters,
        .uses_compatibility_snow_heat_conduction = !has_snow_thermal_record,
        .snow_vapor_parameters = snow_vapor_parameters,
        .uses_compatibility_snow_vapor = !has_snow_vapor_record,
        .snow_vapor_diffusion_parameters = snow_vapor_diffusion_parameters,
        .uses_compatibility_snow_vapor_diffusion = !has_snow_vapor_transport_record,
        .surface_gas_resistance_parameters = surface_gas_resistance_parameters,
        .uses_compatibility_surface_gas_resistance = !has_surface_gas_resistance_record,
        .surface_runoff_parameters = surface_runoff_parameters,
        .uses_compatibility_surface_runoff = !has_surface_runoff_record,
        .rainfall_impact_parameters = rainfall_impact_parameters,
        .uses_compatibility_rainfall_impact = !has_rainfall_impact_record,
        .surface_aerodynamic_parameters = surface_aerodynamic_parameters,
        .uses_compatibility_surface_aerodynamics = !has_surface_aerodynamics_record,
        .ground_air_parameters = ground_air_parameters,
        .uses_compatibility_ground_air = !has_ground_air_record,
        .canopy_surface_exchange_parameters = canopy_surface_exchange_parameters,
        .canopy_sensible_surface_resistance_h_per_m = canopy_sensible_surface_resistance_h_per_m,
        .canopy_latent_surface_resistance_h_per_m = canopy_latent_surface_resistance_h_per_m,
        .uses_compatibility_canopy_surface_exchange = !has_canopy_surface_exchange_record,
        .canopy_ammonia_exchange_parameters = canopy_ammonia_exchange_parameters,
        .uses_compatibility_canopy_ammonia_exchange = !has_canopy_ammonia_exchange_record,
        .root_axes_per_plant = root_axes_per_plant,
        .uses_compatibility_plant_structure = !has_plant_structure_record,
        .canopy_layer_count = canopy_layer_count,
        .uses_compatibility_canopy_layers = !has_canopy_layers_record,
        .canopy_discretization = canopy_discretization,
        .uses_compatibility_canopy_discretization = !has_canopy_discretization_record,
        .stalk_volume_m3_per_g_c = stalk_volume_m3_per_g_c,
        .uses_compatibility_canopy_geometry = !has_canopy_geometry_record,
        .standing_dead_partition_parameters = standing_dead_partition_parameters,
        .uses_compatibility_standing_dead_partition = !has_standing_dead_partition_record,
        .plant_heat_water_parameters = plant_heat_water_parameters,
        .uses_compatibility_plant_heat_water = !has_plant_heat_water_record,
        .plant_geometry_parameters = plant_geometry_parameters,
        .uses_compatibility_plant_geometry = !has_plant_geometry_record,
        .phenology_initialization_parameters = phenology_initialization_parameters,
        .uses_compatibility_phenology_initialization = !has_phenology_initialization_record,
        .root_initialization_parameters = root_initialization_parameters,
        .uses_compatibility_root_initialization = !has_root_initialization_record,
        .root_morphology_parameters = root_morphology_parameters,
        .uses_compatibility_root_morphology = !has_root_morphology_record,
        .standing_dead_sapwood_thickness_m = standing_dead_sapwood_thickness_m,
        .standing_dead_dry_volume_heat_capacity_mj_per_m3_k = standing_dead_dry_volume_heat_capacity_mj_per_m3_k,
        .standing_dead_emissivity = standing_dead_emissivity,
        .standing_dead_activation_heat_capacity_mj_per_m2_k = standing_dead_activation_heat_capacity_mj_per_m2_k,
        .standing_dead_effective_heat_capacity_floor_mj_per_m2_k = standing_dead_effective_heat_capacity_floor_mj_per_m2_k,
        .uses_compatibility_standing_dead_energy = !has_standing_dead_energy_record,
        .woody_optics_parameters = woody_optics_parameters,
        .uses_compatibility_woody_optics = !has_woody_optics_record,
        .canopy_retention_parameters = canopy_retention_parameters,
        .uses_compatibility_canopy_retention = !has_canopy_retention_record,
        .shoot_control_parameters = shoot_control_parameters,
        .uses_compatibility_shoot_control_parameters = !has_shoot_controls_record,
        .c4_carbon_parameters = c4_carbon_parameters,
        .uses_source_c4_carbon_parameters = !has_c4_carbon_record,
        .thermal_acclimation_parameters = thermal_acclimation_parameters,
        .uses_compatibility_thermal_acclimation_parameters = !has_thermal_controls_record,
        .canopy_stress_parameters = canopy_stress_parameters,
        .uses_compatibility_canopy_stress_parameters = !has_canopy_stress_record,
        .phenology_parameters = phenology_parameters,
        .uses_compatibility_phenology_parameters = !has_phenology_controls_record,
        .plant_pool_parameters = plant_pool_parameters,
        .dynamic_plant_salts = dynamic_plant_salts,
        .uses_compatibility_plant_pool_parameters = !has_plant_pool_controls_record,
        .seed_set_parameters = seed_set_parameters,
        .uses_compatibility_seed_set_parameters = !has_seed_set_controls_record,
        .root_gas_parameters = root_gas_parameters,
        .uses_compatibility_root_gas_parameters = !has_root_gas_record and !has_root_gas_transport_record,
        .root_nutrient_parameters = root_nutrient_parameters,
        .uses_compatibility_root_nutrient_parameters = !has_root_nutrients_record,
        .root_salt_parameters = root_salt_parameters,
        .uses_compatibility_root_salt_parameters = !has_root_salts_record,
        .root_mycorrhizal_exchange_parameters = root_mycorrhizal_exchange_parameters,
        .uses_compatibility_root_mycorrhizal_exchange_parameters = !has_root_mycorrhizal_exchange_record,
        .root_exudation_parameters = root_exudation_parameters,
        .uses_compatibility_root_exudation_parameters = !has_root_exudation_record,
        .root_porosity_parameters = root_porosity_parameters,
        .uses_compatibility_root_porosity_parameters = !has_root_porosity_record,
        .root_metabolism_parameters = root_metabolism_parameters,
        .uses_compatibility_root_metabolism_parameters = !has_root_metabolism_record,
        .organ_partition_parameters = organ_partition_parameters,
        .uses_compatibility_organ_partition_parameters = !has_organ_partition_record,
        .shoot_metabolism_parameters = shoot_metabolism_parameters,
        .uses_compatibility_shoot_metabolism_parameters = !has_shoot_metabolism_record,
        .shoot_node_growth_parameters = shoot_node_growth_parameters,
        .uses_source_shoot_node_growth_parameters = !has_shoot_node_growth_record,
        .branch_mobile_exchange_parameters = branch_mobile_exchange_parameters,
        .uses_source_branch_mobile_exchange_parameters = !has_branch_mobile_exchange_record,
        .symbiotic_fixation_parameters = symbiotic_fixation_parameters,
        .uses_source_symbiotic_fixation_parameters = !has_symbiotic_fixation_record,
        .plant_fire_combustion_parameters = plant_fire_combustion_parameters,
        .uses_source_plant_fire_combustion_parameters = !has_plant_fire_combustion_record,
        .soil_fire_combustion_parameters = soil_fire_combustion_parameters,
        .uses_source_soil_fire_combustion_parameters = !has_soil_fire_combustion_record,
        .shoot_root_exchange_parameters = shoot_root_exchange_parameters,
        .uses_compatibility_shoot_root_exchange_parameters = !has_shoot_root_exchange_record,
        .storage_remobilization_parameters = storage_remobilization_parameters,
        .uses_compatibility_storage_remobilization_parameters = !has_storage_remobilization_record,
        .plant_nutrient_initialization = plant_nutrient_initialization,
        .uses_compatibility_plant_nutrient_initialization = !has_plant_nutrients_record,
        .microbial_substrate_count = microbial_substrate_count,
        .microbial_population_count = microbial_population_count,
        .uses_compatibility_microbial_dimensions = !has_microbial_dimensions_record,
        .organic_initialization_file = organic_initialization_file,
        .surface_gas_parameter_file = surface_gas_parameter_file,
        .soil_nitrogen_parameter_file = soil_nitrogen_parameter_file,
        .chemistry_initialization = chemistry_initialization,
        .chemistry_primary_initialization = chemistry_primary_initialization,
        .chemistry_reaction_file = chemistry_reaction_file,
        .fertilizer_nitrogen_molar_mass_g_per_mol = fertilizer_nitrogen_molar_mass_g_per_mol,
        .uses_compatibility_fertilizer_units = !has_fertilizer_units_record and chemistry_primary_initialization == null,
        .execution_repeat_count = execution_repeat_count.?,
        .scenarios = owned_scenarios,
        .scenes = owned_scenes,
    };
}

fn parseSoilPhaseHeatParameters(tokens: anytype) !soil_process_science.RuntimeParameters {
    const result: soil_process_science.RuntimeParameters = .{
        .vapor_equilibrium = .{
            .vapor_density_temperature_coefficient = try nextFloat(tokens, "vapor density temperature coefficient"),
            .molecular_weight_ratio = try nextFloat(tokens, "vapor molecular-weight ratio"),
            .clausius_clapeyron_coefficient_k = try nextFloat(tokens, "Clausius-Clapeyron coefficient in K"),
            .reference_inverse_temperature_per_k = try nextFloat(tokens, "reference inverse temperature per K"),
            .water_molar_mass_g_per_mol = try nextFloat(tokens, "water molar mass in g mol-1"),
            .gas_constant_j_per_mol_k = try nextFloat(tokens, "gas constant in J mol-1 K-1"),
            .latent_heat_of_vaporization_mj_per_m3 = try nextFloat(tokens, "latent heat of vaporization in MJ m-3"),
        },
        .freeze_thaw = .{
            .freezing_potential_numerator_k_mpa = try nextFloat(tokens, "freezing potential numerator in K MPa"),
            .latent_heat_of_fusion_mj_per_m3 = try nextFloat(tokens, "latent heat of fusion in MJ m-3"),
            .ice_density_megagrams_per_m3 = try nextFloat(tokens, "ice density in Mg m-3"),
            .heat_capacity_temperature_feedback_per_k = try nextFloat(tokens, "heat-capacity temperature feedback per K"),
            .pure_water_freezing_temperature_k = try nextFloat(tokens, "pure-water freezing temperature in K"),
        },
        .heat_turbulence = .{
            .water_fraction_threshold = try nextFloat(tokens, "water convection threshold"),
            .air_fraction_threshold = try nextFloat(tokens, "air convection threshold"),
            .water_rayleigh_coefficient = try nextFloat(tokens, "water Rayleigh coefficient"),
            .air_rayleigh_coefficient = try nextFloat(tokens, "air Rayleigh coefficient"),
            .water_nusselt_denominator = try nextFloat(tokens, "water Nusselt denominator"),
            .air_nusselt_denominator = try nextFloat(tokens, "air Nusselt denominator"),
            .maximum_rayleigh_number = try nextFloat(tokens, "maximum Rayleigh number"),
        },
        .liquid_water_heat_capacity_mj_per_m3_k = try nextFloat(tokens, "liquid-water heat capacity in MJ m-3 K-1"),
        .ice_heat_capacity_mj_per_m3_k = try nextFloat(tokens, "ice heat capacity in MJ m-3 K-1"),
    };
    try soil_process_science.validate(result);
    return result;
}

fn parseSoilSolverParameters(tokens: anytype) !soil_solver_properties.RuntimeParameters {
    const result: soil_solver_properties.RuntimeParameters = .{
        .retention = .{
            .saturation_water_potential_mpa = try nextFloat(tokens, "saturation water potential in MPa"),
            .minimum_water_potential_mpa = try nextFloat(tokens, "minimum water potential in MPa"),
            .organic_soil_threshold_g_per_megagram = try nextFloat(tokens, "organic soil threshold in g Mg-1"),
            .mineral_field_capacity_intercept = try nextFloat(tokens, "mineral field-capacity intercept"),
            .mineral_field_capacity_sand_coefficient = try nextFloat(tokens, "mineral field-capacity sand coefficient"),
            .mineral_field_capacity_clay_coefficient = try nextFloat(tokens, "mineral field-capacity clay coefficient"),
            .mineral_field_capacity_organic_coefficient_per_g_per_megagram = try nextFloat(tokens, "mineral field-capacity organic coefficient per g Mg-1"),
            .mineral_wilting_point_intercept = try nextFloat(tokens, "mineral wilting-point intercept"),
            .mineral_wilting_point_clay_coefficient = try nextFloat(tokens, "mineral wilting-point clay coefficient"),
            .mineral_wilting_point_organic_coefficient_per_g_per_megagram = try nextFloat(tokens, "mineral wilting-point organic coefficient per g Mg-1"),
            .organic_bulk_density_threshold_1_megagrams_per_m3 = try nextFloat(tokens, "first organic bulk-density threshold in Mg m-3"),
            .organic_bulk_density_threshold_2_megagrams_per_m3 = try nextFloat(tokens, "second organic bulk-density threshold in Mg m-3"),
            .organic_field_capacity_1 = try nextFloat(tokens, "first organic field capacity"),
            .organic_field_capacity_2 = try nextFloat(tokens, "second organic field capacity"),
            .organic_field_capacity_3 = try nextFloat(tokens, "third organic field capacity"),
            .organic_wilting_point_1 = try nextFloat(tokens, "first organic wilting point"),
            .organic_wilting_point_2 = try nextFloat(tokens, "second organic wilting point"),
            .organic_wilting_point_3 = try nextFloat(tokens, "third organic wilting point"),
            .maximum_field_capacity_fraction_of_porosity = try nextFloat(tokens, "maximum field-capacity fraction of porosity"),
            .maximum_wilting_point_fraction_of_field_capacity = try nextFloat(tokens, "maximum wilting-point fraction of field capacity"),
            .saturation_to_field_shape = try nextFloat(tokens, "saturation-to-field retention shape"),
            .below_wilting_shape = try nextFloat(tokens, "below-wilting retention shape"),
        },
        .mualem_van_genuchten_fit_max_iterations = 80,
        .hydraulic_conductivity_class_count = try nextUnsigned(tokens, "hydraulic conductivity class count"),
        .pore_interaction_exponent = try nextFloat(tokens, "hydraulic pore-interaction exponent"),
        .air_entry_fraction_of_vertical_saturated_conductivity = try nextFloat(tokens, "air-entry conductivity fraction"),
        .mineral_saturated_conductivity_scale_m2_per_h_mpa = try nextFloat(tokens, "mineral saturated-conductivity scale in m2 h-1 MPa-1"),
        .mineral_reference_water_potential_mpa_magnitude = try nextFloat(tokens, "mineral conductivity reference water-potential magnitude in MPa"),
        .organic_saturated_conductivity_intercept_m2_per_h_mpa = try nextFloat(tokens, "organic saturated-conductivity intercept in m2 h-1 MPa-1"),
        .organic_saturated_conductivity_scale_m2_per_h_mpa = try nextFloat(tokens, "organic saturated-conductivity scale in m2 h-1 MPa-1"),
        .organic_saturated_conductivity_bulk_density_base = try nextFloat(tokens, "organic saturated-conductivity bulk-density base"),
        .profile_derivation = .{
            .particulate_carbon_fraction = try nextFloat(tokens, "missing particulate-carbon fraction"),
            .organic_nitrogen_maximum_fraction_of_carbon = try nextFloat(tokens, "maximum organic-nitrogen fraction of carbon"),
            .organic_nitrogen_scale_g_per_megagram = try nextFloat(tokens, "organic-nitrogen derivation scale in g Mg-1"),
            .organic_nitrogen_reference_carbon_g_per_megagram = try nextFloat(tokens, "organic-nitrogen reference carbon in g Mg-1"),
            .organic_nitrogen_exponent = try nextFloat(tokens, "organic-nitrogen derivation exponent"),
            .organic_phosphorus_maximum_fraction_of_carbon = try nextFloat(tokens, "maximum organic-phosphorus fraction of carbon"),
            .organic_phosphorus_scale_g_per_megagram = try nextFloat(tokens, "organic-phosphorus derivation scale in g Mg-1"),
            .organic_phosphorus_reference_carbon_g_per_megagram = try nextFloat(tokens, "organic-phosphorus reference carbon in g Mg-1"),
            .organic_phosphorus_exponent = try nextFloat(tokens, "organic-phosphorus derivation exponent"),
            .cec_conversion_mol_per_megagram_per_cmol_per_kg = try nextFloat(tokens, "CEC conversion from cmol kg-1 to mol Mg-1"),
            .organic_matter_per_carbon = try nextFloat(tokens, "organic-matter to organic-carbon ratio"),
            .organic_matter_cec_cmol_per_kg = try nextFloat(tokens, "organic-matter CEC in cmol kg-1"),
            .clay_cec_cmol_per_kg = try nextFloat(tokens, "clay CEC in cmol kg-1"),
            .silt_cec_cmol_per_kg = try nextFloat(tokens, "silt CEC in cmol kg-1"),
            .sand_cec_cmol_per_kg = try nextFloat(tokens, "sand CEC in cmol kg-1"),
            .minimum_ammonium_g_per_megagram = try nextFloat(tokens, "minimum initial ammonium in g Mg-1"),
            .minimum_calcium_g_per_megagram = try nextFloat(tokens, "minimum initial calcium in g Mg-1"),
        },
    };
    try result.retention.validate();
    try result.profile_derivation.validate();
    if (result.mualem_van_genuchten_fit_max_iterations == 0 or result.hydraulic_conductivity_class_count == 0 or !std.math.isFinite(result.pore_interaction_exponent) or result.pore_interaction_exponent <= 0 or !std.math.isFinite(result.air_entry_fraction_of_vertical_saturated_conductivity) or result.air_entry_fraction_of_vertical_saturated_conductivity < 0 or result.air_entry_fraction_of_vertical_saturated_conductivity > 1 or !std.math.isFinite(result.mineral_saturated_conductivity_scale_m2_per_h_mpa) or result.mineral_saturated_conductivity_scale_m2_per_h_mpa < 0 or !std.math.isFinite(result.mineral_reference_water_potential_mpa_magnitude) or result.mineral_reference_water_potential_mpa_magnitude <= 0 or !std.math.isFinite(result.organic_saturated_conductivity_intercept_m2_per_h_mpa) or result.organic_saturated_conductivity_intercept_m2_per_h_mpa < 0 or !std.math.isFinite(result.organic_saturated_conductivity_scale_m2_per_h_mpa) or result.organic_saturated_conductivity_scale_m2_per_h_mpa < 0 or !std.math.isFinite(result.organic_saturated_conductivity_bulk_density_base) or result.organic_saturated_conductivity_bulk_density_base <= 0) return error.InvalidSoilSolverRuntimeParameters;
    return result;
}

fn parseScene(allocator: std.mem.Allocator, tokens: anytype) !SceneFiles {
    var scene: SceneFiles = undefined;
    var allocated: usize = 0;
    errdefer {
        if (allocated >= 1) allocator.free(scene.weather_grid_file);
        if (allocated >= 2) allocator.free(scene.options);
        if (allocated >= 3) allocator.free(scene.land_management);
        if (allocated >= 4) allocator.free(scene.plant_management);
        var index: usize = 0;
        while (index + 4 < allocated) : (index += 1) allocator.free(scene.output_editors[index]);
    }
    const weather_grid_record = try next(tokens, "weather_grid record");
    if (!std.ascii.eqlIgnoreCase(weather_grid_record, "weather_grid"))
        return error.MissingWeatherGridRecord;
    scene.weather_grid_file = try allocator.dupe(
        u8,
        try next(tokens, "weather-grid mapping filename"),
    );
    allocated += 1;
    scene.options = try allocator.dupe(u8, try next(tokens, "options filename"));
    allocated += 1;
    scene.land_management = try allocator.dupe(u8, try next(tokens, "land-management filename"));
    allocated += 1;
    scene.plant_management = try allocator.dupe(u8, try next(tokens, "plant-management filename"));
    allocated += 1;
    for (&scene.output_editors) |*name| {
        name.* = try allocator.dupe(u8, try next(tokens, "output-editor filename"));
        allocated += 1;
    }
    return scene;
}

fn freeScene(allocator: std.mem.Allocator, scene: SceneFiles) void {
    allocator.free(scene.weather_grid_file);
    allocator.free(scene.options);
    allocator.free(scene.land_management);
    allocator.free(scene.plant_management);
    for (scene.output_editors) |name| allocator.free(name);
}

fn next(tokens: anytype, comptime field: []const u8) ![]const u8 {
    return tokens.next() orelse {
        std.log.err("runscript ended while reading {s}", .{field});
        return error.UnexpectedEndOfRunscript;
    };
}

fn validatePlantNutrientsRecordArity(body: []const u8) !void {
    var records = delimited_input.records(body);
    while (records.next()) |record| {
        var fields = delimited_input.recordTokens(record);
        const tag = fields.next() orelse continue;
        if (!std.ascii.eqlIgnoreCase(tag, "plant_nutrients")) continue;
        var value_count: usize = 0;
        while (fields.next() != null) value_count += 1;
        if (value_count < 7) return error.ShortPlantNutrientsRecord;
        if (value_count > 7) return error.LongPlantNutrientsRecord;
    }
}

fn nextUnsigned(tokens: anytype, comptime field: []const u8) !usize {
    const token = try next(tokens, field);
    return std.fmt.parseUnsigned(usize, token, 10) catch {
        std.log.err("invalid unsigned integer for {s}: '{s}'", .{ field, token });
        return error.InvalidRunscriptInteger;
    };
}

fn nextUnsignedType(comptime T: type, tokens: anytype, comptime field: []const u8) !T {
    const token = try next(tokens, field);
    return std.fmt.parseUnsigned(T, token, 10) catch {
        std.log.err("invalid unsigned integer for {s}: '{s}'", .{ field, token });
        return error.InvalidRunscriptInteger;
    };
}

fn nextFloat(tokens: anytype, comptime field: []const u8) !f64 {
    const token = try next(tokens, field);
    const value = std.fmt.parseFloat(f64, token) catch {
        std.log.warn("invalid floating-point value for {s}: '{s}'", .{ field, token });
        return error.InvalidRunscriptFloat;
    };
    return value;
}

fn nextBool(tokens: anytype, comptime field: []const u8) !bool {
    const token = try next(tokens, field);
    return delimited_input.parseYesNo(token) catch {
        std.log.err("invalid YES/NO value for {s}: '{s}'", .{ field, token });
        return error.InvalidRunscriptBoolean;
    };
}

fn rejectEmptyExplicitFields(body: []const u8) !void {
    var lines = std.mem.splitScalar(u8, body, '\n');
    while (lines.next()) |line| if (hasEmptyExplicitField(line))
        return error.EmptyRunscriptRecordValue;
}

fn hasEmptyExplicitField(line: []const u8) bool {
    const content = if (std.mem.indexOfScalar(u8, line, '#')) |comment|
        line[0..comment]
    else
        line;
    const trimmed = std.mem.trim(u8, content, " \r");
    if (trimmed.len == 0) return false;

    var field_start: usize = 0;
    var saw_explicit_delimiter = false;
    for (trimmed, 0..) |byte, index| {
        if (byte != ',' and byte != '|' and byte != '\t') continue;
        if (std.mem.trim(u8, trimmed[field_start..index], " \r").len == 0)
            return true;
        field_start = index + 1;
        saw_explicit_delimiter = true;
    }
    return saw_explicit_delimiter and
        std.mem.trim(u8, trimmed[field_start..], " \r").len == 0;
}

fn extractInputBody(source: []const u8) []const u8 {
    var lines = std.mem.splitScalar(u8, source, '\n');
    var offset: usize = 0;
    while (lines.next()) |line_with_cr| {
        const line = std.mem.trim(u8, line_with_cr, " \t\r");
        if (std.mem.indexOf(u8, line, "<<") != null) {
            const marker_start = std.mem.indexOf(u8, line, "<<").? + 2;
            const marker_tail = std.mem.trim(u8, line[marker_start..], " \t\r");
            var marker_tokens = std.mem.tokenizeAny(u8, marker_tail, " \t\r");
            const marker = marker_tokens.next() orelse return source;
            const body_start = @min(offset + line_with_cr.len + 1, source.len);
            var body_lines = std.mem.splitScalar(u8, source[body_start..], '\n');
            var body_length: usize = 0;
            while (body_lines.next()) |body_line_with_cr| {
                const body_line = std.mem.trim(u8, body_line_with_cr, " \t\r");
                if (std.mem.eql(u8, body_line, marker)) return source[body_start .. body_start + body_length];
                body_length += body_line_with_cr.len + 1;
            }
            return source[body_start..];
        }
        offset += line_with_cr.len + 1;
    }
    return source;
}

test "parse shell-wrapped self-contained runscript" {
    const source = "#!/bin/sh\nmodel << eor\n1,1,1,1\ngrid_inputs grid_inputs\n1 1\n1 1\nweather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10\n0 0\neor\n";
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(usize, 1), try script.gridCellCount());
    try std.testing.expectEqual(@as(usize, 1), script.scenarios.len);
    try std.testing.expectEqual(@as(usize, 1), script.scenes.len);
    try std.testing.expectEqualStrings("weather", script.scenes[0].weather_grid_file);
    try std.testing.expectEqualStrings("o10", script.scenes[0].output_editors[9]);
}

test "parse mixed-case stop record" {
    const source = "1 1 1 1\ngrid_inputs grid_inputs\n1 1\n1 1\nweather_grid weather options no plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10\n0 0\n";
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(usize, 1), script.scenes.len);
}

test "plant nutrients record enforces seven compulsory values on one line" {
    try validatePlantNutrientsRecordArity(
        "plant_nutrients,0,0,0,1,0.75,0.8,0.9 # metres\n",
    );
    try std.testing.expectError(
        error.ShortPlantNutrientsRecord,
        validatePlantNutrientsRecordArity(
            "plant_nutrients|0|0|0|1|0.75|0.8\n",
        ),
    );
    try std.testing.expectError(
        error.LongPlantNutrientsRecord,
        validatePlantNutrientsRecordArity(
            "plant_nutrients 0 0 0 1 0.75 0.8 0.9 1.0\n",
        ),
    );
}

test "global site topography and weather positional schema is rejected" {
    const old_source = "1 1 1 1\nsite\ntopography\n1 1\n1 1\nweather options no plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10\n0 0\n";
    try std.testing.expectError(
        error.MissingGridInputsRecord,
        parse(std.testing.allocator, old_source),
    );

    const global_weather_source = "1 1 1 1\ngrid_inputs cells.txt\n1 1\n1 1\nweather options no plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10\n0 0\n";
    try std.testing.expectError(
        error.MissingWeatherGridRecord,
        parse(std.testing.allocator, global_weather_source),
    );
}

test "parse comma-delimited raw runscript" {
    const source =
        \\1,1,2,3
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(usize, 6), try script.gridCellCount());
}

test "parse pipe-delimited extensionless runscript content" {
    const source =
        \\1|1|1|1|37
        \\grid_inputs grid_inputs
        \\1|1
        \\1|1
        \\weather_grid|weather|options|NO|plants|o1|o2|o3|o4|o5|o6|o7|o8|o9|o10
        \\0|0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqualStrings("options", script.scenes[0].options);
    try std.testing.expectEqual(@as(usize, 37), script.plant_species_count);
    try std.testing.expect(!script.uses_four_value_species_default);
}

test "runscript rejects empty explicit delimiter fields before token collapse" {
    inline for (.{
        "1,1,,1,1\ngrid_inputs grid_inputs\n",
        "1|1| |1|1\ngrid_inputs grid_inputs\n",
        "1\t1\t\t1\t1\ngrid_inputs grid_inputs\n",
        "1,1,1,1, # missing species before comment\ngrid_inputs grid_inputs\n",
    }) |source| try std.testing.expectError(
        error.EmptyRunscriptRecordValue,
        parse(std.testing.allocator, source),
    );
}

test "runscript empty-field preflight handles comments and shell wrappers" {
    try rejectEmptyExplicitFields(
        "1, 1 | 1\t1 # mixed valid delimiters\n# comment only\n",
    );
    try std.testing.expectError(
        error.EmptyRunscriptRecordValue,
        parse(
            std.testing.allocator,
            "#!/bin/sh\nmodel << eor\n1,1,,1,1\neor\n",
        ),
    );
}

test "modern compact header defines horizontal cells vertical cells and species" {
    const source =
        \\99,99,12
        \\geospatial_grid,45,48,-80,-76,1,1
        \\tile_layout,2,3,2
        \\runtime,6,128,1e-9,1e-12,55,0.6
        \\soil_solver,-0.0005,-1.5e12,250000,0.2576,-0.20,0.36,0.60e-6,0.0260,0.50,0.32e-6,0.075,0.195,0.27,0.62,0.71,0.04,0.15,0.22,0.75,0.75,0.5,0.5,100,1.33,0.1,1.54,0.033,0.10,75,1e-15,0.067,0.125,890,10000,0.80,0.0125,120,10000,0.52,10,1.82,200,80,20,5,1,1
        \\soil_process,0.0098,8.96e-2,298.15,1.75,0.66,0.03,0.5e-3,1e-6,0.533,0.0267
        \\macropore_van_genuchten,0,15,2.68,0.5
        \\dual_domain_exchange,3,0.4
        \\frozen_hydraulic_impedance,7
        \\surface_residue_freeze_thaw,0,400,2.5
        \\soil_gas_transport,298.15,1.75,4.68e-2,7.80e-2,6.43e-2,5.57e-2,5.57e-2,6.67e-2,5.57e-2,0.66,1e-12,1e6,18
        \\soil_phase_heat,2.173e-3,0.61,5360,3.661e-3,18,8.3143,2450,9.0959e4,333,0.917,6.2913e-3,273.15,0.333,0.333,13990344827.5862,89223880597.0149,1.09495999584408,4.77823246449888,10000,4.19,1.9274
        \\geothermal,yes,10,1,8.1e-3,2.052e-4
        \\subsurface_state,1e-3,1e-6
        \\soil_geometry,1.82e-6,110000,0.083,1e-9
        \\surface_pond_energy,2.496e-6,0.838e-3
        \\surface_energy,0.97,0.97,0.97,0.07,0.43,0,1,173.15,373.15
        \\surface_runoff,0.005,0.011,1e-3,3600,1e-12
        \\rainfall_impact,8.95,8.44,15.8,5.87,2.5,2,1e-3,5e-4
        \\canopy_ammonia_exchange,0.16,0.10,0.05,2,1e-4,0.1,0.513,0.0171
        \\plant_structure,10
        \\canopy_layers,10
        \\canopy_geometry,4e-6
        \\canopy_discretization,4,4,4
        \\standing_dead_partition,0,0.045,0.660,0.295,0.020,0.010,0.010,0.020,0.0020,0.0010,0.0010,0.0020
        \\plant_initial_heat_water,273.15,2.173e-3,0.61,5360,3.661e-3,-1e-3
        \\plant_initial_geometry,5e-6,2,0.75,3.1416,0.33,4,1e-6,0.05,0.01,3.142
        \\plant_initial_phenology,0.005,1,24,10,15,3,4,5
        \\root_initialization,2.5,25,2.5e-6,-0.01,-0.01,1e-3,1
        \\root_morphology,1e-2,5
        \\standing_dead_energy,0.0025,2.496,0.97,0.838e-3,0.838e-4
        \\woody_optics,0.1,0.1,0.1,0.1
        \\canopy_retention,5e-4,2.5e-4,2.5e-4,2.5e-4,0.65,0.05
        \\shoot_controls,3600,1.56,2.10e5,3.96e5
        \\c4_carbon,1.2,4.8,0.083e9,0.025,1000,0.02,5e-7,5e6,0.5,3
        \\thermal_controls,3,2.5,1.25,5,12.5,15,35,27.5,30,2,0.002,0.005,0.010,0.002
        \\canopy_stress,24,60,0.02,8.3143,710,25.229,62500,197500,222500
        \\phenology_controls,8.3143,710,24.269,60000,197500,218500,0.1,0.25,2,0.667,-0.1,-150,-1.5,3600,1e-6,1e-6
        \\plant_pool_controls,no,1e-15,1e-6,1e-6,1e-6,1e-2,1e-3
        \\seed_set_controls,2.5e-2,0.5e-2,0.1e-2
        \\root_gas,298.15,8.57e-6,6,2.925e-2,0.31,0.516,0.0172,5.56e4,2.667,1e-6,0.5,-13.833,-0.857
        \\root_gas_transport,1.75,4.68e-2,7.80e-2,5.57e-2,6.67e-2,5.57e-2,6.43e-2,4.25e-6,7.08e-6,5.72e-6,4.00e-6,7.34e-6,8.57e-6,7.391e-1,3.156e-2,5.241e-1,2.852e2,3.156e-2,2.925e-2,0.14,0.14,0.23,0.07,0.14,0.31,0.843,0.597,0.897,0.513,0.597,0.516,0.0281,0.0199,0.0299,0.0171,0.0199,0.0172
        \\root_nutrients,298.15,4e-6,6e-6,3e-6,6,0.7,1,1,0.01,0.01,1e-4,31
        \\root_mycorrhizal_exchange,0.05,0.5
        \\root_exudation,1e3,0.1,0.1,1e-3
        \\root_porosity,0.75,0.1,0.01
        \\root_metabolism,0.015,0.025,0.1,0.01,0.010,1.70,0.167,0.333,0.667,0.667,2.5e-5,0.167,8.3143,710,62500,197500,25.216,1e3,0.05,0.10,0.25,1e-3,1,4,0.25,1,2,4,336,2.5,25,0.86,0.75,0.5,480
        \\organ_partition,0.75,0.25,0.0200,0.0500,0.0067,0.0167,0.75,1.50,2,2,1.75,1.50,0.25,0.50,0.67,0.67,0.58,0.50,0.10,0.10
        \\shoot_metabolism,0.015,0.025,0.010,0.333,1.70,0.025,0.1,0.01
        \\shoot_node_growth,2.5,25,-0.333,-0.50,-0.667,0.002,2,2,1e-3,5e-3,5e-2,0.75,0.005,0.001,72,360,1440,720,720,5e-3,5e-3,5e-6,5e-6,5e-5,5e-4
        \\branch_mobile_exchange,0.01,0.01,0.05
        \\symbiotic_fixation,1e-4,0.125,0.010,0.25,1e-4,1e-5,10,1000,1e-2,0.167,0.333,0.333,0.333,0.20,0.10,0.05,0.20,0.10,0.05,0.50,0.25,0.125,0.050,0.025,0.0125
        \\plant_fire_combustion,473.15,2,8.3143,60000,12.028,1000,1000,1000,1000,2500,120000,20.620,1,2.667,0,0.5,2100,10,0.0375,0.0125,0.0743
        \\soil_fire_combustion,473.15,2,8.3143,60000,12.028,120000,20.620,1000,5000,1000,1000,1000,1000,1,2.667,0,0.5,2.8,0.005,0.0375,0.0125,0.0743
        \\shoot_root_exchange,0.05,0.005,0.25,1
        \\storage_remobilization,45.8,138.4,0.015,0.005,0.25,0.25,0.75,0.75,0.100,0.100,0.010,0.100,0.100,0.100,0.050,0.20,0.005,0.020,0.10
        \\plant_nutrients,0,0,0,1,1,1,1
        \\microbial_dimensions,6,7
        \\organic_initialization_file,organic-parameters
        \\surface_gas_parameter_file,surface-gas-parameters
        \\soil_nitrogen_parameter_file,soil-nitrogen-parameters
        \\chemistry_initialization,0.01,1e-8,1.9e-21,6.3e-26,7.5,6.2e-5,4.8e-10
        \\chemistry_units,14,31,27,56,40,24.3,23,39.1,32,35.5,1,1,0.01,1
        \\chemistry_reaction_file,chemistry-reaction-parameters
        \\fertilizer_units,14
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(usize, 12), try script.gridCellCount());
    try std.testing.expectEqual(@as(usize, 12), script.plant_species_count);
    try std.testing.expectEqual(@as(usize, 2), script.tile_row_count);
    try std.testing.expectEqual(@as(usize, 3), script.tile_column_count);
    try std.testing.expectEqual(@as(usize, 2), script.lateral_flow_halo_cell_count);
    try std.testing.expectEqual(@as(usize, 6), script.worker_count);
    try std.testing.expectEqual(@as(usize, 128), script.tile_cell_count);
    try std.testing.expectEqual(@as(u16, 55), script.max_nonlinear_iterations);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), script.picard_relaxation, 1.0e-15);
    try std.testing.expect(!script.uses_compatibility_runtime_controls);
    try std.testing.expectEqual(@as(usize, 10), script.root_axes_per_plant);
    try std.testing.expect(!script.uses_compatibility_plant_structure);
    try std.testing.expect(!script.uses_compatibility_canopy_layers);
    try std.testing.expect(!script.uses_compatibility_canopy_geometry);
    try std.testing.expect(!script.uses_compatibility_canopy_discretization);
}

test "geospatial runscript requires explicit improved water and freeze controls" {
    const source =
        \\1,1,1
        \\geospatial_grid,45,46,-80,-79,1,1
        \\tile_layout,1,1,2
        \\runtime,1,1,1e-9,1e-12,40,0.5
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    try std.testing.expectError(error.MissingSoilSolverRecord, parse(std.testing.allocator, source));
}

test "geospatial runscript requires explicit soil gas transport controls" {
    const source =
        \\1,1,1
        \\geospatial_grid,45,46,-80,-79,1,1
        \\tile_layout,1,1,2
        \\runtime,1,1,1e-9,1e-12,40,0.5
        \\soil_solver,-0.0005,-1.5e12,250000,0.2576,-0.20,0.36,0.60e-6,0.0260,0.50,0.32e-6,0.075,0.195,0.27,0.62,0.71,0.04,0.15,0.22,0.75,0.75,0.5,0.5,100,1.33,0.1,1.54,0.033,0.10,75,1e-15,0.067,0.125,890,10000,0.80,0.0125,120,10000,0.52,10,1.82,200,80,20,5,1,1
        \\soil_process,0.0098,8.96e-2,298.15,1.75,0.66,0.03,0.5e-3,1e-6,0.533,0.0267
        \\macropore_van_genuchten,0,15,2.68,0.5
        \\dual_domain_exchange,3,0.4
        \\frozen_hydraulic_impedance,7
        \\surface_residue_freeze_thaw,0,400,2.5
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    try std.testing.expectError(error.MissingSoilGasTransportRecord, parse(std.testing.allocator, source));
}

test "plant root axis count is a case-insensitive runtime input" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\PlAnT_StRuCtUrE,37
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(usize, 37), script.root_axes_per_plant);
    try std.testing.expect(!script.uses_compatibility_plant_structure);
}

test "snow layer count boundaries and density are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\SnOw_LaYeRs,7,0.02,0.05,0.10,0.20,0.40,0.80,1.60,0.075,0.92,333,2.7185
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(usize, 7), script.snow_layer_bottom_depth_m.len);
    try std.testing.expectApproxEqAbs(@as(f64, 1.6), script.snow_layer_bottom_depth_m[6], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.075), script.initial_snow_density_Mg_per_m3, 1e-15);
    try std.testing.expect(!script.uses_compatibility_snow_layers);
}

test "canopy layer count is a case-insensitive runtime input" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\CaNoPy_LaYeRs|17
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(usize, 17), script.canopy_layer_count);
    try std.testing.expect(!script.uses_compatibility_canopy_layers);
}

test "canopy angular dimensions are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\CaNoPy_DiScReTiZaTiOn|7|9|11
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(usize, 7), script.canopy_discretization.leaf_inclination_class_count);
    try std.testing.expectEqual(@as(usize, 9), script.canopy_discretization.leaf_azimuth_class_count);
    try std.testing.expectEqual(@as(usize, 11), script.canopy_discretization.diffuse_sky_sector_count);
    try std.testing.expect(!script.uses_compatibility_canopy_discretization);
}

test "STARTQ standing dead partition coefficients are runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\StAnDiNg_DeAd_PaRtItIoN|0.1|0.2|0.3|0.4|1|2|3|4|5|6|7|8
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual([4]f64{ 0.1, 0.2, 0.3, 0.4 }, script.standing_dead_partition_parameters.carbon_fraction);
    try std.testing.expectEqual([4]f64{ 1, 2, 3, 4 }, script.standing_dead_partition_parameters.nitrogen_weight);
    try std.testing.expectEqual([4]f64{ 5, 6, 7, 8 }, script.standing_dead_partition_parameters.phosphorus_weight);
    try std.testing.expect(!script.uses_compatibility_standing_dead_partition);
}

test "STARTQ initial plant heat and water coefficients are runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\PlAnT_InItIaL_HeAt_WaTeR|270|0.003|0.7|5000|0.004|-0.02
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(f64, 270), script.plant_heat_water_parameters.kelvin_offset_k);
    try std.testing.expectEqual(@as(f64, -0.02), script.plant_heat_water_parameters.initial_total_water_potential_mpa);
    try std.testing.expect(!script.uses_compatibility_plant_heat_water);
}

test "STARTQ seed and root geometry coefficients are runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\PlAnT_InItIaL_GeOmEtRy|6e-6|2.1|0.8|3.14|0.34|4.1|2e-6|0.06|0.02|3.141
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(f64, 6e-6), script.plant_geometry_parameters.seed_volume_m3_per_g_c);
    try std.testing.expectEqual(@as(f64, 0.02), script.plant_geometry_parameters.root_porosity_floor);
    try std.testing.expect(!script.uses_compatibility_plant_geometry);
}

test "READQ and STARTQ initial phenology coefficients are runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\PlAnT_InItIaL_PhEnOlOgY|0.01|2|36|8|14|4|6|8
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(f64, 0.01), script.phenology_initialization_parameters.perennial_input_scale);
    try std.testing.expectEqual(@as(usize, 36), script.phenology_initialization_parameters.perennial_maximum_concurrently_growing_nodes);
    try std.testing.expectEqual(@as(usize, 8), script.phenology_initialization_parameters.late_maximum_concurrently_growing_nodes);
    try std.testing.expect(!script.uses_compatibility_phenology_initialization);
}

test "STARTQ root initialization coefficients are runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\RoOt_InItIaLiZaTiOn|3|30|4e-6|-0.02|-0.03|0.002|0.9
        \\rOoT_MoRpHoLoGy|0.015|6.5
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(f64, 4e-6), script.root_initialization_parameters.mycorrhizal_radius_m);
    try std.testing.expectEqual(@as(f64, -0.02), script.root_initialization_parameters.initial_total_water_potential_mpa);
    try std.testing.expectEqual(@as(f64, 0.9), script.root_initialization_parameters.initial_water_fraction);
    try std.testing.expect(!script.uses_compatibility_root_initialization);
    try std.testing.expectEqual(@as(f64, 0.015), script.root_morphology_parameters.minimum_average_secondary_length_m);
    try std.testing.expectEqual(@as(f64, 6.5), script.root_morphology_parameters.root_elastic_modulus_mpa);
    try std.testing.expect(!script.uses_compatibility_root_morphology);
}

test "microbial substrate and population dimensions are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\MiCrObIaL_DiMeNsIoNs|11|19
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(usize, 11), script.microbial_substrate_count);
    try std.testing.expectEqual(@as(usize, 19), script.microbial_population_count);
    try std.testing.expect(!script.uses_compatibility_microbial_dimensions);
}

test "legacy microbial dimensions remain explicit compatibility values" {
    const source =
        \\1,1,2
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqual(@as(usize, 6), script.microbial_substrate_count);
    try std.testing.expectEqual(@as(usize, 7), script.microbial_population_count);
    try std.testing.expect(script.uses_compatibility_microbial_dimensions);
}

test "STARTQ canopy geometry is a case-insensitive runtime input" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\CaNoPy_LaYeRs|17
        \\cAnOpY_gEoMeTrY|5.25e-6
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 5.25e-6), script.stalk_volume_m3_per_g_c, 1.0e-15);
    try std.testing.expect(!script.uses_compatibility_canopy_geometry);
}

test "UPTAKE standing-dead energy coefficients are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\StAnDiNg_DeAd_EnErGy|0.003|2.6|0.96|0.0009|0.00009
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 0.003), script.standing_dead_sapwood_thickness_m, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.6), script.standing_dead_dry_volume_heat_capacity_mj_per_m3_k, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.96), script.standing_dead_emissivity, 1.0e-15);
    try std.testing.expect(!script.uses_compatibility_standing_dead_energy);
}

test "HOUR1 woody optics are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\WoOdY_OpTiCs|0.12|0.08|0.2|0.15
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 0.12), script.woody_optics_parameters.stalk_shortwave_albedo, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.15), script.woody_optics_parameters.standing_dead_par_albedo, 1.0e-15);
    try std.testing.expect(!script.uses_compatibility_woody_optics);
}

test "HOUR1 canopy retention controls are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\CaNoPy_ReTeNtIoN|0.0006|0.0003|0.00035|0.0004|0.7|0.06
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 0.0006), script.canopy_retention_parameters.surface_water_capacity_m3_per_m2_by_root_profile[0], 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.7), script.canopy_retention_parameters.low_sun_extinction_per_area_index, 1.0e-15);
    try std.testing.expect(!script.uses_compatibility_canopy_retention);
}

test "STARTQ shoot controls are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\ShOoT_CoNtRoLs|3600|1.56|210000|396000
        \\C4_CaRbOn|1.25|4.75|83000000|0.026|1100|0.03|0.0000006|5100000|0.6|3.2
        \\ThErMaL_CoNtRoLs|3|2.5|1.25|5|12.5|15|35|27.5|30|2|0.002|0.005|0.01|0.002
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_shoot_control_parameters);
    try std.testing.expectEqual(@as(f64, 3600), script.shoot_control_parameters.seconds_per_hour);
    try std.testing.expectEqual(@as(f64, 1.56), script.shoot_control_parameters.co2_to_water_cuticular_resistance_ratio);
    try std.testing.expectEqual(@as(f64, 396000), script.shoot_control_parameters.c4_intercellular_oxygen_umol_per_mol);
    try std.testing.expect(!script.uses_source_c4_carbon_parameters);
    try std.testing.expectEqual(@as(f64, 1.25), script.c4_carbon_parameters.bundle_sheath_water_g_per_g_c);
    try std.testing.expectEqual(@as(f64, 0.026), script.c4_carbon_parameters.decarboxylation_fraction_per_h);
    try std.testing.expectEqual(@as(f64, 0.0000006), script.c4_carbon_parameters.leakage_g_c_per_umol_per_l_g_leaf_c_h);
    try std.testing.expectEqual(@as(f64, 5100000), script.c4_carbon_parameters.mesophyll_feedback_half_saturation_umol_per_l);
    try std.testing.expectEqual(@as(f64, 0.6), script.c4_carbon_parameters.co2_compensation_umol_per_l);
    try std.testing.expectEqual(@as(f64, 3.2), script.c4_carbon_parameters.electron_requirement_umol_e_per_umol_co2);
    try std.testing.expect(!script.uses_compatibility_thermal_acclimation_parameters);
    try std.testing.expectEqual(@as(f64, 35), script.thermal_acclimation_parameters.soybean_c3_seed_set_base_c);
}

test "root gas scientific coefficients are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\RoOt_GaS,300,9e-6,5.5,0.03,0.32,0.52,0.018,56000,2.7,2e-6,0.4,-13.7,-0.8
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_root_gas_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 300), script.root_gas_parameters.reference_temperature_k, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 9e-6), script.root_gas_parameters.oxygen_aqueous_diffusivity_m2_per_h_at_reference, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.7), script.root_gas_parameters.oxygen_to_carbon_respiration_ratio_g_o_per_g_c, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2e-6), script.root_gas_parameters.minimum_soil_water_film_m, 1.0e-15);
}

test "all transported root gas coefficients are runtime inputs" {
    const source =
        \\1 1 2
        \\runtime 1 1 1e-8 1e-11 40 0.5
        \\RoOt_GaS_TrAnSpOrT|1.8
        \\0.041|0.072|0.053|0.064|0.055|0.061
        \\4.1e-6|7.1e-6|5.8e-6|4.2e-6|7.5e-6|9.1e-6
        \\0.71|0.032|0.53|280|0.033|0.031
        \\0.15|0.16|0.24|0.08|0.17|0.32
        \\0.84|0.59|0.90|0.51|0.60|0.52
        \\0.028|0.020|0.030|0.017|0.021|0.018
        \\grid_inputs grid_inputs
        \\1 1
        \\1 1
        \\weather_grid weather options no plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0 0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_root_gas_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 1.8), script.root_gas_parameters.gaseous_diffusivity_temperature_exponent, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.041), script.root_gas_parameters.gaseous_diffusivity_m2_per_h_at_reference[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 9.1e-6), script.root_gas_parameters.aqueous_diffusivity_m2_per_h_at_reference[5], 1e-18);
    try std.testing.expectApproxEqAbs(@as(f64, 280), script.root_gas_parameters.water_to_air_mass_solubility_at_25c[3], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.018), script.root_gas_parameters.solubility_temperature_coefficient_per_c[5], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 9.1e-6), script.root_gas_parameters.oxygen_aqueous_diffusivity_m2_per_h_at_reference, 1e-18);
}

test "root nutrient diffusion coefficients are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\RoOt_NuTrIeNtS|300|4.1e-6|6.2e-6|3.3e-6|5.8|0.68|1.1|1.2|0.011|0.012|0.0002|30.97
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_root_nutrient_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 300), script.root_nutrient_parameters.reference_temperature_k, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 6.2e-6), script.root_nutrient_parameters.aqueous_diffusivity_m2_per_h_at_reference[1], 1.0e-18);
    try std.testing.expectApproxEqAbs(@as(f64, 5.8), script.root_nutrient_parameters.aqueous_diffusivity_temperature_exponent, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.68), script.root_nutrient_parameters.liquid_tortuosity_coefficient, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.012), script.root_nutrient_parameters.phosphorus_inhibition_by_nitrogen_g_n_per_g_c, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0002), script.root_nutrient_parameters.minimum_population_uptake_fraction_multiplier, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 30.97), script.root_nutrient_parameters.phosphorus_molar_mass_g_per_mol, 1.0e-15);
}

test "UPTAKE root salt coefficients are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\RoOt_SaLtS|300|1e-6|2e-6|3e-6|4e-6|5e-6|6e-6|7e-6|8e-6|5.9|1e-5|2e-5|3e-5|4e-5|5e-5|6e-5|7e-5|8e-5
        \\RoOt_MyCoRrHiZaL_ExChAnGe|0.07|0.45
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_root_salt_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 300), script.root_salt_parameters.reference_temperature_k, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 8e-6), script.root_salt_parameters.aqueous_diffusivity_m2_per_h_at_reference[7], 1.0e-18);
    try std.testing.expectApproxEqAbs(@as(f64, 5.9), script.root_salt_parameters.aqueous_diffusivity_temperature_exponent, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 8e-5), script.root_salt_parameters.root_concentration_inhibition_mol_per_m3[7], 1.0e-18);
    try std.testing.expect(!script.uses_compatibility_root_mycorrhizal_exchange_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.07), script.root_mycorrhizal_exchange_parameters.minimum_partner_water_volume_ratio, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.45), script.root_mycorrhizal_exchange_parameters.exchange_fraction_per_h, 1.0e-15);
}

test "UPTAKE root exudation coefficients are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\RoOt_ExUdAtIoN|999|0.11|0.12|0.0011
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_root_exudation_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 999), script.root_exudation_parameters.maximum_root_carbon_concentration_g_c_per_m3, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.11), script.root_exudation_parameters.root_mobile_nitrogen_exchange_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0011), script.root_exudation_parameters.exchange_rate_per_h, 1.0e-15);
}

test "GROSUB root porosity coefficients are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\RoOt_PoRoSiTy|0.74|0.11|0.012
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_root_porosity_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.74), script.root_porosity_parameters.maximum_porosity_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.11), script.root_porosity_parameters.oxygen_stress_induction_fraction_per_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.012), script.root_porosity_parameters.relaxation_fraction_per_h, 1.0e-15);
}

test "GROSUB root metabolism coefficients are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\RoOt_MeTaBoLiSm|0.016|0.026|0.11|0.012|0.009|1.75|0.16|0.34|0.66|0.68|2.6e-5|0.17|8.31|711|62600|197600|25.3|999|0.051|0.101|0.26|0.0011|1.1|4.1|0.26|1.1|2.1|4.1|337|2.6|26|0.87|0.76|0.51|481
        \\OrGaN_PaRtItIoN|0.74|0.26|0.021|0.051|0.007|0.017|0.76|1.51|2.01|2.02|1.76|1.52|0.26|0.51|0.68|0.69|0.59|0.52|0.11|0.12
        \\ShOoT_MeTaBoLiSm|0.016|0.026|0.011|0.34|1.71|0.026|0.11|0.012
        \\ShOoT_NoDe_GrOwTh|2.6|26|-0.334|-0.51|-0.668|0.0021|2.1|2.2|0.0011|0.0051|0.051|0.76|0.0051|0.0011|73|361|1441|721|722|0.0052|0.0053|0.0000052|0.0000053|0.000052|0.00052
        \\BrAnCh_MoBiLe_ExChAnGe|0.012|0.013|0.014
        \\SyMbIoTiC_FiXaTiOn|0.00011|0.126|0.011|0.26|0.00011|0.000011|11|1001|0.011|0.168|0.332|0.334|0.335|0.21|0.11|0.051|0.22|0.12|0.052|0.51|0.26|0.126|0.051|0.026|0.013
        \\PlAnT_FiRe_CoMbUsTiOn|474.15|1.9|8.32|60001|12.029|901|902|903|904|905|120001|20.621|1.1|2.668|0.01|0.51|2101|11|0.0376|0.0126|0.0744
        \\SoIl_FiRe_CoMbUsTiOn|475.15|1.8|8.33|60002|12.030|120002|20.622|1001|5001|1002|1003|1004|1005|1.2|2.669|0.02|0.52|2.9|0.006|0.0377|0.0127|0.0745
        \\ShOoT_RoOt_ExChAnGe|0.051|0.0051|0.26|0.91
        \\StOrAgE_ReMoBiLiZaTiOn|46|139|0.016|0.006|0.24|0.26|0.76|0.74|0.11|0.12|0.013|0.14|0.15|0.16|0.051|0.21|0.0051|0.021|0.11
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_root_metabolism_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.016), script.root_metabolism_parameters.maximum_substrate_respiration_fraction_per_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.026), script.root_metabolism_parameters.substrate_respiration_half_saturation_g_c_per_g_c, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.75), script.root_metabolism_parameters.nitrogen_assimilation_respiration_g_c_per_g_n, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.6), script.root_metabolism_parameters.root_protein_carbon_per_nitrogen_g_c_per_g_n, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.87), script.root_metabolism_parameters.nutrient_uptake_respiration_g_c_per_g_element, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.16), script.root_metabolism_parameters.minimum_carbon_recycling_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.6e-5), script.root_metabolism_parameters.storage_exchange_fraction_per_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 8.31), script.root_metabolism_parameters.maintenance_gas_constant_j_per_mol_k, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.1), script.root_metabolism_parameters.acidity_half_effect_hydrogen_activity_mol_per_m3, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 4.1), script.root_metabolism_parameters.deeper_primary_root_sink_multiplier, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 337), script.root_metabolism_parameters.annual_termination_hours_without_grain_fill, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.76), script.root_metabolism_parameters.evergreen_leafoff_remobilization_start_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.51), script.root_metabolism_parameters.deciduous_leafoff_remobilization_start_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 481), script.root_metabolism_parameters.full_senescence_duration_h, 1.0e-15);
    try std.testing.expect(!script.uses_compatibility_organ_partition_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.74), script.organ_partition_parameters.initial_leaf_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.02), script.organ_partition_parameters.leaf_reduction_by_turnover[3], 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.12), script.organ_partition_parameters.low_reserve_redirect_fraction, 1.0e-15);
    try std.testing.expect(!script.uses_compatibility_shoot_metabolism_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.016), script.shoot_metabolism_parameters.maximum_mobile_carbon_oxidation_per_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.34), script.shoot_metabolism_parameters.minimum_leaf_nutrient_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.71), script.shoot_metabolism_parameters.nitrogen_assimilation_respiration_g_c_per_g_n, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.11), script.shoot_metabolism_parameters.mobile_nitrogen_inhibition_g_n_per_g_c, 1.0e-15);
    try std.testing.expect(!script.uses_source_shoot_node_growth_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, -0.334), script.shoot_node_growth_parameters.leaf_mass_exponent, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.2), script.shoot_node_growth_parameters.minimum_internode_carbon_g_c_per_m2_cell, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.051), script.shoot_node_growth_parameters.branch_reserve_nutrient_exchange_fraction_per_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 73), script.shoot_node_growth_parameters.physiological_maturity_no_fill_h, 1.0e-15);
    try std.testing.expectEqual(@as([4]f64, .{ 361, 1441, 721, 722 }), script.shoot_node_growth_parameters.annual_leafoff_delay_h_by_phenology);
    try std.testing.expectApproxEqAbs(@as(f64, 0.00052), script.shoot_node_growth_parameters.leaf_storage_exchange_fraction_per_h_by_turnover[5], 1.0e-15);
    try std.testing.expect(!script.uses_source_branch_mobile_exchange_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.012), script.branch_mobile_exchange_parameters.carbon_exchange_fraction_per_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.013), script.branch_mobile_exchange_parameters.nutrient_exchange_fraction_per_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.014), script.branch_mobile_exchange_parameters.remobilization_redistribution_fraction_per_h, 1.0e-15);
    try std.testing.expect(!script.uses_source_symbiotic_fixation_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.00011), script.symbiotic_fixation_parameters.initial_bacterial_carbon_g_c_per_m2, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.052), script.symbiotic_fixation_parameters.host_exchange_fraction_per_h_by_fixation_type[5], 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.013), script.symbiotic_fixation_parameters.decomposition_control_ratio_by_fixation_type[5], 1.0e-15);
    try std.testing.expect(!script.uses_source_plant_fire_combustion_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 474.15), script.plant_fire_combustion_parameters.minimum_combustion_temperature_k, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 901), script.plant_fire_combustion_parameters.mobile_and_leaf_specific_combustion_g_c_per_m2_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 903), script.plant_fire_combustion_parameters.root_structural_specific_combustion_g_c_per_m2_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 905), script.plant_fire_combustion_parameters.standing_dead_specific_combustion_g_c_per_m2_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 120001), script.plant_fire_combustion_parameters.charcoal_activation_energy_j_per_mol, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.1), script.plant_fire_combustion_parameters.charcoal_specific_combustion_g_c_per_m2_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.668), script.plant_fire_combustion_parameters.oxygen_g_per_g_combusted_carbon, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.51), script.plant_fire_combustion_parameters.maximum_anaerobic_charcoal_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0744), script.plant_fire_combustion_parameters.methane_combustion_energy_mj_per_g_c, 1.0e-15);
    try std.testing.expect(!script.uses_source_soil_fire_combustion_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 475.15), script.soil_fire_combustion_parameters.minimum_combustion_temperature_k, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 5001), script.soil_fire_combustion_parameters.specific_combustion_by_substrate_g_c_per_m2_h[1], 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), script.soil_fire_combustion_parameters.charcoal_specific_combustion_g_c_per_m2_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.9), script.soil_fire_combustion_parameters.oxygen_half_saturation_g_o_per_m3, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.006), script.soil_fire_combustion_parameters.methane_half_saturation_g_c_per_m3, 1.0e-15);
    try std.testing.expect(!script.uses_compatibility_shoot_root_exchange_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.051), script.shoot_root_exchange_parameters.minimum_partner_structural_ratio, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0051), script.shoot_root_exchange_parameters.minimum_annual_carbon_exchange_fraction_per_h, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.91), script.shoot_root_exchange_parameters.salt_exchange_fraction_per_h, 1.0e-15);
    try std.testing.expect(!script.uses_compatibility_storage_remobilization_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 46), script.storage_remobilization_parameters.remobilization_duration_h[0], 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.16), script.storage_remobilization_parameters.perennial_nutrient_equilibration_fraction_per_h[5], 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.051), script.storage_remobilization_parameters.minimum_mobile_nitrogen_per_carbon_g_n_per_g_c, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.21), script.storage_remobilization_parameters.maximum_mobile_nitrogen_per_carbon_g_n_per_g_c, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0051), script.storage_remobilization_parameters.minimum_mobile_phosphorus_per_carbon_g_p_per_g_c, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.021), script.storage_remobilization_parameters.maximum_mobile_phosphorus_per_carbon_g_p_per_g_c, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.11), script.storage_remobilization_parameters.depleted_storage_threshold_g_c_per_g_root_c, 1.0e-15);
}

test "HFUNC phenology coefficients are case-insensitive runtime inputs" {
    const source =
        \\1,1,9
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\CaNoPy_StReSs|25|61|0.03|8.31|711|25.3|62600|197600|222600
        \\PhEnOlOgY_CoNtRoLs|8.3143|710|24.269|60000|197500|218500|0.1|0.25|2|0.667|-0.1|-150|-1.5|3600|1e-6|1e-6
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_phenology_parameters);
    try std.testing.expect(!script.uses_compatibility_canopy_stress_parameters);
    try std.testing.expectEqual(@as(usize, 9), script.plant_species_count);
    try std.testing.expectEqual(@as(f64, 25), script.canopy_stress_parameters.maximum_chilling_h);
    try std.testing.expectEqual(@as(f64, 61), script.canopy_stress_parameters.heat_accumulation_threshold_c);
    try std.testing.expectEqual(@as(f64, 0.03), script.canopy_stress_parameters.heat_recovery_per_h);
    try std.testing.expectEqual(@as(f64, 62600), script.canopy_stress_parameters.growth_temperature.activation_energy_j_per_mol);
    try std.testing.expectEqual(@as(f64, 222600), script.canopy_stress_parameters.growth_temperature.high_temperature_inactivation_j_per_mol);
    try std.testing.expectApproxEqAbs(@as(f64, 24.269), script.phenology_parameters.arrhenius_log_prefactor, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), script.phenology_parameters.oxygen_stress_exponent, 1.0e-15);
}

test "HFUNC plant pool thresholds are case-insensitive runtime inputs" {
    const source =
        \\1,1,8
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\PlAnT_PoOl_CoNtRoLs|YeS|2e-15|4e-6|2e-6|3e-6|0.02|0.002
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_plant_pool_parameters);
    try std.testing.expect(script.dynamic_plant_salts);
    try std.testing.expectApproxEqAbs(@as(f64, 4e-6), script.plant_pool_parameters.grain_fill_detection_g_c_per_plant, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.02), script.plant_pool_parameters.nitrogen_inhibition_g_n_per_g_c, 1.0e-15);
}

test "GROSUB seed set half saturations are case-insensitive runtime inputs" {
    const source =
        \\1,1,8
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\SeEd_SeT_CoNtRoLs|0.03|0.006|0.0012
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_seed_set_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.03), script.seed_set_parameters.carbon_half_saturation_g_per_g, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0012), script.seed_set_parameters.phosphorus_half_saturation_g_per_g, 1.0e-15);
}

test "plant nutrient initial zone and phosphate fractions are runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\PlAnT_NuTrIeNtS,0.2,0.3,0.4,0.65,0.5,0.6,0.7
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_plant_nutrient_initialization);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), script.plant_nutrient_initialization.initial_ammonium_band_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.65), script.plant_nutrient_initialization.initial_h2po4_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), script.plant_nutrient_initialization.initial_nitrate_band_row_spacing_m, 1.0e-15);
}

test "organic initialization parameter filename is case insensitive and extension agnostic" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\OrGaNiC_InItIaLiZaTiOn_FiLe|organic-parameters
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqualStrings("organic-parameters", script.organic_initialization_file.?);
    try std.testing.expectEqualStrings("grid_inputs", script.grid_input_file);
}

test "surface gas parameter filename is case insensitive and extension agnostic" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\SuRfAcE_GaS_PaRaMeTeR_FiLe|surface-gas.csv
        \\SoIl_NiTrOgEn_PaRaMeTeR_FiLe|soil-nitrogen
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectEqualStrings("surface-gas.csv", script.surface_gas_parameter_file.?);
    try std.testing.expectEqualStrings("soil-nitrogen", script.soil_nitrogen_parameter_file.?);
    try std.testing.expectEqualStrings("grid_inputs", script.grid_input_file);
}

test "STARTE initial phosphate constants are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\ChEmIsTrY_InItIaLiZaTiOn,0.01,1e-8,1.9e-21,6.3e-26,7.5,6.2e-5,4.8e-10
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), script.chemistry_initialization.?.saturated_paste_phosphate_multiplier, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 6.2e-5), script.chemistry_initialization.?.phosphate_dissociation.h2po4_to_hpo4_mol_per_m3, 1.0e-15);
}

test "STARTE molar masses and extract conversions are runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\chemistry_initialization,0.01,1e-8,1.9e-21,6.3e-26,7.5,6.2e-5,4.8e-10
        \\ChEmIsTrY_UnItS,14,31,27,56,40,24.3,23,39.1,32,35.5,1,1,0.01,1
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 24.3), script.chemistry_primary_initialization.?.molar_mass_g_per_mol.magnesium, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), script.chemistry_primary_initialization.?.soil_ammonium_extract_multiplier, 1e-15);
}

test "soil solver scientific coefficients and class count are runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\Soil_Solver,-0.0004,-1e10,200000,0.25,-0.2,0.3,5e-7,0.02,0.4,3e-7,0.08,0.2,0.3,0.6,0.7,0.05,0.16,0.23,0.74,0.72,0.55,0.45,73,1.4,0.12,1.5,0.03,0.11,70,2e-15,0.067,0.125,890,10000,0.8,0.0125,120,10000,0.52,10,1.82,200,80,20,5,1,1
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_soil_solver_parameters);
    try std.testing.expectEqual(@as(usize, 73), script.soil_solver_parameters.hydraulic_conductivity_class_count);
    try std.testing.expectApproxEqAbs(@as(f64, -0.0004), script.soil_solver_parameters.retention.saturation_water_potential_mpa, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.4), script.soil_solver_parameters.pore_interaction_exponent, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 70), script.soil_solver_parameters.organic_saturated_conductivity_scale_m2_per_h_mpa, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 890), script.soil_solver_parameters.profile_derivation.organic_nitrogen_scale_g_per_megagram, 1e-15);
}

test "soil process coefficients are runtime inputs and tag is case insensitive" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\SoIl_PrOcEsS,0.0097,0.08,300,1.8,0.65,0.04,0.0006,1.1e-6,0.52,0.025
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_soil_process_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.08), script.soil_process_parameters.reference_water_vapor_diffusivity_m2_per_h, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.04), script.soil_process_parameters.osmotic_reflection_coefficient, 1e-15);
}

test "macropore van Genuchten shape is a case-insensitive runtime record" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\MaCrOpOrE_VaN_GeNuChTeN|0.01|16.5|2.75|0.45
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), script.soil_process_parameters.macropore_residual_saturation, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 16.5), script.soil_process_parameters.macropore_van_genuchten_alpha_per_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.75), script.soil_process_parameters.macropore_van_genuchten_n, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.45), script.soil_process_parameters.macropore_pore_connectivity, 1e-15);
}

test "dual domain exchange coefficients are runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\DuAl_DoMaIn_ExChAnGe,4.5,0.35
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 4.5), script.soil_process_parameters.dual_domain_geometry_factor, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.35), script.soil_process_parameters.dual_domain_scaling_coefficient, 1e-15);
}

test "frozen hydraulic impedance exponent is a case-insensitive runtime input" {
    const source =
        \\1,1,2
        \\FrOzEn_HyDrAuLiC_ImPeDaNcE|6.25
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectApproxEqAbs(
        @as(f64, 6.25),
        script.soil_process_parameters.frozen_hydraulic_impedance_exponent,
        1e-15,
    );
}

test "surface residue freeze thaw retention is a case-insensitive runtime input" {
    const source =
        \\1,1,2
        \\SuRfAcE_ReSiDuE_FrEeZe_ThAw|0.03|250|2.2
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 0.03), script.soil_process_parameters.surface_residue_residual_water_content_m3_per_m3, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 250), script.soil_process_parameters.surface_residue_van_genuchten_alpha_per_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.2), script.soil_process_parameters.surface_residue_van_genuchten_n, 1e-15);
}

test "REDIST soil geometry coefficients are case-insensitive runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\SoIl_GeOmEtRy|1.9e-6|120000|0.081|2e-9
        \\SuRfAcE_PoNd_EnErGy|2.6e-6|9e-4
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_soil_geometry_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 1.9e-6), script.soil_geometry_parameters.organic_carbon_specific_volume_m3_per_g, 1e-18);
    try std.testing.expectEqual(@as(f64, 120_000), script.soil_geometry_parameters.organic_horizon_threshold_g_c_per_Mg);
    try std.testing.expectApproxEqAbs(@as(f64, 0.081), script.soil_geometry_parameters.ice_to_water_specific_volume_difference, 1e-15);
    try std.testing.expectEqual(@as(f64, 2e-9), script.soil_geometry_parameters.minimum_layer_thickness_m);
    try std.testing.expect(!script.uses_compatibility_surface_pond_energy);
    try std.testing.expectEqual(@as(f64, 2.6e-6), script.surface_pond_dry_organic_heat_capacity_mj_per_g_c_k);
    try std.testing.expectEqual(@as(f64, 9e-4), script.surface_pond_activation_heat_capacity_mj_per_m2_k);
}

test "seven gas diffusivities and transport constants are runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\SoIl_GaS_TrAnSpOrT|299|1.8|0.041|0.072|0.063|0.054|0.055|0.066|0.057|0.64|1e-10|998000|18.1
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_soil_gas_transport_parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 299), script.soil_gas_transport_parameters.reference_temperature_k, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.041), script.soil_gas_transport_parameters.free_air_diffusivity_m2_per_h[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.057), script.soil_gas_transport_parameters.free_air_diffusivity_m2_per_h[6], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 18.1), script.soil_gas_transport_parameters.water_molar_mass_g_per_mol, 1e-15);
}

test "surface runoff controls accept mixed delimiters and case" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\SuRfAcE_RuNoFf|0.006|0.012|0.0009|3601|2e-12
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_surface_runoff);
    try std.testing.expectApproxEqAbs(@as(f64, 0.006), script.surface_runoff_parameters.ground_surface_retention_m3_per_m2, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0009), script.surface_runoff_parameters.maximum_hydraulic_volume_m3, 1e-15);
}

test "rainfall impact coefficients are runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\RaInFaLl_ImPaCt,9,8,16,6,3,1.9,0.002,0.0007
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_rainfall_impact);
    try std.testing.expectApproxEqAbs(@as(f64, 16), script.rainfall_impact_parameters.throughfall_energy_height_coefficient_j_per_mm_sqrt_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.002), script.rainfall_impact_parameters.conductivity_damage_per_j_per_megagram_per_megagram, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0007), script.rainfall_impact_parameters.conductivity_recovery_fraction_per_h, 1e-15);
}

test "soil phase and heat coefficients are runtime inputs" {
    const source =
        \\1,1,2
        \\runtime,1,1,1e-8,1e-11,40,0.5
        \\soil_phase_heat,0.002,0.62,5300,0.0036,18.1,8.31,2400,90000,330,0.92,0.006,273.2,0.3,0.34,1e10,8e10,1.2,1.3,9000,4.2,1.9
        \\GeOtHeRmAl,YeS,12,1.5,0.009,0.00021
        \\SuBsUrFaCe_StAtE,0.002,0.000003
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_soil_phase_heat_parameters);
    try std.testing.expectEqual(@as(f64, 2400), script.soil_phase_heat_parameters.vapor_equilibrium.latent_heat_of_vaporization_mj_per_m3);
    try std.testing.expectEqual(@as(f64, 9000), script.soil_phase_heat_parameters.heat_turbulence.maximum_rayleigh_number);
    try std.testing.expectEqual(@as(f64, 4.2), script.soil_phase_heat_parameters.liquid_water_heat_capacity_mj_per_m3_k);
    try std.testing.expect(!script.uses_compatibility_geothermal_controls);
    try std.testing.expect(script.geothermal_controls.enabled);
    try std.testing.expectEqual(@as(f64, 12), script.geothermal_controls.minimum_source_depth_m);
    try std.testing.expectEqual(@as(f64, 0.00021), script.geothermal_controls.geothermal_flux_mj_per_m2_h);
    try std.testing.expect(!script.uses_compatibility_subsurface_state_controls);
    try std.testing.expectEqual(@as(f64, 0.002), script.water_table_air_fraction_threshold);
    try std.testing.expectEqual(@as(f64, 0.000003), script.active_layer_ice_fraction_threshold);
}

test "surface energy controls are tagged and case insensitive" {
    const source =
        \\1,1,8
        \\Runtime,2,64,1e-8,1e-11,40,0.5
        \\Surface_Energy,0.95,0.98,0.96,0.08,0.5,0.12,0.9,180,380
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 0.95), script.soil_longwave_emissivity, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.98), script.snow_longwave_emissivity, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.96), script.canopy_longwave_emissivity, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.08), script.snow_full_cover_depth_m, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), script.surface_sensible_heat_conductance_mj_per_m2_h_k, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.12), script.surface_latent_heat_conductance_mj_per_m2_h_kpa, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.9), script.surface_vapor_activity_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 180), script.minimum_surface_temperature_k, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 380), script.maximum_surface_temperature_k, 1.0e-15);
    try std.testing.expect(!script.uses_compatibility_surface_energy);
}

test "canopy surface exchange coefficients are runtime and case insensitive" {
    const source =
        \\1,1,8
        \\CaNoPy_SuRfAcE_ExChAnGe,0.0015,0.028,0.0014,0.014,17.9,8.314,4.18
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_canopy_surface_exchange);
    try std.testing.expectEqual(@as(f64, 0.0015), script.canopy_sensible_surface_resistance_h_per_m);
    try std.testing.expectEqual(@as(f64, 0.028), script.canopy_latent_surface_resistance_h_per_m);
    try std.testing.expectEqual(@as(f64, 0.0014), script.canopy_surface_exchange_parameters.minimum_boundary_resistance_h_per_m);
    try std.testing.expectEqual(@as(f64, 0.014), script.canopy_surface_exchange_parameters.maximum_boundary_resistance_h_per_m);
    try std.testing.expectEqual(@as(f64, 17.9), script.canopy_surface_exchange_parameters.water_potential_vapor_coefficient_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 8.314), script.canopy_surface_exchange_parameters.universal_gas_constant_j_per_mol_k);
    try std.testing.expectEqual(@as(f64, 4.18), script.canopy_surface_exchange_parameters.liquid_water_heat_capacity_mj_per_m3_k);
}

test "canopy ammonia exchange coefficients are runtime and case insensitive" {
    const source =
        \\1,1,8
        \\CaNoPy_AmMoNiA_ExChAnGe,0.17,0.11,0.06,2.1,0.00012,0.09,0.52,0.018
        \\grid_inputs grid_inputs
        \\1,1
        \\1,1
        \\weather_grid weather options NO plants o1 o2 o3 o4 o5 o6 o7 o8 o9 o10
        \\0,0
    ;
    var script = try parse(std.testing.allocator, source);
    defer script.deinit();
    try std.testing.expect(!script.uses_compatibility_canopy_ammonia_exchange);
    try std.testing.expectEqual(@as(f64, 0.17), script.canopy_ammonia_exchange_parameters.minimum_canopy_dry_matter_fraction);
    try std.testing.expectEqual(@as(f64, 0.11), script.canopy_ammonia_exchange_parameters.water_potential_dry_matter_increment);
    try std.testing.expectEqual(@as(f64, 0.06), script.canopy_ammonia_exchange_parameters.water_potential_denominator_per_mpa);
    try std.testing.expectEqual(@as(f64, 0.09), script.canopy_ammonia_exchange_parameters.maximum_mobile_nitrogen_transfer_fraction_per_step);
    try std.testing.expectEqual(@as(f64, 0.018), script.canopy_ammonia_exchange_parameters.solubility_temperature_coefficient_per_c);
}
