const std = @import("std");
const CellRange = @import("compute.zig").CellRange;
const GridState = @import("grid.zig").GridState;
const PlantState = @import("grid.zig").PlantState;
const numerics = @import("numerics.zig");
const PlantRootState = @import("plant_root_system.zig").State;
const CanopyState = @import("canopy_photosynthesis.zig").State;
const SoilPropertiesState = @import("soil_solver_properties.zig").State;
const soil_water_solver = @import("soil_water_solver.zig");
const biological_domain_count = @import("plant_root_system.zig").biological_domain_count;

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    species_count: usize,
    soil_layer_count: usize,
    root_domain_count: usize,
    transpiration_loss_m: []f64,
    total_root_water_uptake_m: []f64,
    water_balance_residual_m: []f64,
    root_water_uptake_m: []f64,
    iteration_count: []u16,
    newton_raphson_step_count: []u16,
    picard_step_count: []u16,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, species_count: usize, soil_layer_count: usize) !State {
        if (cell_count == 0 or species_count == 0 or soil_layer_count == 0) return error.InvalidPlantWaterDimensions;
        const plant_count = try std.math.mul(usize, cell_count, species_count);
        const root_count = try std.math.mul(usize, try std.math.mul(usize, plant_count, biological_domain_count), soil_layer_count);
        const transpiration = try allocator.alloc(f64, plant_count);
        errdefer allocator.free(transpiration);
        const total_uptake = try allocator.alloc(f64, plant_count);
        errdefer allocator.free(total_uptake);
        const balance_residual = try allocator.alloc(f64, plant_count);
        errdefer allocator.free(balance_residual);
        const root_uptake = try allocator.alloc(f64, root_count);
        errdefer allocator.free(root_uptake);
        const iterations = try allocator.alloc(u16, plant_count);
        errdefer allocator.free(iterations);
        const newton_steps = try allocator.alloc(u16, plant_count);
        errdefer allocator.free(newton_steps);
        const picard_steps = try allocator.alloc(u16, plant_count);
        @memset(transpiration, 0);
        @memset(total_uptake, 0);
        @memset(balance_residual, 0);
        @memset(root_uptake, 0);
        @memset(iterations, 0);
        @memset(newton_steps, 0);
        @memset(picard_steps, 0);
        return .{ .allocator = allocator, .cell_count = cell_count, .species_count = species_count, .soil_layer_count = soil_layer_count, .root_domain_count = biological_domain_count, .transpiration_loss_m = transpiration, .total_root_water_uptake_m = total_uptake, .water_balance_residual_m = balance_residual, .root_water_uptake_m = root_uptake, .iteration_count = iterations, .newton_raphson_step_count = newton_steps, .picard_step_count = picard_steps };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.picard_step_count);
        self.allocator.free(self.newton_raphson_step_count);
        self.allocator.free(self.iteration_count);
        self.allocator.free(self.root_water_uptake_m);
        self.allocator.free(self.water_balance_residual_m);
        self.allocator.free(self.total_root_water_uptake_m);
        self.allocator.free(self.transpiration_loss_m);
        self.* = undefined;
    }

    pub fn validateFinite(self: State) !void {
        inline for (.{ self.transpiration_loss_m, self.total_root_water_uptake_m, self.water_balance_residual_m, self.root_water_uptake_m }) |values| for (values, 0..) |value, index| if (!std.math.isFinite(value)) {
            std.log.err("non-finite plant water balance: index={d} value={e}", .{ index, value });
            return error.NonFinitePlantWaterBalance;
        };
    }
};

pub fn resetDailyMinimumCanopyWaterPotential(canopy: *CanopyState) void {
    @memset(canopy.plant_minimum_daily_canopy_water_potential_mpa, 0);
}

/// UPTAKE PSILZ: retain the most negative canopy total water potential reached
/// since the DAY/HOUR1 daily reset.
pub fn updateDailyMinimumCanopyWaterPotential(canopy: *CanopyState, plants: *const PlantState) !void {
    if (canopy.plant_minimum_daily_canopy_water_potential_mpa.len != plants.canopy_water_potential_mpa.len) return error.PlantWaterDimensionMismatch;
    for (plants.canopy_water_potential_mpa, canopy.plant_minimum_daily_canopy_water_potential_mpa) |potential_mpa, *minimum_mpa| {
        if (!std.math.isFinite(potential_mpa)) return error.NonFiniteCanopyWaterPotential;
        minimum_mpa.* = @min(minimum_mpa.*, potential_mpa);
    }
}

/// DAY automatic irrigation reads PSILZ for source PFT slot one only.
pub fn publishFirstPlantDailyMinimumByCell(canopy: *const CanopyState, output_mpa_by_cell: []f64) !void {
    if (output_mpa_by_cell.len != canopy.cell_count or canopy.species_count == 0) return error.PlantWaterDimensionMismatch;
    for (output_mpa_by_cell, 0..) |*output, cell| {
        output.* = canopy.plant_minimum_daily_canopy_water_potential_mpa[cell * canopy.species_count];
        if (!std.math.isFinite(output.*)) return error.NonFiniteCanopyWaterPotential;
    }
}

pub const Workspace = struct {
    allocator: std.mem.Allocator,
    plant_count: usize,
    root_count: usize,
    active: []bool,
    dynamically_woody: []bool,
    vascular_growth_habit: []bool,
    cell_area_m2: []f64,
    root_conductance_m_per_h_mpa: []f64,
    root_biome_fraction: []f64,
    soil_path_length_m: []f64,
    root_cylinder_radius_m: []f64,
    root_surface_area_per_radius_m: []f64,
    maximum_uptake_m: []f64,
    maximum_release_m: []f64,
    soil_resistance_mpa_h_per_m: []f64,
    root_resistance_mpa_h_per_m: []f64,
    canopy_water_capacitance_m_per_m2_mpa: []f64,
    transpiration_loss_m: []f64,
    leaf_osmotic_potential_at_zero_total_mpa: []f64,
    plant_population_count: []f64,
    root_porosity_fraction: []f64,
    root_radial_resistivity_mpa_h_per_m3: []f64,
    root_axial_resistivity_mpa_h_per_m2: []f64,
    primary_root_radius_m: []f64,
    secondary_root_radius_m: []f64,
    seeding_depth_m: []f64,
    woody_root_fraction: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, species_count: usize, soil_layer_count: usize) !Workspace {
        if (cell_count == 0 or species_count == 0 or soil_layer_count == 0) return error.InvalidPlantWaterDimensions;
        const plant_count = try std.math.mul(usize, cell_count, species_count);
        const root_count = try std.math.mul(usize, try std.math.mul(usize, plant_count, biological_domain_count), soil_layer_count);
        var result: Workspace = undefined;
        result.allocator = allocator;
        result.plant_count = plant_count;
        result.root_count = root_count;
        result.active = try allocator.alloc(bool, plant_count);
        @memset(result.active, false);
        errdefer allocator.free(result.active);
        result.dynamically_woody = try allocator.alloc(bool, plant_count);
        @memset(result.dynamically_woody, false);
        errdefer allocator.free(result.dynamically_woody);
        result.vascular_growth_habit = try allocator.alloc(bool, plant_count);
        @memset(result.vascular_growth_habit, false);
        errdefer allocator.free(result.vascular_growth_habit);
        var allocated: usize = 0;
        errdefer freeWorkspaceFields(&result, allocated);
        inline for (@typeInfo(Workspace).@"struct".fields) |field| if (field.type == []f64) {
            const count = if (comptime std.mem.eql(u8, field.name, "cell_area_m2")) cell_count else if (comptime std.mem.eql(u8, field.name, "root_conductance_m_per_h_mpa") or std.mem.eql(u8, field.name, "root_biome_fraction") or std.mem.eql(u8, field.name, "root_resistance_mpa_h_per_m") or std.mem.eql(u8, field.name, "soil_path_length_m") or std.mem.eql(u8, field.name, "root_cylinder_radius_m") or std.mem.eql(u8, field.name, "root_surface_area_per_radius_m") or std.mem.startsWith(u8, field.name, "maximum_") or std.mem.startsWith(u8, field.name, "soil_resistance_")) root_count else plant_count;
            @field(result, field.name) = try allocator.alloc(f64, count);
            @memset(@field(result, field.name), 0);
            allocated += 1;
        };
        return result;
    }

    pub fn deinit(self: *Workspace) void {
        inline for (@typeInfo(Workspace).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.allocator.free(self.vascular_growth_habit);
        self.allocator.free(self.dynamically_woody);
        self.allocator.free(self.active);
        self.* = undefined;
    }

    pub fn refreshActive(self: *Workspace, soil_layer_count: usize) !void {
        if (soil_layer_count == 0 or self.root_count != try std.math.mul(usize, try std.math.mul(usize, self.plant_count, biological_domain_count), soil_layer_count)) return error.PlantWaterDimensionMismatch;
        const roots_per_plant = biological_domain_count * soil_layer_count;
        for (0..self.plant_count) |plant| {
            var conductance: f64 = 0;
            for (self.root_conductance_m_per_h_mpa[plant * roots_per_plant ..][0..roots_per_plant]) |value| {
                if (!std.math.isFinite(value) or value < 0) return error.InvalidRootHydraulicInput;
                conductance += value;
            }
            const capacitance = self.canopy_water_capacitance_m_per_m2_mpa[plant];
            const transpiration = self.transpiration_loss_m[plant];
            if (!std.math.isFinite(capacitance) or capacitance < 0 or !std.math.isFinite(transpiration) or transpiration < 0) return error.InvalidPlantWaterRuntimeInput;
            self.active[plant] = conductance > 0 and capacitance > 0;
        }
    }
};

fn freeWorkspaceFields(workspace: *Workspace, allocated_count: usize) void {
    var visited: usize = 0;
    inline for (@typeInfo(Workspace).@"struct".fields) |field| if (field.type == []f64) {
        if (visited < allocated_count) workspace.allocator.free(@field(workspace, field.name));
        visited += 1;
    };
}

pub fn refreshCanopyWorkspace(workspace: *Workspace, canopy: *const CanopyState, plants: *const PlantState, species_count: usize) !void {
    if (species_count == 0 or canopy.species_count != species_count or canopy.plant_branch_offsets.len != workspace.plant_count + 1 or plants.canopy_water_potential_mpa.len != workspace.plant_count or workspace.cell_area_m2.len * species_count != workspace.plant_count) return error.PlantWaterDimensionMismatch;
    for (0..workspace.plant_count) |plant| {
        const cell = plant / species_count;
        const area = workspace.cell_area_m2[cell];
        if (!std.math.isFinite(area) or area <= 0) return error.InvalidCellArea;
        var leaf_and_petiole_carbon_g: f64 = 0;
        var stalk_carbon_g: f64 = 0;
        var sapwood_carbon_g: f64 = 0;
        var stalk_surface_area_m2: f64 = 0;
        const branches = try canopy.branchRange(plant);
        for (branches.first..branches.end) |branch| {
            leaf_and_petiole_carbon_g += canopy.branch_leaf_carbon_g[branch] + canopy.branch_sheath_carbon_g[branch];
            stalk_carbon_g += canopy.branch_stalk_carbon_g[branch];
            sapwood_carbon_g += canopy.branch_sapwood_carbon_g[branch];
            const nodes = try canopy.nodeRange(branch);
            for (nodes.first..nodes.end) |node| {
                const samples = try canopy.sampleRange(node);
                for (samples.first..samples.end) |sample| stalk_surface_area_m2 += canopy.sample_stalk_area_m2[sample];
            }
        }
        const water_m3 = @max(0.0, plants.canopy_water_storage_m_per_m2[plant] * area);
        const hydraulics = try canopyHydraulics(leaf_and_petiole_carbon_g, stalk_carbon_g, stalk_surface_area_m2, water_m3, 0, plants.canopy_water_potential_mpa[plant]);
        workspace.canopy_water_capacitance_m_per_m2_mpa[plant] = try canopyWaterCapacitanceM3PerMpa(hydraulics.active_carbon_g, plants.canopy_water_potential_mpa[plant]) / area;
        workspace.transpiration_loss_m[plant] = @max(0.0, -canopy.plant_transpiration_m3_per_h[plant] / area);
        workspace.woody_root_fraction[plant] = if (workspace.dynamically_woody[plant] and stalk_carbon_g > 0) std.math.pow(f64, std.math.clamp(sapwood_carbon_g / stalk_carbon_g, 0, 1), 0.167) else 1.0;
    }
}

pub fn refreshRootWorkspace(workspace: *Workspace, roots: *PlantRootState, canopy: *const CanopyState, plants: *const PlantState, grid: *const GridState, properties: *const SoilPropertiesState, species_count: usize, biological_domain_count_by_plant: []const u8, root_volume_numerator_m3_per_g_c: f64, root_dry_matter_fraction: f64, root_geometry_pi: f64, morphology_parameters: @import("plant_root_system.zig").MorphologyParameters) !void {
    try morphology_parameters.validate();
    inline for (.{ root_volume_numerator_m3_per_g_c, root_dry_matter_fraction, root_geometry_pi }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootMorphologyInput;
    if (root_volume_numerator_m3_per_g_c <= 0 or root_dry_matter_fraction <= 0 or root_geometry_pi <= 0) return error.InvalidRootMorphologyInput;
    if (species_count == 0 or grid.soil_layer_capacity != roots.soil_layer_count or properties.layer_count != grid.layer_count or workspace.plant_count != roots.plant_count or workspace.plant_count != grid.cell_count * species_count or workspace.root_count != roots.plant_count * biological_domain_count * roots.soil_layer_count or biological_domain_count_by_plant.len != roots.plant_count) return error.PlantWaterDimensionMismatch;
    @memset(workspace.root_conductance_m_per_h_mpa, 0);
    @memset(workspace.soil_path_length_m, 0);
    @memset(workspace.root_cylinder_radius_m, 0);
    @memset(workspace.root_surface_area_per_radius_m, 0);
    @memset(workspace.maximum_uptake_m, 0);
    @memset(workspace.maximum_release_m, 0);
    @memset(workspace.soil_resistance_mpa_h_per_m, 0);
    @memset(workspace.root_resistance_mpa_h_per_m, 0);
    @memset(workspace.root_biome_fraction, 0);

    // FPQ: each PFT's fraction of total root C in the cell/domain/layer.
    for (0..grid.cell_count) |cell| for (0..biological_domain_count) |domain| for (0..grid.active_soil_layer_count[cell]) |layer| {
        var total_carbon_g: f64 = 0;
        for (0..species_count) |species| {
            const plant = cell * species_count + species;
            if (biological_domain_count_by_plant[plant] <= domain) continue;
            total_carbon_g += @max(0.0, roots.total_carbon_g[try roots.layerIndex(plant, domain, layer)]);
        }
        for (0..species_count) |species| {
            const plant = cell * species_count + species;
            if (biological_domain_count_by_plant[plant] < 1 or biological_domain_count_by_plant[plant] > biological_domain_count)
                return error.PlantWaterDimensionMismatch;
            if (biological_domain_count_by_plant[plant] <= domain) continue;
            const root = (plant * biological_domain_count + domain) * roots.soil_layer_count + layer;
            const carbon_g = @max(0.0, roots.total_carbon_g[try roots.layerIndex(plant, domain, layer)]);
            workspace.root_biome_fraction[root] = if (total_carbon_g > 0) carbon_g / total_carbon_g else 1.0;
        }
    };

    for (0..workspace.plant_count) |plant| {
        const population = workspace.plant_population_count[plant];
        if (!std.math.isFinite(population) or population < 0) return error.InvalidPlantPopulation;
        if (population == 0) continue;
        const cell = plant / species_count;
        var canopy_height_m: f64 = 0;
        const branches = try canopy.branchRange(plant);
        for (branches.first..branches.end) |branch| {
            const nodes = try canopy.nodeRange(branch);
            for (nodes.first..nodes.end) |node| canopy_height_m = @max(canopy_height_m, canopy.node_height_m[node]);
        }
        const stalk_factor = 3.75e3 * std.math.pow(f64, @max(0.5, 1.0 + plants.canopy_water_potential_mpa[plant] / 50.0), 4);
        const porosity = workspace.root_porosity_fraction[plant];
        if (!std.math.isFinite(porosity) or porosity < 0 or porosity >= 1) return error.InvalidRootPorosity;
        const volume_per_carbon = root_volume_numerator_m3_per_g_c / (root_dry_matter_fraction * (1.0 - porosity));
        for (0..biological_domain_count_by_plant[plant]) |domain| {
            var deepest_primary_root_m: f64 = 0;
            for (0..roots.root_axis_count) |axis| {
                var has_axis = false;
                for (0..grid.active_soil_layer_count[cell]) |layer| if (roots.axis_primary_count[try roots.layerAxisIndex(plant, domain, layer, axis)] > 0) {
                    has_axis = true;
                    break;
                };
                if (has_axis) deepest_primary_root_m = @max(deepest_primary_root_m, roots.axis_depth_m[try roots.axisIndex(plant, domain, axis)]);
            }
            var layer_top_depth_m: f64 = 0;
            for (0..grid.active_soil_layer_count[cell]) |layer| {
                const soil = cell * roots.soil_layer_count + layer;
                const root = (plant * biological_domain_count + domain) * roots.soil_layer_count + layer;
                const thickness = properties.layer_thickness_m[soil];
                const topology = try roots.refreshLayerMorphology(plant, domain, layer, population, thickness, workspace.woody_root_fraction[plant], porosity, volume_per_carbon, root_geometry_pi, workspace.vascular_growth_habit[plant], morphology_parameters);
                const rooted_fraction = try rootedLayerFraction(layer == 0, layer_top_depth_m, thickness, deepest_primary_root_m, workspace.seeding_depth_m[plant], 0);
                layer_top_depth_m += thickness;
                if (topology.root_length_density_m_per_m3 <= 0 or topology.primary_axis_count <= 0 or topology.secondary_axis_count <= 0 or rooted_fraction <= 0 or grid.matrix_liquid_water_m3[soil] <= 0) continue;
                const matrix_volume = properties.matrix_bulk_volume_m3[soil];
                if (matrix_volume <= 0) continue;
                const water_fraction = grid.matrix_liquid_water_m3[soil] / matrix_volume;
                const conductivity = try soil_water_solver.interpolatedConductivity(properties.matrix_hydraulic_conductivity_m2_per_h_mpa, properties.hydraulic_conductivity_class_count, properties.porosity_fraction[soil], soil, 2, water_fraction);
                if (conductivity <= 0) continue;
                const geometry = try rootUptakeGeometry(topology.root_length_density_m_per_m3, rooted_fraction, matrix_volume / properties.layer_volume_m3[soil], roots.aqueous_volume_m3[try roots.layerIndex(plant, domain, layer)], porosity, population, topology.root_length_m_per_plant, roots.secondary_radius_m[try roots.layerIndex(plant, domain, layer)], thickness);
                workspace.soil_path_length_m[root] = geometry.soil_path_length_m;
                workspace.root_cylinder_radius_m[root] = geometry.effective_radius_m;
                workspace.root_surface_area_per_radius_m[root] = geometry.surface_area_per_radius_m;
                const resistance = try hydraulicResistance(.{
                    .soil_path_length_m = geometry.soil_path_length_m,
                    .root_cylinder_radius_m = geometry.effective_radius_m,
                    .root_surface_area_per_radius_m = geometry.surface_area_per_radius_m,
                    .unsaturated_soil_conductivity_m_per_h_mpa = conductivity,
                    .root_radial_resistivity_mpa_h_per_m3 = workspace.root_radial_resistivity_mpa_h_per_m3[plant],
                    .soil_air_volume_m3 = properties.porosity_fraction[soil] * matrix_volume,
                    .soil_water_volume_m3 = grid.matrix_liquid_water_m3[soil],
                    .secondary_root_radius_m = roots.secondary_radius_m[try roots.layerIndex(plant, domain, layer)],
                    .root_length_m_per_plant = topology.root_length_m_per_plant,
                    .root_axial_resistivity_mpa_h_per_m2 = workspace.root_axial_resistivity_mpa_h_per_m2[plant],
                    .primary_root_depth_m = layer_top_depth_m,
                    .primary_root_radius_m = roots.primary_radius_m[try roots.layerIndex(plant, domain, layer)],
                    .primary_axis_count_per_plant = topology.primary_axis_count / population,
                    .canopy_height_m = 0.8 * canopy_height_m,
                    .stalk_conducting_element_factor = stalk_factor,
                    .secondary_root_length_m = topology.average_secondary_length_m,
                    .secondary_axis_count_per_plant = topology.secondary_axis_count / population,
                });
                workspace.root_conductance_m_per_h_mpa[root] = resistance.conductance_m_per_h_mpa;
                workspace.soil_resistance_mpa_h_per_m[root] = resistance.soil_mpa_h_per_m;
                workspace.root_resistance_mpa_h_per_m[root] = resistance.root_mpa_h_per_m;
                const available_depth_m = grid.matrix_liquid_water_m3[soil] * workspace.root_biome_fraction[root] / workspace.cell_area_m2[cell];
                workspace.maximum_uptake_m[root] = @max(0.0, available_depth_m);
                workspace.maximum_release_m[root] = @max(0.0, available_depth_m);
            }
        }
    }
}

pub const Settings = struct {
    minimum_canopy_water_potential_mpa: f64,
    maximum_canopy_water_potential_mpa: f64,
    solver_options: numerics.SolverOptions,
};

pub const HydraulicPath = struct {
    soil_path_length_m: f64,
    root_cylinder_radius_m: f64,
    root_surface_area_per_radius_m: f64,
    unsaturated_soil_conductivity_m_per_h_mpa: f64,
    root_radial_resistivity_mpa_h_per_m3: f64,
    soil_air_volume_m3: f64,
    soil_water_volume_m3: f64,
    secondary_root_radius_m: f64,
    root_length_m_per_plant: f64,
    root_axial_resistivity_mpa_h_per_m2: f64,
    primary_root_depth_m: f64,
    primary_root_radius_m: f64,
    primary_axis_count_per_plant: f64,
    canopy_height_m: f64,
    stalk_conducting_element_factor: f64,
    secondary_root_length_m: f64,
    secondary_axis_count_per_plant: f64,
};

pub const HydraulicResistance = struct {
    soil_mpa_h_per_m: f64,
    root_mpa_h_per_m: f64,
    total_mpa_h_per_m: f64,
    conductance_m_per_h_mpa: f64,
};

pub const CanopyHydraulics = struct {
    active_carbon_g: f64,
    dry_matter_fraction_g_c_per_g: f64,
    water_capacity_m3: f64,
    dry_heat_capacity_mj_per_k: f64,
    wet_heat_capacity_mj_per_k: f64,
};

pub fn canopyHydraulics(
    leaf_and_petiole_carbon_g: f64,
    stalk_carbon_g: f64,
    stalk_surface_area_m2: f64,
    canopy_water_m3: f64,
    intercepted_water_m3: f64,
    total_water_potential_mpa: f64,
) !CanopyHydraulics {
    inline for (.{ leaf_and_petiole_carbon_g, stalk_carbon_g, stalk_surface_area_m2, canopy_water_m3, intercepted_water_m3, total_water_potential_mpa }) |value| if (!std.math.isFinite(value)) return error.NonFiniteCanopyHydraulicInput;
    if (leaf_and_petiole_carbon_g < 0 or stalk_carbon_g < 0 or stalk_surface_area_m2 < 0 or canopy_water_m3 < 0 or intercepted_water_m3 < 0) return error.InvalidCanopyHydraulicInput;
    const minimum_dry_matter_fraction = 0.16;
    const stalk_sapwood_thickness_m = 0.0025;
    const stalk_volume_m3_per_g_c = 4.0e-6;
    const active_carbon_g = @max(0.0, leaf_and_petiole_carbon_g + @min(stalk_carbon_g, stalk_sapwood_thickness_m * stalk_surface_area_m2 / stalk_volume_m3_per_g_c));
    const absolute_potential = @abs(total_water_potential_mpa);
    const dry_matter_fraction = minimum_dry_matter_fraction + 0.10 * absolute_potential / (0.05 * absolute_potential + 2.0);
    const water_capacity_m3 = 1.0e-6 * active_carbon_g / dry_matter_fraction;
    const dry_heat_capacity = 2.496 * active_carbon_g * stalk_volume_m3_per_g_c;
    const wet_heat_capacity = dry_heat_capacity + 4.19 * (intercepted_water_m3 + canopy_water_m3);
    return .{ .active_carbon_g = active_carbon_g, .dry_matter_fraction_g_c_per_g = dry_matter_fraction, .water_capacity_m3 = water_capacity_m3, .dry_heat_capacity_mj_per_k = dry_heat_capacity, .wet_heat_capacity_mj_per_k = wet_heat_capacity };
}

/// Analytic d(VOLWPZ)/d(PSILT) on the physical PSILT <= 0 branch. This
/// replaces the source's secant RSSZ update inside its repeated sub-hour loop.
pub fn canopyWaterCapacitanceM3PerMpa(active_canopy_carbon_g: f64, total_water_potential_mpa: f64) !f64 {
    if (!std.math.isFinite(active_canopy_carbon_g) or active_canopy_carbon_g < 0 or !std.math.isFinite(total_water_potential_mpa) or total_water_potential_mpa > 0) return error.InvalidCanopyCapacitanceInput;
    const absolute_potential = -total_water_potential_mpa;
    const denominator = 0.05 * absolute_potential + 2.0;
    const dry_matter_fraction = 0.16 + 0.10 * absolute_potential / denominator;
    const dry_matter_derivative_per_mpa = 0.20 / (denominator * denominator);
    const capacitance = 1.0e-6 * active_canopy_carbon_g * dry_matter_derivative_per_mpa / (dry_matter_fraction * dry_matter_fraction);
    if (!std.math.isFinite(capacitance) or capacitance < 0) return error.InvalidCanopyCapacitanceResult;
    return capacitance;
}

pub const CanopyVaporFlux = struct {
    surface_vapor_concentration_m3_per_m3: f64,
    intercepted_evaporation_m3_per_h: f64,
    net_transpiration_into_canopy_m3_per_h: f64,
    latent_heat_flux_mj_per_h: f64,
    convective_water_heat_flux_mj_per_h: f64,
};

pub const RootUptakeGeometry = struct {
    effective_radius_m: f64,
    soil_path_length_m: f64,
    surface_area_per_radius_m: f64,
};

pub fn rootedLayerFraction(
    is_surface_layer: bool,
    layer_top_depth_m: f64,
    layer_thickness_m: f64,
    deepest_primary_root_m: f64,
    seeding_depth_m: f64,
    hypocotyledon_height_m: f64,
) !f64 {
    inline for (.{ layer_top_depth_m, layer_thickness_m, deepest_primary_root_m, seeding_depth_m, hypocotyledon_height_m }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootedLayerInput;
    if (layer_top_depth_m < 0 or layer_thickness_m <= 0 or deepest_primary_root_m < 0 or seeding_depth_m < 0 or hypocotyledon_height_m < 0) return error.InvalidRootedLayerInput;
    if (is_surface_layer) return 1.0;
    const penetration_m = @max(0.0, deepest_primary_root_m - layer_top_depth_m);
    const below_seed_m = @min(layer_thickness_m, penetration_m) - @max(0.0, seeding_depth_m - layer_top_depth_m - hypocotyledon_height_m);
    return @max(0.0, below_seed_m) / layer_thickness_m;
}

/// UPTAKE RRADL, PATH and RTARR. The fallback preserves the source branch for
/// unrooted or zero-fraction layers without evaluating singular expressions.
pub fn rootUptakeGeometry(
    root_length_density_m_per_m3: f64,
    rooted_layer_fraction: f64,
    micropore_volume_fraction: f64,
    aqueous_root_volume_m3: f64,
    root_porosity_fraction: f64,
    plant_population_count: f64,
    root_length_m_per_plant: f64,
    minimum_secondary_radius_m: f64,
    layer_thickness_m: f64,
) !RootUptakeGeometry {
    inline for (.{ root_length_density_m_per_m3, rooted_layer_fraction, micropore_volume_fraction, aqueous_root_volume_m3, root_porosity_fraction, plant_population_count, root_length_m_per_plant, minimum_secondary_radius_m, layer_thickness_m }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootUptakeGeometry;
    if (root_length_density_m_per_m3 < 0 or rooted_layer_fraction < 0 or rooted_layer_fraction > 1 or micropore_volume_fraction < 0 or aqueous_root_volume_m3 < 0 or root_porosity_fraction < 0 or root_porosity_fraction >= 1 or plant_population_count < 0 or root_length_m_per_plant < 0 or minimum_secondary_radius_m <= 0 or layer_thickness_m <= 0) return error.InvalidRootUptakeGeometry;
    if (root_length_density_m_per_m3 <= 0 or rooted_layer_fraction <= 0) return .{
        .effective_radius_m = minimum_secondary_radius_m,
        .soil_path_length_m = layer_thickness_m,
        .surface_area_per_radius_m = 6.283 * root_length_m_per_plant,
    };
    if (micropore_volume_fraction <= 0 or plant_population_count <= 0 or root_length_m_per_plant <= 0) return error.InvalidRootUptakeGeometry;
    const volume_radius = @sqrt((aqueous_root_volume_m3 / (1.0 - root_porosity_fraction)) / (3.1416 * plant_population_count * root_length_m_per_plant));
    return .{
        .effective_radius_m = @max(minimum_secondary_radius_m, volume_radius),
        .soil_path_length_m = 1.0 / @sqrt(3.1416 * (root_length_density_m_per_m3 / rooted_layer_fraction) / micropore_volume_fraction),
        .surface_area_per_radius_m = 6.283 * root_length_m_per_plant / rooted_layer_fraction,
    };
}

pub fn canopyVaporFlux(
    canopy_temperature_k: f64,
    canopy_total_water_potential_mpa: f64,
    aerodynamic_vapor_concentration_m3_per_m3: f64,
    vapor_conductance_m3_per_h: f64,
    intercepted_water_m3: f64,
    minimum_stomatal_resistance_h_per_m: f64,
    cuticular_resistance_h_per_m: f64,
    aerodynamic_plus_boundary_resistance_h_per_m: f64,
    turgor_shape_per_mpa: f64,
    canopy_turgor_potential_mpa: f64,
    water_capacity_deficit_change_m3_per_h: f64,
    latent_heat_of_vaporization_mj_per_m3: f64,
) !CanopyVaporFlux {
    inline for (.{ canopy_temperature_k, canopy_total_water_potential_mpa, aerodynamic_vapor_concentration_m3_per_m3, vapor_conductance_m3_per_h, intercepted_water_m3, minimum_stomatal_resistance_h_per_m, cuticular_resistance_h_per_m, aerodynamic_plus_boundary_resistance_h_per_m, turgor_shape_per_mpa, canopy_turgor_potential_mpa, water_capacity_deficit_change_m3_per_h, latent_heat_of_vaporization_mj_per_m3 }) |value| if (!std.math.isFinite(value)) return error.NonFiniteCanopyVaporInput;
    if (canopy_temperature_k <= 0 or aerodynamic_vapor_concentration_m3_per_m3 < 0 or vapor_conductance_m3_per_h < 0 or intercepted_water_m3 < 0 or minimum_stomatal_resistance_h_per_m < 0 or cuticular_resistance_h_per_m < minimum_stomatal_resistance_h_per_m or aerodynamic_plus_boundary_resistance_h_per_m < 0 or latent_heat_of_vaporization_mj_per_m3 < 0) return error.InvalidCanopyVaporInput;
    const surface_vapor = 2.173e-3 / canopy_temperature_k * 0.61 * @exp(5360.0 * (3.661e-3 - 1.0 / canopy_temperature_k)) * @exp(18.0 * canopy_total_water_potential_mpa / (8.3143 * canopy_temperature_k));
    var vapor_flux = vapor_conductance_m3_per_h * (aerodynamic_vapor_concentration_m3_per_m3 - surface_vapor);
    const intercepted_evaporation = if (vapor_flux > 0) vapor_flux else @max(vapor_flux, -intercepted_water_m3);
    vapor_flux -= intercepted_evaporation;
    const water_stress = @exp(turgor_shape_per_mpa * canopy_turgor_potential_mpa);
    const stomatal_resistance = minimum_stomatal_resistance_h_per_m + (cuticular_resistance_h_per_m - minimum_stomatal_resistance_h_per_m) * water_stress;
    const denominator = aerodynamic_plus_boundary_resistance_h_per_m + stomatal_resistance;
    if (denominator <= 0) return error.InvalidCanopyVaporResistance;
    const transpiration_surface = vapor_flux * (aerodynamic_plus_boundary_resistance_h_per_m + cuticular_resistance_h_per_m) / denominator;
    const transpiration = transpiration_surface + water_capacity_deficit_change_m3_per_h;
    return .{
        .surface_vapor_concentration_m3_per_m3 = surface_vapor,
        .intercepted_evaporation_m3_per_h = intercepted_evaporation,
        .net_transpiration_into_canopy_m3_per_h = transpiration,
        .latent_heat_flux_mj_per_h = (transpiration_surface + intercepted_evaporation) * latent_heat_of_vaporization_mj_per_m3,
        .convective_water_heat_flux_mj_per_h = intercepted_evaporation * 4.19 * canopy_temperature_k,
    };
}

/// UPTAKE RSSX, RSRG, RSR1, RSR2, RSRT and RSRS equations. Counts are
/// expressed per plant, matching the source's RTN*/PP cancellation.
pub fn hydraulicResistance(path: HydraulicPath) !HydraulicResistance {
    inline for (@typeInfo(HydraulicPath).@"struct".fields) |field| if (!std.math.isFinite(@field(path, field.name))) return error.NonFiniteHydraulicPath;
    if (path.soil_path_length_m < 0 or path.root_cylinder_radius_m <= 0 or path.root_surface_area_per_radius_m <= 0 or path.unsaturated_soil_conductivity_m_per_h_mpa <= 0 or path.root_radial_resistivity_mpa_h_per_m3 < 0 or path.soil_air_volume_m3 < 0 or path.soil_water_volume_m3 <= 0 or path.secondary_root_radius_m <= 0 or path.root_length_m_per_plant <= 0 or path.root_axial_resistivity_mpa_h_per_m2 < 0 or path.primary_root_depth_m < 0 or path.primary_root_radius_m <= 0 or path.primary_axis_count_per_plant <= 0 or path.canopy_height_m < 0 or path.stalk_conducting_element_factor <= 0 or path.secondary_root_length_m < 0 or path.secondary_axis_count_per_plant <= 0) return error.InvalidHydraulicPath;
    const soil = @log((path.soil_path_length_m + path.root_cylinder_radius_m) / path.root_cylinder_radius_m) / path.root_surface_area_per_radius_m / path.unsaturated_soil_conductivity_m_per_h_mpa;
    const absorbing_area_m2_per_plant = 6.283 * path.secondary_root_radius_m * path.root_length_m_per_plant;
    const radial = path.root_radial_resistivity_mpa_h_per_m3 / absorbing_area_m2_per_plant * path.soil_air_volume_m3 / path.soil_water_volume_m3;
    const primary_radius_factor = std.math.pow(f64, path.primary_root_radius_m / 0.1e-3, 4);
    const secondary_radius_factor = std.math.pow(f64, path.secondary_root_radius_m / 0.1e-3, 4);
    const primary_axial = path.root_axial_resistivity_mpa_h_per_m2 * path.primary_root_depth_m / (primary_radius_factor * path.primary_axis_count_per_plant) + path.root_axial_resistivity_mpa_h_per_m2 * path.canopy_height_m / (path.stalk_conducting_element_factor * path.primary_axis_count_per_plant);
    const secondary_axial = path.root_axial_resistivity_mpa_h_per_m2 * path.secondary_root_length_m / (secondary_radius_factor * path.secondary_axis_count_per_plant);
    const root = radial + primary_axial + secondary_axial;
    const total = soil + root;
    if (!std.math.isFinite(total) or total <= 0) return error.InvalidHydraulicResistance;
    return .{ .soil_mpa_h_per_m = soil, .root_mpa_h_per_m = root, .total_mpa_h_per_m = total, .conductance_m_per_h_mpa = 1.0 / total };
}

pub const ApplyContext = struct {
    result: *State,
    grid: *const GridState,
    plants: *PlantState,
    soil_total_water_potential_mpa: []const f64,
    active: []const bool,
    root_conductance_m_per_h_mpa: []const f64,
    maximum_uptake_m: []const f64,
    maximum_release_m: []const f64,
    canopy_water_capacitance_m_per_m2_mpa: []const f64,
    transpiration_loss_m: []const f64,
    settings: Settings,
};

const ResidualContext = struct {
    soil_water_potential_mpa: []const f64,
    root_conductance_m_per_h_mpa: []const f64,
    maximum_uptake_m: []const f64,
    maximum_release_m: []const f64,
    previous_canopy_water_potential_mpa: f64,
    canopy_water_capacitance_m_per_m2_mpa: f64,
    transpiration_loss_m: f64,
    active_layer_count: usize = 0,
    layer_capacity: usize = 0,
    root_domain_count: usize = 1,
};

/// Solves each cell/species hydraulic balance directly. The verified uptake.f
/// MXN=200 ceiling is enforced here; there is no surrounding NPH model cycle.
pub fn applyTile(context: *ApplyContext, range: CellRange) !void {
    const result = context.result;
    const plant_count = try std.math.mul(usize, result.cell_count, result.species_count);
    const roots_per_plant = try std.math.mul(usize, result.root_domain_count, result.soil_layer_count);
    const root_count = try std.math.mul(usize, plant_count, roots_per_plant);
    if (range.end > result.cell_count or context.grid.cell_count != result.cell_count or context.plants.cell_count != result.cell_count or context.plants.species_count != result.species_count or context.grid.soil_layer_capacity != result.soil_layer_count or context.plants.soil_layer_count != result.soil_layer_count or context.soil_total_water_potential_mpa.len != context.grid.layer_count or context.active.len != plant_count or context.root_conductance_m_per_h_mpa.len != root_count or context.maximum_uptake_m.len != root_count or context.maximum_release_m.len != root_count or context.canopy_water_capacitance_m_per_m2_mpa.len != plant_count or context.transpiration_loss_m.len != plant_count) return error.PlantWaterDimensionMismatch;
    if (!std.math.isFinite(context.settings.minimum_canopy_water_potential_mpa) or !std.math.isFinite(context.settings.maximum_canopy_water_potential_mpa) or context.settings.minimum_canopy_water_potential_mpa >= context.settings.maximum_canopy_water_potential_mpa or context.settings.maximum_canopy_water_potential_mpa > 0) return error.InvalidCanopyWaterPotentialBounds;

    for (range.first..range.end) |cell| for (0..result.species_count) |species| {
        const plant_index = cell * result.species_count + species;
        const root_base = plant_index * roots_per_plant;
        if (!context.active[plant_index]) {
            @memset(result.root_water_uptake_m[root_base..][0..roots_per_plant], 0);
            result.transpiration_loss_m[plant_index] = 0;
            result.total_root_water_uptake_m[plant_index] = 0;
            result.water_balance_residual_m[plant_index] = 0;
            result.iteration_count[plant_index] = 0;
            result.newton_raphson_step_count[plant_index] = 0;
            result.picard_step_count[plant_index] = 0;
            continue;
        }
        const active_layers = context.grid.active_soil_layer_count[cell];
        const capacitance = context.canopy_water_capacitance_m_per_m2_mpa[plant_index];
        const transpiration = context.transpiration_loss_m[plant_index];
        if (!std.math.isFinite(capacitance) or capacitance <= 0 or !std.math.isFinite(transpiration) or transpiration < 0) return error.InvalidPlantWaterRuntimeInput;
        const residual_context: ResidualContext = .{
            .soil_water_potential_mpa = context.soil_total_water_potential_mpa[cell * result.soil_layer_count ..][0..result.soil_layer_count],
            .root_conductance_m_per_h_mpa = context.root_conductance_m_per_h_mpa[root_base..][0..roots_per_plant],
            .maximum_uptake_m = context.maximum_uptake_m[root_base..][0..roots_per_plant],
            .maximum_release_m = context.maximum_release_m[root_base..][0..roots_per_plant],
            .previous_canopy_water_potential_mpa = context.plants.canopy_water_potential_mpa[plant_index],
            .canopy_water_capacitance_m_per_m2_mpa = capacitance,
            .transpiration_loss_m = transpiration,
            .active_layer_count = active_layers,
            .layer_capacity = result.soil_layer_count,
            .root_domain_count = result.root_domain_count,
        };
        try validateRootInputs(residual_context);
        var options = context.settings.solver_options;
        options.max_iterations = 200;
        options.residual_scale = @max(1.0e-12, transpiration + totalUptake(residual_context, residual_context.previous_canopy_water_potential_mpa));
        const solved = try numerics.newtonPicardFiniteDifference(residual_context, residual, picard, context.settings.minimum_canopy_water_potential_mpa, context.settings.maximum_canopy_water_potential_mpa, residual_context.previous_canopy_water_potential_mpa, options);
        var total_uptake: f64 = 0;
        for (0..result.root_domain_count) |domain| for (0..result.soil_layer_count) |layer| {
            const root = domain * result.soil_layer_count + layer;
            const uptake = if (layer < active_layers) layerUptake(residual_context, root, solved.root) else 0;
            result.root_water_uptake_m[root_base + root] = uptake;
            total_uptake += uptake;
        };
        result.transpiration_loss_m[plant_index] = transpiration;
        result.total_root_water_uptake_m[plant_index] = total_uptake;
        result.water_balance_residual_m[plant_index] = solved.residual;
        result.iteration_count[plant_index] = solved.iterations;
        result.newton_raphson_step_count[plant_index] = solved.newton_raphson_steps;
        result.picard_step_count[plant_index] = solved.picard_steps;
        context.plants.canopy_water_storage_m_per_m2[plant_index] += capacitance * (solved.root - residual_context.previous_canopy_water_potential_mpa);
        context.plants.canopy_water_potential_mpa[plant_index] = solved.root;
    };
}

fn residual(context: ResidualContext, canopy_water_potential_mpa: f64) f64 {
    const storage_change = context.canopy_water_capacitance_m_per_m2_mpa * (canopy_water_potential_mpa - context.previous_canopy_water_potential_mpa);
    return totalUptake(context, canopy_water_potential_mpa) - context.transpiration_loss_m - storage_change;
}

fn picard(context: ResidualContext, canopy_water_potential_mpa: f64) f64 {
    return context.previous_canopy_water_potential_mpa + (totalUptake(context, canopy_water_potential_mpa) - context.transpiration_loss_m) / context.canopy_water_capacitance_m_per_m2_mpa;
}

fn totalUptake(context: ResidualContext, canopy_water_potential_mpa: f64) f64 {
    var total: f64 = 0;
    for (0..context.root_domain_count) |domain| for (0..context.active_layer_count) |layer| {
        total += layerUptake(context, domain * context.layer_capacity + layer, canopy_water_potential_mpa);
    };
    return total;
}

fn layerUptake(context: ResidualContext, root: usize, canopy_water_potential_mpa: f64) f64 {
    const soil_layer = root % context.layer_capacity;
    const unconstrained = context.root_conductance_m_per_h_mpa[root] * (context.soil_water_potential_mpa[soil_layer] - canopy_water_potential_mpa);
    return std.math.clamp(unconstrained, -context.maximum_release_m[root], context.maximum_uptake_m[root]);
}

fn validateRootInputs(context: ResidualContext) !void {
    if (context.active_layer_count == 0 or context.active_layer_count > context.layer_capacity or context.soil_water_potential_mpa.len != context.layer_capacity or context.root_domain_count == 0 or context.root_conductance_m_per_h_mpa.len != context.root_domain_count * context.layer_capacity) return error.PlantWaterDimensionMismatch;
    for (context.root_conductance_m_per_h_mpa, context.maximum_uptake_m, context.maximum_release_m) |conductance, uptake_limit, release_limit| if (!std.math.isFinite(conductance) or conductance < 0 or !std.math.isFinite(uptake_limit) or uptake_limit < 0 or !std.math.isFinite(release_limit) or release_limit < 0) return error.InvalidRootHydraulicInput;
    for (context.soil_water_potential_mpa) |potential| if (!std.math.isFinite(potential) or potential > 0) return error.InvalidRootHydraulicInput;
}

/// Commits the converged hydraulic transaction to STARTQ/UPTAKE root storage.
/// Positive Zig uptake becomes source-sign negative UPWTR (soil water loss).
pub fn commitRootHydraulics(
    result: State,
    roots: *PlantRootState,
    grid: *const GridState,
    soil_total_water_potential_mpa: []const f64,
    plants: *const PlantState,
    active: []const bool,
    cell_area_m2: []const f64,
    soil_resistance_mpa_h_per_m: []const f64,
    root_resistance_mpa_h_per_m: []const f64,
    leaf_osmotic_potential_at_zero_total_mpa: []const f64,
    biological_domain_count_by_plant: []const u8,
) !void {
    const plant_count = try std.math.mul(usize, result.cell_count, result.species_count);
    const root_count = try std.math.mul(usize, try std.math.mul(usize, plant_count, result.root_domain_count), result.soil_layer_count);
    if (roots.plant_count != plant_count or roots.soil_layer_count != result.soil_layer_count or result.root_domain_count != biological_domain_count or grid.cell_count != result.cell_count or soil_total_water_potential_mpa.len != grid.layer_count or plants.canopy_water_potential_mpa.len != plant_count or active.len != plant_count or cell_area_m2.len != result.cell_count or soil_resistance_mpa_h_per_m.len != root_count or root_resistance_mpa_h_per_m.len != root_count or leaf_osmotic_potential_at_zero_total_mpa.len != plant_count or biological_domain_count_by_plant.len != plant_count) return error.PlantWaterDimensionMismatch;
    for (0..plant_count) |plant| {
        if (!active[plant]) continue;
        const cell = plant / result.species_count;
        const area = cell_area_m2[cell];
        if (!std.math.isFinite(area) or area <= 0) return error.InvalidCellArea;
        const biological_domain_count_for_plant = biological_domain_count_by_plant[plant];
        if (biological_domain_count_for_plant < 1 or biological_domain_count_for_plant > result.root_domain_count)
            return error.PlantWaterDimensionMismatch;
        for (0..biological_domain_count_for_plant) |domain| for (0..result.soil_layer_count) |layer| {
            const root = (plant * result.root_domain_count + domain) * result.soil_layer_count + layer;
            const root_state_index = try roots.layerIndex(plant, domain, layer);
            const soil_r = soil_resistance_mpa_h_per_m[root];
            const root_r = root_resistance_mpa_h_per_m[root];
            if (!std.math.isFinite(soil_r) or soil_r < 0 or !std.math.isFinite(root_r) or root_r < 0 or soil_r + root_r <= 0) return error.InvalidRootHydraulicResistance;
            roots.water_uptake_m3_per_h[root_state_index] = -result.root_water_uptake_m[root] * area;
            const soil_potential =
                soil_total_water_potential_mpa[
                    cell * result.soil_layer_count + layer
                ];
            const canopy_potential = plants.canopy_water_potential_mpa[plant];
            const root_potential = @min(0.0, (soil_potential * root_r + canopy_potential * soil_r) / (soil_r + root_r));
            const osmotic = leaf_osmotic_potential_at_zero_total_mpa[plant] + root_potential;
            roots.total_water_potential_mpa[root_state_index] = root_potential;
            roots.osmotic_water_potential_mpa[root_state_index] = osmotic;
            roots.turgor_water_potential_mpa[root_state_index] = @max(0.0, root_potential - osmotic);
        };
    }
    try roots.validateFinite();
}

test "plant water balance converges before uptake MXN ceiling" {
    const context: ResidualContext = .{ .soil_water_potential_mpa = &.{ -0.2, -0.4 }, .root_conductance_m_per_h_mpa = &.{ 0.004, 0.002 }, .maximum_uptake_m = &.{ 0.01, 0.01 }, .maximum_release_m = &.{ 0.01, 0.01 }, .previous_canopy_water_potential_mpa = -1.0, .canopy_water_capacitance_m_per_m2_mpa = 0.002, .transpiration_loss_m = 0.003, .active_layer_count = 2, .layer_capacity = 2 };
    const solved = try numerics.newtonPicardFiniteDifference(context, residual, picard, -10, 0, -1, .{ .absolute_tolerance = 1.0e-10, .relative_tolerance = 1.0e-8, .residual_scale = 0.003, .max_iterations = 200 });
    try std.testing.expect(solved.iterations < 200);
    try std.testing.expect(@abs(solved.residual) < 1.0e-8);
}

test "UPTAKE hydraulic balance sums root and mycorrhizal domains" {
    const context: ResidualContext = .{
        .soil_water_potential_mpa = &.{ -0.2, -0.4 },
        .root_conductance_m_per_h_mpa = &.{ 0.004, 0.002, 0.001, 0.003 },
        .maximum_uptake_m = &.{ 1, 1, 1, 1 },
        .maximum_release_m = &.{ 1, 1, 1, 1 },
        .previous_canopy_water_potential_mpa = -1,
        .canopy_water_capacitance_m_per_m2_mpa = 0.002,
        .transpiration_loss_m = 0,
        .active_layer_count = 2,
        .layer_capacity = 2,
        .root_domain_count = 2,
    };
    try validateRootInputs(context);
    const expected = 0.004 * 0.8 + 0.002 * 0.6 + 0.001 * 0.8 + 0.003 * 0.6;
    try std.testing.expectApproxEqAbs(expected, totalUptake(context, -1), 1e-14);
}

test "UPTAKE soil radial and axial resistances retain source equations" {
    const path: HydraulicPath = .{
        .soil_path_length_m = 0.002,
        .root_cylinder_radius_m = 0.001,
        .root_surface_area_per_radius_m = 0.5,
        .unsaturated_soil_conductivity_m_per_h_mpa = 0.02,
        .root_radial_resistivity_mpa_h_per_m3 = 0.003,
        .soil_air_volume_m3 = 0.6,
        .soil_water_volume_m3 = 0.4,
        .secondary_root_radius_m = 0.0002,
        .root_length_m_per_plant = 2,
        .root_axial_resistivity_mpa_h_per_m2 = 0.004,
        .primary_root_depth_m = 0.3,
        .primary_root_radius_m = 0.0004,
        .primary_axis_count_per_plant = 3,
        .canopy_height_m = 1.2,
        .stalk_conducting_element_factor = 5,
        .secondary_root_length_m = 0.1,
        .secondary_axis_count_per_plant = 8,
    };
    const resistance = try hydraulicResistance(path);
    const expected_soil = @log(3.0) / 0.5 / 0.02;
    const expected_radial = 0.003 / (6.283 * 0.0002 * 2) * 0.6 / 0.4;
    const expected_primary: f64 = 0.004 * 0.3 / (std.math.pow(f64, 4, 4) * 3.0) + 0.004 * 1.2 / (5.0 * 3.0);
    const expected_secondary = 0.004 * 0.1 / (std.math.pow(f64, 2, 4) * 8);
    try std.testing.expectApproxEqAbs(expected_soil, resistance.soil_mpa_h_per_m, 1e-12);
    try std.testing.expectApproxEqAbs(expected_radial + expected_primary + expected_secondary, resistance.root_mpa_h_per_m, 1e-12);
    try std.testing.expectApproxEqAbs(1.0 / resistance.total_mpa_h_per_m, resistance.conductance_m_per_h_mpa, 1e-14);
}

test "UPTAKE canopy water and heat capacity retain source constants" {
    const canopy = try canopyHydraulics(10, 100, 0.01, 2e-5, 3e-5, -1);
    try std.testing.expectApproxEqAbs(16.25, canopy.active_carbon_g, 1e-14);
    const expected_dry_fraction = 0.16 + 0.10 / 2.05;
    try std.testing.expectApproxEqAbs(expected_dry_fraction, canopy.dry_matter_fraction_g_c_per_g, 1e-14);
    try std.testing.expectApproxEqAbs(1e-6 * 16.25 / expected_dry_fraction, canopy.water_capacity_m3, 1e-14);
    try std.testing.expectApproxEqAbs(2.496 * 16.25 * 4e-6 + 4.19 * 5e-5, canopy.wet_heat_capacity_mj_per_k, 1e-14);
    const capacitance = try canopyWaterCapacitanceM3PerMpa(canopy.active_carbon_g, -1);
    const epsilon = 1e-6;
    const wetter = try canopyHydraulics(10, 100, 0.01, 2e-5, 3e-5, -1 + epsilon);
    const d_capacity = (wetter.water_capacity_m3 - canopy.water_capacity_m3) / epsilon;
    try std.testing.expectApproxEqRel(d_capacity, capacitance, 1e-6);
}

test "UPTAKE rooted fraction and radial geometry retain source branches" {
    try std.testing.expectEqual(@as(f64, 1), try rootedLayerFraction(true, 0, 0.02, 0, 0.08, 0));
    try std.testing.expectApproxEqAbs(0.5, try rootedLayerFraction(false, 0.1, 0.2, 0.2, 0.05, 0), 1e-14);
    const geometry = try rootUptakeGeometry(1000, 0.5, 0.4, 1e-6, 0.2, 10, 2, 1e-4, 0.2);
    const expected_radius = @sqrt((1e-6 / 0.8) / (3.1416 * 10 * 2));
    try std.testing.expectApproxEqAbs(@max(1e-4, expected_radius), geometry.effective_radius_m, 1e-14);
    try std.testing.expectApproxEqAbs(1.0 / @sqrt(3.1416 * (1000 / 0.5) / 0.4), geometry.soil_path_length_m, 1e-14);
    try std.testing.expectApproxEqAbs(6.283 * 2 / 0.5, geometry.surface_area_per_radius_m, 1e-14);
    const fallback = try rootUptakeGeometry(0, 0, 0, 0, 0.2, 0, 0, 2.5e-6, 0.3);
    try std.testing.expectEqual(@as(f64, 0.3), fallback.soil_path_length_m);
    try std.testing.expectEqual(@as(f64, 2.5e-6), fallback.effective_radius_m);
}

test "UPTAKE canopy vapor flux limits intercepted evaporation before transpiration" {
    const flux = try canopyVaporFlux(290, -1, 0, 10, 1e-5, 0.01, 0.1, 0.02, -1, 1, 0, 2450);
    try std.testing.expectApproxEqAbs(-1e-5, flux.intercepted_evaporation_m3_per_h, 1e-14);
    try std.testing.expect(flux.net_transpiration_into_canopy_m3_per_h < 0);
    try std.testing.expectApproxEqAbs((flux.net_transpiration_into_canopy_m3_per_h + flux.intercepted_evaporation_m3_per_h) * 2450, flux.latent_heat_flux_mj_per_h, 1e-12);
    try std.testing.expectApproxEqAbs(flux.intercepted_evaporation_m3_per_h * 4.19 * 290, flux.convective_water_heat_flux_mj_per_h, 1e-12);
}

test "UPTAKE commit preserves source sign and resistance-weighted root potential" {
    const config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 40 });
    var grid = try GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    grid.matric_potential_mpa[0] = -0.2;
    grid.matric_potential_mpa[1] = -0.4;
    var plants = try PlantState.init(std.testing.allocator, config);
    defer plants.deinit();
    plants.canopy_water_potential_mpa[0] = -1.0;
    var balance = try State.init(std.testing.allocator, 1, 1, 2);
    defer balance.deinit();
    balance.root_water_uptake_m[0] = 0.001;
    balance.root_water_uptake_m[1] = 0.002;
    balance.root_water_uptake_m[2] = 0.003;
    balance.root_water_uptake_m[3] = 0.004;
    var roots = try PlantRootState.init(std.testing.allocator, 1, 2, 3);
    defer roots.deinit();
    try commitRootHydraulics(balance, &roots, &grid, &.{ -0.25, -0.5 }, &plants, &.{true}, &.{10}, &.{ 1, 1, 1, 1 }, &.{ 3, 3, 3, 3 }, &.{-1.5}, &.{2});
    try std.testing.expectApproxEqAbs(-0.01, roots.water_uptake_m3_per_h[try roots.layerIndex(0, 0, 0)], 1e-14);
    try std.testing.expectApproxEqAbs(-0.4375, roots.total_water_potential_mpa[try roots.layerIndex(0, 0, 0)], 1e-14);
    try std.testing.expectApproxEqAbs(-1.9375, roots.osmotic_water_potential_mpa[try roots.layerIndex(0, 0, 0)], 1e-14);
    try std.testing.expectApproxEqAbs(1.5, roots.turgor_water_potential_mpa[try roots.layerIndex(0, 0, 0)], 1e-14);
}

test "UPTAKE PSILZ retains the daily minimum and resets before the next day" {
    const config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 2 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 40 });
    var plants = try PlantState.init(std.testing.allocator, config);
    defer plants.deinit();
    var canopy = try CanopyState.init(std.testing.allocator, 1, 2, &.{ 1, 1 }, &.{ 1, 1 }, &.{ 1, 1 });
    defer canopy.deinit();
    plants.canopy_water_potential_mpa[0] = -0.5;
    plants.canopy_water_potential_mpa[1] = -1.0;
    try updateDailyMinimumCanopyWaterPotential(&canopy, &plants);
    plants.canopy_water_potential_mpa[0] = -1.5;
    plants.canopy_water_potential_mpa[1] = -0.25;
    try updateDailyMinimumCanopyWaterPotential(&canopy, &plants);
    try std.testing.expectEqualSlices(f64, &.{ -1.5, -1.0 }, canopy.plant_minimum_daily_canopy_water_potential_mpa);
    var first_plant_by_cell = [_]f64{0};
    try publishFirstPlantDailyMinimumByCell(&canopy, &first_plant_by_cell);
    try std.testing.expectEqual(@as(f64, -1.5), first_plant_by_cell[0]);
    resetDailyMinimumCanopyWaterPotential(&canopy);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0 }, canopy.plant_minimum_daily_canopy_water_potential_mpa);
}

test "hourly hydraulic workspace skips plants without live root conductance" {
    const config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 40 });
    var grid = try GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var plants = try PlantState.init(std.testing.allocator, config);
    defer plants.deinit();
    var canopy = try CanopyState.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{1});
    defer canopy.deinit();
    var workspace = try Workspace.init(std.testing.allocator, 1, 1, 2);
    defer workspace.deinit();
    workspace.cell_area_m2[0] = 10;
    try refreshCanopyWorkspace(&workspace, &canopy, &plants, 1);
    try workspace.refreshActive(2);
    try std.testing.expect(!workspace.active[0]);
    var balance = try State.init(std.testing.allocator, 1, 1, 2);
    defer balance.deinit();
    var context: ApplyContext = .{ .result = &balance, .grid = &grid, .plants = &plants, .soil_total_water_potential_mpa = grid.matric_potential_mpa, .active = workspace.active, .root_conductance_m_per_h_mpa = workspace.root_conductance_m_per_h_mpa, .maximum_uptake_m = workspace.maximum_uptake_m, .maximum_release_m = workspace.maximum_release_m, .canopy_water_capacitance_m_per_m2_mpa = workspace.canopy_water_capacitance_m_per_m2_mpa, .transpiration_loss_m = workspace.transpiration_loss_m, .settings = .{ .minimum_canopy_water_potential_mpa = -100, .maximum_canopy_water_potential_mpa = 0, .solver_options = .{ .max_iterations = 40 } } };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expectEqual(@as(u16, 0), balance.iteration_count[0]);
    try std.testing.expectEqual(@as(f64, 0), balance.total_root_water_uptake_m[0]);
}

test "hourly workspace builds live root and mycorrhizal conductance without allocation" {
    const config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 40 });
    var grid = try GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    grid.matrix_liquid_water_m3[0] = 0.2;
    grid.matric_potential_mpa[0] = -0.2;
    var plants = try PlantState.init(std.testing.allocator, config);
    defer plants.deinit();
    plants.canopy_water_potential_mpa[0] = -1;
    var canopy = try CanopyState.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{1});
    defer canopy.deinit();
    canopy.branch_leaf_carbon_g[0] = 10;
    var roots = try PlantRootState.init(std.testing.allocator, 1, 1, 3);
    defer roots.deinit();
    const traits = try @import("plant_traits.zig").parse(@import("test_fixtures.zig").plant_traits_source);
    try roots.initializePlant(0, traits, 0, 0.02, @import("plant_root_system.zig").compatibilityInitializationParameters());
    for (0..biological_domain_count) |domain| {
        const axis = try roots.layerAxisIndex(0, domain, 0, 0);
        roots.axis_primary_count[axis] = 10;
        roots.axis_secondary_count[axis] = 10;
        roots.axis_primary_length_m[axis] = 0.1;
        roots.axis_secondary_length_m[axis] = 1;
        roots.axis_depth_m[try roots.axisIndex(0, domain, 0)] = 0.1;
        roots.total_carbon_g[try roots.layerIndex(0, domain, 0)] = 1;
    }
    var workspace = try Workspace.init(std.testing.allocator, 1, 1, 1);
    defer workspace.deinit();
    workspace.cell_area_m2[0] = 10;
    workspace.plant_population_count[0] = 10;
    workspace.root_porosity_fraction[0] = traits.roots.root_porosity_fraction;
    workspace.root_radial_resistivity_mpa_h_per_m3[0] = traits.roots.radial_resistivity_mpa_h_per_m3;
    workspace.root_axial_resistivity_mpa_h_per_m2[0] = traits.roots.axial_resistivity_mpa_h_per_m2;
    workspace.seeding_depth_m[0] = 0.02;
    workspace.woody_root_fraction[0] = 1;
    const curves: []const @import("soil_water_retention.zig").ResolvedCurve = &.{};
    const one = [_]f64{1};
    const half = [_]f64{0.5};
    const tenth = [_]f64{0.1};
    const conductivity = [_]f64{ 0.01, 0.01, 0.01 };
    const properties: SoilPropertiesState = .{ .allocator = undefined, .layer_count = 1, .hydraulic_conductivity_class_count = 1, .retention_curve = @constCast(curves), .matrix_bulk_volume_m3 = @constCast(&half), .layer_volume_m3 = @constCast(&one), .layer_thickness_m = @constCast(&tenth), .layer_midpoint_depth_m = @constCast(&half), .layer_bottom_depth_m = @constCast(&one), .bulk_density_megagrams_per_m3 = @constCast(&one), .sand_mass_fraction = @constCast(&half), .clay_mass_fraction = @constCast(&half), .sand_mass_Mg = @constCast(&half), .silt_mass_Mg = @constCast(&half), .clay_mass_Mg = @constCast(&half), .total_organic_carbon_g_per_megagram = @constCast(&one), .cation_exchange_capacity_mol_per_Mg = @constCast(&one), .anion_exchange_capacity_mol_per_Mg = @constCast(&one), .cation_exchange_capacity_mol = @constCast(&one), .anion_exchange_capacity_mol = @constCast(&one), .porosity_fraction = @constCast(&half), .matrix_air_entry_water_fraction = @constCast(&half), .saturation_water_potential_mpa = @constCast(&tenth), .rainfall_conductivity_multiplier = @constCast(&one), .matrix_hydraulic_conductivity_m2_per_h_mpa = @constCast(&conductivity) };
    try refreshCanopyWorkspace(&workspace, &canopy, &plants, 1);
    try refreshRootWorkspace(&workspace, &roots, &canopy, &plants, &grid, &properties, 1, &.{2}, 1.0e-6, 0.05, 3.142, @import("plant_root_system.zig").compatibilityMorphologyParameters());
    try workspace.refreshActive(1);
    try std.testing.expect(workspace.active[0]);
    try std.testing.expect(workspace.root_conductance_m_per_h_mpa[0] > 0);
    try std.testing.expect(workspace.root_conductance_m_per_h_mpa[1] > 0);
    try std.testing.expect(workspace.maximum_uptake_m[0] > 0);
}
