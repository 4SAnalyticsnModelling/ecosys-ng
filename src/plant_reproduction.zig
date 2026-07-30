const std = @import("std");
const CellRange = @import("compute.zig").CellRange;
const PlantState = @import("grid.zig").PlantState;
const canopy_module = @import("canopy_photosynthesis.zig");
const stages_module = @import("plant_growth_stages.zig");

pub const Controls = struct {
    allocator: std.mem.Allocator,
    potential_sites_per_g_growth: []f64,
    maximum_seeds_per_site: []f64,
    maximum_individual_seed_carbon_g: []f64,
    chilling_temperature_c: []f64,
    stomatal_turgor_shape: []f64,
    shallow_root_profile: []bool,

    pub fn init(allocator: std.mem.Allocator, plant_count: usize) !Controls {
        if (plant_count == 0) return error.InvalidPlantReproductionDimensions;
        var result: Controls = undefined;
        result.allocator = allocator;
        var allocated: usize = 0;
        errdefer inline for (@typeInfo(Controls).@"struct".fields) |field| if (field.type == []f64 and allocated > 0) {
            allocated -= 1;
            allocator.free(@field(result, field.name));
        };
        inline for (@typeInfo(Controls).@"struct".fields) |field| if (field.type == []f64) {
            @field(result, field.name) = try allocator.alloc(f64, plant_count);
            @memset(@field(result, field.name), 0);
            allocated += 1;
        };
        result.shallow_root_profile = try allocator.alloc(bool, plant_count);
        @memset(result.shallow_root_profile, false);
        return result;
    }

    pub fn deinit(self: *Controls) void {
        self.allocator.free(self.shallow_root_profile);
        inline for (@typeInfo(Controls).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }

    pub fn setPlant(self: *Controls, plant: usize, potential_sites_per_g_growth: f64, maximum_seeds_per_site: f64, maximum_individual_seed_carbon_g: f64, chilling_temperature_c: f64, stomatal_turgor_shape: f64, shallow_root_profile: bool) !void {
        if (plant >= self.potential_sites_per_g_growth.len) return error.PlantReproductionIndexOutOfBounds;
        inline for (.{ potential_sites_per_g_growth, maximum_seeds_per_site, maximum_individual_seed_carbon_g }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPlantReproductionControl;
        // OSTGX is a signed exponential coefficient; supplied PFTs commonly
        // use negative values (for example maize uses -5 MPa-1).
        if (!std.math.isFinite(stomatal_turgor_shape)) return error.InvalidPlantReproductionControl;
        if (!std.math.isFinite(chilling_temperature_c)) return error.InvalidPlantReproductionControl;
        self.potential_sites_per_g_growth[plant] = potential_sites_per_g_growth;
        self.maximum_seeds_per_site[plant] = maximum_seeds_per_site;
        self.maximum_individual_seed_carbon_g[plant] = maximum_individual_seed_carbon_g;
        self.chilling_temperature_c[plant] = chilling_temperature_c;
        self.stomatal_turgor_shape[plant] = stomatal_turgor_shape;
        self.shallow_root_profile[plant] = shallow_root_profile;
    }
};

pub const ApplyContext = struct {
    canopy: *canopy_module.State,
    plants: *const PlantState,
    growth_stages: *const stages_module.State,
    controls: *const Controls,
    active_by_plant: []const bool,
    minimum_turgor_potential_mpa: f64,
    seed_set_parameters: canopy_module.SeedSetParameters,
    structural_presence_threshold_g_per_plant: f64,
    timestep_h: f64,
};

pub fn applyTile(context: *ApplyContext, range: CellRange) !void {
    try context.seed_set_parameters.validate();
    const canopy = context.canopy;
    const plant_count = canopy.plant_branch_offsets.len - 1;
    if (range.end > canopy.cell_count or context.plants.cell_count * context.plants.species_count != plant_count or context.growth_stages.plant_count != plant_count or context.active_by_plant.len != plant_count or context.controls.potential_sites_per_g_growth.len != plant_count) return error.PlantReproductionDimensionMismatch;
    for (range.first..range.end) |cell| for (0..canopy.species_count) |species| {
        const plant = cell * canopy.species_count + species;
        if (!context.active_by_plant[plant]) continue;
        const structural_presence_threshold_g = try plantScaledPresenceThresholdG(
            context.structural_presence_threshold_g_per_plant,
            canopy.plant_population_count[plant],
        );
        const response = try canopy_module.canopyWaterGrowthResponse(context.controls.shallow_root_profile[plant], canopy.plant_canopy_turgor_potential_mpa[plant], context.minimum_turgor_potential_mpa, context.plants.canopy_water_potential_mpa[plant], context.controls.stomatal_turgor_shape[plant]);
        const branches = try canopy.branchRange(plant);
        for (branches.first..branches.end) |branch| {
            const stage = context.growth_stages.branches[branch];
            if (stage.dead) continue;
            const branch_shoot_carbon_g = try branchShootCarbon(canopy, branch);
            canopy.branch_potential_seed_site_count[branch] = try canopy_module.accumulatePotentialSeedSites(canopy.branch_potential_seed_site_count[branch], stage.stem_elongation_start_day != 0, stage.anthesis_day != 0, branch_shoot_carbon_g, canopy.plant_total_shoot_carbon_g[plant], @max(0.0, canopy.plant_shoot_growth_g_c_per_step[plant]), context.controls.potential_sites_per_g_growth[plant], structural_presence_threshold_g);
            const result = try canopy_module.updateSeedNumberAndSize(.{ .anthesis_started = stage.anthesis_day != 0, .grain_fill_started = stage.grain_fill_start_day != 0, .final_seed_number_set = stage.seed_number_set_end_day != 0, .maximum_seed_size_set = stage.seed_size_set_end_day != 0, .mobile_carbon_concentration_g_per_g = canopy.branch_mobile_carbon_concentration_g_per_g[branch], .mobile_nitrogen_concentration_g_per_g = canopy.branch_mobile_nitrogen_concentration_g_per_g[branch], .mobile_phosphorus_concentration_g_per_g = canopy.branch_mobile_phosphorus_concentration_g_per_g[branch], .carbon_half_saturation_g_per_g = context.seed_set_parameters.carbon_half_saturation_g_per_g, .nitrogen_half_saturation_g_per_g = context.seed_set_parameters.nitrogen_half_saturation_g_per_g, .phosphorus_half_saturation_g_per_g = context.seed_set_parameters.phosphorus_half_saturation_g_per_g, .canopy_temperature_c = context.plants.canopy_temperature_k[plant] - 273.15, .chilling_temperature_c = context.controls.chilling_temperature_c[plant], .high_temperature_c = canopy.plant_seed_set_high_temperature_c[plant], .seed_loss_fraction_per_c_h = canopy.plant_seed_set_loss_fraction_per_c_h[plant], .timestep_h = context.timestep_h, .water_growth_fraction = response.growth_fraction, .reproductive_stage_increment = stage.reproductive_stage_increment, .maximum_seeds_per_site = context.controls.maximum_seeds_per_site[plant], .potential_site_count = canopy.branch_potential_seed_site_count[branch], .current_seed_count = canopy.branch_seed_count[branch], .maximum_individual_seed_carbon_g = context.controls.maximum_individual_seed_carbon_g[plant], .current_individual_seed_carbon_g = canopy.branch_individual_seed_carbon_g[branch] });
            canopy.branch_seed_count[branch] = result.seed_count;
            canopy.branch_individual_seed_carbon_g[branch] = result.individual_seed_carbon_g;
        }
    };
}

fn plantScaledPresenceThresholdG(threshold_g_per_plant: f64, plant_population_count: f64) !f64 {
    if (!std.math.isFinite(threshold_g_per_plant) or threshold_g_per_plant < 0 or
        !std.math.isFinite(plant_population_count) or plant_population_count < 0)
    {
        return error.InvalidPlantReproductionPresenceThreshold;
    }
    const threshold_g = threshold_g_per_plant * plant_population_count;
    if (!std.math.isFinite(threshold_g)) return error.NonFinitePlantReproductionPresenceThreshold;
    return threshold_g;
}

test "GROSUB ZEROP scales the per-plant presence threshold by current population" {
    try std.testing.expectEqual(
        @as(f64, 1.0e-15) * @as(f64, 300.0),
        try plantScaledPresenceThresholdG(1.0e-15, 300.0),
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        try plantScaledPresenceThresholdG(1.0e-15, 0),
    );
    try std.testing.expectError(
        error.InvalidPlantReproductionPresenceThreshold,
        plantScaledPresenceThresholdG(1.0e-15, -1),
    );
}

fn branchShootCarbon(canopy: *const canopy_module.State, branch: usize) !f64 {
    if (branch >= canopy.branch_mobile_carbon_g.len) return error.CanopyBranchIndexOutOfBounds;
    var total = canopy.branch_leaf_carbon_g[branch] + canopy.branch_sheath_carbon_g[branch] + canopy.branch_stalk_carbon_g[branch] + canopy.branch_reserve_carbon_g[branch] + canopy.branch_husk_carbon_g[branch] + canopy.branch_ear_carbon_g[branch] + canopy.branch_grain_carbon_g[branch] + canopy.branch_mobile_carbon_g[branch];
    const nodes = try canopy.nodeRange(branch);
    for (nodes.first..nodes.end) |node| total += canopy.node_c3_nonstructural_carbon_g[node] + canopy.node_c4_mesophyll_nonstructural_carbon_g[node] + canopy.node_bundle_sheath_co2_carbon_g[node] + canopy.node_bundle_sheath_bicarbonate_carbon_g[node];
    if (!std.math.isFinite(total) or total < 0) return error.InvalidBranchShootCarbon;
    return total;
}

test "GROSUB reproduction kernel advances arbitrary runtime plants" {
    const allocator = std.testing.allocator;
    var canopy = try canopy_module.State.init(allocator, 1, 1, &.{1}, &.{1}, &.{1});
    defer canopy.deinit();
    var plants = try PlantState.init(allocator, try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 }));
    defer plants.deinit();
    var stages = try stages_module.State.init(allocator, &.{1});
    defer stages.deinit();
    var controls = try Controls.init(allocator, 1);
    defer controls.deinit();
    try controls.setPlant(0, 3, 4, 0.5, 5, 1, false);
    canopy.branch_leaf_carbon_g[0] = 20;
    canopy.plant_total_shoot_carbon_g[0] = 20;
    canopy.plant_shoot_growth_g_c_per_step[0] = 2;
    canopy.branch_mobile_carbon_concentration_g_per_g[0] = 0.1;
    canopy.branch_mobile_nitrogen_concentration_g_per_g[0] = 0.02;
    canopy.branch_mobile_phosphorus_concentration_g_per_g[0] = 0.004;
    canopy.plant_seed_set_high_temperature_c[0] = 40;
    canopy.plant_seed_set_loss_fraction_per_c_h[0] = 0.01;
    canopy.plant_canopy_turgor_potential_mpa[0] = 0.2;
    plants.canopy_water_potential_mpa[0] = -0.5;
    plants.canopy_temperature_k[0] = 298.15;
    stages.branches[0].stem_elongation_start_day = 1;
    var context: ApplyContext = .{ .canopy = &canopy, .plants = &plants, .growth_stages = &stages, .controls = &controls, .active_by_plant = &.{true}, .minimum_turgor_potential_mpa = 0.1, .seed_set_parameters = canopy_module.compatibilitySeedSetParameters(), .structural_presence_threshold_g_per_plant = 1.0e-12, .timestep_h = 1 };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expectApproxEqAbs(@as(f64, 6), canopy.branch_potential_seed_site_count[0], 1.0e-14);
}

test "source negative stomatal turgor coefficient is retained" {
    var controls = try Controls.init(std.testing.allocator, 1);
    defer controls.deinit();
    try controls.setPlant(0, 1.2, 6, 0.2, -1, -5, false);
    try std.testing.expectEqual(@as(f64, -5), controls.stomatal_turgor_shape[0]);
}
