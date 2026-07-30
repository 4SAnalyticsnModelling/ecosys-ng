const std = @import("std");
const canopy_module = @import("canopy_photosynthesis.zig");
const stages_module = @import("plant_growth_stages.zig");
const structure_module = @import("canopy_structure.zig");
const PlantState = @import("grid.zig").PlantState;

pub const Controls = struct {
    allocator: std.mem.Allocator,
    leaf_length_to_width_ratio: []f64,
    stalk_volume_m3_per_g_c: []f64,
    biomass_turnover_type: []u8,
    root_profile_type: []u8,
    annual_growth_habit: []bool,
    stem_angle_sine: []f64,

    pub fn init(allocator: std.mem.Allocator, plant_count: usize) !Controls {
        if (plant_count == 0) return error.InvalidCanopyLayerDimensions;
        const ratio = try allocator.alloc(f64, plant_count);
        errdefer allocator.free(ratio);
        const volume = try allocator.alloc(f64, plant_count);
        errdefer allocator.free(volume);
        const turnover = try allocator.alloc(u8, plant_count);
        errdefer allocator.free(turnover);
        const profile = try allocator.alloc(u8, plant_count);
        errdefer allocator.free(profile);
        const annual = try allocator.alloc(bool, plant_count);
        errdefer allocator.free(annual);
        const stem_angle_sine = try allocator.alloc(f64, plant_count);
        @memset(ratio, 0);
        @memset(volume, 0);
        @memset(turnover, 0);
        @memset(profile, 0);
        @memset(annual, true);
        @memset(stem_angle_sine, 1);
        return .{ .allocator = allocator, .leaf_length_to_width_ratio = ratio, .stalk_volume_m3_per_g_c = volume, .biomass_turnover_type = turnover, .root_profile_type = profile, .annual_growth_habit = annual, .stem_angle_sine = stem_angle_sine };
    }

    pub fn deinit(self: *Controls) void {
        self.allocator.free(self.stem_angle_sine);
        self.allocator.free(self.annual_growth_habit);
        self.allocator.free(self.root_profile_type);
        self.allocator.free(self.biomass_turnover_type);
        self.allocator.free(self.stalk_volume_m3_per_g_c);
        self.allocator.free(self.leaf_length_to_width_ratio);
        self.* = undefined;
    }

    pub fn setPlant(self: *Controls, plant: usize, leaf_length_to_width_ratio: f64, stalk_volume_m3_per_g_c: f64, biomass_turnover_type: u8, root_profile_type: u8, annual_growth_habit: bool, stem_angle_sine: f64) !void {
        if (plant >= self.leaf_length_to_width_ratio.len) return error.CanopyLayerPlantIndexOutOfBounds;
        if (!std.math.isFinite(leaf_length_to_width_ratio) or leaf_length_to_width_ratio < 0 or !std.math.isFinite(stalk_volume_m3_per_g_c) or stalk_volume_m3_per_g_c <= 0 or !std.math.isFinite(stem_angle_sine) or stem_angle_sine < 0 or stem_angle_sine > 1) return error.InvalidCanopyLayerControl;
        self.leaf_length_to_width_ratio[plant] = leaf_length_to_width_ratio;
        self.stalk_volume_m3_per_g_c[plant] = stalk_volume_m3_per_g_c;
        self.biomass_turnover_type[plant] = biomass_turnover_type;
        self.root_profile_type[plant] = root_profile_type;
        self.annual_growth_habit[plant] = annual_growth_habit;
        self.stem_angle_sine[plant] = stem_angle_sine;
    }
};

/// Heap-owned ARLFL/WGLFL/WGLFLN/WGLFLP/ARSTK and ZL state. Node and branch
/// domains follow dynamic shoot topology; layer count is selected at runtime.
pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    species_count: usize,
    layer_count: usize,
    inclination_count: usize,
    azimuth_count: usize,
    node_count: usize,
    branch_count: usize,
    boundary_height_m: []f64,
    node_leaf_area_m2: []f64,
    node_leaf_carbon_g: []f64,
    node_leaf_nitrogen_g: []f64,
    node_leaf_phosphorus_g: []f64,
    node_leaf_projected_surface_m2: []f64,
    branch_stalk_area_m2: []f64,
    branch_stalk_projected_surface_m2: []f64,
    plant_leaf_projected_surface_m2: []f64,
    plant_stalk_projected_surface_m2: []f64,
    plant_standing_dead_area_m2: []f64,
    plant_standing_dead_projected_surface_m2: []f64,
    cell_leaf_area_m2: []f64,
    cell_leaf_carbon_g: []f64,
    cell_stalk_area_m2: []f64,
    cell_standing_dead_area_m2: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, species_count: usize, layer_count: usize, inclination_count: usize, azimuth_count: usize, canopy: *const canopy_module.State) !State {
        if (cell_count == 0 or species_count == 0 or layer_count == 0 or inclination_count == 0 or azimuth_count == 0 or canopy.cell_count != cell_count or canopy.species_count != species_count) return error.InvalidCanopyLayerDimensions;
        const node_count = canopy.node_sample_offsets.len - 1;
        const branch_count = canopy.branch_node_offsets.len - 1;
        var result: State = .{ .allocator = allocator, .cell_count = cell_count, .species_count = species_count, .layer_count = layer_count, .inclination_count = inclination_count, .azimuth_count = azimuth_count, .node_count = node_count, .branch_count = branch_count, .boundary_height_m = undefined, .node_leaf_area_m2 = undefined, .node_leaf_carbon_g = undefined, .node_leaf_nitrogen_g = undefined, .node_leaf_phosphorus_g = undefined, .node_leaf_projected_surface_m2 = undefined, .branch_stalk_area_m2 = undefined, .branch_stalk_projected_surface_m2 = undefined, .plant_leaf_projected_surface_m2 = undefined, .plant_stalk_projected_surface_m2 = undefined, .plant_standing_dead_area_m2 = undefined, .plant_standing_dead_projected_surface_m2 = undefined, .cell_leaf_area_m2 = undefined, .cell_leaf_carbon_g = undefined, .cell_stalk_area_m2 = undefined, .cell_standing_dead_area_m2 = undefined };
        var allocated: usize = 0;
        errdefer inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64 and allocated > 0) {
            allocated -= 1;
            allocator.free(@field(result, field.name));
        };
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            const count = try fieldCount(field.name, cell_count, species_count, layer_count, inclination_count, node_count, branch_count);
            @field(result, field.name) = try allocator.alloc(f64, count);
            @memset(@field(result, field.name), 0);
            allocated += 1;
        };
        for (0..cell_count) |cell| result.boundary_height_m[cell * (layer_count + 1) + layer_count] = 0.01;
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }

    pub fn ensureTopology(self: *State, canopy: *const canopy_module.State) !void {
        const node_count = canopy.node_sample_offsets.len - 1;
        const branch_count = canopy.branch_node_offsets.len - 1;
        if (canopy.cell_count != self.cell_count or canopy.species_count != self.species_count) return error.CanopyLayerTopologyMismatch;
        if (node_count == self.node_count and branch_count == self.branch_count) return;
        var replacement = try State.init(self.allocator, self.cell_count, self.species_count, self.layer_count, self.inclination_count, self.azimuth_count, canopy);
        errdefer replacement.deinit();
        @memcpy(replacement.boundary_height_m, self.boundary_height_m);
        @memcpy(replacement.cell_standing_dead_area_m2, self.cell_standing_dead_area_m2);
        var previous = self.*;
        self.* = replacement;
        previous.deinit();
    }

    pub fn resetCurrentAreas(self: *State) void {
        @memset(self.node_leaf_area_m2, 0);
        @memset(self.node_leaf_carbon_g, 0);
        @memset(self.node_leaf_nitrogen_g, 0);
        @memset(self.node_leaf_phosphorus_g, 0);
        @memset(self.node_leaf_projected_surface_m2, 0);
        @memset(self.branch_stalk_area_m2, 0);
        @memset(self.branch_stalk_projected_surface_m2, 0);
        @memset(self.plant_leaf_projected_surface_m2, 0);
        @memset(self.plant_stalk_projected_surface_m2, 0);
        @memset(self.plant_standing_dead_area_m2, 0);
        @memset(self.plant_standing_dead_projected_surface_m2, 0);
        @memset(self.cell_leaf_area_m2, 0);
        @memset(self.cell_leaf_carbon_g, 0);
        @memset(self.cell_stalk_area_m2, 0);
        @memset(self.cell_standing_dead_area_m2, 0);
    }

    pub fn cellBoundaries(self: State, cell: usize) ![]f64 {
        if (cell >= self.cell_count) return error.CanopyLayerCellIndexOutOfBounds;
        const first = cell * (self.layer_count + 1);
        return self.boundary_height_m[first .. first + self.layer_count + 1];
    }

    pub fn nodeLayerRange(self: State, node: usize) !struct { first: usize, end: usize } {
        if (node >= self.node_count) return error.CanopyLayerNodeIndexOutOfBounds;
        return .{ .first = node * self.layer_count, .end = (node + 1) * self.layer_count };
    }

    pub fn branchLayerRange(self: State, branch: usize) !struct { first: usize, end: usize } {
        if (branch >= self.branch_count) return error.CanopyLayerBranchIndexOutOfBounds;
        return .{ .first = branch * self.layer_count, .end = (branch + 1) * self.layer_count };
    }

    /// Late-GROSUB dead-branch clearing. Removes the branch from every
    /// node/branch, plant, and cell layer aggregate without waiting for the
    /// next HOUR1 distribution refresh.
    pub fn clearDeadBranch(self: *State, canopy: *const canopy_module.State, plant: usize, branch: usize) !void {
        if (plant >= self.cell_count * self.species_count or branch >= self.branch_count or
            canopy.node_sample_offsets.len - 1 != self.node_count)
            return error.CanopyLayerDistributionDimensionMismatch;
        const plant_branches = try canopy.branchRange(plant);
        if (branch < plant_branches.first or branch >= plant_branches.end) return error.CanopyBranchIndexOutOfBounds;
        const cell = plant / self.species_count;
        const nodes = try canopy.nodeRange(branch);
        for (nodes.first..nodes.end) |node| {
            for (0..self.layer_count) |layer| {
                const node_layer = node * self.layer_count + layer;
                const cell_layer = cell * self.layer_count + layer;
                self.cell_leaf_area_m2[cell_layer] = @max(0, self.cell_leaf_area_m2[cell_layer] - self.node_leaf_area_m2[node_layer]);
                self.cell_leaf_carbon_g[cell_layer] = @max(0, self.cell_leaf_carbon_g[cell_layer] - self.node_leaf_carbon_g[node_layer]);
                self.node_leaf_area_m2[node_layer] = 0;
                self.node_leaf_carbon_g[node_layer] = 0;
                self.node_leaf_nitrogen_g[node_layer] = 0;
                self.node_leaf_phosphorus_g[node_layer] = 0;
                for (0..self.inclination_count) |inclination| {
                    const node_surface = node_layer * self.inclination_count + inclination;
                    const plant_surface = (plant * self.layer_count + layer) * self.inclination_count + inclination;
                    self.plant_leaf_projected_surface_m2[plant_surface] = @max(
                        0,
                        self.plant_leaf_projected_surface_m2[plant_surface] - self.node_leaf_projected_surface_m2[node_surface],
                    );
                    self.node_leaf_projected_surface_m2[node_surface] = 0;
                }
            }
        }
        for (0..self.layer_count) |layer| {
            const branch_layer = branch * self.layer_count + layer;
            const cell_layer = cell * self.layer_count + layer;
            self.cell_stalk_area_m2[cell_layer] = @max(0, self.cell_stalk_area_m2[cell_layer] - self.branch_stalk_area_m2[branch_layer]);
            self.branch_stalk_area_m2[branch_layer] = 0;
            for (0..self.inclination_count) |inclination| {
                const branch_surface = branch_layer * self.inclination_count + inclination;
                const plant_surface = (plant * self.layer_count + layer) * self.inclination_count + inclination;
                self.plant_stalk_projected_surface_m2[plant_surface] = @max(
                    0,
                    self.plant_stalk_projected_surface_m2[plant_surface] - self.branch_stalk_projected_surface_m2[branch_surface],
                );
                self.branch_stalk_projected_surface_m2[branch_surface] = 0;
            }
        }
    }

    /// ARLFS-equivalent total live leaf area for one runtime plant population.
    pub fn plantLeafAreaM2(self: State, canopy: *const canopy_module.State, plant: usize) !f64 {
        if (plant >= self.cell_count * self.species_count or canopy.cell_count != self.cell_count or canopy.species_count != self.species_count or canopy.node_sample_offsets.len - 1 != self.node_count) return error.CanopyLayerDistributionDimensionMismatch;
        const branches = try canopy.branchRange(plant);
        var area_m2: f64 = 0;
        for (branches.first..branches.end) |branch| {
            const nodes = try canopy.nodeRange(branch);
            for (nodes.first..nodes.end) |node| {
                const layers = try self.nodeLayerRange(node);
                for (self.node_leaf_area_m2[layers.first..layers.end]) |area| area_m2 += area;
            }
        }
        if (!std.math.isFinite(area_m2) or area_m2 < 0) return error.NonFiniteCanopyLayerDistribution;
        return area_m2;
    }

    /// Publishes HOUR1 ARLFS/AREA into the authoritative runtime plant LAI
    /// owner before inclination and radiation calculations.
    pub fn publishLeafAreaIndex(self: State, canopy: *const canopy_module.State, plants: *PlantState, cell_area_m2: []const f64) !void {
        if (plants.cell_count != self.cell_count or plants.species_count != self.species_count or cell_area_m2.len != self.cell_count)
            return error.CanopyLayerDistributionDimensionMismatch;
        for (0..self.cell_count) |cell| {
            const area_m2 = cell_area_m2[cell];
            if (!std.math.isFinite(area_m2) or area_m2 <= 0) return error.InvalidCanopyCellArea;
            for (0..self.species_count) |species| {
                const plant = cell * self.species_count + species;
                plants.leaf_area_index_m2_m2[plant] = try self.plantLeafAreaM2(canopy, plant) / area_m2;
            }
        }
    }

    pub fn refresh(self: *State, canopy: *canopy_module.State, growth_stages: *const stages_module.State, controls: *const Controls, emerged_by_plant: []const bool, inclination_sine: []const f64, inclination_fraction_by_plant: []const f64, solar_angle_sine_by_cell: []const f64, minimum_area_m2: f64) !void {
        try self.ensureTopology(canopy);
        const plant_count = canopy.plant_branch_offsets.len - 1;
        const inclination_count = inclination_sine.len;
        if (growth_stages.plant_count != plant_count or controls.leaf_length_to_width_ratio.len != plant_count or controls.stem_angle_sine.len != plant_count or emerged_by_plant.len != plant_count or inclination_count != self.inclination_count or inclination_fraction_by_plant.len != plant_count * inclination_count or solar_angle_sine_by_cell.len != self.cell_count) return error.CanopyLayerDistributionDimensionMismatch;
        // HOUR1 adjusts boundaries from the prior GROSUB distribution first.
        for (0..self.cell_count) |cell| {
            const first = cell * self.layer_count;
            var canopy_height_m: f64 = 0;
            for (0..self.species_count) |species| {
                const plant = cell * self.species_count + species;
                const branches = try canopy.branchRange(plant);
                for (branches.first..branches.end) |branch| {
                    const nodes = try canopy.nodeRange(branch);
                    for (nodes.first..nodes.end) |node| canopy_height_m = @max(canopy_height_m, canopy.node_height_m[node] + canopy.node_sheath_height_m[node]);
                }
            }
            try structure_module.redistributeLayerBoundariesEqualArea(self.allocator, canopy_height_m, minimum_area_m2, self.cell_leaf_area_m2[first .. first + self.layer_count], self.cell_stalk_area_m2[first .. first + self.layer_count], self.cell_standing_dead_area_m2[first .. first + self.layer_count], try self.cellBoundaries(cell));
        }
        self.resetCurrentAreas();
        for (0..self.cell_count) |cell| for (0..self.species_count) |species| {
            const solar_angle_sine = solar_angle_sine_by_cell[cell];
            if (!std.math.isFinite(solar_angle_sine) or solar_angle_sine < 0 or solar_angle_sine > 1) return error.InvalidSolarAngle;
            const plant = cell * self.species_count + species;
            const boundaries = try self.cellBoundaries(cell);
            const dead_carbon_g = canopy.plant_standing_dead_carbon_g[plant];
            const dead_population_count = canopy.plant_standing_dead_population_count[plant];
            if (dead_carbon_g > minimum_area_m2 and dead_population_count > minimum_area_m2) {
                var live_canopy_height_m: f64 = 0;
                const all_branches = try canopy.branchRange(plant);
                for (all_branches.first..all_branches.end) |candidate_branch| {
                    const all_nodes = try canopy.nodeRange(candidate_branch);
                    for (all_nodes.first..all_nodes.end) |candidate_node| live_canopy_height_m = @max(live_canopy_height_m, canopy.node_height_m[candidate_node] + canopy.node_sheath_height_m[candidate_node]);
                }
                const dead_height_m = @max(0.01, @max(canopy.plant_standing_dead_height_m[plant], live_canopy_height_m));
                canopy.plant_standing_dead_height_m[plant] = dead_height_m;
                const radius_m = @sqrt(controls.stalk_volume_m3_per_g_c[plant] * (dead_carbon_g / dead_population_count) / (3.1416 * dead_height_m));
                const total_dead_area_m2 = 6.2832 * radius_m * dead_height_m * dead_population_count;
                const denominator_height_m = @min(dead_height_m, boundaries[boundaries.len - 1]);
                for (0..self.layer_count) |layer| {
                    const thickness_m = boundaries[layer + 1] - boundaries[layer];
                    if (boundaries[layer] >= dead_height_m or thickness_m <= 0) continue;
                    const fraction = @min(1.0, (dead_height_m - boundaries[layer]) / thickness_m);
                    const area_m2 = fraction * total_dead_area_m2 * thickness_m / denominator_height_m;
                    const plant_layer = plant * self.layer_count + layer;
                    self.plant_standing_dead_area_m2[plant_layer] = area_m2;
                    self.cell_standing_dead_area_m2[cell * self.layer_count + layer] += area_m2;
                    self.plant_standing_dead_projected_surface_m2[(plant_layer * inclination_count) + inclination_count - 1] = area_m2 / @as(f64, @floatFromInt(self.azimuth_count));
                }
            } else {
                canopy.plant_standing_dead_height_m[plant] = 0;
            }
            if (!emerged_by_plant[plant]) continue;
            const branch_range = try canopy.branchRange(plant);
            if (branch_range.first == branch_range.end) continue;
            const main_branch = (try growth_stages.mainLivingBranch(plant)) orelse continue;
            const main_nodes = try canopy.nodeRange(main_branch);
            for (branch_range.first..branch_range.end) |branch| {
                if (growth_stages.branches[branch].dead) continue;
                const nodes = try canopy.nodeRange(branch);
                var branch_base_height_m: f64 = 0;
                if (controls.biomass_turnover_type[plant] != 0 and controls.root_profile_type[plant] > 1 and branch != main_branch) {
                    const order = growth_stages.branches[branch].branch_order;
                    if (order < main_nodes.end - main_nodes.first) branch_base_height_m = canopy.node_height_m[main_nodes.first + order];
                }
                var branch_tip_height_m = branch_base_height_m;
                for (nodes.first..nodes.end) |node| {
                    const output = try self.nodeLayerRange(node);
                    const area = self.node_leaf_area_m2[output.first..output.end];
                    const carbon = self.node_leaf_carbon_g[output.first..output.end];
                    const nitrogen = self.node_leaf_nitrogen_g[output.first..output.end];
                    const phosphorus = self.node_leaf_phosphorus_g[output.first..output.end];
                    const stalk_height_m = branch_base_height_m + canopy.node_height_m[node];
                    branch_tip_height_m = @max(branch_tip_height_m, stalk_height_m);
                    const canopy_height_m = @max(0.0, boundaries[boundaries.len - 1] - 0.01);
                    _ = try canopy_module.allocateLeafAcrossCanopyLayers(canopy.node_leaf_area_m2[node], canopy.node_leaf_carbon_g[node], canopy.node_leaf_nitrogen_g[node], canopy.node_leaf_phosphorus_g[node], canopy.plant_population_per_m2[plant], controls.leaf_length_to_width_ratio[plant], stalk_height_m, canopy.node_sheath_height_m[node], canopy_height_m, boundaries, inclination_sine, inclination_fraction_by_plant[plant * inclination_count ..][0..inclination_count], .{ .area_m2 = area, .carbon_g = carbon, .nitrogen_g = nitrogen, .phosphorus_g = phosphorus });
                    const cell_first = cell * self.layer_count;
                    for (0..self.layer_count) |layer| {
                        self.cell_leaf_area_m2[cell_first + layer] += area[layer];
                        self.cell_leaf_carbon_g[cell_first + layer] += carbon[layer];
                        if (solar_angle_sine > 0) {
                            for (0..inclination_count) |inclination| {
                                self.node_leaf_projected_surface_m2[((node * self.layer_count + layer) * inclination_count) + inclination] = @max(0.0, inclination_fraction_by_plant[plant * inclination_count + inclination] * area[layer] / @as(f64, @floatFromInt(self.azimuth_count)));
                                self.plant_leaf_projected_surface_m2[((plant * self.layer_count + layer) * inclination_count) + inclination] += self.node_leaf_projected_surface_m2[((node * self.layer_count + layer) * inclination_count) + inclination];
                            }
                        }
                    }
                }
                const branch_layers = try self.branchLayerRange(branch);
                const stalk_area = self.branch_stalk_area_m2[branch_layers.first..branch_layers.end];
                const stalk = try canopy_module.allocateStalkAcrossCanopyLayers(canopy.branch_stalk_carbon_g[branch], canopy.branch_stalk_carbon_g[branch], canopy.plant_population_per_m2[plant], controls.stalk_volume_m3_per_g_c[plant], branch_base_height_m, branch_tip_height_m, controls.annual_growth_habit[plant], false, boundaries, stalk_area);
                canopy.branch_sapwood_carbon_g[branch] = stalk.sapwood_carbon_g;
                const cell_first = cell * self.layer_count;
                const inclination_width_radians = (0.5 * std.math.pi) / @as(f64, @floatFromInt(inclination_count));
                const stalk_inclination = if (branch == main_branch) inclination_count - 1 else @min(inclination_count - 1, @as(usize, @intFromFloat(std.math.asin(controls.stem_angle_sine[plant]) / inclination_width_radians)));
                for (0..self.layer_count) |layer| {
                    self.cell_stalk_area_m2[cell_first + layer] += stalk_area[layer];
                    if (solar_angle_sine > 0) {
                        const branch_surface = ((branch * self.layer_count + layer) * inclination_count) + stalk_inclination;
                        const plant_surface = ((plant * self.layer_count + layer) * inclination_count) + stalk_inclination;
                        self.branch_stalk_projected_surface_m2[branch_surface] = stalk_area[layer] / @as(f64, @floatFromInt(self.azimuth_count));
                        self.plant_stalk_projected_surface_m2[plant_surface] += self.branch_stalk_projected_surface_m2[branch_surface];
                    }
                }
            }
        };
    }

    /// Publishes the runtime layer × inclination × azimuth geometry consumed
    /// by canopy carboxylation. C/N/P follow the exact projected-area share.
    pub fn publishNodeSamples(self: State, canopy: *canopy_module.State) !void {
        if (canopy.cell_count != self.cell_count or canopy.species_count != self.species_count or canopy.node_leaf_area_m2.len != self.node_count) return error.CanopySampleTopologyMismatch;
        const angular_count = try std.math.mul(usize, self.inclination_count, self.azimuth_count);
        const expected_samples = try std.math.mul(usize, self.layer_count, angular_count);
        for (0..self.cell_count) |cell| {
            const boundaries = try self.cellBoundaries(cell);
            for (0..self.species_count) |species| {
                const plant = try canopy.plantIndex(cell, species);
                const branches = try canopy.branchRange(plant);
                for (branches.first..branches.end) |branch| {
                    const nodes = try canopy.nodeRange(branch);
                    for (nodes.first..nodes.end) |node| {
                        const samples = try canopy.sampleRange(node);
                        if (samples.end - samples.first != expected_samples) return error.CanopySampleTopologyMismatch;
                        for (0..self.layer_count) |layer| {
                            const node_layer = node * self.layer_count + layer;
                            const layer_area_m2 = self.node_leaf_area_m2[node_layer];
                            for (0..self.inclination_count) |inclination| {
                                const angular_area_m2 = self.node_leaf_projected_surface_m2[(node_layer * self.inclination_count) + inclination];
                                const mass_fraction = if (layer_area_m2 > 0) angular_area_m2 / layer_area_m2 else 0;
                                for (0..self.azimuth_count) |azimuth| {
                                    const sample = samples.first + layer * angular_count + inclination * self.azimuth_count + azimuth;
                                    canopy.sample_exposed_leaf_area_m2[sample] = angular_area_m2;
                                    canopy.sample_leaf_area_m2[sample] = angular_area_m2;
                                    canopy.sample_leaf_carbon_g[sample] = self.node_leaf_carbon_g[node_layer] * mass_fraction;
                                    canopy.sample_leaf_nitrogen_g[sample] = self.node_leaf_nitrogen_g[node_layer] * mass_fraction;
                                    canopy.sample_leaf_phosphorus_g[sample] = self.node_leaf_phosphorus_g[node_layer] * mass_fraction;
                                    canopy.sample_layer_lower_height_m[sample] = boundaries[layer];
                                    canopy.sample_layer_upper_height_m[sample] = boundaries[layer + 1];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};

fn fieldCount(comptime name: []const u8, cell_count: usize, species_count: usize, layer_count: usize, inclination_count: usize, node_count: usize, branch_count: usize) !usize {
    if (comptime std.mem.eql(u8, name, "boundary_height_m")) return try std.math.mul(usize, cell_count, try std.math.add(usize, layer_count, 1));
    if (comptime std.mem.eql(u8, name, "node_leaf_projected_surface_m2")) return try std.math.mul(usize, try std.math.mul(usize, node_count, layer_count), inclination_count);
    if (comptime std.mem.eql(u8, name, "branch_stalk_projected_surface_m2")) return try std.math.mul(usize, try std.math.mul(usize, branch_count, layer_count), inclination_count);
    if (comptime (std.mem.eql(u8, name, "plant_leaf_projected_surface_m2") or std.mem.eql(u8, name, "plant_stalk_projected_surface_m2"))) return try std.math.mul(usize, try std.math.mul(usize, try std.math.mul(usize, cell_count, species_count), layer_count), inclination_count);
    if (comptime std.mem.eql(u8, name, "plant_standing_dead_area_m2")) return try std.math.mul(usize, try std.math.mul(usize, cell_count, species_count), layer_count);
    if (comptime std.mem.eql(u8, name, "plant_standing_dead_projected_surface_m2")) return try std.math.mul(usize, try std.math.mul(usize, try std.math.mul(usize, cell_count, species_count), layer_count), inclination_count);
    if (comptime std.mem.startsWith(u8, name, "node_")) return try std.math.mul(usize, node_count, layer_count);
    if (comptime std.mem.startsWith(u8, name, "branch_")) return try std.math.mul(usize, branch_count, layer_count);
    return try std.math.mul(usize, cell_count, layer_count);
}

test "GROSUB dead branch clears node branch plant and cell layer aggregates" {
    var canopy = try canopy_module.State.init(std.testing.allocator, 1, 2, &.{ 1, 1 }, &.{ 1, 1 }, &.{ 0, 0 });
    defer canopy.deinit();
    var layers = try State.init(std.testing.allocator, 1, 2, 1, 1, 1, &canopy);
    defer layers.deinit();
    layers.node_leaf_area_m2[0] = 2;
    layers.node_leaf_area_m2[1] = 3;
    layers.node_leaf_carbon_g[0] = 4;
    layers.node_leaf_carbon_g[1] = 5;
    layers.node_leaf_nitrogen_g[0] = 0.4;
    layers.node_leaf_phosphorus_g[0] = 0.04;
    layers.node_leaf_projected_surface_m2[0] = 2;
    layers.node_leaf_projected_surface_m2[1] = 3;
    layers.branch_stalk_area_m2[0] = 1;
    layers.branch_stalk_area_m2[1] = 2;
    layers.branch_stalk_projected_surface_m2[0] = 1;
    layers.branch_stalk_projected_surface_m2[1] = 2;
    layers.plant_leaf_projected_surface_m2[0] = 2;
    layers.plant_leaf_projected_surface_m2[1] = 3;
    layers.plant_stalk_projected_surface_m2[0] = 1;
    layers.plant_stalk_projected_surface_m2[1] = 2;
    layers.cell_leaf_area_m2[0] = 5;
    layers.cell_leaf_carbon_g[0] = 9;
    layers.cell_stalk_area_m2[0] = 3;

    try layers.clearDeadBranch(&canopy, 0, 0);
    try std.testing.expectEqual(@as(f64, 0), layers.node_leaf_area_m2[0]);
    try std.testing.expectEqual(@as(f64, 3), layers.node_leaf_area_m2[1]);
    try std.testing.expectEqual(@as(f64, 3), layers.cell_leaf_area_m2[0]);
    try std.testing.expectEqual(@as(f64, 5), layers.cell_leaf_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 2), layers.cell_stalk_area_m2[0]);
    try std.testing.expectEqual(@as(f64, 0), layers.plant_leaf_projected_surface_m2[0]);
    try std.testing.expectEqual(@as(f64, 3), layers.plant_leaf_projected_surface_m2[1]);

    const config = try @import("config.zig").SimulationConfig.init(
        .{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 2 },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 },
    );
    var plants = try PlantState.init(std.testing.allocator, config);
    defer plants.deinit();
    try layers.publishLeafAreaIndex(&canopy, &plants, &.{2});
    try std.testing.expectEqual(@as(f64, 0), plants.leaf_area_index_m2_m2[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), plants.leaf_area_index_m2_m2[1], 1e-12);
}

test "node samples conserve runtime layer and angular C N P geometry" {
    var canopy = try canopy_module.State.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{12});
    defer canopy.deinit();
    var layers = try State.init(std.testing.allocator, 1, 1, 2, 2, 3, &canopy);
    defer layers.deinit();
    layers.boundary_height_m[0] = 0;
    layers.boundary_height_m[1] = 1;
    layers.boundary_height_m[2] = 2;
    layers.node_leaf_area_m2[0] = 3;
    layers.node_leaf_area_m2[1] = 6;
    layers.node_leaf_carbon_g[0] = 9;
    layers.node_leaf_carbon_g[1] = 12;
    layers.node_leaf_nitrogen_g[0] = 0.9;
    layers.node_leaf_nitrogen_g[1] = 1.2;
    layers.node_leaf_phosphorus_g[0] = 0.09;
    layers.node_leaf_phosphorus_g[1] = 0.12;
    // Values are per azimuth, matching refresh().
    layers.node_leaf_projected_surface_m2[0] = 0.4;
    layers.node_leaf_projected_surface_m2[1] = 0.6;
    layers.node_leaf_projected_surface_m2[2] = 0.5;
    layers.node_leaf_projected_surface_m2[3] = 1.5;
    try layers.publishNodeSamples(&canopy);
    var area_m2: f64 = 0;
    var carbon_g: f64 = 0;
    var nitrogen_g: f64 = 0;
    var phosphorus_g: f64 = 0;
    for (canopy.sample_leaf_area_m2, canopy.sample_leaf_carbon_g, canopy.sample_leaf_nitrogen_g, canopy.sample_leaf_phosphorus_g) |area, carbon, nitrogen, phosphorus| {
        area_m2 += area;
        carbon_g += carbon;
        nitrogen_g += nitrogen;
        phosphorus_g += phosphorus;
    }
    try std.testing.expectApproxEqAbs(@as(f64, 9), area_m2, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 21), carbon_g, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2.1), nitrogen_g, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.21), phosphorus_g, 1e-14);
}

test "runtime canopy layers resynchronize after dynamic branch insertion" {
    const allocator = std.testing.allocator;
    var canopy = try canopy_module.State.init(allocator, 1, 2, &.{ 1, 1 }, &.{ 1, 1 }, &.{ 1, 1 });
    defer canopy.deinit();
    var layers = try State.init(allocator, 1, 2, 17, 4, 4, &canopy);
    defer layers.deinit();
    layers.boundary_height_m[17] = 2.01;
    _ = try canopy.appendBranch(0, &.{ 1, 1, 1 });
    try layers.ensureTopology(&canopy);
    try std.testing.expectEqual(@as(usize, 5), layers.node_count);
    try std.testing.expectEqual(@as(usize, 3), layers.branch_count);
    try std.testing.expectEqual(@as(usize, 5 * 17), layers.node_leaf_area_m2.len);
    try std.testing.expectEqual(@as(f64, 2.01), layers.boundary_height_m[17]);
}

test "GROSUB runtime layer surfaces conserve live and standing dead geometry" {
    const allocator = std.testing.allocator;
    var canopy = try canopy_module.State.init(allocator, 1, 1, &.{1}, &.{1}, &.{1});
    defer canopy.deinit();
    canopy.plant_population_per_m2[0] = 3;
    canopy.plant_population_count[0] = 3;
    canopy.plant_standing_dead_population_count[0] = 3;
    canopy.plant_standing_dead_carbon_g[0] = 9;
    canopy.node_leaf_area_m2[0] = 2;
    canopy.node_leaf_carbon_g[0] = 4;
    canopy.node_leaf_nitrogen_g[0] = 0.2;
    canopy.node_leaf_phosphorus_g[0] = 0.04;
    canopy.node_height_m[0] = 0.5;
    canopy.node_sheath_height_m[0] = 0.2;
    canopy.branch_stalk_carbon_g[0] = 6;
    var stages = try stages_module.State.init(allocator, &.{1});
    defer stages.deinit();
    var controls = try Controls.init(allocator, 1);
    defer controls.deinit();
    try controls.setPlant(0, 2, 4.0e-6, 0, 1, true, 0.5);
    var layers = try State.init(allocator, 1, 1, 13, 4, 4, &canopy);
    defer layers.deinit();
    try layers.refresh(&canopy, &stages, &controls, &.{true}, &.{ 0.1951, 0.5556, 0.8315, 0.9808 }, &.{ 0.1, 0.2, 0.3, 0.4 }, &.{0.5}, 1.0e-12);

    var leaf_area_m2: f64 = 0;
    var leaf_projected_m2: f64 = 0;
    var stalk_area_m2: f64 = 0;
    var stalk_projected_m2: f64 = 0;
    var dead_area_m2: f64 = 0;
    var dead_projected_m2: f64 = 0;
    for (layers.node_leaf_area_m2) |value| leaf_area_m2 += value;
    for (layers.node_leaf_projected_surface_m2) |value| leaf_projected_m2 += value;
    for (layers.branch_stalk_area_m2) |value| stalk_area_m2 += value;
    for (layers.branch_stalk_projected_surface_m2) |value| stalk_projected_m2 += value;
    for (layers.plant_standing_dead_area_m2) |value| dead_area_m2 += value;
    for (layers.plant_standing_dead_projected_surface_m2) |value| dead_projected_m2 += value;
    try std.testing.expectApproxEqAbs(@as(f64, 2), leaf_area_m2, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.25 * leaf_area_m2, leaf_projected_m2, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.25 * stalk_area_m2, stalk_projected_m2, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.25 * dead_area_m2, dead_projected_m2, 1.0e-12);
    var plant_leaf_projected_m2: f64 = 0;
    var plant_stalk_projected_m2: f64 = 0;
    for (layers.plant_leaf_projected_surface_m2) |value| plant_leaf_projected_m2 += value;
    for (layers.plant_stalk_projected_surface_m2) |value| plant_stalk_projected_m2 += value;
    try std.testing.expectApproxEqAbs(leaf_area_m2, try layers.plantLeafAreaM2(&canopy, 0), 1.0e-12);
    try std.testing.expectApproxEqAbs(leaf_projected_m2, plant_leaf_projected_m2, 1.0e-12);
    try std.testing.expectApproxEqAbs(stalk_projected_m2, plant_stalk_projected_m2, 1.0e-12);
    const expected_dead_radius_m = @sqrt(4.0e-6 * (9.0 / 3.0) / (3.1416 * 0.7));
    try std.testing.expectApproxEqAbs(6.2832 * expected_dead_radius_m * 0.7 * 3.0, dead_area_m2, 1.0e-12);
}
