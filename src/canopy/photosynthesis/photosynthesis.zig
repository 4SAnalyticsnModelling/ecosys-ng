const std = @import("std");

const c4_mesophyll_bundle_exchange = @import("c4_mesophyll_bundle_exchange.zig");

const branch_organ_growth_publication = @import("../../plant/growth/branch_organ_growth_publication.zig");

const leaf_node_growth_publication = @import("../leaf/node_growth_publication.zig");

const shoot_recycling_fraction = @import("../../plant/growth/shoot_recycling_fraction.zig");

const reserve_maintenance_respiration = @import("../../plant/growth/reserve_maintenance_respiration.zig");

const shoot_total_senescence_setup = @import("../../plant/growth/shoot_total_senescence_setup.zig");

const node_senescence_remobilization_request = @import("../../plant/growth/node_senescence_remobilization_request.zig");

const c4_leaf_nonstructural_carbon_senescence = @import("../leaf/c4_nonstructural_carbon_senescence.zig");

const node_senescence_cascade_progress = @import("../../plant/growth/node_senescence_cascade_progress.zig");

const perennial_stalk_senescence_setup = @import("../../plant/growth/perennial_stalk_senescence_setup.zig");

const internode_senescence_publication = @import("../sheath/internode_senescence_publication.zig");

const residual_stalk_senescence_request = @import("../../plant/growth/residual_stalk_senescence_request.zig");

const residual_stalk_senescence_publication = @import("../../plant/growth/residual_stalk_senescence_publication.zig");

/// Compact heap-owned canopy topology. Prefix offsets permit different branch,
/// node, layer/orientation sample counts for every runtime plant population.
pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    species_count: usize,
    plant_branch_offsets: []usize,
    branch_node_offsets: []usize,
    node_sample_offsets: []usize,
    branch_c3_feedback_fraction: []f64,
    branch_c4_feedback_fraction: []f64,
    branch_carboxylation_umol_per_s: []f64,
    branch_fixed_carbon_g_c_per_h: []f64,
    branch_shoot_carbohydrate_g_c_per_h: []f64,
    branch_mobile_carbon_g: []f64,
    branch_mobile_nitrogen_g: []f64,
    branch_mobile_phosphorus_g: []f64,
    branch_symbiont_mobile_carbon_g: []f64,
    branch_symbiont_mobile_nitrogen_g: []f64,
    branch_symbiont_mobile_phosphorus_g: []f64,
    branch_symbiont_structural_carbon_g: []f64,
    branch_symbiont_structural_nitrogen_g: []f64,
    branch_symbiont_structural_phosphorus_g: []f64,
    branch_symbiotic_fixed_nitrogen_g_n_per_h: []f64,
    branch_symbiotic_respiration_g_c_per_h: []f64,
    branch_canopy_ammonia_exchange_g_n_per_h: []f64,
    branch_salt_content_by_species_mol: []f64,
    branch_combustion_salt_loss_by_species_mol_per_h: []f64,
    branch_mobile_carbon_concentration_g_per_g: []f64,
    branch_mobile_nitrogen_concentration_g_per_g: []f64,
    branch_mobile_phosphorus_concentration_g_per_g: []f64,
    branch_leaf_carbon_g: []f64,
    branch_leaf_nitrogen_g: []f64,
    branch_leaf_phosphorus_g: []f64,
    branch_sheath_carbon_g: []f64,
    branch_sheath_nitrogen_g: []f64,
    branch_sheath_phosphorus_g: []f64,
    branch_stalk_carbon_g: []f64,
    branch_stalk_nitrogen_g: []f64,
    branch_stalk_phosphorus_g: []f64,
    branch_sapwood_carbon_g: []f64,
    branch_reserve_carbon_g: []f64,
    branch_reserve_nitrogen_g: []f64,
    branch_reserve_phosphorus_g: []f64,
    branch_husk_carbon_g: []f64,
    branch_husk_nitrogen_g: []f64,
    branch_husk_phosphorus_g: []f64,
    branch_ear_carbon_g: []f64,
    branch_ear_nitrogen_g: []f64,
    branch_ear_phosphorus_g: []f64,
    branch_grain_carbon_g: []f64,
    branch_grain_nitrogen_g: []f64,
    branch_grain_phosphorus_g: []f64,
    branch_potential_seed_site_count: []f64,
    branch_seed_count: []f64,
    branch_individual_seed_carbon_g: []f64,
    branch_senescing_stalk_carbon_g: []f64,
    branch_senescing_stalk_nitrogen_g: []f64,
    branch_senescing_stalk_phosphorus_g: []f64,
    branch_leaf_area_m2: []f64,
    node_leaf_area_m2: []f64,
    node_height_m: []f64,
    node_internode_length_m: []f64,
    node_sheath_height_m: []f64,
    node_leaf_carbon_g: []f64,
    node_leaf_protein_g: []f64,
    node_leaf_nitrogen_g: []f64,
    node_leaf_phosphorus_g: []f64,
    node_sheath_carbon_g: []f64,
    node_sheath_protein_g: []f64,
    node_sheath_nitrogen_g: []f64,
    node_sheath_phosphorus_g: []f64,
    node_internode_carbon_g: []f64,
    node_internode_nitrogen_g: []f64,
    node_internode_phosphorus_g: []f64,
    node_c3_nonstructural_carbon_g: []f64,
    node_c4_mesophyll_nonstructural_carbon_g: []f64,
    node_bundle_sheath_co2_carbon_g: []f64,
    node_bundle_sheath_bicarbonate_carbon_g: []f64,
    node_co2_unlimited_carboxylation_umol_per_m2_s: []f64,
    node_co2_limited_carboxylation_umol_per_m2_s: []f64,
    node_co2_compensation_umol_per_l: []f64,
    node_co2_solubility_umol_per_l_per_umol_per_mol: []f64,
    node_carboxylation_half_saturation_umol_per_l: []f64,
    node_light_saturated_electron_transport_umol_per_m2_s: []f64,
    node_carboxylation_umol_co2_per_umol_electron: []f64,
    node_c4_feedback_fraction: []f64,
    node_pep_carboxylase_surface_density_g_per_m2: []f64,
    node_mesophyll_chlorophyll_surface_density_g_per_m2: []f64,
    node_bundle_sheath_co2_limited_carboxylation_umol_per_m2_s: []f64,
    node_bundle_sheath_light_saturated_electron_transport_umol_per_m2_s: []f64,
    node_bundle_sheath_carboxylation_umol_co2_per_umol_electron: []f64,
    sample_exposed_leaf_area_m2: []f64,
    sample_leaf_area_m2: []f64,
    sample_leaf_carbon_g: []f64,
    sample_leaf_nitrogen_g: []f64,
    sample_leaf_phosphorus_g: []f64,
    sample_stalk_area_m2: []f64,
    sample_layer_lower_height_m: []f64,
    sample_layer_upper_height_m: []f64,
    sample_direct_par_umol_per_m2_s: []f64,
    sample_diffuse_par_umol_per_m2_s: []f64,
    sample_direct_transmission_fraction: []f64,
    sample_diffuse_transmission_fraction: []f64,
    sample_intercellular_co2_umol_per_mol: []f64,
    sample_carboxylation_umol_per_s: []f64,
    sample_bundle_sheath_carboxylation_umol_per_s: []f64,
    plant_carboxylation_umol_per_s: []f64,
    plant_gross_primary_productivity_g_c_per_h: []f64,
    plant_minimum_water_vapor_resistance_h_per_m: []f64,
    plant_population_per_m2: []f64,
    plant_population_count: []f64,
    plant_population_change_count: []f64,
    plant_stem_diameter_m: []f64,
    plant_standing_dead_population_count: []f64,
    plant_cuticular_water_vapor_resistance_h_per_m: []f64,
    plant_cuticular_co2_resistance_s_per_m: []f64,
    plant_intercellular_oxygen_umol_per_mol: []f64,
    plant_thermal_adaptation_offset_c: []f64,
    plant_chilling_stress_h: []f64,
    plant_heat_stress_h: []f64,
    plant_uptake_growth_temperature_response: []f64,
    plant_minimum_daily_canopy_water_potential_megapascal: []f64,
    plant_leafout_threshold_c: []f64,
    plant_leafoff_threshold_c: []f64,
    plant_seed_set_high_temperature_c: []f64,
    plant_seed_set_loss_fraction_per_c_h: []f64,
    plant_seed_storage_carbon_g: []f64,
    plant_seed_storage_nitrogen_g: []f64,
    plant_seed_storage_phosphorus_g: []f64,
    plant_standing_dead_carbon_g: []f64,
    plant_standing_dead_nitrogen_g: []f64,
    plant_standing_dead_phosphorus_g: []f64,
    plant_charcoal_carbon_g: []f64,
    plant_charcoal_nitrogen_g: []f64,
    plant_charcoal_phosphorus_g: []f64,
    plant_standing_dead_height_m: []f64,
    plant_standing_dead_carbon_by_kinetic_g: []f64,
    plant_standing_dead_nitrogen_by_kinetic_g: []f64,
    plant_standing_dead_phosphorus_by_kinetic_g: []f64,
    plant_canopy_aerodynamic_temperature_k: []f64,
    plant_canopy_aerodynamic_vapor_pressure_kpa: []f64,
    plant_standing_dead_aerodynamic_temperature_k: []f64,
    plant_standing_dead_aerodynamic_vapor_pressure_kpa: []f64,
    plant_standing_dead_surface_temperature_k: []f64,
    plant_phenology_temperature_k: []f64,
    plant_canopy_osmotic_potential_megapascal: []f64,
    plant_canopy_turgor_potential_megapascal: []f64,
    plant_stored_energy_megajoules: []f64,
    plant_transpiration_m3_per_h: []f64,
    plant_hypocotyledon_height_m: []f64,
    plant_mobile_carbon_g: []f64,
    plant_mobile_nitrogen_g: []f64,
    plant_mobile_phosphorus_g: []f64,
    plant_symbiont_mobile_carbon_g: []f64,
    plant_symbiont_mobile_nitrogen_g: []f64,
    plant_symbiont_mobile_phosphorus_g: []f64,
    plant_salt_content_mol: []f64,
    plant_mobile_carbon_concentration_g_per_g: []f64,
    plant_mobile_nitrogen_concentration_g_per_g: []f64,
    plant_mobile_phosphorus_concentration_g_per_g: []f64,
    plant_symbiont_mobile_carbon_concentration_g_per_g: []f64,
    plant_salt_concentration_mol_per_g_c: []f64,
    plant_nitrogen_phosphorus_fixation_constraint_fraction: []f64,
    plant_leaf_sheath_partition_fraction: []f64,
    plant_total_shoot_carbon_g: []f64,
    plant_previous_total_shoot_carbon_g: []f64,
    plant_shoot_growth_g_c_per_step: []f64,
    plant_combustion_carbon_loss_g_c_per_h: []f64,
    plant_combustion_nitrogen_loss_g_n_per_h: []f64,
    plant_combustion_phosphorus_loss_g_p_per_h: []f64,
    plant_live_combustion_g_c_per_h: []f64,
    plant_standing_dead_combustion_g_c_per_h: []f64,
    plant_fire_carbon_dioxide_emission_g_c_per_h: []f64,
    plant_fire_methane_emission_g_c_per_h: []f64,
    plant_fire_oxygen_consumption_g_o_per_h: []f64,
    plant_fire_charcoal_production_g_c_per_h: []f64,
    plant_fire_heat_release_megajoules_per_h: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, species_count: usize, branch_count_by_plant: []const usize, node_count_by_branch: []const usize, sample_count_by_node: []const usize) !State {
        if (cell_count == 0 or species_count == 0) return error.InvalidCanopyPhotosynthesisDimensions;
        const plant_count = try std.math.mul(usize, cell_count, species_count);
        const plant_kinetic_count = try std.math.mul(usize, plant_count, 4);
        if (branch_count_by_plant.len != plant_count) return error.CanopyBranchCountDimensionMismatch;
        const branch_count = try checkedSum(branch_count_by_plant);
        if (node_count_by_branch.len != branch_count) return error.CanopyNodeCountDimensionMismatch;
        const node_count = try checkedSum(node_count_by_branch);
        if (sample_count_by_node.len != node_count) return error.CanopySampleCountDimensionMismatch;
        const sample_count = try checkedSum(sample_count_by_node);

        var result: State = undefined;
        result.allocator = allocator;
        result.cell_count = cell_count;
        result.species_count = species_count;
        result.plant_branch_offsets = try makeOffsets(allocator, branch_count_by_plant);
        errdefer allocator.free(result.plant_branch_offsets);
        result.branch_node_offsets = try makeOffsets(allocator, node_count_by_branch);
        errdefer allocator.free(result.branch_node_offsets);
        result.node_sample_offsets = try makeOffsets(allocator, sample_count_by_node);
        errdefer allocator.free(result.node_sample_offsets);
        var allocated_f64_fields: usize = 0;
        errdefer freeAllocatedF64Fields(&result, allocated_f64_fields);
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            const count = domainCount(field.name, plant_count, plant_kinetic_count, branch_count, node_count, sample_count);
            @field(result, field.name) = try allocateZeroed(allocator, count);
            allocated_f64_fields += 1;
        };
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.allocator.free(self.node_sample_offsets);
        self.allocator.free(self.branch_node_offsets);
        self.allocator.free(self.plant_branch_offsets);
        self.* = undefined;
    }

    pub fn plantIndex(self: State, cell: usize, species: usize) !usize {
        if (cell >= self.cell_count or species >= self.species_count) return error.CanopyPlantIndexOutOfBounds;
        return cell * self.species_count + species;
    }

    pub fn branchRange(self: State, plant: usize) !Range {
        if (plant + 1 >= self.plant_branch_offsets.len) return error.CanopyPlantIndexOutOfBounds;
        return .{ .first = self.plant_branch_offsets[plant], .end = self.plant_branch_offsets[plant + 1] };
    }

    pub fn clone(self: State) !State {
        const plant_count = self.plant_branch_offsets.len - 1;
        const branch_count = self.branch_node_offsets.len - 1;
        const node_count = self.node_sample_offsets.len - 1;
        const branch_counts = try self.allocator.alloc(usize, plant_count);
        defer self.allocator.free(branch_counts);
        const node_counts = try self.allocator.alloc(usize, branch_count);
        defer self.allocator.free(node_counts);
        const sample_counts = try self.allocator.alloc(usize, node_count);
        defer self.allocator.free(sample_counts);
        for (branch_counts, 0..) |*count, plant| count.* = self.plant_branch_offsets[plant + 1] - self.plant_branch_offsets[plant];
        for (node_counts, 0..) |*count, branch| count.* = self.branch_node_offsets[branch + 1] - self.branch_node_offsets[branch];
        for (sample_counts, 0..) |*count, node| count.* = self.node_sample_offsets[node + 1] - self.node_sample_offsets[node];
        const result = try State.init(self.allocator, self.cell_count, self.species_count, branch_counts, node_counts, sample_counts);
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) @memcpy(@field(result, field.name), @field(self, field.name));
        return result;
    }

    pub fn nodeRange(self: State, branch: usize) !Range {
        if (branch + 1 >= self.branch_node_offsets.len) return error.CanopyBranchIndexOutOfBounds;
        return .{ .first = self.branch_node_offsets[branch], .end = self.branch_node_offsets[branch + 1] };
    }

    pub fn sampleRange(self: State, node: usize) !Range {
        if (node + 1 >= self.node_sample_offsets.len) return error.CanopyNodeIndexOutOfBounds;
        return .{ .first = self.node_sample_offsets[node], .end = self.node_sample_offsets[node + 1] };
    }

    /// Clears all persistent and diagnostic canopy coordinates owned by one
    /// runtime plant while retaining its allocated compact topology. STARTQ
    /// initialization can then reconstruct the new crop in-place without
    /// carrying harvested organ, mobile-pool, salt, or photosynthetic history.
    pub fn clearPlantForReconstruction(self: *State, plant: usize) !void {
        @setEvalBranchQuota(40_000);
        const plant_count = try std.math.mul(usize, self.cell_count, self.species_count);
        if (plant >= plant_count) return error.CanopyPlantIndexOutOfBounds;
        const branches = try self.branchRange(plant);
        const node_first = self.branch_node_offsets[branches.first];
        const node_end = self.branch_node_offsets[branches.end];
        const sample_first = self.node_sample_offsets[node_first];
        const sample_end = self.node_sample_offsets[node_end];
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            const values = @field(self, field.name);
            if (comptime std.mem.indexOf(u8, field.name, "_by_kinetic_") != null) {
                @memset(values[plant * 4 .. (plant + 1) * 4], 0);
            } else if (comptime std.mem.indexOf(u8, field.name, "_by_species_") != null and std.mem.startsWith(u8, field.name, "branch_")) {
                @memset(values[branches.first * 8 .. branches.end * 8], 0);
            } else if (comptime std.mem.startsWith(u8, field.name, "plant_")) {
                values[plant] = 0;
            } else if (comptime std.mem.startsWith(u8, field.name, "branch_")) {
                @memset(values[branches.first..branches.end], 0);
            } else if (comptime std.mem.startsWith(u8, field.name, "node_")) {
                @memset(values[node_first..node_end], 0);
            } else if (comptime std.mem.startsWith(u8, field.name, "sample_")) {
                @memset(values[sample_first..sample_end], 0);
            } else unreachable;
        };
    }

    /// Atomically restores one previously grown plant to STARTQ's one branch,
    /// one node topology. Other plants retain every value and relative order.
    pub fn compactPlantToInitialTopology(self: *State, plant: usize) !void {
        @setEvalBranchQuota(50_000);
        const plant_count = self.plant_branch_offsets.len - 1;
        if (plant >= plant_count) return error.CanopyPlantIndexOutOfBounds;
        const branches = try self.branchRange(plant);
        if (branches.first == branches.end) return error.CanopyInitialTopologyMissingBranch;
        if (branches.end - branches.first == 1 and self.branch_node_offsets[branches.first + 1] - self.branch_node_offsets[branches.first] == 1) return;
        const retained_node = self.branch_node_offsets[branches.first];
        if (retained_node == self.branch_node_offsets[branches.first + 1]) return error.CanopyInitialTopologyMissingNode;
        const branch_remove_first = branches.first + 1;
        const branch_remove_end = branches.end;
        const node_remove_first = retained_node + 1;
        const node_remove_end = self.branch_node_offsets[branches.end];
        const sample_remove_first = self.node_sample_offsets[node_remove_first];
        const sample_remove_end = self.node_sample_offsets[node_remove_end];

        const branch_counts = try self.allocator.alloc(usize, plant_count);
        defer self.allocator.free(branch_counts);
        for (branch_counts, 0..) |*count, index| count.* = if (index == plant) 1 else self.plant_branch_offsets[index + 1] - self.plant_branch_offsets[index];
        const new_branch_count = (self.branch_node_offsets.len - 1) - (branch_remove_end - branch_remove_first);
        const node_counts = try self.allocator.alloc(usize, new_branch_count);
        defer self.allocator.free(node_counts);
        var destination_branch: usize = 0;
        for (0..self.branch_node_offsets.len - 1) |source_branch| {
            if (source_branch >= branch_remove_first and source_branch < branch_remove_end) continue;
            node_counts[destination_branch] = if (source_branch == branches.first) 1 else self.branch_node_offsets[source_branch + 1] - self.branch_node_offsets[source_branch];
            destination_branch += 1;
        }
        const new_node_count = (self.node_sample_offsets.len - 1) - (node_remove_end - node_remove_first);
        const sample_counts = try self.allocator.alloc(usize, new_node_count);
        defer self.allocator.free(sample_counts);
        var destination_node: usize = 0;
        for (0..self.node_sample_offsets.len - 1) |source_node| {
            if (source_node >= node_remove_first and source_node < node_remove_end) continue;
            sample_counts[destination_node] = self.node_sample_offsets[source_node + 1] - self.node_sample_offsets[source_node];
            destination_node += 1;
        }
        var replacement = try State.init(self.allocator, self.cell_count, self.species_count, branch_counts, node_counts, sample_counts);
        errdefer replacement.deinit();
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            if (comptime std.mem.indexOf(u8, field.name, "_by_species_") != null and std.mem.startsWith(u8, field.name, "branch_"))
                copyRemovingRange(@field(replacement, field.name), @field(self, field.name), branch_remove_first * 8, branch_remove_end * 8)
            else if (comptime std.mem.startsWith(u8, field.name, "branch_"))
                copyRemovingRange(@field(replacement, field.name), @field(self, field.name), branch_remove_first, branch_remove_end)
            else if (comptime std.mem.startsWith(u8, field.name, "node_"))
                copyRemovingRange(@field(replacement, field.name), @field(self, field.name), node_remove_first, node_remove_end)
            else if (comptime std.mem.startsWith(u8, field.name, "sample_"))
                copyRemovingRange(@field(replacement, field.name), @field(self, field.name), sample_remove_first, sample_remove_end)
            else
                @memcpy(@field(replacement, field.name), @field(self, field.name));
        };
        var previous = self.*;
        self.* = replacement;
        previous.deinit();
    }

    /// Inserts a branch after the selected plant's existing branches. Growth is
    /// infrequent relative to hourly kernels, so an atomic compact rebuild keeps
    /// iteration contiguous while removing any compile-time branch/node ceiling.
    pub fn appendBranch(self: *State, plant: usize, sample_count_by_new_node: []const usize) !usize {
        @setEvalBranchQuota(50_000);
        if (plant + 1 >= self.plant_branch_offsets.len) return error.CanopyPlantIndexOutOfBounds;
        const plant_count = self.plant_branch_offsets.len - 1;
        const old_branch_count = self.branch_node_offsets.len - 1;
        const old_node_count = self.node_sample_offsets.len - 1;
        const inserted_branch = self.plant_branch_offsets[plant + 1];
        const inserted_node = self.branch_node_offsets[inserted_branch];
        const inserted_sample = self.node_sample_offsets[inserted_node];
        const added_node_count = sample_count_by_new_node.len;
        const added_sample_count = try checkedSum(sample_count_by_new_node);

        const branch_counts = try self.allocator.alloc(usize, plant_count);
        defer self.allocator.free(branch_counts);
        for (branch_counts, 0..) |*count, index| count.* = self.plant_branch_offsets[index + 1] - self.plant_branch_offsets[index] + @intFromBool(index == plant);
        const node_counts = try self.allocator.alloc(usize, try std.math.add(usize, old_branch_count, 1));
        defer self.allocator.free(node_counts);
        for (0..inserted_branch) |branch| node_counts[branch] = self.branch_node_offsets[branch + 1] - self.branch_node_offsets[branch];
        node_counts[inserted_branch] = added_node_count;
        for (inserted_branch..old_branch_count) |branch| node_counts[branch + 1] = self.branch_node_offsets[branch + 1] - self.branch_node_offsets[branch];
        const sample_counts = try self.allocator.alloc(usize, try std.math.add(usize, old_node_count, added_node_count));
        defer self.allocator.free(sample_counts);
        for (0..inserted_node) |node| sample_counts[node] = self.node_sample_offsets[node + 1] - self.node_sample_offsets[node];
        @memcpy(sample_counts[inserted_node .. inserted_node + added_node_count], sample_count_by_new_node);
        for (inserted_node..old_node_count) |node| sample_counts[node + added_node_count] = self.node_sample_offsets[node + 1] - self.node_sample_offsets[node];

        var replacement = try State.init(self.allocator, self.cell_count, self.species_count, branch_counts, node_counts, sample_counts);
        errdefer replacement.deinit();
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            if (comptime std.mem.indexOf(u8, field.name, "_by_species_") != null and std.mem.startsWith(u8, field.name, "branch_"))
                copyAroundInsertion(@field(replacement, field.name), @field(self, field.name), inserted_branch * 8, 8)
            else if (comptime std.mem.startsWith(u8, field.name, "branch_"))
                copyAroundInsertion(@field(replacement, field.name), @field(self, field.name), inserted_branch, 1)
            else if (comptime std.mem.startsWith(u8, field.name, "node_"))
                copyAroundInsertion(@field(replacement, field.name), @field(self, field.name), inserted_node, added_node_count)
            else if (comptime std.mem.startsWith(u8, field.name, "sample_"))
                copyAroundInsertion(@field(replacement, field.name), @field(self, field.name), inserted_sample, added_sample_count)
            else
                @memcpy(@field(replacement, field.name), @field(self, field.name));
        };
        var previous = self.*;
        self.* = replacement;
        previous.deinit();
        return inserted_branch;
    }

    /// Appends one node to a runtime branch through an atomic compact rebuild.
    /// Existing branch, node, and sample values retain their indices except for
    /// nodes and samples following the insertion point, which shift together.
    pub fn appendNode(self: *State, branch: usize, sample_count: usize) !usize {
        @setEvalBranchQuota(50_000);
        if (branch + 1 >= self.branch_node_offsets.len) return error.CanopyBranchIndexOutOfBounds;
        if (sample_count == 0) return error.InvalidCanopySampleCount;
        const plant_count = self.plant_branch_offsets.len - 1;
        const branch_count = self.branch_node_offsets.len - 1;
        const old_node_count = self.node_sample_offsets.len - 1;
        const inserted_node = self.branch_node_offsets[branch + 1];
        const inserted_sample = self.node_sample_offsets[inserted_node];

        const branch_counts = try self.allocator.alloc(usize, plant_count);
        defer self.allocator.free(branch_counts);
        const node_counts = try self.allocator.alloc(usize, branch_count);
        defer self.allocator.free(node_counts);
        const sample_counts = try self.allocator.alloc(usize, try std.math.add(usize, old_node_count, 1));
        defer self.allocator.free(sample_counts);
        for (branch_counts, 0..) |*count, plant| count.* = self.plant_branch_offsets[plant + 1] - self.plant_branch_offsets[plant];
        for (node_counts, 0..) |*count, index| count.* = self.branch_node_offsets[index + 1] - self.branch_node_offsets[index] + @intFromBool(index == branch);
        for (0..inserted_node) |node| sample_counts[node] = self.node_sample_offsets[node + 1] - self.node_sample_offsets[node];
        sample_counts[inserted_node] = sample_count;
        for (inserted_node..old_node_count) |node| sample_counts[node + 1] = self.node_sample_offsets[node + 1] - self.node_sample_offsets[node];

        var replacement = try State.init(self.allocator, self.cell_count, self.species_count, branch_counts, node_counts, sample_counts);
        errdefer replacement.deinit();
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            if (comptime std.mem.startsWith(u8, field.name, "node_"))
                copyAroundInsertion(@field(replacement, field.name), @field(self, field.name), inserted_node, 1)
            else if (comptime std.mem.startsWith(u8, field.name, "sample_"))
                copyAroundInsertion(@field(replacement, field.name), @field(self, field.name), inserted_sample, sample_count)
            else
                @memcpy(@field(replacement, field.name), @field(self, field.name));
        };
        var previous = self.*;
        self.* = replacement;
        previous.deinit();
        return inserted_node;
    }

    pub fn validateFinite(self: State) !void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) for (@field(self, field.name), 0..) |value, index| if (!std.math.isFinite(value)) {
            std.log.err("non-finite canopy photosynthesis state: field={s} index={d} value={e}", .{ field.name, index, value });
            return error.NonFiniteCanopyPhotosynthesisState;
        };
    }
};

pub const PersistentReseedInventories = struct {
    seed_storage: ElementalMass,
    standing_dead: ElementalMass,
    standing_dead_height_m: f64,
    standing_dead_carbon_by_kinetic_g: [4]f64,
    standing_dead_nitrogen_by_kinetic_g: [4]f64,
    standing_dead_phosphorus_by_kinetic_g: [4]f64,
    standing_dead_aerodynamic_temperature_k: f64,
    standing_dead_aerodynamic_vapor_pressure_kpa: f64,
    standing_dead_surface_temperature_k: f64,
};

pub fn capturePersistentReseedInventories(state: *const State, plant: usize) !PersistentReseedInventories {
    if (plant >= state.plant_seed_storage_carbon_g.len) return error.CanopyPlantIndexOutOfBounds;
    const first = plant * 4;
    var result: PersistentReseedInventories = .{
        .seed_storage = .{ .carbon_g = state.plant_seed_storage_carbon_g[plant], .nitrogen_g = state.plant_seed_storage_nitrogen_g[plant], .phosphorus_g = state.plant_seed_storage_phosphorus_g[plant] },
        .standing_dead = .{ .carbon_g = state.plant_standing_dead_carbon_g[plant], .nitrogen_g = state.plant_standing_dead_nitrogen_g[plant], .phosphorus_g = state.plant_standing_dead_phosphorus_g[plant] },
        .standing_dead_height_m = state.plant_standing_dead_height_m[plant],
        .standing_dead_carbon_by_kinetic_g = undefined,
        .standing_dead_nitrogen_by_kinetic_g = undefined,
        .standing_dead_phosphorus_by_kinetic_g = undefined,
        .standing_dead_aerodynamic_temperature_k = state.plant_standing_dead_aerodynamic_temperature_k[plant],
        .standing_dead_aerodynamic_vapor_pressure_kpa = state.plant_standing_dead_aerodynamic_vapor_pressure_kpa[plant],
        .standing_dead_surface_temperature_k = state.plant_standing_dead_surface_temperature_k[plant],
    };
    @memcpy(&result.standing_dead_carbon_by_kinetic_g, state.plant_standing_dead_carbon_by_kinetic_g[first..][0..4]);
    @memcpy(&result.standing_dead_nitrogen_by_kinetic_g, state.plant_standing_dead_nitrogen_by_kinetic_g[first..][0..4]);
    @memcpy(&result.standing_dead_phosphorus_by_kinetic_g, state.plant_standing_dead_phosphorus_by_kinetic_g[first..][0..4]);
    inline for (@typeInfo(PersistentReseedInventories).@"struct".fields) |field| switch (field.type) {
        f64 => if (!std.math.isFinite(@field(result, field.name))) return error.NonFinitePersistentReseedInventory,
        [4]f64 => for (@field(result, field.name)) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPersistentReseedInventory,
        ElementalMass => inline for (@typeInfo(ElementalMass).@"struct".fields) |mass_field| if (!std.math.isFinite(@field(@field(result, field.name), mass_field.name)) or @field(@field(result, field.name), mass_field.name) < 0) return error.InvalidPersistentReseedInventory,
        else => unreachable,
    };
    return result;
}

pub fn restorePersistentReseedInventories(state: *State, plant: usize, inventories: PersistentReseedInventories) !void {
    if (plant >= state.plant_seed_storage_carbon_g.len) return error.CanopyPlantIndexOutOfBounds;
    state.plant_seed_storage_carbon_g[plant] = inventories.seed_storage.carbon_g;
    state.plant_seed_storage_nitrogen_g[plant] = inventories.seed_storage.nitrogen_g;
    state.plant_seed_storage_phosphorus_g[plant] = inventories.seed_storage.phosphorus_g;
    state.plant_standing_dead_carbon_g[plant] = inventories.standing_dead.carbon_g;
    state.plant_standing_dead_nitrogen_g[plant] = inventories.standing_dead.nitrogen_g;
    state.plant_standing_dead_phosphorus_g[plant] = inventories.standing_dead.phosphorus_g;
    state.plant_standing_dead_height_m[plant] = inventories.standing_dead_height_m;
    const first = plant * 4;
    @memcpy(state.plant_standing_dead_carbon_by_kinetic_g[first..][0..4], &inventories.standing_dead_carbon_by_kinetic_g);
    @memcpy(state.plant_standing_dead_nitrogen_by_kinetic_g[first..][0..4], &inventories.standing_dead_nitrogen_by_kinetic_g);
    @memcpy(state.plant_standing_dead_phosphorus_by_kinetic_g[first..][0..4], &inventories.standing_dead_phosphorus_by_kinetic_g);
    state.plant_standing_dead_aerodynamic_temperature_k[plant] = inventories.standing_dead_aerodynamic_temperature_k;
    state.plant_standing_dead_aerodynamic_vapor_pressure_kpa[plant] = inventories.standing_dead_aerodynamic_vapor_pressure_kpa;
    state.plant_standing_dead_surface_temperature_k[plant] = inventories.standing_dead_surface_temperature_k;
}

pub const Range = struct { first: usize, end: usize };

pub const BranchMobilePoolFluxes = struct {
    fixed_carbon_g: f64,
    maintenance_respiration_demand_g_c: f64,
    available_respirable_carbon_g_c: f64,
    growth_and_respiration_g_c: f64,
    nitrogen_assimilation_respiration_g_c: f64,
    assimilated_nitrogen_g: f64,
    canopy_ammonia_exchange_g_n: f64,
    assimilated_phosphorus_g: f64,
};

pub const BranchMobilePools = struct { carbon_g_c: f64, nitrogen_g_n: f64, phosphorus_g_p: f64 };

pub fn previewBranchMobilePools(state: *const State, branch: usize, fluxes: BranchMobilePoolFluxes) !BranchMobilePools {
    if (branch >= state.branch_mobile_carbon_g.len) return error.CanopyBranchIndexOutOfBounds;
    inline for (@typeInfo(BranchMobilePoolFluxes).@"struct".fields) |field| if (!std.math.isFinite(@field(fluxes, field.name))) return error.NonFiniteBranchPoolFlux;
    const carbon = state.branch_mobile_carbon_g[branch] + fluxes.fixed_carbon_g - @min(fluxes.maintenance_respiration_demand_g_c, fluxes.available_respirable_carbon_g_c) - fluxes.growth_and_respiration_g_c - fluxes.nitrogen_assimilation_respiration_g_c;
    const nitrogen = state.branch_mobile_nitrogen_g[branch] - fluxes.assimilated_nitrogen_g + fluxes.canopy_ammonia_exchange_g_n;
    const phosphorus = state.branch_mobile_phosphorus_g[branch] - fluxes.assimilated_phosphorus_g;
    if (!std.math.isFinite(carbon) or !std.math.isFinite(nitrogen) or !std.math.isFinite(phosphorus)) {
        std.log.err("non-finite branch mobile-pool result: branch={d} carbon_g_c={e} nitrogen_g_n={e} phosphorus_g_p={e}", .{ branch, carbon, nitrogen, phosphorus });
        return error.NonFiniteBranchMobilePoolResult;
    }
    if (carbon < -1e-12 or nitrogen < -1e-12 or phosphorus < -1e-12) {
        std.log.err("branch mobile pool exhausted: branch={d} carbon_g={e} nitrogen_g={e} phosphorus_g={e}", .{ branch, carbon, nitrogen, phosphorus });
        return error.BranchMobilePoolExhausted;
    }
    return .{ .carbon_g_c = @max(0, carbon), .nitrogen_g_n = @max(0, nitrogen), .phosphorus_g_p = @max(0, phosphorus) };
}

pub fn updateBranchMobilePools(state: *State, branch: usize, fluxes: BranchMobilePoolFluxes) !void {
    const pools = try previewBranchMobilePools(state, branch, fluxes);
    state.branch_mobile_carbon_g[branch] = pools.carbon_g_c;
    state.branch_mobile_nitrogen_g[branch] = pools.nitrogen_g_n;
    state.branch_mobile_phosphorus_g[branch] = pools.phosphorus_g_p;
}

pub const C4CarbonFluxes = struct {
    mesophyll_to_bundle_sheath_g_c: f64,
    bundle_sheath_decarboxylation_g_c: f64,
    bundle_sheath_co2_leakage_g_c: f64,
};

pub const C4CarbonParameters = struct {
    bundle_sheath_water_g_per_g_c: f64,
    mesophyll_water_g_per_g_c: f64,
    co2_concentration_umol_per_l_per_g_c_per_g_leaf_c: f64,
    decarboxylation_fraction_per_h: f64,
    co2_decarboxylation_inhibition_umol_per_l: f64,
    decarboxylated_co2_fraction: f64,
    leakage_g_c_per_umol_per_l_g_leaf_c_h: f64,
    mesophyll_feedback_half_saturation_umol_per_l: f64,
    co2_compensation_umol_per_l: f64,
    electron_requirement_umol_e_per_umol_co2: f64,

    pub fn validate(self: C4CarbonParameters) !void {
        inline for (@typeInfo(C4CarbonParameters).@"struct".fields) |field| {
            const value = @field(self, field.name);
            if (!std.math.isFinite(value) or value < 0) return error.InvalidC4CarbonParameter;
        }
        if (self.bundle_sheath_water_g_per_g_c <= 0 or self.mesophyll_water_g_per_g_c <= 0 or self.co2_concentration_umol_per_l_per_g_c_per_g_leaf_c <= 0 or self.co2_decarboxylation_inhibition_umol_per_l <= 0 or self.decarboxylated_co2_fraction > 1 or self.mesophyll_feedback_half_saturation_umol_per_l <= 0 or self.co2_compensation_umol_per_l < 0 or self.electron_requirement_umol_e_per_umol_co2 <= 0) return error.InvalidC4CarbonParameter;
    }
};

pub fn sourceC4CarbonParameters() C4CarbonParameters {
    return .{
        .bundle_sheath_water_g_per_g_c = 1.2,
        .mesophyll_water_g_per_g_c = 4.8,
        .co2_concentration_umol_per_l_per_g_c_per_g_leaf_c = 0.083e9,
        .decarboxylation_fraction_per_h = 0.025,
        .co2_decarboxylation_inhibition_umol_per_l = 1000,
        .decarboxylated_co2_fraction = 0.02,
        .leakage_g_c_per_umol_per_l_g_leaf_c_h = 5.0e-7,
        .mesophyll_feedback_half_saturation_umol_per_l = 5.0e6,
        .co2_compensation_umol_per_l = 0.5,
        .electron_requirement_umol_e_per_umol_co2 = 3,
    };
}

/// GROSUB C4 mesophyll↔bundle-sheath transaction for one runtime node.
pub fn advanceC4CarbonPools(state: *State, node: usize, mesophyll_fixation_g_c: f64, bundle_sheath_fixation_g_c: f64, intercellular_co2_umol_per_l: f64, parameters: C4CarbonParameters, timestep_h: f64) !C4CarbonFluxes {
    inline for (.{ mesophyll_fixation_g_c, bundle_sheath_fixation_g_c, intercellular_co2_umol_per_l, timestep_h }) |value| if (!std.math.isFinite(value)) return error.NonFiniteC4CarbonInput;
    try parameters.validate();
    if (node >= state.node_leaf_carbon_g.len) return error.CanopyNodeIndexOutOfBounds;
    if (mesophyll_fixation_g_c < 0 or bundle_sheath_fixation_g_c < 0 or intercellular_co2_umol_per_l < 0 or timestep_h <= 0) return error.InvalidC4CarbonInput;
    const leaf_carbon = state.node_leaf_carbon_g[node];
    if (leaf_carbon <= 0) return error.C4NodeHasNoLeafCarbon;
    const exchange = try c4_mesophyll_bundle_exchange.exchange(.{
        .bundle_sheath_nonstructural_carbon_g_c = state.node_c3_nonstructural_carbon_g[node],
        .mesophyll_nonstructural_carbon_g_c = state.node_c4_mesophyll_nonstructural_carbon_g[node],
        .bundle_sheath_fixation_g_c_per_timestep = bundle_sheath_fixation_g_c,
        .mesophyll_fixation_g_c_per_timestep = mesophyll_fixation_g_c,
        .leaf_carbon_g_c = leaf_carbon,
        .bundle_sheath_water_g_h2o_per_g_c = parameters.bundle_sheath_water_g_per_g_c,
        .mesophyll_water_g_h2o_per_g_c = parameters.mesophyll_water_g_per_g_c,
        .timestep_h = timestep_h,
    });
    var bundle_nonstructural = exchange.bundle_sheath_nonstructural_carbon_g_c;
    const mesophyll_nonstructural = exchange.mesophyll_nonstructural_carbon_g_c;
    const transfer = exchange.mesophyll_to_bundle_sheath_carbon_g_c;
    const bundle_co2_concentration = @max(0.0, parameters.co2_concentration_umol_per_l_per_g_c_per_g_leaf_c * state.node_bundle_sheath_co2_carbon_g[node] / (leaf_carbon * parameters.bundle_sheath_water_g_per_g_c));
    const decarboxylation = parameters.decarboxylation_fraction_per_h * bundle_nonstructural / (1.0 + bundle_co2_concentration / parameters.co2_decarboxylation_inhibition_umol_per_l) * timestep_h;
    bundle_nonstructural -= decarboxylation;
    const decarboxylated_bicarbonate_fraction = 1.0 - parameters.decarboxylated_co2_fraction;
    var co2_carbon = state.node_bundle_sheath_co2_carbon_g[node] + parameters.decarboxylated_co2_fraction * decarboxylation;
    var bicarbonate_carbon = state.node_bundle_sheath_bicarbonate_carbon_g[node] + decarboxylated_bicarbonate_fraction * decarboxylation;
    const leakage = parameters.leakage_g_c_per_umol_per_l_g_leaf_c_h * (bundle_co2_concentration - intercellular_co2_umol_per_l) * leaf_carbon * parameters.bundle_sheath_water_g_per_g_c * timestep_h;
    co2_carbon -= parameters.decarboxylated_co2_fraction * leakage;
    bicarbonate_carbon -= decarboxylated_bicarbonate_fraction * leakage;
    inline for (.{ bundle_nonstructural, mesophyll_nonstructural, co2_carbon, bicarbonate_carbon }) |value| if (value < -1e-12 or !std.math.isFinite(value)) {
        std.log.err("C4 node carbon transaction failed: node={d} bundle_nonstructural_g={e} mesophyll_nonstructural_g={e} co2_carbon_g={e} bicarbonate_carbon_g={e}", .{ node, bundle_nonstructural, mesophyll_nonstructural, co2_carbon, bicarbonate_carbon });
        return error.C4CarbonPoolExhausted;
    };
    state.node_c3_nonstructural_carbon_g[node] = @max(0.0, bundle_nonstructural);
    state.node_c4_mesophyll_nonstructural_carbon_g[node] = @max(0.0, mesophyll_nonstructural);
    state.node_bundle_sheath_co2_carbon_g[node] = @max(0.0, co2_carbon);
    state.node_bundle_sheath_bicarbonate_carbon_g[node] = @max(0.0, bicarbonate_carbon);
    return .{ .mesophyll_to_bundle_sheath_g_c = transfer, .bundle_sheath_decarboxylation_g_c = decarboxylation, .bundle_sheath_co2_leakage_g_c = leakage };
}

pub const LeafGrowth = struct { carbon_g: f64, nitrogen_g: f64, phosphorus_g: f64 };

pub const Organ = enum(u8) { leaf, sheath, stalk, reserve, husk, ear, grain };

pub const organ_count = @typeInfo(Organ).@"enum".fields.len;

pub const OrganGrowth = struct {
    carbon_g: [organ_count]f64,
    nitrogen_g: [organ_count]f64,
    phosphorus_g: [organ_count]f64,
    total_shoot_carbon_production_g: f64,

    pub fn value(self: OrganGrowth, organ: Organ) LeafGrowth {
        const index = @intFromEnum(organ);
        return .{ .carbon_g = self.carbon_g[index], .nitrogen_g = self.nitrogen_g[index], .phosphorus_g = self.phosphorus_g[index] };
    }
};

/// GROSUB PART(1:7), organ growth yields, and organ N:C/P:C ratios.
pub fn calculateOrganGrowth(total_growth_carbon_consumption_g: f64, partition_fraction: [organ_count]f64, carbon_growth_yield_g_per_g_consumed: [organ_count]f64, nitrogen_to_carbon_g_per_g: [organ_count]f64, phosphorus_to_carbon_g_per_g: [organ_count]f64, minimum_leaf_nutrient_fraction: f64, nutrient_growth_constraint: f64, shoot_growth_yield_g_per_g_consumed: f64) !OrganGrowth {
    inline for (.{ total_growth_carbon_consumption_g, minimum_leaf_nutrient_fraction, nutrient_growth_constraint, shoot_growth_yield_g_per_g_consumed }) |value| if (!std.math.isFinite(value)) return error.NonFiniteOrganGrowthInput;
    if (total_growth_carbon_consumption_g < 0 or minimum_leaf_nutrient_fraction < 0 or minimum_leaf_nutrient_fraction > 1 or nutrient_growth_constraint < 0 or nutrient_growth_constraint > 1 or shoot_growth_yield_g_per_g_consumed < 0) return error.InvalidOrganGrowthInput;
    var result: OrganGrowth = .{ .carbon_g = @splat(0), .nitrogen_g = @splat(0), .phosphorus_g = @splat(0), .total_shoot_carbon_production_g = total_growth_carbon_consumption_g * shoot_growth_yield_g_per_g_consumed };
    var partition_sum: f64 = 0;
    for (0..organ_count) |index| {
        inline for (.{ partition_fraction[index], carbon_growth_yield_g_per_g_consumed[index], nitrogen_to_carbon_g_per_g[index], phosphorus_to_carbon_g_per_g[index] }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidOrganGrowthInput;
        partition_sum += partition_fraction[index];
        result.carbon_g[index] = partition_fraction[index] * total_growth_carbon_consumption_g * carbon_growth_yield_g_per_g_consumed[index];
        const nutrient_factor = if (index == @intFromEnum(Organ.leaf)) minimum_leaf_nutrient_fraction + (1.0 - minimum_leaf_nutrient_fraction) * nutrient_growth_constraint else 1.0;
        result.nitrogen_g[index] = result.carbon_g[index] * nitrogen_to_carbon_g_per_g[index] * nutrient_factor;
        result.phosphorus_g[index] = result.carbon_g[index] * phosphorus_to_carbon_g_per_g[index] * nutrient_factor;
        inline for (.{ result.carbon_g[index], result.nitrogen_g[index], result.phosphorus_g[index] }) |value| if (!std.math.isFinite(value)) return error.NonFiniteOrganGrowthResult;
    }
    if (!std.math.isFinite(partition_sum) or @abs(partition_sum - 1.0) > 1e-8) return error.OrganPartitionDoesNotSumToOne;
    if (!std.math.isFinite(result.total_shoot_carbon_production_g)) return error.NonFiniteOrganGrowthResult;
    return result;
}

/// Commits the branch totals updated immediately after organ partitioning.
/// Grain allocation remains a later reserve→grain transaction, as in GROSUB.
pub fn applyBranchOrganGrowth(state: *State, branch: usize, growth: OrganGrowth) !void {
    try validateBranchOrganGrowthTransaction(state, branch, growth);
    try branch_organ_growth_publication.publish(.{
        .leaf_carbon_g_c = state.branch_leaf_carbon_g,
        .sheath_carbon_g_c = state.branch_sheath_carbon_g,
        .stalk_carbon_g_c = state.branch_stalk_carbon_g,
        .reserve_carbon_g_c = state.branch_reserve_carbon_g,
        .husk_carbon_g_c = state.branch_husk_carbon_g,
        .ear_carbon_g_c = state.branch_ear_carbon_g,
        .leaf_nitrogen_g_n = state.branch_leaf_nitrogen_g,
        .sheath_nitrogen_g_n = state.branch_sheath_nitrogen_g,
        .stalk_nitrogen_g_n = state.branch_stalk_nitrogen_g,
        .reserve_nitrogen_g_n = state.branch_reserve_nitrogen_g,
        .husk_nitrogen_g_n = state.branch_husk_nitrogen_g,
        .ear_nitrogen_g_n = state.branch_ear_nitrogen_g,
        .leaf_phosphorus_g_p = state.branch_leaf_phosphorus_g,
        .sheath_phosphorus_g_p = state.branch_sheath_phosphorus_g,
        .stalk_phosphorus_g_p = state.branch_stalk_phosphorus_g,
        .reserve_phosphorus_g_p = state.branch_reserve_phosphorus_g,
        .husk_phosphorus_g_p = state.branch_husk_phosphorus_g,
        .ear_phosphorus_g_p = state.branch_ear_phosphorus_g,
    }, branch, .{
        .carbon_g_c_per_timestep = organGrowthElement(growth.carbon_g),
        .nitrogen_g_n_per_timestep = organGrowthElement(growth.nitrogen_g),
        .phosphorus_g_p_per_timestep = organGrowthElement(growth.phosphorus_g),
    });
}

fn organGrowthElement(values: [organ_count]f64) branch_organ_growth_publication.ElementGrowth {
    return .{
        .leaf = values[@intFromEnum(Organ.leaf)],
        .sheath = values[@intFromEnum(Organ.sheath)],
        .stalk = values[@intFromEnum(Organ.stalk)],
        .reserve = values[@intFromEnum(Organ.reserve)],
        .husk = values[@intFromEnum(Organ.husk)],
        .ear = values[@intFromEnum(Organ.ear)],
    };
}

fn organGrowthMappings() @TypeOf(.{
    .{ Organ.leaf, "branch_leaf_carbon_g", "branch_leaf_nitrogen_g", "branch_leaf_phosphorus_g" },
    .{ Organ.sheath, "branch_sheath_carbon_g", "branch_sheath_nitrogen_g", "branch_sheath_phosphorus_g" },
    .{ Organ.stalk, "branch_stalk_carbon_g", "branch_stalk_nitrogen_g", "branch_stalk_phosphorus_g" },
    .{ Organ.reserve, "branch_reserve_carbon_g", "branch_reserve_nitrogen_g", "branch_reserve_phosphorus_g" },
    .{ Organ.husk, "branch_husk_carbon_g", "branch_husk_nitrogen_g", "branch_husk_phosphorus_g" },
    .{ Organ.ear, "branch_ear_carbon_g", "branch_ear_nitrogen_g", "branch_ear_phosphorus_g" },
}) {
    return .{
        .{ Organ.leaf, "branch_leaf_carbon_g", "branch_leaf_nitrogen_g", "branch_leaf_phosphorus_g" },
        .{ Organ.sheath, "branch_sheath_carbon_g", "branch_sheath_nitrogen_g", "branch_sheath_phosphorus_g" },
        .{ Organ.stalk, "branch_stalk_carbon_g", "branch_stalk_nitrogen_g", "branch_stalk_phosphorus_g" },
        .{ Organ.reserve, "branch_reserve_carbon_g", "branch_reserve_nitrogen_g", "branch_reserve_phosphorus_g" },
        .{ Organ.husk, "branch_husk_carbon_g", "branch_husk_nitrogen_g", "branch_husk_phosphorus_g" },
        .{ Organ.ear, "branch_ear_carbon_g", "branch_ear_nitrogen_g", "branch_ear_phosphorus_g" },
    };
}

pub fn validateBranchOrganGrowthTransaction(state: *const State, branch: usize, growth: OrganGrowth) !void {
    if (branch >= state.branch_leaf_carbon_g.len) return error.CanopyBranchIndexOutOfBounds;
    const mappings = organGrowthMappings();
    // Validate the complete transaction before publishing any organ. This
    // prevents a late overflow from leaving earlier organs committed.
    inline for (mappings) |mapping| {
        const index = @intFromEnum(mapping[0]);
        inline for (.{
            @field(state, mapping[1])[branch] + growth.carbon_g[index],
            @field(state, mapping[2])[branch] + growth.nitrogen_g[index],
            @field(state, mapping[3])[branch] + growth.phosphorus_g[index],
        }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidBranchOrganGrowthTransaction;
    }
}

/// Distributes growth across the latest runtime nodes, replacing the historical
/// modulo-25 ring while retaining its equal allocation and SLA equation.
pub fn distributeLeafGrowth(state: *State, branch: usize, newest_node_within_branch: usize, first_growing_node_within_branch: usize, maximum_concurrently_growing_nodes: usize, growth: LeafGrowth, protein_per_nitrogen_g_per_g_n: f64, protein_per_phosphorus_g_per_g_p: f64, etoliation_factor: f64, base_specific_leaf_area_m2_per_g_c: f64, minimum_leaf_carbon_per_cell_g: f64, plant_density_per_m2: f64, leaf_area_exponent: f64, turgor_expansion_fraction: f64) !void {
    inline for (.{ growth.carbon_g, growth.nitrogen_g, growth.phosphorus_g, protein_per_nitrogen_g_per_g_n, protein_per_phosphorus_g_per_g_p, etoliation_factor, base_specific_leaf_area_m2_per_g_c, minimum_leaf_carbon_per_cell_g, plant_density_per_m2, leaf_area_exponent, turgor_expansion_fraction }) |value| if (!std.math.isFinite(value)) return error.NonFiniteLeafGrowthInput;
    if (growth.carbon_g < 0 or growth.nitrogen_g < 0 or growth.phosphorus_g < 0 or protein_per_nitrogen_g_per_g_n < 0 or protein_per_phosphorus_g_per_g_p < 0 or etoliation_factor < 0 or base_specific_leaf_area_m2_per_g_c < 0 or minimum_leaf_carbon_per_cell_g < 0 or plant_density_per_m2 <= 0 or turgor_expansion_fraction < 0 or maximum_concurrently_growing_nodes == 0) return error.InvalidLeafGrowthInput;
    const nodes = try state.nodeRange(branch);
    if (newest_node_within_branch >= nodes.end - nodes.first or first_growing_node_within_branch > newest_node_within_branch) return error.CanopyNodeIndexOutOfBounds;
    const first = @max(first_growing_node_within_branch, newest_node_within_branch + 1 -| maximum_concurrently_growing_nodes);
    const count = newest_node_within_branch - first + 1;
    const allocation = 1.0 / @as(f64, @floatFromInt(count));
    try leaf_node_growth_publication.publish(.{
        .leaf_carbon_g_c = state.node_leaf_carbon_g[nodes.first..nodes.end],
        .leaf_nitrogen_g_n = state.node_leaf_nitrogen_g[nodes.first..nodes.end],
        .leaf_phosphorus_g_p = state.node_leaf_phosphorus_g[nodes.first..nodes.end],
        .leaf_protein_g = state.node_leaf_protein_g[nodes.first..nodes.end],
        .leaf_area_m2 = state.node_leaf_area_m2[nodes.first..nodes.end],
        .branch_leaf_area_m2 = &state.branch_leaf_area_m2[branch],
    }, .{
        .first_node = first,
        .last_node = newest_node_within_branch,
        .carbon_growth_g_c_per_node = allocation * growth.carbon_g,
        .nitrogen_growth_g_n_per_node = allocation * growth.nitrogen_g,
        .phosphorus_growth_g_p_per_node = allocation * growth.phosphorus_g,
        .protein_per_nitrogen_g_per_g_n = protein_per_nitrogen_g_per_g_n,
        .protein_per_phosphorus_g_per_g_p = protein_per_phosphorus_g_per_g_p,
        .etiolation_factor = etoliation_factor,
        .base_specific_leaf_area_m2_per_g_c = base_specific_leaf_area_m2_per_g_c,
        .minimum_leaf_carbon_g_c = minimum_leaf_carbon_per_cell_g,
        .plant_population = plant_density_per_m2,
        .leaf_mass_exponent = leaf_area_exponent,
        .turgor_expansion_fraction = turgor_expansion_fraction,
    });
}

pub fn distributeSheathGrowth(state: *State, branch: usize, newest_node_within_branch: usize, first_growing_node_within_branch: usize, maximum_concurrently_growing_nodes: usize, growth: LeafGrowth, protein_per_nitrogen_g_per_g_n: f64, protein_per_phosphorus_g_per_g_p: f64, etoliation_factor: f64, base_specific_length_m_per_g_c: f64, minimum_sheath_carbon_per_cell_g: f64, plant_density_per_m2: f64, length_exponent: f64, turgor_expansion_fraction: f64, vertical_projection_fraction: f64) !void {
    inline for (.{ growth.carbon_g, growth.nitrogen_g, growth.phosphorus_g, protein_per_nitrogen_g_per_g_n, protein_per_phosphorus_g_per_g_p, etoliation_factor, base_specific_length_m_per_g_c, minimum_sheath_carbon_per_cell_g, plant_density_per_m2, length_exponent, turgor_expansion_fraction, vertical_projection_fraction }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSheathGrowthInput;
    if (growth.carbon_g < 0 or growth.nitrogen_g < 0 or growth.phosphorus_g < 0 or protein_per_nitrogen_g_per_g_n < 0 or protein_per_phosphorus_g_per_g_p < 0 or etoliation_factor < 0 or base_specific_length_m_per_g_c < 0 or minimum_sheath_carbon_per_cell_g < 0 or plant_density_per_m2 <= 0 or turgor_expansion_fraction < 0 or vertical_projection_fraction < 0 or maximum_concurrently_growing_nodes == 0) return error.InvalidSheathGrowthInput;
    const nodes = try state.nodeRange(branch);
    if (newest_node_within_branch >= nodes.end - nodes.first or first_growing_node_within_branch > newest_node_within_branch) return error.CanopyNodeIndexOutOfBounds;
    const first = @max(first_growing_node_within_branch, newest_node_within_branch + 1 -| maximum_concurrently_growing_nodes);
    const count = newest_node_within_branch - first + 1;
    const allocation = 1.0 / @as(f64, @floatFromInt(count));
    for (first..newest_node_within_branch + 1) |local_node| {
        const node = nodes.first + local_node;
        const carbon = allocation * growth.carbon_g;
        const nitrogen = allocation * growth.nitrogen_g;
        const phosphorus = allocation * growth.phosphorus_g;
        const updated_sheath_carbon = state.node_sheath_carbon_g[node] + carbon;
        state.node_sheath_carbon_g[node] = updated_sheath_carbon;
        state.node_sheath_nitrogen_g[node] += nitrogen;
        state.node_sheath_phosphorus_g[node] += phosphorus;
        state.node_sheath_protein_g[node] += @min(nitrogen * protein_per_nitrogen_g_per_g_n, phosphorus * protein_per_phosphorus_g_per_g_p);
        if (state.node_leaf_carbon_g[node] > 0) {
            const specific_length = etoliation_factor * base_specific_length_m_per_g_c * std.math.pow(f64, @max(minimum_sheath_carbon_per_cell_g, updated_sheath_carbon) / plant_density_per_m2, length_exponent) * turgor_expansion_fraction;
            state.node_sheath_height_m[node] += carbon / plant_density_per_m2 * specific_length * vertical_projection_fraction;
        }
    }
}

pub const StalkGrowthResult = struct { stem_diameter_m: f64 };

pub const LayerLeafOutputs = struct {
    area_m2: []f64,
    carbon_g: []f64,
    nitrogen_g: []f64,
    phosphorus_g: []f64,
};

/// GROSUB ARLFL/WGLFL/WGLFLN/WGLFLP allocation for one leafed node. Layer
/// boundaries and inclination classes are runtime extents; no 25-node ring or
/// fixed JC/N dimensions are retained.
pub fn allocateLeafAcrossCanopyLayers(leaf_area_m2: f64, leaf_carbon_g: f64, leaf_nitrogen_g: f64, leaf_phosphorus_g: f64, population_per_m2: f64, leaf_length_to_width_ratio: f64, stalk_height_m: f64, sheath_height_m: f64, maximum_canopy_height_m: f64, layer_boundary_height_m: []const f64, inclination_sine: []const f64, inclination_fraction: []const f64, outputs: LayerLeafOutputs) !f64 {
    inline for (.{ leaf_area_m2, leaf_carbon_g, leaf_nitrogen_g, leaf_phosphorus_g, population_per_m2, leaf_length_to_width_ratio, stalk_height_m, sheath_height_m, maximum_canopy_height_m }) |value| if (!std.math.isFinite(value)) return error.NonFiniteCanopyLayerAllocationInput;
    if (leaf_area_m2 < 0 or leaf_carbon_g < 0 or leaf_nitrogen_g < 0 or leaf_phosphorus_g < 0 or population_per_m2 <= 0 or leaf_length_to_width_ratio < 0 or stalk_height_m < 0 or sheath_height_m < 0 or maximum_canopy_height_m < 0 or layer_boundary_height_m.len < 2 or inclination_sine.len == 0 or inclination_sine.len != inclination_fraction.len) return error.InvalidCanopyLayerAllocationInput;
    const layer_count = layer_boundary_height_m.len - 1;
    inline for (.{ outputs.area_m2, outputs.carbon_g, outputs.nitrogen_g, outputs.phosphorus_g }) |values| if (values.len != layer_count) return error.CanopyLayerAllocationDimensionMismatch;
    for (1..layer_boundary_height_m.len) |index| if (!std.math.isFinite(layer_boundary_height_m[index - 1]) or !std.math.isFinite(layer_boundary_height_m[index]) or layer_boundary_height_m[index] < layer_boundary_height_m[index - 1]) return error.InvalidCanopyLayerBoundary;
    var inclination_total: f64 = 0;
    for (inclination_sine, inclination_fraction) |sine, fraction| {
        if (!std.math.isFinite(sine) or sine < 0 or sine > 1 or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidCanopyInclinationDistribution;
        inclination_total += fraction;
    }
    if (@abs(inclination_total - 1) > 1.0e-4) return error.CanopyInclinationDistributionDoesNotSumToOne;
    @memset(outputs.area_m2, 0);
    @memset(outputs.carbon_g, 0);
    @memset(outputs.nitrogen_g, 0);
    @memset(outputs.phosphorus_g, 0);
    if (leaf_area_m2 == 0) return @min(maximum_canopy_height_m + 0.01, stalk_height_m + sheath_height_m);

    const leaf_length_m = @sqrt(leaf_length_to_width_ratio * leaf_area_m2 / population_per_m2);
    const leaf_base_height_m = stalk_height_m + sheath_height_m;
    const height_cap_m = maximum_canopy_height_m + 0.01;
    var accumulated_leaf_elevation_m: f64 = 0;
    var highest_leaf_height_m: f64 = 0;
    var inclination = inclination_sine.len;
    while (inclination > 0) {
        inclination -= 1;
        const class_fraction = inclination_fraction[inclination];
        if (class_fraction == 0) continue;
        const elevation_m = inclination_sine[inclination] * class_fraction * leaf_length_m;
        const lower_m = @min(height_cap_m - elevation_m, leaf_base_height_m + accumulated_leaf_elevation_m);
        const upper_m = @min(height_cap_m, lower_m + elevation_m);
        if (upper_m <= lower_m) {
            const layer = containingCanopyLayer(layer_boundary_height_m, @max(0.0, lower_m));
            addLeafLayer(outputs, layer, class_fraction, leaf_area_m2, leaf_carbon_g, leaf_nitrogen_g, leaf_phosphorus_g);
        } else {
            var allocated_fraction: f64 = 0;
            for (0..layer_count) |layer| {
                const overlap_m = @max(0.0, @min(upper_m, layer_boundary_height_m[layer + 1]) - @max(lower_m, layer_boundary_height_m[layer]));
                if (overlap_m == 0) continue;
                const fraction = class_fraction * overlap_m / (upper_m - lower_m);
                addLeafLayer(outputs, layer, fraction, leaf_area_m2, leaf_carbon_g, leaf_nitrogen_g, leaf_phosphorus_g);
                allocated_fraction += fraction;
            }
            if (allocated_fraction < class_fraction) {
                const layer = containingCanopyLayer(layer_boundary_height_m, std.math.clamp(lower_m, layer_boundary_height_m[0], layer_boundary_height_m[layer_count]));
                addLeafLayer(outputs, layer, class_fraction - allocated_fraction, leaf_area_m2, leaf_carbon_g, leaf_nitrogen_g, leaf_phosphorus_g);
            }
        }
        accumulated_leaf_elevation_m += elevation_m;
        highest_leaf_height_m = @max(highest_leaf_height_m, upper_m);
    }
    return highest_leaf_height_m;
}

fn containingCanopyLayer(boundaries: []const f64, height_m: f64) usize {
    for (0..boundaries.len - 1) |layer| if (height_m <= boundaries[layer + 1]) return layer;
    return boundaries.len - 2;
}

fn addLeafLayer(outputs: LayerLeafOutputs, layer: usize, fraction: f64, area_m2: f64, carbon_g: f64, nitrogen_g: f64, phosphorus_g: f64) void {
    outputs.area_m2[layer] += fraction * area_m2;
    outputs.carbon_g[layer] += fraction * carbon_g;
    outputs.nitrogen_g[layer] += fraction * nitrogen_g;
    outputs.phosphorus_g[layer] += fraction * phosphorus_g;
}

pub const StalkLayerAllocation = struct {
    radius_m: f64,
    surface_area_m2: f64,
    sapwood_carbon_g: f64,
};

pub fn accumulatePotentialSeedSites(current_site_count: f64, stem_elongation_started: bool, anthesis_started: bool, branch_shoot_carbon_g: f64, canopy_shoot_carbon_g: f64, canopy_shoot_growth_g_c_per_step: f64, potential_sites_per_g_growth: f64, structural_presence_threshold_g: f64) !f64 {
    inline for (.{ current_site_count, branch_shoot_carbon_g, canopy_shoot_carbon_g, canopy_shoot_growth_g_c_per_step, potential_sites_per_g_growth, structural_presence_threshold_g }) |value| if (!std.math.isFinite(value)) return error.NonFinitePotentialSeedSiteInput;
    if (current_site_count < 0 or branch_shoot_carbon_g < 0 or canopy_shoot_carbon_g < 0 or canopy_shoot_growth_g_c_per_step < 0 or potential_sites_per_g_growth < 0 or structural_presence_threshold_g < 0) return error.InvalidPotentialSeedSiteInput;
    if (!stem_elongation_started or anthesis_started or canopy_shoot_carbon_g <= structural_presence_threshold_g or canopy_shoot_growth_g_c_per_step <= structural_presence_threshold_g) return current_site_count;
    const branch_growth_g_c = canopy_shoot_growth_g_c_per_step * branch_shoot_carbon_g / canopy_shoot_carbon_g;
    return current_site_count + potential_sites_per_g_growth * branch_growth_g_c;
}

pub const SeedSetInputs = struct {
    anthesis_started: bool,
    grain_fill_started: bool,
    final_seed_number_set: bool,
    maximum_seed_size_set: bool,
    mobile_carbon_concentration_g_per_g: f64,
    mobile_nitrogen_concentration_g_per_g: f64,
    mobile_phosphorus_concentration_g_per_g: f64,
    carbon_half_saturation_g_per_g: f64,
    nitrogen_half_saturation_g_per_g: f64,
    phosphorus_half_saturation_g_per_g: f64,
    canopy_temperature_c: f64,
    chilling_temperature_c: f64,
    high_temperature_c: f64,
    seed_loss_fraction_per_c_h: f64,
    timestep_h: f64,
    water_growth_fraction: f64,
    reproductive_stage_increment: f64,
    maximum_seeds_per_site: f64,
    potential_site_count: f64,
    current_seed_count: f64,
    maximum_individual_seed_carbon_g: f64,
    current_individual_seed_carbon_g: f64,
};

pub const SeedSetParameters = struct {
    carbon_half_saturation_g_per_g: f64,
    nitrogen_half_saturation_g_per_g: f64,
    phosphorus_half_saturation_g_per_g: f64,

    pub fn validate(self: SeedSetParameters) !void {
        inline for (@typeInfo(SeedSetParameters).@"struct".fields) |field| {
            const value = @field(self, field.name);
            if (!std.math.isFinite(value) or value <= 0) return error.InvalidSeedSetParameter;
        }
    }
};

pub fn compatibilitySeedSetParameters() SeedSetParameters {
    return .{ .carbon_half_saturation_g_per_g = 2.5e-2, .nitrogen_half_saturation_g_per_g = 0.5e-2, .phosphorus_half_saturation_g_per_g = 0.1e-2 };
}

pub const CanopyWaterGrowthResponse = struct { stomatal_fraction: f64, growth_fraction: f64, turgor_expansion_fraction: f64, water_potential_expansion_fraction: f64 };

pub fn canopyWaterGrowthResponse(shallow_root_profile: bool, canopy_turgor_potential_megapascal: f64, minimum_turgor_potential_megapascal: f64, canopy_total_water_potential_megapascal: f64, stomatal_turgor_shape: f64) !CanopyWaterGrowthResponse {
    const response = try @import("../energy/water_stress_response.zig").calculate(.{
        .root_profile = if (shallow_root_profile) .shallow else .non_shallow,
        .canopy_turgor_potential_megapascal = canopy_turgor_potential_megapascal,
        .minimum_canopy_turgor_potential_megapascal = minimum_turgor_potential_megapascal,
        .canopy_water_potential_megapascal = canopy_total_water_potential_megapascal,
        .stomatal_turgor_shape_per_megapascal = stomatal_turgor_shape,
    });
    return .{
        .stomatal_fraction = response.stomatal_resistance_factor,
        .growth_fraction = response.growth_factor,
        .turgor_expansion_fraction = response.turgor_expansion_factor,
        .water_potential_expansion_fraction = response.water_potential_expansion_factor,
    };
}

pub const SeedSetResult = struct {
    nutrient_set_fraction: f64,
    thermal_loss_fraction: f64,
    seed_count: f64,
    individual_seed_carbon_g: f64,
};

/// GROSUB SET/FGRNX/GRNOB/GRWTB update following anthesis.
pub fn updateSeedNumberAndSize(input: SeedSetInputs) !SeedSetResult {
    inline for (@typeInfo(SeedSetInputs).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(input, field.name))) return error.NonFiniteSeedSetInput;
    inline for (.{ input.mobile_carbon_concentration_g_per_g, input.mobile_nitrogen_concentration_g_per_g, input.mobile_phosphorus_concentration_g_per_g, input.carbon_half_saturation_g_per_g, input.nitrogen_half_saturation_g_per_g, input.phosphorus_half_saturation_g_per_g, input.seed_loss_fraction_per_c_h, input.timestep_h, input.water_growth_fraction, input.reproductive_stage_increment, input.maximum_seeds_per_site, input.potential_site_count, input.current_seed_count, input.maximum_individual_seed_carbon_g, input.current_individual_seed_carbon_g }) |value| if (value < 0) return error.InvalidSeedSetInput;
    if (input.timestep_h <= 0 or input.water_growth_fraction > 1) return error.InvalidSeedSetInput;
    if (!input.anthesis_started or input.maximum_seed_size_set) return .{
        .nutrient_set_fraction = 0,
        .thermal_loss_fraction = 0,
        .seed_count = input.current_seed_count,
        .individual_seed_carbon_g = input.current_individual_seed_carbon_g,
    };
    if (input.carbon_half_saturation_g_per_g + input.mobile_carbon_concentration_g_per_g <= 0 or input.nitrogen_half_saturation_g_per_g + input.mobile_nitrogen_concentration_g_per_g <= 0 or input.phosphorus_half_saturation_g_per_g + input.mobile_phosphorus_concentration_g_per_g <= 0) return error.InvalidSeedSetInput;
    const nutrient_set_fraction = @min(input.mobile_carbon_concentration_g_per_g / (input.mobile_carbon_concentration_g_per_g + input.carbon_half_saturation_g_per_g), @min(input.mobile_nitrogen_concentration_g_per_g / (input.mobile_nitrogen_concentration_g_per_g + input.nitrogen_half_saturation_g_per_g), input.mobile_phosphorus_concentration_g_per_g / (input.mobile_phosphorus_concentration_g_per_g + input.phosphorus_half_saturation_g_per_g)));
    var thermal_loss_fraction: f64 = 0;
    if (!input.grain_fill_started or !input.final_seed_number_set) {
        if (input.canopy_temperature_c < input.chilling_temperature_c)
            thermal_loss_fraction = input.seed_loss_fraction_per_c_h * (input.chilling_temperature_c - input.canopy_temperature_c) * input.timestep_h
        else if (input.canopy_temperature_c > input.high_temperature_c)
            thermal_loss_fraction = input.seed_loss_fraction_per_c_h * (input.canopy_temperature_c - input.high_temperature_c) * input.timestep_h;
    }
    var seed_count = input.current_seed_count;
    if (input.anthesis_started and !input.final_seed_number_set) {
        const nutrient_water_set = nutrient_set_fraction * std.math.pow(f64, input.water_growth_fraction, 0.25);
        const maximum_seed_count = input.maximum_seeds_per_site * input.potential_site_count;
        const candidate_seed_count = @min(maximum_seed_count, seed_count + maximum_seed_count * nutrient_water_set * input.reproductive_stage_increment - thermal_loss_fraction * seed_count);
        if (candidate_seed_count < 0) return error.NegativeSeedCount;
        seed_count = candidate_seed_count;
    }
    var individual_seed_carbon_g = input.current_individual_seed_carbon_g;
    if (input.grain_fill_started and !input.maximum_seed_size_set) {
        const nutrient_water_set = std.math.pow(f64, nutrient_set_fraction * input.water_growth_fraction, 0.25);
        individual_seed_carbon_g = @min(input.maximum_individual_seed_carbon_g, individual_seed_carbon_g + input.maximum_individual_seed_carbon_g * @max(0.5, nutrient_water_set) * input.reproductive_stage_increment);
    }
    return .{ .nutrient_set_fraction = nutrient_set_fraction, .thermal_loss_fraction = thermal_loss_fraction, .seed_count = seed_count, .individual_seed_carbon_g = individual_seed_carbon_g };
}

/// GROSUB RSTK/ARSTKB/WVSTKB and ARSTK allocation for one branch.
pub fn allocateStalkAcrossCanopyLayers(stalk_carbon_g: f64, retained_stalk_carbon_g: f64, population_per_m2: f64, stalk_volume_m3_per_g_c: f64, branch_base_height_m: f64, branch_tip_height_m: f64, annual_growth_habit: bool, specific_internode_length_positive: bool, layer_boundary_height_m: []const f64, layer_stalk_area_m2: []f64) !StalkLayerAllocation {
    inline for (.{ stalk_carbon_g, retained_stalk_carbon_g, population_per_m2, stalk_volume_m3_per_g_c, branch_base_height_m, branch_tip_height_m }) |value| if (!std.math.isFinite(value)) return error.NonFiniteStalkLayerAllocationInput;
    if (stalk_carbon_g < 0 or retained_stalk_carbon_g < 0 or population_per_m2 <= 0 or stalk_volume_m3_per_g_c <= 0 or branch_base_height_m < 0 or branch_tip_height_m < branch_base_height_m or layer_boundary_height_m.len < 2 or layer_stalk_area_m2.len != layer_boundary_height_m.len - 1) return error.InvalidStalkLayerAllocationInput;
    for (1..layer_boundary_height_m.len) |index| if (!std.math.isFinite(layer_boundary_height_m[index - 1]) or !std.math.isFinite(layer_boundary_height_m[index]) or layer_boundary_height_m[index] < layer_boundary_height_m[index - 1]) return error.InvalidCanopyLayerBoundary;
    @memset(layer_stalk_area_m2, 0);
    const stalk_height_m = branch_tip_height_m - branch_base_height_m;
    if (stalk_height_m == 0) return .{ .radius_m = 0, .surface_area_m2 = 0, .sapwood_carbon_g = if (specific_internode_length_positive) 0 else retained_stalk_carbon_g };
    const radius_m = @sqrt(stalk_volume_m3_per_g_c * (stalk_carbon_g / population_per_m2) / (3.1416 * stalk_height_m));
    const surface_area_m2 = 6.2832 * radius_m * stalk_height_m * population_per_m2;
    const sapwood_carbon_g = if (annual_growth_habit)
        stalk_carbon_g
    else blk: {
        const sapwood_thickness_m = @min(1.0e-3, 0.05 * radius_m);
        const sapwood_cross_section_m2 = 3.1416 * (2 * radius_m * sapwood_thickness_m - sapwood_thickness_m * sapwood_thickness_m);
        break :blk sapwood_cross_section_m2 / stalk_volume_m3_per_g_c * stalk_height_m * population_per_m2;
    };
    for (0..layer_stalk_area_m2.len) |layer| {
        const overlap_m = @max(0.0, @min(branch_tip_height_m, layer_boundary_height_m[layer + 1]) - @max(branch_base_height_m, layer_boundary_height_m[layer]));
        layer_stalk_area_m2[layer] = surface_area_m2 * overlap_m / stalk_height_m;
    }
    return .{ .radius_m = radius_m, .surface_area_m2 = surface_area_m2, .sapwood_carbon_g = sapwood_carbon_g };
}

pub fn distributeStalkGrowth(state: *State, branch: usize, first_growing_node_within_branch: usize, last_growing_node_within_branch: usize, growth: LeafGrowth, etoliation_factor: f64, base_specific_internode_length_m_per_g_c: f64, minimum_stalk_carbon_per_cell_g: f64, plant_density_per_m2: f64, length_exponent: f64, turgor_expansion_fraction: f64, vertical_projection_fraction: f64, stalk_volume_m3_per_g_c: f64) !StalkGrowthResult {
    inline for (.{ growth.carbon_g, growth.nitrogen_g, growth.phosphorus_g, etoliation_factor, base_specific_internode_length_m_per_g_c, minimum_stalk_carbon_per_cell_g, plant_density_per_m2, length_exponent, turgor_expansion_fraction, vertical_projection_fraction, stalk_volume_m3_per_g_c }) |value| if (!std.math.isFinite(value)) return error.NonFiniteStalkGrowthInput;
    if (growth.carbon_g < 0 or growth.nitrogen_g < 0 or growth.phosphorus_g < 0 or etoliation_factor < 0 or base_specific_internode_length_m_per_g_c < 0 or minimum_stalk_carbon_per_cell_g < 0 or plant_density_per_m2 <= 0 or turgor_expansion_fraction < 0 or vertical_projection_fraction < 0 or stalk_volume_m3_per_g_c < 0) return error.InvalidStalkGrowthInput;
    const nodes = try state.nodeRange(branch);
    if (last_growing_node_within_branch >= nodes.end - nodes.first or first_growing_node_within_branch > last_growing_node_within_branch) return error.CanopyNodeIndexOutOfBounds;
    const first = first_growing_node_within_branch;
    const count = last_growing_node_within_branch - first + 1;
    const allocation = 1.0 / @as(f64, @floatFromInt(count));
    for (first..last_growing_node_within_branch + 1) |local_node| {
        const node = nodes.first + local_node;
        const carbon = allocation * growth.carbon_g;
        state.node_internode_carbon_g[node] += carbon;
        state.node_internode_nitrogen_g[node] += allocation * growth.nitrogen_g;
        state.node_internode_phosphorus_g[node] += allocation * growth.phosphorus_g;
        const specific_length = @sqrt(etoliation_factor) * base_specific_internode_length_m_per_g_c * std.math.pow(f64, @max(minimum_stalk_carbon_per_cell_g, state.node_internode_carbon_g[node]) / plant_density_per_m2, length_exponent) * turgor_expansion_fraction;
        state.node_internode_length_m[node] += carbon / plant_density_per_m2 * specific_length * vertical_projection_fraction;
        state.node_height_m[node] = state.node_internode_length_m[node] + if (local_node > 0) state.node_height_m[node - 1] else 0;
    }
    const diagnostic_node = nodes.first + first;
    const previous_height = if (first > 0) state.node_height_m[diagnostic_node - 1] else 0;
    const height_difference = state.node_height_m[diagnostic_node] - previous_height;
    var diameter: f64 = 0;
    if (height_difference > 0 and state.node_internode_carbon_g[diagnostic_node] > 0) {
        const radius = @sqrt(stalk_volume_m3_per_g_c * (state.node_internode_carbon_g[diagnostic_node] / plant_density_per_m2) / (3.1416 * height_difference));
        diameter = 2.0 * radius;
        if (first > 1) diameter *= std.math.pow(f64, @as(f64, @floatFromInt(first)), 0.167);
    }
    return .{ .stem_diameter_m = diameter };
}

pub const RecyclingFractions = struct { carbon: f64, nitrogen: f64, phosphorus: f64 };

pub const NodeSenescenceAllocation = struct {
    leaf_present: bool,
    leaf_fraction: f64,
    sheath_fraction: f64,
    leaf_recycled_carbon_g: f64,
    sheath_recycled_carbon_g: f64,
    carbon_recovered_to_mobile_pool_g: f64,
    carbon_respired_g: f64,
    remaining_respiration_demand_g_c: f64,
};

pub fn allocateNodeSenescenceDemand(node_respiration_demand_g_c: f64, leaf_carbon_g: f64, sheath_carbon_g: f64, carbon_recycling_fraction: f64, phenological_senescence_fraction: f64, nonwoody_carbon_fraction: f64) !NodeSenescenceAllocation {
    return allocateNodeSenescenceDemandWithThreshold(node_respiration_demand_g_c, leaf_carbon_g, sheath_carbon_g, carbon_recycling_fraction, phenological_senescence_fraction, nonwoody_carbon_fraction, 0);
}

fn allocateNodeSenescenceDemandWithThreshold(node_respiration_demand_g_c: f64, leaf_carbon_g: f64, sheath_carbon_g: f64, carbon_recycling_fraction: f64, phenological_senescence_fraction: f64, nonwoody_carbon_fraction: f64, leaf_presence_threshold_g_c: f64) !NodeSenescenceAllocation {
    inline for (.{ node_respiration_demand_g_c, leaf_carbon_g, sheath_carbon_g, carbon_recycling_fraction, phenological_senescence_fraction, nonwoody_carbon_fraction }) |value| if (!std.math.isFinite(value)) return error.NonFiniteNodeSenescenceInput;
    if (node_respiration_demand_g_c < 0 or leaf_carbon_g < 0 or sheath_carbon_g < 0 or carbon_recycling_fraction < 0 or carbon_recycling_fraction > 1 or phenological_senescence_fraction < 0 or phenological_senescence_fraction > 1 or nonwoody_carbon_fraction < 0 or nonwoody_carbon_fraction > 1) return error.InvalidNodeSenescenceInput;
    const request = try node_senescence_remobilization_request.calculate(.{
        .leaf_carbon_g_c = &.{leaf_carbon_g},
        .sheath_carbon_g_c = &.{sheath_carbon_g},
        .leaf_nitrogen_g_n = &.{0},
        .leaf_phosphorus_g_p = &.{0},
    }, 0, node_respiration_demand_g_c, leaf_presence_threshold_g_c, .{
        .carbon = carbon_recycling_fraction,
        .nitrogen = 0,
        .phosphorus = 0,
    });
    const recyclable_leaf = request.remobilizable_leaf_carbon_g_c;
    const leaf_present = leaf_carbon_g > leaf_presence_threshold_g_c;
    const leaf_fraction = if (leaf_present) request.leaf_mass_removal_fraction else 1;
    const consumed_leaf_recycling = request.leaf_mass_removal_fraction * recyclable_leaf * nonwoody_carbon_fraction;
    var remaining = @max(0.0, node_respiration_demand_g_c - consumed_leaf_recycling);
    const recyclable_sheath = sheath_carbon_g * carbon_recycling_fraction;
    const sheath_fraction = if (request.sheath_senescence_respiration_g_c_per_timestep > 0 and sheath_carbon_g > 0)
        (if (recyclable_sheath > leaf_presence_threshold_g_c) std.math.clamp(request.sheath_senescence_respiration_g_c_per_timestep / recyclable_sheath, 0, 1) else 1)
    else
        0;
    const consumed_sheath_recycling = sheath_fraction * recyclable_sheath * nonwoody_carbon_fraction;
    remaining = @max(0.0, remaining - consumed_sheath_recycling);
    const recycled_carbon = consumed_leaf_recycling + consumed_sheath_recycling;
    return .{
        .leaf_present = leaf_present,
        .leaf_fraction = leaf_fraction,
        .sheath_fraction = sheath_fraction,
        .leaf_recycled_carbon_g = leaf_fraction * recyclable_leaf,
        .sheath_recycled_carbon_g = sheath_fraction * recyclable_sheath,
        .carbon_recovered_to_mobile_pool_g = recycled_carbon * phenological_senescence_fraction,
        .carbon_respired_g = recycled_carbon * (1.0 - phenological_senescence_fraction),
        .remaining_respiration_demand_g_c = remaining,
    };
}

pub fn recyclingFractions(emerged: bool, mobile_carbon_g_per_g: f64, mobile_nitrogen_g_per_g: f64, mobile_phosphorus_g_per_g: f64, nitrogen_inhibition_g_n_per_g_c: f64, phosphorus_inhibition_g_p_per_g_c: f64, minimum_carbon_recycling_fraction: f64, responsive_carbon_recycling_fraction: f64, maximum_nitrogen_recycling_fraction: f64, maximum_phosphorus_recycling_fraction: f64) !RecyclingFractions {
    const minimum_carbon = [1]f64{minimum_carbon_recycling_fraction};
    const responsive_carbon = [1]f64{responsive_carbon_recycling_fraction};
    const maximum_nitrogen = [1]f64{maximum_nitrogen_recycling_fraction};
    const maximum_phosphorus = [1]f64{maximum_phosphorus_recycling_fraction};
    const result = try shoot_recycling_fraction.calculate(
        emerged,
        0,
        .{
            .mobile_carbon_g_c_per_g_c = mobile_carbon_g_per_g,
            .mobile_nitrogen_g_n_per_g_c = mobile_nitrogen_g_per_g,
            .mobile_phosphorus_g_p_per_g_c = mobile_phosphorus_g_per_g,
        },
        .{
            .nitrogen_g_n_per_g_c = nitrogen_inhibition_g_n_per_g_c,
            .phosphorus_g_p_per_g_c = phosphorus_inhibition_g_p_per_g_c,
        },
        .{
            .minimum_carbon_fraction = &minimum_carbon,
            .responsive_carbon_fraction = &responsive_carbon,
            .maximum_nitrogen_fraction = &maximum_nitrogen,
            .maximum_phosphorus_fraction = &maximum_phosphorus,
        },
    );
    return .{
        .carbon = result.carbon,
        .nitrogen = result.nitrogen,
        .phosphorus = result.phosphorus,
    };
}

pub const KineticFractions = struct {
    carbon: [4]f64,
    nitrogen: [4]f64,
    phosphorus: [4]f64,

    pub fn validate(self: KineticFractions) !void {
        inline for (.{ self.carbon, self.nitrogen, self.phosphorus }) |fractions| {
            var sum: f64 = 0;
            for (fractions) |fraction| {
                if (!std.math.isFinite(fraction) or fraction < 0) return error.InvalidLitterKineticFraction;
                sum += fraction;
            }
            if (@abs(sum - 1.0) > 1e-8) return error.LitterKineticFractionsDoNotSumToOne;
        }
    }
};

pub const SenescenceProducts = struct {
    woody_carbon_g: [4]f64 = @splat(0),
    woody_nitrogen_g: [4]f64 = @splat(0),
    woody_phosphorus_g: [4]f64 = @splat(0),
    nonwoody_carbon_g: [4]f64 = @splat(0),
    nonwoody_nitrogen_g: [4]f64 = @splat(0),
    nonwoody_phosphorus_g: [4]f64 = @splat(0),
    recycled_carbon_g: f64 = 0,
    recycled_nitrogen_g: f64 = 0,
    recycled_phosphorus_g: f64 = 0,
    respired_carbon_g: f64 = 0,
};

pub fn addSenescenceProducts(total: *SenescenceProducts, addition: SenescenceProducts) void {
    for (0..4) |kinetic| {
        total.woody_carbon_g[kinetic] += addition.woody_carbon_g[kinetic];
        total.woody_nitrogen_g[kinetic] += addition.woody_nitrogen_g[kinetic];
        total.woody_phosphorus_g[kinetic] += addition.woody_phosphorus_g[kinetic];
        total.nonwoody_carbon_g[kinetic] += addition.nonwoody_carbon_g[kinetic];
        total.nonwoody_nitrogen_g[kinetic] += addition.nonwoody_nitrogen_g[kinetic];
        total.nonwoody_phosphorus_g[kinetic] += addition.nonwoody_phosphorus_g[kinetic];
    }
    total.recycled_carbon_g += addition.recycled_carbon_g;
    total.recycled_nitrogen_g += addition.recycled_nitrogen_g;
    total.recycled_phosphorus_g += addition.recycled_phosphorus_g;
    total.respired_carbon_g += addition.respired_carbon_g;
}

pub fn commitNodeSenescenceDemand(state: *State, branch: usize, node_within_branch: usize, allocation: NodeSenescenceAllocation, recycling: RecyclingFractions, protein_per_nitrogen_g_per_g_n: f64, protein_per_phosphorus_g_per_g_p: f64, woody_fraction: [2]f64, leaf_woody_nitrogen_fraction: [2]f64, sheath_woody_nitrogen_fraction: [2]f64, leaf_woody_phosphorus_fraction: [2]f64, sheath_woody_phosphorus_fraction: [2]f64, woody_kinetics: KineticFractions, leaf_kinetics: KineticFractions, sheath_kinetics: KineticFractions) !SenescenceProducts {
    inline for (.{ allocation.leaf_fraction, allocation.sheath_fraction, allocation.carbon_recovered_to_mobile_pool_g, allocation.carbon_respired_g, recycling.carbon, recycling.nitrogen, recycling.phosphorus, protein_per_nitrogen_g_per_g_n, protein_per_phosphorus_g_per_g_p }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSenescenceInput;
    if (allocation.leaf_fraction < 0 or allocation.leaf_fraction > 1 or allocation.sheath_fraction < 0 or allocation.sheath_fraction > 1 or allocation.carbon_recovered_to_mobile_pool_g < 0 or allocation.carbon_respired_g < 0 or recycling.carbon < 0 or recycling.carbon > 1 or recycling.nitrogen < 0 or recycling.nitrogen > 1 or recycling.phosphorus < 0 or recycling.phosphorus > 1 or protein_per_nitrogen_g_per_g_n < 0 or protein_per_phosphorus_g_per_g_p < 0) return error.InvalidSenescenceInput;
    inline for (.{ woody_fraction, leaf_woody_nitrogen_fraction, sheath_woody_nitrogen_fraction, leaf_woody_phosphorus_fraction, sheath_woody_phosphorus_fraction }) |fractions| {
        for (fractions) |fraction| if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidWoodyFraction;
        if (@abs(fractions[0] + fractions[1] - 1.0) > 1e-8) return error.InvalidWoodyFraction;
    }
    try woody_kinetics.validate();
    try leaf_kinetics.validate();
    try sheath_kinetics.validate();
    const nodes = try state.nodeRange(branch);
    if (node_within_branch >= nodes.end - nodes.first) return error.CanopyNodeIndexOutOfBounds;
    const node = nodes.first + node_within_branch;
    const leaf_c = state.node_leaf_carbon_g[node];
    const leaf_n = state.node_leaf_nitrogen_g[node];
    const leaf_p = state.node_leaf_phosphorus_g[node];
    const sheath_c = state.node_sheath_carbon_g[node];
    const sheath_n = state.node_sheath_nitrogen_g[node];
    const sheath_p = state.node_sheath_phosphorus_g[node];
    const recycled_leaf_c = if (allocation.leaf_present) leaf_c * recycling.carbon else 0;
    const recycled_leaf_n = if (allocation.leaf_present) leaf_n * (recycling.nitrogen + (1.0 - recycling.nitrogen) * recycling.carbon) else 0;
    const recycled_leaf_p = if (allocation.leaf_present) leaf_p * (recycling.phosphorus + (1.0 - recycling.phosphorus) * recycling.carbon) else 0;
    const recycled_sheath_c = sheath_c * recycling.carbon;
    const recycled_sheath_n = sheath_n * (recycling.nitrogen + (1.0 - recycling.nitrogen) * recycling.carbon);
    const recycled_sheath_p = sheath_p * (recycling.phosphorus + (1.0 - recycling.phosphorus) * recycling.carbon);
    var products: SenescenceProducts = .{
        .recycled_carbon_g = allocation.carbon_recovered_to_mobile_pool_g,
        .recycled_nitrogen_g = allocation.leaf_fraction * recycled_leaf_n * leaf_woody_nitrogen_fraction[1] + allocation.sheath_fraction * recycled_sheath_n * sheath_woody_nitrogen_fraction[1],
        .recycled_phosphorus_g = allocation.leaf_fraction * recycled_leaf_p * leaf_woody_phosphorus_fraction[1] + allocation.sheath_fraction * recycled_sheath_p * sheath_woody_phosphorus_fraction[1],
        .respired_carbon_g = allocation.carbon_respired_g,
    };
    for (0..4) |kinetic| {
        products.woody_carbon_g[kinetic] = woody_kinetics.carbon[kinetic] * woody_fraction[0] * (allocation.leaf_fraction * leaf_c + allocation.sheath_fraction * sheath_c);
        products.woody_nitrogen_g[kinetic] = woody_kinetics.nitrogen[kinetic] * (allocation.leaf_fraction * leaf_n * leaf_woody_nitrogen_fraction[0] + allocation.sheath_fraction * sheath_n * sheath_woody_nitrogen_fraction[0]);
        products.woody_phosphorus_g[kinetic] = woody_kinetics.phosphorus[kinetic] * (allocation.leaf_fraction * leaf_p * leaf_woody_phosphorus_fraction[0] + allocation.sheath_fraction * sheath_p * sheath_woody_phosphorus_fraction[0]);
        products.nonwoody_carbon_g[kinetic] = woody_fraction[1] * (leaf_kinetics.carbon[kinetic] * allocation.leaf_fraction * (leaf_c - recycled_leaf_c) + sheath_kinetics.carbon[kinetic] * allocation.sheath_fraction * (sheath_c - recycled_sheath_c));
        products.nonwoody_nitrogen_g[kinetic] = leaf_woody_nitrogen_fraction[1] * leaf_kinetics.nitrogen[kinetic] * allocation.leaf_fraction * (leaf_n - recycled_leaf_n) + sheath_woody_nitrogen_fraction[1] * sheath_kinetics.nitrogen[kinetic] * allocation.sheath_fraction * (sheath_n - recycled_sheath_n);
        products.nonwoody_phosphorus_g[kinetic] = leaf_woody_phosphorus_fraction[1] * leaf_kinetics.phosphorus[kinetic] * allocation.leaf_fraction * (leaf_p - recycled_leaf_p) + sheath_woody_phosphorus_fraction[1] * sheath_kinetics.phosphorus[kinetic] * allocation.sheath_fraction * (sheath_p - recycled_sheath_p);
    }
    const leaf_publication = try @import("../leaf/senescence_state_publication.zig").calculate(.{
        .branch = .{ .area_m2 = state.branch_leaf_area_m2[branch], .carbon_g_c = state.branch_leaf_carbon_g[branch], .nitrogen_g_n = state.branch_leaf_nitrogen_g[branch], .phosphorus_g_p = state.branch_leaf_phosphorus_g[branch] },
        .node = .{ .area_m2 = state.node_leaf_area_m2[node], .carbon_g_c = leaf_c, .nitrogen_g_n = leaf_n, .phosphorus_g_p = leaf_p },
        .node_protein_g = state.node_leaf_protein_g[node],
        .senescing_snapshot = .{ .area_m2 = state.node_leaf_area_m2[node], .carbon_g_c = leaf_c, .nitrogen_g_n = leaf_n, .phosphorus_g_p = leaf_p },
        .area_removal_fraction = allocation.leaf_fraction,
        .mass_removal_fraction = allocation.leaf_fraction,
        .protein_per_nitrogen_g_per_g_n = protein_per_nitrogen_g_per_g_n,
        .protein_per_phosphorus_g_per_g_p = protein_per_phosphorus_g_per_g_p,
        .branch_mobile_carbon_g_c = state.branch_mobile_carbon_g[branch],
        .branch_mobile_nitrogen_g_n = state.branch_mobile_nitrogen_g[branch],
        .branch_mobile_phosphorus_g_p = state.branch_mobile_phosphorus_g[branch],
        .recycled_carbon_g_c = products.recycled_carbon_g,
        .recycled_nitrogen_g_n = products.recycled_nitrogen_g,
        .recycled_phosphorus_g_p = products.recycled_phosphorus_g,
    });
    const c4_state: c4_leaf_nonstructural_carbon_senescence.State = .{
        .bundle_sheath_carbon_g_c = state.node_c3_nonstructural_carbon_g,
        .mesophyll_carbon_g_c = state.node_c4_mesophyll_nonstructural_carbon_g,
        .foliar_litter_carbon_g_c_by_kinetic_pool = &products.nonwoody_carbon_g,
    };
    const c4_routing: c4_leaf_nonstructural_carbon_senescence.Routing = .{
        .selected_node = node,
        .foliar_litter_kinetic_pool = 1,
    };
    if (allocation.leaf_present)
        try c4_leaf_nonstructural_carbon_senescence.routePartial(c4_state, c4_routing, allocation.leaf_fraction)
    else
        try c4_leaf_nonstructural_carbon_senescence.routeAll(c4_state, c4_routing);

    state.branch_leaf_area_m2[branch] = leaf_publication.branch.area_m2;
    state.branch_leaf_carbon_g[branch] = leaf_publication.branch.carbon_g_c;
    state.branch_leaf_nitrogen_g[branch] = leaf_publication.branch.nitrogen_g_n;
    state.branch_leaf_phosphorus_g[branch] = leaf_publication.branch.phosphorus_g_p;
    state.node_leaf_area_m2[node] = leaf_publication.node.area_m2;
    state.node_leaf_carbon_g[node] = leaf_publication.node.carbon_g_c;
    state.node_leaf_nitrogen_g[node] = leaf_publication.node.nitrogen_g_n;
    state.node_leaf_phosphorus_g[node] = leaf_publication.node.phosphorus_g_p;
    state.node_leaf_protein_g[node] = leaf_publication.node_protein_g;
    state.node_sheath_height_m[node] *= 1.0 - allocation.sheath_fraction;
    state.node_sheath_carbon_g[node] = sheath_c * (1.0 - allocation.sheath_fraction);
    state.node_sheath_nitrogen_g[node] = sheath_n * (1.0 - allocation.sheath_fraction);
    state.node_sheath_phosphorus_g[node] = sheath_p * (1.0 - allocation.sheath_fraction);
    state.node_sheath_protein_g[node] = @max(0.0, state.node_sheath_protein_g[node] - allocation.sheath_fraction * @max(sheath_n * protein_per_nitrogen_g_per_g_n, sheath_p * protein_per_phosphorus_g_per_g_p));
    state.branch_sheath_carbon_g[branch] = @max(0.0, state.branch_sheath_carbon_g[branch] - allocation.sheath_fraction * sheath_c);
    state.branch_sheath_nitrogen_g[branch] = @max(0.0, state.branch_sheath_nitrogen_g[branch] - allocation.sheath_fraction * sheath_n);
    state.branch_sheath_phosphorus_g[branch] = @max(0.0, state.branch_sheath_phosphorus_g[branch] - allocation.sheath_fraction * sheath_p);
    state.branch_mobile_carbon_g[branch] = leaf_publication.branch_mobile_carbon_g_c;
    state.branch_mobile_nitrogen_g[branch] = leaf_publication.branch_mobile_nitrogen_g_n;
    state.branch_mobile_phosphorus_g[branch] = leaf_publication.branch_mobile_phosphorus_g_p;
    return products;
}

pub const InternodeSenescenceResult = struct { fraction: f64, remaining_respiration_demand_g_c: f64, products: SenescenceProducts };

fn commitInternodeSenescenceDemandScaled(state: *State, branch: usize, node_within_branch: usize, respiration_demand_g_c: f64, phenological_senescence_fraction: f64, scaled_recycling: perennial_stalk_senescence_setup.RecyclingFractions, presence_threshold_g_c: f64, woody_carbon_fraction: [2]f64, woody_nitrogen_fraction: [2]f64, woody_phosphorus_fraction: [2]f64, woody_kinetics: KineticFractions, stalk_kinetics: KineticFractions) !InternodeSenescenceResult {
    const nodes = try state.nodeRange(branch);
    if (node_within_branch >= nodes.end - nodes.first) return error.CanopyNodeIndexOutOfBounds;
    const node = nodes.first + node_within_branch;
    var products: SenescenceProducts = .{};
    const reserve_carbon_before = state.branch_reserve_carbon_g[branch];
    const reserve_nitrogen_before = state.branch_reserve_nitrogen_g[branch];
    const reserve_phosphorus_before = state.branch_reserve_phosphorus_g[branch];
    const result = try internode_senescence_publication.publish(.{
        .branch_stalk_carbon_g_c = &state.branch_stalk_carbon_g[branch],
        .branch_stalk_nitrogen_g_n = &state.branch_stalk_nitrogen_g[branch],
        .branch_stalk_phosphorus_g_p = &state.branch_stalk_phosphorus_g[branch],
        .node_height_m = state.node_height_m,
        .internode_length_m = state.node_internode_length_m,
        .internode_carbon_g_c = state.node_internode_carbon_g,
        .internode_nitrogen_g_n = state.node_internode_nitrogen_g,
        .internode_phosphorus_g_p = state.node_internode_phosphorus_g,
        .reserve_carbon_g_c = &state.branch_reserve_carbon_g[branch],
        .reserve_nitrogen_g_n = &state.branch_reserve_nitrogen_g[branch],
        .reserve_phosphorus_g_p = &state.branch_reserve_phosphorus_g[branch],
        .litter = .{
            .woody_carbon_g_c = &products.woody_carbon_g,
            .woody_nitrogen_g_n = &products.woody_nitrogen_g,
            .woody_phosphorus_g_p = &products.woody_phosphorus_g,
            .stalk_carbon_g_c = &products.nonwoody_carbon_g,
            .stalk_nitrogen_g_n = &products.nonwoody_nitrogen_g,
            .stalk_phosphorus_g_p = &products.nonwoody_phosphorus_g,
        },
    }, .{
        .selected_node = node,
        .respiration_demand_g_c_per_timestep = respiration_demand_g_c,
        .presence_threshold_g_c = presence_threshold_g_c,
        .phenological_senescence_fraction = phenological_senescence_fraction,
        .sapwood_recycling = .{ .carbon = scaled_recycling.carbon, .nitrogen = scaled_recycling.nitrogen, .phosphorus = scaled_recycling.phosphorus },
        .woody_fraction = .{ .carbon = woody_carbon_fraction[0], .nitrogen = woody_nitrogen_fraction[0], .phosphorus = woody_phosphorus_fraction[0] },
        .nonwoody_fraction = .{ .carbon = woody_carbon_fraction[1], .nitrogen = woody_nitrogen_fraction[1], .phosphorus = woody_phosphorus_fraction[1] },
        .woody_kinetics = .{ .carbon = &woody_kinetics.carbon, .nitrogen = &woody_kinetics.nitrogen, .phosphorus = &woody_kinetics.phosphorus },
        .stalk_kinetics = .{ .carbon = &stalk_kinetics.carbon, .nitrogen = &stalk_kinetics.nitrogen, .phosphorus = &stalk_kinetics.phosphorus },
    });
    products.recycled_carbon_g = state.branch_reserve_carbon_g[branch] - reserve_carbon_before;
    products.recycled_nitrogen_g = state.branch_reserve_nitrogen_g[branch] - reserve_nitrogen_before;
    products.recycled_phosphorus_g = state.branch_reserve_phosphorus_g[branch] - reserve_phosphorus_before;
    products.respired_carbon_g = respiration_demand_g_c - result.remaining_respiration_g_c_per_timestep - products.recycled_carbon_g;
    return .{ .fraction = result.removal_fraction, .remaining_respiration_demand_g_c = result.remaining_respiration_g_c_per_timestep, .products = products };
}

pub fn commitInternodeSenescenceDemand(state: *State, branch: usize, node_within_branch: usize, respiration_demand_g_c: f64, phenological_senescence_fraction: f64, recycling: RecyclingFractions, woody_carbon_fraction: [2]f64, woody_nitrogen_fraction: [2]f64, woody_phosphorus_fraction: [2]f64, woody_kinetics: KineticFractions, stalk_kinetics: KineticFractions) !InternodeSenescenceResult {
    inline for (.{ respiration_demand_g_c, phenological_senescence_fraction, recycling.carbon, recycling.nitrogen, recycling.phosphorus }) |value| if (!std.math.isFinite(value)) return error.NonFiniteInternodeSenescenceInput;
    if (respiration_demand_g_c < 0 or phenological_senescence_fraction < 0 or phenological_senescence_fraction > 1 or recycling.carbon < 0 or recycling.carbon > 1 or recycling.nitrogen < 0 or recycling.nitrogen > 1 or recycling.phosphorus < 0 or recycling.phosphorus > 1) return error.InvalidInternodeSenescenceInput;
    inline for (.{ woody_carbon_fraction, woody_nitrogen_fraction, woody_phosphorus_fraction }) |fractions| {
        for (fractions) |fraction| if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidWoodyFraction;
        if (@abs(fractions[0] + fractions[1] - 1.0) > 1e-8) return error.InvalidWoodyFraction;
    }
    try woody_kinetics.validate();
    try stalk_kinetics.validate();
    const nodes = try state.nodeRange(branch);
    if (node_within_branch >= nodes.end - nodes.first) return error.CanopyNodeIndexOutOfBounds;
    const node = nodes.first + node_within_branch;
    const stalk_c = state.branch_stalk_carbon_g[branch];
    if (stalk_c <= 0) return .{ .fraction = 0, .remaining_respiration_demand_g_c = respiration_demand_g_c, .products = .{} };
    const sapwood_fraction = std.math.clamp(state.branch_sapwood_carbon_g[branch] / stalk_c, 0, 1);
    const effective_c_recycling = recycling.carbon * sapwood_fraction;
    const effective_n_recycling = recycling.nitrogen * sapwood_fraction;
    const effective_p_recycling = recycling.phosphorus * sapwood_fraction;
    const node_c = state.node_internode_carbon_g[node];
    const node_n = state.node_internode_nitrogen_g[node];
    const node_p = state.node_internode_phosphorus_g[node];
    if (node_c <= 0) return .{ .fraction = 0, .remaining_respiration_demand_g_c = respiration_demand_g_c, .products = .{} };
    const recycled_c = effective_c_recycling * node_c;
    const recycled_n = node_n * (effective_n_recycling + (1.0 - effective_n_recycling) * effective_c_recycling);
    const recycled_p = node_p * (effective_p_recycling + (1.0 - effective_p_recycling) * effective_c_recycling);
    const fraction = if (recycled_c > 0) std.math.clamp(respiration_demand_g_c / recycled_c, 0, 1) else 1;
    const consumed_recycled_c = fraction * recycled_c * woody_carbon_fraction[1];
    var products: SenescenceProducts = .{
        .recycled_carbon_g = consumed_recycled_c * phenological_senescence_fraction,
        .recycled_nitrogen_g = fraction * recycled_n * woody_nitrogen_fraction[1],
        .recycled_phosphorus_g = fraction * recycled_p * woody_phosphorus_fraction[1],
        .respired_carbon_g = consumed_recycled_c * (1.0 - phenological_senescence_fraction),
    };
    for (0..4) |kinetic| {
        products.woody_carbon_g[kinetic] = woody_kinetics.carbon[kinetic] * fraction * node_c * woody_carbon_fraction[0];
        products.woody_nitrogen_g[kinetic] = woody_kinetics.nitrogen[kinetic] * fraction * node_n * woody_nitrogen_fraction[0];
        products.woody_phosphorus_g[kinetic] = woody_kinetics.phosphorus[kinetic] * fraction * node_p * woody_phosphorus_fraction[0];
        products.nonwoody_carbon_g[kinetic] = stalk_kinetics.carbon[kinetic] * fraction * (node_c - recycled_c) * woody_carbon_fraction[1];
        products.nonwoody_nitrogen_g[kinetic] = stalk_kinetics.nitrogen[kinetic] * fraction * (node_n - recycled_n) * woody_nitrogen_fraction[1];
        products.nonwoody_phosphorus_g[kinetic] = stalk_kinetics.phosphorus[kinetic] * fraction * (node_p - recycled_p) * woody_phosphorus_fraction[1];
    }
    state.branch_stalk_carbon_g[branch] = @max(0.0, state.branch_stalk_carbon_g[branch] - fraction * node_c);
    state.branch_stalk_nitrogen_g[branch] = @max(0.0, state.branch_stalk_nitrogen_g[branch] - fraction * node_n);
    state.branch_stalk_phosphorus_g[branch] = @max(0.0, state.branch_stalk_phosphorus_g[branch] - fraction * node_p);
    state.node_height_m[node] = @max(0.0, state.node_height_m[node] - fraction * state.node_internode_length_m[node]);
    state.node_internode_carbon_g[node] *= 1.0 - fraction;
    state.node_internode_nitrogen_g[node] *= 1.0 - fraction;
    state.node_internode_phosphorus_g[node] *= 1.0 - fraction;
    state.node_internode_length_m[node] *= 1.0 - fraction;
    state.branch_reserve_carbon_g[branch] += products.recycled_carbon_g;
    state.branch_reserve_nitrogen_g[branch] += products.recycled_nitrogen_g;
    state.branch_reserve_phosphorus_g[branch] += products.recycled_phosphorus_g;
    return .{ .fraction = fraction, .remaining_respiration_demand_g_c = @max(0.0, respiration_demand_g_c - consumed_recycled_c), .products = products };
}

pub fn commitResidualStalkSenescenceDemand(state: *State, branch: usize, respiration_demand_g_c: f64, phenological_senescence_fraction: f64, recycling: RecyclingFractions, woody_carbon_fraction: [2]f64, woody_nitrogen_fraction: [2]f64, woody_phosphorus_fraction: [2]f64, woody_kinetics: KineticFractions, stalk_kinetics: KineticFractions) !InternodeSenescenceResult {
    inline for (.{ respiration_demand_g_c, phenological_senescence_fraction, recycling.carbon, recycling.nitrogen, recycling.phosphorus }) |value| if (!std.math.isFinite(value)) return error.NonFiniteResidualStalkSenescenceInput;
    if (branch >= state.branch_stalk_carbon_g.len) return error.CanopyBranchIndexOutOfBounds;
    if (respiration_demand_g_c < 0 or phenological_senescence_fraction < 0 or phenological_senescence_fraction > 1) return error.InvalidResidualStalkSenescenceInput;
    inline for (.{ woody_carbon_fraction, woody_nitrogen_fraction, woody_phosphorus_fraction }) |fractions| {
        for (fractions) |fraction_value| if (!std.math.isFinite(fraction_value) or fraction_value < 0 or fraction_value > 1) return error.InvalidWoodyFraction;
        if (@abs(fractions[0] + fractions[1] - 1.0) > 1e-8) return error.InvalidWoodyFraction;
    }
    try woody_kinetics.validate();
    try stalk_kinetics.validate();
    const stalk_c = state.branch_stalk_carbon_g[branch];
    const residual_c = state.branch_senescing_stalk_carbon_g[branch];
    const residual_n = state.branch_senescing_stalk_nitrogen_g[branch];
    const residual_p = state.branch_senescing_stalk_phosphorus_g[branch];
    if (stalk_c <= 0 or residual_c <= 0) return .{ .fraction = 0, .remaining_respiration_demand_g_c = respiration_demand_g_c, .products = .{} };
    const sapwood_fraction = std.math.clamp(state.branch_sapwood_carbon_g[branch] / stalk_c, 0, 1);
    const effective_c = recycling.carbon * sapwood_fraction;
    const effective_n = recycling.nitrogen * sapwood_fraction;
    const effective_p = recycling.phosphorus * sapwood_fraction;
    const recyclable_c = effective_c * residual_c;
    const recyclable_n = residual_n * (effective_n + (1.0 - effective_n) * effective_c);
    const recyclable_p = residual_p * (effective_p + (1.0 - effective_p) * effective_c);
    const fraction = if (recyclable_c > 0) std.math.clamp(respiration_demand_g_c / recyclable_c, 0, 1) else 1;
    const consumed_recyclable_c = fraction * recyclable_c * woody_carbon_fraction[1];
    var products: SenescenceProducts = .{
        .recycled_carbon_g = consumed_recyclable_c * phenological_senescence_fraction,
        .recycled_nitrogen_g = fraction * recyclable_n * woody_nitrogen_fraction[1],
        .recycled_phosphorus_g = fraction * recyclable_p * woody_phosphorus_fraction[1],
        .respired_carbon_g = consumed_recyclable_c * (1.0 - phenological_senescence_fraction),
    };
    for (0..4) |kinetic| {
        products.woody_carbon_g[kinetic] = woody_kinetics.carbon[kinetic] * fraction * residual_c * woody_carbon_fraction[0];
        products.woody_nitrogen_g[kinetic] = woody_kinetics.nitrogen[kinetic] * fraction * residual_n * woody_nitrogen_fraction[0];
        products.woody_phosphorus_g[kinetic] = woody_kinetics.phosphorus[kinetic] * fraction * residual_p * woody_phosphorus_fraction[0];
        products.nonwoody_carbon_g[kinetic] = stalk_kinetics.carbon[kinetic] * fraction * (residual_c - recyclable_c) * woody_carbon_fraction[1];
        products.nonwoody_nitrogen_g[kinetic] = stalk_kinetics.nitrogen[kinetic] * fraction * (residual_n - recyclable_n) * woody_nitrogen_fraction[1];
        products.nonwoody_phosphorus_g[kinetic] = stalk_kinetics.phosphorus[kinetic] * fraction * (residual_p - recyclable_p) * woody_phosphorus_fraction[1];
    }
    state.branch_stalk_carbon_g[branch] = @max(0.0, stalk_c - fraction * residual_c);
    state.branch_stalk_nitrogen_g[branch] = @max(0.0, state.branch_stalk_nitrogen_g[branch] - fraction * residual_n);
    state.branch_stalk_phosphorus_g[branch] = @max(0.0, state.branch_stalk_phosphorus_g[branch] - fraction * residual_p);
    state.branch_senescing_stalk_carbon_g[branch] *= 1.0 - fraction;
    state.branch_senescing_stalk_nitrogen_g[branch] *= 1.0 - fraction;
    state.branch_senescing_stalk_phosphorus_g[branch] *= 1.0 - fraction;
    const nodes = try state.nodeRange(branch);
    var maximum_height_m: f64 = 0;
    for (state.node_height_m[nodes.first..nodes.end]) |height_m| maximum_height_m = @max(maximum_height_m, height_m);
    const reduced_maximum_height_m = maximum_height_m * (1.0 - fraction);
    for (state.node_height_m[nodes.first..nodes.end]) |*height_m| height_m.* = @min(height_m.*, reduced_maximum_height_m);
    state.branch_reserve_carbon_g[branch] += products.recycled_carbon_g;
    state.branch_reserve_nitrogen_g[branch] += products.recycled_nitrogen_g;
    state.branch_reserve_phosphorus_g[branch] += products.recycled_phosphorus_g;
    return .{ .fraction = fraction, .remaining_respiration_demand_g_c = @max(0.0, respiration_demand_g_c - consumed_recyclable_c), .products = products };
}

fn commitResidualStalkSenescenceDemandScaled(state: *State, branch: usize, respiration_demand_g_c: f64, phenological_senescence_fraction: f64, scaled_recycling: perennial_stalk_senescence_setup.RecyclingFractions, presence_threshold_g_c: f64, woody_carbon_fraction: [2]f64, woody_nitrogen_fraction: [2]f64, woody_phosphorus_fraction: [2]f64, woody_kinetics: KineticFractions, stalk_kinetics: KineticFractions) !InternodeSenescenceResult {
    const request = try residual_stalk_senescence_request.calculate(.{
        .carbon = state.branch_senescing_stalk_carbon_g[branch],
        .nitrogen = state.branch_senescing_stalk_nitrogen_g[branch],
        .phosphorus = state.branch_senescing_stalk_phosphorus_g[branch],
    }, .{
        .carbon = scaled_recycling.carbon,
        .nitrogen = scaled_recycling.nitrogen,
        .phosphorus = scaled_recycling.phosphorus,
    }, respiration_demand_g_c, presence_threshold_g_c) orelse return .{
        .fraction = 0,
        .remaining_respiration_demand_g_c = respiration_demand_g_c,
        .products = .{},
    };
    const nodes = try state.nodeRange(branch);
    var products: SenescenceProducts = .{};
    const reserve_carbon_before = state.branch_reserve_carbon_g[branch];
    const reserve_nitrogen_before = state.branch_reserve_nitrogen_g[branch];
    const reserve_phosphorus_before = state.branch_reserve_phosphorus_g[branch];
    const result = try residual_stalk_senescence_publication.publish(.{
        .branch_stalk_carbon_g_c = &state.branch_stalk_carbon_g[branch],
        .branch_stalk_nitrogen_g_n = &state.branch_stalk_nitrogen_g[branch],
        .branch_stalk_phosphorus_g_p = &state.branch_stalk_phosphorus_g[branch],
        .residual_stalk_carbon_g_c = &state.branch_senescing_stalk_carbon_g[branch],
        .residual_stalk_nitrogen_g_n = &state.branch_senescing_stalk_nitrogen_g[branch],
        .residual_stalk_phosphorus_g_p = &state.branch_senescing_stalk_phosphorus_g[branch],
        .node_height_m = state.node_height_m[nodes.first..nodes.end],
        .reserve_carbon_g_c = &state.branch_reserve_carbon_g[branch],
        .reserve_nitrogen_g_n = &state.branch_reserve_nitrogen_g[branch],
        .reserve_phosphorus_g_p = &state.branch_reserve_phosphorus_g[branch],
        .litter = .{
            .woody_carbon_g_c = &products.woody_carbon_g,
            .woody_nitrogen_g_n = &products.woody_nitrogen_g,
            .woody_phosphorus_g_p = &products.woody_phosphorus_g,
            .stalk_carbon_g_c = &products.nonwoody_carbon_g,
            .stalk_nitrogen_g_n = &products.nonwoody_nitrogen_g,
            .stalk_phosphorus_g_p = &products.nonwoody_phosphorus_g,
        },
    }, .{
        .removal_fraction = request.removal_fraction,
        .recyclable = .{ .carbon = request.recyclable.carbon, .nitrogen = request.recyclable.nitrogen, .phosphorus = request.recyclable.phosphorus },
        .respiration_demand_g_c_per_timestep = respiration_demand_g_c,
        .phenological_senescence_fraction = phenological_senescence_fraction,
        .woody_fraction = .{ .carbon = woody_carbon_fraction[0], .nitrogen = woody_nitrogen_fraction[0], .phosphorus = woody_phosphorus_fraction[0] },
        .nonwoody_fraction = .{ .carbon = woody_carbon_fraction[1], .nitrogen = woody_nitrogen_fraction[1], .phosphorus = woody_phosphorus_fraction[1] },
        .woody_kinetics = .{ .carbon = &woody_kinetics.carbon, .nitrogen = &woody_kinetics.nitrogen, .phosphorus = &woody_kinetics.phosphorus },
        .stalk_kinetics = .{ .carbon = &stalk_kinetics.carbon, .nitrogen = &stalk_kinetics.nitrogen, .phosphorus = &stalk_kinetics.phosphorus },
    });
    products.recycled_carbon_g = state.branch_reserve_carbon_g[branch] - reserve_carbon_before;
    products.recycled_nitrogen_g = state.branch_reserve_nitrogen_g[branch] - reserve_nitrogen_before;
    products.recycled_phosphorus_g = state.branch_reserve_phosphorus_g[branch] - reserve_phosphorus_before;
    products.respired_carbon_g = respiration_demand_g_c - result.remaining_respiration_g_c_per_timestep - products.recycled_carbon_g;
    return .{ .fraction = request.removal_fraction, .remaining_respiration_demand_g_c = result.remaining_respiration_g_c_per_timestep, .products = products };
}

pub const ReserveFallbackPolicy = enum { source_compatible, consume_available };

pub const SenescenceLitterParameters = struct {
    woody_carbon_fraction: [2]f64,
    leaf_woody_nitrogen_fraction: [2]f64,
    sheath_woody_nitrogen_fraction: [2]f64,
    stalk_woody_nitrogen_fraction: [2]f64,
    leaf_woody_phosphorus_fraction: [2]f64,
    sheath_woody_phosphorus_fraction: [2]f64,
    stalk_woody_phosphorus_fraction: [2]f64,
    woody_kinetics: KineticFractions,
    leaf_kinetics: KineticFractions,
    sheath_kinetics: KineticFractions,
    stalk_kinetics: KineticFractions,
};

pub const BranchSenescenceRequest = struct {
    total_respiration_demand_g_c: f64,
    phenological_senescence_fraction: f64,
    first_node_within_branch: usize,
    last_node_within_branch: usize,
    node_group_count: usize,
    perennial: bool,
    reserve_fallback_policy: ReserveFallbackPolicy = .source_compatible,
    leaf_presence_threshold_g_c: f64 = 0,
    demand_tolerance_g_c: f64,
};

pub const BranchSenescenceResult = struct {
    remaining_respiration_demand_g_c: f64,
    reserve_carbon_respired_g_c: f64,
    products: SenescenceProducts,
};

/// Executes the GROSUB senescence cascade from progressively older leaf-node
/// groups through reserve C, internodes, and residual stalk. Runtime node
/// offsets replace the source's modulo-25 storage ring.
pub fn commitBranchSenescenceDemand(state: *State, branch: usize, request: BranchSenescenceRequest, recycling: RecyclingFractions, protein_per_nitrogen_g_per_g_n: f64, protein_per_phosphorus_g_per_g_p: f64, litter: SenescenceLitterParameters) !BranchSenescenceResult {
    inline for (.{ request.total_respiration_demand_g_c, request.phenological_senescence_fraction, request.leaf_presence_threshold_g_c, request.demand_tolerance_g_c }) |value| if (!std.math.isFinite(value)) return error.NonFiniteBranchSenescenceInput;
    if (request.total_respiration_demand_g_c < 0 or request.phenological_senescence_fraction < 0 or request.phenological_senescence_fraction > 1 or request.node_group_count == 0 or request.leaf_presence_threshold_g_c < 0 or request.demand_tolerance_g_c < 0) return error.InvalidBranchSenescenceInput;
    const nodes = try state.nodeRange(branch);
    const node_count = nodes.end - nodes.first;
    if (request.first_node_within_branch > request.last_node_within_branch or request.last_node_within_branch >= node_count) return error.CanopyNodeIndexOutOfBounds;
    var result: BranchSenescenceResult = .{ .remaining_respiration_demand_g_c = 0, .reserve_carbon_respired_g_c = 0, .products = .{} };
    const cascade_threshold_g_c = if (request.reserve_fallback_policy == .source_compatible)
        request.leaf_presence_threshold_g_c
    else
        request.demand_tolerance_g_c;
    const group_demand_g_c = try node_senescence_remobilization_request.respirationPerPass(
        request.total_respiration_demand_g_c,
        request.node_group_count,
    );
    for (0..request.node_group_count) |group| {
        var remaining_g_c = group_demand_g_c;
        const group_start = try node_senescence_cascade_progress.firstNodeForPass(request.first_node_within_branch, request.last_node_within_branch, group);
        if (group_start) |first_node| for (first_node..request.last_node_within_branch + 1) |node_within_branch| {
            const node = nodes.first + node_within_branch;
            const allocation = try allocateNodeSenescenceDemandWithThreshold(remaining_g_c, state.node_leaf_carbon_g[node], state.node_sheath_carbon_g[node], recycling.carbon, request.phenological_senescence_fraction, litter.woody_carbon_fraction[1], request.leaf_presence_threshold_g_c);
            const products = try commitNodeSenescenceDemand(state, branch, node_within_branch, allocation, recycling, protein_per_nitrogen_g_per_g_n, protein_per_phosphorus_g_per_g_p, litter.woody_carbon_fraction, litter.leaf_woody_nitrogen_fraction, litter.sheath_woody_nitrogen_fraction, litter.leaf_woody_phosphorus_fraction, litter.sheath_woody_phosphorus_fraction, litter.woody_kinetics, litter.leaf_kinetics, litter.sheath_kinetics);
            addSenescenceProducts(&result.products, products);
            remaining_g_c = allocation.remaining_respiration_demand_g_c;
            if (remaining_g_c <= cascade_threshold_g_c) break;
        };
        if (remaining_g_c > cascade_threshold_g_c) {
            if (request.reserve_fallback_policy == .consume_available) {
                const reserve_g_c = state.branch_reserve_carbon_g[branch];
                const consumed_g_c = @min(reserve_g_c, remaining_g_c);
                state.branch_reserve_carbon_g[branch] -= consumed_g_c;
                result.reserve_carbon_respired_g_c += consumed_g_c;
                remaining_g_c -= consumed_g_c;
            } else {
                const reserve = try node_senescence_cascade_progress.applyReserveFallback(.{
                    .reserve_carbon_g_c = &state.branch_reserve_carbon_g[branch],
                }, remaining_g_c, request.phenological_senescence_fraction);
                result.reserve_carbon_respired_g_c += reserve.reserve_carbon_respired_g_c_per_timestep;
                remaining_g_c = reserve.excess_maintenance_respiration_g_c_per_timestep;
            }
        }
        if (request.perennial and remaining_g_c > cascade_threshold_g_c) {
            const stalk_setup = try perennial_stalk_senescence_setup.prepare(.{
                .is_perennial = true,
                .excess_maintenance_respiration_g_c_per_timestep = remaining_g_c,
                .stalk_carbon_g_c = state.branch_stalk_carbon_g[branch],
                .sapwood_carbon_g_c = state.branch_sapwood_carbon_g[branch],
                .presence_threshold_g_c = request.leaf_presence_threshold_g_c,
                .first_internode = request.first_node_within_branch,
                .last_internode = request.last_node_within_branch,
                .shoot_recycling = .{ .carbon = recycling.carbon, .nitrogen = recycling.nitrogen, .phosphorus = recycling.phosphorus },
            });
            if (stalk_setup) |setup| for (0..setup.last_internode - setup.first_internode + 1) |iteration| {
                const ordinal = perennial_stalk_senescence_setup.descendingInternode(setup, iteration).?;
                const internode = try commitInternodeSenescenceDemandScaled(state, branch, ordinal, remaining_g_c, setup.phenological_senescence_fraction, setup.sapwood_recycling, request.leaf_presence_threshold_g_c, litter.woody_carbon_fraction, litter.stalk_woody_nitrogen_fraction, litter.stalk_woody_phosphorus_fraction, litter.woody_kinetics, litter.stalk_kinetics);
                addSenescenceProducts(&result.products, internode.products);
                remaining_g_c = internode.remaining_respiration_demand_g_c;
                if (remaining_g_c <= cascade_threshold_g_c) break;
            };
            if (stalk_setup != null and remaining_g_c > cascade_threshold_g_c) {
                const setup = stalk_setup.?;
                const residual = try commitResidualStalkSenescenceDemandScaled(state, branch, remaining_g_c, setup.phenological_senescence_fraction, setup.sapwood_recycling, request.leaf_presence_threshold_g_c, litter.woody_carbon_fraction, litter.stalk_woody_nitrogen_fraction, litter.stalk_woody_phosphorus_fraction, litter.woody_kinetics, litter.stalk_kinetics);
                addSenescenceProducts(&result.products, residual.products);
                remaining_g_c = residual.remaining_respiration_demand_g_c;
            }
        }
        result.remaining_respiration_demand_g_c += @max(0.0, remaining_g_c);
    }
    return result;
}

pub const GrainFillResult = struct {
    carbon_translocated_g: f64,
    nitrogen_translocated_g: f64,
    phosphorus_translocated_g: f64,
    maximum_carbon_translocation_g: f64,
};

pub const LeafNutrientRemobilization = struct { nitrogen_g: f64, phosphorus_g: f64 };

/// GROSUB leaf structural nutrient equilibration. The 1e-3 exchange
/// coefficient and coupled 10:1 N:P bounds are retained from the source.
pub fn remobilizeNodeLeafNutrients(state: *State, branch: usize, node_within_branch: usize, exchange_fraction: f64, minimum_leaf_nutrient_fraction: f64, maximum_leaf_nitrogen_per_carbon_g_n_per_g_c: f64, maximum_leaf_phosphorus_per_carbon_g_p_per_g_c: f64, protein_per_nitrogen_g_per_g_n: f64, protein_per_phosphorus_g_per_g_p: f64) !LeafNutrientRemobilization {
    inline for (.{ exchange_fraction, minimum_leaf_nutrient_fraction, maximum_leaf_nitrogen_per_carbon_g_n_per_g_c, maximum_leaf_phosphorus_per_carbon_g_p_per_g_c, protein_per_nitrogen_g_per_g_n, protein_per_phosphorus_g_per_g_p }) |value| if (!std.math.isFinite(value)) return error.NonFiniteLeafNutrientRemobilizationInput;
    if (exchange_fraction < 0 or minimum_leaf_nutrient_fraction < 0 or minimum_leaf_nutrient_fraction > 1 or maximum_leaf_nitrogen_per_carbon_g_n_per_g_c < 0 or maximum_leaf_phosphorus_per_carbon_g_p_per_g_c < 0 or protein_per_nitrogen_g_per_g_n < 0 or protein_per_phosphorus_g_per_g_p < 0) return error.InvalidLeafNutrientRemobilizationInput;
    const nodes = try state.nodeRange(branch);
    if (node_within_branch >= nodes.end - nodes.first) return error.CanopyNodeIndexOutOfBounds;
    const node = nodes.first + node_within_branch;
    const leaf_c = state.node_leaf_carbon_g[node];
    const leaf_n = state.node_leaf_nitrogen_g[node];
    const leaf_p = state.node_leaf_phosphorus_g[node];
    if (leaf_c <= 0) return .{ .nitrogen_g = 0, .phosphorus_g = 0 };
    const total_carbon_g = leaf_c + state.branch_mobile_carbon_g[branch];
    if (total_carbon_g <= 0) return .{ .nitrogen_g = 0, .phosphorus_g = 0 };
    const nitrogen_gradient_g2 = leaf_n * state.branch_mobile_carbon_g[branch] - state.branch_mobile_nitrogen_g[branch] * leaf_c;
    const phosphorus_gradient_g2 = leaf_p * state.branch_mobile_carbon_g[branch] - state.branch_mobile_phosphorus_g[branch] * leaf_c;
    const unconstrained_n = @max(0.0, exchange_fraction * nitrogen_gradient_g2 / total_carbon_g);
    const unconstrained_p = @max(0.0, exchange_fraction * phosphorus_gradient_g2 / total_carbon_g);
    const removable_n = @max(0.0, leaf_n - minimum_leaf_nutrient_fraction * maximum_leaf_nitrogen_per_carbon_g_n_per_g_c * leaf_c);
    const removable_p = @max(0.0, leaf_p - minimum_leaf_nutrient_fraction * maximum_leaf_phosphorus_per_carbon_g_p_per_g_c * leaf_c);
    const base_n = @min(unconstrained_n, removable_n);
    const base_p = @min(unconstrained_p, removable_p);
    const nitrogen_g = @min(leaf_n, @max(base_n, 10.0 * base_p));
    const phosphorus_g = @min(leaf_p, @max(base_p, 0.1 * base_n));
    state.node_leaf_nitrogen_g[node] -= nitrogen_g;
    state.branch_leaf_nitrogen_g[branch] = @max(0.0, state.branch_leaf_nitrogen_g[branch] - nitrogen_g);
    state.branch_mobile_nitrogen_g[branch] += nitrogen_g;
    state.node_leaf_phosphorus_g[node] -= phosphorus_g;
    state.branch_leaf_phosphorus_g[branch] = @max(0.0, state.branch_leaf_phosphorus_g[branch] - phosphorus_g);
    state.branch_mobile_phosphorus_g[branch] += phosphorus_g;
    state.node_leaf_protein_g[node] = @max(0.0, state.node_leaf_protein_g[node] - @max(nitrogen_g * protein_per_nitrogen_g_per_g_n, phosphorus_g * protein_per_phosphorus_g_per_g_p));
    return .{ .nitrogen_g = nitrogen_g, .phosphorus_g = phosphorus_g };
}

/// GROSUB IFLGZ=0 reserve oxidation before leaf/stalk senescence. When shoot
/// remobilization is already enabled (IFLGZ=1), the later senescence cascade
/// owns the reserve transaction instead.
pub fn consumeReserveForRespiration(state: *State, branch: usize, shoot_remobilization_enabled: bool, remaining_respiration_demand_g_c: f64, maximum_nonstructural_carbon_oxidation_per_h: f64, canopy_growth_temperature_factor: f64, timestep_h: f64) !f64 {
    if (branch >= state.branch_reserve_carbon_g.len) return error.CanopyBranchIndexOutOfBounds;
    var demand = [1]f64{remaining_respiration_demand_g_c};
    _ = try reserve_maintenance_respiration.apply(.{
        .reserve_carbon_g_c = state.branch_reserve_carbon_g[branch .. branch + 1],
        .excess_maintenance_respiration_g_c_per_timestep = &demand,
    }, .{
        .branch = 0,
        .shoot_remobilization_status = if (shoot_remobilization_enabled) .active else .not_started,
        .maximum_nonstructural_carbon_oxidation_per_h = maximum_nonstructural_carbon_oxidation_per_h,
        .canopy_growth_temperature_response = canopy_growth_temperature_factor,
        .timestep_h = timestep_h,
    });
    return demand[0];
}

pub const SenescenceDemand = struct {
    phenological_respiration_g_c: f64,
    total_respiration_g_c: f64,
    phenological_fraction: f64,
    node_group_count: usize,
    first_preceding_node: usize,
};

pub fn senescenceDemand(shoot_remobilization_enabled: bool, phenological_remobilization_enabled: bool, perennial: bool, canopy_leaf_area_m2: f64, horizontal_cell_area_m2: f64, leaf_storage_exchange_per_h: f64, branch_leaf_sheath_carbon_g: f64, remobilization_elapsed_h: f64, full_senescence_h: f64, timestep_h: f64, excess_maintenance_respiration_g_c: f64, highest_leaf_ordinal: usize, lowest_leaf_ordinal: usize) !SenescenceDemand {
    const exchange = [1]f64{leaf_storage_exchange_per_h};
    const setup = try shoot_total_senescence_setup.calculate(.{
        .shoot_remobilization = if (shoot_remobilization_enabled) .enabled else .disabled,
        .phenological_remobilization = if (phenological_remobilization_enabled) .enabled else .disabled,
        .perennial_growth_habit = perennial,
        .plant_leaf_area_m2 = canopy_leaf_area_m2,
        .horizontal_cell_area_m2 = horizontal_cell_area_m2,
        .aboveground_turnover_index = 0,
        .leaf_storage_exchange_fraction_per_h_by_turnover = &exchange,
        .branch_leaf_and_sheath_carbon_g_c = branch_leaf_sheath_carbon_g,
        .remobilization_elapsed_h = remobilization_elapsed_h,
        .full_senescence_duration_h = full_senescence_h,
        .timestep_h = timestep_h,
        .excess_maintenance_respiration_g_c_per_timestep = excess_maintenance_respiration_g_c,
        .structural_presence_threshold_g_c = 0,
        .newest_leaf_node = highest_leaf_ordinal,
        .lowest_leaf_node = lowest_leaf_ordinal,
        .runtime_node_count = try std.math.add(usize, highest_leaf_ordinal, 1),
    });
    if (setup == null) return .{
        .phenological_respiration_g_c = 0,
        .total_respiration_g_c = 0,
        .phenological_fraction = 0,
        .node_group_count = (highest_leaf_ordinal - lowest_leaf_ordinal) / 2 + 1,
        .first_preceding_node = lowest_leaf_ordinal -| 1,
    };
    const result = setup.?;
    return .{
        .phenological_respiration_g_c = result.phenological_respiration_g_c_per_timestep,
        .total_respiration_g_c = result.total_senescence_respiration_g_c_per_timestep,
        .phenological_fraction = result.phenological_fraction,
        .node_group_count = result.node_group_count,
        .first_preceding_node = result.first_preceding_node,
    };
}

pub const ReserveExchange = struct { carbon_g: f64, nitrogen_g: f64, phosphorus_g: f64 };

pub const ElementalMass = struct { carbon_g: f64 = 0, nitrogen_g: f64 = 0, phosphorus_g: f64 = 0 };

pub const HarvestProducts = struct { ecosystem_export: ElementalMass = .{}, litter: ElementalMass = .{} };

const HarvestPartition = struct { remaining: ElementalMass, products: HarvestProducts };

pub const SourceOrderMobileRemoval = struct {
    remaining: ElementalMass,
    unclamped_carbon_retention_fraction: f64,
};

/// Exact GROSUB 9194-9217 non-grazing branch mobile retention ratio.
pub fn sourceOrderNonGrazingMobileRetention(
    previous_leaf_sheath_carbon_g_c: f64,
    remaining_leaf_sheath_carbon_g_c: f64,
    plant_presence_threshold_g_c: f64,
) !f64 {
    inline for (.{
        previous_leaf_sheath_carbon_g_c,
        remaining_leaf_sheath_carbon_g_c,
        plant_presence_threshold_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidHarvestMass;
    if (remaining_leaf_sheath_carbon_g_c > previous_leaf_sheath_carbon_g_c)
        return error.InvalidHarvestMass;
    if (previous_leaf_sheath_carbon_g_c <= plant_presence_threshold_g_c) return 0;
    return @max(0, @min(1, remaining_leaf_sheath_carbon_g_c / previous_leaf_sheath_carbon_g_c));
}

/// Exact GROSUB 9221-9238 proportional host/symbiont mobile removal. The
/// unclamped fraction is exposed because source applies it to symbiont
/// structural pools even when requested C exceeds mobile C.
pub fn sourceOrderProportionalMobileRemoval(
    initial: ElementalMass,
    requested_carbon_g_c: f64,
    plant_presence_threshold_g_c: f64,
) !SourceOrderMobileRemoval {
    inline for (.{ initial.carbon_g, initial.nitrogen_g, initial.phosphorus_g, requested_carbon_g_c, plant_presence_threshold_g_c }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidHarvestMass;
    if (initial.carbon_g <= plant_presence_threshold_g_c) return .{
        .remaining = .{},
        .unclamped_carbon_retention_fraction = 0,
    };
    const ratio = 1 - requested_carbon_g_c / initial.carbon_g;
    return .{
        .remaining = .{
            .carbon_g = @max(0, initial.carbon_g - requested_carbon_g_c),
            .nitrogen_g = @max(0, initial.nitrogen_g - requested_carbon_g_c * initial.nitrogen_g / initial.carbon_g),
            .phosphorus_g = @max(0, initial.phosphorus_g - requested_carbon_g_c * initial.phosphorus_g / initial.carbon_g),
        },
        .unclamped_carbon_retention_fraction = ratio,
    };
}

/// Exact GROSUB 9303-9315 C4-intermediate retention selector.
pub fn sourceOrderC4IntermediateRetention(
    is_c4: bool,
    initial_host_mobile_carbon_g_c: f64,
    remaining_host_mobile_carbon_g_c: f64,
    plant_presence_threshold_g_c: f64,
) !f64 {
    inline for (.{ initial_host_mobile_carbon_g_c, remaining_host_mobile_carbon_g_c, plant_presence_threshold_g_c }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidHarvestMass;
    if (remaining_host_mobile_carbon_g_c > initial_host_mobile_carbon_g_c)
        return error.InvalidHarvestMass;
    if (!is_c4 or initial_host_mobile_carbon_g_c <= plant_presence_threshold_g_c) return 1;
    return remaining_host_mobile_carbon_g_c / initial_host_mobile_carbon_g_c;
}

pub const GrazingPools = struct {
    leaf_carbon_g: f64,
    sheath_carbon_g: f64,
    husk_carbon_g: f64,
    ear_carbon_g: f64,
    grain_carbon_g: f64,
    stalk_carbon_g: f64,
    reserve_carbon_g: f64,
};

pub const GrazingAllocation = struct {
    structural_leaf_carbon_g: f64 = 0,
    structural_sheath_carbon_g: f64 = 0,
    husk_carbon_g: f64 = 0,
    ear_carbon_g: f64 = 0,
    grain_carbon_g: f64 = 0,
    stalk_carbon_g: f64 = 0,
    reserve_carbon_g: f64 = 0,
    mobile_carbon_g: f64 = 0,
    symbiont_mobile_carbon_g: f64 = 0,
    unmet_carbon_g: f64 = 0,
};

pub const SourceOrderAdditionalGrazingRemoval = struct {
    next_total_removed_carbon_g_c: f64,
    next_unmet_carbon_g_c: f64,
};

pub const SourceOrderGrazingNodeRemoval = struct {
    remaining_fraction: f64,
    remaining_branch_layer_demand_g_c: f64,
};

/// Exact GROSUB 8864-8897 branch-layer demand and node retention selector.
pub fn sourceOrderGrazingNodeRemoval(
    node_layer_carbon_g_c: f64,
    branch_layer_demand_g_c: f64,
) !SourceOrderGrazingNodeRemoval {
    inline for (.{ node_layer_carbon_g_c, branch_layer_demand_g_c }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidGrazingPool;
    if (branch_layer_demand_g_c == 0) return .{
        .remaining_fraction = 1,
        .remaining_branch_layer_demand_g_c = 0,
    };
    const remaining_fraction = if (node_layer_carbon_g_c > branch_layer_demand_g_c)
        @max(0, @min(1, (node_layer_carbon_g_c - branch_layer_demand_g_c) / node_layer_carbon_g_c))
    else
        1;
    return .{
        .remaining_fraction = remaining_fraction,
        .remaining_branch_layer_demand_g_c = branch_layer_demand_g_c -
            (1 - remaining_fraction) * node_layer_carbon_g_c,
    };
}

/// Exact GROSUB 8864-8870 plant-to-branch-layer grazing allocation.
pub fn sourceOrderBranchLayerLeafDemand(
    plant_leaf_carbon_g_c: f64,
    plant_structural_leaf_removal_g_c: f64,
    branch_layer_leaf_carbon_g_c: f64,
    plant_leaf_presence_threshold_g_c: f64,
) !f64 {
    inline for (.{
        plant_leaf_carbon_g_c,
        plant_structural_leaf_removal_g_c,
        branch_layer_leaf_carbon_g_c,
        plant_leaf_presence_threshold_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidGrazingPool;
    if (plant_leaf_carbon_g_c <= plant_leaf_presence_threshold_g_c) return 0;
    const demand = plant_structural_leaf_removal_g_c *
        @max(0, branch_layer_leaf_carbon_g_c) / plant_leaf_carbon_g_c;
    if (!std.math.isFinite(demand)) return error.NonFiniteGrazingAllocationInput;
    return demand;
}

/// Exact GROSUB 8758-8783 redistribution operand for one nonfoliar organ.
/// Source caps the additional removal by the original pool, not its remainder.
pub fn sourceOrderAdditionalGrazingRemoval(
    organ_pool_carbon_g_c: f64,
    already_removed_carbon_g_c: f64,
    requested_additional_carbon_g_c: f64,
) !SourceOrderAdditionalGrazingRemoval {
    inline for (.{ organ_pool_carbon_g_c, already_removed_carbon_g_c, requested_additional_carbon_g_c }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidGrazingPool;
    const additional_removed_g_c = @min(organ_pool_carbon_g_c, requested_additional_carbon_g_c);
    return .{
        .next_total_removed_carbon_g_c = already_removed_carbon_g_c + additional_removed_g_c,
        .next_unmet_carbon_g_c = @max(0, requested_additional_carbon_g_c - additional_removed_g_c),
    };
}

pub fn allocateGrazingDemand(total_demand_g_c: f64, harvested_leaf_fraction: f64, harvested_nonfoliar_fraction: f64, harvested_woody_fraction: f64, canopy_mobile_carbon_concentration_g_per_g: f64, symbiont_mobile_carbon_concentration_g_per_g: f64, pools: GrazingPools) !GrazingAllocation {
    inline for (@typeInfo(GrazingPools).@"struct".fields) |field| if (!std.math.isFinite(@field(pools, field.name)) or @field(pools, field.name) < 0) return error.InvalidGrazingPool;
    inline for (.{ total_demand_g_c, harvested_leaf_fraction, harvested_nonfoliar_fraction, harvested_woody_fraction, canopy_mobile_carbon_concentration_g_per_g, symbiont_mobile_carbon_concentration_g_per_g }) |value| if (!std.math.isFinite(value)) return error.NonFiniteGrazingAllocationInput;
    if (total_demand_g_c < 0 or harvested_leaf_fraction < 0 or harvested_leaf_fraction > 1 or harvested_nonfoliar_fraction < 0 or harvested_nonfoliar_fraction > 1 or harvested_woody_fraction < 0 or harvested_woody_fraction > 1 or canopy_mobile_carbon_concentration_g_per_g < 0 or symbiont_mobile_carbon_concentration_g_per_g < 0) return error.InvalidGrazingAllocationInput;
    const mobile_fraction = canopy_mobile_carbon_concentration_g_per_g / (1.0 + canopy_mobile_carbon_concentration_g_per_g);
    const symbiont_fraction = symbiont_mobile_carbon_concentration_g_per_g / (1.0 + symbiont_mobile_carbon_concentration_g_per_g);
    var result: GrazingAllocation = .{};
    const requested_leaf = total_demand_g_c * harvested_leaf_fraction;
    const removed_leaf = @min(pools.leaf_carbon_g, requested_leaf);
    result.structural_leaf_carbon_g = removed_leaf * (1.0 - mobile_fraction);
    result.mobile_carbon_g = removed_leaf * mobile_fraction;
    result.symbiont_mobile_carbon_g = removed_leaf * symbiont_fraction;
    var unmet = @max(0.0, requested_leaf - removed_leaf);
    var removed_sheath_total: f64 = 0;

    const nonfoliar_total = pools.sheath_carbon_g + pools.husk_carbon_g + pools.ear_carbon_g + pools.grain_carbon_g;
    const requested_nonfoliar = total_demand_g_c * harvested_nonfoliar_fraction;
    if (nonfoliar_total > 0) {
        var request = requested_nonfoliar * pools.sheath_carbon_g / nonfoliar_total + unmet;
        const removed = @min(pools.sheath_carbon_g, request);
        removed_sheath_total = removed;
        result.structural_sheath_carbon_g = removed * (1.0 - mobile_fraction);
        result.mobile_carbon_g += removed * mobile_fraction;
        result.symbiont_mobile_carbon_g += removed * symbiont_fraction;
        unmet = @max(0.0, request - removed);
        request = requested_nonfoliar * pools.husk_carbon_g / nonfoliar_total + unmet;
        result.husk_carbon_g = @min(pools.husk_carbon_g, request);
        unmet = @max(0.0, request - result.husk_carbon_g);
        request = requested_nonfoliar * pools.ear_carbon_g / nonfoliar_total + unmet;
        result.ear_carbon_g = @min(pools.ear_carbon_g, request);
        unmet = @max(0.0, request - result.ear_carbon_g);
        request = requested_nonfoliar * pools.grain_carbon_g / nonfoliar_total + unmet;
        result.grain_carbon_g = @min(pools.grain_carbon_g, request);
        unmet = @max(0.0, request - result.grain_carbon_g);
    } else {
        unmet += requested_nonfoliar;
    }

    const woody_total = pools.stalk_carbon_g + pools.reserve_carbon_g;
    const requested_woody = total_demand_g_c * harvested_woody_fraction;
    if (woody_total > requested_woody + unmet) {
        var request = requested_woody * pools.stalk_carbon_g / woody_total + unmet;
        result.stalk_carbon_g = @min(pools.stalk_carbon_g, request);
        unmet = @max(0.0, request - result.stalk_carbon_g);
        request = requested_woody * pools.reserve_carbon_g / woody_total + unmet;
        result.reserve_carbon_g = @min(pools.reserve_carbon_g, request);
        unmet = @max(0.0, request - result.reserve_carbon_g);
    } else {
        result.stalk_carbon_g = 0;
        result.reserve_carbon_g = 0;
        unmet = 0;
    }

    if (unmet > 0) {
        var removed = @min(@max(0.0, pools.leaf_carbon_g - removed_leaf), unmet);
        result.structural_leaf_carbon_g += removed * (1.0 - mobile_fraction);
        result.mobile_carbon_g += removed * mobile_fraction;
        result.symbiont_mobile_carbon_g += removed * symbiont_fraction;
        unmet = @max(0.0, unmet - removed);
        if (nonfoliar_total > 0) {
            var request = unmet * pools.sheath_carbon_g / nonfoliar_total;
            removed = @min(@max(0.0, pools.sheath_carbon_g - removed_sheath_total), request);
            removed_sheath_total += removed;
            result.structural_sheath_carbon_g += removed * (1.0 - mobile_fraction);
            result.mobile_carbon_g += removed * mobile_fraction;
            result.symbiont_mobile_carbon_g += removed * symbiont_fraction;
            unmet = @max(0.0, unmet - removed);
            request = unmet * pools.husk_carbon_g / nonfoliar_total;
            removed = @min(@max(0.0, pools.husk_carbon_g - result.husk_carbon_g), request);
            result.husk_carbon_g += removed;
            unmet = @max(0.0, unmet - removed);
            request = unmet * pools.ear_carbon_g / nonfoliar_total;
            removed = @min(@max(0.0, pools.ear_carbon_g - result.ear_carbon_g), request);
            result.ear_carbon_g += removed;
            unmet = @max(0.0, request - removed);
            request = unmet * pools.grain_carbon_g / nonfoliar_total;
            removed = @min(@max(0.0, pools.grain_carbon_g - result.grain_carbon_g), request);
            result.grain_carbon_g += removed;
            unmet = @max(0.0, request - removed);
        }
    }
    result.unmet_carbon_g = unmet;
    return result;
}

pub const ReproductiveRetention = struct { husk_remaining: f64, husk_unexported: f64, ear_remaining: f64, ear_unexported: f64, grain_remaining: f64, grain_unexported: f64 };

/// Operands for the exact GROSUB 9532-9573 reproductive-organ selector.
pub const SourceOrderReproductiveRetentionInput = struct {
    grazing: bool,
    reproductive_organs_reached_by_cut: bool,
    grain_or_pruning: bool,
    thinning_fraction: f64,
    harvested_nonfoliar_fraction: f64,
    total_husk_carbon_g_c: f64,
    total_ear_carbon_g_c: f64,
    total_grain_carbon_g_c: f64,
    grazed_husk_carbon_g_c: f64,
    grazed_ear_carbon_g_c: f64,
    grazed_grain_carbon_g_c: f64,
    plant_presence_threshold_g_c: f64,
};

/// Preserves GROSUB's strict plant-level `ZEROP` gates and branch order.
pub fn sourceOrderReproductiveRetention(input: SourceOrderReproductiveRetentionInput) !ReproductiveRetention {
    inline for (.{
        input.thinning_fraction,
        input.harvested_nonfoliar_fraction,
        input.total_husk_carbon_g_c,
        input.total_ear_carbon_g_c,
        input.total_grain_carbon_g_c,
        input.grazed_husk_carbon_g_c,
        input.grazed_ear_carbon_g_c,
        input.grazed_grain_carbon_g_c,
        input.plant_presence_threshold_g_c,
    }) |value| if (!std.math.isFinite(value)) return error.NonFiniteReproductiveHarvestInput;
    if (input.thinning_fraction < 0 or input.thinning_fraction > 1 or
        input.harvested_nonfoliar_fraction < 0 or input.harvested_nonfoliar_fraction > 1 or
        input.total_husk_carbon_g_c < 0 or input.total_ear_carbon_g_c < 0 or
        input.total_grain_carbon_g_c < 0 or input.grazed_husk_carbon_g_c < 0 or
        input.grazed_ear_carbon_g_c < 0 or input.grazed_grain_carbon_g_c < 0 or
        input.plant_presence_threshold_g_c < 0)
        return error.InvalidReproductiveHarvestInput;

    if (!input.grazing) {
        const reached = input.reproductive_organs_reached_by_cut or input.grain_or_pruning;
        const remaining = if (reached)
            if (input.thinning_fraction == 0)
                1 - input.harvested_nonfoliar_fraction
            else
                1 - input.thinning_fraction
        else
            1 - input.thinning_fraction;
        const unexported = if (reached and input.thinning_fraction != 0)
            1 - input.harvested_nonfoliar_fraction * input.thinning_fraction
        else
            remaining;
        return .{
            .husk_remaining = remaining,
            .husk_unexported = unexported,
            .ear_remaining = remaining,
            .ear_unexported = unexported,
            .grain_remaining = remaining,
            .grain_unexported = unexported,
        };
    }

    const husk_remaining = if (input.total_husk_carbon_g_c > input.plant_presence_threshold_g_c)
        std.math.clamp(1 - input.grazed_husk_carbon_g_c / input.total_husk_carbon_g_c, 0, 1)
    else
        1;
    const ear_remaining = if (input.total_ear_carbon_g_c > input.plant_presence_threshold_g_c)
        std.math.clamp(1 - input.grazed_ear_carbon_g_c / input.total_ear_carbon_g_c, 0, 1)
    else
        1;
    const grain_remaining = if (input.total_grain_carbon_g_c > input.plant_presence_threshold_g_c)
        std.math.clamp(1 - input.grazed_grain_carbon_g_c / input.total_grain_carbon_g_c, 0, 1)
    else
        1;
    return .{
        .husk_remaining = husk_remaining,
        .husk_unexported = husk_remaining,
        .ear_remaining = ear_remaining,
        .ear_unexported = ear_remaining,
        .grain_remaining = grain_remaining,
        .grain_unexported = grain_remaining,
    };
}

pub fn reproductiveRetention(grazing: bool, reproductive_organs_reached_by_cut: bool, grain_or_pruning: bool, thinning_fraction: f64, harvested_nonfoliar_fraction: f64, total_husk_carbon_g: f64, total_ear_carbon_g: f64, total_grain_carbon_g: f64, grazed_husk_carbon_g: f64, grazed_ear_carbon_g: f64, grazed_grain_carbon_g: f64) !ReproductiveRetention {
    inline for (.{ thinning_fraction, harvested_nonfoliar_fraction, total_husk_carbon_g, total_ear_carbon_g, total_grain_carbon_g, grazed_husk_carbon_g, grazed_ear_carbon_g, grazed_grain_carbon_g }) |value| if (!std.math.isFinite(value)) return error.NonFiniteReproductiveHarvestInput;
    if (thinning_fraction < 0 or thinning_fraction > 1 or harvested_nonfoliar_fraction < 0 or harvested_nonfoliar_fraction > 1 or total_husk_carbon_g < 0 or total_ear_carbon_g < 0 or total_grain_carbon_g < 0 or grazed_husk_carbon_g < 0 or grazed_ear_carbon_g < 0 or grazed_grain_carbon_g < 0) return error.InvalidReproductiveHarvestInput;
    if (grazing) return .{
        .husk_remaining = if (total_husk_carbon_g > 0) std.math.clamp(1.0 - grazed_husk_carbon_g / total_husk_carbon_g, 0, 1) else 1,
        .husk_unexported = if (total_husk_carbon_g > 0) std.math.clamp(1.0 - grazed_husk_carbon_g / total_husk_carbon_g, 0, 1) else 1,
        .ear_remaining = if (total_ear_carbon_g > 0) std.math.clamp(1.0 - grazed_ear_carbon_g / total_ear_carbon_g, 0, 1) else 1,
        .ear_unexported = if (total_ear_carbon_g > 0) std.math.clamp(1.0 - grazed_ear_carbon_g / total_ear_carbon_g, 0, 1) else 1,
        .grain_remaining = if (total_grain_carbon_g > 0) std.math.clamp(1.0 - grazed_grain_carbon_g / total_grain_carbon_g, 0, 1) else 1,
        .grain_unexported = if (total_grain_carbon_g > 0) std.math.clamp(1.0 - grazed_grain_carbon_g / total_grain_carbon_g, 0, 1) else 1,
    };
    const reached = reproductive_organs_reached_by_cut or grain_or_pruning;
    const remaining = if (reached) (if (thinning_fraction == 0) 1.0 - harvested_nonfoliar_fraction else 1.0 - thinning_fraction) else 1.0 - thinning_fraction;
    const unexported = if (reached and thinning_fraction != 0) 1.0 - harvested_nonfoliar_fraction * thinning_fraction else remaining;
    return .{ .husk_remaining = remaining, .husk_unexported = unexported, .ear_remaining = remaining, .ear_unexported = unexported, .grain_remaining = remaining, .grain_unexported = unexported };
}

pub const ReproductiveHarvestResult = struct { products: HarvestProducts, harvested_grain: ElementalMass };

pub fn harvestReproductiveOrgans(state: *State, branch: usize, retention: ReproductiveRetention) !ReproductiveHarvestResult {
    if (branch >= state.branch_husk_carbon_g.len) return error.CanopyBranchIndexOutOfBounds;
    const husk = try partitionHarvest(.{ .carbon_g = state.branch_husk_carbon_g[branch], .nitrogen_g = state.branch_husk_nitrogen_g[branch], .phosphorus_g = state.branch_husk_phosphorus_g[branch] }, retention.husk_remaining, retention.husk_unexported);
    const ear = try partitionHarvest(.{ .carbon_g = state.branch_ear_carbon_g[branch], .nitrogen_g = state.branch_ear_nitrogen_g[branch], .phosphorus_g = state.branch_ear_phosphorus_g[branch] }, retention.ear_remaining, retention.ear_unexported);
    const grain_mass: ElementalMass = .{ .carbon_g = state.branch_grain_carbon_g[branch], .nitrogen_g = state.branch_grain_nitrogen_g[branch], .phosphorus_g = state.branch_grain_phosphorus_g[branch] };
    const grain = try partitionHarvest(grain_mass, retention.grain_remaining, retention.grain_unexported);
    state.branch_husk_carbon_g[branch] = husk.remaining.carbon_g;
    state.branch_husk_nitrogen_g[branch] = husk.remaining.nitrogen_g;
    state.branch_husk_phosphorus_g[branch] = husk.remaining.phosphorus_g;
    state.branch_ear_carbon_g[branch] = ear.remaining.carbon_g;
    state.branch_ear_nitrogen_g[branch] = ear.remaining.nitrogen_g;
    state.branch_ear_phosphorus_g[branch] = ear.remaining.phosphorus_g;
    state.branch_grain_carbon_g[branch] = grain.remaining.carbon_g;
    state.branch_grain_nitrogen_g[branch] = grain.remaining.nitrogen_g;
    state.branch_grain_phosphorus_g[branch] = grain.remaining.phosphorus_g;
    state.branch_potential_seed_site_count[branch] *= retention.grain_remaining;
    state.branch_seed_count[branch] *= retention.grain_remaining;
    return .{
        .products = .{
            .ecosystem_export = addElementalMass(addElementalMass(husk.products.ecosystem_export, ear.products.ecosystem_export), grain.products.ecosystem_export),
            .litter = addElementalMass(addElementalMass(husk.products.litter, ear.products.litter), grain.products.litter),
        },
        .harvested_grain = .{ .carbon_g = (1.0 - retention.grain_remaining) * grain_mass.carbon_g, .nitrogen_g = (1.0 - retention.grain_remaining) * grain_mass.nitrogen_g, .phosphorus_g = (1.0 - retention.grain_remaining) * grain_mass.phosphorus_g },
    };
}

fn addElementalMass(left: ElementalMass, right: ElementalMass) ElementalMass {
    return .{ .carbon_g = left.carbon_g + right.carbon_g, .nitrogen_g = left.nitrogen_g + right.nitrogen_g, .phosphorus_g = left.phosphorus_g + right.phosphorus_g };
}

pub fn cuttingHeightForLeafAreaRemoval(removal_fraction: f64, layer_boundary_height_m: []const f64, leaf_area_by_layer_m2: []const f64) !f64 {
    if (!std.math.isFinite(removal_fraction) or removal_fraction < 0 or removal_fraction > 1 or layer_boundary_height_m.len != leaf_area_by_layer_m2.len + 1 or leaf_area_by_layer_m2.len == 0) return error.InvalidLeafAreaHarvestGeometry;
    var total_leaf_area: f64 = 0;
    for (leaf_area_by_layer_m2) |area| {
        if (!std.math.isFinite(area) or area < 0) return error.InvalidLeafAreaHarvestGeometry;
        total_leaf_area += area;
    }
    for (1..layer_boundary_height_m.len) |index| if (!std.math.isFinite(layer_boundary_height_m[index - 1]) or !std.math.isFinite(layer_boundary_height_m[index]) or layer_boundary_height_m[index] < layer_boundary_height_m[index - 1]) return error.InvalidLeafAreaHarvestGeometry;
    const target_remaining_leaf_area = (1.0 - removal_fraction) * total_leaf_area;
    var accumulated_leaf_area: f64 = 0;
    var cutting_height_m: f64 = 0;
    for (leaf_area_by_layer_m2, 0..) |layer_area, layer| {
        const lower = layer_boundary_height_m[layer];
        const upper = layer_boundary_height_m[layer + 1];
        if (upper > lower and layer_area > 0 and accumulated_leaf_area < target_remaining_leaf_area) {
            if (accumulated_leaf_area + layer_area > target_remaining_leaf_area)
                cutting_height_m = lower + (target_remaining_leaf_area - accumulated_leaf_area) / layer_area * (upper - lower)
            else
                cutting_height_m = 0;
            accumulated_leaf_area += layer_area;
        }
    }
    return cutting_height_m;
}

pub const LayerHarvestRetention = struct { remaining_fraction: f64, unexported_fraction: f64, height_below_cut_fraction: f64 };

pub const RemainingNodeLeaf = struct {
    area_m2: f64,
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
    protein_mass_g: f64,
};

/// Exact GROSUB 8962-9051 reconstruction of one node from its runtime layers.
pub fn sourceOrderRemainingNodeLeaf(
    layer_area_m2: []const f64,
    layer_carbon_g_c: []const f64,
    layer_nitrogen_g_n: []const f64,
    layer_phosphorus_g_p: []const f64,
    previous_node_area_m2: f64,
    previous_protein_mass_g: f64,
    plant_presence_threshold: f64,
) !RemainingNodeLeaf {
    const layer_count = layer_area_m2.len;
    if (layer_count == 0 or layer_carbon_g_c.len != layer_count or
        layer_nitrogen_g_n.len != layer_count or layer_phosphorus_g_p.len != layer_count)
        return error.NodeLeafDimensionMismatch;
    inline for (.{ previous_node_area_m2, previous_protein_mass_g, plant_presence_threshold }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidNodeLeafState;
    var result: RemainingNodeLeaf = .{
        .area_m2 = 0,
        .carbon_g_c = 0,
        .nitrogen_g_n = 0,
        .phosphorus_g_p = 0,
        .protein_mass_g = 0,
    };
    for (layer_area_m2, layer_carbon_g_c, layer_nitrogen_g_n, layer_phosphorus_g_p) |area, carbon, nitrogen, phosphorus| {
        inline for (.{ area, carbon, nitrogen, phosphorus }) |value|
            if (!std.math.isFinite(value) or value < 0) return error.InvalidNodeLeafState;
        result.area_m2 += area;
        result.carbon_g_c += carbon;
        result.nitrogen_g_n += nitrogen;
        result.phosphorus_g_p += phosphorus;
    }
    result.protein_mass_g = if (previous_node_area_m2 > plant_presence_threshold)
        previous_protein_mass_g * result.area_m2 / previous_node_area_m2
    else
        0;
    return result;
}

/// Exact GROSUB 8999-9022 retained-plant and retained-ecosystem fractions for
/// sheath/petiole and internode processing after leaf reconstruction.
pub fn sourceOrderNodeOrganRetention(
    grazing: bool,
    no_harvest_kind: bool,
    previous_leaf_carbon_g_c: f64,
    remaining_leaf_carbon_g_c: f64,
    leaf_harvest_fraction: f64,
    nonfoliar_harvest_fraction: f64,
    thinning_fraction: f64,
    plant_presence_threshold_g_c: f64,
) !LayerHarvestRetention {
    inline for (.{
        previous_leaf_carbon_g_c,
        remaining_leaf_carbon_g_c,
        leaf_harvest_fraction,
        nonfoliar_harvest_fraction,
        thinning_fraction,
        plant_presence_threshold_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidLayerHarvestInput;
    if (remaining_leaf_carbon_g_c > previous_leaf_carbon_g_c or leaf_harvest_fraction > 1 or
        nonfoliar_harvest_fraction > 1 or thinning_fraction > 1)
        return error.InvalidLayerHarvestInput;
    if (grazing) return .{ .remaining_fraction = 0, .unexported_fraction = 0, .height_below_cut_fraction = 0 };
    if (previous_leaf_carbon_g_c > plant_presence_threshold_g_c and leaf_harvest_fraction > 0) {
        const retention = @max(0, @min(1, 1 -
            (1 - @max(0, remaining_leaf_carbon_g_c) / previous_leaf_carbon_g_c) *
                nonfoliar_harvest_fraction / leaf_harvest_fraction));
        return .{ .remaining_fraction = retention, .unexported_fraction = retention, .height_below_cut_fraction = 0 };
    }
    if (thinning_fraction == 0) {
        const retention = 1 - nonfoliar_harvest_fraction;
        return .{ .remaining_fraction = retention, .unexported_fraction = retention, .height_below_cut_fraction = 0 };
    }
    const remaining = 1 - thinning_fraction;
    const unexported = if (no_harvest_kind)
        1 - nonfoliar_harvest_fraction * thinning_fraction
    else
        remaining;
    return .{ .remaining_fraction = remaining, .unexported_fraction = unexported, .height_below_cut_fraction = 0 };
}

pub fn layerHarvestRetention(layer_lower_height_m: f64, layer_upper_height_m: f64, cutting_height_m: f64, pruning: bool, no_harvest_kind: bool, thinning_fraction: f64, harvested_fraction: f64) !LayerHarvestRetention {
    inline for (.{ layer_lower_height_m, layer_upper_height_m, cutting_height_m, thinning_fraction, harvested_fraction }) |value| if (!std.math.isFinite(value)) return error.NonFiniteLayerHarvestInput;
    if (layer_upper_height_m < layer_lower_height_m or thinning_fraction < 0 or thinning_fraction > 1 or harvested_fraction < 0 or harvested_fraction > 1) return error.InvalidLayerHarvestInput;
    const height_fraction = if (pruning) 0 else if (layer_upper_height_m > layer_lower_height_m) std.math.clamp(1.0 - (layer_upper_height_m - cutting_height_m) / (layer_upper_height_m - layer_lower_height_m), 0, 1) else 1;
    if (thinning_fraction == 0) {
        const retention = @max(0.0, 1.0 - (1.0 - height_fraction) * harvested_fraction);
        return .{ .remaining_fraction = retention, .unexported_fraction = retention, .height_below_cut_fraction = height_fraction };
    }
    const remaining = @max(0.0, 1.0 - thinning_fraction);
    const unexported = if (no_harvest_kind) 1.0 - (1.0 - height_fraction) * harvested_fraction * thinning_fraction else remaining;
    return .{ .remaining_fraction = remaining, .unexported_fraction = unexported, .height_below_cut_fraction = height_fraction };
}

pub const LayerLeafHarvestProducts = struct {
    foliar: HarvestProducts = .{},
    woody: HarvestProducts = .{},
    removed_leaf_area_m2: f64 = 0,
};

pub const NodeOrganHarvestProducts = struct { nonwoody: HarvestProducts = .{}, woody: HarvestProducts = .{} };

/// Commits GROSUB WGLFL/WGLFLN/WGLFLP/ARLFL removal for one runtime
/// canopy layer sample and reconciles its node and branch aggregates.
pub fn harvestLeafLayerSample(state: *State, branch: usize, node_within_branch: usize, sample_within_node: usize, retention: LayerHarvestRetention, carbon_woody_fraction: [2]f64, nitrogen_woody_fraction: [2]f64, phosphorus_woody_fraction: [2]f64, scale_stalk_area: bool) !LayerLeafHarvestProducts {
    inline for (.{ carbon_woody_fraction, nitrogen_woody_fraction, phosphorus_woody_fraction }) |fractions| {
        for (fractions) |fraction| if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidWoodyFraction;
        if (@abs(fractions[0] + fractions[1] - 1.0) > 1e-8) return error.InvalidWoodyFraction;
    }
    if (!std.math.isFinite(retention.remaining_fraction) or !std.math.isFinite(retention.unexported_fraction) or retention.remaining_fraction < 0 or retention.unexported_fraction < retention.remaining_fraction or retention.unexported_fraction > 1) return error.InvalidHarvestRetentionFraction;
    const nodes = try state.nodeRange(branch);
    if (node_within_branch >= nodes.end - nodes.first) return error.CanopyNodeIndexOutOfBounds;
    const node = nodes.first + node_within_branch;
    const samples = try state.sampleRange(node);
    if (sample_within_node >= samples.end - samples.first) return error.CanopySampleIndexOutOfBounds;
    const sample = samples.first + sample_within_node;
    const initial_mass: ElementalMass = .{ .carbon_g = state.sample_leaf_carbon_g[sample], .nitrogen_g = state.sample_leaf_nitrogen_g[sample], .phosphorus_g = state.sample_leaf_phosphorus_g[sample] };
    inline for (@typeInfo(ElementalMass).@"struct".fields) |field| if (!std.math.isFinite(@field(initial_mass, field.name)) or @field(initial_mass, field.name) < 0) return error.InvalidHarvestMass;
    const initial_area_m2 = state.sample_leaf_area_m2[sample];
    if (!std.math.isFinite(initial_area_m2) or initial_area_m2 < 0) return error.InvalidLeafAreaHarvestGeometry;
    const removed_fraction = 1.0 - retention.remaining_fraction;
    const export_fraction = 1.0 - retention.unexported_fraction;
    const litter_fraction = retention.unexported_fraction - retention.remaining_fraction;
    var result: LayerLeafHarvestProducts = .{ .removed_leaf_area_m2 = removed_fraction * initial_area_m2 };
    inline for (@typeInfo(ElementalMass).@"struct".fields, .{ carbon_woody_fraction, nitrogen_woody_fraction, phosphorus_woody_fraction }) |field, fractions| {
        const initial = @field(initial_mass, field.name);
        @field(result.woody.ecosystem_export, field.name) = export_fraction * initial * fractions[0];
        @field(result.woody.litter, field.name) = litter_fraction * initial * fractions[0];
        @field(result.foliar.ecosystem_export, field.name) = export_fraction * initial * fractions[1];
        @field(result.foliar.litter, field.name) = litter_fraction * initial * fractions[1];
    }
    const old_node_area_m2 = state.node_leaf_area_m2[node];
    state.sample_leaf_area_m2[sample] = retention.remaining_fraction * initial_area_m2;
    state.sample_leaf_carbon_g[sample] = retention.remaining_fraction * initial_mass.carbon_g;
    state.sample_leaf_nitrogen_g[sample] = retention.remaining_fraction * initial_mass.nitrogen_g;
    state.sample_leaf_phosphorus_g[sample] = retention.remaining_fraction * initial_mass.phosphorus_g;
    state.sample_exposed_leaf_area_m2[sample] *= retention.remaining_fraction;
    if (scale_stalk_area) state.sample_stalk_area_m2[sample] *= retention.remaining_fraction;
    state.node_leaf_area_m2[node] = @max(0.0, old_node_area_m2 - result.removed_leaf_area_m2);
    state.node_leaf_carbon_g[node] = @max(0.0, state.node_leaf_carbon_g[node] - removed_fraction * initial_mass.carbon_g);
    state.node_leaf_nitrogen_g[node] = @max(0.0, state.node_leaf_nitrogen_g[node] - removed_fraction * initial_mass.nitrogen_g);
    state.node_leaf_phosphorus_g[node] = @max(0.0, state.node_leaf_phosphorus_g[node] - removed_fraction * initial_mass.phosphorus_g);
    state.node_leaf_protein_g[node] *= if (old_node_area_m2 > 0) state.node_leaf_area_m2[node] / old_node_area_m2 else 0;
    state.branch_leaf_area_m2[branch] = @max(0.0, state.branch_leaf_area_m2[branch] - result.removed_leaf_area_m2);
    state.branch_leaf_carbon_g[branch] = @max(0.0, state.branch_leaf_carbon_g[branch] - removed_fraction * initial_mass.carbon_g);
    state.branch_leaf_nitrogen_g[branch] = @max(0.0, state.branch_leaf_nitrogen_g[branch] - removed_fraction * initial_mass.nitrogen_g);
    state.branch_leaf_phosphorus_g[branch] = @max(0.0, state.branch_leaf_phosphorus_g[branch] - removed_fraction * initial_mass.phosphorus_g);
    return result;
}

/// Commits GROSUB sheath/petiole removal for one node. Direct cutting uses
/// geometric truncation at the cutting plane; pruning, thinning, and grazing
/// scale length with retained biomass as in the source branches.
pub fn harvestNodeSheath(state: *State, branch: usize, node_within_branch: usize, remaining_fraction: f64, unexported_fraction: f64, carbon_woody_fraction: [2]f64, nitrogen_woody_fraction: [2]f64, phosphorus_woody_fraction: [2]f64, use_cutting_plane: bool, cutting_height_m: f64) !NodeOrganHarvestProducts {
    if (!std.math.isFinite(remaining_fraction) or !std.math.isFinite(unexported_fraction) or !std.math.isFinite(cutting_height_m) or remaining_fraction < 0 or unexported_fraction < remaining_fraction or unexported_fraction > 1 or cutting_height_m < 0) return error.InvalidNodeSheathHarvestInput;
    inline for (.{ carbon_woody_fraction, nitrogen_woody_fraction, phosphorus_woody_fraction }) |fractions| {
        for (fractions) |fraction| if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidWoodyFraction;
        if (@abs(fractions[0] + fractions[1] - 1.0) > 1e-8) return error.InvalidWoodyFraction;
    }
    const nodes = try state.nodeRange(branch);
    if (node_within_branch >= nodes.end - nodes.first) return error.CanopyNodeIndexOutOfBounds;
    const node = nodes.first + node_within_branch;
    const initial_mass: ElementalMass = .{ .carbon_g = state.node_sheath_carbon_g[node], .nitrogen_g = state.node_sheath_nitrogen_g[node], .phosphorus_g = state.node_sheath_phosphorus_g[node] };
    inline for (@typeInfo(ElementalMass).@"struct".fields) |field| if (!std.math.isFinite(@field(initial_mass, field.name)) or @field(initial_mass, field.name) < 0) return error.InvalidHarvestMass;
    const export_fraction = 1.0 - unexported_fraction;
    const litter_fraction = unexported_fraction - remaining_fraction;
    var result: NodeOrganHarvestProducts = .{};
    inline for (@typeInfo(ElementalMass).@"struct".fields, .{ carbon_woody_fraction, nitrogen_woody_fraction, phosphorus_woody_fraction }) |field, fractions| {
        const initial = @field(initial_mass, field.name);
        @field(result.woody.ecosystem_export, field.name) = export_fraction * initial * fractions[0];
        @field(result.woody.litter, field.name) = litter_fraction * initial * fractions[0];
        @field(result.nonwoody.ecosystem_export, field.name) = export_fraction * initial * fractions[1];
        @field(result.nonwoody.litter, field.name) = litter_fraction * initial * fractions[1];
    }
    const removed_fraction = 1.0 - remaining_fraction;
    state.node_sheath_carbon_g[node] *= remaining_fraction;
    state.node_sheath_nitrogen_g[node] *= remaining_fraction;
    state.node_sheath_phosphorus_g[node] *= remaining_fraction;
    state.node_sheath_protein_g[node] *= remaining_fraction;
    state.branch_sheath_carbon_g[branch] = @max(0.0, state.branch_sheath_carbon_g[branch] - removed_fraction * initial_mass.carbon_g);
    state.branch_sheath_nitrogen_g[branch] = @max(0.0, state.branch_sheath_nitrogen_g[branch] - removed_fraction * initial_mass.nitrogen_g);
    state.branch_sheath_phosphorus_g[branch] = @max(0.0, state.branch_sheath_phosphorus_g[branch] - removed_fraction * initial_mass.phosphorus_g);
    const initial_height_m = state.node_sheath_height_m[node];
    if (use_cutting_plane and initial_height_m > 0) {
        const fraction_above_cut = std.math.clamp((state.node_height_m[node] + initial_height_m - cutting_height_m) / initial_height_m, 0, 1);
        state.node_sheath_height_m[node] = (1.0 - fraction_above_cut) * initial_height_m;
    } else state.node_sheath_height_m[node] *= remaining_fraction;
    return result;
}

pub fn internodeHarvestRetention(node_height_m: f64, internode_length_m: f64, cutting_height_m: f64, pruning: bool, thinning_fraction: f64, woody_harvest_fraction: f64, grazing: bool, grazed_stalk_carbon_g: f64, total_stalk_carbon_g: f64) !f64 {
    inline for (.{ node_height_m, internode_length_m, cutting_height_m, thinning_fraction, woody_harvest_fraction, grazed_stalk_carbon_g, total_stalk_carbon_g }) |value| if (!std.math.isFinite(value)) return error.NonFiniteInternodeHarvestInput;
    if (node_height_m < 0 or internode_length_m < 0 or cutting_height_m < 0 or thinning_fraction < 0 or thinning_fraction > 1 or woody_harvest_fraction < 0 or woody_harvest_fraction > 1 or grazed_stalk_carbon_g < 0 or total_stalk_carbon_g < 0) return error.InvalidInternodeHarvestInput;
    if (grazing) return if (total_stalk_carbon_g > 0) std.math.clamp(1.0 - grazed_stalk_carbon_g / total_stalk_carbon_g, 0, 1) else 1;
    if (internode_length_m <= 0) return 1;
    const fraction_above_cut = if (pruning) 0 else std.math.clamp((node_height_m - cutting_height_m) / internode_length_m, 0, 1);
    return if (thinning_fraction == 0) @max(0.0, 1.0 - fraction_above_cut * woody_harvest_fraction) else @max(0.0, 1.0 - thinning_fraction);
}

/// Exact GROSUB 9339-9369 branch stalk retained-plant/ecosystem fractions.
pub fn sourceOrderBranchStalkRetention(
    grazing: bool,
    no_harvest_kind: bool,
    pruning: bool,
    maximum_internode_height_m: f64,
    cutting_height_m: f64,
    thinning_fraction: f64,
    woody_harvest_fraction: f64,
    total_stalk_carbon_g_c: f64,
    allocated_woody_demand_g_c: f64,
    removed_stalk_carbon_g_c: f64,
    plant_presence_threshold_g_c: f64,
) !LayerHarvestRetention {
    inline for (.{
        maximum_internode_height_m,
        cutting_height_m,
        thinning_fraction,
        woody_harvest_fraction,
        total_stalk_carbon_g_c,
        allocated_woody_demand_g_c,
        removed_stalk_carbon_g_c,
        plant_presence_threshold_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidInternodeHarvestInput;
    if (thinning_fraction > 1 or woody_harvest_fraction > 1) return error.InvalidInternodeHarvestInput;
    if (grazing) {
        const retention = if (total_stalk_carbon_g_c > plant_presence_threshold_g_c)
            @max(0, @min(1, 1 - (removed_stalk_carbon_g_c + allocated_woody_demand_g_c) / total_stalk_carbon_g_c))
        else
            1;
        return .{ .remaining_fraction = retention, .unexported_fraction = retention, .height_below_cut_fraction = 0 };
    }
    if (maximum_internode_height_m == 0)
        return .{ .remaining_fraction = 1, .unexported_fraction = 1, .height_below_cut_fraction = 1 };
    const height_fraction = if (pruning) 0 else @max(0, @min(1, cutting_height_m / maximum_internode_height_m));
    if (thinning_fraction == 0) {
        const retention = @max(0, 1 - (1 - height_fraction) * woody_harvest_fraction);
        return .{ .remaining_fraction = retention, .unexported_fraction = retention, .height_below_cut_fraction = height_fraction };
    }
    const remaining = @max(0, 1 - thinning_fraction);
    const unexported = if (no_harvest_kind)
        1 - (1 - height_fraction) * woody_harvest_fraction * thinning_fraction
    else
        remaining;
    return .{ .remaining_fraction = remaining, .unexported_fraction = unexported, .height_below_cut_fraction = height_fraction };
}

/// Exact GROSUB 9474-9490 reserve retention after branch stalk publication.
pub fn sourceOrderStalkReserveRetention(
    grazing: bool,
    remaining_stalk_carbon_g_c: f64,
    stalk_retention: LayerHarvestRetention,
    reserve_carbon_g_c: f64,
    removed_reserve_carbon_g_c: f64,
    plant_presence_threshold_g_c: f64,
) !LayerHarvestRetention {
    inline for (.{ remaining_stalk_carbon_g_c, reserve_carbon_g_c, removed_reserve_carbon_g_c, plant_presence_threshold_g_c }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidHarvestMass;
    if (!grazing) {
        if (remaining_stalk_carbon_g_c > plant_presence_threshold_g_c) return stalk_retention;
        return .{ .remaining_fraction = 0, .unexported_fraction = 0, .height_below_cut_fraction = 0 };
    }
    const retention = if (reserve_carbon_g_c > plant_presence_threshold_g_c)
        @max(0, @min(1, 1 - removed_reserve_carbon_g_c / reserve_carbon_g_c))
    else
        0;
    return .{ .remaining_fraction = retention, .unexported_fraction = retention, .height_below_cut_fraction = 0 };
}

/// Applies the source node-stalk retention after branch-level harvested stalk
/// products have been accounted for.
pub fn commitInternodeHarvest(state: *State, branch: usize, node_within_branch: usize, remaining_fraction: f64, direct_cut_without_thinning: bool, cutting_height_m: f64) !void {
    if (!std.math.isFinite(remaining_fraction) or !std.math.isFinite(cutting_height_m) or remaining_fraction < 0 or remaining_fraction > 1 or cutting_height_m < 0) return error.InvalidInternodeHarvestInput;
    const nodes = try state.nodeRange(branch);
    if (node_within_branch >= nodes.end - nodes.first) return error.CanopyNodeIndexOutOfBounds;
    const node = nodes.first + node_within_branch;
    state.node_internode_carbon_g[node] *= remaining_fraction;
    state.node_internode_nitrogen_g[node] *= remaining_fraction;
    state.node_internode_phosphorus_g[node] *= remaining_fraction;
    if (direct_cut_without_thinning) {
        state.node_internode_length_m[node] *= remaining_fraction;
        state.node_height_m[node] = @min(state.node_height_m[node], cutting_height_m);
    }
}

pub fn grazingCarbonDemandGPerH(animal_grazing: bool, grazer_biomass_g_living_mass_per_m2: f64, specific_consumption_g_dry_matter_per_g_living_mass_d: f64, horizontal_cell_area_m2: f64, leaf_plus_stalk_area_m2: f64, canopy_growth_temperature_factor: f64, cell_shoot_carbon_g: f64, landscape_average_shoot_carbon_g: f64) !f64 {
    inline for (.{ grazer_biomass_g_living_mass_per_m2, specific_consumption_g_dry_matter_per_g_living_mass_d, horizontal_cell_area_m2, leaf_plus_stalk_area_m2, canopy_growth_temperature_factor, cell_shoot_carbon_g, landscape_average_shoot_carbon_g }) |value| if (!std.math.isFinite(value)) return error.NonFiniteGrazingDemandInput;
    if (grazer_biomass_g_living_mass_per_m2 < 0 or specific_consumption_g_dry_matter_per_g_living_mass_d < 0 or horizontal_cell_area_m2 <= 0 or leaf_plus_stalk_area_m2 < 0 or canopy_growth_temperature_factor < 0 or cell_shoot_carbon_g < 0 or landscape_average_shoot_carbon_g < 0) return error.InvalidGrazingDemandInput;
    if (landscape_average_shoot_carbon_g == 0) return 0;
    const spatial_share = cell_shoot_carbon_g / landscape_average_shoot_carbon_g;
    return if (animal_grazing)
        grazer_biomass_g_living_mass_per_m2 * specific_consumption_g_dry_matter_per_g_living_mass_d * horizontal_cell_area_m2 * 0.5 / 24.0 * spatial_share
    else
        grazer_biomass_g_living_mass_per_m2 * specific_consumption_g_dry_matter_per_g_living_mass_d * leaf_plus_stalk_area_m2 * 0.5 / 24.0 * canopy_growth_temperature_factor * spatial_share;
}

/// Exact GROSUB 8650-8662 demand gate using plant ZEROP.
pub fn sourceOrderGrazingCarbonDemandGPerH(
    animal_grazing: bool,
    grazer_biomass_g_living_mass_per_m2: f64,
    specific_consumption_g_dry_matter_per_g_living_mass_d: f64,
    horizontal_cell_area_m2: f64,
    leaf_plus_stalk_area_m2: f64,
    canopy_growth_temperature_factor: f64,
    cell_shoot_carbon_g: f64,
    landscape_average_shoot_carbon_g: f64,
    plant_presence_threshold_g_c: f64,
) !f64 {
    if (!std.math.isFinite(plant_presence_threshold_g_c) or plant_presence_threshold_g_c < 0)
        return error.InvalidGrazingDemandInput;
    if (landscape_average_shoot_carbon_g <= plant_presence_threshold_g_c) return 0;
    return grazingCarbonDemandGPerH(
        animal_grazing,
        grazer_biomass_g_living_mass_per_m2,
        specific_consumption_g_dry_matter_per_g_living_mass_d,
        horizontal_cell_area_m2,
        leaf_plus_stalk_area_m2,
        canopy_growth_temperature_factor,
        cell_shoot_carbon_g,
        landscape_average_shoot_carbon_g,
    );
}

fn partitionHarvest(mass: ElementalMass, remaining_fraction: f64, unexported_fraction: f64) !HarvestPartition {
    if (!std.math.isFinite(remaining_fraction) or !std.math.isFinite(unexported_fraction) or remaining_fraction < 0 or unexported_fraction < remaining_fraction or unexported_fraction > 1) return error.InvalidHarvestRetentionFraction;
    inline for (@typeInfo(ElementalMass).@"struct".fields) |field| if (!std.math.isFinite(@field(mass, field.name)) or @field(mass, field.name) < 0) return error.InvalidHarvestMass;
    var result: HarvestPartition = .{ .remaining = .{}, .products = .{} };
    inline for (@typeInfo(ElementalMass).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        @field(result.remaining, field.name) = remaining_fraction * value;
        @field(result.products.ecosystem_export, field.name) = (1.0 - unexported_fraction) * value;
        @field(result.products.litter, field.name) = (unexported_fraction - remaining_fraction) * value;
    }
    return result;
}

pub fn harvestBranchStalkAndReserve(state: *State, branch: usize, stalk_remaining_fraction: f64, stalk_unexported_fraction: f64, reserve_remaining_fraction: f64, reserve_unexported_fraction: f64) !HarvestProducts {
    if (branch >= state.branch_stalk_carbon_g.len) return error.CanopyBranchIndexOutOfBounds;
    const stalk = try partitionHarvest(.{ .carbon_g = state.branch_stalk_carbon_g[branch], .nitrogen_g = state.branch_stalk_nitrogen_g[branch], .phosphorus_g = state.branch_stalk_phosphorus_g[branch] }, stalk_remaining_fraction, stalk_unexported_fraction);
    const reserve = try partitionHarvest(.{ .carbon_g = state.branch_reserve_carbon_g[branch], .nitrogen_g = state.branch_reserve_nitrogen_g[branch], .phosphorus_g = state.branch_reserve_phosphorus_g[branch] }, reserve_remaining_fraction, reserve_unexported_fraction);
    state.branch_stalk_carbon_g[branch] = stalk.remaining.carbon_g;
    state.branch_stalk_nitrogen_g[branch] = stalk.remaining.nitrogen_g;
    state.branch_stalk_phosphorus_g[branch] = stalk.remaining.phosphorus_g;
    state.branch_sapwood_carbon_g[branch] *= stalk_remaining_fraction;
    state.branch_senescing_stalk_carbon_g[branch] *= stalk_remaining_fraction;
    state.branch_senescing_stalk_nitrogen_g[branch] *= stalk_remaining_fraction;
    state.branch_senescing_stalk_phosphorus_g[branch] *= stalk_remaining_fraction;
    state.branch_reserve_carbon_g[branch] = reserve.remaining.carbon_g;
    state.branch_reserve_nitrogen_g[branch] = reserve.remaining.nitrogen_g;
    state.branch_reserve_phosphorus_g[branch] = reserve.remaining.phosphorus_g;
    return .{
        .ecosystem_export = .{ .carbon_g = stalk.products.ecosystem_export.carbon_g + reserve.products.ecosystem_export.carbon_g, .nitrogen_g = stalk.products.ecosystem_export.nitrogen_g + reserve.products.ecosystem_export.nitrogen_g, .phosphorus_g = stalk.products.ecosystem_export.phosphorus_g + reserve.products.ecosystem_export.phosphorus_g },
        .litter = .{ .carbon_g = stalk.products.litter.carbon_g + reserve.products.litter.carbon_g, .nitrogen_g = stalk.products.litter.nitrogen_g + reserve.products.litter.nitrogen_g, .phosphorus_g = stalk.products.litter.phosphorus_g + reserve.products.litter.phosphorus_g },
    };
}

/// Applies GROSUB's mobile-pool retention ratio and the same ratio to every C4
/// intermediate belonging to the branch.
pub fn harvestBranchMobilePools(state: *State, branch: usize, remaining_fraction: f64) !ElementalMass {
    return harvestBranchMobilePoolsWithIntermediateRetention(
        state,
        branch,
        remaining_fraction,
        remaining_fraction,
    );
}

/// Applies independent host-mobile and biochemical-intermediate retention.
/// GROSUB uses the second fraction only for C4 plants with significant
/// pre-harvest host-mobile carbon.
pub fn harvestBranchMobilePoolsWithIntermediateRetention(
    state: *State,
    branch: usize,
    mobile_remaining_fraction: f64,
    intermediate_remaining_fraction: f64,
) !ElementalMass {
    inline for (.{ mobile_remaining_fraction, intermediate_remaining_fraction }) |fraction|
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.InvalidHarvestRetentionFraction;
    const nodes = try state.nodeRange(branch);
    inline for (.{ state.branch_mobile_carbon_g[branch], state.branch_mobile_nitrogen_g[branch], state.branch_mobile_phosphorus_g[branch] }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidHarvestMass;
    for (nodes.first..nodes.end) |node| {
        inline for (.{ "node_c3_nonstructural_carbon_g", "node_c4_mesophyll_nonstructural_carbon_g", "node_bundle_sheath_co2_carbon_g", "node_bundle_sheath_bicarbonate_carbon_g" }) |field_name| {
            const value = @field(state, field_name)[node];
            if (!std.math.isFinite(value) or value < 0) return error.InvalidHarvestMass;
        }
    }
    var removed: ElementalMass = .{
        .carbon_g = (1.0 - mobile_remaining_fraction) * state.branch_mobile_carbon_g[branch],
        .nitrogen_g = (1.0 - mobile_remaining_fraction) * state.branch_mobile_nitrogen_g[branch],
        .phosphorus_g = (1.0 - mobile_remaining_fraction) * state.branch_mobile_phosphorus_g[branch],
    };
    state.branch_mobile_carbon_g[branch] *= mobile_remaining_fraction;
    state.branch_mobile_nitrogen_g[branch] *= mobile_remaining_fraction;
    state.branch_mobile_phosphorus_g[branch] *= mobile_remaining_fraction;
    for (nodes.first..nodes.end) |node| {
        inline for (.{ "node_c3_nonstructural_carbon_g", "node_c4_mesophyll_nonstructural_carbon_g", "node_bundle_sheath_co2_carbon_g", "node_bundle_sheath_bicarbonate_carbon_g" }) |field_name| {
            removed.carbon_g += (1.0 - intermediate_remaining_fraction) * @field(state, field_name)[node];
            @field(state, field_name)[node] *= intermediate_remaining_fraction;
        }
    }
    return removed;
}

/// Conservative pairwise branch-reserve equilibration from the GROSUB main-
/// branch loop. Signed flux is positive from main_branch to other_branch.
pub fn equilibrateBranchReserves(state: *State, main_branch: usize, other_branch: usize, carbon_exchange_per_h: f64, nutrient_exchange_per_h: f64, timestep_h: f64, presence_threshold_g_c: f64) !ReserveExchange {
    inline for (.{ carbon_exchange_per_h, nutrient_exchange_per_h, timestep_h, presence_threshold_g_c }) |value| if (!std.math.isFinite(value)) return error.NonFiniteBranchReserveExchangeInput;
    if (main_branch >= state.branch_reserve_carbon_g.len or other_branch >= state.branch_reserve_carbon_g.len or main_branch == other_branch) return error.CanopyBranchIndexOutOfBounds;
    if (carbon_exchange_per_h < 0 or nutrient_exchange_per_h < 0 or timestep_h <= 0 or presence_threshold_g_c < 0) return error.InvalidBranchReserveExchangeInput;
    const other_sapwood = state.branch_sapwood_carbon_g[other_branch];
    if (other_sapwood <= presence_threshold_g_c) return .{ .carbon_g = 0, .nitrogen_g = 0, .phosphorus_g = 0 };
    const main_sapwood = state.branch_sapwood_carbon_g[main_branch];
    const total_sapwood = main_sapwood + other_sapwood;
    if (total_sapwood <= 0) return error.InvalidBranchSapwoodPool;
    const initial_total_reserve_c = state.branch_reserve_carbon_g[main_branch] + state.branch_reserve_carbon_g[other_branch];
    const carbon_difference = (state.branch_reserve_carbon_g[main_branch] * other_sapwood - state.branch_reserve_carbon_g[other_branch] * main_sapwood) / total_sapwood;
    const carbon_flux = carbon_exchange_per_h * carbon_difference * timestep_h;
    const main_c = state.branch_reserve_carbon_g[main_branch] - carbon_flux;
    const other_c = state.branch_reserve_carbon_g[other_branch] + carbon_flux;
    if (main_c < -1e-12 or other_c < -1e-12) return error.BranchReserveExchangeExhaustedCarbon;
    var nitrogen_flux: f64 = 0;
    var phosphorus_flux: f64 = 0;
    if (initial_total_reserve_c > presence_threshold_g_c) {
        const nitrogen_difference = (state.branch_reserve_nitrogen_g[main_branch] * other_c - state.branch_reserve_nitrogen_g[other_branch] * main_c) / initial_total_reserve_c;
        const phosphorus_difference = (state.branch_reserve_phosphorus_g[main_branch] * other_c - state.branch_reserve_phosphorus_g[other_branch] * main_c) / initial_total_reserve_c;
        nitrogen_flux = nutrient_exchange_per_h * nitrogen_difference * timestep_h;
        phosphorus_flux = nutrient_exchange_per_h * phosphorus_difference * timestep_h;
    }
    const main_n = state.branch_reserve_nitrogen_g[main_branch] - nitrogen_flux;
    const other_n = state.branch_reserve_nitrogen_g[other_branch] + nitrogen_flux;
    const main_p = state.branch_reserve_phosphorus_g[main_branch] - phosphorus_flux;
    const other_p = state.branch_reserve_phosphorus_g[other_branch] + phosphorus_flux;
    if (main_n < -1e-12 or other_n < -1e-12 or main_p < -1e-12 or other_p < -1e-12) return error.BranchReserveExchangeExhaustedNutrient;
    state.branch_reserve_carbon_g[main_branch] = @max(0.0, main_c);
    state.branch_reserve_carbon_g[other_branch] = @max(0.0, other_c);
    state.branch_reserve_nitrogen_g[main_branch] = @max(0.0, main_n);
    state.branch_reserve_nitrogen_g[other_branch] = @max(0.0, other_n);
    state.branch_reserve_phosphorus_g[main_branch] = @max(0.0, main_p);
    state.branch_reserve_phosphorus_g[other_branch] = @max(0.0, other_p);
    return .{ .carbon_g = carbon_flux, .nitrogen_g = nitrogen_flux, .phosphorus_g = phosphorus_flux };
}

pub fn fillGrainFromReserve(state: *State, branch: usize, grain_fill_started: bool, final_seed_count: f64, maximum_seed_carbon_g: f64, grain_fill_g_c_per_seed_h_25c: f64, growth_temperature_factor: f64, timestep_h: f64, minimum_grain_nutrient_fraction: f64, maximum_grain_nitrogen_to_carbon_g_per_g: f64, maximum_grain_phosphorus_to_carbon_g_per_g: f64, reserve_nitrogen_half_saturation_g_per_g_c: f64, reserve_phosphorus_half_saturation_g_per_g_c: f64, grain_precursor_growth: LeafGrowth) !GrainFillResult {
    inline for (.{ final_seed_count, maximum_seed_carbon_g, grain_fill_g_c_per_seed_h_25c, growth_temperature_factor, timestep_h, minimum_grain_nutrient_fraction, maximum_grain_nitrogen_to_carbon_g_per_g, maximum_grain_phosphorus_to_carbon_g_per_g, reserve_nitrogen_half_saturation_g_per_g_c, reserve_phosphorus_half_saturation_g_per_g_c, grain_precursor_growth.carbon_g, grain_precursor_growth.nitrogen_g, grain_precursor_growth.phosphorus_g }) |value| if (!std.math.isFinite(value)) return error.NonFiniteGrainFillInput;
    if (branch >= state.branch_grain_carbon_g.len) return error.CanopyBranchIndexOutOfBounds;
    if (final_seed_count < 0 or maximum_seed_carbon_g < 0 or grain_fill_g_c_per_seed_h_25c < 0 or growth_temperature_factor < 0 or timestep_h <= 0 or minimum_grain_nutrient_fraction < 0 or minimum_grain_nutrient_fraction > 1 or maximum_grain_nitrogen_to_carbon_g_per_g <= 0 or maximum_grain_phosphorus_to_carbon_g_per_g <= 0 or reserve_nitrogen_half_saturation_g_per_g_c < 0 or reserve_phosphorus_half_saturation_g_per_g_c < 0 or grain_precursor_growth.carbon_g < 0 or grain_precursor_growth.nitrogen_g < 0 or grain_precursor_growth.phosphorus_g < 0) return error.InvalidGrainFillInput;
    if (!grain_fill_started) return .{ .carbon_translocated_g = 0, .nitrogen_translocated_g = 0, .phosphorus_translocated_g = 0, .maximum_carbon_translocation_g = 0 };
    const grain_c = state.branch_grain_carbon_g[branch];
    const grain_n = state.branch_grain_nitrogen_g[branch];
    const grain_p = state.branch_grain_phosphorus_g[branch];
    const reserve_c = state.branch_reserve_carbon_g[branch];
    const reserve_n = state.branch_reserve_nitrogen_g[branch];
    const reserve_p = state.branch_reserve_phosphorus_g[branch];
    inline for (.{ grain_c, grain_n, grain_p, reserve_c, reserve_n, reserve_p }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidGrainFillState;
    const sink_capacity = maximum_seed_carbon_g * final_seed_count;
    const maximum_fill = if (grain_c >= sink_capacity) 0 else @max(0.0, grain_fill_g_c_per_seed_h_25c * final_seed_count * @sqrt(growth_temperature_factor) * timestep_h);
    const nutrient_deficient = grain_n < minimum_grain_nutrient_fraction * maximum_grain_nitrogen_to_carbon_g_per_g * grain_c or grain_p < minimum_grain_nutrient_fraction * maximum_grain_phosphorus_to_carbon_g_per_g * grain_c;
    const actual_fill_ceiling = if (nutrient_deficient) 0 else maximum_fill;
    const maximum_carbon_translocation = @min(maximum_fill, reserve_c);
    const carbon_translocated = @min(actual_fill_ceiling, reserve_c);
    const responsive_fraction = 1.0 - minimum_grain_nutrient_fraction;
    var nitrogen_translocated: f64 = 0;
    if (reserve_n > 0) {
        const reserve_constraint = reserve_n / (reserve_n + reserve_nitrogen_half_saturation_g_per_g_c * reserve_c);
        const grain_ratio = minimum_grain_nutrient_fraction + responsive_fraction * std.math.clamp(reserve_constraint, 0, 1);
        // GROSUB 5181--5183 is a three-argument AMIN1.  The reserve term is
        // bounded below independently; the grain deficit is deliberately not.
        nitrogen_translocated = @min(maximum_carbon_translocation * maximum_grain_nitrogen_to_carbon_g_per_g, @min(@max(0.0, reserve_n * grain_ratio), (grain_c + carbon_translocated) * maximum_grain_nitrogen_to_carbon_g_per_g - grain_n));
    }
    var phosphorus_translocated: f64 = 0;
    if (reserve_p > 0) {
        const reserve_constraint = reserve_p / (reserve_p + reserve_phosphorus_half_saturation_g_per_g_c * reserve_c);
        const grain_ratio = minimum_grain_nutrient_fraction + responsive_fraction * std.math.clamp(reserve_constraint, 0, 1);
        // GROSUB 5191--5193 has the same source-ordered three-way minimum.
        phosphorus_translocated = @min(maximum_carbon_translocation * maximum_grain_phosphorus_to_carbon_g_per_g, @min(@max(0.0, reserve_p * grain_ratio), (grain_c + carbon_translocated) * maximum_grain_phosphorus_to_carbon_g_per_g - grain_p));
    }
    nitrogen_translocated = @min(nitrogen_translocated, phosphorus_translocated * maximum_grain_nitrogen_to_carbon_g_per_g / maximum_grain_phosphorus_to_carbon_g_per_g);
    phosphorus_translocated = @min(phosphorus_translocated, nitrogen_translocated * maximum_grain_phosphorus_to_carbon_g_per_g / maximum_grain_nitrogen_to_carbon_g_per_g);
    const next_reserve_c = reserve_c + grain_precursor_growth.carbon_g - carbon_translocated;
    const next_reserve_n = reserve_n + grain_precursor_growth.nitrogen_g - nitrogen_translocated;
    const next_reserve_p = reserve_p + grain_precursor_growth.phosphorus_g - phosphorus_translocated;
    if (next_reserve_c < -1e-12 or next_reserve_n < -1e-12 or next_reserve_p < -1e-12) {
        std.log.err("grain fill exhausted reserve: branch={d} reserve_c_g={e} reserve_n_g={e} reserve_p_g={e}", .{ branch, next_reserve_c, next_reserve_n, next_reserve_p });
        return error.GrainFillExhaustedReserve;
    }
    state.branch_reserve_carbon_g[branch] = @max(0.0, next_reserve_c);
    state.branch_reserve_nitrogen_g[branch] = @max(0.0, next_reserve_n);
    state.branch_reserve_phosphorus_g[branch] = @max(0.0, next_reserve_p);
    state.branch_grain_carbon_g[branch] = grain_c + carbon_translocated;
    state.branch_grain_nitrogen_g[branch] = grain_n + nitrogen_translocated;
    state.branch_grain_phosphorus_g[branch] = grain_p + phosphorus_translocated;
    return .{ .carbon_translocated_g = carbon_translocated, .nitrogen_translocated_g = nitrogen_translocated, .phosphorus_translocated_g = phosphorus_translocated, .maximum_carbon_translocation_g = maximum_carbon_translocation };
}

pub fn senesceLeafAndSheathNode(state: *State, branch: usize, node_within_branch: usize, requested_fraction: f64, recycling: RecyclingFractions, protein_per_nitrogen_g_per_g_n: f64, protein_per_phosphorus_g_per_g_p: f64, woody_fraction: [2]f64, leaf_woody_nitrogen_fraction: [2]f64, sheath_woody_nitrogen_fraction: [2]f64, leaf_woody_phosphorus_fraction: [2]f64, sheath_woody_phosphorus_fraction: [2]f64, woody_kinetics: KineticFractions, leaf_kinetics: KineticFractions, sheath_kinetics: KineticFractions) !SenescenceProducts {
    inline for (.{ requested_fraction, recycling.carbon, recycling.nitrogen, recycling.phosphorus, protein_per_nitrogen_g_per_g_n, protein_per_phosphorus_g_per_g_p }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSenescenceInput;
    if (requested_fraction < 0 or recycling.carbon < 0 or recycling.carbon > 1 or recycling.nitrogen < 0 or recycling.nitrogen > 1 or recycling.phosphorus < 0 or recycling.phosphorus > 1 or protein_per_nitrogen_g_per_g_n < 0 or protein_per_phosphorus_g_per_g_p < 0) return error.InvalidSenescenceInput;
    inline for (.{ woody_fraction, leaf_woody_nitrogen_fraction, sheath_woody_nitrogen_fraction, leaf_woody_phosphorus_fraction, sheath_woody_phosphorus_fraction }) |fractions| for (fractions) |fraction| if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1 or @abs(fractions[0] + fractions[1] - 1.0) > 1e-8) return error.InvalidWoodyFraction;
    try woody_kinetics.validate();
    try leaf_kinetics.validate();
    try sheath_kinetics.validate();
    const nodes = try state.nodeRange(branch);
    if (node_within_branch >= nodes.end - nodes.first) return error.CanopyNodeIndexOutOfBounds;
    const node = nodes.first + node_within_branch;
    const fraction = @min(1.0, requested_fraction);
    var products: SenescenceProducts = .{};

    const leaf_c = state.node_leaf_carbon_g[node];
    const leaf_n = state.node_leaf_nitrogen_g[node];
    const leaf_p = state.node_leaf_phosphorus_g[node];
    const recycled_leaf_c = leaf_c * recycling.carbon;
    const recycled_leaf_n = leaf_n * (recycling.nitrogen + (1.0 - recycling.nitrogen) * recycling.carbon);
    const recycled_leaf_p = leaf_p * (recycling.phosphorus + (1.0 - recycling.phosphorus) * recycling.carbon);
    const sheath_c = state.node_sheath_carbon_g[node];
    const sheath_n = state.node_sheath_nitrogen_g[node];
    const sheath_p = state.node_sheath_phosphorus_g[node];
    const recycled_sheath_c = sheath_c * recycling.carbon;
    const recycled_sheath_n = if (sheath_c > 0) sheath_n * (recycling.nitrogen + (1.0 - recycling.nitrogen) * recycled_sheath_c / sheath_c) else 0;
    const recycled_sheath_p = if (sheath_c > 0) sheath_p * (recycling.phosphorus + (1.0 - recycling.phosphorus) * recycled_sheath_c / sheath_c) else 0;

    for (0..4) |kinetic| {
        products.woody_carbon_g[kinetic] = woody_kinetics.carbon[kinetic] * fraction * (leaf_c + sheath_c) * woody_fraction[0];
        products.woody_nitrogen_g[kinetic] = woody_kinetics.nitrogen[kinetic] * fraction * (leaf_n * leaf_woody_nitrogen_fraction[0] + sheath_n * sheath_woody_nitrogen_fraction[0]);
        products.woody_phosphorus_g[kinetic] = woody_kinetics.phosphorus[kinetic] * fraction * (leaf_p * leaf_woody_phosphorus_fraction[0] + sheath_p * sheath_woody_phosphorus_fraction[0]);
        products.nonwoody_carbon_g[kinetic] = fraction * woody_fraction[1] * (leaf_kinetics.carbon[kinetic] * (leaf_c - recycled_leaf_c) + sheath_kinetics.carbon[kinetic] * (sheath_c - recycled_sheath_c));
        products.nonwoody_nitrogen_g[kinetic] = fraction * (leaf_woody_nitrogen_fraction[1] * leaf_kinetics.nitrogen[kinetic] * (leaf_n - recycled_leaf_n) + sheath_woody_nitrogen_fraction[1] * sheath_kinetics.nitrogen[kinetic] * (sheath_n - recycled_sheath_n));
        products.nonwoody_phosphorus_g[kinetic] = fraction * (leaf_woody_phosphorus_fraction[1] * leaf_kinetics.phosphorus[kinetic] * (leaf_p - recycled_leaf_p) + sheath_woody_phosphorus_fraction[1] * sheath_kinetics.phosphorus[kinetic] * (sheath_p - recycled_sheath_p));
    }
    products.recycled_carbon_g = fraction * woody_fraction[1] * (recycled_leaf_c + recycled_sheath_c);
    products.recycled_nitrogen_g = fraction * (leaf_woody_nitrogen_fraction[1] * recycled_leaf_n + sheath_woody_nitrogen_fraction[1] * recycled_sheath_n);
    products.recycled_phosphorus_g = fraction * (leaf_woody_phosphorus_fraction[1] * recycled_leaf_p + sheath_woody_phosphorus_fraction[1] * recycled_sheath_p);

    const area_removed = fraction * state.node_leaf_area_m2[node];
    state.node_leaf_area_m2[node] -= area_removed;
    state.branch_leaf_area_m2[branch] = @max(0.0, state.branch_leaf_area_m2[branch] - area_removed);
    state.node_leaf_carbon_g[node] = @max(0.0, leaf_c * (1.0 - fraction));
    state.node_leaf_nitrogen_g[node] = @max(0.0, leaf_n * (1.0 - fraction));
    state.node_leaf_phosphorus_g[node] = @max(0.0, leaf_p * (1.0 - fraction));
    state.node_leaf_protein_g[node] = @max(0.0, state.node_leaf_protein_g[node] - fraction * @max(leaf_n * protein_per_nitrogen_g_per_g_n, leaf_p * protein_per_phosphorus_g_per_g_p));
    state.node_sheath_height_m[node] *= 1.0 - fraction;
    state.node_sheath_carbon_g[node] = @max(0.0, sheath_c * (1.0 - fraction));
    state.node_sheath_nitrogen_g[node] = @max(0.0, sheath_n * (1.0 - fraction));
    state.node_sheath_phosphorus_g[node] = @max(0.0, sheath_p * (1.0 - fraction));
    state.node_sheath_protein_g[node] = @max(0.0, state.node_sheath_protein_g[node] - fraction * @max(sheath_n * protein_per_nitrogen_g_per_g_n, sheath_p * protein_per_phosphorus_g_per_g_p));
    state.branch_leaf_carbon_g[branch] = @max(0.0, state.branch_leaf_carbon_g[branch] - fraction * leaf_c);
    state.branch_leaf_nitrogen_g[branch] = @max(0.0, state.branch_leaf_nitrogen_g[branch] - fraction * leaf_n);
    state.branch_leaf_phosphorus_g[branch] = @max(0.0, state.branch_leaf_phosphorus_g[branch] - fraction * leaf_p);
    state.branch_sheath_carbon_g[branch] = @max(0.0, state.branch_sheath_carbon_g[branch] - fraction * sheath_c);
    state.branch_sheath_nitrogen_g[branch] = @max(0.0, state.branch_sheath_nitrogen_g[branch] - fraction * sheath_n);
    state.branch_sheath_phosphorus_g[branch] = @max(0.0, state.branch_sheath_phosphorus_g[branch] - fraction * sheath_p);
    state.branch_mobile_carbon_g[branch] += products.recycled_carbon_g;
    state.branch_mobile_nitrogen_g[branch] += products.recycled_nitrogen_g;
    state.branch_mobile_phosphorus_g[branch] += products.recycled_phosphorus_g;
    state.branch_senescing_stalk_carbon_g[branch] += state.node_internode_carbon_g[node];
    state.branch_senescing_stalk_nitrogen_g[branch] += state.node_internode_nitrogen_g[node];
    state.branch_senescing_stalk_phosphorus_g[branch] += state.node_internode_phosphorus_g[node];
    state.node_internode_carbon_g[node] = 0;
    state.node_internode_nitrogen_g[node] = 0;
    state.node_internode_phosphorus_g[node] = 0;
    state.node_internode_length_m[node] = 0;
    return products;
}

fn checkedSum(counts: []const usize) !usize {
    var total: usize = 0;
    for (counts) |count| total = try std.math.add(usize, total, count);
    return total;
}

fn makeOffsets(allocator: std.mem.Allocator, counts: []const usize) ![]usize {
    const offsets = try allocator.alloc(usize, try std.math.add(usize, counts.len, 1));
    errdefer allocator.free(offsets);
    offsets[0] = 0;
    for (counts, 0..) |count, index| offsets[index + 1] = try std.math.add(usize, offsets[index], count);
    return offsets;
}

fn allocateZeroed(allocator: std.mem.Allocator, count: usize) ![]f64 {
    const values = try allocator.alloc(f64, count);
    @memset(values, 0);
    return values;
}

fn domainCount(comptime field_name: []const u8, plant_count: usize, plant_kinetic_count: usize, branch_count: usize, node_count: usize, sample_count: usize) usize {
    if (comptime std.mem.indexOf(u8, field_name, "_by_kinetic_") != null) return plant_kinetic_count;
    if (comptime std.mem.indexOf(u8, field_name, "_by_species_") != null and std.mem.startsWith(u8, field_name, "branch_")) return branch_count * 8;
    if (comptime std.mem.startsWith(u8, field_name, "plant_")) return plant_count;
    if (comptime std.mem.startsWith(u8, field_name, "branch_")) return branch_count;
    if (comptime std.mem.startsWith(u8, field_name, "node_")) return node_count;
    if (comptime std.mem.startsWith(u8, field_name, "sample_")) return sample_count;
    @compileError("runtime canopy field must use a domain prefix: " ++ field_name);
}

fn freeAllocatedF64Fields(state: *State, allocated_count: usize) void {
    var visited: usize = 0;
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        if (visited < allocated_count) state.allocator.free(@field(state, field.name));
        visited += 1;
    };
}

fn copyAroundInsertion(destination: []f64, source: []const f64, insertion_index: usize, inserted_count: usize) void {
    @memcpy(destination[0..insertion_index], source[0..insertion_index]);
    @memcpy(destination[insertion_index + inserted_count ..], source[insertion_index..]);
}

fn copyRemovingRange(destination: []f64, source: []const f64, first: usize, end: usize) void {
    std.debug.assert(first <= end and end <= source.len and destination.len == source.len - (end - first));
    @memcpy(destination[0..first], source[0..first]);
    @memcpy(destination[first..], source[end..]);
}

test "GROSUB absent leaf routes all C4 carbon without consuming sheath" {
    var state = try State.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{0});
    defer state.deinit();
    state.node_leaf_carbon_g[0] = 1.0e-12;
    state.branch_leaf_carbon_g[0] = 1.0e-12;
    state.node_sheath_carbon_g[0] = 2;
    state.branch_sheath_carbon_g[0] = 2;
    state.node_c3_nonstructural_carbon_g[0] = 0.4;
    state.node_c4_mesophyll_nonstructural_carbon_g[0] = 0.6;
    const allocation = try allocateNodeSenescenceDemandWithThreshold(1, 1.0e-12, 2, 0.5, 0.25, 0.8, 1.0e-9);
    try std.testing.expect(!allocation.leaf_present);
    try std.testing.expectEqual(@as(f64, 0), allocation.sheath_fraction);
    const kinetics: KineticFractions = .{ .carbon = @splat(0.25), .nitrogen = @splat(0.25), .phosphorus = @splat(0.25) };
    const products = try commitNodeSenescenceDemand(&state, 0, 0, allocation, .{ .carbon = 0.5, .nitrogen = 0.6, .phosphorus = 0.7 }, 2, 20, .{ 0.25, 0.75 }, .{ 0.25, 0.75 }, .{ 0.25, 0.75 }, .{ 0.25, 0.75 }, .{ 0.25, 0.75 }, kinetics, kinetics, kinetics);
    try std.testing.expectEqual(@as(f64, 0), state.node_c3_nonstructural_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.node_c4_mesophyll_nonstructural_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.node_leaf_carbon_g[0]);
    var litter_carbon_g_c: f64 = 0;
    for (0..4) |kinetic| litter_carbon_g_c += products.woody_carbon_g[kinetic] + products.nonwoody_carbon_g[kinetic];
    try std.testing.expectApproxEqAbs(@as(f64, 1.000000000001), litter_carbon_g_c, 1.0e-15);
    try std.testing.expectEqual(@as(f64, 2), state.node_sheath_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 2), state.branch_sheath_carbon_g[0]);
}

test "GROSUB descending stalk sweep scales sapwood recycling only once" {
    var state = try State.init(std.testing.allocator, 1, 1, &.{1}, &.{2}, &.{ 0, 0 });
    defer state.deinit();
    state.branch_stalk_carbon_g[0] = 20;
    state.branch_stalk_nitrogen_g[0] = 2;
    state.branch_stalk_phosphorus_g[0] = 0.2;
    state.branch_sapwood_carbon_g[0] = 5;
    @memset(state.node_internode_carbon_g, 4);
    @memset(state.node_internode_nitrogen_g, 0.4);
    @memset(state.node_internode_phosphorus_g, 0.04);
    @memset(state.node_internode_length_m, 1);
    state.node_height_m[0] = 1;
    state.node_height_m[1] = 2;
    const setup = (try perennial_stalk_senescence_setup.prepare(.{
        .is_perennial = true,
        .excess_maintenance_respiration_g_c_per_timestep = 1,
        .stalk_carbon_g_c = 20,
        .sapwood_carbon_g_c = 5,
        .presence_threshold_g_c = 0,
        .first_internode = 0,
        .last_internode = 1,
        .shoot_recycling = .{ .carbon = 1, .nitrogen = 1, .phosphorus = 1 },
    })).?;
    const kinetics: KineticFractions = .{ .carbon = @splat(0.25), .nitrogen = @splat(0.25), .phosphorus = @splat(0.25) };
    const first = try commitInternodeSenescenceDemandScaled(&state, 0, 1, 1, 0, setup.sapwood_recycling, 0, .{ 0.2, 0.8 }, .{ 0.2, 0.8 }, .{ 0.2, 0.8 }, kinetics, kinetics);
    const second = try commitInternodeSenescenceDemandScaled(&state, 0, 0, first.remaining_respiration_demand_g_c, 0, setup.sapwood_recycling, 0, .{ 0.2, 0.8 }, .{ 0.2, 0.8 }, .{ 0.2, 0.8 }, kinetics, kinetics);
    try std.testing.expectEqual(@as(f64, 1), first.fraction);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), second.fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 3.2), state.node_internode_carbon_g[0], 1.0e-15);
}

test {
    // Keeps the extracted tests discoverable by `zig build test`,
    // which only reaches files reachable by import.
    _ = @import("../../validation/canopy_photosynthesis_test.zig");
}
