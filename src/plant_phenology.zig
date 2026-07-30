const std = @import("std");
const CellRange = @import("compute.zig").CellRange;
const PlantState = @import("grid.zig").PlantState;
const PlantRootState = @import("plant_root_system.zig").State;
const execution_calendar_date = @import("execution_calendar_date.zig");

pub const Parameters = struct {
    gas_constant_j_per_mol_k: f64,
    temperature_scale_k: f64,
    arrhenius_log_prefactor: f64,
    activation_energy_j_per_mol: f64,
    low_temperature_inactivation_j_per_mol: f64,
    high_temperature_inactivation_j_per_mol: f64,
    minimum_turgor_potential_mpa: f64,
    oxygen_stress_exponent: f64,
    vegetative_stage_duration: f64,
    reproductive_stage_duration: f64,
    drought_leafout_total_water_potential_mpa: f64,
    nonvascular_leafoff_total_water_potential_mpa: f64,
    vascular_leafoff_total_water_potential_mpa: f64,
    maximum_photoperiod_counter_h: f64,
    emergence_area_threshold_m2_per_plant: f64,
    emergence_root_depth_margin_m: f64,

    pub fn validate(self: Parameters) !void {
        inline for (@typeInfo(Parameters).@"struct".fields) |field| if (!std.math.isFinite(@field(self, field.name))) return error.NonFinitePhenologyParameter;
        if (self.gas_constant_j_per_mol_k <= 0 or self.temperature_scale_k <= 0 or self.activation_energy_j_per_mol <= 0 or self.low_temperature_inactivation_j_per_mol <= 0 or self.high_temperature_inactivation_j_per_mol <= 0 or self.oxygen_stress_exponent <= 0 or self.vegetative_stage_duration <= 0 or self.reproductive_stage_duration <= 0 or self.maximum_photoperiod_counter_h <= 0 or self.emergence_area_threshold_m2_per_plant < 0 or self.emergence_root_depth_margin_m < 0) return error.InvalidPhenologyParameter;
    }
};

/// Runtime defaults reproduce HFUNC's source coefficients without making
/// them compile-time scientific controls.
pub fn compatibilityParameters() Parameters {
    return .{
        .gas_constant_j_per_mol_k = 8.3143,
        .temperature_scale_k = 710,
        .arrhenius_log_prefactor = 24.269,
        .activation_energy_j_per_mol = 60_000,
        .low_temperature_inactivation_j_per_mol = 197_500,
        .high_temperature_inactivation_j_per_mol = 218_500,
        .minimum_turgor_potential_mpa = 0.1,
        .oxygen_stress_exponent = 0.25,
        .vegetative_stage_duration = 2.0,
        .reproductive_stage_duration = 0.667,
        .drought_leafout_total_water_potential_mpa = -0.1,
        .nonvascular_leafoff_total_water_potential_mpa = -150,
        .vascular_leafoff_total_water_potential_mpa = -1.5,
        .maximum_photoperiod_counter_h = 3600,
        .emergence_area_threshold_m2_per_plant = 1.0e-6,
        .emergence_root_depth_margin_m = 1.0e-6,
    };
}

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    species_count: usize,
    active: []bool,
    emerged: []bool,
    lifecycle_initialized: []bool,
    reseed_pending: []bool,
    death_replant_pending: []bool,
    leafout_transition_this_step: []bool,
    replant_day_of_year: []u16,
    replant_year: []u16,
    annual_growth_habit: []bool,
    floral_initiated: []bool,
    node_initiation_rate_at_25c_per_h: []f64,
    leaf_appearance_rate_at_25c_per_h: []f64,
    initiated_node_count: []f64,
    appeared_leaf_count: []f64,
    node_initiation_per_timestep: []f64,
    leaf_appearance_per_timestep: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, species_count: usize) !State {
        if (cell_count == 0 or species_count == 0) return error.InvalidPhenologyDimensions;
        const count = try std.math.mul(usize, cell_count, species_count);
        var result: State = undefined;
        result.allocator = allocator;
        result.cell_count = cell_count;
        result.species_count = species_count;
        result.active = try allocator.alloc(bool, count);
        errdefer allocator.free(result.active);
        result.emerged = try allocator.alloc(bool, count);
        errdefer allocator.free(result.emerged);
        result.lifecycle_initialized = try allocator.alloc(bool, count);
        errdefer allocator.free(result.lifecycle_initialized);
        result.reseed_pending = try allocator.alloc(bool, count);
        errdefer allocator.free(result.reseed_pending);
        result.death_replant_pending = try allocator.alloc(bool, count);
        errdefer allocator.free(result.death_replant_pending);
        result.leafout_transition_this_step = try allocator.alloc(bool, count);
        errdefer allocator.free(result.leafout_transition_this_step);
        result.replant_day_of_year = try allocator.alloc(u16, count);
        errdefer allocator.free(result.replant_day_of_year);
        result.replant_year = try allocator.alloc(u16, count);
        errdefer allocator.free(result.replant_year);
        result.annual_growth_habit = try allocator.alloc(bool, count);
        errdefer allocator.free(result.annual_growth_habit);
        result.floral_initiated = try allocator.alloc(bool, count);
        errdefer allocator.free(result.floral_initiated);
        var allocated: usize = 0;
        errdefer freeF64Allocated(&result, allocated);
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            @field(result, field.name) = try allocator.alloc(f64, count);
            @memset(@field(result, field.name), 0);
            allocated += 1;
        };
        @memset(result.active, false);
        @memset(result.emerged, false);
        @memset(result.lifecycle_initialized, false);
        @memset(result.reseed_pending, false);
        @memset(result.death_replant_pending, false);
        @memset(result.leafout_transition_this_step, false);
        @memset(result.replant_day_of_year, 0);
        @memset(result.replant_year, 0);
        @memset(result.annual_growth_habit, false);
        @memset(result.floral_initiated, false);
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.allocator.free(self.floral_initiated);
        self.allocator.free(self.annual_growth_habit);
        self.allocator.free(self.lifecycle_initialized);
        self.allocator.free(self.reseed_pending);
        self.allocator.free(self.death_replant_pending);
        self.allocator.free(self.leafout_transition_this_step);
        self.allocator.free(self.replant_year);
        self.allocator.free(self.replant_day_of_year);
        self.allocator.free(self.emerged);
        self.allocator.free(self.active);
        self.* = undefined;
    }

    pub fn validateFinite(self: State) !void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) for (@field(self, field.name), 0..) |value, index| if (!std.math.isFinite(value)) {
            std.log.err("non-finite plant phenology: field={s} index={d} value={e}", .{ field.name, index, value });
            return error.NonFinitePlantPhenology;
        };
        if (self.replant_day_of_year.len != self.replant_year.len)
            return error.PlantPhenologyReplantDimensionMismatch;
        for (self.replant_day_of_year, self.replant_year) |day, year| {
            if (day == 0 and year == 0) continue;
            if (day == 0 or year == 0)
                return error.InvalidPlantPhenologyReplantDate;
            const maximum_day: u16 =
                if (execution_calendar_date.isLeapYear(year)) 366 else 365;
            if (day > maximum_day)
                return error.InvalidPlantPhenologyReplantDate;
        }
    }
};

pub const BranchDevelopmentState = struct {
    allocator: std.mem.Allocator,
    branch_count: usize,
    maturity_group: []f64,
    initial_reproductive_stage: []f64,
    final_reproductive_stage: []f64,
    vegetative_stage_change: []f64,
    initial_total_reproductive_node_change: []f64,
    final_total_reproductive_node_change: []f64,
    hours_without_grain_fill: []f64,
    remobilization_progress_h: []f64,
    perennial_node_scaling: []f64,
    maximum_concurrently_growing_nodes: []usize,
    stage_day: []u32,
    leafout_initialization_enabled: []bool,
    dead: []bool,

    pub fn init(allocator: std.mem.Allocator, branch_count: usize) !BranchDevelopmentState {
        if (branch_count == 0) return error.InvalidBranchDevelopmentDimensions;
        var result: BranchDevelopmentState = undefined;
        result.allocator = allocator;
        result.branch_count = branch_count;
        var allocated_f64: usize = 0;
        errdefer inline for (@typeInfo(BranchDevelopmentState).@"struct".fields) |field| if (field.type == []f64 and allocated_f64 > 0) {
            allocated_f64 -= 1;
            allocator.free(@field(result, field.name));
        };
        inline for (@typeInfo(BranchDevelopmentState).@"struct".fields) |field| if (field.type == []f64) {
            @field(result, field.name) = try allocator.alloc(f64, branch_count);
            @memset(@field(result, field.name), 0);
            allocated_f64 += 1;
        };
        result.stage_day = try allocator.alloc(u32, try std.math.mul(usize, branch_count, 10));
        errdefer allocator.free(result.stage_day);
        @memset(result.stage_day, 0);
        result.maximum_concurrently_growing_nodes = try allocator.alloc(usize, branch_count);
        errdefer allocator.free(result.maximum_concurrently_growing_nodes);
        @memset(result.maximum_concurrently_growing_nodes, 0);
        result.leafout_initialization_enabled = try allocator.alloc(bool, branch_count);
        errdefer allocator.free(result.leafout_initialization_enabled);
        @memset(result.leafout_initialization_enabled, false);
        result.dead = try allocator.alloc(bool, branch_count);
        @memset(result.dead, false);
        return result;
    }

    pub fn deinit(self: *BranchDevelopmentState) void {
        self.allocator.free(self.dead);
        self.allocator.free(self.leafout_initialization_enabled);
        self.allocator.free(self.maximum_concurrently_growing_nodes);
        self.allocator.free(self.stage_day);
        inline for (@typeInfo(BranchDevelopmentState).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }

    pub fn validateFinite(self: BranchDevelopmentState) !void {
        inline for (@typeInfo(BranchDevelopmentState).@"struct".fields) |field| if (field.type == []f64) for (@field(self, field.name), 0..) |value, index| if (!std.math.isFinite(value)) {
            std.log.err("non-finite branch development state: field={s} branch={d} value={e}", .{ field.name, index, value });
            return error.NonFiniteBranchDevelopmentState;
        };
    }

    pub fn clone(self: BranchDevelopmentState) !BranchDevelopmentState {
        const result = try BranchDevelopmentState.init(self.allocator, self.branch_count);
        inline for (@typeInfo(BranchDevelopmentState).@"struct".fields) |field| switch (field.type) {
            []f64, []usize, []u32, []bool => @memcpy(@field(result, field.name), @field(self, field.name)),
            else => {},
        };
        return result;
    }

    pub fn insertBranch(self: *BranchDevelopmentState, index: usize) !void {
        if (index > self.branch_count) return error.BranchDevelopmentIndexOutOfBounds;
        var replacement = try BranchDevelopmentState.init(self.allocator, try std.math.add(usize, self.branch_count, 1));
        errdefer replacement.deinit();
        inline for (@typeInfo(BranchDevelopmentState).@"struct".fields) |field| switch (field.type) {
            []f64, []usize, []bool => copySliceAroundInsertion(@field(replacement, field.name), @field(self, field.name), index, 1),
            []u32 => copySliceAroundInsertion(@field(replacement, field.name), @field(self, field.name), index * 10, 10),
            else => {},
        };
        var previous = self.*;
        self.* = replacement;
        previous.deinit();
    }

    pub fn initializeRange(self: *BranchDevelopmentState, first: usize, end: usize, maturity_group: f64, seed_initial_stage: f64, leafout_enabled: bool, perennial_node_scaling: f64, maximum_concurrently_growing_nodes: usize) !void {
        if (first > end or end > self.branch_count) return error.InvalidBranchDevelopmentRange;
        if (!std.math.isFinite(maturity_group) or maturity_group < 0 or !std.math.isFinite(seed_initial_stage) or seed_initial_stage < 0 or !std.math.isFinite(perennial_node_scaling) or perennial_node_scaling < 1 or maximum_concurrently_growing_nodes == 0) return error.InvalidBranchDevelopmentInitialization;
        for (first..end) |branch| {
            self.maturity_group[branch] = maturity_group;
            self.initial_reproductive_stage[branch] = seed_initial_stage;
            self.final_reproductive_stage[branch] = 0;
            self.leafout_initialization_enabled[branch] = leafout_enabled;
            self.perennial_node_scaling[branch] = perennial_node_scaling;
            self.maximum_concurrently_growing_nodes[branch] = maximum_concurrently_growing_nodes;
        }
    }

    pub fn clearRangeForReconstruction(self: *BranchDevelopmentState, first: usize, end: usize) !void {
        if (first > end or end > self.branch_count) return error.InvalidBranchDevelopmentRange;
        inline for (@typeInfo(BranchDevelopmentState).@"struct".fields) |field| switch (field.type) {
            []f64 => @memset(@field(self, field.name)[first..end], 0),
            []usize => @memset(@field(self, field.name)[first..end], 0),
            []bool => @memset(@field(self, field.name)[first..end], false),
            []u32 => @memset(@field(self, field.name)[first * 10 .. end * 10], 0),
            else => {},
        };
    }

    pub fn removeRange(self: *BranchDevelopmentState, first: usize, end: usize) !void {
        if (first > end or end > self.branch_count) return error.InvalidBranchDevelopmentRange;
        if (first == end) return;
        var replacement = try BranchDevelopmentState.init(self.allocator, self.branch_count - (end - first));
        errdefer replacement.deinit();
        inline for (@typeInfo(BranchDevelopmentState).@"struct".fields) |field| switch (field.type) {
            []f64, []usize, []bool => {
                @memcpy(@field(replacement, field.name)[0..first], @field(self, field.name)[0..first]);
                @memcpy(@field(replacement, field.name)[first..], @field(self, field.name)[end..]);
            },
            []u32 => {
                @memcpy(@field(replacement, field.name)[0 .. first * 10], @field(self, field.name)[0 .. first * 10]);
                @memcpy(@field(replacement, field.name)[first * 10 ..], @field(self, field.name)[end * 10 ..]);
            },
            else => {},
        };
        var previous = self.*;
        self.* = replacement;
        previous.deinit();
    }
};

pub const PostCutResetRequest = struct {
    selected_branch: usize,
    main_branch: usize,
    plant_branch_first: usize,
    plant_branch_end: usize,
    biomass_turnover_type: u8,
    root_profile_type: u8,
    grazing: bool,
    canopy_height_m: f64,
    cutting_height_m: f64,
    winter_phenology_type: u8,
    accumulated_vernalization_h: f64,
    required_vernalization_h: f64,
    vernalization_reset_fraction: f64,
    initial_maturity_group: f64,
    current_reproductive_stage_by_branch: []const f64,
    current_day: u32,
};

test "STARTQ branch development initializes runtime ranges and concurrent-node controls" {
    var state = try BranchDevelopmentState.init(std.testing.allocator, 7);
    defer state.deinit();
    try state.initializeRange(2, 7, 8, 1.5, true, 2, 24);
    for (2..7) |branch| {
        try std.testing.expectEqual(@as(f64, 8), state.maturity_group[branch]);
        try std.testing.expectEqual(@as(f64, 1.5), state.initial_reproductive_stage[branch]);
        try std.testing.expectEqual(@as(f64, 0), state.final_reproductive_stage[branch]);
        try std.testing.expect(state.leafout_initialization_enabled[branch]);
        try std.testing.expectEqual(@as(f64, 2), state.perennial_node_scaling[branch]);
        try std.testing.expectEqual(@as(usize, 24), state.maximum_concurrently_growing_nodes[branch]);
    }
}

pub fn resetDevelopmentAfterCut(state: *BranchDevelopmentState, request: PostCutResetRequest) !bool {
    inline for (.{ request.canopy_height_m, request.cutting_height_m, request.accumulated_vernalization_h, request.required_vernalization_h, request.vernalization_reset_fraction, request.initial_maturity_group }) |value| if (!std.math.isFinite(value)) return error.NonFinitePostCutResetInput;
    if (request.plant_branch_first > request.plant_branch_end or request.plant_branch_end > state.branch_count or request.current_reproductive_stage_by_branch.len < request.plant_branch_end or request.selected_branch < request.plant_branch_first or request.selected_branch >= request.plant_branch_end or request.main_branch < request.plant_branch_first or request.main_branch >= request.plant_branch_end or request.canopy_height_m < 0 or request.cutting_height_m < 0 or request.accumulated_vernalization_h < 0 or request.required_vernalization_h < 0 or request.vernalization_reset_fraction < 0 or request.vernalization_reset_fraction > 1) return error.InvalidPostCutResetInput;
    for (request.current_reproductive_stage_by_branch[request.plant_branch_first..request.plant_branch_end]) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePostCutResetInput;
    const herbaceous = request.biomass_turnover_type == 0 or request.root_profile_type <= 1;
    const stage_one_day = state.stage_day[request.selected_branch * 10];
    const developmental_window = (request.winter_phenology_type != 0 and request.accumulated_vernalization_h <= request.vernalization_reset_fraction * request.required_vernalization_h) or (request.winter_phenology_type == 0 and stage_one_day != 0);
    if (!herbaceous or request.grazing or request.canopy_height_m <= request.cutting_height_m or !developmental_window) return false;
    const first = if (request.selected_branch == request.main_branch) request.plant_branch_first else request.selected_branch;
    const end = if (request.selected_branch == request.main_branch) request.plant_branch_end else request.selected_branch + 1;
    for (first..end) |branch| {
        state.maturity_group[branch] = request.initial_maturity_group;
        state.initial_reproductive_stage[branch] = request.current_reproductive_stage_by_branch[branch];
        state.final_reproductive_stage[branch] = 0;
        state.vegetative_stage_change[branch] = 0;
        state.initial_total_reproductive_node_change[branch] = 0;
        state.final_total_reproductive_node_change[branch] = 0;
        state.hours_without_grain_fill[branch] = 0;
        state.stage_day[branch * 10] = request.current_day;
        @memset(state.stage_day[branch * 10 + 1 .. branch * 10 + 10], 0);
        state.leafout_initialization_enabled[branch] = false;
    }
    return true;
}

pub fn terminatePlantBranches(state: *BranchDevelopmentState, plant_branch_first: usize, plant_branch_end: usize, terminate: bool) !void {
    if (plant_branch_first > plant_branch_end or plant_branch_end > state.branch_count) return error.InvalidBranchDevelopmentRange;
    if (terminate) @memset(state.dead[plant_branch_first..plant_branch_end], true);
}

pub const AdvanceContext = struct {
    phenology: *State,
    plants: *const PlantState,
    thermal_acclimation_offset_k: []const f64,
    canopy_turgor_potential_mpa: []const f64,
    root_oxygen_uptake_to_demand_fraction: []const f64,
    emerged: []const bool,
    timestep_hours: f64,
    parameters: Parameters,
};

/// Ports HFUNC node initiation and leaf appearance. Independent species slots
/// are safe for horizontal tile dispatch; branch phenology is tracked separately
/// in the ledger and is not silently represented by these species totals.
pub fn advanceTile(context: *AdvanceContext, range: CellRange) !void {
    try context.parameters.validate();
    const state = context.phenology;
    const count = try std.math.mul(usize, state.cell_count, state.species_count);
    if (range.end > state.cell_count or context.plants.cell_count != state.cell_count or context.plants.species_count != state.species_count or context.thermal_acclimation_offset_k.len != count or context.canopy_turgor_potential_mpa.len != count or context.root_oxygen_uptake_to_demand_fraction.len != count or context.emerged.len != count) return error.PhenologyDimensionMismatch;
    if (!std.math.isFinite(context.timestep_hours) or context.timestep_hours <= 0) return error.InvalidPhenologyTimestep;
    for (range.first..range.end) |cell| for (0..state.species_count) |species| {
        const index = cell * state.species_count + species;
        state.node_initiation_per_timestep[index] = 0;
        state.leaf_appearance_per_timestep[index] = 0;
        if (!state.active[index] or !context.emerged[index]) continue;
        const acclimated_temperature_k = context.plants.canopy_temperature_k[index] + context.thermal_acclimation_offset_k[index];
        const turgor = context.canopy_turgor_potential_mpa[index];
        const oxygen_ratio = context.root_oxygen_uptake_to_demand_fraction[index];
        if (!std.math.isFinite(acclimated_temperature_k) or acclimated_temperature_k <= 0 or !std.math.isFinite(turgor) or !std.math.isFinite(oxygen_ratio) or oxygen_ratio < 0) return error.InvalidPhenologyRuntimeInput;
        const temperature_function = arrheniusTemperatureFunction(acclimated_temperature_k, context.parameters);
        var node_rate = @max(0.0, state.node_initiation_rate_at_25c_per_h[index] * temperature_function * context.timestep_hours);
        var leaf_rate = @max(0.0, state.leaf_appearance_rate_at_25c_per_h[index] * temperature_function * context.timestep_hours);
        if (state.annual_growth_habit[index] and !state.floral_initiated[index]) {
            const water_stress = std.math.clamp(turgor - context.parameters.minimum_turgor_potential_mpa, 0, 1);
            const oxygen_stress = std.math.pow(f64, oxygen_ratio, context.parameters.oxygen_stress_exponent);
            node_rate *= water_stress * oxygen_stress;
            leaf_rate *= water_stress * oxygen_stress;
        }
        state.node_initiation_per_timestep[index] = node_rate;
        state.leaf_appearance_per_timestep[index] = leaf_rate;
        state.initiated_node_count[index] += node_rate;
        state.appeared_leaf_count[index] += leaf_rate;
    };
}

pub fn arrheniusTemperatureFunction(temperature_k: f64, parameters: Parameters) f64 {
    const rt = parameters.gas_constant_j_per_mol_k * temperature_k;
    const scaled_temperature = parameters.temperature_scale_k * temperature_k;
    const inactivation = 1.0 + @exp((parameters.low_temperature_inactivation_j_per_mol - scaled_temperature) / rt) + @exp((scaled_temperature - parameters.high_temperature_inactivation_j_per_mol) / rt);
    return @exp(parameters.arrhenius_log_prefactor - parameters.activation_energy_j_per_mol / rt) / inactivation;
}

pub const NonstructuralConcentrations = struct {
    carbon_g_per_g: f64,
    nitrogen_g_per_g: f64,
    phosphorus_g_per_g: f64,
    salt_mol_per_g_c: f64,
};

/// HFUNC concentration convention: each mobile pool is divided by structural
/// mass plus that same pool, rather than by a shared total dry mass.
pub fn nonstructuralConcentrations(structural_c_g: f64, mobile_c_g: f64, mobile_n_g: f64, mobile_p_g: f64, salt_mol: f64, dynamic_salts: bool, structural_presence_threshold_g: f64) !NonstructuralConcentrations {
    inline for (.{ structural_c_g, mobile_c_g, mobile_n_g, mobile_p_g, salt_mol, structural_presence_threshold_g }) |value| if (!std.math.isFinite(value)) return error.NonFinitePlantPoolInput;
    if (structural_c_g < 0 or mobile_c_g < 0 or mobile_n_g < 0 or mobile_p_g < 0 or salt_mol < 0 or structural_presence_threshold_g < 0) return error.InvalidPlantPoolInput;
    if (structural_c_g <= structural_presence_threshold_g) return .{ .carbon_g_per_g = 1, .nitrogen_g_per_g = 1, .phosphorus_g_per_g = 1, .salt_mol_per_g_c = 0 };
    return .{
        .carbon_g_per_g = @max(0.0, mobile_c_g / (structural_c_g + mobile_c_g)),
        .nitrogen_g_per_g = @max(0.0, mobile_n_g / (structural_c_g + mobile_n_g)),
        .phosphorus_g_per_g = @max(0.0, mobile_p_g / (structural_c_g + mobile_p_g)),
        .salt_mol_per_g_c = if (dynamic_salts) @max(0.0, salt_mol / (structural_c_g + salt_mol)) else 0,
    };
}

pub const ShootBranchDecision = struct {
    initiate: bool,
    branch_ordinal: usize,
    branch_order: usize,
    maturity_group: f64,
};

pub fn shootBranchDecision(initialization_complete: bool, crop_active: bool, plant_density_per_m2: f64, root_turgor_potential_mpa: f64, minimum_root_turgor_potential_mpa: f64, growth_habit: u8, leaf_phenology_type: u8, floral_initiation_day: usize, current_branch_count: usize, reserve_c_g: f64, canopy_nonstructural_c_g_per_g: f64, branching_threshold_g_per_g: f64, candidate_branch_ordinal: usize, main_branch_ordinal: usize, main_branch_node_count: f64, branch_count_ever_initiated: usize, maximum_concurrently_growing_nodes: usize, perennial_node_scaling: f64, seed_node_count: f64, initial_maturity_group: f64) !ShootBranchDecision {
    inline for (.{ plant_density_per_m2, root_turgor_potential_mpa, minimum_root_turgor_potential_mpa, reserve_c_g, canopy_nonstructural_c_g_per_g, branching_threshold_g_per_g, main_branch_node_count, perennial_node_scaling, seed_node_count, initial_maturity_group }) |value| if (!std.math.isFinite(value)) return error.NonFiniteShootBranchInput;
    if (plant_density_per_m2 < 0 or reserve_c_g < 0 or canopy_nonstructural_c_g_per_g < 0 or branching_threshold_g_per_g < 0 or perennial_node_scaling <= 0 or seed_node_count < 0) return error.InvalidShootBranchInput;
    const developmental_window = growth_habit != 0 or (growth_habit == 0 and leaf_phenology_type == 1) or floral_initiation_day == 0;
    const carbon_permits = (current_branch_count == 0 and reserve_c_g > 0) or (canopy_nonstructural_c_g_per_g > branching_threshold_g_per_g and branching_threshold_g_per_g > 0);
    const spacing_permits = candidate_branch_ordinal == main_branch_ordinal or main_branch_node_count > @as(f64, @floatFromInt(branch_count_ever_initiated)) + @as(f64, @floatFromInt(maximum_concurrently_growing_nodes)) / perennial_node_scaling + seed_node_count or (growth_habit == 0 and leaf_phenology_type == 1);
    // IFLGC is the crop-active flag in HFUNC. It is deliberately not the
    // emergence flag: the source may establish topology before shoot emergence.
    const initiate = !initialization_complete and crop_active and plant_density_per_m2 > 0 and root_turgor_potential_mpa > minimum_root_turgor_potential_mpa and developmental_window and carbon_permits and spacing_permits;
    const next_order = branch_count_ever_initiated + 1;
    return .{
        .initiate = initiate,
        .branch_ordinal = candidate_branch_ordinal,
        .branch_order = if (initiate) next_order - 1 else branch_count_ever_initiated,
        .maturity_group = if (initiate and growth_habit == 0) @max(0.0, initial_maturity_group - @as(f64, @floatFromInt(next_order - 1))) else initial_maturity_group,
    };
}

pub fn shouldAddRootAxis(root_turgor_potential_mpa: f64, minimum_root_turgor_potential_mpa: f64, current_root_axis_count: usize, main_branch_node_count: f64, perennial_node_scaling: f64, seed_node_count: f64, reserve_c_g: f64, canopy_nonstructural_c_g_per_g: f64, root_branching_threshold_g_per_g: f64, maximum_root_axis_count: usize) !bool {
    inline for (.{ root_turgor_potential_mpa, minimum_root_turgor_potential_mpa, main_branch_node_count, perennial_node_scaling, seed_node_count, reserve_c_g, canopy_nonstructural_c_g_per_g, root_branching_threshold_g_per_g }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootAxisInput;
    if (perennial_node_scaling <= 0 or seed_node_count < 0 or reserve_c_g < 0 or canopy_nonstructural_c_g_per_g < 0 or root_branching_threshold_g_per_g < 0 or maximum_root_axis_count == 0) return error.InvalidRootAxisInput;
    if (current_root_axis_count >= maximum_root_axis_count or root_turgor_potential_mpa <= minimum_root_turgor_potential_mpa) return false;
    const stage_permits = current_root_axis_count == 0 or main_branch_node_count > @as(f64, @floatFromInt(current_root_axis_count)) / perennial_node_scaling + seed_node_count;
    const carbon_permits = (current_root_axis_count == 0 and reserve_c_g > 0) or (canopy_nonstructural_c_g_per_g > root_branching_threshold_g_per_g and root_branching_threshold_g_per_g > 0);
    return stage_permits and carbon_permits;
}

pub fn activateNextRootAxis(roots: *PlantRootState, plant: usize, root_turgor_potential_mpa: f64, minimum_root_turgor_potential_mpa: f64, main_branch_node_count: f64, perennial_node_scaling: f64, seed_node_count: f64, reserve_c_g: f64, canopy_nonstructural_c_g_per_g: f64, root_branching_threshold_g_per_g: f64) !bool {
    if (plant >= roots.plant_count) return error.PlantRootIndexOutOfBounds;
    const current_count = roots.active_root_axis_count[plant];
    if (!try shouldAddRootAxis(root_turgor_potential_mpa, minimum_root_turgor_potential_mpa, current_count, main_branch_node_count, perennial_node_scaling, seed_node_count, reserve_c_g, canopy_nonstructural_c_g_per_g, root_branching_threshold_g_per_g, roots.root_axis_count)) return false;
    roots.active_root_axis_count[plant] = try std.math.add(usize, current_count, 1);
    roots.roots_dead[plant] = false;
    return true;
}

fn freeF64Allocated(state: *State, count: usize) void {
    var visited: usize = 0;
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        if (visited < count) state.allocator.free(@field(state, field.name));
        visited += 1;
    };
}

fn copySliceAroundInsertion(destination: anytype, source: anytype, insertion_index: usize, inserted_count: usize) void {
    @memcpy(destination[0..insertion_index], source[0..insertion_index]);
    @memcpy(destination[insertion_index + inserted_count ..], source[insertion_index..]);
}

test "phenology temperature function is normalized near 25 C" {
    const parameters = compatibilityParameters();
    try std.testing.expectApproxEqRel(@as(f64, 1), arrheniusTemperatureFunction(298.15, parameters), 0.01);
}

test "HFUNC root axes activate up to runtime capacity without ten-axis ceiling" {
    var roots = try PlantRootState.init(std.testing.allocator, 1, 2, 13);
    defer roots.deinit();
    for (0..13) |axis| {
        const changed = try activateNextRootAxis(&roots, 0, 0.2, 0.1, @floatFromInt(axis + 2), 1, 0, 1, 0.5, 0.1);
        try std.testing.expect(changed);
    }
    try std.testing.expectEqual(@as(usize, 13), roots.active_root_axis_count[0]);
    try std.testing.expect(!(try activateNextRootAxis(&roots, 0, 0.2, 0.1, 100, 1, 0, 1, 0.5, 0.1)));
}

test "GROSUB post-cut main branch reset propagates to all runtime branches" {
    var branches = try BranchDevelopmentState.init(std.testing.allocator, 3);
    defer branches.deinit();
    @memset(branches.maturity_group, 9);
    @memset(branches.final_reproductive_stage, 4);
    @memset(branches.hours_without_grain_fill, 100);
    @memset(branches.remobilization_progress_h, 48);
    @memset(branches.leafout_initialization_enabled, true);
    for (0..3) |branch| {
        branches.stage_day[branch * 10] = 50;
        branches.stage_day[branch * 10 + 3] = 80;
    }
    const reproductive_stage = [_]f64{ 1.25, 2.5, 3.75 };
    const reset = try resetDevelopmentAfterCut(&branches, .{ .selected_branch = 0, .main_branch = 0, .plant_branch_first = 0, .plant_branch_end = 3, .biomass_turnover_type = 0, .root_profile_type = 1, .grazing = false, .canopy_height_m = 2, .cutting_height_m = 0.5, .winter_phenology_type = 0, .accumulated_vernalization_h = 0, .required_vernalization_h = 0, .vernalization_reset_fraction = 0.5, .initial_maturity_group = 2, .current_reproductive_stage_by_branch = &reproductive_stage, .current_day = 120 });
    try std.testing.expect(reset);
    for (0..3) |branch| {
        try std.testing.expectEqual(2, branches.maturity_group[branch]);
        try std.testing.expectEqual(reproductive_stage[branch], branches.initial_reproductive_stage[branch]);
        try std.testing.expectEqual(0, branches.final_reproductive_stage[branch]);
        try std.testing.expectEqual(@as(u32, 120), branches.stage_day[branch * 10]);
        try std.testing.expectEqual(@as(u32, 0), branches.stage_day[branch * 10 + 3]);
        try std.testing.expect(!branches.leafout_initialization_enabled[branch]);
        try std.testing.expectEqual(@as(f64, 48), branches.remobilization_progress_h[branch]);
    }
    try terminatePlantBranches(&branches, 0, 3, true);
    for (branches.dead) |dead| try std.testing.expect(dead);
    try branches.validateFinite();
}

test "GROSUB grazing does not reset post-cut development" {
    var branches = try BranchDevelopmentState.init(std.testing.allocator, 1);
    defer branches.deinit();
    branches.stage_day[0] = 10;
    branches.maturity_group[0] = 7;
    const reproductive_stage = [_]f64{1};
    const reset = try resetDevelopmentAfterCut(&branches, .{ .selected_branch = 0, .main_branch = 0, .plant_branch_first = 0, .plant_branch_end = 1, .biomass_turnover_type = 0, .root_profile_type = 1, .grazing = true, .canopy_height_m = 2, .cutting_height_m = 0.5, .winter_phenology_type = 0, .accumulated_vernalization_h = 0, .required_vernalization_h = 0, .vernalization_reset_fraction = 0.5, .initial_maturity_group = 2, .current_reproductive_stage_by_branch = &reproductive_stage, .current_day = 20 });
    try std.testing.expect(!reset);
    try std.testing.expectEqual(7, branches.maturity_group[0]);
}

test "annual pre-floral rates respond to water and oxygen stress" {
    const allocator = std.testing.allocator;
    const config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var plants = try PlantState.init(allocator, config);
    defer plants.deinit();
    plants.canopy_temperature_k[0] = 298.15;
    var state = try State.init(allocator, 1, 1);
    defer state.deinit();
    state.active[0] = true;
    state.annual_growth_habit[0] = true;
    state.node_initiation_rate_at_25c_per_h[0] = 0.02;
    state.leaf_appearance_rate_at_25c_per_h[0] = 0.03;
    const parameters = compatibilityParameters();
    var context: AdvanceContext = .{ .phenology = &state, .plants = &plants, .thermal_acclimation_offset_k = &.{0}, .canopy_turgor_potential_mpa = &.{0.6}, .root_oxygen_uptake_to_demand_fraction = &.{0.0625}, .emerged = &.{true}, .timestep_hours = 1, .parameters = parameters };
    try advanceTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expectApproxEqRel(@as(f64, 0.005), state.node_initiation_per_timestep[0], 0.02);
    try std.testing.expectApproxEqRel(@as(f64, 0.0075), state.leaf_appearance_per_timestep[0], 0.02);
}

test "HFUNC pool concentrations retain separate denominators" {
    const value = try nonstructuralConcentrations(10, 2, 1, 0.5, 0.25, true, 1e-12);
    try std.testing.expectApproxEqAbs(2.0 / 12.0, value.carbon_g_per_g, 1e-15);
    try std.testing.expectApproxEqAbs(1.0 / 11.0, value.nitrogen_g_per_g, 1e-15);
    try std.testing.expectApproxEqAbs(0.5 / 10.5, value.phosphorus_g_per_g, 1e-15);
    try std.testing.expectApproxEqAbs(0.25 / 10.25, value.salt_mol_per_g_c, 1e-15);
}

test "HFUNC shoot and root branching use runtime limits" {
    const shoot = try shootBranchDecision(false, true, 100, 0.5, 0.1, 1, 0, 0, 0, 1, 0.2, 0.1, 0, 0, 1, 0, 24, 2, 1, 10);
    try std.testing.expect(shoot.initiate);
    try std.testing.expectEqual(@as(usize, 0), shoot.branch_order);
    try std.testing.expect(try shouldAddRootAxis(0.5, 0.1, 12, 20, 2, 1, 1, 0.2, 0.1, 64));
    try std.testing.expect(!(try shouldAddRootAxis(0.5, 0.1, 64, 100, 2, 1, 1, 0.2, 0.1, 64)));
}

test "phenology replant dates preserve DAY modulo-four chronology" {
    var state = try State.init(std.testing.allocator, 1, 1);
    defer state.deinit();

    state.replant_day_of_year[0] = 366;
    state.replant_year[0] = 1900;
    try state.validateFinite();

    state.replant_year[0] = 1901;
    try std.testing.expectError(
        error.InvalidPlantPhenologyReplantDate,
        state.validateFinite(),
    );

    state.replant_day_of_year[0] = 0;
    try std.testing.expectError(
        error.InvalidPlantPhenologyReplantDate,
        state.validateFinite(),
    );

    state.replant_year[0] = 0;
    try state.validateFinite();
}
