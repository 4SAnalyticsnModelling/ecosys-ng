const std = @import("std");

pub const Inputs = struct {
    enabled: []const bool,
    previous_total_nitrous_oxide_demand_g_n: []const f64,
    previous_nitrous_oxide_capacity_g_n: []const f64,
    microbial_active_fraction: []const f64,
    preceding_nitrite_demand_g_n: []const f64,
    preceding_total_nitrite_reduction_g_n: []const f64,
    nitrous_oxide_concentration_g_n_per_m3: []const f64,
    biologically_active_water_m3: []const f64,
    substrate_complex_fraction: []const f64,
    preceding_remaining_doc_g_c: []const f64,
    preceding_nitrite_respiration_g_c: []const f64,
    nitrous_oxide_g_n: []const f64,
    nitrous_oxide_half_saturation_g_n_per_m3: []const f64,
    product_inhibition_g_n_per_m3_step: []const f64,
    carbon_per_nitrous_oxide_n_g_c_per_g_n: []const f64,
    nitrate_supply_capacity_respiration_g_c: []const f64,
    nitrite_supply_capacity_respiration_g_c: []const f64,
    nitrate_limited_respiration_g_c: []const f64,
    nitrite_limited_respiration_g_c: []const f64,
    minimum_competition_fraction: f64,
    negligible_demand_g_n: f64,
    negligible_water_m3: f64,
    nitrous_oxide_demand_multiplier: f64,
    timestep_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    total_competition_fraction: f64,
    competition_fraction: []f64,
    nitrous_oxide_reduction_demand_g_n: []f64,
    raw_capacity_g_n: []f64,
    product_inhibition_fraction: []f64,
    capacity_g_n: []f64,
    remaining_doc_g_c: []f64,
    doc_budget_g_n: []f64,
    nitrous_oxide_supply_g_n: []f64,
    supply_capacity_reduction_g_n: []f64,
    nitrous_oxide_reduction_g_n: []f64,
    supply_capacity_respiration_g_c: []f64,
    nitrous_oxide_reduction_respiration_g_c: []f64,
    total_supply_capacity_denitrification_respiration_g_c: []f64,
    total_limited_denitrification_respiration_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0)
            return error.InvalidNitrousOxideReductionDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
        state.total_competition_fraction = 0;
        var allocated: usize = 0;
        errdefer {
            inline for (@typeInfo(State).@"struct".fields) |field| {
                if (field.type == []f64 and allocated > 0) {
                    allocated -= 1;
                    allocator.free(@field(state, field.name));
                }
            }
        }
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64) {
                @field(state, field.name) = try allocator.alloc(f64, unit_count);
                @memset(@field(state, field.name), 0);
                allocated += 1;
            }
        }
        return state;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

/// NITRO.F 1870--1936. Inputs explicitly carry nitrate/nitrite stage outputs.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.unit_count);
    defer staged.deinit();
    for (0..state.unit_count) |unit| {
        if (!inputs.enabled[unit]) continue;
        const competition = @max(
            inputs.minimum_competition_fraction,
            if (inputs.previous_total_nitrous_oxide_demand_g_n[unit] >
                inputs.negligible_demand_g_n)
                inputs.previous_nitrous_oxide_capacity_g_n[unit] /
                    inputs.previous_total_nitrous_oxide_demand_g_n[unit]
            else
                inputs.microbial_active_fraction[unit],
        );
        staged.total_competition_fraction += competition;
        calculateReduction(&staged, inputs, unit, competition);
    }
    try validateResult(&staged);
    state.total_competition_fraction = staged.total_competition_fraction;
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64)
            @memcpy(@field(state, field.name), @field(staged, field.name));
}

fn calculateReduction(
    state: *State,
    inputs: Inputs,
    unit: usize,
    competition: f64,
) void {
    const demand =
        (inputs.preceding_nitrite_demand_g_n[unit] -
            inputs.preceding_total_nitrite_reduction_g_n[unit]) *
        inputs.nitrous_oxide_demand_multiplier;
    const concentration =
        inputs.nitrous_oxide_concentration_g_n_per_m3[unit];
    const raw_capacity =
        demand * concentration /
        (concentration +
            inputs.nitrous_oxide_half_saturation_g_n_per_m3[unit]);
    const inhibition = if (inputs.biologically_active_water_m3[unit] >
        inputs.negligible_water_m3 and
        inputs.substrate_complex_fraction[unit] > 0)
        1 / (1 + raw_capacity /
            (inputs.product_inhibition_g_n_per_m3_step[unit] *
                inputs.biologically_active_water_m3[unit] *
                inputs.substrate_complex_fraction[unit]))
    else
        0;
    const capacity = raw_capacity * inhibition;
    const remaining_doc = @max(
        0,
        inputs.preceding_remaining_doc_g_c[unit] -
            inputs.preceding_nitrite_respiration_g_c[unit],
    );
    const doc_budget =
        remaining_doc /
        inputs.carbon_per_nitrous_oxide_n_g_c_per_g_n[unit] *
        inputs.timestep_h;
    const supply =
        (inputs.nitrous_oxide_g_n[unit] * inputs.timestep_h +
            inputs.preceding_total_nitrite_reduction_g_n[unit]) *
        competition;
    const supply_capacity = @max(0, @min(supply, capacity));
    const reduction = @max(0, @min(capacity, doc_budget, supply));
    const supply_respiration =
        inputs.carbon_per_nitrous_oxide_n_g_c_per_g_n[unit] *
        supply_capacity;
    const limited_respiration =
        inputs.carbon_per_nitrous_oxide_n_g_c_per_g_n[unit] * reduction;
    state.competition_fraction[unit] = competition;
    state.nitrous_oxide_reduction_demand_g_n[unit] = demand;
    state.raw_capacity_g_n[unit] = raw_capacity;
    state.product_inhibition_fraction[unit] = inhibition;
    state.capacity_g_n[unit] = capacity;
    state.remaining_doc_g_c[unit] = remaining_doc;
    state.doc_budget_g_n[unit] = doc_budget;
    state.nitrous_oxide_supply_g_n[unit] = supply;
    state.supply_capacity_reduction_g_n[unit] = supply_capacity;
    state.nitrous_oxide_reduction_g_n[unit] = reduction;
    state.supply_capacity_respiration_g_c[unit] = supply_respiration;
    state.nitrous_oxide_reduction_respiration_g_c[unit] =
        limited_respiration;
    state.total_supply_capacity_denitrification_respiration_g_c[unit] =
        inputs.nitrate_supply_capacity_respiration_g_c[unit] +
        inputs.nitrite_supply_capacity_respiration_g_c[unit] +
        supply_respiration;
    state.total_limited_denitrification_respiration_g_c[unit] =
        inputs.nitrate_limited_respiration_g_c[unit] +
        inputs.nitrite_limited_respiration_g_c[unit] +
        limited_respiration;
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.enabled.len != n)
        return error.InvalidNitrousOxideReductionDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == []const f64) {
            const values = @field(inputs, field.name);
            if (values.len != n)
                return error.InvalidNitrousOxideReductionDimensions;
            for (values) |value| if (!std.math.isFinite(value) or value < 0)
                return error.InvalidNitrousOxideReductionInput;
        } else if (field.type == f64) {
            const value = @field(inputs, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidNitrousOxideReductionInput;
        }
    }
    for (0..n) |unit| if (inputs.enabled[unit]) {
        if (inputs.preceding_total_nitrite_reduction_g_n[unit] >
            inputs.preceding_nitrite_demand_g_n[unit] or
            inputs.nitrous_oxide_half_saturation_g_n_per_m3[unit] <= 0 or
            inputs.product_inhibition_g_n_per_m3_step[unit] <= 0 or
            inputs.carbon_per_nitrous_oxide_n_g_c_per_g_n[unit] <= 0)
            return error.InvalidNitrousOxideReductionInput;
    };
    if (inputs.nitrous_oxide_demand_multiplier <= 0 or inputs.timestep_h <= 0)
        return error.InvalidNitrousOxideReductionInput;
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (field.type == f64) {
            const value = @field(state.*, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteNitrousOxideReductionResult;
        } else if (field.type == []f64) {
            for (@field(state.*, field.name)) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteNitrousOxideReductionResult;
        }
    }
}

test "NITRO terminal denitrification stage preserves sequential totals" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    const one = [_]f64{1};
    const zero = [_]f64{0};
    try calculate(&state, .{
        .enabled = &.{true},
        .previous_total_nitrous_oxide_demand_g_n = &zero,
        .previous_nitrous_oxide_capacity_g_n = &zero,
        .microbial_active_fraction = &.{0.5},
        .preceding_nitrite_demand_g_n = &.{4},
        .preceding_total_nitrite_reduction_g_n = &one,
        .nitrous_oxide_concentration_g_n_per_m3 = &one,
        .biologically_active_water_m3 = &one,
        .substrate_complex_fraction = &one,
        .preceding_remaining_doc_g_c = &.{10},
        .preceding_nitrite_respiration_g_c = &one,
        .nitrous_oxide_g_n = &one,
        .nitrous_oxide_half_saturation_g_n_per_m3 = &one,
        .product_inhibition_g_n_per_m3_step = &one,
        .carbon_per_nitrous_oxide_n_g_c_per_g_n = &one,
        .nitrate_supply_capacity_respiration_g_c = &.{2},
        .nitrite_supply_capacity_respiration_g_c = &.{3},
        .nitrate_limited_respiration_g_c = &one,
        .nitrite_limited_respiration_g_c = &.{2},
        .minimum_competition_fraction = 0.01,
        .negligible_demand_g_n = 1e-12,
        .negligible_water_m3 = 1e-12,
        .nitrous_oxide_demand_multiplier = 2,
        .timestep_h = 1,
    });
    try std.testing.expectEqual(@as(f64, 6), state.nitrous_oxide_reduction_demand_g_n[0]);
    try std.testing.expectEqual(@as(f64, 0.5), state.competition_fraction[0]);
    try std.testing.expectEqual(
        @as(f64, 5) + state.supply_capacity_respiration_g_c[0],
        state.total_supply_capacity_denitrification_respiration_g_c[0],
    );
}
