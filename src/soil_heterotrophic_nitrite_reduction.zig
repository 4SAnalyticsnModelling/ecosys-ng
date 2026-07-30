const std = @import("std");

pub const Inputs = struct {
    enabled: []const bool,
    previous_total_non_band_nitrite_demand_g_n: []const f64,
    previous_total_band_nitrite_demand_g_n: []const f64,
    previous_non_band_nitrite_capacity_g_n: []const f64,
    previous_band_nitrite_capacity_g_n: []const f64,
    microbial_active_fraction: []const f64,
    non_band_nitrate_fraction: []const f64,
    band_nitrate_fraction: []const f64,
    non_band_nitrite_fraction: []const f64,
    band_nitrite_fraction: []const f64,
    preceding_nitrate_demand_g_n: []const f64,
    preceding_total_nitrate_reduction_g_n: []const f64,
    non_band_nitrate_reduction_g_n: []const f64,
    band_nitrate_reduction_g_n: []const f64,
    non_band_nitrite_concentration_g_n_per_m3: []const f64,
    band_nitrite_concentration_g_n_per_m3: []const f64,
    nitrous_oxide_concentration_g_n_per_m3: []const f64,
    biologically_active_water_m3: []const f64,
    substrate_complex_fraction: []const f64,
    preceding_remaining_doc_g_c: []const f64,
    preceding_nitrate_respiration_g_c: []const f64,
    non_band_nitrite_g_n: []const f64,
    band_nitrite_g_n: []const f64,
    nitrite_half_saturation_g_n_per_m3: []const f64,
    nitrous_oxide_half_saturation_g_n_per_m3: []const f64,
    product_inhibition_g_n_per_m3_step: []const f64,
    carbon_per_nitrite_n_g_c_per_g_n: []const f64,
    minimum_competition_fraction: f64,
    negligible_demand_g_n: f64,
    negligible_water_m3: f64,
    timestep_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    total_non_band_competition_fraction: f64,
    total_band_competition_fraction: f64,
    non_band_competition_fraction: []f64,
    band_competition_fraction: []f64,
    nitrite_reduction_demand_g_n: []f64,
    raw_non_band_capacity_g_n: []f64,
    raw_band_capacity_g_n: []f64,
    raw_total_capacity_g_n: []f64,
    product_inhibition_fraction: []f64,
    non_band_capacity_g_n: []f64,
    band_capacity_g_n: []f64,
    remaining_doc_g_c: []f64,
    doc_nitrite_budget_g_n: []f64,
    non_band_doc_budget_g_n: []f64,
    band_doc_budget_g_n: []f64,
    non_band_nitrite_supply_g_n: []f64,
    band_nitrite_supply_g_n: []f64,
    non_band_supply_capacity_reduction_g_n: []f64,
    band_supply_capacity_reduction_g_n: []f64,
    non_band_nitrite_reduction_g_n: []f64,
    band_nitrite_reduction_g_n: []f64,
    total_supply_capacity_reduction_g_n: []f64,
    total_nitrite_reduction_g_n: []f64,
    supply_capacity_respiration_g_c: []f64,
    nitrite_reduction_respiration_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidNitriteReductionDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
        state.total_non_band_competition_fraction = 0;
        state.total_band_competition_fraction = 0;
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

/// NITRO.F 1771--1868. Inputs explicitly carry the preceding nitrate stage.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.unit_count);
    defer staged.deinit();
    for (0..state.unit_count) |unit| {
        if (!inputs.enabled[unit]) continue;
        const non_band_competition = @max(
            inputs.minimum_competition_fraction,
            if (inputs.previous_total_non_band_nitrite_demand_g_n[unit] >
                inputs.negligible_demand_g_n)
                inputs.previous_non_band_nitrite_capacity_g_n[unit] /
                    inputs.previous_total_non_band_nitrite_demand_g_n[unit]
            else
                inputs.microbial_active_fraction[unit] *
                    inputs.non_band_nitrate_fraction[unit],
        );
        const band_competition = @max(
            inputs.minimum_competition_fraction,
            if (inputs.previous_total_band_nitrite_demand_g_n[unit] >
                inputs.negligible_demand_g_n)
                inputs.previous_band_nitrite_capacity_g_n[unit] /
                    inputs.previous_total_band_nitrite_demand_g_n[unit]
            else
                inputs.microbial_active_fraction[unit] *
                    inputs.band_nitrate_fraction[unit],
        );
        staged.total_non_band_competition_fraction += non_band_competition;
        staged.total_band_competition_fraction += band_competition;
        calculateReduction(&staged, inputs, unit, non_band_competition, band_competition);
    }
    try validateResult(&staged);
    state.total_non_band_competition_fraction =
        staged.total_non_band_competition_fraction;
    state.total_band_competition_fraction =
        staged.total_band_competition_fraction;
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64)
            @memcpy(@field(state, field.name), @field(staged, field.name));
}

fn calculateReduction(
    state: *State,
    inputs: Inputs,
    unit: usize,
    non_band_competition: f64,
    band_competition: f64,
) void {
    const demand = inputs.preceding_nitrate_demand_g_n[unit] -
        inputs.preceding_total_nitrate_reduction_g_n[unit];
    const raw_non_band = rawCapacity(inputs, unit, false, demand);
    const raw_band = rawCapacity(inputs, unit, true, demand);
    const raw_total = raw_non_band + raw_band;
    const inhibition = if (inputs.biologically_active_water_m3[unit] >
        inputs.negligible_water_m3 and
        inputs.substrate_complex_fraction[unit] > 0)
        1 / (1 + raw_total /
            (inputs.product_inhibition_g_n_per_m3_step[unit] *
                inputs.biologically_active_water_m3[unit] *
                inputs.substrate_complex_fraction[unit]))
    else
        0;
    const non_band_capacity = raw_non_band * inhibition;
    const band_capacity = raw_band * inhibition;
    const remaining_doc = @max(
        0,
        inputs.preceding_remaining_doc_g_c[unit] -
            inputs.preceding_nitrate_respiration_g_c[unit],
    );
    const doc_budget =
        remaining_doc / inputs.carbon_per_nitrite_n_g_c_per_g_n[unit] *
        inputs.timestep_h;
    const non_band_doc =
        doc_budget * inputs.non_band_nitrate_fraction[unit];
    const band_doc = doc_budget * inputs.band_nitrate_fraction[unit];
    const non_band_supply =
        (inputs.non_band_nitrite_g_n[unit] * inputs.timestep_h +
            inputs.non_band_nitrate_reduction_g_n[unit]) *
        non_band_competition;
    const band_supply =
        (inputs.band_nitrite_g_n[unit] * inputs.timestep_h +
            inputs.band_nitrate_reduction_g_n[unit]) *
        band_competition;
    const non_band_supply_capacity =
        @max(0, @min(non_band_supply, non_band_capacity));
    const band_supply_capacity =
        @max(0, @min(band_supply, band_capacity));
    const non_band_reduction =
        @max(0, @min(non_band_capacity, non_band_doc, non_band_supply));
    const band_reduction =
        @max(0, @min(band_capacity, band_doc, band_supply));
    publish(state, inputs, unit, .{
        non_band_competition,
        band_competition,
        demand,
        raw_non_band,
        raw_band,
        raw_total,
        inhibition,
        non_band_capacity,
        band_capacity,
        remaining_doc,
        doc_budget,
        non_band_doc,
        band_doc,
        non_band_supply,
        band_supply,
        non_band_supply_capacity,
        band_supply_capacity,
        non_band_reduction,
        band_reduction,
    });
}

fn rawCapacity(inputs: Inputs, unit: usize, band: bool, demand: f64) f64 {
    const nitrite = if (band)
        inputs.band_nitrite_concentration_g_n_per_m3[unit]
    else
        inputs.non_band_nitrite_concentration_g_n_per_m3[unit];
    if (nitrite <= 0) return 0;
    const fraction = if (band)
        inputs.band_nitrite_fraction[unit]
    else
        inputs.non_band_nitrite_fraction[unit];
    return fraction * demand * nitrite /
        (nitrite + inputs.nitrite_half_saturation_g_n_per_m3[unit]) /
        (1 + inputs.nitrous_oxide_concentration_g_n_per_m3[unit] *
            inputs.nitrite_half_saturation_g_n_per_m3[unit] /
            (nitrite *
                inputs.nitrous_oxide_half_saturation_g_n_per_m3[unit]));
}

fn publish(state: *State, inputs: Inputs, unit: usize, value: [19]f64) void {
    state.non_band_competition_fraction[unit] = value[0];
    state.band_competition_fraction[unit] = value[1];
    state.nitrite_reduction_demand_g_n[unit] = value[2];
    state.raw_non_band_capacity_g_n[unit] = value[3];
    state.raw_band_capacity_g_n[unit] = value[4];
    state.raw_total_capacity_g_n[unit] = value[5];
    state.product_inhibition_fraction[unit] = value[6];
    state.non_band_capacity_g_n[unit] = value[7];
    state.band_capacity_g_n[unit] = value[8];
    state.remaining_doc_g_c[unit] = value[9];
    state.doc_nitrite_budget_g_n[unit] = value[10];
    state.non_band_doc_budget_g_n[unit] = value[11];
    state.band_doc_budget_g_n[unit] = value[12];
    state.non_band_nitrite_supply_g_n[unit] = value[13];
    state.band_nitrite_supply_g_n[unit] = value[14];
    state.non_band_supply_capacity_reduction_g_n[unit] = value[15];
    state.band_supply_capacity_reduction_g_n[unit] = value[16];
    state.non_band_nitrite_reduction_g_n[unit] = value[17];
    state.band_nitrite_reduction_g_n[unit] = value[18];
    state.total_supply_capacity_reduction_g_n[unit] = value[15] + value[16];
    state.total_nitrite_reduction_g_n[unit] = value[17] + value[18];
    state.supply_capacity_respiration_g_c[unit] =
        inputs.carbon_per_nitrite_n_g_c_per_g_n[unit] *
        state.total_supply_capacity_reduction_g_n[unit];
    state.nitrite_reduction_respiration_g_c[unit] =
        inputs.carbon_per_nitrite_n_g_c_per_g_n[unit] *
        state.total_nitrite_reduction_g_n[unit];
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.enabled.len != n)
        return error.InvalidNitriteReductionDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == []const f64) {
            const values = @field(inputs, field.name);
            if (values.len != n) return error.InvalidNitriteReductionDimensions;
            for (values) |value| if (!std.math.isFinite(value) or value < 0)
                return error.InvalidNitriteReductionInput;
        } else if (field.type == f64) {
            const value = @field(inputs, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidNitriteReductionInput;
        }
    }
    for (0..n) |unit| if (inputs.enabled[unit]) {
        if (inputs.preceding_total_nitrate_reduction_g_n[unit] >
            inputs.preceding_nitrate_demand_g_n[unit] or
            inputs.nitrite_half_saturation_g_n_per_m3[unit] <= 0 or
            inputs.nitrous_oxide_half_saturation_g_n_per_m3[unit] <= 0 or
            inputs.product_inhibition_g_n_per_m3_step[unit] <= 0 or
            inputs.carbon_per_nitrite_n_g_c_per_g_n[unit] <= 0)
            return error.InvalidNitriteReductionInput;
    };
    if (inputs.timestep_h <= 0) return error.InvalidNitriteReductionInput;
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (field.type == f64) {
            const value = @field(state.*, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteNitriteReductionResult;
        } else if (field.type == []f64) {
            for (@field(state.*, field.name)) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteNitriteReductionResult;
        }
    }
}

test "NITRO nitrite stage consumes nitrate products in source order" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    const one = [_]f64{1};
    const zero = [_]f64{0};
    try calculate(&state, .{
        .enabled = &.{true},
        .previous_total_non_band_nitrite_demand_g_n = &zero,
        .previous_total_band_nitrite_demand_g_n = &zero,
        .previous_non_band_nitrite_capacity_g_n = &zero,
        .previous_band_nitrite_capacity_g_n = &zero,
        .microbial_active_fraction = &one,
        .non_band_nitrate_fraction = &.{0.75},
        .band_nitrate_fraction = &.{0.25},
        .non_band_nitrite_fraction = &.{0.75},
        .band_nitrite_fraction = &.{0.25},
        .preceding_nitrate_demand_g_n = &.{4},
        .preceding_total_nitrate_reduction_g_n = &.{1},
        .non_band_nitrate_reduction_g_n = &.{0.75},
        .band_nitrate_reduction_g_n = &.{0.25},
        .non_band_nitrite_concentration_g_n_per_m3 = &one,
        .band_nitrite_concentration_g_n_per_m3 = &one,
        .nitrous_oxide_concentration_g_n_per_m3 = &zero,
        .biologically_active_water_m3 = &one,
        .substrate_complex_fraction = &one,
        .preceding_remaining_doc_g_c = &.{10},
        .preceding_nitrate_respiration_g_c = &one,
        .non_band_nitrite_g_n = &one,
        .band_nitrite_g_n = &one,
        .nitrite_half_saturation_g_n_per_m3 = &one,
        .nitrous_oxide_half_saturation_g_n_per_m3 = &one,
        .product_inhibition_g_n_per_m3_step = &one,
        .carbon_per_nitrite_n_g_c_per_g_n = &one,
        .minimum_competition_fraction = 0.01,
        .negligible_demand_g_n = 1e-12,
        .negligible_water_m3 = 1e-12,
        .timestep_h = 1,
    });
    try std.testing.expectEqual(@as(f64, 3), state.nitrite_reduction_demand_g_n[0]);
    try std.testing.expect(state.non_band_nitrite_supply_g_n[0] > 0.75);
    try std.testing.expect(state.total_nitrite_reduction_g_n[0] > 0);
}
