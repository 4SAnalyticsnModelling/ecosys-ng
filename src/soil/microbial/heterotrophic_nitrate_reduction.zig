const std = @import("std");

pub const Inputs = struct {
    enabled: []const bool,
    previous_total_non_band_nitrate_demand_g_n: []const f64,
    previous_total_band_nitrate_demand_g_n: []const f64,
    previous_non_band_nitrate_capacity_g_n: []const f64,
    previous_band_nitrate_capacity_g_n: []const f64,
    microbial_active_fraction: []const f64,
    non_band_nitrate_fraction: []const f64,
    band_nitrate_fraction: []const f64,
    oxygen_demand_g_o: []const f64,
    oxygen_uptake_g_o: []const f64,
    non_band_nitrate_concentration_g_n_per_m3: []const f64,
    band_nitrate_concentration_g_n_per_m3: []const f64,
    non_band_nitrite_concentration_g_n_per_m3: []const f64,
    band_nitrite_concentration_g_n_per_m3: []const f64,
    biologically_active_water_m3: []const f64,
    substrate_complex_fraction: []const f64,
    doc_supply_g_c: []const f64,
    doc_competition_fraction: []const f64,
    preceding_aerobic_respiration_g_c: []const f64,
    oxygen_satisfaction_fraction: []const f64,
    non_band_nitrate_g_n: []const f64,
    band_nitrate_g_n: []const f64,
    nitrate_half_saturation_g_n_per_m3: []const f64,
    nitrite_half_saturation_g_n_per_m3: []const f64,
    product_inhibition_g_n_per_m3_step: []const f64,
    carbon_per_nitrate_n_g_c_per_g_n: []const f64,
    nitrate_n_per_unmet_oxygen_g_n_per_g_o: f64,
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
    unmet_oxygen_demand_g_o: []f64,
    nitrate_reduction_demand_g_n: []f64,
    raw_non_band_capacity_g_n: []f64,
    raw_band_capacity_g_n: []f64,
    raw_total_capacity_g_n: []f64,
    product_inhibition_fraction: []f64,
    non_band_capacity_g_n: []f64,
    band_capacity_g_n: []f64,
    remaining_doc_g_c: []f64,
    doc_nitrate_budget_g_n: []f64,
    non_band_doc_budget_g_n: []f64,
    band_doc_budget_g_n: []f64,
    non_band_nitrate_supply_g_n: []f64,
    band_nitrate_supply_g_n: []f64,
    non_band_supply_capacity_reduction_g_n: []f64,
    band_supply_capacity_reduction_g_n: []f64,
    non_band_nitrate_reduction_g_n: []f64,
    band_nitrate_reduction_g_n: []f64,
    total_supply_capacity_reduction_g_n: []f64,
    total_nitrate_reduction_g_n: []f64,
    supply_capacity_respiration_g_c: []f64,
    nitrate_reduction_respiration_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidNitrateReductionDimensions;
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

/// NITRO.F 1667--1769. Historical K<=4,N=2 admission is runtime data.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.unit_count);
    defer staged.deinit();
    for (0..state.unit_count) |unit| {
        if (!inputs.enabled[unit]) continue;
        const non_band_competition = @max(
            inputs.minimum_competition_fraction,
            if (inputs.previous_total_non_band_nitrate_demand_g_n[unit] >
                inputs.negligible_demand_g_n)
                inputs.previous_non_band_nitrate_capacity_g_n[unit] /
                    inputs.previous_total_non_band_nitrate_demand_g_n[unit]
            else
                inputs.microbial_active_fraction[unit] *
                    inputs.non_band_nitrate_fraction[unit],
        );
        const band_competition = @max(
            inputs.minimum_competition_fraction,
            if (inputs.previous_total_band_nitrate_demand_g_n[unit] >
                inputs.negligible_demand_g_n)
                inputs.previous_band_nitrate_capacity_g_n[unit] /
                    inputs.previous_total_band_nitrate_demand_g_n[unit]
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
    const unmet_oxygen = @max(
        0,
        inputs.oxygen_demand_g_o[unit] - inputs.oxygen_uptake_g_o[unit],
    );
    const demand =
        inputs.nitrate_n_per_unmet_oxygen_g_n_per_g_o * unmet_oxygen;
    const raw_non_band = rawCapacity(inputs, unit, false, demand);
    const raw_band = rawCapacity(inputs, unit, true, demand);
    const raw_total = raw_non_band + raw_band;
    const inhibition = if (inputs.biologically_active_water_m3[unit] >
        inputs.negligible_water_m3 and
        inputs.substrate_complex_fraction[unit] > 0)
        1 / (1 +
            raw_total /
                (inputs.product_inhibition_g_n_per_m3_step[unit] *
                    inputs.biologically_active_water_m3[unit] *
                    inputs.substrate_complex_fraction[unit]))
    else
        0;
    const non_band_capacity = raw_non_band * inhibition;
    const band_capacity = raw_band * inhibition;
    const remaining_doc = @max(
        0,
        inputs.doc_supply_g_c[unit] * inputs.doc_competition_fraction[unit] -
            inputs.preceding_aerobic_respiration_g_c[unit] *
                inputs.oxygen_satisfaction_fraction[unit],
    );
    const doc_budget =
        remaining_doc / inputs.carbon_per_nitrate_n_g_c_per_g_n[unit] *
        inputs.timestep_h;
    const non_band_doc =
        doc_budget * inputs.non_band_nitrate_fraction[unit];
    const band_doc = doc_budget * inputs.band_nitrate_fraction[unit];
    const non_band_supply =
        inputs.non_band_nitrate_g_n[unit] * non_band_competition *
        inputs.timestep_h;
    const band_supply =
        inputs.band_nitrate_g_n[unit] * band_competition * inputs.timestep_h;
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
        unmet_oxygen,
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
    const nitrate = if (band)
        inputs.band_nitrate_concentration_g_n_per_m3[unit]
    else
        inputs.non_band_nitrate_concentration_g_n_per_m3[unit];
    if (nitrate <= 0) return 0;
    const nitrite = if (band)
        inputs.band_nitrite_concentration_g_n_per_m3[unit]
    else
        inputs.non_band_nitrite_concentration_g_n_per_m3[unit];
    const fraction = if (band)
        inputs.band_nitrate_fraction[unit]
    else
        inputs.non_band_nitrate_fraction[unit];
    return fraction * demand * nitrate /
        (nitrate + inputs.nitrate_half_saturation_g_n_per_m3[unit]) /
        (1 + nitrite * inputs.nitrate_half_saturation_g_n_per_m3[unit] /
            (nitrate * inputs.nitrite_half_saturation_g_n_per_m3[unit]));
}

fn publish(state: *State, inputs: Inputs, unit: usize, value: [20]f64) void {
    state.non_band_competition_fraction[unit] = value[0];
    state.band_competition_fraction[unit] = value[1];
    state.unmet_oxygen_demand_g_o[unit] = value[2];
    state.nitrate_reduction_demand_g_n[unit] = value[3];
    state.raw_non_band_capacity_g_n[unit] = value[4];
    state.raw_band_capacity_g_n[unit] = value[5];
    state.raw_total_capacity_g_n[unit] = value[6];
    state.product_inhibition_fraction[unit] = value[7];
    state.non_band_capacity_g_n[unit] = value[8];
    state.band_capacity_g_n[unit] = value[9];
    state.remaining_doc_g_c[unit] = value[10];
    state.doc_nitrate_budget_g_n[unit] = value[11];
    state.non_band_doc_budget_g_n[unit] = value[12];
    state.band_doc_budget_g_n[unit] = value[13];
    state.non_band_nitrate_supply_g_n[unit] = value[14];
    state.band_nitrate_supply_g_n[unit] = value[15];
    state.non_band_supply_capacity_reduction_g_n[unit] = value[16];
    state.band_supply_capacity_reduction_g_n[unit] = value[17];
    state.non_band_nitrate_reduction_g_n[unit] = value[18];
    state.band_nitrate_reduction_g_n[unit] = value[19];
    state.total_supply_capacity_reduction_g_n[unit] = value[16] + value[17];
    state.total_nitrate_reduction_g_n[unit] = value[18] + value[19];
    state.supply_capacity_respiration_g_c[unit] =
        inputs.carbon_per_nitrate_n_g_c_per_g_n[unit] *
        state.total_supply_capacity_reduction_g_n[unit];
    state.nitrate_reduction_respiration_g_c[unit] =
        inputs.carbon_per_nitrate_n_g_c_per_g_n[unit] *
        state.total_nitrate_reduction_g_n[unit];
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.enabled.len != n)
        return error.InvalidNitrateReductionDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == []const f64) {
            const values = @field(inputs, field.name);
            if (values.len != n) return error.InvalidNitrateReductionDimensions;
            for (values) |value| if (!std.math.isFinite(value) or value < 0)
                return error.InvalidNitrateReductionInput;
        } else if (field.type == f64) {
            const value = @field(inputs, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidNitrateReductionInput;
        }
    }
    for (0..n) |unit| if (inputs.enabled[unit]) {
        if (inputs.nitrate_half_saturation_g_n_per_m3[unit] <= 0 or
            inputs.nitrite_half_saturation_g_n_per_m3[unit] <= 0 or
            inputs.product_inhibition_g_n_per_m3_step[unit] <= 0 or
            inputs.carbon_per_nitrate_n_g_c_per_g_n[unit] <= 0)
            return error.InvalidNitrateReductionInput;
    };
    if (inputs.timestep_h <= 0)
        return error.InvalidNitrateReductionInput;
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (field.type == f64) {
            const value = @field(state.*, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteNitrateReductionResult;
        } else if (field.type == []f64) {
            for (@field(state.*, field.name)) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteNitrateReductionResult;
        }
    }
}

test "NITRO nitrate stage preserves competition and sequential limits" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    const one = [_]f64{1};
    const zero = [_]f64{0};
    try calculate(&state, .{
        .enabled = &.{true},
        .previous_total_non_band_nitrate_demand_g_n = &one,
        .previous_total_band_nitrate_demand_g_n = &zero,
        .previous_non_band_nitrate_capacity_g_n = &.{0.25},
        .previous_band_nitrate_capacity_g_n = &zero,
        .microbial_active_fraction = &one,
        .non_band_nitrate_fraction = &.{0.75},
        .band_nitrate_fraction = &.{0.25},
        .oxygen_demand_g_o = &.{4},
        .oxygen_uptake_g_o = &zero,
        .non_band_nitrate_concentration_g_n_per_m3 = &one,
        .band_nitrate_concentration_g_n_per_m3 = &one,
        .non_band_nitrite_concentration_g_n_per_m3 = &zero,
        .band_nitrite_concentration_g_n_per_m3 = &zero,
        .biologically_active_water_m3 = &one,
        .substrate_complex_fraction = &one,
        .doc_supply_g_c = &.{10},
        .doc_competition_fraction = &one,
        .preceding_aerobic_respiration_g_c = &zero,
        .oxygen_satisfaction_fraction = &one,
        .non_band_nitrate_g_n = &.{10},
        .band_nitrate_g_n = &.{10},
        .nitrate_half_saturation_g_n_per_m3 = &one,
        .nitrite_half_saturation_g_n_per_m3 = &one,
        .product_inhibition_g_n_per_m3_step = &one,
        .carbon_per_nitrate_n_g_c_per_g_n = &one,
        .nitrate_n_per_unmet_oxygen_g_n_per_g_o = 0.875,
        .minimum_competition_fraction = 0.01,
        .negligible_demand_g_n = 1e-12,
        .negligible_water_m3 = 1e-12,
        .timestep_h = 1,
    });
    try std.testing.expectEqual(@as(f64, 0.25), state.non_band_competition_fraction[0]);
    try std.testing.expectEqual(@as(f64, 0.25), state.band_competition_fraction[0]);
    try std.testing.expect(state.total_nitrate_reduction_g_n[0] > 0);
    try std.testing.expect(state.total_nitrate_reduction_g_n[0] <=
        state.total_supply_capacity_reduction_g_n[0]);
}
