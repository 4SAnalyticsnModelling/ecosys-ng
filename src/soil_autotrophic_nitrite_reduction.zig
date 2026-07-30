const std = @import("std");

pub const Inputs = struct {
    enabled: []const bool,
    previous_total_non_band_nitrite_demand_g_n: []const f64,
    previous_total_band_nitrite_demand_g_n: []const f64,
    previous_non_band_nitrite_capacity_g_n: []const f64,
    previous_band_nitrite_capacity_g_n: []const f64,
    active_biomass_fraction: []const f64,
    non_band_nitrate_fraction: []const f64,
    band_nitrate_fraction: []const f64,
    non_band_nitrite_fraction: []const f64,
    band_nitrite_fraction: []const f64,
    oxygen_demand_g_o: []const f64,
    oxygen_uptake_g_o: []const f64,
    aqueous_carbon_dioxide_limitation_fraction: []const f64,
    non_band_nitrite_concentration_g_n_per_m3: []const f64,
    band_nitrite_concentration_g_n_per_m3: []const f64,
    biologically_active_water_m3: []const f64,
    non_band_nitrite_g_n: []const f64,
    band_nitrite_g_n: []const f64,
    preceding_non_band_ammonia_oxidation_g_n: []const f64,
    preceding_band_ammonia_oxidation_g_n: []const f64,
    non_band_ammonium_competition_fraction: []const f64,
    band_ammonium_competition_fraction: []const f64,
    non_band_ammonium_g_n: []const f64,
    band_ammonium_g_n: []const f64,
    nitrite_half_saturation_g_n_per_m3: []const f64,
    product_inhibition_g_n_per_m3_step: []const f64,
    carbon_efficiency_g_c_per_g_n: []const f64,
    anaerobic_growth_respiration_fraction: []const f64,
    nitrate_n_per_unmet_oxygen_g_n_per_g_o: f64,
    ammonium_supply_per_nitrite_reduction: f64,
    additional_ammonium_oxidation_per_nitrite_reduction: f64,
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
    nitrite_reduction_demand_g_n: []f64,
    raw_non_band_capacity_g_n: []f64,
    raw_band_capacity_g_n: []f64,
    raw_total_capacity_g_n: []f64,
    product_inhibition_fraction: []f64,
    non_band_capacity_g_n: []f64,
    band_capacity_g_n: []f64,
    non_band_nitrite_supply_g_n: []f64,
    band_nitrite_supply_g_n: []f64,
    non_band_ammonium_donor_limit_g_n: []f64,
    band_ammonium_donor_limit_g_n: []f64,
    non_band_nitrite_reduction_g_n: []f64,
    band_nitrite_reduction_g_n: []f64,
    total_nitrite_reduction_g_n: []f64,
    supply_capacity_respiration_g_c: []f64,
    nitrite_reduction_respiration_g_c: []f64,
    nitrate_reduction_g_n: []f64,
    band_nitrate_reduction_g_n: []f64,
    nitrous_oxide_reduction_g_n: []f64,
    updated_non_band_ammonia_oxidation_g_n: []f64,
    updated_band_ammonia_oxidation_g_n: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0)
            return error.InvalidAutotrophicNitriteReductionDimensions;
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

/// NITRO.F 1955--2041. Historical K=5,N=1 admission is runtime data.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.unit_count);
    defer staged.deinit();
    for (0..state.unit_count) |unit| {
        if (!inputs.enabled[unit]) continue;
        const non_band_competition = competition(
            inputs.previous_total_non_band_nitrite_demand_g_n[unit],
            inputs.previous_non_band_nitrite_capacity_g_n[unit],
            inputs.active_biomass_fraction[unit] *
                inputs.non_band_nitrate_fraction[unit],
            inputs,
        );
        const band_competition = competition(
            inputs.previous_total_band_nitrite_demand_g_n[unit],
            inputs.previous_band_nitrite_capacity_g_n[unit],
            inputs.active_biomass_fraction[unit] *
                inputs.band_nitrate_fraction[unit],
            inputs,
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

fn competition(total: f64, capacity: f64, fallback: f64, inputs: Inputs) f64 {
    return @max(
        inputs.minimum_competition_fraction,
        if (total > inputs.negligible_demand_g_n) capacity / total else fallback,
    );
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
        inputs.nitrate_n_per_unmet_oxygen_g_n_per_g_o *
        unmet_oxygen *
        inputs.aqueous_carbon_dioxide_limitation_fraction[unit];
    const raw_non_band =
        inputs.non_band_nitrite_fraction[unit] * demand *
        inputs.non_band_nitrite_concentration_g_n_per_m3[unit] /
        (inputs.non_band_nitrite_concentration_g_n_per_m3[unit] +
            inputs.nitrite_half_saturation_g_n_per_m3[unit]);
    const raw_band =
        inputs.band_nitrite_fraction[unit] * demand *
        inputs.band_nitrite_concentration_g_n_per_m3[unit] /
        (inputs.band_nitrite_concentration_g_n_per_m3[unit] +
            inputs.nitrite_half_saturation_g_n_per_m3[unit]);
    const raw_total = raw_non_band + raw_band;
    const inhibition = if (inputs.biologically_active_water_m3[unit] >
        inputs.negligible_water_m3)
        1 / (1 + raw_total /
            (inputs.product_inhibition_g_n_per_m3_step[unit] *
                inputs.biologically_active_water_m3[unit]))
    else
        0;
    const non_band_capacity = raw_non_band * inhibition;
    const band_capacity = raw_band * inhibition;
    const non_band_supply =
        inputs.non_band_nitrite_g_n[unit] +
        inputs.preceding_non_band_ammonia_oxidation_g_n[unit];
    const band_supply =
        inputs.band_nitrite_g_n[unit] +
        inputs.preceding_band_ammonia_oxidation_g_n[unit];
    const non_band_donor =
        inputs.non_band_ammonium_competition_fraction[unit] *
        inputs.ammonium_supply_per_nitrite_reduction *
        inputs.non_band_ammonium_g_n[unit] *
        inputs.timestep_h;
    const band_donor =
        inputs.band_ammonium_competition_fraction[unit] *
        inputs.ammonium_supply_per_nitrite_reduction *
        inputs.band_ammonium_g_n[unit] *
        inputs.timestep_h;
    const non_band_reduction =
        @max(0, @min(non_band_capacity, non_band_supply, non_band_donor));
    const band_reduction =
        @max(0, @min(band_capacity, band_supply, band_donor));
    const total = non_band_reduction + band_reduction;
    state.non_band_competition_fraction[unit] = non_band_competition;
    state.band_competition_fraction[unit] = band_competition;
    state.unmet_oxygen_demand_g_o[unit] = unmet_oxygen;
    state.nitrite_reduction_demand_g_n[unit] = demand;
    state.raw_non_band_capacity_g_n[unit] = raw_non_band;
    state.raw_band_capacity_g_n[unit] = raw_band;
    state.raw_total_capacity_g_n[unit] = raw_total;
    state.product_inhibition_fraction[unit] = inhibition;
    state.non_band_capacity_g_n[unit] = non_band_capacity;
    state.band_capacity_g_n[unit] = band_capacity;
    state.non_band_nitrite_supply_g_n[unit] = non_band_supply;
    state.band_nitrite_supply_g_n[unit] = band_supply;
    state.non_band_ammonium_donor_limit_g_n[unit] = non_band_donor;
    state.band_ammonium_donor_limit_g_n[unit] = band_donor;
    state.non_band_nitrite_reduction_g_n[unit] = non_band_reduction;
    state.band_nitrite_reduction_g_n[unit] = band_reduction;
    state.total_nitrite_reduction_g_n[unit] = total;
    state.supply_capacity_respiration_g_c[unit] = 0;
    state.nitrite_reduction_respiration_g_c[unit] =
        total * inputs.carbon_efficiency_g_c_per_g_n[unit] *
        inputs.anaerobic_growth_respiration_fraction[unit];
    state.nitrate_reduction_g_n[unit] = 0;
    state.band_nitrate_reduction_g_n[unit] = 0;
    state.nitrous_oxide_reduction_g_n[unit] = 0;
    state.updated_non_band_ammonia_oxidation_g_n[unit] =
        inputs.preceding_non_band_ammonia_oxidation_g_n[unit] +
        inputs.additional_ammonium_oxidation_per_nitrite_reduction *
            non_band_reduction;
    state.updated_band_ammonia_oxidation_g_n[unit] =
        inputs.preceding_band_ammonia_oxidation_g_n[unit] +
        inputs.additional_ammonium_oxidation_per_nitrite_reduction *
            band_reduction;
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.enabled.len != n)
        return error.InvalidAutotrophicNitriteReductionDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == []const f64) {
            const values = @field(inputs, field.name);
            if (values.len != n)
                return error.InvalidAutotrophicNitriteReductionDimensions;
            for (values) |value| if (!std.math.isFinite(value) or value < 0)
                return error.InvalidAutotrophicNitriteReductionInput;
        } else if (field.type == f64) {
            const value = @field(inputs, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidAutotrophicNitriteReductionInput;
        }
    }
    for (0..n) |unit| if (inputs.enabled[unit]) {
        if (inputs.nitrite_half_saturation_g_n_per_m3[unit] <= 0 or
            inputs.product_inhibition_g_n_per_m3_step[unit] <= 0)
            return error.InvalidAutotrophicNitriteReductionInput;
    };
    if (inputs.nitrate_n_per_unmet_oxygen_g_n_per_g_o <= 0 or
        inputs.ammonium_supply_per_nitrite_reduction <= 0 or
        inputs.additional_ammonium_oxidation_per_nitrite_reduction <= 0 or
        inputs.timestep_h <= 0)
        return error.InvalidAutotrophicNitriteReductionInput;
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (field.type == f64) {
            const value = @field(state.*, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteAutotrophicNitriteReductionResult;
        } else if (field.type == []f64) {
            for (@field(state.*, field.name)) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteAutotrophicNitriteReductionResult;
        }
    }
}

test "NITRO autotrophic denitrification retains distinct 3 and 0.333 factors" {
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
        .active_biomass_fraction = &one,
        .non_band_nitrate_fraction = &.{0.5},
        .band_nitrate_fraction = &.{0.5},
        .non_band_nitrite_fraction = &.{0.5},
        .band_nitrite_fraction = &.{0.5},
        .oxygen_demand_g_o = &one,
        .oxygen_uptake_g_o = &zero,
        .aqueous_carbon_dioxide_limitation_fraction = &one,
        .non_band_nitrite_concentration_g_n_per_m3 = &one,
        .band_nitrite_concentration_g_n_per_m3 = &one,
        .biologically_active_water_m3 = &one,
        .non_band_nitrite_g_n = &one,
        .band_nitrite_g_n = &one,
        .preceding_non_band_ammonia_oxidation_g_n = &zero,
        .preceding_band_ammonia_oxidation_g_n = &zero,
        .non_band_ammonium_competition_fraction = &one,
        .band_ammonium_competition_fraction = &one,
        .non_band_ammonium_g_n = &one,
        .band_ammonium_g_n = &one,
        .nitrite_half_saturation_g_n_per_m3 = &one,
        .product_inhibition_g_n_per_m3_step = &one,
        .carbon_efficiency_g_c_per_g_n = &one,
        .anaerobic_growth_respiration_fraction = &one,
        .nitrate_n_per_unmet_oxygen_g_n_per_g_o = 0.875,
        .ammonium_supply_per_nitrite_reduction = 3,
        .additional_ammonium_oxidation_per_nitrite_reduction = 0.333,
        .minimum_competition_fraction = 0.01,
        .negligible_demand_g_n = 1e-12,
        .negligible_water_m3 = 1e-12,
        .timestep_h = 1,
    });
    try std.testing.expectEqual(
        0.333 * state.non_band_nitrite_reduction_g_n[0],
        state.updated_non_band_ammonia_oxidation_g_n[0],
    );
    try std.testing.expectEqual(@as(f64, 0), state.nitrate_reduction_g_n[0]);
    try std.testing.expect(state.non_band_ammonium_donor_limit_g_n[0] >
        state.updated_non_band_ammonia_oxidation_g_n[0]);
}
