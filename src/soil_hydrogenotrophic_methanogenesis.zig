const std = @import("std");

pub const Inputs = struct {
    enabled: []const bool,
    hydrogen_feedback_energy_kj_per_mol: []const f64,
    aqueous_hydrogen_concentration_g_h_per_m3: []const f64,
    aqueous_hydrogen_g_h: []const f64,
    fermentation_hydrogen_production_equivalent_g_c: []const f64,
    temperature_water_response: []const f64,
    combined_nutrient_limitation_fraction: []const f64,
    aqueous_carbon_dioxide_limitation_fraction: []const f64,
    active_biomass_g_c: []const f64,
    specific_carbon_dioxide_reduction_g_c_per_g_c_h: []const f64,
    hydrogen_half_saturation_g_h_per_m3: []const f64,
    reference_energy_yield_kj_per_g_c: []const f64,
    growth_energy_requirement_kj_per_g_c: []const f64,
    minimum_growth_respiration_fraction: []const f64,
    hydrogen_feedback_molar_to_carbon_divisor: f64,
    fermentation_carbon_to_hydrogen_g_h_per_g_c: f64,
    carbon_dioxide_reduction_per_hydrogen_g_c_per_g_h: f64,
    methane_carbon_per_reduced_carbon_g_c_per_g_c: f64,
    timestep_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    hydrogen_feedback_energy_kj_per_g_c: []f64,
    growth_respiration_fraction: []f64,
    hydrogen_unlimited_carbon_dioxide_reduction_g_c: []f64,
    total_available_hydrogen_g_h: []f64,
    hydrogen_monod_fraction: []f64,
    kinetic_limited_carbon_dioxide_reduction_g_c: []f64,
    supply_limited_carbon_dioxide_reduction_g_c: []f64,
    carbon_dioxide_reduction_g_c: []f64,
    methane_production_g_c: []f64,
    oxygen_demand_g_o: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0)
            return error.InvalidHydrogenotrophicMethanogenesisDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
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

/// NITRO.F 1287--1327. Historical K=5,N=5 admission is runtime data. The
/// source conversions 1/12, 0.111 and 1.5 are named runtime coefficients.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.unit_count);
    defer staged.deinit();
    calculateValidated(&staged, inputs);
    try validateResult(&staged);
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64)
            @memcpy(@field(state, field.name), @field(staged, field.name));
}

fn calculateValidated(state: *State, inputs: Inputs) void {
    clear(state);
    for (0..state.unit_count) |unit| {
        if (!inputs.enabled[unit]) continue;
        const feedback_g_c =
            inputs.hydrogen_feedback_energy_kj_per_mol[unit] /
            inputs.hydrogen_feedback_molar_to_carbon_divisor;
        const respiration_fraction = @max(
            inputs.minimum_growth_respiration_fraction[unit],
            @min(
                1,
                1 /
                    (1 +
                        @max(
                            0,
                            inputs.reference_energy_yield_kj_per_g_c[unit] +
                                feedback_g_c,
                        ) /
                            inputs.growth_energy_requirement_kj_per_g_c[unit]),
            ),
        );
        const unlimited =
            inputs.specific_carbon_dioxide_reduction_g_c_per_g_c_h[unit] *
            inputs.temperature_water_response[unit] *
            inputs.combined_nutrient_limitation_fraction[unit] *
            inputs.aqueous_carbon_dioxide_limitation_fraction[unit] *
            inputs.active_biomass_g_c[unit] *
            inputs.timestep_h;
        const total_hydrogen =
            inputs.aqueous_hydrogen_g_h[unit] +
            inputs.fermentation_carbon_to_hydrogen_g_h_per_g_c *
                inputs.fermentation_hydrogen_production_equivalent_g_c[unit];
        const concentration =
            inputs.aqueous_hydrogen_concentration_g_h_per_m3[unit];
        const monod =
            concentration /
            (concentration +
                inputs.hydrogen_half_saturation_g_h_per_m3[unit]);
        const kinetic_limit = unlimited * monod;
        const supply_limit =
            inputs.carbon_dioxide_reduction_per_hydrogen_g_c_per_g_h *
            total_hydrogen *
            inputs.timestep_h;
        const reduction = @max(0, @min(kinetic_limit, supply_limit));

        state.hydrogen_feedback_energy_kj_per_g_c[unit] = feedback_g_c;
        state.growth_respiration_fraction[unit] = respiration_fraction;
        state.hydrogen_unlimited_carbon_dioxide_reduction_g_c[unit] =
            unlimited;
        state.total_available_hydrogen_g_h[unit] = total_hydrogen;
        state.hydrogen_monod_fraction[unit] = monod;
        state.kinetic_limited_carbon_dioxide_reduction_g_c[unit] =
            kinetic_limit;
        state.supply_limited_carbon_dioxide_reduction_g_c[unit] =
            supply_limit;
        state.carbon_dioxide_reduction_g_c[unit] = reduction;
        state.methane_production_g_c[unit] =
            inputs.methane_carbon_per_reduced_carbon_g_c_per_g_c *
            reduction;
        state.oxygen_demand_g_o[unit] = 0;
    }
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state.*, field.name)) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteHydrogenotrophicMethanogenesisResult;
            if (!std.mem.eql(
                u8,
                field.name,
                "hydrogen_feedback_energy_kj_per_g_c",
            ) and value < 0)
                return error.NonFiniteHydrogenotrophicMethanogenesisResult;
        };
}

pub fn sourceEnabled(
    zero_based_complex: usize,
    zero_based_population: usize,
) bool {
    return zero_based_complex == 5 and zero_based_population == 4;
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.enabled.len != n)
        return error.InvalidHydrogenotrophicMethanogenesisDimensions;
    inline for (.{
        inputs.hydrogen_feedback_energy_kj_per_mol,
        inputs.aqueous_hydrogen_concentration_g_h_per_m3,
        inputs.aqueous_hydrogen_g_h,
        inputs.fermentation_hydrogen_production_equivalent_g_c,
        inputs.temperature_water_response,
        inputs.combined_nutrient_limitation_fraction,
        inputs.aqueous_carbon_dioxide_limitation_fraction,
        inputs.active_biomass_g_c,
        inputs.specific_carbon_dioxide_reduction_g_c_per_g_c_h,
        inputs.hydrogen_half_saturation_g_h_per_m3,
        inputs.reference_energy_yield_kj_per_g_c,
        inputs.growth_energy_requirement_kj_per_g_c,
        inputs.minimum_growth_respiration_fraction,
    }) |values| {
        if (values.len != n)
            return error.InvalidHydrogenotrophicMethanogenesisDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidHydrogenotrophicMethanogenesisInput;
    }
    for (0..n) |unit| if (inputs.enabled[unit]) {
        inline for (.{
            inputs.aqueous_hydrogen_concentration_g_h_per_m3[unit],
            inputs.aqueous_hydrogen_g_h[unit],
            inputs.fermentation_hydrogen_production_equivalent_g_c[unit],
            inputs.temperature_water_response[unit],
            inputs.combined_nutrient_limitation_fraction[unit],
            inputs.aqueous_carbon_dioxide_limitation_fraction[unit],
            inputs.active_biomass_g_c[unit],
            inputs.specific_carbon_dioxide_reduction_g_c_per_g_c_h[unit],
            inputs.minimum_growth_respiration_fraction[unit],
        }) |value| if (value < 0)
            return error.InvalidHydrogenotrophicMethanogenesisInput;
        if (inputs.combined_nutrient_limitation_fraction[unit] > 1 or
            inputs.aqueous_carbon_dioxide_limitation_fraction[unit] > 1 or
            inputs.hydrogen_half_saturation_g_h_per_m3[unit] <= 0 or
            inputs.growth_energy_requirement_kj_per_g_c[unit] <= 0 or
            inputs.minimum_growth_respiration_fraction[unit] > 1)
            return error.InvalidHydrogenotrophicMethanogenesisInput;
    };
    inline for (.{
        inputs.hydrogen_feedback_molar_to_carbon_divisor,
        inputs.fermentation_carbon_to_hydrogen_g_h_per_g_c,
        inputs.carbon_dioxide_reduction_per_hydrogen_g_c_per_g_h,
        inputs.methane_carbon_per_reduced_carbon_g_c_per_g_c,
        inputs.timestep_h,
    }) |value| if (!std.math.isFinite(value) or value <= 0)
        return error.InvalidHydrogenotrophicMethanogenesisInput;
}

fn clear(state: *State) void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) @memset(@field(state, field.name), 0);
}

test "source hydrogenotrophic role reproduces NITRO selector" {
    try std.testing.expect(sourceEnabled(5, 4));
    try std.testing.expect(!sourceEnabled(5, 3));
    try std.testing.expect(!sourceEnabled(4, 4));
}

test "fermentation hydrogen and finite supply reproduce NITRO methanogenesis" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, .{
        .enabled = &.{true},
        .hydrogen_feedback_energy_kj_per_mol = &.{12},
        .aqueous_hydrogen_concentration_g_h_per_m3 = &.{1},
        .aqueous_hydrogen_g_h = &.{0.01},
        .fermentation_hydrogen_production_equivalent_g_c = &.{0.9},
        .temperature_water_response = &.{0.5},
        .combined_nutrient_limitation_fraction = &.{0.8},
        .aqueous_carbon_dioxide_limitation_fraction = &.{0.5},
        .active_biomass_g_c = &.{10},
        .specific_carbon_dioxide_reduction_g_c_per_g_c_h = &.{0.2},
        .hydrogen_half_saturation_g_h_per_m3 = &.{1},
        .reference_energy_yield_kj_per_g_c = &.{1},
        .growth_energy_requirement_kj_per_g_c = &.{4},
        .minimum_growth_respiration_fraction = &.{0.2},
        .hydrogen_feedback_molar_to_carbon_divisor = 12,
        .fermentation_carbon_to_hydrogen_g_h_per_g_c = 0.111,
        .carbon_dioxide_reduction_per_hydrogen_g_c_per_g_h = 1.5,
        .methane_carbon_per_reduced_carbon_g_c_per_g_c = 1,
        .timestep_h = 1,
    });
    const hydrogen = 0.01 + 0.111 * 0.9;
    try std.testing.expectApproxEqAbs(hydrogen, state.total_available_hydrogen_g_h[0], 1e-15);
    // VMXA=.4, Monod=.5 -> .2; H supply=1.5*0.1099=.16485.
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), state.hydrogen_unlimited_carbon_dioxide_reduction_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(1.5 * hydrogen, state.carbon_dioxide_reduction_g_c[0], 1e-15);
    try std.testing.expectEqual(state.carbon_dioxide_reduction_g_c[0], state.methane_production_g_c[0]);
    try std.testing.expectEqual(@as(f64, 0), state.oxygen_demand_g_o[0]);
}

test "invalid late hydrogenotroph preserves prior state" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.methane_production_g_c[0] = 7;
    const zeros = [_]f64{ 0, 0 };
    const ones = [_]f64{ 1, 1 };
    try std.testing.expectError(error.InvalidHydrogenotrophicMethanogenesisInput, calculate(&state, .{
        .enabled = &.{ true, true },
        .hydrogen_feedback_energy_kj_per_mol = &zeros,
        .aqueous_hydrogen_concentration_g_h_per_m3 = &zeros,
        .aqueous_hydrogen_g_h = &zeros,
        .fermentation_hydrogen_production_equivalent_g_c = &zeros,
        .temperature_water_response = &ones,
        .combined_nutrient_limitation_fraction = &ones,
        .aqueous_carbon_dioxide_limitation_fraction = &.{ 1, 2 },
        .active_biomass_g_c = &ones,
        .specific_carbon_dioxide_reduction_g_c_per_g_c_h = &ones,
        .hydrogen_half_saturation_g_h_per_m3 = &ones,
        .reference_energy_yield_kj_per_g_c = &ones,
        .growth_energy_requirement_kj_per_g_c = &ones,
        .minimum_growth_respiration_fraction = &ones,
        .hydrogen_feedback_molar_to_carbon_divisor = 12,
        .fermentation_carbon_to_hydrogen_g_h_per_g_c = 0.111,
        .carbon_dioxide_reduction_per_hydrogen_g_c_per_g_h = 1.5,
        .methane_carbon_per_reduced_carbon_g_c_per_g_c = 1,
        .timestep_h = 1,
    }));
    try std.testing.expectEqual(@as(f64, 7), state.methane_production_g_c[0]);
}

test "NITRO 1287-1327 derived overflow preserves hydrogenotroph state" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.methane_production_g_c[0] = 7;
    const zero = [_]f64{0};
    const one = [_]f64{1};
    try std.testing.expectError(
        error.NonFiniteHydrogenotrophicMethanogenesisResult,
        calculate(&state, .{
            .enabled = &.{true},
            .hydrogen_feedback_energy_kj_per_mol = &zero,
            .aqueous_hydrogen_concentration_g_h_per_m3 = &one,
            .aqueous_hydrogen_g_h = &one,
            .fermentation_hydrogen_production_equivalent_g_c = &zero,
            .temperature_water_response = &one,
            .combined_nutrient_limitation_fraction = &one,
            .aqueous_carbon_dioxide_limitation_fraction = &one,
            .active_biomass_g_c = &.{std.math.floatMax(f64)},
            .specific_carbon_dioxide_reduction_g_c_per_g_c_h = &.{
                std.math.floatMax(f64),
            },
            .hydrogen_half_saturation_g_h_per_m3 = &one,
            .reference_energy_yield_kj_per_g_c = &one,
            .growth_energy_requirement_kj_per_g_c = &one,
            .minimum_growth_respiration_fraction = &one,
            .hydrogen_feedback_molar_to_carbon_divisor = 12,
            .fermentation_carbon_to_hydrogen_g_h_per_g_c = 0.111,
            .carbon_dioxide_reduction_per_hydrogen_g_c_per_g_h = 1.5,
            .methane_carbon_per_reduced_carbon_g_c_per_g_c = 1,
            .timestep_h = 1,
        }),
    );
    try std.testing.expectEqual(@as(f64, 7), state.methane_production_g_c[0]);
}
