const std = @import("std");
const nitrogen_parameters = @import("soil_nitrogen_parameters.zig");

pub const ParameterArrays = struct {
    enabled: []bool,
    acetate_product_inhibition_g_c_per_m3: []f64,
    acetate_half_saturation_g_c_per_m3: []f64,
    specific_respiration_rate_g_c_per_g_c_h: []f64,
    reference_energy_yield_kj_per_g_c: []f64,
    growth_energy_requirement_kj_per_g_c: []f64,
    minimum_growth_respiration_fraction: []f64,
};

/// Expands the runtime acetotrophic role and scalar NITRO parameters into the
/// vector-kernel boundary without assuming a fixed population count.
pub fn fillParameterArrays(
    parameters: nitrogen_parameters.AcetotrophicMethanogenesisParameters,
    population_index: []const usize,
    output: ParameterArrays,
) !void {
    const count = population_index.len;
    inline for (@typeInfo(ParameterArrays).@"struct".fields) |field|
        if (@field(output, field.name).len != count)
            return error.InvalidAcetotrophicMethanogenesisDimensions;
    for (population_index, 0..) |population, unit| {
        output.enabled[unit] = population == parameters.population_index;
        output.acetate_product_inhibition_g_c_per_m3[unit] = parameters.acetate_product_inhibition_g_c_per_m3;
        output.acetate_half_saturation_g_c_per_m3[unit] = parameters.acetate_half_saturation_g_c_per_m3;
        output.specific_respiration_rate_g_c_per_g_c_h[unit] = parameters.specific_respiration_rate_g_c_per_g_c_h;
        output.reference_energy_yield_kj_per_g_c[unit] = parameters.reference_energy_yield_kj_per_g_c;
        output.growth_energy_requirement_kj_per_g_c[unit] = parameters.growth_energy_requirement_kj_per_g_c;
        output.minimum_growth_respiration_fraction[unit] = parameters.minimum_growth_respiration_fraction;
    }
}

pub const Inputs = struct {
    enabled: []const bool,
    soil_temperature_k: []const f64,
    aqueous_acetate_concentration_g_c_per_m3: []const f64,
    aqueous_acetate_g_c: []const f64,
    acetate_competition_fraction: []const f64,
    combined_nutrient_limitation_fraction: []const f64,
    water_response: []const f64,
    active_biomass_g_c: []const f64,
    growth_temperature_response: []const f64,
    acetate_product_inhibition_g_c_per_m3: []const f64,
    acetate_half_saturation_g_c_per_m3: []const f64,
    specific_respiration_rate_g_c_per_g_c_h: []const f64,
    reference_energy_yield_kj_per_g_c: []const f64,
    growth_energy_requirement_kj_per_g_c: []const f64,
    minimum_growth_respiration_fraction: []const f64,
    minimum_acetate_concentration_g_c_per_m3: f64,
    gas_constant_kj_per_mol_k: f64,
    feedback_carbon_conversion_g_c_per_mol: f64,
    methane_carbon_yield_g_c_per_g_c_oxidized: f64,
    timestep_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    acetate_feedback_energy_kj_per_mol: []f64,
    acetate_feedback_energy_kj_per_g_c: []f64,
    growth_respiration_fraction: []f64,
    acetate_monod_fraction: []f64,
    substrate_unlimited_oxidation_g_c: []f64,
    microbial_limited_oxidation_g_c: []f64,
    supply_limited_oxidation_g_c: []f64,
    acetate_oxidation_g_c: []f64,
    potential_acetate_demand_g_c: []f64,
    methane_production_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0)
            return error.InvalidAcetotrophicMethanogenesisDimensions;
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

/// NITRO.F 986--1048. The historical N=5 admission is a runtime mask;
/// therefore any user-sized population catalog can designate acetotrophic
/// methanogens without inheriting the seven-population ceiling.
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
        const concentration_ratio =
            @max(
                inputs.minimum_acetate_concentration_g_c_per_m3,
                inputs.aqueous_acetate_concentration_g_c_per_m3[unit],
            ) / inputs.acetate_product_inhibition_g_c_per_m3[unit];
        const feedback_mol =
            inputs.gas_constant_kj_per_mol_k *
            inputs.soil_temperature_k[unit] *
            @log(concentration_ratio);
        const feedback_g_c =
            feedback_mol / inputs.feedback_carbon_conversion_g_c_per_mol;
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
        const acetate_concentration =
            inputs.aqueous_acetate_concentration_g_c_per_m3[unit];
        const monod =
            acetate_concentration /
            (acetate_concentration +
                inputs.acetate_half_saturation_g_c_per_m3[unit]);
        const unlimited = @max(
            0,
            inputs.specific_respiration_rate_g_c_per_g_c_h[unit] *
                inputs.combined_nutrient_limitation_fraction[unit] *
                inputs.water_response[unit] *
                inputs.active_biomass_g_c[unit] *
                inputs.timestep_h,
        );
        const microbial_limit =
            unlimited * monod * inputs.growth_temperature_response[unit];
        const supply_limit = @max(
            0,
            inputs.aqueous_acetate_g_c[unit] *
                inputs.acetate_competition_fraction[unit] *
                respiration_fraction *
                inputs.timestep_h,
        );
        const oxidation = @min(supply_limit, microbial_limit);

        state.acetate_feedback_energy_kj_per_mol[unit] = feedback_mol;
        state.acetate_feedback_energy_kj_per_g_c[unit] = feedback_g_c;
        state.growth_respiration_fraction[unit] = respiration_fraction;
        state.acetate_monod_fraction[unit] = monod;
        state.substrate_unlimited_oxidation_g_c[unit] = unlimited;
        state.microbial_limited_oxidation_g_c[unit] = microbial_limit;
        state.supply_limited_oxidation_g_c[unit] = supply_limit;
        state.acetate_oxidation_g_c[unit] = oxidation;
        state.potential_acetate_demand_g_c[unit] = microbial_limit;
        state.methane_production_g_c[unit] =
            inputs.methane_carbon_yield_g_c_per_g_c_oxidized * oxidation;
    }
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state.*, field.name)) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteAcetotrophicMethanogenesisResult;
            const signed_feedback =
                std.mem.eql(u8, field.name, "acetate_feedback_energy_kj_per_mol") or
                std.mem.eql(u8, field.name, "acetate_feedback_energy_kj_per_g_c");
            if (!signed_feedback and value < 0)
                return error.NonFiniteAcetotrophicMethanogenesisResult;
        };
}

pub fn sourceEnabled(zero_based_population: usize) bool {
    return zero_based_population == 4;
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.enabled.len != n)
        return error.InvalidAcetotrophicMethanogenesisDimensions;
    inline for (.{
        inputs.soil_temperature_k,
        inputs.aqueous_acetate_concentration_g_c_per_m3,
        inputs.aqueous_acetate_g_c,
        inputs.acetate_competition_fraction,
        inputs.combined_nutrient_limitation_fraction,
        inputs.water_response,
        inputs.active_biomass_g_c,
        inputs.growth_temperature_response,
        inputs.acetate_product_inhibition_g_c_per_m3,
        inputs.acetate_half_saturation_g_c_per_m3,
        inputs.specific_respiration_rate_g_c_per_g_c_h,
        inputs.reference_energy_yield_kj_per_g_c,
        inputs.growth_energy_requirement_kj_per_g_c,
        inputs.minimum_growth_respiration_fraction,
    }) |values| {
        if (values.len != n)
            return error.InvalidAcetotrophicMethanogenesisDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidAcetotrophicMethanogenesisInput;
    }
    for (0..n) |unit| if (inputs.enabled[unit]) {
        inline for (.{
            inputs.aqueous_acetate_concentration_g_c_per_m3[unit],
            inputs.aqueous_acetate_g_c[unit],
            inputs.acetate_competition_fraction[unit],
            inputs.combined_nutrient_limitation_fraction[unit],
            inputs.water_response[unit],
            inputs.active_biomass_g_c[unit],
            inputs.growth_temperature_response[unit],
            inputs.specific_respiration_rate_g_c_per_g_c_h[unit],
            inputs.minimum_growth_respiration_fraction[unit],
        }) |value| if (value < 0)
            return error.InvalidAcetotrophicMethanogenesisInput;
        if (inputs.soil_temperature_k[unit] <= 0 or
            inputs.combined_nutrient_limitation_fraction[unit] > 1 or
            inputs.acetate_product_inhibition_g_c_per_m3[unit] <= 0 or
            inputs.acetate_half_saturation_g_c_per_m3[unit] <= 0 or
            inputs.growth_energy_requirement_kj_per_g_c[unit] <= 0 or
            inputs.minimum_growth_respiration_fraction[unit] > 1)
            return error.InvalidAcetotrophicMethanogenesisInput;
    };
    inline for (.{
        inputs.minimum_acetate_concentration_g_c_per_m3,
        inputs.gas_constant_kj_per_mol_k,
        inputs.feedback_carbon_conversion_g_c_per_mol,
        inputs.methane_carbon_yield_g_c_per_g_c_oxidized,
        inputs.timestep_h,
    }) |value| if (!std.math.isFinite(value) or value <= 0)
        return error.InvalidAcetotrophicMethanogenesisInput;
}

fn clear(state: *State) void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) @memset(@field(state, field.name), 0);
}

test "source acetotrophic role reproduces NITRO selector" {
    try std.testing.expect(sourceEnabled(4));
    try std.testing.expect(!sourceEnabled(3));
    try std.testing.expect(!sourceEnabled(5));
}

test "runtime owner expands role and drives source acetotrophic calculation" {
    const parameters = (try nitrogen_parameters.sourceAnaerobicEnergyParameters()).acetotrophic_methanogenesis;
    const populations = [_]usize{ 4, 9 };
    var enabled: [2]bool = undefined;
    var inhibition: [2]f64 = undefined;
    var half_saturation: [2]f64 = undefined;
    var rate: [2]f64 = undefined;
    var reference: [2]f64 = undefined;
    var growth: [2]f64 = undefined;
    var minimum: [2]f64 = undefined;
    try fillParameterArrays(parameters, &populations, .{
        .enabled = &enabled,
        .acetate_product_inhibition_g_c_per_m3 = &inhibition,
        .acetate_half_saturation_g_c_per_m3 = &half_saturation,
        .specific_respiration_rate_g_c_per_g_c_h = &rate,
        .reference_energy_yield_kj_per_g_c = &reference,
        .growth_energy_requirement_kj_per_g_c = &growth,
        .minimum_growth_respiration_fraction = &minimum,
    });
    try std.testing.expectEqual([2]bool{ true, false }, enabled);

    var state = try State.init(std.testing.allocator, populations.len);
    defer state.deinit();
    const ones = [_]f64{ 1, 1 };
    const temperatures = [_]f64{ 300, 300 };
    const acetate_concentration = [_]f64{ 12, 12 };
    const acetate_mass = [_]f64{ 10, 10 };
    try calculate(&state, .{
        .enabled = &enabled,
        .soil_temperature_k = &temperatures,
        .aqueous_acetate_concentration_g_c_per_m3 = &acetate_concentration,
        .aqueous_acetate_g_c = &acetate_mass,
        .acetate_competition_fraction = &ones,
        .combined_nutrient_limitation_fraction = &ones,
        .water_response = &ones,
        .active_biomass_g_c = &ones,
        .growth_temperature_response = &ones,
        .acetate_product_inhibition_g_c_per_m3 = &inhibition,
        .acetate_half_saturation_g_c_per_m3 = &half_saturation,
        .specific_respiration_rate_g_c_per_g_c_h = &rate,
        .reference_energy_yield_kj_per_g_c = &reference,
        .growth_energy_requirement_kj_per_g_c = &growth,
        .minimum_growth_respiration_fraction = &minimum,
        .minimum_acetate_concentration_g_c_per_m3 = 1e-12,
        .gas_constant_kj_per_mol_k = parameters.gas_constant_kj_per_mol_k,
        .feedback_carbon_conversion_g_c_per_mol = parameters.feedback_carbon_conversion_g_c_per_mol,
        .methane_carbon_yield_g_c_per_g_c_oxidized = parameters.methane_carbon_yield_g_c_per_g_c_oxidized,
        .timestep_h = 1,
    });
    try std.testing.expect(state.acetate_oxidation_g_c[0] > 0);
    try std.testing.expectEqual(
        parameters.methane_carbon_yield_g_c_per_g_c_oxidized * state.acetate_oxidation_g_c[0],
        state.methane_production_g_c[0],
    );
    try std.testing.expectEqual(@as(f64, 0), state.acetate_oxidation_g_c[1]);
}

test "acetate thermodynamics kinetics and methane yield reproduce NITRO" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, .{
        .enabled = &.{true},
        .soil_temperature_k = &.{300},
        .aqueous_acetate_concentration_g_c_per_m3 = &.{2},
        .aqueous_acetate_g_c = &.{4},
        .acetate_competition_fraction = &.{0.5},
        .combined_nutrient_limitation_fraction = &.{0.8},
        .water_response = &.{0.5},
        .active_biomass_g_c = &.{10},
        .growth_temperature_response = &.{2},
        .acetate_product_inhibition_g_c_per_m3 = &.{1},
        .acetate_half_saturation_g_c_per_m3 = &.{1},
        .specific_respiration_rate_g_c_per_g_c_h = &.{0.2},
        .reference_energy_yield_kj_per_g_c = &.{1},
        .growth_energy_requirement_kj_per_g_c = &.{4},
        .minimum_growth_respiration_fraction = &.{0.2},
        .minimum_acetate_concentration_g_c_per_m3 = 1e-12,
        .gas_constant_kj_per_mol_k = 8.3143e-3,
        .feedback_carbon_conversion_g_c_per_mol = 24,
        .methane_carbon_yield_g_c_per_g_c_oxidized = 0.5,
        .timestep_h = 1,
    });
    const feedback_mol = 8.3143e-3 * 300 * @log(@as(f64, 2));
    const requirement = @max(
        @as(f64, 0.2),
        @min(@as(f64, 1), 1 / (1 + @max(0, 1 + feedback_mol / 24) / 4)),
    );
    try std.testing.expectApproxEqAbs(feedback_mol, state.acetate_feedback_energy_kj_per_mol[0], 1e-15);
    try std.testing.expectApproxEqAbs(requirement, state.growth_respiration_fraction[0], 1e-15);
    // Unlimited=.8, microbial=.8*(2/3)*2.
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), state.substrate_unlimited_oxidation_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 16.0 / 15.0), state.microbial_limited_oxidation_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@min(2 * requirement, 16.0 / 15.0), state.acetate_oxidation_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(0.5 * state.acetate_oxidation_g_c[0], state.methane_production_g_c[0], 1e-15);
}

test "invalid late methanogen leaves prior result unchanged" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.methane_production_g_c[0] = 7;
    const zeros = [_]f64{ 0, 0 };
    const ones = [_]f64{ 1, 1 };
    try std.testing.expectError(error.InvalidAcetotrophicMethanogenesisInput, calculate(&state, .{
        .enabled = &.{ true, true },
        .soil_temperature_k = &.{ 300, -1 },
        .aqueous_acetate_concentration_g_c_per_m3 = &zeros,
        .aqueous_acetate_g_c = &zeros,
        .acetate_competition_fraction = &ones,
        .combined_nutrient_limitation_fraction = &ones,
        .water_response = &ones,
        .active_biomass_g_c = &ones,
        .growth_temperature_response = &ones,
        .acetate_product_inhibition_g_c_per_m3 = &ones,
        .acetate_half_saturation_g_c_per_m3 = &ones,
        .specific_respiration_rate_g_c_per_g_c_h = &ones,
        .reference_energy_yield_kj_per_g_c = &ones,
        .growth_energy_requirement_kj_per_g_c = &ones,
        .minimum_growth_respiration_fraction = &ones,
        .minimum_acetate_concentration_g_c_per_m3 = 1e-12,
        .gas_constant_kj_per_mol_k = 8.3143e-3,
        .feedback_carbon_conversion_g_c_per_mol = 24,
        .methane_carbon_yield_g_c_per_g_c_oxidized = 0.5,
        .timestep_h = 1,
    }));
    try std.testing.expectEqual(@as(f64, 7), state.methane_production_g_c[0]);
}

test "NITRO 986-1048 derived overflow leaves methanogenesis state unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.methane_production_g_c[0] = 7;
    const one = [_]f64{1};
    try std.testing.expectError(
        error.NonFiniteAcetotrophicMethanogenesisResult,
        calculate(&state, .{
            .enabled = &.{true},
            .soil_temperature_k = &.{300},
            .aqueous_acetate_concentration_g_c_per_m3 = &.{2},
            .aqueous_acetate_g_c = &one,
            .acetate_competition_fraction = &one,
            .combined_nutrient_limitation_fraction = &one,
            .water_response = &one,
            .active_biomass_g_c = &one,
            .growth_temperature_response = &one,
            .acetate_product_inhibition_g_c_per_m3 = &one,
            .acetate_half_saturation_g_c_per_m3 = &one,
            .specific_respiration_rate_g_c_per_g_c_h = &one,
            .reference_energy_yield_kj_per_g_c = &one,
            .growth_energy_requirement_kj_per_g_c = &one,
            .minimum_growth_respiration_fraction = &one,
            .minimum_acetate_concentration_g_c_per_m3 = 1e-12,
            .gas_constant_kj_per_mol_k = std.math.floatMax(f64),
            .feedback_carbon_conversion_g_c_per_mol = 24,
            .methane_carbon_yield_g_c_per_g_c_oxidized = 0.5,
            .timestep_h = 1,
        }),
    );
    try std.testing.expectEqual(@as(f64, 7), state.methane_production_g_c[0]);
}
