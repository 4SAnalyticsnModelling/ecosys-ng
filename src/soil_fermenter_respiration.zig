const std = @import("std");
const nitrogen_parameters = @import("soil_nitrogen_parameters.zig");

pub const ParameterArrays = struct {
    enabled: []bool,
    reference_energy_yield_kj_per_g_c: []f64,
    growth_energy_requirement_kj_per_g_c: []f64,
    minimum_respiration_requirement_g_c_per_g_c: []f64,
    specific_oxidation_rate_g_c_per_g_c_h: []f64,
    dissolved_organic_carbon_half_saturation_g_c_per_m3: []f64,
    acetate_product_inhibition_g_c_per_m3: []f64,
};

/// Expands runtime NITRO fermenter roles and scalar scientific parameters
/// into the existing vector-kernel boundary without a seven-population cap.
pub fn fillParameterArrays(
    parameters: nitrogen_parameters.FermenterRespirationParameters,
    population_index: []const usize,
    output: ParameterArrays,
) !void {
    const count = population_index.len;
    inline for (@typeInfo(ParameterArrays).@"struct".fields) |field|
        if (@field(output, field.name).len != count)
            return error.InvalidFermenterRespirationDimensions;
    for (population_index, 0..) |population, unit| {
        const is_fermenter = population == parameters.fermenter_population_index;
        const is_diazotroph = population == parameters.anaerobic_diazotroph_population_index;
        output.enabled[unit] = is_fermenter or is_diazotroph;
        output.reference_energy_yield_kj_per_g_c[unit] = parameters.reference_energy_yield_kj_per_g_c;
        output.growth_energy_requirement_kj_per_g_c[unit] = if (is_diazotroph)
            parameters.diazotroph_growth_energy_requirement_kj_per_g_c
        else
            parameters.growth_energy_requirement_kj_per_g_c;
        output.minimum_respiration_requirement_g_c_per_g_c[unit] = if (is_diazotroph)
            parameters.diazotroph_minimum_respiration_requirement_g_c_per_g_c
        else
            parameters.minimum_respiration_requirement_g_c_per_g_c;
        output.specific_oxidation_rate_g_c_per_g_c_h[unit] = parameters.specific_oxidation_rate_g_c_per_g_c_h;
        output.dissolved_organic_carbon_half_saturation_g_c_per_m3[unit] = parameters.dissolved_organic_carbon_half_saturation_g_c_per_m3;
        output.acetate_product_inhibition_g_c_per_m3[unit] = parameters.acetate_product_inhibition_g_c_per_m3;
    }
}

pub const Inputs = struct {
    enabled: []const bool,
    aqueous_acetate_concentration_g_c_per_m3: []const f64,
    aqueous_dissolved_organic_carbon_concentration_g_c_per_m3: []const f64,
    aqueous_dissolved_organic_carbon_g_c: []const f64,
    dissolved_organic_carbon_competition_fraction: []const f64,
    combined_nutrient_limitation_fraction: []const f64,
    water_response: []const f64,
    active_biomass_g_c: []const f64,
    growth_temperature_response: []const f64,
    fermentation_oxygen_inhibition_fraction: []const f64,
    soil_temperature_k: []const f64,
    hydrogen_feedback_energy_kj_per_mol: []const f64,
    reference_fermentation_energy_yield_kj_per_g_c: []const f64,
    growth_energy_requirement_kj_per_g_c: []const f64,
    minimum_respiration_requirement_g_c_per_g_c: []const f64,
    specific_oxidation_rate_g_c_per_g_c_h: []const f64,
    dissolved_organic_carbon_half_saturation_g_c_per_m3: []const f64,
    acetate_product_inhibition_g_c_per_m3: []const f64,
    minimum_acetate_concentration_g_c_per_m3: f64,
    gas_constant_kj_per_mol_k: f64,
    acetate_feedback_stoichiometric_exponent: f64,
    feedback_carbon_conversion_g_c_per_mol: f64,
    timestep_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    hydrogen_feedback_kj_per_g_c: []f64,
    acetate_feedback_energy_kj_per_mol: []f64,
    acetate_feedback_kj_per_g_c: []f64,
    combined_product_feedback_kj_per_g_c: []f64,
    respiration_requirement_g_c_per_g_c: []f64,
    dissolved_organic_carbon_substrate_response: []f64,
    substrate_unlimited_respiration_g_c: []f64,
    microbial_limited_respiration_g_c: []f64,
    supply_limited_respiration_g_c: []f64,
    actual_respiration_g_c: []f64,
    potential_dissolved_organic_carbon_demand_g_c: []f64,
    hydrogen_production_equivalent_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidFermenterRespirationDimensions;
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

/// NITRO.F 900--973. Source N=4 and N=7 differ only through their runtime
/// minimum respiration and growth-energy requirements, so no fixed role
/// indexing is required.
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
        const hydrogen_feedback =
            inputs.hydrogen_feedback_energy_kj_per_mol[unit] /
            inputs.feedback_carbon_conversion_g_c_per_mol;
        const acetate_ratio =
            @max(
                inputs.minimum_acetate_concentration_g_c_per_m3,
                inputs.aqueous_acetate_concentration_g_c_per_m3[unit],
            ) / inputs.acetate_product_inhibition_g_c_per_m3[unit];
        const acetate_feedback_mol =
            inputs.gas_constant_kj_per_mol_k *
            inputs.soil_temperature_k[unit] *
            @log(std.math.pow(
                f64,
                acetate_ratio,
                inputs.acetate_feedback_stoichiometric_exponent,
            ));
        const acetate_feedback = acetate_feedback_mol /
            inputs.feedback_carbon_conversion_g_c_per_mol;
        const combined_feedback = hydrogen_feedback + acetate_feedback;
        const respiration_requirement = @max(
            inputs.minimum_respiration_requirement_g_c_per_g_c[unit],
            @min(
                1,
                1 /
                    (1 +
                        @max(
                            0,
                            inputs.reference_fermentation_energy_yield_kj_per_g_c[unit] -
                                combined_feedback,
                        ) /
                            inputs.growth_energy_requirement_kj_per_g_c[unit]),
            ),
        );
        const doc_concentration =
            inputs.aqueous_dissolved_organic_carbon_concentration_g_c_per_m3[unit];
        const substrate_response =
            doc_concentration /
            (doc_concentration +
                inputs.dissolved_organic_carbon_half_saturation_g_c_per_m3[unit]) *
            inputs.fermentation_oxygen_inhibition_fraction[unit];
        const unlimited = @max(
            0,
            inputs.specific_oxidation_rate_g_c_per_g_c_h[unit] *
                inputs.combined_nutrient_limitation_fraction[unit] *
                inputs.water_response[unit] *
                inputs.active_biomass_g_c[unit] *
                inputs.timestep_h,
        );
        const microbial_limit =
            unlimited * substrate_response *
            inputs.growth_temperature_response[unit];
        const supply_limit = @max(
            0,
            inputs.aqueous_dissolved_organic_carbon_g_c[unit] *
                inputs.dissolved_organic_carbon_competition_fraction[unit] *
                respiration_requirement *
                inputs.timestep_h,
        );
        const actual = @min(supply_limit, microbial_limit);

        state.hydrogen_feedback_kj_per_g_c[unit] = hydrogen_feedback;
        state.acetate_feedback_energy_kj_per_mol[unit] =
            acetate_feedback_mol;
        state.acetate_feedback_kj_per_g_c[unit] = acetate_feedback;
        state.combined_product_feedback_kj_per_g_c[unit] =
            combined_feedback;
        state.respiration_requirement_g_c_per_g_c[unit] =
            respiration_requirement;
        state.dissolved_organic_carbon_substrate_response[unit] =
            substrate_response;
        state.substrate_unlimited_respiration_g_c[unit] = unlimited;
        state.microbial_limited_respiration_g_c[unit] = microbial_limit;
        state.supply_limited_respiration_g_c[unit] = supply_limit;
        state.actual_respiration_g_c[unit] = actual;
        state.potential_dissolved_organic_carbon_demand_g_c[unit] =
            microbial_limit;
        state.hydrogen_production_equivalent_g_c[unit] = actual;
    }
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state.*, field.name)) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteFermenterRespirationResult;
            const signed_feedback =
                std.mem.eql(u8, field.name, "hydrogen_feedback_kj_per_g_c") or
                std.mem.eql(u8, field.name, "acetate_feedback_energy_kj_per_mol") or
                std.mem.eql(u8, field.name, "acetate_feedback_kj_per_g_c") or
                std.mem.eql(u8, field.name, "combined_product_feedback_kj_per_g_c");
            if (!signed_feedback and
                value < 0)
                return error.NonFiniteFermenterRespirationResult;
        };
}

/// Exact source fermenter selector for zero-based population N.
pub fn sourceEnabled(zero_based_population: usize) bool {
    return zero_based_population == 3 or zero_based_population == 6;
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.enabled.len != n)
        return error.InvalidFermenterRespirationDimensions;
    inline for (.{
        inputs.aqueous_acetate_concentration_g_c_per_m3,
        inputs.aqueous_dissolved_organic_carbon_concentration_g_c_per_m3,
        inputs.aqueous_dissolved_organic_carbon_g_c,
        inputs.dissolved_organic_carbon_competition_fraction,
        inputs.combined_nutrient_limitation_fraction,
        inputs.water_response,
        inputs.active_biomass_g_c,
        inputs.growth_temperature_response,
        inputs.fermentation_oxygen_inhibition_fraction,
        inputs.soil_temperature_k,
        inputs.hydrogen_feedback_energy_kj_per_mol,
        inputs.reference_fermentation_energy_yield_kj_per_g_c,
        inputs.growth_energy_requirement_kj_per_g_c,
        inputs.minimum_respiration_requirement_g_c_per_g_c,
        inputs.specific_oxidation_rate_g_c_per_g_c_h,
        inputs.dissolved_organic_carbon_half_saturation_g_c_per_m3,
        inputs.acetate_product_inhibition_g_c_per_m3,
    }) |values| {
        if (values.len != n) return error.InvalidFermenterRespirationDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidFermenterRespirationInput;
    }
    for (0..n) |unit| if (inputs.enabled[unit]) {
        inline for (.{
            inputs.aqueous_acetate_concentration_g_c_per_m3[unit],
            inputs.aqueous_dissolved_organic_carbon_concentration_g_c_per_m3[unit],
            inputs.aqueous_dissolved_organic_carbon_g_c[unit],
            inputs.dissolved_organic_carbon_competition_fraction[unit],
            inputs.combined_nutrient_limitation_fraction[unit],
            inputs.water_response[unit],
            inputs.active_biomass_g_c[unit],
            inputs.growth_temperature_response[unit],
            inputs.fermentation_oxygen_inhibition_fraction[unit],
            inputs.reference_fermentation_energy_yield_kj_per_g_c[unit],
            inputs.minimum_respiration_requirement_g_c_per_g_c[unit],
            inputs.specific_oxidation_rate_g_c_per_g_c_h[unit],
        }) |value| if (value < 0)
            return error.InvalidFermenterRespirationInput;
        if (inputs.combined_nutrient_limitation_fraction[unit] > 1 or
            inputs.fermentation_oxygen_inhibition_fraction[unit] > 1 or
            inputs.soil_temperature_k[unit] <= 0 or
            inputs.growth_energy_requirement_kj_per_g_c[unit] <= 0 or
            inputs.minimum_respiration_requirement_g_c_per_g_c[unit] > 1 or
            inputs.dissolved_organic_carbon_half_saturation_g_c_per_m3[unit] <= 0 or
            inputs.acetate_product_inhibition_g_c_per_m3[unit] <= 0)
            return error.InvalidFermenterRespirationInput;
    };
    if (!std.math.isFinite(
        inputs.minimum_acetate_concentration_g_c_per_m3,
    ) or inputs.minimum_acetate_concentration_g_c_per_m3 <= 0 or
        !std.math.isFinite(inputs.gas_constant_kj_per_mol_k) or
        inputs.gas_constant_kj_per_mol_k <= 0 or
        !std.math.isFinite(inputs.acetate_feedback_stoichiometric_exponent) or
        inputs.acetate_feedback_stoichiometric_exponent <= 0 or
        !std.math.isFinite(inputs.feedback_carbon_conversion_g_c_per_mol) or
        inputs.feedback_carbon_conversion_g_c_per_mol <= 0 or
        !std.math.isFinite(inputs.timestep_h) or inputs.timestep_h <= 0)
        return error.InvalidFermenterRespirationInput;
}

fn clear(state: *State) void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) @memset(@field(state, field.name), 0);
}

test "source fermenter roles reproduce NITRO selector" {
    try std.testing.expect(sourceEnabled(3));
    try std.testing.expect(sourceEnabled(6));
    try std.testing.expect(!sourceEnabled(2));
    try std.testing.expect(!sourceEnabled(4));
}

test "runtime owner expands roles and drives source fermenter calculation" {
    const parameters = (try nitrogen_parameters.sourceAnaerobicEnergyParameters()).fermenter;
    const populations = [_]usize{ 3, 6, 11 };
    var enabled: [3]bool = undefined;
    var reference: [3]f64 = undefined;
    var growth: [3]f64 = undefined;
    var minimum: [3]f64 = undefined;
    var rate: [3]f64 = undefined;
    var half_saturation: [3]f64 = undefined;
    var inhibition: [3]f64 = undefined;
    try fillParameterArrays(parameters, &populations, .{
        .enabled = &enabled,
        .reference_energy_yield_kj_per_g_c = &reference,
        .growth_energy_requirement_kj_per_g_c = &growth,
        .minimum_respiration_requirement_g_c_per_g_c = &minimum,
        .specific_oxidation_rate_g_c_per_g_c_h = &rate,
        .dissolved_organic_carbon_half_saturation_g_c_per_m3 = &half_saturation,
        .acetate_product_inhibition_g_c_per_m3 = &inhibition,
    });
    try std.testing.expectEqual([3]bool{ true, true, false }, enabled);
    try std.testing.expectEqual([3]f64{ 0.4, 0.5, 0.4 }, minimum);

    var state = try State.init(std.testing.allocator, populations.len);
    defer state.deinit();
    const ones = [_]f64{ 1, 1, 1 };
    const temperatures = [_]f64{ 300, 300, 300 };
    const acetate = [_]f64{ 12, 12, 12 };
    const doc_concentration = [_]f64{ 12, 12, 12 };
    const doc_mass = [_]f64{ 10, 10, 10 };
    const zeros = [_]f64{ 0, 0, 0 };
    try calculate(&state, .{
        .enabled = &enabled,
        .aqueous_acetate_concentration_g_c_per_m3 = &acetate,
        .aqueous_dissolved_organic_carbon_concentration_g_c_per_m3 = &doc_concentration,
        .aqueous_dissolved_organic_carbon_g_c = &doc_mass,
        .dissolved_organic_carbon_competition_fraction = &ones,
        .combined_nutrient_limitation_fraction = &ones,
        .water_response = &ones,
        .active_biomass_g_c = &ones,
        .growth_temperature_response = &ones,
        .fermentation_oxygen_inhibition_fraction = &ones,
        .soil_temperature_k = &temperatures,
        .hydrogen_feedback_energy_kj_per_mol = &zeros,
        .reference_fermentation_energy_yield_kj_per_g_c = &reference,
        .growth_energy_requirement_kj_per_g_c = &growth,
        .minimum_respiration_requirement_g_c_per_g_c = &minimum,
        .specific_oxidation_rate_g_c_per_g_c_h = &rate,
        .dissolved_organic_carbon_half_saturation_g_c_per_m3 = &half_saturation,
        .acetate_product_inhibition_g_c_per_m3 = &inhibition,
        .minimum_acetate_concentration_g_c_per_m3 = 1e-12,
        .gas_constant_kj_per_mol_k = parameters.gas_constant_kj_per_mol_k,
        .acetate_feedback_stoichiometric_exponent = parameters.acetate_feedback_stoichiometric_exponent,
        .feedback_carbon_conversion_g_c_per_mol = parameters.feedback_carbon_conversion_g_c_per_mol,
        .timestep_h = 1,
    });
    try std.testing.expect(state.actual_respiration_g_c[0] > 0);
    try std.testing.expect(state.actual_respiration_g_c[1] > 0);
    try std.testing.expectEqual(@as(f64, 0), state.actual_respiration_g_c[2]);
}

test "product energy feedback and DOC caps reproduce NITRO fermentation" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, .{
        .enabled = &.{true},
        .aqueous_acetate_concentration_g_c_per_m3 = &.{2},
        .aqueous_dissolved_organic_carbon_concentration_g_c_per_m3 = &.{3},
        .aqueous_dissolved_organic_carbon_g_c = &.{4},
        .dissolved_organic_carbon_competition_fraction = &.{0.5},
        .combined_nutrient_limitation_fraction = &.{0.8},
        .water_response = &.{0.5},
        .active_biomass_g_c = &.{10},
        .growth_temperature_response = &.{2},
        .fermentation_oxygen_inhibition_fraction = &.{0.25},
        .soil_temperature_k = &.{300},
        .hydrogen_feedback_energy_kj_per_mol = &.{72},
        .reference_fermentation_energy_yield_kj_per_g_c = &.{10},
        .growth_energy_requirement_kj_per_g_c = &.{4},
        .minimum_respiration_requirement_g_c_per_g_c = &.{0.2},
        .specific_oxidation_rate_g_c_per_g_c_h = &.{0.2},
        .dissolved_organic_carbon_half_saturation_g_c_per_m3 = &.{1},
        .acetate_product_inhibition_g_c_per_m3 = &.{1},
        .minimum_acetate_concentration_g_c_per_m3 = 1e-12,
        .gas_constant_kj_per_mol_k = 8.3143e-3,
        .acetate_feedback_stoichiometric_exponent = 2,
        .feedback_carbon_conversion_g_c_per_mol = 72,
        .timestep_h = 1,
    });
    const acetate_feedback_mol =
        8.3143e-3 * 300 * @log(@as(f64, 4));
    const combined_feedback = 1 + acetate_feedback_mol / 72;
    const requirement = @max(
        @as(f64, 0.2),
        @min(@as(f64, 1), 1 / (1 + @max(0, 10 - combined_feedback) / 4)),
    );
    try std.testing.expectApproxEqAbs(acetate_feedback_mol, state.acetate_feedback_energy_kj_per_mol[0], 1e-15);
    try std.testing.expectApproxEqAbs(requirement, state.respiration_requirement_g_c_per_g_c[0], 1e-15);
    // Unlimited=.8; substrate=(3/4)*.25; microbial=.8*.1875*2=.3.
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), state.substrate_unlimited_respiration_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), state.microbial_limited_respiration_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@min(2 * requirement, 0.3), state.actual_respiration_g_c[0], 1e-15);
    try std.testing.expectEqual(state.actual_respiration_g_c[0], state.hydrogen_production_equivalent_g_c[0]);
}

test "invalid late fermenter leaves prior result unchanged" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.actual_respiration_g_c[0] = 7;
    const zeros = [_]f64{ 0, 0 };
    const ones = [_]f64{ 1, 1 };
    try std.testing.expectError(error.InvalidFermenterRespirationInput, calculate(&state, .{
        .enabled = &.{ true, true },
        .aqueous_acetate_concentration_g_c_per_m3 = &zeros,
        .aqueous_dissolved_organic_carbon_concentration_g_c_per_m3 = &zeros,
        .aqueous_dissolved_organic_carbon_g_c = &zeros,
        .dissolved_organic_carbon_competition_fraction = &ones,
        .combined_nutrient_limitation_fraction = &ones,
        .water_response = &ones,
        .active_biomass_g_c = &ones,
        .growth_temperature_response = &ones,
        .fermentation_oxygen_inhibition_fraction = &ones,
        .soil_temperature_k = &.{ 300, -1 },
        .hydrogen_feedback_energy_kj_per_mol = &zeros,
        .reference_fermentation_energy_yield_kj_per_g_c = &ones,
        .growth_energy_requirement_kj_per_g_c = &ones,
        .minimum_respiration_requirement_g_c_per_g_c = &ones,
        .specific_oxidation_rate_g_c_per_g_c_h = &ones,
        .dissolved_organic_carbon_half_saturation_g_c_per_m3 = &ones,
        .acetate_product_inhibition_g_c_per_m3 = &ones,
        .minimum_acetate_concentration_g_c_per_m3 = 1e-12,
        .gas_constant_kj_per_mol_k = 8.3143e-3,
        .acetate_feedback_stoichiometric_exponent = 2,
        .feedback_carbon_conversion_g_c_per_mol = 72,
        .timestep_h = 1,
    }));
    try std.testing.expectEqual(@as(f64, 7), state.actual_respiration_g_c[0]);
}

test "NITRO 900-973 derived overflow leaves fermenter state unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.actual_respiration_g_c[0] = 7;
    const one = [_]f64{1};
    try std.testing.expectError(
        error.NonFiniteFermenterRespirationResult,
        calculate(&state, .{
            .enabled = &.{true},
            .aqueous_acetate_concentration_g_c_per_m3 = &.{2},
            .aqueous_dissolved_organic_carbon_concentration_g_c_per_m3 = &one,
            .aqueous_dissolved_organic_carbon_g_c = &one,
            .dissolved_organic_carbon_competition_fraction = &one,
            .combined_nutrient_limitation_fraction = &one,
            .water_response = &one,
            .active_biomass_g_c = &one,
            .growth_temperature_response = &one,
            .fermentation_oxygen_inhibition_fraction = &one,
            .soil_temperature_k = &.{300},
            .hydrogen_feedback_energy_kj_per_mol = &one,
            .reference_fermentation_energy_yield_kj_per_g_c = &one,
            .growth_energy_requirement_kj_per_g_c = &one,
            .minimum_respiration_requirement_g_c_per_g_c = &one,
            .specific_oxidation_rate_g_c_per_g_c_h = &one,
            .dissolved_organic_carbon_half_saturation_g_c_per_m3 = &one,
            .acetate_product_inhibition_g_c_per_m3 = &one,
            .minimum_acetate_concentration_g_c_per_m3 = 1e-12,
            .gas_constant_kj_per_mol_k = std.math.floatMax(f64),
            .acetate_feedback_stoichiometric_exponent = 2,
            .feedback_carbon_conversion_g_c_per_mol = 72,
            .timestep_h = 1,
        }),
    );
    try std.testing.expectEqual(@as(f64, 7), state.actual_respiration_g_c[0]);
}
