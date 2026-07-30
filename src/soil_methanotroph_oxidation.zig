const std = @import("std");
const methane_oxidation = @import("soil_methane_oxidation.zig");

pub const Inputs = struct {
    enabled: []const bool,
    gaseous_methane_g_c: []const f64,
    aqueous_methane_g_c: []const f64,
    gaseous_methane_flux_g_c: []const f64,
    aqueous_methane_flux_g_c: []const f64,
    heterotrophic_methanogenesis_g_c: []const f64,
    autotrophic_methanogenesis_g_c: []const f64,
    water_volume_m3: []const f64,
    air_volume_m3: []const f64,
    methane_solubility_water_to_air: []const f64,
    gas_exchange_rate_per_h: []const f64,
    gas_exchange_enabled: []const bool,
    minimum_gaseous_methane_g_c: []const f64,
    methane_half_saturation_g_c_per_m3: []const f64,
    specific_methane_oxidation_g_c_per_g_c_h: []const f64,
    temperature_water_response: []const f64,
    combined_nutrient_limitation_fraction: []const f64,
    active_biomass_g_c: []const f64,
    biomass_conversion_efficiency_g_c_per_g_c: []const f64,
    growth_respiration_g_c_per_g_c: []const f64,
    oxygen_per_respired_carbon_g_o_per_g_c: f64,
    oxygen_per_oxidized_methane_carbon_g_o_per_g_c: f64,
    timestep_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    methane_unlimited_oxidation_g_c: []f64,
    final_gaseous_methane_g_c: []f64,
    final_aqueous_methane_g_c: []f64,
    methane_oxidation_g_c: []f64,
    growth_respiration_g_c: []f64,
    gas_to_water_exchange_g_c: []f64,
    respiration_oxygen_demand_g_o: []f64,
    total_oxygen_demand_g_o: []f64,
    iterations: []u16,
    newton_raphson_steps: []u16,
    picard_steps: []u16,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidMethanotrophDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
        var allocated_f64: usize = 0;
        var allocated_u16: usize = 0;
        errdefer {
            inline for (@typeInfo(State).@"struct".fields) |field| {
                if (field.type == []f64 and allocated_f64 > 0) {
                    allocated_f64 -= 1;
                    allocator.free(@field(state, field.name));
                } else if (field.type == []u16 and allocated_u16 > 0) {
                    allocated_u16 -= 1;
                    allocator.free(@field(state, field.name));
                }
            }
        }
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64) {
                @field(state, field.name) = try allocator.alloc(f64, unit_count);
                @memset(@field(state, field.name), 0);
                allocated_f64 += 1;
            } else if (field.type == []u16) {
                @field(state, field.name) = try allocator.alloc(u16, unit_count);
                @memset(@field(state, field.name), 0);
                allocated_u16 += 1;
            }
        }
        return state;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64 or field.type == []u16)
                self.allocator.free(@field(self, field.name));
        }
        self.* = undefined;
    }
};

/// NITRO.F 1336--1444. The source NPH×NPT repetition is replaced by the
/// equivalent conservative whole-step methane residual solved with bounded
/// Newton/Picard convergence. `options.gas_max_iterations` is the runtime
/// NPH×NPG ceiling and convergence exits early.
pub fn calculate(
    state: *State,
    inputs: Inputs,
    options: methane_oxidation.Options,
) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.unit_count);
    defer staged.deinit();
    const temporary = try state.allocator.alloc(
        methane_oxidation.Result,
        state.unit_count,
    );
    defer state.allocator.free(temporary);
    const unlimited = try state.allocator.alloc(f64, state.unit_count);
    defer state.allocator.free(unlimited);
    for (0..state.unit_count) |unit| {
        if (!inputs.enabled[unit]) {
            temporary[unit] = zeroResult(
                inputs.gaseous_methane_g_c[unit],
                inputs.aqueous_methane_g_c[unit],
            );
            unlimited[unit] = 0;
            continue;
        }
        unlimited[unit] =
            inputs.specific_methane_oxidation_g_c_per_g_c_h[unit] *
            inputs.temperature_water_response[unit] *
            inputs.combined_nutrient_limitation_fraction[unit] *
            inputs.active_biomass_g_c[unit] *
            inputs.timestep_h;
        temporary[unit] = try methane_oxidation.solve(.{
            .gaseous_methane_g_c = inputs.gaseous_methane_g_c[unit],
            .aqueous_methane_g_c = inputs.aqueous_methane_g_c[unit],
            .gaseous_methane_flux_g_c = inputs.gaseous_methane_flux_g_c[unit],
            .aqueous_methane_flux_g_c = inputs.aqueous_methane_flux_g_c[unit],
            .methanogenesis_g_c = inputs.heterotrophic_methanogenesis_g_c[unit] +
                inputs.autotrophic_methanogenesis_g_c[unit],
            .water_volume_m3 = inputs.water_volume_m3[unit],
            .air_volume_m3 = inputs.air_volume_m3[unit],
            .methane_solubility_water_to_air = inputs.methane_solubility_water_to_air[unit],
            .gas_exchange_rate_per_step = inputs.gas_exchange_rate_per_h[unit] * inputs.timestep_h,
            .gas_exchange_enabled = inputs.gas_exchange_enabled[unit],
            .minimum_gaseous_methane_g_c = inputs.minimum_gaseous_methane_g_c[unit],
            .methane_half_saturation_g_c_per_m3 = inputs.methane_half_saturation_g_c_per_m3[unit],
            .maximum_methane_oxidation_g_c = unlimited[unit],
            .biomass_conversion_efficiency_g_c_per_g_c = inputs.biomass_conversion_efficiency_g_c_per_g_c[unit],
            .growth_respiration_g_c_per_g_c = inputs.growth_respiration_g_c_per_g_c[unit],
        }, options);
    }

    // No live output changes until every unit has converged.
    for (temporary, 0..) |result, unit| {
        staged.methane_unlimited_oxidation_g_c[unit] = unlimited[unit];
        staged.final_gaseous_methane_g_c[unit] =
            result.gaseous_methane_g_c;
        staged.final_aqueous_methane_g_c[unit] =
            result.aqueous_methane_g_c;
        staged.methane_oxidation_g_c[unit] = result.methane_oxidation_g_c;
        staged.growth_respiration_g_c[unit] =
            result.growth_respiration_g_c;
        staged.gas_to_water_exchange_g_c[unit] =
            result.gas_to_water_exchange_g_c;
        staged.respiration_oxygen_demand_g_o[unit] =
            inputs.oxygen_per_respired_carbon_g_o_per_g_c *
            result.growth_respiration_g_c;
        staged.total_oxygen_demand_g_o[unit] =
            staged.respiration_oxygen_demand_g_o[unit] +
            inputs.oxygen_per_oxidized_methane_carbon_g_o_per_g_c *
                result.methane_oxidation_g_c;
        staged.iterations[unit] = result.iterations;
        staged.newton_raphson_steps[unit] = result.newton_raphson_steps;
        staged.picard_steps[unit] = result.picard_steps;
    }
    try validateResult(&staged);
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64 or field.type == []u16)
            @memcpy(@field(state, field.name), @field(staged, field.name));
}

pub fn sourceEnabled(
    zero_based_complex: usize,
    zero_based_population: usize,
) bool {
    return zero_based_complex == 5 and zero_based_population == 2;
}

fn zeroResult(gas_g_c: f64, aqueous_g_c: f64) methane_oxidation.Result {
    return .{
        .gaseous_methane_g_c = gas_g_c,
        .aqueous_methane_g_c = aqueous_g_c,
        .methane_oxidation_g_c = 0,
        .growth_respiration_g_c = 0,
        .gas_to_water_exchange_g_c = 0,
        .oxygen_demand_g_o = 0,
        .iterations = 0,
        .newton_raphson_steps = 0,
        .picard_steps = 0,
        .residual_g_c = 0,
    };
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.enabled.len != n or
        inputs.gas_exchange_enabled.len != n)
        return error.InvalidMethanotrophDimensions;
    inline for (.{
        inputs.gaseous_methane_g_c,
        inputs.aqueous_methane_g_c,
        inputs.gaseous_methane_flux_g_c,
        inputs.aqueous_methane_flux_g_c,
        inputs.heterotrophic_methanogenesis_g_c,
        inputs.autotrophic_methanogenesis_g_c,
        inputs.water_volume_m3,
        inputs.air_volume_m3,
        inputs.methane_solubility_water_to_air,
        inputs.gas_exchange_rate_per_h,
        inputs.minimum_gaseous_methane_g_c,
        inputs.methane_half_saturation_g_c_per_m3,
        inputs.specific_methane_oxidation_g_c_per_g_c_h,
        inputs.temperature_water_response,
        inputs.combined_nutrient_limitation_fraction,
        inputs.active_biomass_g_c,
        inputs.biomass_conversion_efficiency_g_c_per_g_c,
        inputs.growth_respiration_g_c_per_g_c,
    }) |values| {
        if (values.len != n) return error.InvalidMethanotrophDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidMethanotrophInput;
    }
    for (0..n) |unit| if (inputs.enabled[unit]) {
        inline for (.{
            inputs.gaseous_methane_g_c[unit],
            inputs.aqueous_methane_g_c[unit],
            inputs.heterotrophic_methanogenesis_g_c[unit],
            inputs.autotrophic_methanogenesis_g_c[unit],
            inputs.air_volume_m3[unit],
            inputs.methane_solubility_water_to_air[unit],
            inputs.gas_exchange_rate_per_h[unit],
            inputs.minimum_gaseous_methane_g_c[unit],
            inputs.specific_methane_oxidation_g_c_per_g_c_h[unit],
            inputs.temperature_water_response[unit],
            inputs.combined_nutrient_limitation_fraction[unit],
            inputs.active_biomass_g_c[unit],
            inputs.biomass_conversion_efficiency_g_c_per_g_c[unit],
            inputs.growth_respiration_g_c_per_g_c[unit],
        }) |value| if (value < 0) return error.InvalidMethanotrophInput;
        if (inputs.water_volume_m3[unit] <= 0 or
            inputs.methane_half_saturation_g_c_per_m3[unit] <= 0 or
            inputs.combined_nutrient_limitation_fraction[unit] > 1)
            return error.InvalidMethanotrophInput;
    };
    inline for (.{
        inputs.oxygen_per_respired_carbon_g_o_per_g_c,
        inputs.oxygen_per_oxidized_methane_carbon_g_o_per_g_c,
        inputs.timestep_h,
    }) |value| if (!std.math.isFinite(value) or value <= 0)
        return error.InvalidMethanotrophInput;
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state.*, field.name)) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteMethanotrophResult;
            if (!std.mem.eql(u8, field.name, "gas_to_water_exchange_g_c") and
                value < 0)
                return error.NonFiniteMethanotrophResult;
        };
}

fn testOptions() methane_oxidation.Options {
    return .{
        .absolute_tolerance_g_c = 1e-12,
        .relative_tolerance = 1e-10,
        .derivative_floor = 1e-14,
        .picard_relaxation = 0.5,
        .gas_max_iterations = 80,
    };
}

test "source methanotroph role reproduces NITRO selector" {
    try std.testing.expect(sourceEnabled(5, 2));
    try std.testing.expect(!sourceEnabled(5, 1));
    try std.testing.expect(!sourceEnabled(4, 2));
}

test "runtime methanotroph solve converges early and conserves carbon" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, .{
        .enabled = &.{ true, true },
        .gaseous_methane_g_c = &.{ 4, 2 },
        .aqueous_methane_g_c = &.{ 1, 0.5 },
        .gaseous_methane_flux_g_c = &.{ 0.1, 0 },
        .aqueous_methane_flux_g_c = &.{ 0.1, 0 },
        .heterotrophic_methanogenesis_g_c = &.{ 0.1, 0 },
        .autotrophic_methanogenesis_g_c = &.{ 0.1, 0 },
        .water_volume_m3 = &.{ 2, 1 },
        .air_volume_m3 = &.{ 3, 1 },
        .methane_solubility_water_to_air = &.{ 0.03, 0.03 },
        .gas_exchange_rate_per_h = &.{ 0.5, 0.5 },
        .gas_exchange_enabled = &.{ true, true },
        .minimum_gaseous_methane_g_c = &.{ 1e-12, 1e-12 },
        .methane_half_saturation_g_c_per_m3 = &.{ 0.2, 0.2 },
        .specific_methane_oxidation_g_c_per_g_c_h = &.{ 0.1, 0.1 },
        .temperature_water_response = &.{ 0.5, 0.5 },
        .combined_nutrient_limitation_fraction = &.{ 0.8, 0.8 },
        .active_biomass_g_c = &.{ 10, 5 },
        .biomass_conversion_efficiency_g_c_per_g_c = &.{ 0.4, 0.4 },
        .growth_respiration_g_c_per_g_c = &.{ 0.5, 0.5 },
        .oxygen_per_respired_carbon_g_o_per_g_c = 2.667,
        .oxygen_per_oxidized_methane_carbon_g_o_per_g_c = 5.333,
        .timestep_h = 1,
    }, testOptions());
    const before: f64 = 4 + 1 + 0.1 + 0.1 + 0.1 + 0.1;
    try std.testing.expect(state.iterations[0] < 80);
    try std.testing.expectApproxEqAbs(
        before,
        state.final_gaseous_methane_g_c[0] +
            state.final_aqueous_methane_g_c[0] +
            state.methane_oxidation_g_c[0] +
            state.growth_respiration_g_c[0],
        1e-9,
    );
    try std.testing.expectApproxEqAbs(
        2.667 * state.growth_respiration_g_c[0] +
            5.333 * state.methane_oxidation_g_c[0],
        state.total_oxygen_demand_g_o[0],
        1e-14,
    );
}

test "late invalid methanotroph leaves every output unchanged" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.methane_oxidation_g_c[0] = 7;
    const zeros = [_]f64{ 0, 0 };
    const ones = [_]f64{ 1, 1 };
    try std.testing.expectError(error.InvalidMethanotrophInput, calculate(&state, .{
        .enabled = &.{ true, true },
        .gaseous_methane_g_c = &ones,
        .aqueous_methane_g_c = &ones,
        .gaseous_methane_flux_g_c = &zeros,
        .aqueous_methane_flux_g_c = &zeros,
        .heterotrophic_methanogenesis_g_c = &zeros,
        .autotrophic_methanogenesis_g_c = &zeros,
        .water_volume_m3 = &.{ 1, -1 },
        .air_volume_m3 = &ones,
        .methane_solubility_water_to_air = &ones,
        .gas_exchange_rate_per_h = &ones,
        .gas_exchange_enabled = &.{ true, true },
        .minimum_gaseous_methane_g_c = &zeros,
        .methane_half_saturation_g_c_per_m3 = &ones,
        .specific_methane_oxidation_g_c_per_g_c_h = &ones,
        .temperature_water_response = &ones,
        .combined_nutrient_limitation_fraction = &ones,
        .active_biomass_g_c = &ones,
        .biomass_conversion_efficiency_g_c_per_g_c = &ones,
        .growth_respiration_g_c_per_g_c = &ones,
        .oxygen_per_respired_carbon_g_o_per_g_c = 2.667,
        .oxygen_per_oxidized_methane_carbon_g_o_per_g_c = 5.333,
        .timestep_h = 1,
    }, testOptions()));
    try std.testing.expectEqual(@as(f64, 7), state.methane_oxidation_g_c[0]);
}

test "NITRO 1336-1444 diagnostic overflow preserves methanotroph state" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.methane_oxidation_g_c[0] = 7;
    try std.testing.expectError(
        error.NonFiniteMethanotrophResult,
        calculate(&state, .{
            .enabled = &.{true},
            .gaseous_methane_g_c = &.{1},
            .aqueous_methane_g_c = &.{4},
            .gaseous_methane_flux_g_c = &.{0},
            .aqueous_methane_flux_g_c = &.{0},
            .heterotrophic_methanogenesis_g_c = &.{0},
            .autotrophic_methanogenesis_g_c = &.{0},
            .water_volume_m3 = &.{1},
            .air_volume_m3 = &.{1},
            .methane_solubility_water_to_air = &.{0},
            .gas_exchange_rate_per_h = &.{0},
            .gas_exchange_enabled = &.{false},
            .minimum_gaseous_methane_g_c = &.{0},
            .methane_half_saturation_g_c_per_m3 = &.{1},
            .specific_methane_oxidation_g_c_per_g_c_h = &.{1},
            .temperature_water_response = &.{1},
            .combined_nutrient_limitation_fraction = &.{1},
            .active_biomass_g_c = &.{4},
            .biomass_conversion_efficiency_g_c_per_g_c = &.{1},
            .growth_respiration_g_c_per_g_c = &.{1},
            .oxygen_per_respired_carbon_g_o_per_g_c = std.math.floatMax(f64),
            .oxygen_per_oxidized_methane_carbon_g_o_per_g_c = std.math.floatMax(f64),
            .timestep_h = 1,
        }, testOptions()),
    );
    try std.testing.expectEqual(@as(f64, 7), state.methane_oxidation_g_c[0]);
}
