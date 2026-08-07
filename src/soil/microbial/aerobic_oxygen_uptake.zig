const std = @import("std");
const oxygen_solver = @import("../gas/oxygen_solver.zig");

pub const Inputs = struct {
    enabled: []const bool,
    oxygen_demand_g_o: []const f64,
    oxygen_competition_fraction: []const f64,
    gaseous_oxygen_g_o: []const f64,
    aqueous_oxygen_g_o: []const f64,
    gaseous_oxygen_flux_g_o: []const f64,
    aqueous_oxygen_flux_g_o: []const f64,
    water_volume_m3: []const f64,
    air_volume_m3: []const f64,
    oxygen_solubility_water_to_air: []const f64,
    gas_exchange_rate_per_h: []const f64,
    microbial_radius_m: []const f64,
    water_film_thickness_m: []const f64,
    tortuosity: []const f64,
    aqueous_oxygen_diffusivity_m2_per_h: []const f64,
    microbial_count_per_g_c: []const f64,
    active_biomass_g_c: []const f64,
    oxygen_half_saturation_g_o_per_m3: []const f64,
    maximum_aqueous_oxygen_concentration_g_o_per_m3: []const f64,
    timestep_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    uptake_conductance_m3: []f64,
    allocated_final_gaseous_oxygen_g_o: []f64,
    allocated_final_aqueous_oxygen_g_o: []f64,
    oxygen_uptake_g_o: []f64,
    gas_to_water_exchange_g_o: []f64,
    demand_satisfaction_fraction: []f64,
    iterations: []u16,
    newton_raphson_steps: []u16,
    picard_steps: []u16,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidAerobicOxygenDimensions;
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
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64 or field.type == []u16)
                self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

const Temporary = struct {
    conductance_m3: f64,
    result: oxygen_solver.Result,
};

/// NITRO.F 1458--1614. The nested NPH×NPT dissolution/uptake cycle is one
/// conservative whole-step Newton/Picard residual per runtime population.
/// The runtime NPH×NPG value is supplied as `options.gas_max_iterations`.
pub fn calculate(
    state: *State,
    inputs: Inputs,
    options: oxygen_solver.Options,
) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.unit_count);
    defer staged.deinit();
    const temporary = try state.allocator.alloc(Temporary, state.unit_count);
    defer state.allocator.free(temporary);
    for (0..state.unit_count) |unit| {
        if (!inputs.enabled[unit] or inputs.oxygen_demand_g_o[unit] == 0) {
            temporary[unit] = .{
                .conductance_m3 = 0,
                .result = zeroResult(
                    inputs.gaseous_oxygen_g_o[unit] *
                        inputs.oxygen_competition_fraction[unit],
                    inputs.aqueous_oxygen_g_o[unit] *
                        inputs.oxygen_competition_fraction[unit],
                ),
            };
            continue;
        }
        const conductance = try oxygen_solver.uptakeConductance_m3_per_step(.{
            .microbial_radius_m = inputs.microbial_radius_m[unit],
            .water_film_thickness_m = inputs.water_film_thickness_m[unit],
            .tortuosity = inputs.tortuosity[unit],
            .aqueous_oxygen_diffusivity_m2_per_step = inputs.aqueous_oxygen_diffusivity_m2_per_h[unit] *
                inputs.timestep_h,
            .microbial_count_per_g_c = inputs.microbial_count_per_g_c[unit],
            .active_biomass_g_c = inputs.active_biomass_g_c[unit],
        });
        const share = inputs.oxygen_competition_fraction[unit];
        temporary[unit] = .{
            .conductance_m3 = conductance,
            .result = try oxygen_solver.solve(.{
                .allocated_gaseous_oxygen_g_o = inputs.gaseous_oxygen_g_o[unit] * share,
                .allocated_aqueous_oxygen_g_o = inputs.aqueous_oxygen_g_o[unit] * share,
                .allocated_gaseous_flux_g_o = inputs.gaseous_oxygen_flux_g_o[unit] *
                    inputs.timestep_h * share,
                .allocated_aqueous_flux_g_o = inputs.aqueous_oxygen_flux_g_o[unit] *
                    inputs.timestep_h * share,
                .water_volume_m3 = inputs.water_volume_m3[unit],
                .air_volume_m3 = inputs.air_volume_m3[unit],
                .population_allocation_fraction = share,
                .oxygen_solubility_water_to_air = inputs.oxygen_solubility_water_to_air[unit],
                .gas_exchange_rate_per_step = inputs.gas_exchange_rate_per_h[unit] *
                    inputs.timestep_h,
                .uptake_conductance_m3_per_step = conductance,
                .oxygen_half_saturation_g_o_per_m3 = inputs.oxygen_half_saturation_g_o_per_m3[unit],
                .maximum_oxygen_uptake_g_o = inputs.oxygen_demand_g_o[unit] *
                    inputs.timestep_h,
                .maximum_aqueous_oxygen_concentration_g_o_per_m3 = inputs.maximum_aqueous_oxygen_concentration_g_o_per_m3[unit],
            }, options),
        };
    }
    // Publish only after all population solves converge.
    for (temporary, 0..) |item, unit| {
        staged.uptake_conductance_m3[unit] = item.conductance_m3;
        staged.allocated_final_gaseous_oxygen_g_o[unit] =
            item.result.gaseous_oxygen_g_o;
        staged.allocated_final_aqueous_oxygen_g_o[unit] =
            item.result.aqueous_oxygen_g_o;
        staged.oxygen_uptake_g_o[unit] = item.result.oxygen_uptake_g_o;
        staged.gas_to_water_exchange_g_o[unit] =
            item.result.gas_to_water_exchange_g_o;
        staged.demand_satisfaction_fraction[unit] =
            if (!inputs.enabled[unit] or
            inputs.oxygen_demand_g_o[unit] == 0)
                1
            else
                item.result.demand_satisfaction_fraction;
        staged.iterations[unit] = item.result.iterations;
        staged.newton_raphson_steps[unit] =
            item.result.newton_raphson_steps;
        staged.picard_steps[unit] = item.result.picard_steps;
    }
    try validateResult(&staged);
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64 or field.type == []u16)
            @memcpy(@field(state, field.name), @field(staged, field.name));
}

pub fn sourceEnabled(zero_based_population: usize) bool {
    return zero_based_population <= 2 or zero_based_population == 5;
}

fn zeroResult(gas: f64, aqueous: f64) oxygen_solver.Result {
    return .{
        .gaseous_oxygen_g_o = gas,
        .aqueous_oxygen_g_o = aqueous,
        .oxygen_uptake_g_o = 0,
        .gas_to_water_exchange_g_o = 0,
        .demand_satisfaction_fraction = 1,
        .iterations = 0,
        .newton_raphson_steps = 0,
        .picard_steps = 0,
        .residual_g_o = 0,
    };
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.enabled.len != n)
        return error.InvalidAerobicOxygenDimensions;
    inline for (.{
        inputs.oxygen_demand_g_o,
        inputs.oxygen_competition_fraction,
        inputs.gaseous_oxygen_g_o,
        inputs.aqueous_oxygen_g_o,
        inputs.gaseous_oxygen_flux_g_o,
        inputs.aqueous_oxygen_flux_g_o,
        inputs.water_volume_m3,
        inputs.air_volume_m3,
        inputs.oxygen_solubility_water_to_air,
        inputs.gas_exchange_rate_per_h,
        inputs.microbial_radius_m,
        inputs.water_film_thickness_m,
        inputs.tortuosity,
        inputs.aqueous_oxygen_diffusivity_m2_per_h,
        inputs.microbial_count_per_g_c,
        inputs.active_biomass_g_c,
        inputs.oxygen_half_saturation_g_o_per_m3,
        inputs.maximum_aqueous_oxygen_concentration_g_o_per_m3,
    }) |values| {
        if (values.len != n) return error.InvalidAerobicOxygenDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidAerobicOxygenInput;
    }
    for (0..n) |unit| if (inputs.enabled[unit]) {
        inline for (.{
            inputs.oxygen_demand_g_o[unit],
            inputs.oxygen_competition_fraction[unit],
            inputs.gaseous_oxygen_g_o[unit],
            inputs.aqueous_oxygen_g_o[unit],
            inputs.air_volume_m3[unit],
            inputs.oxygen_solubility_water_to_air[unit],
            inputs.gas_exchange_rate_per_h[unit],
            inputs.tortuosity[unit],
            inputs.aqueous_oxygen_diffusivity_m2_per_h[unit],
            inputs.microbial_count_per_g_c[unit],
            inputs.active_biomass_g_c[unit],
        }) |value| if (value < 0) return error.InvalidAerobicOxygenInput;
        if (inputs.oxygen_competition_fraction[unit] <= 0 or
            inputs.oxygen_competition_fraction[unit] > 1 or
            inputs.water_volume_m3[unit] <= 0 or
            inputs.microbial_radius_m[unit] <= 0 or
            inputs.water_film_thickness_m[unit] <= 0 or
            inputs.oxygen_half_saturation_g_o_per_m3[unit] <= 0 or
            inputs.maximum_aqueous_oxygen_concentration_g_o_per_m3[unit] <= 0)
            return error.InvalidAerobicOxygenInput;
    };
    if (!std.math.isFinite(inputs.timestep_h) or inputs.timestep_h <= 0)
        return error.InvalidAerobicOxygenInput;
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state.*, field.name)) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteAerobicOxygenResult;
            if (!std.mem.eql(u8, field.name, "gas_to_water_exchange_g_o") and
                value < 0)
                return error.NonFiniteAerobicOxygenResult;
        };
    for (state.demand_satisfaction_fraction) |fraction|
        if (fraction > 1) return error.InvalidAerobicOxygenSatisfaction;
}

fn testOptions() oxygen_solver.Options {
    return .{
        .absolute_tolerance_g_o = 1e-12,
        .relative_tolerance = 1e-10,
        .derivative_floor = 1e-14,
        .picard_relaxation = 0.5,
        .gas_max_iterations = 80,
    };
}

test "source aerobic roles reproduce NITRO selector" {
    try std.testing.expect(sourceEnabled(0));
    try std.testing.expect(sourceEnabled(2));
    try std.testing.expect(sourceEnabled(5));
    try std.testing.expect(!sourceEnabled(3));
}

test "runtime population oxygen uptake converges and conserves allocation" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, .{
        .enabled = &.{ true, true },
        .oxygen_demand_g_o = &.{ 0.8, 0.4 },
        .oxygen_competition_fraction = &.{ 0.6, 0.4 },
        .gaseous_oxygen_g_o = &.{ 5, 5 },
        .aqueous_oxygen_g_o = &.{ 1, 1 },
        .gaseous_oxygen_flux_g_o = &.{ 0, 0 },
        .aqueous_oxygen_flux_g_o = &.{ 0, 0 },
        .water_volume_m3 = &.{ 2, 2 },
        .air_volume_m3 = &.{ 3, 3 },
        .oxygen_solubility_water_to_air = &.{ 0.03, 0.03 },
        .gas_exchange_rate_per_h = &.{ 0.5, 0.5 },
        .microbial_radius_m = &.{ 1e-6, 1e-6 },
        .water_film_thickness_m = &.{ 2e-6, 2e-6 },
        .tortuosity = &.{ 0.5, 0.5 },
        .aqueous_oxygen_diffusivity_m2_per_h = &.{ 1e-4, 1e-4 },
        .microbial_count_per_g_c = &.{ 1e12, 1e12 },
        .active_biomass_g_c = &.{ 0.01, 0.01 },
        .oxygen_half_saturation_g_o_per_m3 = &.{ 0.1, 0.1 },
        .maximum_aqueous_oxygen_concentration_g_o_per_m3 = &.{ 1, 1 },
        .timestep_h = 1,
    }, testOptions());
    for (0..2) |unit| {
        const share: f64 = if (unit == 0) 0.6 else 0.4;
        try std.testing.expect(state.iterations[unit] < 80);
        try std.testing.expectApproxEqAbs(
            6 * share,
            state.allocated_final_gaseous_oxygen_g_o[unit] +
                state.allocated_final_aqueous_oxygen_g_o[unit] +
                state.oxygen_uptake_g_o[unit],
            1e-9,
        );
        try std.testing.expect(state.demand_satisfaction_fraction[unit] >= 0 and
            state.demand_satisfaction_fraction[unit] <= 1);
    }
}

test "late invalid oxygen population preserves all outputs" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.oxygen_uptake_g_o[0] = 7;
    const zeros = [_]f64{ 0, 0 };
    const ones = [_]f64{ 1, 1 };
    try std.testing.expectError(error.InvalidAerobicOxygenInput, calculate(&state, .{
        .enabled = &.{ true, true },
        .oxygen_demand_g_o = &ones,
        .oxygen_competition_fraction = &.{ 0.5, 2 },
        .gaseous_oxygen_g_o = &ones,
        .aqueous_oxygen_g_o = &ones,
        .gaseous_oxygen_flux_g_o = &zeros,
        .aqueous_oxygen_flux_g_o = &zeros,
        .water_volume_m3 = &ones,
        .air_volume_m3 = &ones,
        .oxygen_solubility_water_to_air = &ones,
        .gas_exchange_rate_per_h = &ones,
        .microbial_radius_m = &ones,
        .water_film_thickness_m = &ones,
        .tortuosity = &ones,
        .aqueous_oxygen_diffusivity_m2_per_h = &ones,
        .microbial_count_per_g_c = &ones,
        .active_biomass_g_c = &ones,
        .oxygen_half_saturation_g_o_per_m3 = &ones,
        .maximum_aqueous_oxygen_concentration_g_o_per_m3 = &ones,
        .timestep_h = 1,
    }, testOptions()));
    try std.testing.expectEqual(@as(f64, 7), state.oxygen_uptake_g_o[0]);
}
