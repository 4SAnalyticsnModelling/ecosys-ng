const std = @import("std");

pub const Metabolism = enum {
    aerobic,
    fermenter,
    methanogen,
};

pub const Inputs = struct {
    metabolism: []const Metabolism,
    autotrophic_complex: []const bool,
    oxygen_unlimited_respiration_g_c: []const f64,
    oxygen_demand_satisfaction_fraction: []const f64,
    respiration_oxygen_demand_g_o: []const f64,
    non_band_oxidation_g: []const f64,
    band_oxidation_g: []const f64,
    fermenter_carbon_dioxide_yield_g_c_per_g_c: f64,
    fermenter_acetate_yield_g_c_per_g_c: f64,
    fermenter_hydrogen_yield_g_h_per_g_c: f64,
    acetotrophic_carbon_dioxide_yield_g_c_per_g_c: f64,
    acetotrophic_methane_yield_g_c_per_g_c: f64,
    hydrogenotrophic_methane_yield_g_c_per_g_c: f64,
    hydrogenotrophic_hydrogen_uptake_g_h_per_g_c: f64,
    carbon_balance_tolerance_g_c: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    oxygen_limited_respiration_g_c: []f64,
    carbon_dioxide_production_g_c: []f64,
    acetate_production_g_c: []f64,
    methane_production_g_c: []f64,
    oxygen_limited_respiration_oxygen_g_o: []f64,
    hydrogen_production_g_h: []f64,
    hydrogen_uptake_g_h: []f64,
    oxygen_limited_non_band_oxidation_g: []f64,
    oxygen_limited_band_oxidation_g: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0)
            return error.InvalidRespirationProductDimensions;
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

/// NITRO.F 1616--1665. Runtime roles replace fixed N indices while retaining
/// the distinct heterotrophic and autotrophic N=5 product branches.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc(
        [9]f64,
        state.unit_count,
    );
    defer state.allocator.free(temporary);
    for (0..state.unit_count) |unit| {
        const unlimited =
            inputs.oxygen_unlimited_respiration_g_c[unit];
        var actual = unlimited;
        var carbon_dioxide: f64 = 0;
        var acetate: f64 = 0;
        var methane: f64 = 0;
        var respiration_oxygen =
            inputs.respiration_oxygen_demand_g_o[unit];
        var hydrogen_production: f64 = 0;
        var hydrogen_uptake: f64 = 0;
        var non_band_oxidation = inputs.non_band_oxidation_g[unit];
        var band_oxidation = inputs.band_oxidation_g[unit];
        switch (inputs.metabolism[unit]) {
            .aerobic => {
                const satisfaction =
                    inputs.oxygen_demand_satisfaction_fraction[unit];
                actual *= satisfaction;
                carbon_dioxide = actual;
                respiration_oxygen *= satisfaction;
                if (inputs.autotrophic_complex[unit]) {
                    non_band_oxidation *= satisfaction;
                    band_oxidation *= satisfaction;
                }
            },
            .fermenter => {
                carbon_dioxide =
                    inputs.fermenter_carbon_dioxide_yield_g_c_per_g_c *
                    actual;
                acetate =
                    inputs.fermenter_acetate_yield_g_c_per_g_c * actual;
                if (!inputs.autotrophic_complex[unit])
                    hydrogen_production =
                        inputs.fermenter_hydrogen_yield_g_h_per_g_c *
                        actual;
            },
            .methanogen => {
                if (inputs.autotrophic_complex[unit]) {
                    methane =
                        inputs.hydrogenotrophic_methane_yield_g_c_per_g_c *
                        actual;
                    hydrogen_uptake =
                        inputs.hydrogenotrophic_hydrogen_uptake_g_h_per_g_c *
                        actual;
                } else {
                    carbon_dioxide =
                        inputs.acetotrophic_carbon_dioxide_yield_g_c_per_g_c *
                        actual;
                    methane =
                        inputs.acetotrophic_methane_yield_g_c_per_g_c *
                        actual;
                }
            },
        }
        const carbon_products = carbon_dioxide + acetate + methane;
        const tolerance =
            inputs.carbon_balance_tolerance_g_c *
            @max(1, actual);
        inline for (.{
            actual,
            carbon_dioxide,
            acetate,
            methane,
            respiration_oxygen,
            hydrogen_production,
            hydrogen_uptake,
            non_band_oxidation,
            band_oxidation,
            carbon_products,
            tolerance,
        }) |value| if (!std.math.isFinite(value) or value < 0)
            return error.NonFiniteRespirationProductResult;
        if (@abs(carbon_products - actual) > tolerance)
            return error.RespirationProductCarbonBalanceFailure;
        temporary[unit] = .{
            actual,
            carbon_dioxide,
            acetate,
            methane,
            respiration_oxygen,
            hydrogen_production,
            hydrogen_uptake,
            non_band_oxidation,
            band_oxidation,
        };
    }
    for (temporary, 0..) |values, unit| {
        state.oxygen_limited_respiration_g_c[unit] = values[0];
        state.carbon_dioxide_production_g_c[unit] = values[1];
        state.acetate_production_g_c[unit] = values[2];
        state.methane_production_g_c[unit] = values[3];
        state.oxygen_limited_respiration_oxygen_g_o[unit] = values[4];
        state.hydrogen_production_g_h[unit] = values[5];
        state.hydrogen_uptake_g_h[unit] = values[6];
        state.oxygen_limited_non_band_oxidation_g[unit] = values[7];
        state.oxygen_limited_band_oxidation_g[unit] = values[8];
    }
}

pub fn sourceMetabolism(zero_based_population: usize) Metabolism {
    return switch (zero_based_population) {
        3, 6 => .fermenter,
        4 => .methanogen,
        else => .aerobic,
    };
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.metabolism.len != n or
        inputs.autotrophic_complex.len != n)
        return error.InvalidRespirationProductDimensions;
    inline for (.{
        inputs.oxygen_unlimited_respiration_g_c,
        inputs.oxygen_demand_satisfaction_fraction,
        inputs.respiration_oxygen_demand_g_o,
        inputs.non_band_oxidation_g,
        inputs.band_oxidation_g,
    }) |values| {
        if (values.len != n) return error.InvalidRespirationProductDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRespirationProductInput;
    }
    for (inputs.oxygen_demand_satisfaction_fraction) |fraction|
        if (fraction > 1) return error.InvalidRespirationProductInput;
    inline for (.{
        inputs.fermenter_carbon_dioxide_yield_g_c_per_g_c,
        inputs.fermenter_acetate_yield_g_c_per_g_c,
        inputs.fermenter_hydrogen_yield_g_h_per_g_c,
        inputs.acetotrophic_carbon_dioxide_yield_g_c_per_g_c,
        inputs.acetotrophic_methane_yield_g_c_per_g_c,
        inputs.hydrogenotrophic_methane_yield_g_c_per_g_c,
        inputs.hydrogenotrophic_hydrogen_uptake_g_h_per_g_c,
        inputs.carbon_balance_tolerance_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidRespirationProductInput;
}

test "source population roles reproduce NITRO product branches" {
    try std.testing.expectEqual(Metabolism.aerobic, sourceMetabolism(0));
    try std.testing.expectEqual(Metabolism.fermenter, sourceMetabolism(3));
    try std.testing.expectEqual(Metabolism.methanogen, sourceMetabolism(4));
    try std.testing.expectEqual(Metabolism.fermenter, sourceMetabolism(6));
}

test "oxygen scaling and all respiration products reproduce NITRO" {
    var state = try State.init(std.testing.allocator, 4);
    defer state.deinit();
    try calculate(&state, .{
        .metabolism = &.{ .aerobic, .fermenter, .methanogen, .methanogen },
        .autotrophic_complex = &.{ true, false, false, true },
        .oxygen_unlimited_respiration_g_c = &.{ 2, 3, 4, 5 },
        .oxygen_demand_satisfaction_fraction = &.{ 0.25, 0, 0, 0 },
        .respiration_oxygen_demand_g_o = &.{ 6, 0, 0, 0 },
        .non_band_oxidation_g = &.{ 8, 0, 0, 0 },
        .band_oxidation_g = &.{ 4, 0, 0, 0 },
        .fermenter_carbon_dioxide_yield_g_c_per_g_c = 0.333,
        .fermenter_acetate_yield_g_c_per_g_c = 0.667,
        .fermenter_hydrogen_yield_g_h_per_g_c = 0.111,
        .acetotrophic_carbon_dioxide_yield_g_c_per_g_c = 0.5,
        .acetotrophic_methane_yield_g_c_per_g_c = 0.5,
        .hydrogenotrophic_methane_yield_g_c_per_g_c = 1,
        .hydrogenotrophic_hydrogen_uptake_g_h_per_g_c = 0.667,
        .carbon_balance_tolerance_g_c = 1e-12,
    });
    try std.testing.expectEqual(@as(f64, 0.5), state.oxygen_limited_respiration_g_c[0]);
    try std.testing.expectEqual(@as(f64, 0.5), state.carbon_dioxide_production_g_c[0]);
    try std.testing.expectEqual(@as(f64, 1.5), state.oxygen_limited_respiration_oxygen_g_o[0]);
    try std.testing.expectEqual(@as(f64, 2), state.oxygen_limited_non_band_oxidation_g[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.999), state.carbon_dioxide_production_g_c[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.001), state.acetate_production_g_c[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.333), state.hydrogen_production_g_h[1], 1e-15);
    try std.testing.expectEqual(@as(f64, 2), state.carbon_dioxide_production_g_c[2]);
    try std.testing.expectEqual(@as(f64, 2), state.methane_production_g_c[2]);
    try std.testing.expectEqual(@as(f64, 5), state.methane_production_g_c[3]);
    try std.testing.expectApproxEqAbs(@as(f64, 3.335), state.hydrogen_uptake_g_h[3], 1e-15);
}

test "invalid product yields fail without replacing prior state" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.carbon_dioxide_production_g_c[0] = 7;
    try std.testing.expectError(error.RespirationProductCarbonBalanceFailure, calculate(&state, .{
        .metabolism = &.{.fermenter},
        .autotrophic_complex = &.{false},
        .oxygen_unlimited_respiration_g_c = &.{1},
        .oxygen_demand_satisfaction_fraction = &.{1},
        .respiration_oxygen_demand_g_o = &.{0},
        .non_band_oxidation_g = &.{0},
        .band_oxidation_g = &.{0},
        .fermenter_carbon_dioxide_yield_g_c_per_g_c = 0.2,
        .fermenter_acetate_yield_g_c_per_g_c = 0.2,
        .fermenter_hydrogen_yield_g_h_per_g_c = 0.111,
        .acetotrophic_carbon_dioxide_yield_g_c_per_g_c = 0.5,
        .acetotrophic_methane_yield_g_c_per_g_c = 0.5,
        .hydrogenotrophic_methane_yield_g_c_per_g_c = 1,
        .hydrogenotrophic_hydrogen_uptake_g_h_per_g_c = 0.667,
        .carbon_balance_tolerance_g_c = 1e-12,
    }));
    try std.testing.expectEqual(@as(f64, 7), state.carbon_dioxide_production_g_c[0]);
}

test "NITRO 1616-1665 product overflow preserves prior state" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.carbon_dioxide_production_g_c[0] = 7;
    try std.testing.expectError(
        error.NonFiniteRespirationProductResult,
        calculate(&state, .{
            .metabolism = &.{.fermenter},
            .autotrophic_complex = &.{false},
            .oxygen_unlimited_respiration_g_c = &.{
                std.math.floatMax(f64),
            },
            .oxygen_demand_satisfaction_fraction = &.{1},
            .respiration_oxygen_demand_g_o = &.{0},
            .non_band_oxidation_g = &.{0},
            .band_oxidation_g = &.{0},
            .fermenter_carbon_dioxide_yield_g_c_per_g_c = std.math.floatMax(f64),
            .fermenter_acetate_yield_g_c_per_g_c = 0,
            .fermenter_hydrogen_yield_g_h_per_g_c = 0,
            .acetotrophic_carbon_dioxide_yield_g_c_per_g_c = 0.5,
            .acetotrophic_methane_yield_g_c_per_g_c = 0.5,
            .hydrogenotrophic_methane_yield_g_c_per_g_c = 1,
            .hydrogenotrophic_hydrogen_uptake_g_h_per_g_c = 0.667,
            .carbon_balance_tolerance_g_c = 1e-12,
        }),
    );
    try std.testing.expectEqual(
        @as(f64, 7),
        state.carbon_dioxide_production_g_c[0],
    );
}
