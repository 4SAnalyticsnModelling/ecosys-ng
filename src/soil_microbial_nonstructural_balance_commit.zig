const std = @import("std");

pub const Inputs = struct {
    structural_fraction_count: usize,
    enabled: []const bool,
    surface_litter_layer: []const bool,
    nonstructural_carbon_g_c: []const f64,
    nonstructural_nitrogen_g_n: []const f64,
    nonstructural_phosphorus_g_p: []const f64,
    carbon_dioxide_emission_g_c: []const f64,
    total_carbon_uptake_g_c: []const f64,
    aerobic_respiration_g_c: []const f64,
    denitrification_respiration_g_c: []const f64,
    nitrogen_fixation_respiration_g_c: []const f64,
    dissolved_organic_nitrogen_uptake_g_n: []const f64,
    dissolved_organic_phosphorus_uptake_g_p: []const f64,
    non_band_ammonium_exchange_g_n: []const f64,
    band_ammonium_exchange_g_n: []const f64,
    non_band_nitrate_exchange_g_n: []const f64,
    band_nitrate_exchange_g_n: []const f64,
    fixed_dinitrogen_g_n: []const f64,
    non_band_dihydrogen_phosphate_exchange_g_p: []const f64,
    band_dihydrogen_phosphate_exchange_g_p: []const f64,
    non_band_hydrogen_phosphate_exchange_g_p: []const f64,
    band_hydrogen_phosphate_exchange_g_p: []const f64,
    surface_ammonium_exchange_g_n: []const f64,
    surface_nitrate_exchange_g_n: []const f64,
    surface_dihydrogen_phosphate_exchange_g_p: []const f64,
    surface_hydrogen_phosphate_exchange_g_p: []const f64,
    structural_carbon_assimilation_g_c: []const f64,
    structural_nitrogen_assimilation_g_n: []const f64,
    structural_phosphorus_assimilation_g_p: []const f64,
    basal_recycled_carbon_g_c: []const f64,
    basal_recycled_nitrogen_g_n: []const f64,
    basal_recycled_phosphorus_g_p: []const f64,
    senescence_recycled_carbon_g_c: []const f64,
    senescence_recycled_nitrogen_g_n: []const f64,
    senescence_recycled_phosphorus_g_p: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    nonstructural_carbon_g_c: []f64,
    nonstructural_nitrogen_g_n: []f64,
    nonstructural_phosphorus_g_p: []f64,
    carbon_dioxide_emission_g_c: []f64,
    net_carbon_uptake_after_respiration_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidNonstructuralBalanceDimensions;
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
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            @field(state, field.name) = try allocator.alloc(f64, unit_count);
            @memset(@field(state, field.name), 0);
            allocated += 1;
        };
        return state;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

/// Exact NITRO.F 3794--3858 nonstructural microbial C/N/P balance.
/// Units must be supplied in enclosing source K, then N order.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([5]f64, state.unit_count);
    defer state.allocator.free(temporary);
    for (0..state.unit_count) |unit| {
        if (!inputs.enabled[unit]) {
            temporary[unit] = .{
                inputs.nonstructural_carbon_g_c[unit],
                inputs.nonstructural_nitrogen_g_n[unit],
                inputs.nonstructural_phosphorus_g_p[unit],
                inputs.carbon_dioxide_emission_g_c[unit],
                0,
            };
            continue;
        }
        const net_carbon = inputs.total_carbon_uptake_g_c[unit] -
            inputs.aerobic_respiration_g_c[unit] -
            inputs.denitrification_respiration_g_c[unit] -
            inputs.nitrogen_fixation_respiration_g_c[unit];
        var carbon = inputs.nonstructural_carbon_g_c[unit];
        var nitrogen = inputs.nonstructural_nitrogen_g_n[unit];
        var phosphorus = inputs.nonstructural_phosphorus_g_p[unit];
        var carbon_dioxide = inputs.carbon_dioxide_emission_g_c[unit] +
            inputs.nitrogen_fixation_respiration_g_c[unit];
        for (0..inputs.structural_fraction_count) |fraction| {
            const item = unit * inputs.structural_fraction_count + fraction;
            carbon = carbon - inputs.structural_carbon_assimilation_g_c[item] +
                inputs.basal_recycled_carbon_g_c[item];
            nitrogen = nitrogen - inputs.structural_nitrogen_assimilation_g_n[item] +
                inputs.basal_recycled_nitrogen_g_n[item] +
                inputs.senescence_recycled_nitrogen_g_n[item];
            phosphorus = phosphorus -
                inputs.structural_phosphorus_assimilation_g_p[item] +
                inputs.basal_recycled_phosphorus_g_p[item] +
                inputs.senescence_recycled_phosphorus_g_p[item];
            carbon_dioxide = carbon_dioxide +
                inputs.senescence_recycled_carbon_g_c[item];
            try validateIntermediate(carbon, nitrogen, phosphorus, carbon_dioxide);
        }
        carbon = carbon + net_carbon;
        nitrogen = nitrogen + inputs.dissolved_organic_nitrogen_uptake_g_n[unit] +
            inputs.non_band_ammonium_exchange_g_n[unit] +
            inputs.band_ammonium_exchange_g_n[unit] +
            inputs.non_band_nitrate_exchange_g_n[unit] +
            inputs.band_nitrate_exchange_g_n[unit] +
            inputs.fixed_dinitrogen_g_n[unit];
        phosphorus = phosphorus +
            inputs.dissolved_organic_phosphorus_uptake_g_p[unit] +
            inputs.non_band_dihydrogen_phosphate_exchange_g_p[unit] +
            inputs.band_dihydrogen_phosphate_exchange_g_p[unit] +
            inputs.non_band_hydrogen_phosphate_exchange_g_p[unit] +
            inputs.band_hydrogen_phosphate_exchange_g_p[unit];
        try validateIntermediate(carbon, nitrogen, phosphorus, carbon_dioxide);
        if (inputs.surface_litter_layer[unit]) {
            nitrogen = nitrogen + inputs.surface_ammonium_exchange_g_n[unit] +
                inputs.surface_nitrate_exchange_g_n[unit];
            phosphorus = phosphorus +
                inputs.surface_dihydrogen_phosphate_exchange_g_p[unit] +
                inputs.surface_hydrogen_phosphate_exchange_g_p[unit];
            try validateIntermediate(carbon, nitrogen, phosphorus, carbon_dioxide);
        }
        temporary[unit] = .{ carbon, nitrogen, phosphorus, carbon_dioxide, net_carbon };
    }
    for (temporary) |values| {
        inline for (values[0..4]) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidNonstructuralBalanceResult;
        if (!std.math.isFinite(values[4])) return error.InvalidNonstructuralBalanceResult;
    }
    for (temporary, 0..) |values, unit| {
        state.nonstructural_carbon_g_c[unit] = values[0];
        state.nonstructural_nitrogen_g_n[unit] = values[1];
        state.nonstructural_phosphorus_g_p[unit] = values[2];
        state.carbon_dioxide_emission_g_c[unit] = values[3];
        state.net_carbon_uptake_after_respiration_g_c[unit] = values[4];
    }
}

fn validateIntermediate(carbon: f64, nitrogen: f64, phosphorus: f64, carbon_dioxide: f64) !void {
    inline for (.{ carbon, nitrogen, phosphorus, carbon_dioxide }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidNonstructuralBalanceResult;
}

fn validate(state: *const State, inputs: Inputs) !void {
    comptime {
        @setEvalBranchQuota(5000);
    }
    if (inputs.structural_fraction_count == 0 or inputs.enabled.len != state.unit_count or
        inputs.surface_litter_layer.len != state.unit_count)
        return error.InvalidNonstructuralBalanceDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        const fraction_field = comptime std.mem.startsWith(u8, field.name, "structural_") or
            std.mem.startsWith(u8, field.name, "basal_recycled_") or
            std.mem.startsWith(u8, field.name, "senescence_recycled_");
        const expected = if (fraction_field)
            state.unit_count * inputs.structural_fraction_count
        else
            state.unit_count;
        if (values.len != expected) return error.InvalidNonstructuralBalanceDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidNonstructuralBalanceInput;
        const signed_exchange = comptime std.mem.indexOf(u8, field.name, "_exchange_") != null;
        if (!signed_exchange) for (values) |value| if (value < 0)
            return error.InvalidNonstructuralBalanceInput;
    };
}

fn fixture() Inputs {
    return .{
        .structural_fraction_count = 2,
        .enabled = &.{true},
        .surface_litter_layer = &.{true},
        .nonstructural_carbon_g_c = &.{10},
        .nonstructural_nitrogen_g_n = &.{3},
        .nonstructural_phosphorus_g_p = &.{2},
        .carbon_dioxide_emission_g_c = &.{1},
        .total_carbon_uptake_g_c = &.{8},
        .aerobic_respiration_g_c = &.{2},
        .denitrification_respiration_g_c = &.{1},
        .nitrogen_fixation_respiration_g_c = &.{1},
        .dissolved_organic_nitrogen_uptake_g_n = &.{0.5},
        .dissolved_organic_phosphorus_uptake_g_p = &.{0.2},
        .non_band_ammonium_exchange_g_n = &.{0.1},
        .band_ammonium_exchange_g_n = &.{0.1},
        .non_band_nitrate_exchange_g_n = &.{0.1},
        .band_nitrate_exchange_g_n = &.{0.1},
        .fixed_dinitrogen_g_n = &.{0.2},
        .non_band_dihydrogen_phosphate_exchange_g_p = &.{0.05},
        .band_dihydrogen_phosphate_exchange_g_p = &.{0.05},
        .non_band_hydrogen_phosphate_exchange_g_p = &.{0.05},
        .band_hydrogen_phosphate_exchange_g_p = &.{0.05},
        .surface_ammonium_exchange_g_n = &.{0.1},
        .surface_nitrate_exchange_g_n = &.{0.1},
        .surface_dihydrogen_phosphate_exchange_g_p = &.{0.05},
        .surface_hydrogen_phosphate_exchange_g_p = &.{0.05},
        .structural_carbon_assimilation_g_c = &.{ 1, 1 },
        .structural_nitrogen_assimilation_g_n = &.{ 0.2, 0.2 },
        .structural_phosphorus_assimilation_g_p = &.{ 0.1, 0.1 },
        .basal_recycled_carbon_g_c = &.{ 0.5, 0.5 },
        .basal_recycled_nitrogen_g_n = &.{ 0.1, 0.1 },
        .basal_recycled_phosphorus_g_p = &.{ 0.05, 0.05 },
        .senescence_recycled_carbon_g_c = &.{ 0.25, 0.25 },
        .senescence_recycled_nitrogen_g_n = &.{ 0.05, 0.05 },
        .senescence_recycled_phosphorus_g_p = &.{ 0.025, 0.025 },
    };
}

test "nonstructural balance preserves uptake respiration transfer and recycling" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectEqual(4, state.net_carbon_uptake_after_respiration_g_c[0]);
    try std.testing.expectApproxEqAbs(13, state.nonstructural_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(2.5, state.carbon_dioxide_emission_g_c[0], 1e-12);
}

test "surface mineral exchanges enter only surface litter balances" {
    var surface = try State.init(std.testing.allocator, 1);
    defer surface.deinit();
    var soil = try State.init(std.testing.allocator, 1);
    defer soil.deinit();
    const surface_inputs = fixture();
    try calculate(&surface, surface_inputs);
    var soil_inputs = fixture();
    soil_inputs.surface_litter_layer = &.{false};
    try calculate(&soil, soil_inputs);
    try std.testing.expectApproxEqAbs(0.2, surface.nonstructural_nitrogen_g_n[0] - soil.nonstructural_nitrogen_g_n[0], 1e-12);
}

test "signed mineralization exchanges can reduce nonstructural nutrient pools" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.non_band_ammonium_exchange_g_n = &.{-0.5};
    try calculate(&state, inputs);
    try std.testing.expect(state.nonstructural_nitrogen_g_n[0] > 0);
}

test "negative final pool fails atomically" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.nonstructural_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.total_carbon_uptake_g_c = &.{0};
    inputs.aerobic_respiration_g_c = &.{20};
    try std.testing.expectError(error.InvalidNonstructuralBalanceResult, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.nonstructural_carbon_g_c[0]);
}

test "disabled source unit retains pools and zeroes unpublished CGROMC" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.enabled = &.{false};
    try calculate(&state, inputs);
    try std.testing.expectEqual(10, state.nonstructural_carbon_g_c[0]);
    try std.testing.expectEqual(3, state.nonstructural_nitrogen_g_n[0]);
    try std.testing.expectEqual(2, state.nonstructural_phosphorus_g_p[0]);
    try std.testing.expectEqual(1, state.carbon_dioxide_emission_g_c[0]);
    try std.testing.expectEqual(0, state.net_carbon_uptake_after_respiration_g_c[0]);
}

test "NITRO 3794-3858 late surface overflow fails atomically" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.nonstructural_carbon_g_c[0] = 7;
    state.nonstructural_phosphorus_g_p[0] = 11;
    var inputs = fixture();
    inputs.surface_dihydrogen_phosphate_exchange_g_p =
        &.{std.math.floatMax(f64)};
    inputs.surface_hydrogen_phosphate_exchange_g_p =
        &.{std.math.floatMax(f64)};
    try std.testing.expectError(
        error.InvalidNonstructuralBalanceResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.nonstructural_carbon_g_c[0]);
    try std.testing.expectEqual(11, state.nonstructural_phosphorus_g_p[0]);
}
