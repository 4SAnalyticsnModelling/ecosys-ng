const std = @import("std");

pub const Inputs = struct {
    nonstructural_carbon_g_c: []const f64,
    nonstructural_nitrogen_g_n: []const f64,
    nonstructural_phosphorus_g_p: []const f64,
    labile_nitrogen_to_carbon_ratio_g_n_per_g_c: []const f64,
    labile_phosphorus_to_carbon_ratio_g_p_per_g_c: []const f64,
    minimum_carbon_recycling_fraction: f64,
    carbon_recycling_response_fraction: f64,
    maximum_nitrogen_recycling_fraction: f64,
    maximum_phosphorus_recycling_fraction: f64,
    negligible_nonstructural_carbon_g_c: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    carbon_nutrient_sufficiency: []f64,
    nitrogen_carbon_sufficiency: []f64,
    phosphorus_carbon_sufficiency: []f64,
    carbon_recycling_fraction: []f64,
    nitrogen_recycling_fraction: []f64,
    phosphorus_recycling_fraction: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidNutrientRecyclingDimensions;
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

/// Exact NITRO.F 2616--2646 C, N, and P recycling controls.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([6]f64, state.unit_count);
    defer state.allocator.free(temporary);
    for (0..state.unit_count) |unit| {
        const carbon = inputs.nonstructural_carbon_g_c[unit];
        const nitrogen = inputs.nonstructural_nitrogen_g_n[unit];
        const phosphorus = inputs.nonstructural_phosphorus_g_p[unit];
        var carbon_sufficiency: f64 = 1;
        var nitrogen_sufficiency: f64 = 0;
        var phosphorus_sufficiency: f64 = 0;
        if (carbon > inputs.negligible_nonstructural_carbon_g_c) {
            const nitrogen_term = nitrogen /
                (nitrogen + carbon *
                    inputs.labile_nitrogen_to_carbon_ratio_g_n_per_g_c[unit]);
            const phosphorus_term = phosphorus /
                (phosphorus + carbon *
                    inputs.labile_phosphorus_to_carbon_ratio_g_p_per_g_c[unit]);
            carbon_sufficiency = @max(0, @min(1, @min(nitrogen_term, phosphorus_term)));
            nitrogen_sufficiency = @max(0, @min(1, carbon / (carbon + nitrogen /
                inputs.labile_nitrogen_to_carbon_ratio_g_n_per_g_c[unit])));
            phosphorus_sufficiency = @max(0, @min(1, carbon / (carbon + phosphorus /
                inputs.labile_phosphorus_to_carbon_ratio_g_p_per_g_c[unit])));
        }
        temporary[unit] = .{
            carbon_sufficiency,
            nitrogen_sufficiency,
            phosphorus_sufficiency,
            inputs.minimum_carbon_recycling_fraction +
                carbon_sufficiency * inputs.carbon_recycling_response_fraction,
            nitrogen_sufficiency * inputs.maximum_nitrogen_recycling_fraction,
            phosphorus_sufficiency * inputs.maximum_phosphorus_recycling_fraction,
        };
        for (temporary[unit]) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteNutrientRecyclingResult;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, unit| @field(state, field.name)[unit] = values[index];
    };
}

fn stateFieldIndex(comptime name: []const u8) usize {
    var index: usize = 0;
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        if (comptime std.mem.eql(u8, field.name, name)) return index;
        index += 1;
    };
    unreachable;
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != n) return error.InvalidNutrientRecyclingDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidNutrientRecyclingInput;
    };
    inline for (.{
        inputs.labile_nitrogen_to_carbon_ratio_g_n_per_g_c,
        inputs.labile_phosphorus_to_carbon_ratio_g_p_per_g_c,
    }) |values| for (values) |value| if (value == 0)
        return error.InvalidNutrientRecyclingInput;
    inline for (.{
        inputs.minimum_carbon_recycling_fraction,
        inputs.carbon_recycling_response_fraction,
        inputs.maximum_nitrogen_recycling_fraction,
        inputs.maximum_phosphorus_recycling_fraction,
        inputs.negligible_nonstructural_carbon_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidNutrientRecyclingInput;
}

fn fixture() Inputs {
    return .{
        .nonstructural_carbon_g_c = &.{10},
        .nonstructural_nitrogen_g_n = &.{2},
        .nonstructural_phosphorus_g_p = &.{1},
        .labile_nitrogen_to_carbon_ratio_g_n_per_g_c = &.{0.2},
        .labile_phosphorus_to_carbon_ratio_g_p_per_g_c = &.{0.1},
        .minimum_carbon_recycling_fraction = 0.1,
        .carbon_recycling_response_fraction = 0.5,
        .maximum_nitrogen_recycling_fraction = 0.8,
        .maximum_phosphorus_recycling_fraction = 0.6,
        .negligible_nonstructural_carbon_g_c = 1e-12,
    };
}

test "nutrient ratios control separate recycling fractions" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(0.5, state.carbon_nutrient_sufficiency[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.5, state.nitrogen_carbon_sufficiency[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.5, state.phosphorus_carbon_sufficiency[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.35, state.carbon_recycling_fraction[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.4, state.nitrogen_recycling_fraction[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.3, state.phosphorus_recycling_fraction[0], 1e-12);
}

test "negligible carbon selects source fallback recycling controls" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.nonstructural_carbon_g_c = &.{0};
    try calculate(&state, inputs);
    try std.testing.expectEqual(1, state.carbon_nutrient_sufficiency[0]);
    try std.testing.expectEqual(0, state.nitrogen_recycling_fraction[0]);
    try std.testing.expectEqual(0, state.phosphorus_recycling_fraction[0]);
}

test "runtime units use independent stoichiometry" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.nonstructural_carbon_g_c = &.{ 10, 10 };
    inputs.nonstructural_nitrogen_g_n = &.{ 2, 0.2 };
    inputs.nonstructural_phosphorus_g_p = &.{ 1, 0.1 };
    inputs.labile_nitrogen_to_carbon_ratio_g_n_per_g_c = &.{ 0.2, 0.2 };
    inputs.labile_phosphorus_to_carbon_ratio_g_p_per_g_c = &.{ 0.1, 0.1 };
    try calculate(&state, inputs);
    try std.testing.expect(state.carbon_recycling_fraction[1] <
        state.carbon_recycling_fraction[0]);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.carbon_recycling_fraction[0] = 7;
    var inputs = fixture();
    inputs.nonstructural_phosphorus_g_p = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidNutrientRecyclingInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.carbon_recycling_fraction[0]);
}

test "NITRO 2616-2646 derived overflow preserves recycling state" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.carbon_recycling_fraction[0] = 7;
    var inputs = fixture();
    inputs.minimum_carbon_recycling_fraction = std.math.floatMax(f64);
    inputs.carbon_recycling_response_fraction = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteNutrientRecyclingResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.carbon_recycling_fraction[0]);
}
