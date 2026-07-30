const std = @import("std");

pub const Inputs = struct {
    enabled: []const bool,
    nonstructural_carbon_g_c: []const f64,
    nonstructural_nitrogen_g_n: []const f64,
    maximum_nitrogen_to_carbon_ratio_g_n_per_g_c: []const f64,
    growth_respiration_g_c: []const f64,
    aqueous_dinitrogen_concentration_g_n_m3: []const f64,
    fixation_yield_g_n_per_g_c: []const f64,
    nonstructural_to_structural_rate_h: f64,
    dinitrogen_half_saturation_g_n_m3: f64,
    negligible_respiration_g_c: f64,
    biochemical_time_fraction_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    nitrogen_deficit_g_n: []f64,
    required_fixation_respiration_g_c: []f64,
    respiration_demand_share_g_c: []f64,
    aqueous_dinitrogen_response: []f64,
    substrate_limited_respiration_g_c: []f64,
    fixation_respiration_g_c: []f64,
    fixed_dinitrogen_g_n: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidNitrogenFixationDimensions;
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

/// Exact NITRO.F 2511--2553 nonsymbiotic aerobic/anaerobic N2 fixation.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([7]f64, state.unit_count);
    defer state.allocator.free(temporary);
    for (0..state.unit_count) |unit| {
        if (!inputs.enabled[unit]) {
            temporary[unit] = @splat(0);
            continue;
        }
        const carbon = @max(0, inputs.nonstructural_carbon_g_c[unit]);
        const nitrogen = @max(0, inputs.nonstructural_nitrogen_g_n[unit]);
        const deficit = @max(0, carbon *
            inputs.maximum_nitrogen_to_carbon_ratio_g_n_per_g_c[unit] - nitrogen);
        const required_respiration = deficit / inputs.fixation_yield_g_n_per_g_c[unit];
        var demand_share: f64 = 0;
        var nitrogen_response: f64 = 0;
        var substrate_limit: f64 = 0;
        var fixation_respiration: f64 = 0;
        if (inputs.growth_respiration_g_c[unit] > inputs.negligible_respiration_g_c) {
            const denominator = inputs.growth_respiration_g_c[unit] + required_respiration;
            demand_share = inputs.growth_respiration_g_c[unit] *
                required_respiration / denominator;
            nitrogen_response = inputs.aqueous_dinitrogen_concentration_g_n_m3[unit] /
                (inputs.aqueous_dinitrogen_concentration_g_n_m3[unit] +
                    inputs.dinitrogen_half_saturation_g_n_m3);
            substrate_limit = inputs.nonstructural_to_structural_rate_h * carbon *
                inputs.biochemical_time_fraction_h;
            fixation_respiration = @min(demand_share * nitrogen_response, substrate_limit);
        }
        temporary[unit] = .{
            deficit,
            required_respiration,
            demand_share,
            nitrogen_response,
            substrate_limit,
            fixation_respiration,
            fixation_respiration * inputs.fixation_yield_g_n_per_g_c[unit],
        };
        for (temporary[unit]) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteNitrogenFixationResult;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, unit| @field(state, field.name)[unit] = values[index];
    };
}

pub fn sourceEnabled(zero_based_complex: usize, zero_based_population: usize) bool {
    return zero_based_complex <= 3 and
        (zero_based_population == 5 or zero_based_population == 6);
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
    if (inputs.enabled.len != n) return error.InvalidNitrogenFixationDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != n) return error.InvalidNitrogenFixationDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidNitrogenFixationInput;
        const clamped = comptime std.mem.eql(u8, field.name, "nonstructural_carbon_g_c") or
            std.mem.eql(u8, field.name, "nonstructural_nitrogen_g_n");
        if (!clamped) for (values) |value| if (value < 0)
            return error.InvalidNitrogenFixationInput;
    };
    for (inputs.fixation_yield_g_n_per_g_c) |value| if (value == 0)
        return error.InvalidNitrogenFixationInput;
    inline for (.{
        inputs.nonstructural_to_structural_rate_h,
        inputs.dinitrogen_half_saturation_g_n_m3,
        inputs.negligible_respiration_g_c,
        inputs.biochemical_time_fraction_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidNitrogenFixationInput;
    if (inputs.dinitrogen_half_saturation_g_n_m3 == 0)
        return error.InvalidNitrogenFixationInput;
}

fn fixture() Inputs {
    return .{
        .enabled = &.{true},
        .nonstructural_carbon_g_c = &.{10},
        .nonstructural_nitrogen_g_n = &.{1},
        .maximum_nitrogen_to_carbon_ratio_g_n_per_g_c = &.{0.2},
        .growth_respiration_g_c = &.{4},
        .aqueous_dinitrogen_concentration_g_n_m3 = &.{3},
        .fixation_yield_g_n_per_g_c = &.{0.5},
        .nonstructural_to_structural_rate_h = 0.1,
        .dinitrogen_half_saturation_g_n_m3 = 1,
        .negligible_respiration_g_c = 1e-12,
        .biochemical_time_fraction_h = 1,
    };
}

test "source roles match aerobic and anaerobic nonsymbiotic fixers" {
    try std.testing.expect(sourceEnabled(3, 5));
    try std.testing.expect(sourceEnabled(0, 6));
    try std.testing.expect(!sourceEnabled(4, 5));
    try std.testing.expect(!sourceEnabled(3, 4));
}

test "fixation preserves respiration sharing nitrogen response and carbon ceiling" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(1, state.nitrogen_deficit_g_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(2, state.required_fixation_respiration_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(1, state.fixation_respiration_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.5, state.fixed_dinitrogen_g_n[0], 1e-12);
}

test "negligible growth respiration disables fixation" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.growth_respiration_g_c = &.{0};
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.fixed_dinitrogen_g_n[0]);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.fixed_dinitrogen_g_n[0] = 7;
    var inputs = fixture();
    inputs.aqueous_dinitrogen_concentration_g_n_m3 = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidNitrogenFixationInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.fixed_dinitrogen_g_n[0]);
}

test "NITRO 2511-2553 derived overflow preserves fixation state" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.fixed_dinitrogen_g_n[0] = 7;
    var inputs = fixture();
    inputs.nonstructural_carbon_g_c = &.{std.math.floatMax(f64)};
    inputs.maximum_nitrogen_to_carbon_ratio_g_n_per_g_c = &.{2};
    try std.testing.expectError(
        error.NonFiniteNitrogenFixationResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.fixed_dinitrogen_g_n[0]);
}
