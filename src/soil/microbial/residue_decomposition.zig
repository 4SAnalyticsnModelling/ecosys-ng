const std = @import("std");

pub const Inputs = struct {
    residue_fraction_count: usize,
    residue_carbon_g_c: []const f64,
    residue_nitrogen_g_n: []const f64,
    residue_phosphorus_g_p: []const f64,
    specific_decomposition_rate_h: []const f64,
    colonized_soil_organic_carbon_g_c: []const f64,
    microbial_activity_respiration_g_c: []const f64,
    substrate_density_response: []const f64,
    dissolved_carbon_product_response: []const f64,
    nitrogen_limitation_fraction: []const f64,
    phosphorus_limitation_fraction: []const f64,
    temperature_response: []const f64,
    negligible_carbon_g_c: f64,
    biochemical_time_fraction_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    item_count: usize,
    nitrogen_to_carbon_ratio_g_n_per_g_c: []f64,
    phosphorus_to_carbon_ratio_g_p_per_g_c: []f64,
    decomposed_carbon_g_c: []f64,
    decomposed_nitrogen_g_n: []f64,
    decomposed_phosphorus_g_p: []f64,

    pub fn init(allocator: std.mem.Allocator, complex_count: usize, fraction_count: usize) !State {
        if (complex_count == 0 or fraction_count == 0)
            return error.InvalidResidueDecompositionDimensions;
        const item_count = std.math.mul(usize, complex_count, fraction_count) catch
            return error.InvalidResidueDecompositionDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.item_count = item_count;
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
            @field(state, field.name) = try allocator.alloc(f64, item_count);
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

/// Exact NITRO.F 3352--3398 microbial-residue C/N/P decomposition.
/// Traversal is complex-major and residue-fraction-minor.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const complex_count = state.item_count / inputs.residue_fraction_count;
    const temporary = try state.allocator.alloc([5]f64, state.item_count);
    defer state.allocator.free(temporary);
    for (0..complex_count) |complex| {
        const colonized = inputs.colonized_soil_organic_carbon_g_c[complex];
        for (0..inputs.residue_fraction_count) |fraction| {
            const item = complex * inputs.residue_fraction_count + fraction;
            if (colonized <= inputs.negligible_carbon_g_c or
                inputs.residue_carbon_g_c[item] <= inputs.negligible_carbon_g_c)
            {
                temporary[item] = @splat(0);
                continue;
            }
            const nitrogen_ratio = @max(0, inputs.residue_nitrogen_g_n[item] / inputs.residue_carbon_g_c[item]);
            const phosphorus_ratio = @max(0, inputs.residue_phosphorus_g_p[item] / inputs.residue_carbon_g_c[item]);
            const activity_limited = inputs.specific_decomposition_rate_h[fraction] *
                inputs.microbial_activity_respiration_g_c[complex] *
                inputs.substrate_density_response[complex] *
                inputs.dissolved_carbon_product_response[complex] *
                inputs.temperature_response[complex] *
                inputs.residue_carbon_g_c[item] / colonized;
            const decomposed_carbon = @max(0, @min(
                inputs.residue_carbon_g_c[item] * inputs.biochemical_time_fraction_h,
                activity_limited,
            ));
            const decomposed_nitrogen = @max(0, @min(
                inputs.residue_nitrogen_g_n[item] * inputs.biochemical_time_fraction_h,
                nitrogen_ratio * decomposed_carbon,
            )) / inputs.nitrogen_limitation_fraction[complex];
            const decomposed_phosphorus = @max(0, @min(
                inputs.residue_phosphorus_g_p[item] * inputs.biochemical_time_fraction_h,
                phosphorus_ratio * decomposed_carbon,
            )) / inputs.phosphorus_limitation_fraction[complex];
            temporary[item] = .{
                nitrogen_ratio,      phosphorus_ratio,      decomposed_carbon,
                decomposed_nitrogen, decomposed_phosphorus,
            };
            for (temporary[item]) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteResidueDecompositionResult;
        }
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, item| @field(state, field.name)[item] = values[index];
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
    if (inputs.residue_fraction_count == 0 or state.item_count % inputs.residue_fraction_count != 0)
        return error.InvalidResidueDecompositionDimensions;
    const complexes = state.item_count / inputs.residue_fraction_count;
    inline for (.{
        inputs.residue_carbon_g_c,     inputs.residue_nitrogen_g_n,
        inputs.residue_phosphorus_g_p,
    }) |values| {
        if (values.len != state.item_count) return error.InvalidResidueDecompositionDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidResidueDecompositionInput;
    }
    if (inputs.specific_decomposition_rate_h.len != inputs.residue_fraction_count)
        return error.InvalidResidueDecompositionDimensions;
    for (inputs.specific_decomposition_rate_h) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidResidueDecompositionInput;
    inline for (.{
        inputs.colonized_soil_organic_carbon_g_c, inputs.microbial_activity_respiration_g_c,
        inputs.substrate_density_response,        inputs.dissolved_carbon_product_response,
        inputs.nitrogen_limitation_fraction,      inputs.phosphorus_limitation_fraction,
        inputs.temperature_response,
    }) |values| {
        if (values.len != complexes) return error.InvalidResidueDecompositionDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidResidueDecompositionInput;
    }
    inline for (.{
        inputs.nitrogen_limitation_fraction, inputs.phosphorus_limitation_fraction,
    }) |values| for (values) |value| if (value == 0)
        return error.InvalidResidueDecompositionInput;
    inline for (.{
        inputs.negligible_carbon_g_c, inputs.biochemical_time_fraction_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidResidueDecompositionInput;
}

fn fixture() Inputs {
    return .{
        .residue_fraction_count = 2,
        .residue_carbon_g_c = &.{ 10, 20 },
        .residue_nitrogen_g_n = &.{ 2, 2 },
        .residue_phosphorus_g_p = &.{ 1, 1 },
        .specific_decomposition_rate_h = &.{ 0.1, 0.05 },
        .colonized_soil_organic_carbon_g_c = &.{10},
        .microbial_activity_respiration_g_c = &.{4},
        .substrate_density_response = &.{0.5},
        .dissolved_carbon_product_response = &.{0.5},
        .nitrogen_limitation_fraction = &.{0.5},
        .phosphorus_limitation_fraction = &.{0.5},
        .temperature_response = &.{1},
        .negligible_carbon_g_c = 1e-12,
        .biochemical_time_fraction_h = 1,
    };
}

test "runtime residue fractions decompose with shared complex controls" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(0.2, state.nitrogen_to_carbon_ratio_g_n_per_g_c[0], 1e-12);
    try std.testing.expect(state.decomposed_carbon_g_c[0] > 0);
    try std.testing.expect(state.decomposed_carbon_g_c[1] > 0);
}

test "residue nutrient decomposition preserves source limitation divisors" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(
        state.decomposed_carbon_g_c[0] * 0.2 / 0.5,
        state.decomposed_nitrogen_g_n[0],
        1e-12,
    );
}

test "absent colonized substrate zeros residue decomposition" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.colonized_soil_organic_carbon_g_c = &.{0};
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.decomposed_carbon_g_c[0]);
}

test "NITRO outer and residue threshold equality zero all outputs" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.residue_carbon_g_c = &.{ 10, 20, 1e-12, 20 };
    inputs.residue_nitrogen_g_n = &.{ 2, 2, 2, 2 };
    inputs.residue_phosphorus_g_p = &.{ 1, 1, 1, 1 };
    inputs.colonized_soil_organic_carbon_g_c = &.{ 1e-12, 10 };
    inputs.microbial_activity_respiration_g_c = &.{ 4, 4 };
    inputs.substrate_density_response = &.{ 0.5, 0.5 };
    inputs.dissolved_carbon_product_response = &.{ 0.5, 0.5 };
    inputs.nitrogen_limitation_fraction = &.{ 0.5, 0.5 };
    inputs.phosphorus_limitation_fraction = &.{ 0.5, 0.5 };
    inputs.temperature_response = &.{ 1, 1 };
    try calculate(&state, inputs);
    for (0..3) |item| inline for (.{
        state.nitrogen_to_carbon_ratio_g_n_per_g_c,
        state.phosphorus_to_carbon_ratio_g_p_per_g_c,
        state.decomposed_carbon_g_c,
        state.decomposed_nitrogen_g_n,
        state.decomposed_phosphorus_g_p,
    }) |values| try std.testing.expectEqual(0, values[item]);
    try std.testing.expect(state.decomposed_carbon_g_c[3] > 0);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.decomposed_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.temperature_response = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidResidueDecompositionInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.decomposed_carbon_g_c[0]);
}

test "NITRO 3352-3398 late derived overflow preserves all residue outputs" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.decomposed_carbon_g_c[0] = 7;
    state.decomposed_phosphorus_g_p[1] = 11;
    var inputs = fixture();
    inputs.specific_decomposition_rate_h =
        &.{ 0.1, std.math.floatMax(f64) };
    inputs.biochemical_time_fraction_h = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteResidueDecompositionResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.decomposed_carbon_g_c[0]);
    try std.testing.expectEqual(11, state.decomposed_phosphorus_g_p[1]);
}
