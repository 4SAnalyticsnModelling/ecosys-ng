const std = @import("std");

pub const Inputs = struct {
    substrate_class_count: usize,
    salt_chemistry_enabled: bool,
    substrate_carbon_g_c: []const f64,
    substrate_nitrogen_g_n: []const f64,
    substrate_phosphorus_g_p: []const f64,
    colonizing_active_biomass_g_c: []const f64,
    specific_decomposition_rate_per_g_c_h: []const f64,
    colonized_soil_organic_carbon_g_c: []const f64,
    microbial_activity_respiration_g_c: []const f64,
    substrate_density_response: []const f64,
    dissolved_carbon_product_response: []const f64,
    nitrogen_limitation_fraction: []const f64,
    phosphorus_limitation_fraction: []const f64,
    temperature_response: []const f64,
    negligible_substrate_carbon_g_c: f64,
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
    lignin_decomposition_for_acidity_g_c: f64,

    pub fn init(allocator: std.mem.Allocator, complex_count: usize, substrate_class_count: usize) !State {
        if (complex_count == 0 or substrate_class_count == 0)
            return error.InvalidSolidDecompositionDimensions;
        const item_count = std.math.mul(usize, complex_count, substrate_class_count) catch
            return error.InvalidSolidDecompositionDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.item_count = item_count;
        state.lignin_decomposition_for_acidity_g_c = 0;
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

/// Exact NITRO.F 3242--3298 solid substrate C/N/P decomposition.
/// Traversal is complex-major and structural-fraction-minor.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const complex_count = state.item_count / inputs.substrate_class_count;
    const temporary = try state.allocator.alloc([5]f64, state.item_count);
    defer state.allocator.free(temporary);
    var lignin_acidity: f64 = 0;
    for (0..complex_count) |complex| {
        const colonized = inputs.colonized_soil_organic_carbon_g_c[complex];
        if (colonized <= inputs.negligible_substrate_carbon_g_c) {
            for (0..inputs.substrate_class_count) |substrate|
                temporary[complex * inputs.substrate_class_count + substrate] = @splat(0);
            continue;
        }
        for (0..inputs.substrate_class_count) |substrate| {
            const item = complex * inputs.substrate_class_count + substrate;
            if (inputs.substrate_carbon_g_c[item] <= inputs.negligible_substrate_carbon_g_c) {
                temporary[item] = @splat(0);
                continue;
            }
            const nitrogen_ratio = @max(0, inputs.substrate_nitrogen_g_n[item] / inputs.substrate_carbon_g_c[item]);
            const phosphorus_ratio = @max(0, inputs.substrate_phosphorus_g_p[item] / inputs.substrate_carbon_g_c[item]);
            const activity_limited =
                inputs.specific_decomposition_rate_per_g_c_h[item] *
                inputs.microbial_activity_respiration_g_c[complex] *
                inputs.substrate_density_response[complex] *
                inputs.dissolved_carbon_product_response[complex] *
                inputs.temperature_response[complex] *
                inputs.colonizing_active_biomass_g_c[item] / colonized;
            const decomposed_carbon = @max(0, @min(
                inputs.colonizing_active_biomass_g_c[item] *
                    inputs.biochemical_time_fraction_h,
                activity_limited,
            ));
            const decomposed_nitrogen = @max(0, @min(
                inputs.substrate_nitrogen_g_n[item] * inputs.biochemical_time_fraction_h,
                nitrogen_ratio * decomposed_carbon,
            )) / inputs.nitrogen_limitation_fraction[complex];
            const decomposed_phosphorus = @max(0, @min(
                inputs.substrate_phosphorus_g_p[item] * inputs.biochemical_time_fraction_h,
                phosphorus_ratio * decomposed_carbon,
            )) / inputs.phosphorus_limitation_fraction[complex];
            if (inputs.salt_chemistry_enabled and
                sourceContributesToLigninAcidity(complex, substrate))
            {
                lignin_acidity = lignin_acidity + decomposed_carbon;
                if (!std.math.isFinite(lignin_acidity))
                    return error.NonFiniteSolidDecompositionResult;
            }
            temporary[item] = .{
                nitrogen_ratio,      phosphorus_ratio,      decomposed_carbon,
                decomposed_nitrogen, decomposed_phosphorus,
            };
            for (temporary[item]) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteSolidDecompositionResult;
        }
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, item| @field(state, field.name)[item] = values[index];
    };
    state.lignin_decomposition_for_acidity_g_c = lignin_acidity;
}

pub fn sourceContributesToLigninAcidity(
    zero_based_complex: usize,
    zero_based_substrate: usize,
) bool {
    return zero_based_complex <= 2 and zero_based_substrate == 3;
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
    if (inputs.substrate_class_count == 0 or state.item_count % inputs.substrate_class_count != 0)
        return error.InvalidSolidDecompositionDimensions;
    const complexes = state.item_count / inputs.substrate_class_count;
    inline for (.{
        inputs.substrate_carbon_g_c,                  inputs.substrate_nitrogen_g_n,
        inputs.substrate_phosphorus_g_p,              inputs.colonizing_active_biomass_g_c,
        inputs.specific_decomposition_rate_per_g_c_h,
    }) |values| {
        if (values.len != state.item_count) return error.InvalidSolidDecompositionDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSolidDecompositionInput;
    }
    inline for (.{
        inputs.colonized_soil_organic_carbon_g_c, inputs.microbial_activity_respiration_g_c,
        inputs.substrate_density_response,        inputs.dissolved_carbon_product_response,
        inputs.nitrogen_limitation_fraction,      inputs.phosphorus_limitation_fraction,
        inputs.temperature_response,
    }) |values| {
        if (values.len != complexes) return error.InvalidSolidDecompositionDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSolidDecompositionInput;
    }
    inline for (.{
        inputs.nitrogen_limitation_fraction, inputs.phosphorus_limitation_fraction,
    }) |values| for (values) |value| if (value == 0)
        return error.InvalidSolidDecompositionInput;
    inline for (.{
        inputs.negligible_substrate_carbon_g_c, inputs.biochemical_time_fraction_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSolidDecompositionInput;
}

fn fixture() Inputs {
    return .{
        .substrate_class_count = 5,
        .salt_chemistry_enabled = true,
        .substrate_carbon_g_c = &.{ 10, 10, 10, 10, 10 },
        .substrate_nitrogen_g_n = &.{ 2, 2, 2, 2, 2 },
        .substrate_phosphorus_g_p = &.{ 1, 1, 1, 1, 1 },
        .colonizing_active_biomass_g_c = &.{ 2, 2, 2, 2, 2 },
        .specific_decomposition_rate_per_g_c_h = &.{ 0.1, 0.1, 0.1, 0.1, 0.1 },
        .colonized_soil_organic_carbon_g_c = &.{10},
        .microbial_activity_respiration_g_c = &.{4},
        .substrate_density_response = &.{0.5},
        .dissolved_carbon_product_response = &.{0.5},
        .nitrogen_limitation_fraction = &.{0.5},
        .phosphorus_limitation_fraction = &.{0.5},
        .temperature_response = &.{1},
        .negligible_substrate_carbon_g_c = 1e-12,
        .biochemical_time_fraction_h = 1,
    };
}

test "solid substrate decomposition preserves C N P equations" {
    var state = try State.init(std.testing.allocator, 1, 5);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(0.2, state.nitrogen_to_carbon_ratio_g_n_per_g_c[0], 1e-12);
    try std.testing.expect(state.decomposed_carbon_g_c[0] > 0);
    try std.testing.expectApproxEqAbs(
        state.decomposed_carbon_g_c[0] * 0.2 / 0.5,
        state.decomposed_nitrogen_g_n[0],
        1e-12,
    );
}

test "salt lignin acidity uses exact source class gate" {
    var state = try State.init(std.testing.allocator, 1, 5);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expect(sourceContributesToLigninAcidity(2, 3));
    try std.testing.expect(!sourceContributesToLigninAcidity(3, 3));
    try std.testing.expectApproxEqAbs(
        state.decomposed_carbon_g_c[3],
        state.lignin_decomposition_for_acidity_g_c,
        1e-12,
    );
}

test "negligible substrate produces zero decomposition" {
    var state = try State.init(std.testing.allocator, 1, 5);
    defer state.deinit();
    var inputs = fixture();
    inputs.substrate_carbon_g_c = &.{ 0, 10, 10, 10, 10 };
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.decomposed_carbon_g_c[0]);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1, 5);
    defer state.deinit();
    state.decomposed_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.temperature_response = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidSolidDecompositionInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.decomposed_carbon_g_c[0]);
}

test "NITRO outer colonized threshold equality zeros every structural output" {
    var state = try State.init(std.testing.allocator, 1, 5);
    defer state.deinit();
    var inputs = fixture();
    inputs.colonized_soil_organic_carbon_g_c = &.{1e-12};
    try calculate(&state, inputs);
    inline for (.{
        state.nitrogen_to_carbon_ratio_g_n_per_g_c,
        state.phosphorus_to_carbon_ratio_g_p_per_g_c,
        state.decomposed_carbon_g_c,
        state.decomposed_nitrogen_g_n,
        state.decomposed_phosphorus_g_p,
    }) |values| for (values) |value| try std.testing.expectEqual(0, value);
    try std.testing.expectEqual(0, state.lignin_decomposition_for_acidity_g_c);
}

test "NITRO 3242-3298 late derived overflow preserves all outputs" {
    var state = try State.init(std.testing.allocator, 1, 5);
    defer state.deinit();
    state.decomposed_carbon_g_c[0] = 7;
    state.decomposed_phosphorus_g_p[4] = 11;
    state.lignin_decomposition_for_acidity_g_c = 13;
    var inputs = fixture();
    inputs.specific_decomposition_rate_per_g_c_h =
        &.{ 0.1, 0.1, 0.1, 0.1, std.math.floatMax(f64) };
    inputs.biochemical_time_fraction_h = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSolidDecompositionResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.decomposed_carbon_g_c[0]);
    try std.testing.expectEqual(11, state.decomposed_phosphorus_g_p[4]);
    try std.testing.expectEqual(13, state.lignin_decomposition_for_acidity_g_c);
}
