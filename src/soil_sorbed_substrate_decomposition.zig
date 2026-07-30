const std = @import("std");

pub const Inputs = struct {
    sorbed_carbon_g_c: []const f64,
    sorbed_nitrogen_g_n: []const f64,
    sorbed_phosphorus_g_p: []const f64,
    sorbed_acetate_g_c: []const f64,
    colonized_soil_organic_carbon_g_c: []const f64,
    microbial_activity_respiration_g_c: []const f64,
    substrate_density_response: []const f64,
    dissolved_carbon_product_response: []const f64,
    nitrogen_limitation_fraction: []const f64,
    phosphorus_limitation_fraction: []const f64,
    temperature_response: []const f64,
    sorbed_carbon_decomposition_rate_h: f64,
    sorbed_acetate_decomposition_rate_h: f64,
    negligible_carbon_g_c: f64,
    biochemical_time_fraction_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    complex_count: usize,
    nitrogen_to_carbon_ratio_g_n_per_g_c: []f64,
    phosphorus_to_carbon_ratio_g_p_per_g_c: []f64,
    decomposed_carbon_g_c: []f64,
    decomposed_nitrogen_g_n: []f64,
    decomposed_phosphorus_g_p: []f64,
    decomposed_acetate_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, complex_count: usize) !State {
        if (complex_count == 0) return error.InvalidSorbedDecompositionDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.complex_count = complex_count;
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
            @field(state, field.name) = try allocator.alloc(f64, complex_count);
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

/// Exact NITRO.F 3399--3458 sorbed C/N/P and acetate decomposition.
/// Complex traversal and floating-point operation order follow the source.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([6]f64, state.complex_count);
    defer state.allocator.free(temporary);
    for (0..state.complex_count) |complex| {
        if (inputs.colonized_soil_organic_carbon_g_c[complex] <=
            inputs.negligible_carbon_g_c or
            inputs.sorbed_carbon_g_c[complex] <= inputs.negligible_carbon_g_c)
        {
            temporary[complex] = @splat(0);
            continue;
        }
        const carbon = inputs.sorbed_carbon_g_c[complex];
        const nitrogen_ratio = @max(0, inputs.sorbed_nitrogen_g_n[complex] / carbon);
        const phosphorus_ratio = @max(0, inputs.sorbed_phosphorus_g_p[complex] / carbon);
        const colonized = inputs.colonized_soil_organic_carbon_g_c[complex];
        const activity_limited_carbon = inputs.sorbed_carbon_decomposition_rate_h *
            inputs.microbial_activity_respiration_g_c[complex] *
            inputs.substrate_density_response[complex] *
            inputs.dissolved_carbon_product_response[complex] *
            inputs.temperature_response[complex] *
            carbon / colonized;
        const decomposed_carbon = @max(0, @min(
            carbon * inputs.biochemical_time_fraction_h,
            activity_limited_carbon,
        ));
        const decomposed_nitrogen = @max(0, @min(
            inputs.sorbed_nitrogen_g_n[complex] * inputs.biochemical_time_fraction_h,
            nitrogen_ratio * decomposed_carbon,
        )) / inputs.nitrogen_limitation_fraction[complex];
        const decomposed_phosphorus = @max(0, @min(
            inputs.sorbed_phosphorus_g_p[complex] * inputs.biochemical_time_fraction_h,
            phosphorus_ratio * decomposed_carbon,
        )) / inputs.phosphorus_limitation_fraction[complex];
        // Source intentionally excludes OQCI from sorbed acetate decomposition.
        const activity_limited_acetate = inputs.sorbed_acetate_decomposition_rate_h *
            inputs.microbial_activity_respiration_g_c[complex] *
            inputs.substrate_density_response[complex] *
            inputs.temperature_response[complex] *
            inputs.sorbed_acetate_g_c[complex] / colonized;
        const decomposed_acetate = @max(0, @min(
            inputs.sorbed_acetate_g_c[complex] * inputs.biochemical_time_fraction_h,
            activity_limited_acetate,
        ));
        temporary[complex] = .{
            nitrogen_ratio,      phosphorus_ratio,      decomposed_carbon,
            decomposed_nitrogen, decomposed_phosphorus, decomposed_acetate,
        };
        for (temporary[complex]) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteSorbedDecompositionResult;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, complex| @field(state, field.name)[complex] = values[index];
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
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != state.complex_count) return error.InvalidSorbedDecompositionDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSorbedDecompositionInput;
    };
    inline for (.{
        inputs.nitrogen_limitation_fraction, inputs.phosphorus_limitation_fraction,
    }) |values| for (values) |value| if (value == 0)
        return error.InvalidSorbedDecompositionInput;
    inline for (.{
        inputs.sorbed_carbon_decomposition_rate_h,
        inputs.sorbed_acetate_decomposition_rate_h,
        inputs.negligible_carbon_g_c,
        inputs.biochemical_time_fraction_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSorbedDecompositionInput;
}

fn fixture() Inputs {
    return .{
        .sorbed_carbon_g_c = &.{10},
        .sorbed_nitrogen_g_n = &.{2},
        .sorbed_phosphorus_g_p = &.{1},
        .sorbed_acetate_g_c = &.{4},
        .colonized_soil_organic_carbon_g_c = &.{10},
        .microbial_activity_respiration_g_c = &.{4},
        .substrate_density_response = &.{0.5},
        .dissolved_carbon_product_response = &.{0.25},
        .nitrogen_limitation_fraction = &.{0.5},
        .phosphorus_limitation_fraction = &.{0.5},
        .temperature_response = &.{1},
        .sorbed_carbon_decomposition_rate_h = 0.1,
        .sorbed_acetate_decomposition_rate_h = 0.2,
        .negligible_carbon_g_c = 1e-12,
        .biochemical_time_fraction_h = 1,
    };
}

test "sorbed C N and P decomposition preserves source controls" {
    var state = try State.init(std.testing.allocator, 1);
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

test "sorbed acetate decomposition omits DOC product inhibition" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(0.16, state.decomposed_acetate_g_c[0], 1e-12);
}

test "zero sorbed carbon gates acetate decomposition too" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.sorbed_carbon_g_c = &.{0};
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.decomposed_acetate_g_c[0]);
}

test "NITRO outer and inner threshold equality zero every sorbed output" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.sorbed_carbon_g_c = &.{ 10, 1e-12 };
    inputs.sorbed_nitrogen_g_n = &.{ 2, 2 };
    inputs.sorbed_phosphorus_g_p = &.{ 1, 1 };
    inputs.sorbed_acetate_g_c = &.{ 4, 4 };
    inputs.colonized_soil_organic_carbon_g_c = &.{ 1e-12, 10 };
    inputs.microbial_activity_respiration_g_c = &.{ 4, 4 };
    inputs.substrate_density_response = &.{ 0.5, 0.5 };
    inputs.dissolved_carbon_product_response = &.{ 0.25, 0.25 };
    inputs.nitrogen_limitation_fraction = &.{ 0.5, 0.5 };
    inputs.phosphorus_limitation_fraction = &.{ 0.5, 0.5 };
    inputs.temperature_response = &.{ 1, 1 };
    try calculate(&state, inputs);
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64)
        for (@field(state, field.name)) |value| try std.testing.expectEqual(0, value);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.decomposed_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.temperature_response = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidSorbedDecompositionInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.decomposed_carbon_g_c[0]);
}

test "NITRO 3399-3458 late derived overflow preserves all sorbed outputs" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.decomposed_carbon_g_c[0] = 7;
    state.decomposed_acetate_g_c[0] = 11;
    var inputs = fixture();
    inputs.sorbed_acetate_decomposition_rate_h = std.math.floatMax(f64);
    inputs.biochemical_time_fraction_h = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSorbedDecompositionResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.decomposed_carbon_g_c[0]);
    try std.testing.expectEqual(11, state.decomposed_acetate_g_c[0]);
}
