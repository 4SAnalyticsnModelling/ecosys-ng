const std = @import("std");

pub const Inputs = struct {
    temperature_response: []const f64,
    water_response: []const f64,
    decay_carbon_response: []const f64,
    carbon_recycling_fraction: []const f64,
    nitrogen_recycling_fraction: []const f64,
    phosphorus_recycling_fraction: []const f64,
    basal_decomposition_rate_per_h: []const f64,
    structural_carbon_g_c: []const f64,
    structural_nitrogen_g_n: []const f64,
    structural_phosphorus_g_p: []const f64,
    biochemical_time_fraction_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    component_count: usize,
    decomposition_fraction: []f64,
    decomposed_carbon_g_c: []f64,
    decomposed_nitrogen_g_n: []f64,
    decomposed_phosphorus_g_p: []f64,
    recycled_carbon_g_c: []f64,
    recycled_nitrogen_g_n: []f64,
    recycled_phosphorus_g_p: []f64,
    litterfall_carbon_g_c: []f64,
    litterfall_nitrogen_g_n: []f64,
    litterfall_phosphorus_g_p: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        unit_count: usize,
        component_count: usize,
    ) !State {
        if (unit_count == 0 or component_count == 0)
            return error.InvalidBiomassDecompositionDimensions;
        const item_count = std.math.mul(usize, unit_count, component_count) catch
            return error.InvalidBiomassDecompositionDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
        state.component_count = component_count;
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

/// Exact NITRO.F 2682--2704 structural microbial biomass decomposition.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const item_count = state.unit_count * state.component_count;
    const temporary = try state.allocator.alloc([10]f64, item_count);
    defer state.allocator.free(temporary);

    for (0..state.unit_count) |unit| {
        for (0..state.component_count) |component| {
            const item = unit * state.component_count + component;
            const decomposition_fraction = @sqrt(inputs.temperature_response[unit]) *
                inputs.water_response[unit] *
                inputs.basal_decomposition_rate_per_h[item] *
                inputs.decay_carbon_response[unit] *
                inputs.biochemical_time_fraction_h;
            const decomposed_carbon = @max(
                0,
                inputs.structural_carbon_g_c[item] * decomposition_fraction,
            );
            const decomposed_nitrogen = @max(
                0,
                inputs.structural_nitrogen_g_n[item] * decomposition_fraction,
            );
            const decomposed_phosphorus = @max(
                0,
                inputs.structural_phosphorus_g_p[item] * decomposition_fraction,
            );
            const recycled_carbon =
                decomposed_carbon * inputs.carbon_recycling_fraction[unit];
            const recycled_nitrogen = decomposed_nitrogen *
                (inputs.nitrogen_recycling_fraction[unit] +
                    (1 - inputs.nitrogen_recycling_fraction[unit]) *
                        inputs.carbon_recycling_fraction[unit]);
            const recycled_phosphorus = decomposed_phosphorus *
                (inputs.phosphorus_recycling_fraction[unit] +
                    (1 - inputs.phosphorus_recycling_fraction[unit]) *
                        inputs.carbon_recycling_fraction[unit]);
            temporary[item] = .{
                decomposition_fraction,
                decomposed_carbon,
                decomposed_nitrogen,
                decomposed_phosphorus,
                recycled_carbon,
                recycled_nitrogen,
                recycled_phosphorus,
                decomposed_carbon - recycled_carbon,
                decomposed_nitrogen - recycled_nitrogen,
                decomposed_phosphorus - recycled_phosphorus,
            };
            for (temporary[item]) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteBiomassDecompositionResult;
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
    const item_count = std.math.mul(usize, state.unit_count, state.component_count) catch
        return error.InvalidBiomassDecompositionDimensions;
    inline for (.{
        inputs.temperature_response,
        inputs.water_response,
        inputs.decay_carbon_response,
        inputs.carbon_recycling_fraction,
        inputs.nitrogen_recycling_fraction,
        inputs.phosphorus_recycling_fraction,
    }) |values| {
        if (values.len != state.unit_count)
            return error.InvalidBiomassDecompositionDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidBiomassDecompositionInput;
    }
    inline for (.{
        inputs.carbon_recycling_fraction,
        inputs.nitrogen_recycling_fraction,
        inputs.phosphorus_recycling_fraction,
    }) |values| for (values) |value| if (value > 1)
        return error.InvalidBiomassDecompositionInput;
    inline for (.{
        inputs.basal_decomposition_rate_per_h,
        inputs.structural_carbon_g_c,
        inputs.structural_nitrogen_g_n,
        inputs.structural_phosphorus_g_p,
    }) |values| {
        if (values.len != item_count)
            return error.InvalidBiomassDecompositionDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidBiomassDecompositionInput;
    }
    if (!std.math.isFinite(inputs.biochemical_time_fraction_h) or
        inputs.biochemical_time_fraction_h < 0)
        return error.InvalidBiomassDecompositionInput;
}

fn fixture() Inputs {
    return .{
        .temperature_response = &.{0.64},
        .water_response = &.{0.5},
        .decay_carbon_response = &.{1},
        .carbon_recycling_fraction = &.{0.2},
        .nitrogen_recycling_fraction = &.{0.5},
        .phosphorus_recycling_fraction = &.{0.25},
        .basal_decomposition_rate_per_h = &.{ 0.1, 0.05 },
        .structural_carbon_g_c = &.{ 10, 20 },
        .structural_nitrogen_g_n = &.{ 2, 2 },
        .structural_phosphorus_g_p = &.{ 1, 1 },
        .biochemical_time_fraction_h = 1,
    };
}

test "NITRO 2682-2704 preserves component decomposition order" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(0.04, state.decomposition_fraction[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.4, state.decomposed_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.08, state.recycled_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.32, state.litterfall_carbon_g_c[0], 1e-12);
}

test "NITRO 2682-2704 element partitions close exactly" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    for (0..2) |item| {
        try std.testing.expectApproxEqAbs(
            state.decomposed_carbon_g_c[item],
            state.recycled_carbon_g_c[item] + state.litterfall_carbon_g_c[item],
            1e-12,
        );
        try std.testing.expectApproxEqAbs(
            state.decomposed_nitrogen_g_n[item],
            state.recycled_nitrogen_g_n[item] + state.litterfall_nitrogen_g_n[item],
            1e-12,
        );
        try std.testing.expectApproxEqAbs(
            state.decomposed_phosphorus_g_p[item],
            state.recycled_phosphorus_g_p[item] + state.litterfall_phosphorus_g_p[item],
            1e-12,
        );
    }
}

test "NITRO 2682-2704 derived overflow preserves decomposition state" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.decomposed_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.temperature_response = &.{std.math.floatMax(f64)};
    inputs.water_response = &.{std.math.floatMax(f64)};
    try std.testing.expectError(
        error.NonFiniteBiomassDecompositionResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.decomposed_carbon_g_c[0]);
}
