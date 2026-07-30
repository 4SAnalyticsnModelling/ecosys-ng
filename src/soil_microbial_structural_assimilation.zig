const std = @import("std");

pub const Inputs = struct {
    nonstructural_carbon_g_c: []const f64,
    nonstructural_nitrogen_g_n: []const f64,
    nonstructural_phosphorus_g_p: []const f64,
    biological_activity_fraction: []const f64,
    structural_partition_fraction: []const f64,
    maximum_nitrogen_to_carbon_ratio_g_n_per_g_c: []const f64,
    maximum_phosphorus_to_carbon_ratio_g_p_per_g_c: []const f64,
    nonstructural_to_structural_rate_per_h: f64,
    biochemical_time_fraction_h: f64,
    negligible_nonstructural_carbon_g_c: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    component_count: usize,
    total_carbon_assimilation_g_c: []f64,
    component_carbon_assimilation_g_c: []f64,
    component_nitrogen_assimilation_g_n: []f64,
    component_phosphorus_assimilation_g_p: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        unit_count: usize,
        component_count: usize,
    ) !State {
        if (unit_count == 0 or component_count == 0)
            return error.InvalidStructuralAssimilationDimensions;
        const item_count = std.math.mul(usize, unit_count, component_count) catch
            return error.InvalidStructuralAssimilationDimensions;
        const total_carbon = try allocator.alloc(f64, unit_count);
        errdefer allocator.free(total_carbon);
        const carbon = try allocator.alloc(f64, item_count);
        errdefer allocator.free(carbon);
        const nitrogen = try allocator.alloc(f64, item_count);
        errdefer allocator.free(nitrogen);
        const phosphorus = try allocator.alloc(f64, item_count);
        @memset(total_carbon, 0);
        @memset(carbon, 0);
        @memset(nitrogen, 0);
        @memset(phosphorus, 0);
        return .{
            .allocator = allocator,
            .unit_count = unit_count,
            .component_count = component_count,
            .total_carbon_assimilation_g_c = total_carbon,
            .component_carbon_assimilation_g_c = carbon,
            .component_nitrogen_assimilation_g_n = nitrogen,
            .component_phosphorus_assimilation_g_p = phosphorus,
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.total_carbon_assimilation_g_c);
        self.allocator.free(self.component_carbon_assimilation_g_c);
        self.allocator.free(self.component_nitrogen_assimilation_g_n);
        self.allocator.free(self.component_phosphorus_assimilation_g_p);
        self.* = undefined;
    }
};

/// Exact NITRO.F 2653--2680 nonstructural-to-structural C/N/P assimilation.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const item_count = state.unit_count * state.component_count;
    const total_carbon = try state.allocator.alloc(f64, state.unit_count);
    defer state.allocator.free(total_carbon);
    const temporary = try state.allocator.alloc([3]f64, item_count);
    defer state.allocator.free(temporary);

    for (0..state.unit_count) |unit| {
        const carbon = inputs.nonstructural_carbon_g_c[unit];
        total_carbon[unit] = inputs.biological_activity_fraction[unit] *
            inputs.nonstructural_to_structural_rate_per_h *
            @max(0, carbon) * inputs.biochemical_time_fraction_h;
        if (!std.math.isFinite(total_carbon[unit]) or total_carbon[unit] < 0)
            return error.NonFiniteStructuralAssimilationResult;
        for (0..state.component_count) |component| {
            const item = unit * state.component_count + component;
            const component_carbon =
                inputs.structural_partition_fraction[item] * total_carbon[unit];
            var component_nitrogen: f64 = 0;
            var component_phosphorus: f64 = 0;
            if (carbon > inputs.negligible_nonstructural_carbon_g_c) {
                component_nitrogen = @min(
                    inputs.structural_partition_fraction[item] *
                        @max(0, inputs.nonstructural_nitrogen_g_n[unit]),
                    component_carbon * @min(
                        inputs.maximum_nitrogen_to_carbon_ratio_g_n_per_g_c[item],
                        inputs.nonstructural_nitrogen_g_n[unit] / carbon,
                    ),
                );
                component_phosphorus = @min(
                    inputs.structural_partition_fraction[item] *
                        @max(0, inputs.nonstructural_phosphorus_g_p[unit]),
                    component_carbon * @min(
                        inputs.maximum_phosphorus_to_carbon_ratio_g_p_per_g_c[item],
                        inputs.nonstructural_phosphorus_g_p[unit] / carbon,
                    ),
                );
            }
            temporary[item] = .{
                component_carbon,
                component_nitrogen,
                component_phosphorus,
            };
            for (temporary[item]) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteStructuralAssimilationResult;
        }
    }

    @memcpy(state.total_carbon_assimilation_g_c, total_carbon);
    for (temporary, 0..) |values, item| {
        state.component_carbon_assimilation_g_c[item] = values[0];
        state.component_nitrogen_assimilation_g_n[item] = values[1];
        state.component_phosphorus_assimilation_g_p[item] = values[2];
    }
}

fn validate(state: *const State, inputs: Inputs) !void {
    const item_count = std.math.mul(usize, state.unit_count, state.component_count) catch
        return error.InvalidStructuralAssimilationDimensions;
    inline for (.{
        inputs.nonstructural_carbon_g_c,
        inputs.nonstructural_nitrogen_g_n,
        inputs.nonstructural_phosphorus_g_p,
        inputs.biological_activity_fraction,
    }) |values| {
        if (values.len != state.unit_count)
            return error.InvalidStructuralAssimilationDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidStructuralAssimilationInput;
    }
    inline for (.{
        inputs.structural_partition_fraction,
        inputs.maximum_nitrogen_to_carbon_ratio_g_n_per_g_c,
        inputs.maximum_phosphorus_to_carbon_ratio_g_p_per_g_c,
    }) |values| {
        if (values.len != item_count)
            return error.InvalidStructuralAssimilationDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidStructuralAssimilationInput;
    }
    for (inputs.biological_activity_fraction) |value| if (value < 0)
        return error.InvalidStructuralAssimilationInput;
    inline for (.{
        inputs.nonstructural_to_structural_rate_per_h,
        inputs.biochemical_time_fraction_h,
        inputs.negligible_nonstructural_carbon_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidStructuralAssimilationInput;
}

fn fixture() Inputs {
    return .{
        .nonstructural_carbon_g_c = &.{10},
        .nonstructural_nitrogen_g_n = &.{2},
        .nonstructural_phosphorus_g_p = &.{1},
        .biological_activity_fraction = &.{0.5},
        .structural_partition_fraction = &.{ 0.6, 0.4 },
        .maximum_nitrogen_to_carbon_ratio_g_n_per_g_c = &.{ 0.15, 0.1 },
        .maximum_phosphorus_to_carbon_ratio_g_p_per_g_c = &.{ 0.08, 0.05 },
        .nonstructural_to_structural_rate_per_h = 0.2,
        .biochemical_time_fraction_h = 1,
        .negligible_nonstructural_carbon_g_c = 1e-12,
    };
}

test "NITRO 2653-2680 preserves component assimilation caps" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(1, state.total_carbon_assimilation_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.6, state.component_carbon_assimilation_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.09, state.component_nitrogen_assimilation_g_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.048, state.component_phosphorus_assimilation_g_p[0], 1e-12);
}

test "NITRO 2653-2680 strict carbon threshold zeros nutrient assimilation" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.nonstructural_carbon_g_c = &.{1e-12};
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.component_nitrogen_assimilation_g_n[0]);
    try std.testing.expectEqual(0, state.component_phosphorus_assimilation_g_p[1]);
}

test "NITRO 2653-2680 derived overflow preserves assimilation state" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.total_carbon_assimilation_g_c[0] = 7;
    var inputs = fixture();
    inputs.biological_activity_fraction = &.{std.math.floatMax(f64)};
    inputs.nonstructural_to_structural_rate_per_h = 2;
    try std.testing.expectError(
        error.NonFiniteStructuralAssimilationResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.total_carbon_assimilation_g_c[0]);
}
