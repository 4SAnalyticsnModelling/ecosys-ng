const std = @import("std");

pub const Inputs = struct {
    substrate_class_count: usize,
    substrate_carbon_g_c: []const f64,
    colonized_substrate_carbon_g_c: []const f64,
    colonization_rate_per_g_respired_carbon: []const f64,
    microbial_activity_respiration_g_c: []const f64,
    negligible_substrate_carbon_g_c: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    complex_count: usize,
    item_count: usize,
    total_substrate_carbon_g_c: []f64,
    total_colonized_carbon_before_growth_g_c: []f64,
    activity_driven_colonization_g_c: []f64,
    colonized_substrate_carbon_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, complex_count: usize, class_count: usize) !State {
        if (complex_count == 0 or class_count == 0)
            return error.InvalidLitterColonizationDimensions;
        const items = std.math.mul(usize, complex_count, class_count) catch
            return error.InvalidLitterColonizationDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.complex_count = complex_count;
        state.item_count = items;
        state.total_substrate_carbon_g_c = try allocator.alloc(f64, complex_count);
        errdefer allocator.free(state.total_substrate_carbon_g_c);
        state.total_colonized_carbon_before_growth_g_c = try allocator.alloc(f64, complex_count);
        errdefer allocator.free(state.total_colonized_carbon_before_growth_g_c);
        state.activity_driven_colonization_g_c = try allocator.alloc(f64, complex_count);
        errdefer allocator.free(state.activity_driven_colonization_g_c);
        state.colonized_substrate_carbon_g_c = try allocator.alloc(f64, items);
        @memset(state.total_substrate_carbon_g_c, 0);
        @memset(state.total_colonized_carbon_before_growth_g_c, 0);
        @memset(state.activity_driven_colonization_g_c, 0);
        @memset(state.colonized_substrate_carbon_g_c, 0);
        return state;
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.total_substrate_carbon_g_c);
        self.allocator.free(self.total_colonized_carbon_before_growth_g_c);
        self.allocator.free(self.activity_driven_colonization_g_c);
        self.allocator.free(self.colonized_substrate_carbon_g_c);
        self.* = undefined;
    }
};

/// Exact NITRO.F 3860--3902 activity-driven colonization of new litter.
/// Traversal follows K, then structural fraction M in both source passes.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const totals = try state.allocator.alloc([3]f64, state.complex_count);
    defer state.allocator.free(totals);
    const colonized = try state.allocator.alloc(f64, state.item_count);
    defer state.allocator.free(colonized);
    for (0..state.complex_count) |complex| {
        var total_carbon: f64 = 0;
        var total_colonized: f64 = 0;
        for (0..inputs.substrate_class_count) |class| {
            const item = complex * inputs.substrate_class_count + class;
            total_carbon += inputs.substrate_carbon_g_c[item];
            total_colonized += inputs.colonized_substrate_carbon_g_c[item];
        }
        const active = total_carbon > inputs.negligible_substrate_carbon_g_c;
        const colonization = if (active)
            inputs.colonization_rate_per_g_respired_carbon[complex] *
                @max(0, inputs.microbial_activity_respiration_g_c[complex])
        else
            0;
        totals[complex] = .{ total_carbon, total_colonized, colonization };
        for (0..inputs.substrate_class_count) |class| {
            const item = complex * inputs.substrate_class_count + class;
            colonized[item] = if (active)
                @min(inputs.substrate_carbon_g_c[item], inputs.colonized_substrate_carbon_g_c[item] +
                    colonization * inputs.substrate_carbon_g_c[item] / total_carbon)
            else
                @min(inputs.substrate_carbon_g_c[item], inputs.colonized_substrate_carbon_g_c[item]);
        }
    }
    for (totals) |values| inline for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLitterColonizationResult;
    for (colonized) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidLitterColonizationResult;
    for (totals, 0..) |values, complex| {
        state.total_substrate_carbon_g_c[complex] = values[0];
        state.total_colonized_carbon_before_growth_g_c[complex] = values[1];
        state.activity_driven_colonization_g_c[complex] = values[2];
    }
    @memcpy(state.colonized_substrate_carbon_g_c, colonized);
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.substrate_class_count == 0 or
        state.item_count != state.complex_count * inputs.substrate_class_count)
        return error.InvalidLitterColonizationDimensions;
    inline for (.{
        inputs.substrate_carbon_g_c, inputs.colonized_substrate_carbon_g_c,
    }) |values| {
        if (values.len != state.item_count) return error.InvalidLitterColonizationDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLitterColonizationInput;
    }
    if (inputs.colonization_rate_per_g_respired_carbon.len != state.complex_count or
        inputs.microbial_activity_respiration_g_c.len != state.complex_count)
        return error.InvalidLitterColonizationDimensions;
    for (inputs.colonization_rate_per_g_respired_carbon) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLitterColonizationInput;
    for (inputs.microbial_activity_respiration_g_c) |value|
        if (!std.math.isFinite(value))
            return error.InvalidLitterColonizationInput;
    if (!std.math.isFinite(inputs.negligible_substrate_carbon_g_c) or
        inputs.negligible_substrate_carbon_g_c < 0)
        return error.InvalidLitterColonizationInput;
}

fn fixture() Inputs {
    return .{
        .substrate_class_count = 3,
        .substrate_carbon_g_c = &.{ 10, 20, 30 },
        .colonized_substrate_carbon_g_c = &.{ 1, 2, 3 },
        .colonization_rate_per_g_respired_carbon = &.{0.5},
        .microbial_activity_respiration_g_c = &.{6},
        .negligible_substrate_carbon_g_c = 1e-12,
    };
}

test "activity colonizes substrate classes proportional to carbon stock" {
    var state = try State.init(std.testing.allocator, 1, 3);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectEqual(60, state.total_substrate_carbon_g_c[0]);
    try std.testing.expectEqual(6, state.total_colonized_carbon_before_growth_g_c[0]);
    try std.testing.expectEqual(3, state.activity_driven_colonization_g_c[0]);
    try std.testing.expectApproxEqAbs(1.5, state.colonized_substrate_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(4.5, state.colonized_substrate_carbon_g_c[2], 1e-12);
}

test "colonization cannot exceed substrate carbon" {
    var state = try State.init(std.testing.allocator, 1, 3);
    defer state.deinit();
    var inputs = fixture();
    inputs.colonization_rate_per_g_respired_carbon = &.{100};
    try calculate(&state, inputs);
    try std.testing.expectEqual(10, state.colonized_substrate_carbon_g_c[0]);
    try std.testing.expectEqual(30, state.colonized_substrate_carbon_g_c[2]);
}

test "zero total carbon only clamps existing colonization" {
    var state = try State.init(std.testing.allocator, 1, 3);
    defer state.deinit();
    var inputs = fixture();
    inputs.substrate_carbon_g_c = &.{ 0, 0, 0 };
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.colonized_substrate_carbon_g_c[0]);
    try std.testing.expectEqual(0, state.activity_driven_colonization_g_c[0]);
}

test "threshold equality skips DOSAK and clamps every structural fraction" {
    var state = try State.init(std.testing.allocator, 1, 3);
    defer state.deinit();
    var inputs = fixture();
    inputs.substrate_carbon_g_c = &.{ 1e-12, 0, 0 };
    inputs.colonized_substrate_carbon_g_c = &.{ 2e-12, 0, 0 };
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.activity_driven_colonization_g_c[0]);
    try std.testing.expectEqual(1e-12, state.colonized_substrate_carbon_g_c[0]);
}

test "negative microbial activity is zero floored by source AMAX" {
    var state = try State.init(std.testing.allocator, 1, 3);
    defer state.deinit();
    var inputs = fixture();
    inputs.microbial_activity_respiration_g_c = &.{-6};
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.activity_driven_colonization_g_c[0]);
    try std.testing.expectEqual(1, state.colonized_substrate_carbon_g_c[0]);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1, 3);
    defer state.deinit();
    state.colonized_substrate_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.microbial_activity_respiration_g_c = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidLitterColonizationInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.colonized_substrate_carbon_g_c[0]);
}

test "NITRO 3860-3902 late derived overflow leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1, 3);
    defer state.deinit();
    state.colonized_substrate_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.colonization_rate_per_g_respired_carbon =
        &.{std.math.floatMax(f64)};
    inputs.microbial_activity_respiration_g_c =
        &.{std.math.floatMax(f64)};
    try std.testing.expectError(
        error.InvalidLitterColonizationResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.colonized_substrate_carbon_g_c[0]);
}
