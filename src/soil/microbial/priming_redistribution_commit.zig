const std = @import("std");

pub const Inputs = struct {
    complex_count: usize,
    population_count: usize,
    biomass_fraction_count: usize,
    microbial_activity_respiration_g_c: []const f64,
    activity_transfer_g_c: []const f64,
    dissolved_organic_carbon_g_c: []const f64,
    dissolved_carbon_transfer_g_c: []const f64,
    dissolved_organic_nitrogen_g_n: []const f64,
    dissolved_nitrogen_transfer_g_n: []const f64,
    dissolved_organic_phosphorus_g_p: []const f64,
    dissolved_phosphorus_transfer_g_p: []const f64,
    acetate_g_c: []const f64,
    acetate_transfer_g_c: []const f64,
    microbial_carbon_g_c: []const f64,
    microbial_carbon_transfer_g_c: []const f64,
    microbial_nitrogen_g_n: []const f64,
    microbial_nitrogen_transfer_g_n: []const f64,
    microbial_phosphorus_g_p: []const f64,
    microbial_phosphorus_transfer_g_p: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    backing: []f64,
    complex_count: usize,
    microbial_item_count: usize,
    microbial_activity_respiration_g_c: []f64,
    dissolved_organic_carbon_g_c: []f64,
    dissolved_organic_nitrogen_g_n: []f64,
    dissolved_organic_phosphorus_g_p: []f64,
    acetate_g_c: []f64,
    microbial_carbon_g_c: []f64,
    microbial_nitrogen_g_n: []f64,
    microbial_phosphorus_g_p: []f64,
    total_microbial_activity_respiration_g_c: f64,

    pub fn init(
        allocator: std.mem.Allocator,
        complex_count: usize,
        population_count: usize,
        biomass_fraction_count: usize,
    ) !State {
        if (complex_count == 0 or population_count == 0 or biomass_fraction_count == 0)
            return error.InvalidPrimingCommitDimensions;
        const cp = std.math.mul(usize, complex_count, population_count) catch
            return error.InvalidPrimingCommitDimensions;
        const items = std.math.mul(usize, cp, biomass_fraction_count) catch
            return error.InvalidPrimingCommitDimensions;
        const complex_values = std.math.mul(usize, complex_count, 5) catch
            return error.InvalidPrimingCommitDimensions;
        const item_values = std.math.mul(usize, items, 3) catch
            return error.InvalidPrimingCommitDimensions;
        const count = std.math.add(usize, complex_values, item_values) catch
            return error.InvalidPrimingCommitDimensions;
        const backing = try allocator.alloc(f64, count);
        @memset(backing, 0);
        var offset: usize = 0;
        const activity = take(backing, &offset, complex_count);
        const carbon = take(backing, &offset, complex_count);
        const nitrogen = take(backing, &offset, complex_count);
        const phosphorus = take(backing, &offset, complex_count);
        const acetate = take(backing, &offset, complex_count);
        const microbial_carbon = take(backing, &offset, items);
        const microbial_nitrogen = take(backing, &offset, items);
        const microbial_phosphorus = take(backing, &offset, items);
        return .{
            .allocator = allocator,
            .backing = backing,
            .complex_count = complex_count,
            .microbial_item_count = items,
            .microbial_activity_respiration_g_c = activity,
            .dissolved_organic_carbon_g_c = carbon,
            .dissolved_organic_nitrogen_g_n = nitrogen,
            .dissolved_organic_phosphorus_g_p = phosphorus,
            .acetate_g_c = acetate,
            .microbial_carbon_g_c = microbial_carbon,
            .microbial_nitrogen_g_n = microbial_nitrogen,
            .microbial_phosphorus_g_p = microbial_phosphorus,
            .total_microbial_activity_respiration_g_c = 0,
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.backing);
        self.* = undefined;
    }
};

/// Exact NITRO.F 3140--3171 priming-transfer commit and activity aggregation.
/// Publication staging follows K, then N, then M source traversal.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc(f64, state.backing.len);
    defer state.allocator.free(temporary);
    var offset: usize = 0;
    const activity = take(temporary, &offset, state.complex_count);
    const carbon = take(temporary, &offset, state.complex_count);
    const nitrogen = take(temporary, &offset, state.complex_count);
    const phosphorus = take(temporary, &offset, state.complex_count);
    const acetate = take(temporary, &offset, state.complex_count);
    const microbial_carbon = take(temporary, &offset, state.microbial_item_count);
    const microbial_nitrogen = take(temporary, &offset, state.microbial_item_count);
    const microbial_phosphorus = take(temporary, &offset, state.microbial_item_count);

    // TOQCK is reset before the K loop and accumulated immediately after ROQCK.
    var total_activity: f64 = 0;
    for (0..state.complex_count) |complex| {
        activity[complex] = try checkedPoolAdd(
            inputs.microbial_activity_respiration_g_c[complex],
            inputs.activity_transfer_g_c[complex],
        );
        total_activity = total_activity + activity[complex];
        if (!std.math.isFinite(total_activity)) return error.InvalidPrimingCommitResult;
        carbon[complex] = try checkedPoolAdd(
            inputs.dissolved_organic_carbon_g_c[complex],
            inputs.dissolved_carbon_transfer_g_c[complex],
        );
        nitrogen[complex] = try checkedPoolAdd(
            inputs.dissolved_organic_nitrogen_g_n[complex],
            inputs.dissolved_nitrogen_transfer_g_n[complex],
        );
        phosphorus[complex] = try checkedPoolAdd(
            inputs.dissolved_organic_phosphorus_g_p[complex],
            inputs.dissolved_phosphorus_transfer_g_p[complex],
        );
        acetate[complex] = try checkedPoolAdd(
            inputs.acetate_g_c[complex],
            inputs.acetate_transfer_g_c[complex],
        );
        for (0..inputs.population_count) |population| {
            for (0..inputs.biomass_fraction_count) |fraction| {
                const item = (complex * inputs.population_count + population) *
                    inputs.biomass_fraction_count + fraction;
                microbial_carbon[item] = try checkedPoolAdd(
                    inputs.microbial_carbon_g_c[item],
                    inputs.microbial_carbon_transfer_g_c[item],
                );
                microbial_nitrogen[item] = try checkedPoolAdd(
                    inputs.microbial_nitrogen_g_n[item],
                    inputs.microbial_nitrogen_transfer_g_n[item],
                );
                microbial_phosphorus[item] = try checkedPoolAdd(
                    inputs.microbial_phosphorus_g_p[item],
                    inputs.microbial_phosphorus_transfer_g_p[item],
                );
            }
        }
    }
    @memcpy(state.backing, temporary);
    state.total_microbial_activity_respiration_g_c = total_activity;
}

fn checkedPoolAdd(pool: f64, transfer: f64) !f64 {
    const result = pool + transfer;
    if (!std.math.isFinite(result) or result < 0)
        return error.InvalidPrimingCommitResult;
    return result;
}

fn take(backing: []f64, offset: *usize, count: usize) []f64 {
    const result = backing[offset.* .. offset.* + count];
    offset.* += count;
    return result;
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.complex_count != state.complex_count or inputs.population_count == 0 or
        inputs.biomass_fraction_count == 0)
        return error.InvalidPrimingCommitDimensions;
    const cp = std.math.mul(usize, inputs.complex_count, inputs.population_count) catch
        return error.InvalidPrimingCommitDimensions;
    const items = std.math.mul(usize, cp, inputs.biomass_fraction_count) catch
        return error.InvalidPrimingCommitDimensions;
    if (items != state.microbial_item_count) return error.InvalidPrimingCommitDimensions;
    inline for (.{
        inputs.microbial_activity_respiration_g_c, inputs.activity_transfer_g_c,
        inputs.dissolved_organic_carbon_g_c,       inputs.dissolved_carbon_transfer_g_c,
        inputs.dissolved_organic_nitrogen_g_n,     inputs.dissolved_nitrogen_transfer_g_n,
        inputs.dissolved_organic_phosphorus_g_p,   inputs.dissolved_phosphorus_transfer_g_p,
        inputs.acetate_g_c,                        inputs.acetate_transfer_g_c,
    }) |values| {
        if (values.len != state.complex_count) return error.InvalidPrimingCommitDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidPrimingCommitInput;
    }
    inline for (.{
        inputs.microbial_carbon_g_c,     inputs.microbial_carbon_transfer_g_c,
        inputs.microbial_nitrogen_g_n,   inputs.microbial_nitrogen_transfer_g_n,
        inputs.microbial_phosphorus_g_p, inputs.microbial_phosphorus_transfer_g_p,
    }) |values| {
        if (values.len != items) return error.InvalidPrimingCommitDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidPrimingCommitInput;
    }
}

fn fixture() Inputs {
    return .{
        .complex_count = 2,
        .population_count = 1,
        .biomass_fraction_count = 2,
        .microbial_activity_respiration_g_c = &.{ 4, 2 },
        .activity_transfer_g_c = &.{ -1, 1 },
        .dissolved_organic_carbon_g_c = &.{ 8, 4 },
        .dissolved_carbon_transfer_g_c = &.{ -2, 2 },
        .dissolved_organic_nitrogen_g_n = &.{ 2, 1 },
        .dissolved_nitrogen_transfer_g_n = &.{ -0.5, 0.5 },
        .dissolved_organic_phosphorus_g_p = &.{ 1, 0.5 },
        .dissolved_phosphorus_transfer_g_p = &.{ -0.25, 0.25 },
        .acetate_g_c = &.{ 3, 1 },
        .acetate_transfer_g_c = &.{ -0.5, 0.5 },
        .microbial_carbon_g_c = &.{ 4, 2, 2, 1 },
        .microbial_carbon_transfer_g_c = &.{ -1, -0.5, 1, 0.5 },
        .microbial_nitrogen_g_n = &.{ 2, 1, 1, 0.5 },
        .microbial_nitrogen_transfer_g_n = &.{ -0.5, -0.25, 0.5, 0.25 },
        .microbial_phosphorus_g_p = &.{ 1, 0.5, 0.5, 0.25 },
        .microbial_phosphorus_transfer_g_p = &.{ -0.25, -0.125, 0.25, 0.125 },
    };
}

test "commit applies all priming accumulators and totals activity" {
    var state = try State.init(std.testing.allocator, 2, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(3, state.microbial_activity_respiration_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(3, state.microbial_activity_respiration_g_c[1], 1e-12);
    try std.testing.expectApproxEqAbs(6, state.total_microbial_activity_respiration_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(6, state.dissolved_organic_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(3, state.microbial_carbon_g_c[0], 1e-12);
}

test "conservative deltas preserve every published pool total" {
    var state = try State.init(std.testing.allocator, 2, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    inline for (.{
        .{ state.microbial_activity_respiration_g_c, @as(f64, 6) },
        .{ state.dissolved_organic_carbon_g_c, @as(f64, 12) },
        .{ state.dissolved_organic_nitrogen_g_n, @as(f64, 3) },
        .{ state.dissolved_organic_phosphorus_g_p, @as(f64, 1.5) },
        .{ state.acetate_g_c, @as(f64, 4) },
        .{ state.microbial_carbon_g_c, @as(f64, 9) },
        .{ state.microbial_nitrogen_g_n, @as(f64, 4.5) },
        .{ state.microbial_phosphorus_g_p, @as(f64, 2.25) },
    }) |expected| {
        var total: f64 = 0;
        for (expected[0]) |value| total += value;
        try std.testing.expectApproxEqAbs(expected[1], total, 1e-12);
    }
}

test "negative proposed pool fails atomically" {
    var state = try State.init(std.testing.allocator, 2, 1, 2);
    defer state.deinit();
    state.dissolved_organic_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.dissolved_carbon_transfer_g_c = &.{ -9, 9 };
    try std.testing.expectError(error.InvalidPrimingCommitResult, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.dissolved_organic_carbon_g_c[0]);
}

test "non-finite transfer fails before publication" {
    var state = try State.init(std.testing.allocator, 2, 1, 2);
    defer state.deinit();
    state.microbial_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.microbial_carbon_transfer_g_c = &.{ std.math.nan(f64), 0, 0, 0 };
    try std.testing.expectError(error.InvalidPrimingCommitInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.microbial_carbon_g_c[0]);
}

test "late microbial overflow leaves the entire priming publication unchanged" {
    var state = try State.init(std.testing.allocator, 2, 1, 2);
    defer state.deinit();
    state.microbial_activity_respiration_g_c[0] = 7;
    state.microbial_phosphorus_g_p[3] = 11;
    state.total_microbial_activity_respiration_g_c = 13;
    var inputs = fixture();
    inputs.microbial_phosphorus_g_p =
        &.{ 1, 0.5, 0.5, std.math.floatMax(f64) };
    inputs.microbial_phosphorus_transfer_g_p =
        &.{ -0.25, -0.125, 0.25, std.math.floatMax(f64) };
    try std.testing.expectError(
        error.InvalidPrimingCommitResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.microbial_activity_respiration_g_c[0]);
    try std.testing.expectEqual(11, state.microbial_phosphorus_g_p[3]);
    try std.testing.expectEqual(13, state.total_microbial_activity_respiration_g_c);
}
