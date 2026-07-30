const std = @import("std");

pub const Inputs = struct {
    dissolved_organic_carbon_g_c: []const f64,
    dissolved_organic_nitrogen_g_n: []const f64,
    dissolved_organic_phosphorus_g_p: []const f64,
    dissolved_acetate_g_c: []const f64,
    sorbed_organic_carbon_g_c: []const f64,
    sorbed_organic_nitrogen_g_n: []const f64,
    sorbed_organic_phosphorus_g_p: []const f64,
    sorbed_acetate_g_c: []const f64,
    organic_carbon_sorption_g_c: []const f64,
    organic_nitrogen_sorption_g_n: []const f64,
    organic_phosphorus_sorption_g_p: []const f64,
    acetate_sorption_g_c: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    backing: []f64,
    complex_count: usize,
    dissolved_organic_carbon_g_c: []f64,
    dissolved_organic_nitrogen_g_n: []f64,
    dissolved_organic_phosphorus_g_p: []f64,
    dissolved_acetate_g_c: []f64,
    sorbed_organic_carbon_g_c: []f64,
    sorbed_organic_nitrogen_g_n: []f64,
    sorbed_organic_phosphorus_g_p: []f64,
    sorbed_acetate_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, complex_count: usize) !State {
        if (complex_count == 0) return error.InvalidSorptionCommitDimensions;
        const count = std.math.mul(usize, complex_count, 8) catch
            return error.InvalidSorptionCommitDimensions;
        const backing = try allocator.alloc(f64, count);
        @memset(backing, 0);
        var offset: usize = 0;
        return .{
            .allocator = allocator,
            .backing = backing,
            .complex_count = complex_count,
            .dissolved_organic_carbon_g_c = take(backing, &offset, complex_count),
            .dissolved_organic_nitrogen_g_n = take(backing, &offset, complex_count),
            .dissolved_organic_phosphorus_g_p = take(backing, &offset, complex_count),
            .dissolved_acetate_g_c = take(backing, &offset, complex_count),
            .sorbed_organic_carbon_g_c = take(backing, &offset, complex_count),
            .sorbed_organic_nitrogen_g_n = take(backing, &offset, complex_count),
            .sorbed_organic_phosphorus_g_p = take(backing, &offset, complex_count),
            .sorbed_acetate_g_c = take(backing, &offset, complex_count),
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.backing);
        self.* = undefined;
    }
};

/// Exact NITRO.F 3671--3685 signed dissolved/sorbed exchange commit.
/// Per K, all dissolved subtractions precede all sorbed additions.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([8]f64, state.complex_count);
    defer state.allocator.free(temporary);
    for (0..state.complex_count) |complex| {
        temporary[complex] = .{
            inputs.dissolved_organic_carbon_g_c[complex] -
                inputs.organic_carbon_sorption_g_c[complex],
            inputs.dissolved_organic_nitrogen_g_n[complex] -
                inputs.organic_nitrogen_sorption_g_n[complex],
            inputs.dissolved_organic_phosphorus_g_p[complex] -
                inputs.organic_phosphorus_sorption_g_p[complex],
            inputs.dissolved_acetate_g_c[complex] - inputs.acetate_sorption_g_c[complex],
            inputs.sorbed_organic_carbon_g_c[complex] +
                inputs.organic_carbon_sorption_g_c[complex],
            inputs.sorbed_organic_nitrogen_g_n[complex] +
                inputs.organic_nitrogen_sorption_g_n[complex],
            inputs.sorbed_organic_phosphorus_g_p[complex] +
                inputs.organic_phosphorus_sorption_g_p[complex],
            inputs.sorbed_acetate_g_c[complex] + inputs.acetate_sorption_g_c[complex],
        };
    }
    for (temporary) |values| inline for (values) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSorptionCommitResult;
    for (temporary, 0..) |values, complex| {
        state.dissolved_organic_carbon_g_c[complex] = values[0];
        state.dissolved_organic_nitrogen_g_n[complex] = values[1];
        state.dissolved_organic_phosphorus_g_p[complex] = values[2];
        state.dissolved_acetate_g_c[complex] = values[3];
        state.sorbed_organic_carbon_g_c[complex] = values[4];
        state.sorbed_organic_nitrogen_g_n[complex] = values[5];
        state.sorbed_organic_phosphorus_g_p[complex] = values[6];
        state.sorbed_acetate_g_c[complex] = values[7];
    }
}

fn take(backing: []f64, offset: *usize, count: usize) []f64 {
    const result = backing[offset.* .. offset.* + count];
    offset.* += count;
    return result;
}

fn validate(state: *const State, inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const values = @field(inputs, field.name);
        if (values.len != state.complex_count) return error.InvalidSorptionCommitDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidSorptionCommitInput;
        const signed = comptime std.mem.endsWith(u8, field.name, "_sorption_g_c") or
            std.mem.endsWith(u8, field.name, "_sorption_g_n") or
            std.mem.endsWith(u8, field.name, "_sorption_g_p");
        if (!signed) for (values) |value| if (value < 0)
            return error.InvalidSorptionCommitInput;
    }
}

fn fixture() Inputs {
    return .{
        .dissolved_organic_carbon_g_c = &.{ 10, 10 },
        .dissolved_organic_nitrogen_g_n = &.{ 2, 2 },
        .dissolved_organic_phosphorus_g_p = &.{ 1, 1 },
        .dissolved_acetate_g_c = &.{ 4, 4 },
        .sorbed_organic_carbon_g_c = &.{ 2, 2 },
        .sorbed_organic_nitrogen_g_n = &.{ 0.4, 0.4 },
        .sorbed_organic_phosphorus_g_p = &.{ 0.2, 0.2 },
        .sorbed_acetate_g_c = &.{ 1, 1 },
        .organic_carbon_sorption_g_c = &.{ 3, -1 },
        .organic_nitrogen_sorption_g_n = &.{ 0.6, -0.2 },
        .organic_phosphorus_sorption_g_p = &.{ 0.3, -0.1 },
        .acetate_sorption_g_c = &.{ 1, -0.5 },
    };
}

test "signed sorption commits adsorption and desorption" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectEqual(7, state.dissolved_organic_carbon_g_c[0]);
    try std.testing.expectEqual(5, state.sorbed_organic_carbon_g_c[0]);
    try std.testing.expectEqual(11, state.dissolved_organic_carbon_g_c[1]);
    try std.testing.expectEqual(1, state.sorbed_organic_carbon_g_c[1]);
}

test "every sorption commit conserves dissolved plus sorbed C N P and acetate" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    inline for (.{
        .{
            state.dissolved_organic_carbon_g_c,
            state.sorbed_organic_carbon_g_c,
            @as(f64, 12),
        },
        .{
            state.dissolved_organic_nitrogen_g_n,
            state.sorbed_organic_nitrogen_g_n,
            @as(f64, 2.4),
        },
        .{
            state.dissolved_organic_phosphorus_g_p,
            state.sorbed_organic_phosphorus_g_p,
            @as(f64, 1.2),
        },
        .{ state.dissolved_acetate_g_c, state.sorbed_acetate_g_c, @as(f64, 5) },
    }) |element| for (0..2) |complex|
        try std.testing.expectApproxEqAbs(
            element[2],
            element[0][complex] + element[1][complex],
            1e-12,
        );
}

test "desorption may consume the entire sorbed pool but not overdraw it" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.organic_carbon_sorption_g_c = &.{ 3, -2 };
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.sorbed_organic_carbon_g_c[1]);
}

test "overdraw fails atomically" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.dissolved_organic_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.organic_carbon_sorption_g_c = &.{ 11, -1 };
    try std.testing.expectError(error.InvalidSorptionCommitResult, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.dissolved_organic_carbon_g_c[0]);
}

test "NITRO 3671-3685 late sorbed overflow fails atomically" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.dissolved_organic_carbon_g_c[0] = 7;
    state.sorbed_acetate_g_c[1] = 11;
    var inputs = fixture();
    inputs.dissolved_acetate_g_c = &.{ 4, std.math.floatMax(f64) };
    inputs.sorbed_acetate_g_c = &.{ 1, std.math.floatMax(f64) };
    inputs.acetate_sorption_g_c = &.{ 1, std.math.floatMax(f64) };
    try std.testing.expectError(
        error.InvalidSorptionCommitResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.dissolved_organic_carbon_g_c[0]);
    try std.testing.expectEqual(11, state.sorbed_acetate_g_c[1]);
}
