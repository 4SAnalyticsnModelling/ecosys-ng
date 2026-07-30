const std = @import("std");

pub const Inputs = struct {
    residue_fraction_count: usize,
    residue_carbon_g_c: []const f64,
    residue_nitrogen_g_n: []const f64,
    residue_phosphorus_g_p: []const f64,
    decomposed_residue_carbon_g_c: []const f64,
    decomposed_residue_nitrogen_g_n: []const f64,
    decomposed_residue_phosphorus_g_p: []const f64,
    sorbed_carbon_g_c: []const f64,
    sorbed_nitrogen_g_n: []const f64,
    sorbed_phosphorus_g_p: []const f64,
    sorbed_acetate_g_c: []const f64,
    decomposed_sorbed_carbon_g_c: []const f64,
    decomposed_sorbed_nitrogen_g_n: []const f64,
    decomposed_sorbed_phosphorus_g_p: []const f64,
    decomposed_sorbed_acetate_g_c: []const f64,
    dissolved_organic_carbon_g_c: []const f64,
    dissolved_organic_nitrogen_g_n: []const f64,
    dissolved_organic_phosphorus_g_p: []const f64,
    dissolved_acetate_g_c: []const f64,
    chemodenitrification_dissolved_nitrogen_g_n: f64,
    chemodenitrification_distribution_fraction: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    backing: []f64,
    residue_item_count: usize,
    complex_count: usize,
    residue_carbon_g_c: []f64,
    residue_nitrogen_g_n: []f64,
    residue_phosphorus_g_p: []f64,
    sorbed_carbon_g_c: []f64,
    sorbed_nitrogen_g_n: []f64,
    sorbed_phosphorus_g_p: []f64,
    sorbed_acetate_g_c: []f64,
    dissolved_organic_carbon_g_c: []f64,
    dissolved_organic_nitrogen_g_n: []f64,
    dissolved_organic_phosphorus_g_p: []f64,
    dissolved_acetate_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, complex_count: usize, fraction_count: usize) !State {
        if (complex_count == 0 or fraction_count == 0)
            return error.InvalidResidueProductCommitDimensions;
        const items = std.math.mul(usize, complex_count, fraction_count) catch
            return error.InvalidResidueProductCommitDimensions;
        const residue_values = std.math.mul(usize, items, 3) catch
            return error.InvalidResidueProductCommitDimensions;
        const complex_values = std.math.mul(usize, complex_count, 8) catch
            return error.InvalidResidueProductCommitDimensions;
        const backing = try allocator.alloc(f64, residue_values + complex_values);
        @memset(backing, 0);
        var offset: usize = 0;
        return .{
            .allocator = allocator,
            .backing = backing,
            .residue_item_count = items,
            .complex_count = complex_count,
            .residue_carbon_g_c = take(backing, &offset, items),
            .residue_nitrogen_g_n = take(backing, &offset, items),
            .residue_phosphorus_g_p = take(backing, &offset, items),
            .sorbed_carbon_g_c = take(backing, &offset, complex_count),
            .sorbed_nitrogen_g_n = take(backing, &offset, complex_count),
            .sorbed_phosphorus_g_p = take(backing, &offset, complex_count),
            .sorbed_acetate_g_c = take(backing, &offset, complex_count),
            .dissolved_organic_carbon_g_c = take(backing, &offset, complex_count),
            .dissolved_organic_nitrogen_g_n = take(backing, &offset, complex_count),
            .dissolved_organic_phosphorus_g_p = take(backing, &offset, complex_count),
            .dissolved_acetate_g_c = take(backing, &offset, complex_count),
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.backing);
        self.* = undefined;
    }
};

/// Exact NITRO.F 3605--3629 residue and sorbed decomposition-product commit.
/// For each K, residue fraction M updates precede sorbed-product updates.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const residue_results = try state.allocator.alloc([3]f64, state.residue_item_count);
    defer state.allocator.free(residue_results);
    const pool_results = try state.allocator.alloc([8]f64, state.complex_count);
    defer state.allocator.free(pool_results);
    for (0..state.complex_count) |complex| {
        var dissolved_carbon = inputs.dissolved_organic_carbon_g_c[complex];
        var dissolved_nitrogen = inputs.dissolved_organic_nitrogen_g_n[complex];
        var dissolved_phosphorus = inputs.dissolved_organic_phosphorus_g_p[complex];
        for (0..inputs.residue_fraction_count) |fraction| {
            const item = complex * inputs.residue_fraction_count + fraction;
            const carbon = inputs.residue_carbon_g_c[item] -
                inputs.decomposed_residue_carbon_g_c[item];
            const nitrogen = inputs.residue_nitrogen_g_n[item] -
                inputs.decomposed_residue_nitrogen_g_n[item];
            const phosphorus = inputs.residue_phosphorus_g_p[item] -
                inputs.decomposed_residue_phosphorus_g_p[item];
            inline for (.{ carbon, nitrogen, phosphorus }) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.InvalidResidueProductCommitResult;
            residue_results[item] = .{ carbon, nitrogen, phosphorus };
            dissolved_carbon += inputs.decomposed_residue_carbon_g_c[item];
            dissolved_nitrogen += inputs.decomposed_residue_nitrogen_g_n[item];
            dissolved_phosphorus += inputs.decomposed_residue_phosphorus_g_p[item];
            inline for (.{ dissolved_carbon, dissolved_nitrogen, dissolved_phosphorus }) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.InvalidResidueProductCommitResult;
        }
        dissolved_carbon = dissolved_carbon +
            inputs.decomposed_sorbed_carbon_g_c[complex];
        dissolved_nitrogen = dissolved_nitrogen +
            inputs.decomposed_sorbed_nitrogen_g_n[complex] +
            inputs.chemodenitrification_dissolved_nitrogen_g_n *
                inputs.chemodenitrification_distribution_fraction[complex];
        dissolved_phosphorus = dissolved_phosphorus +
            inputs.decomposed_sorbed_phosphorus_g_p[complex];
        const dissolved_acetate = inputs.dissolved_acetate_g_c[complex] +
            inputs.decomposed_sorbed_acetate_g_c[complex];
        pool_results[complex] = .{
            inputs.sorbed_carbon_g_c[complex] - inputs.decomposed_sorbed_carbon_g_c[complex],
            inputs.sorbed_nitrogen_g_n[complex] - inputs.decomposed_sorbed_nitrogen_g_n[complex],
            inputs.sorbed_phosphorus_g_p[complex] -
                inputs.decomposed_sorbed_phosphorus_g_p[complex],
            inputs.sorbed_acetate_g_c[complex] -
                inputs.decomposed_sorbed_acetate_g_c[complex],
            dissolved_carbon,
            dissolved_nitrogen,
            dissolved_phosphorus,
            dissolved_acetate,
        };
    }
    for (pool_results) |values| inline for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidResidueProductCommitResult;
    for (residue_results, 0..) |values, item| {
        state.residue_carbon_g_c[item] = values[0];
        state.residue_nitrogen_g_n[item] = values[1];
        state.residue_phosphorus_g_p[item] = values[2];
    }
    for (pool_results, 0..) |values, complex| {
        state.sorbed_carbon_g_c[complex] = values[0];
        state.sorbed_nitrogen_g_n[complex] = values[1];
        state.sorbed_phosphorus_g_p[complex] = values[2];
        state.sorbed_acetate_g_c[complex] = values[3];
        state.dissolved_organic_carbon_g_c[complex] = values[4];
        state.dissolved_organic_nitrogen_g_n[complex] = values[5];
        state.dissolved_organic_phosphorus_g_p[complex] = values[6];
        state.dissolved_acetate_g_c[complex] = values[7];
    }
}

fn take(backing: []f64, offset: *usize, count: usize) []f64 {
    const result = backing[offset.* .. offset.* + count];
    offset.* += count;
    return result;
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.residue_fraction_count == 0 or
        state.residue_item_count != state.complex_count * inputs.residue_fraction_count)
        return error.InvalidResidueProductCommitDimensions;
    inline for (.{
        inputs.residue_carbon_g_c,              inputs.residue_nitrogen_g_n,
        inputs.residue_phosphorus_g_p,          inputs.decomposed_residue_carbon_g_c,
        inputs.decomposed_residue_nitrogen_g_n, inputs.decomposed_residue_phosphorus_g_p,
    }) |values| {
        if (values.len != state.residue_item_count)
            return error.InvalidResidueProductCommitDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidResidueProductCommitInput;
    }
    inline for (.{
        inputs.sorbed_carbon_g_c,                          inputs.sorbed_nitrogen_g_n,
        inputs.sorbed_phosphorus_g_p,                      inputs.sorbed_acetate_g_c,
        inputs.decomposed_sorbed_carbon_g_c,               inputs.decomposed_sorbed_nitrogen_g_n,
        inputs.decomposed_sorbed_phosphorus_g_p,           inputs.decomposed_sorbed_acetate_g_c,
        inputs.dissolved_organic_carbon_g_c,               inputs.dissolved_organic_nitrogen_g_n,
        inputs.dissolved_organic_phosphorus_g_p,           inputs.dissolved_acetate_g_c,
        inputs.chemodenitrification_distribution_fraction,
    }) |values| {
        if (values.len != state.complex_count)
            return error.InvalidResidueProductCommitDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidResidueProductCommitInput;
    }
    if (!std.math.isFinite(inputs.chemodenitrification_dissolved_nitrogen_g_n) or
        inputs.chemodenitrification_dissolved_nitrogen_g_n < 0)
        return error.InvalidResidueProductCommitInput;
}

fn fixture() Inputs {
    return .{
        .residue_fraction_count = 2,
        .residue_carbon_g_c = &.{ 10, 20, 30, 40 },
        .residue_nitrogen_g_n = &.{ 2, 4, 6, 8 },
        .residue_phosphorus_g_p = &.{ 1, 2, 3, 4 },
        .decomposed_residue_carbon_g_c = &.{ 1, 2, 3, 4 },
        .decomposed_residue_nitrogen_g_n = &.{ 0.2, 0.4, 0.6, 0.8 },
        .decomposed_residue_phosphorus_g_p = &.{ 0.1, 0.2, 0.3, 0.4 },
        .sorbed_carbon_g_c = &.{ 5, 5 },
        .sorbed_nitrogen_g_n = &.{ 1, 1 },
        .sorbed_phosphorus_g_p = &.{ 0.5, 0.5 },
        .sorbed_acetate_g_c = &.{ 2, 2 },
        .decomposed_sorbed_carbon_g_c = &.{ 0.5, 0.5 },
        .decomposed_sorbed_nitrogen_g_n = &.{ 0.1, 0.1 },
        .decomposed_sorbed_phosphorus_g_p = &.{ 0.05, 0.05 },
        .decomposed_sorbed_acetate_g_c = &.{ 0.2, 0.2 },
        .dissolved_organic_carbon_g_c = &.{ 5, 5 },
        .dissolved_organic_nitrogen_g_n = &.{ 1, 1 },
        .dissolved_organic_phosphorus_g_p = &.{ 0.5, 0.5 },
        .dissolved_acetate_g_c = &.{ 2, 2 },
        .chemodenitrification_dissolved_nitrogen_g_n = 1,
        .chemodenitrification_distribution_fraction = &.{ 0.25, 0.75 },
    };
}

test "commit moves residue and sorbed products to dissolved pools" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectEqual(9, state.residue_carbon_g_c[0]);
    try std.testing.expectEqual(4.5, state.sorbed_carbon_g_c[0]);
    try std.testing.expectApproxEqAbs(8.5, state.dissolved_organic_carbon_g_c[0], 1e-12);
}

test "chemodenitrification DON follows residue distribution" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(1.95, state.dissolved_organic_nitrogen_g_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(3.25, state.dissolved_organic_nitrogen_g_n[1], 1e-12);
}

test "element totals close except explicit chemodenitrification addition" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    var total_don: f64 = 0;
    for (state.dissolved_organic_nitrogen_g_n) |value| total_don += value;
    try std.testing.expectApproxEqAbs(5.2, total_don, 1e-12);
}

test "residue and sorbed publication conserves C N P and acetate" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    const inputs = fixture();
    try calculate(&state, inputs);
    inline for (.{
        .{
            inputs.residue_carbon_g_c,
            state.residue_carbon_g_c,
            inputs.sorbed_carbon_g_c,
            state.sorbed_carbon_g_c,
            inputs.dissolved_organic_carbon_g_c,
            state.dissolved_organic_carbon_g_c,
            @as(f64, 0),
        },
        .{
            inputs.residue_nitrogen_g_n,
            state.residue_nitrogen_g_n,
            inputs.sorbed_nitrogen_g_n,
            state.sorbed_nitrogen_g_n,
            inputs.dissolved_organic_nitrogen_g_n,
            state.dissolved_organic_nitrogen_g_n,
            inputs.chemodenitrification_dissolved_nitrogen_g_n,
        },
        .{
            inputs.residue_phosphorus_g_p,
            state.residue_phosphorus_g_p,
            inputs.sorbed_phosphorus_g_p,
            state.sorbed_phosphorus_g_p,
            inputs.dissolved_organic_phosphorus_g_p,
            state.dissolved_organic_phosphorus_g_p,
            @as(f64, 0),
        },
    }) |element| {
        var initial: f64 = element[6];
        var final: f64 = 0;
        for (element[0]) |value| initial += value;
        for (element[2]) |value| initial += value;
        for (element[4]) |value| initial += value;
        for (element[1]) |value| final += value;
        for (element[3]) |value| final += value;
        for (element[5]) |value| final += value;
        try std.testing.expectApproxEqAbs(initial, final, 1e-12);
    }
    var initial_acetate: f64 = 0;
    var final_acetate: f64 = 0;
    for (inputs.sorbed_acetate_g_c) |value| initial_acetate += value;
    for (inputs.dissolved_acetate_g_c) |value| initial_acetate += value;
    for (state.sorbed_acetate_g_c) |value| final_acetate += value;
    for (state.dissolved_acetate_g_c) |value| final_acetate += value;
    try std.testing.expectApproxEqAbs(initial_acetate, final_acetate, 1e-12);
}

test "overdraw fails atomically" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    state.residue_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.decomposed_sorbed_carbon_g_c = &.{ 6, 0.5 };
    try std.testing.expectError(error.InvalidResidueProductCommitResult, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.residue_carbon_g_c[0]);
}

test "NITRO 3605-3629 late dissolved overflow fails atomically" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    state.residue_carbon_g_c[0] = 7;
    state.dissolved_organic_nitrogen_g_n[1] = 11;
    var inputs = fixture();
    inputs.chemodenitrification_dissolved_nitrogen_g_n =
        std.math.floatMax(f64);
    inputs.chemodenitrification_distribution_fraction =
        &.{ 0.25, std.math.floatMax(f64) };
    try std.testing.expectError(
        error.InvalidResidueProductCommitResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.residue_carbon_g_c[0]);
    try std.testing.expectEqual(11, state.dissolved_organic_nitrogen_g_n[1]);
}
