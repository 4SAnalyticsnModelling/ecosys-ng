const std = @import("std");

pub const Inputs = struct {
    complex_count: usize,
    population_count: usize,
    residue_fraction_count: usize,
    dissolved_organic_carbon_g_c: []const f64,
    dissolved_organic_nitrogen_g_n: []const f64,
    dissolved_organic_phosphorus_g_p: []const f64,
    dissolved_acetate_g_c: []const f64,
    dissolved_organic_carbon_uptake_g_c: []const f64,
    dissolved_organic_nitrogen_uptake_g_n: []const f64,
    dissolved_organic_phosphorus_uptake_g_p: []const f64,
    acetate_uptake_g_c: []const f64,
    fermentation_acetate_production_g_c: []const f64,
    residue_carbon_g_c: []const f64,
    residue_nitrogen_g_n: []const f64,
    residue_phosphorus_g_p: []const f64,
    basal_litterfall_carbon_g_c: []const f64,
    basal_litterfall_nitrogen_g_n: []const f64,
    basal_litterfall_phosphorus_g_p: []const f64,
    redistributed_autotrophic_litter_carbon_g_c: []const f64,
    redistributed_autotrophic_litter_nitrogen_g_n: []const f64,
    redistributed_autotrophic_litter_phosphorus_g_p: []const f64,
    senescence_litterfall_carbon_g_c: []const f64,
    senescence_litterfall_nitrogen_g_n: []const f64,
    senescence_litterfall_phosphorus_g_p: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    backing: []f64,
    complex_count: usize,
    residue_item_count: usize,
    dissolved_organic_carbon_g_c: []f64,
    dissolved_organic_nitrogen_g_n: []f64,
    dissolved_organic_phosphorus_g_p: []f64,
    dissolved_acetate_g_c: []f64,
    residue_carbon_g_c: []f64,
    residue_nitrogen_g_n: []f64,
    residue_phosphorus_g_p: []f64,

    pub fn init(allocator: std.mem.Allocator, complex_count: usize, fraction_count: usize) !State {
        if (complex_count == 0 or fraction_count == 0)
            return error.InvalidUptakeLitterCommitDimensions;
        const residue_items = std.math.mul(usize, complex_count, fraction_count) catch
            return error.InvalidUptakeLitterCommitDimensions;
        const dissolved_values = std.math.mul(usize, complex_count, 4) catch
            return error.InvalidUptakeLitterCommitDimensions;
        const residue_values = std.math.mul(usize, residue_items, 3) catch
            return error.InvalidUptakeLitterCommitDimensions;
        const backing = try allocator.alloc(f64, dissolved_values + residue_values);
        @memset(backing, 0);
        var offset: usize = 0;
        return .{
            .allocator = allocator,
            .backing = backing,
            .complex_count = complex_count,
            .residue_item_count = residue_items,
            .dissolved_organic_carbon_g_c = take(backing, &offset, complex_count),
            .dissolved_organic_nitrogen_g_n = take(backing, &offset, complex_count),
            .dissolved_organic_phosphorus_g_p = take(backing, &offset, complex_count),
            .dissolved_acetate_g_c = take(backing, &offset, complex_count),
            .residue_carbon_g_c = take(backing, &offset, residue_items),
            .residue_nitrogen_g_n = take(backing, &offset, residue_items),
            .residue_phosphorus_g_p = take(backing, &offset, residue_items),
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.backing);
        self.* = undefined;
    }
};

/// Exact NITRO.F 3631--3669 population uptake and microbial-litter commit.
/// Traversal follows K, then N, then residue fraction M.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const pool_results = try state.allocator.alloc([4]f64, state.complex_count);
    defer state.allocator.free(pool_results);
    const residue_results = try state.allocator.alloc([3]f64, state.residue_item_count);
    defer state.allocator.free(residue_results);
    for (0..state.complex_count) |complex| {
        for (0..inputs.residue_fraction_count) |fraction| {
            const residue_item = complex * inputs.residue_fraction_count + fraction;
            residue_results[residue_item] = .{
                inputs.residue_carbon_g_c[residue_item],
                inputs.residue_nitrogen_g_n[residue_item],
                inputs.residue_phosphorus_g_p[residue_item],
            };
        }
        var carbon = inputs.dissolved_organic_carbon_g_c[complex];
        var nitrogen = inputs.dissolved_organic_nitrogen_g_n[complex];
        var phosphorus = inputs.dissolved_organic_phosphorus_g_p[complex];
        var acetate = inputs.dissolved_acetate_g_c[complex];
        for (0..inputs.population_count) |population| {
            const population_item = complex * inputs.population_count + population;
            carbon = carbon - inputs.dissolved_organic_carbon_uptake_g_c[population_item];
            nitrogen = nitrogen - inputs.dissolved_organic_nitrogen_uptake_g_n[population_item];
            phosphorus = phosphorus -
                inputs.dissolved_organic_phosphorus_uptake_g_p[population_item];
            acetate = acetate - inputs.acetate_uptake_g_c[population_item] +
                inputs.fermentation_acetate_production_g_c[population_item];
            for (0..inputs.residue_fraction_count) |fraction| {
                const residue_item = complex * inputs.residue_fraction_count + fraction;
                const flux_item = (complex * inputs.population_count + population) *
                    inputs.residue_fraction_count + fraction;
                residue_results[residue_item][0] =
                    residue_results[residue_item][0] +
                    inputs.basal_litterfall_carbon_g_c[flux_item] +
                    inputs.redistributed_autotrophic_litter_carbon_g_c[flux_item] +
                    inputs.senescence_litterfall_carbon_g_c[flux_item];
                residue_results[residue_item][1] =
                    residue_results[residue_item][1] +
                    inputs.basal_litterfall_nitrogen_g_n[flux_item] +
                    inputs.redistributed_autotrophic_litter_nitrogen_g_n[flux_item] +
                    inputs.senescence_litterfall_nitrogen_g_n[flux_item];
                residue_results[residue_item][2] =
                    residue_results[residue_item][2] +
                    inputs.basal_litterfall_phosphorus_g_p[flux_item] +
                    inputs.redistributed_autotrophic_litter_phosphorus_g_p[flux_item] +
                    inputs.senescence_litterfall_phosphorus_g_p[flux_item];
                inline for (residue_results[residue_item]) |value|
                    if (!std.math.isFinite(value) or value < 0)
                        return error.InvalidUptakeLitterCommitResult;
            }
            inline for (.{ carbon, nitrogen, phosphorus, acetate }) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.InvalidUptakeLitterCommitResult;
        }
        pool_results[complex] = .{ carbon, nitrogen, phosphorus, acetate };
    }
    for (pool_results) |values| inline for (values) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidUptakeLitterCommitResult;
    for (residue_results) |values| inline for (values) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidUptakeLitterCommitResult;
    for (pool_results, 0..) |values, complex| {
        state.dissolved_organic_carbon_g_c[complex] = values[0];
        state.dissolved_organic_nitrogen_g_n[complex] = values[1];
        state.dissolved_organic_phosphorus_g_p[complex] = values[2];
        state.dissolved_acetate_g_c[complex] = values[3];
    }
    for (residue_results, 0..) |values, item| {
        state.residue_carbon_g_c[item] = values[0];
        state.residue_nitrogen_g_n[item] = values[1];
        state.residue_phosphorus_g_p[item] = values[2];
    }
}

fn take(backing: []f64, offset: *usize, count: usize) []f64 {
    const result = backing[offset.* .. offset.* + count];
    offset.* += count;
    return result;
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.complex_count != state.complex_count or inputs.population_count == 0 or
        inputs.residue_fraction_count == 0 or
        state.residue_item_count != state.complex_count * inputs.residue_fraction_count)
        return error.InvalidUptakeLitterCommitDimensions;
    const population_items = state.complex_count * inputs.population_count;
    inline for (.{
        inputs.dissolved_organic_carbon_g_c,     inputs.dissolved_organic_nitrogen_g_n,
        inputs.dissolved_organic_phosphorus_g_p, inputs.dissolved_acetate_g_c,
    }) |values| if (values.len != state.complex_count)
        return error.InvalidUptakeLitterCommitDimensions;
    inline for (.{
        inputs.dissolved_organic_carbon_uptake_g_c,
        inputs.dissolved_organic_nitrogen_uptake_g_n,
        inputs.dissolved_organic_phosphorus_uptake_g_p,
        inputs.acetate_uptake_g_c,
        inputs.fermentation_acetate_production_g_c,
    }) |values| if (values.len != population_items)
        return error.InvalidUptakeLitterCommitDimensions;
    inline for (.{
        inputs.residue_carbon_g_c,     inputs.residue_nitrogen_g_n,
        inputs.residue_phosphorus_g_p,
    }) |values| if (values.len != state.residue_item_count)
        return error.InvalidUptakeLitterCommitDimensions;
    const flux_items = population_items * inputs.residue_fraction_count;
    inline for (.{
        inputs.basal_litterfall_carbon_g_c,                   inputs.basal_litterfall_nitrogen_g_n,
        inputs.basal_litterfall_phosphorus_g_p,               inputs.redistributed_autotrophic_litter_carbon_g_c,
        inputs.redistributed_autotrophic_litter_nitrogen_g_n, inputs.redistributed_autotrophic_litter_phosphorus_g_p,
        inputs.senescence_litterfall_carbon_g_c,              inputs.senescence_litterfall_nitrogen_g_n,
        inputs.senescence_litterfall_phosphorus_g_p,
    }) |values| if (values.len != flux_items)
        return error.InvalidUptakeLitterCommitDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64)
        for (@field(inputs, field.name)) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidUptakeLitterCommitInput;
}

fn fixture() Inputs {
    return .{
        .complex_count = 1,
        .population_count = 2,
        .residue_fraction_count = 2,
        .dissolved_organic_carbon_g_c = &.{10},
        .dissolved_organic_nitrogen_g_n = &.{4},
        .dissolved_organic_phosphorus_g_p = &.{2},
        .dissolved_acetate_g_c = &.{5},
        .dissolved_organic_carbon_uptake_g_c = &.{ 1, 2 },
        .dissolved_organic_nitrogen_uptake_g_n = &.{ 0.5, 0.5 },
        .dissolved_organic_phosphorus_uptake_g_p = &.{ 0.25, 0.25 },
        .acetate_uptake_g_c = &.{ 0.5, 1 },
        .fermentation_acetate_production_g_c = &.{ 1, 0 },
        .residue_carbon_g_c = &.{ 10, 20 },
        .residue_nitrogen_g_n = &.{ 2, 4 },
        .residue_phosphorus_g_p = &.{ 1, 2 },
        .basal_litterfall_carbon_g_c = &.{ 1, 2, 3, 4 },
        .basal_litterfall_nitrogen_g_n = &.{ 0.1, 0.2, 0.3, 0.4 },
        .basal_litterfall_phosphorus_g_p = &.{ 0.05, 0.1, 0.15, 0.2 },
        .redistributed_autotrophic_litter_carbon_g_c = &.{ 0.5, 1, 1.5, 2 },
        .redistributed_autotrophic_litter_nitrogen_g_n = &.{ 0.05, 0.1, 0.15, 0.2 },
        .redistributed_autotrophic_litter_phosphorus_g_p = &.{ 0.025, 0.05, 0.075, 0.1 },
        .senescence_litterfall_carbon_g_c = &.{ 0.25, 0.5, 0.75, 1 },
        .senescence_litterfall_nitrogen_g_n = &.{ 0.025, 0.05, 0.075, 0.1 },
        .senescence_litterfall_phosphorus_g_p = &.{ 0.0125, 0.025, 0.0375, 0.05 },
    };
}

test "population uptake and fermentation update dissolved pools" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectEqual(7, state.dissolved_organic_carbon_g_c[0]);
    try std.testing.expectEqual(3, state.dissolved_organic_nitrogen_g_n[0]);
    try std.testing.expectApproxEqAbs(4.5, state.dissolved_acetate_g_c[0], 1e-12);
}

test "all three litter sources accumulate by residue fraction" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(17, state.residue_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(30.5, state.residue_carbon_g_c[1], 1e-12);
}

test "runtime populations contribute independently" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.population_count = 1;
    inputs.dissolved_organic_carbon_uptake_g_c = &.{1};
    inputs.dissolved_organic_nitrogen_uptake_g_n = &.{0.5};
    inputs.dissolved_organic_phosphorus_uptake_g_p = &.{0.25};
    inputs.acetate_uptake_g_c = &.{0.5};
    inputs.fermentation_acetate_production_g_c = &.{1};
    inputs.basal_litterfall_carbon_g_c = &.{ 1, 2 };
    inputs.basal_litterfall_nitrogen_g_n = &.{ 0.1, 0.2 };
    inputs.basal_litterfall_phosphorus_g_p = &.{ 0.05, 0.1 };
    inputs.redistributed_autotrophic_litter_carbon_g_c = &.{ 0.5, 1 };
    inputs.redistributed_autotrophic_litter_nitrogen_g_n = &.{ 0.05, 0.1 };
    inputs.redistributed_autotrophic_litter_phosphorus_g_p = &.{ 0.025, 0.05 };
    inputs.senescence_litterfall_carbon_g_c = &.{ 0.25, 0.5 };
    inputs.senescence_litterfall_nitrogen_g_n = &.{ 0.025, 0.05 };
    inputs.senescence_litterfall_phosphorus_g_p = &.{ 0.0125, 0.025 };
    try calculate(&state, inputs);
    try std.testing.expectEqual(9, state.dissolved_organic_carbon_g_c[0]);
}

test "uptake and litter publication closes every explicit source and sink" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    const inputs = fixture();
    try calculate(&state, inputs);
    inline for (.{
        .{
            inputs.dissolved_organic_carbon_g_c[0],
            inputs.dissolved_organic_carbon_uptake_g_c,
            @as([]const f64, &.{}),
            state.dissolved_organic_carbon_g_c[0],
        },
        .{
            inputs.dissolved_organic_nitrogen_g_n[0],
            inputs.dissolved_organic_nitrogen_uptake_g_n,
            @as([]const f64, &.{}),
            state.dissolved_organic_nitrogen_g_n[0],
        },
        .{
            inputs.dissolved_organic_phosphorus_g_p[0],
            inputs.dissolved_organic_phosphorus_uptake_g_p,
            @as([]const f64, &.{}),
            state.dissolved_organic_phosphorus_g_p[0],
        },
        .{
            inputs.dissolved_acetate_g_c[0],
            inputs.acetate_uptake_g_c,
            inputs.fermentation_acetate_production_g_c,
            state.dissolved_acetate_g_c[0],
        },
    }) |pool| {
        var expected = pool[0];
        for (pool[1]) |value| expected -= value;
        for (pool[2]) |value| expected += value;
        try std.testing.expectApproxEqAbs(expected, pool[3], 1e-12);
    }
}

test "uptake overdraw fails atomically" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.dissolved_organic_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.dissolved_organic_carbon_uptake_g_c = &.{ 11, 2 };
    try std.testing.expectError(error.InvalidUptakeLitterCommitResult, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.dissolved_organic_carbon_g_c[0]);
}

test "NITRO 3631-3669 late litter overflow fails atomically" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.dissolved_organic_carbon_g_c[0] = 7;
    state.residue_phosphorus_g_p[1] = 11;
    var inputs = fixture();
    inputs.senescence_litterfall_phosphorus_g_p =
        &.{ 0.0125, 0.025, 0.0375, std.math.floatMax(f64) };
    inputs.redistributed_autotrophic_litter_phosphorus_g_p =
        &.{ 0.025, 0.05, 0.075, std.math.floatMax(f64) };
    try std.testing.expectError(
        error.InvalidUptakeLitterCommitResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.dissolved_organic_carbon_g_c[0]);
    try std.testing.expectEqual(11, state.residue_phosphorus_g_p[1]);
}
