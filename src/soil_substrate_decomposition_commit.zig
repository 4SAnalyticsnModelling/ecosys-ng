const std = @import("std");

pub const Inputs = struct {
    complex_count: usize,
    substrate_class_count: usize,
    surface_litter_layer: bool,
    substrate_carbon_g_c: []const f64,
    colonized_substrate_carbon_g_c: []const f64,
    substrate_nitrogen_g_n: []const f64,
    substrate_phosphorus_g_p: []const f64,
    decomposed_carbon_g_c: []const f64,
    decomposed_nitrogen_g_n: []const f64,
    decomposed_phosphorus_g_p: []const f64,
    dissolved_product_carbon_g_c: []const f64,
    dissolved_product_nitrogen_g_n: []const f64,
    dissolved_product_phosphorus_g_p: []const f64,
    humified_product_carbon_g_c: []const f64,
    humified_product_nitrogen_g_n: []const f64,
    humified_product_phosphorus_g_p: []const f64,
    dissolved_organic_carbon_g_c: []const f64,
    dissolved_organic_nitrogen_g_n: []const f64,
    dissolved_organic_phosphorus_g_p: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    backing: []f64,
    item_count: usize,
    complex_count: usize,
    substrate_carbon_g_c: []f64,
    colonized_substrate_carbon_g_c: []f64,
    substrate_nitrogen_g_n: []f64,
    substrate_phosphorus_g_p: []f64,
    dissolved_organic_carbon_g_c: []f64,
    dissolved_organic_nitrogen_g_n: []f64,
    dissolved_organic_phosphorus_g_p: []f64,
    particulate_carbon_addition_g_c: f64,
    colonized_particulate_carbon_addition_g_c: f64,
    particulate_nitrogen_addition_g_n: f64,
    particulate_phosphorus_addition_g_p: f64,
    particulate_target_is_surface_soil_layer: bool,

    pub fn init(allocator: std.mem.Allocator, complex_count: usize, class_count: usize) !State {
        if (complex_count == 0 or class_count == 0)
            return error.InvalidSubstrateCommitDimensions;
        const items = std.math.mul(usize, complex_count, class_count) catch
            return error.InvalidSubstrateCommitDimensions;
        const item_values = std.math.mul(usize, items, 4) catch
            return error.InvalidSubstrateCommitDimensions;
        const complex_values = std.math.mul(usize, complex_count, 3) catch
            return error.InvalidSubstrateCommitDimensions;
        const count = std.math.add(usize, item_values, complex_values) catch
            return error.InvalidSubstrateCommitDimensions;
        const backing = try allocator.alloc(f64, count);
        @memset(backing, 0);
        var offset: usize = 0;
        return .{
            .allocator = allocator,
            .backing = backing,
            .item_count = items,
            .complex_count = complex_count,
            .substrate_carbon_g_c = take(backing, &offset, items),
            .colonized_substrate_carbon_g_c = take(backing, &offset, items),
            .substrate_nitrogen_g_n = take(backing, &offset, items),
            .substrate_phosphorus_g_p = take(backing, &offset, items),
            .dissolved_organic_carbon_g_c = take(backing, &offset, complex_count),
            .dissolved_organic_nitrogen_g_n = take(backing, &offset, complex_count),
            .dissolved_organic_phosphorus_g_p = take(backing, &offset, complex_count),
            .particulate_carbon_addition_g_c = 0,
            .colonized_particulate_carbon_addition_g_c = 0,
            .particulate_nitrogen_addition_g_n = 0,
            .particulate_phosphorus_addition_g_p = 0,
            .particulate_target_is_surface_soil_layer = false,
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.backing);
        self.* = undefined;
    }
};

/// Exact NITRO.F 3560--3603 substrate-product state redistribution.
/// Traversal follows K, then structural fraction M.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const item_results = try state.allocator.alloc([4]f64, state.item_count);
    defer state.allocator.free(item_results);
    const dissolved_results = try state.allocator.alloc([3]f64, state.complex_count);
    defer state.allocator.free(dissolved_results);
    var particulate_carbon: f64 = 0;
    var particulate_nitrogen: f64 = 0;
    var particulate_phosphorus: f64 = 0;
    for (0..state.complex_count) |complex| {
        var dissolved_carbon = inputs.dissolved_organic_carbon_g_c[complex];
        var dissolved_nitrogen = inputs.dissolved_organic_nitrogen_g_n[complex];
        var dissolved_phosphorus = inputs.dissolved_organic_phosphorus_g_p[complex];
        for (0..inputs.substrate_class_count) |class| {
            const item = complex * inputs.substrate_class_count + class;
            const carbon = inputs.substrate_carbon_g_c[item] -
                inputs.decomposed_carbon_g_c[item];
            const colonized = inputs.colonized_substrate_carbon_g_c[item] -
                inputs.decomposed_carbon_g_c[item];
            const nitrogen = inputs.substrate_nitrogen_g_n[item] -
                inputs.decomposed_nitrogen_g_n[item];
            const phosphorus = inputs.substrate_phosphorus_g_p[item] -
                inputs.decomposed_phosphorus_g_p[item];
            inline for (.{ carbon, colonized, nitrogen, phosphorus }) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.InvalidSubstrateCommitResult;
            item_results[item] = .{ carbon, colonized, nitrogen, phosphorus };
            dissolved_carbon += inputs.dissolved_product_carbon_g_c[item];
            dissolved_nitrogen += inputs.dissolved_product_nitrogen_g_n[item];
            dissolved_phosphorus += inputs.dissolved_product_phosphorus_g_p[item];
            particulate_carbon += inputs.humified_product_carbon_g_c[item];
            particulate_nitrogen += inputs.humified_product_nitrogen_g_n[item];
            particulate_phosphorus += inputs.humified_product_phosphorus_g_p[item];
            inline for (.{
                dissolved_carbon,   dissolved_nitrogen,   dissolved_phosphorus,
                particulate_carbon, particulate_nitrogen, particulate_phosphorus,
            }) |value| if (!std.math.isFinite(value) or value < 0)
                return error.InvalidSubstrateCommitResult;
        }
        dissolved_results[complex] = .{
            dissolved_carbon, dissolved_nitrogen, dissolved_phosphorus,
        };
    }
    for (dissolved_results) |values| inline for (values) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSubstrateCommitResult;
    for (item_results, 0..) |values, item| {
        state.substrate_carbon_g_c[item] = values[0];
        state.colonized_substrate_carbon_g_c[item] = values[1];
        state.substrate_nitrogen_g_n[item] = values[2];
        state.substrate_phosphorus_g_p[item] = values[3];
    }
    for (dissolved_results, 0..) |values, complex| {
        state.dissolved_organic_carbon_g_c[complex] = values[0];
        state.dissolved_organic_nitrogen_g_n[complex] = values[1];
        state.dissolved_organic_phosphorus_g_p[complex] = values[2];
    }
    state.particulate_carbon_addition_g_c = particulate_carbon;
    state.colonized_particulate_carbon_addition_g_c = particulate_carbon;
    state.particulate_nitrogen_addition_g_n = particulate_nitrogen;
    state.particulate_phosphorus_addition_g_p = particulate_phosphorus;
    state.particulate_target_is_surface_soil_layer = inputs.surface_litter_layer;
}

fn take(backing: []f64, offset: *usize, count: usize) []f64 {
    const result = backing[offset.* .. offset.* + count];
    offset.* += count;
    return result;
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.complex_count != state.complex_count or inputs.substrate_class_count == 0 or
        inputs.complex_count * inputs.substrate_class_count != state.item_count)
        return error.InvalidSubstrateCommitDimensions;
    inline for (.{
        inputs.substrate_carbon_g_c,            inputs.colonized_substrate_carbon_g_c,
        inputs.substrate_nitrogen_g_n,          inputs.substrate_phosphorus_g_p,
        inputs.decomposed_carbon_g_c,           inputs.decomposed_nitrogen_g_n,
        inputs.decomposed_phosphorus_g_p,       inputs.dissolved_product_carbon_g_c,
        inputs.dissolved_product_nitrogen_g_n,  inputs.dissolved_product_phosphorus_g_p,
        inputs.humified_product_carbon_g_c,     inputs.humified_product_nitrogen_g_n,
        inputs.humified_product_phosphorus_g_p,
    }) |values| {
        if (values.len != state.item_count) return error.InvalidSubstrateCommitDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSubstrateCommitInput;
    }
    inline for (.{
        inputs.dissolved_organic_carbon_g_c,     inputs.dissolved_organic_nitrogen_g_n,
        inputs.dissolved_organic_phosphorus_g_p,
    }) |values| {
        if (values.len != state.complex_count) return error.InvalidSubstrateCommitDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSubstrateCommitInput;
    }
}

fn fixture() Inputs {
    return .{
        .complex_count = 2,
        .substrate_class_count = 2,
        .surface_litter_layer = false,
        .substrate_carbon_g_c = &.{ 10, 20, 30, 40 },
        .colonized_substrate_carbon_g_c = &.{ 5, 10, 15, 20 },
        .substrate_nitrogen_g_n = &.{ 2, 4, 6, 8 },
        .substrate_phosphorus_g_p = &.{ 1, 2, 3, 4 },
        .decomposed_carbon_g_c = &.{ 1, 2, 3, 4 },
        .decomposed_nitrogen_g_n = &.{ 0.2, 0.4, 0.6, 0.8 },
        .decomposed_phosphorus_g_p = &.{ 0.1, 0.2, 0.3, 0.4 },
        .dissolved_product_carbon_g_c = &.{ 0.5, 1, 1.5, 2 },
        .dissolved_product_nitrogen_g_n = &.{ 0.1, 0.2, 0.3, 0.4 },
        .dissolved_product_phosphorus_g_p = &.{ 0.05, 0.1, 0.15, 0.2 },
        .humified_product_carbon_g_c = &.{ 0.5, 1, 1.5, 2 },
        .humified_product_nitrogen_g_n = &.{ 0.1, 0.2, 0.3, 0.4 },
        .humified_product_phosphorus_g_p = &.{ 0.05, 0.1, 0.15, 0.2 },
        .dissolved_organic_carbon_g_c = &.{ 5, 5 },
        .dissolved_organic_nitrogen_g_n = &.{ 1, 1 },
        .dissolved_organic_phosphorus_g_p = &.{ 0.5, 0.5 },
    };
}

test "commit subtracts decomposition and adds dissolved products" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectEqual(9, state.substrate_carbon_g_c[0]);
    try std.testing.expectEqual(4, state.colonized_substrate_carbon_g_c[0]);
    try std.testing.expectApproxEqAbs(6.5, state.dissolved_organic_carbon_g_c[0], 1e-12);
}

test "humified products aggregate to particulate target" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(5, state.particulate_carbon_addition_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(
        state.particulate_carbon_addition_g_c,
        state.colonized_particulate_carbon_addition_g_c,
        1e-12,
    );
}

test "surface litter redirects particulate target to surface soil" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.surface_litter_layer = true;
    try calculate(&state, inputs);
    try std.testing.expect(state.particulate_target_is_surface_soil_layer);
}

test "structural decomposition routing closes carbon nitrogen and phosphorus" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    const inputs = fixture();
    try calculate(&state, inputs);
    inline for (.{
        .{
            inputs.substrate_carbon_g_c,
            state.substrate_carbon_g_c,
            inputs.dissolved_organic_carbon_g_c,
            state.dissolved_organic_carbon_g_c,
            state.particulate_carbon_addition_g_c,
        },
        .{
            inputs.substrate_nitrogen_g_n,
            state.substrate_nitrogen_g_n,
            inputs.dissolved_organic_nitrogen_g_n,
            state.dissolved_organic_nitrogen_g_n,
            state.particulate_nitrogen_addition_g_n,
        },
        .{
            inputs.substrate_phosphorus_g_p,
            state.substrate_phosphorus_g_p,
            inputs.dissolved_organic_phosphorus_g_p,
            state.dissolved_organic_phosphorus_g_p,
            state.particulate_phosphorus_addition_g_p,
        },
    }) |element| {
        var initial: f64 = 0;
        var final: f64 = element[4];
        for (element[0]) |value| initial += value;
        for (element[2]) |value| initial += value;
        for (element[1]) |value| final += value;
        for (element[3]) |value| final += value;
        try std.testing.expectApproxEqAbs(initial, final, 1e-12);
    }
}

test "overdraw fails atomically" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    state.substrate_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.decomposed_carbon_g_c = &.{ 6, 2, 3, 4 };
    try std.testing.expectError(error.InvalidSubstrateCommitResult, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.substrate_carbon_g_c[0]);
}

test "NITRO 3560-3603 late running-total overflow fails atomically" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    state.substrate_carbon_g_c[0] = 7;
    state.dissolved_organic_phosphorus_g_p[1] = 11;
    state.particulate_phosphorus_addition_g_p = 13;
    var inputs = fixture();
    inputs.humified_product_phosphorus_g_p =
        &.{ 0.05, 0.1, std.math.floatMax(f64), std.math.floatMax(f64) };
    try std.testing.expectError(
        error.InvalidSubstrateCommitResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.substrate_carbon_g_c[0]);
    try std.testing.expectEqual(11, state.dissolved_organic_phosphorus_g_p[1]);
    try std.testing.expectEqual(13, state.particulate_phosphorus_addition_g_p);
}
