const std = @import("std");

pub const Inputs = struct {
    complex_count: usize,
    population_count: usize,
    biomass_fraction_count: usize,
    microbial_residue_carbon_by_complex_g_c: []const f64,
    total_microbial_residue_carbon_g_c: f64,
    autotrophic_litterfall_carbon_g_c: []const f64,
    autotrophic_senescence_litterfall_carbon_g_c: []const f64,
    autotrophic_litterfall_nitrogen_g_n: []const f64,
    autotrophic_senescence_litterfall_nitrogen_g_n: []const f64,
    autotrophic_litterfall_phosphorus_g_p: []const f64,
    autotrophic_senescence_litterfall_phosphorus_g_p: []const f64,
    zero_residue_fallback_complex: usize,
    negligible_residue_carbon_g_c: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    backing: []f64,
    complex_count: usize,
    source_item_count: usize,
    destination_item_count: usize,
    destination_fraction: []f64,
    transferred_carbon_g_c: []f64,
    transferred_nitrogen_g_n: []f64,
    transferred_phosphorus_g_p: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        complex_count: usize,
        population_count: usize,
        biomass_fraction_count: usize,
    ) !State {
        if (complex_count == 0 or population_count == 0 or biomass_fraction_count == 0)
            return error.InvalidAutotrophicLitterDimensions;
        const source_items = std.math.mul(usize, population_count, biomass_fraction_count) catch
            return error.InvalidAutotrophicLitterDimensions;
        const destination_items = std.math.mul(usize, complex_count, source_items) catch
            return error.InvalidAutotrophicLitterDimensions;
        const transfer_values = std.math.mul(usize, destination_items, 3) catch
            return error.InvalidAutotrophicLitterDimensions;
        const count = std.math.add(usize, complex_count, transfer_values) catch
            return error.InvalidAutotrophicLitterDimensions;
        const backing = try allocator.alloc(f64, count);
        @memset(backing, 0);
        return .{
            .allocator = allocator,
            .backing = backing,
            .complex_count = complex_count,
            .source_item_count = source_items,
            .destination_item_count = destination_items,
            .destination_fraction = backing[0..complex_count],
            .transferred_carbon_g_c = backing[complex_count .. complex_count + destination_items],
            .transferred_nitrogen_g_n = backing[complex_count + destination_items .. complex_count + 2 * destination_items],
            .transferred_phosphorus_g_p = backing[complex_count + 2 * destination_items ..],
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.backing);
        self.* = undefined;
    }
};

/// Exact NITRO.F 3527--3558 redistribution of K=5 autotrophic litterfall.
/// Traversal follows K, then N, then M.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc(f64, state.backing.len);
    defer state.allocator.free(temporary);
    const fractions = temporary[0..state.complex_count];
    const carbon = temporary[state.complex_count .. state.complex_count + state.destination_item_count];
    const nitrogen = temporary[state.complex_count + state.destination_item_count .. state.complex_count + 2 * state.destination_item_count];
    const phosphorus = temporary[state.complex_count + 2 * state.destination_item_count ..];
    for (0..state.complex_count) |complex| {
        fractions[complex] = if (inputs.total_microbial_residue_carbon_g_c >
            inputs.negligible_residue_carbon_g_c)
            inputs.microbial_residue_carbon_by_complex_g_c[complex] /
                inputs.total_microbial_residue_carbon_g_c
        else if (complex == inputs.zero_residue_fallback_complex) 1 else 0;
        if (!std.math.isFinite(fractions[complex]) or fractions[complex] < 0)
            return error.InvalidAutotrophicLitterDistribution;
        for (0..inputs.population_count) |population| {
            for (0..inputs.biomass_fraction_count) |fraction| {
                const source_item = population * inputs.biomass_fraction_count + fraction;
                const destination = complex * state.source_item_count + source_item;
                carbon[destination] =
                    (inputs.autotrophic_litterfall_carbon_g_c[source_item] +
                        inputs.autotrophic_senescence_litterfall_carbon_g_c[source_item]) *
                    fractions[complex];
                nitrogen[destination] =
                    (inputs.autotrophic_litterfall_nitrogen_g_n[source_item] +
                        inputs.autotrophic_senescence_litterfall_nitrogen_g_n[source_item]) *
                    fractions[complex];
                phosphorus[destination] =
                    (inputs.autotrophic_litterfall_phosphorus_g_p[source_item] +
                        inputs.autotrophic_senescence_litterfall_phosphorus_g_p[source_item]) *
                    fractions[complex];
                inline for (.{ carbon[destination], nitrogen[destination], phosphorus[destination] }) |value|
                    if (!std.math.isFinite(value) or value < 0)
                        return error.NonFiniteAutotrophicLitterDistribution;
            }
        }
    }
    @memcpy(state.backing, temporary);
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.complex_count != state.complex_count or inputs.population_count == 0 or
        inputs.biomass_fraction_count == 0 or
        inputs.zero_residue_fallback_complex >= state.complex_count)
        return error.InvalidAutotrophicLitterDimensions;
    const source_items = std.math.mul(usize, inputs.population_count, inputs.biomass_fraction_count) catch return error.InvalidAutotrophicLitterDimensions;
    if (source_items != state.source_item_count or
        inputs.microbial_residue_carbon_by_complex_g_c.len != state.complex_count)
        return error.InvalidAutotrophicLitterDimensions;
    inline for (.{
        inputs.autotrophic_litterfall_carbon_g_c,
        inputs.autotrophic_senescence_litterfall_carbon_g_c,
        inputs.autotrophic_litterfall_nitrogen_g_n,
        inputs.autotrophic_senescence_litterfall_nitrogen_g_n,
        inputs.autotrophic_litterfall_phosphorus_g_p,
        inputs.autotrophic_senescence_litterfall_phosphorus_g_p,
    }) |values| if (values.len != source_items) return error.InvalidAutotrophicLitterDimensions;
    inline for (.{
        inputs.microbial_residue_carbon_by_complex_g_c,
        inputs.autotrophic_litterfall_carbon_g_c,
        inputs.autotrophic_senescence_litterfall_carbon_g_c,
        inputs.autotrophic_litterfall_nitrogen_g_n,
        inputs.autotrophic_senescence_litterfall_nitrogen_g_n,
        inputs.autotrophic_litterfall_phosphorus_g_p,
        inputs.autotrophic_senescence_litterfall_phosphorus_g_p,
    }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidAutotrophicLitterInput;
    inline for (.{
        inputs.total_microbial_residue_carbon_g_c,
        inputs.negligible_residue_carbon_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidAutotrophicLitterInput;
}

fn fixture() Inputs {
    return .{
        .complex_count = 4,
        .population_count = 1,
        .biomass_fraction_count = 2,
        .microbial_residue_carbon_by_complex_g_c = &.{ 1, 2, 3, 4 },
        .total_microbial_residue_carbon_g_c = 10,
        .autotrophic_litterfall_carbon_g_c = &.{ 2, 4 },
        .autotrophic_senescence_litterfall_carbon_g_c = &.{ 1, 2 },
        .autotrophic_litterfall_nitrogen_g_n = &.{ 0.2, 0.4 },
        .autotrophic_senescence_litterfall_nitrogen_g_n = &.{ 0.1, 0.2 },
        .autotrophic_litterfall_phosphorus_g_p = &.{ 0.1, 0.2 },
        .autotrophic_senescence_litterfall_phosphorus_g_p = &.{ 0.05, 0.1 },
        .zero_residue_fallback_complex = 3,
        .negligible_residue_carbon_g_c = 1e-12,
    };
}

test "autotrophic litter follows microbial residue distribution" {
    var state = try State.init(std.testing.allocator, 4, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(0.1, state.destination_fraction[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.4, state.destination_fraction[3], 1e-12);
    try std.testing.expectApproxEqAbs(0.3, state.transferred_carbon_g_c[0], 1e-12);
}

test "redistribution conserves every source product" {
    var state = try State.init(std.testing.allocator, 4, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    inline for (.{
        .{ state.transferred_carbon_g_c, @as(f64, 3), @as(f64, 6) },
        .{ state.transferred_nitrogen_g_n, @as(f64, 0.3), @as(f64, 0.6) },
        .{ state.transferred_phosphorus_g_p, @as(f64, 0.15), @as(f64, 0.3) },
    }) |element| for (0..2) |source_item| {
        var total: f64 = 0;
        for (0..4) |complex| total += element[0][complex * 2 + source_item];
        const expected = if (source_item == 0) element[1] else element[2];
        try std.testing.expectApproxEqAbs(expected, total, 1e-12);
    };
}

test "zero residue falls back to source K equals three" {
    var state = try State.init(std.testing.allocator, 4, 1, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.total_microbial_residue_carbon_g_c = 0;
    try calculate(&state, inputs);
    try std.testing.expectEqual(1, state.destination_fraction[3]);
    try std.testing.expectEqual(3, state.transferred_carbon_g_c[6]);
}

test "NITRO TORC threshold equality uses only fallback complex" {
    var state = try State.init(std.testing.allocator, 4, 1, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.total_microbial_residue_carbon_g_c = 1e-12;
    try calculate(&state, inputs);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0, 1 }, state.destination_fraction);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 4, 1, 2);
    defer state.deinit();
    state.transferred_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.total_microbial_residue_carbon_g_c = std.math.nan(f64);
    try std.testing.expectError(error.InvalidAutotrophicLitterInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.transferred_carbon_g_c[0]);
}

test "NITRO 3527-3558 late overflow preserves all redistribution outputs" {
    var state = try State.init(std.testing.allocator, 4, 1, 2);
    defer state.deinit();
    state.destination_fraction[0] = 7;
    state.transferred_phosphorus_g_p[7] = 11;
    var inputs = fixture();
    inputs.autotrophic_litterfall_phosphorus_g_p =
        &.{ 0.1, std.math.floatMax(f64) };
    inputs.autotrophic_senescence_litterfall_phosphorus_g_p =
        &.{ 0.05, std.math.floatMax(f64) };
    try std.testing.expectError(
        error.NonFiniteAutotrophicLitterDistribution,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.destination_fraction[0]);
    try std.testing.expectEqual(11, state.transferred_phosphorus_g_p[7]);
}
