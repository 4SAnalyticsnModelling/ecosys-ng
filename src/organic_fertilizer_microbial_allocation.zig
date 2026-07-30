const std = @import("std");

pub const ElementPool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const Inputs = struct {
    amendment: ElementPool,
    kinetic_fraction_by_class: []const f64,
    heterotroph_fraction_by_population: []const f64,
    nitrogen_to_carbon_by_population_class: []const f64,
    phosphorus_to_carbon_by_population_class: []const f64,
    autotroph_fraction_by_population: []const f64,
};

pub const Outputs = struct {
    /// Flattened heterotroph-population-major, then kinetic class.
    heterotroph_delta: []ElementPool,
    /// Flattened kinetic-class-major, then autotroph population.
    autotroph_delta: []ElementPool,
};

/// HOUR1 lines 797--825. Runtime extents replace source fixed 7x3 and 7
/// dimensions without changing N -> M -> NN traversal or addition order.
pub fn compute(inputs: Inputs, outputs: Outputs) !ElementPool {
    try validate(inputs, outputs);
    @memset(outputs.heterotroph_delta, .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 });
    @memset(outputs.autotroph_delta, .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 });

    const class_count = inputs.kinetic_fraction_by_class.len;
    const autotroph_count = inputs.autotroph_fraction_by_population.len;
    var allocated: ElementPool = .{
        .carbon_g_c = 0,
        .nitrogen_g_n = 0,
        .phosphorus_g_p = 0,
    };
    for (inputs.heterotroph_fraction_by_population, 0..) |population_fraction, population| {
        for (inputs.kinetic_fraction_by_class, 0..) |kinetic_fraction, class| {
            const index = population * class_count + class;
            const carbon = @max(0, @min(
                inputs.amendment.carbon_g_c *
                    kinetic_fraction *
                    population_fraction,
                inputs.amendment.carbon_g_c - allocated.carbon_g_c,
            ));
            const nitrogen = @max(0, @min(
                carbon * inputs.nitrogen_to_carbon_by_population_class[index],
                inputs.amendment.nitrogen_g_n - allocated.nitrogen_g_n,
            ));
            const phosphorus = @max(0, @min(
                carbon * inputs.phosphorus_to_carbon_by_population_class[index],
                inputs.amendment.phosphorus_g_p - allocated.phosphorus_g_p,
            ));
            outputs.heterotroph_delta[index] = .{
                .carbon_g_c = carbon,
                .nitrogen_g_n = nitrogen,
                .phosphorus_g_p = phosphorus,
            };
            allocated.carbon_g_c += carbon;
            allocated.nitrogen_g_n += nitrogen;
            allocated.phosphorus_g_p += phosphorus;

            for (inputs.autotroph_fraction_by_population, 0..) |autotroph_fraction, autotroph| {
                const autotroph_index = class * autotroph_count + autotroph;
                const carbon_delta = carbon * autotroph_fraction;
                const nitrogen_delta = nitrogen * autotroph_fraction;
                const phosphorus_delta = phosphorus * autotroph_fraction;
                outputs.autotroph_delta[autotroph_index].carbon_g_c += carbon_delta;
                outputs.autotroph_delta[autotroph_index].nitrogen_g_n += nitrogen_delta;
                outputs.autotroph_delta[autotroph_index].phosphorus_g_p += phosphorus_delta;
                allocated.carbon_g_c += carbon_delta;
                allocated.nitrogen_g_n += nitrogen_delta;
                allocated.phosphorus_g_p += phosphorus_delta;
            }
            try finitePool(allocated);
        }
    }
    return allocated;
}

fn validate(inputs: Inputs, outputs: Outputs) !void {
    try finiteNonnegativePool(inputs.amendment);
    const class_count = inputs.kinetic_fraction_by_class.len;
    const heterotroph_count = inputs.heterotroph_fraction_by_population.len;
    const autotroph_count = inputs.autotroph_fraction_by_population.len;
    if (class_count == 0 or heterotroph_count == 0 or autotroph_count == 0)
        return error.ZeroOrganicMicrobialAllocationExtent;
    const heterotroph_size = try std.math.mul(usize, heterotroph_count, class_count);
    const autotroph_size = try std.math.mul(usize, class_count, autotroph_count);
    if (inputs.nitrogen_to_carbon_by_population_class.len != heterotroph_size or
        inputs.phosphorus_to_carbon_by_population_class.len != heterotroph_size or
        outputs.heterotroph_delta.len != heterotroph_size or
        outputs.autotroph_delta.len != autotroph_size)
        return error.OrganicMicrobialAllocationDimensionMismatch;
    inline for (.{
        inputs.kinetic_fraction_by_class,
        inputs.heterotroph_fraction_by_population,
        inputs.nitrogen_to_carbon_by_population_class,
        inputs.phosphorus_to_carbon_by_population_class,
        inputs.autotroph_fraction_by_population,
    }) |values| for (values) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteOrganicMicrobialAllocationInput;
        if (value < 0) return error.InvalidOrganicMicrobialAllocationInput;
    };
}

fn finiteNonnegativePool(pool: ElementPool) !void {
    try finitePool(pool);
    inline for (@typeInfo(ElementPool).@"struct".fields) |field|
        if (@field(pool, field.name) < 0)
            return error.InvalidOrganicMicrobialAllocationInput;
}

fn finitePool(pool: ElementPool) !void {
    inline for (@typeInfo(ElementPool).@"struct".fields) |field|
        if (!std.math.isFinite(@field(pool, field.name)))
            return error.NonFiniteOrganicMicrobialAllocationResult;
}

test "runtime microbial allocation preserves N M NN source traversal" {
    var heterotroph: [4]ElementPool = undefined;
    var autotroph: [4]ElementPool = undefined;
    const allocated = try compute(.{
        .amendment = .{ .carbon_g_c = 100, .nitrogen_g_n = 20, .phosphorus_g_p = 10 },
        .kinetic_fraction_by_class = &.{ 0.1, 0.2 },
        .heterotroph_fraction_by_population = &.{ 0.5, 0.25 },
        .nitrogen_to_carbon_by_population_class = &.{ 0.1, 0.2, 0.3, 0.4 },
        .phosphorus_to_carbon_by_population_class = &.{ 0.01, 0.02, 0.03, 0.04 },
        .autotroph_fraction_by_population = &.{ 0.1, 0.2 },
    }, .{ .heterotroph_delta = &heterotroph, .autotroph_delta = &autotroph });
    try std.testing.expectEqual(@as(f64, 5), heterotroph[0].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 10), heterotroph[1].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 2.5), heterotroph[2].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 5), heterotroph[3].carbon_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 0.75), autotroph[0].carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), autotroph[1].carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 29.25), allocated.carbon_g_c, 1e-15);
}

test "runtime microbial allocation supports nonlegacy extents" {
    var heterotroph: [5]ElementPool = undefined;
    var autotroph: [15]ElementPool = undefined;
    _ = try compute(.{
        .amendment = .{ .carbon_g_c = 10, .nitrogen_g_n = 2, .phosphorus_g_p = 1 },
        .kinetic_fraction_by_class = &.{ 0.01, 0.02, 0.03, 0.04, 0.05 },
        .heterotroph_fraction_by_population = &.{0.5},
        .nitrogen_to_carbon_by_population_class = &.{ 0.1, 0.1, 0.1, 0.1, 0.1 },
        .phosphorus_to_carbon_by_population_class = &.{ 0.01, 0.01, 0.01, 0.01, 0.01 },
        .autotroph_fraction_by_population = &.{ 0.1, 0.2, 0.3 },
    }, .{ .heterotroph_delta = &heterotroph, .autotroph_delta = &autotroph });
}

test "dimension mismatch leaves caller output values unchanged" {
    var heterotroph = [_]ElementPool{.{ .carbon_g_c = 41, .nitrogen_g_n = 0, .phosphorus_g_p = 0 }};
    var autotroph = [_]ElementPool{.{ .carbon_g_c = 42, .nitrogen_g_n = 0, .phosphorus_g_p = 0 }};
    try std.testing.expectError(
        error.OrganicMicrobialAllocationDimensionMismatch,
        compute(.{
            .amendment = .{ .carbon_g_c = 1, .nitrogen_g_n = 1, .phosphorus_g_p = 1 },
            .kinetic_fraction_by_class = &.{ 0.1, 0.2 },
            .heterotroph_fraction_by_population = &.{0.5},
            .nitrogen_to_carbon_by_population_class = &.{ 0.1, 0.1 },
            .phosphorus_to_carbon_by_population_class = &.{ 0.01, 0.01 },
            .autotroph_fraction_by_population = &.{0.1},
        }, .{ .heterotroph_delta = &heterotroph, .autotroph_delta = &autotroph }),
    );
    try std.testing.expectEqual(@as(f64, 41), heterotroph[0].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 42), autotroph[0].carbon_g_c);
}
