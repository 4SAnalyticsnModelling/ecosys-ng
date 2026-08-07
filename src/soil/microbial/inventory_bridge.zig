const std = @import("std");
const microbial = @import("state.zig");
const organic = @import("../organic/initialization.zig");

/// Publishes runtime microbial biomass into the fixed soil-organic inventory
/// used by mass-balance reconstruction. Runtime dimensions larger than the
/// representable legacy role axes fail before any mirror pool is changed.
pub fn publishToOrganic(
    source: *const microbial.State,
    destination: *organic.State,
) !void {
    if (source.cell_count == 0 or source.layer_count == 0 or
        source.cell_count * source.layer_count != destination.layer_count)
        return error.SoilMicrobialInventoryBridgeDimensionMismatch;
    if (source.substrate_count > organic.microbial_substrate_count or
        source.population_count > organic.microbial_population_count)
        return error.SoilMicrobialInventoryExceedsOrganicMirror;

    const runtime_pool_count = try std.math.mul(
        usize,
        try std.math.mul(
            usize,
            destination.layer_count,
            source.substrate_count,
        ),
        source.population_count,
    );
    if (source.nonstructural.len != runtime_pool_count or
        source.structural.len != try std.math.mul(usize, runtime_pool_count, 2))
        return error.SoilMicrobialInventoryBridgeDimensionMismatch;

    for (source.nonstructural) |pool| try validatePool(pool);
    for (source.structural) |pool| try validatePool(pool);

    for (0..destination.layer_count) |layer| {
        for (0..source.substrate_count) |substrate| {
            for (0..source.population_count) |population| {
                const runtime_index = try source.populationIndex(
                    layer / source.layer_count,
                    layer % source.layer_count,
                    substrate,
                    population,
                );
                const mirror_index =
                    ((layer * organic.microbial_substrate_count + substrate) *
                        organic.microbial_population_count + population) *
                    organic.kinetic_fraction_count;
                destination.microbial[mirror_index] = toOrganic(
                    source.structural[runtime_index * 2],
                );
                destination.microbial[mirror_index + 1] = toOrganic(
                    source.structural[runtime_index * 2 + 1],
                );
                destination.microbial[mirror_index + 2] = toOrganic(
                    source.nonstructural[runtime_index],
                );
            }
        }
    }
}

/// Publishes the fixed soil-organic microbial mirror back into runtime
/// microbial biomass after processes such as pond transfer modify the mirror.
/// The complete mapped source is validated before runtime state is changed.
pub fn publishFromOrganic(
    source: *const organic.State,
    destination: *microbial.State,
) !void {
    try validateMappingDimensions(destination, source);

    for (0..source.layer_count) |layer| {
        for (0..destination.substrate_count) |substrate| {
            for (0..destination.population_count) |population| {
                const first =
                    ((layer * organic.microbial_substrate_count + substrate) *
                        organic.microbial_population_count + population) *
                    organic.kinetic_fraction_count;
                for (source.microbial[first .. first + organic.kinetic_fraction_count]) |pool|
                    try validateOrganicPool(pool);
            }
        }
    }

    for (0..source.layer_count) |layer| {
        for (0..destination.substrate_count) |substrate| {
            for (0..destination.population_count) |population| {
                const runtime_index = try destination.populationIndex(
                    layer / destination.layer_count,
                    layer % destination.layer_count,
                    substrate,
                    population,
                );
                const first =
                    ((layer * organic.microbial_substrate_count + substrate) *
                        organic.microbial_population_count + population) *
                    organic.kinetic_fraction_count;
                destination.structural[runtime_index * 2] = toRuntime(
                    source.microbial[first],
                );
                destination.structural[runtime_index * 2 + 1] = toRuntime(
                    source.microbial[first + 1],
                );
                destination.nonstructural[runtime_index] = toRuntime(
                    source.microbial[first + 2],
                );
            }
        }
    }
}

fn validateMappingDimensions(
    runtime: *const microbial.State,
    mirror: *const organic.State,
) !void {
    if (runtime.cell_count == 0 or runtime.layer_count == 0 or
        runtime.cell_count * runtime.layer_count != mirror.layer_count)
        return error.SoilMicrobialInventoryBridgeDimensionMismatch;
    if (runtime.substrate_count > organic.microbial_substrate_count or
        runtime.population_count > organic.microbial_population_count)
        return error.SoilMicrobialInventoryExceedsOrganicMirror;
    const runtime_pool_count = try std.math.mul(
        usize,
        try std.math.mul(usize, mirror.layer_count, runtime.substrate_count),
        runtime.population_count,
    );
    if (runtime.nonstructural.len != runtime_pool_count or
        runtime.structural.len != try std.math.mul(usize, runtime_pool_count, 2))
        return error.SoilMicrobialInventoryBridgeDimensionMismatch;
}

fn validatePool(pool: microbial.ElementalPool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoilMicrobialInventoryBridgePool;
}

fn toOrganic(pool: microbial.ElementalPool) organic.ElementPool {
    return .{
        .carbon_g_c = pool.carbon_g_c,
        .nitrogen_g_n = pool.nitrogen_g_n,
        .phosphorus_g_p = pool.phosphorus_g_p,
    };
}

fn validateOrganicPool(pool: organic.ElementPool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoilMicrobialInventoryBridgePool;
}

fn toRuntime(pool: organic.ElementPool) microbial.ElementalPool {
    return .{
        .carbon_g_c = pool.carbon_g_c,
        .nitrogen_g_n = pool.nitrogen_g_n,
        .phosphorus_g_p = pool.phosphorus_g_p,
    };
}

const ElementTotals = struct {
    carbon_g_c: f64 = 0,
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
};

fn runtimeTotals(state: *const microbial.State) ElementTotals {
    var totals: ElementTotals = .{};
    for (state.structural) |pool| addPool(&totals, pool);
    for (state.nonstructural) |pool| addPool(&totals, pool);
    return totals;
}

fn mappedOrganicTotals(
    state: *const organic.State,
    substrate_count: usize,
    population_count: usize,
) ElementTotals {
    var totals: ElementTotals = .{};
    for (0..state.layer_count) |layer| for (0..substrate_count) |substrate| for (0..population_count) |population| {
        const first = ((layer * organic.microbial_substrate_count + substrate) *
            organic.microbial_population_count + population) *
            organic.kinetic_fraction_count;
        for (state.microbial[first .. first + organic.kinetic_fraction_count]) |pool|
            addPool(&totals, pool);
    };
    return totals;
}

fn addPool(totals: *ElementTotals, pool: anytype) void {
    totals.carbon_g_c += pool.carbon_g_c;
    totals.nitrogen_g_n += pool.nitrogen_g_n;
    totals.phosphorus_g_p += pool.phosphorus_g_p;
}

test "microbial inventory publication preserves mapped C N P and is idempotent" {
    var source = try microbial.State.init(std.testing.allocator, 1, 2, 2, 3);
    defer source.deinit();
    var destination = try organic.State.init(std.testing.allocator, 2);
    defer destination.deinit();

    for (source.nonstructural, 0..) |*pool, index| pool.* = .{
        .carbon_g_c = @floatFromInt(index + 1),
        .nitrogen_g_n = @as(f64, @floatFromInt(index + 1)) / 10,
        .phosphorus_g_p = @as(f64, @floatFromInt(index + 1)) / 100,
    };
    for (source.structural, 0..) |*pool, index| pool.* = .{
        .carbon_g_c = @floatFromInt(index + 21),
        .nitrogen_g_n = @as(f64, @floatFromInt(index + 21)) / 10,
        .phosphorus_g_p = @as(f64, @floatFromInt(index + 21)) / 100,
    };
    @memset(destination.microbial, .{
        .carbon_g_c = 999,
        .nitrogen_g_n = 99,
        .phosphorus_g_p = 9,
    });
    const unmapped_index =
        ((organic.microbial_substrate_count - 1) *
            organic.microbial_population_count) *
        organic.kinetic_fraction_count;

    try publishToOrganic(&source, &destination);
    const runtime_totals = runtimeTotals(&source);
    const mirror_totals = mappedOrganicTotals(&destination, 2, 3);
    try std.testing.expectApproxEqAbs(
        runtime_totals.carbon_g_c,
        mirror_totals.carbon_g_c,
        1.0e-12,
    );
    try std.testing.expectApproxEqAbs(
        runtime_totals.nitrogen_g_n,
        mirror_totals.nitrogen_g_n,
        1.0e-13,
    );
    try std.testing.expectApproxEqAbs(
        runtime_totals.phosphorus_g_p,
        mirror_totals.phosphorus_g_p,
        1.0e-14,
    );
    try std.testing.expectEqual(@as(f64, 999), destination.microbial[unmapped_index].carbon_g_c);
    const first_publication = try std.testing.allocator.dupe(
        organic.ElementPool,
        destination.microbial,
    );
    defer std.testing.allocator.free(first_publication);
    try publishToOrganic(&source, &destination);
    try std.testing.expectEqualSlices(
        organic.ElementPool,
        first_publication,
        destination.microbial,
    );
}

test "oversized runtime microbial axes fail before mirror mutation" {
    var source = try microbial.State.init(
        std.testing.allocator,
        1,
        1,
        organic.microbial_substrate_count + 1,
        organic.microbial_population_count,
    );
    defer source.deinit();
    var destination = try organic.State.init(std.testing.allocator, 1);
    defer destination.deinit();
    destination.microbial[0].carbon_g_c = 17;

    try std.testing.expectError(
        error.SoilMicrobialInventoryExceedsOrganicMirror,
        publishToOrganic(&source, &destination),
    );
    try std.testing.expectEqual(@as(f64, 17), destination.microbial[0].carbon_g_c);
}

test "invalid late runtime pool fails before mirror mutation" {
    var source = try microbial.State.init(std.testing.allocator, 1, 1, 1, 1);
    defer source.deinit();
    var destination = try organic.State.init(std.testing.allocator, 1);
    defer destination.deinit();
    destination.microbial[0].carbon_g_c = 23;
    source.structural[1].carbon_g_c = std.math.nan(f64);

    try std.testing.expectError(
        error.InvalidSoilMicrobialInventoryBridgePool,
        publishToOrganic(&source, &destination),
    );
    try std.testing.expectEqual(@as(f64, 23), destination.microbial[0].carbon_g_c);
}

test "reverse microbial publication preserves mapped C N P and is idempotent" {
    var source = try organic.State.init(std.testing.allocator, 2);
    defer source.deinit();
    var destination = try microbial.State.init(std.testing.allocator, 1, 2, 2, 3);
    defer destination.deinit();

    for (0..source.layer_count) |layer| for (0..destination.substrate_count) |substrate| for (0..destination.population_count) |population| {
        const first = ((layer * organic.microbial_substrate_count + substrate) *
            organic.microbial_population_count + population) *
            organic.kinetic_fraction_count;
        for (source.microbial[first .. first + organic.kinetic_fraction_count], 0..) |*pool, fraction| {
            const ordinal = 1 + first + fraction;
            pool.* = .{
                .carbon_g_c = @floatFromInt(ordinal),
                .nitrogen_g_n = @as(f64, @floatFromInt(ordinal)) / 10,
                .phosphorus_g_p = @as(f64, @floatFromInt(ordinal)) / 100,
            };
        }
    };
    @memset(destination.nonstructural, .{
        .carbon_g_c = 999,
        .nitrogen_g_n = 99,
        .phosphorus_g_p = 9,
    });
    @memset(destination.structural, .{
        .carbon_g_c = 999,
        .nitrogen_g_n = 99,
        .phosphorus_g_p = 9,
    });

    try publishFromOrganic(&source, &destination);
    const mirror_totals = mappedOrganicTotals(&source, 2, 3);
    const runtime_totals = runtimeTotals(&destination);
    try std.testing.expectApproxEqAbs(mirror_totals.carbon_g_c, runtime_totals.carbon_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(mirror_totals.nitrogen_g_n, runtime_totals.nitrogen_g_n, 1.0e-13);
    try std.testing.expectApproxEqAbs(mirror_totals.phosphorus_g_p, runtime_totals.phosphorus_g_p, 1.0e-14);
    const nonstructural = try std.testing.allocator.dupe(microbial.ElementalPool, destination.nonstructural);
    defer std.testing.allocator.free(nonstructural);
    const structural = try std.testing.allocator.dupe(microbial.ElementalPool, destination.structural);
    defer std.testing.allocator.free(structural);
    try publishFromOrganic(&source, &destination);
    try std.testing.expectEqualSlices(microbial.ElementalPool, nonstructural, destination.nonstructural);
    try std.testing.expectEqualSlices(microbial.ElementalPool, structural, destination.structural);
}

test "reverse publication rejects oversized runtime axes atomically" {
    var source = try organic.State.init(std.testing.allocator, 1);
    defer source.deinit();
    var destination = try microbial.State.init(
        std.testing.allocator,
        1,
        1,
        organic.microbial_substrate_count,
        organic.microbial_population_count + 1,
    );
    defer destination.deinit();
    destination.nonstructural[0].carbon_g_c = 29;

    try std.testing.expectError(
        error.SoilMicrobialInventoryExceedsOrganicMirror,
        publishFromOrganic(&source, &destination),
    );
    try std.testing.expectEqual(@as(f64, 29), destination.nonstructural[0].carbon_g_c);
}

test "invalid late organic mirror pool leaves runtime inventory unchanged" {
    var source = try organic.State.init(std.testing.allocator, 1);
    defer source.deinit();
    var destination = try microbial.State.init(std.testing.allocator, 1, 1, 1, 1);
    defer destination.deinit();
    destination.nonstructural[0].carbon_g_c = 31;
    source.microbial[2].phosphorus_g_p = std.math.nan(f64);

    try std.testing.expectError(
        error.InvalidSoilMicrobialInventoryBridgePool,
        publishFromOrganic(&source, &destination),
    );
    try std.testing.expectEqual(@as(f64, 31), destination.nonstructural[0].carbon_g_c);
}
