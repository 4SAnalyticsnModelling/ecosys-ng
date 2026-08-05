const std = @import("std");
const destination_transfer = @import("redist_pond_root_nodule_transfer.zig");

pub const NodulePool = destination_transfer.NodulePool;
pub const Topology = destination_transfer.Topology;
pub const Pools = destination_transfer.Pools;
pub const nodule_pool_count = destination_transfer.nodule_pool_count;
pub const unit = destination_transfer.unit;

fn rootIndex(t: Topology, species: usize, layer: usize, axis: usize) usize {
    return (species * t.layer_count + layer) * t.maximum_axis_count + axis;
}
fn noduleIndex(t: Topology, species: usize, layer: usize) usize {
    return species * t.layer_count + layer;
}
fn poolIndex(pools: Pools, pool: NodulePool, species: usize, layer: usize) usize {
    const nodule_len = pools.topology.species_count * pools.topology.layer_count;
    return @intFromEnum(pool) * nodule_len + noduleIndex(pools.topology, species, layer);
}
fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}
fn eligible(t: Topology, species: usize, source: usize, destination: usize) bool {
    return t.active_root_mass_g_c[rootIndex(t, species, source, 0)] > t.minimum_root_mass_g_c[species] and
        t.active_root_mass_g_c[rootIndex(t, species, destination, 0)] > t.minimum_root_mass_g_c[species];
}

/// Direct translation of REDIST 9284--9292, closing the source plant guards.
pub fn reduce(source_layer: usize, destination_layer: usize, remaining_fraction: f64, pools: Pools) !void {
    const t = pools.topology;
    const root_len = t.species_count * t.layer_count * t.maximum_axis_count;
    const nodule_len = t.species_count * t.layer_count;
    if (t.layer_count == 0 or (t.species_count > 0 and t.maximum_axis_count == 0) or
        source_layer >= t.layer_count or destination_layer >= t.layer_count or source_layer == destination_layer or
        t.minimum_root_mass_g_c.len != t.species_count or t.active_root_mass_g_c.len != root_len or
        pools.amounts.len != nodule_pool_count * nodule_len)
        return error.PondSourceRootNoduleReductionDimensionMismatch;
    if (!finiteSlice(t.minimum_root_mass_g_c) or !finiteSlice(t.active_root_mass_g_c) or
        !finiteSlice(pools.amounts) or !std.math.isFinite(remaining_fraction) or
        remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceRootNoduleReductionInput;
    if (source_layer == 0) return;

    for (0..t.species_count) |species| {
        if (!eligible(t, species, source_layer, destination_layer)) continue;
        inline for (std.meta.fields(NodulePool)) |field| {
            const pool: NodulePool = @enumFromInt(field.value);
            const index = poolIndex(pools, pool, species, source_layer);
            if (!std.math.isFinite(remaining_fraction * pools.amounts[index]))
                return error.NonFinitePondSourceRootNoduleReductionResult;
        }
    }
    for (0..t.species_count) |species| {
        if (!eligible(t, species, source_layer, destination_layer)) continue;
        inline for (std.meta.fields(NodulePool)) |field| {
            const pool: NodulePool = @enumFromInt(field.value);
            const index = poolIndex(pools, pool, species, source_layer);
            pools.amounts[index] = remaining_fraction * pools.amounts[index];
        }
    }
}

test "REDIST source nodules scale exact pools for eligible runtime species" {
    const root_len = 2 * 3;
    const nodule_len = 2 * 3;
    var root_mass: [root_len]f64 = @splat(0);
    var amounts: [nodule_pool_count * nodule_len]f64 = @splat(0);
    const topology = Topology{ .layer_count = 3, .species_count = 2, .maximum_axis_count = 1, .minimum_root_mass_g_c = &.{ 0.1, 0.1 }, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    root_mass[rootIndex(topology, 0, 1, 0)] = 1;
    root_mass[rootIndex(topology, 0, 2, 0)] = 1;
    root_mass[rootIndex(topology, 1, 1, 0)] = 0.1;
    root_mass[rootIndex(topology, 1, 2, 0)] = 1;
    inline for (std.meta.fields(NodulePool)) |field| {
        const pool: NodulePool = @enumFromInt(field.value);
        amounts[poolIndex(pools, pool, 0, 1)] = 4;
        amounts[poolIndex(pools, pool, 1, 1)] = 8;
    }
    try reduce(1, 2, 0.25, pools);
    inline for (std.meta.fields(NodulePool)) |field| {
        const pool: NodulePool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 1), amounts[poolIndex(pools, pool, 0, 1)]);
        try std.testing.expectEqual(@as(f64, 8), amounts[poolIndex(pools, pool, 1, 1)]);
    }
}

test "REDIST source nodule NP zero is a zero-trip loop" {
    var root_mass: [0]f64 = .{};
    var amounts: [0]f64 = .{};
    const topology = Topology{ .layer_count = 3, .species_count = 0, .maximum_axis_count = 0, .minimum_root_mass_g_c = &.{}, .active_root_mass_g_c = &root_mass };
    try reduce(1, 2, 0, .{ .topology = topology, .amounts = &amounts });
}

test "REDIST source nodule layer zero is excluded" {
    var root_mass: [2]f64 = @splat(1);
    var amounts: [nodule_pool_count * 2]f64 = @splat(1);
    const topology = Topology{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 1, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try reduce(0, 1, 0, .{ .topology = topology, .amounts = &amounts });
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST source nodule validation is atomic" {
    var root_mass: [3]f64 = .{ 0, 1, 1 };
    var amounts: [nodule_pool_count * 3]f64 = @splat(1);
    const topology = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    amounts[poolIndex(pools, .nonstructural_phosphorus, 0, 1)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceRootNoduleReductionInput, reduce(1, 2, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[poolIndex(pools, .bacterial_carbon, 0, 1)]);
}
