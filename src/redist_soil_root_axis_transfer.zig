const std = @import("std");
const definitions = @import("redist_pond_root_axis_transfer.zig");

pub const RootAxisPool = definitions.RootAxisPool;
pub const Unit = definitions.Unit;
pub const unit = definitions.unit;
pub const Topology = definitions.Topology;
pub const Pools = definitions.Pools;
pub const root_axis_pool_count = definitions.root_axis_pool_count;

fn rootIndex(t: Topology, species: usize, layer: usize, axis: usize) usize {
    return (species * t.layer_count + layer) * t.maximum_axis_count + axis;
}
fn poolIndex(p: Pools, pool: RootAxisPool, species: usize, layer: usize, axis: usize) usize {
    const topology_len = p.topology.species_count * p.topology.layer_count * p.topology.maximum_axis_count;
    return @intFromEnum(pool) * topology_len + rootIndex(p.topology, species, layer, axis);
}
fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 10483--10531 under the shared plant guards.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, minimum_volume_m3: f64, layer_volume_m3: []const f64, pools: Pools) !void {
    const t = pools.topology;
    const topology_len = t.species_count * t.layer_count * t.maximum_axis_count;
    if (t.layer_count == 0 or (t.species_count > 0 and t.maximum_axis_count == 0) or source_layer >= t.layer_count or
        destination_layer >= t.layer_count or source_layer == destination_layer or layer_volume_m3.len != t.layer_count or
        t.active_axis_count.len != t.species_count or t.minimum_root_mass_g_c.len != t.species_count or
        t.active_root_mass_g_c.len != topology_len or pools.amounts.len != root_axis_pool_count * topology_len)
        return error.SoilRootAxisTransferDimensionMismatch;
    for (t.active_axis_count) |count| if (count > t.maximum_axis_count) return error.SoilRootAxisTransferDimensionMismatch;
    if (!finiteSlice(layer_volume_m3) or !finiteSlice(t.minimum_root_mass_g_c) or !finiteSlice(t.active_root_mass_g_c) or
        !finiteSlice(pools.amounts) or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1 or
        !std.math.isFinite(minimum_volume_m3) or minimum_volume_m3 < 0)
        return error.InvalidSoilRootAxisTransferInput;
    for (layer_volume_m3) |volume_m3| if (volume_m3 < 0) return error.InvalidSoilRootAxisTransferInput;
    for (t.minimum_root_mass_g_c) |mass_g_c| if (mass_g_c < 0) return error.InvalidSoilRootAxisTransferInput;
    for (t.active_root_mass_g_c) |mass_g_c| if (mass_g_c < 0) return error.InvalidSoilRootAxisTransferInput;
    if (source_layer == 0 or layer_volume_m3[source_layer] <= minimum_volume_m3 or layer_volume_m3[destination_layer] <= minimum_volume_m3) return;

    for (0..t.species_count) |species| {
        const eligible = t.active_root_mass_g_c[rootIndex(t, species, source_layer, 0)] > t.minimum_root_mass_g_c[species] and
            t.active_root_mass_g_c[rootIndex(t, species, destination_layer, 0)] > t.minimum_root_mass_g_c[species];
        if (!eligible) continue;
        for (0..t.active_axis_count[species]) |axis| inline for (std.meta.fields(RootAxisPool)) |field| {
            const pool: RootAxisPool = @enumFromInt(field.value);
            const source = poolIndex(pools, pool, species, source_layer, axis);
            const destination = poolIndex(pools, pool, species, destination_layer, axis);
            const moved = fraction * pools.amounts[source];
            if (!std.math.isFinite(moved) or !std.math.isFinite(pools.amounts[destination] + moved) or
                !std.math.isFinite(pools.amounts[source] - moved)) return error.NonFiniteSoilRootAxisTransferResult;
        };
    }
    for (0..t.species_count) |species| {
        const eligible = t.active_root_mass_g_c[rootIndex(t, species, source_layer, 0)] > t.minimum_root_mass_g_c[species] and
            t.active_root_mass_g_c[rootIndex(t, species, destination_layer, 0)] > t.minimum_root_mass_g_c[species];
        if (!eligible) continue;
        for (0..t.active_axis_count[species]) |axis| inline for (std.meta.fields(RootAxisPool)) |field| {
            const pool: RootAxisPool = @enumFromInt(field.value);
            const source = poolIndex(pools, pool, species, source_layer, axis);
            const destination = poolIndex(pools, pool, species, destination_layer, axis);
            const moved = fraction * pools.amounts[source];
            pools.amounts[destination] = pools.amounts[destination] + moved;
            pools.amounts[source] = pools.amounts[source] - moved;
        };
    }
}

test "REDIST soil root axes support arbitrary runtime species and MY" {
    const len = 2 * 3 * 2;
    var root_mass: [len]f64 = @splat(0);
    var amounts: [root_axis_pool_count * len]f64 = @splat(0);
    const t = Topology{ .layer_count = 3, .species_count = 2, .maximum_axis_count = 2, .active_axis_count = &.{ 2, 1 }, .minimum_root_mass_g_c = &.{ 0.1, 0.1 }, .active_root_mass_g_c = &root_mass };
    const p = Pools{ .topology = t, .amounts = &amounts };
    const volumes = [_]f64{ 1, 1, 1 };
    root_mass[rootIndex(t, 0, 1, 0)] = 1;
    root_mass[rootIndex(t, 0, 2, 0)] = 1;
    root_mass[rootIndex(t, 1, 1, 0)] = 0.1;
    root_mass[rootIndex(t, 1, 2, 0)] = 1;
    inline for (std.meta.fields(RootAxisPool)) |field| {
        const pool: RootAxisPool = @enumFromInt(field.value);
        amounts[poolIndex(p, pool, 0, 1, 1)] = 4;
        amounts[poolIndex(p, pool, 0, 2, 1)] = 2;
        amounts[poolIndex(p, pool, 1, 1, 0)] = 8;
    }
    try transfer(1, 2, 0.25, 0, &volumes, p);
    inline for (std.meta.fields(RootAxisPool)) |field| {
        const pool: RootAxisPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 3), amounts[poolIndex(p, pool, 0, 1, 1)]);
        try std.testing.expectEqual(@as(f64, 3), amounts[poolIndex(p, pool, 0, 2, 1)]);
        try std.testing.expectEqual(@as(f64, 8), amounts[poolIndex(p, pool, 1, 1, 0)]);
    }
}
test "REDIST soil root axis elemental pools conserve C N and P" {
    const len = 1 * 3 * 1;
    var root_mass: [len]f64 = @splat(1);
    var amounts: [root_axis_pool_count * len]f64 = @splat(2);
    const t = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .active_axis_count = &.{1}, .minimum_root_mass_g_c = &.{0}, .active_root_mass_g_c = &root_mass };
    const p = Pools{ .topology = t, .amounts = &amounts };
    const v = [_]f64{ 1, 1, 1 };
    try transfer(1, 2, 0.4, 0, &v, p);
    inline for (.{ RootAxisPool.nonstructural_carbon, RootAxisPool.nonstructural_nitrogen, RootAxisPool.nonstructural_phosphorus }) |pool| try std.testing.expectEqual(@as(f64, 4), amounts[poolIndex(p, pool, 0, 1, 0)] + amounts[poolIndex(p, pool, 0, 2, 0)]);
}
test "REDIST soil root axis validation is atomic" {
    const len = 1 * 3 * 1;
    var root_mass: [len]f64 = @splat(1);
    var amounts: [root_axis_pool_count * len]f64 = @splat(1);
    const t = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .active_axis_count = &.{1}, .minimum_root_mass_g_c = &.{0}, .active_root_mass_g_c = &root_mass };
    const p = Pools{ .topology = t, .amounts = &amounts };
    const v = [_]f64{ 1, 1, 1 };
    amounts[poolIndex(p, .average_secondary_length, 0, 1, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilRootAxisTransferInput, transfer(1, 2, 0.5, 0, &v, p));
    try std.testing.expectEqual(@as(f64, 1), amounts[poolIndex(p, .nonstructural_carbon, 0, 2, 0)]);
}
