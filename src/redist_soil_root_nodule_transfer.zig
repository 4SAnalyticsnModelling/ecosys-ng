const std = @import("std");
const definitions = @import("redist_pond_root_nodule_transfer.zig");

pub const NodulePool = definitions.NodulePool;
pub const Unit = definitions.Unit;
pub const unit = definitions.unit;
pub const Topology = definitions.Topology;
pub const Pools = definitions.Pools;
pub const nodule_pool_count = definitions.nodule_pool_count;

fn rootIndex(t: Topology, species: usize, layer: usize, axis: usize) usize {
    return (species * t.layer_count + layer) * t.maximum_axis_count + axis;
}
fn poolIndex(p: Pools, pool: NodulePool, species: usize, layer: usize) usize {
    const species_layer_len = p.topology.species_count * p.topology.layer_count;
    return @intFromEnum(pool) * species_layer_len + species * p.topology.layer_count + layer;
}
fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 10535--10555, closing the plant and organic guards.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, minimum_volume_m3: f64, layer_volume_m3: []const f64, pools: Pools) !void {
    const t = pools.topology;
    const root_len = t.species_count * t.layer_count * t.maximum_axis_count;
    const nodule_len = t.species_count * t.layer_count;
    if (t.layer_count == 0 or (t.species_count > 0 and t.maximum_axis_count == 0) or source_layer >= t.layer_count or
        destination_layer >= t.layer_count or source_layer == destination_layer or layer_volume_m3.len != t.layer_count or
        t.minimum_root_mass_g_c.len != t.species_count or t.active_root_mass_g_c.len != root_len or
        pools.amounts.len != nodule_pool_count * nodule_len) return error.SoilRootNoduleTransferDimensionMismatch;
    if (!finiteSlice(layer_volume_m3) or !finiteSlice(t.minimum_root_mass_g_c) or !finiteSlice(t.active_root_mass_g_c) or
        !finiteSlice(pools.amounts) or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1 or
        !std.math.isFinite(minimum_volume_m3) or minimum_volume_m3 < 0)
        return error.InvalidSoilRootNoduleTransferInput;
    for (layer_volume_m3) |volume_m3| if (volume_m3 < 0) return error.InvalidSoilRootNoduleTransferInput;
    for (t.minimum_root_mass_g_c) |mass_g_c| if (mass_g_c < 0) return error.InvalidSoilRootNoduleTransferInput;
    for (t.active_root_mass_g_c) |mass_g_c| if (mass_g_c < 0) return error.InvalidSoilRootNoduleTransferInput;
    if (source_layer == 0 or layer_volume_m3[source_layer] <= minimum_volume_m3 or layer_volume_m3[destination_layer] <= minimum_volume_m3) return;

    for (0..t.species_count) |species| {
        const eligible = t.active_root_mass_g_c[rootIndex(t, species, source_layer, 0)] > t.minimum_root_mass_g_c[species] and
            t.active_root_mass_g_c[rootIndex(t, species, destination_layer, 0)] > t.minimum_root_mass_g_c[species];
        if (!eligible) continue;
        inline for (std.meta.fields(NodulePool)) |field| {
            const pool: NodulePool = @enumFromInt(field.value);
            const source = poolIndex(pools, pool, species, source_layer);
            const destination = poolIndex(pools, pool, species, destination_layer);
            const moved = fraction * pools.amounts[source];
            if (!std.math.isFinite(moved) or !std.math.isFinite(pools.amounts[destination] + moved) or
                !std.math.isFinite(pools.amounts[source] - moved)) return error.NonFiniteSoilRootNoduleTransferResult;
        }
    }
    for (0..t.species_count) |species| {
        const eligible = t.active_root_mass_g_c[rootIndex(t, species, source_layer, 0)] > t.minimum_root_mass_g_c[species] and
            t.active_root_mass_g_c[rootIndex(t, species, destination_layer, 0)] > t.minimum_root_mass_g_c[species];
        if (!eligible) continue;
        inline for (std.meta.fields(NodulePool)) |field| {
            const pool: NodulePool = @enumFromInt(field.value);
            const source = poolIndex(pools, pool, species, source_layer);
            const destination = poolIndex(pools, pool, species, destination_layer);
            const moved = fraction * pools.amounts[source];
            pools.amounts[destination] = pools.amounts[destination] + moved;
            pools.amounts[source] = pools.amounts[source] - moved;
        }
    }
}

test "REDIST soil root nodules support arbitrary runtime species" {
    const root_len = 2 * 3;
    const nodule_len = 2 * 3;
    var root_mass: [root_len]f64 = @splat(0);
    var amounts: [nodule_pool_count * nodule_len]f64 = @splat(0);
    const t = Topology{ .layer_count = 3, .species_count = 2, .maximum_axis_count = 1, .minimum_root_mass_g_c = &.{ 0.1, 0.1 }, .active_root_mass_g_c = &root_mass };
    const p = Pools{ .topology = t, .amounts = &amounts };
    const volumes = [_]f64{ 1, 1, 1 };
    root_mass[rootIndex(t, 0, 1, 0)] = 1;
    root_mass[rootIndex(t, 0, 2, 0)] = 1;
    root_mass[rootIndex(t, 1, 1, 0)] = 0.1;
    root_mass[rootIndex(t, 1, 2, 0)] = 1;
    inline for (std.meta.fields(NodulePool)) |field| {
        const pool: NodulePool = @enumFromInt(field.value);
        amounts[poolIndex(p, pool, 0, 1)] = 4;
        amounts[poolIndex(p, pool, 0, 2)] = 2;
        amounts[poolIndex(p, pool, 1, 1)] = 8;
    }
    try transfer(1, 2, 0.25, 0, &volumes, p);
    inline for (std.meta.fields(NodulePool)) |field| {
        const pool: NodulePool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 3), amounts[poolIndex(p, pool, 0, 1)]);
        try std.testing.expectEqual(@as(f64, 3), amounts[poolIndex(p, pool, 0, 2)]);
        try std.testing.expectEqual(@as(f64, 8), amounts[poolIndex(p, pool, 1, 1)]);
    }
}
test "REDIST soil root nodules conserve C N and P" {
    const len = 1 * 3;
    var root_mass: [len]f64 = @splat(1);
    var amounts: [nodule_pool_count * len]f64 = @splat(2);
    const t = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .minimum_root_mass_g_c = &.{0}, .active_root_mass_g_c = &root_mass };
    const p = Pools{ .topology = t, .amounts = &amounts };
    const v = [_]f64{ 1, 1, 1 };
    try transfer(1, 2, 0.4, 0, &v, p);
    inline for (std.meta.fields(NodulePool)) |field| {
        const pool: NodulePool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 4), amounts[poolIndex(p, pool, 0, 1)] + amounts[poolIndex(p, pool, 0, 2)]);
    }
}
test "REDIST soil root nodule validation is atomic" {
    const len = 1 * 3;
    var root_mass: [len]f64 = @splat(1);
    var amounts: [nodule_pool_count * len]f64 = @splat(1);
    const t = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .minimum_root_mass_g_c = &.{0}, .active_root_mass_g_c = &root_mass };
    const p = Pools{ .topology = t, .amounts = &amounts };
    const v = [_]f64{ 1, 1, 1 };
    amounts[poolIndex(p, .nonstructural_phosphorus, 0, 1)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilRootNoduleTransferInput, transfer(1, 2, 0.5, 0, &v, p));
    try std.testing.expectEqual(@as(f64, 1), amounts[poolIndex(p, .bacterial_carbon, 0, 2)]);
}
