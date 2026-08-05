const std = @import("std");

/// Exact REDIST 8885--8896 order for each plant species NZ.
pub const NodulePool = enum(u8) {
    bacterial_carbon,
    bacterial_nitrogen,
    bacterial_phosphorus,
    nonstructural_carbon,
    nonstructural_nitrogen,
    nonstructural_phosphorus,
};

pub const nodule_pool_count = std.meta.fields(NodulePool).len;
pub const Unit = enum { g_c, g_n, g_p };

pub fn unit(pool: NodulePool) Unit {
    return switch (pool) {
        .bacterial_carbon, .nonstructural_carbon => .g_c,
        .bacterial_nitrogen, .nonstructural_nitrogen => .g_n,
        .bacterial_phosphorus, .nonstructural_phosphorus => .g_p,
    };
}

pub const Topology = struct {
    layer_count: usize,
    species_count: usize, // NP
    maximum_axis_count: usize,
    minimum_root_mass_g_c: []const f64, // ZEROP(NZ)
    active_root_mass_g_c: []const f64, // WTRTL

    fn rootIndex(self: Topology, species: usize, layer: usize, axis: usize) usize {
        return (species * self.layer_count + layer) * self.maximum_axis_count + axis;
    }

    fn noduleIndex(self: Topology, species: usize, layer: usize) usize {
        return species * self.layer_count + layer;
    }
};

pub const Pools = struct {
    topology: Topology,
    /// Pool-major over species/layer storage.
    amounts: []f64,

    fn index(self: Pools, pool: NodulePool, species: usize, layer: usize) usize {
        const nodule_storage_len = self.topology.species_count * self.topology.layer_count;
        return @intFromEnum(pool) * nodule_storage_len + self.topology.noduleIndex(species, layer);
    }
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

fn eligible(topology: Topology, species: usize, source_layer: usize, destination_layer: usize) bool {
    return topology.active_root_mass_g_c[topology.rootIndex(species, source_layer, 0)] >
        topology.minimum_root_mass_g_c[species] and
        topology.active_root_mass_g_c[topology.rootIndex(species, destination_layer, 0)] >
            topology.minimum_root_mass_g_c[species];
}

/// Direct translation of REDIST 8885--8899, closing the shared plant guards.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, pools: Pools) !void {
    const topology = pools.topology;
    const root_storage_len = topology.species_count * topology.layer_count * topology.maximum_axis_count;
    const nodule_storage_len = topology.species_count * topology.layer_count;
    if (topology.layer_count == 0 or (topology.species_count > 0 and topology.maximum_axis_count == 0) or
        source_layer >= topology.layer_count or destination_layer >= topology.layer_count or
        source_layer == destination_layer or topology.minimum_root_mass_g_c.len != topology.species_count or
        topology.active_root_mass_g_c.len != root_storage_len or pools.amounts.len != nodule_pool_count * nodule_storage_len)
        return error.PondRootNoduleTransferDimensionMismatch;
    if (!finiteSlice(topology.minimum_root_mass_g_c) or !finiteSlice(topology.active_root_mass_g_c) or
        !finiteSlice(pools.amounts) or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondRootNoduleTransferInput;
    if (source_layer == 0) return;

    for (0..topology.species_count) |species| {
        if (!eligible(topology, species, source_layer, destination_layer)) continue;
        inline for (std.meta.fields(NodulePool)) |field| {
            const pool: NodulePool = @enumFromInt(field.value);
            const source = pools.index(pool, species, source_layer);
            const destination = pools.index(pool, species, destination_layer);
            if (!std.math.isFinite(pools.amounts[destination] + fraction * pools.amounts[source]))
                return error.NonFinitePondRootNoduleTransferResult;
        }
    }
    for (0..topology.species_count) |species| {
        if (!eligible(topology, species, source_layer, destination_layer)) continue;
        inline for (std.meta.fields(NodulePool)) |field| {
            const pool: NodulePool = @enumFromInt(field.value);
            const source = pools.index(pool, species, source_layer);
            const destination = pools.index(pool, species, destination_layer);
            pools.amounts[destination] += fraction * pools.amounts[source];
        }
    }
}

test "REDIST root nodules transfer per eligible runtime species" {
    const root_len = 2 * 3;
    const nodule_len = 2 * 3;
    var root_mass: [root_len]f64 = @splat(0);
    var amounts: [nodule_pool_count * nodule_len]f64 = @splat(0);
    const topology = Topology{ .layer_count = 3, .species_count = 2, .maximum_axis_count = 1, .minimum_root_mass_g_c = &.{ 0.1, 0.1 }, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    root_mass[topology.rootIndex(0, 1, 0)] = 1;
    root_mass[topology.rootIndex(0, 2, 0)] = 1;
    root_mass[topology.rootIndex(1, 1, 0)] = 0.1;
    root_mass[topology.rootIndex(1, 2, 0)] = 1;
    inline for (std.meta.fields(NodulePool)) |field| {
        const pool: NodulePool = @enumFromInt(field.value);
        amounts[pools.index(pool, 0, 1)] = 2;
        amounts[pools.index(pool, 0, 2)] = 1;
        amounts[pools.index(pool, 1, 2)] = 3;
    }
    try transfer(1, 2, 0.25, pools);
    inline for (std.meta.fields(NodulePool)) |field| {
        const pool: NodulePool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 1.5), amounts[pools.index(pool, 0, 2)]);
        try std.testing.expectEqual(@as(f64, 3), amounts[pools.index(pool, 1, 2)]);
    }
}

test "REDIST root nodule units preserve C N P order" {
    try std.testing.expectEqual(Unit.g_c, unit(.bacterial_carbon));
    try std.testing.expectEqual(Unit.g_n, unit(.bacterial_nitrogen));
    try std.testing.expectEqual(Unit.g_p, unit(.nonstructural_phosphorus));
}

test "REDIST root nodule source layer zero is excluded" {
    var root_mass: [2]f64 = @splat(1);
    var amounts: [nodule_pool_count * 2]f64 = @splat(1);
    const topology = Topology{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 1, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try transfer(0, 1, 0.5, .{ .topology = topology, .amounts = &amounts });
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST root nodule runtime topology is validated" {
    var root_mass: [2]f64 = @splat(1);
    var amounts: [nodule_pool_count * 2 - 1]f64 = @splat(1);
    const topology = Topology{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 1, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try std.testing.expectError(error.PondRootNoduleTransferDimensionMismatch, transfer(0, 1, 0.5, .{ .topology = topology, .amounts = &amounts }));
}

test "REDIST root nodule validation is atomic" {
    var root_mass: [3]f64 = .{ 0, 1, 1 };
    var amounts: [nodule_pool_count * 3]f64 = @splat(1);
    const topology = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    amounts[pools.index(.nonstructural_phosphorus, 0, 1)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondRootNoduleTransferInput, transfer(1, 2, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.bacterial_carbon, 0, 2)]);
}
