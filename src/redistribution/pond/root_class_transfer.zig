const std = @import("std");

/// Exact REDIST 8830--8847 order for each root class NR.
pub const RootPool = enum(u8) {
    primary_root_carbon,
    primary_root_nitrogen,
    primary_root_phosphorus,
    secondary_root_carbon,
    secondary_root_nitrogen,
    secondary_root_phosphorus,
    primary_root_length,
    secondary_root_length,
    secondary_root_axis_count,
};

pub const root_pool_count = std.meta.fields(RootPool).len;
pub const Unit = enum { g_c, g_n, g_p, m, count };

pub fn unit(pool: RootPool) Unit {
    return switch (pool) {
        .primary_root_carbon, .secondary_root_carbon => .g_c,
        .primary_root_nitrogen, .secondary_root_nitrogen => .g_n,
        .primary_root_phosphorus, .secondary_root_phosphorus => .g_p,
        .primary_root_length, .secondary_root_length => .m,
        .secondary_root_axis_count => .count,
    };
}

pub const Topology = struct {
    layer_count: usize,
    species_count: usize, // NP
    maximum_axis_count: usize,
    maximum_root_class_count: usize,
    active_axis_count: []const usize, // MY(NZ)
    active_root_class_count: []const usize, // NRT(NZ)
    minimum_root_mass_g_c: []const f64, // ZEROP(NZ)
    /// Species-major, layer, axis. WTRTL.
    active_root_mass_g_c: []const f64,

    fn rootIndex(self: Topology, species: usize, layer: usize, axis: usize) usize {
        return (species * self.layer_count + layer) * self.maximum_axis_count + axis;
    }

    fn classIndex(self: Topology, species: usize, layer: usize, axis: usize, root_class: usize) usize {
        return self.rootIndex(species, layer, axis) * self.maximum_root_class_count + root_class;
    }
};

pub const Pools = struct {
    topology: Topology,
    /// Pool-major over species/layer/axis/root-class storage.
    amounts: []f64,

    fn index(self: Pools, pool: RootPool, species: usize, layer: usize, axis: usize, root_class: usize) usize {
        const class_storage_len = self.topology.species_count * self.topology.layer_count *
            self.topology.maximum_axis_count * self.topology.maximum_root_class_count;
        return @intFromEnum(pool) * class_storage_len +
            self.topology.classIndex(species, layer, axis, root_class);
    }
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8829--8848 under the existing plant guards.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, pools: Pools) !void {
    const topology = pools.topology;
    const root_storage_len = topology.species_count * topology.layer_count * topology.maximum_axis_count;
    const class_storage_len = root_storage_len * topology.maximum_root_class_count;
    if (topology.layer_count == 0 or (topology.species_count > 0 and topology.maximum_axis_count == 0) or
        source_layer >= topology.layer_count or
        destination_layer >= topology.layer_count or source_layer == destination_layer or
        topology.active_axis_count.len != topology.species_count or
        topology.active_root_class_count.len != topology.species_count or
        topology.minimum_root_mass_g_c.len != topology.species_count or
        topology.active_root_mass_g_c.len != root_storage_len or pools.amounts.len != root_pool_count * class_storage_len)
        return error.PondRootClassTransferDimensionMismatch;
    for (topology.active_axis_count) |count| if (count > topology.maximum_axis_count)
        return error.PondRootClassTransferDimensionMismatch;
    for (topology.active_root_class_count) |count| if (count > topology.maximum_root_class_count)
        return error.PondRootClassTransferDimensionMismatch;
    if (!finiteSlice(topology.minimum_root_mass_g_c) or !finiteSlice(topology.active_root_mass_g_c) or
        !finiteSlice(pools.amounts) or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondRootClassTransferInput;
    if (source_layer == 0) return;

    for (0..topology.species_count) |species| {
        const eligible = topology.active_root_mass_g_c[topology.rootIndex(species, source_layer, 0)] >
            topology.minimum_root_mass_g_c[species] and
            topology.active_root_mass_g_c[topology.rootIndex(species, destination_layer, 0)] >
                topology.minimum_root_mass_g_c[species];
        if (!eligible) continue;
        for (0..topology.active_axis_count[species]) |axis| {
            for (0..topology.active_root_class_count[species]) |root_class| {
                inline for (std.meta.fields(RootPool)) |field| {
                    const pool: RootPool = @enumFromInt(field.value);
                    const source = pools.index(pool, species, source_layer, axis, root_class);
                    const destination = pools.index(pool, species, destination_layer, axis, root_class);
                    if (!std.math.isFinite(pools.amounts[destination] + fraction * pools.amounts[source]))
                        return error.NonFinitePondRootClassTransferResult;
                }
            }
        }
    }
    for (0..topology.species_count) |species| {
        const eligible = topology.active_root_mass_g_c[topology.rootIndex(species, source_layer, 0)] >
            topology.minimum_root_mass_g_c[species] and
            topology.active_root_mass_g_c[topology.rootIndex(species, destination_layer, 0)] >
                topology.minimum_root_mass_g_c[species];
        if (!eligible) continue;
        for (0..topology.active_axis_count[species]) |axis| {
            for (0..topology.active_root_class_count[species]) |root_class| {
                inline for (std.meta.fields(RootPool)) |field| {
                    const pool: RootPool = @enumFromInt(field.value);
                    const source = pools.index(pool, species, source_layer, axis, root_class);
                    const destination = pools.index(pool, species, destination_layer, axis, root_class);
                    pools.amounts[destination] += fraction * pools.amounts[source];
                }
            }
        }
    }
}

test "REDIST root class pools follow runtime species axis and class counts" {
    const root_len = 1 * 3 * 2;
    const class_len = root_len * 2;
    var root_mass: [root_len]f64 = @splat(0);
    var amounts: [root_pool_count * class_len]f64 = @splat(0);
    const topology = Topology{
        .layer_count = 3,
        .species_count = 1,
        .maximum_axis_count = 2,
        .maximum_root_class_count = 2,
        .active_axis_count = &.{2},
        .active_root_class_count = &.{2},
        .minimum_root_mass_g_c = &.{0.1},
        .active_root_mass_g_c = &root_mass,
    };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    root_mass[topology.rootIndex(0, 1, 0)] = 1;
    root_mass[topology.rootIndex(0, 2, 0)] = 1;
    inline for (std.meta.fields(RootPool)) |field| {
        const pool: RootPool = @enumFromInt(field.value);
        for (0..2) |axis| for (0..2) |root_class| {
            amounts[pools.index(pool, 0, 1, axis, root_class)] = 2;
            amounts[pools.index(pool, 0, 2, axis, root_class)] = 1;
        };
    }
    try transfer(1, 2, 0.25, pools);
    inline for (std.meta.fields(RootPool)) |field| {
        const pool: RootPool = @enumFromInt(field.value);
        for (0..2) |axis| for (0..2) |root_class|
            try std.testing.expectEqual(@as(f64, 1.5), amounts[pools.index(pool, 0, 2, axis, root_class)]);
    }
}

test "REDIST root class transfer obeys strict active-root eligibility" {
    var root_mass: [3]f64 = .{ 0, 0.1, 1 };
    var amounts: [root_pool_count * 3]f64 = @splat(1);
    const topology = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .maximum_root_class_count = 1, .active_axis_count = &.{1}, .active_root_class_count = &.{1}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try transfer(1, 2, 0.5, .{ .topology = topology, .amounts = &amounts });
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST root class runtime counts are validated" {
    var root_mass: [2]f64 = @splat(1);
    var amounts: [root_pool_count * 2]f64 = @splat(1);
    const topology = Topology{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 1, .maximum_root_class_count = 1, .active_axis_count = &.{1}, .active_root_class_count = &.{2}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try std.testing.expectError(error.PondRootClassTransferDimensionMismatch, transfer(0, 1, 0.5, .{ .topology = topology, .amounts = &amounts }));
}

test "REDIST root class zero loop counts are valid" {
    var root_mass: [2]f64 = @splat(1);
    var amounts: [root_pool_count * 2]f64 = @splat(1);
    const topology = Topology{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 1, .maximum_root_class_count = 1, .active_axis_count = &.{0}, .active_root_class_count = &.{0}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try transfer(1, 0, 0.5, .{ .topology = topology, .amounts = &amounts });
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST root class validation is atomic" {
    var root_mass: [3]f64 = .{ 0, 1, 1 };
    var amounts: [root_pool_count * 3]f64 = @splat(1);
    const topology = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .maximum_root_class_count = 1, .active_axis_count = &.{1}, .active_root_class_count = &.{1}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    amounts[pools.index(.secondary_root_axis_count, 0, 1, 0, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondRootClassTransferInput, transfer(1, 2, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.primary_root_carbon, 0, 2, 0, 0)]);
}
