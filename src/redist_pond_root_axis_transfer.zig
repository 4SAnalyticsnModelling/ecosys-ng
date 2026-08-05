const std = @import("std");

/// Exact REDIST 8849--8880 assignment order for each root axis N.
pub const RootAxisPool = enum(u8) {
    nonstructural_carbon,
    nonstructural_nitrogen,
    nonstructural_phosphorus,
    active_root_carbon,
    actual_root_carbon,
    protein_carbon,
    primary_axis_count,
    secondary_axis_count,
    root_length_per_plant,
    root_length_density,
    gaseous_volume,
    aqueous_volume,
    primary_radius,
    secondary_radius,
    surface_area_per_plant,
    average_secondary_length,
};

pub const root_axis_pool_count = std.meta.fields(RootAxisPool).len;
pub const Unit = enum { g_c, g_n, g_p, count, m, m_per_m3, m2, m3 };

pub fn unit(pool: RootAxisPool) Unit {
    return switch (pool) {
        .nonstructural_carbon, .active_root_carbon, .actual_root_carbon, .protein_carbon => .g_c,
        .nonstructural_nitrogen => .g_n,
        .nonstructural_phosphorus => .g_p,
        .primary_axis_count, .secondary_axis_count => .count,
        .root_length_per_plant, .primary_radius, .secondary_radius, .average_secondary_length => .m,
        .root_length_density => .m_per_m3,
        .surface_area_per_plant => .m2,
        .gaseous_volume, .aqueous_volume => .m3,
    };
}

pub const Topology = struct {
    layer_count: usize,
    species_count: usize,
    maximum_axis_count: usize,
    active_axis_count: []const usize, // MY(NZ)
    minimum_root_mass_g_c: []const f64, // ZEROP(NZ)
    active_root_mass_g_c: []const f64, // WTRTL

    fn index(self: Topology, species: usize, layer: usize, axis: usize) usize {
        return (species * self.layer_count + layer) * self.maximum_axis_count + axis;
    }
};

pub const Pools = struct {
    topology: Topology,
    amounts: []f64,

    fn index(self: Pools, pool: RootAxisPool, species: usize, layer: usize, axis: usize) usize {
        const topology_len = self.topology.species_count * self.topology.layer_count * self.topology.maximum_axis_count;
        return @intFromEnum(pool) * topology_len + self.topology.index(species, layer, axis);
    }
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

fn eligible(topology: Topology, species: usize, source_layer: usize, destination_layer: usize) bool {
    return topology.active_root_mass_g_c[topology.index(species, source_layer, 0)] >
        topology.minimum_root_mass_g_c[species] and
        topology.active_root_mass_g_c[topology.index(species, destination_layer, 0)] >
            topology.minimum_root_mass_g_c[species];
}

/// Direct translation of REDIST 8849--8880 under the existing plant guards.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, pools: Pools) !void {
    const topology = pools.topology;
    const topology_len = topology.species_count * topology.layer_count * topology.maximum_axis_count;
    if (topology.layer_count == 0 or (topology.species_count > 0 and topology.maximum_axis_count == 0) or
        source_layer >= topology.layer_count or destination_layer >= topology.layer_count or
        source_layer == destination_layer or topology.active_axis_count.len != topology.species_count or
        topology.minimum_root_mass_g_c.len != topology.species_count or
        topology.active_root_mass_g_c.len != topology_len or pools.amounts.len != root_axis_pool_count * topology_len)
        return error.PondRootAxisTransferDimensionMismatch;
    for (topology.active_axis_count) |count| if (count > topology.maximum_axis_count)
        return error.PondRootAxisTransferDimensionMismatch;
    if (!finiteSlice(topology.minimum_root_mass_g_c) or !finiteSlice(topology.active_root_mass_g_c) or
        !finiteSlice(pools.amounts) or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondRootAxisTransferInput;
    if (source_layer == 0) return;

    for (0..topology.species_count) |species| {
        if (!eligible(topology, species, source_layer, destination_layer)) continue;
        for (0..topology.active_axis_count[species]) |axis| {
            inline for (std.meta.fields(RootAxisPool)) |field| {
                const pool: RootAxisPool = @enumFromInt(field.value);
                const source = pools.index(pool, species, source_layer, axis);
                const destination = pools.index(pool, species, destination_layer, axis);
                if (!std.math.isFinite(pools.amounts[destination] + fraction * pools.amounts[source]))
                    return error.NonFinitePondRootAxisTransferResult;
            }
        }
    }
    for (0..topology.species_count) |species| {
        if (!eligible(topology, species, source_layer, destination_layer)) continue;
        for (0..topology.active_axis_count[species]) |axis| {
            inline for (std.meta.fields(RootAxisPool)) |field| {
                const pool: RootAxisPool = @enumFromInt(field.value);
                const source = pools.index(pool, species, source_layer, axis);
                const destination = pools.index(pool, species, destination_layer, axis);
                pools.amounts[destination] += fraction * pools.amounts[source];
            }
        }
    }
}

test "REDIST per-root-axis pools follow runtime species and axis topology" {
    const topology_len = 2 * 3 * 2;
    var root_mass: [topology_len]f64 = @splat(0);
    var amounts: [root_axis_pool_count * topology_len]f64 = @splat(0);
    const topology = Topology{ .layer_count = 3, .species_count = 2, .maximum_axis_count = 2, .active_axis_count = &.{ 2, 1 }, .minimum_root_mass_g_c = &.{ 0.1, 0.1 }, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    root_mass[topology.index(0, 1, 0)] = 1;
    root_mass[topology.index(0, 2, 0)] = 1;
    root_mass[topology.index(1, 1, 0)] = 0.1;
    root_mass[topology.index(1, 2, 0)] = 1;
    inline for (std.meta.fields(RootAxisPool)) |field| {
        const pool: RootAxisPool = @enumFromInt(field.value);
        for (0..2) |axis| {
            amounts[pools.index(pool, 0, 1, axis)] = 2;
            amounts[pools.index(pool, 0, 2, axis)] = 1;
        }
        amounts[pools.index(pool, 1, 2, 0)] = 3;
    }
    try transfer(1, 2, 0.25, pools);
    inline for (std.meta.fields(RootAxisPool)) |field| {
        const pool: RootAxisPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 1.5), amounts[pools.index(pool, 0, 2, 0)]);
        try std.testing.expectEqual(@as(f64, 1.5), amounts[pools.index(pool, 0, 2, 1)]);
        try std.testing.expectEqual(@as(f64, 3), amounts[pools.index(pool, 1, 2, 0)]);
    }
}

test "REDIST root axis pool unit map covers geometry and nutrients" {
    try std.testing.expectEqual(Unit.g_c, unit(.active_root_carbon));
    try std.testing.expectEqual(Unit.g_n, unit(.nonstructural_nitrogen));
    try std.testing.expectEqual(Unit.m_per_m3, unit(.root_length_density));
    try std.testing.expectEqual(Unit.m3, unit(.aqueous_volume));
    try std.testing.expectEqual(Unit.m2, unit(.surface_area_per_plant));
}

test "REDIST root axis runtime topology is validated" {
    var root_mass: [2]f64 = @splat(1);
    var amounts: [root_axis_pool_count * 2]f64 = @splat(1);
    const topology = Topology{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 1, .active_axis_count = &.{2}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try std.testing.expectError(error.PondRootAxisTransferDimensionMismatch, transfer(0, 1, 0.5, .{ .topology = topology, .amounts = &amounts }));
}

test "REDIST root axis zero MY is a valid empty loop" {
    var root_mass: [2]f64 = @splat(1);
    var amounts: [root_axis_pool_count * 2]f64 = @splat(1);
    const topology = Topology{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 1, .active_axis_count = &.{0}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try transfer(1, 0, 0.5, .{ .topology = topology, .amounts = &amounts });
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST root axis validation is atomic" {
    var root_mass: [3]f64 = .{ 0, 1, 1 };
    var amounts: [root_axis_pool_count * 3]f64 = @splat(1);
    const topology = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .active_axis_count = &.{1}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    amounts[pools.index(.average_secondary_length, 0, 1, 0)] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidPondRootAxisTransferInput, transfer(1, 2, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.nonstructural_carbon, 0, 2, 0)]);
}
