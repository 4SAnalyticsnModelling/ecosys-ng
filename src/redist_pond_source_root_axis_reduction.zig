const std = @import("std");
const destination_transfer = @import("redist_pond_root_axis_transfer.zig");

pub const RootAxisPool = destination_transfer.RootAxisPool;
pub const Topology = destination_transfer.Topology;
pub const Pools = destination_transfer.Pools;
pub const root_axis_pool_count = destination_transfer.root_axis_pool_count;
pub const unit = destination_transfer.unit;

fn topologyIndex(t: Topology, species: usize, layer: usize, axis: usize) usize {
    return (species * t.layer_count + layer) * t.maximum_axis_count + axis;
}
fn poolIndex(pools: Pools, pool: RootAxisPool, species: usize, layer: usize, axis: usize) usize {
    const topology_len = pools.topology.species_count * pools.topology.layer_count * pools.topology.maximum_axis_count;
    return @intFromEnum(pool) * topology_len + topologyIndex(pools.topology, species, layer, axis);
}
fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}
fn eligible(t: Topology, species: usize, source: usize, destination: usize) bool {
    return t.active_root_mass_g_c[topologyIndex(t, species, source, 0)] > t.minimum_root_mass_g_c[species] and
        t.active_root_mass_g_c[topologyIndex(t, species, destination, 0)] > t.minimum_root_mass_g_c[species];
}

/// Direct translation of REDIST 9264--9280 under the open plant guards.
pub fn reduce(source_layer: usize, destination_layer: usize, remaining_fraction: f64, pools: Pools) !void {
    const t = pools.topology;
    const topology_len = t.species_count * t.layer_count * t.maximum_axis_count;
    if (t.layer_count == 0 or (t.species_count > 0 and t.maximum_axis_count == 0) or
        source_layer >= t.layer_count or destination_layer >= t.layer_count or source_layer == destination_layer or
        t.active_axis_count.len != t.species_count or t.minimum_root_mass_g_c.len != t.species_count or
        t.active_root_mass_g_c.len != topology_len or pools.amounts.len != root_axis_pool_count * topology_len)
        return error.PondSourceRootAxisReductionDimensionMismatch;
    for (t.active_axis_count) |count| if (count > t.maximum_axis_count)
        return error.PondSourceRootAxisReductionDimensionMismatch;
    if (!finiteSlice(t.minimum_root_mass_g_c) or !finiteSlice(t.active_root_mass_g_c) or
        !finiteSlice(pools.amounts) or !std.math.isFinite(remaining_fraction) or
        remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceRootAxisReductionInput;
    if (source_layer == 0) return;

    for (0..t.species_count) |species| {
        if (!eligible(t, species, source_layer, destination_layer)) continue;
        for (0..t.active_axis_count[species]) |axis| inline for (std.meta.fields(RootAxisPool)) |field| {
            const pool: RootAxisPool = @enumFromInt(field.value);
            const index = poolIndex(pools, pool, species, source_layer, axis);
            if (!std.math.isFinite(remaining_fraction * pools.amounts[index]))
                return error.NonFinitePondSourceRootAxisReductionResult;
        };
    }
    for (0..t.species_count) |species| {
        if (!eligible(t, species, source_layer, destination_layer)) continue;
        for (0..t.active_axis_count[species]) |axis| inline for (std.meta.fields(RootAxisPool)) |field| {
            const pool: RootAxisPool = @enumFromInt(field.value);
            const index = poolIndex(pools, pool, species, source_layer, axis);
            pools.amounts[index] = remaining_fraction * pools.amounts[index];
        };
    }
}

test "REDIST source root-axis aggregates scale exact pools and runtime axes" {
    const topology_len = 1 * 3 * 2;
    var root_mass: [topology_len]f64 = @splat(0);
    var amounts: [root_axis_pool_count * topology_len]f64 = @splat(0);
    const topology = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 2, .active_axis_count = &.{2}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    root_mass[topologyIndex(topology, 0, 1, 0)] = 1;
    root_mass[topologyIndex(topology, 0, 2, 0)] = 1;
    inline for (std.meta.fields(RootAxisPool)) |field| {
        const pool: RootAxisPool = @enumFromInt(field.value);
        amounts[poolIndex(pools, pool, 0, 1, 0)] = 4;
        amounts[poolIndex(pools, pool, 0, 1, 1)] = 8;
    }
    try reduce(1, 2, 0.25, pools);
    inline for (std.meta.fields(RootAxisPool)) |field| {
        const pool: RootAxisPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 1), amounts[poolIndex(pools, pool, 0, 1, 0)]);
        try std.testing.expectEqual(@as(f64, 2), amounts[poolIndex(pools, pool, 0, 1, 1)]);
    }
}

test "REDIST source root-axis NP and MY zero are zero-trip" {
    var no_root_mass: [0]f64 = .{};
    var no_amounts: [0]f64 = .{};
    const no_species = Topology{ .layer_count = 3, .species_count = 0, .maximum_axis_count = 0, .active_axis_count = &.{}, .minimum_root_mass_g_c = &.{}, .active_root_mass_g_c = &no_root_mass };
    try reduce(1, 2, 0, .{ .topology = no_species, .amounts = &no_amounts });
    var root_mass: [3]f64 = .{ 0, 1, 1 };
    var amounts: [root_axis_pool_count * 3]f64 = @splat(1);
    const no_axes = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .active_axis_count = &.{0}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try reduce(1, 2, 0, .{ .topology = no_axes, .amounts = &amounts });
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST source root-axis validation is atomic" {
    var root_mass: [3]f64 = .{ 0, 1, 1 };
    var amounts: [root_axis_pool_count * 3]f64 = @splat(1);
    const topology = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .active_axis_count = &.{1}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    amounts[poolIndex(pools, .average_secondary_length, 0, 1, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceRootAxisReductionInput, reduce(1, 2, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[poolIndex(pools, .nonstructural_carbon, 0, 1, 0)]);
}
