const std = @import("std");
const destination_transfer = @import("root_class_transfer.zig");

pub const RootPool = destination_transfer.RootPool;
pub const Topology = destination_transfer.Topology;
pub const Pools = destination_transfer.Pools;
pub const root_pool_count = destination_transfer.root_pool_count;
pub const unit = destination_transfer.unit;

fn rootIndex(t: Topology, species: usize, layer: usize, axis: usize) usize {
    return (species * t.layer_count + layer) * t.maximum_axis_count + axis;
}
fn classIndex(t: Topology, species: usize, layer: usize, axis: usize, root_class: usize) usize {
    return rootIndex(t, species, layer, axis) * t.maximum_root_class_count + root_class;
}
fn poolIndex(pools: Pools, pool: RootPool, species: usize, layer: usize, axis: usize, root_class: usize) usize {
    const class_len = pools.topology.species_count * pools.topology.layer_count *
        pools.topology.maximum_axis_count * pools.topology.maximum_root_class_count;
    return @intFromEnum(pool) * class_len + classIndex(pools.topology, species, layer, axis, root_class);
}
fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}
fn eligible(t: Topology, species: usize, source: usize, destination: usize) bool {
    return t.active_root_mass_g_c[rootIndex(t, species, source, 0)] > t.minimum_root_mass_g_c[species] and
        t.active_root_mass_g_c[rootIndex(t, species, destination, 0)] > t.minimum_root_mass_g_c[species];
}

/// Direct translation of REDIST 9253--9263 under the open plant guards.
pub fn reduce(source_layer: usize, destination_layer: usize, remaining_fraction: f64, pools: Pools) !void {
    const t = pools.topology;
    const root_len = t.species_count * t.layer_count * t.maximum_axis_count;
    const class_len = root_len * t.maximum_root_class_count;
    if (t.layer_count == 0 or (t.species_count > 0 and t.maximum_axis_count == 0) or
        source_layer >= t.layer_count or destination_layer >= t.layer_count or source_layer == destination_layer or
        t.active_axis_count.len != t.species_count or t.active_root_class_count.len != t.species_count or
        t.minimum_root_mass_g_c.len != t.species_count or t.active_root_mass_g_c.len != root_len or
        pools.amounts.len != root_pool_count * class_len)
        return error.PondSourceRootClassReductionDimensionMismatch;
    for (t.active_axis_count) |count| if (count > t.maximum_axis_count)
        return error.PondSourceRootClassReductionDimensionMismatch;
    for (t.active_root_class_count) |count| if (count > t.maximum_root_class_count)
        return error.PondSourceRootClassReductionDimensionMismatch;
    if (!finiteSlice(t.minimum_root_mass_g_c) or !finiteSlice(t.active_root_mass_g_c) or
        !finiteSlice(pools.amounts) or !std.math.isFinite(remaining_fraction) or
        remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceRootClassReductionInput;
    if (source_layer == 0) return;

    for (0..t.species_count) |species| {
        if (!eligible(t, species, source_layer, destination_layer)) continue;
        for (0..t.active_axis_count[species]) |axis| for (0..t.active_root_class_count[species]) |root_class| inline for (std.meta.fields(RootPool)) |field| {
            const pool: RootPool = @enumFromInt(field.value);
            const index = poolIndex(pools, pool, species, source_layer, axis, root_class);
            if (!std.math.isFinite(remaining_fraction * pools.amounts[index]))
                return error.NonFinitePondSourceRootClassReductionResult;
        };
    }
    for (0..t.species_count) |species| {
        if (!eligible(t, species, source_layer, destination_layer)) continue;
        for (0..t.active_axis_count[species]) |axis| for (0..t.active_root_class_count[species]) |root_class| inline for (std.meta.fields(RootPool)) |field| {
            const pool: RootPool = @enumFromInt(field.value);
            const index = poolIndex(pools, pool, species, source_layer, axis, root_class);
            pools.amounts[index] = remaining_fraction * pools.amounts[index];
        };
    }
}

test "REDIST source root classes scale runtime axes classes and exact pools" {
    const root_len = 1 * 3 * 2;
    const class_len = root_len * 2;
    var root_mass: [root_len]f64 = @splat(0);
    var amounts: [root_pool_count * class_len]f64 = @splat(0);
    const topology = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 2, .maximum_root_class_count = 2, .active_axis_count = &.{2}, .active_root_class_count = &.{2}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    root_mass[rootIndex(topology, 0, 1, 0)] = 1;
    root_mass[rootIndex(topology, 0, 2, 0)] = 1;
    inline for (std.meta.fields(RootPool)) |field| {
        const pool: RootPool = @enumFromInt(field.value);
        for (0..2) |axis| for (0..2) |root_class| {
            amounts[poolIndex(pools, pool, 0, 1, axis, root_class)] = 4;
        };
    }
    try reduce(1, 2, 0.25, pools);
    inline for (std.meta.fields(RootPool)) |field| {
        const pool: RootPool = @enumFromInt(field.value);
        for (0..2) |axis| for (0..2) |root_class| {
            try std.testing.expectEqual(@as(f64, 1), amounts[poolIndex(pools, pool, 0, 1, axis, root_class)]);
        };
    }
}

test "REDIST source root classes exclude layer zero" {
    var root_mass: [2]f64 = @splat(1);
    var amounts: [root_pool_count * 2]f64 = @splat(1);
    const topology = Topology{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 1, .maximum_root_class_count = 1, .active_axis_count = &.{1}, .active_root_class_count = &.{1}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try reduce(0, 1, 0, .{ .topology = topology, .amounts = &amounts });
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST source root class zero runtime counts are zero-trip loops" {
    var root_mass: [3]f64 = .{ 0, 1, 1 };
    var amounts: [root_pool_count * 3]f64 = @splat(1);
    const topology = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .maximum_root_class_count = 1, .active_axis_count = &.{0}, .active_root_class_count = &.{0}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try reduce(1, 2, 0, .{ .topology = topology, .amounts = &amounts });
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST source root class zero species and root-class extent are valid" {
    var no_root_mass: [0]f64 = .{};
    var no_amounts: [0]f64 = .{};
    const no_species = Topology{ .layer_count = 3, .species_count = 0, .maximum_axis_count = 0, .maximum_root_class_count = 0, .active_axis_count = &.{}, .active_root_class_count = &.{}, .minimum_root_mass_g_c = &.{}, .active_root_mass_g_c = &no_root_mass };
    try reduce(1, 2, 0, .{ .topology = no_species, .amounts = &no_amounts });

    var root_mass: [3]f64 = .{ 0, 1, 1 };
    const no_classes = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .maximum_root_class_count = 0, .active_axis_count = &.{1}, .active_root_class_count = &.{0}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try reduce(1, 2, 0, .{ .topology = no_classes, .amounts = &no_amounts });
}

test "REDIST source root class validation is atomic" {
    var root_mass: [3]f64 = .{ 0, 1, 1 };
    var amounts: [root_pool_count * 3]f64 = @splat(1);
    const topology = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .maximum_root_class_count = 1, .active_axis_count = &.{1}, .active_root_class_count = &.{1}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    amounts[poolIndex(pools, .secondary_root_axis_count, 0, 1, 0, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceRootClassReductionInput, reduce(1, 2, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[poolIndex(pools, .primary_root_carbon, 0, 1, 0, 0)]);
}
