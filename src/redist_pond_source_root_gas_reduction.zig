const std = @import("std");
const destination_transfer = @import("redist_pond_root_gas_transfer.zig");

pub const GasPool = destination_transfer.GasPool;
pub const Topology = destination_transfer.Topology;
pub const Pools = destination_transfer.Pools;
pub const gas_pool_count = destination_transfer.gas_pool_count;
pub const unit = destination_transfer.unit;

fn rootIndex(t: Topology, species: usize, layer: usize, axis: usize) usize {
    return (species * t.layer_count + layer) * t.maximum_axis_count + axis;
}

fn poolIndex(pools: Pools, pool: GasPool, species: usize, layer: usize, axis: usize) usize {
    const topology_len = pools.topology.species_count * pools.topology.layer_count * pools.topology.maximum_axis_count;
    return @intFromEnum(pool) * topology_len + rootIndex(pools.topology, species, layer, axis);
}

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 9230--9249 through source pond root gases.
pub fn reduce(source_layer: usize, destination_layer: usize, remaining_fraction: f64, pools: Pools) !void {
    const t = pools.topology;
    const topology_len = t.species_count * t.layer_count * t.maximum_axis_count;
    if (t.layer_count == 0 or (t.species_count != 0 and t.maximum_axis_count == 0) or
        source_layer >= t.layer_count or destination_layer >= t.layer_count or source_layer == destination_layer or
        t.active_axis_count.len != t.species_count or t.minimum_root_mass_g_c.len != t.species_count or
        t.active_root_mass_g_c.len != topology_len or pools.amounts.len != gas_pool_count * topology_len)
        return error.PondSourceRootGasReductionDimensionMismatch;
    for (t.active_axis_count) |count| if (count > t.maximum_axis_count)
        return error.PondSourceRootGasReductionDimensionMismatch;
    if (!finiteSlice(t.minimum_root_mass_g_c) or !finiteSlice(t.active_root_mass_g_c) or
        !finiteSlice(pools.amounts) or !std.math.isFinite(remaining_fraction) or
        remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceRootGasReductionInput;
    if (source_layer == 0) return;

    for (0..t.species_count) |species| {
        const eligible = t.active_root_mass_g_c[rootIndex(t, species, source_layer, 0)] > t.minimum_root_mass_g_c[species] and
            t.active_root_mass_g_c[rootIndex(t, species, destination_layer, 0)] > t.minimum_root_mass_g_c[species];
        if (!eligible) continue;
        for (0..t.active_axis_count[species]) |axis| inline for (std.meta.fields(GasPool)) |field| {
            const pool: GasPool = @enumFromInt(field.value);
            const index = poolIndex(pools, pool, species, source_layer, axis);
            if (!std.math.isFinite(remaining_fraction * pools.amounts[index]))
                return error.NonFinitePondSourceRootGasReductionResult;
        };
    }
    for (0..t.species_count) |species| {
        const eligible = t.active_root_mass_g_c[rootIndex(t, species, source_layer, 0)] > t.minimum_root_mass_g_c[species] and
            t.active_root_mass_g_c[rootIndex(t, species, destination_layer, 0)] > t.minimum_root_mass_g_c[species];
        if (!eligible) continue;
        for (0..t.active_axis_count[species]) |axis| inline for (std.meta.fields(GasPool)) |field| {
            const pool: GasPool = @enumFromInt(field.value);
            const index = poolIndex(pools, pool, species, source_layer, axis);
            pools.amounts[index] = remaining_fraction * pools.amounts[index];
        };
    }
}

test "REDIST source root gas accepts zero-trip NP and MY loops" {
    var empty_amounts: [0]f64 = .{};
    try reduce(0, 1, 0.5, .{ .topology = .{ .layer_count = 2, .species_count = 0, .maximum_axis_count = 0, .active_axis_count = &.{}, .minimum_root_mass_g_c = &.{}, .active_root_mass_g_c = &.{} }, .amounts = &empty_amounts });

    var root_mass: [2]f64 = @splat(1);
    var amounts: [gas_pool_count * 2]f64 = @splat(1);
    try reduce(1, 0, 0.5, .{ .topology = .{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 1, .active_axis_count = &.{0}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass }, .amounts = &amounts });
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST source root gases scale runtime eligible species and axes" {
    const topology_len = 2 * 3 * 2;
    var root_mass: [topology_len]f64 = @splat(0);
    var amounts: [gas_pool_count * topology_len]f64 = @splat(0);
    const topology = Topology{ .layer_count = 3, .species_count = 2, .maximum_axis_count = 2, .active_axis_count = &.{ 2, 1 }, .minimum_root_mass_g_c = &.{ 0.1, 0.1 }, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    root_mass[rootIndex(topology, 0, 1, 0)] = 1;
    root_mass[rootIndex(topology, 0, 2, 0)] = 1;
    root_mass[rootIndex(topology, 1, 1, 0)] = 0.1;
    root_mass[rootIndex(topology, 1, 2, 0)] = 1;
    inline for (std.meta.fields(GasPool)) |field| {
        const pool: GasPool = @enumFromInt(field.value);
        amounts[poolIndex(pools, pool, 0, 1, 0)] = 4;
        amounts[poolIndex(pools, pool, 0, 1, 1)] = 8;
        amounts[poolIndex(pools, pool, 1, 1, 0)] = 12;
    }
    try reduce(1, 2, 0.25, pools);
    inline for (std.meta.fields(GasPool)) |field| {
        const pool: GasPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 1), amounts[poolIndex(pools, pool, 0, 1, 0)]);
        try std.testing.expectEqual(@as(f64, 2), amounts[poolIndex(pools, pool, 0, 1, 1)]);
        try std.testing.expectEqual(@as(f64, 12), amounts[poolIndex(pools, pool, 1, 1, 0)]);
    }
}

test "REDIST source root gas layer zero is excluded" {
    var root_mass: [2]f64 = @splat(1);
    var amounts: [gas_pool_count * 2]f64 = @splat(1);
    const topology = Topology{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 1, .active_axis_count = &.{1}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try reduce(0, 1, 0, .{ .topology = topology, .amounts = &amounts });
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST source root gas topology validation is atomic" {
    var root_mass: [3]f64 = .{ 0, 1, 1 };
    var amounts: [gas_pool_count * 3]f64 = @splat(1);
    const topology = Topology{ .layer_count = 3, .species_count = 1, .maximum_axis_count = 1, .active_axis_count = &.{1}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    amounts[poolIndex(pools, .aqueous_hydrogen, 0, 1, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceRootGasReductionInput, reduce(1, 2, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[poolIndex(pools, .gaseous_carbon_dioxide, 0, 1, 0)]);
}
