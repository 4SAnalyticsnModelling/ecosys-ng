const std = @import("std");

/// Declaration order matches REDIST 8802--8825 for each root axis N.
pub const GasPool = enum(u8) {
    gaseous_carbon_dioxide,
    gaseous_oxygen,
    gaseous_methane,
    gaseous_nitrous_oxide,
    gaseous_ammonia,
    gaseous_hydrogen,
    aqueous_carbon_dioxide,
    aqueous_oxygen,
    aqueous_methane,
    aqueous_nitrous_oxide,
    aqueous_ammonia,
    aqueous_hydrogen,
};

pub const gas_pool_count = std.meta.fields(GasPool).len;

pub const ElementUnit = enum { g_c, g_o, g_n, g_h };

pub fn unit(pool: GasPool) ElementUnit {
    return switch (pool) {
        .gaseous_carbon_dioxide, .gaseous_methane, .aqueous_carbon_dioxide, .aqueous_methane => .g_c,
        .gaseous_oxygen, .aqueous_oxygen => .g_o,
        .gaseous_nitrous_oxide, .gaseous_ammonia, .aqueous_nitrous_oxide, .aqueous_ammonia => .g_n,
        .gaseous_hydrogen, .aqueous_hydrogen => .g_h,
    };
}

pub const Topology = struct {
    layer_count: usize,
    species_count: usize, // NP
    maximum_axis_count: usize,
    active_axis_count: []const usize, // MY(NZ)
    minimum_root_mass_g_c: []const f64, // ZEROP(NZ)
    /// Species-major, then layer, then root axis. WTRTL.
    active_root_mass_g_c: []const f64,

    fn rootIndex(self: Topology, species: usize, layer: usize, axis: usize) usize {
        return (species * self.layer_count + layer) * self.maximum_axis_count + axis;
    }
};

pub const Pools = struct {
    topology: Topology,
    /// Pool-major over the topology's root storage.
    amounts: []f64,

    fn index(self: Pools, pool: GasPool, species: usize, layer: usize, axis: usize) usize {
        const root_storage_len = self.topology.species_count * self.topology.layer_count *
            self.topology.maximum_axis_count;
        return @intFromEnum(pool) * root_storage_len + self.topology.rootIndex(species, layer, axis);
    }
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8794--8825 through the pond root-gas block.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, pools: Pools) !void {
    const topology = pools.topology;
    const root_storage_len = topology.species_count * topology.layer_count * topology.maximum_axis_count;
    if (topology.layer_count == 0 or (topology.species_count != 0 and topology.maximum_axis_count == 0) or
        source_layer >= topology.layer_count or destination_layer >= topology.layer_count or
        source_layer == destination_layer or topology.active_axis_count.len != topology.species_count or
        topology.minimum_root_mass_g_c.len != topology.species_count or
        topology.active_root_mass_g_c.len != root_storage_len or pools.amounts.len != gas_pool_count * root_storage_len)
        return error.PondRootGasTransferDimensionMismatch;
    for (topology.active_axis_count) |count| if (count > topology.maximum_axis_count)
        return error.PondRootGasTransferDimensionMismatch;
    if (!finiteSlice(topology.minimum_root_mass_g_c) or !finiteSlice(topology.active_root_mass_g_c) or
        !finiteSlice(pools.amounts) or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondRootGasTransferInput;
    if (source_layer == 0) return;

    for (0..topology.species_count) |species| {
        const source_root = topology.active_root_mass_g_c[topology.rootIndex(species, source_layer, 0)];
        const destination_root = topology.active_root_mass_g_c[topology.rootIndex(species, destination_layer, 0)];
        if (source_root > topology.minimum_root_mass_g_c[species] and
            destination_root > topology.minimum_root_mass_g_c[species])
        {
            for (0..topology.active_axis_count[species]) |axis| {
                inline for (std.meta.fields(GasPool)) |field| {
                    const pool: GasPool = @enumFromInt(field.value);
                    const source = pools.index(pool, species, source_layer, axis);
                    const destination = pools.index(pool, species, destination_layer, axis);
                    if (!std.math.isFinite(pools.amounts[destination] + fraction * pools.amounts[source]))
                        return error.NonFinitePondRootGasTransferResult;
                }
            }
        }
    }
    for (0..topology.species_count) |species| {
        const source_root = topology.active_root_mass_g_c[topology.rootIndex(species, source_layer, 0)];
        const destination_root = topology.active_root_mass_g_c[topology.rootIndex(species, destination_layer, 0)];
        if (source_root > topology.minimum_root_mass_g_c[species] and
            destination_root > topology.minimum_root_mass_g_c[species])
        {
            for (0..topology.active_axis_count[species]) |axis| {
                inline for (std.meta.fields(GasPool)) |field| {
                    const pool: GasPool = @enumFromInt(field.value);
                    const source = pools.index(pool, species, source_layer, axis);
                    const destination = pools.index(pool, species, destination_layer, axis);
                    pools.amounts[destination] += fraction * pools.amounts[source];
                }
            }
        }
    }
}

test "REDIST pond root gases use runtime species and axis counts" {
    const topology_shape = 2 * 3 * 2;
    var root_mass: [topology_shape]f64 = @splat(0);
    var amounts: [gas_pool_count * topology_shape]f64 = @splat(0);
    const topology = Topology{
        .layer_count = 3,
        .species_count = 2,
        .maximum_axis_count = 2,
        .active_axis_count = &.{ 2, 0 },
        .minimum_root_mass_g_c = &.{ 0.1, 0.1 },
        .active_root_mass_g_c = &root_mass,
    };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    root_mass[topology.rootIndex(0, 1, 0)] = 1;
    root_mass[topology.rootIndex(0, 2, 0)] = 1;
    root_mass[topology.rootIndex(1, 1, 0)] = 0.05;
    root_mass[topology.rootIndex(1, 2, 0)] = 1;
    inline for (std.meta.fields(GasPool)) |field| {
        const pool: GasPool = @enumFromInt(field.value);
        amounts[pools.index(pool, 0, 1, 0)] = 2;
        amounts[pools.index(pool, 0, 1, 1)] = 4;
        amounts[pools.index(pool, 0, 2, 0)] = 1;
        amounts[pools.index(pool, 0, 2, 1)] = 1;
        amounts[pools.index(pool, 1, 2, 0)] = 3;
    }
    try transfer(1, 2, 0.25, pools);
    inline for (std.meta.fields(GasPool)) |field| {
        const pool: GasPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 1.5), amounts[pools.index(pool, 0, 2, 0)]);
        try std.testing.expectEqual(@as(f64, 2), amounts[pools.index(pool, 0, 2, 1)]);
        try std.testing.expectEqual(@as(f64, 3), amounts[pools.index(pool, 1, 2, 0)]);
    }
}

test "REDIST pond root gas source layer zero is excluded" {
    var root_mass: [4]f64 = @splat(1);
    var amounts: [gas_pool_count * 4]f64 = @splat(1);
    const topology = Topology{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 2, .active_axis_count = &.{1}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try transfer(0, 1, 0.5, .{ .topology = topology, .amounts = &amounts });
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST pond root gas runtime topology is validated" {
    var root_mass: [4]f64 = @splat(1);
    var amounts: [gas_pool_count * 4]f64 = @splat(1);
    const topology = Topology{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 2, .active_axis_count = &.{3}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    try std.testing.expectError(error.PondRootGasTransferDimensionMismatch, transfer(0, 1, 0.5, .{ .topology = topology, .amounts = &amounts }));
}

test "REDIST pond root gas validation is atomic" {
    var root_mass: [2]f64 = @splat(1);
    var amounts: [gas_pool_count * 2]f64 = @splat(1);
    const topology = Topology{ .layer_count = 2, .species_count = 1, .maximum_axis_count = 1, .active_axis_count = &.{1}, .minimum_root_mass_g_c = &.{0.1}, .active_root_mass_g_c = &root_mass };
    const pools = Pools{ .topology = topology, .amounts = &amounts };
    amounts[pools.index(.aqueous_hydrogen, 0, 0, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondRootGasTransferInput, transfer(0, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.gaseous_carbon_dioxide, 0, 1, 0)]);
}
