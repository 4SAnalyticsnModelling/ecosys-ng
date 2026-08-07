const std = @import("std");

pub const Topology = struct {
    layer_count: usize,
    organic_fraction_count: usize = 6, // K=0..5 for OM*
    organic_group_count: usize = 7, // N=1..7
    organic_component_count: usize = 3, // M=1..3
    residue_fraction_count: usize = 5, // K=0..4
    residue_component_count: usize = 2, // M=1..2 for OR*
    soluble_component_count: usize = 5, // M=1..5 for OS*

    fn organicIndex(self: Topology, fraction: usize, group: usize, component: usize, layer: usize) usize {
        return (((fraction * self.organic_group_count + group) * self.organic_component_count + component) * self.layer_count) + layer;
    }
    fn residueIndex(self: Topology, fraction: usize, component: usize, layer: usize) usize {
        return ((fraction * self.residue_component_count + component) * self.layer_count) + layer;
    }
    fn fractionIndex(self: Topology, fraction: usize, layer: usize) usize {
        return fraction * self.layer_count + layer;
    }
    fn solubleIndex(self: Topology, fraction: usize, component: usize, layer: usize) usize {
        return ((fraction * self.soluble_component_count + component) * self.layer_count) + layer;
    }
};

pub const Pools = struct {
    topology: Topology,
    total_organic_carbon_g_c: []f64, // ORGC
    organic_carbon_g_c: []f64, // OMC
    organic_nitrogen_g_n: []f64, // OMN
    organic_phosphorus_g_p: []f64, // OMP
    residue_carbon_g_c: []f64, // ORC
    residue_nitrogen_g_n: []f64, // ORN
    residue_phosphorus_g_p: []f64, // ORP
    humus_carbon_g_c: []f64, // OHC
    humus_nitrogen_g_n: []f64, // OHN
    humus_phosphorus_g_p: []f64, // OHP
    humus_acetate_g_c: []f64, // OHA
    soluble_carbon_g_c: []f64, // OSC
    soluble_acetate_g_c: []f64, // OSA
    soluble_nitrogen_g_n: []f64, // OSN
    soluble_phosphorus_g_p: []f64, // OSP
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

fn validateMove(destination: f64, fraction: f64, source: f64) !void {
    if (!std.math.isFinite(destination + fraction * source))
        return error.NonFinitePondSolidOrganicTransferResult;
}

fn move(values: []f64, destination: usize, fraction: f64, source: usize) void {
    values[destination] = values[destination] + fraction * values[source];
}

fn validateDimensions(pools: Pools) !void {
    const t = pools.topology;
    if (t.layer_count == 0 or t.organic_fraction_count != 6 or t.organic_group_count != 7 or
        t.organic_component_count != 3 or t.residue_fraction_count != 5 or
        t.residue_component_count != 2 or t.soluble_component_count != 5)
        return error.PondSolidOrganicTransferDimensionMismatch;
    const organic_len = t.organic_fraction_count * t.organic_group_count * t.organic_component_count * t.layer_count;
    const residue_len = t.residue_fraction_count * t.residue_component_count * t.layer_count;
    const fraction_len = t.residue_fraction_count * t.layer_count;
    const soluble_len = t.residue_fraction_count * t.soluble_component_count * t.layer_count;
    if (pools.total_organic_carbon_g_c.len != t.layer_count or pools.organic_carbon_g_c.len != organic_len or
        pools.organic_nitrogen_g_n.len != organic_len or pools.organic_phosphorus_g_p.len != organic_len or
        pools.residue_carbon_g_c.len != residue_len or pools.residue_nitrogen_g_n.len != residue_len or
        pools.residue_phosphorus_g_p.len != residue_len or pools.humus_carbon_g_c.len != fraction_len or
        pools.humus_nitrogen_g_n.len != fraction_len or pools.humus_phosphorus_g_p.len != fraction_len or
        pools.humus_acetate_g_c.len != fraction_len or pools.soluble_carbon_g_c.len != soluble_len or
        pools.soluble_acetate_g_c.len != soluble_len or pools.soluble_nitrogen_g_n.len != soluble_len or
        pools.soluble_phosphorus_g_p.len != soluble_len)
        return error.PondSolidOrganicTransferDimensionMismatch;
}

/// Direct translation of REDIST 9039--9080, including both enclosing gates.
pub fn transfer(
    surface_reappearance_flag: u8,
    source_layer: usize,
    destination_layer: usize,
    fraction: f64,
    pools: Pools,
) !void {
    try validateDimensions(pools);
    const t = pools.topology;
    if (source_layer >= t.layer_count or destination_layer >= t.layer_count or source_layer == destination_layer)
        return error.PondSolidOrganicTransferDimensionMismatch;
    inline for (.{ pools.total_organic_carbon_g_c, pools.organic_carbon_g_c, pools.organic_nitrogen_g_n, pools.organic_phosphorus_g_p, pools.residue_carbon_g_c, pools.residue_nitrogen_g_n, pools.residue_phosphorus_g_p, pools.humus_carbon_g_c, pools.humus_nitrogen_g_n, pools.humus_phosphorus_g_p, pools.humus_acetate_g_c, pools.soluble_carbon_g_c, pools.soluble_acetate_g_c, pools.soluble_nitrogen_g_n, pools.soluble_phosphorus_g_p }) |values|
        if (!finiteSlice(values)) return error.InvalidPondSolidOrganicTransferInput;
    if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondSolidOrganicTransferInput;
    if (fraction != 1.0 or surface_reappearance_flag != 0) return;

    try validateMove(pools.total_organic_carbon_g_c[destination_layer], fraction, pools.total_organic_carbon_g_c[source_layer]);
    for (0..t.organic_fraction_count) |k| for (0..t.organic_group_count) |n| for (0..t.organic_component_count) |m| {
        const source = t.organicIndex(k, n, m, source_layer);
        const destination = t.organicIndex(k, n, m, destination_layer);
        try validateMove(pools.organic_carbon_g_c[destination], fraction, pools.organic_carbon_g_c[source]);
        try validateMove(pools.organic_nitrogen_g_n[destination], fraction, pools.organic_nitrogen_g_n[source]);
        try validateMove(pools.organic_phosphorus_g_p[destination], fraction, pools.organic_phosphorus_g_p[source]);
    };
    for (0..t.residue_fraction_count) |k| {
        for (0..t.residue_component_count) |m| {
            const source = t.residueIndex(k, m, source_layer);
            const destination = t.residueIndex(k, m, destination_layer);
            try validateMove(pools.residue_carbon_g_c[destination], fraction, pools.residue_carbon_g_c[source]);
            try validateMove(pools.residue_nitrogen_g_n[destination], fraction, pools.residue_nitrogen_g_n[source]);
            try validateMove(pools.residue_phosphorus_g_p[destination], fraction, pools.residue_phosphorus_g_p[source]);
        }
        const source_fraction = t.fractionIndex(k, source_layer);
        const destination_fraction = t.fractionIndex(k, destination_layer);
        inline for (.{ pools.humus_carbon_g_c, pools.humus_nitrogen_g_n, pools.humus_phosphorus_g_p, pools.humus_acetate_g_c }) |values|
            try validateMove(values[destination_fraction], fraction, values[source_fraction]);
        for (0..t.soluble_component_count) |m| {
            const source = t.solubleIndex(k, m, source_layer);
            const destination = t.solubleIndex(k, m, destination_layer);
            inline for (.{ pools.soluble_carbon_g_c, pools.soluble_acetate_g_c, pools.soluble_nitrogen_g_n, pools.soluble_phosphorus_g_p }) |values|
                try validateMove(values[destination], fraction, values[source]);
        }
    }

    move(pools.total_organic_carbon_g_c, destination_layer, fraction, source_layer);
    for (0..t.organic_fraction_count) |k| for (0..t.organic_group_count) |n| for (0..t.organic_component_count) |m| {
        const source = t.organicIndex(k, n, m, source_layer);
        const destination = t.organicIndex(k, n, m, destination_layer);
        move(pools.organic_carbon_g_c, destination, fraction, source);
        move(pools.organic_nitrogen_g_n, destination, fraction, source);
        move(pools.organic_phosphorus_g_p, destination, fraction, source);
    };
    for (0..t.residue_fraction_count) |k| {
        for (0..t.residue_component_count) |m| {
            const source = t.residueIndex(k, m, source_layer);
            const destination = t.residueIndex(k, m, destination_layer);
            move(pools.residue_carbon_g_c, destination, fraction, source);
            move(pools.residue_nitrogen_g_n, destination, fraction, source);
            move(pools.residue_phosphorus_g_p, destination, fraction, source);
        }
        const source_fraction = t.fractionIndex(k, source_layer);
        const destination_fraction = t.fractionIndex(k, destination_layer);
        inline for (.{ pools.humus_carbon_g_c, pools.humus_nitrogen_g_n, pools.humus_phosphorus_g_p, pools.humus_acetate_g_c }) |values|
            move(values, destination_fraction, fraction, source_fraction);
        for (0..t.soluble_component_count) |m| {
            const source = t.solubleIndex(k, m, source_layer);
            const destination = t.solubleIndex(k, m, destination_layer);
            inline for (.{ pools.soluble_carbon_g_c, pools.soluble_acetate_g_c, pools.soluble_nitrogen_g_n, pools.soluble_phosphorus_g_p }) |values|
                move(values, destination, fraction, source);
        }
    }
}

const Fixture = struct {
    total: [2]f64 = .{ 2, 1 },
    organic_c: [252]f64 = @splat(1),
    organic_n: [252]f64 = @splat(1),
    organic_p: [252]f64 = @splat(1),
    residue_c: [20]f64 = @splat(1),
    residue_n: [20]f64 = @splat(1),
    residue_p: [20]f64 = @splat(1),
    humus_c: [10]f64 = @splat(1),
    humus_n: [10]f64 = @splat(1),
    humus_p: [10]f64 = @splat(1),
    humus_a: [10]f64 = @splat(1),
    soluble_c: [50]f64 = @splat(1),
    soluble_a: [50]f64 = @splat(1),
    soluble_n: [50]f64 = @splat(1),
    soluble_p: [50]f64 = @splat(1),
    fn pools(self: *Fixture) Pools {
        return .{ .topology = .{ .layer_count = 2 }, .total_organic_carbon_g_c = &self.total, .organic_carbon_g_c = &self.organic_c, .organic_nitrogen_g_n = &self.organic_n, .organic_phosphorus_g_p = &self.organic_p, .residue_carbon_g_c = &self.residue_c, .residue_nitrogen_g_n = &self.residue_n, .residue_phosphorus_g_p = &self.residue_p, .humus_carbon_g_c = &self.humus_c, .humus_nitrogen_g_n = &self.humus_n, .humus_phosphorus_g_p = &self.humus_p, .humus_acetate_g_c = &self.humus_a, .soluble_carbon_g_c = &self.soluble_c, .soluble_acetate_g_c = &self.soluble_a, .soluble_nitrogen_g_n = &self.soluble_n, .soluble_phosphorus_g_p = &self.soluble_p };
    }
};

test "REDIST full-fraction solid organic pools transfer through all fixed loops" {
    var fixture = Fixture{};
    try transfer(0, 0, 1, 1, fixture.pools());
    try std.testing.expectEqual(@as(f64, 3), fixture.total[1]);
    try std.testing.expectEqual(@as(f64, 2), fixture.organic_c[1]);
    try std.testing.expectEqual(@as(f64, 2), fixture.organic_p[251]);
    try std.testing.expectEqual(@as(f64, 2), fixture.residue_n[19]);
    try std.testing.expectEqual(@as(f64, 2), fixture.humus_a[9]);
    try std.testing.expectEqual(@as(f64, 2), fixture.soluble_p[49]);
}

test "REDIST solid organic transfer obeys both gates" {
    var fixture = Fixture{};
    try transfer(1, 0, 1, 1, fixture.pools());
    try std.testing.expectEqual(@as(f64, 1), fixture.total[1]);
    try transfer(0, 0, 1, 0.999999999, fixture.pools());
    try std.testing.expectEqual(@as(f64, 1), fixture.total[1]);
}

test "REDIST solid organic zero category count is invalid not zero-trip" {
    var fixture = Fixture{};
    var pools = fixture.pools();
    pools.topology.organic_group_count = 0;
    try std.testing.expectError(error.PondSolidOrganicTransferDimensionMismatch, transfer(0, 0, 1, 1, pools));
}

test "REDIST solid organic validation is atomic" {
    var fixture = Fixture{};
    fixture.soluble_p[48] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSolidOrganicTransferInput, transfer(0, 0, 1, 1, fixture.pools()));
    try std.testing.expectEqual(@as(f64, 1), fixture.total[1]);
}
