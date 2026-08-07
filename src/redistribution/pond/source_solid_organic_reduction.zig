const std = @import("std");
const destination_transfer = @import("solid_organic_transfer.zig");

pub const Topology = destination_transfer.Topology;
pub const Pools = destination_transfer.Pools;

fn organicIndex(t: Topology, k: usize, n: usize, m: usize, layer: usize) usize {
    return (((k * t.organic_group_count + n) * t.organic_component_count + m) * t.layer_count) + layer;
}
fn residueIndex(t: Topology, k: usize, m: usize, layer: usize) usize {
    return ((k * t.residue_component_count + m) * t.layer_count) + layer;
}
fn fractionIndex(t: Topology, k: usize, layer: usize) usize {
    return k * t.layer_count + layer;
}
fn solubleIndex(t: Topology, k: usize, m: usize, layer: usize) usize {
    return ((k * t.soluble_component_count + m) * t.layer_count) + layer;
}
fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}
fn validateScale(value: f64, remaining_fraction: f64) !void {
    if (!std.math.isFinite(remaining_fraction * value)) return error.NonFinitePondSourceSolidOrganicReductionResult;
}
fn scale(values: []f64, at: usize, remaining_fraction: f64) void {
    values[at] = remaining_fraction * values[at];
}

fn validateDimensions(pools: Pools) !void {
    const t = pools.topology;
    if (t.layer_count == 0 or t.organic_fraction_count != 6 or t.organic_group_count != 7 or
        t.organic_component_count != 3 or t.residue_fraction_count != 5 or
        t.residue_component_count != 2 or t.soluble_component_count != 5)
        return error.PondSourceSolidOrganicReductionDimensionMismatch;
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
        return error.PondSourceSolidOrganicReductionDimensionMismatch;
}

/// Direct translation of REDIST 9374--9400 under the enclosing `FX == 1.0`.
pub fn reduce(surface_reappearance_flag: u8, source_layer: usize, redistribution_fraction: f64, remaining_fraction: f64, pools: Pools) !void {
    try validateDimensions(pools);
    const t = pools.topology;
    if (source_layer >= t.layer_count) return error.PondSourceSolidOrganicReductionDimensionMismatch;
    inline for (.{ pools.total_organic_carbon_g_c, pools.organic_carbon_g_c, pools.organic_nitrogen_g_n, pools.organic_phosphorus_g_p, pools.residue_carbon_g_c, pools.residue_nitrogen_g_n, pools.residue_phosphorus_g_p, pools.humus_carbon_g_c, pools.humus_nitrogen_g_n, pools.humus_phosphorus_g_p, pools.humus_acetate_g_c, pools.soluble_carbon_g_c, pools.soluble_acetate_g_c, pools.soluble_nitrogen_g_n, pools.soluble_phosphorus_g_p }) |values|
        if (!finiteSlice(values)) return error.InvalidPondSourceSolidOrganicReductionInput;
    if (!std.math.isFinite(redistribution_fraction) or redistribution_fraction < 0 or redistribution_fraction > 1 or
        !std.math.isFinite(remaining_fraction) or remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceSolidOrganicReductionInput;
    if (redistribution_fraction != 1.0 or surface_reappearance_flag != 0) return;

    try validateScale(pools.total_organic_carbon_g_c[source_layer], remaining_fraction);
    for (0..t.organic_fraction_count) |k| for (0..t.organic_group_count) |n| for (0..t.organic_component_count) |m| {
        const at = organicIndex(t, k, n, m, source_layer);
        inline for (.{ pools.organic_carbon_g_c, pools.organic_nitrogen_g_n, pools.organic_phosphorus_g_p }) |values|
            try validateScale(values[at], remaining_fraction);
    };
    for (0..t.residue_fraction_count) |k| {
        for (0..t.residue_component_count) |m| {
            const at = residueIndex(t, k, m, source_layer);
            inline for (.{ pools.residue_carbon_g_c, pools.residue_nitrogen_g_n, pools.residue_phosphorus_g_p }) |values|
                try validateScale(values[at], remaining_fraction);
        }
        const fraction_at = fractionIndex(t, k, source_layer);
        inline for (.{ pools.humus_carbon_g_c, pools.humus_nitrogen_g_n, pools.humus_phosphorus_g_p, pools.humus_acetate_g_c }) |values|
            try validateScale(values[fraction_at], remaining_fraction);
        for (0..t.soluble_component_count) |m| {
            const at = solubleIndex(t, k, m, source_layer);
            inline for (.{ pools.soluble_carbon_g_c, pools.soluble_acetate_g_c, pools.soluble_nitrogen_g_n, pools.soluble_phosphorus_g_p }) |values|
                try validateScale(values[at], remaining_fraction);
        }
    }

    scale(pools.total_organic_carbon_g_c, source_layer, remaining_fraction);
    for (0..t.organic_fraction_count) |k| for (0..t.organic_group_count) |n| for (0..t.organic_component_count) |m| {
        const at = organicIndex(t, k, n, m, source_layer);
        inline for (.{ pools.organic_carbon_g_c, pools.organic_nitrogen_g_n, pools.organic_phosphorus_g_p }) |values| scale(values, at, remaining_fraction);
    };
    for (0..t.residue_fraction_count) |k| {
        for (0..t.residue_component_count) |m| {
            const at = residueIndex(t, k, m, source_layer);
            inline for (.{ pools.residue_carbon_g_c, pools.residue_nitrogen_g_n, pools.residue_phosphorus_g_p }) |values| scale(values, at, remaining_fraction);
        }
        const fraction_at = fractionIndex(t, k, source_layer);
        inline for (.{ pools.humus_carbon_g_c, pools.humus_nitrogen_g_n, pools.humus_phosphorus_g_p, pools.humus_acetate_g_c }) |values| scale(values, fraction_at, remaining_fraction);
        for (0..t.soluble_component_count) |m| {
            const at = solubleIndex(t, k, m, source_layer);
            inline for (.{ pools.soluble_carbon_g_c, pools.soluble_acetate_g_c, pools.soluble_nitrogen_g_n, pools.soluble_phosphorus_g_p }) |values| scale(values, at, remaining_fraction);
        }
    }
}

const Fixture = struct {
    total: [2]f64 = @splat(4),
    organic_c: [252]f64 = @splat(4),
    organic_n: [252]f64 = @splat(4),
    organic_p: [252]f64 = @splat(4),
    residue_c: [20]f64 = @splat(4),
    residue_n: [20]f64 = @splat(4),
    residue_p: [20]f64 = @splat(4),
    humus_c: [10]f64 = @splat(4),
    humus_n: [10]f64 = @splat(4),
    humus_p: [10]f64 = @splat(4),
    humus_a: [10]f64 = @splat(4),
    soluble_c: [50]f64 = @splat(4),
    soluble_a: [50]f64 = @splat(4),
    soluble_n: [50]f64 = @splat(4),
    soluble_p: [50]f64 = @splat(4),
    fn pools(self: *Fixture) Pools {
        return .{ .topology = .{ .layer_count = 2 }, .total_organic_carbon_g_c = &self.total, .organic_carbon_g_c = &self.organic_c, .organic_nitrogen_g_n = &self.organic_n, .organic_phosphorus_g_p = &self.organic_p, .residue_carbon_g_c = &self.residue_c, .residue_nitrogen_g_n = &self.residue_n, .residue_phosphorus_g_p = &self.residue_p, .humus_carbon_g_c = &self.humus_c, .humus_nitrogen_g_n = &self.humus_n, .humus_phosphorus_g_p = &self.humus_p, .humus_acetate_g_c = &self.humus_a, .soluble_carbon_g_c = &self.soluble_c, .soluble_acetate_g_c = &self.soluble_a, .soluble_nitrogen_g_n = &self.soluble_n, .soluble_phosphorus_g_p = &self.soluble_p };
    }
};

test "REDIST source solid organic scales every fixed-loop pool and source layer zero" {
    var fixture = Fixture{};
    try reduce(0, 0, 1, 0.25, fixture.pools());
    try std.testing.expectEqual(@as(f64, 1), fixture.total[0]);
    try std.testing.expectEqual(@as(f64, 1), fixture.organic_p[250]);
    try std.testing.expectEqual(@as(f64, 1), fixture.residue_n[18]);
    try std.testing.expectEqual(@as(f64, 1), fixture.humus_a[8]);
    try std.testing.expectEqual(@as(f64, 1), fixture.soluble_p[48]);
    try std.testing.expectEqual(@as(f64, 4), fixture.total[1]);
}

test "REDIST source solid organic obeys exact redistribution and reappearance guards" {
    var fixture = Fixture{};
    try reduce(1, 0, 1, 0, fixture.pools());
    try reduce(0, 0, 0.999999999, 0, fixture.pools());
    try std.testing.expectEqual(@as(f64, 4), fixture.total[0]);
}

test "REDIST source solid organic validation is atomic" {
    var fixture = Fixture{};
    fixture.soluble_p[48] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceSolidOrganicReductionInput, reduce(0, 0, 1, 0.5, fixture.pools()));
    try std.testing.expectEqual(@as(f64, 4), fixture.total[0]);
}
