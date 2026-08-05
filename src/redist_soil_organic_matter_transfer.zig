const std = @import("std");

pub const Topology = struct {
    layer_count: usize,
    organic_fraction_count: usize = 6,
    organic_group_count: usize = 7,
    organic_component_count: usize = 3,
    residue_fraction_count: usize = 5,
    residue_component_count: usize = 2,
    soluble_component_count: usize = 5,
};
pub const Pools = struct {
    topology: Topology,
    layer_volume_m3: []const f64,
    preferential_fraction: []const f64,
    organic_carbon_g_c: []f64,
    organic_nitrogen_g_n: []f64,
    organic_phosphorus_g_p: []f64,
    residue_carbon_g_c: []f64,
    residue_nitrogen_g_n: []f64,
    residue_phosphorus_g_p: []f64,
    dissolved_carbon_g_c: []f64,
    dissolved_nitrogen_g_n: []f64,
    dissolved_phosphorus_g_p: []f64,
    dissolved_acetate_g_c: []f64,
    macropore_dissolved_carbon_g_c: []f64,
    macropore_dissolved_nitrogen_g_n: []f64,
    macropore_dissolved_phosphorus_g_p: []f64,
    macropore_dissolved_acetate_g_c: []f64,
    humus_carbon_g_c: []f64,
    humus_nitrogen_g_n: []f64,
    humus_phosphorus_g_p: []f64,
    humus_acetate_g_c: []f64,
    soluble_carbon_g_c: []f64,
    soluble_acetate_g_c: []f64,
    soluble_nitrogen_g_n: []f64,
    soluble_phosphorus_g_p: []f64,
};
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
fn validateMove(values: []const f64, source: usize, destination: usize, fraction: f64) !void {
    const moved = fraction * values[source];
    if (!std.math.isFinite(moved) or !std.math.isFinite(values[destination] + moved) or !std.math.isFinite(values[source] - moved)) return error.NonFiniteSoilOrganicMatterTransferResult;
}
fn move(values: []f64, source: usize, destination: usize, fraction: f64) void {
    const moved = fraction * values[source];
    values[destination] = values[destination] + moved;
    values[source] = values[source] - moved;
}

fn validateDimensions(p: Pools) !void {
    const t = p.topology;
    if (t.layer_count == 0 or t.organic_fraction_count != 6 or t.organic_group_count != 7 or t.organic_component_count != 3 or t.residue_fraction_count != 5 or t.residue_component_count != 2 or t.soluble_component_count != 5) return error.SoilOrganicMatterTransferDimensionMismatch;
    const organic_len = 6 * 7 * 3 * t.layer_count;
    const residue_len = 5 * 2 * t.layer_count;
    const fraction_len = 5 * t.layer_count;
    const soluble_len = 5 * 5 * t.layer_count;
    if (p.layer_volume_m3.len != t.layer_count or p.preferential_fraction.len != t.layer_count) return error.SoilOrganicMatterTransferDimensionMismatch;
    inline for (.{ p.organic_carbon_g_c, p.organic_nitrogen_g_n, p.organic_phosphorus_g_p }) |v| if (v.len != organic_len) return error.SoilOrganicMatterTransferDimensionMismatch;
    inline for (.{ p.residue_carbon_g_c, p.residue_nitrogen_g_n, p.residue_phosphorus_g_p }) |v| if (v.len != residue_len) return error.SoilOrganicMatterTransferDimensionMismatch;
    inline for (.{ p.dissolved_carbon_g_c, p.dissolved_nitrogen_g_n, p.dissolved_phosphorus_g_p, p.dissolved_acetate_g_c, p.macropore_dissolved_carbon_g_c, p.macropore_dissolved_nitrogen_g_n, p.macropore_dissolved_phosphorus_g_p, p.macropore_dissolved_acetate_g_c, p.humus_carbon_g_c, p.humus_nitrogen_g_n, p.humus_phosphorus_g_p, p.humus_acetate_g_c }) |v| if (v.len != fraction_len) return error.SoilOrganicMatterTransferDimensionMismatch;
    inline for (.{ p.soluble_carbon_g_c, p.soluble_acetate_g_c, p.soluble_nitrogen_g_n, p.soluble_phosphorus_g_p }) |v| if (v.len != soluble_len) return error.SoilOrganicMatterTransferDimensionMismatch;
}

/// Direct translation of REDIST 10322--10401.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, minimum_volume_m3: f64, zero_tolerance: f64, p: Pools) !void {
    try validateDimensions(p);
    const t = p.topology;
    if (source_layer >= t.layer_count or destination_layer >= t.layer_count or source_layer == destination_layer) return error.SoilOrganicMatterTransferDimensionMismatch;
    inline for (.{ p.layer_volume_m3, p.preferential_fraction, p.organic_carbon_g_c, p.organic_nitrogen_g_n, p.organic_phosphorus_g_p, p.residue_carbon_g_c, p.residue_nitrogen_g_n, p.residue_phosphorus_g_p, p.dissolved_carbon_g_c, p.dissolved_nitrogen_g_n, p.dissolved_phosphorus_g_p, p.dissolved_acetate_g_c, p.macropore_dissolved_carbon_g_c, p.macropore_dissolved_nitrogen_g_n, p.macropore_dissolved_phosphorus_g_p, p.macropore_dissolved_acetate_g_c, p.humus_carbon_g_c, p.humus_nitrogen_g_n, p.humus_phosphorus_g_p, p.humus_acetate_g_c, p.soluble_carbon_g_c, p.soluble_acetate_g_c, p.soluble_nitrogen_g_n, p.soluble_phosphorus_g_p }) |v| if (!finiteSlice(v)) return error.InvalidSoilOrganicMatterTransferInput;
    if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1 or
        !std.math.isFinite(minimum_volume_m3) or minimum_volume_m3 < 0 or
        !std.math.isFinite(zero_tolerance) or zero_tolerance < 0)
        return error.InvalidSoilOrganicMatterTransferInput;
    if (source_layer == 0 or p.layer_volume_m3[source_layer] <= minimum_volume_m3 or p.layer_volume_m3[destination_layer] <= minimum_volume_m3) return;
    const move_macro = p.preferential_fraction[source_layer] > zero_tolerance and p.preferential_fraction[destination_layer] > zero_tolerance;
    for (0..6) |k| for (0..7) |n| for (0..3) |m| {
        const s = organicIndex(t, k, n, m, source_layer);
        const d = organicIndex(t, k, n, m, destination_layer);
        inline for (.{ p.organic_carbon_g_c, p.organic_nitrogen_g_n, p.organic_phosphorus_g_p }) |v| try validateMove(v, s, d, fraction);
    };
    for (0..5) |k| {
        for (0..2) |m| {
            const s = residueIndex(t, k, m, source_layer);
            const d = residueIndex(t, k, m, destination_layer);
            inline for (.{ p.residue_carbon_g_c, p.residue_nitrogen_g_n, p.residue_phosphorus_g_p }) |v| try validateMove(v, s, d, fraction);
        }
        const s = fractionIndex(t, k, source_layer);
        const d = fractionIndex(t, k, destination_layer);
        inline for (.{ p.dissolved_carbon_g_c, p.dissolved_nitrogen_g_n, p.dissolved_phosphorus_g_p, p.dissolved_acetate_g_c }) |v| try validateMove(v, s, d, fraction);
        if (move_macro) inline for (.{ p.macropore_dissolved_carbon_g_c, p.macropore_dissolved_nitrogen_g_n, p.macropore_dissolved_phosphorus_g_p, p.macropore_dissolved_acetate_g_c }) |v| try validateMove(v, s, d, fraction);
        inline for (.{ p.humus_carbon_g_c, p.humus_nitrogen_g_n, p.humus_phosphorus_g_p, p.humus_acetate_g_c }) |v| try validateMove(v, s, d, fraction);
        for (0..5) |m| {
            const ss = solubleIndex(t, k, m, source_layer);
            const dd = solubleIndex(t, k, m, destination_layer);
            inline for (.{ p.soluble_carbon_g_c, p.soluble_acetate_g_c, p.soluble_nitrogen_g_n, p.soluble_phosphorus_g_p }) |v| try validateMove(v, ss, dd, fraction);
        }
    }
    for (0..6) |k| for (0..7) |n| for (0..3) |m| {
        const s = organicIndex(t, k, n, m, source_layer);
        const d = organicIndex(t, k, n, m, destination_layer);
        inline for (.{ p.organic_carbon_g_c, p.organic_nitrogen_g_n, p.organic_phosphorus_g_p }) |v| move(v, s, d, fraction);
    };
    for (0..5) |k| {
        for (0..2) |m| {
            const s = residueIndex(t, k, m, source_layer);
            const d = residueIndex(t, k, m, destination_layer);
            inline for (.{ p.residue_carbon_g_c, p.residue_nitrogen_g_n, p.residue_phosphorus_g_p }) |v| move(v, s, d, fraction);
        }
        const s = fractionIndex(t, k, source_layer);
        const d = fractionIndex(t, k, destination_layer);
        inline for (.{ p.dissolved_carbon_g_c, p.dissolved_nitrogen_g_n, p.dissolved_phosphorus_g_p, p.dissolved_acetate_g_c }) |v| move(v, s, d, fraction);
        if (move_macro) inline for (.{ p.macropore_dissolved_carbon_g_c, p.macropore_dissolved_nitrogen_g_n, p.macropore_dissolved_phosphorus_g_p, p.macropore_dissolved_acetate_g_c }) |v| move(v, s, d, fraction);
        inline for (.{ p.humus_carbon_g_c, p.humus_nitrogen_g_n, p.humus_phosphorus_g_p, p.humus_acetate_g_c }) |v| move(v, s, d, fraction);
        for (0..5) |m| {
            const ss = solubleIndex(t, k, m, source_layer);
            const dd = solubleIndex(t, k, m, destination_layer);
            inline for (.{ p.soluble_carbon_g_c, p.soluble_acetate_g_c, p.soluble_nitrogen_g_n, p.soluble_phosphorus_g_p }) |v| move(v, ss, dd, fraction);
        }
    }
}

const Fixture = struct {
    volume: [2]f64 = .{ 1, 1 },
    pref: [2]f64 = .{ 1, 1 },
    organic: [3][252]f64 = .{@as([252]f64, @splat(2))} ** 3,
    residue: [3][20]f64 = .{@as([20]f64, @splat(2))} ** 3,
    dissolved: [4][10]f64 = .{@as([10]f64, @splat(2))} ** 4,
    macro: [4][10]f64 = .{@as([10]f64, @splat(2))} ** 4,
    humus: [4][10]f64 = .{@as([10]f64, @splat(2))} ** 4,
    soluble: [4][50]f64 = .{@as([50]f64, @splat(2))} ** 4,
    fn pools(s: *Fixture) Pools {
        return .{ .topology = .{ .layer_count = 2 }, .layer_volume_m3 = &s.volume, .preferential_fraction = &s.pref, .organic_carbon_g_c = &s.organic[0], .organic_nitrogen_g_n = &s.organic[1], .organic_phosphorus_g_p = &s.organic[2], .residue_carbon_g_c = &s.residue[0], .residue_nitrogen_g_n = &s.residue[1], .residue_phosphorus_g_p = &s.residue[2], .dissolved_carbon_g_c = &s.dissolved[0], .dissolved_nitrogen_g_n = &s.dissolved[1], .dissolved_phosphorus_g_p = &s.dissolved[2], .dissolved_acetate_g_c = &s.dissolved[3], .macropore_dissolved_carbon_g_c = &s.macro[0], .macropore_dissolved_nitrogen_g_n = &s.macro[1], .macropore_dissolved_phosphorus_g_p = &s.macro[2], .macropore_dissolved_acetate_g_c = &s.macro[3], .humus_carbon_g_c = &s.humus[0], .humus_nitrogen_g_n = &s.humus[1], .humus_phosphorus_g_p = &s.humus[2], .humus_acetate_g_c = &s.humus[3], .soluble_carbon_g_c = &s.soluble[0], .soluble_acetate_g_c = &s.soluble[1], .soluble_nitrogen_g_n = &s.soluble[2], .soluble_phosphorus_g_p = &s.soluble[3] };
    }
};
test "REDIST soil organic transfers all fixed topologies and conserves elements" {
    var f = Fixture{};
    try transfer(1, 0, 0.25, 0, 0, f.pools());
    try std.testing.expectEqual(@as(f64, 1.5), f.organic[0][1]);
    try std.testing.expectEqual(@as(f64, 2.5), f.organic[0][0]);
    try std.testing.expectEqual(@as(f64, 1.5), f.macro[3][9]);
    try std.testing.expectEqual(@as(f64, 4), f.soluble[3][48] + f.soluble[3][49]);
}
test "REDIST soil organic VOLX and FHOL gates are exact" {
    var f = Fixture{};
    f.pref[0] = 0;
    try transfer(1, 0, 0.5, 0, 0, f.pools());
    try std.testing.expectEqual(@as(f64, 2), f.macro[0][1]);
    try std.testing.expectEqual(@as(f64, 1), f.humus[0][1]);
}
test "REDIST soil organic validation is atomic" {
    var f = Fixture{};
    f.soluble[3][49] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilOrganicMatterTransferInput, transfer(1, 0, 0.5, 0, 0, f.pools()));
    try std.testing.expectEqual(@as(f64, 2), f.organic[0][1]);
}
