const std = @import("std");
const partition_module = @import("combustion_partition.zig");

pub const Partition = struct { charcoal_fraction: f64, aerobic_fraction: f64, anaerobic_fraction: f64, mineral_n_fraction: f64, gaseous_n_fraction: f64, mineral_p_fraction: f64, gaseous_p_fraction: f64 };
pub const Rates = struct { layer_count: usize, carbon_g_c: []const f64, nitrogen_g_n: []const f64, phosphorus_g_p: []const f64 };
pub const State = struct {
    /// OSC(5,K,L), K=0..4, fraction-major then layer.
    charcoal_carbon_g_c: []f64,
    ammonium_g_n: []f64,
    dihydrogen_phosphate_g_p: []f64,
    accumulators: *partition_module.Accumulators,
};
fn rateIndex(layer_count: usize, k: usize, n: usize, m: usize, layer: usize) usize {
    return (((k * 7 + n) * 3 + m) * layer_count) + layer;
}
fn charcoalIndex(layer_count: usize, k: usize, layer: usize) usize {
    return k * layer_count + layer;
}
fn finiteSlice(v: []const f64) bool {
    for (v) |x| if (!std.math.isFinite(x)) return false;
    return true;
}
fn finitePartition(p: Partition) bool {
    inline for (std.meta.fields(Partition)) |f| if (!std.math.isFinite(@field(p, f.name))) return false;
    return true;
}

/// Direct translation of REDIST 10739--10761 for one active combustion layer.
pub fn combust(layer: usize, partition: Partition, rates: Rates, state: State) !void {
    const len = 6 * 7 * 3 * rates.layer_count;
    if (rates.layer_count == 0 or layer >= rates.layer_count or rates.carbon_g_c.len != len or rates.nitrogen_g_n.len != len or rates.phosphorus_g_p.len != len or state.charcoal_carbon_g_c.len != 5 * rates.layer_count or state.ammonium_g_n.len != rates.layer_count or state.dihydrogen_phosphate_g_p.len != rates.layer_count) return error.FireMicrobialBiomassDimensionMismatch;
    if (!finitePartition(partition) or !finiteSlice(rates.carbon_g_c) or !finiteSlice(rates.nitrogen_g_n) or !finiteSlice(rates.phosphorus_g_p) or !finiteSlice(state.charcoal_carbon_g_c) or !finiteSlice(state.ammonium_g_n) or !finiteSlice(state.dihydrogen_phosphate_g_p)) return error.InvalidFireMicrobialBiomassInput;
    inline for (std.meta.fields(Partition)) |field| {
        const value = @field(partition, field.name);
        if (value < 0 or value > 1) return error.InvalidFireMicrobialBiomassInput;
    }
    if ((partition.charcoal_fraction < 1 and @abs(partition.aerobic_fraction + partition.anaerobic_fraction - 1) > 1e-12) or
        @abs(partition.mineral_n_fraction + partition.gaseous_n_fraction - 1) > 1e-12 or
        @abs(partition.mineral_p_fraction + partition.gaseous_p_fraction - 1) > 1e-12)
        return error.InvalidFireMicrobialBiomassInput;
    inline for (.{ rates.carbon_g_c, rates.nitrogen_g_n, rates.phosphorus_g_p, state.charcoal_carbon_g_c, state.ammonium_g_n, state.dihydrogen_phosphate_g_p }) |values|
        for (values) |value| if (value < 0) return error.InvalidFireMicrobialBiomassInput;
    var next_charcoal: [5]f64 = undefined;
    for (0..5) |k| next_charcoal[k] = state.charcoal_carbon_g_c[charcoalIndex(rates.layer_count, k, layer)];
    var next_nh4 = state.ammonium_g_n[layer];
    var next_po4 = state.dihydrogen_phosphate_g_p[layer];
    var next = state.accumulators.*;
    for (0..6) |k| {
        if (layer == 0 and (k == 3 or k == 4)) continue;
        for (0..7) |n| for (0..3) |m| {
            const at = rateIndex(rates.layer_count, k, n, m, layer);
            const c = rates.carbon_g_c[at];
            const nitrogen = rates.nitrogen_g_n[at];
            const phosphorus = rates.phosphorus_g_p[at];
            const kk = if (k <= 4) k else 1;
            next_charcoal[kk] = next_charcoal[kk] + c * partition.charcoal_fraction;
            next_nh4 = next_nh4 + nitrogen * partition.mineral_n_fraction;
            next_po4 = next_po4 + phosphorus * partition.mineral_p_fraction;
            next.carbon_dioxide_emission_g_c = next.carbon_dioxide_emission_g_c + c * (1.0 - partition.charcoal_fraction) * partition.aerobic_fraction;
            next.methane_emission_g_c = next.methane_emission_g_c + c * (1.0 - partition.charcoal_fraction) * partition.anaerobic_fraction;
            next.charcoal_carbon_g_c = next.charcoal_carbon_g_c + c * partition.charcoal_fraction;
            next.gaseous_nitrogen_g_n = next.gaseous_nitrogen_g_n + nitrogen * partition.gaseous_n_fraction;
            next.gaseous_phosphorus_g_p = next.gaseous_phosphorus_g_p + phosphorus * partition.gaseous_p_fraction;
            next.mineral_nitrogen_g_n = next.mineral_nitrogen_g_n + nitrogen * partition.mineral_n_fraction;
            next.mineral_phosphorus_g_p = next.mineral_phosphorus_g_p + phosphorus * partition.mineral_p_fraction;
            inline for (.{ next_charcoal[kk], next_nh4, next_po4, next.carbon_dioxide_emission_g_c, next.methane_emission_g_c, next.charcoal_carbon_g_c, next.gaseous_nitrogen_g_n, next.gaseous_phosphorus_g_p, next.mineral_nitrogen_g_n, next.mineral_phosphorus_g_p }) |x| if (!std.math.isFinite(x)) return error.NonFiniteFireMicrobialBiomassResult;
        };
    }
    for (0..5) |k| state.charcoal_carbon_g_c[charcoalIndex(rates.layer_count, k, layer)] = next_charcoal[k];
    state.ammonium_g_n[layer] = next_nh4;
    state.dihydrogen_phosphate_g_p[layer] = next_po4;
    state.accumulators.* = next;
}

test "REDIST fire microbial biomass preserves loops remap and C N P conservation" {
    const layers = 2;
    var c: [6 * 7 * 3 * layers]f64 = @splat(0);
    var n: [c.len]f64 = @splat(0);
    var p: [c.len]f64 = @splat(0);
    for (0..6) |k| for (0..7) |group| for (0..3) |m| {
        const at = rateIndex(layers, k, group, m, 1);
        c[at] = 1;
        n[at] = 2;
        p[at] = 3;
    };
    var charcoal: [5 * layers]f64 = @splat(0);
    var nh4: [layers]f64 = @splat(0);
    var po4: [layers]f64 = @splat(0);
    var a = partition_module.Accumulators{ .carbon_dioxide_emission_g_c = 0, .methane_emission_g_c = 0, .charcoal_carbon_g_c = 0, .gaseous_nitrogen_g_n = 0, .gaseous_phosphorus_g_p = 0, .mineral_nitrogen_g_n = 0, .mineral_phosphorus_g_p = 0 };
    try combust(1, .{ .charcoal_fraction = 0.25, .aerobic_fraction = 0.6, .anaerobic_fraction = 0.4, .mineral_n_fraction = 0.7, .gaseous_n_fraction = 0.3, .mineral_p_fraction = 0.8, .gaseous_p_fraction = 0.2 }, .{ .layer_count = layers, .carbon_g_c = &c, .nitrogen_g_n = &n, .phosphorus_g_p = &p }, .{ .charcoal_carbon_g_c = &charcoal, .ammonium_g_n = &nh4, .dihydrogen_phosphate_g_p = &po4, .accumulators = &a });
    const count: f64 = 126;
    try std.testing.expectApproxEqAbs(count, a.charcoal_carbon_g_c + a.carbon_dioxide_emission_g_c + a.methane_emission_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(2 * count, a.mineral_nitrogen_g_n + a.gaseous_nitrogen_g_n, 1e-12);
    try std.testing.expectApproxEqAbs(3 * count, a.mineral_phosphorus_g_p + a.gaseous_phosphorus_g_p, 1e-12);
    try std.testing.expectEqual(@as(f64, 42 * 0.25), charcoal[charcoalIndex(layers, 1, 1)]);
}
test "REDIST fire microbial biomass surface exclusions are exact" {
    var c: [6 * 7 * 3]f64 = @splat(1);
    var n: [c.len]f64 = @splat(0);
    var p: [c.len]f64 = @splat(0);
    var charcoal: [5]f64 = @splat(0);
    var mineral: [1]f64 = @splat(0);
    var a = partition_module.Accumulators{ .carbon_dioxide_emission_g_c = 0, .methane_emission_g_c = 0, .charcoal_carbon_g_c = 0, .gaseous_nitrogen_g_n = 0, .gaseous_phosphorus_g_p = 0, .mineral_nitrogen_g_n = 0, .mineral_phosphorus_g_p = 0 };
    try combust(0, .{ .charcoal_fraction = 1, .aerobic_fraction = 0, .anaerobic_fraction = 0, .mineral_n_fraction = 1, .gaseous_n_fraction = 0, .mineral_p_fraction = 1, .gaseous_p_fraction = 0 }, .{ .layer_count = 1, .carbon_g_c = &c, .nitrogen_g_n = &n, .phosphorus_g_p = &p }, .{ .charcoal_carbon_g_c = &charcoal, .ammonium_g_n = &mineral, .dihydrogen_phosphate_g_p = &mineral, .accumulators = &a });
    try std.testing.expectEqual(@as(f64, 84), a.charcoal_carbon_g_c);
}
test "REDIST fire microbial biomass validation is atomic" {
    var c: [6 * 7 * 3]f64 = @splat(1);
    var n: [c.len]f64 = @splat(1);
    var p: [c.len]f64 = @splat(1);
    p[p.len - 1] = std.math.inf(f64);
    var charcoal: [5]f64 = @splat(0);
    var mineral: [1]f64 = @splat(0);
    var a = partition_module.Accumulators{ .carbon_dioxide_emission_g_c = 0, .methane_emission_g_c = 0, .charcoal_carbon_g_c = 0, .gaseous_nitrogen_g_n = 0, .gaseous_phosphorus_g_p = 0, .mineral_nitrogen_g_n = 0, .mineral_phosphorus_g_p = 0 };
    try std.testing.expectError(error.InvalidFireMicrobialBiomassInput, combust(0, .{ .charcoal_fraction = 1, .aerobic_fraction = 0, .anaerobic_fraction = 0, .mineral_n_fraction = 1, .gaseous_n_fraction = 0, .mineral_p_fraction = 1, .gaseous_p_fraction = 0 }, .{ .layer_count = 1, .carbon_g_c = &c, .nitrogen_g_n = &n, .phosphorus_g_p = &p }, .{ .charcoal_carbon_g_c = &charcoal, .ammonium_g_n = &mineral, .dihydrogen_phosphate_g_p = &mineral, .accumulators = &a }));
    try std.testing.expectEqual(@as(f64, 0), a.charcoal_carbon_g_c);
}
