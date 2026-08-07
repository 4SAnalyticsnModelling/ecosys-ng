const std = @import("std");
const partition_module = @import("combustion_partition.zig");

pub const Partition = @import("microbial_biomass_combustion.zig").Partition;
pub const Rates = struct {
    layer_count: usize,
    dissolved_organic_carbon_g_c: []const f64,
    dissolved_organic_acetate_g_c: []const f64,
    dissolved_organic_nitrogen_g_n: []const f64,
    dissolved_organic_phosphorus_g_p: []const f64,
};
pub const State = struct {
    /// OSC(5,K,L), K=0..4, fraction-major then layer.
    charcoal_carbon_g_c: []f64,
    ammonium_g_n: []f64,
    dihydrogen_phosphate_g_p: []f64,
    accumulators: *partition_module.Accumulators,
};

fn index(layer_count: usize, k: usize, layer: usize) usize {
    return k * layer_count + layer;
}

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

fn validPartition(partition: Partition) bool {
    inline for (std.meta.fields(Partition)) |field| {
        const value = @field(partition, field.name);
        if (!std.math.isFinite(value) or value < 0 or value > 1) return false;
    }
    return (partition.charcoal_fraction >= 1 or @abs(partition.aerobic_fraction + partition.anaerobic_fraction - 1) <= 1e-12) and
        @abs(partition.mineral_n_fraction + partition.gaseous_n_fraction - 1) <= 1e-12 and
        @abs(partition.mineral_p_fraction + partition.gaseous_p_fraction - 1) <= 1e-12;
}

/// Direct translation of REDIST 10800--10812 for one active combustion layer.
pub fn combust(layer: usize, partition: Partition, rates: Rates, state: State) !void {
    const rate_len = 5 * rates.layer_count;
    if (rates.layer_count == 0 or layer >= rates.layer_count or
        rates.dissolved_organic_carbon_g_c.len != rate_len or
        rates.dissolved_organic_acetate_g_c.len != rate_len or
        rates.dissolved_organic_nitrogen_g_n.len != rate_len or
        rates.dissolved_organic_phosphorus_g_p.len != rate_len or
        state.charcoal_carbon_g_c.len != rate_len or
        state.ammonium_g_n.len != rates.layer_count or
        state.dihydrogen_phosphate_g_p.len != rates.layer_count)
        return error.FireDissolvedOrganicDimensionMismatch;
    if (!validPartition(partition) or
        !finiteSlice(rates.dissolved_organic_carbon_g_c) or
        !finiteSlice(rates.dissolved_organic_acetate_g_c) or
        !finiteSlice(rates.dissolved_organic_nitrogen_g_n) or
        !finiteSlice(rates.dissolved_organic_phosphorus_g_p) or
        !finiteSlice(state.charcoal_carbon_g_c) or
        !finiteSlice(state.ammonium_g_n) or
        !finiteSlice(state.dihydrogen_phosphate_g_p))
        return error.InvalidFireDissolvedOrganicInput;
    inline for (.{ rates.dissolved_organic_carbon_g_c, rates.dissolved_organic_acetate_g_c, rates.dissolved_organic_nitrogen_g_n, rates.dissolved_organic_phosphorus_g_p, state.charcoal_carbon_g_c, state.ammonium_g_n, state.dihydrogen_phosphate_g_p }) |values|
        for (values) |value| if (value < 0) return error.InvalidFireDissolvedOrganicInput;

    var next_charcoal: [5]f64 = undefined;
    for (0..5) |k| next_charcoal[k] = state.charcoal_carbon_g_c[index(rates.layer_count, k, layer)];
    var next_nh4 = state.ammonium_g_n[layer];
    var next_po4 = state.dihydrogen_phosphate_g_p[layer];
    var next = state.accumulators.*;
    for (0..5) |k| {
        const at = index(rates.layer_count, k, layer);
        const carbon = rates.dissolved_organic_carbon_g_c[at] + rates.dissolved_organic_acetate_g_c[at];
        const nitrogen = rates.dissolved_organic_nitrogen_g_n[at];
        const phosphorus = rates.dissolved_organic_phosphorus_g_p[at];
        next_charcoal[k] = next_charcoal[k] + carbon * partition.charcoal_fraction;
        next_nh4 = next_nh4 + nitrogen * partition.mineral_n_fraction;
        next_po4 = next_po4 + phosphorus * partition.mineral_p_fraction;
        next.carbon_dioxide_emission_g_c = next.carbon_dioxide_emission_g_c + carbon * (1.0 - partition.charcoal_fraction) * partition.aerobic_fraction;
        next.methane_emission_g_c = next.methane_emission_g_c + carbon * (1.0 - partition.charcoal_fraction) * partition.anaerobic_fraction;
        next.charcoal_carbon_g_c = next.charcoal_carbon_g_c + carbon * partition.charcoal_fraction;
        next.gaseous_nitrogen_g_n = next.gaseous_nitrogen_g_n + nitrogen * partition.gaseous_n_fraction;
        next.gaseous_phosphorus_g_p = next.gaseous_phosphorus_g_p + phosphorus * partition.gaseous_p_fraction;
        next.mineral_nitrogen_g_n = next.mineral_nitrogen_g_n + nitrogen * partition.mineral_n_fraction;
        next.mineral_phosphorus_g_p = next.mineral_phosphorus_g_p + phosphorus * partition.mineral_p_fraction;
        inline for (.{ next_charcoal[k], next_nh4, next_po4, next.carbon_dioxide_emission_g_c, next.methane_emission_g_c, next.charcoal_carbon_g_c, next.gaseous_nitrogen_g_n, next.gaseous_phosphorus_g_p, next.mineral_nitrogen_g_n, next.mineral_phosphorus_g_p }) |value|
            if (!std.math.isFinite(value)) return error.NonFiniteFireDissolvedOrganicResult;
    }
    for (0..5) |k| state.charcoal_carbon_g_c[index(rates.layer_count, k, layer)] = next_charcoal[k];
    state.ammonium_g_n[layer] = next_nh4;
    state.dihydrogen_phosphate_g_p[layer] = next_po4;
    state.accumulators.* = next;
}

fn zeroAccumulators() partition_module.Accumulators {
    return .{ .carbon_dioxide_emission_g_c = 0, .methane_emission_g_c = 0, .charcoal_carbon_g_c = 0, .gaseous_nitrogen_g_n = 0, .gaseous_phosphorus_g_p = 0, .mineral_nitrogen_g_n = 0, .mineral_phosphorus_g_p = 0 };
}

test "REDIST dissolved organic combustion conserves C N P across runtime layers" {
    const layers = 2;
    var doc: [5 * layers]f64 = @splat(0);
    var acetate: [doc.len]f64 = @splat(0);
    var nitrogen: [doc.len]f64 = @splat(0);
    var phosphorus: [doc.len]f64 = @splat(0);
    for (0..5) |k| {
        const at = index(layers, k, 1);
        doc[at] = 1;
        acetate[at] = 2;
        nitrogen[at] = 4;
        phosphorus[at] = 6;
    }
    var charcoal: [doc.len]f64 = @splat(0);
    var ammonium: [layers]f64 = @splat(0);
    var phosphate: [layers]f64 = @splat(0);
    var accumulators = zeroAccumulators();
    try combust(1, .{ .charcoal_fraction = 0.25, .aerobic_fraction = 0.6, .anaerobic_fraction = 0.4, .mineral_n_fraction = 0.7, .gaseous_n_fraction = 0.3, .mineral_p_fraction = 0.8, .gaseous_p_fraction = 0.2 }, .{ .layer_count = layers, .dissolved_organic_carbon_g_c = &doc, .dissolved_organic_acetate_g_c = &acetate, .dissolved_organic_nitrogen_g_n = &nitrogen, .dissolved_organic_phosphorus_g_p = &phosphorus }, .{ .charcoal_carbon_g_c = &charcoal, .ammonium_g_n = &ammonium, .dihydrogen_phosphate_g_p = &phosphate, .accumulators = &accumulators });
    try std.testing.expectApproxEqAbs(15, accumulators.charcoal_carbon_g_c + accumulators.carbon_dioxide_emission_g_c + accumulators.methane_emission_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(20, accumulators.mineral_nitrogen_g_n + accumulators.gaseous_nitrogen_g_n, 1e-12);
    try std.testing.expectApproxEqAbs(30, accumulators.mineral_phosphorus_g_p + accumulators.gaseous_phosphorus_g_p, 1e-12);
    try std.testing.expectEqual(@as(f64, 0), charcoal[0]);
    for (0..5) |k| try std.testing.expectEqual(@as(f64, 0.75), charcoal[index(layers, k, 1)]);
}

test "REDIST dissolved organic validation failure is atomic" {
    var doc: [5]f64 = @splat(1);
    var acetate: [5]f64 = @splat(1);
    var nitrogen: [5]f64 = @splat(1);
    var phosphorus: [5]f64 = @splat(1);
    phosphorus[4] = std.math.inf(f64);
    var charcoal: [5]f64 = @splat(4);
    var ammonium = [1]f64{5};
    var phosphate = [1]f64{6};
    var accumulators = zeroAccumulators();
    try std.testing.expectError(error.InvalidFireDissolvedOrganicInput, combust(0, .{ .charcoal_fraction = 0.25, .aerobic_fraction = 0.6, .anaerobic_fraction = 0.4, .mineral_n_fraction = 0.7, .gaseous_n_fraction = 0.3, .mineral_p_fraction = 0.8, .gaseous_p_fraction = 0.2 }, .{ .layer_count = 1, .dissolved_organic_carbon_g_c = &doc, .dissolved_organic_acetate_g_c = &acetate, .dissolved_organic_nitrogen_g_n = &nitrogen, .dissolved_organic_phosphorus_g_p = &phosphorus }, .{ .charcoal_carbon_g_c = &charcoal, .ammonium_g_n = &ammonium, .dihydrogen_phosphate_g_p = &phosphate, .accumulators = &accumulators }));
    try std.testing.expectEqual(@as(f64, 4), charcoal[0]);
    try std.testing.expectEqual(@as(f64, 5), ammonium[0]);
    try std.testing.expectEqual(@as(f64, 0), accumulators.charcoal_carbon_g_c);
}
