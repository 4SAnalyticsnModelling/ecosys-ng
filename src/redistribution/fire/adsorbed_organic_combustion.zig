const std = @import("std");
const dissolved = @import("dissolved_organic_combustion.zig");
const partition_module = @import("combustion_partition.zig");

pub const Partition = dissolved.Partition;
pub const Rates = struct {
    layer_count: usize,
    adsorbed_organic_carbon_g_c: []const f64,
    adsorbed_organic_acetate_g_c: []const f64,
    adsorbed_organic_nitrogen_g_n: []const f64,
    adsorbed_organic_phosphorus_g_p: []const f64,
};
pub const State = dissolved.State;

/// Direct translation of REDIST 10825--10837 for one active combustion layer.
/// The legacy equations and K=0..4 operation order are identical to the
/// dissolved-organic combustion kernel, but consume the adsorbed pools.
pub fn combust(layer: usize, partition: Partition, rates: Rates, state: State) !void {
    dissolved.combust(layer, partition, .{
        .layer_count = rates.layer_count,
        .dissolved_organic_carbon_g_c = rates.adsorbed_organic_carbon_g_c,
        .dissolved_organic_acetate_g_c = rates.adsorbed_organic_acetate_g_c,
        .dissolved_organic_nitrogen_g_n = rates.adsorbed_organic_nitrogen_g_n,
        .dissolved_organic_phosphorus_g_p = rates.adsorbed_organic_phosphorus_g_p,
    }, state) catch |err| switch (err) {
        error.FireDissolvedOrganicDimensionMismatch => return error.FireAdsorbedOrganicDimensionMismatch,
        error.InvalidFireDissolvedOrganicInput => return error.InvalidFireAdsorbedOrganicInput,
        error.NonFiniteFireDissolvedOrganicResult => return error.NonFiniteFireAdsorbedOrganicResult,
    };
}

fn zeroAccumulators() partition_module.Accumulators {
    return .{ .carbon_dioxide_emission_g_c = 0, .methane_emission_g_c = 0, .charcoal_carbon_g_c = 0, .gaseous_nitrogen_g_n = 0, .gaseous_phosphorus_g_p = 0, .mineral_nitrogen_g_n = 0, .mineral_phosphorus_g_p = 0 };
}

test "REDIST adsorbed organic combustion conserves C N P in exact K order" {
    const layers = 2;
    var carbon: [5 * layers]f64 = @splat(0);
    var acetate: [carbon.len]f64 = @splat(0);
    var nitrogen: [carbon.len]f64 = @splat(0);
    var phosphorus: [carbon.len]f64 = @splat(0);
    for (0..5) |k| {
        const at = k * layers + 1;
        carbon[at] = 1;
        acetate[at] = 2;
        nitrogen[at] = 4;
        phosphorus[at] = 6;
    }
    var charcoal: [carbon.len]f64 = @splat(0);
    var ammonium: [layers]f64 = @splat(0);
    var phosphate: [layers]f64 = @splat(0);
    var accumulators = zeroAccumulators();
    try combust(1, .{ .charcoal_fraction = 0.25, .aerobic_fraction = 0.6, .anaerobic_fraction = 0.4, .mineral_n_fraction = 0.7, .gaseous_n_fraction = 0.3, .mineral_p_fraction = 0.8, .gaseous_p_fraction = 0.2 }, .{ .layer_count = layers, .adsorbed_organic_carbon_g_c = &carbon, .adsorbed_organic_acetate_g_c = &acetate, .adsorbed_organic_nitrogen_g_n = &nitrogen, .adsorbed_organic_phosphorus_g_p = &phosphorus }, .{ .charcoal_carbon_g_c = &charcoal, .ammonium_g_n = &ammonium, .dihydrogen_phosphate_g_p = &phosphate, .accumulators = &accumulators });
    try std.testing.expectApproxEqAbs(15, accumulators.charcoal_carbon_g_c + accumulators.carbon_dioxide_emission_g_c + accumulators.methane_emission_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(20, accumulators.mineral_nitrogen_g_n + accumulators.gaseous_nitrogen_g_n, 1e-12);
    try std.testing.expectApproxEqAbs(30, accumulators.mineral_phosphorus_g_p + accumulators.gaseous_phosphorus_g_p, 1e-12);
    try std.testing.expectEqual(@as(f64, 0), charcoal[0]);
    for (0..5) |k| try std.testing.expectEqual(@as(f64, 0.75), charcoal[k * layers + 1]);
}

test "REDIST adsorbed organic validation failure is atomic" {
    var carbon: [5]f64 = @splat(1);
    var acetate: [5]f64 = @splat(1);
    var nitrogen: [5]f64 = @splat(1);
    var phosphorus: [5]f64 = @splat(1);
    phosphorus[4] = std.math.inf(f64);
    var charcoal: [5]f64 = @splat(4);
    var ammonium = [1]f64{5};
    var phosphate = [1]f64{6};
    var accumulators = zeroAccumulators();
    try std.testing.expectError(error.InvalidFireAdsorbedOrganicInput, combust(0, .{ .charcoal_fraction = 0.25, .aerobic_fraction = 0.6, .anaerobic_fraction = 0.4, .mineral_n_fraction = 0.7, .gaseous_n_fraction = 0.3, .mineral_p_fraction = 0.8, .gaseous_p_fraction = 0.2 }, .{ .layer_count = 1, .adsorbed_organic_carbon_g_c = &carbon, .adsorbed_organic_acetate_g_c = &acetate, .adsorbed_organic_nitrogen_g_n = &nitrogen, .adsorbed_organic_phosphorus_g_p = &phosphorus }, .{ .charcoal_carbon_g_c = &charcoal, .ammonium_g_n = &ammonium, .dihydrogen_phosphate_g_p = &phosphate, .accumulators = &accumulators }));
    try std.testing.expectEqual(@as(f64, 4), charcoal[0]);
    try std.testing.expectEqual(@as(f64, 5), ammonium[0]);
    try std.testing.expectEqual(@as(f64, 0), accumulators.charcoal_carbon_g_c);
}
