const std = @import("std");
const snow = @import("snow_solute_transport.zig");

pub const Report = struct {
    transfers: usize,
    maximum_fraction: f64,
};

/// REDIST snow relayering for runtime layer counts. Every physical constituent,
/// sensible energy, and tracked solute moves by the same source fraction. The
/// complete candidate state is validated before atomic publication.
pub fn apply(allocator: std.mem.Allocator, state: *snow.State, negligible_volume_m3: f64) !Report {
    if (!std.math.isFinite(negligible_volume_m3) or negligible_volume_m3 < 0) return error.InvalidSnowRelayeringTolerance;
    const solid = try allocator.dupe(f64, state.solid_snow_water_equivalent_m3);
    defer allocator.free(solid);
    const liquid = try allocator.dupe(f64, state.liquid_water_volume_m3);
    defer allocator.free(liquid);
    const vapor = try allocator.dupe(f64, state.vapor_water_equivalent_m3);
    defer allocator.free(vapor);
    const ice = try allocator.dupe(f64, state.ice_volume_m3);
    defer allocator.free(ice);
    const temperature = try allocator.dupe(f64, state.temperature_k);
    defer allocator.free(temperature);
    const heat_capacity = try allocator.dupe(f64, state.heat_capacity_megajoules_per_k);
    defer allocator.free(heat_capacity);
    const amount_g = try allocator.dupe(f64, state.amount_g);
    defer allocator.free(amount_g);
    var report: Report = .{ .transfers = 0, .maximum_fraction = 0 };

    for (0..state.cell_count) |cell| {
        for (0..state.layer_capacity -| 1) |layer| {
            const current = cell * state.layer_capacity + layer;
            const lower = current + 1;
            const density = state.snow_density_megagrams_per_m3[current];
            const lower_density = state.snow_density_megagrams_per_m3[lower];
            inline for (.{ density, lower_density, state.target_layer_volume_m3[current], solid[current], liquid[current], vapor[current], ice[current], temperature[current], temperature[lower], heat_capacity[current], heat_capacity[lower] }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSnowRelayeringState;
            if (density <= 0 or lower_density <= 0 or state.target_layer_volume_m3[current] < 0 or solid[current] < 0 or liquid[current] < 0 or vapor[current] < 0 or ice[current] < 0 or temperature[current] <= 0 or temperature[lower] <= 0 or heat_capacity[current] < 0 or heat_capacity[lower] < 0) return error.InvalidSnowRelayeringState;
            const current_volume_m3 = solid[current] / density + liquid[current] + ice[current];
            const lower_volume_m3 = solid[lower] / lower_density + liquid[lower] + ice[lower];
            if (current_volume_m3 <= negligible_volume_m3) continue;
            const volume_deficit_m3 = state.target_layer_volume_m3[current] - current_volume_m3;
            var source: usize = undefined;
            var destination: usize = undefined;
            var fraction: f64 = 0;
            if (volume_deficit_m3 > negligible_volume_m3 and lower_volume_m3 > negligible_volume_m3) {
                source = lower;
                destination = current;
                fraction = @min(1, volume_deficit_m3 / lower_volume_m3);
            } else if (volume_deficit_m3 < -negligible_volume_m3 and current_volume_m3 >= state.target_layer_volume_m3[current]) {
                source = current;
                destination = lower;
                fraction = @min(1, -volume_deficit_m3 / current_volume_m3);
            } else continue;
            if (!std.math.isFinite(fraction) or fraction <= 0 or fraction > 1) return error.InvalidSnowRelayeringFraction;
            const retained = 1 - fraction;
            const source_energy_megajoules = heat_capacity[source] * temperature[source];
            const destination_energy_megajoules = heat_capacity[destination] * temperature[destination];
            solid[destination] += fraction * solid[source];
            liquid[destination] += fraction * liquid[source];
            vapor[destination] += fraction * vapor[source];
            ice[destination] += fraction * ice[source];
            solid[source] *= retained;
            liquid[source] *= retained;
            vapor[source] *= retained;
            ice[source] *= retained;
            heat_capacity[destination] = 2.095 * solid[destination] + 4.19 * (liquid[destination] + vapor[destination]) + 1.9274 * ice[destination];
            heat_capacity[source] = 2.095 * solid[source] + 4.19 * (liquid[source] + vapor[source]) + 1.9274 * ice[source];
            if (heat_capacity[destination] > negligible_volume_m3) temperature[destination] = (destination_energy_megajoules + fraction * source_energy_megajoules) / heat_capacity[destination] else temperature[destination] = temperature[source];
            if (heat_capacity[source] > negligible_volume_m3) temperature[source] = retained * source_energy_megajoules / heat_capacity[source] else temperature[source] = temperature[destination];
            for (0..snow.species_count) |species| {
                const source_amount = source * snow.species_count + species;
                const destination_amount = destination * snow.species_count + species;
                amount_g[destination_amount] += fraction * amount_g[source_amount];
                amount_g[source_amount] *= retained;
            }
            report.transfers += 1;
            report.maximum_fraction = @max(report.maximum_fraction, fraction);
        }
    }
    inline for (.{ solid, liquid, vapor, ice, temperature, heat_capacity, amount_g }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSnowRelayeringResult;
    @memcpy(state.solid_snow_water_equivalent_m3, solid);
    @memcpy(state.liquid_water_volume_m3, liquid);
    @memcpy(state.vapor_water_equivalent_m3, vapor);
    @memcpy(state.ice_volume_m3, ice);
    @memcpy(state.temperature_k, temperature);
    @memcpy(state.heat_capacity_megajoules_per_k, heat_capacity);
    @memcpy(state.amount_g, amount_g);
    for (state.active, solid, liquid, vapor, ice) |*active, solid_m3, liquid_m3, vapor_m3, ice_m3| active.* = solid_m3 + liquid_m3 + vapor_m3 + ice_m3 > negligible_volume_m3;
    state.refreshAllGeometry();
    return report;
}

test "runtime relayering conserves physical content solutes and energy" {
    var state = try snow.State.init(std.testing.allocator, 1, 3);
    defer state.deinit();
    try state.initializePhysicalState(&.{0.12}, &.{1}, &.{268}, &.{ 0.05, 0.10, 0.20 }, 0.05);
    // Force the first layer above its runtime target and give all layers solute.
    state.solid_snow_water_equivalent_m3[0] += 0.003;
    state.amount_g[0] = 2;
    state.amount_g[snow.species_count] = 3;
    state.refreshAllGeometry();
    var water_before: f64 = 0;
    var energy_before: f64 = 0;
    for (state.solid_snow_water_equivalent_m3, state.liquid_water_volume_m3, state.vapor_water_equivalent_m3, state.ice_volume_m3, state.heat_capacity_megajoules_per_k, state.temperature_k) |solid_m3, liquid_m3, vapor_m3, ice_m3, capacity, temperature_k| {
        water_before += solid_m3 + liquid_m3 + vapor_m3 + ice_m3;
        energy_before += capacity * temperature_k;
    }
    const report = try apply(std.testing.allocator, &state, 1e-15);
    var water_after: f64 = 0;
    var energy_after: f64 = 0;
    for (state.solid_snow_water_equivalent_m3, state.liquid_water_volume_m3, state.vapor_water_equivalent_m3, state.ice_volume_m3, state.heat_capacity_megajoules_per_k, state.temperature_k) |solid_m3, liquid_m3, vapor_m3, ice_m3, capacity, temperature_k| {
        water_after += solid_m3 + liquid_m3 + vapor_m3 + ice_m3;
        energy_after += capacity * temperature_k;
    }
    try std.testing.expect(report.transfers > 0);
    try std.testing.expectApproxEqAbs(water_before, water_after, 1e-14);
    try std.testing.expectApproxEqAbs(energy_before, energy_after, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5), state.amount_g[0] + state.amount_g[snow.species_count] + state.amount_g[2 * snow.species_count], 1e-14);
}
