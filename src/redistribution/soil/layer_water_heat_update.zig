const std = @import("std");
const balance_module = @import("layer_balance_initialization.zig");

pub const LayerState = struct {
    temperature_k: f64,
    heat_capacity_megajoules_k: f64,
    matrix_heat_capacity_megajoules_k: f64,
    matrix_water_m3: f64,
    water_vapor_m3: f64,
    mobile_matrix_water_m3: f64,
    matrix_ice_m3: f64,
    macropore_water_m3: f64,
    macropore_ice_m3: f64,
    previous_matrix_water_m3: f64,
    previous_macropore_water_m3: f64,
    previous_matrix_ice_m3: f64,
    previous_macropore_ice_m3: f64,
    matrix_air_capacity_m3: f64,
    macropore_air_capacity_m3: f64,
    air_filled_volume_m3: f64,
    bulk_density_megagrams_m3: f64,
    layer_thickness_m: f64,
    temperature_c: f64,
    maximum_temperature_c: f64,
    minimum_temperature_c: f64,
    water_change_m3: f64,
    ice_change_m3: f64,
};
pub const LayerFlux = struct {
    matrix_water_m3: f64,
    vapor_m3: f64,
    mobile_water_m3: f64,
    evaporation_condensation_m3: f64,
    matrix_freeze_thaw_m3: f64,
    matrix_macropore_exchange_m3: f64,
    root_water_uptake_m3: f64,
    subsurface_water_input_m3: f64,
    macropore_water_m3: f64,
    macropore_freeze_thaw_m3: f64,
    conductive_heat_megajoules: f64,
    matrix_latent_heat_megajoules: f64,
    vapor_latent_heat_megajoules: f64,
    root_convective_heat_megajoules: f64,
    subsurface_convective_heat_megajoules: f64,
    combustion_heat_megajoules: f64,
};
pub const Parameters = struct {
    negligible_water_m3: f64,
    ice_density_megagrams_m3: f64,
    substeps_per_hour_reciprocal: f64,
    liquid_heat_capacity_megajoules_m3_k: f64,
    ice_heat_capacity_megajoules_m3_k: f64,
    minimum_heat_capacity_megajoules_k: f64,
    minimum_layer_thickness_m: f64,
};
pub const Result = struct { balance: balance_module.SoilLayerBalance, heat_input_megajoules: f64 };

fn valid(value: anytype) bool {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |f| if (!std.math.isFinite(@field(value, f.name))) return false;
    return true;
}

/// Exact runtime-slice translation of redist.f lines 5943--6015.
pub fn updateLayers(states: []LayerState, fluxes: []const LayerFlux, initial_balance: balance_module.SoilLayerBalance, initial_heat_input_megajoules: f64, top_boundary_temperature_k: f64, p: Parameters) !Result {
    if (states.len == 0 or states.len != fluxes.len or !valid(initial_balance) or !std.math.isFinite(initial_heat_input_megajoules) or !std.math.isFinite(top_boundary_temperature_k) or !valid(p) or p.ice_density_megagrams_m3 == 0)
        return error.InvalidSoilLayerUpdateInput;
    var balance = initial_balance;
    var heat_input = initial_heat_input_megajoules;
    for (states, fluxes, 0..) |*s, f, layer| {
        if (!valid(s.*) or !valid(f)) return error.InvalidSoilLayerUpdateInput;
        const old_heat_capacity = s.heat_capacity_megajoules_k;
        s.matrix_water_m3 = s.matrix_water_m3 + f.matrix_water_m3 + f.evaporation_condensation_m3 + f.matrix_freeze_thaw_m3 + f.matrix_macropore_exchange_m3 + f.root_water_uptake_m3 + f.subsurface_water_input_m3;
        s.water_vapor_m3 = s.water_vapor_m3 + f.vapor_m3 - f.evaporation_condensation_m3;
        if (s.matrix_water_m3 > p.negligible_water_m3) {
            s.mobile_matrix_water_m3 = s.mobile_matrix_water_m3 + f.mobile_water_m3 + f.evaporation_condensation_m3 + f.matrix_freeze_thaw_m3 + f.matrix_macropore_exchange_m3 + f.root_water_uptake_m3 * @max(0.0, s.mobile_matrix_water_m3 / s.matrix_water_m3) + f.subsurface_water_input_m3;
            s.mobile_matrix_water_m3 = @min(s.matrix_water_m3, s.mobile_matrix_water_m3 + 2.88e-3 * (s.matrix_water_m3 - s.mobile_matrix_water_m3) * p.substeps_per_hour_reciprocal);
        } else s.mobile_matrix_water_m3 = 0;
        s.matrix_ice_m3 = s.matrix_ice_m3 - f.matrix_freeze_thaw_m3 / p.ice_density_megagrams_m3;
        s.macropore_water_m3 = s.macropore_water_m3 + f.macropore_water_m3 - f.matrix_macropore_exchange_m3 + f.macropore_freeze_thaw_m3;
        s.macropore_ice_m3 = s.macropore_ice_m3 - f.macropore_freeze_thaw_m3 / p.ice_density_megagrams_m3;
        s.water_change_m3 = s.previous_matrix_water_m3 + s.previous_macropore_water_m3 - s.matrix_water_m3 - s.macropore_water_m3;
        s.ice_change_m3 = s.previous_matrix_ice_m3 + s.previous_macropore_ice_m3 - s.matrix_ice_m3 - s.macropore_ice_m3;
        s.air_filled_volume_m3 = if (s.bulk_density_megagrams_m3 > 0) @max(0, s.matrix_air_capacity_m3 - s.matrix_water_m3 - s.matrix_ice_m3 + s.macropore_air_capacity_m3 - s.macropore_water_m3 - s.macropore_ice_m3) else 0;
        const energy = old_heat_capacity * s.temperature_k;
        s.heat_capacity_megajoules_k = s.matrix_heat_capacity_megajoules_k + p.liquid_heat_capacity_megajoules_m3_k * (s.matrix_water_m3 + s.water_vapor_m3 + s.macropore_water_m3) + p.ice_heat_capacity_megajoules_m3_k * (s.matrix_ice_m3 + s.macropore_ice_m3);
        balance.total_heat_capacity_megajoules_k += s.heat_capacity_megajoules_k;
        balance.total_matrix_heat_capacity_megajoules_k += s.matrix_heat_capacity_megajoules_k;
        balance.total_matrix_water_m3 += s.matrix_water_m3;
        balance.total_water_vapor_m3 += s.water_vapor_m3;
        balance.total_macropore_water_m3 += s.macropore_water_m3;
        balance.total_matrix_ice_m3 += s.matrix_ice_m3;
        balance.total_macropore_ice_m3 += s.macropore_ice_m3;
        balance.total_energy_megajoules += energy;
        if (s.heat_capacity_megajoules_k > p.minimum_heat_capacity_megajoules_k and s.layer_thickness_m > p.minimum_layer_thickness_m) {
            s.temperature_k = (energy + f.conductive_heat_megajoules + f.matrix_latent_heat_megajoules + f.vapor_latent_heat_megajoules + f.root_convective_heat_megajoules + f.subsurface_convective_heat_megajoules + f.combustion_heat_megajoules) / s.heat_capacity_megajoules_k;
        } else {
            const replacement = if (layer == 0) top_boundary_temperature_k else states[layer - 1].temperature_k;
            heat_input += (s.temperature_k - replacement) * s.heat_capacity_megajoules_k;
            s.temperature_k = replacement;
        }
        heat_input += f.combustion_heat_megajoules;
        s.temperature_c = s.temperature_k - 273.15;
        s.maximum_temperature_c = @max(s.maximum_temperature_c, s.temperature_c);
        s.minimum_temperature_c = @min(s.minimum_temperature_c, s.temperature_c);
        if (!valid(s.*) or !valid(balance) or !std.math.isFinite(heat_input)) return error.NonFiniteSoilLayerUpdate;
    }
    return .{ .balance = balance, .heat_input_megajoules = heat_input };
}

test "REDIST soil layer update traverses runtime layers and fallback order" {
    var states = [_]LayerState{ std.mem.zeroes(LayerState), std.mem.zeroes(LayerState) };
    for (&states) |*s| {
        s.temperature_k = 280;
        s.heat_capacity_megajoules_k = 1;
        s.matrix_heat_capacity_megajoules_k = 1;
        s.layer_thickness_m = 0;
        s.maximum_temperature_c = -100;
        s.minimum_temperature_c = 100;
    }
    const fluxes = [_]LayerFlux{ std.mem.zeroes(LayerFlux), std.mem.zeroes(LayerFlux) };
    _ = try updateLayers(&states, &fluxes, balance_module.initialize(), 0, 275, .{ .negligible_water_m3 = 0, .ice_density_megagrams_m3 = 0.9, .substeps_per_hour_reciprocal = 1, .liquid_heat_capacity_megajoules_m3_k = 4.19, .ice_heat_capacity_megajoules_m3_k = 1.9274, .minimum_heat_capacity_megajoules_k = 0, .minimum_layer_thickness_m = 1 });
    try std.testing.expectEqual(@as(f64, 275), states[0].temperature_k);
    try std.testing.expectEqual(@as(f64, 275), states[1].temperature_k);
}

test "REDIST soil layer update preserves water and heat source order" {
    var state = std.mem.zeroes(LayerState);
    state.temperature_k = 280;
    state.heat_capacity_megajoules_k = 2;
    state.matrix_heat_capacity_megajoules_k = 1;
    state.matrix_water_m3 = 1;
    state.mobile_matrix_water_m3 = 0.5;
    state.layer_thickness_m = 1;
    state.bulk_density_megagrams_m3 = 1;
    state.matrix_air_capacity_m3 = 5;
    state.maximum_temperature_c = -100;
    state.minimum_temperature_c = 100;
    var flux = std.mem.zeroes(LayerFlux);
    flux.matrix_water_m3 = 1;
    flux.combustion_heat_megajoules = 2;
    var states = [_]LayerState{state};
    const fluxes = [_]LayerFlux{flux};
    const result = try updateLayers(&states, &fluxes, balance_module.initialize(), 0, 270, .{ .negligible_water_m3 = 0, .ice_density_megagrams_m3 = 1, .substeps_per_hour_reciprocal = 1, .liquid_heat_capacity_megajoules_m3_k = 4, .ice_heat_capacity_megajoules_m3_k = 2, .minimum_heat_capacity_megajoules_k = 0, .minimum_layer_thickness_m = 0 });
    try std.testing.expectEqual(@as(f64, 2), result.heat_input_megajoules);
}
