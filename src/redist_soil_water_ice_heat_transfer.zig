const std = @import("std");

pub const Context = struct {
    source_layer: usize,
    source_bulk_density_megagrams_m3: f64,
    destination_bulk_density_megagrams_m3: f64,
    redistribution_fraction: f64, // FWO=FX
    pond_water_volume_before_redistribution_m3: f64, // XVOLWP
    retained_surface_water_volume_m3: f64, // VOLWD
    minimum_heat_capacity_megajoules_per_k: f64, // VHCPRX
    current_layer_temperature_k: f64, // TKS(L)
};

pub const Layer = struct {
    water_volume_m3: f64, // VOLW
    vapor_volume_m3: f64, // VOLV
    ice_volume_m3: f64, // VOLI
    total_water_equivalent_volume_m3: f64, // VOLY
    auxiliary_water_volume_m3: f64, // VOLWX
    mineral_heat_capacity_megajoules_per_k: f64, // VHCM
    total_heat_capacity_megajoules_per_k: f64, // VHCP
    temperature_k: f64, // TKS
    temperature_c: f64, // TCS
};

pub const State = struct { source: Layer, destination: Layer };

fn finiteStruct(value: anytype) bool {
    inline for (std.meta.fields(@TypeOf(value))) |field|
        if (field.type == f64 and !std.math.isFinite(@field(value, field.name))) return false;
    return true;
}

fn move(source: f64, destination: f64, amount: f64) !struct { source: f64, destination: f64 } {
    const next_destination = destination + amount;
    const next_source = source - amount;
    if (!std.math.isFinite(next_destination) or !std.math.isFinite(next_source))
        return error.NonFiniteSoilWaterIceHeatTransferResult;
    return .{ .source = next_source, .destination = next_destination };
}

/// Direct translation of REDIST 9625--9666 under the positive-`BKDS` soil gate.
pub fn transfer(context: Context, state: *State) !void {
    if (!finiteStruct(context) or !finiteStruct(state.source) or !finiteStruct(state.destination) or
        context.redistribution_fraction < 0 or context.redistribution_fraction > 1)
        return error.InvalidSoilWaterIceHeatTransferInput;
    if (context.source_bulk_density_megagrams_m3 <= 0 or context.destination_bulk_density_megagrams_m3 <= 0) return;

    const water_fraction = context.redistribution_fraction; // FWO=FX.
    const moved_water_m3 = if (context.source_layer == 0)
        water_fraction * @max(0.0, context.pond_water_volume_before_redistribution_m3 - context.retained_surface_water_volume_m3)
    else
        water_fraction * state.source.water_volume_m3;
    const moved_vapor_m3 = if (context.source_layer == 0) 0.0 else water_fraction * state.source.vapor_volume_m3;
    const moved_ice_m3 = water_fraction * state.source.ice_volume_m3;
    const moved_total_water_m3 = water_fraction * state.source.total_water_equivalent_volume_m3;
    const moved_auxiliary_water_m3 = water_fraction * state.source.auxiliary_water_volume_m3;
    const moved_mineral_heat_capacity_megajoules_per_k = water_fraction * state.source.mineral_heat_capacity_megajoules_per_k;
    inline for (.{ moved_water_m3, moved_vapor_m3, moved_ice_m3, moved_total_water_m3, moved_auxiliary_water_m3, moved_mineral_heat_capacity_megajoules_per_k }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilWaterIceHeatTransferResult;

    var next = state.*;
    inline for (.{
        .{ "water_volume_m3", moved_water_m3 },
        .{ "vapor_volume_m3", moved_vapor_m3 },
        .{ "ice_volume_m3", moved_ice_m3 },
        .{ "total_water_equivalent_volume_m3", moved_total_water_m3 },
        .{ "auxiliary_water_volume_m3", moved_auxiliary_water_m3 },
        .{ "mineral_heat_capacity_megajoules_per_k", moved_mineral_heat_capacity_megajoules_per_k },
    }) |field_and_amount| {
        const result = try move(@field(state.source, field_and_amount[0]), @field(state.destination, field_and_amount[0]), field_and_amount[1]);
        @field(next.source, field_and_amount[0]) = result.source;
        @field(next.destination, field_and_amount[0]) = result.destination;
    }
    const moved_heat_capacity_megajoules_per_k = moved_mineral_heat_capacity_megajoules_per_k +
        4.19 * (moved_water_m3 + moved_vapor_m3) + 1.9274 * moved_ice_m3;
    const moved_energy_megajoules = state.source.temperature_k * moved_heat_capacity_megajoules_per_k;
    const destination_energy_megajoules = state.destination.total_heat_capacity_megajoules_per_k * state.destination.temperature_k + moved_energy_megajoules;
    const source_energy_megajoules = state.source.total_heat_capacity_megajoules_per_k * state.source.temperature_k - moved_energy_megajoules;
    next.destination.total_heat_capacity_megajoules_per_k = state.destination.total_heat_capacity_megajoules_per_k + moved_heat_capacity_megajoules_per_k;
    next.source.total_heat_capacity_megajoules_per_k = state.source.total_heat_capacity_megajoules_per_k - moved_heat_capacity_megajoules_per_k;
    next.destination.temperature_k = if (next.destination.total_heat_capacity_megajoules_per_k > context.minimum_heat_capacity_megajoules_per_k)
        destination_energy_megajoules / next.destination.total_heat_capacity_megajoules_per_k
    else
        context.current_layer_temperature_k;
    next.destination.temperature_c = next.destination.temperature_k - 273.15;
    next.source.temperature_k = if (next.source.total_heat_capacity_megajoules_per_k > context.minimum_heat_capacity_megajoules_per_k)
        source_energy_megajoules / next.source.total_heat_capacity_megajoules_per_k
    else
        context.current_layer_temperature_k;
    next.source.temperature_c = next.source.temperature_k - 273.15;
    if (!finiteStruct(next.source) or !finiteStruct(next.destination)) return error.NonFiniteSoilWaterIceHeatTransferResult;
    state.* = next;
}

fn fixture() State {
    return .{
        .source = .{ .water_volume_m3 = 2, .vapor_volume_m3 = 1, .ice_volume_m3 = 0, .total_water_equivalent_volume_m3 = 3, .auxiliary_water_volume_m3 = 1, .mineral_heat_capacity_megajoules_per_k = 2, .total_heat_capacity_megajoules_per_k = 20, .temperature_k = 280, .temperature_c = 6.85 },
        .destination = .{ .water_volume_m3 = 10, .vapor_volume_m3 = 2, .ice_volume_m3 = 0, .total_water_equivalent_volume_m3 = 12, .auxiliary_water_volume_m3 = 2, .mineral_heat_capacity_megajoules_per_k = 4, .total_heat_capacity_megajoules_per_k = 40, .temperature_k = 300, .temperature_c = 26.85 },
    };
}

test "REDIST subsurface water vapor and heat transfer preserves order and energy" {
    var state = fixture();
    const initial_energy = state.source.total_heat_capacity_megajoules_per_k * state.source.temperature_k + state.destination.total_heat_capacity_megajoules_per_k * state.destination.temperature_k;
    try transfer(.{ .source_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .redistribution_fraction = 0.25, .pond_water_volume_before_redistribution_m3 = 0, .retained_surface_water_volume_m3 = 0, .minimum_heat_capacity_megajoules_per_k = 0.01, .current_layer_temperature_k = 290 }, &state);
    try std.testing.expectEqual(@as(f64, 1.5), state.source.water_volume_m3);
    try std.testing.expectEqual(@as(f64, 0.75), state.source.vapor_volume_m3);
    try std.testing.expectEqual(@as(f64, 10.5), state.destination.water_volume_m3);
    const final_energy = state.source.total_heat_capacity_megajoules_per_k * state.source.temperature_k + state.destination.total_heat_capacity_megajoules_per_k * state.destination.temperature_k;
    try std.testing.expectApproxEqAbs(initial_energy, final_energy, 1e-10);
}

test "REDIST surface source uses pond reserve and moves no vapor" {
    var state = fixture();
    try transfer(.{ .source_layer = 0, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .redistribution_fraction = 0.5, .pond_water_volume_before_redistribution_m3 = 1.5, .retained_surface_water_volume_m3 = 0.5, .minimum_heat_capacity_megajoules_per_k = 0.01, .current_layer_temperature_k = 290 }, &state);
    try std.testing.expectEqual(@as(f64, 1.5), state.source.water_volume_m3);
    try std.testing.expectEqual(@as(f64, 10.5), state.destination.water_volume_m3);
    try std.testing.expectEqual(@as(f64, 1), state.source.vapor_volume_m3);
    try std.testing.expectEqual(@as(f64, 2), state.destination.vapor_volume_m3);
}

test "REDIST heat-capacity fallback and invalid input are atomic" {
    var state = fixture();
    try transfer(.{ .source_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .redistribution_fraction = 0, .pond_water_volume_before_redistribution_m3 = 0, .retained_surface_water_volume_m3 = 0, .minimum_heat_capacity_megajoules_per_k = 100, .current_layer_temperature_k = 290 }, &state);
    try std.testing.expectEqual(@as(f64, 290), state.source.temperature_k);
    try std.testing.expectEqual(@as(f64, 290), state.destination.temperature_k);
    const before = state;
    try std.testing.expectError(error.InvalidSoilWaterIceHeatTransferInput, transfer(.{ .source_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .redistribution_fraction = std.math.nan(f64), .pond_water_volume_before_redistribution_m3 = 0, .retained_surface_water_volume_m3 = 0, .minimum_heat_capacity_megajoules_per_k = 0, .current_layer_temperature_k = 290 }, &state));
    try std.testing.expectEqualDeep(before, state);
}
