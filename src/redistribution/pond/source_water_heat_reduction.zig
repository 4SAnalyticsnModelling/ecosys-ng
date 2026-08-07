const std = @import("std");

pub const Inputs = struct {
    source_layer: usize, // L0
    destination_layer: usize, // L1
    remaining_fraction: f64, // FY
    macropore_water_volume_m3: []const f64, // VOLWH
    macropore_ice_volume_m3: []const f64, // VOLIH
    minimum_heat_capacity_megajoules_per_k: f64, // VHCPRX
};

pub const State = struct {
    water_volume_m3: []f64, // VOLW
    ice_volume_m3: []f64, // VOLI
    pore_volume_m3: []f64, // VOLP
    air_volume_m3: []f64, // VOLA
    total_pore_volume_m3: []f64, // VOLY
    previous_water_volume_m3: []f64, // VOLWX
    source_energy_megajoules: *f64, // ENGY0
    mineral_heat_capacity_megajoules_per_k: []f64, // VHCM
    total_heat_capacity_megajoules_per_k: []f64, // VHCP
    temperature_k: []f64, // TKS
    temperature_c: []f64, // TCS
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 9086--9107: reduce pond source water and heat.
pub fn apply(inputs: Inputs, state: State) !void {
    const len = state.water_volume_m3.len;
    if (len == 0 or inputs.source_layer >= len or inputs.destination_layer >= len or
        inputs.source_layer == inputs.destination_layer or inputs.macropore_water_volume_m3.len != len or
        inputs.macropore_ice_volume_m3.len != len or state.ice_volume_m3.len != len or
        state.pore_volume_m3.len != len or state.air_volume_m3.len != len or
        state.total_pore_volume_m3.len != len or state.previous_water_volume_m3.len != len or
        state.mineral_heat_capacity_megajoules_per_k.len != len or state.total_heat_capacity_megajoules_per_k.len != len or
        state.temperature_k.len != len or state.temperature_c.len != len)
        return error.PondSourceWaterHeatReductionDimensionMismatch;
    inline for (.{ inputs.macropore_water_volume_m3, inputs.macropore_ice_volume_m3, state.water_volume_m3, state.ice_volume_m3, state.pore_volume_m3, state.air_volume_m3, state.total_pore_volume_m3, state.previous_water_volume_m3, state.mineral_heat_capacity_megajoules_per_k, state.total_heat_capacity_megajoules_per_k, state.temperature_k, state.temperature_c }) |values|
        if (!finiteSlice(values)) return error.InvalidPondSourceWaterHeatReductionInput;
    inline for (.{ inputs.remaining_fraction, inputs.minimum_heat_capacity_megajoules_per_k, state.source_energy_megajoules.* }) |value|
        if (!std.math.isFinite(value)) return error.InvalidPondSourceWaterHeatReductionInput;
    if (inputs.remaining_fraction < 0 or inputs.remaining_fraction > 1 or inputs.minimum_heat_capacity_megajoules_per_k < 0)
        return error.InvalidPondSourceWaterHeatReductionInput;

    const source = inputs.source_layer;
    const water_volume_m3 = inputs.remaining_fraction * state.water_volume_m3[source];
    const ice_volume_m3 = inputs.remaining_fraction * state.ice_volume_m3[source];
    const pore_volume_m3 = inputs.remaining_fraction * state.pore_volume_m3[source];
    const air_volume_m3 = inputs.remaining_fraction * state.air_volume_m3[source];
    const total_pore_volume_m3 = inputs.remaining_fraction * state.total_pore_volume_m3[source];
    const source_energy_megajoules = inputs.remaining_fraction * state.source_energy_megajoules.*;
    const mineral_heat_capacity_megajoules_per_k = inputs.remaining_fraction * state.mineral_heat_capacity_megajoules_per_k[source];
    const total_heat_capacity_megajoules_per_k = if (source != 0)
        mineral_heat_capacity_megajoules_per_k +
            4.19 * (water_volume_m3 + inputs.macropore_water_volume_m3[source]) +
            1.9274 * (ice_volume_m3 + inputs.macropore_ice_volume_m3[source])
    else
        mineral_heat_capacity_megajoules_per_k + 4.19 * water_volume_m3 + 1.9274 * ice_volume_m3;
    const temperature_k = if (total_heat_capacity_megajoules_per_k > inputs.minimum_heat_capacity_megajoules_per_k)
        source_energy_megajoules / total_heat_capacity_megajoules_per_k
    else
        state.temperature_k[inputs.destination_layer];
    const temperature_c = temperature_k - 273.15;
    inline for (.{ water_volume_m3, ice_volume_m3, pore_volume_m3, air_volume_m3, total_pore_volume_m3, source_energy_megajoules, mineral_heat_capacity_megajoules_per_k, total_heat_capacity_megajoules_per_k, temperature_k, temperature_c }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePondSourceWaterHeatReductionResult;

    state.water_volume_m3[source] = water_volume_m3;
    state.ice_volume_m3[source] = ice_volume_m3;
    state.pore_volume_m3[source] = pore_volume_m3;
    state.air_volume_m3[source] = air_volume_m3;
    state.total_pore_volume_m3[source] = total_pore_volume_m3;
    state.previous_water_volume_m3[source] = state.water_volume_m3[source];
    state.source_energy_megajoules.* = source_energy_megajoules;
    state.mineral_heat_capacity_megajoules_per_k[source] = mineral_heat_capacity_megajoules_per_k;
    state.total_heat_capacity_megajoules_per_k[source] = total_heat_capacity_megajoules_per_k;
    state.temperature_k[source] = temperature_k;
    state.temperature_c[source] = temperature_c;
}

const Fixture = struct {
    macro_water: [2]f64 = .{ 0.2, 0.3 },
    macro_ice: [2]f64 = .{ 0.02, 0.03 },
    water: [2]f64 = .{ 2, 1 },
    ice: [2]f64 = .{ 0.4, 0.2 },
    pore: [2]f64 = .{ 3, 2 },
    air: [2]f64 = .{ 0.6, 0.5 },
    total_pore: [2]f64 = .{ 3.6, 2.5 },
    previous_water: [2]f64 = .{ 1.8, 0.8 },
    mineral_heat: [2]f64 = .{ 4, 3 },
    total_heat: [2]f64 = .{ 14, 9 },
    temperature_k: [2]f64 = .{ 280, 270 },
    temperature_c: [2]f64 = .{ 6.85, -3.15 },
    energy: f64 = 3920,
    fn state(self: *Fixture) State {
        return .{ .water_volume_m3 = &self.water, .ice_volume_m3 = &self.ice, .pore_volume_m3 = &self.pore, .air_volume_m3 = &self.air, .total_pore_volume_m3 = &self.total_pore, .previous_water_volume_m3 = &self.previous_water, .source_energy_megajoules = &self.energy, .mineral_heat_capacity_megajoules_per_k = &self.mineral_heat, .total_heat_capacity_megajoules_per_k = &self.total_heat, .temperature_k = &self.temperature_k, .temperature_c = &self.temperature_c };
    }
};

test "REDIST subsurface source reduction includes macropore heat capacity" {
    var fixture = Fixture{};
    const inputs = Inputs{ .source_layer = 1, .destination_layer = 0, .remaining_fraction = 0.25, .macropore_water_volume_m3 = &fixture.macro_water, .macropore_ice_volume_m3 = &fixture.macro_ice, .minimum_heat_capacity_megajoules_per_k = 0.01 };
    fixture.energy = 9 * 270;
    try apply(inputs, fixture.state());
    try std.testing.expectEqual(@as(f64, 0.25), fixture.water[1]);
    try std.testing.expectEqual(@as(f64, 0.05), fixture.ice[1]);
    try std.testing.expectEqual(@as(f64, 0.5), fixture.pore[1]);
    try std.testing.expectEqual(@as(f64, 0.125), fixture.air[1]);
    try std.testing.expectEqual(@as(f64, 0.625), fixture.total_pore[1]);
    try std.testing.expectEqual(@as(f64, 0.75), fixture.mineral_heat[1]);
    const heat = 0.75 + 4.19 * 0.55 + 1.9274 * 0.08;
    try std.testing.expectApproxEqAbs(heat, fixture.total_heat[1], 1.0e-12);
    try std.testing.expectApproxEqAbs((0.25 * 9 * 270) / heat, fixture.temperature_k[1], 1.0e-12);
}

test "REDIST layer zero source excludes macropore heat capacity" {
    var fixture = Fixture{};
    const inputs = Inputs{ .source_layer = 0, .destination_layer = 1, .remaining_fraction = 0.5, .macropore_water_volume_m3 = &fixture.macro_water, .macropore_ice_volume_m3 = &fixture.macro_ice, .minimum_heat_capacity_megajoules_per_k = 0.01 };
    try apply(inputs, fixture.state());
    const heat = 2 + 4.19 * 1 + 1.9274 * 0.2;
    try std.testing.expectApproxEqAbs(heat, fixture.total_heat[0], 1.0e-12);
}

test "REDIST exhausted source inherits destination temperature" {
    var fixture = Fixture{};
    const inputs = Inputs{ .source_layer = 0, .destination_layer = 1, .remaining_fraction = 0, .macropore_water_volume_m3 = &fixture.macro_water, .macropore_ice_volume_m3 = &fixture.macro_ice, .minimum_heat_capacity_megajoules_per_k = 0.01 };
    try apply(inputs, fixture.state());
    try std.testing.expectEqual(@as(f64, 270), fixture.temperature_k[0]);
    try std.testing.expectApproxEqAbs(@as(f64, -3.15), fixture.temperature_c[0], 1.0e-12);
}

test "REDIST source heat validation is atomic" {
    var fixture = Fixture{};
    fixture.macro_ice[1] = std.math.inf(f64);
    const inputs = Inputs{ .source_layer = 1, .destination_layer = 0, .remaining_fraction = 0.5, .macropore_water_volume_m3 = &fixture.macro_water, .macropore_ice_volume_m3 = &fixture.macro_ice, .minimum_heat_capacity_megajoules_per_k = 0.01 };
    try std.testing.expectError(error.InvalidPondSourceWaterHeatReductionInput, apply(inputs, fixture.state()));
    try std.testing.expectEqual(@as(f64, 1), fixture.water[1]);
}
