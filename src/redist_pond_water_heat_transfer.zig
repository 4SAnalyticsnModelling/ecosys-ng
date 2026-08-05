const std = @import("std");

pub const Inputs = struct {
    source_layer: usize, // L0
    destination_layer: usize, // L1
    redistribution_pass: usize, // NN
    boundary_flag: u8, // IFLGL(L,NN)
    fraction: f64, // FX
    bulk_density_megagrams_per_m3: []const f64, // BKDS
    macropore_water_volume_m3: []const f64, // VOLWH
    macropore_ice_volume_m3: []const f64, // VOLIH
    minimum_heat_capacity_megajoules_per_k: f64, // VHCPRX
    zero_tolerance: f64, // ZERO
};

pub const State = struct {
    disturbance_flag: *u8, // IFLGS
    water_volume_m3: []f64, // VOLW
    ice_volume_m3: []f64, // VOLI
    pore_volume_m3: []f64, // VOLP
    air_volume_m3: []f64, // VOLA
    total_pore_volume_m3: []f64, // VOLY
    previous_water_volume_m3: []f64, // VOLWX
    mineral_heat_capacity_megajoules_per_k: []f64, // VHCM
    total_heat_capacity_megajoules_per_k: []f64, // VHCP
    temperature_k: []f64, // TKS
    temperature_c: []f64, // TCS
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8523--8588, through pond water, ice, and heat.
/// Returns whether the pond-material branch was selected.
pub fn apply(inputs: Inputs, state: State) !bool {
    const len = state.water_volume_m3.len;
    if (len == 0 or inputs.source_layer >= len or inputs.destination_layer >= len or
        inputs.source_layer == inputs.destination_layer or inputs.bulk_density_megagrams_per_m3.len != len or
        inputs.macropore_water_volume_m3.len != len or inputs.macropore_ice_volume_m3.len != len or
        state.ice_volume_m3.len != len or state.pore_volume_m3.len != len or state.air_volume_m3.len != len or
        state.total_pore_volume_m3.len != len or state.previous_water_volume_m3.len != len or
        state.mineral_heat_capacity_megajoules_per_k.len != len or state.total_heat_capacity_megajoules_per_k.len != len or
        state.temperature_k.len != len or state.temperature_c.len != len)
        return error.PondWaterHeatTransferDimensionMismatch;
    inline for (.{ inputs.bulk_density_megagrams_per_m3, inputs.macropore_water_volume_m3, inputs.macropore_ice_volume_m3, state.water_volume_m3, state.ice_volume_m3, state.pore_volume_m3, state.air_volume_m3, state.total_pore_volume_m3, state.previous_water_volume_m3, state.mineral_heat_capacity_megajoules_per_k, state.total_heat_capacity_megajoules_per_k, state.temperature_k, state.temperature_c }) |values|
        if (!finiteSlice(values)) return error.InvalidPondWaterHeatTransferInput;
    inline for (.{ inputs.fraction, inputs.minimum_heat_capacity_megajoules_per_k, inputs.zero_tolerance }) |value|
        if (!std.math.isFinite(value)) return error.InvalidPondWaterHeatTransferInput;
    if (inputs.fraction < 0 or inputs.fraction > 1 or inputs.minimum_heat_capacity_megajoules_per_k < 0 or
        inputs.zero_tolerance < 0)
        return error.InvalidPondWaterHeatTransferInput;
    if (inputs.fraction <= 0) return false;

    const source = inputs.source_layer;
    const destination = inputs.destination_layer;
    const is_pond_transfer = (inputs.bulk_density_megagrams_per_m3[source] <= inputs.zero_tolerance and
        inputs.bulk_density_megagrams_per_m3[destination] <= inputs.zero_tolerance) or
        (inputs.redistribution_pass > 1 and inputs.boundary_flag == 1);
    if (!is_pond_transfer) {
        state.disturbance_flag.* = 1;
        return false;
    }

    const water_volume_m3 = state.water_volume_m3[destination] +
        inputs.fraction * state.water_volume_m3[source];
    const ice_volume_m3 = state.ice_volume_m3[destination] +
        inputs.fraction * state.ice_volume_m3[source];
    const pore_volume_m3 = state.pore_volume_m3[destination] +
        inputs.fraction * state.pore_volume_m3[source];
    const air_volume_m3 = state.air_volume_m3[destination] +
        inputs.fraction * state.air_volume_m3[source];
    const total_pore_volume_m3 = state.total_pore_volume_m3[destination] +
        inputs.fraction * state.total_pore_volume_m3[source];
    const destination_energy_megajoules = state.total_heat_capacity_megajoules_per_k[destination] *
        state.temperature_k[destination];
    const source_energy_megajoules = state.total_heat_capacity_megajoules_per_k[source] * state.temperature_k[source];
    const combined_energy_megajoules = destination_energy_megajoules + inputs.fraction * source_energy_megajoules;
    const mineral_heat_capacity_megajoules_per_k = state.mineral_heat_capacity_megajoules_per_k[destination] +
        inputs.fraction * state.mineral_heat_capacity_megajoules_per_k[source];
    const total_heat_capacity_megajoules_per_k = mineral_heat_capacity_megajoules_per_k +
        4.19 * (water_volume_m3 + inputs.macropore_water_volume_m3[destination]) +
        1.9274 * (ice_volume_m3 + inputs.macropore_ice_volume_m3[destination]);
    const temperature_k = if (total_heat_capacity_megajoules_per_k > inputs.minimum_heat_capacity_megajoules_per_k)
        combined_energy_megajoules / total_heat_capacity_megajoules_per_k
    else
        state.temperature_k[source];
    const temperature_c = temperature_k - 273.15;
    inline for (.{ water_volume_m3, ice_volume_m3, pore_volume_m3, air_volume_m3, total_pore_volume_m3, combined_energy_megajoules, mineral_heat_capacity_megajoules_per_k, total_heat_capacity_megajoules_per_k, temperature_k, temperature_c }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePondWaterHeatTransferResult;

    state.disturbance_flag.* = 1;
    state.water_volume_m3[destination] = water_volume_m3;
    state.ice_volume_m3[destination] = ice_volume_m3;
    state.pore_volume_m3[destination] = pore_volume_m3;
    state.air_volume_m3[destination] = air_volume_m3;
    state.total_pore_volume_m3[destination] = total_pore_volume_m3;
    state.previous_water_volume_m3[destination] = state.water_volume_m3[destination];
    state.mineral_heat_capacity_megajoules_per_k[destination] = mineral_heat_capacity_megajoules_per_k;
    state.total_heat_capacity_megajoules_per_k[destination] = total_heat_capacity_megajoules_per_k;
    state.temperature_k[destination] = temperature_k;
    state.temperature_c[destination] = temperature_c;
    return true;
}

const Fixture = struct {
    density: [2]f64 = .{ 0, 0 },
    macro_water: [2]f64 = .{ 0.1, 0.2 },
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
    disturbance: u8 = 0,

    fn inputs(self: *Fixture, fraction: f64) Inputs {
        return .{
            .source_layer = 0,
            .destination_layer = 1,
            .redistribution_pass = 1,
            .boundary_flag = 0,
            .fraction = fraction,
            .bulk_density_megagrams_per_m3 = &self.density,
            .macropore_water_volume_m3 = &self.macro_water,
            .macropore_ice_volume_m3 = &self.macro_ice,
            .minimum_heat_capacity_megajoules_per_k = 0.01,
            .zero_tolerance = 1.0e-12,
        };
    }

    fn state(self: *Fixture) State {
        return .{
            .disturbance_flag = &self.disturbance,
            .water_volume_m3 = &self.water,
            .ice_volume_m3 = &self.ice,
            .pore_volume_m3 = &self.pore,
            .air_volume_m3 = &self.air,
            .total_pore_volume_m3 = &self.total_pore,
            .previous_water_volume_m3 = &self.previous_water,
            .mineral_heat_capacity_megajoules_per_k = &self.mineral_heat,
            .total_heat_capacity_megajoules_per_k = &self.total_heat,
            .temperature_k = &self.temperature_k,
            .temperature_c = &self.temperature_c,
        };
    }
};

test "REDIST pond water ice pore volume and heat move in source order" {
    var fixture = Fixture{};
    try std.testing.expect(try apply(fixture.inputs(0.25), fixture.state()));
    try std.testing.expectEqual(@as(u8, 1), fixture.disturbance);
    try std.testing.expectEqual(@as(f64, 1.5), fixture.water[1]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), fixture.ice[1], 1.0e-14);
    try std.testing.expectEqual(@as(f64, 2.75), fixture.pore[1]);
    try std.testing.expectEqual(@as(f64, 0.65), fixture.air[1]);
    try std.testing.expectEqual(@as(f64, 3.4), fixture.total_pore[1]);
    try std.testing.expectEqual(fixture.water[1], fixture.previous_water[1]);
    try std.testing.expectEqual(@as(f64, 4), fixture.mineral_heat[1]);
    const expected_heat = 4 + 4.19 * 1.7 + 1.9274 * 0.33;
    const expected_temperature = (9 * 270 + 0.25 * 14 * 280) / expected_heat;
    try std.testing.expectApproxEqAbs(expected_heat, fixture.total_heat[1], 1.0e-12);
    try std.testing.expectApproxEqAbs(expected_temperature, fixture.temperature_k[1], 1.0e-12);
    try std.testing.expectApproxEqAbs(expected_temperature - 273.15, fixture.temperature_c[1], 1.0e-12);
}

test "REDIST boundary pass selects pond transfer across density classes" {
    var fixture = Fixture{};
    fixture.density[1] = 1.2;
    var inputs = fixture.inputs(0.5);
    inputs.redistribution_pass = 3;
    inputs.boundary_flag = 1;
    try std.testing.expect(try apply(inputs, fixture.state()));
    try std.testing.expectEqual(@as(f64, 2), fixture.water[1]);
}

test "REDIST soil branch only marks disturbance in this pond kernel" {
    var fixture = Fixture{};
    fixture.density = .{ 1.1, 1.2 };
    try std.testing.expect(!try apply(fixture.inputs(0.25), fixture.state()));
    try std.testing.expectEqual(@as(u8, 1), fixture.disturbance);
    try std.testing.expectEqual(@as(f64, 1), fixture.water[1]);
}

test "REDIST zero fraction leaves state unchanged" {
    var fixture = Fixture{};
    try std.testing.expect(!try apply(fixture.inputs(0), fixture.state()));
    try std.testing.expectEqual(@as(u8, 0), fixture.disturbance);
    try std.testing.expectEqual(@as(f64, 1), fixture.water[1]);
}

test "REDIST pond transfer validates atomically before mutation" {
    var fixture = Fixture{};
    fixture.macro_water[1] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondWaterHeatTransferInput, apply(fixture.inputs(0.25), fixture.state()));
    try std.testing.expectEqual(@as(u8, 0), fixture.disturbance);
    try std.testing.expectEqual(@as(f64, 1), fixture.water[1]);
}
