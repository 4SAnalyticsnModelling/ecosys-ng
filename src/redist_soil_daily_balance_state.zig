const std = @import("std");

pub const PhysicalState = struct {
    heat_capacity_megajoules_per_k: f64, // VHCP
    temperature_k: f64, // TKS
    sand_megagrams: f64, // SAND
    silt_megagrams: f64, // SILT
    clay_megagrams: f64, // CLAY
};

pub const GasState = struct {
    co2_gas_g_c: f64,
    co2_micropore_g_c: f64,
    co2_macropore_g_c: f64,
    root_co2_g_c: f64,
    ch4_gas_g_c: f64,
    ch4_micropore_g_c: f64,
    ch4_macropore_g_c: f64,
    root_ch4_g_c: f64,
    h2_gas_g_h: f64,
    h2_micropore_g_h: f64,
    h2_macropore_g_h: f64,
    root_h2_g_h: f64,
    o2_gas_g_o: f64,
    o2_micropore_g_o: f64,
    o2_macropore_g_o: f64,
    root_o2_g_o: f64,
    n2_gas_g_n: f64,
    n2_micropore_g_n: f64,
    n2_macropore_g_n: f64,
    root_n2o_g_n: f64, // TLN2OP
    n2o_gas_g_n: f64,
    n2o_micropore_g_n: f64,
    n2o_macropore_g_n: f64,
    root_nh3_g_n: f64,
    nh3_gas_g_n: f64,
};

pub const MineralNitrogenState = struct {
    nh4_micropore_g_n: f64,
    nh4_macropore_g_n: f64,
    nh4_band_micropore_g_n: f64,
    nh4_band_macropore_g_n: f64,
    nh3_micropore_g_n: f64,
    nh3_macropore_g_n: f64,
    nh3_band_micropore_g_n: f64,
    nh3_band_macropore_g_n: f64,
    nh4_exchange_nonband_mol: f64,
    nh4_exchange_band_mol: f64,
    nh4_fertilizer_a_mol: f64,
    urea_fertilizer_a_mol: f64,
    nh3_fertilizer_a_mol: f64,
    nh4_fertilizer_b_mol: f64,
    urea_fertilizer_b_mol: f64,
    nh3_fertilizer_b_mol: f64,
    no3_micropore_g_n: f64,
    no3_macropore_g_n: f64,
    no3_band_micropore_g_n: f64,
    no3_band_macropore_g_n: f64,
    no2_micropore_g_n: f64,
    no2_macropore_g_n: f64,
    no2_band_micropore_g_n: f64,
    no2_band_macropore_g_n: f64,
    no3_fertilizer_a_mol: f64,
    no3_fertilizer_b_mol: f64,
};

pub const PhosphorusState = struct {
    h2po4_micropore_g_p: f64,
    h2po4_macropore_g_p: f64,
    h2po4_band_micropore_g_p: f64,
    h2po4_band_macropore_g_p: f64,
    hpo4_micropore_g_p: f64,
    hpo4_macropore_g_p: f64,
    hpo4_band_micropore_g_p: f64,
    hpo4_band_macropore_g_p: f64,
    hpo4_exchange_nonband_mol: f64,
    h2po4_exchange_nonband_mol: f64,
    hpo4_exchange_band_mol: f64,
    h2po4_exchange_band_mol: f64,
    alpo4_nonband_mol: f64,
    fepo4_nonband_mol: f64,
    cahpo4_nonband_mol: f64,
    alpo4_band_mol: f64,
    fepo4_band_mol: f64,
    cahpo4_band_mol: f64,
    cah2po4_nonband_mol: f64,
    cah2po4_band_mol: f64,
    apatite_nonband_mol: f64,
    apatite_band_mol: f64,
};

pub const LayerState = struct {
    physical: PhysicalState,
    gases: GasState,
    nitrogen: MineralNitrogenState,
    phosphorus: PhosphorusState,
};

pub const DailyBalance = struct {
    landscape_heat_megajoules: f64,
    landscape_sediment_megagrams: f64,
    landscape_carbon_g_c: f64,
    grid_cell_carbon_g_c: f64,
    landscape_hydrogen_g_h: f64,
    landscape_oxygen_g_o: f64,
    landscape_nitrogen_g_n: f64,
    landscape_ammoniacal_n_g_n: f64,
    grid_cell_ammoniacal_n_g_n: f64,
    landscape_oxidized_n_g_n: f64,
    grid_cell_oxidized_n_g_n: f64,
    landscape_phosphate_g_p: f64,
    grid_cell_soluble_phosphate_g_p: f64,
    grid_cell_exchangeable_phosphate_g_p: f64,
    grid_cell_precipitated_phosphate_g_p: f64,
};

fn allFinite(value: anytype) bool {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        const field_value = @field(value, field.name);
        switch (@typeInfo(@TypeOf(field_value))) {
            .float => if (!std.math.isFinite(field_value)) return false,
            .@"struct" => if (!allFinite(field_value)) return false,
            else => @compileError("daily balance state supports only floats and nested structs"),
        }
    }
    return true;
}

/// Direct translation of REDIST 6685--6760 in ascending runtime-layer order.
pub fn accumulateLayers(
    balance: *DailyBalance,
    layers: []const LayerState,
    hour_of_day: u8,
    forcing_zone_index: usize,
    reference_forcing_zone_index: usize,
) !void {
    if (layers.len == 0) return error.SoilDailyBalanceDimensionMismatch;
    if (!allFinite(balance.*)) return error.InvalidSoilDailyBalance;
    for (layers) |layer| if (!allFinite(layer)) return error.InvalidSoilDailyLayer;
    if (hour_of_day != 24 or forcing_zone_index != reference_forcing_zone_index) return;

    var next = balance.*;
    for (layers) |layer| {
        const p = layer.physical;
        const g = layer.gases;
        const n = layer.nitrogen;
        const phosphorus = layer.phosphorus;
        next.landscape_heat_megajoules = next.landscape_heat_megajoules + p.heat_capacity_megajoules_per_k * p.temperature_k;
        const sediment_megagrams = p.sand_megagrams + p.silt_megagrams + p.clay_megagrams;
        next.landscape_sediment_megagrams = next.landscape_sediment_megagrams + sediment_megagrams;
        const carbon_g_c = g.co2_gas_g_c + g.co2_micropore_g_c + g.co2_macropore_g_c + g.root_co2_g_c +
            g.ch4_gas_g_c + g.ch4_micropore_g_c + g.ch4_macropore_g_c + g.root_ch4_g_c;
        next.landscape_carbon_g_c = next.landscape_carbon_g_c + carbon_g_c;
        next.grid_cell_carbon_g_c = next.grid_cell_carbon_g_c + carbon_g_c;
        const hydrogen_g_h = g.h2_gas_g_h + g.h2_micropore_g_h + g.h2_macropore_g_h + g.root_h2_g_h;
        next.landscape_hydrogen_g_h = next.landscape_hydrogen_g_h + hydrogen_g_h;
        const oxygen_g_o = g.o2_gas_g_o + g.o2_micropore_g_o + g.o2_macropore_g_o + g.root_o2_g_o;
        next.landscape_oxygen_g_o = next.landscape_oxygen_g_o + oxygen_g_o;
        const nitrogen_g_n = g.n2_gas_g_n + g.n2_micropore_g_n + g.n2_macropore_g_n + g.root_n2o_g_n +
            g.n2o_gas_g_n + g.n2o_micropore_g_n + g.n2o_macropore_g_n + g.root_nh3_g_n + g.nh3_gas_g_n;
        next.landscape_nitrogen_g_n = next.landscape_nitrogen_g_n + nitrogen_g_n;

        const ammoniacal_solute_g_n = n.nh4_micropore_g_n + n.nh4_macropore_g_n + n.nh4_band_micropore_g_n + n.nh4_band_macropore_g_n +
            n.nh3_micropore_g_n + n.nh3_macropore_g_n + n.nh3_band_micropore_g_n + n.nh3_band_macropore_g_n;
        const ammoniacal_exchange_g_n = 14.0 * (n.nh4_exchange_nonband_mol + n.nh4_exchange_band_mol);
        const ammoniacal_fertilizer_g_n = 14.0 * (n.nh4_fertilizer_a_mol + n.urea_fertilizer_a_mol + n.nh3_fertilizer_a_mol +
            n.nh4_fertilizer_b_mol + n.urea_fertilizer_b_mol + n.nh3_fertilizer_b_mol);
        next.landscape_ammoniacal_n_g_n = next.landscape_ammoniacal_n_g_n + ammoniacal_solute_g_n + ammoniacal_exchange_g_n + ammoniacal_fertilizer_g_n;
        next.grid_cell_ammoniacal_n_g_n = next.grid_cell_ammoniacal_n_g_n + ammoniacal_solute_g_n + ammoniacal_exchange_g_n;

        const oxidized_solute_g_n = n.no3_micropore_g_n + n.no3_macropore_g_n + n.no3_band_micropore_g_n + n.no3_band_macropore_g_n +
            n.no2_micropore_g_n + n.no2_macropore_g_n + n.no2_band_micropore_g_n + n.no2_band_macropore_g_n;
        const oxidized_fertilizer_g_n = 14.0 * (n.no3_fertilizer_a_mol + n.no3_fertilizer_b_mol);
        next.landscape_oxidized_n_g_n = next.landscape_oxidized_n_g_n + oxidized_solute_g_n + oxidized_fertilizer_g_n;
        next.grid_cell_oxidized_n_g_n = next.grid_cell_oxidized_n_g_n + oxidized_solute_g_n;

        const soluble_phosphate_g_p = phosphorus.h2po4_micropore_g_p + phosphorus.h2po4_macropore_g_p + phosphorus.h2po4_band_micropore_g_p + phosphorus.h2po4_band_macropore_g_p +
            phosphorus.hpo4_micropore_g_p + phosphorus.hpo4_macropore_g_p + phosphorus.hpo4_band_micropore_g_p + phosphorus.hpo4_band_macropore_g_p;
        const exchangeable_phosphate_g_p = 31.0 * (phosphorus.hpo4_exchange_nonband_mol + phosphorus.h2po4_exchange_nonband_mol +
            phosphorus.hpo4_exchange_band_mol + phosphorus.h2po4_exchange_band_mol);
        const precipitated_phosphate_g_p = 31.0 * (phosphorus.alpo4_nonband_mol + phosphorus.fepo4_nonband_mol + phosphorus.cahpo4_nonband_mol +
            phosphorus.alpo4_band_mol + phosphorus.fepo4_band_mol + phosphorus.cahpo4_band_mol) +
            62.0 * (phosphorus.cah2po4_nonband_mol + phosphorus.cah2po4_band_mol) +
            93.0 * (phosphorus.apatite_nonband_mol + phosphorus.apatite_band_mol);
        next.landscape_phosphate_g_p = next.landscape_phosphate_g_p + soluble_phosphate_g_p + exchangeable_phosphate_g_p + precipitated_phosphate_g_p;
        next.grid_cell_soluble_phosphate_g_p = next.grid_cell_soluble_phosphate_g_p + soluble_phosphate_g_p;
        next.grid_cell_exchangeable_phosphate_g_p = next.grid_cell_exchangeable_phosphate_g_p + exchangeable_phosphate_g_p;
        next.grid_cell_precipitated_phosphate_g_p = next.grid_cell_precipitated_phosphate_g_p + precipitated_phosphate_g_p;
        if (!allFinite(next)) return error.NonFiniteSoilDailyBalance;
    }
    balance.* = next;
}

fn unitLayer() LayerState {
    var layer: LayerState = undefined;
    inline for (@typeInfo(LayerState).@"struct".fields) |group_field| {
        const group = &@field(layer, group_field.name);
        inline for (@typeInfo(@TypeOf(group.*)).@"struct".fields) |field| @field(group, field.name) = 1.0;
    }
    return layer;
}

test "REDIST daily balance maps complete physical gas N and P census" {
    var balance = std.mem.zeroes(DailyBalance);
    const layers = [_]LayerState{unitLayer()};
    try accumulateLayers(&balance, &layers, 24, 1, 1);
    try std.testing.expectEqual(@as(f64, 1.0), balance.landscape_heat_megajoules);
    try std.testing.expectEqual(@as(f64, 3.0), balance.landscape_sediment_megagrams);
    try std.testing.expectEqual(@as(f64, 8.0), balance.landscape_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 4.0), balance.landscape_hydrogen_g_h);
    try std.testing.expectEqual(@as(f64, 4.0), balance.landscape_oxygen_g_o);
    try std.testing.expectEqual(@as(f64, 9.0), balance.landscape_nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 120.0), balance.landscape_ammoniacal_n_g_n);
    try std.testing.expectEqual(@as(f64, 36.0), balance.grid_cell_ammoniacal_n_g_n);
    try std.testing.expectEqual(@as(f64, 36.0), balance.landscape_oxidized_n_g_n);
    try std.testing.expectEqual(@as(f64, 8.0), balance.grid_cell_oxidized_n_g_n);
    try std.testing.expectEqual(@as(f64, 628.0), balance.landscape_phosphate_g_p);
    try std.testing.expectEqual(@as(f64, 124.0), balance.grid_cell_exchangeable_phosphate_g_p);
    try std.testing.expectEqual(@as(f64, 496.0), balance.grid_cell_precipitated_phosphate_g_p);
}

test "REDIST daily balance requires both hour and reference-zone gates" {
    var balance = std.mem.zeroes(DailyBalance);
    const layers = [_]LayerState{unitLayer()};
    try accumulateLayers(&balance, &layers, 23, 1, 1);
    try std.testing.expectEqual(std.mem.zeroes(DailyBalance), balance);
    try accumulateLayers(&balance, &layers, 24, 2, 1);
    try std.testing.expectEqual(std.mem.zeroes(DailyBalance), balance);
}

test "REDIST daily balance runtime layers remain independent and ordered" {
    var balance = std.mem.zeroes(DailyBalance);
    var layers = [_]LayerState{ unitLayer(), unitLayer() };
    layers[0].physical.heat_capacity_megajoules_per_k = 2.0;
    layers[0].physical.temperature_k = 3.0;
    layers[1].physical.heat_capacity_megajoules_per_k = 4.0;
    layers[1].physical.temperature_k = 5.0;
    try accumulateLayers(&balance, &layers, 24, 0, 0);
    try std.testing.expectEqual(@as(f64, 26.0), balance.landscape_heat_megajoules);
    try std.testing.expectEqual(@as(f64, 16.0), balance.landscape_carbon_g_c);
}

test "REDIST daily balance rejects dimensions invalid values and overflow" {
    var balance = std.mem.zeroes(DailyBalance);
    const no_layers: [0]LayerState = .{};
    try std.testing.expectError(error.SoilDailyBalanceDimensionMismatch, accumulateLayers(&balance, &no_layers, 24, 0, 0));
    var layers = [_]LayerState{unitLayer()};
    layers[0].gases.n2o_gas_g_n = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilDailyLayer, accumulateLayers(&balance, &layers, 24, 0, 0));
    layers[0] = unitLayer();
    balance.landscape_heat_megajoules = std.math.floatMax(f64);
    layers[0].physical.heat_capacity_megajoules_per_k = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSoilDailyBalance, accumulateLayers(&balance, &layers, 24, 0, 0));
}
