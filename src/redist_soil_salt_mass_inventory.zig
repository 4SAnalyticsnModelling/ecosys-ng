const std = @import("std");
const micro = @import("redist_soil_salt_micropore_update.zig");
const macro = @import("redist_soil_salt_macropore_update.zig");
const exchange_update = @import("redist_soil_exchange_site_update.zig");
const phosphate_update = @import("redist_soil_phosphate_precipitate_update.zig");
const mineral_update = @import("redist_soil_precipitate_silicate_update.zig");

pub const FertilizerPools = struct {
    nh3_a_mol: f64,
    urea_a_mol: f64,
    no3_a_mol: f64,
    nh4_a_mol: f64,
    nh3_b_mol: f64,
    urea_b_mol: f64,
    no3_b_mol: f64,
    nh4_b_mol: f64,
};

pub const LayerInventory = struct {
    micropore: micro.SaltMicroporePools,
    macropore: macro.SaltMacroporePools,
    fertilizer: FertilizerPools,
    exchange: exchange_update.ExchangeSitePools,
    carbonate_exchange: macro.CarbonateExchangeSite,
    phosphate_precipitates: phosphate_update.PhosphatePrecipitatePools,
    hydroxide_precipitates: mineral_update.HydroxidePrecipitatePools,
    silicates: mineral_update.SilicatePools,
};

pub const MassInventory = struct {
    landscape_phosphate_g_p: f64, // TLPO4
    landscape_ion_mol: f64, // TION
    grid_cell_ion_mol: f64, // UION
};

fn finiteStruct(value: anytype) bool {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name))) return false;
    return true;
}

fn validateLayer(layer: LayerInventory) bool {
    return finiteStruct(layer.micropore) and finiteStruct(layer.macropore) and
        finiteStruct(layer.fertilizer) and finiteStruct(layer.exchange) and
        finiteStruct(layer.carbonate_exchange) and finiteStruct(layer.phosphate_precipitates) and
        finiteStruct(layer.hydroxide_precipitates) and finiteStruct(layer.silicates);
}

/// Direct translation of REDIST 7211--7287. The ion total intentionally uses
/// the legacy atom-count convention, not electrical charge equivalents.
pub fn accumulateLayers(
    inventory: *MassInventory,
    layers: []const LayerInventory,
    hour_of_day: u8,
    forcing_zone_index: usize,
    reference_forcing_zone_index: usize,
) !void {
    if (layers.len == 0) return error.SoilSaltInventoryDimensionMismatch;
    if (!finiteStruct(inventory.*)) return error.InvalidSoilSaltInventory;
    for (layers) |layer| if (!validateLayer(layer)) return error.InvalidSoilSaltInventoryLayer;
    if (hour_of_day != 24 or forcing_zone_index != reference_forcing_zone_index) return;

    var next = inventory.*;
    for (layers) |layer| {
        const s = layer.micropore;
        const h = layer.macropore;
        const f = layer.fertilizer;
        const x = layer.exchange;
        const p = layer.phosphate_precipitates;
        const q = layer.silicates;
        const z = layer.hydroxide_precipitates;
        const phosphate_g_p = 31.0 * (s.h0po4_mol + s.h3po4_mol + s.fe1p_mol + s.fe2p_mol + s.ca0p_mol + s.ca1p_mol + s.ca2p_mol + s.mg1p_mol +
            s.h0pob_mol + s.h3pob_mol + s.fe1pb_mol + s.fe2pb_mol + s.ca0pb_mol + s.ca1pb_mol + s.ca2pb_mol + s.mg1pb_mol +
            h.h0po4_mol + h.h3po4_mol + h.fe1p_mol + h.fe2p_mol + h.ca0p_mol + h.ca1p_mol + h.ca2p_mol + h.mg1p_mol +
            h.h0pob_mol + h.h3pob_mol + h.fe1pb_mol + h.fe2pb_mol + h.ca0pb_mol + h.ca1pb_mol + h.ca2pb_mol + h.mg1pb_mol);
        next.landscape_phosphate_g_p = next.landscape_phosphate_g_p + phosphate_g_p;

        const micropore_ion_mol = s.al_mol + s.fe_mol + s.h_mol + s.ca_mol + s.megagrams_mol + s.na_mol + s.ka_mol + s.oh_mol +
            s.so4_mol + s.cl_mol + s.co3_mol + s.h0po4_mol + s.h0pob_mol +
            2.0 * (s.hco3_mol + s.aloh1_mol + s.als_mol + s.feoh1_mol + s.fes_mol + s.cao_mol + s.cac_mol + s.cas_mol + s.mgo_mol + s.mgc_mol + s.mgs_mol + s.nac_mol + s.nas_mol + s.kas_mol + s.ca0p_mol + s.ca0pb_mol) +
            3.0 * (s.aloh2_mol + s.feoh2_mol + s.cah_mol + s.mgh_mol + s.fe1p_mol + s.ca1p_mol + s.mg1p_mol + s.fe1pb_mol + s.ca1pb_mol + s.mg1pb_mol) +
            4.0 * (s.aloh3_mol + s.feoh3_mol + s.h3po4_mol + s.fe2p_mol + s.ca2p_mol + s.h3pob_mol + s.fe2pb_mol + s.ca2pb_mol + s.hysi_mol) +
            5.0 * (s.aloh4_mol + s.feoh4_mol);
        const macropore_ion_mol = h.al_mol + h.fe_mol + h.h_mol + h.ca_mol + h.megagrams_mol + h.na_mol + h.ka_mol + h.oh_mol +
            h.so4_mol + h.cl_mol + h.co3_mol + h.h0po4_mol + h.h0pob_mol +
            2.0 * (h.hco3_mol + h.aloh1_mol + h.als_mol + h.feoh1_mol + h.fes_mol + h.cao_mol + h.cac_mol + h.cas_mol + h.mgo_mol + h.mgc_mol + h.mgs_mol + h.nac_mol + h.nas_mol + h.kas_mol + h.ca0p_mol + h.ca0pb_mol) +
            3.0 * (h.aloh2_mol + h.feoh2_mol + h.cah_mol + h.mgh_mol + h.fe1p_mol + h.ca1p_mol + h.mg1p_mol + h.fe1pb_mol + h.ca1pb_mol + h.mg1pb_mol) +
            4.0 * (h.aloh3_mol + h.feoh3_mol + h.h3po4_mol + h.fe2p_mol + h.ca2p_mol + h.h3pob_mol + h.fe2pb_mol + h.ca2pb_mol) +
            5.0 * (h.aloh4_mol + h.feoh4_mol);
        const fertilizer_ion_mol = f.nh3_a_mol + f.urea_a_mol + f.no3_a_mol + f.nh3_b_mol + f.urea_b_mol + f.no3_b_mol +
            2.0 * (f.nh4_a_mol + f.nh4_b_mol);
        const exchange_ion_mol = x.h_mol + x.al_mol + x.fe_mol + x.ca_mol + x.megagrams_mol + x.na_mol + x.ka_mol + layer.carbonate_exchange.hco3_mol +
            x.oh0_mol + x.oh0_band_mol +
            2.0 * (x.nh4_nonband_mol + x.nh4_band_mol + x.oh1_mol + x.oh1_band_mol) +
            3.0 * (x.oh2_mol + x.oh2_band_mol + x.hpo4_mol + x.hpo4_band_mol) +
            4.0 * (x.h2po4_mol + x.h2po4_band_mol);
        const precipitate_ion_mol = 2.0 * (z.caco_mol + z.caso_mol + p.alpo4_nonband_mol + p.fepo4_nonband_mol + p.alpo4_band_mol + p.fepo4_band_mol) +
            3.0 * (p.cahpo4_nonband_mol + p.cahpo4_band_mol) +
            4.0 * (z.aloh_mol + z.feoh_mol) +
            7.0 * (p.cah2po4_nonband_mol + p.cah2po4_band_mol) +
            9.0 * (p.apatite_nonband_mol + p.apatite_band_mol) +
            q.al_mol + q.fe_mol + q.ca_mol + q.megagrams_mol + q.na_mol + q.ka_mol +
            q.al2_mol + q.fe2_mol + q.ca2_mol + q.mg2_mol + q.na2_mol + q.ka2_mol;
        const layer_ion_mol = micropore_ion_mol + macropore_ion_mol + fertilizer_ion_mol + exchange_ion_mol + precipitate_ion_mol;
        next.landscape_ion_mol = next.landscape_ion_mol + layer_ion_mol;
        next.grid_cell_ion_mol = next.grid_cell_ion_mol + layer_ion_mol;
        if (!finiteStruct(next)) return error.NonFiniteSoilSaltInventory;
    }
    inventory.* = next;
}

fn unitStruct(comptime T: type) T {
    var value: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field| @field(value, field.name) = 1.0;
    return value;
}

fn unitLayer() LayerInventory {
    return .{
        .micropore = unitStruct(micro.SaltMicroporePools),
        .macropore = unitStruct(macro.SaltMacroporePools),
        .fertilizer = unitStruct(FertilizerPools),
        .exchange = unitStruct(exchange_update.ExchangeSitePools),
        .carbonate_exchange = unitStruct(macro.CarbonateExchangeSite),
        .phosphate_precipitates = unitStruct(phosphate_update.PhosphatePrecipitatePools),
        .hydroxide_precipitates = unitStruct(mineral_update.HydroxidePrecipitatePools),
        .silicates = unitStruct(mineral_update.SilicatePools),
    };
}

test "REDIST soil salt inventory preserves complete atom-count weighting" {
    var inventory = std.mem.zeroes(MassInventory);
    const layers = [_]LayerInventory{unitLayer()};
    try accumulateLayers(&inventory, &layers, 24, 1, 1);
    try std.testing.expectEqual(@as(f64, 992), inventory.landscape_phosphate_g_p);
    try std.testing.expectEqual(@as(f64, 356), inventory.landscape_ion_mol);
    try std.testing.expectEqual(@as(f64, 356), inventory.grid_cell_ion_mol);
}

test "REDIST soil salt inventory requires both daily gates" {
    var inventory = std.mem.zeroes(MassInventory);
    const layers = [_]LayerInventory{unitLayer()};
    try accumulateLayers(&inventory, &layers, 23, 1, 1);
    try std.testing.expectEqual(std.mem.zeroes(MassInventory), inventory);
    try accumulateLayers(&inventory, &layers, 24, 2, 1);
    try std.testing.expectEqual(std.mem.zeroes(MassInventory), inventory);
}

test "REDIST soil salt inventory runtime layers accumulate independently" {
    var inventory = std.mem.zeroes(MassInventory);
    const layers = [_]LayerInventory{ unitLayer(), unitLayer() };
    try accumulateLayers(&inventory, &layers, 24, 0, 0);
    try std.testing.expectEqual(@as(f64, 1984), inventory.landscape_phosphate_g_p);
    try std.testing.expectEqual(@as(f64, 712), inventory.landscape_ion_mol);
}

test "REDIST soil salt inventory rejects dimensions invalid input and overflow" {
    var inventory = std.mem.zeroes(MassInventory);
    const no_layers: [0]LayerInventory = .{};
    try std.testing.expectError(error.SoilSaltInventoryDimensionMismatch, accumulateLayers(&inventory, &no_layers, 24, 0, 0));
    var layers = [_]LayerInventory{unitLayer()};
    layers[0].micropore.al_mol = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilSaltInventoryLayer, accumulateLayers(&inventory, &layers, 24, 0, 0));
    layers[0] = unitLayer();
    inventory.landscape_ion_mol = std.math.floatMax(f64);
    layers[0].micropore.al_mol = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSoilSaltInventory, accumulateLayers(&inventory, &layers, 24, 0, 0));
}
