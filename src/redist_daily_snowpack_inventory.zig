const std = @import("std");
const snow_salt = @import("redist_snow_redistribution_salt_update.zig");

pub const SnowLayer = struct {
    solid_water_equivalent_m3: f64,
    liquid_water_m3: f64,
    vapor_water_equivalent_m3: f64,
    ice_volume_m3: f64,
    heat_capacity_megajoules_k: f64,
    temperature_k: f64,
    co2_c_g: f64,
    ch4_c_g: f64,
    oxygen_g: f64,
    n2_n_g: f64,
    n2o_n_g: f64,
    nh4_n_g: f64,
    nh3_n_g: f64,
    no3_n_g: f64,
    hpo4_p_g: f64,
    h2po4_p_g: f64,
    ions: snow_salt.SnowIonPools,
};

pub const Inventory = struct {
    landscape_water_m3: f64 = 0.0,
    cell_water_m3: f64 = 0.0,
    landscape_heat_megajoules: f64 = 0.0,
    landscape_co2_c_g: f64 = 0.0,
    cell_co2_c_g: f64 = 0.0,
    landscape_oxygen_g: f64 = 0.0,
    landscape_n2_n_g: f64 = 0.0,
    landscape_nh4_n_g: f64 = 0.0,
    landscape_no3_n_g: f64 = 0.0,
    landscape_po4_p_g: f64 = 0.0,
    landscape_ion_inventory_mol: f64 = 0.0,
};

pub fn isDailyFinalSubstep(hour: u8, current_substep: u32, final_substep: u32) !bool {
    if (hour < 1 or hour > 24 or final_substep == 0 or current_substep > final_substep)
        return error.InvalidSnowInventorySchedule;
    return hour == 24 and current_substep == final_substep;
}

fn validateLayer(layer: SnowLayer) !void {
    inline for (@typeInfo(SnowLayer).@"struct".fields) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(layer, field.name)))
            return error.InvalidSnowInventoryLayer;
    }
    inline for (@typeInfo(snow_salt.SnowIonPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(layer.ions, field.name)))
            return error.InvalidSnowInventoryLayer;
}

/// REDIST lines 5485--5500, retaining the original charge-class grouping.
pub fn stoichiometricIonInventory(ions: snow_salt.SnowIonPools) !f64 {
    inline for (@typeInfo(snow_salt.SnowIonPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(ions, field.name))) return error.InvalidSnowInventoryLayer;
    const value = ions.al_mol + ions.fe_mol + ions.h_mol + ions.ca_mol +
        ions.megagrams_mol + ions.na_mol + ions.ka_mol + ions.oh_mol + ions.so4_mol +
        ions.cl_mol + ions.co3_mol + ions.h0p_mol +
        2.0 * (ions.hco3_mol + ions.aloh1_mol + ions.als_mol + ions.feoh1_mol +
            ions.fes_mol + ions.cao_mol + ions.cac_mol + ions.cas_mol +
            ions.mgo_mol + ions.mgc_mol + ions.mgs_mol + ions.nac_mol +
            ions.nas_mol + ions.kas_mol + ions.ca0p_mol) +
        3.0 * (ions.aloh2_mol + ions.feoh2_mol + ions.cah_mol + ions.mgh_mol +
            ions.fe1p_mol + ions.ca1p_mol + ions.mg1p_mol) +
        4.0 * (ions.aloh3_mol + ions.feoh3_mol + ions.h3p_mol +
            ions.fe2p_mol + ions.ca2p_mol) +
        5.0 * (ions.aloh4_mol + ions.feoh4_mol);
    if (!std.math.isFinite(value)) return error.NonFiniteSnowInventory;
    return value;
}

/// Exact standalone translation of REDIST lines 5470--5518.
/// `ice_density_megagrams_m3` retains the legacy DENSI arithmetic exactly.
pub fn accumulate(
    initial: Inventory,
    layers: []const SnowLayer,
    ice_density_megagrams_m3: f64,
    salt_mode: snow_salt.SaltSimulationMode,
) !Inventory {
    if (!std.math.isFinite(ice_density_megagrams_m3)) return error.InvalidSnowIceDensity;
    var result = initial;
    inline for (@typeInfo(Inventory).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name))) return error.InvalidSnowInventory;

    for (layers) |layer| {
        try validateLayer(layer);
        const water = layer.solid_water_equivalent_m3 + layer.liquid_water_m3 +
            layer.vapor_water_equivalent_m3 + layer.ice_volume_m3 * ice_density_megagrams_m3;
        result.landscape_water_m3 += water;
        result.cell_water_m3 += water;
        const energy = layer.heat_capacity_megajoules_k * layer.temperature_k;
        result.landscape_heat_megajoules += energy;
        result.landscape_co2_c_g += layer.co2_c_g + layer.ch4_c_g;
        result.cell_co2_c_g += layer.co2_c_g + layer.ch4_c_g;
        result.landscape_oxygen_g += layer.oxygen_g;
        result.landscape_n2_n_g += layer.n2_n_g + layer.n2o_n_g;
        result.landscape_nh4_n_g += layer.nh4_n_g + layer.nh3_n_g;
        result.landscape_no3_n_g += layer.no3_n_g;
        result.landscape_po4_p_g += layer.hpo4_p_g + layer.h2po4_p_g;
        if (salt_mode == .dynamic)
            result.landscape_ion_inventory_mol += try stoichiometricIonInventory(layer.ions);
        inline for (@typeInfo(Inventory).@"struct".fields) |field|
            if (!std.math.isFinite(@field(result, field.name)))
                return error.NonFiniteSnowInventory;
    }
    return result;
}

fn zeroLayer() SnowLayer {
    return std.mem.zeroes(SnowLayer);
}

test "REDIST daily snow inventory preserves source arithmetic and layer order" {
    var layers = [_]SnowLayer{ zeroLayer(), zeroLayer() };
    layers[0].solid_water_equivalent_m3 = 1.0;
    layers[0].liquid_water_m3 = 2.0;
    layers[0].vapor_water_equivalent_m3 = 3.0;
    layers[0].ice_volume_m3 = 4.0;
    layers[0].heat_capacity_megajoules_k = 2.0;
    layers[0].temperature_k = 5.0;
    layers[0].co2_c_g = 1.0;
    layers[0].ch4_c_g = 2.0;
    layers[0].oxygen_g = 3.0;
    layers[0].n2_n_g = 4.0;
    layers[0].n2o_n_g = 5.0;
    layers[0].nh4_n_g = 6.0;
    layers[0].nh3_n_g = 7.0;
    layers[0].no3_n_g = 8.0;
    layers[0].hpo4_p_g = 9.0;
    layers[0].h2po4_p_g = 10.0;
    layers[1].solid_water_equivalent_m3 = 1.0;
    const result = try accumulate(.{}, &layers, 0.5, .static_equilibrium);
    try std.testing.expectEqual(@as(f64, 9.0), result.landscape_water_m3);
    try std.testing.expectEqual(result.landscape_water_m3, result.cell_water_m3);
    try std.testing.expectEqual(@as(f64, 10.0), result.landscape_heat_megajoules);
    try std.testing.expectEqual(@as(f64, 3.0), result.landscape_co2_c_g);
    try std.testing.expectEqual(@as(f64, 9.0), result.landscape_n2_n_g);
    try std.testing.expectEqual(@as(f64, 13.0), result.landscape_nh4_n_g);
    try std.testing.expectEqual(@as(f64, 19.0), result.landscape_po4_p_g);
}

test "REDIST snow ion inventory applies every charge class" {
    var ions = std.mem.zeroes(snow_salt.SnowIonPools);
    ions.al_mol = 1.0;
    ions.hco3_mol = 1.0;
    ions.aloh2_mol = 1.0;
    ions.aloh3_mol = 1.0;
    ions.aloh4_mol = 1.0;
    try std.testing.expectEqual(@as(f64, 15.0), try stoichiometricIonInventory(ions));
}

test "REDIST daily snow inventory gate is exact" {
    try std.testing.expect(try isDailyFinalSubstep(24, 4, 4));
    try std.testing.expect(!try isDailyFinalSubstep(23, 4, 4));
    try std.testing.expect(!try isDailyFinalSubstep(24, 3, 4));
    try std.testing.expectError(error.InvalidSnowInventorySchedule, isDailyFinalSubstep(25, 4, 4));
}

test "REDIST daily snow inventory rejects invalid and overflow" {
    var layer = zeroLayer();
    layer.temperature_k = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSnowInventoryLayer,
        accumulate(.{}, &.{layer}, 0.9, .dynamic),
    );
    layer = zeroLayer();
    layer.solid_water_equivalent_m3 = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSnowInventory,
        accumulate(.{ .landscape_water_m3 = std.math.floatMax(f64) }, &.{layer}, 0.9, .static_equilibrium),
    );
}
