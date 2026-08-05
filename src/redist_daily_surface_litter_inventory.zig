const std = @import("std");

pub const OrganicMatter = struct {
    carbon_g: f64,
    colonized_carbon_g: f64,
    nitrogen_g: f64,
    colonized_nitrogen_g: f64,
    phosphorus_g: f64,
    colonized_phosphorus_g: f64,
};

pub const SurfacePools = struct {
    heat_megajoules: f64,
    aqueous_co2_c_g: f64,
    aqueous_ch4_c_g: f64,
    gaseous_co2_c_g: f64,
    gaseous_ch4_c_g: f64,
    aqueous_h2_h_g: f64,
    gaseous_h2_h_g: f64,
    aqueous_oxygen_g: f64,
    gaseous_oxygen_g: f64,
    aqueous_n2_n_g: f64,
    aqueous_n2o_n_g: f64,
    gaseous_n2_n_g: f64,
    gaseous_n2o_n_g: f64,
    aqueous_nh4_n_g: f64,
    aqueous_nh3_n_g: f64,
    gaseous_nh3_n_g: f64,
    exchangeable_nh4_mol: f64,
    fertilizer_nh4_mol: f64,
    fertilizer_urea_mol: f64,
    fertilizer_nh3_mol: f64,
    aqueous_no3_n_g: f64,
    aqueous_no2_n_g: f64,
    fertilizer_no3_mol: f64,
    aqueous_hpo4_p_g: f64,
    aqueous_h2po4_p_g: f64,
    exchangeable_hpo4_mol: f64,
    exchangeable_h2po4_mol: f64,
    precipitated_alpo4_mol: f64,
    precipitated_fepo4_mol: f64,
    precipitated_cahpo4_mol: f64,
    precipitated_cah2po4_mol: f64,
    precipitated_apatite_mol: f64,
};

pub const Parameters = struct { nitrogen_g_mol: f64, phosphorus_g_mol: f64 };

pub const Inventory = struct {
    landscape_litter_carbon_g: f64 = 0,
    cell_litter_carbon_g: f64 = 0,
    landscape_litter_nitrogen_g: f64 = 0,
    cell_litter_nitrogen_g: f64 = 0,
    landscape_litter_phosphorus_g: f64 = 0,
    cell_litter_phosphorus_g: f64 = 0,
    landscape_heat_megajoules: f64 = 0,
    landscape_co2_c_g: f64 = 0,
    cell_co2_c_g: f64 = 0,
    landscape_h2_h_g: f64 = 0,
    landscape_oxygen_g: f64 = 0,
    landscape_n2_n_g: f64 = 0,
    landscape_nh4_n_g: f64 = 0,
    cell_nh4_n_g: f64 = 0,
    landscape_no3_n_g: f64 = 0,
    cell_no3_n_g: f64 = 0,
    landscape_po4_p_g: f64 = 0,
    cell_aqueous_po4_p_g: f64 = 0,
    cell_exchangeable_po4_p_g: f64 = 0,
    cell_precipitated_po4_p_g: f64 = 0,
};

pub fn isDailyFinalSubstep(hour: u8, current_substep: u32, final_substep: u32) !bool {
    if (hour < 1 or hour > 24 or final_substep == 0 or current_substep > final_substep)
        return error.InvalidDailySurfaceInventorySchedule;
    return hour == 24 and current_substep == final_substep;
}

fn finiteStruct(value: anytype) bool {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name))) return false;
    return true;
}

/// Exact standalone translation of REDIST lines 5637--5704.
pub fn accumulate(initial: Inventory, organic: OrganicMatter, pools: SurfacePools, parameters: Parameters) !Inventory {
    if (!finiteStruct(initial) or !finiteStruct(organic) or !finiteStruct(pools) or
        !finiteStruct(parameters) or parameters.nitrogen_g_mol <= 0 or parameters.phosphorus_g_mol <= 0)
        return error.InvalidDailySurfaceInventoryInput;
    var r = initial;
    const carbon = organic.carbon_g + organic.colonized_carbon_g;
    r.landscape_litter_carbon_g += carbon;
    r.cell_litter_carbon_g += carbon;
    const nitrogen = organic.nitrogen_g + organic.colonized_nitrogen_g;
    r.landscape_litter_nitrogen_g += nitrogen;
    r.cell_litter_nitrogen_g += nitrogen;
    const phosphorus = organic.phosphorus_g + organic.colonized_phosphorus_g;
    r.landscape_litter_phosphorus_g += phosphorus;
    r.cell_litter_phosphorus_g += phosphorus;
    r.landscape_heat_megajoules += pools.heat_megajoules;
    const carbon_gas = pools.aqueous_co2_c_g + pools.aqueous_ch4_c_g + pools.gaseous_co2_c_g + pools.gaseous_ch4_c_g;
    r.landscape_co2_c_g += carbon_gas;
    r.cell_co2_c_g += carbon_gas;
    r.landscape_h2_h_g += pools.aqueous_h2_h_g + pools.gaseous_h2_h_g;
    r.landscape_oxygen_g += pools.aqueous_oxygen_g + pools.gaseous_oxygen_g;
    r.landscape_n2_n_g += pools.aqueous_n2_n_g + pools.aqueous_n2o_n_g + pools.gaseous_n2_n_g + pools.gaseous_n2o_n_g;
    const aqueous_nh4 = pools.aqueous_nh4_n_g + pools.aqueous_nh3_n_g + pools.gaseous_nh3_n_g;
    const exchangeable_nh4 = parameters.nitrogen_g_mol * pools.exchangeable_nh4_mol;
    const fertilizer_nh4 = parameters.nitrogen_g_mol * (pools.fertilizer_nh4_mol + pools.fertilizer_urea_mol + pools.fertilizer_nh3_mol);
    r.landscape_nh4_n_g += aqueous_nh4 + exchangeable_nh4 + fertilizer_nh4;
    r.cell_nh4_n_g += aqueous_nh4 + exchangeable_nh4;
    const aqueous_no3 = pools.aqueous_no3_n_g + pools.aqueous_no2_n_g;
    const fertilizer_no3 = parameters.nitrogen_g_mol * pools.fertilizer_no3_mol;
    r.landscape_no3_n_g += aqueous_no3 + fertilizer_no3;
    r.cell_no3_n_g += aqueous_no3;
    const aqueous_po4 = pools.aqueous_hpo4_p_g + pools.aqueous_h2po4_p_g;
    const exchangeable_po4 = parameters.phosphorus_g_mol * (pools.exchangeable_hpo4_mol + pools.exchangeable_h2po4_mol);
    const precipitated_po4 = parameters.phosphorus_g_mol * (pools.precipitated_alpo4_mol + pools.precipitated_fepo4_mol + pools.precipitated_cahpo4_mol) +
        (2.0 * parameters.phosphorus_g_mol) * pools.precipitated_cah2po4_mol +
        (3.0 * parameters.phosphorus_g_mol) * pools.precipitated_apatite_mol;
    r.landscape_po4_p_g += aqueous_po4 + exchangeable_po4 + precipitated_po4;
    r.cell_aqueous_po4_p_g += aqueous_po4;
    r.cell_exchangeable_po4_p_g += exchangeable_po4;
    r.cell_precipitated_po4_p_g += precipitated_po4;
    if (!finiteStruct(r)) return error.NonFiniteDailySurfaceInventory;
    return r;
}

test "REDIST daily surface inventory preserves carrier inclusion and exclusions" {
    var p = std.mem.zeroes(SurfacePools);
    p.aqueous_nh4_n_g = 1;
    p.aqueous_nh3_n_g = 2;
    p.gaseous_nh3_n_g = 3;
    p.exchangeable_nh4_mol = 1;
    p.fertilizer_nh4_mol = 1;
    p.fertilizer_urea_mol = 1;
    p.fertilizer_nh3_mol = 1;
    p.aqueous_no3_n_g = 4;
    p.aqueous_no2_n_g = 5;
    p.fertilizer_no3_mol = 1;
    p.aqueous_hpo4_p_g = 6;
    p.aqueous_h2po4_p_g = 7;
    p.exchangeable_hpo4_mol = 1;
    p.exchangeable_h2po4_mol = 1;
    p.precipitated_alpo4_mol = 1;
    p.precipitated_fepo4_mol = 1;
    p.precipitated_cahpo4_mol = 1;
    p.precipitated_cah2po4_mol = 1;
    p.precipitated_apatite_mol = 1;
    const r = try accumulate(.{}, std.mem.zeroes(OrganicMatter), p, .{ .nitrogen_g_mol = 14, .phosphorus_g_mol = 31 });
    try std.testing.expectEqual(@as(f64, 20), r.cell_nh4_n_g);
    try std.testing.expectEqual(@as(f64, 62), r.landscape_nh4_n_g);
    try std.testing.expectEqual(@as(f64, 9), r.cell_no3_n_g);
    try std.testing.expectEqual(@as(f64, 23), r.landscape_no3_n_g);
    try std.testing.expectEqual(@as(f64, 13), r.cell_aqueous_po4_p_g);
    try std.testing.expectEqual(@as(f64, 62), r.cell_exchangeable_po4_p_g);
    try std.testing.expectEqual(@as(f64, 248), r.cell_precipitated_po4_p_g);
}

test "REDIST daily surface inventory gate and overflow fail fast" {
    try std.testing.expect(try isDailyFinalSubstep(24, 4, 4));
    try std.testing.expect(!try isDailyFinalSubstep(23, 4, 4));
    var initial = Inventory{};
    initial.landscape_heat_megajoules = std.math.floatMax(f64);
    var pools = std.mem.zeroes(SurfacePools);
    pools.heat_megajoules = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteDailySurfaceInventory, accumulate(initial, std.mem.zeroes(OrganicMatter), pools, .{ .nitrogen_g_mol = 14, .phosphorus_g_mol = 31 }));
}

test "REDIST daily surface inventory maps organic heat and gas totals" {
    var pools = std.mem.zeroes(SurfacePools);
    pools.heat_megajoules = 1;
    pools.aqueous_co2_c_g = 2;
    pools.aqueous_ch4_c_g = 3;
    pools.gaseous_co2_c_g = 4;
    pools.gaseous_ch4_c_g = 5;
    pools.aqueous_h2_h_g = 6;
    pools.gaseous_h2_h_g = 7;
    pools.aqueous_oxygen_g = 8;
    pools.gaseous_oxygen_g = 9;
    pools.aqueous_n2_n_g = 10;
    pools.aqueous_n2o_n_g = 11;
    pools.gaseous_n2_n_g = 12;
    pools.gaseous_n2o_n_g = 13;
    const r = try accumulate(.{}, .{ .carbon_g = 1, .colonized_carbon_g = 2, .nitrogen_g = 3, .colonized_nitrogen_g = 4, .phosphorus_g = 5, .colonized_phosphorus_g = 6 }, pools, .{ .nitrogen_g_mol = 14, .phosphorus_g_mol = 31 });
    try std.testing.expectEqual(@as(f64, 3), r.landscape_litter_carbon_g);
    try std.testing.expectEqual(@as(f64, 7), r.cell_litter_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 11), r.landscape_litter_phosphorus_g);
    try std.testing.expectEqual(@as(f64, 14), r.landscape_co2_c_g);
    try std.testing.expectEqual(@as(f64, 13), r.landscape_h2_h_g);
    try std.testing.expectEqual(@as(f64, 17), r.landscape_oxygen_g);
    try std.testing.expectEqual(@as(f64, 46), r.landscape_n2_n_g);
}
