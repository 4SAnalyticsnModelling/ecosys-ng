const std = @import("std");

/// Internal fire products accumulated across all combusted soil layers.
pub const CombustionTotals = struct {
    organic_carbon_retained_g_c: f64, // COM
    carbon_dioxide_emitted_g_c: f64, // COX
    methane_emitted_g_c: f64, // CHX
    gaseous_nitrogen_emitted_g_n: f64, // ZOX
    gaseous_phosphorus_emitted_g_p: f64, // POX
    mineral_nitrogen_retained_g_n: f64, // Z4M
    mineral_phosphorus_retained_g_p: f64, // P4M
};

/// Cumulative landscape-owned inventory and surface-exchange ledgers.
pub const LandscapeState = struct {
    surface_gas_nitrogen_input_g_n: f64, // ZN2GIN; emission is negative input
    surface_phosphorus_input_g_p: f64, // TPIN; emission is negative input
    litter_carbon_g_c: f64, // TLRSDC
    ammonium_g_n: f64, // TLNH4
    phosphate_g_p: f64, // TLPO4
};

/// Cell-owned signed boundary fluxes for the current model step.
pub const CellBoundaryState = struct {
    fire_soil_organic_carbon_loss_g_c_step: f64, // VCOXFS
    carbon_dioxide_input_g_c_step: f64, // VCO2G; emission is negative input
    methane_input_g_c_step: f64, // VCH4G; emission is negative input
    nitrogen_oxide_input_g_n_step: f64, // VNOXG; emission is negative input
    phosphorus_oxide_input_g_p_step: f64, // VPOXG; emission is negative input
};

fn finiteStruct(value: anytype) bool {
    inline for (std.meta.fields(@TypeOf(value))) |field|
        if (!std.math.isFinite(@field(value, field.name))) return false;
    return true;
}

/// Direct translation of REDIST 10919--10928.
/// Landscape inventories and cell boundary fluxes retain distinct ownership.
pub fn publish(totals: CombustionTotals, landscape: *LandscapeState, cell: *CellBoundaryState) !void {
    if (!finiteStruct(totals) or !finiteStruct(landscape.*) or !finiteStruct(cell.*))
        return error.InvalidFireBoundaryLedgerInput;
    inline for (std.meta.fields(CombustionTotals)) |field|
        if (@field(totals, field.name) < 0) return error.InvalidFireBoundaryLedgerInput;

    var next_landscape = landscape.*;
    var next_cell = cell.*;
    next_landscape.surface_gas_nitrogen_input_g_n = next_landscape.surface_gas_nitrogen_input_g_n - totals.gaseous_nitrogen_emitted_g_n;
    next_landscape.surface_phosphorus_input_g_p = next_landscape.surface_phosphorus_input_g_p - totals.gaseous_phosphorus_emitted_g_p;
    next_landscape.litter_carbon_g_c = next_landscape.litter_carbon_g_c + totals.organic_carbon_retained_g_c;
    next_landscape.ammonium_g_n = next_landscape.ammonium_g_n + totals.mineral_nitrogen_retained_g_n;
    next_landscape.phosphate_g_p = next_landscape.phosphate_g_p + totals.mineral_phosphorus_retained_g_p;
    next_cell.fire_soil_organic_carbon_loss_g_c_step = next_cell.fire_soil_organic_carbon_loss_g_c_step + totals.organic_carbon_retained_g_c;
    next_cell.carbon_dioxide_input_g_c_step = next_cell.carbon_dioxide_input_g_c_step - totals.carbon_dioxide_emitted_g_c;
    next_cell.methane_input_g_c_step = next_cell.methane_input_g_c_step - totals.methane_emitted_g_c;
    next_cell.nitrogen_oxide_input_g_n_step = next_cell.nitrogen_oxide_input_g_n_step - totals.gaseous_nitrogen_emitted_g_n;
    next_cell.phosphorus_oxide_input_g_p_step = next_cell.phosphorus_oxide_input_g_p_step - totals.gaseous_phosphorus_emitted_g_p;
    if (!finiteStruct(next_landscape) or !finiteStruct(next_cell))
        return error.NonFiniteFireBoundaryLedgerResult;
    landscape.* = next_landscape;
    cell.* = next_cell;
}

test "REDIST fire boundary ledgers preserve exact signs and ownership" {
    const totals = CombustionTotals{
        .organic_carbon_retained_g_c = 1,
        .carbon_dioxide_emitted_g_c = 2,
        .methane_emitted_g_c = 3,
        .gaseous_nitrogen_emitted_g_n = 4,
        .gaseous_phosphorus_emitted_g_p = 5,
        .mineral_nitrogen_retained_g_n = 6,
        .mineral_phosphorus_retained_g_p = 7,
    };
    var landscape = LandscapeState{ .surface_gas_nitrogen_input_g_n = 10, .surface_phosphorus_input_g_p = 20, .litter_carbon_g_c = 30, .ammonium_g_n = 40, .phosphate_g_p = 50 };
    var cell = CellBoundaryState{ .fire_soil_organic_carbon_loss_g_c_step = 60, .carbon_dioxide_input_g_c_step = 70, .methane_input_g_c_step = 80, .nitrogen_oxide_input_g_n_step = 90, .phosphorus_oxide_input_g_p_step = 100 };
    try publish(totals, &landscape, &cell);
    try std.testing.expectEqual(@as(f64, 6), landscape.surface_gas_nitrogen_input_g_n);
    try std.testing.expectEqual(@as(f64, 15), landscape.surface_phosphorus_input_g_p);
    try std.testing.expectEqual(@as(f64, 31), landscape.litter_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 46), landscape.ammonium_g_n);
    try std.testing.expectEqual(@as(f64, 57), landscape.phosphate_g_p);
    try std.testing.expectEqual(@as(f64, 61), cell.fire_soil_organic_carbon_loss_g_c_step);
    try std.testing.expectEqual(@as(f64, 68), cell.carbon_dioxide_input_g_c_step);
    try std.testing.expectEqual(@as(f64, 77), cell.methane_input_g_c_step);
    try std.testing.expectEqual(@as(f64, 86), cell.nitrogen_oxide_input_g_n_step);
    try std.testing.expectEqual(@as(f64, 95), cell.phosphorus_oxide_input_g_p_step);
}

test "REDIST fire boundary ledger late overflow is atomic" {
    const totals = CombustionTotals{ .organic_carbon_retained_g_c = 1, .carbon_dioxide_emitted_g_c = 2, .methane_emitted_g_c = 3, .gaseous_nitrogen_emitted_g_n = 4, .gaseous_phosphorus_emitted_g_p = 5, .mineral_nitrogen_retained_g_n = 6, .mineral_phosphorus_retained_g_p = std.math.floatMax(f64) };
    var landscape = LandscapeState{ .surface_gas_nitrogen_input_g_n = 10, .surface_phosphorus_input_g_p = 20, .litter_carbon_g_c = 30, .ammonium_g_n = 40, .phosphate_g_p = std.math.floatMax(f64) };
    var cell = CellBoundaryState{ .fire_soil_organic_carbon_loss_g_c_step = 60, .carbon_dioxide_input_g_c_step = 70, .methane_input_g_c_step = 80, .nitrogen_oxide_input_g_n_step = 90, .phosphorus_oxide_input_g_p_step = 100 };
    try std.testing.expectError(error.NonFiniteFireBoundaryLedgerResult, publish(totals, &landscape, &cell));
    try std.testing.expectEqual(@as(f64, 10), landscape.surface_gas_nitrogen_input_g_n);
    try std.testing.expectEqual(std.math.floatMax(f64), landscape.phosphate_g_p);
    try std.testing.expectEqual(@as(f64, 60), cell.fire_soil_organic_carbon_loss_g_c_step);
}
