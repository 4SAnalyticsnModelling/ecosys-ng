const std = @import("std");

pub const LayerFluxes = struct {
    co2_transformation_g_c: f64, // TRCO2
    water_transformation_mol: f64, // TRH2O
    band_water_transformation_mol: f64, // TBH2O
    nh4_nonband_g_n: f64, // TRN4S
    nh3_nonband_g_n: f64, // TRN3S
    nh4_band_g_n: f64, // TRN4B
    nh3_band_g_n: f64, // TRN3B
    no3_nonband_g_n: f64, // TRNO3
    no3_band_g_n: f64, // TRNOB
    hpo4_nonband_g_p: f64, // TRH1P
    h2po4_nonband_g_p: f64, // TRH2P
    hpo4_band_g_p: f64, // TRH1B
    h2po4_band_g_p: f64, // TRH2B
    hydrogen_exchange_mol: f64, // XZHYS, subtracted
    aluminum_uptake_mol: f64,
    iron_uptake_mol: f64,
    calcium_uptake_mol: f64,
    magnesium_uptake_mol: f64,
    sodium_uptake_mol: f64,
    potassium_uptake_mol: f64,
    sulfate_uptake_mol: f64,
    chloride_uptake_mol: f64,
    aluminum_senescence_mol: f64,
    iron_senescence_mol: f64,
    calcium_senescence_mol: f64,
    magnesium_senescence_mol: f64,
    sodium_senescence_mol: f64,
    potassium_senescence_mol: f64,
    sulfate_senescence_mol: f64,
    chloride_senescence_mol: f64,
};

pub const BoundaryInventory = struct {
    ion_input_mol: f64, // TIONIN
    ion_output_mol: f64, // TIONOU
};

pub const SaltEquilibriumMode = enum { static, dynamic };

fn finiteStruct(value: anytype) bool {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name))) return false;
    return true;
}

/// Direct translation of REDIST 7294--7306 in ascending runtime-layer order.
pub fn accumulateLayers(
    inventory: *BoundaryInventory,
    fluxes_by_layer: []const LayerFluxes,
    mode: SaltEquilibriumMode,
) !void {
    if (fluxes_by_layer.len == 0) return error.SoilSaltBoundaryDimensionMismatch;
    if (!finiteStruct(inventory.*)) return error.InvalidSoilSaltBoundaryInventory;
    for (fluxes_by_layer) |fluxes|
        if (!finiteStruct(fluxes)) return error.InvalidSoilSaltBoundaryFlux;
    if (mode == .static) return;

    var next = inventory.*;
    for (fluxes_by_layer) |f| {
        const output_mol = f.co2_transformation_g_c / 12.0 +
            2.0 * f.water_transformation_mol + 2.0 * f.band_water_transformation_mol +
            (2.0 * f.nh4_nonband_g_n + f.nh3_nonband_g_n) / 14.0 +
            (2.0 * f.nh4_band_g_n + f.nh3_band_g_n) / 14.0 +
            (f.no3_nonband_g_n + f.no3_band_g_n) / 14.0 +
            (2.0 * f.hpo4_nonband_g_p + 3.0 * f.h2po4_nonband_g_p) / 31.0 +
            (2.0 * f.hpo4_band_g_p + 3.0 * f.h2po4_band_g_p) / 31.0 -
            f.hydrogen_exchange_mol + f.aluminum_uptake_mol + f.iron_uptake_mol +
            f.calcium_uptake_mol + f.magnesium_uptake_mol + f.sodium_uptake_mol +
            f.potassium_uptake_mol + f.sulfate_uptake_mol + f.chloride_uptake_mol;
        next.ion_output_mol = next.ion_output_mol + output_mol;
        next.ion_input_mol = next.ion_input_mol + f.aluminum_senescence_mol +
            f.iron_senescence_mol + f.calcium_senescence_mol + f.magnesium_senescence_mol +
            f.sodium_senescence_mol + f.potassium_senescence_mol +
            f.sulfate_senescence_mol + f.chloride_senescence_mol;
        if (!finiteStruct(next)) return error.NonFiniteSoilSaltBoundaryInventory;
    }
    inventory.* = next;
}

fn unitFluxes() LayerFluxes {
    var value: LayerFluxes = undefined;
    inline for (@typeInfo(LayerFluxes).@"struct".fields) |field| @field(value, field.name) = 1.0;
    return value;
}

test "REDIST soil salt boundary preserves every conversion and sign" {
    var inventory = std.mem.zeroes(BoundaryInventory);
    const fluxes = [_]LayerFluxes{unitFluxes()};
    try accumulateLayers(&inventory, &fluxes, .dynamic);
    const expected_output = @as(f64, 1.0) / 12.0 + 4.0 + 3.0 / 14.0 + 3.0 / 14.0 +
        2.0 / 14.0 + 5.0 / 31.0 + 5.0 / 31.0 - 1.0 + 8.0;
    try std.testing.expectApproxEqAbs(expected_output, inventory.ion_output_mol, 8.0 * std.math.floatEps(f64));
    try std.testing.expectEqual(@as(f64, 8.0), inventory.ion_input_mol);
}

test "REDIST soil salt boundary runtime layers accumulate in source order" {
    var inventory = BoundaryInventory{ .ion_input_mol = 0, .ion_output_mol = 1 };
    var fluxes = [_]LayerFluxes{ std.mem.zeroes(LayerFluxes), std.mem.zeroes(LayerFluxes) };
    fluxes[0].co2_transformation_g_c = 12;
    fluxes[1].co2_transformation_g_c = 24;
    try accumulateLayers(&inventory, &fluxes, .dynamic);
    try std.testing.expectEqual(@as(f64, 4), inventory.ion_output_mol);
}

test "REDIST soil salt boundary static mode and failures" {
    var inventory = BoundaryInventory{ .ion_input_mol = 2, .ion_output_mol = 3 };
    const fluxes = [_]LayerFluxes{unitFluxes()};
    try accumulateLayers(&inventory, &fluxes, .static);
    try std.testing.expectEqual(@as(f64, 2), inventory.ion_input_mol);
    try std.testing.expectEqual(@as(f64, 3), inventory.ion_output_mol);

    const no_fluxes: [0]LayerFluxes = .{};
    try std.testing.expectError(error.SoilSaltBoundaryDimensionMismatch, accumulateLayers(&inventory, &no_fluxes, .dynamic));
    var invalid = fluxes;
    invalid[0].chloride_uptake_mol = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilSaltBoundaryFlux, accumulateLayers(&inventory, &invalid, .dynamic));
    inventory.ion_output_mol = std.math.floatMax(f64);
    var overflow = fluxes;
    overflow[0].aluminum_uptake_mol = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSoilSaltBoundaryInventory, accumulateLayers(&inventory, &overflow, .dynamic));
}
