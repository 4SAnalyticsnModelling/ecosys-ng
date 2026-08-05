const std = @import("std");

/// Direct translation of STARTE.F lines 1907--1910. `XHC1` is `mol Mg-1`,
/// `BKVL(0)` is Mg, and the published `XHC(0)` inventory is mol.
pub fn publish(
    dynamic_salt_enabled: bool,
    occupied_carboxyl_hydrogen_mol_per_megagram: f64,
    surface_litter_mass_megagrams: f64,
) !f64 {
    if (!dynamic_salt_enabled) return 0.0;
    if (!std.math.isFinite(occupied_carboxyl_hydrogen_mol_per_megagram) or
        occupied_carboxyl_hydrogen_mol_per_megagram < 0)
        return error.InvalidCarboxylHydrogenConcentration;
    if (!std.math.isFinite(surface_litter_mass_megagrams) or surface_litter_mass_megagrams < 0)
        return error.InvalidSurfaceLitterMass;
    const inventory_mol = occupied_carboxyl_hydrogen_mol_per_megagram * surface_litter_mass_megagrams;
    if (!std.math.isFinite(inventory_mol)) return error.InvalidCarboxylHydrogenInventory;
    return inventory_mol;
}

test "STARTE publishes converged carboxyl hydrogen as extensive inventory" {
    try std.testing.expectEqual(@as(f64, 6), try publish(true, 2, 3));
}

test "STARTE static salt branch forces exact zero before dormant validation" {
    try std.testing.expectEqual(@as(f64, 0), try publish(false, std.math.nan(f64), std.math.nan(f64)));
}

test "STARTE carboxyl hydrogen publication rejects overflow" {
    try std.testing.expectError(error.InvalidCarboxylHydrogenInventory, publish(true, std.math.floatMax(f64), 2));
}
