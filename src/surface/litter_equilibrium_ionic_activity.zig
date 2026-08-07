const std = @import("std");

pub const Species = struct {
    aluminum: f64,
    iron: f64,
    phosphate: f64,
    calcium: f64,
    magnesium: f64,
    aluminum_monohydroxide: f64,
    iron_monohydroxide: f64,
    iron_dihydrogen_phosphate: f64,
    sulfate: f64,
    carbonate: f64,
    hydrogen_phosphate: f64,
    ammonium: f64,
    hydrogen: f64,
    sodium: f64,
    potassium: f64,
    aluminum_dihydroxide: f64,
    iron_dihydroxide: f64,
    aluminum_sulfate: f64,
    iron_sulfate: f64,
    calcium_hydroxide: f64,
    calcium_bicarbonate: f64,
    magnesium_hydroxide: f64,
    magnesium_bicarbonate: f64,
    iron_hydrogen_phosphate: f64,
    calcium_dihydrogen_phosphate: f64,
    hydroxide: f64,
    bicarbonate: f64,
    chloride: f64,
    aluminum_tetrahydroxide: f64,
    iron_tetrahydroxide: f64,
    sodium_carbonate: f64,
    sodium_sulfate: f64,
    potassium_sulfate: f64,
    dihydrogen_phosphate: f64,
    calcium_phosphate: f64,
};

pub const ChargeTotals = struct {
    trivalent_cations_mol_per_m3: f64,
    trivalent_anions_mol_per_m3: f64,
    divalent_cations_mol_per_m3: f64,
    divalent_anions_mol_per_m3: f64,
    monovalent_cations_mol_per_m3: f64,
    monovalent_anions_mol_per_m3: f64,
};

pub const Coefficients = struct { monovalent: f64, divalent: f64, trivalent: f64 };
pub const Activities = struct {
    hydrogen: f64,
    hydroxide: f64,
    ammonium: f64,
    aluminum: f64,
    iron: f64,
    calcium: f64,
    magnesium: f64,
    sodium: f64,
    potassium: f64,
};
pub const Result = struct {
    charge_totals: ChargeTotals,
    ionic_strength_mol_per_l: f64,
    activity_expression: f64,
    coefficients: Coefficients,
    activities_mol_per_m3: Activities,
};

/// Direct translation of one pass through `starte.f` lines 1787--1810.
pub fn calculate(dynamic_salt_enabled: bool, species: Species, surface_litter_ph: f64, water_ion_product: f64) !?Result {
    if (!dynamic_salt_enabled) return null;
    inline for (@typeInfo(Species).@"struct".fields) |field| {
        const value = @field(species, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSurfaceLitterSpecies;
    }
    if (!std.math.isFinite(surface_litter_ph) or !std.math.isFinite(water_ion_product) or water_ion_product <= 0)
        return error.InvalidSurfaceLitterActivityParameter;

    const totals: ChargeTotals = .{
        .trivalent_cations_mol_per_m3 = species.aluminum + species.iron,
        .trivalent_anions_mol_per_m3 = species.phosphate,
        .divalent_cations_mol_per_m3 = species.calcium + species.magnesium + species.aluminum_monohydroxide + species.iron_monohydroxide + species.iron_dihydrogen_phosphate,
        .divalent_anions_mol_per_m3 = species.sulfate + species.carbonate + species.hydrogen_phosphate,
        .monovalent_cations_mol_per_m3 = species.ammonium + species.hydrogen + species.sodium + species.potassium + species.aluminum_dihydroxide + species.iron_dihydroxide + species.aluminum_sulfate + species.iron_sulfate + species.calcium_hydroxide + species.calcium_bicarbonate + species.magnesium_hydroxide + species.magnesium_bicarbonate + species.iron_hydrogen_phosphate + species.calcium_dihydrogen_phosphate,
        .monovalent_anions_mol_per_m3 = species.hydroxide + species.bicarbonate + species.chloride + species.aluminum_tetrahydroxide + species.iron_tetrahydroxide + species.sodium_carbonate + species.sodium_sulfate + species.potassium_sulfate + species.dihydrogen_phosphate + species.calcium_phosphate,
    };
    const ionic_strength = 1.0e-3 * (9.0 * (totals.trivalent_cations_mol_per_m3 + totals.trivalent_anions_mol_per_m3) + 4.0 * (totals.divalent_cations_mol_per_m3 + totals.divalent_anions_mol_per_m3) + totals.monovalent_cations_mol_per_m3 + totals.monovalent_anions_mol_per_m3);
    const strength_root = @sqrt(ionic_strength);
    const activity_expression = strength_root / (1.0 + strength_root) - 0.20 * ionic_strength;
    const coefficients: Coefficients = .{
        .monovalent = @min(1.0, std.math.pow(f64, 10.0, -0.509 * 1.0 * activity_expression)),
        .divalent = @min(1.0, std.math.pow(f64, 10.0, -0.509 * 4.0 * activity_expression)),
        .trivalent = @min(1.0, std.math.pow(f64, 10.0, -0.509 * 9.0 * activity_expression)),
    };
    const hydrogen = std.math.pow(f64, 10.0, -(surface_litter_ph - 3.0));
    const activities: Activities = .{
        .hydrogen = hydrogen,
        .hydroxide = water_ion_product / hydrogen,
        .ammonium = species.ammonium * coefficients.monovalent,
        .aluminum = species.aluminum * coefficients.trivalent,
        .iron = species.iron * coefficients.trivalent,
        .calcium = species.calcium * coefficients.divalent,
        .magnesium = species.magnesium * coefficients.divalent,
        .sodium = species.sodium * coefficients.monovalent,
        .potassium = species.potassium * coefficients.monovalent,
    };
    inline for (@typeInfo(Activities).@"struct".fields) |field|
        if (!std.math.isFinite(@field(activities, field.name))) return error.InvalidSurfaceLitterActivityResult;
    return .{ .charge_totals = totals, .ionic_strength_mol_per_l = ionic_strength, .activity_expression = activity_expression, .coefficients = coefficients, .activities_mol_per_m3 = activities };
}

fn filled(value: f64) Species {
    var result: Species = undefined;
    inline for (@typeInfo(Species).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

test "STARTE surface litter ionic activity preserves charge-class source order" {
    const species = filled(1);
    const result = (try calculate(true, species, 6, 1.0e-14)).?;
    try std.testing.expectEqual(@as(f64, 2), result.charge_totals.trivalent_cations_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), result.charge_totals.trivalent_anions_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 5), result.charge_totals.divalent_cations_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 3), result.charge_totals.divalent_anions_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 14), result.charge_totals.monovalent_cations_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 10), result.charge_totals.monovalent_anions_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 0.083), result.ionic_strength_mol_per_l);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0e-3), result.activities_mol_per_m3.hydrogen, 1e-18);
}

test "STARTE activity coefficients retain one cap at zero ionic strength" {
    const result = (try calculate(true, filled(0), 7, 1.0e-14)).?;
    try std.testing.expectEqual(@as(f64, 1), result.coefficients.monovalent);
    try std.testing.expectEqual(@as(f64, 1), result.coefficients.divalent);
    try std.testing.expectEqual(@as(f64, 1), result.coefficients.trivalent);
}

test "STARTE inactive surface litter activity ignores invalid dormant inputs" {
    try std.testing.expectEqual(@as(?Result, null), try calculate(false, filled(std.math.nan(f64)), std.math.nan(f64), std.math.nan(f64)));
}
