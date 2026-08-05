const std = @import("std");

pub const ProfileSource = enum {
    atmospheric_equilibrium,
    supplied_profile,
};

pub const Control = struct {
    profile_source: ProfileSource, // DATA(20) == 'NO' selects equilibrium
    gas_initialization_index: usize, // IGO
};

/// Published aqueous concentrations. Nitrogen values are `mol N m-3` and
/// phosphorus values are `mol P m-3`.
pub const Concentrations = struct {
    ammonium_mol_n_per_m3: f64, // CN4U
    ammonia_mol_n_per_m3: f64, // CN3U
    nitrate_mol_n_per_m3: f64, // CNOU
    dihydrogen_phosphate_mol_p_per_m3: f64, // CH2PU
    hydrogen_phosphate_mol_p_per_m3: f64, // CH1PU
};

/// Zone fractions are used independently exactly as supplied by STARTE.F;
/// this routine does not require paired fractions to sum to one.
pub const ZoneFractions = struct {
    ammonium_non_band: f64, // VLNH4
    ammonium_band: f64, // VLNHB
    nitrate_non_band: f64, // VLNO3
    nitrate_band: f64, // VLNOB
    phosphate_non_band: f64, // VLPO4
    phosphate_band: f64, // VLPOB
};

/// Mineral nutrient inventories are grams of N or P per modeled layer.
pub const Inventories = struct {
    ammonium_non_band_g_n: f64, // ZNH4S
    ammonium_band_g_n: f64, // ZNH4B
    ammonia_non_band_g_n: f64, // ZNH3S
    ammonia_band_g_n: f64, // ZNH3B
    nitrate_non_band_g_n: f64, // ZNO3S
    nitrate_band_g_n: f64, // ZNO3B
    nitrite_non_band_g_n: f64, // ZNO2S; source initializes zero
    nitrite_band_g_n: f64, // ZNO2B; source initializes zero
    dihydrogen_phosphate_non_band_g_p: f64, // H2PO4
    dihydrogen_phosphate_band_g_p: f64, // H2POB
    hydrogen_phosphate_non_band_g_p: f64, // H1PO4
    hydrogen_phosphate_band_g_p: f64, // H1POB
};

/// Direct translation of STARTE.F lines 1443--1454 inside the enclosing
/// `DATA(20) == 'NO' .AND. IGO == 0` branch.
pub fn initialize(
    control: Control,
    concentrations: Concentrations,
    aqueous_volume_m3: f64,
    fractions: ZoneFractions,
) !?Inventories {
    if (control.profile_source != .atmospheric_equilibrium or
        control.gas_initialization_index != 0) return null;
    inline for (@typeInfo(Concentrations).@"struct".fields) |field| {
        const value = @field(concentrations, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidMineralNutrientConcentration;
    }
    if (!std.math.isFinite(aqueous_volume_m3) or aqueous_volume_m3 < 0)
        return error.InvalidAqueousVolume;
    inline for (@typeInfo(ZoneFractions).@"struct".fields) |field| {
        const value = @field(fractions, field.name);
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return error.InvalidMineralNutrientZoneFraction;
    }

    return .{
        .ammonium_non_band_g_n = concentrations.ammonium_mol_n_per_m3 * aqueous_volume_m3 * fractions.ammonium_non_band * 14.0,
        .ammonium_band_g_n = concentrations.ammonium_mol_n_per_m3 * aqueous_volume_m3 * fractions.ammonium_band * 14.0,
        .ammonia_non_band_g_n = concentrations.ammonia_mol_n_per_m3 * aqueous_volume_m3 * fractions.ammonium_non_band * 14.0,
        .ammonia_band_g_n = concentrations.ammonia_mol_n_per_m3 * aqueous_volume_m3 * fractions.ammonium_band * 14.0,
        .nitrate_non_band_g_n = concentrations.nitrate_mol_n_per_m3 * aqueous_volume_m3 * fractions.nitrate_non_band * 14.0,
        .nitrate_band_g_n = concentrations.nitrate_mol_n_per_m3 * aqueous_volume_m3 * fractions.nitrate_band * 14.0,
        .nitrite_non_band_g_n = 0.0,
        .nitrite_band_g_n = 0.0,
        .dihydrogen_phosphate_non_band_g_p = concentrations.dihydrogen_phosphate_mol_p_per_m3 * aqueous_volume_m3 * fractions.phosphate_non_band * 31.0,
        .dihydrogen_phosphate_band_g_p = concentrations.dihydrogen_phosphate_mol_p_per_m3 * aqueous_volume_m3 * fractions.phosphate_band * 31.0,
        .hydrogen_phosphate_non_band_g_p = concentrations.hydrogen_phosphate_mol_p_per_m3 * aqueous_volume_m3 * fractions.phosphate_non_band * 31.0,
        .hydrogen_phosphate_band_g_p = concentrations.hydrogen_phosphate_mol_p_per_m3 * aqueous_volume_m3 * fractions.phosphate_band * 31.0,
    };
}

test "STARTE mineral nutrient zones preserve source conversion and independent fractions" {
    const result = (try initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, .{
        .ammonium_mol_n_per_m3 = 2,
        .ammonia_mol_n_per_m3 = 3,
        .nitrate_mol_n_per_m3 = 4,
        .dihydrogen_phosphate_mol_p_per_m3 = 5,
        .hydrogen_phosphate_mol_p_per_m3 = 6,
    }, 10, .{
        .ammonium_non_band = 0.6,
        .ammonium_band = 0.7,
        .nitrate_non_band = 0.4,
        .nitrate_band = 0.2,
        .phosphate_non_band = 0.3,
        .phosphate_band = 0.8,
    })).?;
    try std.testing.expectEqual(@as(f64, 168), result.ammonium_non_band_g_n);
    try std.testing.expectEqual(@as(f64, 196), result.ammonium_band_g_n);
    try std.testing.expectEqual(@as(f64, 224), result.nitrate_non_band_g_n);
    try std.testing.expectEqual(@as(f64, 465), result.dihydrogen_phosphate_non_band_g_p);
    try std.testing.expectEqual(@as(f64, 1488), result.hydrogen_phosphate_band_g_p);
    try std.testing.expectEqual(@as(f64, 0), result.nitrite_non_band_g_n);
    try std.testing.expectEqual(@as(f64, 0), result.nitrite_band_g_n);
}

test "STARTE inactive mineral nutrient initialization ignores invalid dormant input" {
    const invalid = std.math.nan(f64);
    try std.testing.expectEqual(@as(?Inventories, null), try initialize(.{ .profile_source = .supplied_profile, .gas_initialization_index = 0 }, .{
        .ammonium_mol_n_per_m3 = invalid,
        .ammonia_mol_n_per_m3 = invalid,
        .nitrate_mol_n_per_m3 = invalid,
        .dihydrogen_phosphate_mol_p_per_m3 = invalid,
        .hydrogen_phosphate_mol_p_per_m3 = invalid,
    }, invalid, .{
        .ammonium_non_band = invalid,
        .ammonium_band = invalid,
        .nitrate_non_band = invalid,
        .nitrate_band = invalid,
        .phosphate_non_band = invalid,
        .phosphate_band = invalid,
    }));
}

test "STARTE mineral nutrient zones reject invalid late fraction before publication" {
    try std.testing.expectError(error.InvalidMineralNutrientZoneFraction, initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, .{
        .ammonium_mol_n_per_m3 = 1,
        .ammonia_mol_n_per_m3 = 1,
        .nitrate_mol_n_per_m3 = 1,
        .dihydrogen_phosphate_mol_p_per_m3 = 1,
        .hydrogen_phosphate_mol_p_per_m3 = 1,
    }, 1, .{
        .ammonium_non_band = 0.5,
        .ammonium_band = 0.5,
        .nitrate_non_band = 0.5,
        .nitrate_band = 0.5,
        .phosphate_non_band = 0.5,
        .phosphate_band = std.math.nan(f64),
    }));
}
