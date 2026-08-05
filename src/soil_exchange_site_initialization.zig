const std = @import("std");

pub const ProfileSource = enum { atmospheric_equilibrium, supplied_profile };
pub const Control = struct { profile_source: ProfileSource, gas_initialization_index: usize };

/// Exchange coefficients are `mol Mg-1` and are multiplied by soil mass.
pub const Coefficients = struct {
    ammonium: f64,
    hydrogen: f64,
    aluminum: f64,
    iron: f64,
    calcium: f64,
    magnesium: f64,
    sodium: f64,
    potassium: f64,
    bicarbonate: f64,
    hydroxyl_site_0: f64,
    hydroxyl_site_1: f64,
    hydroxyl_site_2: f64,
    hydrogen_phosphate: f64,
    dihydrogen_phosphate: f64,
};

pub const PhosphateZoneInventories = struct {
    hydroxyl_site_0_mol: f64,
    hydroxyl_site_1_mol: f64,
    hydroxyl_site_2_mol: f64,
    hydrogen_phosphate_mol: f64,
    dihydrogen_phosphate_mol: f64,
};

pub const Inventories = struct {
    ammonium_non_band_mol: f64,
    ammonium_band_mol: f64,
    hydrogen_mol: f64,
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    bicarbonate_mol: f64,
    phosphate_non_band: PhosphateZoneInventories,
    phosphate_band: PhosphateZoneInventories,
};

/// Direct translation of STARTE.F lines 1613--1632 inside the enclosing
/// `DATA(20) == 'NO' .AND. IGO == 0` branch.
pub fn initialize(
    control: Control,
    coefficients: Coefficients,
    soil_mass_megagrams: f64,
    ammonium_non_band_fraction: f64,
    ammonium_band_fraction: f64,
    phosphate_non_band_fraction: f64,
    phosphate_band_fraction: f64,
) !?Inventories {
    if (control.profile_source != .atmospheric_equilibrium or
        control.gas_initialization_index != 0) return null;
    if (!std.math.isFinite(soil_mass_megagrams) or soil_mass_megagrams < 0)
        return error.InvalidSoilMass;
    inline for (.{ ammonium_non_band_fraction, ammonium_band_fraction, phosphate_non_band_fraction, phosphate_band_fraction }) |fraction| {
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.InvalidExchangeZoneFraction;
    }
    inline for (@typeInfo(Coefficients).@"struct".fields) |field| {
        const value = @field(coefficients, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidExchangeCoefficient;
    }
    const result: Inventories = .{
        .ammonium_non_band_mol = coefficients.ammonium * soil_mass_megagrams * ammonium_non_band_fraction,
        .ammonium_band_mol = coefficients.ammonium * soil_mass_megagrams * ammonium_band_fraction,
        .hydrogen_mol = coefficients.hydrogen * soil_mass_megagrams,
        .aluminum_mol = coefficients.aluminum * soil_mass_megagrams,
        .iron_mol = coefficients.iron * soil_mass_megagrams,
        .calcium_mol = coefficients.calcium * soil_mass_megagrams,
        .magnesium_mol = coefficients.magnesium * soil_mass_megagrams,
        .sodium_mol = coefficients.sodium * soil_mass_megagrams,
        .potassium_mol = coefficients.potassium * soil_mass_megagrams,
        .bicarbonate_mol = coefficients.bicarbonate * soil_mass_megagrams,
        .phosphate_non_band = zoneInventories(coefficients, soil_mass_megagrams, phosphate_non_band_fraction),
        .phosphate_band = zoneInventories(coefficients, soil_mass_megagrams, phosphate_band_fraction),
    };
    inline for (@typeInfo(Inventories).@"struct".fields) |field| {
        if (@typeInfo(field.type) == .float and !std.math.isFinite(@field(result, field.name)))
            return error.InvalidExchangeInventory;
    }
    inline for (@typeInfo(PhosphateZoneInventories).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result.phosphate_non_band, field.name)) or
            !std.math.isFinite(@field(result.phosphate_band, field.name)))
            return error.InvalidExchangeInventory;
    }
    return result;
}

fn zoneInventories(coefficients: Coefficients, soil_mass_megagrams: f64, fraction: f64) PhosphateZoneInventories {
    return .{
        .hydroxyl_site_0_mol = coefficients.hydroxyl_site_0 * soil_mass_megagrams * fraction,
        .hydroxyl_site_1_mol = coefficients.hydroxyl_site_1 * soil_mass_megagrams * fraction,
        .hydroxyl_site_2_mol = coefficients.hydroxyl_site_2 * soil_mass_megagrams * fraction,
        .hydrogen_phosphate_mol = coefficients.hydrogen_phosphate * soil_mass_megagrams * fraction,
        .dihydrogen_phosphate_mol = coefficients.dihydrogen_phosphate * soil_mass_megagrams * fraction,
    };
}

fn filled(value: f64) Coefficients {
    var result: Coefficients = undefined;
    inline for (@typeInfo(Coefficients).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

test "STARTE exchange inventories preserve selective zone fractions" {
    var coefficients = filled(1);
    coefficients.ammonium = 2;
    coefficients.calcium = 3;
    coefficients.hydrogen_phosphate = 4;
    const result = (try initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, coefficients, 10, 0.2, 0.7, 0.3, 0.8)).?;
    try std.testing.expectEqual(@as(f64, 4), result.ammonium_non_band_mol);
    try std.testing.expectEqual(@as(f64, 14), result.ammonium_band_mol);
    try std.testing.expectEqual(@as(f64, 30), result.calcium_mol);
    try std.testing.expectEqual(@as(f64, 12), result.phosphate_non_band.hydrogen_phosphate_mol);
    try std.testing.expectEqual(@as(f64, 32), result.phosphate_band.hydrogen_phosphate_mol);
}

test "STARTE exchange zone fractions need not sum to one" {
    const result = (try initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, filled(1), 2, 0.8, 0.7, 0.9, 0.6)).?;
    try std.testing.expectEqual(@as(f64, 1.6), result.ammonium_non_band_mol);
    try std.testing.expectEqual(@as(f64, 1.4), result.ammonium_band_mol);
}

test "STARTE inactive exchange initialization ignores invalid dormant input" {
    try std.testing.expectEqual(@as(?Inventories, null), try initialize(.{ .profile_source = .supplied_profile, .gas_initialization_index = 0 }, filled(std.math.nan(f64)), std.math.nan(f64), std.math.nan(f64), std.math.nan(f64), std.math.nan(f64), std.math.nan(f64)));
}
