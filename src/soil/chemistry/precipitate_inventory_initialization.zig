const std = @import("std");

pub const ProfileSource = enum { atmospheric_equilibrium, supplied_profile };
pub const Control = struct { profile_source: ProfileSource, gas_initialization_index: usize };

/// Precipitate coefficients are `mol Mg-1`.
pub const Coefficients = struct {
    aluminum_hydroxide: f64,
    iron_hydroxide: f64,
    calcium_carbonate: f64,
    calcium_sulfate: f64,
    aluminum_phosphate: f64,
    iron_phosphate: f64,
    calcium_hydrogen_phosphate: f64,
    calcium_dihydrogen_phosphate: f64,
};

pub const PhosphateZone = struct {
    aluminum_phosphate_mol: f64,
    iron_phosphate_mol: f64,
    calcium_hydrogen_phosphate_mol: f64,
    calcium_dihydrogen_phosphate_mol: f64,
    hydroxyapatite_mol: f64,
};

pub const Inventories = struct {
    aluminum_hydroxide_mol: f64,
    iron_hydroxide_mol: f64,
    calcium_carbonate_mol: f64,
    calcium_sulfate_mol: f64,
    phosphate_non_band: PhosphateZone,
    phosphate_band: PhosphateZone,
    electrical_conductivity: f64, // ECND; exact source initialization
};

fn phosphateZone(coefficients: Coefficients, soil_mass_megagrams: f64, fraction: f64) PhosphateZone {
    return .{
        .aluminum_phosphate_mol = coefficients.aluminum_phosphate * soil_mass_megagrams * fraction,
        .iron_phosphate_mol = coefficients.iron_phosphate * soil_mass_megagrams * fraction,
        .calcium_hydrogen_phosphate_mol = coefficients.calcium_hydrogen_phosphate * soil_mass_megagrams * fraction,
        .calcium_dihydrogen_phosphate_mol = coefficients.calcium_dihydrogen_phosphate * soil_mass_megagrams * fraction,
        .hydroxyapatite_mol = 0.0,
    };
}

/// Direct translation of `starte.f` lines 1672--1686 inside the enclosing
/// `DATA(20) == 'NO' .AND. IGO == 0` branch.
pub fn initialize(control: Control, coefficients: Coefficients, soil_mass_megagrams: f64, phosphate_non_band_fraction: f64, phosphate_band_fraction: f64) !?Inventories {
    if (control.profile_source != .atmospheric_equilibrium or control.gas_initialization_index != 0) return null;
    if (!std.math.isFinite(soil_mass_megagrams) or soil_mass_megagrams < 0) return error.InvalidSoilMass;
    inline for (.{ phosphate_non_band_fraction, phosphate_band_fraction }) |fraction|
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidPrecipitateZoneFraction;
    inline for (@typeInfo(Coefficients).@"struct".fields) |field| {
        const value = @field(coefficients, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidPrecipitateCoefficient;
    }
    const result: Inventories = .{
        .aluminum_hydroxide_mol = coefficients.aluminum_hydroxide * soil_mass_megagrams,
        .iron_hydroxide_mol = coefficients.iron_hydroxide * soil_mass_megagrams,
        .calcium_carbonate_mol = coefficients.calcium_carbonate * soil_mass_megagrams,
        .calcium_sulfate_mol = coefficients.calcium_sulfate * soil_mass_megagrams,
        .phosphate_non_band = phosphateZone(coefficients, soil_mass_megagrams, phosphate_non_band_fraction),
        .phosphate_band = phosphateZone(coefficients, soil_mass_megagrams, phosphate_band_fraction),
        .electrical_conductivity = 0.0,
    };
    inline for (.{ result.aluminum_hydroxide_mol, result.iron_hydroxide_mol, result.calcium_carbonate_mol, result.calcium_sulfate_mol }) |value|
        if (!std.math.isFinite(value)) return error.InvalidPrecipitateInventory;
    inline for (@typeInfo(PhosphateZone).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.phosphate_non_band, field.name)) or !std.math.isFinite(@field(result.phosphate_band, field.name))) return error.InvalidPrecipitateInventory;
    return result;
}

fn filled(value: f64) Coefficients {
    var result: Coefficients = undefined;
    inline for (@typeInfo(Coefficients).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

test "STARTE precipitates preserve whole-soil and independent zone scaling" {
    var coefficients = filled(1);
    coefficients.calcium_carbonate = 2;
    coefficients.calcium_hydrogen_phosphate = 3;
    const result = (try initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, coefficients, 10, 0.2, 0.7)).?;
    try std.testing.expectEqual(@as(f64, 20), result.calcium_carbonate_mol);
    try std.testing.expectEqual(@as(f64, 6), result.phosphate_non_band.calcium_hydrogen_phosphate_mol);
    try std.testing.expectEqual(@as(f64, 21), result.phosphate_band.calcium_hydrogen_phosphate_mol);
    try std.testing.expectEqual(@as(f64, 0), result.phosphate_non_band.hydroxyapatite_mol);
    try std.testing.expectEqual(@as(f64, 0), result.phosphate_band.hydroxyapatite_mol);
    try std.testing.expectEqual(@as(f64, 0), result.electrical_conductivity);
}

test "STARTE precipitate zone fractions need not sum to one" {
    const result = (try initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, filled(1), 2, 0.8, 0.7)).?;
    try std.testing.expectEqual(@as(f64, 1.6), result.phosphate_non_band.aluminum_phosphate_mol);
    try std.testing.expectEqual(@as(f64, 1.4), result.phosphate_band.aluminum_phosphate_mol);
}

test "STARTE inactive precipitate initialization ignores invalid dormant input" {
    const nan = std.math.nan(f64);
    try std.testing.expectEqual(@as(?Inventories, null), try initialize(.{ .profile_source = .supplied_profile, .gas_initialization_index = 0 }, filled(nan), nan, nan, nan));
}
