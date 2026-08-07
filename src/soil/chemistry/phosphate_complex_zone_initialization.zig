const std = @import("std");

pub const ProfileSource = enum { atmospheric_equilibrium, supplied_profile };
pub const Control = struct { profile_source: ProfileSource, gas_initialization_index: usize };

/// Phosphate-bearing aqueous species in `mol P m-3`.
pub const Concentrations = struct {
    phosphate: f64, // CH0PU
    phosphoric_acid: f64, // CH3PU
    iron_hydrogen_phosphate: f64, // CF1PU
    iron_dihydrogen_phosphate: f64, // CF2PU
    calcium_phosphate: f64, // CC0PU
    calcium_hydrogen_phosphate: f64, // CC1PU
    calcium_dihydrogen_phosphate: f64, // CC2PU
    magnesium_hydrogen_phosphate: f64, // CM1PU
};

/// Zone inventories retain source mole units; STARTE.F applies no phosphorus
/// atomic-mass conversion to these complex species.
pub const ZoneInventories = Concentrations;

pub const Inventories = struct {
    non_band: ZoneInventories,
    band: ZoneInventories,
    fixed_hydrogen_mol: f64, // ZHYSI
};

/// Direct translation of `starte.f` lines 1502--1518 inside the enclosing
/// `DATA(20) == 'NO' .AND. IGO == 0` branch.
pub fn initialize(
    control: Control,
    concentrations: Concentrations,
    aqueous_volume_m3: f64,
    phosphate_non_band_fraction: f64,
    phosphate_band_fraction: f64,
) !?Inventories {
    if (control.profile_source != .atmospheric_equilibrium or
        control.gas_initialization_index != 0) return null;
    if (!std.math.isFinite(aqueous_volume_m3) or aqueous_volume_m3 < 0)
        return error.InvalidAqueousVolume;
    inline for (.{ phosphate_non_band_fraction, phosphate_band_fraction }) |fraction| {
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.InvalidPhosphateZoneFraction;
    }
    var non_band: ZoneInventories = undefined;
    var band: ZoneInventories = undefined;
    inline for (@typeInfo(Concentrations).@"struct".fields) |field| {
        const value = @field(concentrations, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidPhosphateComplexConcentration;
        @field(non_band, field.name) = value * aqueous_volume_m3 * phosphate_non_band_fraction;
        @field(band, field.name) = value * aqueous_volume_m3 * phosphate_band_fraction;
        if (!std.math.isFinite(@field(non_band, field.name)) or
            !std.math.isFinite(@field(band, field.name)))
            return error.InvalidPhosphateComplexInventory;
    }
    const fixed_hydrogen_mol = 1.0e-3 * aqueous_volume_m3;
    if (!std.math.isFinite(fixed_hydrogen_mol)) return error.InvalidPhosphateComplexInventory;
    return .{ .non_band = non_band, .band = band, .fixed_hydrogen_mol = fixed_hydrogen_mol };
}

fn filled(value: f64) Concentrations {
    var result: Concentrations = undefined;
    inline for (@typeInfo(Concentrations).@"struct".fields) |field| @field(result, field.name) = value;
    return result;
}

test "STARTE phosphate complexes preserve zone multiplication without atomic mass factor" {
    var concentrations = filled(1);
    concentrations.phosphate = 2;
    concentrations.magnesium_hydrogen_phosphate = 3;
    const result = (try initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, concentrations, 10, 0.25, 0.75)).?;
    try std.testing.expectEqual(@as(f64, 5), result.non_band.phosphate);
    try std.testing.expectEqual(@as(f64, 15), result.band.phosphate);
    try std.testing.expectEqual(@as(f64, 7.5), result.non_band.magnesium_hydrogen_phosphate);
    try std.testing.expectEqual(@as(f64, 22.5), result.band.magnesium_hydrogen_phosphate);
    try std.testing.expectEqual(@as(f64, 0.01), result.fixed_hydrogen_mol);
}

test "STARTE phosphate complex zone fractions remain independent" {
    const result = (try initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, filled(1), 2, 0.8, 0.7)).?;
    try std.testing.expectEqual(@as(f64, 1.6), result.non_band.calcium_phosphate);
    try std.testing.expectEqual(@as(f64, 1.4), result.band.calcium_phosphate);
}

test "STARTE inactive phosphate complex initialization ignores invalid dormant input" {
    try std.testing.expectEqual(@as(?Inventories, null), try initialize(.{ .profile_source = .supplied_profile, .gas_initialization_index = 0 }, filled(std.math.nan(f64)), std.math.nan(f64), std.math.nan(f64), std.math.nan(f64)));
}
