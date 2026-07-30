const std = @import("std");

/// Elemental mass changes consumed by redistribution.  The molecular species
/// remain explicit so a gram of N or P cannot accidentally be treated as a
/// mole in a later transport kernel.
pub const MolarChanges = struct {
    carbon_dioxide_mol: f64,
    ammonium_mol_n: f64,
    ammonia_mol_n: f64,
    nitrate_mol_n: f64,
    hydrogen_phosphate_mol_p: f64,
    dihydrogen_phosphate_mol_p: f64,
};

pub const ElementMassChanges = struct {
    carbon_dioxide_g_c: f64,
    ammonium_g_n: f64,
    ammonia_g_n: f64,
    nitrate_g_n: f64,
    hydrogen_phosphate_g_p: f64,
    dihydrogen_phosphate_g_p: f64,
};

pub const carbon_molar_mass_g_per_mol = 12.0;
pub const nitrogen_molar_mass_g_per_mol = 14.0;
pub const phosphorus_molar_mass_g_per_mol = 31.0;

/// Implements the final SOLUTE.F REDIST conversion.  These are elemental
/// masses (C, N, and P), not whole-ion molecular masses.
pub fn toElementMass(changes: MolarChanges) !ElementMassChanges {
    inline for (@typeInfo(MolarChanges).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(changes, field.name))) {
            return error.NonFiniteSoluteMolarChange;
        }
    }
    return .{
        .carbon_dioxide_g_c = changes.carbon_dioxide_mol * carbon_molar_mass_g_per_mol,
        .ammonium_g_n = changes.ammonium_mol_n * nitrogen_molar_mass_g_per_mol,
        .ammonia_g_n = changes.ammonia_mol_n * nitrogen_molar_mass_g_per_mol,
        .nitrate_g_n = changes.nitrate_mol_n * nitrogen_molar_mass_g_per_mol,
        .hydrogen_phosphate_g_p = changes.hydrogen_phosphate_mol_p * phosphorus_molar_mass_g_per_mol,
        .dihydrogen_phosphate_g_p = changes.dihydrogen_phosphate_mol_p * phosphorus_molar_mass_g_per_mol,
    };
}

test "SOLUTE elemental conversions retain the Fortran factors and signs" {
    const result = try toElementMass(.{
        .carbon_dioxide_mol = 2,
        .ammonium_mol_n = -3,
        .ammonia_mol_n = 4,
        .nitrate_mol_n = 5,
        .hydrogen_phosphate_mol_p = -6,
        .dihydrogen_phosphate_mol_p = 7,
    });
    try std.testing.expectEqual(@as(f64, 24), result.carbon_dioxide_g_c);
    try std.testing.expectEqual(@as(f64, -42), result.ammonium_g_n);
    try std.testing.expectEqual(@as(f64, 56), result.ammonia_g_n);
    try std.testing.expectEqual(@as(f64, 70), result.nitrate_g_n);
    try std.testing.expectEqual(@as(f64, -186), result.hydrogen_phosphate_g_p);
    try std.testing.expectEqual(@as(f64, 217), result.dihydrogen_phosphate_g_p);
}

test "non-finite molar changes fail before redistribution" {
    try std.testing.expectError(error.NonFiniteSoluteMolarChange, toElementMass(.{
        .carbon_dioxide_mol = std.math.nan(f64),
        .ammonium_mol_n = 0,
        .ammonia_mol_n = 0,
        .nitrate_mol_n = 0,
        .hydrogen_phosphate_mol_p = 0,
        .dihydrogen_phosphate_mol_p = 0,
    }));
}
