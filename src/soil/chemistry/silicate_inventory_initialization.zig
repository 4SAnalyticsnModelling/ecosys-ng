const std = @import("std");

pub const ProfileSource = enum { atmospheric_equilibrium, supplied_profile };
pub const Control = struct { profile_source: ProfileSource, gas_initialization_index: usize };

pub const Texture = struct {
    sand_g_per_megagram: f64, // CSAND
    silt_g_per_megagram: f64, // CSILT
    clay_g_per_megagram: f64, // CCLAY
};

pub const ElementInventories = struct {
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
};

pub const Inventories = struct {
    mineral_surface_area_m2: f64, // SSAL
    soil_silicates: ElementInventories, // Q*SI
    ground_rock_silicates: ElementInventories, // Q*SIF
};

fn zeroElements() ElementInventories {
    return .{ .aluminum_mol = 0, .iron_mol = 0, .calcium_mol = 0, .magnesium_mol = 0, .sodium_mol = 0, .potassium_mol = 0 };
}

/// Direct translation of `starte.f` lines 1642--1662 inside the enclosing
/// `DATA(20) == 'NO' .AND. IGO == 0` branch.
pub fn initialize(control: Control, texture: Texture, soil_mass_megagrams: f64, ph: f64) !?Inventories {
    if (control.profile_source != .atmospheric_equilibrium or
        control.gas_initialization_index != 0) return null;
    inline for (@typeInfo(Texture).@"struct".fields) |field| {
        const value = @field(texture, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSoilTexture;
    }
    if (!std.math.isFinite(soil_mass_megagrams) or soil_mass_megagrams < 0) return error.InvalidSoilMass;
    if (!std.math.isFinite(ph) or ph < 0 or ph > 14) return error.InvalidSoilPh;

    const surface_area_m2 = (0.3 * texture.sand_g_per_megagram + 2.2 * texture.silt_g_per_megagram +
        8.0 * texture.clay_g_per_megagram) * soil_mass_megagrams;
    const acid_insensitive_mol = 0.167 * surface_area_m2 * 5.0e3;
    if (!std.math.isFinite(surface_area_m2) or !std.math.isFinite(acid_insensitive_mol))
        return error.InvalidSilicateInventory;
    const base_silicate_mol = if (ph > 4.5) acid_insensitive_mol else 0.0;
    return .{
        .mineral_surface_area_m2 = surface_area_m2,
        .soil_silicates = .{
            .aluminum_mol = acid_insensitive_mol,
            .iron_mol = acid_insensitive_mol,
            .calcium_mol = base_silicate_mol,
            .magnesium_mol = base_silicate_mol,
            .sodium_mol = base_silicate_mol,
            .potassium_mol = base_silicate_mol,
        },
        .ground_rock_silicates = zeroElements(),
    };
}

test "STARTE silicate initialization preserves texture operation order and alkaline branch" {
    const result = (try initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, .{
        .sand_g_per_megagram = 10,
        .silt_g_per_megagram = 20,
        .clay_g_per_megagram = 30,
    }, 2, 4.6)).?;
    const expected_area = (0.3 * 10.0 + 2.2 * 20.0 + 8.0 * 30.0) * 2.0;
    const expected_mol = 0.167 * expected_area * 5.0e3;
    try std.testing.expectEqual(expected_area, result.mineral_surface_area_m2);
    try std.testing.expectEqual(expected_mol, result.soil_silicates.aluminum_mol);
    try std.testing.expectEqual(expected_mol, result.soil_silicates.potassium_mol);
    try std.testing.expectEqual(@as(f64, 0), result.ground_rock_silicates.aluminum_mol);
}

test "STARTE pH threshold is strict and suppresses four base silicates" {
    const result = (try initialize(.{ .profile_source = .atmospheric_equilibrium, .gas_initialization_index = 0 }, .{ .sand_g_per_megagram = 1, .silt_g_per_megagram = 1, .clay_g_per_megagram = 1 }, 1, 4.5)).?;
    try std.testing.expect(result.soil_silicates.aluminum_mol > 0);
    try std.testing.expect(result.soil_silicates.iron_mol > 0);
    try std.testing.expectEqual(@as(f64, 0), result.soil_silicates.calcium_mol);
    try std.testing.expectEqual(@as(f64, 0), result.soil_silicates.magnesium_mol);
    try std.testing.expectEqual(@as(f64, 0), result.soil_silicates.sodium_mol);
    try std.testing.expectEqual(@as(f64, 0), result.soil_silicates.potassium_mol);
}

test "STARTE inactive silicate initialization ignores invalid dormant input" {
    const nan = std.math.nan(f64);
    try std.testing.expectEqual(@as(?Inventories, null), try initialize(.{ .profile_source = .supplied_profile, .gas_initialization_index = 0 }, .{ .sand_g_per_megagram = nan, .silt_g_per_megagram = nan, .clay_g_per_megagram = nan }, nan, nan));
}
