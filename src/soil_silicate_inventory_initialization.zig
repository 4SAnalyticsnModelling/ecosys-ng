const std = @import("std");

pub const Texture = struct {
    sand_g_per_Mg: f64,
    silt_g_per_Mg: f64,
    clay_g_per_Mg: f64,
};

pub const SurfaceAreaCoefficients = struct {
    sand_m2_per_g: f64,
    silt_m2_per_g: f64,
    clay_m2_per_g: f64,
};

pub const Parameters = struct {
    surface_area: SurfaceAreaCoefficients,
    elemental_inventory_fraction: f64,
    silicate_inventory_mol_per_m2: f64,
    base_cation_activation_ph: f64,
};

pub const ElementalSilicates = struct {
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
};

pub const Result = struct {
    natural_silicate_surface_area_m2: f64,
    natural: ElementalSilicates,
    ground: ElementalSilicates,
};

/// Direct translation of STARTE lines 1642--1662. Natural Al and Fe
/// inventories are always initialized; Ca, Mg, Na, and K use the strict
/// source pH gate. Ground-rock inventories start at zero.
pub fn calculate(
    texture: Texture,
    soil_mass_Mg: f64,
    soil_ph: f64,
    parameters: Parameters,
) !Result {
    inline for (@typeInfo(Texture).@"struct".fields) |field| {
        const value = @field(texture, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSilicateInventoryInput;
        if (value < 0) return error.InvalidSilicateInventoryInput;
    }
    inline for (@typeInfo(SurfaceAreaCoefficients).@"struct".fields) |field| {
        const value = @field(parameters.surface_area, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSilicateInventoryParameter;
        if (value < 0) return error.InvalidSilicateInventoryParameter;
    }
    inline for (.{
        soil_mass_Mg,
        soil_ph,
        parameters.elemental_inventory_fraction,
        parameters.silicate_inventory_mol_per_m2,
        parameters.base_cation_activation_ph,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteSilicateInventoryInput;
    if (soil_mass_Mg < 0 or
        parameters.elemental_inventory_fraction < 0 or
        parameters.silicate_inventory_mol_per_m2 < 0)
        return error.InvalidSilicateInventoryInput;

    const surface_area_m2 =
        (parameters.surface_area.sand_m2_per_g * texture.sand_g_per_Mg +
            parameters.surface_area.silt_m2_per_g * texture.silt_g_per_Mg +
            parameters.surface_area.clay_m2_per_g * texture.clay_g_per_Mg) *
        soil_mass_Mg;
    const elemental_inventory_mol =
        parameters.elemental_inventory_fraction *
        surface_area_m2 *
        parameters.silicate_inventory_mol_per_m2;
    const base_cation_inventory_mol =
        if (soil_ph > parameters.base_cation_activation_ph)
            elemental_inventory_mol
        else
            0;
    const result: Result = .{
        .natural_silicate_surface_area_m2 = surface_area_m2,
        .natural = .{
            .aluminum_mol = elemental_inventory_mol,
            .iron_mol = elemental_inventory_mol,
            .calcium_mol = base_cation_inventory_mol,
            .magnesium_mol = base_cation_inventory_mol,
            .sodium_mol = base_cation_inventory_mol,
            .potassium_mol = base_cation_inventory_mol,
        },
        .ground = std.mem.zeroes(ElementalSilicates),
    };
    if (!std.math.isFinite(result.natural_silicate_surface_area_m2))
        return error.NonFiniteSilicateInventoryResult;
    inline for (.{ result.natural, result.ground }) |silicates|
        inline for (@typeInfo(ElementalSilicates).@"struct".fields) |field|
            if (!std.math.isFinite(@field(silicates, field.name)))
                return error.NonFiniteSilicateInventoryResult;
    return result;
}

fn sourceParameters() Parameters {
    return .{
        .surface_area = .{
            .sand_m2_per_g = 0.3,
            .silt_m2_per_g = 2.2,
            .clay_m2_per_g = 8.0,
        },
        .elemental_inventory_fraction = 0.167,
        .silicate_inventory_mol_per_m2 = 5.0e3,
        .base_cation_activation_ph = 4.5,
    };
}

test "STARTE silicate inventories preserve texture and pH source order" {
    const result = try calculate(.{
        .sand_g_per_Mg = 10,
        .silt_g_per_Mg = 20,
        .clay_g_per_Mg = 30,
    }, 2, 6, sourceParameters());
    const expected_area = (0.3 * 10 + 2.2 * 20 + 8.0 * 30) * 2;
    const expected_inventory = 0.167 * expected_area * 5.0e3;
    try std.testing.expectEqual(
        expected_area,
        result.natural_silicate_surface_area_m2,
    );
    inline for (@typeInfo(ElementalSilicates).@"struct".fields) |field| {
        try std.testing.expectEqual(
            expected_inventory,
            @field(result.natural, field.name),
        );
        try std.testing.expectEqual(@as(f64, 0), @field(result.ground, field.name));
    }
}

test "STARTE silicate base-cation gate is strict at threshold pH" {
    const result = try calculate(.{
        .sand_g_per_Mg = 1,
        .silt_g_per_Mg = 1,
        .clay_g_per_Mg = 1,
    }, 1, 4.5, sourceParameters());
    try std.testing.expect(result.natural.aluminum_mol > 0);
    try std.testing.expect(result.natural.iron_mol > 0);
    try std.testing.expectEqual(@as(f64, 0), result.natural.calcium_mol);
    try std.testing.expectEqual(@as(f64, 0), result.natural.potassium_mol);
}

test "STARTE silicate initialization rejects invalid texture" {
    try std.testing.expectError(
        error.InvalidSilicateInventoryInput,
        calculate(.{
            .sand_g_per_Mg = -1,
            .silt_g_per_Mg = 1,
            .clay_g_per_Mg = 1,
        }, 1, 7, sourceParameters()),
    );
}
