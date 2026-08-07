const std = @import("std");

pub const AqueousConcentrations = struct {
    ammonium_mol_n_per_m3: f64,
    ammonia_mol_n_per_m3: f64,
    nitrate_mol_n_per_m3: f64,
    dihydrogen_phosphate_mol_p_per_m3: f64,
    hydrogen_phosphate_mol_p_per_m3: f64,
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
    sodium_mol_per_m3: f64,
    potassium_mol_per_m3: f64,
};

pub const Parameters = struct {
    nitrogen_molar_mass_g_per_mol: f64,
    water_ion_product_mol2_per_m6: f64,
};

pub const Inventories = struct {
    ammonium_g_n: f64,
    /// STARTE line 1949 omits the nitrogen molar-mass conversion used for NH4.
    ammonia_source_mol_n: f64,
    /// STARTE line 1950 omits the nitrogen molar-mass conversion used for NH4.
    nitrate_source_mol_n: f64,
    nitrite_g_n: f64,
    dihydrogen_phosphate_mol_p: f64,
    hydrogen_phosphate_mol_p: f64,
    band_ammonium_g_n: f64,
    band_ammonia_g_n: f64,
    band_nitrate_g_n: f64,
    band_nitrite_g_n: f64,
    band_dihydrogen_phosphate_mol_p: f64,
    band_hydrogen_phosphate_mol_p: f64,
    hydrogen_mol: f64,
    hydroxide_mol: f64,
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
};

/// Direct translation of `starte.f` lines 1948--1970 for one surface litter cell.
/// Concentrations are multiplied by runtime field-capacity water volume.
pub fn initialize(
    water_capacity_m3: f64,
    surface_ph: f64,
    concentrations: AqueousConcentrations,
    parameters: Parameters,
) !Inventories {
    inline for (.{
        water_capacity_m3,
        surface_ph,
        concentrations.ammonium_mol_n_per_m3,
        concentrations.ammonia_mol_n_per_m3,
        concentrations.nitrate_mol_n_per_m3,
        concentrations.dihydrogen_phosphate_mol_p_per_m3,
        concentrations.hydrogen_phosphate_mol_p_per_m3,
        concentrations.aluminum_mol_per_m3,
        concentrations.iron_mol_per_m3,
        concentrations.calcium_mol_per_m3,
        concentrations.magnesium_mol_per_m3,
        concentrations.sodium_mol_per_m3,
        concentrations.potassium_mol_per_m3,
        parameters.nitrogen_molar_mass_g_per_mol,
        parameters.water_ion_product_mol2_per_m6,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterSoluteInput;
    }
    if (water_capacity_m3 < 0 or
        surface_ph < 0 or surface_ph > 14 or
        parameters.nitrogen_molar_mass_g_per_mol <= 0 or
        parameters.water_ion_product_mol2_per_m6 <= 0)
        return error.InvalidSurfaceLitterSoluteInput;
    inline for (.{
        concentrations.ammonium_mol_n_per_m3,
        concentrations.ammonia_mol_n_per_m3,
        concentrations.nitrate_mol_n_per_m3,
        concentrations.dihydrogen_phosphate_mol_p_per_m3,
        concentrations.hydrogen_phosphate_mol_p_per_m3,
        concentrations.aluminum_mol_per_m3,
        concentrations.iron_mol_per_m3,
        concentrations.calcium_mol_per_m3,
        concentrations.magnesium_mol_per_m3,
        concentrations.sodium_mol_per_m3,
        concentrations.potassium_mol_per_m3,
    }) |value| {
        if (value < 0) return error.InvalidSurfaceLitterSoluteInput;
    }

    const hydrogen_mol_per_m3 = std.math.pow(f64, 10, -(surface_ph - 3));
    const hydroxide_mol_per_m3 =
        parameters.water_ion_product_mol2_per_m6 / hydrogen_mol_per_m3;
    const result: Inventories = .{
        .ammonium_g_n = concentrations.ammonium_mol_n_per_m3 *
            water_capacity_m3 * parameters.nitrogen_molar_mass_g_per_mol,
        .ammonia_source_mol_n = concentrations.ammonia_mol_n_per_m3 *
            water_capacity_m3,
        .nitrate_source_mol_n = concentrations.nitrate_mol_n_per_m3 *
            water_capacity_m3,
        .nitrite_g_n = 0,
        .dihydrogen_phosphate_mol_p = concentrations.dihydrogen_phosphate_mol_p_per_m3 * water_capacity_m3,
        .hydrogen_phosphate_mol_p = concentrations.hydrogen_phosphate_mol_p_per_m3 * water_capacity_m3,
        .band_ammonium_g_n = 0,
        .band_ammonia_g_n = 0,
        .band_nitrate_g_n = 0,
        .band_nitrite_g_n = 0,
        .band_dihydrogen_phosphate_mol_p = 0,
        .band_hydrogen_phosphate_mol_p = 0,
        .hydrogen_mol = hydrogen_mol_per_m3 * water_capacity_m3,
        .hydroxide_mol = hydroxide_mol_per_m3 * water_capacity_m3,
        .aluminum_mol = concentrations.aluminum_mol_per_m3 * water_capacity_m3,
        .iron_mol = concentrations.iron_mol_per_m3 * water_capacity_m3,
        .calcium_mol = concentrations.calcium_mol_per_m3 * water_capacity_m3,
        .magnesium_mol = concentrations.magnesium_mol_per_m3 * water_capacity_m3,
        .sodium_mol = concentrations.sodium_mol_per_m3 * water_capacity_m3,
        .potassium_mol = concentrations.potassium_mol_per_m3 * water_capacity_m3,
    };
    inline for (std.meta.fields(Inventories)) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceLitterSoluteResult;
    }
    return result;
}

fn exampleConcentrations() AqueousConcentrations {
    return .{
        .ammonium_mol_n_per_m3 = 3,
        .ammonia_mol_n_per_m3 = 4,
        .nitrate_mol_n_per_m3 = 5,
        .dihydrogen_phosphate_mol_p_per_m3 = 6,
        .hydrogen_phosphate_mol_p_per_m3 = 7,
        .aluminum_mol_per_m3 = 8,
        .iron_mol_per_m3 = 9,
        .calcium_mol_per_m3 = 10,
        .magnesium_mol_per_m3 = 11,
        .sodium_mol_per_m3 = 12,
        .potassium_mol_per_m3 = 13,
    };
}

test "STARTE surface primary solutes preserve source operation order and scaling" {
    const result = try initialize(
        2,
        6,
        exampleConcentrations(),
        .{
            .nitrogen_molar_mass_g_per_mol = 14,
            .water_ion_product_mol2_per_m6 = 1.0e-8,
        },
    );
    try std.testing.expectEqual(@as(f64, 84), result.ammonium_g_n);
    try std.testing.expectEqual(@as(f64, 8), result.ammonia_source_mol_n);
    try std.testing.expectEqual(@as(f64, 10), result.nitrate_source_mol_n);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0e-3), result.hydrogen_mol, 1.0e-18);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0e-5), result.hydroxide_mol, 1.0e-20);
    try std.testing.expectEqual(@as(f64, 26), result.potassium_mol);
}

test "STARTE surface band inventories begin at zero" {
    const result = try initialize(
        1,
        7,
        exampleConcentrations(),
        .{
            .nitrogen_molar_mass_g_per_mol = 14,
            .water_ion_product_mol2_per_m6 = 1.0e-8,
        },
    );
    try std.testing.expectEqual(@as(f64, 0), result.band_ammonium_g_n);
    try std.testing.expectEqual(@as(f64, 0), result.band_nitrite_g_n);
    try std.testing.expectEqual(@as(f64, 0), result.band_hydrogen_phosphate_mol_p);
}

test "STARTE surface solute initialization rejects invalid input" {
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterSoluteInput,
        initialize(
            1,
            std.math.nan(f64),
            exampleConcentrations(),
            .{
                .nitrogen_molar_mass_g_per_mol = 14,
                .water_ion_product_mol2_per_m6 = 1.0e-8,
            },
        ),
    );
}
