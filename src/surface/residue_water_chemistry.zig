const std = @import("std");

pub const ActivityCoefficients = struct {
    monovalent: f64,
    divalent: f64,
    trivalent: f64,
};

pub const NutrientMasses = struct {
    ammonium_g_n: f64,
    ammonia_g_n: f64,
    nitrate_g_n: f64,
    nitrite_g_n: f64,
    dihydrogen_phosphate_g_p: f64,
    hydrogen_phosphate_g_p: f64,
};

pub const NutrientConcentrations = struct {
    ammonium_g_n_m3: f64,
    ammonia_g_n_m3: f64,
    nitrate_g_n_m3: f64,
    nitrite_g_n_m3: f64,
    dihydrogen_phosphate_g_p_m3: f64,
    hydrogen_phosphate_g_p_m3: f64,
};

pub const WaterPotentials = struct {
    matric_megapascal: f64,
    osmotic_megapascal: f64,
    gravitational_megapascal: f64,
    total_megapascal: f64,
};

pub const Inputs = struct {
    dry_residue_volume_m3: f64,
    residue_water_volume_m3: f64,
    volume_threshold_m3: f64,
    improved_matric_potential_megapascal: f64,
    residue_temperature_k: f64,
    total_ion_activity_mol_m3: f64,
    surface_elevation_m: f64,
    reference_depth_below_surface_m: f64,
    residue_thickness_m: f64,
    wet_activity_coefficients: ActivityCoefficients,
    nutrient_masses: NutrientMasses,
    underlying_soil_potentials: WaterPotentials,
};

pub const Result = struct {
    potentials: WaterPotentials,
    activity_coefficients: ActivityCoefficients,
    nutrients: NutrientConcentrations,
};

pub const CalculationError = error{
    NonFiniteInput,
    NegativeVolume,
    InvalidMatricPotential,
    InvalidTemperature,
    NegativeIonActivity,
    InvalidActivityCoefficient,
    NonFiniteResult,
};

/// Translates `hour1.f` lines 4488--4525 while replacing the preceding residue
/// PSISM calculation with `improved_matric_potential_megapascal`.
pub fn calculate(inputs: Inputs) CalculationError!Result {
    inline for (std.meta.fields(Inputs)) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) {
            return error.NonFiniteInput;
        }
    }
    const volumes = [_]f64{
        inputs.dry_residue_volume_m3,
        inputs.residue_water_volume_m3,
        inputs.volume_threshold_m3,
        inputs.residue_thickness_m,
        inputs.reference_depth_below_surface_m,
    };
    for (volumes) |volume| if (volume < 0.0) return error.NegativeVolume;
    if (inputs.improved_matric_potential_megapascal > 0.0) return error.InvalidMatricPotential;
    if (inputs.residue_temperature_k <= 0.0) return error.InvalidTemperature;
    if (inputs.total_ion_activity_mol_m3 < 0.0) return error.NegativeIonActivity;
    inline for (std.meta.fields(ActivityCoefficients)) |field| {
        const value = @field(inputs.wet_activity_coefficients, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (value < 0.0 or value > 1.0) return error.InvalidActivityCoefficient;
    }
    inline for (std.meta.fields(NutrientMasses)) |field| {
        if (!std.math.isFinite(@field(inputs.nutrient_masses, field.name))) {
            return error.NonFiniteInput;
        }
    }
    inline for (std.meta.fields(WaterPotentials)) |field| {
        if (!std.math.isFinite(@field(inputs.underlying_soil_potentials, field.name))) {
            return error.NonFiniteInput;
        }
    }

    const has_wet_residue =
        inputs.dry_residue_volume_m3 > inputs.volume_threshold_m3 and
        inputs.residue_water_volume_m3 > inputs.volume_threshold_m3;
    if (!has_wet_residue) {
        return .{
            .potentials = inputs.underlying_soil_potentials,
            .activity_coefficients = .{
                .monovalent = 1.0,
                .divalent = 1.0,
                .trivalent = 1.0,
            },
            .nutrients = std.mem.zeroes(NutrientConcentrations),
        };
    }

    const osmotic_megapascal =
        -8.3143e-6 * inputs.residue_temperature_k * inputs.total_ion_activity_mol_m3;
    const residue_midpoint_elevation_m =
        inputs.surface_elevation_m -
        inputs.reference_depth_below_surface_m +
        0.5 * inputs.residue_thickness_m;
    const gravitational_megapascal = 0.0098 * residue_midpoint_elevation_m;
    const potentials = WaterPotentials{
        .matric_megapascal = inputs.improved_matric_potential_megapascal,
        .osmotic_megapascal = osmotic_megapascal,
        .gravitational_megapascal = gravitational_megapascal,
        .total_megapascal = @min(
            0.0,
            inputs.improved_matric_potential_megapascal + osmotic_megapascal + gravitational_megapascal,
        ),
    };
    const nutrients = NutrientConcentrations{
        .ammonium_g_n_m3 = @max(0.0, inputs.nutrient_masses.ammonium_g_n / inputs.residue_water_volume_m3),
        .ammonia_g_n_m3 = @max(0.0, inputs.nutrient_masses.ammonia_g_n / inputs.residue_water_volume_m3),
        .nitrate_g_n_m3 = @max(0.0, inputs.nutrient_masses.nitrate_g_n / inputs.residue_water_volume_m3),
        .nitrite_g_n_m3 = @max(0.0, inputs.nutrient_masses.nitrite_g_n / inputs.residue_water_volume_m3),
        .dihydrogen_phosphate_g_p_m3 = @max(
            0.0,
            inputs.nutrient_masses.dihydrogen_phosphate_g_p /
                inputs.residue_water_volume_m3,
        ),
        .hydrogen_phosphate_g_p_m3 = @max(
            0.0,
            inputs.nutrient_masses.hydrogen_phosphate_g_p /
                inputs.residue_water_volume_m3,
        ),
    };
    inline for (std.meta.fields(WaterPotentials)) |field| {
        if (!std.math.isFinite(@field(potentials, field.name))) return error.NonFiniteResult;
    }
    inline for (std.meta.fields(NutrientConcentrations)) |field| {
        if (!std.math.isFinite(@field(nutrients, field.name))) return error.NonFiniteResult;
    }
    return .{
        .potentials = potentials,
        .activity_coefficients = inputs.wet_activity_coefficients,
        .nutrients = nutrients,
    };
}

test "wet residue calculates elevation-aware potentials and nutrient concentrations" {
    const result = try calculate(.{
        .dry_residue_volume_m3 = 1.0,
        .residue_water_volume_m3 = 2.0,
        .volume_threshold_m3 = 0.0,
        .improved_matric_potential_megapascal = -0.2,
        .residue_temperature_k = 280.0,
        .total_ion_activity_mol_m3 = 1.0,
        .surface_elevation_m = 100.0,
        .reference_depth_below_surface_m = 1.0,
        .residue_thickness_m = 0.2,
        .wet_activity_coefficients = .{
            .monovalent = 0.9,
            .divalent = 0.8,
            .trivalent = 0.7,
        },
        .nutrient_masses = .{
            .ammonium_g_n = 4.0,
            .ammonia_g_n = 2.0,
            .nitrate_g_n = 6.0,
            .nitrite_g_n = 1.0,
            .dihydrogen_phosphate_g_p = 3.0,
            .hydrogen_phosphate_g_p = 5.0,
        },
        .underlying_soil_potentials = std.mem.zeroes(WaterPotentials),
    });
    try std.testing.expectEqual(@as(f64, 2.0), result.nutrients.ammonium_g_n_m3);
    try std.testing.expect(result.potentials.gravitational_megapascal > 0.0);
    try std.testing.expectEqual(@as(f64, 0.9), result.activity_coefficients.monovalent);
}

test "dry residue inherits soil potentials and zeros chemistry" {
    const inherited = WaterPotentials{
        .matric_megapascal = -0.3,
        .osmotic_megapascal = -0.1,
        .gravitational_megapascal = 0.05,
        .total_megapascal = -0.35,
    };
    const result = try calculate(.{
        .dry_residue_volume_m3 = 0.0,
        .residue_water_volume_m3 = 0.0,
        .volume_threshold_m3 = 0.0,
        .improved_matric_potential_megapascal = 0.0,
        .residue_temperature_k = 273.15,
        .total_ion_activity_mol_m3 = 0.0,
        .surface_elevation_m = 0.0,
        .reference_depth_below_surface_m = 0.0,
        .residue_thickness_m = 0.0,
        .wet_activity_coefficients = std.mem.zeroes(ActivityCoefficients),
        .nutrient_masses = std.mem.zeroes(NutrientMasses),
        .underlying_soil_potentials = inherited,
    });
    try std.testing.expectEqual(inherited, result.potentials);
    try std.testing.expectEqual(
        std.mem.zeroes(NutrientConcentrations),
        result.nutrients,
    );
    try std.testing.expectEqual(@as(f64, 1.0), result.activity_coefficients.trivalent);
}
