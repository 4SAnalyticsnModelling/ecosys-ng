const std = @import("std");

pub const Inputs = struct {
    organic_carbon_g: f64,
    organic_residue_g: f64,
    bulk_volume_m3: f64,
    minimum_bulk_volume_m3: f64,
    sand_mg: f64,
    silt_mg: f64,
    clay_mg: f64,
    bulk_density_mg_m3: f64,
    surface_litter_volume_m3: f64,
    grid_cell_area_m2: f64,
    stalk_area_m2: f64,
    canopy_height_m: f64,
    mass_threshold_mg: f64,
    bulk_density_threshold_mg_m3: f64,
};

pub const Result = struct {
    adjusted_bulk_volume_m3: f64,
    organic_carbon_mass_fraction: f64,
    organic_residue_mass_fraction: f64,
    organic_carbon_concentration_g_mg: f64,
    sand_mass_fraction: f64,
    silt_mass_fraction: f64,
    clay_mass_fraction: f64,
    characteristic_particle_size_um: f64,
    particle_roughness_m: f64,
    soil_surface_roughness_m: f64,
    total_surface_roughness_m: f64,
};

pub const CalculationError = error{
    NonFiniteInput,
    NegativeInput,
    InvalidGridCellArea,
    InvalidCanopyHeight,
    NonFiniteResult,
};

/// Translates HOUR1 lines 2927-2967 in legacy operation order.
pub fn calculate(inputs: Inputs) CalculationError!Result {
    inline for (std.meta.fields(Inputs)) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (value < 0.0) return error.NegativeInput;
    }
    if (inputs.grid_cell_area_m2 <= 0.0) return error.InvalidGridCellArea;

    const organic_carbon_g = @max(0.0, inputs.organic_carbon_g);
    const organic_residue_g = @max(0.0, inputs.organic_residue_g);
    const adjusted_bulk_volume_m3 = @max(
        inputs.bulk_volume_m3,
        inputs.minimum_bulk_volume_m3 + 1.82e-6 * organic_carbon_g,
    );
    const normalized_mass_mg = inputs.sand_mg + inputs.silt_mg +
        inputs.clay_mg + 1.82e-6 * organic_carbon_g;

    var organic_carbon_fraction: f64 = 0.0;
    var organic_residue_fraction: f64 = 0.0;
    var organic_carbon_concentration_g_mg: f64 = 0.0;
    var sand_fraction: f64 = 0.0;
    var silt_fraction: f64 = 1.0;
    var clay_fraction: f64 = 0.0;
    if (normalized_mass_mg > inputs.mass_threshold_mg) {
        organic_carbon_fraction = 1.82e-6 * organic_carbon_g / normalized_mass_mg;
        organic_residue_fraction = 1.82e-6 * organic_residue_g / normalized_mass_mg;
        organic_carbon_concentration_g_mg = 0.55e6 * organic_carbon_fraction;
        sand_fraction = inputs.sand_mg / normalized_mass_mg;
        silt_fraction = inputs.silt_mg / normalized_mass_mg;
        clay_fraction = inputs.clay_mg / normalized_mass_mg;
    }

    const humified_organic_fraction = organic_carbon_fraction - organic_residue_fraction;
    var characteristic_particle_size_um: f64 = 0.0;
    var particle_roughness_m: f64 = 0.0;
    var soil_surface_roughness_m: f64 = 0.001;
    if (inputs.bulk_density_mg_m3 > inputs.bulk_density_threshold_mg_m3) {
        characteristic_particle_size_um = 1.0 * clay_fraction +
            10.0 * silt_fraction + 100.0 * sand_fraction +
            10.0 * humified_organic_fraction + 100.0 * organic_residue_fraction;
        particle_roughness_m = 0.041 *
            std.math.pow(f64, 1.0e-6 * characteristic_particle_size_um, 0.167);
        soil_surface_roughness_m = 0.01;
    }
    if (inputs.canopy_height_m == 0.0 and inputs.stalk_area_m2 > 0.0) {
        return error.InvalidCanopyHeight;
    }
    const total_surface_roughness_m = particle_roughness_m + soil_surface_roughness_m +
        1.0 * @min(
            soil_surface_roughness_m,
            inputs.surface_litter_volume_m3 / inputs.grid_cell_area_m2,
        ) +
        0.1 * inputs.stalk_area_m2 / inputs.grid_cell_area_m2 *
            soil_surface_roughness_m / @max(soil_surface_roughness_m, inputs.canopy_height_m);

    const result = Result{
        .adjusted_bulk_volume_m3 = adjusted_bulk_volume_m3,
        .organic_carbon_mass_fraction = organic_carbon_fraction,
        .organic_residue_mass_fraction = organic_residue_fraction,
        .organic_carbon_concentration_g_mg = organic_carbon_concentration_g_mg,
        .sand_mass_fraction = sand_fraction,
        .silt_mass_fraction = silt_fraction,
        .clay_mass_fraction = clay_fraction,
        .characteristic_particle_size_um = characteristic_particle_size_um,
        .particle_roughness_m = particle_roughness_m,
        .soil_surface_roughness_m = soil_surface_roughness_m,
        .total_surface_roughness_m = total_surface_roughness_m,
    };
    inline for (std.meta.fields(Result)) |field| {
        if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteResult;
    }
    return result;
}

test "disturbed mineral surface preserves composition and roughness equations" {
    const result = try calculate(.{
        .organic_carbon_g = 100_000.0,
        .organic_residue_g = 20_000.0,
        .bulk_volume_m3 = 2.0,
        .minimum_bulk_volume_m3 = 1.0,
        .sand_mg = 1.0,
        .silt_mg = 2.0,
        .clay_mg = 1.0,
        .bulk_density_mg_m3 = 1.0,
        .surface_litter_volume_m3 = 0.5,
        .grid_cell_area_m2 = 100.0,
        .stalk_area_m2 = 2.0,
        .canopy_height_m = 1.0,
        .mass_threshold_mg = 1.0e-12,
        .bulk_density_threshold_mg_m3 = 1.0e-12,
    });
    const normalized_mass = 1.0 + 2.0 + 1.0 + 1.82e-6 * 100_000.0;
    try std.testing.expectEqual(1.0 / normalized_mass, result.sand_mass_fraction);
    try std.testing.expect(result.total_surface_roughness_m > result.soil_surface_roughness_m);
}

test "empty surface uses legacy silt default and minimum roughness" {
    const result = try calculate(.{
        .organic_carbon_g = 0.0,
        .organic_residue_g = 0.0,
        .bulk_volume_m3 = 0.0,
        .minimum_bulk_volume_m3 = 0.0,
        .sand_mg = 0.0,
        .silt_mg = 0.0,
        .clay_mg = 0.0,
        .bulk_density_mg_m3 = 0.0,
        .surface_litter_volume_m3 = 0.0,
        .grid_cell_area_m2 = 1.0,
        .stalk_area_m2 = 0.0,
        .canopy_height_m = 0.0,
        .mass_threshold_mg = 0.0,
        .bulk_density_threshold_mg_m3 = 0.0,
    });
    try std.testing.expectEqual(@as(f64, 1.0), result.silt_mass_fraction);
    try std.testing.expectEqual(@as(f64, 0.001), result.total_surface_roughness_m);
}
