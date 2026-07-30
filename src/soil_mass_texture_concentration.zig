const std = @import("std");

pub const SaltMode = enum { static_equilibrium, dynamic_transport };

pub const Inputs = struct {
    bulk_density_mg_m3: f64,
    rock_excluded_volume_m3: f64,
    organic_carbon_g: f64,
    sand_mg: f64,
    silt_mg: f64,
    clay_mg: f64,
    soil_mass_threshold_mg: f64,
};

pub const Result = struct {
    soil_mass_mg: f64,
    organic_carbon_concentration_g_mg: f64,
    sand_mass_fraction_mg_mg: f64,
    silt_mass_fraction_mg_mg: f64,
    clay_mass_fraction_mg_mg: f64,
    dynamic_salt_surface_area_m2_m3: ?f64,
};

pub const CalculationError = error{
    NonFiniteInput,
    NegativeInput,
    InvalidTextureFraction,
    NonFiniteResult,
};

/// Translates HOUR1 lines 3656-3671 in source operation order.
pub fn calculate(
    salt_mode: SaltMode,
    inputs: Inputs,
) CalculationError!Result {
    inline for (std.meta.fields(Inputs)) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (value < 0.0) return error.NegativeInput;
    }

    const soil_mass_mg = inputs.bulk_density_mg_m3 * inputs.rock_excluded_volume_m3;
    var organic_carbon_concentration_g_mg: f64 = 0.0;
    var sand_mass_fraction_mg_mg: f64 = 0.0;
    var silt_mass_fraction_mg_mg: f64 = 0.0;
    var clay_mass_fraction_mg_mg: f64 = 0.0;
    if (soil_mass_mg > inputs.soil_mass_threshold_mg) {
        organic_carbon_concentration_g_mg =
            @min(0.55e6, inputs.organic_carbon_g / soil_mass_mg);
        sand_mass_fraction_mg_mg = inputs.sand_mg / soil_mass_mg;
        silt_mass_fraction_mg_mg = inputs.silt_mg / soil_mass_mg;
        clay_mass_fraction_mg_mg = inputs.clay_mg / soil_mass_mg;
        if (sand_mass_fraction_mg_mg > 1.0 or
            silt_mass_fraction_mg_mg > 1.0 or
            clay_mass_fraction_mg_mg > 1.0 or
            sand_mass_fraction_mg_mg +
                silt_mass_fraction_mg_mg +
                clay_mass_fraction_mg_mg > 1.0 + 1.0e-12)
        {
            return error.InvalidTextureFraction;
        }
    }

    const dynamic_salt_surface_area_m2_m3 = if (salt_mode == .dynamic_transport)
        (0.3 * sand_mass_fraction_mg_mg +
            2.2 * silt_mass_fraction_mg_mg +
            8.0 * clay_mass_fraction_mg_mg) * inputs.bulk_density_mg_m3
    else
        null;
    const result = Result{
        .soil_mass_mg = soil_mass_mg,
        .organic_carbon_concentration_g_mg = organic_carbon_concentration_g_mg,
        .sand_mass_fraction_mg_mg = sand_mass_fraction_mg_mg,
        .silt_mass_fraction_mg_mg = silt_mass_fraction_mg_mg,
        .clay_mass_fraction_mg_mg = clay_mass_fraction_mg_mg,
        .dynamic_salt_surface_area_m2_m3 = dynamic_salt_surface_area_m2_m3,
    };
    inline for (std.meta.fields(Result)) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(result, field.name))) {
            return error.NonFiniteResult;
        }
    }
    if (dynamic_salt_surface_area_m2_m3) |surface_area| {
        if (!std.math.isFinite(surface_area)) return error.NonFiniteResult;
    }
    return result;
}

test "positive soil mass produces concentrations and dynamic surface area" {
    const result = try calculate(.dynamic_transport, .{
        .bulk_density_mg_m3 = 1.25,
        .rock_excluded_volume_m3 = 100.0,
        .organic_carbon_g = 1_000_000.0,
        .sand_mg = 50.0,
        .silt_mg = 40.0,
        .clay_mg = 30.0,
        .soil_mass_threshold_mg = 1.0e-12,
    });
    try std.testing.expectEqual(@as(f64, 125.0), result.soil_mass_mg);
    try std.testing.expectEqual(@as(f64, 0.4), result.sand_mass_fraction_mg_mg);
    try std.testing.expect(result.dynamic_salt_surface_area_m2_m3.? > 0.0);
}

test "negligible mass uses zero concentrations and static salt has no surface area update" {
    const result = try calculate(.static_equilibrium, .{
        .bulk_density_mg_m3 = 0.0,
        .rock_excluded_volume_m3 = 0.0,
        .organic_carbon_g = 0.0,
        .sand_mg = 0.0,
        .silt_mg = 0.0,
        .clay_mg = 0.0,
        .soil_mass_threshold_mg = 0.0,
    });
    try std.testing.expectEqual(@as(f64, 0.0), result.organic_carbon_concentration_g_mg);
    try std.testing.expect(result.dynamic_salt_surface_area_m2_m3 == null);
}
