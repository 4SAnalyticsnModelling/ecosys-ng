const std = @import("std");

pub const Inputs = struct {
    matric_potential_mpa: f64,
    soil_temperature_k: f64,
    total_ion_activity_mol_m3: f64,
    surface_elevation_m: f64,
    layer_midpoint_depth_below_surface_m: f64,
};

pub const Result = struct {
    osmotic_potential_mpa: f64,
    gravitational_potential_mpa: f64,
    total_potential_mpa: f64,
};

pub const CalculationError = error{
    NonFiniteInput,
    InvalidMatricPotential,
    InvalidSoilTemperature,
    NegativeIonActivity,
    NegativeLayerDepth,
    NonFiniteResult,
};

/// Translates HOUR1 lines 4181-4184 while making the elevation reference
/// explicit. Matric potential is supplied by the Mualem-van Genuchten solver.
pub fn calculate(inputs: Inputs) CalculationError!Result {
    inline for (std.meta.fields(Inputs)) |field| {
        if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteInput;
    }
    if (inputs.matric_potential_mpa > 0.0) return error.InvalidMatricPotential;
    if (inputs.soil_temperature_k <= 0.0) return error.InvalidSoilTemperature;
    if (inputs.total_ion_activity_mol_m3 < 0.0) return error.NegativeIonActivity;
    if (inputs.layer_midpoint_depth_below_surface_m < 0.0) return error.NegativeLayerDepth;

    const osmotic_potential_mpa =
        -8.3143e-6 * inputs.soil_temperature_k * inputs.total_ion_activity_mol_m3;
    const layer_midpoint_elevation_m =
        inputs.surface_elevation_m - inputs.layer_midpoint_depth_below_surface_m;
    const gravitational_potential_mpa = 0.0098 * layer_midpoint_elevation_m;
    const total_potential_mpa = @min(
        0.0,
        inputs.matric_potential_mpa +
            osmotic_potential_mpa +
            gravitational_potential_mpa,
    );
    const result = Result{
        .osmotic_potential_mpa = osmotic_potential_mpa,
        .gravitational_potential_mpa = gravitational_potential_mpa,
        .total_potential_mpa = total_potential_mpa,
    };
    inline for (std.meta.fields(Result)) |field| {
        if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteResult;
    }
    return result;
}

test "water potential includes absolute layer elevation" {
    const result = try calculate(.{
        .matric_potential_mpa = -0.1,
        .soil_temperature_k = 280.0,
        .total_ion_activity_mol_m3 = 10.0,
        .surface_elevation_m = 100.0,
        .layer_midpoint_depth_below_surface_m = 2.0,
    });
    const expected_osmotic = -8.3143e-6 * 280.0 * 10.0;
    const expected_gravitational = 0.0098 * (100.0 - 2.0);
    try std.testing.expectEqual(expected_osmotic, result.osmotic_potential_mpa);
    try std.testing.expectApproxEqAbs(
        expected_gravitational,
        result.gravitational_potential_mpa,
        1.0e-15,
    );
    try std.testing.expectEqual(@as(f64, 0.0), result.total_potential_mpa);
}

test "negative absolute elevation lowers total potential" {
    const result = try calculate(.{
        .matric_potential_mpa = -0.1,
        .soil_temperature_k = 273.15,
        .total_ion_activity_mol_m3 = 0.0,
        .surface_elevation_m = -10.0,
        .layer_midpoint_depth_below_surface_m = 5.0,
    });
    try std.testing.expectEqual(@as(f64, -0.147), result.gravitational_potential_mpa);
    try std.testing.expectEqual(@as(f64, -0.247), result.total_potential_mpa);
}
