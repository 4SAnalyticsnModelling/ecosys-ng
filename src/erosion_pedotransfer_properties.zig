const std = @import("std");

pub const DisturbanceMode = enum {
    no_effects,
    freeze_thaw,
    freeze_thaw_and_erosion,
    freeze_thaw_and_organic_matter_change,
    freeze_thaw_erosion_and_organic_matter_change,

    fn enabled(self: DisturbanceMode) bool {
        return self != .no_effects;
    }
};

pub const Inputs = struct {
    characteristic_particle_size_um: f64,
    silt_mass_fraction: f64,
    organic_residue_mass_fraction: f64,
    humified_organic_mass_fraction: f64,
    clay_mass_fraction: f64,
    organic_carbon_mass_fraction: f64,
    root_length_density_m_m3: f64,
    water_viscosity_megagrams_m_s: f64,
    surface_temperature_c: f64,
};

pub const Properties = struct {
    runoff_capacity_coefficient: f64,
    runoff_capacity_exponent: f64,
    rainfall_detachability_g_j: f64,
    soil_cohesion_kpa: f64,
    detachability_coefficient: f64,
    particle_density_megagrams_m3: f64,
    water_viscosity_megagrams_m_s: f64,
    sediment_settling_rate_m_h: f64,
};

pub const CalculationError = error{
    NonFiniteInput,
    InvalidParticleSize,
    InvalidMassFraction,
    InvalidRootLengthDensity,
    InvalidWaterViscosity,
    InvalidSurfaceTemperature,
    NonFiniteResult,
    InvalidResult,
};

/// Translates HOUR1 lines 2986-3011 in source operation order.
pub fn calculate(
    mode: DisturbanceMode,
    inputs: Inputs,
) CalculationError!?Properties {
    if (!mode.enabled()) return null;
    inline for (std.meta.fields(Inputs)) |field| {
        if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteInput;
    }
    if (inputs.characteristic_particle_size_um < 0.0) return error.InvalidParticleSize;
    const fractions = [_]f64{
        inputs.silt_mass_fraction,
        inputs.organic_residue_mass_fraction,
        inputs.humified_organic_mass_fraction,
        inputs.clay_mass_fraction,
        inputs.organic_carbon_mass_fraction,
    };
    for (fractions) |fraction| {
        if (fraction < 0.0 or fraction > 1.0) return error.InvalidMassFraction;
    }
    if (inputs.silt_mass_fraction + inputs.organic_residue_mass_fraction > 1.0) {
        return error.InvalidMassFraction;
    }
    if (inputs.root_length_density_m_m3 < 0.0) return error.InvalidRootLengthDensity;
    if (inputs.water_viscosity_megagrams_m_s <= 0.0) return error.InvalidWaterViscosity;
    if (inputs.surface_temperature_c <= -273.15) return error.InvalidSurfaceTemperature;

    const runoff_capacity_coefficient = std.math.pow(
        f64,
        (inputs.characteristic_particle_size_um + 5.0) / 0.32,
        -0.6,
    );
    const runoff_capacity_exponent = std.math.pow(
        f64,
        (inputs.characteristic_particle_size_um + 5.0) / 300.0,
        0.25,
    );
    const rainfall_detachability_g_j = 1.0e-6 *
        (1.5 + 2.5 *
            (1.0 - inputs.silt_mass_fraction - inputs.organic_residue_mass_fraction));
    const soil_cohesion_kpa = 2.0 +
        5.0 * inputs.organic_residue_mass_fraction +
        5.0 * inputs.humified_organic_mass_fraction +
        5.0 * inputs.clay_mass_fraction +
        5.0 * (1.0 - @exp(-1.0e-5 * inputs.root_length_density_m_m3));
    const detachability_coefficient = 0.79 * @exp(-0.85 * soil_cohesion_kpa);
    const particle_density_megagrams_m3 = 1.30 * inputs.organic_carbon_mass_fraction +
        2.66 * (1.0 - inputs.organic_carbon_mass_fraction);
    const water_viscosity_megagrams_m_s = inputs.water_viscosity_megagrams_m_s *
        @exp(0.533 - 0.0267 * inputs.surface_temperature_c);
    const sediment_settling_rate_m_h = 3.6e3 * 9.8 *
        (particle_density_megagrams_m3 - 1.0) *
        std.math.pow(f64, 1.0e-6 * inputs.characteristic_particle_size_um, 2.0) /
        (18.0 * water_viscosity_megagrams_m_s);

    const properties = Properties{
        .runoff_capacity_coefficient = runoff_capacity_coefficient,
        .runoff_capacity_exponent = runoff_capacity_exponent,
        .rainfall_detachability_g_j = rainfall_detachability_g_j,
        .soil_cohesion_kpa = soil_cohesion_kpa,
        .detachability_coefficient = detachability_coefficient,
        .particle_density_megagrams_m3 = particle_density_megagrams_m3,
        .water_viscosity_megagrams_m_s = water_viscosity_megagrams_m_s,
        .sediment_settling_rate_m_h = sediment_settling_rate_m_h,
    };
    inline for (std.meta.fields(Properties)) |field| {
        const value = @field(properties, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteResult;
        if (value < 0.0) return error.InvalidResult;
    }
    return properties;
}

test "enabled disturbance preserves erosion pedotransfer equations" {
    const properties = (try calculate(.freeze_thaw_and_erosion, .{
        .characteristic_particle_size_um = 50.0,
        .silt_mass_fraction = 0.3,
        .organic_residue_mass_fraction = 0.05,
        .humified_organic_mass_fraction = 0.1,
        .clay_mass_fraction = 0.2,
        .organic_carbon_mass_fraction = 0.15,
        .root_length_density_m_m3 = 10_000.0,
        .water_viscosity_megagrams_m_s = 1.0e-6,
        .surface_temperature_c = 20.0,
    })).?;
    const expected_detachability = 1.0e-6 * (1.5 + 2.5 * (1.0 - 0.3 - 0.05));
    try std.testing.expectApproxEqAbs(
        expected_detachability,
        properties.rainfall_detachability_g_j,
        1.0e-20,
    );
    try std.testing.expect(properties.sediment_settling_rate_m_h > 0.0);
}

test "no disturbance returns no recalculation" {
    const result = try calculate(.no_effects, undefined);
    try std.testing.expect(result == null);
}
