const std = @import("std");

pub const ReferenceDiffusivities = struct {
    carbon_dioxide_gas_m2_h: f64,
    methane_gas_m2_h: f64,
    oxygen_gas_m2_h: f64,
    nitrogen_gas_m2_h: f64,
    nitrous_oxide_gas_m2_h: f64,
    ammonia_gas_m2_h: f64,
    hydrogen_gas_m2_h: f64,
    methane_aqueous_m2_h: f64,
    oxygen_aqueous_m2_h: f64,
    nitrogen_aqueous_m2_h: f64,
    ammonia_aqueous_m2_h: f64,
    hydrogen_aqueous_m2_h: f64,
    nitrate_aqueous_m2_h: f64,
    nitrous_oxide_aqueous_m2_h: f64,
    dihydrogen_phosphate_aqueous_m2_h: f64,
    dissolved_organic_carbon_m2_h: f64,
    dissolved_organic_nitrogen_m2_h: f64,
    dissolved_organic_phosphorus_m2_h: f64,
    acetate_aqueous_m2_h: f64,
};

pub const AdjustedDiffusivities = ReferenceDiffusivities;

pub const Result = struct {
    gaseous_temperature_factor: f64,
    aqueous_temperature_factor: f64,
    nitrogen_diffusion_temperature_factor: f64,
    diffusivities: AdjustedDiffusivities,
};

pub const CalculationError = error{
    NonFiniteInput,
    NonPositiveTemperature,
    NegativeReferenceDiffusivity,
    NonFiniteResult,
};

/// Translates `hour1.f` lines 4628--4649 for one surface-residue layer.
pub fn calculate(
    residue_temperature_k: f64,
    reference: ReferenceDiffusivities,
) CalculationError!Result {
    if (!std.math.isFinite(residue_temperature_k)) return error.NonFiniteInput;
    if (residue_temperature_k <= 0.0) return error.NonPositiveTemperature;
    inline for (std.meta.fields(ReferenceDiffusivities)) |field| {
        const value = @field(reference, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (value < 0.0) return error.NegativeReferenceDiffusivity;
    }

    const relative_temperature = residue_temperature_k / 298.15;
    const gaseous_temperature_factor = std.math.pow(f64, relative_temperature, 1.75);
    const aqueous_temperature_factor = std.math.pow(f64, relative_temperature, 6.0);

    // Preserve the assignment order in HOUR1: nitrogen diffusion receives
    // the aqueous factor before individual gas and solute diffusivities.
    const nitrogen_diffusion_temperature_factor = aqueous_temperature_factor;
    const diffusivities = AdjustedDiffusivities{
        .carbon_dioxide_gas_m2_h = reference.carbon_dioxide_gas_m2_h *
            gaseous_temperature_factor,
        .methane_gas_m2_h = reference.methane_gas_m2_h * gaseous_temperature_factor,
        .oxygen_gas_m2_h = reference.oxygen_gas_m2_h * gaseous_temperature_factor,
        .nitrogen_gas_m2_h = reference.nitrogen_gas_m2_h * gaseous_temperature_factor,
        .nitrous_oxide_gas_m2_h = reference.nitrous_oxide_gas_m2_h *
            gaseous_temperature_factor,
        .ammonia_gas_m2_h = reference.ammonia_gas_m2_h * gaseous_temperature_factor,
        .hydrogen_gas_m2_h = reference.hydrogen_gas_m2_h * gaseous_temperature_factor,
        .methane_aqueous_m2_h = reference.methane_aqueous_m2_h *
            aqueous_temperature_factor,
        .oxygen_aqueous_m2_h = reference.oxygen_aqueous_m2_h *
            aqueous_temperature_factor,
        .nitrogen_aqueous_m2_h = reference.nitrogen_aqueous_m2_h *
            aqueous_temperature_factor,
        .ammonia_aqueous_m2_h = reference.ammonia_aqueous_m2_h *
            aqueous_temperature_factor,
        .hydrogen_aqueous_m2_h = reference.hydrogen_aqueous_m2_h *
            aqueous_temperature_factor,
        .nitrate_aqueous_m2_h = reference.nitrate_aqueous_m2_h *
            aqueous_temperature_factor,
        .nitrous_oxide_aqueous_m2_h = reference.nitrous_oxide_aqueous_m2_h *
            aqueous_temperature_factor,
        .dihydrogen_phosphate_aqueous_m2_h = reference.dihydrogen_phosphate_aqueous_m2_h *
            aqueous_temperature_factor,
        .dissolved_organic_carbon_m2_h = reference.dissolved_organic_carbon_m2_h *
            aqueous_temperature_factor,
        .dissolved_organic_nitrogen_m2_h = reference.dissolved_organic_nitrogen_m2_h *
            aqueous_temperature_factor,
        .dissolved_organic_phosphorus_m2_h = reference.dissolved_organic_phosphorus_m2_h *
            aqueous_temperature_factor,
        .acetate_aqueous_m2_h = reference.acetate_aqueous_m2_h *
            aqueous_temperature_factor,
    };

    if (!std.math.isFinite(gaseous_temperature_factor) or
        !std.math.isFinite(aqueous_temperature_factor))
    {
        return error.NonFiniteResult;
    }
    inline for (std.meta.fields(AdjustedDiffusivities)) |field| {
        if (!std.math.isFinite(@field(diffusivities, field.name))) {
            return error.NonFiniteResult;
        }
    }
    return .{
        .gaseous_temperature_factor = gaseous_temperature_factor,
        .aqueous_temperature_factor = aqueous_temperature_factor,
        .nitrogen_diffusion_temperature_factor = nitrogen_diffusion_temperature_factor,
        .diffusivities = diffusivities,
    };
}

fn testReferenceDiffusivities() ReferenceDiffusivities {
    return .{
        .carbon_dioxide_gas_m2_h = 4.68e-2,
        .methane_gas_m2_h = 7.80e-2,
        .oxygen_gas_m2_h = 6.43e-2,
        .nitrogen_gas_m2_h = 5.57e-2,
        .nitrous_oxide_gas_m2_h = 5.57e-2,
        .ammonia_gas_m2_h = 6.67e-2,
        .hydrogen_gas_m2_h = 5.57e-2,
        .methane_aqueous_m2_h = 7.08e-6,
        .oxygen_aqueous_m2_h = 8.57e-6,
        .nitrogen_aqueous_m2_h = 7.34e-6,
        .ammonia_aqueous_m2_h = 4.00e-6,
        .hydrogen_aqueous_m2_h = 7.34e-6,
        .nitrate_aqueous_m2_h = 6.00e-6,
        .nitrous_oxide_aqueous_m2_h = 5.72e-6,
        .dihydrogen_phosphate_aqueous_m2_h = 3.00e-6,
        .dissolved_organic_carbon_m2_h = 1.00e-8,
        .dissolved_organic_nitrogen_m2_h = 1.00e-8,
        .dissolved_organic_phosphorus_m2_h = 1.00e-8,
        .acetate_aqueous_m2_h = 3.64e-6,
    };
}

test "reference temperature preserves every residue diffusivity" {
    const reference = testReferenceDiffusivities();
    const result = try calculate(298.15, reference);

    try std.testing.expectEqual(@as(f64, 1.0), result.gaseous_temperature_factor);
    try std.testing.expectEqual(@as(f64, 1.0), result.aqueous_temperature_factor);
    try std.testing.expectEqual(reference, result.diffusivities);
}

test "gas and aqueous diffusivities use their distinct temperature exponents" {
    const reference = testReferenceDiffusivities();
    const temperature_k = 280.0;
    const result = try calculate(temperature_k, reference);
    const expected_gas_factor = std.math.pow(f64, temperature_k / 298.15, 1.75);
    const expected_aqueous_factor = std.math.pow(f64, temperature_k / 298.15, 6.0);

    try std.testing.expectApproxEqRel(
        reference.oxygen_gas_m2_h * expected_gas_factor,
        result.diffusivities.oxygen_gas_m2_h,
        1.0e-14,
    );
    try std.testing.expectApproxEqRel(
        reference.nitrate_aqueous_m2_h * expected_aqueous_factor,
        result.diffusivities.nitrate_aqueous_m2_h,
        1.0e-14,
    );
    try std.testing.expectEqual(
        result.aqueous_temperature_factor,
        result.nitrogen_diffusion_temperature_factor,
    );
}
