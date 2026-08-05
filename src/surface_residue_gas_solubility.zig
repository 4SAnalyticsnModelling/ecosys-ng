//! **A5 DISPOSITION: never production bound.**
//!
//! Legacy source: `ecosys_f77/hour1.f` lines 4598--4619. This module is an exact,
//! tested, source-order translation of that range and is imported only by
//! `src/root.zig`, so it is reachable by `zig build test` and unreachable
//! from `executeHourlyScience`. That is intentional and must stay that way.
//!
//! Classification: architecturally superseded. Production stores the same physics in a
//! different representation, deliberately, with the deviation recorded in
//! `docs/model_changes.md`. Binding this kernel would reintroduce the
//! formulation the project chose to leave behind, as a second writer.
//!
//! Superseded by: the same `gas_transport.surfaceSolubilityWaterToAir` owner.
//!
//! Surface-litter twin of `soil_water_gas_solubility`; identical argument.
//!
//! Do not bind this module, and do not delete it: the tests are a source-order
//! comparison oracle for the range above. Full argument and the census behind
//! it: `docs/traceability/hour1_2039_5200_binding_survey.md`.
const std = @import("std");

pub const ReferenceSolubilityRatios = struct {
    carbon_dioxide: f64,
    methane: f64,
    oxygen: f64,
    nitrogen: f64,
    nitrous_oxide: f64,
    ammonia: f64,
    hydrogen: f64,
};

pub const SolubilityRatios = ReferenceSolubilityRatios;

pub const AqueousGasMasses = struct {
    carbon_dioxide_g: f64,
    methane_g: f64,
    oxygen_g: f64,
    nitrogen_g: f64,
    nitrous_oxide_g: f64,
    hydrogen_g: f64,
};

pub const AqueousGasConcentrations = struct {
    carbon_dioxide_g_m3: f64,
    methane_g_m3: f64,
    oxygen_g_m3: f64,
    nitrogen_g_m3: f64,
    nitrous_oxide_g_m3: f64,
    hydrogen_g_m3: f64,
};

pub const Result = struct {
    solubility_ratios: SolubilityRatios,
    aqueous_concentrations: AqueousGasConcentrations,
};

pub const CalculationError = error{
    NonFiniteInput,
    NegativeReferenceSolubility,
    TemperatureBelowAbsoluteZero,
    NegativeWaterVolume,
    NonFiniteResult,
};

/// Translates HOUR1 lines 4598-4619 for one surface-residue layer.
pub fn calculate(
    reference: ReferenceSolubilityRatios,
    residue_temperature_c: f64,
    aqueous_masses: AqueousGasMasses,
    water_volume_m3: f64,
    water_volume_threshold_m3: f64,
) CalculationError!Result {
    inline for (std.meta.fields(ReferenceSolubilityRatios)) |field| {
        const value = @field(reference, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (value < 0.0) return error.NegativeReferenceSolubility;
    }
    inline for (std.meta.fields(AqueousGasMasses)) |field| {
        if (!std.math.isFinite(@field(aqueous_masses, field.name))) {
            return error.NonFiniteInput;
        }
    }
    if (!std.math.isFinite(residue_temperature_c) or
        !std.math.isFinite(water_volume_m3) or
        !std.math.isFinite(water_volume_threshold_m3))
    {
        return error.NonFiniteInput;
    }
    if (residue_temperature_c < -273.15) return error.TemperatureBelowAbsoluteZero;
    if (water_volume_m3 < 0.0 or water_volume_threshold_m3 < 0.0) {
        return error.NegativeWaterVolume;
    }

    // Preserve the legacy statement order: temperature-adjust every
    // gas-water partition ratio before evaluating aqueous concentrations.
    const solubility_ratios = SolubilityRatios{
        .carbon_dioxide = reference.carbon_dioxide *
            @exp(0.843 - 0.0281 * residue_temperature_c),
        .methane = reference.methane *
            @exp(0.597 - 0.0199 * residue_temperature_c),
        .oxygen = reference.oxygen *
            @exp(0.516 - 0.0172 * residue_temperature_c),
        .nitrogen = reference.nitrogen *
            @exp(0.456 - 0.0152 * residue_temperature_c),
        .nitrous_oxide = reference.nitrous_oxide *
            @exp(0.897 - 0.0299 * residue_temperature_c),
        .ammonia = reference.ammonia *
            @exp(0.513 - 0.0171 * residue_temperature_c),
        .hydrogen = reference.hydrogen *
            @exp(0.597 - 0.0199 * residue_temperature_c),
    };

    const aqueous_concentrations = if (water_volume_m3 > water_volume_threshold_m3)
        AqueousGasConcentrations{
            .carbon_dioxide_g_m3 = @max(0.0, aqueous_masses.carbon_dioxide_g / water_volume_m3),
            .methane_g_m3 = @max(0.0, aqueous_masses.methane_g / water_volume_m3),
            .oxygen_g_m3 = @max(0.0, aqueous_masses.oxygen_g / water_volume_m3),
            .nitrogen_g_m3 = @max(0.0, aqueous_masses.nitrogen_g / water_volume_m3),
            .nitrous_oxide_g_m3 = @max(0.0, aqueous_masses.nitrous_oxide_g / water_volume_m3),
            .hydrogen_g_m3 = @max(0.0, aqueous_masses.hydrogen_g / water_volume_m3),
        }
    else
        std.mem.zeroes(AqueousGasConcentrations);

    inline for (std.meta.fields(SolubilityRatios)) |field| {
        if (!std.math.isFinite(@field(solubility_ratios, field.name))) {
            return error.NonFiniteResult;
        }
    }
    inline for (std.meta.fields(AqueousGasConcentrations)) |field| {
        if (!std.math.isFinite(@field(aqueous_concentrations, field.name))) {
            return error.NonFiniteResult;
        }
    }
    return .{
        .solubility_ratios = solubility_ratios,
        .aqueous_concentrations = aqueous_concentrations,
    };
}

test "temperature-adjusted ratios and wet residue concentrations preserve HOUR1 equations" {
    const result = try calculate(.{
        .carbon_dioxide = 0.7391,
        .methane = 0.03156,
        .oxygen = 0.02925,
        .nitrogen = 0.015,
        .nitrous_oxide = 0.5,
        .ammonia = 1.0,
        .hydrogen = 0.02,
    }, 20.0, .{
        .carbon_dioxide_g = 4.0,
        .methane_g = 2.0,
        .oxygen_g = 8.0,
        .nitrogen_g = 10.0,
        .nitrous_oxide_g = 1.0,
        .hydrogen_g = -0.5,
    }, 2.0, 0.0);

    try std.testing.expectApproxEqRel(
        0.7391 * @exp(0.843 - 0.0281 * 20.0),
        result.solubility_ratios.carbon_dioxide,
        1.0e-14,
    );
    try std.testing.expectApproxEqRel(
        0.02 * @exp(0.597 - 0.0199 * 20.0),
        result.solubility_ratios.hydrogen,
        1.0e-14,
    );
    try std.testing.expectEqual(@as(f64, 2.0), result.aqueous_concentrations.carbon_dioxide_g_m3);
    try std.testing.expectEqual(@as(f64, 0.0), result.aqueous_concentrations.hydrogen_g_m3);
}

test "water volume at threshold produces zero aqueous concentrations" {
    const result = try calculate(
        .{
            .carbon_dioxide = 1.0,
            .methane = 1.0,
            .oxygen = 1.0,
            .nitrogen = 1.0,
            .nitrous_oxide = 1.0,
            .ammonia = 1.0,
            .hydrogen = 1.0,
        },
        0.0,
        std.mem.zeroes(AqueousGasMasses),
        0.0,
        0.0,
    );
    try std.testing.expectEqual(
        std.mem.zeroes(AqueousGasConcentrations),
        result.aqueous_concentrations,
    );
}
