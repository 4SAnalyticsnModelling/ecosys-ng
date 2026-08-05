//! **A5 DISPOSITION: never production bound.**
//!
//! Legacy source: `ecosys_f77/hour1.f` lines 3688--3716. This module is an exact,
//! tested, source-order translation of that range and is imported only by
//! `src/root.zig`, so it is reachable by `zig build test` and unreachable
//! from `executeHourlyScience`. That is intentional and must stay that way.
//!
//! Classification: architecturally superseded. Production stores the same physics in a
//! different representation, deliberately, with the deviation recorded in
//! `docs/model_changes.md`. Binding this kernel would reintroduce the
//! formulation the project chose to leave behind, as a second writer.
//!
//! Superseded by: `soil_solver_properties` and `soil_water_heat_step`, both
//! verified production bound (`soil_solver_properties` at
//! `ecosys_ng.zig:7249` via `State.initMapped`,
//! `soil_water_heat_step.advanceMappedDeferred` at
//! `ecosys_ng.zig:4279`).
//!
//! **Correction, 2026-08-05.** This line originally also named
//! `mineral_layer_phase_initialization`. It should not have:
//! `tools/a5_owner_call_path_check.ps1` shows that module is referenced by
//! nothing in `src/` except `root.zig`, so it is itself unbound and cannot
//! supersede anything. Citing an unbound module as the production owner is the
//! exact error that makes a "do not bind" disposition suppress a real gap.
//! The two owners named above are sufficient and were re-verified.
//!
//! Those three carry the dual-domain micropore/macropore split that the
//! legacy single-domain porosity and particle-density parameters cannot express.
//! The legacy form computes one porosity per layer; production needs a
//! micropore capacity, a macropore capacity, and their independent phase
//! fractions, and derives the de Vries thermal properties from all of them.
//!
//! Do not bind this module, and do not delete it: the tests are a source-order
//! comparison oracle for the range above. Full argument and the census behind
//! it: `docs/traceability/hour1_2039_5200_binding_survey.md`.
const std = @import("std");

pub const LayerPosition = enum { surface_mineral, other };

pub const Inputs = struct {
    soil_mass_megagrams: f64,
    organic_carbon_concentration_g_per_megagram: f64,
    silt_concentration_g_per_megagram: f64,
    clay_concentration_g_per_megagram: f64,
    sand_concentration_g_per_megagram: f64,
    bulk_density_megagrams_m3: f64,
    micropore_volume_fraction: f64,
    rock_volume_fraction: f64,
    existing_porosity_m3_m3: f64,
    existing_solid_heat_capacity_megajoules_k: f64,
    effective_soil_volume_m3: f64,
    total_layer_volume_m3: f64,
    macropore_volume_fraction: f64,
    positive_mass_threshold_megagrams: f64,
    layer_position: LayerPosition,
};

pub const Result = struct {
    organic_matter_mass_fraction: f64,
    particle_density_megagrams_m3: f64,
    organic_volume_fraction_m3_m3: f64,
    mineral_volume_fraction_m3_m3: f64,
    sand_volume_fraction_m3_m3: f64,
    solid_heat_capacity_megajoules_k: f64,
    weighted_thermal_conductivity_m_megajoules_h_k: f64,
    thermal_conductivity_weight: f64,
    porosity_m3_m3: f64,
    micropore_capacity_m3: f64,
    macropore_capacity_m3: f64,
};

pub const CalculationError = error{
    NonFiniteInput,
    NegativeInput,
    InvalidFraction,
    InvalidPorosity,
    NonFiniteResult,
};

/// Translates HOUR1 lines 3688-3716 in source operation order. The positive
/// mass path deliberately retains the existing VHCM value because its legacy
/// assignment is commented out; the zero-mass path clears it.
pub fn calculate(inputs: Inputs) CalculationError!Result {
    inline for (std.meta.fields(Inputs)) |field| {
        if (field.type == f64) {
            const value = @field(inputs, field.name);
            if (!std.math.isFinite(value)) return error.NonFiniteInput;
            if (value < 0.0) return error.NegativeInput;
        }
    }
    if (inputs.micropore_volume_fraction > 1.0 or
        inputs.rock_volume_fraction > 1.0 or
        inputs.existing_porosity_m3_m3 > 1.0 or
        inputs.macropore_volume_fraction > 1.0)
    {
        return error.InvalidFraction;
    }

    var organic_matter_fraction: f64 = 0.0;
    var particle_density_megagrams_m3: f64 = 0.0;
    var organic_volume_fraction: f64 = 0.0;
    var mineral_volume_fraction: f64 = 0.0;
    var sand_volume_fraction: f64 = 0.0;
    var solid_heat_capacity_megajoules_k = inputs.existing_solid_heat_capacity_megajoules_k;
    var weighted_thermal_conductivity: f64 = 0.0;
    var thermal_conductivity_weight: f64 = 0.0;
    var porosity_m3_m3: f64 = 1.0;

    if (inputs.soil_mass_megagrams > inputs.positive_mass_threshold_megagrams) {
        organic_matter_fraction = @max(
            0.0,
            @min(1.0, 1.82e-6 * inputs.organic_carbon_concentration_g_per_megagram),
        );
        particle_density_megagrams_m3 =
            1.30 * organic_matter_fraction + 2.66 * (1.0 - organic_matter_fraction);
        const bulk_to_particle_ratio = inputs.bulk_density_megagrams_m3 / particle_density_megagrams_m3;
        const composition_total = organic_matter_fraction +
            inputs.silt_concentration_g_per_megagram +
            inputs.clay_concentration_g_per_megagram +
            inputs.sand_concentration_g_per_megagram;
        const normalization = if (composition_total == 0.0)
            1.0
        else
            @min(1.0, 1.0 / composition_total);
        organic_volume_fraction = organic_matter_fraction * normalization * bulk_to_particle_ratio;
        mineral_volume_fraction =
            (inputs.silt_concentration_g_per_megagram + inputs.clay_concentration_g_per_megagram) *
            normalization * bulk_to_particle_ratio;
        sand_volume_fraction =
            inputs.sand_concentration_g_per_megagram * normalization * bulk_to_particle_ratio;
        weighted_thermal_conductivity =
            (1.253 * organic_volume_fraction * 9.050e-4 +
                0.514 * mineral_volume_fraction * 1.056e-2 +
                0.386 * sand_volume_fraction * 2.112e-2) *
            inputs.micropore_volume_fraction +
            0.514 * inputs.rock_volume_fraction * 1.056e-2;
        thermal_conductivity_weight =
            (1.253 * organic_volume_fraction +
                0.514 * mineral_volume_fraction +
                0.386 * sand_volume_fraction) *
            inputs.micropore_volume_fraction +
            0.514 * inputs.rock_volume_fraction;
        const calculated_porosity = 1.0 - bulk_to_particle_ratio;
        porosity_m3_m3 = if (inputs.layer_position == .surface_mineral)
            @max(inputs.existing_porosity_m3_m3, calculated_porosity)
        else
            calculated_porosity;
        if (porosity_m3_m3 < 0.0 or porosity_m3_m3 > 1.0) {
            return error.InvalidPorosity;
        }
    } else {
        solid_heat_capacity_megajoules_k = 0.0;
    }

    const result = Result{
        .organic_matter_mass_fraction = organic_matter_fraction,
        .particle_density_megagrams_m3 = particle_density_megagrams_m3,
        .organic_volume_fraction_m3_m3 = organic_volume_fraction,
        .mineral_volume_fraction_m3_m3 = mineral_volume_fraction,
        .sand_volume_fraction_m3_m3 = sand_volume_fraction,
        .solid_heat_capacity_megajoules_k = solid_heat_capacity_megajoules_k,
        .weighted_thermal_conductivity_m_megajoules_h_k = weighted_thermal_conductivity,
        .thermal_conductivity_weight = thermal_conductivity_weight,
        .porosity_m3_m3 = porosity_m3_m3,
        .micropore_capacity_m3 = porosity_m3_m3 * inputs.effective_soil_volume_m3,
        .macropore_capacity_m3 = inputs.macropore_volume_fraction * inputs.total_layer_volume_m3,
    };
    inline for (std.meta.fields(Result)) |field| {
        if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteResult;
    }
    return result;
}

test "positive soil mass calculates thermal properties and retains VHCM" {
    const result = try calculate(.{
        .soil_mass_megagrams = 100.0,
        .organic_carbon_concentration_g_per_megagram = 100_000.0,
        .silt_concentration_g_per_megagram = 0.3,
        .clay_concentration_g_per_megagram = 0.2,
        .sand_concentration_g_per_megagram = 0.4,
        .bulk_density_megagrams_m3 = 1.2,
        .micropore_volume_fraction = 0.8,
        .rock_volume_fraction = 0.1,
        .existing_porosity_m3_m3 = 0.6,
        .existing_solid_heat_capacity_megajoules_k = 12.0,
        .effective_soil_volume_m3 = 80.0,
        .total_layer_volume_m3 = 100.0,
        .macropore_volume_fraction = 0.05,
        .positive_mass_threshold_megagrams = 0.0,
        .layer_position = .surface_mineral,
    });
    try std.testing.expectEqual(@as(f64, 12.0), result.solid_heat_capacity_megajoules_k);
    try std.testing.expectEqual(@as(f64, 0.6), result.porosity_m3_m3);
    try std.testing.expectEqual(@as(f64, 48.0), result.micropore_capacity_m3);
}

test "zero soil mass clears thermal state and uses unit porosity" {
    const result = try calculate(.{
        .soil_mass_megagrams = 0.0,
        .organic_carbon_concentration_g_per_megagram = 0.0,
        .silt_concentration_g_per_megagram = 0.0,
        .clay_concentration_g_per_megagram = 0.0,
        .sand_concentration_g_per_megagram = 0.0,
        .bulk_density_megagrams_m3 = 0.0,
        .micropore_volume_fraction = 0.0,
        .rock_volume_fraction = 0.0,
        .existing_porosity_m3_m3 = 0.0,
        .existing_solid_heat_capacity_megajoules_k = 9.0,
        .effective_soil_volume_m3 = 2.0,
        .total_layer_volume_m3 = 3.0,
        .macropore_volume_fraction = 0.0,
        .positive_mass_threshold_megagrams = 0.0,
        .layer_position = .other,
    });
    try std.testing.expectEqual(@as(f64, 0.0), result.solid_heat_capacity_megajoules_k);
    try std.testing.expectEqual(@as(f64, 1.0), result.porosity_m3_m3);
    try std.testing.expectEqual(@as(f64, 2.0), result.micropore_capacity_m3);
}
