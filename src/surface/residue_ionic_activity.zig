//! **A5 DISPOSITION: never production bound.**
//!
//! Legacy source: `ecosys_f77/hour1.f` lines 4414--4454. This module is an exact,
//! tested, source-order translation of that range and is imported only by
//! `src/module_index.zig`, so it is reachable by `zig build test` and unreachable
//! from `executeHourlyScience`. That is intentional and must stay that way.
//!
//! Classification: architecturally superseded. Production stores the same physics in a
//! different representation, deliberately, with the deviation recorded in
//! `docs/model_changes.md`. Binding this kernel would reintroduce the
//! formulation the project chose to leave behind, as a second writer.
//!
//! Superseded by: the same `solute_activity_coefficients` / `solute_ion_activities` path.
//!
//! This is the surface-litter twin of `soil_ionic_strength_conductivity`
//! and is superseded for the same reason. EXEC-004 settled that litter chemistry
//! state is owned by the runtime litter SOLUTE state, not by a HOUR1 recompute.
//!
//! Do not bind this module, and do not delete it: the tests are a source-order
//! comparison oracle for the range above. Full argument and the census behind
//! it: `docs/traceability/hour1_2039_5200_binding_survey.md`.
const std = @import("std");

/// Species amounts follow `hour1.f` lines 4415--4431 source order.
pub const ChargeInputs = struct {
    trivalent_cations_mol: []const f64,
    trivalent_anions_mol: []const f64,
    divalent_cations_mol: []const f64,
    divalent_anions_mol: []const f64,
    divalent_anion_phosphorus_g_p: []const f64,
    monovalent_cation_nitrogen_g_n: []const f64,
    monovalent_cations_mol: []const f64,
    monovalent_anion_nitrogen_g_n: []const f64,
    monovalent_anions_mol: []const f64,
    monovalent_anion_phosphorus_g_p: []const f64,
    monovalent_calcium_phosphate_mol: []const f64,
    neutral_species_mol: []const f64,
};

pub const ActivityCoefficients = struct {
    monovalent: f64,
    divalent: f64,
    trivalent: f64,
};

pub const Result = struct {
    ionic_strength: f64,
    activity_expression: f64,
    activity_coefficients: ActivityCoefficients,
    total_ion_activity_mol_m3: f64,
};

pub const CalculationError = error{
    DimensionMismatch,
    NonFiniteInput,
    NegativeInput,
    InvalidWaterVolume,
    NonFiniteResult,
};

/// Translates `hour1.f` lines 4414--4454 for surface residue water.
pub fn calculate(
    water_volume_m3: f64,
    water_volume_threshold_m3: f64,
    inputs: ChargeInputs,
) CalculationError!Result {
    if (!std.math.isFinite(water_volume_m3) or
        !std.math.isFinite(water_volume_threshold_m3))
    {
        return error.NonFiniteInput;
    }
    if (water_volume_m3 < 0.0 or water_volume_threshold_m3 < 0.0) {
        return error.InvalidWaterVolume;
    }
    try validateDimensions(inputs);
    inline for (std.meta.fields(ChargeInputs)) |field| {
        for (@field(inputs, field.name)) |amount| {
            if (!std.math.isFinite(amount)) return error.NonFiniteInput;
            if (amount < 0.0) return error.NegativeInput;
        }
    }
    if (water_volume_m3 <= water_volume_threshold_m3) {
        return .{
            .ionic_strength = 0.0,
            .activity_expression = 0.0,
            .activity_coefficients = .{
                .monovalent = 1.0,
                .divalent = 1.0,
                .trivalent = 1.0,
            },
            .total_ion_activity_mol_m3 = 0.0,
        };
    }

    const cations3 = orderedSum(inputs.trivalent_cations_mol);
    const anions3 = orderedSum(inputs.trivalent_anions_mol);
    const cations2 = orderedSum(inputs.divalent_cations_mol);
    var anions2 = orderedSum(inputs.divalent_anions_mol);
    anions2 = anions2 + orderedSum(inputs.divalent_anion_phosphorus_g_p) / 31.0;
    var cations1 = orderedSum(inputs.monovalent_cation_nitrogen_g_n) / 14.0;
    cations1 = cations1 + orderedSum(inputs.monovalent_cations_mol);
    var anions1 = orderedSum(inputs.monovalent_anion_nitrogen_g_n) / 14.0;
    anions1 = anions1 + orderedSum(inputs.monovalent_anions_mol);
    anions1 = anions1 + orderedSum(inputs.monovalent_anion_phosphorus_g_p) / 31.0;
    anions1 = anions1 + orderedSum(inputs.monovalent_calcium_phosphate_mol);
    const neutral = orderedSum(inputs.neutral_species_mol);

    const ionic_strength = @max(
        0.0,
        0.5e-3 *
            (9.0 * (cations3 + anions3) +
                4.0 * (cations2 + anions2) +
                cations1 +
                anions1) /
            water_volume_m3,
    );
    const root_strength = @sqrt(ionic_strength);
    const activity_expression =
        root_strength / (1.0 + root_strength) - 0.20 * ionic_strength;
    const coefficients = ActivityCoefficients{
        .monovalent = @min(
            1.0,
            std.math.pow(f64, 10.0, -0.509 * 1.0 * activity_expression),
        ),
        .divalent = @min(
            1.0,
            std.math.pow(f64, 10.0, -0.509 * 4.0 * activity_expression),
        ),
        .trivalent = @min(
            1.0,
            std.math.pow(f64, 10.0, -0.509 * 9.0 * activity_expression),
        ),
    };
    const total_activity = @max(
        0.0,
        ((cations3 + anions3) * coefficients.trivalent +
            (cations2 + anions2) * coefficients.divalent +
            (cations1 + anions1) * coefficients.monovalent +
            neutral) /
            water_volume_m3,
    );
    const result = Result{
        .ionic_strength = ionic_strength,
        .activity_expression = activity_expression,
        .activity_coefficients = coefficients,
        .total_ion_activity_mol_m3 = total_activity,
    };
    if (!std.math.isFinite(ionic_strength) or
        !std.math.isFinite(activity_expression) or
        !std.math.isFinite(total_activity) or
        !std.math.isFinite(coefficients.monovalent) or
        !std.math.isFinite(coefficients.divalent) or
        !std.math.isFinite(coefficients.trivalent))
    {
        return error.NonFiniteResult;
    }
    return result;
}

fn validateDimensions(inputs: ChargeInputs) CalculationError!void {
    if (inputs.trivalent_cations_mol.len != 2 or
        inputs.trivalent_anions_mol.len != 1 or
        inputs.divalent_cations_mol.len != 5 or
        inputs.divalent_anions_mol.len != 2 or
        inputs.divalent_anion_phosphorus_g_p.len != 1 or
        inputs.monovalent_cation_nitrogen_g_n.len != 1 or
        inputs.monovalent_cations_mol.len != 13 or
        inputs.monovalent_anion_nitrogen_g_n.len != 1 or
        inputs.monovalent_anions_mol.len != 8 or
        inputs.monovalent_anion_phosphorus_g_p.len != 1 or
        inputs.monovalent_calcium_phosphate_mol.len != 1 or
        inputs.neutral_species_mol.len != 9)
    {
        return error.DimensionMismatch;
    }
}

fn orderedSum(values: []const f64) f64 {
    var total: f64 = 0.0;
    for (values) |value| total = total + value;
    return total;
}

test "surface residue monovalent ions determine ionic activity" {
    const z1 = [_]f64{0.0};
    const z2 = [_]f64{ 0.0, 0.0 };
    const z5 = [_]f64{0.0} ** 5;
    const z8 = [_]f64{0.0} ** 8;
    const z9 = [_]f64{0.0} ** 9;
    var c13 = [_]f64{0.0} ** 13;
    c13[0] = 1.0;
    var a8 = z8;
    a8[0] = 1.0;
    const result = try calculate(1.0, 0.0, .{
        .trivalent_cations_mol = &z2,
        .trivalent_anions_mol = &z1,
        .divalent_cations_mol = &z5,
        .divalent_anions_mol = &z2,
        .divalent_anion_phosphorus_g_p = &z1,
        .monovalent_cation_nitrogen_g_n = &z1,
        .monovalent_cations_mol = &c13,
        .monovalent_anion_nitrogen_g_n = &z1,
        .monovalent_anions_mol = &a8,
        .monovalent_anion_phosphorus_g_p = &z1,
        .monovalent_calcium_phosphate_mol = &z1,
        .neutral_species_mol = &z9,
    });
    try std.testing.expectEqual(@as(f64, 0.001), result.ionic_strength);
    try std.testing.expect(result.total_ion_activity_mol_m3 > 0.0);
}

test "dry residue uses unit activity coefficients" {
    const z1 = [_]f64{0.0};
    const z2 = [_]f64{ 0.0, 0.0 };
    const z5 = [_]f64{0.0} ** 5;
    const z8 = [_]f64{0.0} ** 8;
    const z9 = [_]f64{0.0} ** 9;
    const z13 = [_]f64{0.0} ** 13;
    const result = try calculate(0.0, 0.0, .{
        .trivalent_cations_mol = &z2,
        .trivalent_anions_mol = &z1,
        .divalent_cations_mol = &z5,
        .divalent_anions_mol = &z2,
        .divalent_anion_phosphorus_g_p = &z1,
        .monovalent_cation_nitrogen_g_n = &z1,
        .monovalent_cations_mol = &z13,
        .monovalent_anion_nitrogen_g_n = &z1,
        .monovalent_anions_mol = &z8,
        .monovalent_anion_phosphorus_g_p = &z1,
        .monovalent_calcium_phosphate_mol = &z1,
        .neutral_species_mol = &z9,
    });
    try std.testing.expectEqual(@as(f64, 1.0), result.activity_coefficients.trivalent);
    try std.testing.expectEqual(@as(f64, 0.0), result.total_ion_activity_mol_m3);
}
