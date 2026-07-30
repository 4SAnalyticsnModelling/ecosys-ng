const std = @import("std");

pub const Inputs = struct {
    retained_pre_calcite_calcium_concentration_mol_per_m3: f64,
    retained_pre_calcite_calcium_activity_mol_per_m3: f64,
    calcium_transformation_mol_per_m3_step: f64,
    sulfate_concentration_mol_s_per_m3: f64,
    sulfate_transformation_mol_s_per_m3_step: f64,
    gypsum_concentration_mol_mineral_per_m3: f64,
    divalent_activity_coefficient: f64,
    gypsum_solubility_product: f64,
    substrate_limit_fraction_per_step: f64,
    maximum_precipitation_mol_mineral_per_m3_step: f64,
};

pub const Result = struct {
    provisional_sulfate_concentration_mol_s_per_m3: f64,
    provisional_sulfate_activity_mol_s_per_m3: f64,
    substrate_limit_mol_mineral_per_m3_step: f64,
    equilibrium_calcium_activity_mol_per_m3: f64,
    precipitation_mol_mineral_per_m3_step: f64,
    calcium_transformation_mol_per_m3_step: f64,
    sulfate_transformation_mol_s_per_m3_step: f64,
    gypsum_concentration_mol_mineral_per_m3: f64,
};

/// Direct source-order translation of SOLUTE.F lines 5067--5075.
///
/// The source reuses Ca concentration/activity computed before the immediately
/// preceding calcite update. Those retained values are explicit inputs here.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);

    // SOLUTE.F 5067--5072.
    const provisional_sulfate =
        inputs.sulfate_concentration_mol_s_per_m3 +
        inputs.sulfate_transformation_mol_s_per_m3_step;
    if (!std.math.isFinite(provisional_sulfate))
        return error.NonFiniteSurfaceLitterGypsumCorrectionResult;
    if (provisional_sulfate <= 0)
        return error.NonPositiveSurfaceLitterProvisionalSulfate;

    const sulfate_activity =
        provisional_sulfate * inputs.divalent_activity_coefficient;
    const substrate_limit =
        inputs.substrate_limit_fraction_per_step *
        @min(
            inputs.retained_pre_calcite_calcium_concentration_mol_per_m3,
            provisional_sulfate,
        );
    const equilibrium_calcium_activity =
        inputs.gypsum_solubility_product / sulfate_activity;
    const maximum =
        inputs.maximum_precipitation_mol_mineral_per_m3_step;
    const precipitation = @max(
        -@max(0.0, inputs.gypsum_concentration_mol_mineral_per_m3),
        -maximum,
        @min(
            maximum,
            substrate_limit,
            (inputs.retained_pre_calcite_calcium_activity_mol_per_m3 -
                equilibrium_calcium_activity) /
                inputs.divalent_activity_coefficient,
        ),
    );

    // SOLUTE.F 5073--5075.
    const result: Result = .{
        .provisional_sulfate_concentration_mol_s_per_m3 = provisional_sulfate,
        .provisional_sulfate_activity_mol_s_per_m3 = sulfate_activity,
        .substrate_limit_mol_mineral_per_m3_step = substrate_limit,
        .equilibrium_calcium_activity_mol_per_m3 = equilibrium_calcium_activity,
        .precipitation_mol_mineral_per_m3_step = precipitation,
        .calcium_transformation_mol_per_m3_step = inputs.calcium_transformation_mol_per_m3_step -
            precipitation,
        .sulfate_transformation_mol_s_per_m3_step = inputs.sulfate_transformation_mol_s_per_m3_step -
            precipitation,
        .gypsum_concentration_mol_mineral_per_m3 = inputs.gypsum_concentration_mol_mineral_per_m3 +
            precipitation,
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidSurfaceLitterGypsumCorrectionInput;
    }
    if (inputs.retained_pre_calcite_calcium_concentration_mol_per_m3 < 0 or
        inputs.retained_pre_calcite_calcium_activity_mol_per_m3 < 0 or
        inputs.sulfate_concentration_mol_s_per_m3 < 0 or
        inputs.gypsum_concentration_mol_mineral_per_m3 < 0 or
        inputs.divalent_activity_coefficient <= 0 or
        inputs.gypsum_solubility_product <= 0 or
        inputs.substrate_limit_fraction_per_step < 0 or
        inputs.substrate_limit_fraction_per_step > 1 or
        inputs.maximum_precipitation_mol_mineral_per_m3_step < 0)
    {
        return error.InvalidSurfaceLitterGypsumCorrectionInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceLitterGypsumCorrectionResult;
    }
    if (result.provisional_sulfate_concentration_mol_s_per_m3 <= 0 or
        result.provisional_sulfate_activity_mol_s_per_m3 <= 0 or
        result.substrate_limit_mol_mineral_per_m3_step < 0 or
        result.equilibrium_calcium_activity_mol_per_m3 < 0 or
        result.gypsum_concentration_mol_mineral_per_m3 < 0)
    {
        return error.InvalidSurfaceLitterGypsumCorrectionResult;
    }
}

fn testInputs() Inputs {
    return .{
        .retained_pre_calcite_calcium_concentration_mol_per_m3 = 0.6,
        .retained_pre_calcite_calcium_activity_mol_per_m3 = 0.3,
        .calcium_transformation_mol_per_m3_step = 0.05,
        .sulfate_concentration_mol_s_per_m3 = 0.5,
        .sulfate_transformation_mol_s_per_m3_step = 0.1,
        .gypsum_concentration_mol_mineral_per_m3 = 0.2,
        .divalent_activity_coefficient = 0.5,
        .gypsum_solubility_product = 0.03,
        .substrate_limit_fraction_per_step = 0.2,
        .maximum_precipitation_mol_mineral_per_m3_step = 0.1,
    };
}

test "SOLUTE gypsum correction preserves every source expression" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const sulfate =
        inputs.sulfate_concentration_mol_s_per_m3 +
        inputs.sulfate_transformation_mol_s_per_m3_step;
    const sulfate_activity =
        sulfate * inputs.divalent_activity_coefficient;
    const limit =
        inputs.substrate_limit_fraction_per_step *
        @min(
            inputs.retained_pre_calcite_calcium_concentration_mol_per_m3,
            sulfate,
        );
    const equilibrium =
        inputs.gypsum_solubility_product / sulfate_activity;
    const maximum =
        inputs.maximum_precipitation_mol_mineral_per_m3_step;
    const rate = @max(
        -@max(0.0, inputs.gypsum_concentration_mol_mineral_per_m3),
        -maximum,
        @min(
            maximum,
            limit,
            (inputs.retained_pre_calcite_calcium_activity_mol_per_m3 -
                equilibrium) /
                inputs.divalent_activity_coefficient,
        ),
    );

    try std.testing.expectEqual(
        sulfate,
        result.provisional_sulfate_concentration_mol_s_per_m3,
    );
    try std.testing.expectEqual(
        sulfate_activity,
        result.provisional_sulfate_activity_mol_s_per_m3,
    );
    try std.testing.expectEqual(
        equilibrium,
        result.equilibrium_calcium_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        rate,
        result.precipitation_mol_mineral_per_m3_step,
    );
}

test "gypsum correction conserves calcium and sulfate" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const rate = result.precipitation_mol_mineral_per_m3_step;

    try std.testing.expectApproxEqAbs(
        inputs.calcium_transformation_mol_per_m3_step,
        result.calcium_transformation_mol_per_m3_step + rate,
        1.0e-16,
    );
    try std.testing.expectApproxEqAbs(
        inputs.sulfate_transformation_mol_s_per_m3_step,
        result.sulfate_transformation_mol_s_per_m3_step + rate,
        1.0e-16,
    );
    try std.testing.expectApproxEqAbs(
        inputs.gypsum_concentration_mol_mineral_per_m3,
        result.gypsum_concentration_mol_mineral_per_m3 - rate,
        1.0e-16,
    );
}

test "gypsum retains pre-calcite calcium state for its rate" {
    const inputs = testInputs();
    const baseline = try calculateSourceOrder(inputs);
    var changed_transformation = inputs;
    changed_transformation.calcium_transformation_mol_per_m3_step = -0.4;
    const changed = try calculateSourceOrder(changed_transformation);

    try std.testing.expectEqual(
        baseline.precipitation_mol_mineral_per_m3_step,
        changed.precipitation_mol_mineral_per_m3_step,
    );
    try std.testing.expectEqual(
        baseline.substrate_limit_mol_mineral_per_m3_step,
        changed.substrate_limit_mol_mineral_per_m3_step,
    );
}

test "gypsum dissolution cannot exceed existing solid" {
    var inputs = testInputs();
    inputs.gypsum_concentration_mol_mineral_per_m3 = 0.02;
    inputs.gypsum_solubility_product = 10;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectEqual(
        -inputs.gypsum_concentration_mol_mineral_per_m3,
        result.precipitation_mol_mineral_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.gypsum_concentration_mol_mineral_per_m3,
    );
}

test "gypsum correction rejects invalid state and overflow" {
    var inputs = testInputs();
    inputs.divalent_activity_coefficient = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterGypsumCorrectionInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.sulfate_transformation_mol_s_per_m3_step = -1;
    try std.testing.expectError(
        error.NonPositiveSurfaceLitterProvisionalSulfate,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.sulfate_concentration_mol_s_per_m3 =
        std.math.floatMax(f64);
    inputs.sulfate_transformation_mol_s_per_m3_step =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterGypsumCorrectionResult,
        calculateSourceOrder(inputs),
    );
}
