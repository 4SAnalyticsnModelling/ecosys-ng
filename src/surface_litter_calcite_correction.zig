const std = @import("std");

pub const Inputs = struct {
    calcium_concentration_mol_per_m3: f64,
    calcium_transformation_mol_per_m3_step: f64,
    carbonate_concentration_mol_c_per_m3: f64,
    carbonate_transformation_mol_c_per_m3_step: f64,
    hydroxide_activity_mol_per_m3: f64,
    calcite_concentration_mol_mineral_per_m3: f64,
    divalent_activity_coefficient: f64,
    calcite_solubility_product: f64,
    hydroxide_inhibition_constant_mol_per_m3: f64,
    substrate_limit_fraction_per_step: f64,
    maximum_precipitation_mol_mineral_per_m3_step: f64,
};

pub const Result = struct {
    maximum_dissolution_mol_mineral_per_m3_step: f64,
    provisional_calcium_concentration_mol_per_m3: f64,
    provisional_carbonate_concentration_mol_c_per_m3: f64,
    provisional_calcium_activity_mol_per_m3: f64,
    provisional_carbonate_activity_mol_c_per_m3: f64,
    substrate_limit_mol_mineral_per_m3_step: f64,
    equilibrium_calcium_activity_mol_per_m3: f64,
    precipitation_mol_mineral_per_m3_step: f64,
    calcium_transformation_mol_per_m3_step: f64,
    carbonate_transformation_mol_c_per_m3_step: f64,
    calcite_concentration_mol_mineral_per_m3: f64,
};

/// Direct source-order translation of SOLUTE.F lines 5051--5063.
///
/// The source applies hydroxide inhibition only to calcite dissolution.
/// Positive rates precipitate one Ca and one carbonate per calcite.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);

    // SOLUTE.F 5051--5057.
    const maximum_dissolution =
        inputs.maximum_precipitation_mol_mineral_per_m3_step /
        (1.0 +
            inputs.hydroxide_activity_mol_per_m3 /
                inputs.hydroxide_inhibition_constant_mol_per_m3);
    const provisional_calcium =
        inputs.calcium_concentration_mol_per_m3 +
        inputs.calcium_transformation_mol_per_m3_step;
    const provisional_carbonate =
        inputs.carbonate_concentration_mol_c_per_m3 +
        inputs.carbonate_transformation_mol_c_per_m3_step;
    if (!std.math.isFinite(provisional_calcium) or
        !std.math.isFinite(provisional_carbonate))
    {
        return error.NonFiniteSurfaceLitterCalciteCorrectionResult;
    }
    if (provisional_calcium < 0)
        return error.NegativeSurfaceLitterProvisionalCalcium;
    if (provisional_carbonate <= 0)
        return error.NonPositiveSurfaceLitterProvisionalCarbonate;

    const calcium_activity =
        provisional_calcium * inputs.divalent_activity_coefficient;
    const carbonate_activity =
        provisional_carbonate * inputs.divalent_activity_coefficient;
    const substrate_limit =
        inputs.substrate_limit_fraction_per_step *
        @min(provisional_calcium, provisional_carbonate);

    // SOLUTE.F 5058--5060.
    const equilibrium_calcium_activity =
        inputs.calcite_solubility_product / carbonate_activity;
    const maximum_precipitation =
        inputs.maximum_precipitation_mol_mineral_per_m3_step;
    const precipitation = @max(
        -@max(0.0, inputs.calcite_concentration_mol_mineral_per_m3),
        -maximum_dissolution,
        @min(
            maximum_precipitation,
            substrate_limit,
            (calcium_activity - equilibrium_calcium_activity) /
                inputs.divalent_activity_coefficient,
        ),
    );

    // SOLUTE.F 5061--5063.
    const result: Result = .{
        .maximum_dissolution_mol_mineral_per_m3_step = maximum_dissolution,
        .provisional_calcium_concentration_mol_per_m3 = provisional_calcium,
        .provisional_carbonate_concentration_mol_c_per_m3 = provisional_carbonate,
        .provisional_calcium_activity_mol_per_m3 = calcium_activity,
        .provisional_carbonate_activity_mol_c_per_m3 = carbonate_activity,
        .substrate_limit_mol_mineral_per_m3_step = substrate_limit,
        .equilibrium_calcium_activity_mol_per_m3 = equilibrium_calcium_activity,
        .precipitation_mol_mineral_per_m3_step = precipitation,
        .calcium_transformation_mol_per_m3_step = inputs.calcium_transformation_mol_per_m3_step -
            precipitation,
        .carbonate_transformation_mol_c_per_m3_step = inputs.carbonate_transformation_mol_c_per_m3_step -
            precipitation,
        .calcite_concentration_mol_mineral_per_m3 = inputs.calcite_concentration_mol_mineral_per_m3 +
            precipitation,
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidSurfaceLitterCalciteCorrectionInput;
    }
    if (inputs.calcium_concentration_mol_per_m3 < 0 or
        inputs.carbonate_concentration_mol_c_per_m3 < 0 or
        inputs.hydroxide_activity_mol_per_m3 < 0 or
        inputs.calcite_concentration_mol_mineral_per_m3 < 0 or
        inputs.divalent_activity_coefficient <= 0 or
        inputs.calcite_solubility_product <= 0 or
        inputs.hydroxide_inhibition_constant_mol_per_m3 <= 0 or
        inputs.substrate_limit_fraction_per_step < 0 or
        inputs.substrate_limit_fraction_per_step > 1 or
        inputs.maximum_precipitation_mol_mineral_per_m3_step < 0)
    {
        return error.InvalidSurfaceLitterCalciteCorrectionInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceLitterCalciteCorrectionResult;
    }
    if (result.maximum_dissolution_mol_mineral_per_m3_step < 0 or
        result.provisional_calcium_concentration_mol_per_m3 < 0 or
        result.provisional_carbonate_concentration_mol_c_per_m3 <= 0 or
        result.provisional_calcium_activity_mol_per_m3 < 0 or
        result.provisional_carbonate_activity_mol_c_per_m3 <= 0 or
        result.substrate_limit_mol_mineral_per_m3_step < 0 or
        result.equilibrium_calcium_activity_mol_per_m3 < 0 or
        result.calcite_concentration_mol_mineral_per_m3 < 0)
    {
        return error.InvalidSurfaceLitterCalciteCorrectionResult;
    }
}

fn testInputs() Inputs {
    return .{
        .calcium_concentration_mol_per_m3 = 0.5,
        .calcium_transformation_mol_per_m3_step = 0.1,
        .carbonate_concentration_mol_c_per_m3 = 0.4,
        .carbonate_transformation_mol_c_per_m3_step = 0.1,
        .hydroxide_activity_mol_per_m3 = 0.2,
        .calcite_concentration_mol_mineral_per_m3 = 0.3,
        .divalent_activity_coefficient = 0.5,
        .calcite_solubility_product = 0.025,
        .hydroxide_inhibition_constant_mol_per_m3 = 0.4,
        .substrate_limit_fraction_per_step = 0.2,
        .maximum_precipitation_mol_mineral_per_m3_step = 0.1,
    };
}

test "SOLUTE calcite correction preserves every source expression" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const maximum_dissolution =
        inputs.maximum_precipitation_mol_mineral_per_m3_step /
        (1.0 +
            inputs.hydroxide_activity_mol_per_m3 /
                inputs.hydroxide_inhibition_constant_mol_per_m3);
    const calcium =
        inputs.calcium_concentration_mol_per_m3 +
        inputs.calcium_transformation_mol_per_m3_step;
    const carbonate =
        inputs.carbonate_concentration_mol_c_per_m3 +
        inputs.carbonate_transformation_mol_c_per_m3_step;
    const calcium_activity = calcium * inputs.divalent_activity_coefficient;
    const carbonate_activity =
        carbonate * inputs.divalent_activity_coefficient;
    const limit =
        inputs.substrate_limit_fraction_per_step *
        @min(calcium, carbonate);
    const equilibrium =
        inputs.calcite_solubility_product / carbonate_activity;
    const rate = @max(
        -@max(0.0, inputs.calcite_concentration_mol_mineral_per_m3),
        -maximum_dissolution,
        @min(
            inputs.maximum_precipitation_mol_mineral_per_m3_step,
            limit,
            (calcium_activity - equilibrium) /
                inputs.divalent_activity_coefficient,
        ),
    );

    try std.testing.expectEqual(
        maximum_dissolution,
        result.maximum_dissolution_mol_mineral_per_m3_step,
    );
    try std.testing.expectEqual(
        calcium_activity,
        result.provisional_calcium_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        carbonate_activity,
        result.provisional_carbonate_activity_mol_c_per_m3,
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

test "calcite correction conserves calcium and carbonate" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const rate = result.precipitation_mol_mineral_per_m3_step;

    try std.testing.expectApproxEqAbs(
        inputs.calcium_transformation_mol_per_m3_step,
        result.calcium_transformation_mol_per_m3_step + rate,
        1.0e-16,
    );
    try std.testing.expectApproxEqAbs(
        inputs.carbonate_transformation_mol_c_per_m3_step,
        result.carbonate_transformation_mol_c_per_m3_step + rate,
        1.0e-16,
    );
    try std.testing.expectApproxEqAbs(
        inputs.calcite_concentration_mol_mineral_per_m3,
        result.calcite_concentration_mol_mineral_per_m3 - rate,
        1.0e-16,
    );
}

test "hydroxide inhibition applies only to calcite dissolution" {
    var inputs = testInputs();
    inputs.hydroxide_activity_mol_per_m3 = 0.8;
    inputs.hydroxide_inhibition_constant_mol_per_m3 = 0.2;
    inputs.calcite_solubility_product = 10;
    var result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        -result.maximum_dissolution_mol_mineral_per_m3_step,
        result.precipitation_mol_mineral_per_m3_step,
    );

    inputs.calcite_solubility_product = 1.0e-6;
    result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        inputs.maximum_precipitation_mol_mineral_per_m3_step,
        result.precipitation_mol_mineral_per_m3_step,
    );
}

test "calcite correction rejects invalid state and overflow" {
    var inputs = testInputs();
    inputs.divalent_activity_coefficient = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterCalciteCorrectionInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.carbonate_transformation_mol_c_per_m3_step = -1;
    try std.testing.expectError(
        error.NonPositiveSurfaceLitterProvisionalCarbonate,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.calcium_concentration_mol_per_m3 = std.math.floatMax(f64);
    inputs.calcium_transformation_mol_per_m3_step =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterCalciteCorrectionResult,
        calculateSourceOrder(inputs),
    );
}
