const std = @import("std");

pub const Inputs = struct {
    aluminum_concentration_mol_per_m3: f64,
    aluminum_transformation_mol_per_m3_step: f64,
    hydroxide_transformation_mol_per_m3_step: f64,
    hydroxide_activity_mol_per_m3: f64,
    gibbsite_concentration_mol_mineral_per_m3: f64,
    trivalent_activity_coefficient: f64,
    gibbsite_solubility_product: f64,
    substrate_limit_fraction_per_step: f64,
    maximum_precipitation_mol_mineral_per_m3_step: f64,
};

pub const Result = struct {
    provisional_aluminum_concentration_mol_per_m3: f64,
    provisional_aluminum_activity_mol_per_m3: f64,
    substrate_limit_mol_mineral_per_m3_step: f64,
    equilibrium_aluminum_activity_mol_per_m3: f64,
    precipitation_mol_mineral_per_m3_step: f64,
    aluminum_transformation_mol_per_m3_step: f64,
    hydroxide_transformation_mol_per_m3_step: f64,
    gibbsite_concentration_mol_mineral_per_m3: f64,
};

/// Direct source-order translation of SOLUTE.F lines 5027--5035.
///
/// Positive precipitation removes one aqueous Al and three OH per gibbsite;
/// negative precipitation dissolves the existing mineral.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);

    // SOLUTE.F 5027--5032.
    const provisional_aluminum =
        inputs.aluminum_concentration_mol_per_m3 +
        inputs.aluminum_transformation_mol_per_m3_step;
    if (!std.math.isFinite(provisional_aluminum))
        return error.NonFiniteSurfaceLitterGibbsiteCorrectionResult;
    if (provisional_aluminum < 0)
        return error.NegativeSurfaceLitterProvisionalAluminum;

    const aluminum_activity =
        provisional_aluminum * inputs.trivalent_activity_coefficient;
    const substrate_limit =
        inputs.substrate_limit_fraction_per_step * provisional_aluminum;
    const equilibrium_aluminum_activity =
        inputs.gibbsite_solubility_product /
        std.math.pow(f64, inputs.hydroxide_activity_mol_per_m3, 3.0);
    const maximum =
        inputs.maximum_precipitation_mol_mineral_per_m3_step;
    const precipitation = @max(
        -@max(0.0, inputs.gibbsite_concentration_mol_mineral_per_m3),
        -maximum,
        @min(
            maximum,
            substrate_limit,
            (aluminum_activity - equilibrium_aluminum_activity) /
                inputs.trivalent_activity_coefficient,
        ),
    );

    // SOLUTE.F 5033--5035.
    const result: Result = .{
        .provisional_aluminum_concentration_mol_per_m3 = provisional_aluminum,
        .provisional_aluminum_activity_mol_per_m3 = aluminum_activity,
        .substrate_limit_mol_mineral_per_m3_step = substrate_limit,
        .equilibrium_aluminum_activity_mol_per_m3 = equilibrium_aluminum_activity,
        .precipitation_mol_mineral_per_m3_step = precipitation,
        .aluminum_transformation_mol_per_m3_step = inputs.aluminum_transformation_mol_per_m3_step -
            precipitation,
        .hydroxide_transformation_mol_per_m3_step = inputs.hydroxide_transformation_mol_per_m3_step -
            3.0 * precipitation,
        .gibbsite_concentration_mol_mineral_per_m3 = inputs.gibbsite_concentration_mol_mineral_per_m3 +
            precipitation,
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidSurfaceLitterGibbsiteCorrectionInput;
    }
    if (inputs.aluminum_concentration_mol_per_m3 < 0 or
        inputs.hydroxide_activity_mol_per_m3 <= 0 or
        inputs.gibbsite_concentration_mol_mineral_per_m3 < 0 or
        inputs.trivalent_activity_coefficient <= 0 or
        inputs.gibbsite_solubility_product <= 0 or
        inputs.substrate_limit_fraction_per_step < 0 or
        inputs.substrate_limit_fraction_per_step > 1 or
        inputs.maximum_precipitation_mol_mineral_per_m3_step < 0)
    {
        return error.InvalidSurfaceLitterGibbsiteCorrectionInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceLitterGibbsiteCorrectionResult;
    }
    if (result.provisional_aluminum_concentration_mol_per_m3 < 0 or
        result.provisional_aluminum_activity_mol_per_m3 < 0 or
        result.substrate_limit_mol_mineral_per_m3_step < 0 or
        result.equilibrium_aluminum_activity_mol_per_m3 < 0 or
        result.gibbsite_concentration_mol_mineral_per_m3 < 0)
    {
        return error.InvalidSurfaceLitterGibbsiteCorrectionResult;
    }
}

fn testInputs() Inputs {
    return .{
        .aluminum_concentration_mol_per_m3 = 0.3,
        .aluminum_transformation_mol_per_m3_step = 0.1,
        .hydroxide_transformation_mol_per_m3_step = 0.2,
        .hydroxide_activity_mol_per_m3 = 0.5,
        .gibbsite_concentration_mol_mineral_per_m3 = 0.2,
        .trivalent_activity_coefficient = 0.5,
        .gibbsite_solubility_product = 0.001,
        .substrate_limit_fraction_per_step = 0.2,
        .maximum_precipitation_mol_mineral_per_m3_step = 0.05,
    };
}

test "SOLUTE gibbsite correction preserves every source expression" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const provisional =
        inputs.aluminum_concentration_mol_per_m3 +
        inputs.aluminum_transformation_mol_per_m3_step;
    const activity = provisional * inputs.trivalent_activity_coefficient;
    const limit = inputs.substrate_limit_fraction_per_step * provisional;
    const equilibrium = inputs.gibbsite_solubility_product /
        std.math.pow(f64, inputs.hydroxide_activity_mol_per_m3, 3.0);
    const maximum =
        inputs.maximum_precipitation_mol_mineral_per_m3_step;
    const rate = @max(
        -@max(0.0, inputs.gibbsite_concentration_mol_mineral_per_m3),
        -maximum,
        @min(
            maximum,
            limit,
            (activity - equilibrium) /
                inputs.trivalent_activity_coefficient,
        ),
    );

    try std.testing.expectEqual(
        provisional,
        result.provisional_aluminum_concentration_mol_per_m3,
    );
    try std.testing.expectEqual(
        activity,
        result.provisional_aluminum_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        equilibrium,
        result.equilibrium_aluminum_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        rate,
        result.precipitation_mol_mineral_per_m3_step,
    );
    try std.testing.expectEqual(
        inputs.aluminum_transformation_mol_per_m3_step - rate,
        result.aluminum_transformation_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        inputs.hydroxide_transformation_mol_per_m3_step - 3.0 * rate,
        result.hydroxide_transformation_mol_per_m3_step,
    );
}

test "gibbsite correction conserves aluminum and hydroxide stoichiometry" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const rate = result.precipitation_mol_mineral_per_m3_step;

    try std.testing.expectApproxEqAbs(
        inputs.aluminum_transformation_mol_per_m3_step,
        result.aluminum_transformation_mol_per_m3_step + rate,
        1.0e-16,
    );
    try std.testing.expectApproxEqAbs(
        inputs.hydroxide_transformation_mol_per_m3_step,
        result.hydroxide_transformation_mol_per_m3_step + 3.0 * rate,
        1.0e-16,
    );
    try std.testing.expectApproxEqAbs(
        inputs.gibbsite_concentration_mol_mineral_per_m3,
        result.gibbsite_concentration_mol_mineral_per_m3 - rate,
        1.0e-16,
    );
}

test "gibbsite dissolution cannot exceed existing solid" {
    var inputs = testInputs();
    inputs.gibbsite_concentration_mol_mineral_per_m3 = 0.02;
    inputs.gibbsite_solubility_product = 1;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectEqual(
        -inputs.gibbsite_concentration_mol_mineral_per_m3,
        result.precipitation_mol_mineral_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.gibbsite_concentration_mol_mineral_per_m3,
    );
}

test "gibbsite correction rejects invalid state and overflow" {
    var inputs = testInputs();
    inputs.hydroxide_activity_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterGibbsiteCorrectionInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.aluminum_transformation_mol_per_m3_step = -1;
    try std.testing.expectError(
        error.NegativeSurfaceLitterProvisionalAluminum,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.aluminum_concentration_mol_per_m3 = std.math.floatMax(f64);
    inputs.aluminum_transformation_mol_per_m3_step =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterGibbsiteCorrectionResult,
        calculateSourceOrder(inputs),
    );
}
