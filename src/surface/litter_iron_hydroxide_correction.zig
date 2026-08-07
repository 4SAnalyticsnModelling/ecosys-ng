const std = @import("std");

pub const Inputs = struct {
    iron_concentration_mol_per_m3: f64,
    iron_transformation_mol_per_m3_step: f64,
    hydroxide_transformation_mol_per_m3_step: f64,
    hydroxide_activity_mol_per_m3: f64,
    iron_hydroxide_concentration_mol_mineral_per_m3: f64,
    trivalent_activity_coefficient: f64,
    iron_hydroxide_solubility_product: f64,
    substrate_limit_fraction_per_step: f64,
    maximum_precipitation_mol_mineral_per_m3_step: f64,
};

pub const Result = struct {
    provisional_iron_concentration_mol_per_m3: f64,
    provisional_iron_activity_mol_per_m3: f64,
    substrate_limit_mol_mineral_per_m3_step: f64,
    equilibrium_iron_activity_mol_per_m3: f64,
    precipitation_mol_mineral_per_m3_step: f64,
    iron_transformation_mol_per_m3_step: f64,
    hydroxide_transformation_mol_per_m3_step: f64,
    iron_hydroxide_concentration_mol_mineral_per_m3: f64,
};

/// Direct source-order translation of SOLUTE.F lines 5039--5047.
///
/// Positive precipitation removes one aqueous Fe and three OH per Fe(OH)3;
/// negative precipitation dissolves the existing mineral.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);

    // SOLUTE.F 5039--5044.
    const provisional_iron =
        inputs.iron_concentration_mol_per_m3 +
        inputs.iron_transformation_mol_per_m3_step;
    if (!std.math.isFinite(provisional_iron))
        return error.NonFiniteSurfaceLitterIronHydroxideCorrectionResult;
    if (provisional_iron < 0)
        return error.NegativeSurfaceLitterProvisionalIron;

    const iron_activity =
        provisional_iron * inputs.trivalent_activity_coefficient;
    const substrate_limit =
        inputs.substrate_limit_fraction_per_step * provisional_iron;
    const equilibrium_iron_activity =
        inputs.iron_hydroxide_solubility_product /
        std.math.pow(f64, inputs.hydroxide_activity_mol_per_m3, 3.0);
    const maximum =
        inputs.maximum_precipitation_mol_mineral_per_m3_step;
    const precipitation = @max(
        -@max(
            0.0,
            inputs.iron_hydroxide_concentration_mol_mineral_per_m3,
        ),
        -maximum,
        @min(
            maximum,
            substrate_limit,
            (iron_activity - equilibrium_iron_activity) /
                inputs.trivalent_activity_coefficient,
        ),
    );

    // SOLUTE.F 5045--5047.
    const result: Result = .{
        .provisional_iron_concentration_mol_per_m3 = provisional_iron,
        .provisional_iron_activity_mol_per_m3 = iron_activity,
        .substrate_limit_mol_mineral_per_m3_step = substrate_limit,
        .equilibrium_iron_activity_mol_per_m3 = equilibrium_iron_activity,
        .precipitation_mol_mineral_per_m3_step = precipitation,
        .iron_transformation_mol_per_m3_step = inputs.iron_transformation_mol_per_m3_step -
            precipitation,
        .hydroxide_transformation_mol_per_m3_step = inputs.hydroxide_transformation_mol_per_m3_step -
            3.0 * precipitation,
        .iron_hydroxide_concentration_mol_mineral_per_m3 = inputs.iron_hydroxide_concentration_mol_mineral_per_m3 +
            precipitation,
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidSurfaceLitterIronHydroxideCorrectionInput;
    }
    if (inputs.iron_concentration_mol_per_m3 < 0 or
        inputs.hydroxide_activity_mol_per_m3 <= 0 or
        inputs.iron_hydroxide_concentration_mol_mineral_per_m3 < 0 or
        inputs.trivalent_activity_coefficient <= 0 or
        inputs.iron_hydroxide_solubility_product <= 0 or
        inputs.substrate_limit_fraction_per_step < 0 or
        inputs.substrate_limit_fraction_per_step > 1 or
        inputs.maximum_precipitation_mol_mineral_per_m3_step < 0)
    {
        return error.InvalidSurfaceLitterIronHydroxideCorrectionInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceLitterIronHydroxideCorrectionResult;
    }
    if (result.provisional_iron_concentration_mol_per_m3 < 0 or
        result.provisional_iron_activity_mol_per_m3 < 0 or
        result.substrate_limit_mol_mineral_per_m3_step < 0 or
        result.equilibrium_iron_activity_mol_per_m3 < 0 or
        result.iron_hydroxide_concentration_mol_mineral_per_m3 < 0)
    {
        return error.InvalidSurfaceLitterIronHydroxideCorrectionResult;
    }
}

fn testInputs() Inputs {
    return .{
        .iron_concentration_mol_per_m3 = 0.3,
        .iron_transformation_mol_per_m3_step = 0.1,
        .hydroxide_transformation_mol_per_m3_step = 0.2,
        .hydroxide_activity_mol_per_m3 = 0.5,
        .iron_hydroxide_concentration_mol_mineral_per_m3 = 0.2,
        .trivalent_activity_coefficient = 0.5,
        .iron_hydroxide_solubility_product = 0.001,
        .substrate_limit_fraction_per_step = 0.2,
        .maximum_precipitation_mol_mineral_per_m3_step = 0.05,
    };
}

test "SOLUTE iron hydroxide correction preserves every source expression" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const provisional =
        inputs.iron_concentration_mol_per_m3 +
        inputs.iron_transformation_mol_per_m3_step;
    const activity = provisional * inputs.trivalent_activity_coefficient;
    const limit = inputs.substrate_limit_fraction_per_step * provisional;
    const equilibrium = inputs.iron_hydroxide_solubility_product /
        std.math.pow(f64, inputs.hydroxide_activity_mol_per_m3, 3.0);
    const maximum =
        inputs.maximum_precipitation_mol_mineral_per_m3_step;
    const rate = @max(
        -@max(
            0.0,
            inputs.iron_hydroxide_concentration_mol_mineral_per_m3,
        ),
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
        result.provisional_iron_concentration_mol_per_m3,
    );
    try std.testing.expectEqual(
        activity,
        result.provisional_iron_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        equilibrium,
        result.equilibrium_iron_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        rate,
        result.precipitation_mol_mineral_per_m3_step,
    );
    try std.testing.expectEqual(
        inputs.iron_transformation_mol_per_m3_step - rate,
        result.iron_transformation_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        inputs.hydroxide_transformation_mol_per_m3_step - 3.0 * rate,
        result.hydroxide_transformation_mol_per_m3_step,
    );
}

test "iron hydroxide correction conserves iron and hydroxide" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const rate = result.precipitation_mol_mineral_per_m3_step;

    try std.testing.expectApproxEqAbs(
        inputs.iron_transformation_mol_per_m3_step,
        result.iron_transformation_mol_per_m3_step + rate,
        1.0e-16,
    );
    try std.testing.expectApproxEqAbs(
        inputs.hydroxide_transformation_mol_per_m3_step,
        result.hydroxide_transformation_mol_per_m3_step + 3.0 * rate,
        1.0e-16,
    );
    try std.testing.expectApproxEqAbs(
        inputs.iron_hydroxide_concentration_mol_mineral_per_m3,
        result.iron_hydroxide_concentration_mol_mineral_per_m3 - rate,
        1.0e-16,
    );
}

test "iron hydroxide dissolution cannot exceed existing solid" {
    var inputs = testInputs();
    inputs.iron_hydroxide_concentration_mol_mineral_per_m3 = 0.02;
    inputs.iron_hydroxide_solubility_product = 1;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectEqual(
        -inputs.iron_hydroxide_concentration_mol_mineral_per_m3,
        result.precipitation_mol_mineral_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.iron_hydroxide_concentration_mol_mineral_per_m3,
    );
}

test "iron hydroxide correction rejects invalid state and overflow" {
    var inputs = testInputs();
    inputs.hydroxide_activity_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterIronHydroxideCorrectionInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.iron_transformation_mol_per_m3_step = -1;
    try std.testing.expectError(
        error.NegativeSurfaceLitterProvisionalIron,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.iron_concentration_mol_per_m3 = std.math.floatMax(f64);
    inputs.iron_transformation_mol_per_m3_step =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterIronHydroxideCorrectionResult,
        calculateSourceOrder(inputs),
    );
}
