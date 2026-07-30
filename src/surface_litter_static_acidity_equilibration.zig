const std = @import("std");

pub const SaltEquilibriumMode = enum {
    static_concentrations,
    dynamic_equilibria,
};

pub const Inputs = struct {
    mode: SaltEquilibriumMode,
    soil_ph: f64,
    mol_per_liter_to_mol_per_m3: f64,
    monovalent_activity_coefficient: f64,
    water_dissociation_product_mol2_per_m6: f64,
    existing_hydrogen_activity_mol_per_m3: f64,
    existing_hydroxide_activity_mol_per_m3: f64,
    hydrogen_transformation_mol_per_m3_step: f64,
    hydroxide_transformation_mol_per_m3_step: f64,
};

pub const StaticEquilibration = struct {
    target_hydrogen_activity_mol_per_m3: f64,
    target_hydroxide_activity_mol_per_m3: f64,
    hydrogen_equilibration_mol_per_m3_step: f64,
    hydroxide_equilibration_mol_per_m3_step: f64,
};

pub const Result = union(enum) {
    dynamic_equilibria_skipped,
    static: StaticEquilibration,
};

/// Direct source-order translation of SOLUTE.F lines 5100--5113.
///
/// SOLUTE.F 5076--5090 only conditionally prints mineral diagnostics and
/// mutates no scientific state. Structured results from the preceding mineral
/// kernels replace those scheduled console writes; this numerical kernel
/// preserves the following static-salt branch.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);

    // SOLUTE.F 5100.
    if (inputs.mode == .dynamic_equilibria)
        return .dynamic_equilibria_skipped;

    // SOLUTE.F 5102--5105.
    const target_hydrogen_activity =
        std.math.pow(f64, 10.0, -inputs.soil_ph) *
        inputs.mol_per_liter_to_mol_per_m3 *
        inputs.monovalent_activity_coefficient;
    const target_hydroxide_activity =
        inputs.water_dissociation_product_mol2_per_m6 /
        target_hydrogen_activity;
    const hydrogen_equilibration =
        target_hydrogen_activity -
        inputs.existing_hydrogen_activity_mol_per_m3 -
        inputs.hydrogen_transformation_mol_per_m3_step;
    const hydroxide_equilibration =
        target_hydroxide_activity -
        inputs.existing_hydroxide_activity_mol_per_m3 -
        inputs.hydroxide_transformation_mol_per_m3_step;

    const result: StaticEquilibration = .{
        .target_hydrogen_activity_mol_per_m3 = target_hydrogen_activity,
        .target_hydroxide_activity_mol_per_m3 = target_hydroxide_activity,
        .hydrogen_equilibration_mol_per_m3_step = hydrogen_equilibration,
        .hydroxide_equilibration_mol_per_m3_step = hydroxide_equilibration,
    };
    try validateResult(result);
    return .{ .static = result };
}

fn validateInputs(inputs: Inputs) !void {
    inline for (.{
        inputs.soil_ph,
        inputs.mol_per_liter_to_mol_per_m3,
        inputs.monovalent_activity_coefficient,
        inputs.water_dissociation_product_mol2_per_m6,
        inputs.existing_hydrogen_activity_mol_per_m3,
        inputs.existing_hydroxide_activity_mol_per_m3,
        inputs.hydrogen_transformation_mol_per_m3_step,
        inputs.hydroxide_transformation_mol_per_m3_step,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.InvalidSurfaceLitterStaticAcidityInput;
    }
    if (inputs.mol_per_liter_to_mol_per_m3 <= 0 or
        inputs.monovalent_activity_coefficient <= 0 or
        inputs.water_dissociation_product_mol2_per_m6 <= 0 or
        inputs.existing_hydrogen_activity_mol_per_m3 < 0 or
        inputs.existing_hydroxide_activity_mol_per_m3 < 0)
    {
        return error.InvalidSurfaceLitterStaticAcidityInput;
    }
}

fn validateResult(result: StaticEquilibration) !void {
    inline for (@typeInfo(StaticEquilibration).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceLitterStaticAcidityResult;
    }
    if (result.target_hydrogen_activity_mol_per_m3 <= 0 or
        result.target_hydroxide_activity_mol_per_m3 <= 0)
    {
        return error.InvalidSurfaceLitterStaticAcidityResult;
    }
}

fn testInputs() Inputs {
    return .{
        .mode = .static_concentrations,
        .soil_ph = 7,
        .mol_per_liter_to_mol_per_m3 = 1.0e3,
        .monovalent_activity_coefficient = 0.8,
        .water_dissociation_product_mol2_per_m6 = 1.0e-8,
        .existing_hydrogen_activity_mol_per_m3 = 2.0e-5,
        .existing_hydroxide_activity_mol_per_m3 = 3.0e-5,
        .hydrogen_transformation_mol_per_m3_step = 1.0e-5,
        .hydroxide_transformation_mol_per_m3_step = -1.0e-5,
    };
}

test "SOLUTE static acidity preserves every source expression" {
    const inputs = testInputs();
    const result = (try calculateSourceOrder(inputs)).static;
    const target_hydrogen =
        std.math.pow(f64, 10.0, -inputs.soil_ph) *
        inputs.mol_per_liter_to_mol_per_m3 *
        inputs.monovalent_activity_coefficient;
    const target_hydroxide =
        inputs.water_dissociation_product_mol2_per_m6 /
        target_hydrogen;

    try std.testing.expectEqual(
        target_hydrogen,
        result.target_hydrogen_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        target_hydroxide,
        result.target_hydroxide_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        target_hydrogen -
            inputs.existing_hydrogen_activity_mol_per_m3 -
            inputs.hydrogen_transformation_mol_per_m3_step,
        result.hydrogen_equilibration_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        target_hydroxide -
            inputs.existing_hydroxide_activity_mol_per_m3 -
            inputs.hydroxide_transformation_mol_per_m3_step,
        result.hydroxide_equilibration_mol_per_m3_step,
    );
}

test "static acidity fluxes reach target activities and water product" {
    const inputs = testInputs();
    const result = (try calculateSourceOrder(inputs)).static;
    const reached_hydrogen =
        inputs.existing_hydrogen_activity_mol_per_m3 +
        inputs.hydrogen_transformation_mol_per_m3_step +
        result.hydrogen_equilibration_mol_per_m3_step;
    const reached_hydroxide =
        inputs.existing_hydroxide_activity_mol_per_m3 +
        inputs.hydroxide_transformation_mol_per_m3_step +
        result.hydroxide_equilibration_mol_per_m3_step;

    try std.testing.expectApproxEqAbs(
        result.target_hydrogen_activity_mol_per_m3,
        reached_hydrogen,
        1.0e-20,
    );
    try std.testing.expectApproxEqAbs(
        result.target_hydroxide_activity_mol_per_m3,
        reached_hydroxide,
        1.0e-20,
    );
    try std.testing.expectApproxEqAbs(
        inputs.water_dissociation_product_mol2_per_m6,
        result.target_hydrogen_activity_mol_per_m3 *
            result.target_hydroxide_activity_mol_per_m3,
        4.0e-24,
    );
}

test "dynamic salt mode skips static acidity equilibration" {
    var inputs = testInputs();
    inputs.mode = .dynamic_equilibria;

    try std.testing.expectEqual(
        Result.dynamic_equilibria_skipped,
        try calculateSourceOrder(inputs),
    );
}

test "static acidity rejects invalid input and pH overflow" {
    var inputs = testInputs();
    inputs.monovalent_activity_coefficient = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticAcidityInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.existing_hydrogen_activity_mol_per_m3 = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticAcidityInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.soil_ph = -400;
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterStaticAcidityResult,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.soil_ph = 400;
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterStaticAcidityResult,
        calculateSourceOrder(inputs),
    );
}
