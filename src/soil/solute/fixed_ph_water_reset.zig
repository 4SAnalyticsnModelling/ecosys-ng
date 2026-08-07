const std = @import("std");

pub const Inputs = struct {
    prescribed_soil_ph: f64,
    monovalent_activity_coefficient: f64,
    water_activity_product_mol2_per_m6: f64,
    current_hydrogen_activity_mol_per_m3: f64,
    current_hydroxide_activity_mol_per_m3: f64,
    assembled_hydrogen_change_mol_per_m3_step: f64,
    assembled_hydroxide_change_mol_per_m3_step: f64,
};

pub const Result = struct {
    target_hydrogen_activity_mol_per_m3: f64,
    target_hydroxide_activity_mol_per_m3: f64,
    hydrogen_reset_mol_per_m3_step: f64,
    hydroxide_reset_mol_per_m3_step: f64,
};

/// Exact source-order comparator for the prescribed-pH water reset in
/// SOLUTE.F lines 3640--3650.
///
/// The source treats `10^(-pH) * 1000` as an H concentration and multiplies
/// it by the monovalent activity coefficient. Consequently, the target
/// activity represents the prescribed pH only when that coefficient is one.
/// This convention is retained here solely for traceability.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteFixedPhWaterResetInput;
    }
    if (inputs.monovalent_activity_coefficient <= 0 or
        inputs.water_activity_product_mol2_per_m6 <= 0 or
        inputs.current_hydrogen_activity_mol_per_m3 < 0 or
        inputs.current_hydroxide_activity_mol_per_m3 < 0)
        return error.InvalidFixedPhWaterResetInput;

    // SOLUTE.F 3647--3650. Preserve the source statement and operation order.
    const target_hydrogen_activity =
        std.math.pow(f64, 10.0, -inputs.prescribed_soil_ph) *
        1.0e3 *
        inputs.monovalent_activity_coefficient;
    if (!std.math.isFinite(target_hydrogen_activity) or
        target_hydrogen_activity <= 0)
        return error.InvalidFixedPhWaterResetTarget;
    const target_hydroxide_activity =
        inputs.water_activity_product_mol2_per_m6 /
        target_hydrogen_activity;
    const hydrogen_reset =
        target_hydrogen_activity -
        inputs.current_hydrogen_activity_mol_per_m3 -
        inputs.assembled_hydrogen_change_mol_per_m3_step;
    const hydroxide_reset =
        target_hydroxide_activity -
        inputs.current_hydroxide_activity_mol_per_m3 -
        inputs.assembled_hydroxide_change_mol_per_m3_step;

    const result: Result = .{
        .target_hydrogen_activity_mol_per_m3 = target_hydrogen_activity,
        .target_hydroxide_activity_mol_per_m3 = target_hydroxide_activity,
        .hydrogen_reset_mol_per_m3_step = hydrogen_reset,
        .hydroxide_reset_mol_per_m3_step = hydroxide_reset,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteFixedPhWaterResetResult;
    }
    return result;
}

test "fixed-pH water reset preserves every source operation" {
    const inputs: Inputs = .{
        .prescribed_soil_ph = 6.5,
        .monovalent_activity_coefficient = 0.8,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
        .current_hydrogen_activity_mol_per_m3 = 2.0e-4,
        .current_hydroxide_activity_mol_per_m3 = 3.0e-5,
        .assembled_hydrogen_change_mol_per_m3_step = 2.0e-5,
        .assembled_hydroxide_change_mol_per_m3_step = -1.0e-5,
    };
    const result = try calculateSourceOrder(inputs);
    const expected_hydrogen =
        std.math.pow(f64, 10.0, -inputs.prescribed_soil_ph) *
        1.0e3 *
        inputs.monovalent_activity_coefficient;
    const expected_hydroxide =
        inputs.water_activity_product_mol2_per_m6 / expected_hydrogen;

    try std.testing.expectEqual(
        expected_hydrogen,
        result.target_hydrogen_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        expected_hydroxide,
        result.target_hydroxide_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        expected_hydrogen -
            inputs.current_hydrogen_activity_mol_per_m3 -
            inputs.assembled_hydrogen_change_mol_per_m3_step,
        result.hydrogen_reset_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        expected_hydroxide -
            inputs.current_hydroxide_activity_mol_per_m3 -
            inputs.assembled_hydroxide_change_mol_per_m3_step,
        result.hydroxide_reset_mol_per_m3_step,
    );
}

test "fixed-pH water reset closes both targets and water activity product" {
    const inputs: Inputs = .{
        .prescribed_soil_ph = 7.2,
        .monovalent_activity_coefficient = 0.65,
        .water_activity_product_mol2_per_m6 = 1.3e-8,
        .current_hydrogen_activity_mol_per_m3 = 4.0e-5,
        .current_hydroxide_activity_mol_per_m3 = 1.0e-4,
        .assembled_hydrogen_change_mol_per_m3_step = -3.0e-6,
        .assembled_hydroxide_change_mol_per_m3_step = 7.0e-6,
    };
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectApproxEqRel(
        inputs.water_activity_product_mol2_per_m6,
        result.target_hydrogen_activity_mol_per_m3 *
            result.target_hydroxide_activity_mol_per_m3,
        2 * std.math.floatEps(f64),
    );
    try std.testing.expectApproxEqAbs(
        result.target_hydrogen_activity_mol_per_m3,
        inputs.current_hydrogen_activity_mol_per_m3 +
            inputs.assembled_hydrogen_change_mol_per_m3_step +
            result.hydrogen_reset_mol_per_m3_step,
        1.0e-18,
    );
    try std.testing.expectApproxEqAbs(
        result.target_hydroxide_activity_mol_per_m3,
        inputs.current_hydroxide_activity_mol_per_m3 +
            inputs.assembled_hydroxide_change_mol_per_m3_step +
            result.hydroxide_reset_mol_per_m3_step,
        1.0e-18,
    );
}

test "fixed-pH source convention shifts activity pH by coefficient" {
    const result = try calculateSourceOrder(.{
        .prescribed_soil_ph = 6.0,
        .monovalent_activity_coefficient = 0.5,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
        .current_hydrogen_activity_mol_per_m3 = 0,
        .current_hydroxide_activity_mol_per_m3 = 0,
        .assembled_hydrogen_change_mol_per_m3_step = 0,
        .assembled_hydroxide_change_mol_per_m3_step = 0,
    });
    const standard_activity_target_mol_per_m3 = 1.0e-3;
    const implied_activity_ph =
        -@log10(result.target_hydrogen_activity_mol_per_m3 * 1.0e-3);

    try std.testing.expectApproxEqAbs(
        5.0e-4,
        result.target_hydrogen_activity_mol_per_m3,
        1.0e-18,
    );
    try std.testing.expectApproxEqAbs(6.301029995663981, implied_activity_ph, 1.0e-15);
    try std.testing.expectApproxEqAbs(
        5.0e-4,
        standard_activity_target_mol_per_m3 -
            result.target_hydrogen_activity_mol_per_m3,
        1.0e-18,
    );
}

test "fixed-pH water reset rejects invalid domains and non-finite state" {
    const valid: Inputs = .{
        .prescribed_soil_ph = 7.0,
        .monovalent_activity_coefficient = 0.8,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
        .current_hydrogen_activity_mol_per_m3 = 1.0e-4,
        .current_hydroxide_activity_mol_per_m3 = 1.0e-4,
        .assembled_hydrogen_change_mol_per_m3_step = 0,
        .assembled_hydroxide_change_mol_per_m3_step = 0,
    };
    var invalid = valid;
    invalid.monovalent_activity_coefficient = 0;
    try std.testing.expectError(
        error.InvalidFixedPhWaterResetInput,
        calculateSourceOrder(invalid),
    );
    invalid = valid;
    invalid.water_activity_product_mol2_per_m6 = 0;
    try std.testing.expectError(
        error.InvalidFixedPhWaterResetInput,
        calculateSourceOrder(invalid),
    );
    invalid = valid;
    invalid.current_hydrogen_activity_mol_per_m3 = -1;
    try std.testing.expectError(
        error.InvalidFixedPhWaterResetInput,
        calculateSourceOrder(invalid),
    );
    invalid = valid;
    invalid.assembled_hydroxide_change_mol_per_m3_step = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteFixedPhWaterResetInput,
        calculateSourceOrder(invalid),
    );
    invalid = valid;
    invalid.prescribed_soil_ph = 400;
    try std.testing.expectError(
        error.InvalidFixedPhWaterResetTarget,
        calculateSourceOrder(invalid),
    );
}
