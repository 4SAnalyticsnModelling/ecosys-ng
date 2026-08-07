const std = @import("std");

pub const Inputs = struct {
    litter_water_volume_m3: f64,
    phosphorus_molar_mass_g_p_per_mol: f64,
    hydrogen_phosphate_net_change_g_p_per_step: f64,
    dihydrogen_phosphate_net_change_g_p_per_step: f64,
    hydrogen_phosphate_inventory_g_p: f64,
    dihydrogen_phosphate_inventory_g_p: f64,
    aqueous_phosphorus_floor_mol_p_per_m3: f64,
};

pub const Result = struct {
    phosphorus_water_scale_g_p_m3_per_mol: f64,
    hydrogen_phosphate_input_mol_p_per_m3_step: f64,
    dihydrogen_phosphate_input_mol_p_per_m3_step: f64,
    aqueous_hydrogen_phosphate_mol_p_per_m3: f64,
    aqueous_dihydrogen_phosphate_mol_p_per_m3: f64,
    hydrogen_phosphate_floor_was_applied: bool,
    dihydrogen_phosphate_floor_was_applied: bool,
};

/// Direct source-order translation of SOLUTE.F lines 4149--4161.
///
/// The caller selects one runtime horizontal cell already admitted by the
/// surrounding wet-litter branch. No inventory is mutated by this kernel.
pub fn calculate(inputs: Inputs) !Result {
    try validateInputs(inputs);

    // SOLUTE.F 4157--4161.
    const phosphorus_water_scale =
        inputs.phosphorus_molar_mass_g_p_per_mol *
        inputs.litter_water_volume_m3;
    if (!std.math.isFinite(phosphorus_water_scale) or
        phosphorus_water_scale <= 0)
    {
        return error.NonFiniteSurfaceLitterPhosphateConcentration;
    }
    const hydrogen_phosphate_input =
        inputs.hydrogen_phosphate_net_change_g_p_per_step /
        phosphorus_water_scale;
    const dihydrogen_phosphate_input =
        inputs.dihydrogen_phosphate_net_change_g_p_per_step /
        phosphorus_water_scale;
    const unconstrained_hydrogen_phosphate =
        inputs.hydrogen_phosphate_inventory_g_p /
        phosphorus_water_scale +
        hydrogen_phosphate_input;
    const unconstrained_dihydrogen_phosphate =
        inputs.dihydrogen_phosphate_inventory_g_p /
        phosphorus_water_scale +
        dihydrogen_phosphate_input;
    inline for (.{
        hydrogen_phosphate_input,
        dihydrogen_phosphate_input,
        unconstrained_hydrogen_phosphate,
        unconstrained_dihydrogen_phosphate,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterPhosphateConcentration;
    }

    const result: Result = .{
        .phosphorus_water_scale_g_p_m3_per_mol = phosphorus_water_scale,
        .hydrogen_phosphate_input_mol_p_per_m3_step = hydrogen_phosphate_input,
        .dihydrogen_phosphate_input_mol_p_per_m3_step = dihydrogen_phosphate_input,
        .aqueous_hydrogen_phosphate_mol_p_per_m3 = @max(
            inputs.aqueous_phosphorus_floor_mol_p_per_m3,
            unconstrained_hydrogen_phosphate,
        ),
        .aqueous_dihydrogen_phosphate_mol_p_per_m3 = @max(
            inputs.aqueous_phosphorus_floor_mol_p_per_m3,
            unconstrained_dihydrogen_phosphate,
        ),
        .hydrogen_phosphate_floor_was_applied = unconstrained_hydrogen_phosphate <
            inputs.aqueous_phosphorus_floor_mol_p_per_m3,
        .dihydrogen_phosphate_floor_was_applied = unconstrained_dihydrogen_phosphate <
            inputs.aqueous_phosphorus_floor_mol_p_per_m3,
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidSurfaceLitterPhosphateConcentrationInput;
    }
    if (inputs.litter_water_volume_m3 <= 0 or
        inputs.phosphorus_molar_mass_g_p_per_mol <= 0 or
        inputs.hydrogen_phosphate_inventory_g_p < 0 or
        inputs.dihydrogen_phosphate_inventory_g_p < 0 or
        inputs.aqueous_phosphorus_floor_mol_p_per_m3 < 0)
    {
        return error.InvalidSurfaceLitterPhosphateConcentrationInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (.{
        result.phosphorus_water_scale_g_p_m3_per_mol,
        result.hydrogen_phosphate_input_mol_p_per_m3_step,
        result.dihydrogen_phosphate_input_mol_p_per_m3_step,
        result.aqueous_hydrogen_phosphate_mol_p_per_m3,
        result.aqueous_dihydrogen_phosphate_mol_p_per_m3,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterPhosphateConcentration;
    }
    if (result.aqueous_hydrogen_phosphate_mol_p_per_m3 < 0 or
        result.aqueous_dihydrogen_phosphate_mol_p_per_m3 < 0)
    {
        return error.InvalidSurfaceLitterPhosphateConcentrationResult;
    }
}

fn testInputs() Inputs {
    return .{
        .litter_water_volume_m3 = 2,
        .phosphorus_molar_mass_g_p_per_mol = 31,
        .hydrogen_phosphate_net_change_g_p_per_step = 6.2,
        .dihydrogen_phosphate_net_change_g_p_per_step = -3.1,
        .hydrogen_phosphate_inventory_g_p = 55.8,
        .dihydrogen_phosphate_inventory_g_p = 34.1,
        .aqueous_phosphorus_floor_mol_p_per_m3 = 1.0e-20,
    };
}

test "SOLUTE surface phosphate reconstruction preserves every source expression" {
    const inputs = testInputs();
    const result = try calculate(inputs);

    const expected_scale =
        inputs.phosphorus_molar_mass_g_p_per_mol *
        inputs.litter_water_volume_m3;
    const expected_hydrogen_input =
        inputs.hydrogen_phosphate_net_change_g_p_per_step /
        expected_scale;
    const expected_dihydrogen_input =
        inputs.dihydrogen_phosphate_net_change_g_p_per_step /
        expected_scale;
    const expected_hydrogen = @max(
        inputs.aqueous_phosphorus_floor_mol_p_per_m3,
        inputs.hydrogen_phosphate_inventory_g_p /
            expected_scale +
            expected_hydrogen_input,
    );
    const expected_dihydrogen = @max(
        inputs.aqueous_phosphorus_floor_mol_p_per_m3,
        inputs.dihydrogen_phosphate_inventory_g_p /
            expected_scale +
            expected_dihydrogen_input,
    );

    try std.testing.expectEqual(
        expected_scale,
        result.phosphorus_water_scale_g_p_m3_per_mol,
    );
    try std.testing.expectEqual(
        expected_hydrogen_input,
        result.hydrogen_phosphate_input_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        expected_dihydrogen_input,
        result.dihydrogen_phosphate_input_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        expected_hydrogen,
        result.aqueous_hydrogen_phosphate_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        expected_dihydrogen,
        result.aqueous_dihydrogen_phosphate_mol_p_per_m3,
    );
}

test "surface phosphate reconstruction recovers valid extensive masses" {
    const inputs = testInputs();
    const result = try calculate(inputs);

    try std.testing.expectApproxEqAbs(
        inputs.hydrogen_phosphate_inventory_g_p +
            inputs.hydrogen_phosphate_net_change_g_p_per_step,
        result.aqueous_hydrogen_phosphate_mol_p_per_m3 *
            result.phosphorus_water_scale_g_p_m3_per_mol,
        1.0e-14,
    );
    try std.testing.expectApproxEqAbs(
        inputs.dihydrogen_phosphate_inventory_g_p +
            inputs.dihydrogen_phosphate_net_change_g_p_per_step,
        result.aqueous_dihydrogen_phosphate_mol_p_per_m3 *
            result.phosphorus_water_scale_g_p_m3_per_mol,
        1.0e-14,
    );
}

test "surface phosphate depletion floors are explicit" {
    var inputs = testInputs();
    inputs.hydrogen_phosphate_net_change_g_p_per_step = -100;
    inputs.dihydrogen_phosphate_net_change_g_p_per_step = -100;
    const result = try calculate(inputs);

    try std.testing.expect(result.hydrogen_phosphate_floor_was_applied);
    try std.testing.expect(result.dihydrogen_phosphate_floor_was_applied);
    try std.testing.expectEqual(
        inputs.aqueous_phosphorus_floor_mol_p_per_m3,
        result.aqueous_hydrogen_phosphate_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        inputs.aqueous_phosphorus_floor_mol_p_per_m3,
        result.aqueous_dihydrogen_phosphate_mol_p_per_m3,
    );
}

test "surface phosphate reconstruction rejects invalid input and overflow" {
    var inputs = testInputs();
    inputs.litter_water_volume_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterPhosphateConcentrationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.hydrogen_phosphate_net_change_g_p_per_step = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterPhosphateConcentrationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.dihydrogen_phosphate_inventory_g_p = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterPhosphateConcentrationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.phosphorus_molar_mass_g_p_per_mol = std.math.floatMax(f64);
    inputs.litter_water_volume_m3 = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterPhosphateConcentration,
        calculate(inputs),
    );
}
