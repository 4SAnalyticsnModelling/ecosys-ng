const std = @import("std");

pub const PhosphateState = struct {
    hydrogen_phosphate_concentration_mol_p_per_m3: f64,
    hydrogen_phosphate_activity_mol_p_per_m3: f64,
    dihydrogen_phosphate_concentration_mol_p_per_m3: f64,
    dihydrogen_phosphate_activity_mol_p_per_m3: f64,
    hydrogen_activity_mol_per_m3: f64,
};

pub const Kinetics = struct {
    substrate_limit_fraction_per_step: f64,
    maximum_association_mol_p_per_m3_step: f64,
};

pub const Inputs = struct {
    state: PhosphateState,
    dihydrogen_phosphate_dissociation_mol_per_m3: f64,
    divalent_activity_coefficient: f64,
    kinetics: Kinetics,
};

pub const Result = struct {
    dissociation_limit_mol_p_per_m3_step: f64,
    association_limit_mol_p_per_m3_step: f64,
    equilibrium_hydrogen_phosphate_activity_mol_p_per_m3: f64,
    dihydrogen_phosphate_association_mol_p_per_m3_step: f64,
};

/// Direct source-order translation of SOLUTE.F lines 4924--4926.
///
/// Positive output associates HPO4 with H into H2PO4; negative output
/// dissociates H2PO4. The pure kernel operates on one surface-litter cell.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const state = inputs.state;
    const kinetics = inputs.kinetics;

    // SOLUTE.F 4924--4926. The static branch limits association with HPO4
    // alone; it does not include hydrogen concentration in that ceiling.
    const equilibrium_hydrogen_phosphate =
        inputs.dihydrogen_phosphate_dissociation_mol_per_m3 *
        state.dihydrogen_phosphate_activity_mol_p_per_m3 /
        state.hydrogen_activity_mol_per_m3;
    const dissociation_limit =
        kinetics.substrate_limit_fraction_per_step *
        state.dihydrogen_phosphate_concentration_mol_p_per_m3;
    const association_limit =
        kinetics.substrate_limit_fraction_per_step *
        state.hydrogen_phosphate_concentration_mol_p_per_m3;
    const maximum = kinetics.maximum_association_mol_p_per_m3_step;
    const association = @max(
        -maximum,
        -dissociation_limit,
        @min(
            maximum,
            association_limit,
            (state.hydrogen_phosphate_activity_mol_p_per_m3 -
                equilibrium_hydrogen_phosphate) /
                inputs.divalent_activity_coefficient,
        ),
    );

    const result: Result = .{
        .dissociation_limit_mol_p_per_m3_step = dissociation_limit,
        .association_limit_mol_p_per_m3_step = association_limit,
        .equilibrium_hydrogen_phosphate_activity_mol_p_per_m3 = equilibrium_hydrogen_phosphate,
        .dihydrogen_phosphate_association_mol_p_per_m3_step = association,
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    try validateNonnegativeStruct(inputs.state);
    try validateNonnegativeStruct(inputs.kinetics);
    if (inputs.state.hydrogen_activity_mol_per_m3 <= 0 or
        !std.math.isFinite(
            inputs.dihydrogen_phosphate_dissociation_mol_per_m3,
        ) or
        inputs.dihydrogen_phosphate_dissociation_mol_per_m3 <= 0 or
        !std.math.isFinite(inputs.divalent_activity_coefficient) or
        inputs.divalent_activity_coefficient <= 0 or
        inputs.kinetics.substrate_limit_fraction_per_step > 1)
    {
        return error.InvalidSurfaceLitterStaticPhosphateEquilibriumInput;
    }
}

fn validateNonnegativeStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        const field_value = @field(value, field.name);
        if (!std.math.isFinite(field_value) or field_value < 0)
            return error.InvalidSurfaceLitterStaticPhosphateEquilibriumInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceLitterStaticPhosphateEquilibriumResult;
    }
}

fn testInputs() Inputs {
    return .{
        .state = .{
            .hydrogen_phosphate_concentration_mol_p_per_m3 = 0.8,
            .hydrogen_phosphate_activity_mol_p_per_m3 = 0.48,
            .dihydrogen_phosphate_concentration_mol_p_per_m3 = 0.9,
            .dihydrogen_phosphate_activity_mol_p_per_m3 = 0.72,
            .hydrogen_activity_mol_per_m3 = 0.4,
        },
        .dihydrogen_phosphate_dissociation_mol_per_m3 = 0.1,
        .divalent_activity_coefficient = 0.6,
        .kinetics = .{
            .substrate_limit_fraction_per_step = 0.2,
            .maximum_association_mol_p_per_m3_step = 0.15,
        },
    };
}

test "SOLUTE static phosphate equilibrium preserves every source expression" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const state = inputs.state;
    const kinetics = inputs.kinetics;
    const expected_equilibrium =
        inputs.dihydrogen_phosphate_dissociation_mol_per_m3 *
        state.dihydrogen_phosphate_activity_mol_p_per_m3 /
        state.hydrogen_activity_mol_per_m3;
    const expected_dissociation =
        kinetics.substrate_limit_fraction_per_step *
        state.dihydrogen_phosphate_concentration_mol_p_per_m3;
    const expected_association =
        kinetics.substrate_limit_fraction_per_step *
        state.hydrogen_phosphate_concentration_mol_p_per_m3;
    const expected_rate = @max(
        -kinetics.maximum_association_mol_p_per_m3_step,
        -expected_dissociation,
        @min(
            kinetics.maximum_association_mol_p_per_m3_step,
            expected_association,
            (state.hydrogen_phosphate_activity_mol_p_per_m3 -
                expected_equilibrium) /
                inputs.divalent_activity_coefficient,
        ),
    );

    try std.testing.expectEqual(
        expected_equilibrium,
        result.equilibrium_hydrogen_phosphate_activity_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        expected_dissociation,
        result.dissociation_limit_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        expected_association,
        result.association_limit_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        expected_rate,
        result.dihydrogen_phosphate_association_mol_p_per_m3_step,
    );
}

test "static phosphate association ceiling uses hydrogen phosphate alone" {
    var inputs = testInputs();
    inputs.state.hydrogen_phosphate_concentration_mol_p_per_m3 = 0.05;
    inputs.state.hydrogen_phosphate_activity_mol_p_per_m3 = 10;
    inputs.state.hydrogen_activity_mol_per_m3 = 1.0e-12;
    inputs.dihydrogen_phosphate_dissociation_mol_per_m3 = 1.0e-14;
    inputs.state.dihydrogen_phosphate_activity_mol_p_per_m3 = 0.01;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectEqual(
        inputs.kinetics.substrate_limit_fraction_per_step *
            inputs.state.hydrogen_phosphate_concentration_mol_p_per_m3,
        result.association_limit_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        result.association_limit_mol_p_per_m3_step,
        result.dihydrogen_phosphate_association_mol_p_per_m3_step,
    );
}

test "static phosphate dissociation respects inventory ceiling" {
    var inputs = testInputs();
    inputs.state.hydrogen_phosphate_activity_mol_p_per_m3 = 0;
    inputs.state.dihydrogen_phosphate_concentration_mol_p_per_m3 = 0.1;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectEqual(
        -result.dissociation_limit_mol_p_per_m3_step,
        result.dihydrogen_phosphate_association_mol_p_per_m3_step,
    );
}

test "static phosphate equilibrium rejects invalid input and overflow" {
    var inputs = testInputs();
    inputs.state.hydrogen_activity_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticPhosphateEquilibriumInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.state.hydrogen_phosphate_concentration_mol_p_per_m3 = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticPhosphateEquilibriumInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.divalent_activity_coefficient = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticPhosphateEquilibriumInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.dihydrogen_phosphate_dissociation_mol_per_m3 =
        std.math.floatMax(f64);
    inputs.state.dihydrogen_phosphate_activity_mol_p_per_m3 =
        std.math.floatMax(f64);
    inputs.state.hydrogen_activity_mol_per_m3 =
        std.math.floatMin(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterStaticPhosphateEquilibriumResult,
        calculateSourceOrder(inputs),
    );
}
