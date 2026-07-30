const std = @import("std");

pub const Concentrations = struct {
    ammonium_mol_n_per_m3: f64,
    ammonia_mol_n_per_m3: f64,
    hydrogen_mol_per_m3: f64,
};

pub const Activities = struct {
    ammonium_mol_n_per_m3: f64,
    ammonia_mol_n_per_m3: f64,
    hydrogen_mol_per_m3: f64,
};

pub const Kinetics = struct {
    substrate_limit_fraction_per_step: f64,
    maximum_association_mol_n_per_m3_step: f64,
};

pub const Inputs = struct {
    concentrations: Concentrations,
    activities: Activities,
    ammonium_dissociation_mol_per_m3: f64,
    kinetics: Kinetics,
};

pub const Result = struct {
    dissociation_limit_mol_n_per_m3_step: f64,
    association_limit_mol_n_per_m3_step: f64,
    equilibrium_ammonia_activity_mol_n_per_m3: f64,
    association_mol_n_per_m3_step: f64,
};

/// Direct source-order translation of SOLUTE.F lines 4906--4910.
///
/// Positive output associates NH3 + H into NH4; negative output dissociates
/// NH4. The pure kernel operates on one surface-litter cell.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const concentration = inputs.concentrations;
    const activity = inputs.activities;
    const kinetics = inputs.kinetics;

    // SOLUTE.F 4906--4910.
    const dissociation_limit =
        kinetics.substrate_limit_fraction_per_step *
        concentration.ammonium_mol_n_per_m3;
    const association_limit =
        kinetics.substrate_limit_fraction_per_step *
        @min(
            concentration.hydrogen_mol_per_m3,
            concentration.ammonia_mol_n_per_m3,
        );
    const equilibrium_ammonia =
        inputs.ammonium_dissociation_mol_per_m3 *
        activity.ammonium_mol_n_per_m3 /
        activity.hydrogen_mol_per_m3;
    const maximum = kinetics.maximum_association_mol_n_per_m3_step;
    const association = @max(
        -maximum,
        -dissociation_limit,
        @min(
            maximum,
            association_limit,
            activity.ammonia_mol_n_per_m3 - equilibrium_ammonia,
        ),
    );

    const result: Result = .{
        .dissociation_limit_mol_n_per_m3_step = dissociation_limit,
        .association_limit_mol_n_per_m3_step = association_limit,
        .equilibrium_ammonia_activity_mol_n_per_m3 = equilibrium_ammonia,
        .association_mol_n_per_m3_step = association,
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    try validateNonnegativeStruct(inputs.concentrations);
    try validateNonnegativeStruct(inputs.activities);
    try validateNonnegativeStruct(inputs.kinetics);
    if (inputs.activities.hydrogen_mol_per_m3 <= 0 or
        !std.math.isFinite(inputs.ammonium_dissociation_mol_per_m3) or
        inputs.ammonium_dissociation_mol_per_m3 <= 0 or
        inputs.kinetics.substrate_limit_fraction_per_step > 1)
    {
        return error.InvalidSurfaceLitterStaticAmmoniumEquilibriumInput;
    }
}

fn validateNonnegativeStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        const field_value = @field(value, field.name);
        if (!std.math.isFinite(field_value) or field_value < 0)
            return error.InvalidSurfaceLitterStaticAmmoniumEquilibriumInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceLitterStaticAmmoniumEquilibriumResult;
    }
}

fn testInputs() Inputs {
    return .{
        .concentrations = .{
            .ammonium_mol_n_per_m3 = 0.8,
            .ammonia_mol_n_per_m3 = 0.6,
            .hydrogen_mol_per_m3 = 0.5,
        },
        .activities = .{
            .ammonium_mol_n_per_m3 = 0.64,
            .ammonia_mol_n_per_m3 = 0.48,
            .hydrogen_mol_per_m3 = 0.4,
        },
        .ammonium_dissociation_mol_per_m3 = 0.1,
        .kinetics = .{
            .substrate_limit_fraction_per_step = 0.2,
            .maximum_association_mol_n_per_m3_step = 0.15,
        },
    };
}

test "SOLUTE static ammonium equilibrium preserves every source expression" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const concentration = inputs.concentrations;
    const activity = inputs.activities;
    const kinetics = inputs.kinetics;
    const expected_dissociation =
        kinetics.substrate_limit_fraction_per_step *
        concentration.ammonium_mol_n_per_m3;
    const expected_association =
        kinetics.substrate_limit_fraction_per_step *
        @min(
            concentration.hydrogen_mol_per_m3,
            concentration.ammonia_mol_n_per_m3,
        );
    const expected_equilibrium =
        inputs.ammonium_dissociation_mol_per_m3 *
        activity.ammonium_mol_n_per_m3 /
        activity.hydrogen_mol_per_m3;
    const expected_rate = @max(
        -kinetics.maximum_association_mol_n_per_m3_step,
        -expected_dissociation,
        @min(
            kinetics.maximum_association_mol_n_per_m3_step,
            expected_association,
            activity.ammonia_mol_n_per_m3 - expected_equilibrium,
        ),
    );

    try std.testing.expectEqual(
        expected_dissociation,
        result.dissociation_limit_mol_n_per_m3_step,
    );
    try std.testing.expectEqual(
        expected_association,
        result.association_limit_mol_n_per_m3_step,
    );
    try std.testing.expectEqual(
        expected_equilibrium,
        result.equilibrium_ammonia_activity_mol_n_per_m3,
    );
    try std.testing.expectEqual(
        expected_rate,
        result.association_mol_n_per_m3_step,
    );
}

test "static ammonium association respects substrate ceiling" {
    const result = try calculateSourceOrder(testInputs());

    try std.testing.expectEqual(
        result.association_limit_mol_n_per_m3_step,
        result.association_mol_n_per_m3_step,
    );
}

test "static ammonium dissociation respects kinetic and inventory ceilings" {
    var inputs = testInputs();
    inputs.activities.ammonia_mol_n_per_m3 = 0;
    var result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        -inputs.kinetics.maximum_association_mol_n_per_m3_step,
        result.association_mol_n_per_m3_step,
    );

    inputs.concentrations.ammonium_mol_n_per_m3 = 0.1;
    result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        -result.dissociation_limit_mol_n_per_m3_step,
        result.association_mol_n_per_m3_step,
    );
}

test "static ammonium equilibrium accepts zero available substrates" {
    var inputs = testInputs();
    inputs.concentrations.ammonia_mol_n_per_m3 = 0;
    inputs.concentrations.hydrogen_mol_per_m3 = 0;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectEqual(
        @as(f64, 0),
        result.association_limit_mol_n_per_m3_step,
    );
    try std.testing.expect(result.association_mol_n_per_m3_step <= 0);
}

test "static ammonium equilibrium rejects invalid input and overflow" {
    var inputs = testInputs();
    inputs.activities.hydrogen_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticAmmoniumEquilibriumInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.concentrations.ammonium_mol_n_per_m3 = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticAmmoniumEquilibriumInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.kinetics.substrate_limit_fraction_per_step = 1.1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticAmmoniumEquilibriumInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.ammonium_dissociation_mol_per_m3 = std.math.floatMax(f64);
    inputs.activities.ammonium_mol_n_per_m3 = std.math.floatMax(f64);
    inputs.activities.hydrogen_mol_per_m3 = std.math.floatMin(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterStaticAmmoniumEquilibriumResult,
        calculateSourceOrder(inputs),
    );
}
