const std = @import("std");

pub const SpeciesState = struct {
    concentration_mol_per_m3: f64,
    activity_mol_per_m3: f64,
};

pub const State = struct {
    hydrogen_concentration_mol_per_m3: f64,
    hydrogen_activity_mol_per_m3: f64,
    ammonium: SpeciesState,
    ammonia: SpeciesState,
    carbonate: SpeciesState,
    bicarbonate: SpeciesState,
    carbon_dioxide: SpeciesState,
    hydrogen_phosphate: SpeciesState,
    dihydrogen_phosphate: SpeciesState,
};

pub const EquilibriumConstants = struct {
    ammonium_dissociation_mol_per_m3: f64,
    bicarbonate_dissociation_mol_per_m3: f64,
    carbon_dioxide_dissociation_mol_per_m3: f64,
    dihydrogen_phosphate_dissociation_mol_per_m3: f64,
};

pub const Kinetics = struct {
    ammonium_substrate_limit_fraction_per_step: f64,
    general_substrate_limit_fraction_per_step: f64,
    maximum_ammonium_association_mol_per_m3_step: f64,
    maximum_general_association_mol_per_m3_step: f64,
};

pub const Inputs = struct {
    state: State,
    equilibrium: EquilibriumConstants,
    kinetics: Kinetics,
    monovalent_activity_coefficient: f64,
    divalent_activity_coefficient: f64,
    /// SOLUTE line 4538 uses carboxyl `XMIN`, not bicarbonate `XMINN`.
    retained_carboxyl_substrate_limit_source_units: f64,
};

pub const ReactionDiagnostics = struct {
    calculated_dissociation_limit_mol_per_m3_step: f64,
    association_limit_mol_per_m3_step: f64,
    equilibrium_free_activity_mol_per_m3: f64,
    association_mol_per_m3_step: f64,
};

pub const Result = struct {
    ammonium: ReactionDiagnostics,
    bicarbonate: ReactionDiagnostics,
    carbon_dioxide: ReactionDiagnostics,
    dihydrogen_phosphate: ReactionDiagnostics,
    hydrogen_after_ammonium_mol_per_m3: f64,
    hydrogen_activity_after_ammonium_mol_per_m3: f64,
    bicarbonate_source_dissociation_bound: f64,
};

/// Direct source-order translation of SOLUTE.F lines 4521--4565.
///
/// Positive rates associate the named conjugate base with H. The NH4 rate
/// updates H activity before the following three equilibrium targets, exactly
/// as in the surface source block.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const state = inputs.state;
    const kinetics = inputs.kinetics;
    const general_fraction =
        kinetics.general_substrate_limit_fraction_per_step;

    // SOLUTE.F 4523--4529.
    const ammonium_dissociation =
        kinetics.ammonium_substrate_limit_fraction_per_step *
        state.ammonium.concentration_mol_per_m3;
    const ammonium_association =
        kinetics.ammonium_substrate_limit_fraction_per_step *
        @min(
            state.hydrogen_concentration_mol_per_m3,
            state.ammonia.concentration_mol_per_m3,
        );
    const ammonia_equilibrium =
        inputs.equilibrium.ammonium_dissociation_mol_per_m3 *
        state.ammonium.activity_mol_per_m3 /
        state.hydrogen_activity_mol_per_m3;
    const ammonium_rate = bounded(
        state.ammonia.activity_mol_per_m3 - ammonia_equilibrium,
        ammonium_dissociation,
        ammonium_association,
        kinetics.maximum_ammonium_association_mol_per_m3_step,
    );
    const hydrogen_after =
        state.hydrogen_concentration_mol_per_m3 - ammonium_rate;
    const hydrogen_activity_after =
        hydrogen_after * inputs.monovalent_activity_coefficient;
    if (!std.math.isFinite(hydrogen_activity_after) or
        hydrogen_activity_after <= 0)
    {
        return error.InvalidSurfaceLitterPostAmmoniumHydrogenActivity;
    }

    // SOLUTE.F 4535--4539. The calculated bicarbonate dissociation limit is
    // retained diagnostically, but line 4538 reads the preceding `XMIN`.
    const bicarbonate_dissociation =
        general_fraction * state.bicarbonate.concentration_mol_per_m3;
    const bicarbonate_association =
        general_fraction *
        @min(
            state.hydrogen_concentration_mol_per_m3,
            state.carbonate.concentration_mol_per_m3,
        );
    const carbonate_equilibrium =
        inputs.equilibrium.bicarbonate_dissociation_mol_per_m3 *
        state.bicarbonate.activity_mol_per_m3 /
        hydrogen_activity_after;
    const bicarbonate_rate = bounded(
        (state.carbonate.activity_mol_per_m3 - carbonate_equilibrium) /
            inputs.divalent_activity_coefficient,
        inputs.retained_carboxyl_substrate_limit_source_units,
        bicarbonate_association,
        kinetics.maximum_general_association_mol_per_m3_step,
    );

    // SOLUTE.F 4543--4547.
    const carbon_dioxide_dissociation =
        general_fraction * state.carbon_dioxide.concentration_mol_per_m3;
    const carbon_dioxide_association =
        general_fraction *
        @min(
            state.hydrogen_concentration_mol_per_m3,
            state.bicarbonate.concentration_mol_per_m3,
        );
    const bicarbonate_equilibrium =
        inputs.equilibrium.carbon_dioxide_dissociation_mol_per_m3 *
        state.carbon_dioxide.activity_mol_per_m3 /
        hydrogen_activity_after;
    const carbon_dioxide_rate = bounded(
        (state.bicarbonate.activity_mol_per_m3 - bicarbonate_equilibrium) /
            inputs.monovalent_activity_coefficient,
        carbon_dioxide_dissociation,
        carbon_dioxide_association,
        kinetics.maximum_general_association_mol_per_m3_step,
    );

    // SOLUTE.F 4561--4565.
    const phosphate_dissociation =
        general_fraction *
        state.dihydrogen_phosphate.concentration_mol_per_m3;
    const phosphate_association =
        general_fraction *
        @min(
            state.hydrogen_concentration_mol_per_m3,
            state.hydrogen_phosphate.concentration_mol_per_m3,
        );
    const hydrogen_phosphate_equilibrium =
        inputs.equilibrium.dihydrogen_phosphate_dissociation_mol_per_m3 *
        state.dihydrogen_phosphate.activity_mol_per_m3 /
        hydrogen_activity_after;
    const phosphate_rate = bounded(
        (state.hydrogen_phosphate.activity_mol_per_m3 -
            hydrogen_phosphate_equilibrium) /
            inputs.divalent_activity_coefficient,
        phosphate_dissociation,
        phosphate_association,
        kinetics.maximum_general_association_mol_per_m3_step,
    );

    const result: Result = .{
        .ammonium = .{
            .calculated_dissociation_limit_mol_per_m3_step = ammonium_dissociation,
            .association_limit_mol_per_m3_step = ammonium_association,
            .equilibrium_free_activity_mol_per_m3 = ammonia_equilibrium,
            .association_mol_per_m3_step = ammonium_rate,
        },
        .bicarbonate = .{
            .calculated_dissociation_limit_mol_per_m3_step = bicarbonate_dissociation,
            .association_limit_mol_per_m3_step = bicarbonate_association,
            .equilibrium_free_activity_mol_per_m3 = carbonate_equilibrium,
            .association_mol_per_m3_step = bicarbonate_rate,
        },
        .carbon_dioxide = .{
            .calculated_dissociation_limit_mol_per_m3_step = carbon_dioxide_dissociation,
            .association_limit_mol_per_m3_step = carbon_dioxide_association,
            .equilibrium_free_activity_mol_per_m3 = bicarbonate_equilibrium,
            .association_mol_per_m3_step = carbon_dioxide_rate,
        },
        .dihydrogen_phosphate = .{
            .calculated_dissociation_limit_mol_per_m3_step = phosphate_dissociation,
            .association_limit_mol_per_m3_step = phosphate_association,
            .equilibrium_free_activity_mol_per_m3 = hydrogen_phosphate_equilibrium,
            .association_mol_per_m3_step = phosphate_rate,
        },
        .hydrogen_after_ammonium_mol_per_m3 = hydrogen_after,
        .hydrogen_activity_after_ammonium_mol_per_m3 = hydrogen_activity_after,
        .bicarbonate_source_dissociation_bound = inputs.retained_carboxyl_substrate_limit_source_units,
    };
    try validateResult(result);
    return result;
}

fn bounded(
    driving: f64,
    dissociation_limit: f64,
    association_limit: f64,
    maximum: f64,
) f64 {
    return @max(
        -maximum,
        -dissociation_limit,
        @min(maximum, association_limit, driving),
    );
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (field.type == SpeciesState) {
            const species_state = @field(inputs.state, field.name);
            inline for (@typeInfo(SpeciesState).@"struct".fields) |member| {
                const value = @field(species_state, member.name);
                if (!std.math.isFinite(value) or value < 0)
                    return error.InvalidSurfaceLitterAcidBaseReactionInput;
            }
        } else {
            const value = @field(inputs.state, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidSurfaceLitterAcidBaseReactionInput;
        }
    }
    inline for (@typeInfo(EquilibriumConstants).@"struct".fields) |field| {
        const value = @field(inputs.equilibrium, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterAcidBaseReactionInput;
    }
    inline for (@typeInfo(Kinetics).@"struct".fields) |field| {
        const value = @field(inputs.kinetics, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterAcidBaseReactionInput;
    }
    if (inputs.state.hydrogen_activity_mol_per_m3 <= 0 or
        !std.math.isFinite(inputs.monovalent_activity_coefficient) or
        inputs.monovalent_activity_coefficient <= 0 or
        !std.math.isFinite(inputs.divalent_activity_coefficient) or
        inputs.divalent_activity_coefficient <= 0 or
        !std.math.isFinite(inputs.retained_carboxyl_substrate_limit_source_units) or
        inputs.retained_carboxyl_substrate_limit_source_units < 0 or
        inputs.kinetics.ammonium_substrate_limit_fraction_per_step > 1 or
        inputs.kinetics.general_substrate_limit_fraction_per_step > 1)
    {
        return error.InvalidSurfaceLitterAcidBaseReactionInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (.{
        result.ammonium,
        result.bicarbonate,
        result.carbon_dioxide,
        result.dihydrogen_phosphate,
    }) |reaction| {
        inline for (@typeInfo(ReactionDiagnostics).@"struct".fields) |field| {
            if (!std.math.isFinite(@field(reaction, field.name)))
                return error.NonFiniteSurfaceLitterAcidBaseReactionResult;
        }
    }
    inline for (.{
        result.hydrogen_after_ammonium_mol_per_m3,
        result.hydrogen_activity_after_ammonium_mol_per_m3,
        result.bicarbonate_source_dissociation_bound,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterAcidBaseReactionResult;
    }
}

fn species(concentration: f64, activity: f64) SpeciesState {
    return .{
        .concentration_mol_per_m3 = concentration,
        .activity_mol_per_m3 = activity,
    };
}

fn testInputs() Inputs {
    return .{
        .state = .{
            .hydrogen_concentration_mol_per_m3 = 0.8,
            .hydrogen_activity_mol_per_m3 = 0.64,
            .ammonium = species(0.6, 0.48),
            .ammonia = species(0.5, 0.5),
            .carbonate = species(0.4, 0.24),
            .bicarbonate = species(0.3, 0.24),
            .carbon_dioxide = species(0.2, 0.2),
            .hydrogen_phosphate = species(0.35, 0.21),
            .dihydrogen_phosphate = species(0.25, 0.2),
        },
        .equilibrium = .{
            .ammonium_dissociation_mol_per_m3 = 0.2,
            .bicarbonate_dissociation_mol_per_m3 = 0.15,
            .carbon_dioxide_dissociation_mol_per_m3 = 0.1,
            .dihydrogen_phosphate_dissociation_mol_per_m3 = 0.12,
        },
        .kinetics = .{
            .ammonium_substrate_limit_fraction_per_step = 0.4,
            .general_substrate_limit_fraction_per_step = 0.3,
            .maximum_ammonium_association_mol_per_m3_step = 0.2,
            .maximum_general_association_mol_per_m3_step = 0.1,
        },
        .monovalent_activity_coefficient = 0.8,
        .divalent_activity_coefficient = 0.6,
        .retained_carboxyl_substrate_limit_source_units = 0.07,
    };
}

test "SOLUTE surface acid-base reactions preserve every source expression" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const state = inputs.state;
    const k = inputs.kinetics;
    const equilibrium = inputs.equilibrium;

    const ammonium_dissociation =
        k.ammonium_substrate_limit_fraction_per_step *
        state.ammonium.concentration_mol_per_m3;
    const ammonium_association =
        k.ammonium_substrate_limit_fraction_per_step *
        @min(
            state.hydrogen_concentration_mol_per_m3,
            state.ammonia.concentration_mol_per_m3,
        );
    const ammonia_equilibrium =
        equilibrium.ammonium_dissociation_mol_per_m3 *
        state.ammonium.activity_mol_per_m3 /
        state.hydrogen_activity_mol_per_m3;
    const ammonium_rate = @max(
        -k.maximum_ammonium_association_mol_per_m3_step,
        -ammonium_dissociation,
        @min(
            k.maximum_ammonium_association_mol_per_m3_step,
            ammonium_association,
            state.ammonia.activity_mol_per_m3 - ammonia_equilibrium,
        ),
    );
    const hydrogen_after =
        state.hydrogen_concentration_mol_per_m3 - ammonium_rate;
    const hydrogen_activity_after =
        hydrogen_after * inputs.monovalent_activity_coefficient;
    const carbonate_equilibrium =
        equilibrium.bicarbonate_dissociation_mol_per_m3 *
        state.bicarbonate.activity_mol_per_m3 /
        hydrogen_activity_after;

    try std.testing.expectEqual(
        ammonium_dissociation,
        result.ammonium.calculated_dissociation_limit_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        ammonium_association,
        result.ammonium.association_limit_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        ammonia_equilibrium,
        result.ammonium.equilibrium_free_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        ammonium_rate,
        result.ammonium.association_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        hydrogen_activity_after,
        result.hydrogen_activity_after_ammonium_mol_per_m3,
    );
    try std.testing.expectEqual(
        carbonate_equilibrium,
        result.bicarbonate.equilibrium_free_activity_mol_per_m3,
    );

    const bicarbonate_expected = bounded(
        (state.carbonate.activity_mol_per_m3 - carbonate_equilibrium) /
            inputs.divalent_activity_coefficient,
        inputs.retained_carboxyl_substrate_limit_source_units,
        k.general_substrate_limit_fraction_per_step *
            @min(
                state.hydrogen_concentration_mol_per_m3,
                state.carbonate.concentration_mol_per_m3,
            ),
        k.maximum_general_association_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        bicarbonate_expected,
        result.bicarbonate.association_mol_per_m3_step,
    );
}

test "surface ammonium association updates later equilibrium hydrogen activity" {
    var inputs = testInputs();
    inputs.state = .{
        .hydrogen_concentration_mol_per_m3 = 1,
        .hydrogen_activity_mol_per_m3 = 1,
        .ammonium = species(1, 1),
        .ammonia = species(1, 1),
        .carbonate = species(1, 1),
        .bicarbonate = species(1, 1),
        .carbon_dioxide = species(1, 1),
        .hydrogen_phosphate = species(1, 1),
        .dihydrogen_phosphate = species(1, 1),
    };
    inputs.equilibrium.ammonium_dissociation_mol_per_m3 = 0.5;
    inputs.equilibrium.bicarbonate_dissociation_mol_per_m3 = 1;
    inputs.kinetics.ammonium_substrate_limit_fraction_per_step = 0.5;
    inputs.kinetics.general_substrate_limit_fraction_per_step = 1;
    inputs.kinetics.maximum_ammonium_association_mol_per_m3_step = 1;
    inputs.kinetics.maximum_general_association_mol_per_m3_step = 1;
    inputs.monovalent_activity_coefficient = 1;
    inputs.divalent_activity_coefficient = 1;
    inputs.retained_carboxyl_substrate_limit_source_units = 1;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectEqual(
        @as(f64, 0.5),
        result.ammonium.association_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.5),
        result.hydrogen_activity_after_ammonium_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 2),
        result.bicarbonate.equilibrium_free_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, -1),
        result.bicarbonate.association_mol_per_m3_step,
    );
}

test "surface bicarbonate lower bound retains carboxyl XMIN" {
    var inputs = testInputs();
    inputs.state.hydrogen_concentration_mol_per_m3 = 1;
    inputs.state.hydrogen_activity_mol_per_m3 = 1;
    inputs.state.ammonium = species(1, 1);
    inputs.state.ammonia = species(1, 1);
    inputs.state.carbonate = species(0, 0);
    inputs.state.bicarbonate = species(0.01, 0.01);
    inputs.equilibrium.ammonium_dissociation_mol_per_m3 = 1;
    inputs.equilibrium.bicarbonate_dissociation_mol_per_m3 = 100;
    inputs.kinetics.ammonium_substrate_limit_fraction_per_step = 1;
    inputs.kinetics.general_substrate_limit_fraction_per_step = 0.1;
    inputs.kinetics.maximum_ammonium_association_mol_per_m3_step = 1;
    inputs.kinetics.maximum_general_association_mol_per_m3_step = 1;
    inputs.monovalent_activity_coefficient = 1;
    inputs.divalent_activity_coefficient = 1;
    inputs.retained_carboxyl_substrate_limit_source_units = 0.1;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectEqual(
        @as(f64, 0.001),
        result.bicarbonate.calculated_dissociation_limit_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, -0.1),
        result.bicarbonate.association_mol_per_m3_step,
    );
    try std.testing.expect(
        -result.bicarbonate.association_mol_per_m3_step >
            inputs.state.bicarbonate.concentration_mol_per_m3,
    );
}

test "surface acid-base reactions reject zero post-ammonium hydrogen" {
    var inputs = testInputs();
    inputs.state.hydrogen_concentration_mol_per_m3 = 1;
    inputs.state.hydrogen_activity_mol_per_m3 = 1;
    inputs.state.ammonium = species(1, 1);
    inputs.state.ammonia = species(1, 2);
    inputs.equilibrium.ammonium_dissociation_mol_per_m3 = 1;
    inputs.kinetics.ammonium_substrate_limit_fraction_per_step = 1;
    inputs.kinetics.maximum_ammonium_association_mol_per_m3_step = 1;
    inputs.monovalent_activity_coefficient = 1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterPostAmmoniumHydrogenActivity,
        calculateSourceOrder(inputs),
    );
}

test "surface acid-base reactions reject invalid input and overflow" {
    var inputs = testInputs();
    inputs.state.bicarbonate.concentration_mol_per_m3 = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterAcidBaseReactionInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.divalent_activity_coefficient = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterAcidBaseReactionInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.equilibrium.bicarbonate_dissociation_mol_per_m3 =
        std.math.floatMax(f64);
    inputs.state.bicarbonate.activity_mol_per_m3 =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterAcidBaseReactionResult,
        calculateSourceOrder(inputs),
    );
}
