const std = @import("std");

pub const SaltEquilibriumMode = enum {
    static_concentrations,
    dynamic_equilibria,
};

pub const InitialIonTransformations = struct {
    hydrogen_mol_per_m3_step: f64,
    hydroxide_mol_per_m3_step: f64,
    aluminum_mol_per_m3_step: f64,
    iron_mol_per_m3_step: f64,
    calcium_mol_per_m3_step: f64,
};

pub const AssociationRates = struct {
    carbon_dioxide_mol_c_per_m3_step: f64,
    bicarbonate_mol_c_per_m3_step: f64,
    dihydrogen_phosphate_mol_p_per_m3_step: f64,
    ammonium_mol_n_per_m3_step: f64,
};

pub const PhosphateMineralRates = struct {
    aluminum_phosphate_mol_mineral_per_m3_step: f64,
    iron_phosphate_mol_mineral_per_m3_step: f64,
    dicalcium_phosphate_mol_mineral_per_m3_step: f64,
    hydroxyapatite_mol_mineral_per_m3_step: f64,
    monocalcium_phosphate_mol_mineral_per_m3_step: f64,
};

pub const Inputs = struct {
    mode: SaltEquilibriumMode,
    initial: InitialIonTransformations,
    external_hydrogen_mol_per_m3_step: f64,
    carboxyl_hydrogen_exchange_mol_per_Mg_step: f64,
    litter_mass_per_water_volume_Mg_per_m3: f64,
    association: AssociationRates,
    phosphate_minerals: PhosphateMineralRates,
};

pub const DynamicTransformations = struct {
    hydrogen_mol_per_m3_step: f64,
    hydroxide_mol_per_m3_step: f64,
    aluminum_mol_per_m3_step: f64,
    iron_mol_per_m3_step: f64,
    calcium_mol_per_m3_step: f64,
    carbonate_mol_c_per_m3_step: f64,
    bicarbonate_mol_c_per_m3_step: f64,
    carbon_dioxide_mol_c_per_m3_step: f64,
    sulfate_mol_s_per_m3_step: f64,
    water_mol_per_m3_step: f64,
};

pub const Result = union(enum) {
    static_concentrations_skipped,
    dynamic: DynamicTransformations,
};

/// Direct source-order translation of SOLUTE.F lines 5011--5023.
///
/// The source enters this block only for dynamically solved salt equilibria.
/// This pure one-cell assembly does not apply the subsequent mineral updates.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);

    // SOLUTE.F 5011.
    if (inputs.mode == .static_concentrations)
        return .static_concentrations_skipped;

    const initial = inputs.initial;
    const association = inputs.association;
    const minerals = inputs.phosphate_minerals;

    // SOLUTE.F 5012--5023, preserved in assignment order.
    const hydrogen =
        initial.hydrogen_mol_per_m3_step +
        inputs.external_hydrogen_mol_per_m3_step -
        inputs.carboxyl_hydrogen_exchange_mol_per_Mg_step *
            inputs.litter_mass_per_water_volume_Mg_per_m3 -
        association.carbon_dioxide_mol_c_per_m3_step -
        association.bicarbonate_mol_c_per_m3_step +
        2.0 *
            (minerals.aluminum_phosphate_mol_mineral_per_m3_step +
                minerals.iron_phosphate_mol_mineral_per_m3_step) +
        minerals.dicalcium_phosphate_mol_mineral_per_m3_step +
        6.0 * minerals.hydroxyapatite_mol_mineral_per_m3_step -
        association.dihydrogen_phosphate_mol_p_per_m3_step -
        association.ammonium_mol_n_per_m3_step;
    const hydroxide =
        initial.hydroxide_mol_per_m3_step -
        minerals.hydroxyapatite_mol_mineral_per_m3_step;
    const aluminum =
        initial.aluminum_mol_per_m3_step -
        minerals.aluminum_phosphate_mol_mineral_per_m3_step;
    const iron =
        initial.iron_mol_per_m3_step -
        minerals.iron_phosphate_mol_mineral_per_m3_step;
    const calcium =
        initial.calcium_mol_per_m3_step -
        minerals.dicalcium_phosphate_mol_mineral_per_m3_step -
        5.0 * minerals.hydroxyapatite_mol_mineral_per_m3_step -
        minerals.monocalcium_phosphate_mol_mineral_per_m3_step;
    const carbonate = -association.bicarbonate_mol_c_per_m3_step;
    const bicarbonate =
        -association.carbon_dioxide_mol_c_per_m3_step +
        association.bicarbonate_mol_c_per_m3_step;
    const carbon_dioxide = association.carbon_dioxide_mol_c_per_m3_step;
    const sulfate = 0.0;
    const water = association.carbon_dioxide_mol_c_per_m3_step;

    const result: DynamicTransformations = .{
        .hydrogen_mol_per_m3_step = hydrogen,
        .hydroxide_mol_per_m3_step = hydroxide,
        .aluminum_mol_per_m3_step = aluminum,
        .iron_mol_per_m3_step = iron,
        .calcium_mol_per_m3_step = calcium,
        .carbonate_mol_c_per_m3_step = carbonate,
        .bicarbonate_mol_c_per_m3_step = bicarbonate,
        .carbon_dioxide_mol_c_per_m3_step = carbon_dioxide,
        .sulfate_mol_s_per_m3_step = sulfate,
        .water_mol_per_m3_step = water,
    };
    try validateResult(result);
    return .{ .dynamic = result };
}

fn validateInputs(inputs: Inputs) !void {
    try validateFiniteStruct(inputs.initial);
    try validateFiniteStruct(inputs.association);
    try validateFiniteStruct(inputs.phosphate_minerals);
    inline for (.{
        inputs.external_hydrogen_mol_per_m3_step,
        inputs.carboxyl_hydrogen_exchange_mol_per_Mg_step,
        inputs.litter_mass_per_water_volume_Mg_per_m3,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.InvalidSurfaceLitterDynamicSaltAssemblyInput;
    }
    if (inputs.litter_mass_per_water_volume_Mg_per_m3 <= 0)
        return error.InvalidSurfaceLitterDynamicSaltAssemblyInput;
}

fn validateFiniteStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(value, field.name)))
            return error.InvalidSurfaceLitterDynamicSaltAssemblyInput;
    }
}

fn validateResult(result: DynamicTransformations) !void {
    inline for (@typeInfo(DynamicTransformations).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceLitterDynamicSaltAssemblyResult;
    }
}

fn testInputs() Inputs {
    return .{
        .mode = .dynamic_equilibria,
        .initial = .{
            .hydrogen_mol_per_m3_step = 0.10,
            .hydroxide_mol_per_m3_step = 0.20,
            .aluminum_mol_per_m3_step = 0.30,
            .iron_mol_per_m3_step = 0.40,
            .calcium_mol_per_m3_step = 0.50,
        },
        .external_hydrogen_mol_per_m3_step = 0.60,
        .carboxyl_hydrogen_exchange_mol_per_Mg_step = 0.07,
        .litter_mass_per_water_volume_Mg_per_m3 = 2,
        .association = .{
            .carbon_dioxide_mol_c_per_m3_step = 0.08,
            .bicarbonate_mol_c_per_m3_step = 0.09,
            .dihydrogen_phosphate_mol_p_per_m3_step = 0.10,
            .ammonium_mol_n_per_m3_step = 0.11,
        },
        .phosphate_minerals = .{
            .aluminum_phosphate_mol_mineral_per_m3_step = 0.012,
            .iron_phosphate_mol_mineral_per_m3_step = 0.013,
            .dicalcium_phosphate_mol_mineral_per_m3_step = 0.014,
            .hydroxyapatite_mol_mineral_per_m3_step = 0.015,
            .monocalcium_phosphate_mol_mineral_per_m3_step = 0.016,
        },
    };
}

test "SOLUTE dynamic salt assembly preserves every source equation" {
    const inputs = testInputs();
    const result = (try calculateSourceOrder(inputs)).dynamic;
    const initial = inputs.initial;
    const association = inputs.association;
    const minerals = inputs.phosphate_minerals;

    try std.testing.expectEqual(
        initial.hydrogen_mol_per_m3_step +
            inputs.external_hydrogen_mol_per_m3_step -
            inputs.carboxyl_hydrogen_exchange_mol_per_Mg_step *
                inputs.litter_mass_per_water_volume_Mg_per_m3 -
            association.carbon_dioxide_mol_c_per_m3_step -
            association.bicarbonate_mol_c_per_m3_step +
            2.0 *
                (minerals.aluminum_phosphate_mol_mineral_per_m3_step +
                    minerals.iron_phosphate_mol_mineral_per_m3_step) +
            minerals.dicalcium_phosphate_mol_mineral_per_m3_step +
            6.0 * minerals.hydroxyapatite_mol_mineral_per_m3_step -
            association.dihydrogen_phosphate_mol_p_per_m3_step -
            association.ammonium_mol_n_per_m3_step,
        result.hydrogen_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        initial.hydroxide_mol_per_m3_step -
            minerals.hydroxyapatite_mol_mineral_per_m3_step,
        result.hydroxide_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        initial.calcium_mol_per_m3_step -
            minerals.dicalcium_phosphate_mol_mineral_per_m3_step -
            5.0 * minerals.hydroxyapatite_mol_mineral_per_m3_step -
            minerals.monocalcium_phosphate_mol_mineral_per_m3_step,
        result.calcium_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        association.carbon_dioxide_mol_c_per_m3_step,
        result.water_mol_per_m3_step,
    );
}

test "dynamic salt carbon transformations close exactly" {
    const result = (try calculateSourceOrder(testInputs())).dynamic;
    const carbon_sum =
        result.carbonate_mol_c_per_m3_step +
        result.bicarbonate_mol_c_per_m3_step +
        result.carbon_dioxide_mol_c_per_m3_step;

    try std.testing.expectApproxEqAbs(@as(f64, 0), carbon_sum, 1.0e-16);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.sulfate_mol_s_per_m3_step,
    );
}

test "static salt concentrations skip dynamic transformation assembly" {
    var inputs = testInputs();
    inputs.mode = .static_concentrations;

    try std.testing.expectEqual(
        Result.static_concentrations_skipped,
        try calculateSourceOrder(inputs),
    );
}

test "dynamic salt assembly accepts signed reversible rates" {
    var inputs = testInputs();
    inline for (@typeInfo(AssociationRates).@"struct".fields) |field|
        @field(inputs.association, field.name) *= -1;
    inline for (@typeInfo(PhosphateMineralRates).@"struct".fields) |field|
        @field(inputs.phosphate_minerals, field.name) *= -1;
    const result = (try calculateSourceOrder(inputs)).dynamic;

    try std.testing.expect(result.carbon_dioxide_mol_c_per_m3_step < 0);
    try std.testing.expect(result.water_mol_per_m3_step < 0);
}

test "dynamic salt assembly rejects invalid input and overflow" {
    var inputs = testInputs();
    inputs.litter_mass_per_water_volume_Mg_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterDynamicSaltAssemblyInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.association.ammonium_mol_n_per_m3_step = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterDynamicSaltAssemblyInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.litter_mass_per_water_volume_Mg_per_m3 =
        std.math.floatMax(f64);
    inputs.carboxyl_hydrogen_exchange_mol_per_Mg_step =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterDynamicSaltAssemblyResult,
        calculateSourceOrder(inputs),
    );
}
