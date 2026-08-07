const std = @import("std");

pub const AqueousState = struct {
    hydrogen_activity_mol_per_m3: f64,
    hydroxide_activity_mol_per_m3: f64,
    aluminum_concentration_mol_per_m3: f64,
    aluminum_activity_mol_per_m3: f64,
    iron_concentration_mol_per_m3: f64,
    iron_activity_mol_per_m3: f64,
    calcium_concentration_mol_per_m3: f64,
    calcium_activity_mol_per_m3: f64,
    hydrogen_phosphate_concentration_mol_p_per_m3: f64,
    hydrogen_phosphate_activity_mol_p_per_m3: f64,
    dihydrogen_phosphate_concentration_mol_p_per_m3: f64,
    dihydrogen_phosphate_activity_mol_p_per_m3: f64,
};

pub const MineralState = struct {
    aluminum_phosphate_mol_per_m3: f64,
    iron_phosphate_mol_per_m3: f64,
    dicalcium_phosphate_mol_per_m3: f64,
    hydroxyapatite_mol_per_m3: f64,
    monocalcium_phosphate_mol_per_m3: f64,
};

pub const SolubilityProducts = struct {
    aluminum_phosphate: f64,
    iron_phosphate: f64,
    dicalcium_phosphate: f64,
    hydroxyapatite: f64,
    monocalcium_phosphate: f64,
};

pub const DissociationConstants = struct {
    dihydrogen_phosphate_mol_per_m3: f64,
    hydrogen_phosphate_mol_per_m3: f64,
};

pub const Kinetics = struct {
    substrate_limit_fraction_per_step: f64,
    maximum_phosphate_precipitation_mol_per_m3_step: f64,
    maximum_apatite_precipitation_mol_per_m3_step: f64,
    maximum_monocalcium_dissolution_mol_per_m3_step: f64,
};

pub const Inputs = struct {
    aqueous: AqueousState,
    minerals: MineralState,
    solubility_products: SolubilityProducts,
    dissociation: DissociationConstants,
    kinetics: Kinetics,
    monovalent_activity_coefficient: f64,
    divalent_activity_coefficient: f64,
};

pub const MineralReaction = struct {
    substrate_limit_mol_per_m3_step: f64,
    equilibrium_phosphate_activity_mol_per_m3: f64,
    precipitation_mol_per_m3_step: f64,
};

pub const Result = struct {
    aluminum_phosphate: MineralReaction,
    iron_phosphate: MineralReaction,
    dicalcium_phosphate: MineralReaction,
    hydroxyapatite: MineralReaction,
    monocalcium_phosphate: MineralReaction,
};

/// Direct source-order translation of SOLUTE.F lines 4569--4607.
///
/// Positive rates precipitate solid; negative rates dissolve it. This pure
/// one-cell kernel does not allocate or mutate authoritative inventories.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const aqueous = inputs.aqueous;
    const minerals = inputs.minerals;
    const products = inputs.solubility_products;
    const dissociation = inputs.dissociation;
    const kinetics = inputs.kinetics;
    const fraction = kinetics.substrate_limit_fraction_per_step;
    const phosphate_maximum =
        kinetics.maximum_phosphate_precipitation_mol_per_m3_step;

    // SOLUTE.F 4571--4574.
    const aluminum_limit = fraction *
        @min(
            aqueous.aluminum_concentration_mol_per_m3,
            aqueous.dihydrogen_phosphate_concentration_mol_p_per_m3,
        );
    const aluminum_equilibrium =
        products.aluminum_phosphate *
        std.math.pow(f64, aqueous.hydrogen_activity_mol_per_m3, 2.0) /
        (dissociation.dihydrogen_phosphate_mol_per_m3 *
            dissociation.hydrogen_phosphate_mol_per_m3 *
            aqueous.aluminum_activity_mol_per_m3);
    const aluminum_rate = boundedRate(
        (aqueous.dihydrogen_phosphate_activity_mol_p_per_m3 -
            aluminum_equilibrium) /
            inputs.monovalent_activity_coefficient,
        minerals.aluminum_phosphate_mol_per_m3,
        phosphate_maximum,
        phosphate_maximum,
        aluminum_limit,
    );

    // SOLUTE.F 4580--4583.
    const iron_limit = fraction *
        @min(
            aqueous.iron_concentration_mol_per_m3,
            aqueous.dihydrogen_phosphate_concentration_mol_p_per_m3,
        );
    const iron_equilibrium =
        products.iron_phosphate *
        std.math.pow(f64, aqueous.hydrogen_activity_mol_per_m3, 2.0) /
        (dissociation.dihydrogen_phosphate_mol_per_m3 *
            dissociation.hydrogen_phosphate_mol_per_m3 *
            aqueous.iron_activity_mol_per_m3);
    const iron_rate = boundedRate(
        (aqueous.dihydrogen_phosphate_activity_mol_p_per_m3 -
            iron_equilibrium) /
            inputs.monovalent_activity_coefficient,
        minerals.iron_phosphate_mol_per_m3,
        phosphate_maximum,
        phosphate_maximum,
        iron_limit,
    );

    // SOLUTE.F 4589--4592. The source limits this HPO4 mineral with H2PO4
    // concentration; preserve that operand instead of substituting HPO4.
    const dicalcium_limit = fraction *
        @min(
            aqueous.calcium_concentration_mol_per_m3,
            aqueous.dihydrogen_phosphate_concentration_mol_p_per_m3,
        );
    const dicalcium_equilibrium =
        products.dicalcium_phosphate /
        aqueous.calcium_activity_mol_per_m3;
    const dicalcium_rate = boundedRate(
        (aqueous.hydrogen_phosphate_activity_mol_p_per_m3 -
            dicalcium_equilibrium) /
            inputs.divalent_activity_coefficient,
        minerals.dicalcium_phosphate_mol_per_m3,
        phosphate_maximum,
        phosphate_maximum,
        dicalcium_limit,
    );

    // SOLUTE.F 4596--4599 retains source exponent 0.333.
    const apatite_limit = fraction *
        @min(
            aqueous.calcium_concentration_mol_per_m3 / 5.0,
            aqueous.dihydrogen_phosphate_concentration_mol_p_per_m3 / 3.0,
        );
    const apatite_equilibrium = std.math.pow(
        f64,
        products.hydroxyapatite *
            std.math.pow(f64, aqueous.hydrogen_activity_mol_per_m3, 6.0) /
            (std.math.pow(f64, aqueous.calcium_activity_mol_per_m3, 5.0) *
                aqueous.hydroxide_activity_mol_per_m3 *
                std.math.pow(
                    f64,
                    dissociation.dihydrogen_phosphate_mol_per_m3 *
                        dissociation.hydrogen_phosphate_mol_per_m3,
                    3.0,
                )),
        0.333,
    );
    const apatite_maximum =
        kinetics.maximum_apatite_precipitation_mol_per_m3_step;
    const apatite_rate = boundedRate(
        (aqueous.dihydrogen_phosphate_activity_mol_p_per_m3 -
            apatite_equilibrium) /
            inputs.monovalent_activity_coefficient,
        minerals.hydroxyapatite_mol_per_m3,
        apatite_maximum,
        apatite_maximum,
        apatite_limit,
    );

    // SOLUTE.F 4604--4607. Dissolution uses TPZ while precipitation retains
    // TPD, including the source's asymmetric kinetic ceilings.
    const monocalcium_limit = fraction *
        @min(
            aqueous.calcium_concentration_mol_per_m3,
            aqueous.dihydrogen_phosphate_concentration_mol_p_per_m3,
        );
    const monocalcium_equilibrium = std.math.pow(
        f64,
        products.monocalcium_phosphate /
            aqueous.calcium_activity_mol_per_m3,
        0.5,
    );
    const monocalcium_rate = boundedRate(
        (aqueous.dihydrogen_phosphate_activity_mol_p_per_m3 -
            monocalcium_equilibrium) /
            inputs.monovalent_activity_coefficient,
        minerals.monocalcium_phosphate_mol_per_m3,
        kinetics.maximum_monocalcium_dissolution_mol_per_m3_step,
        phosphate_maximum,
        monocalcium_limit,
    );

    const result: Result = .{
        .aluminum_phosphate = reaction(
            aluminum_limit,
            aluminum_equilibrium,
            aluminum_rate,
        ),
        .iron_phosphate = reaction(
            iron_limit,
            iron_equilibrium,
            iron_rate,
        ),
        .dicalcium_phosphate = reaction(
            dicalcium_limit,
            dicalcium_equilibrium,
            dicalcium_rate,
        ),
        .hydroxyapatite = reaction(
            apatite_limit,
            apatite_equilibrium,
            apatite_rate,
        ),
        .monocalcium_phosphate = reaction(
            monocalcium_limit,
            monocalcium_equilibrium,
            monocalcium_rate,
        ),
    };
    try validateResult(result);
    return result;
}

fn boundedRate(
    driving: f64,
    solid_mol_per_m3: f64,
    maximum_dissolution: f64,
    maximum_precipitation: f64,
    substrate_limit: f64,
) f64 {
    return @max(
        -@max(0.0, solid_mol_per_m3),
        -maximum_dissolution,
        @min(maximum_precipitation, substrate_limit, driving),
    );
}

fn reaction(
    substrate_limit: f64,
    equilibrium_activity: f64,
    rate: f64,
) MineralReaction {
    return .{
        .substrate_limit_mol_per_m3_step = substrate_limit,
        .equilibrium_phosphate_activity_mol_per_m3 = equilibrium_activity,
        .precipitation_mol_per_m3_step = rate,
    };
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(AqueousState).@"struct".fields) |field| {
        const value = @field(inputs.aqueous, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterPhosphateMineralInput;
    }
    inline for (@typeInfo(MineralState).@"struct".fields) |field| {
        const value = @field(inputs.minerals, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterPhosphateMineralInput;
    }
    inline for (@typeInfo(SolubilityProducts).@"struct".fields) |field| {
        const value = @field(inputs.solubility_products, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidSurfaceLitterPhosphateMineralParameter;
    }
    inline for (@typeInfo(DissociationConstants).@"struct".fields) |field| {
        const value = @field(inputs.dissociation, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidSurfaceLitterPhosphateMineralParameter;
    }
    inline for (@typeInfo(Kinetics).@"struct".fields) |field| {
        const value = @field(inputs.kinetics, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterPhosphateMineralParameter;
    }
    if (inputs.aqueous.hydrogen_activity_mol_per_m3 <= 0 or
        inputs.aqueous.hydroxide_activity_mol_per_m3 <= 0 or
        inputs.aqueous.aluminum_activity_mol_per_m3 <= 0 or
        inputs.aqueous.iron_activity_mol_per_m3 <= 0 or
        inputs.aqueous.calcium_activity_mol_per_m3 <= 0 or
        !std.math.isFinite(inputs.monovalent_activity_coefficient) or
        inputs.monovalent_activity_coefficient <= 0 or
        !std.math.isFinite(inputs.divalent_activity_coefficient) or
        inputs.divalent_activity_coefficient <= 0 or
        inputs.kinetics.substrate_limit_fraction_per_step > 1)
    {
        return error.InvalidSurfaceLitterPhosphateMineralInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        const mineral = @field(result, field.name);
        inline for (@typeInfo(MineralReaction).@"struct".fields) |member| {
            if (!std.math.isFinite(@field(mineral, member.name)))
                return error.NonFiniteSurfaceLitterPhosphateMineralResult;
        }
    }
}

fn testInputs() Inputs {
    return .{
        .aqueous = .{
            .hydrogen_activity_mol_per_m3 = 0.4,
            .hydroxide_activity_mol_per_m3 = 0.25,
            .aluminum_concentration_mol_per_m3 = 0.7,
            .aluminum_activity_mol_per_m3 = 0.35,
            .iron_concentration_mol_per_m3 = 0.6,
            .iron_activity_mol_per_m3 = 0.3,
            .calcium_concentration_mol_per_m3 = 1.2,
            .calcium_activity_mol_per_m3 = 0.72,
            .hydrogen_phosphate_concentration_mol_p_per_m3 = 0.8,
            .hydrogen_phosphate_activity_mol_p_per_m3 = 0.48,
            .dihydrogen_phosphate_concentration_mol_p_per_m3 = 0.9,
            .dihydrogen_phosphate_activity_mol_p_per_m3 = 0.72,
        },
        .minerals = .{
            .aluminum_phosphate_mol_per_m3 = 0.2,
            .iron_phosphate_mol_per_m3 = 0.25,
            .dicalcium_phosphate_mol_per_m3 = 0.3,
            .hydroxyapatite_mol_per_m3 = 0.35,
            .monocalcium_phosphate_mol_per_m3 = 0.4,
        },
        .solubility_products = .{
            .aluminum_phosphate = 0.02,
            .iron_phosphate = 0.03,
            .dicalcium_phosphate = 0.04,
            .hydroxyapatite = 0.05,
            .monocalcium_phosphate = 0.06,
        },
        .dissociation = .{
            .dihydrogen_phosphate_mol_per_m3 = 0.2,
            .hydrogen_phosphate_mol_per_m3 = 0.15,
        },
        .kinetics = .{
            .substrate_limit_fraction_per_step = 0.4,
            .maximum_phosphate_precipitation_mol_per_m3_step = 0.1,
            .maximum_apatite_precipitation_mol_per_m3_step = 0.08,
            .maximum_monocalcium_dissolution_mol_per_m3_step = 0.06,
        },
        .monovalent_activity_coefficient = 0.8,
        .divalent_activity_coefficient = 0.6,
    };
}

test "SOLUTE surface phosphate minerals preserve every source expression" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const aqueous = inputs.aqueous;
    const products = inputs.solubility_products;
    const dissociation = inputs.dissociation;
    const fraction = inputs.kinetics.substrate_limit_fraction_per_step;
    const phosphate_maximum =
        inputs.kinetics.maximum_phosphate_precipitation_mol_per_m3_step;

    const aluminum_limit = fraction *
        @min(
            aqueous.aluminum_concentration_mol_per_m3,
            aqueous.dihydrogen_phosphate_concentration_mol_p_per_m3,
        );
    const aluminum_equilibrium =
        products.aluminum_phosphate *
        std.math.pow(f64, aqueous.hydrogen_activity_mol_per_m3, 2.0) /
        (dissociation.dihydrogen_phosphate_mol_per_m3 *
            dissociation.hydrogen_phosphate_mol_per_m3 *
            aqueous.aluminum_activity_mol_per_m3);
    const aluminum_rate = @max(
        -@max(0.0, inputs.minerals.aluminum_phosphate_mol_per_m3),
        -phosphate_maximum,
        @min(
            phosphate_maximum,
            aluminum_limit,
            (aqueous.dihydrogen_phosphate_activity_mol_p_per_m3 -
                aluminum_equilibrium) /
                inputs.monovalent_activity_coefficient,
        ),
    );
    try std.testing.expectEqual(
        aluminum_limit,
        result.aluminum_phosphate.substrate_limit_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        aluminum_equilibrium,
        result.aluminum_phosphate.equilibrium_phosphate_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        aluminum_rate,
        result.aluminum_phosphate.precipitation_mol_per_m3_step,
    );

    const iron_equilibrium =
        products.iron_phosphate *
        std.math.pow(f64, aqueous.hydrogen_activity_mol_per_m3, 2.0) /
        (dissociation.dihydrogen_phosphate_mol_per_m3 *
            dissociation.hydrogen_phosphate_mol_per_m3 *
            aqueous.iron_activity_mol_per_m3);
    try std.testing.expectEqual(
        iron_equilibrium,
        result.iron_phosphate.equilibrium_phosphate_activity_mol_per_m3,
    );

    const dicalcium_equilibrium =
        products.dicalcium_phosphate /
        aqueous.calcium_activity_mol_per_m3;
    try std.testing.expectEqual(
        dicalcium_equilibrium,
        result.dicalcium_phosphate.equilibrium_phosphate_activity_mol_per_m3,
    );

    const apatite_equilibrium = std.math.pow(
        f64,
        products.hydroxyapatite *
            std.math.pow(f64, aqueous.hydrogen_activity_mol_per_m3, 6.0) /
            (std.math.pow(f64, aqueous.calcium_activity_mol_per_m3, 5.0) *
                aqueous.hydroxide_activity_mol_per_m3 *
                std.math.pow(
                    f64,
                    dissociation.dihydrogen_phosphate_mol_per_m3 *
                        dissociation.hydrogen_phosphate_mol_per_m3,
                    3.0,
                )),
        0.333,
    );
    try std.testing.expectEqual(
        apatite_equilibrium,
        result.hydroxyapatite.equilibrium_phosphate_activity_mol_per_m3,
    );

    const monocalcium_equilibrium = std.math.pow(
        f64,
        products.monocalcium_phosphate /
            aqueous.calcium_activity_mol_per_m3,
        0.5,
    );
    try std.testing.expectEqual(
        monocalcium_equilibrium,
        result.monocalcium_phosphate.equilibrium_phosphate_activity_mol_per_m3,
    );
}

test "surface mineral rates preserve source substrate and kinetic ceilings" {
    var inputs = testInputs();
    inputs.aqueous.dihydrogen_phosphate_activity_mol_p_per_m3 = 100;
    var result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        inputs.kinetics.maximum_phosphate_precipitation_mol_per_m3_step,
        result.aluminum_phosphate.precipitation_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        inputs.kinetics.maximum_apatite_precipitation_mol_per_m3_step,
        result.hydroxyapatite.precipitation_mol_per_m3_step,
    );

    inputs.aqueous.dihydrogen_phosphate_activity_mol_p_per_m3 = 0;
    inputs.solubility_products.monocalcium_phosphate = 100;
    inputs.minerals.monocalcium_phosphate_mol_per_m3 = 1;
    result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        -inputs.kinetics.maximum_monocalcium_dissolution_mol_per_m3_step,
        result.monocalcium_phosphate.precipitation_mol_per_m3_step,
    );
}

test "surface dicalcium source limit uses dihydrogen phosphate" {
    var inputs = testInputs();
    inputs.aqueous.hydrogen_phosphate_concentration_mol_p_per_m3 = 10;
    inputs.aqueous.hydrogen_phosphate_activity_mol_p_per_m3 = 10;
    inputs.aqueous.dihydrogen_phosphate_concentration_mol_p_per_m3 = 0.01;
    inputs.kinetics.substrate_limit_fraction_per_step = 0.5;
    inputs.kinetics.maximum_phosphate_precipitation_mol_per_m3_step = 1;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectEqual(
        @as(f64, 0.005),
        result.dicalcium_phosphate.substrate_limit_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.005),
        result.dicalcium_phosphate.precipitation_mol_per_m3_step,
    );
}

test "surface apatite target retains source 0.333 exponent" {
    var inputs = testInputs();
    inputs.solubility_products.hydroxyapatite =
        std.math.pow(f64, inputs.aqueous.calcium_activity_mol_per_m3, 5.0) *
        inputs.aqueous.hydroxide_activity_mol_per_m3 *
        std.math.pow(
            f64,
            inputs.dissociation.dihydrogen_phosphate_mol_per_m3 *
                inputs.dissociation.hydrogen_phosphate_mol_per_m3,
            3.0,
        ) *
        8.0 /
        std.math.pow(f64, inputs.aqueous.hydrogen_activity_mol_per_m3, 6.0);
    const result = try calculateSourceOrder(inputs);
    const source = std.math.pow(f64, 8.0, 0.333);
    const exact_cube_root = std.math.cbrt(@as(f64, 8.0));

    try std.testing.expectEqual(
        source,
        result.hydroxyapatite.equilibrium_phosphate_activity_mol_per_m3,
    );
    try std.testing.expect(source != exact_cube_root);
}

test "surface phosphate minerals reject invalid input and overflow" {
    var inputs = testInputs();
    inputs.minerals.hydroxyapatite_mol_per_m3 = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterPhosphateMineralInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.dissociation.hydrogen_phosphate_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterPhosphateMineralParameter,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.aqueous.hydrogen_activity_mol_per_m3 =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterPhosphateMineralResult,
        calculateSourceOrder(inputs),
    );
}
