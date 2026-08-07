const std = @import("std");

pub const AqueousState = struct {
    hydrogen_activity_mol_per_m3: f64,
    hydroxide_activity_mol_per_m3: f64,
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

/// Runtime equilibrium products retain the source parameter dimensions.
/// Each expression below reduces its product to an activity in mol m^-3.
pub const EquilibriumProducts = struct {
    aluminum_hydroxide: f64,
    aluminum_phosphate_hydrogen_form: f64,
    iron_hydroxide: f64,
    iron_phosphate_hydrogen_form: f64,
    dicalcium_phosphate: f64,
    hydroxyapatite_hydrogen_form: f64,
    monocalcium_phosphate: f64,
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
    equilibrium_products: EquilibriumProducts,
    kinetics: Kinetics,
    monovalent_activity_coefficient: f64,
    divalent_activity_coefficient: f64,
};

pub const MetalPhosphateReaction = struct {
    substrate_limit_mol_p_per_m3_step: f64,
    equilibrium_metal_activity_mol_per_m3: f64,
    equilibrium_dihydrogen_phosphate_activity_mol_p_per_m3: f64,
    precipitation_mol_p_per_m3_step: f64,
};

pub const MineralReaction = struct {
    substrate_limit_mol_p_per_m3_step: f64,
    equilibrium_phosphate_activity_mol_p_per_m3: f64,
    precipitation_mol_p_per_m3_step: f64,
};

pub const Result = struct {
    aluminum_phosphate: MetalPhosphateReaction,
    iron_phosphate: MetalPhosphateReaction,
    dicalcium_phosphate: MineralReaction,
    hydroxyapatite: MineralReaction,
    monocalcium_phosphate: MineralReaction,
};

/// Direct source-order translation of SOLUTE.F lines 4684--4719.
///
/// Positive rates precipitate and negative rates dissolve. The caller resolves
/// surface-litter cell index `(0, NY, NX)` and applies returned rates to state.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const aqueous = inputs.aqueous;
    const minerals = inputs.minerals;
    const products = inputs.equilibrium_products;
    const kinetics = inputs.kinetics;
    const fraction = kinetics.substrate_limit_fraction_per_step;
    const phosphate_maximum =
        kinetics.maximum_phosphate_precipitation_mol_per_m3_step;

    // SOLUTE.F 4684--4688. The source does not include dissolved aluminum in
    // XMIN; the equilibrium Al activity is instead derived from OH activity.
    const aluminum_limit =
        fraction * aqueous.dihydrogen_phosphate_concentration_mol_p_per_m3;
    const aluminum_equilibrium = products.aluminum_hydroxide /
        std.math.pow(f64, aqueous.hydroxide_activity_mol_per_m3, 3.0);
    const aluminum_phosphate_equilibrium =
        products.aluminum_phosphate_hydrogen_form *
        std.math.pow(f64, aqueous.hydrogen_activity_mol_per_m3, 2.0) /
        aluminum_equilibrium;
    const aluminum_rate = boundedRate(
        (aqueous.dihydrogen_phosphate_activity_mol_p_per_m3 -
            aluminum_phosphate_equilibrium) /
            inputs.monovalent_activity_coefficient,
        minerals.aluminum_phosphate_mol_per_m3,
        phosphate_maximum,
        phosphate_maximum,
        aluminum_limit,
    );

    // SOLUTE.F 4692--4696 likewise derives equilibrium Fe from OH activity.
    const iron_limit =
        fraction * aqueous.dihydrogen_phosphate_concentration_mol_p_per_m3;
    const iron_equilibrium = products.iron_hydroxide /
        std.math.pow(f64, aqueous.hydroxide_activity_mol_per_m3, 3.0);
    const iron_phosphate_equilibrium =
        products.iron_phosphate_hydrogen_form *
        std.math.pow(f64, aqueous.hydrogen_activity_mol_per_m3, 2.0) /
        iron_equilibrium;
    const iron_rate = boundedRate(
        (aqueous.dihydrogen_phosphate_activity_mol_p_per_m3 -
            iron_phosphate_equilibrium) /
            inputs.monovalent_activity_coefficient,
        minerals.iron_phosphate_mol_per_m3,
        phosphate_maximum,
        phosphate_maximum,
        iron_limit,
    );

    // SOLUTE.F 4700--4705.
    const dicalcium_limit = fraction *
        @min(
            aqueous.calcium_concentration_mol_per_m3,
            aqueous.hydrogen_phosphate_concentration_mol_p_per_m3,
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

    // SOLUTE.F 4709--4712 retains both the unscaled Ca/H2PO4 substrate
    // minimum and the literal source exponent 0.333.
    const apatite_limit = fraction *
        @min(
            aqueous.calcium_concentration_mol_per_m3,
            aqueous.dihydrogen_phosphate_concentration_mol_p_per_m3,
        );
    const apatite_equilibrium = std.math.pow(
        f64,
        products.hydroxyapatite_hydrogen_form *
            std.math.pow(f64, aqueous.hydrogen_activity_mol_per_m3, 7.0) /
            std.math.pow(f64, aqueous.calcium_activity_mol_per_m3, 5.0),
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

    // SOLUTE.F 4716--4719. Dissolution uses TPZ while precipitation uses TPD.
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
        .aluminum_phosphate = .{
            .substrate_limit_mol_p_per_m3_step = aluminum_limit,
            .equilibrium_metal_activity_mol_per_m3 = aluminum_equilibrium,
            .equilibrium_dihydrogen_phosphate_activity_mol_p_per_m3 = aluminum_phosphate_equilibrium,
            .precipitation_mol_p_per_m3_step = aluminum_rate,
        },
        .iron_phosphate = .{
            .substrate_limit_mol_p_per_m3_step = iron_limit,
            .equilibrium_metal_activity_mol_per_m3 = iron_equilibrium,
            .equilibrium_dihydrogen_phosphate_activity_mol_p_per_m3 = iron_phosphate_equilibrium,
            .precipitation_mol_p_per_m3_step = iron_rate,
        },
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
    driving_mol_per_m3_step: f64,
    solid_mol_per_m3: f64,
    maximum_dissolution_mol_per_m3_step: f64,
    maximum_precipitation_mol_per_m3_step: f64,
    substrate_limit_mol_per_m3_step: f64,
) f64 {
    return @max(
        -@max(0.0, solid_mol_per_m3),
        -maximum_dissolution_mol_per_m3_step,
        @min(
            maximum_precipitation_mol_per_m3_step,
            substrate_limit_mol_per_m3_step,
            driving_mol_per_m3_step,
        ),
    );
}

fn reaction(
    substrate_limit_mol_p_per_m3_step: f64,
    equilibrium_phosphate_activity_mol_p_per_m3: f64,
    precipitation_mol_p_per_m3_step: f64,
) MineralReaction {
    return .{
        .substrate_limit_mol_p_per_m3_step = substrate_limit_mol_p_per_m3_step,
        .equilibrium_phosphate_activity_mol_p_per_m3 = equilibrium_phosphate_activity_mol_p_per_m3,
        .precipitation_mol_p_per_m3_step = precipitation_mol_p_per_m3_step,
    };
}

fn validateInputs(inputs: Inputs) !void {
    try validateNonnegativeStruct(inputs.aqueous);
    try validateNonnegativeStruct(inputs.minerals);
    try validatePositiveStruct(inputs.equilibrium_products);
    try validateNonnegativeStruct(inputs.kinetics);
    if (inputs.aqueous.hydrogen_activity_mol_per_m3 <= 0 or
        inputs.aqueous.hydroxide_activity_mol_per_m3 <= 0 or
        inputs.aqueous.calcium_activity_mol_per_m3 <= 0 or
        !std.math.isFinite(inputs.monovalent_activity_coefficient) or
        inputs.monovalent_activity_coefficient <= 0 or
        !std.math.isFinite(inputs.divalent_activity_coefficient) or
        inputs.divalent_activity_coefficient <= 0 or
        inputs.kinetics.substrate_limit_fraction_per_step > 1)
    {
        return error.InvalidSurfaceLitterStaticPhosphateMineralInput;
    }
}

fn validateNonnegativeStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        const field_value = @field(value, field.name);
        if (!std.math.isFinite(field_value) or field_value < 0)
            return error.InvalidSurfaceLitterStaticPhosphateMineralInput;
    }
}

fn validatePositiveStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        const field_value = @field(value, field.name);
        if (!std.math.isFinite(field_value) or field_value <= 0)
            return error.InvalidSurfaceLitterStaticPhosphateMineralParameter;
    }
}

fn validateResult(result: Result) !void {
    try validateFiniteStruct(result.aluminum_phosphate);
    try validateFiniteStruct(result.iron_phosphate);
    try validateFiniteStruct(result.dicalcium_phosphate);
    try validateFiniteStruct(result.hydroxyapatite);
    try validateFiniteStruct(result.monocalcium_phosphate);
}

fn validateFiniteStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(value, field.name)))
            return error.NonFiniteSurfaceLitterStaticPhosphateMineralResult;
    }
}

fn testInputs() Inputs {
    return .{
        .aqueous = .{
            .hydrogen_activity_mol_per_m3 = 0.4,
            .hydroxide_activity_mol_per_m3 = 0.5,
            .calcium_concentration_mol_per_m3 = 1.2,
            .calcium_activity_mol_per_m3 = 0.6,
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
        .equilibrium_products = .{
            .aluminum_hydroxide = 0.001,
            .aluminum_phosphate_hydrogen_form = 0.02,
            .iron_hydroxide = 0.002,
            .iron_phosphate_hydrogen_form = 0.03,
            .dicalcium_phosphate = 0.12,
            .hydroxyapatite_hydrogen_form = 0.01,
            .monocalcium_phosphate = 0.06,
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

test "SOLUTE static phosphate minerals preserve every source expression" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const aqueous = inputs.aqueous;
    const products = inputs.equilibrium_products;
    const fraction = inputs.kinetics.substrate_limit_fraction_per_step;

    const expected_aluminum_activity = products.aluminum_hydroxide /
        std.math.pow(f64, aqueous.hydroxide_activity_mol_per_m3, 3.0);
    const expected_aluminum_phosphate =
        products.aluminum_phosphate_hydrogen_form *
        std.math.pow(f64, aqueous.hydrogen_activity_mol_per_m3, 2.0) /
        expected_aluminum_activity;
    try std.testing.expectEqual(
        fraction * aqueous.dihydrogen_phosphate_concentration_mol_p_per_m3,
        result.aluminum_phosphate.substrate_limit_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        expected_aluminum_activity,
        result.aluminum_phosphate.equilibrium_metal_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        expected_aluminum_phosphate,
        result.aluminum_phosphate
            .equilibrium_dihydrogen_phosphate_activity_mol_p_per_m3,
    );

    const expected_iron_activity = products.iron_hydroxide /
        std.math.pow(f64, aqueous.hydroxide_activity_mol_per_m3, 3.0);
    const expected_iron_phosphate =
        products.iron_phosphate_hydrogen_form *
        std.math.pow(f64, aqueous.hydrogen_activity_mol_per_m3, 2.0) /
        expected_iron_activity;
    try std.testing.expectEqual(
        expected_iron_activity,
        result.iron_phosphate.equilibrium_metal_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        expected_iron_phosphate,
        result.iron_phosphate
            .equilibrium_dihydrogen_phosphate_activity_mol_p_per_m3,
    );

    try std.testing.expectEqual(
        fraction *
            @min(
                aqueous.calcium_concentration_mol_per_m3,
                aqueous.hydrogen_phosphate_concentration_mol_p_per_m3,
            ),
        result.dicalcium_phosphate.substrate_limit_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        products.dicalcium_phosphate /
            aqueous.calcium_activity_mol_per_m3,
        result.dicalcium_phosphate
            .equilibrium_phosphate_activity_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        fraction *
            @min(
                aqueous.calcium_concentration_mol_per_m3,
                aqueous.dihydrogen_phosphate_concentration_mol_p_per_m3,
            ),
        result.hydroxyapatite.substrate_limit_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        std.math.pow(
            f64,
            products.monocalcium_phosphate /
                aqueous.calcium_activity_mol_per_m3,
            0.5,
        ),
        result.monocalcium_phosphate
            .equilibrium_phosphate_activity_mol_p_per_m3,
    );
}

test "static mineral rates preserve source kinetic and substrate ceilings" {
    const result = try calculateSourceOrder(testInputs());

    try std.testing.expectEqual(
        0.1,
        result.aluminum_phosphate.precipitation_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        0.1,
        result.iron_phosphate.precipitation_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        0.1,
        result.dicalcium_phosphate.precipitation_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        0.08,
        result.hydroxyapatite.precipitation_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        0.1,
        result.monocalcium_phosphate.precipitation_mol_p_per_m3_step,
    );
}

test "static mineral dissolution respects solid and asymmetric limits" {
    var inputs = testInputs();
    inputs.aqueous.dihydrogen_phosphate_activity_mol_p_per_m3 = 0;
    inputs.aqueous.hydrogen_phosphate_activity_mol_p_per_m3 = 0;
    inputs.minerals.aluminum_phosphate_mol_per_m3 = 0.03;
    inputs.minerals.monocalcium_phosphate_mol_per_m3 = 0.4;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectEqual(
        -0.03,
        result.aluminum_phosphate.precipitation_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        -inputs.kinetics.maximum_monocalcium_dissolution_mol_per_m3_step,
        result.monocalcium_phosphate.precipitation_mol_p_per_m3_step,
    );
}

test "static apatite equilibrium retains source exponent 0.333" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const inner = inputs.equilibrium_products
        .hydroxyapatite_hydrogen_form *
        std.math.pow(
            f64,
            inputs.aqueous.hydrogen_activity_mol_per_m3,
            7.0,
        ) /
        std.math.pow(
            f64,
            inputs.aqueous.calcium_activity_mol_per_m3,
            5.0,
        );
    const source_expected = std.math.pow(f64, inner, 0.333);

    try std.testing.expectEqual(
        source_expected,
        result.hydroxyapatite.equilibrium_phosphate_activity_mol_p_per_m3,
    );
    try std.testing.expect(
        source_expected != std.math.pow(f64, inner, 1.0 / 3.0),
    );
}

test "static phosphate minerals reject invalid input and overflow" {
    var inputs = testInputs();
    inputs.aqueous.hydroxide_activity_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticPhosphateMineralInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.kinetics.substrate_limit_fraction_per_step = 1.1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticPhosphateMineralInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.equilibrium_products.dicalcium_phosphate = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticPhosphateMineralParameter,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.aqueous.hydrogen_activity_mol_per_m3 =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterStaticPhosphateMineralResult,
        calculateSourceOrder(inputs),
    );
}
