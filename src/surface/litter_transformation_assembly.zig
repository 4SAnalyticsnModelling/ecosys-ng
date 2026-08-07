const std = @import("std");

pub const PhosphateMineralRates = struct {
    aluminum_phosphate_mol_mineral_per_m3_step: f64,
    iron_phosphate_mol_mineral_per_m3_step: f64,
    dicalcium_phosphate_mol_mineral_per_m3_step: f64,
    monocalcium_phosphate_mol_mineral_per_m3_step: f64,
    hydroxyapatite_mol_mineral_per_m3_step: f64,
};

pub const CationExchangeRates = struct {
    hydrogen_mol_per_megagram_step: f64,
    aluminum_mol_per_megagram_step: f64,
    iron_mol_per_megagram_step: f64,
    calcium_mol_per_megagram_step: f64,
    magnesium_mol_per_megagram_step: f64,
    sodium_mol_per_megagram_step: f64,
    potassium_mol_per_megagram_step: f64,
};

pub const Inputs = struct {
    ammonium_association_mol_n_per_m3_step: f64,
    ammonium_exchange_mol_n_per_megagram_step: f64,
    dihydrogen_phosphate_association_mol_p_per_m3_step: f64,
    phosphate_minerals: PhosphateMineralRates,
    cation_exchange: CationExchangeRates,
    litter_mass_per_water_volume_megagrams_per_m3: f64,
};

pub const NitrogenTransformations = struct {
    ammonium_mol_n_per_m3_step: f64,
    ammonia_mol_n_per_m3_step: f64,
};

pub const PhosphorusTransformations = struct {
    hydrogen_phosphate_mol_p_per_m3_step: f64,
    dihydrogen_phosphate_mol_p_per_m3_step: f64,
};

pub const IonTransformations = struct {
    hydrogen_mol_per_m3_step: f64,
    hydroxide_mol_per_m3_step: f64,
    aluminum_mol_per_m3_step: f64,
    iron_mol_per_m3_step: f64,
    calcium_mol_per_m3_step: f64,
    magnesium_mol_per_m3_step: f64,
    sodium_mol_per_m3_step: f64,
    potassium_mol_per_m3_step: f64,
};

pub const Result = struct {
    nitrogen: NitrogenTransformations,
    phosphorus: PhosphorusTransformations,
    ions: IonTransformations,
};

/// Direct source-order translation of SOLUTE.F lines 4987--4998.
///
/// This pure one-cell assembly converts exchange rates from mol Mg^-1 per
/// step to mol m^-3 per step using runtime litter mass-to-water density.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const density = inputs.litter_mass_per_water_volume_megagrams_per_m3;
    const minerals = inputs.phosphate_minerals;
    const exchange = inputs.cation_exchange;

    // SOLUTE.F 4987--4998, preserved in assignment order.
    const ammonium =
        inputs.ammonium_association_mol_n_per_m3_step -
        inputs.ammonium_exchange_mol_n_per_megagram_step * density;
    const ammonia = -inputs.ammonium_association_mol_n_per_m3_step;
    const hydrogen_phosphate =
        -inputs.dihydrogen_phosphate_association_mol_p_per_m3_step;
    const dihydrogen_phosphate =
        inputs.dihydrogen_phosphate_association_mol_p_per_m3_step -
        minerals.aluminum_phosphate_mol_mineral_per_m3_step -
        minerals.iron_phosphate_mol_mineral_per_m3_step -
        minerals.dicalcium_phosphate_mol_mineral_per_m3_step -
        2.0 * minerals.monocalcium_phosphate_mol_mineral_per_m3_step -
        3.0 * minerals.hydroxyapatite_mol_mineral_per_m3_step;
    const hydrogen = -exchange.hydrogen_mol_per_megagram_step * density;
    const hydroxide = 0.0;
    const aluminum = -exchange.aluminum_mol_per_megagram_step * density;
    const iron = -exchange.iron_mol_per_megagram_step * density;
    const calcium = -exchange.calcium_mol_per_megagram_step * density;
    const magnesium = -exchange.magnesium_mol_per_megagram_step * density;
    const sodium = -exchange.sodium_mol_per_megagram_step * density;
    const potassium = -exchange.potassium_mol_per_megagram_step * density;

    const result: Result = .{
        .nitrogen = .{
            .ammonium_mol_n_per_m3_step = ammonium,
            .ammonia_mol_n_per_m3_step = ammonia,
        },
        .phosphorus = .{
            .hydrogen_phosphate_mol_p_per_m3_step = hydrogen_phosphate,
            .dihydrogen_phosphate_mol_p_per_m3_step = dihydrogen_phosphate,
        },
        .ions = .{
            .hydrogen_mol_per_m3_step = hydrogen,
            .hydroxide_mol_per_m3_step = hydroxide,
            .aluminum_mol_per_m3_step = aluminum,
            .iron_mol_per_m3_step = iron,
            .calcium_mol_per_m3_step = calcium,
            .magnesium_mol_per_m3_step = magnesium,
            .sodium_mol_per_m3_step = sodium,
            .potassium_mol_per_m3_step = potassium,
        },
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (.{
        inputs.ammonium_association_mol_n_per_m3_step,
        inputs.ammonium_exchange_mol_n_per_megagram_step,
        inputs.dihydrogen_phosphate_association_mol_p_per_m3_step,
        inputs.litter_mass_per_water_volume_megagrams_per_m3,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.InvalidSurfaceLitterTransformationAssemblyInput;
    }
    if (inputs.litter_mass_per_water_volume_megagrams_per_m3 <= 0)
        return error.InvalidSurfaceLitterTransformationAssemblyInput;
    try validateFiniteStruct(inputs.phosphate_minerals);
    try validateFiniteStruct(inputs.cation_exchange);
}

fn validateFiniteStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(value, field.name)))
            return error.InvalidSurfaceLitterTransformationAssemblyInput;
    }
}

fn validateResult(result: Result) !void {
    try validateFiniteResultStruct(result.nitrogen);
    try validateFiniteResultStruct(result.phosphorus);
    try validateFiniteResultStruct(result.ions);
}

fn validateFiniteResultStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(value, field.name)))
            return error.NonFiniteSurfaceLitterTransformationAssemblyResult;
    }
}

fn testInputs() Inputs {
    return .{
        .ammonium_association_mol_n_per_m3_step = 0.4,
        .ammonium_exchange_mol_n_per_megagram_step = 0.1,
        .dihydrogen_phosphate_association_mol_p_per_m3_step = 0.3,
        .phosphate_minerals = .{
            .aluminum_phosphate_mol_mineral_per_m3_step = 0.01,
            .iron_phosphate_mol_mineral_per_m3_step = 0.02,
            .dicalcium_phosphate_mol_mineral_per_m3_step = 0.03,
            .monocalcium_phosphate_mol_mineral_per_m3_step = 0.04,
            .hydroxyapatite_mol_mineral_per_m3_step = 0.05,
        },
        .cation_exchange = .{
            .hydrogen_mol_per_megagram_step = 0.06,
            .aluminum_mol_per_megagram_step = 0.07,
            .iron_mol_per_megagram_step = 0.08,
            .calcium_mol_per_megagram_step = 0.09,
            .magnesium_mol_per_megagram_step = 0.10,
            .sodium_mol_per_megagram_step = 0.11,
            .potassium_mol_per_megagram_step = 0.12,
        },
        .litter_mass_per_water_volume_megagrams_per_m3 = 2,
    };
}

test "SOLUTE surface transformation assembly preserves every source equation" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const density = inputs.litter_mass_per_water_volume_megagrams_per_m3;
    const minerals = inputs.phosphate_minerals;
    const exchange = inputs.cation_exchange;

    try std.testing.expectEqual(
        inputs.ammonium_association_mol_n_per_m3_step -
            inputs.ammonium_exchange_mol_n_per_megagram_step * density,
        result.nitrogen.ammonium_mol_n_per_m3_step,
    );
    try std.testing.expectEqual(
        -inputs.ammonium_association_mol_n_per_m3_step,
        result.nitrogen.ammonia_mol_n_per_m3_step,
    );
    try std.testing.expectEqual(
        -inputs.dihydrogen_phosphate_association_mol_p_per_m3_step,
        result.phosphorus.hydrogen_phosphate_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        inputs.dihydrogen_phosphate_association_mol_p_per_m3_step -
            minerals.aluminum_phosphate_mol_mineral_per_m3_step -
            minerals.iron_phosphate_mol_mineral_per_m3_step -
            minerals.dicalcium_phosphate_mol_mineral_per_m3_step -
            2.0 * minerals.monocalcium_phosphate_mol_mineral_per_m3_step -
            3.0 * minerals.hydroxyapatite_mol_mineral_per_m3_step,
        result.phosphorus.dihydrogen_phosphate_mol_p_per_m3_step,
    );
    try std.testing.expectEqual(
        -exchange.calcium_mol_per_megagram_step * density,
        result.ions.calcium_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.ions.hydroxide_mol_per_m3_step,
    );
}

test "surface transformation assembly closes nitrogen and phosphorus" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const density = inputs.litter_mass_per_water_volume_megagrams_per_m3;
    const nitrogen_with_exchange =
        result.nitrogen.ammonium_mol_n_per_m3_step +
        result.nitrogen.ammonia_mol_n_per_m3_step +
        inputs.ammonium_exchange_mol_n_per_megagram_step * density;
    const minerals = inputs.phosphate_minerals;
    const phosphorus_with_solids =
        result.phosphorus.hydrogen_phosphate_mol_p_per_m3_step +
        result.phosphorus.dihydrogen_phosphate_mol_p_per_m3_step +
        minerals.aluminum_phosphate_mol_mineral_per_m3_step +
        minerals.iron_phosphate_mol_mineral_per_m3_step +
        minerals.dicalcium_phosphate_mol_mineral_per_m3_step +
        2.0 * minerals.monocalcium_phosphate_mol_mineral_per_m3_step +
        3.0 * minerals.hydroxyapatite_mol_mineral_per_m3_step;

    try std.testing.expectApproxEqAbs(
        @as(f64, 0),
        nitrogen_with_exchange,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0),
        phosphorus_with_solids,
        1.0e-15,
    );
}

test "surface transformation assembly preserves signed reversibility" {
    var inputs = testInputs();
    inputs.ammonium_association_mol_n_per_m3_step *= -1;
    inputs.dihydrogen_phosphate_association_mol_p_per_m3_step *= -1;
    inline for (@typeInfo(PhosphateMineralRates).@"struct".fields) |field|
        @field(inputs.phosphate_minerals, field.name) *= -1;
    inline for (@typeInfo(CationExchangeRates).@"struct".fields) |field|
        @field(inputs.cation_exchange, field.name) *= -1;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expect(result.nitrogen.ammonia_mol_n_per_m3_step > 0);
    try std.testing.expect(result.ions.calcium_mol_per_m3_step > 0);
}

test "surface transformation assembly rejects invalid input and overflow" {
    var inputs = testInputs();
    inputs.litter_mass_per_water_volume_megagrams_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterTransformationAssemblyInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.phosphate_minerals.dicalcium_phosphate_mol_mineral_per_m3_step =
        std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterTransformationAssemblyInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.litter_mass_per_water_volume_megagrams_per_m3 =
        std.math.floatMax(f64);
    inputs.cation_exchange.calcium_mol_per_megagram_step =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterTransformationAssemblyResult,
        calculateSourceOrder(inputs),
    );
}
