const std = @import("std");

pub const AqueousNutrientTotals = struct {
    ammonium_mol_n_per_step: f64,
    ammonia_mol_n_per_step: f64,
    nitrate_mol_n_per_step: f64,
    hydrogen_phosphate_mol_p_per_step: f64,
    dihydrogen_phosphate_mol_p_per_step: f64,
};

pub const ExchangeTotals = struct {
    ammonium_mol_n_per_step: f64,
    hydrogen_mol_per_step: f64,
    aluminum_mol_per_step: f64,
    iron_mol_per_step: f64,
    calcium_mol_per_step: f64,
    magnesium_mol_per_step: f64,
    sodium_mol_per_step: f64,
    potassium_mol_per_step: f64,
};

pub const PhosphateMineralTotals = struct {
    aluminum_phosphate_mol_mineral_per_step: f64,
    iron_phosphate_mol_mineral_per_step: f64,
    dicalcium_phosphate_mol_mineral_per_step: f64,
    hydroxyapatite_mol_mineral_per_step: f64,
    monocalcium_phosphate_mol_mineral_per_step: f64,
};

pub const AqueousCationTotals = struct {
    aluminum_mol_per_step: f64,
    iron_mol_per_step: f64,
    calcium_mol_per_step: f64,
    magnesium_mol_per_step: f64,
    sodium_mol_per_step: f64,
    potassium_mol_per_step: f64,
};

pub const FertilizerInventories = struct {
    ammonium_mol_n: f64,
    ammonia_mol_n: f64,
    urea_mol_n: f64,
    nitrate_mol_n: f64,
};

pub const State = struct {
    aqueous_nutrients: AqueousNutrientTotals,
    exchange: ExchangeTotals,
    phosphate_minerals: PhosphateMineralTotals,
    aqueous_cations: AqueousCationTotals,
    fertilizer: FertilizerInventories,
};

pub const AqueousNutrientRates = struct {
    ammonium_mol_n_per_m3_step: f64,
    ammonia_mol_n_per_m3_step: f64,
    hydrogen_phosphate_mol_p_per_m3_step: f64,
    dihydrogen_phosphate_mol_p_per_m3_step: f64,
};

pub const ExchangeRates = struct {
    ammonium_mol_n_per_megagram_step: f64,
    hydrogen_mol_per_megagram_step: f64,
    aluminum_mol_per_megagram_step: f64,
    iron_mol_per_megagram_step: f64,
    calcium_mol_per_megagram_step: f64,
    magnesium_mol_per_megagram_step: f64,
    sodium_mol_per_megagram_step: f64,
    potassium_mol_per_megagram_step: f64,
};

pub const PhosphateMineralRates = struct {
    aluminum_phosphate_mol_mineral_per_m3_step: f64,
    iron_phosphate_mol_mineral_per_m3_step: f64,
    dicalcium_phosphate_mol_mineral_per_m3_step: f64,
    hydroxyapatite_mol_mineral_per_m3_step: f64,
    monocalcium_phosphate_mol_mineral_per_m3_step: f64,
};

pub const AqueousCationRates = struct {
    aluminum_mol_per_m3_step: f64,
    iron_mol_per_m3_step: f64,
    calcium_mol_per_m3_step: f64,
    magnesium_mol_per_m3_step: f64,
    sodium_mol_per_m3_step: f64,
    potassium_mol_per_m3_step: f64,
};

pub const FertilizerDissolution = struct {
    ammonium_mol_n_per_step: f64,
    ammonia_mol_n_per_step: f64,
    urea_mol_n_per_step: f64,
    nitrate_mol_n_per_step: f64,
};

pub const Inputs = struct {
    existing: State,
    aqueous_nutrient_rates: AqueousNutrientRates,
    exchange_rates: ExchangeRates,
    phosphate_mineral_rates: PhosphateMineralRates,
    aqueous_cation_rates: AqueousCationRates,
    fertilizer_dissolution: FertilizerDissolution,
    litter_water_volume_m3: f64,
    litter_dry_mass_megagrams: f64,
};

/// Direct source-order translation of SOLUTE.F lines 5186--5215.
///
/// Aqueous/mineral concentration rates use water volume, exchange rates use
/// dry mass, and fertilizer rates transfer existing molar inventories.
pub fn calculateSourceOrder(inputs: Inputs) !State {
    try validateInputs(inputs);
    var state = inputs.existing;
    const water_volume = inputs.litter_water_volume_m3;
    const dry_mass = inputs.litter_dry_mass_megagrams;
    const nutrient = inputs.aqueous_nutrient_rates;
    const exchange = inputs.exchange_rates;
    const mineral = inputs.phosphate_mineral_rates;
    const cation = inputs.aqueous_cation_rates;
    const fertilizer = inputs.fertilizer_dissolution;

    // SOLUTE.F 5186--5189.
    state.aqueous_nutrients.ammonium_mol_n_per_step +=
        nutrient.ammonium_mol_n_per_m3_step * water_volume;
    state.aqueous_nutrients.ammonia_mol_n_per_step +=
        nutrient.ammonia_mol_n_per_m3_step * water_volume;
    state.aqueous_nutrients.hydrogen_phosphate_mol_p_per_step +=
        nutrient.hydrogen_phosphate_mol_p_per_m3_step * water_volume;
    state.aqueous_nutrients.dihydrogen_phosphate_mol_p_per_step +=
        nutrient.dihydrogen_phosphate_mol_p_per_m3_step * water_volume;

    // SOLUTE.F 5190--5197.
    state.exchange.ammonium_mol_n_per_step +=
        exchange.ammonium_mol_n_per_megagram_step * dry_mass;
    state.exchange.hydrogen_mol_per_step +=
        exchange.hydrogen_mol_per_megagram_step * dry_mass;
    state.exchange.aluminum_mol_per_step +=
        exchange.aluminum_mol_per_megagram_step * dry_mass;
    state.exchange.iron_mol_per_step +=
        exchange.iron_mol_per_megagram_step * dry_mass;
    state.exchange.calcium_mol_per_step +=
        exchange.calcium_mol_per_megagram_step * dry_mass;
    state.exchange.magnesium_mol_per_step +=
        exchange.magnesium_mol_per_megagram_step * dry_mass;
    state.exchange.sodium_mol_per_step +=
        exchange.sodium_mol_per_megagram_step * dry_mass;
    state.exchange.potassium_mol_per_step +=
        exchange.potassium_mol_per_megagram_step * dry_mass;

    // SOLUTE.F 5198--5202.
    state.phosphate_minerals.aluminum_phosphate_mol_mineral_per_step +=
        mineral.aluminum_phosphate_mol_mineral_per_m3_step *
        water_volume;
    state.phosphate_minerals.iron_phosphate_mol_mineral_per_step +=
        mineral.iron_phosphate_mol_mineral_per_m3_step * water_volume;
    state.phosphate_minerals.dicalcium_phosphate_mol_mineral_per_step +=
        mineral.dicalcium_phosphate_mol_mineral_per_m3_step *
        water_volume;
    state.phosphate_minerals.hydroxyapatite_mol_mineral_per_step +=
        mineral.hydroxyapatite_mol_mineral_per_m3_step * water_volume;
    state.phosphate_minerals.monocalcium_phosphate_mol_mineral_per_step +=
        mineral.monocalcium_phosphate_mol_mineral_per_m3_step *
        water_volume;

    // SOLUTE.F 5203--5208.
    state.aqueous_cations.aluminum_mol_per_step +=
        cation.aluminum_mol_per_m3_step * water_volume;
    state.aqueous_cations.iron_mol_per_step +=
        cation.iron_mol_per_m3_step * water_volume;
    state.aqueous_cations.calcium_mol_per_step +=
        cation.calcium_mol_per_m3_step * water_volume;
    state.aqueous_cations.magnesium_mol_per_step +=
        cation.magnesium_mol_per_m3_step * water_volume;
    state.aqueous_cations.sodium_mol_per_step +=
        cation.sodium_mol_per_m3_step * water_volume;
    state.aqueous_cations.potassium_mol_per_step +=
        cation.potassium_mol_per_m3_step * water_volume;

    // SOLUTE.F 5209--5212.
    state.fertilizer.ammonium_mol_n -= fertilizer.ammonium_mol_n_per_step;
    state.fertilizer.ammonia_mol_n -= fertilizer.ammonia_mol_n_per_step;
    state.fertilizer.urea_mol_n -= fertilizer.urea_mol_n_per_step;
    state.fertilizer.nitrate_mol_n -= fertilizer.nitrate_mol_n_per_step;

    // SOLUTE.F 5213--5215.
    state.aqueous_nutrients.ammonium_mol_n_per_step +=
        fertilizer.ammonium_mol_n_per_step;
    state.aqueous_nutrients.ammonia_mol_n_per_step +=
        fertilizer.ammonia_mol_n_per_step +
        fertilizer.urea_mol_n_per_step;
    state.aqueous_nutrients.nitrate_mol_n_per_step +=
        fertilizer.nitrate_mol_n_per_step;

    try validateResult(state);
    return state;
}

fn validateInputs(inputs: Inputs) !void {
    try validateState(
        inputs.existing,
        error.InvalidSurfaceLitterNutrientAccumulationInput,
    );
    try validateFiniteStruct(
        inputs.aqueous_nutrient_rates,
        error.InvalidSurfaceLitterNutrientAccumulationInput,
    );
    try validateFiniteStruct(
        inputs.exchange_rates,
        error.InvalidSurfaceLitterNutrientAccumulationInput,
    );
    try validateFiniteStruct(
        inputs.phosphate_mineral_rates,
        error.InvalidSurfaceLitterNutrientAccumulationInput,
    );
    try validateFiniteStruct(
        inputs.aqueous_cation_rates,
        error.InvalidSurfaceLitterNutrientAccumulationInput,
    );
    try validateNonnegativeStruct(
        inputs.fertilizer_dissolution,
        error.InvalidSurfaceLitterNutrientAccumulationInput,
    );
    inline for (.{
        inputs.litter_water_volume_m3,
        inputs.litter_dry_mass_megagrams,
    }) |value| {
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterNutrientAccumulationInput;
    }
    if (inputs.fertilizer_dissolution.ammonium_mol_n_per_step >
        inputs.existing.fertilizer.ammonium_mol_n or
        inputs.fertilizer_dissolution.ammonia_mol_n_per_step >
            inputs.existing.fertilizer.ammonia_mol_n or
        inputs.fertilizer_dissolution.urea_mol_n_per_step >
            inputs.existing.fertilizer.urea_mol_n or
        inputs.fertilizer_dissolution.nitrate_mol_n_per_step >
            inputs.existing.fertilizer.nitrate_mol_n)
    {
        return error.SurfaceLitterFertilizerDissolutionExceedsInventory;
    }
}

fn validateResult(state: State) !void {
    try validateState(
        state,
        error.NonFiniteSurfaceLitterNutrientAccumulationResult,
    );
}

fn validateState(state: State, failure: anyerror) !void {
    try validateFiniteStruct(state.aqueous_nutrients, failure);
    try validateFiniteStruct(state.exchange, failure);
    try validateFiniteStruct(state.phosphate_minerals, failure);
    try validateFiniteStruct(state.aqueous_cations, failure);
    try validateNonnegativeStruct(state.fertilizer, failure);
}

fn validateFiniteStruct(value: anytype, failure: anyerror) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(value, field.name)))
            return failure;
    }
}

fn validateNonnegativeStruct(value: anytype, failure: anyerror) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        const field_value = @field(value, field.name);
        if (!std.math.isFinite(field_value) or field_value < 0)
            return failure;
    }
}

fn filled(comptime T: type, start: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields, 0..) |field, index|
        @field(result, field.name) = start + @as(f64, @floatFromInt(index));
    return result;
}

fn testInputs() Inputs {
    return .{
        .existing = .{
            .aqueous_nutrients = filled(AqueousNutrientTotals, 1),
            .exchange = filled(ExchangeTotals, 10),
            .phosphate_minerals = filled(PhosphateMineralTotals, 20),
            .aqueous_cations = filled(AqueousCationTotals, 30),
            .fertilizer = .{
                .ammonium_mol_n = 2,
                .ammonia_mol_n = 3,
                .urea_mol_n = 4,
                .nitrate_mol_n = 5,
            },
        },
        .aqueous_nutrient_rates = .{
            .ammonium_mol_n_per_m3_step = 0.1,
            .ammonia_mol_n_per_m3_step = 0.2,
            .hydrogen_phosphate_mol_p_per_m3_step = 0.3,
            .dihydrogen_phosphate_mol_p_per_m3_step = 0.4,
        },
        .exchange_rates = filled(ExchangeRates, 0.01),
        .phosphate_mineral_rates = filled(PhosphateMineralRates, 0.02),
        .aqueous_cation_rates = filled(AqueousCationRates, 0.03),
        .fertilizer_dissolution = .{
            .ammonium_mol_n_per_step = 0.5,
            .ammonia_mol_n_per_step = 0.6,
            .urea_mol_n_per_step = 0.7,
            .nitrate_mol_n_per_step = 0.8,
        },
        .litter_water_volume_m3 = 2,
        .litter_dry_mass_megagrams = 3,
    };
}

test "SOLUTE nutrient accumulation preserves source dimensional updates" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const initial = inputs.existing;
    const volume = inputs.litter_water_volume_m3;
    const mass = inputs.litter_dry_mass_megagrams;

    try std.testing.expectEqual(
        initial.aqueous_nutrients.hydrogen_phosphate_mol_p_per_step +
            inputs.aqueous_nutrient_rates
                .hydrogen_phosphate_mol_p_per_m3_step * volume,
        result.aqueous_nutrients.hydrogen_phosphate_mol_p_per_step,
    );
    try std.testing.expectEqual(
        initial.exchange.aluminum_mol_per_step +
            inputs.exchange_rates.aluminum_mol_per_megagram_step * mass,
        result.exchange.aluminum_mol_per_step,
    );
    try std.testing.expectEqual(
        initial.phosphate_minerals
            .hydroxyapatite_mol_mineral_per_step +
            inputs.phosphate_mineral_rates
                .hydroxyapatite_mol_mineral_per_m3_step * volume,
        result.phosphate_minerals
            .hydroxyapatite_mol_mineral_per_step,
    );
    try std.testing.expectEqual(
        initial.aqueous_cations.calcium_mol_per_step +
            inputs.aqueous_cation_rates.calcium_mol_per_m3_step * volume,
        result.aqueous_cations.calcium_mol_per_step,
    );
}

test "fertilizer inventory transfer conserves nitrogen" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const initial = inputs.existing;
    const volume = inputs.litter_water_volume_m3;
    const ammonium_fertilizer_gain =
        result.aqueous_nutrients.ammonium_mol_n_per_step -
        initial.aqueous_nutrients.ammonium_mol_n_per_step -
        inputs.aqueous_nutrient_rates.ammonium_mol_n_per_m3_step * volume;
    const ammonia_fertilizer_gain =
        result.aqueous_nutrients.ammonia_mol_n_per_step -
        initial.aqueous_nutrients.ammonia_mol_n_per_step -
        inputs.aqueous_nutrient_rates.ammonia_mol_n_per_m3_step * volume;

    try std.testing.expectApproxEqAbs(
        initial.fertilizer.ammonium_mol_n -
            result.fertilizer.ammonium_mol_n,
        ammonium_fertilizer_gain,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        initial.fertilizer.ammonia_mol_n -
            result.fertilizer.ammonia_mol_n +
            initial.fertilizer.urea_mol_n -
            result.fertilizer.urea_mol_n,
        ammonia_fertilizer_gain,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        initial.fertilizer.nitrate_mol_n -
            result.fertilizer.nitrate_mol_n,
        result.aqueous_nutrients.nitrate_mol_n_per_step -
            initial.aqueous_nutrients.nitrate_mol_n_per_step,
        1.0e-15,
    );
}

test "zero litter dimensions retain only fertilizer transfers" {
    var inputs = testInputs();
    inputs.litter_water_volume_m3 = 0;
    inputs.litter_dry_mass_megagrams = 0;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectEqualDeep(
        inputs.existing.exchange,
        result.exchange,
    );
    try std.testing.expectEqualDeep(
        inputs.existing.phosphate_minerals,
        result.phosphate_minerals,
    );
    try std.testing.expectEqualDeep(
        inputs.existing.aqueous_cations,
        result.aqueous_cations,
    );
    try std.testing.expect(
        result.aqueous_nutrients.ammonium_mol_n_per_step >
            inputs.existing.aqueous_nutrients.ammonium_mol_n_per_step,
    );
}

test "nutrient accumulation rejects invalid inventory and overflow" {
    var inputs = testInputs();
    inputs.litter_dry_mass_megagrams = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterNutrientAccumulationInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.fertilizer_dissolution.ammonium_mol_n_per_step = 3;
    try std.testing.expectError(
        error.SurfaceLitterFertilizerDissolutionExceedsInventory,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.exchange_rates.calcium_mol_per_megagram_step = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterNutrientAccumulationInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.litter_water_volume_m3 = std.math.floatMax(f64);
    inputs.aqueous_cation_rates.calcium_mol_per_m3_step =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterNutrientAccumulationResult,
        calculateSourceOrder(inputs),
    );
}
