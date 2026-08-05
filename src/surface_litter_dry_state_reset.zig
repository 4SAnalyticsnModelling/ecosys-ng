const std = @import("std");

pub const Inputs = struct {
    pH: f64,
    water_dissociation_product_mol2_per_m6: f64,
    mol_per_liter_to_mol_per_m3: f64,
    dry_state_mass_per_water_sentinel_megagrams_per_m3: f64,
};

pub const FertilizerDissolutionRates = struct {
    ammonium_mol_n_per_step: f64,
    ammonia_mol_n_per_step: f64,
    urea_mol_n_per_step: f64,
    nitrate_mol_n_per_step: f64,
};

pub const MineralReactionRates = struct {
    aluminum_hydroxide_mol_per_m3_step: f64,
    iron_hydroxide_mol_per_m3_step: f64,
    calcite_mol_per_m3_step: f64,
    gypsum_mol_per_m3_step: f64,
    aluminum_phosphate_mol_per_m3_step: f64,
    iron_phosphate_mol_per_m3_step: f64,
    dicalcium_phosphate_mol_per_m3_step: f64,
    hydroxyapatite_mol_per_m3_step: f64,
    monocalcium_phosphate_mol_per_m3_step: f64,
};

pub const NitrogenState = struct {
    ammonium_input_g_n_per_step: f64,
    ammonia_input_g_n_per_step: f64,
    ammonium_concentration_mol_n_per_m3: f64,
    ammonia_concentration_mol_n_per_m3: f64,
    exchangeable_ammonium_mol_n_per_megagram: f64,
    ammonium_exchange_mol_n_per_megagram_step: f64,
    ammonium_association_mol_n_per_m3_step: f64,
};

pub const PhosphorusState = struct {
    hydrogen_phosphate_input_mol_p_per_m3_step: f64,
    dihydrogen_phosphate_input_mol_p_per_m3_step: f64,
    hydrogen_phosphate_concentration_mol_p_per_m3: f64,
    dihydrogen_phosphate_concentration_mol_p_per_m3: f64,
    dihydrogen_phosphate_association_mol_p_per_m3_step: f64,
};

pub const CationExchangeRates = struct {
    hydrogen_mol_per_megagram_step: f64,
    aluminum_mol_per_megagram_step: f64,
    iron_mol_per_megagram_step: f64,
    calcium_mol_per_megagram_step: f64,
    magnesium_mol_per_megagram_step: f64,
    sodium_mol_per_megagram_step: f64,
    potassium_mol_per_megagram_step: f64,
    carboxyl_hydrogen_mol_per_megagram_step: f64,
};

pub const CarbonateReactionRates = struct {
    bicarbonate_association_mol_c_per_m3_step: f64,
    carbon_dioxide_association_mol_c_per_m3_step: f64,
};

pub const SolidConcentrations = struct {
    aluminum_hydroxide_mol_per_m3: f64,
    iron_hydroxide_mol_per_m3: f64,
    calcite_mol_per_m3: f64,
    gypsum_mol_per_m3: f64,
};

pub const Result = struct {
    litter_dry_mass_megagrams: f64,
    hydrogen_hydroxide_equilibration_mol_per_m3_step: f64,
    external_hydrogen_mol_per_m3_step: f64,
    hydrogen_activity_mol_per_m3: f64,
    hydroxide_activity_mol_per_m3: f64,
    fertilizer: FertilizerDissolutionRates,
    minerals: MineralReactionRates,
    nitrogen: NitrogenState,
    phosphorus: PhosphorusState,
    cation_exchange: CationExchangeRates,
    carbonate: CarbonateReactionRates,
    solids: SolidConcentrations,
    litter_mass_per_water_volume_megagrams_per_m3: f64,
};

/// Direct source-order translation of SOLUTE.F lines 4931--4976.
///
/// The enclosing wet-litter test has failed. This pure one-cell initializer
/// returns the explicit dry-state values used by subsequent transformations.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);

    // SOLUTE.F 4931--4935.
    const litter_dry_mass_megagrams = 0.0;
    const hydrogen_hydroxide_equilibration = 0.0;
    const external_hydrogen = 0.0;
    const hydrogen_activity = std.math.pow(f64, 10.0, -inputs.pH) *
        inputs.mol_per_liter_to_mol_per_m3;
    const hydroxide_activity =
        inputs.water_dissociation_product_mol2_per_m6 /
        hydrogen_activity;

    // SOLUTE.F 4936--4975. Declare every scalar in source assignment order;
    // domain structs are assembled only after the final assignment.
    const broadcast_ammonium_dissolution = 0.0;
    const broadcast_ammonia_dissolution = 0.0;
    const urea_dissolution = 0.0;
    const nitrate_dissolution = 0.0;
    const aluminum_hydroxide_reaction = 0.0;
    const iron_hydroxide_reaction = 0.0;
    const calcite_reaction = 0.0;
    const gypsum_reaction = 0.0;
    const aluminum_phosphate_reaction = 0.0;
    const iron_phosphate_reaction = 0.0;
    const dicalcium_phosphate_reaction = 0.0;
    const hydroxyapatite_reaction = 0.0;
    const monocalcium_phosphate_reaction = 0.0;
    const ammonium_input = 0.0;
    const ammonia_input = 0.0;
    const ammonium_concentration = 0.0;
    const ammonia_concentration = 0.0;
    const exchangeable_ammonium = 0.0;
    const hydrogen_phosphate_input = 0.0;
    const dihydrogen_phosphate_input = 0.0;
    const hydrogen_phosphate_concentration = 0.0;
    const dihydrogen_phosphate_concentration = 0.0;
    const ammonium_exchange = 0.0;
    const ammonium_association = 0.0;
    const bicarbonate_association = 0.0;
    const carbon_dioxide_association = 0.0;
    const dihydrogen_phosphate_association = 0.0;
    const hydrogen_exchange = 0.0;
    const aluminum_exchange = 0.0;
    const iron_exchange = 0.0;
    const calcium_exchange = 0.0;
    const magnesium_exchange = 0.0;
    const sodium_exchange = 0.0;
    const potassium_exchange = 0.0;
    const carboxyl_hydrogen_exchange = 0.0;
    const aluminum_hydroxide_concentration = 0.0;
    const iron_hydroxide_concentration = 0.0;
    const calcite_concentration = 0.0;
    const gypsum_concentration = 0.0;
    const litter_mass_per_water =
        inputs.dry_state_mass_per_water_sentinel_megagrams_per_m3;

    const fertilizer: FertilizerDissolutionRates = .{
        .ammonium_mol_n_per_step = broadcast_ammonium_dissolution,
        .ammonia_mol_n_per_step = broadcast_ammonia_dissolution,
        .urea_mol_n_per_step = urea_dissolution,
        .nitrate_mol_n_per_step = nitrate_dissolution,
    };
    const minerals: MineralReactionRates = .{
        .aluminum_hydroxide_mol_per_m3_step = aluminum_hydroxide_reaction,
        .iron_hydroxide_mol_per_m3_step = iron_hydroxide_reaction,
        .calcite_mol_per_m3_step = calcite_reaction,
        .gypsum_mol_per_m3_step = gypsum_reaction,
        .aluminum_phosphate_mol_per_m3_step = aluminum_phosphate_reaction,
        .iron_phosphate_mol_per_m3_step = iron_phosphate_reaction,
        .dicalcium_phosphate_mol_per_m3_step = dicalcium_phosphate_reaction,
        .hydroxyapatite_mol_per_m3_step = hydroxyapatite_reaction,
        .monocalcium_phosphate_mol_per_m3_step = monocalcium_phosphate_reaction,
    };
    const nitrogen: NitrogenState = .{
        .ammonium_input_g_n_per_step = ammonium_input,
        .ammonia_input_g_n_per_step = ammonia_input,
        .ammonium_concentration_mol_n_per_m3 = ammonium_concentration,
        .ammonia_concentration_mol_n_per_m3 = ammonia_concentration,
        .exchangeable_ammonium_mol_n_per_megagram = exchangeable_ammonium,
        .ammonium_exchange_mol_n_per_megagram_step = ammonium_exchange,
        .ammonium_association_mol_n_per_m3_step = ammonium_association,
    };
    const phosphorus: PhosphorusState = .{
        .hydrogen_phosphate_input_mol_p_per_m3_step = hydrogen_phosphate_input,
        .dihydrogen_phosphate_input_mol_p_per_m3_step = dihydrogen_phosphate_input,
        .hydrogen_phosphate_concentration_mol_p_per_m3 = hydrogen_phosphate_concentration,
        .dihydrogen_phosphate_concentration_mol_p_per_m3 = dihydrogen_phosphate_concentration,
        .dihydrogen_phosphate_association_mol_p_per_m3_step = dihydrogen_phosphate_association,
    };
    const carbonate: CarbonateReactionRates = .{
        .bicarbonate_association_mol_c_per_m3_step = bicarbonate_association,
        .carbon_dioxide_association_mol_c_per_m3_step = carbon_dioxide_association,
    };
    const cation_exchange: CationExchangeRates = .{
        .hydrogen_mol_per_megagram_step = hydrogen_exchange,
        .aluminum_mol_per_megagram_step = aluminum_exchange,
        .iron_mol_per_megagram_step = iron_exchange,
        .calcium_mol_per_megagram_step = calcium_exchange,
        .magnesium_mol_per_megagram_step = magnesium_exchange,
        .sodium_mol_per_megagram_step = sodium_exchange,
        .potassium_mol_per_megagram_step = potassium_exchange,
        .carboxyl_hydrogen_mol_per_megagram_step = carboxyl_hydrogen_exchange,
    };
    const solids: SolidConcentrations = .{
        .aluminum_hydroxide_mol_per_m3 = aluminum_hydroxide_concentration,
        .iron_hydroxide_mol_per_m3 = iron_hydroxide_concentration,
        .calcite_mol_per_m3 = calcite_concentration,
        .gypsum_mol_per_m3 = gypsum_concentration,
    };
    const result: Result = .{
        .litter_dry_mass_megagrams = litter_dry_mass_megagrams,
        .hydrogen_hydroxide_equilibration_mol_per_m3_step = hydrogen_hydroxide_equilibration,
        .external_hydrogen_mol_per_m3_step = external_hydrogen,
        .hydrogen_activity_mol_per_m3 = hydrogen_activity,
        .hydroxide_activity_mol_per_m3 = hydroxide_activity,
        .fertilizer = fertilizer,
        .minerals = minerals,
        .nitrogen = nitrogen,
        .phosphorus = phosphorus,
        .cation_exchange = cation_exchange,
        .carbonate = carbonate,
        .solids = solids,
        .litter_mass_per_water_volume_megagrams_per_m3 = litter_mass_per_water,
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidSurfaceLitterDryStateInput;
    }
    if (inputs.water_dissociation_product_mol2_per_m6 <= 0 or
        inputs.mol_per_liter_to_mol_per_m3 <= 0 or
        inputs.dry_state_mass_per_water_sentinel_megagrams_per_m3 <= 0)
    {
        return error.InvalidSurfaceLitterDryStateInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (.{
        result.litter_dry_mass_megagrams,
        result.hydrogen_hydroxide_equilibration_mol_per_m3_step,
        result.external_hydrogen_mol_per_m3_step,
        result.hydrogen_activity_mol_per_m3,
        result.hydroxide_activity_mol_per_m3,
        result.litter_mass_per_water_volume_megagrams_per_m3,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterDryStateResult;
    }
    if (result.hydrogen_activity_mol_per_m3 <= 0 or
        result.hydroxide_activity_mol_per_m3 <= 0)
    {
        return error.InvalidSurfaceLitterDryStateResult;
    }
}

fn testInputs() Inputs {
    return .{
        .pH = 7,
        .water_dissociation_product_mol2_per_m6 = 1.0e-8,
        .mol_per_liter_to_mol_per_m3 = 1.0e3,
        .dry_state_mass_per_water_sentinel_megagrams_per_m3 = 1,
    };
}

test "SOLUTE dry surface reset preserves pH activity source expressions" {
    const inputs = testInputs();
    const result = try calculateSourceOrder(inputs);
    const expected_hydrogen =
        std.math.pow(f64, 10.0, -inputs.pH) *
        inputs.mol_per_liter_to_mol_per_m3;

    try std.testing.expectEqual(
        expected_hydrogen,
        result.hydrogen_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        inputs.water_dissociation_product_mol2_per_m6 /
            expected_hydrogen,
        result.hydroxide_activity_mol_per_m3,
    );
    try std.testing.expectEqual(
        inputs.water_dissociation_product_mol2_per_m6,
        result.hydrogen_activity_mol_per_m3 *
            result.hydroxide_activity_mol_per_m3,
    );
}

test "dry surface reset clears every source state group" {
    const result = try calculateSourceOrder(testInputs());

    try std.testing.expectEqual(@as(f64, 0), result.litter_dry_mass_megagrams);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.hydrogen_hydroxide_equilibration_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.external_hydrogen_mol_per_m3_step,
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(FertilizerDissolutionRates),
        result.fertilizer,
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(MineralReactionRates),
        result.minerals,
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(NitrogenState),
        result.nitrogen,
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(PhosphorusState),
        result.phosphorus,
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(CationExchangeRates),
        result.cation_exchange,
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(CarbonateReactionRates),
        result.carbonate,
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(SolidConcentrations),
        result.solids,
    );
}

test "dry surface reset uses runtime mass-to-water sentinel" {
    var inputs = testInputs();
    inputs.dry_state_mass_per_water_sentinel_megagrams_per_m3 = 2.5;
    const result = try calculateSourceOrder(inputs);

    try std.testing.expectEqual(
        @as(f64, 2.5),
        result.litter_mass_per_water_volume_megagrams_per_m3,
    );
}

test "dry surface reset rejects invalid inputs and pH overflow" {
    var inputs = testInputs();
    inputs.water_dissociation_product_mol2_per_m6 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterDryStateInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.pH = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterDryStateInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.pH = -400;
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterDryStateResult,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.pH = 400;
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterDryStateResult,
        calculateSourceOrder(inputs),
    );
}
