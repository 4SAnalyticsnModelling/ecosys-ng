const std = @import("std");
const aqueous_network = @import("solute_aqueous_network.zig");

const minimum_concentration_mol_per_m3 = 1.0e-32;

/// Direct translation of the SOLUTE.F lines 158--161 layer admission gate.
/// DLYR(3,L,NY,NX) is vertical layer thickness (m); VOLW is water volume
/// (m3). Equality to either threshold remains inactive, matching `.GT.`.
pub fn isReactiveSoilLayer(
    vertical_layer_thickness_m: f64,
    water_volume_m3: f64,
    minimum_layer_thickness_m: f64,
    minimum_water_volume_m3: f64,
) !bool {
    inline for (.{
        vertical_layer_thickness_m,
        water_volume_m3,
        minimum_layer_thickness_m,
        minimum_water_volume_m3,
    }) |value| {
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoluteLayerAdmissionInput;
    }
    return vertical_layer_thickness_m > minimum_layer_thickness_m and
        water_volume_m3 > minimum_water_volume_m3;
}

test "SOLUTE lines 158-161 require both layer extents above thresholds" {
    try std.testing.expect(try isReactiveSoilLayer(0.2, 0.03, 0.01, 1.0e-12));
    try std.testing.expect(!try isReactiveSoilLayer(0.01, 0.03, 0.01, 1.0e-12));
    try std.testing.expect(!try isReactiveSoilLayer(0.2, 1.0e-12, 0.01, 1.0e-12));
    try std.testing.expect(!try isReactiveSoilLayer(0.005, 0, 0.01, 1.0e-12));
}

test "SOLUTE layer admission rejects non-finite runtime geometry" {
    try std.testing.expectError(
        error.InvalidSoluteLayerAdmissionInput,
        isReactiveSoilLayer(std.math.nan(f64), 1, 0, 0),
    );
}

pub const NitrateZoneInput = struct {
    nitrate_g_n: f64,
    water_volume_m3: f64,
};

pub const NitrateConcentrations = struct {
    non_band_mol_n_per_m3: f64,
    band_mol_n_per_m3: f64,
};

pub const ReactionIterationMode = enum {
    restricted_chemistry,
    dynamic_salt_chemistry,
};

pub const NonBandAmmoniumInputs = struct {
    water_volume_m3: f64,
    soil_mass_megagrams: f64,
    soluble_ammonium_g_n: f64,
    soluble_ammonia_g_n: f64,
    exchangeable_ammonium_mol_n: f64,
    fertilizer_ammonium_input_g_n_per_step: f64,
    root_ammonium_uptake_g_n_per_step: f64,
    ammonium_reaction_source_mol_n_per_step: f64,
    root_ammonia_uptake_g_n_per_step: f64,
    ammonia_reaction_source_mol_n_per_step: f64,
};

pub const BandAmmoniumInputs = struct {
    water_volume_m3: f64,
    soil_mass_megagrams: f64,
    soluble_ammonium_g_n: f64,
    soluble_ammonia_g_n: f64,
    exchangeable_ammonium_mol_n: f64,
    fertilizer_ammonium_input_g_n_per_step: f64,
    root_ammonium_uptake_g_n_per_step: f64,
    ammonium_reaction_source_first_mol_n_per_step: f64,
    ammonium_reaction_source_second_mol_n_per_step: f64,
    root_ammonia_uptake_g_n_per_step: f64,
    ammonia_reaction_source_first_mol_n_per_step: f64,
    ammonia_reaction_source_second_mol_n_per_step: f64,
};

pub const AmmoniumZoneConcentrations = struct {
    ammonium_source_mol_n_per_m3_iteration: f64,
    ammonia_source_mol_n_per_m3_iteration: f64,
    ammonium_mol_n_per_m3: f64,
    ammonia_mol_n_per_m3: f64,
    exchangeable_ammonium_mol_n_per_megagram: f64,
    soil_mass_per_water_volume_megagrams_per_m3: f64,
};

pub const AmmoniumConcentrations = struct {
    non_band: AmmoniumZoneConcentrations,
    band: AmmoniumZoneConcentrations,
};

/// Direct source-order translation of SOLUTE.F lines 367--414.
/// Nitrogen amounts stored in g N are converted using 14 g N mol-1. Reaction
/// sources are already mol N step-1. The returned source increments and
/// concentrations are mol N m-3 per nonlinear iteration and mol N m-3.
pub fn initializeAmmoniumConcentrations(
    non_band: NonBandAmmoniumInputs,
    band: BandAmmoniumInputs,
    mode: ReactionIterationMode,
    maximum_equilibrium_iterations: usize,
    minimum_water_volume_m3: f64,
    concentration_floor_mol_per_m3: f64,
) !AmmoniumConcentrations {
    try validateAmmoniumInputs(
        non_band,
        band,
        maximum_equilibrium_iterations,
        minimum_water_volume_m3,
        concentration_floor_mol_per_m3,
    );
    const reaction_iterations: f64 = switch (mode) {
        .dynamic_salt_chemistry => @floatFromInt(maximum_equilibrium_iterations),
        .restricted_chemistry => 1.0,
    };

    const non_band_result = if (non_band.water_volume_m3 > minimum_water_volume_m3)
        initializeNonBandAmmonium(
            non_band,
            reaction_iterations,
            concentration_floor_mol_per_m3,
        )
    else
        zeroAmmoniumZone();
    const band_result = if (band.water_volume_m3 > minimum_water_volume_m3)
        initializeBandAmmonium(
            band,
            reaction_iterations,
            concentration_floor_mol_per_m3,
        )
    else
        zeroAmmoniumZone();
    return .{ .non_band = non_band_result, .band = band_result };
}

fn initializeNonBandAmmonium(
    input: NonBandAmmoniumInputs,
    reaction_iterations: f64,
    minimum_concentration: f64,
) AmmoniumZoneConcentrations {
    const nitrogen_water_volume_g_n_m3 = 14.0 * input.water_volume_m3;
    const ammonium_source =
        (input.fertilizer_ammonium_input_g_n_per_step -
            input.root_ammonium_uptake_g_n_per_step +
            14.0 * input.ammonium_reaction_source_mol_n_per_step) /
        (reaction_iterations * nitrogen_water_volume_g_n_m3);
    const ammonia_source =
        (-input.root_ammonia_uptake_g_n_per_step +
            14.0 * input.ammonia_reaction_source_mol_n_per_step) /
        (reaction_iterations * nitrogen_water_volume_g_n_m3);
    return .{
        .ammonium_source_mol_n_per_m3_iteration = ammonium_source,
        .ammonia_source_mol_n_per_m3_iteration = ammonia_source,
        .ammonium_mol_n_per_m3 = @max(
            minimum_concentration,
            input.soluble_ammonium_g_n / nitrogen_water_volume_g_n_m3 +
                ammonium_source,
        ),
        .ammonia_mol_n_per_m3 = @max(
            minimum_concentration,
            input.soluble_ammonia_g_n / nitrogen_water_volume_g_n_m3 +
                ammonia_source,
        ),
        .exchangeable_ammonium_mol_n_per_megagram = @max(
            minimum_concentration,
            input.exchangeable_ammonium_mol_n / input.soil_mass_megagrams,
        ),
        .soil_mass_per_water_volume_megagrams_per_m3 = input.soil_mass_megagrams / input.water_volume_m3,
    };
}

fn initializeBandAmmonium(
    input: BandAmmoniumInputs,
    reaction_iterations: f64,
    minimum_concentration: f64,
) AmmoniumZoneConcentrations {
    const nitrogen_water_volume_g_n_m3 = 14.0 * input.water_volume_m3;
    const ammonium_source =
        (input.fertilizer_ammonium_input_g_n_per_step -
            input.root_ammonium_uptake_g_n_per_step +
            14.0 * (input.ammonium_reaction_source_first_mol_n_per_step +
                input.ammonium_reaction_source_second_mol_n_per_step)) /
        (reaction_iterations * nitrogen_water_volume_g_n_m3);
    const ammonia_source =
        (-input.root_ammonia_uptake_g_n_per_step +
            14.0 * (input.ammonia_reaction_source_first_mol_n_per_step +
                input.ammonia_reaction_source_second_mol_n_per_step)) /
        (reaction_iterations * nitrogen_water_volume_g_n_m3);
    return .{
        .ammonium_source_mol_n_per_m3_iteration = ammonium_source,
        .ammonia_source_mol_n_per_m3_iteration = ammonia_source,
        .ammonium_mol_n_per_m3 = @max(
            minimum_concentration,
            input.soluble_ammonium_g_n / nitrogen_water_volume_g_n_m3 +
                ammonium_source,
        ),
        .ammonia_mol_n_per_m3 = @max(
            minimum_concentration,
            input.soluble_ammonia_g_n / nitrogen_water_volume_g_n_m3 +
                ammonia_source,
        ),
        .exchangeable_ammonium_mol_n_per_megagram = @max(
            minimum_concentration,
            input.exchangeable_ammonium_mol_n / input.soil_mass_megagrams,
        ),
        .soil_mass_per_water_volume_megagrams_per_m3 = input.soil_mass_megagrams / input.water_volume_m3,
    };
}

fn zeroAmmoniumZone() AmmoniumZoneConcentrations {
    return .{
        .ammonium_source_mol_n_per_m3_iteration = 0,
        .ammonia_source_mol_n_per_m3_iteration = 0,
        .ammonium_mol_n_per_m3 = 0,
        .ammonia_mol_n_per_m3 = 0,
        .exchangeable_ammonium_mol_n_per_megagram = 0,
        .soil_mass_per_water_volume_megagrams_per_m3 = 0,
    };
}

fn validateAmmoniumInputs(
    non_band: NonBandAmmoniumInputs,
    band: BandAmmoniumInputs,
    maximum_equilibrium_iterations: usize,
    minimum_water_volume_m3: f64,
    concentration_floor_mol_per_m3: f64,
) !void {
    if (maximum_equilibrium_iterations == 0 or
        !std.math.isFinite(minimum_water_volume_m3) or
        minimum_water_volume_m3 < 0 or
        !std.math.isFinite(concentration_floor_mol_per_m3) or
        concentration_floor_mol_per_m3 <= 0)
        return error.InvalidAmmoniumConcentrationInput;
    inline for (@typeInfo(NonBandAmmoniumInputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(non_band, field.name)))
            return error.InvalidAmmoniumConcentrationInput;
    inline for (@typeInfo(BandAmmoniumInputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(band, field.name)))
            return error.InvalidAmmoniumConcentrationInput;
    if (non_band.water_volume_m3 < 0 or band.water_volume_m3 < 0 or
        (non_band.water_volume_m3 > minimum_water_volume_m3 and
            non_band.soil_mass_megagrams <= 0) or
        (band.water_volume_m3 > minimum_water_volume_m3 and
            band.soil_mass_megagrams <= 0))
        return error.InvalidAmmoniumConcentrationInput;
}

/// SOLUTE lines 526--542. The source stores nitrate as g N and divides by
/// 14 g N mol-1 before converting each runtime water-volume zone.
pub fn initializeNitrateConcentrations(
    non_band: NitrateZoneInput,
    band: NitrateZoneInput,
    minimum_water_volume_m3: f64,
    nitrogen_molar_mass_g_per_mol: f64,
    concentration_floor_mol_n_per_m3: f64,
) !NitrateConcentrations {
    if (!std.math.isFinite(minimum_water_volume_m3) or
        minimum_water_volume_m3 < 0 or
        !std.math.isFinite(nitrogen_molar_mass_g_per_mol) or
        nitrogen_molar_mass_g_per_mol <= 0 or
        !std.math.isFinite(concentration_floor_mol_n_per_m3) or
        concentration_floor_mol_n_per_m3 <= 0)
        return error.InvalidNitrateConcentrationInput;
    return .{
        .non_band_mol_n_per_m3 = try nitrateConcentration(
            non_band,
            minimum_water_volume_m3,
            nitrogen_molar_mass_g_per_mol,
            concentration_floor_mol_n_per_m3,
        ),
        .band_mol_n_per_m3 = try nitrateConcentration(
            band,
            minimum_water_volume_m3,
            nitrogen_molar_mass_g_per_mol,
            concentration_floor_mol_n_per_m3,
        ),
    };
}

fn nitrateConcentration(
    input: NitrateZoneInput,
    minimum_water_volume_m3: f64,
    nitrogen_molar_mass_g_per_mol: f64,
    concentration_floor_mol_n_per_m3: f64,
) !f64 {
    if (!std.math.isFinite(input.nitrate_g_n) or
        !std.math.isFinite(input.water_volume_m3) or
        input.water_volume_m3 < 0)
        return error.InvalidNitrateConcentrationInput;
    if (input.water_volume_m3 <= minimum_water_volume_m3) return 0;
    return @max(
        concentration_floor_mol_n_per_m3,
        input.nitrate_g_n /
            (nitrogen_molar_mass_g_per_mol * input.water_volume_m3),
    );
}

fn testNonBandAmmoniumInputs() NonBandAmmoniumInputs {
    return .{
        .water_volume_m3 = 2,
        .soil_mass_megagrams = 4,
        .soluble_ammonium_g_n = 28,
        .soluble_ammonia_g_n = 14,
        .exchangeable_ammonium_mol_n = 8,
        .fertilizer_ammonium_input_g_n_per_step = 14,
        .root_ammonium_uptake_g_n_per_step = 7,
        .ammonium_reaction_source_mol_n_per_step = 0.5,
        .root_ammonia_uptake_g_n_per_step = 7,
        .ammonia_reaction_source_mol_n_per_step = 1,
    };
}

fn testBandAmmoniumInputs() BandAmmoniumInputs {
    return .{
        .water_volume_m3 = 1,
        .soil_mass_megagrams = 3,
        .soluble_ammonium_g_n = 14,
        .soluble_ammonia_g_n = 7,
        .exchangeable_ammonium_mol_n = 6,
        .fertilizer_ammonium_input_g_n_per_step = 28,
        .root_ammonium_uptake_g_n_per_step = 14,
        .ammonium_reaction_source_first_mol_n_per_step = 0.25,
        .ammonium_reaction_source_second_mol_n_per_step = 0.75,
        .root_ammonia_uptake_g_n_per_step = 7,
        .ammonia_reaction_source_first_mol_n_per_step = 0.25,
        .ammonia_reaction_source_second_mol_n_per_step = 0.25,
    };
}

test "SOLUTE 367-414 preserves dynamic-salt zone calculation order" {
    const result = try initializeAmmoniumConcentrations(
        testNonBandAmmoniumInputs(),
        testBandAmmoniumInputs(),
        .dynamic_salt_chemistry,
        60,
        1.0e-12,
        1.0e-32,
    );
    const expected_non_band_ammonium_source =
        (14.0 - 7.0 + 14.0 * 0.5) / (60.0 * (14.0 * 2.0));
    const expected_band_ammonium_source =
        (28.0 - 14.0 + 14.0 * (0.25 + 0.75)) /
        (60.0 * (14.0 * 1.0));
    try std.testing.expectEqual(
        expected_non_band_ammonium_source,
        result.non_band.ammonium_source_mol_n_per_m3_iteration,
    );
    try std.testing.expectEqual(
        28.0 / (14.0 * 2.0) + expected_non_band_ammonium_source,
        result.non_band.ammonium_mol_n_per_m3,
    );
    try std.testing.expectEqual(
        expected_band_ammonium_source,
        result.band.ammonium_source_mol_n_per_m3_iteration,
    );
    try std.testing.expectEqual(
        @as(f64, 2),
        result.non_band.exchangeable_ammonium_mol_n_per_megagram,
    );
    try std.testing.expectEqual(
        @as(f64, 2),
        result.non_band.soil_mass_per_water_volume_megagrams_per_m3,
    );
}

test "SOLUTE 367-414 restricted chemistry uses one reaction iteration" {
    const result = try initializeAmmoniumConcentrations(
        testNonBandAmmoniumInputs(),
        testBandAmmoniumInputs(),
        .restricted_chemistry,
        60,
        1.0e-12,
        1.0e-32,
    );
    try std.testing.expectEqual(
        (14.0 - 7.0 + 14.0 * 0.5) / (14.0 * 2.0),
        result.non_band.ammonium_source_mol_n_per_m3_iteration,
    );
}

test "SOLUTE 389-414 dry zones reset every derived value" {
    var non_band = testNonBandAmmoniumInputs();
    var band = testBandAmmoniumInputs();
    non_band.water_volume_m3 = 1.0e-13;
    non_band.soil_mass_megagrams = 0;
    band.water_volume_m3 = 0;
    band.soil_mass_megagrams = 0;
    const result = try initializeAmmoniumConcentrations(
        non_band,
        band,
        .dynamic_salt_chemistry,
        60,
        1.0e-12,
        1.0e-32,
    );
    try std.testing.expectEqualDeep(zeroAmmoniumZone(), result.non_band);
    try std.testing.expectEqualDeep(zeroAmmoniumZone(), result.band);
}

test "SOLUTE ammonium initialization rejects active zero soil mass" {
    var non_band = testNonBandAmmoniumInputs();
    non_band.soil_mass_megagrams = 0;
    try std.testing.expectError(
        error.InvalidAmmoniumConcentrationInput,
        initializeAmmoniumConcentrations(
            non_band,
            testBandAmmoniumInputs(),
            .dynamic_salt_chemistry,
            60,
            1.0e-12,
            1.0e-32,
        ),
    );
}

pub const CationExchangeAmounts = struct {
    capacity_mol_charge: f64,
    hydrogen_mol: f64,
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    carboxyl_bound_hydrogen_mol: f64,
};

pub const CationExchangeConcentrations = struct {
    capacity_mol_charge_per_megagram: f64,
    hydrogen_mol_per_megagram: f64,
    aluminum_mol_per_megagram: f64,
    iron_mol_per_megagram: f64,
    calcium_mol_per_megagram: f64,
    magnesium_mol_per_megagram: f64,
    sodium_mol_per_megagram: f64,
    potassium_mol_per_megagram: f64,
    carboxyl_bound_hydrogen_mol_per_megagram: f64,
    total_carboxyl_sites_mol_per_megagram: f64,
    soil_mass_per_water_volume_megagrams_per_m3: f64,
};

/// SOLUTE lines 544--578. `organic_carbon_g_c` is converted to Mg C by
/// 1e-6 before applying the carboxyl-site density in mol Mg-C-1.
pub fn initializeCationExchangeConcentrations(
    amounts: CationExchangeAmounts,
    soil_mass_megagrams: f64,
    water_volume_m3: f64,
    carboxyl_sites_mol_per_megagram_c: f64,
    organic_carbon_g_c: f64,
    minimum_active_extent: f64,
) !CationExchangeConcentrations {
    inline for (@typeInfo(CationExchangeAmounts).@"struct".fields) |field|
        if (!std.math.isFinite(@field(amounts, field.name)))
            return error.InvalidCationExchangeInitialization;
    inline for (.{
        soil_mass_megagrams,
        water_volume_m3,
        carboxyl_sites_mol_per_megagram_c,
        organic_carbon_g_c,
        minimum_active_extent,
    }) |value| try validateFiniteNonnegative(value);

    // Legacy ZEROS(NY,NX) is one cell-specific numerical threshold applied
    // contextually to both the Mg soil mass and m3 water volume.
    if (soil_mass_megagrams <= minimum_active_extent or
        water_volume_m3 <= minimum_active_extent)
        return zeroCationExchangeConcentrations();

    return .{
        .capacity_mol_charge_per_megagram = concentration(amounts.capacity_mol_charge, soil_mass_megagrams),
        .hydrogen_mol_per_megagram = concentration(amounts.hydrogen_mol, soil_mass_megagrams),
        .aluminum_mol_per_megagram = concentration(amounts.aluminum_mol, soil_mass_megagrams),
        .iron_mol_per_megagram = concentration(amounts.iron_mol, soil_mass_megagrams),
        .calcium_mol_per_megagram = concentration(amounts.calcium_mol, soil_mass_megagrams),
        .magnesium_mol_per_megagram = concentration(amounts.magnesium_mol, soil_mass_megagrams),
        .sodium_mol_per_megagram = concentration(amounts.sodium_mol, soil_mass_megagrams),
        .potassium_mol_per_megagram = concentration(amounts.potassium_mol, soil_mass_megagrams),
        .carboxyl_bound_hydrogen_mol_per_megagram = concentration(
            amounts.carboxyl_bound_hydrogen_mol,
            soil_mass_megagrams,
        ),
        .total_carboxyl_sites_mol_per_megagram = concentration(
            carboxyl_sites_mol_per_megagram_c * 1.0e-6 * organic_carbon_g_c,
            soil_mass_megagrams,
        ),
        .soil_mass_per_water_volume_megagrams_per_m3 = soil_mass_megagrams / water_volume_m3,
    };
}

fn concentration(amount_mol: f64, soil_mass_megagrams: f64) f64 {
    return @max(minimum_concentration_mol_per_m3, amount_mol / soil_mass_megagrams);
}

fn zeroCationExchangeConcentrations() CationExchangeConcentrations {
    return .{
        .capacity_mol_charge_per_megagram = 0,
        .hydrogen_mol_per_megagram = 0,
        .aluminum_mol_per_megagram = 0,
        .iron_mol_per_megagram = 0,
        .calcium_mol_per_megagram = 0,
        .magnesium_mol_per_megagram = 0,
        .sodium_mol_per_megagram = 0,
        .potassium_mol_per_megagram = 0,
        .carboxyl_bound_hydrogen_mol_per_megagram = 0,
        .total_carboxyl_sites_mol_per_megagram = 0,
        .soil_mass_per_water_volume_megagrams_per_m3 = 1,
    };
}

fn validateFiniteNonnegative(value: f64) !void {
    if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSoluteConcentrationInitialization;
}

pub const FreeIonAmounts = struct {
    hydrogen_mol: f64,
    hydroxide_mol: f64,
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    sulfate_mol: f64,
    chloride_mol: f64,
    carbonate_mol: f64,
    bicarbonate_mol: f64,
    carbon_dioxide_g_c: f64,
};

pub const DissolvedComplexAmounts = struct {
    aluminum_hydroxide_1_mol: f64,
    aluminum_hydroxide_2_mol: f64,
    aluminum_hydroxide_3_mol: f64,
    aluminum_hydroxide_4_mol: f64,
    aluminum_sulfate_mol: f64,
    iron_hydroxide_1_mol: f64,
    iron_hydroxide_2_mol: f64,
    iron_hydroxide_3_mol: f64,
    iron_hydroxide_4_mol: f64,
    iron_sulfate_mol: f64,
    calcium_hydroxide_mol: f64,
    calcium_carbonate_mol: f64,
    calcium_bicarbonate_mol: f64,
    calcium_sulfate_mol: f64,
    magnesium_hydroxide_mol: f64,
    magnesium_carbonate_mol: f64,
    magnesium_bicarbonate_mol: f64,
    magnesium_sulfate_mol: f64,
    sodium_carbonate_mol: f64,
    sodium_sulfate_mol: f64,
    potassium_sulfate_mol: f64,
    hydrogen_silicate_mol: f64,
};

pub const DynamicSaltAqueousInputs = struct {
    free_ions: FreeIonAmounts,
    dissolved_complexes: DissolvedComplexAmounts,
    water_volume_m3: f64,
    minimum_water_volume_m3: f64,
    dry_soil_ph: f64,
    water_activity_product_mol2_per_m6: f64,
    carbon_molar_mass_g_per_mol: f64,
};

pub const RestrictedSaltAmounts = struct {
    hydrogen_mol: f64,
    hydroxide_mol: f64,
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    sulfate_mol: f64,
    chloride_mol: f64,
    aluminum_hydroxide_2_mol: f64,
    iron_hydroxide_2_mol: f64,
};

pub const RestrictedSaltInputs = struct {
    amounts: RestrictedSaltAmounts,
    water_volume_m3: f64,
    exchange_soil_mass_megagrams: f64,
    cation_exchange_capacity_mol_charge: f64,
    carbon_dioxide_g_c_per_m3: f64,
    bicarbonate_dissociation_constant_mol_per_m3: f64,
    carbonate_dissociation_constant_mol2_per_m6: f64,
    minimum_exchange_soil_mass_megagrams: f64,
    concentration_floor_mol_per_m3: f64,
    carbon_molar_mass_g_per_mol: f64,
};

pub const RestrictedSaltResult = struct {
    aqueous: aqueous_network.State,
    cation_exchange_capacity_mol_charge_per_megagram: f64,
};

/// Direct source-order concentration initialization for SOLUTE.F 2914--2938.
/// Lines 2914--2923 duplicate the same CEC assignment; the returned scalar is
/// their common result. Existing nutrient-zone fields in `current` are kept.
pub fn initializeRestrictedSaltConcentrations(
    current: aqueous_network.State,
    inputs: RestrictedSaltInputs,
) !RestrictedSaltResult {
    inline for (@typeInfo(RestrictedSaltAmounts).@"struct".fields) |field| {
        const value = @field(inputs.amounts, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRestrictedSaltInput;
    }
    inline for (@typeInfo(RestrictedSaltInputs).@"struct".fields) |field| {
        if (field.type == RestrictedSaltAmounts) continue;
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRestrictedSaltInput;
    }
    if (inputs.water_volume_m3 <= 0 or
        inputs.concentration_floor_mol_per_m3 <= 0 or
        inputs.carbon_molar_mass_g_per_mol <= 0)
        return error.InvalidRestrictedSaltInput;

    const cec = if (inputs.exchange_soil_mass_megagrams >
        inputs.minimum_exchange_soil_mass_megagrams)
        @max(
            inputs.concentration_floor_mol_per_m3,
            inputs.cation_exchange_capacity_mol_charge /
                inputs.exchange_soil_mass_megagrams,
        )
    else
        0;
    const a = inputs.amounts;
    const v = inputs.water_volume_m3;
    var next = current;
    inline for (.{
        .{ "hydrogen", "hydrogen_mol" },
        .{ "hydroxide", "hydroxide_mol" },
        .{ "aluminum", "aluminum_mol" },
        .{ "iron", "iron_mol" },
        .{ "calcium", "calcium_mol" },
        .{ "magnesium", "magnesium_mol" },
        .{ "sodium", "sodium_mol" },
        .{ "potassium", "potassium_mol" },
        .{ "sulfate", "sulfate_mol" },
        .{ "chloride", "chloride_mol" },
    }) |names|
        @field(next, names[0]) = @max(
            inputs.concentration_floor_mol_per_m3,
            @field(a, names[1]) / v,
        );
    next.carbon_dioxide = @max(
        inputs.concentration_floor_mol_per_m3,
        inputs.carbon_dioxide_g_c_per_m3 / inputs.carbon_molar_mass_g_per_mol,
    );
    next.bicarbonate = @max(
        inputs.concentration_floor_mol_per_m3,
        next.carbon_dioxide * inputs.bicarbonate_dissociation_constant_mol_per_m3 /
            next.hydrogen,
    );
    next.carbonate = @max(
        inputs.concentration_floor_mol_per_m3,
        next.carbon_dioxide * inputs.carbonate_dissociation_constant_mol2_per_m6 /
            std.math.pow(f64, next.hydrogen, 2),
    );
    next.aluminum_hydroxide_2 = @max(
        inputs.concentration_floor_mol_per_m3,
        a.aluminum_hydroxide_2_mol / v,
    );
    next.iron_hydroxide_2 = @max(
        inputs.concentration_floor_mol_per_m3,
        a.iron_hydroxide_2_mol / v,
    );
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)) or
            @field(next, field.name) < 0)
            return error.InvalidRestrictedSaltConcentration;
    return .{
        .aqueous = next,
        .cation_exchange_capacity_mol_charge_per_megagram = cec,
    };
}

/// Direct translation of SOLUTE lines 646--719. Nutrient-zone ammonium,
/// ammonia, and nitrate fields are separate owners and remain unchanged.
pub fn applyDynamicSaltAqueousConcentrations(
    state: *aqueous_network.State,
    inputs: DynamicSaltAqueousInputs,
) !void {
    try validateDynamicSaltAqueousInputs(inputs);
    var staged = state.*;
    if (inputs.water_volume_m3 > inputs.minimum_water_volume_m3) {
        const volume = inputs.water_volume_m3;
        const free = inputs.free_ions;
        staged.hydrogen = amountPerWaterVolume(free.hydrogen_mol, volume);
        staged.hydroxide = amountPerWaterVolume(free.hydroxide_mol, volume);
        staged.aluminum = amountPerWaterVolume(free.aluminum_mol, volume);
        staged.iron = amountPerWaterVolume(free.iron_mol, volume);
        staged.calcium = amountPerWaterVolume(free.calcium_mol, volume);
        staged.magnesium = amountPerWaterVolume(free.magnesium_mol, volume);
        staged.sodium = amountPerWaterVolume(free.sodium_mol, volume);
        staged.potassium = amountPerWaterVolume(free.potassium_mol, volume);
        staged.sulfate = amountPerWaterVolume(free.sulfate_mol, volume);
        staged.chloride = amountPerWaterVolume(free.chloride_mol, volume);
        staged.carbonate = amountPerWaterVolume(free.carbonate_mol, volume);
        staged.bicarbonate = amountPerWaterVolume(free.bicarbonate_mol, volume);
        // SOLUTE.F line 658: preserve CO2S / (12.0 * VOLW), including the
        // source multiplication-before-division floating-point order.
        staged.carbon_dioxide = @max(
            minimum_concentration_mol_per_m3,
            free.carbon_dioxide_g_c /
                (inputs.carbon_molar_mass_g_per_mol * volume),
        );

        const complexes = inputs.dissolved_complexes;
        staged.aluminum_hydroxide_1 = amountPerWaterVolume(
            complexes.aluminum_hydroxide_1_mol,
            volume,
        );
        staged.aluminum_hydroxide_2 = amountPerWaterVolume(
            complexes.aluminum_hydroxide_2_mol,
            volume,
        );
        staged.aluminum_hydroxide_3 = amountPerWaterVolume(
            complexes.aluminum_hydroxide_3_mol,
            volume,
        );
        staged.aluminum_hydroxide_4 = amountPerWaterVolume(
            complexes.aluminum_hydroxide_4_mol,
            volume,
        );
        staged.aluminum_sulfate = amountPerWaterVolume(
            complexes.aluminum_sulfate_mol,
            volume,
        );
        staged.iron_hydroxide_1 = amountPerWaterVolume(
            complexes.iron_hydroxide_1_mol,
            volume,
        );
        staged.iron_hydroxide_2 = amountPerWaterVolume(
            complexes.iron_hydroxide_2_mol,
            volume,
        );
        staged.iron_hydroxide_3 = amountPerWaterVolume(
            complexes.iron_hydroxide_3_mol,
            volume,
        );
        staged.iron_hydroxide_4 = amountPerWaterVolume(
            complexes.iron_hydroxide_4_mol,
            volume,
        );
        staged.iron_sulfate = amountPerWaterVolume(
            complexes.iron_sulfate_mol,
            volume,
        );
        staged.calcium_hydroxide = amountPerWaterVolume(
            complexes.calcium_hydroxide_mol,
            volume,
        );
        staged.calcium_carbonate = amountPerWaterVolume(
            complexes.calcium_carbonate_mol,
            volume,
        );
        staged.calcium_bicarbonate = amountPerWaterVolume(
            complexes.calcium_bicarbonate_mol,
            volume,
        );
        staged.calcium_sulfate = amountPerWaterVolume(
            complexes.calcium_sulfate_mol,
            volume,
        );
        staged.magnesium_hydroxide = amountPerWaterVolume(
            complexes.magnesium_hydroxide_mol,
            volume,
        );
        staged.magnesium_carbonate = amountPerWaterVolume(
            complexes.magnesium_carbonate_mol,
            volume,
        );
        staged.magnesium_bicarbonate = amountPerWaterVolume(
            complexes.magnesium_bicarbonate_mol,
            volume,
        );
        staged.magnesium_sulfate = amountPerWaterVolume(
            complexes.magnesium_sulfate_mol,
            volume,
        );
        staged.sodium_carbonate = amountPerWaterVolume(
            complexes.sodium_carbonate_mol,
            volume,
        );
        staged.sodium_sulfate = amountPerWaterVolume(
            complexes.sodium_sulfate_mol,
            volume,
        );
        staged.potassium_sulfate = amountPerWaterVolume(
            complexes.potassium_sulfate_mol,
            volume,
        );
        staged.hydrogen_silicate = amountPerWaterVolume(
            complexes.hydrogen_silicate_mol,
            volume,
        );
    } else {
        staged.hydrogen =
            std.math.pow(f64, 10, -inputs.dry_soil_ph) * 1.0e3;
        staged.hydroxide =
            inputs.water_activity_product_mol2_per_m6 / staged.hydrogen;
        inline for (.{
            "aluminum",
            "iron",
            "calcium",
            "magnesium",
            "sodium",
            "potassium",
            "sulfate",
            "chloride",
            "carbonate",
            "bicarbonate",
            "carbon_dioxide",
            "aluminum_hydroxide_1",
            "aluminum_hydroxide_2",
            "aluminum_hydroxide_3",
            "aluminum_hydroxide_4",
            "aluminum_sulfate",
            "iron_hydroxide_1",
            "iron_hydroxide_2",
            "iron_hydroxide_3",
            "iron_hydroxide_4",
            "iron_sulfate",
            "calcium_hydroxide",
            "calcium_carbonate",
            "calcium_bicarbonate",
            "calcium_sulfate",
            "magnesium_hydroxide",
            "magnesium_carbonate",
            "magnesium_bicarbonate",
            "magnesium_sulfate",
            "sodium_carbonate",
            "sodium_sulfate",
            "potassium_sulfate",
            "hydrogen_silicate",
        }) |field_name| @field(staged, field_name) = 0;
    }
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field| {
        const value = @field(staged, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDynamicSaltAqueousConcentration;
    }
    state.* = staged;
}

fn validateDynamicSaltAqueousInputs(inputs: DynamicSaltAqueousInputs) !void {
    inline for (@typeInfo(FreeIonAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.free_ions, field.name)))
            return error.InvalidDynamicSaltAqueousInput;
    }
    inline for (@typeInfo(DissolvedComplexAmounts).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.dissolved_complexes, field.name)))
            return error.InvalidDynamicSaltAqueousInput;
    }
    if (!std.math.isFinite(inputs.water_volume_m3) or
        inputs.water_volume_m3 < 0 or
        !std.math.isFinite(inputs.minimum_water_volume_m3) or
        inputs.minimum_water_volume_m3 < 0 or
        !std.math.isFinite(inputs.dry_soil_ph) or
        !std.math.isFinite(inputs.water_activity_product_mol2_per_m6) or
        inputs.water_activity_product_mol2_per_m6 <= 0 or
        !std.math.isFinite(inputs.carbon_molar_mass_g_per_mol) or
        inputs.carbon_molar_mass_g_per_mol <= 0)
    {
        return error.InvalidDynamicSaltAqueousInput;
    }
}

fn amountPerWaterVolume(amount_mol: f64, water_volume_m3: f64) f64 {
    return @max(
        minimum_concentration_mol_per_m3,
        amount_mol / water_volume_m3,
    );
}

fn dynamicSaltInputs() DynamicSaltAqueousInputs {
    return .{
        .free_ions = .{
            .hydrogen_mol = 2,
            .hydroxide_mol = 4,
            .aluminum_mol = 6,
            .iron_mol = 8,
            .calcium_mol = 10,
            .magnesium_mol = 12,
            .sodium_mol = 14,
            .potassium_mol = 16,
            .sulfate_mol = 18,
            .chloride_mol = 20,
            .carbonate_mol = 22,
            .bicarbonate_mol = 24,
            .carbon_dioxide_g_c = 48,
        },
        .dissolved_complexes = .{
            .aluminum_hydroxide_1_mol = 2,
            .aluminum_hydroxide_2_mol = 4,
            .aluminum_hydroxide_3_mol = 6,
            .aluminum_hydroxide_4_mol = 8,
            .aluminum_sulfate_mol = 10,
            .iron_hydroxide_1_mol = 12,
            .iron_hydroxide_2_mol = 14,
            .iron_hydroxide_3_mol = 16,
            .iron_hydroxide_4_mol = 18,
            .iron_sulfate_mol = 20,
            .calcium_hydroxide_mol = 22,
            .calcium_carbonate_mol = 24,
            .calcium_bicarbonate_mol = 26,
            .calcium_sulfate_mol = 28,
            .magnesium_hydroxide_mol = 30,
            .magnesium_carbonate_mol = 32,
            .magnesium_bicarbonate_mol = 34,
            .magnesium_sulfate_mol = 36,
            .sodium_carbonate_mol = 38,
            .sodium_sulfate_mol = 40,
            .potassium_sulfate_mol = 42,
            .hydrogen_silicate_mol = 44,
        },
        .water_volume_m3 = 2,
        .minimum_water_volume_m3 = 1.0e-12,
        .dry_soil_ph = 7,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
        .carbon_molar_mass_g_per_mol = 12,
    };
}

test "SOLUTE nitrate concentrations preserve source conversion and dry branch" {
    const result = try initializeNitrateConcentrations(
        .{ .nitrate_g_n = 28, .water_volume_m3 = 2 },
        .{ .nitrate_g_n = -7, .water_volume_m3 = 1 },
        0,
        14,
        1.0e-32,
    );
    try std.testing.expectEqual(@as(f64, 1), result.non_band_mol_n_per_m3);
    try std.testing.expectEqual(
        minimum_concentration_mol_per_m3,
        result.band_mol_n_per_m3,
    );
    const dry = try initializeNitrateConcentrations(
        .{ .nitrate_g_n = 28, .water_volume_m3 = 0 },
        .{ .nitrate_g_n = 28, .water_volume_m3 = 1.0e-9 },
        1.0e-9,
        14,
        1.0e-32,
    );
    try std.testing.expectEqual(@as(f64, 0), dry.non_band_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 0), dry.band_mol_n_per_m3);
}

test "SOLUTE nitrate conversion uses compulsory runtime constants" {
    const result = try initializeNitrateConcentrations(
        .{ .nitrate_g_n = 30, .water_volume_m3 = 2 },
        .{ .nitrate_g_n = -1, .water_volume_m3 = 1 },
        0,
        15,
        2.0e-30,
    );
    try std.testing.expectEqual(@as(f64, 1), result.non_band_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 2.0e-30), result.band_mol_n_per_m3);
    try std.testing.expectError(
        error.InvalidNitrateConcentrationInput,
        initializeNitrateConcentrations(
            .{ .nitrate_g_n = 1, .water_volume_m3 = 1 },
            .{ .nitrate_g_n = 1, .water_volume_m3 = 1 },
            0,
            0,
            1.0e-32,
        ),
    );
}

test "SOLUTE cation exchange setup preserves source order and units" {
    const result = try initializeCationExchangeConcentrations(
        .{
            .capacity_mol_charge = 20,
            .hydrogen_mol = 2,
            .aluminum_mol = 4,
            .iron_mol = 6,
            .calcium_mol = 8,
            .magnesium_mol = 10,
            .sodium_mol = 12,
            .potassium_mol = 14,
            .carboxyl_bound_hydrogen_mol = 16,
        },
        2,
        0.5,
        5,
        4.0e8,
        0,
    );
    try std.testing.expectEqual(@as(f64, 10), result.capacity_mol_charge_per_megagram);
    try std.testing.expectEqual(@as(f64, 1), result.hydrogen_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 8), result.carboxyl_bound_hydrogen_mol_per_megagram);
    try std.testing.expectApproxEqAbs(
        @as(f64, 1000),
        result.total_carboxyl_sites_mol_per_megagram,
        1.0e-12,
    );
    try std.testing.expectEqual(@as(f64, 4), result.soil_mass_per_water_volume_megagrams_per_m3);
}

test "SOLUTE cation exchange dry branch returns source defaults" {
    const zero = std.mem.zeroes(CationExchangeAmounts);
    const result = try initializeCationExchangeConcentrations(
        zero,
        0,
        0,
        0,
        0,
        1.0e-12,
    );
    try std.testing.expectEqual(@as(f64, 0), result.capacity_mol_charge_per_megagram);
    try std.testing.expectEqual(@as(f64, 1), result.soil_mass_per_water_volume_megagrams_per_m3);
}

test "SOLUTE wet dynamic salt initialization preserves units and nutrient owners" {
    var state = std.mem.zeroes(aqueous_network.State);
    state.ammonium_non_band = 101;
    state.ammonia_non_band = 102;
    state.ammonium_band = 103;
    state.ammonia_band = 104;
    state.nitrate_non_band = 105;
    state.nitrate_band = 106;
    var inputs = dynamicSaltInputs();
    inputs.free_ions.iron_mol = -1;
    try applyDynamicSaltAqueousConcentrations(&state, inputs);

    try std.testing.expectEqual(@as(f64, 1), state.hydrogen);
    try std.testing.expectEqual(
        minimum_concentration_mol_per_m3,
        state.iron,
    );
    try std.testing.expectEqual(@as(f64, 2), state.carbon_dioxide);
    try std.testing.expectEqual(@as(f64, 1), state.aluminum_hydroxide_1);
    try std.testing.expectEqual(@as(f64, 22), state.hydrogen_silicate);
    try std.testing.expectEqual(@as(f64, 101), state.ammonium_non_band);
    try std.testing.expectEqual(@as(f64, 106), state.nitrate_band);
}

test "SOLUTE line 658 preserves carbon dioxide denominator operation order" {
    var state = std.mem.zeroes(aqueous_network.State);
    var inputs = dynamicSaltInputs();
    inputs.water_volume_m3 = 10;
    inputs.free_ions.carbon_dioxide_g_c = 1;
    inputs.carbon_molar_mass_g_per_mol = 13;
    try applyDynamicSaltAqueousConcentrations(&state, inputs);

    const source_order = 1.0 / (13.0 * 10.0);
    try std.testing.expectEqual(source_order, state.carbon_dioxide);
}

test "SOLUTE line 658 rejects invalid runtime carbon molar mass" {
    var state = std.mem.zeroes(aqueous_network.State);
    var inputs = dynamicSaltInputs();
    inputs.carbon_molar_mass_g_per_mol = 0;
    try std.testing.expectError(
        error.InvalidDynamicSaltAqueousInput,
        applyDynamicSaltAqueousConcentrations(&state, inputs),
    );
}

test "SOLUTE dry dynamic salt initialization uses pH and clears only salts" {
    var state: aqueous_network.State = undefined;
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field|
        @field(state, field.name) = 9;
    var inputs = dynamicSaltInputs();
    inputs.water_volume_m3 = inputs.minimum_water_volume_m3;
    try applyDynamicSaltAqueousConcentrations(&state, inputs);

    try std.testing.expectApproxEqAbs(@as(f64, 1.0e-4), state.hydrogen, 1.0e-18);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0e-4), state.hydroxide, 1.0e-18);
    try std.testing.expectEqual(@as(f64, 0), state.aluminum);
    try std.testing.expectEqual(@as(f64, 0), state.hydrogen_silicate);
    try std.testing.expectEqual(@as(f64, 9), state.ammonium_band);
    try std.testing.expectEqual(@as(f64, 9), state.nitrate_non_band);
}

test "SOLUTE dynamic salt initialization rejects late invalid input atomically" {
    var state: aqueous_network.State = undefined;
    inline for (@typeInfo(aqueous_network.State).@"struct".fields) |field|
        @field(state, field.name) = 3;
    const before = state;
    var inputs = dynamicSaltInputs();
    inputs.dissolved_complexes.potassium_sulfate_mol =
        std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidDynamicSaltAqueousInput,
        applyDynamicSaltAqueousConcentrations(&state, inputs),
    );
    try std.testing.expectEqualDeep(before, state);
}

test "SOLUTE restricted salt concentrations preserve source subset and order" {
    var current = std.mem.zeroes(aqueous_network.State);
    current.ammonium_non_band = 9;
    current.nitrate_band = 8;
    const inputs = RestrictedSaltInputs{
        .amounts = .{
            .hydrogen_mol = 2,
            .hydroxide_mol = 4,
            .aluminum_mol = 6,
            .iron_mol = 8,
            .calcium_mol = 10,
            .magnesium_mol = 12,
            .sodium_mol = 14,
            .potassium_mol = 16,
            .sulfate_mol = 18,
            .chloride_mol = 20,
            .aluminum_hydroxide_2_mol = 6,
            .iron_hydroxide_2_mol = 8,
        },
        .water_volume_m3 = 2,
        .exchange_soil_mass_megagrams = 2,
        .cation_exchange_capacity_mol_charge = 10,
        .carbon_dioxide_g_c_per_m3 = 24,
        .bicarbonate_dissociation_constant_mol_per_m3 = 0.5,
        .carbonate_dissociation_constant_mol2_per_m6 = 0.25,
        .minimum_exchange_soil_mass_megagrams = 1.0e-12,
        .concentration_floor_mol_per_m3 = 1.0e-32,
        .carbon_molar_mass_g_per_mol = 12,
    };
    const result = try initializeRestrictedSaltConcentrations(current, inputs);
    try std.testing.expectEqual(@as(f64, 5), result.cation_exchange_capacity_mol_charge_per_megagram);
    try std.testing.expectEqual(@as(f64, 1), result.aqueous.hydrogen);
    try std.testing.expectEqual(@as(f64, 2), result.aqueous.carbon_dioxide);
    try std.testing.expectEqual(@as(f64, 1), result.aqueous.bicarbonate);
    try std.testing.expectEqual(@as(f64, 0.5), result.aqueous.carbonate);
    try std.testing.expectEqual(@as(f64, 3), result.aqueous.aluminum_hydroxide_2);
    try std.testing.expectEqual(@as(f64, 4), result.aqueous.iron_hydroxide_2);
    try std.testing.expectEqual(@as(f64, 9), result.aqueous.ammonium_non_band);
    try std.testing.expectEqual(@as(f64, 8), result.aqueous.nitrate_band);

    var threshold = inputs;
    threshold.exchange_soil_mass_megagrams = threshold.minimum_exchange_soil_mass_megagrams;
    const inactive = try initializeRestrictedSaltConcentrations(current, threshold);
    try std.testing.expectEqual(
        @as(f64, 0),
        inactive.cation_exchange_capacity_mol_charge_per_megagram,
    );
}
