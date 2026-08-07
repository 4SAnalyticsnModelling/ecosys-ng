const std = @import("std");
const activity_coefficients = @import("activity_coefficients.zig");

/// Extensive aqueous inventories read by the restricted-salt branch.
/// Every field except CO2 is mol per runtime soil layer and grid cell.
pub const ExtensiveAqueousState = struct {
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
    carbon_dioxide_g_c_per_m3: f64,
};

/// Concentrations already reconstructed for one fertilizer placement zone.
pub const ZoneConcentrations = struct {
    hpo4_mol_p_per_m3: f64,
    h2po4_mol_p_per_m3: f64,
    ammonium_mol_n_per_m3: f64,
    ammonia_mol_n_per_m3: f64,
};

pub const CarbonateEquilibrium = struct {
    bicarbonate_dissociation_constant_mol_per_m3: f64,
    carbonate_dissociation_constant_mol2_per_m6: f64,
};

pub const Inputs = struct {
    extensive: ExtensiveAqueousState,
    non_band: ZoneConcentrations,
    band: ZoneConcentrations,
    matrix_water_volume_m3: f64,
    soil_mass_megagrams: f64,
    minimum_active_soil_mass_megagrams: f64,
    cation_exchange_capacity_mol_charge: f64,
    /// Runtime replacement for source `ZEROC`. The source reuses this numeric
    /// floor for both mol m-3 and mol Mg-1 values.
    minimum_positive_value: f64,
    carbonate_equilibrium: CarbonateEquilibrium,
    activity_coefficients: activity_coefficients.Result,
};

pub const AqueousConcentrations = struct {
    hydrogen_mol_per_m3: f64,
    hydroxide_mol_per_m3: f64,
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
    sodium_mol_per_m3: f64,
    potassium_mol_per_m3: f64,
    sulfate_mol_per_m3: f64,
    chloride_mol_per_m3: f64,
    carbon_dioxide_mol_per_m3: f64,
    bicarbonate_mol_per_m3: f64,
    carbonate_mol_per_m3: f64,
    aluminum_hydroxide_2_mol_per_m3: f64,
    iron_hydroxide_2_mol_per_m3: f64,
};

pub const AqueousActivities = struct {
    hydrogen_mol_per_m3: f64,
    hydroxide_mol_per_m3: f64,
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
    sodium_mol_per_m3: f64,
    potassium_mol_per_m3: f64,
    sulfate_mol_per_m3: f64,
    carbon_dioxide_mol_per_m3: f64,
    bicarbonate_mol_per_m3: f64,
    carbonate_mol_per_m3: f64,
    aluminum_hydroxide_2_mol_per_m3: f64,
    iron_hydroxide_2_mol_per_m3: f64,
};

pub const ZoneActivities = struct {
    hpo4_mol_p_per_m3: f64,
    h2po4_mol_p_per_m3: f64,
    ammonium_mol_n_per_m3: f64,
    ammonia_mol_n_per_m3: f64,
};

pub const Result = struct {
    cation_exchange_capacity_mol_charge_per_megagram: f64,
    concentrations: AqueousConcentrations,
    activities: AqueousActivities,
    non_band_activities: ZoneActivities,
    band_activities: ZoneActivities,
};

/// Direct source-order translation of SOLUTE.F lines 2914--2960.
///
/// This pure cell/layer kernel neither owns grid dimensions nor mutates
/// chemistry state. Callers allocate runtime-sized state and publish the
/// returned values only after this function succeeds.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validate(inputs);
    const floor = inputs.minimum_positive_value;
    const water_volume_m3 = inputs.matrix_water_volume_m3;

    // SOLUTE.F repeats this complete CCEC branch twice at lines 2914--2923.
    // The duplicate is a numerical no-op and is collapsed without reordering.
    const cation_exchange_capacity = if (inputs.soil_mass_megagrams >
        inputs.minimum_active_soil_mass_megagrams)
        @max(
            floor,
            inputs.cation_exchange_capacity_mol_charge /
                inputs.soil_mass_megagrams,
        )
    else
        0;

    const extensive = inputs.extensive;
    const hydrogen = @max(floor, extensive.hydrogen_mol / water_volume_m3);
    const hydroxide = @max(floor, extensive.hydroxide_mol / water_volume_m3);
    const aluminum = @max(floor, extensive.aluminum_mol / water_volume_m3);
    const iron = @max(floor, extensive.iron_mol / water_volume_m3);
    const calcium = @max(floor, extensive.calcium_mol / water_volume_m3);
    const magnesium = @max(floor, extensive.magnesium_mol / water_volume_m3);
    const sodium = @max(floor, extensive.sodium_mol / water_volume_m3);
    const potassium = @max(floor, extensive.potassium_mol / water_volume_m3);
    const sulfate = @max(floor, extensive.sulfate_mol / water_volume_m3);
    const chloride = @max(floor, extensive.chloride_mol / water_volume_m3);

    // HOUR1.F has already converted CCO2S to g C m-3, so SOLUTE.F divides
    // only by the carbon molar mass here.
    const carbon_dioxide = @max(
        floor,
        extensive.carbon_dioxide_g_c_per_m3 / 12.0,
    );
    const equilibrium = inputs.carbonate_equilibrium;
    const bicarbonate = @max(
        floor,
        carbon_dioxide *
            equilibrium.bicarbonate_dissociation_constant_mol_per_m3 /
            hydrogen,
    );
    const carbonate = @max(
        floor,
        carbon_dioxide *
            equilibrium.carbonate_dissociation_constant_mol2_per_m6 /
            (hydrogen * hydrogen),
    );
    const aluminum_hydroxide_2 = @max(
        floor,
        extensive.aluminum_hydroxide_2_mol / water_volume_m3,
    );
    const iron_hydroxide_2 = @max(
        floor,
        extensive.iron_hydroxide_2_mol / water_volume_m3,
    );
    const concentrations = AqueousConcentrations{
        .hydrogen_mol_per_m3 = hydrogen,
        .hydroxide_mol_per_m3 = hydroxide,
        .aluminum_mol_per_m3 = aluminum,
        .iron_mol_per_m3 = iron,
        .calcium_mol_per_m3 = calcium,
        .magnesium_mol_per_m3 = magnesium,
        .sodium_mol_per_m3 = sodium,
        .potassium_mol_per_m3 = potassium,
        .sulfate_mol_per_m3 = sulfate,
        .chloride_mol_per_m3 = chloride,
        .carbon_dioxide_mol_per_m3 = carbon_dioxide,
        .bicarbonate_mol_per_m3 = bicarbonate,
        .carbonate_mol_per_m3 = carbonate,
        .aluminum_hydroxide_2_mol_per_m3 = aluminum_hydroxide_2,
        .iron_hydroxide_2_mol_per_m3 = iron_hydroxide_2,
    };

    const coefficients = inputs.activity_coefficients;
    const monovalent = coefficients.monovalent_activity_coefficient;
    const divalent = coefficients.divalent_activity_coefficient;
    const trivalent = coefficients.trivalent_activity_coefficient;
    const result = Result{
        .cation_exchange_capacity_mol_charge_per_megagram = cation_exchange_capacity,
        .concentrations = concentrations,
        .activities = .{
            .hydrogen_mol_per_m3 = hydrogen * monovalent,
            .hydroxide_mol_per_m3 = hydroxide * monovalent,
            .aluminum_mol_per_m3 = aluminum * trivalent,
            .iron_mol_per_m3 = iron * trivalent,
            .calcium_mol_per_m3 = calcium * divalent,
            .magnesium_mol_per_m3 = magnesium * divalent,
            .sodium_mol_per_m3 = sodium * monovalent,
            .potassium_mol_per_m3 = potassium * monovalent,
            .sulfate_mol_per_m3 = sulfate * divalent,
            .carbon_dioxide_mol_per_m3 = carbon_dioxide,
            .bicarbonate_mol_per_m3 = bicarbonate * monovalent,
            .carbonate_mol_per_m3 = carbonate * divalent,
            .aluminum_hydroxide_2_mol_per_m3 = aluminum_hydroxide_2 * monovalent,
            .iron_hydroxide_2_mol_per_m3 = iron_hydroxide_2 * monovalent,
        },
        .non_band_activities = zoneActivities(
            inputs.non_band,
            monovalent,
            divalent,
        ),
        .band_activities = zoneActivities(
            inputs.band,
            monovalent,
            divalent,
        ),
    };
    try validateResult(result);
    return result;
}

fn zoneActivities(
    concentrations: ZoneConcentrations,
    monovalent: f64,
    divalent: f64,
) ZoneActivities {
    return .{
        .hpo4_mol_p_per_m3 = concentrations.hpo4_mol_p_per_m3 * divalent,
        .h2po4_mol_p_per_m3 = concentrations.h2po4_mol_p_per_m3 * monovalent,
        .ammonium_mol_n_per_m3 = concentrations.ammonium_mol_n_per_m3 *
            monovalent,
        .ammonia_mol_n_per_m3 = concentrations.ammonia_mol_n_per_m3,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(ExtensiveAqueousState).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.extensive, field.name)))
            return error.InvalidFixedPhInitialization;
    inline for (@typeInfo(ZoneConcentrations).@"struct".fields) |field| {
        const non_band = @field(inputs.non_band, field.name);
        const band = @field(inputs.band, field.name);
        if (!std.math.isFinite(non_band) or non_band < 0 or
            !std.math.isFinite(band) or band < 0)
            return error.InvalidFixedPhInitialization;
    }
    inline for (.{
        inputs.matrix_water_volume_m3,
        inputs.minimum_positive_value,
    }) |positive| {
        if (!std.math.isFinite(positive) or positive <= 0)
            return error.InvalidFixedPhInitialization;
    }
    inline for (.{
        inputs.soil_mass_megagrams,
        inputs.minimum_active_soil_mass_megagrams,
    }) |nonnegative| {
        if (!std.math.isFinite(nonnegative) or nonnegative < 0)
            return error.InvalidFixedPhInitialization;
    }
    if (!std.math.isFinite(inputs.cation_exchange_capacity_mol_charge))
        return error.InvalidFixedPhInitialization;
    inline for (@typeInfo(CarbonateEquilibrium).@"struct".fields) |field| {
        const value = @field(inputs.carbonate_equilibrium, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhInitialization;
    }
    inline for (.{
        inputs.activity_coefficients.monovalent_activity_coefficient,
        inputs.activity_coefficients.divalent_activity_coefficient,
        inputs.activity_coefficients.trivalent_activity_coefficient,
    }) |coefficient| {
        if (!std.math.isFinite(coefficient) or coefficient <= 0)
            return error.InvalidFixedPhInitialization;
    }
}

fn validateResult(result: Result) !void {
    if (!std.math.isFinite(result.cation_exchange_capacity_mol_charge_per_megagram))
        return error.NonFiniteFixedPhInitialization;
    inline for (@typeInfo(AqueousConcentrations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.concentrations, field.name)))
            return error.NonFiniteFixedPhInitialization;
    inline for (@typeInfo(AqueousActivities).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.activities, field.name)))
            return error.NonFiniteFixedPhInitialization;
    inline for (@typeInfo(ZoneActivities).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result.non_band_activities, field.name)) or
            !std.math.isFinite(@field(result.band_activities, field.name)))
            return error.NonFiniteFixedPhInitialization;
    }
}

fn validInputs() Inputs {
    return .{
        .extensive = .{
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
            .aluminum_hydroxide_2_mol = 22,
            .iron_hydroxide_2_mol = 24,
            .carbon_dioxide_g_c_per_m3 = 48,
        },
        .non_band = .{
            .hpo4_mol_p_per_m3 = 8,
            .h2po4_mol_p_per_m3 = 10,
            .ammonium_mol_n_per_m3 = 16,
            .ammonia_mol_n_per_m3 = 20,
        },
        .band = .{
            .hpo4_mol_p_per_m3 = 12,
            .h2po4_mol_p_per_m3 = 14,
            .ammonium_mol_n_per_m3 = 18,
            .ammonia_mol_n_per_m3 = 22,
        },
        .matrix_water_volume_m3 = 2,
        .soil_mass_megagrams = 4,
        .minimum_active_soil_mass_megagrams = 1,
        .cation_exchange_capacity_mol_charge = 8,
        .minimum_positive_value = 0.01,
        .carbonate_equilibrium = .{
            .bicarbonate_dissociation_constant_mol_per_m3 = 0.5,
            .carbonate_dissociation_constant_mol2_per_m6 = 0.25,
        },
        .activity_coefficients = .{
            .ionic_strength_mol_per_l = 0,
            .monovalent_activity_coefficient = 0.5,
            .divalent_activity_coefficient = 0.25,
            .trivalent_activity_coefficient = 0.125,
            .total_ion_activity_mol_per_m3 = 0,
            .electrical_conductivity_dS_per_m = 0,
        },
    };
}

test "fixed-pH initialization preserves source concentration and activity order" {
    const result = try calculateSourceOrder(validInputs());
    try std.testing.expectEqual(
        @as(f64, 2),
        result.cation_exchange_capacity_mol_charge_per_megagram,
    );
    try std.testing.expectEqual(@as(f64, 1), result.concentrations.hydrogen_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 10), result.concentrations.chloride_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4), result.concentrations.carbon_dioxide_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 2), result.concentrations.bicarbonate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), result.concentrations.carbonate_mol_per_m3);
    try std.testing.expectEqual(
        @as(f64, 11),
        result.concentrations.aluminum_hydroxide_2_mol_per_m3,
    );
    try std.testing.expectEqual(@as(f64, 0.5), result.activities.hydrogen_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 0.375), result.activities.aluminum_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1.25), result.activities.calcium_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4), result.activities.carbon_dioxide_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 5.5), result.activities.aluminum_hydroxide_2_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 2), result.non_band_activities.hpo4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 5), result.non_band_activities.h2po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 8), result.non_band_activities.ammonium_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 20), result.non_band_activities.ammonia_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 3), result.band_activities.hpo4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 7), result.band_activities.h2po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 9), result.band_activities.ammonium_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 22), result.band_activities.ammonia_mol_n_per_m3);
}

test "fixed-pH initialization floors signed source inventories" {
    var inputs = validInputs();
    inputs.extensive.hydrogen_mol = -2;
    inputs.extensive.carbon_dioxide_g_c_per_m3 = -48;
    inputs.cation_exchange_capacity_mol_charge = -8;
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        inputs.minimum_positive_value,
        result.cation_exchange_capacity_mol_charge_per_megagram,
    );
    try std.testing.expectEqual(
        inputs.minimum_positive_value,
        result.concentrations.hydrogen_mol_per_m3,
    );
    try std.testing.expectEqual(
        inputs.minimum_positive_value,
        result.concentrations.carbon_dioxide_mol_per_m3,
    );
    try std.testing.expect(std.math.isFinite(
        result.concentrations.carbonate_mol_per_m3,
    ));
}

test "fixed-pH initialization preserves inactive-soil CEC zero branch" {
    var inputs = validInputs();
    inputs.soil_mass_megagrams = inputs.minimum_active_soil_mass_megagrams;
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.cation_exchange_capacity_mol_charge_per_megagram,
    );
}

test "fixed-pH initialization rejects invalid runtime state" {
    var inputs = validInputs();
    inputs.matrix_water_volume_m3 = 0;
    try std.testing.expectError(
        error.InvalidFixedPhInitialization,
        calculateSourceOrder(inputs),
    );
    inputs = validInputs();
    inputs.band.ammonia_mol_n_per_m3 = -1;
    try std.testing.expectError(
        error.InvalidFixedPhInitialization,
        calculateSourceOrder(inputs),
    );
    inputs = validInputs();
    inputs.activity_coefficients.divalent_activity_coefficient =
        std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidFixedPhInitialization,
        calculateSourceOrder(inputs),
    );
}
