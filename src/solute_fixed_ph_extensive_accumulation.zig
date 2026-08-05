const std = @import("std");
const cation_exchange = @import("solute_cation_exchange.zig");
const transformations = @import("solute_fixed_ph_transformation_assembly.zig");
const water_reset = @import("solute_fixed_ph_water_reset.zig");

pub const Geometry = struct {
    shared_water_volume_m3: f64,
    non_band_ammonium_water_volume_m3: f64,
    band_ammonium_water_volume_m3: f64,
    non_band_phosphate_water_volume_m3: f64,
    band_phosphate_water_volume_m3: f64,
    shared_soil_mass_megagrams: f64,
    non_band_ammonium_soil_mass_megagrams: f64,
    band_ammonium_soil_mass_megagrams: f64,
};

pub const WaterResetChanges = struct {
    hydrogen_mol_per_m3_step: f64,
    hydroxide_mol_per_m3_step: f64,
};

pub const AmmoniumTotals = struct {
    non_band_ammonium_mol_n_per_step: f64,
    band_ammonium_mol_n_per_step: f64,
    non_band_ammonia_mol_n_per_step: f64,
    band_ammonia_mol_n_per_step: f64,
};

pub const SharedAqueousTotals = struct {
    hydrogen_mol_per_step: f64,
    hydroxide_mol_per_step: f64,
    aluminum_mol_per_step: f64,
    iron_mol_per_step: f64,
    calcium_mol_per_step: f64,
    magnesium_mol_per_step: f64,
    sodium_mol_per_step: f64,
    potassium_mol_per_step: f64,
};

pub const CationExchangeTotals = struct {
    non_band_ammonium_mol_n_per_step: f64,
    band_ammonium_mol_n_per_step: f64,
    hydrogen_mol_per_step: f64,
    aluminum_mol_per_step: f64,
    iron_mol_per_step: f64,
    calcium_mol_per_step: f64,
    magnesium_mol_per_step: f64,
    sodium_mol_per_step: f64,
    potassium_mol_per_step: f64,
};

pub const SurfacePhosphateTotals = struct {
    hydroxyl_site_mol_per_step: f64,
    protonated_site_mol_per_step: f64,
    adsorbed_hpo4_mol_p_per_step: f64,
    adsorbed_h2po4_mol_p_per_step: f64,
};

pub const PhosphateMineralTotals = struct {
    aluminum_phosphate_mol_mineral_per_step: f64,
    iron_phosphate_mol_mineral_per_step: f64,
    dicalcium_phosphate_mol_mineral_per_step: f64,
    hydroxyapatite_mol_mineral_per_step: f64,
    monocalcium_phosphate_mol_mineral_per_step: f64,
};

pub const PhosphateZoneTotals = struct {
    dissolved_hpo4_mol_p_per_step: f64,
    dissolved_h2po4_mol_p_per_step: f64,
    surface: SurfacePhosphateTotals,
    minerals: PhosphateMineralTotals,
};

pub const Totals = struct {
    ammonium: AmmoniumTotals,
    shared_aqueous: SharedAqueousTotals,
    water_equilibration_hydrogen_mol_per_step: f64,
    cation_exchange: CationExchangeTotals,
    non_band_phosphate: PhosphateZoneTotals,
    band_phosphate: PhosphateZoneTotals,
};

pub const Inputs = struct {
    previous_totals: Totals,
    geometry: Geometry,
    transformation_rates: transformations.Result,
    water_reset_changes: WaterResetChanges,
    cation_exchange_rates_mol_per_megagram_step: cation_exchange.Cations,
};

/// Exact source-order comparator for SOLUTE.F lines 3659--3725.
///
/// The source accumulates accepted concentration-space reaction changes into
/// extensive totals for REDIST. Aqueous, phosphate-surface, and mineral
/// changes use their owning water volume; Gapon changes use their owning soil
/// mass. This pure function stages all additions before returning.
pub fn accumulateSourceOrder(inputs: Inputs) !Totals {
    try validateInputs(inputs);
    const geometry = inputs.geometry;
    const rates = inputs.transformation_rates;
    const exchange = inputs.cation_exchange_rates_mol_per_megagram_step;
    var result = inputs.previous_totals;

    // SOLUTE.F 3682--3725. Preserve source statement and operation order.
    result.ammonium.non_band_ammonium_mol_n_per_step +=
        rates.ammonium.non_band_ammonium_mol_n_per_m3_step *
        geometry.non_band_ammonium_water_volume_m3;
    result.ammonium.band_ammonium_mol_n_per_step +=
        rates.ammonium.band_ammonium_mol_n_per_m3_step *
        geometry.band_ammonium_water_volume_m3;
    result.ammonium.non_band_ammonia_mol_n_per_step +=
        rates.ammonium.non_band_ammonia_mol_n_per_m3_step *
        geometry.non_band_ammonium_water_volume_m3;
    result.ammonium.band_ammonia_mol_n_per_step +=
        rates.ammonium.band_ammonia_mol_n_per_m3_step *
        geometry.band_ammonium_water_volume_m3;
    result.shared_aqueous.hydrogen_mol_per_step +=
        (rates.shared_aqueous.hydrogen_mol_per_m3_step +
            inputs.water_reset_changes.hydrogen_mol_per_m3_step) *
        geometry.shared_water_volume_m3;
    result.shared_aqueous.hydroxide_mol_per_step +=
        (rates.shared_aqueous.hydroxide_mol_per_m3_step +
            inputs.water_reset_changes.hydroxide_mol_per_m3_step) *
        geometry.shared_water_volume_m3;
    result.shared_aqueous.aluminum_mol_per_step +=
        rates.shared_aqueous.aluminum_mol_per_m3_step *
        geometry.shared_water_volume_m3;
    result.shared_aqueous.iron_mol_per_step +=
        rates.shared_aqueous.iron_mol_per_m3_step *
        geometry.shared_water_volume_m3;
    result.shared_aqueous.calcium_mol_per_step +=
        rates.shared_aqueous.calcium_mol_per_m3_step *
        geometry.shared_water_volume_m3;
    result.shared_aqueous.magnesium_mol_per_step +=
        rates.shared_aqueous.magnesium_mol_per_m3_step *
        geometry.shared_water_volume_m3;
    result.shared_aqueous.sodium_mol_per_step +=
        rates.shared_aqueous.sodium_mol_per_m3_step *
        geometry.shared_water_volume_m3;
    result.shared_aqueous.potassium_mol_per_step +=
        rates.shared_aqueous.potassium_mol_per_m3_step *
        geometry.shared_water_volume_m3;
    result.water_equilibration_hydrogen_mol_per_step +=
        inputs.water_reset_changes.hydrogen_mol_per_m3_step *
        geometry.shared_water_volume_m3;

    accumulatePhosphateDissolved(
        &result.non_band_phosphate,
        rates.non_band_phosphate,
        geometry.non_band_phosphate_water_volume_m3,
    );
    accumulatePhosphateDissolved(
        &result.band_phosphate,
        rates.band_phosphate,
        geometry.band_phosphate_water_volume_m3,
    );
    accumulateCationExchange(&result.cation_exchange, exchange, geometry);
    accumulatePhosphateSurface(
        &result.non_band_phosphate,
        rates.non_band_phosphate,
        geometry.non_band_phosphate_water_volume_m3,
    );
    accumulatePhosphateSurface(
        &result.band_phosphate,
        rates.band_phosphate,
        geometry.band_phosphate_water_volume_m3,
    );
    accumulatePhosphateMinerals(
        &result.non_band_phosphate,
        rates.non_band_phosphate,
        geometry.non_band_phosphate_water_volume_m3,
    );
    accumulatePhosphateMinerals(
        &result.band_phosphate,
        rates.band_phosphate,
        geometry.band_phosphate_water_volume_m3,
    );

    try validateTotals(result, error.NonFiniteFixedPhExtensiveResult);
    return result;
}

fn accumulatePhosphateDissolved(
    total: *PhosphateZoneTotals,
    rates: transformations.PhosphateZoneTransformations,
    water_volume_m3: f64,
) void {
    total.dissolved_hpo4_mol_p_per_step +=
        rates.dissolved_hpo4_mol_p_per_m3_step * water_volume_m3;
    total.dissolved_h2po4_mol_p_per_step +=
        rates.dissolved_h2po4_mol_p_per_m3_step * water_volume_m3;
}

fn accumulateCationExchange(
    total: *CationExchangeTotals,
    rates: cation_exchange.Cations,
    geometry: Geometry,
) void {
    total.non_band_ammonium_mol_n_per_step +=
        rates.ammonium_non_band * geometry.non_band_ammonium_soil_mass_megagrams;
    total.band_ammonium_mol_n_per_step +=
        rates.ammonium_band * geometry.band_ammonium_soil_mass_megagrams;
    total.hydrogen_mol_per_step += rates.hydrogen * geometry.shared_soil_mass_megagrams;
    total.aluminum_mol_per_step += rates.aluminum * geometry.shared_soil_mass_megagrams;
    total.iron_mol_per_step += rates.iron * geometry.shared_soil_mass_megagrams;
    total.calcium_mol_per_step += rates.calcium * geometry.shared_soil_mass_megagrams;
    total.magnesium_mol_per_step += rates.magnesium * geometry.shared_soil_mass_megagrams;
    total.sodium_mol_per_step += rates.sodium * geometry.shared_soil_mass_megagrams;
    total.potassium_mol_per_step += rates.potassium * geometry.shared_soil_mass_megagrams;
}

fn accumulatePhosphateSurface(
    total: *PhosphateZoneTotals,
    rates: transformations.PhosphateZoneTransformations,
    water_volume_m3: f64,
) void {
    total.surface.hydroxyl_site_mol_per_step +=
        rates.surface.hydroxyl_site_source_change_per_step * water_volume_m3;
    total.surface.protonated_site_mol_per_step +=
        rates.surface.protonated_site_source_change_per_step * water_volume_m3;
    total.surface.adsorbed_hpo4_mol_p_per_step +=
        rates.surface.adsorbed_hpo4_source_change_per_step * water_volume_m3;
    total.surface.adsorbed_h2po4_mol_p_per_step +=
        rates.surface.adsorbed_h2po4_source_change_per_step * water_volume_m3;
}

fn accumulatePhosphateMinerals(
    total: *PhosphateZoneTotals,
    rates: transformations.PhosphateZoneTransformations,
    water_volume_m3: f64,
) void {
    total.minerals.aluminum_phosphate_mol_mineral_per_step +=
        rates.solids.aluminum_phosphate_mol_per_m3_step * water_volume_m3;
    total.minerals.iron_phosphate_mol_mineral_per_step +=
        rates.solids.iron_phosphate_mol_per_m3_step * water_volume_m3;
    total.minerals.dicalcium_phosphate_mol_mineral_per_step +=
        rates.solids.dicalcium_phosphate_mol_per_m3_step * water_volume_m3;
    total.minerals.hydroxyapatite_mol_mineral_per_step +=
        rates.solids.hydroxyapatite_mol_per_m3_step * water_volume_m3;
    total.minerals.monocalcium_phosphate_mol_mineral_per_step +=
        rates.solids.monocalcium_phosphate_mol_per_m3_step * water_volume_m3;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Geometry).@"struct".fields) |field| {
        const value = @field(inputs.geometry, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhExtensiveGeometry;
    }
    try validateTotals(
        inputs.previous_totals,
        error.NonFiniteFixedPhExtensiveInput,
    );
    try validateFiniteStruct(
        inputs.water_reset_changes,
        error.NonFiniteFixedPhExtensiveInput,
    );
    try validateFiniteStruct(
        inputs.cation_exchange_rates_mol_per_megagram_step,
        error.NonFiniteFixedPhExtensiveInput,
    );
    try validateTransformationRates(inputs.transformation_rates);
}

fn validateTransformationRates(rates: transformations.Result) !void {
    try validateFiniteStruct(rates.ammonium, error.NonFiniteFixedPhExtensiveInput);
    try validateFiniteStruct(
        rates.shared_aqueous,
        error.NonFiniteFixedPhExtensiveInput,
    );
    try validatePhosphateZoneRates(rates.non_band_phosphate);
    try validatePhosphateZoneRates(rates.band_phosphate);
}

fn validatePhosphateZoneRates(
    rates: transformations.PhosphateZoneTransformations,
) !void {
    if (!std.math.isFinite(rates.dissolved_hpo4_mol_p_per_m3_step) or
        !std.math.isFinite(rates.dissolved_h2po4_mol_p_per_m3_step))
        return error.NonFiniteFixedPhExtensiveInput;
    try validateFiniteStruct(
        rates.surface,
        error.NonFiniteFixedPhExtensiveInput,
    );
    try validateFiniteStruct(
        rates.solids,
        error.NonFiniteFixedPhExtensiveInput,
    );
}

fn validateTotals(totals: Totals, failure: anyerror) !void {
    try validateFiniteStruct(totals.ammonium, failure);
    try validateFiniteStruct(totals.shared_aqueous, failure);
    if (!std.math.isFinite(totals.water_equilibration_hydrogen_mol_per_step))
        return failure;
    try validateFiniteStruct(totals.cation_exchange, failure);
    try validatePhosphateTotals(totals.non_band_phosphate, failure);
    try validatePhosphateTotals(totals.band_phosphate, failure);
}

fn validatePhosphateTotals(
    totals: PhosphateZoneTotals,
    failure: anyerror,
) !void {
    if (!std.math.isFinite(totals.dissolved_hpo4_mol_p_per_step) or
        !std.math.isFinite(totals.dissolved_h2po4_mol_p_per_step))
        return failure;
    try validateFiniteStruct(totals.surface, failure);
    try validateFiniteStruct(totals.minerals, failure);
}

fn validateFiniteStruct(values: anytype, failure: anyerror) !void {
    inline for (@typeInfo(@TypeOf(values)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(values, field.name))) return failure;
}

fn filledTotals(value: f64) Totals {
    return .{
        .ammonium = filled(AmmoniumTotals, value),
        .shared_aqueous = filled(SharedAqueousTotals, value),
        .water_equilibration_hydrogen_mol_per_step = value,
        .cation_exchange = filled(CationExchangeTotals, value),
        .non_band_phosphate = .{
            .dissolved_hpo4_mol_p_per_step = value,
            .dissolved_h2po4_mol_p_per_step = value,
            .surface = filled(SurfacePhosphateTotals, value),
            .minerals = filled(PhosphateMineralTotals, value),
        },
        .band_phosphate = .{
            .dissolved_hpo4_mol_p_per_step = value,
            .dissolved_h2po4_mol_p_per_step = value,
            .surface = filled(SurfacePhosphateTotals, value),
            .minerals = filled(PhosphateMineralTotals, value),
        },
    };
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn testGeometry() Geometry {
    return .{
        .shared_water_volume_m3 = 2,
        .non_band_ammonium_water_volume_m3 = 3,
        .band_ammonium_water_volume_m3 = 5,
        .non_band_phosphate_water_volume_m3 = 7,
        .band_phosphate_water_volume_m3 = 11,
        .shared_soil_mass_megagrams = 13,
        .non_band_ammonium_soil_mass_megagrams = 17,
        .band_ammonium_soil_mass_megagrams = 19,
    };
}

fn testRates() transformations.Result {
    return .{
        .ammonium = .{
            .non_band_ammonium_mol_n_per_m3_step = 1,
            .band_ammonium_mol_n_per_m3_step = 2,
            .non_band_ammonia_mol_n_per_m3_step = 3,
            .band_ammonia_mol_n_per_m3_step = 4,
        },
        .shared_aqueous = .{
            .hydrogen_mol_per_m3_step = 5,
            .hydroxide_mol_per_m3_step = 6,
            .aluminum_mol_per_m3_step = 7,
            .iron_mol_per_m3_step = 8,
            .calcium_mol_per_m3_step = 9,
            .magnesium_mol_per_m3_step = 10,
            .sodium_mol_per_m3_step = 11,
            .potassium_mol_per_m3_step = 12,
        },
        .non_band_phosphate = testPhosphateRates(13),
        .band_phosphate = testPhosphateRates(24),
    };
}

fn testPhosphateRates(start: f64) transformations.PhosphateZoneTransformations {
    return .{
        .dissolved_hpo4_mol_p_per_m3_step = start,
        .dissolved_h2po4_mol_p_per_m3_step = start + 1,
        .surface = .{
            .hydroxyl_site_source_change_per_step = start + 2,
            .protonated_site_source_change_per_step = start + 3,
            .adsorbed_hpo4_source_change_per_step = start + 4,
            .adsorbed_h2po4_source_change_per_step = start + 5,
        },
        .solids = .{
            .aluminum_phosphate_mol_per_m3_step = start + 6,
            .iron_phosphate_mol_per_m3_step = start + 7,
            .dicalcium_phosphate_mol_per_m3_step = start + 8,
            .hydroxyapatite_mol_per_m3_step = start + 9,
            .monocalcium_phosphate_mol_per_m3_step = start + 10,
        },
    };
}

fn testCationRates() cation_exchange.Cations {
    return .{
        .ammonium_non_band = 1,
        .ammonium_band = 2,
        .hydrogen = 3,
        .aluminum = 4,
        .iron = 5,
        .calcium = 6,
        .magnesium = 7,
        .sodium = 8,
        .potassium = 9,
    };
}

test "fixed-pH extensive accumulation matches every source statement" {
    const inputs: Inputs = .{
        .previous_totals = filledTotals(0.25),
        .geometry = testGeometry(),
        .transformation_rates = testRates(),
        .water_reset_changes = .{
            .hydrogen_mol_per_m3_step = 0.5,
            .hydroxide_mol_per_m3_step = -0.75,
        },
        .cation_exchange_rates_mol_per_megagram_step = testCationRates(),
    };
    const actual = try accumulateSourceOrder(inputs);
    const g = inputs.geometry;
    const r = inputs.transformation_rates;
    const x = inputs.cation_exchange_rates_mol_per_megagram_step;
    var expected = inputs.previous_totals;
    expected.ammonium.non_band_ammonium_mol_n_per_step +=
        r.ammonium.non_band_ammonium_mol_n_per_m3_step *
        g.non_band_ammonium_water_volume_m3;
    expected.ammonium.band_ammonium_mol_n_per_step +=
        r.ammonium.band_ammonium_mol_n_per_m3_step *
        g.band_ammonium_water_volume_m3;
    expected.ammonium.non_band_ammonia_mol_n_per_step +=
        r.ammonium.non_band_ammonia_mol_n_per_m3_step *
        g.non_band_ammonium_water_volume_m3;
    expected.ammonium.band_ammonia_mol_n_per_step +=
        r.ammonium.band_ammonia_mol_n_per_m3_step *
        g.band_ammonium_water_volume_m3;
    expected.shared_aqueous.hydrogen_mol_per_step +=
        (r.shared_aqueous.hydrogen_mol_per_m3_step +
            inputs.water_reset_changes.hydrogen_mol_per_m3_step) *
        g.shared_water_volume_m3;
    expected.shared_aqueous.hydroxide_mol_per_step +=
        (r.shared_aqueous.hydroxide_mol_per_m3_step +
            inputs.water_reset_changes.hydroxide_mol_per_m3_step) *
        g.shared_water_volume_m3;
    expected.shared_aqueous.aluminum_mol_per_step +=
        r.shared_aqueous.aluminum_mol_per_m3_step * g.shared_water_volume_m3;
    expected.shared_aqueous.iron_mol_per_step +=
        r.shared_aqueous.iron_mol_per_m3_step * g.shared_water_volume_m3;
    expected.shared_aqueous.calcium_mol_per_step +=
        r.shared_aqueous.calcium_mol_per_m3_step * g.shared_water_volume_m3;
    expected.shared_aqueous.magnesium_mol_per_step +=
        r.shared_aqueous.magnesium_mol_per_m3_step * g.shared_water_volume_m3;
    expected.shared_aqueous.sodium_mol_per_step +=
        r.shared_aqueous.sodium_mol_per_m3_step * g.shared_water_volume_m3;
    expected.shared_aqueous.potassium_mol_per_step +=
        r.shared_aqueous.potassium_mol_per_m3_step * g.shared_water_volume_m3;
    expected.water_equilibration_hydrogen_mol_per_step +=
        inputs.water_reset_changes.hydrogen_mol_per_m3_step *
        g.shared_water_volume_m3;
    accumulatePhosphateDissolved(
        &expected.non_band_phosphate,
        r.non_band_phosphate,
        g.non_band_phosphate_water_volume_m3,
    );
    accumulatePhosphateDissolved(
        &expected.band_phosphate,
        r.band_phosphate,
        g.band_phosphate_water_volume_m3,
    );
    accumulateCationExchange(&expected.cation_exchange, x, g);
    accumulatePhosphateSurface(
        &expected.non_band_phosphate,
        r.non_band_phosphate,
        g.non_band_phosphate_water_volume_m3,
    );
    accumulatePhosphateSurface(
        &expected.band_phosphate,
        r.band_phosphate,
        g.band_phosphate_water_volume_m3,
    );
    accumulatePhosphateMinerals(
        &expected.non_band_phosphate,
        r.non_band_phosphate,
        g.non_band_phosphate_water_volume_m3,
    );
    accumulatePhosphateMinerals(
        &expected.band_phosphate,
        r.band_phosphate,
        g.band_phosphate_water_volume_m3,
    );

    try std.testing.expectEqualDeep(expected, actual);
}

test "fixed-pH extensive accumulation closes pH reset and water product" {
    var rates = testRates();
    rates.shared_aqueous.hydrogen_mol_per_m3_step = 2.0e-5;
    rates.shared_aqueous.hydroxide_mol_per_m3_step = -1.0e-5;
    const reset = try water_reset.calculateSourceOrder(.{
        .prescribed_soil_ph = 6.5,
        .monovalent_activity_coefficient = 0.8,
        .water_activity_product_mol2_per_m6 = 1.0e-8,
        .current_hydrogen_activity_mol_per_m3 = 2.0e-4,
        .current_hydroxide_activity_mol_per_m3 = 3.0e-5,
        .assembled_hydrogen_change_mol_per_m3_step = rates.shared_aqueous.hydrogen_mol_per_m3_step,
        .assembled_hydroxide_change_mol_per_m3_step = rates.shared_aqueous.hydroxide_mol_per_m3_step,
    });
    const geometry = testGeometry();
    const totals = try accumulateSourceOrder(.{
        .previous_totals = filledTotals(0),
        .geometry = geometry,
        .transformation_rates = rates,
        .water_reset_changes = .{
            .hydrogen_mol_per_m3_step = reset.hydrogen_reset_mol_per_m3_step,
            .hydroxide_mol_per_m3_step = reset.hydroxide_reset_mol_per_m3_step,
        },
        .cation_exchange_rates_mol_per_megagram_step = testCationRates(),
    });

    try std.testing.expectApproxEqAbs(
        reset.target_hydrogen_activity_mol_per_m3 -
            2.0e-4,
        totals.shared_aqueous.hydrogen_mol_per_step /
            geometry.shared_water_volume_m3,
        1.0e-18,
    );
    try std.testing.expectApproxEqAbs(
        reset.target_hydroxide_activity_mol_per_m3 -
            3.0e-5,
        totals.shared_aqueous.hydroxide_mol_per_step /
            geometry.shared_water_volume_m3,
        1.0e-18,
    );
    try std.testing.expectEqual(
        reset.hydrogen_reset_mol_per_m3_step *
            geometry.shared_water_volume_m3,
        totals.water_equilibration_hydrogen_mol_per_step,
    );
}

test "fixed-pH extensive scaling conserves exchange and phosphate transfers" {
    const geometry: Geometry = .{
        .shared_water_volume_m3 = 2,
        .non_band_ammonium_water_volume_m3 = 3,
        .band_ammonium_water_volume_m3 = 5,
        .non_band_phosphate_water_volume_m3 = 7,
        .band_phosphate_water_volume_m3 = 11,
        .shared_soil_mass_megagrams = 13,
        .non_band_ammonium_soil_mass_megagrams = 17,
        .band_ammonium_soil_mass_megagrams = 19,
    };
    const exchange: cation_exchange.Cations = .{
        .ammonium_non_band = 0.02,
        .ammonium_band = -0.03,
        .hydrogen = 0.04,
        .aluminum = -0.05,
        .iron = 0.06,
        .calcium = -0.07,
        .magnesium = 0.08,
        .sodium = -0.09,
        .potassium = 0.1,
    };
    const assembled = try transformations.assembleSourceOrder(.{
        .geometry = .{
            .shared_soil_mass_per_water_volume_megagrams_per_m3 = geometry.shared_soil_mass_megagrams / geometry.shared_water_volume_m3,
            .non_band_ammonium_soil_mass_per_water_volume_megagrams_per_m3 = geometry.non_band_ammonium_soil_mass_megagrams /
                geometry.non_band_ammonium_water_volume_m3,
            .band_ammonium_soil_mass_per_water_volume_megagrams_per_m3 = geometry.band_ammonium_soil_mass_megagrams /
                geometry.band_ammonium_water_volume_m3,
        },
        .ammonium = .{
            .non_band_association_mol_n_per_m3_step = 0,
            .band_association_mol_n_per_m3_step = 0,
        },
        .non_band_phosphate = .{
            .hpo4_hydrogen_association_mol_p_per_m3_step = 0.11,
            .surface = .{
                .h2po4_with_protonated_site_source_extent_per_step = 0.12,
                .h2po4_with_hydroxyl_site_source_extent_per_step = -0.03,
                .hpo4_with_hydroxyl_site_source_extent_per_step = 0.04,
            },
            .minerals = .{
                .aluminum_phosphate_mol_per_m3_step = 0.05,
                .iron_phosphate_mol_per_m3_step = -0.01,
                .dicalcium_phosphate_mol_per_m3_step = 0.02,
                .hydroxyapatite_mol_per_m3_step = 0.03,
                .monocalcium_phosphate_mol_per_m3_step = -0.04,
            },
        },
        .band_phosphate = .{
            .hpo4_hydrogen_association_mol_p_per_m3_step = 0,
            .surface = std.mem.zeroes(transformations.SurfacePhosphateRates),
            .minerals = std.mem.zeroes(transformations.PhosphateMineralRates),
        },
        .cation_exchange_mol_per_megagram_step = exchange,
    });
    const totals = try accumulateSourceOrder(.{
        .previous_totals = filledTotals(0),
        .geometry = geometry,
        .transformation_rates = assembled,
        .water_reset_changes = .{
            .hydrogen_mol_per_m3_step = 0,
            .hydroxide_mol_per_m3_step = 0,
        },
        .cation_exchange_rates_mol_per_megagram_step = exchange,
    });

    try std.testing.expectApproxEqAbs(
        0,
        totals.ammonium.non_band_ammonium_mol_n_per_step +
            totals.cation_exchange.non_band_ammonium_mol_n_per_step,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        0,
        totals.ammonium.band_ammonium_mol_n_per_step +
            totals.cation_exchange.band_ammonium_mol_n_per_step,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        0,
        totals.shared_aqueous.calcium_mol_per_step +
            totals.cation_exchange.calcium_mol_per_step,
        1.0e-15,
    );
    const phosphate = totals.non_band_phosphate;
    const phosphorus_residual =
        phosphate.dissolved_hpo4_mol_p_per_step +
        phosphate.dissolved_h2po4_mol_p_per_step +
        phosphate.surface.adsorbed_hpo4_mol_p_per_step +
        phosphate.surface.adsorbed_h2po4_mol_p_per_step +
        phosphate.minerals.aluminum_phosphate_mol_mineral_per_step +
        phosphate.minerals.iron_phosphate_mol_mineral_per_step +
        phosphate.minerals.dicalcium_phosphate_mol_mineral_per_step +
        3.0 * phosphate.minerals.hydroxyapatite_mol_mineral_per_step +
        2.0 * phosphate.minerals.monocalcium_phosphate_mol_mineral_per_step;
    try std.testing.expectApproxEqAbs(0, phosphorus_residual, 1.0e-15);
}

test "fixed-pH extensive accumulation preserves zero zones and rejects invalid state" {
    var inputs: Inputs = .{
        .previous_totals = filledTotals(0),
        .geometry = testGeometry(),
        .transformation_rates = testRates(),
        .water_reset_changes = .{
            .hydrogen_mol_per_m3_step = 0,
            .hydroxide_mol_per_m3_step = 0,
        },
        .cation_exchange_rates_mol_per_megagram_step = testCationRates(),
    };
    inputs.geometry.non_band_ammonium_water_volume_m3 = 0;
    inputs.geometry.non_band_phosphate_water_volume_m3 = 0;
    inputs.geometry.non_band_ammonium_soil_mass_megagrams = 0;
    const zero_zone = try accumulateSourceOrder(inputs);
    try std.testing.expectEqual(
        @as(f64, 0),
        zero_zone.ammonium.non_band_ammonium_mol_n_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        zero_zone.non_band_phosphate.dissolved_hpo4_mol_p_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        zero_zone.cation_exchange.non_band_ammonium_mol_n_per_step,
    );

    inputs.geometry.shared_water_volume_m3 = -1;
    try std.testing.expectError(
        error.InvalidFixedPhExtensiveGeometry,
        accumulateSourceOrder(inputs),
    );
    inputs.geometry = testGeometry();
    inputs.previous_totals.shared_aqueous.iron_mol_per_step = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteFixedPhExtensiveInput,
        accumulateSourceOrder(inputs),
    );
    inputs.previous_totals = filledTotals(0);
    inputs.transformation_rates.shared_aqueous.iron_mol_per_m3_step =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteFixedPhExtensiveResult,
        accumulateSourceOrder(inputs),
    );
}
