const std = @import("std");
const aqueous_network = @import("solute_aqueous_network.zig");
const phosphate_network = @import("solute_phosphate_network.zig");
const cation_exchange = @import("solute_cation_exchange.zig");
const geochemistry = @import("solute_geochemistry_network.zig");

pub const SharedTotals = struct {
    aqueous_mol_per_m3: aqueous_network.Transformations,
    water_mol_per_m3: f64,
    water_balance_mol_per_m3: f64,
};

pub const SharedWaterVolumes = struct {
    non_band_ammonium_water_m3: f64,
    band_ammonium_water_m3: f64,
    shared_water_m3: f64,
};

pub const ExchangeTotals = struct {
    cations_mol_per_megagram: cation_exchange.Cations,
    carboxyl_hydrogen_mol_per_megagram: f64,
};

pub const ExchangeSoilMasses = struct {
    non_band_ammonium_soil_megagrams: f64,
    band_ammonium_soil_megagrams: f64,
    shared_exchange_soil_megagrams: f64,
};

/// Accumulates the shared-aqueous portion of SOLUTE.F 2601--2642 in source
/// order. Inputs and totals are one runtime-selected cell/layer; no grid size
/// is compiled into this scalar kernel.
pub fn accumulateShared(
    totals: *SharedTotals,
    iteration: aqueous_network.Transformations,
    water_change_mol_per_m3: f64,
) !void {
    if (!std.math.isFinite(water_change_mol_per_m3) or
        !std.math.isFinite(totals.water_mol_per_m3) or
        !std.math.isFinite(totals.water_balance_mol_per_m3))
        return error.NonFiniteSoluteIterationTotal;
    inline for (@typeInfo(aqueous_network.Transformations).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(totals.aqueous_mol_per_m3, field.name)) or
            !std.math.isFinite(@field(iteration, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    }

    var next = totals.*;
    inline for (.{
        "ammonium_non_band",
        "ammonium_band",
        "ammonia_non_band",
        "ammonia_band",
        "hydrogen",
        "hydroxide",
        "aluminum",
        "iron",
        "calcium",
        "magnesium",
        "sodium",
        "potassium",
        "sulfate",
        "carbonate",
        "bicarbonate",
        "carbon_dioxide",
    }) |name|
        @field(next.aqueous_mol_per_m3, name) += @field(iteration, name);
    next.water_mol_per_m3 += water_change_mol_per_m3;
    inline for (.{
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
    }) |name|
        @field(next.aqueous_mol_per_m3, name) += @field(iteration, name);
    inline for (@typeInfo(aqueous_network.Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next.aqueous_mol_per_m3, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    if (!std.math.isFinite(next.water_mol_per_m3))
        return error.NonFiniteSoluteIterationTotal;
    totals.* = next;
}

/// Converts shared concentration totals to amounts in SOLUTE.F 2763--2801.
/// Hydrogen silicate is deliberately excluded here; the source scales it at
/// line 2822 after both phosphate-zone blocks.
pub fn scaleSharedToAmounts(
    totals: *SharedTotals,
    volumes: SharedWaterVolumes,
) !void {
    inline for (@typeInfo(SharedWaterVolumes).@"struct".fields) |field| {
        const value = @field(volumes, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoluteTransformationVolume;
    }
    inline for (@typeInfo(aqueous_network.Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals.aqueous_mol_per_m3, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    if (!std.math.isFinite(totals.water_mol_per_m3) or
        !std.math.isFinite(totals.water_balance_mol_per_m3))
        return error.NonFiniteSoluteIterationTotal;

    var next = totals.*;
    inline for (.{ "ammonium_non_band", "ammonia_non_band" }) |name|
        @field(next.aqueous_mol_per_m3, name) *=
            volumes.non_band_ammonium_water_m3;
    inline for (.{ "ammonium_band", "ammonia_band" }) |name|
        @field(next.aqueous_mol_per_m3, name) *=
            volumes.band_ammonium_water_m3;
    inline for (.{
        "hydrogen",
        "hydroxide",
        "aluminum",
        "iron",
        "calcium",
        "magnesium",
        "sodium",
        "potassium",
        "sulfate",
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
    }) |name|
        @field(next.aqueous_mol_per_m3, name) *= volumes.shared_water_m3;
    next.water_balance_mol_per_m3 *= volumes.shared_water_m3;
    next.water_mol_per_m3 *= volumes.shared_water_m3;
    inline for (@typeInfo(aqueous_network.Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next.aqueous_mol_per_m3, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    if (!std.math.isFinite(next.water_mol_per_m3) or
        !std.math.isFinite(next.water_balance_mol_per_m3))
        return error.NonFiniteSoluteIterationTotal;
    totals.* = next;
}

/// Converts one zone's phosphate aqueous/pair totals in SOLUTE.F 2802--2821.
/// Invoke for non-band then band with the corresponding runtime water volume.
pub fn scalePhosphateAqueousToAmounts(
    totals: *phosphate_network.Transformations,
    zone_water_m3: f64,
) !void {
    if (!std.math.isFinite(zone_water_m3) or zone_water_m3 < 0)
        return error.InvalidSoluteTransformationVolume;
    inline for (@typeInfo(phosphate_network.Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals.*, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    var next = totals.*;
    inline for (.{
        "dissolved_po4_mol_p_per_m3",
        "dissolved_hpo4_mol_p_per_m3",
        "dissolved_h2po4_mol_p_per_m3",
        "dissolved_h3po4_mol_p_per_m3",
        "iron_hpo4_pair_mol_per_m3",
        "iron_h2po4_pair_mol_per_m3",
        "calcium_po4_pair_mol_per_m3",
        "calcium_hpo4_pair_mol_per_m3",
        "calcium_h2po4_pair_mol_per_m3",
        "magnesium_hpo4_pair_mol_per_m3",
    }) |name|
        @field(next, name) *= zone_water_m3;
    totals.* = next;
}

/// Performs the separately ordered H-silicate conversion at SOLUTE.F 2822.
pub fn scaleHydrogenSilicateToAmount(
    totals: *SharedTotals,
    shared_water_m3: f64,
) !void {
    if (!std.math.isFinite(shared_water_m3) or shared_water_m3 < 0)
        return error.InvalidSoluteTransformationVolume;
    const current = totals.aqueous_mol_per_m3.hydrogen_silicate;
    if (!std.math.isFinite(current))
        return error.NonFiniteSoluteIterationTotal;
    const next = current * shared_water_m3;
    if (!std.math.isFinite(next))
        return error.NonFiniteSoluteIterationTotal;
    totals.aqueous_mol_per_m3.hydrogen_silicate = next;
}

/// Accumulates one phosphate zone from SOLUTE.F 2643--2662. Invoke first for
/// non-band and then for band to preserve source order.
pub fn accumulatePhosphateAqueous(
    totals: *phosphate_network.Transformations,
    iteration: phosphate_network.Transformations,
) !void {
    inline for (@typeInfo(phosphate_network.Transformations).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(totals.*, field.name)) or
            !std.math.isFinite(@field(iteration, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    }
    var next = totals.*;
    inline for (.{
        "dissolved_po4_mol_p_per_m3",
        "dissolved_hpo4_mol_p_per_m3",
        "dissolved_h2po4_mol_p_per_m3",
        "dissolved_h3po4_mol_p_per_m3",
        "iron_hpo4_pair_mol_per_m3",
        "iron_h2po4_pair_mol_per_m3",
        "calcium_po4_pair_mol_per_m3",
        "calcium_hpo4_pair_mol_per_m3",
        "calcium_h2po4_pair_mol_per_m3",
        "magnesium_hpo4_pair_mol_per_m3",
    }) |name|
        @field(next, name) += @field(iteration, name);
    inline for (@typeInfo(phosphate_network.Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    totals.* = next;
}

/// Accumulates SOLUTE.F 2663--2672 in exact cation then carboxyl order.
pub fn accumulateExchange(
    totals: *ExchangeTotals,
    cation_iteration_mol_per_megagram: cation_exchange.Cations,
    carboxyl_hydrogen_iteration_mol_per_megagram: f64,
) !void {
    if (!std.math.isFinite(totals.carboxyl_hydrogen_mol_per_megagram) or
        !std.math.isFinite(carboxyl_hydrogen_iteration_mol_per_megagram))
        return error.NonFiniteSoluteIterationTotal;
    inline for (@typeInfo(cation_exchange.Cations).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(totals.cations_mol_per_megagram, field.name)) or
            !std.math.isFinite(@field(cation_iteration_mol_per_megagram, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    }
    var next = totals.*;
    inline for (.{
        "ammonium_non_band",
        "ammonium_band",
        "hydrogen",
        "aluminum",
        "iron",
        "calcium",
        "magnesium",
        "sodium",
        "potassium",
    }) |name|
        @field(next.cations_mol_per_megagram, name) +=
            @field(cation_iteration_mol_per_megagram, name);
    next.carboxyl_hydrogen_mol_per_megagram +=
        carboxyl_hydrogen_iteration_mol_per_megagram;
    inline for (@typeInfo(cation_exchange.Cations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next.cations_mol_per_megagram, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    if (!std.math.isFinite(next.carboxyl_hydrogen_mol_per_megagram))
        return error.NonFiniteSoluteIterationTotal;
    totals.* = next;
}

/// Converts exchange concentrations to amounts in SOLUTE.F 2823--2832.
/// NH4 zones retain their distinct carriers; all other cations and carboxyl H
/// use the shared runtime exchange-soil mass.
pub fn scaleExchangeToAmounts(
    totals: *ExchangeTotals,
    soil: ExchangeSoilMasses,
) !void {
    inline for (@typeInfo(ExchangeSoilMasses).@"struct".fields) |field| {
        const value = @field(soil, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoluteTransformationMass;
    }
    if (!std.math.isFinite(totals.carboxyl_hydrogen_mol_per_megagram))
        return error.NonFiniteSoluteIterationTotal;
    inline for (@typeInfo(cation_exchange.Cations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals.cations_mol_per_megagram, field.name)))
            return error.NonFiniteSoluteIterationTotal;

    var next = totals.*;
    next.cations_mol_per_megagram.ammonium_non_band *=
        soil.non_band_ammonium_soil_megagrams;
    next.cations_mol_per_megagram.ammonium_band *=
        soil.band_ammonium_soil_megagrams;
    inline for (.{
        "hydrogen",
        "aluminum",
        "iron",
        "calcium",
        "magnesium",
        "sodium",
        "potassium",
    }) |name|
        @field(next.cations_mol_per_megagram, name) *= soil.shared_exchange_soil_megagrams;
    next.carboxyl_hydrogen_mol_per_megagram *= soil.shared_exchange_soil_megagrams;
    inline for (@typeInfo(cation_exchange.Cations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next.cations_mol_per_megagram, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    if (!std.math.isFinite(next.carboxyl_hydrogen_mol_per_megagram))
        return error.NonFiniteSoluteIterationTotal;
    totals.* = next;
}

/// Accumulates one zone's five surface-site totals from SOLUTE.F 2673--2682.
/// Invoke for non-band and then band to preserve source order.
pub fn accumulatePhosphateSurface(
    totals: *phosphate_network.Transformations,
    iteration: phosphate_network.Transformations,
) !void {
    inline for (@typeInfo(phosphate_network.Transformations).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(totals.*, field.name)) or
            !std.math.isFinite(@field(iteration, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    }
    var next = totals.*;
    inline for (.{
        "deprotonated_site_mol_per_megagram",
        "hydroxyl_site_mol_per_megagram",
        "protonated_site_mol_per_megagram",
        "adsorbed_hpo4_mol_p_per_megagram",
        "adsorbed_h2po4_mol_p_per_megagram",
    }) |name|
        @field(next, name) += @field(iteration, name);
    inline for (@typeInfo(phosphate_network.Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    totals.* = next;
}

/// Converts one zone's five phosphate surface totals in SOLUTE.F 2835--2844.
/// Invoke for non-band then band with its runtime water volume.
pub fn scalePhosphateSurfaceToAmounts(
    totals: *phosphate_network.Transformations,
    zone_water_m3: f64,
) !void {
    if (!std.math.isFinite(zone_water_m3) or zone_water_m3 < 0)
        return error.InvalidSoluteTransformationVolume;
    inline for (@typeInfo(phosphate_network.Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals.*, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    var next = totals.*;
    inline for (.{
        "deprotonated_site_mol_per_megagram",
        "hydroxyl_site_mol_per_megagram",
        "protonated_site_mol_per_megagram",
        "adsorbed_hpo4_mol_p_per_megagram",
        "adsorbed_h2po4_mol_p_per_megagram",
    }) |name|
        @field(next, name) *= zone_water_m3;
    totals.* = next;
}

/// Accumulates mineral-completion then silicate totals from SOLUTE.F
/// 2683--2698 in exact source order.
pub fn accumulateGeochemistrySolids(
    totals: *geochemistry.Transformations,
    iteration: geochemistry.Transformations,
) !void {
    inline for (@typeInfo(geochemistry.Transformations).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(totals.*, field.name)) or
            !std.math.isFinite(@field(iteration, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    }
    var next = totals.*;
    inline for (.{
        "gibbsite_solid_mol_per_m3",
        "iron_hydroxide_solid_mol_per_m3",
        "calcite_solid_mol_per_m3",
        "gypsum_solid_mol_per_m3",
        "aluminum_natural_silicate_mol_per_m3",
        "iron_natural_silicate_mol_per_m3",
        "calcium_natural_silicate_mol_per_m3",
        "magnesium_natural_silicate_mol_per_m3",
        "sodium_natural_silicate_mol_per_m3",
        "potassium_natural_silicate_mol_per_m3",
        "aluminum_ground_silicate_mol_per_m3",
        "iron_ground_silicate_mol_per_m3",
        "calcium_ground_silicate_mol_per_m3",
        "magnesium_ground_silicate_mol_per_m3",
        "sodium_ground_silicate_mol_per_m3",
        "potassium_ground_silicate_mol_per_m3",
    }) |name|
        @field(next, name) += @field(iteration, name);
    inline for (@typeInfo(geochemistry.Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    totals.* = next;
}

/// Converts completion-mineral and silicate totals in SOLUTE.F 2845--2860
/// with the shared runtime water volume.
pub fn scaleGeochemistrySolidsToAmounts(
    totals: *geochemistry.Transformations,
    shared_water_m3: f64,
) !void {
    if (!std.math.isFinite(shared_water_m3) or shared_water_m3 < 0)
        return error.InvalidSoluteTransformationVolume;
    inline for (@typeInfo(geochemistry.Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals.*, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    var next = totals.*;
    inline for (.{
        "gibbsite_solid_mol_per_m3",
        "iron_hydroxide_solid_mol_per_m3",
        "calcite_solid_mol_per_m3",
        "gypsum_solid_mol_per_m3",
        "aluminum_natural_silicate_mol_per_m3",
        "iron_natural_silicate_mol_per_m3",
        "calcium_natural_silicate_mol_per_m3",
        "magnesium_natural_silicate_mol_per_m3",
        "sodium_natural_silicate_mol_per_m3",
        "potassium_natural_silicate_mol_per_m3",
        "aluminum_ground_silicate_mol_per_m3",
        "iron_ground_silicate_mol_per_m3",
        "calcium_ground_silicate_mol_per_m3",
        "magnesium_ground_silicate_mol_per_m3",
        "sodium_ground_silicate_mol_per_m3",
        "potassium_ground_silicate_mol_per_m3",
    }) |name|
        @field(next, name) *= shared_water_m3;
    totals.* = next;
}

/// Accumulates one zone's five phosphate-mineral totals from SOLUTE.F
/// 2699--2708. Invoke for non-band then band to retain source order.
pub fn accumulatePhosphateMinerals(
    totals: *phosphate_network.Transformations,
    iteration: phosphate_network.Transformations,
) !void {
    inline for (@typeInfo(phosphate_network.Transformations).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(totals.*, field.name)) or
            !std.math.isFinite(@field(iteration, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    }
    var next = totals.*;
    inline for (.{
        "aluminum_phosphate_solid_mol_per_m3",
        "iron_phosphate_solid_mol_per_m3",
        "dicalcium_phosphate_solid_mol_per_m3",
        "hydroxyapatite_solid_mol_per_m3",
        "monocalcium_phosphate_solid_mol_per_m3",
    }) |name|
        @field(next, name) += @field(iteration, name);
    inline for (@typeInfo(phosphate_network.Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    totals.* = next;
}

/// Converts one zone's five phosphate-mineral totals in SOLUTE.F 2861--2870.
/// Invoke for non-band then band with the corresponding runtime water volume.
pub fn scalePhosphateMineralsToAmounts(
    totals: *phosphate_network.Transformations,
    zone_water_m3: f64,
) !void {
    if (!std.math.isFinite(zone_water_m3) or zone_water_m3 < 0)
        return error.InvalidSoluteTransformationVolume;
    inline for (@typeInfo(phosphate_network.Transformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals.*, field.name)))
            return error.NonFiniteSoluteIterationTotal;
    var next = totals.*;
    inline for (.{
        "aluminum_phosphate_solid_mol_per_m3",
        "iron_phosphate_solid_mol_per_m3",
        "dicalcium_phosphate_solid_mol_per_m3",
        "hydroxyapatite_solid_mol_per_m3",
        "monocalcium_phosphate_solid_mol_per_m3",
    }) |name|
        @field(next, name) *= zone_water_m3;
    totals.* = next;
}

fn filled(value: f64) aqueous_network.Transformations {
    var result: aqueous_network.Transformations = undefined;
    inline for (@typeInfo(aqueous_network.Transformations).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn filledPhosphate(value: f64) phosphate_network.Transformations {
    var result: phosphate_network.Transformations = undefined;
    inline for (@typeInfo(phosphate_network.Transformations).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn filledGeochemistry(value: f64) geochemistry.Transformations {
    var result: geochemistry.Transformations = undefined;
    inline for (@typeInfo(geochemistry.Transformations).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

test "shared iteration totals preserve source accumulation membership" {
    var totals = SharedTotals{
        .aqueous_mol_per_m3 = filled(1),
        .water_mol_per_m3 = 2,
        .water_balance_mol_per_m3 = 3,
    };
    const iteration = filled(0.25);
    try accumulateShared(&totals, iteration, 0.5);

    try std.testing.expectEqual(@as(f64, 1.25), totals.aqueous_mol_per_m3.hydrogen);
    try std.testing.expectEqual(@as(f64, 1.25), totals.aqueous_mol_per_m3.potassium_sulfate);
    try std.testing.expectEqual(@as(f64, 1.25), totals.aqueous_mol_per_m3.hydrogen_silicate);
    try std.testing.expectEqual(@as(f64, 2.5), totals.water_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.aqueous_mol_per_m3.chloride);
    try std.testing.expectEqual(@as(f64, 1), totals.aqueous_mol_per_m3.nitrate_non_band);
    try std.testing.expectEqual(@as(f64, 1), totals.aqueous_mol_per_m3.nitrate_band);
}

test "shared amount conversion preserves source carrier membership" {
    var totals = SharedTotals{
        .aqueous_mol_per_m3 = filled(1),
        .water_mol_per_m3 = 2,
        .water_balance_mol_per_m3 = 3,
    };
    try scaleSharedToAmounts(&totals, .{
        .non_band_ammonium_water_m3 = 4,
        .band_ammonium_water_m3 = 5,
        .shared_water_m3 = 6,
    });
    try std.testing.expectEqual(@as(f64, 4), totals.aqueous_mol_per_m3.ammonium_non_band);
    try std.testing.expectEqual(@as(f64, 5), totals.aqueous_mol_per_m3.ammonia_band);
    try std.testing.expectEqual(@as(f64, 6), totals.aqueous_mol_per_m3.hydrogen);
    try std.testing.expectEqual(@as(f64, 6), totals.aqueous_mol_per_m3.potassium_sulfate);
    try std.testing.expectEqual(@as(f64, 1), totals.aqueous_mol_per_m3.hydrogen_silicate);
    try std.testing.expectEqual(@as(f64, 12), totals.water_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 18), totals.water_balance_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.aqueous_mol_per_m3.chloride);
}

test "phosphate aqueous and later silicate conversions retain source order" {
    var phosphate = filledPhosphate(1);
    try scalePhosphateAqueousToAmounts(&phosphate, 4);
    try std.testing.expectEqual(@as(f64, 4), phosphate.dissolved_po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 4), phosphate.magnesium_hpo4_pair_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), phosphate.deprotonated_site_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 1), phosphate.aluminum_phosphate_solid_mol_per_m3);

    var shared = SharedTotals{
        .aqueous_mol_per_m3 = filled(1),
        .water_mol_per_m3 = 0,
        .water_balance_mol_per_m3 = 0,
    };
    try scaleSharedToAmounts(&shared, .{
        .non_band_ammonium_water_m3 = 2,
        .band_ammonium_water_m3 = 2,
        .shared_water_m3 = 3,
    });
    try std.testing.expectEqual(@as(f64, 1), shared.aqueous_mol_per_m3.hydrogen_silicate);
    try scaleHydrogenSilicateToAmount(&shared, 3);
    try std.testing.expectEqual(@as(f64, 3), shared.aqueous_mol_per_m3.hydrogen_silicate);
}

test "phosphate aqueous totals include exactly ten source pools per zone" {
    var totals = filledPhosphate(1);
    try accumulatePhosphateAqueous(&totals, filledPhosphate(0.25));

    try std.testing.expectEqual(@as(f64, 1.25), totals.dissolved_po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 1.25), totals.iron_h2po4_pair_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1.25), totals.magnesium_hpo4_pair_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.deprotonated_site_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 1), totals.adsorbed_h2po4_mol_p_per_megagram);
    try std.testing.expectEqual(@as(f64, 1), totals.aluminum_phosphate_solid_mol_per_m3);
}

test "exchange totals preserve cation then carboxyl source membership" {
    var totals = ExchangeTotals{
        .cations_mol_per_megagram = .{
            .ammonium_non_band = 1,
            .ammonium_band = 1,
            .hydrogen = 1,
            .aluminum = 1,
            .iron = 1,
            .calcium = 1,
            .magnesium = 1,
            .sodium = 1,
            .potassium = 1,
        },
        .carboxyl_hydrogen_mol_per_megagram = 2,
    };
    const changes = cation_exchange.Cations{
        .ammonium_non_band = 0.1,
        .ammonium_band = 0.2,
        .hydrogen = 0.3,
        .aluminum = 0.4,
        .iron = 0.5,
        .calcium = 0.6,
        .magnesium = 0.7,
        .sodium = 0.8,
        .potassium = 0.9,
    };
    try accumulateExchange(&totals, changes, 0.25);
    try std.testing.expectEqual(@as(f64, 1.1), totals.cations_mol_per_megagram.ammonium_non_band);
    try std.testing.expectEqual(@as(f64, 1.9), totals.cations_mol_per_megagram.potassium);
    try std.testing.expectEqual(@as(f64, 2.25), totals.carboxyl_hydrogen_mol_per_megagram);
}

test "exchange amount conversion uses three source carriers" {
    var totals = ExchangeTotals{
        .cations_mol_per_megagram = .{
            .ammonium_non_band = 1,
            .ammonium_band = 1,
            .hydrogen = 1,
            .aluminum = 1,
            .iron = 1,
            .calcium = 1,
            .magnesium = 1,
            .sodium = 1,
            .potassium = 1,
        },
        .carboxyl_hydrogen_mol_per_megagram = 1,
    };
    try scaleExchangeToAmounts(&totals, .{
        .non_band_ammonium_soil_megagrams = 2,
        .band_ammonium_soil_megagrams = 3,
        .shared_exchange_soil_megagrams = 4,
    });
    try std.testing.expectEqual(@as(f64, 2), totals.cations_mol_per_megagram.ammonium_non_band);
    try std.testing.expectEqual(@as(f64, 3), totals.cations_mol_per_megagram.ammonium_band);
    try std.testing.expectEqual(@as(f64, 4), totals.cations_mol_per_megagram.hydrogen);
    try std.testing.expectEqual(@as(f64, 4), totals.cations_mol_per_megagram.potassium);
    try std.testing.expectEqual(@as(f64, 4), totals.carboxyl_hydrogen_mol_per_megagram);
}

test "phosphate surface totals include exactly five source pools per zone" {
    var totals = filledPhosphate(1);
    try accumulatePhosphateSurface(&totals, filledPhosphate(0.25));

    try std.testing.expectEqual(@as(f64, 1.25), totals.deprotonated_site_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 1.25), totals.protonated_site_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 1.25), totals.adsorbed_h2po4_mol_p_per_megagram);
    try std.testing.expectEqual(@as(f64, 1), totals.dissolved_h2po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.hydroxyapatite_solid_mol_per_m3);
}

test "phosphate surface amount conversion leaves other zone pools untouched" {
    var totals = filledPhosphate(1);
    try scalePhosphateSurfaceToAmounts(&totals, 4);
    try std.testing.expectEqual(@as(f64, 4), totals.deprotonated_site_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 4), totals.adsorbed_h2po4_mol_p_per_megagram);
    try std.testing.expectEqual(@as(f64, 1), totals.dissolved_h2po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.aluminum_phosphate_solid_mol_per_m3);
}

test "geochemistry totals preserve completion then silicate membership" {
    var totals = filledGeochemistry(1);
    try accumulateGeochemistrySolids(&totals, filledGeochemistry(0.25));

    try std.testing.expectEqual(@as(f64, 1.25), totals.gibbsite_solid_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1.25), totals.gypsum_solid_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1.25), totals.aluminum_natural_silicate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1.25), totals.potassium_ground_silicate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.dissolved_hydrogen_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.dissolved_hydrogen_silicate_mol_per_m3);
}

test "geochemistry solid amount conversion leaves dissolved totals untouched" {
    var totals = filledGeochemistry(1);
    try scaleGeochemistrySolidsToAmounts(&totals, 4);
    try std.testing.expectEqual(@as(f64, 4), totals.gibbsite_solid_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4), totals.gypsum_solid_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4), totals.aluminum_natural_silicate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4), totals.potassium_ground_silicate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.dissolved_aluminum_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.dissolved_hydrogen_silicate_mol_per_m3);
}

test "phosphate mineral totals include exactly five source pools per zone" {
    var totals = filledPhosphate(1);
    try accumulatePhosphateMinerals(&totals, filledPhosphate(0.25));

    try std.testing.expectEqual(@as(f64, 1.25), totals.aluminum_phosphate_solid_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1.25), totals.dicalcium_phosphate_solid_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1.25), totals.monocalcium_phosphate_solid_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.dissolved_h2po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.adsorbed_h2po4_mol_p_per_megagram);
}

test "phosphate mineral amount conversion leaves other zone pools untouched" {
    var totals = filledPhosphate(1);
    try scalePhosphateMineralsToAmounts(&totals, 4);
    try std.testing.expectEqual(@as(f64, 4), totals.aluminum_phosphate_solid_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 4), totals.monocalcium_phosphate_solid_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.dissolved_h2po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 1), totals.adsorbed_h2po4_mol_p_per_megagram);
}
