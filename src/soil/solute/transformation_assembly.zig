const std = @import("std");
const aqueous_rates = @import("aqueous_reaction_rates.zig");
const cation_exchange = @import("cation_exchange.zig");
const phosphate_exchange = @import("phosphate_exchange.zig");
const phosphate_network = @import("phosphate_network.zig");
const water_equilibrium = @import("water_equilibrium.zig");

pub const ZonePhosphateFluxes = struct {
    h2po4_dissociation_mol_p_per_m3: f64,
    exchange: phosphate_exchange.RestrictedFlux,
    minerals: phosphate_network.MineralFluxes,
};

pub const SoilMassWaterRatios = struct {
    shared_exchange_megagrams_per_m3: f64,
    non_band_ammonium_exchange_megagrams_per_m3: f64,
    band_ammonium_exchange_megagrams_per_m3: f64,
};

pub const Inputs = struct {
    ammonium: aqueous_rates.RestrictedAmmoniumFluxes,
    cation_exchange_mol_per_megagram: cation_exchange.Cations,
    non_band_phosphate: ZonePhosphateFluxes,
    band_phosphate: ZonePhosphateFluxes,
    soil_mass_water: SoilMassWaterRatios,
};

pub const AmmoniumTransformations = struct {
    ammonium_non_band_mol_n_per_m3: f64,
    ammonium_band_mol_n_per_m3: f64,
    ammonia_non_band_mol_n_per_m3: f64,
    ammonia_band_mol_n_per_m3: f64,
};

pub const ZonePhosphateTransformations = struct {
    hpo4_mol_p_per_m3: f64,
    h2po4_mol_p_per_m3: f64,
    hydroxyl_site_mol_per_megagram: f64,
    protonated_site_mol_per_megagram: f64,
    adsorbed_hpo4_mol_p_per_megagram: f64,
    adsorbed_h2po4_mol_p_per_megagram: f64,
};

pub const SharedCationTransformations = struct {
    hydrogen_mol_per_m3: f64,
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
    magnesium_mol_per_m3: f64,
    sodium_mol_per_m3: f64,
    potassium_mol_per_m3: f64,
    hydroxide_mol_per_m3: f64,
};

pub const Result = struct {
    ammonium: AmmoniumTransformations,
    non_band_phosphate: ZonePhosphateTransformations,
    band_phosphate: ZonePhosphateTransformations,
    shared_cations: SharedCationTransformations,
};

pub const SharedAmounts = struct {
    ammonium_non_band_mol_n: f64,
    ammonium_band_mol_n: f64,
    ammonia_non_band_mol_n: f64,
    ammonia_band_mol_n: f64,
    hydrogen_mol: f64,
    hydroxide_mol: f64,
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    water_balance_mol: f64,
};

pub const SharedWaterVolumes = struct {
    non_band_ammonium_water_m3: f64,
    band_ammonium_water_m3: f64,
    shared_water_m3: f64,
};

pub const ZoneAmounts = struct {
    hpo4_mol_p: f64,
    h2po4_mol_p: f64,
    hydroxyl_site_mol: f64,
    protonated_site_mol: f64,
    adsorbed_hpo4_mol_p: f64,
    adsorbed_h2po4_mol_p: f64,
    aluminum_phosphate_mol: f64,
    iron_phosphate_mol: f64,
    dicalcium_phosphate_mol: f64,
    hydroxyapatite_mol: f64,
    monocalcium_phosphate_mol: f64,
};

pub const DomainAmounts = struct {
    cation_exchange_mol: cation_exchange.Cations,
    non_band_phosphate: ZoneAmounts,
    band_phosphate: ZoneAmounts,
};

pub const DomainCarriers = struct {
    non_band_phosphate_water_m3: f64,
    band_phosphate_water_m3: f64,
    non_band_ammonium_exchange_soil_megagrams: f64,
    band_ammonium_exchange_soil_megagrams: f64,
    shared_exchange_soil_megagrams: f64,
};

/// Direct coupled transformation assembly from SOLUTE.F lines 3603--3628.
pub fn assembleSourceOrder(inputs: Inputs) !Result {
    try validate(inputs);
    const non_band = assemblePhosphateZone(inputs.non_band_phosphate);
    const band = assemblePhosphateZone(inputs.band_phosphate);
    const exchange = inputs.cation_exchange_mol_per_megagram;
    const ratios = inputs.soil_mass_water;
    const result = Result{
        .ammonium = .{
            .ammonium_non_band_mol_n_per_m3 = inputs.ammonium.non_band_association_mol_n_per_m3 - exchange.ammonium_non_band * ratios.non_band_ammonium_exchange_megagrams_per_m3,
            .ammonium_band_mol_n_per_m3 = inputs.ammonium.band_association_mol_n_per_m3 - exchange.ammonium_band * ratios.band_ammonium_exchange_megagrams_per_m3,
            .ammonia_non_band_mol_n_per_m3 = -inputs.ammonium.non_band_association_mol_n_per_m3,
            .ammonia_band_mol_n_per_m3 = -inputs.ammonium.band_association_mol_n_per_m3,
        },
        .non_band_phosphate = non_band,
        .band_phosphate = band,
        .shared_cations = .{
            .hydrogen_mol_per_m3 = -exchange.hydrogen * ratios.shared_exchange_megagrams_per_m3,
            .aluminum_mol_per_m3 = -exchange.aluminum * ratios.shared_exchange_megagrams_per_m3,
            .iron_mol_per_m3 = -exchange.iron * ratios.shared_exchange_megagrams_per_m3,
            .calcium_mol_per_m3 = -exchange.calcium * ratios.shared_exchange_megagrams_per_m3,
            .magnesium_mol_per_m3 = -exchange.magnesium * ratios.shared_exchange_megagrams_per_m3,
            .sodium_mol_per_m3 = -exchange.sodium * ratios.shared_exchange_megagrams_per_m3,
            .potassium_mol_per_m3 = -exchange.potassium * ratios.shared_exchange_megagrams_per_m3,
            .hydroxide_mol_per_m3 = 0,
        },
    };
    return result;
}

/// Accumulates the shared aqueous amounts in SOLUTE.F lines 3682--3694.
/// The pH-reset increments are added to H/OH, while only the H reset is also
/// recorded in the source `TBH2O` water-balance carrier.
pub fn accumulateSharedAmountsSourceOrder(
    totals: *SharedAmounts,
    transformations: Result,
    fixed_ph_reset: water_equilibrium.FixedPhReset,
    volumes: SharedWaterVolumes,
) !void {
    inline for (@typeInfo(SharedAmounts).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals.*, field.name)))
            return error.InvalidSoilTransformationAmount;
    inline for (@typeInfo(SharedWaterVolumes).@"struct".fields) |field| {
        const value = @field(volumes, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoilTransformationAmount;
    }
    inline for (@typeInfo(water_equilibrium.FixedPhReset).@"struct".fields) |field|
        if (!std.math.isFinite(@field(fixed_ph_reset, field.name)))
            return error.InvalidSoilTransformationAmount;

    var next = totals.*;
    next.ammonium_non_band_mol_n += transformations.ammonium.ammonium_non_band_mol_n_per_m3 * volumes.non_band_ammonium_water_m3;
    next.ammonium_band_mol_n += transformations.ammonium.ammonium_band_mol_n_per_m3 * volumes.band_ammonium_water_m3;
    next.ammonia_non_band_mol_n += transformations.ammonium.ammonia_non_band_mol_n_per_m3 * volumes.non_band_ammonium_water_m3;
    next.ammonia_band_mol_n += transformations.ammonium.ammonia_band_mol_n_per_m3 * volumes.band_ammonium_water_m3;
    next.hydrogen_mol += (transformations.shared_cations.hydrogen_mol_per_m3 + fixed_ph_reset.hydrogen_reset_mol_per_m3) * volumes.shared_water_m3;
    next.hydroxide_mol += (transformations.shared_cations.hydroxide_mol_per_m3 + fixed_ph_reset.hydroxide_reset_mol_per_m3) * volumes.shared_water_m3;
    next.aluminum_mol += transformations.shared_cations.aluminum_mol_per_m3 * volumes.shared_water_m3;
    next.iron_mol += transformations.shared_cations.iron_mol_per_m3 * volumes.shared_water_m3;
    next.calcium_mol += transformations.shared_cations.calcium_mol_per_m3 * volumes.shared_water_m3;
    next.magnesium_mol += transformations.shared_cations.magnesium_mol_per_m3 * volumes.shared_water_m3;
    next.sodium_mol += transformations.shared_cations.sodium_mol_per_m3 * volumes.shared_water_m3;
    next.potassium_mol += transformations.shared_cations.potassium_mol_per_m3 * volumes.shared_water_m3;
    next.water_balance_mol += fixed_ph_reset.hydrogen_reset_mol_per_m3 * volumes.shared_water_m3;
    inline for (@typeInfo(SharedAmounts).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.InvalidSoilTransformationAmount;
    totals.* = next;
}

/// Accumulates phosphate, exchange, and mineral amounts in the exact category
/// order and with the carriers used by SOLUTE.F lines 3695--3725.
pub fn accumulateDomainAmountsSourceOrder(
    totals: *DomainAmounts,
    transformations: Result,
    raw_fluxes: Inputs,
    carriers: DomainCarriers,
) !void {
    try validate(raw_fluxes);
    inline for (@typeInfo(cation_exchange.Cations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(totals.cation_exchange_mol, field.name)))
            return error.InvalidSoilTransformationAmount;
    inline for (.{ totals.non_band_phosphate, totals.band_phosphate }) |zone|
        inline for (@typeInfo(ZoneAmounts).@"struct".fields) |field|
            if (!std.math.isFinite(@field(zone, field.name)))
                return error.InvalidSoilTransformationAmount;
    inline for (@typeInfo(DomainCarriers).@"struct".fields) |field| {
        const value = @field(carriers, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSoilTransformationAmount;
    }

    var next = totals.*;
    accumulateZoneAqueous(&next.non_band_phosphate, transformations.non_band_phosphate, carriers.non_band_phosphate_water_m3);
    accumulateZoneAqueous(&next.band_phosphate, transformations.band_phosphate, carriers.band_phosphate_water_m3);
    next.cation_exchange_mol.ammonium_non_band += raw_fluxes.cation_exchange_mol_per_megagram.ammonium_non_band * carriers.non_band_ammonium_exchange_soil_megagrams;
    next.cation_exchange_mol.ammonium_band += raw_fluxes.cation_exchange_mol_per_megagram.ammonium_band * carriers.band_ammonium_exchange_soil_megagrams;
    inline for (.{ "hydrogen", "aluminum", "iron", "calcium", "magnesium", "sodium", "potassium" }) |name|
        @field(next.cation_exchange_mol, name) += @field(raw_fluxes.cation_exchange_mol_per_megagram, name) * carriers.shared_exchange_soil_megagrams;
    accumulateZoneSites(&next.non_band_phosphate, transformations.non_band_phosphate, carriers.non_band_phosphate_water_m3);
    accumulateZoneSites(&next.band_phosphate, transformations.band_phosphate, carriers.band_phosphate_water_m3);
    accumulateZoneMinerals(&next.non_band_phosphate, raw_fluxes.non_band_phosphate.minerals, carriers.non_band_phosphate_water_m3);
    accumulateZoneMinerals(&next.band_phosphate, raw_fluxes.band_phosphate.minerals, carriers.band_phosphate_water_m3);

    inline for (@typeInfo(cation_exchange.Cations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next.cation_exchange_mol, field.name)))
            return error.InvalidSoilTransformationAmount;
    inline for (.{ next.non_band_phosphate, next.band_phosphate }) |zone|
        inline for (@typeInfo(ZoneAmounts).@"struct".fields) |field|
            if (!std.math.isFinite(@field(zone, field.name)))
                return error.InvalidSoilTransformationAmount;
    totals.* = next;
}

fn accumulateZoneAqueous(
    totals: *ZoneAmounts,
    transformations: ZonePhosphateTransformations,
    water_m3: f64,
) void {
    totals.hpo4_mol_p += transformations.hpo4_mol_p_per_m3 * water_m3;
    totals.h2po4_mol_p += transformations.h2po4_mol_p_per_m3 * water_m3;
}

fn accumulateZoneSites(
    totals: *ZoneAmounts,
    transformations: ZonePhosphateTransformations,
    water_m3: f64,
) void {
    totals.hydroxyl_site_mol += transformations.hydroxyl_site_mol_per_megagram * water_m3;
    totals.protonated_site_mol += transformations.protonated_site_mol_per_megagram * water_m3;
    totals.adsorbed_hpo4_mol_p += transformations.adsorbed_hpo4_mol_p_per_megagram * water_m3;
    totals.adsorbed_h2po4_mol_p += transformations.adsorbed_h2po4_mol_p_per_megagram * water_m3;
}

fn accumulateZoneMinerals(
    totals: *ZoneAmounts,
    minerals: phosphate_network.MineralFluxes,
    water_m3: f64,
) void {
    totals.aluminum_phosphate_mol += minerals.aluminum_phosphate_mol_per_m3 * water_m3;
    totals.iron_phosphate_mol += minerals.iron_phosphate_mol_per_m3 * water_m3;
    totals.dicalcium_phosphate_mol += minerals.dicalcium_phosphate_mol_per_m3 * water_m3;
    totals.hydroxyapatite_mol += minerals.hydroxyapatite_mol_per_m3 * water_m3;
    totals.monocalcium_phosphate_mol += minerals.monocalcium_phosphate_mol_per_m3 * water_m3;
}

fn assemblePhosphateZone(fluxes: ZonePhosphateFluxes) ZonePhosphateTransformations {
    const surface = fluxes.exchange;
    const minerals = fluxes.minerals;
    return .{
        .hpo4_mol_p_per_m3 = -fluxes.h2po4_dissociation_mol_p_per_m3 - surface.hpo4_with_hydroxyl_site_mol_p_per_megagram,
        .h2po4_mol_p_per_m3 = fluxes.h2po4_dissociation_mol_p_per_m3 - surface.h2po4_with_protonated_site_mol_p_per_megagram - surface.h2po4_with_hydroxyl_site_mol_p_per_megagram - minerals.aluminum_phosphate_mol_per_m3 - minerals.iron_phosphate_mol_per_m3 - minerals.dicalcium_phosphate_mol_per_m3 - 2 * minerals.monocalcium_phosphate_mol_per_m3 - 3 * minerals.hydroxyapatite_mol_per_m3,
        .hydroxyl_site_mol_per_megagram = -surface.h2po4_with_hydroxyl_site_mol_p_per_megagram - surface.hpo4_with_hydroxyl_site_mol_p_per_megagram,
        .protonated_site_mol_per_megagram = -surface.h2po4_with_protonated_site_mol_p_per_megagram,
        .adsorbed_hpo4_mol_p_per_megagram = surface.hpo4_with_hydroxyl_site_mol_p_per_megagram,
        .adsorbed_h2po4_mol_p_per_megagram = surface.h2po4_with_protonated_site_mol_p_per_megagram + surface.h2po4_with_hydroxyl_site_mol_p_per_megagram,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(aqueous_rates.RestrictedAmmoniumFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.ammonium, field.name))) return error.InvalidSoilTransformationFlux;
    inline for (@typeInfo(cation_exchange.Cations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.cation_exchange_mol_per_megagram, field.name))) return error.InvalidSoilTransformationFlux;
    inline for (@typeInfo(SoilMassWaterRatios).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs.soil_mass_water, field.name)) or @field(inputs.soil_mass_water, field.name) < 0) return error.InvalidSoilTransformationFlux;
    inline for (.{ inputs.non_band_phosphate, inputs.band_phosphate }) |zone| {
        if (!std.math.isFinite(zone.h2po4_dissociation_mol_p_per_m3)) return error.InvalidSoilTransformationFlux;
        inline for (@typeInfo(phosphate_exchange.RestrictedFlux).@"struct".fields) |field|
            if (!std.math.isFinite(@field(zone.exchange, field.name))) return error.InvalidSoilTransformationFlux;
        inline for (@typeInfo(phosphate_network.MineralFluxes).@"struct".fields) |field|
            if (!std.math.isFinite(@field(zone.minerals, field.name))) return error.InvalidSoilTransformationFlux;
    }
}

test "SOLUTE 3603-3628 assembles coupled transformations in source order" {
    const cations = cation_exchange.Cations{ .ammonium_non_band = 1, .ammonium_band = 2, .hydrogen = 3, .aluminum = 4, .iron = 5, .calcium = 6, .magnesium = 7, .sodium = 8, .potassium = 9 };
    const surface = phosphate_exchange.RestrictedFlux{ .h2po4_with_protonated_site_mol_p_per_megagram = 2, .h2po4_with_hydroxyl_site_mol_p_per_megagram = 3, .hpo4_with_hydroxyl_site_mol_p_per_megagram = 5 };
    const minerals = phosphate_network.MineralFluxes{ .aluminum_phosphate_mol_per_m3 = 7, .iron_phosphate_mol_per_m3 = 11, .dicalcium_phosphate_mol_per_m3 = 13, .hydroxyapatite_mol_per_m3 = 17, .monocalcium_phosphate_mol_per_m3 = 19 };
    const result = try assembleSourceOrder(.{
        .ammonium = .{ .non_band_association_mol_n_per_m3 = 10, .band_association_mol_n_per_m3 = 20 },
        .cation_exchange_mol_per_megagram = cations,
        .non_band_phosphate = .{ .h2po4_dissociation_mol_p_per_m3 = 23, .exchange = surface, .minerals = minerals },
        .band_phosphate = .{ .h2po4_dissociation_mol_p_per_m3 = 29, .exchange = surface, .minerals = minerals },
        .soil_mass_water = .{ .shared_exchange_megagrams_per_m3 = 2, .non_band_ammonium_exchange_megagrams_per_m3 = 3, .band_ammonium_exchange_megagrams_per_m3 = 4 },
    });
    try std.testing.expectEqual(@as(f64, 7), result.ammonium.ammonium_non_band_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, 12), result.ammonium.ammonium_band_mol_n_per_m3);
    try std.testing.expectEqual(@as(f64, -28), result.non_band_phosphate.hpo4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, -102), result.non_band_phosphate.h2po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, -8), result.non_band_phosphate.hydroxyl_site_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 5), result.band_phosphate.adsorbed_h2po4_mol_p_per_megagram);
    try std.testing.expectEqual(@as(f64, -12), result.shared_cations.calcium_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 0), result.shared_cations.hydroxide_mol_per_m3);
}

test "SOLUTE 3682-3694 accumulates shared amounts with distinct water carriers" {
    var transformations = std.mem.zeroes(Result);
    transformations.ammonium = .{
        .ammonium_non_band_mol_n_per_m3 = 2,
        .ammonium_band_mol_n_per_m3 = 3,
        .ammonia_non_band_mol_n_per_m3 = -2,
        .ammonia_band_mol_n_per_m3 = -3,
    };
    transformations.shared_cations = .{
        .hydrogen_mol_per_m3 = 4,
        .hydroxide_mol_per_m3 = 5,
        .aluminum_mol_per_m3 = 6,
        .iron_mol_per_m3 = 7,
        .calcium_mol_per_m3 = 8,
        .magnesium_mol_per_m3 = 9,
        .sodium_mol_per_m3 = 10,
        .potassium_mol_per_m3 = 11,
    };
    var totals: SharedAmounts = undefined;
    inline for (@typeInfo(SharedAmounts).@"struct".fields) |field|
        @field(totals, field.name) = 1;
    try accumulateSharedAmountsSourceOrder(&totals, transformations, .{
        .target_hydrogen_activity_mol_per_m3 = 1,
        .target_hydroxide_activity_mol_per_m3 = 1,
        .hydrogen_reset_mol_per_m3 = 0.5,
        .hydroxide_reset_mol_per_m3 = 0.25,
    }, .{
        .non_band_ammonium_water_m3 = 10,
        .band_ammonium_water_m3 = 20,
        .shared_water_m3 = 30,
    });
    try std.testing.expectEqual(@as(f64, 21), totals.ammonium_non_band_mol_n);
    try std.testing.expectEqual(@as(f64, 61), totals.ammonium_band_mol_n);
    try std.testing.expectEqual(@as(f64, -19), totals.ammonia_non_band_mol_n);
    try std.testing.expectEqual(@as(f64, -59), totals.ammonia_band_mol_n);
    try std.testing.expectEqual(@as(f64, 136), totals.hydrogen_mol);
    try std.testing.expectEqual(@as(f64, 158.5), totals.hydroxide_mol);
    try std.testing.expectEqual(@as(f64, 181), totals.aluminum_mol);
    try std.testing.expectEqual(@as(f64, 16), totals.water_balance_mol);
}

test "SOLUTE 3695-3725 preserves phosphate water and exchange soil carriers" {
    var transformations = std.mem.zeroes(Result);
    transformations.non_band_phosphate.hpo4_mol_p_per_m3 = 2;
    transformations.non_band_phosphate.hydroxyl_site_mol_per_megagram = 3;
    transformations.band_phosphate.hpo4_mol_p_per_m3 = 4;
    var raw_fluxes = std.mem.zeroes(Inputs);
    raw_fluxes.cation_exchange_mol_per_megagram.ammonium_non_band = 5;
    raw_fluxes.cation_exchange_mol_per_megagram.ammonium_band = 6;
    raw_fluxes.cation_exchange_mol_per_megagram.hydrogen = 7;
    raw_fluxes.non_band_phosphate.minerals.aluminum_phosphate_mol_per_m3 = 8;
    raw_fluxes.band_phosphate.minerals.iron_phosphate_mol_per_m3 = 9;
    var totals = std.mem.zeroes(DomainAmounts);
    try accumulateDomainAmountsSourceOrder(&totals, transformations, raw_fluxes, .{
        .non_band_phosphate_water_m3 = 10,
        .band_phosphate_water_m3 = 20,
        .non_band_ammonium_exchange_soil_megagrams = 30,
        .band_ammonium_exchange_soil_megagrams = 40,
        .shared_exchange_soil_megagrams = 50,
    });
    try std.testing.expectEqual(@as(f64, 20), totals.non_band_phosphate.hpo4_mol_p);
    try std.testing.expectEqual(@as(f64, 30), totals.non_band_phosphate.hydroxyl_site_mol);
    try std.testing.expectEqual(@as(f64, 80), totals.band_phosphate.hpo4_mol_p);
    try std.testing.expectEqual(@as(f64, 150), totals.cation_exchange_mol.ammonium_non_band);
    try std.testing.expectEqual(@as(f64, 240), totals.cation_exchange_mol.ammonium_band);
    try std.testing.expectEqual(@as(f64, 350), totals.cation_exchange_mol.hydrogen);
    try std.testing.expectEqual(@as(f64, 80), totals.non_band_phosphate.aluminum_phosphate_mol);
    try std.testing.expectEqual(@as(f64, 180), totals.band_phosphate.iron_phosphate_mol);
}
