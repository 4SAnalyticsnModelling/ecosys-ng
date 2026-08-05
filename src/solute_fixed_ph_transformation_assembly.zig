const std = @import("std");
const cation_exchange = @import("solute_cation_exchange.zig");

pub const Geometry = struct {
    shared_soil_mass_per_water_volume_megagrams_per_m3: f64,
    non_band_ammonium_soil_mass_per_water_volume_megagrams_per_m3: f64,
    band_ammonium_soil_mass_per_water_volume_megagrams_per_m3: f64,
};

pub const AmmoniumRates = struct {
    non_band_association_mol_n_per_m3_step: f64,
    band_association_mol_n_per_m3_step: f64,
};

pub const SurfacePhosphateRates = struct {
    h2po4_with_protonated_site_source_extent_per_step: f64,
    h2po4_with_hydroxyl_site_source_extent_per_step: f64,
    hpo4_with_hydroxyl_site_source_extent_per_step: f64,
};

pub const PhosphateMineralRates = struct {
    aluminum_phosphate_mol_per_m3_step: f64,
    iron_phosphate_mol_per_m3_step: f64,
    dicalcium_phosphate_mol_per_m3_step: f64,
    hydroxyapatite_mol_per_m3_step: f64,
    monocalcium_phosphate_mol_per_m3_step: f64,
};

pub const PhosphateZoneRates = struct {
    hpo4_hydrogen_association_mol_p_per_m3_step: f64,
    surface: SurfacePhosphateRates,
    minerals: PhosphateMineralRates,
};

pub const Inputs = struct {
    geometry: Geometry,
    ammonium: AmmoniumRates,
    non_band_phosphate: PhosphateZoneRates,
    band_phosphate: PhosphateZoneRates,
    cation_exchange_mol_per_megagram_step: cation_exchange.Cations,
};

pub const AmmoniumTransformations = struct {
    non_band_ammonium_mol_n_per_m3_step: f64,
    non_band_ammonia_mol_n_per_m3_step: f64,
    band_ammonium_mol_n_per_m3_step: f64,
    band_ammonia_mol_n_per_m3_step: f64,
};

pub const SurfacePhosphateTransformations = struct {
    hydroxyl_site_source_change_per_step: f64,
    protonated_site_source_change_per_step: f64,
    adsorbed_hpo4_source_change_per_step: f64,
    adsorbed_h2po4_source_change_per_step: f64,
};

pub const PhosphateZoneTransformations = struct {
    dissolved_hpo4_mol_p_per_m3_step: f64,
    dissolved_h2po4_mol_p_per_m3_step: f64,
    surface: SurfacePhosphateTransformations,
    solids: PhosphateMineralRates,
};

pub const SharedAqueousTransformations = struct {
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
    ammonium: AmmoniumTransformations,
    non_band_phosphate: PhosphateZoneTransformations,
    band_phosphate: PhosphateZoneTransformations,
    shared_aqueous: SharedAqueousTransformations,
};

/// Direct source-order translation of SOLUTE.F lines 3590--3628.
///
/// Surface extents retain the source's mixed numerical basis. This pure
/// comparator does not apply the production soil-mass:water correction.
pub fn assembleSourceOrder(inputs: Inputs) !Result {
    try validate(inputs);
    const exchange = inputs.cation_exchange_mol_per_megagram_step;
    const shared_density =
        inputs.geometry.shared_soil_mass_per_water_volume_megagrams_per_m3;
    const result = Result{
        .ammonium = .{
            .non_band_ammonium_mol_n_per_m3_step = inputs.ammonium.non_band_association_mol_n_per_m3_step -
                exchange.ammonium_non_band *
                    inputs.geometry
                        .non_band_ammonium_soil_mass_per_water_volume_megagrams_per_m3,
            .band_ammonium_mol_n_per_m3_step = inputs.ammonium.band_association_mol_n_per_m3_step -
                exchange.ammonium_band *
                    inputs.geometry
                        .band_ammonium_soil_mass_per_water_volume_megagrams_per_m3,
            .non_band_ammonia_mol_n_per_m3_step = -inputs.ammonium.non_band_association_mol_n_per_m3_step,
            .band_ammonia_mol_n_per_m3_step = -inputs.ammonium.band_association_mol_n_per_m3_step,
        },
        .non_band_phosphate = assemblePhosphateZone(inputs.non_band_phosphate),
        .band_phosphate = assemblePhosphateZone(inputs.band_phosphate),
        .shared_aqueous = .{
            .hydrogen_mol_per_m3_step = -exchange.hydrogen * shared_density,
            .aluminum_mol_per_m3_step = -exchange.aluminum * shared_density,
            .iron_mol_per_m3_step = -exchange.iron * shared_density,
            .calcium_mol_per_m3_step = -exchange.calcium * shared_density,
            .magnesium_mol_per_m3_step = -exchange.magnesium * shared_density,
            .sodium_mol_per_m3_step = -exchange.sodium * shared_density,
            .potassium_mol_per_m3_step = -exchange.potassium * shared_density,
            .hydroxide_mol_per_m3_step = 0,
        },
    };
    try validateResult(result);
    return result;
}

fn assemblePhosphateZone(
    rates: PhosphateZoneRates,
) PhosphateZoneTransformations {
    const surface = rates.surface;
    const minerals = rates.minerals;
    return .{
        .dissolved_hpo4_mol_p_per_m3_step = -rates.hpo4_hydrogen_association_mol_p_per_m3_step -
            surface.hpo4_with_hydroxyl_site_source_extent_per_step,
        .dissolved_h2po4_mol_p_per_m3_step = rates.hpo4_hydrogen_association_mol_p_per_m3_step -
            surface.h2po4_with_protonated_site_source_extent_per_step -
            surface.h2po4_with_hydroxyl_site_source_extent_per_step -
            minerals.aluminum_phosphate_mol_per_m3_step -
            minerals.iron_phosphate_mol_per_m3_step -
            minerals.dicalcium_phosphate_mol_per_m3_step -
            2.0 * minerals.monocalcium_phosphate_mol_per_m3_step -
            3.0 * minerals.hydroxyapatite_mol_per_m3_step,
        .surface = .{
            .hydroxyl_site_source_change_per_step = -surface.h2po4_with_hydroxyl_site_source_extent_per_step -
                surface.hpo4_with_hydroxyl_site_source_extent_per_step,
            .protonated_site_source_change_per_step = -surface.h2po4_with_protonated_site_source_extent_per_step,
            .adsorbed_hpo4_source_change_per_step = surface.hpo4_with_hydroxyl_site_source_extent_per_step,
            .adsorbed_h2po4_source_change_per_step = surface.h2po4_with_protonated_site_source_extent_per_step +
                surface.h2po4_with_hydroxyl_site_source_extent_per_step,
        },
        .solids = minerals,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Geometry).@"struct".fields) |field| {
        const value = @field(inputs.geometry, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhTransformationGeometry;
    }
    try validateSignedStruct(inputs.ammonium);
    try validatePhosphateZone(inputs.non_band_phosphate);
    try validatePhosphateZone(inputs.band_phosphate);
    try validateSignedStruct(inputs.cation_exchange_mol_per_megagram_step);
}

fn validatePhosphateZone(zone: PhosphateZoneRates) !void {
    if (!std.math.isFinite(
        zone.hpo4_hydrogen_association_mol_p_per_m3_step,
    ))
        return error.NonFiniteFixedPhTransformationRate;
    try validateSignedStruct(zone.surface);
    try validateSignedStruct(zone.minerals);
}

fn validateSignedStruct(values: anytype) !void {
    inline for (@typeInfo(@TypeOf(values)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(values, field.name)))
            return error.NonFiniteFixedPhTransformationRate;
}

fn validateResult(result: Result) !void {
    try validateSignedStruct(result.ammonium);
    try validatePhosphateResult(result.non_band_phosphate);
    try validatePhosphateResult(result.band_phosphate);
    try validateSignedStruct(result.shared_aqueous);
}

fn validatePhosphateResult(result: PhosphateZoneTransformations) !void {
    if (!std.math.isFinite(result.dissolved_hpo4_mol_p_per_m3_step) or
        !std.math.isFinite(result.dissolved_h2po4_mol_p_per_m3_step))
        return error.NonFiniteFixedPhTransformationResult;
    try validateSignedStruct(result.surface);
    try validateSignedStruct(result.solids);
}

fn zeroCations() cation_exchange.Cations {
    return std.mem.zeroes(cation_exchange.Cations);
}

fn zeroZoneRates() PhosphateZoneRates {
    return .{
        .hpo4_hydrogen_association_mol_p_per_m3_step = 0,
        .surface = std.mem.zeroes(SurfacePhosphateRates),
        .minerals = std.mem.zeroes(PhosphateMineralRates),
    };
}

fn validInputs() Inputs {
    return .{
        .geometry = .{
            .shared_soil_mass_per_water_volume_megagrams_per_m3 = 1.4,
            .non_band_ammonium_soil_mass_per_water_volume_megagrams_per_m3 = 1.2,
            .band_ammonium_soil_mass_per_water_volume_megagrams_per_m3 = 0.8,
        },
        .ammonium = .{
            .non_band_association_mol_n_per_m3_step = 0.11,
            .band_association_mol_n_per_m3_step = -0.07,
        },
        .non_band_phosphate = .{
            .hpo4_hydrogen_association_mol_p_per_m3_step = 0.01,
            .surface = .{
                .h2po4_with_protonated_site_source_extent_per_step = 0.02,
                .h2po4_with_hydroxyl_site_source_extent_per_step = 0.03,
                .hpo4_with_hydroxyl_site_source_extent_per_step = 0.04,
            },
            .minerals = .{
                .aluminum_phosphate_mol_per_m3_step = 0.05,
                .iron_phosphate_mol_per_m3_step = 0.06,
                .dicalcium_phosphate_mol_per_m3_step = 0.07,
                .hydroxyapatite_mol_per_m3_step = 0.08,
                .monocalcium_phosphate_mol_per_m3_step = 0.09,
            },
        },
        .band_phosphate = .{
            .hpo4_hydrogen_association_mol_p_per_m3_step = 0.11,
            .surface = .{
                .h2po4_with_protonated_site_source_extent_per_step = 0.12,
                .h2po4_with_hydroxyl_site_source_extent_per_step = 0.13,
                .hpo4_with_hydroxyl_site_source_extent_per_step = 0.14,
            },
            .minerals = .{
                .aluminum_phosphate_mol_per_m3_step = 0.15,
                .iron_phosphate_mol_per_m3_step = 0.16,
                .dicalcium_phosphate_mol_per_m3_step = 0.17,
                .hydroxyapatite_mol_per_m3_step = 0.18,
                .monocalcium_phosphate_mol_per_m3_step = 0.19,
            },
        },
        .cation_exchange_mol_per_megagram_step = .{
            .ammonium_non_band = 0.02,
            .ammonium_band = -0.03,
            .hydrogen = 0.04,
            .aluminum = -0.05,
            .iron = 0.06,
            .calcium = -0.07,
            .magnesium = 0.08,
            .sodium = -0.09,
            .potassium = 0.1,
        },
    };
}

test "fixed-pH transformation assembly matches every source equation" {
    const inputs = validInputs();
    const result = try assembleSourceOrder(inputs);
    const n = inputs.non_band_phosphate;
    const b = inputs.band_phosphate;
    const x = inputs.cation_exchange_mol_per_megagram_step;
    const shared_density =
        inputs.geometry.shared_soil_mass_per_water_volume_megagrams_per_m3;
    const expected = Result{
        .ammonium = .{
            .non_band_ammonium_mol_n_per_m3_step = inputs.ammonium.non_band_association_mol_n_per_m3_step -
                x.ammonium_non_band *
                    inputs.geometry
                        .non_band_ammonium_soil_mass_per_water_volume_megagrams_per_m3,
            .band_ammonium_mol_n_per_m3_step = inputs.ammonium.band_association_mol_n_per_m3_step -
                x.ammonium_band *
                    inputs.geometry
                        .band_ammonium_soil_mass_per_water_volume_megagrams_per_m3,
            .non_band_ammonia_mol_n_per_m3_step = -inputs.ammonium.non_band_association_mol_n_per_m3_step,
            .band_ammonia_mol_n_per_m3_step = -inputs.ammonium.band_association_mol_n_per_m3_step,
        },
        .non_band_phosphate = .{
            .dissolved_hpo4_mol_p_per_m3_step = -n.hpo4_hydrogen_association_mol_p_per_m3_step -
                n.surface.hpo4_with_hydroxyl_site_source_extent_per_step,
            .dissolved_h2po4_mol_p_per_m3_step = n.hpo4_hydrogen_association_mol_p_per_m3_step -
                n.surface.h2po4_with_protonated_site_source_extent_per_step -
                n.surface.h2po4_with_hydroxyl_site_source_extent_per_step -
                n.minerals.aluminum_phosphate_mol_per_m3_step -
                n.minerals.iron_phosphate_mol_per_m3_step -
                n.minerals.dicalcium_phosphate_mol_per_m3_step -
                2.0 * n.minerals.monocalcium_phosphate_mol_per_m3_step -
                3.0 * n.minerals.hydroxyapatite_mol_per_m3_step,
            .surface = .{
                .hydroxyl_site_source_change_per_step = -n.surface
                    .h2po4_with_hydroxyl_site_source_extent_per_step -
                    n.surface.hpo4_with_hydroxyl_site_source_extent_per_step,
                .protonated_site_source_change_per_step = -n.surface
                    .h2po4_with_protonated_site_source_extent_per_step,
                .adsorbed_hpo4_source_change_per_step = n.surface.hpo4_with_hydroxyl_site_source_extent_per_step,
                .adsorbed_h2po4_source_change_per_step = n.surface
                    .h2po4_with_protonated_site_source_extent_per_step +
                    n.surface
                        .h2po4_with_hydroxyl_site_source_extent_per_step,
            },
            .solids = n.minerals,
        },
        .band_phosphate = .{
            .dissolved_hpo4_mol_p_per_m3_step = -b.hpo4_hydrogen_association_mol_p_per_m3_step -
                b.surface.hpo4_with_hydroxyl_site_source_extent_per_step,
            .dissolved_h2po4_mol_p_per_m3_step = b.hpo4_hydrogen_association_mol_p_per_m3_step -
                b.surface.h2po4_with_protonated_site_source_extent_per_step -
                b.surface.h2po4_with_hydroxyl_site_source_extent_per_step -
                b.minerals.aluminum_phosphate_mol_per_m3_step -
                b.minerals.iron_phosphate_mol_per_m3_step -
                b.minerals.dicalcium_phosphate_mol_per_m3_step -
                2.0 * b.minerals.monocalcium_phosphate_mol_per_m3_step -
                3.0 * b.minerals.hydroxyapatite_mol_per_m3_step,
            .surface = .{
                .hydroxyl_site_source_change_per_step = -b.surface
                    .h2po4_with_hydroxyl_site_source_extent_per_step -
                    b.surface.hpo4_with_hydroxyl_site_source_extent_per_step,
                .protonated_site_source_change_per_step = -b.surface
                    .h2po4_with_protonated_site_source_extent_per_step,
                .adsorbed_hpo4_source_change_per_step = b.surface.hpo4_with_hydroxyl_site_source_extent_per_step,
                .adsorbed_h2po4_source_change_per_step = b.surface
                    .h2po4_with_protonated_site_source_extent_per_step +
                    b.surface
                        .h2po4_with_hydroxyl_site_source_extent_per_step,
            },
            .solids = b.minerals,
        },
        .shared_aqueous = .{
            .hydrogen_mol_per_m3_step = -x.hydrogen * shared_density,
            .aluminum_mol_per_m3_step = -x.aluminum * shared_density,
            .iron_mol_per_m3_step = -x.iron * shared_density,
            .calcium_mol_per_m3_step = -x.calcium * shared_density,
            .magnesium_mol_per_m3_step = -x.magnesium * shared_density,
            .sodium_mol_per_m3_step = -x.sodium * shared_density,
            .potassium_mol_per_m3_step = -x.potassium * shared_density,
            .hydroxide_mol_per_m3_step = 0,
        },
    };
    try std.testing.expectEqualDeep(expected, result);
}

test "fixed-pH assembly closes internal N P sites and cation exchange" {
    var inputs = validInputs();
    inputs.non_band_phosphate.minerals =
        std.mem.zeroes(PhosphateMineralRates);
    inputs.band_phosphate.minerals =
        std.mem.zeroes(PhosphateMineralRates);
    const result = try assembleSourceOrder(inputs);
    const non_band_n =
        result.ammonium.non_band_ammonium_mol_n_per_m3_step +
        result.ammonium.non_band_ammonia_mol_n_per_m3_step +
        inputs.cation_exchange_mol_per_megagram_step.ammonium_non_band *
            inputs.geometry
                .non_band_ammonium_soil_mass_per_water_volume_megagrams_per_m3;
    const band_n =
        result.ammonium.band_ammonium_mol_n_per_m3_step +
        result.ammonium.band_ammonia_mol_n_per_m3_step +
        inputs.cation_exchange_mol_per_megagram_step.ammonium_band *
            inputs.geometry
                .band_ammonium_soil_mass_per_water_volume_megagrams_per_m3;
    try std.testing.expectApproxEqAbs(@as(f64, 0), non_band_n, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0), band_n, 1e-15);
    try expectSourcePhosphateClosure(result.non_band_phosphate);
    try expectSourcePhosphateClosure(result.band_phosphate);
    const density =
        inputs.geometry.shared_soil_mass_per_water_volume_megagrams_per_m3;
    inline for (.{
        .{ result.shared_aqueous.hydrogen_mol_per_m3_step, inputs.cation_exchange_mol_per_megagram_step.hydrogen },
        .{ result.shared_aqueous.aluminum_mol_per_m3_step, inputs.cation_exchange_mol_per_megagram_step.aluminum },
        .{ result.shared_aqueous.iron_mol_per_m3_step, inputs.cation_exchange_mol_per_megagram_step.iron },
        .{ result.shared_aqueous.calcium_mol_per_m3_step, inputs.cation_exchange_mol_per_megagram_step.calcium },
        .{ result.shared_aqueous.magnesium_mol_per_m3_step, inputs.cation_exchange_mol_per_megagram_step.magnesium },
        .{ result.shared_aqueous.sodium_mol_per_m3_step, inputs.cation_exchange_mol_per_megagram_step.sodium },
        .{ result.shared_aqueous.potassium_mol_per_m3_step, inputs.cation_exchange_mol_per_megagram_step.potassium },
    }) |pair| try std.testing.expectApproxEqAbs(
        @as(f64, 0),
        pair[0] + pair[1] * density,
        1e-15,
    );
}

fn expectSourcePhosphateClosure(
    zone: PhosphateZoneTransformations,
) !void {
    const surface_p =
        zone.surface.adsorbed_hpo4_source_change_per_step +
        zone.surface.adsorbed_h2po4_source_change_per_step;
    const solid_p =
        zone.solids.aluminum_phosphate_mol_per_m3_step +
        zone.solids.iron_phosphate_mol_per_m3_step +
        zone.solids.dicalcium_phosphate_mol_per_m3_step +
        3.0 * zone.solids.hydroxyapatite_mol_per_m3_step +
        2.0 * zone.solids.monocalcium_phosphate_mol_per_m3_step;
    try std.testing.expectApproxEqAbs(
        @as(f64, 0),
        zone.dissolved_hpo4_mol_p_per_m3_step +
            zone.dissolved_h2po4_mol_p_per_m3_step +
            surface_p + solid_p,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0),
        zone.surface.hydroxyl_site_source_change_per_step +
            zone.surface.protonated_site_source_change_per_step +
            surface_p,
        1e-15,
    );
}

test "fixed-pH assembly exposes omitted dissolved mineral metal changes" {
    var inputs = validInputs();
    inputs.ammonium = std.mem.zeroes(AmmoniumRates);
    inputs.cation_exchange_mol_per_megagram_step = zeroCations();
    inputs.non_band_phosphate = zeroZoneRates();
    inputs.band_phosphate = zeroZoneRates();
    inputs.non_band_phosphate.minerals
        .aluminum_phosphate_mol_per_m3_step = 0.1;
    const result = try assembleSourceOrder(inputs);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.shared_aqueous.aluminum_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.1),
        result.non_band_phosphate.solids
            .aluminum_phosphate_mol_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.1),
        result.shared_aqueous.aluminum_mol_per_m3_step +
            result.non_band_phosphate.solids
                .aluminum_phosphate_mol_per_m3_step,
    );
    try expectSourcePhosphateClosure(result.non_band_phosphate);
}

test "fixed-pH assembly propagates inactive zero rates and rejects NaN" {
    var inputs = validInputs();
    inputs.ammonium = std.mem.zeroes(AmmoniumRates);
    inputs.non_band_phosphate = zeroZoneRates();
    inputs.band_phosphate = zeroZoneRates();
    inputs.cation_exchange_mol_per_megagram_step = zeroCations();
    const result = try assembleSourceOrder(inputs);
    try std.testing.expectEqualDeep(
        std.mem.zeroes(AmmoniumTransformations),
        result.ammonium,
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(SharedAqueousTransformations),
        result.shared_aqueous,
    );

    inputs.ammonium.band_association_mol_n_per_m3_step =
        std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteFixedPhTransformationRate,
        assembleSourceOrder(inputs),
    );
}
