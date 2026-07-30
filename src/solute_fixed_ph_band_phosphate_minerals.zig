const std = @import("std");
const non_band_minerals = @import("solute_fixed_ph_phosphate_minerals.zig");
const phosphate_network = @import("solute_phosphate_network.zig");
const phosphate_precipitation = @import("solute_phosphate_precipitation.zig");

pub const ZoneGeometry = struct {
    band_water_volume_m3: f64,
    minimum_active_water_volume_m3: f64,
};

pub const SharedConcentrations = struct {
    aluminum_mol_per_m3: f64,
    iron_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
};

pub const BandPhosphateState = struct {
    hpo4_concentration_mol_p_per_m3: f64,
    h2po4_concentration_mol_p_per_m3: f64,
    hpo4_activity_mol_p_per_m3: f64,
    h2po4_activity_mol_p_per_m3: f64,
    aluminum_phosphate_solid_mol_per_m3: f64,
    iron_phosphate_solid_mol_per_m3: f64,
    dicalcium_phosphate_solid_mol_per_m3: f64,
    hydroxyapatite_solid_mol_per_m3: f64,
    monocalcium_phosphate_solid_mol_per_m3: f64,
};

pub const Inputs = struct {
    geometry: ZoneGeometry,
    shared_concentrations: SharedConcentrations,
    shared_activities: non_band_minerals.SharedActivities,
    zone: BandPhosphateState,
    products: non_band_minerals.EquilibriumProducts,
    kinetics: non_band_minerals.KineticControls,
    coefficients: non_band_minerals.ActivityCoefficients,
};

pub const Result = struct {
    status: non_band_minerals.ZoneStatus,
    equilibrium_activities: non_band_minerals.EquilibriumActivities,
    extents: phosphate_network.MineralFluxes,
};

/// Direct source-order translation of SOLUTE.F lines 3180--3235.
///
/// Positive extents precipitate and negative extents dissolve. Runtime grid
/// and layer allocation remains caller-owned; this pure compatibility
/// calculation is not bound to the production mineral solver.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validate(inputs);
    if (inputs.geometry.band_water_volume_m3 <=
        inputs.geometry.minimum_active_water_volume_m3)
    {
        return zeroResult();
    }

    const concentrations = inputs.shared_concentrations;
    const a = inputs.shared_activities;
    const z = inputs.zone;
    const p = inputs.products;
    const k = inputs.kinetics;
    const g1 = inputs.coefficients.monovalent;
    const g2 = inputs.coefficients.divalent;
    const fraction = k.substrate_limit_fraction;

    // Band variscite: SOLUTE.F 3202--3206.
    const aluminum_limit =
        fraction * @min(
            concentrations.aluminum_mol_per_m3,
            z.h2po4_concentration_mol_p_per_m3,
        );
    const aluminum_equilibrium =
        p.aluminum_hydroxide_solubility_product /
        std.math.pow(f64, a.hydroxide_mol_per_m3, 3);
    const variscite_h2po4_equilibrium =
        p.aluminum_phosphate_h2po4_product *
        std.math.pow(f64, a.hydrogen_mol_per_m3, 2) /
        aluminum_equilibrium;
    const variscite_extent = sourceExtent(
        z.aluminum_phosphate_solid_mol_per_m3,
        z.h2po4_activity_mol_p_per_m3,
        variscite_h2po4_equilibrium,
        g1,
        aluminum_limit,
        k.maximum_phosphate_precipitation_mol_per_m3_step,
        k.maximum_phosphate_precipitation_mol_per_m3_step,
    );

    // Band strengite: SOLUTE.F 3210--3214.
    const iron_limit =
        fraction * @min(
            concentrations.iron_mol_per_m3,
            z.h2po4_concentration_mol_p_per_m3,
        );
    const iron_equilibrium =
        p.iron_hydroxide_solubility_product /
        std.math.pow(f64, a.hydroxide_mol_per_m3, 3);
    const strengite_h2po4_equilibrium =
        p.iron_phosphate_h2po4_product *
        std.math.pow(f64, a.hydrogen_mol_per_m3, 2) /
        iron_equilibrium;
    const strengite_extent = sourceExtent(
        z.iron_phosphate_solid_mol_per_m3,
        z.h2po4_activity_mol_p_per_m3,
        strengite_h2po4_equilibrium,
        g1,
        iron_limit,
        k.maximum_phosphate_precipitation_mol_per_m3_step,
        k.maximum_phosphate_precipitation_mol_per_m3_step,
    );

    // Band dicalcium phosphate: SOLUTE.F 3218--3221.
    const dicalcium_limit =
        fraction * @min(
            concentrations.calcium_mol_per_m3,
            z.hpo4_concentration_mol_p_per_m3,
        );
    const dicalcium_hpo4_equilibrium =
        p.dicalcium_phosphate_solubility_product / a.calcium_mol_per_m3;
    const dicalcium_extent = sourceExtent(
        z.dicalcium_phosphate_solid_mol_per_m3,
        z.hpo4_activity_mol_p_per_m3,
        dicalcium_hpo4_equilibrium,
        g2,
        dicalcium_limit,
        k.maximum_phosphate_precipitation_mol_per_m3_step,
        k.maximum_general_solute_reaction_mol_per_m3_step,
    );

    // Band hydroxyapatite: SOLUTE.F 3225--3228 retains exponent 0.333 and
    // uses TPZ as both precipitation and dissolution ceiling.
    const h2po4_calcium_limit =
        fraction * @min(
            concentrations.calcium_mol_per_m3,
            z.h2po4_concentration_mol_p_per_m3,
        );
    const hydroxyapatite_h2po4_equilibrium = std.math.pow(
        f64,
        p.hydroxyapatite_h2po4_product *
            std.math.pow(f64, a.hydrogen_mol_per_m3, 7) /
            std.math.pow(f64, a.calcium_mol_per_m3, 5),
        0.333,
    );
    const hydroxyapatite_extent = sourceExtent(
        z.hydroxyapatite_solid_mol_per_m3,
        z.h2po4_activity_mol_p_per_m3,
        hydroxyapatite_h2po4_equilibrium,
        g1,
        h2po4_calcium_limit,
        k.maximum_hydroxide_mineral_reaction_mol_per_m3_step,
        k.maximum_hydroxide_mineral_reaction_mol_per_m3_step,
    );

    // Band monocalcium phosphate: SOLUTE.F 3232--3235.
    const monocalcium_h2po4_equilibrium = std.math.pow(
        f64,
        p.monocalcium_phosphate_solubility_product /
            a.calcium_mol_per_m3,
        0.5,
    );
    const monocalcium_extent = sourceExtent(
        z.monocalcium_phosphate_solid_mol_per_m3,
        z.h2po4_activity_mol_p_per_m3,
        monocalcium_h2po4_equilibrium,
        g1,
        h2po4_calcium_limit,
        k.maximum_phosphate_precipitation_mol_per_m3_step,
        k.maximum_hydroxide_mineral_reaction_mol_per_m3_step,
    );

    const result = Result{
        .status = .active,
        .equilibrium_activities = .{
            .aluminum_mol_per_m3 = aluminum_equilibrium,
            .variscite_h2po4_mol_p_per_m3 = variscite_h2po4_equilibrium,
            .iron_mol_per_m3 = iron_equilibrium,
            .strengite_h2po4_mol_p_per_m3 = strengite_h2po4_equilibrium,
            .dicalcium_hpo4_mol_p_per_m3 = dicalcium_hpo4_equilibrium,
            .hydroxyapatite_h2po4_mol_p_per_m3 = hydroxyapatite_h2po4_equilibrium,
            .monocalcium_h2po4_mol_p_per_m3 = monocalcium_h2po4_equilibrium,
        },
        .extents = .{
            .aluminum_phosphate_mol_per_m3 = variscite_extent,
            .iron_phosphate_mol_per_m3 = strengite_extent,
            .dicalcium_phosphate_mol_per_m3 = dicalcium_extent,
            .hydroxyapatite_mol_per_m3 = hydroxyapatite_extent,
            .monocalcium_phosphate_mol_per_m3 = monocalcium_extent,
        },
    };
    try validateResult(result);
    return result;
}

fn sourceExtent(
    solid_mol_per_m3: f64,
    phosphate_activity_mol_p_per_m3: f64,
    equilibrium_activity_mol_p_per_m3: f64,
    activity_coefficient: f64,
    precipitation_inventory_limit_mol_per_m3: f64,
    precipitation_maximum_mol_per_m3_step: f64,
    dissolution_maximum_mol_per_m3_step: f64,
) f64 {
    return @max(
        -@max(0.0, solid_mol_per_m3),
        -dissolution_maximum_mol_per_m3_step,
        @min(
            precipitation_maximum_mol_per_m3_step,
            precipitation_inventory_limit_mol_per_m3,
            (phosphate_activity_mol_p_per_m3 -
                equilibrium_activity_mol_p_per_m3) /
                activity_coefficient,
        ),
    );
}

fn zeroResult() Result {
    return .{
        .status = .dry,
        .equilibrium_activities = std.mem.zeroes(non_band_minerals.EquilibriumActivities),
        .extents = std.mem.zeroes(phosphate_network.MineralFluxes),
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(ZoneGeometry).@"struct".fields) |field| {
        const value = @field(inputs.geometry, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhBandMineralInput;
    }
    inline for (@typeInfo(SharedConcentrations).@"struct".fields) |field| {
        const value = @field(inputs.shared_concentrations, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhBandMineralInput;
    }
    inline for (@typeInfo(non_band_minerals.SharedActivities).@"struct".fields) |field| {
        const value = @field(inputs.shared_activities, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFixedPhBandMineralInput;
    }
    inline for (@typeInfo(BandPhosphateState).@"struct".fields) |field| {
        const value = @field(inputs.zone, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhBandMineralInput;
    }
    inline for (@typeInfo(non_band_minerals.EquilibriumProducts).@"struct".fields) |field| {
        const value = @field(inputs.products, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFixedPhBandMineralInput;
    }
    if (!std.math.isFinite(inputs.kinetics.substrate_limit_fraction) or
        inputs.kinetics.substrate_limit_fraction < 0 or
        inputs.kinetics.substrate_limit_fraction > 1)
        return error.InvalidFixedPhBandMineralInput;
    inline for (@typeInfo(non_band_minerals.KineticControls).@"struct".fields) |field| {
        if (comptime std.mem.eql(
            u8,
            field.name,
            "substrate_limit_fraction",
        )) continue;
        const value = @field(inputs.kinetics, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhBandMineralInput;
    }
    inline for (@typeInfo(non_band_minerals.ActivityCoefficients).@"struct".fields) |field| {
        const value = @field(inputs.coefficients, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFixedPhBandMineralInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(non_band_minerals.EquilibriumActivities).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.equilibrium_activities, field.name)))
            return error.NonFiniteFixedPhBandMineralResult;
    inline for (@typeInfo(phosphate_network.MineralFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.extents, field.name)))
            return error.NonFiniteFixedPhBandMineralResult;
}

fn validInputs() Inputs {
    return .{
        .geometry = .{
            .band_water_volume_m3 = 0.6,
            .minimum_active_water_volume_m3 = 0.1,
        },
        .shared_concentrations = .{
            .aluminum_mol_per_m3 = 0.9,
            .iron_mol_per_m3 = 0.7,
            .calcium_mol_per_m3 = 0.5,
        },
        .shared_activities = .{
            .hydrogen_mol_per_m3 = 1.1,
            .hydroxide_mol_per_m3 = 0.9,
            .calcium_mol_per_m3 = 1.2,
        },
        .zone = .{
            .hpo4_concentration_mol_p_per_m3 = 0.6,
            .h2po4_concentration_mol_p_per_m3 = 0.8,
            .hpo4_activity_mol_p_per_m3 = 0.3,
            .h2po4_activity_mol_p_per_m3 = 0.56,
            .aluminum_phosphate_solid_mol_per_m3 = 0.1,
            .iron_phosphate_solid_mol_per_m3 = 0.2,
            .dicalcium_phosphate_solid_mol_per_m3 = 0.3,
            .hydroxyapatite_solid_mol_per_m3 = 0.4,
            .monocalcium_phosphate_solid_mol_per_m3 = 0.5,
        },
        .products = .{
            .aluminum_hydroxide_solubility_product = 0.2,
            .aluminum_phosphate_h2po4_product = 0.001,
            .iron_hydroxide_solubility_product = 0.3,
            .iron_phosphate_h2po4_product = 0.001,
            .dicalcium_phosphate_solubility_product = 0.02,
            .hydroxyapatite_h2po4_product = 0.00001,
            .monocalcium_phosphate_solubility_product = 0.012,
        },
        .kinetics = .{
            .substrate_limit_fraction = 0.4,
            .maximum_phosphate_precipitation_mol_per_m3_step = 0.25,
            .maximum_hydroxide_mineral_reaction_mol_per_m3_step = 0.3,
            .maximum_general_solute_reaction_mol_per_m3_step = 0.2,
            .maximum_apatite_reaction_mol_per_m3_step = 0.35,
        },
        .coefficients = .{
            .monovalent = 0.7,
            .divalent = 0.5,
        },
    };
}

test "fixed-pH band minerals match every source equation exactly" {
    const inputs = validInputs();
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(non_band_minerals.ZoneStatus.active, result.status);
    const c = inputs.shared_concentrations;
    const a = inputs.shared_activities;
    const z = inputs.zone;
    const p = inputs.products;
    const k = inputs.kinetics;

    const aluminum_equilibrium =
        p.aluminum_hydroxide_solubility_product /
        std.math.pow(f64, a.hydroxide_mol_per_m3, 3);
    const variscite_equilibrium =
        p.aluminum_phosphate_h2po4_product *
        std.math.pow(f64, a.hydrogen_mol_per_m3, 2) /
        aluminum_equilibrium;
    const iron_equilibrium =
        p.iron_hydroxide_solubility_product /
        std.math.pow(f64, a.hydroxide_mol_per_m3, 3);
    const strengite_equilibrium =
        p.iron_phosphate_h2po4_product *
        std.math.pow(f64, a.hydrogen_mol_per_m3, 2) /
        iron_equilibrium;
    const dicalcium_equilibrium =
        p.dicalcium_phosphate_solubility_product / a.calcium_mol_per_m3;
    const hydroxyapatite_equilibrium = std.math.pow(
        f64,
        p.hydroxyapatite_h2po4_product *
            std.math.pow(f64, a.hydrogen_mol_per_m3, 7) /
            std.math.pow(f64, a.calcium_mol_per_m3, 5),
        0.333,
    );
    const monocalcium_equilibrium = std.math.pow(
        f64,
        p.monocalcium_phosphate_solubility_product / a.calcium_mol_per_m3,
        0.5,
    );
    try std.testing.expectEqual(
        aluminum_equilibrium,
        result.equilibrium_activities.aluminum_mol_per_m3,
    );
    try std.testing.expectEqual(
        variscite_equilibrium,
        result.equilibrium_activities.variscite_h2po4_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        iron_equilibrium,
        result.equilibrium_activities.iron_mol_per_m3,
    );
    try std.testing.expectEqual(
        strengite_equilibrium,
        result.equilibrium_activities.strengite_h2po4_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        dicalcium_equilibrium,
        result.equilibrium_activities.dicalcium_hpo4_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        hydroxyapatite_equilibrium,
        result.equilibrium_activities.hydroxyapatite_h2po4_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        monocalcium_equilibrium,
        result.equilibrium_activities.monocalcium_h2po4_mol_p_per_m3,
    );

    try std.testing.expectEqual(
        sourceExtent(
            z.aluminum_phosphate_solid_mol_per_m3,
            z.h2po4_activity_mol_p_per_m3,
            variscite_equilibrium,
            inputs.coefficients.monovalent,
            k.substrate_limit_fraction *
                @min(c.aluminum_mol_per_m3, z.h2po4_concentration_mol_p_per_m3),
            k.maximum_phosphate_precipitation_mol_per_m3_step,
            k.maximum_phosphate_precipitation_mol_per_m3_step,
        ),
        result.extents.aluminum_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(
        sourceExtent(
            z.iron_phosphate_solid_mol_per_m3,
            z.h2po4_activity_mol_p_per_m3,
            strengite_equilibrium,
            inputs.coefficients.monovalent,
            k.substrate_limit_fraction *
                @min(c.iron_mol_per_m3, z.h2po4_concentration_mol_p_per_m3),
            k.maximum_phosphate_precipitation_mol_per_m3_step,
            k.maximum_phosphate_precipitation_mol_per_m3_step,
        ),
        result.extents.iron_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 0.2),
        result.extents.dicalcium_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 0.2),
        result.extents.hydroxyapatite_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 0.2),
        result.extents.monocalcium_phosphate_mol_per_m3,
    );
}

test "fixed-pH band minerals preserve source dissolution ceilings" {
    var inputs = validInputs();
    inputs.zone.hpo4_activity_mol_p_per_m3 = 0;
    inputs.zone.h2po4_activity_mol_p_per_m3 = 0;
    inputs.zone.aluminum_phosphate_solid_mol_per_m3 = 1;
    inputs.zone.iron_phosphate_solid_mol_per_m3 = 1;
    inputs.zone.dicalcium_phosphate_solid_mol_per_m3 = 1;
    inputs.zone.hydroxyapatite_solid_mol_per_m3 = 1;
    inputs.zone.monocalcium_phosphate_solid_mol_per_m3 = 1;
    inputs.products.aluminum_phosphate_h2po4_product = 100;
    inputs.products.iron_phosphate_h2po4_product = 100;
    inputs.products.dicalcium_phosphate_solubility_product = 100;
    inputs.products.hydroxyapatite_h2po4_product = 100;
    inputs.products.monocalcium_phosphate_solubility_product = 100;
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        -inputs.kinetics.maximum_phosphate_precipitation_mol_per_m3_step,
        result.extents.aluminum_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(
        -inputs.kinetics.maximum_phosphate_precipitation_mol_per_m3_step,
        result.extents.iron_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(
        -inputs.kinetics.maximum_general_solute_reaction_mol_per_m3_step,
        result.extents.dicalcium_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(
        -inputs.kinetics.maximum_hydroxide_mineral_reaction_mol_per_m3_step,
        result.extents.hydroxyapatite_mol_per_m3,
    );
    try std.testing.expectEqual(
        -inputs.kinetics.maximum_hydroxide_mineral_reaction_mol_per_m3_step,
        result.extents.monocalcium_phosphate_mol_per_m3,
    );
}

test "fixed-pH band water gate clears every mineral result" {
    var inputs = validInputs();
    inputs.geometry.band_water_volume_m3 =
        inputs.geometry.minimum_active_water_volume_m3;
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(non_band_minerals.ZoneStatus.dry, result.status);
    try std.testing.expectEqualDeep(
        std.mem.zeroes(non_band_minerals.EquilibriumActivities),
        result.equilibrium_activities,
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(phosphate_network.MineralFluxes),
        result.extents,
    );
}

test "fixed-pH band mineral source bounds expose phosphate stoichiometry overdraw" {
    const source_hydroxyapatite = sourceExtent(0, 1, 0, 1, 0.4 * @min(0.5, 0.8), 1, 1);
    const bounded_hydroxyapatite = try phosphate_precipitation.calculateExtent(
        .{
            .dissolved_cation_mol_per_m3 = 0.5,
            .dissolved_phosphate_mol_p_per_m3 = 0.8,
            .precipitate_mol_per_m3 = 0,
        },
        1,
        0,
        .{
            .cation_mol_per_mol_precipitate = 5,
            .phosphorus_mol_per_mol_precipitate = 3,
        },
        .{
            .substrate_limit_fraction = 0.4,
            .maximum_precipitation_mol_per_m3_step = 1,
            .maximum_dissolution_mol_per_m3_step = 1,
            .phosphate_activity_coefficient = 1,
        },
    );
    try std.testing.expectEqual(@as(f64, 0.2), source_hydroxyapatite);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.04),
        bounded_hydroxyapatite,
        1e-15,
    );

    const source_monocalcium = sourceExtent(0, 1, 0, 1, 0.4 * @min(0.5, 0.8), 1, 1);
    const bounded_monocalcium = try phosphate_precipitation.calculateExtent(
        .{
            .dissolved_cation_mol_per_m3 = 0.5,
            .dissolved_phosphate_mol_p_per_m3 = 0.8,
            .precipitate_mol_per_m3 = 0,
        },
        1,
        0,
        .{
            .cation_mol_per_mol_precipitate = 1,
            .phosphorus_mol_per_mol_precipitate = 2,
        },
        .{
            .substrate_limit_fraction = 0.4,
            .maximum_precipitation_mol_per_m3_step = 1,
            .maximum_dissolution_mol_per_m3_step = 1,
            .phosphate_activity_coefficient = 1,
        },
    );
    try std.testing.expectEqual(@as(f64, 0.2), source_monocalcium);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.16),
        bounded_monocalcium,
        1e-15,
    );
}

test "fixed-pH band minerals reject invalid runtime state" {
    var inputs = validInputs();
    inputs.shared_activities.calcium_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidFixedPhBandMineralInput,
        calculateSourceOrder(inputs),
    );
    inputs = validInputs();
    inputs.shared_concentrations.iron_mol_per_m3 = -1;
    try std.testing.expectError(
        error.InvalidFixedPhBandMineralInput,
        calculateSourceOrder(inputs),
    );
    inputs = validInputs();
    inputs.coefficients.monovalent = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidFixedPhBandMineralInput,
        calculateSourceOrder(inputs),
    );
}
