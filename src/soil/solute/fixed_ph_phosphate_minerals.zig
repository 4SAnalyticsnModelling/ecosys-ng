const std = @import("std");
const phosphate_network = @import("phosphate_network.zig");
const phosphate_precipitation = @import("phosphate_precipitation.zig");

pub const ZoneGeometry = struct {
    non_band_water_volume_m3: f64,
    minimum_active_water_volume_m3: f64,
};

pub const SharedActivities = struct {
    hydrogen_mol_per_m3: f64,
    hydroxide_mol_per_m3: f64,
    calcium_mol_per_m3: f64,
};

pub const NonBandPhosphateState = struct {
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

/// Runtime equilibrium products retain the already-combined constants used by
/// the source statements so the compatibility path does not reassociate them.
pub const EquilibriumProducts = struct {
    aluminum_hydroxide_solubility_product: f64,
    aluminum_phosphate_h2po4_product: f64,
    iron_hydroxide_solubility_product: f64,
    iron_phosphate_h2po4_product: f64,
    dicalcium_phosphate_solubility_product: f64,
    hydroxyapatite_h2po4_product: f64,
    monocalcium_phosphate_solubility_product: f64,
};

pub const KineticControls = struct {
    substrate_limit_fraction: f64,
    maximum_phosphate_precipitation_mol_per_m3_step: f64,
    maximum_hydroxide_mineral_reaction_mol_per_m3_step: f64,
    maximum_general_solute_reaction_mol_per_m3_step: f64,
    maximum_apatite_reaction_mol_per_m3_step: f64,
};

pub const ActivityCoefficients = struct {
    monovalent: f64,
    divalent: f64,
};

pub const Inputs = struct {
    geometry: ZoneGeometry,
    shared_activities: SharedActivities,
    zone: NonBandPhosphateState,
    products: EquilibriumProducts,
    kinetics: KineticControls,
    coefficients: ActivityCoefficients,
};

pub const ZoneStatus = enum {
    dry,
    active,
};

pub const EquilibriumActivities = struct {
    aluminum_mol_per_m3: f64,
    variscite_h2po4_mol_p_per_m3: f64,
    iron_mol_per_m3: f64,
    strengite_h2po4_mol_p_per_m3: f64,
    dicalcium_hpo4_mol_p_per_m3: f64,
    hydroxyapatite_h2po4_mol_p_per_m3: f64,
    monocalcium_h2po4_mol_p_per_m3: f64,
};

pub const Result = struct {
    status: ZoneStatus,
    equilibrium_activities: EquilibriumActivities,
    extents: phosphate_network.MineralFluxes,
};

/// Direct source-order translation of SOLUTE.F lines 2987--3041.
///
/// Positive extents precipitate and negative extents dissolve. This
/// diagnostic kernel deliberately retains the source's phosphate-only
/// precipitation inventory bounds. It neither mutates state nor participates
/// in the production mineral solver.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validate(inputs);
    if (inputs.geometry.non_band_water_volume_m3 <=
        inputs.geometry.minimum_active_water_volume_m3)
    {
        return zeroResult();
    }

    const a = inputs.shared_activities;
    const z = inputs.zone;
    const p = inputs.products;
    const k = inputs.kinetics;
    const g1 = inputs.coefficients.monovalent;
    const g2 = inputs.coefficients.divalent;

    // Variscite: SOLUTE.F 2991--2995.
    const phosphate_limit_h2po4 =
        k.substrate_limit_fraction * z.h2po4_concentration_mol_p_per_m3;
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
        phosphate_limit_h2po4,
        k.maximum_phosphate_precipitation_mol_per_m3_step,
        k.maximum_phosphate_precipitation_mol_per_m3_step,
    );

    // Strengite: SOLUTE.F 3004--3008.
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
        phosphate_limit_h2po4,
        k.maximum_phosphate_precipitation_mol_per_m3_step,
        k.maximum_phosphate_precipitation_mol_per_m3_step,
    );

    // Dicalcium phosphate: SOLUTE.F 3017--3020.
    const phosphate_limit_hpo4 =
        k.substrate_limit_fraction * z.hpo4_concentration_mol_p_per_m3;
    const dicalcium_hpo4_equilibrium =
        p.dicalcium_phosphate_solubility_product / a.calcium_mol_per_m3;
    const dicalcium_extent = sourceExtent(
        z.dicalcium_phosphate_solid_mol_per_m3,
        z.hpo4_activity_mol_p_per_m3,
        dicalcium_hpo4_equilibrium,
        g2,
        phosphate_limit_hpo4,
        k.maximum_phosphate_precipitation_mol_per_m3_step,
        k.maximum_general_solute_reaction_mol_per_m3_step,
    );

    // Hydroxyapatite: SOLUTE.F 3024--3027 retains exponent 0.333.
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
        phosphate_limit_h2po4,
        k.maximum_apatite_reaction_mol_per_m3_step,
        k.maximum_apatite_reaction_mol_per_m3_step,
    );

    // Monocalcium phosphate: SOLUTE.F 3038--3041.
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
        phosphate_limit_h2po4,
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
        .equilibrium_activities = std.mem.zeroes(EquilibriumActivities),
        .extents = std.mem.zeroes(phosphate_network.MineralFluxes),
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(ZoneGeometry).@"struct".fields) |field| {
        const value = @field(inputs.geometry, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhPhosphateMineralInput;
    }
    inline for (@typeInfo(SharedActivities).@"struct".fields) |field| {
        const value = @field(inputs.shared_activities, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFixedPhPhosphateMineralInput;
    }
    inline for (@typeInfo(NonBandPhosphateState).@"struct".fields) |field| {
        const value = @field(inputs.zone, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhPhosphateMineralInput;
    }
    inline for (@typeInfo(EquilibriumProducts).@"struct".fields) |field| {
        const value = @field(inputs.products, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFixedPhPhosphateMineralInput;
    }
    if (!std.math.isFinite(inputs.kinetics.substrate_limit_fraction) or
        inputs.kinetics.substrate_limit_fraction < 0 or
        inputs.kinetics.substrate_limit_fraction > 1)
        return error.InvalidFixedPhPhosphateMineralInput;
    inline for (@typeInfo(KineticControls).@"struct".fields) |field| {
        if (comptime std.mem.eql(
            u8,
            field.name,
            "substrate_limit_fraction",
        )) continue;
        const value = @field(inputs.kinetics, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhPhosphateMineralInput;
    }
    inline for (@typeInfo(ActivityCoefficients).@"struct".fields) |field| {
        const value = @field(inputs.coefficients, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFixedPhPhosphateMineralInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(EquilibriumActivities).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.equilibrium_activities, field.name)))
            return error.NonFiniteFixedPhPhosphateMineralResult;
    inline for (@typeInfo(phosphate_network.MineralFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.extents, field.name)))
            return error.NonFiniteFixedPhPhosphateMineralResult;
}

fn validInputs() Inputs {
    return .{
        .geometry = .{
            .non_band_water_volume_m3 = 0.7,
            .minimum_active_water_volume_m3 = 0.1,
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
            .aluminum_phosphate_h2po4_product = 0.4,
            .iron_hydroxide_solubility_product = 0.3,
            .iron_phosphate_h2po4_product = 0.2,
            .dicalcium_phosphate_solubility_product = 0.2,
            .hydroxyapatite_h2po4_product = 0.01,
            .monocalcium_phosphate_solubility_product = 0.12,
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

test "fixed-pH non-band minerals match every source equation exactly" {
    const inputs = validInputs();
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(ZoneStatus.active, result.status);

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

    const h2_limit =
        k.substrate_limit_fraction * z.h2po4_concentration_mol_p_per_m3;
    const hpo4_limit =
        k.substrate_limit_fraction * z.hpo4_concentration_mol_p_per_m3;
    try std.testing.expectEqual(
        sourceExtent(
            z.aluminum_phosphate_solid_mol_per_m3,
            z.h2po4_activity_mol_p_per_m3,
            variscite_equilibrium,
            inputs.coefficients.monovalent,
            h2_limit,
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
            h2_limit,
            k.maximum_phosphate_precipitation_mol_per_m3_step,
            k.maximum_phosphate_precipitation_mol_per_m3_step,
        ),
        result.extents.iron_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 0.24),
        result.extents.dicalcium_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(
        sourceExtent(
            z.hydroxyapatite_solid_mol_per_m3,
            z.h2po4_activity_mol_p_per_m3,
            hydroxyapatite_equilibrium,
            inputs.coefficients.monovalent,
            h2_limit,
            k.maximum_apatite_reaction_mol_per_m3_step,
            k.maximum_apatite_reaction_mol_per_m3_step,
        ),
        result.extents.hydroxyapatite_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 0.25),
        result.extents.monocalcium_phosphate_mol_per_m3,
    );
    try std.testing.expectEqual(@as(f64, 0.24), hpo4_limit);
}

test "fixed-pH non-band minerals preserve distinct dissolution ceilings" {
    var inputs = validInputs();
    inputs.zone.hpo4_activity_mol_p_per_m3 = 0;
    inputs.zone.h2po4_activity_mol_p_per_m3 = 0;
    inputs.products.aluminum_phosphate_h2po4_product = 100;
    inputs.products.iron_phosphate_h2po4_product = 100;
    inputs.products.dicalcium_phosphate_solubility_product = 100;
    inputs.products.hydroxyapatite_h2po4_product = 100;
    inputs.products.monocalcium_phosphate_solubility_product = 100;
    inputs.zone.aluminum_phosphate_solid_mol_per_m3 = 1;
    inputs.zone.iron_phosphate_solid_mol_per_m3 = 1;
    inputs.zone.dicalcium_phosphate_solid_mol_per_m3 = 1;
    inputs.zone.hydroxyapatite_solid_mol_per_m3 = 1;
    inputs.zone.monocalcium_phosphate_solid_mol_per_m3 = 1;
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
        -inputs.kinetics.maximum_apatite_reaction_mol_per_m3_step,
        result.extents.hydroxyapatite_mol_per_m3,
    );
    try std.testing.expectEqual(
        -inputs.kinetics.maximum_hydroxide_mineral_reaction_mol_per_m3_step,
        result.extents.monocalcium_phosphate_mol_per_m3,
    );
}

test "fixed-pH non-band mineral gate clears every dry-zone result" {
    var inputs = validInputs();
    inputs.geometry.non_band_water_volume_m3 =
        inputs.geometry.minimum_active_water_volume_m3;
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(ZoneStatus.dry, result.status);
    try std.testing.expectEqualDeep(
        std.mem.zeroes(EquilibriumActivities),
        result.equilibrium_activities,
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(phosphate_network.MineralFluxes),
        result.extents,
    );
}

test "fixed-pH non-band minerals reject invalid state before evaluation" {
    var inputs = validInputs();
    inputs.shared_activities.hydroxide_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidFixedPhPhosphateMineralInput,
        calculateSourceOrder(inputs),
    );
    inputs = validInputs();
    inputs.zone.monocalcium_phosphate_solid_mol_per_m3 = -1;
    try std.testing.expectError(
        error.InvalidFixedPhPhosphateMineralInput,
        calculateSourceOrder(inputs),
    );
    inputs = validInputs();
    inputs.kinetics.maximum_apatite_reaction_mol_per_m3_step =
        std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidFixedPhPhosphateMineralInput,
        calculateSourceOrder(inputs),
    );
}

test "fixed-pH source phosphate-only bounds expose metal and stoichiometry overdraw" {
    const source_metal_omission = sourceExtent(
        0,
        1,
        0,
        1,
        0.4 * 0.8,
        1,
        1,
    );
    const metal_bounded = try phosphate_precipitation.calculateExtent(
        .{
            .dissolved_cation_mol_per_m3 = 0.001,
            .dissolved_phosphate_mol_p_per_m3 = 0.8,
            .precipitate_mol_per_m3 = 0,
        },
        1,
        0,
        .{
            .cation_mol_per_mol_precipitate = 1,
            .phosphorus_mol_per_mol_precipitate = 1,
        },
        .{
            .substrate_limit_fraction = 0.4,
            .maximum_precipitation_mol_per_m3_step = 1,
            .maximum_dissolution_mol_per_m3_step = 1,
            .phosphate_activity_coefficient = 1,
        },
    );
    try std.testing.expectEqual(@as(f64, 0.32), source_metal_omission);
    try std.testing.expectEqual(@as(f64, 0.0004), metal_bounded);

    const source_two_phosphorus_omission = sourceExtent(
        0,
        1,
        0,
        1,
        0.4 * 0.8,
        1,
        1,
    );
    const stoichiometry_bounded = try phosphate_precipitation.calculateExtent(
        .{
            .dissolved_cation_mol_per_m3 = 10,
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
    try std.testing.expectEqual(
        @as(f64, 0.32),
        source_two_phosphorus_omission,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.16),
        stoichiometry_bounded,
        1e-15,
    );
}
