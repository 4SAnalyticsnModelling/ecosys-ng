const std = @import("std");
const non_band_exchange = @import("solute_fixed_ph_phosphate_exchange.zig");

pub const ZoneGeometry = struct {
    band_water_volume_m3: f64,
    minimum_active_water_volume_m3: f64,
};

pub const KineticControls = struct {
    primary_substrate_limit_fraction: f64,
    hydroxyl_h2po4_substrate_limit_fraction: f64,
    /// Source `TADA` mixes aqueous and surface concentration bases.
    maximum_anion_exchange_source_extent_per_step: f64,
    maximum_general_reaction_mol_per_m3_step: f64,
};

pub const Inputs = struct {
    geometry: ZoneGeometry,
    capacity: non_band_exchange.ExchangeCapacity,
    aqueous: non_band_exchange.AqueousPhosphate,
    sites: non_band_exchange.SurfaceSites,
    shared_activities: non_band_exchange.SharedActivities,
    products: non_band_exchange.EquilibriumProducts,
    kinetics: KineticControls,
    coefficients: non_band_exchange.ActivityCoefficients,
};

pub const Result = struct {
    status: non_band_exchange.ZoneStatus,
    equilibrium_activities: non_band_exchange.EquilibriumActivities,
    extents: non_band_exchange.Extents,
};

/// Direct source-order translation of SOLUTE.F lines 3243--3345.
///
/// The band protonated-site equation deliberately uses the configured water
/// activity product, while the hydroxyl-site H2PO4 equation uses the separate
/// source `FIONX` limiter. This pure comparator is not production-bound.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validate(inputs);
    if (inputs.geometry.band_water_volume_m3 <=
        inputs.geometry.minimum_active_water_volume_m3)
    {
        return zeroResult(.dry);
    }

    const exchange_is_active =
        inputs.capacity.anion_exchange_capacity_mol >
        inputs.capacity.minimum_active_capacity_mol;
    var equilibrium =
        std.mem.zeroes(non_band_exchange.EquilibriumActivities);
    var extents = std.mem.zeroes(non_band_exchange.Extents);
    if (exchange_is_active) {
        try validateActiveSites(inputs.sites);
        const a = inputs.aqueous;
        const s = inputs.sites;
        const shared = inputs.shared_activities;
        const p = inputs.products;
        const maximum =
            inputs.kinetics.maximum_anion_exchange_source_extent_per_step;
        const primary = inputs.kinetics.primary_substrate_limit_fraction;
        const hydroxyl_fraction =
            inputs.kinetics.hydroxyl_h2po4_substrate_limit_fraction;

        // H2PO4 with protonated site: SOLUTE.F 3270--3275.
        equilibrium.h2po4_protonated_site_product =
            p.h2po4_exchange_constant *
            p.water_activity_product_mol2_per_m6;
        equilibrium.h2po4_at_protonated_site_mol_p_per_m3 =
            equilibrium.h2po4_protonated_site_product *
            s.adsorbed_h2po4_mol_p_per_megagram /
            s.protonated_site_mol_per_megagram;
        extents.h2po4_with_protonated_site_source_extent_per_step =
            sourceBoundedExtent(
                a.h2po4_activity_mol_p_per_m3,
                equilibrium.h2po4_at_protonated_site_mol_p_per_m3,
                inputs.coefficients.monovalent,
                primary * s.adsorbed_h2po4_mol_p_per_megagram,
                primary * @min(
                    a.h2po4_concentration_mol_p_per_m3,
                    s.protonated_site_mol_per_megagram,
                ),
                maximum,
            );

        // H2PO4 with hydroxyl site: SOLUTE.F 3276--3281 uses FIONX.
        equilibrium.h2po4_hydroxyl_site_product =
            p.h2po4_exchange_constant *
            shared.hydroxide_mol_per_m3;
        equilibrium.h2po4_at_hydroxyl_site_mol_p_per_m3 =
            equilibrium.h2po4_hydroxyl_site_product *
            s.adsorbed_h2po4_mol_p_per_megagram /
            s.hydroxyl_site_mol_per_megagram;
        extents.h2po4_with_hydroxyl_site_source_extent_per_step =
            sourceBoundedExtent(
                a.h2po4_activity_mol_p_per_m3,
                equilibrium.h2po4_at_hydroxyl_site_mol_p_per_m3,
                inputs.coefficients.monovalent,
                hydroxyl_fraction * s.adsorbed_h2po4_mol_p_per_megagram,
                hydroxyl_fraction * @min(
                    a.h2po4_concentration_mol_p_per_m3,
                    s.hydroxyl_site_mol_per_megagram,
                ),
                maximum,
            );

        // HPO4 with hydroxyl site: SOLUTE.F 3297--3302.
        equilibrium.hpo4_hydroxyl_site_product =
            p.hpo4_exchange_constant *
            p.water_activity_product_mol2_per_m6 /
            p.h2po4_dissociation_constant_mol_per_m3;
        equilibrium.hpo4_at_hydroxyl_site_mol_p_per_m3 =
            equilibrium.hpo4_hydroxyl_site_product *
            s.adsorbed_hpo4_mol_p_per_megagram /
            s.hydroxyl_site_mol_per_megagram;
        extents.hpo4_with_hydroxyl_site_source_extent_per_step =
            sourceBoundedExtent(
                a.hpo4_activity_mol_p_per_m3,
                equilibrium.hpo4_at_hydroxyl_site_mol_p_per_m3,
                inputs.coefficients.divalent,
                primary * s.adsorbed_hpo4_mol_p_per_megagram,
                primary * @min(
                    a.hpo4_concentration_mol_p_per_m3,
                    s.hydroxyl_site_mol_per_megagram,
                ),
                maximum,
            );
    }

    // H2PO4 <-> H + HPO4 remains outside the capacity gate: 3331--3335.
    const primary = inputs.kinetics.primary_substrate_limit_fraction;
    equilibrium.hpo4_from_h2po4_mol_p_per_m3 =
        inputs.products.h2po4_dissociation_constant_mol_per_m3 *
        inputs.aqueous.h2po4_activity_mol_p_per_m3 /
        inputs.shared_activities.hydrogen_mol_per_m3;
    extents.h2po4_to_hydrogen_plus_hpo4_mol_p_per_m3_step =
        sourceBoundedExtent(
            inputs.aqueous.hpo4_activity_mol_p_per_m3,
            equilibrium.hpo4_from_h2po4_mol_p_per_m3,
            inputs.coefficients.divalent,
            primary * inputs.aqueous.h2po4_concentration_mol_p_per_m3,
            primary * inputs.aqueous.hpo4_concentration_mol_p_per_m3,
            inputs.kinetics.maximum_general_reaction_mol_per_m3_step,
        );

    const result = Result{
        .status = if (exchange_is_active) .exchange_active else .aqueous_only,
        .equilibrium_activities = equilibrium,
        .extents = extents,
    };
    try validateResult(result);
    return result;
}

fn sourceBoundedExtent(
    current_activity: f64,
    equilibrium_activity: f64,
    activity_coefficient: f64,
    reverse_inventory_limit: f64,
    forward_inventory_limit: f64,
    maximum: f64,
) f64 {
    return @max(
        -maximum,
        -reverse_inventory_limit,
        @min(
            maximum,
            forward_inventory_limit,
            (current_activity - equilibrium_activity) / activity_coefficient,
        ),
    );
}

fn zeroResult(status: non_band_exchange.ZoneStatus) Result {
    return .{
        .status = status,
        .equilibrium_activities = std.mem.zeroes(non_band_exchange.EquilibriumActivities),
        .extents = std.mem.zeroes(non_band_exchange.Extents),
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(ZoneGeometry).@"struct".fields) |field|
        try finiteNonnegative(@field(inputs.geometry, field.name));
    inline for (@typeInfo(non_band_exchange.ExchangeCapacity).@"struct".fields) |field|
        try finiteNonnegative(@field(inputs.capacity, field.name));
    inline for (@typeInfo(non_band_exchange.AqueousPhosphate).@"struct".fields) |field|
        try finiteNonnegative(@field(inputs.aqueous, field.name));
    inline for (@typeInfo(non_band_exchange.SurfaceSites).@"struct".fields) |field|
        try finiteNonnegative(@field(inputs.sites, field.name));
    inline for (@typeInfo(non_band_exchange.SharedActivities).@"struct".fields) |field| {
        const value = @field(inputs.shared_activities, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFixedPhBandExchangeInput;
    }
    inline for (@typeInfo(non_band_exchange.EquilibriumProducts).@"struct".fields) |field| {
        const value = @field(inputs.products, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFixedPhBandExchangeInput;
    }
    inline for (.{
        inputs.kinetics.primary_substrate_limit_fraction,
        inputs.kinetics.hydroxyl_h2po4_substrate_limit_fraction,
    }) |fraction| {
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.InvalidFixedPhBandExchangeInput;
    }
    inline for (.{
        inputs.kinetics.maximum_anion_exchange_source_extent_per_step,
        inputs.kinetics.maximum_general_reaction_mol_per_m3_step,
    }) |maximum| try finiteNonnegative(maximum);
    inline for (@typeInfo(non_band_exchange.ActivityCoefficients).@"struct".fields) |field| {
        const value = @field(inputs.coefficients, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFixedPhBandExchangeInput;
    }
}

fn finiteNonnegative(value: f64) !void {
    if (!std.math.isFinite(value) or value < 0)
        return error.InvalidFixedPhBandExchangeInput;
}

fn validateActiveSites(sites: non_band_exchange.SurfaceSites) !void {
    if (sites.hydroxyl_site_mol_per_megagram <= 0 or
        sites.protonated_site_mol_per_megagram <= 0)
        return error.InvalidActiveBandPhosphateExchangeSite;
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(non_band_exchange.EquilibriumActivities).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.equilibrium_activities, field.name)))
            return error.NonFiniteFixedPhBandExchangeResult;
    inline for (@typeInfo(non_band_exchange.Extents).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.extents, field.name)))
            return error.NonFiniteFixedPhBandExchangeResult;
}

fn validInputs() Inputs {
    return .{
        .geometry = .{
            .band_water_volume_m3 = 0.6,
            .minimum_active_water_volume_m3 = 0.1,
        },
        .capacity = .{
            .anion_exchange_capacity_mol = 1,
            .minimum_active_capacity_mol = 0.1,
        },
        .aqueous = .{
            .hpo4_concentration_mol_p_per_m3 = 0.1,
            .h2po4_concentration_mol_p_per_m3 = 0.2,
            .hpo4_activity_mol_p_per_m3 = 0.06,
            .h2po4_activity_mol_p_per_m3 = 0.16,
        },
        .sites = .{
            .hydroxyl_site_mol_per_megagram = 0.4,
            .protonated_site_mol_per_megagram = 0.3,
            .adsorbed_hpo4_mol_p_per_megagram = 0.04,
            .adsorbed_h2po4_mol_p_per_megagram = 0.05,
        },
        .shared_activities = .{
            .hydrogen_mol_per_m3 = 0.08,
            .hydroxide_mol_per_m3 = 1.25e-7,
        },
        .products = .{
            .h2po4_exchange_constant = 0.2,
            .hpo4_exchange_constant = 0.3,
            .water_activity_product_mol2_per_m6 = 1.0e-8,
            .h2po4_dissociation_constant_mol_per_m3 = 1.0e-3,
        },
        .kinetics = .{
            .primary_substrate_limit_fraction = 0.5,
            .hydroxyl_h2po4_substrate_limit_fraction = 0.25,
            .maximum_anion_exchange_source_extent_per_step = 1,
            .maximum_general_reaction_mol_per_m3_step = 0.2,
        },
        .coefficients = .{
            .monovalent = 0.8,
            .divalent = 0.6,
        },
    };
}

test "fixed-pH band exchange matches every source equation exactly" {
    const inputs = validInputs();
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        non_band_exchange.ZoneStatus.exchange_active,
        result.status,
    );
    const a = inputs.aqueous;
    const s = inputs.sites;
    const shared = inputs.shared_activities;
    const p = inputs.products;
    const primary = inputs.kinetics.primary_substrate_limit_fraction;
    const secondary =
        inputs.kinetics.hydroxyl_h2po4_substrate_limit_fraction;
    const maximum =
        inputs.kinetics.maximum_anion_exchange_source_extent_per_step;

    const protonated_product =
        p.h2po4_exchange_constant *
        p.water_activity_product_mol2_per_m6;
    const protonated_equilibrium =
        protonated_product * s.adsorbed_h2po4_mol_p_per_megagram /
        s.protonated_site_mol_per_megagram;
    const hydroxyl_product =
        p.h2po4_exchange_constant * shared.hydroxide_mol_per_m3;
    const hydroxyl_equilibrium =
        hydroxyl_product * s.adsorbed_h2po4_mol_p_per_megagram /
        s.hydroxyl_site_mol_per_megagram;
    const hpo4_product =
        p.hpo4_exchange_constant *
        p.water_activity_product_mol2_per_m6 /
        p.h2po4_dissociation_constant_mol_per_m3;
    const hpo4_equilibrium =
        hpo4_product * s.adsorbed_hpo4_mol_p_per_megagram /
        s.hydroxyl_site_mol_per_megagram;
    const dissociation_equilibrium =
        p.h2po4_dissociation_constant_mol_per_m3 *
        a.h2po4_activity_mol_p_per_m3 /
        shared.hydrogen_mol_per_m3;
    try std.testing.expectEqual(
        protonated_product,
        result.equilibrium_activities.h2po4_protonated_site_product,
    );
    try std.testing.expectEqual(
        protonated_equilibrium,
        result.equilibrium_activities.h2po4_at_protonated_site_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        hydroxyl_product,
        result.equilibrium_activities.h2po4_hydroxyl_site_product,
    );
    try std.testing.expectEqual(
        hydroxyl_equilibrium,
        result.equilibrium_activities.h2po4_at_hydroxyl_site_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        hpo4_product,
        result.equilibrium_activities.hpo4_hydroxyl_site_product,
    );
    try std.testing.expectEqual(
        hpo4_equilibrium,
        result.equilibrium_activities.hpo4_at_hydroxyl_site_mol_p_per_m3,
    );
    try std.testing.expectEqual(
        dissociation_equilibrium,
        result.equilibrium_activities.hpo4_from_h2po4_mol_p_per_m3,
    );

    try std.testing.expectEqual(
        sourceBoundedExtent(
            a.h2po4_activity_mol_p_per_m3,
            protonated_equilibrium,
            inputs.coefficients.monovalent,
            primary * s.adsorbed_h2po4_mol_p_per_megagram,
            primary * @min(
                a.h2po4_concentration_mol_p_per_m3,
                s.protonated_site_mol_per_megagram,
            ),
            maximum,
        ),
        result.extents.h2po4_with_protonated_site_source_extent_per_step,
    );
    try std.testing.expectEqual(
        sourceBoundedExtent(
            a.h2po4_activity_mol_p_per_m3,
            hydroxyl_equilibrium,
            inputs.coefficients.monovalent,
            secondary * s.adsorbed_h2po4_mol_p_per_megagram,
            secondary * @min(
                a.h2po4_concentration_mol_p_per_m3,
                s.hydroxyl_site_mol_per_megagram,
            ),
            maximum,
        ),
        result.extents.h2po4_with_hydroxyl_site_source_extent_per_step,
    );
    try std.testing.expectEqual(
        sourceBoundedExtent(
            a.hpo4_activity_mol_p_per_m3,
            hpo4_equilibrium,
            inputs.coefficients.divalent,
            primary * s.adsorbed_hpo4_mol_p_per_megagram,
            primary * @min(
                a.hpo4_concentration_mol_p_per_m3,
                s.hydroxyl_site_mol_per_megagram,
            ),
            maximum,
        ),
        result.extents.hpo4_with_hydroxyl_site_source_extent_per_step,
    );
    try std.testing.expectEqual(
        sourceBoundedExtent(
            a.hpo4_activity_mol_p_per_m3,
            dissociation_equilibrium,
            inputs.coefficients.divalent,
            primary * a.h2po4_concentration_mol_p_per_m3,
            primary * a.hpo4_concentration_mol_p_per_m3,
            inputs.kinetics.maximum_general_reaction_mol_per_m3_step,
        ),
        result.extents.h2po4_to_hydrogen_plus_hpo4_mol_p_per_m3_step,
    );
}

test "fixed-pH band capacity gate retains aqueous association" {
    var inputs = validInputs();
    inputs.capacity.anion_exchange_capacity_mol =
        inputs.capacity.minimum_active_capacity_mol;
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        non_band_exchange.ZoneStatus.aqueous_only,
        result.status,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.extents.h2po4_with_protonated_site_source_extent_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.extents.h2po4_with_hydroxyl_site_source_extent_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.05),
        result.extents.h2po4_to_hydrogen_plus_hpo4_mol_p_per_m3_step,
    );
}

test "fixed-pH band dry gate clears every exchange result" {
    var inputs = validInputs();
    inputs.geometry.band_water_volume_m3 =
        inputs.geometry.minimum_active_water_volume_m3;
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(non_band_exchange.ZoneStatus.dry, result.status);
    try std.testing.expectEqualDeep(
        std.mem.zeroes(non_band_exchange.EquilibriumActivities),
        result.equilibrium_activities,
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(non_band_exchange.Extents),
        result.extents,
    );
}

test "fixed-pH band exchange preserves independent H2PO4 fractions" {
    var inputs = validInputs();
    inputs.aqueous.h2po4_activity_mol_p_per_m3 = 0;
    inputs.products.h2po4_exchange_constant = 1.0e8;
    inputs.kinetics.primary_substrate_limit_fraction = 1;
    inputs.kinetics.hydroxyl_h2po4_substrate_limit_fraction = 0.25;
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        -inputs.sites.adsorbed_h2po4_mol_p_per_megagram,
        result.extents.h2po4_with_protonated_site_source_extent_per_step,
    );
    try std.testing.expectEqual(
        -0.25 * inputs.sites.adsorbed_h2po4_mol_p_per_megagram,
        result.extents.h2po4_with_hydroxyl_site_source_extent_per_step,
    );
}

test "fixed-pH band water product can reverse the non-band exchange sign" {
    var band_inputs = validInputs();
    band_inputs.shared_activities.hydroxide_mol_per_m3 = 0.01;
    band_inputs.aqueous.h2po4_activity_mol_p_per_m3 = 1.0e-5;
    band_inputs.kinetics.primary_substrate_limit_fraction = 1;
    const band = try calculateSourceOrder(band_inputs);
    const non_band = try non_band_exchange.calculateSourceOrder(.{
        .geometry = .{
            .non_band_water_volume_m3 = band_inputs.geometry.band_water_volume_m3,
            .minimum_active_water_volume_m3 = band_inputs.geometry.minimum_active_water_volume_m3,
        },
        .capacity = band_inputs.capacity,
        .aqueous = band_inputs.aqueous,
        .sites = band_inputs.sites,
        .shared_activities = band_inputs.shared_activities,
        .products = band_inputs.products,
        .kinetics = .{
            .substrate_limit_fraction = 1,
            .maximum_anion_exchange_source_extent_per_step = band_inputs.kinetics.maximum_anion_exchange_source_extent_per_step,
            .maximum_general_reaction_mol_per_m3_step = band_inputs.kinetics.maximum_general_reaction_mol_per_m3_step,
        },
        .coefficients = band_inputs.coefficients,
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.2499583333333334e-5),
        band.extents.h2po4_with_protonated_site_source_extent_per_step,
        1e-19,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, -2.0833333333333333e-5),
        non_band.extents.h2po4_with_protonated_site_source_extent_per_step,
        1e-19,
    );
}

test "fixed-pH band exchange rejects invalid runtime input" {
    var inputs = validInputs();
    inputs.sites.protonated_site_mol_per_megagram = 0;
    try std.testing.expectError(
        error.InvalidActiveBandPhosphateExchangeSite,
        calculateSourceOrder(inputs),
    );
    inputs = validInputs();
    inputs.kinetics.hydroxyl_h2po4_substrate_limit_fraction = 1.1;
    try std.testing.expectError(
        error.InvalidFixedPhBandExchangeInput,
        calculateSourceOrder(inputs),
    );
    inputs = validInputs();
    inputs.coefficients.monovalent = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidFixedPhBandExchangeInput,
        calculateSourceOrder(inputs),
    );
}
