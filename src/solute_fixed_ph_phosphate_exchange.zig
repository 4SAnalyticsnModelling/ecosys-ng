const std = @import("std");

pub const ZoneGeometry = struct {
    non_band_water_volume_m3: f64,
    minimum_active_water_volume_m3: f64,
};

pub const ExchangeCapacity = struct {
    anion_exchange_capacity_mol: f64,
    minimum_active_capacity_mol: f64,
};

pub const AqueousPhosphate = struct {
    hpo4_concentration_mol_p_per_m3: f64,
    h2po4_concentration_mol_p_per_m3: f64,
    hpo4_activity_mol_p_per_m3: f64,
    h2po4_activity_mol_p_per_m3: f64,
};

pub const SurfaceSites = struct {
    hydroxyl_site_mol_per_Mg: f64,
    protonated_site_mol_per_Mg: f64,
    adsorbed_hpo4_mol_p_per_Mg: f64,
    adsorbed_h2po4_mol_p_per_Mg: f64,
};

pub const SharedActivities = struct {
    hydrogen_mol_per_m3: f64,
    hydroxide_mol_per_m3: f64,
};

pub const EquilibriumProducts = struct {
    h2po4_exchange_constant: f64,
    hpo4_exchange_constant: f64,
    water_activity_product_mol2_per_m6: f64,
    h2po4_dissociation_constant_mol_per_m3: f64,
};

pub const KineticControls = struct {
    substrate_limit_fraction: f64,
    /// Source `TADA` is combined with both mol m-3 aqueous concentrations and
    /// mol Mg-1 surface concentrations, so no coherent unit is assigned here.
    maximum_anion_exchange_source_extent_per_step: f64,
    maximum_general_reaction_mol_per_m3_step: f64,
};

pub const ActivityCoefficients = struct {
    monovalent: f64,
    divalent: f64,
};

pub const Inputs = struct {
    geometry: ZoneGeometry,
    capacity: ExchangeCapacity,
    aqueous: AqueousPhosphate,
    sites: SurfaceSites,
    shared_activities: SharedActivities,
    products: EquilibriumProducts,
    kinetics: KineticControls,
    coefficients: ActivityCoefficients,
};

pub const ZoneStatus = enum {
    dry,
    aqueous_only,
    exchange_active,
};

pub const EquilibriumActivities = struct {
    h2po4_protonated_site_product: f64,
    h2po4_at_protonated_site_mol_p_per_m3: f64,
    h2po4_hydroxyl_site_product: f64,
    h2po4_at_hydroxyl_site_mol_p_per_m3: f64,
    hpo4_hydroxyl_site_product: f64,
    hpo4_at_hydroxyl_site_mol_p_per_m3: f64,
    hpo4_from_h2po4_mol_p_per_m3: f64,
};

/// The three surface values retain source-native numerical units because
/// SOLUTE combines mol m-3 solution and mol Mg-1 site terms. The aqueous
/// association extent has the declared mol P m-3 step-1 unit.
pub const Extents = struct {
    h2po4_with_protonated_site_source_extent_per_step: f64,
    h2po4_with_hydroxyl_site_source_extent_per_step: f64,
    hpo4_with_hydroxyl_site_source_extent_per_step: f64,
    h2po4_to_hydrogen_plus_hpo4_mol_p_per_m3_step: f64,
};

pub const Result = struct {
    status: ZoneStatus,
    equilibrium_activities: EquilibriumActivities,
    extents: Extents,
};

/// Direct source-order translation of SOLUTE.F lines 3050--3168.
///
/// Positive surface extents adsorb phosphate; negative values desorb.
/// Positive aqueous extent associates H + HPO4 to H2PO4; negative values
/// dissociate H2PO4. The function is pure and is not production-bound.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validate(inputs);
    if (inputs.geometry.non_band_water_volume_m3 <=
        inputs.geometry.minimum_active_water_volume_m3)
    {
        return zeroResult(.dry);
    }

    const exchange_is_active =
        inputs.capacity.anion_exchange_capacity_mol >
        inputs.capacity.minimum_active_capacity_mol;
    var equilibrium = std.mem.zeroes(EquilibriumActivities);
    var extents = std.mem.zeroes(Extents);
    if (exchange_is_active) {
        try validateActiveSites(inputs.sites);
        const a = inputs.aqueous;
        const s = inputs.sites;
        const shared = inputs.shared_activities;
        const p = inputs.products;
        const k = inputs.kinetics;
        const maximum = k.maximum_anion_exchange_source_extent_per_step;
        const fraction = k.substrate_limit_fraction;

        // H2PO4 with protonated site: SOLUTE.F 3080--3085.
        const adsorbed_h2po4_limit =
            fraction * s.adsorbed_h2po4_mol_p_per_Mg;
        const protonated_site_limit =
            fraction * @min(
                a.h2po4_concentration_mol_p_per_m3,
                s.protonated_site_mol_per_Mg,
            );
        equilibrium.h2po4_protonated_site_product =
            p.h2po4_exchange_constant *
            shared.hydrogen_mol_per_m3 *
            shared.hydroxide_mol_per_m3;
        equilibrium.h2po4_at_protonated_site_mol_p_per_m3 =
            equilibrium.h2po4_protonated_site_product *
            s.adsorbed_h2po4_mol_p_per_Mg /
            s.protonated_site_mol_per_Mg;
        extents.h2po4_with_protonated_site_source_extent_per_step =
            sourceBoundedExtent(
                a.h2po4_activity_mol_p_per_m3,
                equilibrium.h2po4_at_protonated_site_mol_p_per_m3,
                inputs.coefficients.monovalent,
                adsorbed_h2po4_limit,
                protonated_site_limit,
                maximum,
            );

        // H2PO4 with hydroxyl site: SOLUTE.F 3086--3091.
        const hydroxyl_site_h2po4_limit =
            fraction * @min(
                a.h2po4_concentration_mol_p_per_m3,
                s.hydroxyl_site_mol_per_Mg,
            );
        equilibrium.h2po4_hydroxyl_site_product =
            p.h2po4_exchange_constant *
            shared.hydroxide_mol_per_m3;
        equilibrium.h2po4_at_hydroxyl_site_mol_p_per_m3 =
            equilibrium.h2po4_hydroxyl_site_product *
            s.adsorbed_h2po4_mol_p_per_Mg /
            s.hydroxyl_site_mol_per_Mg;
        extents.h2po4_with_hydroxyl_site_source_extent_per_step =
            sourceBoundedExtent(
                a.h2po4_activity_mol_p_per_m3,
                equilibrium.h2po4_at_hydroxyl_site_mol_p_per_m3,
                inputs.coefficients.monovalent,
                adsorbed_h2po4_limit,
                hydroxyl_site_h2po4_limit,
                maximum,
            );

        // HPO4 with hydroxyl site: SOLUTE.F 3111--3116.
        const adsorbed_hpo4_limit =
            fraction * s.adsorbed_hpo4_mol_p_per_Mg;
        const hydroxyl_site_hpo4_limit =
            fraction * @min(
                a.hpo4_concentration_mol_p_per_m3,
                s.hydroxyl_site_mol_per_Mg,
            );
        equilibrium.hpo4_hydroxyl_site_product =
            p.hpo4_exchange_constant *
            p.water_activity_product_mol2_per_m6 /
            p.h2po4_dissociation_constant_mol_per_m3;
        equilibrium.hpo4_at_hydroxyl_site_mol_p_per_m3 =
            equilibrium.hpo4_hydroxyl_site_product *
            s.adsorbed_hpo4_mol_p_per_Mg /
            s.hydroxyl_site_mol_per_Mg;
        extents.hpo4_with_hydroxyl_site_source_extent_per_step =
            sourceBoundedExtent(
                a.hpo4_activity_mol_p_per_m3,
                equilibrium.hpo4_at_hydroxyl_site_mol_p_per_m3,
                inputs.coefficients.divalent,
                adsorbed_hpo4_limit,
                hydroxyl_site_hpo4_limit,
                maximum,
            );
    }

    // H2PO4 <-> H + HPO4 is outside the capacity gate: lines 3154--3158.
    const fraction = inputs.kinetics.substrate_limit_fraction;
    equilibrium.hpo4_from_h2po4_mol_p_per_m3 =
        inputs.products.h2po4_dissociation_constant_mol_per_m3 *
        inputs.aqueous.h2po4_activity_mol_p_per_m3 /
        inputs.shared_activities.hydrogen_mol_per_m3;
    extents.h2po4_to_hydrogen_plus_hpo4_mol_p_per_m3_step =
        sourceBoundedExtent(
            inputs.aqueous.hpo4_activity_mol_p_per_m3,
            equilibrium.hpo4_from_h2po4_mol_p_per_m3,
            inputs.coefficients.divalent,
            fraction * inputs.aqueous.h2po4_concentration_mol_p_per_m3,
            fraction * inputs.aqueous.hpo4_concentration_mol_p_per_m3,
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

fn zeroResult(status: ZoneStatus) Result {
    return .{
        .status = status,
        .equilibrium_activities = std.mem.zeroes(EquilibriumActivities),
        .extents = std.mem.zeroes(Extents),
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(ZoneGeometry).@"struct".fields) |field| {
        const value = @field(inputs.geometry, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhPhosphateExchangeInput;
    }
    inline for (@typeInfo(ExchangeCapacity).@"struct".fields) |field| {
        const value = @field(inputs.capacity, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhPhosphateExchangeInput;
    }
    inline for (@typeInfo(AqueousPhosphate).@"struct".fields) |field| {
        const value = @field(inputs.aqueous, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhPhosphateExchangeInput;
    }
    inline for (@typeInfo(SurfaceSites).@"struct".fields) |field| {
        const value = @field(inputs.sites, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFixedPhPhosphateExchangeInput;
    }
    inline for (@typeInfo(SharedActivities).@"struct".fields) |field| {
        const value = @field(inputs.shared_activities, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFixedPhPhosphateExchangeInput;
    }
    inline for (@typeInfo(EquilibriumProducts).@"struct".fields) |field| {
        const value = @field(inputs.products, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFixedPhPhosphateExchangeInput;
    }
    if (!std.math.isFinite(inputs.kinetics.substrate_limit_fraction) or
        inputs.kinetics.substrate_limit_fraction < 0 or
        inputs.kinetics.substrate_limit_fraction > 1)
        return error.InvalidFixedPhPhosphateExchangeInput;
    inline for (.{
        inputs.kinetics.maximum_anion_exchange_source_extent_per_step,
        inputs.kinetics.maximum_general_reaction_mol_per_m3_step,
    }) |maximum| {
        if (!std.math.isFinite(maximum) or maximum < 0)
            return error.InvalidFixedPhPhosphateExchangeInput;
    }
    inline for (@typeInfo(ActivityCoefficients).@"struct".fields) |field| {
        const value = @field(inputs.coefficients, field.name);
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidFixedPhPhosphateExchangeInput;
    }
}

fn validateActiveSites(sites: SurfaceSites) !void {
    if (sites.hydroxyl_site_mol_per_Mg <= 0 or
        sites.protonated_site_mol_per_Mg <= 0)
        return error.InvalidActivePhosphateExchangeSite;
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(EquilibriumActivities).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.equilibrium_activities, field.name)))
            return error.NonFiniteFixedPhPhosphateExchangeResult;
    inline for (@typeInfo(Extents).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.extents, field.name)))
            return error.NonFiniteFixedPhPhosphateExchangeResult;
}

fn validInputs() Inputs {
    return .{
        .geometry = .{
            .non_band_water_volume_m3 = 0.7,
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
            .hydroxyl_site_mol_per_Mg = 0.4,
            .protonated_site_mol_per_Mg = 0.3,
            .adsorbed_hpo4_mol_p_per_Mg = 0.04,
            .adsorbed_h2po4_mol_p_per_Mg = 0.05,
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
            .substrate_limit_fraction = 0.5,
            .maximum_anion_exchange_source_extent_per_step = 1,
            .maximum_general_reaction_mol_per_m3_step = 0.2,
        },
        .coefficients = .{
            .monovalent = 0.8,
            .divalent = 0.6,
        },
    };
}

test "fixed-pH non-band exchange matches every source equation exactly" {
    const inputs = validInputs();
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(ZoneStatus.exchange_active, result.status);
    const a = inputs.aqueous;
    const s = inputs.sites;
    const shared = inputs.shared_activities;
    const p = inputs.products;
    const fraction = inputs.kinetics.substrate_limit_fraction;
    const maximum =
        inputs.kinetics.maximum_anion_exchange_source_extent_per_step;

    const protonated_product =
        p.h2po4_exchange_constant *
        shared.hydrogen_mol_per_m3 *
        shared.hydroxide_mol_per_m3;
    const protonated_equilibrium =
        protonated_product * s.adsorbed_h2po4_mol_p_per_Mg /
        s.protonated_site_mol_per_Mg;
    const hydroxyl_product =
        p.h2po4_exchange_constant * shared.hydroxide_mol_per_m3;
    const hydroxyl_equilibrium =
        hydroxyl_product * s.adsorbed_h2po4_mol_p_per_Mg /
        s.hydroxyl_site_mol_per_Mg;
    const hpo4_product =
        p.hpo4_exchange_constant *
        p.water_activity_product_mol2_per_m6 /
        p.h2po4_dissociation_constant_mol_per_m3;
    const hpo4_equilibrium =
        hpo4_product * s.adsorbed_hpo4_mol_p_per_Mg /
        s.hydroxyl_site_mol_per_Mg;
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
            fraction * s.adsorbed_h2po4_mol_p_per_Mg,
            fraction * @min(
                a.h2po4_concentration_mol_p_per_m3,
                s.protonated_site_mol_per_Mg,
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
            fraction * s.adsorbed_h2po4_mol_p_per_Mg,
            fraction * @min(
                a.h2po4_concentration_mol_p_per_m3,
                s.hydroxyl_site_mol_per_Mg,
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
            fraction * s.adsorbed_hpo4_mol_p_per_Mg,
            fraction * @min(
                a.hpo4_concentration_mol_p_per_m3,
                s.hydroxyl_site_mol_per_Mg,
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
            fraction * a.h2po4_concentration_mol_p_per_m3,
            fraction * a.hpo4_concentration_mol_p_per_m3,
            inputs.kinetics.maximum_general_reaction_mol_per_m3_step,
        ),
        result.extents.h2po4_to_hydrogen_plus_hpo4_mol_p_per_m3_step,
    );
}

test "fixed-pH capacity gate retains aqueous phosphate association" {
    var inputs = validInputs();
    inputs.capacity.anion_exchange_capacity_mol =
        inputs.capacity.minimum_active_capacity_mol;
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(ZoneStatus.aqueous_only, result.status);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.extents.h2po4_with_protonated_site_source_extent_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.extents.h2po4_with_hydroxyl_site_source_extent_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.extents.hpo4_with_hydroxyl_site_source_extent_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.05),
        result.extents.h2po4_to_hydrogen_plus_hpo4_mol_p_per_m3_step,
    );
}

test "fixed-pH dry-zone gate clears exchange and aqueous association" {
    var inputs = validInputs();
    inputs.geometry.non_band_water_volume_m3 =
        inputs.geometry.minimum_active_water_volume_m3;
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(ZoneStatus.dry, result.status);
    try std.testing.expectEqualDeep(
        std.mem.zeroes(EquilibriumActivities),
        result.equilibrium_activities,
    );
    try std.testing.expectEqualDeep(std.mem.zeroes(Extents), result.extents);
}

test "fixed-pH exchange preserves independent source desorption bounds" {
    var inputs = validInputs();
    inputs.kinetics.substrate_limit_fraction = 1;
    inputs.aqueous.h2po4_activity_mol_p_per_m3 = 0;
    inputs.aqueous.hpo4_activity_mol_p_per_m3 = 0;
    inputs.products.h2po4_exchange_constant = 1.0e8;
    inputs.products.hpo4_exchange_constant = 1.0e8;
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(
        -inputs.sites.adsorbed_h2po4_mol_p_per_Mg,
        result.extents.h2po4_with_protonated_site_source_extent_per_step,
    );
    try std.testing.expectEqual(
        -inputs.sites.adsorbed_h2po4_mol_p_per_Mg,
        result.extents.h2po4_with_hydroxyl_site_source_extent_per_step,
    );
    try std.testing.expectEqual(
        -inputs.sites.adsorbed_hpo4_mol_p_per_Mg,
        result.extents.hpo4_with_hydroxyl_site_source_extent_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.extents.h2po4_to_hydrogen_plus_hpo4_mol_p_per_m3_step,
    );
}

test "fixed-pH exchange rejects invalid active denominators" {
    var inputs = validInputs();
    inputs.sites.hydroxyl_site_mol_per_Mg = 0;
    try std.testing.expectError(
        error.InvalidActivePhosphateExchangeSite,
        calculateSourceOrder(inputs),
    );
    inputs = validInputs();
    inputs.shared_activities.hydrogen_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidFixedPhPhosphateExchangeInput,
        calculateSourceOrder(inputs),
    );
    inputs = validInputs();
    inputs.coefficients.divalent = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidFixedPhPhosphateExchangeInput,
        calculateSourceOrder(inputs),
    );
}
