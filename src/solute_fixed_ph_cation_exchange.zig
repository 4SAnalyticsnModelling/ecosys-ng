const std = @import("std");
const cation_exchange = @import("solute_cation_exchange.zig");

pub const Cations = cation_exchange.Cations;
pub const Selectivity = cation_exchange.Selectivity;

pub const Geometry = struct {
    cation_exchange_capacity_mol: f64,
    minimum_active_capacity_mol: f64,
    soil_mass_per_water_volume_megagrams_per_m3: f64,
    ammonium_non_band_fraction: f64,
    ammonium_band_fraction: f64,
};

pub const ExchangeState = struct {
    cation_exchange_capacity_mol_charge_per_megagram: f64,
    aqueous_concentration_mol_per_m3: Cations,
    aqueous_activity_mol_per_m3: Cations,
    exchange_concentration_mol_per_megagram: Cations,
};

pub const KineticControls = struct {
    /// Runtime replacement for restricted-salt `FIONN`.
    substrate_limit_fraction: f64,
    /// Runtime replacement for restricted-salt `TADC`.
    maximum_adsorption_mol_charge_per_m3_step: f64,
    /// Runtime replacement for source `ZEROC`.
    minimum_activity_mol_per_m3: f64,
};

pub const Inputs = struct {
    geometry: Geometry,
    state: ExchangeState,
    selectivity: Selectivity,
    kinetics: KineticControls,
};

pub const Status = enum {
    capacity_inactive,
    equilibrium_inactive,
    active,
};

pub const ActivityRoots = struct {
    aluminum_mol_per_m3_cuberoot: f64,
    iron_mol_per_m3_cuberoot: f64,
    calcium_mol_per_m3_squareroot: f64,
    magnesium_mol_per_m3_squareroot: f64,
};

pub const Normalization = struct {
    calcium_basis_mol_charge_per_megagram: f64,
    equilibrium_total_mol_charge_per_megagram: f64,
    current_total_mol_charge_per_megagram: f64,
    equilibrium_scale: f64,
    current_scale: f64,
};

pub const Result = struct {
    status: Status,
    activity_roots: ActivityRoots,
    normalization: Normalization,
    equilibrium_charge_mol_per_megagram: Cations,
    current_charge_mol_per_megagram: Cations,
    raw_charge_change_mol_per_megagram_step: Cations,
    ion_change_mol_per_megagram_step: Cations,
    source_substrate_scaling_m3_per_megagram_step: f64,
    maximum_charge_change_mol_per_megagram_step: f64,
    raw_charge_sum_mol_per_megagram_step: f64,
    raw_charge_magnitude_mol_per_megagram_step: f64,
};

/// Direct source-order translation of SOLUTE.F lines 3348--3530.
///
/// This pure restricted-salt comparator preserves the extensive-capacity
/// gate and uses `FIONN`/`TADC`. It does not alter production Gapon callers.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validate(inputs);
    if (inputs.geometry.cation_exchange_capacity_mol <=
        inputs.geometry.minimum_active_capacity_mol)
    {
        return zeroResult(.capacity_inactive, std.mem.zeroes(ActivityRoots));
    }

    const roots = activityRoots(inputs);
    const capacity =
        inputs.state.cation_exchange_capacity_mol_charge_per_megagram;
    if (roots.calcium_mol_per_m3_squareroot <= 0 or capacity <= 0)
        return zeroResult(.equilibrium_inactive, roots);

    const normalization_and_targets =
        try calculateNormalizedCharges(inputs, roots);
    const source_substrate_scaling =
        inputs.kinetics.substrate_limit_fraction /
        inputs.geometry.soil_mass_per_water_volume_megagrams_per_m3;
    const maximum_charge_change =
        inputs.kinetics.maximum_adsorption_mol_charge_per_m3_step /
        inputs.geometry.soil_mass_per_water_volume_megagrams_per_m3;
    var raw = calculateRawChargeChanges(
        inputs,
        normalization_and_targets.equilibrium,
        normalization_and_targets.current,
        source_substrate_scaling,
        maximum_charge_change,
    );
    const raw_sum = sumCharge(raw);
    const raw_magnitude = sumAbsoluteCharge(raw);
    if (raw_magnitude == 0)
        return error.ZeroGaponChargeRedistributionMagnitude;

    const raw_before_closure = raw;
    closeChargeAndConvertToIonMoles(&raw, raw_sum, raw_magnitude);
    const result = Result{
        .status = .active,
        .activity_roots = roots,
        .normalization = normalization_and_targets.normalization,
        .equilibrium_charge_mol_per_megagram = normalization_and_targets.equilibrium,
        .current_charge_mol_per_megagram = normalization_and_targets.current,
        .raw_charge_change_mol_per_megagram_step = raw_before_closure,
        .ion_change_mol_per_megagram_step = raw,
        .source_substrate_scaling_m3_per_megagram_step = source_substrate_scaling,
        .maximum_charge_change_mol_per_megagram_step = maximum_charge_change,
        .raw_charge_sum_mol_per_megagram_step = raw_sum,
        .raw_charge_magnitude_mol_per_megagram_step = raw_magnitude,
    };
    try validateResult(result);
    return result;
}

const NormalizedCharges = struct {
    normalization: Normalization,
    equilibrium: Cations,
    current: Cations,
};

fn calculateNormalizedCharges(
    inputs: Inputs,
    roots: ActivityRoots,
) !NormalizedCharges {
    const activity = inputs.state.aqueous_activity_mol_per_m3;
    const s = inputs.selectivity;
    const calcium_root = roots.calcium_mol_per_m3_squareroot;
    const capacity =
        inputs.state.cation_exchange_capacity_mol_charge_per_megagram;
    var denominator: f64 = 1.0;
    denominator += s.calcium_ammonium *
        activity.ammonium_non_band / calcium_root *
        inputs.geometry.ammonium_non_band_fraction;
    denominator += s.calcium_ammonium *
        activity.ammonium_band / calcium_root *
        inputs.geometry.ammonium_band_fraction;
    denominator += s.calcium_hydrogen * activity.hydrogen / calcium_root;
    denominator += s.calcium_aluminum_and_iron *
        roots.aluminum_mol_per_m3_cuberoot / calcium_root * 3.0;
    denominator += s.calcium_aluminum_and_iron *
        roots.iron_mol_per_m3_cuberoot / calcium_root * 3.0;
    denominator += s.calcium_magnesium *
        roots.magnesium_mol_per_m3_squareroot / calcium_root * 2.0;
    denominator += s.calcium_sodium * activity.sodium / calcium_root;
    denominator += s.calcium_potassium * activity.potassium / calcium_root;
    const calcium_basis = capacity / denominator;
    if (!std.math.isFinite(calcium_basis) or calcium_basis < 0)
        return error.NonFiniteFixedPhCationExchangeResult;

    var equilibrium = sourceEquilibriumCharge(
        inputs,
        roots,
        calcium_basis,
    );
    const equilibrium_total = sumCharge(equilibrium);
    const current_total = sourceCurrentChargeTotal(
        inputs.state.exchange_concentration_mol_per_megagram,
    );
    const equilibrium_scale =
        if (equilibrium_total > 0) capacity / equilibrium_total else 0;
    const current_scale =
        if (current_total > 0) capacity / current_total else 0;
    scaleCations(&equilibrium, equilibrium_scale);
    const current = sourceCurrentCharge(
        inputs.state.exchange_concentration_mol_per_megagram,
        current_scale,
    );
    try validateCations(equilibrium);
    try validateCations(current);
    return .{
        .normalization = .{
            .calcium_basis_mol_charge_per_megagram = calcium_basis,
            .equilibrium_total_mol_charge_per_megagram = equilibrium_total,
            .current_total_mol_charge_per_megagram = current_total,
            .equilibrium_scale = equilibrium_scale,
            .current_scale = current_scale,
        },
        .equilibrium = equilibrium,
        .current = current,
    };
}

fn sourceEquilibriumCharge(
    inputs: Inputs,
    roots: ActivityRoots,
    calcium_basis: f64,
) Cations {
    const activity = inputs.state.aqueous_activity_mol_per_m3;
    const s = inputs.selectivity;
    const calcium_root = roots.calcium_mol_per_m3_squareroot;
    return .{
        .ammonium_non_band = calcium_basis *
            activity.ammonium_non_band / calcium_root *
            s.calcium_ammonium *
            inputs.geometry.ammonium_non_band_fraction,
        .ammonium_band = calcium_basis *
            activity.ammonium_band / calcium_root *
            s.calcium_ammonium *
            inputs.geometry.ammonium_band_fraction,
        .hydrogen = calcium_basis * activity.hydrogen / calcium_root *
            s.calcium_hydrogen,
        .aluminum = calcium_basis *
            roots.aluminum_mol_per_m3_cuberoot / calcium_root *
            s.calcium_aluminum_and_iron * 3.0,
        .iron = calcium_basis *
            roots.iron_mol_per_m3_cuberoot / calcium_root *
            s.calcium_aluminum_and_iron * 3.0,
        .calcium = calcium_basis * calcium_root / calcium_root * 2.0,
        .magnesium = calcium_basis *
            roots.magnesium_mol_per_m3_squareroot / calcium_root *
            s.calcium_magnesium * 2.0,
        .sodium = calcium_basis * activity.sodium / calcium_root *
            s.calcium_sodium,
        .potassium = calcium_basis * activity.potassium / calcium_root *
            s.calcium_potassium,
    };
}

fn sourceCurrentChargeTotal(exchange: Cations) f64 {
    return exchange.ammonium_non_band +
        exchange.ammonium_band +
        exchange.hydrogen +
        exchange.aluminum * 3.0 +
        exchange.iron * 3.0 +
        exchange.calcium * 2.0 +
        exchange.magnesium * 2.0 +
        exchange.sodium +
        exchange.potassium;
}

fn sourceCurrentCharge(exchange: Cations, scale: f64) Cations {
    return .{
        .ammonium_non_band = scale * exchange.ammonium_non_band,
        .ammonium_band = scale * exchange.ammonium_band,
        .hydrogen = scale * exchange.hydrogen,
        .aluminum = scale * exchange.aluminum * 3.0,
        .iron = scale * exchange.iron * 3.0,
        .calcium = scale * exchange.calcium * 2.0,
        .magnesium = scale * exchange.magnesium * 2.0,
        .sodium = scale * exchange.sodium,
        .potassium = scale * exchange.potassium,
    };
}

fn calculateRawChargeChanges(
    inputs: Inputs,
    equilibrium: Cations,
    current: Cations,
    source_substrate_scaling: f64,
    maximum: f64,
) Cations {
    const exchange = inputs.state.exchange_concentration_mol_per_megagram;
    const aqueous = inputs.state.aqueous_concentration_mol_per_m3;
    return .{
        .ammonium_non_band = sourceBoundedChargeChange(
            equilibrium.ammonium_non_band - current.ammonium_non_band,
            source_substrate_scaling *
                @min(exchange.ammonium_non_band, aqueous.ammonium_non_band),
            maximum,
        ) * inputs.geometry.ammonium_non_band_fraction,
        .ammonium_band = sourceBoundedChargeChange(
            equilibrium.ammonium_band - current.ammonium_band,
            source_substrate_scaling *
                @min(exchange.ammonium_band, aqueous.ammonium_band),
            maximum,
        ) * inputs.geometry.ammonium_band_fraction,
        .hydrogen = sourceBoundedChargeChange(
            equilibrium.hydrogen - current.hydrogen,
            source_substrate_scaling *
                @min(exchange.hydrogen, aqueous.hydrogen),
            maximum,
        ),
        .aluminum = sourceBoundedChargeChange(
            equilibrium.aluminum - current.aluminum,
            source_substrate_scaling * 3.0 *
                @min(exchange.aluminum, aqueous.aluminum),
            maximum,
        ),
        .iron = sourceBoundedChargeChange(
            equilibrium.iron - current.iron,
            source_substrate_scaling * 3.0 *
                @min(exchange.iron, aqueous.iron),
            maximum,
        ),
        .calcium = sourceBoundedChargeChange(
            equilibrium.calcium - current.calcium,
            source_substrate_scaling * 2.0 *
                @min(exchange.calcium, aqueous.calcium),
            maximum,
        ),
        .magnesium = sourceBoundedChargeChange(
            equilibrium.magnesium - current.magnesium,
            source_substrate_scaling * 2.0 *
                @min(exchange.magnesium, aqueous.magnesium),
            maximum,
        ),
        .sodium = sourceBoundedChargeChange(
            equilibrium.sodium - current.sodium,
            source_substrate_scaling *
                @min(exchange.sodium, aqueous.sodium),
            maximum,
        ),
        .potassium = sourceBoundedChargeChange(
            equilibrium.potassium - current.potassium,
            source_substrate_scaling *
                @min(exchange.potassium, aqueous.potassium),
            maximum,
        ),
    };
}

fn sourceBoundedChargeChange(
    driving_change: f64,
    substrate_limit: f64,
    maximum: f64,
) f64 {
    return @max(
        -maximum,
        -substrate_limit,
        @min(maximum, substrate_limit, driving_change),
    );
}

fn activityRoots(inputs: Inputs) ActivityRoots {
    const activity = inputs.state.aqueous_activity_mol_per_m3;
    const floor = inputs.kinetics.minimum_activity_mol_per_m3;
    return .{
        .aluminum_mol_per_m3_cuberoot = std.math.pow(f64, @max(floor, activity.aluminum), 0.333),
        .iron_mol_per_m3_cuberoot = std.math.pow(f64, @max(floor, activity.iron), 0.333),
        .calcium_mol_per_m3_squareroot = @sqrt(@max(floor, activity.calcium)),
        .magnesium_mol_per_m3_squareroot = @sqrt(@max(floor, activity.magnesium)),
    };
}

fn closeChargeAndConvertToIonMoles(
    changes: *Cations,
    raw_sum: f64,
    raw_magnitude: f64,
) void {
    inline for (@typeInfo(Cations).@"struct".fields) |field|
        @field(changes.*, field.name) -=
            raw_sum * @abs(@field(changes.*, field.name)) / raw_magnitude;
    changes.aluminum /= 3.0;
    changes.iron /= 3.0;
    changes.calcium /= 2.0;
    changes.magnesium /= 2.0;
}

fn sumCharge(values: Cations) f64 {
    return values.ammonium_non_band +
        values.ammonium_band +
        values.hydrogen +
        values.aluminum +
        values.iron +
        values.calcium +
        values.magnesium +
        values.sodium +
        values.potassium;
}

fn sumAbsoluteCharge(values: Cations) f64 {
    return @abs(values.ammonium_non_band) +
        @abs(values.ammonium_band) +
        @abs(values.hydrogen) +
        @abs(values.aluminum) +
        @abs(values.iron) +
        @abs(values.calcium) +
        @abs(values.magnesium) +
        @abs(values.sodium) +
        @abs(values.potassium);
}

fn scaleCations(values: *Cations, scale: f64) void {
    inline for (@typeInfo(Cations).@"struct".fields) |field|
        @field(values.*, field.name) *= scale;
}

fn zeroResult(status: Status, roots: ActivityRoots) Result {
    return .{
        .status = status,
        .activity_roots = roots,
        .normalization = std.mem.zeroes(Normalization),
        .equilibrium_charge_mol_per_megagram = std.mem.zeroes(Cations),
        .current_charge_mol_per_megagram = std.mem.zeroes(Cations),
        .raw_charge_change_mol_per_megagram_step = std.mem.zeroes(Cations),
        .ion_change_mol_per_megagram_step = std.mem.zeroes(Cations),
        .source_substrate_scaling_m3_per_megagram_step = 0,
        .maximum_charge_change_mol_per_megagram_step = 0,
        .raw_charge_sum_mol_per_megagram_step = 0,
        .raw_charge_magnitude_mol_per_megagram_step = 0,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Geometry).@"struct".fields) |field|
        try finiteNonnegative(@field(inputs.geometry, field.name));
    if (inputs.geometry.soil_mass_per_water_volume_megagrams_per_m3 <= 0 or
        inputs.geometry.ammonium_non_band_fraction > 1 or
        inputs.geometry.ammonium_band_fraction > 1)
        return error.InvalidFixedPhCationExchangeInput;
    try finiteNonnegative(
        inputs.state.cation_exchange_capacity_mol_charge_per_megagram,
    );
    try validateCations(inputs.state.aqueous_concentration_mol_per_m3);
    try validateCations(inputs.state.aqueous_activity_mol_per_m3);
    try validateCations(inputs.state.exchange_concentration_mol_per_megagram);
    inline for (@typeInfo(Selectivity).@"struct".fields) |field|
        try finiteNonnegative(@field(inputs.selectivity, field.name));
    if (!std.math.isFinite(inputs.kinetics.substrate_limit_fraction) or
        inputs.kinetics.substrate_limit_fraction < 0 or
        inputs.kinetics.substrate_limit_fraction > 1 or
        !std.math.isFinite(
            inputs.kinetics.maximum_adsorption_mol_charge_per_m3_step,
        ) or
        inputs.kinetics.maximum_adsorption_mol_charge_per_m3_step < 0 or
        !std.math.isFinite(inputs.kinetics.minimum_activity_mol_per_m3) or
        inputs.kinetics.minimum_activity_mol_per_m3 <= 0)
        return error.InvalidFixedPhCationExchangeInput;
}

fn validateCations(values: Cations) !void {
    inline for (@typeInfo(Cations).@"struct".fields) |field|
        try finiteNonnegative(@field(values, field.name));
}

fn finiteNonnegative(value: f64) !void {
    if (!std.math.isFinite(value) or value < 0)
        return error.InvalidFixedPhCationExchangeInput;
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(ActivityRoots).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.activity_roots, field.name)))
            return error.NonFiniteFixedPhCationExchangeResult;
    inline for (@typeInfo(Normalization).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.normalization, field.name)))
            return error.NonFiniteFixedPhCationExchangeResult;
    try validateCations(result.equilibrium_charge_mol_per_megagram);
    try validateCations(result.current_charge_mol_per_megagram);
    inline for (.{
        result.raw_charge_change_mol_per_megagram_step,
        result.ion_change_mol_per_megagram_step,
    }) |values| inline for (@typeInfo(Cations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(values, field.name)))
            return error.NonFiniteFixedPhCationExchangeResult;
}

fn validInputs() Inputs {
    return .{
        .geometry = .{
            .cation_exchange_capacity_mol = 2,
            .minimum_active_capacity_mol = 0.1,
            .soil_mass_per_water_volume_megagrams_per_m3 = 1.4,
            .ammonium_non_band_fraction = 0.7,
            .ammonium_band_fraction = 0.3,
        },
        .state = .{
            .cation_exchange_capacity_mol_charge_per_megagram = 1.7,
            .aqueous_concentration_mol_per_m3 = .{
                .ammonium_non_band = 0.4,
                .ammonium_band = 0.2,
                .hydrogen = 0.03,
                .aluminum = 0.04,
                .iron = 0.02,
                .calcium = 0.5,
                .magnesium = 0.25,
                .sodium = 0.12,
                .potassium = 0.08,
            },
            .aqueous_activity_mol_per_m3 = .{
                .ammonium_non_band = 0.32,
                .ammonium_band = 0.14,
                .hydrogen = 0.02,
                .aluminum = 0.027,
                .iron = 0.008,
                .calcium = 0.36,
                .magnesium = 0.16,
                .sodium = 0.09,
                .potassium = 0.05,
            },
            .exchange_concentration_mol_per_megagram = .{
                .ammonium_non_band = 0.2,
                .ammonium_band = 0.1,
                .hydrogen = 0.03,
                .aluminum = 0.04,
                .iron = 0.02,
                .calcium = 0.25,
                .magnesium = 0.1,
                .sodium = 0.08,
                .potassium = 0.06,
            },
        },
        .selectivity = .{
            .calcium_ammonium = 1.1,
            .calcium_hydrogen = 0.9,
            .calcium_aluminum_and_iron = 1.2,
            .calcium_magnesium = 0.8,
            .calcium_sodium = 1.3,
            .calcium_potassium = 1.4,
        },
        .kinetics = .{
            .substrate_limit_fraction = 0.35,
            .maximum_adsorption_mol_charge_per_m3_step = 0.07,
            .minimum_activity_mol_per_m3 = 1e-32,
        },
    };
}

test "restricted-salt Gapon matches source roots normalization and rates" {
    const inputs = validInputs();
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(Status.active, result.status);
    const floor = inputs.kinetics.minimum_activity_mol_per_m3;
    try std.testing.expectEqual(
        std.math.pow(
            f64,
            @max(floor, inputs.state.aqueous_activity_mol_per_m3.aluminum),
            0.333,
        ),
        result.activity_roots.aluminum_mol_per_m3_cuberoot,
    );
    try std.testing.expectEqual(
        @sqrt(@max(
            floor,
            inputs.state.aqueous_activity_mol_per_m3.calcium,
        )),
        result.activity_roots.calcium_mol_per_m3_squareroot,
    );
    try std.testing.expectApproxEqAbs(
        inputs.state.cation_exchange_capacity_mol_charge_per_megagram,
        sumCharge(result.equilibrium_charge_mol_per_megagram),
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        inputs.state.cation_exchange_capacity_mol_charge_per_megagram,
        sumCharge(result.current_charge_mol_per_megagram),
        1e-15,
    );
    const final_charge =
        result.ion_change_mol_per_megagram_step.ammonium_non_band +
        result.ion_change_mol_per_megagram_step.ammonium_band +
        result.ion_change_mol_per_megagram_step.hydrogen +
        3.0 * result.ion_change_mol_per_megagram_step.aluminum +
        3.0 * result.ion_change_mol_per_megagram_step.iron +
        2.0 * result.ion_change_mol_per_megagram_step.calcium +
        2.0 * result.ion_change_mol_per_megagram_step.magnesium +
        result.ion_change_mol_per_megagram_step.sodium +
        result.ion_change_mol_per_megagram_step.potassium;
    try std.testing.expectApproxEqAbs(@as(f64, 0), final_charge, 1e-15);
}

test "restricted-salt Gapon equals earlier source block when controls match" {
    const inputs = validInputs();
    const restricted = try calculateSourceOrder(inputs);
    const earlier = try cation_exchange.calculateSourceOrder(.{
        .cation_exchange_capacity_mol_charge_per_megagram = inputs.state.cation_exchange_capacity_mol_charge_per_megagram,
        .aqueous_concentration_mol_per_m3 = inputs.state.aqueous_concentration_mol_per_m3,
        .aqueous_activity_mol_per_m3 = inputs.state.aqueous_activity_mol_per_m3,
        .exchange_concentration_mol_per_megagram = inputs.state.exchange_concentration_mol_per_megagram,
        .ammonium_non_band_fraction = inputs.geometry.ammonium_non_band_fraction,
        .ammonium_band_fraction = inputs.geometry.ammonium_band_fraction,
        .soil_mass_per_water_volume_megagrams_per_m3 = inputs.geometry.soil_mass_per_water_volume_megagrams_per_m3,
    }, .{
        .selectivity = inputs.selectivity,
        .substrate_limit_fraction = inputs.kinetics.substrate_limit_fraction,
        .maximum_adsorption_mol_charge_per_m3_step = inputs.kinetics.maximum_adsorption_mol_charge_per_m3_step,
    }, .{
        .minimum_activity_mol_per_m3 = inputs.kinetics.minimum_activity_mol_per_m3,
    });
    try std.testing.expectEqualDeep(
        earlier,
        restricted.ion_change_mol_per_megagram_step,
    );
}

test "restricted-salt Gapon preserves both source gates" {
    var inputs = validInputs();
    inputs.geometry.cation_exchange_capacity_mol =
        inputs.geometry.minimum_active_capacity_mol;
    var result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(Status.capacity_inactive, result.status);
    try std.testing.expectEqualDeep(
        std.mem.zeroes(Cations),
        result.ion_change_mol_per_megagram_step,
    );

    inputs = validInputs();
    inputs.state.cation_exchange_capacity_mol_charge_per_megagram = 0;
    result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(Status.equilibrium_inactive, result.status);
    try std.testing.expect(
        result.activity_roots.calcium_mol_per_m3_squareroot > 0,
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(Cations),
        result.ion_change_mol_per_megagram_step,
    );
}

test "restricted-salt Gapon fails before source zero-magnitude division" {
    var inputs = validInputs();
    inputs.state.aqueous_concentration_mol_per_m3 = std.mem.zeroes(Cations);
    inputs.state.aqueous_activity_mol_per_m3 = std.mem.zeroes(Cations);
    inputs.state.aqueous_activity_mol_per_m3.calcium = 1;
    inputs.state.exchange_concentration_mol_per_megagram =
        std.mem.zeroes(Cations);
    inputs.state.exchange_concentration_mol_per_megagram.calcium = 0.5;
    inputs.state.cation_exchange_capacity_mol_charge_per_megagram = 1;
    inputs.selectivity = std.mem.zeroes(Selectivity);
    try std.testing.expectError(
        error.ZeroGaponChargeRedistributionMagnitude,
        calculateSourceOrder(inputs),
    );
}

test "restricted-salt Gapon rejects invalid runtime controls" {
    var inputs = validInputs();
    inputs.kinetics.minimum_activity_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidFixedPhCationExchangeInput,
        calculateSourceOrder(inputs),
    );
    inputs = validInputs();
    inputs.geometry.ammonium_band_fraction = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidFixedPhCationExchangeInput,
        calculateSourceOrder(inputs),
    );
}
