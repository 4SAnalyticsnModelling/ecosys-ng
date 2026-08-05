const std = @import("std");

pub const Cations = struct {
    ammonium: f64,
    hydrogen: f64,
    aluminum: f64,
    iron: f64,
    calcium: f64,
    magnesium: f64,
    sodium: f64,
    potassium: f64,
};

pub const ActivityRoots = struct {
    aluminum_mol_per_m3_cuberoot: f64,
    iron_mol_per_m3_cuberoot: f64,
    calcium_mol_per_m3_squareroot: f64,
    magnesium_mol_per_m3_squareroot: f64,
};

pub const Selectivity = struct {
    calcium_ammonium: f64,
    calcium_hydrogen: f64,
    calcium_aluminum_and_iron: f64,
    calcium_magnesium: f64,
    calcium_sodium: f64,
    calcium_potassium: f64,
};

pub const Inputs = struct {
    cation_exchange_capacity_mol_charge_per_megagram: f64,
    aqueous_concentration_mol_per_m3: Cations,
    aqueous_activity_mol_per_m3: Cations,
    exchange_concentration_mol_per_megagram: Cations,
    activity_roots: ActivityRoots,
    selectivity: Selectivity,
    litter_mass_per_water_volume_megagrams_per_m3: f64,
    substrate_limit_fraction_per_step: f64,
    maximum_cation_adsorption_mol_charge_per_m3_step: f64,
    activation_threshold: f64,
    /// SOLUTE.F 4857--4861 consumes RXNBQ without assigning it in this
    /// surface branch. Requiring retained state prevents implicit storage.
    retained_band_ammonium_raw_change_mol_charge_per_megagram_step: f64,
};

pub const Normalization = struct {
    calcium_basis_mol_per_megagram: f64,
    equilibrium_total_mol_per_megagram: f64,
    current_total_mol_charge_per_megagram: f64,
    equilibrium_scale: f64,
    current_scale: f64,
};

pub const ActiveResult = struct {
    normalization: Normalization,
    equilibrium_target_mol_per_megagram: Cations,
    current_charge_mol_per_megagram: Cations,
    raw_charge_change_mol_per_megagram_step: Cations,
    ion_change_mol_per_megagram_step: Cations,
    retained_band_ammonium_raw_change_mol_charge_per_megagram_step: f64,
    discarded_band_ammonium_change_mol_per_megagram_step: f64,
    substrate_scaling_m3_per_megagram_step: f64,
    maximum_charge_change_mol_per_megagram_step: f64,
    raw_charge_sum_mol_per_megagram_step: f64,
    raw_charge_magnitude_mol_per_megagram_step: f64,
};

pub const Result = union(enum) {
    exchange_inactive,
    active: ActiveResult,
};

const NormalizedTargets = struct {
    normalization: Normalization,
    equilibrium: Cations,
    current: Cations,
};

/// Direct source-order translation of SOLUTE.F lines 4768--4891.
///
/// The caller resolves one surface-litter cell and applies returned ion-mole
/// changes. The kernel is pure and preserves the source's charge closure and
/// retained band-ammonium dependency.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const calcium_root = inputs.activity_roots
        .calcium_mol_per_m3_squareroot;
    const capacity = inputs.cation_exchange_capacity_mol_charge_per_megagram;

    // SOLUTE.F 4768 and 4882--4891.
    if (calcium_root <= inputs.activation_threshold or
        capacity <= inputs.activation_threshold)
    {
        return .exchange_inactive;
    }

    const targets = try normalizedTargets(inputs);

    // SOLUTE.F 4831--4832.
    const substrate_scaling =
        inputs.substrate_limit_fraction_per_step /
        inputs.litter_mass_per_water_volume_megagrams_per_m3;
    const maximum =
        inputs.maximum_cation_adsorption_mol_charge_per_m3_step /
        inputs.litter_mass_per_water_volume_megagrams_per_m3;
    if (!std.math.isFinite(substrate_scaling) or
        !std.math.isFinite(maximum))
    {
        return error.NonFiniteSurfaceLitterStaticCationExchangeResult;
    }

    // SOLUTE.F 4833--4856.
    var charge_change = rawChargeChanges(
        inputs,
        targets.equilibrium,
        targets.current,
        substrate_scaling,
        maximum,
    );
    const unclosed_charge_change = charge_change;
    const retained_band =
        inputs.retained_band_ammonium_raw_change_mol_charge_per_megagram_step;

    // SOLUTE.F 4857--4859 preserves the exact summation order.
    const raw_sum = sum(charge_change) + retained_band;
    const raw_magnitude = sumAbsolute(charge_change) + @abs(retained_band);
    if (!std.math.isFinite(raw_sum) or !std.math.isFinite(raw_magnitude))
        return error.NonFiniteSurfaceLitterStaticCationExchangeResult;
    if (raw_magnitude == 0)
        return error.ZeroSurfaceLitterStaticCationExchangeClosureMagnitude;

    // SOLUTE.F 4860--4868.
    const band_after_closure =
        retained_band - raw_sum * @abs(retained_band) / raw_magnitude;
    closeChargeAndConvertToIonMoles(
        &charge_change,
        raw_sum,
        raw_magnitude,
    );

    const active: ActiveResult = .{
        .normalization = targets.normalization,
        .equilibrium_target_mol_per_megagram = targets.equilibrium,
        .current_charge_mol_per_megagram = targets.current,
        .raw_charge_change_mol_per_megagram_step = unclosed_charge_change,
        .ion_change_mol_per_megagram_step = charge_change,
        .retained_band_ammonium_raw_change_mol_charge_per_megagram_step = retained_band,
        .discarded_band_ammonium_change_mol_per_megagram_step = band_after_closure,
        .substrate_scaling_m3_per_megagram_step = substrate_scaling,
        .maximum_charge_change_mol_per_megagram_step = maximum,
        .raw_charge_sum_mol_per_megagram_step = raw_sum,
        .raw_charge_magnitude_mol_per_megagram_step = raw_magnitude,
    };
    try validateActiveResult(active);
    return .{ .active = active };
}

fn normalizedTargets(inputs: Inputs) !NormalizedTargets {
    const activity = inputs.aqueous_activity_mol_per_m3;
    const roots = inputs.activity_roots;
    const selectivity = inputs.selectivity;
    const calcium_root = roots.calcium_mol_per_m3_squareroot;
    const capacity = inputs.cation_exchange_capacity_mol_charge_per_megagram;

    // SOLUTE.F 4769--4776.
    const calcium_basis = capacity /
        (1.0 +
            selectivity.calcium_ammonium * activity.ammonium / calcium_root +
            selectivity.calcium_hydrogen * activity.hydrogen / calcium_root +
            selectivity.calcium_aluminum_and_iron *
                roots.aluminum_mol_per_m3_cuberoot / calcium_root * 3.0 +
            selectivity.calcium_aluminum_and_iron *
                roots.iron_mol_per_m3_cuberoot / calcium_root * 3.0 +
            selectivity.calcium_magnesium *
                roots.magnesium_mol_per_m3_squareroot / calcium_root * 2.0 +
            selectivity.calcium_sodium * activity.sodium / calcium_root +
            selectivity.calcium_potassium * activity.potassium / calcium_root);
    if (!std.math.isFinite(calcium_basis) or calcium_basis < 0)
        return error.NonFiniteSurfaceLitterStaticCationExchangeResult;

    // SOLUTE.F 4777--4786. Multivalent equilibrium targets are unweighted,
    // while the current exchange total is explicitly charge weighted.
    var equilibrium: Cations = .{
        .ammonium = calcium_basis * activity.ammonium / calcium_root *
            selectivity.calcium_ammonium,
        .hydrogen = calcium_basis * activity.hydrogen / calcium_root *
            selectivity.calcium_hydrogen,
        .aluminum = calcium_basis *
            roots.aluminum_mol_per_m3_cuberoot / calcium_root *
            selectivity.calcium_aluminum_and_iron,
        .iron = calcium_basis *
            roots.iron_mol_per_m3_cuberoot / calcium_root *
            selectivity.calcium_aluminum_and_iron,
        .calcium = calcium_basis * calcium_root / calcium_root,
        .magnesium = calcium_basis *
            roots.magnesium_mol_per_m3_squareroot / calcium_root *
            selectivity.calcium_magnesium,
        .sodium = calcium_basis * activity.sodium / calcium_root *
            selectivity.calcium_sodium,
        .potassium = calcium_basis * activity.potassium / calcium_root *
            selectivity.calcium_potassium,
    };
    const equilibrium_total = sum(equilibrium);
    const current_total =
        currentChargeTotal(inputs.exchange_concentration_mol_per_megagram);

    // SOLUTE.F 4787--4796.
    const equilibrium_scale =
        if (equilibrium_total > 0) capacity / equilibrium_total else 0;
    const current_scale =
        if (current_total > 0) capacity / current_total else 0;

    // SOLUTE.F 4797--4812.
    scale(&equilibrium, equilibrium_scale);
    const current =
        currentCharge(inputs.exchange_concentration_mol_per_megagram, current_scale);
    try validateNonnegativeCations(equilibrium);
    try validateNonnegativeCations(current);
    return .{
        .normalization = .{
            .calcium_basis_mol_per_megagram = calcium_basis,
            .equilibrium_total_mol_per_megagram = equilibrium_total,
            .current_total_mol_charge_per_megagram = current_total,
            .equilibrium_scale = equilibrium_scale,
            .current_scale = current_scale,
        },
        .equilibrium = equilibrium,
        .current = current,
    };
}

fn rawChargeChanges(
    inputs: Inputs,
    equilibrium: Cations,
    current: Cations,
    substrate_scaling: f64,
    maximum: f64,
) Cations {
    const exchange = inputs.exchange_concentration_mol_per_megagram;
    const aqueous = inputs.aqueous_concentration_mol_per_m3;
    return .{
        .ammonium = boundedChargeChange(
            equilibrium.ammonium - current.ammonium,
            substrate_scaling * @min(exchange.ammonium, aqueous.ammonium),
            maximum,
        ),
        .hydrogen = boundedChargeChange(
            equilibrium.hydrogen - current.hydrogen,
            substrate_scaling * @min(exchange.hydrogen, aqueous.hydrogen),
            maximum,
        ),
        .aluminum = boundedChargeChange(
            equilibrium.aluminum - current.aluminum,
            substrate_scaling * 3.0 *
                @min(exchange.aluminum, aqueous.aluminum),
            maximum,
        ),
        .iron = boundedChargeChange(
            equilibrium.iron - current.iron,
            substrate_scaling * 3.0 *
                @min(exchange.iron, aqueous.iron),
            maximum,
        ),
        .calcium = boundedChargeChange(
            equilibrium.calcium - current.calcium,
            substrate_scaling * 2.0 *
                @min(exchange.calcium, aqueous.calcium),
            maximum,
        ),
        .magnesium = boundedChargeChange(
            equilibrium.magnesium - current.magnesium,
            substrate_scaling * 2.0 *
                @min(exchange.magnesium, aqueous.magnesium),
            maximum,
        ),
        .sodium = boundedChargeChange(
            equilibrium.sodium - current.sodium,
            substrate_scaling * @min(exchange.sodium, aqueous.sodium),
            maximum,
        ),
        .potassium = boundedChargeChange(
            equilibrium.potassium - current.potassium,
            substrate_scaling * @min(exchange.potassium, aqueous.potassium),
            maximum,
        ),
    };
}

fn boundedChargeChange(
    driving_mol_per_megagram_step: f64,
    substrate_limit_mol_per_megagram_step: f64,
    maximum_mol_per_megagram_step: f64,
) f64 {
    return @max(
        -maximum_mol_per_megagram_step,
        -substrate_limit_mol_per_megagram_step,
        @min(
            maximum_mol_per_megagram_step,
            substrate_limit_mol_per_megagram_step,
            driving_mol_per_megagram_step,
        ),
    );
}

fn currentChargeTotal(exchange: Cations) f64 {
    return exchange.ammonium +
        exchange.hydrogen +
        exchange.aluminum * 3.0 +
        exchange.iron * 3.0 +
        exchange.calcium * 2.0 +
        exchange.magnesium * 2.0 +
        exchange.sodium +
        exchange.potassium;
}

fn currentCharge(exchange: Cations, factor: f64) Cations {
    return .{
        .ammonium = factor * exchange.ammonium,
        .hydrogen = factor * exchange.hydrogen,
        .aluminum = factor * exchange.aluminum * 3.0,
        .iron = factor * exchange.iron * 3.0,
        .calcium = factor * exchange.calcium * 2.0,
        .magnesium = factor * exchange.magnesium * 2.0,
        .sodium = factor * exchange.sodium,
        .potassium = factor * exchange.potassium,
    };
}

fn closeChargeAndConvertToIonMoles(
    changes: *Cations,
    raw_sum: f64,
    raw_magnitude: f64,
) void {
    inline for (@typeInfo(Cations).@"struct".fields) |field| {
        @field(changes.*, field.name) -=
            raw_sum * @abs(@field(changes.*, field.name)) / raw_magnitude;
    }
    changes.aluminum /= 3.0;
    changes.iron /= 3.0;
    changes.calcium /= 2.0;
    changes.magnesium /= 2.0;
}

fn sum(values: Cations) f64 {
    return values.ammonium +
        values.hydrogen +
        values.aluminum +
        values.iron +
        values.calcium +
        values.magnesium +
        values.sodium +
        values.potassium;
}

fn sumAbsolute(values: Cations) f64 {
    return @abs(values.ammonium) +
        @abs(values.hydrogen) +
        @abs(values.aluminum) +
        @abs(values.iron) +
        @abs(values.calcium) +
        @abs(values.magnesium) +
        @abs(values.sodium) +
        @abs(values.potassium);
}

fn scale(values: *Cations, factor: f64) void {
    inline for (@typeInfo(Cations).@"struct".fields) |field|
        @field(values.*, field.name) *= factor;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Cations).@"struct".fields) |field| {
        inline for (.{
            @field(inputs.aqueous_concentration_mol_per_m3, field.name),
            @field(inputs.aqueous_activity_mol_per_m3, field.name),
            @field(inputs.exchange_concentration_mol_per_megagram, field.name),
        }) |value| {
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidSurfaceLitterStaticCationExchangeInput;
        }
    }
    try validateNonnegativeStruct(inputs.activity_roots);
    try validateNonnegativeStruct(inputs.selectivity);
    inline for (.{
        inputs.cation_exchange_capacity_mol_charge_per_megagram,
        inputs.substrate_limit_fraction_per_step,
        inputs.maximum_cation_adsorption_mol_charge_per_m3_step,
        inputs.activation_threshold,
    }) |value| {
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterStaticCationExchangeInput;
    }
    if (!std.math.isFinite(inputs.litter_mass_per_water_volume_megagrams_per_m3) or
        inputs.litter_mass_per_water_volume_megagrams_per_m3 <= 0 or
        inputs.substrate_limit_fraction_per_step > 1 or
        !std.math.isFinite(
            inputs
                .retained_band_ammonium_raw_change_mol_charge_per_megagram_step,
        ))
    {
        return error.InvalidSurfaceLitterStaticCationExchangeInput;
    }
}

fn validateNonnegativeStruct(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| {
        const field_value = @field(value, field.name);
        if (!std.math.isFinite(field_value) or field_value < 0)
            return error.InvalidSurfaceLitterStaticCationExchangeInput;
    }
}

fn validateNonnegativeCations(values: Cations) !void {
    inline for (@typeInfo(Cations).@"struct".fields) |field| {
        const value = @field(values, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.NonFiniteSurfaceLitterStaticCationExchangeResult;
    }
}

fn validateActiveResult(result: ActiveResult) !void {
    try validateNonnegativeCations(result.equilibrium_target_mol_per_megagram);
    try validateNonnegativeCations(result.current_charge_mol_per_megagram);
    inline for (@typeInfo(Cations).@"struct".fields) |field| {
        inline for (.{
            @field(result.raw_charge_change_mol_per_megagram_step, field.name),
            @field(result.ion_change_mol_per_megagram_step, field.name),
        }) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteSurfaceLitterStaticCationExchangeResult;
        }
    }
    try validateNonnegativeStruct(result.normalization);
    inline for (.{
        result.retained_band_ammonium_raw_change_mol_charge_per_megagram_step,
        result.discarded_band_ammonium_change_mol_per_megagram_step,
        result.substrate_scaling_m3_per_megagram_step,
        result.maximum_charge_change_mol_per_megagram_step,
        result.raw_charge_sum_mol_per_megagram_step,
        result.raw_charge_magnitude_mol_per_megagram_step,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterStaticCationExchangeResult;
    }
}

fn testInputs() Inputs {
    return .{
        .cation_exchange_capacity_mol_charge_per_megagram = 2,
        .aqueous_concentration_mol_per_m3 = .{
            .ammonium = 0.7,
            .hydrogen = 0.4,
            .aluminum = 0.3,
            .iron = 0.2,
            .calcium = 0.8,
            .magnesium = 0.5,
            .sodium = 0.6,
            .potassium = 0.25,
        },
        .aqueous_activity_mol_per_m3 = .{
            .ammonium = 0.63,
            .hydrogen = 0.36,
            .aluminum = 0.21,
            .iron = 0.14,
            .calcium = 0.56,
            .magnesium = 0.35,
            .sodium = 0.54,
            .potassium = 0.225,
        },
        .exchange_concentration_mol_per_megagram = .{
            .ammonium = 0.3,
            .hydrogen = 0.2,
            .aluminum = 0.04,
            .iron = 0.03,
            .calcium = 0.4,
            .magnesium = 0.2,
            .sodium = 0.15,
            .potassium = 0.1,
        },
        .activity_roots = .{
            .aluminum_mol_per_m3_cuberoot = 0.8,
            .iron_mol_per_m3_cuberoot = 0.7,
            .calcium_mol_per_m3_squareroot = 0.75,
            .magnesium_mol_per_m3_squareroot = 0.6,
        },
        .selectivity = .{
            .calcium_ammonium = 1.1,
            .calcium_hydrogen = 0.9,
            .calcium_aluminum_and_iron = 1.2,
            .calcium_magnesium = 0.8,
            .calcium_sodium = 1.3,
            .calcium_potassium = 1.4,
        },
        .litter_mass_per_water_volume_megagrams_per_m3 = 1.4,
        .substrate_limit_fraction_per_step = 0.2,
        .maximum_cation_adsorption_mol_charge_per_m3_step = 0.1,
        .activation_threshold = 1.0e-15,
        .retained_band_ammonium_raw_change_mol_charge_per_megagram_step = 0,
    };
}

test "SOLUTE static Gapon equilibrium preserves source expressions" {
    const inputs = testInputs();
    const active = (try calculateSourceOrder(inputs)).active;
    const activity = inputs.aqueous_activity_mol_per_m3;
    const roots = inputs.activity_roots;
    const selectivity = inputs.selectivity;
    const calcium_root = roots.calcium_mol_per_m3_squareroot;
    const basis = inputs.cation_exchange_capacity_mol_charge_per_megagram /
        (1.0 +
            selectivity.calcium_ammonium * activity.ammonium / calcium_root +
            selectivity.calcium_hydrogen * activity.hydrogen / calcium_root +
            selectivity.calcium_aluminum_and_iron *
                roots.aluminum_mol_per_m3_cuberoot / calcium_root * 3.0 +
            selectivity.calcium_aluminum_and_iron *
                roots.iron_mol_per_m3_cuberoot / calcium_root * 3.0 +
            selectivity.calcium_magnesium *
                roots.magnesium_mol_per_m3_squareroot / calcium_root * 2.0 +
            selectivity.calcium_sodium * activity.sodium / calcium_root +
            selectivity.calcium_potassium * activity.potassium / calcium_root);
    const raw_aluminum = basis *
        roots.aluminum_mol_per_m3_cuberoot / calcium_root *
        selectivity.calcium_aluminum_and_iron;

    try std.testing.expectEqual(
        basis,
        active.normalization.calcium_basis_mol_per_megagram,
    );
    try std.testing.expectEqual(
        active.normalization.equilibrium_scale * raw_aluminum,
        active.equilibrium_target_mol_per_megagram.aluminum,
    );
    try std.testing.expectApproxEqAbs(
        inputs.cation_exchange_capacity_mol_charge_per_megagram,
        sum(active.equilibrium_target_mol_per_megagram),
        1.0e-14,
    );
    try std.testing.expectApproxEqAbs(
        inputs.cation_exchange_capacity_mol_charge_per_megagram,
        sum(active.current_charge_mol_per_megagram),
        1.0e-14,
    );
}

test "static Gapon rates preserve nested bounds and charge closure" {
    const inputs = testInputs();
    const active = (try calculateSourceOrder(inputs)).active;
    const scaling = inputs.substrate_limit_fraction_per_step /
        inputs.litter_mass_per_water_volume_megagrams_per_m3;
    const maximum =
        inputs.maximum_cation_adsorption_mol_charge_per_m3_step /
        inputs.litter_mass_per_water_volume_megagrams_per_m3;
    const aluminum_limit = scaling * 3.0 *
        @min(
            inputs.exchange_concentration_mol_per_megagram.aluminum,
            inputs.aqueous_concentration_mol_per_m3.aluminum,
        );
    const expected_raw_aluminum = @max(
        -maximum,
        -aluminum_limit,
        @min(
            maximum,
            aluminum_limit,
            active.equilibrium_target_mol_per_megagram.aluminum -
                active.current_charge_mol_per_megagram.aluminum,
        ),
    );
    try std.testing.expectEqual(
        expected_raw_aluminum,
        active.raw_charge_change_mol_per_megagram_step.aluminum,
    );

    const change = active.ion_change_mol_per_megagram_step;
    const closed_charge = change.ammonium +
        change.hydrogen +
        3.0 * change.aluminum +
        3.0 * change.iron +
        2.0 * change.calcium +
        2.0 * change.magnesium +
        change.sodium +
        change.potassium +
        active.discarded_band_ammonium_change_mol_per_megagram_step;
    try std.testing.expectApproxEqAbs(@as(f64, 0), closed_charge, 1.0e-14);
}

test "static source closure exposes retained band ammonium state" {
    var inputs = testInputs();
    inputs.retained_band_ammonium_raw_change_mol_charge_per_megagram_step = 0.05;
    const with_band = (try calculateSourceOrder(inputs)).active;
    inputs.retained_band_ammonium_raw_change_mol_charge_per_megagram_step = 0;
    const without_band = (try calculateSourceOrder(inputs)).active;

    try std.testing.expect(
        with_band.ion_change_mol_per_megagram_step.calcium !=
            without_band.ion_change_mol_per_megagram_step.calcium,
    );
    try std.testing.expectEqual(
        @as(f64, 0.05),
        with_band
            .retained_band_ammonium_raw_change_mol_charge_per_megagram_step,
    );
}

test "static Gapon inactive and zero closure states are explicit" {
    var inputs = testInputs();
    inputs.cation_exchange_capacity_mol_charge_per_megagram =
        inputs.activation_threshold;
    try std.testing.expectEqual(
        Result.exchange_inactive,
        try calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.activity_roots.calcium_mol_per_m3_squareroot =
        inputs.activation_threshold;
    try std.testing.expectEqual(
        Result.exchange_inactive,
        try calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.aqueous_concentration_mol_per_m3 = std.mem.zeroes(Cations);
    inputs.exchange_concentration_mol_per_megagram = std.mem.zeroes(Cations);
    try std.testing.expectError(
        error.ZeroSurfaceLitterStaticCationExchangeClosureMagnitude,
        calculateSourceOrder(inputs),
    );
}

test "static Gapon exchange rejects invalid input and overflow" {
    var inputs = testInputs();
    inputs.exchange_concentration_mol_per_megagram.calcium = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticCationExchangeInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.substrate_limit_fraction_per_step = 1.1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticCationExchangeInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.retained_band_ammonium_raw_change_mol_charge_per_megagram_step =
        std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterStaticCationExchangeInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.litter_mass_per_water_volume_megagrams_per_m3 =
        std.math.floatMin(f64);
    inputs.maximum_cation_adsorption_mol_charge_per_m3_step =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterStaticCationExchangeResult,
        calculateSourceOrder(inputs),
    );
}
