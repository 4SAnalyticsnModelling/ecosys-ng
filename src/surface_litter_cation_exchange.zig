const std = @import("std");
const initialization = @import("surface_litter_cation_exchange_initialization.zig");

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

pub const Selectivity = struct {
    calcium_ammonium: f64,
    calcium_hydrogen: f64,
    calcium_aluminum_and_iron: f64,
    calcium_magnesium: f64,
    calcium_sodium: f64,
    calcium_potassium: f64,
};

pub const Inputs = struct {
    cation_exchange_capacity_mol_charge_per_Mg: f64,
    aqueous_concentration_mol_per_m3: Cations,
    aqueous_activity_mol_per_m3: Cations,
    exchange_concentration_mol_per_Mg: Cations,
    activity_roots: initialization.MultivalentActivityRoots,
    selectivity: Selectivity,
    litter_mass_per_water_volume_Mg_per_m3: f64,
    substrate_limit_fraction_per_step: f64,
    maximum_cation_adsorption_mol_charge_per_m3_step: f64,
    /// SOLUTE surface lines 4449--4453 read `RXNBQ` without assigning it in
    /// the surface block. The comparator requires the retained soil/band value
    /// explicitly rather than reproducing implicit or uninitialized storage.
    retained_band_ammonium_raw_charge_change_mol_per_Mg_step: f64,
};

pub const Normalization = struct {
    calcium_basis_mol_per_Mg: f64,
    equilibrium_total_mol_per_Mg: f64,
    current_total_mol_charge_per_Mg: f64,
    equilibrium_scale: f64,
    current_scale: f64,
};

pub const ActiveResult = struct {
    normalization: Normalization,
    /// The source treats these unweighted ion targets as charge during rates.
    equilibrium_target_mol_per_Mg: Cations,
    current_charge_mol_per_Mg: Cations,
    raw_charge_change_mol_per_Mg_step: Cations,
    retained_band_ammonium_raw_charge_change_mol_per_Mg_step: f64,
    ion_change_mol_per_Mg_step: Cations,
    discarded_band_ammonium_change_mol_per_Mg_step: f64,
    substrate_scaling_m3_per_Mg_step: f64,
    maximum_charge_change_mol_per_Mg_step: f64,
    raw_charge_sum_mol_per_Mg_step: f64,
    raw_charge_magnitude_mol_per_Mg_step: f64,
};

pub const Result = union(enum) {
    equilibrium_inactive,
    active: ActiveResult,
};

const NormalizedTargets = struct {
    normalization: Normalization,
    equilibrium: Cations,
    current: Cations,
};

/// Direct source-order translation of SOLUTE.F lines 4360--4481.
///
/// This surface formulation is intentionally separate from the earlier soil
/// Gapon kernel: the legacy surface equilibrium omits multivalent charge
/// factors in `XALQ..XMGQ`, and its final closure includes a retained band-NH4
/// rate that is not otherwise part of the surface system.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validateInputs(inputs);
    const calcium_root = inputs.activity_roots.calcium_mol_per_m3_squareroot;
    const capacity = inputs.cation_exchange_capacity_mol_charge_per_Mg;
    if (calcium_root <= 0 or capacity <= 0)
        return .equilibrium_inactive;

    const targets = try normalizedTargets(inputs);
    const substrate_scaling =
        inputs.substrate_limit_fraction_per_step /
        inputs.litter_mass_per_water_volume_Mg_per_m3;
    const maximum =
        inputs.maximum_cation_adsorption_mol_charge_per_m3_step /
        inputs.litter_mass_per_water_volume_Mg_per_m3;
    if (!std.math.isFinite(substrate_scaling) or
        !std.math.isFinite(maximum))
    {
        return error.NonFiniteSurfaceLitterCationExchangeResult;
    }

    var raw = rawChargeChanges(
        inputs,
        targets.equilibrium,
        targets.current,
        substrate_scaling,
        maximum,
    );
    const retained_band =
        inputs.retained_band_ammonium_raw_charge_change_mol_per_Mg_step;
    const raw_sum = sum(raw) + retained_band;
    const raw_magnitude = sumAbsolute(raw) + @abs(retained_band);
    if (!std.math.isFinite(raw_sum) or !std.math.isFinite(raw_magnitude))
        return error.NonFiniteSurfaceLitterCationExchangeResult;
    if (raw_magnitude == 0)
        return error.ZeroSurfaceLitterCationExchangeClosureMagnitude;

    const band_after_closure =
        retained_band - raw_sum * @abs(retained_band) / raw_magnitude;
    closeChargeAndConvertToIonMoles(&raw, raw_sum, raw_magnitude);
    const active: ActiveResult = .{
        .normalization = targets.normalization,
        .equilibrium_target_mol_per_Mg = targets.equilibrium,
        .current_charge_mol_per_Mg = targets.current,
        .raw_charge_change_mol_per_Mg_step = rawChargeChanges(
            inputs,
            targets.equilibrium,
            targets.current,
            substrate_scaling,
            maximum,
        ),
        .retained_band_ammonium_raw_charge_change_mol_per_Mg_step = retained_band,
        .ion_change_mol_per_Mg_step = raw,
        .discarded_band_ammonium_change_mol_per_Mg_step = band_after_closure,
        .substrate_scaling_m3_per_Mg_step = substrate_scaling,
        .maximum_charge_change_mol_per_Mg_step = maximum,
        .raw_charge_sum_mol_per_Mg_step = raw_sum,
        .raw_charge_magnitude_mol_per_Mg_step = raw_magnitude,
    };
    try validateActiveResult(active);
    return .{ .active = active };
}

fn normalizedTargets(inputs: Inputs) !NormalizedTargets {
    const activity = inputs.aqueous_activity_mol_per_m3;
    const roots = inputs.activity_roots;
    const selectivity = inputs.selectivity;
    const calcium_root = roots.calcium_mol_per_m3_squareroot;
    const capacity = inputs.cation_exchange_capacity_mol_charge_per_Mg;

    // SOLUTE.F 4361--4368.
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
        return error.NonFiniteSurfaceLitterCationExchangeResult;

    // SOLUTE.F 4369--4378. Unlike the soil block, the surface source omits
    // 3, 2, and 2 from the Al, Fe, Ca, and Mg equilibrium target terms.
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
        currentChargeTotal(inputs.exchange_concentration_mol_per_Mg);
    const equilibrium_scale =
        if (equilibrium_total > 0) capacity / equilibrium_total else 0;
    const current_scale =
        if (current_total > 0) capacity / current_total else 0;
    scale(&equilibrium, equilibrium_scale);
    const current =
        currentCharge(inputs.exchange_concentration_mol_per_Mg, current_scale);
    try validateCations(equilibrium);
    try validateCations(current);
    return .{
        .normalization = .{
            .calcium_basis_mol_per_Mg = calcium_basis,
            .equilibrium_total_mol_per_Mg = equilibrium_total,
            .current_total_mol_charge_per_Mg = current_total,
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
    const exchange = inputs.exchange_concentration_mol_per_Mg;
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

fn boundedChargeChange(driving: f64, substrate_limit: f64, maximum: f64) f64 {
    // SOLUTE.F 4425--4448 nested AMAX1/AMIN1 order.
    return @max(-maximum, -substrate_limit, @min(maximum, substrate_limit, driving));
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

fn currentCharge(exchange: Cations, current_scale: f64) Cations {
    return .{
        .ammonium = current_scale * exchange.ammonium,
        .hydrogen = current_scale * exchange.hydrogen,
        .aluminum = current_scale * exchange.aluminum * 3.0,
        .iron = current_scale * exchange.iron * 3.0,
        .calcium = current_scale * exchange.calcium * 2.0,
        .magnesium = current_scale * exchange.magnesium * 2.0,
        .sodium = current_scale * exchange.sodium,
        .potassium = current_scale * exchange.potassium,
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
            @field(inputs.exchange_concentration_mol_per_Mg, field.name),
        }) |value| {
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidSurfaceLitterCationExchangeInput;
        }
    }
    inline for (@typeInfo(initialization.MultivalentActivityRoots).@"struct".fields) |field| {
        const value = @field(inputs.activity_roots, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterCationExchangeInput;
    }
    inline for (@typeInfo(Selectivity).@"struct".fields) |field| {
        const value = @field(inputs.selectivity, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterCationExchangeInput;
    }
    if (!std.math.isFinite(inputs.cation_exchange_capacity_mol_charge_per_Mg) or
        inputs.cation_exchange_capacity_mol_charge_per_Mg < 0 or
        !std.math.isFinite(inputs.litter_mass_per_water_volume_Mg_per_m3) or
        inputs.litter_mass_per_water_volume_Mg_per_m3 <= 0 or
        !std.math.isFinite(inputs.substrate_limit_fraction_per_step) or
        inputs.substrate_limit_fraction_per_step < 0 or
        inputs.substrate_limit_fraction_per_step > 1 or
        !std.math.isFinite(inputs.maximum_cation_adsorption_mol_charge_per_m3_step) or
        inputs.maximum_cation_adsorption_mol_charge_per_m3_step < 0 or
        !std.math.isFinite(inputs.retained_band_ammonium_raw_charge_change_mol_per_Mg_step))
    {
        return error.InvalidSurfaceLitterCationExchangeInput;
    }
}

fn validateCations(values: Cations) !void {
    inline for (@typeInfo(Cations).@"struct".fields) |field| {
        const value = @field(values, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.NonFiniteSurfaceLitterCationExchangeResult;
    }
}

fn validateActiveResult(result: ActiveResult) !void {
    try validateCations(result.equilibrium_target_mol_per_Mg);
    try validateCations(result.current_charge_mol_per_Mg);
    inline for (@typeInfo(Cations).@"struct".fields) |field| {
        inline for (.{
            @field(result.raw_charge_change_mol_per_Mg_step, field.name),
            @field(result.ion_change_mol_per_Mg_step, field.name),
        }) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteSurfaceLitterCationExchangeResult;
        }
    }
    inline for (@typeInfo(Normalization).@"struct".fields) |field| {
        const value = @field(result.normalization, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.NonFiniteSurfaceLitterCationExchangeResult;
    }
    inline for (.{
        result.retained_band_ammonium_raw_charge_change_mol_per_Mg_step,
        result.discarded_band_ammonium_change_mol_per_Mg_step,
        result.substrate_scaling_m3_per_Mg_step,
        result.maximum_charge_change_mol_per_Mg_step,
        result.raw_charge_sum_mol_per_Mg_step,
        result.raw_charge_magnitude_mol_per_Mg_step,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterCationExchangeResult;
    }
}

fn testInputs() Inputs {
    return .{
        .cation_exchange_capacity_mol_charge_per_Mg = 2,
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
        .exchange_concentration_mol_per_Mg = .{
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
        .litter_mass_per_water_volume_Mg_per_m3 = 1.4,
        .substrate_limit_fraction_per_step = 0.2,
        .maximum_cation_adsorption_mol_charge_per_m3_step = 0.1,
        .retained_band_ammonium_raw_charge_change_mol_per_Mg_step = 0,
    };
}

test "surface Gapon equilibrium preserves source unweighted multivalent targets" {
    const inputs = testInputs();
    const active = (try calculateSourceOrder(inputs)).active;
    const activity = inputs.aqueous_activity_mol_per_m3;
    const roots = inputs.activity_roots;
    const s = inputs.selectivity;
    const calcium_root = roots.calcium_mol_per_m3_squareroot;
    const basis = inputs.cation_exchange_capacity_mol_charge_per_Mg /
        (1.0 +
            s.calcium_ammonium * activity.ammonium / calcium_root +
            s.calcium_hydrogen * activity.hydrogen / calcium_root +
            s.calcium_aluminum_and_iron *
                roots.aluminum_mol_per_m3_cuberoot / calcium_root * 3.0 +
            s.calcium_aluminum_and_iron *
                roots.iron_mol_per_m3_cuberoot / calcium_root * 3.0 +
            s.calcium_magnesium *
                roots.magnesium_mol_per_m3_squareroot / calcium_root * 2.0 +
            s.calcium_sodium * activity.sodium / calcium_root +
            s.calcium_potassium * activity.potassium / calcium_root);
    const raw_aluminum = basis * roots.aluminum_mol_per_m3_cuberoot /
        calcium_root * s.calcium_aluminum_and_iron;
    const raw_calcium = basis * calcium_root / calcium_root;
    const scale_factor = active.normalization.equilibrium_scale;

    try std.testing.expectEqual(
        basis,
        active.normalization.calcium_basis_mol_per_Mg,
    );
    try std.testing.expectEqual(
        scale_factor * raw_aluminum,
        active.equilibrium_target_mol_per_Mg.aluminum,
    );
    try std.testing.expectEqual(
        scale_factor * raw_calcium,
        active.equilibrium_target_mol_per_Mg.calcium,
    );
    try std.testing.expectApproxEqAbs(
        inputs.cation_exchange_capacity_mol_charge_per_Mg,
        sum(active.equilibrium_target_mol_per_Mg),
        1.0e-14,
    );
    try std.testing.expectApproxEqAbs(
        inputs.cation_exchange_capacity_mol_charge_per_Mg,
        sum(active.current_charge_mol_per_Mg),
        1.0e-14,
    );
}

test "surface and soil Gapon equations differ only by controlled source formulation" {
    const soil_exchange = @import("solute_cation_exchange.zig");
    var inputs = testInputs();
    const activity = inputs.aqueous_activity_mol_per_m3;
    inputs.activity_roots = .{
        .aluminum_mol_per_m3_cuberoot = std.math.pow(f64, activity.aluminum, 0.333),
        .iron_mol_per_m3_cuberoot = std.math.pow(f64, activity.iron, 0.333),
        .calcium_mol_per_m3_squareroot = std.math.pow(f64, activity.calcium, 0.500),
        .magnesium_mol_per_m3_squareroot = std.math.pow(f64, activity.magnesium, 0.500),
    };
    const surface = (try calculateSourceOrder(inputs)).active;
    const soil = try soil_exchange.calculateSourceOrder(.{
        .cation_exchange_capacity_mol_charge_per_Mg = inputs.cation_exchange_capacity_mol_charge_per_Mg,
        .aqueous_concentration_mol_per_m3 = .{
            .ammonium_non_band = inputs.aqueous_concentration_mol_per_m3.ammonium,
            .ammonium_band = 0,
            .hydrogen = inputs.aqueous_concentration_mol_per_m3.hydrogen,
            .aluminum = inputs.aqueous_concentration_mol_per_m3.aluminum,
            .iron = inputs.aqueous_concentration_mol_per_m3.iron,
            .calcium = inputs.aqueous_concentration_mol_per_m3.calcium,
            .magnesium = inputs.aqueous_concentration_mol_per_m3.magnesium,
            .sodium = inputs.aqueous_concentration_mol_per_m3.sodium,
            .potassium = inputs.aqueous_concentration_mol_per_m3.potassium,
        },
        .aqueous_activity_mol_per_m3 = .{
            .ammonium_non_band = activity.ammonium,
            .ammonium_band = 0,
            .hydrogen = activity.hydrogen,
            .aluminum = activity.aluminum,
            .iron = activity.iron,
            .calcium = activity.calcium,
            .magnesium = activity.magnesium,
            .sodium = activity.sodium,
            .potassium = activity.potassium,
        },
        .exchange_concentration_mol_per_Mg = .{
            .ammonium_non_band = inputs.exchange_concentration_mol_per_Mg.ammonium,
            .ammonium_band = 0,
            .hydrogen = inputs.exchange_concentration_mol_per_Mg.hydrogen,
            .aluminum = inputs.exchange_concentration_mol_per_Mg.aluminum,
            .iron = inputs.exchange_concentration_mol_per_Mg.iron,
            .calcium = inputs.exchange_concentration_mol_per_Mg.calcium,
            .magnesium = inputs.exchange_concentration_mol_per_Mg.magnesium,
            .sodium = inputs.exchange_concentration_mol_per_Mg.sodium,
            .potassium = inputs.exchange_concentration_mol_per_Mg.potassium,
        },
        .ammonium_non_band_fraction = 1,
        .ammonium_band_fraction = 0,
        .soil_mass_per_water_volume_Mg_per_m3 = inputs.litter_mass_per_water_volume_Mg_per_m3,
    }, .{
        .selectivity = .{
            .calcium_ammonium = inputs.selectivity.calcium_ammonium,
            .calcium_hydrogen = inputs.selectivity.calcium_hydrogen,
            .calcium_aluminum_and_iron = inputs.selectivity.calcium_aluminum_and_iron,
            .calcium_magnesium = inputs.selectivity.calcium_magnesium,
            .calcium_sodium = inputs.selectivity.calcium_sodium,
            .calcium_potassium = inputs.selectivity.calcium_potassium,
        },
        .substrate_limit_fraction = inputs.substrate_limit_fraction_per_step,
        .maximum_adsorption_mol_charge_per_m3_step = inputs.maximum_cation_adsorption_mol_charge_per_m3_step,
    }, .{ .minimum_activity_mol_per_m3 = 1.0e-32 });

    try std.testing.expect(
        @abs(surface.ion_change_mol_per_Mg_step.ammonium -
            soil.ammonium_non_band) > 1.0e-3,
    );
    try std.testing.expect(
        @abs(surface.equilibrium_target_mol_per_Mg.aluminum -
            3.0 * surface.normalization.equilibrium_scale *
                surface.normalization.calcium_basis_mol_per_Mg *
                inputs.activity_roots.aluminum_mol_per_m3_cuberoot /
                inputs.activity_roots.calcium_mol_per_m3_squareroot *
                inputs.selectivity.calcium_aluminum_and_iron) > 1.0e-3,
    );
}

test "surface Gapon rates preserve nested limits and charge closure" {
    const inputs = testInputs();
    const active = (try calculateSourceOrder(inputs)).active;
    const scaling = inputs.substrate_limit_fraction_per_step /
        inputs.litter_mass_per_water_volume_Mg_per_m3;
    const maximum =
        inputs.maximum_cation_adsorption_mol_charge_per_m3_step /
        inputs.litter_mass_per_water_volume_Mg_per_m3;
    const aluminum_limit = scaling * 3.0 *
        @min(
            inputs.exchange_concentration_mol_per_Mg.aluminum,
            inputs.aqueous_concentration_mol_per_m3.aluminum,
        );
    const expected_raw_aluminum = @max(
        -maximum,
        -aluminum_limit,
        @min(
            maximum,
            aluminum_limit,
            active.equilibrium_target_mol_per_Mg.aluminum -
                active.current_charge_mol_per_Mg.aluminum,
        ),
    );
    try std.testing.expectEqual(
        expected_raw_aluminum,
        active.raw_charge_change_mol_per_Mg_step.aluminum,
    );

    const change = active.ion_change_mol_per_Mg_step;
    const charge = change.ammonium +
        change.hydrogen +
        3.0 * change.aluminum +
        3.0 * change.iron +
        2.0 * change.calcium +
        2.0 * change.magnesium +
        change.sodium +
        change.potassium +
        active.discarded_band_ammonium_change_mol_per_Mg_step;
    try std.testing.expectApproxEqAbs(@as(f64, 0), charge, 1.0e-14);
}

test "surface source closure exposes retained subsurface band ammonium" {
    var inputs = testInputs();
    inputs.retained_band_ammonium_raw_charge_change_mol_per_Mg_step = 0.05;
    const with_band = (try calculateSourceOrder(inputs)).active;
    inputs.retained_band_ammonium_raw_charge_change_mol_per_Mg_step = 0;
    const without_band = (try calculateSourceOrder(inputs)).active;

    try std.testing.expect(
        with_band.ion_change_mol_per_Mg_step.calcium !=
            without_band.ion_change_mol_per_Mg_step.calcium,
    );
    const surface = with_band.ion_change_mol_per_Mg_step;
    const surface_charge = surface.ammonium +
        surface.hydrogen +
        3.0 * surface.aluminum +
        3.0 * surface.iron +
        2.0 * surface.calcium +
        2.0 * surface.magnesium +
        surface.sodium +
        surface.potassium;
    try std.testing.expectApproxEqAbs(
        -with_band.discarded_band_ammonium_change_mol_per_Mg_step,
        surface_charge,
        1.0e-14,
    );
}

test "surface Gapon inactive and zero-magnitude states are explicit" {
    var inputs = testInputs();
    inputs.cation_exchange_capacity_mol_charge_per_Mg = 0;
    try std.testing.expectEqual(
        Result.equilibrium_inactive,
        try calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.aqueous_concentration_mol_per_m3 = std.mem.zeroes(Cations);
    inputs.exchange_concentration_mol_per_Mg = std.mem.zeroes(Cations);
    inputs.retained_band_ammonium_raw_charge_change_mol_per_Mg_step = 0;
    try std.testing.expectError(
        error.ZeroSurfaceLitterCationExchangeClosureMagnitude,
        calculateSourceOrder(inputs),
    );
}

test "surface Gapon exchange rejects invalid input and derived overflow" {
    var inputs = testInputs();
    inputs.exchange_concentration_mol_per_Mg.calcium = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterCationExchangeInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.retained_band_ammonium_raw_charge_change_mol_per_Mg_step =
        std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterCationExchangeInput,
        calculateSourceOrder(inputs),
    );

    inputs = testInputs();
    inputs.litter_mass_per_water_volume_Mg_per_m3 =
        std.math.floatMin(f64);
    inputs.maximum_cation_adsorption_mol_charge_per_m3_step =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterCationExchangeResult,
        calculateSourceOrder(inputs),
    );
}
