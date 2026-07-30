const std = @import("std");

pub const Geometry = struct {
    non_band_water_volume_m3: f64,
    band_water_volume_m3: f64,
    minimum_active_water_volume_m3: f64,
};

pub const AqueousState = struct {
    ammonium_non_band_concentration_mol_n_per_m3: f64,
    ammonia_non_band_concentration_mol_n_per_m3: f64,
    ammonium_band_concentration_mol_n_per_m3: f64,
    ammonia_band_concentration_mol_n_per_m3: f64,
    ammonium_non_band_activity_mol_n_per_m3: f64,
    ammonia_non_band_activity_mol_n_per_m3: f64,
    ammonium_band_activity_mol_n_per_m3: f64,
    ammonia_band_activity_mol_n_per_m3: f64,
    hydrogen_activity_mol_per_m3: f64,
};

pub const KineticControls = struct {
    ammonium_substrate_limit_fraction: f64,
    maximum_ammonium_association_mol_n_per_m3_step: f64,
    ammonium_dissociation_constant_mol_per_m3: f64,
};

pub const Inputs = struct {
    geometry: Geometry,
    aqueous: AqueousState,
    kinetics: KineticControls,
    /// Explicit value of the source scalar `XMIN` retained from a preceding
    /// reaction. SOLUTE.F line 3578 uses it instead of the new band `XMINN`.
    retained_preceding_reaction_reverse_limit_per_step: f64,
};

pub const ZoneStatus = enum {
    dry,
    active,
};

pub const ZoneResult = struct {
    status: ZoneStatus,
    equilibrium_ammonia_activity_mol_n_per_m3: f64,
    fresh_ammonium_dissociation_limit_mol_n_per_m3_step: f64,
    ammonia_association_limit_mol_n_per_m3_step: f64,
    applied_reverse_limit_per_step: f64,
    ammonia_hydrogen_to_ammonium_mol_n_per_m3_step: f64,
};

pub const Result = struct {
    non_band: ZoneResult,
    band: ZoneResult,
};

/// Direct source-order translation of SOLUTE.F lines 3547--3588.
///
/// Positive rates associate NH3 + H into NH4. The band path intentionally
/// applies the retained preceding `XMIN`, while exposing its unused fresh
/// NH4 limit. This comparator is pure and not production-bound.
pub fn calculateSourceOrder(inputs: Inputs) !Result {
    try validate(inputs);
    const non_band =
        if (inputs.geometry.non_band_water_volume_m3 >
        inputs.geometry.minimum_active_water_volume_m3)
            try calculateNonBand(inputs)
        else
            zeroZone();
    const band =
        if (inputs.geometry.band_water_volume_m3 >
        inputs.geometry.minimum_active_water_volume_m3)
            try calculateBand(inputs)
        else
            zeroZone();
    return .{
        .non_band = non_band,
        .band = band,
    };
}

fn calculateNonBand(inputs: Inputs) !ZoneResult {
    const fraction = inputs.kinetics.ammonium_substrate_limit_fraction;
    const aqueous = inputs.aqueous;
    const reverse_limit =
        fraction *
        aqueous.ammonium_non_band_concentration_mol_n_per_m3;
    const forward_limit =
        fraction *
        aqueous.ammonia_non_band_concentration_mol_n_per_m3;
    const equilibrium_ammonia =
        inputs.kinetics.ammonium_dissociation_constant_mol_per_m3 *
        aqueous.ammonium_non_band_activity_mol_n_per_m3 /
        aqueous.hydrogen_activity_mol_per_m3;
    const rate = sourceBoundedRate(
        aqueous.ammonia_non_band_activity_mol_n_per_m3 -
            equilibrium_ammonia,
        reverse_limit,
        forward_limit,
        inputs.kinetics.maximum_ammonium_association_mol_n_per_m3_step,
    );
    return validatedZone(.{
        .status = .active,
        .equilibrium_ammonia_activity_mol_n_per_m3 = equilibrium_ammonia,
        .fresh_ammonium_dissociation_limit_mol_n_per_m3_step = reverse_limit,
        .ammonia_association_limit_mol_n_per_m3_step = forward_limit,
        .applied_reverse_limit_per_step = reverse_limit,
        .ammonia_hydrogen_to_ammonium_mol_n_per_m3_step = rate,
    });
}

fn calculateBand(inputs: Inputs) !ZoneResult {
    const fraction = inputs.kinetics.ammonium_substrate_limit_fraction;
    const aqueous = inputs.aqueous;
    const fresh_reverse_limit =
        fraction * aqueous.ammonium_band_concentration_mol_n_per_m3;
    const forward_limit =
        fraction * aqueous.ammonia_band_concentration_mol_n_per_m3;
    const equilibrium_ammonia =
        inputs.kinetics.ammonium_dissociation_constant_mol_per_m3 *
        aqueous.ammonium_band_activity_mol_n_per_m3 /
        aqueous.hydrogen_activity_mol_per_m3;
    const applied_reverse_limit =
        inputs.retained_preceding_reaction_reverse_limit_per_step;
    const rate = sourceBoundedRate(
        aqueous.ammonia_band_activity_mol_n_per_m3 -
            equilibrium_ammonia,
        applied_reverse_limit,
        forward_limit,
        inputs.kinetics.maximum_ammonium_association_mol_n_per_m3_step,
    );
    return validatedZone(.{
        .status = .active,
        .equilibrium_ammonia_activity_mol_n_per_m3 = equilibrium_ammonia,
        .fresh_ammonium_dissociation_limit_mol_n_per_m3_step = fresh_reverse_limit,
        .ammonia_association_limit_mol_n_per_m3_step = forward_limit,
        .applied_reverse_limit_per_step = applied_reverse_limit,
        .ammonia_hydrogen_to_ammonium_mol_n_per_m3_step = rate,
    });
}

fn sourceBoundedRate(
    driving_activity_mol_n_per_m3: f64,
    reverse_limit_per_step: f64,
    forward_limit_mol_n_per_m3_step: f64,
    maximum_mol_n_per_m3_step: f64,
) f64 {
    return @max(
        -maximum_mol_n_per_m3_step,
        -reverse_limit_per_step,
        @min(
            maximum_mol_n_per_m3_step,
            forward_limit_mol_n_per_m3_step,
            driving_activity_mol_n_per_m3,
        ),
    );
}

fn zeroZone() ZoneResult {
    return .{
        .status = .dry,
        .equilibrium_ammonia_activity_mol_n_per_m3 = 0,
        .fresh_ammonium_dissociation_limit_mol_n_per_m3_step = 0,
        .ammonia_association_limit_mol_n_per_m3_step = 0,
        .applied_reverse_limit_per_step = 0,
        .ammonia_hydrogen_to_ammonium_mol_n_per_m3_step = 0,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Geometry).@"struct".fields) |field|
        try finiteNonnegative(@field(inputs.geometry, field.name));
    inline for (@typeInfo(AqueousState).@"struct".fields) |field|
        try finiteNonnegative(@field(inputs.aqueous, field.name));
    if (inputs.aqueous.hydrogen_activity_mol_per_m3 <= 0)
        return error.InvalidFixedPhAmmoniumAssociationInput;
    inline for (@typeInfo(KineticControls).@"struct".fields) |field|
        try finiteNonnegative(@field(inputs.kinetics, field.name));
    if (inputs.kinetics.ammonium_substrate_limit_fraction > 1)
        return error.InvalidFixedPhAmmoniumAssociationInput;
    try finiteNonnegative(
        inputs.retained_preceding_reaction_reverse_limit_per_step,
    );
}

fn finiteNonnegative(value: f64) !void {
    if (!std.math.isFinite(value) or value < 0)
        return error.InvalidFixedPhAmmoniumAssociationInput;
}

fn validatedZone(zone: ZoneResult) !ZoneResult {
    inline for (@typeInfo(ZoneResult).@"struct".fields) |field| {
        if (field.type == ZoneStatus) continue;
        if (!std.math.isFinite(@field(zone, field.name)))
            return error.NonFiniteFixedPhAmmoniumAssociationResult;
    }
    return zone;
}

fn validInputs() Inputs {
    return .{
        .geometry = .{
            .non_band_water_volume_m3 = 0.7,
            .band_water_volume_m3 = 0.3,
            .minimum_active_water_volume_m3 = 0.1,
        },
        .aqueous = .{
            .ammonium_non_band_concentration_mol_n_per_m3 = 0.3,
            .ammonia_non_band_concentration_mol_n_per_m3 = 0.4,
            .ammonium_band_concentration_mol_n_per_m3 = 0.5,
            .ammonia_band_concentration_mol_n_per_m3 = 0.6,
            .ammonium_non_band_activity_mol_n_per_m3 = 0.24,
            .ammonia_non_band_activity_mol_n_per_m3 = 0.32,
            .ammonium_band_activity_mol_n_per_m3 = 0.4,
            .ammonia_band_activity_mol_n_per_m3 = 0.48,
            .hydrogen_activity_mol_per_m3 = 0.05,
        },
        .kinetics = .{
            .ammonium_substrate_limit_fraction = 0.4,
            .maximum_ammonium_association_mol_n_per_m3_step = 0.2,
            .ammonium_dissociation_constant_mol_per_m3 = 0.1,
        },
        .retained_preceding_reaction_reverse_limit_per_step = 0.03,
    };
}

test "fixed-pH ammonium association matches every source equation" {
    const inputs = validInputs();
    const result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(ZoneStatus.active, result.non_band.status);
    try std.testing.expectEqual(ZoneStatus.active, result.band.status);

    const non_band_equilibrium =
        inputs.kinetics.ammonium_dissociation_constant_mol_per_m3 *
        inputs.aqueous.ammonium_non_band_activity_mol_n_per_m3 /
        inputs.aqueous.hydrogen_activity_mol_per_m3;
    const band_equilibrium =
        inputs.kinetics.ammonium_dissociation_constant_mol_per_m3 *
        inputs.aqueous.ammonium_band_activity_mol_n_per_m3 /
        inputs.aqueous.hydrogen_activity_mol_per_m3;
    try std.testing.expectEqual(
        non_band_equilibrium,
        result.non_band.equilibrium_ammonia_activity_mol_n_per_m3,
    );
    try std.testing.expectEqual(
        band_equilibrium,
        result.band.equilibrium_ammonia_activity_mol_n_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 0.12),
        result.non_band
            .fresh_ammonium_dissociation_limit_mol_n_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.2),
        result.band.fresh_ammonium_dissociation_limit_mol_n_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, -0.12),
        result.non_band
            .ammonia_hydrogen_to_ammonium_mol_n_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, -0.03),
        result.band.ammonia_hydrogen_to_ammonium_mol_n_per_m3_step,
    );
}

test "fixed-pH ammonium association gates each runtime zone independently" {
    var inputs = validInputs();
    inputs.geometry.non_band_water_volume_m3 =
        inputs.geometry.minimum_active_water_volume_m3;
    var result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(ZoneStatus.dry, result.non_band.status);
    try std.testing.expectEqual(ZoneStatus.active, result.band.status);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.non_band
            .ammonia_hydrogen_to_ammonium_mol_n_per_m3_step,
    );

    inputs = validInputs();
    inputs.geometry.band_water_volume_m3 =
        inputs.geometry.minimum_active_water_volume_m3;
    result = try calculateSourceOrder(inputs);
    try std.testing.expectEqual(ZoneStatus.active, result.non_band.status);
    try std.testing.expectEqual(ZoneStatus.dry, result.band.status);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.band.ammonia_hydrogen_to_ammonium_mol_n_per_m3_step,
    );
}

test "fixed-pH band association exposes retained reverse-bound dependency" {
    const inputs = validInputs();
    const source = try calculateSourceOrder(inputs);
    var fresh_inputs = inputs;
    fresh_inputs.retained_preceding_reaction_reverse_limit_per_step =
        inputs.kinetics.ammonium_substrate_limit_fraction *
        inputs.aqueous.ammonium_band_concentration_mol_n_per_m3;
    const fresh = try calculateSourceOrder(fresh_inputs);
    try std.testing.expectEqual(
        @as(f64, -0.03),
        source.band.ammonia_hydrogen_to_ammonium_mol_n_per_m3_step,
    );
    try std.testing.expectEqual(
        @as(f64, -0.2),
        fresh.band.ammonia_hydrogen_to_ammonium_mol_n_per_m3_step,
    );
}

test "fixed-pH ammonium association rejects invalid runtime state" {
    var inputs = validInputs();
    inputs.aqueous.hydrogen_activity_mol_per_m3 = 0;
    try std.testing.expectError(
        error.InvalidFixedPhAmmoniumAssociationInput,
        calculateSourceOrder(inputs),
    );
    inputs = validInputs();
    inputs.retained_preceding_reaction_reverse_limit_per_step =
        std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidFixedPhAmmoniumAssociationInput,
        calculateSourceOrder(inputs),
    );
}
