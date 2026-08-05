const std = @import("std");

pub const Inputs = struct {
    litter_water_volume_m3: f64,
    volumetric_litter_water_m3_per_m3: f64,
    litter_dry_mass_megagrams: f64,
    topsoil_dry_mass_megagrams: f64,
    ordinary_organic_carbon_g_c: f64,
    charcoal_carbon_g_c: f64,
    carboxyl_sites_mol_per_megagram_c: f64,
    active_water_volume_threshold_m3: f64,
    active_litter_dry_mass_threshold_megagrams: f64,
    cation_exchange_capacity_floor_mol_charge_per_megagram_litter: f64,
};

pub const WetState = struct {
    volumetric_litter_water_m3_per_m3: f64,
    litter_dry_mass_megagrams: f64,
    litter_mass_per_water_volume_megagrams_per_m3: f64,
    topsoil_dry_mass_megagrams: f64,
    cation_exchange_capacity_mol_charge_per_megagram_litter: f64,
};

pub const Result = union(enum) {
    inactive,
    wet: WetState,
};

/// Direct source-order translation of SOLUTE.F lines 4018--4036.
///
/// The caller owns the runtime horizontal-cell traversal. This scalar kernel
/// first applies the source's strict litter-water gate, then reconstructs the
/// wet-litter geometry and carboxyl-derived cation-exchange capacity without
/// mutating model state.
pub fn calculate(inputs: Inputs) !Result {
    try validateInputs(inputs);
    if (inputs.litter_water_volume_m3 <=
        inputs.active_water_volume_threshold_m3)
        return .inactive;

    // SOLUTE.F 4027--4030. Preserve assignment and division order.
    const volumetric_water = inputs.volumetric_litter_water_m3_per_m3;
    const litter_mass = inputs.litter_dry_mass_megagrams;
    const litter_mass_per_water =
        litter_mass / inputs.litter_water_volume_m3;
    const topsoil_mass = inputs.topsoil_dry_mass_megagrams;

    // SOLUTE.F 4031--4036. The dry-mass gate is strict and the source floor
    // applies only inside its admitted branch.
    const cation_exchange_capacity = if (litter_mass >
        inputs.active_litter_dry_mass_threshold_megagrams)
        @max(
            inputs.cation_exchange_capacity_floor_mol_charge_per_megagram_litter,
            inputs.carboxyl_sites_mol_per_megagram_c *
                1.0e-6 *
                (inputs.ordinary_organic_carbon_g_c +
                    inputs.charcoal_carbon_g_c) /
                litter_mass,
        )
    else
        0.0;

    const result: WetState = .{
        .volumetric_litter_water_m3_per_m3 = volumetric_water,
        .litter_dry_mass_megagrams = litter_mass,
        .litter_mass_per_water_volume_megagrams_per_m3 = litter_mass_per_water,
        .topsoil_dry_mass_megagrams = topsoil_mass,
        .cation_exchange_capacity_mol_charge_per_megagram_litter = cation_exchange_capacity,
    };
    inline for (@typeInfo(WetState).@"struct".fields) |field| {
        const value = @field(result, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterReactionInitialization;
    }
    return .{ .wet = result };
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterReactionInitializationInput;
    }
    if (inputs.cation_exchange_capacity_floor_mol_charge_per_megagram_litter == 0)
        return error.InvalidSurfaceLitterReactionInitializationInput;
}

fn testInputs() Inputs {
    return .{
        .litter_water_volume_m3 = 0.25,
        .volumetric_litter_water_m3_per_m3 = 0.4,
        .litter_dry_mass_megagrams = 0.5,
        .topsoil_dry_mass_megagrams = 10,
        .ordinary_organic_carbon_g_c = 200_000,
        .charcoal_carbon_g_c = 50_000,
        .carboxyl_sites_mol_per_megagram_c = 250,
        // Source ZEROS2 supplies the same numeric threshold to water and
        // mass comparisons; separate typed fields make those units explicit.
        .active_water_volume_threshold_m3 = 1.0e-12,
        .active_litter_dry_mass_threshold_megagrams = 1.0e-12,
        .cation_exchange_capacity_floor_mol_charge_per_megagram_litter = 1.0e-32,
    };
}

test "SOLUTE wet surface-litter initialization preserves every source expression" {
    const inputs = testInputs();
    const result = try calculate(inputs);
    const wet = switch (result) {
        .inactive => return error.ExpectedWetSurfaceLitter,
        .wet => |value| value,
    };

    const expected_mass_per_water =
        inputs.litter_dry_mass_megagrams / inputs.litter_water_volume_m3;
    const expected_capacity = @max(
        inputs.cation_exchange_capacity_floor_mol_charge_per_megagram_litter,
        inputs.carboxyl_sites_mol_per_megagram_c *
            1.0e-6 *
            (inputs.ordinary_organic_carbon_g_c +
                inputs.charcoal_carbon_g_c) /
            inputs.litter_dry_mass_megagrams,
    );
    try std.testing.expectEqual(
        inputs.volumetric_litter_water_m3_per_m3,
        wet.volumetric_litter_water_m3_per_m3,
    );
    try std.testing.expectEqual(
        inputs.litter_dry_mass_megagrams,
        wet.litter_dry_mass_megagrams,
    );
    try std.testing.expectEqual(
        expected_mass_per_water,
        wet.litter_mass_per_water_volume_megagrams_per_m3,
    );
    try std.testing.expectEqual(
        inputs.topsoil_dry_mass_megagrams,
        wet.topsoil_dry_mass_megagrams,
    );
    try std.testing.expectEqual(
        expected_capacity,
        wet.cation_exchange_capacity_mol_charge_per_megagram_litter,
    );
}

test "surface-litter CEC closes carboxyl charge on the litter-mass basis" {
    const inputs = testInputs();
    const result = try calculate(inputs);
    const wet = result.wet;
    const represented_carboxyl_charge_mol =
        inputs.carboxyl_sites_mol_per_megagram_c *
        1.0e-6 *
        (inputs.ordinary_organic_carbon_g_c +
            inputs.charcoal_carbon_g_c);
    const reconstructed_charge_mol =
        wet.cation_exchange_capacity_mol_charge_per_megagram_litter *
        wet.litter_dry_mass_megagrams;
    try std.testing.expectApproxEqAbs(
        represented_carboxyl_charge_mol,
        reconstructed_charge_mol,
        1.0e-14,
    );
}

test "surface-litter water and dry-mass gates retain strict source boundaries" {
    var inputs = testInputs();
    inputs.litter_water_volume_m3 =
        inputs.active_water_volume_threshold_m3;
    switch (try calculate(inputs)) {
        .inactive => {},
        .wet => return error.ExpectedInactiveSurfaceLitter,
    }

    inputs = testInputs();
    inputs.litter_dry_mass_megagrams =
        inputs.active_litter_dry_mass_threshold_megagrams;
    const threshold_wet = (try calculate(inputs)).wet;
    try std.testing.expectEqual(
        @as(f64, 0),
        threshold_wet.cation_exchange_capacity_mol_charge_per_megagram_litter,
    );

    inputs = testInputs();
    inputs.ordinary_organic_carbon_g_c = 0;
    inputs.charcoal_carbon_g_c = 0;
    const zero_carbon_wet = (try calculate(inputs)).wet;
    try std.testing.expectEqual(
        inputs.cation_exchange_capacity_floor_mol_charge_per_megagram_litter,
        zero_carbon_wet.cation_exchange_capacity_mol_charge_per_megagram_litter,
    );
}

test "surface-litter initialization rejects invalid inputs and overflow" {
    var inputs = testInputs();
    inputs.litter_water_volume_m3 = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterReactionInitializationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.ordinary_organic_carbon_g_c = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterReactionInitializationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.cation_exchange_capacity_floor_mol_charge_per_megagram_litter = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterReactionInitializationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.ordinary_organic_carbon_g_c = std.math.floatMax(f64);
    inputs.charcoal_carbon_g_c = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterReactionInitialization,
        calculate(inputs),
    );
}
