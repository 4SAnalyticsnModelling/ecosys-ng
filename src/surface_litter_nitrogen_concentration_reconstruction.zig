const std = @import("std");

pub const Inputs = struct {
    litter_water_volume_m3: f64,
    nitrogen_molar_mass_g_n_per_mol: f64,
    ammonium_net_change_g_n_per_step: f64,
    ammonium_fertilizer_dissolution_mol_n_per_step: f64,
    ammonia_from_urea_hydrolysis_mol_n_per_step: f64,
    aqueous_ammonium_inventory_g_n: f64,
    aqueous_ammonia_inventory_g_n: f64,
    exchangeable_ammonium_inventory_mol_n: f64,
    litter_dry_mass_megagrams: f64,
    active_litter_dry_mass_threshold_megagrams: f64,
    aqueous_nitrogen_mass_floor_g_n: f64,
    exchangeable_ammonium_floor_mol_n_per_megagram: f64,
};

pub const Result = struct {
    nitrogen_water_scale_g_n_m3_per_mol: f64,
    ammonium_input_g_n_per_step: f64,
    ammonia_input_g_n_per_step: f64,
    aqueous_ammonium_mol_n_per_m3: f64,
    aqueous_ammonia_mol_n_per_m3: f64,
    exchangeable_ammonium_mol_n_per_megagram: f64,
    aqueous_ammonium_floor_was_applied: bool,
    aqueous_ammonia_floor_was_applied: bool,
    exchangeable_ammonium_floor_was_applied: bool,
};

/// Direct source-order translation of SOLUTE.F lines 4127--4147.
///
/// The caller selects one runtime horizontal cell already admitted by the
/// surrounding wet-litter branch. No inventory is mutated by this kernel.
pub fn calculate(inputs: Inputs) !Result {
    try validateInputs(inputs);

    // SOLUTE.F 4138--4142.
    const nitrogen_water_scale =
        inputs.nitrogen_molar_mass_g_n_per_mol *
        inputs.litter_water_volume_m3;
    const ammonium_input =
        inputs.ammonium_net_change_g_n_per_step +
        inputs.nitrogen_molar_mass_g_n_per_mol *
            inputs.ammonium_fertilizer_dissolution_mol_n_per_step;
    const ammonia_input =
        inputs.nitrogen_molar_mass_g_n_per_mol *
        inputs.ammonia_from_urea_hydrolysis_mol_n_per_step;
    const aqueous_ammonium_mass =
        inputs.aqueous_ammonium_inventory_g_n + ammonium_input;
    const aqueous_ammonia_mass =
        inputs.aqueous_ammonia_inventory_g_n + ammonia_input;
    try validateIntermediate(.{
        nitrogen_water_scale,
        ammonium_input,
        ammonia_input,
        aqueous_ammonium_mass,
        aqueous_ammonia_mass,
    });

    const ammonium_concentration =
        @max(inputs.aqueous_nitrogen_mass_floor_g_n, aqueous_ammonium_mass) /
        nitrogen_water_scale;
    const ammonia_concentration =
        @max(inputs.aqueous_nitrogen_mass_floor_g_n, aqueous_ammonia_mass) /
        nitrogen_water_scale;

    // SOLUTE.F 4143--4147. Threshold equality selects the dry branch.
    var exchangeable_ammonium_concentration: f64 = 0;
    var exchangeable_floor_was_applied = false;
    if (inputs.litter_dry_mass_megagrams >
        inputs.active_litter_dry_mass_threshold_megagrams)
    {
        const unconstrained_exchangeable =
            inputs.exchangeable_ammonium_inventory_mol_n /
            inputs.litter_dry_mass_megagrams;
        exchangeable_ammonium_concentration = @max(
            inputs.exchangeable_ammonium_floor_mol_n_per_megagram,
            unconstrained_exchangeable,
        );
        exchangeable_floor_was_applied =
            unconstrained_exchangeable <
            inputs.exchangeable_ammonium_floor_mol_n_per_megagram;
    }

    const result: Result = .{
        .nitrogen_water_scale_g_n_m3_per_mol = nitrogen_water_scale,
        .ammonium_input_g_n_per_step = ammonium_input,
        .ammonia_input_g_n_per_step = ammonia_input,
        .aqueous_ammonium_mol_n_per_m3 = ammonium_concentration,
        .aqueous_ammonia_mol_n_per_m3 = ammonia_concentration,
        .exchangeable_ammonium_mol_n_per_megagram = exchangeable_ammonium_concentration,
        .aqueous_ammonium_floor_was_applied = aqueous_ammonium_mass <
            inputs.aqueous_nitrogen_mass_floor_g_n,
        .aqueous_ammonia_floor_was_applied = aqueous_ammonia_mass <
            inputs.aqueous_nitrogen_mass_floor_g_n,
        .exchangeable_ammonium_floor_was_applied = exchangeable_floor_was_applied,
    };
    try validateResult(result);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidSurfaceLitterNitrogenConcentrationInput;
    }
    inline for (.{
        inputs.ammonium_fertilizer_dissolution_mol_n_per_step,
        inputs.ammonia_from_urea_hydrolysis_mol_n_per_step,
        inputs.aqueous_ammonium_inventory_g_n,
        inputs.aqueous_ammonia_inventory_g_n,
        inputs.exchangeable_ammonium_inventory_mol_n,
        inputs.litter_dry_mass_megagrams,
        inputs.active_litter_dry_mass_threshold_megagrams,
        inputs.aqueous_nitrogen_mass_floor_g_n,
        inputs.exchangeable_ammonium_floor_mol_n_per_megagram,
    }) |value| {
        if (value < 0)
            return error.InvalidSurfaceLitterNitrogenConcentrationInput;
    }
    if (inputs.litter_water_volume_m3 <= 0 or
        inputs.nitrogen_molar_mass_g_n_per_mol <= 0)
    {
        return error.InvalidSurfaceLitterNitrogenConcentrationInput;
    }
}

fn validateIntermediate(values: anytype) !void {
    inline for (values) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterNitrogenConcentration;
    }
}

fn validateResult(result: Result) !void {
    inline for (.{
        result.nitrogen_water_scale_g_n_m3_per_mol,
        result.ammonium_input_g_n_per_step,
        result.ammonia_input_g_n_per_step,
        result.aqueous_ammonium_mol_n_per_m3,
        result.aqueous_ammonia_mol_n_per_m3,
        result.exchangeable_ammonium_mol_n_per_megagram,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterNitrogenConcentration;
    }
    if (result.nitrogen_water_scale_g_n_m3_per_mol <= 0 or
        result.aqueous_ammonium_mol_n_per_m3 < 0 or
        result.aqueous_ammonia_mol_n_per_m3 < 0 or
        result.exchangeable_ammonium_mol_n_per_megagram < 0)
    {
        return error.InvalidSurfaceLitterNitrogenConcentrationResult;
    }
}

fn testInputs() Inputs {
    return .{
        .litter_water_volume_m3 = 2,
        .nitrogen_molar_mass_g_n_per_mol = 14,
        .ammonium_net_change_g_n_per_step = -1,
        .ammonium_fertilizer_dissolution_mol_n_per_step = 0.5,
        .ammonia_from_urea_hydrolysis_mol_n_per_step = 0.25,
        .aqueous_ammonium_inventory_g_n = 22,
        .aqueous_ammonia_inventory_g_n = 10.5,
        .exchangeable_ammonium_inventory_mol_n = 3,
        .litter_dry_mass_megagrams = 1.5,
        .active_litter_dry_mass_threshold_megagrams = 1.0e-12,
        .aqueous_nitrogen_mass_floor_g_n = 1.0e-20,
        .exchangeable_ammonium_floor_mol_n_per_megagram = 1.0e-20,
    };
}

test "SOLUTE surface nitrogen reconstruction preserves every source expression" {
    const inputs = testInputs();
    const result = try calculate(inputs);

    const expected_scale =
        inputs.nitrogen_molar_mass_g_n_per_mol *
        inputs.litter_water_volume_m3;
    const expected_ammonium_input =
        inputs.ammonium_net_change_g_n_per_step +
        inputs.nitrogen_molar_mass_g_n_per_mol *
            inputs.ammonium_fertilizer_dissolution_mol_n_per_step;
    const expected_ammonia_input =
        inputs.nitrogen_molar_mass_g_n_per_mol *
        inputs.ammonia_from_urea_hydrolysis_mol_n_per_step;
    const expected_ammonium =
        @max(
            inputs.aqueous_nitrogen_mass_floor_g_n,
            inputs.aqueous_ammonium_inventory_g_n +
                expected_ammonium_input,
        ) / expected_scale;
    const expected_ammonia =
        @max(
            inputs.aqueous_nitrogen_mass_floor_g_n,
            inputs.aqueous_ammonia_inventory_g_n +
                expected_ammonia_input,
        ) / expected_scale;
    const expected_exchangeable = @max(
        inputs.exchangeable_ammonium_floor_mol_n_per_megagram,
        inputs.exchangeable_ammonium_inventory_mol_n /
            inputs.litter_dry_mass_megagrams,
    );

    try std.testing.expectEqual(
        expected_scale,
        result.nitrogen_water_scale_g_n_m3_per_mol,
    );
    try std.testing.expectEqual(
        expected_ammonium_input,
        result.ammonium_input_g_n_per_step,
    );
    try std.testing.expectEqual(
        expected_ammonia_input,
        result.ammonia_input_g_n_per_step,
    );
    try std.testing.expectEqual(
        expected_ammonium,
        result.aqueous_ammonium_mol_n_per_m3,
    );
    try std.testing.expectEqual(
        expected_ammonia,
        result.aqueous_ammonia_mol_n_per_m3,
    );
    try std.testing.expectEqual(
        expected_exchangeable,
        result.exchangeable_ammonium_mol_n_per_megagram,
    );
}

test "surface nitrogen reconstruction recovers valid extensive masses" {
    const inputs = testInputs();
    const result = try calculate(inputs);

    try std.testing.expectApproxEqAbs(
        inputs.aqueous_ammonium_inventory_g_n +
            result.ammonium_input_g_n_per_step,
        result.aqueous_ammonium_mol_n_per_m3 *
            result.nitrogen_water_scale_g_n_m3_per_mol,
        1.0e-14,
    );
    try std.testing.expectApproxEqAbs(
        inputs.aqueous_ammonia_inventory_g_n +
            result.ammonia_input_g_n_per_step,
        result.aqueous_ammonia_mol_n_per_m3 *
            result.nitrogen_water_scale_g_n_m3_per_mol,
        1.0e-14,
    );
    try std.testing.expectApproxEqAbs(
        inputs.exchangeable_ammonium_inventory_mol_n,
        result.exchangeable_ammonium_mol_n_per_megagram *
            inputs.litter_dry_mass_megagrams,
        1.0e-14,
    );
}

test "surface nitrogen floors and dry-mass threshold are explicit" {
    var inputs = testInputs();
    inputs.ammonium_net_change_g_n_per_step = -100;
    inputs.ammonium_fertilizer_dissolution_mol_n_per_step = 0;
    inputs.aqueous_ammonia_inventory_g_n = 0;
    inputs.ammonia_from_urea_hydrolysis_mol_n_per_step = 0;
    inputs.exchangeable_ammonium_inventory_mol_n = 0;
    const floored = try calculate(inputs);
    try std.testing.expect(floored.aqueous_ammonium_floor_was_applied);
    try std.testing.expect(floored.aqueous_ammonia_floor_was_applied);
    try std.testing.expect(floored.exchangeable_ammonium_floor_was_applied);

    inputs = testInputs();
    inputs.litter_dry_mass_megagrams =
        inputs.active_litter_dry_mass_threshold_megagrams;
    const dry = try calculate(inputs);
    try std.testing.expectEqual(
        @as(f64, 0),
        dry.exchangeable_ammonium_mol_n_per_megagram,
    );
    try std.testing.expect(!dry.exchangeable_ammonium_floor_was_applied);
}

test "surface nitrogen reconstruction rejects invalid input and overflow" {
    var inputs = testInputs();
    inputs.litter_water_volume_m3 = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterNitrogenConcentrationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.ammonium_net_change_g_n_per_step = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterNitrogenConcentrationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.aqueous_ammonium_inventory_g_n = -1;
    try std.testing.expectError(
        error.InvalidSurfaceLitterNitrogenConcentrationInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.nitrogen_molar_mass_g_n_per_mol = std.math.floatMax(f64);
    inputs.litter_water_volume_m3 = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterNitrogenConcentration,
        calculate(inputs),
    );
}
