const std = @import("std");

pub const Inputs = struct {
    biologically_active_water_volume_m3: f64,
    active_biomass_respiration_g_c_per_step: f64,
    step_duration_h: f64,
    minimum_urea_half_saturation_mol_n_per_Mg: f64,
    microbial_activity_inhibition_g_c_per_m3_h: f64,
    initial_urease_inhibition_fraction: f64,
    current_urease_inhibition_fraction: f64,
    inhibition_decline_fraction_per_step: f64,
    urea_fertilizer_mol_n: f64,
    litter_dry_mass_Mg: f64,
    microbial_temperature_response_fraction: f64,
    specific_urea_hydrolysis_mol_n_per_g_c: f64,
    active_water_volume_threshold_m3: f64,
    active_inhibition_threshold_fraction: f64,
    active_urea_threshold_mol_n: f64,
    active_litter_dry_mass_threshold_Mg: f64,
};

pub const Result = struct {
    microbial_activity_g_c_per_m3_h: f64,
    effective_urea_half_saturation_mol_n_per_Mg: f64,
    updated_urease_inhibition_fraction: f64,
    urea_concentration_mol_n_per_Mg: ?f64,
    urea_substrate_response_fraction: ?f64,
    urea_hydrolysis_mol_n_per_step: f64,
};

/// Direct source-order translation of SOLUTE.F lines 4038--4096.
///
/// The caller selects one runtime horizontal cell. This pure kernel stages
/// the source mutation of current urease inhibition and its urea-hydrolysis
/// flux, leaving publication to a later authoritative owner.
pub fn calculate(inputs: Inputs) !Result {
    try validateInputs(inputs);

    // SOLUTE.F 4050--4055.
    const microbial_activity = if (inputs.biologically_active_water_volume_m3 >
        inputs.active_water_volume_threshold_m3)
    wet: {
        const active_water_time_m3_h =
            inputs.biologically_active_water_volume_m3 *
            inputs.step_duration_h;
        if (!std.math.isFinite(active_water_time_m3_h))
            return error.NonFiniteSurfaceLitterUreaHydrolysis;
        break :wet @min(
            0.1e+06,
            inputs.active_biomass_respiration_g_c_per_step /
                active_water_time_m3_h,
        );
    } else 0.1e+06;
    const effective_half_saturation =
        inputs.minimum_urea_half_saturation_mol_n_per_Mg *
        (1.0 + microbial_activity /
            inputs.microbial_activity_inhibition_g_c_per_m3_h);

    // SOLUTE.F 4064--4071. Every right-hand occurrence reads the preceding
    // inhibition value before the one staged assignment.
    const updated_inhibition =
        if (inputs.initial_urease_inhibition_fraction >
        inputs.active_inhibition_threshold_fraction and
        inputs.current_urease_inhibition_fraction >
            inputs.active_inhibition_threshold_fraction)
            inputs.current_urease_inhibition_fraction -
                inputs.inhibition_decline_fraction_per_step *
                    inputs.current_urease_inhibition_fraction *
                    @max(
                        inputs.inhibition_decline_fraction_per_step,
                        1.0 - inputs.current_urease_inhibition_fraction /
                            inputs.initial_urease_inhibition_fraction,
                    )
        else
            0.0;

    var urea_concentration: ?f64 = null;
    var substrate_response: ?f64 = null;
    var hydrolysis: f64 = 0;
    // SOLUTE.F 4088--4096. Threshold equality selects the inactive branch.
    if (inputs.urea_fertilizer_mol_n > inputs.active_urea_threshold_mol_n and
        inputs.litter_dry_mass_Mg >
            inputs.active_litter_dry_mass_threshold_Mg)
    {
        const concentration =
            inputs.urea_fertilizer_mol_n / inputs.litter_dry_mass_Mg;
        const response =
            concentration / (concentration + effective_half_saturation);
        const potential_hydrolysis =
            inputs.specific_urea_hydrolysis_mol_n_per_g_c *
            inputs.active_biomass_respiration_g_c_per_step *
            response *
            inputs.microbial_temperature_response_fraction *
            (1.0 - updated_inhibition);
        if (!std.math.isFinite(potential_hydrolysis))
            return error.NonFiniteSurfaceLitterUreaHydrolysis;
        hydrolysis = @min(
            inputs.urea_fertilizer_mol_n,
            potential_hydrolysis,
        );
        urea_concentration = concentration;
        substrate_response = response;
    }

    const result: Result = .{
        .microbial_activity_g_c_per_m3_h = microbial_activity,
        .effective_urea_half_saturation_mol_n_per_Mg = effective_half_saturation,
        .updated_urease_inhibition_fraction = updated_inhibition,
        .urea_concentration_mol_n_per_Mg = urea_concentration,
        .urea_substrate_response_fraction = substrate_response,
        .urea_hydrolysis_mol_n_per_step = hydrolysis,
    };
    try validateResult(result, inputs.urea_fertilizer_mol_n);
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSurfaceLitterUreaHydrolysisInput;
    }
    if (inputs.step_duration_h == 0 or
        inputs.minimum_urea_half_saturation_mol_n_per_Mg == 0 or
        inputs.microbial_activity_inhibition_g_c_per_m3_h == 0 or
        inputs.initial_urease_inhibition_fraction > 1 or
        inputs.current_urease_inhibition_fraction > 1 or
        inputs.inhibition_decline_fraction_per_step > 1 or
        inputs.microbial_temperature_response_fraction > 1 or
        inputs.active_inhibition_threshold_fraction > 1)
    {
        return error.InvalidSurfaceLitterUreaHydrolysisInput;
    }
}

fn validateResult(result: Result, available_urea_mol_n: f64) !void {
    inline for (.{
        result.microbial_activity_g_c_per_m3_h,
        result.effective_urea_half_saturation_mol_n_per_Mg,
        result.updated_urease_inhibition_fraction,
        result.urea_hydrolysis_mol_n_per_step,
    }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterUreaHydrolysis;
    }
    inline for (.{
        result.urea_concentration_mol_n_per_Mg,
        result.urea_substrate_response_fraction,
    }) |optional| if (optional) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceLitterUreaHydrolysis;
    };
    if (result.updated_urease_inhibition_fraction < 0 or
        result.updated_urease_inhibition_fraction > 1 or
        result.urea_hydrolysis_mol_n_per_step < 0 or
        result.urea_hydrolysis_mol_n_per_step > available_urea_mol_n)
    {
        return error.InvalidSurfaceLitterUreaHydrolysisResult;
    }
}

fn testInputs() Inputs {
    return .{
        .biologically_active_water_volume_m3 = 2,
        .active_biomass_respiration_g_c_per_step = 30,
        .step_duration_h = 0.5,
        .minimum_urea_half_saturation_mol_n_per_Mg = 4,
        .microbial_activity_inhibition_g_c_per_m3_h = 10,
        .initial_urease_inhibition_fraction = 0.8,
        .current_urease_inhibition_fraction = 0.6,
        .inhibition_decline_fraction_per_step = 0.1,
        .urea_fertilizer_mol_n = 12,
        .litter_dry_mass_Mg = 3,
        .microbial_temperature_response_fraction = 0.7,
        .specific_urea_hydrolysis_mol_n_per_g_c = 0.2,
        .active_water_volume_threshold_m3 = 1.0e-12,
        .active_inhibition_threshold_fraction = 1.0e-12,
        .active_urea_threshold_mol_n = 1.0e-12,
        .active_litter_dry_mass_threshold_Mg = 1.0e-12,
    };
}

test "SOLUTE surface urea hydrolysis preserves every source expression" {
    const inputs = testInputs();
    const result = try calculate(inputs);

    const expected_activity = @min(
        0.1e+06,
        inputs.active_biomass_respiration_g_c_per_step /
            (inputs.biologically_active_water_volume_m3 *
                inputs.step_duration_h),
    );
    const expected_half_saturation =
        inputs.minimum_urea_half_saturation_mol_n_per_Mg *
        (1.0 + expected_activity /
            inputs.microbial_activity_inhibition_g_c_per_m3_h);
    const expected_inhibition =
        inputs.current_urease_inhibition_fraction -
        inputs.inhibition_decline_fraction_per_step *
            inputs.current_urease_inhibition_fraction *
            @max(
                inputs.inhibition_decline_fraction_per_step,
                1.0 - inputs.current_urease_inhibition_fraction /
                    inputs.initial_urease_inhibition_fraction,
            );
    const expected_concentration =
        inputs.urea_fertilizer_mol_n / inputs.litter_dry_mass_Mg;
    const expected_response = expected_concentration /
        (expected_concentration + expected_half_saturation);
    const expected_hydrolysis = @min(
        inputs.urea_fertilizer_mol_n,
        inputs.specific_urea_hydrolysis_mol_n_per_g_c *
            inputs.active_biomass_respiration_g_c_per_step *
            expected_response *
            inputs.microbial_temperature_response_fraction *
            (1.0 - expected_inhibition),
    );

    try std.testing.expectEqual(
        expected_activity,
        result.microbial_activity_g_c_per_m3_h,
    );
    try std.testing.expectEqual(
        expected_half_saturation,
        result.effective_urea_half_saturation_mol_n_per_Mg,
    );
    try std.testing.expectEqual(
        expected_inhibition,
        result.updated_urease_inhibition_fraction,
    );
    try std.testing.expectEqual(
        expected_concentration,
        result.urea_concentration_mol_n_per_Mg.?,
    );
    try std.testing.expectEqual(
        expected_response,
        result.urea_substrate_response_fraction.?,
    );
    try std.testing.expectEqual(
        expected_hydrolysis,
        result.urea_hydrolysis_mol_n_per_step,
    );
}

test "zero active water uses the source microbial-activity cap" {
    var inputs = testInputs();
    inputs.biologically_active_water_volume_m3 = 0;
    const result = try calculate(inputs);
    try std.testing.expectEqual(
        @as(f64, 100_000),
        result.microbial_activity_g_c_per_m3_h,
    );
    try std.testing.expect(result.urea_hydrolysis_mol_n_per_step > 0);
}

test "surface urea thresholds are strict and donor transfer conserves nitrogen" {
    var inputs = testInputs();
    const active = try calculate(inputs);
    const remaining_urea_mol_n =
        inputs.urea_fertilizer_mol_n -
        active.urea_hydrolysis_mol_n_per_step;
    try std.testing.expectApproxEqAbs(
        inputs.urea_fertilizer_mol_n,
        remaining_urea_mol_n +
            active.urea_hydrolysis_mol_n_per_step,
        1.0e-14,
    );

    inputs.current_urease_inhibition_fraction =
        inputs.active_inhibition_threshold_fraction;
    const inactive_inhibition = try calculate(inputs);
    try std.testing.expectEqual(
        @as(f64, 0),
        inactive_inhibition.updated_urease_inhibition_fraction,
    );

    inputs = testInputs();
    inputs.urea_fertilizer_mol_n = inputs.active_urea_threshold_mol_n;
    const inactive_urea = try calculate(inputs);
    try std.testing.expect(inactive_urea.urea_concentration_mol_n_per_Mg == null);
    try std.testing.expectEqual(
        @as(f64, 0),
        inactive_urea.urea_hydrolysis_mol_n_per_step,
    );

    inputs = testInputs();
    inputs.litter_dry_mass_Mg =
        inputs.active_litter_dry_mass_threshold_Mg;
    const inactive_mass = try calculate(inputs);
    try std.testing.expect(inactive_mass.urea_substrate_response_fraction == null);
    try std.testing.expectEqual(
        @as(f64, 0),
        inactive_mass.urea_hydrolysis_mol_n_per_step,
    );
}

test "surface urea hydrolysis rejects invalid inputs and overflow" {
    var inputs = testInputs();
    inputs.step_duration_h = 0;
    try std.testing.expectError(
        error.InvalidSurfaceLitterUreaHydrolysisInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.current_urease_inhibition_fraction = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSurfaceLitterUreaHydrolysisInput,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.minimum_urea_half_saturation_mol_n_per_Mg =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterUreaHydrolysis,
        calculate(inputs),
    );

    inputs = testInputs();
    inputs.active_biomass_respiration_g_c_per_step =
        std.math.floatMax(f64);
    inputs.specific_urea_hydrolysis_mol_n_per_g_c =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterUreaHydrolysis,
        calculate(inputs),
    );
}
