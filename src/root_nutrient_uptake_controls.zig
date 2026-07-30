const std = @import("std");

pub const Admission = struct {
    layer_thickness_m: f64,
    minimum_layer_thickness_m: f64,
    root_length_density_m_per_m3: f64,
    root_aqueous_volume_m3: f64,
    soil_water_content_m3_per_m3: f64,
    negligible_root_amount: f64,
    negligible_soil_water_content_m3_per_m3: f64,
};

pub const Protein = struct {
    protein_content_g: f64,
    root_carbon_mass_g_c: f64,
    maximum_protein_concentration_g_per_g_c: f64,
    reference_protein_concentration_g_per_g_c: f64,
};

pub const Respiration = struct {
    unconstrained_uptake_respiration_g_c_per_step: f64,
    nonstructural_carbon_pool_g_c: f64,
    biological_timestep_h_per_step: f64,
};

pub const NutrientFeedback = struct {
    nonstructural_carbon_g_c_per_g_c: f64,
    nonstructural_nitrogen_g_n_per_g_c: f64,
    nonstructural_phosphorus_g_p_per_g_c: f64,
    nutrient_uptake_strategy_code: u8,
    vascular_carbon_g_c: f64,
    total_root_carbon_g_c: f64,
    vascular_fraction_threshold: f64,
    nitrogen_inhibition_by_nitrogen_g_n_per_g_c: f64,
    nitrogen_inhibition_by_phosphorus_g_n_per_g_c: f64,
    phosphorus_inhibition_by_phosphorus_g_p_per_g_c: f64,
    phosphorus_inhibition_by_nitrogen_g_p_per_g_c: f64,
    significant_concentration_g_per_g_c: f64,
};

pub const Water = struct {
    root_water_uptake_m3_per_step: f64,
    plant_population: f64,
    gas_flux_timestep_h_per_step: f64,
};

pub const UptakeAccumulators = struct {
    oxygen: f64 = 0,
    ammonium_non_band: f64 = 0,
    ammonium_band: f64 = 0,
    nitrate_non_band: f64 = 0,
    nitrate_band: f64 = 0,
    phosphate_non_band: f64 = 0,
    phosphate_band: f64 = 0,
    phosphate_1_non_band: f64 = 0,
    phosphate_1_band: f64 = 0,
};

pub const Active = struct {
    accumulators: UptakeAccumulators,
    protein_concentration_g_per_g_c: f64,
    relative_protein_uptake_capacity: f64,
    carbon_respiration_limitation: f64,
    nitrogen_feedback_limitation: f64,
    phosphorus_feedback_limitation: f64,
    water_uptake_per_plant_m3_per_step: f64,
    water_uptake_for_gas_step_m3_per_step: f64,
};

pub const Outcome = union(enum) {
    inactive,
    active: Active,
};

/// UPTAKE.F 1662--1740. Evaluates one root/mycorrhizal layer admission and
/// the protein, respiration, nutrient-feedback, and water-uptake controls.
pub fn calculate(
    admission: Admission,
    protein: Protein,
    respiration: Respiration,
    feedback: NutrientFeedback,
    water: Water,
) !Outcome {
    try validateAdmission(admission);
    if (!(admission.layer_thickness_m > admission.minimum_layer_thickness_m and
        admission.root_length_density_m_per_m3 >
            admission.negligible_root_amount and
        admission.root_aqueous_volume_m3 >
            admission.negligible_root_amount and
        admission.soil_water_content_m3_per_m3 >
            admission.negligible_soil_water_content_m3_per_m3))
        return .inactive;
    try validateActive(protein, respiration, feedback, water);

    const protein_concentration =
        if (protein.root_carbon_mass_g_c > admission.negligible_root_amount)
            @min(
                protein.maximum_protein_concentration_g_per_g_c,
                protein.protein_content_g / protein.root_carbon_mass_g_c,
            )
        else
            protein.maximum_protein_concentration_g_per_g_c;
    const relative_protein_capacity =
        if (protein.root_carbon_mass_g_c > admission.negligible_root_amount)
            protein_concentration /
                protein.reference_protein_concentration_g_per_g_c
        else
            1;
    const carbon_limitation =
        if (respiration.unconstrained_uptake_respiration_g_c_per_step >
        admission.negligible_root_amount)
            @max(
                0,
                @min(
                    1,
                    respiration.nonstructural_carbon_pool_g_c *
                        respiration.biological_timestep_h_per_step /
                        respiration.unconstrained_uptake_respiration_g_c_per_step,
                ),
            )
        else
            0;

    var nitrogen_limitation: f64 = 1;
    var phosphorus_limitation: f64 = 1;
    if (feedback.nonstructural_carbon_g_c_per_g_c >
        feedback.significant_concentration_g_per_g_c and
        (feedback.nutrient_uptake_strategy_code != 0 or
            feedback.vascular_carbon_g_c <
                feedback.vascular_fraction_threshold *
                    feedback.total_root_carbon_g_c))
    {
        nitrogen_limitation = @min(
            feedback.nonstructural_carbon_g_c_per_g_c /
                (feedback.nonstructural_carbon_g_c_per_g_c +
                    feedback.nonstructural_nitrogen_g_n_per_g_c /
                        feedback.nitrogen_inhibition_by_nitrogen_g_n_per_g_c),
            feedback.nonstructural_phosphorus_g_p_per_g_c /
                (feedback.nonstructural_phosphorus_g_p_per_g_c +
                    feedback.nonstructural_nitrogen_g_n_per_g_c /
                        feedback.nitrogen_inhibition_by_phosphorus_g_n_per_g_c),
        );
        phosphorus_limitation = @min(
            feedback.nonstructural_carbon_g_c_per_g_c /
                (feedback.nonstructural_carbon_g_c_per_g_c +
                    feedback.nonstructural_phosphorus_g_p_per_g_c /
                        feedback.phosphorus_inhibition_by_phosphorus_g_p_per_g_c),
            feedback.nonstructural_nitrogen_g_n_per_g_c /
                (feedback.nonstructural_nitrogen_g_n_per_g_c +
                    feedback.nonstructural_phosphorus_g_p_per_g_c /
                        feedback.phosphorus_inhibition_by_nitrogen_g_p_per_g_c),
        );
    }
    const uptake_per_plant = @max(
        0,
        -water.root_water_uptake_m3_per_step / water.plant_population,
    );
    const uptake_for_gas_step =
        uptake_per_plant * water.gas_flux_timestep_h_per_step;
    inline for (.{
        protein_concentration,
        relative_protein_capacity,
        carbon_limitation,
        nitrogen_limitation,
        phosphorus_limitation,
        uptake_per_plant,
        uptake_for_gas_step,
    }) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteRootNutrientUptakeControlResult;
    return .{ .active = .{
        .accumulators = .{},
        .protein_concentration_g_per_g_c = protein_concentration,
        .relative_protein_uptake_capacity = relative_protein_capacity,
        .carbon_respiration_limitation = carbon_limitation,
        .nitrogen_feedback_limitation = nitrogen_limitation,
        .phosphorus_feedback_limitation = phosphorus_limitation,
        .water_uptake_per_plant_m3_per_step = uptake_per_plant,
        .water_uptake_for_gas_step_m3_per_step = uptake_for_gas_step,
    } };
}

fn validateAdmission(admission: Admission) !void {
    inline for (@typeInfo(Admission).@"struct".fields) |field|
        if (!std.math.isFinite(@field(admission, field.name)))
            return error.InvalidRootNutrientUptakeAdmission;
    if (admission.layer_thickness_m < 0 or
        admission.minimum_layer_thickness_m < 0 or
        admission.root_length_density_m_per_m3 < 0 or
        admission.root_aqueous_volume_m3 < 0 or
        admission.soil_water_content_m3_per_m3 < 0 or
        admission.negligible_root_amount < 0 or
        admission.negligible_soil_water_content_m3_per_m3 < 0)
        return error.InvalidRootNutrientUptakeAdmission;
}

fn validateActive(
    protein: Protein,
    respiration: Respiration,
    feedback: NutrientFeedback,
    water: Water,
) !void {
    inline for (@typeInfo(Protein).@"struct".fields) |field|
        if (!std.math.isFinite(@field(protein, field.name)))
            return error.InvalidRootNutrientUptakeControlInput;
    inline for (@typeInfo(Respiration).@"struct".fields) |field|
        if (!std.math.isFinite(@field(respiration, field.name)))
            return error.InvalidRootNutrientUptakeControlInput;
    inline for (@typeInfo(NutrientFeedback).@"struct".fields) |field| {
        if (field.type == u8) continue;
        if (!std.math.isFinite(@field(feedback, field.name)))
            return error.InvalidRootNutrientUptakeControlInput;
    }
    inline for (@typeInfo(Water).@"struct".fields) |field|
        if (!std.math.isFinite(@field(water, field.name)))
            return error.InvalidRootNutrientUptakeControlInput;
    if (protein.protein_content_g < 0 or protein.root_carbon_mass_g_c < 0 or
        protein.maximum_protein_concentration_g_per_g_c < 0 or
        protein.reference_protein_concentration_g_per_g_c <= 0 or
        respiration.nonstructural_carbon_pool_g_c < 0 or
        respiration.biological_timestep_h_per_step < 0 or
        feedback.vascular_fraction_threshold < 0 or
        feedback.nitrogen_inhibition_by_nitrogen_g_n_per_g_c <= 0 or
        feedback.nitrogen_inhibition_by_phosphorus_g_n_per_g_c <= 0 or
        feedback.phosphorus_inhibition_by_phosphorus_g_p_per_g_c <= 0 or
        feedback.phosphorus_inhibition_by_nitrogen_g_p_per_g_c <= 0 or
        feedback.significant_concentration_g_per_g_c < 0 or
        water.plant_population <= 0 or
        water.gas_flux_timestep_h_per_step < 0)
        return error.InvalidRootNutrientUptakeControlInput;
}

fn activeAdmission() Admission {
    return .{
        .layer_thickness_m = 0.1,
        .minimum_layer_thickness_m = 0.01,
        .root_length_density_m_per_m3 = 10,
        .root_aqueous_volume_m3 = 0.2,
        .soil_water_content_m3_per_m3 = 0.3,
        .negligible_root_amount = 1e-12,
        .negligible_soil_water_content_m3_per_m3 = 1e-12,
    };
}

fn proteinInputs() Protein {
    return .{
        .protein_content_g = 2,
        .root_carbon_mass_g_c = 20,
        .maximum_protein_concentration_g_per_g_c = 0.08,
        .reference_protein_concentration_g_per_g_c = 0.05,
    };
}

fn respirationInputs() Respiration {
    return .{
        .unconstrained_uptake_respiration_g_c_per_step = 4,
        .nonstructural_carbon_pool_g_c = 1,
        .biological_timestep_h_per_step = 2,
    };
}

fn feedbackInputs() NutrientFeedback {
    return .{
        .nonstructural_carbon_g_c_per_g_c = 0.2,
        .nonstructural_nitrogen_g_n_per_g_c = 0.04,
        .nonstructural_phosphorus_g_p_per_g_c = 0.01,
        .nutrient_uptake_strategy_code = 1,
        .vascular_carbon_g_c = 1,
        .total_root_carbon_g_c = 10,
        .vascular_fraction_threshold = 1e-3,
        .nitrogen_inhibition_by_nitrogen_g_n_per_g_c = 0.02,
        .nitrogen_inhibition_by_phosphorus_g_n_per_g_c = 0.01,
        .phosphorus_inhibition_by_phosphorus_g_p_per_g_c = 0.005,
        .phosphorus_inhibition_by_nitrogen_g_p_per_g_c = 0.002,
        .significant_concentration_g_per_g_c = 1e-12,
    };
}

fn waterInputs() Water {
    return .{
        .root_water_uptake_m3_per_step = -2,
        .plant_population = 4,
        .gas_flux_timestep_h_per_step = 0.25,
    };
}

test "UPTAKE active root nutrient controls preserve source order" {
    const outcome = try calculate(
        activeAdmission(),
        proteinInputs(),
        respirationInputs(),
        feedbackInputs(),
        waterInputs(),
    );
    const active = outcome.active;
    try std.testing.expectEqual(@as(f64, 0.08), active.protein_concentration_g_per_g_c);
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.6),
        active.relative_protein_uptake_capacity,
        3e-16,
    );
    try std.testing.expectEqual(@as(f64, 0.5), active.carbon_respiration_limitation);
    try std.testing.expect(active.nitrogen_feedback_limitation >= 0);
    try std.testing.expect(active.phosphorus_feedback_limitation >= 0);
    try std.testing.expectEqual(@as(f64, 0.5), active.water_uptake_per_plant_m3_per_step);
    try std.testing.expectEqual(@as(f64, 0.125), active.water_uptake_for_gas_step_m3_per_step);
    try std.testing.expectEqual(@as(f64, 0), active.accumulators.oxygen);
}

test "failed strict admission gate returns inactive without downstream division" {
    var admission = activeAdmission();
    admission.layer_thickness_m = admission.minimum_layer_thickness_m;
    var water = waterInputs();
    water.plant_population = 0;
    const outcome = try calculate(
        admission,
        proteinInputs(),
        respirationInputs(),
        feedbackInputs(),
        water,
    );
    try std.testing.expect(outcome == .inactive);
}

test "source protein respiration and feedback fallback branches are retained" {
    var protein = proteinInputs();
    protein.root_carbon_mass_g_c = 0;
    var respiration = respirationInputs();
    respiration.unconstrained_uptake_respiration_g_c_per_step = 0;
    var feedback = feedbackInputs();
    feedback.nonstructural_carbon_g_c_per_g_c = 0;
    const outcome = try calculate(
        activeAdmission(),
        protein,
        respiration,
        feedback,
        waterInputs(),
    );
    const active = outcome.active;
    try std.testing.expectEqual(
        protein.maximum_protein_concentration_g_per_g_c,
        active.protein_concentration_g_per_g_c,
    );
    try std.testing.expectEqual(@as(f64, 1), active.relative_protein_uptake_capacity);
    try std.testing.expectEqual(@as(f64, 0), active.carbon_respiration_limitation);
    try std.testing.expectEqual(@as(f64, 1), active.nitrogen_feedback_limitation);
    try std.testing.expectEqual(@as(f64, 1), active.phosphorus_feedback_limitation);
}
