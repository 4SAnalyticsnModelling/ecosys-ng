//! Tests for `plant_root_metabolism.zig`.
//!
//! Extracted verbatim so the module beside it contains only the model
//! code. Tests that use private declarations of that module stay there,
//! since a sibling file can only reach `pub` declarations.

const NutrientResult = @import("../plant/root/plant_root_nutrient_uptake.zig").Result;
const RootState = @import("../plant/root/plant_root_system.zig").State;
const root_domain_count = @import("../plant/root/plant_root_system.zig").biological_domain_count;
const std = @import("std");
const plant_root_metabolism = @import("../plant/root/plant_root_metabolism.zig");
test "negative primary growth removes concurrent mycorrhizal structure and mobile pools" {
    const kinetics: plant_root_metabolism.RootLitterFractions = .{
        .woody_carbon = .{ 0.1, 0.2, 0.3, 0.4 },
        .woody_nitrogen = .{ 0.1, 0.2, 0.3, 0.4 },
        .woody_phosphorus = .{ 0.1, 0.2, 0.3, 0.4 },
        .nonwoody_carbon = .{ 0.4, 0.3, 0.2, 0.1 },
        .nonwoody_nitrogen = .{ 0.4, 0.3, 0.2, 0.1 },
        .nonwoody_phosphorus = .{ 0.4, 0.3, 0.2, 0.1 },
    };
    const result = try plant_root_metabolism.mycorrhizalLossWithSecondaryRoots(
        -2,
        8,
        20,
        1.0e-12,
        .{
            .structural_carbon_g_c = 4,
            .structural_nitrogen_g_n = 2,
            .structural_phosphorus_g_p = 1,
            .length_m = 12,
            .mobile_carbon_g_c = 5,
            .mobile_nitrogen_g_n = 3,
            .mobile_phosphorus_g_p = 2,
        },
        .{ .{ 0.6, 0.4 }, .{ 0.5, 0.5 }, .{ 0.25, 0.75 } },
        kinetics,
    );
    try std.testing.expectApproxEqAbs(0.25, result.structural_loss_fraction, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.1, result.mobile_loss_fraction, 1.0e-12);
    try std.testing.expectApproxEqAbs(3, result.remaining.structural_carbon_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(9, result.remaining.length_m, 1.0e-12);
    try std.testing.expectApproxEqAbs(4.5, result.remaining.mobile_carbon_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.06, result.litter.woody_carbon_g_c[0], 1.0e-12);
    try std.testing.expectApproxEqAbs(0.36, result.litter.nonwoody_carbon_g_c[0], 1.0e-12);

    const unchanged = try plant_root_metabolism.mycorrhizalLossWithSecondaryRoots(
        0.2,
        0,
        0,
        1.0e-12,
        result.remaining,
        .{ .{ 0.6, 0.4 }, .{ 0.5, 0.5 }, .{ 0.25, 0.75 } },
        kinetics,
    );
    try std.testing.expectEqual(@as(f64, 0), unchanged.structural_loss_fraction);
    try std.testing.expectEqual(result.remaining, unchanged.remaining);
}

test "live GROSUB mycorrhizal loss follows current then upper host-root deficits" {
    var roots = try RootState.init(std.testing.allocator, 1, 2, 1);
    defer roots.deinit();
    var workspace = try plant_root_metabolism.AxisWorkspace.init(std.testing.allocator, 1, 2);
    defer workspace.deinit();
    try workspace.resetAxes(1);
    workspace.primary_deficit_active[0] = true;
    workspace.primary_deficit_absorption[0] = try plant_root_metabolism.absorbPrimaryDeficitFromSecondaryRoots(
        3,
        0,
        0,
        .{ .carbon_g_c = 2, .nitrogen_g_n = 0, .phosphorus_g_p = 0, .length_m = 2 },
        .{ .carbon_g_c = 4, .nitrogen_g_n = 0, .phosphorus_g_p = 0, .length_m = 4 },
    );
    const upper_axis = try roots.layerAxisIndex(0, 1, 0, 0);
    const current_axis = try roots.layerAxisIndex(0, 1, 1, 0);
    const upper_root = try roots.layerIndex(0, 1, 0);
    const current_root = try roots.layerIndex(0, 1, 1);
    roots.axis_secondary_carbon_g[upper_axis] = 8;
    roots.axis_secondary_carbon_g[current_axis] = 4;
    roots.axis_secondary_length_m[upper_axis] = 16;
    roots.axis_secondary_length_m[current_axis] = 8;
    roots.mobile_carbon_g[upper_root] = 8;
    roots.mobile_carbon_g[current_root] = 10;
    const unit = [_]f64{1} ** 4;
    const zero = [_]f64{0} ** 4;
    const litter = try plant_root_metabolism.commitMycorrhizalLossWithSecondaryRoots(
        &roots,
        0,
        1,
        workspace,
        1,
        .{ 10, 8 },
        1.0e-12,
        .{ .{ 0, 1 }, .{ 0, 1 }, .{ 0, 1 } },
        .{
            .woody_carbon = zero,
            .woody_nitrogen = zero,
            .woody_phosphorus = zero,
            .nonwoody_carbon = unit,
            .nonwoody_nitrogen = unit,
            .nonwoody_phosphorus = unit,
        },
    );
    try std.testing.expectEqual(@as(f64, 0), roots.axis_secondary_carbon_g[current_axis]);
    try std.testing.expectEqual(@as(f64, 0), roots.axis_secondary_length_m[current_axis]);
    try std.testing.expectApproxEqAbs(@as(f64, 16.0 / 3.0), roots.axis_secondary_carbon_g[upper_axis], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 32.0 / 3.0), roots.axis_secondary_length_m[upper_axis], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 7), roots.mobile_carbon_g[current_root], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 7), roots.mobile_carbon_g[upper_root], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 7), litter.current.nonwoody_carbon_g_c[0], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 11.0 / 3.0), litter.upper.nonwoody_carbon_g_c[0], 1.0e-12);
}

test "GROSUB root respiration assembly preserves RCO2T and RCO2TM equations" {
    const components: plant_root_metabolism.Components = .{
        .maintenance_demand_g_c = 3,
        .substrate_respiration_actual_g_c = 2,
        .substrate_respiration_oxygen_unlimited_g_c = 4,
        .growth_respiration_actual_g_c = 0.5,
        .growth_respiration_oxygen_unlimited_g_c = 0.8,
        .senescence_respiration_actual_g_c = 0.2,
        .senescence_respiration_oxygen_unlimited_g_c = 0.3,
        .nitrogen_assimilation_respiration_actual_g_c = 0.1,
        .nitrogen_assimilation_respiration_oxygen_unlimited_g_c = 0.15,
    };
    const result = try plant_root_metabolism.assemble(components);
    try std.testing.expectApproxEqAbs(@as(f64, 2.8), result.actual_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 4.25), result.oxygen_unlimited_g_c, 1.0e-12);
    try std.testing.expectEqual(result.actual_g_c, result.carbon_unlimited_g_c);
}

test "GROSUB secondary-root metabolism preserves source equations" {
    const parameters: plant_root_metabolism.SecondaryRootParameters = .{
        .maximum_substrate_respiration_fraction_per_h = 0.015,
        .substrate_respiration_half_saturation_g_c_per_g_c = 0.025,
        .nitrogen_feedback_half_saturation_g_n_per_g_c = 0.1,
        .phosphorus_feedback_half_saturation_g_p_per_g_c = 0.01,
        .maintenance_respiration_g_c_per_g_n_h = 0.010,
        .nitrogen_assimilation_respiration_g_c_per_g_n = 1.70,
        .minimum_carbon_recycling_fraction = 0.167,
        .responsive_carbon_recycling_fraction = 0.333,
        .maximum_nitrogen_recycling_fraction = 0.667,
        .maximum_phosphorus_recycling_fraction = 0.667,
        .storage_exchange_fraction_per_h = 2.5e-5,
        .nonwoody_root_fraction_exponent = 0.167,
        .maintenance_gas_constant_j_per_mol_k = 8.3143,
        .maintenance_enthalpy_j_per_mol_k = 710,
        .maintenance_activation_energy_j_per_mol = 62500,
        .maintenance_low_temperature_inactivation_energy_j_per_mol = 197500,
        .maintenance_normalization_log_intercept = 25.216,
        .maximum_maintenance_temperature_response = 1.0e3,
        .shallow_root_water_response_per_megapascal = 0.05,
        .deep_root_water_response_per_megapascal = 0.10,
        .maintenance_water_response_exponent = 0.25,
        .root_penetration_reference_radius_m = 1.0e-3,
        .acidity_half_effect_hydrogen_activity_mol_per_m3 = 1,
        .maximum_acidity_enhancement = 4,
        .shallow_primary_root_sink_multiplier = 0.25,
        .intermediate_primary_root_sink_multiplier = 1,
        .deep_primary_root_sink_multiplier = 2,
        .deeper_primary_root_sink_multiplier = 4,
        .annual_termination_hours_without_grain_fill = 336,
        .root_protein_carbon_per_nitrogen_g_c_per_g_n = 2.5,
        .root_protein_carbon_per_phosphorus_g_c_per_g_p = 25,
        .nutrient_uptake_respiration_g_c_per_g_element = 0.86,
        .evergreen_leafoff_remobilization_start_fraction = 0.75,
        .deciduous_leafoff_remobilization_start_fraction = 0.5,
        .full_senescence_duration_h = 480,
    };
    const inputs: plant_root_metabolism.SecondaryRootInputs = .{
        .mobile_carbon_g_c = 0.2,
        .nonstructural_nitrogen_g_n = 0.04,
        .nonstructural_phosphorus_g_p = 0.004,
        .root_carbon_g_c = 2,
        .root_nitrogen_g_n = 0.04,
        .root_nitrogen_to_carbon_ratio_g_n_per_g_c = 0.02,
        .root_phosphorus_to_carbon_ratio_g_p_per_g_c = 0.002,
        .root_growth_yield_g_c_per_g_c = 0.8,
        .active_root_fraction = 0.5,
        .biological_timestep_h = 1,
        .substrate_temperature_response = 0.9,
        .maintenance_temperature_response = 0.8,
        .acidity_response = 0.75,
        .substrate_feedback = 0.6,
        .oxygen_limitation = 0.7,
        .substrate_water_response = 0.5,
        .maintenance_water_response = 0.5,
    };
    const result = try plant_root_metabolism.secondaryRootMetabolism(parameters, inputs);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0 / 3.0), result.nutrient_feedback, 1.0e-12);
    const rco2rm = 0.015 * 0.5 * 0.2 * 0.9 * (2.0 / 3.0) * 0.6 * 0.5 * 0.1 / 0.125;
    const rmncr = 0.010 * 0.04 * 0.8 * 0.75 * 0.5;
    try std.testing.expectApproxEqAbs(rco2rm, result.substrate_respiration_oxygen_unlimited_g_c_per_h, 1.0e-12);
    try std.testing.expectApproxEqAbs(rmncr, result.maintenance_respiration_g_c_per_h, 1.0e-12);
    try std.testing.expectApproxEqAbs(1.70 * result.nitrogen_growth_actual_g_n_per_h, result.nitrogen_assimilation_respiration_actual_g_c_per_h, 1.0e-12);
    var singular = inputs;
    singular.root_growth_yield_g_c_per_g_c = 1;
    try std.testing.expectError(error.InvalidSecondaryRootMetabolismInput, plant_root_metabolism.secondaryRootMetabolism(parameters, singular));
}

test "GROSUB secondary-root entry and respiration water selectors preserve source gates" {
    try std.testing.expect(plant_root_metabolism.secondaryRootAxisActive(2, 2, false));
    try std.testing.expect(!plant_root_metabolism.secondaryRootAxisActive(3, 2, false));
    try std.testing.expect(!plant_root_metabolism.secondaryRootAxisActive(2, 2, true));
    try std.testing.expect(plant_root_metabolism.rootRespirationActive(true, true));
    try std.testing.expect(plant_root_metabolism.rootRespirationActive(false, false));
    try std.testing.expect(!plant_root_metabolism.rootRespirationActive(false, true));

    const evergreen_deep = try plant_root_metabolism.sourceRootRespirationWaterResponses(2, 0, 0.25, 0.70);
    try std.testing.expectEqual(@as(f64, 0.25), evergreen_deep.substrate);
    try std.testing.expectEqual(@as(f64, 0.70), evergreen_deep.maintenance);
    const drought_deciduous = try plant_root_metabolism.sourceRootRespirationWaterResponses(2, 2, 0.25, 0.70);
    try std.testing.expectEqual(@as(f64, 0.25), drought_deciduous.substrate);
    try std.testing.expectEqual(@as(f64, 0.25), drought_deciduous.maintenance);
}

test "runtime root metabolism plant parameters reject invalid dimensions and codes" {
    const valid: plant_root_metabolism.RuntimePlantParameters = .{
        .root_profile_type = 2,
        .mycorrhizal_type = 2,
        .growth_habit = 1,
        .leaf_phenology_type = 0,
        .root_growth_yield_g_c_per_g_c = 0.7,
        .root_nitrogen_to_carbon_g_n_per_g_c = 0.03,
        .root_phosphorus_to_carbon_g_p_per_g_c = 0.004,
        .stalk_nitrogen_to_carbon_g_n_per_g_c = 0.01,
        .stalk_phosphorus_to_carbon_g_p_per_g_c = 0.001,
        .primary_root_radius_m = 0.001,
        .secondary_root_radius_m = 0.0002,
        .primary_specific_length_m_per_g_c = 10,
        .secondary_specific_length_m_per_g_c = 100,
        .secondary_root_branching_per_m = 20,
        .shoot_root_equilibration_fraction_per_h = 0.1,
    };
    try valid.validate();
    var invalid = valid;
    invalid.root_profile_type = 4;
    try std.testing.expectError(error.InvalidRootMetabolismPlantCode, invalid.validate());
    invalid = valid;
    invalid.secondary_specific_length_m_per_g_c = 0;
    try std.testing.expectError(error.InvalidRootMetabolismPlantParameter, invalid.validate());
}

test "STARTQ CNRTS and CPRTS yield-scaled ratios bound GROSUB root growth respiration" {
    const respiration_fraction = 0.3;
    const growth_yield = 0.8;
    const active_fraction = 0.5;
    const nitrogen_ratio = 0.02;
    const phosphorus_ratio = 0.002;
    const nitrogen = 0.04;
    const phosphorus = 0.003;
    const result = try plant_root_metabolism.nutrientLimitedRootGrowthRespiration(
        nitrogen,
        phosphorus,
        active_fraction,
        respiration_fraction,
        growth_yield,
        nitrogen_ratio,
        phosphorus_ratio,
    );
    const source_n_limit = nitrogen * respiration_fraction * active_fraction / (nitrogen_ratio * growth_yield);
    const source_p_limit = phosphorus * respiration_fraction * active_fraction / (phosphorus_ratio * growth_yield);
    try std.testing.expect(source_p_limit < source_n_limit);
    try std.testing.expectApproxEqAbs(source_p_limit, result, 1.0e-15);
    // Converting source growth respiration back through DMRTD and DMRT
    // consumes exactly the limiting active phosphorus inventory.
    const structural_growth_g_c = result / respiration_fraction * growth_yield;
    try std.testing.expectApproxEqAbs(phosphorus * active_fraction, structural_growth_g_c * phosphorus_ratio, 1.0e-15);
}

test "GROSUB primary-root bottom cap uses current axis maintenance demand" {
    const result = try plant_root_metabolism.primaryRootMetabolism(plant_root_metabolism.compatibilitySecondaryRootParameters(), .{
        .shared = .{
            .mobile_carbon_g_c = 1,
            .nonstructural_nitrogen_g_n = 0.1,
            .nonstructural_phosphorus_g_p = 0.01,
            .root_carbon_g_c = 2,
            .root_nitrogen_g_n = 0.02,
            .root_nitrogen_to_carbon_ratio_g_n_per_g_c = 0.02,
            .root_phosphorus_to_carbon_ratio_g_p_per_g_c = 0.002,
            .root_growth_yield_g_c_per_g_c = 0.8,
            .active_root_fraction = 1,
            .biological_timestep_h = 1,
            .substrate_temperature_response = 1,
            .maintenance_temperature_response = 1,
            .acidity_response = 1,
            .substrate_feedback = 1,
            .oxygen_limitation = 1,
            .substrate_water_response = 1,
            .maintenance_water_response = 1,
        },
        .primary_tip_at_or_below_profile_bottom = true,
    });
    try std.testing.expectApproxEqAbs(result.maintenance_respiration_g_c_per_h, result.substrate_respiration_oxygen_unlimited_g_c_per_h, 1e-15);
    try std.testing.expectEqual(@as(f64, 0), result.growth_respiration_actual_g_c_per_h);
}

test "GROSUB primary-root plant_root_metabolism.commit updates axis and shared pools atomically" {
    var roots = try RootState.init(std.testing.allocator, 1, 1, 1);
    defer roots.deinit();
    roots.mobile_carbon_g[0] = 10;
    roots.mobile_nitrogen_g[0] = 1;
    roots.mobile_phosphorus_g[0] = 0.1;
    roots.axis_primary_carbon_g[0] = 2;
    roots.axis_primary_nitrogen_g[0] = 0.2;
    roots.axis_primary_phosphorus_g[0] = 0.02;
    roots.axis_primary_length_m[0] = 4;
    roots.axis_depth_m[0] = 1;
    const senescence: plant_root_metabolism.SecondaryRootSenescence = .{
        .respiration_oxygen_unlimited_g_c_per_h = 0.1,
        .respiration_actual_g_c_per_h = 0.1,
        .phenological_senescence_g_c_per_h = 0,
        .senesced_fraction = 0.25,
        .recyclable_carbon_g_c = 0.8,
        .recyclable_nitrogen_g_n = 0.08,
        .recyclable_phosphorus_g_p = 0.008,
    };
    try plant_root_metabolism.commitPrimaryRoot(&roots, 0, 0, 0, .{
        .metabolism = std.mem.zeroes(plant_root_metabolism.SecondaryRootResult),
        .senescence = senescence,
        .primary_specific_length_m_per_g_c = 10,
        .root_extension_water_response = 1,
        .nonwoody_carbon_fraction = 0.6,
        .nonwoody_nitrogen_fraction = 0.7,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 0,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 0,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 10.02), roots.mobile_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.014), roots.mobile_nitrogen_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.102), roots.mobile_phosphorus_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), roots.axis_primary_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), roots.axis_primary_length_m[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), roots.axis_depth_m[0], 1e-12);
    const mobile_before = roots.mobile_carbon_g[0];
    try std.testing.expectError(error.PlantRootIndexOutOfBounds, plant_root_metabolism.commitPrimaryRoot(&roots, 0, 0, 2, .{
        .metabolism = std.mem.zeroes(plant_root_metabolism.SecondaryRootResult),
        .senescence = senescence,
        .primary_specific_length_m_per_g_c = 10,
        .root_extension_water_response = 1,
        .nonwoody_carbon_fraction = 0.6,
        .nonwoody_nitrogen_fraction = 0.7,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 0,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 0,
    }));
    try std.testing.expectEqual(mobile_before, roots.mobile_carbon_g[0]);
}

test "GROSUB primary-root respiration is allocated across traversed layers" {
    var roots = try RootState.init(std.testing.allocator, 1, 3, 1);
    defer roots.deinit();
    const indices = [_]usize{
        try roots.layerIndex(0, 0, 0),
        try roots.layerIndex(0, 0, 1),
        try roots.layerIndex(0, 0, 2),
    };
    var fractions: [3]f64 = undefined;
    try plant_root_metabolism.allocatePrimaryRootRespiration(&roots, &indices, &.{ 0.2, 0.3, 0.1 }, &fractions, 1.1, 0.1, true, .{
        .actual_g_c = 2,
        .oxygen_unlimited_g_c = 3,
        .carbon_unlimited_g_c = 4,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), fractions[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), fractions[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), fractions[2], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), roots.actual_respiration_g_c_per_h[indices[0]], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), roots.actual_respiration_g_c_per_h[indices[2]], 1e-15);
    var actual_total: f64 = 0;
    for (indices) |root| actual_total += roots.actual_respiration_g_c_per_h[root];
    try std.testing.expectApproxEqAbs(@as(f64, 2), actual_total, 1e-15);
}

test "GROSUB primary root depth retracts under negative net growth" {
    const positive = try plant_root_metabolism.primaryRootLengthChange(0.2, 0.2, 4, 1.1, 0.1, 10, 2, 0.5);
    try std.testing.expectApproxEqAbs(0.5, positive, 1e-15);
    const retraction = try plant_root_metabolism.primaryRootLengthChange(0.02, -0.8, 4, 1.1, 0.1, 10, 2, 0.5);
    // Gross extension is 0.05 m; proportional withdrawal is -0.20 m.
    try std.testing.expectApproxEqAbs(-0.15, retraction, 1e-15);
}

test "GROSUB secondary-root recycling and senescence preserve source branches" {
    const parameters = plant_root_metabolism.compatibilitySecondaryRootParameters();
    const recycling = try plant_root_metabolism.secondaryRootRecyclingFractions(true, 0.2, 0.04, 0.004, parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.167 + 0.333 * (2.0 / 3.0)), recycling.carbon, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.667 / 3.0), recycling.nitrogen, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.667 / 3.0), recycling.phosphorus, 1.0e-12);
    const senescence = try plant_root_metabolism.secondaryRootSenescence(.{
        .oxygen_unlimited_substrate_minus_maintenance_g_c_per_h = -0.3,
        .actual_substrate_minus_maintenance_g_c_per_h = -0.4,
        .root_carbon_g_c = 1,
        .root_nitrogen_g_n = 0.02,
        .root_phosphorus_g_p = 0.002,
        .oxygen_limitation = 0.5,
        .phenological_remobilization_enabled = true,
        .root_remobilization_enabled = true,
        .storage_exchange_fraction_per_h = 0.01,
        .remobilization_elapsed_h = 50,
        .full_senescence_h = 100,
        .biological_timestep_h = 1,
        .structural_presence_threshold_g_c = 1.0e-12,
    }, recycling);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), senescence.respiration_oxygen_unlimited_g_c_per_h, 1.0e-12);
    try std.testing.expectApproxEqAbs(senescence.recyclable_carbon_g_c * 0.5 + 0.005, senescence.respiration_actual_g_c_per_h, 1.0e-12);
    try std.testing.expectApproxEqAbs(senescence.respiration_actual_g_c_per_h / senescence.recyclable_carbon_g_c, senescence.senesced_fraction, 1.0e-12);
}

test "GROSUB primary senescence excludes secondary phenological remobilization" {
    const recycling: plant_root_metabolism.RecyclingFractions = .{ .carbon = 0.5, .nitrogen = 0.5, .phosphorus = 0.5 };
    const inputs: plant_root_metabolism.SecondaryRootSenescenceInputs = .{
        .oxygen_unlimited_substrate_minus_maintenance_g_c_per_h = 0,
        .actual_substrate_minus_maintenance_g_c_per_h = 0,
        .root_carbon_g_c = 10,
        .root_nitrogen_g_n = 1,
        .root_phosphorus_g_p = 0.1,
        .oxygen_limitation = 1,
        .phenological_remobilization_enabled = true,
        .root_remobilization_enabled = true,
        .storage_exchange_fraction_per_h = 0.1,
        .remobilization_elapsed_h = 10,
        .full_senescence_h = 10,
        .biological_timestep_h = 1,
        .structural_presence_threshold_g_c = 1e-12,
    };
    const secondary = try plant_root_metabolism.secondaryRootSenescence(inputs, recycling);
    const primary = try plant_root_metabolism.primaryRootSenescence(inputs, recycling);
    try std.testing.expectEqual(@as(f64, 1), secondary.phenological_senescence_g_c_per_h);
    try std.testing.expectEqual(@as(f64, 0), primary.phenological_senescence_g_c_per_h);
    try std.testing.expectEqual(@as(f64, 0), primary.respiration_actual_g_c_per_h);
}

test "GROSUB root wood composition preserves FWODR and weighted growth ratios" {
    const composition = try plant_root_metabolism.rootWoodComposition(true, true, 8, 2, 0.01, 0.03, 0.001, 0.003, 1.0e-12, 0.167);
    const nonwoody = std.math.pow(f64, 0.25, 0.167);
    try std.testing.expectApproxEqAbs(nonwoody, composition.carbon_fraction[1], 1.0e-12);
    try std.testing.expectEqual(composition.carbon_fraction, composition.nitrogen_fraction);
    try std.testing.expectApproxEqAbs((1 - nonwoody) * 0.01 + nonwoody * 0.03, composition.growth_nitrogen_to_carbon_g_n_per_g_c, 1.0e-12);
    const herbaceous = try plant_root_metabolism.rootWoodComposition(false, false, 0, 0, 0.01, 0.03, 0.001, 0.003, 1.0e-12, 0.167);
    try std.testing.expectEqual([2]f64{ 0, 1 }, herbaceous.carbon_fraction);
}

test "GROSUB root environment preserves TFN6 FPH and water responses" {
    const parameters = plant_root_metabolism.compatibilitySecondaryRootParameters();
    const response = try plant_root_metabolism.rootEnvironmentResponses(parameters, 298.15, 0, 7, -1, 0.5, 0, 0.2, 0.5e-3, true);
    const adjusted_temperature_k = 298.15;
    const rtk = 8.3143 * adjusted_temperature_k;
    const expected_temperature = @min(@as(f64, 1.0e3), std.math.exp(25.216 - @as(f64, 62500) / rtk) / (1 + std.math.exp((197500 - 710 * adjusted_temperature_k) / rtk)));
    try std.testing.expectApproxEqAbs(expected_temperature, response.maintenance_temperature, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0001), response.acidity, 1.0e-12);
    try std.testing.expectApproxEqAbs(std.math.exp(@as(f64, -0.05)), response.growth_water, 1.0e-12);
    try std.testing.expectApproxEqAbs(std.math.pow(f64, response.growth_water, 0.25), response.maintenance_water, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), response.scaled_penetration_resistance_megapascal, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.45), response.extension_water, 1.0e-12);
}

test "GROSUB next lower root layer skips thin layers but retains the bottom" {
    const thickness_m = [_]f64{ 0.1, 1e-8, 0.2, 0 };
    try std.testing.expectEqual(@as(usize, 2), try plant_root_metabolism.nextLowerRootLayer(&thickness_m, 0, 1e-6));
    try std.testing.expectEqual(@as(usize, 3), try plant_root_metabolism.nextLowerRootLayer(&thickness_m, 2, 1e-6));
    try std.testing.expectError(error.NoLowerRootLayer, plant_root_metabolism.nextLowerRootLayer(&thickness_m, 3, 1e-6));
    try std.testing.expectError(error.InvalidRootLayerThickness, plant_root_metabolism.nextLowerRootLayer(&.{ 0.1, std.math.nan(f64) }, 0, 1e-6));
}

test "GROSUB root axis sink strengths retain series equation and runtime axes" {
    const parameters = plant_root_metabolism.compatibilitySecondaryRootParameters();
    const first = try plant_root_metabolism.rootAxisSinkStrength(parameters, .{
        .root_profile_type = 2,
        .primary_axis_count_multiplier = 3,
        .primary_root_radius_m = 2e-3,
        .primary_root_depth_from_canopy_m = 0.4,
        .secondary_root_depth_from_canopy_m = 0.3,
        .secondary_axis_count = 8,
        .secondary_root_radius_m = 1e-3,
        .average_secondary_root_length_m = 0.2,
        .primary_biological_domain = true,
    });
    const primary_series = 3 * std.math.pow(f64, 2e-3, 2) / 0.3;
    const secondary_parallel = 8 * std.math.pow(f64, 1e-3, 2) / 0.2;
    try std.testing.expectApproxEqAbs(2 * 3 * std.math.pow(f64, 2e-3, 2) / 0.4, first.primary_m, 1.0e-15);
    try std.testing.expectApproxEqAbs(primary_series * secondary_parallel / (primary_series + secondary_parallel), first.secondary_m, 1.0e-15);
    var strengths = [_]plant_root_metabolism.RootAxisSinkStrength{first} ** 7;
    strengths[6] = .{ .primary_m = 0, .secondary_m = first.secondary_m };
    var primary_fractions: [7]f64 = undefined;
    var secondary_fractions: [7]f64 = undefined;
    const total = try plant_root_metabolism.normalizeRootAxisSinkFractions(&strengths, &primary_fractions, &secondary_fractions, 1e-20);
    try std.testing.expect(total > 0);
    var fraction_sum: f64 = 0;
    for (primary_fractions, secondary_fractions) |primary, secondary| fraction_sum += primary + secondary;
    try std.testing.expectApproxEqAbs(@as(f64, 1), fraction_sum, 1.0e-12);
}

test "GROSUB source sink comparator gates primary tips and uses rooted midpoint depth" {
    const parameters = plant_root_metabolism.compatibilitySecondaryRootParameters();
    const shared: plant_root_metabolism.SourceOrderRootAxisSinkInputs = .{
        .root_profile_type = 2,
        .primary_axis_count_multiplier = 3,
        .primary_root_radius_m = 2e-3,
        .primary_root_depth_from_surface_m = 0.45,
        .layer_top_depth_m = 0.2,
        .layer_thickness_m = 0.2,
        .secondary_root_origin_offset_m = 0.05,
        .seeding_depth_m = 0.1,
        .hypocotyledon_height_m = 0.02,
        .canopy_height_m = 0.3,
        .secondary_axis_count = 8,
        .secondary_root_radius_m = 1e-3,
        .average_secondary_root_length_m = 0.2,
        .negligible_sink_m = 1e-20,
        .primary_biological_domain = true,
    };
    const outside_tip_layer = try plant_root_metabolism.sourceOrderRootAxisSinkStrength(parameters, shared);
    try std.testing.expectEqual(@as(f64, 0), outside_tip_layer.primary_m);
    const rooted_length_m = 0.2;
    const secondary_depth_m = 0.2 + 0.5 * rooted_length_m + 0.3;
    const primary_series = 3 * std.math.pow(f64, 2e-3, 2) / secondary_depth_m;
    const secondary_parallel = 8 * std.math.pow(f64, 1e-3, 2) / 0.2;
    try std.testing.expectApproxEqAbs(primary_series * secondary_parallel / (primary_series + secondary_parallel), outside_tip_layer.secondary_m, 1e-15);

    var mycorrhiza = shared;
    mycorrhiza.primary_biological_domain = false;
    const mycorrhizal = try plant_root_metabolism.sourceOrderRootAxisSinkStrength(parameters, mycorrhiza);
    try std.testing.expectEqual(@as(f64, 0), mycorrhizal.primary_m);
    try std.testing.expectApproxEqAbs(secondary_parallel, mycorrhizal.secondary_m, 1e-15);
}

test "root metabolism grid workspace is runtime sized and cell independent" {
    var workspace = try plant_root_metabolism.GridWorkspace.init(std.testing.allocator, 4, 19, 7);
    defer workspace.deinit();
    try std.testing.expectEqual(@as(usize, 4), workspace.per_cell.len);
    for (workspace.per_cell) |cell| {
        try std.testing.expectEqual(@as(usize, 19), cell.sink_strengths.len);
        try std.testing.expectEqual(@as(usize, 7), cell.primary_respiration_allocation_fractions.len);
    }
    workspace.per_cell[0].primary_active[18] = true;
    try workspace.per_cell[0].beginPlantHour(19);
    try workspace.per_cell[0].markPrimaryProcessed(0, 18);
    try workspace.per_cell[0].resetAxes(19);
    try std.testing.expect(!workspace.per_cell[0].primary_active[18]);
    try std.testing.expect(try workspace.per_cell[0].primaryWasProcessed(0, 18));
    try workspace.per_cell[0].beginPlantHour(19);
    try std.testing.expect(!try workspace.per_cell[0].primaryWasProcessed(0, 18));
    try std.testing.expect(workspace.per_cell[0].sink_strengths.ptr != workspace.per_cell[1].sink_strengths.ptr);
    try std.testing.expectError(error.RootMetabolismWorkspaceCapacityExceeded, workspace.per_cell[0].resetAxes(20));
}

test "staged primary and secondary axes plant_root_metabolism.commit one shared mobile-pool transaction" {
    var roots = try RootState.init(std.testing.allocator, 1, 1, 2);
    defer roots.deinit();
    var workspace = try plant_root_metabolism.AxisWorkspace.init(std.testing.allocator, 2, 1);
    defer workspace.deinit();
    try workspace.resetAxes(2);
    roots.mobile_carbon_g[0] = 10;
    roots.mobile_nitrogen_g[0] = 1;
    roots.mobile_phosphorus_g[0] = 0.1;
    roots.axis_primary_carbon_g[0] = 1;
    roots.axis_secondary_carbon_g[1] = 1;
    workspace.primary_active[0] = true;
    workspace.secondary_active[1] = true;
    workspace.primary_metabolism[0].root_growth_actual_g_c_per_h = 0.2;
    workspace.primary_metabolism[0].growth_and_respiration_carbon_actual_g_c_per_h = 0.25;
    workspace.primary_metabolism[0].nitrogen_growth_actual_g_n_per_h = 0.02;
    workspace.primary_metabolism[0].phosphorus_growth_actual_g_p_per_h = 0.002;
    workspace.secondary_metabolism[1].root_growth_actual_g_c_per_h = 0.3;
    workspace.secondary_metabolism[1].growth_and_respiration_carbon_actual_g_c_per_h = 0.375;
    workspace.secondary_metabolism[1].nitrogen_growth_actual_g_n_per_h = 0.03;
    workspace.secondary_metabolism[1].phosphorus_growth_actual_g_p_per_h = 0.003;
    const parameters: plant_root_metabolism.StagedLayerCommitParameters = .{
        .primary_specific_length_m_per_g_c = 10,
        .secondary_specific_length_m_per_g_c = 20,
        .plant_population_count = 1,
        .seeding_depth_m = 0,
        .current_layer_bottom_depth_m = 2,
        .next_layer_thickness_m = 0,
        .extension_presence_threshold_m = 0,
        .root_extension_water_response = 0.5,
        .nonwoody_carbon_fraction = 1,
        .nonwoody_nitrogen_fraction = 1,
        .nonwoody_phosphorus_fraction = 1,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 0,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 0,
    };
    try plant_root_metabolism.commitStagedLayerAxes(&roots, 0, 0, 0, &workspace, 2, parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 9.375), roots.mobile_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.95), roots.mobile_nitrogen_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.095), roots.mobile_phosphorus_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), roots.axis_primary_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.3), roots.axis_secondary_carbon_g[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), roots.axis_primary_length_m[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 3), roots.axis_secondary_length_m[1], 1e-15);

    const mobile_before = roots.mobile_carbon_g[0];
    const primary_before = roots.axis_primary_carbon_g[0];
    try workspace.resetAxes(2);
    workspace.primary_active[0] = true;
    workspace.primary_metabolism[0].growth_and_respiration_carbon_actual_g_c_per_h = 20;
    try std.testing.expectError(error.StagedRootCommitWouldOverdrawPool, plant_root_metabolism.commitStagedLayerAxes(&roots, 0, 0, 0, &workspace, 2, parameters));
    try std.testing.expectEqual(mobile_before, roots.mobile_carbon_g[0]);
    try std.testing.expectEqual(primary_before, roots.axis_primary_carbon_g[0]);
}

test "GROSUB primary crossing routes growth and mobile pools to the next runtime layer" {
    const placement = try plant_root_metabolism.primaryRootExtensionPlacement(2, 0.9, 1.0, 0.5);
    try std.testing.expect(placement.crosses_into_next_layer);
    try std.testing.expectEqual(@as(f64, 0.5), placement.extension_m);

    var roots = try RootState.init(std.testing.allocator, 1, 2, 1);
    defer roots.deinit();
    var workspace = try plant_root_metabolism.AxisWorkspace.init(std.testing.allocator, 1, 2);
    defer workspace.deinit();
    try workspace.resetAxes(1);
    const current_root = try roots.layerIndex(0, 0, 0);
    const next_root = try roots.layerIndex(0, 0, 1);
    const current_axis = try roots.layerAxisIndex(0, 0, 0, 0);
    const next_axis = try roots.layerAxisIndex(0, 0, 1, 0);
    roots.mobile_carbon_g[current_root] = 8;
    roots.mobile_nitrogen_g[current_root] = 1;
    roots.mobile_phosphorus_g[current_root] = 0.1;
    roots.total_water_potential_megapascal[current_root] = -0.4;
    roots.osmotic_water_potential_megapascal[current_root] = -0.8;
    roots.turgor_water_potential_megapascal[current_root] = 0.4;
    roots.primary_radius_m[current_root] = 0.001;
    roots.axis_primary_carbon_g[current_axis] = 1;
    roots.axis_primary_nitrogen_g[current_axis] = 0.1;
    roots.axis_primary_phosphorus_g[current_axis] = 0.01;
    roots.axis_depth_m[try roots.axisIndex(0, 0, 0)] = 0.9;
    workspace.primary_active[0] = true;
    workspace.primary_sink_fractions[0] = 0.25;
    workspace.primary_metabolism[0].root_growth_actual_g_c_per_h = 0.2;
    workspace.primary_metabolism[0].nitrogen_growth_actual_g_n_per_h = 0.02;
    workspace.primary_metabolism[0].phosphorus_growth_actual_g_p_per_h = 0.002;
    try plant_root_metabolism.commitStagedLayerAxes(&roots, 0, 0, 0, &workspace, 1, .{
        .primary_specific_length_m_per_g_c = 10,
        .secondary_specific_length_m_per_g_c = 20,
        .plant_population_count = 1,
        .seeding_depth_m = 0,
        .current_layer_bottom_depth_m = 1,
        .next_layer_thickness_m = 0.5,
        .extension_presence_threshold_m = 0,
        .root_extension_water_response = 1,
        .nonwoody_carbon_fraction = 1,
        .nonwoody_nitrogen_fraction = 1,
        .nonwoody_phosphorus_fraction = 1,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 0,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 0,
    });
    try std.testing.expectEqual(@as(f64, 1), roots.axis_primary_carbon_g[current_axis]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), roots.axis_primary_carbon_g[next_axis], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), roots.axis_primary_length_m[next_axis], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.4), roots.axis_depth_m[try roots.axisIndex(0, 0, 0)], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 6), roots.mobile_carbon_g[current_root], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2), roots.mobile_carbon_g[next_root], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.735), roots.mobile_nitrogen_g[current_root], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.245), roots.mobile_nitrogen_g[next_root], 1e-15);
    try std.testing.expectEqual(roots.total_water_potential_megapascal[current_root], roots.total_water_potential_megapascal[next_root]);
    try std.testing.expectEqual(roots.primary_radius_m[current_root], roots.primary_radius_m[next_root]);
}

test "GROSUB primary crossing requires extension above ZEROP" {
    const below = try plant_root_metabolism.sourceOrderPrimaryRootExtensionPlacement(5e-7, 0.9999998, 1, 0.5, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 5e-7), below.extension_m, 1e-18);
    try std.testing.expect(!below.crosses_into_next_layer);
    const above = try plant_root_metabolism.sourceOrderPrimaryRootExtensionPlacement(2e-6, 0.9999998, 1, 0.5, 1e-6);
    try std.testing.expect(above.crosses_into_next_layer);
}

test "GROSUB negative primary growth consumes secondary roots in tip then upper layer" {
    const result = try plant_root_metabolism.absorbPrimaryDeficitFromSecondaryRoots(
        3.5,
        0.35,
        0.035,
        .{ .carbon_g_c = 2, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.02, .length_m = 4 },
        .{ .carbon_g_c = 3, .nitrogen_g_n = 0.3, .phosphorus_g_p = 0.03, .length_m = 6 },
    );
    try std.testing.expectEqual(@as(f64, 0), result.current.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), result.current.length_m);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), result.upper.carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 3), result.upper.length_m, 1e-15);
    try std.testing.expectEqual(@as(f64, 0), result.residual_carbon_deficit_g_c);
    try std.testing.expectEqual(@as(f64, 0), result.residual_nitrogen_deficit_g_n);
    try std.testing.expectEqual(@as(f64, 0), result.residual_phosphorus_deficit_g_p);

    const exhausted = try plant_root_metabolism.absorbPrimaryDeficitFromSecondaryRoots(
        7,
        0.7,
        0.07,
        .{ .carbon_g_c = 2, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.02, .length_m = 4 },
        .{ .carbon_g_c = 3, .nitrogen_g_n = 0.3, .phosphorus_g_p = 0.03, .length_m = 6 },
    );
    try std.testing.expectEqual(@as(f64, 2), exhausted.residual_carbon_deficit_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), exhausted.residual_nitrogen_deficit_g_n, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.02), exhausted.residual_phosphorus_deficit_g_p, 1e-15);
}

test "live staged GROSUB plant_root_metabolism.commit absorbs primary senescence from secondary layers first" {
    var roots = try RootState.init(std.testing.allocator, 1, 2, 1);
    defer roots.deinit();
    var workspace = try plant_root_metabolism.AxisWorkspace.init(std.testing.allocator, 1, 2);
    defer workspace.deinit();
    try workspace.resetAxes(1);
    const upper_axis = try roots.layerAxisIndex(0, 0, 0, 0);
    const tip_axis = try roots.layerAxisIndex(0, 0, 1, 0);
    roots.axis_primary_carbon_g[tip_axis] = 2;
    roots.axis_primary_nitrogen_g[tip_axis] = 0.2;
    roots.axis_primary_phosphorus_g[tip_axis] = 0.02;
    roots.axis_secondary_carbon_g[tip_axis] = 0.6;
    roots.axis_secondary_nitrogen_g[tip_axis] = 0.06;
    roots.axis_secondary_phosphorus_g[tip_axis] = 0.006;
    roots.axis_secondary_length_m[tip_axis] = 3;
    roots.axis_secondary_carbon_g[upper_axis] = 0.6;
    roots.axis_secondary_nitrogen_g[upper_axis] = 0.06;
    roots.axis_secondary_phosphorus_g[upper_axis] = 0.006;
    roots.axis_secondary_length_m[upper_axis] = 3;
    roots.axis_depth_m[try roots.axisIndex(0, 0, 0)] = 1.5;
    workspace.primary_active[0] = true;
    workspace.primary_senescence[0].senesced_fraction = 0.5;
    try plant_root_metabolism.commitStagedLayerAxes(&roots, 0, 0, 1, &workspace, 1, .{
        .primary_specific_length_m_per_g_c = 10,
        .secondary_specific_length_m_per_g_c = 20,
        .plant_population_count = 1,
        .seeding_depth_m = 0,
        .current_layer_bottom_depth_m = 2,
        .next_layer_thickness_m = 0,
        .extension_presence_threshold_m = 0,
        .root_extension_water_response = 1,
        .nonwoody_carbon_fraction = 1,
        .nonwoody_nitrogen_fraction = 1,
        .nonwoody_phosphorus_fraction = 1,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 0,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 0,
    });
    try std.testing.expectEqual(@as(f64, 2), roots.axis_primary_carbon_g[tip_axis]);
    try std.testing.expectEqual(@as(f64, 0), roots.axis_secondary_carbon_g[tip_axis]);
    try std.testing.expectEqual(@as(f64, 0), roots.axis_secondary_length_m[tip_axis]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), roots.axis_secondary_carbon_g[upper_axis], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), roots.axis_secondary_length_m[upper_axis], 1e-15);
    try std.testing.expect(workspace.primary_deficit_active[0]);
}

test "STOMATE annual termination feedback is shared with root metabolism" {
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), try plant_root_metabolism.annualTerminationFeedback(0, 168, 336), 1.0e-12);
    try std.testing.expectEqual(@as(f64, 1), try plant_root_metabolism.annualTerminationFeedback(1, 168, 336));
    try std.testing.expectEqual(@as(f64, 0), try plant_root_metabolism.annualTerminationFeedback(0, 400, 336));
}

test "GROSUB secondary-root litter allocation and plant_root_metabolism.commit are atomic" {
    const senescence: plant_root_metabolism.SecondaryRootSenescence = .{
        .respiration_oxygen_unlimited_g_c_per_h = 0.1,
        .respiration_actual_g_c_per_h = 0.1,
        .phenological_senescence_g_c_per_h = 0,
        .senesced_fraction = 0.25,
        .recyclable_carbon_g_c = 0.8,
        .recyclable_nitrogen_g_n = 0.08,
        .recyclable_phosphorus_g_p = 0.008,
    };
    const quarter = [_]f64{0.25} ** 4;
    const litter = try plant_root_metabolism.secondaryRootLitter(senescence, 2, 0.2, 0.02, .{ 0.4, 0.6 }, .{ 0.3, 0.7 }, .{ 0.2, 0.8 }, .{
        .woody_carbon = quarter,
        .woody_nitrogen = quarter,
        .woody_phosphorus = quarter,
        .nonwoody_carbon = quarter,
        .nonwoody_nitrogen = quarter,
        .nonwoody_phosphorus = quarter,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), litter.woody_carbon_g_c[0], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.045), litter.nonwoody_carbon_g_c[0], 1.0e-12);

    var roots = try RootState.init(std.testing.allocator, 1, 1, 1);
    defer roots.deinit();
    roots.mobile_carbon_g[0] = 10;
    roots.mobile_nitrogen_g[0] = 1;
    roots.mobile_phosphorus_g[0] = 0.1;
    roots.axis_secondary_carbon_g[0] = 2;
    roots.axis_secondary_nitrogen_g[0] = 0.2;
    roots.axis_secondary_phosphorus_g[0] = 0.02;
    roots.axis_secondary_length_m[0] = 4;
    const metabolism = try plant_root_metabolism.secondaryRootMetabolism(plant_root_metabolism.compatibilitySecondaryRootParameters(), .{
        .mobile_carbon_g_c = 0.2,
        .nonstructural_nitrogen_g_n = 0.04,
        .nonstructural_phosphorus_g_p = 0.004,
        .root_carbon_g_c = 2,
        .root_nitrogen_g_n = 0.04,
        .root_nitrogen_to_carbon_ratio_g_n_per_g_c = 0.02,
        .root_phosphorus_to_carbon_ratio_g_p_per_g_c = 0.002,
        .root_growth_yield_g_c_per_g_c = 0.8,
        .active_root_fraction = 0.5,
        .biological_timestep_h = 1,
        .substrate_temperature_response = 0.9,
        .maintenance_temperature_response = 0.8,
        .acidity_response = 0.75,
        .substrate_feedback = 0.6,
        .oxygen_limitation = 0.7,
        .substrate_water_response = 0.5,
        .maintenance_water_response = 0.5,
    });
    const commit_inputs: plant_root_metabolism.SecondaryRootCommitInputs = .{
        .metabolism = metabolism,
        .senescence = senescence,
        .root_specific_length_m_per_g_c = 10,
        .root_extension_water_response = 0.8,
        .nonwoody_carbon_fraction = 0.6,
        .nonwoody_nitrogen_fraction = 0.7,
        .nonwoody_phosphorus_fraction = 0.8,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 2,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 20,
    };
    try plant_root_metabolism.commitSecondaryRoot(&roots, 0, 0, commit_inputs);
    try std.testing.expect(roots.axis_secondary_carbon_g[0] < 2);
    try std.testing.expect(roots.actual_respiration_g_c_per_h[0] > 0);
    const mobile_before_failure = roots.mobile_carbon_g[0];
    roots.mobile_carbon_g[0] = 0;
    var failing_inputs = commit_inputs;
    failing_inputs.senescence.recyclable_carbon_g_c = 0;
    failing_inputs.senescence.recyclable_nitrogen_g_n = 0;
    failing_inputs.senescence.recyclable_phosphorus_g_p = 0;
    try std.testing.expectError(error.SecondaryRootCommitWouldOverdrawPool, plant_root_metabolism.commitSecondaryRoot(&roots, 0, 0, failing_inputs));
    try std.testing.expectEqual(@as(f64, 0), roots.mobile_carbon_g[0]);
    try std.testing.expect(mobile_before_failure > 0);
}

test "GROSUB nutrient uptake respiration retains 0.86 coefficient and three limits" {
    const result = NutrientResult{ .demand_g_element = 1, .uptake_g_element = 0.5, .oxygen_unlimited_uptake_g_element = 0.75, .carbon_unlimited_uptake_g_element = 1, .available_g_element = 2 };
    const respiration = try plant_root_metabolism.nutrientUptakeRespiration(&([_]NutrientResult{result} ** 8), 0.86);
    try std.testing.expectApproxEqAbs(@as(f64, 3.44), respiration.actual_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5.16), respiration.oxygen_unlimited_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 6.88), respiration.carbon_unlimited_g_c, 1.0e-12);
}

test "root respiration plant_root_metabolism.commit is conservative and rollback safe" {
    var roots = try RootState.init(std.testing.allocator, 1, 1, 1);
    defer roots.deinit();
    roots.mobile_carbon_g[0] = 5;
    try plant_root_metabolism.commit(&roots, 0, .{ .actual_g_c = 2, .oxygen_unlimited_g_c = 3, .carbon_unlimited_g_c = 4 });
    try std.testing.expectEqual(@as(f64, 3), roots.mobile_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 2), roots.actual_respiration_g_c_per_h[0]);
    try std.testing.expectError(error.InsufficientRootMobileCarbonForRespiration, plant_root_metabolism.commit(&roots, 0, .{ .actual_g_c = 4, .oxygen_unlimited_g_c = 4, .carbon_unlimited_g_c = 4 }));
    try std.testing.expectEqual(@as(f64, 3), roots.mobile_carbon_g[0]);
}

test "GROSUB primary-root axis scaling retains the larger carbon basis" {
    const decay_limited = try plant_root_metabolism.primaryRootAxisScaling(4, 2, 1, 1);
    const expected_retained = 0.999992087 * 4.0;
    try std.testing.expectApproxEqAbs(expected_retained, decay_limited.retained_root_carbon_g_c_per_plant, 1.0e-15);
    try std.testing.expectApproxEqAbs(std.math.pow(f64, expected_retained, 0.667), decay_limited.primary_axis_count_multiplier, 1.0e-15);

    const biomass_limited = try plant_root_metabolism.primaryRootAxisScaling(1, 18, 3, 1);
    try std.testing.expectEqual(@as(f64, 6), biomass_limited.retained_root_carbon_g_c_per_plant);
    try std.testing.expectApproxEqAbs(std.math.pow(f64, 6, 0.667) * 3, biomass_limited.primary_axis_count_multiplier, 1.0e-15);
}

test "GROSUB primary-root axis scaling clears state without population" {
    const result = try plant_root_metabolism.primaryRootAxisScaling(12, 30, 0, 1);
    try std.testing.expectEqual(@as(f64, 0), result.retained_root_carbon_g_c_per_plant);
    try std.testing.expectEqual(@as(f64, 0), result.primary_axis_count_multiplier);
}

test "GROSUB primary-root axis scaling rejects invalid state" {
    try std.testing.expectError(error.InvalidPrimaryRootAxisScalingInput, plant_root_metabolism.primaryRootAxisScaling(-1, 1, 1, 1));
    try std.testing.expectError(error.InvalidPrimaryRootAxisScalingInput, plant_root_metabolism.primaryRootAxisScaling(1, 1, 1, 0));
    try std.testing.expectError(error.InvalidPrimaryRootAxisScalingInput, plant_root_metabolism.primaryRootAxisScaling(1, std.math.nan(f64), 1, 1));
}
