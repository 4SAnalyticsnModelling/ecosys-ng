//! Tests for `photosynthesis.zig`.
//!
//! Extracted verbatim so the module beside it contains only the model
//! code. Tests that use private declarations of that module stay there,
//! since a sibling file can only reach `pub` declarations.

const branch_organ_growth_publication = @import("../plant/growth/branch_organ_growth_publication.zig");
const c4_leaf_nonstructural_carbon_senescence = @import("../canopy/leaf/c4_nonstructural_carbon_senescence.zig");
const c4_mesophyll_bundle_exchange = @import("../canopy/photosynthesis/c4_mesophyll_bundle_exchange.zig");
const internode_senescence_publication = @import("../canopy/sheath/internode_senescence_publication.zig");
const leaf_node_growth_publication = @import("../canopy/leaf/node_growth_publication.zig");
const node_senescence_cascade_progress = @import("../plant/growth/node_senescence_cascade_progress.zig");
const node_senescence_remobilization_request = @import("../plant/growth/node_senescence_remobilization_request.zig");
const perennial_stalk_senescence_setup = @import("../plant/growth/perennial_stalk_senescence_setup.zig");
const reserve_maintenance_respiration = @import("../plant/growth/reserve_maintenance_respiration.zig");
const residual_stalk_senescence_publication = @import("../plant/growth/residual_stalk_senescence_publication.zig");
const residual_stalk_senescence_request = @import("../plant/growth/residual_stalk_senescence_request.zig");
const shoot_recycling_fraction = @import("../plant/growth/shoot_recycling_fraction.zig");
const shoot_total_senescence_setup = @import("../plant/growth/shoot_total_senescence_setup.zig");
const std = @import("std");
const canopy_photosynthesis = @import("../canopy/photosynthesis/photosynthesis.zig");
test "reseed reconstruction preserves seed and standing-dead inventories" {
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{0});
    defer state.deinit();
    state.plant_seed_storage_carbon_g[0] = 7;
    state.plant_seed_storage_nitrogen_g[0] = 0.7;
    state.plant_seed_storage_phosphorus_g[0] = 0.07;
    state.plant_standing_dead_carbon_g[0] = 10;
    state.plant_standing_dead_nitrogen_g[0] = 1;
    state.plant_standing_dead_phosphorus_g[0] = 0.1;
    state.plant_standing_dead_height_m[0] = 1.5;
    state.plant_standing_dead_aerodynamic_temperature_k[0] = 280;
    state.plant_standing_dead_aerodynamic_vapor_pressure_kpa[0] = 0.8;
    state.plant_standing_dead_surface_temperature_k[0] = 281;
    state.plant_standing_dead_carbon_by_kinetic_g[0..4].* = .{ 1, 2, 3, 4 };
    state.plant_standing_dead_nitrogen_by_kinetic_g[0..4].* = .{ 0.1, 0.2, 0.3, 0.4 };
    state.plant_standing_dead_phosphorus_by_kinetic_g[0..4].* = .{ 0.01, 0.02, 0.03, 0.04 };
    const carried = try canopy_photosynthesis.capturePersistentReseedInventories(&state, 0);
    try state.clearPlantForReconstruction(0);
    try canopy_photosynthesis.restorePersistentReseedInventories(&state, 0, carried);
    try std.testing.expectApproxEqAbs(@as(f64, 7), state.plant_seed_storage_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 10), state.plant_standing_dead_carbon_g[0], 1e-12);
    try std.testing.expectEqual([4]f64{ 1, 2, 3, 4 }, state.plant_standing_dead_carbon_by_kinetic_g[0..4].*);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), state.plant_standing_dead_height_m[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 281), state.plant_standing_dead_surface_temperature_k[0], 1e-12);
}

test "GROSUB zero-area leaf still advances canopy height to capped leaf base" {
    var area = [_]f64{ 9, 9 };
    var carbon = [_]f64{ 9, 9 };
    var nitrogen = [_]f64{ 9, 9 };
    var phosphorus = [_]f64{ 9, 9 };
    const height_m = try canopy_photosynthesis.allocateLeafAcrossCanopyLayers(0, 0, 0, 0, 100, 2, 0.8, 0.4, 1, &.{ 0, 0.5, 1.01 }, &.{ 0.2, 0.5, 0.8, 1 }, &.{ 0.1, 0.2, 0.3, 0.4 }, .{ .area_m2 = &area, .carbon_g = &carbon, .nitrogen_g = &nitrogen, .phosphorus_g = &phosphorus });
    try std.testing.expectEqual(@as(f64, 1.01), height_m);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0 }, &area);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0 }, &carbon);
}

test "production canopy water response retains source values above one" {
    const response = try canopy_photosynthesis.canopyWaterGrowthResponse(false, 0.5, -1.0, 2.0, 2.0);
    try std.testing.expectApproxEqAbs(std.math.exp(1.0), response.stomatal_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(std.math.exp(0.2), response.growth_fraction, 1.0e-15);
    try std.testing.expect(response.water_potential_expansion_fraction > 1);
}

test "canopy topology accepts arbitrary runtime branch node and sample counts" {
    const branch_counts = [_]usize{ 2, 1, 0, 3, 1, 2 };
    const node_counts = [_]usize{ 1, 3, 2, 0, 1, 2, 2, 1, 2 };
    const sample_counts = [_]usize{ 16, 3, 8, 1, 4, 2, 9, 5, 1, 7, 3, 2, 6, 1 };
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 2, 3, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    try std.testing.expectEqual(6, state.plant_carboxylation_umol_per_s.len);
    try std.testing.expectEqual(9, state.branch_c3_feedback_fraction.len);
    try std.testing.expectEqual(14, state.node_co2_limited_carboxylation_umol_per_m2_s.len);
    try std.testing.expectEqual(68, state.sample_carboxylation_umol_per_s.len);
    try std.testing.expectEqual(canopy_photosynthesis.Range{ .first = 3, .end = 3 }, try state.branchRange(2));
    try std.testing.expectEqual(canopy_photosynthesis.Range{ .first = 7, .end = 9 }, try state.nodeRange(5));
    try state.validateFinite();
}

test "replant canopy reconstruction clears all compact domains for only one plant" {
    const branch_counts = [_]usize{ 1, 2 };
    const node_counts = [_]usize{ 1, 1, 2 };
    const sample_counts = [_]usize{ 2, 1, 2, 1 };
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 2, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    @memset(state.plant_mobile_carbon_g, 1);
    @memset(state.plant_standing_dead_carbon_by_kinetic_g, 2);
    @memset(state.branch_leaf_carbon_g, 3);
    @memset(state.branch_salt_content_by_species_mol, 4);
    @memset(state.node_leaf_carbon_g, 5);
    @memset(state.sample_leaf_carbon_g, 6);

    try state.clearPlantForReconstruction(1);

    try std.testing.expectEqual(@as(f64, 1), state.plant_mobile_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.plant_mobile_carbon_g[1]);
    try std.testing.expectEqual(@as(f64, 2), state.plant_standing_dead_carbon_by_kinetic_g[3]);
    for (state.plant_standing_dead_carbon_by_kinetic_g[4..8]) |value| try std.testing.expectEqual(@as(f64, 0), value);
    try std.testing.expectEqual(@as(f64, 3), state.branch_leaf_carbon_g[0]);
    for (state.branch_leaf_carbon_g[1..3]) |value| try std.testing.expectEqual(@as(f64, 0), value);
    try std.testing.expectEqual(@as(f64, 4), state.branch_salt_content_by_species_mol[7]);
    for (state.branch_salt_content_by_species_mol[8..24]) |value| try std.testing.expectEqual(@as(f64, 0), value);
    try std.testing.expectEqual(@as(f64, 5), state.node_leaf_carbon_g[0]);
    for (state.node_leaf_carbon_g[1..4]) |value| try std.testing.expectEqual(@as(f64, 0), value);
    try std.testing.expectEqual(@as(f64, 6), state.sample_leaf_carbon_g[1]);
    for (state.sample_leaf_carbon_g[2..6]) |value| try std.testing.expectEqual(@as(f64, 0), value);
}

test "replant topology compaction retains one branch and node without changing neighbors" {
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 2, &.{ 1, 1 }, &.{ 1, 1 }, &.{ 2, 3 });
    defer state.deinit();
    _ = try state.appendNode(0, 4);
    _ = try state.appendBranch(0, &.{5});
    @memset(state.plant_mobile_carbon_g, 7);
    state.branch_leaf_carbon_g[state.branch_leaf_carbon_g.len - 1] = 11;
    const neighbor_branch_before = state.branch_leaf_carbon_g[(try state.branchRange(1)).first];
    const neighbor_sample_count = (try state.sampleRange((try state.nodeRange((try state.branchRange(1)).first)).first)).end -
        (try state.sampleRange((try state.nodeRange((try state.branchRange(1)).first)).first)).first;

    try state.compactPlantToInitialTopology(0);

    const first = try state.branchRange(0);
    const neighbor = try state.branchRange(1);
    try std.testing.expectEqual(@as(usize, 1), first.end - first.first);
    try std.testing.expectEqual(@as(usize, 1), (try state.nodeRange(first.first)).end - (try state.nodeRange(first.first)).first);
    try std.testing.expectEqual(@as(usize, 1), neighbor.end - neighbor.first);
    try std.testing.expectEqual(neighbor_sample_count, (try state.sampleRange((try state.nodeRange(neighbor.first)).first)).end - (try state.sampleRange((try state.nodeRange(neighbor.first)).first)).first);
    try std.testing.expectEqual(neighbor_branch_before, state.branch_leaf_carbon_g[neighbor.first]);
    try std.testing.expectEqual(@as(f64, 7), state.plant_mobile_carbon_g[1]);
}

test "GROSUB runtime canopy-layer leaf allocation conserves area and C N P" {
    var area = [_]f64{0} ** 3;
    var carbon = [_]f64{0} ** 3;
    var nitrogen = [_]f64{0} ** 3;
    var phosphorus = [_]f64{0} ** 3;
    const height = try canopy_photosynthesis.allocateLeafAcrossCanopyLayers(2, 4, 0.4, 0.08, 100, 2, 0.2, 0.05, 2, &.{ 0, 0.25, 0.5, 2 }, &.{ 0.2, 0.5, 0.8, 1.0 }, &.{ 0.1, 0.2, 0.3, 0.4 }, .{ .area_m2 = &area, .carbon_g = &carbon, .nitrogen_g = &nitrogen, .phosphorus_g = &phosphorus });
    try std.testing.expect(height > 0.25);
    var area_sum: f64 = 0;
    var carbon_sum: f64 = 0;
    var nitrogen_sum: f64 = 0;
    var phosphorus_sum: f64 = 0;
    for (area, carbon, nitrogen, phosphorus) |a, c, n, p| {
        area_sum += a;
        carbon_sum += c;
        nitrogen_sum += n;
        phosphorus_sum += p;
    }
    try std.testing.expectApproxEqAbs(@as(f64, 2), area_sum, 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 4), carbon_sum, 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), nitrogen_sum, 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.08), phosphorus_sum, 1.0e-14);
}

test "GROSUB stalk geometry and runtime layer allocation conserve surface area" {
    var layer_area = [_]f64{0} ** 3;
    const result = try canopy_photosynthesis.allocateStalkAcrossCanopyLayers(20, 3, 10, 1.0e-6, 0.1, 1.6, true, false, &.{ 0, 0.5, 1.0, 2.0 }, &layer_area);
    try std.testing.expect(result.radius_m > 0);
    try std.testing.expectEqual(@as(f64, 20), result.sapwood_carbon_g);
    var sum: f64 = 0;
    for (layer_area) |area| sum += area;
    try std.testing.expectApproxEqAbs(result.surface_area_m2, sum, 1.0e-14);
    var no_height_area = [_]f64{1} ** 3;
    const no_height = try canopy_photosynthesis.allocateStalkAcrossCanopyLayers(20, 3, 10, 1.0e-6, 0.1, 0.1, false, false, &.{ 0, 0.5, 1.0, 2.0 }, &no_height_area);
    try std.testing.expectEqual(@as(f64, 3), no_height.sapwood_carbon_g);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, &no_height_area);
}

test "GROSUB potential sites and post-anthesis seed number retain source limitations" {
    const sites = try canopy_photosynthesis.accumulatePotentialSeedSites(10, true, false, 20, 100, 5, 3, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 13), sites, 1.0e-15);
    const result = try canopy_photosynthesis.updateSeedNumberAndSize(.{ .anthesis_started = true, .grain_fill_started = true, .final_seed_number_set = false, .maximum_seed_size_set = false, .mobile_carbon_concentration_g_per_g = 0.2, .mobile_nitrogen_concentration_g_per_g = 0.1, .mobile_phosphorus_concentration_g_per_g = 0.05, .carbon_half_saturation_g_per_g = 0.2, .nitrogen_half_saturation_g_per_g = 0.1, .phosphorus_half_saturation_g_per_g = 0.05, .canopy_temperature_c = 42, .chilling_temperature_c = 5, .high_temperature_c = 40, .seed_loss_fraction_per_c_h = 0.01, .timestep_h = 1, .water_growth_fraction = 1, .reproductive_stage_increment = 0.1, .maximum_seeds_per_site = 4, .potential_site_count = 13, .current_seed_count = 10, .maximum_individual_seed_carbon_g = 0.5, .current_individual_seed_carbon_g = 0.1 });
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), result.nutrient_set_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.02), result.thermal_loss_fraction, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 12.4), result.seed_count, 1.0e-14);
    try std.testing.expect(result.individual_seed_carbon_g > 0.1);
}

test "GROSUB end of maximum seed-size setting closes the whole seed-set block" {
    const result = try canopy_photosynthesis.updateSeedNumberAndSize(.{
        .anthesis_started = true,
        .grain_fill_started = true,
        .final_seed_number_set = false,
        .maximum_seed_size_set = true,
        .mobile_carbon_concentration_g_per_g = 0,
        .mobile_nitrogen_concentration_g_per_g = 0,
        .mobile_phosphorus_concentration_g_per_g = 0,
        .carbon_half_saturation_g_per_g = 0,
        .nitrogen_half_saturation_g_per_g = 0,
        .phosphorus_half_saturation_g_per_g = 0,
        .canopy_temperature_c = 50,
        .chilling_temperature_c = 5,
        .high_temperature_c = 40,
        .seed_loss_fraction_per_c_h = 1,
        .timestep_h = 1,
        .water_growth_fraction = 1,
        .reproductive_stage_increment = 1,
        .maximum_seeds_per_site = 10,
        .potential_site_count = 10,
        .current_seed_count = 7,
        .maximum_individual_seed_carbon_g = 1,
        .current_individual_seed_carbon_g = 0.4,
    });
    try std.testing.expectEqual(@as(f64, 0), result.nutrient_set_fraction);
    try std.testing.expectEqual(@as(f64, 0), result.thermal_loss_fraction);
    try std.testing.expectEqual(@as(f64, 7), result.seed_count);
    try std.testing.expectEqual(@as(f64, 0.4), result.individual_seed_carbon_g);
}

test "GROSUB impossible negative seed count fails instead of silently clamping" {
    try std.testing.expectError(error.NegativeSeedCount, canopy_photosynthesis.updateSeedNumberAndSize(.{
        .anthesis_started = true,
        .grain_fill_started = false,
        .final_seed_number_set = false,
        .maximum_seed_size_set = false,
        .mobile_carbon_concentration_g_per_g = 0,
        .mobile_nitrogen_concentration_g_per_g = 0,
        .mobile_phosphorus_concentration_g_per_g = 0,
        .carbon_half_saturation_g_per_g = 1,
        .nitrogen_half_saturation_g_per_g = 1,
        .phosphorus_half_saturation_g_per_g = 1,
        .canopy_temperature_c = 100,
        .chilling_temperature_c = 5,
        .high_temperature_c = 40,
        .seed_loss_fraction_per_c_h = 1,
        .timestep_h = 1,
        .water_growth_fraction = 1,
        .reproductive_stage_increment = 0,
        .maximum_seeds_per_site = 1,
        .potential_site_count = 1,
        .current_seed_count = 1,
        .maximum_individual_seed_carbon_g = 1,
        .current_individual_seed_carbon_g = 0,
    }));
}

test "branch insertion atomically preserves compact topology state" {
    const branch_counts = [_]usize{ 1, 1 };
    const node_counts = [_]usize{ 1, 1 };
    const sample_counts = [_]usize{ 2, 1 };
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 2, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.branch_c3_feedback_fraction[0] = 0.4;
    state.branch_c3_feedback_fraction[1] = 0.8;
    state.node_co2_limited_carboxylation_umol_per_m2_s[1] = 12;
    state.sample_carboxylation_umol_per_s[2] = 7;
    const new_samples = [_]usize{ 3, 4 };
    try std.testing.expectEqual(@as(usize, 1), try state.appendBranch(0, &new_samples));
    try std.testing.expectEqual(canopy_photosynthesis.Range{ .first = 0, .end = 2 }, try state.branchRange(0));
    try std.testing.expectEqual(canopy_photosynthesis.Range{ .first = 2, .end = 3 }, try state.branchRange(1));
    try std.testing.expectEqual(@as(usize, 3), state.branch_c3_feedback_fraction.len);
    try std.testing.expectEqual(@as(usize, 4), state.node_co2_limited_carboxylation_umol_per_m2_s.len);
    try std.testing.expectEqual(@as(usize, 10), state.sample_carboxylation_umol_per_s.len);
    try std.testing.expectEqual(0.4, state.branch_c3_feedback_fraction[0]);
    try std.testing.expectEqual(0.8, state.branch_c3_feedback_fraction[2]);
    try std.testing.expectEqual(12, state.node_co2_limited_carboxylation_umol_per_m2_s[3]);
    try std.testing.expectEqual(7, state.sample_carboxylation_umol_per_s[9]);
}

test "node insertion preserves following runtime node and sample state" {
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 2, &.{ 1, 1 }, &.{ 1, 1 }, &.{ 2, 1 });
    defer state.deinit();
    state.branch_leaf_carbon_g[1] = 4;
    state.node_leaf_carbon_g[0] = 1;
    state.node_leaf_carbon_g[1] = 2;
    state.sample_carboxylation_umol_per_s[2] = 7;
    try std.testing.expectEqual(@as(usize, 1), try state.appendNode(0, 3));
    try std.testing.expectEqual(canopy_photosynthesis.Range{ .first = 0, .end = 2 }, try state.nodeRange(0));
    try std.testing.expectEqual(canopy_photosynthesis.Range{ .first = 2, .end = 3 }, try state.nodeRange(1));
    try std.testing.expectEqual(@as(usize, 6), state.sample_carboxylation_umol_per_s.len);
    try std.testing.expectEqual(@as(f64, 1), state.node_leaf_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.node_leaf_carbon_g[1]);
    try std.testing.expectEqual(@as(f64, 2), state.node_leaf_carbon_g[2]);
    try std.testing.expectEqual(@as(f64, 7), state.sample_carboxylation_umol_per_s[5]);
    try std.testing.expectEqual(@as(f64, 4), state.branch_leaf_carbon_g[1]);
}

test "GROSUB branch pools and C4 transfer conserve internal carbon" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{1};
    const sample_counts = [_]usize{0};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.branch_mobile_carbon_g[0] = 2;
    state.branch_mobile_nitrogen_g[0] = 0.2;
    state.branch_mobile_phosphorus_g[0] = 0.03;
    try canopy_photosynthesis.updateBranchMobilePools(&state, 0, .{ .fixed_carbon_g = 1, .maintenance_respiration_demand_g_c = 0.4, .available_respirable_carbon_g_c = 0.25, .growth_and_respiration_g_c = 0.5, .nitrogen_assimilation_respiration_g_c = 0.1, .assimilated_nitrogen_g = 0.02, .canopy_ammonia_exchange_g_n = 0.005, .assimilated_phosphorus_g = 0.01 });
    try std.testing.expectApproxEqAbs(2.15, state.branch_mobile_carbon_g[0], 1e-14);
    try std.testing.expectApproxEqAbs(0.185, state.branch_mobile_nitrogen_g[0], 1e-14);
    try std.testing.expectApproxEqAbs(0.02, state.branch_mobile_phosphorus_g[0], 1e-14);

    state.node_leaf_carbon_g[0] = 10;
    state.node_c3_nonstructural_carbon_g[0] = 0.5;
    state.node_c4_mesophyll_nonstructural_carbon_g[0] = 0.8;
    state.node_bundle_sheath_co2_carbon_g[0] = 1e-5;
    state.node_bundle_sheath_bicarbonate_carbon_g[0] = 2e-5;
    const before = state.node_c3_nonstructural_carbon_g[0] + state.node_c4_mesophyll_nonstructural_carbon_g[0] + state.node_bundle_sheath_co2_carbon_g[0] + state.node_bundle_sheath_bicarbonate_carbon_g[0];
    var parameters = canopy_photosynthesis.sourceC4CarbonParameters();
    parameters.decarboxylated_co2_fraction = 0.7;
    const expected_exchange = try c4_mesophyll_bundle_exchange.exchange(.{
        .bundle_sheath_nonstructural_carbon_g_c = 0.5,
        .mesophyll_nonstructural_carbon_g_c = 0.8,
        .bundle_sheath_fixation_g_c_per_timestep = 0.01,
        .mesophyll_fixation_g_c_per_timestep = 0.02,
        .leaf_carbon_g_c = 10,
        .bundle_sheath_water_g_h2o_per_g_c = parameters.bundle_sheath_water_g_per_g_c,
        .mesophyll_water_g_h2o_per_g_c = parameters.mesophyll_water_g_per_g_c,
        .timestep_h = 0.1,
    });
    const fluxes = try canopy_photosynthesis.advanceC4CarbonPools(&state, 0, 0.02, 0.01, 10, parameters, 0.1);
    try std.testing.expectEqual(expected_exchange.mesophyll_to_bundle_sheath_carbon_g_c, fluxes.mesophyll_to_bundle_sheath_g_c);
    const after = state.node_c3_nonstructural_carbon_g[0] + state.node_c4_mesophyll_nonstructural_carbon_g[0] + state.node_bundle_sheath_co2_carbon_g[0] + state.node_bundle_sheath_bicarbonate_carbon_g[0];
    try std.testing.expectApproxEqAbs(before + 0.02 - 0.01 - fluxes.bundle_sheath_co2_leakage_g_c, after, 1e-13);
}

test "GROSUB leaf growth uses all runtime concurrent nodes without a ring" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{30};
    const sample_counts = [_]usize{0} ** 30;
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    try canopy_photosynthesis.distributeLeafGrowth(&state, 0, 29, 0, 5, .{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.02 }, 2, 10, 1, 0.02, 0.01, 1, 0, 1);
    for (0..25) |node| try std.testing.expectEqual(0, state.node_leaf_carbon_g[node]);
    for (25..30) |node| try std.testing.expectApproxEqAbs(0.2, state.node_leaf_carbon_g[node], 1e-15);
    try std.testing.expectApproxEqAbs(0.02, state.branch_leaf_area_m2[0], 1e-15);
}

test "GROSUB organ partition and branch commit preserve all C N P products" {
    const partition = [canopy_photosynthesis.organ_count]f64{ 0.2, 0.1, 0.2, 0.15, 0.1, 0.1, 0.15 };
    const yields = [canopy_photosynthesis.organ_count]f64{ 0.8, 0.7, 0.75, 0.9, 0.7, 0.65, 0.85 };
    const n_to_c = [canopy_photosynthesis.organ_count]f64{ 0.04, 0.02, 0.01, 0.015, 0.012, 0.02, 0.015 };
    const p_to_c = [canopy_photosynthesis.organ_count]f64{ 0.004, 0.002, 0.001, 0.0015, 0.0012, 0.002, 0.0015 };
    const growth = try canopy_photosynthesis.calculateOrganGrowth(10, partition, yields, n_to_c, p_to_c, 0.25, 0.5, 0.8);
    try std.testing.expectApproxEqAbs(1.6, growth.value(.leaf).carbon_g, 1e-15);
    try std.testing.expectApproxEqAbs(1.6 * 0.04 * 0.625, growth.value(.leaf).nitrogen_g, 1e-15);
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{0};
    const sample_counts = [_]usize{};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    try canopy_photosynthesis.applyBranchOrganGrowth(&state, 0, growth);
    try std.testing.expectApproxEqAbs(growth.value(.leaf).carbon_g, state.branch_leaf_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(growth.value(.ear).phosphorus_g, state.branch_ear_phosphorus_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(growth.value(.reserve).carbon_g, state.branch_reserve_carbon_g[0], 1e-15);
    try std.testing.expectEqual(0, state.branch_grain_carbon_g[0]);
}

test "organ growth publication is atomic when a destination overflows" {
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &.{1}, &.{0}, &.{});
    defer state.deinit();
    state.branch_leaf_carbon_g[0] = 2;
    state.branch_sheath_carbon_g[0] = std.math.floatMax(f64);
    var growth: canopy_photosynthesis.OrganGrowth = .{ .carbon_g = @splat(0), .nitrogen_g = @splat(0), .phosphorus_g = @splat(0), .total_shoot_carbon_production_g = 2 };
    growth.carbon_g[@intFromEnum(canopy_photosynthesis.Organ.leaf)] = 1;
    growth.carbon_g[@intFromEnum(canopy_photosynthesis.Organ.sheath)] = std.math.floatMax(f64);
    try std.testing.expectError(error.InvalidBranchOrganGrowthTransaction, canopy_photosynthesis.applyBranchOrganGrowth(&state, 0, growth));
    try std.testing.expectEqual(@as(f64, 2), state.branch_leaf_carbon_g[0]);
    try std.testing.expectEqual(std.math.floatMax(f64), state.branch_sheath_carbon_g[0]);
}

test "GROSUB sheath and stalk growth retain runtime node geometry" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{5};
    const sample_counts = [_]usize{0} ** 5;
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    @memset(state.node_leaf_carbon_g, 1);
    try canopy_photosynthesis.distributeSheathGrowth(&state, 0, 4, 1, 3, .{ .carbon_g = 0.6, .nitrogen_g = 0.06, .phosphorus_g = 0.012 }, 2, 10, 1.5, 0.4, 0.01, 2, 0, 0.8, 0.5);
    try std.testing.expectEqual(0, state.node_sheath_carbon_g[1]);
    for (2..5) |node| try std.testing.expectApproxEqAbs(0.2, state.node_sheath_carbon_g[node], 1e-15);
    const stalk = try canopy_photosynthesis.distributeStalkGrowth(&state, 0, 2, 4, .{ .carbon_g = 0.9, .nitrogen_g = 0.09, .phosphorus_g = 0.009 }, 1, 0.5, 0.01, 2, 0, 1, 0.8, 1e-5);
    try std.testing.expect(stalk.stem_diameter_m > 0);
    try std.testing.expect(state.node_height_m[4] > state.node_height_m[3]);
    try std.testing.expectApproxEqAbs(0.9, state.node_internode_carbon_g[2] + state.node_internode_carbon_g[3] + state.node_internode_carbon_g[4], 1e-15);
}

test "GROSUB stalk allocation accepts an exact runtime window wider than twenty five nodes" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{32};
    const sample_counts = [_]usize{0} ** 32;
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();

    _ = try canopy_photosynthesis.distributeStalkGrowth(
        &state,
        0,
        1,
        30,
        .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 },
        1,
        0.5,
        0.01,
        2,
        0,
        1,
        0.8,
        1e-5,
    );
    try std.testing.expectEqual(@as(f64, 0), state.node_internode_carbon_g[0]);
    for (1..31) |node| try std.testing.expectApproxEqAbs(@as(f64, 0.1), state.node_internode_carbon_g[node], 1e-15);
    try std.testing.expectEqual(@as(f64, 0), state.node_internode_carbon_g[31]);
}

test "GROSUB recycling and leaf sheath senescence conserve C N P" {
    const recycling = try canopy_photosynthesis.recyclingFractions(true, 0.2, 0.02, 0.002, 0.025, 0.0025, 0.1, 0.5, 0.8, 0.7);
    try std.testing.expect(recycling.carbon >= 0.1 and recycling.carbon <= 0.6);
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{1};
    const sample_counts = [_]usize{0};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.node_leaf_carbon_g[0] = 4;
    state.node_leaf_nitrogen_g[0] = 0.4;
    state.node_leaf_phosphorus_g[0] = 0.04;
    state.node_leaf_protein_g[0] = 0.3;
    state.node_leaf_area_m2[0] = 2;
    state.node_sheath_carbon_g[0] = 2;
    state.node_sheath_nitrogen_g[0] = 0.1;
    state.node_sheath_phosphorus_g[0] = 0.02;
    state.node_sheath_protein_g[0] = 0.08;
    state.node_sheath_height_m[0] = 0.5;
    state.node_internode_carbon_g[0] = 0.7;
    state.node_internode_nitrogen_g[0] = 0.01;
    state.node_internode_phosphorus_g[0] = 0.001;
    state.branch_leaf_carbon_g[0] = 4;
    state.branch_leaf_nitrogen_g[0] = 0.4;
    state.branch_leaf_phosphorus_g[0] = 0.04;
    state.branch_sheath_carbon_g[0] = 2;
    state.branch_sheath_nitrogen_g[0] = 0.1;
    state.branch_sheath_phosphorus_g[0] = 0.02;
    state.branch_leaf_area_m2[0] = 2;
    const kinetics: canopy_photosynthesis.KineticFractions = .{ .carbon = @splat(0.25), .nitrogen = @splat(0.25), .phosphorus = @splat(0.25) };
    const products = try canopy_photosynthesis.senesceLeafAndSheathNode(&state, 0, 0, 0.25, recycling, 0.5, 5, .{ 0.2, 0.8 }, .{ 0.1, 0.9 }, .{ 0.3, 0.7 }, .{ 0.15, 0.85 }, .{ 0.25, 0.75 }, kinetics, kinetics, kinetics);
    var litter_c: f64 = 0;
    var litter_n: f64 = 0;
    var litter_p: f64 = 0;
    for (0..4) |index| {
        litter_c += products.woody_carbon_g[index] + products.nonwoody_carbon_g[index];
        litter_n += products.woody_nitrogen_g[index] + products.nonwoody_nitrogen_g[index];
        litter_p += products.woody_phosphorus_g[index] + products.nonwoody_phosphorus_g[index];
    }
    try std.testing.expectApproxEqAbs(6.0, state.node_leaf_carbon_g[0] + state.node_sheath_carbon_g[0] + litter_c + products.recycled_carbon_g, 1e-13);
    try std.testing.expectApproxEqAbs(0.5, state.node_leaf_nitrogen_g[0] + state.node_sheath_nitrogen_g[0] + litter_n + products.recycled_nitrogen_g, 1e-13);
    try std.testing.expectApproxEqAbs(0.06, state.node_leaf_phosphorus_g[0] + state.node_sheath_phosphorus_g[0] + litter_p + products.recycled_phosphorus_g, 1e-13);
    try std.testing.expectEqual(0, state.node_internode_carbon_g[0]);
    try std.testing.expectEqual(0.7, state.branch_senescing_stalk_carbon_g[0]);
}

test "GROSUB node senescence splits demand and phenological carbon recovery" {
    const allocation = try canopy_photosynthesis.allocateNodeSenescenceDemand(1, 3, 1, 0.5, 0.25, 0.8);
    try std.testing.expectApproxEqAbs(0.5, allocation.leaf_fraction, 1e-15);
    try std.testing.expectApproxEqAbs(0.5, allocation.sheath_fraction, 1e-15);
    try std.testing.expectApproxEqAbs(0.2, allocation.carbon_recovered_to_mobile_pool_g, 1e-15);
    try std.testing.expectApproxEqAbs(0.6, allocation.carbon_respired_g, 1e-15);
    try std.testing.expectApproxEqAbs(0.2, allocation.remaining_respiration_demand_g_c, 1e-15);
}

test "GROSUB node senescence commit conserves structural mobile litter and respired elements" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{1};
    const sample_counts = [_]usize{0};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.node_leaf_carbon_g[0] = 3;
    state.node_leaf_nitrogen_g[0] = 0.3;
    state.node_leaf_phosphorus_g[0] = 0.03;
    state.node_leaf_area_m2[0] = 1.5;
    state.node_sheath_carbon_g[0] = 1;
    state.node_sheath_nitrogen_g[0] = 0.1;
    state.node_sheath_phosphorus_g[0] = 0.01;
    state.node_c3_nonstructural_carbon_g[0] = 0.1;
    state.node_c4_mesophyll_nonstructural_carbon_g[0] = 0.2;
    state.branch_leaf_carbon_g[0] = 3;
    state.branch_leaf_nitrogen_g[0] = 0.3;
    state.branch_leaf_phosphorus_g[0] = 0.03;
    state.branch_sheath_carbon_g[0] = 1;
    state.branch_sheath_nitrogen_g[0] = 0.1;
    state.branch_sheath_phosphorus_g[0] = 0.01;
    state.branch_leaf_area_m2[0] = 1.5;
    const allocation = try canopy_photosynthesis.allocateNodeSenescenceDemand(1, 3, 1, 0.5, 0.25, 0.8);
    const kinetics: canopy_photosynthesis.KineticFractions = .{ .carbon = @splat(0.25), .nitrogen = @splat(0.25), .phosphorus = @splat(0.25) };
    const products = try canopy_photosynthesis.commitNodeSenescenceDemand(&state, 0, 0, allocation, .{ .carbon = 0.5, .nitrogen = 0.6, .phosphorus = 0.7 }, 0.5, 5, .{ 0.2, 0.8 }, .{ 0.1, 0.9 }, .{ 0.3, 0.7 }, .{ 0.15, 0.85 }, .{ 0.25, 0.75 }, kinetics, kinetics, kinetics);
    var litter_c: f64 = 0;
    var litter_n: f64 = 0;
    var litter_p: f64 = 0;
    for (0..4) |index| {
        litter_c += products.woody_carbon_g[index] + products.nonwoody_carbon_g[index];
        litter_n += products.woody_nitrogen_g[index] + products.nonwoody_nitrogen_g[index];
        litter_p += products.woody_phosphorus_g[index] + products.nonwoody_phosphorus_g[index];
    }
    const remaining_c = state.node_leaf_carbon_g[0] + state.node_sheath_carbon_g[0] + state.node_c3_nonstructural_carbon_g[0] + state.node_c4_mesophyll_nonstructural_carbon_g[0];
    try std.testing.expectApproxEqAbs(4.3, remaining_c + litter_c + products.recycled_carbon_g + products.respired_carbon_g, 1e-13);
    try std.testing.expectApproxEqAbs(0.4, state.node_leaf_nitrogen_g[0] + state.node_sheath_nitrogen_g[0] + litter_n + products.recycled_nitrogen_g, 1e-13);
    try std.testing.expectApproxEqAbs(0.04, state.node_leaf_phosphorus_g[0] + state.node_sheath_phosphorus_g[0] + litter_p + products.recycled_phosphorus_g, 1e-13);
}

test "GROSUB internode senescence conserves stalk reserve litter and respired C N P" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{1};
    const sample_counts = [_]usize{0};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.branch_stalk_carbon_g[0] = 10;
    state.branch_stalk_nitrogen_g[0] = 1;
    state.branch_stalk_phosphorus_g[0] = 0.1;
    state.branch_sapwood_carbon_g[0] = 5;
    state.node_internode_carbon_g[0] = 4;
    state.node_internode_nitrogen_g[0] = 0.4;
    state.node_internode_phosphorus_g[0] = 0.04;
    state.node_internode_length_m[0] = 1;
    state.node_height_m[0] = 1;
    const kinetics: canopy_photosynthesis.KineticFractions = .{ .carbon = @splat(0.25), .nitrogen = @splat(0.25), .phosphorus = @splat(0.25) };
    const result = try canopy_photosynthesis.commitInternodeSenescenceDemand(&state, 0, 0, 0.6, 0.25, .{ .carbon = 0.6, .nitrogen = 0.7, .phosphorus = 0.8 }, .{ 0.25, 0.75 }, .{ 0.25, 0.75 }, .{ 0.25, 0.75 }, kinetics, kinetics);
    var litter_c: f64 = 0;
    var litter_n: f64 = 0;
    var litter_p: f64 = 0;
    for (0..4) |kinetic| {
        litter_c += result.products.woody_carbon_g[kinetic] + result.products.nonwoody_carbon_g[kinetic];
        litter_n += result.products.woody_nitrogen_g[kinetic] + result.products.nonwoody_nitrogen_g[kinetic];
        litter_p += result.products.woody_phosphorus_g[kinetic] + result.products.nonwoody_phosphorus_g[kinetic];
    }
    try std.testing.expectApproxEqAbs(0.5, result.fraction, 1e-15);
    try std.testing.expectApproxEqAbs(4.0, state.node_internode_carbon_g[0] + litter_c + result.products.recycled_carbon_g + result.products.respired_carbon_g, 1e-13);
    try std.testing.expectApproxEqAbs(0.4, state.node_internode_nitrogen_g[0] + litter_n + result.products.recycled_nitrogen_g, 1e-13);
    try std.testing.expectApproxEqAbs(0.04, state.node_internode_phosphorus_g[0] + litter_p + result.products.recycled_phosphorus_g, 1e-13);
    try std.testing.expectApproxEqAbs(8.0, state.branch_stalk_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(0.15, result.remaining_respiration_demand_g_c, 1e-15);
}

test "GROSUB residual stalk senescence conserves elements and lowers canopy height" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{2};
    const sample_counts = [_]usize{ 0, 0 };
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.branch_stalk_carbon_g[0] = 10;
    state.branch_stalk_nitrogen_g[0] = 1;
    state.branch_stalk_phosphorus_g[0] = 0.1;
    state.branch_sapwood_carbon_g[0] = 5;
    state.branch_senescing_stalk_carbon_g[0] = 4;
    state.branch_senescing_stalk_nitrogen_g[0] = 0.4;
    state.branch_senescing_stalk_phosphorus_g[0] = 0.04;
    state.node_height_m[0] = 1;
    state.node_height_m[1] = 2;
    const kinetics: canopy_photosynthesis.KineticFractions = .{ .carbon = @splat(0.25), .nitrogen = @splat(0.25), .phosphorus = @splat(0.25) };
    const result = try canopy_photosynthesis.commitResidualStalkSenescenceDemand(&state, 0, 0.6, 0.25, .{ .carbon = 0.6, .nitrogen = 0.7, .phosphorus = 0.8 }, .{ 0.25, 0.75 }, .{ 0.25, 0.75 }, .{ 0.25, 0.75 }, kinetics, kinetics);
    var litter_c: f64 = 0;
    var litter_n: f64 = 0;
    var litter_p: f64 = 0;
    for (0..4) |kinetic| {
        litter_c += result.products.woody_carbon_g[kinetic] + result.products.nonwoody_carbon_g[kinetic];
        litter_n += result.products.woody_nitrogen_g[kinetic] + result.products.nonwoody_nitrogen_g[kinetic];
        litter_p += result.products.woody_phosphorus_g[kinetic] + result.products.nonwoody_phosphorus_g[kinetic];
    }
    try std.testing.expectApproxEqAbs(4.0, state.branch_senescing_stalk_carbon_g[0] + litter_c + result.products.recycled_carbon_g + result.products.respired_carbon_g, 1e-13);
    try std.testing.expectApproxEqAbs(0.4, state.branch_senescing_stalk_nitrogen_g[0] + litter_n + result.products.recycled_nitrogen_g, 1e-13);
    try std.testing.expectApproxEqAbs(0.04, state.branch_senescing_stalk_phosphorus_g[0] + litter_p + result.products.recycled_phosphorus_g, 1e-13);
    try std.testing.expectApproxEqAbs(1.0, state.node_height_m[1], 1e-15);
    try std.testing.expectApproxEqAbs(8.0, state.branch_stalk_carbon_g[0], 1e-15);
}

test "GROSUB leaf nutrient equilibration preserves N P and source coupling" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{1};
    const sample_counts = [_]usize{0};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.node_leaf_carbon_g[0] = 10;
    state.node_leaf_nitrogen_g[0] = 1;
    state.node_leaf_phosphorus_g[0] = 0.1;
    state.node_leaf_protein_g[0] = 2;
    state.branch_leaf_nitrogen_g[0] = 1;
    state.branch_leaf_phosphorus_g[0] = 0.1;
    state.branch_mobile_carbon_g[0] = 10;
    const initial_n = state.node_leaf_nitrogen_g[0] + state.branch_mobile_nitrogen_g[0];
    const initial_p = state.node_leaf_phosphorus_g[0] + state.branch_mobile_phosphorus_g[0];
    const flux = try canopy_photosynthesis.remobilizeNodeLeafNutrients(&state, 0, 0, 1.0e-3, 0.1, 0.08, 0.008, 2, 20);
    try std.testing.expectApproxEqAbs(0.0005, flux.nitrogen_g, 1e-15);
    try std.testing.expectApproxEqAbs(0.00005, flux.phosphorus_g, 1e-15);
    try std.testing.expectApproxEqAbs(initial_n, state.node_leaf_nitrogen_g[0] + state.branch_mobile_nitrogen_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(initial_p, state.node_leaf_phosphorus_g[0] + state.branch_mobile_phosphorus_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(1.999, state.node_leaf_protein_g[0], 1e-15);
}

test "GROSUB branch senescence cascade exits after reserve satisfies node remainder" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{1};
    const sample_counts = [_]usize{0};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.node_leaf_carbon_g[0] = 2;
    state.node_leaf_nitrogen_g[0] = 0.2;
    state.node_leaf_phosphorus_g[0] = 0.02;
    state.node_leaf_area_m2[0] = 1;
    state.branch_leaf_carbon_g[0] = 2;
    state.branch_leaf_nitrogen_g[0] = 0.2;
    state.branch_leaf_phosphorus_g[0] = 0.02;
    state.branch_leaf_area_m2[0] = 1;
    state.branch_reserve_carbon_g[0] = 0.1;
    const kinetics: canopy_photosynthesis.KineticFractions = .{ .carbon = @splat(0.25), .nitrogen = @splat(0.25), .phosphorus = @splat(0.25) };
    const litter: canopy_photosynthesis.SenescenceLitterParameters = .{
        .woody_carbon_fraction = .{ 0.25, 0.75 },
        .leaf_woody_nitrogen_fraction = .{ 0.25, 0.75 },
        .sheath_woody_nitrogen_fraction = .{ 0.25, 0.75 },
        .stalk_woody_nitrogen_fraction = .{ 0.25, 0.75 },
        .leaf_woody_phosphorus_fraction = .{ 0.25, 0.75 },
        .sheath_woody_phosphorus_fraction = .{ 0.25, 0.75 },
        .stalk_woody_phosphorus_fraction = .{ 0.25, 0.75 },
        .woody_kinetics = kinetics,
        .leaf_kinetics = kinetics,
        .sheath_kinetics = kinetics,
        .stalk_kinetics = kinetics,
    };
    const result = try canopy_photosynthesis.commitBranchSenescenceDemand(&state, 0, .{ .total_respiration_demand_g_c = 0.3, .phenological_senescence_fraction = 0.5, .first_node_within_branch = 0, .last_node_within_branch = 0, .node_group_count = 1, .perennial = false, .reserve_fallback_policy = .consume_available, .demand_tolerance_g_c = 1e-12 }, .{ .carbon = 0.5, .nitrogen = 0.6, .phosphorus = 0.7 }, 2, 20, litter);
    var litter_c: f64 = 0;
    for (0..4) |kinetic| litter_c += result.products.woody_carbon_g[kinetic] + result.products.nonwoody_carbon_g[kinetic];
    try std.testing.expectEqual(0, result.remaining_respiration_demand_g_c);
    try std.testing.expectApproxEqAbs(0.075, result.reserve_carbon_respired_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(2.1, state.node_leaf_carbon_g[0] + state.branch_reserve_carbon_g[0] + litter_c + result.products.recycled_carbon_g + result.products.respired_carbon_g + result.reserve_carbon_respired_g_c, 1e-13);
}

test "GROSUB cascade does not repeat final node and reserve equality is strict" {
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{0});
    defer state.deinit();
    state.node_leaf_carbon_g[0] = 10;
    state.branch_leaf_carbon_g[0] = 10;
    state.branch_reserve_carbon_g[0] = 2.25;
    const kinetics: canopy_photosynthesis.KineticFractions = .{ .carbon = @splat(0.25), .nitrogen = @splat(0.25), .phosphorus = @splat(0.25) };
    const litter: canopy_photosynthesis.SenescenceLitterParameters = .{
        .woody_carbon_fraction = .{ 0.25, 0.75 },
        .leaf_woody_nitrogen_fraction = .{ 0.25, 0.75 },
        .sheath_woody_nitrogen_fraction = .{ 0.25, 0.75 },
        .stalk_woody_nitrogen_fraction = .{ 0.25, 0.75 },
        .leaf_woody_phosphorus_fraction = .{ 0.25, 0.75 },
        .sheath_woody_phosphorus_fraction = .{ 0.25, 0.75 },
        .stalk_woody_phosphorus_fraction = .{ 0.25, 0.75 },
        .woody_kinetics = kinetics,
        .leaf_kinetics = kinetics,
        .sheath_kinetics = kinetics,
        .stalk_kinetics = kinetics,
    };
    const result = try canopy_photosynthesis.commitBranchSenescenceDemand(&state, 0, .{
        .total_respiration_demand_g_c = 3,
        .phenological_senescence_fraction = 0,
        .first_node_within_branch = 0,
        .last_node_within_branch = 0,
        .node_group_count = 3,
        .perennial = false,
        .demand_tolerance_g_c = 0,
    }, .{ .carbon = 1, .nitrogen = 1, .phosphorus = 1 }, 0, 0, litter);
    try std.testing.expectEqual(@as(f64, 9), state.node_leaf_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 1), result.remaining_respiration_demand_g_c);
    try std.testing.expectEqual(@as(f64, 1), state.branch_reserve_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 1.25), result.reserve_carbon_respired_g_c);
}

test "GROSUB reserve to grain fill conserves precursor and translocated C N P" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{0};
    const sample_counts = [_]usize{};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.branch_reserve_carbon_g[0] = 2;
    state.branch_reserve_nitrogen_g[0] = 0.2;
    state.branch_reserve_phosphorus_g[0] = 0.02;
    state.branch_grain_carbon_g[0] = 1;
    state.branch_grain_nitrogen_g[0] = 0.05;
    state.branch_grain_phosphorus_g[0] = 0.005;
    const precursor: canopy_photosynthesis.LeafGrowth = .{ .carbon_g = 0.15, .nitrogen_g = 0.003, .phosphorus_g = 0.0003 };
    const result = try canopy_photosynthesis.fillGrainFromReserve(&state, 0, true, 100, 0.1, 0.01, 1, 1, 0.5, 0.04, 0.004, 0.02, 0.002, precursor);
    try std.testing.expectApproxEqAbs(1.0, result.carbon_translocated_g, 1e-15);
    try std.testing.expectApproxEqAbs(3.0 + precursor.carbon_g, state.branch_reserve_carbon_g[0] + state.branch_grain_carbon_g[0], 1e-14);
    try std.testing.expectApproxEqAbs(0.25 + precursor.nitrogen_g, state.branch_reserve_nitrogen_g[0] + state.branch_grain_nitrogen_g[0], 1e-14);
    try std.testing.expectApproxEqAbs(0.025 + precursor.phosphorus_g, state.branch_reserve_phosphorus_g[0] + state.branch_grain_phosphorus_g[0], 1e-14);
}

test "GROSUB grain nutrient fill uses three-way minimum rather than reserve-deficit maximum" {
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &.{1}, &.{0}, &.{});
    defer state.deinit();
    state.branch_reserve_carbon_g[0] = 10;
    state.branch_reserve_nitrogen_g[0] = 0.1;
    state.branch_reserve_phosphorus_g[0] = 0.1;
    state.branch_grain_carbon_g[0] = 1;
    state.branch_grain_nitrogen_g[0] = 0;
    state.branch_grain_phosphorus_g[0] = 0;

    const result = try canopy_photosynthesis.fillGrainFromReserve(
        &state,
        0,
        true,
        1,
        10,
        1,
        1,
        1,
        0,
        1,
        1,
        0,
        0,
        .{ .carbon_g = 0, .nitrogen_g = 0, .phosphorus_g = 0 },
    );
    // The grain deficit is 2 g for both nutrients, but each reserve contains
    // only 0.1 g. The former production max selected the deficit and overdrawn
    // the reserve; the source minimum selects the available reserve term.
    try std.testing.expectApproxEqAbs(@as(f64, 1), result.carbon_translocated_g, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), result.nitrogen_translocated_g, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), result.phosphorus_translocated_g, 1e-15);
}

test "GROSUB grain fill rejects non-finite authoritative pool state" {
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &.{1}, &.{0}, &.{});
    defer state.deinit();
    state.branch_grain_carbon_g[0] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidGrainFillState, canopy_photosynthesis.fillGrainFromReserve(
        &state,
        0,
        true,
        1,
        1,
        1,
        1,
        1,
        0.5,
        0.04,
        0.004,
        0.02,
        0.002,
        .{ .carbon_g = 0, .nitrogen_g = 0, .phosphorus_g = 0 },
    ));
}

test "GROSUB reserve respiration demand and phenological senescence are exact" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{0};
    const sample_counts = [_]usize{};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.branch_reserve_carbon_g[0] = 4;
    const remaining = try canopy_photosynthesis.consumeReserveForRespiration(&state, 0, false, 1, 0.1, 0.5, 1);
    try std.testing.expectApproxEqAbs(0.8, remaining, 1e-15);
    try std.testing.expectApproxEqAbs(3.8, state.branch_reserve_carbon_g[0], 1e-15);
    const remobilizing_remaining = try canopy_photosynthesis.consumeReserveForRespiration(&state, 0, true, 1, 0.1, 0.5, 1);
    try std.testing.expectEqual(@as(f64, 1), remobilizing_remaining);
    try std.testing.expectApproxEqAbs(3.8, state.branch_reserve_carbon_g[0], 1e-15);
    const demand = try canopy_photosynthesis.senescenceDemand(true, true, true, 6, 1, 0.01, 20, 50, 100, 1, 0.2, 10, 2);
    try std.testing.expectApproxEqAbs(0.3, demand.phenological_respiration_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(0.5, demand.total_respiration_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(0.6, demand.phenological_fraction, 1e-15);
    try std.testing.expectEqual(@as(usize, 5), demand.node_group_count);
}

test "GROSUB interbranch reserve exchange is pairwise conservative" {
    const branch_counts = [_]usize{2};
    const node_counts = [_]usize{ 0, 0 };
    const sample_counts = [_]usize{};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.branch_sapwood_carbon_g[0] = 2;
    state.branch_sapwood_carbon_g[1] = 1;
    state.branch_reserve_carbon_g[0] = 3;
    state.branch_reserve_nitrogen_g[0] = 0.3;
    state.branch_reserve_phosphorus_g[0] = 0.03;
    const flux = try canopy_photosynthesis.equilibrateBranchReserves(&state, 0, 1, 0.5, 0.5, 1, 0);
    try std.testing.expectApproxEqAbs(0.5, flux.carbon_g, 1e-15);
    try std.testing.expectApproxEqAbs(3.0, state.branch_reserve_carbon_g[0] + state.branch_reserve_carbon_g[1], 1e-15);
    try std.testing.expectApproxEqAbs(0.3, state.branch_reserve_nitrogen_g[0] + state.branch_reserve_nitrogen_g[1], 1e-15);
    try std.testing.expectApproxEqAbs(0.03, state.branch_reserve_phosphorus_g[0] + state.branch_reserve_phosphorus_g[1], 1e-15);
}

test "GROSUB harvest retention partitions export litter and remaining pools" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{2};
    const sample_counts = [_]usize{ 0, 0 };
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.branch_stalk_carbon_g[0] = 10;
    state.branch_stalk_nitrogen_g[0] = 1;
    state.branch_stalk_phosphorus_g[0] = 0.1;
    state.branch_sapwood_carbon_g[0] = 4;
    state.branch_senescing_stalk_carbon_g[0] = 2;
    state.branch_reserve_carbon_g[0] = 5;
    state.branch_reserve_nitrogen_g[0] = 0.5;
    state.branch_reserve_phosphorus_g[0] = 0.05;
    const products = try canopy_photosynthesis.harvestBranchStalkAndReserve(&state, 0, 0.4, 0.7, 0.2, 0.8);
    try std.testing.expectApproxEqAbs(15.0, state.branch_stalk_carbon_g[0] + state.branch_reserve_carbon_g[0] + products.ecosystem_export.carbon_g + products.litter.carbon_g, 1e-14);
    try std.testing.expectApproxEqAbs(1.5, state.branch_stalk_nitrogen_g[0] + state.branch_reserve_nitrogen_g[0] + products.ecosystem_export.nitrogen_g + products.litter.nitrogen_g, 1e-14);
    try std.testing.expectApproxEqAbs(1.6, state.branch_sapwood_carbon_g[0], 1e-15);
    state.branch_mobile_carbon_g[0] = 2;
    state.branch_mobile_nitrogen_g[0] = 0.2;
    state.branch_mobile_phosphorus_g[0] = 0.02;
    state.node_c3_nonstructural_carbon_g[0] = 0.4;
    state.node_c4_mesophyll_nonstructural_carbon_g[1] = 0.6;
    const removed = try canopy_photosynthesis.harvestBranchMobilePools(&state, 0, 0.25);
    try std.testing.expectApproxEqAbs(2.25, removed.carbon_g, 1e-15);
    try std.testing.expectApproxEqAbs(0.5, state.branch_mobile_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(0.1, state.node_c3_nonstructural_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(0.15, state.node_c4_mesophyll_nonstructural_carbon_g[1], 1e-15);
    try std.testing.expectApproxEqAbs(3.0, removed.carbon_g + state.branch_mobile_carbon_g[0] + state.node_c3_nonstructural_carbon_g[0] + state.node_c4_mesophyll_nonstructural_carbon_g[1], 1e-15);
}

test "GROSUB branch mobile removal preserves source ratios and thresholds" {
    try std.testing.expectEqual(
        @as(f64, 0.25),
        try canopy_photosynthesis.sourceOrderNonGrazingMobileRetention(8, 2, 1.0e-12),
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        try canopy_photosynthesis.sourceOrderNonGrazingMobileRetention(1.0e-12, 0.5e-12, 1.0e-12),
    );
    const result = try canopy_photosynthesis.sourceOrderProportionalMobileRemoval(
        .{ .carbon_g = 4, .nitrogen_g = 2, .phosphorus_g = 1 },
        1,
        1.0e-12,
    );
    try std.testing.expectEqual(@as(f64, 3), result.remaining.carbon_g);
    try std.testing.expectEqual(@as(f64, 1.5), result.remaining.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 0.75), result.remaining.phosphorus_g);
    try std.testing.expectEqual(@as(f64, 0.75), result.unclamped_carbon_retention_fraction);
}

test "GROSUB symbiont structural scaling exposes source overdraw" {
    const result = try canopy_photosynthesis.sourceOrderProportionalMobileRemoval(
        .{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 },
        2,
        1.0e-12,
    );
    try std.testing.expectEqual(@as(f64, 0), result.remaining.carbon_g);
    try std.testing.expectEqual(@as(f64, -1), result.unclamped_carbon_retention_fraction);
}

test "GROSUB C4 intermediate retention uses host mobile carbon ratio" {
    try std.testing.expectEqual(@as(f64, 0.25), try canopy_photosynthesis.sourceOrderC4IntermediateRetention(true, 4, 1, 1.0e-12));
    try std.testing.expectEqual(@as(f64, 1), try canopy_photosynthesis.sourceOrderC4IntermediateRetention(false, 4, 1, 1.0e-12));
    try std.testing.expectEqual(@as(f64, 1), try canopy_photosynthesis.sourceOrderC4IntermediateRetention(true, 1.0e-12, 0, 1.0e-12));
}

test "GROSUB layer leaf harvest reconciles sample node branch and products" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{1};
    const sample_counts = [_]usize{1};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.sample_leaf_area_m2[0] = 2;
    state.sample_exposed_leaf_area_m2[0] = 1;
    state.sample_stalk_area_m2[0] = 2;
    state.sample_leaf_carbon_g[0] = 4;
    state.sample_leaf_nitrogen_g[0] = 0.4;
    state.sample_leaf_phosphorus_g[0] = 0.04;
    state.node_leaf_area_m2[0] = 2;
    state.node_leaf_carbon_g[0] = 4;
    state.node_leaf_nitrogen_g[0] = 0.4;
    state.node_leaf_phosphorus_g[0] = 0.04;
    state.node_leaf_protein_g[0] = 1;
    state.branch_leaf_area_m2[0] = 2;
    state.branch_leaf_carbon_g[0] = 4;
    state.branch_leaf_nitrogen_g[0] = 0.4;
    state.branch_leaf_phosphorus_g[0] = 0.04;
    const products = try canopy_photosynthesis.harvestLeafLayerSample(&state, 0, 0, 0, .{ .remaining_fraction = 0.5, .unexported_fraction = 0.8, .height_below_cut_fraction = 0.5 }, .{ 0.25, 0.75 }, .{ 0.2, 0.8 }, .{ 0.1, 0.9 }, true);
    const exported_c = products.foliar.ecosystem_export.carbon_g + products.woody.ecosystem_export.carbon_g;
    const litter_c = products.foliar.litter.carbon_g + products.woody.litter.carbon_g;
    const exported_n = products.foliar.ecosystem_export.nitrogen_g + products.woody.ecosystem_export.nitrogen_g;
    const litter_n = products.foliar.litter.nitrogen_g + products.woody.litter.nitrogen_g;
    const exported_p = products.foliar.ecosystem_export.phosphorus_g + products.woody.ecosystem_export.phosphorus_g;
    const litter_p = products.foliar.litter.phosphorus_g + products.woody.litter.phosphorus_g;
    try std.testing.expectApproxEqAbs(4.0, state.sample_leaf_carbon_g[0] + exported_c + litter_c, 1e-15);
    try std.testing.expectApproxEqAbs(0.4, state.sample_leaf_nitrogen_g[0] + exported_n + litter_n, 1e-15);
    try std.testing.expectApproxEqAbs(0.04, state.sample_leaf_phosphorus_g[0] + exported_p + litter_p, 1e-15);
    try std.testing.expectEqual(state.sample_leaf_carbon_g[0], state.node_leaf_carbon_g[0]);
    try std.testing.expectEqual(state.node_leaf_carbon_g[0], state.branch_leaf_carbon_g[0]);
    try std.testing.expectApproxEqAbs(0.5, state.node_leaf_protein_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(1.0, state.sample_stalk_area_m2[0], 1e-15);
}

test "GROSUB node sheath harvest conserves elements and truncates at cutting plane" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{1};
    const sample_counts = [_]usize{0};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.node_height_m[0] = 1;
    state.node_sheath_height_m[0] = 2;
    state.node_sheath_carbon_g[0] = 2;
    state.node_sheath_nitrogen_g[0] = 0.2;
    state.node_sheath_phosphorus_g[0] = 0.02;
    state.node_sheath_protein_g[0] = 0.5;
    state.branch_sheath_carbon_g[0] = 2;
    state.branch_sheath_nitrogen_g[0] = 0.2;
    state.branch_sheath_phosphorus_g[0] = 0.02;
    const products = try canopy_photosynthesis.harvestNodeSheath(&state, 0, 0, 0.5, 0.8, .{ 0.25, 0.75 }, .{ 0.2, 0.8 }, .{ 0.1, 0.9 }, true, 2);
    const exported_c = products.nonwoody.ecosystem_export.carbon_g + products.woody.ecosystem_export.carbon_g;
    const litter_c = products.nonwoody.litter.carbon_g + products.woody.litter.carbon_g;
    try std.testing.expectApproxEqAbs(2.0, state.node_sheath_carbon_g[0] + exported_c + litter_c, 1e-15);
    try std.testing.expectEqual(state.node_sheath_carbon_g[0], state.branch_sheath_carbon_g[0]);
    try std.testing.expectApproxEqAbs(1.0, state.node_sheath_height_m[0], 1e-15);
    try std.testing.expectApproxEqAbs(0.25, state.node_sheath_protein_g[0], 1e-15);
}

test "GROSUB remaining node leaf rebuilds runtime layers and scales protein" {
    const result = try canopy_photosynthesis.sourceOrderRemainingNodeLeaf(
        &.{ 1, 2, 3 },
        &.{ 4, 5, 6 },
        &.{ 0.4, 0.5, 0.6 },
        &.{ 0.04, 0.05, 0.06 },
        12,
        8,
        1.0e-12,
    );
    try std.testing.expectEqual(@as(f64, 6), result.area_m2);
    try std.testing.expectEqual(@as(f64, 15), result.carbon_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), result.nitrogen_g_n, 1.0e-15);
    try std.testing.expectEqual(@as(f64, 4), result.protein_mass_g);
}

test "GROSUB node organ retention preserves leaf coupling and kind zero litter" {
    const coupled = try canopy_photosynthesis.sourceOrderNodeOrganRetention(false, false, 10, 5, 0.5, 0.25, 0, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 0.75), coupled.remaining_fraction);
    try std.testing.expectEqual(coupled.remaining_fraction, coupled.unexported_fraction);
    const kind_zero = try canopy_photosynthesis.sourceOrderNodeOrganRetention(false, true, 0, 0, 0.5, 0.25, 0.4, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 0.6), kind_zero.remaining_fraction);
    try std.testing.expectEqual(@as(f64, 0.9), kind_zero.unexported_fraction);
    const equality = try canopy_photosynthesis.sourceOrderNodeOrganRetention(false, false, 1.0e-12, 0, 0.5, 0.25, 0, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 0.75), equality.remaining_fraction);
}

test "GROSUB internode harvest retains source geometry and runtime node state" {
    try std.testing.expectApproxEqAbs(0.5, try canopy_photosynthesis.internodeHarvestRetention(3, 2, 2, false, 0, 1, false, 0, 10), 1e-15);
    try std.testing.expectApproxEqAbs(0.8, try canopy_photosynthesis.internodeHarvestRetention(3, 2, 2, false, 0.2, 1, false, 0, 10), 1e-15);
    try std.testing.expectApproxEqAbs(0.7, try canopy_photosynthesis.internodeHarvestRetention(3, 2, 2, false, 0, 1, true, 3, 10), 1e-15);
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{1};
    const sample_counts = [_]usize{0};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.node_height_m[0] = 3;
    state.node_internode_length_m[0] = 2;
    state.node_internode_carbon_g[0] = 4;
    state.node_internode_nitrogen_g[0] = 0.4;
    state.node_internode_phosphorus_g[0] = 0.04;
    try canopy_photosynthesis.commitInternodeHarvest(&state, 0, 0, 0.5, true, 2);
    try std.testing.expectApproxEqAbs(2.0, state.node_internode_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(1.0, state.node_internode_length_m[0], 1e-15);
    try std.testing.expectApproxEqAbs(2.0, state.node_height_m[0], 1e-15);
}

test "GROSUB branch stalk retention preserves kind zero and grazing operands" {
    const kind_zero = try canopy_photosynthesis.sourceOrderBranchStalkRetention(false, true, false, 4, 2, 0.5, 0.8, 10, 0, 0, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 0.5), kind_zero.remaining_fraction);
    try std.testing.expectEqual(@as(f64, 0.8), kind_zero.unexported_fraction);
    const grazing = try canopy_photosynthesis.sourceOrderBranchStalkRetention(true, false, false, 4, 2, 0, 0.8, 10, 2, 3, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 0.5), grazing.remaining_fraction);
}

test "GROSUB stalk reserve is discarded when no stalk remains" {
    const prior: canopy_photosynthesis.LayerHarvestRetention = .{ .remaining_fraction = 0.5, .unexported_fraction = 0.8, .height_below_cut_fraction = 0 };
    const retained = try canopy_photosynthesis.sourceOrderStalkReserveRetention(false, 2, prior, 4, 0, 1.0e-12);
    try std.testing.expectEqualDeep(prior, retained);
    const discarded = try canopy_photosynthesis.sourceOrderStalkReserveRetention(false, 1.0e-12, prior, 4, 0, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 0), discarded.remaining_fraction);
}

test "GROSUB cutting height layer retention and grazing demand retain equations" {
    const boundaries = [_]f64{ 0, 1, 2 };
    const leaf_area = [_]f64{ 2, 2 };
    try std.testing.expectApproxEqAbs(1.5, try canopy_photosynthesis.cuttingHeightForLeafAreaRemoval(0.25, &boundaries, &leaf_area), 1e-15);
    const direct_cut = try canopy_photosynthesis.layerHarvestRetention(1, 2, 1.5, false, false, 0, 0.8);
    try std.testing.expectApproxEqAbs(0.5, direct_cut.height_below_cut_fraction, 1e-15);
    try std.testing.expectApproxEqAbs(0.6, direct_cut.remaining_fraction, 1e-15);
    const thinned = try canopy_photosynthesis.layerHarvestRetention(1, 2, 1.5, false, true, 0.25, 0.8);
    try std.testing.expectApproxEqAbs(0.75, thinned.remaining_fraction, 1e-15);
    try std.testing.expectApproxEqAbs(0.9, thinned.unexported_fraction, 1e-15);
    try std.testing.expectApproxEqAbs(1.0 / 12.0, try canopy_photosynthesis.grazingCarbonDemandGPerH(true, 10, 0.1, 2, 0, 1, 20, 10), 1e-15);
    try std.testing.expectApproxEqAbs(1.0 / 12.0, try canopy_photosynthesis.grazingCarbonDemandGPerH(false, 10, 0.1, 1, 4, 0.5, 20, 10), 1e-15);
}

test "GROSUB grazing cascade carries unmet demand through organ order" {
    const allocation = try canopy_photosynthesis.allocateGrazingDemand(1, 0.5, 0.3, 0.2, 0.25, 0.1, .{ .leaf_carbon_g = 2, .sheath_carbon_g = 1, .husk_carbon_g = 1, .ear_carbon_g = 1, .grain_carbon_g = 1, .stalk_carbon_g = 2, .reserve_carbon_g = 2 });
    const physical_removal = allocation.structural_leaf_carbon_g + allocation.structural_sheath_carbon_g + allocation.husk_carbon_g + allocation.ear_carbon_g + allocation.grain_carbon_g + allocation.stalk_carbon_g + allocation.reserve_carbon_g + allocation.mobile_carbon_g;
    try std.testing.expectApproxEqAbs(1.0, physical_removal, 1e-14);
    try std.testing.expectEqual(0, allocation.unmet_carbon_g);
    try std.testing.expect(allocation.symbiont_mobile_carbon_g > 0);
}

test "GROSUB grazing demand retains ZEROP gate and animal insect drivers" {
    const animal = try canopy_photosynthesis.sourceOrderGrazingCarbonDemandGPerH(true, 2, 3, 10, 20, 0.5, 5, 10, 1.0e-12);
    const insect = try canopy_photosynthesis.sourceOrderGrazingCarbonDemandGPerH(false, 2, 3, 10, 20, 0.5, 5, 10, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.625), animal, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.625), insect, 1.0e-15);
    try std.testing.expectEqual(@as(f64, 0), try canopy_photosynthesis.sourceOrderGrazingCarbonDemandGPerH(true, 2, 3, 10, 20, 0.5, 5, 1.0e-12, 1.0e-12));
}

test "GROSUB additional nonfoliar redistribution exposes source overdraw" {
    const source = try canopy_photosynthesis.sourceOrderAdditionalGrazingRemoval(1, 0.8, 0.5);
    try std.testing.expectEqual(@as(f64, 1.3), source.next_total_removed_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), source.next_unmet_carbon_g_c);
}

test "GROSUB top-down grazing node selector characterizes source skip" {
    const partial = try canopy_photosynthesis.sourceOrderGrazingNodeRemoval(4, 1);
    try std.testing.expectEqual(@as(f64, 0.75), partial.remaining_fraction);
    try std.testing.expectEqual(@as(f64, 0), partial.remaining_branch_layer_demand_g_c);

    const skipped = try canopy_photosynthesis.sourceOrderGrazingNodeRemoval(1, 1);
    try std.testing.expectEqual(@as(f64, 1), skipped.remaining_fraction);
    try std.testing.expectEqual(@as(f64, 1), skipped.remaining_branch_layer_demand_g_c);
}

test "GROSUB branch-layer leaf demand retains ZEROP2 gate" {
    try std.testing.expectEqual(
        @as(f64, 0.5),
        try canopy_photosynthesis.sourceOrderBranchLayerLeafDemand(4, 2, 1, 1.0e-12),
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        try canopy_photosynthesis.sourceOrderBranchLayerLeafDemand(1.0e-12, 2, 1, 1.0e-12),
    );
}

test "GROSUB reproductive harvest conserves export litter and retained C N P" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{0};
    const sample_counts = [_]usize{};
    var state = try canopy_photosynthesis.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    state.branch_husk_carbon_g[0] = 2;
    state.branch_husk_nitrogen_g[0] = 0.2;
    state.branch_husk_phosphorus_g[0] = 0.02;
    state.branch_ear_carbon_g[0] = 3;
    state.branch_ear_nitrogen_g[0] = 0.3;
    state.branch_ear_phosphorus_g[0] = 0.03;
    state.branch_grain_carbon_g[0] = 5;
    state.branch_grain_nitrogen_g[0] = 0.5;
    state.branch_grain_phosphorus_g[0] = 0.05;
    state.branch_potential_seed_site_count[0] = 100;
    state.branch_seed_count[0] = 80;
    state.branch_individual_seed_carbon_g[0] = 0.1;
    const retention = try canopy_photosynthesis.reproductiveRetention(false, true, false, 0.2, 0.5, 2, 3, 5, 0, 0, 0);
    const result = try canopy_photosynthesis.harvestReproductiveOrgans(&state, 0, retention);
    const retained_c = state.branch_husk_carbon_g[0] + state.branch_ear_carbon_g[0] + state.branch_grain_carbon_g[0];
    try std.testing.expectApproxEqAbs(10.0, retained_c + result.products.ecosystem_export.carbon_g + result.products.litter.carbon_g, 1e-14);
    try std.testing.expectApproxEqAbs(1.0, result.harvested_grain.carbon_g, 1e-14);
    try std.testing.expectApproxEqAbs(80, state.branch_potential_seed_site_count[0], 1e-14);
    try std.testing.expectApproxEqAbs(64, state.branch_seed_count[0], 1e-14);
    try std.testing.expectEqual(0.1, state.branch_individual_seed_carbon_g[0]);
}

test "GROSUB reproductive retention preserves non-grazing source branches" {
    const direct_cut = try canopy_photosynthesis.sourceOrderReproductiveRetention(.{
        .grazing = false,
        .reproductive_organs_reached_by_cut = true,
        .grain_or_pruning = false,
        .thinning_fraction = 0,
        .harvested_nonfoliar_fraction = 0.6,
        .total_husk_carbon_g_c = 2,
        .total_ear_carbon_g_c = 3,
        .total_grain_carbon_g_c = 5,
        .grazed_husk_carbon_g_c = 0,
        .grazed_ear_carbon_g_c = 0,
        .grazed_grain_carbon_g_c = 0,
        .plant_presence_threshold_g_c = 1.0e-12,
    });
    try std.testing.expectEqual(@as(f64, 0.4), direct_cut.grain_remaining);
    try std.testing.expectEqual(direct_cut.grain_remaining, direct_cut.grain_unexported);

    var thinned = canopy_photosynthesis.SourceOrderReproductiveRetentionInput{
        .grazing = false,
        .reproductive_organs_reached_by_cut = false,
        .grain_or_pruning = true,
        .thinning_fraction = 0.25,
        .harvested_nonfoliar_fraction = 0.6,
        .total_husk_carbon_g_c = 2,
        .total_ear_carbon_g_c = 3,
        .total_grain_carbon_g_c = 5,
        .grazed_husk_carbon_g_c = 0,
        .grazed_ear_carbon_g_c = 0,
        .grazed_grain_carbon_g_c = 0,
        .plant_presence_threshold_g_c = 1.0e-12,
    };
    const grain_harvest = try canopy_photosynthesis.sourceOrderReproductiveRetention(thinned);
    try std.testing.expectEqual(@as(f64, 0.75), grain_harvest.grain_remaining);
    try std.testing.expectEqual(@as(f64, 0.85), grain_harvest.grain_unexported);

    thinned.grain_or_pruning = false;
    const above_cut = try canopy_photosynthesis.sourceOrderReproductiveRetention(thinned);
    try std.testing.expectEqual(@as(f64, 0.75), above_cut.grain_remaining);
    try std.testing.expectEqual(above_cut.grain_remaining, above_cut.grain_unexported);
}

test "GROSUB reproductive grazing uses strict plant presence threshold" {
    const retention = try canopy_photosynthesis.sourceOrderReproductiveRetention(.{
        .grazing = true,
        .reproductive_organs_reached_by_cut = false,
        .grain_or_pruning = false,
        .thinning_fraction = 0,
        .harvested_nonfoliar_fraction = 0,
        .total_husk_carbon_g_c = 1.0e-12,
        .total_ear_carbon_g_c = 2.0e-12,
        .total_grain_carbon_g_c = 4,
        .grazed_husk_carbon_g_c = 1.0e-12,
        .grazed_ear_carbon_g_c = 1.0e-12,
        .grazed_grain_carbon_g_c = 6,
        .plant_presence_threshold_g_c = 1.0e-12,
    });
    try std.testing.expectEqual(@as(f64, 1), retention.husk_remaining);
    try std.testing.expectEqual(@as(f64, 0.5), retention.ear_remaining);
    try std.testing.expectEqual(@as(f64, 0), retention.grain_remaining);
    try std.testing.expectEqual(retention.grain_remaining, retention.grain_unexported);
}
