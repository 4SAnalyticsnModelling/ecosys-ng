//! Tests for `plant_harvest_runtime.zig`.
//!
//! Extracted verbatim so the module beside it contains only the model
//! code. Tests that use private declarations of that module stay there,
//! since a sibling file can only reach `pub` declarations.

const builtin = @import("builtin");
const canopy = @import("../canopy/photosynthesis/photosynthesis.zig");
const canopy_biochemistry = @import("../canopy/photosynthesis/biochemistry.zig");
const canopy_layers = @import("../canopy/radiation/layer_distribution.zig");
const canopy_structure = @import("../canopy/morphology/structure.zig");
const carbon_exchange = @import("../canopy/photosynthesis/carbon_exchange.zig");
const dormancy = @import("../plant/lifecycle/dormancy.zig");
const grazing_manure = @import("../management/grazing_manure.zig");
const grid_module = @import("../state/grid.zig");
const growth_stages = @import("../plant/lifecycle/growth_stages.zig");
const litter_partition = @import("../plant/partition/litter.zig");
const management = @import("../management/plant_management.zig");
const phenology = @import("../plant/lifecycle/phenology.zig");
const root_disturbance = @import("../plant/root/plant_root_disturbance.zig");
const root_litter_ledger = @import("../plant/root/plant_root_litter_ledger.zig");
const root_litterfall = @import("../plant/root/plant_root_litterfall.zig");
const root_system = @import("../plant/root/plant_root_system.zig");
const shoot_litter_bridge = @import("../plant/growth/shoot_litter_bridge.zig");
const soil_organic = @import("../soil/organic/initialization.zig");
const spring_reproductive_litterfall = @import("../plant/growth/spring_reproductive_litterfall.zig");
const std = @import("std");
const surface_nutrients = @import("../soil/biogeochemistry/organic_matter_fire_exchange.zig");
const symbiotic_fixation = @import("../canopy/symbiosis/plant_symbiotic_fixation.zig");
const plant_harvest_runtime = @import("../management/plant_harvest_runtime.zig");
test "GROSUB forest self thinning retains exact density law and runtime units" {
    try std.testing.expectEqual(@as(f64, 0), try plant_harvest_runtime.forestSelfThinningFraction(0, 10));
    try std.testing.expectEqual(@as(f64, 0), try plant_harvest_runtime.forestSelfThinningFraction(0.25, 0.05));
    // At 0.25 m PPQ is exactly 0.1 plants m-2.
    try std.testing.expectApproxEqAbs(@as(f64, 0.09), try plant_harvest_runtime.forestSelfThinningFraction(0.25, 1), 1.0e-15);
    try std.testing.expectError(error.InvalidForestSelfThinningInput, plant_harvest_runtime.forestSelfThinningFraction(-0.1, 1));
}

test "GROSUB first substep resets only hourly disturbance products" {
    const current: plant_harvest_runtime.HourlyDisturbanceReset = .{
        .previous_cumulative_harvest_carbon_g_c = 3,
        .manure_organic_carbon_g_c = .{ 1, 2, 3, 4 },
        .manure_organic_nitrogen_g_n = .{ 5, 6, 7, 8 },
        .manure_organic_phosphorus_g_p = .{ 9, 10, 11, 12 },
        .manure_inorganic_nitrogen_g_n = 13,
        .manure_inorganic_phosphorus_g_p = 14,
    };
    const unchanged = try plant_harvest_runtime.sourceOrderHourlyDisturbanceReset(false, 20, current);
    try std.testing.expectEqualDeep(current, unchanged);
    const reset = try plant_harvest_runtime.sourceOrderHourlyDisturbanceReset(true, 20, current);
    try std.testing.expectEqual(@as(f64, 20), reset.previous_cumulative_harvest_carbon_g_c);
    try std.testing.expectEqual([4]f64{ 0, 0, 0, 0 }, reset.manure_organic_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), reset.manure_inorganic_nitrogen_g_n);
}

test "GROSUB forest self thinning selector preserves monthly noon and event gates" {
    try std.testing.expect(try plant_harvest_runtime.sourceOrderForestSelfThinningIsEnabled(30, 12, 12.75, 1, 2, -1));
    try std.testing.expect(try plant_harvest_runtime.sourceOrderForestSelfThinningIsEnabled(360, 12, 12.75, 1, 2, 4));
    try std.testing.expect(try plant_harvest_runtime.sourceOrderForestSelfThinningIsEnabled(30, 12, 12.75, 1, 2, 6));
    try std.testing.expect(!try plant_harvest_runtime.sourceOrderForestSelfThinningIsEnabled(30, 12, 12.75, 1, 2, 0));
    try std.testing.expect(!try plant_harvest_runtime.sourceOrderForestSelfThinningIsEnabled(29, 12, 12.75, 1, 2, -1));
    try std.testing.expect(!try plant_harvest_runtime.sourceOrderForestSelfThinningIsEnabled(30, 13, 12.75, 1, 2, -1));
    try std.testing.expect(!try plant_harvest_runtime.sourceOrderForestSelfThinningIsEnabled(30, 12, 12.75, 0, 2, -1));
    try std.testing.expect(!try plant_harvest_runtime.sourceOrderForestSelfThinningIsEnabled(30, 12, 12.75, 1, 1, -1));
}

test "GROSUB pruning multiplies the persistent canopy clumping factor" {
    try std.testing.expectApproxEqAbs(@as(f64, 0.48), try plant_harvest_runtime.prunedClumpingFactor(0.8, 0.6), 1.0e-15);
    try std.testing.expectError(error.InvalidPruningClumpingFraction, plant_harvest_runtime.prunedClumpingFactor(0.8, -0.1));
}

test "GROSUB negative HVST interpolates combined canopy leaf area across runtime layers" {
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.5),
        try plant_harvest_runtime.cuttingHeightFromLeafAreaRemoval(0.5, &.{ 0, 1, 3 }, &.{ 2, 4 }, 1.0e-12),
        1.0e-15,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        try plant_harvest_runtime.cuttingHeightFromLeafAreaRemoval(1, &.{ 0, 1, 3 }, &.{ 2, 4 }, 1.0e-12),
    );
    try std.testing.expectError(
        error.InvalidLeafAreaHarvestGeometry,
        plant_harvest_runtime.cuttingHeightFromLeafAreaRemoval(1.1, &.{ 0, 1 }, &.{1}, 1.0e-12),
    );
}

test "GROSUB cutting height uses authoritative combined canopy leaf area" {
    const exact = try plant_harvest_runtime.sourceOrderCuttingHeightFromLeafAreaRemoval(
        0.5,
        8,
        &.{ 0, 1, 3 },
        &.{ 2, 4 },
        1.0e-12,
    );
    try std.testing.expectEqual(@as(f64, 2), exact);
    const recomputed = try plant_harvest_runtime.cuttingHeightFromLeafAreaRemoval(0.5, &.{ 0, 1, 3 }, &.{ 2, 4 }, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 1.5), recomputed);
}

test "GROSUB aboveground disturbance dispatch and population update are exact" {
    try std.testing.expect(try plant_harvest_runtime.sourceOrderAbovegroundDisturbanceIsEnabled(4, 3, 12.5));
    try std.testing.expect(try plant_harvest_runtime.sourceOrderAbovegroundDisturbanceIsEnabled(6, 3, 12.5));
    try std.testing.expect(try plant_harvest_runtime.sourceOrderAbovegroundDisturbanceIsEnabled(2, 12, 12.5));
    try std.testing.expect(!try plant_harvest_runtime.sourceOrderAbovegroundDisturbanceIsEnabled(2, 11, 12.5));
    try std.testing.expect(!try plant_harvest_runtime.sourceOrderAbovegroundDisturbanceIsEnabled(-1, 12, 12.5));

    const current: plant_harvest_runtime.PopulationAfterDisturbance = .{
        .living_population_per_m2 = 4,
        .living_population_count = 40,
        .standing_dead_population_count = 10,
    };
    const thinned = try plant_harvest_runtime.sourceOrderPopulationAfterDisturbance(false, 0.25, current, 7, 10);
    try std.testing.expectEqual(@as(f64, 3), thinned.living_population_per_m2);
    try std.testing.expectEqual(@as(f64, 30), thinned.living_population_count);
    try std.testing.expectEqual(@as(f64, 7.5), thinned.standing_dead_population_count);
    const reseeded = try plant_harvest_runtime.sourceOrderPopulationAfterDisturbance(true, 0.25, current, 7, 10);
    try std.testing.expectEqual(@as(f64, 7), reseeded.living_population_per_m2);
    try std.testing.expectEqual(@as(f64, 70), reseeded.living_population_count);
    try std.testing.expectEqual(@as(f64, 70), reseeded.standing_dead_population_count);
}

test "GROSUB harvest litter uses organ-specific runtime kinetics conservatively" {
    const foliar: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 0, 1, 0, 0 },
        .phosphorus = .{ 0, 0, 1, 0 },
    };
    const nonfoliar: litter_partition.ElementFractions = .{
        .carbon = .{ 0, 1, 0, 0 },
        .nitrogen = .{ 0, 0, 1, 0 },
        .phosphorus = .{ 0, 0, 0, 1 },
    };
    const woody: litter_partition.ElementFractions = .{
        .carbon = .{ 0, 0, 1, 0 },
        .nitrogen = .{ 0, 0, 0, 1 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const products: plant_harvest_runtime.ProductLedger = .{
        .foliar = .{ .litter = .{ .carbon_g = 2, .nitrogen_g = 0.2, .phosphorus_g = 0.02 } },
        .nonfoliar = .{ .litter = .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 } },
        .woody = .{ .litter = .{ .carbon_g = 5, .nitrogen_g = 0.5, .phosphorus_g = 0.05 } },
    };
    const result = try plant_harvest_runtime.harvestLitterToKinetics(products, foliar, foliar, nonfoliar, woody);
    try std.testing.expectEqual([4]f64{ 2, 3, 0, 0 }, result.nonwoody_carbon_g);
    try std.testing.expectEqual([4]f64{ 0, 0, 5, 0 }, result.woody_carbon_g);
    var carbon_g_c: f64 = 0;
    var nitrogen_g_n: f64 = 0;
    var phosphorus_g_p: f64 = 0;
    for (0..4) |kinetic| {
        carbon_g_c += result.nonwoody_carbon_g[kinetic] + result.woody_carbon_g[kinetic];
        nitrogen_g_n += result.nonwoody_nitrogen_g[kinetic] + result.woody_nitrogen_g[kinetic];
        phosphorus_g_p += result.nonwoody_phosphorus_g[kinetic] + result.woody_phosphorus_g[kinetic];
    }
    try std.testing.expectApproxEqAbs(@as(f64, 10), carbon_g_c, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), nitrogen_g_n, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), phosphorus_g_p, 1.0e-15);
}

test "GROSUB automatic deciduous annual harvest fires once at reproductive turnover" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{0}, &.{});
    defer state.deinit();
    var layers = try canopy_layers.State.init(std.testing.allocator, 1, 1, 1, 1, 1, &state);
    defer layers.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    var plant_state = try phenology.State.init(std.testing.allocator, 1, 1);
    defer plant_state.deinit();
    plant_state.active[0] = true;
    plant_state.lifecycle_initialized[0] = true;
    var growth = try growth_stages.State.init(std.testing.allocator, &.{1});
    defer growth.deinit();
    state.branch_grain_carbon_g[0] = 5;
    state.branch_grain_nitrogen_g[0] = 0.5;
    state.branch_grain_phosphorus_g[0] = 0.05;
    state.plant_population_per_m2[0] = 2;
    const science = [_]plant_harvest_runtime.ScienceParameters{.{ .carbon_woody_fraction = .{ 0, 1 }, .leaf_nitrogen_woody_fraction = .{ 0, 1 }, .sheath_nitrogen_woody_fraction = .{ 0, 1 }, .leaf_phosphorus_woody_fraction = .{ 0, 1 }, .sheath_phosphorus_woody_fraction = .{ 0, 1 } }};
    var ledgers = [_]plant_harvest_runtime.ProductLedger{.{}};
    const population = [_]f64{4};
    const area = [_]f64{2};
    const root_woody = [_]f64{0};
    var automatic_dates = [_]management.PackedDate{.{ .day = 1, .month = 1, .year = 0 }};
    var context: plant_harvest_runtime.Context = .{
        .canopy_state = &state,
        .canopy_layer_state = &layers,
        .branch_development = &development,
        .science_by_plant = &science,
        .products_by_plant = &ledgers,
        .leaf_area_presence_tolerance_m2 = 1e-12,
        .plant_phenology = &plant_state,
        .growth_stages = &growth,
        .reseed_population_per_m2_by_plant = &population,
        .cell_area_m2_by_cell = &area,
        .root_woody_fraction_by_plant = &root_woody,
        .automatic_harvest_date_by_plant = &automatic_dates,
    };
    const current_date: management.PackedDate = .{ .day = 17, .month = 9, .year = 2004 };
    try std.testing.expectEqual(@as(usize, 1), try plant_harvest_runtime.applyAutomaticSelfSeedingHarvests(&context, &.{true}, &.{0}, &.{1}, current_date));
    try std.testing.expectApproxEqAbs(@as(f64, 5), state.plant_seed_storage_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 4), state.plant_population_per_m2[0], 1e-12);
    try std.testing.expect(plant_state.reseed_pending[0]);
    try std.testing.expectEqual(current_date, automatic_dates[0]);
    try std.testing.expectEqual(@as(usize, 0), try plant_harvest_runtime.applyAutomaticSelfSeedingHarvests(&context, &.{true}, &.{0}, &.{1}, current_date));
}

test "GROSUB perennial start-of-season residue is conserved before reconstruction" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{1});
    defer state.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    var partitions = try litter_partition.State.init(std.testing.allocator, 1);
    defer partitions.deinit();
    @memset(partitions.by_plant_and_organ, .{ .carbon = .{ 0.1, 0.2, 0.3, 0.4 }, .nitrogen = .{ 0.1, 0.2, 0.3, 0.4 }, .phosphorus = .{ 0.1, 0.2, 0.3, 0.4 } });
    state.sample_leaf_area_m2[0] = 1;
    state.sample_leaf_carbon_g[0] = 2;
    state.sample_leaf_nitrogen_g[0] = 0.2;
    state.sample_leaf_phosphorus_g[0] = 0.02;
    state.node_leaf_area_m2[0] = 1;
    state.node_leaf_carbon_g[0] = 2;
    state.node_leaf_nitrogen_g[0] = 0.2;
    state.node_leaf_phosphorus_g[0] = 0.02;
    state.branch_leaf_area_m2[0] = 1;
    state.branch_leaf_carbon_g[0] = 2;
    state.branch_leaf_nitrogen_g[0] = 0.2;
    state.branch_leaf_phosphorus_g[0] = 0.02;
    state.node_sheath_carbon_g[0] = 1;
    state.node_sheath_nitrogen_g[0] = 0.1;
    state.node_sheath_phosphorus_g[0] = 0.01;
    state.branch_sheath_carbon_g[0] = 1;
    state.branch_sheath_nitrogen_g[0] = 0.1;
    state.branch_sheath_phosphorus_g[0] = 0.01;
    state.branch_husk_carbon_g[0] = 1;
    state.branch_ear_carbon_g[0] = 1;
    state.branch_grain_carbon_g[0] = 1;
    state.branch_stalk_carbon_g[0] = 4;
    state.branch_stalk_nitrogen_g[0] = 0.4;
    state.branch_stalk_phosphorus_g[0] = 0.04;
    state.branch_reserve_carbon_g[0] = 2;
    state.branch_reserve_nitrogen_g[0] = 0.2;
    state.branch_reserve_phosphorus_g[0] = 0.02;
    const science = [_]plant_harvest_runtime.ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.25, 0.75 }, .sheath_nitrogen_woody_fraction = .{ 0.25, 0.75 }, .leaf_phosphorus_woody_fraction = .{ 0.25, 0.75 }, .sheath_phosphorus_woody_fraction = .{ 0.25, 0.75 } }};
    var ledgers = [_]plant_harvest_runtime.ProductLedger{.{}};
    var context: plant_harvest_runtime.Context = .{ .canopy_state = &state, .branch_development = &development, .science_by_plant = &science, .products_by_plant = &ledgers, .leaf_area_presence_tolerance_m2 = 1e-12, .root_litter_partition = &partitions };
    var before_failure = try state.clone();
    defer before_failure.deinit();
    const products_before_failure = ledgers;
    partitions.by_plant_and_organ[@intFromEnum(litter_partition.Organ.stalk)].carbon[3] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidPlantLitterFraction, plant_harvest_runtime.applyStartOfSeasonResidue(&context, 0, 0, 1));
    inline for (@typeInfo(canopy.State).@"struct".fields) |field| if (field.type == []f64)
        try std.testing.expectEqualSlices(f64, @field(before_failure, field.name), @field(state, field.name));
    try std.testing.expectEqualDeep(products_before_failure, ledgers);
    partitions.by_plant_and_organ[@intFromEnum(litter_partition.Organ.stalk)].carbon = .{ 0.1, 0.2, 0.3, 0.4 };
    try plant_harvest_runtime.applyStartOfSeasonResidue(&context, 0, 0, 1);
    try std.testing.expectApproxEqAbs(@as(f64, 4), state.plant_standing_dead_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2), state.plant_seed_storage_carbon_g[0], 1e-12);
    var direct_litter_carbon_g_c: f64 = 0;
    for (ledgers[0].direct_litter.nonwoody_carbon_g) |value| direct_litter_carbon_g_c += value;
    const litter_carbon_g_c = ledgers[0].foliar.litter.carbon_g + ledgers[0].nonfoliar.litter.carbon_g + ledgers[0].woody.litter.carbon_g + direct_litter_carbon_g_c;
    try std.testing.expectApproxEqAbs(@as(f64, 6), litter_carbon_g_c, 1e-12);
    var standing_kinetic_carbon_g_c: f64 = 0;
    for (state.plant_standing_dead_carbon_by_kinetic_g[0..4]) |value| standing_kinetic_carbon_g_c += value;
    try std.testing.expectApproxEqAbs(@as(f64, 4), standing_kinetic_carbon_g_c, 1e-12);
    try std.testing.expectEqual(@as(f64, 0), state.branch_stalk_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_reserve_carbon_g[0]);
}

test "GROSUB whole-plant death conserves shoot storage symbiont and standing dead carbon" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{1});
    defer state.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    var partitions = try litter_partition.State.init(std.testing.allocator, 1);
    defer partitions.deinit();
    @memset(partitions.by_plant_and_organ, .{
        .carbon = .{ 0.1, 0.2, 0.3, 0.4 },
        .nitrogen = .{ 0.1, 0.2, 0.3, 0.4 },
        .phosphorus = .{ 0.1, 0.2, 0.3, 0.4 },
    });
    state.plant_seed_storage_carbon_g[0] = 4;
    state.branch_stalk_carbon_g[0] = 3;
    state.branch_reserve_carbon_g[0] = 2;
    state.branch_mobile_carbon_g[0] = 1;
    state.node_c4_mesophyll_nonstructural_carbon_g[0] = 1;
    state.branch_symbiont_mobile_carbon_g[0] = 0.5;
    state.branch_symbiont_structural_carbon_g[0] = 0.5;
    state.branch_husk_carbon_g[0] = 1;
    const science = [_]plant_harvest_runtime.ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.25, 0.75 }, .sheath_nitrogen_woody_fraction = .{ 0.25, 0.75 }, .leaf_phosphorus_woody_fraction = .{ 0.25, 0.75 }, .sheath_phosphorus_woody_fraction = .{ 0.25, 0.75 } }};
    var ledgers = [_]plant_harvest_runtime.ProductLedger{.{}};
    var context: plant_harvest_runtime.Context = .{
        .canopy_state = &state,
        .branch_development = &development,
        .science_by_plant = &science,
        .products_by_plant = &ledgers,
        .leaf_area_presence_tolerance_m2 = 1e-12,
        .root_litter_partition = &partitions,
    };

    try plant_harvest_runtime.applyWholePlantMortalityResidue(&context, 0);
    var direct_litter_carbon_g_c: f64 = 0;
    var direct_woody_carbon_g_c: f64 = 0;
    var direct_nonwoody_carbon_g_c: f64 = 0;
    for (ledgers[0].direct_litter.woody_carbon_g, ledgers[0].direct_litter.nonwoody_carbon_g) |woody, nonwoody| {
        direct_woody_carbon_g_c += woody;
        direct_nonwoody_carbon_g_c += nonwoody;
        direct_litter_carbon_g_c += woody + nonwoody;
    }
    const litter_carbon_g_c = direct_litter_carbon_g_c + ledgers[0].nonstructural.litter.carbon_g +
        ledgers[0].foliar.litter.carbon_g +
        ledgers[0].nonfoliar.litter.carbon_g +
        ledgers[0].woody.litter.carbon_g;
    try std.testing.expectApproxEqAbs(@as(f64, 5), state.plant_standing_dead_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 8), litter_carbon_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 13), state.plant_standing_dead_carbon_g[0] + litter_carbon_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), direct_woody_carbon_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), direct_nonwoody_carbon_g_c, 1e-12);
    try std.testing.expectEqual(@as(f64, 0), state.plant_seed_storage_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_stalk_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_reserve_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_mobile_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.node_c4_mesophyll_nonstructural_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_symbiont_mobile_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_symbiont_structural_carbon_g[0]);
}

test "GROSUB natural winter-annual branch death recovers mobile reserve and grain" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{0});
    defer state.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    var partitions = try litter_partition.State.init(std.testing.allocator, 1);
    defer partitions.deinit();
    @memset(partitions.by_plant_and_organ, .{
        .carbon = .{ 0.1, 0.2, 0.3, 0.4 },
        .nitrogen = .{ 0.1, 0.2, 0.3, 0.4 },
        .phosphorus = .{ 0.1, 0.2, 0.3, 0.4 },
    });
    state.branch_mobile_carbon_g[0] = 2;
    state.node_c4_mesophyll_nonstructural_carbon_g[0] = 1;
    state.branch_reserve_carbon_g[0] = 3;
    state.branch_grain_carbon_g[0] = 4;
    state.branch_stalk_carbon_g[0] = 5;
    state.branch_symbiont_structural_carbon_g[0] = 1;
    const science = [_]plant_harvest_runtime.ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.25, 0.75 }, .sheath_nitrogen_woody_fraction = .{ 0.25, 0.75 }, .leaf_phosphorus_woody_fraction = .{ 0.25, 0.75 }, .sheath_phosphorus_woody_fraction = .{ 0.25, 0.75 } }};
    var ledgers = [_]plant_harvest_runtime.ProductLedger{.{}};
    var context: plant_harvest_runtime.Context = .{ .canopy_state = &state, .branch_development = &development, .science_by_plant = &science, .products_by_plant = &ledgers, .leaf_area_presence_tolerance_m2 = 1e-12, .root_litter_partition = &partitions };

    try plant_harvest_runtime.applyNaturalDeadBranchResidue(&context, 0, 0, true);
    try std.testing.expectApproxEqAbs(@as(f64, 10), state.plant_seed_storage_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5), state.plant_standing_dead_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), ledgers[0].foliar.litter.carbon_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 16), state.plant_seed_storage_carbon_g[0] + state.plant_standing_dead_carbon_g[0] + ledgers[0].foliar.litter.carbon_g, 1e-12);
    try std.testing.expectEqual(@as(f64, 0), state.branch_mobile_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.node_c4_mesophyll_nonstructural_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_grain_carbon_g[0]);
}

test "non-grazing runtime callback conserves branch carbon through cutting" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{1};
    const sample_counts = [_]usize{1};
    var canopy_state = try canopy.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer canopy_state.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    canopy_state.sample_leaf_area_m2[0] = 2;
    canopy_state.sample_layer_lower_height_m[0] = 0;
    canopy_state.sample_layer_upper_height_m[0] = 2;
    canopy_state.sample_leaf_carbon_g[0] = 4;
    canopy_state.sample_leaf_nitrogen_g[0] = 0.4;
    canopy_state.sample_leaf_phosphorus_g[0] = 0.04;
    canopy_state.node_leaf_area_m2[0] = 2;
    canopy_state.node_leaf_carbon_g[0] = 4;
    canopy_state.node_leaf_nitrogen_g[0] = 0.4;
    canopy_state.node_leaf_phosphorus_g[0] = 0.04;
    canopy_state.branch_leaf_area_m2[0] = 2;
    canopy_state.branch_leaf_carbon_g[0] = 4;
    canopy_state.branch_leaf_nitrogen_g[0] = 0.4;
    canopy_state.branch_leaf_phosphorus_g[0] = 0.04;
    canopy_state.node_height_m[0] = 2;
    canopy_state.node_sheath_height_m[0] = 1;
    canopy_state.node_sheath_carbon_g[0] = 2;
    canopy_state.node_sheath_nitrogen_g[0] = 0.2;
    canopy_state.node_sheath_phosphorus_g[0] = 0.02;
    canopy_state.branch_sheath_carbon_g[0] = 2;
    canopy_state.branch_sheath_nitrogen_g[0] = 0.2;
    canopy_state.branch_sheath_phosphorus_g[0] = 0.02;
    canopy_state.node_internode_length_m[0] = 2;
    canopy_state.node_internode_carbon_g[0] = 4;
    canopy_state.node_internode_nitrogen_g[0] = 0.4;
    canopy_state.node_internode_phosphorus_g[0] = 0.04;
    canopy_state.branch_stalk_carbon_g[0] = 4;
    canopy_state.branch_stalk_nitrogen_g[0] = 0.4;
    canopy_state.branch_stalk_phosphorus_g[0] = 0.04;
    canopy_state.branch_reserve_carbon_g[0] = 1;
    canopy_state.branch_reserve_nitrogen_g[0] = 0.1;
    canopy_state.branch_reserve_phosphorus_g[0] = 0.01;
    canopy_state.branch_mobile_carbon_g[0] = 2;
    canopy_state.branch_mobile_nitrogen_g[0] = 0.2;
    canopy_state.branch_mobile_phosphorus_g[0] = 0.02;
    const science = [_]plant_harvest_runtime.ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .sheath_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .leaf_phosphorus_woody_fraction = .{ 0.1, 0.9 }, .sheath_phosphorus_woody_fraction = .{ 0.1, 0.9 } }};
    var ledgers = [_]plant_harvest_runtime.ProductLedger{.{}};
    var context: plant_harvest_runtime.Context = .{ .canopy_state = &canopy_state, .branch_development = &development, .science_by_plant = &science, .products_by_plant = &ledgers, .leaf_area_presence_tolerance_m2 = 1.0e-12 };
    const event: management.HarvestEvent = .{ .date = .{ .day = 1, .month = 1, .year = 9999 }, .kind = .above_ground, .termination = .retain, .cutting_height_m_or_lai_fraction = 1, .thinning_fraction_or_consumption_rate = 0, .harvested_fraction = .{ .leaf = 1, .nonfoliar = 1, .woody = 1, .standing_dead = 0 }, .ecosystem_export_fraction = .{ .leaf = 0.8, .nonfoliar = 0.8, .woody = 0.8, .standing_dead = 0 } };
    try plant_harvest_runtime.applyEvent(&context, 0, event);
    const remaining_c = canopy_state.branch_leaf_carbon_g[0] + canopy_state.branch_sheath_carbon_g[0] + canopy_state.branch_stalk_carbon_g[0] + canopy_state.branch_reserve_carbon_g[0] + canopy_state.branch_mobile_carbon_g[0];
    const products_c = ledgers[0].nonstructural.ecosystem_export.carbon_g + ledgers[0].nonstructural.litter.carbon_g + ledgers[0].foliar.ecosystem_export.carbon_g + ledgers[0].foliar.litter.carbon_g + ledgers[0].nonfoliar.ecosystem_export.carbon_g + ledgers[0].nonfoliar.litter.carbon_g + ledgers[0].woody.ecosystem_export.carbon_g + ledgers[0].woody.litter.carbon_g;
    try std.testing.expectApproxEqAbs(13.0, remaining_c + products_c, 1e-12);
    try std.testing.expectApproxEqAbs(2.0, canopy_state.branch_leaf_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(1.0, canopy_state.branch_sheath_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(2.0, canopy_state.branch_stalk_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(1.0, canopy_state.branch_mobile_carbon_g[0], 1e-15);

    canopy_state.branch_grain_carbon_g[0] = 3;
    canopy_state.branch_grain_nitrogen_g[0] = 0.3;
    canopy_state.branch_grain_phosphorus_g[0] = 0.03;
    var high_cut = event;
    high_cut.cutting_height_m_or_lai_fraction = 3;
    try plant_harvest_runtime.applyEvent(&context, 0, high_cut);
    try std.testing.expectApproxEqAbs(
        3,
        canopy_state.branch_grain_carbon_g[0],
        1e-15,
    );
}

test "runtime callback rejects grazing approximation" {
    const event: management.HarvestEvent = .{ .date = .{ .day = 1, .month = 1, .year = 9999 }, .kind = .animal_grazing, .termination = .retain, .cutting_height_m_or_lai_fraction = 0, .thinning_fraction_or_consumption_rate = 0, .harvested_fraction = .{ .leaf = 0, .nonfoliar = 0, .woody = 0, .standing_dead = 0 }, .ecosystem_export_fraction = .{ .leaf = 0, .nonfoliar = 0, .woody = 0, .standing_dead = 0 } };
    var context: plant_harvest_runtime.Context = undefined;
    try std.testing.expectError(error.GrazingRequiresDemandDrivenKernel, plant_harvest_runtime.applyEvent(&context, 0, event));
}

test "GROSUB grazing distributes host and symbiont pools by branch leaf sheath mass" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{2}, &.{ 1, 1 }, &.{ 1, 1 });
    defer state.deinit();
    var layers = try canopy_layers.State.init(std.testing.allocator, 1, 1, 1, 1, 1, &state);
    defer layers.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 2);
    defer development.deinit();
    for (0..2) |branch| {
        const leaf_g_c: f64 = if (branch == 0) 3 else 1;
        state.sample_leaf_area_m2[branch] = leaf_g_c;
        state.sample_exposed_leaf_area_m2[branch] = leaf_g_c;
        state.sample_leaf_carbon_g[branch] = leaf_g_c;
        state.node_leaf_area_m2[branch] = leaf_g_c;
        state.node_leaf_carbon_g[branch] = leaf_g_c;
        state.branch_leaf_area_m2[branch] = leaf_g_c;
        state.branch_leaf_carbon_g[branch] = leaf_g_c;
        layers.node_leaf_area_m2[branch] = leaf_g_c;
        layers.node_leaf_carbon_g[branch] = leaf_g_c;
        state.branch_mobile_carbon_g[branch] = if (branch == 0) 1 else 3;
        state.branch_symbiont_mobile_carbon_g[branch] = 1;
        state.branch_symbiont_structural_carbon_g[branch] = 2;
    }
    state.plant_total_shoot_carbon_g[0] = 14;
    state.plant_uptake_growth_temperature_response[0] = 1;
    state.plant_mobile_carbon_concentration_g_per_g[0] = 1;
    state.plant_symbiont_mobile_carbon_concentration_g_per_g[0] = 1;
    const science = [_]plant_harvest_runtime.ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .sheath_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .leaf_phosphorus_woody_fraction = .{ 0.1, 0.9 }, .sheath_phosphorus_woody_fraction = .{ 0.1, 0.9 } }};
    var ledgers = [_]plant_harvest_runtime.ProductLedger{.{}};
    var context: plant_harvest_runtime.Context = .{ .canopy_state = &state, .canopy_layer_state = &layers, .branch_development = &development, .science_by_plant = &science, .products_by_plant = &ledgers, .leaf_area_presence_tolerance_m2 = 1e-12 };
    const removed_g_c = try plant_harvest_runtime.applyGrazingEvent(&context, 0, .{
        .date = .{ .day = 1, .month = 1, .year = 9999 },
        .kind = .animal_grazing,
        .termination = .retain,
        .cutting_height_m_or_lai_fraction = 96,
        .thinning_fraction_or_consumption_rate = 1,
        .harvested_fraction = .{ .leaf = 1, .nonfoliar = 0, .woody = 0, .standing_dead = 0 },
        .ecosystem_export_fraction = .{ .leaf = 1, .nonfoliar = 1, .woody = 1, .standing_dead = 0 },
    }, 14, 1);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), state.branch_mobile_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2.75), state.branch_mobile_carbon_g[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), state.branch_symbiont_mobile_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.75), state.branch_symbiont_mobile_carbon_g[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), state.branch_symbiont_structural_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), state.branch_symbiont_structural_carbon_g[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5), removed_g_c, 1e-12);
}

test "grazing publication commits manure mineral nutrients and export without examples" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{0}, &.{});
    defer state.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    var partitions = try litter_partition.State.init(std.testing.allocator, 1);
    defer partitions.deinit();
    const uniform: litter_partition.ElementFractions = .{
        .carbon = .{ 0.25, 0.25, 0.25, 0.25 },
        .nitrogen = .{ 0.25, 0.25, 0.25, 0.25 },
        .phosphorus = .{ 0.25, 0.25, 0.25, 0.25 },
    };
    @memset(partitions.by_plant_and_organ, uniform);
    const config = try @import("../core/config.zig").SimulationConfig.init(
        .{ .lon_count = 1, .lat_count = 1, .soil_layers = 1, .plant_populations = 1 },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-12, .max_nonlinear_iterations = 10 },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var surface = try soil_organic.State.init(std.testing.allocator, 1);
    defer surface.deinit();
    var nutrients = try surface_nutrients.State.init(std.testing.allocator, 1, soil_organic.substrate_count);
    defer nutrients.deinit();
    var daily_manure_c = [_]f64{0};
    var daily_manure_n = [_]f64{0};
    var daily_manure_p = [_]f64{0};
    var shoot_litter_carbon_g_c = [_]f64{0};
    var shoot_litter_nitrogen_g_n = [_]f64{0};
    var shoot_litter_phosphorus_g_p = [_]f64{0};
    var hourly_manure_products = [_]grazing_manure.Products{.{}};
    const science = [_]plant_harvest_runtime.ScienceParameters{.{ .carbon_woody_fraction = .{ 0, 1 }, .leaf_nitrogen_woody_fraction = .{ 0, 1 }, .sheath_nitrogen_woody_fraction = .{ 0, 1 }, .leaf_phosphorus_woody_fraction = .{ 0, 1 }, .sheath_phosphorus_woody_fraction = .{ 0, 1 } }};
    var ledgers = [_]plant_harvest_runtime.ProductLedger{.{}};
    ledgers[0].direct_litter.nonwoody_carbon_g[0] = 2;
    ledgers[0].direct_litter.nonwoody_nitrogen_g[0] = 0.2;
    ledgers[0].direct_litter.nonwoody_phosphorus_g[0] = 0.02;
    ledgers[0].manure = try grazing_manure.partition(.animal_grazing, .{ .carbon_g = 10, .nitrogen_g = 1, .phosphorus_g = 0.2 });
    const expected_hourly_manure = ledgers[0].manure;
    ledgers[0].standing_dead_export = .{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 };
    var context: plant_harvest_runtime.Context = .{
        .canopy_state = &state,
        .branch_development = &development,
        .science_by_plant = &science,
        .products_by_plant = &ledgers,
        .root_litter_partition = &partitions,
        .surface_organic_state = &surface,
        .surface_nutrient_state = &nutrients,
        .daily_manure_carbon_input_g_c = &daily_manure_c,
        .daily_manure_nitrogen_input_g_n = &daily_manure_n,
        .daily_manure_phosphorus_input_g_p = &daily_manure_p,
        .hourly_manure_products_by_plant = &hourly_manure_products,
        .shoot_litter_carbon_g_c_by_plant = &shoot_litter_carbon_g_c,
        .shoot_litter_nitrogen_g_n_by_plant = &shoot_litter_nitrogen_g_n,
        .shoot_litter_phosphorus_g_p_by_plant = &shoot_litter_phosphorus_g_p,
        .grid = &grid,
        .leaf_area_presence_tolerance_m2 = 1e-12,
    };
    const exported = try plant_harvest_runtime.publishPlantProducts(&context, 0);
    try std.testing.expectApproxEqAbs(@as(f64, 1), exported.carbon_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2), shoot_litter_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), shoot_litter_nitrogen_g_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.02), shoot_litter_phosphorus_g_p[0], 1e-12);
    var manure: canopy.ElementalMass = .{};
    for (0..4) |fraction| {
        const index = (2 * soil_organic.structural_fraction_count) + fraction;
        manure.carbon_g += surface.structural[index].carbon_g_c;
        manure.nitrogen_g += surface.structural[index].nitrogen_g_n;
        manure.phosphorus_g += surface.structural[index].phosphorus_g_p;
    }
    try std.testing.expectApproxEqAbs(@as(f64, 10), manure.carbon_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), manure.nitrogen_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), manure.phosphorus_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5 / 14.0), nutrients.pending_surface_ammonium_mol_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1 / 31.0), nutrients.pending_surface_phosphate_mol_p[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 10), daily_manure_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), daily_manure_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), daily_manure_p[0], 1e-12);
    try std.testing.expectEqual(
        expected_hourly_manure,
        hourly_manure_products[0],
    );
    try std.testing.expectEqual(plant_harvest_runtime.ProductLedger{}, ledgers[0]);

    const surface_before = try std.testing.allocator.dupe(soil_organic.ElementPool, surface.structural);
    defer std.testing.allocator.free(surface_before);
    shoot_litter_phosphorus_g_p[0] = std.math.floatMax(f64);
    ledgers[0].direct_litter.nonwoody_phosphorus_g[0] = std.math.floatMax(f64);
    try std.testing.expectError(
        error.InvalidPlantHarvestProduct,
        plant_harvest_runtime.publishPlantProducts(&context, 0),
    );
    try std.testing.expectEqualSlices(soil_organic.ElementPool, surface_before, surface.structural);
    try std.testing.expectEqual(
        std.math.floatMax(f64),
        ledgers[0].direct_litter.nonwoody_phosphorus_g[0],
    );
}

test "source-order tillage population reduction scales all population fields" {
    const result = try plant_harvest_runtime.sourceOrderTillagePopulationReduction(.{
        .hour_of_day = 12,
        .local_solar_noon_h = 12.75,
        .biomass_turnover_type = 0,
        .root_profile_type = 2,
        .current_day_of_year = 151,
        .current_year = 2001,
        .planting_day_of_year = 100,
        .planting_year = 2001,
        .tillage_code = 8,
        .is_first_plant_population = true,
        .remaining_fraction = 0.25,
        .zero_population_threshold = 1.0e-12,
        .state = .{
            .living_population_per_m2 = 8,
            .living_population_count = 40,
            .standing_dead_population_count = 12,
            .canopy_radiation_fraction = 0.8,
        },
    });
    try std.testing.expect(result.applied);
    try std.testing.expect(result.clear_leaf_sheath_and_sapwood_totals);
    try std.testing.expect(!result.terminate_living_branches);
    try std.testing.expectEqual(@as(f64, 2), result.state.living_population_per_m2);
    try std.testing.expectEqual(@as(f64, 10), result.state.living_population_count);
    try std.testing.expectEqual(@as(f64, 3), result.state.standing_dead_population_count);
    try std.testing.expectEqual(@as(f64, 0.2), result.state.canopy_radiation_fraction);
}

test "source-order tillage population selector preserves crop and date gates" {
    const base: plant_harvest_runtime.SourceOrderTillagePopulationInput = .{
        .hour_of_day = 11,
        .local_solar_noon_h = 11.4,
        .biomass_turnover_type = 1,
        .root_profile_type = 1,
        .current_day_of_year = 200,
        .current_year = 2000,
        .planting_day_of_year = 100,
        .planting_year = 2001,
        .tillage_code = 15,
        .is_first_plant_population = true,
        .remaining_fraction = 0,
        .zero_population_threshold = 1.0e-9,
        .state = .{
            .living_population_per_m2 = 1,
            .living_population_count = 1,
            .standing_dead_population_count = 1,
            .canopy_radiation_fraction = 1,
        },
    };
    try std.testing.expect(!(try plant_harvest_runtime.sourceOrderTillagePopulationReduction(base)).applied);

    var second_population = base;
    second_population.is_first_plant_population = false;
    const source_date_result = try plant_harvest_runtime.sourceOrderTillagePopulationReduction(second_population);
    try std.testing.expect(source_date_result.applied);
    try std.testing.expect(source_date_result.terminate_living_branches);

    var planting_date = second_population;
    planting_date.current_day_of_year = planting_date.planting_day_of_year;
    planting_date.current_year = planting_date.planting_year;
    try std.testing.expect(!(try plant_harvest_runtime.sourceOrderTillagePopulationReduction(planting_date)).applied);
}

test "source-order tillage branch litter conserves every element" {
    const kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 0.1, 0.2, 0.3, 0.4 },
        .nitrogen = .{ 0.4, 0.3, 0.2, 0.1 },
        .phosphorus = .{ 0.25, 0.25, 0.25, 0.25 },
    };
    const composition: plant_harvest_runtime.TillageElementComposition = .{
        .carbon = .{ 0.25, 0.75 },
        .nitrogen = .{ 0.4, 0.6 },
        .phosphorus = .{ 0.5, 0.5 },
    };
    const pools: plant_harvest_runtime.SourceOrderTillageBranchPools = .{
        .host_mobile = .{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 },
        .symbiont_mobile = .{ .carbon_g = 2, .nitrogen_g = 0.2, .phosphorus_g = 0.02 },
        .c4_mobile_carbon_g_c = 0.5,
        .stalk_reserve = .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 },
        .leaf = .{ .carbon_g = 4, .nitrogen_g = 0.4, .phosphorus_g = 0.04 },
        .symbiont_structural = .{ .carbon_g = 5, .nitrogen_g = 0.5, .phosphorus_g = 0.05 },
        .sheath = .{ .carbon_g = 6, .nitrogen_g = 0.6, .phosphorus_g = 0.06 },
        .husk = .{ .carbon_g = 7, .nitrogen_g = 0.7, .phosphorus_g = 0.07 },
        .ear = .{ .carbon_g = 8, .nitrogen_g = 0.8, .phosphorus_g = 0.08 },
        .grain = .{ .carbon_g = 9, .nitrogen_g = 0.9, .phosphorus_g = 0.09 },
        .stalk = .{ .carbon_g = 10, .nitrogen_g = 1, .phosphorus_g = 0.1 },
    };
    const result = try plant_harvest_runtime.sourceOrderTillageBranchLitter(.{
        .remaining_fraction = 0.4,
        .winter_annual = true,
        .pools = pools,
        .leaf_composition = composition,
        .sheath_composition = composition,
        .stalk_composition = composition,
        .nonstructural_kinetics = kinetics,
        .foliar_kinetics = kinetics,
        .nonfoliar_kinetics = kinetics,
        .stalk_kinetics = kinetics,
        .coarse_wood_kinetics = kinetics,
    });
    var litter_carbon_g_c: f64 = 0;
    var litter_nitrogen_g_n: f64 = 0;
    var litter_phosphorus_g_p: f64 = 0;
    for (0..litter_partition.kinetic_component_count) |kinetic| {
        litter_carbon_g_c += result.litter.woody_carbon_g[kinetic] + result.litter.nonwoody_carbon_g[kinetic];
        litter_nitrogen_g_n += result.litter.woody_nitrogen_g[kinetic] + result.litter.nonwoody_nitrogen_g[kinetic];
        litter_phosphorus_g_p += result.litter.woody_phosphorus_g[kinetic] + result.litter.nonwoody_phosphorus_g[kinetic];
    }
    try std.testing.expectApproxEqAbs(0.6 * 55.5, litter_carbon_g_c + result.seasonal_storage.carbon_g, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.6 * 5.5, litter_nitrogen_g_n + result.seasonal_storage.nitrogen_g, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.6 * 0.55, litter_phosphorus_g_p + result.seasonal_storage.phosphorus_g, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.6 * pools.grain.carbon_g, result.seasonal_storage.carbon_g, 1.0e-12);
}

test "source-order tillage routes non-winter grain to nonwoody litter" {
    const one_hot: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const zero_mass: canopy.ElementalMass = .{};
    const result = try plant_harvest_runtime.sourceOrderTillageBranchLitter(.{
        .remaining_fraction = 0.25,
        .winter_annual = false,
        .pools = .{
            .host_mobile = zero_mass,
            .symbiont_mobile = zero_mass,
            .c4_mobile_carbon_g_c = 0,
            .stalk_reserve = zero_mass,
            .leaf = zero_mass,
            .symbiont_structural = zero_mass,
            .sheath = zero_mass,
            .husk = zero_mass,
            .ear = zero_mass,
            .grain = .{ .carbon_g = 4, .nitrogen_g = 0.4, .phosphorus_g = 0.04 },
            .stalk = zero_mass,
        },
        .leaf_composition = .{ .carbon = .{ 0, 1 }, .nitrogen = .{ 0, 1 }, .phosphorus = .{ 0, 1 } },
        .sheath_composition = .{ .carbon = .{ 0, 1 }, .nitrogen = .{ 0, 1 }, .phosphorus = .{ 0, 1 } },
        .stalk_composition = .{ .carbon = .{ 0, 1 }, .nitrogen = .{ 0, 1 }, .phosphorus = .{ 0, 1 } },
        .nonstructural_kinetics = one_hot,
        .foliar_kinetics = one_hot,
        .nonfoliar_kinetics = one_hot,
        .stalk_kinetics = one_hot,
        .coarse_wood_kinetics = one_hot,
    });
    try std.testing.expectEqual(@as(f64, 3), result.litter.nonwoody_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), result.seasonal_storage.carbon_g);
}

test "source-order tillage retains complete branch state over runtime extents" {
    const mass: canopy.ElementalMass = .{ .carbon_g = 8, .nitrogen_g = 4, .phosphorus_g = 2 };
    var scalar: plant_harvest_runtime.SourceOrderTillageBranchScalarState = .{
        .host_mobile = mass,
        .c4_mobile_carbon_g_c = 8,
        .symbiont_mobile = mass,
        .total_shoot = mass,
        .leaf = mass,
        .symbiont_structural = mass,
        .sheath = mass,
        .stalk = mass,
        .sapwood_carbon_g_c = 8,
        .stalk_reserve = mass,
        .husk = mass,
        .ear = mass,
        .grain = mass,
        .potential_seed_site_count = 8,
        .seed_count = 8,
        .individual_seed_carbon_g_c = 3,
        .leaf_area_m2 = 8,
        .stalk_total = mass,
    };
    var node_values: [16][3]f64 = @splat(.{ 2, 4, 6 });
    var sample_values: [4][2]f64 = @splat(.{ 10, 12 });
    const result = try plant_harvest_runtime.sourceOrderRetainTillageBranchState(&scalar, .{
        .c3_mobile_carbon_g_c = &node_values[0],
        .c4_mobile_carbon_g_c = &node_values[1],
        .carbon_dioxide_g_c = &node_values[2],
        .bicarbonate_g_c = &node_values[3],
        .leaf_area_m2 = &node_values[4],
        .growing_leaf_carbon_g_c = &node_values[5],
        .senescing_leaf_carbon_g_c = &node_values[6],
        .growing_sheath_carbon_g_c = &node_values[7],
        .senescing_sheath_carbon_g_c = &node_values[8],
        .growing_node_carbon_g_c = &node_values[9],
        .growing_leaf_nitrogen_g_n = &node_values[10],
        .growing_sheath_nitrogen_g_n = &node_values[11],
        .growing_node_nitrogen_g_n = &node_values[12],
        .growing_leaf_phosphorus_g_p = &node_values[13],
        .growing_sheath_phosphorus_g_p = &node_values[14],
        .growing_node_phosphorus_g_p = &node_values[15],
    }, .{
        .leaf_area_m2 = &sample_values[0],
        .growing_leaf_carbon_g_c = &sample_values[1],
        .growing_leaf_nitrogen_g_n = &sample_values[2],
        .growing_leaf_phosphorus_g_p = &sample_values[3],
    }, 0.25);

    try std.testing.expectEqual(@as(f64, 2), scalar.host_mobile.carbon_g);
    try std.testing.expectEqual(@as(f64, 3), scalar.individual_seed_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 4), result.leaf_sheath_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 2), result.sapwood_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 2), node_values[0][0]);
    try std.testing.expectEqual(@as(f64, 1), node_values[0][1]);
    try std.testing.expectEqual(@as(f64, 0.5), node_values[4][0]);
    try std.testing.expectEqual(@as(f64, 2.5), sample_values[0][0]);
}

test "source-order tillage standing dead repartitions and conserves C N P" {
    const stalk_kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 0.1, 0.2, 0.3, 0.4 },
        .nitrogen = .{ 0.4, 0.3, 0.2, 0.1 },
        .phosphorus = .{ 0.25, 0.25, 0.25, 0.25 },
    };
    const coarse_kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 0.4, 0.3, 0.2, 0.1 },
        .nitrogen = .{ 0.1, 0.2, 0.3, 0.4 },
        .phosphorus = .{ 0.25, 0.25, 0.25, 0.25 },
    };
    const result = try plant_harvest_runtime.sourceOrderTillageStandingDead(.{
        .remaining_fraction = 0.25,
        .standing_dead_by_source_component = .{
            .{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 },
            .{ .carbon_g = 2, .nitrogen_g = 0.2, .phosphorus_g = 0.02 },
            .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 },
            .{ .carbon_g = 4, .nitrogen_g = 0.4, .phosphorus_g = 0.04 },
        },
        .composition = .{
            .carbon = .{ 0.6, 0.4 },
            .nitrogen = .{ 0.7, 0.3 },
            .phosphorus = .{ 0.8, 0.2 },
        },
        .stalk_kinetics = stalk_kinetics,
        .coarse_wood_kinetics = coarse_kinetics,
    });
    var litter_carbon_g_c: f64 = 0;
    var litter_nitrogen_g_n: f64 = 0;
    var litter_phosphorus_g_p: f64 = 0;
    for (0..litter_partition.kinetic_component_count) |kinetic| {
        litter_carbon_g_c += result.litter.woody_carbon_g[kinetic] + result.litter.nonwoody_carbon_g[kinetic];
        litter_nitrogen_g_n += result.litter.woody_nitrogen_g[kinetic] + result.litter.nonwoody_nitrogen_g[kinetic];
        litter_phosphorus_g_p += result.litter.woody_phosphorus_g[kinetic] + result.litter.nonwoody_phosphorus_g[kinetic];
    }
    try std.testing.expectApproxEqAbs(7.5, litter_carbon_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.75, litter_nitrogen_g_n, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.075, litter_phosphorus_g_p, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 1), result.remaining_by_source_component[3].carbon_g);
    try std.testing.expectApproxEqAbs(
        0.75 * 0.4 * 10 * 0.6,
        result.litter.woody_carbon_g[0],
        1.0e-12,
    );
}

test "source-order tillage termination sets all death flags and harvest date" {
    const prior: plant_harvest_runtime.SourceOrderTillageTerminationState = .{
        .roots_dead = false,
        .shoots_dead = false,
        .plant_dead = false,
        .harvest_termination_code = 0,
        .harvest_day_of_year = 0,
        .harvest_year = 0,
    };
    const alive = try plant_harvest_runtime.sourceOrderTillageTermination(.{
        .living_population_count = 1.1e-6,
        .zero_population_threshold = 1.0e-6,
        .current_day_of_year = 240,
        .current_year = 2004,
        .state = prior,
    });
    try std.testing.expect(!alive.terminated);
    try std.testing.expectEqualDeep(prior, alive.state);

    const terminated = try plant_harvest_runtime.sourceOrderTillageTermination(.{
        .living_population_count = 1.0e-6,
        .zero_population_threshold = 1.0e-6,
        .current_day_of_year = 240,
        .current_year = 2004,
        .state = prior,
    });
    try std.testing.expect(terminated.terminated);
    try std.testing.expect(terminated.state.roots_dead);
    try std.testing.expect(terminated.state.shoots_dead);
    try std.testing.expect(terminated.state.plant_dead);
    try std.testing.expectEqual(@as(u8, 1), terminated.state.harvest_termination_code);
    try std.testing.expectEqual(@as(u16, 240), terminated.state.harvest_day_of_year);
    try std.testing.expectEqual(@as(u32, 2004), terminated.state.harvest_year);
}

test "source-order standing dead harvest separates harvest litter and retention" {
    const result = try plant_harvest_runtime.sourceOrderStandingDeadHarvest(.{
        .harvest_code = 0,
        .hour_of_day = 12,
        .local_solar_noon_h = 12.8,
        .thinning_fraction_or_specific_consumption_rate = 0.5,
        .standing_dead_removal_fraction = 0.4,
        .grazer_live_mass_g_per_m2 = 0,
        .animal_accessible_area_m2 = 0,
        .insect_accessible_area_m2 = 0,
        .standing_dead_presence_threshold_g_c = 1.0e-12,
        .standing_dead_by_component = .{
            .{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 },
            .{ .carbon_g = 2, .nitrogen_g = 0.2, .phosphorus_g = 0.02 },
            .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 },
            .{ .carbon_g = 4, .nitrogen_g = 0.4, .phosphorus_g = 0.04 },
            .{ .carbon_g = 5, .nitrogen_g = 0.5, .phosphorus_g = 0.05 },
        },
    });
    try std.testing.expectEqual(@as(f64, 0.5), result.retained_fraction);
    try std.testing.expectEqual(@as(f64, 0.8), result.harvested_fraction);
    try std.testing.expectApproxEqAbs(3, result.harvested.carbon_g, 1.0e-14);
    try std.testing.expectApproxEqAbs(4.5, result.returned_to_litter.carbon_g, 1.0e-14);
    try std.testing.expectEqual(@as(f64, 2.5), result.remaining_by_component[4].carbon_g);
    try std.testing.expectApproxEqAbs(
        15,
        result.harvested.carbon_g + result.returned_to_litter.carbon_g +
            result.remaining_by_component[0].carbon_g + result.remaining_by_component[1].carbon_g +
            result.remaining_by_component[2].carbon_g + result.remaining_by_component[3].carbon_g +
            result.remaining_by_component[4].carbon_g,
        1.0e-14,
    );
}

test "source-order standing dead grazing selects animal and insect areas" {
    const base: plant_harvest_runtime.SourceOrderStandingDeadHarvestInput = .{
        .harvest_code = 4,
        .hour_of_day = 3,
        .local_solar_noon_h = 12,
        .thinning_fraction_or_specific_consumption_rate = 1,
        .standing_dead_removal_fraction = 0.5,
        .grazer_live_mass_g_per_m2 = 96,
        .animal_accessible_area_m2 = 2,
        .insect_accessible_area_m2 = 5,
        .standing_dead_presence_threshold_g_c = 1,
        .standing_dead_by_component = .{
            .{ .carbon_g = 2 },
            .{ .carbon_g = 2 },
            .{ .carbon_g = 2 },
            .{ .carbon_g = 2 },
            .{ .carbon_g = 2 },
        },
    };
    const animal = try plant_harvest_runtime.sourceOrderStandingDeadHarvest(base);
    try std.testing.expectEqual(@as(f64, 0.8), animal.retained_fraction);
    var insect_input = base;
    insect_input.harvest_code = 6;
    const insect = try plant_harvest_runtime.sourceOrderStandingDeadHarvest(insect_input);
    try std.testing.expectEqual(@as(f64, 0.5), insect.retained_fraction);
    var absent = base;
    absent.standing_dead_presence_threshold_g_c = 10;
    try std.testing.expectEqual(
        @as(f64, 1),
        (try plant_harvest_runtime.sourceOrderStandingDeadHarvest(absent)).retained_fraction,
    );
}

test "source-order harvest residue routes all five components" {
    const harvested = [plant_harvest_runtime.harvest_product_component_count]canopy.ElementalMass{
        .{ .carbon_g = 10, .nitrogen_g = 1, .phosphorus_g = 0.1 },
        .{ .carbon_g = 20, .nitrogen_g = 2, .phosphorus_g = 0.2 },
        .{ .carbon_g = 30, .nitrogen_g = 3, .phosphorus_g = 0.3 },
        .{ .carbon_g = 40, .nitrogen_g = 4, .phosphorus_g = 0.4 },
        .{ .carbon_g = 50, .nitrogen_g = 5, .phosphorus_g = 0.5 },
    };
    const residue = try plant_harvest_runtime.sourceOrderHarvestResidueRouting(.{
        .harvest_code = 0,
        .harvested_by_component = harvested,
        .harvested_grain = .{},
        .ecosystem_export_fraction = .{ 0.1, 0.2, 0.3, 0.4 },
    });
    try std.testing.expectEqual(@as(f64, 9), residue[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 18), residue[1].carbon_g);
    try std.testing.expectEqual(@as(f64, 24), residue[2].carbon_g);
    try std.testing.expectEqual(@as(f64, 28), residue[3].carbon_g);
    try std.testing.expectEqual(@as(f64, 30), residue[4].carbon_g);

    const grazing = try plant_harvest_runtime.sourceOrderHarvestResidueRouting(.{
        .harvest_code = 6,
        .harvested_by_component = harvested,
        .harvested_grain = .{},
        .ecosystem_export_fraction = .{ 0.1, 0.2, 0.3, 0.4 },
    });
    try std.testing.expectEqualDeep(residue, grazing);
}

test "source-order grain harvest subtracts only exported grain from component two" {
    const harvested = [plant_harvest_runtime.harvest_product_component_count]canopy.ElementalMass{
        .{ .carbon_g = 10, .nitrogen_g = 1, .phosphorus_g = 0.1 },
        .{ .carbon_g = 20, .nitrogen_g = 2, .phosphorus_g = 0.2 },
        .{ .carbon_g = 30, .nitrogen_g = 3, .phosphorus_g = 0.3 },
        .{ .carbon_g = 40, .nitrogen_g = 4, .phosphorus_g = 0.4 },
        .{ .carbon_g = 50, .nitrogen_g = 5, .phosphorus_g = 0.5 },
    };
    const residue = try plant_harvest_runtime.sourceOrderHarvestResidueRouting(.{
        .harvest_code = 1,
        .harvested_by_component = harvested,
        .harvested_grain = .{ .carbon_g = 10, .nitrogen_g = 1, .phosphorus_g = 0.1 },
        .ecosystem_export_fraction = .{ 0.9, 0.5, 0.8, 0.7 },
    });
    try std.testing.expectEqualDeep(harvested[0], residue[0]);
    try std.testing.expectEqualDeep(harvested[1], residue[1]);
    try std.testing.expectEqual(@as(f64, 25), residue[2].carbon_g);
    try std.testing.expectEqualDeep(harvested[3], residue[3]);
    try std.testing.expectEqualDeep(harvested[4], residue[4]);
}

test "source-order disturbance totals route ordinary export and reseed storage" {
    const harvested: [plant_harvest_runtime.harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 2, .nitrogen_g = 1, .phosphorus_g = 0.5 });
    const residue: [plant_harvest_runtime.harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 1, .nitrogen_g = 0.5, .phosphorus_g = 0.25 });
    const direct_litter: [plant_harvest_runtime.harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 0.2, .nitrogen_g = 0.1, .phosphorus_g = 0.05 });
    const ordinary = try plant_harvest_runtime.sourceOrderTotalDisturbanceRemoval(.{
        .harvest_code = 2,
        .terminate_and_reseed = false,
        .grazer_growth_yield = 0,
        .grazer_respiration_fraction = 0,
        .harvested_by_component = harvested,
        .residue_by_component = residue,
        .direct_litter_by_component = direct_litter,
    });
    try std.testing.expectEqual(@as(f64, 10), ordinary.harvested_total.carbon_g);
    try std.testing.expectEqual(@as(f64, 5), ordinary.residue_total.carbon_g);
    try std.testing.expectEqual(@as(f64, 1), ordinary.direct_litter_total.carbon_g);
    try std.testing.expectEqual(@as(f64, 5), ordinary.plant_ecosystem_removal.carbon_g);
    try std.testing.expectEqual(@as(f64, -5), ordinary.net_biome_production_carbon_change_g_c_per_h);

    const reseed = try plant_harvest_runtime.sourceOrderTotalDisturbanceRemoval(.{
        .harvest_code = 2,
        .terminate_and_reseed = true,
        .grazer_growth_yield = 0,
        .grazer_respiration_fraction = 0,
        .harvested_by_component = harvested,
        .residue_by_component = residue,
        .direct_litter_by_component = direct_litter,
    });
    try std.testing.expectEqual(@as(f64, 5), reseed.reseed_storage_addition.carbon_g);
    try std.testing.expectEqual(@as(f64, 0), reseed.plant_ecosystem_removal.carbon_g);
    try std.testing.expectEqual(@as(f64, 0), reseed.grid_ecosystem_removal.carbon_g);
}

test "source-order grazing totals split growth and respiration carbon" {
    const harvested: [plant_harvest_runtime.harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 2, .nitrogen_g = 1, .phosphorus_g = 0.5 });
    const residue: [plant_harvest_runtime.harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 1, .nitrogen_g = 0.5, .phosphorus_g = 0.25 });
    const result = try plant_harvest_runtime.sourceOrderTotalDisturbanceRemoval(.{
        .harvest_code = 4,
        .terminate_and_reseed = false,
        .grazer_growth_yield = 0.4,
        .grazer_respiration_fraction = 0.6,
        .harvested_by_component = harvested,
        .residue_by_component = residue,
        .direct_litter_by_component = @splat(.{}),
    });
    try std.testing.expectEqual(@as(f64, 2), result.plant_ecosystem_removal.carbon_g);
    try std.testing.expectEqual(@as(f64, 2.5), result.plant_ecosystem_removal.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 1.25), result.plant_ecosystem_removal.phosphorus_g);
    try std.testing.expectEqual(@as(f64, -3), result.plant_total_respiration_change_g_c_per_h);
    try std.testing.expectEqual(@as(f64, -3), result.plant_actual_respiration_change_g_c_per_h);
    try std.testing.expectEqual(@as(f64, -3), result.ecosystem_respiration_change_g_c_per_h);
    try std.testing.expectEqual(@as(f64, -3), result.autotrophic_respiration_change_g_c_per_h);
}

test "source-order aboveground harvest litter preserves herbaceous and woody routing" {
    const residue = [plant_harvest_runtime.harvest_product_component_count]canopy.ElementalMass{
        .{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 },
        .{ .carbon_g = 2, .nitrogen_g = 0.2, .phosphorus_g = 0.02 },
        .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 },
        .{ .carbon_g = 4, .nitrogen_g = 0.4, .phosphorus_g = 0.04 },
        .{ .carbon_g = 5, .nitrogen_g = 0.5, .phosphorus_g = 0.05 },
    };
    const direct: [plant_harvest_runtime.harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 });
    const one_hot: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const common: plant_harvest_runtime.SourceOrderAbovegroundHarvestLitterInput = .{
        .harvest_code = 2,
        .biomass_turnover_type = 0,
        .root_profile_type = 2,
        .residue_by_component = residue,
        .direct_litter_by_component = direct,
        .woody_composition = .{
            .carbon = .{ 0.25, 0.75 },
            .nitrogen = .{ 0.25, 0.75 },
            .phosphorus = .{ 0.25, 0.75 },
        },
        .nonstructural_kinetics = one_hot,
        .foliar_kinetics = one_hot,
        .nonfoliar_kinetics = one_hot,
        .stalk_kinetics = one_hot,
        .coarse_wood_kinetics = one_hot,
    };
    const herbaceous = try plant_harvest_runtime.sourceOrderAbovegroundHarvestLitter(common);
    try std.testing.expectEqual(@as(f64, 20), herbaceous.litter.nonwoody_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), herbaceous.litter.woody_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), herbaceous.standing_dead_addition.woody_carbon_g[0]);

    var woody_input = common;
    woody_input.biomass_turnover_type = 2;
    const woody = try plant_harvest_runtime.sourceOrderAbovegroundHarvestLitter(woody_input);
    try std.testing.expectEqual(@as(f64, 15.75), woody.litter.nonwoody_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 2.25), woody.litter.woody_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 2), woody.standing_dead_addition.woody_carbon_g[0]);
}

test "source-order grazing ledgers preserve five components and manure partition" {
    const residue: [plant_harvest_runtime.harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 1, .nitrogen_g = 0.2, .phosphorus_g = 0.04 });
    const direct = residue;
    const current: plant_harvest_runtime.SourceOrderGrazingLitterLedgerState = .{
        .hourly_litter = .{},
        .cumulative_litter = .{ .carbon_g = 20, .nitrogen_g = 10, .phosphorus_g = 1 },
        .cumulative_aboveground_litter = .{ .carbon_g = 4, .nitrogen_g = 3, .phosphorus_g = 0.3 },
        .surface_litter_carbon_g_c = 7,
        .accumulated_application = .{},
    };
    const animal = try plant_harvest_runtime.sourceOrderGrazingLitterLedgers(4, residue, direct, current);
    try std.testing.expectEqual(@as(f64, 10), animal.returned_mass.carbon_g);
    try std.testing.expectEqual(@as(f64, 2), animal.returned_mass.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 14), animal.state.cumulative_aboveground_litter.carbon_g);
    try std.testing.expectEqual(@as(f64, 14), animal.state.cumulative_aboveground_litter.nitrogen_g);
    try std.testing.expectApproxEqAbs(1.8, animal.state.cumulative_aboveground_litter.phosphorus_g, 1.0e-14);
    try std.testing.expectEqual(@as(f64, 17), animal.state.surface_litter_carbon_g_c);
    try std.testing.expectApproxEqAbs(0.36, animal.manure.organic_by_biochemical_fraction[0].carbon_g, 1.0e-14);
    try std.testing.expectEqual(@as(f64, 1), animal.manure.inorganic_nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 0.2), animal.manure.inorganic_phosphorus_g_p);

    const insect = try plant_harvest_runtime.sourceOrderGrazingLitterLedgers(6, residue, direct, current);
    try std.testing.expectApproxEqAbs(1.38, insect.manure.organic_by_biochemical_fraction[0].carbon_g, 1.0e-14);
}

test "source-order population thresholds preserve runtime scaling and evaluation order" {
    const result = try plant_harvest_runtime.sourceOrderPopulationScaledNumericalThresholds(
        250,
        50,
        1.0e-15,
        1.0e-6,
    );
    try std.testing.expectEqual(@as(f64, 2.5e-13), result.plant_mass_presence_g);
    try std.testing.expectEqual(@as(f64, 5.0e-15), result.plant_mass_density_g_m2);
    try std.testing.expectEqual(@as(f64, 2.5e-4), result.plant_flux_presence_g_per_step);

    const extinct = try plant_harvest_runtime.sourceOrderPopulationScaledNumericalThresholds(0, 50, 1.0e-15, 1.0e-6);
    try std.testing.expectEqual(@as(f64, 0), extinct.plant_mass_presence_g);
    try std.testing.expectEqual(@as(f64, 0), extinct.plant_mass_density_g_m2);
    try std.testing.expectEqual(@as(f64, 0), extinct.plant_flux_presence_g_per_step);

    try std.testing.expectError(
        error.InvalidPopulationScaledNumericalThresholdInput,
        plant_harvest_runtime.sourceOrderPopulationScaledNumericalThresholds(1, 0, 1.0e-15, 1.0e-6),
    );
    try std.testing.expectError(
        error.InvalidPopulationScaledNumericalThresholdInput,
        plant_harvest_runtime.sourceOrderPopulationScaledNumericalThresholds(1, 1, -1, 1.0e-6),
    );
}

test "source-order dead branch reset preserves selector and runtime branch loop" {
    const stale: plant_harvest_runtime.SourceOrderDeadBranchPhenologyState = .{
        .dead = true,
        .maturity_group_node_count = 9,
        .initiated_node_count = 8,
        .nodes_at_floral_initiation = 7,
        .nodes_at_anthesis = 6,
        .appeared_leaf_count = 5,
        .leaves_at_floral_initiation = 4,
        .current_leaf_ordinal = 3,
        .current_growing_leaf_ordinal = 2,
        .normalized_vegetative_node_change = 1,
        .normalized_reproductive_node_change = 1,
        .accumulated_leafout_h = 1,
        .accumulated_leafoff_h = 1,
        .lengthening_photoperiod_h = 1,
        .shortening_photoperiod_h = 1,
        .time_since_germination_h = 1,
        .hours_without_grain_fill = 1,
        .carbon_fixation_feedback = 0.5,
        .carbon_fixation_feedback_previous = 0.5,
        .leafout_initialization_enabled = false,
        .emergence_initialization_disabled = false,
        .leafoff_enabled = false,
        .remobilization_enabled = false,
        .hours_after_maturity_h = 1,
        .new_branch_count = 3,
        .stage_day_of_year = @splat(100),
    };
    var branches = [_]plant_harvest_runtime.SourceOrderDeadBranchPhenologyState{ stale, stale };
    branches[1].dead = false;
    const live_before = branches[1];
    const applied = try plant_harvest_runtime.sourceOrderResetDeadBranchPhenology(&branches, .{
        .first_living_branch_emergence_day_of_year = 20,
        .perennial_growth_habit = false,
        .current_day_of_year = 200,
        .current_year = 2025,
        .harvest_day_of_year = 200,
        .harvest_year = 2025,
        .hour_of_day = 12,
        .local_solar_noon_h = 12.8,
        .initial_maturity_group_node_count = 4,
        .initial_node_count = 0.25,
    });
    try std.testing.expect(applied);
    try std.testing.expectEqual(@as(f64, 4), branches[0].maturity_group_node_count);
    try std.testing.expectEqual(@as(f64, 0.25), branches[0].initiated_node_count);
    try std.testing.expectEqual(@as(f64, 0.25), branches[0].nodes_at_floral_initiation);
    try std.testing.expectEqual(@as(usize, 1), branches[0].current_leaf_ordinal);
    try std.testing.expectEqual(@as(f64, 1), branches[0].carbon_fixation_feedback);
    try std.testing.expect(branches[0].leafout_initialization_enabled);
    try std.testing.expect(branches[0].emergence_initialization_disabled);
    try std.testing.expectEqual([_]u16{0} ** 10, branches[0].stage_day_of_year);
    try std.testing.expectEqualDeep(live_before, branches[1]);

    var unselected = [_]plant_harvest_runtime.SourceOrderDeadBranchPhenologyState{stale};
    const skipped = try plant_harvest_runtime.sourceOrderResetDeadBranchPhenology(&unselected, .{
        .first_living_branch_emergence_day_of_year = 20,
        .perennial_growth_habit = false,
        .current_day_of_year = 199,
        .current_year = 2025,
        .harvest_day_of_year = 200,
        .harvest_year = 2025,
        .hour_of_day = 12,
        .local_solar_noon_h = 12,
        .initial_maturity_group_node_count = 4,
        .initial_node_count = 0.25,
    });
    try std.testing.expect(!skipped);
    try std.testing.expectEqualDeep(stale, unselected[0]);
}

test "source-order dead branch litterfall conserves runtime branch components" {
    const kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const composition: plant_harvest_runtime.TillageElementComposition = .{
        .carbon = .{ 0.25, 0.75 },
        .nitrogen = .{ 0.25, 0.75 },
        .phosphorus = .{ 0.25, 0.75 },
    };
    const common: plant_harvest_runtime.SourceOrderDeadBranchLitterInput = .{
        .annual_growth_habit = false,
        .deciduous_phenology = false,
        .pools = .{
            .bacterial_nonstructural = .{ .carbon_g = 1, .nitrogen_g = 1, .phosphorus_g = 1 },
            .bacterial_structural = .{ .carbon_g = 2, .nitrogen_g = 2, .phosphorus_g = 2 },
            .leaf = .{ .carbon_g = 4, .nitrogen_g = 4, .phosphorus_g = 4 },
            .sheath = .{ .carbon_g = 8, .nitrogen_g = 8, .phosphorus_g = 8 },
            .husk = .{ .carbon_g = 16, .nitrogen_g = 16, .phosphorus_g = 16 },
            .ear = .{ .carbon_g = 32, .nitrogen_g = 32, .phosphorus_g = 32 },
            .grain = .{ .carbon_g = 64, .nitrogen_g = 64, .phosphorus_g = 64 },
            .stalk = .{ .carbon_g = 128, .nitrogen_g = 128, .phosphorus_g = 128 },
        },
        .leaf_woody_fraction = composition,
        .sheath_woody_fraction = composition,
        .nonstructural_kinetics = kinetics,
        .foliar_kinetics = kinetics,
        .nonfoliar_kinetics = kinetics,
        .stalk_kinetics = kinetics,
        .coarse_wood_kinetics = kinetics,
    };
    const ordinary = try plant_harvest_runtime.sourceOrderDeadBranchLitterfall(common);
    try std.testing.expectEqual(@as(f64, 124), ordinary.nonwoody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 3), ordinary.woody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 128), ordinary.standing_dead_stalk[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 0), ordinary.seasonal_storage_addition.carbon_g);

    var winter_input = common;
    winter_input.annual_growth_habit = true;
    winter_input.deciduous_phenology = true;
    const winter = try plant_harvest_runtime.sourceOrderDeadBranchLitterfall(winter_input);
    try std.testing.expectEqual(@as(f64, 60), winter.nonwoody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 64), winter.seasonal_storage_addition.carbon_g);
    try std.testing.expectEqual(@as(f64, 255), winter.nonwoody_litter[0].carbon_g +
        winter.woody_litter[0].carbon_g +
        winter.standing_dead_stalk[0].carbon_g +
        winter.seasonal_storage_addition.carbon_g);
}

test "source-order dead branch storage recovery preserves six assignments" {
    const result = try plant_harvest_runtime.sourceOrderDeadBranchStorageRecovery(.{
        .current_seasonal_storage = .{
            .carbon_g = 100,
            .nitrogen_g = 20,
            .phosphorus_g = 5,
        },
        .branch_mobile = .{
            .carbon_g = 7,
            .nitrogen_g = 3,
            .phosphorus_g = 0.4,
        },
        .c4_intermediate_carbon_g_c = 11,
        .stalk_reserve = .{
            .carbon_g = 13,
            .nitrogen_g = 2,
            .phosphorus_g = 0.6,
        },
    });
    try std.testing.expectEqual(@as(f64, 131), result.seasonal_storage.carbon_g);
    try std.testing.expectEqual(@as(f64, 25), result.seasonal_storage.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 6), result.seasonal_storage.phosphorus_g);
    try std.testing.expectEqual(@as(f64, 31), result.recovered.carbon_g);
    try std.testing.expectEqual(@as(f64, 5), result.recovered.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 1), result.recovered.phosphorus_g);

    try std.testing.expectError(
        error.InvalidDeadBranchStorageRecoveryInput,
        plant_harvest_runtime.sourceOrderDeadBranchStorageRecovery(.{
            .current_seasonal_storage = .{},
            .branch_mobile = .{},
            .c4_intermediate_carbon_g_c = -1,
            .stalk_reserve = .{},
        }),
    );
}

test "source-order dead branch canopy reset preserves node zero exceptions" {
    var scalar: plant_harvest_runtime.SourceOrderDeadBranchScalarResetState = undefined;
    inline for (@typeInfo(plant_harvest_runtime.SourceOrderDeadBranchScalarResetState).@"struct".fields) |field| {
        if (field.type == f64) {
            @field(scalar, field.name) = 1;
        } else {
            @field(scalar, field.name) = .{
                .carbon_g = 1,
                .nitrogen_g = 1,
                .phosphorus_g = 1,
            };
        }
    }
    const node: plant_harvest_runtime.SourceOrderDeadBranchNodeResetState = .{
        .bundle_sheath_mobile_carbon_g_c = 1,
        .mesophyll_mobile_carbon_g_c = 1,
        .bundle_sheath_co2_carbon_g_c = 1,
        .bundle_sheath_bicarbonate_carbon_g_c = 1,
        .leaf_area_m2 = 1,
        .node_height_m = 1,
        .node_height_previous_m = 1,
        .sheath_height_m = 1,
        .leaf = .{ .carbon_g = 1, .nitrogen_g = 1, .phosphorus_g = 1 },
        .leaf_protein_g = 1,
        .sheath = .{ .carbon_g = 1, .nitrogen_g = 1, .phosphorus_g = 1 },
        .sheath_protein_g = 1,
        .stalk = .{ .carbon_g = 1, .nitrogen_g = 1, .phosphorus_g = 1 },
    };
    var nodes = [_]plant_harvest_runtime.SourceOrderDeadBranchNodeResetState{ node, node };
    var node_layers = [_]plant_harvest_runtime.SourceOrderDeadBranchLayerResetState{
        .{
            .leaf_area_m2 = 2,
            .leaf = .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 },
            .projected_leaf_surface_m2 = @splat(4),
        },
        .{
            .leaf_area_m2 = 5,
            .leaf = .{ .carbon_g = 7, .nitrogen_g = 0.7, .phosphorus_g = 0.07 },
            .projected_leaf_surface_m2 = @splat(8),
        },
    };
    var canopy_area = [_]f64{10};
    var canopy_carbon = [_]f64{20};
    var stalk_area = [_]f64{6};
    var stalk_surface = [_][4]f64{@splat(9)};
    try plant_harvest_runtime.sourceOrderResetDeadBranchCanopy(.{
        .scalar = &scalar,
        .nodes = &nodes,
        .node_layers = &node_layers,
        .canopy_leaf_area_m2_by_layer = &canopy_area,
        .canopy_leaf_carbon_g_c_by_layer = &canopy_carbon,
        .branch_stalk_area_m2_by_layer = &stalk_area,
        .branch_projected_stalk_surface_m2 = &stalk_surface,
    });
    try std.testing.expectEqual(@as(f64, 3), canopy_area[0]);
    try std.testing.expectEqual(@as(f64, 10), canopy_carbon[0]);
    try std.testing.expectEqual(@as(f64, 0), scalar.host_mobile.carbon_g);
    try std.testing.expectEqual(@as(f64, 1), nodes[0].bundle_sheath_mobile_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), nodes[1].bundle_sheath_mobile_carbon_g_c);
    try std.testing.expectEqual([_]f64{4} ** 4, node_layers[0].projected_leaf_surface_m2);
    try std.testing.expectEqual([_]f64{0} ** 4, node_layers[1].projected_leaf_surface_m2);
    try std.testing.expectEqual(@as(f64, 0), stalk_area[0]);
    try std.testing.expectEqual([_]f64{0} ** 4, stalk_surface[0]);
}

test "source-order whole plant termination preserves winter annual reseed branch" {
    const current: plant_harvest_runtime.SourceOrderWholePlantTerminationState = .{
        .shoot_alive = true,
        .root_alive = true,
        .total_node_count = 12,
        .hours_below_leaf_turgor_threshold_h = 8,
        .main_stalk_diameter_m = 0.03,
        .branch_count = 3,
        .living_population_per_m2 = 20,
        .living_population_count = 200,
        .hypocotyl_height_m = 0.1,
    };
    const partial = try plant_harvest_runtime.sourceOrderWholePlantTermination(&.{ true, false, true }, false, current);
    try std.testing.expect(!partial.all_branches_dead);
    try std.testing.expectEqual(@as(usize, 2), partial.dead_branch_count);
    try std.testing.expectEqualDeep(current, partial.state);

    const winter = try plant_harvest_runtime.sourceOrderWholePlantTermination(&.{ true, true, true }, true, current);
    try std.testing.expect(winter.all_branches_dead);
    try std.testing.expectEqual(@as(usize, 1), winter.state.branch_count);
    try std.testing.expectEqual(@as(f64, 20), winter.state.living_population_per_m2);
    try std.testing.expectEqual(@as(f64, 200), winter.state.living_population_count);
    try std.testing.expect(!winter.state.shoot_alive);
    try std.testing.expect(!winter.state.root_alive);
    try std.testing.expectEqual(@as(f64, 0), winter.state.hypocotyl_height_m);

    const ordinary = try plant_harvest_runtime.sourceOrderWholePlantTermination(&.{ true, true, true }, false, current);
    try std.testing.expectEqual(@as(usize, 0), ordinary.state.branch_count);
    try std.testing.expectEqual(@as(f64, 0), ordinary.state.living_population_per_m2);
    try std.testing.expectEqual(@as(f64, 0), ordinary.state.living_population_count);
}

test "source-order dead root litterfall conserves runtime domains and axes" {
    const kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const mobile = [_]canopy.ElementalMass{
        .{ .carbon_g = 1, .nitrogen_g = 1, .phosphorus_g = 1 },
        .{ .carbon_g = 2, .nitrogen_g = 2, .phosphorus_g = 2 },
    };
    const axis: plant_harvest_runtime.SourceOrderDeadRootAxisPools = .{
        .primary = .{ .carbon_g = 1, .nitrogen_g = 1, .phosphorus_g = 1 },
        .secondary = .{ .carbon_g = 2, .nitrogen_g = 2, .phosphorus_g = 2 },
    };
    const structural = [_]plant_harvest_runtime.SourceOrderDeadRootAxisPools{ axis, axis, axis, axis };
    const base: plant_harvest_runtime.SourceOrderDeadRootLitterInput = .{
        .roots_dead = true,
        .root_domain_count = 2,
        .soil_layer_count = 1,
        .root_axis_count = 2,
        .mobile_by_domain_layer = &mobile,
        .structural_by_domain_layer_axis = &structural,
        .root_woody_fraction = .{
            .carbon = .{ 0.25, 0.75 },
            .nitrogen = .{ 0.25, 0.75 },
            .phosphorus = .{ 0.25, 0.75 },
        },
        .mobile_kinetics = kinetics,
        .fine_root_kinetics = kinetics,
        .coarse_root_kinetics = kinetics,
    };
    const dead = try plant_harvest_runtime.sourceOrderDeadRootLitterfall(std.testing.allocator, base);
    defer std.testing.allocator.free(dead);
    try std.testing.expectEqual(@as(f64, 3), dead[0].woody[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 12), dead[0].nonwoody[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 15), dead[0].woody[0].nitrogen_g +
        dead[0].nonwoody[0].nitrogen_g);

    var live_input = base;
    live_input.roots_dead = false;
    const live = try plant_harvest_runtime.sourceOrderDeadRootLitterfall(std.testing.allocator, live_input);
    defer std.testing.allocator.free(live);
    try std.testing.expectEqual(@as(f64, 0), live[0].woody[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 0), live[0].nonwoody[0].carbon_g);
}

test "source-order dead root gas release clears both phases conservatively" {
    const gaseous: plant_harvest_runtime.SourceOrderRootGasInventory = .{
        .carbon_dioxide_carbon_g_c = 1,
        .oxygen_g_o = 2,
        .methane_carbon_g_c = 3,
        .nitrous_oxide_nitrogen_g_n = 4,
        .ammonia_nitrogen_g_n = 5,
        .hydrogen_g_h = 6,
    };
    const aqueous: plant_harvest_runtime.SourceOrderRootGasInventory = .{
        .carbon_dioxide_carbon_g_c = 10,
        .oxygen_g_o = 20,
        .methane_carbon_g_c = 30,
        .nitrous_oxide_nitrogen_g_n = 40,
        .ammonia_nitrogen_g_n = 50,
        .hydrogen_g_h = 60,
    };
    var phases = [_]plant_harvest_runtime.SourceOrderRootGasPhases{
        .{ .gaseous = gaseous, .aqueous = aqueous },
        .{ .gaseous = gaseous, .aqueous = aqueous },
    };
    const loss = try plant_harvest_runtime.sourceOrderReleaseDeadRootGases(
        true,
        &phases,
        std.mem.zeroes(plant_harvest_runtime.SourceOrderRootGasInventory),
    );
    try std.testing.expectEqual(@as(f64, -22), loss.carbon_dioxide_carbon_g_c);
    try std.testing.expectEqual(@as(f64, -44), loss.oxygen_g_o);
    try std.testing.expectEqual(@as(f64, -132), loss.hydrogen_g_h);
    try std.testing.expectEqualDeep(
        std.mem.zeroes(plant_harvest_runtime.SourceOrderRootGasPhases),
        phases[0],
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(plant_harvest_runtime.SourceOrderRootGasPhases),
        phases[1],
    );

    var living = [_]plant_harvest_runtime.SourceOrderRootGasPhases{.{ .gaseous = gaseous, .aqueous = aqueous }};
    const unchanged_loss = try plant_harvest_runtime.sourceOrderReleaseDeadRootGases(
        false,
        &living,
        gaseous,
    );
    try std.testing.expectEqualDeep(gaseous, unchanged_loss);
    try std.testing.expectEqualDeep(gaseous, living[0].gaseous);
    try std.testing.expectEqualDeep(aqueous, living[0].aqueous);
}

test "source-order dead root reset preserves runtime extents and radius defaults" {
    const mass: canopy.ElementalMass = .{
        .carbon_g = 1,
        .nitrogen_g = 1,
        .phosphorus_g = 1,
    };
    const axis_layer_value: plant_harvest_runtime.SourceOrderDeadRootAxisLayerState = .{
        .primary = mass,
        .secondary = mass,
        .primary_length_m = 1,
        .secondary_length_m = 1,
        .secondary_axis_count = 1,
    };
    var axis_layers = [_]plant_harvest_runtime.SourceOrderDeadRootAxisLayerState{
        axis_layer_value,
        axis_layer_value,
        axis_layer_value,
        axis_layer_value,
    };
    var domain_axes = [_]plant_harvest_runtime.SourceOrderDeadRootDomainAxisState{
        .{ .primary_total = mass },
        .{ .primary_total = mass },
    };
    const layer_value: plant_harvest_runtime.SourceOrderDeadRootDomainLayerState = .{
        .mobile = mass,
        .active_root_carbon_g_c = 1,
        .actual_root_carbon_g_c = 1,
        .root_protein_g = 1,
        .primary_axis_count = 1,
        .total_axis_count = 1,
        .root_length_per_plant_m = 1,
        .root_length_density_m_m3 = 1,
        .gaseous_volume_m3 = 1,
        .aqueous_volume_m3 = 1,
        .root_surface_area_per_plant_m2 = 1,
        .primary_radius_m = 1,
        .secondary_radius_m = 1,
        .average_secondary_root_length_m = 1,
    };
    var domain_layers = [_]plant_harvest_runtime.SourceOrderDeadRootDomainLayerState{
        layer_value,
        layer_value,
    };
    try plant_harvest_runtime.sourceOrderResetDeadRootState(true, .{
        .root_domain_count = 1,
        .soil_layer_count = 2,
        .root_axis_count = 2,
        .axis_layer = &axis_layers,
        .domain_axis = &domain_axes,
        .domain_layer = &domain_layers,
        .initial_primary_radius_m_by_domain = &.{0.01},
        .initial_secondary_radius_m_by_domain = &.{0.002},
        .initial_average_secondary_root_length_m = 0.1,
    });
    for (axis_layers) |axis| {
        try std.testing.expectEqual(@as(f64, 0), axis.primary.carbon_g);
        try std.testing.expectEqual(@as(f64, 0), axis.secondary_length_m);
        try std.testing.expectEqual(@as(f64, 0), axis.secondary_axis_count);
    }
    for (domain_axes) |axis|
        try std.testing.expectEqual(@as(f64, 0), axis.primary_total.carbon_g);
    for (domain_layers) |layer| {
        try std.testing.expectEqual(@as(f64, 0), layer.mobile.carbon_g);
        try std.testing.expectEqual(@as(f64, 0.01), layer.primary_radius_m);
        try std.testing.expectEqual(@as(f64, 0.002), layer.secondary_radius_m);
        try std.testing.expectEqual(@as(f64, 0.1), layer.average_secondary_root_length_m);
        try std.testing.expectEqual(@as(f64, 0), layer.root_surface_area_per_plant_m2);
    }
}

test "source-order dead nodule litterfall uses only first root domain" {
    const kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const pools: plant_harvest_runtime.SourceOrderDeadNoduleLayerPools = .{
        .structural = .{ .carbon_g = 10, .nitrogen_g = 5, .phosphorus_g = 1 },
        .mobile = .{ .carbon_g = 2, .nitrogen_g = 1, .phosphorus_g = 0.2 },
    };
    var layers = [_]plant_harvest_runtime.SourceOrderDeadNoduleLayerPools{ pools, pools };
    const litter = try plant_harvest_runtime.sourceOrderDeadNoduleLitterfall(std.testing.allocator, .{
        .roots_dead = true,
        .nitrogen_fixation_enabled = true,
        .root_domain_count = 3,
        .layer_pools = &layers,
        .structural_kinetics = kinetics,
        .mobile_kinetics = kinetics,
    });
    defer std.testing.allocator.free(litter);
    try std.testing.expectEqual(@as(f64, 12), litter[0][0].carbon_g);
    try std.testing.expectEqual(@as(f64, 6), litter[1][0].nitrogen_g);
    try std.testing.expectEqual(@as(f64, 0), layers[0].structural.carbon_g);
    try std.testing.expectEqual(@as(f64, 0), layers[1].mobile.carbon_g);

    var nonfixing = [_]plant_harvest_runtime.SourceOrderDeadNoduleLayerPools{pools};
    const empty = try plant_harvest_runtime.sourceOrderDeadNoduleLitterfall(std.testing.allocator, .{
        .roots_dead = true,
        .nitrogen_fixation_enabled = false,
        .root_domain_count = 2,
        .layer_pools = &nonfixing,
        .structural_kinetics = kinetics,
        .mobile_kinetics = kinetics,
    });
    defer std.testing.allocator.free(empty);
    try std.testing.expectEqual(@as(f64, 0), empty[0][0].carbon_g);
    try std.testing.expectEqualDeep(pools, nonfixing[0]);
}

test "source-order dead root depth reset preserves axis-major domain order" {
    var deepest_by_axis = [_]usize{ 8, 9 };
    var depths = [_]f64{ 1, 2, 3, 4, 5, 6 };
    const mass: canopy.ElementalMass = .{
        .carbon_g = 1,
        .nitrogen_g = 0.1,
        .phosphorus_g = 0.01,
    };
    var totals = [_]canopy.ElementalMass{mass} ** 6;
    var deepest_active: usize = 9;
    var active_axes: usize = 2;
    try plant_harvest_runtime.sourceOrderResetDeadRootDepth(true, 3, 0.04, .{
        .root_domain_count = 3,
        .root_axis_count = 2,
        .deepest_layer_by_axis = &deepest_by_axis,
        .primary_depth_from_surface_m_by_axis_domain = &depths,
        .primary_total_by_axis_domain = &totals,
        .deepest_active_root_layer = &deepest_active,
        .active_root_axis_count = &active_axes,
    });
    try std.testing.expectEqual([_]usize{ 3, 3 }, deepest_by_axis);
    try std.testing.expectEqual([_]f64{0.04} ** 6, depths);
    for (totals) |total|
        try std.testing.expectEqual(@as(f64, 0), total.carbon_g);
    try std.testing.expectEqual(@as(usize, 3), deepest_active);
    try std.testing.expectEqual(@as(usize, 0), active_axes);
}

test "source-order complete death shoot litterfall conserves storage and branches" {
    const kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const composition: plant_harvest_runtime.TillageElementComposition = .{
        .carbon = .{ 0.25, 0.75 },
        .nitrogen = .{ 0.25, 0.75 },
        .phosphorus = .{ 0.25, 0.75 },
    };
    const carbon = struct {
        fn mass(value: f64) canopy.ElementalMass {
            return .{ .carbon_g = value, .nitrogen_g = 0, .phosphorus_g = 0 };
        }
    }.mass;
    const branch: plant_harvest_runtime.SourceOrderCompleteDeathBranchPools = .{
        .host_mobile = carbon(1),
        .symbiont_mobile = carbon(2),
        .c4_intermediate_carbon_g_c = 3,
        .leaf = carbon(4),
        .symbiont_structural = carbon(5),
        .sheath = carbon(6),
        .husk = carbon(7),
        .ear = carbon(8),
        .grain = carbon(9),
        .stalk = carbon(10),
        .stalk_reserve = carbon(11),
    };
    const base: plant_harvest_runtime.SourceOrderCompleteDeathShootInput = .{
        .shoot_dead = true,
        .roots_dead = true,
        .perennial_growth_habit = true,
        .deciduous_phenology = true,
        .seasonal_storage = carbon(8),
        .branches = &.{branch},
        .root_woody_fraction = composition,
        .leaf_woody_fraction = composition,
        .sheath_woody_fraction = composition,
        .nonstructural_kinetics = kinetics,
        .foliar_kinetics = kinetics,
        .nonfoliar_kinetics = kinetics,
        .stalk_kinetics = kinetics,
        .coarse_wood_kinetics = kinetics,
    };
    const perennial = try plant_harvest_runtime.sourceOrderCompleteDeathShootLitterfall(base);
    try std.testing.expect(perennial.plant_death_initialized);
    try std.testing.expectEqual(@as(f64, 2), perennial.planting_layer_woody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 6), perennial.planting_layer_nonwoody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 42.5), perennial.surface_nonwoody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 2.5), perennial.surface_woody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 21), perennial.standing_dead_stalk[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 0), perennial.seasonal_storage.carbon_g);

    var winter_input = base;
    winter_input.perennial_growth_habit = false;
    const winter = try plant_harvest_runtime.sourceOrderCompleteDeathShootLitterfall(winter_input);
    try std.testing.expect(!winter.plant_death_initialized);
    try std.testing.expectEqual(@as(f64, 33.5), winter.surface_nonwoody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 17), winter.seasonal_storage.carbon_g);
    try std.testing.expectEqual(@as(f64, 74), winter.surface_nonwoody_litter[0].carbon_g +
        winter.surface_woody_litter[0].carbon_g +
        winter.standing_dead_stalk[0].carbon_g +
        winter.seasonal_storage.carbon_g);
}

test "source-order complete death root litterfall conserves domains axes and layers" {
    const kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const mobile = [_]canopy.ElementalMass{
        .{ .carbon_g = 1, .nitrogen_g = 2, .phosphorus_g = 3 },
        .{ .carbon_g = 10, .nitrogen_g = 20, .phosphorus_g = 30 },
        .{ .carbon_g = 4, .nitrogen_g = 5, .phosphorus_g = 6 },
        .{ .carbon_g = 40, .nitrogen_g = 50, .phosphorus_g = 60 },
    };
    const axis: plant_harvest_runtime.SourceOrderDeadRootAxisPools = .{
        .primary = .{ .carbon_g = 2, .nitrogen_g = 2, .phosphorus_g = 2 },
        .secondary = .{ .carbon_g = 2, .nitrogen_g = 2, .phosphorus_g = 2 },
    };
    const structural = [_]plant_harvest_runtime.SourceOrderDeadRootAxisPools{
        axis, axis, axis, axis, axis, axis, axis, axis,
    };
    const base: plant_harvest_runtime.SourceOrderCompleteDeathRootInput = .{
        .shoot_dead = true,
        .roots_dead = true,
        .root_domain_count = 2,
        .soil_layer_count = 2,
        .root_axis_count = 2,
        .mobile_by_domain_layer = &mobile,
        .structural_by_domain_layer_axis = &structural,
        .root_woody_fraction = .{
            .carbon = .{ 0.25, 0.75 },
            .nitrogen = .{ 0.25, 0.75 },
            .phosphorus = .{ 0.25, 0.75 },
        },
        .nonstructural_kinetics = kinetics,
        .fine_root_kinetics = kinetics,
        .coarse_root_kinetics = kinetics,
    };
    const dead = try plant_harvest_runtime.sourceOrderCompleteDeathRootLitterfall(
        std.testing.allocator,
        base,
    );
    defer std.testing.allocator.free(dead);
    try std.testing.expectEqual(@as(f64, 4), dead[0].woody[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 17), dead[0].nonwoody[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 62), dead[1].nonwoody[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 21), dead[0].woody[0].carbon_g +
        dead[0].nonwoody[0].carbon_g);

    var partial_death = base;
    partial_death.shoot_dead = false;
    const live = try plant_harvest_runtime.sourceOrderCompleteDeathRootLitterfall(
        std.testing.allocator,
        partial_death,
    );
    defer std.testing.allocator.free(live);
    try std.testing.expectEqual(@as(f64, 0), live[0].nonwoody[0].carbon_g);
}

test "source-order complete death branch reset clears every runtime branch field" {
    const mass: canopy.ElementalMass = .{
        .carbon_g = 1,
        .nitrogen_g = 2,
        .phosphorus_g = 3,
    };
    const populated: plant_harvest_runtime.SourceOrderCompleteDeathBranchState = .{
        .host_mobile = mass,
        .c4_intermediate_carbon_g_c = 4,
        .symbiont_mobile = mass,
        .shoot = mass,
        .leaf = mass,
        .nodule = mass,
        .sheath = mass,
        .stalk = mass,
        .stalk_volume_m3 = 5,
        .reserve = mass,
        .husk = mass,
        .ear = mass,
        .grain = mass,
        .leaf_starch_carbon_g_c = 6,
        .stalk_extra = mass,
    };
    var partial = [_]plant_harvest_runtime.SourceOrderCompleteDeathBranchState{populated};
    try plant_harvest_runtime.sourceOrderResetCompleteDeathBranches(true, false, &partial);
    try std.testing.expectEqual(@as(f64, 1), partial[0].grain.carbon_g);

    var dead = [_]plant_harvest_runtime.SourceOrderCompleteDeathBranchState{ populated, populated };
    try plant_harvest_runtime.sourceOrderResetCompleteDeathBranches(true, true, &dead);
    for (dead) |branch| {
        inline for (@typeInfo(plant_harvest_runtime.SourceOrderCompleteDeathBranchState).@"struct".fields) |field| {
            const value = @field(branch, field.name);
            if (field.type == f64) {
                try std.testing.expectEqual(@as(f64, 0), value);
            } else inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element|
                try std.testing.expectEqual(@as(f64, 0), @field(value, element.name));
        }
    }
}

test "source-order complete death root reset preserves layer domain axis extents" {
    const mass: canopy.ElementalMass = .{
        .carbon_g = 1,
        .nitrogen_g = 2,
        .phosphorus_g = 3,
    };
    const axis: plant_harvest_runtime.SourceOrderDeadRootAxisLayerState = .{
        .primary = mass,
        .secondary = mass,
        .primary_length_m = 4,
        .secondary_length_m = 5,
        .secondary_axis_count = 6,
    };
    var mobile = [_]canopy.ElementalMass{ mass, mass, mass, mass };
    var structural = [_]plant_harvest_runtime.SourceOrderDeadRootAxisLayerState{
        axis, axis, axis, axis, axis, axis, axis, axis,
    };
    var totals = [_]plant_harvest_runtime.SourceOrderDeadRootDomainAxisState{
        .{ .primary_total = mass },
        .{ .primary_total = mass },
        .{ .primary_total = mass },
        .{ .primary_total = mass },
    };
    const state: plant_harvest_runtime.SourceOrderCompleteDeathRootResetState = .{
        .root_domain_count = 2,
        .soil_layer_count = 2,
        .root_axis_count = 2,
        .mobile_by_domain_layer = &mobile,
        .structural_by_domain_layer_axis = &structural,
        .primary_total_by_domain_axis = &totals,
    };
    try plant_harvest_runtime.sourceOrderResetCompleteDeathRoots(true, false, state);
    try std.testing.expectEqual(@as(f64, 1), mobile[0].carbon_g);

    try plant_harvest_runtime.sourceOrderResetCompleteDeathRoots(true, true, state);
    for (mobile) |value|
        try std.testing.expectEqual(@as(f64, 0), value.carbon_g);
    for (structural) |value| {
        try std.testing.expectEqual(@as(f64, 0), value.primary.carbon_g);
        try std.testing.expectEqual(@as(f64, 0), value.secondary.nitrogen_g);
        try std.testing.expectEqual(@as(f64, 0), value.primary_length_m);
        try std.testing.expectEqual(@as(f64, 0), value.secondary_length_m);
        try std.testing.expectEqual(@as(f64, 0), value.secondary_axis_count);
    }
    for (totals) |value|
        try std.testing.expectEqual(@as(f64, 0), value.primary_total.phosphorus_g);
}

test "source-order dead perennial reseed preserves next-day year rollover and gates" {
    const same_year = try plant_harvest_runtime.sourceOrderScheduleDeadPerennialReseed(
        true,
        false,
        200,
        365,
        2001,
    );
    try std.testing.expect(same_year.plant_death_flag);
    try std.testing.expectEqual(@as(u16, 201), same_year.reseed_date.?.day_of_year);
    try std.testing.expectEqual(@as(u32, 2001), same_year.reseed_date.?.year);

    const rollover = try plant_harvest_runtime.sourceOrderScheduleDeadPerennialReseed(
        true,
        false,
        366,
        366,
        2004,
    );
    try std.testing.expectEqual(@as(u16, 1), rollover.reseed_date.?.day_of_year);
    try std.testing.expectEqual(@as(u32, 2005), rollover.reseed_date.?.year);

    const annual = try plant_harvest_runtime.sourceOrderScheduleDeadPerennialReseed(
        false,
        false,
        100,
        365,
        2001,
    );
    try std.testing.expect(annual.reseed_date == null);
    try std.testing.expect(!annual.plant_death_flag);
    const terminated = try plant_harvest_runtime.sourceOrderScheduleDeadPerennialReseed(
        true,
        true,
        100,
        365,
        2001,
    );
    try std.testing.expect(terminated.reseed_date == null);
}

test "source-order soil plant exchange separates uptake fixation and NPP ledgers" {
    const result = try plant_harvest_runtime.sourceOrderAccumulateSoilPlantExchange(.{
        .organic_carbon_exchange_g_c_step = 1,
        .organic_nitrogen_exchange_g_n_step = 2,
        .ammonium_uptake_g_n_step = 3,
        .nitrate_uptake_g_n_step = 4,
        .root_fixation_g_n_step = 5,
        .canopy_fixation_g_n_step = 6,
        .organic_phosphorus_exchange_g_p_step = 7,
        .dihydrogen_phosphate_uptake_g_p_step = 8,
        .hydrogen_phosphate_uptake_g_p_step = 9,
        .cumulative_soil_exchange = .{
            .carbon_g = 10,
            .nitrogen_g = 20,
            .phosphorus_g = 30,
        },
        .cumulative_fixation_g_n = 40,
        .cumulative_plant_carbon_g_c = 50,
        .cumulative_respired_carbon_g_c = -12,
    });
    try std.testing.expectEqual(@as(f64, 1), result.hourly_net_exchange.carbon_g);
    try std.testing.expectEqual(@as(f64, 14), result.hourly_net_exchange.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 24), result.hourly_net_exchange.phosphorus_g);
    try std.testing.expectEqual(@as(f64, 11), result.cumulative_soil_exchange.carbon_g);
    try std.testing.expectEqual(@as(f64, 29), result.cumulative_soil_exchange.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 54), result.cumulative_soil_exchange.phosphorus_g);
    try std.testing.expectEqual(@as(f64, 51), result.cumulative_fixation_g_n);
    try std.testing.expectEqual(
        @as(f64, 38),
        result.cumulative_net_primary_productivity_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 5),
        result.hourly_net_exchange.nitrogen_g -
            (result.cumulative_soil_exchange.nitrogen_g - 20),
    );
}

test "source-order standing dead geometry aggregates runtime components and layers" {
    const components = [_]canopy.ElementalMass{
        .{ .carbon_g = 1.1416, .nitrogen_g = 2, .phosphorus_g = 3 },
        .{ .carbon_g = 2, .nitrogen_g = 4, .phosphorus_g = 5 },
    };
    const edges = [_]f64{ 0, 0.5, 1 };
    const result = try plant_harvest_runtime.sourceOrderStandingDeadGeometry(std.testing.allocator, .{
        .components = &components,
        .negligible_mass_g_c = 1.0e-12,
        .standing_dead_population_count = 1,
        .previous_height_m = 0.75,
        .canopy_height_m = 1,
        .stalk_volume_per_carbon_m3_g_c = 1,
        .canopy_layer_edges_m = &edges,
    });
    defer result.deinit(std.testing.allocator);
    try std.testing.expectApproxEqAbs(@as(f64, 3.1416), result.total_mass.carbon_g, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 6), result.total_mass.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 8), result.total_mass.phosphorus_g);
    try std.testing.expectApproxEqAbs(@as(f64, 1), result.height_m, 1.0e-12);
    try std.testing.expectApproxEqAbs(
        @as(f64, 6.2832),
        result.total_surface_area_m2,
        1.0e-12,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 3.1416),
        result.layer_surface_area_m2[0],
        1.0e-12,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.7854),
        result.projected_surface_area_m2[1],
        1.0e-12,
    );

    const empty = try plant_harvest_runtime.sourceOrderStandingDeadGeometry(std.testing.allocator, .{
        .components = &.{.{}},
        .negligible_mass_g_c = 1.0e-12,
        .standing_dead_population_count = 1,
        .previous_height_m = 1,
        .canopy_height_m = 1,
        .stalk_volume_per_carbon_m3_g_c = 1,
        .canopy_layer_edges_m = &edges,
    });
    defer empty.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(f64, 0), empty.height_m);
    try std.testing.expectEqual(@as(f64, 0), empty.layer_surface_area_m2[1]);
}

test "source-order fire inventory preserves plant layer domain and axis sums" {
    const carbon = struct {
        fn pools(value: f64) plant_harvest_runtime.SourceOrderDeadRootAxisPools {
            return .{
                .primary = .{ .carbon_g = value },
                .secondary = .{ .carbon_g = value },
            };
        }
    }.pools;
    const structural_a = [_]plant_harvest_runtime.SourceOrderDeadRootAxisPools{
        carbon(1), carbon(1), carbon(1), carbon(1),
        carbon(1), carbon(1), carbon(1), carbon(1),
    };
    const structural_b = [_]plant_harvest_runtime.SourceOrderDeadRootAxisPools{ carbon(2), carbon(3) };
    const shoot_a: plant_harvest_runtime.SourceOrderFireShootCarbon = .{
        .canopy_nonstructural_g_c = 1,
        .leaf_g_c = 2,
        .sheath_g_c = 3,
        .stalk_g_c = 4,
        .reserve_g_c = 5,
        .husk_g_c = 6,
        .ear_g_c = 7,
        .grain_g_c = 8,
        .symbiont_nonstructural_g_c = 9,
        .symbiont_biomass_g_c = 10,
        .seasonal_storage_g_c = 11,
        .standing_dead_g_c = 12,
    };
    var shoot_b = shoot_a;
    shoot_b.symbiont_nonstructural_g_c = 90;
    shoot_b.symbiont_biomass_g_c = 100;
    const plants = [_]plant_harvest_runtime.SourceOrderFirePlantInventory{
        .{
            .shoot = shoot_a,
            .canopy_symbiont_included = true,
            .root_domain_count = 2,
            .root_axis_count = 2,
            .nodule_nonstructural_by_layer_g_c = &.{ 1, 2 },
            .nodule_biomass_by_layer_g_c = &.{ 3, 4 },
            .root_nonstructural_by_layer_domain_g_c = &.{ 1, 2, 3, 4 },
            .root_structural_by_layer_domain_axis = &structural_a,
        },
        .{
            .shoot = shoot_b,
            .canopy_symbiont_included = false,
            .root_domain_count = 1,
            .root_axis_count = 1,
            .nodule_nonstructural_by_layer_g_c = &.{ 10, 20 },
            .nodule_biomass_by_layer_g_c = &.{ 30, 40 },
            .root_nonstructural_by_layer_domain_g_c = &.{ 10, 20 },
            .root_structural_by_layer_domain_axis = &structural_b,
        },
    };
    const optional = try plant_harvest_runtime.sourceOrderAggregateFireCarbonInventory(
        std.testing.allocator,
        true,
        2,
        &plants,
    );
    const result = optional.?;
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(f64, 2), result.shoot.canopy_nonstructural_g_c);
    try std.testing.expectEqual(@as(f64, 9), result.shoot.symbiont_nonstructural_g_c);
    try std.testing.expectEqual(@as(f64, 22), result.shoot.seasonal_storage_g_c);
    try std.testing.expectEqual(@as(f64, 13), result.layers[0].root_nonstructural_g_c);
    try std.testing.expectEqual(@as(f64, 12), result.layers[0].root_structural_g_c);
    try std.testing.expectEqual(@as(f64, 11), result.layers[0].nodule_nonstructural_g_c);
    try std.testing.expectEqual(@as(f64, 44), result.layers[1].nodule_biomass_g_c);
    try std.testing.expect((try plant_harvest_runtime.sourceOrderAggregateFireCarbonInventory(
        std.testing.allocator,
        false,
        0,
        &.{},
    )) == null);
}

test "source-order combustion rates preserve independent temperature gates" {
    const specific: plant_harvest_runtime.SourceOrderCombustionSpecificRates = .{
        .living_nonstructural_and_leaf_g_c_m2_h = 1,
        .living_sheath_g_c_m2_h = 2,
        .living_stalk_g_c_m2_h = 3,
        .living_reproductive_g_c_m2_h = 4,
        .standing_dead_g_c_m2_h = 5,
    };
    const living = try plant_harvest_runtime.sourceOrderCombustionRates(
        600,
        500,
        550,
        0.5,
        2,
        3,
        specific,
    );
    try std.testing.expectEqual(@as(f64, 0.5), living.living_temperature_fraction);
    try std.testing.expectEqual(@as(f64, 3), living.living_nonstructural_and_leaf_g_c_step);
    try std.testing.expectEqual(@as(f64, 6), living.living_sheath_g_c_step);
    try std.testing.expectEqual(@as(f64, 9), living.living_stalk_g_c_step);
    try std.testing.expectEqual(@as(f64, 12), living.living_reproductive_g_c_step);
    try std.testing.expectEqual(@as(f64, 0), living.standing_dead_g_c_step);

    const dead = try plant_harvest_runtime.sourceOrderCombustionRates(
        500,
        600,
        550,
        0.5,
        2,
        3,
        specific,
    );
    try std.testing.expectEqual(@as(f64, 0), dead.living_temperature_fraction);
    try std.testing.expectEqual(@as(f64, 0.5), dead.standing_dead_temperature_fraction);
    try std.testing.expectEqual(@as(f64, 15), dead.standing_dead_g_c_step);

    const cold = try plant_harvest_runtime.sourceOrderCombustionRates(
        550,
        550,
        550,
        0.5,
        2,
        3,
        specific,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        cold.living_nonstructural_and_leaf_g_c_step,
    );
    try std.testing.expectEqual(@as(f64, 0), cold.standing_dead_g_c_step);
}

test "source-order shoot combustion fractions preserve pool gates and routing" {
    const totals: plant_harvest_runtime.SourceOrderShootCombustionTotals = .{
        .canopy_nonstructural_g_c = 10,
        .leaf_g_c = 10,
        .sheath_g_c = 10,
        .stalk_g_c = 10,
        .husk_g_c = 10,
        .ear_g_c = 10,
        .grain_g_c = 10,
        .symbiont_nonstructural_g_c = 10,
        .symbiont_biomass_g_c = 10,
        .standing_dead_g_c = 10,
    };
    const rates: plant_harvest_runtime.SourceOrderCombustionRates = .{
        .living_temperature_fraction = 1,
        .standing_dead_temperature_fraction = 1,
        .living_nonstructural_and_leaf_g_c_step = 1,
        .living_sheath_g_c_step = 2,
        .living_stalk_g_c_step = 3,
        .living_reproductive_g_c_step = 4,
        .standing_dead_g_c_step = 5,
    };
    const woody = try plant_harvest_runtime.sourceOrderShootCombustionFractions(
        totals,
        rates,
        0,
        1,
        2,
    );
    try std.testing.expectEqual(@as(f64, 0.1), woody.canopy_nonstructural);
    try std.testing.expectEqual(@as(f64, 0.3), woody.sheath);
    try std.testing.expectEqual(@as(f64, 0.4), woody.stalk);
    try std.testing.expectEqual(woody.stalk, woody.reserve);
    try std.testing.expectEqual(@as(f64, 0.2), woody.ear);
    try std.testing.expectEqual(@as(f64, 0.5), woody.standing_dead);

    const herbaceous = try plant_harvest_runtime.sourceOrderShootCombustionFractions(
        totals,
        rates,
        0,
        0,
        2,
    );
    try std.testing.expectEqual(@as(f64, 0.2), herbaceous.sheath);
    try std.testing.expectEqual(@as(f64, 0.2), herbaceous.stalk);
    var sparse = totals;
    sparse.leaf_g_c = 1.0e-12;
    const gated = try plant_harvest_runtime.sourceOrderShootCombustionFractions(
        sparse,
        rates,
        1.0e-12,
        0,
        2,
    );
    try std.testing.expectEqual(@as(f64, 0), gated.leaf);
    inline for (@typeInfo(plant_harvest_runtime.SourceOrderShootCombustionFractions).@"struct".fields) |field| {
        const value = @field(woody, field.name);
        try std.testing.expect(value >= 0 and value <= 1);
    }
}

test "source-order branch combustion conserves C N P and signed ledgers" {
    const mass: canopy.ElementalMass = .{
        .carbon_g = 1,
        .nitrogen_g = 2,
        .phosphorus_g = 3,
    };
    const branch: plant_harvest_runtime.SourceOrderShootCombustionBranchPools = .{
        .canopy_nonstructural = mass,
        .leaf = mass,
        .sheath = mass,
        .stalk = mass,
        .reserve = mass,
        .husk = mass,
        .ear = mass,
        .grain = mass,
        .symbiont_nonstructural = mass,
        .symbiont_biomass = mass,
    };
    const fractions: plant_harvest_runtime.SourceOrderShootCombustionFractions = .{
        .canopy_nonstructural = 0.5,
        .leaf = 0.5,
        .sheath = 0.5,
        .stalk = 0.5,
        .reserve = 0.5,
        .husk = 0.5,
        .ear = 0.5,
        .grain = 0.5,
        .symbiont_nonstructural = 0.5,
        .symbiont_biomass = 0.5,
        .standing_dead = 0.5,
    };
    const result = try plant_harvest_runtime.sourceOrderShootCombustionLosses(
        std.testing.allocator,
        &.{ branch, branch },
        fractions,
        1,
        .{ .carbon_g = 100, .nitrogen_g = 200, .phosphorus_g = 300 },
    );
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 2), result.branches.len);
    try std.testing.expectEqual(@as(f64, 5), result.branches[0].total_combusted.carbon_g);
    try std.testing.expectEqual(@as(f64, 10), result.branches[0].total_combusted.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 15), result.branches[0].total_combusted.phosphorus_g);
    try std.testing.expectEqual(@as(f64, 11), result.cumulative_canopy_combustion_g_c);
    try std.testing.expectEqual(@as(f64, 90), result.disturbance_emission_ledger.carbon_g);
    try std.testing.expectEqual(@as(f64, 180), result.disturbance_emission_ledger.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 270), result.disturbance_emission_ledger.phosphorus_g);
    try std.testing.expectEqual(
        @as(f64, 10),
        100 - result.disturbance_emission_ledger.carbon_g,
    );
}

test "source-order shoot salt combustion conserves all eight species" {
    const branch: plant_harvest_runtime.SourceOrderShootSaltInventory = .{
        .aluminum_mol = 1,
        .iron_mol = 2,
        .calcium_mol = 3,
        .magnesium_mol = 4,
        .sodium_mol = 5,
        .potassium_mol = 6,
        .sulfate_mol = 7,
        .chloride_mol = 8,
    };
    const optional = try plant_harvest_runtime.sourceOrderShootSaltCombustion(
        std.testing.allocator,
        true,
        0.25,
        &.{ branch, branch },
    );
    const result = optional.?;
    defer std.testing.allocator.free(result);
    try std.testing.expectEqual(@as(usize, 2), result.len);
    try std.testing.expectEqual(@as(f64, 0.25), result[0].combusted.aluminum_mol);
    try std.testing.expectEqual(@as(f64, 2), result[0].combusted.chloride_mol);
    try std.testing.expectEqual(@as(f64, 6), result[0].remaining.chloride_mol);
    inline for (@typeInfo(plant_harvest_runtime.SourceOrderShootSaltInventory).@"struct".fields) |field| {
        try std.testing.expectEqual(
            @field(branch, field.name),
            @field(result[0].combusted, field.name) +
                @field(result[0].remaining, field.name),
        );
    }
    try std.testing.expect((try plant_harvest_runtime.sourceOrderShootSaltCombustion(
        std.testing.allocator,
        false,
        2,
        &.{},
    )) == null);
}

test "source-order uncombusted shoot state conserves pools and runtime topology" {
    const pool: canopy.ElementalMass = .{
        .carbon_g = 10,
        .nitrogen_g = 20,
        .phosphorus_g = 30,
    };
    const burned_mass: canopy.ElementalMass = .{
        .carbon_g = 2,
        .nitrogen_g = 4,
        .phosphorus_g = 6,
    };
    const pools: plant_harvest_runtime.SourceOrderShootCombustionBranchPools = .{
        .canopy_nonstructural = pool,
        .leaf = pool,
        .sheath = pool,
        .stalk = pool,
        .reserve = pool,
        .husk = pool,
        .ear = pool,
        .grain = pool,
        .symbiont_nonstructural = pool,
        .symbiont_biomass = pool,
    };
    const burned: plant_harvest_runtime.SourceOrderShootCombustionBranchPools = .{
        .canopy_nonstructural = burned_mass,
        .leaf = burned_mass,
        .sheath = burned_mass,
        .stalk = burned_mass,
        .reserve = burned_mass,
        .husk = burned_mass,
        .ear = burned_mass,
        .grain = burned_mass,
        .symbiont_nonstructural = burned_mass,
        .symbiont_biomass = burned_mass,
    };
    var nodes = [_]plant_harvest_runtime.SourceOrderShootCombustionNodeState{.{
        .leaf_area_m2 = 4,
        .sheath_height_m = 6,
        .green_leaf = pool,
        .senescent_leaf_carbon_g_c = 8,
        .green_sheath = pool,
        .senescent_sheath_carbon_g_c = 10,
        .node = pool,
    }};
    var layers = [_]plant_harvest_runtime.SourceOrderShootCombustionNodeLayerState{
        .{ .leaf_area_m2 = 2, .green_leaf = pool },
        .{ .leaf_area_m2 = 4, .green_leaf = pool },
    };
    var branches = [_]plant_harvest_runtime.SourceOrderUncombustedBranchState{.{
        .pools = pools,
        .c4_intermediate_carbon_g_c = 1,
        .total_shoot = .{},
        .leaf_area_m2 = 12,
        .nodes = &nodes,
        .node_layers = &layers,
        .canopy_layer_count = 2,
    }};
    var fractions = std.mem.zeroes(plant_harvest_runtime.SourceOrderShootCombustionFractions);
    fractions.leaf = 0.25;
    fractions.sheath = 0.5;
    fractions.stalk = 0.75;
    try plant_harvest_runtime.sourceOrderApplyUncombustedShootState(&branches, &.{burned}, fractions);
    try std.testing.expectEqual(@as(f64, 8), branches[0].pools.leaf.carbon_g);
    try std.testing.expectEqual(@as(f64, 65), branches[0].total_shoot.carbon_g);
    try std.testing.expectEqual(@as(f64, 128), branches[0].total_shoot.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 9), branches[0].leaf_area_m2);
    try std.testing.expectEqual(@as(f64, 3), nodes[0].leaf_area_m2);
    try std.testing.expectEqual(@as(f64, 3), nodes[0].sheath_height_m);
    try std.testing.expectEqual(@as(f64, 2.5), nodes[0].node.carbon_g);
    try std.testing.expectEqual(@as(f64, 1.5), layers[0].leaf_area_m2);
    try std.testing.expectEqual(@as(f64, 7.5), layers[1].green_leaf.carbon_g);
}

test "source-order standing dead combustion conserves components and ledgers" {
    const components = [_]canopy.ElementalMass{
        .{ .carbon_g = 1, .nitrogen_g = 2, .phosphorus_g = 3 },
        .{ .carbon_g = 2, .nitrogen_g = 4, .phosphorus_g = 6 },
        .{ .carbon_g = 3, .nitrogen_g = 6, .phosphorus_g = 9 },
        .{ .carbon_g = 4, .nitrogen_g = 8, .phosphorus_g = 12 },
    };
    const result = try plant_harvest_runtime.sourceOrderStandingDeadCombustion(
        std.testing.allocator,
        &components,
        0.25,
        10,
        .{ .carbon_g = 100, .nitrogen_g = 200, .phosphorus_g = 300 },
    );
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 4), result.components.len);
    try std.testing.expectEqual(@as(f64, 0.25), result.components[0].combusted.carbon_g);
    try std.testing.expectEqual(@as(f64, 3), result.components[3].remaining.carbon_g);
    try std.testing.expectEqual(
        @as(f64, 12.5),
        result.cumulative_standing_dead_combustion_g_c,
    );
    try std.testing.expectEqual(@as(f64, 97.5), result.disturbance_emission_ledger.carbon_g);
    try std.testing.expectEqual(@as(f64, 195), result.disturbance_emission_ledger.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 292.5), result.disturbance_emission_ledger.phosphorus_g);
    for (components, result.components) |before, after| {
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
            try std.testing.expectEqual(
                @field(before, field.name),
                @field(after.combusted, field.name) +
                    @field(after.remaining, field.name),
            );
    }
}

test "source-order charcoal combustion preserves response fraction and ledgers" {
    const result = try plant_harvest_runtime.sourceOrderCharcoalCombustion(
        700,
        0.5,
        2,
        3,
        4,
        24,
        1.0e-12,
        .{ .carbon_g = 10, .nitrogen_g = 2, .phosphorus_g = 1 },
        7,
        .{ .carbon_g = 100, .nitrogen_g = 20, .phosphorus_g = 10 },
        100,
        3,
    );
    try std.testing.expectEqual(@as(f64, 0.5), result.temperature_response);
    try std.testing.expectEqual(@as(f64, 12), result.potential_combustion_g_c_step);
    try std.testing.expectEqual(@as(f64, 0.5), result.combustion_fraction);
    try std.testing.expectEqual(@as(f64, 5), result.combusted.carbon_g);
    try std.testing.expectEqual(@as(f64, 5), result.remaining.carbon_g);
    try std.testing.expectEqual(@as(f64, 12), result.cumulative_standing_dead_combustion_g_c);
    try std.testing.expectEqual(@as(f64, 95), result.disturbance_emission_ledger.carbon_g);
    try std.testing.expectEqual(@as(f64, 19), result.disturbance_emission_ledger.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 9.5), result.disturbance_emission_ledger.phosphorus_g);
    try std.testing.expectEqual(@as(f64, 115), result.grid_total_combustion_g_c);
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        try std.testing.expectEqual(
            @field(result.combusted, field.name) +
                @field(result.remaining, field.name),
            @field(canopy.ElementalMass{
                .carbon_g = 10,
                .nitrogen_g = 2,
                .phosphorus_g = 1,
            }, field.name),
        );
}

test "source-order no-combustion branch resets every combustion rate" {
    const result = try plant_harvest_runtime.sourceOrderResetNoCombustion(std.testing.allocator, 3, true);
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(f64, 0), result.canopy_combustion_g_c_step);
    try std.testing.expectEqual(@as(f64, 0), result.standing_dead_combustion_g_c_step);
    try std.testing.expectEqual(@as(f64, 0), result.canopy_temperature_response);
    try std.testing.expectEqual(@as(f64, 0), result.standing_dead_temperature_response);
    try std.testing.expectEqual(@as(usize, 3), result.branch_combustion.len);
    try std.testing.expectEqual(@as(usize, 5), result.standing_dead_combustion.len);
    for (result.branch_combustion) |branch|
        inline for (@typeInfo(plant_harvest_runtime.SourceOrderShootCombustionBranchPools).@"struct".fields) |field|
            try std.testing.expectEqual(
                std.mem.zeroes(canopy.ElementalMass),
                @field(branch, field.name),
            );
    const salts = result.branch_salt_combustion.?;
    try std.testing.expectEqual(@as(usize, 3), salts.len);
    for (salts) |salt|
        inline for (@typeInfo(plant_harvest_runtime.SourceOrderShootSaltInventory).@"struct".fields) |field|
            try std.testing.expectEqual(@as(f64, 0), @field(salt, field.name));
    for (result.standing_dead_combustion) |component|
        try std.testing.expectEqual(std.mem.zeroes(canopy.ElementalMass), component);

    const without_salt = try plant_harvest_runtime.sourceOrderResetNoCombustion(
        std.testing.allocator,
        0,
        false,
    );
    defer without_salt.deinit(std.testing.allocator);
    try std.testing.expect(without_salt.branch_salt_combustion == null);
}

test "source-order root and surface storage combustion preserves source fractions" {
    const specific: plant_harvest_runtime.SourceOrderCombustionSpecificRates = .{
        .living_nonstructural_and_leaf_g_c_m2_h = 1,
        .living_sheath_g_c_m2_h = 2,
        .living_stalk_g_c_m2_h = 3,
        .living_reproductive_g_c_m2_h = 4,
        .standing_dead_g_c_m2_h = 5,
    };
    const result = try plant_harvest_runtime.sourceOrderRootStorageCombustion(.{
        .soil_temperature_k = 600,
        .minimum_combustion_temperature_k = 500,
        .maximum_temperature_response = 0.5,
        .surface_area_m2 = 2,
        .biological_timestep_h = 1,
        .specific_rates = specific,
        .totals = .{
            .root_nonstructural_g_c = 4,
            .active_root_g_c = 12,
            .nodule_nonstructural_g_c = 8,
            .nodule_biomass_g_c = 16,
        },
        .negligible_carbon_g_c = 1.0e-12,
        .is_surface_layer = true,
        .storage = .{ .carbon_g = 6, .nitrogen_g = 3, .phosphorus_g = 1.5 },
        .preceding_layer_combustion_g_c = 10,
        .preceding_disturbance_emission_ledger = .{ .carbon_g = 100, .nitrogen_g = 50, .phosphorus_g = 25 },
    });
    try std.testing.expectEqual(@as(f64, 0.5), result.temperature_response);
    try std.testing.expectEqual(@as(f64, 1), result.potential_rates.root_nonstructural_g_c_step);
    try std.testing.expectEqual(@as(f64, 2), result.potential_rates.nodule_biomass_g_c_step);
    try std.testing.expectEqual(@as(f64, 3), result.potential_rates.active_root_g_c_step);
    try std.testing.expectEqual(@as(f64, 0.25), result.fractions.root_nonstructural);
    try std.testing.expectEqual(@as(f64, 0.25), result.fractions.active_root);
    try std.testing.expectEqual(@as(f64, 0.125), result.fractions.nodule_nonstructural);
    try std.testing.expectEqual(@as(f64, 0.125), result.fractions.nodule_biomass);
    try std.testing.expectEqual(result.fractions.active_root, result.fractions.storage);
    try std.testing.expectEqual(@as(f64, 1.5), result.storage_combusted.carbon_g);
    try std.testing.expectEqual(@as(f64, 4.5), result.storage_remaining.carbon_g);
    try std.testing.expectEqual(@as(f64, 11.5), result.layer_combustion_g_c);
    try std.testing.expectEqual(@as(f64, 98.5), result.disturbance_emission_ledger.carbon_g);
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        try std.testing.expectEqual(
            @field(result.storage_combusted, field.name) +
                @field(result.storage_remaining, field.name),
            @field(canopy.ElementalMass{
                .carbon_g = 6,
                .nitrogen_g = 3,
                .phosphorus_g = 1.5,
            }, field.name),
        );

    const subsurface = try plant_harvest_runtime.sourceOrderRootStorageCombustion(.{
        .soil_temperature_k = 600,
        .minimum_combustion_temperature_k = 500,
        .maximum_temperature_response = 0.5,
        .surface_area_m2 = 2,
        .biological_timestep_h = 1,
        .specific_rates = specific,
        .totals = .{
            .root_nonstructural_g_c = 4,
            .active_root_g_c = 12,
            .nodule_nonstructural_g_c = 8,
            .nodule_biomass_g_c = 16,
        },
        .negligible_carbon_g_c = 1.0e-12,
        .is_surface_layer = false,
        .storage = .{ .carbon_g = 6, .nitrogen_g = 3, .phosphorus_g = 1.5 },
        .preceding_layer_combustion_g_c = 10,
        .preceding_disturbance_emission_ledger = .{ .carbon_g = 100, .nitrogen_g = 50, .phosphorus_g = 25 },
    });
    try std.testing.expectEqual(@as(f64, 0), subsurface.fractions.storage);
    try std.testing.expectEqual(@as(f64, 6), subsurface.storage_remaining.carbon_g);
    try std.testing.expectEqual(@as(f64, 10), subsurface.layer_combustion_g_c);
}

test "source-order root-domain combustion conserves runtime domains and axes" {
    var axes = [_]plant_harvest_runtime.SourceOrderRootCombustionAxisState{
        .{
            .primary = .{ .carbon_g = 8, .nitrogen_g = 4, .phosphorus_g = 2 },
            .secondary = .{ .carbon_g = 4, .nitrogen_g = 2, .phosphorus_g = 1 },
            .whole_primary = .{ .carbon_g = 16, .nitrogen_g = 8, .phosphorus_g = 4 },
            .primary_length_m = 10,
            .secondary_length_m = 6,
            .secondary_root_number = 2,
        },
        .{
            .primary = .{ .carbon_g = 12, .nitrogen_g = 6, .phosphorus_g = 3 },
            .secondary = .{ .carbon_g = 8, .nitrogen_g = 4, .phosphorus_g = 2 },
            .whole_primary = .{ .carbon_g = 20, .nitrogen_g = 10, .phosphorus_g = 5 },
            .primary_length_m = 14,
            .secondary_length_m = 8,
            .secondary_root_number = 4,
        },
    };
    var domains = [_]plant_harvest_runtime.SourceOrderRootCombustionDomainState{.{
        .nonstructural = .{ .carbon_g = 10, .nitrogen_g = 5, .phosphorus_g = 2 },
        .salts = .{
            .aluminum_mol = 1,
            .iron_mol = 2,
            .calcium_mol = 3,
            .magnesium_mol = 4,
            .sodium_mol = 5,
            .potassium_mol = 6,
            .sulfate_mol = 7,
            .chloride_mol = 8,
        },
        .active_root_carbon_g_c = 100,
        .root_density_g_c_m3 = 20,
        .root_surface_area_m2 = 30,
        .primary_root_number = 4,
        .root_length_m = 5,
        .root_length_growth_m_step = 6,
        .root_depth_growth_m_step = 7,
        .root_volume_growth_m3_step = 8,
        .root_volume_m3 = 9,
        .root_area_m2 = 10,
        .axes = &axes,
    }};
    const result = try plant_harvest_runtime.sourceOrderApplyRootDomainCombustion(
        std.testing.allocator,
        &domains,
        0.2,
        0.25,
        true,
        5,
        .{ .carbon_g = 100, .nitrogen_g = 50, .phosphorus_g = 25 },
    );
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 1), result.losses.len);
    try std.testing.expectEqual(@as(f64, 2), result.losses[0].nonstructural.carbon_g);
    try std.testing.expectEqual(@as(f64, 8), domains[0].nonstructural.carbon_g);
    try std.testing.expectEqual(@as(f64, 0.2), result.losses[0].salts.?.aluminum_mol);
    try std.testing.expectEqual(@as(f64, 0.8), domains[0].salts.aluminum_mol);
    try std.testing.expectEqual(@as(f64, 8), result.losses[0].structural.carbon_g);
    try std.testing.expectEqual(@as(f64, 6), axes[0].primary.carbon_g);
    try std.testing.expectEqual(@as(f64, 3), axes[0].secondary.carbon_g);
    try std.testing.expectEqual(@as(f64, 75), domains[0].active_root_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 7.5), axes[0].primary_length_m);
    try std.testing.expectEqual(@as(f64, 15), result.layer_combustion_g_c);
    try std.testing.expectEqual(@as(f64, 90), result.disturbance_emission_ledger.carbon_g);
    try std.testing.expectEqual(@as(f64, 45), result.disturbance_emission_ledger.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 22.6), result.disturbance_emission_ledger.phosphorus_g);
}

test "source-order root-nodule combustion conserves pools and source ledgers" {
    const result = try plant_harvest_runtime.sourceOrderRootNoduleCombustion(
        .{ .carbon_g = 8, .nitrogen_g = 4, .phosphorus_g = 2 },
        .{ .carbon_g = 12, .nitrogen_g = 6, .phosphorus_g = 3 },
        0.25,
        0.5,
        10,
        100,
        .{ .carbon_g = 50, .nitrogen_g = 25, .phosphorus_g = 12 },
    );
    try std.testing.expectEqual(@as(f64, 2), result.nonstructural_combusted.carbon_g);
    try std.testing.expectEqual(@as(f64, 6), result.nonstructural_remaining.carbon_g);
    try std.testing.expectEqual(@as(f64, 6), result.biomass_combusted.carbon_g);
    try std.testing.expectEqual(@as(f64, 6), result.biomass_remaining.carbon_g);
    try std.testing.expectEqual(@as(f64, 18), result.layer_plant_combustion_g_c);
    try std.testing.expectEqual(@as(f64, 118), result.grid_layer_combustion_g_c);
    try std.testing.expectEqual(@as(f64, 42), result.disturbance_emission_ledger.carbon_g);
    try std.testing.expectEqual(@as(f64, 21), result.disturbance_emission_ledger.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 10), result.disturbance_emission_ledger.phosphorus_g);
    inline for (.{ .{
        result.nonstructural_combusted,
        result.nonstructural_remaining,
        canopy.ElementalMass{ .carbon_g = 8, .nitrogen_g = 4, .phosphorus_g = 2 },
    }, .{
        result.biomass_combusted,
        result.biomass_remaining,
        canopy.ElementalMass{ .carbon_g = 12, .nitrogen_g = 6, .phosphorus_g = 3 },
    } }) |pools|
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
            try std.testing.expectEqual(
                @field(pools[2], field.name),
                @field(pools[0], field.name) + @field(pools[1], field.name),
            );
}

test "source-order cold-soil reset preserves heterogeneous runtime topology" {
    const axis_counts = [_]usize{ 2, 0, 1 };
    const result = try plant_harvest_runtime.sourceOrderResetColdSoilCombustion(
        std.testing.allocator,
        &axis_counts,
        true,
        true,
    );
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(f64, 0), result.layer_plant_combustion_g_c);
    try std.testing.expectEqual(std.mem.zeroes(canopy.ElementalMass), result.storage_combustion.?);
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 2, 3 }, result.axis_offsets);
    try std.testing.expectEqual(@as(usize, 3), result.domains.len);
    try std.testing.expectEqual(@as(usize, 3), result.axes.len);
    for (result.domains) |domain| {
        try std.testing.expectEqual(std.mem.zeroes(canopy.ElementalMass), domain.nonstructural);
        try std.testing.expectEqual(std.mem.zeroes(canopy.ElementalMass), domain.structural);
        const salts = domain.salts.?;
        inline for (@typeInfo(plant_harvest_runtime.SourceOrderShootSaltInventory).@"struct".fields) |field|
            try std.testing.expectEqual(@as(f64, 0), @field(salts, field.name));
    }
    for (result.axes) |axis| {
        try std.testing.expectEqual(std.mem.zeroes(canopy.ElementalMass), axis.primary);
        try std.testing.expectEqual(std.mem.zeroes(canopy.ElementalMass), axis.secondary);
    }
    try std.testing.expectEqual(
        std.mem.zeroes(canopy.ElementalMass),
        result.nodule_nonstructural_combustion,
    );
    try std.testing.expectEqual(
        std.mem.zeroes(canopy.ElementalMass),
        result.nodule_biomass_combustion,
    );

    const subsurface = try plant_harvest_runtime.sourceOrderResetColdSoilCombustion(
        std.testing.allocator,
        &.{},
        false,
        false,
    );
    defer subsurface.deinit(std.testing.allocator);
    try std.testing.expect(subsurface.storage_combustion == null);
    try std.testing.expectEqualSlices(usize, &.{0}, subsurface.axis_offsets);
}

test "source-order dormant-seed activation uses first qualifying branch" {
    const activation = (try plant_harvest_runtime.sourceOrderDormantSeedActivation(
        true,
        true,
        &.{
            .{
                .leafout_disabled = true,
                .accumulated_leafout_h = 100,
                .required_leafout_h = 10,
            },
            .{
                .leafout_disabled = false,
                .accumulated_leafout_h = 9,
                .required_leafout_h = 10,
            },
            .{
                .leafout_disabled = false,
                .accumulated_leafout_h = 10,
                .required_leafout_h = 10,
            },
            .{
                .leafout_disabled = false,
                .accumulated_leafout_h = 20,
                .required_leafout_h = 10,
            },
        },
        150,
        2025,
        0.02,
    )).?;
    try std.testing.expectEqual(@as(usize, 2), activation.qualifying_branch_index);
    try std.testing.expectEqual(@as(u16, 150), activation.planting_day_of_year);
    try std.testing.expectEqual(@as(i32, 2025), activation.planting_year);
    try std.testing.expectEqual(@as(f64, 0.025), activation.seeding_depth_m);
    try std.testing.expect(!activation.initialization_pending);
}

test "source-order dormant-seed activation preserves iteration and pending gates" {
    const ready = [_]plant_harvest_runtime.SourceOrderDormantSeedBranch{.{
        .leafout_disabled = false,
        .accumulated_leafout_h = 1,
        .required_leafout_h = 1,
    }};
    try std.testing.expect((try plant_harvest_runtime.sourceOrderDormantSeedActivation(
        false,
        true,
        &ready,
        1,
        2025,
        0,
    )) == null);
    try std.testing.expect((try plant_harvest_runtime.sourceOrderDormantSeedActivation(
        true,
        false,
        &ready,
        1,
        2025,
        0,
    )) == null);
    try std.testing.expectError(
        error.InvalidDormantSeedLeafoutHours,
        plant_harvest_runtime.sourceOrderDormantSeedActivation(
            true,
            true,
            &.{.{
                .leafout_disabled = false,
                .accumulated_leafout_h = std.math.nan(f64),
                .required_leafout_h = 1,
            }},
            1,
            2025,
            0,
        ),
    );
}

test "source-order litterfall accumulation preserves position fraction layer order" {
    var carbon: [20]f64 = undefined;
    var nitrogen: [20]f64 = undefined;
    var phosphorus: [20]f64 = undefined;
    for (0..20) |index| {
        carbon[index] = @floatFromInt(index + 1);
        nitrogen[index] = carbon[index] * 0.1;
        phosphorus[index] = carbon[index] * 0.01;
    }
    const result = try plant_harvest_runtime.sourceOrderAccumulateLitterfall(std.testing.allocator, .{
        .carbon_g_c_by_position_fraction_layer = &carbon,
        .nitrogen_g_n_by_position_fraction_layer = &nitrogen,
        .phosphorus_g_p_by_position_fraction_layer = &phosphorus,
        .layer_count_including_surface = 2,
        .preceding_cumulative_surface_litter = .{ .carbon_g = 5, .nitrogen_g = 0.5, .phosphorus_g = 0.05 },
        .preceding_hourly_litter = .{ .carbon_g = 7, .nitrogen_g = 0.7, .phosphorus_g = 0.07 },
        .preceding_cumulative_litter = .{ .carbon_g = 11, .nitrogen_g = 1.1, .phosphorus_g = 0.11 },
        .preceding_layer_carbon_g_c = &.{ 13, 17 },
    });
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(f64, 105), result.cumulative_surface_litter.carbon_g);
    try std.testing.expectApproxEqAbs(
        @as(f64, 10.5),
        result.cumulative_surface_litter.nitrogen_g,
        1.0e-14,
    );
    try std.testing.expectEqual(@as(f64, 217), result.hourly_litter.carbon_g);
    try std.testing.expectEqual(@as(f64, 221), result.cumulative_litter.carbon_g);
    try std.testing.expectEqualSlices(f64, &.{ 113, 127 }, result.layer_carbon_g_c);
}

test "source-order litterfall accumulation rejects a late invalid value atomically" {
    var values = [_]f64{0} ** 10;
    values[9] = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidLitterfallAccumulationInput,
        plant_harvest_runtime.sourceOrderAccumulateLitterfall(std.testing.allocator, .{
            .carbon_g_c_by_position_fraction_layer = &values,
            .nitrogen_g_n_by_position_fraction_layer = &([_]f64{0} ** 10),
            .phosphorus_g_p_by_position_fraction_layer = &([_]f64{0} ** 10),
            .layer_count_including_surface = 1,
            .preceding_cumulative_surface_litter = .{},
            .preceding_hourly_litter = .{},
            .preceding_cumulative_litter = .{},
            .preceding_layer_carbon_g_c = &.{0},
        }),
    );
}
