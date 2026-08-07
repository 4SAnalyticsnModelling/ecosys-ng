//! Tests for `landscape_mass_inventory.zig`.
//!
//! Extracted verbatim so the module beside it contains only the model
//! code. Tests that use private declarations of that module stay there,
//! since a sibling file can only reach `pub` declarations.

const audit = @import("mass_balance_audit.zig");
const canopy_retention = @import("../canopy/energy/precipitation_retention.zig");
const gas = @import("../soil/gas/transport.zig");
const grid_module = @import("../state/grid.zig");
const litter_chemistry = @import("../surface/litter_chemistry.zig");
const litter_fertilizer = @import("../surface/litter_fertilizer.zig");
const mineral_fertilizer = @import("../management/mineral_fertilizer_inventory.zig");
const mineral_nitrogen = @import("../soil/biogeochemistry/mineral_nitrogen_transport.zig");
const nitrogen_fertilizer = @import("../management/fertilizer_nitrogen_inventory.zig");
const organic = @import("../soil/organic/initialization.zig");
const organic_transport = @import("../soil/organic/transport.zig");
const plant_roots = @import("../plant/root/plant_root_system.zig");
const snow = @import("../soil/solute/snow_solute_transport.zig");
const soil_chemistry = @import("../soil/solute/chemistry_state.zig");
const solute_species = @import("../soil/solute/transport_species.zig");
const solute_transport = @import("../soil/solute/transport.zig");
const std = @import("std");
const surface_precipitation = @import("../surface/precipitation.zig");
const zone_classification = @import("../soil/solute/charge_classification.zig");
const landscape_mass_inventory = @import("landscape_mass_inventory.zig");
test "snow inventory rejects hidden non-finite inactive-layer mass" {
    var state = try snow.State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.amount_g[snow.species_count] = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteSnowInventory,
        landscape_mass_inventory.aggregateSnow(&state, 0.92),
    );
}

test "REDIST soil inventory uses runtime active layers and all gas phases" {
    const config = try @import("../core/config.zig").SimulationConfig.init(
        .{
            .lon_count = 2,
            .lat_count = 1,
            .soil_layers = 3,
            .plant_populations = 1,
        },
        .{ .worker_threads = 1, .tile_cells = 2 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    grid.active_soil_layer_count[0] = 2;
    grid.active_soil_layer_count[1] = 1;
    for (0..grid.layer_count) |index| {
        grid.matrix_liquid_water_m3[index] = 1;
        grid.macropore_liquid_water_m3[index] = 2;
        grid.matrix_ice_water_m3[index] = 3;
        grid.macropore_ice_water_m3[index] = 4;
        grid.water_vapor_volume_m3[index] = 5;
        grid.soil_temperature_k[index] = 250;
    }
    var gas_state = try gas.State.init(std.testing.allocator, grid.layer_count);
    defer gas_state.deinit();
    for (0..grid.layer_count) |layer| {
        const first = layer * gas.species_count;
        inline for (.{ gas_state.gaseous_mass_g, gas_state.dissolved_mass_g, gas_state.macropore_dissolved_mass_g, gas_state.band_dissolved_mass_g }) |phase| {
            phase[first + @intFromEnum(gas.Species.carbon_dioxide)] = 1;
            phase[first + @intFromEnum(gas.Species.methane)] = 2;
            phase[first + @intFromEnum(gas.Species.oxygen)] = 3;
            phase[first + @intFromEnum(gas.Species.nitrogen)] = 4;
            phase[first + @intFromEnum(gas.Species.nitrous_oxide)] = 5;
            phase[first + @intFromEnum(gas.Species.ammonia)] = 6;
        }
    }
    const dry_heat_capacity = [_]f64{2} ** 6;
    const volume = [_]f64{0.5} ** 6;
    const inventory = try landscape_mass_inventory.aggregateSoilPhysicalAndGas(
        &grid,
        &dry_heat_capacity,
        &volume,
        4,
        1.5,
        333,
        &gas_state,
    );
    const active: f64 = 3;
    try std.testing.expectEqual(active * 15, inventory.water_m3);
    // HEAT-001 resolution A, second layer: enthalpy relative to LIQUID water at
    // 0 K. Frozen carriers do not sit exactly `L` below liquid at the same
    // temperature; they follow the liquid branch to the melting point and the
    // ice branch back down, giving `C_l*Tm - L + C_i*(T - Tm)` per cubic metre.
    //
    // Note this test's heat capacities are the deliberately non-physical
    // `C_l = 4`, `C_i = 1.5`, which is useful here: the correction
    // `(C_l - C_i)*Tm` is then `2.5*273.15 = 682.875` per cubic metre rather
    // than the production `617.5`, so the expectation would not accidentally
    // agree if the code hardcoded production constants instead of using the
    // passed-in ones.
    const frozen_water_equivalent_m3: f64 = 3 + 4;
    try std.testing.expectEqual(
        active * ((2 * 0.5 + 4 * (1 + 2 + 5)) * 250 +
            (4 * 273.15 - 333 + 1.5 * (250 - 273.15)) * frozen_water_equivalent_m3),
        inventory.heat_megajoules,
    );
    try std.testing.expectEqual(active * 4 * 3, inventory.carbon_dioxide_carbon_g);
    try std.testing.expectEqual(active * 4 * 3, inventory.oxygen_g);
    try std.testing.expectEqual(active * 4 * 9, inventory.dinitrogen_nitrogen_g);
    try std.testing.expectEqual(active * 4 * 6, inventory.ammonium_nitrogen_g);
}

test "REDIST root gas inventory includes both root phases" {
    var roots = try plant_roots.State.init(std.testing.allocator, 1, 1, 1);
    defer roots.deinit();
    roots.gaseous_carbon_dioxide_g_c[0] = 1;
    roots.aqueous_carbon_dioxide_g_c[0] = 2;
    roots.gaseous_methane_g_c[0] = 3;
    roots.aqueous_methane_g_c[0] = 4;
    roots.gaseous_oxygen_g_o[0] = 5;
    roots.aqueous_oxygen_g_o[0] = 6;
    roots.gaseous_nitrous_oxide_g_n[0] = 7;
    roots.aqueous_nitrous_oxide_g_n[0] = 8;
    roots.gaseous_ammonia_g_n[0] = 9;
    roots.aqueous_ammonia_g_n[0] = 10;

    const inventory = try landscape_mass_inventory.aggregateRootGas(&roots);
    try std.testing.expectEqual(@as(f64, 10), inventory.carbon_dioxide_carbon_g);
    try std.testing.expectEqual(@as(f64, 11), inventory.oxygen_g);
    try std.testing.expectEqual(@as(f64, 15), inventory.dinitrogen_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 19), inventory.ammonium_nitrogen_g);
}

test "REDIST surface organic category split excludes humus microbial pool" {
    var state = try organic.State.init(std.testing.allocator, 2);
    defer state.deinit();
    // One counted microbial pool and one excluded K=4 humus microbial pool.
    state.microbial[0] = .{
        .carbon_g_c = 1,
        .nitrogen_g_n = 2,
        .phosphorus_g_p = 3,
    };
    const humus_microbial =
        4 * organic.microbial_population_count * organic.kinetic_fraction_count;
    state.microbial[humus_microbial] = .{
        .carbon_g_c = 100,
        .nitrogen_g_n = 100,
        .phosphorus_g_p = 100,
    };
    state.residue[0] = .{
        .carbon_g_c = 4,
        .nitrogen_g_n = 5,
        .phosphorus_g_p = 6,
    };
    state.dissolved[0] = .{
        .carbon_g_c = 7,
        .nitrogen_g_n = 8,
        .phosphorus_g_p = 9,
    };
    state.adsorbed[0] = .{
        .carbon_g_c = 10,
        .nitrogen_g_n = 11,
        .phosphorus_g_p = 12,
    };
    state.dissolved_acetate_carbon_g_c[0] = 13;
    state.adsorbed_acetate_carbon_g_c[0] = 14;
    // Fifth structural fraction is persistent charcoal and must be included.
    state.structural[organic.structural_fraction_count - 1] = .{
        .carbon_g_c = 15,
        .nitrogen_g_n = 16,
        .phosphorus_g_p = 17,
    };
    state.colonized_structural_carbon_g_c[
        organic.structural_fraction_count - 1
    ] = 9; // subset diagnostic; must not be counted a second time.
    const inventory = try landscape_mass_inventory.aggregateSurfaceOrganic(&state);
    try std.testing.expectEqual(@as(f64, 64), inventory.residue_carbon_g);
    try std.testing.expectEqual(@as(f64, 42), inventory.residue_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 47), inventory.residue_phosphorus_g);
    try std.testing.expectEqual(@as(f64, 0), inventory.organic_carbon_g);
}

test "REDIST soil organic split assigns only K=4 to humus" {
    const config = try @import("../core/config.zig").SimulationConfig.init(
        .{
            .lon_count = 1,
            .lat_count = 1,
            .soil_layers = 2,
            .plant_populations = 1,
        },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    grid.active_soil_layer_count[0] = 1;
    var state = try organic.State.init(std.testing.allocator, grid.layer_count);
    defer state.deinit();

    state.microbial[0] = .{
        .carbon_g_c = 1,
        .nitrogen_g_n = 2,
        .phosphorus_g_p = 3,
    };
    const humus_microbial =
        4 * organic.microbial_population_count * organic.kinetic_fraction_count;
    state.microbial[humus_microbial] = .{
        .carbon_g_c = 4,
        .nitrogen_g_n = 5,
        .phosphorus_g_p = 6,
    };
    const humus_mobile = 4;
    state.dissolved[humus_mobile] = .{
        .carbon_g_c = 7,
        .nitrogen_g_n = 8,
        .phosphorus_g_p = 9,
    };
    state.dissolved_acetate_carbon_g_c[humus_mobile] = 10;
    state.structural[organic.structural_fraction_count - 1] = .{
        .carbon_g_c = 11,
        .nitrogen_g_n = 12,
        .phosphorus_g_p = 13,
    };
    const humus_charcoal =
        4 * organic.structural_fraction_count +
        organic.structural_fraction_count - 1;
    state.structural[humus_charcoal] = .{
        .carbon_g_c = 14,
        .nitrogen_g_n = 15,
        .phosphorus_g_p = 16,
    };
    // Inactive capacity must not enter the profile inventory.
    const inactive_first =
        organic.microbial_substrate_count *
        organic.microbial_population_count *
        organic.kinetic_fraction_count;
    state.microbial[inactive_first].carbon_g_c = 1000;

    const inventory = try landscape_mass_inventory.aggregateSoilOrganic(&state, &grid);
    try std.testing.expectEqual(@as(f64, 12), inventory.residue_carbon_g);
    try std.testing.expectEqual(@as(f64, 14), inventory.residue_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 16), inventory.residue_phosphorus_g);
    try std.testing.expectEqual(@as(f64, 35), inventory.organic_carbon_g);
    try std.testing.expectEqual(@as(f64, 28), inventory.organic_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 31), inventory.organic_phosphorus_g);
}

test "REDIST surface chemistry retains N P and TION stoichiometry" {
    var chemistry = try litter_chemistry.State.init(std.testing.allocator, 1);
    defer chemistry.deinit();
    var fertilizer = try litter_fertilizer.State.init(std.testing.allocator, 1);
    defer fertilizer.deinit();
    const cell = &chemistry.cells[0];
    cell.ammonium_mol_per_m3 = 1;
    cell.ammonia_mol_per_m3 = 2;
    cell.nitrate_mol_per_m3 = 3;
    cell.hpo4_mol_p_per_m3 = 4;
    cell.h2po4_mol_p_per_m3 = 5;
    cell.calcium_mol_per_m3 = 6;
    cell.bicarbonate_mol_per_m3 = 7;
    cell.exchange.ammonium_mol_per_megagram = 8;
    cell.exchange.calcium_mol_per_megagram = 9;
    cell.phosphate_surface.adsorbed_hpo4_mol_p_per_megagram = 10;
    cell.phosphate_surface.adsorbed_h2po4_mol_p_per_megagram = 11;
    cell.phosphate_minerals.monocalcium_phosphate_mol_per_m3 = 12;
    cell.phosphate_minerals.hydroxyapatite_mol_per_m3 = 13;
    cell.salt_minerals.gypsum_mol_per_m3 = 14;
    fertilizer.cells[0].ammonium_mol_n = 15;
    fertilizer.cells[0].ammonia_mol_n = 16;
    fertilizer.cells[0].urea_mol_n = 17;
    fertilizer.cells[0].nitrate_mol_n = 18;

    const inventory = try landscape_mass_inventory.aggregateSurfaceChemistry(
        &chemistry,
        &fertilizer,
        &.{19},
        &.{2},
        &.{3},
        12,
        14,
        31,
    );
    try std.testing.expectEqual(
        @as(f64, 14 * (2 * (1 + 2) + 3 * 8 + 15 + 16 + 17)),
        inventory.ammonium_nitrogen_g,
    );
    try std.testing.expectEqual(
        @as(f64, 14 * (2 * 3 + 18) + 19),
        inventory.nitrate_nitrogen_g,
    );
    try std.testing.expectEqual(
        @as(f64, 31 * (2 * (4 + 5 + 2 * 12 + 3 * 13) + 3 * (10 + 11))),
        inventory.phosphate_phosphorus_g,
    );
    try std.testing.expectEqual(
        @as(f64, 12 * 2 * 7),
        inventory.carbon_dioxide_carbon_g,
    );
    const expected_ions =
        2 * ((6 + 2 * 7 + 3 * 4 + 4 * 5) +
            (2 * 14 + 7 * 12 + 9 * 13)) +
        3 * (2 * 8 + 9 + 3 * 10 + 4 * 11) +
        (2 * 15 + 16 + 17 + 18);
    try std.testing.expectEqual(
        @as(f64, expected_ions),
        inventory.ion_inventory_mol,
    );
}

test "storage publication cannot overwrite cumulative EXEC ledgers" {
    var totals = std.mem.zeroes(audit.Totals);
    totals.landscape_area_m2 = 10;
    totals.cumulative_rain_m3 = 17;
    totals.cumulative_nitrogen_output_g = 19;
    try landscape_mass_inventory.publishStorage(&totals, .{
        .water_m3 = 1,
        .organic_carbon_g = 2,
        .ion_inventory_mol = 3,
    });
    try std.testing.expectEqual(@as(f64, 1), totals.water_storage_m3);
    try std.testing.expectEqual(@as(f64, 2), totals.organic_carbon_g);
    try std.testing.expectEqual(@as(f64, 3), totals.ion_inventory_mol);
    try std.testing.expectEqual(@as(f64, 17), totals.cumulative_rain_m3);
    try std.testing.expectEqual(
        @as(f64, 19),
        totals.cumulative_nitrogen_output_g,
    );
}

test "REDIST surface physical inventory includes vapor and four gas phases" {
    var surface = try surface_precipitation.RuntimeState.init(
        std.testing.allocator,
        2,
    );
    defer surface.deinit();
    const config = try @import("../core/config.zig").SimulationConfig.init(
        .{
            .lon_count = 2,
            .lat_count = 1,
            .soil_layers = 1,
            .plant_populations = 1,
        },
        .{ .worker_threads = 1, .tile_cells = 2 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var gas_state = try gas.State.init(std.testing.allocator, 2);
    defer gas_state.deinit();
    var organic_state = try organic.State.init(std.testing.allocator, 2);
    defer organic_state.deinit();
    surface.litter_water_m3[0] = 1;
    surface.litter_water_m3[1] = 2;
    grid.surface_temperature_k[0] = 250;
    grid.surface_temperature_k[1] = 300;
    gas_state.water_vapor_mol[0] = 10;
    gas_state.water_vapor_mol[1] = 20;
    organic_state.structural[0].carbon_g_c = 100;
    organic_state.structural[
        organic.substrate_count * organic.structural_fraction_count
    ].carbon_g_c = 200;
    for (0..2) |cell| {
        const first = cell * gas.species_count;
        inline for (.{ gas_state.gaseous_mass_g, gas_state.dissolved_mass_g, gas_state.macropore_dissolved_mass_g, gas_state.band_dissolved_mass_g }) |phase| {
            phase[first + @intFromEnum(gas.Species.carbon_dioxide)] = 1;
            phase[first + @intFromEnum(gas.Species.methane)] = 2;
            phase[first + @intFromEnum(gas.Species.oxygen)] = 3;
            phase[first + @intFromEnum(gas.Species.nitrogen)] = 4;
            phase[first + @intFromEnum(gas.Species.nitrous_oxide)] = 5;
            phase[first + @intFromEnum(gas.Species.ammonia)] = 6;
        }
    }
    const ice = [_]f64{ 0.5, 0.25 };
    const parameters: landscape_mass_inventory.SurfacePhysicalParameters = .{
        .dry_organic_heat_capacity_megajoules_per_g_c_k = 2.5e-6,
        .liquid_water_heat_capacity_megajoules_per_m3_k = 4.19,
        .ice_heat_capacity_megajoules_per_m3_k = 1.9274,
        .water_molar_mass_g_per_mol = 18,
        .liquid_water_density_g_per_m3 = 1e6,
    };
    const inventory = try landscape_mass_inventory.aggregateSurfacePhysicalAndGas(
        &surface,
        &ice,
        &grid,
        &gas_state,
        &organic_state,
        parameters,
    );
    const vapor0 = 10.0 * 18.0 / 1e6;
    const vapor1 = 20.0 * 18.0 / 1e6;
    try std.testing.expectApproxEqAbs(
        1 + 2 + 0.5 + 0.25 + vapor0 + vapor1,
        inventory.water_m3,
        1e-14,
    );
    // HEAT-001 resolution A, second layer. Surface litter/pond ice is a frozen
    // carrier and carries `C_l*Tm - L + C_i*(T - Tm)` per cubic metre, not
    // `C_i*T - L`. The two cells sit at 250 K and 300 K, so this expectation
    // also pins that the ice branch is evaluated at each cell's own
    // temperature rather than at the melting point.
    const frozen_enthalpy_at = struct {
        fn f(temperature_k: f64) f64 {
            return 4.19 * 273.15 - 333 + 1.9274 * (temperature_k - 273.15);
        }
    }.f;
    const expected_heat =
        (2.5e-6 * 100 + 4.19 * (1 + vapor0)) * 250 +
        (2.5e-6 * 200 + 4.19 * (2 + vapor1)) * 300 +
        frozen_enthalpy_at(250) * 0.5 +
        frozen_enthalpy_at(300) * 0.25;
    try std.testing.expectApproxEqAbs(expected_heat, inventory.heat_megajoules, 1e-10);
    try std.testing.expectEqual(@as(f64, 24), inventory.carbon_dioxide_carbon_g);
    try std.testing.expectEqual(@as(f64, 24), inventory.oxygen_g);
    try std.testing.expectEqual(@as(f64, 72), inventory.dinitrogen_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 48), inventory.ammonium_nitrogen_g);
}

test "EXTRACT canopy inventory supports more than five runtime species" {
    const config = try @import("../core/config.zig").SimulationConfig.init(
        .{
            .lon_count = 2,
            .lat_count = 1,
            .soil_layers = 1,
            .plant_populations = 7,
        },
        .{ .worker_threads = 3, .tile_cells = 2 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var plants = try grid_module.PlantState.init(std.testing.allocator, config);
    defer plants.deinit();
    var retention = try canopy_retention.State.init(
        std.testing.allocator,
        2,
        7,
    );
    defer retention.deinit();
    for (0..14) |plant| {
        plants.canopy_water_storage_m_per_m2[plant] = 0.001;
        retention.living_surface_water_m3[plant] = 0.01;
        retention.standing_dead_surface_water_m3[plant] = 0.02;
        retention.previous_water_energy_megajoules[plant] = 3;
    }
    const inventory = try landscape_mass_inventory.aggregateCanopyWaterAndHeat(
        &plants,
        &retention,
        &.{ 10, 20 },
    );
    try std.testing.expectApproxEqAbs(
        7 * (0.001 * 10 + 0.01 + 0.02) +
            7 * (0.001 * 20 + 0.01 + 0.02),
        inventory.water_m3,
        1e-14,
    );
    try std.testing.expectEqual(@as(f64, 42), inventory.heat_megajoules);
}

test "REDIST profile mineral nitrogen counts each runtime owner once" {
    const config = try @import("../core/config.zig").SimulationConfig.init(
        .{
            .lon_count = 2,
            .lat_count = 1,
            .soil_layers = 2,
            .plant_populations = 1,
        },
        .{ .worker_threads = 2, .tile_cells = 2 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    grid.active_soil_layer_count[0] = 1;
    grid.active_soil_layer_count[1] = 1;

    var transport = try mineral_nitrogen.State.init(
        std.testing.allocator,
        grid.layer_count,
    );
    defer transport.deinit();
    var chemistry = try soil_chemistry.State.init(
        std.testing.allocator,
        grid.layer_count,
    );
    defer chemistry.deinit();
    var fertilizer = try nitrogen_fertilizer.State.init(
        std.testing.allocator,
        grid.cell_count,
        grid.soil_layer_capacity,
    );
    defer fertilizer.deinit();

    // Active profile cells are 0 and 2. Matrix and macropore are both
    // authoritative and must be included.
    for ([_]usize{ 0, 2 }) |profile_cell| {
        const matrix = try transport.matrix.cellAmounts(profile_cell);
        const macropore = try transport.macropore.cellAmounts(profile_cell);
        matrix[@intFromEnum(mineral_nitrogen.Species.ammonium_non_band)] = 1;
        macropore[@intFromEnum(mineral_nitrogen.Species.ammonia_band)] = 2;
        matrix[@intFromEnum(mineral_nitrogen.Species.nitrate_non_band)] = 3;
        macropore[@intFromEnum(mineral_nitrogen.Species.nitrite_band)] = 4;
        chemistry.cation_exchange_mol_per_megagram[profile_cell]
            .ammonium_non_band = 0.5;
        fertilizer.soil[profile_cell].broadcast_ammonium_mol_n = 5;
        fertilizer.soil[profile_cell].banded_urea_mol_n = 6;
        fertilizer.soil[profile_cell].broadcast_nitrate_mol_n = 7;
    }
    // Inactive allocated capacity must never enter an authoritative total.
    (try transport.matrix.cellAmounts(1))[0] = 1_000_000;
    chemistry.cation_exchange_mol_per_megagram[3].ammonium_band = 1_000_000;
    fertilizer.soil[1].broadcast_ammonium_mol_n = 1_000_000;

    const inventory = try landscape_mass_inventory.aggregateProfileMineralNitrogen(
        &grid,
        &transport,
        &chemistry,
        &fertilizer,
        &.{ 2, 2, 2, 2 },
        14,
    );
    // Per active cell: NHx=(1+2)+(0.5*2)+(5+6)=15 mol N;
    // NOx=(3+4)+7=14 mol N.
    try std.testing.expectEqual(@as(f64, 2 * 15 * 14), inventory.ammonium_nitrogen_g);
    try std.testing.expectEqual(@as(f64, 2 * 14 * 14), inventory.nitrate_nitrogen_g);
    // Per cell: exchange NH4 2*1 + dry NH4 2*5 + urea 6 + nitrate 7.
    try std.testing.expectEqual(@as(f64, 2 * 25), inventory.ion_inventory_mol);
}

test "REDIST profile phosphorus and ions include matrix macropore and immobile owners" {
    const config = try @import("../core/config.zig").SimulationConfig.init(
        .{
            .lon_count = 2,
            .lat_count = 1,
            .soil_layers = 2,
            .plant_populations = 1,
        },
        .{ .worker_threads = 2, .tile_cells = 2 },
        .{
            .relative_tolerance = 1e-8,
            .absolute_tolerance = 1e-11,
            .max_nonlinear_iterations = 20,
        },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    grid.active_soil_layer_count[0] = 1;
    grid.active_soil_layer_count[1] = 1;

    var micropore = try solute_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
        solute_species.AqueousSpecies.count,
    );
    defer micropore.deinit();
    var macropore = try solute_transport.State.init(
        std.testing.allocator,
        grid.layer_count,
        solute_species.AqueousSpecies.count,
    );
    defer macropore.deinit();
    var chemistry = try soil_chemistry.State.init(
        std.testing.allocator,
        grid.layer_count,
    );
    defer chemistry.deinit();
    var pending = try mineral_fertilizer.State.init(
        std.testing.allocator,
        grid.cell_count,
        grid.soil_layer_capacity,
    );
    defer pending.deinit();

    for ([_]usize{ 0, 2 }) |profile_cell| {
        const matrix = try micropore.cellAmounts(profile_cell);
        const macro = try macropore.cellAmounts(profile_cell);
        matrix[solute_species.index(.aluminum)] = 1;
        matrix[solute_species.index(.bicarbonate)] = 2;
        matrix[solute_species.index(.aluminum_hydroxide_2)] = 3;
        matrix[solute_species.index(.aluminum_hydroxide_3)] = 4;
        matrix[solute_species.index(.aluminum_hydroxide_4)] = 5;
        matrix[solute_species.index(.non_band_phosphate)] = 6;
        matrix[solute_species.index(.band_iron_hpo4)] = 7;
        macro[solute_species.index(.band_phosphoric_acid)] = 8;

        const non_band = &chemistry.non_band_phosphate[profile_cell];
        non_band.dissolved_hpo4_mol_p_per_m3 = 1;
        non_band.adsorbed_hpo4_mol_p_per_megagram = 2;
        non_band.deprotonated_site_mol_per_megagram = 1;
        non_band.monocalcium_phosphate_solid_mol_per_m3 = 1;
        const band = &chemistry.band_phosphate[profile_cell];
        band.dissolved_h2po4_mol_p_per_m3 = 2;
        band.adsorbed_h2po4_mol_p_per_megagram = 1;
        band.hydroxyl_site_mol_per_megagram = 1;
        band.hydroxyapatite_solid_mol_per_m3 = 2;

        chemistry.cation_exchange_mol_per_megagram[profile_cell].hydrogen = 1;
        chemistry.cation_exchange_mol_per_megagram[profile_cell].calcium = 2;
        chemistry.carboxyl_bound_hydrogen_mol_per_megagram[profile_cell] = 0.5;
        chemistry.geochemistry_solids[profile_cell]
            .calcite_solid_mol_per_m3 = 1;
        chemistry.geochemistry_solids[profile_cell]
            .gibbsite_solid_mol_per_m3 = 2;
        chemistry.geochemistry_solids[profile_cell]
            .aluminum_natural_silicate_mol_per_m3 = 3;

        pending.soil[profile_cell].broadcast_monocalcium_phosphate_mol = 1;
        pending.soil[profile_cell].hydroxyapatite_mol = 2;
        pending.soil[profile_cell].calcite_mol = 3;
        pending.soil[profile_cell].aluminum_ground_silicate_mol = 4;
    }
    // Capacity layers 1 and 3 are not part of the runtime profile.
    (try micropore.cellAmounts(1))[solute_species.index(.band_phosphate)] =
        1_000_000;
    chemistry.non_band_phosphate[3]
        .hydroxyapatite_solid_mol_per_m3 = 1_000_000;
    pending.soil[1].hydroxyapatite_mol = 1_000_000;

    const inventory = try landscape_mass_inventory.aggregateProfilePhosphorusAndIons(
        &grid,
        &micropore,
        &macropore,
        &chemistry,
        &pending,
        &.{ 3, 3, 3, 3 },
        &.{ 2, 2, 2, 2 },
        .{
            .ammonium_non_band = 0.75,
            .ammonium_band = 0.25,
            .nitrate_non_band = 0.75,
            .nitrate_band = 0.25,
            .phosphate_non_band = 0.75,
            .phosphate_band = 0.25,
        },
        12,
        31,
    );
    try std.testing.expectEqual(
        @as(f64, 2 * 45.25 * 31),
        inventory.phosphate_phosphorus_g,
    );
    try std.testing.expectEqual(
        @as(f64, 2 * 250.5),
        inventory.ion_inventory_mol,
    );
    try std.testing.expectEqual(
        @as(f64, 2 * 96),
        inventory.carbon_dioxide_carbon_g,
    );
}

test "REDIST dry surface mineral fertilizer remains in landscape storage" {
    var state = try mineral_fertilizer.State.init(
        std.testing.allocator,
        2,
        3,
    );
    defer state.deinit();
    state.surface[0].broadcast_monocalcium_phosphate_mol = 2;
    state.surface[0].banded_monocalcium_phosphate_mol = 3;
    state.surface[0].hydroxyapatite_mol = 4;
    state.surface[0].calcite_mol = 5;
    state.surface[0].gypsum_mol = 6;
    state.surface[0].aluminum_ground_silicate_mol = 7;
    state.surface[1].potassium_ground_silicate_mol = 8;

    const inventory = try landscape_mass_inventory.aggregatePendingSurfaceMinerals(&state, 12, 31);
    try std.testing.expectEqual(
        @as(f64, (2 * (2 + 3) + 3 * 4) * 31),
        inventory.phosphate_phosphorus_g,
    );
    try std.testing.expectEqual(
        @as(f64, 5 * 12),
        inventory.carbon_dioxide_carbon_g,
    );
    try std.testing.expectEqual(
        @as(f64, 7 * (2 + 3) + 9 * 4 + 2 * (5 + 6) + 7 + 8),
        inventory.ion_inventory_mol,
    );
}
