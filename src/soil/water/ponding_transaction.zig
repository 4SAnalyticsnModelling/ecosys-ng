const properties_module = @import("solver_properties.zig");
const organic_module = @import("../organic/initialization.zig");
const chemistry_module = @import("../solute/chemistry_state.zig");
const reactive_module = @import("../nutrients/reactive_nitrogen_state.zig");
const gas_module = @import("../gas/transport.zig");
const nitrogen_fertilizer_module = @import("../../management/fertilizer_nitrogen_inventory.zig");
const mineral_fertilizer_module = @import("../../management/mineral_fertilizer_inventory.zig");
const roots_module = @import("../../plant/root/plant_root_system.zig");
const grid_module = @import("../../state/grid.zig");
const thermal_module = @import("../heat/thermal.zig");
const mineral_remap = @import("../nutrients/mineral_layer_remap.zig");
const organic_remap = @import("../organic/layer_remap.zig");
const chemistry_remap = @import("../chemistry/layer_remap.zig");
const nitrite_remap = @import("../microbial/nitrite_layer_remap.zig");
const gas_remap = @import("../gas/layer_remap.zig");
const fertilizer_remap = @import("../nutrients/fertilizer_layer_remap.zig");
const root_remap = @import("../../plant/root/plant_root_layer_remap.zig");
const water_heat_remap = @import("heat_layer_remap.zig");
const pond_transition = @import("../../redistribution/pond/layer_transition.zig");

pub const Owners = struct {
    properties: *properties_module.State,
    organic: *organic_module.State,
    chemistry: *chemistry_module.State,
    reactive_nitrogen: *reactive_module.State,
    gas: *gas_module.State,
    nitrogen_fertilizer: *nitrogen_fertilizer_module.State,
    mineral_fertilizer: *mineral_fertilizer_module.State,
    roots: *roots_module.State,
    grid: *grid_module.GridState,
    thermal: *thermal_module.State,
};

pub const Inputs = struct {
    cell: usize,
    species_count: usize,
    plant_population_count: []const f64,
    root_structural_presence_g_per_plant: f64,
    source_layer: usize,
    destination_layer: usize,
    fraction: f64,
    source_soil_mass_before_megagrams: f64,
    destination_soil_mass_before_megagrams: f64,
    source_soil_mass_after_megagrams: f64,
    destination_soil_mass_after_megagrams: f64,
    source_water_before: chemistry_remap.ZoneWaterVolumes,
    destination_water_before: chemistry_remap.ZoneWaterVolumes,
    source_water_after: chemistry_remap.ZoneWaterVolumes,
    destination_water_after: chemistry_remap.ZoneWaterVolumes,
    dynamic_salts: bool,
    water_heat_parameters: water_heat_remap.Parameters,
};

/// Atomic boundary for all currently translated REDIST pond-content owners.
/// Every owner validates against the same pre/post carrier state before the
/// first commit. After that boundary, each commit is deterministic and cannot
/// fail unless an owner was mutated concurrently.
pub fn apply(owners: Owners, inputs: Inputs) !void {
    const layer_capacity = owners.nitrogen_fertilizer.layer_capacity;
    if (inputs.cell >= owners.nitrogen_fertilizer.cell_count or inputs.source_layer >= layer_capacity or inputs.destination_layer >= layer_capacity) return error.SoilPondingTransactionDimensionMismatch;
    const source = inputs.cell * layer_capacity + inputs.source_layer;
    const destination = inputs.cell * layer_capacity + inputs.destination_layer;
    if (source >= owners.properties.layer_count or source >= owners.organic.layer_count or source >= owners.chemistry.cell_count or source >= owners.reactive_nitrogen.layer_count or source >= owners.gas.cell_count or destination >= owners.properties.layer_count or destination >= owners.organic.layer_count or destination >= owners.chemistry.cell_count or destination >= owners.reactive_nitrogen.layer_count or destination >= owners.gas.cell_count) return error.SoilPondingTransactionDimensionMismatch;

    try water_heat_remap.validateLayerFraction(owners.grid, owners.thermal, source, destination, inputs.fraction, inputs.water_heat_parameters);
    try mineral_remap.validateLayerFraction(owners.properties, source, destination, inputs.fraction, inputs.source_soil_mass_after_megagrams, inputs.destination_soil_mass_after_megagrams);
    try organic_remap.validateLayerFraction(owners.organic, source, destination, inputs.fraction);
    try chemistry_remap.validateSolidLayerFraction(owners.chemistry, source, destination, inputs.source_soil_mass_before_megagrams, inputs.destination_soil_mass_before_megagrams, inputs.source_water_before.shared_m3, inputs.destination_water_before.shared_m3, inputs.source_soil_mass_after_megagrams, inputs.destination_soil_mass_after_megagrams, inputs.source_water_after.shared_m3, inputs.destination_water_after.shared_m3, inputs.fraction);
    try chemistry_remap.validateAqueousLayerFraction(owners.chemistry, source, destination, inputs.source_water_before, inputs.destination_water_before, inputs.source_water_after, inputs.destination_water_after, inputs.dynamic_salts, inputs.fraction);
    try nitrite_remap.validateLayerFraction(owners.reactive_nitrogen, source, destination, inputs.destination_water_before.nitrate_band_m3, inputs.fraction);
    try gas_remap.validateLayerFraction(
        owners.gas,
        source,
        destination,
        inputs.destination_water_before.ammonium_band_m3,
        inputs.fraction,
    );
    try fertilizer_remap.validateCellLayerFraction(owners.nitrogen_fertilizer, owners.mineral_fertilizer, inputs.cell, inputs.source_layer, inputs.destination_layer, inputs.fraction);
    try root_remap.validatePondedCellLayerFraction(
        owners.roots,
        inputs.cell,
        inputs.species_count,
        inputs.plant_population_count,
        inputs.root_structural_presence_g_per_plant,
        inputs.source_layer,
        inputs.destination_layer,
        inputs.fraction,
    );

    water_heat_remap.transferLayerFraction(owners.grid, owners.thermal, source, destination, inputs.fraction, inputs.water_heat_parameters) catch unreachable;
    mineral_remap.transferLayerFraction(owners.properties, source, destination, inputs.fraction, inputs.source_soil_mass_after_megagrams, inputs.destination_soil_mass_after_megagrams) catch unreachable;
    organic_remap.transferLayerFraction(owners.organic, source, destination, inputs.fraction) catch unreachable;
    chemistry_remap.transferSolidLayerFraction(owners.chemistry, source, destination, inputs.source_soil_mass_before_megagrams, inputs.destination_soil_mass_before_megagrams, inputs.source_water_before.shared_m3, inputs.destination_water_before.shared_m3, inputs.source_soil_mass_after_megagrams, inputs.destination_soil_mass_after_megagrams, inputs.source_water_after.shared_m3, inputs.destination_water_after.shared_m3, inputs.fraction) catch unreachable;
    chemistry_remap.transferAqueousLayerFraction(owners.chemistry, source, destination, inputs.source_water_before, inputs.destination_water_before, inputs.source_water_after, inputs.destination_water_after, inputs.dynamic_salts, inputs.fraction) catch unreachable;
    nitrite_remap.transferLayerFraction(owners.reactive_nitrogen, source, destination, inputs.destination_water_before.nitrate_band_m3, inputs.fraction) catch unreachable;
    gas_remap.transferLayerFraction(
        owners.gas,
        source,
        destination,
        inputs.destination_water_before.ammonium_band_m3,
        inputs.fraction,
    ) catch unreachable;
    fertilizer_remap.transferCellLayerFraction(owners.nitrogen_fertilizer, owners.mineral_fertilizer, inputs.cell, inputs.source_layer, inputs.destination_layer, inputs.fraction) catch unreachable;
    root_remap.transferPondedCellLayerFraction(
        owners.roots,
        inputs.cell,
        inputs.species_count,
        inputs.plant_population_count,
        inputs.root_structural_presence_g_per_plant,
        inputs.source_layer,
        inputs.destination_layer,
        inputs.fraction,
    ) catch unreachable;
}

/// Dispatches a validated selector result when both endpoints are runtime soil
/// layers. Surface-pond endpoints require the separate surface-owner adapter
/// and are rejected explicitly instead of aliasing the topsoil layer.
pub fn applySelectedSoilTransition(owners: Owners, base_inputs: Inputs, transition: pond_transition.Transition) !void {
    const source_layer = switch (transition.source) {
        .soil_layer => |layer| layer,
        .surface_pond => return error.SurfacePondAdapterRequired,
    };
    const destination_layer = switch (transition.destination) {
        .soil_layer => |layer| layer,
        .surface_pond => return error.SurfacePondAdapterRequired,
    };
    if (source_layer == destination_layer or transition.transfer_fraction <= 0) return error.InvalidSelectedPondTransition;
    var inputs = base_inputs;
    inputs.source_layer = source_layer;
    inputs.destination_layer = destination_layer;
    inputs.fraction = transition.transfer_fraction;
    try apply(owners, inputs);
}

test "REDIST combined pool transaction rejects late root failure before earlier owners mutate" {
    const std = @import("std");
    var properties: properties_module.State = undefined;
    properties.layer_count = 2;
    var sand_mass_megagrams = [_]f64{ 6, 2 };
    var silt_mass_megagrams = [_]f64{ 3, 5 };
    var clay_mass_megagrams = [_]f64{ 1, 3 };
    var cation_exchange_capacity_mol = [_]f64{ 100, 200 };
    var anion_exchange_capacity_mol = [_]f64{ 10, 20 };
    var sand_mass_fraction = [_]f64{ 0.6, 0.2 };
    var clay_mass_fraction = [_]f64{ 0.1, 0.3 };
    var cation_exchange_capacity_mol_per_megagram = [_]f64{ 10, 20 };
    var anion_exchange_capacity_mol_per_megagram = [_]f64{ 1, 2 };
    properties.sand_mass_megagrams = &sand_mass_megagrams;
    properties.silt_mass_megagrams = &silt_mass_megagrams;
    properties.clay_mass_megagrams = &clay_mass_megagrams;
    properties.cation_exchange_capacity_mol = &cation_exchange_capacity_mol;
    properties.anion_exchange_capacity_mol = &anion_exchange_capacity_mol;
    properties.sand_mass_fraction = &sand_mass_fraction;
    properties.clay_mass_fraction = &clay_mass_fraction;
    properties.cation_exchange_capacity_mol_per_megagram = &cation_exchange_capacity_mol_per_megagram;
    properties.anion_exchange_capacity_mol_per_megagram = &anion_exchange_capacity_mol_per_megagram;
    var organic = try organic_module.State.init(std.testing.allocator, 2);
    defer organic.deinit();
    var chemistry = try chemistry_module.State.init(std.testing.allocator, 2);
    defer chemistry.deinit();
    var reactive = try reactive_module.State.init(std.testing.allocator, 2, 1);
    defer reactive.deinit();
    var gas = try gas_module.State.init(std.testing.allocator, 2);
    defer gas.deinit();
    var nitrogen_fertilizer = try nitrogen_fertilizer_module.State.init(std.testing.allocator, 1, 2);
    defer nitrogen_fertilizer.deinit();
    var mineral_fertilizer = try mineral_fertilizer_module.State.init(std.testing.allocator, 1, 2);
    defer mineral_fertilizer.deinit();
    var roots = try roots_module.State.init(std.testing.allocator, 1, 2, 1);
    defer roots.deinit();
    const config = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 1, .lat_count = 1, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 4 });
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var thermal: thermal_module.State = undefined;
    var layer_volume_m3 = [_]f64{ 1, 1 };
    var dry_heat_capacity_megajoules_per_m3_k = [_]f64{ 1, 1 };
    var total_heat_capacity_megajoules_per_m3_k = [_]f64{ 1, 1 };
    var porosity_fraction = [_]f64{ 0.5, 0.5 };
    thermal.layer_volume_m3 = &layer_volume_m3;
    thermal.dry_solid_heat_capacity_megajoules_per_m3_k = &dry_heat_capacity_megajoules_per_m3_k;
    thermal.total_heat_capacity_megajoules_per_m3_k = &total_heat_capacity_megajoules_per_m3_k;
    thermal.porosity_fraction = &porosity_fraction;
    grid.matrix_liquid_water_m3[0] = 0.2;
    grid.matrix_liquid_water_m3[1] = 0.2;
    organic.microbial[0].carbon_g_c = 5;
    chemistry.aqueous[0].calcium = 4;
    reactive.non_band_nitrite_g_n[0] = 3;
    gas.gaseous_mass_g[0] = 2;
    nitrogen_fertilizer.soil[0].broadcast_ammonium_mol_n = 1;
    roots.total_carbon_g[try roots.layerIndex(0, 0, 0)] = 2;
    roots.total_carbon_g[try roots.layerIndex(0, 0, 1)] = 2;
    roots.gaseous_hydrogen_g_h[try roots.layerIndex(0, 1, 0)] = std.math.nan(f64);
    const water: chemistry_remap.ZoneWaterVolumes = .{ .shared_m3 = 1, .ammonium_non_band_m3 = 1, .ammonium_band_m3 = 0, .nitrate_non_band_m3 = 1, .nitrate_band_m3 = 0, .phosphate_non_band_m3 = 1, .phosphate_band_m3 = 0 };
    const inputs: Inputs = .{ .cell = 0, .species_count = 1, .plant_population_count = &.{5}, .root_structural_presence_g_per_plant = 0.1, .source_layer = 0, .destination_layer = 1, .fraction = 0.5, .source_soil_mass_before_megagrams = 10, .destination_soil_mass_before_megagrams = 10, .source_soil_mass_after_megagrams = 5, .destination_soil_mass_after_megagrams = 15, .source_water_before = water, .destination_water_before = water, .source_water_after = water, .destination_water_after = water, .dynamic_salts = false, .water_heat_parameters = .{ .liquid_water_heat_capacity_megajoules_per_m3_k = 4.19, .ice_heat_capacity_megajoules_per_m3_k = 1.9274, .minimum_heat_capacity_megajoules_per_k = 0 } };
    const owners: Owners = .{ .properties = &properties, .organic = &organic, .chemistry = &chemistry, .reactive_nitrogen = &reactive, .gas = &gas, .nitrogen_fertilizer = &nitrogen_fertilizer, .mineral_fertilizer = &mineral_fertilizer, .roots = &roots, .grid = &grid, .thermal = &thermal };
    try std.testing.expectError(error.InvalidRootLayerExtensiveState, apply(owners, inputs));
    try std.testing.expectEqual(@as(f64, 5), organic.microbial[0].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 4), chemistry.aqueous[0].calcium);
    try std.testing.expectEqual(@as(f64, 3), reactive.non_band_nitrite_g_n[0]);
    try std.testing.expectEqual(@as(f64, 2), gas.gaseous_mass_g[0]);
    try std.testing.expectEqual(@as(f64, 1), nitrogen_fertilizer.soil[0].broadcast_ammonium_mol_n);
    try std.testing.expectEqual(@as(f64, 0.2), grid.matrix_liquid_water_m3[0]);
    roots.gaseous_hydrogen_g_h[try roots.layerIndex(0, 1, 0)] = 0;
    try apply(owners, inputs);
    try std.testing.expectEqual(@as(f64, 2.5), organic.microbial[0].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 2), chemistry.aqueous[0].calcium);
    try std.testing.expectEqual(@as(f64, 1.5), reactive.non_band_nitrite_g_n[0]);
    try std.testing.expectEqual(@as(f64, 1), gas.gaseous_mass_g[0]);
    try std.testing.expectEqual(@as(f64, 0.5), nitrogen_fertilizer.soil[0].broadcast_ammonium_mol_n);
    try std.testing.expectEqual(@as(f64, 1), roots.total_carbon_g[try roots.layerIndex(0, 0, 0)]);
    try std.testing.expectEqual(@as(f64, 0.1), grid.matrix_liquid_water_m3[0]);
    try std.testing.expectEqual(@as(f64, 3), sand_mass_megagrams[0]);
    try std.testing.expectEqual(@as(f64, 5), sand_mass_megagrams[1]);
}

test "selected surface endpoint cannot be silently aliased to topsoil" {
    const std = @import("std");
    const owners: Owners = undefined;
    const inputs: Inputs = undefined;
    const transition: pond_transition.Transition = .{
        .source = .surface_pond,
        .destination = .{ .soil_layer = 0 },
        .transfer_fraction = 0.5,
        .boundary_change_m = -0.01,
        .next_first_active_layer = 0,
    };
    try std.testing.expectError(error.SurfacePondAdapterRequired, applySelectedSoilTransition(owners, inputs, transition));
}
