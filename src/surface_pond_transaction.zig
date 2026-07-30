const inventory = @import("surface_pond_inventory_transfer.zig");
const chemistry = @import("surface_pond_chemistry_transfer.zig");
const surface_chemistry_module = @import("surface_litter_chemistry.zig");
const soil_chemistry_module = @import("solute_chemistry_state.zig");
const water_heat = @import("surface_pond_water_heat_transfer.zig");
const geometry_module = @import("soil_layer_geometry.zig");

pub const Owners = struct {
    inventories: inventory.Owners,
    surface_chemistry: *surface_chemistry_module.State,
    soil_chemistry: *soil_chemistry_module.State,
    water_heat: water_heat.Owners,
    soil_geometry: *geometry_module.State,
};

pub const Inputs = struct {
    cell: usize,
    destination_soil_layer: usize,
    fraction: f64,
    chemistry_carriers: chemistry.CarrierVolumes,
    dynamic_salts: bool,
    water_heat_parameters: water_heat.Parameters,
    geometry_changes: geometry_module.DisturbanceChanges,
    minimum_soil_layer_thickness_m: f64,
};

/// Atomic REDIST surface-pond inventory boundary. Both distinct owner schemas
/// validate before either commits.
pub fn transferSurfaceFractionToSoil(owners: Owners, inputs: Inputs) !void {
    const destination = inputs.cell * owners.inventories.soil_nitrogen_fertilizer.layer_capacity + inputs.destination_soil_layer;
    const surface_organic_carbon_g_c = try owners.inventories.surface_organic.totalCarbon_g_c(inputs.cell);
    try inventory.validateSurfaceFractionToSoil(owners.inventories, .{ .cell = inputs.cell, .destination_soil_layer = inputs.destination_soil_layer, .fraction = inputs.fraction });
    try chemistry.validateSurfaceFractionToSoil(owners.surface_chemistry, owners.soil_chemistry, inputs.cell, destination, inputs.chemistry_carriers, inputs.dynamic_salts, inputs.fraction);
    try water_heat.validateSurfaceFractionToSoil(owners.water_heat, .{ .cell = inputs.cell, .destination_soil_layer = inputs.destination_soil_layer, .fraction = inputs.fraction, .surface_organic_carbon_g_c = surface_organic_carbon_g_c, .parameters = inputs.water_heat_parameters });
    try geometry_module.validateDisturbances(owners.soil_geometry, inputs.geometry_changes, inputs.minimum_soil_layer_thickness_m);
    inventory.transferSurfaceFractionToSoil(owners.inventories, .{ .cell = inputs.cell, .destination_soil_layer = inputs.destination_soil_layer, .fraction = inputs.fraction }) catch unreachable;
    chemistry.transferSurfaceFractionToSoil(owners.surface_chemistry, owners.soil_chemistry, inputs.cell, destination, inputs.chemistry_carriers, inputs.dynamic_salts, inputs.fraction) catch unreachable;
    water_heat.transferSurfaceFractionToSoil(owners.water_heat, .{ .cell = inputs.cell, .destination_soil_layer = inputs.destination_soil_layer, .fraction = inputs.fraction, .surface_organic_carbon_g_c = surface_organic_carbon_g_c, .parameters = inputs.water_heat_parameters }) catch unreachable;
    geometry_module.applyDisturbances(owners.soil_geometry, inputs.geometry_changes, inputs.minimum_soil_layer_thickness_m) catch unreachable;
}

test "late surface chemistry failure leaves organic gas and fertilizer inventories unchanged" {
    const std = @import("std");
    const organic_module = @import("soil_organic_initialization.zig");
    const gas_module = @import("gas_transport.zig");
    const surface_fertilizer_module = @import("surface_litter_fertilizer.zig");
    const soil_fertilizer_module = @import("fertilizer_nitrogen_inventory.zig");
    const mineral_fertilizer_module = @import("mineral_fertilizer_inventory.zig");
    var surface_organic = try organic_module.State.init(std.testing.allocator, 1);
    defer surface_organic.deinit();
    var soil_organic = try organic_module.State.init(std.testing.allocator, 1);
    defer soil_organic.deinit();
    var surface_gas = try gas_module.State.init(std.testing.allocator, 1);
    defer surface_gas.deinit();
    var soil_gas = try gas_module.State.init(std.testing.allocator, 1);
    defer soil_gas.deinit();
    var surface_n = try surface_fertilizer_module.State.init(std.testing.allocator, 1);
    defer surface_n.deinit();
    var soil_n = try soil_fertilizer_module.State.init(std.testing.allocator, 1, 1);
    defer soil_n.deinit();
    var mineral = try mineral_fertilizer_module.State.init(std.testing.allocator, 1, 1);
    defer mineral.deinit();
    var surface_chemistry = try surface_chemistry_module.State.init(std.testing.allocator, 1);
    defer surface_chemistry.deinit();
    var soil_chemistry = try soil_chemistry_module.State.init(std.testing.allocator, 1);
    defer soil_chemistry.deinit();
    const config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 4 });
    var grid = try @import("grid.zig").GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var surface_geometry = try @import("surface_litter_geometry_step.zig").State.init(std.testing.allocator, 1);
    defer surface_geometry.deinit();
    var soil_geometry = try geometry_module.State.init(std.testing.allocator, 1, 1);
    defer soil_geometry.deinit();
    try geometry_module.initializeCell(&soil_geometry, 0, 0, &.{0.1}, 0, 1e-6);
    var thermal: @import("soil_thermal.zig").State = undefined;
    var layer_volume = [_]f64{1};
    var dry_capacity = [_]f64{1};
    var total_capacity = [_]f64{1};
    var porosity = [_]f64{0.5};
    thermal.layer_volume_m3 = &layer_volume;
    thermal.dry_solid_heat_capacity_mj_per_m3_k = &dry_capacity;
    thermal.total_heat_capacity_mj_per_m3_k = &total_capacity;
    thermal.porosity_fraction = &porosity;
    var surface_water = [_]f64{0.2};
    var surface_ice = [_]f64{0};
    surface_geometry.expanded_total_volume_m3[0] = 0.2;
    surface_geometry.pore_volume_m3[0] = 0.1;
    surface_geometry.air_volume_m3[0] = 0.05;
    surface_organic.microbial[0].carbon_g_c = 8;
    surface_gas.gaseous_mass_g[0] = 6;
    surface_n.cells[0].ammonium_mol_n = 4;
    surface_chemistry.cells[0].potassium_mol_per_m3 = std.math.nan(f64);
    const inventories: inventory.Owners = .{ .surface_organic = &surface_organic, .soil_organic = &soil_organic, .surface_gas = &surface_gas, .soil_gas = &soil_gas, .surface_nitrogen_fertilizer = &surface_n, .soil_nitrogen_fertilizer = &soil_n, .mineral_fertilizer = &mineral };
    const carriers: chemistry.CarrierVolumes = .{ .surface_water_before_m3 = 1, .soil_shared_water_before_m3 = 1, .soil_phosphate_non_band_water_before_m3 = 1, .surface_water_after_m3 = 0.5, .soil_shared_water_after_m3 = 1.5, .soil_phosphate_non_band_water_after_m3 = 1.5, .surface_dry_mass_before_Mg = 1, .soil_dry_mass_before_Mg = 1, .surface_dry_mass_after_Mg = 0.5, .soil_dry_mass_after_Mg = 1.5 };
    const water_heat_owners: water_heat.Owners = .{ .surface_liquid_water_m3 = &surface_water, .surface_ice_m3 = &surface_ice, .surface_temperature_k = grid.surface_temperature_k, .surface_geometry = &surface_geometry, .grid = &grid, .soil_thermal = &thermal };
    const zero_geometry_change = [_]f64{ 0, 0 };
    try std.testing.expectError(error.InvalidSurfacePondChemistryState, transferSurfaceFractionToSoil(.{ .inventories = inventories, .surface_chemistry = &surface_chemistry, .soil_chemistry = &soil_chemistry, .water_heat = water_heat_owners, .soil_geometry = &soil_geometry }, .{ .cell = 0, .destination_soil_layer = 0, .fraction = 0.5, .chemistry_carriers = carriers, .dynamic_salts = false, .water_heat_parameters = .{ .dry_organic_heat_capacity_mj_per_g_c_k = 2.496e-6, .liquid_water_heat_capacity_mj_per_m3_k = 4.19, .ice_heat_capacity_mj_per_m3_k = 1.9274, .minimum_heat_capacity_mj_per_k = 0 }, .geometry_changes = .{ .pond_m = &zero_geometry_change, .freeze_thaw_m = &zero_geometry_change, .erosion_m = &zero_geometry_change, .organic_carbon_m = &zero_geometry_change }, .minimum_soil_layer_thickness_m = 1e-6 }));
    try std.testing.expectEqual(@as(f64, 8), surface_organic.microbial[0].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 6), surface_gas.gaseous_mass_g[0]);
    try std.testing.expectEqual(@as(f64, 4), surface_n.cells[0].ammonium_mol_n);
    try std.testing.expectEqual(@as(f64, 0), soil_organic.microbial[0].carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0.2), surface_water[0]);
    try std.testing.expectEqual(@as(f64, 0.1), soil_geometry.layer_thickness_m[0]);
}
