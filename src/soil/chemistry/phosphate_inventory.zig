const std = @import("std");
const Grid = @import("../../state/grid.zig");
const SoilChemistry = @import("../solute/chemistry_state.zig");
const SurfaceChemistry = @import("../../surface/litter_chemistry.zig");
const Phosphate = @import("../solute/phosphate_network.zig");

/// REDIST UPP4 inventory. Chemistry stores mineral concentrations by water
/// volume, so each fertilizer zone is weighted by its runtime water fraction
/// before applying the exact 1/1/1/3/2 phosphorus stoichiometry.
pub fn precipitatedPhosphorus_g_p(
    grid: *const Grid.GridState,
    soil: *const SoilChemistry.State,
    surface: *const SurfaceChemistry.State,
    soil_water_m3: []const f64,
    litter_water_m3: []const f64,
    cell: usize,
    non_band_fraction: f64,
    band_fraction: f64,
    phosphorus_molar_mass_g_per_mol: f64,
) !f64 {
    if (cell >= grid.cell_count or soil.cell_count != grid.layer_count or surface.cells.len != grid.cell_count or soil_water_m3.len != grid.layer_count or litter_water_m3.len != grid.cell_count) return error.PhosphateInventoryDimensionMismatch;
    inline for (.{ non_band_fraction, band_fraction, phosphorus_molar_mass_g_per_mol }) |value| if (!std.math.isFinite(value)) return error.InvalidPhosphateInventory;
    if (non_band_fraction < 0 or band_fraction < 0 or @abs(non_band_fraction + band_fraction - 1) > 1e-12 or phosphorus_molar_mass_g_per_mol <= 0) return error.InvalidPhosphateInventory;
    var phosphorus_mol: f64 = 0;
    const surface_water = litter_water_m3[cell];
    if (!std.math.isFinite(surface_water) or surface_water < 0) return error.InvalidPhosphateInventory;
    const surface_minerals = surface.cells[cell].phosphate_minerals;
    phosphorus_mol += surface_water * (surface_minerals.aluminum_phosphate_mol_per_m3 + surface_minerals.iron_phosphate_mol_per_m3 + surface_minerals.dicalcium_phosphate_mol_per_m3 + 3 * surface_minerals.hydroxyapatite_mol_per_m3 + 2 * surface_minerals.monocalcium_phosphate_mol_per_m3);
    for (0..grid.active_soil_layer_count[cell]) |local_layer| {
        const layer = try grid.layerIndex(cell, local_layer);
        const water = soil_water_m3[layer];
        if (!std.math.isFinite(water) or water < 0) return error.InvalidPhosphateInventory;
        phosphorus_mol += water * (non_band_fraction * solidPhosphorus_mol_per_m3(soil.non_band_phosphate[layer]) + band_fraction * solidPhosphorus_mol_per_m3(soil.band_phosphate[layer]));
    }
    const result = phosphorus_mol * phosphorus_molar_mass_g_per_mol;
    if (!std.math.isFinite(result) or result < 0) return error.InvalidPhosphateInventory;
    return result;
}

fn solidPhosphorus_mol_per_m3(state: Phosphate.State) f64 {
    return state.aluminum_phosphate_solid_mol_per_m3 +
        state.iron_phosphate_solid_mol_per_m3 +
        state.dicalcium_phosphate_solid_mol_per_m3 +
        3 * state.hydroxyapatite_solid_mol_per_m3 +
        2 * state.monocalcium_phosphate_solid_mol_per_m3;
}

/// REDIST UPX4 inventory: adsorbed HPO4 plus H2PO4 on litter and soil
/// phosphate surfaces. Surface and soil concentrations are both mol P/Mg,
/// but use their distinct runtime dry-mass owners.
pub fn exchangeablePhosphorus_g_p(
    grid: *const Grid.GridState,
    soil: *const SoilChemistry.State,
    surface: *const SurfaceChemistry.State,
    matrix_bulk_volume_m3: []const f64,
    bulk_density_megagrams_per_m3: []const f64,
    litter_dry_mass_megagrams: []const f64,
    cell: usize,
    non_band_fraction: f64,
    band_fraction: f64,
    phosphorus_molar_mass_g_per_mol: f64,
) !f64 {
    if (cell >= grid.cell_count or soil.cell_count != grid.layer_count or surface.cells.len != grid.cell_count or matrix_bulk_volume_m3.len != grid.layer_count or bulk_density_megagrams_per_m3.len != grid.layer_count or litter_dry_mass_megagrams.len != grid.cell_count) return error.PhosphateInventoryDimensionMismatch;
    inline for (.{ non_band_fraction, band_fraction, phosphorus_molar_mass_g_per_mol }) |value| if (!std.math.isFinite(value)) return error.InvalidPhosphateInventory;
    if (non_band_fraction < 0 or band_fraction < 0 or @abs(non_band_fraction + band_fraction - 1) > 1e-12 or phosphorus_molar_mass_g_per_mol <= 0) return error.InvalidPhosphateInventory;
    const litter_mass = litter_dry_mass_megagrams[cell];
    if (!std.math.isFinite(litter_mass) or litter_mass < 0) return error.InvalidPhosphateInventory;
    const litter_sites = surface.cells[cell].phosphate_surface;
    var phosphorus_mol = litter_mass * (litter_sites.adsorbed_hpo4_mol_p_per_megagram + litter_sites.adsorbed_h2po4_mol_p_per_megagram);
    for (0..grid.active_soil_layer_count[cell]) |local_layer| {
        const layer = try grid.layerIndex(cell, local_layer);
        const volume = matrix_bulk_volume_m3[layer];
        const density = bulk_density_megagrams_per_m3[layer];
        if (!std.math.isFinite(volume) or volume < 0 or !std.math.isFinite(density) or density < 0) return error.InvalidPhosphateInventory;
        const non_band = soil.non_band_phosphate[layer];
        const band = soil.band_phosphate[layer];
        phosphorus_mol += volume * density * (non_band_fraction * (non_band.adsorbed_hpo4_mol_p_per_megagram + non_band.adsorbed_h2po4_mol_p_per_megagram) + band_fraction * (band.adsorbed_hpo4_mol_p_per_megagram + band.adsorbed_h2po4_mol_p_per_megagram));
    }
    const result = phosphorus_mol * phosphorus_molar_mass_g_per_mol;
    if (!std.math.isFinite(result) or result < 0) return error.InvalidPhosphateInventory;
    return result;
}

/// REDIST UPO4: soluble HPO4 plus H2PO4 in surface litter and both runtime
/// soil fertilizer zones.
pub fn solublePhosphorus_g_p(
    grid: *const Grid.GridState,
    soil: *const SoilChemistry.State,
    surface: *const SurfaceChemistry.State,
    soil_water_m3: []const f64,
    litter_water_m3: []const f64,
    cell: usize,
    non_band_fraction: f64,
    band_fraction: f64,
    phosphorus_molar_mass_g_per_mol: f64,
) !f64 {
    if (cell >= grid.cell_count or soil.cell_count != grid.layer_count or surface.cells.len != grid.cell_count or soil_water_m3.len != grid.layer_count or litter_water_m3.len != grid.cell_count) return error.PhosphateInventoryDimensionMismatch;
    inline for (.{ non_band_fraction, band_fraction, phosphorus_molar_mass_g_per_mol }) |value| if (!std.math.isFinite(value)) return error.InvalidPhosphateInventory;
    if (non_band_fraction < 0 or band_fraction < 0 or @abs(non_band_fraction + band_fraction - 1) > 1e-12 or phosphorus_molar_mass_g_per_mol <= 0) return error.InvalidPhosphateInventory;
    const litter_water = litter_water_m3[cell];
    if (!std.math.isFinite(litter_water) or litter_water < 0) return error.InvalidPhosphateInventory;
    const litter = surface.cells[cell];
    var phosphorus_mol = litter_water * (litter.hpo4_mol_p_per_m3 + litter.h2po4_mol_p_per_m3);
    for (0..grid.active_soil_layer_count[cell]) |local_layer| {
        const layer = try grid.layerIndex(cell, local_layer);
        const water = soil_water_m3[layer];
        if (!std.math.isFinite(water) or water < 0) return error.InvalidPhosphateInventory;
        const non_band = soil.non_band_phosphate[layer];
        const band = soil.band_phosphate[layer];
        phosphorus_mol += water * (non_band_fraction * (non_band.dissolved_hpo4_mol_p_per_m3 + non_band.dissolved_h2po4_mol_p_per_m3) + band_fraction * (band.dissolved_hpo4_mol_p_per_m3 + band.dissolved_h2po4_mol_p_per_m3));
    }
    const result = phosphorus_mol * phosphorus_molar_mass_g_per_mol;
    if (!std.math.isFinite(result) or result < 0) return error.InvalidPhosphateInventory;
    return result;
}

test "UPP4 applies exact phosphate mineral stoichiometry across runtime zones and surface" {
    const config = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 1, .lat_count = 1, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try Grid.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var soil = try SoilChemistry.State.init(std.testing.allocator, grid.layer_count);
    defer soil.deinit();
    var surface = try SurfaceChemistry.State.init(std.testing.allocator, grid.cell_count);
    defer surface.deinit();
    soil.non_band_phosphate[0].aluminum_phosphate_solid_mol_per_m3 = 2;
    soil.band_phosphate[0].hydroxyapatite_solid_mol_per_m3 = 4;
    soil.non_band_phosphate[1].monocalcium_phosphate_solid_mol_per_m3 = 3;
    surface.cells[0].phosphate_minerals.iron_phosphate_mol_per_m3 = 5;
    const result = try precipitatedPhosphorus_g_p(&grid, &soil, &surface, &.{ 10, 20 }, &.{2}, 0, 0.75, 0.25, 31);
    // surface 10 mol P; layer 0: 10*(1.5+3); layer 1: 20*4.5 mol P
    try std.testing.expectEqual(@as(f64, 145 * 31), result);
}

test "UPX4 uses distinct litter and soil masses with runtime phosphate zones" {
    const config = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 1, .lat_count = 1, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try Grid.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var soil = try SoilChemistry.State.init(std.testing.allocator, grid.layer_count);
    defer soil.deinit();
    var surface = try SurfaceChemistry.State.init(std.testing.allocator, grid.cell_count);
    defer surface.deinit();
    surface.cells[0].phosphate_surface.adsorbed_hpo4_mol_p_per_megagram = 2;
    soil.non_band_phosphate[0].adsorbed_h2po4_mol_p_per_megagram = 4;
    soil.band_phosphate[1].adsorbed_hpo4_mol_p_per_megagram = 8;
    const result = try exchangeablePhosphorus_g_p(&grid, &soil, &surface, &.{ 10, 20 }, &.{ 2, 3 }, &.{5}, 0, 0.75, 0.25, 31);
    // litter 10; layer 0 60; layer 1 120 = 190 mol P.
    try std.testing.expectEqual(@as(f64, 190 * 31), result);
}

test "UPO4 sums surface and runtime-zone soluble phosphate inventories" {
    const config = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 1, .lat_count = 1, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try Grid.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var soil = try SoilChemistry.State.init(std.testing.allocator, grid.layer_count);
    defer soil.deinit();
    var surface = try SurfaceChemistry.State.init(std.testing.allocator, grid.cell_count);
    defer surface.deinit();
    surface.cells[0].h2po4_mol_p_per_m3 = 2;
    soil.non_band_phosphate[0].dissolved_hpo4_mol_p_per_m3 = 4;
    soil.band_phosphate[1].dissolved_h2po4_mol_p_per_m3 = 8;
    const result = try solublePhosphorus_g_p(&grid, &soil, &surface, &.{ 10, 20 }, &.{5}, 0, 0.75, 0.25, 31);
    try std.testing.expectEqual(@as(f64, 80 * 31), result);
}
