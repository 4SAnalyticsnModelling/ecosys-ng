const std = @import("std");
const GridState = @import("grid.zig").GridState;
const SoilHydrology = @import("soil_hydrology.zig").SoilHydrology;
const CellRange = @import("compute.zig").CellRange;
const SoilCatalogEntry = @import("soil_catalog.zig").Entry;

/// Copies one fully resolved soil profile into one horizontal grid cell. The
/// caller controls profile sharing and tile scheduling; this kernel has no file
/// I/O or global state and is directly portable to a future device backend.
pub fn initializeCellHydrology(grid: *GridState, cell_index: usize, hydrology: SoilHydrology) !void {
    if (cell_index >= grid.cell_count) return error.GridIndexOutOfBounds;
    if (hydrology.layer_count == 0) return error.NoSoilLayers;
    if (hydrology.layer_count > grid.soil_layer_capacity) return error.SoilLayerCountMismatch;
    grid.active_soil_layer_count[cell_index] = hydrology.layer_count;
    for (0..hydrology.layer_count) |layer| {
        const grid_index = try grid.layerIndex(cell_index, layer);
        grid.liquid_water_m3[grid_index] = hydrology.matrix_water_volume_m3[layer] + hydrology.macropore_water_volume_m3[layer];
        grid.ice_water_m3[grid_index] = hydrology.matrix_ice_volume_m3[layer] + hydrology.macropore_ice_volume_m3[layer];
        grid.matrix_liquid_water_m3[grid_index] = hydrology.matrix_water_volume_m3[layer];
        grid.macropore_liquid_water_m3[grid_index] = hydrology.macropore_water_volume_m3[layer];
        grid.matrix_ice_water_m3[grid_index] = hydrology.matrix_ice_volume_m3[layer];
        grid.macropore_ice_water_m3[grid_index] = hydrology.macropore_ice_volume_m3[layer];
        grid.matrix_pore_capacity_m3[grid_index] = hydrology.matrix_pore_volume_m3[layer];
        grid.macropore_pore_capacity_m3[grid_index] = hydrology.macropore_volume_m3[layer];
        grid.matrix_air_volume_m3[grid_index] = hydrology.matrix_air_volume_m3[layer];
        grid.macropore_air_volume_m3[grid_index] = hydrology.macropore_air_volume_m3[layer];
        grid.air_volume_m3[grid_index] = hydrology.air_volume_m3[layer];
    }
    try grid.validateFinite();
}

pub const UniformHydrologyContext = struct {
    grid: *GridState,
    hydrology: SoilHydrology,
    reference_cell_area_m2: f64,
    horizontal_cell_width_m: []const f64,
    vertical_cell_width_m: []const f64,
};

/// Tile kernel used when a landscape unit shares one resolved soil profile.
/// More complex domains invoke the same kernel once per profile-backed unit.
pub fn initializeUniformHydrologyTile(context: *UniformHydrologyContext, range: CellRange) !void {
    if (!std.math.isFinite(context.reference_cell_area_m2) or context.reference_cell_area_m2 <= 0.0) return error.InvalidReferenceCellArea;
    if (context.horizontal_cell_width_m.len != context.grid.cell_count or
        context.vertical_cell_width_m.len != context.grid.cell_count)
        return error.CellGeometryDimensionMismatch;
    for (range.first..range.end) |cell_index| {
        const cell_area_m2 = context.horizontal_cell_width_m[cell_index] *
            context.vertical_cell_width_m[cell_index];
        const volume_scale = cell_area_m2 / context.reference_cell_area_m2;
        try initializeCellHydrologyScaled(context.grid, cell_index, context.hydrology, volume_scale);
    }
}

fn initializeCellHydrologyScaled(grid: *GridState, cell_index: usize, hydrology: SoilHydrology, volume_scale: f64) !void {
    if (!std.math.isFinite(volume_scale) or volume_scale <= 0.0) return error.InvalidCellVolumeScale;
    if (cell_index >= grid.cell_count) return error.GridIndexOutOfBounds;
    if (hydrology.layer_count > grid.soil_layer_capacity) return error.SoilLayerCountMismatch;
    grid.active_soil_layer_count[cell_index] = hydrology.layer_count;
    for (0..hydrology.layer_count) |layer| {
        const grid_index = try grid.layerIndex(cell_index, layer);
        grid.liquid_water_m3[grid_index] = volume_scale * (hydrology.matrix_water_volume_m3[layer] + hydrology.macropore_water_volume_m3[layer]);
        grid.ice_water_m3[grid_index] = volume_scale * (hydrology.matrix_ice_volume_m3[layer] + hydrology.macropore_ice_volume_m3[layer]);
        grid.matrix_liquid_water_m3[grid_index] = volume_scale * hydrology.matrix_water_volume_m3[layer];
        grid.macropore_liquid_water_m3[grid_index] = volume_scale * hydrology.macropore_water_volume_m3[layer];
        grid.matrix_ice_water_m3[grid_index] = volume_scale * hydrology.matrix_ice_volume_m3[layer];
        grid.macropore_ice_water_m3[grid_index] = volume_scale * hydrology.macropore_ice_volume_m3[layer];
        grid.matrix_pore_capacity_m3[grid_index] = volume_scale * hydrology.matrix_pore_volume_m3[layer];
        grid.macropore_pore_capacity_m3[grid_index] = volume_scale * hydrology.macropore_volume_m3[layer];
        grid.matrix_air_volume_m3[grid_index] = volume_scale * hydrology.matrix_air_volume_m3[layer];
        grid.macropore_air_volume_m3[grid_index] = volume_scale * hydrology.macropore_air_volume_m3[layer];
        grid.air_volume_m3[grid_index] = volume_scale * hydrology.air_volume_m3[layer];
    }
}

pub const MappedHydrologyContext = struct {
    grid: *GridState,
    catalog_entries: []const SoilCatalogEntry,
    catalog_index_by_cell: []const usize,
    horizontal_cell_width_m: []const f64,
    vertical_cell_width_m: []const f64,
};

pub fn initializeMappedHydrologyTile(context: *MappedHydrologyContext, range: CellRange) !void {
    if (context.horizontal_cell_width_m.len != context.grid.cell_count or
        context.vertical_cell_width_m.len != context.grid.cell_count)
        return error.CellGeometryDimensionMismatch;
    for (range.first..range.end) |cell_index| {
        if (cell_index >= context.catalog_index_by_cell.len) return error.SoilCatalogMapOutOfBounds;
        const catalog_index = context.catalog_index_by_cell[cell_index];
        if (catalog_index >= context.catalog_entries.len) return error.SoilCatalogMapOutOfBounds;
        const cell_area_m2 = context.horizontal_cell_width_m[cell_index] *
            context.vertical_cell_width_m[cell_index];
        try initializeCellHydrologyScaled(context.grid, cell_index, context.catalog_entries[catalog_index].hydrology_per_m2, cell_area_m2);
    }
}

test "resolved hydrology populates heap-backed grid state" {
    const allocator = std.testing.allocator;
    const source = try @import("test_fixtures.zig").soilProfileSource(allocator, @typeInfo(@import("soil_profile.zig").LayerProperty).@"enum".fields.len);
    defer allocator.free(source);
    var profile = try @import("soil_profile.zig").parsePhysicalProfile(allocator, source);
    defer profile.deinit();
    var material = try @import("soil_initialization.zig").SoilMaterial.init(allocator, profile, @import("soil_profile_derivation.zig").compatibilityParameters());
    defer material.deinit();
    var hydrology = try SoilHydrology.init(allocator, profile, material, 1.0, @import("soil_water_retention.zig").compatibilityParameters());
    defer hydrology.deinit();
    const config = try @import("config.zig").SimulationConfig.init(
        .{ .lon_count = 2, .lat_count = 1, .soil_layers = profile.total_layer_count, .plant_populations = 8 },
        .{ .worker_threads = 1, .tile_cells = 2 },
        .{ .relative_tolerance = 1.0e-8, .absolute_tolerance = 1.0e-11, .max_nonlinear_iterations = 40 },
    );
    var grid = try GridState.init(allocator, config);
    defer grid.deinit();
    try initializeCellHydrology(&grid, 1, hydrology);
    const first_layer = try grid.layerIndex(1, 0);
    try std.testing.expectApproxEqAbs(
        hydrology.matrix_water_volume_m3[0] + hydrology.macropore_water_volume_m3[0],
        grid.liquid_water_m3[first_layer],
        1.0e-12,
    );
    try std.testing.expectApproxEqAbs(hydrology.matrix_water_volume_m3[0], grid.matrix_liquid_water_m3[first_layer], 1.0e-12);
    try std.testing.expectApproxEqAbs(hydrology.macropore_water_volume_m3[0], grid.macropore_liquid_water_m3[first_layer], 1.0e-12);
    try std.testing.expectApproxEqAbs(hydrology.air_volume_m3[0], grid.air_volume_m3[first_layer], 1.0e-12);
}
