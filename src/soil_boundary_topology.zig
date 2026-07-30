const std = @import("std");
const grid_module = @import("grid.zig");
const site_module = @import("site.zig");
const terrain_module = @import("terrain_hydrology.zig");

pub const Direction = enum { north, east, south, west, lower };

pub const Face = struct {
    cell_index: usize,
    layer_index: usize,
    direction: Direction,
    direction_sign: f64,
    directional_layer_width_m: f64,
    slope_sine: f64,
    natural_water_table_distance_m: f64,
    natural_exchange_fraction: f64,
    artificial_water_table_distance_m: f64,
    artificial_exchange_fraction: f64,
    surface_runoff_fraction: f64,
    is_lower_boundary: bool,
};

/// Heap-owned perimeter and profile-bottom topology. Its size follows the
/// runtime grid and each cell's active layer count; no historical JX/JY/JZ
/// ceiling survives in this representation.
pub const State = struct {
    allocator: std.mem.Allocator,
    faces: []Face,
    water_table_mode: []u8,
    natural_water_table_depth_m: []f64,
    internal_water_table_depth_m: []f64,
    active_layer_depth_m: []f64,
    artificial_water_table_depth_m: []f64,
    natural_water_table_surface_slope: []f64,
    artificial_water_table_surface_slope: []f64,

    pub fn initMapped(
        allocator: std.mem.Allocator,
        grid: *const grid_module.GridState,
        terrain: *const terrain_module.State,
        columns: usize,
        rows: usize,
        horizontal_cell_width_m: []const f64,
        vertical_cell_width_m: []const f64,
        site_by_cell: []const site_module.Site,
    ) !State {
        if (columns == 0 or rows == 0 or
            try std.math.mul(usize, columns, rows) != grid.cell_count or
            terrain.columns != columns or terrain.rows != rows or
            horizontal_cell_width_m.len != grid.cell_count or
            vertical_cell_width_m.len != grid.cell_count or
            site_by_cell.len != grid.cell_count)
            return error.SoilBoundaryTopologyDimensionMismatch;
        var face_count: usize = grid.cell_count;
        for (0..grid.cell_count) |cell| {
            const row = cell / columns;
            const column = cell % columns;
            var lateral_edges: usize = 0;
            if (row == 0) lateral_edges += 1;
            if (column + 1 == columns) lateral_edges += 1;
            if (row + 1 == rows) lateral_edges += 1;
            if (column == 0) lateral_edges += 1;
            face_count = try std.math.add(usize, face_count, try std.math.mul(usize, lateral_edges, grid.active_soil_layer_count[cell]));
        }
        const faces = try allocator.alloc(Face, face_count);
        errdefer allocator.free(faces);
        const natural_depth = try allocator.alloc(f64, grid.cell_count);
        errdefer allocator.free(natural_depth);
        const artificial_depth = try allocator.alloc(f64, grid.cell_count);
        errdefer allocator.free(artificial_depth);
        const internal_depth = try allocator.alloc(f64, grid.cell_count);
        errdefer allocator.free(internal_depth);
        const active_layer_depth = try allocator.alloc(f64, grid.cell_count);
        errdefer allocator.free(active_layer_depth);
        const water_table_mode = try allocator.alloc(u8, grid.cell_count);
        errdefer allocator.free(water_table_mode);
        const natural_slope = try allocator.alloc(f64, grid.cell_count);
        errdefer allocator.free(natural_slope);
        const artificial_slope = try allocator.alloc(f64, grid.cell_count);
        errdefer allocator.free(artificial_slope);
        for (0..grid.cell_count) |cell| {
            const site = site_by_cell[cell];
            const elevation_adjustment_m = terrain.minimum_surface_elevation_m - terrain.relative_surface_elevation_m[cell];
            natural_depth[cell] = site.initial_water_table_depth_m - elevation_adjustment_m * (1 - site.natural_water_table_surface_slope);
            internal_depth[cell] = natural_depth[cell];
            active_layer_depth[cell] = 9999;
            artificial_depth[cell] = if (site.artificial_water_table_depth_m) |depth| @max(0, depth - elevation_adjustment_m * (1 - site.artificial_water_table_surface_slope.?)) else 0;
            water_table_mode[cell] = site.water_table_mode;
            natural_slope[cell] = site.natural_water_table_surface_slope;
            artificial_slope[cell] = site.artificial_water_table_surface_slope orelse 0;
        }
        var next: usize = 0;
        for (0..grid.cell_count) |cell| {
            const site = &site_by_cell[cell];
            const row = cell / columns;
            const column = cell % columns;
            const active_layers = grid.active_soil_layer_count[cell];
            if (active_layers == 0 or active_layers > grid.soil_layer_capacity) return error.InvalidActiveSoilLayerCount;
            for (0..active_layers) |layer| {
                const layer_index = try grid.layerIndex(cell, layer);
                const x_width_m = horizontal_cell_width_m[cell];
                const y_width_m = vertical_cell_width_m[cell];
                if (row == 0) appendLateral(faces, &next, cell, layer_index, .north, 0, y_width_m, terrain.north_south_slope_m_per_m[cell], site);
                if (column + 1 == columns) appendLateral(faces, &next, cell, layer_index, .east, 1, x_width_m, terrain.east_west_slope_m_per_m[cell], site);
                if (row + 1 == rows) appendLateral(faces, &next, cell, layer_index, .south, 2, y_width_m, terrain.north_south_slope_m_per_m[cell], site);
                if (column == 0) appendLateral(faces, &next, cell, layer_index, .west, 3, x_width_m, terrain.east_west_slope_m_per_m[cell], site);
            }
            const bottom_layer_index = try grid.layerIndex(cell, active_layers - 1);
            faces[next] = .{
                .cell_index = cell,
                .layer_index = bottom_layer_index,
                .direction = .lower,
                .direction_sign = -1,
                .directional_layer_width_m = 1,
                .slope_sine = 1,
                .natural_water_table_distance_m = 1,
                .natural_exchange_fraction = site.lower_boundary_exchange_fraction,
                .artificial_water_table_distance_m = 0,
                .artificial_exchange_fraction = 0,
                .surface_runoff_fraction = 0,
                .is_lower_boundary = true,
            };
            next += 1;
        }
        std.debug.assert(next == faces.len);
        return .{ .allocator = allocator, .faces = faces, .water_table_mode = water_table_mode, .natural_water_table_depth_m = natural_depth, .internal_water_table_depth_m = internal_depth, .active_layer_depth_m = active_layer_depth, .artificial_water_table_depth_m = artificial_depth, .natural_water_table_surface_slope = natural_slope, .artificial_water_table_surface_slope = artificial_slope };
    }

    /// HOUR1 DPTHT refresh. A saturated zone must remain continuous downward
    /// until the prescribed natural table is reached; its upper boundary is
    /// interpolated from THETS (air-entry water content) in the layer above.
    pub fn refreshInternalWaterTable(self: *State, grid: *const grid_module.GridState, matrix_bulk_volume_m3: []const f64, porosity_fraction: []const f64, air_entry_water_fraction: []const f64, layer_thickness_m: []const f64, layer_midpoint_depth_m: []const f64, layer_bottom_depth_m: []const f64, minimum_air_filled_porosity: f64, minimum_frozen_pore_fraction: f64) !void {
        const layers = grid.layer_count;
        if (matrix_bulk_volume_m3.len != layers or porosity_fraction.len != layers or air_entry_water_fraction.len != layers or layer_thickness_m.len != layers or layer_midpoint_depth_m.len != layers or layer_bottom_depth_m.len != layers or !std.math.isFinite(minimum_air_filled_porosity) or minimum_air_filled_porosity < 0 or !std.math.isFinite(minimum_frozen_pore_fraction) or minimum_frozen_pore_fraction < 0) return error.SoilBoundaryTopologyDimensionMismatch;
        for (0..grid.cell_count) |cell| {
            const active_layers = grid.active_soil_layer_count[cell];
            self.active_layer_depth_m[cell] = 9999;
            for (0..active_layers) |local_layer| {
                const layer = try grid.layerIndex(cell, local_layer);
                const pore_volume_m3 = grid.matrix_pore_capacity_m3[layer] + grid.macropore_pore_capacity_m3[layer];
                const ice_volume_m3 = grid.matrix_ice_water_m3[layer] + grid.macropore_ice_water_m3[layer];
                if (pore_volume_m3 <= 0 or ice_volume_m3 < minimum_frozen_pore_fraction * pore_volume_m3) continue;
                var lower_layers_frozen = true;
                var lower_local = local_layer + 1;
                while (lower_local < active_layers) : (lower_local += 1) {
                    const lower = try grid.layerIndex(cell, lower_local);
                    const lower_pore_m3 = grid.matrix_pore_capacity_m3[lower] + grid.macropore_pore_capacity_m3[lower];
                    const lower_ice_m3 = grid.matrix_ice_water_m3[lower] + grid.macropore_ice_water_m3[lower];
                    if (lower_pore_m3 > 0 and lower_ice_m3 < minimum_frozen_pore_fraction * lower_pore_m3) {
                        lower_layers_frozen = false;
                        break;
                    }
                }
                if (!lower_layers_frozen) continue;
                self.active_layer_depth_m[cell] = layer_bottom_depth_m[layer] - layer_thickness_m[layer] * std.math.clamp(ice_volume_m3 / pore_volume_m3, 0, 1);
                break;
            }
            if (self.water_table_mode[cell] == 0) continue;
            var found = false;
            for (0..active_layers) |local_layer| {
                const layer = try grid.layerIndex(cell, local_layer);
                const total_pore_volume_m3 = grid.matrix_pore_capacity_m3[layer] + grid.macropore_pore_capacity_m3[layer];
                const air_fraction = if (matrix_bulk_volume_m3[layer] > 0) grid.air_volume_m3[layer] / matrix_bulk_volume_m3[layer] else 0;
                if (total_pore_volume_m3 <= 0 or (air_fraction >= minimum_air_filled_porosity and local_layer + 1 != active_layers)) continue;
                var continuous = true;
                if (layer_midpoint_depth_m[layer] < self.natural_water_table_depth_m[cell]) {
                    var lower_local = local_layer + 1;
                    while (lower_local < active_layers) : (lower_local += 1) {
                        const lower = try grid.layerIndex(cell, lower_local);
                        const lower_air_fraction = if (matrix_bulk_volume_m3[lower] > 0) grid.air_volume_m3[lower] / matrix_bulk_volume_m3[lower] else 0;
                        if (lower_air_fraction >= minimum_air_filled_porosity and lower_local + 1 != active_layers) {
                            continuous = false;
                            break;
                        }
                        if (layer_midpoint_depth_m[lower] >= self.natural_water_table_depth_m[cell]) break;
                    }
                }
                if (!continuous) continue;
                if (local_layer == 0) {
                    self.internal_water_table_depth_m[cell] = 0;
                } else {
                    const above = try grid.layerIndex(cell, local_layer - 1);
                    const denominator = porosity_fraction[above] - air_entry_water_fraction[above];
                    const water_fraction = grid.matrix_liquid_water_m3[above] / matrix_bulk_volume_m3[above];
                    const saturated_fraction = if (denominator > 0) std.math.clamp((water_fraction - air_entry_water_fraction[above]) / denominator, 0, 1) else 0;
                    self.internal_water_table_depth_m[cell] = layer_bottom_depth_m[above] - layer_thickness_m[above] * saturated_fraction;
                }
                found = true;
                break;
            }
            if (!found) self.internal_water_table_depth_m[cell] = self.natural_water_table_depth_m[cell];
            if (!std.math.isFinite(self.internal_water_table_depth_m[cell])) return error.NonFiniteInternalWaterTableDepth;
        }
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.faces);
        self.allocator.free(self.artificial_water_table_surface_slope);
        self.allocator.free(self.natural_water_table_surface_slope);
        self.allocator.free(self.water_table_mode);
        self.allocator.free(self.artificial_water_table_depth_m);
        self.allocator.free(self.internal_water_table_depth_m);
        self.allocator.free(self.active_layer_depth_m);
        self.allocator.free(self.natural_water_table_depth_m);
        self.* = undefined;
    }
};

fn appendLateral(faces: []Face, next: *usize, cell_index: usize, layer_index: usize, direction: Direction, site_index: usize, directional_layer_width_m: f64, slope_sine: f64, site: *const site_module.Site) void {
    faces[next.*] = .{
        .cell_index = cell_index,
        .layer_index = layer_index,
        .direction = direction,
        // Exact WATSUB XN convention: east/south=-1, west/north=+1.
        .direction_sign = if (direction == .east or direction == .south) -1 else 1,
        .directional_layer_width_m = directional_layer_width_m,
        .slope_sine = slope_sine,
        .natural_water_table_distance_m = site.natural_water_table_distance_m[site_index],
        .natural_exchange_fraction = site.natural_subsurface_exchange_fraction[site_index],
        .artificial_water_table_distance_m = site.artificial_water_table_distance_m[site_index],
        .artificial_exchange_fraction = site.artificial_subsurface_exchange_fraction[site_index],
        .surface_runoff_fraction = site.surface_runoff_boundary_fraction[site_index],
        .is_lower_boundary = false,
    };
    next.* += 1;
}

test "runtime boundary topology preserves READI compass order and WATSUB signs" {
    const cfg = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 2, .grid_rows = 1, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    var site = try site_module.parse(std.testing.allocator, @import("test_fixtures.zig").site_source_2x1, 2, 1);
    defer site.deinit();
    const topography_source = "1 1 2 1 0 0 0 0\nsoil\n";
    var topography = try @import("topography.zig").parse(std.testing.allocator, topography_source);
    defer topography.deinit();
    var terrain = try terrain_module.State.initMapped(std.testing.allocator, topography, &.{ 0, 0 }, &.{ 1, 1 }, &.{ 1, 1 }, 2, 1);
    defer terrain.deinit();
    var sites = [_]site_module.Site{ site, site };
    sites[1].initial_water_table_depth_m = 4.0;
    sites[1].water_table_mode = 1;
    sites[1].natural_water_table_distance_m[1] = 25;
    sites[1].natural_subsurface_exchange_fraction[1] = 0.75;
    var state = try State.initMapped(std.testing.allocator, &grid, &terrain, 2, 1, &.{ 1, 1 }, &.{ 1, 1 }, &sites);
    defer state.deinit();
    // Two cells x two layers: N/S on both cells, W on first, E on second,
    // plus one bottom face per cell.
    try std.testing.expectEqual(@as(usize, 14), state.faces.len);
    const east_expected_depth = sites[1].initial_water_table_depth_m -
        (terrain.minimum_surface_elevation_m - terrain.relative_surface_elevation_m[1]) *
            (1 - sites[1].natural_water_table_surface_slope);
    try std.testing.expectEqual(east_expected_depth, state.natural_water_table_depth_m[1]);
    try std.testing.expectEqual(@as(u8, 1), state.water_table_mode[1]);
    try std.testing.expectEqual(Direction.north, state.faces[0].direction);
    try std.testing.expectEqual(@as(f64, 1), state.faces[0].direction_sign);
    try std.testing.expectEqual(@as(f64, 10), state.faces[0].natural_water_table_distance_m);
    var east_found = false;
    var lower_count: usize = 0;
    for (state.faces) |face| {
        if (face.direction == .east) {
            east_found = true;
            try std.testing.expectEqual(@as(f64, -1), face.direction_sign);
            try std.testing.expectEqual(@as(f64, 25), face.natural_water_table_distance_m);
            try std.testing.expectEqual(@as(f64, 0.75), face.natural_exchange_fraction);
        }
        if (face.is_lower_boundary) {
            lower_count += 1;
            try std.testing.expectEqual(site.lower_boundary_exchange_fraction, face.natural_exchange_fraction);
        }
    }
    try std.testing.expect(east_found);
    try std.testing.expectEqual(@as(usize, 2), lower_count);
    const expected_depth = site.initial_water_table_depth_m - (terrain.minimum_surface_elevation_m - terrain.relative_surface_elevation_m[0]) * (1 - site.natural_water_table_surface_slope);
    try std.testing.expectApproxEqAbs(expected_depth, state.natural_water_table_depth_m[0], 1e-12);
}

test "internal water table refresh interpolates HOUR1 DPTHT from runtime saturation" {
    const cfg = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    var site = try site_module.parse(std.testing.allocator, @import("test_fixtures.zig").site_source, 1, 1);
    defer site.deinit();
    var topography = try @import("topography.zig").parse(std.testing.allocator, "1 1 1 1 0 0 0 0\nsoil\n");
    defer topography.deinit();
    var terrain = try terrain_module.State.initMapped(std.testing.allocator, topography, &.{0}, &.{1}, &.{1}, 1, 1);
    defer terrain.deinit();
    var state = try State.initMapped(std.testing.allocator, &grid, &terrain, 1, 1, &.{1}, &.{1}, &.{site});
    defer state.deinit();
    state.water_table_mode[0] = 1;
    grid.matrix_pore_capacity_m3[0] = 0.5;
    grid.matrix_pore_capacity_m3[1] = 0.5;
    grid.matrix_liquid_water_m3[0] = 0.4;
    grid.matrix_liquid_water_m3[1] = 0.5;
    grid.matrix_air_volume_m3[0] = 0.1;
    grid.matrix_air_volume_m3[1] = 0;
    grid.air_volume_m3[0] = 0.1;
    grid.air_volume_m3[1] = 0;
    try state.refreshInternalWaterTable(&grid, &.{ 1, 1 }, &.{ 0.5, 0.5 }, &.{ 0.3, 0.3 }, &.{ 0.1, 0.1 }, &.{ 0.05, 0.15 }, &.{ 0.1, 0.2 }, 1e-3, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), state.internal_water_table_depth_m[0], 1e-12);
    grid.matrix_ice_water_m3[1] = 0.25;
    try state.refreshInternalWaterTable(&grid, &.{ 1, 1 }, &.{ 0.5, 0.5 }, &.{ 0.3, 0.3 }, &.{ 0.1, 0.1 }, &.{ 0.05, 0.15 }, &.{ 0.1, 0.2 }, 1e-3, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 0.15), state.active_layer_depth_m[0], 1e-12);
}
