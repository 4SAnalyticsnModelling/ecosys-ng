const std = @import("std");
const terrain_module = @import("../state/terrain_hydrology.zig");
const spatial_grid = @import("../state/spatial_grid.zig");
const lateral_store = @import("../state/tile_lateral_contribution_store.zig");

pub const Parameters = struct {
    ground_surface_retention_m3_per_m2: f64,
    runoff_roughness_h_per_m_one_third: f64,
    maximum_hydraulic_volume_m3: f64 = 1.0e-3,
    manning_time_conversion_s_per_h: f64 = 3.6e3,
    negligible_water_m3: f64 = 1.0e-12,
};

pub const BoundaryFractions = struct {
    north: []const f64,
    east: []const f64,
    south: []const f64,
    west: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    excess_surface_water_m3: []f64,
    excess_surface_ice_m3: []f64,
    runoff_velocity_m_per_s: []f64,
    total_runoff_m3_per_step: []f64,
    east_runoff_m3_per_step: []f64,
    west_runoff_m3_per_step: []f64,
    south_runoff_m3_per_step: []f64,
    north_runoff_m3_per_step: []f64,
    water_change_m3: []f64,
    exported_water_m3: []f64,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize) !State {
        if (cell_count == 0) return error.ZeroSurfaceRunoffCellCount;
        var result: State = undefined;
        result.allocator = allocator;
        result.cell_count = cell_count;
        var allocated: usize = 0;
        errdefer inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64 and allocated > 0) {
            allocated -= 1;
            allocator.free(@field(result, field.name));
        };
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            @field(result, field.name) = try allocator.alloc(f64, cell_count);
            @memset(@field(result, field.name), 0);
            allocated += 1;
        };
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

/// WATSUB `XVOLW/XVOLI`, Manning `QRM/QRV`, and directional `QRMN`.
/// Every face is generated from the converged hourly surface state and then
/// reduced once; no full ecosystem sub-hour cycle is repeated.
pub fn route(
    state: *State,
    columns: usize,
    rows: usize,
    terrain: *const terrain_module.State,
    cell_area_m2: []const f64,
    surface_water_m3: []f64,
    surface_ice_m3: []const f64,
    litter_retention_capacity_m3: []const f64,
    lateral_connection_mode_by_cell: []const u8,
    boundary_fractions: BoundaryFractions,
    parameters: Parameters,
) !void {
    try calculateFluxes(
        state,
        columns,
        rows,
        terrain,
        cell_area_m2,
        surface_water_m3,
        surface_ice_m3,
        litter_retention_capacity_m3,
        lateral_connection_mode_by_cell,
        boundary_fractions,
        parameters,
    );
    try commitWaterChanges(surface_water_m3, state.water_change_m3);
}

/// Calculates the immutable-source runoff snapshot and signed changes without
/// mutating surface water. This phase can be followed either by the resident
/// commit or by the two-pass Morton contribution transaction.
pub fn calculateFluxes(
    state: *State,
    columns: usize,
    rows: usize,
    terrain: *const terrain_module.State,
    cell_area_m2: []const f64,
    surface_water_m3: []f64,
    surface_ice_m3: []const f64,
    litter_retention_capacity_m3: []const f64,
    lateral_connection_mode_by_cell: []const u8,
    boundary_fractions: BoundaryFractions,
    parameters: Parameters,
) !void {
    const count = try std.math.mul(usize, columns, rows);
    if (count != state.cell_count or terrain.columns != columns or terrain.rows != rows or cell_area_m2.len != count or surface_water_m3.len != count or surface_ice_m3.len != count or litter_retention_capacity_m3.len != count or lateral_connection_mode_by_cell.len != count or boundary_fractions.north.len != count or boundary_fractions.east.len != count or boundary_fractions.south.len != count or boundary_fractions.west.len != count) return error.SurfaceRunoffDimensionMismatch;
    for (lateral_connection_mode_by_cell) |mode| if (mode != 1 and mode != 3) return error.InvalidLateralConnectionMode;
    try validateParameters(parameters, boundary_fractions);
    inline for (.{ state.excess_surface_water_m3, state.excess_surface_ice_m3, state.runoff_velocity_m_per_s, state.total_runoff_m3_per_step, state.east_runoff_m3_per_step, state.west_runoff_m3_per_step, state.south_runoff_m3_per_step, state.north_runoff_m3_per_step, state.water_change_m3, state.exported_water_m3 }) |values| @memset(values, 0);

    for (0..count) |cell| {
        const water = surface_water_m3[cell];
        const ice = surface_ice_m3[cell];
        const retention = litter_retention_capacity_m3[cell];
        const area = cell_area_m2[cell];
        inline for (.{ water, ice, retention, area }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSurfaceRunoffState;
        if (area <= 0) return error.InvalidSurfaceRunoffState;
        const total = water + ice;
        if (total <= parameters.negligible_water_m3) continue;
        const retained_water = water / total * retention;
        const retained_ice = ice / total * retention;
        state.excess_surface_water_m3[cell] = @max(0, water - retained_water);
        state.excess_surface_ice_m3[cell] = @max(0, ice - retained_ice);
        const excess_total = state.excess_surface_water_m3[cell] + state.excess_surface_ice_m3[cell];
        const ground_retention = parameters.ground_surface_retention_m3_per_m2 * area;
        if (excess_total <= ground_retention or state.excess_surface_water_m3[cell] <= parameters.negligible_water_m3) continue;
        const hydraulic_volume_m3 = @min(
            parameters.maximum_hydraulic_volume_m3,
            (excess_total - ground_retention) * state.excess_surface_water_m3[cell] / excess_total,
        );
        const hydraulic_depth_m = hydraulic_volume_m3 / area;
        const velocity_m_per_s = std.math.pow(f64, hydraulic_depth_m, 0.67) * @sqrt(terrain.slope_m_per_m[cell]) / parameters.runoff_roughness_h_per_m_one_third;
        state.runoff_velocity_m_per_s[cell] = velocity_m_per_s;
        state.total_runoff_m3_per_step[cell] = @min(hydraulic_volume_m3, velocity_m_per_s * hydraulic_depth_m * terrain.flow_width_m[cell] * parameters.manning_time_conversion_s_per_h);
    }

    for (0..rows) |row| for (0..columns) |column| {
        const source = row * columns + column;
        const available = state.total_runoff_m3_per_step[source];
        if (available <= parameters.negligible_water_m3) continue;
        const source_surface_m = terrain.current_surface_elevation_m[source] +
            (state.excess_surface_water_m3[source] + state.excess_surface_ice_m3[source]) / cell_area_m2[source];
        if (terrain.runoff_to_east[source]) state.east_runoff_m3_per_step[source] = try directionalFlux(source, if (column + 1 < columns) source + 1 else null, source_surface_m, available, terrain.east_west_runoff_fraction[source], boundary_fractions.east[source], lateral_connection_mode_by_cell, terrain, cell_area_m2, state);
        if (terrain.runoff_to_west[source]) state.west_runoff_m3_per_step[source] = try directionalFlux(source, if (column > 0) source - 1 else null, source_surface_m, available, terrain.east_west_runoff_fraction[source], boundary_fractions.west[source], lateral_connection_mode_by_cell, terrain, cell_area_m2, state);
        if (terrain.runoff_to_south[source]) state.south_runoff_m3_per_step[source] = try directionalFlux(source, if (row + 1 < rows) source + columns else null, source_surface_m, available, terrain.north_south_runoff_fraction[source], boundary_fractions.south[source], lateral_connection_mode_by_cell, terrain, cell_area_m2, state);
        if (terrain.runoff_to_north[source]) state.north_runoff_m3_per_step[source] = try directionalFlux(source, if (row > 0) source - columns else null, source_surface_m, available, terrain.north_south_runoff_fraction[source], boundary_fractions.north[source], lateral_connection_mode_by_cell, terrain, cell_area_m2, state);
    };

    for (0..rows) |row| for (0..columns) |column| {
        const source = row * columns + column;
        const fluxes = [_]f64{ state.east_runoff_m3_per_step[source], state.west_runoff_m3_per_step[source], state.south_runoff_m3_per_step[source], state.north_runoff_m3_per_step[source] };
        const outgoing = fluxes[0] + fluxes[1] + fluxes[2] + fluxes[3];
        if (outgoing > surface_water_m3[source] + 64 * std.math.floatEps(f64) * @max(1, surface_water_m3[source])) return error.SurfaceRunoffExceedsAvailableWater;
        state.water_change_m3[source] -= outgoing;
        if (column + 1 < columns) state.water_change_m3[source + 1] += fluxes[0] else state.exported_water_m3[source] += fluxes[0];
        if (column > 0) state.water_change_m3[source - 1] += fluxes[1] else state.exported_water_m3[source] += fluxes[1];
        if (row + 1 < rows) state.water_change_m3[source + columns] += fluxes[2] else state.exported_water_m3[source] += fluxes[2];
        if (row > 0) state.water_change_m3[source - columns] += fluxes[3] else state.exported_water_m3[source] += fluxes[3];
    };
    try validateWaterChanges(surface_water_m3, state.water_change_m3);
}

pub fn commitWaterChanges(
    surface_water_m3: []f64,
    water_change_m3: []const f64,
) !void {
    try validateWaterChanges(surface_water_m3, water_change_m3);
    for (surface_water_m3, water_change_m3) |*water, change|
        water.* = @max(0, water.* + change);
}

fn validateWaterChanges(
    surface_water_m3: []const f64,
    water_change_m3: []const f64,
) !void {
    if (surface_water_m3.len != water_change_m3.len)
        return error.SurfaceRunoffDimensionMismatch;
    for (surface_water_m3, water_change_m3) |water, change|
        if (!std.math.isFinite(water) or !std.math.isFinite(change) or
            water + change < -1e-12)
            return error.InvalidSurfaceRunoffCandidate;
}

pub const lateral_component_count: usize = 2;
pub const water_change_component: usize = 0;
pub const boundary_export_component: usize = 1;

/// Converts the immutable directional runoff snapshot for one source tile
/// into signed endpoint records. Each source cell belongs to exactly one
/// tile, so outgoing water is emitted once; neighboring receipts are carried
/// by the durable contribution sidecar until the destination tile commits.
pub fn appendOwnedTileContributions(
    allocator: std.mem.Allocator,
    plan: spatial_grid.TilePlan,
    tile_index: usize,
    state: *const State,
    contributions: *std.ArrayList(lateral_store.Contribution),
) !void {
    if (state.cell_count != plan.lat_count * plan.lon_count)
        return error.SurfaceRunoffDimensionMismatch;
    const owned_cells = try plan.ownedCells(tile_index);
    for (owned_cells) |source| {
        const row = source / plan.lon_count;
        const column = source % plan.lon_count;
        const fluxes = [_]f64{
            state.east_runoff_m3_per_step[source],
            state.west_runoff_m3_per_step[source],
            state.south_runoff_m3_per_step[source],
            state.north_runoff_m3_per_step[source],
        };
        var outgoing_m3: f64 = 0;
        for (fluxes) |flux_m3| {
            if (!std.math.isFinite(flux_m3) or flux_m3 < 0)
                return error.InvalidSurfaceRunoffFlux;
            outgoing_m3 += flux_m3;
        }
        if (outgoing_m3 > 0) try contributions.append(allocator, .{
            .target_cell = source,
            .component = water_change_component,
            .delta = -outgoing_m3,
        });
        const destinations = [_]?usize{
            if (column + 1 < plan.lon_count) source + 1 else null,
            if (column > 0) source - 1 else null,
            if (row + 1 < plan.lat_count)
                source + plan.lon_count
            else
                null,
            if (row > 0) source - plan.lon_count else null,
        };
        for (fluxes, destinations) |flux_m3, destination| {
            if (flux_m3 == 0) continue;
            try contributions.append(allocator, .{
                .target_cell = destination orelse source,
                .component = if (destination == null)
                    boundary_export_component
                else
                    water_change_component,
                .delta = flux_m3,
            });
        }
    }
}

/// Commits only the destination tile's gathered runoff components.
pub fn commitOwnedTileContributions(
    plan: spatial_grid.TilePlan,
    tile_index: usize,
    state: *State,
    surface_water_m3: []f64,
    gathered: []const f64,
) !void {
    if (surface_water_m3.len != state.cell_count or
        gathered.len != state.cell_count * lateral_component_count)
        return error.SurfaceRunoffDimensionMismatch;
    const owned_cells = try plan.ownedCells(tile_index);
    for (owned_cells) |cell| {
        const water_change_m3 =
            gathered[cell * lateral_component_count + water_change_component];
        const exported_water_m3 =
            gathered[cell * lateral_component_count + boundary_export_component];
        if (!std.math.isFinite(water_change_m3) or
            !std.math.isFinite(exported_water_m3) or
            exported_water_m3 < 0 or
            surface_water_m3[cell] + water_change_m3 < -1e-12)
            return error.InvalidSurfaceRunoffCandidate;
    }
    for (owned_cells) |cell| {
        const water_change_m3 =
            gathered[cell * lateral_component_count + water_change_component];
        const exported_water_m3 =
            gathered[cell * lateral_component_count + boundary_export_component];
        state.water_change_m3[cell] = water_change_m3;
        state.exported_water_m3[cell] = exported_water_m3;
        surface_water_m3[cell] =
            @max(0, surface_water_m3[cell] + water_change_m3);
    }
}

fn directionalFlux(source: usize, destination: ?usize, source_surface_m: f64, available_m3: f64, direction_fraction: f64, boundary_fraction: f64, lateral_connection_mode_by_cell: []const u8, terrain: *const terrain_module.State, cell_area_m2: []const f64, state: *const State) !f64 {
    if (destination) |target| {
        if (lateral_connection_mode_by_cell[source] != 1 or
            lateral_connection_mode_by_cell[target] != 1) return 0;
        const destination_surface_m = terrain.current_surface_elevation_m[target] +
            (state.excess_surface_water_m3[target] + state.excess_surface_ice_m3[target]) / cell_area_m2[target];
        if (source_surface_m <= destination_surface_m) return 0;
        const equilibrium_m3 = @max(0, 0.5 * (source_surface_m - destination_surface_m) * cell_area_m2[source] * cell_area_m2[target] / (cell_area_m2[source] + cell_area_m2[target]));
        return @min(equilibrium_m3, available_m3) * direction_fraction;
    }
    return available_m3 * direction_fraction * boundary_fraction;
}

fn validateParameters(parameters: Parameters, boundaries: BoundaryFractions) !void {
    inline for (@typeInfo(Parameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters, field.name)) or @field(parameters, field.name) < 0) return error.InvalidSurfaceRunoffParameter;
    if (parameters.runoff_roughness_h_per_m_one_third <= 0 or parameters.maximum_hydraulic_volume_m3 <= 0 or parameters.manning_time_conversion_s_per_h <= 0) return error.InvalidSurfaceRunoffParameter;
    inline for (.{ boundaries.north, boundaries.east, boundaries.south, boundaries.west }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidSurfaceRunoffBoundary;
}

test "internal Manning runoff conserves water on a runtime grid" {
    const topography = @import("../state/topography.zig");
    var units = [_]topography.LandscapeUnit{.{ .west_column = 1, .north_row = 1, .east_column = 2, .south_row = 1, .compass_aspect_degrees = 90, .geometric_aspect_degrees = 90, .slope_degrees = 5, .unused_slope_input = 0, .initial_snowpack_depth_m = 0, .soil_profile_file = "soil" }};
    var terrain = try terrain_module.State.initMapped(std.testing.allocator, .{ .allocator = std.testing.allocator, .units = &units }, &.{ 0, 0 }, &.{ 1, 1 }, &.{ 1, 1 }, 2, 1);
    defer terrain.deinit();
    // Ensure an unambiguous eastward surface gradient for this carrier test.
    terrain.runoff_to_east[0] = true;
    terrain.runoff_to_west[0] = false;
    terrain.east_west_runoff_fraction[0] = 1;
    terrain.north_south_runoff_fraction[0] = 0;
    try terrain.bindInitialSurfaceElevations(&.{ 1, 0 });
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var water = [_]f64{ 0.02, 0.01 };
    const initial = water[0] + water[1];
    try route(&state, 2, 1, &terrain, &.{ 1, 1 }, &water, &.{ 0, 0 }, &.{ 0.005, 0.005 }, &.{ 1, 1 }, .{ .north = &.{ 0, 0 }, .east = &.{ 0, 0 }, .south = &.{ 0, 0 }, .west = &.{ 0, 0 } }, .{ .ground_surface_retention_m3_per_m2 = 0, .runoff_roughness_h_per_m_one_third = 0.1 });
    try std.testing.expect(state.east_runoff_m3_per_step[0] > 0);
    try std.testing.expectApproxEqAbs(initial, water[0] + water[1], 1e-15);
    water = .{ 0.02, 0.01 };
    try route(&state, 2, 1, &terrain, &.{ 1, 1 }, &water, &.{ 0, 0 }, &.{ 0.005, 0.005 }, &.{ 3, 1 }, .{ .north = &.{ 0, 0 }, .east = &.{ 0, 0 }, .south = &.{ 0, 0 }, .west = &.{ 0, 0 } }, .{ .ground_surface_retention_m3_per_m2 = 0, .runoff_roughness_h_per_m_one_third = 0.1 });
    try std.testing.expectEqual(@as(f64, 0), state.east_runoff_m3_per_step[0]);
    try std.testing.expectEqual(@as(f64, 0.02), water[0]);
    try std.testing.expectEqual(@as(f64, 0.01), water[1]);
}

test "open runoff boundary records exported water" {
    const topography = @import("../state/topography.zig");
    var units = [_]topography.LandscapeUnit{.{ .west_column = 1, .north_row = 1, .east_column = 1, .south_row = 1, .compass_aspect_degrees = 90, .geometric_aspect_degrees = 90, .slope_degrees = 5, .unused_slope_input = 0, .initial_snowpack_depth_m = 0, .soil_profile_file = "soil" }};
    var terrain = try terrain_module.State.initMapped(std.testing.allocator, .{ .allocator = std.testing.allocator, .units = &units }, &.{0}, &.{1}, &.{1}, 1, 1);
    defer terrain.deinit();
    terrain.runoff_to_east[0] = true;
    terrain.east_west_runoff_fraction[0] = 1;
    terrain.north_south_runoff_fraction[0] = 0;
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var water = [_]f64{0.02};
    try route(&state, 1, 1, &terrain, &.{1}, &water, &.{0}, &.{0.005}, &.{3}, .{ .north = &.{0}, .east = &.{0.5}, .south = &.{0}, .west = &.{0} }, .{ .ground_surface_retention_m3_per_m2 = 0, .runoff_roughness_h_per_m_one_third = 0.1 });
    try std.testing.expect(state.exported_water_m3[0] > 0);
    try std.testing.expectApproxEqAbs(@as(f64, 0.02), water[0] + state.exported_water_m3[0], 1e-15);
}

test "Morton two-pass runoff matches resident conservative commit across tiles" {
    const topography = @import("../state/topography.zig");
    var units = [_]topography.LandscapeUnit{.{
        .west_column = 1,
        .north_row = 1,
        .east_column = 4,
        .south_row = 1,
        .compass_aspect_degrees = 90,
        .geometric_aspect_degrees = 90,
        .slope_degrees = 5,
        .unused_slope_input = 0,
        .initial_snowpack_depth_m = 0,
        .soil_profile_file = "soil",
    }};
    var terrain = try terrain_module.State.initMapped(
        std.testing.allocator,
        .{ .allocator = std.testing.allocator, .units = &units },
        &.{ 0, 0, 0, 0 },
        &.{ 1, 1, 1, 1 },
        &.{ 1, 1, 1, 1 },
        4,
        1,
    );
    defer terrain.deinit();
    for (0..4) |cell| {
        terrain.runoff_to_east[cell] = cell + 1 < 4;
        terrain.runoff_to_west[cell] = false;
        terrain.runoff_to_north[cell] = false;
        terrain.runoff_to_south[cell] = false;
        terrain.east_west_runoff_fraction[cell] = 1;
        terrain.north_south_runoff_fraction[cell] = 0;
        terrain.initial_surface_elevation_m[cell] =
            @as(f64, @floatFromInt(4 - cell));
        terrain.current_surface_elevation_m[cell] =
            terrain.initial_surface_elevation_m[cell];
    }
    var state = try State.init(std.testing.allocator, 4);
    defer state.deinit();
    const initial_water = [_]f64{ 0.04, 0.03, 0.02, 0.01 };
    var resident_water = initial_water;
    try calculateFluxes(
        &state,
        4,
        1,
        &terrain,
        &.{ 1, 1, 1, 1 },
        &resident_water,
        &.{ 0, 0, 0, 0 },
        &.{ 0, 0, 0, 0 },
        &.{ 1, 1, 1, 1 },
        .{
            .north = &.{ 0, 0, 0, 0 },
            .east = &.{ 0, 0, 0, 0 },
            .south = &.{ 0, 0, 0, 0 },
            .west = &.{ 0, 0, 0, 0 },
        },
        .{
            .ground_surface_retention_m3_per_m2 = 0,
            .runoff_roughness_h_per_m_one_third = 0.1,
        },
    );
    try commitWaterChanges(&resident_water, state.water_change_m3);

    var plan = try spatial_grid.TilePlan.init(
        std.testing.allocator,
        1,
        4,
        1,
        2,
        2,
    );
    defer plan.deinit();
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const store = try lateral_store.FileStore.init(
        std.testing.allocator,
        std.testing.io,
        temporary.dir,
        4096,
        30,
        31,
    );
    for (plan.tiles, 0..) |_, tile_index| {
        var contributions =
            std.ArrayList(lateral_store.Contribution).empty;
        defer contributions.deinit(std.testing.allocator);
        try appendOwnedTileContributions(
            std.testing.allocator,
            plan,
            tile_index,
            &state,
            &contributions,
        );
        try store.saveSourceTile(
            plan,
            tile_index,
            lateral_component_count,
            contributions.items,
        );
    }
    try store.publish(plan);
    var tiled_water = initial_water;
    @memset(state.water_change_m3, 0);
    @memset(state.exported_water_m3, 0);
    const gathered = try std.testing.allocator.alloc(
        f64,
        4 * lateral_component_count,
    );
    defer std.testing.allocator.free(gathered);
    @memset(gathered, 0);
    for (plan.tiles, 0..) |_, tile_index| {
        try store.gatherOwnedTile(
            plan,
            tile_index,
            lateral_component_count,
            gathered,
        );
        try commitOwnedTileContributions(
            plan,
            tile_index,
            &state,
            &tiled_water,
            gathered,
        );
    }
    try std.testing.expectEqualSlices(f64, &resident_water, &tiled_water);
    var initial_total_m3: f64 = 0;
    var tiled_total_m3: f64 = 0;
    for (initial_water) |water_m3| initial_total_m3 += water_m3;
    for (tiled_water) |water_m3| tiled_total_m3 += water_m3;
    try std.testing.expectApproxEqAbs(
        initial_total_m3,
        tiled_total_m3,
        1e-15,
    );
}
