const std = @import("std");
const grid_module = @import("grid.zig");
const snow_module = @import("snow_solute_transport.zig");
const solute = @import("solute_transport.zig");

pub const State = struct {
    allocator: std.mem.Allocator,
    columns: usize,
    rows: usize,
    soil_layer_capacity: usize,
    snow_layer_capacity: usize,
    micropore_water_volume_m3: []f64,
    macropore_water_volume_m3: []f64,
    matrix_air_volume_m3: []f64,
    macropore_air_volume_m3: []f64,
    air_volume_m3: []f64,
    water_vapor_volume_m3: []f64,
    /// Soil-cell-major +x,+y,+z water flux; positive moves toward + direction.
    micropore_face_flux_m3_per_step: []f64,
    macropore_face_flux_m3_per_step: []f64,
    vapor_face_flux_m3_per_step: []f64,
    heat_face_flux_megajoules_per_step: []f64,
    /// Positive values leave the modeled domain; negative values recharge it.
    micropore_external_water_flux_m3_per_step: []f64,
    macropore_external_water_flux_m3_per_step: []f64,
    /// Accepted whole-hour REDIST `DVOLI` owners. Positive values are freezing.
    matrix_ice_volume_change_m3_per_step: []f64,
    macropore_ice_volume_change_m3_per_step: []f64,
    total_ice_volume_change_m3_per_step: []f64,
    /// Positive artificial-drain discharge, aggregated by horizontal cell.
    artificial_drainage_outflow_m3_per_step: []f64,
    runoff_total_m3_per_step: []f64,
    runoff_east_m3_per_step: []f64,
    runoff_west_m3_per_step: []f64,
    runoff_south_m3_per_step: []f64,
    runoff_north_m3_per_step: []f64,
    snow_transfer_total_m3_per_step: []f64,
    snow_transfer_east_m3_per_step: []f64,
    snow_transfer_west_m3_per_step: []f64,
    snow_transfer_south_m3_per_step: []f64,
    snow_transfer_north_m3_per_step: []f64,
    /// Total snow+liquid+ice carrier volume (`VOLSL(1)`) at the surface.
    snow_surface_carrier_volume_m3: []f64,
    snow_liquid_water_volume_m3: []f64,
    snow_downward_water_flux_m3_per_step: []f64,
    snow_to_litter_water_flux_m3_per_step: []f64,
    snow_to_soil_micropore_flux_m3_per_step: []f64,
    snow_to_soil_macropore_flux_m3_per_step: []f64,

    pub fn init(allocator: std.mem.Allocator, columns: usize, rows: usize, soil_layers: usize, snow_layers: usize) !State {
        if (columns == 0 or rows == 0 or soil_layers == 0 or snow_layers == 0) return error.ZeroTransportHydrologyDimension;
        const cells = try std.math.mul(usize, columns, rows);
        const soil_components = try std.math.mul(usize, cells, soil_layers);
        const snow_components = try std.math.mul(usize, cells, snow_layers);
        var result: State = undefined;
        result.allocator = allocator;
        result.columns = columns;
        result.rows = rows;
        result.soil_layer_capacity = soil_layers;
        result.snow_layer_capacity = snow_layers;
        result.micropore_water_volume_m3 = try zeroes(allocator, soil_components);
        errdefer allocator.free(result.micropore_water_volume_m3);
        result.macropore_water_volume_m3 = try zeroes(allocator, soil_components);
        errdefer allocator.free(result.macropore_water_volume_m3);
        result.matrix_air_volume_m3 = try zeroes(allocator, soil_components);
        errdefer allocator.free(result.matrix_air_volume_m3);
        result.macropore_air_volume_m3 = try zeroes(allocator, soil_components);
        errdefer allocator.free(result.macropore_air_volume_m3);
        result.air_volume_m3 = try zeroes(allocator, soil_components);
        errdefer allocator.free(result.air_volume_m3);
        result.water_vapor_volume_m3 = try zeroes(allocator, soil_components);
        errdefer allocator.free(result.water_vapor_volume_m3);
        result.micropore_face_flux_m3_per_step = try zeroes(allocator, try std.math.mul(usize, soil_components, 3));
        errdefer allocator.free(result.micropore_face_flux_m3_per_step);
        result.macropore_face_flux_m3_per_step = try zeroes(allocator, try std.math.mul(usize, soil_components, 3));
        errdefer allocator.free(result.macropore_face_flux_m3_per_step);
        result.vapor_face_flux_m3_per_step = try zeroes(allocator, try std.math.mul(usize, soil_components, 3));
        errdefer allocator.free(result.vapor_face_flux_m3_per_step);
        result.heat_face_flux_megajoules_per_step = try zeroes(allocator, try std.math.mul(usize, soil_components, 3));
        errdefer allocator.free(result.heat_face_flux_megajoules_per_step);
        result.micropore_external_water_flux_m3_per_step = try zeroes(allocator, soil_components);
        errdefer allocator.free(result.micropore_external_water_flux_m3_per_step);
        result.macropore_external_water_flux_m3_per_step = try zeroes(allocator, soil_components);
        errdefer allocator.free(result.macropore_external_water_flux_m3_per_step);
        result.matrix_ice_volume_change_m3_per_step = try zeroes(allocator, soil_components);
        errdefer allocator.free(result.matrix_ice_volume_change_m3_per_step);
        result.macropore_ice_volume_change_m3_per_step = try zeroes(allocator, soil_components);
        errdefer allocator.free(result.macropore_ice_volume_change_m3_per_step);
        result.total_ice_volume_change_m3_per_step = try zeroes(allocator, soil_components);
        errdefer allocator.free(result.total_ice_volume_change_m3_per_step);
        result.artificial_drainage_outflow_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.artificial_drainage_outflow_m3_per_step);
        result.runoff_total_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.runoff_total_m3_per_step);
        result.runoff_east_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.runoff_east_m3_per_step);
        result.runoff_west_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.runoff_west_m3_per_step);
        result.runoff_south_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.runoff_south_m3_per_step);
        result.runoff_north_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.runoff_north_m3_per_step);
        result.snow_transfer_total_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.snow_transfer_total_m3_per_step);
        result.snow_transfer_east_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.snow_transfer_east_m3_per_step);
        result.snow_transfer_west_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.snow_transfer_west_m3_per_step);
        result.snow_transfer_south_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.snow_transfer_south_m3_per_step);
        result.snow_transfer_north_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.snow_transfer_north_m3_per_step);
        result.snow_surface_carrier_volume_m3 = try zeroes(allocator, cells);
        errdefer allocator.free(result.snow_surface_carrier_volume_m3);
        result.snow_liquid_water_volume_m3 = try zeroes(allocator, snow_components);
        errdefer allocator.free(result.snow_liquid_water_volume_m3);
        result.snow_downward_water_flux_m3_per_step = try zeroes(allocator, snow_components);
        errdefer allocator.free(result.snow_downward_water_flux_m3_per_step);
        result.snow_to_litter_water_flux_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.snow_to_litter_water_flux_m3_per_step);
        result.snow_to_soil_micropore_flux_m3_per_step = try zeroes(allocator, cells);
        errdefer allocator.free(result.snow_to_soil_micropore_flux_m3_per_step);
        result.snow_to_soil_macropore_flux_m3_per_step = try zeroes(allocator, cells);
        return result;
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.snow_to_soil_macropore_flux_m3_per_step);
        self.allocator.free(self.snow_to_soil_micropore_flux_m3_per_step);
        self.allocator.free(self.snow_to_litter_water_flux_m3_per_step);
        self.allocator.free(self.snow_downward_water_flux_m3_per_step);
        self.allocator.free(self.snow_liquid_water_volume_m3);
        self.allocator.free(self.snow_surface_carrier_volume_m3);
        self.allocator.free(self.snow_transfer_north_m3_per_step);
        self.allocator.free(self.snow_transfer_south_m3_per_step);
        self.allocator.free(self.snow_transfer_west_m3_per_step);
        self.allocator.free(self.snow_transfer_east_m3_per_step);
        self.allocator.free(self.snow_transfer_total_m3_per_step);
        self.allocator.free(self.runoff_north_m3_per_step);
        self.allocator.free(self.runoff_south_m3_per_step);
        self.allocator.free(self.runoff_west_m3_per_step);
        self.allocator.free(self.runoff_east_m3_per_step);
        self.allocator.free(self.runoff_total_m3_per_step);
        self.allocator.free(self.artificial_drainage_outflow_m3_per_step);
        self.allocator.free(self.total_ice_volume_change_m3_per_step);
        self.allocator.free(self.macropore_ice_volume_change_m3_per_step);
        self.allocator.free(self.matrix_ice_volume_change_m3_per_step);
        self.allocator.free(self.macropore_external_water_flux_m3_per_step);
        self.allocator.free(self.micropore_external_water_flux_m3_per_step);
        self.allocator.free(self.macropore_face_flux_m3_per_step);
        self.allocator.free(self.heat_face_flux_megajoules_per_step);
        self.allocator.free(self.vapor_face_flux_m3_per_step);
        self.allocator.free(self.micropore_face_flux_m3_per_step);
        self.allocator.free(self.air_volume_m3);
        self.allocator.free(self.water_vapor_volume_m3);
        self.allocator.free(self.macropore_air_volume_m3);
        self.allocator.free(self.matrix_air_volume_m3);
        self.allocator.free(self.macropore_water_volume_m3);
        self.allocator.free(self.micropore_water_volume_m3);
        self.* = undefined;
    }

    pub fn syncStorage(self: *State, grid: *const grid_module.GridState, snow: *const snow_module.State) !void {
        if (grid.cell_count != self.columns * self.rows or grid.soil_layer_capacity != self.soil_layer_capacity or snow.cell_count != grid.cell_count or snow.layer_capacity != self.snow_layer_capacity) return error.TransportHydrologyDimensionMismatch;
        @memcpy(self.micropore_water_volume_m3, grid.matrix_liquid_water_m3);
        @memcpy(self.macropore_water_volume_m3, grid.macropore_liquid_water_m3);
        @memcpy(self.matrix_air_volume_m3, grid.matrix_air_volume_m3);
        @memcpy(self.macropore_air_volume_m3, grid.macropore_air_volume_m3);
        @memcpy(self.air_volume_m3, grid.air_volume_m3);
        @memcpy(self.water_vapor_volume_m3, grid.water_vapor_volume_m3);
        @memcpy(self.snow_liquid_water_volume_m3, snow.liquid_water_volume_m3);
        try self.validateFinite();
    }

    pub fn validateFinite(self: *const State) !void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            for (@field(self, field.name), 0..) |value, index| {
                if (!std.math.isFinite(value)) {
                    std.log.err("non-finite transport hydrology: field={s} index={d} value={e}", .{ field.name, index, value });
                    return error.NonFiniteTransportHydrology;
                }
                const signed_flux_field = comptime std.mem.eql(u8, field.name, "micropore_face_flux_m3_per_step") or std.mem.eql(u8, field.name, "macropore_face_flux_m3_per_step") or std.mem.eql(u8, field.name, "vapor_face_flux_m3_per_step") or std.mem.eql(u8, field.name, "heat_face_flux_megajoules_per_step") or std.mem.eql(u8, field.name, "micropore_external_water_flux_m3_per_step") or std.mem.eql(u8, field.name, "macropore_external_water_flux_m3_per_step") or std.mem.endsWith(u8, field.name, "ice_volume_change_m3_per_step");
                if (!signed_flux_field and value < 0) {
                    std.log.err("negative transport hydrology: field={s} index={d} value={e}", .{ field.name, index, value });
                    return error.NegativeTransportHydrologyValue;
                }
            }
        };
    }

    pub fn soilIndex(self: *const State, column: usize, row: usize, layer: usize) !usize {
        if (column >= self.columns or row >= self.rows or layer >= self.soil_layer_capacity) return error.TransportHydrologyIndexOutOfBounds;
        return (row * self.columns + column) * self.soil_layer_capacity + layer;
    }
};

pub const SoilFaces = struct {
    allocator: std.mem.Allocator,
    /// 0=x, 1=y, 2=z; parallel to both pore-domain face arrays.
    direction_axis: []u2,
    micropore_faces: []solute.Face,
    macropore_faces: []solute.Face,
    micropore_water_flux_m3_per_step: []f64,
    macropore_water_flux_m3_per_step: []f64,
    vapor_flux_m3_per_step: []f64,
    heat_flux_megajoules_per_step: []f64,

    pub fn deinit(self: *SoilFaces) void {
        self.allocator.free(self.macropore_water_flux_m3_per_step);
        self.allocator.free(self.heat_flux_megajoules_per_step);
        self.allocator.free(self.vapor_flux_m3_per_step);
        self.allocator.free(self.micropore_water_flux_m3_per_step);
        self.allocator.free(self.macropore_faces);
        self.allocator.free(self.micropore_faces);
        self.allocator.free(self.direction_axis);
        self.* = undefined;
    }
};

/// Builds the active runtime +x/+y/+z topology shared by water and solutes.
pub fn buildSoilFaces(allocator: std.mem.Allocator, hydrology: *const State, grid: *const grid_module.GridState) !SoilFaces {
    const modes = try allocator.alloc(u8, grid.cell_count);
    defer allocator.free(modes);
    @memset(modes, 1);
    return buildSoilFacesMapped(allocator, hydrology, grid, modes);
}

pub fn buildSoilFacesMapped(allocator: std.mem.Allocator, hydrology: *const State, grid: *const grid_module.GridState, lateral_connection_mode_by_cell: []const u8) !SoilFaces {
    if (grid.cell_count != hydrology.columns * hydrology.rows or grid.soil_layer_capacity != hydrology.soil_layer_capacity or lateral_connection_mode_by_cell.len != grid.cell_count) return error.TransportHydrologyDimensionMismatch;
    for (lateral_connection_mode_by_cell) |mode| if (mode != 1 and mode != 3) return error.InvalidLateralConnectionMode;
    var count: usize = 0;
    for (0..hydrology.rows) |row| for (0..hydrology.columns) |column| {
        const horizontal = row * hydrology.columns + column;
        for (0..grid.active_soil_layer_count[horizontal]) |layer| {
            if (column + 1 < hydrology.columns and lateral_connection_mode_by_cell[horizontal] == 1 and lateral_connection_mode_by_cell[horizontal + 1] == 1 and layer < grid.active_soil_layer_count[horizontal + 1]) count = try std.math.add(usize, count, 1);
            if (row + 1 < hydrology.rows and lateral_connection_mode_by_cell[horizontal] == 1 and lateral_connection_mode_by_cell[horizontal + hydrology.columns] == 1 and layer < grid.active_soil_layer_count[horizontal + hydrology.columns]) count = try std.math.add(usize, count, 1);
            if (layer + 1 < grid.active_soil_layer_count[horizontal]) count = try std.math.add(usize, count, 1);
        }
    };
    const micropore_faces = try allocator.alloc(solute.Face, count);
    errdefer allocator.free(micropore_faces);
    const direction_axis = try allocator.alloc(u2, count);
    errdefer allocator.free(direction_axis);
    const macropore_faces = try allocator.alloc(solute.Face, count);
    errdefer allocator.free(macropore_faces);
    const micropore = try allocator.alloc(f64, count);
    errdefer allocator.free(micropore);
    const macropore = try allocator.alloc(f64, count);
    errdefer allocator.free(macropore);
    const vapor = try allocator.alloc(f64, count);
    errdefer allocator.free(vapor);
    const heat = try allocator.alloc(f64, count);
    errdefer allocator.free(heat);
    var face_index: usize = 0;
    for (0..hydrology.rows) |row| for (0..hydrology.columns) |column| {
        const horizontal = row * hydrology.columns + column;
        for (0..grid.active_soil_layer_count[horizontal]) |layer| {
            const first = try hydrology.soilIndex(column, row, layer);
            if (column + 1 < hydrology.columns and lateral_connection_mode_by_cell[horizontal] == 1 and lateral_connection_mode_by_cell[horizontal + 1] == 1 and layer < grid.active_soil_layer_count[horizontal + 1]) appendFace(hydrology, first, try hydrology.soilIndex(column + 1, row, layer), 0, direction_axis, micropore_faces, macropore_faces, micropore, macropore, vapor, heat, &face_index);
            if (row + 1 < hydrology.rows and lateral_connection_mode_by_cell[horizontal] == 1 and lateral_connection_mode_by_cell[horizontal + hydrology.columns] == 1 and layer < grid.active_soil_layer_count[horizontal + hydrology.columns]) appendFace(hydrology, first, try hydrology.soilIndex(column, row + 1, layer), 1, direction_axis, micropore_faces, macropore_faces, micropore, macropore, vapor, heat, &face_index);
            if (layer + 1 < grid.active_soil_layer_count[horizontal]) appendFace(hydrology, first, try hydrology.soilIndex(column, row, layer + 1), 2, direction_axis, micropore_faces, macropore_faces, micropore, macropore, vapor, heat, &face_index);
        }
    };
    std.debug.assert(face_index == count);
    return .{ .allocator = allocator, .direction_axis = direction_axis, .micropore_faces = micropore_faces, .macropore_faces = macropore_faces, .micropore_water_flux_m3_per_step = micropore, .macropore_water_flux_m3_per_step = macropore, .vapor_flux_m3_per_step = vapor, .heat_flux_megajoules_per_step = heat };
}

/// Rebinds the compact active-face cache after an authoritative hydrology
/// checkpoint owner is restored. Topology remains runtime-defined, while the
/// accepted previous-step water, vapor, and heat fluxes come from `State`.
pub fn refreshSoilFacesFromHydrology(
    faces: *SoilFaces,
    hydrology: *const State,
) !void {
    const count = faces.direction_axis.len;
    if (faces.micropore_faces.len != count or
        faces.macropore_faces.len != count or
        faces.micropore_water_flux_m3_per_step.len != count or
        faces.macropore_water_flux_m3_per_step.len != count or
        faces.vapor_flux_m3_per_step.len != count or
        faces.heat_flux_megajoules_per_step.len != count)
        return error.SoilFaceCacheDimensionMismatch;
    const soil_layer_count = try std.math.mul(
        usize,
        try std.math.mul(usize, hydrology.columns, hydrology.rows),
        hydrology.soil_layer_capacity,
    );
    for (0..count) |face_index| {
        const micropore_face = &faces.micropore_faces[face_index];
        const macropore_face = &faces.macropore_faces[face_index];
        const axis: usize = faces.direction_axis[face_index];
        if (axis > 2 or
            micropore_face.first_cell >= soil_layer_count or
            micropore_face.second_cell >= soil_layer_count or
            micropore_face.first_cell != macropore_face.first_cell or
            micropore_face.second_cell != macropore_face.second_cell)
            return error.InvalidSoilFaceTopology;
        const state_index = micropore_face.first_cell * 3 + axis;
        const micropore_flux =
            hydrology.micropore_face_flux_m3_per_step[state_index];
        const macropore_flux =
            hydrology.macropore_face_flux_m3_per_step[state_index];
        const vapor_flux = hydrology.vapor_face_flux_m3_per_step[state_index];
        const heat_flux = hydrology.heat_face_flux_megajoules_per_step[state_index];
        inline for (.{
            micropore_flux,
            macropore_flux,
            vapor_flux,
            heat_flux,
        }) |value| if (!std.math.isFinite(value))
            return error.NonFiniteSoilFaceCache;
        micropore_face.water_flux_m3_per_step = micropore_flux;
        macropore_face.water_flux_m3_per_step = macropore_flux;
        faces.micropore_water_flux_m3_per_step[face_index] = micropore_flux;
        faces.macropore_water_flux_m3_per_step[face_index] = macropore_flux;
        faces.vapor_flux_m3_per_step[face_index] = vapor_flux;
        faces.heat_flux_megajoules_per_step[face_index] = heat_flux;
    }
}

fn appendFace(hydrology: *const State, first: usize, second: usize, direction: u2, direction_axis: []u2, micropore_faces: []solute.Face, macropore_faces: []solute.Face, micropore: []f64, macropore: []f64, vapor: []f64, heat: []f64, face_index: *usize) void {
    const micro_flux = hydrology.micropore_face_flux_m3_per_step[first * 3 + direction];
    const macro_flux = hydrology.macropore_face_flux_m3_per_step[first * 3 + direction];
    micropore_faces[face_index.*] = .{ .first_cell = first, .second_cell = second, .water_flux_m3_per_step = micro_flux };
    macropore_faces[face_index.*] = .{ .first_cell = first, .second_cell = second, .water_flux_m3_per_step = macro_flux };
    direction_axis[face_index.*] = direction;
    micropore[face_index.*] = micro_flux;
    macropore[face_index.*] = macro_flux;
    vapor[face_index.*] = hydrology.vapor_face_flux_m3_per_step[first * 3 + direction];
    heat[face_index.*] = hydrology.heat_face_flux_megajoules_per_step[first * 3 + direction];
    face_index.* += 1;
}

fn zeroes(allocator: std.mem.Allocator, count: usize) ![]f64 {
    const values = try allocator.alloc(f64, count);
    @memset(values, 0);
    return values;
}

test "transport hydrology retains separate pore domains and builds runtime faces" {
    const config = try @import("config.zig").SimulationConfig.init(.{ .lon_count = 2, .lat_count = 2, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 4 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    @memset(grid.matrix_liquid_water_m3, 1);
    @memset(grid.macropore_liquid_water_m3, 0.25);
    @memset(grid.air_volume_m3, 0.5);
    var snow = try snow_module.State.init(std.testing.allocator, 4, 3);
    defer snow.deinit();
    var hydrology = try State.init(std.testing.allocator, 2, 2, 2, 3);
    defer hydrology.deinit();
    try hydrology.syncStorage(&grid, &snow);
    hydrology.micropore_face_flux_m3_per_step[0] = 0.1;
    hydrology.macropore_face_flux_m3_per_step[0] = 0.02;
    var faces = try buildSoilFaces(std.testing.allocator, &hydrology, &grid);
    defer faces.deinit();
    try std.testing.expectEqual(@as(usize, 12), faces.micropore_faces.len);
    try std.testing.expectEqual(@as(usize, 12), faces.macropore_faces.len);
    try std.testing.expectEqual(@as(f64, 0.1), faces.micropore_water_flux_m3_per_step[0]);
    try std.testing.expectEqual(@as(f64, 0.02), faces.macropore_water_flux_m3_per_step[0]);
    hydrology.micropore_face_flux_m3_per_step[0] = -0.35;
    hydrology.macropore_face_flux_m3_per_step[0] = 0.07;
    hydrology.vapor_face_flux_m3_per_step[0] = 0.004;
    hydrology.heat_face_flux_megajoules_per_step[0] = -1.25;
    try refreshSoilFacesFromHydrology(&faces, &hydrology);
    try std.testing.expectEqual(
        @as(f64, -0.35),
        faces.micropore_faces[0].water_flux_m3_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.07),
        faces.macropore_water_flux_m3_per_step[0],
    );
    try std.testing.expectEqual(@as(f64, 0.004), faces.vapor_flux_m3_per_step[0]);
    try std.testing.expectEqual(@as(f64, -1.25), faces.heat_flux_megajoules_per_step[0]);
    try std.testing.expectEqual(@as(f64, 1), hydrology.micropore_water_volume_m3[0]);
    var standalone_faces = try buildSoilFacesMapped(
        std.testing.allocator,
        &hydrology,
        &grid,
        &.{ 3, 3, 3, 3 },
    );
    defer standalone_faces.deinit();
    try std.testing.expectEqual(@as(usize, 4), standalone_faces.micropore_faces.len);
}
