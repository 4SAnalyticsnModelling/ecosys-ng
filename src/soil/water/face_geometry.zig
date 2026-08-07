const std = @import("std");
const grid_module = @import("../../state/grid.zig");
const hydrology_module = @import("../../transport/hydrology.zig");

pub const State = struct {
    allocator: std.mem.Allocator,
    source_path_length_m: []f64,
    destination_path_length_m: []f64,
    face_area_m2: []f64,
    macropore_hydraulic_conductance_m_per_h_megapascal: []f64,

    /// Reconstructs STARTS `DLYR`, `AREA`, and `DIST` geometry for the shared
    /// +x/+y/+z face topology. AREA is taken from the first/source cell, as in
    /// `XDPTH(N,...)=AREA(N,N3,N2,N1)/DIST(...)`.
    pub fn initMapped(
        allocator: std.mem.Allocator,
        grid: *const grid_module.GridState,
        faces: *const hydrology_module.SoilFaces,
        layer_thickness_m: []const f64,
        horizontal_cell_width_m: []const f64,
        vertical_cell_width_m: []const f64,
    ) !State {
        const count = faces.micropore_faces.len;
        if (faces.macropore_faces.len != count or faces.direction_axis.len != count or
            layer_thickness_m.len != grid.layer_count or
            horizontal_cell_width_m.len != grid.cell_count or
            vertical_cell_width_m.len != grid.cell_count)
            return error.SoilFaceGeometryDimensionMismatch;
        var result: State = undefined;
        result.allocator = allocator;
        result.source_path_length_m = try allocator.alloc(f64, count);
        errdefer allocator.free(result.source_path_length_m);
        result.destination_path_length_m = try allocator.alloc(f64, count);
        errdefer allocator.free(result.destination_path_length_m);
        result.face_area_m2 = try allocator.alloc(f64, count);
        errdefer allocator.free(result.face_area_m2);
        result.macropore_hydraulic_conductance_m_per_h_megapascal = try allocator.alloc(f64, count);
        errdefer allocator.free(result.macropore_hydraulic_conductance_m_per_h_megapascal);
        @memset(result.macropore_hydraulic_conductance_m_per_h_megapascal, 0);
        try result.refreshMapped(grid, faces, layer_thickness_m, horizontal_cell_width_m, vertical_cell_width_m);
        try result.validateFinite();
        return result;
    }

    /// Recomputes all DLYR-dependent face paths and areas without allocation.
    /// Validation is completed for the entire topology before cached values
    /// are changed, so it can participate in a wider atomic geometry commit.
    pub fn refreshMapped(
        self: *State,
        grid: *const grid_module.GridState,
        faces: *const hydrology_module.SoilFaces,
        layer_thickness_m: []const f64,
        horizontal_cell_width_m: []const f64,
        vertical_cell_width_m: []const f64,
    ) !void {
        try self.validateMapped(grid, faces, layer_thickness_m, horizontal_cell_width_m, vertical_cell_width_m);
        for (faces.micropore_faces, faces.direction_axis, 0..) |face, axis, face_index| {
            const source = dimensions(grid, face.first_cell, layer_thickness_m, horizontal_cell_width_m, vertical_cell_width_m) catch unreachable;
            const destination = dimensions(grid, face.second_cell, layer_thickness_m, horizontal_cell_width_m, vertical_cell_width_m) catch unreachable;
            switch (axis) {
                0 => {
                    self.source_path_length_m[face_index] = source.x_width_m;
                    self.destination_path_length_m[face_index] = destination.x_width_m;
                    self.face_area_m2[face_index] = source.layer_thickness_m * source.y_width_m;
                },
                1 => {
                    self.source_path_length_m[face_index] = source.y_width_m;
                    self.destination_path_length_m[face_index] = destination.y_width_m;
                    self.face_area_m2[face_index] = source.layer_thickness_m * source.x_width_m;
                },
                2 => {
                    self.source_path_length_m[face_index] = source.layer_thickness_m;
                    self.destination_path_length_m[face_index] = destination.layer_thickness_m;
                    self.face_area_m2[face_index] = source.x_width_m * source.y_width_m;
                },
                else => unreachable,
            }
        }
    }

    pub fn validateMapped(
        self: *const State,
        grid: *const grid_module.GridState,
        faces: *const hydrology_module.SoilFaces,
        layer_thickness_m: []const f64,
        horizontal_cell_width_m: []const f64,
        vertical_cell_width_m: []const f64,
    ) !void {
        const count = faces.micropore_faces.len;
        if (faces.macropore_faces.len != count or faces.direction_axis.len != count or
            self.source_path_length_m.len != count or self.destination_path_length_m.len != count or
            self.face_area_m2.len != count or layer_thickness_m.len != grid.layer_count or
            horizontal_cell_width_m.len != grid.cell_count or
            vertical_cell_width_m.len != grid.cell_count)
            return error.SoilFaceGeometryDimensionMismatch;
        for (faces.micropore_faces, faces.direction_axis) |face, axis| {
            if (face.first_cell >= grid.layer_count or face.second_cell >= grid.layer_count or axis > 2) return error.InvalidSoilFaceTopology;
            _ = try dimensions(grid, face.first_cell, layer_thickness_m, horizontal_cell_width_m, vertical_cell_width_m);
            _ = try dimensions(grid, face.second_cell, layer_thickness_m, horizontal_cell_width_m, vertical_cell_width_m);
        }
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.macropore_hydraulic_conductance_m_per_h_megapascal);
        self.allocator.free(self.face_area_m2);
        self.allocator.free(self.destination_path_length_m);
        self.allocator.free(self.source_path_length_m);
        self.* = undefined;
    }

    pub fn validateFinite(self: *const State) !void {
        inline for (.{ self.source_path_length_m, self.destination_path_length_m, self.face_area_m2 }) |values| for (values) |value| if (!std.math.isFinite(value) or value <= 0) return error.InvalidSoilFaceGeometry;
        inline for (.{self.macropore_hydraulic_conductance_m_per_h_megapascal}) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSoilFaceGeometry;
    }
};

const Dimensions = struct { x_width_m: f64, y_width_m: f64, layer_thickness_m: f64 };

fn dimensions(grid: *const grid_module.GridState, layer_cell: usize, layer_thickness_m: []const f64, horizontal_cell_width_m: []const f64, vertical_cell_width_m: []const f64) !Dimensions {
    const horizontal_cell = layer_cell / grid.soil_layer_capacity;
    if (horizontal_cell >= grid.cell_count) return error.InvalidSoilFaceTopology;
    if (horizontal_cell_width_m.len != grid.cell_count or vertical_cell_width_m.len != grid.cell_count) return error.CellGeometryDimensionMismatch;
    const result: Dimensions = .{ .x_width_m = horizontal_cell_width_m[horizontal_cell], .y_width_m = vertical_cell_width_m[horizontal_cell], .layer_thickness_m = layer_thickness_m[layer_cell] };
    inline for (.{ result.x_width_m, result.y_width_m, result.layer_thickness_m }) |value| if (!std.math.isFinite(value) or value <= 0) return error.InvalidSoilFaceGeometry;
    return result;
}

test "mapped faces reproduce STARTS source AREA and full DLYR paths" {
    const cfg = try @import("../../core/config.zig").SimulationConfig.init(.{ .lon_count = 2, .lat_count = 1, .soil_layers = 2, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 2 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var grid = try grid_module.GridState.init(std.testing.allocator, cfg);
    defer grid.deinit();
    @memset(grid.active_soil_layer_count, 2);
    var snow = try @import("../solute/snow_solute_transport.zig").State.init(std.testing.allocator, 2, 1);
    defer snow.deinit();
    var hydrology = try hydrology_module.State.init(std.testing.allocator, 2, 1, 2, 1);
    defer hydrology.deinit();
    var faces = try hydrology_module.buildSoilFaces(std.testing.allocator, &hydrology, &grid);
    defer faces.deinit();
    var geometry = try State.initMapped(std.testing.allocator, &grid, &faces, &.{ 0.1, 0.2, 0.1, 0.2 }, &.{ 3, 5 }, &.{ 4, 4 });
    defer geometry.deinit();
    // First face is +x from cell 0 layer 0: AREA(1)=DLYR(3)*DLYR(2).
    try std.testing.expectEqual(@as(f64, 3), geometry.source_path_length_m[0]);
    try std.testing.expectEqual(@as(f64, 5), geometry.destination_path_length_m[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), geometry.face_area_m2[0], 1e-15);
    // Second face is +z inside cell 0: AREA(3)=DLYR(1)*DLYR(2).
    try std.testing.expectEqual(@as(f64, 0.1), geometry.source_path_length_m[1]);
    try std.testing.expectEqual(@as(f64, 0.2), geometry.destination_path_length_m[1]);
    try std.testing.expectEqual(@as(f64, 12), geometry.face_area_m2[1]);
    const updated_thickness = [_]f64{ 0.3, 0.2, 0.1, 0.2 };
    try geometry.refreshMapped(&grid, &faces, &updated_thickness, &.{ 3, 5 }, &.{ 4, 4 });
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), geometry.face_area_m2[0], 1e-15);
    try std.testing.expectEqual(@as(f64, 0.3), geometry.source_path_length_m[1]);
}
