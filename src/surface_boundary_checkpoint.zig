const std = @import("std");
const GroundAir = @import("ground_air_exchange.zig").State;
const SurfaceAerodynamics = @import("surface_aerodynamics.zig").State;

pub const ground_air_field_count: usize = 4;
pub const aerodynamic_field_count: usize = 5;

pub const View = struct {
    ground_air: *const GroundAir,
    surface_aerodynamics: *const SurfaceAerodynamics,
};

pub const Snapshot = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    ground_air_fields: [ground_air_field_count][]f64,
    ground_air_iteration_count: []u16,
    aerodynamic_fields: [aerodynamic_field_count][]f64,

    pub fn deinit(self: *Snapshot) void {
        for (self.aerodynamic_fields) |values| self.allocator.free(values);
        self.allocator.free(self.ground_air_iteration_count);
        for (self.ground_air_fields) |values| self.allocator.free(values);
        self.* = undefined;
    }

    pub fn validate(self: Snapshot) !void {
        if (self.cell_count == 0 or
            self.ground_air_iteration_count.len != self.cell_count)
            return error.SurfaceBoundaryCheckpointDimensionMismatch;
        for (self.ground_air_fields) |values| {
            if (values.len != self.cell_count)
                return error.SurfaceBoundaryCheckpointDimensionMismatch;
            try validateFinite(values);
        }
        for (self.aerodynamic_fields) |values| {
            if (values.len != self.cell_count)
                return error.SurfaceBoundaryCheckpointDimensionMismatch;
            try validateFinite(values);
        }
        for (0..self.cell_count) |cell| {
            if (self.ground_air_fields[0][cell] <= 0 or
                self.ground_air_fields[1][cell] < 0 or
                self.ground_air_fields[2][cell] <= 0 or
                self.ground_air_fields[3][cell] <= 0 or
                self.aerodynamic_fields[0][cell] < 0 or
                self.aerodynamic_fields[1][cell] <= 0 or
                self.aerodynamic_fields[2][cell] <= 0 or
                self.aerodynamic_fields[3][cell] < 0 or
                self.aerodynamic_fields[4][cell] < 0)
                return error.InvalidSurfaceBoundaryCheckpointState;
        }
    }

    pub fn restoreInto(
        self: Snapshot,
        ground_air: *GroundAir,
        surface_aerodynamics: *SurfaceAerodynamics,
    ) !void {
        try self.validate();
        try validateTargetDimensions(.{
            .ground_air = ground_air,
            .surface_aerodynamics = surface_aerodynamics,
        });
        if (ground_air.cell_count != self.cell_count)
            return error.SurfaceBoundaryCheckpointDimensionMismatch;
        inline for (groundAirFields(ground_air), self.ground_air_fields) |target, source|
            @memcpy(target, source);
        @memcpy(ground_air.iteration_count, self.ground_air_iteration_count);
        inline for (
            aerodynamicFields(surface_aerodynamics),
            self.aerodynamic_fields,
        ) |target, source| @memcpy(target, source);
    }
};

pub fn write(writer: anytype, view: View) !void {
    try validateView(view);
    try writer.writeInt(u64, @intCast(view.ground_air.cell_count), .little);
    inline for (groundAirConstFields(view.ground_air)) |values|
        try writeF64Slice(writer, values);
    for (view.ground_air.iteration_count) |value|
        try writer.writeInt(u16, value, .little);
    inline for (aerodynamicConstFields(view.surface_aerodynamics)) |values|
        try writeF64Slice(writer, values);
}

pub fn read(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    expected_cell_count: usize,
) !Snapshot {
    if (try reader.takeInt(u64, .little) != expected_cell_count)
        return error.SurfaceBoundaryCheckpointDimensionMismatch;
    var result: Snapshot = .{
        .allocator = allocator,
        .cell_count = expected_cell_count,
        .ground_air_fields = undefined,
        .ground_air_iteration_count = undefined,
        .aerodynamic_fields = undefined,
    };
    var ground_allocated: usize = 0;
    var aerodynamic_allocated: usize = 0;
    var iterations_allocated = false;
    errdefer {
        for (result.aerodynamic_fields[0..aerodynamic_allocated]) |values|
            allocator.free(values);
        if (iterations_allocated) allocator.free(result.ground_air_iteration_count);
        for (result.ground_air_fields[0..ground_allocated]) |values|
            allocator.free(values);
    }
    for (&result.ground_air_fields) |*values| {
        values.* = try allocator.alloc(f64, expected_cell_count);
        ground_allocated += 1;
        try readF64Slice(reader, values.*);
    }
    result.ground_air_iteration_count =
        try allocator.alloc(u16, expected_cell_count);
    iterations_allocated = true;
    for (result.ground_air_iteration_count) |*value|
        value.* = try reader.takeInt(u16, .little);
    for (&result.aerodynamic_fields) |*values| {
        values.* = try allocator.alloc(f64, expected_cell_count);
        aerodynamic_allocated += 1;
        try readF64Slice(reader, values.*);
    }
    try result.validate();
    return result;
}

pub fn validateView(view: View) !void {
    try validateTargetDimensions(view);
    const cells = view.ground_air.cell_count;
    inline for (groundAirConstFields(view.ground_air)) |values|
        try validateFinite(values);
    inline for (aerodynamicConstFields(view.surface_aerodynamics)) |values|
        try validateFinite(values);
    for (0..cells) |cell| {
        if (view.ground_air.temperature_k[cell] <= 0 or
            view.ground_air.vapor_volume_fraction[cell] < 0 or
            view.ground_air.heat_capacity_megajoules_per_k[cell] <= 0 or
            view.ground_air.air_volume_m3[cell] <= 0 or
            view.surface_aerodynamics.zero_plane_displacement_m[cell] < 0 or
            view.surface_aerodynamics.effective_roughness_height_m[cell] <= 0 or
            view.surface_aerodynamics.wind_reference_height_m[cell] <= 0 or
            view.surface_aerodynamics.bulk_richardson_coefficient_k[cell] < 0 or
            view.surface_aerodynamics.isothermal_aerodynamic_resistance_h_per_m[cell] < 0)
            return error.InvalidSurfaceBoundaryCheckpointState;
    }
}

pub fn validateTargetDimensions(view: View) !void {
    const cells = view.ground_air.cell_count;
    if (cells == 0 or view.surface_aerodynamics.cell_count != cells or
        view.ground_air.iteration_count.len != cells)
        return error.SurfaceBoundaryCheckpointDimensionMismatch;
    inline for (groundAirConstFields(view.ground_air)) |values|
        if (values.len != cells)
            return error.SurfaceBoundaryCheckpointDimensionMismatch;
    inline for (aerodynamicConstFields(view.surface_aerodynamics)) |values|
        if (values.len != cells)
            return error.SurfaceBoundaryCheckpointDimensionMismatch;
}

fn groundAirConstFields(state: *const GroundAir) [ground_air_field_count][]const f64 {
    return .{
        state.temperature_k,
        state.vapor_volume_fraction,
        state.heat_capacity_megajoules_per_k,
        state.air_volume_m3,
    };
}

fn groundAirFields(state: *GroundAir) [ground_air_field_count][]f64 {
    return .{
        state.temperature_k,
        state.vapor_volume_fraction,
        state.heat_capacity_megajoules_per_k,
        state.air_volume_m3,
    };
}

fn aerodynamicConstFields(
    state: *const SurfaceAerodynamics,
) [aerodynamic_field_count][]const f64 {
    return .{
        state.zero_plane_displacement_m,
        state.effective_roughness_height_m,
        state.wind_reference_height_m,
        state.bulk_richardson_coefficient_k,
        state.isothermal_aerodynamic_resistance_h_per_m,
    };
}

fn aerodynamicFields(
    state: *SurfaceAerodynamics,
) [aerodynamic_field_count][]f64 {
    return .{
        state.zero_plane_displacement_m,
        state.effective_roughness_height_m,
        state.wind_reference_height_m,
        state.bulk_richardson_coefficient_k,
        state.isothermal_aerodynamic_resistance_h_per_m,
    };
}

fn writeF64Slice(writer: anytype, values: []const f64) !void {
    for (values) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceBoundaryCheckpoint;
        try writer.writeInt(u64, @bitCast(value), .little);
    }
}

fn readF64Slice(reader: *std.Io.Reader, values: []f64) !void {
    for (values) |*value| {
        value.* = @bitCast(try reader.takeInt(u64, .little));
        if (!std.math.isFinite(value.*))
            return error.NonFiniteSurfaceBoundaryCheckpoint;
    }
}

fn validateFinite(values: []const f64) !void {
    for (values) |value| if (!std.math.isFinite(value))
        return error.NonFiniteSurfaceBoundaryCheckpoint;
}
