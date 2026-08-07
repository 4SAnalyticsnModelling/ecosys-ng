const std = @import("std");
const Geometry = @import("../../soil/profile/layer_geometry.zig").State;
const Hydrology = @import("../../transport/hydrology.zig").State;
const Surface = @import("../../surface/precipitation.zig").RuntimeState;
const Erosion = @import("../../soil/profile/erosion.zig").RuntimeState;
const Climate = @import("../input/climate_change.zig").State;
const ErodedMinerals = @import("../../soil/profile/erosion_mineral_bridge.zig").State;
const Runtime = @import("soil_runtime_checkpoint.zig");
const SurfaceBoundary = @import("surface_boundary_checkpoint.zig");
const SurfaceLitterGeometry = @import("../../surface/litter_geometry_step.zig").State;

const magic = "ECOSGEOM";
const version: u32 = 12;

pub const View = struct {
    geometry: *const Geometry,
    hydrology: *const Hydrology,
    surface: *const Surface,
    erosion: *const Erosion,
    climate: *const Climate,
    eroded_minerals: *const ErodedMinerals,
    runtime: ?Runtime.View,
    surface_boundary: ?SurfaceBoundary.View,
    surface_litter_geometry: *const SurfaceLitterGeometry,
    surface_litter_ice_m3: []const f64,
    delayed_live_canopy_combustion_heat_megajoules: []const f64,
    delayed_standing_dead_combustion_heat_megajoules: []const f64,
    delayed_subsurface_combustion_heat_megajoules: []const f64,
    delayed_surface_combustion_heat_megajoules: []const f64,
};

pub const Limits = struct {
    maximum_columns: usize,
    maximum_rows: usize,
    maximum_soil_layers: usize,
    maximum_snow_layers: usize,
    maximum_plants: usize,
};

pub const Owned = struct {
    geometry: Geometry,
    hydrology: Hydrology,
    surface: Surface,
    erosion: Erosion,
    climate: Climate,
    eroded_minerals: ErodedMinerals,
    runtime: Runtime.Snapshot,
    surface_boundary: SurfaceBoundary.Snapshot,
    surface_litter_geometry: SurfaceLitterGeometry,
    surface_litter_ice_m3: []f64,
    allocator: std.mem.Allocator,
    delayed_live_canopy_combustion_heat_megajoules: []f64,
    delayed_standing_dead_combustion_heat_megajoules: []f64,
    delayed_subsurface_combustion_heat_megajoules: []f64,
    delayed_surface_combustion_heat_megajoules: []f64,

    pub fn deinit(self: *Owned) void {
        self.surface_litter_geometry.deinit();
        self.surface_boundary.deinit();
        self.allocator.free(self.surface_litter_ice_m3);
        self.runtime.deinit();
        self.allocator.free(self.delayed_surface_combustion_heat_megajoules);
        self.allocator.free(self.delayed_subsurface_combustion_heat_megajoules);
        self.allocator.free(self.delayed_standing_dead_combustion_heat_megajoules);
        self.allocator.free(self.delayed_live_canopy_combustion_heat_megajoules);
        self.eroded_minerals.deinit();
        self.erosion.deinit();
        self.surface.deinit();
        self.hydrology.deinit();
        self.geometry.deinit();
        self.* = undefined;
    }
};

pub fn write(writer: anytype, view: View) !void {
    try validate(view);
    const runtime = view.runtime orelse
        return error.MissingSoilRuntimeCheckpointState;
    const surface_boundary = view.surface_boundary orelse
        return error.MissingSurfaceBoundaryCheckpointState;
    try writer.writeAll(magic);
    try writer.writeInt(u32, version, .little);
    try writer.writeInt(u64, @intCast(view.hydrology.columns), .little);
    try writer.writeInt(u64, @intCast(view.hydrology.rows), .little);
    try writer.writeInt(u64, @intCast(view.hydrology.soil_layer_capacity), .little);
    try writer.writeInt(u64, @intCast(view.hydrology.snow_layer_capacity), .little);
    try writer.writeInt(u64, @intCast(view.delayed_live_canopy_combustion_heat_megajoules.len), .little);
    try writeUsizeSlice(writer, view.geometry.first_active_layer);
    try writeUsizeSlice(writer, view.geometry.active_layer_count);
    inline for (@typeInfo(Geometry).@"struct".fields) |field| {
        if (field.type == []f64) try writeF64Slice(writer, @field(view.geometry, field.name));
    }
    inline for (@typeInfo(Hydrology).@"struct".fields) |field| {
        if (field.type == []f64) try writeF64Slice(writer, @field(view.hydrology, field.name));
    }
    inline for (@typeInfo(Surface).@"struct".fields) |field| switch (field.type) {
        []f64 => try writeF64Slice(writer, @field(view.surface, field.name)),
        []bool => try writeBoolSlice(writer, @field(view.surface, field.name)),
        else => {},
    };
    inline for (@typeInfo(SurfaceLitterGeometry).@"struct".fields) |field|
        if (field.type == []f64)
            try writeF64Slice(
                writer,
                @field(view.surface_litter_geometry, field.name),
            );
    try writeF64Slice(writer, view.erosion.surface_sediment_megagrams);
    try writeF64Slice(writer, view.erosion.surface_soil_mass_megagrams);
    try writeBoolSlice(writer, view.erosion.surface_soil_mass_initialized);
    for (view.climate.modifiers) |modifier| inline for (@typeInfo(@TypeOf(modifier)).@"struct".fields) |field| {
        const value = @field(modifier, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteSoilGeometryCheckpoint;
        try writer.writeInt(u64, @bitCast(value), .little);
    };
    try writer.writeByte(@intFromBool(view.eroded_minerals.initialized));
    try writeF64Slice(writer, view.eroded_minerals.workspace.pools);
    try writeF64Slice(writer, view.eroded_minerals.workspace.exported);
    try writeF64Slice(writer, view.delayed_live_canopy_combustion_heat_megajoules);
    try writeF64Slice(writer, view.delayed_standing_dead_combustion_heat_megajoules);
    try writeF64Slice(writer, view.delayed_subsurface_combustion_heat_megajoules);
    try writeF64Slice(writer, view.delayed_surface_combustion_heat_megajoules);
    try writeF64Slice(writer, view.surface_litter_ice_m3);
    try Runtime.write(writer, runtime);
    try SurfaceBoundary.write(writer, surface_boundary);
}

pub fn read(allocator: std.mem.Allocator, reader: *std.Io.Reader, limits: Limits) !Owned {
    if (limits.maximum_columns == 0 or limits.maximum_rows == 0 or limits.maximum_soil_layers == 0 or limits.maximum_snow_layers == 0) return error.InvalidSoilGeometryCheckpointLimits;
    if (!std.mem.eql(u8, try reader.takeArray(magic.len), magic)) return error.InvalidSoilGeometryCheckpointMagic;
    const file_version = try reader.takeInt(u32, .little);
    if (file_version != version)
        return error.UnsupportedSoilGeometryCheckpointVersion;
    const columns = try bounded(reader, limits.maximum_columns, error.SoilGeometryCheckpointColumnLimitExceeded);
    const rows = try bounded(reader, limits.maximum_rows, error.SoilGeometryCheckpointRowLimitExceeded);
    const soil_layers = try bounded(reader, limits.maximum_soil_layers, error.SoilGeometryCheckpointSoilLayerLimitExceeded);
    const snow_layers = try bounded(reader, limits.maximum_snow_layers, error.SoilGeometryCheckpointSnowLayerLimitExceeded);
    const plants = try bounded(reader, limits.maximum_plants, error.SoilGeometryCheckpointPlantLimitExceeded);
    if (columns == 0 or rows == 0 or soil_layers == 0 or snow_layers == 0) return error.InvalidSoilGeometryCheckpointDimensions;
    const cells = try std.math.mul(usize, columns, rows);
    var geometry = try Geometry.init(allocator, cells, soil_layers);
    errdefer geometry.deinit();
    var hydrology = try Hydrology.init(allocator, columns, rows, soil_layers, snow_layers);
    errdefer hydrology.deinit();
    var surface = try Surface.init(allocator, cells);
    errdefer surface.deinit();
    var surface_litter_geometry =
        try SurfaceLitterGeometry.init(allocator, cells);
    errdefer surface_litter_geometry.deinit();
    var erosion = try Erosion.init(allocator, columns, rows);
    errdefer erosion.deinit();
    var climate: Climate = .{};
    var eroded_minerals = try ErodedMinerals.init(allocator, cells);
    errdefer eroded_minerals.deinit();
    const layer_cells = try std.math.mul(usize, cells, soil_layers);
    const delayed_live_heat = try allocator.alloc(f64, plants);
    errdefer allocator.free(delayed_live_heat);
    const delayed_dead_heat = try allocator.alloc(f64, plants);
    errdefer allocator.free(delayed_dead_heat);
    const delayed_subsurface_heat = try allocator.alloc(f64, layer_cells);
    errdefer allocator.free(delayed_subsurface_heat);
    const delayed_surface_heat = try allocator.alloc(f64, cells);
    errdefer allocator.free(delayed_surface_heat);
    const surface_litter_ice_m3 = try allocator.alloc(f64, cells);
    errdefer allocator.free(surface_litter_ice_m3);
    try readUsizeSlice(reader, geometry.first_active_layer);
    try readUsizeSlice(reader, geometry.active_layer_count);
    inline for (@typeInfo(Geometry).@"struct".fields) |field| {
        if (field.type == []f64) try readF64Slice(reader, @field(geometry, field.name));
    }
    inline for (@typeInfo(Hydrology).@"struct".fields) |field| {
        if (field.type == []f64)
            try readF64Slice(reader, @field(hydrology, field.name));
    }
    inline for (@typeInfo(Surface).@"struct".fields) |field| switch (field.type) {
        []f64 => try readF64Slice(reader, @field(surface, field.name)),
        []bool => try readBoolSlice(reader, @field(surface, field.name)),
        else => {},
    };
    inline for (@typeInfo(SurfaceLitterGeometry).@"struct".fields) |field|
        if (field.type == []f64)
            try readF64Slice(
                reader,
                @field(surface_litter_geometry, field.name),
            );
    try readF64Slice(reader, erosion.surface_sediment_megagrams);
    try readF64Slice(reader, erosion.surface_soil_mass_megagrams);
    try readBoolSlice(reader, erosion.surface_soil_mass_initialized);
    for (&climate.modifiers) |*modifier| inline for (@typeInfo(@TypeOf(modifier.*)).@"struct".fields) |field| {
        @field(modifier.*, field.name) = @bitCast(try reader.takeInt(u64, .little));
        if (!std.math.isFinite(@field(modifier.*, field.name))) return error.NonFiniteSoilGeometryCheckpoint;
    };
    eroded_minerals.initialized = switch (try reader.takeByte()) {
        0 => false,
        1 => true,
        else => return error.InvalidSoilGeometryCheckpointBoolean,
    };
    try readF64Slice(reader, eroded_minerals.workspace.pools);
    try readF64Slice(reader, eroded_minerals.workspace.exported);
    try readF64Slice(reader, delayed_live_heat);
    try readF64Slice(reader, delayed_dead_heat);
    try readF64Slice(reader, delayed_subsurface_heat);
    try readF64Slice(reader, delayed_surface_heat);
    try readF64Slice(reader, surface_litter_ice_m3);
    var runtime = try Runtime.read(allocator, reader, layer_cells);
    errdefer runtime.deinit();
    var surface_boundary = try SurfaceBoundary.read(allocator, reader, cells);
    errdefer surface_boundary.deinit();
    if (reader.peekByte()) |_| {
        return error.TrailingSoilGeometryCheckpointData;
    } else |err| switch (err) {
        error.EndOfStream => {},
        else => return err,
    }
    const result = Owned{ .geometry = geometry, .hydrology = hydrology, .surface = surface, .erosion = erosion, .climate = climate, .eroded_minerals = eroded_minerals, .runtime = runtime, .surface_boundary = surface_boundary, .surface_litter_geometry = surface_litter_geometry, .surface_litter_ice_m3 = surface_litter_ice_m3, .allocator = allocator, .delayed_live_canopy_combustion_heat_megajoules = delayed_live_heat, .delayed_standing_dead_combustion_heat_megajoules = delayed_dead_heat, .delayed_subsurface_combustion_heat_megajoules = delayed_subsurface_heat, .delayed_surface_combustion_heat_megajoules = delayed_surface_heat };
    try validate(.{ .geometry = &result.geometry, .hydrology = &result.hydrology, .surface = &result.surface, .erosion = &result.erosion, .climate = &result.climate, .eroded_minerals = &result.eroded_minerals, .runtime = null, .surface_boundary = null, .surface_litter_geometry = &result.surface_litter_geometry, .surface_litter_ice_m3 = result.surface_litter_ice_m3, .delayed_live_canopy_combustion_heat_megajoules = result.delayed_live_canopy_combustion_heat_megajoules, .delayed_standing_dead_combustion_heat_megajoules = result.delayed_standing_dead_combustion_heat_megajoules, .delayed_subsurface_combustion_heat_megajoules = result.delayed_subsurface_combustion_heat_megajoules, .delayed_surface_combustion_heat_megajoules = result.delayed_surface_combustion_heat_megajoules });
    return result;
}

fn validate(view: View) !void {
    const geometry = view.geometry;
    const hydrology = view.hydrology;
    const cells = std.math.mul(usize, hydrology.columns, hydrology.rows) catch return error.InvalidSoilGeometryCheckpointDimensions;
    if (hydrology.columns == 0 or hydrology.rows == 0 or hydrology.soil_layer_capacity == 0 or hydrology.snow_layer_capacity == 0 or geometry.cell_count != cells or geometry.layer_capacity != hydrology.soil_layer_capacity or view.surface.cell_count != cells or view.surface_litter_geometry.cell_count != cells or view.erosion.cell_count != cells or view.delayed_live_canopy_combustion_heat_megajoules.len == 0 or view.delayed_standing_dead_combustion_heat_megajoules.len != view.delayed_live_canopy_combustion_heat_megajoules.len or view.delayed_subsurface_combustion_heat_megajoules.len != try std.math.mul(usize, cells, hydrology.soil_layer_capacity) or view.delayed_surface_combustion_heat_megajoules.len != cells or view.surface_litter_ice_m3.len != cells) return error.InvalidSoilGeometryCheckpointDimensions;
    for (0..cells) |cell| {
        const first = geometry.first_active_layer[cell];
        const count = geometry.active_layer_count[cell];
        if (count == 0 or first >= geometry.layer_capacity or count > geometry.layer_capacity - first) return error.InvalidCheckpointActiveSoilLayerRange;
        const boundary_base = cell * (geometry.layer_capacity + 1);
        const layer_base = cell * geometry.layer_capacity;
        for (first..first + count) |layer| {
            const top = geometry.boundary_depth_m[boundary_base + layer];
            const bottom = geometry.boundary_depth_m[boundary_base + layer + 1];
            const top_without_freeze = geometry.boundary_depth_without_freeze_m[boundary_base + layer];
            const bottom_without_freeze = geometry.boundary_depth_without_freeze_m[boundary_base + layer + 1];
            if (bottom <= top or bottom_without_freeze <= top_without_freeze or geometry.layer_thickness_m[layer_base + layer] <= 0) return error.InvalidCheckpointSoilLayerGeometry;
        }
    }
    for (view.climate.modifiers) |modifier| inline for (@typeInfo(@TypeOf(modifier)).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(modifier, field.name))) return error.NonFiniteSoilGeometryCheckpoint;
    };
    if (view.eroded_minerals.workspace.cell_count != cells or view.eroded_minerals.workspace.component_count != @import("../../soil/profile/erosion_mineral_bridge.zig").component_count) return error.InvalidSoilGeometryCheckpointDimensions;
    try validateFinite(view.eroded_minerals.workspace.pools);
    try validateFinite(view.eroded_minerals.workspace.exported);
    try validateFinite(view.delayed_live_canopy_combustion_heat_megajoules);
    try validateFinite(view.delayed_standing_dead_combustion_heat_megajoules);
    try validateFinite(view.delayed_subsurface_combustion_heat_megajoules);
    try validateFinite(view.delayed_surface_combustion_heat_megajoules);
    try validateFinite(view.surface_litter_ice_m3);
    inline for (@typeInfo(SurfaceLitterGeometry).@"struct".fields) |field|
        if (field.type == []f64)
            try validateFinite(@field(view.surface_litter_geometry, field.name));
    for (0..cells) |cell| {
        const litter = view.surface_litter_geometry;
        inline for (.{
            litter.water_retention_capacity_m3[cell],
            litter.dry_litter_volume_m3[cell],
            litter.expanded_total_volume_m3[cell],
            litter.dry_mass_megagrams[cell],
            litter.pore_volume_m3[cell],
            litter.air_volume_m3[cell],
            litter.porosity_m3_per_m3[cell],
            litter.field_capacity_m3_per_m3[cell],
            litter.wilting_point_m3_per_m3[cell],
        }) |value| if (value < 0)
            return error.InvalidCheckpointSurfaceLitterGeometry;
        if (litter.air_volume_m3[cell] > litter.pore_volume_m3[cell] or
            litter.pore_volume_m3[cell] >
                litter.expanded_total_volume_m3[cell] or
            litter.porosity_m3_per_m3[cell] > 1 or
            litter.field_capacity_m3_per_m3[cell] >
                litter.porosity_m3_per_m3[cell] or
            litter.wilting_point_m3_per_m3[cell] >
                litter.field_capacity_m3_per_m3[cell])
            return error.InvalidCheckpointSurfaceLitterGeometry;
    }
    inline for (.{ view.delayed_live_canopy_combustion_heat_megajoules, view.delayed_standing_dead_combustion_heat_megajoules, view.delayed_subsurface_combustion_heat_megajoules, view.delayed_surface_combustion_heat_megajoules }) |values| for (values) |value| if (value < 0) return error.InvalidCheckpointDelayedCombustionHeat;
    for (view.surface_litter_ice_m3) |value| if (value < 0)
        return error.InvalidCheckpointSurfaceLitterIce;
    for (view.eroded_minerals.workspace.pools) |value| if (value < -1e-14) return error.InvalidSoilGeometryCheckpointInventory;
    for (view.eroded_minerals.workspace.exported) |value| if (value < -1e-14) return error.InvalidSoilGeometryCheckpointInventory;
    inline for (@typeInfo(Geometry).@"struct".fields) |field| {
        if (field.type == []f64) try validateFinite(@field(geometry, field.name));
    }
    try hydrology.validateFinite();
    inline for (@typeInfo(Surface).@"struct".fields) |field| if (field.type == []f64) try validateFinite(@field(view.surface, field.name));
    if (view.runtime) |runtime| try Runtime.validateView(runtime);
    if (view.surface_boundary) |surface_boundary|
        try SurfaceBoundary.validateView(surface_boundary);
    try validateFinite(view.erosion.surface_sediment_megagrams);
    try validateFinite(view.erosion.surface_soil_mass_megagrams);
    for (0..cells) |cell| {
        if (view.erosion.surface_sediment_megagrams[cell] < -1e-14 or view.erosion.surface_soil_mass_megagrams[cell] < -1e-14) return error.InvalidSoilGeometryCheckpointInventory;
        if (view.erosion.surface_soil_mass_initialized[cell] and view.erosion.surface_soil_mass_megagrams[cell] <= 0) return error.InvalidCheckpointSurfaceSoilMass;
    }
    inline for (.{ hydrology.micropore_water_volume_m3, hydrology.macropore_water_volume_m3, hydrology.matrix_air_volume_m3, hydrology.macropore_air_volume_m3, hydrology.air_volume_m3, hydrology.water_vapor_volume_m3, hydrology.snow_surface_carrier_volume_m3, hydrology.snow_liquid_water_volume_m3 }) |values| {
        for (values) |value| if (value < -1e-14) return error.InvalidSoilGeometryCheckpointInventory;
    }
}

fn validateFinite(values: []const f64) !void {
    for (values) |value| if (!std.math.isFinite(value)) return error.NonFiniteSoilGeometryCheckpoint;
}

fn bounded(reader: *std.Io.Reader, limit: usize, comptime too_large: anyerror) !usize {
    const value = try reader.takeInt(u64, .little);
    if (value > limit or value > std.math.maxInt(usize)) return too_large;
    return @intCast(value);
}

fn writeUsizeSlice(writer: anytype, values: []const usize) !void {
    for (values) |value| try writer.writeInt(u64, @intCast(value), .little);
}

fn readUsizeSlice(reader: *std.Io.Reader, values: []usize) !void {
    for (values) |*value| {
        const stored = try reader.takeInt(u64, .little);
        if (stored > std.math.maxInt(usize)) return error.InvalidCheckpointInteger;
        value.* = @intCast(stored);
    }
}

fn writeF64Slice(writer: anytype, values: []const f64) !void {
    for (values) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteSoilGeometryCheckpoint;
        try writer.writeInt(u64, @bitCast(value), .little);
    }
}

fn readF64Slice(reader: *std.Io.Reader, values: []f64) !void {
    for (values) |*value| {
        value.* = @bitCast(try reader.takeInt(u64, .little));
        if (!std.math.isFinite(value.*)) return error.NonFiniteSoilGeometryCheckpoint;
    }
}
fn writeBoolSlice(writer: anytype, values: []const bool) !void {
    for (values) |value| try writer.writeByte(@intFromBool(value));
}
fn readBoolSlice(reader: *std.Io.Reader, values: []bool) !void {
    for (values) |*value| value.* = switch (try reader.takeByte()) {
        0 => false,
        1 => true,
        else => return error.InvalidSoilGeometryCheckpointBoolean,
    };
}

const TestRuntime = struct {
    allocator: std.mem.Allocator,
    properties: @import("../../soil/water/solver_properties.zig").State,
    thermal: @import("../../soil/heat/thermal.zig").State,
    solver_fields: [Runtime.solver_field_count][]f64,
    thermal_fields: [Runtime.thermal_field_count][]f64,

    fn init(
        allocator: std.mem.Allocator,
        cell_count: usize,
        layer_capacity: usize,
    ) !TestRuntime {
        const layer_count = try std.math.mul(usize, cell_count, layer_capacity);
        var result: TestRuntime = undefined;
        result.allocator = allocator;
        var solver_allocated: usize = 0;
        var thermal_allocated: usize = 0;
        errdefer {
            for (result.thermal_fields[0..thermal_allocated]) |values|
                allocator.free(values);
            for (result.solver_fields[0..solver_allocated]) |values|
                allocator.free(values);
        }
        for (&result.solver_fields) |*values| {
            values.* = try allocator.alloc(f64, layer_count);
            @memset(values.*, 0);
            solver_allocated += 1;
        }
        for (&result.thermal_fields) |*values| {
            values.* = try allocator.alloc(f64, layer_count);
            @memset(values.*, 0);
            thermal_allocated += 1;
        }
        result.properties = undefined;
        result.properties.layer_count = layer_count;
        inline for (@typeInfo(Runtime.SolverField).@"enum".fields, 0..) |field, index|
            @field(result.properties, field.name) = result.solver_fields[index];
        result.thermal = undefined;
        result.thermal.cell_count = cell_count;
        result.thermal.soil_layer_capacity = layer_capacity;
        inline for (@typeInfo(Runtime.ThermalField).@"enum".fields, 0..) |field, index|
            @field(result.thermal, field.name) = result.thermal_fields[index];
        @memset(result.properties.matrix_bulk_volume_m3, 0.5);
        @memset(result.properties.layer_volume_m3, 1);
        @memset(result.properties.layer_thickness_m, 1);
        @memset(result.properties.layer_midpoint_depth_m, 0.5);
        @memset(result.properties.layer_bottom_depth_m, 1);
        @memset(result.properties.bulk_density_megagrams_per_m3, 1);
        @memset(result.properties.porosity_fraction, 0.5);
        @memset(result.properties.rainfall_conductivity_multiplier, 1);
        @memset(result.thermal.layer_volume_m3, 1);
        @memset(result.thermal.layer_thickness_m, 1);
        @memset(result.thermal.porosity_fraction, 0.5);
        @memset(result.thermal.dry_solid_heat_capacity_megajoules_per_m3_k, 1);
        @memset(result.thermal.solid_thermal_conductivity_numerator_m_megajoules_per_h_k, 1);
        @memset(result.thermal.solid_thermal_conductivity_denominator, 1);
        @memset(result.thermal.total_heat_capacity_megajoules_per_m3_k, 2);
        @memset(result.thermal.thermal_conductivity_m_megajoules_per_h_k, 0.5);
        return result;
    }

    fn deinit(self: *TestRuntime) void {
        for (self.thermal_fields) |values| self.allocator.free(values);
        for (self.solver_fields) |values| self.allocator.free(values);
        self.* = undefined;
    }

    fn view(self: *const TestRuntime) Runtime.View {
        return .{
            .soil_properties = &self.properties,
            .soil_thermal = &self.thermal,
        };
    }
};

const TestSurfaceBoundary = struct {
    allocator: std.mem.Allocator,
    ground_air: @import("../../surface/ground_air_exchange.zig").State,
    aerodynamics: @import("../../surface/aerodynamics.zig").State,
    ground_fields: [SurfaceBoundary.ground_air_field_count][]f64,
    iterations: []u16,
    aerodynamic_fields: [SurfaceBoundary.aerodynamic_field_count][]f64,

    fn init(
        allocator: std.mem.Allocator,
        cell_count: usize,
    ) !TestSurfaceBoundary {
        var result: TestSurfaceBoundary = undefined;
        result.allocator = allocator;
        var ground_allocated: usize = 0;
        var aerodynamic_allocated: usize = 0;
        var iterations_allocated = false;
        errdefer {
            for (result.aerodynamic_fields[0..aerodynamic_allocated]) |values|
                allocator.free(values);
            if (iterations_allocated) allocator.free(result.iterations);
            for (result.ground_fields[0..ground_allocated]) |values|
                allocator.free(values);
        }
        for (&result.ground_fields) |*values| {
            values.* = try allocator.alloc(f64, cell_count);
            @memset(values.*, 1);
            ground_allocated += 1;
        }
        result.iterations = try allocator.alloc(u16, cell_count);
        @memset(result.iterations, 0);
        iterations_allocated = true;
        for (&result.aerodynamic_fields) |*values| {
            values.* = try allocator.alloc(f64, cell_count);
            @memset(values.*, 1);
            aerodynamic_allocated += 1;
        }
        result.ground_air = undefined;
        result.ground_air.cell_count = cell_count;
        result.ground_air.temperature_k = result.ground_fields[0];
        result.ground_air.vapor_volume_fraction = result.ground_fields[1];
        result.ground_air.heat_capacity_megajoules_per_k = result.ground_fields[2];
        result.ground_air.air_volume_m3 = result.ground_fields[3];
        result.ground_air.iteration_count = result.iterations;
        result.aerodynamics = undefined;
        result.aerodynamics.cell_count = cell_count;
        result.aerodynamics.zero_plane_displacement_m =
            result.aerodynamic_fields[0];
        result.aerodynamics.effective_roughness_height_m =
            result.aerodynamic_fields[1];
        result.aerodynamics.wind_reference_height_m =
            result.aerodynamic_fields[2];
        result.aerodynamics.bulk_richardson_coefficient_k =
            result.aerodynamic_fields[3];
        result.aerodynamics.isothermal_aerodynamic_resistance_h_per_m =
            result.aerodynamic_fields[4];
        @memset(result.ground_air.temperature_k, 280);
        @memset(result.ground_air.vapor_volume_fraction, 0.005);
        @memset(result.aerodynamics.zero_plane_displacement_m, 0);
        return result;
    }

    fn deinit(self: *TestSurfaceBoundary) void {
        for (self.aerodynamic_fields) |values| self.allocator.free(values);
        self.allocator.free(self.iterations);
        for (self.ground_fields) |values| self.allocator.free(values);
        self.* = undefined;
    }

    fn view(self: *const TestSurfaceBoundary) SurfaceBoundary.View {
        return .{
            .ground_air = &self.ground_air,
            .surface_aerodynamics = &self.aerodynamics,
        };
    }
};

test "soil geometry checkpoint round trips arbitrary grid layers and flux ledgers" {
    var geometry = try Geometry.init(std.testing.allocator, 6, 7);
    defer geometry.deinit();
    geometry.first_active_layer[5] = 2;
    geometry.active_layer_count[5] = 4;
    geometry.boundary_depth_m[5 * 8 + 3] = 3.25;
    geometry.layer_thickness_m[5 * 7 + 2] = 1.25;
    var hydrology = try Hydrology.init(std.testing.allocator, 3, 2, 7, 4);
    defer hydrology.deinit();
    var surface = try Surface.init(std.testing.allocator, 6);
    defer surface.deinit();
    var surface_litter_geometry =
        try SurfaceLitterGeometry.init(std.testing.allocator, 6);
    defer surface_litter_geometry.deinit();
    var erosion = try Erosion.init(std.testing.allocator, 3, 2);
    defer erosion.deinit();
    var climate: Climate = .{};
    var eroded_minerals = try ErodedMinerals.init(std.testing.allocator, 6);
    defer eroded_minerals.deinit();
    var runtime = try TestRuntime.init(std.testing.allocator, 6, 7);
    defer runtime.deinit();
    var surface_boundary =
        try TestSurfaceBoundary.init(std.testing.allocator, 6);
    defer surface_boundary.deinit();
    surface_boundary.ground_air.temperature_k[5] = 267.25;
    surface_boundary.ground_air.vapor_volume_fraction[5] = 0.0027;
    surface_boundary.aerodynamics.effective_roughness_height_m[5] = 0.13;
    runtime.properties.matrix_bulk_volume_m3[41] = 0.73;
    runtime.properties.layer_volume_m3[41] = 1.25;
    runtime.properties.layer_thickness_m[41] = 1.25;
    runtime.properties.bulk_density_megagrams_per_m3[41] = 1.37;
    runtime.properties.porosity_fraction[41] = 0.42;
    runtime.thermal.layer_volume_m3[41] = 0.44;
    runtime.thermal.layer_thickness_m[41] = 0.044;
    runtime.thermal.porosity_fraction[41] = 0.42;
    runtime.thermal.dry_solid_heat_capacity_megajoules_per_m3_k[41] = 1.67;
    runtime.thermal.total_heat_capacity_megajoules_per_m3_k[41] = 2.91;
    runtime.thermal.thermal_conductivity_m_megajoules_per_h_k[41] = 0.0042;
    hydrology.micropore_water_volume_m3[41] = 2.5;
    hydrology.heat_face_flux_megajoules_per_step[125] = -3.5;
    hydrology.snow_liquid_water_volume_m3[23] = 0.4;
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    surface.solid_snow_water_equivalent_m3[5] = 0.7;
    surface_litter_geometry.expanded_total_volume_m3[5] = 0.41;
    surface_litter_geometry.pore_volume_m3[5] = 0.30;
    surface_litter_geometry.air_volume_m3[5] = 0.17;
    surface_litter_geometry.porosity_m3_per_m3[5] = 0.63;
    surface_litter_geometry.field_capacity_m3_per_m3[5] = 0.31;
    surface_litter_geometry.wilting_point_m3_per_m3[5] = 0.12;
    erosion.surface_sediment_megagrams[5] = 0.25;
    erosion.surface_soil_mass_megagrams[5] = 4.5;
    erosion.surface_soil_mass_initialized[5] = true;
    climate.modifiers[2].precipitation = 1.25;
    eroded_minerals.initialized = true;
    eroded_minerals.workspace.pools[5] = 2.75;
    var live_heat = [_]f64{0} ** 12;
    var dead_heat = [_]f64{0} ** 12;
    var subsurface_heat = [_]f64{0} ** 42;
    var surface_heat = [_]f64{0} ** 6;
    var surface_ice = [_]f64{0} ** 6;
    live_heat[11] = 1.25;
    dead_heat[10] = 2.5;
    subsurface_heat[41] = 3.75;
    surface_heat[5] = 5;
    surface_ice[5] = 0.375;
    try write(&bytes.writer, .{ .geometry = &geometry, .hydrology = &hydrology, .surface = &surface, .erosion = &erosion, .climate = &climate, .eroded_minerals = &eroded_minerals, .runtime = runtime.view(), .surface_boundary = surface_boundary.view(), .surface_litter_geometry = &surface_litter_geometry, .surface_litter_ice_m3 = &surface_ice, .delayed_live_canopy_combustion_heat_megajoules = &live_heat, .delayed_standing_dead_combustion_heat_megajoules = &dead_heat, .delayed_subsurface_combustion_heat_megajoules = &subsurface_heat, .delayed_surface_combustion_heat_megajoules = &surface_heat });
    var reader: std.Io.Reader = .fixed(bytes.written());
    var restored = try read(std.testing.allocator, &reader, .{ .maximum_columns = 10, .maximum_rows = 10, .maximum_soil_layers = 20, .maximum_snow_layers = 10, .maximum_plants = 20 });
    defer restored.deinit();
    try std.testing.expectEqualSlices(usize, geometry.first_active_layer, restored.geometry.first_active_layer);
    try std.testing.expectEqualSlices(f64, geometry.boundary_depth_m, restored.geometry.boundary_depth_m);
    try std.testing.expectEqualSlices(f64, hydrology.heat_face_flux_megajoules_per_step, restored.hydrology.heat_face_flux_megajoules_per_step);
    try std.testing.expectEqualSlices(f64, hydrology.snow_liquid_water_volume_m3, restored.hydrology.snow_liquid_water_volume_m3);
    try std.testing.expectEqualSlices(f64, surface.solid_snow_water_equivalent_m3, restored.surface.solid_snow_water_equivalent_m3);
    try std.testing.expectEqualSlices(f64, surface_litter_geometry.expanded_total_volume_m3, restored.surface_litter_geometry.expanded_total_volume_m3);
    try std.testing.expectEqualSlices(f64, surface_litter_geometry.air_volume_m3, restored.surface_litter_geometry.air_volume_m3);
    try std.testing.expectEqualSlices(f64, erosion.surface_sediment_megagrams, restored.erosion.surface_sediment_megagrams);
    try std.testing.expectEqualSlices(f64, erosion.surface_soil_mass_megagrams, restored.erosion.surface_soil_mass_megagrams);
    try std.testing.expectEqualSlices(bool, erosion.surface_soil_mass_initialized, restored.erosion.surface_soil_mass_initialized);
    try std.testing.expectEqual(climate.modifiers[2], restored.climate.modifiers[2]);
    try std.testing.expectEqualSlices(f64, eroded_minerals.workspace.pools, restored.eroded_minerals.workspace.pools);
    try std.testing.expectEqualSlices(f64, &live_heat, restored.delayed_live_canopy_combustion_heat_megajoules);
    try std.testing.expectEqualSlices(f64, &dead_heat, restored.delayed_standing_dead_combustion_heat_megajoules);
    try std.testing.expectEqualSlices(f64, &subsurface_heat, restored.delayed_subsurface_combustion_heat_megajoules);
    try std.testing.expectEqualSlices(f64, &surface_heat, restored.delayed_surface_combustion_heat_megajoules);
    try std.testing.expectEqualSlices(f64, &surface_ice, restored.surface_litter_ice_m3);
    try std.testing.expectEqual(
        surface_boundary.ground_air.temperature_k[5],
        restored.surface_boundary.ground_air_fields[0][5],
    );
    try std.testing.expectEqual(
        surface_boundary.aerodynamics.effective_roughness_height_m[5],
        restored.surface_boundary.aerodynamic_fields[1][5],
    );
    var target_surface_boundary =
        try TestSurfaceBoundary.init(std.testing.allocator, 6);
    defer target_surface_boundary.deinit();
    @memset(target_surface_boundary.aerodynamics.wind_reference_height_m, 0);
    @memset(
        target_surface_boundary.aerodynamics
            .isothermal_aerodynamic_resistance_h_per_m,
        0,
    );
    try restored.surface_boundary.restoreInto(
        &target_surface_boundary.ground_air,
        &target_surface_boundary.aerodynamics,
    );
    try std.testing.expectEqual(
        surface_boundary.ground_air.temperature_k[5],
        target_surface_boundary.ground_air.temperature_k[5],
    );
    try std.testing.expectEqual(
        surface_boundary.aerodynamics.effective_roughness_height_m[5],
        target_surface_boundary.aerodynamics.effective_roughness_height_m[5],
    );
    try std.testing.expectEqual(
        runtime.properties.matrix_bulk_volume_m3[41],
        restored.runtime.solver_fields[
            @intFromEnum(Runtime.SolverField.matrix_bulk_volume_m3)
        ][41],
    );
    try std.testing.expectEqual(
        runtime.thermal.dry_solid_heat_capacity_megajoules_per_m3_k[41],
        restored.runtime.thermal_fields[
            @intFromEnum(
                Runtime.ThermalField.dry_solid_heat_capacity_megajoules_per_m3_k,
            )
        ][41],
    );
    var target_runtime = try TestRuntime.init(std.testing.allocator, 6, 7);
    defer target_runtime.deinit();
    try restored.runtime.restoreInto(
        &target_runtime.properties,
        &target_runtime.thermal,
    );
    try std.testing.expectEqual(
        runtime.properties.bulk_density_megagrams_per_m3[41],
        target_runtime.properties.bulk_density_megagrams_per_m3[41],
    );
    try std.testing.expectEqual(
        runtime.thermal.total_heat_capacity_megajoules_per_m3_k[41],
        target_runtime.thermal.total_heat_capacity_megajoules_per_m3_k[41],
    );
    try std.testing.expectEqual(
        runtime.thermal.layer_thickness_m[41],
        target_runtime.thermal.layer_thickness_m[41],
    );
    try std.testing.expect(
        target_runtime.thermal.layer_thickness_m[41] !=
            geometry.layer_thickness_m[41],
    );
}

test "soil geometry checkpoint applies limits before allocation" {
    var geometry = try Geometry.init(std.testing.allocator, 6, 7);
    defer geometry.deinit();
    var hydrology = try Hydrology.init(std.testing.allocator, 3, 2, 7, 4);
    defer hydrology.deinit();
    var surface = try Surface.init(std.testing.allocator, 6);
    defer surface.deinit();
    var surface_litter_geometry =
        try SurfaceLitterGeometry.init(std.testing.allocator, 6);
    defer surface_litter_geometry.deinit();
    var erosion = try Erosion.init(std.testing.allocator, 3, 2);
    defer erosion.deinit();
    var climate: Climate = .{};
    var eroded_minerals = try ErodedMinerals.init(std.testing.allocator, 6);
    defer eroded_minerals.deinit();
    var runtime = try TestRuntime.init(std.testing.allocator, 6, 7);
    defer runtime.deinit();
    var surface_boundary =
        try TestSurfaceBoundary.init(std.testing.allocator, 6);
    defer surface_boundary.deinit();
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    const plant_heat = [_]f64{0} ** 12;
    const subsurface_heat = [_]f64{0} ** 42;
    const surface_heat = [_]f64{0} ** 6;
    try write(&bytes.writer, .{ .geometry = &geometry, .hydrology = &hydrology, .surface = &surface, .erosion = &erosion, .climate = &climate, .eroded_minerals = &eroded_minerals, .runtime = runtime.view(), .surface_boundary = surface_boundary.view(), .surface_litter_geometry = &surface_litter_geometry, .surface_litter_ice_m3 = &surface_heat, .delayed_live_canopy_combustion_heat_megajoules = &plant_heat, .delayed_standing_dead_combustion_heat_megajoules = &plant_heat, .delayed_subsurface_combustion_heat_megajoules = &subsurface_heat, .delayed_surface_combustion_heat_megajoules = &surface_heat });
    var reader: std.Io.Reader = .fixed(bytes.written());
    try std.testing.expectError(error.SoilGeometryCheckpointSoilLayerLimitExceeded, read(std.testing.allocator, &reader, .{ .maximum_columns = 3, .maximum_rows = 2, .maximum_soil_layers = 6, .maximum_snow_layers = 4, .maximum_plants = 12 }));
}

test "soil geometry checkpoint rejects corruption and nonfinite state" {
    var geometry = try Geometry.init(std.testing.allocator, 1, 1);
    defer geometry.deinit();
    var hydrology = try Hydrology.init(std.testing.allocator, 1, 1, 1, 1);
    defer hydrology.deinit();
    var surface = try Surface.init(std.testing.allocator, 1);
    defer surface.deinit();
    var surface_litter_geometry =
        try SurfaceLitterGeometry.init(std.testing.allocator, 1);
    defer surface_litter_geometry.deinit();
    var erosion = try Erosion.init(std.testing.allocator, 1, 1);
    defer erosion.deinit();
    var climate: Climate = .{};
    var eroded_minerals = try ErodedMinerals.init(std.testing.allocator, 1);
    defer eroded_minerals.deinit();
    var runtime = try TestRuntime.init(std.testing.allocator, 1, 1);
    defer runtime.deinit();
    var surface_boundary =
        try TestSurfaceBoundary.init(std.testing.allocator, 1);
    defer surface_boundary.deinit();
    geometry.boundary_depth_m[0] = std.math.nan(f64);
    var invalid: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer invalid.deinit();
    const heat = [_]f64{0};
    try std.testing.expectError(error.NonFiniteSoilGeometryCheckpoint, write(&invalid.writer, .{ .geometry = &geometry, .hydrology = &hydrology, .surface = &surface, .erosion = &erosion, .climate = &climate, .eroded_minerals = &eroded_minerals, .runtime = runtime.view(), .surface_boundary = surface_boundary.view(), .surface_litter_geometry = &surface_litter_geometry, .surface_litter_ice_m3 = &heat, .delayed_live_canopy_combustion_heat_megajoules = &heat, .delayed_standing_dead_combustion_heat_megajoules = &heat, .delayed_subsurface_combustion_heat_megajoules = &heat, .delayed_surface_combustion_heat_megajoules = &heat }));
    geometry.boundary_depth_m[0] = 0;
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try write(&bytes.writer, .{ .geometry = &geometry, .hydrology = &hydrology, .surface = &surface, .erosion = &erosion, .climate = &climate, .eroded_minerals = &eroded_minerals, .runtime = runtime.view(), .surface_boundary = surface_boundary.view(), .surface_litter_geometry = &surface_litter_geometry, .surface_litter_ice_m3 = &heat, .delayed_live_canopy_combustion_heat_megajoules = &heat, .delayed_standing_dead_combustion_heat_megajoules = &heat, .delayed_subsurface_combustion_heat_megajoules = &heat, .delayed_surface_combustion_heat_megajoules = &heat });
    try bytes.writer.writeByte(0xff);
    var reader: std.Io.Reader = .fixed(bytes.written());
    try std.testing.expectError(error.TrailingSoilGeometryCheckpointData, read(std.testing.allocator, &reader, .{ .maximum_columns = 1, .maximum_rows = 1, .maximum_soil_layers = 1, .maximum_snow_layers = 1, .maximum_plants = 1 }));
}
