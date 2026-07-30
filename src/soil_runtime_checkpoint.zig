const std = @import("std");
const SoilProperties = @import("soil_solver_properties.zig").State;
const SoilThermal = @import("soil_thermal.zig").State;

pub const solver_field_count = @typeInfo(SolverField).@"enum".fields.len;
pub const thermal_field_count = @typeInfo(ThermalField).@"enum".fields.len;

pub const SolverField = enum {
    matrix_bulk_volume_m3,
    layer_volume_m3,
    layer_thickness_m,
    layer_midpoint_depth_m,
    layer_bottom_depth_m,
    bulk_density_megagrams_per_m3,
    sand_mass_fraction,
    clay_mass_fraction,
    sand_mass_Mg,
    silt_mass_Mg,
    clay_mass_Mg,
    total_organic_carbon_g_per_megagram,
    cation_exchange_capacity_mol_per_Mg,
    anion_exchange_capacity_mol_per_Mg,
    cation_exchange_capacity_mol,
    anion_exchange_capacity_mol,
    porosity_fraction,
    rainfall_conductivity_multiplier,
};

pub const ThermalField = enum {
    layer_volume_m3,
    layer_thickness_m,
    porosity_fraction,
    dry_solid_heat_capacity_mj_per_m3_k,
    solid_thermal_conductivity_numerator_m_mj_per_h_k,
    solid_thermal_conductivity_denominator,
    total_heat_capacity_mj_per_m3_k,
    thermal_conductivity_m_mj_per_h_k,
};

pub const View = struct {
    soil_properties: *const SoilProperties,
    soil_thermal: *const SoilThermal,
};

pub const Snapshot = struct {
    allocator: std.mem.Allocator,
    layer_count: usize,
    solver_fields: [solver_field_count][]f64,
    thermal_fields: [thermal_field_count][]f64,

    pub fn deinit(self: *Snapshot) void {
        for (self.thermal_fields) |values| self.allocator.free(values);
        for (self.solver_fields) |values| self.allocator.free(values);
        self.* = undefined;
    }

    pub fn restoreInto(
        self: Snapshot,
        soil_properties: *SoilProperties,
        soil_thermal: *SoilThermal,
    ) !void {
        try validateDimensions(
            .{ .soil_properties = soil_properties, .soil_thermal = soil_thermal },
            self.layer_count,
        );
        try self.validate();
        inline for (@typeInfo(SolverField).@"enum".fields, 0..) |field, index|
            @memcpy(
                mutableSolverValues(soil_properties, @enumFromInt(field.value)),
                self.solver_fields[index],
            );
        inline for (@typeInfo(ThermalField).@"enum".fields, 0..) |field, index|
            @memcpy(
                mutableThermalValues(soil_thermal, @enumFromInt(field.value)),
                self.thermal_fields[index],
            );
        try validateView(.{
            .soil_properties = soil_properties,
            .soil_thermal = soil_thermal,
        });
    }

    pub fn validate(self: Snapshot) !void {
        if (self.layer_count == 0)
            return error.SoilRuntimeCheckpointDimensionMismatch;
        for (self.solver_fields) |values| {
            if (values.len != self.layer_count)
                return error.SoilRuntimeCheckpointDimensionMismatch;
            try validateFinite(values);
        }
        for (self.thermal_fields) |values| {
            if (values.len != self.layer_count)
                return error.SoilRuntimeCheckpointDimensionMismatch;
            try validateFinite(values);
        }
        const solver_volume = self.solver_fields[
            @intFromEnum(SolverField.layer_volume_m3)
        ];
        const matrix_volume = self.solver_fields[
            @intFromEnum(SolverField.matrix_bulk_volume_m3)
        ];
        const density = self.solver_fields[
            @intFromEnum(SolverField.bulk_density_megagrams_per_m3)
        ];
        const solver_porosity = self.solver_fields[
            @intFromEnum(SolverField.porosity_fraction)
        ];
        const solver_thickness = self.solver_fields[
            @intFromEnum(SolverField.layer_thickness_m)
        ];
        const thermal_volume = self.thermal_fields[
            @intFromEnum(ThermalField.layer_volume_m3)
        ];
        const thermal_thickness = self.thermal_fields[
            @intFromEnum(ThermalField.layer_thickness_m)
        ];
        const thermal_porosity = self.thermal_fields[
            @intFromEnum(ThermalField.porosity_fraction)
        ];
        const dry_heat_capacity = self.thermal_fields[
            @intFromEnum(ThermalField.dry_solid_heat_capacity_mj_per_m3_k)
        ];
        const total_heat_capacity = self.thermal_fields[
            @intFromEnum(ThermalField.total_heat_capacity_mj_per_m3_k)
        ];
        const conductivity = self.thermal_fields[
            @intFromEnum(ThermalField.thermal_conductivity_m_mj_per_h_k)
        ];
        for (0..self.layer_count) |index| {
            if (solver_volume[index] <= 0 or matrix_volume[index] < 0 or
                matrix_volume[index] > solver_volume[index] or
                density[index] < 0 or solver_porosity[index] < 0 or
                solver_porosity[index] > 1 or solver_thickness[index] <= 0 or
                thermal_volume[index] <= 0 or thermal_thickness[index] <= 0 or
                thermal_porosity[index] < 0 or thermal_porosity[index] > 1 or
                dry_heat_capacity[index] < 0 or total_heat_capacity[index] < 0 or
                conductivity[index] < 0)
                return error.InvalidSoilRuntimeCheckpointState;
        }
    }
};

pub fn write(writer: anytype, view: View) !void {
    try validateView(view);
    try writer.writeInt(u64, @intCast(view.soil_properties.layer_count), .little);
    inline for (@typeInfo(SolverField).@"enum".fields) |field|
        try writeF64Slice(
            writer,
            solverValues(view.soil_properties, @enumFromInt(field.value)),
        );
    inline for (@typeInfo(ThermalField).@"enum".fields) |field|
        try writeF64Slice(
            writer,
            thermalValues(view.soil_thermal, @enumFromInt(field.value)),
        );
}

pub fn read(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    expected_layer_count: usize,
) !Snapshot {
    const layer_count_u64 = try reader.takeInt(u64, .little);
    if (layer_count_u64 != expected_layer_count)
        return error.SoilRuntimeCheckpointDimensionMismatch;
    var result: Snapshot = .{
        .allocator = allocator,
        .layer_count = expected_layer_count,
        .solver_fields = undefined,
        .thermal_fields = undefined,
    };
    var solver_allocated: usize = 0;
    var thermal_allocated: usize = 0;
    errdefer {
        for (result.thermal_fields[0..thermal_allocated]) |values|
            allocator.free(values);
        for (result.solver_fields[0..solver_allocated]) |values|
            allocator.free(values);
    }
    for (&result.solver_fields) |*values| {
        values.* = try allocator.alloc(f64, expected_layer_count);
        solver_allocated += 1;
        try readF64Slice(reader, values.*);
    }
    for (&result.thermal_fields) |*values| {
        values.* = try allocator.alloc(f64, expected_layer_count);
        thermal_allocated += 1;
        try readF64Slice(reader, values.*);
    }
    try result.validate();
    return result;
}

pub fn validateView(view: View) !void {
    const layer_count = view.soil_properties.layer_count;
    try validateDimensions(view, layer_count);
    inline for (@typeInfo(SolverField).@"enum".fields) |field|
        try validateFinite(
            solverValues(view.soil_properties, @enumFromInt(field.value)),
        );
    inline for (@typeInfo(ThermalField).@"enum".fields) |field|
        try validateFinite(
            thermalValues(view.soil_thermal, @enumFromInt(field.value)),
        );
    for (view.soil_properties.layer_volume_m3, 0..) |volume, index| {
        const matrix_volume = view.soil_properties.matrix_bulk_volume_m3[index];
        const density = view.soil_properties.bulk_density_megagrams_per_m3[index];
        const porosity = view.soil_properties.porosity_fraction[index];
        if (volume <= 0 or matrix_volume < 0 or matrix_volume > volume or
            density < 0 or porosity < 0 or porosity > 1 or
            view.soil_properties.layer_thickness_m[index] <= 0 or
            view.soil_thermal.layer_volume_m3[index] <= 0 or
            view.soil_thermal.layer_thickness_m[index] <= 0 or
            view.soil_thermal.porosity_fraction[index] < 0 or
            view.soil_thermal.porosity_fraction[index] > 1 or
            view.soil_thermal.dry_solid_heat_capacity_mj_per_m3_k[index] < 0 or
            view.soil_thermal.total_heat_capacity_mj_per_m3_k[index] < 0 or
            view.soil_thermal.thermal_conductivity_m_mj_per_h_k[index] < 0)
            return error.InvalidSoilRuntimeCheckpointState;
    }
}

fn validateDimensions(view: View, layer_count: usize) !void {
    if (layer_count == 0 or
        view.soil_thermal.cell_count == 0 or
        view.soil_thermal.soil_layer_capacity == 0 or
        view.soil_thermal.cell_count * view.soil_thermal.soil_layer_capacity !=
            layer_count)
        return error.SoilRuntimeCheckpointDimensionMismatch;
    inline for (@typeInfo(SolverField).@"enum".fields) |field|
        if (solverValues(view.soil_properties, @enumFromInt(field.value)).len !=
            layer_count)
            return error.SoilRuntimeCheckpointDimensionMismatch;
    inline for (@typeInfo(ThermalField).@"enum".fields) |field|
        if (thermalValues(view.soil_thermal, @enumFromInt(field.value)).len !=
            layer_count)
            return error.SoilRuntimeCheckpointDimensionMismatch;
}

fn solverValues(state: *const SoilProperties, field: SolverField) []const f64 {
    return switch (field) {
        inline else => |selected| @field(state, @tagName(selected)),
    };
}

fn mutableSolverValues(state: *SoilProperties, field: SolverField) []f64 {
    return switch (field) {
        inline else => |selected| @field(state, @tagName(selected)),
    };
}

fn thermalValues(state: *const SoilThermal, field: ThermalField) []const f64 {
    return switch (field) {
        inline else => |selected| @field(state, @tagName(selected)),
    };
}

fn mutableThermalValues(state: *SoilThermal, field: ThermalField) []f64 {
    return switch (field) {
        inline else => |selected| @field(state, @tagName(selected)),
    };
}

fn writeF64Slice(writer: anytype, values: []const f64) !void {
    for (values) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSoilRuntimeCheckpoint;
        try writer.writeInt(u64, @bitCast(value), .little);
    }
}

fn readF64Slice(reader: *std.Io.Reader, values: []f64) !void {
    for (values) |*value| {
        value.* = @bitCast(try reader.takeInt(u64, .little));
        if (!std.math.isFinite(value.*))
            return error.NonFiniteSoilRuntimeCheckpoint;
    }
}

fn validateFinite(values: []const f64) !void {
    for (values) |value| if (!std.math.isFinite(value))
        return error.NonFiniteSoilRuntimeCheckpoint;
}
