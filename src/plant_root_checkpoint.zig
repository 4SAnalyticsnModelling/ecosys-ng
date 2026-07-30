const std = @import("std");
const RootState = @import("plant_root_system.zig").State;
const magic = "ECOSROOT";
const version: u32 = 3;

pub const Limits = struct { maximum_plants: usize, maximum_soil_layers: usize, maximum_root_axes: usize };

pub fn write(writer: anytype, state: RootState) !void {
    try validate(state);
    try writer.writeAll(magic);
    try writer.writeInt(u32, version, .little);
    try writer.writeInt(u64, @intCast(state.plant_count), .little);
    try writer.writeInt(u64, @intCast(state.soil_layer_count), .little);
    try writer.writeInt(u64, @intCast(state.root_axis_count), .little);
    inline for (@typeInfo(RootState).@"struct".fields) |field| switch (field.type) {
        []f64 => try writeF64Slice(writer, @field(state, field.name)),
        []usize => try writeUsizeSlice(writer, @field(state, field.name)),
        []bool => try writeBoolSlice(writer, @field(state, field.name)),
        else => {},
    };
}

pub fn read(allocator: std.mem.Allocator, reader: *std.Io.Reader, limits: Limits) !RootState {
    if (limits.maximum_plants == 0 or limits.maximum_soil_layers == 0 or limits.maximum_root_axes == 0) return error.InvalidPlantRootCheckpointLimits;
    if (!std.mem.eql(u8, try reader.takeArray(magic.len), magic)) return error.InvalidPlantRootCheckpointMagic;
    if (try reader.takeInt(u32, .little) != version) return error.UnsupportedPlantRootCheckpointVersion;
    const plants = try bounded(reader, limits.maximum_plants, error.PlantRootCheckpointPlantLimitExceeded);
    const layers = try bounded(reader, limits.maximum_soil_layers, error.PlantRootCheckpointLayerLimitExceeded);
    const axes = try bounded(reader, limits.maximum_root_axes, error.PlantRootCheckpointAxisLimitExceeded);
    if (plants == 0 or layers == 0 or axes == 0) return error.InvalidPlantRootCheckpointDimensions;
    var state = try RootState.init(allocator, plants, layers, axes);
    errdefer state.deinit();
    inline for (@typeInfo(RootState).@"struct".fields) |field| switch (field.type) {
        []f64 => try readF64Slice(reader, @field(state, field.name)),
        []usize => try readUsizeSlice(reader, @field(state, field.name)),
        []bool => try readBoolSlice(reader, @field(state, field.name)),
        else => {},
    };
    if (reader.peekByte()) |_| return error.TrailingPlantRootCheckpointData else |err| switch (err) {
        error.EndOfStream => {},
        else => return err,
    }
    try validate(state);
    return state;
}

fn validate(state: RootState) !void {
    if (state.plant_count == 0 or state.soil_layer_count == 0 or state.root_axis_count == 0 or state.planting_layer_by_plant.len != state.plant_count or state.active_root_axis_count.len != state.plant_count or state.roots_dead.len != state.plant_count) return error.InvalidPlantRootCheckpointDimensions;
    for (state.planting_layer_by_plant) |layer| if (layer >= state.soil_layer_count) return error.InvalidCheckpointPlantingLayer;
    for (state.active_root_axis_count) |count| if (count > state.root_axis_count) return error.InvalidCheckpointActiveRootAxisCount;
    try state.validateFinite();
}
fn bounded(reader: *std.Io.Reader, limit: usize, comptime too_large: anyerror) !usize {
    const value = try reader.takeInt(u64, .little);
    if (value > limit or value > std.math.maxInt(usize)) return too_large;
    return @intCast(value);
}
fn writeF64Slice(writer: anytype, values: []const f64) !void {
    for (values) |value| {
        if (!std.math.isFinite(value)) return error.NonFinitePlantRootCheckpoint;
        try writer.writeInt(u64, @bitCast(value), .little);
    }
}
fn readF64Slice(reader: *std.Io.Reader, values: []f64) !void {
    for (values) |*value| {
        value.* = @bitCast(try reader.takeInt(u64, .little));
        if (!std.math.isFinite(value.*)) return error.NonFinitePlantRootCheckpoint;
    }
}
fn writeUsizeSlice(writer: anytype, values: []const usize) !void {
    for (values) |value| try writer.writeInt(u64, @intCast(value), .little);
}
fn readUsizeSlice(reader: *std.Io.Reader, values: []usize) !void {
    for (values) |*value| {
        const stored = try reader.takeInt(u64, .little);
        if (stored > std.math.maxInt(usize)) return error.PlantRootCheckpointIntegerOverflow;
        value.* = @intCast(stored);
    }
}
fn writeBoolSlice(writer: anytype, values: []const bool) !void {
    for (values) |value| try writer.writeByte(@intFromBool(value));
}
fn readBoolSlice(reader: *std.Io.Reader, values: []bool) !void {
    for (values) |*value| value.* = switch (try reader.takeByte()) {
        0 => false,
        1 => true,
        else => return error.InvalidPlantRootCheckpointBoolean,
    };
}

test "root checkpoint round trip preserves runtime plants layers axes and every field" {
    var source = try RootState.init(std.testing.allocator, 7, 4, 12);
    defer source.deinit();
    source.planting_layer_by_plant[6] = 3;
    source.active_root_axis_count[6] = 11;
    source.roots_dead[6] = false;
    source.total_carbon_g[source.total_carbon_g.len - 1] = 42;
    source.axis_secondary_phosphorus_g[source.axis_secondary_phosphorus_g.len - 1] = 0.25;
    source.salt_content_mol[source.salt_content_mol.len - 1] = 3.5;
    source.exudate_carbon_exchange_g_c_per_h[source.exudate_carbon_exchange_g_c_per_h.len - 1] = -0.1;
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try write(&bytes.writer, source);
    var reader: std.Io.Reader = .fixed(bytes.written());
    var restored = try read(std.testing.allocator, &reader, .{ .maximum_plants = 20, .maximum_soil_layers = 30, .maximum_root_axes = 40 });
    defer restored.deinit();
    try std.testing.expectEqual(@as(usize, 7), restored.plant_count);
    try std.testing.expectEqual(@as(usize, 12), restored.root_axis_count);
    try std.testing.expectEqual(@as(usize, 11), restored.active_root_axis_count[6]);
    try std.testing.expect(!restored.roots_dead[6]);
    try std.testing.expectEqualSlices(f64, source.total_carbon_g, restored.total_carbon_g);
    try std.testing.expectEqualSlices(f64, source.axis_secondary_phosphorus_g, restored.axis_secondary_phosphorus_g);
    try std.testing.expectEqualSlices(f64, source.salt_content_mol, restored.salt_content_mol);
    try std.testing.expectEqualSlices(f64, source.exudate_carbon_exchange_g_c_per_h, restored.exudate_carbon_exchange_g_c_per_h);
}

test "root checkpoint enforces runtime axis limit before allocation" {
    var source = try RootState.init(std.testing.allocator, 1, 2, 12);
    defer source.deinit();
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try write(&bytes.writer, source);
    var reader: std.Io.Reader = .fixed(bytes.written());
    try std.testing.expectError(error.PlantRootCheckpointAxisLimitExceeded, read(std.testing.allocator, &reader, .{ .maximum_plants = 1, .maximum_soil_layers = 2, .maximum_root_axes = 10 }));
}
