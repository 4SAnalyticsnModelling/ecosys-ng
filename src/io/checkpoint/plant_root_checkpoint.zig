const std = @import("std");
const RootState = @import("../../plant/root/plant_root_system.zig").State;
const magic = "ECOSROOT";
/// Version 7 adds `retained_root_carbon_g_c_per_plant` (GROSUB 507 `WTRTA`),
/// the persistent per-plant recurrence that `BIND-GROSUB-506` introduced.
/// Because `WTRTA` is a recurrence on its own previous value, it MUST be
/// serialized or a checkpoint-resume run diverges from a continuous one, which
/// is exactly the Wave 2 `RESTART-EQUIVALENCE` obligation.
const version: u32 = 7;
const without_retained_root_carbon_version: u32 = 6;
const previous_version: u32 = 5;
const compressed_porosity_version: u32 = 4;

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
    const stored_version = try reader.takeInt(u32, .little);
    if (stored_version != version and
        stored_version != without_retained_root_carbon_version and
        stored_version != previous_version and
        stored_version != compressed_porosity_version)
        return error.UnsupportedPlantRootCheckpointVersion;
    const plants = try bounded(reader, limits.maximum_plants, error.PlantRootCheckpointPlantLimitExceeded);
    const layers = try bounded(reader, limits.maximum_soil_layers, error.PlantRootCheckpointLayerLimitExceeded);
    const axes = try bounded(reader, limits.maximum_root_axes, error.PlantRootCheckpointAxisLimitExceeded);
    if (plants == 0 or layers == 0 or axes == 0) return error.InvalidPlantRootCheckpointDimensions;
    var state = try RootState.init(allocator, plants, layers, axes);
    errdefer state.deinit();
    inline for (@typeInfo(RootState).@"struct".fields) |field| switch (field.type) {
        []f64 => if (comptime std.mem.eql(u8, field.name, "retained_root_carbon_g_c_per_plant")) {
            // Absent before version 7. A pre-7 checkpoint has no `WTRTA`
            // history to restore, so it stays at the `State.init` zero and the
            // recurrence re-seeds itself on the first resumed hour from the
            // `WTRT/PP` branch of GROSUB 507--508, which is the same branch a
            // fresh planting takes. That is a defined restart, not silent
            // corruption, and it is stated here so nobody reads a pre-7 resume
            // as evidence for or against restart equivalence.
            if (stored_version >= version) try readF64Slice(reader, @field(state, field.name));
        } else if (stored_version == compressed_porosity_version and
            (comptime std.mem.eql(u8, field.name, "current_porosity_fraction_by_domain") or
                std.mem.eql(u8, field.name, "initial_porosity_fraction_by_domain")))
        {
            const per_plant = try allocator.alloc(f64, plants);
            defer allocator.free(per_plant);
            try readF64Slice(reader, per_plant);
            for (per_plant, 0..) |value, plant| {
                for (0..@import("../../plant/root/plant_root_system.zig").biological_domain_count) |domain| {
                    @field(state, field.name)[try state.domainIndex(plant, domain)] = value;
                }
            }
        } else try readF64Slice(reader, @field(state, field.name)),
        []usize => if (stored_version < version and (comptime isRootedLayerBound(field.name))) {
            @memcpy(@field(state, field.name), state.planting_layer_by_plant);
        } else try readUsizeSlice(reader, @field(state, field.name)),
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
    source.current_deepest_rooted_layer_by_plant[6] = 3;
    source.next_deepest_rooted_layer_by_plant[6] = 3;
    try source.includeNextDeepestRootedLayer(5, 2);
    source.active_root_axis_count[6] = 11;
    source.roots_dead[6] = false;
    source.total_carbon_g[source.total_carbon_g.len - 1] = 42;
    source.axis_secondary_phosphorus_g[source.axis_secondary_phosphorus_g.len - 1] = 0.25;
    source.salt_content_mol[source.salt_content_mol.len - 1] = 3.5;
    source.exudate_carbon_exchange_g_c_per_h[source.exudate_carbon_exchange_g_c_per_h.len - 1] = -0.1;
    source.current_porosity_fraction_by_domain[try source.domainIndex(6, 0)] = 0.21;
    source.current_porosity_fraction_by_domain[try source.domainIndex(6, 1)] = 0.47;
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
    try std.testing.expectEqualSlices(f64, source.current_porosity_fraction_by_domain, restored.current_porosity_fraction_by_domain);
    try std.testing.expectEqualSlices(usize, source.current_deepest_rooted_layer_by_plant, restored.current_deepest_rooted_layer_by_plant);
    try std.testing.expectEqualSlices(usize, source.next_deepest_rooted_layer_by_plant, restored.next_deepest_rooted_layer_by_plant);
}

fn isRootedLayerBound(comptime name: []const u8) bool {
    return std.mem.eql(u8, name, "current_deepest_rooted_layer_by_plant") or
        std.mem.eql(u8, name, "next_deepest_rooted_layer_by_plant");
}

fn writePriorVersionForTest(writer: anytype, state: RootState, stored_version: u32) !void {
    try validate(state);
    try writer.writeAll(magic);
    if (stored_version != compressed_porosity_version and stored_version != previous_version)
        return error.UnsupportedPlantRootCheckpointVersion;
    try writer.writeInt(u32, stored_version, .little);
    try writer.writeInt(u64, @intCast(state.plant_count), .little);
    try writer.writeInt(u64, @intCast(state.soil_layer_count), .little);
    try writer.writeInt(u64, @intCast(state.root_axis_count), .little);
    inline for (@typeInfo(RootState).@"struct".fields) |field| switch (field.type) {
        []f64 => if (comptime std.mem.eql(u8, field.name, "retained_root_carbon_g_c_per_plant")) {
            // Version 7 field. A pre-7 stream does not contain it, so this
            // fixture writer must not emit it either or the reflective reader
            // walk desynchronises and the read reports
            // `TrailingPlantRootCheckpointData` on a stream that is actually
            // well formed. Found by lane A7b against the in-flight tree, and
            // it is the correct failure: these two fixtures are precisely the
            // regression guard for reader/writer field-set drift, so they
            // caught a real asymmetry rather than a spurious one.
        } else if (stored_version == compressed_porosity_version and
            (comptime std.mem.eql(u8, field.name, "current_porosity_fraction_by_domain") or
                std.mem.eql(u8, field.name, "initial_porosity_fraction_by_domain")))
        {
            for (0..state.plant_count) |plant|
                try writeF64Slice(writer, @field(state, field.name)[try state.domainIndex(plant, 0)..][0..1]);
        } else try writeF64Slice(writer, @field(state, field.name)),
        []usize => if (!(comptime isRootedLayerBound(field.name)))
            try writeUsizeSlice(writer, @field(state, field.name)),
        []bool => try writeBoolSlice(writer, @field(state, field.name)),
        else => {},
    };
}

test "version four root checkpoint expands plant porosity into every biological domain" {
    var source = try RootState.init(std.testing.allocator, 2, 2, 3);
    defer source.deinit();
    source.current_porosity_fraction_by_domain[try source.domainIndex(0, 0)] = 0.23;
    source.initial_porosity_fraction_by_domain[try source.domainIndex(0, 0)] = 0.19;
    source.current_porosity_fraction_by_domain[try source.domainIndex(0, 1)] = 0.61;
    source.initial_porosity_fraction_by_domain[try source.domainIndex(0, 1)] = 0.62;
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try writePriorVersionForTest(&bytes.writer, source, compressed_porosity_version);
    var reader: std.Io.Reader = .fixed(bytes.written());
    var restored = try read(std.testing.allocator, &reader, .{ .maximum_plants = 2, .maximum_soil_layers = 2, .maximum_root_axes = 3 });
    defer restored.deinit();
    try std.testing.expectEqual(@as(f64, 0.23), restored.current_porosity_fraction_by_domain[try restored.domainIndex(0, 0)]);
    try std.testing.expectEqual(@as(f64, 0.23), restored.current_porosity_fraction_by_domain[try restored.domainIndex(0, 1)]);
    try std.testing.expectEqual(@as(f64, 0.19), restored.initial_porosity_fraction_by_domain[try restored.domainIndex(0, 1)]);
}

test "version five root checkpoint initializes NI and NIX from planting layer" {
    var source = try RootState.init(std.testing.allocator, 2, 4, 3);
    defer source.deinit();
    source.planting_layer_by_plant[0] = 1;
    source.planting_layer_by_plant[1] = 2;
    source.current_deepest_rooted_layer_by_plant[0] = 3;
    source.next_deepest_rooted_layer_by_plant[0] = 3;
    source.current_deepest_rooted_layer_by_plant[1] = 3;
    source.next_deepest_rooted_layer_by_plant[1] = 3;
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try writePriorVersionForTest(&bytes.writer, source, previous_version);
    var reader: std.Io.Reader = .fixed(bytes.written());
    var restored = try read(std.testing.allocator, &reader, .{ .maximum_plants = 2, .maximum_soil_layers = 4, .maximum_root_axes = 3 });
    defer restored.deinit();
    try std.testing.expectEqualSlices(usize, source.planting_layer_by_plant, restored.current_deepest_rooted_layer_by_plant);
    try std.testing.expectEqualSlices(usize, source.planting_layer_by_plant, restored.next_deepest_rooted_layer_by_plant);
}

test "BIND-GROSUB-506 version seven round trips the WTRTA recurrence and pre-seven defaults it" {
    // Restart-equivalence evidence for the persistent `WTRTA` field added with
    // `BIND-GROSUB-506`. `WTRTA` is a recurrence on its OWN previous value
    // (`grosub.f` 507), so unlike a per-hour flux it cannot be rebuilt from the
    // resumed state, and an unserialized copy would make a resumed run diverge
    // from a continuous one. That is the Wave 2 `RESTART-EQUIVALENCE`
    // obligation, so it is proven here rather than assumed.
    var source = try RootState.init(std.testing.allocator, 3, 2, 4);
    defer source.deinit();
    source.retained_root_carbon_g_c_per_plant[0] = 0.125;
    source.retained_root_carbon_g_c_per_plant[2] = 17.5;
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try write(&bytes.writer, source);
    var reader: std.Io.Reader = .fixed(bytes.written());
    var restored = try read(std.testing.allocator, &reader, .{ .maximum_plants = 3, .maximum_soil_layers = 2, .maximum_root_axes = 4 });
    defer restored.deinit();
    try std.testing.expectEqualSlices(
        f64,
        source.retained_root_carbon_g_c_per_plant,
        restored.retained_root_carbon_g_c_per_plant,
    );

    // A pre-7 stream carries no `WTRTA` history. It must default to zero
    // rather than desynchronise the reflective walk: on the first resumed hour
    // the recurrence re-seeds itself from the `WTRT/PP` branch of 507--508,
    // which is the same branch a fresh planting takes. Stated as a defined
    // restart so nobody reads a pre-7 resume as restart-equivalence evidence.
    var legacy_bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer legacy_bytes.deinit();
    try writePriorVersionForTest(&legacy_bytes.writer, source, previous_version);
    var legacy_reader: std.Io.Reader = .fixed(legacy_bytes.written());
    var legacy_restored = try read(std.testing.allocator, &legacy_reader, .{ .maximum_plants = 3, .maximum_soil_layers = 2, .maximum_root_axes = 4 });
    defer legacy_restored.deinit();
    for (legacy_restored.retained_root_carbon_g_c_per_plant) |value|
        try std.testing.expectEqual(@as(f64, 0), value);
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
