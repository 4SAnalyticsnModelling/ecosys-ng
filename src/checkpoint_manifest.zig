const std = @import("std");
const execution_calendar_date = @import("execution_calendar_date.zig");

const magic = "ECOSBNDL";
const version: u32 = 2;
const atomic_replace_max_attempts: u8 = 20;
const atomic_replace_retry_delay_ms: i64 = 25;

/// Authoritative state owners that must all belong to one committed restart.
/// A manifest is accepted only when every section occurs exactly once.
pub const Section = enum(u8) {
    grid_and_plants,
    plant_metadata,
    plant_development,
    plant_roots,
    plant_canopy,
    soil_biogeochemistry,
    soil_organic_matter,
    solute_gas_and_snow_transport,
    soil_geometry_and_hydrology,
    landscape_mass_balance,
};

pub const section_count = @typeInfo(Section).@"enum".fields.len;

pub const SimulationInstant = struct {
    year: i32,
    day_of_year: u16,
    hour: u8,
    execution_iteration: u64 = 0,
    scenario_index: u64 = 0,
    scenario_iteration: u64 = 0,
    scene_index: u64 = 0,
    completed_scene_hours: u64 = 0,
};

pub const RuntimeShape = struct {
    columns: usize,
    rows: usize,
    soil_layers: usize,
    snow_layers: usize,
    plant_species_per_cell: usize,
    root_axes_per_plant: usize,
};

pub const Entry = struct {
    section: Section,
    file_name: []u8,
    byte_length: u64,
    checksum: u64,
};

pub const Manifest = struct {
    allocator: std.mem.Allocator,
    generation: u64,
    instant: SimulationInstant,
    shape: RuntimeShape,
    entries: []Entry,

    pub fn deinit(self: *Manifest) void {
        for (self.entries) |entry_value| self.allocator.free(entry_value.file_name);
        self.allocator.free(self.entries);
        self.* = undefined;
    }

    pub fn entry(self: Manifest, section: Section) ?Entry {
        for (self.entries) |candidate| if (candidate.section == section) return candidate;
        return null;
    }
};

pub const Limits = struct {
    maximum_columns: usize,
    maximum_rows: usize,
    maximum_soil_layers: usize,
    maximum_snow_layers: usize,
    maximum_plant_species_per_cell: usize,
    maximum_root_axes_per_plant: usize,
    maximum_file_name_bytes: usize = 255,
};

pub const EntrySource = struct {
    section: Section,
    file_name: []const u8,
    bytes: []const u8,
};

pub const EntryDescriptor = struct {
    section: Section,
    file_name: []const u8,
    byte_length: u64,
    checksum: u64,
};

pub fn build(allocator: std.mem.Allocator, generation: u64, instant: SimulationInstant, shape: RuntimeShape, sources: []const EntrySource) !Manifest {
    const descriptors = try allocator.alloc(EntryDescriptor, sources.len);
    defer allocator.free(descriptors);
    for (sources, descriptors) |source, *descriptor| descriptor.* = .{ .section = source.section, .file_name = source.file_name, .byte_length = source.bytes.len, .checksum = checksum(source.bytes) };
    return buildFromDescriptors(allocator, generation, instant, shape, descriptors);
}

pub fn buildFromDescriptors(allocator: std.mem.Allocator, generation: u64, instant: SimulationInstant, shape: RuntimeShape, descriptors: []const EntryDescriptor) !Manifest {
    if (descriptors.len != section_count) return error.IncompleteCheckpointBundle;
    const entries = try allocator.alloc(Entry, descriptors.len);
    errdefer allocator.free(entries);
    var initialized: usize = 0;
    errdefer for (entries[0..initialized]) |entry| allocator.free(entry.file_name);
    for (descriptors, entries) |descriptor, *destination| {
        const file_name = try allocator.dupe(u8, descriptor.file_name);
        destination.* = .{ .section = descriptor.section, .file_name = file_name, .byte_length = descriptor.byte_length, .checksum = descriptor.checksum };
        initialized += 1;
    }
    const result = Manifest{ .allocator = allocator, .generation = generation, .instant = instant, .shape = shape, .entries = entries };
    try validate(result);
    return result;
}

pub fn write(writer: anytype, manifest: Manifest) !void {
    try validate(manifest);
    try writer.writeAll(magic);
    try writer.writeInt(u32, version, .little);
    try writer.writeInt(u64, manifest.generation, .little);
    try writer.writeInt(i32, manifest.instant.year, .little);
    try writer.writeInt(u16, manifest.instant.day_of_year, .little);
    try writer.writeByte(manifest.instant.hour);
    try writer.writeInt(u64, manifest.instant.execution_iteration, .little);
    try writer.writeInt(u64, manifest.instant.scenario_index, .little);
    try writer.writeInt(u64, manifest.instant.scenario_iteration, .little);
    try writer.writeInt(u64, manifest.instant.scene_index, .little);
    try writer.writeInt(u64, manifest.instant.completed_scene_hours, .little);
    inline for (@typeInfo(RuntimeShape).@"struct".fields) |field| try writer.writeInt(u64, @intCast(@field(manifest.shape, field.name)), .little);
    try writer.writeInt(u16, @intCast(manifest.entries.len), .little);
    for (manifest.entries) |entry_value| {
        try writer.writeByte(@intFromEnum(entry_value.section));
        try writer.writeInt(u32, @intCast(entry_value.file_name.len), .little);
        try writer.writeAll(entry_value.file_name);
        try writer.writeInt(u64, entry_value.byte_length, .little);
        try writer.writeInt(u64, entry_value.checksum, .little);
    }
}

/// Publishes the commit marker last. Section files must already have been
/// written, synchronized, and renamed to their manifest names. A crash before
/// this atomic replacement leaves the previous generation authoritative.
pub fn publishAtomic(allocator: std.mem.Allocator, io: std.Io, directory: std.Io.Dir, manifest_file_name: []const u8, manifest: Manifest, buffer_bytes: usize) !void {
    if (!safeFileName(manifest_file_name)) return error.InvalidCheckpointFileName;
    if (buffer_bytes == 0) return error.InvalidCheckpointManifestBufferSize;
    try validate(manifest);
    const buffer = try allocator.alloc(u8, buffer_bytes);
    defer allocator.free(buffer);
    var atomic_file = try directory.createFileAtomic(io, manifest_file_name, .{ .replace = true });
    defer atomic_file.deinit(io);
    var file_writer = atomic_file.file.writerStreaming(io, buffer);
    try write(&file_writer.interface, manifest);
    try file_writer.interface.flush();
    try atomic_file.file.sync(io);
    var attempt: u8 = 1;
    while (true) {
        atomic_file.replace(io) catch |err| {
            if (!retryAtomicReplace(err, attempt)) return err;
            try std.Io.sleep(io, .fromMilliseconds(atomic_replace_retry_delay_ms), .awake);
            attempt += 1;
            continue;
        };
        break;
    }
}

fn retryAtomicReplace(err: anyerror, attempt: u8) bool {
    // Zig 0.16 currently reports Windows STATUS_SHARING_VIOLATION from rename
    // as Unexpected. The retry is deliberately restricted to that transient
    // class and a short bounded budget; all other storage errors fail fast.
    return err == error.Unexpected and attempt < atomic_replace_max_attempts;
}

pub fn read(allocator: std.mem.Allocator, reader: *std.Io.Reader, limits: Limits) !Manifest {
    try validateLimits(limits);
    if (!std.mem.eql(u8, try reader.takeArray(magic.len), magic)) return error.InvalidCheckpointManifestMagic;
    if (try reader.takeInt(u32, .little) != version) return error.UnsupportedCheckpointManifestVersion;
    const generation = try reader.takeInt(u64, .little);
    const instant: SimulationInstant = .{
        .year = try reader.takeInt(i32, .little),
        .day_of_year = try reader.takeInt(u16, .little),
        .hour = try reader.takeByte(),
        .execution_iteration = try reader.takeInt(u64, .little),
        .scenario_index = try reader.takeInt(u64, .little),
        .scenario_iteration = try reader.takeInt(u64, .little),
        .scene_index = try reader.takeInt(u64, .little),
        .completed_scene_hours = try reader.takeInt(u64, .little),
    };
    const shape: RuntimeShape = .{
        .columns = try bounded(reader, limits.maximum_columns, error.CheckpointColumnLimitExceeded),
        .rows = try bounded(reader, limits.maximum_rows, error.CheckpointRowLimitExceeded),
        .soil_layers = try bounded(reader, limits.maximum_soil_layers, error.CheckpointSoilLayerLimitExceeded),
        .snow_layers = try bounded(reader, limits.maximum_snow_layers, error.CheckpointSnowLayerLimitExceeded),
        .plant_species_per_cell = try bounded(reader, limits.maximum_plant_species_per_cell, error.CheckpointPlantSpeciesLimitExceeded),
        .root_axes_per_plant = try bounded(reader, limits.maximum_root_axes_per_plant, error.CheckpointRootAxisLimitExceeded),
    };
    const stored_entry_count = try reader.takeInt(u16, .little);
    if (stored_entry_count != section_count) return error.IncompleteCheckpointBundle;
    const entries = try allocator.alloc(Entry, stored_entry_count);
    errdefer allocator.free(entries);
    var initialized: usize = 0;
    errdefer for (entries[0..initialized]) |entry_value| allocator.free(entry_value.file_name);
    for (entries) |*entry_value| {
        const raw_section = try reader.takeByte();
        const section = std.enums.fromInt(Section, raw_section) orelse return error.InvalidCheckpointSection;
        const file_name_length = try reader.takeInt(u32, .little);
        if (file_name_length == 0 or file_name_length > limits.maximum_file_name_bytes) return error.InvalidCheckpointFileName;
        const file_name = try allocator.alloc(u8, file_name_length);
        errdefer allocator.free(file_name);
        try reader.readSliceAll(file_name);
        entry_value.* = .{ .section = section, .file_name = file_name, .byte_length = try reader.takeInt(u64, .little), .checksum = try reader.takeInt(u64, .little) };
        initialized += 1;
    }
    if (reader.peekByte()) |_| {
        return error.TrailingCheckpointManifestData;
    } else |err| switch (err) {
        error.EndOfStream => {},
        else => return err,
    }
    const result = Manifest{ .allocator = allocator, .generation = generation, .instant = instant, .shape = shape, .entries = entries };
    try validate(result);
    return result;
}

/// Call before parsing a section. This prevents a valid but stale or partially
/// replaced checkpoint file from being mixed into the committed generation.
pub fn verifySectionBytes(manifest: Manifest, section: Section, bytes: []const u8) !void {
    const expected = manifest.entry(section) orelse return error.MissingCheckpointSection;
    if (bytes.len != expected.byte_length) return error.CheckpointSectionLengthMismatch;
    if (checksum(bytes) != expected.checksum) return error.CheckpointSectionChecksumMismatch;
}

pub fn checksum(bytes: []const u8) u64 {
    return std.hash.Wyhash.hash(0x45434f5359534e47, bytes);
}

fn validate(manifest: Manifest) !void {
    if (manifest.generation == 0) return error.InvalidCheckpointGeneration;
    if (manifest.instant.year <= 0 or manifest.instant.year > 9999 or
        manifest.instant.hour > 23)
        return error.InvalidCheckpointInstant;
    _ = execution_calendar_date.fromDayOfYear(
        manifest.instant.day_of_year,
        @intCast(manifest.instant.year),
    ) catch return error.InvalidCheckpointInstant;
    inline for (@typeInfo(RuntimeShape).@"struct".fields) |field| if (@field(manifest.shape, field.name) == 0) return error.InvalidCheckpointRuntimeShape;
    if (manifest.entries.len != section_count) return error.IncompleteCheckpointBundle;
    var seen = [_]bool{false} ** section_count;
    for (manifest.entries) |entry_value| {
        const index = @intFromEnum(entry_value.section);
        if (seen[index]) return error.DuplicateCheckpointSection;
        seen[index] = true;
        if (!safeFileName(entry_value.file_name)) return error.InvalidCheckpointFileName;
    }
    for (seen) |present| if (!present) return error.IncompleteCheckpointBundle;
}

fn safeFileName(name: []const u8) bool {
    if (name.len == 0 or std.mem.eql(u8, name, ".") or std.mem.eql(u8, name, "..")) return false;
    if (std.mem.indexOf(u8, name, "..") != null) return false;
    if (name[0] == ' ' or name[name.len - 1] == ' ' or name[name.len - 1] == '.') return false;
    for (name) |byte| {
        if (byte == 0 or byte < 0x20 or
            std.mem.indexOfScalar(u8, "/\\<>:\"|?*", byte) != null)
            return false;
    }
    return true;
}

fn validateLimits(limits: Limits) !void {
    inline for (@typeInfo(Limits).@"struct".fields) |field| if (@field(limits, field.name) == 0) return error.InvalidCheckpointManifestLimits;
}

fn bounded(reader: *std.Io.Reader, limit: usize, comptime too_large: anyerror) !usize {
    const value = try reader.takeInt(u64, .little);
    if (value > limit or value > std.math.maxInt(usize)) return too_large;
    return @intCast(value);
}

fn testSources() [section_count]EntrySource {
    var result: [section_count]EntrySource = undefined;
    inline for (@typeInfo(Section).@"enum".fields, 0..) |field, index| result[index] = .{ .section = @enumFromInt(field.value), .file_name = field.name ++ ".bin", .bytes = field.name };
    return result;
}

const test_shape: RuntimeShape = .{ .columns = 3, .rows = 2, .soil_layers = 7, .snow_layers = 4, .plant_species_per_cell = 12, .root_axes_per_plant = 16 };
const test_instant: SimulationInstant = .{ .year = 2001, .day_of_year = 123, .hour = 17, .execution_iteration = 2, .scenario_index = 3, .scenario_iteration = 4, .scene_index = 5, .completed_scene_hours = 876 };
const test_limits: Limits = .{ .maximum_columns = 10, .maximum_rows = 10, .maximum_soil_layers = 20, .maximum_snow_layers = 10, .maximum_plant_species_per_cell = 30, .maximum_root_axes_per_plant = 40 };

test "checkpoint manifest round trips complete arbitrary runtime generation" {
    const sources = testSources();
    var manifest = try build(std.testing.allocator, 42, test_instant, test_shape, &sources);
    defer manifest.deinit();
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try write(&bytes.writer, manifest);
    var reader: std.Io.Reader = .fixed(bytes.written());
    var restored = try read(std.testing.allocator, &reader, test_limits);
    defer restored.deinit();
    try std.testing.expectEqual(@as(u64, 42), restored.generation);
    try std.testing.expectEqual(test_instant, restored.instant);
    try std.testing.expectEqual(@as(usize, 12), restored.shape.plant_species_per_cell);
    try verifySectionBytes(restored, .plant_roots, "plant_roots");
}

test "checkpoint instant preserves DAY modulo-four chronology" {
    const sources = testSources();
    var leap_instant = test_instant;
    leap_instant.year = 1900;
    leap_instant.day_of_year = 366;
    var manifest = try build(
        std.testing.allocator,
        1,
        leap_instant,
        test_shape,
        &sources,
    );
    defer manifest.deinit();

    leap_instant.year = 1901;
    try std.testing.expectError(
        error.InvalidCheckpointInstant,
        build(std.testing.allocator, 1, leap_instant, test_shape, &sources),
    );
    leap_instant.year = -1;
    leap_instant.day_of_year = 1;
    try std.testing.expectError(
        error.InvalidCheckpointInstant,
        build(std.testing.allocator, 1, leap_instant, test_shape, &sources),
    );
}

test "checkpoint manifest rejects duplicate missing unsafe and corrupt sections" {
    var sources = testSources();
    sources[1].section = sources[0].section;
    try std.testing.expectError(error.DuplicateCheckpointSection, build(std.testing.allocator, 1, test_instant, test_shape, &sources));
    sources = testSources();
    sources[0].file_name = "../grid.bin";
    try std.testing.expectError(error.InvalidCheckpointFileName, build(std.testing.allocator, 1, test_instant, test_shape, &sources));
    try std.testing.expectError(error.IncompleteCheckpointBundle, build(std.testing.allocator, 1, test_instant, test_shape, sources[0 .. sources.len - 1]));
    sources = testSources();
    var manifest = try build(std.testing.allocator, 1, test_instant, test_shape, &sources);
    defer manifest.deinit();
    try std.testing.expectError(error.CheckpointSectionChecksumMismatch, verifySectionBytes(manifest, .plant_roots, "plant_rootx"));
}

test "checkpoint manifest rejects nonportable file names before I/O" {
    inline for (.{
        "",
        ".",
        "..",
        "../grid.bin",
        "grid/section.bin",
        "grid\\section.bin",
        " grid.bin",
        "grid.bin ",
        "grid.",
        "grid:section.bin",
        "grid|section.bin",
        "grid?section.bin",
        "grid*section.bin",
        "grid<section.bin",
        "grid>section.bin",
        "grid\"section.bin",
    }) |name| {
        var sources = testSources();
        sources[0].file_name = name;
        try std.testing.expectError(
            error.InvalidCheckpointFileName,
            build(std.testing.allocator, 1, test_instant, test_shape, &sources),
        );
    }
}

test "checkpoint manifest accepts simple extensionless and binary file names" {
    inline for (.{ "grid_state", "grid-state.bin", "grid.state.bin" }) |name| {
        var sources = testSources();
        sources[0].file_name = name;
        var manifest = try build(
            std.testing.allocator,
            1,
            test_instant,
            test_shape,
            &sources,
        );
        defer manifest.deinit();
        try std.testing.expectEqualStrings(name, manifest.entries[0].file_name);
    }
}

test "checkpoint manifest enforces runtime limits before entry allocation" {
    const sources = testSources();
    var manifest = try build(std.testing.allocator, 1, test_instant, test_shape, &sources);
    defer manifest.deinit();
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try write(&bytes.writer, manifest);
    var reader: std.Io.Reader = .fixed(bytes.written());
    var limits = test_limits;
    limits.maximum_plant_species_per_cell = 5;
    try std.testing.expectError(error.CheckpointPlantSpeciesLimitExceeded, read(std.testing.allocator, &reader, limits));
}

test "checkpoint manifest is published as one atomic commit marker" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const sources = testSources();
    var first = try build(std.testing.allocator, 1, test_instant, test_shape, &sources);
    defer first.deinit();
    try publishAtomic(std.testing.allocator, std.testing.io, temporary.dir, "restart.manifest", first, 128);
    var second = try build(std.testing.allocator, 2, test_instant, test_shape, &sources);
    defer second.deinit();
    try publishAtomic(std.testing.allocator, std.testing.io, temporary.dir, "restart.manifest", second, 128);
    const bytes = try temporary.dir.readFileAlloc(std.testing.io, "restart.manifest", std.testing.allocator, .limited(16 * 1024));
    defer std.testing.allocator.free(bytes);
    var reader: std.Io.Reader = .fixed(bytes);
    var restored = try read(std.testing.allocator, &reader, test_limits);
    defer restored.deinit();
    try std.testing.expectEqual(@as(u64, 2), restored.generation);
}

test "atomic checkpoint replacement retries only bounded transient sharing failures" {
    try std.testing.expect(retryAtomicReplace(error.Unexpected, 1));
    try std.testing.expect(retryAtomicReplace(error.Unexpected, atomic_replace_max_attempts - 1));
    try std.testing.expect(!retryAtomicReplace(error.Unexpected, atomic_replace_max_attempts));
    try std.testing.expect(!retryAtomicReplace(error.AccessDenied, 1));
}
