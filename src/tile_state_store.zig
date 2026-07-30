const std = @import("std");
const spatial_grid = @import("spatial_grid.zig");

const magic = "ECOSTILE";
const format_version: u32 = 6;
const morton_order_tag: u8 = 1;
const generation_manifest_name = "generation.morton.manifest";
const generation_manifest_magic = "ECOGEN01";
const generation_manifest_version: u32 = 3;

const ScalarType = enum(u8) {
    float64 = 1,
    unsigned64 = 2,
};

pub const Field = struct {
    name: []const u8 = "unnamed",
    values: []const f64 = &.{},
    unsigned_values: ?[]const usize = null,
    components_per_cell: usize,
};

pub const MutableField = struct {
    name: []const u8 = "unnamed",
    values: []f64 = &.{},
    unsigned_values: ?[]usize = null,
    components_per_cell: usize,
};

/// One binary file per tile. The file name is keyed by the tile's Morton
/// index, and values inside the file are also laid out in Morton order. This
/// keeps an out-of-core load/compute/commit transaction contiguous on disk.
pub const FileStore = struct {
    allocator: std.mem.Allocator,
    io: std.Io,
    directory: std.Io.Dir,
    buffer_byte_count: usize,
    generation_index: u64,

    pub fn init(
        allocator: std.mem.Allocator,
        io: std.Io,
        directory: std.Io.Dir,
        buffer_byte_count: usize,
    ) !FileStore {
        if (buffer_byte_count == 0) return error.InvalidTileIoBufferSize;
        return .{
            .allocator = allocator,
            .io = io,
            .directory = directory,
            .buffer_byte_count = buffer_byte_count,
            .generation_index = 0,
        };
    }

    /// Creates a store whose identity participates in source/destination
    /// generation validation. Generation indices must increase monotonically
    /// in the simulation driver; the binary records remain keyed by Morton
    /// tile index within the generation directory.
    pub fn initGeneration(
        allocator: std.mem.Allocator,
        io: std.Io,
        directory: std.Io.Dir,
        buffer_byte_count: usize,
        generation_index: u64,
    ) !FileStore {
        var store = try init(allocator, io, directory, buffer_byte_count);
        store.generation_index = generation_index;
        return store;
    }

    pub fn saveOwnedFields(
        self: FileStore,
        tile: spatial_grid.Tile,
        grid_row_count: usize,
        grid_column_count: usize,
        fields: []const Field,
    ) !void {
        const file_name = try tileFileName(self.allocator, tile.z_order_index);
        defer self.allocator.free(file_name);
        const buffer = try self.allocator.alloc(u8, self.buffer_byte_count);
        defer self.allocator.free(buffer);
        var atomic_file = try self.directory.createFileAtomic(
            self.io,
            file_name,
            .{ .replace = true },
        );
        defer atomic_file.deinit(self.io);
        var writer = atomic_file.file.writerStreaming(self.io, buffer);
        try writeOwnedFieldsForGeneration(
            self.allocator,
            &writer.interface,
            self.generation_index,
            tile,
            grid_row_count,
            grid_column_count,
            fields,
        );
        try writer.interface.flush();
        try atomic_file.file.sync(self.io);
        try atomic_file.replace(self.io);
    }

    pub fn loadTileFieldsInto(
        self: FileStore,
        plan: spatial_grid.TilePlan,
        destination_tile: spatial_grid.Tile,
        fields: []const MutableField,
    ) !void {
        // Tile files contain authoritative owned interiors only. Assemble the
        // requested tile and its halo from every intersecting owner, avoiding
        // duplicated mutable halo records that could become stale.
        const buffer = try self.allocator.alloc(u8, self.buffer_byte_count);
        defer self.allocator.free(buffer);
        for (plan.tiles) |source_tile| {
            if (!ownedIntersectsLoaded(source_tile, destination_tile)) continue;
            const file_name = try tileFileName(
                self.allocator,
                source_tile.z_order_index,
            );
            defer self.allocator.free(file_name);
            var file = try self.directory.openFile(self.io, file_name, .{});
            defer file.close(self.io);
            var reader = file.readerStreaming(self.io, buffer);
            try readOwnedFieldsIntoForGeneration(
                self.allocator,
                &reader.interface,
                self.generation_index,
                source_tile,
                destination_tile,
                plan.grid_row_count,
                plan.grid_column_count,
                fields,
            );
        }
    }

    /// Publishes the generation only after every owned tile file is durable.
    /// The manifest is the sole completion marker and is atomically replaced.
    pub fn publishGeneration(
        self: FileStore,
        plan: spatial_grid.TilePlan,
        fields: []const Field,
    ) !void {
        try validateFieldSchema(fields);
        const buffer = try self.allocator.alloc(u8, self.buffer_byte_count);
        defer self.allocator.free(buffer);
        var atomic_file = try self.directory.createFileAtomic(
            self.io,
            generation_manifest_name,
            .{ .replace = true },
        );
        defer atomic_file.deinit(self.io);
        var writer = atomic_file.file.writerStreaming(self.io, buffer);
        try writeGenerationManifest(
            &writer.interface,
            self.generation_index,
            plan,
            fields,
        );
        try writer.interface.flush();
        try atomic_file.file.sync(self.io);
        try atomic_file.replace(self.io);
    }

    /// Requires a complete, shape-compatible manifest and every referenced
    /// Morton tile file before any science kernel is allowed to run.
    pub fn validateGeneration(
        self: FileStore,
        plan: spatial_grid.TilePlan,
        fields: []const MutableField,
    ) !void {
        try validateMutableFieldSchema(fields);
        const buffer = try self.allocator.alloc(u8, self.buffer_byte_count);
        defer self.allocator.free(buffer);
        var manifest_file = try self.directory.openFile(
            self.io,
            generation_manifest_name,
            .{},
        );
        defer manifest_file.close(self.io);
        var reader = manifest_file.readerStreaming(self.io, buffer);
        try readAndValidateGenerationManifest(
            &reader.interface,
            self.generation_index,
            plan,
            fields,
        );
        for (plan.tiles) |tile| {
            const file_name = try tileFileName(
                self.allocator,
                tile.z_order_index,
            );
            defer self.allocator.free(file_name);
            var tile_file = try self.directory.openFile(
                self.io,
                file_name,
                .{},
            );
            defer tile_file.close(self.io);
            var tile_reader = tile_file.readerStreaming(self.io, buffer);
            try validateOwnedFieldsForGeneration(
                self.allocator,
                &tile_reader.interface,
                self.generation_index,
                tile,
                plan.grid_row_count,
                plan.grid_column_count,
                fields,
            );
        }
    }

    pub fn requireUnpublished(self: FileStore) !void {
        var file = self.directory.openFile(
            self.io,
            generation_manifest_name,
            .{},
        ) catch |err| switch (err) {
            error.FileNotFound => return,
            else => return err,
        };
        file.close(self.io);
        return error.TileDestinationGenerationAlreadyPublished;
    }
};

fn validateFieldSchema(fields: []const Field) !void {
    if (fields.len == 0 or fields.len > std.math.maxInt(u32))
        return error.InvalidTileFieldCount;
    for (fields) |field| {
        if (field.name.len == 0 or field.name.len > std.math.maxInt(u32) or
            field.components_per_cell == 0)
            return error.TileFieldDimensionMismatch;
        _ = try fieldScalarType(field);
    }
    try ensureUniqueFieldNames(fields);
}

fn validateMutableFieldSchema(fields: []const MutableField) !void {
    if (fields.len == 0 or fields.len > std.math.maxInt(u32))
        return error.InvalidTileFieldCount;
    for (fields) |field| {
        if (field.name.len == 0 or field.name.len > std.math.maxInt(u32) or
            field.components_per_cell == 0)
            return error.TileFieldDimensionMismatch;
        _ = try mutableFieldScalarType(field);
    }
    try ensureUniqueFieldNames(fields);
}

fn fieldScalarType(field: Field) !ScalarType {
    if (field.unsigned_values != null) {
        if (field.values.len != 0)
            return error.AmbiguousTileFieldScalarType;
        return .unsigned64;
    }
    return .float64;
}

fn mutableFieldScalarType(field: MutableField) !ScalarType {
    if (field.unsigned_values != null) {
        if (field.values.len != 0)
            return error.AmbiguousTileFieldScalarType;
        return .unsigned64;
    }
    return .float64;
}

fn fieldValueCount(field: Field) usize {
    return if (field.unsigned_values) |values|
        values.len
    else
        field.values.len;
}

fn mutableFieldValueCount(field: MutableField) usize {
    return if (field.unsigned_values) |values|
        values.len
    else
        field.values.len;
}

fn ensureUniqueFieldNames(fields: anytype) !void {
    for (fields, 0..) |field, index| {
        for (fields[0..index]) |previous| {
            if (std.mem.eql(u8, field.name, previous.name))
                return error.DuplicateTileFieldName;
        }
    }
}

fn writeFieldName(writer: *std.Io.Writer, name: []const u8) !void {
    if (name.len == 0 or name.len > std.math.maxInt(u32))
        return error.InvalidTileFieldName;
    try writer.writeInt(u32, @intCast(name.len), .little);
    try writer.writeAll(name);
}

fn validateFieldName(
    reader: *std.Io.Reader,
    expected: []const u8,
) !void {
    const byte_count = try reader.takeInt(u32, .little);
    if (byte_count != expected.len)
        return error.TileStateFieldNameMismatch;
    if (!std.mem.eql(u8, try reader.take(byte_count), expected))
        return error.TileStateFieldNameMismatch;
}

fn writeGenerationManifest(
    writer: *std.Io.Writer,
    generation_index: u64,
    plan: spatial_grid.TilePlan,
    fields: []const Field,
) !void {
    try writer.writeAll(generation_manifest_magic);
    try writer.writeInt(u32, generation_manifest_version, .little);
    try writer.writeInt(u64, generation_index, .little);
    try writer.writeByte(morton_order_tag);
    try writer.writeInt(u64, @intCast(plan.grid_row_count), .little);
    try writer.writeInt(u64, @intCast(plan.grid_column_count), .little);
    try writer.writeInt(u32, @intCast(fields.len), .little);
    for (fields) |field| {
        try writeFieldName(writer, field.name);
        try writer.writeByte(@intFromEnum(try fieldScalarType(field)));
        try writer.writeInt(
            u64,
            @intCast(field.components_per_cell),
            .little,
        );
    }
    try writer.writeInt(u64, @intCast(plan.tiles.len), .little);
    for (plan.tiles) |tile| {
        try writer.writeInt(u64, tile.z_order_index, .little);
        inline for (.{
            tile.owned_north_row,
            tile.owned_west_column,
            tile.owned_south_row_exclusive,
            tile.owned_east_column_exclusive,
            tile.loaded_north_row,
            tile.loaded_west_column,
            tile.loaded_south_row_exclusive,
            tile.loaded_east_column_exclusive,
        }) |value| try writer.writeInt(u64, @intCast(value), .little);
    }
}

fn readAndValidateGenerationManifest(
    reader: *std.Io.Reader,
    generation_index: u64,
    plan: spatial_grid.TilePlan,
    fields: []const MutableField,
) !void {
    const file_magic = try reader.takeArray(generation_manifest_magic.len);
    if (!std.mem.eql(u8, file_magic, generation_manifest_magic))
        return error.InvalidTileGenerationManifestMagic;
    if (try reader.takeInt(u32, .little) != generation_manifest_version)
        return error.UnsupportedTileGenerationManifestVersion;
    if (try reader.takeInt(u64, .little) != generation_index)
        return error.TileGenerationManifestIndexMismatch;
    if (try reader.takeByte() != morton_order_tag)
        return error.UnsupportedTileTraversalOrder;
    if (try reader.takeInt(u64, .little) != plan.grid_row_count or
        try reader.takeInt(u64, .little) != plan.grid_column_count)
        return error.TileGenerationManifestGridMismatch;
    if (try reader.takeInt(u32, .little) != fields.len)
        return error.TileStateFieldCountMismatch;
    for (fields) |field| {
        try validateFieldName(reader, field.name);
        if (try reader.takeByte() !=
            @intFromEnum(try mutableFieldScalarType(field)))
            return error.TileStateFieldScalarTypeMismatch;
        if (try reader.takeInt(u64, .little) !=
            field.components_per_cell)
            return error.TileStateFieldShapeMismatch;
    }
    if (try reader.takeInt(u64, .little) != plan.tiles.len)
        return error.TileGenerationManifestTileCountMismatch;
    for (plan.tiles) |tile| {
        if (try reader.takeInt(u64, .little) != tile.z_order_index)
            return error.TileGenerationManifestOrderMismatch;
        inline for (.{
            tile.owned_north_row,
            tile.owned_west_column,
            tile.owned_south_row_exclusive,
            tile.owned_east_column_exclusive,
            tile.loaded_north_row,
            tile.loaded_west_column,
            tile.loaded_south_row_exclusive,
            tile.loaded_east_column_exclusive,
        }) |expected| if (try reader.takeInt(u64, .little) != expected)
            return error.TileGenerationManifestBoundsMismatch;
    }
    if (reader.peekByte()) |_| return error.TrailingTileGenerationManifestData else |err| switch (err) {
        error.EndOfStream => {},
        else => return err,
    }
}

fn tileFileName(allocator: std.mem.Allocator, z_order_index: u64) ![]u8 {
    return std.fmt.allocPrint(
        allocator,
        "tile_{x:0>16}.morton.bin",
        .{z_order_index},
    );
}

/// Writes loaded tile fields in global Morton/Z-order. Field order is a
/// caller-owned schema; every field contains one scalar per global grid cell.
pub fn writeOwnedFields(
    allocator: std.mem.Allocator,
    writer: anytype,
    tile: spatial_grid.Tile,
    grid_row_count: usize,
    grid_column_count: usize,
    fields: []const Field,
) !void {
    return writeOwnedFieldsForGeneration(
        allocator,
        writer,
        0,
        tile,
        grid_row_count,
        grid_column_count,
        fields,
    );
}

fn writeOwnedFieldsForGeneration(
    allocator: std.mem.Allocator,
    writer: anytype,
    generation_index: u64,
    tile: spatial_grid.Tile,
    grid_row_count: usize,
    grid_column_count: usize,
    fields: []const Field,
) !void {
    const grid_cell_count = try std.math.mul(usize, grid_row_count, grid_column_count);
    if (fields.len == 0 or fields.len > std.math.maxInt(u32))
        return error.InvalidTileFieldCount;
    try ensureUniqueFieldNames(fields);
    for (fields) |field| {
        _ = try fieldScalarType(field);
        if (field.name.len == 0 or field.name.len > std.math.maxInt(u32) or
            field.components_per_cell == 0 or
            fieldValueCount(field) !=
                try std.math.mul(
                    usize,
                    grid_cell_count,
                    field.components_per_cell,
                ))
            return error.TileFieldDimensionMismatch;
    }
    const indices = try tile.ownedCellIndicesZOrder(
        allocator,
        grid_row_count,
        grid_column_count,
    );
    defer allocator.free(indices);

    try writer.writeAll(magic);
    try writer.writeInt(u32, format_version, .little);
    try writer.writeInt(u64, generation_index, .little);
    // Both the serial tile schedule and the records within each tile use
    // Morton order. These tags make that persistent contract explicit.
    try writer.writeByte(morton_order_tag);
    try writer.writeByte(morton_order_tag);
    try writer.writeInt(u64, @intCast(grid_row_count), .little);
    try writer.writeInt(u64, @intCast(grid_column_count), .little);
    try writer.writeInt(u64, tile.z_order_index, .little);
    inline for (.{
        tile.owned_north_row,
        tile.owned_west_column,
        tile.owned_south_row_exclusive,
        tile.owned_east_column_exclusive,
        tile.loaded_north_row,
        tile.loaded_west_column,
        tile.loaded_south_row_exclusive,
        tile.loaded_east_column_exclusive,
    }) |value| try writer.writeInt(u64, @intCast(value), .little);
    try writer.writeInt(u32, @intCast(fields.len), .little);
    for (fields) |field| {
        try writeFieldName(writer, field.name);
        try writer.writeByte(@intFromEnum(try fieldScalarType(field)));
        try writer.writeInt(
            u64,
            @intCast(field.components_per_cell),
            .little,
        );
    }
    try writer.writeInt(u64, @intCast(indices.len), .little);
    for (fields) |field| for (indices) |cell| {
        const first = cell * field.components_per_cell;
        if (field.unsigned_values) |values| {
            for (values[first..][0..field.components_per_cell]) |value|
                try writer.writeInt(u64, @intCast(value), .little);
        } else {
            for (field.values[first..][0..field.components_per_cell]) |value| {
                if (!std.math.isFinite(value))
                    return error.NonFiniteTileState;
                try writer.writeInt(u64, @bitCast(value), .little);
            }
        }
    };
}

/// Loads the halo and owned interior into preallocated global fields. The
/// caller later commits only the owned interior after its parallel tile kernel.
pub fn readOwnedFieldsInto(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    source_tile: spatial_grid.Tile,
    destination_tile: spatial_grid.Tile,
    grid_row_count: usize,
    grid_column_count: usize,
    fields: []const MutableField,
) !void {
    return readOwnedFieldsIntoForGeneration(
        allocator,
        reader,
        0,
        source_tile,
        destination_tile,
        grid_row_count,
        grid_column_count,
        fields,
    );
}

fn readOwnedFieldsIntoForGeneration(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    generation_index: u64,
    source_tile: spatial_grid.Tile,
    destination_tile: spatial_grid.Tile,
    grid_row_count: usize,
    grid_column_count: usize,
    fields: []const MutableField,
) !void {
    return consumeOwnedFieldsForGeneration(
        allocator,
        reader,
        generation_index,
        source_tile,
        destination_tile,
        grid_row_count,
        grid_column_count,
        fields,
        true,
    );
}

fn validateOwnedFieldsForGeneration(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    generation_index: u64,
    source_tile: spatial_grid.Tile,
    grid_row_count: usize,
    grid_column_count: usize,
    fields: []const MutableField,
) !void {
    return consumeOwnedFieldsForGeneration(
        allocator,
        reader,
        generation_index,
        source_tile,
        source_tile,
        grid_row_count,
        grid_column_count,
        fields,
        false,
    );
}

fn consumeOwnedFieldsForGeneration(
    allocator: std.mem.Allocator,
    reader: *std.Io.Reader,
    generation_index: u64,
    source_tile: spatial_grid.Tile,
    destination_tile: spatial_grid.Tile,
    grid_row_count: usize,
    grid_column_count: usize,
    fields: []const MutableField,
    apply_values: bool,
) !void {
    const grid_cell_count = try std.math.mul(usize, grid_row_count, grid_column_count);
    if (fields.len == 0 or fields.len > std.math.maxInt(u32))
        return error.InvalidTileFieldCount;
    try ensureUniqueFieldNames(fields);
    for (fields) |field| {
        _ = try mutableFieldScalarType(field);
        if (field.name.len == 0 or field.name.len > std.math.maxInt(u32) or
            field.components_per_cell == 0 or
            mutableFieldValueCount(field) !=
                try std.math.mul(
                    usize,
                    grid_cell_count,
                    field.components_per_cell,
                ))
            return error.TileFieldDimensionMismatch;
    }
    const file_magic = try reader.takeArray(magic.len);
    if (!std.mem.eql(u8, file_magic, magic)) return error.InvalidTileStateMagic;
    if (try reader.takeInt(u32, .little) != format_version)
        return error.UnsupportedTileStateVersion;
    if (try reader.takeInt(u64, .little) != generation_index)
        return error.TileStateGenerationMismatch;
    if (try reader.takeByte() != morton_order_tag)
        return error.UnsupportedTileTraversalOrder;
    if (try reader.takeByte() != morton_order_tag)
        return error.UnsupportedTileCellOrder;
    if (try reader.takeInt(u64, .little) != grid_row_count or
        try reader.takeInt(u64, .little) != grid_column_count or
        try reader.takeInt(u64, .little) != source_tile.z_order_index)
        return error.TileStateGridMismatch;
    inline for (.{
        source_tile.owned_north_row,
        source_tile.owned_west_column,
        source_tile.owned_south_row_exclusive,
        source_tile.owned_east_column_exclusive,
        source_tile.loaded_north_row,
        source_tile.loaded_west_column,
        source_tile.loaded_south_row_exclusive,
        source_tile.loaded_east_column_exclusive,
    }) |expected| if (try reader.takeInt(u64, .little) != expected)
        return error.TileStateBoundsMismatch;
    if (try reader.takeInt(u32, .little) != fields.len)
        return error.TileStateFieldCountMismatch;
    for (fields) |field| {
        try validateFieldName(reader, field.name);
        if (try reader.takeByte() !=
            @intFromEnum(try mutableFieldScalarType(field)))
            return error.TileStateFieldScalarTypeMismatch;
        if (try reader.takeInt(u64, .little) !=
            field.components_per_cell)
            return error.TileStateFieldShapeMismatch;
    }
    const indices = try source_tile.ownedCellIndicesZOrder(
        allocator,
        grid_row_count,
        grid_column_count,
    );
    defer allocator.free(indices);
    if (try reader.takeInt(u64, .little) != indices.len)
        return error.TileStateCellCountMismatch;
    for (fields) |field| for (indices) |cell| {
        const first = cell * field.components_per_cell;
        const row = cell / grid_column_count;
        const column = cell % grid_column_count;
        const inside_loaded =
            row >= destination_tile.loaded_north_row and
            row < destination_tile.loaded_south_row_exclusive and
            column >= destination_tile.loaded_west_column and
            column < destination_tile.loaded_east_column_exclusive;
        if (field.unsigned_values) |values| {
            for (values[first..][0..field.components_per_cell]) |*value| {
                const decoded = try reader.takeInt(u64, .little);
                const runtime_value = std.math.cast(usize, decoded) orelse
                    return error.TileUnsignedValueExceedsRuntime;
                if (apply_values and inside_loaded)
                    value.* = runtime_value;
            }
        } else {
            for (field.values[first..][0..field.components_per_cell]) |*value| {
                const decoded: f64 = @bitCast(
                    try reader.takeInt(u64, .little),
                );
                if (!std.math.isFinite(decoded))
                    return error.NonFiniteTileState;
                if (apply_values and inside_loaded) value.* = decoded;
            }
        }
    };
    if (reader.peekByte()) |_| return error.TrailingTileStateData else |err| switch (err) {
        error.EndOfStream => {},
        else => return err,
    }
}

fn ownedIntersectsLoaded(
    source: spatial_grid.Tile,
    destination: spatial_grid.Tile,
) bool {
    return source.owned_north_row < destination.loaded_south_row_exclusive and
        source.owned_south_row_exclusive > destination.loaded_north_row and
        source.owned_west_column < destination.loaded_east_column_exclusive and
        source.owned_east_column_exclusive > destination.loaded_west_column;
}

test "binary tile state round trips owned cells in Morton order" {
    var plan = try spatial_grid.TilePlan.init(std.testing.allocator, 6, 7, 2, 3, 2);
    defer plan.deinit();
    const tile = plan.tiles[3];
    const cells = 6 * 7;
    const source_a = try std.testing.allocator.alloc(f64, cells);
    defer std.testing.allocator.free(source_a);
    const source_b = try std.testing.allocator.alloc(f64, cells);
    defer std.testing.allocator.free(source_b);
    for (source_a, source_b, 0..) |*a, *b, cell| {
        a.* = @floatFromInt(cell + 1);
        b.* = -a.*;
    }
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try writeOwnedFields(
        std.testing.allocator,
        &bytes.writer,
        tile,
        6,
        7,
        &.{
            .{ .name = "source_a", .values = source_a, .components_per_cell = 1 },
            .{ .name = "source_b", .values = source_b, .components_per_cell = 1 },
        },
    );
    const destination_a = try std.testing.allocator.alloc(f64, cells);
    defer std.testing.allocator.free(destination_a);
    const destination_b = try std.testing.allocator.alloc(f64, cells);
    defer std.testing.allocator.free(destination_b);
    @memset(destination_a, 0);
    @memset(destination_b, 0);
    var reader: std.Io.Reader = .fixed(bytes.written());
    try readOwnedFieldsInto(
        std.testing.allocator,
        &reader,
        tile,
        tile,
        6,
        7,
        &.{
            .{ .name = "source_a", .values = destination_a, .components_per_cell = 1 },
            .{ .name = "source_b", .values = destination_b, .components_per_cell = 1 },
        },
    );
    const indices = try tile.ownedCellIndicesZOrder(std.testing.allocator, 6, 7);
    defer std.testing.allocator.free(indices);
    for (indices) |cell| {
        try std.testing.expectEqual(source_a[cell], destination_a[cell]);
        try std.testing.expectEqual(source_b[cell], destination_b[cell]);
    }
}

test "binary tile state rejects non-Morton traversal metadata" {
    var plan = try spatial_grid.TilePlan.init(std.testing.allocator, 1, 1, 1, 1, 2);
    defer plan.deinit();
    const source = [_]f64{1};
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try writeOwnedFields(
        std.testing.allocator,
        &bytes.writer,
        plan.tiles[0],
        1,
        1,
        &.{.{ .values = &source, .components_per_cell = 1 }},
    );
    const encoded = bytes.written();
    encoded[magic.len + @sizeOf(u32) + @sizeOf(u64)] = 0;
    var destination = [_]f64{0};
    var reader: std.Io.Reader = .fixed(encoded);
    try std.testing.expectError(
        error.UnsupportedTileTraversalOrder,
        readOwnedFieldsInto(
            std.testing.allocator,
            &reader,
            plan.tiles[0],
            plan.tiles[0],
            1,
            1,
            &.{.{ .values = &destination, .components_per_cell = 1 }},
        ),
    );
}

test "binary tile state rejects reordered same-shaped fields" {
    var plan = try spatial_grid.TilePlan.init(
        std.testing.allocator,
        1,
        1,
        1,
        1,
        2,
    );
    defer plan.deinit();
    const first = [_]f64{1};
    const second = [_]f64{2};
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try writeOwnedFields(
        std.testing.allocator,
        &bytes.writer,
        plan.tiles[0],
        1,
        1,
        &.{
            .{ .name = "water_m3", .values = &first, .components_per_cell = 1 },
            .{ .name = "temperature_k", .values = &second, .components_per_cell = 1 },
        },
    );
    var destination_first = [_]f64{0};
    var destination_second = [_]f64{0};
    var reader: std.Io.Reader = .fixed(bytes.written());
    try std.testing.expectError(
        error.TileStateFieldNameMismatch,
        readOwnedFieldsInto(
            std.testing.allocator,
            &reader,
            plan.tiles[0],
            plan.tiles[0],
            1,
            1,
            &.{
                .{ .name = "temperature_k", .values = &destination_second, .components_per_cell = 1 },
                .{ .name = "water_m3", .values = &destination_first, .components_per_cell = 1 },
            },
        ),
    );
}

test "binary tile state preserves unsigned topology and rejects type drift" {
    var plan = try spatial_grid.TilePlan.init(
        std.testing.allocator,
        1,
        3,
        1,
        3,
        2,
    );
    defer plan.deinit();
    const source = [_]usize{ 1, 4, 9 };
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try writeOwnedFields(
        std.testing.allocator,
        &bytes.writer,
        plan.tiles[0],
        1,
        3,
        &.{.{
            .name = "active_soil_layer_count",
            .unsigned_values = &source,
            .components_per_cell = 1,
        }},
    );
    var destination = [_]usize{ 0, 0, 0 };
    var reader: std.Io.Reader = .fixed(bytes.written());
    try readOwnedFieldsInto(
        std.testing.allocator,
        &reader,
        plan.tiles[0],
        plan.tiles[0],
        1,
        3,
        &.{.{
            .name = "active_soil_layer_count",
            .unsigned_values = &destination,
            .components_per_cell = 1,
        }},
    );
    try std.testing.expectEqualSlices(usize, &source, &destination);

    var incorrect_float_destination = [_]f64{ 0, 0, 0 };
    var type_drift_reader: std.Io.Reader = .fixed(bytes.written());
    try std.testing.expectError(
        error.TileStateFieldScalarTypeMismatch,
        readOwnedFieldsInto(
            std.testing.allocator,
            &type_drift_reader,
            plan.tiles[0],
            plan.tiles[0],
            1,
            3,
            &.{.{
                .name = "active_soil_layer_count",
                .values = &incorrect_float_destination,
                .components_per_cell = 1,
            }},
        ),
    );
}

test "file store atomically persists one Morton binary file per tile" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    var plan = try spatial_grid.TilePlan.init(std.testing.allocator, 4, 5, 2, 2, 2);
    defer plan.deinit();
    const tile = plan.tiles[2];
    const source = try std.testing.allocator.alloc(f64, 4 * 5);
    defer std.testing.allocator.free(source);
    for (source, 0..) |*value, cell| value.* = @floatFromInt(cell + 11);
    const destination = try std.testing.allocator.alloc(f64, 4 * 5);
    defer std.testing.allocator.free(destination);
    @memset(destination, 0);

    const store = try FileStore.init(
        std.testing.allocator,
        std.testing.io,
        temporary.dir,
        256,
    );
    for (plan.tiles) |source_tile| try store.saveOwnedFields(
        source_tile,
        4,
        5,
        &.{.{ .values = source, .components_per_cell = 1 }},
    );
    try store.loadTileFieldsInto(
        plan,
        tile,
        &.{.{ .values = destination, .components_per_cell = 1 }},
    );

    const loaded_cells = try tile.loadedCellIndicesZOrder(
        std.testing.allocator,
        4,
        5,
    );
    defer std.testing.allocator.free(loaded_cells);
    for (loaded_cells) |cell|
        try std.testing.expectEqual(source[cell], destination[cell]);

    const file_name = try tileFileName(std.testing.allocator, tile.z_order_index);
    defer std.testing.allocator.free(file_name);
    var file = try temporary.dir.openFile(std.testing.io, file_name, .{});
    file.close(std.testing.io);
}

test "generation validation rejects a published manifest with a missing tile" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    var plan = try spatial_grid.TilePlan.init(
        std.testing.allocator,
        2,
        4,
        2,
        2,
        2,
    );
    defer plan.deinit();
    var values = [_]f64{1} ** 8;
    const store = try FileStore.initGeneration(
        std.testing.allocator,
        std.testing.io,
        temporary.dir,
        256,
        9,
    );
    try store.saveOwnedFields(
        plan.tiles[0],
        2,
        4,
        &.{.{ .values = &values, .components_per_cell = 1 }},
    );
    try store.publishGeneration(
        plan,
        &.{.{ .values = &values, .components_per_cell = 1 }},
    );
    try std.testing.expectError(
        error.FileNotFound,
        store.validateGeneration(
            plan,
            &.{.{ .values = &values, .components_per_cell = 1 }},
        ),
    );
}

test "generation validation rejects runtime field shape drift" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    var plan = try spatial_grid.TilePlan.init(
        std.testing.allocator,
        1,
        1,
        1,
        1,
        2,
    );
    defer plan.deinit();
    var values = [_]f64{ 1, 2 };
    const store = try FileStore.initGeneration(
        std.testing.allocator,
        std.testing.io,
        temporary.dir,
        256,
        3,
    );
    try store.saveOwnedFields(
        plan.tiles[0],
        1,
        1,
        &.{.{ .values = values[0..1], .components_per_cell = 1 }},
    );
    try store.publishGeneration(
        plan,
        &.{.{ .values = values[0..1], .components_per_cell = 1 }},
    );
    try std.testing.expectError(
        error.TileStateFieldShapeMismatch,
        store.validateGeneration(
            plan,
            &.{.{ .values = &values, .components_per_cell = 2 }},
        ),
    );
}

test "generation validation consumes every tile record without mutating state" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    var plan = try spatial_grid.TilePlan.init(
        std.testing.allocator,
        2,
        3,
        1,
        2,
        2,
    );
    defer plan.deinit();
    const source = [_]f64{ 11, 12, 13, 14, 15, 16 };
    var destination = [_]f64{-99} ** source.len;
    const store = try FileStore.initGeneration(
        std.testing.allocator,
        std.testing.io,
        temporary.dir,
        256,
        27,
    );
    for (plan.tiles) |tile| try store.saveOwnedFields(
        tile,
        plan.grid_row_count,
        plan.grid_column_count,
        &.{.{
            .name = "surface_temperature_k",
            .values = &source,
            .components_per_cell = 1,
        }},
    );
    try store.publishGeneration(
        plan,
        &.{.{
            .name = "surface_temperature_k",
            .values = &source,
            .components_per_cell = 1,
        }},
    );
    try store.validateGeneration(
        plan,
        &.{.{
            .name = "surface_temperature_k",
            .values = &destination,
            .components_per_cell = 1,
        }},
    );
    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{-99} ** source.len),
        &destination,
    );

    const corrupt_name = try tileFileName(
        std.testing.allocator,
        plan.tiles[plan.tiles.len - 1].z_order_index,
    );
    defer std.testing.allocator.free(corrupt_name);
    var corrupt_file = try temporary.dir.createFile(
        std.testing.io,
        corrupt_name,
        .{},
    );
    try corrupt_file.writeStreamingAll(std.testing.io, "ECOSTILE");
    corrupt_file.close(std.testing.io);
    try std.testing.expectError(
        error.EndOfStream,
        store.validateGeneration(
            plan,
            &.{.{
                .name = "surface_temperature_k",
                .values = &destination,
                .components_per_cell = 1,
            }},
        ),
    );
}

test "tile load rejects a stale file from another generation" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    var plan = try spatial_grid.TilePlan.init(
        std.testing.allocator,
        1,
        1,
        1,
        1,
        2,
    );
    defer plan.deinit();
    const source = [_]f64{12};
    var destination = [_]f64{0};
    const generation_one = try FileStore.initGeneration(
        std.testing.allocator,
        std.testing.io,
        temporary.dir,
        256,
        1,
    );
    try generation_one.saveOwnedFields(
        plan.tiles[0],
        1,
        1,
        &.{.{ .values = &source, .components_per_cell = 1 }},
    );
    const generation_two = try FileStore.initGeneration(
        std.testing.allocator,
        std.testing.io,
        temporary.dir,
        256,
        2,
    );
    try std.testing.expectError(
        error.TileStateGenerationMismatch,
        generation_two.loadTileFieldsInto(
            plan,
            plan.tiles[0],
            &.{.{ .values = &destination, .components_per_cell = 1 }},
        ),
    );
}

test "tile halo is reconstructed from the neighboring authoritative interior" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    var plan = try spatial_grid.TilePlan.init(
        std.testing.allocator,
        4,
        6,
        2,
        2,
        2,
    );
    defer plan.deinit();
    const source = try std.testing.allocator.alloc(f64, 4 * 6);
    defer std.testing.allocator.free(source);
    @memset(source, 1);
    const destination = try std.testing.allocator.alloc(f64, 4 * 6);
    defer std.testing.allocator.free(destination);
    @memset(destination, 0);
    const store = try FileStore.init(
        std.testing.allocator,
        std.testing.io,
        temporary.dir,
        256,
    );
    for (plan.tiles) |tile| try store.saveOwnedFields(
        tile,
        4,
        6,
        &.{.{ .values = source, .components_per_cell = 1 }},
    );

    const target = plan.tiles[0];
    var neighbor: ?spatial_grid.Tile = null;
    for (plan.tiles) |tile| {
        if (tile.owned_west_column == target.owned_east_column_exclusive and
            tile.owned_north_row == target.owned_north_row)
        {
            neighbor = tile;
            break;
        }
    }
    const east_neighbor = neighbor orelse return error.MissingTestNeighbor;
    const neighbor_cell =
        east_neighbor.owned_north_row * 6 + east_neighbor.owned_west_column;
    source[neighbor_cell] = 99;
    try store.saveOwnedFields(
        east_neighbor,
        4,
        6,
        &.{.{ .values = source, .components_per_cell = 1 }},
    );
    try store.loadTileFieldsInto(
        plan,
        target,
        &.{.{ .values = destination, .components_per_cell = 1 }},
    );
    try std.testing.expectEqual(@as(f64, 99), destination[neighbor_cell]);
}

test "Morton tile fields preserve runtime components per grid cell" {
    var plan = try spatial_grid.TilePlan.init(
        std.testing.allocator,
        3,
        4,
        2,
        2,
        2,
    );
    defer plan.deinit();
    const tile = plan.tiles[0];
    const components_per_cell: usize = 7;
    const value_count = 3 * 4 * components_per_cell;
    const source = try std.testing.allocator.alloc(f64, value_count);
    defer std.testing.allocator.free(source);
    for (source, 0..) |*value, index| value.* = @floatFromInt(index + 1);
    const destination = try std.testing.allocator.alloc(f64, value_count);
    defer std.testing.allocator.free(destination);
    @memset(destination, 0);
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try writeOwnedFields(
        std.testing.allocator,
        &bytes.writer,
        tile,
        3,
        4,
        &.{.{
            .values = source,
            .components_per_cell = components_per_cell,
        }},
    );
    var reader: std.Io.Reader = .fixed(bytes.written());
    try readOwnedFieldsInto(
        std.testing.allocator,
        &reader,
        tile,
        tile,
        3,
        4,
        &.{.{
            .values = destination,
            .components_per_cell = components_per_cell,
        }},
    );
    const loaded_cells = try tile.ownedCellIndicesZOrder(
        std.testing.allocator,
        3,
        4,
    );
    defer std.testing.allocator.free(loaded_cells);
    for (loaded_cells) |cell| {
        const first = cell * components_per_cell;
        try std.testing.expectEqualSlices(
            f64,
            source[first..][0..components_per_cell],
            destination[first..][0..components_per_cell],
        );
    }
}
