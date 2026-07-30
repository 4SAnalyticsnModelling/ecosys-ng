const std = @import("std");
const output_record = @import("output_record.zig");
const output_stream_rotation = @import("output_stream_rotation.zig");

/// Runtime-sized, fixed-memory storage for one output editor family.
/// Disabled families allocate nothing. Enabled families own one reusable
/// scientific row and one bounded rotating stream per grid cell.
pub const Bank = struct {
    allocator: std.mem.Allocator,
    cell_count: usize,
    variable_count: usize,
    values: []f64,
    streams: []output_stream_rotation.RotatingStream,
    enabled: bool,

    pub fn init(
        allocator: std.mem.Allocator,
        io: std.Io,
        directory: std.Io.Dir,
        enabled: bool,
        cell_count: usize,
        variable_count: usize,
        stream_buffer_bytes: usize,
        delimiter: output_record.Delimiter,
    ) !Bank {
        if (cell_count == 0 or variable_count == 0) return error.InvalidOutputStreamBankDimensions;
        if (enabled and stream_buffer_bytes == 0)
            return error.InvalidOutputBufferSize;
        if (!enabled) return .{
            .allocator = allocator,
            .cell_count = cell_count,
            .variable_count = variable_count,
            .values = @constCast(&.{}),
            .streams = @constCast(&.{}),
            .enabled = false,
        };
        const value_count = std.math.mul(
            usize,
            cell_count,
            variable_count,
        ) catch return error.OutputStreamBankSizeOverflow;
        const values = try allocator.alloc(f64, value_count);
        errdefer allocator.free(values);
        @memset(values, 0);
        const streams = try allocator.alloc(output_stream_rotation.RotatingStream, cell_count);
        var initialized: usize = 0;
        errdefer {
            for (streams[0..initialized]) |*stream| stream.deinit();
            allocator.free(streams);
        }
        for (streams) |*stream| {
            stream.* = try output_stream_rotation.RotatingStream.init(allocator, io, directory, stream_buffer_bytes, delimiter);
            initialized += 1;
        }
        return .{ .allocator = allocator, .cell_count = cell_count, .variable_count = variable_count, .values = values, .streams = streams, .enabled = true };
    }

    pub fn deinit(self: *Bank) void {
        if (self.enabled) {
            for (self.streams) |*stream| stream.deinit();
            self.allocator.free(self.streams);
            self.allocator.free(self.values);
        }
        self.* = undefined;
    }

    pub fn row(self: *Bank, cell: usize) ![]f64 {
        if (!self.enabled) return error.OutputStreamBankDisabled;
        if (cell >= self.cell_count) return error.OutputStreamBankCellOutOfBounds;
        return self.values[cell * self.variable_count ..][0..self.variable_count];
    }

    pub fn finish(self: *Bank) !void {
        if (!self.enabled) return;
        for (self.streams) |*stream| try stream.finish();
    }
};

test "disabled output bank owns no per-cell storage" {
    var bank = try Bank.init(std.testing.allocator, std.testing.io, undefined, false, 1_000_000, 500, 64, .comma);
    defer bank.deinit();
    try std.testing.expect(!bank.enabled);
    try std.testing.expectEqual(@as(usize, 0), bank.values.len);
    try std.testing.expectEqual(@as(usize, 0), bank.streams.len);
    try std.testing.expectError(error.OutputStreamBankDisabled, bank.row(0));
    try bank.finish();
}

test "enabled output bank exposes stable runtime rows" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    var bank = try Bank.init(std.testing.allocator, std.testing.io, temporary.dir, true, 3, 4, 32, .pipe);
    defer bank.deinit();
    const second = try bank.row(1);
    second[2] = 7;
    try std.testing.expectEqual(@as(f64, 7), bank.values[6]);
    try std.testing.expectError(error.OutputStreamBankCellOutOfBounds, bank.row(3));
    try bank.finish();
}

test "enabled output bank rejects invalid runtime size before allocation" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    try std.testing.expectError(
        error.InvalidOutputBufferSize,
        Bank.init(
            std.testing.allocator,
            std.testing.io,
            temporary.dir,
            true,
            1,
            1,
            0,
            .tab,
        ),
    );
    try std.testing.expectError(
        error.OutputStreamBankSizeOverflow,
        Bank.init(
            std.testing.allocator,
            std.testing.io,
            temporary.dir,
            true,
            std.math.maxInt(usize),
            2,
            64,
            .tab,
        ),
    );
}

test "disabled output bank permits zero buffer without allocating" {
    var bank = try Bank.init(
        std.testing.allocator,
        std.testing.io,
        undefined,
        false,
        2,
        3,
        0,
        .tab,
    );
    defer bank.deinit();
    try std.testing.expectEqual(@as(usize, 0), bank.values.len);
    try std.testing.expectEqual(@as(usize, 0), bank.streams.len);
}
