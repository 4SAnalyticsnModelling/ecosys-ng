const std = @import("std");

pub const Progress = struct {
    execution_day: u64,
    year: i32,
};

/// Structured ecosys-ng equivalent of EXEC format 212. This is a progress
/// diagnostic, not a scientific output record.
pub fn write(writer: *std.Io.Writer, progress: Progress) !void {
    if (progress.execution_day == 0)
        return error.InvalidExecutionProgressDay;
    try writer.print(
        "now_executing_day={d}\tyear={d}\n",
        .{ progress.execution_day, progress.year },
    );
}

test "EXEC cadence progress preserves day and year ordering" {
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try write(&bytes.writer, .{
        .execution_day = 81,
        .year = 1998,
    });
    try std.testing.expectEqualStrings(
        "now_executing_day=81\tyear=1998\n",
        bytes.written(),
    );
}

test "signed runtime year remains explicit" {
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try write(&bytes.writer, .{
        .execution_day = std.math.maxInt(u64),
        .year = -1,
    });
    try std.testing.expectEqualStrings(
        "now_executing_day=18446744073709551615\tyear=-1\n",
        bytes.written(),
    );
}

test "invalid day writes no partial progress message" {
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try std.testing.expectError(
        error.InvalidExecutionProgressDay,
        write(&bytes.writer, .{
            .execution_day = 0,
            .year = 2001,
        }),
    );
    try std.testing.expectEqual(@as(usize, 0), bytes.written().len);
}
