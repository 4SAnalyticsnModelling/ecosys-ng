const std = @import("std");
const execution_calendar_date = @import("execution_calendar_date.zig");

pub const Delimiter = enum {
    comma,
    tab,
    space,
    pipe,

    pub fn byte(self: Delimiter) u8 {
        return switch (self) {
            .comma => ',',
            .tab => '\t',
            .space => ' ',
            .pipe => '|',
        };
    }
};

pub const Variable = struct {
    name: []const u8,
    unit: []const u8,
};

pub const Timestamp = struct {
    year: i32,
    day_of_year: u16,
    month: u8,
    day: u8,
    hour: u8,
};

pub const Record = struct {
    timestamp: Timestamp,
    grid_column: usize,
    grid_row: usize,
    values: []const f64,
};

/// Writes a stable, allocation-free heading for only the variables selected
/// by the runtime output editor.
pub fn writeHeader(writer: *std.Io.Writer, variables: []const Variable, enabled: []const bool, delimiter: Delimiter) !void {
    if (variables.len != enabled.len) return error.OutputSelectionDimensionMismatch;
    const separator = delimiter.byte();
    const fixed = [_][]const u8{ "year", "day_of_year", "month", "day", "hour", "grid_column", "grid_row" };
    for (fixed, 0..) |heading, index| {
        if (index != 0) try writer.writeByte(separator);
        try writer.writeAll(heading);
    }
    for (variables, enabled) |variable, selected| {
        if (!selected) continue;
        try validateLabel(variable.name, separator);
        try validateLabel(variable.unit, separator);
        try writer.writeByte(separator);
        try writer.print("{s}[{s}]", .{ variable.name, variable.unit });
    }
    try writer.writeByte('\n');
}

/// Streams one selected record without assembling a second output row. Any
/// NaN or infinity aborts before that value can silently enter an output file.
pub fn writeRecord(writer: *std.Io.Writer, record: Record, enabled: []const bool, delimiter: Delimiter) !void {
    if (record.values.len != enabled.len) return error.OutputSelectionDimensionMismatch;
    try validateTimestamp(record.timestamp);
    if (record.grid_column == 0 or record.grid_row == 0)
        return error.InvalidOutputGridCoordinate;
    for (record.values, enabled) |value, selected|
        if (selected and !std.math.isFinite(value))
            return error.NonFiniteOutputValue;
    const separator = delimiter.byte();
    try writer.print("{d}{c}{d}{c}{d}{c}{d}{c}{d}{c}{d}{c}{d}", .{ record.timestamp.year, separator, record.timestamp.day_of_year, separator, record.timestamp.month, separator, record.timestamp.day, separator, record.timestamp.hour, separator, record.grid_column, separator, record.grid_row });
    for (record.values, enabled) |value, selected| {
        if (!selected) continue;
        try writer.writeByte(separator);
        try writer.print("{e}", .{value});
    }
    try writer.writeByte('\n');
}

/// FOUTS-compatible per-cell/year stem, without a compile-time grid width.
pub fn buildCellFileName(allocator: std.mem.Allocator, grid_column: usize, grid_row: usize, year: i32, editor_name: []const u8) ![]u8 {
    if (grid_column == 0 or grid_row == 0 or year <= 0 or year > 9999 or
        !safeEditorName(editor_name))
        return error.InvalidOutputFileName;
    if (hasTextExtension(editor_name))
        return std.fmt.allocPrint(allocator, "{d:0>2}{d:0>2}0{d:0>4}{s}", .{ grid_column, grid_row, @as(u32, @intCast(year)), editor_name });
    return std.fmt.allocPrint(allocator, "{d:0>2}{d:0>2}0{d:0>4}{s}.txt", .{ grid_column, grid_row, @as(u32, @intCast(year)), editor_name });
}

/// FOUTP-compatible per-plant/year stem. Species is runtime-sized and
/// one-based, so populations are not restricted to the source model's five.
pub fn buildPlantFileName(allocator: std.mem.Allocator, grid_column: usize, grid_row: usize, one_based_species: usize, year: i32, editor_name: []const u8) ![]u8 {
    if (grid_column == 0 or grid_row == 0 or one_based_species == 0 or
        year <= 0 or year > 9999 or !safeEditorName(editor_name))
        return error.InvalidOutputFileName;
    if (hasTextExtension(editor_name))
        return std.fmt.allocPrint(allocator, "{d:0>2}{d:0>2}{d}{d:0>4}{s}", .{ grid_column, grid_row, one_based_species, @as(u32, @intCast(year)), editor_name });
    return std.fmt.allocPrint(allocator, "{d:0>2}{d:0>2}{d}{d:0>4}{s}.txt", .{ grid_column, grid_row, one_based_species, @as(u32, @intCast(year)), editor_name });
}

fn hasTextExtension(name: []const u8) bool {
    return name.len > 4 and std.ascii.eqlIgnoreCase(name[name.len - 4 ..], ".txt");
}

fn safeEditorName(name: []const u8) bool {
    if (name.len == 0 or std.mem.indexOf(u8, name, "..") != null) return false;
    if (name[0] == ' ' or name[name.len - 1] == ' ' or name[name.len - 1] == '.')
        return false;
    for (name) |byte| {
        if (byte == 0 or byte < 0x20 or
            std.mem.indexOfScalar(u8, "/\\<>:\"|?*", byte) != null)
            return false;
    }
    return true;
}

fn validateTimestamp(timestamp: Timestamp) !void {
    if (timestamp.year <= 0 or timestamp.year > 9999 or timestamp.month == 0 or
        timestamp.month > 12 or timestamp.hour > 23)
        return error.InvalidOutputTimestamp;
    const expected = execution_calendar_date.dayOfYear(.{
        .day = timestamp.day,
        .month = timestamp.month,
        .year = @intCast(timestamp.year),
    }) catch return error.InvalidOutputTimestamp;
    if (expected != timestamp.day_of_year)
        return error.InvalidOutputTimestamp;
}

fn validateLabel(label: []const u8, delimiter: u8) !void {
    if (label.len == 0 or std.mem.indexOfScalar(u8, label, delimiter) != null or std.mem.indexOfAny(u8, label, "\r\n[]") != null) return error.InvalidOutputLabel;
}

test "selected output record streams headings units and finite values" {
    const variables = [_]Variable{
        .{ .name = "runoff", .unit = "mm" },
        .{ .name = "soil_temperature", .unit = "degC" },
        .{ .name = "water_table_depth", .unit = "m" },
    };
    const enabled = [_]bool{ true, false, true };
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try writeHeader(&bytes.writer, &variables, &enabled, .pipe);
    try writeRecord(&bytes.writer, .{ .timestamp = .{ .year = 2001, .day_of_year = 32, .month = 2, .day = 1, .hour = 5 }, .grid_column = 12, .grid_row = 3, .values = &.{ 1.25, 99, 2.5 } }, &enabled, .pipe);
    try std.testing.expectEqualStrings("year|day_of_year|month|day|hour|grid_column|grid_row|runoff[mm]|water_table_depth[m]\n2001|32|2|1|5|12|3|1.25e0|2.5e0\n", bytes.written());
}

test "output writer rejects nonfinite selected values" {
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try std.testing.expectError(error.NonFiniteOutputValue, writeRecord(&bytes.writer, .{ .timestamp = .{ .year = 2001, .day_of_year = 1, .month = 1, .day = 1, .hour = 0 }, .grid_column = 1, .grid_row = 1, .values = &.{std.math.nan(f64)} }, &.{true}, .comma));
}

test "cell output name retains source coordinate and year convention" {
    const name = try buildCellFileName(std.testing.allocator, 3, 7, 2001, "hourly_water.txt");
    defer std.testing.allocator.free(name);
    try std.testing.expectEqualStrings("030702001hourly_water.txt", name);
}

test "extensionless editor input still produces a text output filename" {
    const name = try buildCellFileName(std.testing.allocator, 1, 2, 2001, "hourly_water");
    defer std.testing.allocator.free(name);
    try std.testing.expectEqualStrings("010202001hourly_water.txt", name);
}

test "plant output name preserves source convention without a species limit" {
    const ordinary = try buildPlantFileName(std.testing.allocator, 3, 7, 1, 2001, "hourly_carbon.txt");
    defer std.testing.allocator.free(ordinary);
    try std.testing.expectEqualStrings("030712001hourly_carbon.txt", ordinary);
    const many_species = try buildPlantFileName(std.testing.allocator, 3, 7, 127, 2001, "hourly_carbon.txt");
    defer std.testing.allocator.free(many_species);
    try std.testing.expectEqualStrings("03071272001hourly_carbon.txt", many_species);
    try std.testing.expectError(error.InvalidOutputFileName, buildPlantFileName(std.testing.allocator, 3, 7, 0, 2001, "hourly_carbon.txt"));
}

test "output timestamps reject impossible or inconsistent calendar values" {
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    inline for (.{
        Timestamp{ .year = 2001, .day_of_year = 60, .month = 2, .day = 29, .hour = 0 },
        Timestamp{ .year = 2001, .day_of_year = 121, .month = 4, .day = 31, .hour = 0 },
        Timestamp{ .year = 2001, .day_of_year = 2, .month = 1, .day = 1, .hour = 0 },
        Timestamp{ .year = 0, .day_of_year = 1, .month = 1, .day = 1, .hour = 0 },
    }) |timestamp| try std.testing.expectError(
        error.InvalidOutputTimestamp,
        writeRecord(
            &bytes.writer,
            .{
                .timestamp = timestamp,
                .grid_column = 1,
                .grid_row = 1,
                .values = &.{1},
            },
            &.{true},
            .tab,
        ),
    );
    try validateTimestamp(.{
        .year = 2000,
        .day_of_year = 60,
        .month = 2,
        .day = 29,
        .hour = 23,
    });
    try validateTimestamp(.{
        .year = 1900,
        .day_of_year = 60,
        .month = 2,
        .day = 29,
        .hour = 23,
    });
}

test "output file builders reject invalid coordinates and editor names" {
    inline for (.{
        "",
        "../hourly",
        "subdir/hourly",
        "subdir\\hourly",
        " hourly",
        "hourly ",
        "hourly.",
        "hourly:data",
        "hourly|data",
        "hourly?data",
    }) |name| try std.testing.expectError(
        error.InvalidOutputFileName,
        buildCellFileName(std.testing.allocator, 1, 1, 2001, name),
    );
    try std.testing.expectError(
        error.InvalidOutputFileName,
        buildCellFileName(std.testing.allocator, 0, 1, 2001, "hourly"),
    );
    try std.testing.expectError(
        error.InvalidOutputFileName,
        buildPlantFileName(std.testing.allocator, 1, 0, 1, 2001, "hourly"),
    );
}
