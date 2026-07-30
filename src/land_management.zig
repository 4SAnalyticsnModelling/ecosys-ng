const std = @import("std");
const delimited_input = @import("delimited_input.zig");
const Domain = @import("runscript.zig").Domain;

pub const Unit = struct {
    west_column: usize,
    north_row: usize,
    east_column: usize,
    south_row: usize,
    tillage_file: []const u8,
    fertilizer_file: []const u8,
    irrigation_file: []const u8,
};

pub const Assignments = struct {
    allocator: std.mem.Allocator,
    units: []Unit,

    pub fn deinit(self: *Assignments) void {
        for (self.units) |unit| {
            self.allocator.free(unit.tillage_file);
            self.allocator.free(unit.fertilizer_file);
            self.allocator.free(unit.irrigation_file);
        }
        self.allocator.free(self.units);
        self.* = undefined;
    }

    pub fn buildCellUnitMap(self: Assignments, allocator: std.mem.Allocator, domain: Domain) ![]usize {
        const columns = try domain.columns();
        const rows = try domain.rows();
        const map = try allocator.alloc(usize, try std.math.mul(usize, columns, rows));
        errdefer allocator.free(map);
        @memset(map, std.math.maxInt(usize));
        for (self.units, 0..) |unit, unit_index| {
            const west = @max(unit.west_column, domain.west_column);
            const east = @min(unit.east_column, domain.east_column);
            const north = @max(unit.north_row, domain.north_row);
            const south = @min(unit.south_row, domain.south_row);
            if (west > east or north > south) continue;
            for (north..south + 1) |row| for (west..east + 1) |column| {
                const index = (row - domain.north_row) * columns + column - domain.west_column;
                if (map[index] != std.math.maxInt(usize)) return error.OverlappingLandManagementUnits;
                map[index] = unit_index;
            };
        }
        for (map) |unit_index| if (unit_index == std.math.maxInt(usize)) return error.IncompleteLandManagementCoverage;
        return map;
    }
};

pub fn parse(allocator: std.mem.Allocator, source: []const u8) !Assignments {
    var records = std.mem.splitScalar(u8, source, '\n');
    var units: std.ArrayList(Unit) = .empty;
    defer {
        for (units.items) |unit| freeUnit(allocator, unit);
        units.deinit(allocator);
    }
    while (try nextRecord(&records)) |bounds_record| {
        if (hasEmptyExplicitField(bounds_record)) return error.EmptyLandManagementRecordValue;
        var bounds = delimited_input.recordTokens(bounds_record);
        const west = try integer(&bounds);
        const north = try integer(&bounds);
        const east = try integer(&bounds);
        const south = try integer(&bounds);
        if (bounds.next() != null) return error.TrailingLandManagementBounds;
        if (west == 0 or north == 0 or east < west or south < north) return error.InvalidLandManagementUnit;
        const files_record = try nextRecord(&records) orelse return error.MissingLandManagementFiles;
        if (hasEmptyExplicitField(files_record)) return error.EmptyLandManagementRecordValue;
        var names = delimited_input.recordTokens(files_record);
        const tillage = try duplicateNext(allocator, &names);
        errdefer allocator.free(tillage);
        const fertilizer = try duplicateNext(allocator, &names);
        errdefer allocator.free(fertilizer);
        const irrigation = try duplicateNext(allocator, &names);
        errdefer allocator.free(irrigation);
        if (names.next() != null) return error.TrailingLandManagementFiles;
        try units.append(allocator, .{
            .west_column = west,
            .north_row = north,
            .east_column = east,
            .south_row = south,
            .tillage_file = tillage,
            .fertilizer_file = fertilizer,
            .irrigation_file = irrigation,
        });
    }
    if (units.items.len == 0) return error.EmptyLandManagement;
    return .{ .allocator = allocator, .units = try units.toOwnedSlice(allocator) };
}

fn integer(tokens: anytype) !usize {
    return std.fmt.parseUnsigned(usize, tokens.next() orelse return error.IncompleteLandManagementBounds, 10);
}

fn duplicateNext(allocator: std.mem.Allocator, tokens: anytype) ![]u8 {
    return allocator.dupe(u8, tokens.next() orelse return error.IncompleteLandManagementFiles);
}

fn nextRecord(records: anytype) !?[]const u8 {
    while (records.next()) |record| {
        var tokens = delimited_input.recordTokens(record);
        if (tokens.next() != null) return record;
    }
    return null;
}

fn hasEmptyExplicitField(record: []const u8) bool {
    const content = if (std.mem.indexOfScalar(u8, record, '#')) |comment|
        record[0..comment]
    else
        record;
    const trimmed = std.mem.trim(u8, content, " \r");
    if (trimmed.len == 0) return false;

    var field_start: usize = 0;
    var saw_explicit_delimiter = false;
    for (trimmed, 0..) |byte, index| {
        if (byte != ',' and byte != '|' and byte != '\t') continue;
        if (std.mem.trim(u8, trimmed[field_start..index], " \r").len == 0)
            return true;
        field_start = index + 1;
        saw_explicit_delimiter = true;
    }
    return saw_explicit_delimiter and
        std.mem.trim(u8, trimmed[field_start..], " \r").len == 0;
}

fn freeUnit(allocator: std.mem.Allocator, unit: Unit) void {
    allocator.free(unit.tillage_file);
    allocator.free(unit.fertilizer_file);
    allocator.free(unit.irrigation_file);
}

test "parse and map self-contained land management" {
    const allocator = std.testing.allocator;
    const source = @import("test_fixtures.zig").land_management_source;
    var assignments = try parse(allocator, source);
    defer assignments.deinit();
    try std.testing.expectEqualStrings("tillage", assignments.units[0].tillage_file);
    try std.testing.expectEqualStrings("fertilizer", assignments.units[0].fertilizer_file);
    try std.testing.expectEqualStrings("NO", assignments.units[0].irrigation_file);
    const map = try assignments.buildCellUnitMap(allocator, .{ .west_column = 1, .north_row = 1, .east_column = 1, .south_row = 1 });
    defer allocator.free(map);
    try std.testing.expectEqual(@as(usize, 0), map[0]);
}

test "land management accepts pipe and tab delimiters" {
    const allocator = std.testing.allocator;
    var assignments = try parse(allocator, "1|1|2|1\ntillage\tfertilizer|irrigation\n");
    defer assignments.deinit();
    try std.testing.expectEqual(@as(usize, 1), assignments.units.len);
}

test "land management comments do not displace compulsory record pairs" {
    const allocator = std.testing.allocator;
    var assignments = try parse(
        allocator,
        \\# Each unit has one bounds line and one three-file line.
        \\1,1,1,1 # western unit
        \\   # Filenames may be extensionless.
        \\tillage_a|fertilizer_a|NO # irrigation disabled
        \\# A second valid unit may reuse or select other files.
        \\2 1 2 1
        \\tillage_b fertilizer_b irrigation_b
        \\# End of assignments.
        ,
    );
    defer assignments.deinit();
    try std.testing.expectEqual(@as(usize, 2), assignments.units.len);
    try std.testing.expectEqualStrings(
        "tillage_a",
        assignments.units[0].tillage_file,
    );
    try std.testing.expectEqualStrings(
        "NO",
        assignments.units[0].irrigation_file,
    );
    try std.testing.expectEqualStrings(
        "irrigation_b",
        assignments.units[1].irrigation_file,
    );
}

test "land management comments cannot hide short or long records" {
    try std.testing.expectError(
        error.IncompleteLandManagementBounds,
        parse(
            std.testing.allocator,
            "1 1 1 # missing south row\none two three\n",
        ),
    );
    try std.testing.expectError(
        error.TrailingLandManagementFiles,
        parse(
            std.testing.allocator,
            "1 1 1 1\none two three four # extra filename\n",
        ),
    );
}

test "land management rejects empty explicit delimiter fields" {
    try std.testing.expectError(
        error.EmptyLandManagementRecordValue,
        parse(std.testing.allocator, "1,,1,1\none two three\n"),
    );
    try std.testing.expectError(
        error.EmptyLandManagementRecordValue,
        parse(std.testing.allocator, "1 1 1 1\none| |three\n"),
    );
    try std.testing.expectError(
        error.EmptyLandManagementRecordValue,
        parse(std.testing.allocator, "1 1 1 1\none two,\n"),
    );
}
