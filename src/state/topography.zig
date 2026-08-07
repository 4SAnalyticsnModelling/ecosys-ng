const std = @import("std");
const delimited_input = @import("../io/input/delimited_input.zig");
const Domain = @import("../driver/runscript.zig").Domain;

pub const LandscapeUnit = struct {
    west_column: usize,
    north_row: usize,
    east_column: usize,
    south_row: usize,
    compass_aspect_degrees: f64,
    geometric_aspect_degrees: f64,
    slope_degrees: f64,
    unused_slope_input: f64,
    initial_snowpack_depth_m: f64,
    soil_profile_file: []const u8,

    pub fn contains(self: LandscapeUnit, column: usize, row: usize) bool {
        return column >= self.west_column and column <= self.east_column and
            row >= self.north_row and row <= self.south_row;
    }
};

pub const Topography = struct {
    allocator: std.mem.Allocator,
    units: []LandscapeUnit,

    pub fn deinit(self: *Topography) void {
        for (self.units) |unit| self.allocator.free(unit.soil_profile_file);
        self.allocator.free(self.units);
        self.* = undefined;
    }

    pub fn validateCoverage(self: Topography, domain: Domain, allocator: std.mem.Allocator) !void {
        const map = try self.buildCellUnitMap(domain, allocator);
        allocator.free(map);
    }

    /// Resolves exactly one landscape unit at a global one-based grid
    /// coordinate. Separate per-cell topography files still obey the same
    /// overlap and coverage invariants as a single domain-wide file.
    pub fn unitForCell(self: Topography, column: usize, row: usize) !usize {
        if (column == 0 or row == 0) return error.InvalidTopographyCoordinate;
        var found: ?usize = null;
        for (self.units, 0..) |unit, index| {
            if (!unit.contains(column, row)) continue;
            if (found != null) return error.OverlappingTopographyUnits;
            found = index;
        }
        return found orelse error.IncompleteTopographyCoverage;
    }

    /// Returns one landscape-unit index per selected cell in row-major order.
    /// Overlaps and gaps are fatal because either would make state provenance
    /// ambiguous.
    pub fn buildCellUnitMap(self: Topography, domain: Domain, allocator: std.mem.Allocator) ![]usize {
        const columns = try domain.columns();
        const rows = try domain.rows();
        const cell_count = try std.math.mul(usize, columns, rows);
        const unit_by_cell = try allocator.alloc(usize, cell_count);
        errdefer allocator.free(unit_by_cell);
        @memset(unit_by_cell, std.math.maxInt(usize));
        for (self.units, 0..) |unit, unit_index| {
            const west = @max(unit.west_column, domain.west_column);
            const east = @min(unit.east_column, domain.east_column);
            const north = @max(unit.north_row, domain.north_row);
            const south = @min(unit.south_row, domain.south_row);
            if (west > east or north > south) continue;
            var column = west;
            while (column <= east) : (column += 1) {
                var row = north;
                while (row <= south) : (row += 1) {
                    const index = (row - domain.north_row) * columns + (column - domain.west_column);
                    if (unit_by_cell[index] != std.math.maxInt(usize)) return error.OverlappingTopographyUnits;
                    unit_by_cell[index] = unit_index;
                }
            }
        }
        for (unit_by_cell, 0..) |unit_index, index| {
            if (unit_index == std.math.maxInt(usize)) {
                std.log.err("topography does not cover selected domain cell index={d}", .{index});
                return error.IncompleteTopographyCoverage;
            }
        }
        return unit_by_cell;
    }

    pub fn commonSoilProfileFile(self: Topography, unit_by_cell: []const usize) ![]const u8 {
        if (unit_by_cell.len == 0) return error.EmptyGrid;
        const first_name = self.units[unit_by_cell[0]].soil_profile_file;
        for (unit_by_cell[1..]) |unit_index| {
            if (!std.mem.eql(u8, first_name, self.units[unit_index].soil_profile_file)) {
                return error.MultipleSoilProfilesRequireCatalog;
            }
        }
        return first_name;
    }
};

pub fn parse(allocator: std.mem.Allocator, source: []const u8) !Topography {
    var records = std.mem.splitScalar(u8, source, '\n');
    var units: std.ArrayList(LandscapeUnit) = .empty;
    defer {
        for (units.items) |unit| allocator.free(unit.soil_profile_file);
        units.deinit(allocator);
    }
    while (nextNonemptyRecord(&records)) |record| {
        if (hasEmptyExplicitField(record)) return error.EmptyTopographyRecordValue;
        var fields = delimited_input.recordTokens(record);
        var unit: LandscapeUnit = undefined;
        unit.west_column = try number(usize, &fields);
        unit.north_row = try number(usize, &fields);
        unit.east_column = try number(usize, &fields);
        unit.south_row = try number(usize, &fields);
        unit.compass_aspect_degrees = try number(f64, &fields);
        unit.slope_degrees = try number(f64, &fields);
        unit.unused_slope_input = try number(f64, &fields);
        unit.initial_snowpack_depth_m = try number(f64, &fields);
        if (fields.next() != null) return error.TrailingTopographyRecordData;
        if (unit.west_column == 0 or unit.north_row == 0 or unit.east_column < unit.west_column or unit.south_row < unit.north_row) return error.InvalidLandscapeUnit;
        if (!std.math.isFinite(unit.compass_aspect_degrees) or !std.math.isFinite(unit.slope_degrees) or !std.math.isFinite(unit.initial_snowpack_depth_m)) return error.NonFiniteTopography;
        unit.geometric_aspect_degrees = 450.0 - unit.compass_aspect_degrees;
        if (unit.geometric_aspect_degrees >= 360.0) unit.geometric_aspect_degrees -= 360.0;
        const soil_record = nextNonemptyRecord(&records) orelse return error.MissingSoilProfileFilename;
        if (hasEmptyExplicitField(soil_record)) return error.EmptyTopographyRecordValue;
        var soil_fields = delimited_input.recordTokens(soil_record);
        unit.soil_profile_file = try allocator.dupe(u8, soil_fields.next() orelse return error.MissingSoilProfileFilename);
        errdefer allocator.free(unit.soil_profile_file);
        if (soil_fields.next() != null) return error.InvalidSoilProfileFilenameRecord;
        try units.append(allocator, unit);
    }
    if (units.items.len == 0) return error.EmptyTopography;
    return .{ .allocator = allocator, .units = try units.toOwnedSlice(allocator) };
}

fn nextNonemptyRecord(records: anytype) ?[]const u8 {
    while (records.next()) |record| {
        var fields = delimited_input.recordTokens(record);
        if (fields.next() != null) return record;
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

fn number(comptime T: type, tokens: anytype) !T {
    const text = tokens.next() orelse return error.IncompleteTopographyRecord;
    return switch (@typeInfo(T)) {
        .float => std.fmt.parseFloat(T, text),
        .int => std.fmt.parseInt(T, text, 10),
        else => @compileError("unsupported topography number type"),
    };
}

test "parse self-contained topography and validate selected-domain coverage" {
    const source = "1 1 2 2 94.6 2 0 0\nsoil_profile\n";
    var topography = try parse(std.testing.allocator, source);
    defer topography.deinit();
    try std.testing.expectEqual(@as(usize, 1), topography.units.len);
    try std.testing.expectEqualStrings("soil_profile", topography.units[0].soil_profile_file);
    try std.testing.expectApproxEqAbs(@as(f64, 355.4), topography.units[0].geometric_aspect_degrees, 1.0e-12);
    try topography.validateCoverage(.{
        .west_column = 1,
        .north_row = 1,
        .east_column = 1,
        .south_row = 1,
    }, std.testing.allocator);
    const full_domain = Domain{ .west_column = 1, .north_row = 1, .east_column = 2, .south_row = 2 };
    const unit_by_cell = try topography.buildCellUnitMap(full_domain, std.testing.allocator);
    defer std.testing.allocator.free(unit_by_cell);
    try std.testing.expectEqual(@as(usize, 4), unit_by_cell.len);
    try std.testing.expectEqualStrings("soil_profile", try topography.commonSoilProfileFile(unit_by_cell));
}

test "single-cell lookup rejects gaps and overlaps" {
    var topography = try parse(
        std.testing.allocator,
        "1 1 2 1 90 1 0 0\nsoil_a\n2 1 3 1 180 2 0 0\nsoil_b\n",
    );
    defer topography.deinit();
    try std.testing.expectEqual(@as(usize, 0), try topography.unitForCell(1, 1));
    try std.testing.expectError(error.OverlappingTopographyUnits, topography.unitForCell(2, 1));
    try std.testing.expectError(error.IncompleteTopographyCoverage, topography.unitForCell(4, 1));
}

test "topography accepts whole-line and trailing comments without joining records" {
    var topography = try parse(
        std.testing.allocator,
        \\# Landscape unit followed by its compulsory soil-profile record.
        \\1,1,1,1,90,2,0,0 # north-facing test cell
        \\# The filename may have no extension.
        \\soil_profile_a # cell-specific soil input
        \\   # comment with leading whitespace
        \\2|1|2|1|180|3|0|0
        \\soil_profile_b
        ,
    );
    defer topography.deinit();
    try std.testing.expectEqual(@as(usize, 2), topography.units.len);
    try std.testing.expectEqualStrings(
        "soil_profile_a",
        topography.units[0].soil_profile_file,
    );
    try std.testing.expectEqualStrings(
        "soil_profile_b",
        topography.units[1].soil_profile_file,
    );
}

test "topography comments cannot hide missing or extra required values" {
    try std.testing.expectError(
        error.IncompleteTopographyRecord,
        parse(
            std.testing.allocator,
            "1 1 1 1 90 2 0 # missing snow depth\nsoil_profile\n",
        ),
    );
    try std.testing.expectError(
        error.TrailingTopographyRecordData,
        parse(
            std.testing.allocator,
            "1 1 1 1 90 2 0 0 99 # extra value\nsoil_profile\n",
        ),
    );
}

test "topography records reject empty explicit delimiter fields" {
    try std.testing.expectError(
        error.EmptyTopographyRecordValue,
        parse(std.testing.allocator, "1,,1,1,90,2,0,0\nsoil_profile\n"),
    );
    try std.testing.expectError(
        error.EmptyTopographyRecordValue,
        parse(std.testing.allocator, "1 1 1 1 90 2 0 0\n,soil_profile\n"),
    );
    try std.testing.expectError(
        error.EmptyTopographyRecordValue,
        parse(std.testing.allocator, "1,1,1,1,90,2,0,0\nsoil_profile,\n"),
    );
}
