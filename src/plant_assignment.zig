const std = @import("std");
const delimited_input = @import("delimited_input.zig");
const Domain = @import("runscript.zig").Domain;

pub const SpeciesAssignment = struct {
    species_file: []const u8,
    management_file: []const u8,
};

pub const Unit = struct {
    west_column: usize,
    north_row: usize,
    east_column: usize,
    south_row: usize,
    species: []SpeciesAssignment,
};

pub const Assignments = struct {
    allocator: std.mem.Allocator,
    units: []Unit,

    pub fn deinit(self: *Assignments) void {
        for (self.units) |unit| {
            for (unit.species) |species| {
                self.allocator.free(species.species_file);
                self.allocator.free(species.management_file);
            }
            self.allocator.free(unit.species);
        }
        self.allocator.free(self.units);
        self.* = undefined;
    }

    pub fn buildCellUnitMap(self: Assignments, allocator: std.mem.Allocator, domain: Domain, species_capacity: usize) ![]usize {
        const columns = try domain.columns();
        const rows = try domain.rows();
        const map = try allocator.alloc(usize, try std.math.mul(usize, columns, rows));
        errdefer allocator.free(map);
        @memset(map, std.math.maxInt(usize));
        for (self.units, 0..) |unit, unit_index| {
            if (unit.species.len > species_capacity) {
                std.log.err("plant assignment exceeds runtime species capacity: assigned={d} capacity={d} unit={d}", .{ unit.species.len, species_capacity, unit_index });
                return error.PlantSpeciesCapacityExceeded;
            }
            const west = @max(unit.west_column, domain.west_column);
            const east = @min(unit.east_column, domain.east_column);
            const north = @max(unit.north_row, domain.north_row);
            const south = @min(unit.south_row, domain.south_row);
            if (west > east or north > south) continue;
            for (north..south + 1) |row| for (west..east + 1) |column| {
                const index = (row - domain.north_row) * columns + (column - domain.west_column);
                if (map[index] != std.math.maxInt(usize)) return error.OverlappingPlantAssignmentUnits;
                map[index] = unit_index;
            };
        }
        for (map) |unit_index| if (unit_index == std.math.maxInt(usize)) return error.IncompletePlantAssignmentCoverage;
        return map;
    }
};

pub fn parse(allocator: std.mem.Allocator, source: []const u8, ecosystem_type: i32) !Assignments {
    if (ecosystem_type < 0 or ecosystem_type > 99) return error.InvalidEcosystemTypeForSpeciesFilename;
    var records = std.mem.splitScalar(u8, source, '\n');
    var units: std.ArrayList(Unit) = .empty;
    defer {
        for (units.items) |unit| {
            for (unit.species) |species| {
                allocator.free(species.species_file);
                allocator.free(species.management_file);
            }
            allocator.free(unit.species);
        }
        units.deinit(allocator);
    }
    while (try nextRecordOrNull(&records)) |header_record| {
        var header = delimited_input.recordTokens(header_record);
        const west = try number(usize, &header);
        const north = try number(usize, &header);
        const east = try number(usize, &header);
        const south = try number(usize, &header);
        const species_count = try number(usize, &header);
        if (header.next() != null) return error.TrailingPlantAssignmentHeaderData;
        if (west == 0 or north == 0 or east < west or south < north) return error.InvalidPlantAssignmentUnit;
        const species = try allocator.alloc(SpeciesAssignment, species_count);
        var allocated: usize = 0;
        errdefer {
            for (species[0..allocated]) |entry| {
                allocator.free(entry.species_file);
                allocator.free(entry.management_file);
            }
            allocator.free(species);
        }
        if (species_count > 0) {
            const names_record = try nextRecordOrNull(&records) orelse return error.MissingPlantAssignmentNames;
            var names = delimited_input.recordTokens(names_record);
            while (allocated < species_count) : (allocated += 1) {
                const species_base = names.next() orelse return error.IncompletePlantAssignmentNames;
                const management_name = names.next() orelse return error.IncompletePlantAssignmentNames;
                species[allocated] = .{
                    .species_file = try resolvedSpeciesName(allocator, species_base, ecosystem_type),
                    .management_file = allocator.dupe(u8, management_name) catch |err| {
                        allocator.free(species[allocated].species_file);
                        return err;
                    },
                };
            }
            if (names.next() != null) return error.TrailingPlantAssignmentNames;
        }
        try units.append(allocator, .{
            .west_column = west,
            .north_row = north,
            .east_column = east,
            .south_row = south,
            .species = species,
        });
    }
    if (units.items.len == 0) return error.EmptyPlantAssignments;
    return .{ .allocator = allocator, .units = try units.toOwnedSlice(allocator) };
}

fn resolvedSpeciesName(allocator: std.mem.Allocator, base: []const u8, ecosystem_type: i32) ![]u8 {
    if (ecosystem_type == 0) return allocator.dupe(u8, base);
    if (base.len < 4) return error.SpeciesFilenameTooShort;
    return if (ecosystem_type < 10)
        std.fmt.allocPrint(allocator, "{s}0{d}", .{ base[0..4], ecosystem_type })
    else
        std.fmt.allocPrint(allocator, "{s}{d}", .{ base[0..4], ecosystem_type });
}

fn nextRecordOrNull(records: anytype) !?[]const u8 {
    while (records.next()) |record| {
        if (hasEmptyExplicitField(record))
            return error.EmptyPlantAssignmentRecordValue;
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

fn number(comptime T: type, tokens: anytype) !T {
    return std.fmt.parseInt(T, tokens.next() orelse return error.IncompletePlantAssignmentHeader, 10);
}

test "parse runtime-sized self-contained plant assignments" {
    const allocator = std.testing.allocator;
    const source = @import("test_fixtures.zig").plant_assignment_source;
    var assignments = try parse(allocator, source, 33);
    defer assignments.deinit();
    try std.testing.expectEqual(@as(usize, 1), assignments.units[0].species.len);
    try std.testing.expectEqualStrings("maiz33", assignments.units[0].species[0].species_file);
    try std.testing.expectEqualStrings("management", assignments.units[0].species[0].management_file);
    const map = try assignments.buildCellUnitMap(allocator, .{ .west_column = 1, .north_row = 1, .east_column = 1, .south_row = 1 }, 17);
    defer allocator.free(map);
    try std.testing.expectEqual(@as(usize, 0), map[0]);
}

test "plant assignment comments preserve runtime-sized record pairs" {
    const allocator = std.testing.allocator;
    var assignments = try parse(
        allocator,
        \\# Bounds and runtime species count.
        \\1,1,1,1,2 # two species in the western cell
        \\   # Exactly two species/management filename pairs follow.
        \\maize|maize_management|soybean|soybean_management
        \\# A neighboring cell can have a different runtime count.
        \\2 1 2 1 1
        \\wheat wheat_management # extensionless inputs
        \\# End of assignments.
    ,
        0,
    );
    defer assignments.deinit();
    try std.testing.expectEqual(@as(usize, 2), assignments.units.len);
    try std.testing.expectEqual(
        @as(usize, 2),
        assignments.units[0].species.len,
    );
    try std.testing.expectEqualStrings(
        "soybean",
        assignments.units[0].species[1].species_file,
    );
    try std.testing.expectEqual(
        @as(usize, 1),
        assignments.units[1].species.len,
    );
}

test "plant assignment comments cannot hide short or long records" {
    try std.testing.expectError(
        error.IncompletePlantAssignmentHeader,
        parse(
            std.testing.allocator,
            "1 1 1 1 # missing runtime species count\n",
            0,
        ),
    );
    try std.testing.expectError(
        error.TrailingPlantAssignmentNames,
        parse(
            std.testing.allocator,
            "1 1 1 1 1\nspecies management extra # too many names\n",
            0,
        ),
    );
}

test "plant assignment rejects empty explicit delimiter fields" {
    try std.testing.expectError(
        error.EmptyPlantAssignmentRecordValue,
        parse(std.testing.allocator, "1,1,1,1,1\nsp,,mgmt\n", 0),
    );
    try std.testing.expectError(
        error.EmptyPlantAssignmentRecordValue,
        parse(
            std.testing.allocator,
            "1 1 1 1 1\nsp| |mgmt\n",
            0,
        ),
    );
    try std.testing.expectError(
        error.EmptyPlantAssignmentRecordValue,
        parse(std.testing.allocator, "1 1 1 1 1\nsp mgmt,\n", 0),
    );
}
