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
    /// Site longitude in degrees east, from the site file's
    /// `longitude_degrees_east`. This is the physical coordinate the user
    /// entered, not a grid index: an output row must be locatable on the earth
    /// without knowing the grid layout.
    longitude_degrees_east: f64,
    /// Site latitude in degrees north, from `latitude_degrees_north`.
    latitude_degrees_north: f64,
    values: []const f64,
};

/// Writes a stable, allocation-free heading for only the variables selected
/// by the runtime output editor.
pub fn writeHeader(writer: *std.Io.Writer, variables: []const Variable, enabled: []const bool, delimiter: Delimiter) !void {
    if (variables.len != enabled.len) return error.OutputSelectionDimensionMismatch;
    const separator = delimiter.byte();
    const fixed = [_][]const u8{ "year", "day_of_year", "month", "day", "hour", "longitude", "latitude" };
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
    // Physical coordinate ranges, matching `site.zig`'s own validation. A grid
    // index could only be zero-checked; a real coordinate can be checked against
    // the earth, which also catches a caller still passing an index.
    if (!std.math.isFinite(record.longitude_degrees_east) or
        record.longitude_degrees_east < -180 or record.longitude_degrees_east > 180 or
        !std.math.isFinite(record.latitude_degrees_north) or
        record.latitude_degrees_north < -90 or record.latitude_degrees_north > 90)
        return error.InvalidOutputSiteCoordinate;
    for (record.values, enabled) |value, selected|
        if (selected and !std.math.isFinite(value))
            return error.NonFiniteOutputValue;
    const separator = delimiter.byte();
    try writer.print("{d}{c}{d}{c}{d}{c}{d}{c}{d}{c}{d}{c}{d}", .{ record.timestamp.year, separator, record.timestamp.day_of_year, separator, record.timestamp.month, separator, record.timestamp.day, separator, record.timestamp.hour, separator, record.longitude_degrees_east, separator, record.latitude_degrees_north });
    for (record.values, enabled) |value, selected| {
        if (!selected) continue;
        try writer.writeByte(separator);
        try writer.print("{e}", .{value});
    }
    try writer.writeByte('\n');
}

/// Subject of an output file: the whole soil/ecosystem column, or one plant
/// species. The source model encoded this as a digit in the file name, `0` for
/// FOUTS soil output and `1..5` for FOUTP per-plant output.
pub const Subject = union(enum) {
    /// Whole-column soil and ecosystem output, the source model's `0`.
    soil_or_eco,
    /// One plant species, named by its species input file rather than numbered.
    species: []const u8,
};

/// Builds a self-describing output file name:
///
///     lat_<latitude>_lon_<longitude>_<subject>_<year>_<editor>.txt
///
/// for example `lat_45.30_lon_-75.70_soil_or_eco_1998_f25ed1.txt`.
///
/// This replaces the source model's positional stem (`010101998f25ed1.txt`), whose
/// leading digits were grid column, grid row and a species digit. That encoding
/// could not be read without knowing the grid layout, capped species at one digit,
/// and gave no hint of the site's real location. Every
/// (latitude, longitude, subject, year, editor) combination gets its own file, so
/// the name is a complete key for the row set it contains.
///
/// Coordinates are printed to two decimal places, which distinguishes sites about
/// a kilometre apart and keeps the name short. A site file supplying more precision
/// than that would collide, which is why `latitude_degrees_north` and
/// `longitude_degrees_east` are validated to real ranges by the caller.
pub fn buildOutputFileName(
    allocator: std.mem.Allocator,
    latitude_degrees_north: f64,
    longitude_degrees_east: f64,
    subject: Subject,
    year: i32,
    editor_name: []const u8,
) ![]u8 {
    if (year <= 0 or year > 9999 or !safeEditorName(editor_name))
        return error.InvalidOutputFileName;
    if (!std.math.isFinite(latitude_degrees_north) or
        latitude_degrees_north < -90 or latitude_degrees_north > 90 or
        !std.math.isFinite(longitude_degrees_east) or
        longitude_degrees_east < -180 or longitude_degrees_east > 180)
        return error.InvalidOutputFileName;
    const subject_text = switch (subject) {
        .soil_or_eco => "soil_or_eco",
        .species => |name| name,
    };
    // A species name becomes a path component, so it needs the same safety check
    // as the editor name: no separators, no parent traversal, no empty name.
    if (!safeEditorName(subject_text)) return error.InvalidOutputFileName;
    const stem = trimTextExtension(editor_name);
    if (stem.len == 0) return error.InvalidOutputFileName;
    return std.fmt.allocPrint(
        allocator,
        "lat_{d:.2}_lon_{d:.2}_{s}_{d:0>4}_{s}.txt",
        .{ latitude_degrees_north, longitude_degrees_east, subject_text, @as(u32, @intCast(year)), stem },
    );
}

/// Drops a trailing `.txt` so an editor name that already carries one does not
/// produce `..._f25ed1.txt.txt`.
fn trimTextExtension(name: []const u8) []const u8 {
    return if (hasTextExtension(name)) name[0 .. name.len - 4] else name;
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
    try writeRecord(&bytes.writer, .{ .timestamp = .{ .year = 2001, .day_of_year = 32, .month = 2, .day = 1, .hour = 5 }, .longitude_degrees_east = -75.7, .latitude_degrees_north = 45.3, .values = &.{ 1.25, 99, 2.5 } }, &enabled, .pipe);
    try std.testing.expectEqualStrings("year|day_of_year|month|day|hour|longitude|latitude|runoff[mm]|water_table_depth[m]\n2001|32|2|1|5|-75.7|45.3|1.25e0|2.5e0\n", bytes.written());
}

test "SI unit labels with spaces stream under tab and are rejected under space" {
    const variables = [_]Variable{
        .{ .name = "litter_water_vapor_density", .unit = "g m-3" },
        .{ .name = "nitrous_oxide_emission", .unit = "g N m-2 h-1" },
    };
    const enabled = [_]bool{ true, true };
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try writeHeader(&bytes.writer, &variables, &enabled, .tab);
    try std.testing.expectEqualStrings(
        "year\tday_of_year\tmonth\tday\thour\tlongitude\tlatitude\tlitter_water_vapor_density[g m-3]\tnitrous_oxide_emission[g N m-2 h-1]\n",
        bytes.written(),
    );

    // A space-delimited stream cannot carry a unit that contains spaces, so the
    // heading is rejected instead of producing ambiguous columns.
    var space_bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer space_bytes.deinit();
    try std.testing.expectError(
        error.InvalidOutputLabel,
        writeHeader(&space_bytes.writer, &variables, &enabled, .space),
    );
}

test "output writer rejects nonfinite selected values" {
    var bytes: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer bytes.deinit();
    try std.testing.expectError(error.NonFiniteOutputValue, writeRecord(&bytes.writer, .{ .timestamp = .{ .year = 2001, .day_of_year = 1, .month = 1, .day = 1, .hour = 0 }, .longitude_degrees_east = -75.7, .latitude_degrees_north = 45.3, .values = &.{std.math.nan(f64)} }, &.{true}, .comma));
}

test "output name is self describing: lat, lon, subject, year, editor" {
    const name = try buildOutputFileName(std.testing.allocator, 45.3, -75.7, .soil_or_eco, 1998, "f25ed1");
    defer std.testing.allocator.free(name);
    // Replaces the source model's positional stem `010101998f25ed1.txt`, whose
    // leading digits were grid column, grid row and a species digit.
    try std.testing.expectEqualStrings("lat_45.30_lon_-75.70_soil_or_eco_1998_f25ed1.txt", name);
}

test "an editor name that already ends in .txt does not double the extension" {
    const name = try buildOutputFileName(std.testing.allocator, 45.3, -75.7, .soil_or_eco, 1998, "f25ed1.txt");
    defer std.testing.allocator.free(name);
    try std.testing.expectEqualStrings("lat_45.30_lon_-75.70_soil_or_eco_1998_f25ed1.txt", name);
}

test "plant output names carry the species name rather than a digit" {
    const maize = try buildOutputFileName(std.testing.allocator, 45.3, -75.7, .{ .species = "maize" }, 1998, "f25ch1");
    defer std.testing.allocator.free(maize);
    try std.testing.expectEqualStrings("lat_45.30_lon_-75.70_maize_1998_f25ch1.txt", maize);
    // Species are no longer limited to one digit's worth of populations, and two
    // species at the same site and year get distinct files.
    const soybean = try buildOutputFileName(std.testing.allocator, 45.3, -75.7, .{ .species = "soybean" }, 1998, "f25ch1");
    defer std.testing.allocator.free(soybean);
    try std.testing.expectEqualStrings("lat_45.30_lon_-75.70_soybean_1998_f25ch1.txt", soybean);
}

test "southern and western sites keep their coordinate signs" {
    const name = try buildOutputFileName(std.testing.allocator, -33.87, 151.21, .soil_or_eco, 2001, "hourly_water");
    defer std.testing.allocator.free(name);
    try std.testing.expectEqualStrings("lat_-33.87_lon_151.21_soil_or_eco_2001_hourly_water.txt", name);
}

test "each subject, year and editor combination yields a distinct file" {
    // The name is the complete key for the row set it contains, so varying any
    // one component must change it.
    const base = try buildOutputFileName(std.testing.allocator, 45.3, -75.7, .soil_or_eco, 1998, "f25ch1");
    defer std.testing.allocator.free(base);
    inline for (.{
        .{ 45.3, -75.7, Subject.soil_or_eco, 1999, "f25ch1" },
        .{ 45.3, -75.7, Subject.soil_or_eco, 1998, "f25wh1" },
        .{ 46.0, -75.7, Subject.soil_or_eco, 1998, "f25ch1" },
        .{ 45.3, -74.0, Subject.soil_or_eco, 1998, "f25ch1" },
    }) |variant| {
        const other = try buildOutputFileName(std.testing.allocator, variant[0], variant[1], variant[2], variant[3], variant[4]);
        defer std.testing.allocator.free(other);
        try std.testing.expect(!std.mem.eql(u8, base, other));
    }
}

test "output file names reject unsafe editor names, species and coordinates" {
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
    }) |name| {
        // Unsafe as an editor name.
        try std.testing.expectError(
            error.InvalidOutputFileName,
            buildOutputFileName(std.testing.allocator, 45.3, -75.7, .soil_or_eco, 2001, name),
        );
        // And equally unsafe as a species name, since both become path components.
        try std.testing.expectError(
            error.InvalidOutputFileName,
            buildOutputFileName(std.testing.allocator, 45.3, -75.7, .{ .species = name }, 2001, "hourly"),
        );
    }
    // Coordinates must be real, which also catches a caller still passing an index
    // for longitude only if it is out of range; the year bounds are unchanged.
    inline for (.{
        .{ 91.0, 0.0 },
        .{ -91.0, 0.0 },
        .{ 0.0, 181.0 },
        .{ 0.0, -181.0 },
    }) |pair| try std.testing.expectError(
        error.InvalidOutputFileName,
        buildOutputFileName(std.testing.allocator, pair[0], pair[1], .soil_or_eco, 2001, "hourly"),
    );
    try std.testing.expectError(
        error.InvalidOutputFileName,
        buildOutputFileName(std.testing.allocator, 45.3, -75.7, .soil_or_eco, 0, "hourly"),
    );
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
                .longitude_degrees_east = -75.7,
                .latitude_degrees_north = 45.3,
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

