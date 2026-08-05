//! Parser for the legacy-to-modern column map data file.
//!
//! The map is a data file, not code, because the mapping is evidence: each entry
//! cites the `foutp.f`/`fouts.f` heading slot and the `outsh.f`/`outsd.f`
//! expression that produced the legacy value, and carries the explicit unit
//! conversion factor plus the justification for its tolerance. A hardcoded table
//! could not be reviewed against the Fortran without recompiling.
//!
//! Format, one record per line, `|`-separated, `#` starts a comment:
//!
//! ```text
//! stream <id> <cadence> <legacy_relative_path> <modern_relative_path>
//! column <legacy_heading> | <modern_heading> | <factor> | <offset> | <abs> | <rel> | <justification>
//! ```
//!
//! A `column` line before any `stream` line is malformed. A `stream` with no
//! `column` lines is malformed, because it would compare nothing while
//! appearing to succeed.

const std = @import("std");
const comparison = @import("legacy_comparison.zig");

pub const Error = comparison.Error;

pub const ColumnMapping = struct {
    legacy_heading: []const u8,
    modern_heading: []const u8,
    conversion: comparison.Conversion,
    tolerance: comparison.Tolerance,
};

pub const StreamMapping = struct {
    id: []const u8,
    cadence: comparison.Cadence,
    legacy_path: []const u8,
    modern_path: []const u8,
    columns: []ColumnMapping,
};

pub const Map = struct {
    allocator: std.mem.Allocator,
    /// Owns every string referenced by the streams and columns.
    text: []u8,
    streams: []StreamMapping,

    pub fn deinit(self: *Map) void {
        for (self.streams) |stream| self.allocator.free(stream.columns);
        self.allocator.free(self.streams);
        self.allocator.free(self.text);
        self.* = undefined;
    }

    pub fn columnCount(self: Map) usize {
        var total: usize = 0;
        for (self.streams) |stream| total += stream.columns.len;
        return total;
    }
};

fn nextField(iterator: *std.mem.SplitIterator(u8, .scalar)) ![]const u8 {
    const field = iterator.next() orelse return Error.MalformedDataFile;
    return std.mem.trim(u8, field, " \t\r");
}

/// Parse `source`. The returned map borrows from an owned copy of `source`, so
/// the caller may free their own buffer immediately.
pub fn parse(allocator: std.mem.Allocator, source: []const u8) !Map {
    const owned = try allocator.dupe(u8, source);
    errdefer allocator.free(owned);

    var streams: std.ArrayList(StreamMapping) = .empty;
    errdefer {
        for (streams.items) |stream| allocator.free(stream.columns);
        streams.deinit(allocator);
    }
    var columns: std.ArrayList(ColumnMapping) = .empty;
    errdefer columns.deinit(allocator);
    var have_stream = false;

    var lines = std.mem.tokenizeAny(u8, owned, "\n");
    while (lines.next()) |raw_line| {
        const line = std.mem.trim(u8, raw_line, " \t\r");
        if (line.len == 0 or line[0] == '#') continue;

        if (std.mem.startsWith(u8, line, "stream ")) {
            if (have_stream) {
                if (columns.items.len == 0) return Error.MalformedDataFile;
                streams.items[streams.items.len - 1].columns =
                    try columns.toOwnedSlice(allocator);
            }
            var fields = std.mem.tokenizeAny(u8, line["stream ".len..], " \t");
            const id = fields.next() orelse return Error.MalformedDataFile;
            const cadence_text = fields.next() orelse return Error.MalformedDataFile;
            const cadence: comparison.Cadence =
                if (std.mem.eql(u8, cadence_text, "hourly"))
                    .hourly
                else if (std.mem.eql(u8, cadence_text, "daily"))
                    .daily
                else
                    return Error.MalformedDataFile;
            const legacy_path = fields.next() orelse return Error.MalformedDataFile;
            const modern_path = fields.next() orelse return Error.MalformedDataFile;
            if (fields.next() != null) return Error.MalformedDataFile;
            try streams.append(allocator, .{
                .id = id,
                .cadence = cadence,
                .legacy_path = legacy_path,
                .modern_path = modern_path,
                .columns = &.{},
            });
            have_stream = true;
            continue;
        }

        if (std.mem.startsWith(u8, line, "column ")) {
            if (!have_stream) return Error.MalformedDataFile;
            var fields = std.mem.splitScalar(u8, line["column ".len..], '|');
            const legacy_heading = try nextField(&fields);
            const modern_heading = try nextField(&fields);
            const factor_text = try nextField(&fields);
            const offset_text = try nextField(&fields);
            const absolute_text = try nextField(&fields);
            const relative_text = try nextField(&fields);
            const justification = try nextField(&fields);
            if (fields.next() != null) return Error.MalformedDataFile;
            if (legacy_heading.len == 0 or modern_heading.len == 0)
                return Error.MalformedDataFile;
            const mapping: ColumnMapping = .{
                .legacy_heading = legacy_heading,
                .modern_heading = modern_heading,
                .conversion = .{
                    .factor = std.fmt.parseFloat(f64, factor_text) catch
                        return Error.MalformedDataFile,
                    .offset = std.fmt.parseFloat(f64, offset_text) catch
                        return Error.MalformedDataFile,
                },
                .tolerance = .{
                    .absolute = std.fmt.parseFloat(f64, absolute_text) catch
                        return Error.MalformedDataFile,
                    .relative = std.fmt.parseFloat(f64, relative_text) catch
                        return Error.MalformedDataFile,
                    .justification = justification,
                },
            };
            // Reject a tolerance that cannot be justified at parse time rather
            // than discovering it mid-comparison.
            try mapping.tolerance.validate();
            try columns.append(allocator, mapping);
            continue;
        }
        return Error.MalformedDataFile;
    }

    if (!have_stream) return Error.MalformedDataFile;
    if (columns.items.len == 0) return Error.MalformedDataFile;
    streams.items[streams.items.len - 1].columns =
        try columns.toOwnedSlice(allocator);

    return .{
        .allocator = allocator,
        .text = owned,
        .streams = try streams.toOwnedSlice(allocator),
    };
}

// ---------------------------------------------------------------------- tests

const example =
    "# comment line\n" ++
    "stream soil_hourly_carbon hourly ottawa_run/01998f25ch1 010101998f25ch1.txt\n" ++
    "column SOIL_CO2_FLUX | carbon_dioxide_emission[umol m-2 s-1] | 1.0 | 0.0 | 1e-6 | 1e-6 | seven significant digits in E16.7E3\n" ++
    "column O2_1 | dissolved_oxygen_concentration_layer_1[g O2 m-3 water] | 1.0 | 0.0 | 1e-6 | 1e-6 | same units, no conversion\n" ++
    "\n" ++
    "stream soil_daily_carbon daily ottawa_run/01998f25cd1 010101998f25cd1.txt\n" ++
    "column HUMUS_C | humus_carbon[g C m-2] | 1.0 | 0.0 | 1e-5 | 1e-7 | large pool, absolute bound scaled to magnitude\n";

test "column map parses streams, cadence, conversions and tolerances" {
    var map = try parse(std.testing.allocator, example);
    defer map.deinit();
    try std.testing.expectEqual(@as(usize, 2), map.streams.len);
    try std.testing.expectEqualStrings("soil_hourly_carbon", map.streams[0].id);
    try std.testing.expectEqual(comparison.Cadence.hourly, map.streams[0].cadence);
    try std.testing.expectEqual(@as(usize, 2), map.streams[0].columns.len);
    try std.testing.expectEqual(comparison.Cadence.daily, map.streams[1].cadence);
    try std.testing.expectEqual(@as(usize, 1), map.streams[1].columns.len);
    try std.testing.expectEqual(@as(usize, 3), map.columnCount());
    try std.testing.expectEqualStrings(
        "carbon_dioxide_emission[umol m-2 s-1]",
        map.streams[0].columns[0].modern_heading,
    );
}

test "a column line before any stream line is malformed" {
    try std.testing.expectError(Error.MalformedDataFile, parse(
        std.testing.allocator,
        "column A | b[m] | 1.0 | 0.0 | 1e-6 | 1e-6 | reason\n",
    ));
}

test "a stream with no columns is malformed, since it would compare nothing" {
    try std.testing.expectError(Error.MalformedDataFile, parse(
        std.testing.allocator,
        "stream empty hourly legacy modern.txt\n",
    ));
    try std.testing.expectError(Error.MalformedDataFile, parse(
        std.testing.allocator,
        "stream a hourly l m.txt\n" ++
            "column A | b[m] | 1.0 | 0.0 | 1e-6 | 1e-6 | reason\n" ++
            "stream b hourly l2 m2.txt\n",
    ));
}

test "an unjustified tolerance in the data file is rejected at parse time" {
    try std.testing.expectError(comparison.Error.UnjustifiedTolerance, parse(
        std.testing.allocator,
        "stream a hourly l m.txt\n" ++
            "column A | b[m] | 1.0 | 0.0 | 1e-6 | 1e-6 |    \n",
    ));
}

test "an unknown cadence, unknown directive or short column line is malformed" {
    try std.testing.expectError(Error.MalformedDataFile, parse(
        std.testing.allocator,
        "stream a weekly l m.txt\ncolumn A | b[m] | 1 | 0 | 1e-6 | 1e-6 | r\n",
    ));
    try std.testing.expectError(Error.MalformedDataFile, parse(
        std.testing.allocator,
        "banana a hourly l m.txt\n",
    ));
    try std.testing.expectError(Error.MalformedDataFile, parse(
        std.testing.allocator,
        "stream a hourly l m.txt\ncolumn A | b[m] | 1 | 0 | 1e-6\n",
    ));
}

test "an empty data file is malformed rather than a silent zero-column pass" {
    try std.testing.expectError(Error.MalformedDataFile, parse(
        std.testing.allocator,
        "# nothing but comments\n\n",
    ));
}
