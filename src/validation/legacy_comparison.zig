//! Divergence-reporting core for legacy-versus-ecosys-ng output comparison.
//!
//! This module is deliberately free of any front-end concern: it parses the two
//! table dialects, aligns rows on an absolute simulated-hour ordinal, and
//! reports per-column divergence including the **earliest** divergent
//! timestamp, which is the field that localizes a bug. `tools/compare_legacy.zig`
//! is one front end; lane A8's initial-state inventory comparison is intended to
//! be the second, reusing `Divergence` and `compareColumn` unchanged.
//!
//! Non-negotiable behaviour: a missing, truncated, header-only or
//! column-mismatched stream is an error, never a silent pass. A harness that
//! passes when it measured nothing is worse than no harness.

const std = @import("std");

/// A simulated instant, as carried by both dialects.
///
/// `hour` is 1..24 in the legacy dialect and 1..23 plus a 0 that belongs to the
/// following calendar day in the modern dialect. `ordinal` erases that
/// difference; see `Timestamp.ordinal`.
pub const Timestamp = struct {
    year: u16,
    day_of_year: u16,
    /// 0 means "hour 24 of the previous calendar day" in the modern dialect and
    /// is normalized during parsing, so a stored 0 only occurs for daily rows.
    hour: u8,

    /// Absolute hour ordinal within the simulated year.
    ///
    /// Legacy day 1 hour 24 and modern day 2 hour 0 are the same instant, and
    /// both map to 1*24+24 == 2*24 == 48. Daily rows carry hour 0 with the
    /// calendar day they summarize, giving day_of_year*24.
    pub fn ordinal(self: Timestamp) i64 {
        return @as(i64, self.year) * 24 * 400 +
            @as(i64, self.day_of_year) * 24 + @as(i64, self.hour);
    }

    pub fn format(self: Timestamp, writer: *std.Io.Writer) std.Io.Writer.Error!void {
        try writer.print(
            "year{d}-day{d}-hour{d}",
            .{ self.year, self.day_of_year, self.hour },
        );
    }
};

/// Whether a stream carries an hour column.
pub const Cadence = enum { hourly, daily };

pub const Error = error{
    LegacyStreamMissing,
    LegacyStreamEmpty,
    LegacyHeadingUnparsable,
    LegacyRowTruncated,
    LegacyValueUnparsable,
    ModernStreamMissing,
    ModernStreamEmpty,
    ModernHeadingMissing,
    ModernRowTruncated,
    ModernValueUnparsable,
    LegacyColumnNotFound,
    ModernColumnNotFound,
    NoAlignedRows,
    DuplicateOrdinal,
    UnjustifiedTolerance,
    MalformedDataFile,
};

/// A parsed output table: owned headings plus a dense row-major value block.
pub const Table = struct {
    allocator: std.mem.Allocator,
    /// Backing store for every heading string, freed as one block.
    heading_bytes: []u8,
    headings: [][]const u8,
    timestamps: []Timestamp,
    values: []f64,
    row_count: usize,
    column_count: usize,

    pub fn deinit(self: *Table) void {
        self.allocator.free(self.heading_bytes);
        self.allocator.free(self.headings);
        self.allocator.free(self.timestamps);
        self.allocator.free(self.values);
        self.* = undefined;
    }

    pub fn columnIndex(self: Table, name: []const u8) ?usize {
        for (self.headings, 0..) |heading, index| {
            if (std.mem.eql(u8, heading, name)) return index;
        }
        return null;
    }

    pub fn at(self: Table, row: usize, column: usize) f64 {
        return self.values[row * self.column_count + column];
    }
};

/// Days elapsed before the first of `month`, using the source's modulo-four
/// leap chronology (`fouts.f` retains exactly that rule, so the comparison
/// harness must too rather than inventing a Gregorian correction).
fn dayOfYear(year: u16, month: u8, day: u8) !u16 {
    if (month < 1 or month > 12) return Error.LegacyRowTruncated;
    const leap = (year % 4) == 0;
    const cumulative = [12]u16{ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
    var result: u16 = cumulative[month - 1] + day;
    if (leap and month > 2) result += 1;
    return result;
}

/// Parse a legacy fixed-width table as written by `outsh.f`/`outsd.f`.
///
/// Heading line is Fortran `(A12,2A8,50A16)` for hourly and `(A12,A8,50A16)`
/// for daily, so the scientific headings begin at byte 28 and 20 respectively
/// and occupy 16 bytes each. Data rows are
/// `(A16,F8.3,4X,A8,I8,50E16.7E3)` / `(A16,F8.3,4X,A8,50E16.7E3)`; the fields
/// are wide enough that whitespace splitting is unambiguous, and splitting is
/// robust against the gfortran column drift seen for `E16.7E3` overflow.
pub fn parseLegacy(
    allocator: std.mem.Allocator,
    text: []const u8,
    cadence: Cadence,
) !Table {
    if (text.len == 0) return Error.LegacyStreamEmpty;
    var lines = std.mem.tokenizeAny(u8, text, "\r\n");
    const heading_line = lines.next() orelse return Error.LegacyStreamEmpty;
    const heading_start: usize = switch (cadence) {
        .hourly => 28,
        .daily => 20,
    };
    if (heading_line.len <= heading_start) return Error.LegacyHeadingUnparsable;

    var heading_bytes: std.ArrayList(u8) = .empty;
    errdefer heading_bytes.deinit(allocator);
    var spans: std.ArrayList([2]usize) = .empty;
    defer spans.deinit(allocator);

    var cursor = heading_start;
    while (cursor < heading_line.len) : (cursor += 16) {
        const end = @min(cursor + 16, heading_line.len);
        const name = std.mem.trim(u8, heading_line[cursor..end], " \t");
        if (name.len == 0) continue;
        const start = heading_bytes.items.len;
        try heading_bytes.appendSlice(allocator, name);
        try spans.append(allocator, .{ start, heading_bytes.items.len });
    }
    if (spans.items.len == 0) return Error.LegacyHeadingUnparsable;
    const column_count = spans.items.len;

    // Leading tokens per row: tag, DOY, DATE, and HOUR for hourly streams.
    const leading_tokens: usize = switch (cadence) {
        .hourly => 4,
        .daily => 3,
    };

    var timestamps: std.ArrayList(Timestamp) = .empty;
    errdefer timestamps.deinit(allocator);
    var values: std.ArrayList(f64) = .empty;
    errdefer values.deinit(allocator);

    while (lines.next()) |line| {
        if (std.mem.trim(u8, line, " \t").len == 0) continue;
        var tokens = std.mem.tokenizeAny(u8, line, " \t");
        _ = tokens.next() orelse return Error.LegacyRowTruncated; // cell tag
        _ = tokens.next() orelse return Error.LegacyRowTruncated; // fractional DOY
        const date = tokens.next() orelse return Error.LegacyRowTruncated;
        if (date.len != 8) return Error.LegacyRowTruncated;
        const day = std.fmt.parseInt(u8, date[0..2], 10) catch
            return Error.LegacyRowTruncated;
        const month = std.fmt.parseInt(u8, date[2..4], 10) catch
            return Error.LegacyRowTruncated;
        const year = std.fmt.parseInt(u16, date[4..8], 10) catch
            return Error.LegacyRowTruncated;
        var hour: u8 = 0;
        if (cadence == .hourly) {
            const hour_text = tokens.next() orelse return Error.LegacyRowTruncated;
            hour = std.fmt.parseInt(u8, hour_text, 10) catch
                return Error.LegacyRowTruncated;
        }
        std.debug.assert(leading_tokens >= 3);
        try timestamps.append(allocator, .{
            .year = year,
            .day_of_year = try dayOfYear(year, month, day),
            .hour = hour,
        });
        var seen: usize = 0;
        while (seen < column_count) : (seen += 1) {
            const token = tokens.next() orelse return Error.LegacyRowTruncated;
            const parsed = std.fmt.parseFloat(f64, token) catch
                return Error.LegacyValueUnparsable;
            try values.append(allocator, parsed);
        }
    }
    if (timestamps.items.len == 0) return Error.LegacyStreamEmpty;

    const owned_bytes = try heading_bytes.toOwnedSlice(allocator);
    errdefer allocator.free(owned_bytes);
    const headings = try allocator.alloc([]const u8, column_count);
    for (spans.items, 0..) |span, index| {
        headings[index] = owned_bytes[span[0]..span[1]];
    }
    const row_count = timestamps.items.len;
    const owned_timestamps = try timestamps.toOwnedSlice(allocator);
    errdefer allocator.free(owned_timestamps);
    return .{
        .allocator = allocator,
        .heading_bytes = owned_bytes,
        .headings = headings,
        .timestamps = owned_timestamps,
        .values = try values.toOwnedSlice(allocator),
        .row_count = row_count,
        .column_count = column_count,
    };
}

/// Parse a modern tab-delimited stream as written by `stream.zig`.
///
/// The leading index columns are `year day_of_year month day hour lon lat` for
/// hourly streams and the same without `hour` for daily streams; every
/// remaining heading is a scientific `name[unit]` column enforced by
/// `unit_convention.zig`. Modern hour 0 belongs to the previous
/// calendar day, so it is normalized to hour 24 of `day_of_year - 1` here and
/// nowhere else.
pub fn parseModern(allocator: std.mem.Allocator, text: []const u8) !Table {
    if (text.len == 0) return Error.ModernStreamEmpty;
    var lines = std.mem.tokenizeAny(u8, text, "\r\n");
    const heading_line = lines.next() orelse return Error.ModernStreamEmpty;

    var all_headings: std.ArrayList([]const u8) = .empty;
    defer all_headings.deinit(allocator);
    var heading_fields = std.mem.splitScalar(u8, heading_line, '\t');
    while (heading_fields.next()) |field| {
        const name = std.mem.trim(u8, field, " \t\r");
        if (name.len == 0) continue;
        try all_headings.append(allocator, name);
    }
    if (all_headings.items.len == 0) return Error.ModernHeadingMissing;

    var year_column: ?usize = null;
    var day_column: ?usize = null;
    var hour_column: ?usize = null;
    var first_scientific: usize = 0;
    for (all_headings.items, 0..) |heading, index| {
        if (std.mem.eql(u8, heading, "year")) year_column = index;
        if (std.mem.eql(u8, heading, "day_of_year")) day_column = index;
        if (std.mem.eql(u8, heading, "hour")) hour_column = index;
        if (std.mem.indexOfScalar(u8, heading, '[') != null) {
            first_scientific = index;
            break;
        }
    }
    if (year_column == null or day_column == null) return Error.ModernHeadingMissing;
    if (first_scientific == 0) return Error.ModernHeadingMissing;
    const column_count = all_headings.items.len - first_scientific;

    var heading_bytes: std.ArrayList(u8) = .empty;
    errdefer heading_bytes.deinit(allocator);
    var spans: std.ArrayList([2]usize) = .empty;
    defer spans.deinit(allocator);
    for (all_headings.items[first_scientific..]) |heading| {
        const start = heading_bytes.items.len;
        try heading_bytes.appendSlice(allocator, heading);
        try spans.append(allocator, .{ start, heading_bytes.items.len });
    }

    var timestamps: std.ArrayList(Timestamp) = .empty;
    errdefer timestamps.deinit(allocator);
    var values: std.ArrayList(f64) = .empty;
    errdefer values.deinit(allocator);

    while (lines.next()) |line| {
        if (std.mem.trim(u8, line, " \t").len == 0) continue;
        var fields = std.mem.splitScalar(u8, line, '\t');
        var row: std.ArrayList([]const u8) = .empty;
        defer row.deinit(allocator);
        while (fields.next()) |field| {
            try row.append(allocator, std.mem.trim(u8, field, " \t\r"));
        }
        if (row.items.len < all_headings.items.len) return Error.ModernRowTruncated;
        const year = std.fmt.parseInt(u16, row.items[year_column.?], 10) catch
            return Error.ModernRowTruncated;
        var day_of_year = std.fmt.parseInt(u16, row.items[day_column.?], 10) catch
            return Error.ModernRowTruncated;
        var hour: u8 = 0;
        if (hour_column) |index| {
            hour = std.fmt.parseInt(u8, row.items[index], 10) catch
                return Error.ModernRowTruncated;
            if (hour == 0) {
                if (day_of_year == 0) return Error.ModernRowTruncated;
                day_of_year -= 1;
                hour = 24;
            }
        }
        try timestamps.append(allocator, .{
            .year = year,
            .day_of_year = day_of_year,
            .hour = hour,
        });
        for (row.items[first_scientific..][0..column_count]) |field| {
            const parsed = std.fmt.parseFloat(f64, field) catch
                return Error.ModernValueUnparsable;
            try values.append(allocator, parsed);
        }
    }
    if (timestamps.items.len == 0) return Error.ModernStreamEmpty;

    const owned_bytes = try heading_bytes.toOwnedSlice(allocator);
    errdefer allocator.free(owned_bytes);
    const headings = try allocator.alloc([]const u8, column_count);
    for (spans.items, 0..) |span, index| {
        headings[index] = owned_bytes[span[0]..span[1]];
    }
    const row_count = timestamps.items.len;
    const owned_timestamps = try timestamps.toOwnedSlice(allocator);
    errdefer allocator.free(owned_timestamps);
    return .{
        .allocator = allocator,
        .heading_bytes = owned_bytes,
        .headings = headings,
        .timestamps = owned_timestamps,
        .values = try values.toOwnedSlice(allocator),
        .row_count = row_count,
        .column_count = column_count,
    };
}

/// Explicit unit reconciliation. Never compare raw numbers across differing
/// units; the factor and offset are recorded per column in the column-map data
/// file, together with the reason they are what they are.
pub const Conversion = struct {
    factor: f64 = 1.0,
    offset: f64 = 0.0,

    pub fn apply(self: Conversion, legacy_value: f64) f64 {
        return legacy_value * self.factor + self.offset;
    }
};

/// Per-column tolerance, with a mandatory justification. A tolerance without a
/// justification is a hidden failure, so `Tolerance.validate` rejects it.
pub const Tolerance = struct {
    absolute: f64,
    relative: f64,
    justification: []const u8,

    pub fn validate(self: Tolerance) !void {
        if (std.mem.trim(u8, self.justification, " \t").len == 0)
            return Error.UnjustifiedTolerance;
        if (!std.math.isFinite(self.absolute) or self.absolute < 0)
            return Error.UnjustifiedTolerance;
        if (!std.math.isFinite(self.relative) or self.relative < 0)
            return Error.UnjustifiedTolerance;
    }

    /// A divergence counts only when it breaks both bounds, so that a tiny
    /// absolute difference on a large value and a large relative difference on
    /// a near-zero value are both accepted deliberately rather than by luck.
    pub fn exceeded(self: Tolerance, absolute: f64, relative: f64) bool {
        return absolute > self.absolute and relative > self.relative;
    }
};

/// The verdict for one aligned column.
pub const Divergence = struct {
    compared_rows: usize = 0,
    max_absolute: f64 = 0.0,
    max_relative: f64 = 0.0,
    max_absolute_at: ?Timestamp = null,
    /// Signed `modern - legacy` at the row where `max_absolute` occurred.
    /// Magnitude alone cannot tell a lane whether the modern model runs hot or
    /// cold, which is the first thing an attribution needs, so the sign is
    /// carried explicitly rather than reconstructed.
    max_absolute_signed: f64 = 0.0,
    earliest_divergent: ?Timestamp = null,
    earliest_divergent_absolute: f64 = 0.0,
    earliest_divergent_relative: f64 = 0.0,
    /// Signed `modern - legacy` at the earliest divergent row.
    earliest_divergent_signed: f64 = 0.0,
    divergent_rows: usize = 0,
    nonfinite_rows: usize = 0,

    pub fn within(self: Divergence) bool {
        return self.divergent_rows == 0 and self.nonfinite_rows == 0;
    }

    /// `-1`, `0` or `+1` for the signed divergence at the maximum. A report
    /// reader wants the direction, not the magnitude, at a glance.
    pub fn signCharacter(self: Divergence) u8 {
        if (self.max_absolute_signed > 0) return '+';
        if (self.max_absolute_signed < 0) return '-';
        return '0';
    }
};

fn relativeDifference(expected: f64, actual: f64) f64 {
    const scale = @max(@abs(expected), @abs(actual));
    if (scale == 0.0) return 0.0;
    return @abs(expected - actual) / scale;
}

/// Row alignment on the absolute hour ordinal.
pub const Alignment = struct {
    /// Ordinals present in both tables, ascending.
    legacy_rows: []usize,
    modern_rows: []usize,
    legacy_only: usize,
    modern_only: usize,
    allocator: std.mem.Allocator,

    pub fn deinit(self: *Alignment) void {
        self.allocator.free(self.legacy_rows);
        self.allocator.free(self.modern_rows);
        self.* = undefined;
    }

    pub fn count(self: Alignment) usize {
        return self.legacy_rows.len;
    }
};

pub fn alignRows(
    allocator: std.mem.Allocator,
    legacy: Table,
    modern: Table,
) !Alignment {
    var legacy_index: std.AutoHashMapUnmanaged(i64, usize) = .empty;
    defer legacy_index.deinit(allocator);
    for (legacy.timestamps, 0..) |timestamp, row| {
        const entry = try legacy_index.getOrPut(allocator, timestamp.ordinal());
        if (entry.found_existing) return Error.DuplicateOrdinal;
        entry.value_ptr.* = row;
    }
    var legacy_rows: std.ArrayList(usize) = .empty;
    errdefer legacy_rows.deinit(allocator);
    var modern_rows: std.ArrayList(usize) = .empty;
    errdefer modern_rows.deinit(allocator);
    var matched: std.AutoHashMapUnmanaged(i64, void) = .empty;
    defer matched.deinit(allocator);
    var modern_only: usize = 0;
    for (modern.timestamps, 0..) |timestamp, row| {
        const ordinal = timestamp.ordinal();
        if (legacy_index.get(ordinal)) |legacy_row| {
            if ((try matched.getOrPut(allocator, ordinal)).found_existing)
                return Error.DuplicateOrdinal;
            try legacy_rows.append(allocator, legacy_row);
            try modern_rows.append(allocator, row);
        } else modern_only += 1;
    }
    if (legacy_rows.items.len == 0) return Error.NoAlignedRows;
    return .{
        .allocator = allocator,
        .legacy_only = legacy.row_count - legacy_rows.items.len,
        .modern_only = modern_only,
        .legacy_rows = try legacy_rows.toOwnedSlice(allocator),
        .modern_rows = try modern_rows.toOwnedSlice(allocator),
    };
}

/// Compare one legacy column against one modern column over aligned rows.
///
/// Atomic in the sense that matters here: it either produces a complete verdict
/// for the column or returns an error, and it never mutates either table.
pub fn compareColumn(
    legacy: Table,
    legacy_column: usize,
    modern: Table,
    modern_column: usize,
    conversion: Conversion,
    tolerance: Tolerance,
    alignment: Alignment,
) !Divergence {
    try tolerance.validate();
    if (legacy_column >= legacy.column_count) return Error.LegacyColumnNotFound;
    if (modern_column >= modern.column_count) return Error.ModernColumnNotFound;
    var result: Divergence = .{};
    for (alignment.legacy_rows, alignment.modern_rows) |legacy_row, modern_row| {
        const expected = conversion.apply(legacy.at(legacy_row, legacy_column));
        const actual = modern.at(modern_row, modern_column);
        const timestamp = modern.timestamps[modern_row];
        if (!std.math.isFinite(expected) or !std.math.isFinite(actual)) {
            result.nonfinite_rows += 1;
            if (result.earliest_divergent == null) result.earliest_divergent = timestamp;
            continue;
        }
        result.compared_rows += 1;
        const absolute = @abs(expected - actual);
        const relative = relativeDifference(expected, actual);
        if (absolute > result.max_absolute) {
            result.max_absolute = absolute;
            result.max_absolute_at = timestamp;
            result.max_absolute_signed = actual - expected;
        }
        result.max_relative = @max(result.max_relative, relative);
        if (tolerance.exceeded(absolute, relative)) {
            result.divergent_rows += 1;
            if (result.earliest_divergent == null) {
                result.earliest_divergent = timestamp;
                result.earliest_divergent_absolute = absolute;
                result.earliest_divergent_relative = relative;
                result.earliest_divergent_signed = actual - expected;
            }
        }
    }
    return result;
}

// ---------------------------------------------------------------------- tests

test "legacy hourly heading offset and value parse match the outsh.f format" {
    const text =
        "    DOY     DATE    HOUR    SOIL_CO2_FLUX   ECO_CO2_FLUX    \n" ++
        "010101998f25ch1    0.042    01011998       1 -0.8703152E+001 -0.8703152E+001\n" ++
        "010101998f25ch1    0.083    01011998       2 -0.6444155E+001 -0.6444155E+001\n";
    var table = try parseLegacy(std.testing.allocator, text, .hourly);
    defer table.deinit();
    try std.testing.expectEqual(@as(usize, 2), table.column_count);
    try std.testing.expectEqualStrings("SOIL_CO2_FLUX", table.headings[0]);
    try std.testing.expectEqualStrings("ECO_CO2_FLUX", table.headings[1]);
    try std.testing.expectEqual(@as(usize, 2), table.row_count);
    try std.testing.expectApproxEqAbs(-8.703152, table.at(0, 0), 1e-12);
    try std.testing.expectEqual(@as(u16, 1), table.timestamps[0].day_of_year);
    try std.testing.expectEqual(@as(u8, 1), table.timestamps[0].hour);
}

test "legacy daily heading offset is 20 bytes per the outsd.f format" {
    const text =
        "    DOY     DATE    LIT+POM_C       HUMUS_C     \n" ++
        "010101998f25cd1    1.000    01011998  0.1604418E+004  0.4541760E+004\n";
    var table = try parseLegacy(std.testing.allocator, text, .daily);
    defer table.deinit();
    try std.testing.expectEqualStrings("LIT+POM_C", table.headings[0]);
    try std.testing.expectEqualStrings("HUMUS_C", table.headings[1]);
    try std.testing.expectApproxEqAbs(1604.418, table.at(0, 0), 1e-9);
    try std.testing.expectEqual(@as(u8, 0), table.timestamps[0].hour);
}

test "legacy leap chronology follows the source modulo-four rule" {
    const text =
        "    DOY     DATE    HOUR    X               \n" ++
        "tag    60.000    01032000       1  0.1000000E+001\n";
    var table = try parseLegacy(std.testing.allocator, text, .hourly);
    defer table.deinit();
    // 2000 is a modulo-four leap year in the source chronology, so 1 March is
    // day 61, not day 60.
    try std.testing.expectEqual(@as(u16, 61), table.timestamps[0].day_of_year);
}

test "modern hour zero normalizes onto the previous calendar day" {
    const text =
        "year\tday_of_year\tmonth\tday\thour\tlon\tlat\tcarbon_dioxide_emission[umol m-2 s-1]\n" ++
        "1998\t1\t1\t1\t23\t1\t1\t1.5e0\n" ++
        "1998\t2\t1\t2\t0\t1\t1\t2.5e0\n";
    var table = try parseModern(std.testing.allocator, text);
    defer table.deinit();
    try std.testing.expectEqual(@as(usize, 1), table.column_count);
    try std.testing.expectEqualStrings(
        "carbon_dioxide_emission[umol m-2 s-1]",
        table.headings[0],
    );
    try std.testing.expectEqual(@as(u16, 1), table.timestamps[1].day_of_year);
    try std.testing.expectEqual(@as(u8, 24), table.timestamps[1].hour);
}

test "legacy hour twenty-four and modern hour zero share one ordinal" {
    const legacy: Timestamp = .{ .year = 1998, .day_of_year = 1, .hour = 24 };
    const modern: Timestamp = .{ .year = 1998, .day_of_year = 1, .hour = 24 };
    try std.testing.expectEqual(legacy.ordinal(), modern.ordinal());
    const next: Timestamp = .{ .year = 1998, .day_of_year = 2, .hour = 1 };
    try std.testing.expect(next.ordinal() > legacy.ordinal());
}

test "empty and header-only streams are errors, never a silent pass" {
    try std.testing.expectError(
        Error.LegacyStreamEmpty,
        parseLegacy(std.testing.allocator, "", .hourly),
    );
    try std.testing.expectError(
        Error.LegacyStreamEmpty,
        parseLegacy(
            std.testing.allocator,
            "    DOY     DATE    HOUR    X               \n",
            .hourly,
        ),
    );
    try std.testing.expectError(
        Error.ModernStreamEmpty,
        parseModern(std.testing.allocator, ""),
    );
    try std.testing.expectError(
        Error.ModernStreamEmpty,
        parseModern(
            std.testing.allocator,
            "year\tday_of_year\thour\tlon\tlat\tx[m]\n",
        ),
    );
}

test "a truncated legacy row is rejected rather than padded" {
    const text =
        "    DOY     DATE    HOUR    A               B               \n" ++
        "tag    0.042    01011998       1  0.1000000E+001\n";
    try std.testing.expectError(
        Error.LegacyRowTruncated,
        parseLegacy(std.testing.allocator, text, .hourly),
    );
}

test "a tolerance without justification is rejected" {
    const bad: Tolerance = .{ .absolute = 1e-9, .relative = 1e-9, .justification = "  " };
    try std.testing.expectError(Error.UnjustifiedTolerance, bad.validate());
    const good: Tolerance = .{
        .absolute = 1e-9,
        .relative = 1e-9,
        .justification = "gfortran E16.7E3 carries seven significant digits",
    };
    try good.validate();
}

test "divergence reports the earliest divergent timestamp, not just the largest" {
    const legacy_text =
        "    DOY     DATE    HOUR    X               \n" ++
        "tag    0.042    01011998       1  0.1000000E+001\n" ++
        "tag    0.083    01011998       2  0.1000000E+001\n" ++
        "tag    0.125    01011998       3  0.1000000E+001\n";
    const modern_text =
        "year\tday_of_year\tmonth\tday\thour\tlon\tlat\tx[m]\n" ++
        "1998\t1\t1\t1\t1\t1\t1\t1.0e0\n" ++
        "1998\t1\t1\t1\t2\t1\t1\t1.5e0\n" ++
        "1998\t1\t1\t1\t3\t1\t1\t9.0e0\n";
    var legacy = try parseLegacy(std.testing.allocator, legacy_text, .hourly);
    defer legacy.deinit();
    var modern = try parseModern(std.testing.allocator, modern_text);
    defer modern.deinit();
    var alignment = try alignRows(std.testing.allocator, legacy, modern);
    defer alignment.deinit();
    try std.testing.expectEqual(@as(usize, 3), alignment.count());
    const verdict = try compareColumn(
        legacy,
        0,
        modern,
        0,
        .{},
        .{ .absolute = 1e-6, .relative = 1e-6, .justification = "test" },
        alignment,
    );
    try std.testing.expectEqual(@as(usize, 2), verdict.divergent_rows);
    try std.testing.expectEqual(@as(u8, 2), verdict.earliest_divergent.?.hour);
    try std.testing.expectEqual(@as(u8, 3), verdict.max_absolute_at.?.hour);
    try std.testing.expectApproxEqAbs(@as(f64, 8.0), verdict.max_absolute, 1e-12);
    try std.testing.expect(!verdict.within());
}

test "unit conversion is applied to the legacy side and agreement is then exact" {
    // Legacy hourly water is millimetres via *1000.0 in outsh.f from metres.
    const legacy_text =
        "    DOY     DATE    HOUR    EVAPN           \n" ++
        "tag    0.042    01011998       1  0.2000000E+001\n";
    const modern_text =
        "year\tday_of_year\tmonth\tday\thour\tlon\tlat\tevapotranspiration[m]\n" ++
        "1998\t1\t1\t1\t1\t1\t1\t2.0e-3\n";
    var legacy = try parseLegacy(std.testing.allocator, legacy_text, .hourly);
    defer legacy.deinit();
    var modern = try parseModern(std.testing.allocator, modern_text);
    defer modern.deinit();
    var alignment = try alignRows(std.testing.allocator, legacy, modern);
    defer alignment.deinit();
    const raw = try compareColumn(
        legacy,
        0,
        modern,
        0,
        .{},
        .{ .absolute = 1e-9, .relative = 1e-9, .justification = "test" },
        alignment,
    );
    try std.testing.expect(!raw.within());
    const converted = try compareColumn(
        legacy,
        0,
        modern,
        0,
        .{ .factor = 1.0e-3 },
        .{ .absolute = 1e-12, .relative = 1e-12, .justification = "millimetre to metre" },
        alignment,
    );
    try std.testing.expect(converted.within());
}

test "non-overlapping timestamps are an error, not an empty pass" {
    const legacy_text =
        "    DOY     DATE    HOUR    X               \n" ++
        "tag    0.042    01011998       1  0.1000000E+001\n";
    const modern_text =
        "year\tday_of_year\tmonth\tday\thour\tlon\tlat\tx[m]\n" ++
        "1998\t200\t7\t19\t5\t1\t1\t1.0e0\n";
    var legacy = try parseLegacy(std.testing.allocator, legacy_text, .hourly);
    defer legacy.deinit();
    var modern = try parseModern(std.testing.allocator, modern_text);
    defer modern.deinit();
    try std.testing.expectError(
        Error.NoAlignedRows,
        alignRows(std.testing.allocator, legacy, modern),
    );
}

test "a partial modern stream still aligns and reports its shorter overlap" {
    // The Ottawa example currently stops in a coupled surface-litter phosphate
    // equilibrium tail, so the modern stream is a prefix of the legacy one. The
    // harness must produce a useful partial report for exactly this case.
    const legacy_text =
        "    DOY     DATE    HOUR    X               \n" ++
        "tag    0.042    01011998       1  0.1000000E+001\n" ++
        "tag    0.083    01011998       2  0.2000000E+001\n" ++
        "tag    0.125    01011998       3  0.3000000E+001\n";
    const modern_text =
        "year\tday_of_year\tmonth\tday\thour\tlon\tlat\tx[m]\n" ++
        "1998\t1\t1\t1\t1\t1\t1\t1.0e0\n" ++
        "1998\t1\t1\t1\t2\t1\t1\t2.0e0\n";
    var legacy = try parseLegacy(std.testing.allocator, legacy_text, .hourly);
    defer legacy.deinit();
    var modern = try parseModern(std.testing.allocator, modern_text);
    defer modern.deinit();
    var alignment = try alignRows(std.testing.allocator, legacy, modern);
    defer alignment.deinit();
    try std.testing.expectEqual(@as(usize, 2), alignment.count());
    try std.testing.expectEqual(@as(usize, 1), alignment.legacy_only);
    try std.testing.expectEqual(@as(usize, 0), alignment.modern_only);
    const verdict = try compareColumn(
        legacy,
        0,
        modern,
        0,
        .{},
        .{ .absolute = 1e-9, .relative = 1e-9, .justification = "test" },
        alignment,
    );
    try std.testing.expect(verdict.within());
    try std.testing.expectEqual(@as(usize, 2), verdict.compared_rows);
}

test "column lookup is by exact heading text including the unit bracket" {
    const modern_text =
        "year\tday_of_year\tmonth\tday\thour\tlon\tlat\tsoil_temperature_layer_1[degC]\n" ++
        "1998\t1\t1\t1\t1\t1\t1\t3.0e0\n";
    var modern = try parseModern(std.testing.allocator, modern_text);
    defer modern.deinit();
    try std.testing.expectEqual(
        @as(?usize, 0),
        modern.columnIndex("soil_temperature_layer_1[degC]"),
    );
    try std.testing.expectEqual(
        @as(?usize, null),
        modern.columnIndex("soil_temperature_layer_1"),
    );
}

test "a non-finite value counts as divergence and sets the earliest timestamp" {
    const legacy_text =
        "    DOY     DATE    HOUR    X               \n" ++
        "tag    0.042    01011998       1  0.1000000E+001\n" ++
        "tag    0.083    01011998       2  0.1000000E+001\n";
    const modern_text =
        "year\tday_of_year\tmonth\tday\thour\tlon\tlat\tx[m]\n" ++
        "1998\t1\t1\t1\t1\t1\t1\t1.0e0\n" ++
        "1998\t1\t1\t1\t2\t1\t1\tnan\n";
    var legacy = try parseLegacy(std.testing.allocator, legacy_text, .hourly);
    defer legacy.deinit();
    var modern = try parseModern(std.testing.allocator, modern_text);
    defer modern.deinit();
    var alignment = try alignRows(std.testing.allocator, legacy, modern);
    defer alignment.deinit();
    const verdict = try compareColumn(
        legacy,
        0,
        modern,
        0,
        .{},
        .{ .absolute = 1e9, .relative = 1e9, .justification = "wide on purpose" },
        alignment,
    );
    try std.testing.expectEqual(@as(usize, 1), verdict.nonfinite_rows);
    try std.testing.expect(!verdict.within());
    try std.testing.expectEqual(@as(u8, 2), verdict.earliest_divergent.?.hour);
}

test "divergence order does not depend on row order of the modern table" {
    const legacy_text =
        "    DOY     DATE    HOUR    X               \n" ++
        "tag    0.042    01011998       1  0.1000000E+001\n" ++
        "tag    0.083    01011998       2  0.2000000E+001\n" ++
        "tag    0.125    01011998       3  0.3000000E+001\n";
    const ascending =
        "year\tday_of_year\tmonth\tday\thour\tlon\tlat\tx[m]\n" ++
        "1998\t1\t1\t1\t1\t1\t1\t1.0e0\n" ++
        "1998\t1\t1\t1\t2\t1\t1\t2.5e0\n" ++
        "1998\t1\t1\t1\t3\t1\t1\t3.0e0\n";
    const shuffled =
        "year\tday_of_year\tmonth\tday\thour\tlon\tlat\tx[m]\n" ++
        "1998\t1\t1\t1\t3\t1\t1\t3.0e0\n" ++
        "1998\t1\t1\t1\t1\t1\t1\t1.0e0\n" ++
        "1998\t1\t1\t1\t2\t1\t1\t2.5e0\n";
    var legacy = try parseLegacy(std.testing.allocator, legacy_text, .hourly);
    defer legacy.deinit();
    const tolerance: Tolerance = .{
        .absolute = 1e-9,
        .relative = 1e-9,
        .justification = "test",
    };
    var first_max: f64 = 0;
    var first_hour: u8 = 0;
    for ([_][]const u8{ ascending, shuffled }, 0..) |text, index| {
        var modern = try parseModern(std.testing.allocator, text);
        defer modern.deinit();
        var alignment = try alignRows(std.testing.allocator, legacy, modern);
        defer alignment.deinit();
        const verdict = try compareColumn(
            legacy,
            0,
            modern,
            0,
            .{},
            tolerance,
            alignment,
        );
        if (index == 0) {
            first_max = verdict.max_absolute;
            first_hour = verdict.earliest_divergent.?.hour;
        } else {
            try std.testing.expectApproxEqAbs(first_max, verdict.max_absolute, 0);
            try std.testing.expectEqual(first_hour, verdict.earliest_divergent.?.hour);
        }
    }
}
