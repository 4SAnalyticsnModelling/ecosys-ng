const std = @import("std");
const execution_calendar_date = @import("execution_calendar_date.zig");

pub const Error = error{
    InvalidLegacyMonth,
    InvalidLegacyDay,
};

/// Legacy `BLOCKDATA001` payload copied from `validation/legacy_ottawa_gfortran_16_1/source/BLOCKDATA001.f`.
///
/// `ICOR(M)` is indexed by month with Fortran 1-based indexing and values
/// used only in legacy day-index reconstruction logic.
pub const icor: [12]i16 = .{ 1, -1, 0, 0, 1, 1, 2, 3, 3, 4, 4, 5 };

/// Legacy calendar reconstruction used by translated readers/output selection paths.
///
/// Inputs use 1-based month/day. This intentionally preserves the historical
/// arithmetic shape rather than re-deriving a cleaner equivalent in
/// `execution_calendar_date.dayOfYear`.
pub fn dayIndex(month: u8, day: u8, year: u16) !u16 {
    if (month == 0 or month > 12) return Error.InvalidLegacyMonth;
    const validation_year = if (year == 0) 2000 else year;
    _ = execution_calendar_date.dayOfYear(.{
        .day = day,
        .month = month,
        .year = validation_year,
    }) catch return Error.InvalidLegacyDay;

    if (month == 1) return day;

    const leap: i16 = if (execution_calendar_date.isLeapYear(year) and month > 2) 1 else 0;
    const result: i16 = 30 * (@as(i16, month) - 1) +
        icor[month - 1] + @as(i16, @intCast(day)) + leap;
    if (result <= 0) return Error.InvalidLegacyDay;
    return @intCast(result);
}

test "legacy ICOR payload is byte-for-byte source data" {
    try std.testing.expectEqualSlices(i16, &.{ 1, -1, 0, 0, 1, 1, 2, 3, 3, 4, 4, 5 }, &icor);
}

test "legacy day index preserves source month-to-index arithmetic" {
    try std.testing.expectEqual(@as(u16, 1), try dayIndex(1, 1, 1999));
    try std.testing.expectEqual(@as(u16, 30), try dayIndex(2, 1, 1999));
    try std.testing.expectEqual(@as(u16, 57), try dayIndex(2, 28, 1999));
    try std.testing.expectEqual(@as(u16, 61), try dayIndex(3, 1, 1999));
    try std.testing.expectEqual(@as(u16, 365), try dayIndex(12, 30, 1999));
    try std.testing.expectEqual(@as(u16, 62), try dayIndex(3, 1, 2000));
}

test "legacy month/day conversion rejects invalid month and zero day" {
    try std.testing.expectError(Error.InvalidLegacyMonth, dayIndex(13, 1, 1999));
    try std.testing.expectError(Error.InvalidLegacyMonth, dayIndex(0, 1, 1999));
    try std.testing.expectError(Error.InvalidLegacyDay, dayIndex(1, 0, 1999));
    try std.testing.expectError(Error.InvalidLegacyDay, dayIndex(2, 29, 1999));
    try std.testing.expectEqual(@as(u16, 58), try dayIndex(2, 29, 1900));
}
