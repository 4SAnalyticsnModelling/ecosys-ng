const std = @import("std");
const execution_calendar_date = @import("../../driver/execution_calendar_date.zig");

pub const Season = enum(u2) {
    winter = 0,
    spring = 1,
    summer = 2,
    autumn = 3,
};

/// Exact WTHR four-period selector:
/// winter 335..year-end plus 1..59, spring 60..151,
/// summer 152..243, autumn 244..334.
pub fn season(
    day_of_year: u16,
    execution_year: u16,
) !Season {
    _ = execution_calendar_date.fromDayOfYear(
        day_of_year,
        execution_year,
    ) catch return error.InvalidClimateForcingCalendarDate;
    if (day_of_year > 334 or day_of_year <= 59) return .winter;
    if (day_of_year <= 151) return .spring;
    if (day_of_year <= 243) return .summer;
    return .autumn;
}

pub fn index(day_of_year: u16, execution_year: u16) !usize {
    return @intFromEnum(try season(day_of_year, execution_year));
}

/// Selects one of the four compulsory runtime climate records while retaining
/// strict arity validation.
pub fn select(
    comptime T: type,
    records: []const T,
    day_of_year: u16,
    execution_year: u16,
) !*const T {
    if (records.len != 4)
        return error.ClimateForcingPeriodArityMismatch;
    return &records[try index(day_of_year, execution_year)];
}

test "WTHR seasonal boundaries are exact and winter wraps year end" {
    const cases = [_]struct { day: u16, expected: Season }{
        .{ .day = 1, .expected = .winter },
        .{ .day = 59, .expected = .winter },
        .{ .day = 60, .expected = .spring },
        .{ .day = 151, .expected = .spring },
        .{ .day = 152, .expected = .summer },
        .{ .day = 243, .expected = .summer },
        .{ .day = 244, .expected = .autumn },
        .{ .day = 334, .expected = .autumn },
        .{ .day = 335, .expected = .winter },
        .{ .day = 365, .expected = .winter },
    };
    for (cases) |case|
        try std.testing.expectEqual(
            case.expected,
            try season(case.day, 1999),
        );
    try std.testing.expectEqual(Season.winter, try season(366, 1900));
}

test "runtime record selection enforces exactly four periods" {
    const records = [_]u32{ 10, 20, 30, 40 };
    try std.testing.expectEqual(
        @as(u32, 10),
        (try select(u32, &records, 10, 1999)).*,
    );
    try std.testing.expectEqual(
        @as(u32, 20),
        (try select(u32, &records, 100, 1999)).*,
    );
    try std.testing.expectEqual(
        @as(u32, 30),
        (try select(u32, &records, 200, 1999)).*,
    );
    try std.testing.expectEqual(
        @as(u32, 40),
        (try select(u32, &records, 300, 1999)).*,
    );
    try std.testing.expectError(
        error.ClimateForcingPeriodArityMismatch,
        select(u32, records[0..3], 1, 1999),
    );
}

test "invalid dates fail immediately under DAY chronology" {
    try std.testing.expectError(
        error.InvalidClimateForcingCalendarDate,
        season(0, 1999),
    );
    try std.testing.expectError(
        error.InvalidClimateForcingCalendarDate,
        season(366, 1901),
    );
    try std.testing.expectError(
        error.InvalidClimateForcingCalendarDate,
        season(367, 1900),
    );
    try std.testing.expectError(
        error.InvalidClimateForcingCalendarDate,
        season(1, 0),
    );
}
