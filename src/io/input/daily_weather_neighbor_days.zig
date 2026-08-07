const std = @import("std");
const execution_calendar_date = @import("../../driver/execution_calendar_date.zig");

pub const Inputs = struct {
    current_day_of_year: u16,
    scene_begin_day_of_year: u16,
    execution_year: u16,
};

pub const NeighborDays = struct {
    previous_day_of_year: u16,
    current_day_of_year: u16,
    next_day_of_year: u16,
};

/// Exact DAY I2/I/I3 selection from day.f:257-263.
///
/// At day one, the preliminary previous record is the previous year's final
/// day. If day one is also scene begin, the following source statement
/// overrides it with the current record. At current-year end, next repeats
/// current.
pub fn select(inputs: Inputs) !NeighborDays {
    if (inputs.execution_year == 0 or
        (inputs.execution_year == 1 and inputs.current_day_of_year == 1))
        return error.InvalidDailyWeatherNeighborDay;
    _ = execution_calendar_date.fromDayOfYear(
        inputs.current_day_of_year,
        inputs.execution_year,
    ) catch return error.InvalidDailyWeatherNeighborDay;
    _ = execution_calendar_date.fromDayOfYear(
        inputs.scene_begin_day_of_year,
        inputs.execution_year,
    ) catch return error.InvalidDailyWeatherNeighborDay;
    const current_year_day_count: u16 =
        if (execution_calendar_date.isLeapYear(inputs.execution_year)) 366 else 365;
    const previous_year_day_count: u16 =
        if (inputs.execution_year > 1 and
        execution_calendar_date.isLeapYear(inputs.execution_year - 1))
            366
        else
            365;

    var previous = if (inputs.current_day_of_year == 1)
        previous_year_day_count
    else
        inputs.current_day_of_year - 1;
    if (inputs.current_day_of_year == inputs.scene_begin_day_of_year)
        previous = inputs.current_day_of_year;
    const next = if (inputs.current_day_of_year == current_year_day_count)
        inputs.current_day_of_year
    else
        inputs.current_day_of_year + 1;
    return .{
        .previous_day_of_year = previous,
        .current_day_of_year = inputs.current_day_of_year,
        .next_day_of_year = next,
    };
}

test "ordinary day selects immediate previous and next records" {
    const result = try select(.{
        .current_day_of_year = 100,
        .scene_begin_day_of_year = 1,
        .execution_year = 1999,
    });
    try std.testing.expectEqual(@as(u16, 99), result.previous_day_of_year);
    try std.testing.expectEqual(@as(u16, 100), result.current_day_of_year);
    try std.testing.expectEqual(@as(u16, 101), result.next_day_of_year);
}

test "day one uses previous year final day when scene began earlier" {
    const result = try select(.{
        .current_day_of_year = 1,
        .scene_begin_day_of_year = 10,
        .execution_year = 1901,
    });
    try std.testing.expectEqual(@as(u16, 366), result.previous_day_of_year);
    try std.testing.expectEqual(@as(u16, 2), result.next_day_of_year);
}

test "scene begin overrides day-one previous-year selection" {
    const result = try select(.{
        .current_day_of_year = 1,
        .scene_begin_day_of_year = 1,
        .execution_year = 1901,
    });
    try std.testing.expectEqual(@as(u16, 1), result.previous_day_of_year);
}

test "scene begin and year end repeat the current record independently" {
    const result = try select(.{
        .current_day_of_year = 366,
        .scene_begin_day_of_year = 366,
        .execution_year = 1900,
    });
    try std.testing.expectEqual(@as(u16, 366), result.previous_day_of_year);
    try std.testing.expectEqual(@as(u16, 366), result.next_day_of_year);
}

test "invalid runtime dates fail immediately" {
    try std.testing.expectError(
        error.InvalidDailyWeatherNeighborDay,
        select(.{
            .current_day_of_year = 366,
            .scene_begin_day_of_year = 1,
            .execution_year = 1901,
        }),
    );
    try std.testing.expectError(
        error.InvalidDailyWeatherNeighborDay,
        select(.{
            .current_day_of_year = 1,
            .scene_begin_day_of_year = 1,
            .execution_year = 1,
        }),
    );
}
