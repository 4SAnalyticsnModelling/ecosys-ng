const std = @import("std");

pub const State = struct {
    reported_day: i64,
    previous_day: i64,
    management_event_count: u64,
    year_transition_count: u64,
};

pub const Inputs = struct {
    current_day: i64,
    requested_reported_day: i64,
    days_in_current_year: i64,
};

/// Atomic translation of the EXEC tail. A negative requested day remains
/// relative to the runtime year length; otherwise EXEC publishes current_day.
/// Daily management and year-transition counters clear only after the next
/// complete state has been validated.
pub fn advance(state: *State, inputs: Inputs) !void {
    if (inputs.current_day < 0)
        return error.InvalidExecutionCurrentDay;
    if (inputs.days_in_current_year <= 0)
        return error.InvalidExecutionYearLength;

    const reported_day = if (inputs.requested_reported_day < 0)
        std.math.add(
            i64,
            inputs.days_in_current_year,
            inputs.requested_reported_day,
        ) catch return error.ExecutionReportedDayOverflow
    else
        inputs.current_day;
    if (reported_day < 0 or reported_day > inputs.days_in_current_year)
        return error.ExecutionReportedDayOutOfRange;

    state.* = .{
        .reported_day = reported_day,
        .previous_day = inputs.current_day,
        .management_event_count = 0,
        .year_transition_count = 0,
    };
}

test "negative IDAYR remains relative to runtime common and leap years" {
    var state: State = .{
        .reported_day = 50,
        .previous_day = 49,
        .management_event_count = 8,
        .year_transition_count = 2,
    };
    try advance(&state, .{
        .current_day = 100,
        .requested_reported_day = -2,
        .days_in_current_year = 365,
    });
    try std.testing.expectEqual(@as(i64, 363), state.reported_day);
    try std.testing.expectEqual(@as(i64, 100), state.previous_day);
    try std.testing.expectEqual(@as(u64, 0), state.management_event_count);
    try std.testing.expectEqual(@as(u64, 0), state.year_transition_count);

    state.management_event_count = 1;
    try advance(&state, .{
        .current_day = 101,
        .requested_reported_day = -2,
        .days_in_current_year = 366,
    });
    try std.testing.expectEqual(@as(i64, 364), state.reported_day);
}

test "nonnegative IDAYR publishes current execution day" {
    var state: State = .{
        .reported_day = 1,
        .previous_day = 1,
        .management_event_count = 5,
        .year_transition_count = 3,
    };
    try advance(&state, .{
        .current_day = 207,
        .requested_reported_day = 9,
        .days_in_current_year = 365,
    });
    try std.testing.expectEqual(@as(i64, 207), state.reported_day);
    try std.testing.expectEqual(@as(i64, 207), state.previous_day);
}

test "invalid relative day cannot partially clear daily counters" {
    var state: State = .{
        .reported_day = 70,
        .previous_day = 69,
        .management_event_count = 12,
        .year_transition_count = 4,
    };
    const before = state;
    try std.testing.expectError(
        error.ExecutionReportedDayOutOfRange,
        advance(&state, .{
            .current_day = 100,
            .requested_reported_day = -400,
            .days_in_current_year = 365,
        }),
    );
    try std.testing.expectEqualDeep(before, state);
}

test "invalid current day rolls back state" {
    var state: State = .{
        .reported_day = 2,
        .previous_day = 1,
        .management_event_count = 2,
        .year_transition_count = 1,
    };
    const before = state;
    try std.testing.expectError(
        error.InvalidExecutionCurrentDay,
        advance(&state, .{
            .current_day = -1,
            .requested_reported_day = -1,
            .days_in_current_year = 365,
        }),
    );
    try std.testing.expectEqualDeep(before, state);
}
