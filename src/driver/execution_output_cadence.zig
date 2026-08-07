const std = @import("std");

pub const Schedule = struct {
    interval_days: u64,
};

/// Exact positive-day equivalent of EXEC `(I / IOUT) * IOUT == I`.
/// Modulo avoids overflow in the multiply-back expression.
pub fn isDue(execution_day: u64, schedule: Schedule) !bool {
    if (execution_day == 0)
        return error.InvalidExecutionOutputDay;
    if (schedule.interval_days == 0)
        return error.InvalidExecutionOutputInterval;
    return execution_day % schedule.interval_days == 0;
}

/// Validates an entire runtime sequence before publishing any due flags.
/// This is useful when constructing an out-of-core scene output schedule.
pub fn fillDueFlags(
    execution_days: []const u64,
    schedule: Schedule,
    due: []bool,
) !void {
    if (execution_days.len != due.len)
        return error.ExecutionOutputCadenceDimensionMismatch;
    if (schedule.interval_days == 0)
        return error.InvalidExecutionOutputInterval;
    for (execution_days) |day| if (day == 0)
        return error.InvalidExecutionOutputDay;
    for (execution_days, due) |day, *flag|
        flag.* = day % schedule.interval_days == 0;
}

test "EXEC cadence selects exact positive multiples" {
    const schedule: Schedule = .{ .interval_days = 10 };
    try std.testing.expect(!try isDue(1, schedule));
    try std.testing.expect(!try isDue(9, schedule));
    try std.testing.expect(try isDue(10, schedule));
    try std.testing.expect(try isDue(20, schedule));
    try std.testing.expect(!try isDue(21, schedule));
}

test "modulo cadence remains safe at maximum runtime day" {
    const maximum = std.math.maxInt(u64);
    try std.testing.expect(try isDue(maximum, .{ .interval_days = maximum }));
    try std.testing.expect(!try isDue(
        maximum - 1,
        .{ .interval_days = maximum },
    ));
}

test "invalid cadence controls fail immediately" {
    try std.testing.expectError(
        error.InvalidExecutionOutputDay,
        isDue(0, .{ .interval_days = 1 }),
    );
    try std.testing.expectError(
        error.InvalidExecutionOutputInterval,
        isDue(1, .{ .interval_days = 0 }),
    );
}

test "runtime schedule fill rolls back on invalid late day" {
    var due = [_]bool{ true, false, true, false };
    const before = due;
    try std.testing.expectError(
        error.InvalidExecutionOutputDay,
        fillDueFlags(
            &.{ 10, 20, 0, 40 },
            .{ .interval_days = 10 },
            &due,
        ),
    );
    try std.testing.expectEqualSlices(bool, &before, &due);
    try fillDueFlags(
        &.{ 9, 10, 19, 20 },
        .{ .interval_days = 10 },
        &due,
    );
    try std.testing.expectEqualSlices(
        bool,
        &[_]bool{ false, true, false, true },
        &due,
    );
}
