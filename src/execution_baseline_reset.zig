const std = @import("std");

pub const ResetReason = packed struct {
    execution_begin: bool = false,
    scene_start: bool = false,
    after_scene_end: bool = false,

    pub fn any(self: ResetReason) bool {
        return self.execution_begin or self.scene_start or
            self.after_scene_end;
    }
};

pub const Schedule = struct {
    execution_begin_day: u64,
    scene_start_day: u64,
    scene_last_day: u64,
};

/// Resolves the three exact EXEC reset comparisons without allowing ILAST+1
/// to wrap at runtime.
pub fn reason(current_day: u64, schedule: Schedule) !ResetReason {
    const after_scene_end = std.math.add(
        u64,
        schedule.scene_last_day,
        1,
    ) catch return error.ExecutionResetDayOverflow;
    return .{
        .execution_begin = current_day == schedule.execution_begin_day,
        .scene_start = current_day == schedule.scene_start_day,
        .after_scene_end = current_day == after_scene_end,
    };
}

/// Atomically replaces every runtime balance-domain baseline when EXEC would
/// reset TLW/TLH/TLO/TLC/TLN/TLP/TLI. This must run before the cadence audit
/// for the same day.
pub fn resetIfScheduled(
    baseline_by_domain: []f64,
    current_balance_by_domain: []const f64,
    current_day: u64,
    schedule: Schedule,
) !ResetReason {
    if (baseline_by_domain.len == 0 or
        baseline_by_domain.len != current_balance_by_domain.len)
        return error.ExecutionBaselineDimensionMismatch;
    const reset_reason = try reason(current_day, schedule);
    if (!reset_reason.any()) return reset_reason;

    // Validate the complete source before mutation so a non-finite late
    // domain cannot leave a partially reset baseline.
    for (current_balance_by_domain) |balance| {
        if (!std.math.isFinite(balance))
            return error.NonFiniteExecutionBalance;
    }
    @memcpy(baseline_by_domain, current_balance_by_domain);
    return reset_reason;
}

test "EXEC resets before audit at each exact boundary" {
    const schedule: Schedule = .{
        .execution_begin_day = 2,
        .scene_start_day = 10,
        .scene_last_day = 20,
    };
    const current = [_]f64{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    inline for (.{ @as(u64, 2), @as(u64, 10), @as(u64, 21) }) |day| {
        var baseline = [_]f64{0} ** current.len;
        const reset_reason = try resetIfScheduled(
            &baseline,
            &current,
            day,
            schedule,
        );
        try std.testing.expect(reset_reason.any());
        try std.testing.expectEqualSlices(f64, &current, &baseline);
    }
}

test "overlapping EXEC reset reasons remain explicit" {
    const reset_reason = try reason(5, .{
        .execution_begin_day = 5,
        .scene_start_day = 5,
        .scene_last_day = 4,
    });
    try std.testing.expect(reset_reason.execution_begin);
    try std.testing.expect(reset_reason.scene_start);
    try std.testing.expect(reset_reason.after_scene_end);
}

test "ordinary day leaves baseline byte-for-byte unchanged" {
    var baseline = [_]f64{ 10, 20, 30 };
    const before = baseline;
    const reset_reason = try resetIfScheduled(
        &baseline,
        &.{ 90, 80, 70 },
        8,
        .{
            .execution_begin_day = 1,
            .scene_start_day = 2,
            .scene_last_day = 20,
        },
    );
    try std.testing.expect(!reset_reason.any());
    try std.testing.expectEqualSlices(f64, &before, &baseline);
}

test "invalid late domain and end-day overflow roll back baseline" {
    var baseline = [_]f64{ 10, 20, 30 };
    const before = baseline;
    try std.testing.expectError(
        error.NonFiniteExecutionBalance,
        resetIfScheduled(
            &baseline,
            &.{ 1, 2, std.math.nan(f64) },
            4,
            .{
                .execution_begin_day = 4,
                .scene_start_day = 8,
                .scene_last_day = 9,
            },
        ),
    );
    try std.testing.expectEqualSlices(f64, &before, &baseline);
    try std.testing.expectError(
        error.ExecutionResetDayOverflow,
        resetIfScheduled(
            &baseline,
            &.{ 1, 2, 3 },
            4,
            .{
                .execution_begin_day = 4,
                .scene_start_day = 8,
                .scene_last_day = std.math.maxInt(u64),
            },
        ),
    );
    try std.testing.expectEqualSlices(f64, &before, &baseline);
}
