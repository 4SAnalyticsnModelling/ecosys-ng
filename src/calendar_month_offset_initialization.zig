const std = @import("std");

/// Exact translation of `BLOCKDATA001.f` line 3.
///
/// The caller owns and runtime-allocates the twelve Gregorian month offsets.
/// Values convert the source's `30 * month + offset` approximation into its
/// calendar boundary convention; they are dimensionless day counts.
pub fn initialize(month_offset_days: []i8) !void {
    if (month_offset_days.len != 12)
        return error.CalendarMonthOffsetDimensionMismatch;
    const source_offsets = [_]i8{ 1, -1, 0, 0, 1, 1, 2, 3, 3, 4, 4, 5 };
    @memcpy(month_offset_days, &source_offsets);
}

test "BLOCKDATA001 month offsets preserve exact source order" {
    var offsets = [_]i8{99} ** 12;
    try initialize(&offsets);
    try std.testing.expectEqualSlices(
        i8,
        &.{ 1, -1, 0, 0, 1, 1, 2, 3, 3, 4, 4, 5 },
        &offsets,
    );
}

test "month offset dimension failure is atomic" {
    var offsets = [_]i8{99} ** 11;
    try std.testing.expectError(
        error.CalendarMonthOffsetDimensionMismatch,
        initialize(&offsets),
    );
    const expected = [_]i8{99} ** 11;
    try std.testing.expectEqualSlices(i8, &expected, &offsets);
}
