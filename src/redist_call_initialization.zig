const std = @import("std");

pub const GridWindow = struct {
    first_column: usize,
    last_column_inclusive: usize,
    first_row: usize,
    last_row_inclusive: usize,
};

/// Exact meaningful REDIST line 204 call-local reset.
///
/// Traceability: `DORGE(NY,NX)=0` in the source's column-outer, row-inner
/// traversal. Validation completes before any caller-owned value changes.
pub fn resetErodedOrganicCarbon(
    eroded_organic_carbon_g_c_by_cell: []f64,
    column_count: usize,
    row_count: usize,
    window: GridWindow,
) !void {
    if (column_count == 0 or row_count == 0)
        return error.InvalidRedistGridDimensions;
    const cell_count = std.math.mul(usize, column_count, row_count) catch
        return error.InvalidRedistGridDimensions;
    if (eroded_organic_carbon_g_c_by_cell.len != cell_count or
        window.first_column > window.last_column_inclusive or
        window.first_row > window.last_row_inclusive or
        window.last_column_inclusive >= column_count or
        window.last_row_inclusive >= row_count)
        return error.InvalidRedistGridWindow;
    for (eroded_organic_carbon_g_c_by_cell) |value|
        if (!std.math.isFinite(value))
            return error.InvalidErodedOrganicCarbonAccumulator;

    for (window.first_column..window.last_column_inclusive + 1) |column| {
        for (window.first_row..window.last_row_inclusive + 1) |row|
            eroded_organic_carbon_g_c_by_cell[row * column_count + column] = 0;
    }
}

test "REDIST call reset preserves column outer row inner runtime window" {
    var carbon = [_]f64{
        1, 2,  3,  4,
        5, 6,  7,  8,
        9, 10, 11, 12,
    };
    try resetErodedOrganicCarbon(&carbon, 4, 3, .{
        .first_column = 1,
        .last_column_inclusive = 2,
        .first_row = 0,
        .last_row_inclusive = 1,
    });
    try std.testing.expectEqualSlices(
        f64,
        &.{ 1, 0, 0, 4, 5, 0, 0, 8, 9, 10, 11, 12 },
        &carbon,
    );
}

test "REDIST call reset rejects invalid late state atomically" {
    var carbon = [_]f64{ 1, 2, 3, std.math.nan(f64) };
    try std.testing.expectError(
        error.InvalidErodedOrganicCarbonAccumulator,
        resetErodedOrganicCarbon(&carbon, 2, 2, .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        }),
    );
    try std.testing.expectEqual(@as(f64, 1), carbon[0]);
    try std.testing.expect(std.math.isNan(carbon[3]));
}

test "REDIST call reset rejects reversed and out-of-range windows" {
    var carbon = [_]f64{ 1, 2, 3, 4 };
    try std.testing.expectError(
        error.InvalidRedistGridWindow,
        resetErodedOrganicCarbon(&carbon, 2, 2, .{
            .first_column = 1,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 1,
        }),
    );
    try std.testing.expectError(
        error.InvalidRedistGridWindow,
        resetErodedOrganicCarbon(&carbon, 2, 2, .{
            .first_column = 0,
            .last_column_inclusive = 2,
            .first_row = 0,
            .last_row_inclusive = 1,
        }),
    );
    try std.testing.expectEqualSlices(f64, &.{ 1, 2, 3, 4 }, &carbon);
}
