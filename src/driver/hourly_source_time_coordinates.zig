const std = @import("std");

pub const Inputs = struct {
    source_day_of_year: u16,
    source_hour_one_through_twenty_four: u8,
    source_subhourly_iteration: u32,
};

pub const Result = struct {
    source_hour_real: f32,
    fractional_day_of_year: f32,
    source_subhourly_iteration_real: f32,
};

/// `hour1.f` lines 130--132. Implicit Fortran default-REAL temporaries remain
/// binary32 and preserve source assignment and operation order.
pub fn compute(inputs: Inputs) !Result {
    if (inputs.source_day_of_year == 0 or
        inputs.source_day_of_year > 366 or
        inputs.source_hour_one_through_twenty_four == 0 or
        inputs.source_hour_one_through_twenty_four > 24 or
        inputs.source_subhourly_iteration == 0)
        return error.InvalidHourlySourceTimeCoordinate;

    var result: Result = undefined;
    result.source_hour_real =
        @floatFromInt(inputs.source_hour_one_through_twenty_four);
    result.fractional_day_of_year =
        @as(f32, @floatFromInt(inputs.source_day_of_year - 1)) +
        result.source_hour_real / 24;
    result.source_subhourly_iteration_real =
        @floatFromInt(inputs.source_subhourly_iteration);
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteHourlySourceTimeCoordinate;
    return result;
}

test "HOUR1 source time coordinates preserve one-based hour convention" {
    const result = try compute(.{
        .source_day_of_year = 100,
        .source_hour_one_through_twenty_four = 6,
        .source_subhourly_iteration = 3,
    });
    try std.testing.expectEqual(@as(f32, 6), result.source_hour_real);
    try std.testing.expectEqual(@as(f32, 99.25), result.fractional_day_of_year);
    try std.testing.expectEqual(
        @as(f32, 3),
        result.source_subhourly_iteration_real,
    );
}

test "source hour twenty four reaches integer day coordinate" {
    const result = try compute(.{
        .source_day_of_year = 366,
        .source_hour_one_through_twenty_four = 24,
        .source_subhourly_iteration = 1,
    });
    try std.testing.expectEqual(@as(f32, 366), result.fractional_day_of_year);
}

test "invalid source coordinates fail explicitly" {
    try std.testing.expectError(
        error.InvalidHourlySourceTimeCoordinate,
        compute(.{
            .source_day_of_year = 1,
            .source_hour_one_through_twenty_four = 0,
            .source_subhourly_iteration = 1,
        }),
    );
    try std.testing.expectError(
        error.InvalidHourlySourceTimeCoordinate,
        compute(.{
            .source_day_of_year = 367,
            .source_hour_one_through_twenty_four = 1,
            .source_subhourly_iteration = 1,
        }),
    );
}
