const std = @import("std");

pub const TemporalResolution = enum(u8) {
    daily = 1,
    hourly = 2,
};

pub const Sources = struct {
    first: TemporalResolution,
    second: TemporalResolution,
};

pub const Decision = struct {
    resolution: TemporalResolution,
    source_index: usize,
};

/// Exact DAY weather source decision:
/// IGO==0 or I<=ILAST selects IWTHR(1), otherwise IWTHR(2).
pub fn select(
    first_execution_scene: bool,
    execution_day: u64,
    first_source_last_day: u64,
    sources: Sources,
) !Decision {
    if (execution_day == 0)
        return error.InvalidWeatherSourceExecutionDay;
    const first = first_execution_scene or
        execution_day <= first_source_last_day;
    return .{
        .resolution = if (first) sources.first else sources.second,
        .source_index = if (first) 0 else 1,
    };
}

pub fn resolutionFromCode(code: u8) !TemporalResolution {
    return switch (code) {
        1 => .daily,
        2 => .hourly,
        else => error.InvalidWeatherTemporalResolution,
    };
}

test "first execution scene always selects first weather source" {
    const decision = try select(true, 500, 100, .{
        .first = .daily,
        .second = .hourly,
    });
    try std.testing.expectEqual(@as(usize, 0), decision.source_index);
    try std.testing.expectEqual(TemporalResolution.daily, decision.resolution);
}

test "later execution switches only after exact ILAST boundary" {
    const sources: Sources = .{
        .first = .daily,
        .second = .hourly,
    };
    var decision = try select(false, 100, 100, sources);
    try std.testing.expectEqual(@as(usize, 0), decision.source_index);
    decision = try select(false, 101, 100, sources);
    try std.testing.expectEqual(@as(usize, 1), decision.source_index);
    try std.testing.expectEqual(TemporalResolution.hourly, decision.resolution);
}

test "runtime weather codes reject unsupported temporal modes" {
    try std.testing.expectEqual(
        TemporalResolution.daily,
        try resolutionFromCode(1),
    );
    try std.testing.expectEqual(
        TemporalResolution.hourly,
        try resolutionFromCode(2),
    );
    try std.testing.expectError(
        error.InvalidWeatherTemporalResolution,
        resolutionFromCode(0),
    );
    try std.testing.expectError(
        error.InvalidWeatherTemporalResolution,
        resolutionFromCode(3),
    );
}

test "zero execution day cannot choose an out-of-core stream" {
    try std.testing.expectError(
        error.InvalidWeatherSourceExecutionDay,
        select(false, 0, 100, .{
            .first = .daily,
            .second = .hourly,
        }),
    );
}
