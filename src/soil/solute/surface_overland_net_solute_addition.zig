const std = @import("std");

pub const runoff_species_count = 42;
pub const snow_species_count = 41;
pub const side_count = 2;
pub const Direction = enum { east_west, north_south, vertical };

pub const Inputs = struct {
    is_surface_layer: bool,
    direction: Direction,
    /// Side-major `[side][species]`, forward then reverse.
    runoff_boundary_amount_per_step: []const f64,
    snow_boundary_amount_per_step: []const f64,
};

pub const Totals = struct {
    runoff_amount_per_step: []f64,
    snow_amount_per_step: []f64,
};

/// Compatibility translation of TRNSFRS.F lines 8529--8637.
/// For each NN side in order, all 42 runoff species are added first, followed
/// by all 41 snow-transfer species. Vertical and non-surface calls do nothing.
pub fn add(inputs: Inputs, totals: Totals) !bool {
    if (!inputs.is_surface_layer or inputs.direction == .vertical) return false;
    try validate(inputs, totals);

    for (0..runoff_species_count) |species| {
        const after_forward = totals.runoff_amount_per_step[species] +
            inputs.runoff_boundary_amount_per_step[species];
        const after_reverse = after_forward +
            inputs.runoff_boundary_amount_per_step[runoff_species_count + species];
        if (!std.math.isFinite(after_forward) or !std.math.isFinite(after_reverse))
            return error.NonFiniteSurfaceOverlandNetSoluteResult;
    }
    for (0..snow_species_count) |species| {
        const after_forward = totals.snow_amount_per_step[species] +
            inputs.snow_boundary_amount_per_step[species];
        const after_reverse = after_forward +
            inputs.snow_boundary_amount_per_step[snow_species_count + species];
        if (!std.math.isFinite(after_forward) or !std.math.isFinite(after_reverse))
            return error.NonFiniteSurfaceOverlandNetSoluteResult;
    }
    for (0..side_count) |side| {
        for (0..runoff_species_count) |species|
            totals.runoff_amount_per_step[species] = totals.runoff_amount_per_step[species] +
                inputs.runoff_boundary_amount_per_step[side * runoff_species_count + species];
        for (0..snow_species_count) |species|
            totals.snow_amount_per_step[species] = totals.snow_amount_per_step[species] +
                inputs.snow_boundary_amount_per_step[side * snow_species_count + species];
    }
    return true;
}

fn validate(inputs: Inputs, totals: Totals) !void {
    if (inputs.runoff_boundary_amount_per_step.len != side_count * runoff_species_count or
        inputs.snow_boundary_amount_per_step.len != side_count * snow_species_count or
        totals.runoff_amount_per_step.len != runoff_species_count or totals.snow_amount_per_step.len != snow_species_count)
        return error.SurfaceOverlandNetSoluteDimensionMismatch;
    for (inputs.runoff_boundary_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfaceOverlandNetSoluteInput;
    for (inputs.snow_boundary_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfaceOverlandNetSoluteInput;
    for (totals.runoff_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfaceOverlandNetSoluteInput;
    for (totals.snow_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfaceOverlandNetSoluteInput;
}

test "TRNSFRS adds forward then reverse runoff and snow boundary fluxes" {
    var runoff_boundary = [_]f64{0} ** (side_count * runoff_species_count);
    var snow_boundary = [_]f64{0} ** (side_count * snow_species_count);
    @memset(runoff_boundary[0..runoff_species_count], 2);
    @memset(runoff_boundary[runoff_species_count..], -0.5);
    @memset(snow_boundary[0..snow_species_count], 3);
    @memset(snow_boundary[snow_species_count..], -1);
    var runoff_total = [_]f64{1} ** runoff_species_count;
    var snow_total = [_]f64{1} ** snow_species_count;
    try std.testing.expect(try add(.{ .is_surface_layer = true, .direction = .east_west, .runoff_boundary_amount_per_step = &runoff_boundary, .snow_boundary_amount_per_step = &snow_boundary }, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total }));
    try std.testing.expectEqual(@as(f64, 2.5), runoff_total[0]);
    try std.testing.expectEqual(@as(f64, 2.5), runoff_total[41]);
    try std.testing.expectEqual(@as(f64, 3), snow_total[0]);
    try std.testing.expectEqual(@as(f64, 3), snow_total[40]);
}

test "TRNSFRS retains different runoff and snow species topology" {
    const runoff_boundary = [_]f64{1} ** (side_count * runoff_species_count);
    const snow_boundary = [_]f64{2} ** (side_count * snow_species_count);
    var runoff_total = [_]f64{0} ** runoff_species_count;
    var snow_total = [_]f64{0} ** snow_species_count;
    _ = try add(.{ .is_surface_layer = true, .direction = .north_south, .runoff_boundary_amount_per_step = &runoff_boundary, .snow_boundary_amount_per_step = &snow_boundary }, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total });
    try std.testing.expectEqual(@as(f64, 2), runoff_total[41]);
    try std.testing.expectEqual(@as(f64, 4), snow_total[40]);
}

test "TRNSFRS vertical direction leaves both totals untouched" {
    const runoff_boundary = [_]f64{};
    const snow_boundary = [_]f64{};
    var runoff_total = [_]f64{3} ** runoff_species_count;
    var snow_total = [_]f64{4} ** snow_species_count;
    try std.testing.expect(!try add(.{ .is_surface_layer = true, .direction = .vertical, .runoff_boundary_amount_per_step = &runoff_boundary, .snow_boundary_amount_per_step = &snow_boundary }, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total }));
    try std.testing.expectEqual(@as(f64, 3), runoff_total[0]);
    try std.testing.expectEqual(@as(f64, 4), snow_total[0]);
}

test "TRNSFRS non-surface guard leaves invalid dormant topology unexamined" {
    const empty = [_]f64{};
    try std.testing.expect(!try add(.{ .is_surface_layer = false, .direction = .east_west, .runoff_boundary_amount_per_step = &empty, .snow_boundary_amount_per_step = &empty }, .{ .runoff_amount_per_step = @constCast(&empty), .snow_amount_per_step = @constCast(&empty) }));
}

test "runtime topology mismatch fails atomically" {
    const short_runoff = [_]f64{1} ** (side_count * runoff_species_count - 1);
    const snow_boundary = [_]f64{2} ** (side_count * snow_species_count);
    var runoff_total = [_]f64{3} ** runoff_species_count;
    var snow_total = [_]f64{4} ** snow_species_count;
    try std.testing.expectError(error.SurfaceOverlandNetSoluteDimensionMismatch, add(.{ .is_surface_layer = true, .direction = .east_west, .runoff_boundary_amount_per_step = &short_runoff, .snow_boundary_amount_per_step = &snow_boundary }, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total }));
    try std.testing.expectEqual(@as(f64, 3), runoff_total[0]);
}

test "late invalid snow flux leaves runoff total atomic" {
    const runoff_boundary = [_]f64{1} ** (side_count * runoff_species_count);
    var snow_boundary = [_]f64{2} ** (side_count * snow_species_count);
    snow_boundary[snow_boundary.len - 1] = std.math.inf(f64);
    var runoff_total = [_]f64{3} ** runoff_species_count;
    var snow_total = [_]f64{4} ** snow_species_count;
    try std.testing.expectError(error.NonFiniteSurfaceOverlandNetSoluteInput, add(.{ .is_surface_layer = true, .direction = .east_west, .runoff_boundary_amount_per_step = &runoff_boundary, .snow_boundary_amount_per_step = &snow_boundary }, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total }));
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** runoff_species_count), &runoff_total);
}

test "cumulative two-side overflow fails atomically" {
    var runoff_boundary = [_]f64{0} ** (side_count * runoff_species_count);
    const snow_boundary = [_]f64{0} ** (side_count * snow_species_count);
    runoff_boundary[0] = std.math.floatMax(f64) * 0.75;
    runoff_boundary[runoff_species_count] = std.math.floatMax(f64) * 0.75;
    var runoff_total = [_]f64{0} ** runoff_species_count;
    var snow_total = [_]f64{0} ** snow_species_count;
    try std.testing.expectError(error.NonFiniteSurfaceOverlandNetSoluteResult, add(.{ .is_surface_layer = true, .direction = .east_west, .runoff_boundary_amount_per_step = &runoff_boundary, .snow_boundary_amount_per_step = &snow_boundary }, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total }));
    try std.testing.expectEqual(@as(f64, 0), runoff_total[0]);
}
