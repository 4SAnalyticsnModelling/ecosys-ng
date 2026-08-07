const std = @import("std");

/// Exact grosub.f lines 303--304 canopy-height rollover (ZCX=ZC; ZC=0).
///
/// Both slices use the production cell-major, species-minor runtime plant
/// order. Validation precedes mutation so a corrupt late plant cannot expose
/// a mixture of previous and current hourly canopy heights.
pub fn rollover(
    previous_canopy_height_m_by_plant: []f64,
    current_canopy_height_m_by_plant: []f64,
) !void {
    if (current_canopy_height_m_by_plant.len == 0 or
        previous_canopy_height_m_by_plant.len !=
            current_canopy_height_m_by_plant.len)
        return error.CanopyHeightRolloverDimensionMismatch;
    for (current_canopy_height_m_by_plant) |height_m| {
        if (!std.math.isFinite(height_m))
            return error.NonFiniteCanopyHeight;
        if (height_m < 0)
            return error.NegativeCanopyHeight;
    }
    for (previous_canopy_height_m_by_plant, current_canopy_height_m_by_plant) |*previous_height_m, *current_height_m| {
        previous_height_m.* = current_height_m.*;
        current_height_m.* = 0;
    }
}

test "GROSUB rolls arbitrary runtime plant heights in source order" {
    const allocator = std.testing.allocator;
    const plant_count = 41;
    const previous = try allocator.alloc(f64, plant_count);
    defer allocator.free(previous);
    const current = try allocator.alloc(f64, plant_count);
    defer allocator.free(current);
    @memset(previous, -1);
    for (current, 0..) |*height_m, plant|
        height_m.* = @as(f64, @floatFromInt(plant)) * 0.125;

    try rollover(previous, current);

    for (previous, current, 0..) |previous_height_m, current_height_m, plant| {
        try std.testing.expectEqual(
            @as(f64, @floatFromInt(plant)) * 0.125,
            previous_height_m,
        );
        try std.testing.expectEqual(@as(f64, 0), current_height_m);
    }
}

test "invalid late canopy height leaves both generations unchanged" {
    var previous = [_]f64{ 1, 2, 3 };
    var current = [_]f64{ 4, 5, std.math.nan(f64) };
    const previous_before = previous;
    const current_before = current;

    try std.testing.expectError(
        error.NonFiniteCanopyHeight,
        rollover(&previous, &current),
    );

    try std.testing.expectEqualSlices(f64, &previous_before, &previous);
    try std.testing.expectEqualSlices(
        u8,
        std.mem.asBytes(&current_before),
        std.mem.asBytes(&current),
    );
}

test "dimension mismatch and negative height fail explicitly" {
    var previous = [_]f64{ 1, 2 };
    var short_current = [_]f64{3};
    try std.testing.expectError(
        error.CanopyHeightRolloverDimensionMismatch,
        rollover(&previous, &short_current),
    );
    var negative_current = [_]f64{ -0.1, 0 };
    try std.testing.expectError(
        error.NegativeCanopyHeight,
        rollover(&previous, &negative_current),
    );
}
