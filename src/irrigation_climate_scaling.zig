const std = @import("std");

pub const Totals = struct {
    surface_before_m: f64,
    subsurface_before_m: f64,
    surface_after_m: f64,
    subsurface_after_m: f64,
};

/// Applies WTHR's active seasonal TDIRI multiplier to both irrigation
/// pathways. Validation of all runtime cells and products precedes mutation,
/// so chemistry can subsequently use the scaled carrier without partial mass.
pub fn apply(
    surface_irrigation_depth_m_by_cell: []f64,
    subsurface_irrigation_depth_m_by_cell: []f64,
    irrigation_multiplier_by_cell: []const f64,
) !Totals {
    const cell_count = surface_irrigation_depth_m_by_cell.len;
    if (cell_count == 0 or
        subsurface_irrigation_depth_m_by_cell.len != cell_count or
        irrigation_multiplier_by_cell.len != cell_count)
        return error.IrrigationClimateScalingDimensionMismatch;

    var totals: Totals = .{
        .surface_before_m = 0,
        .subsurface_before_m = 0,
        .surface_after_m = 0,
        .subsurface_after_m = 0,
    };
    for (
        surface_irrigation_depth_m_by_cell,
        subsurface_irrigation_depth_m_by_cell,
        irrigation_multiplier_by_cell,
    ) |surface_m, subsurface_m, multiplier| {
        inline for (.{ surface_m, subsurface_m, multiplier }) |value|
            if (!std.math.isFinite(value))
                return error.NonFiniteIrrigationClimateScalingInput;
        if (surface_m < 0 or subsurface_m < 0 or multiplier < 0)
            return error.InvalidIrrigationClimateScalingInput;
        const scaled_surface_m = surface_m * multiplier;
        const scaled_subsurface_m = subsurface_m * multiplier;
        inline for (.{ scaled_surface_m, scaled_subsurface_m }) |value|
            if (!std.math.isFinite(value))
                return error.IrrigationClimateScalingOverflow;
        totals.surface_before_m += surface_m;
        totals.subsurface_before_m += subsurface_m;
        totals.surface_after_m += scaled_surface_m;
        totals.subsurface_after_m += scaled_subsurface_m;
        inline for (@typeInfo(Totals).@"struct".fields) |field|
            if (!std.math.isFinite(@field(totals, field.name)))
                return error.IrrigationClimateScalingOverflow;
    }

    for (
        surface_irrigation_depth_m_by_cell,
        subsurface_irrigation_depth_m_by_cell,
        irrigation_multiplier_by_cell,
    ) |*surface_m, *subsurface_m, multiplier| {
        surface_m.* *= multiplier;
        subsurface_m.* *= multiplier;
    }
    return totals;
}

test "seasonal irrigation multiplier scales both runtime pathways" {
    var surface = [_]f64{ 0.001, 0.002, 0 };
    var subsurface = [_]f64{ 0.004, 0, 0.006 };
    const totals = try apply(
        &surface,
        &subsurface,
        &.{ 2, 0.5, 1.5 },
    );
    try std.testing.expectApproxEqAbs(@as(f64, 0.003), totals.surface_before_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.010), totals.subsurface_before_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.003), totals.surface_after_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.017), totals.subsurface_after_m, 1e-15);
    try std.testing.expectEqualSlices(
        f64,
        &[_]f64{ 0.002, 0.001, 0 },
        &surface,
    );
    for (subsurface, [_]f64{ 0.008, 0, 0.009 }) |actual, expected|
        try std.testing.expectApproxEqAbs(expected, actual, 1e-15);
}

test "zero climate multiplier disables both irrigation carriers" {
    var surface = [_]f64{0.001};
    var subsurface = [_]f64{0.002};
    _ = try apply(&surface, &subsurface, &.{0});
    try std.testing.expectEqual(@as(f64, 0), surface[0]);
    try std.testing.expectEqual(@as(f64, 0), subsurface[0]);
}

test "invalid late cell rolls back every irrigation carrier" {
    var surface = [_]f64{ 1, 2, 3 };
    var subsurface = [_]f64{ 4, 5, 6 };
    const surface_before = surface;
    const subsurface_before = subsurface;
    try std.testing.expectError(
        error.NonFiniteIrrigationClimateScalingInput,
        apply(
            &surface,
            &subsurface,
            &.{ 1, 2, std.math.nan(f64) },
        ),
    );
    try std.testing.expectEqualSlices(f64, &surface_before, &surface);
    try std.testing.expectEqualSlices(f64, &subsurface_before, &subsurface);
}

test "overflowing product rolls back runtime arrays" {
    var surface = [_]f64{std.math.floatMax(f64)};
    var subsurface = [_]f64{1};
    try std.testing.expectError(
        error.IrrigationClimateScalingOverflow,
        apply(&surface, &subsurface, &.{2}),
    );
    try std.testing.expectEqual(std.math.floatMax(f64), surface[0]);
    try std.testing.expectEqual(@as(f64, 1), subsurface[0]);
}
