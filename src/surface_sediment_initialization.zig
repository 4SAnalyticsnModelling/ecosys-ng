const std = @import("std");

pub const DisturbanceMode = enum {
    no_profile_change,
    freeze_thaw,
    freeze_thaw_and_erosion,
    freeze_thaw_and_organic_matter_change,
    freeze_thaw_erosion_and_organic_matter_change,

    fn includesErosion(self: DisturbanceMode) bool {
        return self == .freeze_thaw_and_erosion or
            self == .freeze_thaw_erosion_and_organic_matter_change;
    }
};

/// Direct translation of STARTS lines 677--679 over runtime grid cells.
/// Surface-water sediment inventory is extensive Mg. Non-erosion modes retain
/// their initialized values, matching the source branch.
pub fn initialize(
    surface_sediment_Mg: []f64,
    disturbance_modes: []const DisturbanceMode,
) !void {
    if (surface_sediment_Mg.len == 0 or
        surface_sediment_Mg.len != disturbance_modes.len)
        return error.SurfaceSedimentInitializationDimensionMismatch;
    for (surface_sediment_Mg, disturbance_modes) |sediment_Mg, mode| {
        if (!mode.includesErosion() and !std.math.isFinite(sediment_Mg))
            return error.NonFiniteRetainedSurfaceSediment;
        if (!mode.includesErosion() and sediment_Mg < 0)
            return error.InvalidRetainedSurfaceSediment;
    }
    for (surface_sediment_Mg, disturbance_modes) |*sediment_Mg, mode| {
        if (mode.includesErosion()) sediment_Mg.* = 0;
    }
}

test "STARTS surface sediment reset follows erosion disturbance modes" {
    var sediment_Mg = [_]f64{ 1, 2, 3, 4, 5 };
    try initialize(&sediment_Mg, &.{
        .no_profile_change,
        .freeze_thaw,
        .freeze_thaw_and_erosion,
        .freeze_thaw_and_organic_matter_change,
        .freeze_thaw_erosion_and_organic_matter_change,
    });
    try std.testing.expectEqualSlices(
        f64,
        &.{ 1, 2, 0, 4, 0 },
        &sediment_Mg,
    );
}

test "STARTS surface sediment initialization supports runtime cell count" {
    var sediment_Mg = [_]f64{9} ** 7;
    const modes = [_]DisturbanceMode{.freeze_thaw_and_erosion} ** sediment_Mg.len;
    try initialize(&sediment_Mg, &modes);
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0, 0, 0, 0, 0, 0, 0 },
        &sediment_Mg,
    );
}

test "STARTS surface sediment preflight prevents partial reset" {
    var sediment_Mg = [_]f64{ 5, std.math.nan(f64), 7 };
    try std.testing.expectError(
        error.NonFiniteRetainedSurfaceSediment,
        initialize(&sediment_Mg, &.{
            .freeze_thaw_and_erosion,
            .freeze_thaw,
            .freeze_thaw_and_erosion,
        }),
    );
    try std.testing.expectEqual(@as(f64, 5), sediment_Mg[0]);
    try std.testing.expect(std.math.isNan(sediment_Mg[1]));
    try std.testing.expectEqual(@as(f64, 7), sediment_Mg[2]);
}
