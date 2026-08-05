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
    surface_sediment_megagrams: []f64,
    disturbance_modes: []const DisturbanceMode,
) !void {
    if (surface_sediment_megagrams.len == 0 or
        surface_sediment_megagrams.len != disturbance_modes.len)
        return error.SurfaceSedimentInitializationDimensionMismatch;
    for (surface_sediment_megagrams, disturbance_modes) |sediment_megagrams, mode| {
        if (!mode.includesErosion() and !std.math.isFinite(sediment_megagrams))
            return error.NonFiniteRetainedSurfaceSediment;
        if (!mode.includesErosion() and sediment_megagrams < 0)
            return error.InvalidRetainedSurfaceSediment;
    }
    for (surface_sediment_megagrams, disturbance_modes) |*sediment_megagrams, mode| {
        if (mode.includesErosion()) sediment_megagrams.* = 0;
    }
}

test "STARTS surface sediment reset follows erosion disturbance modes" {
    var sediment_megagrams = [_]f64{ 1, 2, 3, 4, 5 };
    try initialize(&sediment_megagrams, &.{
        .no_profile_change,
        .freeze_thaw,
        .freeze_thaw_and_erosion,
        .freeze_thaw_and_organic_matter_change,
        .freeze_thaw_erosion_and_organic_matter_change,
    });
    try std.testing.expectEqualSlices(
        f64,
        &.{ 1, 2, 0, 4, 0 },
        &sediment_megagrams,
    );
}

test "STARTS surface sediment initialization supports runtime cell count" {
    var sediment_megagrams = [_]f64{9} ** 7;
    const modes = [_]DisturbanceMode{.freeze_thaw_and_erosion} ** sediment_megagrams.len;
    try initialize(&sediment_megagrams, &modes);
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0, 0, 0, 0, 0, 0, 0 },
        &sediment_megagrams,
    );
}

test "STARTS surface sediment preflight prevents partial reset" {
    var sediment_megagrams = [_]f64{ 5, std.math.nan(f64), 7 };
    try std.testing.expectError(
        error.NonFiniteRetainedSurfaceSediment,
        initialize(&sediment_megagrams, &.{
            .freeze_thaw_and_erosion,
            .freeze_thaw,
            .freeze_thaw_and_erosion,
        }),
    );
    try std.testing.expectEqual(@as(f64, 5), sediment_megagrams[0]);
    try std.testing.expect(std.math.isNan(sediment_megagrams[1]));
    try std.testing.expectEqual(@as(f64, 7), sediment_megagrams[2]);
}
