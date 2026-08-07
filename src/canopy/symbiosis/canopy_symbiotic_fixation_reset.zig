const std = @import("std");

/// Exact grosub.f line 302 hourly reset of UPNFC (g N h-1).
///
/// The slice is runtime plant-sized across all grid cells and species. The
/// subsequent canopy symbiont kernels accumulate accepted fixation into the
/// same owner during this biological hour.
pub fn resetHourly(
    canopy_symbiotic_nitrogen_fixation_g_n_per_h_by_plant: []f64,
) !void {
    if (canopy_symbiotic_nitrogen_fixation_g_n_per_h_by_plant.len == 0)
        return error.EmptyCanopySymbioticFixationState;
    @memset(canopy_symbiotic_nitrogen_fixation_g_n_per_h_by_plant, 0);
}

test "GROSUB resets canopy fixation for arbitrary runtime plant count" {
    const allocator = std.testing.allocator;
    const plant_count = 37;
    const fixation = try allocator.alloc(f64, plant_count);
    defer allocator.free(fixation);
    for (fixation, 0..) |*value, plant| value.* = @floatFromInt(plant + 1);

    try resetHourly(fixation);

    for (fixation) |value|
        try std.testing.expectEqual(@as(f64, 0), value);
}

test "empty runtime plant domain fails explicitly" {
    try std.testing.expectError(
        error.EmptyCanopySymbioticFixationState,
        resetHourly(&.{}),
    );
}
