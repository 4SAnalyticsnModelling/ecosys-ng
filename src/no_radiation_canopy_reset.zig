pub const SpeciesDiagnostics = struct {
    live_shortwave_mj_h: []f64,
    standing_dead_shortwave_mj_h: []f64,
    live_par_umol_s: []f64,
    standing_dead_par_umol_s: []f64,
};

pub const Result = struct {
    ground_shortwave_mj_h: f64 = 0,
};

/// HOUR1 control lines 1763--1767 and executable lines 1768--1774. Resets
/// ground radiation, then runtime-species canopy diagnostics in source order.
pub fn apply(diagnostics: SpeciesDiagnostics) !Result {
    const species_count = diagnostics.live_shortwave_mj_h.len;
    if (species_count == 0 or
        diagnostics.standing_dead_shortwave_mj_h.len != species_count or
        diagnostics.live_par_umol_s.len != species_count or
        diagnostics.standing_dead_par_umol_s.len != species_count)
        return error.NoRadiationCanopyDimensionMismatch;
    const result: Result = .{};
    for (0..species_count) |species| {
        diagnostics.live_shortwave_mj_h[species] = 0.0;
        diagnostics.standing_dead_shortwave_mj_h[species] = 0.0;
        diagnostics.live_par_umol_s[species] = 0.0;
        diagnostics.standing_dead_par_umol_s[species] = 0.0;
    }
    return result;
}

test "no radiation resets ground and runtime species diagnostics" {
    var live_sw = [_]f64{ 1, 2, 3, 4 };
    var dead_sw = [_]f64{ 5, 6, 7, 8 };
    var live_par = [_]f64{ 9, 10, 11, 12 };
    var dead_par = [_]f64{ 13, 14, 15, 16 };
    const result = try apply(.{
        .live_shortwave_mj_h = &live_sw,
        .standing_dead_shortwave_mj_h = &dead_sw,
        .live_par_umol_s = &live_par,
        .standing_dead_par_umol_s = &dead_par,
    });
    try std.testing.expectEqual(
        @as(f64, 0),
        result.ground_shortwave_mj_h,
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0, 0, 0, 0 },
        &live_sw,
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 0, 0, 0, 0 },
        &dead_par,
    );
}

test "dimension mismatch leaves diagnostics unchanged" {
    var live_sw = [_]f64{42};
    var empty = [_]f64{};
    try std.testing.expectError(
        error.NoRadiationCanopyDimensionMismatch,
        apply(.{
            .live_shortwave_mj_h = &live_sw,
            .standing_dead_shortwave_mj_h = &empty,
            .live_par_umol_s = &live_sw,
            .standing_dead_par_umol_s = &live_sw,
        }),
    );
    try std.testing.expectEqual(@as(f64, 42), live_sw[0]);
}
const std = @import("std");
