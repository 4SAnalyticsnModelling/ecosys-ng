const std = @import("std");

pub const species_count = 50;
pub const Direction = enum { east_west, north_south, vertical };
pub const CellConnectivity = enum { interconnected, standalone };

/// Compatibility translation of TRNSFRS.F lines 8006--8057.
/// This ELSE belongs to `NCN(M2,M1) != 3 OR N == 3`: only a standalone
/// horizontal boundary clears the instantaneous 50-species micropore flux.
pub fn resetIfRequired(
    connectivity: CellConnectivity,
    direction: Direction,
    boundary_flux_amount_per_step: []f64,
) !bool {
    if (connectivity != .standalone or direction == .vertical) return false;
    if (boundary_flux_amount_per_step.len != species_count)
        return error.StandaloneExternalMicroporeFluxDimensionMismatch;
    for (boundary_flux_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteStandaloneExternalMicroporeFluxInput;
    @memset(boundary_flux_amount_per_step, 0);
    return true;
}

test "TRNSFRS standalone east-west boundary clears all 50 instantaneous species" {
    var flux = [_]f64{3} ** species_count;
    try std.testing.expect(try resetIfRequired(.standalone, .east_west, &flux));
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** species_count), &flux);
}

test "TRNSFRS standalone north-south boundary includes H4SiO4 and both P pools" {
    var flux = [_]f64{3} ** species_count;
    try std.testing.expect(try resetIfRequired(.standalone, .north_south, &flux));
    try std.testing.expectEqual(@as(f64, 0), flux[33]);
    try std.testing.expectEqual(@as(f64, 0), flux[34]);
    try std.testing.expectEqual(@as(f64, 0), flux[49]);
}

test "TRNSFRS vertical exception leaves standalone flux untouched" {
    var flux = [_]f64{3} ** species_count;
    try std.testing.expect(!try resetIfRequired(.standalone, .vertical, &flux));
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** species_count), &flux);
}

test "TRNSFRS interconnected horizontal boundary leaves flux untouched" {
    var flux = [_]f64{3} ** species_count;
    try std.testing.expect(!try resetIfRequired(.interconnected, .east_west, &flux));
    try std.testing.expectEqual(@as(f64, 3), flux[49]);
}

test "active reset validates atomically before clearing" {
    var flux = [_]f64{3} ** species_count;
    flux[49] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteStandaloneExternalMicroporeFluxInput, resetIfRequired(.standalone, .north_south, &flux));
    try std.testing.expectEqual(@as(f64, 3), flux[0]);
}
