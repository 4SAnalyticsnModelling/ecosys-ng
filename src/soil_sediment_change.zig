const std = @import("std");

/// Publishes REDIST `TSEDER` from the authoritative extensive mineral owner.
/// Positive is deposition and negative is erosion.
pub fn publishAcceptedNetSedimentMg(hour_start_surface_soil_mass_Mg: []const f64, final_surface_soil_mass_Mg: []const f64, net_sediment_Mg: []f64) !void {
    if (hour_start_surface_soil_mass_Mg.len != final_surface_soil_mass_Mg.len or net_sediment_Mg.len != final_surface_soil_mass_Mg.len) return error.NetSedimentLedgerDimensionMismatch;
    for (hour_start_surface_soil_mass_Mg, final_surface_soil_mass_Mg) |initial_mass_Mg, final_mass_Mg| {
        const change_Mg = final_mass_Mg - initial_mass_Mg;
        if (!std.math.isFinite(initial_mass_Mg) or initial_mass_Mg <= 0 or !std.math.isFinite(final_mass_Mg) or final_mass_Mg <= 0 or !std.math.isFinite(change_Mg)) return error.InvalidNetSedimentLedgerState;
    }
    for (hour_start_surface_soil_mass_Mg, final_surface_soil_mass_Mg, net_sediment_Mg) |initial_mass_Mg, final_mass_Mg, *change_Mg| change_Mg.* = final_mass_Mg - initial_mass_Mg;
}

test "REDIST TSEDER retains deposition sign and publishes atomically" {
    var change = [_]f64{ 8, 9, 10 };
    try publishAcceptedNetSedimentMg(&.{ 10, 20, 30 }, &.{ 9, 22, 30 }, &change);
    try std.testing.expectEqualSlices(f64, &.{ -1, 2, 0 }, &change);
    change = .{ 8, 9, 10 };
    try std.testing.expectError(error.InvalidNetSedimentLedgerState, publishAcceptedNetSedimentMg(&.{ 10, 20, 30 }, &.{ 9, std.math.nan(f64), 30 }, &change));
    try std.testing.expectEqualSlices(f64, &.{ 8, 9, 10 }, &change);
}
