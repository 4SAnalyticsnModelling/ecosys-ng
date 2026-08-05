const std = @import("std");

pub const snow_drift_species_count = 41;
pub const BoundaryAxis = enum { east_west, north_south };

pub const BoundaryFluxes = struct {
    axis: BoundaryAxis,
    /// Legacy side 2; always present for the selected internal direction.
    forward_side_flux_mol_per_step: []f64,
    /// Legacy side 1; null at a physical backward edge.
    backward_side_flux_mol_per_step: ?[]f64,
};

/// Exact inactive-snow-drift branch from TRNSFRS.F lines 4754--4838.
/// Instantaneous compact 41-field fluxes are zeroed; accumulated REDIST state
/// is not modified. All present dimensions validate before either mutation.
pub fn reset(fluxes: BoundaryFluxes) !void {
    _ = fluxes.axis;
    if (fluxes.forward_side_flux_mol_per_step.len != snow_drift_species_count)
        return error.InactiveSnowDriftBoundaryFluxDimensionMismatch;
    if (fluxes.backward_side_flux_mol_per_step) |backward|
        if (backward.len != snow_drift_species_count)
            return error.InactiveSnowDriftBoundaryFluxDimensionMismatch;
    @memset(fluxes.forward_side_flux_mol_per_step, 0);
    if (fluxes.backward_side_flux_mol_per_step) |backward| @memset(backward, 0);
}

test "TRNSFRS inactive snow drift clears both valid compact sides" {
    var forward = [_]f64{9} ** snow_drift_species_count;
    var backward = [_]f64{7} ** snow_drift_species_count;
    try reset(.{ .axis = .east_west, .forward_side_flux_mol_per_step = &forward, .backward_side_flux_mol_per_step = &backward });
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** snow_drift_species_count), &forward);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** snow_drift_species_count), &backward);
}

test "physical backward edge clears only forward compact side" {
    var forward = [_]f64{9} ** snow_drift_species_count;
    try reset(.{ .axis = .north_south, .forward_side_flux_mol_per_step = &forward, .backward_side_flux_mol_per_step = null });
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** snow_drift_species_count), &forward);
}

test "late backward dimension failure leaves forward side atomic" {
    var forward = [_]f64{9} ** snow_drift_species_count;
    var short_backward = [_]f64{7} ** (snow_drift_species_count - 1);
    try std.testing.expectError(
        error.InactiveSnowDriftBoundaryFluxDimensionMismatch,
        reset(.{ .axis = .east_west, .forward_side_flux_mol_per_step = &forward, .backward_side_flux_mol_per_step = &short_backward }),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** snow_drift_species_count), &forward);
}

test "inactive snow reset rejects phantom H4SiO4 slot" {
    var long_forward = [_]f64{9} ** (snow_drift_species_count + 1);
    try std.testing.expectError(
        error.InactiveSnowDriftBoundaryFluxDimensionMismatch,
        reset(.{ .axis = .north_south, .forward_side_flux_mol_per_step = &long_forward, .backward_side_flux_mol_per_step = null }),
    );
}
