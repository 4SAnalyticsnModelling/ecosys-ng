const std = @import("std");

pub const runoff_species_count = 42;

pub const BoundaryAxis = enum {
    east_west,
    north_south,
};

pub const BoundaryFluxes = struct {
    axis: BoundaryAxis,
    /// Legacy side 2 at the forward neighbor; always present for this block.
    forward_side_flux_amount_per_step: []f64,
    /// Legacy side 1 at the backward neighbor; absent at a physical edge.
    backward_side_flux_amount_per_step: ?[]f64,
};

/// Exact source-order inactive-runoff branch from TRNSFRS.F lines 4368--4454.
/// Only instantaneous 42-field boundary fluxes are zeroed; REDIST accumulated
/// fluxes are deliberately untouched. Both present dimensions validate first.
pub fn reset(fluxes: BoundaryFluxes) !void {
    _ = fluxes.axis;
    if (fluxes.forward_side_flux_amount_per_step.len != runoff_species_count)
        return error.InactiveRunoffBoundaryFluxDimensionMismatch;
    if (fluxes.backward_side_flux_amount_per_step) |backward|
        if (backward.len != runoff_species_count)
            return error.InactiveRunoffBoundaryFluxDimensionMismatch;

    @memset(fluxes.forward_side_flux_amount_per_step, 0);
    if (fluxes.backward_side_flux_amount_per_step) |backward| @memset(backward, 0);
}

test "TRNSFRS inactive runoff clears both valid boundary sides" {
    var forward = [_]f64{9} ** runoff_species_count;
    var backward = [_]f64{7} ** runoff_species_count;
    try reset(.{ .axis = .east_west, .forward_side_flux_amount_per_step = &forward, .backward_side_flux_amount_per_step = &backward });
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** runoff_species_count), &forward);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** runoff_species_count), &backward);
}

test "physical backward edge clears only the always-present forward side" {
    var forward = [_]f64{9} ** runoff_species_count;
    try reset(.{ .axis = .north_south, .forward_side_flux_amount_per_step = &forward, .backward_side_flux_amount_per_step = null });
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** runoff_species_count), &forward);
}

test "late backward dimension failure leaves forward side atomic" {
    var forward = [_]f64{9} ** runoff_species_count;
    var short_backward = [_]f64{7} ** (runoff_species_count - 1);
    try std.testing.expectError(
        error.InactiveRunoffBoundaryFluxDimensionMismatch,
        reset(.{ .axis = .east_west, .forward_side_flux_amount_per_step = &forward, .backward_side_flux_amount_per_step = &short_backward }),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** runoff_species_count), &forward);
    try std.testing.expectEqualSlices(f64, &([_]f64{7} ** (runoff_species_count - 1)), &short_backward);
}

test "long forward topology is rejected" {
    var long_forward = [_]f64{9} ** (runoff_species_count + 1);
    try std.testing.expectError(
        error.InactiveRunoffBoundaryFluxDimensionMismatch,
        reset(.{ .axis = .north_south, .forward_side_flux_amount_per_step = &long_forward, .backward_side_flux_amount_per_step = null }),
    );
}
