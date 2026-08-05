const std = @import("std");
const species_registry = @import("solute_transport_species.zig");

pub const species_count = species_registry.AqueousSpecies.count - 1;

/// Exact source-order translation of TRNSFRS.F lines 1420--1468.
///
/// Whole-step and substep irrigation fluxes are extensive mol step-1. The
/// complete 49-field non-band/band source axis is traversed in source order.
pub fn publish(
    whole_step_flux_mol: []const f64,
    substep_flux_mol: []f64,
    substep_fraction: f64,
) !void {
    if (whole_step_flux_mol.len != species_count or substep_flux_mol.len != species_count)
        return error.SubsurfaceIrrigationSubstepDimensionMismatch;
    if (!std.math.isFinite(substep_fraction) or substep_fraction < 0 or substep_fraction > 1)
        return error.InvalidSoluteSubstepFraction;

    for (whole_step_flux_mol) |flux_mol| {
        if (!std.math.isFinite(flux_mol) or !std.math.isFinite(flux_mol * substep_fraction))
            return error.NonFiniteSubsurfaceIrrigationFlux;
    }
    for (whole_step_flux_mol, substep_flux_mol) |flux_mol, *scaled_mol| {
        scaled_mol.* = flux_mol * substep_fraction;
    }
}

test "TRNSFRS subsurface irrigation scales 49 source fluxes in order" {
    var whole_step: [species_count]f64 = undefined;
    var substep = [_]f64{-1} ** species_count;
    for (&whole_step, 0..) |*flux_mol, index| flux_mol.* = @floatFromInt(index + 1);
    try publish(&whole_step, &substep, 0.25);
    try std.testing.expectEqual(@as(f64, 0.25), substep[0]);
    try std.testing.expectEqual(
        whole_step[species_count - 1] * 0.25,
        substep[species_count - 1],
    );
}

test "subsurface irrigation substep late failure is atomic" {
    var whole_step = [_]f64{2} ** species_count;
    var substep = [_]f64{7} ** species_count;
    whole_step[species_count - 1] = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteSubsurfaceIrrigationFlux,
        publish(&whole_step, &substep, 0.5),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{7} ** species_count), &substep);
}

test "subsurface irrigation substep topology is exact" {
    const whole_step = [_]f64{2} ** (species_count - 1);
    var substep = [_]f64{7} ** species_count;
    try std.testing.expectError(
        error.SubsurfaceIrrigationSubstepDimensionMismatch,
        publish(&whole_step, &substep, 0.5),
    );
}
