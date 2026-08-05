const std = @import("std");
const species_registry = @import("solute_transport_species.zig");

pub const surface_species_count = species_registry.index(.band_phosphate);

/// Exact source-order translation of TRNSFRS.F lines 1091--1132.
///
/// Whole-step and substep boundary fluxes are extensive mol step-1. Unlike
/// the preceding snowpack block, this surface-litter block explicitly
/// includes H4SiO4 at source line 1124.
pub fn publish(
    whole_step_flux_mol: []const f64,
    substep_flux_mol: []f64,
    substep_fraction: f64,
) !void {
    if (whole_step_flux_mol.len != surface_species_count or
        substep_flux_mol.len != surface_species_count)
        return error.SurfaceLitterSoluteSubstepDimensionMismatch;
    if (!std.math.isFinite(substep_fraction) or
        substep_fraction < 0 or substep_fraction > 1)
        return error.InvalidSoluteSubstepFraction;

    for (whole_step_flux_mol) |flux_mol| {
        if (!std.math.isFinite(flux_mol) or
            !std.math.isFinite(flux_mol * substep_fraction))
            return error.NonFiniteSurfaceLitterSoluteFlux;
    }
    for (whole_step_flux_mol, substep_flux_mol) |flux_mol, *scaled_mol|
        scaled_mol.* = flux_mol * substep_fraction;
}

test "TRNSFRS surface litter scales all 42 source species in order" {
    var whole_step: [surface_species_count]f64 = undefined;
    var substep = [_]f64{-1} ** surface_species_count;
    for (&whole_step, 0..) |*flux_mol, index| flux_mol.* = @floatFromInt(index + 1);
    try publish(&whole_step, &substep, 0.25);
    for (whole_step, substep) |whole_mol, scaled_mol|
        try std.testing.expectEqual(whole_mol * 0.25, scaled_mol);
}

test "TRNSFRS surface litter explicitly publishes hydrogen silicate" {
    var whole_step = [_]f64{0} ** surface_species_count;
    var substep = [_]f64{-1} ** surface_species_count;
    whole_step[species_registry.index(.hydrogen_silicate)] = 8;
    try publish(&whole_step, &substep, 0.5);
    try std.testing.expectEqual(
        @as(f64, 4),
        substep[species_registry.index(.hydrogen_silicate)],
    );
}

test "surface litter substep publication failure is atomic" {
    var whole_step = [_]f64{2} ** surface_species_count;
    var substep = [_]f64{5} ** surface_species_count;
    whole_step[surface_species_count - 1] = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterSoluteFlux,
        publish(&whole_step, &substep, 0.5),
    );
    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{5} ** surface_species_count),
        &substep,
    );
}

test "surface litter substep dimensions are exact" {
    const whole_step = [_]f64{2} ** (surface_species_count - 1);
    var substep = [_]f64{5} ** surface_species_count;
    try std.testing.expectError(
        error.SurfaceLitterSoluteSubstepDimensionMismatch,
        publish(&whole_step, &substep, 0.5),
    );
}
