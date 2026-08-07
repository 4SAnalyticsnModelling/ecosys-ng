const std = @import("std");
const species_registry = @import("transport_species.zig");

pub const topsoil_species_count = species_registry.AqueousSpecies.count;

/// Exact source-order translation of TRNSFRS.F lines 1133--1182.
///
/// Whole-step and substep topsoil boundary fluxes are extensive mol step-1.
/// The topology includes 42 non-band species followed by the eight fertilizer-
/// band phosphate species, exactly matching the source statement sequence.
pub fn publish(
    whole_step_flux_mol: []const f64,
    substep_flux_mol: []f64,
    substep_fraction: f64,
) !void {
    if (whole_step_flux_mol.len != topsoil_species_count or
        substep_flux_mol.len != topsoil_species_count)
        return error.TopsoilSoluteSubstepDimensionMismatch;
    if (!std.math.isFinite(substep_fraction) or
        substep_fraction < 0 or substep_fraction > 1)
        return error.InvalidSoluteSubstepFraction;

    for (whole_step_flux_mol) |flux_mol| {
        if (!std.math.isFinite(flux_mol) or
            !std.math.isFinite(flux_mol * substep_fraction))
            return error.NonFiniteTopsoilSoluteFlux;
    }
    for (whole_step_flux_mol, substep_flux_mol) |flux_mol, *scaled_mol|
        scaled_mol.* = flux_mol * substep_fraction;
}

test "TRNSFRS topsoil scales all 50 source species in order" {
    var whole_step: [topsoil_species_count]f64 = undefined;
    var substep = [_]f64{-1} ** topsoil_species_count;
    for (&whole_step, 0..) |*flux_mol, index| flux_mol.* = @floatFromInt(index + 1);
    try publish(&whole_step, &substep, 0.25);
    for (whole_step, substep) |whole_mol, scaled_mol|
        try std.testing.expectEqual(whole_mol * 0.25, scaled_mol);
}

test "TRNSFRS topsoil includes hydrogen silicate and band phosphate" {
    var whole_step = [_]f64{0} ** topsoil_species_count;
    var substep = [_]f64{-1} ** topsoil_species_count;
    whole_step[species_registry.index(.hydrogen_silicate)] = 8;
    whole_step[species_registry.index(.band_phosphate)] = 10;
    whole_step[species_registry.index(.band_magnesium_hpo4)] = 12;
    try publish(&whole_step, &substep, 0.5);
    try std.testing.expectEqual(@as(f64, 4), substep[species_registry.index(.hydrogen_silicate)]);
    try std.testing.expectEqual(@as(f64, 5), substep[species_registry.index(.band_phosphate)]);
    try std.testing.expectEqual(@as(f64, 6), substep[species_registry.index(.band_magnesium_hpo4)]);
}

test "topsoil substep publication failure is atomic" {
    var whole_step = [_]f64{2} ** topsoil_species_count;
    var substep = [_]f64{5} ** topsoil_species_count;
    whole_step[topsoil_species_count - 1] = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteTopsoilSoluteFlux,
        publish(&whole_step, &substep, 0.5),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{5} ** topsoil_species_count), &substep);
}

test "topsoil substep dimensions are exact" {
    const whole_step = [_]f64{2} ** (topsoil_species_count - 1);
    var substep = [_]f64{5} ** topsoil_species_count;
    try std.testing.expectError(
        error.TopsoilSoluteSubstepDimensionMismatch,
        publish(&whole_step, &substep, 0.5),
    );
}
