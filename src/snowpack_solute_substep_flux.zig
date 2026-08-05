const std = @import("std");
const species_registry = @import("solute_transport_species.zig");

pub const snow_species_count = species_registry.index(.band_phosphate) - 1;

/// Exact source-order translation of TRNSFRS.F lines 1050--1090.
///
/// Whole-step and substep fluxes are mol step-1. The runtime fraction is the
/// legacy XNPH. The dedicated snow axis has 41 actual fields: TRNSFRS has no
/// snow boundary statement for H4SiO4 and no placeholder slot is synthesized.
pub fn publish(
    whole_step_flux_mol: []const f64,
    substep_flux_mol: []f64,
    substep_fraction: f64,
) !void {
    if (whole_step_flux_mol.len != snow_species_count or
        substep_flux_mol.len != snow_species_count)
        return error.SnowpackSoluteSubstepDimensionMismatch;
    if (!std.math.isFinite(substep_fraction) or
        substep_fraction < 0 or substep_fraction > 1)
        return error.InvalidSoluteSubstepFraction;

    for (whole_step_flux_mol) |flux_mol| {
        if (!std.math.isFinite(flux_mol))
            return error.NonFiniteSnowpackSoluteFlux;
        const scaled = flux_mol * substep_fraction;
        if (!std.math.isFinite(scaled))
            return error.NonFiniteSnowpackSoluteFlux;
    }
    for (whole_step_flux_mol, substep_flux_mol) |flux_mol, *scaled_mol| {
        scaled_mol.* = flux_mol * substep_fraction;
    }
}

test "TRNSFRS snowpack boundary scales 41 source species in order" {
    var whole_step: [snow_species_count]f64 = undefined;
    var substep = [_]f64{-7} ** snow_species_count;
    for (&whole_step, 0..) |*flux_mol, index| flux_mol.* = @floatFromInt(index + 1);
    try publish(&whole_step, &substep, 0.25);
    try std.testing.expectEqual(@as(f64, 0.25), substep[species_registry.index(.aluminum)]);
    try std.testing.expectEqual(@as(f64, 10.25), substep[snow_species_count - 1]);
}

test "TRNSFRS snowpack boundary has no hydrogen silicate placeholder" {
    var whole_step = [_]f64{2} ** snow_species_count;
    var substep = [_]f64{5} ** snow_species_count;
    try publish(&whole_step, &substep, 0.5);
    try std.testing.expectEqual(@as(usize, 41), snow_species_count);
    try std.testing.expectEqual(@as(f64, 1), substep[snow_species_count - 1]);
}

test "snowpack substep publication failure is atomic" {
    var whole_step = [_]f64{2} ** snow_species_count;
    var substep = [_]f64{5} ** snow_species_count;
    whole_step[20] = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteSnowpackSoluteFlux,
        publish(&whole_step, &substep, 0.5),
    );
    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{5} ** snow_species_count),
        &substep,
    );
}

test "snowpack substep fraction is runtime validated" {
    const whole_step = [_]f64{2} ** snow_species_count;
    var substep = [_]f64{5} ** snow_species_count;
    try std.testing.expectError(
        error.InvalidSoluteSubstepFraction,
        publish(&whole_step, &substep, 1.1),
    );
}
