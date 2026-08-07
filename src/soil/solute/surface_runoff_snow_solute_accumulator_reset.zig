const std = @import("std");

/// Runoff has 33 salt/complex species, H4SiO4, and eight phosphate species.
pub const runoff_species_count = 42;

/// Snow redistribution omits H4SiO4: 33 salt/complex species and eight
/// phosphate species are stored contiguously in source order.
pub const snow_redistribution_species_count = 41;

pub const Accumulators = struct {
    /// Net surface-runoff solute fluxes, mol step-1.
    runoff_flux_mol_per_step: []f64,
    /// Net snow-redistribution solute fluxes, mol step-1.
    snow_redistribution_flux_mol_per_step: []f64,
};

/// Exact source-order translation of TRNSFRS.F lines 1626--1708.
///
/// The enclosing sub-hour, cell, and dynamic-salt loops remain caller-owned.
/// Both dimensions are checked before either source group is mutated.
pub fn reset(accumulators: Accumulators) !void {
    if (accumulators.runoff_flux_mol_per_step.len != runoff_species_count)
        return error.RunoffSoluteAccumulatorDimensionMismatch;
    if (accumulators.snow_redistribution_flux_mol_per_step.len != snow_redistribution_species_count)
        return error.SnowRedistributionSoluteAccumulatorDimensionMismatch;

    @memset(accumulators.runoff_flux_mol_per_step, 0);
    @memset(accumulators.snow_redistribution_flux_mol_per_step, 0);
}

test "TRNSFRS resets runoff then compact snow redistribution accumulators" {
    var runoff = [_]f64{3} ** runoff_species_count;
    var snow = [_]f64{5} ** snow_redistribution_species_count;

    try reset(.{
        .runoff_flux_mol_per_step = &runoff,
        .snow_redistribution_flux_mol_per_step = &snow,
    });

    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** runoff_species_count), &runoff);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** snow_redistribution_species_count), &snow);
}

test "TRNSFRS runoff retains the H4SiO4 slot omitted by snow redistribution" {
    try std.testing.expectEqual(@as(usize, 42), runoff_species_count);
    try std.testing.expectEqual(@as(usize, 41), snow_redistribution_species_count);
}

test "dimension validation is atomic across both accumulator groups" {
    var runoff = [_]f64{3} ** runoff_species_count;
    var short_snow = [_]f64{5} ** (snow_redistribution_species_count - 1);

    try std.testing.expectError(
        error.SnowRedistributionSoluteAccumulatorDimensionMismatch,
        reset(.{
            .runoff_flux_mol_per_step = &runoff,
            .snow_redistribution_flux_mol_per_step = &short_snow,
        }),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** runoff_species_count), &runoff);
    try std.testing.expectEqualSlices(f64, &([_]f64{5} ** (snow_redistribution_species_count - 1)), &short_snow);
}

test "long runoff topology is rejected before snow mutation" {
    var long_runoff = [_]f64{3} ** (runoff_species_count + 1);
    var snow = [_]f64{5} ** snow_redistribution_species_count;

    try std.testing.expectError(
        error.RunoffSoluteAccumulatorDimensionMismatch,
        reset(.{
            .runoff_flux_mol_per_step = &long_runoff,
            .snow_redistribution_flux_mol_per_step = &snow,
        }),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{5} ** snow_redistribution_species_count), &snow);
}
