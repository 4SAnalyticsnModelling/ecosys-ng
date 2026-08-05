const std = @import("std");
const species_registry = @import("solute_transport_species.zig");

pub const irrigation_concentration_count = species_registry.index(.band_phosphate) - 1;
pub const flux_species_count = species_registry.AqueousSpecies.count - 1;

/// Exact source-order translation of TRNSFRS.F lines 1351--1399.
///
/// Water is m3 step-1, concentrations are mol m-3, and output is mol step-1.
/// H4SiO4 has no source statement and therefore no array slot. Each band-
/// phosphate output reuses its corresponding non-band irrigation concentration.
pub fn calculate(
    subsurface_irrigation_water_m3_per_step: f64,
    irrigation_concentration_mol_per_m3: []const f64,
    non_band_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
    output_flux_mol_per_step: []f64,
) !void {
    if (irrigation_concentration_mol_per_m3.len != irrigation_concentration_count or
        output_flux_mol_per_step.len != flux_species_count)
        return error.SubsurfaceIrrigationSoluteDimensionMismatch;
    inline for (.{ subsurface_irrigation_water_m3_per_step, non_band_phosphate_fraction, band_phosphate_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSubsurfaceIrrigationInput;
    inline for (.{ non_band_phosphate_fraction, band_phosphate_fraction }) |fraction|
        if (fraction < 0 or fraction > 1) return error.InvalidPhosphateZoneFraction;

    for (0..flux_species_count) |species| {
        const concentration_index = concentrationIndex(species);
        const concentration = irrigation_concentration_mol_per_m3[concentration_index];
        if (!std.math.isFinite(concentration)) return error.NonFiniteSubsurfaceIrrigationInput;
        const flux = sourceFlux(species, subsurface_irrigation_water_m3_per_step, concentration, non_band_phosphate_fraction, band_phosphate_fraction);
        if (!std.math.isFinite(flux)) return error.NonFiniteSubsurfaceIrrigationFlux;
    }
    for (0..flux_species_count) |species| {
        const concentration = irrigation_concentration_mol_per_m3[concentrationIndex(species)];
        output_flux_mol_per_step[species] = sourceFlux(species, subsurface_irrigation_water_m3_per_step, concentration, non_band_phosphate_fraction, band_phosphate_fraction);
    }
}

fn concentrationIndex(species: usize) usize {
    if (species < irrigation_concentration_count) return species;
    return species - irrigation_concentration_count + species_registry.index(.non_band_phosphate) - 1;
}

fn sourceFlux(species: usize, water: f64, concentration: f64, non_band_fraction: f64, band_fraction: f64) f64 {
    const base = water * concentration;
    if (species >= irrigation_concentration_count) return base * band_fraction;
    if (species >= species_registry.index(.non_band_phosphate) - 1) return base * non_band_fraction;
    return base;
}

test "TRNSFRS subsurface irrigation preserves species and zone source order" {
    var concentrations: [irrigation_concentration_count]f64 = undefined;
    var flux = [_]f64{-1} ** flux_species_count;
    for (&concentrations, 0..) |*value, index| value.* = @floatFromInt(index + 1);
    try calculate(2, &concentrations, 0.7, 0.3, &flux);
    try std.testing.expectEqual(@as(f64, 2), flux[0]);
    const phosphate_index = species_registry.index(.non_band_phosphate) - 1;
    const phosphate_concentration = concentrations[phosphate_index];
    try std.testing.expectEqual(2 * phosphate_concentration * 0.7, flux[phosphate_index]);
    try std.testing.expectEqual(2 * phosphate_concentration * 0.3, flux[irrigation_concentration_count]);
}

test "subsurface irrigation late failure is atomic" {
    var concentrations = [_]f64{2} ** irrigation_concentration_count;
    var flux = [_]f64{9} ** flux_species_count;
    concentrations[irrigation_concentration_count - 1] = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteSubsurfaceIrrigationInput,
        calculate(1, &concentrations, 0.7, 0.3, &flux),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** flux_species_count), &flux);
}

test "subsurface irrigation topology is exact" {
    const concentrations = [_]f64{2} ** (irrigation_concentration_count - 1);
    var flux = [_]f64{9} ** flux_species_count;
    try std.testing.expectError(
        error.SubsurfaceIrrigationSoluteDimensionMismatch,
        calculate(1, &concentrations, 0.7, 0.3, &flux),
    );
}
