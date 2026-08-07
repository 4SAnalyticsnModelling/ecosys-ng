const std = @import("std");
const species_registry = @import("transport_species.zig");

pub const snowpack_species_count = species_registry.index(.band_phosphate) - 1;
pub const non_band_soil_species_count = species_registry.index(.band_phosphate);
pub const topsoil_species_count = species_registry.AqueousSpecies.count;

pub const State = struct {
    snowpack_flux_mol_per_step: []f64,
    surface_litter_flux_mol_per_step: []f64,
    topsoil_flux_mol_per_step: []f64,
};

/// Exact source-order translation of TRNSFRS.F lines 893--1025.
///
/// The caller owns the source `ELSE` admission at line 892. Snowpack and
/// snowpack contains 41 species because H4SiO4 is absent. Surface litter has
/// 42 non-band species; topsoil also has eight band-phosphate species. All
/// values are extensive mol step-1.
pub fn reset(state: State) !void {
    if (state.snowpack_flux_mol_per_step.len != snowpack_species_count or
        state.surface_litter_flux_mol_per_step.len != non_band_soil_species_count or
        state.topsoil_flux_mol_per_step.len != topsoil_species_count)
        return error.InactiveAtmosphericSoluteFluxDimensionMismatch;

    for (state.snowpack_flux_mol_per_step) |*flux_mol| flux_mol.* = 0;
    for (state.surface_litter_flux_mol_per_step) |*flux_mol| flux_mol.* = 0;
    for (state.topsoil_flux_mol_per_step) |*flux_mol| flux_mol.* = 0;
}

test "TRNSFRS inactive atmosphere clears all three destinations in source order" {
    var snow = [_]f64{2} ** snowpack_species_count;
    var surface = [_]f64{3} ** non_band_soil_species_count;
    var topsoil = [_]f64{5} ** topsoil_species_count;
    try reset(.{
        .snowpack_flux_mol_per_step = &snow,
        .surface_litter_flux_mol_per_step = &surface,
        .topsoil_flux_mol_per_step = &topsoil,
    });
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** snowpack_species_count), &snow);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** non_band_soil_species_count), &surface);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** topsoil_species_count), &topsoil);
}

test "TRNSFRS inactive atmosphere clears topsoil band phosphate" {
    var snow = [_]f64{0} ** snowpack_species_count;
    var surface = [_]f64{0} ** non_band_soil_species_count;
    var topsoil = [_]f64{0} ** topsoil_species_count;
    for (species_registry.index(.band_phosphate)..topsoil_species_count) |index|
        topsoil[index] = @floatFromInt(index + 1);
    try reset(.{
        .snowpack_flux_mol_per_step = &snow,
        .surface_litter_flux_mol_per_step = &surface,
        .topsoil_flux_mol_per_step = &topsoil,
    });
    for (species_registry.index(.band_phosphate)..topsoil_species_count) |index|
        try std.testing.expectEqual(@as(f64, 0), topsoil[index]);
}

test "inactive atmospheric flux dimension failure is atomic" {
    var snow = [_]f64{2} ** snowpack_species_count;
    var surface = [_]f64{3} ** (non_band_soil_species_count - 1);
    var topsoil = [_]f64{5} ** topsoil_species_count;
    try std.testing.expectError(
        error.InactiveAtmosphericSoluteFluxDimensionMismatch,
        reset(.{
            .snowpack_flux_mol_per_step = &snow,
            .surface_litter_flux_mol_per_step = &surface,
            .topsoil_flux_mol_per_step = &topsoil,
        }),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{2} ** snowpack_species_count), &snow);
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** (non_band_soil_species_count - 1)), &surface);
    try std.testing.expectEqualSlices(f64, &([_]f64{5} ** topsoil_species_count), &topsoil);
}
