const std = @import("std");
const species_registry = @import("transport_species.zig");

pub const surface_species_count = species_registry.index(.band_phosphate);
pub const topsoil_species_count = species_registry.AqueousSpecies.count;

pub const State = struct {
    surface_litter_flux_mol_per_step: []f64,
    topsoil_flux_mol_per_step: []f64,
};

/// Exact source-order translation of TRNSFRS.F lines 536--627.
///
/// The enclosing snow-route branch is caller-owned. Surface litter has only
/// the 42 non-band species represented at L=0; top mineral soil additionally
/// clears all eight fertilizer-band phosphate species. Fluxes are mol step-1.
pub fn reset(state: State) !void {
    if (state.surface_litter_flux_mol_per_step.len != surface_species_count or
        state.topsoil_flux_mol_per_step.len != topsoil_species_count)
        return error.SnowRouteSoilSoluteFluxDimensionMismatch;

    for (state.surface_litter_flux_mol_per_step) |*flux_mol| flux_mol.* = 0;
    for (state.topsoil_flux_mol_per_step) |*flux_mol| flux_mol.* = 0;
}

test "TRNSFRS snow route clears surface and topsoil species in source order" {
    var surface = [_]f64{3} ** surface_species_count;
    var topsoil = [_]f64{5} ** topsoil_species_count;
    try reset(.{
        .surface_litter_flux_mol_per_step = &surface,
        .topsoil_flux_mol_per_step = &topsoil,
    });
    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{0} ** surface_species_count),
        &surface,
    );
    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{0} ** topsoil_species_count),
        &topsoil,
    );
}

test "TRNSFRS snow route reset includes every band phosphate species" {
    var surface = [_]f64{0} ** surface_species_count;
    var topsoil = [_]f64{0} ** topsoil_species_count;
    for (species_registry.index(.band_phosphate)..topsoil_species_count) |index|
        topsoil[index] = @floatFromInt(index + 1);
    try reset(.{
        .surface_litter_flux_mol_per_step = &surface,
        .topsoil_flux_mol_per_step = &topsoil,
    });
    for (species_registry.index(.band_phosphate)..topsoil_species_count) |index|
        try std.testing.expectEqual(@as(f64, 0), topsoil[index]);
}

test "snow route solute reset dimension failure is atomic" {
    var short_surface = [_]f64{3} ** (surface_species_count - 1);
    var topsoil = [_]f64{5} ** topsoil_species_count;
    try std.testing.expectError(
        error.SnowRouteSoilSoluteFluxDimensionMismatch,
        reset(.{
            .surface_litter_flux_mol_per_step = &short_surface,
            .topsoil_flux_mol_per_step = &topsoil,
        }),
    );
    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{3} ** (surface_species_count - 1)),
        &short_surface,
    );
    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{5} ** topsoil_species_count),
        &topsoil,
    );
}
