const std = @import("std");
const species_registry = @import("solute_transport_species.zig");

/// TRNSFRS snow-input fields omit H4SiO4, which is present in the contiguous
/// soil aqueous-species registry immediately before the phosphate fields.
pub const snow_input_species_count = species_registry.index(.band_phosphate) - 1;

/// Exact source-order translation of TRNSFRS.F lines 666--706.
///
/// The caller owns the lines 633--634 bare-soil precipitation admission.
/// Snowpack inputs are extensive mol step-1. H4SiO4 and band phosphate species
/// are not represented in this source array and therefore are excluded.
pub fn reset(snowpack_input_flux_mol_per_step: []f64) !void {
    if (snowpack_input_flux_mol_per_step.len != snow_input_species_count)
        return error.BareSoilSnowSoluteFluxDimensionMismatch;
    for (snowpack_input_flux_mol_per_step) |*flux_mol| flux_mol.* = 0;
}

test "TRNSFRS bare soil branch clears every snow input in source order" {
    var flux = [_]f64{4} ** snow_input_species_count;
    try reset(&flux);
    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{0} ** snow_input_species_count),
        &flux,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        flux[snow_input_species_count - 1],
    );
}

test "bare soil snow reset rejects short snow input axis atomically" {
    var flux = [_]f64{4} ** (snow_input_species_count - 1);
    try std.testing.expectError(
        error.BareSoilSnowSoluteFluxDimensionMismatch,
        reset(&flux),
    );
    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{4} ** (snow_input_species_count - 1)),
        &flux,
    );
}

test "snow input axis omits silicate and band phosphate" {
    try std.testing.expectEqual(@as(usize, 41), snow_input_species_count);
    try std.testing.expectEqual(
        snow_input_species_count + 1,
        species_registry.index(.band_phosphate),
    );
}
