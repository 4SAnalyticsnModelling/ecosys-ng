const std = @import("std");

pub const micropore_species_count = 50;
pub const macropore_species_count = 49;
const pre_phosphate_species_count = 33;

pub const Inputs = struct {
    micropore_water_m3: f64,
    macropore_water_m3: f64,
    minimum_water_m3: f64,
    maximum_transfer_volume_fraction: f64,
    soil_bulk_volume_m3: f64,
    substep_fraction: f64,
    nonband_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
    /// Micropore free PO4/H3PO4 are g P; remaining fields are mol.
    micropore_inventory_amount: []const f64,
    /// Compact 33 salt/complex plus 16 phosphate inventories, mol.
    macropore_inventory_mol: []const f64,
};

/// Exact compatibility translation of TRNSFRS.F lines 3379--3599.
/// Positive flux is from macropore toward micropore. H4SiO4 is absent from the
/// compact 49-field result because no macropore H4SiO4 inventory exists.
pub fn calculate(inputs: Inputs, flux_amount_per_step: []f64) !void {
    if (inputs.micropore_inventory_amount.len != micropore_species_count or
        inputs.macropore_inventory_mol.len != macropore_species_count or
        flux_amount_per_step.len != macropore_species_count)
        return error.SurfacePoreDiffusiveSoluteDimensionMismatch;
    try validate(inputs);
    if (!(inputs.macropore_water_m3 > inputs.minimum_water_m3)) {
        @memset(flux_amount_per_step, 0);
        return;
    }

    const transferable_macropore_water_m3 = @min(
        inputs.maximum_transfer_volume_fraction * inputs.soil_bulk_volume_m3,
        inputs.macropore_water_m3,
    );
    const combined_water_m3 = inputs.micropore_water_m3 + transferable_macropore_water_m3;
    if (!(combined_water_m3 > 0) or !std.math.isFinite(combined_water_m3))
        return error.InvalidSurfacePoreDiffusiveSoluteWaterVolume;
    for (0..macropore_species_count) |species| {
        const value = fluxForSpecies(inputs, species, transferable_macropore_water_m3, combined_water_m3);
        if (!std.math.isFinite(value)) return error.NonFiniteSurfacePoreDiffusiveSoluteResult;
    }
    for (0..macropore_species_count) |species|
        flux_amount_per_step[species] =
            fluxForSpecies(inputs, species, transferable_macropore_water_m3, combined_water_m3);
}

fn fluxForSpecies(inputs: Inputs, species: usize, transferable_water_m3: f64, combined_water_m3: f64) f64 {
    const micropore_species = if (species < pre_phosphate_species_count) species else species + 1;
    var flux = inputs.substep_fraction *
        (@max(0, inputs.macropore_inventory_mol[species]) * inputs.micropore_water_m3 -
            @max(0, inputs.micropore_inventory_amount[micropore_species]) * transferable_water_m3) /
        combined_water_m3;
    if (species >= 33 and species < 41) flux *= inputs.nonband_phosphate_fraction;
    if (species >= 41) flux *= inputs.band_phosphate_fraction;
    return flux;
}

fn validate(inputs: Inputs) !void {
    inline for (.{ inputs.micropore_water_m3, inputs.macropore_water_m3, inputs.minimum_water_m3, inputs.maximum_transfer_volume_fraction, inputs.soil_bulk_volume_m3, inputs.substep_fraction, inputs.nonband_phosphate_fraction, inputs.band_phosphate_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfacePoreDiffusiveSoluteInput;
    if (inputs.micropore_water_m3 < 0 or inputs.macropore_water_m3 < 0 or inputs.minimum_water_m3 < 0 or
        inputs.maximum_transfer_volume_fraction < 0 or inputs.soil_bulk_volume_m3 < 0 or
        inputs.substep_fraction < 0 or inputs.nonband_phosphate_fraction < 0 or
        inputs.nonband_phosphate_fraction > 1 or inputs.band_phosphate_fraction < 0 or
        inputs.band_phosphate_fraction > 1)
        return error.InvalidSurfacePoreDiffusiveSoluteInput;
    for (inputs.micropore_inventory_amount) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfacePoreDiffusiveSoluteInput;
    for (inputs.macropore_inventory_mol) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfacePoreDiffusiveSoluteInput;
}

fn fixture(micropore: []const f64, macropore: []const f64) Inputs {
    return .{ .micropore_water_m3 = 4, .macropore_water_m3 = 3, .minimum_water_m3 = 0.01, .maximum_transfer_volume_fraction = 0.5, .soil_bulk_volume_m3 = 4, .substep_fraction = 0.25, .nonband_phosphate_fraction = 0.2, .band_phosphate_fraction = 0.8, .micropore_inventory_amount = micropore, .macropore_inventory_mol = macropore };
}

test "TRNSFRS applies capped transferable volume and compact source mapping" {
    const micropore = [_]f64{2} ** micropore_species_count;
    const macropore = [_]f64{8} ** macropore_species_count;
    var flux = [_]f64{0} ** macropore_species_count;
    try calculate(fixture(&micropore, &macropore), &flux);
    // VOLWHS=min(0.5*4,3)=2; 0.25*(8*4-2*2)/(4+2)=7/6.
    try std.testing.expectApproxEqAbs(@as(f64, 7.0 / 6.0), flux[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 7.0 / 30.0), flux[33], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 14.0 / 15.0), flux[41], 1e-15);
}

test "negative concentration contrast produces micropore-to-macropore flux" {
    const micropore = [_]f64{10} ** micropore_species_count;
    const macropore = [_]f64{1} ** macropore_species_count;
    var flux = [_]f64{0} ** macropore_species_count;
    try calculate(fixture(&micropore, &macropore), &flux);
    try std.testing.expect(flux[0] < 0);
}

test "inactive macropore water explicitly zeroes all compact fluxes" {
    const micropore = [_]f64{2} ** micropore_species_count;
    const macropore = [_]f64{8} ** macropore_species_count;
    var flux = [_]f64{9} ** macropore_species_count;
    var inputs = fixture(&micropore, &macropore);
    inputs.macropore_water_m3 = 0;
    try calculate(inputs, &flux);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** macropore_species_count), &flux);
}

test "late invalid micropore inventory leaves output atomic" {
    var micropore = [_]f64{2} ** micropore_species_count;
    micropore[micropore_species_count - 1] = std.math.inf(f64);
    const macropore = [_]f64{8} ** macropore_species_count;
    var flux = [_]f64{9} ** macropore_species_count;
    try std.testing.expectError(error.NonFiniteSurfacePoreDiffusiveSoluteInput, calculate(fixture(&micropore, &macropore), &flux));
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** macropore_species_count), &flux);
}
