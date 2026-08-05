const std = @import("std");

pub const micropore_species_count = 50;
pub const macropore_species_count = 49;
const pre_phosphate_without_silicate_count = 33;

pub const Inputs = struct {
    /// Positive is macropore to micropore, negative is the reverse, m3 step-1.
    macropore_to_micropore_water_flux_m3_per_step: f64,
    micropore_water_m3: f64,
    macropore_water_m3: f64,
    minimum_water_m3: f64,
    maximum_convective_fraction: f64,
    nonband_phosphate_fraction: f64,
    band_phosphate_fraction: f64,
    /// Micropore free PO4/H3PO4 are g P; remaining fields are mol.
    micropore_inventory_amount: []const f64,
    /// 33 salt/complex (no H4SiO4) plus 16 phosphate inventories, mol.
    macropore_inventory_mol: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 3149--3350.
///
/// The legacy positive/zero branches leave `RFLHYS` stale because macropores
/// have no H4SiO4 field. ecosys-ng explicitly publishes zero for that unsafe
/// missing-source case while preserving every scientific routing equation.
pub fn calculate(inputs: Inputs, flux_amount_per_step: []f64) !void {
    if (inputs.micropore_inventory_amount.len != micropore_species_count or
        inputs.macropore_inventory_mol.len != macropore_species_count or
        flux_amount_per_step.len != micropore_species_count)
        return error.SurfacePoreConvectiveSoluteDimensionMismatch;
    try validate(inputs);

    const water_flux = inputs.macropore_to_micropore_water_flux_m3_per_step;
    const fraction = if (water_flux > 0)
        if (inputs.macropore_water_m3 > inputs.minimum_water_m3)
            std.math.clamp(water_flux / inputs.macropore_water_m3, 0, inputs.maximum_convective_fraction)
        else
            inputs.maximum_convective_fraction
    else if (water_flux < 0)
        if (inputs.micropore_water_m3 > inputs.minimum_water_m3)
            std.math.clamp(water_flux / inputs.micropore_water_m3, -inputs.maximum_convective_fraction, 0)
        else
            -inputs.maximum_convective_fraction
    else
        0;

    for (0..micropore_species_count) |species| {
        const value = sourceFlux(inputs, species, fraction, water_flux > 0);
        if (!std.math.isFinite(value)) return error.NonFiniteSurfacePoreConvectiveSoluteResult;
    }
    for (0..micropore_species_count) |species|
        flux_amount_per_step[species] = sourceFlux(inputs, species, fraction, water_flux > 0);
}

fn sourceFlux(inputs: Inputs, species: usize, fraction: f64, macropore_donor: bool) f64 {
    if (!macropore_donor) {
        const partition = phosphatePartition(inputs, species, fraction);
        return @max(0, inputs.micropore_inventory_amount[species]) * partition;
    }
    if (species == 33) return 0;
    const macropore_species = if (species < 33) species else species - 1;
    const partition = phosphatePartition(inputs, species, fraction);
    return @max(0, inputs.macropore_inventory_mol[macropore_species]) * partition;
}

fn phosphatePartition(inputs: Inputs, species: usize, fraction: f64) f64 {
    if (species >= 34 and species < 42) return fraction * inputs.nonband_phosphate_fraction;
    if (species >= 42) return fraction * inputs.band_phosphate_fraction;
    return fraction;
}

fn validate(inputs: Inputs) !void {
    inline for (.{ inputs.macropore_to_micropore_water_flux_m3_per_step, inputs.micropore_water_m3, inputs.macropore_water_m3, inputs.minimum_water_m3, inputs.maximum_convective_fraction, inputs.nonband_phosphate_fraction, inputs.band_phosphate_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfacePoreConvectiveSoluteInput;
    if (inputs.micropore_water_m3 < 0 or inputs.macropore_water_m3 < 0 or inputs.minimum_water_m3 < 0 or
        inputs.maximum_convective_fraction < 0 or inputs.nonband_phosphate_fraction < 0 or
        inputs.nonband_phosphate_fraction > 1 or inputs.band_phosphate_fraction < 0 or
        inputs.band_phosphate_fraction > 1)
        return error.InvalidSurfacePoreConvectiveSoluteInput;
    for (inputs.micropore_inventory_amount) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfacePoreConvectiveSoluteInput;
    for (inputs.macropore_inventory_mol) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSurfacePoreConvectiveSoluteInput;
}

fn fixture(micropore: []const f64, macropore: []const f64, water_flux: f64) Inputs {
    return .{ .macropore_to_micropore_water_flux_m3_per_step = water_flux, .micropore_water_m3 = 4, .macropore_water_m3 = 2, .minimum_water_m3 = 0.01, .maximum_convective_fraction = 0.75, .nonband_phosphate_fraction = 0.25, .band_phosphate_fraction = 0.5, .micropore_inventory_amount = micropore, .macropore_inventory_mol = macropore };
}

test "TRNSFRS positive flow maps compact macropore topology and partitions P" {
    const micropore = [_]f64{99} ** micropore_species_count;
    const macropore = [_]f64{8} ** macropore_species_count;
    var flux = [_]f64{9} ** micropore_species_count;
    try calculate(fixture(&micropore, &macropore, 1), &flux);
    try std.testing.expectEqual(@as(f64, 4), flux[0]);
    try std.testing.expectEqual(@as(f64, 0), flux[33]);
    try std.testing.expectEqual(@as(f64, 1), flux[34]);
    try std.testing.expectEqual(@as(f64, 2), flux[42]);
}

test "TRNSFRS negative flow uses full micropore topology and sign" {
    var micropore = [_]f64{8} ** micropore_species_count;
    micropore[33] = 12;
    const macropore = [_]f64{99} ** macropore_species_count;
    var flux = [_]f64{0} ** micropore_species_count;
    try calculate(fixture(&micropore, &macropore, -2), &flux);
    try std.testing.expectEqual(@as(f64, -4), flux[0]);
    try std.testing.expectEqual(@as(f64, -6), flux[33]);
    try std.testing.expectEqual(@as(f64, -1), flux[34]);
    try std.testing.expectEqual(@as(f64, -2), flux[42]);
}

test "zero exchange explicitly clears legacy stale H4SiO4 and all fluxes" {
    const micropore = [_]f64{8} ** micropore_species_count;
    const macropore = [_]f64{8} ** macropore_species_count;
    var flux = [_]f64{9} ** micropore_species_count;
    try calculate(fixture(&micropore, &macropore, 0), &flux);
    try std.testing.expectEqualSlices(f64, &([_]f64{0} ** micropore_species_count), &flux);
}

test "late invalid compact donor leaves output atomic" {
    const micropore = [_]f64{8} ** micropore_species_count;
    var macropore = [_]f64{8} ** macropore_species_count;
    macropore[macropore_species_count - 1] = std.math.inf(f64);
    var flux = [_]f64{9} ** micropore_species_count;
    try std.testing.expectError(error.NonFiniteSurfacePoreConvectiveSoluteInput, calculate(fixture(&micropore, &macropore, 1), &flux));
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** micropore_species_count), &flux);
}
