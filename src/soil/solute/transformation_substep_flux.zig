const std = @import("std");
const species_registry = @import("transport_species.zig");

pub const species_count = species_registry.AqueousSpecies.count;

/// Exact source-order translation of TRNSFRS.F lines 1277--1327.
///
/// Arrays are active-mineral-layer-major, then the complete 50-species
/// registry. SOLUTE transformations and published transport transformations
/// are extensive mol step-1. Legacy chloride is forced to zero at each layer.
pub fn publish(
    layer_count: usize,
    solute_transformation_mol_per_step: []const f64,
    transport_transformation_mol_per_step: []f64,
    substep_fraction: f64,
) !void {
    const extent = std.math.mul(usize, layer_count, species_count) catch
        return error.SoilSoluteTransformationDimensionOverflow;
    if (solute_transformation_mol_per_step.len != extent or
        transport_transformation_mol_per_step.len != extent)
        return error.SoilSoluteTransformationDimensionMismatch;
    if (!std.math.isFinite(substep_fraction) or substep_fraction < 0 or substep_fraction > 1)
        return error.InvalidSoluteSubstepFraction;

    for (0..layer_count) |layer| {
        const start = layer * species_count;
        for (solute_transformation_mol_per_step[start..][0..species_count], 0..) |transformation_mol, species| {
            if (species == species_registry.index(.chloride)) continue;
            if (!std.math.isFinite(transformation_mol) or
                !std.math.isFinite(-transformation_mol * substep_fraction))
                return error.NonFiniteSoilSoluteTransformation;
        }
    }
    for (0..layer_count) |layer| {
        const start = layer * species_count;
        for (solute_transformation_mol_per_step[start..][0..species_count], transport_transformation_mol_per_step[start..][0..species_count], 0..) |transformation_mol, *published_mol, species| {
            published_mol.* = if (species == species_registry.index(.chloride))
                0
            else
                -transformation_mol * substep_fraction;
        }
    }
}

test "TRNSFRS soil transformation accepts an empty active-layer range" {
    try publish(0, &.{}, &.{}, 0.5);
}

test "TRNSFRS soil transformations preserve runtime layer-major signs and order" {
    const layer_count = 3;
    var transformations: [layer_count * species_count]f64 = undefined;
    var published = [_]f64{99} ** (layer_count * species_count);
    for (&transformations, 0..) |*value, index| value.* = @floatFromInt(index + 1);
    try publish(layer_count, &transformations, &published, 0.25);
    for (0..layer_count) |layer| {
        const base = layer * species_count;
        try std.testing.expectEqual(-transformations[base] * 0.25, published[base]);
        try std.testing.expectEqual(
            -transformations[base + species_registry.index(.band_magnesium_hpo4)] * 0.25,
            published[base + species_registry.index(.band_magnesium_hpo4)],
        );
    }
}

test "TRNSFRS soil chloride transformation is forced zero per layer" {
    const layer_count = 2;
    var transformations = [_]f64{2} ** (layer_count * species_count);
    var published = [_]f64{99} ** (layer_count * species_count);
    for (0..layer_count) |layer|
        transformations[layer * species_count + species_registry.index(.chloride)] = std.math.nan(f64);
    try publish(layer_count, &transformations, &published, 0.5);
    for (0..layer_count) |layer|
        try std.testing.expectEqual(@as(f64, 0), published[layer * species_count + species_registry.index(.chloride)]);
}

test "soil transformation late-layer failure is atomic" {
    const layer_count = 2;
    var transformations = [_]f64{2} ** (layer_count * species_count);
    var published = [_]f64{99} ** (layer_count * species_count);
    transformations[transformations.len - 1] = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteSoilSoluteTransformation,
        publish(layer_count, &transformations, &published, 0.5),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{99} ** published.len), &published);
}

test "soil transformation topology follows runtime layers" {
    const transformations = [_]f64{2} ** species_count;
    var published = [_]f64{99} ** species_count;
    try std.testing.expectError(
        error.SoilSoluteTransformationDimensionMismatch,
        publish(2, &transformations, &published, 0.5),
    );
}
