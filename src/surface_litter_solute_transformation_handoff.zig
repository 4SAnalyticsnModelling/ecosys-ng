const std = @import("std");
const species_registry = @import("solute_transport_species.zig");

pub const primary_species_count = species_registry.index(.aluminum_hydroxide_1);

/// Exact source-order translation of TRNSFRS.F lines 1253--1276.
///
/// Inventories and transformations are extensive mol cell-1 and mol step-1.
/// The source derives all twelve RZ* values first, then subtracts them from
/// all twelve transport inventories. Chloride and carbonate are explicitly
/// zero and do not consume their supplied SOLUTE transformation slots.
pub fn apply(
    transport_inventory_mol: []f64,
    solute_transformation_mol_per_step: []const f64,
    transport_transformation_mol_per_step: []f64,
    substep_fraction: f64,
) !void {
    if (transport_inventory_mol.len != primary_species_count or
        solute_transformation_mol_per_step.len != primary_species_count or
        transport_transformation_mol_per_step.len != primary_species_count)
        return error.SurfaceLitterTransformationDimensionMismatch;
    if (!std.math.isFinite(substep_fraction) or substep_fraction < 0 or substep_fraction > 1)
        return error.InvalidSoluteSubstepFraction;

    for (0..primary_species_count) |index| {
        if (!std.math.isFinite(transport_inventory_mol[index]))
            return error.NonFiniteSurfaceLitterTransformationInput;
        const reaction_mol = sourceReaction(index, solute_transformation_mol_per_step, substep_fraction) catch |err| return err;
        if (!std.math.isFinite(transport_inventory_mol[index] - reaction_mol))
            return error.NonFiniteSurfaceLitterTransformationResult;
    }
    for (0..primary_species_count) |index|
        transport_transformation_mol_per_step[index] = try sourceReaction(index, solute_transformation_mol_per_step, substep_fraction);
    for (0..primary_species_count) |index|
        transport_inventory_mol[index] = transport_inventory_mol[index] - transport_transformation_mol_per_step[index];
}

fn sourceReaction(index: usize, transformations: []const f64, fraction: f64) !f64 {
    if (index == species_registry.index(.chloride) or index == species_registry.index(.carbonate)) return 0;
    if (!std.math.isFinite(transformations[index]))
        return error.NonFiniteSurfaceLitterTransformationInput;
    const reaction_mol = -transformations[index] * fraction;
    if (!std.math.isFinite(reaction_mol))
        return error.NonFiniteSurfaceLitterTransformationResult;
    return reaction_mol;
}

test "TRNSFRS surface transformation preserves two-pass signs and order" {
    var inventory = [_]f64{10} ** primary_species_count;
    var transformations: [primary_species_count]f64 = undefined;
    var reaction = [_]f64{99} ** primary_species_count;
    for (&transformations, 0..) |*value, index| value.* = @floatFromInt(index + 1);
    try apply(&inventory, &transformations, &reaction, 0.25);
    try std.testing.expectEqual(@as(f64, -0.25), reaction[species_registry.index(.aluminum)]);
    try std.testing.expectEqual(@as(f64, 10.25), inventory[species_registry.index(.aluminum)]);
    try std.testing.expectEqual(@as(f64, -3), reaction[species_registry.index(.bicarbonate)]);
    try std.testing.expectEqual(@as(f64, 13), inventory[species_registry.index(.bicarbonate)]);
}

test "TRNSFRS surface chloride and carbonate transformations are forced zero" {
    var inventory = [_]f64{10} ** primary_species_count;
    var transformations = [_]f64{2} ** primary_species_count;
    var reaction = [_]f64{99} ** primary_species_count;
    transformations[species_registry.index(.chloride)] = std.math.nan(f64);
    transformations[species_registry.index(.carbonate)] = std.math.inf(f64);
    try apply(&inventory, &transformations, &reaction, 0.5);
    inline for (.{ species_registry.index(.chloride), species_registry.index(.carbonate) }) |index| {
        try std.testing.expectEqual(@as(f64, 0), reaction[index]);
        try std.testing.expectEqual(@as(f64, 10), inventory[index]);
    }
}

test "surface transformation late failure is atomic" {
    var inventory = [_]f64{10} ** primary_species_count;
    var transformations = [_]f64{2} ** primary_species_count;
    var reaction = [_]f64{99} ** primary_species_count;
    transformations[primary_species_count - 1] = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterTransformationInput,
        apply(&inventory, &transformations, &reaction, 0.5),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{10} ** primary_species_count), &inventory);
    try std.testing.expectEqualSlices(f64, &([_]f64{99} ** primary_species_count), &reaction);
}

test "surface transformation topology is exact" {
    var inventory = [_]f64{10} ** (primary_species_count - 1);
    const transformations = [_]f64{2} ** primary_species_count;
    var reaction = [_]f64{99} ** primary_species_count;
    try std.testing.expectError(
        error.SurfaceLitterTransformationDimensionMismatch,
        apply(&inventory, &transformations, &reaction, 0.5),
    );
}
