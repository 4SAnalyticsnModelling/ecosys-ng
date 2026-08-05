const std = @import("std");

/// Source order is 34 salt/complex species (including H4SiO4), eight
/// non-band phosphate species, then eight band phosphate species.
pub const species_per_layer = 50;

pub const ActiveLayerRange = struct {
    first_layer_index: usize,
    end_layer_index: usize,
};

pub const Handoff = struct {
    total_layer_count: usize,
    active_layers: ActiveLayerRange,
    /// Transport inventory before the current substep, mol layer-1.
    transport_inventory_mol_per_layer: []f64,
    /// SOLUTE transformation for the current substep, mol layer-1 step-1.
    transformation_mol_per_layer_step: []const f64,
};

/// Exact source-order translation of TRNSFRS.F lines 1880--1930.
///
/// For each active layer, the SOLUTE transformation is subtracted from every
/// transport inventory in the compact source topology. All inputs and results
/// are validated before mutation so a late invalid species cannot contaminate
/// preceding state.
pub fn apply(handoff: Handoff) !void {
    const range = handoff.active_layers;
    if (range.first_layer_index > range.end_layer_index or
        range.end_layer_index > handoff.total_layer_count)
        return error.InvalidActiveSoilLayerRange;
    const required_len = std.math.mul(
        usize,
        handoff.total_layer_count,
        species_per_layer,
    ) catch return error.SoilSoluteTransformationDimensionOverflow;
    if (handoff.transport_inventory_mol_per_layer.len != required_len or
        handoff.transformation_mol_per_layer_step.len != required_len)
        return error.SoilSoluteTransformationDimensionMismatch;

    const start = range.first_layer_index * species_per_layer;
    const end = range.end_layer_index * species_per_layer;
    for (start..end) |index| {
        const inventory_mol = handoff.transport_inventory_mol_per_layer[index];
        const transformation_mol = handoff.transformation_mol_per_layer_step[index];
        if (!std.math.isFinite(inventory_mol) or !std.math.isFinite(transformation_mol))
            return error.NonFiniteSoilSoluteTransformationInput;
        if (!std.math.isFinite(inventory_mol - transformation_mol))
            return error.NonFiniteSoilSoluteTransformationResult;
    }
    for (start..end) |index| {
        handoff.transport_inventory_mol_per_layer[index] =
            handoff.transport_inventory_mol_per_layer[index] -
            handoff.transformation_mol_per_layer_step[index];
    }
}

test "TRNSFRS subtracts all 50 transformations in source order per active layer" {
    const layer_count = 2;
    var inventory = [_]f64{100} ** (layer_count * species_per_layer);
    var transformations: [layer_count * species_per_layer]f64 = undefined;
    for (&transformations, 0..) |*value, index| value.* = @floatFromInt(index + 1);
    try apply(.{
        .total_layer_count = layer_count,
        .active_layers = .{ .first_layer_index = 1, .end_layer_index = 2 },
        .transport_inventory_mol_per_layer = &inventory,
        .transformation_mol_per_layer_step = &transformations,
    });
    try std.testing.expectEqual(@as(f64, 100), inventory[species_per_layer - 1]);
    try std.testing.expectEqual(@as(f64, 49), inventory[species_per_layer]);
    try std.testing.expectEqual(@as(f64, 0), inventory[layer_count * species_per_layer - 1]);
}

test "TRNSFRS soil handoff preserves an empty active-layer range" {
    var inventory = [_]f64{10} ** species_per_layer;
    const transformations = [_]f64{2} ** species_per_layer;
    try apply(.{
        .total_layer_count = 1,
        .active_layers = .{ .first_layer_index = 1, .end_layer_index = 1 },
        .transport_inventory_mol_per_layer = &inventory,
        .transformation_mol_per_layer_step = &transformations,
    });
    try std.testing.expectEqualSlices(f64, &([_]f64{10} ** species_per_layer), &inventory);
}

test "soil transformation late non-finite result is atomic" {
    var inventory = [_]f64{10} ** species_per_layer;
    var transformations = [_]f64{2} ** species_per_layer;
    transformations[species_per_layer - 1] = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteSoilSoluteTransformationInput,
        apply(.{
            .total_layer_count = 1,
            .active_layers = .{ .first_layer_index = 0, .end_layer_index = 1 },
            .transport_inventory_mol_per_layer = &inventory,
            .transformation_mol_per_layer_step = &transformations,
        }),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{10} ** species_per_layer), &inventory);
}

test "soil transformation topology rejects a missing field" {
    var inventory = [_]f64{10} ** (species_per_layer - 1);
    const transformations = [_]f64{2} ** species_per_layer;
    try std.testing.expectError(
        error.SoilSoluteTransformationDimensionMismatch,
        apply(.{
            .total_layer_count = 1,
            .active_layers = .{ .first_layer_index = 0, .end_layer_index = 1 },
            .transport_inventory_mol_per_layer = &inventory,
            .transformation_mol_per_layer_step = &transformations,
        }),
    );
}

test "invalid soil layer range is rejected atomically" {
    var inventory = [_]f64{10} ** species_per_layer;
    const transformations = [_]f64{2} ** species_per_layer;
    try std.testing.expectError(
        error.InvalidActiveSoilLayerRange,
        apply(.{
            .total_layer_count = 1,
            .active_layers = .{ .first_layer_index = 1, .end_layer_index = 0 },
            .transport_inventory_mol_per_layer = &inventory,
            .transformation_mol_per_layer_step = &transformations,
        }),
    );
    try std.testing.expectEqual(@as(f64, 10), inventory[0]);
}
