const std = @import("std");
const species_registry = @import("solute_transport_species.zig");

pub const species_count = species_registry.AqueousSpecies.count;

/// Exact source-order translation of TRNSFRS.F lines 1513--1562.
///
/// The caller selects one runtime mineral layer from the enclosing source
/// loop. Current and working inventories contain the complete 50-species
/// micropore registry and use extensive mol layer-1.
pub fn capture(
    current_inventory_mol: []const f64,
    transport_inventory_mol: []f64,
) !void {
    if (current_inventory_mol.len != species_count or transport_inventory_mol.len != species_count)
        return error.SoilMicroporeSoluteSnapshotDimensionMismatch;
    for (current_inventory_mol) |amount_mol| {
        if (!std.math.isFinite(amount_mol)) return error.NonFiniteSoilMicroporeSoluteInventory;
    }
    for (current_inventory_mol, transport_inventory_mol) |amount_mol, *snapshot_mol|
        snapshot_mol.* = amount_mol;
}

test "TRNSFRS micropore snapshot preserves all 50 source species" {
    var current: [species_count]f64 = undefined;
    var snapshot = [_]f64{-1} ** species_count;
    for (&current, 0..) |*amount_mol, index| amount_mol.* = @floatFromInt(index + 1);
    try capture(&current, &snapshot);
    try std.testing.expectEqualSlices(f64, &current, &snapshot);
}

test "TRNSFRS micropore snapshot includes silicate and band endpoints" {
    var current = [_]f64{0} ** species_count;
    var snapshot = [_]f64{-1} ** species_count;
    current[species_registry.index(.hydrogen_silicate)] = 3;
    current[species_registry.index(.band_phosphate)] = 5;
    current[species_registry.index(.band_magnesium_hpo4)] = 7;
    try capture(&current, &snapshot);
    try std.testing.expectEqual(@as(f64, 3), snapshot[species_registry.index(.hydrogen_silicate)]);
    try std.testing.expectEqual(@as(f64, 5), snapshot[species_registry.index(.band_phosphate)]);
    try std.testing.expectEqual(@as(f64, 7), snapshot[species_registry.index(.band_magnesium_hpo4)]);
}

test "micropore snapshot late failure is atomic" {
    var current = [_]f64{2} ** species_count;
    var snapshot = [_]f64{9} ** species_count;
    current[species_count - 1] = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteSoilMicroporeSoluteInventory,
        capture(&current, &snapshot),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** species_count), &snapshot);
}

test "micropore snapshot topology is exact" {
    const current = [_]f64{2} ** (species_count - 1);
    var snapshot = [_]f64{9} ** species_count;
    try std.testing.expectError(
        error.SoilMicroporeSoluteSnapshotDimensionMismatch,
        capture(&current, &snapshot),
    );
}
