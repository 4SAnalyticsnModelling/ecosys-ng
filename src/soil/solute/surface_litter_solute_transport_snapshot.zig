const std = @import("std");
const species_registry = @import("transport_species.zig");

pub const surface_species_count = species_registry.index(.band_phosphate);

/// Exact source-order translation of TRNSFRS.F lines 374--415.
///
/// Values are extensive mol cell-1. Band-phosphate species begin after this
/// surface-litter registry and are intentionally absent, matching the source
/// arrays whose band variants have no layer-zero element.
pub fn capture(
    dynamic_salt_enabled: bool,
    current_inventory_mol: []const f64,
    transport_inventory_mol: []f64,
) !bool {
    if (current_inventory_mol.len != surface_species_count or
        transport_inventory_mol.len != surface_species_count)
        return error.SurfaceLitterSoluteSnapshotDimensionMismatch;
    if (!dynamic_salt_enabled) return false;

    for (current_inventory_mol) |amount_mol| {
        if (!std.math.isFinite(amount_mol))
            return error.NonFiniteSurfaceLitterSoluteInventory;
    }
    for (current_inventory_mol, transport_inventory_mol) |amount_mol, *snapshot_mol|
        snapshot_mol.* = amount_mol;
    return true;
}

test "TRNSFRS surface litter snapshot preserves all source species order" {
    var current: [surface_species_count]f64 = undefined;
    var snapshot = [_]f64{-1} ** surface_species_count;
    for (&current, 0..) |*amount_mol, index| amount_mol.* = @floatFromInt(index + 1);

    try std.testing.expect(try capture(true, &current, &snapshot));
    try std.testing.expectEqualSlices(f64, &current, &snapshot);
    try std.testing.expectEqual(
        @as(f64, 1),
        snapshot[species_registry.index(.aluminum)],
    );
    try std.testing.expectEqual(
        @as(f64, 42),
        snapshot[species_registry.index(.non_band_magnesium_hpo4)],
    );
}

test "TRNSFRS fixed salt mode retains transport snapshot" {
    const current = [_]f64{2} ** surface_species_count;
    var snapshot = [_]f64{7} ** surface_species_count;
    try std.testing.expect(!try capture(false, &current, &snapshot));
    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{7} ** surface_species_count),
        &snapshot,
    );
}

test "surface litter solute snapshot failure is atomic" {
    var current = [_]f64{2} ** surface_species_count;
    var snapshot = [_]f64{7} ** surface_species_count;
    current[17] = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterSoluteInventory,
        capture(true, &current, &snapshot),
    );
    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{7} ** surface_species_count),
        &snapshot,
    );
}

test "surface litter transport registry excludes band phosphate" {
    try std.testing.expectEqual(@as(usize, 42), surface_species_count);
    try std.testing.expectEqual(
        surface_species_count,
        species_registry.index(.band_phosphate),
    );
}
