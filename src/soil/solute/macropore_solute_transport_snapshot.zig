const std = @import("std");

pub const pre_phosphate_species_count = 33;
pub const phosphate_species_count = 16;
pub const macropore_species_count = pre_phosphate_species_count + phosphate_species_count;

/// Exact source-order translation of TRNSFRS.F lines 1563--1611.
///
/// The compact macropore topology is 33 salt/complex fields followed directly
/// by eight non-band and eight band phosphate fields. H4SiO4 has no macropore
/// field and therefore no phantom slot is introduced. Units are mol layer-1.
pub fn capture(
    current_inventory_mol: []const f64,
    transport_inventory_mol: []f64,
) !void {
    if (current_inventory_mol.len != macropore_species_count or
        transport_inventory_mol.len != macropore_species_count)
        return error.SoilMacroporeSoluteSnapshotDimensionMismatch;
    for (current_inventory_mol) |amount_mol| {
        if (!std.math.isFinite(amount_mol)) return error.NonFiniteSoilMacroporeSoluteInventory;
    }
    for (current_inventory_mol, transport_inventory_mol) |amount_mol, *snapshot_mol|
        snapshot_mol.* = amount_mol;
}

test "TRNSFRS macropore snapshot preserves all 49 contiguous source fields" {
    var current: [macropore_species_count]f64 = undefined;
    var snapshot = [_]f64{-1} ** macropore_species_count;
    for (&current, 0..) |*amount_mol, index| amount_mol.* = @floatFromInt(index + 1);
    try capture(&current, &snapshot);
    try std.testing.expectEqualSlices(f64, &current, &snapshot);
}

test "TRNSFRS macropore phosphate begins immediately after potassium sulfate" {
    var current = [_]f64{0} ** macropore_species_count;
    var snapshot = [_]f64{-1} ** macropore_species_count;
    current[pre_phosphate_species_count - 1] = 3;
    current[pre_phosphate_species_count] = 5;
    current[macropore_species_count - 1] = 7;
    try capture(&current, &snapshot);
    try std.testing.expectEqual(@as(f64, 3), snapshot[32]);
    try std.testing.expectEqual(@as(f64, 5), snapshot[33]);
    try std.testing.expectEqual(@as(f64, 7), snapshot[48]);
}

test "macropore snapshot late failure is atomic" {
    var current = [_]f64{2} ** macropore_species_count;
    var snapshot = [_]f64{9} ** macropore_species_count;
    current[macropore_species_count - 1] = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteSoilMacroporeSoluteInventory,
        capture(&current, &snapshot),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{9} ** macropore_species_count), &snapshot);
}

test "macropore snapshot rejects phantom 50th species" {
    const current = [_]f64{2} ** (macropore_species_count + 1);
    var snapshot = [_]f64{9} ** macropore_species_count;
    try std.testing.expectError(
        error.SoilMacroporeSoluteSnapshotDimensionMismatch,
        capture(&current, &snapshot),
    );
}
