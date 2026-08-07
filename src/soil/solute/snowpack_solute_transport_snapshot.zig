const std = @import("std");
const species_registry = @import("transport_species.zig");

pub const snow_species_count = species_registry.index(.band_phosphate) - 1;

/// Exact source-order translation of TRNSFRS.F lines 1202--1244.
///
/// Arrays are runtime snow-layer-major, then species, and contain extensive
/// mol layer-1. TRNSFRS has no snow H4SiO4 array, so the runtime axis is the
/// 33 salt fields followed immediately by the eight non-band phosphate fields.
pub fn capture(
    layer_count: usize,
    current_inventory_mol: []const f64,
    transport_inventory_mol: []f64,
) !void {
    const extent = std.math.mul(usize, layer_count, snow_species_count) catch
        return error.SnowpackSoluteSnapshotDimensionOverflow;
    if (current_inventory_mol.len != extent or transport_inventory_mol.len != extent)
        return error.SnowpackSoluteSnapshotDimensionMismatch;

    for (0..layer_count) |layer| {
        const start = layer * snow_species_count;
        for (current_inventory_mol[start..][0..snow_species_count]) |amount_mol| {
            if (!std.math.isFinite(amount_mol))
                return error.NonFiniteSnowpackSoluteInventory;
        }
    }
    for (0..layer_count) |layer| {
        const start = layer * snow_species_count;
        for (current_inventory_mol[start..][0..snow_species_count], transport_inventory_mol[start..][0..snow_species_count]) |amount_mol, *snapshot_mol| {
            snapshot_mol.* = amount_mol;
        }
    }
}

test "TRNSFRS snow snapshot accepts the Fortran zero-trip layer loop" {
    try capture(0, &.{}, &.{});
}

test "TRNSFRS snow snapshot preserves runtime layer-major source order" {
    const layer_count = 3;
    var current: [layer_count * snow_species_count]f64 = undefined;
    var snapshot = [_]f64{-1} ** (layer_count * snow_species_count);
    for (&current, 0..) |*amount_mol, index| amount_mol.* = @floatFromInt(index + 1);
    try capture(layer_count, &current, &snapshot);
    for (0..layer_count) |layer| {
        const start = layer * snow_species_count;
        for (0..snow_species_count) |species| {
            try std.testing.expectEqual(current[start + species], snapshot[start + species]);
        }
    }
}

test "snowpack snapshot late failure is atomic" {
    const layer_count = 2;
    var current = [_]f64{2} ** (layer_count * snow_species_count);
    var snapshot = [_]f64{7} ** (layer_count * snow_species_count);
    current[current.len - 1] = std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteSnowpackSoluteInventory,
        capture(layer_count, &current, &snapshot),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{7} ** snapshot.len), &snapshot);
}

test "snowpack snapshot dimensions follow runtime layer count" {
    const current = [_]f64{2} ** snow_species_count;
    var snapshot = [_]f64{7} ** snow_species_count;
    try std.testing.expectError(
        error.SnowpackSoluteSnapshotDimensionMismatch,
        capture(2, &current, &snapshot),
    );
}
