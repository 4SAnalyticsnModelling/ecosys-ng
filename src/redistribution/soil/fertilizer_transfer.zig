const std = @import("std");
const pond_fertilizer = @import("../pond/fertilizer_transfer.zig");

pub const FertilizerPool = pond_fertilizer.FertilizerPool;
pub const Pools = pond_fertilizer.Pools;
pub const fertilizer_pool_count = pond_fertilizer.fertilizer_pool_count;

fn index(pools: Pools, pool: FertilizerPool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}

/// Direct translation of REDIST 9670--9693 under the positive-`BKDS` soil gate.
pub fn transfer(
    current_layer: usize,
    source_layer: usize,
    destination_layer: usize,
    source_bulk_density_megagrams_m3: f64,
    destination_bulk_density_megagrams_m3: f64,
    redistribution_fraction: f64,
    pools: Pools,
) !void {
    if (pools.layer_count == 0 or current_layer >= pools.layer_count or source_layer >= pools.layer_count or
        destination_layer >= pools.layer_count or source_layer == destination_layer or
        pools.amounts_mol_n.len != fertilizer_pool_count * pools.layer_count)
        return error.SoilFertilizerTransferDimensionMismatch;
    if (!std.math.isFinite(source_bulk_density_megagrams_m3) or !std.math.isFinite(destination_bulk_density_megagrams_m3) or
        !std.math.isFinite(redistribution_fraction) or redistribution_fraction < 0 or redistribution_fraction > 1)
        return error.InvalidSoilFertilizerTransferInput;
    for (pools.amounts_mol_n) |amount| if (!std.math.isFinite(amount))
        return error.InvalidSoilFertilizerTransferInput;
    if (source_bulk_density_megagrams_m3 <= 0 or destination_bulk_density_megagrams_m3 <= 0) return;

    inline for (std.meta.fields(FertilizerPool)) |field| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        const moved_mol_n = @min(
            redistribution_fraction * pools.amounts_mol_n[index(pools, pool, current_layer)],
            pools.amounts_mol_n[index(pools, pool, source_layer)],
        );
        const next_destination = pools.amounts_mol_n[index(pools, pool, destination_layer)] + moved_mol_n;
        const next_source = pools.amounts_mol_n[index(pools, pool, source_layer)] - moved_mol_n;
        if (!std.math.isFinite(moved_mol_n) or !std.math.isFinite(next_destination) or !std.math.isFinite(next_source))
            return error.NonFiniteSoilFertilizerTransferResult;
    }
    inline for (std.meta.fields(FertilizerPool)) |field| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        const current = index(pools, pool, current_layer);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved_mol_n = @min(redistribution_fraction * pools.amounts_mol_n[current], pools.amounts_mol_n[source]);
        pools.amounts_mol_n[destination] = pools.amounts_mol_n[destination] + moved_mol_n;
        pools.amounts_mol_n[source] = pools.amounts_mol_n[source] - moved_mol_n;
    }
}

test "REDIST soil fertilizer uses current amount cap in exact pool order" {
    var amounts: [fertilizer_pool_count * 3]f64 = @splat(0);
    const pools = Pools{ .layer_count = 3, .amounts_mol_n = &amounts };
    inline for (std.meta.fields(FertilizerPool), 0..) |field, ordinal| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        amounts[index(pools, pool, 0)] = @floatFromInt(4 * (ordinal + 1));
        amounts[index(pools, pool, 1)] = @floatFromInt(ordinal + 1);
        amounts[index(pools, pool, 2)] = 10;
    }
    try transfer(0, 1, 2, 1, 1, 0.5, pools);
    inline for (std.meta.fields(FertilizerPool), 0..) |field, ordinal| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        const source_before: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(@as(f64, 0), amounts[index(pools, pool, 1)]);
        try std.testing.expectEqual(10 + source_before, amounts[index(pools, pool, 2)]);
    }
}

test "REDIST soil fertilizer conserves all eight pools" {
    var amounts: [fertilizer_pool_count * 3]f64 = @splat(0);
    const pools = Pools{ .layer_count = 3, .amounts_mol_n = &amounts };
    inline for (std.meta.fields(FertilizerPool)) |field| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        amounts[index(pools, pool, 0)] = 2;
        amounts[index(pools, pool, 1)] = 4;
        amounts[index(pools, pool, 2)] = 8;
    }
    try transfer(0, 1, 2, 1, 1, 0.5, pools);
    inline for (std.meta.fields(FertilizerPool)) |field| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 12), amounts[index(pools, pool, 1)] + amounts[index(pools, pool, 2)]);
    }
}

test "REDIST soil fertilizer validation is atomic" {
    var amounts: [fertilizer_pool_count * 3]f64 = @splat(1);
    const pools = Pools{ .layer_count = 3, .amounts_mol_n = &amounts };
    amounts[index(pools, .banded_nitrate, 1)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilFertilizerTransferInput, transfer(0, 1, 2, 1, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .broadcast_ammonium, 1)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .broadcast_ammonium, 2)]);
}
