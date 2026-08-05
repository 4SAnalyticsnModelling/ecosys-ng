const std = @import("std");
const destination_transfer = @import("redist_pond_fertilizer_transfer.zig");

pub const FertilizerPool = destination_transfer.FertilizerPool;
pub const Pools = destination_transfer.Pools;
pub const fertilizer_pool_count = destination_transfer.fertilizer_pool_count;

fn poolIndex(pools: Pools, pool: FertilizerPool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}
fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 9331--9338 under the enclosing `FX == 1.0`.
pub fn reduce(source_layer: usize, redistribution_fraction: f64, remaining_fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or
        pools.amounts_mol_n.len != fertilizer_pool_count * pools.layer_count)
        return error.PondSourceFertilizerReductionDimensionMismatch;
    if (!finiteSlice(pools.amounts_mol_n) or !std.math.isFinite(redistribution_fraction) or
        redistribution_fraction < 0 or redistribution_fraction > 1 or !std.math.isFinite(remaining_fraction) or
        remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceFertilizerReductionInput;
    if (redistribution_fraction != 1.0) return;

    inline for (std.meta.fields(FertilizerPool)) |field| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        const index = poolIndex(pools, pool, source_layer);
        if (!std.math.isFinite(remaining_fraction * pools.amounts_mol_n[index]))
            return error.NonFinitePondSourceFertilizerReductionResult;
    }
    inline for (std.meta.fields(FertilizerPool)) |field| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        const index = poolIndex(pools, pool, source_layer);
        pools.amounts_mol_n[index] = remaining_fraction * pools.amounts_mol_n[index];
    }
}

test "REDIST source fertilizer scales exact broadcast and banded order" {
    var amounts: [fertilizer_pool_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .amounts_mol_n = &amounts };
    inline for (std.meta.fields(FertilizerPool), 0..) |field, ordinal| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        amounts[poolIndex(pools, pool, 0)] = @floatFromInt(4 * (ordinal + 1));
    }
    try reduce(0, 1, 0.25, pools);
    inline for (std.meta.fields(FertilizerPool), 0..) |field, ordinal| {
        const pool: FertilizerPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, @floatFromInt(ordinal + 1)), amounts[poolIndex(pools, pool, 0)]);
    }
}

test "REDIST source fertilizer permits layer zero and exhaustion" {
    var amounts: [fertilizer_pool_count]f64 = @splat(1);
    const pools = Pools{ .layer_count = 1, .amounts_mol_n = &amounts };
    try reduce(0, 1, 0, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 0), amount);
}

test "REDIST source fertilizer requires exact full redistribution fraction" {
    var amounts: [fertilizer_pool_count]f64 = @splat(1);
    const pools = Pools{ .layer_count = 1, .amounts_mol_n = &amounts };
    try reduce(0, 0.999999999, 0, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST source fertilizer validation is atomic" {
    var amounts: [fertilizer_pool_count]f64 = @splat(1);
    const pools = Pools{ .layer_count = 1, .amounts_mol_n = &amounts };
    amounts[poolIndex(pools, .banded_nitrate, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceFertilizerReductionInput, reduce(0, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[poolIndex(pools, .broadcast_ammonium, 0)]);
}
