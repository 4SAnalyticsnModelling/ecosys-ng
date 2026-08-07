const std = @import("std");
const destination_transfer = @import("sediment_transfer.zig");

pub const SedimentPool = destination_transfer.SedimentPool;
pub const Pools = destination_transfer.Pools;
pub const sediment_pool_count = destination_transfer.sediment_pool_count;
pub const unit = destination_transfer.unit;

fn poolIndex(pools: Pools, pool: SedimentPool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}
fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 9296--9327: full-fraction source sediment.
pub fn reduce(source_layer: usize, redistribution_fraction: f64, remaining_fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or
        pools.amounts.len != sediment_pool_count * pools.layer_count)
        return error.PondSourceSedimentReductionDimensionMismatch;
    if (!finiteSlice(pools.amounts) or !std.math.isFinite(redistribution_fraction) or
        redistribution_fraction < 0 or redistribution_fraction > 1 or !std.math.isFinite(remaining_fraction) or
        remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceSedimentReductionInput;
    if (redistribution_fraction != 1.0 or source_layer == 0) return;

    inline for (std.meta.fields(SedimentPool)) |field| {
        const pool: SedimentPool = @enumFromInt(field.value);
        const index = poolIndex(pools, pool, source_layer);
        if (!std.math.isFinite(remaining_fraction * pools.amounts[index]))
            return error.NonFinitePondSourceSedimentReductionResult;
    }
    inline for (std.meta.fields(SedimentPool)) |field| {
        const pool: SedimentPool = @enumFromInt(field.value);
        const index = poolIndex(pools, pool, source_layer);
        pools.amounts[index] = remaining_fraction * pools.amounts[index];
    }
}

test "REDIST full-fraction source sediment scales exact 29-pool order" {
    var amounts: [sediment_pool_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .amounts = &amounts };
    inline for (std.meta.fields(SedimentPool), 0..) |field, ordinal| {
        const pool: SedimentPool = @enumFromInt(field.value);
        amounts[poolIndex(pools, pool, 1)] = @floatFromInt(4 * (ordinal + 1));
    }
    try reduce(1, 1, 0.25, pools);
    inline for (std.meta.fields(SedimentPool), 0..) |field, ordinal| {
        const pool: SedimentPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, @floatFromInt(ordinal + 1)), amounts[poolIndex(pools, pool, 1)]);
    }
}

test "REDIST source sediment requires exact full fraction and nonzero layer" {
    var amounts: [sediment_pool_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts = &amounts };
    try reduce(1, 0.999999999, 0, pools);
    try reduce(0, 1, 0, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST source sediment runtime dimensions are exact" {
    var amounts: [sediment_pool_count * 2 - 1]f64 = @splat(0);
    try std.testing.expectError(error.PondSourceSedimentReductionDimensionMismatch, reduce(1, 1, 0.5, .{ .layer_count = 2, .amounts = &amounts }));
}

test "REDIST source sediment validation is atomic" {
    var amounts: [sediment_pool_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts = &amounts };
    amounts[poolIndex(pools, .fertilizer_potassium_silicate, 1)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceSedimentReductionInput, reduce(1, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[poolIndex(pools, .sand, 1)]);
}
