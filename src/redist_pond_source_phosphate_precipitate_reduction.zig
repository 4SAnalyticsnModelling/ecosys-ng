const std = @import("std");
const destination_transfer = @import("redist_pond_phosphate_precipitate_transfer.zig");

pub const PrecipitatePool = destination_transfer.PrecipitatePool;
pub const Pools = destination_transfer.Pools;
pub const precipitate_pool_count = destination_transfer.precipitate_pool_count;

fn index(pools: Pools, pool: PrecipitatePool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}

/// Direct translation of REDIST 9360--9369 under the enclosing `FX == 1.0`.
pub fn reduce(source_layer: usize, redistribution_fraction: f64, remaining_fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or
        pools.amounts_mol_p.len != precipitate_pool_count * pools.layer_count)
        return error.PondSourcePhosphatePrecipitateReductionDimensionMismatch;
    for (pools.amounts_mol_p) |amount| if (!std.math.isFinite(amount))
        return error.InvalidPondSourcePhosphatePrecipitateReductionInput;
    if (!std.math.isFinite(redistribution_fraction) or redistribution_fraction < 0 or redistribution_fraction > 1 or
        !std.math.isFinite(remaining_fraction) or remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourcePhosphatePrecipitateReductionInput;
    if (redistribution_fraction != 1.0) return;

    inline for (std.meta.fields(PrecipitatePool)) |field| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        if (!std.math.isFinite(remaining_fraction * pools.amounts_mol_p[index(pools, pool, source_layer)]))
            return error.NonFinitePondSourcePhosphatePrecipitateReductionResult;
    }
    inline for (std.meta.fields(PrecipitatePool)) |field| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        pools.amounts_mol_p[source] = remaining_fraction * pools.amounts_mol_p[source];
    }
}

test "REDIST source phosphate precipitates scale exact non-band then band order" {
    var amounts: [precipitate_pool_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .amounts_mol_p = &amounts };
    inline for (std.meta.fields(PrecipitatePool), 0..) |field, ordinal| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        amounts[index(pools, pool, 0)] = @floatFromInt(4 * (ordinal + 1));
        amounts[index(pools, pool, 1)] = 99;
    }
    try reduce(0, 1, 0.25, pools);
    inline for (std.meta.fields(PrecipitatePool), 0..) |field, ordinal| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, @floatFromInt(ordinal + 1)), amounts[index(pools, pool, 0)]);
        try std.testing.expectEqual(@as(f64, 99), amounts[index(pools, pool, 1)]);
    }
}

test "REDIST source phosphate precipitates permit layer zero and require exact full fraction" {
    var amounts: [precipitate_pool_count]f64 = @splat(1);
    const pools = Pools{ .layer_count = 1, .amounts_mol_p = &amounts };
    try reduce(0, 0.999999999, 0, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
    try reduce(0, 1, 0, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 0), amount);
}

test "REDIST source phosphate precipitate validation is atomic" {
    var amounts: [precipitate_pool_count]f64 = @splat(1);
    const pools = Pools{ .layer_count = 1, .amounts_mol_p = &amounts };
    amounts[index(pools, .banded_monocalcium_phosphate, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourcePhosphatePrecipitateReductionInput, reduce(0, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .nonband_aluminum_phosphate, 0)]);
}
