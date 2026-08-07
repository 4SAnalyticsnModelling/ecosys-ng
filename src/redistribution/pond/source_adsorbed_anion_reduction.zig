const std = @import("std");
const destination_transfer = @import("adsorbed_anion_transfer.zig");

pub const AdsorbedAnionPool = destination_transfer.AdsorbedAnionPool;
pub const Pools = destination_transfer.Pools;
pub const adsorbed_anion_pool_count = destination_transfer.adsorbed_anion_pool_count;

fn index(pools: Pools, pool: AdsorbedAnionPool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}

/// Direct translation of REDIST 9347--9356 under the enclosing `FX == 1.0`.
pub fn reduce(source_layer: usize, redistribution_fraction: f64, remaining_fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or
        pools.amounts_mol.len != adsorbed_anion_pool_count * pools.layer_count)
        return error.PondSourceAdsorbedAnionReductionDimensionMismatch;
    for (pools.amounts_mol) |amount| if (!std.math.isFinite(amount))
        return error.InvalidPondSourceAdsorbedAnionReductionInput;
    if (!std.math.isFinite(redistribution_fraction) or redistribution_fraction < 0 or redistribution_fraction > 1 or
        !std.math.isFinite(remaining_fraction) or remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceAdsorbedAnionReductionInput;
    if (redistribution_fraction != 1.0) return;

    inline for (std.meta.fields(AdsorbedAnionPool)) |field| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        if (!std.math.isFinite(remaining_fraction * pools.amounts_mol[index(pools, pool, source_layer)]))
            return error.NonFinitePondSourceAdsorbedAnionReductionResult;
    }
    inline for (std.meta.fields(AdsorbedAnionPool)) |field| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        pools.amounts_mol[source] = remaining_fraction * pools.amounts_mol[source];
    }
}

test "REDIST source adsorbed anions scale exact non-band then band order" {
    var amounts: [adsorbed_anion_pool_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .amounts_mol = &amounts };
    inline for (std.meta.fields(AdsorbedAnionPool), 0..) |field, ordinal| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        amounts[index(pools, pool, 0)] = @floatFromInt(4 * (ordinal + 1));
        amounts[index(pools, pool, 1)] = 99;
    }
    try reduce(0, 1, 0.25, pools);
    inline for (std.meta.fields(AdsorbedAnionPool), 0..) |field, ordinal| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, @floatFromInt(ordinal + 1)), amounts[index(pools, pool, 0)]);
        try std.testing.expectEqual(@as(f64, 99), amounts[index(pools, pool, 1)]);
    }
}

test "REDIST source adsorbed anions permit layer zero and require exact full fraction" {
    var amounts: [adsorbed_anion_pool_count]f64 = @splat(1);
    const pools = Pools{ .layer_count = 1, .amounts_mol = &amounts };
    try reduce(0, 0.999999999, 0, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
    try reduce(0, 1, 0, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 0), amount);
}

test "REDIST source adsorbed anion validation is atomic" {
    var amounts: [adsorbed_anion_pool_count]f64 = @splat(1);
    const pools = Pools{ .layer_count = 1, .amounts_mol = &amounts };
    amounts[index(pools, .banded_dihydrogen_phosphate, 0)] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidPondSourceAdsorbedAnionReductionInput, reduce(0, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .nonband_deprotonated_site, 0)]);
}
