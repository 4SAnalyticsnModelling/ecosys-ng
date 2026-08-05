const std = @import("std");

/// Exact REDIST 8993--9012 order: non-band, then band.
pub const AdsorbedAnionPool = enum(u8) {
    nonband_deprotonated_site,
    nonband_hydroxyl_site,
    nonband_protonated_site,
    nonband_hydrogen_phosphate,
    nonband_dihydrogen_phosphate,
    banded_deprotonated_site,
    banded_hydroxyl_site,
    banded_protonated_site,
    banded_hydrogen_phosphate,
    banded_dihydrogen_phosphate,
};

pub const adsorbed_anion_pool_count = std.meta.fields(AdsorbedAnionPool).len;

pub const Pools = struct {
    layer_count: usize,
    /// Pool-major extensive adsorption-site or phosphate amounts, mol.
    amounts_mol: []f64,

    fn index(self: Pools, pool: AdsorbedAnionPool, layer: usize) usize {
        return @intFromEnum(pool) * self.layer_count + layer;
    }
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8993--9012 under `FX == 1.0`.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.amounts_mol.len != adsorbed_anion_pool_count * pools.layer_count)
        return error.PondAdsorbedAnionTransferDimensionMismatch;
    if (!finiteSlice(pools.amounts_mol) or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondAdsorbedAnionTransferInput;
    if (fraction != 1.0) return;

    inline for (std.meta.fields(AdsorbedAnionPool)) |field| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        const source = pools.index(pool, source_layer);
        const destination = pools.index(pool, destination_layer);
        if (!std.math.isFinite(pools.amounts_mol[destination] + fraction * pools.amounts_mol[source]))
            return error.NonFinitePondAdsorbedAnionTransferResult;
    }
    inline for (std.meta.fields(AdsorbedAnionPool)) |field| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        const source = pools.index(pool, source_layer);
        const destination = pools.index(pool, destination_layer);
        pools.amounts_mol[destination] = pools.amounts_mol[destination] + fraction * pools.amounts_mol[source];
    }
}

test "REDIST full-fraction adsorbed anions transfer exact non-band band order" {
    var amounts: [adsorbed_anion_pool_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .amounts_mol = &amounts };
    inline for (std.meta.fields(AdsorbedAnionPool), 0..) |field, ordinal| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        amounts[pools.index(pool, 0)] = @floatFromInt(ordinal + 1);
        amounts[pools.index(pool, 1)] = @floatFromInt((ordinal + 1) * 10);
    }
    try transfer(0, 1, 1, pools);
    inline for (std.meta.fields(AdsorbedAnionPool), 0..) |field, ordinal| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        const factor: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(11 * factor, amounts[pools.index(pool, 1)]);
        try std.testing.expectEqual(factor, amounts[pools.index(pool, 0)]);
    }
}

test "REDIST adsorbed anions require exact full fraction" {
    var amounts: [adsorbed_anion_pool_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts_mol = &amounts };
    try transfer(0, 1, 0.999999999, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST adsorbed anion runtime dimensions are exact" {
    var amounts: [adsorbed_anion_pool_count * 2 - 1]f64 = @splat(0);
    try std.testing.expectError(
        error.PondAdsorbedAnionTransferDimensionMismatch,
        transfer(0, 1, 1, .{ .layer_count = 2, .amounts_mol = &amounts }),
    );
}

test "REDIST adsorbed anion validation is atomic" {
    var amounts: [adsorbed_anion_pool_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts_mol = &amounts };
    amounts[pools.index(.banded_dihydrogen_phosphate, 0)] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidPondAdsorbedAnionTransferInput, transfer(0, 1, 1, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.nonband_deprotonated_site, 1)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.banded_dihydrogen_phosphate, 1)]);
}
