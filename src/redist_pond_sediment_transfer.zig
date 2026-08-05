const std = @import("std");

/// Exact REDIST 8905--8962 sediment assignment order.
pub const SedimentPool = enum(u8) {
    sand,
    silt,
    clay,
    cation_exchange_capacity,
    anion_exchange_capacity,
    exchangeable_hydrogen,
    exchangeable_aluminum,
    exchangeable_iron,
    exchangeable_calcium,
    exchangeable_magnesium,
    exchangeable_sodium,
    exchangeable_potassium,
    exchangeable_hydroxyl,
    aluminum_hydroxide_precipitate,
    iron_hydroxide_precipitate,
    calcium_carbonate_precipitate,
    calcium_sulfate_precipitate,
    aluminum_silicate,
    iron_silicate,
    calcium_silicate,
    magnesium_silicate,
    sodium_silicate,
    potassium_silicate,
    fertilizer_aluminum_silicate,
    fertilizer_iron_silicate,
    fertilizer_calcium_silicate,
    fertilizer_magnesium_silicate,
    fertilizer_sodium_silicate,
    fertilizer_potassium_silicate,
};

pub const sediment_pool_count = std.meta.fields(SedimentPool).len;
pub const Unit = enum { Mg, mol };

pub fn unit(pool: SedimentPool) Unit {
    return switch (pool) {
        .sand, .silt, .clay => .Mg,
        else => .mol,
    };
}

pub const Pools = struct {
    layer_count: usize,
    /// Pool-major: `pool * layer_count + layer`.
    amounts: []f64,

    fn index(self: Pools, pool: SedimentPool, layer: usize) usize {
        return @intFromEnum(pool) * self.layer_count + layer;
    }
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8903--8963: full-fraction pond sediment.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.amounts.len != sediment_pool_count * pools.layer_count)
        return error.PondSedimentTransferDimensionMismatch;
    if (!finiteSlice(pools.amounts) or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondSedimentTransferInput;
    if (fraction != 1.0 or source_layer == 0) return;

    inline for (std.meta.fields(SedimentPool)) |field| {
        const pool: SedimentPool = @enumFromInt(field.value);
        const source = pools.index(pool, source_layer);
        const destination = pools.index(pool, destination_layer);
        if (!std.math.isFinite(pools.amounts[destination] + fraction * pools.amounts[source]))
            return error.NonFinitePondSedimentTransferResult;
    }
    inline for (std.meta.fields(SedimentPool)) |field| {
        const pool: SedimentPool = @enumFromInt(field.value);
        const source = pools.index(pool, source_layer);
        const destination = pools.index(pool, destination_layer);
        pools.amounts[destination] = pools.amounts[destination] + fraction * pools.amounts[source];
    }
}

test "REDIST full-fraction pond sediment transfers all pools in exact order" {
    var amounts: [sediment_pool_count * 3]f64 = @splat(0);
    const pools = Pools{ .layer_count = 3, .amounts = &amounts };
    inline for (std.meta.fields(SedimentPool), 0..) |field, ordinal| {
        const pool: SedimentPool = @enumFromInt(field.value);
        amounts[pools.index(pool, 1)] = @floatFromInt(ordinal + 1);
        amounts[pools.index(pool, 2)] = @floatFromInt((ordinal + 1) * 10);
    }
    try transfer(1, 2, 1.0, pools);
    inline for (std.meta.fields(SedimentPool), 0..) |field, ordinal| {
        const pool: SedimentPool = @enumFromInt(field.value);
        const factor: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(11 * factor, amounts[pools.index(pool, 2)]);
        try std.testing.expectEqual(factor, amounts[pools.index(pool, 1)]);
    }
}

test "REDIST pond sediment requires exact full fraction" {
    var amounts: [sediment_pool_count * 3]f64 = @splat(1);
    const pools = Pools{ .layer_count = 3, .amounts = &amounts };
    try transfer(1, 2, 0.999999999, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST pond sediment excludes source layer zero" {
    var amounts: [sediment_pool_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts = &amounts };
    try transfer(0, 1, 1, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 1), amount);
}

test "REDIST pond sediment units preserve Mg capitalization" {
    try std.testing.expectEqual(Unit.Mg, unit(.sand));
    try std.testing.expectEqual(Unit.Mg, unit(.clay));
    try std.testing.expectEqual(Unit.mol, unit(.cation_exchange_capacity));
    try std.testing.expectEqual(Unit.mol, unit(.fertilizer_potassium_silicate));
}

test "REDIST pond sediment runtime dimensions are exact" {
    var amounts: [sediment_pool_count * 2 - 1]f64 = @splat(0);
    try std.testing.expectError(
        error.PondSedimentTransferDimensionMismatch,
        transfer(0, 1, 1, .{ .layer_count = 2, .amounts = &amounts }),
    );
}

test "REDIST pond sediment validation is atomic" {
    var amounts: [sediment_pool_count * 3]f64 = @splat(1);
    const pools = Pools{ .layer_count = 3, .amounts = &amounts };
    amounts[pools.index(.fertilizer_potassium_silicate, 1)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSedimentTransferInput, transfer(1, 2, 1, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.sand, 2)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[pools.index(.fertilizer_potassium_silicate, 2)]);
}
