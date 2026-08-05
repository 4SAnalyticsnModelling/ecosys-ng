const std = @import("std");

/// Exact REDIST 10063--10083 order. Each value is an extensive elemental mass.
pub const GaseousPool = enum(u8) {
    carbon_dioxide_g_c, // CO2G
    methane_g_c, // CH4G
    oxygen_g_o, // OXYG
    dinitrogen_g_n, // Z2GG
    nitrous_oxide_g_n, // Z2OG
    ammonia_g_n, // ZNH3G
    hydrogen_g_h, // H2GG
};
pub const gaseous_pool_count = std.meta.fields(GaseousPool).len;

pub const Pools = struct {
    layer_count: usize,
    /// Pool-major amounts in g C, g O, g N, or g H as named.
    amounts_g_element: []f64,
};
fn index(pools: Pools, pool: GaseousPool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}

/// Direct translation of REDIST 10063--10084 under the positive-`BKDS` soil gate.
pub fn transfer(source_layer: usize, destination_layer: usize, source_bulk_density_megagrams_m3: f64, destination_bulk_density_megagrams_m3: f64, water_fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.amounts_g_element.len != gaseous_pool_count * pools.layer_count)
        return error.SoilGaseousGasTransferDimensionMismatch;
    if (!std.math.isFinite(source_bulk_density_megagrams_m3) or !std.math.isFinite(destination_bulk_density_megagrams_m3) or
        !std.math.isFinite(water_fraction) or water_fraction < 0 or water_fraction > 1)
        return error.InvalidSoilGaseousGasTransferInput;
    for (pools.amounts_g_element) |amount| if (!std.math.isFinite(amount)) return error.InvalidSoilGaseousGasTransferInput;
    if (source_bulk_density_megagrams_m3 <= 0 or destination_bulk_density_megagrams_m3 <= 0) return;

    inline for (std.meta.fields(GaseousPool)) |field| {
        const pool: GaseousPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = water_fraction * pools.amounts_g_element[source];
        if (!std.math.isFinite(moved) or !std.math.isFinite(pools.amounts_g_element[destination] + moved) or
            !std.math.isFinite(pools.amounts_g_element[source] - moved)) return error.NonFiniteSoilGaseousGasTransferResult;
    }
    inline for (std.meta.fields(GaseousPool)) |field| {
        const pool: GaseousPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = water_fraction * pools.amounts_g_element[source];
        pools.amounts_g_element[destination] = pools.amounts_g_element[destination] + moved;
        pools.amounts_g_element[source] = pools.amounts_g_element[source] - moved;
    }
}

test "REDIST gaseous gases transfer exact CO2-through-H2 order" {
    var amounts: [gaseous_pool_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .amounts_g_element = &amounts };
    inline for (std.meta.fields(GaseousPool), 0..) |field, ordinal| {
        const pool: GaseousPool = @enumFromInt(field.value);
        amounts[index(pools, pool, 0)] = @floatFromInt(4 * (ordinal + 1));
        amounts[index(pools, pool, 1)] = 10;
    }
    try transfer(0, 1, 1, 1, 0.25, pools);
    inline for (std.meta.fields(GaseousPool), 0..) |field, ordinal| {
        const pool: GaseousPool = @enumFromInt(field.value);
        const moved: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(3 * moved, amounts[index(pools, pool, 0)]);
        try std.testing.expectEqual(10 + moved, amounts[index(pools, pool, 1)]);
    }
}

test "REDIST gaseous gas transfer conserves every elemental pool" {
    var amounts: [gaseous_pool_count * 2]f64 = @splat(2);
    const pools = Pools{ .layer_count = 2, .amounts_g_element = &amounts };
    try transfer(0, 1, 1, 1, 0.4, pools);
    inline for (std.meta.fields(GaseousPool)) |field| {
        const pool: GaseousPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 4), amounts[index(pools, pool, 0)] + amounts[index(pools, pool, 1)]);
    }
}

test "REDIST gaseous gas validation is atomic" {
    var amounts: [gaseous_pool_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts_g_element = &amounts };
    amounts[index(pools, .hydrogen_g_h, 0)] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilGaseousGasTransferInput, transfer(0, 1, 1, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .carbon_dioxide_g_c, 0)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .carbon_dioxide_g_c, 1)]);
}
