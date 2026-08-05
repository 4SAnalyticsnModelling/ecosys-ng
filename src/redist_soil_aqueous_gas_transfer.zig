const std = @import("std");

/// Exact REDIST 10089--10106 order. Each value is an extensive elemental mass.
pub const AqueousGasPool = enum(u8) {
    carbon_dioxide_g_c, // CO2S
    methane_g_c, // CH4S
    oxygen_g_o, // OXYS
    dinitrogen_g_n, // Z2GS
    nitrous_oxide_g_n, // Z2OS
    hydrogen_g_h, // H2GS
};
pub const aqueous_gas_pool_count = std.meta.fields(AqueousGasPool).len;

pub const Pools = struct {
    layer_count: usize,
    /// Pool-major amounts in g C, g O, g N, or g H as named.
    amounts_g_element: []f64,
};
fn index(pools: Pools, pool: AqueousGasPool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}

/// Direct translation of REDIST 10088--10106 under the exact `L0 != 0` gate.
pub fn transfer(source_layer: usize, destination_layer: usize, water_fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.amounts_g_element.len != aqueous_gas_pool_count * pools.layer_count)
        return error.SoilAqueousGasTransferDimensionMismatch;
    if (!std.math.isFinite(water_fraction) or water_fraction < 0 or water_fraction > 1)
        return error.InvalidSoilAqueousGasTransferInput;
    for (pools.amounts_g_element) |amount| if (!std.math.isFinite(amount)) return error.InvalidSoilAqueousGasTransferInput;
    if (source_layer == 0) return;

    inline for (std.meta.fields(AqueousGasPool)) |field| {
        const pool: AqueousGasPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = water_fraction * pools.amounts_g_element[source];
        if (!std.math.isFinite(moved) or !std.math.isFinite(pools.amounts_g_element[destination] + moved) or
            !std.math.isFinite(pools.amounts_g_element[source] - moved)) return error.NonFiniteSoilAqueousGasTransferResult;
    }
    inline for (std.meta.fields(AqueousGasPool)) |field| {
        const pool: AqueousGasPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = water_fraction * pools.amounts_g_element[source];
        pools.amounts_g_element[destination] = pools.amounts_g_element[destination] + moved;
        pools.amounts_g_element[source] = pools.amounts_g_element[source] - moved;
    }
}

test "REDIST aqueous gases transfer exact CO2-through-H2 order" {
    var amounts: [aqueous_gas_pool_count * 3]f64 = @splat(0);
    const pools = Pools{ .layer_count = 3, .amounts_g_element = &amounts };
    inline for (std.meta.fields(AqueousGasPool), 0..) |field, ordinal| {
        const pool: AqueousGasPool = @enumFromInt(field.value);
        amounts[index(pools, pool, 2)] = @floatFromInt(4 * (ordinal + 1));
        amounts[index(pools, pool, 1)] = 10;
    }
    try transfer(2, 1, 0.25, pools);
    inline for (std.meta.fields(AqueousGasPool), 0..) |field, ordinal| {
        const pool: AqueousGasPool = @enumFromInt(field.value);
        const moved: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(3 * moved, amounts[index(pools, pool, 2)]);
        try std.testing.expectEqual(10 + moved, amounts[index(pools, pool, 1)]);
    }
}

test "REDIST aqueous gas L0 gate and conservation are exact" {
    var amounts: [aqueous_gas_pool_count * 3]f64 = @splat(2);
    const pools = Pools{ .layer_count = 3, .amounts_g_element = &amounts };
    try transfer(0, 1, 1, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 2), amount);
    try transfer(2, 1, 0.4, pools);
    inline for (std.meta.fields(AqueousGasPool)) |field| {
        const pool: AqueousGasPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 4), amounts[index(pools, pool, 2)] + amounts[index(pools, pool, 1)]);
    }
}

test "REDIST aqueous gas validation is atomic" {
    var amounts: [aqueous_gas_pool_count * 3]f64 = @splat(1);
    const pools = Pools{ .layer_count = 3, .amounts_g_element = &amounts };
    amounts[index(pools, .hydrogen_g_h, 2)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilAqueousGasTransferInput, transfer(2, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .carbon_dioxide_g_c, 2)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .carbon_dioxide_g_c, 1)]);
}
