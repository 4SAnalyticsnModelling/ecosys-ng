const std = @import("std");

/// Exact REDIST 10302--10316 order; amounts are grams of the named element.
pub const MacroporeAqueousGasPool = enum(u8) {
    carbon_dioxide_g_c, // CO2SH
    methane_g_c, // CH4SH
    oxygen_g_o, // OXYSH
    dinitrogen_g_n, // Z2GSH
    nitrous_oxide_g_n, // Z2OSH
};
pub const pool_count = std.meta.fields(MacroporeAqueousGasPool).len;

pub const Pools = struct {
    layer_count: usize,
    amounts_g_element: []f64,
    preferential_fraction: []const f64, // FHOL, dimensionless.
};
fn index(pools: Pools, pool: MacroporeAqueousGasPool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}

/// Direct translation of REDIST 10302--10318 under `L0 != 0` and positive-`FHOL`.
pub fn transfer(source_layer: usize, destination_layer: usize, preferential_transfer_fraction: f64, zero_tolerance: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.amounts_g_element.len != pool_count * pools.layer_count or
        pools.preferential_fraction.len != pools.layer_count)
        return error.SoilMacroporeAqueousGasTransferDimensionMismatch;
    if (!std.math.isFinite(preferential_transfer_fraction) or preferential_transfer_fraction < 0 or preferential_transfer_fraction > 1 or
        !std.math.isFinite(zero_tolerance) or zero_tolerance < 0)
        return error.InvalidSoilMacroporeAqueousGasTransferInput;
    for (pools.amounts_g_element) |amount| if (!std.math.isFinite(amount)) return error.InvalidSoilMacroporeAqueousGasTransferInput;
    for (pools.preferential_fraction) |fraction| if (!std.math.isFinite(fraction)) return error.InvalidSoilMacroporeAqueousGasTransferInput;
    if (source_layer == 0 or pools.preferential_fraction[source_layer] <= zero_tolerance or
        pools.preferential_fraction[destination_layer] <= zero_tolerance) return;

    inline for (std.meta.fields(MacroporeAqueousGasPool)) |field| {
        const pool: MacroporeAqueousGasPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = preferential_transfer_fraction * pools.amounts_g_element[source];
        if (!std.math.isFinite(moved) or !std.math.isFinite(pools.amounts_g_element[destination] + moved) or
            !std.math.isFinite(pools.amounts_g_element[source] - moved)) return error.NonFiniteSoilMacroporeAqueousGasTransferResult;
    }
    inline for (std.meta.fields(MacroporeAqueousGasPool)) |field| {
        const pool: MacroporeAqueousGasPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = preferential_transfer_fraction * pools.amounts_g_element[source];
        pools.amounts_g_element[destination] = pools.amounts_g_element[destination] + moved;
        pools.amounts_g_element[source] = pools.amounts_g_element[source] - moved;
    }
}

test "REDIST macropore aqueous gases transfer exact CO2-through-N2O order" {
    var amounts: [pool_count * 3]f64 = @splat(0);
    const fractions = [_]f64{ 0, 0.2, 0.3 };
    const pools = Pools{ .layer_count = 3, .amounts_g_element = &amounts, .preferential_fraction = &fractions };
    inline for (std.meta.fields(MacroporeAqueousGasPool), 0..) |field, ordinal| {
        const pool: MacroporeAqueousGasPool = @enumFromInt(field.value);
        amounts[index(pools, pool, 2)] = @floatFromInt(4 * (ordinal + 1));
        amounts[index(pools, pool, 1)] = 10;
    }
    try transfer(2, 1, 0.25, 0, pools);
    inline for (std.meta.fields(MacroporeAqueousGasPool), 0..) |field, ordinal| {
        const pool: MacroporeAqueousGasPool = @enumFromInt(field.value);
        const moved: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(3 * moved, amounts[index(pools, pool, 2)]);
        try std.testing.expectEqual(10 + moved, amounts[index(pools, pool, 1)]);
    }
}

test "REDIST macropore aqueous gas gates and conservation are exact" {
    var amounts: [pool_count * 3]f64 = @splat(2);
    const fractions = [_]f64{ 0, 0.2, 0.3 };
    const pools = Pools{ .layer_count = 3, .amounts_g_element = &amounts, .preferential_fraction = &fractions };
    try transfer(2, 1, 0.4, 0, pools);
    inline for (std.meta.fields(MacroporeAqueousGasPool)) |field| {
        const pool: MacroporeAqueousGasPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 4), amounts[index(pools, pool, 2)] + amounts[index(pools, pool, 1)]);
    }
}

test "REDIST macropore aqueous gas validation is atomic" {
    var amounts: [pool_count * 3]f64 = @splat(1);
    const fractions = [_]f64{ 0, 0.2, 0.3 };
    const pools = Pools{ .layer_count = 3, .amounts_g_element = &amounts, .preferential_fraction = &fractions };
    amounts[index(pools, .nitrous_oxide_g_n, 2)] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilMacroporeAqueousGasTransferInput, transfer(2, 1, 0.5, 0, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .carbon_dioxide_g_c, 2)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .carbon_dioxide_g_c, 1)]);
}
