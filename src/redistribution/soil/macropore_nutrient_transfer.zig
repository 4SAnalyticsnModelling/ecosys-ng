const std = @import("std");

/// Exact REDIST 10111--10146 order; amounts are grams of N or P as named.
pub const MacroporeNutrientPool = enum(u8) {
    nonband_ammonium_g_n,
    nonband_ammonia_g_n,
    nonband_nitrate_g_n,
    nonband_nitrite_g_n,
    banded_ammonium_g_n,
    banded_ammonia_g_n,
    banded_nitrate_g_n,
    banded_nitrite_g_n,
    nonband_hydrogen_phosphate_g_p,
    nonband_dihydrogen_phosphate_g_p,
    banded_hydrogen_phosphate_g_p,
    banded_dihydrogen_phosphate_g_p,
};
pub const macropore_nutrient_pool_count = std.meta.fields(MacroporeNutrientPool).len;

pub const Pools = struct {
    layer_count: usize,
    /// Pool-major extensive amounts, g N or g P as named.
    amounts_g_element: []f64,
    /// FHOL by layer, dimensionless preferential-domain fraction.
    preferential_fraction: []const f64,
};
fn index(pools: Pools, pool: MacroporeNutrientPool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}

/// Direct translation of REDIST 10110--10146 inside the exact `L0 != 0` gate.
pub fn transfer(source_layer: usize, destination_layer: usize, preferential_transfer_fraction: f64, zero_tolerance: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.amounts_g_element.len != macropore_nutrient_pool_count * pools.layer_count or
        pools.preferential_fraction.len != pools.layer_count)
        return error.SoilMacroporeNutrientTransferDimensionMismatch;
    if (!std.math.isFinite(preferential_transfer_fraction) or preferential_transfer_fraction < 0 or preferential_transfer_fraction > 1 or
        !std.math.isFinite(zero_tolerance) or zero_tolerance < 0)
        return error.InvalidSoilMacroporeNutrientTransferInput;
    for (pools.amounts_g_element) |amount| if (!std.math.isFinite(amount)) return error.InvalidSoilMacroporeNutrientTransferInput;
    for (pools.preferential_fraction) |fraction| if (!std.math.isFinite(fraction)) return error.InvalidSoilMacroporeNutrientTransferInput;
    if (source_layer == 0 or pools.preferential_fraction[destination_layer] <= zero_tolerance or
        pools.preferential_fraction[source_layer] <= zero_tolerance) return;

    inline for (std.meta.fields(MacroporeNutrientPool)) |field| {
        const pool: MacroporeNutrientPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = preferential_transfer_fraction * pools.amounts_g_element[source];
        if (!std.math.isFinite(moved) or !std.math.isFinite(pools.amounts_g_element[destination] + moved) or
            !std.math.isFinite(pools.amounts_g_element[source] - moved)) return error.NonFiniteSoilMacroporeNutrientTransferResult;
    }
    inline for (std.meta.fields(MacroporeNutrientPool)) |field| {
        const pool: MacroporeNutrientPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = preferential_transfer_fraction * pools.amounts_g_element[source];
        pools.amounts_g_element[destination] = pools.amounts_g_element[destination] + moved;
        pools.amounts_g_element[source] = pools.amounts_g_element[source] - moved;
    }
}

test "REDIST macropore nutrients transfer exact nonband band and phosphate order" {
    var amounts: [macropore_nutrient_pool_count * 3]f64 = @splat(0);
    const fractions = [_]f64{ 0, 0.2, 0.3 };
    const pools = Pools{ .layer_count = 3, .amounts_g_element = &amounts, .preferential_fraction = &fractions };
    inline for (std.meta.fields(MacroporeNutrientPool), 0..) |field, ordinal| {
        const pool: MacroporeNutrientPool = @enumFromInt(field.value);
        amounts[index(pools, pool, 2)] = @floatFromInt(4 * (ordinal + 1));
        amounts[index(pools, pool, 1)] = 10;
    }
    try transfer(2, 1, 0.25, 0, pools);
    inline for (std.meta.fields(MacroporeNutrientPool), 0..) |field, ordinal| {
        const pool: MacroporeNutrientPool = @enumFromInt(field.value);
        const moved: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(3 * moved, amounts[index(pools, pool, 2)]);
        try std.testing.expectEqual(10 + moved, amounts[index(pools, pool, 1)]);
    }
}

test "REDIST macropore nutrient FHOL gates and conservation are exact" {
    var amounts: [macropore_nutrient_pool_count * 3]f64 = @splat(2);
    const fractions = [_]f64{ 0, 0, 0.3 };
    const pools = Pools{ .layer_count = 3, .amounts_g_element = &amounts, .preferential_fraction = &fractions };
    try transfer(2, 1, 1, 0, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 2), amount);
}

test "REDIST macropore nutrient validation is atomic" {
    var amounts: [macropore_nutrient_pool_count * 3]f64 = @splat(1);
    const fractions = [_]f64{ 0, 0.2, 0.3 };
    const pools = Pools{ .layer_count = 3, .amounts_g_element = &amounts, .preferential_fraction = &fractions };
    amounts[index(pools, .banded_dihydrogen_phosphate_g_p, 2)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilMacroporeNutrientTransferInput, transfer(2, 1, 0.5, 0, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .nonband_ammonium_g_n, 2)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .nonband_ammonium_g_n, 1)]);
}
