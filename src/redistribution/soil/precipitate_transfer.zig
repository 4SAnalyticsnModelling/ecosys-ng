const std = @import("std");

/// Exact REDIST 9981--10058 order.
pub const PrecipitatePool = enum(u8) {
    aluminum_hydroxide,
    iron_hydroxide,
    calcium_carbonate,
    calcium_sulfate,
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
    nonband_aluminum_phosphate,
    nonband_iron_phosphate,
    nonband_dicalcium_phosphate,
    nonband_hydroxyapatite,
    nonband_monocalcium_phosphate,
    banded_aluminum_phosphate,
    banded_iron_phosphate,
    banded_dicalcium_phosphate,
    banded_hydroxyapatite,
    banded_monocalcium_phosphate,
};

pub const precipitate_pool_count = std.meta.fields(PrecipitatePool).len;
pub const Unit = enum { mol, mol_p };
pub fn unit(pool: PrecipitatePool) Unit {
    return switch (pool) {
        .nonband_aluminum_phosphate,
        .nonband_iron_phosphate,
        .nonband_dicalcium_phosphate,
        .nonband_hydroxyapatite,
        .nonband_monocalcium_phosphate,
        .banded_aluminum_phosphate,
        .banded_iron_phosphate,
        .banded_dicalcium_phosphate,
        .banded_hydroxyapatite,
        .banded_monocalcium_phosphate,
        => .mol_p,
        else => .mol,
    };
}

pub const Pools = struct {
    layer_count: usize,
    /// Pool-major extensive amounts, in the unit returned by `unit`.
    amounts: []f64,
};
fn index(pools: Pools, pool: PrecipitatePool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}

/// Direct translation of REDIST 9980--10059 inside the `L0 > L1` branch.
pub fn transfer(source_layer: usize, destination_layer: usize, source_bulk_density_megagrams_m3: f64, destination_bulk_density_megagrams_m3: f64, redistribution_fraction: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.amounts.len != precipitate_pool_count * pools.layer_count)
        return error.SoilPrecipitateTransferDimensionMismatch;
    if (!std.math.isFinite(source_bulk_density_megagrams_m3) or !std.math.isFinite(destination_bulk_density_megagrams_m3) or
        !std.math.isFinite(redistribution_fraction) or redistribution_fraction < 0 or redistribution_fraction > 1)
        return error.InvalidSoilPrecipitateTransferInput;
    for (pools.amounts) |amount| if (!std.math.isFinite(amount)) return error.InvalidSoilPrecipitateTransferInput;
    if (source_bulk_density_megagrams_m3 <= 0 or destination_bulk_density_megagrams_m3 <= 0 or source_layer <= destination_layer) return;

    const precipitate_fraction = redistribution_fraction; // FPO=FX.
    inline for (std.meta.fields(PrecipitatePool)) |field| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = precipitate_fraction * pools.amounts[source];
        if (!std.math.isFinite(moved) or !std.math.isFinite(pools.amounts[destination] + moved) or
            !std.math.isFinite(pools.amounts[source] - moved)) return error.NonFiniteSoilPrecipitateTransferResult;
    }
    inline for (std.meta.fields(PrecipitatePool)) |field| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = precipitate_fraction * pools.amounts[source];
        pools.amounts[destination] = pools.amounts[destination] + moved;
        pools.amounts[source] = pools.amounts[source] - moved;
    }
}

test "REDIST soil precipitates transfer exact PALOH-through-PCPMB order" {
    var amounts: [precipitate_pool_count * 3]f64 = @splat(0);
    const pools = Pools{ .layer_count = 3, .amounts = &amounts };
    inline for (std.meta.fields(PrecipitatePool), 0..) |field, ordinal| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        amounts[index(pools, pool, 2)] = @floatFromInt(4 * (ordinal + 1));
        amounts[index(pools, pool, 1)] = 10;
    }
    try transfer(2, 1, 1, 1, 0.25, pools);
    inline for (std.meta.fields(PrecipitatePool), 0..) |field, ordinal| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        const moved: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(3 * moved, amounts[index(pools, pool, 2)]);
        try std.testing.expectEqual(10 + moved, amounts[index(pools, pool, 1)]);
    }
    try std.testing.expectEqual(Unit.mol, unit(.fertilizer_potassium_silicate));
    try std.testing.expectEqual(Unit.mol_p, unit(.banded_monocalcium_phosphate));
}

test "REDIST soil precipitate direction gate and conservation are exact" {
    var amounts: [precipitate_pool_count * 3]f64 = @splat(2);
    const pools = Pools{ .layer_count = 3, .amounts = &amounts };
    try transfer(1, 2, 1, 1, 1, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 2), amount);
    try transfer(2, 1, 1, 1, 0.4, pools);
    inline for (std.meta.fields(PrecipitatePool)) |field| {
        const pool: PrecipitatePool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 4), amounts[index(pools, pool, 2)] + amounts[index(pools, pool, 1)]);
    }
}

test "REDIST soil precipitate validation is atomic" {
    var amounts: [precipitate_pool_count * 3]f64 = @splat(1);
    const pools = Pools{ .layer_count = 3, .amounts = &amounts };
    amounts[index(pools, .banded_monocalcium_phosphate, 2)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilPrecipitateTransferInput, transfer(2, 1, 1, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .aluminum_hydroxide, 2)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .aluminum_hydroxide, 1)]);
}
