const std = @import("std");

/// Exact REDIST 9907--9940 order, including exchange capacity first.
pub const AdsorbedCationPool = enum(u8) {
    cation_exchange_capacity_mol_charge, // XCEC
    nonband_ammonium_mol_n, // XN4
    banded_ammonium_mol_n, // XNB
    hydrogen_mol, // XHY
    aluminum_mol, // XAL
    iron_mol, // XFE
    calcium_mol, // XCA
    magnesium_mol, // XMG
    sodium_mol, // XNA
    potassium_mol, // XKA
    bicarbonate_mol, // XHC
};

pub const adsorbed_cation_pool_count = std.meta.fields(AdsorbedCationPool).len;

pub const Pools = struct {
    layer_count: usize,
    /// Pool-major extensive amounts, mol charge, mol N, or mol as named.
    amounts: []f64,
};

fn index(pools: Pools, pool: AdsorbedCationPool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}

/// Direct translation of REDIST 9905--9940 under the positive-`BKDS` soil gate.
pub fn transfer(
    source_layer: usize,
    destination_layer: usize,
    source_bulk_density_megagrams_m3: f64,
    destination_bulk_density_megagrams_m3: f64,
    redistribution_fraction: f64,
    pools: Pools,
) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.amounts.len != adsorbed_cation_pool_count * pools.layer_count)
        return error.SoilAdsorbedCationTransferDimensionMismatch;
    if (!std.math.isFinite(source_bulk_density_megagrams_m3) or !std.math.isFinite(destination_bulk_density_megagrams_m3) or
        !std.math.isFinite(redistribution_fraction) or redistribution_fraction < 0 or redistribution_fraction > 1)
        return error.InvalidSoilAdsorbedCationTransferInput;
    for (pools.amounts) |amount| if (!std.math.isFinite(amount))
        return error.InvalidSoilAdsorbedCationTransferInput;
    if (source_bulk_density_megagrams_m3 <= 0 or destination_bulk_density_megagrams_m3 <= 0 or source_layer <= destination_layer) return;

    const cation_fraction = redistribution_fraction; // FCO=FX.
    inline for (std.meta.fields(AdsorbedCationPool)) |field| {
        const pool: AdsorbedCationPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = cation_fraction * pools.amounts[source];
        if (!std.math.isFinite(moved) or !std.math.isFinite(pools.amounts[destination] + moved) or
            !std.math.isFinite(pools.amounts[source] - moved))
            return error.NonFiniteSoilAdsorbedCationTransferResult;
    }
    inline for (std.meta.fields(AdsorbedCationPool)) |field| {
        const pool: AdsorbedCationPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = cation_fraction * pools.amounts[source];
        pools.amounts[destination] = pools.amounts[destination] + moved;
        pools.amounts[source] = pools.amounts[source] - moved;
    }
}

test "REDIST adsorbed cations transfer exact XCEC-through-XHC order" {
    var amounts: [adsorbed_cation_pool_count * 3]f64 = @splat(0);
    const pools = Pools{ .layer_count = 3, .amounts = &amounts };
    inline for (std.meta.fields(AdsorbedCationPool), 0..) |field, ordinal| {
        const pool: AdsorbedCationPool = @enumFromInt(field.value);
        amounts[index(pools, pool, 2)] = @floatFromInt(4 * (ordinal + 1));
        amounts[index(pools, pool, 1)] = 10;
    }
    try transfer(2, 1, 1, 1, 0.25, pools);
    inline for (std.meta.fields(AdsorbedCationPool), 0..) |field, ordinal| {
        const pool: AdsorbedCationPool = @enumFromInt(field.value);
        const moved: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(3 * moved, amounts[index(pools, pool, 2)]);
        try std.testing.expectEqual(10 + moved, amounts[index(pools, pool, 1)]);
    }
}

test "REDIST adsorbed cation direction gate and conservation are exact" {
    var amounts: [adsorbed_cation_pool_count * 3]f64 = @splat(2);
    const pools = Pools{ .layer_count = 3, .amounts = &amounts };
    try transfer(1, 2, 1, 1, 1, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 2), amount);
    try transfer(2, 1, 1, 1, 0.4, pools);
    inline for (std.meta.fields(AdsorbedCationPool)) |field| {
        const pool: AdsorbedCationPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 4), amounts[index(pools, pool, 2)] + amounts[index(pools, pool, 1)]);
    }
}

test "REDIST adsorbed cation validation is atomic" {
    var amounts: [adsorbed_cation_pool_count * 3]f64 = @splat(1);
    const pools = Pools{ .layer_count = 3, .amounts = &amounts };
    amounts[index(pools, .bicarbonate_mol, 2)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilAdsorbedCationTransferInput, transfer(2, 1, 1, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .cation_exchange_capacity_mol_charge, 2)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .cation_exchange_capacity_mol_charge, 1)]);
}
