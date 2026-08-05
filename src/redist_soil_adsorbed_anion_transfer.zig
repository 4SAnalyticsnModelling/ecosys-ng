const std = @import("std");

/// Exact REDIST 9944--9976 order, including exchange capacity first.
pub const AdsorbedAnionPool = enum(u8) {
    anion_exchange_capacity_mol_charge, // XAEC
    nonband_deprotonated_site_mol, // XOH0
    nonband_hydroxyl_site_mol, // XOH1
    nonband_protonated_site_mol, // XOH2
    nonband_hydrogen_phosphate_mol, // XH1P
    nonband_dihydrogen_phosphate_mol, // XH2P
    banded_deprotonated_site_mol, // XOH0B
    banded_hydroxyl_site_mol, // XOH1B
    banded_protonated_site_mol, // XOH2B
    banded_hydrogen_phosphate_mol, // XH1PB
    banded_dihydrogen_phosphate_mol, // XH2PB
};

pub const adsorbed_anion_pool_count = std.meta.fields(AdsorbedAnionPool).len;

pub const Pools = struct {
    layer_count: usize,
    /// Pool-major extensive amounts, mol charge or mol as named.
    amounts: []f64,
};

fn index(pools: Pools, pool: AdsorbedAnionPool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}

/// Direct translation of REDIST 9943--9976 inside the `L0 > L1` branch.
pub fn transfer(
    source_layer: usize,
    destination_layer: usize,
    source_bulk_density_megagrams_m3: f64,
    destination_bulk_density_megagrams_m3: f64,
    redistribution_fraction: f64,
    pools: Pools,
) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.amounts.len != adsorbed_anion_pool_count * pools.layer_count)
        return error.SoilAdsorbedAnionTransferDimensionMismatch;
    if (!std.math.isFinite(source_bulk_density_megagrams_m3) or !std.math.isFinite(destination_bulk_density_megagrams_m3) or
        !std.math.isFinite(redistribution_fraction) or redistribution_fraction < 0 or redistribution_fraction > 1)
        return error.InvalidSoilAdsorbedAnionTransferInput;
    for (pools.amounts) |amount| if (!std.math.isFinite(amount))
        return error.InvalidSoilAdsorbedAnionTransferInput;
    if (source_bulk_density_megagrams_m3 <= 0 or destination_bulk_density_megagrams_m3 <= 0 or source_layer <= destination_layer) return;

    const anion_fraction = redistribution_fraction; // FAO=FX.
    inline for (std.meta.fields(AdsorbedAnionPool)) |field| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = anion_fraction * pools.amounts[source];
        if (!std.math.isFinite(moved) or !std.math.isFinite(pools.amounts[destination] + moved) or
            !std.math.isFinite(pools.amounts[source] - moved))
            return error.NonFiniteSoilAdsorbedAnionTransferResult;
    }
    inline for (std.meta.fields(AdsorbedAnionPool)) |field| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = anion_fraction * pools.amounts[source];
        pools.amounts[destination] = pools.amounts[destination] + moved;
        pools.amounts[source] = pools.amounts[source] - moved;
    }
}

test "REDIST adsorbed anions transfer exact XAEC-through-XH2PB order" {
    var amounts: [adsorbed_anion_pool_count * 3]f64 = @splat(0);
    const pools = Pools{ .layer_count = 3, .amounts = &amounts };
    inline for (std.meta.fields(AdsorbedAnionPool), 0..) |field, ordinal| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        amounts[index(pools, pool, 2)] = @floatFromInt(4 * (ordinal + 1));
        amounts[index(pools, pool, 1)] = 10;
    }
    try transfer(2, 1, 1, 1, 0.25, pools);
    inline for (std.meta.fields(AdsorbedAnionPool), 0..) |field, ordinal| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        const moved: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(3 * moved, amounts[index(pools, pool, 2)]);
        try std.testing.expectEqual(10 + moved, amounts[index(pools, pool, 1)]);
    }
}

test "REDIST adsorbed anion direction gate and conservation are exact" {
    var amounts: [adsorbed_anion_pool_count * 3]f64 = @splat(2);
    const pools = Pools{ .layer_count = 3, .amounts = &amounts };
    try transfer(1, 2, 1, 1, 1, pools);
    for (amounts) |amount| try std.testing.expectEqual(@as(f64, 2), amount);
    try transfer(2, 1, 1, 1, 0.4, pools);
    inline for (std.meta.fields(AdsorbedAnionPool)) |field| {
        const pool: AdsorbedAnionPool = @enumFromInt(field.value);
        try std.testing.expectEqual(@as(f64, 4), amounts[index(pools, pool, 2)] + amounts[index(pools, pool, 1)]);
    }
}

test "REDIST adsorbed anion validation is atomic" {
    var amounts: [adsorbed_anion_pool_count * 3]f64 = @splat(1);
    const pools = Pools{ .layer_count = 3, .amounts = &amounts };
    amounts[index(pools, .banded_dihydrogen_phosphate_mol, 2)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilAdsorbedAnionTransferInput, transfer(2, 1, 1, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .anion_exchange_capacity_mol_charge, 2)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .anion_exchange_capacity_mol_charge, 1)]);
}
