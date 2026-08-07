const std = @import("std");

/// Exact REDIST 9697--9714 order: four N pools, then two P pools.
pub const NutrientSolutePool = enum(u8) {
    ammonium_mol_n, // ZNH4S
    ammonia_mol_n, // ZNH3S
    nitrate_mol_n, // ZNO3S
    nitrite_mol_n, // ZNO2S
    hydrogen_phosphate_mol_p, // H1PO4
    dihydrogen_phosphate_mol_p, // H2PO4
};

pub const nutrient_solute_pool_count = std.meta.fields(NutrientSolutePool).len;

pub const Pools = struct {
    layer_count: usize,
    /// Pool-major extensive amounts; units are mol N or mol P as named by the pool.
    amounts_mol_element: []f64,
};

fn index(pools: Pools, pool: NutrientSolutePool, layer: usize) usize {
    return @intFromEnum(pool) * pools.layer_count + layer;
}

/// Direct translation of REDIST 9697--9714 under the positive-`BKDS` soil gate.
pub fn transfer(
    source_layer: usize,
    destination_layer: usize,
    source_bulk_density_megagrams_m3: f64,
    destination_bulk_density_megagrams_m3: f64,
    water_fraction: f64,
    pools: Pools,
) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.amounts_mol_element.len != nutrient_solute_pool_count * pools.layer_count)
        return error.SoilNutrientSoluteTransferDimensionMismatch;
    if (!std.math.isFinite(source_bulk_density_megagrams_m3) or !std.math.isFinite(destination_bulk_density_megagrams_m3) or
        !std.math.isFinite(water_fraction) or water_fraction < 0 or water_fraction > 1)
        return error.InvalidSoilNutrientSoluteTransferInput;
    for (pools.amounts_mol_element) |amount| if (!std.math.isFinite(amount))
        return error.InvalidSoilNutrientSoluteTransferInput;
    if (source_bulk_density_megagrams_m3 <= 0 or destination_bulk_density_megagrams_m3 <= 0) return;

    inline for (std.meta.fields(NutrientSolutePool)) |field| {
        const pool: NutrientSolutePool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = water_fraction * pools.amounts_mol_element[source];
        if (!std.math.isFinite(moved) or !std.math.isFinite(pools.amounts_mol_element[destination] + moved) or
            !std.math.isFinite(pools.amounts_mol_element[source] - moved))
            return error.NonFiniteSoilNutrientSoluteTransferResult;
    }
    inline for (std.meta.fields(NutrientSolutePool)) |field| {
        const pool: NutrientSolutePool = @enumFromInt(field.value);
        const source = index(pools, pool, source_layer);
        const destination = index(pools, pool, destination_layer);
        const moved = water_fraction * pools.amounts_mol_element[source];
        pools.amounts_mol_element[destination] = pools.amounts_mol_element[destination] + moved;
        pools.amounts_mol_element[source] = pools.amounts_mol_element[source] - moved;
    }
}

test "REDIST soil nutrient solutes transfer exact N then P order" {
    var amounts: [nutrient_solute_pool_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .amounts_mol_element = &amounts };
    inline for (std.meta.fields(NutrientSolutePool), 0..) |field, ordinal| {
        const pool: NutrientSolutePool = @enumFromInt(field.value);
        amounts[index(pools, pool, 0)] = @floatFromInt(4 * (ordinal + 1));
        amounts[index(pools, pool, 1)] = 10;
    }
    try transfer(0, 1, 1, 1, 0.25, pools);
    inline for (std.meta.fields(NutrientSolutePool), 0..) |field, ordinal| {
        const pool: NutrientSolutePool = @enumFromInt(field.value);
        const moved: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(3 * moved, amounts[index(pools, pool, 0)]);
        try std.testing.expectEqual(10 + moved, amounts[index(pools, pool, 1)]);
    }
}

test "REDIST soil nutrient solutes conserve every element pool" {
    var amounts: [nutrient_solute_pool_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .amounts_mol_element = &amounts };
    inline for (std.meta.fields(NutrientSolutePool), 0..) |field, ordinal| {
        const pool: NutrientSolutePool = @enumFromInt(field.value);
        amounts[index(pools, pool, 0)] = @floatFromInt(ordinal + 1);
        amounts[index(pools, pool, 1)] = 10;
    }
    try transfer(0, 1, 1, 1, 0.4, pools);
    inline for (std.meta.fields(NutrientSolutePool), 0..) |field, ordinal| {
        const pool: NutrientSolutePool = @enumFromInt(field.value);
        try std.testing.expectApproxEqAbs(@as(f64, @floatFromInt(ordinal + 11)), amounts[index(pools, pool, 0)] + amounts[index(pools, pool, 1)], 1e-14);
    }
}

test "REDIST soil nutrient solute validation is atomic" {
    var amounts: [nutrient_solute_pool_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .amounts_mol_element = &amounts };
    amounts[index(pools, .dihydrogen_phosphate_mol_p, 0)] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilNutrientSoluteTransferInput, transfer(0, 1, 1, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .ammonium_mol_n, 0)]);
    try std.testing.expectEqual(@as(f64, 1), amounts[index(pools, .ammonium_mol_n, 1)]);
}
