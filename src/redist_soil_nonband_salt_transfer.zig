const std = @import("std");
const extended_definition = @import("redist_pond_extended_salt_transfer.zig");

pub const BaseSpecies = enum(u8) { hydrogen, hydroxide, aluminum, iron, calcium, magnesium, sodium, potassium };
pub const base_species_count = std.meta.fields(BaseSpecies).len;
pub const ExtendedSpecies = extended_definition.Species;
pub const extended_species_count = extended_definition.species_count;
pub const ExtendedUnit = extended_definition.Unit;
pub const extendedUnit = extended_definition.unit;

pub const Pools = struct {
    layer_count: usize,
    /// Base-species-major extensive amounts, mol.
    base_amounts_mol: []f64,
    /// Extended-species-major amounts: mol except H0PO4/H3PO4, which are g P.
    extended_amounts: []f64,
};

fn baseIndex(pools: Pools, species: BaseSpecies, layer: usize) usize {
    return @intFromEnum(species) * pools.layer_count + layer;
}
fn extendedIndex(pools: Pools, species: ExtendedSpecies, layer: usize) usize {
    return @intFromEnum(species) * pools.layer_count + layer;
}
fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}
fn validateMove(source: f64, destination: f64, water_fraction: f64) !void {
    const moved = water_fraction * source;
    if (!std.math.isFinite(moved) or !std.math.isFinite(destination + moved) or !std.math.isFinite(source - moved))
        return error.NonFiniteSoilNonbandSaltTransferResult;
}
fn move(values: []f64, source: usize, destination: usize, water_fraction: f64) void {
    const moved = water_fraction * values[source];
    values[destination] = values[destination] + moved;
    values[source] = values[source] - moved;
}

/// Direct translation of REDIST 9718--9845 under the positive-`BKDS` soil gate.
pub fn transfer(
    salts_enabled: bool,
    source_layer: usize,
    destination_layer: usize,
    source_bulk_density_megagrams_m3: f64,
    destination_bulk_density_megagrams_m3: f64,
    water_fraction: f64,
    pools: Pools,
) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.base_amounts_mol.len != base_species_count * pools.layer_count or
        pools.extended_amounts.len != extended_species_count * pools.layer_count)
        return error.SoilNonbandSaltTransferDimensionMismatch;
    if (!finiteSlice(pools.base_amounts_mol) or !finiteSlice(pools.extended_amounts) or
        !std.math.isFinite(source_bulk_density_megagrams_m3) or !std.math.isFinite(destination_bulk_density_megagrams_m3) or
        !std.math.isFinite(water_fraction) or water_fraction < 0 or water_fraction > 1)
        return error.InvalidSoilNonbandSaltTransferInput;
    if (source_bulk_density_megagrams_m3 <= 0 or destination_bulk_density_megagrams_m3 <= 0) return;

    inline for (std.meta.fields(BaseSpecies)) |field| {
        const species: BaseSpecies = @enumFromInt(field.value);
        try validateMove(pools.base_amounts_mol[baseIndex(pools, species, source_layer)], pools.base_amounts_mol[baseIndex(pools, species, destination_layer)], water_fraction);
    }
    if (salts_enabled) inline for (std.meta.fields(ExtendedSpecies)) |field| {
        const species: ExtendedSpecies = @enumFromInt(field.value);
        try validateMove(pools.extended_amounts[extendedIndex(pools, species, source_layer)], pools.extended_amounts[extendedIndex(pools, species, destination_layer)], water_fraction);
    };
    inline for (std.meta.fields(BaseSpecies)) |field| {
        const species: BaseSpecies = @enumFromInt(field.value);
        move(pools.base_amounts_mol, baseIndex(pools, species, source_layer), baseIndex(pools, species, destination_layer), water_fraction);
    }
    if (salts_enabled) inline for (std.meta.fields(ExtendedSpecies)) |field| {
        const species: ExtendedSpecies = @enumFromInt(field.value);
        move(pools.extended_amounts, extendedIndex(pools, species, source_layer), extendedIndex(pools, species, destination_layer), water_fraction);
    };
}

test "REDIST non-band salts transfer base order and obey ISALTG" {
    var base: [base_species_count * 2]f64 = @splat(0);
    var extended: [extended_species_count * 2]f64 = @splat(4);
    const pools = Pools{ .layer_count = 2, .base_amounts_mol = &base, .extended_amounts = &extended };
    inline for (std.meta.fields(BaseSpecies), 0..) |field, ordinal| {
        const species: BaseSpecies = @enumFromInt(field.value);
        base[baseIndex(pools, species, 0)] = @floatFromInt(4 * (ordinal + 1));
        base[baseIndex(pools, species, 1)] = 10;
    }
    try transfer(false, 0, 1, 1, 1, 0.25, pools);
    inline for (std.meta.fields(BaseSpecies), 0..) |field, ordinal| {
        const species: BaseSpecies = @enumFromInt(field.value);
        const moved: f64 = @floatFromInt(ordinal + 1);
        try std.testing.expectEqual(3 * moved, base[baseIndex(pools, species, 0)]);
        try std.testing.expectEqual(10 + moved, base[baseIndex(pools, species, 1)]);
    }
    for (extended) |amount| try std.testing.expectEqual(@as(f64, 4), amount);
}

test "REDIST extended non-band salts preserve first-to-last order and conservation" {
    var base: [base_species_count * 2]f64 = @splat(2);
    var extended: [extended_species_count * 2]f64 = @splat(0);
    const pools = Pools{ .layer_count = 2, .base_amounts_mol = &base, .extended_amounts = &extended };
    inline for (std.meta.fields(ExtendedSpecies), 0..) |field, ordinal| {
        const species: ExtendedSpecies = @enumFromInt(field.value);
        extended[extendedIndex(pools, species, 0)] = @floatFromInt(ordinal + 1);
        extended[extendedIndex(pools, species, 1)] = 10;
    }
    try transfer(true, 0, 1, 1, 1, 0.4, pools);
    inline for (std.meta.fields(ExtendedSpecies), 0..) |field, ordinal| {
        const species: ExtendedSpecies = @enumFromInt(field.value);
        try std.testing.expectApproxEqAbs(@as(f64, @floatFromInt(ordinal + 11)), extended[extendedIndex(pools, species, 0)] + extended[extendedIndex(pools, species, 1)], 1e-14);
    }
    try std.testing.expectEqual(ExtendedUnit.g_p, extendedUnit(.phosphate));
    try std.testing.expectEqual(ExtendedUnit.mol, extendedUnit(.magnesium_phosphate_1));
}

test "REDIST non-band salt validation is atomic" {
    var base: [base_species_count * 2]f64 = @splat(1);
    var extended: [extended_species_count * 2]f64 = @splat(1);
    const pools = Pools{ .layer_count = 2, .base_amounts_mol = &base, .extended_amounts = &extended };
    extended[extendedIndex(pools, .magnesium_phosphate_1, 0)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilNonbandSaltTransferInput, transfer(true, 0, 1, 1, 1, 0.5, pools));
    try std.testing.expectEqual(@as(f64, 1), base[baseIndex(pools, .hydrogen, 0)]);
    try std.testing.expectEqual(@as(f64, 1), base[baseIndex(pools, .hydrogen, 1)]);
}
