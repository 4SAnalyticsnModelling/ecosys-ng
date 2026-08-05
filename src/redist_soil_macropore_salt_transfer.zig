const std = @import("std");
const extended_definition = @import("redist_pond_extended_salt_transfer.zig");

pub const BaseSpecies = enum(u8) { hydrogen, hydroxide, aluminum, iron, calcium, magnesium, sodium, potassium };
pub const base_species_count = std.meta.fields(BaseSpecies).len;
pub const ExtendedSpecies = extended_definition.Species;
pub const extended_species_count = extended_definition.species_count;
pub const BandedPhosphateSpecies = enum(u8) {
    phosphate,
    phosphoric_acid,
    iron_phosphate_1,
    iron_phosphate_2,
    calcium_phosphate_0,
    calcium_phosphate_1,
    calcium_phosphate_2,
    magnesium_phosphate_1,
};
pub const banded_phosphate_species_count = std.meta.fields(BandedPhosphateSpecies).len;

pub const Pools = struct {
    layer_count: usize,
    base_amounts_mol: []f64,
    extended_amounts_mol: []f64,
    banded_phosphate_amounts_mol: []f64,
    preferential_fraction: []const f64, // FHOL, dimensionless.
};
fn index(layer_count: usize, species: anytype, layer: usize) usize {
    return @intFromEnum(species) * layer_count + layer;
}
fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}
fn validateMove(values: []const f64, source: usize, destination: usize, fraction: f64) !void {
    const moved = fraction * values[source];
    if (!std.math.isFinite(moved) or !std.math.isFinite(values[destination] + moved) or
        !std.math.isFinite(values[source] - moved)) return error.NonFiniteSoilMacroporeSaltTransferResult;
}
fn move(values: []f64, source: usize, destination: usize, fraction: f64) void {
    const moved = fraction * values[source];
    values[destination] = values[destination] + moved;
    values[source] = values[source] - moved;
}

/// Direct translation of REDIST 10150--10298 inside the `L0 != 0` and positive-`FHOL` gates.
pub fn transfer(salts_enabled: bool, source_layer: usize, destination_layer: usize, preferential_transfer_fraction: f64, zero_tolerance: f64, pools: Pools) !void {
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.base_amounts_mol.len != base_species_count * pools.layer_count or
        pools.extended_amounts_mol.len != extended_species_count * pools.layer_count or
        pools.banded_phosphate_amounts_mol.len != banded_phosphate_species_count * pools.layer_count or
        pools.preferential_fraction.len != pools.layer_count)
        return error.SoilMacroporeSaltTransferDimensionMismatch;
    if (!finiteSlice(pools.base_amounts_mol) or !finiteSlice(pools.extended_amounts_mol) or
        !finiteSlice(pools.banded_phosphate_amounts_mol) or !finiteSlice(pools.preferential_fraction) or
        !std.math.isFinite(preferential_transfer_fraction) or preferential_transfer_fraction < 0 or preferential_transfer_fraction > 1 or
        !std.math.isFinite(zero_tolerance) or zero_tolerance < 0)
        return error.InvalidSoilMacroporeSaltTransferInput;
    if (source_layer == 0 or pools.preferential_fraction[source_layer] <= zero_tolerance or
        pools.preferential_fraction[destination_layer] <= zero_tolerance) return;

    inline for (std.meta.fields(BaseSpecies)) |field| {
        const species: BaseSpecies = @enumFromInt(field.value);
        try validateMove(pools.base_amounts_mol, index(pools.layer_count, species, source_layer), index(pools.layer_count, species, destination_layer), preferential_transfer_fraction);
    }
    if (salts_enabled) {
        inline for (std.meta.fields(ExtendedSpecies)) |field| {
            const species: ExtendedSpecies = @enumFromInt(field.value);
            try validateMove(pools.extended_amounts_mol, index(pools.layer_count, species, source_layer), index(pools.layer_count, species, destination_layer), preferential_transfer_fraction);
        }
        inline for (std.meta.fields(BandedPhosphateSpecies)) |field| {
            const species: BandedPhosphateSpecies = @enumFromInt(field.value);
            try validateMove(pools.banded_phosphate_amounts_mol, index(pools.layer_count, species, source_layer), index(pools.layer_count, species, destination_layer), preferential_transfer_fraction);
        }
    }
    inline for (std.meta.fields(BaseSpecies)) |field| {
        const species: BaseSpecies = @enumFromInt(field.value);
        move(pools.base_amounts_mol, index(pools.layer_count, species, source_layer), index(pools.layer_count, species, destination_layer), preferential_transfer_fraction);
    }
    if (salts_enabled) {
        inline for (std.meta.fields(ExtendedSpecies)) |field| {
            const species: ExtendedSpecies = @enumFromInt(field.value);
            move(pools.extended_amounts_mol, index(pools.layer_count, species, source_layer), index(pools.layer_count, species, destination_layer), preferential_transfer_fraction);
        }
        inline for (std.meta.fields(BandedPhosphateSpecies)) |field| {
            const species: BandedPhosphateSpecies = @enumFromInt(field.value);
            move(pools.banded_phosphate_amounts_mol, index(pools.layer_count, species, source_layer), index(pools.layer_count, species, destination_layer), preferential_transfer_fraction);
        }
    }
}

const Fixture = struct {
    base: [base_species_count * 3]f64 = @splat(2),
    extended: [extended_species_count * 3]f64 = @splat(2),
    banded: [banded_phosphate_species_count * 3]f64 = @splat(2),
    fractions: [3]f64 = .{ 0, 0.2, 0.3 },
    fn pools(self: *Fixture) Pools {
        return .{ .layer_count = 3, .base_amounts_mol = &self.base, .extended_amounts_mol = &self.extended, .banded_phosphate_amounts_mol = &self.banded, .preferential_fraction = &self.fractions };
    }
};

test "REDIST macropore salts transfer exact base extended and banded order" {
    var fixture = Fixture{};
    try transfer(true, 2, 1, 0.25, 0, fixture.pools());
    inline for (.{ fixture.base, fixture.extended, fixture.banded }) |amounts| {
        var pool: usize = 0;
        while (pool < amounts.len / 3) : (pool += 1) {
            try std.testing.expectEqual(@as(f64, 1.5), amounts[pool * 3 + 2]);
            try std.testing.expectEqual(@as(f64, 2.5), amounts[pool * 3 + 1]);
        }
    }
}

test "REDIST macropore salt ISALTG and FHOL gates are exact" {
    var fixture = Fixture{};
    try transfer(false, 2, 1, 0.5, 0, fixture.pools());
    try std.testing.expectEqual(@as(f64, 1), fixture.base[index(3, BaseSpecies.hydrogen, 2)]);
    try std.testing.expectEqual(@as(f64, 2), fixture.extended[index(3, ExtendedSpecies.sulfate, 2)]);
    try std.testing.expectEqual(@as(f64, 2), fixture.banded[index(3, BandedPhosphateSpecies.magnesium_phosphate_1, 2)]);
}

test "REDIST macropore salt validation is atomic" {
    var fixture = Fixture{};
    fixture.banded[index(3, BandedPhosphateSpecies.magnesium_phosphate_1, 2)] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilMacroporeSaltTransferInput, transfer(true, 2, 1, 0.5, 0, fixture.pools()));
    try std.testing.expectEqual(@as(f64, 2), fixture.base[index(3, BaseSpecies.hydrogen, 2)]);
    try std.testing.expectEqual(@as(f64, 2), fixture.base[index(3, BaseSpecies.hydrogen, 1)]);
}
