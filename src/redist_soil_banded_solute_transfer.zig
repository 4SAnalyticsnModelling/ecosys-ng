const std = @import("std");
const definitions = @import("redist_pond_band_solute_transfer.zig");

pub const Pools = definitions.Pools;
pub const BandVolumes = definitions.BandVolumes;

pub const Context = struct {
    source_layer: usize,
    destination_layer: usize,
    source_bulk_density_megagrams_m3: f64,
    destination_bulk_density_megagrams_m3: f64,
    water_fraction: f64, // FWO
    salts_enabled: bool, // ISALTG != 0
    zero_tolerance_m3: f64, // ZERO
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}
fn validateMove(values: []const f64, source: usize, destination: usize, fraction: f64) !void {
    const moved = fraction * values[source];
    if (!std.math.isFinite(moved) or !std.math.isFinite(values[destination] + moved) or
        !std.math.isFinite(values[source] - moved)) return error.NonFiniteSoilBandedSoluteTransferResult;
}
fn move(values: []f64, source: usize, destination: usize, fraction: f64) void {
    const moved = fraction * values[source];
    values[destination] = values[destination] + moved;
    values[source] = values[source] - moved;
}

/// Direct translation of REDIST 9849--9901 under the positive-`BKDS` soil gate.
pub fn transfer(context: Context, volumes: BandVolumes, pools: Pools) !void {
    const layer_count = pools.ammonium_g_n.len;
    if (layer_count == 0 or context.source_layer >= layer_count or context.destination_layer >= layer_count or
        context.source_layer == context.destination_layer)
        return error.SoilBandedSoluteTransferDimensionMismatch;
    inline for (.{ pools.ammonium_g_n, pools.ammonia_g_n, pools.nitrate_g_n, pools.nitrite_g_n, pools.hydrogen_phosphate_g_p, pools.dihydrogen_phosphate_g_p, pools.phosphate_g_p, pools.phosphoric_acid_g_p, pools.iron_phosphate_1_mol, pools.iron_phosphate_2_mol, pools.calcium_phosphate_0_mol, pools.calcium_phosphate_1_mol, pools.calcium_phosphate_2_mol, pools.magnesium_phosphate_1_mol, volumes.ammonium_m3, volumes.nitrate_m3, volumes.phosphate_m3 }) |values| {
        if (values.len != layer_count) return error.SoilBandedSoluteTransferDimensionMismatch;
        if (!finiteSlice(values)) return error.InvalidSoilBandedSoluteTransferInput;
    }
    if (!std.math.isFinite(context.source_bulk_density_megagrams_m3) or !std.math.isFinite(context.destination_bulk_density_megagrams_m3) or
        !std.math.isFinite(context.water_fraction) or context.water_fraction < 0 or context.water_fraction > 1 or
        !std.math.isFinite(context.zero_tolerance_m3) or context.zero_tolerance_m3 < 0)
        return error.InvalidSoilBandedSoluteTransferInput;
    if (context.source_bulk_density_megagrams_m3 <= 0 or context.destination_bulk_density_megagrams_m3 <= 0 or
        context.source_layer == 0) return;

    const move_ammonium = volumes.ammonium_m3[context.destination_layer] > context.zero_tolerance_m3;
    const move_nitrate = volumes.nitrate_m3[context.destination_layer] > context.zero_tolerance_m3;
    const move_phosphate = volumes.phosphate_m3[context.destination_layer] > context.zero_tolerance_m3;
    if (move_ammonium) inline for (.{ pools.ammonium_g_n, pools.ammonia_g_n }) |values|
        try validateMove(values, context.source_layer, context.destination_layer, context.water_fraction);
    if (move_nitrate) inline for (.{ pools.nitrate_g_n, pools.nitrite_g_n }) |values|
        try validateMove(values, context.source_layer, context.destination_layer, context.water_fraction);
    if (move_phosphate) inline for (.{ pools.hydrogen_phosphate_g_p, pools.dihydrogen_phosphate_g_p }) |values|
        try validateMove(values, context.source_layer, context.destination_layer, context.water_fraction);
    if (context.salts_enabled and move_phosphate) inline for (.{ pools.phosphate_g_p, pools.phosphoric_acid_g_p, pools.iron_phosphate_1_mol, pools.iron_phosphate_2_mol, pools.calcium_phosphate_0_mol, pools.calcium_phosphate_1_mol, pools.calcium_phosphate_2_mol, pools.magnesium_phosphate_1_mol }) |values|
        try validateMove(values, context.source_layer, context.destination_layer, context.water_fraction);

    if (move_ammonium) inline for (.{ pools.ammonium_g_n, pools.ammonia_g_n }) |values|
        move(values, context.source_layer, context.destination_layer, context.water_fraction);
    if (move_nitrate) inline for (.{ pools.nitrate_g_n, pools.nitrite_g_n }) |values|
        move(values, context.source_layer, context.destination_layer, context.water_fraction);
    if (move_phosphate) inline for (.{ pools.hydrogen_phosphate_g_p, pools.dihydrogen_phosphate_g_p }) |values|
        move(values, context.source_layer, context.destination_layer, context.water_fraction);
    if (context.salts_enabled and move_phosphate) inline for (.{ pools.phosphate_g_p, pools.phosphoric_acid_g_p, pools.iron_phosphate_1_mol, pools.iron_phosphate_2_mol, pools.calcium_phosphate_0_mol, pools.calcium_phosphate_1_mol, pools.calcium_phosphate_2_mol, pools.magnesium_phosphate_1_mol }) |values|
        move(values, context.source_layer, context.destination_layer, context.water_fraction);
}

const Fixture = struct {
    storage: [14][3]f64 = .{.{ 0, 1, 2 }} ** 14,
    ammonium_volume: [3]f64 = .{ 0, 1, 0 },
    nitrate_volume: [3]f64 = .{ 0, 1, 0 },
    phosphate_volume: [3]f64 = .{ 0, 1, 0 },
    fn pools(self: *Fixture) Pools {
        return .{ .ammonium_g_n = &self.storage[0], .ammonia_g_n = &self.storage[1], .nitrate_g_n = &self.storage[2], .nitrite_g_n = &self.storage[3], .hydrogen_phosphate_g_p = &self.storage[4], .dihydrogen_phosphate_g_p = &self.storage[5], .phosphate_g_p = &self.storage[6], .phosphoric_acid_g_p = &self.storage[7], .iron_phosphate_1_mol = &self.storage[8], .iron_phosphate_2_mol = &self.storage[9], .calcium_phosphate_0_mol = &self.storage[10], .calcium_phosphate_1_mol = &self.storage[11], .calcium_phosphate_2_mol = &self.storage[12], .magnesium_phosphate_1_mol = &self.storage[13] };
    }
    fn volumes(self: *Fixture) BandVolumes {
        return .{ .ammonium_m3 = &self.ammonium_volume, .nitrate_m3 = &self.nitrate_volume, .phosphate_m3 = &self.phosphate_volume };
    }
};
fn testContext() Context {
    return .{ .source_layer = 2, .destination_layer = 1, .source_bulk_density_megagrams_m3 = 1, .destination_bulk_density_megagrams_m3 = 1, .water_fraction = 0.25, .salts_enabled = true, .zero_tolerance_m3 = 1e-12 };
}

test "REDIST soil banded solutes transfer exact independently gated order" {
    var fixture = Fixture{};
    try transfer(testContext(), fixture.volumes(), fixture.pools());
    for (fixture.storage) |pool| {
        try std.testing.expectEqual(@as(f64, 1.5), pool[1]);
        try std.testing.expectEqual(@as(f64, 1.5), pool[2]);
    }
}

test "REDIST soil banded solute gates are independent and conservative" {
    var fixture = Fixture{};
    fixture.nitrate_volume[1] = 0;
    fixture.phosphate_volume[1] = 0;
    try transfer(testContext(), fixture.volumes(), fixture.pools());
    try std.testing.expectEqual(@as(f64, 1.5), fixture.storage[0][1]);
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[2][1]);
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[13][1]);
    try std.testing.expectEqual(@as(f64, 3), fixture.storage[0][1] + fixture.storage[0][2]);
}

test "REDIST soil banded solute validation is atomic" {
    var fixture = Fixture{};
    fixture.storage[13][2] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidSoilBandedSoluteTransferInput, transfer(testContext(), fixture.volumes(), fixture.pools()));
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[0][1]);
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[13][1]);
}
