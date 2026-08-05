const std = @import("std");
const destination_transfer = @import("redist_pond_band_solute_transfer.zig");

pub const Pools = destination_transfer.Pools;
pub const BandVolumes = destination_transfer.BandVolumes;

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

fn scaled(fraction: f64, value: f64) !f64 {
    const result = fraction * value;
    if (!std.math.isFinite(result)) return error.NonFinitePondSourceBandSoluteReductionResult;
    return result;
}

/// Direct translation of REDIST 9167--9191: source band solutes.
pub fn reduce(source_layer: usize, destination_layer: usize, remaining_fraction: f64, salts_enabled: bool, zero_tolerance: f64, volumes: BandVolumes, pools: Pools) !void {
    const len = pools.ammonium_g_n.len;
    if (len == 0 or source_layer >= len or destination_layer >= len or source_layer == destination_layer or
        pools.ammonia_g_n.len != len or pools.nitrate_g_n.len != len or pools.nitrite_g_n.len != len or
        pools.hydrogen_phosphate_g_p.len != len or pools.dihydrogen_phosphate_g_p.len != len or
        pools.phosphate_g_p.len != len or pools.phosphoric_acid_g_p.len != len or
        pools.iron_phosphate_1_mol.len != len or pools.iron_phosphate_2_mol.len != len or
        pools.calcium_phosphate_0_mol.len != len or pools.calcium_phosphate_1_mol.len != len or
        pools.calcium_phosphate_2_mol.len != len or pools.magnesium_phosphate_1_mol.len != len or
        volumes.ammonium_m3.len != len or volumes.nitrate_m3.len != len or volumes.phosphate_m3.len != len)
        return error.PondSourceBandSoluteReductionDimensionMismatch;
    inline for (.{ pools.ammonium_g_n, pools.ammonia_g_n, pools.nitrate_g_n, pools.nitrite_g_n, pools.hydrogen_phosphate_g_p, pools.dihydrogen_phosphate_g_p, pools.phosphate_g_p, pools.phosphoric_acid_g_p, pools.iron_phosphate_1_mol, pools.iron_phosphate_2_mol, pools.calcium_phosphate_0_mol, pools.calcium_phosphate_1_mol, pools.calcium_phosphate_2_mol, pools.magnesium_phosphate_1_mol, volumes.ammonium_m3, volumes.nitrate_m3, volumes.phosphate_m3 }) |values|
        if (!finiteSlice(values)) return error.InvalidPondSourceBandSoluteReductionInput;
    if (!std.math.isFinite(remaining_fraction) or remaining_fraction < 0 or remaining_fraction > 1 or
        !std.math.isFinite(zero_tolerance) or zero_tolerance < 0)
        return error.InvalidPondSourceBandSoluteReductionInput;
    if (source_layer == 0) return;

    const move_ammonium = volumes.ammonium_m3[destination_layer] > zero_tolerance;
    const move_nitrate = volumes.nitrate_m3[destination_layer] > zero_tolerance;
    const move_phosphate = volumes.phosphate_m3[destination_layer] > zero_tolerance;
    const move_extended = salts_enabled and move_phosphate;
    const ammonium = if (move_ammonium) try scaled(remaining_fraction, pools.ammonium_g_n[source_layer]) else 0;
    const ammonia = if (move_ammonium) try scaled(remaining_fraction, pools.ammonia_g_n[source_layer]) else 0;
    const nitrate = if (move_nitrate) try scaled(remaining_fraction, pools.nitrate_g_n[source_layer]) else 0;
    const nitrite = if (move_nitrate) try scaled(remaining_fraction, pools.nitrite_g_n[source_layer]) else 0;
    const hydrogen_phosphate = if (move_phosphate) try scaled(remaining_fraction, pools.hydrogen_phosphate_g_p[source_layer]) else 0;
    const dihydrogen_phosphate = if (move_phosphate) try scaled(remaining_fraction, pools.dihydrogen_phosphate_g_p[source_layer]) else 0;
    const phosphate = if (move_extended) try scaled(remaining_fraction, pools.phosphate_g_p[source_layer]) else 0;
    const phosphoric_acid = if (move_extended) try scaled(remaining_fraction, pools.phosphoric_acid_g_p[source_layer]) else 0;
    const iron_phosphate_1 = if (move_extended) try scaled(remaining_fraction, pools.iron_phosphate_1_mol[source_layer]) else 0;
    const iron_phosphate_2 = if (move_extended) try scaled(remaining_fraction, pools.iron_phosphate_2_mol[source_layer]) else 0;
    const calcium_phosphate_0 = if (move_extended) try scaled(remaining_fraction, pools.calcium_phosphate_0_mol[source_layer]) else 0;
    const calcium_phosphate_1 = if (move_extended) try scaled(remaining_fraction, pools.calcium_phosphate_1_mol[source_layer]) else 0;
    const calcium_phosphate_2 = if (move_extended) try scaled(remaining_fraction, pools.calcium_phosphate_2_mol[source_layer]) else 0;
    const magnesium_phosphate_1 = if (move_extended) try scaled(remaining_fraction, pools.magnesium_phosphate_1_mol[source_layer]) else 0;

    if (move_ammonium) {
        pools.ammonium_g_n[source_layer] = ammonium;
        pools.ammonia_g_n[source_layer] = ammonia;
    }
    if (move_nitrate) {
        pools.nitrate_g_n[source_layer] = nitrate;
        pools.nitrite_g_n[source_layer] = nitrite;
    }
    if (move_phosphate) {
        pools.hydrogen_phosphate_g_p[source_layer] = hydrogen_phosphate;
        pools.dihydrogen_phosphate_g_p[source_layer] = dihydrogen_phosphate;
    }
    if (move_extended) {
        pools.phosphate_g_p[source_layer] = phosphate;
        pools.phosphoric_acid_g_p[source_layer] = phosphoric_acid;
        pools.iron_phosphate_1_mol[source_layer] = iron_phosphate_1;
        pools.iron_phosphate_2_mol[source_layer] = iron_phosphate_2;
        pools.calcium_phosphate_0_mol[source_layer] = calcium_phosphate_0;
        pools.calcium_phosphate_1_mol[source_layer] = calcium_phosphate_1;
        pools.calcium_phosphate_2_mol[source_layer] = calcium_phosphate_2;
        pools.magnesium_phosphate_1_mol[source_layer] = magnesium_phosphate_1;
    }
}

const Fixture = struct {
    storage: [14][3]f64 = .{.{ 0, 4, 8 }} ** 14,
    ammonium_volume: [3]f64 = .{ 0, 1, 1 },
    nitrate_volume: [3]f64 = .{ 0, 1, 1 },
    phosphate_volume: [3]f64 = .{ 0, 1, 1 },
    fn pools(self: *Fixture) Pools {
        return .{ .ammonium_g_n = &self.storage[0], .ammonia_g_n = &self.storage[1], .nitrate_g_n = &self.storage[2], .nitrite_g_n = &self.storage[3], .hydrogen_phosphate_g_p = &self.storage[4], .dihydrogen_phosphate_g_p = &self.storage[5], .phosphate_g_p = &self.storage[6], .phosphoric_acid_g_p = &self.storage[7], .iron_phosphate_1_mol = &self.storage[8], .iron_phosphate_2_mol = &self.storage[9], .calcium_phosphate_0_mol = &self.storage[10], .calcium_phosphate_1_mol = &self.storage[11], .calcium_phosphate_2_mol = &self.storage[12], .magnesium_phosphate_1_mol = &self.storage[13] };
    }
    fn volumes(self: *Fixture) BandVolumes {
        return .{ .ammonium_m3 = &self.ammonium_volume, .nitrate_m3 = &self.nitrate_volume, .phosphate_m3 = &self.phosphate_volume };
    }
};

test "REDIST source band solutes scale every enabled gated pool" {
    var fixture = Fixture{};
    try reduce(2, 1, 0.25, true, 1.0e-12, fixture.volumes(), fixture.pools());
    for (fixture.storage) |pool| try std.testing.expectEqual(@as(f64, 2), pool[2]);
}

test "REDIST source band gates are independent and source zero is skipped" {
    var fixture = Fixture{};
    fixture.nitrate_volume[1] = 0;
    fixture.phosphate_volume[1] = 0;
    try reduce(2, 1, 0.5, true, 1.0e-12, fixture.volumes(), fixture.pools());
    try std.testing.expectEqual(@as(f64, 4), fixture.storage[0][2]);
    try std.testing.expectEqual(@as(f64, 8), fixture.storage[2][2]);
    try reduce(0, 1, 0, true, 1.0e-12, fixture.volumes(), fixture.pools());
    try std.testing.expectEqual(@as(f64, 0), fixture.storage[0][0]);
}

test "REDIST source band validation is atomic" {
    var fixture = Fixture{};
    fixture.storage[13][2] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceBandSoluteReductionInput, reduce(2, 1, 0.5, true, 1.0e-12, fixture.volumes(), fixture.pools()));
    try std.testing.expectEqual(@as(f64, 8), fixture.storage[0][2]);
}
