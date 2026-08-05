const std = @import("std");

pub const Pools = struct {
    ammonium_g_n: []f64, // ZNH4B
    ammonia_g_n: []f64, // ZNH3B
    nitrate_g_n: []f64, // ZNO3B
    nitrite_g_n: []f64, // ZNO2B
    hydrogen_phosphate_g_p: []f64, // H1POB
    dihydrogen_phosphate_g_p: []f64, // H2POB
    phosphate_g_p: []f64, // H0POB
    phosphoric_acid_g_p: []f64, // H3POB
    iron_phosphate_1_mol: []f64, // ZFE1PB
    iron_phosphate_2_mol: []f64, // ZFE2PB
    calcium_phosphate_0_mol: []f64, // ZCA0PB
    calcium_phosphate_1_mol: []f64, // ZCA1PB
    calcium_phosphate_2_mol: []f64, // ZCA2PB
    magnesium_phosphate_1_mol: []f64, // ZMG1PB
};

pub const BandVolumes = struct {
    ammonium_m3: []const f64, // VLNHB
    nitrate_m3: []const f64, // VLNOB
    phosphate_m3: []const f64, // VLPOB
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

fn moved(destination: f64, fraction: f64, source: f64) !f64 {
    const result = destination + fraction * source;
    if (!std.math.isFinite(result)) return error.NonFinitePondBandSoluteTransferResult;
    return result;
}

/// Direct translation of REDIST 8696--8734: pond band solutes.
pub fn transfer(
    source_layer: usize,
    destination_layer: usize,
    fraction: f64,
    salts_enabled: bool,
    zero_tolerance: f64,
    volumes: BandVolumes,
    pools: Pools,
) !void {
    const len = pools.ammonium_g_n.len;
    if (len == 0 or source_layer >= len or destination_layer >= len or source_layer == destination_layer or
        pools.ammonia_g_n.len != len or pools.nitrate_g_n.len != len or pools.nitrite_g_n.len != len or
        pools.hydrogen_phosphate_g_p.len != len or pools.dihydrogen_phosphate_g_p.len != len or
        pools.phosphate_g_p.len != len or pools.phosphoric_acid_g_p.len != len or
        pools.iron_phosphate_1_mol.len != len or pools.iron_phosphate_2_mol.len != len or
        pools.calcium_phosphate_0_mol.len != len or pools.calcium_phosphate_1_mol.len != len or
        pools.calcium_phosphate_2_mol.len != len or pools.magnesium_phosphate_1_mol.len != len or
        volumes.ammonium_m3.len != len or volumes.nitrate_m3.len != len or volumes.phosphate_m3.len != len)
        return error.PondBandSoluteTransferDimensionMismatch;
    inline for (.{ pools.ammonium_g_n, pools.ammonia_g_n, pools.nitrate_g_n, pools.nitrite_g_n, pools.hydrogen_phosphate_g_p, pools.dihydrogen_phosphate_g_p, pools.phosphate_g_p, pools.phosphoric_acid_g_p, pools.iron_phosphate_1_mol, pools.iron_phosphate_2_mol, pools.calcium_phosphate_0_mol, pools.calcium_phosphate_1_mol, pools.calcium_phosphate_2_mol, pools.magnesium_phosphate_1_mol, volumes.ammonium_m3, volumes.nitrate_m3, volumes.phosphate_m3 }) |values|
        if (!finiteSlice(values)) return error.InvalidPondBandSoluteTransferInput;
    if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1 or
        !std.math.isFinite(zero_tolerance) or zero_tolerance < 0)
        return error.InvalidPondBandSoluteTransferInput;
    if (source_layer == 0) return;

    const move_ammonium = volumes.ammonium_m3[destination_layer] > zero_tolerance;
    const move_nitrate = volumes.nitrate_m3[destination_layer] > zero_tolerance;
    const move_phosphate = volumes.phosphate_m3[destination_layer] > zero_tolerance;

    const ammonium = if (move_ammonium) try moved(pools.ammonium_g_n[destination_layer], fraction, pools.ammonium_g_n[source_layer]) else 0;
    const ammonia = if (move_ammonium) try moved(pools.ammonia_g_n[destination_layer], fraction, pools.ammonia_g_n[source_layer]) else 0;
    const nitrate = if (move_nitrate) try moved(pools.nitrate_g_n[destination_layer], fraction, pools.nitrate_g_n[source_layer]) else 0;
    const nitrite = if (move_nitrate) try moved(pools.nitrite_g_n[destination_layer], fraction, pools.nitrite_g_n[source_layer]) else 0;
    const hydrogen_phosphate = if (move_phosphate) try moved(pools.hydrogen_phosphate_g_p[destination_layer], fraction, pools.hydrogen_phosphate_g_p[source_layer]) else 0;
    const dihydrogen_phosphate = if (move_phosphate) try moved(pools.dihydrogen_phosphate_g_p[destination_layer], fraction, pools.dihydrogen_phosphate_g_p[source_layer]) else 0;
    const move_extended = salts_enabled and move_phosphate;
    const phosphate = if (move_extended) try moved(pools.phosphate_g_p[destination_layer], fraction, pools.phosphate_g_p[source_layer]) else 0;
    const phosphoric_acid = if (move_extended) try moved(pools.phosphoric_acid_g_p[destination_layer], fraction, pools.phosphoric_acid_g_p[source_layer]) else 0;
    const iron_phosphate_1 = if (move_extended) try moved(pools.iron_phosphate_1_mol[destination_layer], fraction, pools.iron_phosphate_1_mol[source_layer]) else 0;
    const iron_phosphate_2 = if (move_extended) try moved(pools.iron_phosphate_2_mol[destination_layer], fraction, pools.iron_phosphate_2_mol[source_layer]) else 0;
    const calcium_phosphate_0 = if (move_extended) try moved(pools.calcium_phosphate_0_mol[destination_layer], fraction, pools.calcium_phosphate_0_mol[source_layer]) else 0;
    const calcium_phosphate_1 = if (move_extended) try moved(pools.calcium_phosphate_1_mol[destination_layer], fraction, pools.calcium_phosphate_1_mol[source_layer]) else 0;
    const calcium_phosphate_2 = if (move_extended) try moved(pools.calcium_phosphate_2_mol[destination_layer], fraction, pools.calcium_phosphate_2_mol[source_layer]) else 0;
    const magnesium_phosphate_1 = if (move_extended) try moved(pools.magnesium_phosphate_1_mol[destination_layer], fraction, pools.magnesium_phosphate_1_mol[source_layer]) else 0;

    if (move_ammonium) {
        pools.ammonium_g_n[destination_layer] = ammonium;
        pools.ammonia_g_n[destination_layer] = ammonia;
    }
    if (move_nitrate) {
        pools.nitrate_g_n[destination_layer] = nitrate;
        pools.nitrite_g_n[destination_layer] = nitrite;
    }
    if (move_phosphate) {
        pools.hydrogen_phosphate_g_p[destination_layer] = hydrogen_phosphate;
        pools.dihydrogen_phosphate_g_p[destination_layer] = dihydrogen_phosphate;
    }
    if (move_extended) {
        pools.phosphate_g_p[destination_layer] = phosphate;
        pools.phosphoric_acid_g_p[destination_layer] = phosphoric_acid;
        pools.iron_phosphate_1_mol[destination_layer] = iron_phosphate_1;
        pools.iron_phosphate_2_mol[destination_layer] = iron_phosphate_2;
        pools.calcium_phosphate_0_mol[destination_layer] = calcium_phosphate_0;
        pools.calcium_phosphate_1_mol[destination_layer] = calcium_phosphate_1;
        pools.calcium_phosphate_2_mol[destination_layer] = calcium_phosphate_2;
        pools.magnesium_phosphate_1_mol[destination_layer] = magnesium_phosphate_1;
    }
}

const Fixture = struct {
    storage: [14][3]f64 = .{.{ 0, 1, 2 }} ** 14,
    ammonium_volume: [3]f64 = .{ 0, 1, 0 },
    nitrate_volume: [3]f64 = .{ 0, 1, 0 },
    phosphate_volume: [3]f64 = .{ 0, 1, 0 },

    fn pools(self: *Fixture) Pools {
        return .{
            .ammonium_g_n = &self.storage[0],
            .ammonia_g_n = &self.storage[1],
            .nitrate_g_n = &self.storage[2],
            .nitrite_g_n = &self.storage[3],
            .hydrogen_phosphate_g_p = &self.storage[4],
            .dihydrogen_phosphate_g_p = &self.storage[5],
            .phosphate_g_p = &self.storage[6],
            .phosphoric_acid_g_p = &self.storage[7],
            .iron_phosphate_1_mol = &self.storage[8],
            .iron_phosphate_2_mol = &self.storage[9],
            .calcium_phosphate_0_mol = &self.storage[10],
            .calcium_phosphate_1_mol = &self.storage[11],
            .calcium_phosphate_2_mol = &self.storage[12],
            .magnesium_phosphate_1_mol = &self.storage[13],
        };
    }

    fn volumes(self: *Fixture) BandVolumes {
        return .{ .ammonium_m3 = &self.ammonium_volume, .nitrate_m3 = &self.nitrate_volume, .phosphate_m3 = &self.phosphate_volume };
    }
};

test "REDIST pond band solutes transfer all independently gated pools" {
    var fixture = Fixture{};
    try transfer(2, 1, 0.25, true, 1.0e-12, fixture.volumes(), fixture.pools());
    for (fixture.storage) |pool| try std.testing.expectEqual(@as(f64, 1.5), pool[1]);
}

test "REDIST pond band volume gates are independent" {
    var fixture = Fixture{};
    fixture.nitrate_volume[1] = 0;
    fixture.phosphate_volume[1] = 0;
    try transfer(2, 1, 0.5, true, 1.0e-12, fixture.volumes(), fixture.pools());
    try std.testing.expectEqual(@as(f64, 2), fixture.storage[0][1]);
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[2][1]);
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[4][1]);
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[13][1]);
}

test "REDIST disabled salts retain extended phosphate pools" {
    var fixture = Fixture{};
    try transfer(2, 1, 0.5, false, 1.0e-12, fixture.volumes(), fixture.pools());
    try std.testing.expectEqual(@as(f64, 2), fixture.storage[4][1]);
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[6][1]);
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[13][1]);
}

test "REDIST source layer zero skips all band transfer" {
    var fixture = Fixture{};
    try transfer(0, 1, 0.5, true, 1.0e-12, fixture.volumes(), fixture.pools());
    for (fixture.storage) |pool| try std.testing.expectEqual(@as(f64, 1), pool[1]);
}

test "REDIST pond band validation is atomic" {
    var fixture = Fixture{};
    fixture.storage[13][2] = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidPondBandSoluteTransferInput,
        transfer(2, 1, 0.5, true, 1.0e-12, fixture.volumes(), fixture.pools()),
    );
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[0][1]);
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[13][1]);
}
