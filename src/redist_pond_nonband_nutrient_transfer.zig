const std = @import("std");

pub const Pools = struct {
    ammonium_g_n: []f64, // ZNH4S
    ammonia_g_n: []f64, // ZNH3S
    nitrate_g_n: []f64, // ZNO3S
    nitrite_g_n: []f64, // ZNO2S
    hydrogen_phosphate_g_p: []f64, // H1PO4
    dihydrogen_phosphate_g_p: []f64, // H2PO4
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8592--8603: pond non-band N and P solutes.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, pools: Pools) !void {
    const len = pools.ammonium_g_n.len;
    if (len == 0 or source_layer >= len or destination_layer >= len or source_layer == destination_layer or
        pools.ammonia_g_n.len != len or pools.nitrate_g_n.len != len or pools.nitrite_g_n.len != len or
        pools.hydrogen_phosphate_g_p.len != len or pools.dihydrogen_phosphate_g_p.len != len)
        return error.PondNonbandNutrientTransferDimensionMismatch;
    inline for (.{ pools.ammonium_g_n, pools.ammonia_g_n, pools.nitrate_g_n, pools.nitrite_g_n, pools.hydrogen_phosphate_g_p, pools.dihydrogen_phosphate_g_p }) |values|
        if (!finiteSlice(values)) return error.InvalidPondNonbandNutrientTransferInput;
    if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondNonbandNutrientTransferInput;

    const ammonium_g_n = pools.ammonium_g_n[destination_layer] + fraction * pools.ammonium_g_n[source_layer];
    const ammonia_g_n = pools.ammonia_g_n[destination_layer] + fraction * pools.ammonia_g_n[source_layer];
    const nitrate_g_n = pools.nitrate_g_n[destination_layer] + fraction * pools.nitrate_g_n[source_layer];
    const nitrite_g_n = pools.nitrite_g_n[destination_layer] + fraction * pools.nitrite_g_n[source_layer];
    const hydrogen_phosphate_g_p = pools.hydrogen_phosphate_g_p[destination_layer] +
        fraction * pools.hydrogen_phosphate_g_p[source_layer];
    const dihydrogen_phosphate_g_p = pools.dihydrogen_phosphate_g_p[destination_layer] +
        fraction * pools.dihydrogen_phosphate_g_p[source_layer];
    inline for (.{ ammonium_g_n, ammonia_g_n, nitrate_g_n, nitrite_g_n, hydrogen_phosphate_g_p, dihydrogen_phosphate_g_p }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePondNonbandNutrientTransferResult;

    pools.ammonium_g_n[destination_layer] = ammonium_g_n;
    pools.ammonia_g_n[destination_layer] = ammonia_g_n;
    pools.nitrate_g_n[destination_layer] = nitrate_g_n;
    pools.nitrite_g_n[destination_layer] = nitrite_g_n;
    pools.hydrogen_phosphate_g_p[destination_layer] = hydrogen_phosphate_g_p;
    pools.dihydrogen_phosphate_g_p[destination_layer] = dihydrogen_phosphate_g_p;
}

const Fixture = struct {
    ammonium: [3]f64 = .{ 4, 1, 0 },
    ammonia: [3]f64 = .{ 2, 3, 0 },
    nitrate: [3]f64 = .{ 8, 2, 0 },
    nitrite: [3]f64 = .{ 6, 1.5, 0 },
    hydrogen_phosphate: [3]f64 = .{ 10, 4, 0 },
    dihydrogen_phosphate: [3]f64 = .{ 12, 5, 0 },

    fn pools(self: *Fixture) Pools {
        return .{
            .ammonium_g_n = &self.ammonium,
            .ammonia_g_n = &self.ammonia,
            .nitrate_g_n = &self.nitrate,
            .nitrite_g_n = &self.nitrite,
            .hydrogen_phosphate_g_p = &self.hydrogen_phosphate,
            .dihydrogen_phosphate_g_p = &self.dihydrogen_phosphate,
        };
    }
};

test "REDIST pond non-band nutrient pools transfer in source order" {
    var fixture = Fixture{};
    try transfer(0, 1, 0.25, fixture.pools());
    try std.testing.expectEqual(@as(f64, 2), fixture.ammonium[1]);
    try std.testing.expectEqual(@as(f64, 3.5), fixture.ammonia[1]);
    try std.testing.expectEqual(@as(f64, 4), fixture.nitrate[1]);
    try std.testing.expectEqual(@as(f64, 3), fixture.nitrite[1]);
    try std.testing.expectEqual(@as(f64, 6.5), fixture.hydrogen_phosphate[1]);
    try std.testing.expectEqual(@as(f64, 8), fixture.dihydrogen_phosphate[1]);
    try std.testing.expectEqual(@as(f64, 4), fixture.ammonium[0]);
}

test "REDIST zero fraction preserves non-band nutrient destination" {
    var fixture = Fixture{};
    try transfer(0, 1, 0, fixture.pools());
    try std.testing.expectEqual(@as(f64, 1), fixture.ammonium[1]);
    try std.testing.expectEqual(@as(f64, 5), fixture.dihydrogen_phosphate[1]);
}

test "REDIST pond non-band nutrient dimensions are runtime validated" {
    var fixture = Fixture{};
    var pools = fixture.pools();
    pools.nitrate_g_n = fixture.nitrate[0..2];
    try std.testing.expectError(
        error.PondNonbandNutrientTransferDimensionMismatch,
        transfer(0, 1, 0.5, pools),
    );
}

test "REDIST pond non-band nutrient transfer validates atomically" {
    var fixture = Fixture{};
    fixture.nitrite[0] = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidPondNonbandNutrientTransferInput,
        transfer(0, 1, 0.5, fixture.pools()),
    );
    try std.testing.expectEqual(@as(f64, 1), fixture.ammonium[1]);
    try std.testing.expectEqual(@as(f64, 5), fixture.dihydrogen_phosphate[1]);
}
