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

/// Direct translation of REDIST 9111--9116: scale source pond N/P solutes by FY.
pub fn reduce(source_layer: usize, remaining_fraction: f64, pools: Pools) !void {
    const len = pools.ammonium_g_n.len;
    if (len == 0 or source_layer >= len or pools.ammonia_g_n.len != len or pools.nitrate_g_n.len != len or
        pools.nitrite_g_n.len != len or pools.hydrogen_phosphate_g_p.len != len or
        pools.dihydrogen_phosphate_g_p.len != len)
        return error.PondSourceNutrientReductionDimensionMismatch;
    inline for (.{ pools.ammonium_g_n, pools.ammonia_g_n, pools.nitrate_g_n, pools.nitrite_g_n, pools.hydrogen_phosphate_g_p, pools.dihydrogen_phosphate_g_p }) |values|
        if (!finiteSlice(values)) return error.InvalidPondSourceNutrientReductionInput;
    if (!std.math.isFinite(remaining_fraction) or remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceNutrientReductionInput;

    const ammonium_g_n = remaining_fraction * pools.ammonium_g_n[source_layer];
    const ammonia_g_n = remaining_fraction * pools.ammonia_g_n[source_layer];
    const nitrate_g_n = remaining_fraction * pools.nitrate_g_n[source_layer];
    const nitrite_g_n = remaining_fraction * pools.nitrite_g_n[source_layer];
    const hydrogen_phosphate_g_p = remaining_fraction * pools.hydrogen_phosphate_g_p[source_layer];
    const dihydrogen_phosphate_g_p = remaining_fraction * pools.dihydrogen_phosphate_g_p[source_layer];
    inline for (.{ ammonium_g_n, ammonia_g_n, nitrate_g_n, nitrite_g_n, hydrogen_phosphate_g_p, dihydrogen_phosphate_g_p }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePondSourceNutrientReductionResult;

    pools.ammonium_g_n[source_layer] = ammonium_g_n;
    pools.ammonia_g_n[source_layer] = ammonia_g_n;
    pools.nitrate_g_n[source_layer] = nitrate_g_n;
    pools.nitrite_g_n[source_layer] = nitrite_g_n;
    pools.hydrogen_phosphate_g_p[source_layer] = hydrogen_phosphate_g_p;
    pools.dihydrogen_phosphate_g_p[source_layer] = dihydrogen_phosphate_g_p;
}

const Fixture = struct {
    ammonium: [2]f64 = .{ 4, 1 },
    ammonia: [2]f64 = .{ 8, 2 },
    nitrate: [2]f64 = .{ 12, 3 },
    nitrite: [2]f64 = .{ 16, 4 },
    hydrogen_phosphate: [2]f64 = .{ 20, 5 },
    dihydrogen_phosphate: [2]f64 = .{ 24, 6 },
    fn pools(self: *Fixture) Pools {
        return .{ .ammonium_g_n = &self.ammonium, .ammonia_g_n = &self.ammonia, .nitrate_g_n = &self.nitrate, .nitrite_g_n = &self.nitrite, .hydrogen_phosphate_g_p = &self.hydrogen_phosphate, .dihydrogen_phosphate_g_p = &self.dihydrogen_phosphate };
    }
};

test "REDIST source pond nutrients scale in exact N then P order" {
    var fixture = Fixture{};
    try reduce(0, 0.25, fixture.pools());
    try std.testing.expectEqual(@as(f64, 1), fixture.ammonium[0]);
    try std.testing.expectEqual(@as(f64, 2), fixture.ammonia[0]);
    try std.testing.expectEqual(@as(f64, 3), fixture.nitrate[0]);
    try std.testing.expectEqual(@as(f64, 4), fixture.nitrite[0]);
    try std.testing.expectEqual(@as(f64, 5), fixture.hydrogen_phosphate[0]);
    try std.testing.expectEqual(@as(f64, 6), fixture.dihydrogen_phosphate[0]);
}

test "REDIST source pond nutrient reduction permits layer zero and exhaustion" {
    var fixture = Fixture{};
    try reduce(0, 0, fixture.pools());
    try std.testing.expectEqual(@as(f64, 0), fixture.ammonium[0]);
    try std.testing.expectEqual(@as(f64, 0), fixture.dihydrogen_phosphate[0]);
}

test "REDIST source pond nutrient runtime dimensions must agree" {
    var fixture = Fixture{};
    var pools = fixture.pools();
    pools.nitrate_g_n = fixture.nitrate[0..1];
    try std.testing.expectError(error.PondSourceNutrientReductionDimensionMismatch, reduce(0, 0.5, pools));
}

test "REDIST source pond nutrient validation is atomic" {
    var fixture = Fixture{};
    fixture.dihydrogen_phosphate[0] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceNutrientReductionInput, reduce(0, 0.5, fixture.pools()));
    try std.testing.expectEqual(@as(f64, 4), fixture.ammonium[0]);
}
