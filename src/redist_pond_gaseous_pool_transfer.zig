const std = @import("std");

pub const Pools = struct {
    carbon_dioxide_g_c: []f64, // CO2G
    methane_g_c: []f64, // CH4G
    oxygen_g_o: []f64, // OXYG
    dinitrogen_g_n: []f64, // Z2GG
    nitrous_oxide_g_n: []f64, // Z2OG
    ammonia_g_n: []f64, // ZNH3G
    hydrogen_g_h: []f64, // H2GG
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8738--8752, including the enclosing `L0 != 0`.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, pools: Pools) !void {
    const len = pools.carbon_dioxide_g_c.len;
    if (len == 0 or source_layer >= len or destination_layer >= len or source_layer == destination_layer or
        pools.methane_g_c.len != len or pools.oxygen_g_o.len != len or pools.dinitrogen_g_n.len != len or
        pools.nitrous_oxide_g_n.len != len or pools.ammonia_g_n.len != len or pools.hydrogen_g_h.len != len)
        return error.PondGaseousPoolTransferDimensionMismatch;
    inline for (.{ pools.carbon_dioxide_g_c, pools.methane_g_c, pools.oxygen_g_o, pools.dinitrogen_g_n, pools.nitrous_oxide_g_n, pools.ammonia_g_n, pools.hydrogen_g_h }) |values|
        if (!finiteSlice(values)) return error.InvalidPondGaseousPoolTransferInput;
    if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondGaseousPoolTransferInput;
    if (source_layer == 0) return;

    const carbon_dioxide_g_c = pools.carbon_dioxide_g_c[destination_layer] + fraction * pools.carbon_dioxide_g_c[source_layer];
    const methane_g_c = pools.methane_g_c[destination_layer] + fraction * pools.methane_g_c[source_layer];
    const oxygen_g_o = pools.oxygen_g_o[destination_layer] + fraction * pools.oxygen_g_o[source_layer];
    const dinitrogen_g_n = pools.dinitrogen_g_n[destination_layer] + fraction * pools.dinitrogen_g_n[source_layer];
    const nitrous_oxide_g_n = pools.nitrous_oxide_g_n[destination_layer] + fraction * pools.nitrous_oxide_g_n[source_layer];
    const ammonia_g_n = pools.ammonia_g_n[destination_layer] + fraction * pools.ammonia_g_n[source_layer];
    const hydrogen_g_h = pools.hydrogen_g_h[destination_layer] + fraction * pools.hydrogen_g_h[source_layer];
    inline for (.{ carbon_dioxide_g_c, methane_g_c, oxygen_g_o, dinitrogen_g_n, nitrous_oxide_g_n, ammonia_g_n, hydrogen_g_h }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePondGaseousPoolTransferResult;

    pools.carbon_dioxide_g_c[destination_layer] = carbon_dioxide_g_c;
    pools.methane_g_c[destination_layer] = methane_g_c;
    pools.oxygen_g_o[destination_layer] = oxygen_g_o;
    pools.dinitrogen_g_n[destination_layer] = dinitrogen_g_n;
    pools.nitrous_oxide_g_n[destination_layer] = nitrous_oxide_g_n;
    pools.ammonia_g_n[destination_layer] = ammonia_g_n;
    pools.hydrogen_g_h[destination_layer] = hydrogen_g_h;
}

const Fixture = struct {
    storage: [7][3]f64 = .{.{ 0, 1, 2 }} ** 7,

    fn pools(self: *Fixture) Pools {
        return .{
            .carbon_dioxide_g_c = &self.storage[0],
            .methane_g_c = &self.storage[1],
            .oxygen_g_o = &self.storage[2],
            .dinitrogen_g_n = &self.storage[3],
            .nitrous_oxide_g_n = &self.storage[4],
            .ammonia_g_n = &self.storage[5],
            .hydrogen_g_h = &self.storage[6],
        };
    }
};

test "REDIST pond gaseous pools transfer in exact source order" {
    var fixture = Fixture{};
    try transfer(2, 1, 0.25, fixture.pools());
    for (fixture.storage) |pool| {
        try std.testing.expectEqual(@as(f64, 1.5), pool[1]);
        try std.testing.expectEqual(@as(f64, 2), pool[2]);
    }
}

test "REDIST pond gaseous pools remain when source is layer zero" {
    var fixture = Fixture{};
    try transfer(0, 1, 0.5, fixture.pools());
    for (fixture.storage) |pool| try std.testing.expectEqual(@as(f64, 1), pool[1]);
}

test "REDIST pond gaseous pools use runtime dimensions" {
    var fixture = Fixture{};
    var pools = fixture.pools();
    pools.oxygen_g_o = fixture.storage[2][0..2];
    try std.testing.expectError(
        error.PondGaseousPoolTransferDimensionMismatch,
        transfer(2, 1, 0.5, pools),
    );
}

test "REDIST pond gaseous pool validation is atomic" {
    var fixture = Fixture{};
    fixture.storage[6][2] = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidPondGaseousPoolTransferInput,
        transfer(2, 1, 0.5, fixture.pools()),
    );
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[0][1]);
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[6][1]);
}
