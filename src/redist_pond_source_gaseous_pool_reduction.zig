const std = @import("std");
const destination_transfer = @import("redist_pond_gaseous_pool_transfer.zig");

pub const Pools = destination_transfer.Pools;

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 9195--9202, including the open `L0 != 0` guard.
pub fn reduce(source_layer: usize, remaining_fraction: f64, pools: Pools) !void {
    const len = pools.carbon_dioxide_g_c.len;
    if (len == 0 or source_layer >= len or pools.methane_g_c.len != len or pools.oxygen_g_o.len != len or
        pools.dinitrogen_g_n.len != len or pools.nitrous_oxide_g_n.len != len or
        pools.ammonia_g_n.len != len or pools.hydrogen_g_h.len != len)
        return error.PondSourceGaseousPoolReductionDimensionMismatch;
    inline for (.{ pools.carbon_dioxide_g_c, pools.methane_g_c, pools.oxygen_g_o, pools.dinitrogen_g_n, pools.nitrous_oxide_g_n, pools.ammonia_g_n, pools.hydrogen_g_h }) |values|
        if (!finiteSlice(values)) return error.InvalidPondSourceGaseousPoolReductionInput;
    if (!std.math.isFinite(remaining_fraction) or remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceGaseousPoolReductionInput;
    if (source_layer == 0) return;

    const carbon_dioxide = remaining_fraction * pools.carbon_dioxide_g_c[source_layer];
    const methane = remaining_fraction * pools.methane_g_c[source_layer];
    const oxygen = remaining_fraction * pools.oxygen_g_o[source_layer];
    const dinitrogen = remaining_fraction * pools.dinitrogen_g_n[source_layer];
    const nitrous_oxide = remaining_fraction * pools.nitrous_oxide_g_n[source_layer];
    const ammonia = remaining_fraction * pools.ammonia_g_n[source_layer];
    const hydrogen = remaining_fraction * pools.hydrogen_g_h[source_layer];
    inline for (.{ carbon_dioxide, methane, oxygen, dinitrogen, nitrous_oxide, ammonia, hydrogen }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePondSourceGaseousPoolReductionResult;

    pools.carbon_dioxide_g_c[source_layer] = carbon_dioxide;
    pools.methane_g_c[source_layer] = methane;
    pools.oxygen_g_o[source_layer] = oxygen;
    pools.dinitrogen_g_n[source_layer] = dinitrogen;
    pools.nitrous_oxide_g_n[source_layer] = nitrous_oxide;
    pools.ammonia_g_n[source_layer] = ammonia;
    pools.hydrogen_g_h[source_layer] = hydrogen;
}

const Fixture = struct {
    storage: [7][2]f64 = .{.{ 4, 8 }} ** 7,
    fn pools(self: *Fixture) Pools {
        return .{ .carbon_dioxide_g_c = &self.storage[0], .methane_g_c = &self.storage[1], .oxygen_g_o = &self.storage[2], .dinitrogen_g_n = &self.storage[3], .nitrous_oxide_g_n = &self.storage[4], .ammonia_g_n = &self.storage[5], .hydrogen_g_h = &self.storage[6] };
    }
};

test "REDIST source gaseous pools scale in exact elemental order" {
    var fixture = Fixture{};
    try reduce(1, 0.25, fixture.pools());
    for (fixture.storage) |pool| try std.testing.expectEqual(@as(f64, 2), pool[1]);
}

test "REDIST source gaseous pools exclude layer zero" {
    var fixture = Fixture{};
    try reduce(0, 0, fixture.pools());
    for (fixture.storage) |pool| try std.testing.expectEqual(@as(f64, 4), pool[0]);
}

test "REDIST source gaseous pool runtime dimensions must agree" {
    var fixture = Fixture{};
    var pools = fixture.pools();
    pools.oxygen_g_o = fixture.storage[2][0..1];
    try std.testing.expectError(error.PondSourceGaseousPoolReductionDimensionMismatch, reduce(1, 0.5, pools));
}

test "REDIST source gaseous pool validation is atomic" {
    var fixture = Fixture{};
    fixture.storage[6][1] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceGaseousPoolReductionInput, reduce(1, 0.5, fixture.pools()));
    try std.testing.expectEqual(@as(f64, 8), fixture.storage[0][1]);
}
