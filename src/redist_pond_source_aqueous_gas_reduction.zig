const std = @import("std");
const destination_transfer = @import("redist_pond_aqueous_gas_transfer.zig");

pub const Pools = destination_transfer.Pools;

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 9206--9211: source aqueous gases scaled by FY.
pub fn reduce(source_layer: usize, remaining_fraction: f64, pools: Pools) !void {
    const len = pools.carbon_dioxide_g_c.len;
    if (len == 0 or source_layer >= len or pools.methane_g_c.len != len or pools.oxygen_g_o.len != len or
        pools.dinitrogen_g_n.len != len or pools.nitrous_oxide_g_n.len != len or pools.hydrogen_g_h.len != len)
        return error.PondSourceAqueousGasReductionDimensionMismatch;
    inline for (.{ pools.carbon_dioxide_g_c, pools.methane_g_c, pools.oxygen_g_o, pools.dinitrogen_g_n, pools.nitrous_oxide_g_n, pools.hydrogen_g_h }) |values|
        if (!finiteSlice(values)) return error.InvalidPondSourceAqueousGasReductionInput;
    if (!std.math.isFinite(remaining_fraction) or remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceAqueousGasReductionInput;

    const carbon_dioxide = remaining_fraction * pools.carbon_dioxide_g_c[source_layer];
    const methane = remaining_fraction * pools.methane_g_c[source_layer];
    const oxygen = remaining_fraction * pools.oxygen_g_o[source_layer];
    const dinitrogen = remaining_fraction * pools.dinitrogen_g_n[source_layer];
    const nitrous_oxide = remaining_fraction * pools.nitrous_oxide_g_n[source_layer];
    const hydrogen = remaining_fraction * pools.hydrogen_g_h[source_layer];
    inline for (.{ carbon_dioxide, methane, oxygen, dinitrogen, nitrous_oxide, hydrogen }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePondSourceAqueousGasReductionResult;

    pools.carbon_dioxide_g_c[source_layer] = carbon_dioxide;
    pools.methane_g_c[source_layer] = methane;
    pools.oxygen_g_o[source_layer] = oxygen;
    pools.dinitrogen_g_n[source_layer] = dinitrogen;
    pools.nitrous_oxide_g_n[source_layer] = nitrous_oxide;
    pools.hydrogen_g_h[source_layer] = hydrogen;
}

const Fixture = struct {
    storage: [6][2]f64 = .{.{ 4, 8 }} ** 6,
    fn pools(self: *Fixture) Pools {
        return .{ .carbon_dioxide_g_c = &self.storage[0], .methane_g_c = &self.storage[1], .oxygen_g_o = &self.storage[2], .dinitrogen_g_n = &self.storage[3], .nitrous_oxide_g_n = &self.storage[4], .hydrogen_g_h = &self.storage[5] };
    }
};

test "REDIST source aqueous gases scale in exact elemental order" {
    var fixture = Fixture{};
    try reduce(1, 0.25, fixture.pools());
    for (fixture.storage) |pool| try std.testing.expectEqual(@as(f64, 2), pool[1]);
}

test "REDIST source aqueous gases permit layer zero and exhaustion" {
    var fixture = Fixture{};
    try reduce(0, 0, fixture.pools());
    for (fixture.storage) |pool| try std.testing.expectEqual(@as(f64, 0), pool[0]);
}

test "REDIST source aqueous gas runtime dimensions must agree" {
    var fixture = Fixture{};
    var pools = fixture.pools();
    pools.dinitrogen_g_n = fixture.storage[3][0..1];
    try std.testing.expectError(error.PondSourceAqueousGasReductionDimensionMismatch, reduce(1, 0.5, pools));
}

test "REDIST source aqueous gas validation is atomic" {
    var fixture = Fixture{};
    fixture.storage[5][0] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceAqueousGasReductionInput, reduce(0, 0.5, fixture.pools()));
    try std.testing.expectEqual(@as(f64, 4), fixture.storage[0][0]);
}
