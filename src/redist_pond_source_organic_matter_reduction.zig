const std = @import("std");
const destination_transfer = @import("redist_pond_organic_matter_transfer.zig");

pub const Pools = destination_transfer.Pools;
pub const organic_fraction_count = destination_transfer.organic_fraction_count;

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 9215--9226: source dissolved organic matter.
pub fn reduce(surface_reappearance_flag: u8, source_layer: usize, remaining_fraction: f64, pools: Pools) !void {
    const expected_len = organic_fraction_count * pools.layer_count;
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or
        pools.micropore_doc_g_c.len != expected_len or pools.micropore_don_g_n.len != expected_len or
        pools.micropore_dop_g_p.len != expected_len or pools.micropore_acetate_g_c.len != expected_len or
        pools.macropore_doc_g_c.len != expected_len or pools.macropore_don_g_n.len != expected_len or
        pools.macropore_dop_g_p.len != expected_len or pools.macropore_acetate_g_c.len != expected_len)
        return error.PondSourceOrganicMatterReductionDimensionMismatch;
    inline for (.{ pools.micropore_doc_g_c, pools.micropore_don_g_n, pools.micropore_dop_g_p, pools.micropore_acetate_g_c, pools.macropore_doc_g_c, pools.macropore_don_g_n, pools.macropore_dop_g_p, pools.macropore_acetate_g_c }) |values|
        if (!finiteSlice(values)) return error.InvalidPondSourceOrganicMatterReductionInput;
    if (!std.math.isFinite(remaining_fraction) or remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceOrganicMatterReductionInput;
    if (surface_reappearance_flag != 0) return;

    for (0..organic_fraction_count) |fraction| {
        const index = fraction * pools.layer_count + source_layer;
        inline for (.{ pools.micropore_doc_g_c, pools.micropore_don_g_n, pools.micropore_dop_g_p, pools.micropore_acetate_g_c, pools.macropore_doc_g_c, pools.macropore_don_g_n, pools.macropore_dop_g_p, pools.macropore_acetate_g_c }) |values|
            if (!std.math.isFinite(remaining_fraction * values[index]))
                return error.NonFinitePondSourceOrganicMatterReductionResult;
    }
    for (0..organic_fraction_count) |fraction| {
        const index = fraction * pools.layer_count + source_layer;
        pools.micropore_doc_g_c[index] = remaining_fraction * pools.micropore_doc_g_c[index];
        pools.micropore_don_g_n[index] = remaining_fraction * pools.micropore_don_g_n[index];
        pools.micropore_dop_g_p[index] = remaining_fraction * pools.micropore_dop_g_p[index];
        pools.micropore_acetate_g_c[index] = remaining_fraction * pools.micropore_acetate_g_c[index];
        pools.macropore_doc_g_c[index] = remaining_fraction * pools.macropore_doc_g_c[index];
        pools.macropore_don_g_n[index] = remaining_fraction * pools.macropore_don_g_n[index];
        pools.macropore_dop_g_p[index] = remaining_fraction * pools.macropore_dop_g_p[index];
        pools.macropore_acetate_g_c[index] = remaining_fraction * pools.macropore_acetate_g_c[index];
    }
}

const Fixture = struct {
    storage: [8][organic_fraction_count * 2]f64 = .{.{ 4, 8 } ** organic_fraction_count} ** 8,
    fn pools(self: *Fixture) Pools {
        return .{ .layer_count = 2, .micropore_doc_g_c = &self.storage[0], .micropore_don_g_n = &self.storage[1], .micropore_dop_g_p = &self.storage[2], .micropore_acetate_g_c = &self.storage[3], .macropore_doc_g_c = &self.storage[4], .macropore_don_g_n = &self.storage[5], .macropore_dop_g_p = &self.storage[6], .macropore_acetate_g_c = &self.storage[7] };
    }
};

test "REDIST source dissolved organics scale five fractions and eight pools" {
    var fixture = Fixture{};
    try reduce(0, 1, 0.25, fixture.pools());
    for (fixture.storage) |pool| for (0..organic_fraction_count) |fraction|
        try std.testing.expectEqual(@as(f64, 2), pool[fraction * 2 + 1]);
}

test "REDIST source dissolved organic guard suppresses reduction" {
    var fixture = Fixture{};
    try reduce(1, 1, 0, fixture.pools());
    for (fixture.storage) |pool| try std.testing.expectEqual(@as(f64, 8), pool[1]);
}

test "REDIST source dissolved organic runtime dimensions are exact" {
    var fixture = Fixture{};
    var pools = fixture.pools();
    pools.macropore_dop_g_p = fixture.storage[6][0 .. organic_fraction_count * 2 - 1];
    try std.testing.expectError(error.PondSourceOrganicMatterReductionDimensionMismatch, reduce(0, 1, 0.5, pools));
}

test "REDIST source dissolved organic validation is atomic" {
    var fixture = Fixture{};
    fixture.storage[7][8] = std.math.inf(f64);
    try std.testing.expectError(error.InvalidPondSourceOrganicMatterReductionInput, reduce(0, 1, 0.5, fixture.pools()));
    try std.testing.expectEqual(@as(f64, 8), fixture.storage[0][1]);
}
