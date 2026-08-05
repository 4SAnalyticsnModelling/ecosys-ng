const std = @import("std");

/// The legacy science defines K=0..4 as five organic matter fractions.
pub const organic_fraction_count: usize = 5;

pub const Pools = struct {
    layer_count: usize,
    micropore_doc_g_c: []f64, // OQC
    micropore_don_g_n: []f64, // OQN
    micropore_dop_g_p: []f64, // OQP
    micropore_acetate_g_c: []f64, // OQA
    macropore_doc_g_c: []f64, // OQCH
    macropore_don_g_n: []f64, // OQNH
    macropore_dop_g_p: []f64, // OQPH
    macropore_acetate_g_c: []f64, // OQAH

    fn index(self: Pools, fraction: usize, layer: usize) usize {
        return fraction * self.layer_count + layer;
    }
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

fn validateMoved(destination: f64, fraction: f64, source: f64) !void {
    if (!std.math.isFinite(destination + fraction * source))
        return error.NonFinitePondOrganicMatterTransferResult;
}

/// Direct translation of REDIST 8771--8790, including `IFLGL(L,3)==0`.
pub fn transfer(
    surface_reappearance_flag: u8,
    source_layer: usize,
    destination_layer: usize,
    redistribution_fraction: f64,
    pools: Pools,
) !void {
    const expected_len = organic_fraction_count * pools.layer_count;
    if (pools.layer_count == 0 or source_layer >= pools.layer_count or destination_layer >= pools.layer_count or
        source_layer == destination_layer or pools.micropore_doc_g_c.len != expected_len or
        pools.micropore_don_g_n.len != expected_len or pools.micropore_dop_g_p.len != expected_len or
        pools.micropore_acetate_g_c.len != expected_len or pools.macropore_doc_g_c.len != expected_len or
        pools.macropore_don_g_n.len != expected_len or pools.macropore_dop_g_p.len != expected_len or
        pools.macropore_acetate_g_c.len != expected_len)
        return error.PondOrganicMatterTransferDimensionMismatch;
    inline for (.{ pools.micropore_doc_g_c, pools.micropore_don_g_n, pools.micropore_dop_g_p, pools.micropore_acetate_g_c, pools.macropore_doc_g_c, pools.macropore_don_g_n, pools.macropore_dop_g_p, pools.macropore_acetate_g_c }) |values|
        if (!finiteSlice(values)) return error.InvalidPondOrganicMatterTransferInput;
    if (!std.math.isFinite(redistribution_fraction) or redistribution_fraction < 0 or redistribution_fraction > 1)
        return error.InvalidPondOrganicMatterTransferInput;
    if (surface_reappearance_flag != 0) return;

    for (0..organic_fraction_count) |organic_fraction| {
        const source = pools.index(organic_fraction, source_layer);
        const destination = pools.index(organic_fraction, destination_layer);
        try validateMoved(pools.micropore_doc_g_c[destination], redistribution_fraction, pools.micropore_doc_g_c[source]);
        try validateMoved(pools.micropore_don_g_n[destination], redistribution_fraction, pools.micropore_don_g_n[source]);
        try validateMoved(pools.micropore_dop_g_p[destination], redistribution_fraction, pools.micropore_dop_g_p[source]);
        try validateMoved(pools.micropore_acetate_g_c[destination], redistribution_fraction, pools.micropore_acetate_g_c[source]);
        try validateMoved(pools.macropore_doc_g_c[destination], redistribution_fraction, pools.macropore_doc_g_c[source]);
        try validateMoved(pools.macropore_don_g_n[destination], redistribution_fraction, pools.macropore_don_g_n[source]);
        try validateMoved(pools.macropore_dop_g_p[destination], redistribution_fraction, pools.macropore_dop_g_p[source]);
        try validateMoved(pools.macropore_acetate_g_c[destination], redistribution_fraction, pools.macropore_acetate_g_c[source]);
    }
    for (0..organic_fraction_count) |organic_fraction| {
        const source = pools.index(organic_fraction, source_layer);
        const destination = pools.index(organic_fraction, destination_layer);
        pools.micropore_doc_g_c[destination] += redistribution_fraction * pools.micropore_doc_g_c[source];
        pools.micropore_don_g_n[destination] += redistribution_fraction * pools.micropore_don_g_n[source];
        pools.micropore_dop_g_p[destination] += redistribution_fraction * pools.micropore_dop_g_p[source];
        pools.micropore_acetate_g_c[destination] += redistribution_fraction * pools.micropore_acetate_g_c[source];
        pools.macropore_doc_g_c[destination] += redistribution_fraction * pools.macropore_doc_g_c[source];
        pools.macropore_don_g_n[destination] += redistribution_fraction * pools.macropore_don_g_n[source];
        pools.macropore_dop_g_p[destination] += redistribution_fraction * pools.macropore_dop_g_p[source];
        pools.macropore_acetate_g_c[destination] += redistribution_fraction * pools.macropore_acetate_g_c[source];
    }
}

const Fixture = struct {
    storage: [8][organic_fraction_count * 3]f64 = .{.{0} ** (organic_fraction_count * 3)} ** 8,

    fn initialize(self: *Fixture) void {
        for (0..organic_fraction_count) |fraction| {
            for (&self.storage) |*pool| {
                pool[fraction * 3] = @floatFromInt(2 * (fraction + 1));
                pool[fraction * 3 + 1] = @floatFromInt(fraction + 1);
            }
        }
    }

    fn pools(self: *Fixture) Pools {
        return .{
            .layer_count = 3,
            .micropore_doc_g_c = &self.storage[0],
            .micropore_don_g_n = &self.storage[1],
            .micropore_dop_g_p = &self.storage[2],
            .micropore_acetate_g_c = &self.storage[3],
            .macropore_doc_g_c = &self.storage[4],
            .macropore_don_g_n = &self.storage[5],
            .macropore_dop_g_p = &self.storage[6],
            .macropore_acetate_g_c = &self.storage[7],
        };
    }
};

test "REDIST pond organic matter transfers five fractions in exact pool order" {
    var fixture = Fixture{};
    fixture.initialize();
    try transfer(0, 0, 1, 0.25, fixture.pools());
    for (0..organic_fraction_count) |fraction| {
        const expected: f64 = 1.5 * @as(f64, @floatFromInt(fraction + 1));
        for (fixture.storage) |pool| try std.testing.expectEqual(expected, pool[fraction * 3 + 1]);
    }
}

test "REDIST pond organic matter is suppressed by surface reappearance" {
    var fixture = Fixture{};
    fixture.initialize();
    try transfer(1, 0, 1, 0.5, fixture.pools());
    for (fixture.storage) |pool| try std.testing.expectEqual(@as(f64, 1), pool[1]);
}

test "REDIST pond organic matter runtime dimensions are exact" {
    var fixture = Fixture{};
    fixture.initialize();
    var pools = fixture.pools();
    pools.macropore_dop_g_p = fixture.storage[6][0 .. organic_fraction_count * 3 - 1];
    try std.testing.expectError(
        error.PondOrganicMatterTransferDimensionMismatch,
        transfer(0, 0, 1, 0.5, pools),
    );
}

test "REDIST pond organic matter validation is atomic" {
    var fixture = Fixture{};
    fixture.initialize();
    fixture.storage[7][12] = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidPondOrganicMatterTransferInput,
        transfer(0, 0, 1, 0.5, fixture.pools()),
    );
    try std.testing.expectEqual(@as(f64, 1), fixture.storage[0][1]);
    try std.testing.expectEqual(@as(f64, 5), fixture.storage[7][13]);
}
