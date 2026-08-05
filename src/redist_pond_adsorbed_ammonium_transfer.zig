const std = @import("std");

pub const Pools = struct {
    nonband_ammonium_mol_n: []f64, // XN4
    banded_ammonium_mol_n: []f64, // XNB
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8986--8989 under `FX == 1.0`.
pub fn transfer(source_layer: usize, destination_layer: usize, fraction: f64, pools: Pools) !void {
    const len = pools.nonband_ammonium_mol_n.len;
    if (len == 0 or source_layer >= len or destination_layer >= len or source_layer == destination_layer or
        pools.banded_ammonium_mol_n.len != len)
        return error.PondAdsorbedAmmoniumTransferDimensionMismatch;
    if (!finiteSlice(pools.nonband_ammonium_mol_n) or !finiteSlice(pools.banded_ammonium_mol_n) or
        !std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
        return error.InvalidPondAdsorbedAmmoniumTransferInput;
    if (fraction != 1.0) return;

    const nonband_ammonium_mol_n = pools.nonband_ammonium_mol_n[destination_layer] +
        fraction * pools.nonband_ammonium_mol_n[source_layer];
    const banded_ammonium_mol_n = pools.banded_ammonium_mol_n[destination_layer] +
        fraction * pools.banded_ammonium_mol_n[source_layer];
    if (!std.math.isFinite(nonband_ammonium_mol_n) or !std.math.isFinite(banded_ammonium_mol_n))
        return error.NonFinitePondAdsorbedAmmoniumTransferResult;

    pools.nonband_ammonium_mol_n[destination_layer] = nonband_ammonium_mol_n;
    pools.banded_ammonium_mol_n[destination_layer] = banded_ammonium_mol_n;
}

test "REDIST full-fraction adsorbed ammonium transfers non-band then band" {
    var nonband = [_]f64{ 2, 1 };
    var banded = [_]f64{ 4, 3 };
    try transfer(0, 1, 1, .{ .nonband_ammonium_mol_n = &nonband, .banded_ammonium_mol_n = &banded });
    try std.testing.expectEqual(@as(f64, 3), nonband[1]);
    try std.testing.expectEqual(@as(f64, 7), banded[1]);
    try std.testing.expectEqual(@as(f64, 2), nonband[0]);
    try std.testing.expectEqual(@as(f64, 4), banded[0]);
}

test "REDIST adsorbed ammonium requires exact full fraction" {
    var nonband = [_]f64{ 2, 1 };
    var banded = [_]f64{ 4, 3 };
    try transfer(0, 1, 0.999999999, .{ .nonband_ammonium_mol_n = &nonband, .banded_ammonium_mol_n = &banded });
    try std.testing.expectEqual(@as(f64, 1), nonband[1]);
    try std.testing.expectEqual(@as(f64, 3), banded[1]);
}

test "REDIST adsorbed ammonium runtime dimensions must agree" {
    var nonband = [_]f64{ 2, 1 };
    var banded = [_]f64{4};
    try std.testing.expectError(
        error.PondAdsorbedAmmoniumTransferDimensionMismatch,
        transfer(0, 1, 1, .{ .nonband_ammonium_mol_n = &nonband, .banded_ammonium_mol_n = &banded }),
    );
}

test "REDIST adsorbed ammonium validation is atomic" {
    var nonband = [_]f64{ 2, 1 };
    var banded = [_]f64{ std.math.inf(f64), 3 };
    try std.testing.expectError(
        error.InvalidPondAdsorbedAmmoniumTransferInput,
        transfer(0, 1, 1, .{ .nonband_ammonium_mol_n = &nonband, .banded_ammonium_mol_n = &banded }),
    );
    try std.testing.expectEqual(@as(f64, 1), nonband[1]);
    try std.testing.expectEqual(@as(f64, 3), banded[1]);
}
