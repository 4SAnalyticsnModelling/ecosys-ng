const std = @import("std");
const destination_transfer = @import("adsorbed_ammonium_transfer.zig");

pub const Pools = destination_transfer.Pools;

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 9342--9343 under the enclosing `FX == 1.0`.
pub fn reduce(source_layer: usize, redistribution_fraction: f64, remaining_fraction: f64, pools: Pools) !void {
    const layer_count = pools.nonband_ammonium_mol_n.len;
    if (layer_count == 0 or source_layer >= layer_count or pools.banded_ammonium_mol_n.len != layer_count)
        return error.PondSourceAdsorbedAmmoniumReductionDimensionMismatch;
    if (!finiteSlice(pools.nonband_ammonium_mol_n) or !finiteSlice(pools.banded_ammonium_mol_n) or
        !std.math.isFinite(redistribution_fraction) or redistribution_fraction < 0 or redistribution_fraction > 1 or
        !std.math.isFinite(remaining_fraction) or remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPondSourceAdsorbedAmmoniumReductionInput;
    if (redistribution_fraction != 1.0) return;

    const nonband_ammonium_mol_n = remaining_fraction * pools.nonband_ammonium_mol_n[source_layer];
    const banded_ammonium_mol_n = remaining_fraction * pools.banded_ammonium_mol_n[source_layer];
    if (!std.math.isFinite(nonband_ammonium_mol_n) or !std.math.isFinite(banded_ammonium_mol_n))
        return error.NonFinitePondSourceAdsorbedAmmoniumReductionResult;

    pools.nonband_ammonium_mol_n[source_layer] = nonband_ammonium_mol_n;
    pools.banded_ammonium_mol_n[source_layer] = banded_ammonium_mol_n;
}

test "REDIST source adsorbed ammonium scales non-band then band" {
    var nonband = [_]f64{ 8, 3 };
    var banded = [_]f64{ 12, 5 };
    try reduce(0, 1, 0.25, .{ .nonband_ammonium_mol_n = &nonband, .banded_ammonium_mol_n = &banded });
    try std.testing.expectEqual(@as(f64, 2), nonband[0]);
    try std.testing.expectEqual(@as(f64, 3), banded[0]);
    try std.testing.expectEqual(@as(f64, 3), nonband[1]);
    try std.testing.expectEqual(@as(f64, 5), banded[1]);
}

test "REDIST source adsorbed ammonium permits layer zero and exhaustion" {
    var nonband = [_]f64{2};
    var banded = [_]f64{4};
    try reduce(0, 1, 0, .{ .nonband_ammonium_mol_n = &nonband, .banded_ammonium_mol_n = &banded });
    try std.testing.expectEqual(@as(f64, 0), nonband[0]);
    try std.testing.expectEqual(@as(f64, 0), banded[0]);
}

test "REDIST source adsorbed ammonium requires exact full fraction" {
    var nonband = [_]f64{2};
    var banded = [_]f64{4};
    try reduce(0, 0.999999999, 0, .{ .nonband_ammonium_mol_n = &nonband, .banded_ammonium_mol_n = &banded });
    try std.testing.expectEqual(@as(f64, 2), nonband[0]);
    try std.testing.expectEqual(@as(f64, 4), banded[0]);
}

test "REDIST source adsorbed ammonium validation is atomic" {
    var nonband = [_]f64{2};
    var banded = [_]f64{std.math.inf(f64)};
    try std.testing.expectError(
        error.InvalidPondSourceAdsorbedAmmoniumReductionInput,
        reduce(0, 1, 0.5, .{ .nonband_ammonium_mol_n = &nonband, .banded_ammonium_mol_n = &banded }),
    );
    try std.testing.expectEqual(@as(f64, 2), nonband[0]);
}
