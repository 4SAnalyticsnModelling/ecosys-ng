const std = @import("std");

/// Dissolved organic micropore pools for one layer, K=0-4 (g C, N, P).
pub const DissolvedOrganicMicropores = struct {
    /// OQC[K=0-4]. DOC in micropores (g C).
    oqc: [5]f64,
    /// OQN[K=0-4]. DON in micropores (g N).
    oqn: [5]f64,
    /// OQP[K=0-4]. DOP in micropores (g P).
    oqp: [5]f64,
};

/// Root-soil nonstructural C/N/P exchange fluxes from extract.f (g C,N,P step-1).
pub const ExudationFluxes = struct {
    /// TDFOMC[K=0-4].
    tdfomc: [5]f64,
    /// TDFOMN[K=0-4].
    tdfomn: [5]f64,
    /// TDFOMP[K=0-4].
    tdfomp: [5]f64,
};

pub const RuntimePools = struct { carbon_g: []f64, nitrogen_g: []f64, phosphorus_g: []f64 };
pub const RuntimeFluxes = struct { carbon_g: []const f64, nitrogen_g: []const f64, phosphorus_g: []const f64 };

/// Runtime layer-major/class-minor translation of redist.f lines 6115--6119.
pub fn applyLayers(layer_count: usize, class_count: usize, pools: RuntimePools, fluxes: RuntimeFluxes) !void {
    if (layer_count == 0 or class_count == 0) return error.InvalidExudationDimensions;
    const count = std.math.mul(usize, layer_count, class_count) catch return error.ExudationDimensionOverflow;
    inline for (.{ pools.carbon_g.len, pools.nitrogen_g.len, pools.phosphorus_g.len, fluxes.carbon_g.len, fluxes.nitrogen_g.len, fluxes.phosphorus_g.len }) |len|
        if (len != count) return error.ExudationDimensionMismatch;
    for (0..layer_count) |layer| for (0..class_count) |class| {
        const i = layer * class_count + class;
        inline for (.{ pools.carbon_g[i], pools.nitrogen_g[i], pools.phosphorus_g[i] }) |v| if (!std.math.isFinite(v)) return error.InvalidExudationPool;
        inline for (.{ fluxes.carbon_g[i], fluxes.nitrogen_g[i], fluxes.phosphorus_g[i] }) |v| if (!std.math.isFinite(v)) return error.InvalidExudationFlux;
        pools.carbon_g[i] = pools.carbon_g[i] + fluxes.carbon_g[i];
        pools.nitrogen_g[i] = pools.nitrogen_g[i] + fluxes.nitrogen_g[i];
        pools.phosphorus_g[i] = pools.phosphorus_g[i] + fluxes.phosphorus_g[i];
        inline for (.{ pools.carbon_g[i], pools.nitrogen_g[i], pools.phosphorus_g[i] }) |v| if (!std.math.isFinite(v)) return error.NonFiniteExudationPool;
    };
}

/// Direct translation of redist.f lines 6115--6119 (inner body of DO 125 L loop).
pub fn apply(pools: DissolvedOrganicMicropores, fluxes: ExudationFluxes) !DissolvedOrganicMicropores {
    for (0..5) |k| {
        if (!std.math.isFinite(pools.oqc[k]) or !std.math.isFinite(pools.oqn[k]) or !std.math.isFinite(pools.oqp[k])) return error.InvalidExudationPool;
        if (!std.math.isFinite(fluxes.tdfomc[k])) return error.InvalidExudationFlux;
        if (!std.math.isFinite(fluxes.tdfomn[k])) return error.InvalidExudationFlux;
        if (!std.math.isFinite(fluxes.tdfomp[k])) return error.InvalidExudationFlux;
    }

    var new = pools;
    for (0..5) |k| {
        new.oqc[k] = pools.oqc[k] + fluxes.tdfomc[k];
        new.oqn[k] = pools.oqn[k] + fluxes.tdfomn[k];
        new.oqp[k] = pools.oqp[k] + fluxes.tdfomp[k];
    }
    for (0..5) |k| if (!std.math.isFinite(new.oqc[k]) or !std.math.isFinite(new.oqn[k]) or !std.math.isFinite(new.oqp[k])) return error.NonFiniteExudationPool;
    return new;
}

test "REDIST root exudation adds to OQC across all K fractions" {
    var fluxes = std.mem.zeroes(ExudationFluxes);
    fluxes.tdfomc[0] = 1.0;
    fluxes.tdfomc[4] = 0.5;
    const result = try apply(std.mem.zeroes(DissolvedOrganicMicropores), fluxes);
    try std.testing.expectApproxEqRel(@as(f64, 1.0), result.oqc[0], 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 0.5), result.oqc[4], 1.0e-15);
}

test "REDIST root exudation OQN and OQP are independent of OQC" {
    var fluxes = std.mem.zeroes(ExudationFluxes);
    fluxes.tdfomn[2] = 0.3;
    fluxes.tdfomp[3] = 0.1;
    const result = try apply(std.mem.zeroes(DissolvedOrganicMicropores), fluxes);
    try std.testing.expectApproxEqRel(@as(f64, 0.3), result.oqn[2], 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 0.1), result.oqp[3], 1.0e-15);
    try std.testing.expectEqual(@as(f64, 0.0), result.oqc[2]);
}

test "REDIST root exudation rejects non-finite flux" {
    var bad = std.mem.zeroes(ExudationFluxes);
    bad.tdfomn[1] = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidExudationFlux,
        apply(std.mem.zeroes(DissolvedOrganicMicropores), bad),
    );
}

test "REDIST root exudation supports runtime layers and classes" {
    var c = [_]f64{0} ** 6;
    var n = [_]f64{0} ** 6;
    var p = [_]f64{0} ** 6;
    const cf = [_]f64{1} ** 6;
    const nf = [_]f64{2} ** 6;
    const pf = [_]f64{-1} ** 6;
    try applyLayers(2, 3, .{ .carbon_g = &c, .nitrogen_g = &n, .phosphorus_g = &p }, .{ .carbon_g = &cf, .nitrogen_g = &nf, .phosphorus_g = &pf });
    try std.testing.expectEqual(@as(f64, 1), c[5]);
    try std.testing.expectEqual(@as(f64, 2), n[5]);
    try std.testing.expectEqual(@as(f64, -1), p[5]);
}

test "REDIST root exudation rejects non-carbon pool overflow" {
    var pools = std.mem.zeroes(DissolvedOrganicMicropores);
    pools.oqn[0] = std.math.floatMax(f64);
    var fluxes = std.mem.zeroes(ExudationFluxes);
    fluxes.tdfomn[0] = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteExudationPool, apply(pools, fluxes));
}
