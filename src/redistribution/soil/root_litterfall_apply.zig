const std = @import("std");

/// SOC/SON/SOP pools for one soil layer, K=0-1 (litter types), M=1-5 (size classes).
pub const SocLayer = struct {
    /// OSC[K=0-1][M=1-5] (g C).
    osc: [2][5]f64,
    /// OSN[K=0-1][M=1-5] (g N).
    osn: [2][5]f64,
    /// OSP[K=0-1][M=1-5] (g P).
    osp: [2][5]f64,
};

/// Root litterfall fluxes from extract.f for one soil layer (g C,N,P step-1).
pub const RootLitterfallFluxes = struct {
    /// CSNT[K=0-1][M=1-5].
    csnt: [2][5]f64,
    /// ZSNT[K=0-1][M=1-5].
    zsnt: [2][5]f64,
    /// PSNT[K=0-1][M=1-5].
    psnt: [2][5]f64,
};

/// Direct translation of redist.f lines 6057--6069 (inner body of DO 125 L loop).
///
/// OSA (colonised SOC) is commented out in source (line 6060); not updated here.
pub fn apply(layer: SocLayer, fluxes: RootLitterfallFluxes) !SocLayer {
    for (0..2) |k| {
        for (0..5) |m| {
            if (!std.math.isFinite(layer.osc[k][m])) return error.InvalidSocPool;
            if (!std.math.isFinite(fluxes.csnt[k][m])) return error.InvalidRootLitterfallFlux;
        }
    }

    var new = layer;
    for (0..2) |k| {
        for (0..5) |m| {
            new.osc[k][m] = layer.osc[k][m] + fluxes.csnt[k][m];
            new.osn[k][m] = layer.osn[k][m] + fluxes.zsnt[k][m];
            new.osp[k][m] = layer.osp[k][m] + fluxes.psnt[k][m];
        }
    }
    for (0..2) |k| {
        for (0..5) |m| {
            if (!std.math.isFinite(new.osc[k][m])) return error.NonFiniteSocPool;
        }
    }
    return new;
}

test "REDIST root litterfall OSC accumulates across K and M dimensions" {
    var fluxes = std.mem.zeroes(RootLitterfallFluxes);
    fluxes.csnt[0][0] = 1.0;
    fluxes.csnt[1][4] = 2.0;
    const result = try apply(std.mem.zeroes(SocLayer), fluxes);
    try std.testing.expectApproxEqRel(@as(f64, 1.0), result.osc[0][0], 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.osc[1][4], 1.0e-15);
}

test "REDIST root litterfall OSN and OSP update independently from OSC" {
    var fluxes = std.mem.zeroes(RootLitterfallFluxes);
    fluxes.zsnt[0][2] = 0.5;
    fluxes.psnt[1][3] = 0.3;
    const result = try apply(std.mem.zeroes(SocLayer), fluxes);
    try std.testing.expectApproxEqRel(@as(f64, 0.5), result.osn[0][2], 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 0.3), result.osp[1][3], 1.0e-15);
    try std.testing.expectEqual(@as(f64, 0.0), result.osc[0][2]);
}

test "REDIST root litterfall rejects non-finite flux" {
    var bad = std.mem.zeroes(RootLitterfallFluxes);
    bad.csnt[1][2] = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidRootLitterfallFlux,
        apply(std.mem.zeroes(SocLayer), bad),
    );
}
