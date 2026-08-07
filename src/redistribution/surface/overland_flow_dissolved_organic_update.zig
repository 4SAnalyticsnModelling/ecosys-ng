const std = @import("std");

/// Litter dissolved organic pool for one fraction K (g).
pub const OrganicPool = struct {
    /// OQC(K,0). DOC content.
    doc_g: f64,
    /// OQN(K,0). DON content.
    don_g: f64,
    /// OQP(K,0). DOP content.
    dop_g: f64,
    /// OQA(K,0). Acetate content.
    acetate_g: f64,
};

/// Overland flow organic runoff flux for one fraction K (g step-1).
pub const OverlandFlux = struct {
    /// TOCQRS(K). DOC overland flux.
    doc_g: f64,
    /// TONQRS(K). DON overland flux.
    don_g: f64,
    /// TOPQRS(K). DOP overland flux.
    dop_g: f64,
    /// TOAQRS(K). Acetate overland flux.
    acetate_g: f64,
};

/// Direct translation of redist.f lines 5065--5075 (K=0,2 loop body).
///
/// Caller must check `ABS(TQR) > ZEROS` before invoking.
/// Updates OQC/OQN/OQP/OQA in each of the three organic fractions.
pub fn update(
    fluxes: [3]OverlandFlux,
    pools: [3]OrganicPool,
) ![3]OrganicPool {
    for (0..3) |k| {
        inline for (@typeInfo(OverlandFlux).@"struct".fields) |field|
            if (!std.math.isFinite(@field(fluxes[k], field.name)))
                return error.InvalidOverlandOrganicFlux;
        inline for (@typeInfo(OrganicPool).@"struct".fields) |field|
            if (!std.math.isFinite(@field(pools[k], field.name)))
                return error.InvalidOverlandOrganicPool;
    }

    var out: [3]OrganicPool = undefined;
    for (0..3) |k| {
        out[k] = .{
            .doc_g = pools[k].doc_g + fluxes[k].doc_g,
            .don_g = pools[k].don_g + fluxes[k].don_g,
            .dop_g = pools[k].dop_g + fluxes[k].dop_g,
            .acetate_g = pools[k].acetate_g + fluxes[k].acetate_g,
        };
        inline for (@typeInfo(OrganicPool).@"struct".fields) |field|
            if (!std.math.isFinite(@field(out[k], field.name)))
                return error.NonFiniteOverlandOrganicPool;
    }
    return out;
}

test "REDIST overland organic DOC accumulates into all three fractions" {
    const fluxes: [3]OverlandFlux = .{
        .{ .doc_g = 1.0, .don_g = 0.0, .dop_g = 0.0, .acetate_g = 0.0 },
        .{ .doc_g = 2.0, .don_g = 0.0, .dop_g = 0.0, .acetate_g = 0.0 },
        .{ .doc_g = 3.0, .don_g = 0.0, .dop_g = 0.0, .acetate_g = 0.0 },
    };
    const pools = [3]OrganicPool{
        std.mem.zeroes(OrganicPool),
        std.mem.zeroes(OrganicPool),
        std.mem.zeroes(OrganicPool),
    };
    const result = try update(fluxes, pools);
    try std.testing.expectEqual(@as(f64, 1.0), result[0].doc_g);
    try std.testing.expectEqual(@as(f64, 2.0), result[1].doc_g);
    try std.testing.expectEqual(@as(f64, 3.0), result[2].doc_g);
}

test "REDIST overland organic fractions are independent" {
    const fluxes: [3]OverlandFlux = .{
        std.mem.zeroes(OverlandFlux),
        .{ .doc_g = 5.0, .don_g = 0.5, .dop_g = 0.05, .acetate_g = 0.5 },
        std.mem.zeroes(OverlandFlux),
    };
    const pools = [3]OrganicPool{
        .{ .doc_g = 10.0, .don_g = 0.0, .dop_g = 0.0, .acetate_g = 0.0 },
        std.mem.zeroes(OrganicPool),
        std.mem.zeroes(OrganicPool),
    };
    const result = try update(fluxes, pools);
    try std.testing.expectEqual(@as(f64, 10.0), result[0].doc_g);
    try std.testing.expectEqual(@as(f64, 5.0), result[1].doc_g);
    try std.testing.expectEqual(@as(f64, 0.0), result[2].doc_g);
}

test "REDIST overland organic rejects non-finite flux" {
    const bad: [3]OverlandFlux = .{
        .{ .doc_g = std.math.nan(f64), .don_g = 0.0, .dop_g = 0.0, .acetate_g = 0.0 },
        std.mem.zeroes(OverlandFlux),
        std.mem.zeroes(OverlandFlux),
    };
    const pools = [3]OrganicPool{ std.mem.zeroes(OrganicPool), std.mem.zeroes(OrganicPool), std.mem.zeroes(OrganicPool) };
    try std.testing.expectError(error.InvalidOverlandOrganicFlux, update(bad, pools));
}
