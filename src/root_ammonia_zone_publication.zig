const std = @import("std");

pub const ZoneTotals = struct {
    non_band_ammonia_g_n_per_step: f64,
    band_ammonia_g_n_per_step: f64,
};

pub const Contributions = struct {
    /// Root and mycorrhizal contributions for one soil layer, in source order.
    non_band_ammonia_g_n_per_step: []const f64,
    band_ammonia_g_n_per_step: []const f64,
};

/// Direct translation of EXTRACT lines 796--797 for one runtime soil layer.
/// Unlike the coalesced production diagnostic, this compatibility owner keeps
/// non-band and band NH3 exchange distinct through publication.
pub fn accumulateLayer(
    preceding: ZoneTotals,
    contributions: Contributions,
) !ZoneTotals {
    if (contributions.non_band_ammonia_g_n_per_step.len == 0 or
        contributions.non_band_ammonia_g_n_per_step.len !=
            contributions.band_ammonia_g_n_per_step.len)
        return error.RootAmmoniaZoneDimensionMismatch;
    inline for (std.meta.fields(ZoneTotals)) |field|
        if (!std.math.isFinite(@field(preceding, field.name)))
            return error.NonFiniteRootAmmoniaZoneTotal;

    var next = preceding;
    for (
        contributions.non_band_ammonia_g_n_per_step,
        contributions.band_ammonia_g_n_per_step,
    ) |non_band, band| {
        if (!std.math.isFinite(non_band) or !std.math.isFinite(band))
            return error.NonFiniteRootAmmoniaZoneContribution;
        next.non_band_ammonia_g_n_per_step += non_band;
        next.band_ammonia_g_n_per_step += band;
        if (!std.math.isFinite(next.non_band_ammonia_g_n_per_step) or
            !std.math.isFinite(next.band_ammonia_g_n_per_step))
            return error.NonFiniteRootAmmoniaZoneResult;
    }
    return next;
}

test "EXTRACT ammonia publication retains separate non-band and band totals" {
    const result = try accumulateLayer(.{
        .non_band_ammonia_g_n_per_step = 1,
        .band_ammonia_g_n_per_step = 10,
    }, .{
        .non_band_ammonia_g_n_per_step = &.{ 2, 3, 4 },
        .band_ammonia_g_n_per_step = &.{ 20, 30, 40 },
    });
    try std.testing.expectEqual(
        @as(f64, 10),
        result.non_band_ammonia_g_n_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 100),
        result.band_ammonia_g_n_per_step,
    );
}

test "EXTRACT ammonia publication preserves signed source exchanges" {
    const result = try accumulateLayer(.{
        .non_band_ammonia_g_n_per_step = 5,
        .band_ammonia_g_n_per_step = 7,
    }, .{
        .non_band_ammonia_g_n_per_step = &.{ -2, 1 },
        .band_ammonia_g_n_per_step = &.{ 3, -4 },
    });
    try std.testing.expectEqual(
        @as(f64, 4),
        result.non_band_ammonia_g_n_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 6),
        result.band_ammonia_g_n_per_step,
    );
}

test "EXTRACT ammonia publication rejects late non-finite contribution" {
    try std.testing.expectError(
        error.NonFiniteRootAmmoniaZoneContribution,
        accumulateLayer(.{
            .non_band_ammonia_g_n_per_step = 5,
            .band_ammonia_g_n_per_step = 7,
        }, .{
            .non_band_ammonia_g_n_per_step = &.{ 1, std.math.nan(f64) },
            .band_ammonia_g_n_per_step = &.{ 2, 3 },
        }),
    );
}
