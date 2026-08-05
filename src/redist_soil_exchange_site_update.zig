const std = @import("std");

/// Soil exchangeable cation and anion/phosphate pools (mol). 19 fields, REDIST 6320-6352.
pub const ExchangeSitePools = struct {
    // Non-band (lines 6320-6328)
    nh4_nonband_mol: f64, // XN4
    nh4_band_mol: f64, // XNB
    h_mol: f64, // XHY
    al_mol: f64, // XAL
    fe_mol: f64, // XFE
    ca_mol: f64, // XCA
    megagrams_mol: f64, // XMG
    na_mol: f64, // XNA
    ka_mol: f64, // XKA
    // Anion sites non-band (lines 6343-6347)
    oh0_mol: f64, // XOH0
    oh1_mol: f64, // XOH1
    oh2_mol: f64, // XOH2
    hpo4_mol: f64, // XH1P
    h2po4_mol: f64, // XH2P
    // Anion sites band (lines 6348-6352)
    oh0_band_mol: f64, // XOH0B
    oh1_band_mol: f64, // XOH1B
    oh2_band_mol: f64, // XOH2B
    hpo4_band_mol: f64, // XH1PB
    h2po4_band_mol: f64, // XH2PB
};

/// Exchange transformation fluxes (mol step-1). 19 fields matching ExchangeSitePools.
pub const ExchangeSiteFluxes = struct {
    nh4_nonband_mol: f64, // TRXN4
    nh4_band_mol: f64, // TRXNB
    h_mol: f64, // TRXHY
    al_mol: f64, // TRXAL
    fe_mol: f64, // TRXFE
    ca_mol: f64, // TRXCA
    megagrams_mol: f64, // TRXMG
    na_mol: f64, // TRXNA
    ka_mol: f64, // TRXKA
    oh0_mol: f64, // TRXH0
    oh1_mol: f64, // TRXH1
    oh2_mol: f64, // TRXH2
    hpo4_mol: f64, // TRX1P
    h2po4_mol: f64, // TRX2P
    oh0_band_mol: f64, // TRBH0
    oh1_band_mol: f64, // TRBH1
    oh2_band_mol: f64, // TRBH2
    hpo4_band_mol: f64, // TRB1P
    h2po4_band_mol: f64, // TRB2P
};

/// Direct translation of REDIST lines 6320--6352 (inner body of DO 125 L loop).
pub fn update(pools: ExchangeSitePools, fluxes: ExchangeSiteFluxes) !ExchangeSitePools {
    inline for (@typeInfo(ExchangeSitePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(pools, field.name)))
            return error.InvalidExchangeSitePool;
    inline for (@typeInfo(ExchangeSiteFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(fluxes, field.name)))
            return error.InvalidExchangeSiteFlux;

    const result = ExchangeSitePools{
        .nh4_nonband_mol = pools.nh4_nonband_mol + fluxes.nh4_nonband_mol,
        .nh4_band_mol = pools.nh4_band_mol + fluxes.nh4_band_mol,
        .h_mol = pools.h_mol + fluxes.h_mol,
        .al_mol = pools.al_mol + fluxes.al_mol,
        .fe_mol = pools.fe_mol + fluxes.fe_mol,
        .ca_mol = pools.ca_mol + fluxes.ca_mol,
        .megagrams_mol = pools.megagrams_mol + fluxes.megagrams_mol,
        .na_mol = pools.na_mol + fluxes.na_mol,
        .ka_mol = pools.ka_mol + fluxes.ka_mol,
        .oh0_mol = pools.oh0_mol + fluxes.oh0_mol,
        .oh1_mol = pools.oh1_mol + fluxes.oh1_mol,
        .oh2_mol = pools.oh2_mol + fluxes.oh2_mol,
        .hpo4_mol = pools.hpo4_mol + fluxes.hpo4_mol,
        .h2po4_mol = pools.h2po4_mol + fluxes.h2po4_mol,
        .oh0_band_mol = pools.oh0_band_mol + fluxes.oh0_band_mol,
        .oh1_band_mol = pools.oh1_band_mol + fluxes.oh1_band_mol,
        .oh2_band_mol = pools.oh2_band_mol + fluxes.oh2_band_mol,
        .hpo4_band_mol = pools.hpo4_band_mol + fluxes.hpo4_band_mol,
        .h2po4_band_mol = pools.h2po4_band_mol + fluxes.h2po4_band_mol,
    };
    inline for (@typeInfo(ExchangeSitePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteExchangeSitePool;
    return result;
}

/// Applies REDIST exchange reactions to every runtime-configured soil layer.
pub fn updateLayers(
    pools_by_layer: []ExchangeSitePools,
    fluxes_by_layer: []const ExchangeSiteFluxes,
) !void {
    if (pools_by_layer.len == 0 or fluxes_by_layer.len != pools_by_layer.len)
        return error.ExchangeSiteDimensionMismatch;

    for (pools_by_layer, fluxes_by_layer) |*pools, fluxes|
        pools.* = try update(pools.*, fluxes);
}

test "REDIST soil exchange site cations non-band and band are independent" {
    var f = std.mem.zeroes(ExchangeSiteFluxes);
    f.nh4_nonband_mol = 1.0;
    f.nh4_band_mol = 2.0;
    const result = try update(std.mem.zeroes(ExchangeSitePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 1.0), result.nh4_nonband_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.nh4_band_mol, 1.0e-15);
}

test "REDIST soil exchange site anion sites non-band and band are independent" {
    var f = std.mem.zeroes(ExchangeSiteFluxes);
    f.oh1_mol = 0.5;
    f.oh1_band_mol = 0.3;
    f.hpo4_band_mol = 0.1;
    const result = try update(std.mem.zeroes(ExchangeSitePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 0.5), result.oh1_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 0.3), result.oh1_band_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 0.1), result.hpo4_band_mol, 1.0e-15);
}

test "REDIST soil exchange site rejects non-finite flux" {
    var bad = std.mem.zeroes(ExchangeSiteFluxes);
    bad.h2po4_band_mol = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidExchangeSiteFlux,
        update(std.mem.zeroes(ExchangeSitePools), bad),
    );
}

test "REDIST soil exchange site maps all source transformations in order" {
    var fluxes: ExchangeSiteFluxes = undefined;
    inline for (@typeInfo(ExchangeSiteFluxes).@"struct".fields, 1..) |field, value|
        @field(fluxes, field.name) = @floatFromInt(value);

    const result = try update(std.mem.zeroes(ExchangeSitePools), fluxes);
    inline for (@typeInfo(ExchangeSitePools).@"struct".fields, 1..) |field, expected|
        try std.testing.expectEqual(@as(f64, @floatFromInt(expected)), @field(result, field.name));
}

test "REDIST soil exchange site runtime layers remain independent" {
    var pools = [_]ExchangeSitePools{
        std.mem.zeroes(ExchangeSitePools),
        std.mem.zeroes(ExchangeSitePools),
    };
    var fluxes = [_]ExchangeSiteFluxes{
        std.mem.zeroes(ExchangeSiteFluxes),
        std.mem.zeroes(ExchangeSiteFluxes),
    };
    fluxes[0].al_mol = 2.0;
    fluxes[1].h2po4_band_mol = 3.0;

    try updateLayers(&pools, &fluxes);

    try std.testing.expectEqual(@as(f64, 2.0), pools[0].al_mol);
    try std.testing.expectEqual(@as(f64, 0.0), pools[1].al_mol);
    try std.testing.expectEqual(@as(f64, 0.0), pools[0].h2po4_band_mol);
    try std.testing.expectEqual(@as(f64, 3.0), pools[1].h2po4_band_mol);
}

test "REDIST soil exchange site runtime layers reject dimensions and overflow" {
    var pools = [_]ExchangeSitePools{std.mem.zeroes(ExchangeSitePools)};
    const no_fluxes: [0]ExchangeSiteFluxes = .{};
    try std.testing.expectError(
        error.ExchangeSiteDimensionMismatch,
        updateLayers(&pools, &no_fluxes),
    );

    pools[0].ca_mol = std.math.floatMax(f64);
    var fluxes = std.mem.zeroes(ExchangeSiteFluxes);
    fluxes.ca_mol = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteExchangeSitePool, update(pools[0], fluxes));
}
