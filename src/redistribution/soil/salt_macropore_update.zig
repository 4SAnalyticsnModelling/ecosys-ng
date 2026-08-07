const std = @import("std");

/// Soil salt macropore pools (mol). 50 fields, REDIST 7114-7162.
/// Gated on ISALTG != 0 by the caller.
pub const SaltMacroporePools = struct {
    al_mol: f64, // ZALH
    fe_mol: f64, // ZFEH
    h_mol: f64, // ZHYH
    ca_mol: f64, // ZCCH
    megagrams_mol: f64, // ZMAH
    na_mol: f64, // ZNAH
    ka_mol: f64, // ZKAH
    oh_mol: f64, // ZOHH
    so4_mol: f64, // ZSO4H
    cl_mol: f64, // ZCLH
    co3_mol: f64, // ZCO3H
    hco3_mol: f64, // ZHCO3H
    aloh1_mol: f64, // ZALO1H
    aloh2_mol: f64, // ZALO2H
    aloh3_mol: f64, // ZALO3H
    aloh4_mol: f64, // ZALO4H
    als_mol: f64, // ZALSH
    feoh1_mol: f64, // ZFEO1H
    feoh2_mol: f64, // ZFEO2H
    feoh3_mol: f64, // ZFEO3H
    feoh4_mol: f64, // ZFEO4H
    fes_mol: f64, // ZFESH
    cao_mol: f64, // ZCAOH
    cac_mol: f64, // ZCACH
    cah_mol: f64, // ZCAHH
    cas_mol: f64, // ZCASH
    mgo_mol: f64, // ZMGOH
    mgc_mol: f64, // ZMGCH
    mgh_mol: f64, // ZMGHH
    mgs_mol: f64, // ZMGSH
    nac_mol: f64, // ZNACH
    nas_mol: f64, // ZNASH
    kas_mol: f64, // ZKASH
    h0po4_mol: f64, // H0PO4H
    h3po4_mol: f64, // H3PO4H
    fe1p_mol: f64, // ZFE1PH
    fe2p_mol: f64, // ZFE2PH
    ca0p_mol: f64, // ZCA0PH
    ca1p_mol: f64, // ZCA1PH
    ca2p_mol: f64, // ZCA2PH
    mg1p_mol: f64, // ZMG1PH
    h0pob_mol: f64, // H0POBH
    h3pob_mol: f64, // H3POBH
    fe1pb_mol: f64, // ZFE1BH
    fe2pb_mol: f64, // ZFE2BH
    ca0pb_mol: f64, // ZCA0BH
    ca1pb_mol: f64, // ZCA1BH
    ca2pb_mol: f64, // ZCA2BH
    mg1pb_mol: f64, // ZMG1BH
};

/// Macropore transport and exchange sign: pool += transport - exchange.
/// Each species has: macropore_transport_mol (T*FHS/T*FHB) and exchange_mol (X*FXS/X*FXB).
pub const SaltMacroporeFluxes = struct {
    al_macropore_mol: f64,
    al_exchange_mol: f64,
    fe_macropore_mol: f64,
    fe_exchange_mol: f64,
    h_macropore_mol: f64,
    h_exchange_mol: f64,
    ca_macropore_mol: f64,
    ca_exchange_mol: f64,
    megagrams_macropore_mol: f64,
    megagrams_exchange_mol: f64,
    na_macropore_mol: f64,
    na_exchange_mol: f64,
    ka_macropore_mol: f64,
    ka_exchange_mol: f64,
    oh_macropore_mol: f64,
    oh_exchange_mol: f64,
    so4_macropore_mol: f64,
    so4_exchange_mol: f64,
    cl_macropore_mol: f64,
    cl_exchange_mol: f64,
    co3_macropore_mol: f64,
    co3_exchange_mol: f64,
    hco3_macropore_mol: f64,
    hco3_exchange_mol: f64,
    aloh1_macropore_mol: f64,
    aloh1_exchange_mol: f64,
    aloh2_macropore_mol: f64,
    aloh2_exchange_mol: f64,
    aloh3_macropore_mol: f64,
    aloh3_exchange_mol: f64,
    aloh4_macropore_mol: f64,
    aloh4_exchange_mol: f64,
    als_macropore_mol: f64,
    als_exchange_mol: f64,
    feoh1_macropore_mol: f64,
    feoh1_exchange_mol: f64,
    feoh2_macropore_mol: f64,
    feoh2_exchange_mol: f64,
    feoh3_macropore_mol: f64,
    feoh3_exchange_mol: f64,
    feoh4_macropore_mol: f64,
    feoh4_exchange_mol: f64,
    fes_macropore_mol: f64,
    fes_exchange_mol: f64,
    cao_macropore_mol: f64,
    cao_exchange_mol: f64,
    cac_macropore_mol: f64,
    cac_exchange_mol: f64,
    cah_macropore_mol: f64,
    cah_exchange_mol: f64,
    cas_macropore_mol: f64,
    cas_exchange_mol: f64,
    mgo_macropore_mol: f64,
    mgo_exchange_mol: f64,
    mgc_macropore_mol: f64,
    mgc_exchange_mol: f64,
    mgh_macropore_mol: f64,
    mgh_exchange_mol: f64,
    mgs_macropore_mol: f64,
    mgs_exchange_mol: f64,
    nac_macropore_mol: f64,
    nac_exchange_mol: f64,
    nas_macropore_mol: f64,
    nas_exchange_mol: f64,
    kas_macropore_mol: f64,
    kas_exchange_mol: f64,
    h0po4_macropore_mol: f64,
    h0po4_exchange_mol: f64,
    h3po4_macropore_mol: f64,
    h3po4_exchange_mol: f64,
    fe1p_macropore_mol: f64,
    fe1p_exchange_mol: f64,
    fe2p_macropore_mol: f64,
    fe2p_exchange_mol: f64,
    ca0p_macropore_mol: f64,
    ca0p_exchange_mol: f64,
    ca1p_macropore_mol: f64,
    ca1p_exchange_mol: f64,
    ca2p_macropore_mol: f64,
    ca2p_exchange_mol: f64,
    mg1p_macropore_mol: f64,
    mg1p_exchange_mol: f64,
    h0pob_macropore_mol: f64,
    h0pob_exchange_mol: f64,
    h3pob_macropore_mol: f64,
    h3pob_exchange_mol: f64,
    fe1pb_macropore_mol: f64,
    fe1pb_exchange_mol: f64,
    fe2pb_macropore_mol: f64,
    fe2pb_exchange_mol: f64,
    ca0pb_macropore_mol: f64,
    ca0pb_exchange_mol: f64,
    ca1pb_macropore_mol: f64,
    ca1pb_exchange_mol: f64,
    ca2pb_macropore_mol: f64,
    ca2pb_exchange_mol: f64,
    mg1pb_macropore_mol: f64,
    mg1pb_exchange_mol: f64,
};

pub const CarbonateExchangeSite = struct {
    hco3_mol: f64, // XHC
};

pub const CarbonateExchangeFlux = struct {
    hco3_mol: f64, // TRXHC
};

pub const LayerResult = struct {
    macropore: SaltMacroporePools,
    carbonate_exchange: CarbonateExchangeSite,
};

pub const SaltEquilibriumMode = enum { static, dynamic };

/// Direct translation of redist.f lines 7114--7162.
/// All macropore pools follow: pool += transport - exchange.
pub fn update(pools: SaltMacroporePools, f: SaltMacroporeFluxes) !SaltMacroporePools {
    inline for (@typeInfo(SaltMacroporePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(pools, field.name)))
            return error.InvalidSaltMacroporePool;
    inline for (@typeInfo(SaltMacroporeFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(f, field.name)))
            return error.InvalidSaltMacroporeFlux;

    const r = SaltMacroporePools{
        .al_mol = pools.al_mol + f.al_macropore_mol - f.al_exchange_mol,
        .fe_mol = pools.fe_mol + f.fe_macropore_mol - f.fe_exchange_mol,
        .h_mol = pools.h_mol + f.h_macropore_mol - f.h_exchange_mol,
        .ca_mol = pools.ca_mol + f.ca_macropore_mol - f.ca_exchange_mol,
        .megagrams_mol = pools.megagrams_mol + f.megagrams_macropore_mol - f.megagrams_exchange_mol,
        .na_mol = pools.na_mol + f.na_macropore_mol - f.na_exchange_mol,
        .ka_mol = pools.ka_mol + f.ka_macropore_mol - f.ka_exchange_mol,
        .oh_mol = pools.oh_mol + f.oh_macropore_mol - f.oh_exchange_mol,
        .so4_mol = pools.so4_mol + f.so4_macropore_mol - f.so4_exchange_mol,
        .cl_mol = pools.cl_mol + f.cl_macropore_mol - f.cl_exchange_mol,
        .co3_mol = pools.co3_mol + f.co3_macropore_mol - f.co3_exchange_mol,
        .hco3_mol = pools.hco3_mol + f.hco3_macropore_mol - f.hco3_exchange_mol,
        .aloh1_mol = pools.aloh1_mol + f.aloh1_macropore_mol - f.aloh1_exchange_mol,
        .aloh2_mol = pools.aloh2_mol + f.aloh2_macropore_mol - f.aloh2_exchange_mol,
        .aloh3_mol = pools.aloh3_mol + f.aloh3_macropore_mol - f.aloh3_exchange_mol,
        .aloh4_mol = pools.aloh4_mol + f.aloh4_macropore_mol - f.aloh4_exchange_mol,
        .als_mol = pools.als_mol + f.als_macropore_mol - f.als_exchange_mol,
        .feoh1_mol = pools.feoh1_mol + f.feoh1_macropore_mol - f.feoh1_exchange_mol,
        .feoh2_mol = pools.feoh2_mol + f.feoh2_macropore_mol - f.feoh2_exchange_mol,
        .feoh3_mol = pools.feoh3_mol + f.feoh3_macropore_mol - f.feoh3_exchange_mol,
        .feoh4_mol = pools.feoh4_mol + f.feoh4_macropore_mol - f.feoh4_exchange_mol,
        .fes_mol = pools.fes_mol + f.fes_macropore_mol - f.fes_exchange_mol,
        .cao_mol = pools.cao_mol + f.cao_macropore_mol - f.cao_exchange_mol,
        .cac_mol = pools.cac_mol + f.cac_macropore_mol - f.cac_exchange_mol,
        .cah_mol = pools.cah_mol + f.cah_macropore_mol - f.cah_exchange_mol,
        .cas_mol = pools.cas_mol + f.cas_macropore_mol - f.cas_exchange_mol,
        .mgo_mol = pools.mgo_mol + f.mgo_macropore_mol - f.mgo_exchange_mol,
        .mgc_mol = pools.mgc_mol + f.mgc_macropore_mol - f.mgc_exchange_mol,
        .mgh_mol = pools.mgh_mol + f.mgh_macropore_mol - f.mgh_exchange_mol,
        .mgs_mol = pools.mgs_mol + f.mgs_macropore_mol - f.mgs_exchange_mol,
        .nac_mol = pools.nac_mol + f.nac_macropore_mol - f.nac_exchange_mol,
        .nas_mol = pools.nas_mol + f.nas_macropore_mol - f.nas_exchange_mol,
        .kas_mol = pools.kas_mol + f.kas_macropore_mol - f.kas_exchange_mol,
        .h0po4_mol = pools.h0po4_mol + f.h0po4_macropore_mol - f.h0po4_exchange_mol,
        .h3po4_mol = pools.h3po4_mol + f.h3po4_macropore_mol - f.h3po4_exchange_mol,
        .fe1p_mol = pools.fe1p_mol + f.fe1p_macropore_mol - f.fe1p_exchange_mol,
        .fe2p_mol = pools.fe2p_mol + f.fe2p_macropore_mol - f.fe2p_exchange_mol,
        .ca0p_mol = pools.ca0p_mol + f.ca0p_macropore_mol - f.ca0p_exchange_mol,
        .ca1p_mol = pools.ca1p_mol + f.ca1p_macropore_mol - f.ca1p_exchange_mol,
        .ca2p_mol = pools.ca2p_mol + f.ca2p_macropore_mol - f.ca2p_exchange_mol,
        .mg1p_mol = pools.mg1p_mol + f.mg1p_macropore_mol - f.mg1p_exchange_mol,
        .h0pob_mol = pools.h0pob_mol + f.h0pob_macropore_mol - f.h0pob_exchange_mol,
        .h3pob_mol = pools.h3pob_mol + f.h3pob_macropore_mol - f.h3pob_exchange_mol,
        .fe1pb_mol = pools.fe1pb_mol + f.fe1pb_macropore_mol - f.fe1pb_exchange_mol,
        .fe2pb_mol = pools.fe2pb_mol + f.fe2pb_macropore_mol - f.fe2pb_exchange_mol,
        .ca0pb_mol = pools.ca0pb_mol + f.ca0pb_macropore_mol - f.ca0pb_exchange_mol,
        .ca1pb_mol = pools.ca1pb_mol + f.ca1pb_macropore_mol - f.ca1pb_exchange_mol,
        .ca2pb_mol = pools.ca2pb_mol + f.ca2pb_macropore_mol - f.ca2pb_exchange_mol,
        .mg1pb_mol = pools.mg1pb_mol + f.mg1pb_macropore_mol - f.mg1pb_exchange_mol,
    };
    inline for (@typeInfo(SaltMacroporePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(r, field.name)))
            return error.NonFiniteSaltMacroporePool;
    return r;
}

/// Direct translation of REDIST 7114--7163 including the XHC exchange pool.
pub fn updateLayer(
    pools: SaltMacroporePools,
    carbonate_exchange: CarbonateExchangeSite,
    fluxes: SaltMacroporeFluxes,
    carbonate_flux: CarbonateExchangeFlux,
) !LayerResult {
    if (!std.math.isFinite(carbonate_exchange.hco3_mol))
        return error.InvalidCarbonateExchangePool;
    if (!std.math.isFinite(carbonate_flux.hco3_mol))
        return error.InvalidCarbonateExchangeFlux;
    const macropore = try update(pools, fluxes);
    const new_exchange = CarbonateExchangeSite{
        .hco3_mol = carbonate_exchange.hco3_mol + carbonate_flux.hco3_mol,
    };
    if (!std.math.isFinite(new_exchange.hco3_mol))
        return error.NonFiniteCarbonateExchangePool;
    return .{ .macropore = macropore, .carbonate_exchange = new_exchange };
}

pub fn updateLayers(
    pools_by_layer: []SaltMacroporePools,
    carbonate_exchange_by_layer: []CarbonateExchangeSite,
    fluxes_by_layer: []const SaltMacroporeFluxes,
    carbonate_fluxes_by_layer: []const CarbonateExchangeFlux,
    mode: SaltEquilibriumMode,
) !void {
    if (pools_by_layer.len == 0 or
        carbonate_exchange_by_layer.len != pools_by_layer.len or
        fluxes_by_layer.len != pools_by_layer.len or
        carbonate_fluxes_by_layer.len != pools_by_layer.len)
        return error.SaltMacroporeDimensionMismatch;
    if (mode == .static) return;

    for (
        pools_by_layer,
        carbonate_exchange_by_layer,
        fluxes_by_layer,
        carbonate_fluxes_by_layer,
    ) |*pools, *carbonate_exchange, fluxes, carbonate_flux| {
        const result = try updateLayer(pools.*, carbonate_exchange.*, fluxes, carbonate_flux);
        pools.* = result.macropore;
        carbonate_exchange.* = result.carbonate_exchange;
    }
}

test "REDIST soil salt macropore transport adds and exchange subtracts" {
    var f = std.mem.zeroes(SaltMacroporeFluxes);
    f.ca_macropore_mol = 3.0;
    f.ca_exchange_mol = 1.0;
    const result = try update(std.mem.zeroes(SaltMacroporePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.ca_mol, 1.0e-15);
}

test "REDIST soil salt macropore non-band and band phosphate are independent" {
    var f = std.mem.zeroes(SaltMacroporeFluxes);
    f.ca1p_macropore_mol = 1.0;
    f.ca1pb_macropore_mol = 2.0;
    const result = try update(std.mem.zeroes(SaltMacroporePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 1.0), result.ca1p_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.ca1pb_mol, 1.0e-15);
}

test "REDIST soil salt macropore rejects non-finite flux" {
    var bad = std.mem.zeroes(SaltMacroporeFluxes);
    bad.feoh3_macropore_mol = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSaltMacroporeFlux,
        update(std.mem.zeroes(SaltMacroporePools), bad),
    );
}

test "REDIST soil salt macropore includes XHC transformation" {
    const result = try updateLayer(
        std.mem.zeroes(SaltMacroporePools),
        .{ .hco3_mol = 2.0 },
        std.mem.zeroes(SaltMacroporeFluxes),
        .{ .hco3_mol = 3.0 },
    );
    try std.testing.expectEqual(@as(f64, 5.0), result.carbonate_exchange.hco3_mol);
}

test "REDIST soil salt macropore dynamic runtime layers remain independent" {
    var pools = [_]SaltMacroporePools{
        std.mem.zeroes(SaltMacroporePools),
        std.mem.zeroes(SaltMacroporePools),
    };
    var exchange = [_]CarbonateExchangeSite{ .{ .hco3_mol = 0 }, .{ .hco3_mol = 0 } };
    var fluxes = [_]SaltMacroporeFluxes{
        std.mem.zeroes(SaltMacroporeFluxes),
        std.mem.zeroes(SaltMacroporeFluxes),
    };
    fluxes[0].al_macropore_mol = 3.0;
    fluxes[0].al_exchange_mol = 1.0;
    fluxes[1].ca2pb_macropore_mol = 5.0;
    const exchange_fluxes = [_]CarbonateExchangeFlux{ .{ .hco3_mol = 7 }, .{ .hco3_mol = 11 } };
    try updateLayers(&pools, &exchange, &fluxes, &exchange_fluxes, .dynamic);
    try std.testing.expectEqual(@as(f64, 2.0), pools[0].al_mol);
    try std.testing.expectEqual(@as(f64, 5.0), pools[1].ca2pb_mol);
    try std.testing.expectEqual(@as(f64, 7.0), exchange[0].hco3_mol);
    try std.testing.expectEqual(@as(f64, 11.0), exchange[1].hco3_mol);
}

test "REDIST soil salt macropore static gate and failures" {
    var pools = [_]SaltMacroporePools{std.mem.zeroes(SaltMacroporePools)};
    var exchange = [_]CarbonateExchangeSite{.{ .hco3_mol = 2 }};
    var fluxes = [_]SaltMacroporeFluxes{std.mem.zeroes(SaltMacroporeFluxes)};
    fluxes[0].al_macropore_mol = 3;
    const exchange_fluxes = [_]CarbonateExchangeFlux{.{ .hco3_mol = 4 }};
    try updateLayers(&pools, &exchange, &fluxes, &exchange_fluxes, .static);
    try std.testing.expectEqual(@as(f64, 0), pools[0].al_mol);
    try std.testing.expectEqual(@as(f64, 2), exchange[0].hco3_mol);

    const no_fluxes: [0]SaltMacroporeFluxes = .{};
    try std.testing.expectError(
        error.SaltMacroporeDimensionMismatch,
        updateLayers(&pools, &exchange, &no_fluxes, &exchange_fluxes, .dynamic),
    );
    try std.testing.expectError(
        error.NonFiniteCarbonateExchangePool,
        updateLayer(
            pools[0],
            .{ .hco3_mol = std.math.floatMax(f64) },
            fluxes[0],
            .{ .hco3_mol = std.math.floatMax(f64) },
        ),
    );
}
