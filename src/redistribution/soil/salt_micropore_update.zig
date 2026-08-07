const std = @import("std");

/// Soil salt solution micropore pools (mol). 50 fields, REDIST 6956-7062.
/// Gated on ISALTG != 0 by the caller.
pub const SaltMicroporePools = struct {
    // Base cations (lines 6956-6977)
    h_mol: f64, // ZHY
    oh_mol: f64, // ZOH
    al_mol: f64, // ZAL
    fe_mol: f64, // ZFE
    ca_mol: f64, // ZCA
    megagrams_mol: f64, // ZMG
    na_mol: f64, // ZNA
    ka_mol: f64, // ZKA
    // Anions (lines 6978-6987)
    so4_mol: f64, // ZSO4
    cl_mol: f64, // ZCL
    co3_mol: f64, // ZCO3
    hco3_mol: f64, // ZHCO3
    // Al complexes (lines 6988-6997)
    aloh1_mol: f64, // ZALOH1
    aloh2_mol: f64, // ZALOH2
    aloh3_mol: f64, // ZALOH3
    aloh4_mol: f64, // ZALOH4
    als_mol: f64, // ZALS
    // Fe complexes (lines 6998-7007)
    feoh1_mol: f64, // ZFEOH1
    feoh2_mol: f64, // ZFEOH2
    feoh3_mol: f64, // ZFEOH3
    feoh4_mol: f64, // ZFEOH4
    fes_mol: f64, // ZFES
    // Ca complexes (lines 7008-7015)
    cao_mol: f64, // ZCAO
    cac_mol: f64, // ZCAC
    cah_mol: f64, // ZCAH
    cas_mol: f64, // ZCAS
    // Mg complexes (lines 7016-7023)
    mgo_mol: f64, // ZMGO
    mgc_mol: f64, // ZMGC
    mgh_mol: f64, // ZMGH
    mgs_mol: f64, // ZMGS
    // Na/K complexes (lines 7024-7029)
    nac_mol: f64, // ZNAC
    nas_mol: f64, // ZNAS
    kas_mol: f64, // ZKAS
    // Silicate (line 7030)
    hysi_mol: f64, // ZHYSI
    // Phosphate non-band (lines 7031-7046)
    h0po4_mol: f64, // H0PO4
    h3po4_mol: f64, // H3PO4
    fe1p_mol: f64, // ZFE1P
    fe2p_mol: f64, // ZFE2P
    ca0p_mol: f64, // ZCA0P
    ca1p_mol: f64, // ZCA1P
    ca2p_mol: f64, // ZCA2P
    mg1p_mol: f64, // ZMG1P
    // Phosphate band (lines 7047-7062)
    h0pob_mol: f64, // H0POB
    h3pob_mol: f64, // H3POB
    fe1pb_mol: f64, // ZFE1PB
    fe2pb_mol: f64, // ZFE2PB
    ca0pb_mol: f64, // ZCA0PB
    ca1pb_mol: f64, // ZCA1PB
    ca2pb_mol: f64, // ZCA2PB
    mg1pb_mol: f64, // ZMG1PB
};

/// Per-species flux bundle for salt micropore transport. Sizes vary:
///  Base cations/anions: transport + subsurface + pore_exchange + [root_uptake] + [senescence]
///  Complexes: transformation + transport + pore_exchange (no root/senes)
///  Special: ZALOH2 and ZFEOH2 have an additional -desorption term (TRXAL2/TRXFE2).
pub const SaltMicroporeFluxes = struct {
    // H (line 6956)
    h_transport_mol: f64, // THYFLS
    h_subsurface_mol: f64, // RHYFLU
    h_exchange_mol: f64, // XHYFXS
    // OH (line 6958)
    oh_transport_mol: f64,
    oh_subsurface_mol: f64,
    oh_exchange_mol: f64,
    // Al (line 6960)
    al_transport_mol: f64,
    al_subsurface_mol: f64,
    al_exchange_mol: f64,
    al_root_uptake_mol: f64, // TUPZAL (subtracted)
    al_senescence_mol: f64, // ALSNT
    // Fe (line 6963)
    fe_transport_mol: f64,
    fe_subsurface_mol: f64,
    fe_exchange_mol: f64,
    fe_root_uptake_mol: f64,
    fe_senescence_mol: f64, // FESNT
    // Ca (line 6966)
    ca_transport_mol: f64,
    ca_subsurface_mol: f64,
    ca_exchange_mol: f64,
    ca_root_uptake_mol: f64,
    ca_senescence_mol: f64, // CASNT
    // Mg (line 6969)
    megagrams_transport_mol: f64,
    megagrams_subsurface_mol: f64,
    megagrams_exchange_mol: f64,
    megagrams_root_uptake_mol: f64,
    megagrams_senescence_mol: f64, // GMSNT
    // Na (line 6972)
    na_transport_mol: f64,
    na_subsurface_mol: f64,
    na_exchange_mol: f64,
    na_root_uptake_mol: f64,
    na_senescence_mol: f64, // ANSNT
    // Ka (line 6975)
    ka_transport_mol: f64,
    ka_subsurface_mol: f64,
    ka_exchange_mol: f64,
    ka_root_uptake_mol: f64,
    ka_senescence_mol: f64, // AKSNT
    // SO4 (line 6978)
    so4_transform_mol: f64, // TRSO4
    so4_transport_mol: f64,
    so4_subsurface_mol: f64,
    so4_exchange_mol: f64,
    so4_root_uptake_mol: f64,
    so4_senescence_mol: f64, // SOSNT
    // Cl (line 6981)
    cl_transport_mol: f64,
    cl_subsurface_mol: f64,
    cl_exchange_mol: f64,
    cl_root_uptake_mol: f64,
    cl_senescence_mol: f64, // CLSNT
    // CO3 (line 6984)
    co3_transform_mol: f64, // TRCO3
    co3_transport_mol: f64,
    co3_exchange_mol: f64,
    // HCO3 (line 6986)
    hco3_transform_mol: f64, // TRHCO
    hco3_transport_mol: f64,
    hco3_exchange_mol: f64,
    // ALOH1 (line 6988)
    aloh1_transform_mol: f64,
    aloh1_transport_mol: f64,
    aloh1_exchange_mol: f64,
    // ALOH2 (line 6990): also has -TRXAL2 (desorption)
    aloh2_transform_mol: f64,
    aloh2_transport_mol: f64,
    aloh2_exchange_mol: f64,
    aloh2_desorption_mol: f64, // TRXAL2 (subtracted)
    // ALOH3 (line 6992)
    aloh3_transform_mol: f64,
    aloh3_transport_mol: f64,
    aloh3_exchange_mol: f64,
    // ALOH4 (line 6994)
    aloh4_transform_mol: f64,
    aloh4_transport_mol: f64,
    aloh4_exchange_mol: f64,
    // ALS (line 6996)
    als_transform_mol: f64,
    als_transport_mol: f64,
    als_exchange_mol: f64,
    // FEOH1 (line 6998)
    feoh1_transform_mol: f64,
    feoh1_transport_mol: f64,
    feoh1_exchange_mol: f64,
    // FEOH2 (line 7000): also has -TRXFE2 (desorption)
    feoh2_transform_mol: f64,
    feoh2_transport_mol: f64,
    feoh2_exchange_mol: f64,
    feoh2_desorption_mol: f64, // TRXFE2 (subtracted)
    // FEOH3 (line 7002)
    feoh3_transform_mol: f64,
    feoh3_transport_mol: f64,
    feoh3_exchange_mol: f64,
    // FEOH4 (line 7004)
    feoh4_transform_mol: f64,
    feoh4_transport_mol: f64,
    feoh4_exchange_mol: f64,
    // FES (line 7006)
    fes_transform_mol: f64,
    fes_transport_mol: f64,
    fes_exchange_mol: f64,
    // CAO (line 7008)
    cao_transform_mol: f64,
    cao_transport_mol: f64,
    cao_exchange_mol: f64,
    // CAC (line 7010)
    cac_transform_mol: f64,
    cac_transport_mol: f64,
    cac_exchange_mol: f64,
    // CAH (line 7012)
    cah_transform_mol: f64,
    cah_transport_mol: f64,
    cah_exchange_mol: f64,
    // CAS (line 7014)
    cas_transform_mol: f64,
    cas_transport_mol: f64,
    cas_exchange_mol: f64,
    // MGO (line 7016)
    mgo_transform_mol: f64,
    mgo_transport_mol: f64,
    mgo_exchange_mol: f64,
    // MGC (line 7018)
    mgc_transform_mol: f64,
    mgc_transport_mol: f64,
    mgc_exchange_mol: f64,
    // MGH (line 7020)
    mgh_transform_mol: f64,
    mgh_transport_mol: f64,
    mgh_exchange_mol: f64,
    // MGS (line 7022)
    mgs_transform_mol: f64,
    mgs_transport_mol: f64,
    mgs_exchange_mol: f64,
    // NAC (line 7024)
    nac_transform_mol: f64,
    nac_transport_mol: f64,
    nac_exchange_mol: f64,
    // NAS (line 7026)
    nas_transform_mol: f64,
    nas_transport_mol: f64,
    nas_exchange_mol: f64,
    // KAS (line 7028)
    kas_transform_mol: f64,
    kas_transport_mol: f64,
    kas_exchange_mol: f64,
    // HYSI (line 7030)
    hysi_transform_mol: f64, // TRHYSI
    hysi_transport_mol: f64, // THYSIS (no exchange term)
    // H0PO4 (line 7031)
    h0po4_transform_mol: f64,
    h0po4_transport_mol: f64,
    h0po4_exchange_mol: f64,
    // H3PO4 (line 7033)
    h3po4_transform_mol: f64,
    h3po4_transport_mol: f64,
    h3po4_exchange_mol: f64,
    // ZFE1P (line 7035)
    fe1p_transform_mol: f64,
    fe1p_transport_mol: f64,
    fe1p_exchange_mol: f64,
    // ZFE2P (line 7037)
    fe2p_transform_mol: f64,
    fe2p_transport_mol: f64,
    fe2p_exchange_mol: f64,
    // ZCA0P (line 7039)
    ca0p_transform_mol: f64,
    ca0p_transport_mol: f64,
    ca0p_exchange_mol: f64,
    // ZCA1P (line 7041)
    ca1p_transform_mol: f64,
    ca1p_transport_mol: f64,
    ca1p_exchange_mol: f64,
    // ZCA2P (line 7043)
    ca2p_transform_mol: f64,
    ca2p_transport_mol: f64,
    ca2p_exchange_mol: f64,
    // ZMG1P (line 7045)
    mg1p_transform_mol: f64,
    mg1p_transport_mol: f64,
    mg1p_exchange_mol: f64,
    // Band phosphate: T*FLB + X*FXB only (lines 7047-7062)
    h0pob_transform_mol: f64,
    h0pob_transport_mol: f64,
    h0pob_exchange_mol: f64,
    h3pob_transform_mol: f64,
    h3pob_transport_mol: f64,
    h3pob_exchange_mol: f64,
    fe1pb_transform_mol: f64,
    fe1pb_transport_mol: f64,
    fe1pb_exchange_mol: f64,
    fe2pb_transform_mol: f64,
    fe2pb_transport_mol: f64,
    fe2pb_exchange_mol: f64,
    ca0pb_transform_mol: f64,
    ca0pb_transport_mol: f64,
    ca0pb_exchange_mol: f64,
    ca1pb_transform_mol: f64,
    ca1pb_transport_mol: f64,
    ca1pb_exchange_mol: f64,
    ca2pb_transform_mol: f64,
    ca2pb_transport_mol: f64,
    ca2pb_exchange_mol: f64,
    mg1pb_transform_mol: f64,
    mg1pb_transport_mol: f64,
    mg1pb_exchange_mol: f64,
};

pub const Result = struct {
    pools: SaltMicroporePools,
    /// Increment to UCO2S = 12*(ZCO3 + ZHCO3), mol→g C (line 7063).
    uco2s_increment_gC: f64,
};

pub const SaltEquilibriumMode = enum {
    static,
    dynamic,
};

/// Direct translation of redist.f lines 6956--7063.
/// Caller must gate on ISALTG != 0.
pub fn update(pools: SaltMicroporePools, f: SaltMicroporeFluxes) !Result {
    inline for (@typeInfo(SaltMicroporePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(pools, field.name)))
            return error.InvalidSaltMicroporePool;
    inline for (@typeInfo(SaltMicroporeFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(f, field.name)))
            return error.InvalidSaltMicroporeFlux;

    const r = SaltMicroporePools{
        .h_mol = pools.h_mol + f.h_transport_mol + f.h_subsurface_mol + f.h_exchange_mol,
        .oh_mol = pools.oh_mol + f.oh_transport_mol + f.oh_subsurface_mol + f.oh_exchange_mol,
        .al_mol = pools.al_mol + f.al_transport_mol + f.al_subsurface_mol + f.al_exchange_mol - f.al_root_uptake_mol + f.al_senescence_mol,
        .fe_mol = pools.fe_mol + f.fe_transport_mol + f.fe_subsurface_mol + f.fe_exchange_mol - f.fe_root_uptake_mol + f.fe_senescence_mol,
        .ca_mol = pools.ca_mol + f.ca_transport_mol + f.ca_subsurface_mol + f.ca_exchange_mol - f.ca_root_uptake_mol + f.ca_senescence_mol,
        .megagrams_mol = pools.megagrams_mol + f.megagrams_transport_mol + f.megagrams_subsurface_mol + f.megagrams_exchange_mol - f.megagrams_root_uptake_mol + f.megagrams_senescence_mol,
        .na_mol = pools.na_mol + f.na_transport_mol + f.na_subsurface_mol + f.na_exchange_mol - f.na_root_uptake_mol + f.na_senescence_mol,
        .ka_mol = pools.ka_mol + f.ka_transport_mol + f.ka_subsurface_mol + f.ka_exchange_mol - f.ka_root_uptake_mol + f.ka_senescence_mol,
        .so4_mol = pools.so4_mol + f.so4_transform_mol + f.so4_transport_mol + f.so4_subsurface_mol + f.so4_exchange_mol - f.so4_root_uptake_mol + f.so4_senescence_mol,
        .cl_mol = pools.cl_mol + f.cl_transport_mol + f.cl_subsurface_mol + f.cl_exchange_mol - f.cl_root_uptake_mol + f.cl_senescence_mol,
        .co3_mol = pools.co3_mol + f.co3_transform_mol + f.co3_transport_mol + f.co3_exchange_mol,
        .hco3_mol = pools.hco3_mol + f.hco3_transform_mol + f.hco3_transport_mol + f.hco3_exchange_mol,
        .aloh1_mol = pools.aloh1_mol + f.aloh1_transform_mol + f.aloh1_transport_mol + f.aloh1_exchange_mol,
        .aloh2_mol = pools.aloh2_mol + f.aloh2_transform_mol + f.aloh2_transport_mol + f.aloh2_exchange_mol - f.aloh2_desorption_mol,
        .aloh3_mol = pools.aloh3_mol + f.aloh3_transform_mol + f.aloh3_transport_mol + f.aloh3_exchange_mol,
        .aloh4_mol = pools.aloh4_mol + f.aloh4_transform_mol + f.aloh4_transport_mol + f.aloh4_exchange_mol,
        .als_mol = pools.als_mol + f.als_transform_mol + f.als_transport_mol + f.als_exchange_mol,
        .feoh1_mol = pools.feoh1_mol + f.feoh1_transform_mol + f.feoh1_transport_mol + f.feoh1_exchange_mol,
        .feoh2_mol = pools.feoh2_mol + f.feoh2_transform_mol + f.feoh2_transport_mol + f.feoh2_exchange_mol - f.feoh2_desorption_mol,
        .feoh3_mol = pools.feoh3_mol + f.feoh3_transform_mol + f.feoh3_transport_mol + f.feoh3_exchange_mol,
        .feoh4_mol = pools.feoh4_mol + f.feoh4_transform_mol + f.feoh4_transport_mol + f.feoh4_exchange_mol,
        .fes_mol = pools.fes_mol + f.fes_transform_mol + f.fes_transport_mol + f.fes_exchange_mol,
        .cao_mol = pools.cao_mol + f.cao_transform_mol + f.cao_transport_mol + f.cao_exchange_mol,
        .cac_mol = pools.cac_mol + f.cac_transform_mol + f.cac_transport_mol + f.cac_exchange_mol,
        .cah_mol = pools.cah_mol + f.cah_transform_mol + f.cah_transport_mol + f.cah_exchange_mol,
        .cas_mol = pools.cas_mol + f.cas_transform_mol + f.cas_transport_mol + f.cas_exchange_mol,
        .mgo_mol = pools.mgo_mol + f.mgo_transform_mol + f.mgo_transport_mol + f.mgo_exchange_mol,
        .mgc_mol = pools.mgc_mol + f.mgc_transform_mol + f.mgc_transport_mol + f.mgc_exchange_mol,
        .mgh_mol = pools.mgh_mol + f.mgh_transform_mol + f.mgh_transport_mol + f.mgh_exchange_mol,
        .mgs_mol = pools.mgs_mol + f.mgs_transform_mol + f.mgs_transport_mol + f.mgs_exchange_mol,
        .nac_mol = pools.nac_mol + f.nac_transform_mol + f.nac_transport_mol + f.nac_exchange_mol,
        .nas_mol = pools.nas_mol + f.nas_transform_mol + f.nas_transport_mol + f.nas_exchange_mol,
        .kas_mol = pools.kas_mol + f.kas_transform_mol + f.kas_transport_mol + f.kas_exchange_mol,
        .hysi_mol = pools.hysi_mol + f.hysi_transform_mol + f.hysi_transport_mol,
        .h0po4_mol = pools.h0po4_mol + f.h0po4_transform_mol + f.h0po4_transport_mol + f.h0po4_exchange_mol,
        .h3po4_mol = pools.h3po4_mol + f.h3po4_transform_mol + f.h3po4_transport_mol + f.h3po4_exchange_mol,
        .fe1p_mol = pools.fe1p_mol + f.fe1p_transform_mol + f.fe1p_transport_mol + f.fe1p_exchange_mol,
        .fe2p_mol = pools.fe2p_mol + f.fe2p_transform_mol + f.fe2p_transport_mol + f.fe2p_exchange_mol,
        .ca0p_mol = pools.ca0p_mol + f.ca0p_transform_mol + f.ca0p_transport_mol + f.ca0p_exchange_mol,
        .ca1p_mol = pools.ca1p_mol + f.ca1p_transform_mol + f.ca1p_transport_mol + f.ca1p_exchange_mol,
        .ca2p_mol = pools.ca2p_mol + f.ca2p_transform_mol + f.ca2p_transport_mol + f.ca2p_exchange_mol,
        .mg1p_mol = pools.mg1p_mol + f.mg1p_transform_mol + f.mg1p_transport_mol + f.mg1p_exchange_mol,
        .h0pob_mol = pools.h0pob_mol + f.h0pob_transform_mol + f.h0pob_transport_mol + f.h0pob_exchange_mol,
        .h3pob_mol = pools.h3pob_mol + f.h3pob_transform_mol + f.h3pob_transport_mol + f.h3pob_exchange_mol,
        .fe1pb_mol = pools.fe1pb_mol + f.fe1pb_transform_mol + f.fe1pb_transport_mol + f.fe1pb_exchange_mol,
        .fe2pb_mol = pools.fe2pb_mol + f.fe2pb_transform_mol + f.fe2pb_transport_mol + f.fe2pb_exchange_mol,
        .ca0pb_mol = pools.ca0pb_mol + f.ca0pb_transform_mol + f.ca0pb_transport_mol + f.ca0pb_exchange_mol,
        .ca1pb_mol = pools.ca1pb_mol + f.ca1pb_transform_mol + f.ca1pb_transport_mol + f.ca1pb_exchange_mol,
        .ca2pb_mol = pools.ca2pb_mol + f.ca2pb_transform_mol + f.ca2pb_transport_mol + f.ca2pb_exchange_mol,
        .mg1pb_mol = pools.mg1pb_mol + f.mg1pb_transform_mol + f.mg1pb_transport_mol + f.mg1pb_exchange_mol,
    };
    inline for (@typeInfo(SaltMicroporePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(r, field.name)))
            return error.NonFiniteSaltMicroporePool;
    const uco2s_increment_gC = 12.0 * (r.co3_mol + r.hco3_mol);
    if (!std.math.isFinite(uco2s_increment_gC))
        return error.NonFiniteSaltMicroporeDiagnostic;
    return Result{
        .pools = r,
        .uco2s_increment_gC = uco2s_increment_gC,
    };
}

/// Applies the complete `ISALTG != 0` micropore block in ascending runtime
/// layer order, including the source-order UCO2S accumulation.
pub fn updateLayers(
    pools_by_layer: []SaltMicroporePools,
    fluxes_by_layer: []const SaltMicroporeFluxes,
    grid_cell_carbon_g_c: *f64,
    mode: SaltEquilibriumMode,
) !void {
    if (pools_by_layer.len == 0 or fluxes_by_layer.len != pools_by_layer.len)
        return error.SaltMicroporeDimensionMismatch;
    if (!std.math.isFinite(grid_cell_carbon_g_c.*))
        return error.InvalidSaltMicroporeDiagnostic;
    if (mode == .static) return;

    var next_carbon_g_c = grid_cell_carbon_g_c.*;
    for (pools_by_layer, fluxes_by_layer) |*pools, fluxes| {
        const result = try update(pools.*, fluxes);
        next_carbon_g_c = next_carbon_g_c + result.uco2s_increment_gC;
        if (!std.math.isFinite(next_carbon_g_c))
            return error.NonFiniteSaltMicroporeDiagnostic;
        pools.* = result.pools;
    }
    grid_cell_carbon_g_c.* = next_carbon_g_c;
}

test "REDIST soil salt micropore base cation transport and subsurface add" {
    var f = std.mem.zeroes(SaltMicroporeFluxes);
    f.ca_transport_mol = 2.0;
    f.ca_subsurface_mol = 1.0;
    const result = try update(std.mem.zeroes(SaltMicroporePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 3.0), result.pools.ca_mol, 1.0e-15);
}

test "REDIST soil salt micropore root uptake subtracts from Al and Fe" {
    var f = std.mem.zeroes(SaltMicroporeFluxes);
    f.al_transport_mol = 5.0;
    f.al_root_uptake_mol = 2.0;
    f.fe_transport_mol = 3.0;
    f.fe_root_uptake_mol = 1.0;
    const result = try update(std.mem.zeroes(SaltMicroporePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 3.0), result.pools.al_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.pools.fe_mol, 1.0e-15);
}

test "REDIST soil salt micropore ALOH2 and FEOH2 subtract desorption" {
    var f = std.mem.zeroes(SaltMicroporeFluxes);
    f.aloh2_transform_mol = 3.0;
    f.aloh2_desorption_mol = 1.0;
    f.feoh2_transform_mol = 4.0;
    f.feoh2_desorption_mol = 2.0;
    const result = try update(std.mem.zeroes(SaltMicroporePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.pools.aloh2_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.pools.feoh2_mol, 1.0e-15);
}

test "REDIST soil salt micropore UCO2S is 12*(CO3+HCO3) of updated pools" {
    var f = std.mem.zeroes(SaltMicroporeFluxes);
    f.co3_transform_mol = 1.0;
    f.hco3_transform_mol = 2.0;
    const result = try update(std.mem.zeroes(SaltMicroporePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 36.0), result.uco2s_increment_gC, 1.0e-15);
}

test "REDIST soil salt micropore SO4 gets transformation and senescence" {
    var f = std.mem.zeroes(SaltMicroporeFluxes);
    f.so4_transform_mol = 1.0;
    f.so4_transport_mol = 2.0;
    f.so4_senescence_mol = 0.5;
    const result = try update(std.mem.zeroes(SaltMicroporePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 3.5), result.pools.so4_mol, 1.0e-15);
}

test "REDIST soil salt micropore rejects non-finite flux" {
    var bad = std.mem.zeroes(SaltMicroporeFluxes);
    bad.ca1p_exchange_mol = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidSaltMicroporeFlux,
        update(std.mem.zeroes(SaltMicroporePools), bad),
    );
}

test "REDIST soil salt micropore dynamic mode traverses runtime layers in order" {
    var pools = [_]SaltMicroporePools{
        std.mem.zeroes(SaltMicroporePools),
        std.mem.zeroes(SaltMicroporePools),
    };
    var fluxes = [_]SaltMicroporeFluxes{
        std.mem.zeroes(SaltMicroporeFluxes),
        std.mem.zeroes(SaltMicroporeFluxes),
    };
    fluxes[0].co3_transform_mol = 1.0;
    fluxes[1].hco3_transform_mol = 2.0;
    fluxes[0].ca_transport_mol = 3.0;
    fluxes[1].megagrams_transport_mol = 4.0;
    var carbon_g_c: f64 = 5.0;
    try updateLayers(&pools, &fluxes, &carbon_g_c, .dynamic);
    try std.testing.expectEqual(@as(f64, 3.0), pools[0].ca_mol);
    try std.testing.expectEqual(@as(f64, 4.0), pools[1].megagrams_mol);
    try std.testing.expectEqual(@as(f64, 41.0), carbon_g_c);
}

test "REDIST soil salt micropore static mode preserves pools and diagnostic" {
    var pools = [_]SaltMicroporePools{std.mem.zeroes(SaltMicroporePools)};
    var fluxes = [_]SaltMicroporeFluxes{std.mem.zeroes(SaltMicroporeFluxes)};
    fluxes[0].ca_transport_mol = 3.0;
    var carbon_g_c: f64 = 5.0;
    try updateLayers(&pools, &fluxes, &carbon_g_c, .static);
    try std.testing.expectEqual(@as(f64, 0.0), pools[0].ca_mol);
    try std.testing.expectEqual(@as(f64, 5.0), carbon_g_c);
}

test "REDIST soil salt micropore runtime traversal rejects dimensions and diagnostic overflow" {
    var pools = [_]SaltMicroporePools{std.mem.zeroes(SaltMicroporePools)};
    const no_fluxes: [0]SaltMicroporeFluxes = .{};
    var carbon_g_c: f64 = 0.0;
    try std.testing.expectError(
        error.SaltMicroporeDimensionMismatch,
        updateLayers(&pools, &no_fluxes, &carbon_g_c, .dynamic),
    );

    const fluxes = [_]SaltMicroporeFluxes{std.mem.zeroes(SaltMicroporeFluxes)};
    pools[0].co3_mol = std.math.floatMax(f64);
    pools[0].hco3_mol = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSaltMicroporeDiagnostic,
        updateLayers(&pools, &fluxes, &carbon_g_c, .dynamic),
    );
}
