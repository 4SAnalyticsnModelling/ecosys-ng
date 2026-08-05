const std = @import("std");

/// Exchangeable cation/anion pools at layer NU (mol).
/// 20 fields in REDIST source order (lines 5232--5251).
pub const ExchangePools = struct {
    /// XN4. Adsorbed NH4 non-band.
    nh4_nonband_mol: f64,
    /// XNB. Adsorbed NH4 band.
    nh4_band_mol: f64,
    /// XHY. Adsorbed H.
    h_mol: f64,
    /// XAL. Adsorbed Al.
    al_mol: f64,
    /// XFE. Adsorbed Fe.
    fe_mol: f64,
    /// XCA. Adsorbed Ca.
    ca_mol: f64,
    /// XMG. Adsorbed Mg.
    megagrams_mol: f64,
    /// XNA. Adsorbed Na.
    na_mol: f64,
    /// XKA. Adsorbed K.
    ka_mol: f64,
    /// XHC. Adsorbed HCO3.
    hco3_mol: f64,
    /// XOH0. Adsorption site R- non-band.
    oh0_nonband_mol: f64,
    /// XOH1. Adsorption site R-OH non-band.
    oh1_nonband_mol: f64,
    /// XOH2. Adsorption site R-OH2 non-band.
    oh2_nonband_mol: f64,
    /// XH1P. Adsorbed HPO4 non-band.
    hpo4_nonband_mol: f64,
    /// XH2P. Adsorbed H2PO4 non-band.
    h2po4_nonband_mol: f64,
    /// XOH0B. Adsorption site R- band.
    oh0_band_mol: f64,
    /// XOH1B. Adsorption site R-OH band.
    oh1_band_mol: f64,
    /// XOH2B. Adsorption site R-OH2 band.
    oh2_band_mol: f64,
    /// XH1PB. Adsorbed HPO4 band.
    hpo4_band_mol: f64,
    /// XH2PB. Adsorbed H2PO4 band.
    h2po4_band_mol: f64,
};

/// Signed exchange erosion increments (mol per model step). Same 20 fields.
pub const ExchangeFluxes = struct {
    nh4_nonband_mol: f64, // TN4ER
    nh4_band_mol: f64, // TNBER
    h_mol: f64, // THYER
    al_mol: f64, // TALER
    fe_mol: f64, // TFEER
    ca_mol: f64, // TCAER
    megagrams_mol: f64, // TMGER
    na_mol: f64, // TNAER
    ka_mol: f64, // TKAER
    hco3_mol: f64, // THCER
    oh0_nonband_mol: f64, // TOH0ER
    oh1_nonband_mol: f64, // TOH1ER
    oh2_nonband_mol: f64, // TOH2ER
    hpo4_nonband_mol: f64, // TH1PER
    h2po4_nonband_mol: f64, // TH2PER
    oh0_band_mol: f64, // TOH0EB
    oh1_band_mol: f64, // TOH1EB
    oh2_band_mol: f64, // TOH2EB
    hpo4_band_mol: f64, // TH1PEB
    h2po4_band_mol: f64, // TH2PEB
};

/// Precipitate and silicate pools at layer NU (mol).
/// 26 fields in REDIST source order (lines 5263--5288).
pub const PrecipitatePools = struct {
    // Hydroxide precipitates
    /// PALOH. Al hydroxide precipitate.
    aloh_mol: f64,
    /// PFEOH. Fe hydroxide precipitate.
    feoh_mol: f64,
    // Salt precipitates
    /// PCACO. CaCO3 precipitate.
    caco3_mol: f64,
    /// PCASO. CaSO4 precipitate.
    caso4_mol: f64,
    // Primary silicates
    /// QALSI. Al primary silicate.
    alsi_mol: f64,
    /// QFESI. Fe primary silicate.
    fesi_mol: f64,
    /// QCASI. Ca primary silicate.
    casi_mol: f64,
    /// QMGSI. Mg primary silicate.
    mgsi_mol: f64,
    /// QNASI. Na primary silicate.
    nasi_mol: f64,
    /// QKASI. K primary silicate.
    kasi_mol: f64,
    // Secondary silicates
    /// QALSIF. Al secondary silicate.
    alsif_mol: f64,
    /// QFESIF. Fe secondary silicate.
    fesif_mol: f64,
    /// QCASIF. Ca secondary silicate.
    casif_mol: f64,
    /// QMGSIF. Mg secondary silicate.
    mgsif_mol: f64,
    /// QNASIF. Na secondary silicate.
    nasif_mol: f64,
    /// QKASIF. K secondary silicate.
    kasif_mol: f64,
    // Phosphate precipitates non-band
    /// PALPO. AlPO4 non-band.
    alpo4_nonband_mol: f64,
    /// PFEPO. FePO4 non-band.
    fepo4_nonband_mol: f64,
    /// PCAPD. CaHPO4 non-band.
    cahpo4_nonband_mol: f64,
    /// PCAPH. Apatite non-band.
    apatite_nonband_mol: f64,
    /// PCAPM. CaH2PO4 non-band.
    cah2po4_nonband_mol: f64,
    // Phosphate precipitates band
    /// PALPB. AlPO4 band.
    alpo4_band_mol: f64,
    /// PFEPB. FePO4 band.
    fepo4_band_mol: f64,
    /// PCPDB. CaHPO4 band.
    cahpo4_band_mol: f64,
    /// PCPHB. Apatite band.
    apatite_band_mol: f64,
    /// PCPMB. CaH2PO4 band.
    cah2po4_band_mol: f64,
};

/// Signed precipitate erosion increments (mol per model step). Same 26 fields.
pub const PrecipitateFluxes = struct {
    aloh_mol: f64, // TALOER
    feoh_mol: f64, // TFEOER
    caco3_mol: f64, // TCACER
    caso4_mol: f64, // TCASER
    alsi_mol: f64, // TQALER
    fesi_mol: f64, // TQFEER
    casi_mol: f64, // TQCAER
    mgsi_mol: f64, // TQMGER
    nasi_mol: f64, // TQNAER
    kasi_mol: f64, // TQKAER
    alsif_mol: f64, // TQALERF
    fesif_mol: f64, // TQFEERF
    casif_mol: f64, // TQCAERF
    mgsif_mol: f64, // TQMGERF
    nasif_mol: f64, // TQNAERF
    kasif_mol: f64, // TQKAERF
    alpo4_nonband_mol: f64, // TALPER
    fepo4_nonband_mol: f64, // TFEPER
    cahpo4_nonband_mol: f64, // TCPDER
    apatite_nonband_mol: f64, // TCPHER
    cah2po4_nonband_mol: f64, // TCPMER
    alpo4_band_mol: f64, // TALPEB
    fepo4_band_mol: f64, // TFEPEB
    cahpo4_band_mol: f64, // TCPDEB
    apatite_band_mol: f64, // TCPHEB
    cah2po4_band_mol: f64, // TCPMEB
};

pub const Result = struct {
    exchange: ExchangePools,
    precipitate: PrecipitatePools,
};

/// Direct translation of REDIST lines 5232--5288 (IERSNG-gated block).
///
/// Caller must preserve the enclosing REDIST gate: IERSNG is 1 or 3 and
/// ABS(TSEDER) > ZEROS. Mode 2 does not execute this erosion block.
pub fn apply(
    exchange: ExchangePools,
    precipitate: PrecipitatePools,
    ef: ExchangeFluxes,
    pf: PrecipitateFluxes,
) !Result {
    inline for (@typeInfo(ExchangePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(exchange, field.name)))
            return error.InvalidExchangePool;
    inline for (@typeInfo(PrecipitatePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(precipitate, field.name)))
            return error.InvalidPrecipitatePool;
    inline for (@typeInfo(ExchangeFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(ef, field.name)))
            return error.InvalidExchangeFlux;
    inline for (@typeInfo(PrecipitateFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(pf, field.name)))
            return error.InvalidPrecipitateFlux;

    const new_exchange = ExchangePools{
        .nh4_nonband_mol = exchange.nh4_nonband_mol + ef.nh4_nonband_mol,
        .nh4_band_mol = exchange.nh4_band_mol + ef.nh4_band_mol,
        .h_mol = exchange.h_mol + ef.h_mol,
        .al_mol = exchange.al_mol + ef.al_mol,
        .fe_mol = exchange.fe_mol + ef.fe_mol,
        .ca_mol = exchange.ca_mol + ef.ca_mol,
        .megagrams_mol = exchange.megagrams_mol + ef.megagrams_mol,
        .na_mol = exchange.na_mol + ef.na_mol,
        .ka_mol = exchange.ka_mol + ef.ka_mol,
        .hco3_mol = exchange.hco3_mol + ef.hco3_mol,
        .oh0_nonband_mol = exchange.oh0_nonband_mol + ef.oh0_nonband_mol,
        .oh1_nonband_mol = exchange.oh1_nonband_mol + ef.oh1_nonband_mol,
        .oh2_nonband_mol = exchange.oh2_nonband_mol + ef.oh2_nonband_mol,
        .hpo4_nonband_mol = exchange.hpo4_nonband_mol + ef.hpo4_nonband_mol,
        .h2po4_nonband_mol = exchange.h2po4_nonband_mol + ef.h2po4_nonband_mol,
        .oh0_band_mol = exchange.oh0_band_mol + ef.oh0_band_mol,
        .oh1_band_mol = exchange.oh1_band_mol + ef.oh1_band_mol,
        .oh2_band_mol = exchange.oh2_band_mol + ef.oh2_band_mol,
        .hpo4_band_mol = exchange.hpo4_band_mol + ef.hpo4_band_mol,
        .h2po4_band_mol = exchange.h2po4_band_mol + ef.h2po4_band_mol,
    };
    const new_precipitate = PrecipitatePools{
        .aloh_mol = precipitate.aloh_mol + pf.aloh_mol,
        .feoh_mol = precipitate.feoh_mol + pf.feoh_mol,
        .caco3_mol = precipitate.caco3_mol + pf.caco3_mol,
        .caso4_mol = precipitate.caso4_mol + pf.caso4_mol,
        .alsi_mol = precipitate.alsi_mol + pf.alsi_mol,
        .fesi_mol = precipitate.fesi_mol + pf.fesi_mol,
        .casi_mol = precipitate.casi_mol + pf.casi_mol,
        .mgsi_mol = precipitate.mgsi_mol + pf.mgsi_mol,
        .nasi_mol = precipitate.nasi_mol + pf.nasi_mol,
        .kasi_mol = precipitate.kasi_mol + pf.kasi_mol,
        .alsif_mol = precipitate.alsif_mol + pf.alsif_mol,
        .fesif_mol = precipitate.fesif_mol + pf.fesif_mol,
        .casif_mol = precipitate.casif_mol + pf.casif_mol,
        .mgsif_mol = precipitate.mgsif_mol + pf.mgsif_mol,
        .nasif_mol = precipitate.nasif_mol + pf.nasif_mol,
        .kasif_mol = precipitate.kasif_mol + pf.kasif_mol,
        .alpo4_nonband_mol = precipitate.alpo4_nonband_mol + pf.alpo4_nonband_mol,
        .fepo4_nonband_mol = precipitate.fepo4_nonband_mol + pf.fepo4_nonband_mol,
        .cahpo4_nonband_mol = precipitate.cahpo4_nonband_mol + pf.cahpo4_nonband_mol,
        .apatite_nonband_mol = precipitate.apatite_nonband_mol + pf.apatite_nonband_mol,
        .cah2po4_nonband_mol = precipitate.cah2po4_nonband_mol + pf.cah2po4_nonband_mol,
        .alpo4_band_mol = precipitate.alpo4_band_mol + pf.alpo4_band_mol,
        .fepo4_band_mol = precipitate.fepo4_band_mol + pf.fepo4_band_mol,
        .cahpo4_band_mol = precipitate.cahpo4_band_mol + pf.cahpo4_band_mol,
        .apatite_band_mol = precipitate.apatite_band_mol + pf.apatite_band_mol,
        .cah2po4_band_mol = precipitate.cah2po4_band_mol + pf.cah2po4_band_mol,
    };

    const result = Result{ .exchange = new_exchange, .precipitate = new_precipitate };
    inline for (@typeInfo(ExchangePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.exchange, field.name)))
            return error.NonFiniteExchangePool;
    inline for (@typeInfo(PrecipitatePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.precipitate, field.name)))
            return error.NonFinitePrecipitatePool;
    return result;
}

test "REDIST erosion exchange NH4 non-band and band accumulate independently" {
    var ef = std.mem.zeroes(ExchangeFluxes);
    ef.nh4_nonband_mol = 1.0;
    ef.nh4_band_mol = 2.0;
    const result = try apply(
        std.mem.zeroes(ExchangePools),
        std.mem.zeroes(PrecipitatePools),
        ef,
        std.mem.zeroes(PrecipitateFluxes),
    );
    try std.testing.expectApproxEqRel(@as(f64, 1.0), result.exchange.nh4_nonband_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.exchange.nh4_band_mol, 1.0e-15);
}

test "REDIST erosion precipitate hydroxides and silicates are independent pools" {
    var pf = std.mem.zeroes(PrecipitateFluxes);
    pf.aloh_mol = 0.3;
    pf.alsi_mol = 0.7;
    const result = try apply(
        std.mem.zeroes(ExchangePools),
        std.mem.zeroes(PrecipitatePools),
        std.mem.zeroes(ExchangeFluxes),
        pf,
    );
    try std.testing.expectApproxEqRel(@as(f64, 0.3), result.precipitate.aloh_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 0.7), result.precipitate.alsi_mol, 1.0e-15);
}

test "REDIST erosion precipitate band phosphate pools" {
    var pf = std.mem.zeroes(PrecipitateFluxes);
    pf.alpo4_band_mol = 0.5;
    pf.apatite_band_mol = 0.25;
    const result = try apply(
        std.mem.zeroes(ExchangePools),
        std.mem.zeroes(PrecipitatePools),
        std.mem.zeroes(ExchangeFluxes),
        pf,
    );
    try std.testing.expectApproxEqRel(@as(f64, 0.5), result.precipitate.alpo4_band_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 0.25), result.precipitate.apatite_band_mol, 1.0e-15);
}

test "REDIST erosion exchange rejects non-finite flux" {
    var ef = std.mem.zeroes(ExchangeFluxes);
    ef.ca_mol = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidExchangeFlux,
        apply(
            std.mem.zeroes(ExchangePools),
            std.mem.zeroes(PrecipitatePools),
            ef,
            std.mem.zeroes(PrecipitateFluxes),
        ),
    );
}

test "REDIST erosion exchange and precipitate map all 46 source fields" {
    var ef = std.mem.zeroes(ExchangeFluxes);
    var exchange_value: f64 = 1.0;
    inline for (@typeInfo(ExchangeFluxes).@"struct".fields) |field| {
        @field(ef, field.name) = exchange_value;
        exchange_value += 1.0;
    }
    var pf = std.mem.zeroes(PrecipitateFluxes);
    var precipitate_value: f64 = 31.0;
    inline for (@typeInfo(PrecipitateFluxes).@"struct".fields) |field| {
        @field(pf, field.name) = precipitate_value;
        precipitate_value += 1.0;
    }
    const result = try apply(
        std.mem.zeroes(ExchangePools),
        std.mem.zeroes(PrecipitatePools),
        ef,
        pf,
    );
    inline for (@typeInfo(ExchangePools).@"struct".fields) |field|
        try std.testing.expectEqual(@field(ef, field.name), @field(result.exchange, field.name));
    inline for (@typeInfo(PrecipitatePools).@"struct".fields) |field|
        try std.testing.expectEqual(@field(pf, field.name), @field(result.precipitate, field.name));
}

test "REDIST erosion exchange and precipitate preserve signed removal" {
    var exchange = std.mem.zeroes(ExchangePools);
    exchange.hpo4_band_mol = 4.0;
    var ef = std.mem.zeroes(ExchangeFluxes);
    ef.hpo4_band_mol = -1.5;
    var precipitate = std.mem.zeroes(PrecipitatePools);
    precipitate.kasif_mol = 3.0;
    var pf = std.mem.zeroes(PrecipitateFluxes);
    pf.kasif_mol = -0.75;
    const result = try apply(exchange, precipitate, ef, pf);
    try std.testing.expectEqual(@as(f64, 2.5), result.exchange.hpo4_band_mol);
    try std.testing.expectEqual(@as(f64, 2.25), result.precipitate.kasif_mol);
}

test "REDIST erosion precipitate rejects non-finite flux" {
    var pf = std.mem.zeroes(PrecipitateFluxes);
    pf.cahpo4_band_mol = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidPrecipitateFlux,
        apply(
            std.mem.zeroes(ExchangePools),
            std.mem.zeroes(PrecipitatePools),
            std.mem.zeroes(ExchangeFluxes),
            pf,
        ),
    );
}

test "REDIST erosion exchange and precipitate reject arithmetic overflow" {
    var exchange = std.mem.zeroes(ExchangePools);
    exchange.oh1_band_mol = std.math.floatMax(f64);
    var ef = std.mem.zeroes(ExchangeFluxes);
    ef.oh1_band_mol = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteExchangePool,
        apply(exchange, std.mem.zeroes(PrecipitatePools), ef, std.mem.zeroes(PrecipitateFluxes)),
    );

    var precipitate = std.mem.zeroes(PrecipitatePools);
    precipitate.apatite_nonband_mol = std.math.floatMax(f64);
    var pf = std.mem.zeroes(PrecipitateFluxes);
    pf.apatite_nonband_mol = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFinitePrecipitatePool,
        apply(std.mem.zeroes(ExchangePools), precipitate, std.mem.zeroes(ExchangeFluxes), pf),
    );
}
