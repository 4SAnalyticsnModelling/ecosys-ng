const std = @import("std");

/// Litter-layer ion exchange adsorption amounts at layer 0 (mol).
pub const ExchangePools = struct {
    /// XN4(0). Exchangeable NH4 (mol N).
    nh4_mol: f64,
    /// XHY(0). Exchangeable H+ (mol).
    h_mol: f64,
    /// XAL(0). Exchangeable Al3+ (mol).
    al_mol: f64,
    /// XFE(0). Exchangeable Fe3+ (mol).
    fe_mol: f64,
    /// XCA(0). Exchangeable Ca2+ (mol).
    ca_mol: f64,
    /// XMG(0). Exchangeable Mg2+ (mol).
    megagrams_mol: f64,
    /// XNA(0). Exchangeable Na+ (mol).
    na_mol: f64,
    /// XKA(0). Exchangeable K+ (mol).
    ka_mol: f64,
    /// XHC(0). Exchangeable HCO3- (mol).
    hco3_mol: f64,
    /// XOH0(0). Adsorption site R- (mol).
    oh0_mol: f64,
    /// XOH1(0). Adsorption site R-OH (mol).
    oh1_mol: f64,
    /// XOH2(0). Adsorption site R-OH2 (mol).
    oh2_mol: f64,
    /// XH1P(0). Exchangeable HPO4 (mol P).
    hpo4_mol: f64,
    /// XH2P(0). Exchangeable H2PO4 (mol P).
    h2po4_mol: f64,
};

/// Litter-layer precipitate amounts at layer 0 (mol P or mol).
pub const PrecipitatePools = struct {
    /// PALOH(0). AlOH precipitate (mol).
    aloh_mol: f64,
    /// PFEOH(0). FeOH precipitate (mol).
    feoh_mol: f64,
    /// PCACO(0). CaCO3 precipitate (mol).
    caco3_mol: f64,
    /// PCASO(0). CaSO4 precipitate (mol).
    caso4_mol: f64,
    /// PALPO(0). AlPO4 precipitate (mol P).
    alpo4_mol: f64,
    /// PFEPO(0). FePO4 precipitate (mol P).
    fepo4_mol: f64,
    /// PCAPD(0). CaHPO4 precipitate (mol P).
    cahpo4_mol: f64,
    /// PCAPH(0). Apatite precipitate (mol P).
    apatite_mol: f64,
    /// PCAPM(0). CaH2PO4 precipitate (mol P).
    cah2po4_mol: f64,
};

/// Signed exchange transformations from `solute.f` (mol per model step).
pub const ExchangeTransformations = struct {
    /// TRXN4(0).
    nh4_mol: f64,
    /// TRXHY(0).
    h_mol: f64,
    /// TRXAL(0).
    al_mol: f64,
    /// TRXFE(0).
    fe_mol: f64,
    /// TRXCA(0).
    ca_mol: f64,
    /// TRXMG(0).
    megagrams_mol: f64,
    /// TRXNA(0).
    na_mol: f64,
    /// TRXKA(0).
    ka_mol: f64,
    /// TRXHC(0).
    hco3_mol: f64,
    /// TRXH0(0).
    oh0_mol: f64,
    /// TRXH1(0).
    oh1_mol: f64,
    /// TRXH2(0).
    oh2_mol: f64,
    /// TRX1P(0).
    hpo4_mol: f64,
    /// TRX2P(0).
    h2po4_mol: f64,
};

/// Signed precipitation transformations from `solute.f` (mol per model step).
pub const PrecipitateTransformations = struct {
    /// TRALOH(0).
    aloh_mol: f64,
    /// TRFEOH(0).
    feoh_mol: f64,
    /// TRCACO(0).
    caco3_mol: f64,
    /// TRCASO(0).
    caso4_mol: f64,
    /// TRALPO(0).
    alpo4_mol: f64,
    /// TRFEPO(0).
    fepo4_mol: f64,
    /// TRCAPD(0).
    cahpo4_mol: f64,
    /// TRCAPH(0).
    apatite_mol: f64,
    /// TRCAPM(0).
    cah2po4_mol: f64,
};

pub const Result = struct {
    exchange: ExchangePools,
    precipitate: PrecipitatePools,
};

/// Direct translation of REDIST lines 4978--5000.
///
/// Applies solute.f exchange and precipitation transformations to litter ion
/// pools.
pub fn apply(
    exchange: ExchangePools,
    precipitate: PrecipitatePools,
    ex_tr: ExchangeTransformations,
    pr_tr: PrecipitateTransformations,
) !Result {
    inline for (@typeInfo(ExchangePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(exchange, field.name)))
            return error.InvalidLitterIonExchangeInput;
    inline for (@typeInfo(PrecipitatePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(precipitate, field.name)))
            return error.InvalidLitterIonExchangeInput;
    inline for (@typeInfo(ExchangeTransformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(ex_tr, field.name)))
            return error.InvalidLitterIonExchangeInput;
    inline for (@typeInfo(PrecipitateTransformations).@"struct".fields) |field|
        if (!std.math.isFinite(@field(pr_tr, field.name)))
            return error.InvalidLitterIonExchangeInput;

    const new_ex = ExchangePools{
        .nh4_mol = exchange.nh4_mol + ex_tr.nh4_mol,
        .h_mol = exchange.h_mol + ex_tr.h_mol,
        .al_mol = exchange.al_mol + ex_tr.al_mol,
        .fe_mol = exchange.fe_mol + ex_tr.fe_mol,
        .ca_mol = exchange.ca_mol + ex_tr.ca_mol,
        .megagrams_mol = exchange.megagrams_mol + ex_tr.megagrams_mol,
        .na_mol = exchange.na_mol + ex_tr.na_mol,
        .ka_mol = exchange.ka_mol + ex_tr.ka_mol,
        .hco3_mol = exchange.hco3_mol + ex_tr.hco3_mol,
        .oh0_mol = exchange.oh0_mol + ex_tr.oh0_mol,
        .oh1_mol = exchange.oh1_mol + ex_tr.oh1_mol,
        .oh2_mol = exchange.oh2_mol + ex_tr.oh2_mol,
        .hpo4_mol = exchange.hpo4_mol + ex_tr.hpo4_mol,
        .h2po4_mol = exchange.h2po4_mol + ex_tr.h2po4_mol,
    };
    const new_pr = PrecipitatePools{
        .aloh_mol = precipitate.aloh_mol + pr_tr.aloh_mol,
        .feoh_mol = precipitate.feoh_mol + pr_tr.feoh_mol,
        .caco3_mol = precipitate.caco3_mol + pr_tr.caco3_mol,
        .caso4_mol = precipitate.caso4_mol + pr_tr.caso4_mol,
        .alpo4_mol = precipitate.alpo4_mol + pr_tr.alpo4_mol,
        .fepo4_mol = precipitate.fepo4_mol + pr_tr.fepo4_mol,
        .cahpo4_mol = precipitate.cahpo4_mol + pr_tr.cahpo4_mol,
        .apatite_mol = precipitate.apatite_mol + pr_tr.apatite_mol,
        .cah2po4_mol = precipitate.cah2po4_mol + pr_tr.cah2po4_mol,
    };

    const result = Result{ .exchange = new_ex, .precipitate = new_pr };
    inline for (@typeInfo(ExchangePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.exchange, field.name)))
            return error.NonFiniteLitterIonExchangePool;
    inline for (@typeInfo(PrecipitatePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.precipitate, field.name)))
            return error.NonFiniteLitterIonExchangePool;
    return result;
}

test "REDIST litter ion exchange accumulates all 14 exchange pools" {
    var ex_tr = std.mem.zeroes(ExchangeTransformations);
    var expected_value: f64 = 1.0;
    inline for (@typeInfo(ExchangeTransformations).@"struct".fields) |field| {
        @field(ex_tr, field.name) = expected_value;
        expected_value += 1.0;
    }
    const result = try apply(
        std.mem.zeroes(ExchangePools),
        std.mem.zeroes(PrecipitatePools),
        ex_tr,
        std.mem.zeroes(PrecipitateTransformations),
    );
    inline for (@typeInfo(ExchangePools).@"struct".fields) |field|
        try std.testing.expectEqual(@field(ex_tr, field.name), @field(result.exchange, field.name));
}

test "REDIST litter precipitate accumulates all 9 precipitate pools" {
    var pr_tr = std.mem.zeroes(PrecipitateTransformations);
    var expected_value: f64 = 1.0;
    inline for (@typeInfo(PrecipitateTransformations).@"struct".fields) |field| {
        @field(pr_tr, field.name) = expected_value;
        expected_value += 1.0;
    }
    const result = try apply(
        std.mem.zeroes(ExchangePools),
        std.mem.zeroes(PrecipitatePools),
        std.mem.zeroes(ExchangeTransformations),
        pr_tr,
    );
    inline for (@typeInfo(PrecipitatePools).@"struct".fields) |field|
        try std.testing.expectEqual(@field(pr_tr, field.name), @field(result.precipitate, field.name));
}

test "REDIST litter ion exchange preserves prior state on zero transformation" {
    var ex = std.mem.zeroes(ExchangePools);
    ex.ca_mol = 7.5;
    const result = try apply(
        ex,
        std.mem.zeroes(PrecipitatePools),
        std.mem.zeroes(ExchangeTransformations),
        std.mem.zeroes(PrecipitateTransformations),
    );
    try std.testing.expectApproxEqRel(@as(f64, 7.5), result.exchange.ca_mol, 1.0e-15);
}

test "REDIST litter ion exchange preserves signed removal increments" {
    var exchange = std.mem.zeroes(ExchangePools);
    exchange.nh4_mol = 4.0;
    var ex_tr = std.mem.zeroes(ExchangeTransformations);
    ex_tr.nh4_mol = -1.5;
    var precipitate = std.mem.zeroes(PrecipitatePools);
    precipitate.alpo4_mol = 3.0;
    var pr_tr = std.mem.zeroes(PrecipitateTransformations);
    pr_tr.alpo4_mol = -0.5;
    const result = try apply(exchange, precipitate, ex_tr, pr_tr);
    try std.testing.expectEqual(@as(f64, 2.5), result.exchange.nh4_mol);
    try std.testing.expectEqual(@as(f64, 2.5), result.precipitate.alpo4_mol);
}

test "REDIST litter ion exchange rejects non-finite transformation" {
    var ex_tr = std.mem.zeroes(ExchangeTransformations);
    ex_tr.oh2_mol = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidLitterIonExchangeInput,
        apply(
            std.mem.zeroes(ExchangePools),
            std.mem.zeroes(PrecipitatePools),
            ex_tr,
            std.mem.zeroes(PrecipitateTransformations),
        ),
    );
}

test "REDIST litter ion exchange rejects arithmetic overflow" {
    var precipitate = std.mem.zeroes(PrecipitatePools);
    precipitate.caso4_mol = std.math.floatMax(f64);
    var pr_tr = std.mem.zeroes(PrecipitateTransformations);
    pr_tr.caso4_mol = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteLitterIonExchangePool,
        apply(std.mem.zeroes(ExchangePools), precipitate, std.mem.zeroes(ExchangeTransformations), pr_tr),
    );
}
