// redist.f lines 8590–9225: simple scalar pool redistribution for every per-layer
// variable that follows the uniform FX/FY pattern: dst += FX*src; src *= FY.
// Water/ice/heat (lines 8564-9107) and plant roots (8794-8898/9228-9270) are
// handled by separate kernels. This file covers N/P solutes, base ions, salt
// ions (ISALTG), band solutes, gaseous gases, aqueous gases, and dissolved
// organic matter.
//
// Mass conservation invariant for every field f:
//   new_dst.f + new_src.f == old_dst.f + old_src.f
//
// Gate conditions (ISALTG≠0, VLNHB>0, L0≠0, etc.) are evaluated by the
// caller; the transfer functions here are always unconditional.

const std = @import("std");

// ─── primitive ──────────────────────────────────────────────────────────────

// Apply the FX/FY transfer to a single scalar pool value pair.
// Returns (new_src, new_dst) where new_dst = dst + fx*src, new_src = src*(1-fx).
pub fn transferScalar(src: f64, dst: f64, fx: f64) struct { src: f64, dst: f64 } {
    return .{
        .src = src * (1.0 - fx),
        .dst = dst + fx * src,
    };
}

// ─── N/P solutes non-band (REDIST 8590–8603 dst, 9111–9116 src) ─────────────

pub const NpSolutesNonBand = struct {
    znh4s: f64, // ZNH4S
    znh3s: f64, // ZNH3S
    zno3s: f64, // ZNO3S
    zno2s: f64, // ZNO2S
    h1po4: f64, // H1PO4
    h2po4: f64, // H2PO4
};

pub const NpSolutesNonBandTransferResult = struct {
    src: NpSolutesNonBand,
    dst: NpSolutesNonBand,
};

pub fn transferNpSolutesNonBand(
    src: NpSolutesNonBand,
    dst: NpSolutesNonBand,
    fx: f64,
) !NpSolutesNonBandTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .znh4s = src.znh4s * fy,
            .znh3s = src.znh3s * fy,
            .zno3s = src.zno3s * fy,
            .zno2s = src.zno2s * fy,
            .h1po4 = src.h1po4 * fy,
            .h2po4 = src.h2po4 * fy,
        },
        .dst = .{
            .znh4s = dst.znh4s + fx * src.znh4s,
            .znh3s = dst.znh3s + fx * src.znh3s,
            .zno3s = dst.zno3s + fx * src.zno3s,
            .zno2s = dst.zno2s + fx * src.zno2s,
            .h1po4 = dst.h1po4 + fx * src.h1po4,
            .h2po4 = dst.h2po4 + fx * src.h2po4,
        },
    };
}

// ─── base ions non-band (REDIST 8607–8622 dst, 9120–9127 src) ───────────────

pub const BaseIonsNonBand = struct {
    zhy: f64, // ZHY  (H+)
    zoh: f64, // ZOH  (OH-)
    zal: f64, // ZAL  (Al3+)
    zfe: f64, // ZFE  (Fe3+)
    zca: f64, // ZCA  (Ca2+)
    zmg: f64, // ZMG  (Mg2+)
    zna: f64, // ZNA  (Na+)
    zka: f64, // ZKA  (K+)
};

pub const BaseIonsNonBandTransferResult = struct {
    src: BaseIonsNonBand,
    dst: BaseIonsNonBand,
};

pub fn transferBaseIonsNonBand(
    src: BaseIonsNonBand,
    dst: BaseIonsNonBand,
    fx: f64,
) !BaseIonsNonBandTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .zhy = src.zhy * fy,
            .zoh = src.zoh * fy,
            .zal = src.zal * fy,
            .zfe = src.zfe * fy,
            .zca = src.zca * fy,
            .zmg = src.zmg * fy,
            .zna = src.zna * fy,
            .zka = src.zka * fy,
        },
        .dst = .{
            .zhy = dst.zhy + fx * src.zhy,
            .zoh = dst.zoh + fx * src.zoh,
            .zal = dst.zal + fx * src.zal,
            .zfe = dst.zfe + fx * src.zfe,
            .zca = dst.zca + fx * src.zca,
            .zmg = dst.zmg + fx * src.zmg,
            .zna = dst.zna + fx * src.zna,
            .zka = dst.zka + fx * src.zka,
        },
    };
}

// ─── salt ions non-band (REDIST 8624–8692 dst, 9129–9163 src; gate: ISALTG≠0) ─

pub const SaltIonsNonBand = struct {
    zso4: f64, // ZSO4
    zcl: f64, // ZCL
    zco3: f64, // ZCO3
    zhco3: f64, // ZHCO3
    zaloh1: f64, // ZALOH1
    zaloh2: f64, // ZALOH2
    zaloh3: f64, // ZALOH3
    zaloh4: f64, // ZALOH4
    zals: f64, // ZALS
    zfeoh1: f64, // ZFEOH1
    zfeoh2: f64, // ZFEOH2
    zfeoh3: f64, // ZFEOH3
    zfeoh4: f64, // ZFEOH4
    zfes: f64, // ZFES
    zcao: f64, // ZCAO
    zcac: f64, // ZCAC
    zcah: f64, // ZCAH
    zcas: f64, // ZCAS
    zmgo: f64, // ZMGO
    zmgc: f64, // ZMGC
    zmgh: f64, // ZMGH
    zmgs: f64, // ZMGS
    znac: f64, // ZNAC
    znas: f64, // ZNAS
    zkas: f64, // ZKAS
    zhysi: f64, // ZHYSI
    h0po4: f64, // H0PO4
    h3po4: f64, // H3PO4
    zfe1p: f64, // ZFE1P
    zfe2p: f64, // ZFE2P
    zca0p: f64, // ZCA0P
    zca1p: f64, // ZCA1P
    zca2p: f64, // ZCA2P
    zmg1p: f64, // ZMG1P
};

pub const SaltIonsNonBandTransferResult = struct {
    src: SaltIonsNonBand,
    dst: SaltIonsNonBand,
};

pub fn transferSaltIonsNonBand(
    src: SaltIonsNonBand,
    dst: SaltIonsNonBand,
    fx: f64,
) !SaltIonsNonBandTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .zso4 = src.zso4 * fy,
            .zcl = src.zcl * fy,
            .zco3 = src.zco3 * fy,
            .zhco3 = src.zhco3 * fy,
            .zaloh1 = src.zaloh1 * fy,
            .zaloh2 = src.zaloh2 * fy,
            .zaloh3 = src.zaloh3 * fy,
            .zaloh4 = src.zaloh4 * fy,
            .zals = src.zals * fy,
            .zfeoh1 = src.zfeoh1 * fy,
            .zfeoh2 = src.zfeoh2 * fy,
            .zfeoh3 = src.zfeoh3 * fy,
            .zfeoh4 = src.zfeoh4 * fy,
            .zfes = src.zfes * fy,
            .zcao = src.zcao * fy,
            .zcac = src.zcac * fy,
            .zcah = src.zcah * fy,
            .zcas = src.zcas * fy,
            .zmgo = src.zmgo * fy,
            .zmgc = src.zmgc * fy,
            .zmgh = src.zmgh * fy,
            .zmgs = src.zmgs * fy,
            .znac = src.znac * fy,
            .znas = src.znas * fy,
            .zkas = src.zkas * fy,
            .zhysi = src.zhysi * fy,
            .h0po4 = src.h0po4 * fy,
            .h3po4 = src.h3po4 * fy,
            .zfe1p = src.zfe1p * fy,
            .zfe2p = src.zfe2p * fy,
            .zca0p = src.zca0p * fy,
            .zca1p = src.zca1p * fy,
            .zca2p = src.zca2p * fy,
            .zmg1p = src.zmg1p * fy,
        },
        .dst = .{
            .zso4 = dst.zso4 + fx * src.zso4,
            .zcl = dst.zcl + fx * src.zcl,
            .zco3 = dst.zco3 + fx * src.zco3,
            .zhco3 = dst.zhco3 + fx * src.zhco3,
            .zaloh1 = dst.zaloh1 + fx * src.zaloh1,
            .zaloh2 = dst.zaloh2 + fx * src.zaloh2,
            .zaloh3 = dst.zaloh3 + fx * src.zaloh3,
            .zaloh4 = dst.zaloh4 + fx * src.zaloh4,
            .zals = dst.zals + fx * src.zals,
            .zfeoh1 = dst.zfeoh1 + fx * src.zfeoh1,
            .zfeoh2 = dst.zfeoh2 + fx * src.zfeoh2,
            .zfeoh3 = dst.zfeoh3 + fx * src.zfeoh3,
            .zfeoh4 = dst.zfeoh4 + fx * src.zfeoh4,
            .zfes = dst.zfes + fx * src.zfes,
            .zcao = dst.zcao + fx * src.zcao,
            .zcac = dst.zcac + fx * src.zcac,
            .zcah = dst.zcah + fx * src.zcah,
            .zcas = dst.zcas + fx * src.zcas,
            .zmgo = dst.zmgo + fx * src.zmgo,
            .zmgc = dst.zmgc + fx * src.zmgc,
            .zmgh = dst.zmgh + fx * src.zmgh,
            .zmgs = dst.zmgs + fx * src.zmgs,
            .znac = dst.znac + fx * src.znac,
            .znas = dst.znas + fx * src.znas,
            .zkas = dst.zkas + fx * src.zkas,
            .zhysi = dst.zhysi + fx * src.zhysi,
            .h0po4 = dst.h0po4 + fx * src.h0po4,
            .h3po4 = dst.h3po4 + fx * src.h3po4,
            .zfe1p = dst.zfe1p + fx * src.zfe1p,
            .zfe2p = dst.zfe2p + fx * src.zfe2p,
            .zca0p = dst.zca0p + fx * src.zca0p,
            .zca1p = dst.zca1p + fx * src.zca1p,
            .zca2p = dst.zca2p + fx * src.zca2p,
            .zmg1p = dst.zmg1p + fx * src.zmg1p,
        },
    };
}

// ─── band N/P solutes (REDIST 8696–8713 dst, 9168–9179 src; gate: L0≠0) ─────
// Sub-gates VLNHB>0 / VLNOB>0 / VLPOB>0 are evaluated by the caller.

pub const NpSolutesBand = struct {
    znh4b: f64, // ZNH4B (gate: VLNHB>0)
    znh3b: f64, // ZNH3B (gate: VLNHB>0)
    zno3b: f64, // ZNO3B (gate: VLNOB>0)
    zno2b: f64, // ZNO2B (gate: VLNOB>0)
    h1pob: f64, // H1POB (gate: VLPOB>0)
    h2pob: f64, // H2POB (gate: VLPOB>0)
};

pub const NpSolutesBandTransferResult = struct {
    src: NpSolutesBand,
    dst: NpSolutesBand,
};

pub fn transferNpSolutesBand(
    src: NpSolutesBand,
    dst: NpSolutesBand,
    fx: f64,
) !NpSolutesBandTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .znh4b = src.znh4b * fy,
            .znh3b = src.znh3b * fy,
            .zno3b = src.zno3b * fy,
            .zno2b = src.zno2b * fy,
            .h1pob = src.h1pob * fy,
            .h2pob = src.h2pob * fy,
        },
        .dst = .{
            .znh4b = dst.znh4b + fx * src.znh4b,
            .znh3b = dst.znh3b + fx * src.znh3b,
            .zno3b = dst.zno3b + fx * src.zno3b,
            .zno2b = dst.zno2b + fx * src.zno2b,
            .h1pob = dst.h1pob + fx * src.h1pob,
            .h2pob = dst.h2pob + fx * src.h2pob,
        },
    };
}

// ─── band salt ions (REDIST 8715–8733 dst, 9180–9190 src; gate: L0≠0 + ISALTG≠0 + VLPOB>0) ─

pub const SaltIonsBand = struct {
    h0pob: f64, // H0POB
    h3pob: f64, // H3POB
    zfe1pb: f64, // ZFE1PB
    zfe2pb: f64, // ZFE2PB
    zca0pb: f64, // ZCA0PB
    zca1pb: f64, // ZCA1PB
    zca2pb: f64, // ZCA2PB
    zmg1pb: f64, // ZMG1PB
};

pub const SaltIonsBandTransferResult = struct {
    src: SaltIonsBand,
    dst: SaltIonsBand,
};

pub fn transferSaltIonsBand(
    src: SaltIonsBand,
    dst: SaltIonsBand,
    fx: f64,
) !SaltIonsBandTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .h0pob = src.h0pob * fy,
            .h3pob = src.h3pob * fy,
            .zfe1pb = src.zfe1pb * fy,
            .zfe2pb = src.zfe2pb * fy,
            .zca0pb = src.zca0pb * fy,
            .zca1pb = src.zca1pb * fy,
            .zca2pb = src.zca2pb * fy,
            .zmg1pb = src.zmg1pb * fy,
        },
        .dst = .{
            .h0pob = dst.h0pob + fx * src.h0pob,
            .h3pob = dst.h3pob + fx * src.h3pob,
            .zfe1pb = dst.zfe1pb + fx * src.zfe1pb,
            .zfe2pb = dst.zfe2pb + fx * src.zfe2pb,
            .zca0pb = dst.zca0pb + fx * src.zca0pb,
            .zca1pb = dst.zca1pb + fx * src.zca1pb,
            .zca2pb = dst.zca2pb + fx * src.zca2pb,
            .zmg1pb = dst.zmg1pb + fx * src.zmg1pb,
        },
    };
}

// ─── gaseous gases (REDIST 8738–8751 dst, 9195–9201 src; gate: L0≠0) ────────

pub const GaseousGases = struct {
    co2g: f64, // CO2G
    ch4g: f64, // CH4G
    oxyg: f64, // OXYG
    z2gg: f64, // Z2GG  (N2)
    z2og: f64, // Z2OG  (N2O)
    znh3g: f64, // ZNH3G (NH3)
    h2gg: f64, // H2GG  (H2)
};

pub const GaseousGasesTransferResult = struct {
    src: GaseousGases,
    dst: GaseousGases,
};

pub fn transferGaseousGases(
    src: GaseousGases,
    dst: GaseousGases,
    fx: f64,
) !GaseousGasesTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .co2g = src.co2g * fy,
            .ch4g = src.ch4g * fy,
            .oxyg = src.oxyg * fy,
            .z2gg = src.z2gg * fy,
            .z2og = src.z2og * fy,
            .znh3g = src.znh3g * fy,
            .h2gg = src.h2gg * fy,
        },
        .dst = .{
            .co2g = dst.co2g + fx * src.co2g,
            .ch4g = dst.ch4g + fx * src.ch4g,
            .oxyg = dst.oxyg + fx * src.oxyg,
            .z2gg = dst.z2gg + fx * src.z2gg,
            .z2og = dst.z2og + fx * src.z2og,
            .znh3g = dst.znh3g + fx * src.znh3g,
            .h2gg = dst.h2gg + fx * src.h2gg,
        },
    };
}

// ─── aqueous gases (REDIST 8756–8767 dst, 9206–9211 src) ────────────────────

pub const AqueousGases = struct {
    co2s: f64, // CO2S
    ch4s: f64, // CH4S
    oxys: f64, // OXYS
    z2gs: f64, // Z2GS  (N2)
    z2os: f64, // Z2OS  (N2O)
    h2gs: f64, // H2GS  (H2)
};

pub const AqueousGasesTransferResult = struct {
    src: AqueousGases,
    dst: AqueousGases,
};

pub fn transferAqueousGases(
    src: AqueousGases,
    dst: AqueousGases,
    fx: f64,
) !AqueousGasesTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .co2s = src.co2s * fy,
            .ch4s = src.ch4s * fy,
            .oxys = src.oxys * fy,
            .z2gs = src.z2gs * fy,
            .z2os = src.z2os * fy,
            .h2gs = src.h2gs * fy,
        },
        .dst = .{
            .co2s = dst.co2s + fx * src.co2s,
            .ch4s = dst.ch4s + fx * src.ch4s,
            .oxys = dst.oxys + fx * src.oxys,
            .z2gs = dst.z2gs + fx * src.z2gs,
            .z2os = dst.z2os + fx * src.z2os,
            .h2gs = dst.h2gs + fx * src.h2gs,
        },
    };
}

// ─── dissolved organic matter — one layer (K) (REDIST 8772–8788 dst, 9216–9224 src) ─
// Fortran DO 7780 K=0,4 loop; call once per K. Gate: IFLGL(L,3)==0 checked by caller.

pub const DissolvedOrganicMatterLayer = struct {
    oqc: f64, // OQC  (C)
    oqn: f64, // OQN  (N)
    oqp: f64, // OQP  (P)
    oqa: f64, // OQA  (charge)
    oqch: f64, // OQCH (C in humus)
    oqnh: f64, // OQNH (N in humus)
    oqph: f64, // OQPH (P in humus)
    oqah: f64, // OQAH (charge in humus)
};

pub const DissolvedOrganicMatterLayerTransferResult = struct {
    src: DissolvedOrganicMatterLayer,
    dst: DissolvedOrganicMatterLayer,
};

pub fn transferDissolvedOrganicMatterLayer(
    src: DissolvedOrganicMatterLayer,
    dst: DissolvedOrganicMatterLayer,
    fx: f64,
) !DissolvedOrganicMatterLayerTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .oqc = src.oqc * fy,
            .oqn = src.oqn * fy,
            .oqp = src.oqp * fy,
            .oqa = src.oqa * fy,
            .oqch = src.oqch * fy,
            .oqnh = src.oqnh * fy,
            .oqph = src.oqph * fy,
            .oqah = src.oqah * fy,
        },
        .dst = .{
            .oqc = dst.oqc + fx * src.oqc,
            .oqn = dst.oqn + fx * src.oqn,
            .oqp = dst.oqp + fx * src.oqp,
            .oqa = dst.oqa + fx * src.oqa,
            .oqch = dst.oqch + fx * src.oqch,
            .oqnh = dst.oqnh + fx * src.oqnh,
            .oqph = dst.oqph + fx * src.oqph,
            .oqah = dst.oqah + fx * src.oqah,
        },
    };
}

// ─── tests ───────────────────────────────────────────────────────────────────

test "transferScalar mass conservation" {
    const s = transferScalar(10.0, 3.0, 0.4);
    try std.testing.expectApproxEqRel(s.dst + s.src, 13.0, 1e-12);
    try std.testing.expectApproxEqRel(s.dst, 7.0, 1e-12);
    try std.testing.expectApproxEqRel(s.src, 6.0, 1e-12);
}

test "transferScalar zero fraction" {
    const s = transferScalar(10.0, 3.0, 0.0);
    try std.testing.expectEqual(s.src, 10.0);
    try std.testing.expectEqual(s.dst, 3.0);
}

test "transferScalar full fraction" {
    const s = transferScalar(10.0, 3.0, 1.0);
    try std.testing.expectApproxEqRel(s.src, 0.0, 1e-12);
    try std.testing.expectApproxEqRel(s.dst, 13.0, 1e-12);
}

test "transferNpSolutesNonBand mass conservation" {
    const src = NpSolutesNonBand{ .znh4s = 5.0, .znh3s = 1.0, .zno3s = 2.0, .zno2s = 0.5, .h1po4 = 3.0, .h2po4 = 0.1 };
    const dst = NpSolutesNonBand{ .znh4s = 1.0, .znh3s = 0.0, .zno3s = 0.5, .zno2s = 0.0, .h1po4 = 0.2, .h2po4 = 0.0 };
    const r = try transferNpSolutesNonBand(src, dst, 0.3);
    try std.testing.expectApproxEqRel(r.src.znh4s + r.dst.znh4s, src.znh4s + dst.znh4s, 1e-12);
    try std.testing.expectApproxEqRel(r.src.h2po4 + r.dst.h2po4, src.h2po4 + dst.h2po4, 1e-12);
}

test "transferNpSolutesNonBand invalid fraction" {
    const zero = NpSolutesNonBand{ .znh4s = 0, .znh3s = 0, .zno3s = 0, .zno2s = 0, .h1po4 = 0, .h2po4 = 0 };
    try std.testing.expectError(error.InvalidFraction, transferNpSolutesNonBand(zero, zero, 1.1));
    try std.testing.expectError(error.InvalidFraction, transferNpSolutesNonBand(zero, zero, -0.01));
}

test "transferBaseIonsNonBand mass conservation" {
    const src = BaseIonsNonBand{ .zhy = 1.0, .zoh = 2.0, .zal = 0.1, .zfe = 0.2, .zca = 4.0, .zmg = 0.5, .zna = 3.0, .zka = 1.5 };
    const dst = BaseIonsNonBand{ .zhy = 0.5, .zoh = 0.0, .zal = 0.0, .zfe = 0.1, .zca = 1.0, .zmg = 0.0, .zna = 0.5, .zka = 0.2 };
    const r = try transferBaseIonsNonBand(src, dst, 0.5);
    try std.testing.expectApproxEqRel(r.src.zca + r.dst.zca, src.zca + dst.zca, 1e-12);
    try std.testing.expectApproxEqRel(r.src.zna + r.dst.zna, src.zna + dst.zna, 1e-12);
}

test "transferSaltIonsNonBand mass conservation" {
    var src = std.mem.zeroes(SaltIonsNonBand);
    var dst = std.mem.zeroes(SaltIonsNonBand);
    src.zso4 = 3.0;
    dst.zso4 = 1.0;
    src.zmg1p = 2.5;
    dst.zmg1p = 0.5;
    const r = try transferSaltIonsNonBand(src, dst, 0.6);
    try std.testing.expectApproxEqRel(r.src.zso4 + r.dst.zso4, src.zso4 + dst.zso4, 1e-12);
    try std.testing.expectApproxEqRel(r.src.zmg1p + r.dst.zmg1p, src.zmg1p + dst.zmg1p, 1e-12);
}

test "transferNpSolutesBand mass conservation" {
    const src = NpSolutesBand{ .znh4b = 2.0, .znh3b = 0.5, .zno3b = 1.0, .zno2b = 0.2, .h1pob = 0.8, .h2pob = 0.1 };
    const dst = NpSolutesBand{ .znh4b = 0.5, .znh3b = 0.0, .zno3b = 0.0, .zno2b = 0.0, .h1pob = 0.1, .h2pob = 0.0 };
    const r = try transferNpSolutesBand(src, dst, 0.25);
    try std.testing.expectApproxEqRel(r.src.znh4b + r.dst.znh4b, src.znh4b + dst.znh4b, 1e-12);
    try std.testing.expectApproxEqRel(r.src.h2pob + r.dst.h2pob, src.h2pob + dst.h2pob, 1e-12);
}

test "transferGaseousGases mass conservation" {
    const src = GaseousGases{ .co2g = 10.0, .ch4g = 2.0, .oxyg = 5.0, .z2gg = 0.5, .z2og = 0.1, .znh3g = 0.3, .h2gg = 0.05 };
    const dst = GaseousGases{ .co2g = 1.0, .ch4g = 0.0, .oxyg = 0.5, .z2gg = 0.0, .z2og = 0.0, .znh3g = 0.0, .h2gg = 0.0 };
    const r = try transferGaseousGases(src, dst, 0.4);
    try std.testing.expectApproxEqRel(r.src.co2g + r.dst.co2g, src.co2g + dst.co2g, 1e-12);
    try std.testing.expectApproxEqRel(r.src.oxyg + r.dst.oxyg, src.oxyg + dst.oxyg, 1e-12);
}

test "transferAqueousGases mass conservation" {
    const src = AqueousGases{ .co2s = 4.0, .ch4s = 1.0, .oxys = 2.0, .z2gs = 0.3, .z2os = 0.1, .h2gs = 0.05 };
    const dst = AqueousGases{ .co2s = 0.5, .ch4s = 0.0, .oxys = 0.2, .z2gs = 0.0, .z2os = 0.0, .h2gs = 0.0 };
    const r = try transferAqueousGases(src, dst, 0.7);
    try std.testing.expectApproxEqRel(r.src.co2s + r.dst.co2s, src.co2s + dst.co2s, 1e-12);
    try std.testing.expectApproxEqRel(r.src.h2gs + r.dst.h2gs, src.h2gs + dst.h2gs, 1e-12);
}

test "transferDissolvedOrganicMatterLayer mass conservation" {
    const src = DissolvedOrganicMatterLayer{ .oqc = 8.0, .oqn = 1.0, .oqp = 0.5, .oqa = 0.2, .oqch = 3.0, .oqnh = 0.4, .oqph = 0.1, .oqah = 0.05 };
    const dst = DissolvedOrganicMatterLayer{ .oqc = 2.0, .oqn = 0.0, .oqp = 0.0, .oqa = 0.0, .oqch = 1.0, .oqnh = 0.0, .oqph = 0.0, .oqah = 0.0 };
    const r = try transferDissolvedOrganicMatterLayer(src, dst, 0.5);
    try std.testing.expectApproxEqRel(r.src.oqc + r.dst.oqc, src.oqc + dst.oqc, 1e-12);
    try std.testing.expectApproxEqRel(r.src.oqah + r.dst.oqah, src.oqah + dst.oqah, 1e-12);
}

test "transferDissolvedOrganicMatterLayer invalid fraction" {
    const zero = std.mem.zeroes(DissolvedOrganicMatterLayer);
    try std.testing.expectError(error.InvalidFraction, transferDissolvedOrganicMatterLayer(zero, zero, -0.1));
}
