// redist.f lines 8903–9080 (dst) and 9296–9400 (src):
// solid/mineral pool redistribution, all gated on FX==1.0 (full-layer transfer).
// Sediment additionally requires L0≠0 (caller enforces).
// Solid organic matter additionally requires IFLGL(L,3)==0 (caller enforces).
// Plant roots (lines 8794-8898 / 9228-9291) are handled by a separate kernel.
//
// Mass conservation invariant: for every field f, new_dst.f + new_src.f == old_dst.f + old_src.f

const std = @import("std");

// ─── sediment + exchange sites + silicates ────────────────────────────────────
// redist.f lines 8904-8963 dst / 9297-9327 src; gate: FX==1.0 AND L0≠0

pub const SedimentLayer = struct {
    sand: f64, // SAND
    silt: f64, // SILT
    clay: f64, // CLAY
    xcec: f64, // XCEC (cation exchange capacity)
    xaec: f64, // XAEC (anion exchange capacity)
    xhy: f64, // XHY  (exchangeable H+)
    xal: f64, // XAL
    xfe: f64, // XFE
    xca: f64, // XCA
    xmg: f64, // XMG
    xna: f64, // XNA
    xka: f64, // XKA
    xhc: f64, // XHC  (humus-C)
    paloh: f64, // PALOH (Al-hydroxide precipitate)
    pfeoh: f64, // PFEOH (Fe-hydroxide precipitate)
    pcaco: f64, // PCACO (calcite)
    pcaso: f64, // PCASO (gypsum)
    qalsi: f64, // QALSI
    qfesi: f64, // QFESI
    qcasi: f64, // QCASI
    qmgsi: f64, // QMGSI
    qnasi: f64, // QNASI
    qkasi: f64, // QKASI
    qalsif: f64, // QALSIF
    qfesif: f64, // QFESIF
    qcasif: f64, // QCASIF
    qmgsif: f64, // QMGSIF
    qnasif: f64, // QNASIF
    qkasif: f64, // QKASIF
};

pub const SedimentLayerTransferResult = struct { src: SedimentLayer, dst: SedimentLayer };

pub fn transferSedimentLayer(
    src: SedimentLayer,
    dst: SedimentLayer,
    fx: f64,
) !SedimentLayerTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    // This kernel is called only when FX==1.0 and L0≠0 (Fortran gate).
    // The math is identical to the scalar pattern regardless.
    return .{
        .src = .{
            .sand = src.sand * fy,
            .silt = src.silt * fy,
            .clay = src.clay * fy,
            .xcec = src.xcec * fy,
            .xaec = src.xaec * fy,
            .xhy = src.xhy * fy,
            .xal = src.xal * fy,
            .xfe = src.xfe * fy,
            .xca = src.xca * fy,
            .xmg = src.xmg * fy,
            .xna = src.xna * fy,
            .xka = src.xka * fy,
            .xhc = src.xhc * fy,
            .paloh = src.paloh * fy,
            .pfeoh = src.pfeoh * fy,
            .pcaco = src.pcaco * fy,
            .pcaso = src.pcaso * fy,
            .qalsi = src.qalsi * fy,
            .qfesi = src.qfesi * fy,
            .qcasi = src.qcasi * fy,
            .qmgsi = src.qmgsi * fy,
            .qnasi = src.qnasi * fy,
            .qkasi = src.qkasi * fy,
            .qalsif = src.qalsif * fy,
            .qfesif = src.qfesif * fy,
            .qcasif = src.qcasif * fy,
            .qmgsif = src.qmgsif * fy,
            .qnasif = src.qnasif * fy,
            .qkasif = src.qkasif * fy,
        },
        .dst = .{
            .sand = dst.sand + fx * src.sand,
            .silt = dst.silt + fx * src.silt,
            .clay = dst.clay + fx * src.clay,
            .xcec = dst.xcec + fx * src.xcec,
            .xaec = dst.xaec + fx * src.xaec,
            .xhy = dst.xhy + fx * src.xhy,
            .xal = dst.xal + fx * src.xal,
            .xfe = dst.xfe + fx * src.xfe,
            .xca = dst.xca + fx * src.xca,
            .xmg = dst.xmg + fx * src.xmg,
            .xna = dst.xna + fx * src.xna,
            .xka = dst.xka + fx * src.xka,
            .xhc = dst.xhc + fx * src.xhc,
            .paloh = dst.paloh + fx * src.paloh,
            .pfeoh = dst.pfeoh + fx * src.pfeoh,
            .pcaco = dst.pcaco + fx * src.pcaco,
            .pcaso = dst.pcaso + fx * src.pcaso,
            .qalsi = dst.qalsi + fx * src.qalsi,
            .qfesi = dst.qfesi + fx * src.qfesi,
            .qcasi = dst.qcasi + fx * src.qcasi,
            .qmgsi = dst.qmgsi + fx * src.qmgsi,
            .qnasi = dst.qnasi + fx * src.qnasi,
            .qkasi = dst.qkasi + fx * src.qkasi,
            .qalsif = dst.qalsif + fx * src.qalsif,
            .qfesif = dst.qfesif + fx * src.qfesif,
            .qcasif = dst.qcasif + fx * src.qcasif,
            .qmgsif = dst.qmgsif + fx * src.qmgsif,
            .qnasif = dst.qnasif + fx * src.qnasif,
            .qkasif = dst.qkasif + fx * src.qkasif,
        },
    };
}

// ─── fertilizer N (REDIST 8967–8982 dst / 9331–9338 src; gate: FX==1.0) ─────

pub const FertilizerNLayer = struct {
    znh4fa: f64, // ZNH4FA (NH4 non-band A)
    znh3fa: f64, // ZNH3FA (NH3 non-band A)
    znhufa: f64, // ZNHUFA (urea non-band A)
    zno3fa: f64, // ZNO3FA (NO3 non-band A)
    znh4fb: f64, // ZNH4FB (NH4 band B)
    znh3fb: f64, // ZNH3FB (NH3 band B)
    znhufb: f64, // ZNHUFB (urea band B)
    zno3fb: f64, // ZNO3FB (NO3 band B)
};

pub const FertilizerNLayerTransferResult = struct { src: FertilizerNLayer, dst: FertilizerNLayer };

pub fn transferFertilizerNLayer(
    src: FertilizerNLayer,
    dst: FertilizerNLayer,
    fx: f64,
) !FertilizerNLayerTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .znh4fa = src.znh4fa * fy,
            .znh3fa = src.znh3fa * fy,
            .znhufa = src.znhufa * fy,
            .zno3fa = src.zno3fa * fy,
            .znh4fb = src.znh4fb * fy,
            .znh3fb = src.znh3fb * fy,
            .znhufb = src.znhufb * fy,
            .zno3fb = src.zno3fb * fy,
        },
        .dst = .{
            .znh4fa = dst.znh4fa + fx * src.znh4fa,
            .znh3fa = dst.znh3fa + fx * src.znh3fa,
            .znhufa = dst.znhufa + fx * src.znhufa,
            .zno3fa = dst.zno3fa + fx * src.zno3fa,
            .znh4fb = dst.znh4fb + fx * src.znh4fb,
            .znh3fb = dst.znh3fb + fx * src.znh3fb,
            .znhufb = dst.znhufb + fx * src.znhufb,
            .zno3fb = dst.zno3fb + fx * src.zno3fb,
        },
    };
}

// ─── adsorbed cations/anions (REDIST 8986–9012 dst / 9342–9356 src; gate: FX==1.0) ─

pub const AdsorbedIonsLayer = struct {
    xn4: f64, // XN4   (NH4+ adsorbed non-band)
    xnb: f64, // XNB   (NH4+ adsorbed band)
    xoh0: f64, // XOH0  (OH- adsorbed)
    xoh1: f64, // XOH1
    xoh2: f64, // XOH2
    xh1p: f64, // XH1P  (H2PO4- adsorbed)
    xh2p: f64, // XH2P  (HPO4-- adsorbed)
    xoh0b: f64, // XOH0B (band)
    xoh1b: f64, // XOH1B
    xoh2b: f64, // XOH2B
    xh1pb: f64, // XH1PB
    xh2pb: f64, // XH2PB
};

pub const AdsorbedIonsLayerTransferResult = struct { src: AdsorbedIonsLayer, dst: AdsorbedIonsLayer };

pub fn transferAdsorbedIonsLayer(
    src: AdsorbedIonsLayer,
    dst: AdsorbedIonsLayer,
    fx: f64,
) !AdsorbedIonsLayerTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .xn4 = src.xn4 * fy,
            .xnb = src.xnb * fy,
            .xoh0 = src.xoh0 * fy,
            .xoh1 = src.xoh1 * fy,
            .xoh2 = src.xoh2 * fy,
            .xh1p = src.xh1p * fy,
            .xh2p = src.xh2p * fy,
            .xoh0b = src.xoh0b * fy,
            .xoh1b = src.xoh1b * fy,
            .xoh2b = src.xoh2b * fy,
            .xh1pb = src.xh1pb * fy,
            .xh2pb = src.xh2pb * fy,
        },
        .dst = .{
            .xn4 = dst.xn4 + fx * src.xn4,
            .xnb = dst.xnb + fx * src.xnb,
            .xoh0 = dst.xoh0 + fx * src.xoh0,
            .xoh1 = dst.xoh1 + fx * src.xoh1,
            .xoh2 = dst.xoh2 + fx * src.xoh2,
            .xh1p = dst.xh1p + fx * src.xh1p,
            .xh2p = dst.xh2p + fx * src.xh2p,
            .xoh0b = dst.xoh0b + fx * src.xoh0b,
            .xoh1b = dst.xoh1b + fx * src.xoh1b,
            .xoh2b = dst.xoh2b + fx * src.xoh2b,
            .xh1pb = dst.xh1pb + fx * src.xh1pb,
            .xh2pb = dst.xh2pb + fx * src.xh2pb,
        },
    };
}

// ─── precipitates (REDIST 9016–9035 dst / 9360–9369 src; gate: FX==1.0) ─────

pub const PrecipitatesLayer = struct {
    palpo: f64, // PALPO (Al-phosphate non-band)
    pfepo: f64, // PFEPO (Fe-phosphate non-band)
    pcapd: f64, // PCAPD (dicalcium phosphate)
    pcaph: f64, // PCAPH (hydroxyapatite)
    pcapm: f64, // PCAPM (monocalcium phosphate)
    palpb: f64, // PALPB (Al-phosphate band)
    pfepb: f64, // PFEPB
    pcpdb: f64, // PCPDB
    pcphb: f64, // PCPHB
    pcpmb: f64, // PCPMB
};

pub const PrecipitatesLayerTransferResult = struct { src: PrecipitatesLayer, dst: PrecipitatesLayer };

pub fn transferPrecipitatesLayer(
    src: PrecipitatesLayer,
    dst: PrecipitatesLayer,
    fx: f64,
) !PrecipitatesLayerTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .palpo = src.palpo * fy,
            .pfepo = src.pfepo * fy,
            .pcapd = src.pcapd * fy,
            .pcaph = src.pcaph * fy,
            .pcapm = src.pcapm * fy,
            .palpb = src.palpb * fy,
            .pfepb = src.pfepb * fy,
            .pcpdb = src.pcpdb * fy,
            .pcphb = src.pcphb * fy,
            .pcpmb = src.pcpmb * fy,
        },
        .dst = .{
            .palpo = dst.palpo + fx * src.palpo,
            .pfepo = dst.pfepo + fx * src.pfepo,
            .pcapd = dst.pcapd + fx * src.pcapd,
            .pcaph = dst.pcaph + fx * src.pcaph,
            .pcapm = dst.pcapm + fx * src.pcapm,
            .palpb = dst.palpb + fx * src.palpb,
            .pfepb = dst.pfepb + fx * src.pfepb,
            .pcpdb = dst.pcpdb + fx * src.pcpdb,
            .pcphb = dst.pcphb + fx * src.pcphb,
            .pcpmb = dst.pcpmb + fx * src.pcpmb,
        },
    };
}

// ─── solid organic matter (REDIST 9039–9079 dst / 9373–9399 src) ─────────────
// Gate: FX==1.0 AND IFLGL(L,3)==0 (caller enforces).
// Array dimensions (0-based in Zig, matching Fortran loop bounds):
//   OMC/OMN/OMP: [M_count=3][N_count=7][K_count=6]  (Fortran M=1..3, N=1..7, K=0..5)
//   ORC/ORN/ORP: [M_count=2][K_count=5]              (Fortran M=1..2, K=0..4)
//   OHC/OHN/OHP/OHA: [K_count=5]                     (Fortran K=0..4)
//   OSC/OSA/OSN/OSP: [M_count=5][K_count=5]          (Fortran M=1..5, K=0..4)

pub const SolidOrganicMatterLayer = struct {
    orgc: f64,
    omc: [3][7][6]f64,
    omn: [3][7][6]f64,
    omp: [3][7][6]f64,
    orc: [2][5]f64,
    orn: [2][5]f64,
    orp: [2][5]f64,
    ohc: [5]f64,
    ohn: [5]f64,
    ohp: [5]f64,
    oha: [5]f64,
    osc: [5][5]f64,
    osa: [5][5]f64,
    osn: [5][5]f64,
    osp: [5][5]f64,
};

pub const SolidOrganicMatterLayerTransferResult = struct {
    src: SolidOrganicMatterLayer,
    dst: SolidOrganicMatterLayer,
};

pub fn transferSolidOrganicMatterLayer(
    src: SolidOrganicMatterLayer,
    dst: SolidOrganicMatterLayer,
    fx: f64,
) !SolidOrganicMatterLayerTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;

    var result: SolidOrganicMatterLayerTransferResult = undefined;
    result.src.orgc = src.orgc * fy;
    result.dst.orgc = dst.orgc + fx * src.orgc;

    // OMC, OMN, OMP: [3][7][6]
    for (0..3) |m| {
        for (0..7) |n| {
            for (0..6) |k| {
                result.src.omc[m][n][k] = src.omc[m][n][k] * fy;
                result.dst.omc[m][n][k] = dst.omc[m][n][k] + fx * src.omc[m][n][k];
                result.src.omn[m][n][k] = src.omn[m][n][k] * fy;
                result.dst.omn[m][n][k] = dst.omn[m][n][k] + fx * src.omn[m][n][k];
                result.src.omp[m][n][k] = src.omp[m][n][k] * fy;
                result.dst.omp[m][n][k] = dst.omp[m][n][k] + fx * src.omp[m][n][k];
            }
        }
    }

    // ORC, ORN, ORP: [2][5]
    for (0..2) |m| {
        for (0..5) |k| {
            result.src.orc[m][k] = src.orc[m][k] * fy;
            result.dst.orc[m][k] = dst.orc[m][k] + fx * src.orc[m][k];
            result.src.orn[m][k] = src.orn[m][k] * fy;
            result.dst.orn[m][k] = dst.orn[m][k] + fx * src.orn[m][k];
            result.src.orp[m][k] = src.orp[m][k] * fy;
            result.dst.orp[m][k] = dst.orp[m][k] + fx * src.orp[m][k];
        }
    }

    // OHC, OHN, OHP, OHA: [5]
    for (0..5) |k| {
        result.src.ohc[k] = src.ohc[k] * fy;
        result.dst.ohc[k] = dst.ohc[k] + fx * src.ohc[k];
        result.src.ohn[k] = src.ohn[k] * fy;
        result.dst.ohn[k] = dst.ohn[k] + fx * src.ohn[k];
        result.src.ohp[k] = src.ohp[k] * fy;
        result.dst.ohp[k] = dst.ohp[k] + fx * src.ohp[k];
        result.src.oha[k] = src.oha[k] * fy;
        result.dst.oha[k] = dst.oha[k] + fx * src.oha[k];
    }

    // OSC, OSA, OSN, OSP: [5][5]
    for (0..5) |m| {
        for (0..5) |k| {
            result.src.osc[m][k] = src.osc[m][k] * fy;
            result.dst.osc[m][k] = dst.osc[m][k] + fx * src.osc[m][k];
            result.src.osa[m][k] = src.osa[m][k] * fy;
            result.dst.osa[m][k] = dst.osa[m][k] + fx * src.osa[m][k];
            result.src.osn[m][k] = src.osn[m][k] * fy;
            result.dst.osn[m][k] = dst.osn[m][k] + fx * src.osn[m][k];
            result.src.osp[m][k] = src.osp[m][k] * fy;
            result.dst.osp[m][k] = dst.osp[m][k] + fx * src.osp[m][k];
        }
    }

    return result;
}

// ─── tests ───────────────────────────────────────────────────────────────────

test "transferSedimentLayer mass conservation" {
    var src = std.mem.zeroes(SedimentLayer);
    var dst = std.mem.zeroes(SedimentLayer);
    src.sand = 10.0;
    dst.sand = 5.0;
    src.clay = 3.0;
    dst.clay = 1.0;
    src.xcec = 2.0;
    dst.xcec = 0.5;
    const r = try transferSedimentLayer(src, dst, 1.0);
    try std.testing.expectApproxEqRel(r.src.sand + r.dst.sand, src.sand + dst.sand, 1e-12);
    try std.testing.expectApproxEqRel(r.src.clay + r.dst.clay, src.clay + dst.clay, 1e-12);
    try std.testing.expectApproxEqRel(r.src.xcec + r.dst.xcec, src.xcec + dst.xcec, 1e-12);
}

test "transferSedimentLayer full transfer FX=1" {
    var src = std.mem.zeroes(SedimentLayer);
    src.sand = 7.0;
    const dst = std.mem.zeroes(SedimentLayer);
    const r = try transferSedimentLayer(src, dst, 1.0);
    try std.testing.expectApproxEqRel(r.src.sand, 0.0, 1e-12);
    try std.testing.expectApproxEqRel(r.dst.sand, 7.0, 1e-12);
}

test "transferSedimentLayer invalid fraction" {
    const s = std.mem.zeroes(SedimentLayer);
    try std.testing.expectError(error.InvalidFraction, transferSedimentLayer(s, s, -0.1));
    try std.testing.expectError(error.InvalidFraction, transferSedimentLayer(s, s, 1.1));
}

test "transferFertilizerNLayer mass conservation" {
    const src = FertilizerNLayer{
        .znh4fa = 5.0,
        .znh3fa = 1.0,
        .znhufa = 0.5,
        .zno3fa = 2.0,
        .znh4fb = 3.0,
        .znh3fb = 0.5,
        .znhufb = 0.2,
        .zno3fb = 1.0,
    };
    const dst = FertilizerNLayer{
        .znh4fa = 1.0,
        .znh3fa = 0.0,
        .znhufa = 0.0,
        .zno3fa = 0.5,
        .znh4fb = 0.5,
        .znh3fb = 0.0,
        .znhufb = 0.0,
        .zno3fb = 0.2,
    };
    const r = try transferFertilizerNLayer(src, dst, 1.0);
    try std.testing.expectApproxEqRel(r.src.znh4fa + r.dst.znh4fa, src.znh4fa + dst.znh4fa, 1e-12);
    try std.testing.expectApproxEqRel(r.src.zno3fb + r.dst.zno3fb, src.zno3fb + dst.zno3fb, 1e-12);
}

test "transferAdsorbedIonsLayer mass conservation" {
    var src = std.mem.zeroes(AdsorbedIonsLayer);
    var dst = std.mem.zeroes(AdsorbedIonsLayer);
    src.xn4 = 4.0;
    dst.xn4 = 1.0;
    src.xh2pb = 0.5;
    dst.xh2pb = 0.1;
    const r = try transferAdsorbedIonsLayer(src, dst, 1.0);
    try std.testing.expectApproxEqRel(r.src.xn4 + r.dst.xn4, src.xn4 + dst.xn4, 1e-12);
    try std.testing.expectApproxEqRel(r.src.xh2pb + r.dst.xh2pb, src.xh2pb + dst.xh2pb, 1e-12);
}

test "transferPrecipitatesLayer mass conservation" {
    var src = std.mem.zeroes(PrecipitatesLayer);
    var dst = std.mem.zeroes(PrecipitatesLayer);
    src.palpo = 2.0;
    dst.palpo = 0.5;
    src.pcpmb = 1.0;
    dst.pcpmb = 0.2;
    const r = try transferPrecipitatesLayer(src, dst, 1.0);
    try std.testing.expectApproxEqRel(r.src.palpo + r.dst.palpo, src.palpo + dst.palpo, 1e-12);
    try std.testing.expectApproxEqRel(r.src.pcpmb + r.dst.pcpmb, src.pcpmb + dst.pcpmb, 1e-12);
}

test "transferSolidOrganicMatterLayer mass conservation ORGC" {
    var src = std.mem.zeroes(SolidOrganicMatterLayer);
    var dst = std.mem.zeroes(SolidOrganicMatterLayer);
    src.orgc = 12.0;
    dst.orgc = 3.0;
    const r = try transferSolidOrganicMatterLayer(src, dst, 1.0);
    try std.testing.expectApproxEqRel(r.src.orgc + r.dst.orgc, src.orgc + dst.orgc, 1e-12);
    try std.testing.expectApproxEqRel(r.src.orgc, 0.0, 1e-12);
    try std.testing.expectApproxEqRel(r.dst.orgc, 15.0, 1e-12);
}

test "transferSolidOrganicMatterLayer mass conservation OMC" {
    var src = std.mem.zeroes(SolidOrganicMatterLayer);
    var dst = std.mem.zeroes(SolidOrganicMatterLayer);
    src.omc[1][3][2] = 5.0;
    dst.omc[1][3][2] = 2.0;
    src.osn[4][3] = 1.5;
    dst.osn[4][3] = 0.3;
    const r = try transferSolidOrganicMatterLayer(src, dst, 1.0);
    const total_omc = src.omc[1][3][2] + dst.omc[1][3][2];
    try std.testing.expectApproxEqRel(r.src.omc[1][3][2] + r.dst.omc[1][3][2], total_omc, 1e-12);
    const total_osn = src.osn[4][3] + dst.osn[4][3];
    try std.testing.expectApproxEqRel(r.src.osn[4][3] + r.dst.osn[4][3], total_osn, 1e-12);
}

test "transferSolidOrganicMatterLayer partial FX=0.6" {
    var src = std.mem.zeroes(SolidOrganicMatterLayer);
    var dst = std.mem.zeroes(SolidOrganicMatterLayer);
    src.orgc = 10.0;
    dst.orgc = 4.0;
    const r = try transferSolidOrganicMatterLayer(src, dst, 0.6);
    try std.testing.expectApproxEqRel(r.src.orgc + r.dst.orgc, src.orgc + dst.orgc, 1e-12);
    try std.testing.expectApproxEqRel(r.src.orgc, 4.0, 1e-12); // 10 * 0.4
    try std.testing.expectApproxEqRel(r.dst.orgc, 10.0, 1e-12); // 4 + 6
}

test "transferSolidOrganicMatterLayer invalid fraction" {
    const s = std.mem.zeroes(SolidOrganicMatterLayer);
    try std.testing.expectError(error.InvalidFraction, transferSolidOrganicMatterLayer(s, s, -0.01));
}
