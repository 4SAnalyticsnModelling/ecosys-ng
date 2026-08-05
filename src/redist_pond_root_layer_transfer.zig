// REDIST lines 8794–8899 (dst) and 9230–9292 (src):
// Plant root and root-nodule pool redistribution during pond/soil relayering.
// Gate enforced by caller: L0≠0, WTRTL(1,L0,NZ,…)>ZEROP, WTRTL(1,L1,NZ,…)>ZEROP.
//
// Fortran loop structure:
//   DO 8900 NZ=1,NP                     ← plant population
//     IF(WTRTL gate)THEN
//       DO 8895 N=1,MY(NZ)              ← root age class
//         [root gases — per N, per layer]
//         DO 8870 NR=1,NRT(NZ)         ← root segment type
//           [root segment pools]
//         [root biomass — per N, per layer]
//       [root nodules — per NZ, per layer]
//
// Caller loops NZ, N, NR; kernel functions take a single scalar layer pair.
// All functions follow the uniform dst += FX*src; src *= FY pattern.

const std = @import("std");

// ─── root gases (lines 8802–8825 dst / 9238–9249 src; per N, per layer) ─────
// Variables: CO2A, OXYA, CH4A, Z2OA, ZH3A, H2GA (aerenchyma)
//            CO2P, OXYP, CH4P, Z2OP, ZH3P, H2GP (prism/root cell)

pub const RootGasesLayer = struct {
    co2a: f64, // CO2A (aerenchyma CO2)
    oxya: f64, // OXYA
    ch4a: f64, // CH4A
    z2oa: f64, // Z2OA (N2O aerenchyma)
    zh3a: f64, // ZH3A (NH3 aerenchyma)
    h2ga: f64, // H2GA (H2 aerenchyma)
    co2p: f64, // CO2P (prism CO2)
    oxyp: f64, // OXYP
    ch4p: f64, // CH4P
    z2op: f64, // Z2OP
    zh3p: f64, // ZH3P
    h2gp: f64, // H2GP
};

pub const RootGasesLayerTransferResult = struct { src: RootGasesLayer, dst: RootGasesLayer };

pub fn transferRootGasesLayer(
    src: RootGasesLayer,
    dst: RootGasesLayer,
    fx: f64,
) !RootGasesLayerTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .co2a = src.co2a * fy,
            .oxya = src.oxya * fy,
            .ch4a = src.ch4a * fy,
            .z2oa = src.z2oa * fy,
            .zh3a = src.zh3a * fy,
            .h2ga = src.h2ga * fy,
            .co2p = src.co2p * fy,
            .oxyp = src.oxyp * fy,
            .ch4p = src.ch4p * fy,
            .z2op = src.z2op * fy,
            .zh3p = src.zh3p * fy,
            .h2gp = src.h2gp * fy,
        },
        .dst = .{
            .co2a = dst.co2a + fx * src.co2a,
            .oxya = dst.oxya + fx * src.oxya,
            .ch4a = dst.ch4a + fx * src.ch4a,
            .z2oa = dst.z2oa + fx * src.z2oa,
            .zh3a = dst.zh3a + fx * src.zh3a,
            .h2ga = dst.h2ga + fx * src.h2ga,
            .co2p = dst.co2p + fx * src.co2p,
            .oxyp = dst.oxyp + fx * src.oxyp,
            .ch4p = dst.ch4p + fx * src.ch4p,
            .z2op = dst.z2op + fx * src.z2op,
            .zh3p = dst.zh3p + fx * src.zh3p,
            .h2gp = dst.h2gp + fx * src.h2gp,
        },
    };
}

// ─── root segment pools (lines 8830–8847 dst / 9254–9261 src; per N, per NR, per layer) ─

pub const RootSegmentLayer = struct {
    wtrt1: f64, // WTRT1  (fine root C)
    wtrt1n: f64, // WTRT1N (fine root N)
    wtrt1p: f64, // WTRT1P (fine root P)
    wtrt2: f64, // WTRT2  (coarse root C)
    wtrt2n: f64, // WTRT2N
    wtrt2p: f64, // WTRT2P
    rtlg1: f64, // RTLG1  (root length 1)
    rtlg2: f64, // RTLG2  (root length 2)
    rtn2: f64, // RTN2   (root N type 2)
};

pub const RootSegmentLayerTransferResult = struct { src: RootSegmentLayer, dst: RootSegmentLayer };

pub fn transferRootSegmentLayer(
    src: RootSegmentLayer,
    dst: RootSegmentLayer,
    fx: f64,
) !RootSegmentLayerTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .wtrt1 = src.wtrt1 * fy,
            .wtrt1n = src.wtrt1n * fy,
            .wtrt1p = src.wtrt1p * fy,
            .wtrt2 = src.wtrt2 * fy,
            .wtrt2n = src.wtrt2n * fy,
            .wtrt2p = src.wtrt2p * fy,
            .rtlg1 = src.rtlg1 * fy,
            .rtlg2 = src.rtlg2 * fy,
            .rtn2 = src.rtn2 * fy,
        },
        .dst = .{
            .wtrt1 = dst.wtrt1 + fx * src.wtrt1,
            .wtrt1n = dst.wtrt1n + fx * src.wtrt1n,
            .wtrt1p = dst.wtrt1p + fx * src.wtrt1p,
            .wtrt2 = dst.wtrt2 + fx * src.wtrt2,
            .wtrt2n = dst.wtrt2n + fx * src.wtrt2n,
            .wtrt2p = dst.wtrt2p + fx * src.wtrt2p,
            .rtlg1 = dst.rtlg1 + fx * src.rtlg1,
            .rtlg2 = dst.rtlg2 + fx * src.rtlg2,
            .rtn2 = dst.rtn2 + fx * src.rtn2,
        },
    };
}

// ─── root biomass (lines 8849–8880 dst / 9263–9279 src; per N, per layer, NOT per NR) ─

pub const RootBiomassLayer = struct {
    cpoolr: f64, // CPOOLR (labile C pool)
    zpoolr: f64, // ZPOOLR (labile N pool)
    ppoolr: f64, // PPOOLR (labile P pool)
    wtrtl: f64, // WTRTL  (live root weight)
    wtrtd: f64, // WTRTD  (dead root weight)
    wsrtl: f64, // WSRTL  (root shoot weight)
    rtn1: f64, // RTN1
    rtnl: f64, // RTNL   (root N live)
    rtlgp: f64, // RTLGP  (root length for P uptake)
    rtdnp: f64, // RTDNP  (root density for N/P)
    rtvlp: f64, // RTVLP  (root volume for P uptake)
    rtvlw: f64, // RTVLW  (root volume for water)
    rrad1: f64, // RRAD1  (root radius 1)
    rrad2: f64, // RRAD2  (root radius 2)
    rtarp: f64, // RTARP  (root area for P uptake)
    rtlga: f64, // RTLGA  (root length accumulated)
};

pub const RootBiomassLayerTransferResult = struct { src: RootBiomassLayer, dst: RootBiomassLayer };

pub fn transferRootBiomassLayer(
    src: RootBiomassLayer,
    dst: RootBiomassLayer,
    fx: f64,
) !RootBiomassLayerTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .cpoolr = src.cpoolr * fy,
            .zpoolr = src.zpoolr * fy,
            .ppoolr = src.ppoolr * fy,
            .wtrtl = src.wtrtl * fy,
            .wtrtd = src.wtrtd * fy,
            .wsrtl = src.wsrtl * fy,
            .rtn1 = src.rtn1 * fy,
            .rtnl = src.rtnl * fy,
            .rtlgp = src.rtlgp * fy,
            .rtdnp = src.rtdnp * fy,
            .rtvlp = src.rtvlp * fy,
            .rtvlw = src.rtvlw * fy,
            .rrad1 = src.rrad1 * fy,
            .rrad2 = src.rrad2 * fy,
            .rtarp = src.rtarp * fy,
            .rtlga = src.rtlga * fy,
        },
        .dst = .{
            .cpoolr = dst.cpoolr + fx * src.cpoolr,
            .zpoolr = dst.zpoolr + fx * src.zpoolr,
            .ppoolr = dst.ppoolr + fx * src.ppoolr,
            .wtrtl = dst.wtrtl + fx * src.wtrtl,
            .wtrtd = dst.wtrtd + fx * src.wtrtd,
            .wsrtl = dst.wsrtl + fx * src.wsrtl,
            .rtn1 = dst.rtn1 + fx * src.rtn1,
            .rtnl = dst.rtnl + fx * src.rtnl,
            .rtlgp = dst.rtlgp + fx * src.rtlgp,
            .rtdnp = dst.rtdnp + fx * src.rtdnp,
            .rtvlp = dst.rtvlp + fx * src.rtvlp,
            .rtvlw = dst.rtvlw + fx * src.rtvlw,
            .rrad1 = dst.rrad1 + fx * src.rrad1,
            .rrad2 = dst.rrad2 + fx * src.rrad2,
            .rtarp = dst.rtarp + fx * src.rtarp,
            .rtlga = dst.rtlga + fx * src.rtlga,
        },
    };
}

// ─── root nodules (lines 8885–8896 dst / 9284–9289 src; per NZ, per layer) ──
// Applied after all N/NR loops; gate is the same WTRTL pair as above.

pub const RootNodulesLayer = struct {
    wtndl: f64, // WTNDL  (nodule C)
    wtndln: f64, // WTNDLN (nodule N)
    wtndlp: f64, // WTNDLP (nodule P)
    cpooln: f64, // CPOOLN (labile C)
    zpooln: f64, // ZPOOLN (labile N)
    ppooln: f64, // PPOOLN (labile P)
};

pub const RootNodulesLayerTransferResult = struct { src: RootNodulesLayer, dst: RootNodulesLayer };

pub fn transferRootNodulesLayer(
    src: RootNodulesLayer,
    dst: RootNodulesLayer,
    fx: f64,
) !RootNodulesLayerTransferResult {
    if (fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    const fy = 1.0 - fx;
    return .{
        .src = .{
            .wtndl = src.wtndl * fy,
            .wtndln = src.wtndln * fy,
            .wtndlp = src.wtndlp * fy,
            .cpooln = src.cpooln * fy,
            .zpooln = src.zpooln * fy,
            .ppooln = src.ppooln * fy,
        },
        .dst = .{
            .wtndl = dst.wtndl + fx * src.wtndl,
            .wtndln = dst.wtndln + fx * src.wtndln,
            .wtndlp = dst.wtndlp + fx * src.wtndlp,
            .cpooln = dst.cpooln + fx * src.cpooln,
            .zpooln = dst.zpooln + fx * src.zpooln,
            .ppooln = dst.ppooln + fx * src.ppooln,
        },
    };
}

// ─── tests ───────────────────────────────────────────────────────────────────

test "transferRootGasesLayer mass conservation" {
    const src = RootGasesLayer{
        .co2a = 5.0,
        .oxya = 2.0,
        .ch4a = 0.5,
        .z2oa = 0.1,
        .zh3a = 0.2,
        .h2ga = 0.05,
        .co2p = 3.0,
        .oxyp = 1.0,
        .ch4p = 0.3,
        .z2op = 0.0,
        .zh3p = 0.1,
        .h2gp = 0.01,
    };
    const dst = RootGasesLayer{
        .co2a = 1.0,
        .oxya = 0.5,
        .ch4a = 0.0,
        .z2oa = 0.0,
        .zh3a = 0.0,
        .h2ga = 0.0,
        .co2p = 0.5,
        .oxyp = 0.0,
        .ch4p = 0.0,
        .z2op = 0.0,
        .zh3p = 0.0,
        .h2gp = 0.0,
    };
    const r = try transferRootGasesLayer(src, dst, 0.4);
    try std.testing.expectApproxEqRel(r.src.co2a + r.dst.co2a, src.co2a + dst.co2a, 1e-12);
    try std.testing.expectApproxEqRel(r.src.h2gp + r.dst.h2gp, src.h2gp + dst.h2gp, 1e-12);
}

test "transferRootGasesLayer zero fraction identity" {
    const src = RootGasesLayer{
        .co2a = 5.0,
        .oxya = 0,
        .ch4a = 0,
        .z2oa = 0,
        .zh3a = 0,
        .h2ga = 0,
        .co2p = 0,
        .oxyp = 0,
        .ch4p = 0,
        .z2op = 0,
        .zh3p = 0,
        .h2gp = 0,
    };
    const dst = std.mem.zeroes(RootGasesLayer);
    const r = try transferRootGasesLayer(src, dst, 0.0);
    try std.testing.expectEqual(r.src.co2a, 5.0);
    try std.testing.expectEqual(r.dst.co2a, 0.0);
}

test "transferRootGasesLayer invalid fraction" {
    const z = std.mem.zeroes(RootGasesLayer);
    try std.testing.expectError(error.InvalidFraction, transferRootGasesLayer(z, z, 1.5));
    try std.testing.expectError(error.InvalidFraction, transferRootGasesLayer(z, z, -0.1));
}

test "transferRootSegmentLayer mass conservation" {
    const src = RootSegmentLayer{
        .wtrt1 = 3.0,
        .wtrt1n = 0.1,
        .wtrt1p = 0.05,
        .wtrt2 = 5.0,
        .wtrt2n = 0.2,
        .wtrt2p = 0.1,
        .rtlg1 = 10.0,
        .rtlg2 = 8.0,
        .rtn2 = 0.15,
    };
    const dst = RootSegmentLayer{
        .wtrt1 = 1.0,
        .wtrt1n = 0.0,
        .wtrt1p = 0.0,
        .wtrt2 = 2.0,
        .wtrt2n = 0.0,
        .wtrt2p = 0.0,
        .rtlg1 = 3.0,
        .rtlg2 = 2.0,
        .rtn2 = 0.0,
    };
    const r = try transferRootSegmentLayer(src, dst, 0.5);
    try std.testing.expectApproxEqRel(r.src.wtrt1 + r.dst.wtrt1, src.wtrt1 + dst.wtrt1, 1e-12);
    try std.testing.expectApproxEqRel(r.src.rtlg2 + r.dst.rtlg2, src.rtlg2 + dst.rtlg2, 1e-12);
}

test "transferRootBiomassLayer mass conservation" {
    var src = std.mem.zeroes(RootBiomassLayer);
    var dst = std.mem.zeroes(RootBiomassLayer);
    src.wtrtl = 8.0;
    dst.wtrtl = 2.0;
    src.rtlga = 50.0;
    dst.rtlga = 10.0;
    const r = try transferRootBiomassLayer(src, dst, 0.6);
    try std.testing.expectApproxEqRel(r.src.wtrtl + r.dst.wtrtl, src.wtrtl + dst.wtrtl, 1e-12);
    try std.testing.expectApproxEqRel(r.src.rtlga + r.dst.rtlga, src.rtlga + dst.rtlga, 1e-12);
}

test "transferRootBiomassLayer full transfer" {
    var src = std.mem.zeroes(RootBiomassLayer);
    src.cpoolr = 4.0;
    const dst = std.mem.zeroes(RootBiomassLayer);
    const r = try transferRootBiomassLayer(src, dst, 1.0);
    try std.testing.expectApproxEqRel(r.src.cpoolr, 0.0, 1e-12);
    try std.testing.expectApproxEqRel(r.dst.cpoolr, 4.0, 1e-12);
}

test "transferRootNodulesLayer mass conservation" {
    const src = RootNodulesLayer{
        .wtndl = 2.0,
        .wtndln = 0.1,
        .wtndlp = 0.05,
        .cpooln = 1.0,
        .zpooln = 0.05,
        .ppooln = 0.02,
    };
    const dst = RootNodulesLayer{
        .wtndl = 0.5,
        .wtndln = 0.0,
        .wtndlp = 0.0,
        .cpooln = 0.2,
        .zpooln = 0.0,
        .ppooln = 0.0,
    };
    const r = try transferRootNodulesLayer(src, dst, 0.3);
    try std.testing.expectApproxEqRel(r.src.wtndl + r.dst.wtndl, src.wtndl + dst.wtndl, 1e-12);
    try std.testing.expectApproxEqRel(r.src.ppooln + r.dst.ppooln, src.ppooln + dst.ppooln, 1e-12);
}

test "transferRootNodulesLayer invalid fraction" {
    const z = std.mem.zeroes(RootNodulesLayer);
    try std.testing.expectError(error.InvalidFraction, transferRootNodulesLayer(z, z, 2.0));
}
