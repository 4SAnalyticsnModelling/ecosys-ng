const std = @import("std");

/// Organic matter C/N/P state grouped by classification.
/// K=0-5 litter types, NO=1-7 organism groups, M=1-3 size classes.
pub const MicrobialPools = struct {
    /// OMC[K=0-5][NO=1-7][M=1-3] (g C).
    omc: [6][7][3]f64,
    /// OMN[K=0-5][NO=1-7][M=1-3] (g N).
    omn: [6][7][3]f64,
    /// OMP[K=0-5][NO=1-7][M=1-3] (g P).
    omp: [6][7][3]f64,
};

/// K=0-4 litter types, M=1-2 size classes.
pub const ResiduePools = struct {
    /// ORC[K=0-4][M=1-2] (g C).
    orc: [5][2]f64,
    /// ORN[K=0-4][M=1-2] (g N).
    orn: [5][2]f64,
    /// ORP[K=0-4][M=1-2] (g P).
    orp: [5][2]f64,
};

/// K=0-4 litter types (scalar per K).
pub const AdsorptionPools = struct {
    /// OHC[K=0-4] (g C).
    ohc: [5]f64,
    /// OHN[K=0-4] (g N).
    ohn: [5]f64,
    /// OHP[K=0-4] (g P).
    ohp: [5]f64,
    /// OHA[K=0-4] (g acetate).
    oha: [5]f64,
};

/// K=0-4 litter types, M=1-5 size classes.
pub const SocPools = struct {
    /// OSC[K=0-4][M=1-5] (g C).
    osc: [5][5]f64,
    /// OSA[K=0-4][M=1-5] (g colonised C).
    osa: [5][5]f64,
    /// OSN[K=0-4][M=1-5] (g N).
    osn: [5][5]f64,
    /// OSP[K=0-4][M=1-5] (g P).
    osp: [5][5]f64,
};

/// Erosion fluxes; same dimensions as pools.
pub const MicrobialFluxes = struct {
    tomcer: [6][7][3]f64,
    tomner: [6][7][3]f64,
    tomper: [6][7][3]f64,
};
pub const ResidueFluxes = struct {
    torcer: [5][2]f64,
    torner: [5][2]f64,
    torper: [5][2]f64,
};
pub const AdsorptionFluxes = struct {
    tohcer: [5]f64,
    tohner: [5]f64,
    tohper: [5]f64,
    tohaer: [5]f64,
};
pub const SocFluxes = struct {
    toscer: [5][5]f64,
    tosaer: [5][5]f64,
    tosner: [5][5]f64,
    tosper: [5][5]f64,
};

pub const Result = struct {
    microbial: MicrobialPools,
    residue: ResiduePools,
    adsorption: AdsorptionPools,
    soc: SocPools,
    /// DORGE accumulator: total C eroded (g).
    dorge: f64,
    /// DORGP accumulator: total P eroded (g).
    dorgp: f64,
};

/// Direct translation of REDIST lines 5300--5342 (IERSNG-gated organic block).
///
/// Caller must preserve the enclosing gate: IERSNG is 1 or 3 and
/// ABS(TSEDER) > ZEROS. Mode 2 does not execute this erosion block.
pub fn apply(
    microbial: MicrobialPools,
    residue: ResiduePools,
    adsorption: AdsorptionPools,
    soc: SocPools,
    mf: MicrobialFluxes,
    rf: ResidueFluxes,
    af: AdsorptionFluxes,
    sf: SocFluxes,
) !Result {
    var dorge: f64 = 0.0;
    var dorgp: f64 = 0.0;

    var new_omc: [6][7][3]f64 = undefined;
    var new_omn: [6][7][3]f64 = undefined;
    var new_omp: [6][7][3]f64 = undefined;

    // REDIST lines 5301-5312: DO K=0,5 / DO NO=1,7 / DO M=1,3
    for (0..6) |k| {
        for (0..7) |no| {
            for (0..3) |m| {
                if (!std.math.isFinite(microbial.omc[k][no][m])) return error.InvalidMicrobialPool;
                if (!std.math.isFinite(microbial.omn[k][no][m]) or
                    !std.math.isFinite(microbial.omp[k][no][m])) return error.InvalidMicrobialPool;
                if (!std.math.isFinite(mf.tomcer[k][no][m]) or
                    !std.math.isFinite(mf.tomner[k][no][m]) or
                    !std.math.isFinite(mf.tomper[k][no][m])) return error.InvalidMicrobialFlux;
                new_omc[k][no][m] = microbial.omc[k][no][m] + mf.tomcer[k][no][m];
                new_omn[k][no][m] = microbial.omn[k][no][m] + mf.tomner[k][no][m];
                new_omp[k][no][m] = microbial.omp[k][no][m] + mf.tomper[k][no][m];
                if (!std.math.isFinite(new_omc[k][no][m]) or
                    !std.math.isFinite(new_omn[k][no][m]) or
                    !std.math.isFinite(new_omp[k][no][m])) return error.NonFiniteMicrobialPool;
                dorge += mf.tomcer[k][no][m];
                dorgp += mf.tomper[k][no][m];
            }
        }
    }

    var new_orc: [5][2]f64 = undefined;
    var new_orn: [5][2]f64 = undefined;
    var new_orp: [5][2]f64 = undefined;
    var new_ohc: [5]f64 = undefined;
    var new_ohn: [5]f64 = undefined;
    var new_ohp: [5]f64 = undefined;
    var new_oha: [5]f64 = undefined;
    var new_osc: [5][5]f64 = undefined;
    var new_osa: [5][5]f64 = undefined;
    var new_osn: [5][5]f64 = undefined;
    var new_osp: [5][5]f64 = undefined;

    // REDIST lines 5313-5342: DO K=0,4
    for (0..5) |k| {
        // DO M=1,2: ORC/ORN/ORP
        for (0..2) |m| {
            if (!std.math.isFinite(residue.orc[k][m])) return error.InvalidResiduePool;
            if (!std.math.isFinite(residue.orn[k][m]) or
                !std.math.isFinite(residue.orp[k][m])) return error.InvalidResiduePool;
            if (!std.math.isFinite(rf.torcer[k][m]) or
                !std.math.isFinite(rf.torner[k][m]) or
                !std.math.isFinite(rf.torper[k][m])) return error.InvalidResidueFlux;
            new_orc[k][m] = residue.orc[k][m] + rf.torcer[k][m];
            new_orn[k][m] = residue.orn[k][m] + rf.torner[k][m];
            new_orp[k][m] = residue.orp[k][m] + rf.torper[k][m];
            if (!std.math.isFinite(new_orc[k][m]) or
                !std.math.isFinite(new_orn[k][m]) or
                !std.math.isFinite(new_orp[k][m])) return error.NonFiniteResiduePool;
            dorge += rf.torcer[k][m];
            dorgp += rf.torper[k][m];
        }
        // OHC/OHN/OHP/OHA scalar per K
        if (!std.math.isFinite(adsorption.ohc[k])) return error.InvalidAdsorptionPool;
        if (!std.math.isFinite(adsorption.ohn[k]) or
            !std.math.isFinite(adsorption.ohp[k]) or
            !std.math.isFinite(adsorption.oha[k])) return error.InvalidAdsorptionPool;
        if (!std.math.isFinite(af.tohcer[k]) or
            !std.math.isFinite(af.tohner[k]) or
            !std.math.isFinite(af.tohper[k]) or
            !std.math.isFinite(af.tohaer[k])) return error.InvalidAdsorptionFlux;
        new_ohc[k] = adsorption.ohc[k] + af.tohcer[k];
        new_ohn[k] = adsorption.ohn[k] + af.tohner[k];
        new_ohp[k] = adsorption.ohp[k] + af.tohper[k];
        new_oha[k] = adsorption.oha[k] + af.tohaer[k];
        if (!std.math.isFinite(new_ohc[k]) or !std.math.isFinite(new_ohn[k]) or
            !std.math.isFinite(new_ohp[k]) or !std.math.isFinite(new_oha[k]))
            return error.NonFiniteAdsorptionPool;
        dorge += af.tohcer[k] + af.tohaer[k];
        dorgp += af.tohper[k];
        // DO M=1,5: OSC/OSA/OSN/OSP
        for (0..5) |m| {
            if (!std.math.isFinite(soc.osc[k][m])) return error.InvalidSocPool;
            if (!std.math.isFinite(soc.osa[k][m]) or
                !std.math.isFinite(soc.osn[k][m]) or
                !std.math.isFinite(soc.osp[k][m])) return error.InvalidSocPool;
            if (!std.math.isFinite(sf.toscer[k][m]) or
                !std.math.isFinite(sf.tosaer[k][m]) or
                !std.math.isFinite(sf.tosner[k][m]) or
                !std.math.isFinite(sf.tosper[k][m])) return error.InvalidSocFlux;
            new_osc[k][m] = soc.osc[k][m] + sf.toscer[k][m];
            new_osa[k][m] = soc.osa[k][m] + sf.tosaer[k][m];
            new_osn[k][m] = soc.osn[k][m] + sf.tosner[k][m];
            new_osp[k][m] = soc.osp[k][m] + sf.tosper[k][m];
            if (!std.math.isFinite(new_osc[k][m]) or !std.math.isFinite(new_osa[k][m]) or
                !std.math.isFinite(new_osn[k][m]) or !std.math.isFinite(new_osp[k][m]))
                return error.NonFiniteSocPool;
            dorge += sf.toscer[k][m];
            dorgp += sf.tosper[k][m];
        }
    }

    if (!std.math.isFinite(dorge) or !std.math.isFinite(dorgp))
        return error.NonFiniteOrganicAccumulator;

    return Result{
        .microbial = MicrobialPools{
            .omc = new_omc,
            .omn = new_omn,
            .omp = new_omp,
        },
        .residue = ResiduePools{
            .orc = new_orc,
            .orn = new_orn,
            .orp = new_orp,
        },
        .adsorption = AdsorptionPools{
            .ohc = new_ohc,
            .ohn = new_ohn,
            .ohp = new_ohp,
            .oha = new_oha,
        },
        .soc = SocPools{
            .osc = new_osc,
            .osa = new_osa,
            .osn = new_osn,
            .osp = new_osp,
        },
        .dorge = dorge,
        .dorgp = dorgp,
    };
}

test "REDIST erosion organic OMC accumulates across all K,NO,M dimensions" {
    var mf = std.mem.zeroes(MicrobialFluxes);
    mf.tomcer[0][0][0] = 1.0;
    mf.tomcer[5][6][2] = 2.0;
    const result = try apply(
        std.mem.zeroes(MicrobialPools),
        std.mem.zeroes(ResiduePools),
        std.mem.zeroes(AdsorptionPools),
        std.mem.zeroes(SocPools),
        mf,
        std.mem.zeroes(ResidueFluxes),
        std.mem.zeroes(AdsorptionFluxes),
        std.mem.zeroes(SocFluxes),
    );
    try std.testing.expectApproxEqRel(@as(f64, 1.0), result.microbial.omc[0][0][0], 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.microbial.omc[5][6][2], 1.0e-15);
}

test "REDIST erosion organic DORGE accumulates OMC + OHC + OHA + OSC" {
    var af = std.mem.zeroes(AdsorptionFluxes);
    af.tohcer[2] = 3.0;
    af.tohaer[2] = 1.0;
    const result = try apply(
        std.mem.zeroes(MicrobialPools),
        std.mem.zeroes(ResiduePools),
        std.mem.zeroes(AdsorptionPools),
        std.mem.zeroes(SocPools),
        std.mem.zeroes(MicrobialFluxes),
        std.mem.zeroes(ResidueFluxes),
        af,
        std.mem.zeroes(SocFluxes),
    );
    // DORGE includes both OHC and OHA
    try std.testing.expectApproxEqRel(@as(f64, 4.0), result.dorge, 1.0e-15);
}

test "REDIST erosion organic ORC uses K=0-4 only" {
    var rf = std.mem.zeroes(ResidueFluxes);
    rf.torcer[4][1] = 5.0;
    const result = try apply(
        std.mem.zeroes(MicrobialPools),
        std.mem.zeroes(ResiduePools),
        std.mem.zeroes(AdsorptionPools),
        std.mem.zeroes(SocPools),
        std.mem.zeroes(MicrobialFluxes),
        rf,
        std.mem.zeroes(AdsorptionFluxes),
        std.mem.zeroes(SocFluxes),
    );
    try std.testing.expectApproxEqRel(@as(f64, 5.0), result.residue.orc[4][1], 1.0e-15);
}

test "REDIST erosion organic rejects non-finite microbial pool" {
    var bad = std.mem.zeroes(MicrobialPools);
    bad.omc[1][2][0] = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidMicrobialPool,
        apply(
            bad,
            std.mem.zeroes(ResiduePools),
            std.mem.zeroes(AdsorptionPools),
            std.mem.zeroes(SocPools),
            std.mem.zeroes(MicrobialFluxes),
            std.mem.zeroes(ResidueFluxes),
            std.mem.zeroes(AdsorptionFluxes),
            std.mem.zeroes(SocFluxes),
        ),
    );
}

test "REDIST erosion organic maps nitrogen phosphorus acetate and SOC fields" {
    var mf = std.mem.zeroes(MicrobialFluxes);
    mf.tomner[5][6][2] = 1.0;
    mf.tomper[5][6][2] = 2.0;
    var rf = std.mem.zeroes(ResidueFluxes);
    rf.torner[4][1] = 3.0;
    rf.torper[4][1] = 4.0;
    var af = std.mem.zeroes(AdsorptionFluxes);
    af.tohner[4] = 5.0;
    af.tohper[4] = 6.0;
    af.tohaer[4] = 7.0;
    var sf = std.mem.zeroes(SocFluxes);
    sf.tosaer[4][4] = 8.0;
    sf.tosner[4][4] = 9.0;
    sf.tosper[4][4] = 10.0;
    const result = try apply(
        std.mem.zeroes(MicrobialPools),
        std.mem.zeroes(ResiduePools),
        std.mem.zeroes(AdsorptionPools),
        std.mem.zeroes(SocPools),
        mf,
        rf,
        af,
        sf,
    );
    try std.testing.expectEqual(@as(f64, 1.0), result.microbial.omn[5][6][2]);
    try std.testing.expectEqual(@as(f64, 2.0), result.microbial.omp[5][6][2]);
    try std.testing.expectEqual(@as(f64, 3.0), result.residue.orn[4][1]);
    try std.testing.expectEqual(@as(f64, 4.0), result.residue.orp[4][1]);
    try std.testing.expectEqual(@as(f64, 7.0), result.adsorption.oha[4]);
    try std.testing.expectEqual(@as(f64, 8.0), result.soc.osa[4][4]);
    try std.testing.expectEqual(@as(f64, 9.0), result.soc.osn[4][4]);
    try std.testing.expectEqual(@as(f64, 10.0), result.soc.osp[4][4]);
    try std.testing.expectEqual(@as(f64, 22.0), result.dorgp);
}

test "REDIST erosion organic rejects non-carbon invalid flux" {
    var sf = std.mem.zeroes(SocFluxes);
    sf.tosner[2][3] = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSocFlux,
        apply(
            std.mem.zeroes(MicrobialPools),
            std.mem.zeroes(ResiduePools),
            std.mem.zeroes(AdsorptionPools),
            std.mem.zeroes(SocPools),
            std.mem.zeroes(MicrobialFluxes),
            std.mem.zeroes(ResidueFluxes),
            std.mem.zeroes(AdsorptionFluxes),
            sf,
        ),
    );
}

test "REDIST erosion organic rejects non-carbon pool overflow" {
    var microbial = std.mem.zeroes(MicrobialPools);
    microbial.omn[0][0][0] = std.math.floatMax(f64);
    var mf = std.mem.zeroes(MicrobialFluxes);
    mf.tomner[0][0][0] = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteMicrobialPool,
        apply(
            microbial,
            std.mem.zeroes(ResiduePools),
            std.mem.zeroes(AdsorptionPools),
            std.mem.zeroes(SocPools),
            mf,
            std.mem.zeroes(ResidueFluxes),
            std.mem.zeroes(AdsorptionFluxes),
            std.mem.zeroes(SocFluxes),
        ),
    );
}
