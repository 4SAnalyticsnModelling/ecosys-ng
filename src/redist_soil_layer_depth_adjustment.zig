const std = @import("std");

/// Per-layer freeze-thaw depth increment (m step⁻¹).
/// Translates REDIST lines 7961–7962:
///   DDLYXF = DVOLI(LX) * DENSJ / (FMPR(LX) * AREA)
/// DENSJ = 1 − DENSI, the density correction for ice volume expansion.
pub fn freezethawIncrement(dvoli: f64, densj: f64, fmpr: f64, area: f64) !f64 {
    if (!std.math.isFinite(dvoli)) return error.InvalidDvoli;
    if (!std.math.isFinite(densj)) return error.InvalidDensj;
    if (!std.math.isFinite(fmpr) or fmpr <= 0.0) return error.InvalidFmpr;
    if (!std.math.isFinite(area) or area <= 0.0) return error.InvalidArea;
    return dvoli * densj / (fmpr * area);
}

/// Per-layer SOC-change depth increment (m step⁻¹).
/// Translates REDIST lines 7898–7899:
///   DDLYXC = 1.82E-06 * DORGC / AREA / ((1 - FHOL) * BKDSI)
/// 1.82E-06 converts g C m⁻² Mg⁻¹ m⁻³ → m of depth per unit SOC change.
pub fn socIncrement(dorgc: f64, area: f64, fhol: f64, bkdsi: f64) !f64 {
    if (!std.math.isFinite(dorgc)) return error.InvalidDorgc;
    if (!std.math.isFinite(area) or area <= 0.0) return error.InvalidArea;
    if (!std.math.isFinite(fhol) or fhol < 0.0 or fhol >= 1.0) return error.InvalidFhol;
    if (!std.math.isFinite(bkdsi) or bkdsi <= 0.0) return error.InvalidBkdsi;
    const denom = (1.0 - fhol) * bkdsi;
    return 1.82e-6 * dorgc / area / denom;
}

/// Column erosion depth increment (m step⁻¹), uniform across all layers.
/// Translates REDIST lines 7845–7856:
///   At NL (bottom mineral layer, if BKDS(NU)>0):
///     DDLYXE  = TSEDER / (AREA * BKVLNU / VOLX_NU)
///   At NUL (litter–soil boundary layer):
///     DDLYXE += TSEDSK / (AREA * BKDSI_NUL)
/// The result is the same for all layers: DDLYX(LX,3) = DDLYR(LX,3) = DDLYXE.
/// Callers must zero DDLYX(:,3)/DDLYR(:,3) when the erosion gate is inactive.
pub fn erosionIncrement(
    tseder: f64,
    tsedsk: f64,
    area: f64,
    bkvlnu_over_volx_nu: f64,
    bkdsi_nul: f64,
    nu_is_soil: bool,
) !f64 {
    if (!std.math.isFinite(tseder)) return error.InvalidTseder;
    if (!std.math.isFinite(tsedsk)) return error.InvalidTsedsk;
    if (!std.math.isFinite(area) or area <= 0.0) return error.InvalidArea;
    var ddlyxe: f64 = 0.0;
    if (nu_is_soil) {
        if (!std.math.isFinite(bkvlnu_over_volx_nu) or bkvlnu_over_volx_nu <= 0.0)
            return error.InvalidBkvlnuRatio;
        ddlyxe = tseder / (area * bkvlnu_over_volx_nu);
    }
    if (!std.math.isFinite(bkdsi_nul) or bkdsi_nul <= 0.0) return error.InvalidBkdsiNul;
    ddlyxe += tsedsk / (area * bkdsi_nul);
    return ddlyxe;
}

/// Per-layer pond depth increment (m step⁻¹) when BKDS(LX) ≤ 0 and BKDS(LX+1) > 0
/// (pond layer immediately above soil, case IFLGL = 2).
/// Translates REDIST line 7752–7753:
///   DDLYXP = DLYR(LX) − (VOLW(LX) + VOLI(LX)) / (FMPR(LX) * AREA)
pub fn pondSoilBoundaryIncrement(
    dlyr: f64,
    volw: f64,
    voli: f64,
    fmpr: f64,
    area: f64,
) !f64 {
    if (!std.math.isFinite(dlyr)) return error.InvalidDlyr;
    if (!std.math.isFinite(volw)) return error.InvalidVolw;
    if (!std.math.isFinite(voli)) return error.InvalidVoli;
    if (!std.math.isFinite(fmpr) or fmpr <= 0.0) return error.InvalidFmpr;
    if (!std.math.isFinite(area) or area <= 0.0) return error.InvalidArea;
    return dlyr - (volw + voli) / (fmpr * area);
}

/// Per-layer pond depth increment (m step⁻¹) when BKDS(LX) ≤ 0 and BKDS(LX+1) ≤ 0
/// (pond on pond, uses initial layer depth DLYRI, REDIST lines 7758–7759):
///   DDLYXP = DLYRI(LX) − (VOLW(LX) + VOLI(LX)) / (FMPR(LX) * AREA)
pub fn pondPondIncrement(
    dlyri: f64,
    volw: f64,
    voli: f64,
    fmpr: f64,
    area: f64,
) !f64 {
    if (!std.math.isFinite(dlyri)) return error.InvalidDlyri;
    if (!std.math.isFinite(volw)) return error.InvalidVolw;
    if (!std.math.isFinite(voli)) return error.InvalidVoli;
    if (!std.math.isFinite(fmpr) or fmpr <= 0.0) return error.InvalidFmpr;
    if (!std.math.isFinite(area) or area <= 0.0) return error.InvalidArea;
    return dlyri - (volw + voli) / (fmpr * area);
}

/// Prefix-sum accumulation of per-layer freeze-thaw DDLYX/DDLYR arrays.
/// Translates REDIST lines 7964–7995 for soil layers only (BKDS(LX) > 0).
/// `ddlyxf_per_layer`: per-layer freeze-thaw increment (already computed by
///   `freezethawIncrement`), length == number of soil layers from NL down to NU.
///   Index 0 = bottom layer (NL), index len-1 = top layer (NU).
/// `ddlyx`, `ddlyr`: output slices of the same length; filled bottom-to-top.
/// After this call the caller is responsible for propagating the surface value
/// one layer above NU (REDIST lines 7973–7976) if needed.
pub fn accumulateFreezethawAdjustments(
    ddlyxf_per_layer: []const f64,
    ddlyx: []f64,
    ddlyr: []f64,
) !void {
    if (ddlyx.len != ddlyxf_per_layer.len) return error.LengthMismatch;
    if (ddlyr.len != ddlyxf_per_layer.len) return error.LengthMismatch;
    const n = ddlyxf_per_layer.len;
    for (0..n) |i| {
        const ddlyxf = ddlyxf_per_layer[i];
        if (!std.math.isFinite(ddlyxf)) return error.NonFiniteIncrement;
        if (i == 0) {
            // bottom layer: NL case — DDLYX(NL,2) = DDLYXF, DDLYR(NL,2) = 0
            ddlyx[0] = ddlyxf;
            ddlyr[0] = 0.0;
        } else {
            // layers above bottom: DDLYX(LX,2) = DDLYXF + DDLYX(LX+1,2)
            //                      DDLYR(LX,2) = DDLYX(LX+1,2)
            ddlyx[i] = ddlyxf + ddlyx[i - 1];
            ddlyr[i] = ddlyx[i - 1];
        }
    }
}

/// Prefix-sum accumulation of per-layer SOC DDLYX/DDLYR arrays.
/// Translates REDIST lines 7896–7948 for soil layers only (BKDS(LX) > 0).
/// Inputs (all bottom-to-top, index 0 = bottom layer NL):
///   `ddlyxc_per_layer`: per-layer SOC increment from `socIncrement`.
///   `is_segment_bottom`: true for layer LX when at NL, or when the layer above
///     (LX+1) is pond (BKDS(LX+1)≤0), or at organic/mineral boundary
///     (CORGC(LX)≥FORGC and CORGC(LX+1)<FORGC). Index 0 always set to true.
///   `dlyr`, `dlyri`: current and initial layer depth (m), same indexing.
/// `ddlyx`, `ddlyr_out`: output slices of the same length; filled bottom-to-top.
/// Propagation of the surface value above NU (REDIST lines 7913–7917) is the
/// caller's responsibility.
pub fn accumulateSocAdjustments(
    ddlyxc_per_layer: []const f64,
    is_segment_bottom: []const bool,
    dlyr: []const f64,
    dlyri: []const f64,
    ddlyx: []f64,
    ddlyr_out: []f64,
) !void {
    const n = ddlyxc_per_layer.len;
    if (is_segment_bottom.len != n) return error.LengthMismatch;
    if (dlyr.len != n) return error.LengthMismatch;
    if (dlyri.len != n) return error.LengthMismatch;
    if (ddlyx.len != n) return error.LengthMismatch;
    if (ddlyr_out.len != n) return error.LengthMismatch;
    for (0..n) |i| {
        const ddlyxc = ddlyxc_per_layer[i];
        if (!std.math.isFinite(ddlyxc)) return error.NonFiniteIncrement;
        if (is_segment_bottom[i]) {
            // NL or organic/mineral or pond boundary: DDLYX=DDLYXC, DDLYR=0
            ddlyx[i] = ddlyxc;
            ddlyr_out[i] = 0.0;
        } else {
            // Middle of continuous mineral segment (prefix sum):
            //   DDLYX(LX,4) = DDLYXC + DDLYX(LX+1,4)
            //   DDLYR(LX,4) = DDLYX(LX+1,4) + (DLYRI(LX) - DLYR(LX))
            ddlyx[i] = ddlyxc + ddlyx[i - 1];
            ddlyr_out[i] = ddlyx[i - 1] + (dlyri[i] - dlyr[i]);
        }
    }
}

test "freezethaw increment: expansion when ice forms" {
    // DVOLI=0.001 m3, DENSJ=0.09 (ice less dense than water), FMPR=0.9, AREA=1.0 m2
    const inc = try freezethawIncrement(0.001, 0.09, 0.9, 1.0);
    try std.testing.expectApproxEqRel(@as(f64, 0.001 * 0.09 / (0.9 * 1.0)), inc, 1.0e-14);
}

test "freezethaw increment: contraction when ice melts" {
    const inc = try freezethawIncrement(-0.002, 0.09, 0.8, 2.0);
    try std.testing.expect(inc < 0.0);
}

test "freezethaw increment: rejects zero fmpr" {
    try std.testing.expectError(error.InvalidFmpr, freezethawIncrement(0.001, 0.09, 0.0, 1.0));
}

test "SOC increment: positive SOC gain deepens layer" {
    // DORGC=100 gC, AREA=1 m2, FHOL=0.05, BKDSI=1.2 Mg/m3
    const inc = try socIncrement(100.0, 1.0, 0.05, 1.2);
    try std.testing.expectApproxEqRel(
        @as(f64, 1.82e-6 * 100.0 / 1.0 / (0.95 * 1.2)),
        inc,
        1.0e-14,
    );
}

test "SOC increment: negative SOC loss shallows layer" {
    const inc = try socIncrement(-50.0, 1.0, 0.0, 1.3);
    try std.testing.expect(inc < 0.0);
}

test "SOC increment: rejects fhol >= 1.0" {
    try std.testing.expectError(error.InvalidFhol, socIncrement(10.0, 1.0, 1.0, 1.2));
}

test "erosion increment: no erosion when tseder=tsedsk=0" {
    const inc = try erosionIncrement(0.0, 0.0, 1.0, 1.5, 1.2, true);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), inc, 1.0e-15);
}

test "erosion increment: surface-only sediment (tsedsk only)" {
    const inc = try erosionIncrement(0.0, 0.01, 1.0, 1.5, 1.2, true);
    const expected = 0.01 / (1.0 * 1.2);
    try std.testing.expectApproxEqRel(expected, inc, 1.0e-14);
}

test "erosion increment: no subsoil contribution when NU is pond" {
    // nu_is_soil=false means the surface layer is pond — TSEDER not applied
    const inc = try erosionIncrement(10.0, 0.01, 1.0, 1.5, 1.2, false);
    const expected = 0.01 / (1.0 * 1.2);
    try std.testing.expectApproxEqRel(expected, inc, 1.0e-14);
}

test "pond soil-boundary increment: water deficit shallows layer" {
    // DLYR=0.1, (VOLW+VOLI)/(FMPR*AREA) = 0.08 → DDLYXP = 0.02
    const inc = try pondSoilBoundaryIncrement(0.1, 0.06, 0.02, 0.9, 1.0);
    try std.testing.expectApproxEqRel(@as(f64, 0.1 - 0.08 / 0.9), inc, 1.0e-14);
}

test "accumulate freeze-thaw: bottom layer has DDLYR=0" {
    const ddlyxf = [_]f64{ 0.005, 0.003, 0.002 }; // index 0 = bottom (NL)
    var ddlyx = [_]f64{ 0.0, 0.0, 0.0 };
    var ddlyr = [_]f64{ 0.0, 0.0, 0.0 };
    try accumulateFreezethawAdjustments(&ddlyxf, &ddlyx, &ddlyr);
    try std.testing.expectApproxEqRel(@as(f64, 0.005), ddlyx[0], 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), ddlyr[0], 1.0e-15);
}

test "accumulate freeze-thaw: prefix sum propagates upward" {
    const ddlyxf = [_]f64{ 0.005, 0.003, 0.002 };
    var ddlyx = [_]f64{ 0.0, 0.0, 0.0 };
    var ddlyr = [_]f64{ 0.0, 0.0, 0.0 };
    try accumulateFreezethawAdjustments(&ddlyxf, &ddlyx, &ddlyr);
    // layer 1: DDLYX = 0.003 + 0.005 = 0.008; DDLYR = 0.005
    try std.testing.expectApproxEqRel(@as(f64, 0.008), ddlyx[1], 1.0e-14);
    try std.testing.expectApproxEqRel(@as(f64, 0.005), ddlyr[1], 1.0e-14);
    // layer 2: DDLYX = 0.002 + 0.008 = 0.010; DDLYR = 0.008
    try std.testing.expectApproxEqRel(@as(f64, 0.010), ddlyx[2], 1.0e-14);
    try std.testing.expectApproxEqRel(@as(f64, 0.008), ddlyr[2], 1.0e-14);
}

test "accumulate SOC: segment bottom resets DDLYR to zero" {
    const ddlyxc = [_]f64{ 0.004, 0.003 };
    const is_bottom = [_]bool{ true, false };
    const dlyr = [_]f64{ 0.10, 0.12 };
    const dlyri = [_]f64{ 0.10, 0.15 };
    var ddlyx = [_]f64{ 0.0, 0.0 };
    var ddlyr_out = [_]f64{ 0.0, 0.0 };
    try accumulateSocAdjustments(&ddlyxc, &is_bottom, &dlyr, &dlyri, &ddlyx, &ddlyr_out);
    try std.testing.expectApproxEqRel(@as(f64, 0.004), ddlyx[0], 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), ddlyr_out[0], 1.0e-15);
    // layer 1 (above): DDLYX = 0.003 + 0.004 = 0.007
    //                  DDLYR = 0.004 + (0.15 - 0.12) = 0.034
    try std.testing.expectApproxEqRel(@as(f64, 0.007), ddlyx[1], 1.0e-14);
    try std.testing.expectApproxEqRel(@as(f64, 0.004 + (0.15 - 0.12)), ddlyr_out[1], 1.0e-14);
}

test "accumulate freeze-thaw: rejects non-finite increment" {
    const ddlyxf = [_]f64{ std.math.nan(f64), 0.003 };
    var ddlyx = [_]f64{ 0.0, 0.0 };
    var ddlyr = [_]f64{ 0.0, 0.0 };
    try std.testing.expectError(
        error.NonFiniteIncrement,
        accumulateFreezethawAdjustments(&ddlyxf, &ddlyx, &ddlyr),
    );
}
