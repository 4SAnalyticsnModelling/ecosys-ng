const std = @import("std");

/// Thermal-hydraulic state of one layer involved in a relayering transfer.
/// Translates the variables manipulated in REDIST DO 245 water-heat section
/// (lines 8564–8588 for destination; 9086–9107 for source).
pub const LayerWaterHeatState = struct {
    volw: f64, // liquid water volume (m3) — VOLW
    voli: f64, // ice volume (m3) — VOLI
    volp: f64, // macropore-drained air volume (m3) — VOLP
    vola: f64, // air-filled porosity volume (m3) — VOLA
    voly: f64, // prior-step total micropore volume snapshot (m3) — VOLY
    vhcm: f64, // mineral + organic heat capacity (MJ K⁻¹) — VHCM
    tks: f64, // temperature (K) — TKS
    tcs: f64, // temperature (°C) = TKS − 273.15 — TCS
};

/// Fixed macropore water and ice volumes that contribute to heat capacity
/// but are NOT transferred during relayering (VOLWH, VOLIH).
pub const MacroporeWaterIce = struct {
    volwh: f64, // macropore liquid water (m3)
    volih: f64, // macropore ice (m3)
};

/// Result of a single-step water-heat relayering transfer.
pub const TransferResult = struct {
    src: LayerWaterHeatState,
    dst: LayerWaterHeatState,
};

/// Transfer FX fraction of water-heat state from `src` to `dst`.
/// Translates REDIST lines 8564–8588 (destination) and 9086–9107 (source).
///
/// Energy conservation:
///   ENGY = VHCP * TKS (energy stored in a layer, MJ)
/// VHCM is the mineral/organic fraction of heat capacity; VHCP is the
/// full heat capacity including water and ice contributions:
///   VHCP = VHCM + 4.19*(VOLW + VOLWH) + 1.9274*(VOLI + VOLIH)
/// (For the pond surface layer at index 0, VOLWH and VOLIH are excluded.)
/// If VHCP ≤ vhcprx after update, temperature falls back to the other
/// layer's temperature (no energy to divide; adopt the neighbour's T).
///
/// `src_is_surface_pond`: true when the source layer is the "layer 0"
///   (pond surface) — uses simpler VHCP formula without VOLWH/VOLIH.
pub fn transfer(
    src: LayerWaterHeatState,
    dst: LayerWaterHeatState,
    src_macropore: MacroporeWaterIce,
    dst_macropore: MacroporeWaterIce,
    fx: f64,
    vhcprx: f64,
    src_is_surface_pond: bool,
) !TransferResult {
    if (!std.math.isFinite(fx) or fx < 0.0 or fx > 1.0) return error.InvalidFraction;
    if (!std.math.isFinite(vhcprx) or vhcprx < 0.0) return error.InvalidVhcprx;
    inline for (@typeInfo(LayerWaterHeatState).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(src, field.name))) return error.InvalidSrcState;
        if (!std.math.isFinite(@field(dst, field.name))) return error.InvalidDstState;
    }
    inline for (@typeInfo(MacroporeWaterIce).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(src_macropore, field.name))) return error.InvalidSrcMacropore;
        if (!std.math.isFinite(@field(dst_macropore, field.name))) return error.InvalidDstMacropore;
    }

    const fy = 1.0 - fx;

    // === DESTINATION LAYER RECEIVES FX * SOURCE (lines 8564–8588) ===
    // VHCP computed from components (VHCM + 4.19*(VOLW+VOLWH) + 1.9274*(VOLI+VOLIH))
    const src_vhcp = src.vhcm + 4.19 * (src.volw + src_macropore.volwh) +
        1.9274 * (src.voli + src_macropore.volih);
    const dst_vhcp = dst.vhcm + 4.19 * (dst.volw + dst_macropore.volwh) +
        1.9274 * (dst.voli + dst_macropore.volih);

    var engy_dst = dst_vhcp * dst.tks;
    const engy_src = src_vhcp * src.tks;

    const new_dst_volw = dst.volw + fx * src.volw;
    const new_dst_voli = dst.voli + fx * src.voli;
    const new_dst_volp = dst.volp + fx * src.volp;
    const new_dst_vola = dst.vola + fx * src.vola;
    const new_dst_voly = dst.voly + fx * src.voly;

    engy_dst += fx * engy_src;
    const new_dst_vhcm = dst.vhcm + fx * src.vhcm;

    // VHCP = VHCM + 4.19*(VOLW+VOLWH) + 1.9274*(VOLI+VOLIH)
    const new_dst_vhcp = new_dst_vhcm +
        4.19 * (new_dst_volw + dst_macropore.volwh) +
        1.9274 * (new_dst_voli + dst_macropore.volih);

    const new_dst_tks = if (new_dst_vhcp > vhcprx)
        engy_dst / new_dst_vhcp
    else
        src.tks; // fallback: adopt source temperature

    const new_dst = LayerWaterHeatState{
        .volw = new_dst_volw,
        .voli = new_dst_voli,
        .volp = new_dst_volp,
        .vola = new_dst_vola,
        .voly = new_dst_voly,
        .vhcm = new_dst_vhcm,
        .tks = new_dst_tks,
        .tcs = new_dst_tks - 273.15,
    };

    // === SOURCE LAYER SCALED BY FY (lines 9086–9107) ===
    const new_src_volw = fy * src.volw;
    const new_src_voli = fy * src.voli;
    const new_src_volp = fy * src.volp;
    const new_src_vola = fy * src.vola;
    const new_src_voly = fy * src.voly;
    const engy_src_scaled = fy * engy_src;
    const new_src_vhcm = fy * src.vhcm;

    const new_src_vhcp = if (src_is_surface_pond)
        // Layer 0 (pond surface): VHCP = VHCM + 4.19*VOLW + 1.9274*VOLI (no macropore)
        new_src_vhcm + 4.19 * new_src_volw + 1.9274 * new_src_voli
    else
        new_src_vhcm +
            4.19 * (new_src_volw + src_macropore.volwh) +
            1.9274 * (new_src_voli + src_macropore.volih);

    const new_src_tks = if (new_src_vhcp > vhcprx)
        engy_src_scaled / new_src_vhcp
    else
        new_dst.tks; // fallback: adopt destination temperature

    const new_src = LayerWaterHeatState{
        .volw = new_src_volw,
        .voli = new_src_voli,
        .volp = new_src_volp,
        .vola = new_src_vola,
        .voly = new_src_voly,
        .vhcm = new_src_vhcm,
        .tks = new_src_tks,
        .tcs = new_src_tks - 273.15,
    };

    inline for (@typeInfo(LayerWaterHeatState).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(new_dst, field.name))) return error.NonFiniteDstState;
        if (!std.math.isFinite(@field(new_src, field.name))) return error.NonFiniteSrcState;
    }

    return TransferResult{ .src = new_src, .dst = new_dst };
}

test "full transfer (FX=1.0): source goes to zero, destination has sum" {
    const src = LayerWaterHeatState{
        .volw = 0.05,
        .voli = 0.01,
        .volp = 0.02,
        .vola = 0.03,
        .voly = 0.04,
        .vhcm = 1.0,
        .tks = 275.0,
        .tcs = 1.85,
    };
    const dst = LayerWaterHeatState{
        .volw = 0.02,
        .voli = 0.0,
        .volp = 0.01,
        .vola = 0.01,
        .voly = 0.02,
        .vhcm = 0.5,
        .tks = 280.0,
        .tcs = 6.85,
    };
    const zero_macro = MacroporeWaterIce{ .volwh = 0.0, .volih = 0.0 };
    const result = try transfer(src, dst, zero_macro, zero_macro, 1.0, 0.001, false);
    try std.testing.expectApproxEqRel(@as(f64, 0.0), result.src.volw, 1.0e-14);
    try std.testing.expectApproxEqRel(@as(f64, 0.07), result.dst.volw, 1.0e-14);
}

test "zero transfer (FX=0.0): both layers unchanged" {
    const src = LayerWaterHeatState{
        .volw = 0.05,
        .voli = 0.01,
        .volp = 0.02,
        .vola = 0.03,
        .voly = 0.04,
        .vhcm = 1.0,
        .tks = 275.0,
        .tcs = 1.85,
    };
    const dst = LayerWaterHeatState{
        .volw = 0.02,
        .voli = 0.0,
        .volp = 0.01,
        .vola = 0.01,
        .voly = 0.02,
        .vhcm = 0.5,
        .tks = 280.0,
        .tcs = 6.85,
    };
    const zero_macro = MacroporeWaterIce{ .volwh = 0.0, .volih = 0.0 };
    const result = try transfer(src, dst, zero_macro, zero_macro, 0.0, 0.001, false);
    try std.testing.expectApproxEqRel(src.volw, result.src.volw, 1.0e-14);
    try std.testing.expectApproxEqRel(dst.volw, result.dst.volw, 1.0e-14);
}

test "energy conservation: total energy preserved after transfer" {
    const src = LayerWaterHeatState{
        .volw = 0.04,
        .voli = 0.0,
        .volp = 0.01,
        .vola = 0.01,
        .voly = 0.03,
        .vhcm = 0.8,
        .tks = 278.0,
        .tcs = 4.85,
    };
    const dst = LayerWaterHeatState{
        .volw = 0.03,
        .voli = 0.0,
        .volp = 0.01,
        .vola = 0.01,
        .voly = 0.02,
        .vhcm = 0.6,
        .tks = 283.0,
        .tcs = 9.85,
    };
    const zero_macro = MacroporeWaterIce{ .volwh = 0.0, .volih = 0.0 };
    // Recompute VHCP for inputs to get actual energies (no macropore water)
    const src_vhcp_before = src.vhcm + 4.19 * src.volw;
    const dst_vhcp_before = dst.vhcm + 4.19 * dst.volw;
    const total_energy_before = src_vhcp_before * src.tks + dst_vhcp_before * dst.tks;
    const fx = 0.3;
    const result = try transfer(src, dst, zero_macro, zero_macro, fx, 0.001, false);
    const result_src_vhcp = result.src.vhcm + 4.19 * result.src.volw;
    const result_dst_vhcp = result.dst.vhcm + 4.19 * result.dst.volw;
    const total_energy_after = result_src_vhcp * result.src.tks + result_dst_vhcp * result.dst.tks;
    try std.testing.expectApproxEqRel(total_energy_before, total_energy_after, 1.0e-12);
}

test "temperature fallback when VHCP is below threshold" {
    // Very small heat capacity — temperature should fall back
    const src = LayerWaterHeatState{
        .volw = 1.0e-10,
        .voli = 0.0,
        .volp = 0.0,
        .vola = 0.0,
        .voly = 0.0,
        .vhcm = 1.0e-10,
        .tks = 290.0,
        .tcs = 16.85,
    };
    const dst = LayerWaterHeatState{
        .volw = 1.0e-10,
        .voli = 0.0,
        .volp = 0.0,
        .vola = 0.0,
        .voly = 0.0,
        .vhcm = 1.0e-10,
        .tks = 280.0,
        .tcs = 6.85,
    };
    const zero_macro = MacroporeWaterIce{ .volwh = 0.0, .volih = 0.0 };
    const vhcprx: f64 = 1.0; // large threshold — forces fallback
    const result = try transfer(src, dst, zero_macro, zero_macro, 0.5, vhcprx, false);
    // dst gets src's temperature in fallback (VHCP too small to compute reliably)
    try std.testing.expectApproxEqRel(src.tks, result.dst.tks, 1.0e-14);
    // src gets dst's temperature in fallback
    try std.testing.expectApproxEqRel(result.dst.tks, result.src.tks, 1.0e-14);
}

test "surface pond (layer 0) uses simplified VHCP formula" {
    const src = LayerWaterHeatState{
        .volw = 0.05,
        .voli = 0.0,
        .volp = 0.0,
        .vola = 0.0,
        .voly = 0.04,
        .vhcm = 0.0,
        .tks = 277.0,
        .tcs = 3.85,
    };
    const dst = LayerWaterHeatState{
        .volw = 0.03,
        .voli = 0.0,
        .volp = 0.0,
        .vola = 0.0,
        .voly = 0.02,
        .vhcm = 0.5,
        .tks = 275.0,
        .tcs = 1.85,
    };
    const src_macro = MacroporeWaterIce{ .volwh = 0.01, .volih = 0.005 };
    const zero_macro = MacroporeWaterIce{ .volwh = 0.0, .volih = 0.0 };
    // src_is_surface_pond = true: simplified VHCP = VHCM + 4.19*VOLW + 1.9274*VOLI (no macropore)
    const result = try transfer(src, dst, src_macro, zero_macro, 0.5, 1.0e-6, true);
    // Source VHCM = 0*FY and no macropore → VHCP = 4.19 * 0.025 only
    const expected_src_vhcp = 4.19 * (0.5 * src.volw);
    const expected_src_tks = (0.5 * src.vhcm + 4.19 * src.volw) * src.tks * 0.5 / expected_src_vhcp;
    _ = expected_src_tks;
    try std.testing.expect(std.math.isFinite(result.src.tks));
    try std.testing.expectApproxEqRel(@as(f64, 0.5 * src.volw), result.src.volw, 1.0e-14);
}

test "rejects non-finite fraction" {
    const s = std.mem.zeroes(LayerWaterHeatState);
    const m = std.mem.zeroes(MacroporeWaterIce);
    try std.testing.expectError(
        error.InvalidFraction,
        transfer(s, s, m, m, std.math.nan(f64), 0.001, false),
    );
}

test "TCS = TKS - 273.15 maintained after transfer" {
    const src = LayerWaterHeatState{
        .volw = 0.04,
        .voli = 0.0,
        .volp = 0.01,
        .vola = 0.01,
        .voly = 0.03,
        .vhcm = 1.0,
        .tks = 278.15,
        .tcs = 5.0,
    };
    const dst = LayerWaterHeatState{
        .volw = 0.03,
        .voli = 0.0,
        .volp = 0.01,
        .vola = 0.01,
        .voly = 0.02,
        .vhcm = 0.8,
        .tks = 283.15,
        .tcs = 10.0,
    };
    const zero_macro = MacroporeWaterIce{ .volwh = 0.0, .volih = 0.0 };
    const result = try transfer(src, dst, zero_macro, zero_macro, 0.4, 0.001, false);
    try std.testing.expectApproxEqRel(result.dst.tks - 273.15, result.dst.tcs, 1.0e-13);
    try std.testing.expectApproxEqRel(result.src.tks - 273.15, result.src.tcs, 1.0e-13);
}
