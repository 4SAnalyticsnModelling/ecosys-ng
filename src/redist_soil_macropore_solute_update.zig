const std = @import("std");

/// Soil macropore solute pools for one layer, stored as grams of the tracked
/// element (C, O, N, H, or P). 18 fields, REDIST 6414--6431.
pub const MacroporeSolutePools = struct {
    // Gases (lines 6414-6419)
    co2_g: f64, // CO2SH
    ch4_g: f64, // CH4SH
    o2_g: f64, // OXYSH
    n2_g: f64, // Z2GSH
    n2o_g: f64, // Z2OSH
    h2_g: f64, // H2GSH
    // N/P non-band (lines 6420-6425)
    nh4_g: f64, // ZNH4SH
    nh3_g: f64, // ZNH3SH
    no3_g: f64, // ZNO3SH
    no2_g: f64, // ZNO2SH
    hpo4_g: f64, // H1PO4H
    h2po4_g: f64, // H2PO4H
    // N/P band (lines 6426-6431)
    nh4_band_g: f64, // ZNH4BH
    nh3_band_g: f64, // ZNH3BH
    no3_band_g: f64, // ZNO3BH
    no2_band_g: f64, // ZNO2BH
    hpo4_band_g: f64, // H1POBH
    h2po4_band_g: f64, // H2POBH
};

/// Macropore flux terms (macropore transport - micropore exchange, g step-1).
/// Each: T*FH* (macropore net) minus X*FX* (micro-macro exchange, positive = micro receives).
pub const MacroporeSoluteFluxes = struct {
    // Gases
    co2_macropore_g: f64, // TCOFHS
    co2_exchange_g: f64, // XCOFXS (subtracted)
    ch4_macropore_g: f64, // TCHFHS
    ch4_exchange_g: f64, // XCHFXS (subtracted)
    o2_macropore_g: f64, // TOXFHS
    o2_exchange_g: f64, // XOXFXS (subtracted)
    n2_macropore_g: f64, // TNGFHS
    n2_exchange_g: f64, // XNGFXS (subtracted)
    n2o_macropore_g: f64, // TN2FHS
    n2o_exchange_g: f64, // XN2FXS (subtracted)
    h2_macropore_g: f64, // THGFHS
    h2_exchange_g: f64, // XHGFXS (subtracted)
    // N/P non-band
    nh4_macropore_g: f64, // TN4FHS
    nh4_exchange_g: f64, // XN4FXW (subtracted)
    nh3_macropore_g: f64, // TN3FHS
    nh3_exchange_g: f64, // XN3FXW (subtracted)
    no3_macropore_g: f64, // TNOFHS
    no3_exchange_g: f64, // XNOFXW (subtracted)
    no2_macropore_g: f64, // TNXFHS
    no2_exchange_g: f64, // XNXFXS (subtracted)
    hpo4_macropore_g: f64, // TP1FHS
    hpo4_exchange_g: f64, // XH1PXS (subtracted)
    h2po4_macropore_g: f64, // TPOFHS
    h2po4_exchange_g: f64, // XH2PXS (subtracted)
    // N/P band
    nh4_band_macropore_g: f64, // TN4FHB
    nh4_band_exchange_g: f64, // XN4FXB (subtracted)
    nh3_band_macropore_g: f64, // TN3FHB
    nh3_band_exchange_g: f64, // XN3FXB (subtracted)
    no3_band_macropore_g: f64, // TNOFHB
    no3_band_exchange_g: f64, // XNOFXB (subtracted)
    no2_band_macropore_g: f64, // TNXFHB
    no2_band_exchange_g: f64, // XNXFXB (subtracted)
    hpo4_band_macropore_g: f64, // TH1BHB
    hpo4_band_exchange_g: f64, // XH1BXB (subtracted)
    h2po4_band_macropore_g: f64, // TH2BHB
    h2po4_band_exchange_g: f64, // XH2BXB (subtracted)
};

/// Direct translation of REDIST lines 6414--6431 (inner body of DO 125 L loop).
pub fn update(pools: MacroporeSolutePools, f: MacroporeSoluteFluxes) !MacroporeSolutePools {
    inline for (@typeInfo(MacroporeSolutePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(pools, field.name)))
            return error.InvalidMacroporeSolutePool;
    inline for (@typeInfo(MacroporeSoluteFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(f, field.name)))
            return error.InvalidMacroporeSoluteFlux;

    const result = MacroporeSolutePools{
        .co2_g = pools.co2_g + f.co2_macropore_g - f.co2_exchange_g,
        .ch4_g = pools.ch4_g + f.ch4_macropore_g - f.ch4_exchange_g,
        .o2_g = pools.o2_g + f.o2_macropore_g - f.o2_exchange_g,
        .n2_g = pools.n2_g + f.n2_macropore_g - f.n2_exchange_g,
        .n2o_g = pools.n2o_g + f.n2o_macropore_g - f.n2o_exchange_g,
        .h2_g = pools.h2_g + f.h2_macropore_g - f.h2_exchange_g,
        .nh4_g = pools.nh4_g + f.nh4_macropore_g - f.nh4_exchange_g,
        .nh3_g = pools.nh3_g + f.nh3_macropore_g - f.nh3_exchange_g,
        .no3_g = pools.no3_g + f.no3_macropore_g - f.no3_exchange_g,
        .no2_g = pools.no2_g + f.no2_macropore_g - f.no2_exchange_g,
        .hpo4_g = pools.hpo4_g + f.hpo4_macropore_g - f.hpo4_exchange_g,
        .h2po4_g = pools.h2po4_g + f.h2po4_macropore_g - f.h2po4_exchange_g,
        .nh4_band_g = pools.nh4_band_g + f.nh4_band_macropore_g - f.nh4_band_exchange_g,
        .nh3_band_g = pools.nh3_band_g + f.nh3_band_macropore_g - f.nh3_band_exchange_g,
        .no3_band_g = pools.no3_band_g + f.no3_band_macropore_g - f.no3_band_exchange_g,
        .no2_band_g = pools.no2_band_g + f.no2_band_macropore_g - f.no2_band_exchange_g,
        .hpo4_band_g = pools.hpo4_band_g + f.hpo4_band_macropore_g - f.hpo4_band_exchange_g,
        .h2po4_band_g = pools.h2po4_band_g + f.h2po4_band_macropore_g - f.h2po4_band_exchange_g,
    };
    inline for (@typeInfo(MacroporeSolutePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteMacroporeSolutePool;
    return result;
}

/// Runtime NU..NL traversal in source layer order.
pub fn updateLayers(pools: []MacroporeSolutePools, fluxes: []const MacroporeSoluteFluxes) !void {
    if (pools.len == 0 or pools.len != fluxes.len)
        return error.MacroporeSoluteDimensionMismatch;
    for (pools, fluxes) |*layer_pools, layer_fluxes|
        layer_pools.* = try update(layer_pools.*, layer_fluxes);
}

test "REDIST soil macropore CO2 transport adds and exchange subtracts" {
    var f = std.mem.zeroes(MacroporeSoluteFluxes);
    f.co2_macropore_g = 3.0;
    f.co2_exchange_g = 1.0;
    const result = try update(std.mem.zeroes(MacroporeSolutePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.co2_g, 1.0e-15);
}

test "REDIST soil macropore non-band and band N pools are independent" {
    var f = std.mem.zeroes(MacroporeSoluteFluxes);
    f.nh4_macropore_g = 1.0;
    f.nh4_band_macropore_g = 2.0;
    const result = try update(std.mem.zeroes(MacroporeSolutePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 1.0), result.nh4_g, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.nh4_band_g, 1.0e-15);
}

test "REDIST soil macropore rejects non-finite flux" {
    var bad = std.mem.zeroes(MacroporeSoluteFluxes);
    bad.n2o_exchange_g = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidMacroporeSoluteFlux,
        update(std.mem.zeroes(MacroporeSolutePools), bad),
    );
}

test "REDIST soil macropore runtime layers preserve gas nutrient and band mappings" {
    var pools = [_]MacroporeSolutePools{ std.mem.zeroes(MacroporeSolutePools), std.mem.zeroes(MacroporeSolutePools) };
    var fluxes = [_]MacroporeSoluteFluxes{ std.mem.zeroes(MacroporeSoluteFluxes), std.mem.zeroes(MacroporeSoluteFluxes) };
    fluxes[0].ch4_macropore_g = 10;
    fluxes[0].ch4_exchange_g = 1;
    fluxes[0].n2o_macropore_g = 11;
    fluxes[0].n2o_exchange_g = 2;
    fluxes[1].h2po4_macropore_g = 12;
    fluxes[1].h2po4_exchange_g = 3;
    fluxes[1].hpo4_band_macropore_g = 13;
    fluxes[1].hpo4_band_exchange_g = 4;
    fluxes[1].no2_band_macropore_g = 14;
    fluxes[1].no2_band_exchange_g = 5;
    try updateLayers(&pools, &fluxes);
    try std.testing.expectEqual(@as(f64, 9), pools[0].ch4_g);
    try std.testing.expectEqual(@as(f64, 9), pools[0].n2o_g);
    try std.testing.expectEqual(@as(f64, 9), pools[1].h2po4_g);
    try std.testing.expectEqual(@as(f64, 9), pools[1].hpo4_band_g);
    try std.testing.expectEqual(@as(f64, 9), pools[1].no2_band_g);
}

test "REDIST soil macropore preserves signed exchange and rejects overflow" {
    var fluxes = std.mem.zeroes(MacroporeSoluteFluxes);
    fluxes.o2_exchange_g = -2;
    const signed = try update(std.mem.zeroes(MacroporeSolutePools), fluxes);
    try std.testing.expectEqual(@as(f64, 2), signed.o2_g);

    var pools = std.mem.zeroes(MacroporeSolutePools);
    pools.nh3_band_g = std.math.floatMax(f64);
    fluxes = std.mem.zeroes(MacroporeSoluteFluxes);
    fluxes.nh3_band_macropore_g = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteMacroporeSolutePool, update(pools, fluxes));
}

test "REDIST soil macropore rejects empty or mismatched runtime layers" {
    var pools: [1]MacroporeSolutePools = .{std.mem.zeroes(MacroporeSolutePools)};
    const no_fluxes: [0]MacroporeSoluteFluxes = .{};
    try std.testing.expectError(error.MacroporeSoluteDimensionMismatch, updateLayers(&pools, &no_fluxes));
}
