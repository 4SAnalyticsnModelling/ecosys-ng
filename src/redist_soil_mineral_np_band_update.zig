const std = @import("std");

/// Soil mineral N/P band micropore pools (g N or P) for one layer.
pub const SoilNPBandPools = struct {
    /// ZNH3B. NH3 in band micropores (g N).
    nh3_g: f64,
    /// ZNH4B. NH4 in band micropores (g N).
    nh4_g: f64,
    /// ZNO3B. NO3 in band micropores (g N).
    no3_g: f64,
    /// ZNO2B. NO2 in band micropores (g N).
    no2_g: f64,
    /// H1POB. HPO4 in band micropores (g P).
    hpo4_g: f64,
    /// H2POB. H2PO4 in band micropores (g P).
    h2po4_g: f64,
};

/// N/P band flux terms (g N or P step-1).
pub const SoilNPBandFluxes = struct {
    // NH3 (line 6275)
    nh3_transport_g: f64, // TN3FLB
    nh3_dissolution_g: f64, // XNBDFG
    nh3_transformation_g: f64, // TRN3B
    nh3_root_uptake_loss_g: f64, // TUPN3B (subtracted)
    nh3_subsurface_g: f64, // RN3FBU
    nh3_pore_exchange_g: f64, // XN3FXB
    nh3_bubble_g: f64, // XNBBBL
    // NH4 (line 6278)
    nh4_transport_g: f64, // TN4FLB
    nh4_adsorption_g: f64, // XNH4B
    nh4_transformation_g: f64, // TRN4B
    nh4_root_uptake_loss_g: f64, // TUPNHB (subtracted)
    nh4_subsurface_g: f64, // RN4FBU
    nh4_pore_exchange_g: f64, // XN4FXB
    // NO3 (line 6281)
    no3_transport_g: f64, // TNOFLB
    no3_nitrification_g: f64, // XNO3B
    no3_transformation_g: f64, // TRNOB
    no3_root_uptake_loss_g: f64, // TUPNOB (subtracted)
    no3_subsurface_g: f64, // RNOFBU
    no3_pore_exchange_g: f64, // XNOFXB
    // NO2 (line 6284)
    no2_transport_g: f64, // TNXFLB
    no2_nitrification_g: f64, // XNO2B
    no2_transformation_g: f64, // TRN2B
    no2_pore_exchange_g: f64, // XNXFXB
    // HPO4 (line 6286)
    hpo4_transport_g: f64, // TH1BFB
    hpo4_desorption_g: f64, // XH1BS
    hpo4_transformation_g: f64, // TRH1B
    hpo4_root_uptake_loss_g: f64, // TUPH1B (subtracted)
    hpo4_subsurface_g: f64, // RH1BBU
    hpo4_pore_exchange_g: f64, // XH1BXB
    // H2PO4 (line 6289)
    h2po4_transport_g: f64, // TH2BFB
    h2po4_desorption_g: f64, // XH2BS
    h2po4_transformation_g: f64, // TRH2B
    h2po4_root_uptake_loss_g: f64, // TUPH2B (subtracted)
    h2po4_subsurface_g: f64, // RH2BBU
    h2po4_pore_exchange_g: f64, // XH2BXB
};

/// Layer accumulator increments from band N chemistry (lines 6292-6295).
pub const LayerAccumulators = struct {
    /// RCO2O + RCH4O: respiration C added to THRE.
    thre_increment_g: f64,
    /// XN2GS: bubble N2 escape added to UN2GS.
    un2gs_increment_g: f64,
    /// RN2G: denitrification N2 production added to UN2GG and HN2GG.
    rn2g_g: f64,
};

pub const Result = struct {
    band: SoilNPBandPools,
    thre_inc: f64,
    un2gs_inc: f64,
    rn2g: f64,
};

/// Grid-cell diagnostics accumulated in ascending soil-layer order.
pub const CellDiagnostics = struct {
    thre_g: f64,
    un2gs_g: f64,
    un2gg_g: f64,
    hn2gg_g: f64,
};

/// Direct translation of REDIST lines 6275--6295 (inner body of DO 125 L loop).
pub fn update(
    band: SoilNPBandPools,
    f: SoilNPBandFluxes,
    acc: LayerAccumulators,
) !Result {
    inline for (@typeInfo(SoilNPBandPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(band, field.name)))
            return error.InvalidSoilNPBandPool;
    inline for (@typeInfo(SoilNPBandFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(f, field.name)))
            return error.InvalidSoilNPBandFlux;
    if (!std.math.isFinite(acc.thre_increment_g) or
        !std.math.isFinite(acc.un2gs_increment_g) or
        !std.math.isFinite(acc.rn2g_g))
        return error.InvalidLayerAccumulator;

    const new_band = SoilNPBandPools{
        .nh3_g = band.nh3_g + f.nh3_transport_g + f.nh3_dissolution_g + f.nh3_transformation_g - f.nh3_root_uptake_loss_g + f.nh3_subsurface_g + f.nh3_pore_exchange_g + f.nh3_bubble_g,
        .nh4_g = band.nh4_g + f.nh4_transport_g + f.nh4_adsorption_g + f.nh4_transformation_g - f.nh4_root_uptake_loss_g + f.nh4_subsurface_g + f.nh4_pore_exchange_g,
        .no3_g = band.no3_g + f.no3_transport_g + f.no3_nitrification_g + f.no3_transformation_g - f.no3_root_uptake_loss_g + f.no3_subsurface_g + f.no3_pore_exchange_g,
        .no2_g = band.no2_g + f.no2_transport_g + f.no2_nitrification_g + f.no2_transformation_g + f.no2_pore_exchange_g,
        .hpo4_g = band.hpo4_g + f.hpo4_transport_g + f.hpo4_desorption_g + f.hpo4_transformation_g - f.hpo4_root_uptake_loss_g + f.hpo4_subsurface_g + f.hpo4_pore_exchange_g,
        .h2po4_g = band.h2po4_g + f.h2po4_transport_g + f.h2po4_desorption_g + f.h2po4_transformation_g - f.h2po4_root_uptake_loss_g + f.h2po4_subsurface_g + f.h2po4_pore_exchange_g,
    };
    inline for (@typeInfo(SoilNPBandPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(new_band, field.name)))
            return error.NonFiniteSoilNPBandPool;

    return Result{
        .band = new_band,
        .thre_inc = acc.thre_increment_g,
        .un2gs_inc = acc.un2gs_increment_g,
        .rn2g = acc.rn2g_g,
    };
}

/// Applies the REDIST DO 125 body over runtime soil layers and accumulates
/// cell diagnostics in the same order as the source loop.
pub fn updateLayers(
    bands_by_layer: []SoilNPBandPools,
    fluxes_by_layer: []const SoilNPBandFluxes,
    increments_by_layer: []const LayerAccumulators,
    diagnostics: *CellDiagnostics,
) !void {
    if (bands_by_layer.len == 0 or
        fluxes_by_layer.len != bands_by_layer.len or
        increments_by_layer.len != bands_by_layer.len)
        return error.SoilNPBandDimensionMismatch;
    inline for (@typeInfo(CellDiagnostics).@"struct".fields) |field|
        if (!std.math.isFinite(@field(diagnostics, field.name)))
            return error.InvalidCellDiagnostic;

    for (bands_by_layer, fluxes_by_layer, increments_by_layer) |*band, fluxes, increments| {
        const result = try update(band.*, fluxes, increments);
        band.* = result.band;
        diagnostics.thre_g = diagnostics.thre_g + result.thre_inc;
        diagnostics.un2gs_g = diagnostics.un2gs_g + result.un2gs_inc;
        diagnostics.un2gg_g = diagnostics.un2gg_g + result.rn2g;
        diagnostics.hn2gg_g = diagnostics.hn2gg_g + result.rn2g;
        inline for (@typeInfo(CellDiagnostics).@"struct".fields) |field|
            if (!std.math.isFinite(@field(diagnostics, field.name)))
                return error.NonFiniteCellDiagnostic;
    }
}

test "REDIST soil band NH4 transport and adsorption add independently" {
    var f = std.mem.zeroes(SoilNPBandFluxes);
    f.nh4_transport_g = 1.0;
    f.nh4_adsorption_g = 0.5;
    const result = try update(
        std.mem.zeroes(SoilNPBandPools),
        f,
        LayerAccumulators{ .thre_increment_g = 0.0, .un2gs_increment_g = 0.0, .rn2g_g = 0.0 },
    );
    try std.testing.expectApproxEqRel(@as(f64, 1.5), result.band.nh4_g, 1.0e-15);
}

test "REDIST soil band NO2 has no subsurface or root-uptake term" {
    var f = std.mem.zeroes(SoilNPBandFluxes);
    f.no2_nitrification_g = 2.0;
    f.no2_transformation_g = 1.0;
    const result = try update(
        std.mem.zeroes(SoilNPBandPools),
        f,
        LayerAccumulators{ .thre_increment_g = 0.0, .un2gs_increment_g = 0.0, .rn2g_g = 0.0 },
    );
    try std.testing.expectApproxEqRel(@as(f64, 3.0), result.band.no2_g, 1.0e-15);
}

test "REDIST soil band layer accumulators pass through unchanged" {
    const acc = LayerAccumulators{ .thre_increment_g = 5.0, .un2gs_increment_g = 2.0, .rn2g_g = 1.0 };
    const result = try update(std.mem.zeroes(SoilNPBandPools), std.mem.zeroes(SoilNPBandFluxes), acc);
    try std.testing.expectEqual(@as(f64, 5.0), result.thre_inc);
    try std.testing.expectEqual(@as(f64, 1.0), result.rn2g);
}

test "REDIST soil band rejects non-finite flux" {
    var bad = std.mem.zeroes(SoilNPBandFluxes);
    bad.hpo4_transformation_g = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidSoilNPBandFlux,
        update(
            std.mem.zeroes(SoilNPBandPools),
            bad,
            LayerAccumulators{ .thre_increment_g = 0.0, .un2gs_increment_g = 0.0, .rn2g_g = 0.0 },
        ),
    );
}

test "REDIST soil band runtime layers preserve updates and diagnostic order" {
    var bands = [_]SoilNPBandPools{
        std.mem.zeroes(SoilNPBandPools),
        std.mem.zeroes(SoilNPBandPools),
    };
    var fluxes = [_]SoilNPBandFluxes{
        std.mem.zeroes(SoilNPBandFluxes),
        std.mem.zeroes(SoilNPBandFluxes),
    };
    fluxes[0].nh3_transport_g = 2.0;
    fluxes[1].hpo4_desorption_g = 3.0;
    const increments = [_]LayerAccumulators{
        .{ .thre_increment_g = 1.0e16, .un2gs_increment_g = 2.0, .rn2g_g = 3.0 },
        .{ .thre_increment_g = -1.0e16, .un2gs_increment_g = 5.0, .rn2g_g = 7.0 },
    };
    var diagnostics = CellDiagnostics{
        .thre_g = 1.0,
        .un2gs_g = 11.0,
        .un2gg_g = 13.0,
        .hn2gg_g = 17.0,
    };

    try updateLayers(&bands, &fluxes, &increments, &diagnostics);

    try std.testing.expectEqual(@as(f64, 2.0), bands[0].nh3_g);
    try std.testing.expectEqual(@as(f64, 3.0), bands[1].hpo4_g);
    try std.testing.expectEqual(@as(f64, 0.0), diagnostics.thre_g);
    try std.testing.expectEqual(@as(f64, 18.0), diagnostics.un2gs_g);
    try std.testing.expectEqual(@as(f64, 23.0), diagnostics.un2gg_g);
    try std.testing.expectEqual(@as(f64, 27.0), diagnostics.hn2gg_g);
}

test "REDIST soil band runtime layers reject dimensions and diagnostic overflow" {
    var bands = [_]SoilNPBandPools{std.mem.zeroes(SoilNPBandPools)};
    const empty_fluxes: [0]SoilNPBandFluxes = .{};
    const increments = [_]LayerAccumulators{.{
        .thre_increment_g = 0.0,
        .un2gs_increment_g = 0.0,
        .rn2g_g = 0.0,
    }};
    var diagnostics = std.mem.zeroes(CellDiagnostics);
    try std.testing.expectError(
        error.SoilNPBandDimensionMismatch,
        updateLayers(&bands, &empty_fluxes, &increments, &diagnostics),
    );

    const fluxes = [_]SoilNPBandFluxes{std.mem.zeroes(SoilNPBandFluxes)};
    const overflowing = [_]LayerAccumulators{.{
        .thre_increment_g = std.math.floatMax(f64),
        .un2gs_increment_g = 0.0,
        .rn2g_g = 0.0,
    }};
    diagnostics.thre_g = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteCellDiagnostic,
        updateLayers(&bands, &fluxes, &overflowing, &diagnostics),
    );
}
