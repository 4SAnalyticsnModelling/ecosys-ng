const std = @import("std");

/// Dissolved gas pools at the reference mineral layer NU (g).
pub const SoilSurfaceGasPools = struct {
    /// CO2S(NU). Aqueous CO2 (g C).
    co2_g: f64,
    /// CH4S(NU). Aqueous CH4 (g C).
    ch4_g: f64,
    /// OXYS(NU). Aqueous O2 (g O).
    o2_g: f64,
    /// Z2GS(NU). Aqueous N2 (g N).
    n2_g: f64,
    /// Z2OS(NU). Aqueous N2O (g N).
    n2o_g: f64,
    /// ZNH3S(NU). Aqueous NH3 (g N).
    nh3_g: f64,
    /// ZNH3B(NU). Band NH3 (g N).
    nh3_band_g: f64,
    /// H2GS(NU). Aqueous H2 (g H).
    h2_g: f64,
};

/// Signed surface gas diffusion increments from `trnsfr.f`, expressed as grams
/// of the named element per model step. Positive values enter layer NU.
pub const SurfaceDiffusionFluxes = struct {
    /// XCODFS. CO2 surface diffusion.
    co2_g: f64,
    /// XCHDFS. CH4 surface diffusion.
    ch4_g: f64,
    /// XOXDFS. O2 surface diffusion.
    o2_g: f64,
    /// XNGDFS. N2 surface diffusion.
    n2_g: f64,
    /// XN2DFS. N2O surface diffusion.
    n2o_g: f64,
    /// XN3DFS. NH3 surface diffusion.
    nh3_g: f64,
    /// XNBDFS. NH4 (converted to band gas) surface diffusion.
    nh3_band_g: f64,
    /// XHGDFS. H2 surface diffusion.
    h2_g: f64,
};

/// Hourly cumulative accumulators updated alongside surface gas.
pub const GasAccumulators = struct {
    /// THRE. Total heterotrophic respiration (g C).
    thre_g: f64,
    /// UN2GG. N2 cumulative (g N).
    un2gg_g: f64,
    /// HN2GG. Hourly N2 (g N).
    hn2gg_g: f64,
};

/// Litter-layer N2 transformation rate used for accumulator updates.
pub const LitterN2Rate = struct {
    /// RCO2O(0). CO2 respiration rate in litter (g C step-1).
    co2_respiration_g: f64,
    /// RCH4O(0). CH4 respiration rate in litter (g C step-1).
    ch4_respiration_g: f64,
    /// RN2G(0). N2 production rate in litter (g N step-1).
    n2_production_g: f64,
};

pub const Result = struct {
    pools: SoilSurfaceGasPools,
    accumulators: GasAccumulators,
};

/// Direct translation of redist.f lines 4907--4917.
///
/// Adds surface-diffusion fluxes to dissolved gas pools at soil reference layer
/// NU and updates THRE/UN2GG/HN2GG from litter-layer transformation rates.
pub fn update(
    pools: SoilSurfaceGasPools,
    fluxes: SurfaceDiffusionFluxes,
    accumulators: GasAccumulators,
    litter: LitterN2Rate,
) !Result {
    inline for (@typeInfo(SoilSurfaceGasPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(pools, field.name)))
            return error.InvalidSoilSurfaceGasPoolInput;
    inline for (@typeInfo(SurfaceDiffusionFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(fluxes, field.name)))
            return error.InvalidSoilSurfaceGasFluxInput;
    inline for (@typeInfo(GasAccumulators).@"struct".fields) |field|
        if (!std.math.isFinite(@field(accumulators, field.name)))
            return error.InvalidSoilSurfaceGasFluxInput;
    inline for (@typeInfo(LitterN2Rate).@"struct".fields) |field|
        if (!std.math.isFinite(@field(litter, field.name)))
            return error.InvalidSoilSurfaceGasFluxInput;

    const updated_pools = SoilSurfaceGasPools{
        .co2_g = pools.co2_g + fluxes.co2_g,
        .ch4_g = pools.ch4_g + fluxes.ch4_g,
        .o2_g = pools.o2_g + fluxes.o2_g,
        .n2_g = pools.n2_g + fluxes.n2_g,
        .n2o_g = pools.n2o_g + fluxes.n2o_g,
        .nh3_g = pools.nh3_g + fluxes.nh3_g,
        .nh3_band_g = pools.nh3_band_g + fluxes.nh3_band_g,
        .h2_g = pools.h2_g + fluxes.h2_g,
    };
    const updated_acc = GasAccumulators{
        .thre_g = accumulators.thre_g + litter.co2_respiration_g + litter.ch4_respiration_g,
        .un2gg_g = accumulators.un2gg_g + litter.n2_production_g,
        .hn2gg_g = accumulators.hn2gg_g + litter.n2_production_g,
    };

    const result = Result{ .pools = updated_pools, .accumulators = updated_acc };
    inline for (@typeInfo(SoilSurfaceGasPools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.pools, field.name)))
            return error.NonFiniteSoilSurfaceGasPool;
    inline for (@typeInfo(GasAccumulators).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.accumulators, field.name)))
            return error.NonFiniteSoilSurfaceGasAccumulator;
    return result;
}

test "REDIST soil surface gas pools accumulate surface diffusion fluxes" {
    const result = try update(
        std.mem.zeroes(SoilSurfaceGasPools),
        .{ .co2_g = 1.0, .ch4_g = 2.0, .o2_g = 3.0, .n2_g = 4.0, .n2o_g = 5.0, .nh3_g = 6.0, .nh3_band_g = 7.0, .h2_g = 8.0 },
        std.mem.zeroes(GasAccumulators),
        std.mem.zeroes(LitterN2Rate),
    );
    try std.testing.expectEqual(@as(f64, 1.0), result.pools.co2_g);
    try std.testing.expectEqual(@as(f64, 2.0), result.pools.ch4_g);
    try std.testing.expectEqual(@as(f64, 3.0), result.pools.o2_g);
    try std.testing.expectEqual(@as(f64, 4.0), result.pools.n2_g);
    try std.testing.expectEqual(@as(f64, 5.0), result.pools.n2o_g);
    try std.testing.expectEqual(@as(f64, 6.0), result.pools.nh3_g);
    try std.testing.expectEqual(@as(f64, 7.0), result.pools.nh3_band_g);
    try std.testing.expectEqual(@as(f64, 8.0), result.pools.h2_g);
}

test "REDIST soil surface gas THRE sums CO2 and CH4 respiration" {
    const result = try update(
        std.mem.zeroes(SoilSurfaceGasPools),
        std.mem.zeroes(SurfaceDiffusionFluxes),
        std.mem.zeroes(GasAccumulators),
        .{ .co2_respiration_g = 1.5, .ch4_respiration_g = 0.5, .n2_production_g = 0.0 },
    );
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.accumulators.thre_g, 1.0e-15);
}

test "REDIST soil surface gas UN2GG and HN2GG get the same N2 increment" {
    const result = try update(
        std.mem.zeroes(SoilSurfaceGasPools),
        std.mem.zeroes(SurfaceDiffusionFluxes),
        std.mem.zeroes(GasAccumulators),
        .{ .co2_respiration_g = 0.0, .ch4_respiration_g = 0.0, .n2_production_g = 3.0 },
    );
    try std.testing.expectEqual(result.accumulators.un2gg_g, result.accumulators.hn2gg_g);
    try std.testing.expectApproxEqRel(@as(f64, 3.0), result.accumulators.un2gg_g, 1.0e-15);
}

test "REDIST soil surface gas preserves signed outward diffusion" {
    const pools = SoilSurfaceGasPools{ .co2_g = 10.0, .ch4_g = 10.0, .o2_g = 10.0, .n2_g = 10.0, .n2o_g = 10.0, .nh3_g = 10.0, .nh3_band_g = 10.0, .h2_g = 10.0 };
    const fluxes = SurfaceDiffusionFluxes{ .co2_g = -1.0, .ch4_g = -2.0, .o2_g = -3.0, .n2_g = -4.0, .n2o_g = -5.0, .nh3_g = -6.0, .nh3_band_g = -7.0, .h2_g = -8.0 };
    const result = try update(pools, fluxes, std.mem.zeroes(GasAccumulators), std.mem.zeroes(LitterN2Rate));
    try std.testing.expectEqual(@as(f64, 9.0), result.pools.co2_g);
    try std.testing.expectEqual(@as(f64, 2.0), result.pools.h2_g);
}

test "REDIST soil surface gas retains prior accumulator values" {
    const accumulators = GasAccumulators{ .thre_g = 10.0, .un2gg_g = 20.0, .hn2gg_g = 30.0 };
    const litter = LitterN2Rate{ .co2_respiration_g = 1.0, .ch4_respiration_g = 2.0, .n2_production_g = 3.0 };
    const result = try update(std.mem.zeroes(SoilSurfaceGasPools), std.mem.zeroes(SurfaceDiffusionFluxes), accumulators, litter);
    try std.testing.expectEqual(@as(f64, 13.0), result.accumulators.thre_g);
    try std.testing.expectEqual(@as(f64, 23.0), result.accumulators.un2gg_g);
    try std.testing.expectEqual(@as(f64, 33.0), result.accumulators.hn2gg_g);
}

test "REDIST soil surface gas rejects invalid flux and transformation" {
    var fluxes = std.mem.zeroes(SurfaceDiffusionFluxes);
    fluxes.n2o_g = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidSoilSurfaceGasFluxInput,
        update(std.mem.zeroes(SoilSurfaceGasPools), fluxes, std.mem.zeroes(GasAccumulators), std.mem.zeroes(LitterN2Rate)),
    );

    var litter = std.mem.zeroes(LitterN2Rate);
    litter.n2_production_g = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidSoilSurfaceGasFluxInput,
        update(std.mem.zeroes(SoilSurfaceGasPools), std.mem.zeroes(SurfaceDiffusionFluxes), std.mem.zeroes(GasAccumulators), litter),
    );
}

test "REDIST soil surface gas rejects accumulator overflow" {
    var accumulators = std.mem.zeroes(GasAccumulators);
    accumulators.thre_g = std.math.floatMax(f64);
    var litter = std.mem.zeroes(LitterN2Rate);
    litter.co2_respiration_g = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSoilSurfaceGasAccumulator,
        update(std.mem.zeroes(SoilSurfaceGasPools), std.mem.zeroes(SurfaceDiffusionFluxes), accumulators, litter),
    );
}
