const std = @import("std");

/// Net erosion flux of mineral fractions (Mg step-1) at layer NU.
pub const MineralErosionFlux = struct {
    /// TSEDER. Total sediment.
    total_sediment_megagrams: f64,
    /// TSANER. Sand fraction.
    sand_megagrams: f64,
    /// TSILER. Silt fraction.
    silt_megagrams: f64,
    /// TCLAER. Clay fraction.
    clay_megagrams: f64,
    /// TCECER. Cation exchange capacity increment.
    cec_mol: f64,
    /// TAECER. Anion exchange capacity increment.
    aec_mol: f64,
};

/// Net erosion flux of fertilizer pools (mol step-1).
pub const FertilizerErosionFlux = struct {
    // Non-band
    /// TNH4ER. NH4 non-band.
    nh4_nonband_mol: f64,
    /// TNH3ER. NH3 non-band.
    nh3_nonband_mol: f64,
    /// TNHUER. Urea non-band.
    urea_nonband_mol: f64,
    /// TNO3ER. NO3 non-band.
    no3_nonband_mol: f64,
    // Band
    /// TNH4EB. NH4 band.
    nh4_band_mol: f64,
    /// TNH3EB. NH3 band.
    nh3_band_mol: f64,
    /// TNHUEB. Urea band.
    urea_band_mol: f64,
    /// TNO3EB. NO3 band.
    no3_band_mol: f64,
};

/// Mineral layer state at layer NU (Mg or mol).
pub const MineralLayerState = struct {
    /// TSED. Cumulative sediment (Mg).
    total_sediment_megagrams: f64,
    /// SAND(NU). Sand (Mg).
    sand_megagrams: f64,
    /// SILT(NU). Silt (Mg).
    silt_megagrams: f64,
    /// CLAY(NU). Clay (Mg).
    clay_megagrams: f64,
    /// XCEC(NU). CEC (mol).
    cec_mol: f64,
    /// XAEC(NU). AEC (mol).
    aec_mol: f64,
};

/// Fertilizer pool state at layer NU (mol N).
pub const FertilizerLayerState = struct {
    /// ZNH4FA(NU). NH4 non-band.
    nh4_nonband_mol: f64,
    /// ZNH3FA(NU). NH3 non-band.
    nh3_nonband_mol: f64,
    /// ZNHUFA(NU). Urea non-band.
    urea_nonband_mol: f64,
    /// ZNO3FA(NU). NO3 non-band.
    no3_nonband_mol: f64,
    /// ZNH4FB(NU). NH4 band.
    nh4_band_mol: f64,
    /// ZNH3FB(NU). NH3 band.
    nh3_band_mol: f64,
    /// ZNHUFB(NU). Urea band.
    urea_band_mol: f64,
    /// ZNO3FB(NU). NO3 band.
    no3_band_mol: f64,
};

pub const Result = struct {
    mineral: MineralLayerState,
    fertilizer: FertilizerLayerState,
};

/// Runtime disturbance option corresponding to the complete IERSNG domain.
pub const DisturbanceMode = enum {
    no_profile_effects,
    freeze_thaw,
    freeze_thaw_and_erosion,
    freeze_thaw_and_organic_change,
    freeze_thaw_erosion_and_organic_change,
};

/// Exact REDIST lines 5191--5192 erosion predicate.
pub fn erosionIsActive(
    mode: DisturbanceMode,
    total_sediment_flux_megagrams: f64,
    negligible_sediment_megagrams: f64,
) !bool {
    if (!std.math.isFinite(total_sediment_flux_megagrams) or
        !std.math.isFinite(negligible_sediment_megagrams) or negligible_sediment_megagrams < 0.0)
        return error.InvalidErosionGateInput;
    const erosion_mode = mode == .freeze_thaw_and_erosion or
        mode == .freeze_thaw_erosion_and_organic_change;
    return erosion_mode and @abs(total_sediment_flux_megagrams) > negligible_sediment_megagrams;
}

/// Direct translation of REDIST lines 5193--5219.
///
/// Caller must check IERSNG is 1 or 3 and ABS(TSEDER) > ZEROS.
pub fn apply(
    mineral: MineralLayerState,
    fertilizer: FertilizerLayerState,
    mf: MineralErosionFlux,
    ff: FertilizerErosionFlux,
) !Result {
    inline for (@typeInfo(MineralLayerState).@"struct".fields) |field|
        if (!std.math.isFinite(@field(mineral, field.name)))
            return error.InvalidErosionMineralState;
    inline for (@typeInfo(FertilizerLayerState).@"struct".fields) |field|
        if (!std.math.isFinite(@field(fertilizer, field.name)))
            return error.InvalidErosionFertilizerState;
    inline for (@typeInfo(MineralErosionFlux).@"struct".fields) |field|
        if (!std.math.isFinite(@field(mf, field.name)))
            return error.InvalidErosionMineralFlux;
    inline for (@typeInfo(FertilizerErosionFlux).@"struct".fields) |field|
        if (!std.math.isFinite(@field(ff, field.name)))
            return error.InvalidErosionFertilizerFlux;

    const new_mineral = MineralLayerState{
        .total_sediment_megagrams = mineral.total_sediment_megagrams + mf.total_sediment_megagrams,
        .sand_megagrams = mineral.sand_megagrams + mf.sand_megagrams,
        .silt_megagrams = mineral.silt_megagrams + mf.silt_megagrams,
        .clay_megagrams = mineral.clay_megagrams + mf.clay_megagrams,
        .cec_mol = mineral.cec_mol + mf.cec_mol,
        .aec_mol = mineral.aec_mol + mf.aec_mol,
    };
    const new_fertilizer = FertilizerLayerState{
        .nh4_nonband_mol = fertilizer.nh4_nonband_mol + ff.nh4_nonband_mol,
        .nh3_nonband_mol = fertilizer.nh3_nonband_mol + ff.nh3_nonband_mol,
        .urea_nonband_mol = fertilizer.urea_nonband_mol + ff.urea_nonband_mol,
        .no3_nonband_mol = fertilizer.no3_nonband_mol + ff.no3_nonband_mol,
        .nh4_band_mol = fertilizer.nh4_band_mol + ff.nh4_band_mol,
        .nh3_band_mol = fertilizer.nh3_band_mol + ff.nh3_band_mol,
        .urea_band_mol = fertilizer.urea_band_mol + ff.urea_band_mol,
        .no3_band_mol = fertilizer.no3_band_mol + ff.no3_band_mol,
    };

    const result = Result{ .mineral = new_mineral, .fertilizer = new_fertilizer };
    inline for (@typeInfo(MineralLayerState).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.mineral, field.name)))
            return error.NonFiniteErosionMineralState;
    inline for (@typeInfo(FertilizerLayerState).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.fertilizer, field.name)))
            return error.NonFiniteErosionFertilizerState;
    return result;
}

test "REDIST erosion mineral sediment and texture accumulate separately" {
    var mf = std.mem.zeroes(MineralErosionFlux);
    mf.total_sediment_megagrams = 0.5;
    mf.sand_megagrams = 0.3;
    mf.clay_megagrams = 0.2;
    const result = try apply(
        std.mem.zeroes(MineralLayerState),
        std.mem.zeroes(FertilizerLayerState),
        mf,
        std.mem.zeroes(FertilizerErosionFlux),
    );
    try std.testing.expectApproxEqRel(@as(f64, 0.5), result.mineral.total_sediment_megagrams, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 0.3), result.mineral.sand_megagrams, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 0.2), result.mineral.clay_megagrams, 1.0e-15);
}

test "REDIST erosion fertilizer non-band and band are independent" {
    var ff = std.mem.zeroes(FertilizerErosionFlux);
    ff.nh4_nonband_mol = 1.0;
    ff.nh4_band_mol = 2.0;
    const result = try apply(
        std.mem.zeroes(MineralLayerState),
        std.mem.zeroes(FertilizerLayerState),
        std.mem.zeroes(MineralErosionFlux),
        ff,
    );
    try std.testing.expectApproxEqRel(@as(f64, 1.0), result.fertilizer.nh4_nonband_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.fertilizer.nh4_band_mol, 1.0e-15);
}

test "REDIST erosion mineral rejects non-finite flux" {
    var bad = std.mem.zeroes(MineralErosionFlux);
    bad.sand_megagrams = std.math.inf(f64);
    try std.testing.expectError(
        error.InvalidErosionMineralFlux,
        apply(
            std.mem.zeroes(MineralLayerState),
            std.mem.zeroes(FertilizerLayerState),
            bad,
            std.mem.zeroes(FertilizerErosionFlux),
        ),
    );
}

test "REDIST erosion gate requires erosion mode and strict threshold excess" {
    try std.testing.expect(try erosionIsActive(.freeze_thaw_and_erosion, -2.0, 1.0));
    try std.testing.expect(try erosionIsActive(.freeze_thaw_erosion_and_organic_change, 2.0, 1.0));
    try std.testing.expect(!try erosionIsActive(.freeze_thaw, 2.0, 1.0));
    try std.testing.expect(!try erosionIsActive(.freeze_thaw_and_organic_change, 2.0, 1.0));
    try std.testing.expect(!try erosionIsActive(.freeze_thaw_and_erosion, 1.0, 1.0));
}

test "REDIST erosion gate rejects invalid threshold" {
    try std.testing.expectError(
        error.InvalidErosionGateInput,
        erosionIsActive(.freeze_thaw_and_erosion, 1.0, -1.0),
    );
}

test "REDIST erosion mineral and fertilizer map every source field" {
    var mf = std.mem.zeroes(MineralErosionFlux);
    var mineral_value: f64 = 1.0;
    inline for (@typeInfo(MineralErosionFlux).@"struct".fields) |field| {
        @field(mf, field.name) = mineral_value;
        mineral_value += 1.0;
    }
    var ff = std.mem.zeroes(FertilizerErosionFlux);
    var fertilizer_value: f64 = 11.0;
    inline for (@typeInfo(FertilizerErosionFlux).@"struct".fields) |field| {
        @field(ff, field.name) = fertilizer_value;
        fertilizer_value += 1.0;
    }
    const result = try apply(
        std.mem.zeroes(MineralLayerState),
        std.mem.zeroes(FertilizerLayerState),
        mf,
        ff,
    );
    inline for (@typeInfo(MineralLayerState).@"struct".fields) |field|
        try std.testing.expectEqual(@field(mf, field.name), @field(result.mineral, field.name));
    inline for (@typeInfo(FertilizerLayerState).@"struct".fields) |field|
        try std.testing.expectEqual(@field(ff, field.name), @field(result.fertilizer, field.name));
}

test "REDIST erosion fertilizer rejects arithmetic overflow" {
    var fertilizer = std.mem.zeroes(FertilizerLayerState);
    fertilizer.urea_band_mol = std.math.floatMax(f64);
    var ff = std.mem.zeroes(FertilizerErosionFlux);
    ff.urea_band_mol = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteErosionFertilizerState,
        apply(std.mem.zeroes(MineralLayerState), fertilizer, std.mem.zeroes(MineralErosionFlux), ff),
    );
}
