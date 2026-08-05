const std = @import("std");

/// Source-default hourly fertilizer dissolution rate constants from
/// solute.f lines 120--122.
pub const HourlyDissolutionRates = struct {
    /// Legacy SPNH4H. Specific hourly NH4 fertilizer dissolution rate (h-1).
    ammonium_per_h: f64,
    /// Legacy SPNH3H. Specific hourly NH3 fertilizer dissolution rate (h-1).
    ammonia_per_h: f64,
    /// Legacy SPNO3H. Specific hourly NO3 fertilizer dissolution rate (h-1).
    nitrate_per_h: f64,
    /// Legacy SPNHUH. Specific hourly urea hydrolysis rate (mol N g C-1 h-1).
    urea_hydrolysis_per_h: f64,
};

/// Source-default urea hydrolysis kinetic parameters from solute.f line 130.
pub const UreaHydrolysisParameters = struct {
    /// Legacy DUKM. Minimum Michaelis constant for urea hydrolysis
    /// (mol N Mg-1).
    minimum_half_saturation_mol_n_per_megagram: f64,
    /// Legacy DUKI. Michaelis-Menten inhibition constant for microbial
    /// activity effects on urea hydrolysis (g C m-3 h-1).
    microbial_activity_inhibition_g_c_per_m3_h: f64,
};

/// Legacy RNHUI. Inhibition decline rate constants (h-1) indexed by urea
/// formulation type from solute.f line 132.
///   index 0 = urine (fast release)
///   index 1 = normal release
///   index 2 = slow release
pub const inhibition_decline_rates_per_h: [3]f64 = .{
    5.0e-2,
    1.0e-2,
    0.5e-2,
};

/// Runtime urea formulation corresponding to legacy IUTYP. The explicit enum
/// prevents an invalid integer from indexing the three-entry RNHUI table.
pub const UreaFormulation = enum(u2) {
    urine_fast_release = 0,
    normal_release = 1,
    slow_release = 2,
};

/// Legacy ZEROC. Numerical floor used throughout SOLUTE to guard
/// against zero-division and log-of-zero conditions (solute.f line 131).
pub const numerical_floor: f64 = 1.0e-32;

/// Per-step fertilizer dissolution rate constants scaled from hourly values
/// by the solver timestep. These are outputs of `scale`.
pub const PerStepDissolutionRates = struct {
    /// Legacy SPNH4. NH4 fertilizer dissolution rate per timestep (step-1).
    ammonium_per_step: f64,
    /// Legacy SPNH3. NH3 fertilizer dissolution rate per timestep (step-1).
    ammonia_per_step: f64,
    /// Legacy SPNO3. NO3 fertilizer dissolution rate per timestep (step-1).
    nitrate_per_step: f64,
    /// Legacy SPNHU. Urea hydrolysis rate per timestep
    /// (mol N g C-1 step-1).
    urea_hydrolysis_per_step: f64,
};

/// Direct translation of SOLUTE lines 143--146.
///
/// Scales hourly fertilizer dissolution rate constants by the solver
/// timestep `timestep_h` (legacy XNFH, h step-1). Source-order preserved:
/// SPNH4, SPNH3, SPNO3, SPNHU.
pub fn scale(
    hourly: HourlyDissolutionRates,
    timestep_h: f64,
) !PerStepDissolutionRates {
    if (!std.math.isFinite(timestep_h) or timestep_h <= 0)
        return error.InvalidFertilizerDissolutionTimestep;
    inline for (@typeInfo(HourlyDissolutionRates).@"struct".fields) |field| {
        const value = @field(hourly, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidFertilizerDissolutionRate;
    }

    const result = PerStepDissolutionRates{
        .ammonium_per_step = hourly.ammonium_per_h * timestep_h,
        .ammonia_per_step = hourly.ammonia_per_h * timestep_h,
        .nitrate_per_step = hourly.nitrate_per_h * timestep_h,
        .urea_hydrolysis_per_step = hourly.urea_hydrolysis_per_h * timestep_h,
    };
    inline for (@typeInfo(PerStepDissolutionRates).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteFertilizerDissolutionRate;
    return result;
}

/// Direct source-order translation of SOLUTE.F line 157:
/// `RNHUX = RNHUI(IUTYP(NY,NX)) * XNFH`.
///
/// The caller supplies the hourly rates so production runs remain fully
/// runtime-parameterized. Returned units are step-1.
pub fn inhibitionDeclinePerStep(
    hourly_rates_per_h: [3]f64,
    formulation: UreaFormulation,
    timestep_h: f64,
) !f64 {
    if (!std.math.isFinite(timestep_h) or timestep_h <= 0)
        return error.InvalidFertilizerDissolutionTimestep;
    for (hourly_rates_per_h) |rate_per_h| {
        if (!std.math.isFinite(rate_per_h) or rate_per_h < 0)
            return error.InvalidUreaInhibitionDeclineRate;
    }

    const formulation_index: usize = @intFromEnum(formulation);
    const decline_per_step =
        hourly_rates_per_h[formulation_index] * timestep_h;
    if (!std.math.isFinite(decline_per_step))
        return error.NonFiniteUreaInhibitionDeclineRate;
    return decline_per_step;
}

const source_hourly: HourlyDissolutionRates = .{
    .ammonium_per_h = 1.0,
    .ammonia_per_h = 1.0,
    .nitrate_per_h = 1.0,
    .urea_hydrolysis_per_h = 3.0e-2,
};

test "SOLUTE fertilizer dissolution preserves source-order timestep scaling" {
    const result = try scale(source_hourly, 1.0);
    try std.testing.expectEqual(
        source_hourly.ammonium_per_h,
        result.ammonium_per_step,
    );
    try std.testing.expectEqual(
        source_hourly.ammonia_per_h,
        result.ammonia_per_step,
    );
    try std.testing.expectEqual(
        source_hourly.nitrate_per_h,
        result.nitrate_per_step,
    );
    try std.testing.expectEqual(
        source_hourly.urea_hydrolysis_per_h,
        result.urea_hydrolysis_per_step,
    );
}

test "SOLUTE fertilizer dissolution applies runtime timestep" {
    const result = try scale(source_hourly, 0.5);
    try std.testing.expectEqual(
        @as(f64, 0.5),
        result.ammonium_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.5),
        result.ammonia_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 0.5),
        result.nitrate_per_step,
    );
    try std.testing.expectApproxEqRel(
        3.0e-2 * 0.5,
        result.urea_hydrolysis_per_step,
        1.0e-15,
    );
}

test "SOLUTE fertilizer dissolution rejects invalid timestep" {
    try std.testing.expectError(
        error.InvalidFertilizerDissolutionTimestep,
        scale(source_hourly, 0.0),
    );
    try std.testing.expectError(
        error.InvalidFertilizerDissolutionTimestep,
        scale(source_hourly, -1.0),
    );
    try std.testing.expectError(
        error.InvalidFertilizerDissolutionTimestep,
        scale(source_hourly, std.math.nan(f64)),
    );
}

test "SOLUTE fertilizer dissolution rejects negative hourly rates" {
    var hourly = source_hourly;
    hourly.ammonium_per_h = -0.1;
    try std.testing.expectError(
        error.InvalidFertilizerDissolutionRate,
        scale(hourly, 1.0),
    );
}

test "SOLUTE line 157 selects formulation before timestep scaling" {
    const timestep_h = 0.25;
    try std.testing.expectEqual(
        @as(f64, 5.0e-2 * timestep_h),
        try inhibitionDeclinePerStep(
            inhibition_decline_rates_per_h,
            .urine_fast_release,
            timestep_h,
        ),
    );
    try std.testing.expectEqual(
        @as(f64, 1.0e-2 * timestep_h),
        try inhibitionDeclinePerStep(
            inhibition_decline_rates_per_h,
            .normal_release,
            timestep_h,
        ),
    );
    try std.testing.expectEqual(
        @as(f64, 0.5e-2 * timestep_h),
        try inhibitionDeclinePerStep(
            inhibition_decline_rates_per_h,
            .slow_release,
            timestep_h,
        ),
    );
}

test "SOLUTE line 157 rejects invalid runtime rates" {
    var hourly_rates = inhibition_decline_rates_per_h;
    hourly_rates[1] = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidUreaInhibitionDeclineRate,
        inhibitionDeclinePerStep(hourly_rates, .normal_release, 1),
    );
}
