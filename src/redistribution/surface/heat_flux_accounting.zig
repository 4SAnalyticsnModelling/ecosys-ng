const std = @import("std");

pub const Inputs = struct {
    /// TKAM. Air temperature (K).
    air_temperature_k: f64,
    /// PRECA. Rain plus surface irrigation rate (m3 h-1).
    rain_irrigation_m3_per_h: f64,
    /// PRECW. Snowfall rate (m3 h-1).
    snowfall_m3_per_h: f64,
    /// PRECU. Subsurface irrigation rate (m3 h-1).
    subsurface_irrigation_m3_per_h: f64,
    /// HEATH. Net surface heat exchange from watsub.f (MJ step-1).
    surface_heat_exchange_megajoules: f64,
    /// THFLXC. Net canopy heat exchange from watsub.f (MJ step-1).
    canopy_heat_exchange_megajoules: f64,
    /// XNFH. Timestep (h step-1).
    timestep_h: f64,
};

/// Per-step surface heat flux increments. Snow-layer latent heat terms
/// are accumulated by the caller using the slice passed to `compute`.
pub const Increments = struct {
    /// Increment to add to HEATIN: precipitation + surface/canopy + snowpack.
    heatin_megajoules: f64,
    /// Magnitude to subtract from HEATOU: subsurface irrigation heat out.
    heatou_decrement_megajoules: f64,
};

// Heat capacity coefficients (MJ m-3 K-1).
const liquid_water_heat_capacity: f64 = 4.19;
const solid_snow_heat_capacity: f64 = 2.095;

/// Direct translation of redist.f lines 4431--4437.
///
/// Accumulates surface boundary heat fluxes into HEATIN and HEATOU.
/// Source-order preserved:
///   1. Precipitation sensible heat: (4.19*TKAM*PRECA + 2.095*TKAM*PRECW)*XNFH
///   2. Surface and canopy: HEATH + THFLXC
///   3. Snowpack latent heat: sum over layers of XHFLF0(L) + XHFLV0(L)
///   4. Subsurface irrigation heat out: 4.19*TKAM*PRECU*XNFH
///
/// The snow layer loop (DO 5150 L=1,JS) is realized here by summing
/// `snow_freeze_thaw_heat_megajoules` and `snow_evap_condensation_heat_megajoules` slices.
pub fn compute(
    inputs: Inputs,
    snow_freeze_thaw_heat_megajoules: []const f64,
    snow_evap_condensation_heat_megajoules: []const f64,
) !Increments {
    if (snow_freeze_thaw_heat_megajoules.len != snow_evap_condensation_heat_megajoules.len)
        return error.SnowLatentHeatSliceLengthMismatch;
    if (!std.math.isFinite(inputs.timestep_h) or inputs.timestep_h <= 0)
        return error.InvalidSurfaceHeatFluxTimestep;
    if (!std.math.isFinite(inputs.air_temperature_k) or inputs.air_temperature_k <= 0 or
        !std.math.isFinite(inputs.rain_irrigation_m3_per_h) or
        !std.math.isFinite(inputs.snowfall_m3_per_h) or
        !std.math.isFinite(inputs.subsurface_irrigation_m3_per_h) or
        !std.math.isFinite(inputs.surface_heat_exchange_megajoules) or
        !std.math.isFinite(inputs.canopy_heat_exchange_megajoules))
        return error.InvalidSurfaceHeatFluxInput;
    if (inputs.rain_irrigation_m3_per_h < 0 or
        inputs.snowfall_m3_per_h < 0 or
        inputs.subsurface_irrigation_m3_per_h < 0)
    {
        return error.NegativeSurfaceWaterInputRate;
    }
    for (snow_freeze_thaw_heat_megajoules) |v|
        if (!std.math.isFinite(v)) return error.InvalidSnowLatentHeatFlux;
    for (snow_evap_condensation_heat_megajoules) |v|
        if (!std.math.isFinite(v)) return error.InvalidSnowLatentHeatFlux;

    // Line 4431-4432: precipitation sensible heat.
    const precip_heat = (liquid_water_heat_capacity * inputs.air_temperature_k *
        inputs.rain_irrigation_m3_per_h +
        solid_snow_heat_capacity * inputs.air_temperature_k *
            inputs.snowfall_m3_per_h) * inputs.timestep_h;

    // Lines 4433--4436: reproduce the two source statement updates and then
    // each loop iteration. Do not pre-sum these signed fluxes: that changes
    // the floating-point expression tree relative to HEATIN.
    var heatin_megajoules = precip_heat;
    heatin_megajoules = heatin_megajoules + inputs.surface_heat_exchange_megajoules +
        inputs.canopy_heat_exchange_megajoules;
    for (snow_freeze_thaw_heat_megajoules, snow_evap_condensation_heat_megajoules) |freeze, evap| {
        heatin_megajoules = heatin_megajoules + freeze + evap;
    }

    // Line 4437: subsurface irrigation heat out.
    const irrigation_heat_out = liquid_water_heat_capacity *
        inputs.air_temperature_k *
        inputs.subsurface_irrigation_m3_per_h *
        inputs.timestep_h;

    const result = Increments{
        .heatin_megajoules = heatin_megajoules,
        .heatou_decrement_megajoules = irrigation_heat_out,
    };
    inline for (@typeInfo(Increments).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteSurfaceHeatFluxIncrement;
    return result;
}

const default_inputs = Inputs{
    .air_temperature_k = 278.0,
    .rain_irrigation_m3_per_h = 0.001,
    .snowfall_m3_per_h = 0.0005,
    .subsurface_irrigation_m3_per_h = 0.0002,
    .surface_heat_exchange_megajoules = 0.5,
    .canopy_heat_exchange_megajoules = 0.1,
    .timestep_h = 1.0,
};

test "REDIST surface heat flux preserves source-order precipitation heat formula" {
    const result = try compute(default_inputs, &.{}, &.{});
    const expected_precip = (liquid_water_heat_capacity * 278.0 * 0.001 +
        solid_snow_heat_capacity * 278.0 * 0.0005) * 1.0;
    const expected_heatin = expected_precip + 0.5 + 0.1;
    try std.testing.expectApproxEqRel(expected_heatin, result.heatin_megajoules, 1.0e-14);
}

test "REDIST surface heat flux accumulates snow layer latent heat" {
    const freeze = [_]f64{ 0.05, 0.03 };
    const evap = [_]f64{ 0.02, 0.01 };
    const result = try compute(default_inputs, &freeze, &evap);
    const expected_precip = (liquid_water_heat_capacity * 278.0 * 0.001 +
        solid_snow_heat_capacity * 278.0 * 0.0005) * 1.0;
    const expected_heatin = expected_precip + 0.5 + 0.1 +
        0.05 + 0.02 + 0.03 + 0.01;
    try std.testing.expectApproxEqRel(expected_heatin, result.heatin_megajoules, 1.0e-14);
}

test "REDIST surface heat flux computes subsurface irrigation heat out" {
    const result = try compute(default_inputs, &.{}, &.{});
    const expected_out = liquid_water_heat_capacity * 278.0 * 0.0002 * 1.0;
    try std.testing.expectApproxEqRel(expected_out, result.heatou_decrement_megajoules, 1.0e-14);
}

test "REDIST surface heat flux applies runtime timestep" {
    var inp = default_inputs;
    inp.timestep_h = 0.5;
    const full = try compute(default_inputs, &.{}, &.{});
    const half = try compute(inp, &.{}, &.{});
    // Precipitation and irrigation terms scale with timestep; surface/canopy do not.
    const precip_full = (liquid_water_heat_capacity * 278.0 * 0.001 +
        solid_snow_heat_capacity * 278.0 * 0.0005) * 1.0;
    const precip_half = precip_full * 0.5;
    const surface_canopy: f64 = 0.5 + 0.1;
    try std.testing.expectApproxEqRel(
        precip_half + surface_canopy,
        half.heatin_megajoules,
        1.0e-14,
    );
    _ = full;
}

test "REDIST surface heat flux rejects mismatched snow layer slices" {
    try std.testing.expectError(
        error.SnowLatentHeatSliceLengthMismatch,
        compute(default_inputs, &.{ 0.1, 0.2 }, &.{0.1}),
    );
}

test "REDIST surface heat flux rejects negative precipitation and irrigation rates" {
    inline for (.{
        "rain_irrigation_m3_per_h",
        "snowfall_m3_per_h",
        "subsurface_irrigation_m3_per_h",
    }) |field_name| {
        var inputs = default_inputs;
        @field(inputs, field_name) = -1.0e-6;
        try std.testing.expectError(
            error.NegativeSurfaceWaterInputRate,
            compute(inputs, &.{}, &.{}),
        );
    }
}
