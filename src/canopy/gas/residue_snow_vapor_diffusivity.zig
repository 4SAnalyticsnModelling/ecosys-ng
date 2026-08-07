const std = @import("std");

pub const Temperatures = struct {
    canopy_air_k: f64,
    surface_residue_k: f64,
    snow_layers_k: []const f64,
};

pub const Diffusivities = struct {
    canopy_air_m2_h: f64,
    surface_residue_m2_h: f64,
};

pub const Result = struct {
    temperature_factor_canopy: f64,
    temperature_factor_surface_residue: f64,
    diffusivities: Diffusivities,
};

pub const CalculationError = error{
    SnowLayerCountMismatch,
    NonFiniteInput,
    NonPositiveTemperature,
    NegativeReferenceDiffusivity,
    NonFiniteResult,
};

/// Translates `hour1.f` lines 4688--4696 for one grid cell.
///
/// `snow_vapor_diffusivity_m2_h` is caller-owned runtime memory with one
/// element per active snow layer.
pub fn calculate(
    reference_vapor_diffusivity_m2_h: f64,
    temperatures: Temperatures,
    snow_vapor_diffusivity_m2_h: []f64,
) CalculationError!Result {
    if (temperatures.snow_layers_k.len != snow_vapor_diffusivity_m2_h.len) {
        return error.SnowLayerCountMismatch;
    }
    if (!std.math.isFinite(reference_vapor_diffusivity_m2_h) or
        !std.math.isFinite(temperatures.canopy_air_k) or
        !std.math.isFinite(temperatures.surface_residue_k))
    {
        return error.NonFiniteInput;
    }
    if (reference_vapor_diffusivity_m2_h < 0.0) {
        return error.NegativeReferenceDiffusivity;
    }
    if (temperatures.canopy_air_k <= 0.0 or temperatures.surface_residue_k <= 0.0) {
        return error.NonPositiveTemperature;
    }
    for (temperatures.snow_layers_k) |temperature_k| {
        if (!std.math.isFinite(temperature_k)) return error.NonFiniteInput;
        if (temperature_k <= 0.0) return error.NonPositiveTemperature;
    }

    const canopy_factor = std.math.pow(f64, temperatures.canopy_air_k / 298.15, 1.75);
    const canopy_diffusivity = reference_vapor_diffusivity_m2_h * canopy_factor;
    const residue_factor = std.math.pow(f64, temperatures.surface_residue_k / 298.15, 1.75);
    const residue_diffusivity = reference_vapor_diffusivity_m2_h * residue_factor;

    if (!std.math.isFinite(canopy_diffusivity) or !std.math.isFinite(residue_diffusivity)) {
        return error.NonFiniteResult;
    }
    for (temperatures.snow_layers_k, snow_vapor_diffusivity_m2_h) |
        temperature_k,
        *diffusivity_m2_h,
    | {
        const snow_factor = std.math.pow(f64, temperature_k / 298.15, 1.75);
        const snow_diffusivity = reference_vapor_diffusivity_m2_h * snow_factor;
        if (!std.math.isFinite(snow_diffusivity)) return error.NonFiniteResult;
        diffusivity_m2_h.* = snow_diffusivity;
    }

    return .{
        .temperature_factor_canopy = canopy_factor,
        .temperature_factor_surface_residue = residue_factor,
        .diffusivities = .{
            .canopy_air_m2_h = canopy_diffusivity,
            .surface_residue_m2_h = residue_diffusivity,
        },
    };
}

test "canopy residue and runtime snow layers use legacy temperature exponent" {
    const reference_m2_h = 8.96e-2;
    const snow_temperatures_k = [_]f64{ 273.15, 268.15, 260.0 };
    var snow_diffusivities_m2_h: [snow_temperatures_k.len]f64 = undefined;
    const result = try calculate(reference_m2_h, .{
        .canopy_air_k = 300.0,
        .surface_residue_k = 285.0,
        .snow_layers_k = &snow_temperatures_k,
    }, &snow_diffusivities_m2_h);

    try std.testing.expectApproxEqRel(
        reference_m2_h * std.math.pow(f64, 300.0 / 298.15, 1.75),
        result.diffusivities.canopy_air_m2_h,
        1.0e-14,
    );
    try std.testing.expectApproxEqRel(
        reference_m2_h * std.math.pow(f64, 285.0 / 298.15, 1.75),
        result.diffusivities.surface_residue_m2_h,
        1.0e-14,
    );
    for (snow_temperatures_k, snow_diffusivities_m2_h) |temperature_k, actual_m2_h| {
        try std.testing.expectApproxEqRel(
            reference_m2_h * std.math.pow(f64, temperature_k / 298.15, 1.75),
            actual_m2_h,
            1.0e-14,
        );
    }
}

test "invalid snow temperature fails before output mutation" {
    const snow_temperatures_k = [_]f64{ 270.0, std.math.nan(f64) };
    var snow_diffusivities_m2_h = [_]f64{ 11.0, 12.0 };

    try std.testing.expectError(error.NonFiniteInput, calculate(8.96e-2, .{
        .canopy_air_k = 290.0,
        .surface_residue_k = 280.0,
        .snow_layers_k = &snow_temperatures_k,
    }, &snow_diffusivities_m2_h));
    try std.testing.expectEqualSlices(f64, &.{ 11.0, 12.0 }, &snow_diffusivities_m2_h);
}
