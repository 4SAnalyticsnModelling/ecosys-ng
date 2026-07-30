const std = @import("std");

pub const Inputs = struct {
    maximum_observed_inner_iterations: usize,
    maximum_inner_iterations: usize,
    canopy_air_temperature_k: f64,
    canopy_air_vapor_concentration_m3_per_m3: f64,
    canopy_surface_temperature_k: f64,
};

pub const Result = struct {
    canopy_air_temperature_k: f64,
    canopy_air_vapor_concentration_m3_per_m3: f64,
    canopy_surface_temperature_k: f64,
    canopy_surface_temperature_c: f64,
};

/// UPTAKE.F 1358--1362 publishes the converged final canopy state. Unlike the
/// legacy 1363--1412 branch, iteration exhaustion is an explicit fatal error
/// and cannot silently replace the failed solve with manufactured defaults.
pub fn select(inputs: Inputs) !Result {
    inline for (.{
        inputs.canopy_air_temperature_k,
        inputs.canopy_air_vapor_concentration_m3_per_m3,
        inputs.canopy_surface_temperature_k,
    }) |value|
        if (!std.math.isFinite(value))
            return error.InvalidFinalCanopyStateInput;
    if (inputs.maximum_inner_iterations == 0 or
        inputs.canopy_air_temperature_k <= 0 or
        inputs.canopy_air_vapor_concentration_m3_per_m3 < 0 or
        inputs.canopy_surface_temperature_k <= 0)
        return error.InvalidFinalCanopyStateInput;
    if (inputs.maximum_observed_inner_iterations >=
        inputs.maximum_inner_iterations)
        return error.CanopyWaterEnergyNonConvergence;
    const surface_temperature_c =
        inputs.canopy_surface_temperature_k - 273.15;
    if (!std.math.isFinite(surface_temperature_c))
        return error.NonFiniteFinalCanopyStateResult;
    return .{
        .canopy_air_temperature_k = inputs.canopy_air_temperature_k,
        .canopy_air_vapor_concentration_m3_per_m3 = inputs.canopy_air_vapor_concentration_m3_per_m3,
        .canopy_surface_temperature_k = inputs.canopy_surface_temperature_k,
        .canopy_surface_temperature_c = surface_temperature_c,
    };
}

test "UPTAKE converged final canopy state is published exactly" {
    const result = try select(.{
        .maximum_observed_inner_iterations = 12,
        .maximum_inner_iterations = 200,
        .canopy_air_temperature_k = 299,
        .canopy_air_vapor_concentration_m3_per_m3 = 0.012,
        .canopy_surface_temperature_k = 300,
    });
    try std.testing.expectEqual(@as(f64, 299), result.canopy_air_temperature_k);
    try std.testing.expectEqual(@as(f64, 0.012), result.canopy_air_vapor_concentration_m3_per_m3);
    try std.testing.expectEqual(@as(f64, 300), result.canopy_surface_temperature_k);
    try std.testing.expectApproxEqAbs(@as(f64, 26.85), result.canopy_surface_temperature_c, 1e-13);
}

test "runtime iteration exhaustion fails instead of applying legacy defaults" {
    try std.testing.expectError(
        error.CanopyWaterEnergyNonConvergence,
        select(.{
            .maximum_observed_inner_iterations = 200,
            .maximum_inner_iterations = 200,
            .canopy_air_temperature_k = 299,
            .canopy_air_vapor_concentration_m3_per_m3 = 0.012,
            .canopy_surface_temperature_k = 300,
        }),
    );
}

test "non-finite final state fails explicitly" {
    try std.testing.expectError(
        error.InvalidFinalCanopyStateInput,
        select(.{
            .maximum_observed_inner_iterations = 12,
            .maximum_inner_iterations = 200,
            .canopy_air_temperature_k = std.math.nan(f64),
            .canopy_air_vapor_concentration_m3_per_m3 = 0.012,
            .canopy_surface_temperature_k = 300,
        }),
    );
}
