const std = @import("std");

pub const Inputs = struct {
    canopy_surface_temperature_k: f64,
    canopy_total_water_potential_mpa: f64,
    canopy_air_vapor_concentration_m3_per_m3: f64,
    latent_surface_conductance_m3_per_step: f64,
    intercepted_water_volume_m3: f64,
    evaporation_limit_per_step: f64,
    adjusted_sensible_surface_resistance_h_per_m: f64,
    latent_surface_resistance_h_per_m: f64,
    stomatal_resistance_h_per_m: f64,
    canopy_water_capacity_change_m3_per_step: f64,
    latent_heat_of_vaporization_megajoules_per_m3: f64,
    sensible_surface_conductance_megajoules_per_k_step: f64,
    canopy_to_surface_temperature_difference_k: f64,
};

pub const Result = struct {
    surface_vapor_concentration_m3_per_m3: f64,
    residual_vapor_flux_m3_per_step: f64,
    intercepted_evaporation_m3_per_step: f64,
    potential_transpiration_m3_per_step: f64,
    net_transpiration_m3_per_step: f64,
    latent_heat_flux_megajoules_per_step: f64,
    intercepted_water_convective_heat_flux_megajoules_per_step: f64,
    sensible_heat_flux_megajoules_per_step: f64,
};

/// UPTAKE.F 999--1014. Preserves the source evaporation branch and the
/// subsequent stomatal, latent, convective-water, and sensible flux order.
pub fn calculate(inputs: Inputs) !Result {
    try validate(inputs);
    const surface_vapor =
        2.173e-3 / inputs.canopy_surface_temperature_k *
        0.61 *
        @exp(5360.0 *
            (3.661e-3 - 1.0 / inputs.canopy_surface_temperature_k)) *
        @exp(18.0 * inputs.canopy_total_water_potential_mpa /
            (8.3143 * inputs.canopy_surface_temperature_k));
    var residual_vapor =
        inputs.latent_surface_conductance_m3_per_step *
        (inputs.canopy_air_vapor_concentration_m3_per_m3 - surface_vapor);
    const intercepted_evaporation = if (residual_vapor > 0) evaporation: {
        const evaporation = residual_vapor;
        residual_vapor = 0;
        break :evaporation evaporation;
    } else evaporation: {
        const evaporation = @max(
            residual_vapor,
            -@max(
                0,
                inputs.intercepted_water_volume_m3 *
                    inputs.evaporation_limit_per_step,
            ),
        );
        residual_vapor -= evaporation;
        break :evaporation evaporation;
    };
    const transpiration_denominator =
        inputs.adjusted_sensible_surface_resistance_h_per_m +
        inputs.stomatal_resistance_h_per_m;
    if (transpiration_denominator == 0)
        return error.SingularCanopyTranspirationResistance;
    const potential_transpiration =
        residual_vapor *
        (inputs.adjusted_sensible_surface_resistance_h_per_m +
            inputs.latent_surface_resistance_h_per_m) /
        transpiration_denominator;
    const net_transpiration =
        potential_transpiration +
        inputs.canopy_water_capacity_change_m3_per_step;
    const latent_heat =
        (potential_transpiration + intercepted_evaporation) *
        inputs.latent_heat_of_vaporization_megajoules_per_m3;
    const convective_water_heat =
        intercepted_evaporation * 4.19 *
        inputs.canopy_surface_temperature_k;
    const sensible_heat =
        inputs.sensible_surface_conductance_megajoules_per_k_step *
        inputs.canopy_to_surface_temperature_difference_k;
    const result = Result{
        .surface_vapor_concentration_m3_per_m3 = surface_vapor,
        .residual_vapor_flux_m3_per_step = residual_vapor,
        .intercepted_evaporation_m3_per_step = intercepted_evaporation,
        .potential_transpiration_m3_per_step = potential_transpiration,
        .net_transpiration_m3_per_step = net_transpiration,
        .latent_heat_flux_megajoules_per_step = latent_heat,
        .intercepted_water_convective_heat_flux_megajoules_per_step = convective_water_heat,
        .sensible_heat_flux_megajoules_per_step = sensible_heat,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteCanopyVaporHeatFluxResult;
    return result;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidCanopyVaporHeatFluxInput;
    if (inputs.canopy_surface_temperature_k <= 0 or
        inputs.canopy_air_vapor_concentration_m3_per_m3 < 0 or
        inputs.latent_surface_conductance_m3_per_step < 0 or
        inputs.intercepted_water_volume_m3 < 0 or
        inputs.evaporation_limit_per_step < 0 or
        inputs.adjusted_sensible_surface_resistance_h_per_m < 0 or
        inputs.latent_surface_resistance_h_per_m < 0 or
        inputs.stomatal_resistance_h_per_m < 0 or
        inputs.latent_heat_of_vaporization_megajoules_per_m3 <= 0 or
        inputs.sensible_surface_conductance_megajoules_per_k_step < 0)
        return error.InvalidCanopyVaporHeatFluxInput;
}

fn sourceInputs() Inputs {
    return .{
        .canopy_surface_temperature_k = 298.15,
        .canopy_total_water_potential_mpa = -1,
        .canopy_air_vapor_concentration_m3_per_m3 = 0.02,
        .latent_surface_conductance_m3_per_step = 2,
        .intercepted_water_volume_m3 = 0.003,
        .evaporation_limit_per_step = 0.25,
        .adjusted_sensible_surface_resistance_h_per_m = 0.2,
        .latent_surface_resistance_h_per_m = 0.3,
        .stomatal_resistance_h_per_m = 0.8,
        .canopy_water_capacity_change_m3_per_step = 0.001,
        .latent_heat_of_vaporization_megajoules_per_m3 = 2450,
        .sensible_surface_conductance_megajoules_per_k_step = 3,
        .canopy_to_surface_temperature_difference_k = 2,
    };
}

test "positive vapor flux is entirely intercepted evaporation" {
    const inputs = sourceInputs();
    const result = try calculate(inputs);
    const surface_vapor =
        2.173e-3 / inputs.canopy_surface_temperature_k *
        0.61 *
        @exp(5360.0 *
            (3.661e-3 - 1.0 / inputs.canopy_surface_temperature_k)) *
        @exp(18.0 * inputs.canopy_total_water_potential_mpa /
            (8.3143 * inputs.canopy_surface_temperature_k));
    const initial_flux =
        inputs.latent_surface_conductance_m3_per_step *
        (inputs.canopy_air_vapor_concentration_m3_per_m3 - surface_vapor);
    try std.testing.expect(initial_flux > 0);
    try std.testing.expectEqual(surface_vapor, result.surface_vapor_concentration_m3_per_m3);
    try std.testing.expectEqual(initial_flux, result.intercepted_evaporation_m3_per_step);
    try std.testing.expectEqual(@as(f64, 0), result.residual_vapor_flux_m3_per_step);
    try std.testing.expectEqual(@as(f64, 0), result.potential_transpiration_m3_per_step);
    try std.testing.expectEqual(@as(f64, 0.001), result.net_transpiration_m3_per_step);
    try std.testing.expectEqual(@as(f64, 6), result.sensible_heat_flux_megajoules_per_step);
}

test "negative vapor flux limits intercepted loss then passes residual to transpiration" {
    var inputs = sourceInputs();
    inputs.canopy_air_vapor_concentration_m3_per_m3 = 0;
    inputs.latent_surface_conductance_m3_per_step = 100;
    const result = try calculate(inputs);
    const maximum_intercepted_loss =
        inputs.intercepted_water_volume_m3 *
        inputs.evaporation_limit_per_step;
    try std.testing.expectEqual(
        -maximum_intercepted_loss,
        result.intercepted_evaporation_m3_per_step,
    );
    try std.testing.expect(result.residual_vapor_flux_m3_per_step < 0);
    try std.testing.expectEqual(
        result.residual_vapor_flux_m3_per_step *
            (inputs.adjusted_sensible_surface_resistance_h_per_m +
                inputs.latent_surface_resistance_h_per_m) /
            (inputs.adjusted_sensible_surface_resistance_h_per_m +
                inputs.stomatal_resistance_h_per_m),
        result.potential_transpiration_m3_per_step,
    );
}

test "zero transpiration resistance fails explicitly" {
    var inputs = sourceInputs();
    inputs.adjusted_sensible_surface_resistance_h_per_m = 0;
    inputs.stomatal_resistance_h_per_m = 0;
    try std.testing.expectError(
        error.SingularCanopyTranspirationResistance,
        calculate(inputs),
    );
}
