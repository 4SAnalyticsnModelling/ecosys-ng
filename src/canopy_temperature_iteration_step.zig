const std = @import("std");

pub const TemperatureDirection = enum {
    non_increasing,
    increasing,
};

pub const Inputs = struct {
    net_radiation_flux_megajoules_per_step: f64,
    latent_heat_flux_megajoules_per_step: f64,
    sensible_heat_flux_megajoules_per_step: f64,
    intercepted_water_convective_heat_flux_megajoules_per_step: f64,
    precipitation_convective_heat_flux_megajoules_per_step: f64,
    previous_wet_canopy_heat_capacity_megajoules_per_k: f64,
    intercepted_evaporation_m3_per_step: f64,
    foliar_water_retention_m3_per_step: f64,
    minimum_canopy_heat_capacity_megajoules_per_k: f64,
    previous_temperature_estimate_k: f64,
    canopy_surface_temperature_k: f64,
    temperature_step_fraction: f64,
    temperature_step_scale: f64,
    previous_direction: TemperatureDirection,
};

pub const Result = struct {
    storage_heat_flux_megajoules_per_step: f64,
    canopy_heat_capacity_megajoules_per_k: f64,
    previous_temperature_estimate_k: f64,
    equilibrium_temperature_k: f64,
    canopy_surface_temperature_k: f64,
    temperature_step_scale: f64,
    direction: TemperatureDirection,
};

/// UPTAKE.F 1026 and 1041--1057. One source-ordered canopy-temperature
/// iteration step; the caller owns the runtime Newton/Picard loop.
pub fn calculate(inputs: Inputs) !Result {
    try validate(inputs);
    const storage_heat =
        inputs.net_radiation_flux_megajoules_per_step +
        inputs.latent_heat_flux_megajoules_per_step +
        inputs.sensible_heat_flux_megajoules_per_step +
        inputs.intercepted_water_convective_heat_flux_megajoules_per_step +
        inputs.precipitation_convective_heat_flux_megajoules_per_step;
    const heat_capacity =
        inputs.previous_wet_canopy_heat_capacity_megajoules_per_k +
        4.19 * (inputs.intercepted_evaporation_m3_per_step +
            inputs.foliar_water_retention_m3_per_step);

    var previous_temperature = inputs.previous_temperature_estimate_k;
    var equilibrium_temperature = previous_temperature;
    var surface_temperature = previous_temperature;
    if (heat_capacity > inputs.minimum_canopy_heat_capacity_megajoules_per_k) {
        previous_temperature = inputs.canopy_surface_temperature_k;
        equilibrium_temperature =
            (previous_temperature *
                inputs.previous_wet_canopy_heat_capacity_megajoules_per_k +
                storage_heat) /
            heat_capacity;
        surface_temperature =
            previous_temperature +
            inputs.temperature_step_scale *
                inputs.temperature_step_fraction *
                (equilibrium_temperature - previous_temperature);
    }

    var step_scale = inputs.temperature_step_scale;
    if ((inputs.previous_direction == .non_increasing and
        surface_temperature > previous_temperature) or
        (inputs.previous_direction == .increasing and
            surface_temperature < previous_temperature))
        step_scale = 0.5 * step_scale;
    const direction: TemperatureDirection =
        if (surface_temperature > previous_temperature)
            .increasing
        else
            .non_increasing;
    const result = Result{
        .storage_heat_flux_megajoules_per_step = storage_heat,
        .canopy_heat_capacity_megajoules_per_k = heat_capacity,
        .previous_temperature_estimate_k = previous_temperature,
        .equilibrium_temperature_k = equilibrium_temperature,
        .canopy_surface_temperature_k = surface_temperature,
        .temperature_step_scale = step_scale,
        .direction = direction,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (field.type == TemperatureDirection) continue;
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteCanopyTemperatureIterationResult;
    }
    return result;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == TemperatureDirection) continue;
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidCanopyTemperatureIterationInput;
    }
    if (inputs.previous_wet_canopy_heat_capacity_megajoules_per_k < 0 or
        inputs.minimum_canopy_heat_capacity_megajoules_per_k < 0 or
        inputs.previous_temperature_estimate_k <= 0 or
        inputs.canopy_surface_temperature_k <= 0 or
        inputs.temperature_step_fraction < 0 or
        inputs.temperature_step_scale < 0)
        return error.InvalidCanopyTemperatureIterationInput;
}

fn sourceInputs() Inputs {
    return .{
        .net_radiation_flux_megajoules_per_step = 5,
        .latent_heat_flux_megajoules_per_step = -2,
        .sensible_heat_flux_megajoules_per_step = -1,
        .intercepted_water_convective_heat_flux_megajoules_per_step = 0.2,
        .precipitation_convective_heat_flux_megajoules_per_step = 0.3,
        .previous_wet_canopy_heat_capacity_megajoules_per_k = 10,
        .intercepted_evaporation_m3_per_step = 0.01,
        .foliar_water_retention_m3_per_step = 0.02,
        .minimum_canopy_heat_capacity_megajoules_per_k = 0.1,
        .previous_temperature_estimate_k = 297,
        .canopy_surface_temperature_k = 300,
        .temperature_step_fraction = 0.5,
        .temperature_step_scale = 1,
        .previous_direction = .non_increasing,
    };
}

test "UPTAKE canopy temperature step preserves source operation order" {
    const inputs = sourceInputs();
    const result = try calculate(inputs);
    const storage_heat = 5.0 - 2.0 - 1.0 + 0.2 + 0.3;
    const heat_capacity = 10.0 + 4.19 * (0.01 + 0.02);
    const equilibrium = (300.0 * 10.0 + storage_heat) / heat_capacity;
    const surface = 300.0 + 1.0 * 0.5 * (equilibrium - 300.0);
    try std.testing.expectEqual(storage_heat, result.storage_heat_flux_megajoules_per_step);
    try std.testing.expectEqual(heat_capacity, result.canopy_heat_capacity_megajoules_per_k);
    try std.testing.expectEqual(@as(f64, 300), result.previous_temperature_estimate_k);
    try std.testing.expectEqual(equilibrium, result.equilibrium_temperature_k);
    try std.testing.expectApproxEqAbs(
        surface,
        result.canopy_surface_temperature_k,
        1e-13,
    );
    try std.testing.expectEqual(TemperatureDirection.non_increasing, result.direction);
}

test "direction reversal halves the next source step scale" {
    var inputs = sourceInputs();
    inputs.previous_direction = .increasing;
    const result = try calculate(inputs);
    try std.testing.expect(
        result.canopy_surface_temperature_k < result.previous_temperature_estimate_k,
    );
    try std.testing.expectEqual(@as(f64, 0.5), result.temperature_step_scale);
}

test "heat capacity threshold retains incoming previous estimate" {
    var inputs = sourceInputs();
    inputs.previous_wet_canopy_heat_capacity_megajoules_per_k = 0.01;
    inputs.intercepted_evaporation_m3_per_step = 0;
    inputs.foliar_water_retention_m3_per_step = 0;
    const result = try calculate(inputs);
    try std.testing.expectEqual(@as(f64, 297), result.previous_temperature_estimate_k);
    try std.testing.expectEqual(@as(f64, 297), result.equilibrium_temperature_k);
    try std.testing.expectEqual(@as(f64, 297), result.canopy_surface_temperature_k);
}
