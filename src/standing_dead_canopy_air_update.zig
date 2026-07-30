const std = @import("std");

pub const Inputs = struct {
    dry_heat_capacity_mj_per_k: f64,
    minimum_air_solver_heat_capacity_mj_per_k: f64,
    combined_area_radiation_fraction: f64,
    minimum_combined_area_radiation_fraction: f64,
    atmosphere_minus_canopy_air_temperature_k: f64,
    sensible_boundary_conductance_mj_per_m_k_step: f64,
    total_aerodynamic_resistance_h_per_m: f64,
    ambient_vapor_volume_fraction: f64,
    canopy_air_vapor_volume_fraction: f64,
    latent_boundary_conductance_m2_per_step: f64,
    canopy_air_vapor_capacity_m3: f64,
    surface_sensible_heat_mj_per_step: f64,
    surface_transpiration_m3_per_step: f64,
    surface_evaporation_m3_per_step: f64,
    ground_sensible_conductance_mj_per_k_step: f64,
    canopy_air_temperature_k: f64,
    ground_air_temperature_k: f64,
    ground_vapor_conductance_m3_per_step: f64,
    ground_air_vapor_volume_fraction: f64,
    lateral_sensible_heat_rate_mj_per_h: f64,
    lateral_vapor_rate_m3_per_h: f64,
    timestep_h: f64,
    canopy_air_heat_capacity_mj_per_k: f64,
    intercepted_water_m3: f64,
    ambient_air_temperature_k: f64,
};

pub const Result = struct {
    atmosphere_sensible_heat_mj_per_step: f64,
    atmosphere_vapor_flux_m3_per_step: f64,
    net_sensible_heat_mj_per_step: f64,
    net_vapor_flux_m3_per_step: f64,
    ground_sensible_heat_mj_per_step: f64,
    ground_vapor_flux_m3_per_step: f64,
    lateral_sensible_heat_mj_per_step: f64,
    lateral_vapor_flux_m3_per_step: f64,
    canopy_air_temperature_k: f64,
    saturated_vapor_volume_fraction: f64,
    canopy_air_vapor_volume_fraction: f64,
    intercepted_water_m3: f64,
};

/// UPTAKE.F 4198--4225 standing-dead canopy-air temperature, vapor, and
/// intercepted-water update.
pub fn compute(inputs: Inputs) !Result {
    try validateAlwaysUsed(inputs);
    var result = std.mem.zeroes(Result);
    if (inputs.dry_heat_capacity_mj_per_k >
        inputs.minimum_air_solver_heat_capacity_mj_per_k and
        inputs.combined_area_radiation_fraction >
            inputs.minimum_combined_area_radiation_fraction)
    {
        try validateAdmitted(inputs);
        result.atmosphere_sensible_heat_mj_per_step =
            inputs.atmosphere_minus_canopy_air_temperature_k *
            inputs.sensible_boundary_conductance_mj_per_m_k_step /
            inputs.total_aerodynamic_resistance_h_per_m *
            inputs.combined_area_radiation_fraction;
        const atmosphere_minus_canopy_vapor =
            inputs.ambient_vapor_volume_fraction -
            inputs.canopy_air_vapor_volume_fraction;
        result.atmosphere_vapor_flux_m3_per_step =
            atmosphere_minus_canopy_vapor *
            @min(
                inputs.latent_boundary_conductance_m2_per_step /
                    inputs.total_aerodynamic_resistance_h_per_m,
                inputs.canopy_air_vapor_capacity_m3 *
                    inputs.combined_area_radiation_fraction,
            );
        result.net_sensible_heat_mj_per_step =
            result.atmosphere_sensible_heat_mj_per_step -
            inputs.surface_sensible_heat_mj_per_step;
        result.net_vapor_flux_m3_per_step =
            result.atmosphere_vapor_flux_m3_per_step -
            inputs.surface_transpiration_m3_per_step -
            inputs.surface_evaporation_m3_per_step;
        result.ground_sensible_heat_mj_per_step =
            inputs.ground_sensible_conductance_mj_per_k_step *
            (inputs.canopy_air_temperature_k - inputs.ground_air_temperature_k) *
            inputs.combined_area_radiation_fraction;
        result.ground_vapor_flux_m3_per_step =
            inputs.ground_vapor_conductance_m3_per_step *
            (inputs.canopy_air_vapor_volume_fraction -
                inputs.ground_air_vapor_volume_fraction) *
            inputs.combined_area_radiation_fraction;
        result.lateral_sensible_heat_mj_per_step =
            inputs.lateral_sensible_heat_rate_mj_per_h *
            inputs.combined_area_radiation_fraction * inputs.timestep_h;
        result.lateral_vapor_flux_m3_per_step =
            inputs.lateral_vapor_rate_m3_per_h *
            inputs.combined_area_radiation_fraction * inputs.timestep_h;
        result.canopy_air_temperature_k =
            inputs.canopy_air_temperature_k +
            (result.net_sensible_heat_mj_per_step -
                result.ground_sensible_heat_mj_per_step +
                result.lateral_sensible_heat_mj_per_step) /
                (inputs.canopy_air_heat_capacity_mj_per_k *
                    inputs.combined_area_radiation_fraction);
        result.saturated_vapor_volume_fraction =
            saturationVaporVolumeFraction(result.canopy_air_temperature_k);
        result.canopy_air_vapor_volume_fraction = @max(
            0,
            @min(
                result.saturated_vapor_volume_fraction,
                inputs.canopy_air_vapor_volume_fraction +
                    (result.net_vapor_flux_m3_per_step -
                        result.ground_vapor_flux_m3_per_step +
                        result.lateral_vapor_flux_m3_per_step) /
                        (inputs.canopy_air_vapor_capacity_m3 *
                            inputs.combined_area_radiation_fraction),
            ),
        );
    } else {
        result.canopy_air_temperature_k = inputs.ambient_air_temperature_k;
        result.saturated_vapor_volume_fraction =
            saturationVaporVolumeFraction(result.canopy_air_temperature_k);
        result.canopy_air_vapor_volume_fraction =
            inputs.ambient_vapor_volume_fraction;
    }
    result.intercepted_water_m3 =
        inputs.intercepted_water_m3 + inputs.surface_evaporation_m3_per_step;
    try validateResult(result);
    return result;
}

pub fn computeRuntimeSpecies(
    inputs: []const Inputs,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.StandingDeadCanopyAirDimensionMismatch;
    for (inputs, scratch) |species_inputs, *candidate|
        candidate.* = try compute(species_inputs);
    @memcpy(destination, scratch);
}

fn saturationVaporVolumeFraction(temperature_k: f64) f64 {
    return 2.173e-3 / temperature_k * 0.61 *
        @exp(5360 * (3.661e-3 - 1 / temperature_k));
}

fn validateAlwaysUsed(inputs: Inputs) !void {
    inline for (.{
        inputs.dry_heat_capacity_mj_per_k,
        inputs.minimum_air_solver_heat_capacity_mj_per_k,
        inputs.combined_area_radiation_fraction,
        inputs.minimum_combined_area_radiation_fraction,
        inputs.surface_evaporation_m3_per_step,
        inputs.intercepted_water_m3,
        inputs.ambient_air_temperature_k,
        inputs.ambient_vapor_volume_fraction,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteStandingDeadCanopyAirInput;
    if (inputs.dry_heat_capacity_mj_per_k < 0 or
        inputs.minimum_air_solver_heat_capacity_mj_per_k < 0 or
        inputs.combined_area_radiation_fraction < 0 or
        inputs.combined_area_radiation_fraction > 1 or
        inputs.minimum_combined_area_radiation_fraction < 0 or
        inputs.intercepted_water_m3 < 0 or
        inputs.ambient_air_temperature_k <= 0 or
        inputs.ambient_vapor_volume_fraction < 0)
        return error.InvalidStandingDeadCanopyAirInput;
}

fn validateAdmitted(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteStandingDeadCanopyAirInput;
    if (inputs.total_aerodynamic_resistance_h_per_m <= 0 or
        inputs.sensible_boundary_conductance_mj_per_m_k_step < 0 or
        inputs.latent_boundary_conductance_m2_per_step < 0 or
        inputs.canopy_air_vapor_capacity_m3 <= 0 or
        inputs.ground_sensible_conductance_mj_per_k_step < 0 or
        inputs.canopy_air_temperature_k <= 0 or
        inputs.ground_air_temperature_k <= 0 or
        inputs.ground_vapor_conductance_m3_per_step < 0 or
        inputs.ground_air_vapor_volume_fraction < 0 or
        inputs.timestep_h < 0 or inputs.canopy_air_heat_capacity_mj_per_k <= 0)
        return error.InvalidStandingDeadCanopyAirInput;
}

fn validateResult(result: Result) !void {
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteStandingDeadCanopyAirResult;
    if (result.canopy_air_temperature_k <= 0 or
        result.saturated_vapor_volume_fraction < 0 or
        result.canopy_air_vapor_volume_fraction < 0 or
        result.intercepted_water_m3 < 0)
        return error.InvalidStandingDeadCanopyAirResult;
}

fn testInputs() Inputs {
    return .{
        .dry_heat_capacity_mj_per_k = 2,
        .minimum_air_solver_heat_capacity_mj_per_k = 0.1,
        .combined_area_radiation_fraction = 0.5,
        .minimum_combined_area_radiation_fraction = 1e-3,
        .atmosphere_minus_canopy_air_temperature_k = 2,
        .sensible_boundary_conductance_mj_per_m_k_step = 4,
        .total_aerodynamic_resistance_h_per_m = 2,
        .ambient_vapor_volume_fraction = 0.02,
        .canopy_air_vapor_volume_fraction = 0.01,
        .latent_boundary_conductance_m2_per_step = 2,
        .canopy_air_vapor_capacity_m3 = 4,
        .surface_sensible_heat_mj_per_step = 0.5,
        .surface_transpiration_m3_per_step = 0,
        .surface_evaporation_m3_per_step = -0.001,
        .ground_sensible_conductance_mj_per_k_step = 0.1,
        .canopy_air_temperature_k = 280,
        .ground_air_temperature_k = 278,
        .ground_vapor_conductance_m3_per_step = 0.2,
        .ground_air_vapor_volume_fraction = 0.008,
        .lateral_sensible_heat_rate_mj_per_h = 0.4,
        .lateral_vapor_rate_m3_per_h = 0.002,
        .timestep_h = 0.25,
        .canopy_air_heat_capacity_mj_per_k = 2,
        .intercepted_water_m3 = 0.01,
        .ambient_air_temperature_k = 282,
    };
}

test "admitted standing-dead canopy air update preserves source order" {
    const result = try compute(testInputs());
    try std.testing.expectEqual(@as(f64, 2), result.atmosphere_sensible_heat_mj_per_step);
    try std.testing.expectEqual(@as(f64, 1.5), result.net_sensible_heat_mj_per_step);
    try std.testing.expectEqual(@as(f64, 0.1), result.ground_sensible_heat_mj_per_step);
    try std.testing.expectEqual(@as(f64, 0.05), result.lateral_sensible_heat_mj_per_step);
    try std.testing.expectEqual(@as(f64, 281.45), result.canopy_air_temperature_k);
    try std.testing.expect(result.canopy_air_vapor_volume_fraction <=
        result.saturated_vapor_volume_fraction);
    try std.testing.expectApproxEqAbs(@as(f64, 0.009), result.intercepted_water_m3, 1e-15);
}

test "strict threshold fallback resets fluxes and selects ambient state" {
    var inputs = testInputs();
    inputs.combined_area_radiation_fraction =
        inputs.minimum_combined_area_radiation_fraction;
    inputs.total_aerodynamic_resistance_h_per_m = std.math.nan(f64);
    const result = try compute(inputs);
    try std.testing.expectEqual(@as(f64, 282), result.canopy_air_temperature_k);
    try std.testing.expectEqual(@as(f64, 0.02), result.canopy_air_vapor_volume_fraction);
    try std.testing.expectEqual(@as(f64, 0), result.net_sensible_heat_mj_per_step);
}

test "vapor update is bounded at zero and saturation" {
    var dry = testInputs();
    dry.lateral_vapor_rate_m3_per_h = -100;
    try std.testing.expectEqual(@as(f64, 0), (try compute(dry)).canopy_air_vapor_volume_fraction);
    var wet = testInputs();
    wet.lateral_vapor_rate_m3_per_h = 100;
    const result = try compute(wet);
    try std.testing.expectEqual(result.saturated_vapor_volume_fraction, result.canopy_air_vapor_volume_fraction);
}

test "later invalid runtime species leaves destination unchanged" {
    var invalid = testInputs();
    invalid.canopy_air_heat_capacity_mj_per_k = 0;
    const inputs = [_]Inputs{ testInputs(), invalid };
    var scratch: [2]Result = undefined;
    var destination = [_]Result{ std.mem.zeroes(Result), std.mem.zeroes(Result) };
    destination[0].canopy_air_temperature_k = 241;
    destination[1].canopy_air_temperature_k = 242;
    try std.testing.expectError(
        error.InvalidStandingDeadCanopyAirInput,
        computeRuntimeSpecies(&inputs, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 241), destination[0].canopy_air_temperature_k);
    try std.testing.expectEqual(@as(f64, 242), destination[1].canopy_air_temperature_k);
}
