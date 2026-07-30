const std = @import("std");

pub const Inputs = struct {
    dry_canopy_heat_capacity_mj_per_k: f64,
    minimum_canopy_heat_capacity_mj_per_k: f64,
    canopy_radiation_fraction: f64,
    minimum_active_radiation_fraction: f64,
    ambient_to_canopy_temperature_difference_k: f64,
    sensible_boundary_conductance_mj_per_m_k_step: f64,
    aerodynamic_resistance_h_per_m: f64,
    atmospheric_vapor_concentration_m3_per_m3: f64,
    canopy_air_vapor_concentration_m3_per_m3: f64,
    latent_boundary_conductance_m2_per_step: f64,
    canopy_air_water_volume_m3: f64,
    canopy_sensible_heat_flux_mj_per_step: f64,
    canopy_transpiration_m3_per_step: f64,
    intercepted_evaporation_m3_per_step: f64,
    ground_sensible_conductance_mj_per_k_step: f64,
    canopy_air_temperature_k: f64,
    ground_air_temperature_k: f64,
    ground_latent_conductance_m3_per_step: f64,
    ground_vapor_concentration_m3_per_m3: f64,
    lateral_sensible_heat_flux_mj_per_h: f64,
    lateral_vapor_flux_m3_per_h: f64,
    energy_timestep_h_per_step: f64,
    canopy_air_heat_capacity_mj_per_k: f64,
    ambient_air_temperature_k: f64,
    intercepted_water_volume_m3: f64,
    canopy_water_volume_m3: f64,
    root_water_uptake_m3_per_step: f64,
};

pub const Result = struct {
    atmosphere_sensible_heat_flux_mj_per_step: f64,
    atmosphere_vapor_flux_m3_per_step: f64,
    net_canopy_sensible_heat_flux_mj_per_step: f64,
    net_canopy_vapor_flux_m3_per_step: f64,
    ground_sensible_heat_flux_mj_per_step: f64,
    ground_vapor_flux_m3_per_step: f64,
    lateral_sensible_heat_flux_mj_per_step: f64,
    lateral_vapor_flux_m3_per_step: f64,
    canopy_air_temperature_k: f64,
    saturated_vapor_concentration_m3_per_m3: f64,
    canopy_air_vapor_concentration_m3_per_m3: f64,
    intercepted_water_volume_m3: f64,
    canopy_water_volume_m3: f64,
};

/// UPTAKE.F 1263--1291. Updates one species canopy-air contribution and its
/// water stores in exact source branch and operation order.
pub fn calculate(inputs: Inputs) !Result {
    try validate(inputs);
    var atmosphere_sensible: f64 = 0;
    var atmosphere_vapor: f64 = 0;
    var net_sensible: f64 = 0;
    var net_vapor: f64 = 0;
    var ground_sensible: f64 = 0;
    var ground_vapor: f64 = 0;
    var lateral_sensible: f64 = 0;
    var lateral_vapor: f64 = 0;
    var air_temperature = inputs.ambient_air_temperature_k;
    var saturated_vapor = inputs.atmospheric_vapor_concentration_m3_per_m3;
    var air_vapor = inputs.atmospheric_vapor_concentration_m3_per_m3;

    if (inputs.dry_canopy_heat_capacity_mj_per_k >
        inputs.minimum_canopy_heat_capacity_mj_per_k and
        inputs.canopy_radiation_fraction >
            inputs.minimum_active_radiation_fraction)
    {
        if (inputs.aerodynamic_resistance_h_per_m == 0 or
            inputs.canopy_air_heat_capacity_mj_per_k == 0 or
            inputs.canopy_air_water_volume_m3 == 0)
            return error.SingularCanopyAirStateUpdate;
        atmosphere_sensible =
            inputs.ambient_to_canopy_temperature_difference_k *
            inputs.sensible_boundary_conductance_mj_per_m_k_step /
            inputs.aerodynamic_resistance_h_per_m *
            inputs.canopy_radiation_fraction;
        const vapor_difference =
            inputs.atmospheric_vapor_concentration_m3_per_m3 -
            inputs.canopy_air_vapor_concentration_m3_per_m3;
        atmosphere_vapor = vapor_difference * @min(
            inputs.latent_boundary_conductance_m2_per_step /
                inputs.aerodynamic_resistance_h_per_m,
            inputs.canopy_air_water_volume_m3 *
                inputs.canopy_radiation_fraction,
        );
        net_sensible =
            atmosphere_sensible -
            inputs.canopy_sensible_heat_flux_mj_per_step;
        net_vapor =
            atmosphere_vapor -
            inputs.canopy_transpiration_m3_per_step -
            inputs.intercepted_evaporation_m3_per_step;
        ground_sensible =
            inputs.ground_sensible_conductance_mj_per_k_step *
            (inputs.canopy_air_temperature_k -
                inputs.ground_air_temperature_k) *
            inputs.canopy_radiation_fraction;
        ground_vapor =
            inputs.ground_latent_conductance_m3_per_step *
            (inputs.canopy_air_vapor_concentration_m3_per_m3 -
                inputs.ground_vapor_concentration_m3_per_m3) *
            inputs.canopy_radiation_fraction;
        lateral_sensible =
            inputs.lateral_sensible_heat_flux_mj_per_h *
            inputs.canopy_radiation_fraction *
            inputs.energy_timestep_h_per_step;
        lateral_vapor =
            inputs.lateral_vapor_flux_m3_per_h *
            inputs.canopy_radiation_fraction *
            inputs.energy_timestep_h_per_step;
        air_temperature =
            inputs.canopy_air_temperature_k +
            (net_sensible - ground_sensible + lateral_sensible) /
                (inputs.canopy_air_heat_capacity_mj_per_k *
                    inputs.canopy_radiation_fraction);
        if (air_temperature <= 0)
            return error.InvalidCanopyAirTemperatureResult;
        saturated_vapor =
            2.173e-3 / air_temperature *
            0.61 *
            @exp(5360.0 * (3.661e-3 - 1.0 / air_temperature));
        air_vapor = @max(
            0,
            @min(
                saturated_vapor,
                inputs.canopy_air_vapor_concentration_m3_per_m3 +
                    (net_vapor - ground_vapor + lateral_vapor) /
                        (inputs.canopy_air_water_volume_m3 *
                            inputs.canopy_radiation_fraction),
            ),
        );
    }
    const intercepted_water =
        inputs.intercepted_water_volume_m3 +
        inputs.intercepted_evaporation_m3_per_step;
    const canopy_water =
        inputs.canopy_water_volume_m3 +
        inputs.canopy_transpiration_m3_per_step -
        inputs.root_water_uptake_m3_per_step;
    const result = Result{
        .atmosphere_sensible_heat_flux_mj_per_step = atmosphere_sensible,
        .atmosphere_vapor_flux_m3_per_step = atmosphere_vapor,
        .net_canopy_sensible_heat_flux_mj_per_step = net_sensible,
        .net_canopy_vapor_flux_m3_per_step = net_vapor,
        .ground_sensible_heat_flux_mj_per_step = ground_sensible,
        .ground_vapor_flux_m3_per_step = ground_vapor,
        .lateral_sensible_heat_flux_mj_per_step = lateral_sensible,
        .lateral_vapor_flux_m3_per_step = lateral_vapor,
        .canopy_air_temperature_k = air_temperature,
        .saturated_vapor_concentration_m3_per_m3 = saturated_vapor,
        .canopy_air_vapor_concentration_m3_per_m3 = air_vapor,
        .intercepted_water_volume_m3 = intercepted_water,
        .canopy_water_volume_m3 = canopy_water,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteCanopyAirStateUpdateResult;
    return result;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidCanopyAirStateUpdateInput;
    if (inputs.dry_canopy_heat_capacity_mj_per_k < 0 or
        inputs.minimum_canopy_heat_capacity_mj_per_k < 0 or
        inputs.canopy_radiation_fraction < 0 or
        inputs.minimum_active_radiation_fraction < 0 or
        inputs.aerodynamic_resistance_h_per_m < 0 or
        inputs.atmospheric_vapor_concentration_m3_per_m3 < 0 or
        inputs.canopy_air_vapor_concentration_m3_per_m3 < 0 or
        inputs.canopy_air_water_volume_m3 < 0 or
        inputs.canopy_air_temperature_k <= 0 or
        inputs.ground_air_temperature_k <= 0 or
        inputs.ground_vapor_concentration_m3_per_m3 < 0 or
        inputs.energy_timestep_h_per_step < 0 or
        inputs.canopy_air_heat_capacity_mj_per_k < 0 or
        inputs.ambient_air_temperature_k <= 0 or
        inputs.intercepted_water_volume_m3 < 0 or
        inputs.canopy_water_volume_m3 < 0)
        return error.InvalidCanopyAirStateUpdateInput;
}

fn sourceInputs() Inputs {
    return .{
        .dry_canopy_heat_capacity_mj_per_k = 2,
        .minimum_canopy_heat_capacity_mj_per_k = 0.1,
        .canopy_radiation_fraction = 0.5,
        .minimum_active_radiation_fraction = 1e-3,
        .ambient_to_canopy_temperature_difference_k = 2,
        .sensible_boundary_conductance_mj_per_m_k_step = 4,
        .aerodynamic_resistance_h_per_m = 2,
        .atmospheric_vapor_concentration_m3_per_m3 = 0.02,
        .canopy_air_vapor_concentration_m3_per_m3 = 0.01,
        .latent_boundary_conductance_m2_per_step = 3,
        .canopy_air_water_volume_m3 = 4,
        .canopy_sensible_heat_flux_mj_per_step = 1,
        .canopy_transpiration_m3_per_step = -0.1,
        .intercepted_evaporation_m3_per_step = -0.02,
        .ground_sensible_conductance_mj_per_k_step = 0.5,
        .canopy_air_temperature_k = 300,
        .ground_air_temperature_k = 298,
        .ground_latent_conductance_m3_per_step = 0.4,
        .ground_vapor_concentration_m3_per_m3 = 0.008,
        .lateral_sensible_heat_flux_mj_per_h = 0.2,
        .lateral_vapor_flux_m3_per_h = 0.01,
        .energy_timestep_h_per_step = 0.5,
        .canopy_air_heat_capacity_mj_per_k = 5,
        .ambient_air_temperature_k = 302,
        .intercepted_water_volume_m3 = 0.2,
        .canopy_water_volume_m3 = 2,
        .root_water_uptake_m3_per_step = -0.15,
    };
}

test "UPTAKE active canopy air update preserves source flux order" {
    const inputs = sourceInputs();
    const result = try calculate(inputs);
    const atmosphere_sensible = 2.0 * 4.0 / 2.0 * 0.5;
    const atmosphere_vapor = 0.01 * @min(3.0 / 2.0, 4.0 * 0.5);
    const net_sensible = atmosphere_sensible - 1.0;
    const net_vapor = atmosphere_vapor - -0.1 - -0.02;
    const ground_sensible = 0.5 * (300.0 - 298.0) * 0.5;
    const ground_vapor = 0.4 * (0.01 - 0.008) * 0.5;
    const lateral_sensible = 0.2 * 0.5 * 0.5;
    const lateral_vapor = 0.01 * 0.5 * 0.5;
    const temperature =
        300.0 + (net_sensible - ground_sensible + lateral_sensible) /
            (5.0 * 0.5);
    try std.testing.expectEqual(atmosphere_sensible, result.atmosphere_sensible_heat_flux_mj_per_step);
    try std.testing.expectEqual(atmosphere_vapor, result.atmosphere_vapor_flux_m3_per_step);
    try std.testing.expectEqual(temperature, result.canopy_air_temperature_k);
    try std.testing.expect(result.canopy_air_vapor_concentration_m3_per_m3 >= 0);
    try std.testing.expect(result.canopy_air_vapor_concentration_m3_per_m3 <= result.saturated_vapor_concentration_m3_per_m3);
    _ = net_vapor;
    _ = ground_vapor;
    _ = lateral_vapor;
}

test "inactive canopy restores atmosphere while still updating water stores" {
    var inputs = sourceInputs();
    inputs.canopy_radiation_fraction = 1e-3;
    const result = try calculate(inputs);
    try std.testing.expectEqual(inputs.ambient_air_temperature_k, result.canopy_air_temperature_k);
    try std.testing.expectEqual(inputs.atmospheric_vapor_concentration_m3_per_m3, result.canopy_air_vapor_concentration_m3_per_m3);
    try std.testing.expectApproxEqAbs(@as(f64, 0.18), result.intercepted_water_volume_m3, 5e-17);
    try std.testing.expectApproxEqAbs(@as(f64, 2.05), result.canopy_water_volume_m3, 5e-16);
}

test "active zero aerodynamic resistance fails explicitly" {
    var inputs = sourceInputs();
    inputs.aerodynamic_resistance_h_per_m = 0;
    try std.testing.expectError(
        error.SingularCanopyAirStateUpdate,
        calculate(inputs),
    );
}
