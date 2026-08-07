const std = @import("std");

pub const Inputs = struct {
    current_canopy_temperature_k: f64,
    previous_canopy_temperature_k: f64,
    temperature_convergence_tolerance_k: f64,
    canopy_total_water_potential_megapascal: f64,
    canopy_gravitational_water_potential_megapascal: f64,
    available_soil_water_m3_by_layer: []const f64,
    adjusted_soil_total_water_potential_mpa_by_layer: []const f64,
    rooted_by_root_layer: []const bool,
    root_mass_fraction_by_root_layer: []const f64,
    soil_plus_root_resistance_mpa_h_per_m_by_root_layer: []const f64,
    root_axis_count: usize,
    first_soil_layer: usize,
    end_soil_layer: usize,
    plant_population: f64,
    water_flux_timestep_h_per_step: f64,
};

pub const Result = struct {
    energy_converged: bool,
    height_adjusted_canopy_water_potential_megapascal: f64,
    total_root_water_uptake_m3_per_step: f64,
};

/// UPTAKE.F 1095--1127. On strict energy convergence, evaluates runtime-sized
/// root-axis/layer water uptake in the source root-major, layer-minor order.
pub fn calculate(
    inputs: Inputs,
    root_water_uptake_m3_per_step_by_root_layer: []f64,
) !Result {
    try validate(inputs, root_water_uptake_m3_per_step_by_root_layer.len);
    if (@abs(inputs.current_canopy_temperature_k -
        inputs.previous_canopy_temperature_k) >=
        inputs.temperature_convergence_tolerance_k)
        return .{
            .energy_converged = false,
            .height_adjusted_canopy_water_potential_megapascal = inputs.canopy_total_water_potential_megapascal -
                inputs.canopy_gravitational_water_potential_megapascal,
            .total_root_water_uptake_m3_per_step = 0,
        };

    var total_uptake: f64 = 0;
    const canopy_potential =
        inputs.canopy_total_water_potential_megapascal -
        inputs.canopy_gravitational_water_potential_megapascal;
    for (0..inputs.root_axis_count) |root_axis| {
        for (inputs.first_soil_layer..inputs.end_soil_layer) |layer| {
            const index = root_axis * inputs.available_soil_water_m3_by_layer.len +
                layer;
            if (inputs.rooted_by_root_layer[index]) {
                const maximum_uptake =
                    inputs.available_soil_water_m3_by_layer[layer] *
                    inputs.root_mass_fraction_by_root_layer[index] *
                    inputs.water_flux_timestep_h_per_step;
                const hydraulic_flux =
                    (canopy_potential -
                        inputs.adjusted_soil_total_water_potential_mpa_by_layer[layer]) /
                    inputs.soil_plus_root_resistance_mpa_h_per_m_by_root_layer[index] *
                    inputs.plant_population *
                    inputs.water_flux_timestep_h_per_step;
                root_water_uptake_m3_per_step_by_root_layer[index] = @max(
                    @min(0, -maximum_uptake),
                    @min(hydraulic_flux, maximum_uptake),
                );
                total_uptake += root_water_uptake_m3_per_step_by_root_layer[index];
            } else {
                root_water_uptake_m3_per_step_by_root_layer[index] = 0;
            }
        }
    }
    if (!std.math.isFinite(total_uptake))
        return error.NonFiniteRootWaterUptakeResult;
    return .{
        .energy_converged = true,
        .height_adjusted_canopy_water_potential_megapascal = canopy_potential,
        .total_root_water_uptake_m3_per_step = total_uptake,
    };
}

fn validate(inputs: Inputs, output_count: usize) !void {
    inline for (.{
        inputs.current_canopy_temperature_k,
        inputs.previous_canopy_temperature_k,
        inputs.temperature_convergence_tolerance_k,
        inputs.canopy_total_water_potential_megapascal,
        inputs.canopy_gravitational_water_potential_megapascal,
        inputs.plant_population,
        inputs.water_flux_timestep_h_per_step,
    }) |value|
        if (!std.math.isFinite(value))
            return error.InvalidRootWaterUptakeInput;
    const layer_count = inputs.available_soil_water_m3_by_layer.len;
    const root_layer_count = std.math.mul(
        usize,
        inputs.root_axis_count,
        layer_count,
    ) catch return error.RootWaterUptakeDimensionMismatch;
    if (inputs.adjusted_soil_total_water_potential_mpa_by_layer.len != layer_count or
        inputs.rooted_by_root_layer.len != root_layer_count or
        inputs.root_mass_fraction_by_root_layer.len != root_layer_count or
        inputs.soil_plus_root_resistance_mpa_h_per_m_by_root_layer.len != root_layer_count or
        output_count != root_layer_count or
        inputs.first_soil_layer > inputs.end_soil_layer or
        inputs.end_soil_layer > layer_count)
        return error.RootWaterUptakeDimensionMismatch;
    if (inputs.current_canopy_temperature_k <= 0 or
        inputs.previous_canopy_temperature_k <= 0 or
        inputs.temperature_convergence_tolerance_k <= 0 or
        inputs.plant_population < 0 or
        inputs.water_flux_timestep_h_per_step < 0)
        return error.InvalidRootWaterUptakeInput;
    for (inputs.available_soil_water_m3_by_layer) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootWaterUptakeInput;
    for (inputs.adjusted_soil_total_water_potential_mpa_by_layer) |value|
        if (!std.math.isFinite(value))
            return error.InvalidRootWaterUptakeInput;
    for (0..root_layer_count) |index| {
        const fraction = inputs.root_mass_fraction_by_root_layer[index];
        const resistance =
            inputs.soil_plus_root_resistance_mpa_h_per_m_by_root_layer[index];
        if (!std.math.isFinite(fraction) or fraction < 0 or
            !std.math.isFinite(resistance) or
            (inputs.rooted_by_root_layer[index] and resistance <= 0))
            return error.InvalidRootWaterUptakeInput;
    }
}

fn sourceInputs() Inputs {
    return .{
        .current_canopy_temperature_k = 300.001,
        .previous_canopy_temperature_k = 300,
        .temperature_convergence_tolerance_k = 0.01,
        .canopy_total_water_potential_megapascal = -1,
        .canopy_gravitational_water_potential_megapascal = -0.02,
        .available_soil_water_m3_by_layer = &.{ 2, 3 },
        .adjusted_soil_total_water_potential_mpa_by_layer = &.{ -0.5, -2 },
        .rooted_by_root_layer = &.{ true, false, true, true },
        .root_mass_fraction_by_root_layer = &.{ 0.2, 0.3, 0.4, 0.1 },
        .soil_plus_root_resistance_mpa_h_per_m_by_root_layer = &.{ 2, 2, 4, 5 },
        .root_axis_count = 2,
        .first_soil_layer = 0,
        .end_soil_layer = 2,
        .plant_population = 2,
        .water_flux_timestep_h_per_step = 0.5,
    };
}

test "UPTAKE root uptake preserves source clamps and accumulation order" {
    const inputs = sourceInputs();
    var uptake = [_]f64{ 9, 9, 9, 9 };
    const result = try calculate(inputs, &uptake);
    const canopy_potential = -1.0 - -0.02;
    const maximum_0 = 2.0 * 0.2 * 0.5;
    const hydraulic_0 = (canopy_potential - -0.5) / 2.0 * 2.0 * 0.5;
    const expected_0 = @max(@min(0, -maximum_0), @min(hydraulic_0, maximum_0));
    const maximum_2 = 2.0 * 0.4 * 0.5;
    const hydraulic_2 = (canopy_potential - -0.5) / 4.0 * 2.0 * 0.5;
    const expected_2 = @max(@min(0, -maximum_2), @min(hydraulic_2, maximum_2));
    const maximum_3 = 3.0 * 0.1 * 0.5;
    const hydraulic_3 = (canopy_potential - -2.0) / 5.0 * 2.0 * 0.5;
    const expected_3 = @max(@min(0, -maximum_3), @min(hydraulic_3, maximum_3));
    try std.testing.expect(result.energy_converged);
    try std.testing.expectEqual(expected_0, uptake[0]);
    try std.testing.expectEqual(@as(f64, 0), uptake[1]);
    try std.testing.expectEqual(expected_2, uptake[2]);
    try std.testing.expectApproxEqAbs(expected_3, uptake[3], 5e-17);
    try std.testing.expectApproxEqAbs(
        expected_0 + expected_2 + expected_3,
        result.total_root_water_uptake_m3_per_step,
        1e-16,
    );
}

test "strict temperature gate leaves uptake workspace untouched" {
    var inputs = sourceInputs();
    inputs.current_canopy_temperature_k = 300.02;
    var uptake = [_]f64{ 9, 9, 9, 9 };
    const result = try calculate(inputs, &uptake);
    try std.testing.expect(!result.energy_converged);
    try std.testing.expectEqualSlices(f64, &.{ 9, 9, 9, 9 }, &uptake);
}

test "rooted zero resistance fails explicitly" {
    var resistance = [_]f64{ 0, 2, 4, 5 };
    var inputs = sourceInputs();
    inputs.soil_plus_root_resistance_mpa_h_per_m_by_root_layer = &resistance;
    var uptake = [_]f64{0} ** 4;
    try std.testing.expectError(
        error.InvalidRootWaterUptakeInput,
        calculate(inputs, &uptake),
    );
}
