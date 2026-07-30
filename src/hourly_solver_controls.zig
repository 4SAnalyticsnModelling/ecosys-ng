const std = @import("std");

pub const Inputs = struct {
    option_water_heat_solute_iterations: u32,
    option_gas_iterations_per_water_iteration: u32,
    fire_active: bool,
    top_layer_heat_capacity_mj_per_k_by_cell: []const f64,
    top_layer_horizontal_area_m2_by_cell: []const f64,
    normal_rainfall_interception_coefficient: f64,
    fire_rainfall_interception_coefficient: f64,
};

pub const Controls = struct {
    water_heat_solute_max_iterations: u32,
    gas_max_iterations: u32,
    litter_max_iterations: u32,
    snowpack_max_iterations: u32,
    litter_under_snow_max_iterations: u32,
    rainfall_interception_coefficient: f64,
    low_top_layer_heat_capacity_present: bool,
};

/// WTHR timestep counts become nonlinear solver ceilings in ecosys-ng; this
/// function never requests a sub-hour whole-model cycle.
pub fn derive(inputs: Inputs) !Controls {
    if (inputs.option_water_heat_solute_iterations == 0 or
        inputs.option_gas_iterations_per_water_iteration == 0)
        return error.InvalidHourlySolverIterationOption;
    if (inputs.top_layer_heat_capacity_mj_per_k_by_cell.len == 0 or
        inputs.top_layer_heat_capacity_mj_per_k_by_cell.len !=
            inputs.top_layer_horizontal_area_m2_by_cell.len)
        return error.HourlySolverControlDimensionMismatch;
    inline for (.{
        inputs.normal_rainfall_interception_coefficient,
        inputs.fire_rainfall_interception_coefficient,
    }) |coefficient| if (!std.math.isFinite(coefficient) or coefficient < 0)
        return error.InvalidRainfallInterceptionCoefficient;

    var low_heat_capacity = false;
    for (
        inputs.top_layer_heat_capacity_mj_per_k_by_cell,
        inputs.top_layer_horizontal_area_m2_by_cell,
    ) |heat_capacity_mj_per_k, area_m2| {
        if (!std.math.isFinite(heat_capacity_mj_per_k) or
            heat_capacity_mj_per_k < 0 or
            !std.math.isFinite(area_m2) or area_m2 <= 0)
            return error.InvalidTopLayerHeatCapacityState;
        low_heat_capacity = low_heat_capacity or
            heat_capacity_mj_per_k < 4.19e-3 * area_m2;
    }

    const water_iterations = if (inputs.fire_active or low_heat_capacity)
        @max(@as(u32, 20), inputs.option_water_heat_solute_iterations)
    else
        inputs.option_water_heat_solute_iterations;
    const gas_iterations = std.math.mul(
        u32,
        water_iterations,
        inputs.option_gas_iterations_per_water_iteration,
    ) catch return error.HourlyGasIterationCeilingOverflow;
    return .{
        .water_heat_solute_max_iterations = water_iterations,
        .gas_max_iterations = gas_iterations,
        .litter_max_iterations = 30,
        .snowpack_max_iterations = 20,
        .litter_under_snow_max_iterations = 10,
        .rainfall_interception_coefficient = if (inputs.fire_active)
            inputs.fire_rainfall_interception_coefficient
        else
            inputs.normal_rainfall_interception_coefficient,
        .low_top_layer_heat_capacity_present = low_heat_capacity,
    };
}

fn baseInputs() Inputs {
    return .{
        .option_water_heat_solute_iterations = 8,
        .option_gas_iterations_per_water_iteration = 3,
        .fire_active = false,
        .top_layer_heat_capacity_mj_per_k_by_cell = &.{ 1, 2 },
        .top_layer_horizontal_area_m2_by_cell = &.{ 10, 20 },
        .normal_rainfall_interception_coefficient = 0.1,
        .fire_rainfall_interception_coefficient = 0.8,
    };
}

test "ordinary hour uses runtime option ceilings and source process defaults" {
    const controls = try derive(baseInputs());
    try std.testing.expectEqual(
        @as(u32, 8),
        controls.water_heat_solute_max_iterations,
    );
    try std.testing.expectEqual(@as(u32, 24), controls.gas_max_iterations);
    try std.testing.expectEqual(@as(u32, 30), controls.litter_max_iterations);
    try std.testing.expectEqual(@as(u32, 20), controls.snowpack_max_iterations);
    try std.testing.expectEqual(
        @as(u32, 10),
        controls.litter_under_snow_max_iterations,
    );
    try std.testing.expectEqual(
        @as(f64, 0.1),
        controls.rainfall_interception_coefficient,
    );
}

test "fire takes source minimum ceiling and fire interception coefficient" {
    var inputs = baseInputs();
    inputs.fire_active = true;
    const controls = try derive(inputs);
    try std.testing.expectEqual(
        @as(u32, 20),
        controls.water_heat_solute_max_iterations,
    );
    try std.testing.expectEqual(@as(u32, 60), controls.gas_max_iterations);
    try std.testing.expectEqual(
        @as(f64, 0.8),
        controls.rainfall_interception_coefficient,
    );
}

test "any runtime cell below source heat-capacity threshold raises ceiling" {
    var inputs = baseInputs();
    inputs.top_layer_heat_capacity_mj_per_k_by_cell =
        &.{ 1, 4.189e-2, 3 };
    inputs.top_layer_horizontal_area_m2_by_cell = &.{ 10, 10, 30 };
    const controls = try derive(inputs);
    try std.testing.expect(controls.low_top_layer_heat_capacity_present);
    try std.testing.expectEqual(
        @as(u32, 20),
        controls.water_heat_solute_max_iterations,
    );
}

test "low top-layer threshold is strict and exact boundary does not raise ceiling" {
    var inputs = baseInputs();
    inputs.top_layer_heat_capacity_mj_per_k_by_cell =
        &.{ 4.19e-3, 4.19e-3 };
    inputs.top_layer_horizontal_area_m2_by_cell =
        &.{ 1.0, 1.0 };
    const controls_equal = try derive(inputs);
    try std.testing.expect(!controls_equal.low_top_layer_heat_capacity_present);
    try std.testing.expectEqual(
        @as(u32, 8),
        controls_equal.water_heat_solute_max_iterations,
    );
}

test "invalid late cell and gas ceiling overflow fail derivation" {
    var inputs = baseInputs();
    inputs.top_layer_heat_capacity_mj_per_k_by_cell =
        &.{ 1, std.math.nan(f64) };
    try std.testing.expectError(
        error.InvalidTopLayerHeatCapacityState,
        derive(inputs),
    );
    inputs = baseInputs();
    inputs.option_water_heat_solute_iterations = std.math.maxInt(u32);
    inputs.option_gas_iterations_per_water_iteration = 2;
    try std.testing.expectError(
        error.HourlyGasIterationCeilingOverflow,
        derive(inputs),
    );
}

test "invalid layer area fails early before gas ceiling math" {
    var inputs = baseInputs();
    inputs.top_layer_horizontal_area_m2_by_cell = &.{ 10, 0 };
    try std.testing.expectError(
        error.InvalidTopLayerHeatCapacityState,
        derive(inputs),
    );
}
