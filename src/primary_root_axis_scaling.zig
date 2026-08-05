const std = @import("std");

pub const State = struct {
    retained_root_carbon_g_c_per_plant: f64,
};

pub const Inputs = struct {
    total_root_carbon_g_c: f64,
    plant_population: f64,
    biological_timestep_h: f64,
    hourly_retention_fraction: f64,
    axis_scaling_exponent: f64,
};

pub const Result = struct {
    primary_root_axis_count_multiplier: f64,
};

/// Exact GROSUB lines 506--512 WTRTA/XRTN1 recurrence.
///
/// WTRTA is retained root carbon per plant (g C plant-1). The source directly
/// multiplies its hourly retention fraction by XNFH; it does not exponentiate
/// the retention fraction by the timestep.
pub fn advance(state: *State, inputs: Inputs) !Result {
    inline for (.{
        state.retained_root_carbon_g_c_per_plant,
        inputs.total_root_carbon_g_c,
        inputs.plant_population,
        inputs.biological_timestep_h,
        inputs.hourly_retention_fraction,
        inputs.axis_scaling_exponent,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidPrimaryRootAxisScalingInput;
    if (inputs.biological_timestep_h == 0 or
        inputs.hourly_retention_fraction > 1 or
        inputs.axis_scaling_exponent == 0)
        return error.InvalidPrimaryRootAxisScalingInput;

    const retained = if (inputs.plant_population > 0)
        @max(
            inputs.hourly_retention_fraction *
                state.retained_root_carbon_g_c_per_plant *
                inputs.biological_timestep_h,
            inputs.total_root_carbon_g_c / inputs.plant_population,
        )
    else
        0;
    const multiplier = @max(
        @as(f64, 1),
        std.math.pow(f64, retained, inputs.axis_scaling_exponent),
    ) * inputs.plant_population;
    if (!std.math.isFinite(retained) or !std.math.isFinite(multiplier))
        return error.PrimaryRootAxisScalingOverflow;

    state.retained_root_carbon_g_c_per_plant = retained;
    return .{ .primary_root_axis_count_multiplier = multiplier };
}

test "GROSUB retained root carbon uses exact recurrence ordering" {
    var state: State = .{ .retained_root_carbon_g_c_per_plant = 8 };
    const result = try advance(&state, .{
        .total_root_carbon_g_c = 6,
        .plant_population = 2,
        .biological_timestep_h = 0.5,
        .hourly_retention_fraction = 0.999992087,
        .axis_scaling_exponent = 0.667,
    });
    const expected_retained = 0.999992087 * 8 * 0.5;
    try std.testing.expectEqual(expected_retained, state.retained_root_carbon_g_c_per_plant);
    try std.testing.expectApproxEqAbs(
        std.math.pow(f64, expected_retained, 0.667) * 2,
        result.primary_root_axis_count_multiplier,
        1.0e-12,
    );
}

test "current root carbon per plant supplies the larger source branch" {
    var state: State = .{ .retained_root_carbon_g_c_per_plant = 1 };
    const result = try advance(&state, .{
        .total_root_carbon_g_c = 20,
        .plant_population = 4,
        .biological_timestep_h = 1,
        .hourly_retention_fraction = 0.999992087,
        .axis_scaling_exponent = 0.667,
    });
    try std.testing.expectEqual(@as(f64, 5), state.retained_root_carbon_g_c_per_plant);
    try std.testing.expectApproxEqAbs(
        std.math.pow(f64, @as(f64, 5), 0.667) * 4,
        result.primary_root_axis_count_multiplier,
        1.0e-12,
    );
}

test "zero population clears retained mass and publishes zero axes" {
    var state: State = .{ .retained_root_carbon_g_c_per_plant = 7 };
    const result = try advance(&state, .{
        .total_root_carbon_g_c = 9,
        .plant_population = 0,
        .biological_timestep_h = 1,
        .hourly_retention_fraction = 0.999992087,
        .axis_scaling_exponent = 0.667,
    });
    try std.testing.expectEqual(@as(f64, 0), state.retained_root_carbon_g_c_per_plant);
    try std.testing.expectEqual(@as(f64, 0), result.primary_root_axis_count_multiplier);
}

test "invalid input leaves retained root state unchanged" {
    var state: State = .{ .retained_root_carbon_g_c_per_plant = 7 };
    try std.testing.expectError(
        error.InvalidPrimaryRootAxisScalingInput,
        advance(&state, .{
            .total_root_carbon_g_c = std.math.nan(f64),
            .plant_population = 2,
            .biological_timestep_h = 1,
            .hourly_retention_fraction = 0.999992087,
            .axis_scaling_exponent = 0.667,
        }),
    );
    try std.testing.expectEqual(@as(f64, 7), state.retained_root_carbon_g_c_per_plant);
}

test "zero timestep and scaling exponent fail before mutation" {
    inline for (.{
        Inputs{ .total_root_carbon_g_c = 1, .plant_population = 1, .biological_timestep_h = 0, .hourly_retention_fraction = 0.999992087, .axis_scaling_exponent = 0.667 },
        Inputs{ .total_root_carbon_g_c = 1, .plant_population = 1, .biological_timestep_h = 1, .hourly_retention_fraction = 0.999992087, .axis_scaling_exponent = 0 },
    }) |inputs| {
        var state: State = .{ .retained_root_carbon_g_c_per_plant = 7 };
        try std.testing.expectError(error.InvalidPrimaryRootAxisScalingInput, advance(&state, inputs));
        try std.testing.expectEqual(@as(f64, 7), state.retained_root_carbon_g_c_per_plant);
    }
}
