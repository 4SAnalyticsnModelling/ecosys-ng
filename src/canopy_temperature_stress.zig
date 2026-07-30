const std = @import("std");

pub const Inputs = struct {
    canopy_surface_temperature_c: f64,
    chilling_threshold_c: f64,
    accumulated_chilling_h: f64,
    accumulated_heating_degree_h: f64,
    biological_timestep_h_per_step: f64,
    maximum_chilling_h: f64,
    chilling_accumulation_rate: f64,
    heating_threshold_c: f64,
    heating_recovery_rate_per_h: f64,
};

pub const Result = struct {
    accumulated_chilling_h: f64,
    accumulated_heating_degree_h: f64,
};

/// UPTAKE.F 1601--1610. Updates canopy chilling and heating stress in the
/// exact source branch order using runtime thresholds and rates.
pub fn calculate(inputs: Inputs) !Result {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidCanopyTemperatureStressInput;
    if (inputs.accumulated_chilling_h < 0 or
        inputs.accumulated_heating_degree_h < 0 or
        inputs.biological_timestep_h_per_step < 0 or
        inputs.maximum_chilling_h < 0 or
        inputs.chilling_accumulation_rate < 0 or
        inputs.heating_recovery_rate_per_h < 0)
        return error.InvalidCanopyTemperatureStressInput;
    const chilling =
        if (inputs.canopy_surface_temperature_c < inputs.chilling_threshold_c)
            @min(
                inputs.maximum_chilling_h,
                inputs.accumulated_chilling_h +
                    inputs.chilling_accumulation_rate *
                        inputs.biological_timestep_h_per_step,
            )
        else
            @max(
                0,
                inputs.accumulated_chilling_h -
                    inputs.chilling_accumulation_rate *
                        inputs.biological_timestep_h_per_step,
            );
    const heating =
        if (inputs.canopy_surface_temperature_c > inputs.heating_threshold_c)
            inputs.accumulated_heating_degree_h +
                (inputs.canopy_surface_temperature_c -
                    inputs.heating_threshold_c) *
                    inputs.biological_timestep_h_per_step
        else
            @max(
                0,
                inputs.accumulated_heating_degree_h -
                    inputs.heating_recovery_rate_per_h *
                        inputs.biological_timestep_h_per_step,
            );
    if (!std.math.isFinite(chilling) or !std.math.isFinite(heating))
        return error.NonFiniteCanopyTemperatureStressResult;
    return .{
        .accumulated_chilling_h = chilling,
        .accumulated_heating_degree_h = heating,
    };
}

test "cold canopy accumulates capped chilling and recovers heat" {
    const result = try calculate(.{
        .canopy_surface_temperature_c = 2,
        .chilling_threshold_c = 5,
        .accumulated_chilling_h = 23.5,
        .accumulated_heating_degree_h = 1,
        .biological_timestep_h_per_step = 1,
        .maximum_chilling_h = 24,
        .chilling_accumulation_rate = 1,
        .heating_threshold_c = 60,
        .heating_recovery_rate_per_h = 0.02,
    });
    try std.testing.expectEqual(@as(f64, 24), result.accumulated_chilling_h);
    try std.testing.expectEqual(@as(f64, 0.98), result.accumulated_heating_degree_h);
}

test "hot canopy releases chilling and accumulates degree hours" {
    const result = try calculate(.{
        .canopy_surface_temperature_c = 65,
        .chilling_threshold_c = 5,
        .accumulated_chilling_h = 2,
        .accumulated_heating_degree_h = 1,
        .biological_timestep_h_per_step = 0.5,
        .maximum_chilling_h = 24,
        .chilling_accumulation_rate = 1,
        .heating_threshold_c = 60,
        .heating_recovery_rate_per_h = 0.02,
    });
    try std.testing.expectEqual(@as(f64, 1.5), result.accumulated_chilling_h);
    try std.testing.expectEqual(@as(f64, 3.5), result.accumulated_heating_degree_h);
}

test "stress recovery retains zero floors" {
    const result = try calculate(.{
        .canopy_surface_temperature_c = 20,
        .chilling_threshold_c = 5,
        .accumulated_chilling_h = 0.2,
        .accumulated_heating_degree_h = 0.001,
        .biological_timestep_h_per_step = 1,
        .maximum_chilling_h = 24,
        .chilling_accumulation_rate = 1,
        .heating_threshold_c = 60,
        .heating_recovery_rate_per_h = 0.02,
    });
    try std.testing.expectEqual(@as(f64, 0), result.accumulated_chilling_h);
    try std.testing.expectEqual(@as(f64, 0), result.accumulated_heating_degree_h);
}
