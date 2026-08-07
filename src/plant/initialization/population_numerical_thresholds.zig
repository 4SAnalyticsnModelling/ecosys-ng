const std = @import("std");

pub const Inputs = struct {
    plant_population_count: f64,
    cell_area_m2: f64,
    primary_threshold_per_plant: f64,
    secondary_threshold_per_plant: f64,
};

pub const Thresholds = struct {
    population_scaled_primary: f64,
    area_normalized_primary_per_m2: f64,
    population_scaled_secondary: f64,
};

pub const CalculationError = error{
    NonFiniteInput,
    NegativePopulation,
    NonPositiveCellArea,
    NegativeThreshold,
    NonFiniteResult,
};

/// Translates `startq.f` lines 886--888 for one plant species.
pub fn calculate(inputs: Inputs) CalculationError!Thresholds {
    inline for (std.meta.fields(Inputs)) |field| {
        if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteInput;
    }
    if (inputs.plant_population_count < 0.0) return error.NegativePopulation;
    if (inputs.cell_area_m2 <= 0.0) return error.NonPositiveCellArea;
    if (inputs.primary_threshold_per_plant < 0.0 or
        inputs.secondary_threshold_per_plant < 0.0)
    {
        return error.NegativeThreshold;
    }

    const population_scaled_primary =
        inputs.primary_threshold_per_plant * inputs.plant_population_count;
    const area_normalized_primary_per_m2 =
        inputs.primary_threshold_per_plant * inputs.plant_population_count /
        inputs.cell_area_m2;
    const population_scaled_secondary =
        inputs.secondary_threshold_per_plant * inputs.plant_population_count;
    const result = Thresholds{
        .population_scaled_primary = population_scaled_primary,
        .area_normalized_primary_per_m2 = area_normalized_primary_per_m2,
        .population_scaled_secondary = population_scaled_secondary,
    };
    inline for (std.meta.fields(Thresholds)) |field| {
        if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteResult;
    }
    return result;
}

test "thresholds scale with population and cell area in STARTQ order" {
    const result = try calculate(.{
        .plant_population_count = 200.0,
        .cell_area_m2 = 50.0,
        .primary_threshold_per_plant = 1.0e-12,
        .secondary_threshold_per_plant = 1.0e-8,
    });
    try std.testing.expectEqual(@as(f64, 2.0e-10), result.population_scaled_primary);
    try std.testing.expectEqual(@as(f64, 4.0e-12), result.area_normalized_primary_per_m2);
    try std.testing.expectEqual(@as(f64, 2.0e-6), result.population_scaled_secondary);
}

test "zero population produces finite zero thresholds" {
    const result = try calculate(.{
        .plant_population_count = 0.0,
        .cell_area_m2 = 10.0,
        .primary_threshold_per_plant = 1.0e-12,
        .secondary_threshold_per_plant = 1.0e-8,
    });
    try std.testing.expectEqual(std.mem.zeroes(Thresholds), result);
}
