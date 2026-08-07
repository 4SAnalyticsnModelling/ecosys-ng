const std = @import("std");

pub const Inputs = struct {
    current_canopy_water_volume_m3: f64,
    committed_canopy_water_volume_m3: f64,
    current_net_transpiration_m3_per_step: f64,
    current_root_uptake_m3_per_step: f64,
    horizontal_cell_area_m2: f64,
    water_storage_secant_resistance_mpa_per_m: f64,
    flux_secant_resistance_mpa_per_m: f64,
    current_canopy_water_potential_megapascal: f64,
    current_iteration: usize,
    maximum_iterations: usize,
    water_balance_tolerance_m_per_step: f64,
};

pub const Result = struct {
    canopy_storage_change_m3: f64,
    transpiration_minus_uptake_m3_per_step: f64,
    water_balance_residual_m_per_step: f64,
    potential_change_megapascal: f64,
    previous_net_transpiration_m3_per_step: f64,
    previous_root_uptake_m3_per_step: f64,
    previous_canopy_water_volume_m3: f64,
    previous_canopy_water_potential_megapascal: f64,
    next_canopy_water_potential_megapascal: f64,
    iteration_complete: bool,
};

/// UPTAKE.F 1190--1216 compatibility residual and half-step potential update.
/// Production Newton/Picard orchestration remains outside this kernel.
pub fn calculate(inputs: Inputs) !Result {
    try validate(inputs);
    const storage_change =
        inputs.current_canopy_water_volume_m3 -
        inputs.committed_canopy_water_volume_m3;
    const transpiration_minus_uptake =
        inputs.current_net_transpiration_m3_per_step -
        inputs.current_root_uptake_m3_per_step;
    const residual =
        (transpiration_minus_uptake - storage_change) /
        inputs.horizontal_cell_area_m2;
    const potential_change = @min(
        inputs.water_storage_secant_resistance_mpa_per_m,
        inputs.flux_secant_resistance_mpa_per_m,
    ) * residual;
    const next_potential = @min(
        0.5 * inputs.current_canopy_water_potential_megapascal,
        inputs.current_canopy_water_potential_megapascal +
            0.5 * potential_change,
    );
    const result = Result{
        .canopy_storage_change_m3 = storage_change,
        .transpiration_minus_uptake_m3_per_step = transpiration_minus_uptake,
        .water_balance_residual_m_per_step = residual,
        .potential_change_megapascal = potential_change,
        .previous_net_transpiration_m3_per_step = inputs.current_net_transpiration_m3_per_step,
        .previous_root_uptake_m3_per_step = inputs.current_root_uptake_m3_per_step,
        .previous_canopy_water_volume_m3 = inputs.current_canopy_water_volume_m3,
        .previous_canopy_water_potential_megapascal = inputs.current_canopy_water_potential_megapascal,
        .next_canopy_water_potential_megapascal = next_potential,
        .iteration_complete = inputs.current_iteration >= inputs.maximum_iterations or
            @abs(residual) < inputs.water_balance_tolerance_m_per_step,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (field.type == bool) continue;
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteCanopyWaterBalanceUpdateResult;
    }
    return result;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == usize) continue;
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidCanopyWaterBalanceUpdateInput;
    }
    if (inputs.current_canopy_water_volume_m3 < 0 or
        inputs.committed_canopy_water_volume_m3 < 0 or
        inputs.horizontal_cell_area_m2 <= 0 or
        inputs.water_storage_secant_resistance_mpa_per_m < 0 or
        inputs.flux_secant_resistance_mpa_per_m < 0 or
        inputs.current_iteration == 0 or
        inputs.maximum_iterations == 0 or
        inputs.water_balance_tolerance_m_per_step < 0)
        return error.InvalidCanopyWaterBalanceUpdateInput;
}

fn sourceInputs() Inputs {
    return .{
        .current_canopy_water_volume_m3 = 4,
        .committed_canopy_water_volume_m3 = 3,
        .current_net_transpiration_m3_per_step = -0.5,
        .current_root_uptake_m3_per_step = -2,
        .horizontal_cell_area_m2 = 2,
        .water_storage_secant_resistance_mpa_per_m = 4,
        .flux_secant_resistance_mpa_per_m = 3,
        .current_canopy_water_potential_megapascal = -1,
        .current_iteration = 2,
        .maximum_iterations = 200,
        .water_balance_tolerance_m_per_step = 1e-8,
    };
}

test "UPTAKE water balance residual and potential update preserve source order" {
    const inputs = sourceInputs();
    const result = try calculate(inputs);
    const storage_change = 4.0 - 3.0;
    const flux_difference = -0.5 - -2.0;
    const residual = (flux_difference - storage_change) / 2.0;
    const potential_change = @min(4.0, 3.0) * residual;
    const next_potential = @min(-0.5, -1.0 + 0.5 * potential_change);
    try std.testing.expectEqual(storage_change, result.canopy_storage_change_m3);
    try std.testing.expectEqual(flux_difference, result.transpiration_minus_uptake_m3_per_step);
    try std.testing.expectEqual(residual, result.water_balance_residual_m_per_step);
    try std.testing.expectEqual(potential_change, result.potential_change_megapascal);
    try std.testing.expectEqual(next_potential, result.next_canopy_water_potential_megapascal);
    try std.testing.expect(!result.iteration_complete);
}

test "source upper potential limiter prevents positive relaxation" {
    var inputs = sourceInputs();
    inputs.current_root_uptake_m3_per_step = -10;
    const result = try calculate(inputs);
    try std.testing.expectEqual(
        0.5 * inputs.current_canopy_water_potential_megapascal,
        result.next_canopy_water_potential_megapascal,
    );
}

test "runtime maximum and strict residual tolerance terminate iteration" {
    var inputs = sourceInputs();
    inputs.current_iteration = inputs.maximum_iterations;
    var result = try calculate(inputs);
    try std.testing.expect(result.iteration_complete);

    inputs = sourceInputs();
    inputs.current_canopy_water_volume_m3 = 4.5;
    result = try calculate(inputs);
    try std.testing.expectEqual(@as(f64, 0), result.water_balance_residual_m_per_step);
    try std.testing.expect(result.iteration_complete);
}
