const std = @import("std");

pub const Inputs = struct {
    current_canopy_water_potential_mpa: f64,
    previous_canopy_water_potential_mpa: f64,
    current_canopy_water_volume_m3: f64,
    previous_canopy_water_volume_m3: f64,
    current_net_transpiration_m3_per_step: f64,
    previous_net_transpiration_m3_per_step: f64,
    current_root_uptake_m3_per_step: f64,
    previous_root_uptake_m3_per_step: f64,
    horizontal_cell_area_m2: f64,
    significant_change_threshold: f64,
    fallback_resistance_mpa_per_m: f64,
};

pub const Result = struct {
    absolute_water_potential_change_mpa: f64,
    canopy_water_depth_change_m: f64,
    water_storage_secant_resistance_mpa_per_m: f64,
    maximum_flux_depth_change_m_per_step: f64,
    flux_secant_resistance_mpa_per_m: f64,
};

/// UPTAKE.F 1154--1170 compatibility secants. The production nonlinear owner
/// may use analytic Newton/Picard derivatives, but this preserves the legacy
/// diagnostic calculation and fallback exactly.
pub fn calculate(inputs: Inputs) !Result {
    try validate(inputs);
    const potential_change = @abs(
        inputs.current_canopy_water_potential_mpa -
            inputs.previous_canopy_water_potential_mpa,
    );
    const water_depth_change = @abs(
        inputs.current_canopy_water_volume_m3 -
            inputs.previous_canopy_water_volume_m3,
    ) / inputs.horizontal_cell_area_m2;
    const storage_resistance =
        if (water_depth_change > inputs.significant_change_threshold and
        potential_change > inputs.significant_change_threshold)
            potential_change / water_depth_change
        else
            inputs.fallback_resistance_mpa_per_m;
    const flux_depth_change = @max(
        @abs(inputs.current_net_transpiration_m3_per_step -
            inputs.previous_net_transpiration_m3_per_step),
        @abs(inputs.current_root_uptake_m3_per_step -
            inputs.previous_root_uptake_m3_per_step),
    ) / inputs.horizontal_cell_area_m2;
    const flux_resistance =
        if (flux_depth_change > inputs.significant_change_threshold and
        potential_change > inputs.significant_change_threshold)
            potential_change / flux_depth_change
        else
            inputs.fallback_resistance_mpa_per_m;
    const result = Result{
        .absolute_water_potential_change_mpa = potential_change,
        .canopy_water_depth_change_m = water_depth_change,
        .water_storage_secant_resistance_mpa_per_m = storage_resistance,
        .maximum_flux_depth_change_m_per_step = flux_depth_change,
        .flux_secant_resistance_mpa_per_m = flux_resistance,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteCanopyWaterBalanceSecantResult;
    return result;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidCanopyWaterBalanceSecantInput;
    if (inputs.current_canopy_water_volume_m3 < 0 or
        inputs.previous_canopy_water_volume_m3 < 0 or
        inputs.horizontal_cell_area_m2 <= 0 or
        inputs.significant_change_threshold < 0 or
        inputs.fallback_resistance_mpa_per_m <= 0)
        return error.InvalidCanopyWaterBalanceSecantInput;
}

test "UPTAKE canopy water secants preserve source order" {
    const result = try calculate(.{
        .current_canopy_water_potential_mpa = -1.5,
        .previous_canopy_water_potential_mpa = -1,
        .current_canopy_water_volume_m3 = 4,
        .previous_canopy_water_volume_m3 = 2,
        .current_net_transpiration_m3_per_step = -0.6,
        .previous_net_transpiration_m3_per_step = -0.2,
        .current_root_uptake_m3_per_step = -0.8,
        .previous_root_uptake_m3_per_step = -0.1,
        .horizontal_cell_area_m2 = 2,
        .significant_change_threshold = 1e-12,
        .fallback_resistance_mpa_per_m = 1000,
    });
    try std.testing.expectEqual(@as(f64, 0.5), result.absolute_water_potential_change_mpa);
    try std.testing.expectEqual(@as(f64, 1), result.canopy_water_depth_change_m);
    try std.testing.expectEqual(@as(f64, 0.5), result.water_storage_secant_resistance_mpa_per_m);
    try std.testing.expectApproxEqAbs(@as(f64, 0.35), result.maximum_flux_depth_change_m_per_step, 1e-16);
    try std.testing.expectApproxEqAbs(0.5 / 0.35, result.flux_secant_resistance_mpa_per_m, 1e-15);
}

test "insignificant changes retain runtime compatibility fallback" {
    const result = try calculate(.{
        .current_canopy_water_potential_mpa = -1,
        .previous_canopy_water_potential_mpa = -1,
        .current_canopy_water_volume_m3 = 2,
        .previous_canopy_water_volume_m3 = 2,
        .current_net_transpiration_m3_per_step = -0.2,
        .previous_net_transpiration_m3_per_step = -0.2,
        .current_root_uptake_m3_per_step = -0.1,
        .previous_root_uptake_m3_per_step = -0.1,
        .horizontal_cell_area_m2 = 2,
        .significant_change_threshold = 1e-12,
        .fallback_resistance_mpa_per_m = 1000,
    });
    try std.testing.expectEqual(@as(f64, 1000), result.water_storage_secant_resistance_mpa_per_m);
    try std.testing.expectEqual(@as(f64, 1000), result.flux_secant_resistance_mpa_per_m);
}

test "zero cell area fails explicitly" {
    var inputs = Inputs{
        .current_canopy_water_potential_mpa = -1,
        .previous_canopy_water_potential_mpa = -1,
        .current_canopy_water_volume_m3 = 2,
        .previous_canopy_water_volume_m3 = 2,
        .current_net_transpiration_m3_per_step = 0,
        .previous_net_transpiration_m3_per_step = 0,
        .current_root_uptake_m3_per_step = 0,
        .previous_root_uptake_m3_per_step = 0,
        .horizontal_cell_area_m2 = 1,
        .significant_change_threshold = 0,
        .fallback_resistance_mpa_per_m = 1000,
    };
    inputs.horizontal_cell_area_m2 = 0;
    try std.testing.expectError(
        error.InvalidCanopyWaterBalanceSecantInput,
        calculate(inputs),
    );
}
