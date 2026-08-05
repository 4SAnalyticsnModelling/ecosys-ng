const std = @import("std");

pub const Inputs = struct {
    uppermost_soil_bulk_density_megagrams_m3: f64,
    density_sentinel_threshold_megagrams_m3: f64,
    snowpack_heat_capacity_megajoules_k: f64,
    minimum_snowpack_heat_capacity_megajoules_k: f64,
    snow_surface_roughness_m: f64,
    soil_surface_roughness_m: f64,
    surface_slope: f64,
    cell_area_m2: f64,
    adjusted_water_table_depth_m: f64,
    cumulative_depth_above_uppermost_soil_layer_m: f64,
    cumulative_depth_to_uppermost_soil_layer_bottom_m: f64,
    uppermost_soil_layer_thickness_m: f64,
};

pub const Result = struct {
    selected_surface_roughness_m: f64,
    drainable_surface_water_volume_m3: f64,
    groundwater_storage_volume_m3: f64,
    uppermost_soil_layer_center_depth_m: f64,
    disturbance_flag_after_reset: u8,
};

/// HOUR1 lines 2367--2388. Preserves roughness selection, VOLWD/VOLWG,
/// uppermost-layer center depth, and final disturbance-flag reset order.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    const selected_surface_roughness_m =
        if (inputs.uppermost_soil_bulk_density_megagrams_m3 <
        inputs.density_sentinel_threshold_megagrams_m3 or
        inputs.snowpack_heat_capacity_megajoules_k >
            inputs.minimum_snowpack_heat_capacity_megajoules_k)
            inputs.snow_surface_roughness_m
        else
            inputs.soil_surface_roughness_m;
    const drainable_surface_water_volume_m3 = @max(
        inputs.snow_surface_roughness_m,
        0.112 * selected_surface_roughness_m +
            3.10 * selected_surface_roughness_m *
                selected_surface_roughness_m -
            0.012 * selected_surface_roughness_m * inputs.surface_slope,
    ) * inputs.cell_area_m2;
    const groundwater_storage_volume_m3 = @max(
        drainable_surface_water_volume_m3,
        -(inputs.adjusted_water_table_depth_m -
            inputs.cumulative_depth_above_uppermost_soil_layer_m) *
            inputs.cell_area_m2,
    );
    const uppermost_soil_layer_center_depth_m =
        inputs.cumulative_depth_to_uppermost_soil_layer_bottom_m -
        0.5 * inputs.uppermost_soil_layer_thickness_m;
    const disturbance_flag_after_reset: u8 = 0;
    return .{
        .selected_surface_roughness_m = selected_surface_roughness_m,
        .drainable_surface_water_volume_m3 = drainable_surface_water_volume_m3,
        .groundwater_storage_volume_m3 = groundwater_storage_volume_m3,
        .uppermost_soil_layer_center_depth_m = uppermost_soil_layer_center_depth_m,
        .disturbance_flag_after_reset = disturbance_flag_after_reset,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteSurfaceRoughnessInput;
    inline for (.{
        inputs.snowpack_heat_capacity_megajoules_k,
        inputs.minimum_snowpack_heat_capacity_megajoules_k,
        inputs.snow_surface_roughness_m,
        inputs.soil_surface_roughness_m,
        inputs.cell_area_m2,
        inputs.cumulative_depth_to_uppermost_soil_layer_bottom_m,
        inputs.uppermost_soil_layer_thickness_m,
    }) |value| if (value < 0)
        return error.InvalidSurfaceRoughnessInput;
    if (inputs.cell_area_m2 == 0)
        return error.InvalidSurfaceRoughnessInput;
}

test "snow selects snow roughness and preserves storage calculation order" {
    const result = try compute(.{
        .uppermost_soil_bulk_density_megagrams_m3 = 1,
        .density_sentinel_threshold_megagrams_m3 = 0,
        .snowpack_heat_capacity_megajoules_k = 2,
        .minimum_snowpack_heat_capacity_megajoules_k = 1,
        .snow_surface_roughness_m = 0.02,
        .soil_surface_roughness_m = 0.01,
        .surface_slope = 0.1,
        .cell_area_m2 = 100,
        .adjusted_water_table_depth_m = -0.5,
        .cumulative_depth_above_uppermost_soil_layer_m = 0,
        .cumulative_depth_to_uppermost_soil_layer_bottom_m = 0.2,
        .uppermost_soil_layer_thickness_m = 0.2,
    });
    try std.testing.expectEqual(@as(f64, 0.02), result.selected_surface_roughness_m);
    try std.testing.expect(result.drainable_surface_water_volume_m3 >= 2);
    try std.testing.expectEqual(@as(f64, 50), result.groundwater_storage_volume_m3);
    try std.testing.expectEqual(@as(f64, 0.1), result.uppermost_soil_layer_center_depth_m);
    try std.testing.expectEqual(@as(u8, 0), result.disturbance_flag_after_reset);
}

test "bare soil selects soil roughness" {
    const result = try compute(.{
        .uppermost_soil_bulk_density_megagrams_m3 = 1,
        .density_sentinel_threshold_megagrams_m3 = 0,
        .snowpack_heat_capacity_megajoules_k = 0,
        .minimum_snowpack_heat_capacity_megajoules_k = 1,
        .snow_surface_roughness_m = 0.02,
        .soil_surface_roughness_m = 0.01,
        .surface_slope = 0,
        .cell_area_m2 = 1,
        .adjusted_water_table_depth_m = 1,
        .cumulative_depth_above_uppermost_soil_layer_m = 0,
        .cumulative_depth_to_uppermost_soil_layer_bottom_m = 0.2,
        .uppermost_soil_layer_thickness_m = 0.2,
    });
    try std.testing.expectEqual(@as(f64, 0.01), result.selected_surface_roughness_m);
}
