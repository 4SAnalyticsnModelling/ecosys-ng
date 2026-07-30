const std = @import("std");

pub const Inputs = struct {
    surface_water_volume_m3: f64,
    surface_ice_volume_m3: f64,
    litter_water_capacity_m3: f64,
    cell_area_m2: f64,
    surface_roughness_m: f64,
    negligible_volume_m3: f64,
};

pub const Result = struct {
    excess_water_depth_m: f64,
    excess_ice_depth_m: f64,
    depth_above_roughness_m: f64,
};

/// HOUR1 lines 963--974. Partitions litter retention proportionally between
/// surface water and ice before converting excess volumes to depths.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    const total_water_ice_volume_m3 =
        inputs.surface_water_volume_m3 + inputs.surface_ice_volume_m3;
    var excess_water_depth_m: f64 = 0;
    var excess_ice_depth_m: f64 = 0;
    if (total_water_ice_volume_m3 > inputs.negligible_volume_m3) {
        const retained_water_volume_m3 =
            inputs.surface_water_volume_m3 / total_water_ice_volume_m3 *
            inputs.litter_water_capacity_m3;
        const retained_ice_volume_m3 =
            inputs.surface_ice_volume_m3 / total_water_ice_volume_m3 *
            inputs.litter_water_capacity_m3;
        excess_water_depth_m = @max(
            0,
            inputs.surface_water_volume_m3 - retained_water_volume_m3,
        ) / inputs.cell_area_m2;
        excess_ice_depth_m = @max(
            0,
            inputs.surface_ice_volume_m3 - retained_ice_volume_m3,
        ) / inputs.cell_area_m2;
    }
    const depth_above_roughness_m = @max(
        0,
        excess_water_depth_m + excess_ice_depth_m -
            0.5 * inputs.surface_roughness_m,
    );
    if (!std.math.isFinite(depth_above_roughness_m))
        return error.NonFiniteSurfaceRadiationDepth;
    return .{
        .excess_water_depth_m = excess_water_depth_m,
        .excess_ice_depth_m = excess_ice_depth_m,
        .depth_above_roughness_m = depth_above_roughness_m,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSurfaceRadiationDepthInput;
        if (value < 0)
            return error.InvalidSurfaceRadiationDepthInput;
    }
    if (inputs.cell_area_m2 == 0)
        return error.InvalidSurfaceRadiationCellArea;
}

test "surface water and ice share litter capacity proportionally" {
    const result = try compute(.{
        .surface_water_volume_m3 = 6,
        .surface_ice_volume_m3 = 4,
        .litter_water_capacity_m3 = 5,
        .cell_area_m2 = 10,
        .surface_roughness_m = 0.2,
        .negligible_volume_m3 = 1.0e-12,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), result.excess_water_depth_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), result.excess_ice_depth_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), result.depth_above_roughness_m, 1e-15);
}

test "negligible combined volume produces zero interception depth" {
    const result = try compute(.{
        .surface_water_volume_m3 = 0,
        .surface_ice_volume_m3 = 0,
        .litter_water_capacity_m3 = 2,
        .cell_area_m2 = 10,
        .surface_roughness_m = 0.1,
        .negligible_volume_m3 = 1.0e-12,
    });
    try std.testing.expectEqual(@as(f64, 0), result.excess_water_depth_m);
    try std.testing.expectEqual(@as(f64, 0), result.excess_ice_depth_m);
    try std.testing.expectEqual(@as(f64, 0), result.depth_above_roughness_m);
}

test "zero cell area fails before division" {
    try std.testing.expectError(error.InvalidSurfaceRadiationCellArea, compute(.{
        .surface_water_volume_m3 = 1,
        .surface_ice_volume_m3 = 0,
        .litter_water_capacity_m3 = 0,
        .cell_area_m2 = 0,
        .surface_roughness_m = 0,
        .negligible_volume_m3 = 0,
    }));
}
