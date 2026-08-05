const std = @import("std");

pub const ErosionLayerReset = enum { disabled, enabled };

pub const Inputs = struct {
    first_output_layer: usize, // NUI
    last_layer: *usize, // NL; slices include surface index zero
    minimum_layer_thickness_m: f64, // DLYRM
    erosion_layer_reset: ErosionLayerReset, // IERSNG > 0
    cell_area_m2: []const f64, // AREA(3,...)
    liquid_water_m3: []const f64, // VOLW
    ice_water_equivalent_m3: []const f64, // VOLI
    surface_residual_water_m3: f64, // VOLWRX
    held_liquid_water_m3: []const f64, // VOLWH
    held_ice_water_equivalent_m3: []const f64, // VOLIH
    macropore_air_volume_m3: []const f64, // VOLAH
    micropore_fraction: []const f64, // FMPR
    layer_thickness_m: []f64, // DLYR(3,...)
    cumulative_bottom_depth_m: []f64, // CDPTH
    organic_carbon_g_c: []f64, // ORGC
    organic_nitrogen_g_n: []f64, // ORGN
    charcoal_organic_carbon_g_c: []const f64, // ORGCC
    previous_charcoal_organic_carbon_g_c: []const f64, // ORGCCX
    liquid_water_fraction: []f64, // THETWZ
    ice_water_equivalent_fraction: []f64, // THETIZ
    charcoal_carbon_change_g_c: []f64, // DORGCC
};

fn finite(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 13030--13090. The loop limit is captured
/// before erosion can reduce NL, matching Fortran DO-loop semantics.
pub fn update(allocator: std.mem.Allocator, inputs: Inputs) !void {
    const count = inputs.cell_area_m2.len;
    const original_last_layer = inputs.last_layer.*;
    if (count == 0 or original_last_layer >= count or inputs.first_output_layer > original_last_layer) return error.SoilOutputStateDimensionMismatch;
    inline for (.{ inputs.liquid_water_m3, inputs.ice_water_equivalent_m3, inputs.held_liquid_water_m3, inputs.held_ice_water_equivalent_m3, inputs.macropore_air_volume_m3, inputs.micropore_fraction, inputs.layer_thickness_m, inputs.cumulative_bottom_depth_m, inputs.organic_carbon_g_c, inputs.organic_nitrogen_g_n, inputs.charcoal_organic_carbon_g_c, inputs.previous_charcoal_organic_carbon_g_c, inputs.liquid_water_fraction, inputs.ice_water_equivalent_fraction, inputs.charcoal_carbon_change_g_c }) |values| {
        if (values.len != count) return error.SoilOutputStateDimensionMismatch;
        if (!finite(values)) return error.InvalidSoilOutputStateInput;
    }
    if (!std.math.isFinite(inputs.minimum_layer_thickness_m) or !std.math.isFinite(inputs.surface_residual_water_m3) or inputs.minimum_layer_thickness_m < 0) return error.InvalidSoilOutputStateInput;
    for (inputs.cell_area_m2) |area_m2| if (area_m2 <= 0) return error.InvalidSoilOutputStateArea;

    const staged_thickness = try allocator.dupe(f64, inputs.layer_thickness_m);
    defer allocator.free(staged_thickness);
    const staged_depth = try allocator.dupe(f64, inputs.cumulative_bottom_depth_m);
    defer allocator.free(staged_depth);
    const staged_carbon = try allocator.dupe(f64, inputs.organic_carbon_g_c);
    defer allocator.free(staged_carbon);
    const staged_nitrogen = try allocator.dupe(f64, inputs.organic_nitrogen_g_n);
    defer allocator.free(staged_nitrogen);
    const staged_water = try allocator.dupe(f64, inputs.liquid_water_fraction);
    defer allocator.free(staged_water);
    const staged_ice = try allocator.dupe(f64, inputs.ice_water_equivalent_fraction);
    defer allocator.free(staged_ice);
    const staged_charcoal_change = try allocator.dupe(f64, inputs.charcoal_carbon_change_g_c);
    defer allocator.free(staged_charcoal_change);
    var staged_last_layer = original_last_layer;

    staged_water[0] = @max(0.0, (inputs.liquid_water_m3[0] - inputs.surface_residual_water_m3) / inputs.cell_area_m2[0]);
    staged_ice[0] = @max(0.0, (inputs.ice_water_equivalent_m3[0] - inputs.surface_residual_water_m3) / inputs.cell_area_m2[0]);
    staged_charcoal_change[0] = inputs.charcoal_organic_carbon_g_c[0] - inputs.previous_charcoal_organic_carbon_g_c[0];
    inline for (.{ staged_water[0], staged_ice[0], staged_charcoal_change[0] }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSoilOutputStateResult;

    for (inputs.first_output_layer..original_last_layer + 1) |layer| {
        if (staged_thickness[layer] <= inputs.minimum_layer_thickness_m) staged_thickness[layer] = @max(0.0, (inputs.liquid_water_m3[layer] + inputs.ice_water_equivalent_m3[layer]) / inputs.cell_area_m2[layer]);
        const micropore_volume_m3 = inputs.cell_area_m2[layer] * staged_thickness[layer] * inputs.micropore_fraction[layer];
        const total_effective_volume_m3 = micropore_volume_m3 + inputs.macropore_air_volume_m3[layer];
        if (!std.math.isFinite(micropore_volume_m3) or !std.math.isFinite(total_effective_volume_m3) or total_effective_volume_m3 == 0) return error.InvalidSoilOutputStateVolume;
        staged_water[layer] = (inputs.liquid_water_m3[layer] + @min(inputs.macropore_air_volume_m3[layer], inputs.held_liquid_water_m3[layer])) / total_effective_volume_m3;
        staged_ice[layer] = (inputs.ice_water_equivalent_m3[layer] + @min(inputs.macropore_air_volume_m3[layer], inputs.held_ice_water_equivalent_m3[layer])) / total_effective_volume_m3;
        if (!std.math.isFinite(staged_thickness[layer]) or !std.math.isFinite(staged_water[layer]) or !std.math.isFinite(staged_ice[layer])) return error.NonFiniteSoilOutputStateResult;

        if (inputs.erosion_layer_reset == .enabled and staged_last_layer > 0 and layer == staged_last_layer - 1 and staged_depth[staged_last_layer] - staged_depth[layer] <= inputs.minimum_layer_thickness_m) {
            staged_depth[layer] = staged_depth[layer] + staged_thickness[staged_last_layer];
            staged_depth[staged_last_layer] = staged_depth[layer];
            staged_thickness[staged_last_layer] = 0;
            staged_carbon[staged_last_layer] = 0;
            staged_nitrogen[staged_last_layer] = 0;
            staged_last_layer = layer;
            if (!std.math.isFinite(staged_depth[layer])) return error.NonFiniteSoilOutputStateResult;
        }
        staged_charcoal_change[layer] = inputs.charcoal_organic_carbon_g_c[layer] - inputs.previous_charcoal_organic_carbon_g_c[layer];
        if (!std.math.isFinite(staged_charcoal_change[layer])) return error.NonFiniteSoilOutputStateResult;
    }

    @memcpy(inputs.layer_thickness_m, staged_thickness);
    @memcpy(inputs.cumulative_bottom_depth_m, staged_depth);
    @memcpy(inputs.organic_carbon_g_c, staged_carbon);
    @memcpy(inputs.organic_nitrogen_g_n, staged_nitrogen);
    @memcpy(inputs.liquid_water_fraction, staged_water);
    @memcpy(inputs.ice_water_equivalent_fraction, staged_ice);
    @memcpy(inputs.charcoal_carbon_change_g_c, staged_charcoal_change);
    inputs.last_layer.* = staged_last_layer;
}

fn baseInputs(last_layer: *usize, arrays: *[15][3]f64) Inputs {
    return .{
        .first_output_layer = 1,
        .last_layer = last_layer,
        .minimum_layer_thickness_m = 0.05,
        .erosion_layer_reset = .disabled,
        .cell_area_m2 = &arrays[0],
        .liquid_water_m3 = &arrays[1],
        .ice_water_equivalent_m3 = &arrays[2],
        .surface_residual_water_m3 = 1,
        .held_liquid_water_m3 = &arrays[3],
        .held_ice_water_equivalent_m3 = &arrays[4],
        .macropore_air_volume_m3 = &arrays[5],
        .micropore_fraction = &arrays[6],
        .layer_thickness_m = &arrays[7],
        .cumulative_bottom_depth_m = &arrays[8],
        .organic_carbon_g_c = &arrays[9],
        .organic_nitrogen_g_n = &arrays[10],
        .charcoal_organic_carbon_g_c = &arrays[11],
        .previous_charcoal_organic_carbon_g_c = &arrays[12],
        .liquid_water_fraction = &arrays[13],
        .ice_water_equivalent_fraction = &arrays[14],
        .charcoal_carbon_change_g_c = &arrays[3],
    };
}

test "REDIST soil output state preserves surface and runtime layer equations" {
    var last_layer: usize = 2;
    var arrays: [15][3]f64 = @splat(@splat(0));
    arrays[0] = .{ 2, 10, 10 };
    arrays[1] = .{ 5, 2, 3 };
    arrays[2] = .{ 3, 1, 1 };
    arrays[3] = .{ 0, 0.5, 0.5 };
    arrays[4] = .{ 0, 0.25, 0.25 };
    arrays[5] = .{ 0, 1, 1 };
    arrays[6] = .{ 0, 0.5, 0.5 };
    arrays[7] = .{ 0, 0.2, 0.4 };
    arrays[8] = .{ 0, 0.2, 0.6 };
    arrays[11] = .{ 5, 8, 12 };
    arrays[12] = .{ 2, 3, 4 };
    var charcoal_change = [_]f64{ 0, 0, 0 };
    var inputs = baseInputs(&last_layer, &arrays);
    inputs.charcoal_carbon_change_g_c = &charcoal_change;
    try update(std.testing.allocator, inputs);
    try std.testing.expectEqual(@as(f64, 2), arrays[13][0]);
    try std.testing.expectEqual(@as(f64, 1), arrays[14][0]);
    try std.testing.expectApproxEqAbs(@as(f64, 1.25), arrays[13][1], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.625), arrays[14][1], 1.0e-12);
    try std.testing.expectEqual(@as(f64, 3), charcoal_change[0]);
    try std.testing.expectEqual(@as(f64, 5), charcoal_change[1]);
    try std.testing.expectEqual(@as(f64, 8), charcoal_change[2]);
}

test "REDIST erosion collapse retains original Fortran DO upper bound" {
    var last_layer: usize = 2;
    var arrays: [15][3]f64 = @splat(@splat(1));
    arrays[0] = .{ 1, 10, 10 };
    arrays[1] = .{ 2, 1, 1 };
    arrays[2] = .{ 2, 0, 0 };
    arrays[5] = .{ 1, 1, 1 };
    arrays[6] = .{ 1, 1, 1 };
    arrays[7] = .{ 1, 0.2, 0.01 };
    arrays[8] = .{ 0, 0.2, 0.21 };
    arrays[9] = .{ 1, 2, 3 };
    arrays[10] = .{ 1, 2, 3 };
    arrays[11] = .{ 3, 4, 5 };
    arrays[12] = .{ 1, 1, 1 };
    var charcoal_change = [_]f64{ 0, 0, 0 };
    var inputs = baseInputs(&last_layer, &arrays);
    inputs.erosion_layer_reset = .enabled;
    inputs.charcoal_carbon_change_g_c = &charcoal_change;
    try update(std.testing.allocator, inputs);
    try std.testing.expectEqual(@as(usize, 1), last_layer);
    try std.testing.expectEqual(@as(f64, 0), arrays[9][2]);
    try std.testing.expectEqual(@as(f64, 0), arrays[10][2]);
    try std.testing.expectEqual(@as(f64, 4), charcoal_change[2]);
}
