const std = @import("std");

pub const State = struct {
    cumulative_depth_m: []f64, // CDPTH
    previous_cumulative_depth_m: []const f64, // CDPTHX
    freeze_excluded_cumulative_depth_m: []f64, // CDPTHY
    previous_layer_volume_m3: []const f64, // VOLX
    working_layer_volume_m3: []f64, // VOLY
    pond_profile_flag_by_layer: []u8, // IFLGK
};

pub const Inputs = struct {
    bulk_density_megagrams_per_m3: []const f64, // BKDS
    water_m3: []const f64, // VOLW
    ice_m3: []const f64, // VOLI
    area_m2: []const f64,
    boundary_change_m_by_layer_channel: []const f64, // DDLYR(LX,1..4)
    pond_profile_flag: u8, // IFLGJ
    zero_tolerance: f64,
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8022--8141. Channels are stored layer-major
/// in legacy order: pond, freeze-thaw, erosion, SOC.
pub fn apply(
    disturbance_mode: i8,
    top_layer: usize,
    bottom_layer: usize,
    state: State,
    inputs: Inputs,
) !void {
    const len = state.cumulative_depth_m.len;
    if (len == 0 or top_layer == 0 or top_layer > bottom_layer or bottom_layer >= len or
        state.previous_cumulative_depth_m.len != len or state.freeze_excluded_cumulative_depth_m.len != len or
        state.previous_layer_volume_m3.len != len or state.working_layer_volume_m3.len != len or
        state.pond_profile_flag_by_layer.len != len or inputs.bulk_density_megagrams_per_m3.len != len or
        inputs.water_m3.len != len or inputs.ice_m3.len != len or inputs.area_m2.len != len or
        inputs.boundary_change_m_by_layer_channel.len != len * 4)
        return error.SoilProfileDepthApplyDimensionMismatch;
    inline for (.{ state.cumulative_depth_m, state.previous_cumulative_depth_m, state.freeze_excluded_cumulative_depth_m, state.previous_layer_volume_m3, state.working_layer_volume_m3, inputs.bulk_density_megagrams_per_m3, inputs.water_m3, inputs.ice_m3, inputs.area_m2, inputs.boundary_change_m_by_layer_channel }) |values|
        if (!finiteSlice(values)) return error.InvalidSoilProfileDepthApplyInput;
    if (!std.math.isFinite(inputs.zero_tolerance) or inputs.zero_tolerance < 0)
        return error.InvalidSoilProfileDepthApplyInput;
    if (disturbance_mode < 0) return;

    var layer = bottom_layer;
    while (true) : (layer -= 1) {
        const is_pond = inputs.bulk_density_megagrams_per_m3[layer] <= inputs.zero_tolerance;
        if (is_pond and inputs.area_m2[layer] <= 0) return error.InvalidPondGeometry;
        for (0..4) |channel| {
            const boundary_change = inputs.boundary_change_m_by_layer_channel[layer * 4 + channel];
            if (is_pond) {
                state.cumulative_depth_m[layer] = state.cumulative_depth_m[layer] + boundary_change;
                if (layer == top_layer) {
                    state.cumulative_depth_m[layer - 1] = state.cumulative_depth_m[layer] -
                        (inputs.water_m3[layer] + inputs.ice_m3[layer]) / inputs.area_m2[layer];
                }
                if (channel != 1) {
                    state.freeze_excluded_cumulative_depth_m[layer] =
                        state.freeze_excluded_cumulative_depth_m[layer] + boundary_change;
                    if (layer == top_layer) {
                        state.freeze_excluded_cumulative_depth_m[layer - 1] =
                            state.freeze_excluded_cumulative_depth_m[layer] -
                            (inputs.water_m3[layer] + inputs.ice_m3[layer]) / inputs.area_m2[layer];
                    }
                } else {
                    state.freeze_excluded_cumulative_depth_m[layer] = state.previous_cumulative_depth_m[layer];
                }
            } else {
                state.cumulative_depth_m[layer] = state.cumulative_depth_m[layer] + boundary_change;
                if (layer == top_layer) {
                    state.cumulative_depth_m[layer - 1] = state.cumulative_depth_m[layer - 1] +
                        inputs.boundary_change_m_by_layer_channel[(layer - 1) * 4 + channel];
                }
                if (channel != 1) {
                    state.freeze_excluded_cumulative_depth_m[layer] =
                        state.freeze_excluded_cumulative_depth_m[layer] + boundary_change;
                    if (layer == top_layer) {
                        state.freeze_excluded_cumulative_depth_m[layer - 1] =
                            state.freeze_excluded_cumulative_depth_m[layer - 1] +
                            inputs.boundary_change_m_by_layer_channel[(layer - 1) * 4 + channel];
                    }
                } else {
                    state.freeze_excluded_cumulative_depth_m[layer] = state.previous_cumulative_depth_m[layer];
                }
            }
            if (!std.math.isFinite(state.cumulative_depth_m[layer]) or
                !std.math.isFinite(state.freeze_excluded_cumulative_depth_m[layer]))
                return error.NonFiniteSoilProfileDepthApply;
        }
        state.working_layer_volume_m3[layer] = state.previous_layer_volume_m3[layer];
        state.pond_profile_flag_by_layer[layer] = inputs.pond_profile_flag;
        if (layer == top_layer) break;
    }
    state.working_layer_volume_m3[0] = if (inputs.bulk_density_megagrams_per_m3[top_layer] <= inputs.zero_tolerance)
        inputs.water_m3[0] + inputs.ice_m3[0]
    else
        state.previous_layer_volume_m3[0];
    if (!finiteSlice(state.cumulative_depth_m) or
        !finiteSlice(state.freeze_excluded_cumulative_depth_m) or
        !finiteSlice(state.working_layer_volume_m3))
        return error.NonFiniteSoilProfileDepthApply;
}

test "REDIST soil depth applies four soil channels and resets CDPTHY at freeze" {
    var cumulative = [_]f64{ 0, 10, 0 };
    const previous = [_]f64{ 0, 10, 0 };
    var excluded = [_]f64{ 0, 10, 0 };
    const old_volume = [_]f64{ 3, 4, 0 };
    var working = [_]f64{ 0, 0, 0 };
    var flags = [_]u8{ 0, 0, 0 };
    const density = [_]f64{ 0, 1, 0 };
    const zeros = [_]f64{ 0, 0, 0 };
    const ones = [_]f64{ 1, 1, 1 };
    const changes = [_]f64{ 1, 2, 3, 4, 1, 2, 3, 4, 0, 0, 0, 0 };
    try apply(0, 1, 1, .{ .cumulative_depth_m = &cumulative, .previous_cumulative_depth_m = &previous, .freeze_excluded_cumulative_depth_m = &excluded, .previous_layer_volume_m3 = &old_volume, .working_layer_volume_m3 = &working, .pond_profile_flag_by_layer = &flags }, .{ .bulk_density_megagrams_per_m3 = &density, .water_m3 = &zeros, .ice_m3 = &zeros, .area_m2 = &ones, .boundary_change_m_by_layer_channel = &changes, .pond_profile_flag = 1, .zero_tolerance = 0 });
    try std.testing.expectEqual(@as(f64, 20), cumulative[1]);
    try std.testing.expectEqual(@as(f64, 17), excluded[1]);
    try std.testing.expectEqual(@as(f64, 10), cumulative[0]);
    try std.testing.expectEqual(@as(f64, 8), excluded[0]);
    try std.testing.expectEqual(@as(f64, 4), working[1]);
    try std.testing.expectEqual(@as(f64, 3), working[0]);
}

test "REDIST pond depth resets surface boundary from water depth each channel" {
    var cumulative = [_]f64{ 0, 10, 0 };
    const previous = [_]f64{ 0, 10, 0 };
    var excluded = [_]f64{ 0, 10, 0 };
    const volume = [_]f64{ 0, 4, 0 };
    var working = [_]f64{ 0, 0, 0 };
    var flags = [_]u8{ 0, 0, 0 };
    const density = [_]f64{ 0, 0, 0 };
    const water = [_]f64{ 2, 2, 0 };
    const ice = [_]f64{ 1, 1, 0 };
    const area = [_]f64{ 1, 1, 1 };
    const changes = [_]f64{ 0, 0, 0, 0, 1, 2, 3, 4, 0, 0, 0, 0 };
    try apply(0, 1, 1, .{ .cumulative_depth_m = &cumulative, .previous_cumulative_depth_m = &previous, .freeze_excluded_cumulative_depth_m = &excluded, .previous_layer_volume_m3 = &volume, .working_layer_volume_m3 = &working, .pond_profile_flag_by_layer = &flags }, .{ .bulk_density_megagrams_per_m3 = &density, .water_m3 = &water, .ice_m3 = &ice, .area_m2 = &area, .boundary_change_m_by_layer_channel = &changes, .pond_profile_flag = 1, .zero_tolerance = 0 });
    try std.testing.expectEqual(@as(f64, 20), cumulative[1]);
    try std.testing.expectEqual(@as(f64, 17), cumulative[0]);
    try std.testing.expectEqual(@as(f64, 17), excluded[1]);
    try std.testing.expectEqual(@as(f64, 14), excluded[0]);
    try std.testing.expectEqual(@as(f64, 3), working[0]);
}

test "REDIST soil depth inactive disturbance preserves state" {
    var cumulative = [_]f64{ 1, 2 };
    const previous = [_]f64{ 1, 2 };
    var excluded = [_]f64{ 1, 2 };
    const volume = [_]f64{ 1, 2 };
    var working = [_]f64{ 3, 4 };
    var flags = [_]u8{ 0, 0 };
    const values = [_]f64{ 1, 1 };
    const changes = [_]f64{0} ** 8;
    try apply(-1, 1, 1, .{ .cumulative_depth_m = &cumulative, .previous_cumulative_depth_m = &previous, .freeze_excluded_cumulative_depth_m = &excluded, .previous_layer_volume_m3 = &volume, .working_layer_volume_m3 = &working, .pond_profile_flag_by_layer = &flags }, .{ .bulk_density_megagrams_per_m3 = &values, .water_m3 = &values, .ice_m3 = &values, .area_m2 = &values, .boundary_change_m_by_layer_channel = &changes, .pond_profile_flag = 0, .zero_tolerance = 0 });
    try std.testing.expectEqual(@as(f64, 2), cumulative[1]);
    try std.testing.expectEqual(@as(f64, 4), working[1]);
}

test "REDIST soil depth apply rejects dimensions and pond area" {
    var cumulative = [_]f64{ 0, 0 };
    const values = [_]f64{ 0, 0 };
    var mutable = [_]f64{ 0, 0 };
    var flags = [_]u8{ 0, 0 };
    const changes = [_]f64{0} ** 8;
    const state = State{ .cumulative_depth_m = &cumulative, .previous_cumulative_depth_m = &values, .freeze_excluded_cumulative_depth_m = &mutable, .previous_layer_volume_m3 = &values, .working_layer_volume_m3 = &mutable, .pond_profile_flag_by_layer = &flags };
    try std.testing.expectError(error.SoilProfileDepthApplyDimensionMismatch, apply(0, 1, 2, state, .{ .bulk_density_megagrams_per_m3 = &values, .water_m3 = &values, .ice_m3 = &values, .area_m2 = &values, .boundary_change_m_by_layer_channel = &changes, .pond_profile_flag = 0, .zero_tolerance = 0 }));
    const zero_area = [_]f64{ 0, 0 };
    try std.testing.expectError(error.InvalidPondGeometry, apply(0, 1, 1, state, .{ .bulk_density_megagrams_per_m3 = &values, .water_m3 = &values, .ice_m3 = &values, .area_m2 = &zero_area, .boundary_change_m_by_layer_channel = &changes, .pond_profile_flag = 0, .zero_tolerance = 0 }));
}
