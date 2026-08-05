const std = @import("std");

pub const LayerState = struct {
    water_volume_m3: []const f64, // VOLW
    ice_volume_m3: []const f64, // VOLI
    cumulative_depth_m: []f64, // CDPTH
    layer_thickness_m: []f64, // DLYR(3,...)
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 9404--9409.
pub fn closeEmptySourceLayer(
    redistribution_pass: usize,
    source_layer: usize,
    destination_layer: usize,
    empty_volume_threshold_m3: f64,
    state: LayerState,
) !void {
    const layer_count = state.water_volume_m3.len;
    if (layer_count == 0 or source_layer >= layer_count or destination_layer >= layer_count or
        source_layer == destination_layer or state.ice_volume_m3.len != layer_count or
        state.cumulative_depth_m.len != layer_count or state.layer_thickness_m.len != layer_count)
        return error.PondEmptyLayerClosureDimensionMismatch;
    if (!finiteSlice(state.water_volume_m3) or !finiteSlice(state.ice_volume_m3) or
        !finiteSlice(state.cumulative_depth_m) or !finiteSlice(state.layer_thickness_m) or
        !std.math.isFinite(empty_volume_threshold_m3) or empty_volume_threshold_m3 < 0)
        return error.InvalidPondEmptyLayerClosureInput;
    const combined_source_volume_m3 = state.water_volume_m3[source_layer] + state.ice_volume_m3[source_layer];
    if (!std.math.isFinite(combined_source_volume_m3)) return error.NonFinitePondEmptyLayerClosureResult;
    if (redistribution_pass != 1 or combined_source_volume_m3 > empty_volume_threshold_m3) return;

    state.cumulative_depth_m[destination_layer] = state.cumulative_depth_m[source_layer];
    state.layer_thickness_m[source_layer] = 0;
}

test "REDIST first pass closes source at threshold in assignment order" {
    const water = [_]f64{ 0.3, 4 };
    const ice = [_]f64{ 0.2, 5 };
    var depth = [_]f64{ 1.25, 2.5 };
    var thickness = [_]f64{ 0.4, 0.8 };
    try closeEmptySourceLayer(1, 0, 1, 0.5, .{
        .water_volume_m3 = &water,
        .ice_volume_m3 = &ice,
        .cumulative_depth_m = &depth,
        .layer_thickness_m = &thickness,
    });
    try std.testing.expectEqual(@as(f64, 1.25), depth[1]);
    try std.testing.expectEqual(@as(f64, 0), thickness[0]);
    try std.testing.expectEqual(@as(f64, 0.8), thickness[1]);
}

test "REDIST empty-layer closure obeys exact pass and threshold gates" {
    const water = [_]f64{ 0.3, 4 };
    const ice = [_]f64{ 0.2, 5 };
    var depth = [_]f64{ 1.25, 2.5 };
    var thickness = [_]f64{ 0.4, 0.8 };
    const state = LayerState{ .water_volume_m3 = &water, .ice_volume_m3 = &ice, .cumulative_depth_m = &depth, .layer_thickness_m = &thickness };
    try closeEmptySourceLayer(2, 0, 1, 0.5, state);
    try closeEmptySourceLayer(1, 0, 1, 0.499999999, state);
    try std.testing.expectEqual(@as(f64, 2.5), depth[1]);
    try std.testing.expectEqual(@as(f64, 0.4), thickness[0]);
}

test "REDIST empty-layer closure validation is atomic" {
    const water = [_]f64{ 0, 4 };
    const ice = [_]f64{ 0, 5 };
    var depth = [_]f64{ 1.25, 2.5 };
    var thickness = [_]f64{ 0.4, std.math.nan(f64) };
    try std.testing.expectError(error.InvalidPondEmptyLayerClosureInput, closeEmptySourceLayer(1, 0, 1, 0, .{
        .water_volume_m3 = &water,
        .ice_volume_m3 = &ice,
        .cumulative_depth_m = &depth,
        .layer_thickness_m = &thickness,
    }));
    try std.testing.expectEqual(@as(f64, 2.5), depth[1]);
    try std.testing.expectEqual(@as(f64, 0.4), thickness[0]);
}
