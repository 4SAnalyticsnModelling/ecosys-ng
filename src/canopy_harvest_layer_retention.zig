const std = @import("std");

pub const State = struct {
    plant_retained_fraction_by_layer: []f64,
    harvest_retained_fraction_by_layer: []f64,
};

pub const Inputs = struct {
    canopy_layer_count: usize,
    harvest_type: i32,
    harvest_height_m: f64,
    thinning_fraction: f64,
    leaf_removal_efficiency: f64,
    canopy_layer_bottom_height_m: []const f64,
};

fn validateInputs(inputs: Inputs) !void {
    if (inputs.canopy_layer_count == 0 or inputs.canopy_layer_bottom_height_m.len != inputs.canopy_layer_count + 1) return error.HarvestLayerRetentionDimensionMismatch;
    const scalars = [_]f64{ inputs.harvest_height_m, inputs.thinning_fraction, inputs.leaf_removal_efficiency };
    for (scalars) |value| if (!std.math.isFinite(value)) return error.InvalidHarvestLayerRetentionInput;
    if (inputs.thinning_fraction < 0.0 or inputs.leaf_removal_efficiency < 0.0) return error.InvalidHarvestLayerRetentionInput;
    for (inputs.canopy_layer_bottom_height_m) |value| if (!std.math.isFinite(value)) return error.InvalidHarvestLayerRetentionInput;
}

fn validateState(state: State, count: usize) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        if (values.len != count) return error.HarvestLayerRetentionDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0 or value > 1.0) return error.InvalidHarvestLayerRetentionState;
    }
}

/// Exact GROSUB 8821--8847 top-to-bottom canopy-layer retained fractions.
/// Heights are m and fractions dimensionless. Output index is the logical
/// bottom-to-top canopy-layer index even though traversal is reversed.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateInputs(inputs);
    try validateState(state, inputs.canopy_layer_count);
    try validateState(workspace, inputs.canopy_layer_count);
    @memcpy(workspace.plant_retained_fraction_by_layer, state.plant_retained_fraction_by_layer);
    @memcpy(workspace.harvest_retained_fraction_by_layer, state.harvest_retained_fraction_by_layer);
    var reverse_layer = inputs.canopy_layer_count;
    while (reverse_layer > 0) {
        reverse_layer -= 1;
        const layer = reverse_layer;
        var plant_retained_fraction: f64 = undefined;
        var harvest_retained_fraction: f64 = undefined;
        if (inputs.harvest_type != 4 and inputs.harvest_type != 6) {
            const height_retained_fraction = if (inputs.harvest_type != 3)
                if (inputs.canopy_layer_bottom_height_m[layer + 1] > inputs.canopy_layer_bottom_height_m[layer])
                    @max(0.0, @min(1.0, 1.0 - (inputs.canopy_layer_bottom_height_m[layer + 1] - inputs.harvest_height_m) / (inputs.canopy_layer_bottom_height_m[layer + 1] - inputs.canopy_layer_bottom_height_m[layer])))
                else
                    1.0
            else
                0.0;
            if (inputs.thinning_fraction == 0.0) {
                plant_retained_fraction = @max(0.0, 1.0 - (1.0 - height_retained_fraction) * inputs.leaf_removal_efficiency);
                harvest_retained_fraction = plant_retained_fraction;
            } else {
                plant_retained_fraction = @max(0.0, 1.0 - inputs.thinning_fraction);
                harvest_retained_fraction = if (inputs.harvest_type == 0)
                    1.0 - (1.0 - height_retained_fraction) * inputs.leaf_removal_efficiency * inputs.thinning_fraction
                else
                    plant_retained_fraction;
            }
        } else {
            plant_retained_fraction = 0.0;
            harvest_retained_fraction = 0.0;
        }
        workspace.plant_retained_fraction_by_layer[layer] = plant_retained_fraction;
        workspace.harvest_retained_fraction_by_layer[layer] = harvest_retained_fraction;
    }
    try validateState(workspace, inputs.canopy_layer_count);
    @memcpy(state.plant_retained_fraction_by_layer, workspace.plant_retained_fraction_by_layer);
    @memcpy(state.harvest_retained_fraction_by_layer, workspace.harvest_retained_fraction_by_layer);
}

test "GROSUB canopy harvest layer fractions traverse top down and preserve geometry" {
    var plant = [_]f64{ 1, 1 };
    var harvest = [_]f64{ 1, 1 };
    var wp = [_]f64{0} ** 2;
    var wh = [_]f64{0} ** 2;
    const inputs: Inputs = .{ .canopy_layer_count = 2, .harvest_type = 2, .harvest_height_m = 1.5, .thinning_fraction = 0, .leaf_removal_efficiency = 1, .canopy_layer_bottom_height_m = &.{ 0, 1, 2 } };
    try apply(.{ .plant_retained_fraction_by_layer = &plant, .harvest_retained_fraction_by_layer = &harvest }, .{ .plant_retained_fraction_by_layer = &wp, .harvest_retained_fraction_by_layer = &wh }, inputs);
    try std.testing.expectEqualSlices(f64, &.{ 1, 0.5 }, &plant);
    try std.testing.expectEqualSlices(f64, &.{ 1, 0.5 }, &harvest);
}
