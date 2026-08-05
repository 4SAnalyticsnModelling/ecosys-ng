const std = @import("std");

pub const DepthAdjustmentMode = enum { inactive, limited_by_lower_layer, direct };

pub const TransferPlan = struct {
    depth_change_m: f64,
    unconstrained_depth_change_m: f64,
    source_layer: usize,
    destination_layer: usize,
    fraction: f64,
    adjustment_mode: DepthAdjustmentMode,
};

pub const Inputs = struct {
    heat_capacity_megajoules_per_k: []const f64, // VHCPW
    minimum_heat_capacity_megajoules_per_k: f64, // VHCPWX
    target_volume_m3: []const f64, // VOLSI
    current_volume_m3: []const f64, // VOLSL
    solid_water_equivalent_m3: []const f64, // VOLSSL
    liquid_water_m3: []const f64, // VOLWSL
    ice_m3: []const f64, // VOLISL
    snow_density_megagrams_per_m3: []const f64, // DENSS
    area_by_layer_m2: []const f64, // AREA(3,L)
    surface_area_m2: f64, // AREA(3,NU)
    negligible_volume_m3: f64, // ZEROS2
    depth_tolerance_m: f64, // ZERO
};

fn validateSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 7419--7481. Index zero is retained so Zig
/// indexes equal the original snow-layer numbers; executable layers are 1..JS-1.
pub fn buildPlans(
    inputs: Inputs,
    cumulative_depth_m: []f64,
    layer_depth_m: []f64,
    plans_by_layer: []TransferPlan,
) !void {
    const len = inputs.heat_capacity_megajoules_per_k.len;
    if (len < 3 or inputs.target_volume_m3.len != len or inputs.current_volume_m3.len != len or
        inputs.solid_water_equivalent_m3.len != len or inputs.liquid_water_m3.len != len or
        inputs.ice_m3.len != len or inputs.snow_density_megagrams_per_m3.len != len or
        inputs.area_by_layer_m2.len != len or cumulative_depth_m.len != len or
        layer_depth_m.len != len or plans_by_layer.len != len)
        return error.SnowRelayeringDimensionMismatch;
    if (!validateSlice(inputs.heat_capacity_megajoules_per_k) or !validateSlice(inputs.target_volume_m3) or
        !validateSlice(inputs.current_volume_m3) or !validateSlice(inputs.solid_water_equivalent_m3) or
        !validateSlice(inputs.liquid_water_m3) or !validateSlice(inputs.ice_m3) or
        !validateSlice(inputs.snow_density_megagrams_per_m3) or !validateSlice(inputs.area_by_layer_m2) or
        !validateSlice(cumulative_depth_m) or !validateSlice(layer_depth_m) or
        !std.math.isFinite(inputs.minimum_heat_capacity_megajoules_per_k) or
        !std.math.isFinite(inputs.surface_area_m2) or inputs.surface_area_m2 <= 0 or
        !std.math.isFinite(inputs.negligible_volume_m3) or inputs.negligible_volume_m3 < 0 or
        !std.math.isFinite(inputs.depth_tolerance_m) or inputs.depth_tolerance_m < 0)
        return error.InvalidSnowRelayeringInput;
    if (inputs.heat_capacity_megajoules_per_k[1] <= inputs.minimum_heat_capacity_megajoules_per_k) return;

    const js = len - 1;
    for (1..js) |layer| {
        var depth_change_m: f64 = 0.0;
        var unconstrained_m: f64 = 0.0;
        var mode: DepthAdjustmentMode = .inactive;
        if (inputs.current_volume_m3[layer] > inputs.negligible_volume_m3) {
            const density = inputs.snow_density_megagrams_per_m3[layer];
            if (density <= 0) return error.InvalidSnowDensity;
            unconstrained_m = (inputs.target_volume_m3[layer] -
                inputs.solid_water_equivalent_m3[layer] / density -
                inputs.liquid_water_m3[layer] - inputs.ice_m3[layer]) /
                inputs.surface_area_m2;
            if (unconstrained_m < -inputs.depth_tolerance_m or layer_depth_m[layer + 1] > inputs.depth_tolerance_m) {
                depth_change_m = @min(unconstrained_m, layer_depth_m[layer + 1]);
                mode = .limited_by_lower_layer;
            } else {
                // REDIST recomputes DDLYXS before assigning it in this branch.
                unconstrained_m = (inputs.current_volume_m3[layer] -
                    inputs.solid_water_equivalent_m3[layer] / density -
                    inputs.liquid_water_m3[layer] - inputs.ice_m3[layer]) /
                    inputs.surface_area_m2;
                depth_change_m = unconstrained_m;
                mode = .direct;
            }
        }
        cumulative_depth_m[layer] = cumulative_depth_m[layer] + depth_change_m;
        layer_depth_m[layer] = cumulative_depth_m[layer] - cumulative_depth_m[layer - 1];

        var source = layer;
        var destination = layer;
        var fraction: f64 = 0.0;
        if (@abs(depth_change_m) > inputs.depth_tolerance_m) {
            if (depth_change_m > 0.0) {
                destination = layer;
                source = layer + 1;
                if (depth_change_m < unconstrained_m) {
                    fraction = 1.0;
                } else {
                    if (inputs.current_volume_m3[source] <= 0) return error.InvalidSnowSourceVolume;
                    fraction = @min(1.0, depth_change_m * inputs.area_by_layer_m2[source] / inputs.current_volume_m3[source]);
                }
            } else {
                destination = layer + 1;
                source = layer;
                if (inputs.current_volume_m3[source] < inputs.target_volume_m3[source]) {
                    fraction = 0.0;
                } else {
                    if (inputs.current_volume_m3[source] <= 0) return error.InvalidSnowSourceVolume;
                    fraction = @min(1.0, -depth_change_m * inputs.area_by_layer_m2[source] / inputs.current_volume_m3[source]);
                }
            }
        }
        if (!std.math.isFinite(depth_change_m) or !std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.NonFiniteSnowRelayeringPlan;
        plans_by_layer[layer] = .{
            .depth_change_m = depth_change_m,
            .unconstrained_depth_change_m = unconstrained_m,
            .source_layer = source,
            .destination_layer = destination,
            .fraction = fraction,
            .adjustment_mode = mode,
        };
    }
}

fn baseInputs(values: []const f64) Inputs {
    return .{
        .heat_capacity_megajoules_per_k = values,
        .minimum_heat_capacity_megajoules_per_k = 0,
        .target_volume_m3 = values,
        .current_volume_m3 = values,
        .solid_water_equivalent_m3 = values,
        .liquid_water_m3 = &.{ 0, 0, 0 },
        .ice_m3 = &.{ 0, 0, 0 },
        .snow_density_megagrams_per_m3 = &.{ 1, 1, 1 },
        .area_by_layer_m2 = &.{ 1, 1, 1 },
        .surface_area_m2 = 1,
        .negligible_volume_m3 = 0,
        .depth_tolerance_m = 1e-12,
    };
}

test "REDIST snow relayering positive deficit consumes lower layer" {
    const heat = [_]f64{ 0, 1, 1 };
    var inputs = baseInputs(&heat);
    inputs.target_volume_m3 = &.{ 0, 10, 0 };
    inputs.current_volume_m3 = &.{ 0, 2, 4 };
    inputs.solid_water_equivalent_m3 = &.{ 0, 2, 4 };
    var cumulative = [_]f64{ 0, 2, 5 };
    var depth = [_]f64{ 0, 2, 3 };
    var plans: [3]TransferPlan = undefined;
    try buildPlans(inputs, &cumulative, &depth, &plans);
    try std.testing.expectEqual(@as(f64, 3), plans[1].depth_change_m);
    try std.testing.expectEqual(@as(usize, 2), plans[1].source_layer);
    try std.testing.expectEqual(@as(usize, 1), plans[1].destination_layer);
    try std.testing.expectEqual(@as(f64, 1), plans[1].fraction);
}

test "REDIST snow relayering negative excess computes fractional transfer" {
    const heat = [_]f64{ 0, 1, 1 };
    var inputs = baseInputs(&heat);
    inputs.target_volume_m3 = &.{ 0, 1, 0 };
    inputs.current_volume_m3 = &.{ 0, 3, 0 };
    inputs.solid_water_equivalent_m3 = &.{ 0, 3, 0 };
    var cumulative = [_]f64{ 0, 3, 3 };
    var depth = [_]f64{ 0, 3, 0 };
    var plans: [3]TransferPlan = undefined;
    try buildPlans(inputs, &cumulative, &depth, &plans);
    try std.testing.expectEqual(@as(f64, -2), plans[1].depth_change_m);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0 / 3.0), plans[1].fraction, 1e-15);
    try std.testing.expectEqual(@as(usize, 1), plans[1].source_layer);
    try std.testing.expectEqual(@as(usize, 2), plans[1].destination_layer);
}

test "REDIST snow relayering outer heat gate preserves geometry" {
    const heat = [_]f64{ 0, 1, 1 };
    var inputs = baseInputs(&heat);
    inputs.minimum_heat_capacity_megajoules_per_k = 1;
    var cumulative = [_]f64{ 0, 2, 5 };
    var depth = [_]f64{ 0, 2, 3 };
    var plans: [3]TransferPlan = undefined;
    try buildPlans(inputs, &cumulative, &depth, &plans);
    try std.testing.expectEqual(@as(f64, 2), cumulative[1]);
    try std.testing.expectEqual(@as(f64, 2), depth[1]);
}

test "REDIST snow relayering rejects runtime dimensions and invalid density" {
    const short = [_]f64{ 0, 1 };
    var inputs = baseInputs(&short);
    var cumulative = [_]f64{ 0, 1 };
    var depth = [_]f64{ 0, 1 };
    var plans: [2]TransferPlan = undefined;
    try std.testing.expectError(error.SnowRelayeringDimensionMismatch, buildPlans(inputs, &cumulative, &depth, &plans));

    const heat = [_]f64{ 0, 1, 1 };
    inputs = baseInputs(&heat);
    inputs.snow_density_megagrams_per_m3 = &.{ 1, 0, 1 };
    var cumulative_ok = [_]f64{ 0, 1, 2 };
    var depth_ok = [_]f64{ 0, 1, 1 };
    var plans_ok: [3]TransferPlan = undefined;
    try std.testing.expectError(error.InvalidSnowDensity, buildPlans(inputs, &cumulative_ok, &depth_ok, &plans_ok));
}
