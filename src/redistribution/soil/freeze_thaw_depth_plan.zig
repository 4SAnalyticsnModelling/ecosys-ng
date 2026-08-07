const std = @import("std");

pub const Inputs = struct {
    ice_volume_change_m3: []const f64, // DVOLI
    micropore_fraction: []const f64, // FMPR
    area_m2: []const f64,
    bulk_density_megagrams_per_m3: []const f64, // BKDS
    ice_density_volume_correction: f64, // DENSJ
    negligible_ice_change_m3: f64, // ZEROS
    zero_tolerance: f64, // ZERO
};

/// Direct translation of REDIST 7960--7996 in NL-to-NU order. Returns IFLGM.
pub fn build(
    top_layer: usize,
    bottom_layer: usize,
    inputs: Inputs,
    cumulative_change_m: []f64,
    boundary_change_m: []f64,
) !bool {
    const len = inputs.ice_volume_change_m3.len;
    if (len == 0 or top_layer == 0 or top_layer > bottom_layer or bottom_layer >= len or
        inputs.micropore_fraction.len != len or inputs.area_m2.len != len or
        inputs.bulk_density_megagrams_per_m3.len != len or cumulative_change_m.len != len or
        boundary_change_m.len != len)
        return error.SoilFreezeThawDepthDimensionMismatch;
    inline for (.{ inputs.ice_volume_change_m3, inputs.micropore_fraction, inputs.area_m2, inputs.bulk_density_megagrams_per_m3 }) |values|
        for (values) |value| if (!std.math.isFinite(value)) return error.InvalidSoilFreezeThawDepthInput;
    inline for (.{ inputs.ice_density_volume_correction, inputs.negligible_ice_change_m3, inputs.zero_tolerance }) |value|
        if (!std.math.isFinite(value)) return error.InvalidSoilFreezeThawDepthInput;
    if (inputs.negligible_ice_change_m3 < 0 or inputs.zero_tolerance < 0)
        return error.InvalidSoilFreezeThawDepthInput;

    var changed = false;
    var layer = bottom_layer;
    while (true) : (layer -= 1) {
        if (inputs.bulk_density_megagrams_per_m3[layer] > inputs.zero_tolerance) {
            if (@abs(inputs.ice_volume_change_m3[layer]) > inputs.negligible_ice_change_m3) {
                const fraction = inputs.micropore_fraction[layer];
                const area = inputs.area_m2[layer];
                if (fraction <= 0 or area <= 0) return error.InvalidSoilFreezeThawGeometry;
                const layer_change_m = inputs.ice_volume_change_m3[layer] *
                    inputs.ice_density_volume_correction / (fraction * area);
                changed = true;
                if (layer == bottom_layer) {
                    cumulative_change_m[layer] = layer_change_m;
                    boundary_change_m[layer] = 0.0;
                } else {
                    cumulative_change_m[layer] = layer_change_m + cumulative_change_m[layer + 1];
                    boundary_change_m[layer] = cumulative_change_m[layer + 1];
                    if (layer == top_layer or inputs.bulk_density_megagrams_per_m3[layer - 1] <= inputs.zero_tolerance) {
                        cumulative_change_m[layer - 1] = cumulative_change_m[layer];
                        boundary_change_m[layer - 1] = cumulative_change_m[layer];
                    }
                }
            } else if (layer == bottom_layer) {
                cumulative_change_m[layer] = 0.0;
                boundary_change_m[layer] = 0.0;
            } else {
                cumulative_change_m[layer] = cumulative_change_m[layer + 1];
                boundary_change_m[layer] = cumulative_change_m[layer + 1];
                if (layer == top_layer) {
                    cumulative_change_m[layer - 1] = cumulative_change_m[layer];
                    boundary_change_m[layer - 1] = cumulative_change_m[layer];
                }
            }
            if (!std.math.isFinite(cumulative_change_m[layer]) or !std.math.isFinite(boundary_change_m[layer]))
                return error.NonFiniteSoilFreezeThawDepthPlan;
        }
        if (layer == top_layer) break;
    }
    return changed;
}

test "REDIST freeze thaw accumulates active reverse layers and surface boundary" {
    const ice = [_]f64{ 0, 1, 2, 3, 0 };
    const fraction = [_]f64{ 1, 1, 1, 1, 1 };
    const area = [_]f64{ 1, 1, 1, 1, 1 };
    const density = [_]f64{ 0, 1, 1, 1, 0 };
    var cumulative = [_]f64{ 0, 0, 0, 0, 0 };
    var boundary = [_]f64{ 0, 0, 0, 0, 0 };
    const changed = try build(1, 3, .{ .ice_volume_change_m3 = &ice, .micropore_fraction = &fraction, .area_m2 = &area, .bulk_density_megagrams_per_m3 = &density, .ice_density_volume_correction = 0.1, .negligible_ice_change_m3 = 0, .zero_tolerance = 1e-12 }, &cumulative, &boundary);
    try std.testing.expect(changed);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), cumulative[3], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), cumulative[2], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), cumulative[1], 1e-15);
    try std.testing.expectEqual(cumulative[1], cumulative[0]);
}

test "REDIST freeze thaw inactive layer inherits active lower change" {
    const ice = [_]f64{ 0, 0, 0, 2, 0 };
    const values = [_]f64{ 1, 1, 1, 1, 1 };
    var cumulative = [_]f64{ 0, 0, 0, 0, 0 };
    var boundary = [_]f64{ 0, 0, 0, 0, 0 };
    _ = try build(1, 3, .{ .ice_volume_change_m3 = &ice, .micropore_fraction = &values, .area_m2 = &values, .bulk_density_megagrams_per_m3 = &values, .ice_density_volume_correction = 0.1, .negligible_ice_change_m3 = 0, .zero_tolerance = 0 }, &cumulative, &boundary);
    try std.testing.expectEqual(cumulative[3], cumulative[2]);
    try std.testing.expectEqual(cumulative[2], cumulative[1]);
    try std.testing.expectEqual(cumulative[1], cumulative[0]);
}

test "REDIST freeze thaw reports no change and zeros inactive profile" {
    const ice = [_]f64{ 0, 0, 0, 0 };
    const values = [_]f64{ 1, 1, 1, 1 };
    var cumulative = [_]f64{ 9, 9, 9, 9 };
    var boundary = [_]f64{ 9, 9, 9, 9 };
    const changed = try build(1, 3, .{ .ice_volume_change_m3 = &ice, .micropore_fraction = &values, .area_m2 = &values, .bulk_density_megagrams_per_m3 = &values, .ice_density_volume_correction = 0.1, .negligible_ice_change_m3 = 0, .zero_tolerance = 0 }, &cumulative, &boundary);
    try std.testing.expect(!changed);
    try std.testing.expectEqual(@as(f64, 0), cumulative[3]);
    try std.testing.expectEqual(@as(f64, 0), cumulative[1]);
}

test "REDIST freeze thaw rejects dimensions and invalid geometry" {
    const values = [_]f64{ 1, 1, 1, 1 };
    var cumulative = [_]f64{ 0, 0, 0, 0 };
    var boundary = [_]f64{ 0, 0, 0, 0 };
    const inputs = Inputs{ .ice_volume_change_m3 = &values, .micropore_fraction = &values, .area_m2 = &values, .bulk_density_megagrams_per_m3 = &values, .ice_density_volume_correction = 0.1, .negligible_ice_change_m3 = 0, .zero_tolerance = 0 };
    try std.testing.expectError(error.SoilFreezeThawDepthDimensionMismatch, build(3, 1, inputs, &cumulative, &boundary));
    var fraction = values;
    fraction[3] = 0;
    var bad = inputs;
    bad.micropore_fraction = &fraction;
    try std.testing.expectError(error.InvalidSoilFreezeThawGeometry, build(1, 3, bad, &cumulative, &boundary));
}
