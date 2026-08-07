const std = @import("std");

pub const Inputs = struct {
    bulk_density_megagrams_per_m3: []const f64, // BKDS
    area_m2: []const f64, // AREA(3,LX)
    initial_bulk_density_megagrams_per_m3: []const f64, // BKDSI
    top_layer_bulk_mass_megagrams: f64, // BKVLNU
    top_layer_volume_m3: f64, // VOLX(NU)
    net_sediment_flux_megagrams: f64, // TSEDER
    surface_sediment_flux_megagrams: f64, // TSEDSK
    negligible_sediment_megagrams: f64, // ZEROS
};

/// Direct translation of REDIST 7842--7864. Arrays retain Fortran layer
/// numbers and are filled in the original NL-to-NU reverse order.
pub fn build(
    disturbance_mode: i8,
    top_layer: usize,
    bottom_layer: usize,
    organic_mineral_boundary_layer: usize, // NUL
    inputs: Inputs,
    cumulative_change_m: []f64,
    boundary_change_m: []f64,
) !void {
    const len = inputs.bulk_density_megagrams_per_m3.len;
    if (len == 0 or top_layer > organic_mineral_boundary_layer or
        organic_mineral_boundary_layer > bottom_layer or bottom_layer >= len or
        inputs.area_m2.len != len or inputs.initial_bulk_density_megagrams_per_m3.len != len or
        cumulative_change_m.len != len or boundary_change_m.len != len)
        return error.SoilErosionDepthDimensionMismatch;
    inline for (.{ inputs.bulk_density_megagrams_per_m3, inputs.area_m2, inputs.initial_bulk_density_megagrams_per_m3 }) |values|
        for (values) |value| if (!std.math.isFinite(value)) return error.InvalidSoilErosionDepthInput;
    inline for (.{ inputs.top_layer_bulk_mass_megagrams, inputs.top_layer_volume_m3, inputs.net_sediment_flux_megagrams, inputs.surface_sediment_flux_megagrams, inputs.negligible_sediment_megagrams }) |value|
        if (!std.math.isFinite(value)) return error.InvalidSoilErosionDepthInput;
    if (inputs.negligible_sediment_megagrams < 0 or inputs.bulk_density_megagrams_per_m3[bottom_layer] <= 0)
        return error.InvalidSoilErosionDepthInput;

    const active = (disturbance_mode == 1 or disturbance_mode == 3) and
        (@abs(inputs.net_sediment_flux_megagrams) > inputs.negligible_sediment_megagrams or
            inputs.surface_sediment_flux_megagrams > inputs.negligible_sediment_megagrams);
    var erosion_change_m: f64 = 0.0;
    var layer = bottom_layer;
    while (true) : (layer -= 1) {
        if (inputs.bulk_density_megagrams_per_m3[layer] > 0) {
            if (active) {
                if (layer == bottom_layer) {
                    if (inputs.bulk_density_megagrams_per_m3[top_layer] > 0) {
                        if (inputs.area_m2[layer] <= 0 or inputs.top_layer_volume_m3 <= 0 or
                            inputs.top_layer_bulk_mass_megagrams <= 0)
                            return error.InvalidSoilErosionGeometry;
                        erosion_change_m = inputs.net_sediment_flux_megagrams /
                            (inputs.area_m2[layer] * inputs.top_layer_bulk_mass_megagrams / inputs.top_layer_volume_m3);
                    } else {
                        erosion_change_m = 0.0;
                    }
                }
                if (layer == organic_mineral_boundary_layer) {
                    if (inputs.area_m2[layer] <= 0 or inputs.initial_bulk_density_megagrams_per_m3[layer] <= 0)
                        return error.InvalidSoilErosionGeometry;
                    erosion_change_m = erosion_change_m + inputs.surface_sediment_flux_megagrams /
                        (inputs.area_m2[layer] * inputs.initial_bulk_density_megagrams_per_m3[layer]);
                }
                cumulative_change_m[layer] = erosion_change_m;
                boundary_change_m[layer] = erosion_change_m;
            } else {
                cumulative_change_m[layer] = 0.0;
                boundary_change_m[layer] = 0.0;
            }
            if (!std.math.isFinite(cumulative_change_m[layer]))
                return error.NonFiniteSoilErosionDepthPlan;
        }
        if (layer == top_layer) break;
    }
}

test "REDIST erosion depth retains bottom term then augments at NUL" {
    const density = [_]f64{ 0, 1, 1, 1 };
    const area = [_]f64{ 1, 1, 4, 2 };
    const initial_density = [_]f64{ 1, 1, 2, 1 };
    var cumulative = [_]f64{ 0, 0, 0, 0 };
    var boundary = [_]f64{ 0, 0, 0, 0 };
    try build(3, 1, 3, 2, .{
        .bulk_density_megagrams_per_m3 = &density,
        .area_m2 = &area,
        .initial_bulk_density_megagrams_per_m3 = &initial_density,
        .top_layer_bulk_mass_megagrams = 10,
        .top_layer_volume_m3 = 5,
        .net_sediment_flux_megagrams = 8,
        .surface_sediment_flux_megagrams = 8,
        .negligible_sediment_megagrams = 0,
    }, &cumulative, &boundary);
    try std.testing.expectEqual(@as(f64, 2), cumulative[3]);
    try std.testing.expectEqual(@as(f64, 3), cumulative[2]);
    try std.testing.expectEqual(@as(f64, 3), cumulative[1]);
}

test "REDIST erosion depth inactive gate zeros every soil layer" {
    const values = [_]f64{ 1, 1, 1, 1 };
    var cumulative = [_]f64{ 9, 9, 9, 9 };
    var boundary = [_]f64{ 9, 9, 9, 9 };
    try build(2, 1, 3, 2, .{
        .bulk_density_megagrams_per_m3 = &values,
        .area_m2 = &values,
        .initial_bulk_density_megagrams_per_m3 = &values,
        .top_layer_bulk_mass_megagrams = 1,
        .top_layer_volume_m3 = 1,
        .net_sediment_flux_megagrams = 1,
        .surface_sediment_flux_megagrams = 1,
        .negligible_sediment_megagrams = 0,
    }, &cumulative, &boundary);
    try std.testing.expectEqual(@as(f64, 0), cumulative[1]);
    try std.testing.expectEqual(@as(f64, 0), cumulative[2]);
    try std.testing.expectEqual(@as(f64, 0), cumulative[3]);
}

test "REDIST erosion depth ignores TSEDER when top layer is pond" {
    const density = [_]f64{ 0, 0, 1, 1 };
    const area = [_]f64{ 1, 1, 2, 2 };
    const initial_density = [_]f64{ 1, 1, 2, 1 };
    var cumulative = [_]f64{ 0, 0, 0, 0 };
    var boundary = [_]f64{ 0, 0, 0, 0 };
    try build(1, 1, 3, 2, .{ .bulk_density_megagrams_per_m3 = &density, .area_m2 = &area, .initial_bulk_density_megagrams_per_m3 = &initial_density, .top_layer_bulk_mass_megagrams = 10, .top_layer_volume_m3 = 5, .net_sediment_flux_megagrams = 8, .surface_sediment_flux_megagrams = 4, .negligible_sediment_megagrams = 0 }, &cumulative, &boundary);
    try std.testing.expectEqual(@as(f64, 0), cumulative[3]);
    try std.testing.expectEqual(@as(f64, 1), cumulative[2]);
}

test "REDIST erosion depth rejects dimensions and invalid geometry" {
    const values = [_]f64{ 1, 1, 1, 1 };
    var cumulative = [_]f64{ 0, 0, 0, 0 };
    var boundary = [_]f64{ 0, 0, 0, 0 };
    const inputs = Inputs{ .bulk_density_megagrams_per_m3 = &values, .area_m2 = &values, .initial_bulk_density_megagrams_per_m3 = &values, .top_layer_bulk_mass_megagrams = 1, .top_layer_volume_m3 = 1, .net_sediment_flux_megagrams = 1, .surface_sediment_flux_megagrams = 1, .negligible_sediment_megagrams = 0 };
    try std.testing.expectError(error.SoilErosionDepthDimensionMismatch, build(1, 3, 1, 2, inputs, &cumulative, &boundary));
    var bad_area = values;
    bad_area[3] = 0;
    var bad = inputs;
    bad.area_m2 = &bad_area;
    try std.testing.expectError(error.InvalidSoilErosionGeometry, build(1, 1, 3, 2, bad, &cumulative, &boundary));
}
