const std = @import("std");

pub const Inputs = struct {
    organic_carbon_change_g: []const f64, // DORGC
    area_m2: []const f64,
    macropore_fraction: []const f64, // FHOL
    initial_bulk_density_megagrams_per_m3: []const f64, // BKDSI
    bulk_density_megagrams_per_m3: []const f64, // BKDS
    organic_carbon_concentration: []const f64, // CORGC
    initial_depth_m: []const f64, // DLYRI
    current_depth_m: []const f64, // DLYR
    organic_boundary_concentration: f64, // FORGC
    negligible_carbon_change_g: f64, // ZEROS
    zero_tolerance: f64, // ZERO
};

/// Direct translation of REDIST 7896--7948 in NL-to-NU order. Array indexes
/// equal Fortran layer numbers and include NU-1 and NL+1 boundary slots.
pub fn build(
    disturbance_mode: i8,
    top_layer: usize,
    bottom_layer: usize,
    inputs: Inputs,
    cumulative_change_m: []f64,
    boundary_change_m: []f64,
    organic_boundary_flag: []u8,
    water_boundary_flag: []u8,
) !void {
    const len = inputs.organic_carbon_change_g.len;
    if (len == 0 or top_layer == 0 or top_layer > bottom_layer or bottom_layer + 1 >= len or
        inputs.area_m2.len != len or inputs.macropore_fraction.len != len or
        inputs.initial_bulk_density_megagrams_per_m3.len != len or inputs.bulk_density_megagrams_per_m3.len != len or
        inputs.organic_carbon_concentration.len != len or inputs.initial_depth_m.len != len or
        inputs.current_depth_m.len != len or cumulative_change_m.len != len or boundary_change_m.len != len or
        organic_boundary_flag.len != len or water_boundary_flag.len != len)
        return error.SoilSocDepthDimensionMismatch;
    inline for (.{ inputs.organic_carbon_change_g, inputs.area_m2, inputs.macropore_fraction, inputs.initial_bulk_density_megagrams_per_m3, inputs.bulk_density_megagrams_per_m3, inputs.organic_carbon_concentration, inputs.initial_depth_m, inputs.current_depth_m }) |values|
        for (values) |value| if (!std.math.isFinite(value)) return error.InvalidSoilSocDepthInput;
    inline for (.{ inputs.organic_boundary_concentration, inputs.negligible_carbon_change_g, inputs.zero_tolerance }) |value|
        if (!std.math.isFinite(value)) return error.InvalidSoilSocDepthInput;
    if (inputs.negligible_carbon_change_g < 0 or inputs.zero_tolerance < 0)
        return error.InvalidSoilSocDepthInput;

    var layer = bottom_layer;
    while (true) : (layer -= 1) {
        if (inputs.bulk_density_megagrams_per_m3[layer] > inputs.zero_tolerance) {
            const active = (disturbance_mode == 2 or disturbance_mode == 3) and
                @abs(inputs.organic_carbon_change_g[layer]) > inputs.negligible_carbon_change_g;
            if (active) {
                const area = inputs.area_m2[layer];
                const solid_fraction = 1.0 - inputs.macropore_fraction[layer];
                const density = inputs.initial_bulk_density_megagrams_per_m3[layer];
                if (area <= 0 or solid_fraction <= 0 or density <= 0)
                    return error.InvalidSoilSocGeometry;
                const layer_change_m = 1.82e-6 * inputs.organic_carbon_change_g[layer] / area /
                    (solid_fraction * density);
                const segment_bottom = layer == bottom_layer or
                    inputs.bulk_density_megagrams_per_m3[layer + 1] <= inputs.zero_tolerance or
                    (inputs.organic_carbon_concentration[layer] >= inputs.organic_boundary_concentration and
                        inputs.organic_carbon_concentration[layer + 1] < inputs.organic_boundary_concentration);
                if (segment_bottom) {
                    cumulative_change_m[layer] = layer_change_m;
                    boundary_change_m[layer] = 0.0;
                    organic_boundary_flag[layer] = 1;
                    if (layer > top_layer) water_boundary_flag[layer] = 2;
                } else {
                    cumulative_change_m[layer] = layer_change_m + cumulative_change_m[layer + 1];
                    boundary_change_m[layer] = cumulative_change_m[layer + 1] +
                        inputs.initial_depth_m[layer] - inputs.current_depth_m[layer];
                    organic_boundary_flag[layer] = 1;
                    if (layer == top_layer or inputs.bulk_density_megagrams_per_m3[layer - 1] <= inputs.zero_tolerance) {
                        cumulative_change_m[layer - 1] = cumulative_change_m[layer];
                        boundary_change_m[layer - 1] = cumulative_change_m[layer];
                        organic_boundary_flag[layer - 1] = 1;
                    }
                }
            } else if (layer == bottom_layer) {
                cumulative_change_m[layer] = 0.0;
                boundary_change_m[layer] = 0.0;
                organic_boundary_flag[layer] = 0;
            } else {
                cumulative_change_m[layer] = cumulative_change_m[layer + 1];
                boundary_change_m[layer] = cumulative_change_m[layer + 1];
                organic_boundary_flag[layer] = 0;
            }
            if (!std.math.isFinite(cumulative_change_m[layer]) or !std.math.isFinite(boundary_change_m[layer]))
                return error.NonFiniteSoilSocDepthPlan;
        }
        if (layer == top_layer) break;
    }
}

const Fixture = struct {
    dorgc: [5]f64 = .{ 0, 1, 2, 3, 0 },
    area: [5]f64 = .{ 1, 1, 1, 1, 1 },
    fhol: [5]f64 = .{ 0, 0, 0, 0, 0 },
    density: [5]f64 = .{ 0, 1, 1, 1, 0 },
    corgc: [5]f64 = .{ 0, 2, 2, 2, 0 },
    initial_depth: [5]f64 = .{ 0, 2, 2, 2, 0 },
    current_depth: [5]f64 = .{ 0, 1, 1, 1, 0 },
    cumulative: [5]f64 = .{ 0, 0, 0, 0, 0 },
    boundary: [5]f64 = .{ 0, 0, 0, 0, 0 },
    organic_flag: [5]u8 = .{ 0, 0, 0, 0, 0 },
    water_flag: [5]u8 = .{ 0, 0, 0, 0, 0 },
    fn inputs(self: *Fixture) Inputs {
        return .{ .organic_carbon_change_g = &self.dorgc, .area_m2 = &self.area, .macropore_fraction = &self.fhol, .initial_bulk_density_megagrams_per_m3 = &self.density, .bulk_density_megagrams_per_m3 = &self.density, .organic_carbon_concentration = &self.corgc, .initial_depth_m = &self.initial_depth, .current_depth_m = &self.current_depth, .organic_boundary_concentration = 1, .negligible_carbon_change_g = 0, .zero_tolerance = 1e-12 };
    }
};

test "REDIST SOC depth accumulates active layers and propagates NU boundary" {
    var fixture = Fixture{};
    try build(2, 1, 3, fixture.inputs(), &fixture.cumulative, &fixture.boundary, &fixture.organic_flag, &fixture.water_flag);
    try std.testing.expectApproxEqAbs(@as(f64, 3 * 1.82e-6), fixture.cumulative[3], 1e-18);
    try std.testing.expectApproxEqAbs(@as(f64, 5 * 1.82e-6), fixture.cumulative[2], 1e-18);
    try std.testing.expectApproxEqAbs(@as(f64, 6 * 1.82e-6), fixture.cumulative[1], 1e-18);
    try std.testing.expectEqual(fixture.cumulative[1], fixture.cumulative[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 1 + 5 * 1.82e-6), fixture.boundary[1], 2e-15);
}

test "REDIST SOC organic concentration boundary resets segment" {
    var fixture = Fixture{};
    fixture.corgc[2] = 2;
    fixture.corgc[3] = 0;
    try build(3, 1, 3, fixture.inputs(), &fixture.cumulative, &fixture.boundary, &fixture.organic_flag, &fixture.water_flag);
    try std.testing.expectApproxEqAbs(@as(f64, 2 * 1.82e-6), fixture.cumulative[2], 1e-18);
    try std.testing.expectEqual(@as(f64, 0), fixture.boundary[2]);
    try std.testing.expectEqual(@as(u8, 2), fixture.water_flag[2]);
}

test "REDIST SOC inactive layer inherits active lower cumulative change" {
    var fixture = Fixture{};
    fixture.dorgc[2] = 0;
    try build(2, 1, 3, fixture.inputs(), &fixture.cumulative, &fixture.boundary, &fixture.organic_flag, &fixture.water_flag);
    try std.testing.expectEqual(fixture.cumulative[3], fixture.cumulative[2]);
    try std.testing.expectEqual(@as(u8, 0), fixture.organic_flag[2]);
}

test "REDIST SOC depth rejects dimensions and invalid geometry" {
    var fixture = Fixture{};
    try std.testing.expectError(error.SoilSocDepthDimensionMismatch, build(2, 3, 1, fixture.inputs(), &fixture.cumulative, &fixture.boundary, &fixture.organic_flag, &fixture.water_flag));
    fixture.area[3] = 0;
    try std.testing.expectError(error.InvalidSoilSocGeometry, build(2, 1, 3, fixture.inputs(), &fixture.cumulative, &fixture.boundary, &fixture.organic_flag, &fixture.water_flag));
}
