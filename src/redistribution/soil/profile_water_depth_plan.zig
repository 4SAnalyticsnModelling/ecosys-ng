const std = @import("std");

pub const UnderlyingLayer = enum(u8) { none = 0, pond = 1, soil = 2 };

pub const Inputs = struct {
    bulk_density_megagrams_per_m3: []const f64, // BKDS
    current_depth_m: []const f64, // DLYR(3)
    initial_depth_m: []const f64, // DLYRI(3)
    water_m3: []const f64, // VOLW
    ice_m3: []const f64, // VOLI
    micropore_fraction: []const f64, // FMPR
    area_m2: []const f64, // AREA(3,LX)
    cumulative_depth_m: []const f64, // CDPTH
    zero_tolerance: f64, // ZERO
};

pub const Outputs = struct {
    previous_cumulative_depth_m: []f64, // CDPTHX
    working_cumulative_depth_m: []f64, // CDPTHY
    water_cumulative_change_m: []f64, // DDLYX(:,1)
    water_boundary_change_m: []f64, // DDLYR(:,1)
    water_boundary: []UnderlyingLayer, // IFLGL(:,1)
    erosion_cumulative_change_m: []f64, // DDLYX(:,3), pond branch
    erosion_boundary_change_m: []f64, // DDLYR(:,3), pond branch
    soc_cumulative_change_m: []f64, // DDLYX(:,4), pond branch
    soc_boundary_change_m: []f64, // DDLYR(:,4), pond branch
    organic_boundary: []u8, // IFLGO
};

fn validSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 7731--7824. Arrays retain Fortran layer
/// numbers, including `NU-1` and `NL+1` boundary slots.
pub fn build(
    disturbance_mode: i8,
    top_layer: usize,
    bottom_layer: usize,
    inputs: Inputs,
    outputs: Outputs,
) !bool {
    const len = inputs.bulk_density_megagrams_per_m3.len;
    if (len == 0 or top_layer == 0 or top_layer > bottom_layer or bottom_layer + 1 >= len or
        inputs.current_depth_m.len != len or inputs.initial_depth_m.len != len or
        inputs.water_m3.len != len or inputs.ice_m3.len != len or
        inputs.micropore_fraction.len != len or inputs.area_m2.len != len or
        inputs.cumulative_depth_m.len != len or
        outputs.previous_cumulative_depth_m.len != len or outputs.working_cumulative_depth_m.len != len or
        outputs.water_cumulative_change_m.len != len or outputs.water_boundary_change_m.len != len or
        outputs.water_boundary.len != len or outputs.erosion_cumulative_change_m.len != len or
        outputs.erosion_boundary_change_m.len != len or outputs.soc_cumulative_change_m.len != len or
        outputs.soc_boundary_change_m.len != len or outputs.organic_boundary.len != len)
        return error.SoilProfileDepthDimensionMismatch;
    if (!validSlice(inputs.bulk_density_megagrams_per_m3) or !validSlice(inputs.current_depth_m) or
        !validSlice(inputs.initial_depth_m) or !validSlice(inputs.water_m3) or
        !validSlice(inputs.ice_m3) or !validSlice(inputs.micropore_fraction) or
        !validSlice(inputs.area_m2) or !validSlice(inputs.cumulative_depth_m) or
        !validSlice(outputs.water_cumulative_change_m) or !validSlice(outputs.erosion_cumulative_change_m) or
        !validSlice(outputs.soc_cumulative_change_m) or
        !std.math.isFinite(inputs.zero_tolerance) or inputs.zero_tolerance < 0)
        return error.InvalidSoilProfileDepthInput;
    if (disturbance_mode < 0) return false;

    var pond_encountered = false;
    var layer = bottom_layer;
    while (true) : (layer -= 1) {
        outputs.previous_cumulative_depth_m[layer] = inputs.cumulative_depth_m[layer];
        outputs.working_cumulative_depth_m[layer] = inputs.cumulative_depth_m[layer];
        if (inputs.bulk_density_megagrams_per_m3[layer] <= inputs.zero_tolerance) {
            pond_encountered = true;
            const fmpr = inputs.micropore_fraction[layer];
            const area = inputs.area_m2[layer];
            if (fmpr <= 0 or area <= 0) return error.InvalidPondGeometry;
            var pond_change_m: f64 = undefined;
            if (inputs.bulk_density_megagrams_per_m3[layer + 1] > inputs.zero_tolerance) {
                pond_change_m = inputs.current_depth_m[layer] -
                    (inputs.water_m3[layer] + inputs.ice_m3[layer]) / (fmpr * area);
                outputs.water_cumulative_change_m[layer] = pond_change_m + outputs.water_cumulative_change_m[layer + 1];
                outputs.water_boundary_change_m[layer] = outputs.water_cumulative_change_m[layer + 1];
                outputs.water_boundary[layer] = .soil;
            } else {
                pond_change_m = inputs.initial_depth_m[layer] -
                    (inputs.water_m3[layer] + inputs.ice_m3[layer]) / (fmpr * area);
                const underlying_water_depth_m = (inputs.water_m3[layer + 1] + inputs.ice_m3[layer + 1]) / area;
                if (pond_change_m < -inputs.zero_tolerance or underlying_water_depth_m > inputs.zero_tolerance) {
                    outputs.water_cumulative_change_m[layer] = pond_change_m + outputs.water_cumulative_change_m[layer + 1];
                    outputs.water_boundary_change_m[layer] = @min(outputs.water_cumulative_change_m[layer + 1], underlying_water_depth_m);
                    outputs.water_boundary[layer] = if (underlying_water_depth_m > inputs.zero_tolerance) .pond else .soil;
                } else {
                    pond_change_m = inputs.current_depth_m[layer] -
                        (inputs.water_m3[layer] + inputs.ice_m3[layer]) / (fmpr * area);
                    outputs.water_cumulative_change_m[layer] = pond_change_m + outputs.water_cumulative_change_m[layer + 1];
                    outputs.water_boundary_change_m[layer] = outputs.water_cumulative_change_m[layer + 1];
                    outputs.water_boundary[layer] = .soil;
                }
            }
            if (layer == top_layer) {
                outputs.water_cumulative_change_m[layer - 1] = outputs.water_cumulative_change_m[layer];
                outputs.water_boundary_change_m[layer - 1] = outputs.water_cumulative_change_m[layer];
                outputs.water_boundary[layer - 1] = .pond;
            }
            outputs.erosion_cumulative_change_m[layer] = outputs.erosion_cumulative_change_m[layer + 1];
            outputs.erosion_boundary_change_m[layer] = outputs.erosion_cumulative_change_m[layer + 1];
            outputs.soc_cumulative_change_m[layer] = outputs.soc_cumulative_change_m[layer + 1];
            outputs.soc_boundary_change_m[layer] = outputs.soc_cumulative_change_m[layer + 1];
            outputs.organic_boundary[layer] = 0;
            if (layer == top_layer) {
                outputs.erosion_cumulative_change_m[layer - 1] = outputs.erosion_cumulative_change_m[layer];
                outputs.erosion_boundary_change_m[layer - 1] = outputs.erosion_cumulative_change_m[layer];
                outputs.soc_cumulative_change_m[layer - 1] = outputs.soc_cumulative_change_m[layer];
                outputs.soc_boundary_change_m[layer - 1] = outputs.soc_cumulative_change_m[layer];
                outputs.organic_boundary[layer - 1] = 1;
            }
        } else {
            outputs.water_cumulative_change_m[layer] = outputs.water_cumulative_change_m[layer + 1];
            outputs.water_boundary_change_m[layer] = outputs.water_cumulative_change_m[layer + 1];
            outputs.water_boundary[layer] = .none;
            if (layer == top_layer) {
                outputs.water_cumulative_change_m[layer - 1] = outputs.water_cumulative_change_m[layer];
                outputs.water_boundary_change_m[layer - 1] = outputs.water_cumulative_change_m[layer];
                outputs.water_boundary[layer - 1] = .none;
            }
        }
        if (layer == top_layer) break;
    }
    inline for (.{ outputs.previous_cumulative_depth_m, outputs.working_cumulative_depth_m, outputs.water_cumulative_change_m, outputs.water_boundary_change_m, outputs.erosion_cumulative_change_m, outputs.erosion_boundary_change_m, outputs.soc_cumulative_change_m, outputs.soc_boundary_change_m }) |values|
        if (!validSlice(values)) return error.NonFiniteSoilProfileDepthPlan;
    return pond_encountered;
}

const Fixture = struct {
    density: [4]f64 = .{ 0, 0, 1, 1 },
    current: [4]f64 = .{ 0, 1, 1, 1 },
    initial: [4]f64 = .{ 0, 2, 1, 1 },
    water: [4]f64 = .{ 0, 0.5, 0, 0 },
    ice: [4]f64 = .{ 0, 0, 0, 0 },
    fmpr: [4]f64 = .{ 1, 1, 1, 1 },
    area: [4]f64 = .{ 1, 1, 1, 1 },
    cumulative: [4]f64 = .{ 0, 1, 2, 3 },
    previous: [4]f64 = .{ 0, 0, 0, 0 },
    working: [4]f64 = .{ 0, 0, 0, 0 },
    water_x: [4]f64 = .{ 0, 0, 0, 0 },
    water_r: [4]f64 = .{ 0, 0, 0, 0 },
    boundary: [4]UnderlyingLayer = .{ .none, .none, .none, .none },
    erosion_x: [4]f64 = .{ 0, 0, 0, 0 },
    erosion_r: [4]f64 = .{ 0, 0, 0, 0 },
    soc_x: [4]f64 = .{ 0, 0, 0, 0 },
    soc_r: [4]f64 = .{ 0, 0, 0, 0 },
    organic: [4]u8 = .{ 0, 0, 0, 0 },
    fn inputs(self: *Fixture) Inputs {
        return .{ .bulk_density_megagrams_per_m3 = &self.density, .current_depth_m = &self.current, .initial_depth_m = &self.initial, .water_m3 = &self.water, .ice_m3 = &self.ice, .micropore_fraction = &self.fmpr, .area_m2 = &self.area, .cumulative_depth_m = &self.cumulative, .zero_tolerance = 1e-12 };
    }
    fn outputs(self: *Fixture) Outputs {
        return .{ .previous_cumulative_depth_m = &self.previous, .working_cumulative_depth_m = &self.working, .water_cumulative_change_m = &self.water_x, .water_boundary_change_m = &self.water_r, .water_boundary = &self.boundary, .erosion_cumulative_change_m = &self.erosion_x, .erosion_boundary_change_m = &self.erosion_r, .soc_cumulative_change_m = &self.soc_x, .soc_boundary_change_m = &self.soc_r, .organic_boundary = &self.organic };
    }
};

test "REDIST soil profile pond over soil uses current depth and soil flag" {
    var fixture = Fixture{};
    const pond = try build(0, 1, 2, fixture.inputs(), fixture.outputs());
    try std.testing.expect(pond);
    try std.testing.expectEqual(@as(f64, 0.5), fixture.water_x[1]);
    try std.testing.expectEqual(UnderlyingLayer.soil, fixture.boundary[1]);
    try std.testing.expectEqual(UnderlyingLayer.pond, fixture.boundary[0]);
    try std.testing.expectEqual(@as(f64, 2), fixture.previous[2]);
}

test "REDIST soil profile soil layers inherit lower water adjustment" {
    var fixture = Fixture{};
    fixture.density[1] = 1;
    fixture.water_x[3] = 0.25;
    const pond = try build(0, 1, 2, fixture.inputs(), fixture.outputs());
    try std.testing.expect(!pond);
    try std.testing.expectEqual(@as(f64, 0.25), fixture.water_x[2]);
    try std.testing.expectEqual(@as(f64, 0.25), fixture.water_x[1]);
    try std.testing.expectEqual(UnderlyingLayer.none, fixture.boundary[0]);
}

test "REDIST soil profile negative disturbance gate preserves outputs" {
    var fixture = Fixture{};
    fixture.water_x[1] = 9;
    const pond = try build(-1, 1, 2, fixture.inputs(), fixture.outputs());
    try std.testing.expect(!pond);
    try std.testing.expectEqual(@as(f64, 9), fixture.water_x[1]);
}

test "REDIST soil profile depth plan rejects dimensions and pond geometry" {
    var fixture = Fixture{};
    var outputs = fixture.outputs();
    outputs.water_boundary_change_m = outputs.water_boundary_change_m[0..2];
    try std.testing.expectError(error.SoilProfileDepthDimensionMismatch, build(0, 1, 2, fixture.inputs(), outputs));
    fixture.fmpr[1] = 0;
    try std.testing.expectError(error.InvalidPondGeometry, build(0, 1, 2, fixture.inputs(), fixture.outputs()));
}
