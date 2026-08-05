const std = @import("std");

pub const Inputs = struct {
    initial_depth_m: []const f64, // DLYRI
    bulk_density_megagrams_per_m3: []const f64, // BKDS
    water_boundary_flag: []const u8, // IFLGL(:,1)
    pond_profile_flag: []const u8, // IFLGK
    previous_cumulative_depth_m: []const f64, // CDPTHX
    freeze_excluded_cumulative_depth_m: []const f64, // CDPTHY
    freeze_thaw_changed: bool, // IFLGM == 1
    minimum_layer_depth_m: f64, // DLYRM
    zero_tolerance: f64,
};

pub const Geometry = struct {
    cumulative_depth_m: []f64, // CDPTH
    current_depth_m: []f64, // DLYR(3)
    midpoint_depth_m: []f64, // DPTH
    bottom_depth_from_surface_m: []f64, // CDPTHZ
    midpoint_depth_from_surface_m: []f64, // DPTHZ
    relayering_change_m: []f64, // DDLYRX(1) per layer
    restoration_change_m: []f64, // DDLYRY(L)
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8172--8222 (`NN=1`).
pub fn apply(
    top_layer: usize,
    bottom_layer: usize,
    inputs: Inputs,
    geometry: Geometry,
) !void {
    const len = geometry.cumulative_depth_m.len;
    if (len == 0 or top_layer == 0 or top_layer >= bottom_layer or bottom_layer >= len or
        inputs.initial_depth_m.len != len or inputs.bulk_density_megagrams_per_m3.len != len or
        inputs.water_boundary_flag.len != len or inputs.pond_profile_flag.len != len or
        inputs.previous_cumulative_depth_m.len != len or inputs.freeze_excluded_cumulative_depth_m.len != len or
        geometry.current_depth_m.len != len or geometry.midpoint_depth_m.len != len or
        geometry.bottom_depth_from_surface_m.len != len or geometry.midpoint_depth_from_surface_m.len != len or
        geometry.relayering_change_m.len != len or geometry.restoration_change_m.len != len)
        return error.SoilPrimaryRelayeringDimensionMismatch;
    inline for (.{ inputs.initial_depth_m, inputs.bulk_density_megagrams_per_m3, inputs.previous_cumulative_depth_m, inputs.freeze_excluded_cumulative_depth_m, geometry.cumulative_depth_m, geometry.current_depth_m, geometry.midpoint_depth_m, geometry.bottom_depth_from_surface_m, geometry.midpoint_depth_from_surface_m, geometry.relayering_change_m, geometry.restoration_change_m }) |values|
        if (!finiteSlice(values)) return error.InvalidSoilPrimaryRelayeringInput;
    if (!std.math.isFinite(inputs.minimum_layer_depth_m) or
        !std.math.isFinite(inputs.zero_tolerance) or inputs.zero_tolerance < 0)
        return error.InvalidSoilPrimaryRelayeringInput;

    for (top_layer..bottom_layer) |layer| {
        const depth_before_reset = geometry.cumulative_depth_m[layer] - geometry.cumulative_depth_m[layer - 1];
        var relayering_change: f64 = undefined;
        var restoration_change: f64 = undefined;
        if (inputs.water_boundary_flag[layer] == 0 and inputs.water_boundary_flag[layer + 1] == 1) {
            relayering_change = 0.0;
            restoration_change = if (inputs.bulk_density_megagrams_per_m3[layer] <= inputs.zero_tolerance)
                inputs.initial_depth_m[layer] - depth_before_reset
            else
                geometry.current_depth_m[layer] - depth_before_reset;
        } else if (inputs.water_boundary_flag[layer] == 2 and
            (inputs.water_boundary_flag[layer + 1] == 0 or depth_before_reset <= inputs.initial_depth_m[layer]))
        {
            relayering_change = 0.0;
            restoration_change = 0.0;
        } else if (inputs.freeze_thaw_changed and inputs.bulk_density_megagrams_per_m3[layer] > inputs.zero_tolerance) {
            relayering_change = inputs.freeze_excluded_cumulative_depth_m[layer] -
                inputs.previous_cumulative_depth_m[layer];
            restoration_change = inputs.initial_depth_m[layer] - depth_before_reset;
        } else {
            relayering_change = inputs.initial_depth_m[layer] - depth_before_reset;
            restoration_change = relayering_change;
        }
        if (inputs.bulk_density_megagrams_per_m3[layer] > inputs.zero_tolerance and inputs.pond_profile_flag[layer] == 1)
            relayering_change = geometry.cumulative_depth_m[layer] - inputs.previous_cumulative_depth_m[layer];

        if (geometry.current_depth_m[layer] > inputs.minimum_layer_depth_m or
            inputs.bulk_density_megagrams_per_m3[top_layer] <= inputs.zero_tolerance)
        {
            geometry.cumulative_depth_m[layer] = geometry.cumulative_depth_m[layer] + restoration_change;
        } else {
            geometry.cumulative_depth_m[layer] = geometry.cumulative_depth_m[layer - 1] + restoration_change;
            if (inputs.bulk_density_megagrams_per_m3[layer] > inputs.zero_tolerance) {
                relayering_change = if (layer == bottom_layer - 1)
                    geometry.cumulative_depth_m[layer] - inputs.previous_cumulative_depth_m[layer]
                else
                    inputs.freeze_excluded_cumulative_depth_m[layer] - inputs.previous_cumulative_depth_m[layer];
            }
        }
        geometry.relayering_change_m[layer] = relayering_change;
        geometry.restoration_change_m[layer] = restoration_change;
        geometry.current_depth_m[layer] = geometry.cumulative_depth_m[layer] - geometry.cumulative_depth_m[layer - 1];
        geometry.midpoint_depth_m[layer] = 0.5 * (geometry.cumulative_depth_m[layer] + geometry.cumulative_depth_m[layer - 1]);
        geometry.bottom_depth_from_surface_m[layer] = geometry.cumulative_depth_m[layer] -
            geometry.cumulative_depth_m[top_layer - 1];
        if (layer == bottom_layer - 1) {
            geometry.current_depth_m[layer + 1] = geometry.cumulative_depth_m[layer + 1] - geometry.cumulative_depth_m[layer];
            geometry.midpoint_depth_m[layer + 1] = 0.5 * (geometry.cumulative_depth_m[layer + 1] + geometry.cumulative_depth_m[layer]);
            geometry.bottom_depth_from_surface_m[layer + 1] = geometry.cumulative_depth_m[layer + 1] -
                geometry.cumulative_depth_m[top_layer - 1];
        }
        geometry.midpoint_depth_from_surface_m[layer] = if (layer == top_layer)
            0.5 * geometry.bottom_depth_from_surface_m[layer]
        else
            0.5 * (geometry.bottom_depth_from_surface_m[layer] + geometry.bottom_depth_from_surface_m[layer - 1]);
        inline for (.{ geometry.cumulative_depth_m[layer], geometry.current_depth_m[layer], geometry.midpoint_depth_m[layer], geometry.bottom_depth_from_surface_m[layer], geometry.midpoint_depth_from_surface_m[layer], relayering_change, restoration_change }) |value|
            if (!std.math.isFinite(value)) return error.NonFiniteSoilPrimaryRelayeringGeometry;
    }
}

const Fixture = struct {
    initial: [4]f64 = .{ 0, 2, 2, 2 },
    density: [4]f64 = .{ 0, 1, 1, 1 },
    flags: [4]u8 = .{ 0, 0, 0, 0 },
    pond: [4]u8 = .{ 0, 0, 0, 0 },
    previous: [4]f64 = .{ 0, 1, 3, 6 },
    excluded: [4]f64 = .{ 0, 1, 3, 6 },
    cumulative: [4]f64 = .{ 0, 1, 3, 6 },
    depth: [4]f64 = .{ 0, 1, 2, 3 },
    midpoint: [4]f64 = .{ 0, 0, 0, 0 },
    bottom_surface: [4]f64 = .{ 0, 0, 0, 0 },
    midpoint_surface: [4]f64 = .{ 0, 0, 0, 0 },
    change: [4]f64 = .{ 0, 0, 0, 0 },
    restore: [4]f64 = .{ 0, 0, 0, 0 },
    fn inputs(self: *Fixture) Inputs {
        return .{ .initial_depth_m = &self.initial, .bulk_density_megagrams_per_m3 = &self.density, .water_boundary_flag = &self.flags, .pond_profile_flag = &self.pond, .previous_cumulative_depth_m = &self.previous, .freeze_excluded_cumulative_depth_m = &self.excluded, .freeze_thaw_changed = false, .minimum_layer_depth_m = 0.1, .zero_tolerance = 0 };
    }
    fn geometry(self: *Fixture) Geometry {
        return .{ .cumulative_depth_m = &self.cumulative, .current_depth_m = &self.depth, .midpoint_depth_m = &self.midpoint, .bottom_depth_from_surface_m = &self.bottom_surface, .midpoint_depth_from_surface_m = &self.midpoint_surface, .relayering_change_m = &self.change, .restoration_change_m = &self.restore };
    }
};

test "REDIST primary geometry restores runtime layers and derived depths" {
    var fixture = Fixture{};
    try apply(1, 3, fixture.inputs(), fixture.geometry());
    try std.testing.expectEqual(@as(f64, 2), fixture.cumulative[1]);
    try std.testing.expectEqual(@as(f64, 4), fixture.cumulative[2]);
    try std.testing.expectEqual(@as(f64, 2), fixture.depth[1]);
    try std.testing.expectEqual(@as(f64, 2), fixture.depth[3]);
    try std.testing.expectEqual(@as(f64, 1), fixture.midpoint_surface[1]);
}

test "REDIST primary geometry freeze path uses CDPTHY minus CDPTHX" {
    var fixture = Fixture{};
    fixture.excluded[1] = 1.5;
    var inputs = fixture.inputs();
    inputs.freeze_thaw_changed = true;
    try apply(1, 3, inputs, fixture.geometry());
    try std.testing.expectEqual(@as(f64, 0.5), fixture.change[1]);
}

test "REDIST primary geometry pond boundary special case suppresses relayering" {
    var fixture = Fixture{};
    fixture.flags[1] = 2;
    fixture.flags[2] = 0;
    try apply(1, 3, fixture.inputs(), fixture.geometry());
    try std.testing.expectEqual(@as(f64, 0), fixture.change[1]);
    try std.testing.expectEqual(@as(f64, 0), fixture.restore[1]);
}

test "REDIST primary geometry rejects dimensions and non-finite input" {
    var fixture = Fixture{};
    try std.testing.expectError(error.SoilPrimaryRelayeringDimensionMismatch, apply(3, 1, fixture.inputs(), fixture.geometry()));
    fixture.cumulative[1] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilPrimaryRelayeringInput, apply(1, 3, fixture.inputs(), fixture.geometry()));
}
