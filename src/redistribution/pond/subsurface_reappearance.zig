const std = @import("std");

pub const Inputs = struct {
    top_layer: usize, // NU
    bulk_density_megagrams_per_m3: []const f64, // BKDS
    water_volume_m3: []const f64, // VOLW
    ice_volume_m3: []const f64, // VOLI
    air_volume_m3: []const f64, // VOLA
    area_m2: []const f64, // AREA(3)
    minimum_excess_volume_m3: f64, // ZEROS
    zero_tolerance: f64, // ZERO
};

pub const State = struct {
    layer_depth_m: []f64, // DLYR(3)
    cumulative_depth_m: []f64, // CDPTH
    midpoint_depth_m: []f64, // DPTH
    bottom_depth_from_surface_m: []f64, // CDPTHZ
    midpoint_depth_from_surface_m: []f64, // DPTHZ
    relayering_change_m: *f64, // DDLYRX(5)
    pond_boundary_flag: *u8, // IFLGL(L,5)
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8392--8421 (`NN=5`).
pub fn apply(layer: usize, inputs: Inputs, state: State) !void {
    const len = state.layer_depth_m.len;
    if (len < 2 or layer >= len - 1 or inputs.top_layer == 0 or inputs.top_layer >= len or
        inputs.bulk_density_megagrams_per_m3.len != len or inputs.water_volume_m3.len != len or
        inputs.ice_volume_m3.len != len or inputs.air_volume_m3.len != len or inputs.area_m2.len != len or
        state.cumulative_depth_m.len != len or state.midpoint_depth_m.len != len or
        state.bottom_depth_from_surface_m.len != len or state.midpoint_depth_from_surface_m.len != len)
        return error.PondSubsurfaceReappearanceDimensionMismatch;
    inline for (.{ inputs.bulk_density_megagrams_per_m3, inputs.water_volume_m3, inputs.ice_volume_m3, inputs.air_volume_m3, inputs.area_m2, state.layer_depth_m, state.cumulative_depth_m, state.midpoint_depth_m, state.bottom_depth_from_surface_m, state.midpoint_depth_from_surface_m }) |values|
        if (!finiteSlice(values)) return error.InvalidPondSubsurfaceReappearanceInput;
    inline for (.{ inputs.minimum_excess_volume_m3, inputs.zero_tolerance, state.relayering_change_m.* }) |value|
        if (!std.math.isFinite(value)) return error.InvalidPondSubsurfaceReappearanceInput;
    if (inputs.minimum_excess_volume_m3 < 0 or inputs.zero_tolerance < 0 or inputs.area_m2[layer] <= 0)
        return error.InvalidPondSubsurfaceReappearanceInput;

    if (layer >= inputs.top_layer and inputs.bulk_density_megagrams_per_m3[layer] > inputs.zero_tolerance and
        inputs.bulk_density_megagrams_per_m3[layer + 1] <= inputs.zero_tolerance)
    {
        const excess_water_volume_m3 = @max(0.0, inputs.water_volume_m3[layer] +
            inputs.ice_volume_m3[layer] - inputs.air_volume_m3[layer]);
        if (excess_water_volume_m3 > inputs.minimum_excess_volume_m3) {
            state.relayering_change_m.* = -excess_water_volume_m3 / inputs.area_m2[layer];
            state.pond_boundary_flag.* = 1;
            state.layer_depth_m[layer + 1] -= state.relayering_change_m.*;
            state.cumulative_depth_m[layer] = state.cumulative_depth_m[layer + 1] -
                state.layer_depth_m[layer + 1];
            state.midpoint_depth_m[layer + 1] = 0.5 *
                (state.cumulative_depth_m[layer + 1] + state.cumulative_depth_m[layer]);
            state.bottom_depth_from_surface_m[layer + 1] = state.cumulative_depth_m[layer + 1] -
                state.cumulative_depth_m[inputs.top_layer - 1];
            state.midpoint_depth_from_surface_m[layer + 1] = 0.5 *
                (state.bottom_depth_from_surface_m[layer] + state.bottom_depth_from_surface_m[layer + 1]);
        } else {
            state.relayering_change_m.* = 0.0;
            state.pond_boundary_flag.* = 0;
        }
    } else {
        state.relayering_change_m.* = 0.0;
        state.pond_boundary_flag.* = 0;
    }

    inline for (.{ state.layer_depth_m[layer + 1], state.cumulative_depth_m[layer], state.midpoint_depth_m[layer + 1], state.bottom_depth_from_surface_m[layer + 1], state.midpoint_depth_from_surface_m[layer + 1], state.relayering_change_m.* }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePondSubsurfaceReappearanceState;
}

const Fixture = struct {
    density: [4]f64 = .{ 0, 1.2, 1.3, 0 },
    water: [4]f64 = .{ 0, 0.1, 0.8, 0.2 },
    ice: [4]f64 = .{ 0, 0, 0.2, 0 },
    air: [4]f64 = .{ 0, 0, 0.4, 0 },
    area: [4]f64 = .{ 1, 2, 2, 2 },
    depth: [4]f64 = .{ 0, 0.2, 0.3, 0.1 },
    cumulative: [4]f64 = .{ 0, 1, 2, 3 },
    midpoint: [4]f64 = .{ 0, 0.5, 1.5, 2.5 },
    bottom_from_surface: [4]f64 = .{ 0, 1, 2, 3 },
    midpoint_from_surface: [4]f64 = .{ 0, 0.5, 1.5, 2.5 },
    change: f64 = 9,
    flag: u8 = 8,

    fn inputs(self: *Fixture) Inputs {
        return .{
            .top_layer = 1,
            .bulk_density_megagrams_per_m3 = &self.density,
            .water_volume_m3 = &self.water,
            .ice_volume_m3 = &self.ice,
            .air_volume_m3 = &self.air,
            .area_m2 = &self.area,
            .minimum_excess_volume_m3 = 0.01,
            .zero_tolerance = 1.0e-12,
        };
    }

    fn state(self: *Fixture) State {
        return .{
            .layer_depth_m = &self.depth,
            .cumulative_depth_m = &self.cumulative,
            .midpoint_depth_m = &self.midpoint,
            .bottom_depth_from_surface_m = &self.bottom_from_surface,
            .midpoint_depth_from_surface_m = &self.midpoint_from_surface,
            .relayering_change_m = &self.change,
            .pond_boundary_flag = &self.flag,
        };
    }
};

test "REDIST pond subsurface reappearance expands water layer and geometry" {
    var fixture = Fixture{};
    try apply(2, fixture.inputs(), fixture.state());

    try std.testing.expectApproxEqAbs(@as(f64, -0.3), fixture.change, 1.0e-14);
    try std.testing.expectEqual(@as(u8, 1), fixture.flag);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), fixture.depth[3], 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2.6), fixture.cumulative[2], 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2.8), fixture.midpoint[3], 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 3), fixture.bottom_from_surface[3], 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2.5), fixture.midpoint_from_surface[3], 1.0e-14);
}

test "REDIST pond subsurface reappearance requires excess above threshold" {
    var fixture = Fixture{};
    fixture.air[2] = 1;
    try apply(2, fixture.inputs(), fixture.state());
    try std.testing.expectEqual(@as(f64, 0), fixture.change);
    try std.testing.expectEqual(@as(u8, 0), fixture.flag);
    try std.testing.expectEqual(@as(f64, 0.1), fixture.depth[3]);
}

test "REDIST pond subsurface reappearance requires soil over water" {
    var fixture = Fixture{};
    fixture.density[3] = 0.2;
    try apply(2, fixture.inputs(), fixture.state());
    try std.testing.expectEqual(@as(f64, 0), fixture.change);
    try std.testing.expectEqual(@as(u8, 0), fixture.flag);
}

test "REDIST pond subsurface reappearance rejects zero layer area" {
    var fixture = Fixture{};
    fixture.area[2] = 0;
    try std.testing.expectError(
        error.InvalidPondSubsurfaceReappearanceInput,
        apply(2, fixture.inputs(), fixture.state()),
    );
}
