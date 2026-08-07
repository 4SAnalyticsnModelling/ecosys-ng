const std = @import("std");

pub const Inputs = struct {
    initial_top_layer: usize, // NUI
    bulk_density_megagrams_per_m3: []const f64, // BKDS
    surface_water_volume_m3: f64, // VOLW(0)
    surface_ice_volume_m3: f64, // VOLI(0)
    ponding_capacity_m3: f64, // VOLWD
    litter_water_capacity_m3: f64, // VOLWRX
    dry_litter_volume_m3: f64, // VOLR
    initial_pond_bottom_depth_m: f64, // CDPTHI
    minimum_surface_heat_capacity_megajoules_per_k: f64, // VHCPNX
    surface_area_m2: f64, // AREA(3,0)
    zero_tolerance: f64, // ZERO
};

pub const Geometry = struct {
    layer_depth_m: []f64, // DLYR(3)
    cumulative_depth_m: []f64, // CDPTH
    midpoint_depth_m: []f64, // DPTH
    bottom_depth_from_surface_m: []f64, // CDPTHZ
    midpoint_depth_from_surface_m: []f64, // DPTHZ
};

pub const State = struct {
    top_layer: *usize, // NU
    active_surface_layer: *usize, // NUM
    relayering_change_m: *f64, // DDLYRX(3)
    pond_boundary_flag: *u8, // IFLGL(L,3)
    geometry: Geometry,
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8306--8351 (`NN=3`).
pub fn apply(layer: usize, inputs: Inputs, state: State) !void {
    const len = state.geometry.layer_depth_m.len;
    if (len == 0 or layer >= len or inputs.initial_top_layer >= len or state.top_layer.* >= len or
        state.active_surface_layer.* >= len or state.geometry.cumulative_depth_m.len != len or
        state.geometry.midpoint_depth_m.len != len or state.geometry.bottom_depth_from_surface_m.len != len or
        state.geometry.midpoint_depth_from_surface_m.len != len or inputs.bulk_density_megagrams_per_m3.len != len)
        return error.PondSurfaceReappearanceDimensionMismatch;
    inline for (.{ inputs.bulk_density_megagrams_per_m3, state.geometry.layer_depth_m, state.geometry.cumulative_depth_m, state.geometry.midpoint_depth_m, state.geometry.bottom_depth_from_surface_m, state.geometry.midpoint_depth_from_surface_m }) |values|
        if (!finiteSlice(values)) return error.InvalidPondSurfaceReappearanceInput;
    inline for (.{ inputs.surface_water_volume_m3, inputs.surface_ice_volume_m3, inputs.ponding_capacity_m3, inputs.litter_water_capacity_m3, inputs.dry_litter_volume_m3, inputs.initial_pond_bottom_depth_m, inputs.minimum_surface_heat_capacity_megajoules_per_k, inputs.surface_area_m2, inputs.zero_tolerance, state.relayering_change_m.* }) |value|
        if (!std.math.isFinite(value)) return error.InvalidPondSurfaceReappearanceInput;
    if (inputs.surface_area_m2 <= 0 or inputs.zero_tolerance < 0)
        return error.InvalidPondSurfaceReappearanceInput;

    const excess_pond_volume_m3 = @max(0.0, inputs.surface_water_volume_m3 +
        inputs.surface_ice_volume_m3 - inputs.ponding_capacity_m3);
    const old_top_layer = state.top_layer.*;
    if (layer == old_top_layer and state.geometry.cumulative_depth_m[0] > inputs.initial_pond_bottom_depth_m and
        excess_pond_volume_m3 > inputs.minimum_surface_heat_capacity_megajoules_per_k / 4.19)
    {
        if (inputs.bulk_density_megagrams_per_m3[layer] > inputs.zero_tolerance and
            old_top_layer > inputs.initial_top_layer)
        {
            state.relayering_change_m.* = -excess_pond_volume_m3 / inputs.surface_area_m2;
            state.top_layer.* = inputs.initial_top_layer;
            state.active_surface_layer.* = inputs.initial_top_layer;
            state.pond_boundary_flag.* = 1;
            const pond_layer_depth_m = (@max(0.0, inputs.surface_water_volume_m3 +
                inputs.surface_ice_volume_m3 - inputs.litter_water_capacity_m3) +
                inputs.dry_litter_volume_m3) / inputs.surface_area_m2;
            state.geometry.layer_depth_m[0] = pond_layer_depth_m + state.relayering_change_m.*;
            state.geometry.layer_depth_m[inputs.initial_top_layer] -= state.relayering_change_m.*;

            if (layer > 2) {
                var deeper_layer = layer - 2;
                while (deeper_layer >= inputs.initial_top_layer) {
                    state.geometry.cumulative_depth_m[deeper_layer] =
                        state.geometry.cumulative_depth_m[layer - 1];
                    if (deeper_layer == inputs.initial_top_layer) break;
                    deeper_layer -= 1;
                }
            }
            state.geometry.cumulative_depth_m[0] =
                state.geometry.cumulative_depth_m[inputs.initial_top_layer] -
                state.geometry.layer_depth_m[inputs.initial_top_layer];
            state.geometry.midpoint_depth_m[inputs.initial_top_layer] = 0.5 *
                (state.geometry.cumulative_depth_m[inputs.initial_top_layer] +
                    state.geometry.cumulative_depth_m[0]);
            state.geometry.bottom_depth_from_surface_m[inputs.initial_top_layer] =
                state.geometry.layer_depth_m[inputs.initial_top_layer];
            state.geometry.midpoint_depth_from_surface_m[inputs.initial_top_layer] = 0.5 *
                state.geometry.bottom_depth_from_surface_m[inputs.initial_top_layer];
        } else {
            state.relayering_change_m.* = 0.0;
            state.pond_boundary_flag.* = 0;
        }
    } else {
        state.relayering_change_m.* = 0.0;
        state.pond_boundary_flag.* = 0;
    }

    inline for (.{ state.geometry.layer_depth_m[0], state.geometry.layer_depth_m[state.top_layer.*], state.geometry.cumulative_depth_m[0], state.relayering_change_m.* }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePondSurfaceReappearanceState;
}

const Fixture = struct {
    density: [5]f64 = .{ 0, 1.2, 1.3, 1.4, 1.5 },
    depth: [5]f64 = .{ 0.1, 0.2, 0.3, 0.4, 0.5 },
    cumulative: [5]f64 = .{ 0.3, 1, 2, 3, 4 },
    midpoint: [5]f64 = .{ 0, 0, 0, 0, 0 },
    bottom_from_surface: [5]f64 = .{ 0, 0, 0, 0, 0 },
    midpoint_from_surface: [5]f64 = .{ 0, 0, 0, 0, 0 },
    top_layer: usize = 3,
    active_surface_layer: usize = 3,
    change: f64 = 8,
    flag: u8 = 9,

    fn inputs(self: *Fixture) Inputs {
        return .{
            .initial_top_layer = 1,
            .bulk_density_megagrams_per_m3 = &self.density,
            .surface_water_volume_m3 = 0.8,
            .surface_ice_volume_m3 = 0.2,
            .ponding_capacity_m3 = 0.2,
            .litter_water_capacity_m3 = 0.3,
            .dry_litter_volume_m3 = 0.1,
            .initial_pond_bottom_depth_m = 0.1,
            .minimum_surface_heat_capacity_megajoules_per_k = 0.419,
            .surface_area_m2 = 2,
            .zero_tolerance = 1.0e-12,
        };
    }

    fn state(self: *Fixture) State {
        return .{
            .top_layer = &self.top_layer,
            .active_surface_layer = &self.active_surface_layer,
            .relayering_change_m = &self.change,
            .pond_boundary_flag = &self.flag,
            .geometry = .{
                .layer_depth_m = &self.depth,
                .cumulative_depth_m = &self.cumulative,
                .midpoint_depth_m = &self.midpoint,
                .bottom_depth_from_surface_m = &self.bottom_from_surface,
                .midpoint_depth_from_surface_m = &self.midpoint_from_surface,
            },
        };
    }
};

test "REDIST pond surface reappearance restores initial top layer and geometry" {
    var fixture = Fixture{};
    try apply(3, fixture.inputs(), fixture.state());

    try std.testing.expectEqual(@as(usize, 1), fixture.top_layer);
    try std.testing.expectEqual(@as(usize, 1), fixture.active_surface_layer);
    try std.testing.expectEqual(@as(u8, 1), fixture.flag);
    try std.testing.expectApproxEqAbs(@as(f64, -0.4), fixture.change, 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), fixture.depth[0], 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), fixture.depth[1], 1.0e-14);
    try std.testing.expectEqual(@as(f64, 2), fixture.cumulative[1]);
    try std.testing.expectApproxEqAbs(@as(f64, 1.4), fixture.cumulative[0], 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 1.7), fixture.midpoint[1], 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), fixture.bottom_from_surface[1], 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), fixture.midpoint_from_surface[1], 1.0e-14);
}

test "REDIST pond surface reappearance clears outputs when outer condition fails" {
    var fixture = Fixture{};
    fixture.cumulative[0] = 0.05;
    try apply(3, fixture.inputs(), fixture.state());
    try std.testing.expectEqual(@as(f64, 0), fixture.change);
    try std.testing.expectEqual(@as(u8, 0), fixture.flag);
    try std.testing.expectEqual(@as(usize, 3), fixture.top_layer);
}

test "REDIST pond surface reappearance clears outputs when layer is not soil" {
    var fixture = Fixture{};
    fixture.density[3] = 0;
    try apply(3, fixture.inputs(), fixture.state());
    try std.testing.expectEqual(@as(f64, 0), fixture.change);
    try std.testing.expectEqual(@as(u8, 0), fixture.flag);
}

test "REDIST pond surface reappearance rejects zero surface area" {
    var fixture = Fixture{};
    var inputs = fixture.inputs();
    inputs.surface_area_m2 = 0;
    try std.testing.expectError(error.InvalidPondSurfaceReappearanceInput, apply(3, inputs, fixture.state()));
}
