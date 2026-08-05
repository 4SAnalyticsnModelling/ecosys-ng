const std = @import("std");

pub const Inputs = struct {
    bulk_density_megagrams_per_m3: []const f64, // BKDS
    heat_capacity_megajoules_per_k: []const f64, // VHCP
    minimum_surface_heat_capacity_megajoules_per_k: f64, // VHCPNX
    area_m2: []const f64, // AREA(3)
    micropore_fraction: []const f64, // FMPR
    minimum_layer_depth_m: f64, // DLYRM
    zero_tolerance: f64, // ZERO
};

pub const State = struct {
    top_layer: *usize, // NU
    active_surface_layer: usize, // NUM
    layer_depth_m: []f64, // DLYR(3)
    total_volume_m3: []f64, // VOLT
    micropore_volume_m3: []f64, // VOLX
    relayering_change_m: *f64, // DDLYRX(2)
    pond_boundary_flag: *u8, // IFLGL(L,2)
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8256--8280 (`NN=2`).
///
/// Array indexes retain the Fortran layer numbers. `saved_top_layer` is NUX,
/// captured before the surrounding layer loop begins.
pub fn apply(
    layer: usize,
    saved_top_layer: usize,
    bottom_layer: usize,
    inputs: Inputs,
    state: State,
) !void {
    const len = state.layer_depth_m.len;
    if (len == 0 or saved_top_layer >= len or bottom_layer >= len or layer >= len or
        state.top_layer.* >= len or state.active_surface_layer >= len or
        inputs.bulk_density_megagrams_per_m3.len != len or inputs.heat_capacity_megajoules_per_k.len != len or
        inputs.area_m2.len != len or inputs.micropore_fraction.len != len or
        state.total_volume_m3.len != len or state.micropore_volume_m3.len != len)
        return error.PondSurfaceDisappearanceDimensionMismatch;
    inline for (.{ inputs.bulk_density_megagrams_per_m3, inputs.heat_capacity_megajoules_per_k, inputs.area_m2, inputs.micropore_fraction, state.layer_depth_m, state.total_volume_m3, state.micropore_volume_m3 }) |values|
        if (!finiteSlice(values)) return error.InvalidPondSurfaceDisappearanceInput;
    inline for (.{ inputs.minimum_surface_heat_capacity_megajoules_per_k, inputs.minimum_layer_depth_m, inputs.zero_tolerance, state.relayering_change_m.* }) |value|
        if (!std.math.isFinite(value)) return error.InvalidPondSurfaceDisappearanceInput;
    if (inputs.minimum_layer_depth_m < 0 or inputs.zero_tolerance < 0 or
        inputs.area_m2[saved_top_layer] < 0 or inputs.micropore_fraction[saved_top_layer] < 0)
        return error.InvalidPondSurfaceDisappearanceInput;

    const current_top_layer = state.top_layer.*;
    if ((layer == current_top_layer and inputs.bulk_density_megagrams_per_m3[current_top_layer] <= inputs.zero_tolerance) and
        (inputs.heat_capacity_megajoules_per_k[current_top_layer] <= inputs.minimum_surface_heat_capacity_megajoules_per_k or
            state.active_surface_layer > current_top_layer))
    {
        var deeper_layer = saved_top_layer + 1;
        while (deeper_layer <= bottom_layer) : (deeper_layer += 1) {
            if (state.layer_depth_m[deeper_layer] > inputs.minimum_layer_depth_m) {
                state.top_layer.* = deeper_layer;
                state.relayering_change_m.* = state.layer_depth_m[saved_top_layer];
                state.pond_boundary_flag.* = 1;
                state.layer_depth_m[saved_top_layer] = 0.0;
                if (inputs.bulk_density_megagrams_per_m3[saved_top_layer] <= inputs.zero_tolerance) {
                    state.total_volume_m3[saved_top_layer] = inputs.area_m2[saved_top_layer] *
                        state.layer_depth_m[saved_top_layer];
                    state.micropore_volume_m3[saved_top_layer] = state.total_volume_m3[saved_top_layer] *
                        inputs.micropore_fraction[saved_top_layer];
                }
                break;
            }
        }
    } else {
        state.relayering_change_m.* = 0.0;
        state.pond_boundary_flag.* = 0;
    }

    inline for (.{ state.layer_depth_m[saved_top_layer], state.total_volume_m3[saved_top_layer], state.micropore_volume_m3[saved_top_layer], state.relayering_change_m.* }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePondSurfaceDisappearanceState;
}

const Fixture = struct {
    density: [4]f64 = .{ 0, 0, 1.2, 1.3 },
    heat_capacity: [4]f64 = .{ 0, 0.2, 1, 1 },
    area: [4]f64 = .{ 1, 2, 2, 2 },
    micropore_fraction: [4]f64 = .{ 0.8, 0.75, 0.7, 0.65 },
    depth: [4]f64 = .{ 0, 0.03, 0.001, 0.2 },
    total_volume: [4]f64 = .{ 0, 0.06, 0.002, 0.4 },
    micropore_volume: [4]f64 = .{ 0, 0.045, 0.0014, 0.26 },
    top_layer: usize = 1,
    relayering_change: f64 = -9,
    flag: u8 = 7,

    fn inputs(self: *Fixture) Inputs {
        return .{
            .bulk_density_megagrams_per_m3 = &self.density,
            .heat_capacity_megajoules_per_k = &self.heat_capacity,
            .minimum_surface_heat_capacity_megajoules_per_k = 0.25,
            .area_m2 = &self.area,
            .micropore_fraction = &self.micropore_fraction,
            .minimum_layer_depth_m = 0.01,
            .zero_tolerance = 1.0e-12,
        };
    }

    fn state(self: *Fixture, active_surface_layer: usize) State {
        return .{
            .top_layer = &self.top_layer,
            .active_surface_layer = active_surface_layer,
            .layer_depth_m = &self.depth,
            .total_volume_m3 = &self.total_volume,
            .micropore_volume_m3 = &self.micropore_volume,
            .relayering_change_m = &self.relayering_change,
            .pond_boundary_flag = &self.flag,
        };
    }
};

test "REDIST pond surface disappearance selects first viable deeper layer" {
    var fixture = Fixture{};
    try apply(1, 1, 3, fixture.inputs(), fixture.state(1));

    try std.testing.expectEqual(@as(usize, 3), fixture.top_layer);
    try std.testing.expectEqual(@as(f64, 0.03), fixture.relayering_change);
    try std.testing.expectEqual(@as(u8, 1), fixture.flag);
    try std.testing.expectEqual(@as(f64, 0), fixture.depth[1]);
    try std.testing.expectEqual(@as(f64, 0), fixture.total_volume[1]);
    try std.testing.expectEqual(@as(f64, 0), fixture.micropore_volume[1]);
}

test "REDIST pond surface disappearance false condition clears outputs" {
    var fixture = Fixture{};
    fixture.heat_capacity[1] = 0.3;
    try apply(1, 1, 3, fixture.inputs(), fixture.state(1));

    try std.testing.expectEqual(@as(usize, 1), fixture.top_layer);
    try std.testing.expectEqual(@as(f64, 0), fixture.relayering_change);
    try std.testing.expectEqual(@as(u8, 0), fixture.flag);
    try std.testing.expectEqual(@as(f64, 0.03), fixture.depth[1]);
}

test "REDIST pond surface disappearance preserves outputs when no candidate exists" {
    var fixture = Fixture{};
    fixture.depth[3] = 0.005;
    try apply(1, 1, 3, fixture.inputs(), fixture.state(1));

    try std.testing.expectEqual(@as(usize, 1), fixture.top_layer);
    try std.testing.expectEqual(@as(f64, -9), fixture.relayering_change);
    try std.testing.expectEqual(@as(u8, 7), fixture.flag);
    try std.testing.expectEqual(@as(f64, 0.03), fixture.depth[1]);
}

test "REDIST pond surface disappearance rejects invalid state" {
    var fixture = Fixture{};
    fixture.depth[2] = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidPondSurfaceDisappearanceInput,
        apply(1, 1, 3, fixture.inputs(), fixture.state(1)),
    );
}
