const std = @import("std");

pub const Inputs = struct {
    bulk_density_megagrams_per_m3: []const f64, // BKDS
    area_m2: []const f64, // AREA(3)
    micropore_fraction: []const f64, // FMPR
    maximum_disappearing_depth_m: f64, // ZERO2
    zero_tolerance: f64, // ZERO
};

pub const State = struct {
    layer_depth_m: []f64, // DLYR(3)
    total_volume_m3: []f64, // VOLT
    micropore_volume_m3: []f64, // VOLX
    relayering_change_m: *f64, // DDLYRX(4)
    pond_boundary_flag: *u8, // IFLGL(L,4)
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 8356--8371 (`NN=4`).
pub fn apply(layer: usize, inputs: Inputs, state: State) !void {
    const len = state.layer_depth_m.len;
    if (len < 2 or layer >= len - 1 or inputs.bulk_density_megagrams_per_m3.len != len or
        inputs.area_m2.len != len or inputs.micropore_fraction.len != len or
        state.total_volume_m3.len != len or state.micropore_volume_m3.len != len)
        return error.PondSubsurfaceDisappearanceDimensionMismatch;
    inline for (.{ inputs.bulk_density_megagrams_per_m3, inputs.area_m2, inputs.micropore_fraction, state.layer_depth_m, state.total_volume_m3, state.micropore_volume_m3 }) |values|
        if (!finiteSlice(values)) return error.InvalidPondSubsurfaceDisappearanceInput;
    inline for (.{ inputs.maximum_disappearing_depth_m, inputs.zero_tolerance, state.relayering_change_m.* }) |value|
        if (!std.math.isFinite(value)) return error.InvalidPondSubsurfaceDisappearanceInput;
    if (inputs.maximum_disappearing_depth_m < 0 or inputs.zero_tolerance < 0 or
        inputs.area_m2[layer] < 0 or inputs.micropore_fraction[layer + 1] < 0)
        return error.InvalidPondSubsurfaceDisappearanceInput;

    const subsurface_layer = layer + 1;
    if (inputs.bulk_density_megagrams_per_m3[layer] > inputs.zero_tolerance and
        inputs.bulk_density_megagrams_per_m3[subsurface_layer] <= inputs.zero_tolerance and
        state.layer_depth_m[subsurface_layer] <= inputs.maximum_disappearing_depth_m and
        state.layer_depth_m[subsurface_layer] > 0.0)
    {
        state.relayering_change_m.* = state.layer_depth_m[subsurface_layer];
        state.pond_boundary_flag.* = 1;
        state.layer_depth_m[subsurface_layer] = 0.0;
        state.total_volume_m3[subsurface_layer] = inputs.area_m2[layer] *
            state.layer_depth_m[subsurface_layer];
        state.micropore_volume_m3[subsurface_layer] = state.total_volume_m3[subsurface_layer] *
            inputs.micropore_fraction[subsurface_layer];
    } else {
        state.relayering_change_m.* = 0.0;
        state.pond_boundary_flag.* = 0;
    }
}

const Fixture = struct {
    density: [3]f64 = .{ 0, 1.2, 0 },
    area: [3]f64 = .{ 1, 2, 3 },
    micropore_fraction: [3]f64 = .{ 0.8, 0.7, 0.6 },
    depth: [3]f64 = .{ 0, 0.2, 0.005 },
    total_volume: [3]f64 = .{ 0, 0.4, 0.015 },
    micropore_volume: [3]f64 = .{ 0, 0.28, 0.009 },
    change: f64 = -1,
    flag: u8 = 8,

    fn inputs(self: *Fixture) Inputs {
        return .{
            .bulk_density_megagrams_per_m3 = &self.density,
            .area_m2 = &self.area,
            .micropore_fraction = &self.micropore_fraction,
            .maximum_disappearing_depth_m = 0.01,
            .zero_tolerance = 1.0e-12,
        };
    }

    fn state(self: *Fixture) State {
        return .{
            .layer_depth_m = &self.depth,
            .total_volume_m3 = &self.total_volume,
            .micropore_volume_m3 = &self.micropore_volume,
            .relayering_change_m = &self.change,
            .pond_boundary_flag = &self.flag,
        };
    }
};

test "REDIST pond subsurface disappearance removes shallow positive water layer" {
    var fixture = Fixture{};
    try apply(1, fixture.inputs(), fixture.state());
    try std.testing.expectEqual(@as(f64, 0.005), fixture.change);
    try std.testing.expectEqual(@as(u8, 1), fixture.flag);
    try std.testing.expectEqual(@as(f64, 0), fixture.depth[2]);
    try std.testing.expectEqual(@as(f64, 0), fixture.total_volume[2]);
    try std.testing.expectEqual(@as(f64, 0), fixture.micropore_volume[2]);
}

test "REDIST pond subsurface disappearance requires strictly positive depth" {
    var fixture = Fixture{};
    fixture.depth[2] = 0;
    try apply(1, fixture.inputs(), fixture.state());
    try std.testing.expectEqual(@as(f64, 0), fixture.change);
    try std.testing.expectEqual(@as(u8, 0), fixture.flag);
}

test "REDIST pond subsurface disappearance rejects a soil subsurface layer" {
    var fixture = Fixture{};
    fixture.density[2] = 0.1;
    try apply(1, fixture.inputs(), fixture.state());
    try std.testing.expectEqual(@as(f64, 0), fixture.change);
    try std.testing.expectEqual(@as(u8, 0), fixture.flag);
    try std.testing.expectEqual(@as(f64, 0.005), fixture.depth[2]);
}

test "REDIST pond subsurface disappearance rejects non-finite depth" {
    var fixture = Fixture{};
    fixture.depth[2] = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidPondSubsurfaceDisappearanceInput,
        apply(1, fixture.inputs(), fixture.state()),
    );
}
