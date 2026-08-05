const std = @import("std");

pub const Inputs = struct {
    top_layer: usize, // NU
    saved_top_layer: usize, // NUX
    bulk_density_megagrams_per_m3: []const f64, // BKDS
    layer_depth_m: []const f64, // DLYR(3)
    water_volume_m3: []const f64, // VOLW
    ice_volume_m3: []const f64, // VOLI
    area_m2: []const f64, // AREA(3)
    primary_boundary_flag: []const u8, // IFLGL(:,1)
    surface_disappearance_flag: []const u8, // IFLGL(:,2)
    surface_reappearance_flag: []const u8, // IFLGL(:,3)
    subsurface_disappearance_flag: []const u8, // IFLGL(:,4)
    disturbance_flag: []const u8, // IFLGO
    minimum_layer_depth_m: f64, // DLYRM
    zero_tolerance: f64, // ZERO
};

pub const Transfer = struct {
    source_layer: usize, // L0
    destination_layer: usize, // L1
    state_fraction: f64, // FX
    organic_fraction: f64, // FO
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

fn waterDepthM(inputs: Inputs, layer: usize) !f64 {
    if (inputs.area_m2[layer] <= 0) return error.InvalidLayerTransferArea;
    return (inputs.water_volume_m3[layer] + inputs.ice_volume_m3[layer]) / inputs.area_m2[layer];
}

fn fractionsFromPositiveChange(
    layer: usize,
    relayering_change_m: f64,
    disturbance_depth_change_m: f64,
    inputs: Inputs,
) !Transfer {
    const destination_layer = if (inputs.surface_disappearance_flag[layer] == 0)
        layer
    else
        inputs.top_layer;
    const source_layer = if (inputs.surface_disappearance_flag[layer] == 0)
        layer + 1
    else
        inputs.saved_top_layer;

    if ((inputs.bulk_density_megagrams_per_m3[layer] <= inputs.zero_tolerance and
        inputs.primary_boundary_flag[layer] == 2) or
        (inputs.layer_depth_m[source_layer] <= inputs.minimum_layer_depth_m and
            inputs.disturbance_flag[layer] == 1))
        return .{ .source_layer = source_layer, .destination_layer = destination_layer, .state_fraction = 1, .organic_fraction = 1 };

    var state_fraction: f64 = undefined;
    var organic_fraction: f64 = undefined;
    if (inputs.bulk_density_megagrams_per_m3[source_layer] <= inputs.zero_tolerance) {
        const water_depth_m = try waterDepthM(inputs, source_layer);
        if (water_depth_m > inputs.zero_tolerance) {
            state_fraction = @min(1.0, relayering_change_m / water_depth_m);
            organic_fraction = @min(1.0, (relayering_change_m - disturbance_depth_change_m) / water_depth_m);
        } else {
            state_fraction = 1.0;
            organic_fraction = 1.0;
        }
    } else if (inputs.bulk_density_megagrams_per_m3[destination_layer] <= inputs.zero_tolerance and
        inputs.subsurface_disappearance_flag[destination_layer] == 0)
    {
        state_fraction = 0.0;
        organic_fraction = 0.0;
    } else if (inputs.layer_depth_m[source_layer] > inputs.minimum_layer_depth_m) {
        state_fraction = @min(1.0, relayering_change_m / inputs.layer_depth_m[source_layer]);
        organic_fraction = @min(1.0, (relayering_change_m - disturbance_depth_change_m) /
            inputs.layer_depth_m[source_layer]);
    } else {
        state_fraction = 1.0;
        organic_fraction = 1.0;
    }
    return .{ .source_layer = source_layer, .destination_layer = destination_layer, .state_fraction = state_fraction, .organic_fraction = organic_fraction };
}

fn fractionsFromNegativeChange(
    layer: usize,
    relayering_change_m: f64,
    disturbance_depth_change_m: f64,
    inputs: Inputs,
) !Transfer {
    const destination_layer = if (inputs.surface_reappearance_flag[layer] == 0)
        layer + 1
    else
        inputs.top_layer;
    const source_layer = if (inputs.surface_reappearance_flag[layer] == 0)
        layer
    else
        0;

    var state_fraction: f64 = undefined;
    var organic_fraction: f64 = undefined;
    if (inputs.bulk_density_megagrams_per_m3[source_layer] <= inputs.zero_tolerance or
        inputs.surface_reappearance_flag[layer] == 1)
    {
        const water_depth_m = try waterDepthM(inputs, source_layer);
        if (water_depth_m > inputs.zero_tolerance) {
            state_fraction = @min(1.0, -relayering_change_m / water_depth_m);
            organic_fraction = @min(1.0, -(relayering_change_m - disturbance_depth_change_m) / water_depth_m);
        } else {
            state_fraction = 1.0;
            organic_fraction = 1.0;
        }
    } else if (inputs.layer_depth_m[source_layer] > inputs.minimum_layer_depth_m) {
        state_fraction = @min(1.0, -relayering_change_m / inputs.layer_depth_m[source_layer]);
        organic_fraction = @min(1.0, -(relayering_change_m - disturbance_depth_change_m) /
            inputs.layer_depth_m[source_layer]);
    } else {
        state_fraction = 1.0;
        organic_fraction = 1.0;
    }
    return .{ .source_layer = source_layer, .destination_layer = destination_layer, .state_fraction = state_fraction, .organic_fraction = organic_fraction };
}

/// Direct translation of REDIST 8438--8506.
pub fn calculate(
    layer: usize,
    relayering_change_m: f64,
    disturbance_depth_change_m: f64,
    inputs: Inputs,
) !?Transfer {
    const len = inputs.layer_depth_m.len;
    if (len < 2 or layer >= len - 1 or inputs.top_layer >= len or inputs.saved_top_layer >= len or
        inputs.bulk_density_megagrams_per_m3.len != len or inputs.water_volume_m3.len != len or
        inputs.ice_volume_m3.len != len or inputs.area_m2.len != len or inputs.primary_boundary_flag.len != len or
        inputs.surface_disappearance_flag.len != len or inputs.surface_reappearance_flag.len != len or
        inputs.subsurface_disappearance_flag.len != len or inputs.disturbance_flag.len != len)
        return error.LayerTransferFractionDimensionMismatch;
    inline for (.{ inputs.bulk_density_megagrams_per_m3, inputs.layer_depth_m, inputs.water_volume_m3, inputs.ice_volume_m3, inputs.area_m2 }) |values|
        if (!finiteSlice(values)) return error.InvalidLayerTransferFractionInput;
    inline for (.{ inputs.minimum_layer_depth_m, inputs.zero_tolerance, relayering_change_m, disturbance_depth_change_m }) |value|
        if (!std.math.isFinite(value)) return error.InvalidLayerTransferFractionInput;
    if (inputs.minimum_layer_depth_m < 0 or inputs.zero_tolerance < 0)
        return error.InvalidLayerTransferFractionInput;

    if (@abs(relayering_change_m) <= inputs.zero_tolerance) return null;
    const transfer = if (relayering_change_m > 0)
        try fractionsFromPositiveChange(layer, relayering_change_m, disturbance_depth_change_m, inputs)
    else
        try fractionsFromNegativeChange(layer, relayering_change_m, disturbance_depth_change_m, inputs);
    if (!std.math.isFinite(transfer.state_fraction) or !std.math.isFinite(transfer.organic_fraction))
        return error.NonFiniteLayerTransferFraction;
    return transfer;
}

const Fixture = struct {
    density: [4]f64 = .{ 0, 1, 0, 1 },
    depth: [4]f64 = .{ 0.1, 0.4, 0.2, 0.5 },
    water: [4]f64 = .{ 0.4, 0.1, 0.6, 0.1 },
    ice: [4]f64 = .{ 0, 0, 0.2, 0 },
    area: [4]f64 = .{ 2, 2, 2, 2 },
    flag1: [4]u8 = .{ 0, 0, 0, 0 },
    flag2: [4]u8 = .{ 0, 0, 0, 0 },
    flag3: [4]u8 = .{ 0, 0, 0, 0 },
    flag4: [4]u8 = .{ 0, 0, 0, 0 },
    disturbance: [4]u8 = .{ 0, 0, 0, 0 },

    fn inputs(self: *Fixture) Inputs {
        return .{
            .top_layer = 1,
            .saved_top_layer = 1,
            .bulk_density_megagrams_per_m3 = &self.density,
            .layer_depth_m = &self.depth,
            .water_volume_m3 = &self.water,
            .ice_volume_m3 = &self.ice,
            .area_m2 = &self.area,
            .primary_boundary_flag = &self.flag1,
            .surface_disappearance_flag = &self.flag2,
            .surface_reappearance_flag = &self.flag3,
            .subsurface_disappearance_flag = &self.flag4,
            .disturbance_flag = &self.disturbance,
            .minimum_layer_depth_m = 0.01,
            .zero_tolerance = 1.0e-12,
        };
    }
};

test "REDIST positive transfer uses pond water depth fractions" {
    var fixture = Fixture{};
    const transfer = (try calculate(1, 0.1, 0.02, fixture.inputs())).?;
    try std.testing.expectEqual(@as(usize, 2), transfer.source_layer);
    try std.testing.expectEqual(@as(usize, 1), transfer.destination_layer);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), transfer.state_fraction, 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), transfer.organic_fraction, 1.0e-14);
}

test "REDIST positive soil transfer into water is blocked" {
    var fixture = Fixture{};
    fixture.density[2] = 1;
    fixture.density[1] = 0;
    const transfer = (try calculate(1, 0.1, 0, fixture.inputs())).?;
    try std.testing.expectEqual(@as(f64, 0), transfer.state_fraction);
    try std.testing.expectEqual(@as(f64, 0), transfer.organic_fraction);
}

test "REDIST negative reappearance routes surface pond to restored top" {
    var fixture = Fixture{};
    fixture.flag3[2] = 1;
    const transfer = (try calculate(2, -0.05, 0, fixture.inputs())).?;
    try std.testing.expectEqual(@as(usize, 0), transfer.source_layer);
    try std.testing.expectEqual(@as(usize, 1), transfer.destination_layer);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), transfer.state_fraction, 1.0e-14);
}

test "REDIST disturbance condition transfers full source layer" {
    var fixture = Fixture{};
    fixture.depth[2] = 0.005;
    fixture.disturbance[1] = 1;
    const transfer = (try calculate(1, 0.002, 0.001, fixture.inputs())).?;
    try std.testing.expectEqual(@as(f64, 1), transfer.state_fraction);
    try std.testing.expectEqual(@as(f64, 1), transfer.organic_fraction);
}

test "REDIST negligible layer change has no transfer" {
    var fixture = Fixture{};
    try std.testing.expect((try calculate(1, 1.0e-13, 0, fixture.inputs())) == null);
}
