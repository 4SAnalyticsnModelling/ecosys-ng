const std = @import("std");

pub const State = struct {
    layer_thickness_m: []const f64, // DLYR(3,L)
    temperature_k: []f64, // TKS
    temperature_c: []f64, // TCS
};

fn finiteSlice(values: []const f64) bool {
    for (values) |value| if (!std.math.isFinite(value)) return false;
    return true;
}

/// Direct translation of REDIST 10670--10677. Zig's exclusive upper bound
/// reproduces Fortran `L=NUI,NU-1` after index conversion.
pub fn fill(
    initial_upper_layer: usize,
    current_upper_layer: usize,
    minimum_layer_thickness_m: f64,
    state: State,
) !void {
    const layer_count = state.layer_thickness_m.len;
    if (layer_count == 0 or initial_upper_layer >= layer_count or current_upper_layer >= layer_count or
        state.temperature_k.len != layer_count or state.temperature_c.len != layer_count)
        return error.ThinLayerTemperatureFillDimensionMismatch;
    if (!finiteSlice(state.layer_thickness_m) or !finiteSlice(state.temperature_k) or
        !finiteSlice(state.temperature_c) or !std.math.isFinite(minimum_layer_thickness_m) or
        minimum_layer_thickness_m < 0)
        return error.InvalidThinLayerTemperatureFillInput;
    for (state.layer_thickness_m) |thickness_m|
        if (thickness_m < 0) return error.InvalidThinLayerTemperatureFillInput;
    if (current_upper_layer <= initial_upper_layer) return;

    const reference_temperature_k = state.temperature_k[current_upper_layer];
    const reference_temperature_c = reference_temperature_k - 273.15;
    if (!std.math.isFinite(reference_temperature_c)) return error.NonFiniteThinLayerTemperatureFillResult;
    for (initial_upper_layer..current_upper_layer) |layer| {
        if (state.layer_thickness_m[layer] <= minimum_layer_thickness_m) {
            state.temperature_k[layer] = reference_temperature_k;
            state.temperature_c[layer] = state.temperature_k[layer] - 273.15;
        }
    }
}

test "REDIST thin-layer temperature fill preserves NUI through NU minus one" {
    const thickness = [_]f64{ 1, 0.1, 0.1001, 0.05, 1 };
    var temperature_k = [_]f64{ 270, 271, 272, 273, 290 };
    var temperature_c = [_]f64{ -3.15, -2.15, -1.15, -0.15, 16.85 };
    try fill(1, 4, 0.1, .{ .layer_thickness_m = &thickness, .temperature_k = &temperature_k, .temperature_c = &temperature_c });
    try std.testing.expectEqual(@as(f64, 290), temperature_k[1]);
    try std.testing.expectApproxEqAbs(@as(f64, 16.85), temperature_c[1], 1e-12);
    try std.testing.expectEqual(@as(f64, 272), temperature_k[2]);
    try std.testing.expectEqual(@as(f64, 290), temperature_k[3]);
    try std.testing.expectEqual(@as(f64, 290), temperature_k[4]);
}

test "REDIST thin-layer fill has exact zero-trip semantics" {
    const thickness = [_]f64{ 0, 0 };
    var temperature_k = [_]f64{ 270, 280 };
    var temperature_c = [_]f64{ -3.15, 6.85 };
    try fill(1, 1, 0.1, .{ .layer_thickness_m = &thickness, .temperature_k = &temperature_k, .temperature_c = &temperature_c });
    try std.testing.expectEqual(@as(f64, 270), temperature_k[0]);
}

test "REDIST thin-layer fill validation is atomic" {
    const thickness = [_]f64{ 0, 0 };
    var temperature_k = [_]f64{ 270, std.math.nan(f64) };
    var temperature_c = [_]f64{ -3.15, 6.85 };
    try std.testing.expectError(error.InvalidThinLayerTemperatureFillInput, fill(0, 1, 0.1, .{
        .layer_thickness_m = &thickness,
        .temperature_k = &temperature_k,
        .temperature_c = &temperature_c,
    }));
    try std.testing.expectEqual(@as(f64, 270), temperature_k[0]);
}
