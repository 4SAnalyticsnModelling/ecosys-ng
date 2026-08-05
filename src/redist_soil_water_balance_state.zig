const std = @import("std");

/// Water and ice state for one runtime soil layer (m3).
pub const LayerWaterState = struct {
    micropore_water_m3: f64, // VOLW
    water_vapour_m3: f64, // VOLV
    macropore_water_m3: f64, // VOLWH
    micropore_ice_m3: f64, // VOLI
    macropore_ice_m3: f64, // VOLIH
};

pub const WaterBalance = struct {
    landscape_water_m3: f64, // VOLWSO
    landscape_ice_m3: f64, // VOLISO
    grid_cell_water_m3: f64, // UVOLW
};

fn finiteLayer(layer: LayerWaterState) bool {
    inline for (@typeInfo(LayerWaterState).@"struct".fields) |field|
        if (!std.math.isFinite(@field(layer, field.name))) return false;
    return true;
}

/// Direct translation of REDIST 6677--6684.
///
/// The update is active only for the reference forcing zone (`NFZ == NFH`).
/// Layers are accumulated in ascending runtime order to preserve rounding.
pub fn accumulateLayers(
    balance: *WaterBalance,
    layers: []const LayerWaterState,
    forcing_zone_index: usize,
    reference_forcing_zone_index: usize,
    ice_water_density_ratio: f64,
) !void {
    if (layers.len == 0) return error.SoilWaterBalanceDimensionMismatch;
    inline for (@typeInfo(WaterBalance).@"struct".fields) |field|
        if (!std.math.isFinite(@field(balance, field.name)))
            return error.InvalidSoilWaterBalance;
    if (!std.math.isFinite(ice_water_density_ratio) or ice_water_density_ratio <= 0)
        return error.InvalidIceWaterDensityRatio;
    for (layers) |layer|
        if (!finiteLayer(layer)) return error.InvalidSoilWaterLayer;
    if (forcing_zone_index != reference_forcing_zone_index) return;

    var next = balance.*;
    for (layers) |layer| {
        const water_equivalent_m3 = layer.micropore_water_m3 + layer.water_vapour_m3 +
            layer.macropore_water_m3 +
            (layer.micropore_ice_m3 + layer.macropore_ice_m3) * ice_water_density_ratio;
        next.landscape_water_m3 = next.landscape_water_m3 + water_equivalent_m3;
        next.landscape_ice_m3 = next.landscape_ice_m3 + layer.micropore_ice_m3 + layer.macropore_ice_m3;
        next.grid_cell_water_m3 = next.grid_cell_water_m3 + water_equivalent_m3;
        inline for (@typeInfo(WaterBalance).@"struct".fields) |field|
            if (!std.math.isFinite(@field(next, field.name)))
                return error.NonFiniteSoilWaterBalance;
    }
    balance.* = next;
}

test "REDIST soil water balance preserves source expression order" {
    var balance = std.mem.zeroes(WaterBalance);
    const layers = [_]LayerWaterState{.{
        .micropore_water_m3 = 1.0,
        .water_vapour_m3 = 2.0,
        .macropore_water_m3 = 3.0,
        .micropore_ice_m3 = 4.0,
        .macropore_ice_m3 = 5.0,
    }};
    try accumulateLayers(&balance, &layers, 2, 2, 0.9);
    try std.testing.expectEqual(@as(f64, 14.1), balance.landscape_water_m3);
    try std.testing.expectEqual(@as(f64, 9.0), balance.landscape_ice_m3);
    try std.testing.expectEqual(@as(f64, 14.1), balance.grid_cell_water_m3);
}

test "REDIST soil water balance gate leaves non-reference zone unchanged" {
    var balance = WaterBalance{
        .landscape_water_m3 = 1.0,
        .landscape_ice_m3 = 2.0,
        .grid_cell_water_m3 = 3.0,
    };
    const layers = [_]LayerWaterState{std.mem.zeroes(LayerWaterState)};
    try accumulateLayers(&balance, &layers, 1, 2, 0.9);
    try std.testing.expectEqual(@as(f64, 1.0), balance.landscape_water_m3);
    try std.testing.expectEqual(@as(f64, 2.0), balance.landscape_ice_m3);
    try std.testing.expectEqual(@as(f64, 3.0), balance.grid_cell_water_m3);
}

test "REDIST soil water balance runtime layers retain ascending order" {
    var balance = WaterBalance{
        .landscape_water_m3 = 1.0,
        .landscape_ice_m3 = 0.0,
        .grid_cell_water_m3 = 1.0,
    };
    const layers = [_]LayerWaterState{
        .{ .micropore_water_m3 = 1.0e16, .water_vapour_m3 = 0, .macropore_water_m3 = 0, .micropore_ice_m3 = 0, .macropore_ice_m3 = 0 },
        .{ .micropore_water_m3 = -1.0e16, .water_vapour_m3 = 0, .macropore_water_m3 = 0, .micropore_ice_m3 = 0, .macropore_ice_m3 = 0 },
    };
    try accumulateLayers(&balance, &layers, 0, 0, 0.9);
    try std.testing.expectEqual(@as(f64, 0.0), balance.landscape_water_m3);
    try std.testing.expectEqual(@as(f64, 0.0), balance.grid_cell_water_m3);
}

test "REDIST soil water balance rejects dimensions invalid input and overflow" {
    var balance = std.mem.zeroes(WaterBalance);
    const no_layers: [0]LayerWaterState = .{};
    try std.testing.expectError(error.SoilWaterBalanceDimensionMismatch, accumulateLayers(&balance, &no_layers, 0, 0, 0.9));
    const layers = [_]LayerWaterState{std.mem.zeroes(LayerWaterState)};
    try std.testing.expectError(error.InvalidIceWaterDensityRatio, accumulateLayers(&balance, &layers, 0, 0, 0.0));
    var bad_layers = layers;
    bad_layers[0].water_vapour_m3 = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilWaterLayer, accumulateLayers(&balance, &bad_layers, 0, 0, 0.9));
    balance.landscape_water_m3 = std.math.floatMax(f64);
    var overflow_layers = layers;
    overflow_layers[0].micropore_water_m3 = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSoilWaterBalance, accumulateLayers(&balance, &overflow_layers, 0, 0, 0.9));
}
