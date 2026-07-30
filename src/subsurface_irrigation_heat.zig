const std = @import("std");

pub const Totals = struct {
    water_input_m3: f64,
    heat_input_mj: f64,
};

/// Calculates the accepted subsurface-irrigation water and sensible-heat
/// boundary input. Each runtime soil layer uses the atmospheric temperature
/// of its owning horizontal grid cell, matching the temperature used when
/// the same water is injected into the implicit soil heat residual.
pub fn calculate(
    water_input_m3_by_layer: []const f64,
    atmospheric_temperature_k_by_cell: []const f64,
    soil_layer_capacity: usize,
    liquid_water_heat_capacity_mj_per_m3_k: f64,
) !Totals {
    try validateDimensions(
        water_input_m3_by_layer,
        atmospheric_temperature_k_by_cell,
        soil_layer_capacity,
    );
    if (!std.math.isFinite(liquid_water_heat_capacity_mj_per_m3_k) or
        liquid_water_heat_capacity_mj_per_m3_k <= 0)
        return error.InvalidLiquidWaterHeatCapacity;

    var totals: Totals = .{
        .water_input_m3 = 0,
        .heat_input_mj = 0,
    };
    for (water_input_m3_by_layer, 0..) |water_input_m3, layer| {
        if (!std.math.isFinite(water_input_m3) or water_input_m3 < 0)
            return error.InvalidSubsurfaceIrrigationWater;
        const cell = layer / soil_layer_capacity;
        const atmospheric_temperature_k =
            atmospheric_temperature_k_by_cell[cell];
        if (!std.math.isFinite(atmospheric_temperature_k) or
            atmospheric_temperature_k <= 0)
            return error.InvalidAtmosphericTemperature;
        totals.water_input_m3 += water_input_m3;
        totals.heat_input_mj += liquid_water_heat_capacity_mj_per_m3_k *
            atmospheric_temperature_k * water_input_m3;
        if (!std.math.isFinite(totals.water_input_m3) or
            !std.math.isFinite(totals.heat_input_mj))
            return error.NonFiniteSubsurfaceIrrigationHeat;
    }
    return totals;
}

/// Adds the same layer-resolved sensible heat used by `calculate` to the
/// implicit soil heat source. Validation and accumulation occur before
/// mutation, so invalid input cannot partially modify the destination.
pub fn addToLayerHeatSources(
    heat_source_mj_by_layer: []f64,
    water_input_m3_by_layer: []const f64,
    atmospheric_temperature_k_by_cell: []const f64,
    soil_layer_capacity: usize,
    liquid_water_heat_capacity_mj_per_m3_k: f64,
) !Totals {
    if (heat_source_mj_by_layer.len != water_input_m3_by_layer.len)
        return error.SubsurfaceIrrigationLayerDimensionMismatch;
    for (heat_source_mj_by_layer) |heat_source_mj| {
        if (!std.math.isFinite(heat_source_mj))
            return error.NonFiniteSoilHeatSource;
    }
    const totals = try calculate(
        water_input_m3_by_layer,
        atmospheric_temperature_k_by_cell,
        soil_layer_capacity,
        liquid_water_heat_capacity_mj_per_m3_k,
    );
    for (
        heat_source_mj_by_layer,
        water_input_m3_by_layer,
        0..,
    ) |*heat_source_mj, water_input_m3, layer| {
        const cell = layer / soil_layer_capacity;
        heat_source_mj.* += liquid_water_heat_capacity_mj_per_m3_k *
            atmospheric_temperature_k_by_cell[cell] * water_input_m3;
    }
    return totals;
}

fn validateDimensions(
    water_input_m3_by_layer: []const f64,
    atmospheric_temperature_k_by_cell: []const f64,
    soil_layer_capacity: usize,
) !void {
    if (soil_layer_capacity == 0)
        return error.InvalidSoilLayerCapacity;
    const expected_layer_count = std.math.mul(
        usize,
        atmospheric_temperature_k_by_cell.len,
        soil_layer_capacity,
    ) catch return error.SubsurfaceIrrigationDimensionOverflow;
    if (water_input_m3_by_layer.len != expected_layer_count)
        return error.SubsurfaceIrrigationLayerDimensionMismatch;
}

test "runtime layers use their owning cell temperature and conserve totals" {
    const water_input_m3 = [_]f64{
        1, 2, 0,
        4, 0, 5,
    };
    const temperature_k = [_]f64{ 280, 300 };
    var heat_source_mj = [_]f64{ 10, 20, 30, 40, 50, 60 };
    const heat_capacity_mj_per_m3_k = 4.19;

    const totals = try addToLayerHeatSources(
        &heat_source_mj,
        &water_input_m3,
        &temperature_k,
        3,
        heat_capacity_mj_per_m3_k,
    );
    try std.testing.expectEqual(@as(f64, 12), totals.water_input_m3);
    try std.testing.expectApproxEqAbs(
        4.19 * (280 * 3 + 300 * 9),
        totals.heat_input_mj,
        1e-10,
    );
    try std.testing.expectApproxEqAbs(
        10 + 4.19 * 280,
        heat_source_mj[0],
        1e-12,
    );
    try std.testing.expectEqual(@as(f64, 30), heat_source_mj[2]);
    try std.testing.expectApproxEqAbs(
        60 + 4.19 * 300 * 5,
        heat_source_mj[5],
        1e-12,
    );
}

test "invalid irrigation input leaves every heat source unchanged" {
    var heat_source_mj = [_]f64{ 1, 2 };
    const original = heat_source_mj;
    const water_input_m3 = [_]f64{ 0.5, -0.25 };
    const temperature_k = [_]f64{290};

    try std.testing.expectError(
        error.InvalidSubsurfaceIrrigationWater,
        addToLayerHeatSources(
            &heat_source_mj,
            &water_input_m3,
            &temperature_k,
            2,
            4.19,
        ),
    );
    try std.testing.expectEqualSlices(f64, &original, &heat_source_mj);
}
