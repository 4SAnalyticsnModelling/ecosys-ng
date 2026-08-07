const std = @import("std");

/// 33 salt/complex fields without H4SiO4, then 16 phosphate fields.
pub const exchange_species_per_layer = 49;

pub const LayerExchange = struct {
    total_layer_count: usize,
    surface_layer_index: usize,
    /// Layer-major current pore-domain exchange amount per model step.
    current_flux_amount_per_layer_step: []const f64,
    /// Layer-major REDIST accumulated exchange amount per model step.
    accumulated_flux_amount_per_layer_step: []f64,
};

/// Exact source-order translation of TRNSFRS.F lines 3697--3794.
/// Runtime dimensions replace the fixed legacy layer bounds. The complete
/// selected-layer sum is validated before any accumulator is mutated.
pub fn accumulate(exchange: LayerExchange) !void {
    if (exchange.surface_layer_index >= exchange.total_layer_count)
        return error.SurfacePoreExchangeLayerOutOfBounds;
    const required_len = std.math.mul(
        usize,
        exchange.total_layer_count,
        exchange_species_per_layer,
    ) catch return error.SurfacePoreExchangeAccumulationDimensionOverflow;
    if (exchange.current_flux_amount_per_layer_step.len != required_len or
        exchange.accumulated_flux_amount_per_layer_step.len != required_len)
        return error.SurfacePoreExchangeAccumulationDimensionMismatch;

    const start = exchange.surface_layer_index * exchange_species_per_layer;
    const end = start + exchange_species_per_layer;
    for (start..end) |index| {
        const current = exchange.current_flux_amount_per_layer_step[index];
        const accumulated = exchange.accumulated_flux_amount_per_layer_step[index];
        if (!std.math.isFinite(current) or !std.math.isFinite(accumulated))
            return error.NonFiniteSurfacePoreExchangeAccumulationInput;
        if (!std.math.isFinite(accumulated + current))
            return error.NonFiniteSurfacePoreExchangeAccumulationResult;
    }
    for (start..end) |index|
        exchange.accumulated_flux_amount_per_layer_step[index] +=
            exchange.current_flux_amount_per_layer_step[index];
}

test "TRNSFRS accumulates all 49 source fields at the runtime surface layer" {
    const layers = 3;
    var current = [_]f64{0} ** (layers * exchange_species_per_layer);
    var accumulated = [_]f64{10} ** (layers * exchange_species_per_layer);
    for (current[exchange_species_per_layer .. 2 * exchange_species_per_layer], 0..) |*value, species|
        value.* = @floatFromInt(species + 1);
    try accumulate(.{
        .total_layer_count = layers,
        .surface_layer_index = 1,
        .current_flux_amount_per_layer_step = &current,
        .accumulated_flux_amount_per_layer_step = &accumulated,
    });
    try std.testing.expectEqual(@as(f64, 10), accumulated[0]);
    try std.testing.expectEqual(@as(f64, 11), accumulated[exchange_species_per_layer]);
    try std.testing.expectEqual(@as(f64, 59), accumulated[2 * exchange_species_per_layer - 1]);
    try std.testing.expectEqual(@as(f64, 10), accumulated[2 * exchange_species_per_layer]);
}

test "compact accumulation rejects a phantom H4SiO4 field per layer" {
    const current = [_]f64{1} ** (exchange_species_per_layer + 1);
    var accumulated = [_]f64{2} ** exchange_species_per_layer;
    try std.testing.expectError(
        error.SurfacePoreExchangeAccumulationDimensionMismatch,
        accumulate(.{ .total_layer_count = 1, .surface_layer_index = 0, .current_flux_amount_per_layer_step = &current, .accumulated_flux_amount_per_layer_step = &accumulated }),
    );
}

test "late selected-layer invalid value leaves accumulation atomic" {
    const layers = 2;
    var current = [_]f64{1} ** (layers * exchange_species_per_layer);
    current[2 * exchange_species_per_layer - 1] = std.math.inf(f64);
    var accumulated = [_]f64{10} ** (layers * exchange_species_per_layer);
    try std.testing.expectError(
        error.NonFiniteSurfacePoreExchangeAccumulationInput,
        accumulate(.{ .total_layer_count = layers, .surface_layer_index = 1, .current_flux_amount_per_layer_step = &current, .accumulated_flux_amount_per_layer_step = &accumulated }),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{10} ** (layers * exchange_species_per_layer)), &accumulated);
}

test "runtime surface layer bounds fail before mutation" {
    const current = [_]f64{1} ** exchange_species_per_layer;
    var accumulated = [_]f64{10} ** exchange_species_per_layer;
    try std.testing.expectError(
        error.SurfacePoreExchangeLayerOutOfBounds,
        accumulate(.{ .total_layer_count = 1, .surface_layer_index = 1, .current_flux_amount_per_layer_step = &current, .accumulated_flux_amount_per_layer_step = &accumulated }),
    );
    try std.testing.expectEqual(@as(f64, 10), accumulated[0]);
}
