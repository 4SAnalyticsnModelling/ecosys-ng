const std = @import("std");

pub const species_count = 49;

/// Compatibility translation of TRNSFRS.F lines 6821--6869.
/// Accumulates the compact same-layer macro/micropore exchange topology in
/// exact source order for later REDIST use.
pub fn accumulate(accumulated_amount_per_step: []f64, exchange_amount_per_step: []const f64) !void {
    if (accumulated_amount_per_step.len != species_count or exchange_amount_per_step.len != species_count)
        return error.SoilPoreSoluteExchangeAccumulationDimensionMismatch;

    for (accumulated_amount_per_step, exchange_amount_per_step) |accumulated, exchange| {
        if (!std.math.isFinite(accumulated) or !std.math.isFinite(exchange))
            return error.NonFiniteSoilPoreSoluteExchangeAccumulationInput;
        if (!std.math.isFinite(accumulated + exchange))
            return error.NonFiniteSoilPoreSoluteExchangeAccumulationResult;
    }
    for (accumulated_amount_per_step, exchange_amount_per_step) |*accumulated, exchange|
        accumulated.* = accumulated.* + exchange;
}

test "TRNSFRS accumulates all 49 compact exchange species" {
    var accumulated = [_]f64{2} ** species_count;
    const exchange = [_]f64{0.5} ** species_count;
    try accumulate(&accumulated, &exchange);
    try std.testing.expectEqual(@as(f64, 2.5), accumulated[0]);
    try std.testing.expectEqual(@as(f64, 2.5), accumulated[32]);
    try std.testing.expectEqual(@as(f64, 2.5), accumulated[33]);
    try std.testing.expectEqual(@as(f64, 2.5), accumulated[48]);
}

test "TRNSFRS accumulation retains signed exchange" {
    var accumulated = [_]f64{2} ** species_count;
    const exchange = [_]f64{-3} ** species_count;
    try accumulate(&accumulated, &exchange);
    try std.testing.expectEqualSlices(f64, &([_]f64{-1} ** species_count), &accumulated);
}

test "TRNSFRS repeated substeps preserve left-to-right accumulation" {
    var accumulated = [_]f64{1.0e16} ** species_count;
    const first = [_]f64{-1.0e16} ** species_count;
    const second = [_]f64{1} ** species_count;
    try accumulate(&accumulated, &first);
    try accumulate(&accumulated, &second);
    try std.testing.expectEqual(@as(f64, 1), accumulated[0]);
    try std.testing.expectEqual(@as(f64, 1), accumulated[48]);
}

test "dimension mismatch leaves the accumulator unchanged" {
    var accumulated = [_]f64{3} ** species_count;
    const short_exchange = [_]f64{1} ** (species_count - 1);
    try std.testing.expectError(error.SoilPoreSoluteExchangeAccumulationDimensionMismatch, accumulate(&accumulated, &short_exchange));
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** species_count), &accumulated);
}

test "late invalid exchange leaves every accumulator atomic" {
    var accumulated = [_]f64{3} ** species_count;
    var exchange = [_]f64{1} ** species_count;
    exchange[48] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteSoilPoreSoluteExchangeAccumulationInput, accumulate(&accumulated, &exchange));
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** species_count), &accumulated);
}
