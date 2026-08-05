const std = @import("std");

pub const species_count = 50;

/// Compatibility translation of TRNSFRS.F lines 7956--8005.
/// Accumulates all external micropore boundary species in source order:
/// 34 salts/complexes including H4SiO4, 8 non-band P, then 8 band P.
pub fn accumulate(accumulated_amount_per_step: []f64, boundary_flux_amount_per_step: []const f64) !void {
    if (accumulated_amount_per_step.len != species_count or boundary_flux_amount_per_step.len != species_count)
        return error.ExternalMicroporeSoluteFluxAccumulationDimensionMismatch;
    for (accumulated_amount_per_step, boundary_flux_amount_per_step) |accumulated, flux| {
        if (!std.math.isFinite(accumulated) or !std.math.isFinite(flux))
            return error.NonFiniteExternalMicroporeSoluteFluxAccumulationInput;
        if (!std.math.isFinite(accumulated + flux))
            return error.NonFiniteExternalMicroporeSoluteFluxAccumulationResult;
    }
    for (accumulated_amount_per_step, boundary_flux_amount_per_step) |*accumulated, flux|
        accumulated.* = accumulated.* + flux;
}

test "TRNSFRS accumulates all 50 external micropore species" {
    var accumulated = [_]f64{2} ** species_count;
    const flux = [_]f64{0.5} ** species_count;
    try accumulate(&accumulated, &flux);
    try std.testing.expectEqual(@as(f64, 2.5), accumulated[0]);
    try std.testing.expectEqual(@as(f64, 2.5), accumulated[33]);
    try std.testing.expectEqual(@as(f64, 2.5), accumulated[34]);
    try std.testing.expectEqual(@as(f64, 2.5), accumulated[49]);
}

test "TRNSFRS accumulation retains signed recharge flux" {
    var accumulated = [_]f64{2} ** species_count;
    const flux = [_]f64{-3} ** species_count;
    try accumulate(&accumulated, &flux);
    try std.testing.expectEqualSlices(f64, &([_]f64{-1} ** species_count), &accumulated);
}

test "TRNSFRS repeated substeps preserve left-to-right addition" {
    var accumulated = [_]f64{1.0e16} ** species_count;
    const first = [_]f64{-1.0e16} ** species_count;
    const second = [_]f64{1} ** species_count;
    try accumulate(&accumulated, &first);
    try accumulate(&accumulated, &second);
    try std.testing.expectEqual(@as(f64, 1), accumulated[0]);
    try std.testing.expectEqual(@as(f64, 1), accumulated[49]);
}

test "dimension mismatch leaves every accumulator unchanged" {
    var accumulated = [_]f64{3} ** species_count;
    const short_flux = [_]f64{1} ** (species_count - 1);
    try std.testing.expectError(error.ExternalMicroporeSoluteFluxAccumulationDimensionMismatch, accumulate(&accumulated, &short_flux));
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** species_count), &accumulated);
}

test "late invalid flux keeps accumulation atomic" {
    var accumulated = [_]f64{3} ** species_count;
    var flux = [_]f64{1} ** species_count;
    flux[49] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteExternalMicroporeSoluteFluxAccumulationInput, accumulate(&accumulated, &flux));
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** species_count), &accumulated);
}
