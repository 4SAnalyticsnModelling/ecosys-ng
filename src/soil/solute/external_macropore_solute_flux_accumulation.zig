const std = @import("std");

pub const species_count = 49;

/// Compatibility translation of TRNSFRS.F lines 8423--8471.
/// Preserves compact macropore order: 33 salts/complexes without H4SiO4,
/// followed by 8 non-band and 8 band phosphorus species.
pub fn accumulate(accumulated_amount_per_step: []f64, boundary_flux_amount_per_step: []const f64) !void {
    if (accumulated_amount_per_step.len != species_count or boundary_flux_amount_per_step.len != species_count)
        return error.ExternalMacroporeSoluteFluxAccumulationDimensionMismatch;
    for (accumulated_amount_per_step, boundary_flux_amount_per_step) |accumulated, flux| {
        if (!std.math.isFinite(accumulated) or !std.math.isFinite(flux))
            return error.NonFiniteExternalMacroporeSoluteFluxAccumulationInput;
        if (!std.math.isFinite(accumulated + flux))
            return error.NonFiniteExternalMacroporeSoluteFluxAccumulationResult;
    }
    for (accumulated_amount_per_step, boundary_flux_amount_per_step) |*accumulated, flux|
        accumulated.* = accumulated.* + flux;
}

test "TRNSFRS accumulates all 49 external macropore species" {
    var accumulated = [_]f64{2} ** species_count;
    const flux = [_]f64{0.5} ** species_count;
    try accumulate(&accumulated, &flux);
    try std.testing.expectEqual(@as(f64, 2.5), accumulated[0]);
    try std.testing.expectEqual(@as(f64, 2.5), accumulated[32]);
    try std.testing.expectEqual(@as(f64, 2.5), accumulated[33]);
    try std.testing.expectEqual(@as(f64, 2.5), accumulated[48]);
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
    try std.testing.expectEqual(@as(f64, 1), accumulated[48]);
}

test "dimension mismatch leaves compact accumulator unchanged" {
    var accumulated = [_]f64{3} ** species_count;
    const short_flux = [_]f64{1} ** (species_count - 1);
    try std.testing.expectError(error.ExternalMacroporeSoluteFluxAccumulationDimensionMismatch, accumulate(&accumulated, &short_flux));
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** species_count), &accumulated);
}

test "late invalid flux keeps compact accumulation atomic" {
    var accumulated = [_]f64{3} ** species_count;
    var flux = [_]f64{1} ** species_count;
    flux[48] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteExternalMacroporeSoluteFluxAccumulationInput, accumulate(&accumulated, &flux));
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** species_count), &accumulated);
}
