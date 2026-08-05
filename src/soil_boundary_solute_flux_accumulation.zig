const std = @import("std");

pub const micropore_species_count = 50;
pub const macropore_species_count = 49;

pub const Accumulators = struct {
    micropore_amount_per_step: []f64,
    macropore_amount_per_step: []f64,
};

pub const Fluxes = struct {
    micropore_amount_per_step: []const f64,
    macropore_amount_per_step: []const f64,
};

/// Compatibility translation of TRNSFRS.F lines 6179--6277.
/// The legacy source accumulates all 50 micropore species first, including
/// H4SiO4, and then the compact 49-species macropore topology.
pub fn accumulate(accumulators: Accumulators, fluxes: Fluxes) !void {
    try validateDimensions(accumulators, fluxes);
    for (accumulators.micropore_amount_per_step, fluxes.micropore_amount_per_step) |total, flux| {
        if (!std.math.isFinite(total) or !std.math.isFinite(flux))
            return error.NonFiniteSoilBoundarySoluteFluxAccumulationInput;
        if (!std.math.isFinite(total + flux))
            return error.NonFiniteSoilBoundarySoluteFluxAccumulationResult;
    }
    for (accumulators.macropore_amount_per_step, fluxes.macropore_amount_per_step) |total, flux| {
        if (!std.math.isFinite(total) or !std.math.isFinite(flux))
            return error.NonFiniteSoilBoundarySoluteFluxAccumulationInput;
        if (!std.math.isFinite(total + flux))
            return error.NonFiniteSoilBoundarySoluteFluxAccumulationResult;
    }
    for (accumulators.micropore_amount_per_step, fluxes.micropore_amount_per_step) |*total, flux|
        total.* = total.* + flux;
    for (accumulators.macropore_amount_per_step, fluxes.macropore_amount_per_step) |*total, flux|
        total.* = total.* + flux;
}

fn validateDimensions(accumulators: Accumulators, fluxes: Fluxes) !void {
    if (accumulators.micropore_amount_per_step.len != micropore_species_count or
        fluxes.micropore_amount_per_step.len != micropore_species_count or
        accumulators.macropore_amount_per_step.len != macropore_species_count or
        fluxes.macropore_amount_per_step.len != macropore_species_count)
        return error.SoilBoundarySoluteFluxAccumulationDimensionMismatch;
}

test "TRNSFRS accumulates every micropore and compact macropore species" {
    var micropore_total = [_]f64{1} ** micropore_species_count;
    var macropore_total = [_]f64{2} ** macropore_species_count;
    const micropore_flux = [_]f64{0.5} ** micropore_species_count;
    const macropore_flux = [_]f64{-0.25} ** macropore_species_count;
    try accumulate(.{ .micropore_amount_per_step = &micropore_total, .macropore_amount_per_step = &macropore_total }, .{ .micropore_amount_per_step = &micropore_flux, .macropore_amount_per_step = &macropore_flux });
    try std.testing.expectEqual(@as(f64, 1.5), micropore_total[0]);
    try std.testing.expectEqual(@as(f64, 1.5), micropore_total[49]);
    try std.testing.expectEqual(@as(f64, 1.75), macropore_total[0]);
    try std.testing.expectEqual(@as(f64, 1.75), macropore_total[48]);
}

test "TRNSFRS preserves micropore H4SiO4 slot without macropore counterpart" {
    var micropore_total = [_]f64{0} ** micropore_species_count;
    var macropore_total = [_]f64{0} ** macropore_species_count;
    var micropore_flux = [_]f64{0} ** micropore_species_count;
    const macropore_flux = [_]f64{0} ** macropore_species_count;
    micropore_flux[33] = 4;
    try accumulate(.{ .micropore_amount_per_step = &micropore_total, .macropore_amount_per_step = &macropore_total }, .{ .micropore_amount_per_step = &micropore_flux, .macropore_amount_per_step = &macropore_flux });
    try std.testing.expectEqual(@as(f64, 4), micropore_total[33]);
    try std.testing.expectEqual(@as(usize, 49), macropore_total.len);
}

test "TRNSFRS repeated substep accumulation retains source addition order" {
    var micropore_total = [_]f64{1.0e16} ** micropore_species_count;
    var macropore_total = [_]f64{1.0e16} ** macropore_species_count;
    const first_micro = [_]f64{-1.0e16} ** micropore_species_count;
    const first_macro = [_]f64{-1.0e16} ** macropore_species_count;
    const second_micro = [_]f64{1} ** micropore_species_count;
    const second_macro = [_]f64{1} ** macropore_species_count;
    try accumulate(.{ .micropore_amount_per_step = &micropore_total, .macropore_amount_per_step = &macropore_total }, .{ .micropore_amount_per_step = &first_micro, .macropore_amount_per_step = &first_macro });
    try accumulate(.{ .micropore_amount_per_step = &micropore_total, .macropore_amount_per_step = &macropore_total }, .{ .micropore_amount_per_step = &second_micro, .macropore_amount_per_step = &second_macro });
    try std.testing.expectEqual(@as(f64, 1), micropore_total[0]);
    try std.testing.expectEqual(@as(f64, 1), macropore_total[0]);
}

test "dimension mismatch leaves both accumulator families unchanged" {
    var micropore_total = [_]f64{3} ** micropore_species_count;
    var macropore_total = [_]f64{4} ** macropore_species_count;
    const short_micro = [_]f64{1} ** (micropore_species_count - 1);
    const macropore_flux = [_]f64{1} ** macropore_species_count;
    try std.testing.expectError(error.SoilBoundarySoluteFluxAccumulationDimensionMismatch, accumulate(.{ .micropore_amount_per_step = &micropore_total, .macropore_amount_per_step = &macropore_total }, .{ .micropore_amount_per_step = &short_micro, .macropore_amount_per_step = &macropore_flux }));
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** micropore_species_count), &micropore_total);
    try std.testing.expectEqualSlices(f64, &([_]f64{4} ** macropore_species_count), &macropore_total);
}

test "late invalid macropore flux leaves micropore accumulation atomic" {
    var micropore_total = [_]f64{3} ** micropore_species_count;
    var macropore_total = [_]f64{4} ** macropore_species_count;
    const micropore_flux = [_]f64{1} ** micropore_species_count;
    var macropore_flux = [_]f64{1} ** macropore_species_count;
    macropore_flux[48] = std.math.inf(f64);
    try std.testing.expectError(error.NonFiniteSoilBoundarySoluteFluxAccumulationInput, accumulate(.{ .micropore_amount_per_step = &micropore_total, .macropore_amount_per_step = &macropore_total }, .{ .micropore_amount_per_step = &micropore_flux, .macropore_amount_per_step = &macropore_flux }));
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** micropore_species_count), &micropore_total);
    try std.testing.expectEqualSlices(f64, &([_]f64{4} ** macropore_species_count), &macropore_total);
}
