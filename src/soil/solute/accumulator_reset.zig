const std = @import("std");

/// Micropores retain H4SiO4: 34 salt/complex plus 16 phosphate fields.
pub const micropore_species_per_layer = 50;
/// Macropores omit H4SiO4: 33 salt/complex plus 16 phosphate fields.
pub const macropore_species_per_layer = 49;

pub const ActiveLayerRange = struct {
    /// Zero-based translation of legacy `NU`; inclusive.
    first_layer_index: usize,
    /// Zero-based translation of legacy `NL`; exclusive in Zig.
    end_layer_index: usize,
};

pub const SoilAccumulators = struct {
    total_layer_count: usize,
    active_layers: ActiveLayerRange,
    /// Layer-major micropore net solute fluxes, mol layer-1 step-1.
    micropore_flux_mol_per_layer_step: []f64,
    /// Layer-major macropore net solute fluxes, mol layer-1 step-1.
    macropore_flux_mol_per_layer_step: []f64,
};

/// Exact source-order translation of TRNSFRS.F lines 1758--1857.
///
/// The two compact topologies deliberately differ by the micropore H4SiO4
/// slot. Validation completes before the active layer range is mutated.
pub fn reset(accumulators: SoilAccumulators) !void {
    const range = accumulators.active_layers;
    if (range.first_layer_index > range.end_layer_index or
        range.end_layer_index > accumulators.total_layer_count)
        return error.InvalidActiveSoilLayerRange;

    const micropore_len = std.math.mul(
        usize,
        accumulators.total_layer_count,
        micropore_species_per_layer,
    ) catch return error.SoilSoluteAccumulatorDimensionOverflow;
    const macropore_len = std.math.mul(
        usize,
        accumulators.total_layer_count,
        macropore_species_per_layer,
    ) catch return error.SoilSoluteAccumulatorDimensionOverflow;
    if (accumulators.micropore_flux_mol_per_layer_step.len != micropore_len)
        return error.MicroporeSoluteAccumulatorDimensionMismatch;
    if (accumulators.macropore_flux_mol_per_layer_step.len != macropore_len)
        return error.MacroporeSoluteAccumulatorDimensionMismatch;

    const micropore_start = range.first_layer_index * micropore_species_per_layer;
    const micropore_end = range.end_layer_index * micropore_species_per_layer;
    const macropore_start = range.first_layer_index * macropore_species_per_layer;
    const macropore_end = range.end_layer_index * macropore_species_per_layer;
    @memset(accumulators.micropore_flux_mol_per_layer_step[micropore_start..micropore_end], 0);
    @memset(accumulators.macropore_flux_mol_per_layer_step[macropore_start..macropore_end], 0);
}

test "TRNSFRS resets only the inclusive NU through NL soil layers" {
    const layer_count = 4;
    var micropore = [_]f64{3} ** (layer_count * micropore_species_per_layer);
    var macropore = [_]f64{5} ** (layer_count * macropore_species_per_layer);
    try reset(.{
        .total_layer_count = layer_count,
        .active_layers = .{ .first_layer_index = 1, .end_layer_index = 3 },
        .micropore_flux_mol_per_layer_step = &micropore,
        .macropore_flux_mol_per_layer_step = &macropore,
    });

    try std.testing.expectEqual(@as(f64, 3), micropore[0]);
    try std.testing.expectEqual(@as(f64, 0), micropore[micropore_species_per_layer]);
    try std.testing.expectEqual(@as(f64, 0), micropore[3 * micropore_species_per_layer - 1]);
    try std.testing.expectEqual(@as(f64, 3), micropore[3 * micropore_species_per_layer]);
    try std.testing.expectEqual(@as(f64, 5), macropore[0]);
    try std.testing.expectEqual(@as(f64, 0), macropore[macropore_species_per_layer]);
    try std.testing.expectEqual(@as(f64, 5), macropore[3 * macropore_species_per_layer]);
}

test "TRNSFRS soil reset preserves an empty runtime layer range" {
    var micropore = [_]f64{3} ** micropore_species_per_layer;
    var macropore = [_]f64{5} ** macropore_species_per_layer;
    try reset(.{
        .total_layer_count = 1,
        .active_layers = .{ .first_layer_index = 1, .end_layer_index = 1 },
        .micropore_flux_mol_per_layer_step = &micropore,
        .macropore_flux_mol_per_layer_step = &macropore,
    });
    try std.testing.expectEqual(@as(f64, 3), micropore[0]);
    try std.testing.expectEqual(@as(f64, 5), macropore[0]);
}

test "soil reset preserves exact differing compact topologies" {
    try std.testing.expectEqual(@as(usize, 50), micropore_species_per_layer);
    try std.testing.expectEqual(@as(usize, 49), macropore_species_per_layer);
}

test "late macropore dimension failure leaves micropores unchanged" {
    var micropore = [_]f64{3} ** micropore_species_per_layer;
    var short_macropore = [_]f64{5} ** (macropore_species_per_layer - 1);
    try std.testing.expectError(
        error.MacroporeSoluteAccumulatorDimensionMismatch,
        reset(.{
            .total_layer_count = 1,
            .active_layers = .{ .first_layer_index = 0, .end_layer_index = 1 },
            .micropore_flux_mol_per_layer_step = &micropore,
            .macropore_flux_mol_per_layer_step = &short_macropore,
        }),
    );
    try std.testing.expectEqualSlices(f64, &([_]f64{3} ** micropore_species_per_layer), &micropore);
}

test "invalid active range fails before mutation" {
    var micropore = [_]f64{3} ** micropore_species_per_layer;
    var macropore = [_]f64{5} ** macropore_species_per_layer;
    try std.testing.expectError(
        error.InvalidActiveSoilLayerRange,
        reset(.{
            .total_layer_count = 1,
            .active_layers = .{ .first_layer_index = 1, .end_layer_index = 0 },
            .micropore_flux_mol_per_layer_step = &micropore,
            .macropore_flux_mol_per_layer_step = &macropore,
        }),
    );
    try std.testing.expectEqual(@as(f64, 3), micropore[0]);
    try std.testing.expectEqual(@as(f64, 5), macropore[0]);
}
