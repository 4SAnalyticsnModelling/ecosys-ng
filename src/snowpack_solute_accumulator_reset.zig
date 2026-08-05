const std = @import("std");

/// Snowpack layers store 33 salt/complex species followed by eight phosphate
/// species. H4SiO4 is absent from this compact source topology.
pub const species_per_snow_layer = 41;

pub const SnowpackAccumulators = struct {
    /// Runtime snow-layer count (`JS`). A zero count preserves Fortran's
    /// zero-trip `DO L=1,JS` behavior.
    snow_layer_count: usize,
    /// Layer-major net solute fluxes, mol layer-1 step-1.
    flux_mol_per_layer_step: []f64,
};

/// Exact source-order translation of TRNSFRS.F lines 1712--1754.
///
/// The caller owns the enclosing sub-hour, cell, and dynamic-salt loops.
pub fn reset(accumulators: SnowpackAccumulators) !void {
    const required_len = std.math.mul(
        usize,
        accumulators.snow_layer_count,
        species_per_snow_layer,
    ) catch return error.SnowpackSoluteAccumulatorDimensionOverflow;
    if (accumulators.flux_mol_per_layer_step.len != required_len)
        return error.SnowpackSoluteAccumulatorDimensionMismatch;

    @memset(accumulators.flux_mol_per_layer_step, 0);
}

test "TRNSFRS resets every compact snowpack species in every runtime layer" {
    const layer_count = 3;
    var fluxes = [_]f64{7} ** (layer_count * species_per_snow_layer);

    try reset(.{
        .snow_layer_count = layer_count,
        .flux_mol_per_layer_step = &fluxes,
    });

    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{0} ** (layer_count * species_per_snow_layer)),
        &fluxes,
    );
}

test "TRNSFRS snowpack reset preserves JS zero-trip behavior" {
    var no_fluxes: [0]f64 = .{};
    try reset(.{
        .snow_layer_count = 0,
        .flux_mol_per_layer_step = &no_fluxes,
    });
}

test "snowpack topology omits the H4SiO4 phantom slot" {
    try std.testing.expectEqual(@as(usize, 41), species_per_snow_layer);
}

test "snowpack dimension failure is atomic" {
    const layer_count = 2;
    var short_fluxes = [_]f64{9} ** (layer_count * species_per_snow_layer - 1);
    try std.testing.expectError(
        error.SnowpackSoluteAccumulatorDimensionMismatch,
        reset(.{
            .snow_layer_count = layer_count,
            .flux_mol_per_layer_step = &short_fluxes,
        }),
    );
    try std.testing.expectEqualSlices(
        f64,
        &([_]f64{9} ** (layer_count * species_per_snow_layer - 1)),
        &short_fluxes,
    );
}

test "snowpack dimension multiplication overflow is explicit" {
    var no_fluxes: [0]f64 = .{};
    try std.testing.expectError(
        error.SnowpackSoluteAccumulatorDimensionOverflow,
        reset(.{
            .snow_layer_count = std.math.maxInt(usize),
            .flux_mol_per_layer_step = &no_fluxes,
        }),
    );
}
