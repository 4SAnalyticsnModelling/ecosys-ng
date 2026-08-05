const std = @import("std");

/// Soil phosphate precipitate pools (mol). 10 fields, REDIST 6378-6387.
pub const PhosphatePrecipitatePools = struct {
    alpo4_nonband_mol: f64, // PALPO
    fepo4_nonband_mol: f64, // PFEPO
    cahpo4_nonband_mol: f64, // PCAPD
    apatite_nonband_mol: f64, // PCAPH
    cah2po4_nonband_mol: f64, // PCAPM
    alpo4_band_mol: f64, // PALPB
    fepo4_band_mol: f64, // PFEPB
    cahpo4_band_mol: f64, // PCPDB
    apatite_band_mol: f64, // PCPHB
    cah2po4_band_mol: f64, // PCPMB
};

/// Precipitation-dissolution transformation fluxes (mol step-1). Same 10 fields.
pub const PhosphatePrecipitateFluxes = struct {
    alpo4_nonband_mol: f64, // TRALPO
    fepo4_nonband_mol: f64, // TRFEPO
    cahpo4_nonband_mol: f64, // TRCAPD
    apatite_nonband_mol: f64, // TRCAPH
    cah2po4_nonband_mol: f64, // TRCAPM
    alpo4_band_mol: f64, // TRALPB
    fepo4_band_mol: f64, // TRFEPB
    cahpo4_band_mol: f64, // TRCPDB
    apatite_band_mol: f64, // TRCPHB
    cah2po4_band_mol: f64, // TRCPMB
};

/// Direct translation of REDIST lines 6378--6387 (inner body of DO 125 L loop).
pub fn update(pools: PhosphatePrecipitatePools, fluxes: PhosphatePrecipitateFluxes) !PhosphatePrecipitatePools {
    inline for (@typeInfo(PhosphatePrecipitatePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(pools, field.name)))
            return error.InvalidPhosphatePrecipPool;
    inline for (@typeInfo(PhosphatePrecipitateFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(fluxes, field.name)))
            return error.InvalidPhosphatePrecipFlux;

    const result = PhosphatePrecipitatePools{
        .alpo4_nonband_mol = pools.alpo4_nonband_mol + fluxes.alpo4_nonband_mol,
        .fepo4_nonband_mol = pools.fepo4_nonband_mol + fluxes.fepo4_nonband_mol,
        .cahpo4_nonband_mol = pools.cahpo4_nonband_mol + fluxes.cahpo4_nonband_mol,
        .apatite_nonband_mol = pools.apatite_nonband_mol + fluxes.apatite_nonband_mol,
        .cah2po4_nonband_mol = pools.cah2po4_nonband_mol + fluxes.cah2po4_nonband_mol,
        .alpo4_band_mol = pools.alpo4_band_mol + fluxes.alpo4_band_mol,
        .fepo4_band_mol = pools.fepo4_band_mol + fluxes.fepo4_band_mol,
        .cahpo4_band_mol = pools.cahpo4_band_mol + fluxes.cahpo4_band_mol,
        .apatite_band_mol = pools.apatite_band_mol + fluxes.apatite_band_mol,
        .cah2po4_band_mol = pools.cah2po4_band_mol + fluxes.cah2po4_band_mol,
    };
    inline for (@typeInfo(PhosphatePrecipitatePools).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFinitePhosphatePrecipPool;
    return result;
}

/// Applies precipitation-dissolution updates to all runtime soil layers.
pub fn updateLayers(
    pools_by_layer: []PhosphatePrecipitatePools,
    fluxes_by_layer: []const PhosphatePrecipitateFluxes,
) !void {
    if (pools_by_layer.len == 0 or fluxes_by_layer.len != pools_by_layer.len)
        return error.PhosphatePrecipitateDimensionMismatch;

    for (pools_by_layer, fluxes_by_layer) |*pools, fluxes|
        pools.* = try update(pools.*, fluxes);
}

test "REDIST soil phosphate precipitate non-band and band are independent" {
    var f = std.mem.zeroes(PhosphatePrecipitateFluxes);
    f.alpo4_nonband_mol = 1.0;
    f.alpo4_band_mol = 2.0;
    const result = try update(std.mem.zeroes(PhosphatePrecipitatePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 1.0), result.alpo4_nonband_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 2.0), result.alpo4_band_mol, 1.0e-15);
}

test "REDIST soil phosphate apatite and CaH2PO4 update independently" {
    var f = std.mem.zeroes(PhosphatePrecipitateFluxes);
    f.apatite_nonband_mol = 0.5;
    f.cah2po4_nonband_mol = 0.3;
    const result = try update(std.mem.zeroes(PhosphatePrecipitatePools), f);
    try std.testing.expectApproxEqRel(@as(f64, 0.5), result.apatite_nonband_mol, 1.0e-15);
    try std.testing.expectApproxEqRel(@as(f64, 0.3), result.cah2po4_nonband_mol, 1.0e-15);
}

test "REDIST soil phosphate precipitate rejects non-finite flux" {
    var bad = std.mem.zeroes(PhosphatePrecipitateFluxes);
    bad.fepo4_band_mol = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidPhosphatePrecipFlux,
        update(std.mem.zeroes(PhosphatePrecipitatePools), bad),
    );
}

test "REDIST soil phosphate precipitate maps all source transformations" {
    var fluxes: PhosphatePrecipitateFluxes = undefined;
    inline for (@typeInfo(PhosphatePrecipitateFluxes).@"struct".fields, 1..) |field, value|
        @field(fluxes, field.name) = @floatFromInt(value);

    const result = try update(std.mem.zeroes(PhosphatePrecipitatePools), fluxes);
    inline for (@typeInfo(PhosphatePrecipitatePools).@"struct".fields, 1..) |field, expected|
        try std.testing.expectEqual(@as(f64, @floatFromInt(expected)), @field(result, field.name));
}

test "REDIST soil phosphate precipitate runtime layers are independent" {
    var pools = [_]PhosphatePrecipitatePools{
        std.mem.zeroes(PhosphatePrecipitatePools),
        std.mem.zeroes(PhosphatePrecipitatePools),
    };
    var fluxes = [_]PhosphatePrecipitateFluxes{
        std.mem.zeroes(PhosphatePrecipitateFluxes),
        std.mem.zeroes(PhosphatePrecipitateFluxes),
    };
    fluxes[0].cahpo4_nonband_mol = 2.5;
    fluxes[1].apatite_band_mol = 4.5;

    try updateLayers(&pools, &fluxes);

    try std.testing.expectEqual(@as(f64, 2.5), pools[0].cahpo4_nonband_mol);
    try std.testing.expectEqual(@as(f64, 0.0), pools[1].cahpo4_nonband_mol);
    try std.testing.expectEqual(@as(f64, 0.0), pools[0].apatite_band_mol);
    try std.testing.expectEqual(@as(f64, 4.5), pools[1].apatite_band_mol);
}

test "REDIST soil phosphate precipitate runtime layers reject dimensions and overflow" {
    var pools = [_]PhosphatePrecipitatePools{std.mem.zeroes(PhosphatePrecipitatePools)};
    const no_fluxes: [0]PhosphatePrecipitateFluxes = .{};
    try std.testing.expectError(
        error.PhosphatePrecipitateDimensionMismatch,
        updateLayers(&pools, &no_fluxes),
    );

    pools[0].cah2po4_nonband_mol = std.math.floatMax(f64);
    var fluxes = std.mem.zeroes(PhosphatePrecipitateFluxes);
    fluxes.cah2po4_nonband_mol = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFinitePhosphatePrecipPool, update(pools[0], fluxes));
}
