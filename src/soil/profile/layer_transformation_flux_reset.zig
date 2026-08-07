const std = @import("std");

/// Source-ordered transformation fluxes from `hour1.f` lines 3342--3386.
/// Values retain the legacy mass-per-timestep units for their named species.
pub const TransformationFlux = enum {
    ammonium_non_band,
    ammonia_non_band,
    hydrogen_ion,
    ammonia_gas,
    nitrate,
    nitrite,
    exchangeable_ammonium,
    dihydrogen_phosphate,
    hydrogen_phosphate,
    aluminum_hydroxide_0,
    aluminum_hydroxide_1,
    aluminum_hydroxide_2,
    phosphate_1,
    phosphate_2,
    aluminum_phosphate,
    iron_phosphate,
    calcium_phosphate_dicalcium,
    calcium_phosphate_hydroxy,
    calcium_phosphate_mono,
    aluminum,
    iron,
    calcium,
    magnesium,
    sodium,
    potassium,
    bicarbonate_exchange,
    hydrogen,
    hydroxide,
    water,
    band_water,
    carbonate,
    bicarbonate,
    carbon_dioxide,
    aluminum_hydroxide,
    iron_hydroxide,
    calcium_carbonate,
    calcium_sulfate,
    exchangeable_hydrogen,
    exchangeable_aluminum,
    exchangeable_iron,
    exchangeable_calcium,
    exchangeable_magnesium,
    exchangeable_sodium,
    exchangeable_potassium,
    sulfate,
};

pub const LayerFluxes = struct {
    transformations: []f64,
};

pub const ResetError = error{FluxCountMismatch};

/// Translates `hour1.f` lines 3342--3386 over runtime-allocated soil layers.
pub fn reset(layers: []LayerFluxes) ResetError!void {
    const flux_count = std.meta.fields(TransformationFlux).len;
    for (layers) |layer| {
        if (layer.transformations.len != flux_count) return error.FluxCountMismatch;
    }
    for (layers) |layer| {
        for (layer.transformations) |*flux| flux.* = 0.0;
    }
}

test "transformation reset covers all named fluxes and runtime layers" {
    const allocator = std.testing.allocator;
    const layers = try allocator.alloc(LayerFluxes, 9);
    defer allocator.free(layers);
    for (layers) |*layer| {
        layer.transformations =
            try allocator.alloc(f64, std.meta.fields(TransformationFlux).len);
        @memset(layer.transformations, 6.0);
    }
    defer for (layers) |layer| allocator.free(layer.transformations);

    try reset(layers);

    for (layers) |layer| {
        for (layer.transformations) |flux| {
            try std.testing.expectEqual(@as(f64, 0.0), flux);
        }
    }
}

test "dimension mismatch fails before any layer mutation" {
    var valid = [_]f64{8.0} ** std.meta.fields(TransformationFlux).len;
    var short = [_]f64{9.0};
    var layers = [_]LayerFluxes{
        .{ .transformations = &valid },
        .{ .transformations = &short },
    };
    try std.testing.expectError(error.FluxCountMismatch, reset(&layers));
    try std.testing.expectEqual(@as(f64, 8.0), valid[0]);
}
