const std = @import("std");

pub const LayerPosition = enum { surface, subsurface };

/// Source-ordered band and mineral transformations, HOUR1 3431-3486.
pub const MineralTransformation = enum {
    band_ammonium,
    band_ammonia,
    band_nitrate,
    band_nitrite,
    band_dihydrogen_phosphate,
    band_hydrogen_phosphate,
    aluminum_hydroxide_1,
    aluminum_hydroxide_2,
    aluminum_hydroxide_3,
    aluminum_hydroxide_4,
    aluminum_sulfate,
    iron_hydroxide_1,
    iron_hydroxide_2,
    iron_hydroxide_3,
    iron_hydroxide_4,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
    hydrogen_phosphate_0,
    hydrogen_phosphate_3,
    calcium_phosphate_0,
    iron_phosphate_1,
    iron_phosphate_2,
    calcium_phosphate_1,
    calcium_phosphate_2,
    magnesium_phosphate_1,
    band_hydrogen_phosphate_0,
    band_hydrogen_phosphate_3,
    band_iron_phosphate_1,
    band_iron_phosphate_2,
    band_calcium_phosphate_0,
    band_calcium_phosphate_1,
    band_calcium_phosphate_2,
    band_magnesium_phosphate_1,
    band_exchangeable_nitrogen,
    exchangeable_aluminum,
    exchangeable_iron,
    band_aluminum_hydroxide_0,
    band_aluminum_hydroxide_1,
    band_aluminum_hydroxide_2,
    band_phosphate_1,
    band_phosphate_2,
    band_aluminum_phosphate,
    band_iron_phosphate,
    band_calcium_phosphate_dicalcium,
    band_calcium_phosphate_hydroxy,
    band_calcium_phosphate_mono,
};

/// Source-ordered micropore/macropore exchange fluxes, HOUR1 3487-3553.
pub const ExchangeFlux = enum {
    carbon_dioxide,
    methane,
    oxygen,
    nitrogen,
    nitrous_oxide,
    hydrogen,
    ammonium,
    ammonia,
    nitrate,
    other_nitrogen,
    dihydrogen_phosphate,
    hydrogen_phosphate,
    band_ammonium,
    band_ammonia,
    band_nitrate,
    band_other_nitrogen,
    band_dihydrogen_phosphate,
    band_hydrogen_phosphate,
    aluminum,
    iron,
    hydrogen_ion,
    calcium,
    magnesium,
    sodium,
    potassium,
    hydroxide,
    sulfate,
    chloride,
    carbonate,
    bicarbonate,
    aluminum_hydroxide_1,
    aluminum_hydroxide_2,
    aluminum_hydroxide_3,
    aluminum_hydroxide_4,
    aluminum_sulfate,
    iron_hydroxide_1,
    iron_hydroxide_2,
    iron_hydroxide_3,
    iron_hydroxide_4,
    iron_sulfate,
    calcium_hydroxide,
    calcium_carbonate,
    calcium_bicarbonate,
    calcium_sulfate,
    magnesium_hydroxide,
    magnesium_carbonate,
    magnesium_bicarbonate,
    magnesium_sulfate,
    sodium_carbonate,
    sodium_sulfate,
    potassium_sulfate,
    hydrogen_phosphate_0,
    hydrogen_phosphate_3,
    iron_phosphate_1,
    iron_phosphate_2,
    calcium_phosphate_0,
    calcium_phosphate_1,
    calcium_phosphate_2,
    magnesium_phosphate_1,
    band_hydrogen_phosphate_0,
    band_hydrogen_phosphate_3,
    band_iron_phosphate_1,
    band_iron_phosphate_2,
    band_calcium_phosphate_0,
    band_calcium_phosphate_1,
    band_calcium_phosphate_2,
    band_magnesium_phosphate_1,
};

/// Source-ordered water, heat, gas, and root diagnostics, HOUR1 3560-3573.
pub const BoundaryFlux = enum {
    liquid_water,
    macropore_water,
    liquid_heat,
    vapor_water,
    vapor_heat,
    carbon_dioxide,
    methane,
    oxygen,
    nitrogen,
    nitrous_oxide,
    ammonia,
    nitrogen_balance,
    hydrogen,
    root_length_density,
};

pub const LayerFluxes = struct {
    mineral_transformations: []f64,
    exchange_fluxes: []f64,
    organic_carbon_exchange: []f64, // XOCFXS
    organic_nitrogen_exchange: []f64, // XONFXS
    organic_phosphorus_exchange: []f64, // XOPFXS
    organic_acetate_exchange: []f64, // XOAFXS
    boundary_fluxes: []f64,
};

pub const ResetError = error{DimensionMismatch};

/// Translates HOUR1 lines 3431-3573. Fluxes retain their legacy
/// species-mass-per-timestep units; root length density is m m-3.
pub fn reset(
    layer_position: LayerPosition,
    organic_pool_count: usize,
    fluxes: *LayerFluxes,
) ResetError!void {
    if (layer_position == .surface) return;
    if (fluxes.mineral_transformations.len != std.meta.fields(MineralTransformation).len or
        fluxes.exchange_fluxes.len != std.meta.fields(ExchangeFlux).len or
        fluxes.organic_carbon_exchange.len != organic_pool_count or
        fluxes.organic_nitrogen_exchange.len != organic_pool_count or
        fluxes.organic_phosphorus_exchange.len != organic_pool_count or
        fluxes.organic_acetate_exchange.len != organic_pool_count or
        fluxes.boundary_fluxes.len != std.meta.fields(BoundaryFlux).len)
    {
        return error.DimensionMismatch;
    }

    @memset(fluxes.mineral_transformations, 0.0);
    @memset(fluxes.exchange_fluxes, 0.0);
    for (0..organic_pool_count) |pool| {
        fluxes.organic_carbon_exchange[pool] = 0.0;
        fluxes.organic_nitrogen_exchange[pool] = 0.0;
        fluxes.organic_phosphorus_exchange[pool] = 0.0;
        fluxes.organic_acetate_exchange[pool] = 0.0;
    }
    @memset(fluxes.boundary_fluxes, 0.0);
}

test "subsurface reset clears all named groups and runtime organic pools" {
    const allocator = std.testing.allocator;
    var fluxes = LayerFluxes{
        .mineral_transformations = try allocator.alloc(f64, std.meta.fields(MineralTransformation).len),
        .exchange_fluxes = try allocator.alloc(f64, std.meta.fields(ExchangeFlux).len),
        .organic_carbon_exchange = try allocator.alloc(f64, 7),
        .organic_nitrogen_exchange = try allocator.alloc(f64, 7),
        .organic_phosphorus_exchange = try allocator.alloc(f64, 7),
        .organic_acetate_exchange = try allocator.alloc(f64, 7),
        .boundary_fluxes = try allocator.alloc(f64, std.meta.fields(BoundaryFlux).len),
    };
    defer inline for (std.meta.fields(LayerFluxes)) |field| {
        allocator.free(@field(fluxes, field.name));
    };
    inline for (std.meta.fields(LayerFluxes)) |field| {
        @memset(@field(fluxes, field.name), 9.0);
    }

    try reset(.subsurface, 7, &fluxes);

    inline for (std.meta.fields(LayerFluxes)) |field| {
        for (@field(fluxes, field.name)) |value| {
            try std.testing.expectEqual(@as(f64, 0.0), value);
        }
    }
}

test "dimension mismatch fails before mutation" {
    var short = [_]f64{8.0};
    var fluxes: LayerFluxes = undefined;
    fluxes.mineral_transformations = &short;
    try std.testing.expectError(error.DimensionMismatch, reset(.subsurface, 0, &fluxes));
    try std.testing.expectEqual(@as(f64, 8.0), short[0]);
}
