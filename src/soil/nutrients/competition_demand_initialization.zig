const std = @import("std");

pub const Demand = struct {
    oxygen_g_o_per_h: f64 = 0,
    ammonium_nonband_g_n_per_h: f64 = 0,
    nitrate_nonband_g_n_per_h: f64 = 0,
    nitrite_nonband_g_n_per_h: f64 = 0,
    nitrous_oxide_nonband_g_n_per_h: f64 = 0,
    h2po4_nonband_g_p_per_h: f64 = 0,
    hpo4_nonband_g_p_per_h: f64 = 0,
    chemodenitrification_nitrite_nonband_g_n_per_h: f64 = 0,
    ammonium_band_g_n_per_h: f64 = 0,
    nitrate_band_g_n_per_h: f64 = 0,
    nitrite_band_g_n_per_h: f64 = 0,
    h2po4_band_g_p_per_h: f64 = 0,
    hpo4_band_g_p_per_h: f64 = 0,
    chemodenitrification_nitrite_band_g_n_per_h: f64 = 0,
};

pub const Dimensions = struct {
    column_count: usize,
    row_count: usize,
    soil_layer_count_including_surface: usize,
};

/// Resets competition demands before the first hourly process calculation.
///
/// Traceability: `starts.f` lines 1611--1624 (`ROXYX` through `RVMBC`).
/// Storage is cell-major, then layer, and includes the surface layer at
/// layer index zero. The source quantities are extensive hourly demands.
pub fn initialize(dimensions: Dimensions, demand_by_cell_layer: []Demand) !void {
    const expected_count = try checkedValueCount(dimensions);
    if (demand_by_cell_layer.len != expected_count)
        return error.InvalidSoilCompetitionDemandDimensions;

    @memset(demand_by_cell_layer, Demand{});
}

fn checkedValueCount(dimensions: Dimensions) !usize {
    if (dimensions.column_count == 0 or
        dimensions.row_count == 0 or
        dimensions.soil_layer_count_including_surface == 0)
        return error.InvalidSoilCompetitionDemandDimensions;

    const cell_count = std.math.mul(
        usize,
        dimensions.column_count,
        dimensions.row_count,
    ) catch return error.InvalidSoilCompetitionDemandDimensions;
    return std.math.mul(
        usize,
        cell_count,
        dimensions.soil_layer_count_including_surface,
    ) catch return error.InvalidSoilCompetitionDemandDimensions;
}

test "initialize resets every demand for runtime grid and layer dimensions" {
    var demands = [_]Demand{.{
        .oxygen_g_o_per_h = 1,
        .ammonium_nonband_g_n_per_h = 2,
        .nitrate_nonband_g_n_per_h = 3,
        .nitrite_nonband_g_n_per_h = 4,
        .nitrous_oxide_nonband_g_n_per_h = 5,
        .h2po4_nonband_g_p_per_h = 6,
        .hpo4_nonband_g_p_per_h = 7,
        .chemodenitrification_nitrite_nonband_g_n_per_h = 8,
        .ammonium_band_g_n_per_h = 9,
        .nitrate_band_g_n_per_h = 10,
        .nitrite_band_g_n_per_h = 11,
        .h2po4_band_g_p_per_h = 12,
        .hpo4_band_g_p_per_h = 13,
        .chemodenitrification_nitrite_band_g_n_per_h = 14,
    }} ** 12;

    try initialize(.{
        .column_count = 2,
        .row_count = 2,
        .soil_layer_count_including_surface = 3,
    }, &demands);

    for (demands) |demand|
        try std.testing.expectEqualDeep(Demand{}, demand);
}

test "initialize rejects a short runtime allocation without mutation" {
    var demands = [_]Demand{.{
        .oxygen_g_o_per_h = 7,
    }} ** 5;

    try std.testing.expectError(
        error.InvalidSoilCompetitionDemandDimensions,
        initialize(.{
            .column_count = 2,
            .row_count = 1,
            .soil_layer_count_including_surface = 3,
        }, &demands),
    );
    for (demands) |demand|
        try std.testing.expectEqual(@as(f64, 7), demand.oxygen_g_o_per_h);
}

test "initialize rejects zero and overflowing dimensions" {
    var none: [0]Demand = .{};
    try std.testing.expectError(
        error.InvalidSoilCompetitionDemandDimensions,
        initialize(.{
            .column_count = 0,
            .row_count = 1,
            .soil_layer_count_including_surface = 1,
        }, &none),
    );
    try std.testing.expectError(
        error.InvalidSoilCompetitionDemandDimensions,
        initialize(.{
            .column_count = std.math.maxInt(usize),
            .row_count = 2,
            .soil_layer_count_including_surface = 1,
        }, &none),
    );
}
