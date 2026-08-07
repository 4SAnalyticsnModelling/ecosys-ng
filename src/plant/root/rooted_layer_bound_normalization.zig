const std = @import("std");

pub const Dimensions = struct {
    species_count: usize,
    root_population_count: usize,
};

pub const RootedLayerBounds = struct {
    uppermost_rooted_layer: []usize,
    root_search_layer: []usize,
    deepest_rooted_layer_by_population: []usize,
};

/// `hour1.f` lines 1871--1875. Preserves species-outer then root-population
/// traversal while replacing the source fixed ten root populations.
pub fn apply(
    dimensions: Dimensions,
    uppermost_soil_layer: usize,
    soil_layer_count: usize,
    bounds: RootedLayerBounds,
) !void {
    try validate(
        dimensions,
        uppermost_soil_layer,
        soil_layer_count,
        bounds,
    );
    for (0..dimensions.species_count) |species| {
        bounds.uppermost_rooted_layer[species] = @max(
            bounds.uppermost_rooted_layer[species],
            uppermost_soil_layer,
        );
        bounds.root_search_layer[species] = @max(
            bounds.root_search_layer[species],
            uppermost_soil_layer,
        );
        for (0..dimensions.root_population_count) |root_population| {
            const index =
                species * dimensions.root_population_count + root_population;
            bounds.deepest_rooted_layer_by_population[index] = @max(
                bounds.deepest_rooted_layer_by_population[index],
                uppermost_soil_layer,
            );
        }
    }
}

fn validate(
    dimensions: Dimensions,
    uppermost_soil_layer: usize,
    soil_layer_count: usize,
    bounds: RootedLayerBounds,
) !void {
    if (dimensions.species_count == 0 or
        dimensions.root_population_count == 0)
        return error.ZeroRootedLayerBoundExtent;
    if (soil_layer_count == 0 or uppermost_soil_layer >= soil_layer_count)
        return error.UppermostSoilLayerOutOfRange;
    const population_value_count = try std.math.mul(
        usize,
        dimensions.species_count,
        dimensions.root_population_count,
    );
    if (bounds.uppermost_rooted_layer.len != dimensions.species_count or
        bounds.root_search_layer.len != dimensions.species_count or
        bounds.deepest_rooted_layer_by_population.len != population_value_count)
        return error.RootedLayerBoundDimensionMismatch;
    inline for (.{
        bounds.uppermost_rooted_layer,
        bounds.root_search_layer,
        bounds.deepest_rooted_layer_by_population,
    }) |values| for (values) |value|
        if (value >= soil_layer_count)
            return error.RootedLayerBoundOutOfRange;
}

test "runtime species and root populations are bounded by uppermost soil" {
    var uppermost = [_]usize{ 0, 3 };
    var search = [_]usize{ 1, 4 };
    var deepest = [_]usize{ 0, 1, 2, 3, 4, 0 };
    try apply(.{
        .species_count = 2,
        .root_population_count = 3,
    }, 2, 6, .{
        .uppermost_rooted_layer = &uppermost,
        .root_search_layer = &search,
        .deepest_rooted_layer_by_population = &deepest,
    });
    try std.testing.expectEqualSlices(usize, &.{ 2, 3 }, &uppermost);
    try std.testing.expectEqualSlices(usize, &.{ 2, 4 }, &search);
    try std.testing.expectEqualSlices(
        usize,
        &.{ 2, 2, 2, 3, 4, 2 },
        &deepest,
    );
}

test "dimension mismatch leaves all bounds unchanged" {
    var uppermost = [_]usize{1};
    var search = [_]usize{2};
    var deepest = [_]usize{3};
    try std.testing.expectError(
        error.RootedLayerBoundDimensionMismatch,
        apply(.{
            .species_count = 1,
            .root_population_count = 2,
        }, 0, 4, .{
            .uppermost_rooted_layer = &uppermost,
            .root_search_layer = &search,
            .deepest_rooted_layer_by_population = &deepest,
        }),
    );
    try std.testing.expectEqual(@as(usize, 1), uppermost[0]);
    try std.testing.expectEqual(@as(usize, 2), search[0]);
    try std.testing.expectEqual(@as(usize, 3), deepest[0]);
}
