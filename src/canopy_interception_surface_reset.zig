const std = @import("std");

pub const Dimensions = struct {
    species_count: usize,
    layer_count: usize,
    inclination_count: usize,
};

pub const SurfaceAreas = struct {
    total_leaf_area_m2_per_m2: []f64,
    total_branch_area_m2_per_m2: []f64,
    total_standing_dead_area_m2_per_m2: []f64,
};

/// HOUR1 lines 1149--1155. Resets species -> layer -> inclination surface
/// arrays in the exact source loop and assignment order.
pub fn apply(dimensions: Dimensions, surfaces: SurfaceAreas) !void {
    const expected_count = try validate(dimensions, surfaces);
    _ = expected_count;
    for (0..dimensions.species_count) |species| {
        for (0..dimensions.layer_count) |layer| {
            for (0..dimensions.inclination_count) |inclination| {
                const index = (species * dimensions.layer_count + layer) *
                    dimensions.inclination_count + inclination;
                surfaces.total_leaf_area_m2_per_m2[index] = 0.0;
                surfaces.total_branch_area_m2_per_m2[index] = 0.0;
                surfaces.total_standing_dead_area_m2_per_m2[index] = 0.0;
            }
        }
    }
}

fn validate(dimensions: Dimensions, surfaces: SurfaceAreas) !usize {
    if (dimensions.species_count == 0 or dimensions.layer_count == 0 or
        dimensions.inclination_count == 0)
        return error.ZeroCanopyInterceptionSurfaceExtent;
    const species_layers = try std.math.mul(
        usize,
        dimensions.species_count,
        dimensions.layer_count,
    );
    const expected_count = try std.math.mul(
        usize,
        species_layers,
        dimensions.inclination_count,
    );
    if (surfaces.total_leaf_area_m2_per_m2.len != expected_count or
        surfaces.total_branch_area_m2_per_m2.len != expected_count or
        surfaces.total_standing_dead_area_m2_per_m2.len != expected_count)
        return error.CanopyInterceptionSurfaceDimensionMismatch;
    return expected_count;
}

test "runtime species layer inclination surfaces reset in source order" {
    const dimensions: Dimensions = .{
        .species_count = 2,
        .layer_count = 3,
        .inclination_count = 5,
    };
    var leaf = [_]f64{1} ** 30;
    var branch = [_]f64{2} ** 30;
    var dead = [_]f64{3} ** 30;
    try apply(dimensions, .{
        .total_leaf_area_m2_per_m2 = &leaf,
        .total_branch_area_m2_per_m2 = &branch,
        .total_standing_dead_area_m2_per_m2 = &dead,
    });
    for (leaf) |value| try std.testing.expectEqual(@as(f64, 0), value);
    for (branch) |value| try std.testing.expectEqual(@as(f64, 0), value);
    for (dead) |value| try std.testing.expectEqual(@as(f64, 0), value);
}

test "dimension mismatch leaves all caller arrays unchanged" {
    var leaf = [_]f64{ 41, 41 };
    var branch = [_]f64{ 42, 42 };
    var dead = [_]f64{ 43, 43 };
    try std.testing.expectError(
        error.CanopyInterceptionSurfaceDimensionMismatch,
        apply(.{
            .species_count = 1,
            .layer_count = 1,
            .inclination_count = 3,
        }, .{
            .total_leaf_area_m2_per_m2 = &leaf,
            .total_branch_area_m2_per_m2 = &branch,
            .total_standing_dead_area_m2_per_m2 = &dead,
        }),
    );
    try std.testing.expectEqualSlices(f64, &.{ 41, 41 }, &leaf);
    try std.testing.expectEqualSlices(f64, &.{ 42, 42 }, &branch);
    try std.testing.expectEqualSlices(f64, &.{ 43, 43 }, &dead);
}
