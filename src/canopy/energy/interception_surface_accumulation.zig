const std = @import("std");

pub const Dimensions = struct {
    species_count: usize,
    layer_count: usize,
    inclination_count: usize,
    branch_capacity: usize,
    node_count: usize,
};

pub const Inputs = struct {
    dimensions: Dimensions,
    branch_count_by_species: []const usize,
    layer_bottom_height_m: []const f64,
    snow_depth_m: f64,
    surface_water_ice_depth_m: f64,
    depth_tolerance_m: f64,
    /// Species -> layer -> inclination.
    standing_dead_area_m2_per_m2: []const f64,
    /// Species -> layer -> branch -> inclination.
    branch_area_m2_per_m2: []const f64,
    /// Species -> layer -> branch -> inclination -> node.
    node_leaf_area_m2_per_m2: []const f64,
};

pub const Totals = struct {
    /// Species -> layer -> inclination.
    leaf_area_m2_per_m2: []f64,
    branch_area_m2_per_m2: []f64,
    standing_dead_area_m2_per_m2: []f64,
};

/// `hour1.f` lines 1156--1172. Preserves NZ -> L -> dead -> NB -> N -> K
/// source order. Standing dead remains assigned to the final inclination
/// class, corresponding to source class four.
pub fn apply(inputs: Inputs, totals: Totals) !void {
    try validate(inputs, totals);
    const dimensions = inputs.dimensions;
    for (0..dimensions.species_count) |species| {
        for (0..dimensions.layer_count) |layer| {
            if (inputs.layer_bottom_height_m[layer] >
                inputs.snow_depth_m - inputs.depth_tolerance_m and
                inputs.layer_bottom_height_m[layer] >
                    inputs.surface_water_ice_depth_m - inputs.depth_tolerance_m)
            {
                const dead_index = surfaceIndex(
                    dimensions,
                    species,
                    layer,
                    dimensions.inclination_count - 1,
                );
                totals.standing_dead_area_m2_per_m2[dead_index] +=
                    inputs.standing_dead_area_m2_per_m2[dead_index];
                for (0..inputs.branch_count_by_species[species]) |branch| {
                    for (0..dimensions.inclination_count) |inclination| {
                        const surface_index = surfaceIndex(
                            dimensions,
                            species,
                            layer,
                            inclination,
                        );
                        totals.branch_area_m2_per_m2[surface_index] +=
                            inputs.branch_area_m2_per_m2[
                                branchIndex(
                                    dimensions,
                                    species,
                                    layer,
                                    branch,
                                    inclination,
                                )
                            ];
                        for (0..dimensions.node_count) |node| {
                            totals.leaf_area_m2_per_m2[surface_index] +=
                                inputs.node_leaf_area_m2_per_m2[
                                    nodeIndex(
                                        dimensions,
                                        species,
                                        layer,
                                        branch,
                                        inclination,
                                        node,
                                    )
                                ];
                        }
                    }
                }
            }
        }
    }
}

fn validate(inputs: Inputs, totals: Totals) !void {
    const dimensions = inputs.dimensions;
    inline for (.{
        dimensions.species_count,
        dimensions.layer_count,
        dimensions.inclination_count,
        dimensions.branch_capacity,
        dimensions.node_count,
    }) |extent| if (extent == 0)
        return error.ZeroCanopySurfaceAccumulationExtent;
    if (inputs.branch_count_by_species.len != dimensions.species_count or
        inputs.layer_bottom_height_m.len != dimensions.layer_count)
        return error.CanopySurfaceAccumulationDimensionMismatch;
    for (inputs.branch_count_by_species) |count|
        if (count > dimensions.branch_capacity)
            return error.CanopyBranchCountExceedsCapacity;
    const surface_count = try std.math.mul(
        usize,
        try std.math.mul(usize, dimensions.species_count, dimensions.layer_count),
        dimensions.inclination_count,
    );
    const branch_count = try std.math.mul(
        usize,
        surface_count,
        dimensions.branch_capacity,
    );
    const node_value_count = try std.math.mul(
        usize,
        branch_count,
        dimensions.node_count,
    );
    if (inputs.standing_dead_area_m2_per_m2.len != surface_count or
        inputs.branch_area_m2_per_m2.len != branch_count or
        inputs.node_leaf_area_m2_per_m2.len != node_value_count or
        totals.leaf_area_m2_per_m2.len != surface_count or
        totals.branch_area_m2_per_m2.len != surface_count or
        totals.standing_dead_area_m2_per_m2.len != surface_count)
        return error.CanopySurfaceAccumulationDimensionMismatch;
    inline for (.{
        inputs.snow_depth_m,
        inputs.surface_water_ice_depth_m,
        inputs.depth_tolerance_m,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopySurfaceAccumulationInput;
    inline for (.{
        inputs.layer_bottom_height_m,
        inputs.standing_dead_area_m2_per_m2,
        inputs.branch_area_m2_per_m2,
        inputs.node_leaf_area_m2_per_m2,
        totals.leaf_area_m2_per_m2,
        totals.branch_area_m2_per_m2,
        totals.standing_dead_area_m2_per_m2,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopySurfaceAccumulationInput;
}

fn surfaceIndex(
    dimensions: Dimensions,
    species: usize,
    layer: usize,
    inclination: usize,
) usize {
    return (species * dimensions.layer_count + layer) *
        dimensions.inclination_count + inclination;
}

fn branchIndex(
    dimensions: Dimensions,
    species: usize,
    layer: usize,
    branch: usize,
    inclination: usize,
) usize {
    return ((species * dimensions.layer_count + layer) *
        dimensions.branch_capacity + branch) *
        dimensions.inclination_count + inclination;
}

fn nodeIndex(
    dimensions: Dimensions,
    species: usize,
    layer: usize,
    branch: usize,
    inclination: usize,
    node: usize,
) usize {
    return branchIndex(dimensions, species, layer, branch, inclination) *
        dimensions.node_count + node;
}

test "only exposed layers accumulate dead branch and node surfaces" {
    const dimensions: Dimensions = .{
        .species_count = 1,
        .layer_count = 2,
        .inclination_count = 3,
        .branch_capacity = 2,
        .node_count = 2,
    };
    var dead = [_]f64{0} ** 6;
    dead[5] = 5;
    const branch = [_]f64{1} ** 12;
    const nodes = [_]f64{0.5} ** 24;
    var leaf_total = [_]f64{0} ** 6;
    var branch_total = [_]f64{0} ** 6;
    var dead_total = [_]f64{0} ** 6;
    try apply(.{
        .dimensions = dimensions,
        .branch_count_by_species = &.{2},
        .layer_bottom_height_m = &.{ 0, 2 },
        .snow_depth_m = 1,
        .surface_water_ice_depth_m = 0.5,
        .depth_tolerance_m = 0.1,
        .standing_dead_area_m2_per_m2 = &dead,
        .branch_area_m2_per_m2 = &branch,
        .node_leaf_area_m2_per_m2 = &nodes,
    }, .{
        .leaf_area_m2_per_m2 = &leaf_total,
        .branch_area_m2_per_m2 = &branch_total,
        .standing_dead_area_m2_per_m2 = &dead_total,
    });
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, leaf_total[0..3]);
    try std.testing.expectEqualSlices(f64, &.{ 2, 2, 2 }, leaf_total[3..6]);
    try std.testing.expectEqualSlices(f64, &.{ 2, 2, 2 }, branch_total[3..6]);
    try std.testing.expectEqual(@as(f64, 5), dead_total[5]);
}

test "branch count beyond runtime capacity fails before mutation" {
    var values = [_]f64{42};
    try std.testing.expectError(error.CanopyBranchCountExceedsCapacity, apply(.{
        .dimensions = .{
            .species_count = 1,
            .layer_count = 1,
            .inclination_count = 1,
            .branch_capacity = 1,
            .node_count = 1,
        },
        .branch_count_by_species = &.{2},
        .layer_bottom_height_m = &.{1},
        .snow_depth_m = 0,
        .surface_water_ice_depth_m = 0,
        .depth_tolerance_m = 0,
        .standing_dead_area_m2_per_m2 = &.{0},
        .branch_area_m2_per_m2 = &.{0},
        .node_leaf_area_m2_per_m2 = &.{0},
    }, .{
        .leaf_area_m2_per_m2 = &values,
        .branch_area_m2_per_m2 = &values,
        .standing_dead_area_m2_per_m2 = &values,
    }));
    try std.testing.expectEqual(@as(f64, 42), values[0]);
}
