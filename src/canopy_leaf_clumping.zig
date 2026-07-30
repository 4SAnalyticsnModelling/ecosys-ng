const std = @import("std");

pub const PlantActivity = enum {
    inactive,
    active,
};

pub const Dimensions = struct {
    branch_count: usize,
    node_count: usize,
    canopy_layer_count: usize,
    inclination_count: usize,

    fn leafCount(self: Dimensions) !usize {
        if (self.branch_count == 0 or self.node_count == 0 or
            self.canopy_layer_count == 0 or self.inclination_count == 0)
            return error.InvalidCanopyLeafClumpingDimensions;
        return std.math.mul(usize, self.branch_count, self.node_count) catch
            return error.InvalidCanopyLeafClumpingDimensions;
    }

    fn surfaceCount(self: Dimensions) !usize {
        const leaf_count = try self.leafCount();
        const with_layers = std.math.mul(usize, leaf_count, self.canopy_layer_count) catch
            return error.InvalidCanopyLeafClumpingDimensions;
        return std.math.mul(usize, with_layers, self.inclination_count) catch
            return error.InvalidCanopyLeafClumpingDimensions;
    }
};

pub const Inputs = struct {
    dimensions: Dimensions,
    activity: PlantActivity,
    population_count: f64,
    leaf_area_m2_per_m2_by_branch_node: []const f64,
    leaf_protein_g_per_m2_by_branch_node: []const f64,
    leaf_area_significance_threshold_m2_per_m2: f64,
    unclumped_surface_area_m2_per_m2: []const f64,
    clumping_factor: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    dimensions: Dimensions,
    last_qualifying_leaf_node_one_based_by_branch: []usize,
    clumped_surface_area_m2_per_m2: []f64,

    pub fn init(allocator: std.mem.Allocator, dimensions: Dimensions) !State {
        _ = try dimensions.leafCount();
        const surface_count = try dimensions.surfaceCount();
        const nodes = try allocator.alloc(usize, dimensions.branch_count);
        errdefer allocator.free(nodes);
        @memset(nodes, 0);
        const surface = try allocator.alloc(f64, surface_count);
        @memset(surface, 0);
        return .{
            .allocator = allocator,
            .dimensions = dimensions,
            .last_qualifying_leaf_node_one_based_by_branch = nodes,
            .clumped_surface_area_m2_per_m2 = surface,
        };
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.last_qualifying_leaf_node_one_based_by_branch);
        self.allocator.free(self.clumped_surface_area_m2_per_m2);
        self.* = undefined;
    }
};

/// UPTAKE.F 428--463. The source traverses nodes upward, so KLEAFX retains
/// the last qualifying one-based node despite its historical "lowest" label.
pub fn apply(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    if (inputs.activity == .inactive or inputs.population_count <= 0) return;

    var staged = try State.init(state.allocator, state.dimensions);
    errdefer staged.deinit();
    @memcpy(
        staged.last_qualifying_leaf_node_one_based_by_branch,
        state.last_qualifying_leaf_node_one_based_by_branch,
    );
    for (0..inputs.dimensions.branch_count) |branch| {
        for (0..inputs.dimensions.node_count) |node| {
            const leaf = branch * inputs.dimensions.node_count + node;
            if (inputs.leaf_area_m2_per_m2_by_branch_node[leaf] >
                inputs.leaf_area_significance_threshold_m2_per_m2 and
                inputs.leaf_protein_g_per_m2_by_branch_node[leaf] >
                    inputs.leaf_area_significance_threshold_m2_per_m2)
                staged.last_qualifying_leaf_node_one_based_by_branch[branch] = node + 1;
            var layer = inputs.dimensions.canopy_layer_count;
            while (layer > 0) {
                layer -= 1;
                for (0..inputs.dimensions.inclination_count) |inclination| {
                    const index = surfaceIndex(inputs.dimensions, inclination, layer, node, branch);
                    staged.clumped_surface_area_m2_per_m2[index] =
                        inputs.unclumped_surface_area_m2_per_m2[index] *
                        inputs.clumping_factor;
                }
            }
        }
    }
    for (staged.clumped_surface_area_m2_per_m2) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteCanopyLeafClumping;
    state.deinit();
    state.* = staged;
}

fn surfaceIndex(
    dimensions: Dimensions,
    inclination: usize,
    layer: usize,
    node: usize,
    branch: usize,
) usize {
    return (((branch * dimensions.node_count + node) *
        dimensions.canopy_layer_count + layer) *
        dimensions.inclination_count + inclination);
}

fn validate(state: *const State, inputs: Inputs) !void {
    const leaf_count = try inputs.dimensions.leafCount();
    const surface_count = try inputs.dimensions.surfaceCount();
    if (!std.meta.eql(state.dimensions, inputs.dimensions) or
        state.last_qualifying_leaf_node_one_based_by_branch.len != inputs.dimensions.branch_count or
        state.clumped_surface_area_m2_per_m2.len != surface_count or
        inputs.leaf_area_m2_per_m2_by_branch_node.len != leaf_count or
        inputs.leaf_protein_g_per_m2_by_branch_node.len != leaf_count or
        inputs.unclumped_surface_area_m2_per_m2.len != surface_count)
        return error.InvalidCanopyLeafClumpingDimensions;
    inline for (.{
        inputs.leaf_area_m2_per_m2_by_branch_node,
        inputs.leaf_protein_g_per_m2_by_branch_node,
        inputs.unclumped_surface_area_m2_per_m2,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value)) return error.InvalidCanopyLeafClumpingInput;
    inline for (.{
        inputs.population_count,
        inputs.leaf_area_significance_threshold_m2_per_m2,
        inputs.clumping_factor,
    }) |value| if (!std.math.isFinite(value))
        return error.InvalidCanopyLeafClumpingInput;
    if (inputs.population_count < 0 or
        inputs.leaf_area_significance_threshold_m2_per_m2 < 0 or
        inputs.clumping_factor < 0)
        return error.InvalidCanopyLeafClumpingInput;
}

test "UPTAKE active species retains last qualifying node and clumps runtime axes" {
    const dimensions = Dimensions{
        .branch_count = 2,
        .node_count = 3,
        .canopy_layer_count = 2,
        .inclination_count = 5,
    };
    var state = try State.init(std.testing.allocator, dimensions);
    defer state.deinit();
    var surface: [60]f64 = undefined;
    for (&surface, 0..) |*value, index| value.* = @floatFromInt(index + 1);
    try apply(&state, .{
        .dimensions = dimensions,
        .activity = .active,
        .population_count = 10,
        .leaf_area_m2_per_m2_by_branch_node = &.{ 0.2, 0, 0.4, 0, 0.5, 0.7 },
        .leaf_protein_g_per_m2_by_branch_node = &.{ 1, 1, 1, 1, 0, 2 },
        .leaf_area_significance_threshold_m2_per_m2 = 0.1,
        .unclumped_surface_area_m2_per_m2 = &surface,
        .clumping_factor = 0.25,
    });
    try std.testing.expectEqual([2]usize{ 3, 3 }, state.last_qualifying_leaf_node_one_based_by_branch[0..2].*);
    for (surface, state.clumped_surface_area_m2_per_m2) |source, clumped|
        try std.testing.expectEqual(source * 0.25, clumped);
}

test "inactive or zero-population species leaves prior clumping state unchanged" {
    const dimensions = Dimensions{
        .branch_count = 1,
        .node_count = 1,
        .canopy_layer_count = 1,
        .inclination_count = 1,
    };
    var state = try State.init(std.testing.allocator, dimensions);
    defer state.deinit();
    state.last_qualifying_leaf_node_one_based_by_branch[0] = 8;
    state.clumped_surface_area_m2_per_m2[0] = 9;
    try apply(&state, .{
        .dimensions = dimensions,
        .activity = .active,
        .population_count = 0,
        .leaf_area_m2_per_m2_by_branch_node = &.{1},
        .leaf_protein_g_per_m2_by_branch_node = &.{1},
        .leaf_area_significance_threshold_m2_per_m2 = 0,
        .unclumped_surface_area_m2_per_m2 = &.{4},
        .clumping_factor = 0.5,
    });
    try std.testing.expectEqual(@as(usize, 8), state.last_qualifying_leaf_node_one_based_by_branch[0]);
    try std.testing.expectEqual(@as(f64, 9), state.clumped_surface_area_m2_per_m2[0]);
}

test "late clumped-area overflow leaves canopy state unchanged" {
    const dimensions = Dimensions{
        .branch_count = 1,
        .node_count = 1,
        .canopy_layer_count = 1,
        .inclination_count = 1,
    };
    var state = try State.init(std.testing.allocator, dimensions);
    defer state.deinit();
    state.clumped_surface_area_m2_per_m2[0] = 6;
    try std.testing.expectError(
        error.NonFiniteCanopyLeafClumping,
        apply(&state, .{
            .dimensions = dimensions,
            .activity = .active,
            .population_count = 1,
            .leaf_area_m2_per_m2_by_branch_node = &.{1},
            .leaf_protein_g_per_m2_by_branch_node = &.{1},
            .leaf_area_significance_threshold_m2_per_m2 = 0,
            .unclumped_surface_area_m2_per_m2 = &.{std.math.floatMax(f64)},
            .clumping_factor = 2,
        }),
    );
    try std.testing.expectEqual(@as(f64, 6), state.clumped_surface_area_m2_per_m2[0]);
}
