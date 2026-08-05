const std = @import("std");

pub const Workspace = struct {
    total_co2_fixation_umol_per_s: *f64,
    total_carbohydrate_production_umol_per_s: *f64,
    c3_fixation_umol_per_s_by_node: []f64,
    c4_fixation_umol_per_s_by_node: []f64,
};

/// Exact GROSUB lines 943--948 carboxylation-workspace reset.
///
/// The two node slices replace the source's fixed CH2O3(1:25) and
/// CH2O4(1:25) arrays with caller-allocated runtime node dimensions. These
/// values are temporary fixation rates, not persistent nonstructural carbon.
pub fn reset(workspace: Workspace) !void {
    const node_count = workspace.c3_fixation_umol_per_s_by_node.len;
    if (node_count == 0 or
        workspace.c4_fixation_umol_per_s_by_node.len != node_count)
        return error.CanopyCarboxylationWorkspaceDimensionMismatch;

    workspace.total_co2_fixation_umol_per_s.* = 0;
    workspace.total_carbohydrate_production_umol_per_s.* = 0;
    for (
        workspace.c3_fixation_umol_per_s_by_node,
        workspace.c4_fixation_umol_per_s_by_node,
    ) |*c3_fixation_umol_per_s, *c4_fixation_umol_per_s| {
        c3_fixation_umol_per_s.* = 0;
        c4_fixation_umol_per_s.* = 0;
    }
}

test "GROSUB resets arbitrary runtime carboxylation nodes in K order" {
    const allocator = std.testing.allocator;
    const node_count = 43;
    const c3 = try allocator.alloc(f64, node_count);
    defer allocator.free(c3);
    const c4 = try allocator.alloc(f64, node_count);
    defer allocator.free(c4);
    for (c3, c4, 0..) |*c3_rate, *c4_rate, node| {
        c3_rate.* = @as(f64, @floatFromInt(node)) + 0.25;
        c4_rate.* = @as(f64, @floatFromInt(node)) + 0.75;
    }
    var total_co2: f64 = 123;
    var total_carbohydrate: f64 = 456;

    try reset(.{
        .total_co2_fixation_umol_per_s = &total_co2,
        .total_carbohydrate_production_umol_per_s = &total_carbohydrate,
        .c3_fixation_umol_per_s_by_node = c3,
        .c4_fixation_umol_per_s_by_node = c4,
    });

    try std.testing.expectEqual(@as(f64, 0), total_co2);
    try std.testing.expectEqual(@as(f64, 0), total_carbohydrate);
    for (c3, c4) |c3_rate, c4_rate| {
        try std.testing.expectEqual(@as(f64, 0), c3_rate);
        try std.testing.expectEqual(@as(f64, 0), c4_rate);
    }
}

test "dimension failure preserves all workspace values" {
    var total_co2: f64 = 1;
    var total_carbohydrate: f64 = 2;
    var c3 = [_]f64{ 3, 4 };
    var c4 = [_]f64{5};
    try std.testing.expectError(
        error.CanopyCarboxylationWorkspaceDimensionMismatch,
        reset(.{
            .total_co2_fixation_umol_per_s = &total_co2,
            .total_carbohydrate_production_umol_per_s = &total_carbohydrate,
            .c3_fixation_umol_per_s_by_node = &c3,
            .c4_fixation_umol_per_s_by_node = &c4,
        }),
    );
    try std.testing.expectEqual(@as(f64, 1), total_co2);
    try std.testing.expectEqual(@as(f64, 2), total_carbohydrate);
    try std.testing.expectEqualSlices(f64, &.{ 3, 4 }, &c3);
    try std.testing.expectEqualSlices(f64, &.{5}, &c4);
}
