const std = @import("std");

/// Exact grosub.f lines 957--958 node admission in source K order.
///
/// Leaf areas and the plant-specific ZEROP threshold are m2. Scratch and
/// destination are caller allocated for the runtime node count so a late
/// invalid node cannot partially publish the admission mask.
pub fn computeRuntimeNodes(
    leaf_area_m2_by_node: []const f64,
    leaf_area_presence_threshold_m2: f64,
    scratch: []bool,
    destination: []bool,
) !void {
    const node_count = leaf_area_m2_by_node.len;
    if (node_count == 0 or
        scratch.len != node_count or
        destination.len != node_count)
        return error.CanopyNodeCarboxylationAdmissionDimensionMismatch;
    if (!std.math.isFinite(leaf_area_presence_threshold_m2) or
        leaf_area_presence_threshold_m2 < 0)
        return error.InvalidCanopyNodeLeafAreaThreshold;

    for (leaf_area_m2_by_node, scratch) |leaf_area_m2, *candidate| {
        if (!std.math.isFinite(leaf_area_m2) or leaf_area_m2 < 0)
            return error.InvalidCanopyNodeLeafArea;
        candidate.* = leaf_area_m2 > leaf_area_presence_threshold_m2;
    }
    @memcpy(destination, scratch);
}

test "GROSUB admits arbitrary runtime nodes in source K order" {
    const allocator = std.testing.allocator;
    const node_count = 39;
    const leaf_area = try allocator.alloc(f64, node_count);
    defer allocator.free(leaf_area);
    const scratch = try allocator.alloc(bool, node_count);
    defer allocator.free(scratch);
    const admitted = try allocator.alloc(bool, node_count);
    defer allocator.free(admitted);
    for (leaf_area, 0..) |*area_m2, node|
        area_m2.* = @as(f64, @floatFromInt(node)) * 0.01;

    try computeRuntimeNodes(leaf_area, 0.2, scratch, admitted);

    for (admitted, 0..) |is_admitted, node|
        try std.testing.expectEqual(node > 20, is_admitted);
}

test "source node leaf-area comparison remains strict at ZEROP" {
    const area = [_]f64{ 0.099, 0.1, 0.101 };
    var scratch: [3]bool = undefined;
    var admitted: [3]bool = undefined;
    try computeRuntimeNodes(&area, 0.1, &scratch, &admitted);
    try std.testing.expectEqualSlices(
        bool,
        &.{ false, false, true },
        &admitted,
    );
}

test "invalid late node leaves published admission unchanged" {
    const area = [_]f64{ 1, std.math.nan(f64) };
    var scratch: [2]bool = undefined;
    var admitted = [_]bool{ false, true };
    try std.testing.expectError(
        error.InvalidCanopyNodeLeafArea,
        computeRuntimeNodes(&area, 0, &scratch, &admitted),
    );
    try std.testing.expectEqualSlices(bool, &.{ false, true }, &admitted);
}

test "runtime node topology mismatch fails explicitly" {
    const area = [_]f64{ 1, 2 };
    var short_scratch: [1]bool = undefined;
    var admitted = [_]bool{ false, false };
    try std.testing.expectError(
        error.CanopyNodeCarboxylationAdmissionDimensionMismatch,
        computeRuntimeNodes(&area, 0, &short_scratch, &admitted),
    );
}
