const std = @import("std");
const compute = @import("compute.zig");
const spatial_grid = @import("spatial_grid.zig");

/// Serial out-of-core tile transaction:
/// load owned region plus two-cell halo, run parallel owned-cell kernels, then
/// commit only the owned interior. A failed kernel never commits its tile.
pub fn run(
    executor: compute.CpuExecutor,
    plan: spatial_grid.TilePlan,
    context: anytype,
    comptime loadTile: anytype,
    comptime cellKernel: anytype,
    comptime commitOwnedTile: anytype,
) !void {
    var previous_z_order_index: ?u64 = null;
    for (plan.tiles, 0..) |tile, tile_index| {
        if (previous_z_order_index) |previous| {
            if (tile.z_order_index <= previous)
                return error.TilePlanIsNotStrictMortonOrder;
        }
        try loadTile(context, tile);
        const owned_cells = try plan.ownedCells(tile_index);
        try executor.runCells(owned_cells, context, cellKernel);
        try commitOwnedTile(context, tile);
        previous_z_order_index = tile.z_order_index;
    }
}

const TestContext = struct {
    values: []u8,
    events: []u64,
    next_event: usize = 0,
};

fn recordLoad(context: *TestContext, tile: spatial_grid.Tile) !void {
    context.events[context.next_event] = tile.z_order_index * 2;
    context.next_event += 1;
}

fn markCells(context: *TestContext, cells: []const usize) !void {
    for (cells) |cell| {
        if (context.values[cell] != 0) return error.CellProcessedMoreThanOnce;
        context.values[cell] = 1;
    }
}

fn recordCommit(context: *TestContext, tile: spatial_grid.Tile) !void {
    context.events[context.next_event] = tile.z_order_index * 2 + 1;
    context.next_event += 1;
}

test "tiles are serial transactions and owned cells run exactly once" {
    var plan = try spatial_grid.TilePlan.init(std.testing.allocator, 7, 9, 3, 4, 2);
    defer plan.deinit();
    const values = try std.testing.allocator.alloc(u8, 7 * 9);
    defer std.testing.allocator.free(values);
    @memset(values, 0);
    const events = try std.testing.allocator.alloc(u64, plan.tiles.len * 2);
    defer std.testing.allocator.free(events);
    var context = TestContext{ .values = values, .events = events };
    const executor = try compute.CpuExecutor.init(std.testing.allocator, 4, 64);
    try run(
        executor,
        plan,
        &context,
        recordLoad,
        markCells,
        recordCommit,
    );
    for (values) |value| try std.testing.expectEqual(@as(u8, 1), value);
    for (plan.tiles, 0..) |tile, index| {
        try std.testing.expectEqual(tile.z_order_index * 2, events[index * 2]);
        try std.testing.expectEqual(tile.z_order_index * 2 + 1, events[index * 2 + 1]);
    }
}
