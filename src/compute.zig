const std = @import("std");

pub const CellRange = struct {
    first: usize,
    end: usize,

    pub fn count(self: CellRange) usize {
        return self.end - self.first;
    }
};

pub const CpuExecutor = struct {
    allocator: std.mem.Allocator,
    worker_count: usize,
    tile_cell_count: usize,

    pub fn init(allocator: std.mem.Allocator, worker_count: usize, tile_cell_count: usize) !CpuExecutor {
        if (worker_count == 0) return error.NoWorkerThreads;
        if (tile_cell_count == 0) return error.EmptyTile;
        return .{ .allocator = allocator, .worker_count = worker_count, .tile_cell_count = tile_cell_count };
    }

    /// Parallelizes the complete set of grid cells belonging to one
    /// already-loaded spatial tile. Spatial tile traversal is deliberately
    /// absent here and belongs only to `tile_executor`, which serializes
    /// load/compute/commit transactions.
    pub fn run(self: CpuExecutor, cell_count: usize, context: anytype, comptime kernel: anytype) !void {
        if (cell_count == 0) return error.EmptyGrid;
        try self.runTileParallel(
            context,
            kernel,
            .{ .first = 0, .end = cell_count },
        );
    }

    fn runTileParallel(self: CpuExecutor, context: anytype, comptime kernel: anytype, tile: CellRange) !void {
        const active_worker_count = @min(self.worker_count, tile.count());
        if (active_worker_count == 1) return kernel(context, tile);
        const threads = try self.allocator.alloc(std.Thread, active_worker_count);
        defer self.allocator.free(threads);
        const worker_errors = try self.allocator.alloc(?anyerror, active_worker_count);
        defer self.allocator.free(worker_errors);
        @memset(worker_errors, null);

        const Worker = struct {
            fn entry(
                worker_context: @TypeOf(context),
                errors: []?anyerror,
                tile_range: CellRange,
                worker_index: usize,
                workers: usize,
            ) void {
                const cells_per_worker = ceilDivision(tile_range.count(), workers) catch |err| {
                    errors[worker_index] = err;
                    return;
                };
                const local_first = std.math.mul(usize, worker_index, cells_per_worker) catch |err| {
                    errors[worker_index] = err;
                    return;
                };
                if (local_first >= tile_range.count()) return;
                const range = CellRange{
                    .first = tile_range.first + local_first,
                    .end = @min(tile_range.end, tile_range.first + local_first + cells_per_worker),
                };
                kernel(worker_context, range) catch |err| {
                    errors[worker_index] = err;
                };
            }
        };

        var spawned: usize = 0;
        errdefer for (threads[0..spawned]) |thread| thread.join();
        while (spawned < active_worker_count) : (spawned += 1) {
            threads[spawned] = try std.Thread.spawn(.{}, Worker.entry, .{
                context,
                worker_errors,
                tile,
                spawned,
                active_worker_count,
            });
        }
        for (threads) |thread| thread.join();
        // Disarm spawn-failure cleanup: all thread handles are now consumed.
        spawned = 0;
        for (worker_errors) |maybe_error| if (maybe_error) |err| return err;
    }

    /// Parallelizes an explicit Morton-ordered cell list belonging to one
    /// already-loaded tile. This function never advances to another tile.
    pub fn runCells(
        self: CpuExecutor,
        cells: []const usize,
        context: anytype,
        comptime kernel: anytype,
    ) !void {
        if (cells.len == 0) return error.EmptyTile;
        const active_worker_count = @min(self.worker_count, cells.len);
        if (active_worker_count == 1) return kernel(context, cells);
        const threads = try self.allocator.alloc(std.Thread, active_worker_count);
        defer self.allocator.free(threads);
        const worker_errors = try self.allocator.alloc(?anyerror, active_worker_count);
        defer self.allocator.free(worker_errors);
        @memset(worker_errors, null);
        const Worker = struct {
            fn entry(
                worker_context: @TypeOf(context),
                all_cells: []const usize,
                errors: []?anyerror,
                worker_index: usize,
                workers: usize,
            ) void {
                const per_worker = ceilDivision(all_cells.len, workers) catch |err| {
                    errors[worker_index] = err;
                    return;
                };
                const first = std.math.mul(usize, worker_index, per_worker) catch |err| {
                    errors[worker_index] = err;
                    return;
                };
                if (first >= all_cells.len) return;
                kernel(worker_context, all_cells[first..@min(all_cells.len, first + per_worker)]) catch |err| {
                    errors[worker_index] = err;
                };
            }
        };
        var spawned: usize = 0;
        errdefer for (threads[0..spawned]) |thread| thread.join();
        while (spawned < active_worker_count) : (spawned += 1) {
            threads[spawned] = try std.Thread.spawn(.{}, Worker.entry, .{
                context,
                cells,
                worker_errors,
                spawned,
                active_worker_count,
            });
        }
        for (threads) |thread| thread.join();
        spawned = 0;
        for (worker_errors) |maybe_error| if (maybe_error) |err| return err;
    }

    /// Adapts an existing contiguous `CellRange` kernel to an explicit
    /// non-contiguous Morton-owned cell list. Each invocation contains one
    /// global cell, preventing a range from spanning unowned halo cells.
    pub fn runIndexedCells(
        self: CpuExecutor,
        cells: []const usize,
        context: anytype,
        comptime cell_range_kernel: anytype,
    ) !void {
        const Adapter = struct {
            science_context: @TypeOf(context),

            noinline fn apply(
                adapter: *@This(),
                worker_cells: []const usize,
            ) !void {
                for (worker_cells) |cell| {
                    try @call(.never_inline, cell_range_kernel, .{
                        adapter.science_context,
                        CellRange{
                            .first = cell,
                            .end = try std.math.add(usize, cell, 1),
                        },
                    });
                }
            }
        };
        var adapter = Adapter{ .science_context = context };
        try self.runCells(cells, &adapter, Adapter.apply);
    }

    /// Assigns complete horizontal grid cells to workers while adapting a
    /// cell-layer kernel that consumes a contiguous flattened layer range.
    /// No worker boundary can split the vertical column of one grid cell.
    pub fn runCellLayers(
        self: CpuExecutor,
        cell_count: usize,
        layer_count_per_cell: usize,
        context: anytype,
        comptime layer_kernel: anytype,
    ) !void {
        if (cell_count == 0) return error.EmptyGrid;
        if (layer_count_per_cell == 0) return error.NoSoilLayers;
        const Adapter = struct {
            science_context: @TypeOf(context),
            layers_per_cell: usize,

            fn apply(adapter: *@This(), cells: CellRange) !void {
                const first_layer = try std.math.mul(
                    usize,
                    cells.first,
                    adapter.layers_per_cell,
                );
                const end_layer = try std.math.mul(
                    usize,
                    cells.end,
                    adapter.layers_per_cell,
                );
                try layer_kernel(
                    adapter.science_context,
                    .{ .first = first_layer, .end = end_layer },
                );
            }
        };
        var adapter = Adapter{
            .science_context = context,
            .layers_per_cell = layer_count_per_cell,
        };
        // This is one already-loaded spatial tile. Parallelism is across its
        // grid cells; this method never advances to another tile.
        try self.runTileParallel(
            &adapter,
            Adapter.apply,
            .{ .first = 0, .end = cell_count },
        );
    }

    /// Parallelizes the Morton-ordered owned cells of one loaded spatial tile
    /// while preserving each cell's complete vertical column. Unlike
    /// `runCellLayers`, this accepts non-contiguous global cell indices, which
    /// is required when a rectangular tile occupies only part of each domain
    /// row. Tiles themselves remain serial in `tile_executor`.
    pub fn runIndexedCellLayers(
        self: CpuExecutor,
        cells: []const usize,
        layer_count_per_cell: usize,
        context: anytype,
        comptime cell_layer_kernel: anytype,
    ) !void {
        if (cells.len == 0) return error.EmptyTile;
        if (layer_count_per_cell == 0) return error.NoSoilLayers;
        const Adapter = struct {
            science_context: @TypeOf(context),
            layers_per_cell: usize,

            noinline fn apply(
                adapter: *@This(),
                worker_cells: []const usize,
            ) !void {
                for (worker_cells) |cell| {
                    const first_layer = try std.math.mul(
                        usize,
                        cell,
                        adapter.layers_per_cell,
                    );
                    try @call(.never_inline, cell_layer_kernel, .{
                        adapter.science_context,
                        CellRange{
                            .first = first_layer,
                            .end = try std.math.add(
                                usize,
                                first_layer,
                                adapter.layers_per_cell,
                            ),
                        },
                    });
                }
            }
        };
        var adapter = Adapter{
            .science_context = context,
            .layers_per_cell = layer_count_per_cell,
        };
        try self.runCells(cells, &adapter, Adapter.apply);
    }
};

fn ceilDivision(numerator: usize, denominator: usize) !usize {
    return try std.math.add(usize, numerator / denominator, @intFromBool(numerator % denominator != 0));
}

const FillContext = struct { values: []usize };

fn fillSquares(context: *FillContext, range: CellRange) !void {
    for (range.first..range.end) |index| context.values[index] = try std.math.mul(usize, index, index);
}

test "runtime workers process every cell in one loaded tile exactly" {
    const allocator = std.testing.allocator;
    const values = try allocator.alloc(usize, 10_003);
    defer allocator.free(values);
    @memset(values, std.math.maxInt(usize));
    var context = FillContext{ .values = values };
    const executor = try CpuExecutor.init(allocator, 7, 113);
    try executor.run(values.len, &context, fillSquares);
    for (values, 0..) |value, index| try std.testing.expectEqual(index * index, value);
}

const FailureContext = struct { failing_cell: usize };

fn failInsideTile(context: *FailureContext, range: CellRange) !void {
    if (context.failing_cell >= range.first and context.failing_cell < range.end) return error.DeliberateKernelFailure;
}

test "worker kernel errors are propagated" {
    var context = FailureContext{ .failing_cell = 777 };
    const executor = try CpuExecutor.init(std.testing.allocator, 4, 64);
    try std.testing.expectError(error.DeliberateKernelFailure, executor.run(1000, &context, failInsideTile));
}

test "indexed cell dispatch adapts CellRange kernels without touching halos" {
    const values = try std.testing.allocator.alloc(usize, 12);
    defer std.testing.allocator.free(values);
    @memset(values, std.math.maxInt(usize));
    const morton_owned_cells = [_]usize{ 1, 2, 5, 6, 9, 10 };
    var context = FillContext{ .values = values };
    const executor = try CpuExecutor.init(std.testing.allocator, 3, 6);
    try executor.runIndexedCells(
        &morton_owned_cells,
        &context,
        fillSquares,
    );
    for (values, 0..) |value, cell| {
        if (std.mem.indexOfScalar(
            usize,
            &morton_owned_cells,
            cell,
        ) != null) {
            try std.testing.expectEqual(cell * cell, value);
        } else {
            try std.testing.expectEqual(std.math.maxInt(usize), value);
        }
    }
}

const LayerOwnershipContext = struct {
    owner_by_layer: []usize,
    next_owner: std.atomic.Value(usize) = .init(0),
};

fn markWholeCellLayers(
    context: *LayerOwnershipContext,
    range: CellRange,
) !void {
    const owner = context.next_owner.fetchAdd(1, .monotonic);
    for (range.first..range.end) |layer| {
        if (context.owner_by_layer[layer] != std.math.maxInt(usize))
            return error.LayerProcessedMoreThanOnce;
        context.owner_by_layer[layer] = owner;
    }
}

test "cell-layer dispatch never splits one grid column between workers" {
    const layers_per_cell: usize = 7;
    const cell_count: usize = 13;
    const owners = try std.testing.allocator.alloc(
        usize,
        cell_count * layers_per_cell,
    );
    defer std.testing.allocator.free(owners);
    @memset(owners, std.math.maxInt(usize));
    var context = LayerOwnershipContext{
        .owner_by_layer = owners,
    };
    const executor = try CpuExecutor.init(std.testing.allocator, 4, cell_count);
    try executor.runCellLayers(
        cell_count,
        layers_per_cell,
        &context,
        markWholeCellLayers,
    );
    for (0..cell_count) |cell| {
        const first = cell * layers_per_cell;
        for (owners[first..][0..layers_per_cell]) |owner|
            try std.testing.expectEqual(owners[first], owner);
    }
}

const IndexedLayerContext = struct {
    visits_by_layer: []u8,
    visited_cells: []u8,
};

fn markIndexedCellLayers(
    context: *IndexedLayerContext,
    range: CellRange,
) !void {
    const cell = range.first / 3;
    if (context.visited_cells[cell] != 0) return error.CellProcessedMoreThanOnce;
    context.visited_cells[cell] = 1;
    for (range.first..range.end) |layer| {
        if (context.visits_by_layer[layer] != 0)
            return error.LayerProcessedMoreThanOnce;
        context.visits_by_layer[layer] = 1;
    }
}

test "indexed cell-layer dispatch supports non-contiguous Morton tile cells" {
    const cell_count: usize = 12;
    const layers_per_cell: usize = 3;
    // A two-column rectangular tile inside a four-column global domain.
    const morton_owned_cells = [_]usize{ 1, 2, 5, 6, 9, 10 };
    const visits = try std.testing.allocator.alloc(
        u8,
        cell_count * layers_per_cell,
    );
    defer std.testing.allocator.free(visits);
    const visited_cells = try std.testing.allocator.alloc(u8, cell_count);
    defer std.testing.allocator.free(visited_cells);
    @memset(visits, 0);
    @memset(visited_cells, 0);
    var context = IndexedLayerContext{
        .visits_by_layer = visits,
        .visited_cells = visited_cells,
    };
    const executor = try CpuExecutor.init(std.testing.allocator, 3, 6);
    try executor.runIndexedCellLayers(
        &morton_owned_cells,
        layers_per_cell,
        &context,
        markIndexedCellLayers,
    );
    for (0..cell_count) |cell| {
        const expected: u8 =
            if (std.mem.indexOfScalar(usize, &morton_owned_cells, cell) != null)
                1
            else
                0;
        try std.testing.expectEqual(expected, visited_cells[cell]);
        const first = cell * layers_per_cell;
        for (visits[first..][0..layers_per_cell]) |visit|
            try std.testing.expectEqual(expected, visit);
    }
}
