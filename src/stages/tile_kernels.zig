//! Kernel drivers that fan work across cells, layers and spatial tiles.
//!
//! Extracted verbatim from `ecosys_ng.zig` so the entry point holds only
//! `main`. Declaration bodies are unchanged.

const std = @import("std");
const ecosys = @import("ecosys_ng");
/// Dispatches cell-local science across either the resident domain or the
/// non-contiguous Morton-owned cells of one loaded tile. Complete vertical
/// columns remain indivisible worker units.
pub fn runScienceCellLayers(
    context: anytype,
    kernel_context: anytype,
    comptime kernel: anytype,
) !void {
    const cells = context.active_tile_cells.* orelse
        return error.MissingActiveScienceTile;
    try context.executor.runIndexedCellLayers(
        cells,
        context.grid.soil_layer_capacity,
        kernel_context,
        kernel,
    );
}

pub fn runScienceCells(
    context: anytype,
    kernel_context: anytype,
    comptime kernel: anytype,
) !void {
    const cells = context.active_tile_cells.* orelse
        return error.MissingActiveScienceTile;
    try context.executor.runIndexedCells(cells, kernel_context, kernel);
}

/// Preserves an hourly stage boundary while executing that stage across
/// serial Morton tiles. No owned-cell list is allocated here: TilePlan builds
/// and validates the lists once during runtime initialization.
pub fn runKernelAcrossSerialTiles(
    context: anytype,
    kernel_context: anytype,
    comptime kernel: anytype,
) !void {
    const plan = context.tile_plan;
    for (plan.tiles, 0..) |_, tile_index| {
        try context.executor.runIndexedCells(
            try plan.ownedCells(tile_index),
            kernel_context,
            kernel,
        );
    }
}

/// Initialization and refresh kernels obey the same execution boundary as the
/// hourly science: tiles advance serially, while owned cells within the active
/// tile may run concurrently.
pub fn runKernelAcrossSerialTilePlan(
    executor: ecosys.compute.CpuExecutor,
    plan: *const ecosys.spatial_grid.TilePlan,
    kernel_context: anytype,
    comptime kernel: anytype,
) !void {
    for (plan.tiles, 0..) |_, tile_index| {
        try executor.runIndexedCells(
            try plan.ownedCells(tile_index),
            kernel_context,
            kernel,
        );
    }
}
