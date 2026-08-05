const std = @import("std");

/// Per-timestep plant litterfall fluxes from extract.f (g step-1).
pub const LitterfallFluxes = struct {
    /// ZCSNC. Total plant C litterfall (g C).
    carbon_g: f64,
    /// ZZSNC. Total plant N litterfall (g N).
    nitrogen_g: f64,
    /// ZPSNC. Total plant P litterfall (g P).
    phosphorus_g: f64,
};

/// Cell-level cumulative litterfall accumulators (g).
pub const CellAccumulators = struct {
    /// UXCSN. Cumulative cell C litterfall.
    carbon_g: f64,
    /// UXZSN. Cumulative cell N litterfall.
    nitrogen_g: f64,
    /// UXPSN. Cumulative cell P litterfall.
    phosphorus_g: f64,
};

/// Grid-level running sums (g), updated once per cell per timestep.
pub const GridTotals = struct {
    /// XCSN. Domain C litterfall.
    carbon_g: f64,
    /// XZSN. Domain N litterfall.
    nitrogen_g: f64,
    /// XPSN. Domain P litterfall.
    phosphorus_g: f64,
};

pub const Result = struct {
    cell: CellAccumulators,
    grid: GridTotals,
};

/// Direct translation of REDIST lines 4682--4687.
///
/// Accumulates per-cell and domain C, N, and P litterfall from EXTRACT outputs.
pub fn accumulate(
    fluxes: LitterfallFluxes,
    cell: CellAccumulators,
    grid: GridTotals,
) !Result {
    inline for (@typeInfo(LitterfallFluxes).@"struct".fields) |field|
        if (!std.math.isFinite(@field(fluxes, field.name)))
            return error.InvalidLitterfallFlux;
    inline for (@typeInfo(CellAccumulators).@"struct".fields) |field|
        if (!std.math.isFinite(@field(cell, field.name)))
            return error.InvalidLitterfallAccumulator;
    inline for (@typeInfo(GridTotals).@"struct".fields) |field|
        if (!std.math.isFinite(@field(grid, field.name)))
            return error.InvalidLitterfallAccumulator;

    const result = Result{
        .cell = .{
            .carbon_g = cell.carbon_g + fluxes.carbon_g,
            .nitrogen_g = cell.nitrogen_g + fluxes.nitrogen_g,
            .phosphorus_g = cell.phosphorus_g + fluxes.phosphorus_g,
        },
        .grid = .{
            .carbon_g = grid.carbon_g + fluxes.carbon_g,
            .nitrogen_g = grid.nitrogen_g + fluxes.nitrogen_g,
            .phosphorus_g = grid.phosphorus_g + fluxes.phosphorus_g,
        },
    };
    inline for (@typeInfo(CellAccumulators).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.cell, field.name)))
            return error.NonFiniteLitterfallAccumulator;
    inline for (@typeInfo(GridTotals).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.grid, field.name)))
            return error.NonFiniteLitterfallAccumulator;
    return result;
}

test "REDIST litterfall C accumulates into both cell and grid" {
    const result = try accumulate(
        .{ .carbon_g = 1.5, .nitrogen_g = 0.0, .phosphorus_g = 0.0 },
        .{ .carbon_g = 10.0, .nitrogen_g = 0.0, .phosphorus_g = 0.0 },
        .{ .carbon_g = 100.0, .nitrogen_g = 0.0, .phosphorus_g = 0.0 },
    );
    try std.testing.expectEqual(@as(f64, 11.5), result.cell.carbon_g);
    try std.testing.expectEqual(@as(f64, 101.5), result.grid.carbon_g);
}

test "REDIST litterfall accumulates C, N, and P independently" {
    const result = try accumulate(
        .{ .carbon_g = 2.0, .nitrogen_g = 0.1, .phosphorus_g = 0.01 },
        .{ .carbon_g = 0.0, .nitrogen_g = 0.0, .phosphorus_g = 0.0 },
        .{ .carbon_g = 0.0, .nitrogen_g = 0.0, .phosphorus_g = 0.0 },
    );
    try std.testing.expectEqual(@as(f64, 2.0), result.cell.carbon_g);
    try std.testing.expectEqual(@as(f64, 0.1), result.cell.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 0.01), result.cell.phosphorus_g);
}

test "REDIST litterfall preserves prior cell state" {
    const result = try accumulate(
        .{ .carbon_g = 0.0, .nitrogen_g = 0.0, .phosphorus_g = 0.0 },
        .{ .carbon_g = 5.0, .nitrogen_g = 3.0, .phosphorus_g = 0.5 },
        .{ .carbon_g = 50.0, .nitrogen_g = 30.0, .phosphorus_g = 5.0 },
    );
    try std.testing.expectEqual(@as(f64, 5.0), result.cell.carbon_g);
    try std.testing.expectEqual(@as(f64, 50.0), result.grid.carbon_g);
}

test "REDIST litterfall rejects non-finite flux" {
    try std.testing.expectError(
        error.InvalidLitterfallFlux,
        accumulate(
            .{ .carbon_g = std.math.nan(f64), .nitrogen_g = 0.0, .phosphorus_g = 0.0 },
            .{ .carbon_g = 0.0, .nitrogen_g = 0.0, .phosphorus_g = 0.0 },
            .{ .carbon_g = 0.0, .nitrogen_g = 0.0, .phosphorus_g = 0.0 },
        ),
    );
}

test "REDIST litterfall rejects domain total overflow" {
    try std.testing.expectError(
        error.NonFiniteLitterfallAccumulator,
        accumulate(
            .{ .carbon_g = std.math.floatMax(f64), .nitrogen_g = 0.0, .phosphorus_g = 0.0 },
            .{ .carbon_g = 0.0, .nitrogen_g = 0.0, .phosphorus_g = 0.0 },
            .{ .carbon_g = std.math.floatMax(f64), .nitrogen_g = 0.0, .phosphorus_g = 0.0 },
        ),
    );
}
