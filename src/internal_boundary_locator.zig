const std = @import("std");

pub const GridCoordinate = struct {
    row: usize,
    column: usize,
};

pub const GridBounds = struct {
    north_row: usize,
    south_row: usize,
    west_column: usize,
    east_column: usize,
};

/// Exact legacy loop order: N=1/NN=1 east, N=1/NN=2 west,
/// N=2/NN=1 south, N=2/NN=2 north.
pub const BoundaryDirection = enum {
    east,
    west,
    south,
    north,
};

pub const InternalBoundary = struct {
    source: GridCoordinate,
    destination: GridCoordinate,
    /// Zero-based translation of legacy axis `N`: east-west=0, north-south=1.
    axis_index: usize,
    /// Zero-based translation of legacy side `NN`: forward=0, backward=1.
    side_index: usize,
};

/// Exact runtime-geometry translation of TRNSFRS.F lines 4016--4044.
/// Returns null for a physical domain edge, matching `GO TO 4305`.
pub fn locate(bounds: GridBounds, source: GridCoordinate, direction: BoundaryDirection) !?InternalBoundary {
    if (bounds.north_row > bounds.south_row or bounds.west_column > bounds.east_column)
        return error.InvalidInternalBoundaryBounds;
    if (source.row < bounds.north_row or source.row > bounds.south_row or
        source.column < bounds.west_column or source.column > bounds.east_column)
        return error.InternalBoundarySourceOutOfBounds;

    return switch (direction) {
        .east => if (source.column == bounds.east_column) null else .{
            .source = source,
            .destination = .{ .row = source.row, .column = source.column + 1 },
            .axis_index = 0,
            .side_index = 0,
        },
        .west => if (source.column == bounds.west_column) null else .{
            .source = source,
            .destination = .{ .row = source.row, .column = source.column - 1 },
            .axis_index = 0,
            .side_index = 1,
        },
        .south => if (source.row == bounds.south_row) null else .{
            .source = source,
            .destination = .{ .row = source.row + 1, .column = source.column },
            .axis_index = 1,
            .side_index = 0,
        },
        .north => if (source.row == bounds.north_row) null else .{
            .source = source,
            .destination = .{ .row = source.row - 1, .column = source.column },
            .axis_index = 1,
            .side_index = 1,
        },
    };
}

test "TRNSFRS locates four interior neighbors in exact N NN order" {
    const bounds = GridBounds{ .north_row = 2, .south_row = 4, .west_column = 5, .east_column = 7 };
    const source = GridCoordinate{ .row = 3, .column = 6 };
    const east = (try locate(bounds, source, .east)).?;
    const west = (try locate(bounds, source, .west)).?;
    const south = (try locate(bounds, source, .south)).?;
    const north = (try locate(bounds, source, .north)).?;
    try std.testing.expectEqual(GridCoordinate{ .row = 3, .column = 7 }, east.destination);
    try std.testing.expectEqual(@as(usize, 0), east.axis_index);
    try std.testing.expectEqual(@as(usize, 0), east.side_index);
    try std.testing.expectEqual(GridCoordinate{ .row = 3, .column = 5 }, west.destination);
    try std.testing.expectEqual(@as(usize, 1), west.side_index);
    try std.testing.expectEqual(GridCoordinate{ .row = 4, .column = 6 }, south.destination);
    try std.testing.expectEqual(@as(usize, 1), south.axis_index);
    try std.testing.expectEqual(GridCoordinate{ .row = 2, .column = 6 }, north.destination);
}

test "physical domain edges reproduce legacy boundary skips" {
    const bounds = GridBounds{ .north_row = 2, .south_row = 4, .west_column = 5, .east_column = 7 };
    try std.testing.expect((try locate(bounds, .{ .row = 2, .column = 5 }, .north)) == null);
    try std.testing.expect((try locate(bounds, .{ .row = 2, .column = 5 }, .west)) == null);
    try std.testing.expect((try locate(bounds, .{ .row = 4, .column = 7 }, .south)) == null);
    try std.testing.expect((try locate(bounds, .{ .row = 4, .column = 7 }, .east)) == null);
}

test "invalid runtime bounds fail before coordinate arithmetic" {
    try std.testing.expectError(
        error.InvalidInternalBoundaryBounds,
        locate(.{ .north_row = 3, .south_row = 2, .west_column = 0, .east_column = 1 }, .{ .row = 2, .column = 0 }, .east),
    );
}

test "source outside runtime domain is rejected" {
    try std.testing.expectError(
        error.InternalBoundarySourceOutOfBounds,
        locate(.{ .north_row = 0, .south_row = 1, .west_column = 0, .east_column = 1 }, .{ .row = 2, .column = 0 }, .north),
    );
}
