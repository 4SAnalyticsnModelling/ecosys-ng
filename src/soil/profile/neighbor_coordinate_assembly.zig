const std = @import("std");

pub const Dimensions = struct {
    lon_count: usize,
    lat_count: usize,
    soil_layer_capacity: usize,
};

/// Zero-based logical position, independent of row-major or Morton storage.
pub const Coordinate = struct {
    column: usize,
    row: usize,
    layer: usize,
};

/// Source traversal order from REDIST's `N=1,2,3`.
pub const Axis = enum(u8) {
    column,
    row,
    vertical,

    pub const count: usize = @typeInfo(Axis).@"enum".fields.len;
};

pub const AxisMapping = struct {
    axis: Axis,
    source: Coordinate,
    /// REDIST `N4,N5,N6`: +column, +row, or +layer.
    positive_destination: ?Coordinate,
    /// REDIST `N4B,N5B`, defined only for horizontal axes.
    opposite_destination: ?Coordinate,
};

/// Assembles source and destination coordinates for one soil cell-layer.
///
/// Traceability: REDIST.F lines 1933--1953 (`N1`--`N6`, `N4B`, and `N5B`).
/// Output order is column, row, vertical. Null marks a physical allocation
/// boundary and prevents the source's zero/upper-halo coordinates from being
/// used as array indexes. The vertical opposite coordinate is always null
/// because the source does not assign it in the `N=3` branch.
pub fn assemble(
    dimensions: Dimensions,
    source: Coordinate,
    mappings: []AxisMapping,
) !void {
    try validate(dimensions, source, mappings.len);

    mappings[0] = .{
        .axis = .column,
        .source = source,
        .positive_destination = if (source.column + 1 < dimensions.lon_count)
            .{
                .column = source.column + 1,
                .row = source.row,
                .layer = source.layer,
            }
        else
            null,
        .opposite_destination = if (source.column > 0)
            .{
                .column = source.column - 1,
                .row = source.row,
                .layer = source.layer,
            }
        else
            null,
    };
    mappings[1] = .{
        .axis = .row,
        .source = source,
        .positive_destination = if (source.row + 1 < dimensions.lat_count)
            .{
                .column = source.column,
                .row = source.row + 1,
                .layer = source.layer,
            }
        else
            null,
        .opposite_destination = if (source.row > 0)
            .{
                .column = source.column,
                .row = source.row - 1,
                .layer = source.layer,
            }
        else
            null,
    };
    mappings[2] = .{
        .axis = .vertical,
        .source = source,
        .positive_destination = if (source.layer + 1 < dimensions.soil_layer_capacity)
            .{
                .column = source.column,
                .row = source.row,
                .layer = source.layer + 1,
            }
        else
            null,
        .opposite_destination = null,
    };
}

fn validate(dimensions: Dimensions, source: Coordinate, output_count: usize) !void {
    if (dimensions.lon_count == 0 or
        dimensions.lat_count == 0 or
        dimensions.soil_layer_capacity == 0)
    {
        return error.InvalidSoilNeighborDimensions;
    }
    const cell_count = std.math.mul(
        usize,
        dimensions.lon_count,
        dimensions.lat_count,
    ) catch return error.SoilNeighborDimensionOverflow;
    _ = std.math.mul(
        usize,
        cell_count,
        dimensions.soil_layer_capacity,
    ) catch return error.SoilNeighborDimensionOverflow;
    if (source.column >= dimensions.lon_count or
        source.row >= dimensions.lat_count or
        source.layer >= dimensions.soil_layer_capacity)
    {
        return error.SoilNeighborSourceOutOfBounds;
    }
    if (output_count != Axis.count) return error.SoilNeighborOutputLengthMismatch;
}

test "interior mapping preserves exact REDIST axis and destination order" {
    const dimensions = Dimensions{
        .lon_count = 4,
        .lat_count = 3,
        .soil_layer_capacity = 2,
    };
    const source = Coordinate{ .column = 1, .row = 1, .layer = 0 };
    var mappings: [Axis.count]AxisMapping = undefined;

    try assemble(dimensions, source, &mappings);

    try std.testing.expectEqual(Axis.column, mappings[0].axis);
    try std.testing.expectEqual(source, mappings[0].source);
    try std.testing.expectEqual(
        @as(?Coordinate, .{ .column = 2, .row = 1, .layer = 0 }),
        mappings[0].positive_destination,
    );
    try std.testing.expectEqual(
        @as(?Coordinate, .{ .column = 0, .row = 1, .layer = 0 }),
        mappings[0].opposite_destination,
    );
    try std.testing.expectEqual(Axis.row, mappings[1].axis);
    try std.testing.expectEqual(
        @as(?Coordinate, .{ .column = 1, .row = 2, .layer = 0 }),
        mappings[1].positive_destination,
    );
    try std.testing.expectEqual(
        @as(?Coordinate, .{ .column = 1, .row = 0, .layer = 0 }),
        mappings[1].opposite_destination,
    );
    try std.testing.expectEqual(Axis.vertical, mappings[2].axis);
    try std.testing.expectEqual(
        @as(?Coordinate, .{ .column = 1, .row = 1, .layer = 1 }),
        mappings[2].positive_destination,
    );
    try std.testing.expectEqual(@as(?Coordinate, null), mappings[2].opposite_destination);
}

test "domain edges are explicit and never wrap or underflow" {
    const dimensions = Dimensions{
        .lon_count = 1,
        .lat_count = 1,
        .soil_layer_capacity = 1,
    };
    var mappings: [Axis.count]AxisMapping = undefined;

    try assemble(dimensions, .{ .column = 0, .row = 0, .layer = 0 }, &mappings);

    for (mappings) |mapping| {
        try std.testing.expectEqual(@as(?Coordinate, null), mapping.positive_destination);
        try std.testing.expectEqual(@as(?Coordinate, null), mapping.opposite_destination);
    }
}

test "horizontal adjacency is reciprocal throughout runtime grid" {
    const dimensions = Dimensions{
        .lon_count = 5,
        .lat_count = 4,
        .soil_layer_capacity = 3,
    };
    var mappings: [Axis.count]AxisMapping = undefined;
    var neighbor_mappings: [Axis.count]AxisMapping = undefined;

    for (0..dimensions.soil_layer_capacity) |layer| {
        for (0..dimensions.lat_count) |row| {
            for (0..dimensions.lon_count) |column| {
                const source = Coordinate{ .column = column, .row = row, .layer = layer };
                try assemble(dimensions, source, &mappings);
                inline for (.{ Axis.column, Axis.row }) |axis| {
                    const axis_index = @intFromEnum(axis);
                    if (mappings[axis_index].positive_destination) |positive| {
                        try assemble(dimensions, positive, &neighbor_mappings);
                        try std.testing.expectEqual(
                            @as(?Coordinate, source),
                            neighbor_mappings[axis_index].opposite_destination,
                        );
                    }
                    if (mappings[axis_index].opposite_destination) |opposite| {
                        try assemble(dimensions, opposite, &neighbor_mappings);
                        try std.testing.expectEqual(
                            @as(?Coordinate, source),
                            neighbor_mappings[axis_index].positive_destination,
                        );
                    }
                }
            }
        }
    }
}

test "invalid source and output length fail before output mutation" {
    const dimensions = Dimensions{
        .lon_count = 2,
        .lat_count = 2,
        .soil_layer_capacity = 2,
    };
    const sentinel = AxisMapping{
        .axis = .vertical,
        .source = .{ .column = 9, .row = 9, .layer = 9 },
        .positive_destination = null,
        .opposite_destination = null,
    };
    var mappings = [_]AxisMapping{sentinel} ** Axis.count;
    try std.testing.expectError(
        error.SoilNeighborSourceOutOfBounds,
        assemble(dimensions, .{ .column = 2, .row = 0, .layer = 0 }, &mappings),
    );
    for (mappings) |mapping| try std.testing.expectEqual(sentinel, mapping);

    try std.testing.expectError(
        error.SoilNeighborOutputLengthMismatch,
        assemble(dimensions, .{ .column = 0, .row = 0, .layer = 0 }, mappings[0..2]),
    );
    for (mappings) |mapping| try std.testing.expectEqual(sentinel, mapping);
}

test "overflowing runtime dimensions fail explicitly" {
    const dimensions = Dimensions{
        .lon_count = std.math.maxInt(usize),
        .lat_count = 2,
        .soil_layer_capacity = 1,
    };
    var mappings: [Axis.count]AxisMapping = undefined;
    try std.testing.expectError(
        error.SoilNeighborDimensionOverflow,
        assemble(dimensions, .{ .column = 0, .row = 0, .layer = 0 }, &mappings),
    );
}
