const std = @import("std");

pub const Direction = enum { east_west, north_south, vertical };
pub const Side = enum { forward, reverse };

pub const GridBounds = struct {
    first_column: usize,
    last_column: usize,
    first_row: usize,
    last_row: usize,
};

pub const Coordinate3d = struct { column: isize, row: isize, layer: isize };

pub const BoundaryChargeFractions = struct {
    east: f64,
    west: f64,
    south: f64,
    north: f64,
};

pub const Inputs = struct {
    salt_transport_enabled: bool,
    bounds: GridBounds,
    column: usize,
    row: usize,
    layer: usize,
    bottom_layer: usize,
    direction: Direction,
    side: Side,
    charge_fraction: BoundaryChargeFractions,
};

pub const Boundary = struct {
    current: Coordinate3d,
    destination: Coordinate3d,
    forward_neighbor: Coordinate3d,
    backward_neighbor: Coordinate3d,
    outward_sign: f64,
    charge_fraction: ?f64,
};

/// Compatibility translation of TRNSFRS.F lines 7097--7198.
/// Lines 7077--7085 are control-only CONTINUE labels. Coordinates are signed
/// so zero-based west/north ghost positions remain representable.
pub fn locate(inputs: Inputs) !?Boundary {
    if (!inputs.salt_transport_enabled) return null;
    try validate(inputs);
    const column: isize = @intCast(inputs.column);
    const row: isize = @intCast(inputs.row);
    const layer: isize = @intCast(inputs.layer);
    const current = Coordinate3d{ .column = column, .row = row, .layer = layer };

    return switch (inputs.direction) {
        .east_west => switch (inputs.side) {
            .forward => if (inputs.column == inputs.bounds.last_column)
                Boundary{ .current = current, .destination = .{ .column = column + 1, .row = row, .layer = layer }, .forward_neighbor = .{ .column = column + 1, .row = row, .layer = layer }, .backward_neighbor = .{ .column = column - 1, .row = row, .layer = layer }, .outward_sign = -1, .charge_fraction = inputs.charge_fraction.east }
            else
                null,
            .reverse => if (inputs.column == inputs.bounds.first_column)
                Boundary{ .current = current, .destination = current, .forward_neighbor = .{ .column = column + 1, .row = row, .layer = layer }, .backward_neighbor = .{ .column = column - 1, .row = row, .layer = layer }, .outward_sign = 1, .charge_fraction = inputs.charge_fraction.west }
            else
                null,
        },
        .north_south => switch (inputs.side) {
            .forward => if (inputs.row == inputs.bounds.last_row)
                Boundary{ .current = current, .destination = .{ .column = column, .row = row + 1, .layer = layer }, .forward_neighbor = .{ .column = column, .row = row + 1, .layer = layer }, .backward_neighbor = .{ .column = column, .row = row - 1, .layer = layer }, .outward_sign = -1, .charge_fraction = inputs.charge_fraction.south }
            else
                null,
            .reverse => if (inputs.row == inputs.bounds.first_row)
                Boundary{ .current = current, .destination = current, .forward_neighbor = .{ .column = column, .row = row + 1, .layer = layer }, .backward_neighbor = .{ .column = column, .row = row - 1, .layer = layer }, .outward_sign = 1, .charge_fraction = inputs.charge_fraction.north }
            else
                null,
        },
        .vertical => switch (inputs.side) {
            .forward => if (inputs.layer == inputs.bottom_layer)
                Boundary{ .current = current, .destination = .{ .column = column, .row = row, .layer = layer + 1 }, .forward_neighbor = .{ .column = column, .row = row, .layer = layer + 1 }, .backward_neighbor = current, .outward_sign = -1, .charge_fraction = null }
            else
                null,
            .reverse => null,
        },
    };
}

fn validate(inputs: Inputs) !void {
    if (inputs.bounds.first_column > inputs.bounds.last_column or inputs.bounds.first_row > inputs.bounds.last_row or
        inputs.column < inputs.bounds.first_column or inputs.column > inputs.bounds.last_column or
        inputs.row < inputs.bounds.first_row or inputs.row > inputs.bounds.last_row or inputs.layer > inputs.bottom_layer)
        return error.InvalidExternalSoilBoundaryCoordinate;
    inline for (.{ inputs.charge_fraction.east, inputs.charge_fraction.west, inputs.charge_fraction.south, inputs.charge_fraction.north }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteExternalSoilBoundaryChargeFraction;
}

fn fixture(direction: Direction, side: Side, column: usize, row: usize, layer: usize) Inputs {
    return .{ .salt_transport_enabled = true, .bounds = .{ .first_column = 0, .last_column = 2, .first_row = 0, .last_row = 3 }, .column = column, .row = row, .layer = layer, .bottom_layer = 4, .direction = direction, .side = side, .charge_fraction = .{ .east = 0.1, .west = 0.2, .south = 0.3, .north = 0.4 } };
}

test "TRNSFRS locates east and west external boundaries with legacy destinations" {
    const east = (try locate(fixture(.east_west, .forward, 2, 1, 2))).?;
    const west = (try locate(fixture(.east_west, .reverse, 0, 1, 2))).?;
    try std.testing.expectEqual(@as(isize, 3), east.destination.column);
    try std.testing.expectEqual(@as(f64, -1), east.outward_sign);
    try std.testing.expectEqual(@as(f64, 0.1), east.charge_fraction.?);
    try std.testing.expectEqual(@as(isize, 0), west.destination.column);
    try std.testing.expectEqual(@as(isize, -1), west.backward_neighbor.column);
    try std.testing.expectEqual(@as(f64, 1), west.outward_sign);
}

test "TRNSFRS locates south and north external boundaries" {
    const south = (try locate(fixture(.north_south, .forward, 1, 3, 2))).?;
    const north = (try locate(fixture(.north_south, .reverse, 1, 0, 2))).?;
    try std.testing.expectEqual(@as(isize, 4), south.destination.row);
    try std.testing.expectEqual(@as(f64, 0.3), south.charge_fraction.?);
    try std.testing.expectEqual(@as(isize, 0), north.destination.row);
    try std.testing.expectEqual(@as(f64, 0.4), north.charge_fraction.?);
}

test "TRNSFRS permits only forward bottom vertical boundary" {
    const bottom = (try locate(fixture(.vertical, .forward, 1, 1, 4))).?;
    try std.testing.expectEqual(@as(isize, 5), bottom.destination.layer);
    try std.testing.expect(bottom.charge_fraction == null);
    try std.testing.expect((try locate(fixture(.vertical, .reverse, 1, 1, 4))) == null);
}

test "TRNSFRS skips interior sides and disabled salt transport" {
    try std.testing.expect((try locate(fixture(.east_west, .forward, 1, 1, 2))) == null);
    var disabled = fixture(.east_west, .forward, 2, 1, 2);
    disabled.salt_transport_enabled = false;
    disabled.column = 99;
    disabled.charge_fraction.east = std.math.nan(f64);
    try std.testing.expect((try locate(disabled)) == null);
}

test "invalid runtime coordinate fails before boundary publication" {
    var inputs = fixture(.north_south, .forward, 1, 3, 2);
    inputs.row = 4;
    try std.testing.expectError(error.InvalidExternalSoilBoundaryCoordinate, locate(inputs));
}
