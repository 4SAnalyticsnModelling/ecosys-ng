const std = @import("std");

pub const Coordinate3d = struct { row: usize, column: usize, layer: usize };
pub const GridBounds = struct { north_row: usize, south_row: usize, west_column: usize, east_column: usize };
pub const Direction = enum(u2) { east, south, downward };

pub const ScanState = struct {
    /// Named replacement for legacy `IFLGB`.
    gas_boundary_claimed: bool,
};

pub const Inputs = struct {
    bounds: GridBounds,
    source: Coordinate3d,
    source_layer_count: usize,
    /// Zero-based translation of `NCN`; directions before this are skipped.
    first_direction: Direction,
    source_layer_thickness_m: []const f64,
    minimum_transport_layer_thickness_m: f64,
};

/// Exact TRNSFRS.F line 4847 scan initialization.
pub fn beginScan(state: *ScanState) void {
    state.gas_boundary_claimed = false;
}

/// Exact runtime translation of TRNSFRS.F lines 4848--4889.
/// Returns null for directions excluded by `NCN` or at east/south/bottom edges.
/// If no deeper layer exceeds `DLYRM`, the immediate lower layer is retained,
/// matching the legacy loop's unchanged `N6` behavior.
pub fn locate(inputs: Inputs, direction: Direction) !?Coordinate3d {
    try validate(inputs);
    if (@intFromEnum(direction) < @intFromEnum(inputs.first_direction)) return null;
    return switch (direction) {
        .east => if (inputs.source.column == inputs.bounds.east_column) null else .{
            .row = inputs.source.row,
            .column = inputs.source.column + 1,
            .layer = inputs.source.layer,
        },
        .south => if (inputs.source.row == inputs.bounds.south_row) null else .{
            .row = inputs.source.row + 1,
            .column = inputs.source.column,
            .layer = inputs.source.layer,
        },
        .downward => locateLowerLayer(inputs),
    };
}

fn locateLowerLayer(inputs: Inputs) ?Coordinate3d {
    if (inputs.source.layer + 1 >= inputs.source_layer_count) return null;
    var destination_layer = inputs.source.layer + 1;
    for (inputs.source.layer + 1..inputs.source_layer_count) |layer| {
        if (inputs.source_layer_thickness_m[layer] > inputs.minimum_transport_layer_thickness_m) {
            destination_layer = layer;
            break;
        }
    }
    return .{ .row = inputs.source.row, .column = inputs.source.column, .layer = destination_layer };
}

fn validate(inputs: Inputs) !void {
    if (inputs.bounds.north_row > inputs.bounds.south_row or
        inputs.bounds.west_column > inputs.bounds.east_column)
        return error.InvalidSoilInternalBoundaryBounds;
    if (inputs.source.row < inputs.bounds.north_row or inputs.source.row > inputs.bounds.south_row or
        inputs.source.column < inputs.bounds.west_column or inputs.source.column > inputs.bounds.east_column or
        inputs.source.layer >= inputs.source_layer_count)
        return error.SoilInternalBoundarySourceOutOfBounds;
    if (inputs.source_layer_thickness_m.len != inputs.source_layer_count)
        return error.SoilInternalBoundaryLayerDimensionMismatch;
    if (!std.math.isFinite(inputs.minimum_transport_layer_thickness_m) or
        inputs.minimum_transport_layer_thickness_m < 0)
        return error.InvalidSoilInternalBoundaryThickness;
    for (inputs.source_layer_thickness_m) |thickness_m|
        if (!std.math.isFinite(thickness_m) or thickness_m < 0)
            return error.InvalidSoilInternalBoundaryThickness;
}

fn fixture(thickness: []const f64) Inputs {
    return .{ .bounds = .{ .north_row = 2, .south_row = 4, .west_column = 5, .east_column = 7 }, .source = .{ .row = 3, .column = 6, .layer = 0 }, .source_layer_count = thickness.len, .first_direction = .east, .source_layer_thickness_m = thickness, .minimum_transport_layer_thickness_m = 0.1 };
}

test "TRNSFRS resets IFLGB before runtime layer scan" {
    var state = ScanState{ .gas_boundary_claimed = true };
    beginScan(&state);
    try std.testing.expect(!state.gas_boundary_claimed);
}

test "TRNSFRS locates east south and first sufficiently thick lower layer" {
    const thickness = [_]f64{ 1, 0.01, 0.2, 0.3 };
    const inputs = fixture(&thickness);
    try std.testing.expectEqual(Coordinate3d{ .row = 3, .column = 7, .layer = 0 }, (try locate(inputs, .east)).?);
    try std.testing.expectEqual(Coordinate3d{ .row = 4, .column = 6, .layer = 0 }, (try locate(inputs, .south)).?);
    try std.testing.expectEqual(Coordinate3d{ .row = 3, .column = 6, .layer = 2 }, (try locate(inputs, .downward)).?);
}

test "vertical scan retains immediate lower layer when none exceed DLYRM" {
    const thickness = [_]f64{ 1, 0.01, 0.02 };
    const inputs = fixture(&thickness);
    try std.testing.expectEqual(@as(usize, 1), (try locate(inputs, .downward)).?.layer);
}

test "NCN and physical boundaries produce zero-trip direction skips" {
    const thickness = [_]f64{1};
    var inputs = fixture(&thickness);
    inputs.first_direction = .south;
    try std.testing.expect((try locate(inputs, .east)) == null);
    inputs.source.row = inputs.bounds.south_row;
    try std.testing.expect((try locate(inputs, .south)) == null);
    try std.testing.expect((try locate(inputs, .downward)) == null);
}

test "invalid runtime layer dimensions fail" {
    const thickness = [_]f64{1};
    var inputs = fixture(&thickness);
    inputs.source_layer_count = 2;
    try std.testing.expectError(error.SoilInternalBoundaryLayerDimensionMismatch, locate(inputs, .east));
}
