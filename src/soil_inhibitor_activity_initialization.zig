const std = @import("std");

pub const Activity = struct {
    current_urease_inhibition_fraction: f64 = 0,
    initial_urease_inhibition_fraction: f64 = 0,
    current_nitrification_inhibition_fraction: f64 = 0,
    initial_nitrification_inhibition_fraction: f64 = 0,
};

pub const Dimensions = struct {
    column_count: usize,
    row_count: usize,
    soil_layer_count_including_surface: usize,
};

/// Clears the initial inhibitor activity carried by every soil layer.
///
/// Traceability: STARTS.F lines 1636--1639 (`ZNHUI`, `ZNHU0`, `ZNFNI`,
/// and `ZNFN0`). Storage is cell-major, then layer, including surface layer
/// zero. All fields are dimensionless fractions.
pub fn initialize(dimensions: Dimensions, activity_by_cell_layer: []Activity) !void {
    const expected_count = try checkedValueCount(dimensions);
    if (activity_by_cell_layer.len != expected_count)
        return error.InvalidSoilInhibitorActivityDimensions;

    @memset(activity_by_cell_layer, Activity{});
}

fn checkedValueCount(dimensions: Dimensions) !usize {
    if (dimensions.column_count == 0 or
        dimensions.row_count == 0 or
        dimensions.soil_layer_count_including_surface == 0)
        return error.InvalidSoilInhibitorActivityDimensions;

    const cell_count = std.math.mul(
        usize,
        dimensions.column_count,
        dimensions.row_count,
    ) catch return error.InvalidSoilInhibitorActivityDimensions;
    return std.math.mul(
        usize,
        cell_count,
        dimensions.soil_layer_count_including_surface,
    ) catch return error.InvalidSoilInhibitorActivityDimensions;
}

test "initialize clears inhibitor fractions for runtime grid and layers" {
    var activity = [_]Activity{.{
        .current_urease_inhibition_fraction = 0.1,
        .initial_urease_inhibition_fraction = 0.2,
        .current_nitrification_inhibition_fraction = 0.3,
        .initial_nitrification_inhibition_fraction = 0.4,
    }} ** 8;

    try initialize(.{
        .column_count = 2,
        .row_count = 1,
        .soil_layer_count_including_surface = 4,
    }, &activity);

    for (activity) |layer|
        try std.testing.expectEqualDeep(Activity{}, layer);
}

test "initialize rejects mismatched allocation without mutation" {
    var activity = [_]Activity{.{
        .current_urease_inhibition_fraction = 0.75,
    }} ** 3;

    try std.testing.expectError(
        error.InvalidSoilInhibitorActivityDimensions,
        initialize(.{
            .column_count = 2,
            .row_count = 1,
            .soil_layer_count_including_surface = 2,
        }, &activity),
    );
    for (activity) |layer|
        try std.testing.expectEqual(
            @as(f64, 0.75),
            layer.current_urease_inhibition_fraction,
        );
}

test "initialize rejects zero and overflowing dimensions" {
    var none: [0]Activity = .{};
    try std.testing.expectError(
        error.InvalidSoilInhibitorActivityDimensions,
        initialize(.{
            .column_count = 1,
            .row_count = 0,
            .soil_layer_count_including_surface = 1,
        }, &none),
    );
    try std.testing.expectError(
        error.InvalidSoilInhibitorActivityDimensions,
        initialize(.{
            .column_count = std.math.maxInt(usize),
            .row_count = 2,
            .soil_layer_count_including_surface = 1,
        }, &none),
    );
}
