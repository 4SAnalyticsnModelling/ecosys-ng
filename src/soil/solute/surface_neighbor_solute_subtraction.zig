const std = @import("std");

pub const runoff_species_count = 42;
pub const snow_species_count = 41;

pub const BoundaryState = enum {
    connected,
    blocked,
};

pub const Inputs = struct {
    runoff_boundary: BoundaryState,
    snow_boundary: BoundaryState,
    /// Fluxes at the forward neighbour `(N5,N4)`, in TRNSFRS species order.
    runoff_neighbor_amount_per_step: []const f64,
    snow_neighbor_amount_per_step: []const f64,
};

pub const Totals = struct {
    runoff_amount_per_step: []f64,
    snow_amount_per_step: []f64,
};

/// Compatibility translation of TRNSFRS.F lines 8638--8724.
/// The source independently gates runoff with IFLBM and snow with IFLBMS;
/// an IFLB value of zero is represented by `.connected`.
pub fn subtractForwardNeighbor(inputs: Inputs, totals: Totals) !bool {
    const subtract_runoff = inputs.runoff_boundary == .connected;
    const subtract_snow = inputs.snow_boundary == .connected;
    if (!subtract_runoff and !subtract_snow) return false;

    if (subtract_runoff) try validateRunoff(inputs, totals);
    if (subtract_snow) try validateSnow(inputs, totals);

    if (subtract_runoff) {
        for (0..runoff_species_count) |species| {
            const result = totals.runoff_amount_per_step[species] -
                inputs.runoff_neighbor_amount_per_step[species];
            if (!std.math.isFinite(result)) return error.NonFiniteNeighborSoluteResult;
        }
    }
    if (subtract_snow) {
        for (0..snow_species_count) |species| {
            const result = totals.snow_amount_per_step[species] -
                inputs.snow_neighbor_amount_per_step[species];
            if (!std.math.isFinite(result)) return error.NonFiniteNeighborSoluteResult;
        }
    }

    if (subtract_runoff) {
        for (0..runoff_species_count) |species|
            totals.runoff_amount_per_step[species] = totals.runoff_amount_per_step[species] -
                inputs.runoff_neighbor_amount_per_step[species];
    }
    if (subtract_snow) {
        for (0..snow_species_count) |species|
            totals.snow_amount_per_step[species] = totals.snow_amount_per_step[species] -
                inputs.snow_neighbor_amount_per_step[species];
    }
    return true;
}

fn validateRunoff(inputs: Inputs, totals: Totals) !void {
    if (inputs.runoff_neighbor_amount_per_step.len != runoff_species_count or
        totals.runoff_amount_per_step.len != runoff_species_count)
        return error.NeighborSoluteDimensionMismatch;
    for (inputs.runoff_neighbor_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteNeighborSoluteInput;
    for (totals.runoff_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteNeighborSoluteInput;
}

fn validateSnow(inputs: Inputs, totals: Totals) !void {
    if (inputs.snow_neighbor_amount_per_step.len != snow_species_count or
        totals.snow_amount_per_step.len != snow_species_count)
        return error.NeighborSoluteDimensionMismatch;
    for (inputs.snow_neighbor_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteNeighborSoluteInput;
    for (totals.snow_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteNeighborSoluteInput;
}

test "TRNSFRS subtracts both connected forward-neighbor families" {
    const runoff_neighbor = [_]f64{2} ** runoff_species_count;
    const snow_neighbor = [_]f64{3} ** snow_species_count;
    var runoff_total = [_]f64{10} ** runoff_species_count;
    var snow_total = [_]f64{10} ** snow_species_count;
    try std.testing.expect(try subtractForwardNeighbor(.{
        .runoff_boundary = .connected,
        .snow_boundary = .connected,
        .runoff_neighbor_amount_per_step = &runoff_neighbor,
        .snow_neighbor_amount_per_step = &snow_neighbor,
    }, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total }));
    try std.testing.expectEqual(@as(f64, 8), runoff_total[41]);
    try std.testing.expectEqual(@as(f64, 7), snow_total[40]);
}

test "TRNSFRS independently gates runoff and snow neighbor subtraction" {
    const runoff_neighbor = [_]f64{};
    const snow_neighbor = [_]f64{3} ** snow_species_count;
    var runoff_total = [_]f64{10} ** runoff_species_count;
    var snow_total = [_]f64{10} ** snow_species_count;
    _ = try subtractForwardNeighbor(.{
        .runoff_boundary = .blocked,
        .snow_boundary = .connected,
        .runoff_neighbor_amount_per_step = &runoff_neighbor,
        .snow_neighbor_amount_per_step = &snow_neighbor,
    }, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total });
    try std.testing.expectEqual(@as(f64, 10), runoff_total[0]);
    try std.testing.expectEqual(@as(f64, 7), snow_total[0]);
}

test "blocked forward-neighbor families are a no-op" {
    const runoff_neighbor = [_]f64{};
    const snow_neighbor = [_]f64{};
    var runoff_total = [_]f64{};
    var snow_total = [_]f64{};
    try std.testing.expect(!try subtractForwardNeighbor(.{
        .runoff_boundary = .blocked,
        .snow_boundary = .blocked,
        .runoff_neighbor_amount_per_step = &runoff_neighbor,
        .snow_neighbor_amount_per_step = &snow_neighbor,
    }, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total }));
}

test "neighbor topology mismatch fails atomically" {
    const short_runoff = [_]f64{2} ** (runoff_species_count - 1);
    const snow_neighbor = [_]f64{3} ** snow_species_count;
    var runoff_total = [_]f64{10} ** runoff_species_count;
    var snow_total = [_]f64{10} ** snow_species_count;
    try std.testing.expectError(error.NeighborSoluteDimensionMismatch, subtractForwardNeighbor(.{
        .runoff_boundary = .connected,
        .snow_boundary = .connected,
        .runoff_neighbor_amount_per_step = &short_runoff,
        .snow_neighbor_amount_per_step = &snow_neighbor,
    }, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total }));
    try std.testing.expectEqual(@as(f64, 10), runoff_total[0]);
}

test "late invalid snow neighbor leaves runoff atomic" {
    const runoff_neighbor = [_]f64{2} ** runoff_species_count;
    var snow_neighbor = [_]f64{3} ** snow_species_count;
    snow_neighbor[snow_species_count - 1] = std.math.nan(f64);
    var runoff_total = [_]f64{10} ** runoff_species_count;
    var snow_total = [_]f64{10} ** snow_species_count;
    try std.testing.expectError(error.NonFiniteNeighborSoluteInput, subtractForwardNeighbor(.{
        .runoff_boundary = .connected,
        .snow_boundary = .connected,
        .runoff_neighbor_amount_per_step = &runoff_neighbor,
        .snow_neighbor_amount_per_step = &snow_neighbor,
    }, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total }));
    try std.testing.expectEqual(@as(f64, 10), runoff_total[0]);
}
