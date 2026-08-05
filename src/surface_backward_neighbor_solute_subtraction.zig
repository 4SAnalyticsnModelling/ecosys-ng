const std = @import("std");

pub const runoff_species_count = 42;
pub const snow_species_count = 41;

pub const BoundarySide = enum {
    forward,
    reverse,
};

pub const BackwardNeighbor = struct {
    row_fortran: i32,
    column_fortran: i32,
};

pub const Inputs = struct {
    side: BoundarySide,
    backward_neighbor: BackwardNeighbor,
    /// Fluxes at `(N5B,N4B)`, in their distinct TRNSFRS species orders.
    runoff_neighbor_amount_per_step: []const f64,
    snow_neighbor_amount_per_step: []const f64,
};

pub const Totals = struct {
    runoff_amount_per_step: []f64,
    snow_amount_per_step: []f64,
};

/// Compatibility translation of TRNSFRS.F lines 8725--8809.
/// Unlike the preceding forward-neighbour block, the one source guard applies
/// to both the 42 runoff species and the 41 snow species.
pub fn subtractBackwardNeighbor(inputs: Inputs, totals: Totals) !bool {
    if (inputs.backward_neighbor.row_fortran <= 0 or
        inputs.backward_neighbor.column_fortran <= 0 or
        inputs.side != .forward) return false;
    try validate(inputs, totals);

    for (0..runoff_species_count) |species| {
        const result = totals.runoff_amount_per_step[species] -
            inputs.runoff_neighbor_amount_per_step[species];
        if (!std.math.isFinite(result)) return error.NonFiniteBackwardNeighborSoluteResult;
    }
    for (0..snow_species_count) |species| {
        const result = totals.snow_amount_per_step[species] -
            inputs.snow_neighbor_amount_per_step[species];
        if (!std.math.isFinite(result)) return error.NonFiniteBackwardNeighborSoluteResult;
    }
    for (0..runoff_species_count) |species|
        totals.runoff_amount_per_step[species] = totals.runoff_amount_per_step[species] -
            inputs.runoff_neighbor_amount_per_step[species];
    for (0..snow_species_count) |species|
        totals.snow_amount_per_step[species] = totals.snow_amount_per_step[species] -
            inputs.snow_neighbor_amount_per_step[species];
    return true;
}

fn validate(inputs: Inputs, totals: Totals) !void {
    if (inputs.runoff_neighbor_amount_per_step.len != runoff_species_count or
        inputs.snow_neighbor_amount_per_step.len != snow_species_count or
        totals.runoff_amount_per_step.len != runoff_species_count or
        totals.snow_amount_per_step.len != snow_species_count)
        return error.BackwardNeighborSoluteDimensionMismatch;
    for (inputs.runoff_neighbor_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteBackwardNeighborSoluteInput;
    for (inputs.snow_neighbor_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteBackwardNeighborSoluteInput;
    for (totals.runoff_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteBackwardNeighborSoluteInput;
    for (totals.snow_amount_per_step) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteBackwardNeighborSoluteInput;
}

fn validInputs(runoff: []const f64, snow: []const f64) Inputs {
    return .{
        .side = .forward,
        .backward_neighbor = .{ .row_fortran = 2, .column_fortran = 3 },
        .runoff_neighbor_amount_per_step = runoff,
        .snow_neighbor_amount_per_step = snow,
    };
}

test "TRNSFRS subtracts both backward-neighbor species families" {
    const runoff_neighbor = [_]f64{2} ** runoff_species_count;
    const snow_neighbor = [_]f64{3} ** snow_species_count;
    var runoff_total = [_]f64{10} ** runoff_species_count;
    var snow_total = [_]f64{10} ** snow_species_count;
    try std.testing.expect(try subtractBackwardNeighbor(validInputs(&runoff_neighbor, &snow_neighbor), .{
        .runoff_amount_per_step = &runoff_total,
        .snow_amount_per_step = &snow_total,
    }));
    try std.testing.expectEqual(@as(f64, 8), runoff_total[41]);
    try std.testing.expectEqual(@as(f64, 7), snow_total[40]);
}

test "reverse side skips backward-neighbor subtraction" {
    const runoff_neighbor = [_]f64{};
    const snow_neighbor = [_]f64{};
    var inputs = validInputs(&runoff_neighbor, &snow_neighbor);
    inputs.side = .reverse;
    var runoff_total = [_]f64{10} ** runoff_species_count;
    var snow_total = [_]f64{10} ** snow_species_count;
    try std.testing.expect(!try subtractBackwardNeighbor(inputs, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total }));
    try std.testing.expectEqual(@as(f64, 10), runoff_total[0]);
    try std.testing.expectEqual(@as(f64, 10), snow_total[0]);
}

test "nonpositive Fortran neighbor coordinate skips both families" {
    const runoff_neighbor = [_]f64{};
    const snow_neighbor = [_]f64{};
    var inputs = validInputs(&runoff_neighbor, &snow_neighbor);
    inputs.backward_neighbor.row_fortran = 0;
    var runoff_total = [_]f64{10} ** runoff_species_count;
    var snow_total = [_]f64{10} ** snow_species_count;
    try std.testing.expect(!try subtractBackwardNeighbor(inputs, .{ .runoff_amount_per_step = &runoff_total, .snow_amount_per_step = &snow_total }));
    try std.testing.expectEqual(@as(f64, 10), snow_total[0]);
}

test "backward-neighbor dimension mismatch fails atomically" {
    const short_runoff = [_]f64{2} ** (runoff_species_count - 1);
    const snow_neighbor = [_]f64{3} ** snow_species_count;
    var runoff_total = [_]f64{10} ** runoff_species_count;
    var snow_total = [_]f64{10} ** snow_species_count;
    try std.testing.expectError(error.BackwardNeighborSoluteDimensionMismatch, subtractBackwardNeighbor(validInputs(&short_runoff, &snow_neighbor), .{
        .runoff_amount_per_step = &runoff_total,
        .snow_amount_per_step = &snow_total,
    }));
    try std.testing.expectEqual(@as(f64, 10), runoff_total[0]);
}

test "late nonfinite snow result leaves runoff atomic" {
    const runoff_neighbor = [_]f64{2} ** runoff_species_count;
    var snow_neighbor = [_]f64{3} ** snow_species_count;
    snow_neighbor[snow_species_count - 1] = -std.math.floatMax(f64);
    var runoff_total = [_]f64{10} ** runoff_species_count;
    var snow_total = [_]f64{10} ** snow_species_count;
    snow_total[snow_species_count - 1] = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteBackwardNeighborSoluteResult, subtractBackwardNeighbor(validInputs(&runoff_neighbor, &snow_neighbor), .{
        .runoff_amount_per_step = &runoff_total,
        .snow_amount_per_step = &snow_total,
    }));
    try std.testing.expectEqual(@as(f64, 10), runoff_total[0]);
}
