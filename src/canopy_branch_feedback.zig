const std = @import("std");
const stomatal = @import("canopy_stomatal_resistance.zig");

pub const Inputs = struct {
    phenology_type: u8,
    growth_habit: u8,
    aboveground_turnover_type: u8,
    leafout_accumulated_h: f64,
    leafout_required_h: f64,
    leafoff_accumulated_h: f64,
    leafoff_required_h: f64,
    nonstructural_c_g_per_g_c: f64,
    nonstructural_n_g_per_g_c: f64,
    nonstructural_p_g_per_g_c: f64,
    heat_stress_h: f64,
    dehardening_h: f64,
    hours_without_grain_fill_h: f64,
    annual_termination_hours_without_grain_fill_h: f64,
};

pub const Result = struct {
    c3_feedback_fraction: f64,
    annual_termination_fraction: f64,
    photosynthetically_active: bool,
};

/// STOMATE.F 110--189 branch topology and FDBK/FDBKX calculation, with the
/// inactive source values at lines 631--633 mapped explicitly.
pub fn compute(inputs: Inputs) !Result {
    const feedback = try stomatal.branchFeedback(
        inputs.phenology_type,
        inputs.growth_habit,
        inputs.aboveground_turnover_type,
        inputs.leafout_accumulated_h,
        inputs.leafout_required_h,
        inputs.leafoff_accumulated_h,
        inputs.leafoff_required_h,
        inputs.nonstructural_c_g_per_g_c,
        inputs.nonstructural_n_g_per_g_c,
        inputs.nonstructural_p_g_per_g_c,
        inputs.heat_stress_h,
        inputs.dehardening_h,
        inputs.hours_without_grain_fill_h,
        inputs.annual_termination_hours_without_grain_fill_h,
    );
    return .{
        .c3_feedback_fraction = feedback.c3_fraction,
        .annual_termination_fraction = feedback.annual_termination_fraction,
        .photosynthetically_active = feedback.photosynthetically_active,
    };
}

/// `branch_offsets` assigns one contiguous runtime branch range to every
/// caller-defined cell/species entry.
pub fn computeRuntimeTopology(
    branch_offsets: []const usize,
    inputs: []const Inputs,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len or
        branch_offsets.len == 0)
        return error.CanopyBranchFeedbackDimensionMismatch;
    if (branch_offsets[0] != 0 or
        branch_offsets[branch_offsets.len - 1] != inputs.len)
        return error.InvalidCanopyBranchFeedbackOffsets;
    for (0..branch_offsets.len - 1) |owner|
        if (branch_offsets[owner] > branch_offsets[owner + 1])
            return error.InvalidCanopyBranchFeedbackOffsets;
    for (inputs, scratch) |branch_inputs, *candidate|
        candidate.* = try compute(branch_inputs);
    @memcpy(destination, scratch);
}

fn activeInputs() Inputs {
    return .{
        .phenology_type = 1,
        .growth_habit = 0,
        .aboveground_turnover_type = 2,
        .leafout_accumulated_h = 300,
        .leafout_required_h = 200,
        .leafoff_accumulated_h = 0,
        .leafoff_required_h = 100,
        .nonstructural_c_g_per_g_c = 10,
        .nonstructural_n_g_per_g_c = 0.5,
        .nonstructural_p_g_per_g_c = 0.05,
        .heat_stress_h = 1,
        .dehardening_h = 125,
        .hours_without_grain_fill_h = 168,
        .annual_termination_hours_without_grain_fill_h = 336,
    };
}

test "active branch maps exact source FDBK and FDBKX fields" {
    const result = try compute(activeInputs());
    try std.testing.expect(result.photosynthetically_active);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0 / 12.0), result.c3_feedback_fraction, 1e-15);
    try std.testing.expectEqual(@as(f64, 0.5), result.annual_termination_fraction);
}

test "inactive gate maps source FDBK zero and FDBKX one" {
    var inputs = activeInputs();
    inputs.growth_habit = 1;
    inputs.leafout_accumulated_h = 0;
    inputs.leafoff_accumulated_h = inputs.leafoff_required_h;
    const result = try compute(inputs);
    try std.testing.expect(!result.photosynthetically_active);
    try std.testing.expectEqual(@as(f64, 0), result.c3_feedback_fraction);
    try std.testing.expectEqual(@as(f64, 1), result.annual_termination_fraction);
}

test "runtime topology permits empty and unequal owner ranges" {
    const offsets = [_]usize{ 0, 2, 2, 3 };
    const inputs = [_]Inputs{ activeInputs(), activeInputs(), activeInputs() };
    var scratch: [3]Result = undefined;
    var destination: [3]Result = undefined;
    try computeRuntimeTopology(&offsets, &inputs, &scratch, &destination);
    try std.testing.expectEqual(@as(f64, 0.5), destination[0].annual_termination_fraction);
    try std.testing.expectEqual(@as(f64, 0.5), destination[2].annual_termination_fraction);
}

test "later invalid branch leaves every destination unchanged" {
    const offsets = [_]usize{ 0, 1, 2 };
    var invalid = activeInputs();
    invalid.heat_stress_h = std.math.nan(f64);
    const inputs = [_]Inputs{ activeInputs(), invalid };
    var scratch: [2]Result = undefined;
    var destination = [_]Result{
        .{ .c3_feedback_fraction = 41, .annual_termination_fraction = 1, .photosynthetically_active = true },
        .{ .c3_feedback_fraction = 42, .annual_termination_fraction = 1, .photosynthetically_active = true },
    };
    try std.testing.expectError(
        error.NonFiniteBranchFeedbackInput,
        computeRuntimeTopology(&offsets, &inputs, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 41), destination[0].c3_feedback_fraction);
    try std.testing.expectEqual(@as(f64, 42), destination[1].c3_feedback_fraction);
}
