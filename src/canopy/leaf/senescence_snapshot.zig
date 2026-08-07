const std = @import("std");

pub const NodeState = struct {
    leaf_carbon_g_c: []const f64,
    leaf_nitrogen_g_n: []const f64,
    leaf_phosphorus_g_p: []const f64,
    leaf_area_m2: []const f64,
};

pub const RecyclingFractions = struct {
    carbon: f64,
    nitrogen: f64,
    phosphorus: f64,
};

pub const Snapshot = struct {
    leaf_carbon_g_c: f64,
    leaf_nitrogen_g_n: f64,
    leaf_phosphorus_g_p: f64,
    leaf_area_m2: f64,
    remobilizable_carbon_g_c: f64,
    remobilizable_nitrogen_g_n: f64,
    remobilizable_phosphorus_g_p: f64,
    leaf_remobilization_fraction: f64,
    leaf_area_remobilization_fraction: f64,
};

/// grosub.f lines 2561--2588. Reads the selected logical runtime node, snapshots
/// its leaf state, calculates structurally constrained remobilizable C/N/P, and
/// applies the source leaf-availability cap in exact operation order.
pub fn prepare(
    leaf_remobilization_enabled: bool,
    state: NodeState,
    selected_node: usize,
    leaf_presence_threshold_g_c: f64,
    requested_remobilization_fraction: f64,
    recycling: RecyclingFractions,
) !?Snapshot {
    // IFLGP encloses the snapshot in the source. Returning null prevents stale
    // snapshot values from a preceding timestep being reused.
    if (!leaf_remobilization_enabled) return null;

    const runtime_node_count = state.leaf_carbon_g_c.len;
    inline for (.{ state.leaf_nitrogen_g_n, state.leaf_phosphorus_g_p, state.leaf_area_m2 }) |values|
        if (values.len != runtime_node_count)
            return error.LeafSenescenceSnapshotDimensionMismatch;
    if (selected_node >= runtime_node_count)
        return error.LeafSenescenceSnapshotIndexOutOfBounds;
    inline for (.{ leaf_presence_threshold_g_c, requested_remobilization_fraction }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteLeafSenescenceSnapshotInput;
    if (leaf_presence_threshold_g_c < 0 or requested_remobilization_fraction < 0)
        return error.InvalidLeafSenescenceSnapshotInput;
    inline for (@typeInfo(RecyclingFractions).@"struct".fields) |field| {
        const value = @field(recycling, field.name);
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return error.InvalidLeafSenescenceRecyclingFraction;
    }

    const leaf_carbon_g_c = state.leaf_carbon_g_c[selected_node];
    const leaf_nitrogen_g_n = state.leaf_nitrogen_g_n[selected_node];
    const leaf_phosphorus_g_p = state.leaf_phosphorus_g_p[selected_node];
    const leaf_area_m2 = state.leaf_area_m2[selected_node];
    inline for (.{ leaf_carbon_g_c, leaf_nitrogen_g_n, leaf_phosphorus_g_p, leaf_area_m2 }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLeafSenescenceNodeState;

    var remobilizable_carbon_g_c: f64 = 0;
    var remobilizable_nitrogen_g_n: f64 = 0;
    var remobilizable_phosphorus_g_p: f64 = 0;
    if (leaf_carbon_g_c > leaf_presence_threshold_g_c) {
        remobilizable_carbon_g_c = leaf_carbon_g_c * recycling.carbon;
        remobilizable_nitrogen_g_n = leaf_nitrogen_g_n *
            (recycling.nitrogen + (1.0 - recycling.nitrogen) * recycling.carbon);
        remobilizable_phosphorus_g_p = leaf_phosphorus_g_p *
            (recycling.phosphorus + (1.0 - recycling.phosphorus) * recycling.carbon);
    }

    var leaf_fraction = requested_remobilization_fraction;
    if (requested_remobilization_fraction * leaf_carbon_g_c > leaf_carbon_g_c and
        leaf_carbon_g_c > leaf_presence_threshold_g_c)
        leaf_fraction = @max(0.0, leaf_carbon_g_c / leaf_carbon_g_c);
    // Below the source presence threshold its condition does not cap FSNC. Do
    // not reproduce a resulting negative leaf pool; fail at the source instead.
    if (leaf_fraction > 1)
        return error.LeafSenescenceFractionExceedsAvailableMass;

    inline for (.{ remobilizable_carbon_g_c, remobilizable_nitrogen_g_n, remobilizable_phosphorus_g_p, leaf_fraction }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLeafSenescenceSnapshotResult;
    return .{
        .leaf_carbon_g_c = leaf_carbon_g_c,
        .leaf_nitrogen_g_n = leaf_nitrogen_g_n,
        .leaf_phosphorus_g_p = leaf_phosphorus_g_p,
        .leaf_area_m2 = leaf_area_m2,
        .remobilizable_carbon_g_c = remobilizable_carbon_g_c,
        .remobilizable_nitrogen_g_n = remobilizable_nitrogen_g_n,
        .remobilizable_phosphorus_g_p = remobilizable_phosphorus_g_p,
        .leaf_remobilization_fraction = leaf_fraction,
        .leaf_area_remobilization_fraction = leaf_fraction,
    };
}

fn sampleState() NodeState {
    return .{
        .leaf_carbon_g_c = &.{ 1, 4, 9 },
        .leaf_nitrogen_g_n = &.{ 0.1, 0.4, 0.9 },
        .leaf_phosphorus_g_p = &.{ 0.01, 0.04, 0.09 },
        .leaf_area_m2 = &.{ 0.2, 0.8, 1.8 },
    };
}

const sample_recycling: RecyclingFractions = .{ .carbon = 0.5, .nitrogen = 0.6, .phosphorus = 0.7 };

test "GROSUB snapshots selected logical node and remobilizable elements" {
    const result = (try prepare(true, sampleState(), 1, 1.0e-9, 0.25, sample_recycling)).?;
    try std.testing.expectEqual(@as(f64, 4), result.leaf_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0.8), result.leaf_area_m2);
    try std.testing.expectEqual(@as(f64, 2), result.remobilizable_carbon_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 0.32), result.remobilizable_nitrogen_g_n, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.034), result.remobilizable_phosphorus_g_p, 1.0e-15);
    try std.testing.expectEqual(@as(f64, 0.25), result.leaf_remobilization_fraction);
    try std.testing.expectEqual(result.leaf_remobilization_fraction, result.leaf_area_remobilization_fraction);
}

test "source cap limits above-threshold leaf removal to available mass" {
    const result = (try prepare(true, sampleState(), 2, 1.0e-9, 1.5, sample_recycling)).?;
    try std.testing.expectEqual(@as(f64, 1), result.leaf_remobilization_fraction);
}

test "below-threshold leaf has zero remobilizable pools" {
    const result = (try prepare(true, sampleState(), 0, 1, 0.4, sample_recycling)).?;
    try std.testing.expectEqual(@as(f64, 0), result.remobilizable_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), result.remobilizable_nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 0), result.remobilizable_phosphorus_g_p);
    try std.testing.expectEqual(@as(f64, 0.4), result.leaf_remobilization_fraction);
}

test "unsafe source over-removal below threshold fails instead of going negative" {
    try std.testing.expectError(
        error.LeafSenescenceFractionExceedsAvailableMass,
        prepare(true, sampleState(), 0, 1, 1.1, sample_recycling),
    );
}

test "disabled leaf remobilization returns no stale snapshot" {
    const malformed: NodeState = .{
        .leaf_carbon_g_c = &.{},
        .leaf_nitrogen_g_n = &.{1},
        .leaf_phosphorus_g_p = &.{},
        .leaf_area_m2 = &.{},
    };
    try std.testing.expectEqual(@as(?Snapshot, null), try prepare(false, malformed, 99, -1, -1, sample_recycling));
}

test "runtime dimensions indexes and non-finite state fail explicitly" {
    var malformed = sampleState();
    malformed.leaf_area_m2 = malformed.leaf_area_m2[0..2];
    try std.testing.expectError(error.LeafSenescenceSnapshotDimensionMismatch, prepare(true, malformed, 0, 0, 0.1, sample_recycling));
    try std.testing.expectError(error.LeafSenescenceSnapshotIndexOutOfBounds, prepare(true, sampleState(), 3, 0, 0.1, sample_recycling));
    const non_finite: NodeState = .{
        .leaf_carbon_g_c = &.{std.math.nan(f64)},
        .leaf_nitrogen_g_n = &.{0},
        .leaf_phosphorus_g_p = &.{0},
        .leaf_area_m2 = &.{0},
    };
    try std.testing.expectError(error.InvalidLeafSenescenceNodeState, prepare(true, non_finite, 0, 0, 0.1, sample_recycling));
}
