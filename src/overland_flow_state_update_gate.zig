const std = @import("std");

pub const SaltSimulation = enum { disabled, enabled };

pub const Decision = struct {
    apply_overland_updates: bool,
    apply_salt_updates: bool,
};

/// Exact REDIST.F line 5056 outer predicate and line 5114 nested salt gate,
/// closed at lines 5157--5158. TQR is net runoff volume per model step (m3);
/// ZEROS is the cell-specific negligible-volume threshold (m3).
pub fn decide(
    net_runoff_volume_m3_per_step: f64,
    negligible_volume_m3: f64,
    salt_simulation: SaltSimulation,
) !Decision {
    if (!std.math.isFinite(net_runoff_volume_m3_per_step) or
        !std.math.isFinite(negligible_volume_m3))
        return error.NonFiniteOverlandFlowUpdateGateInput;
    if (negligible_volume_m3 < 0) return error.InvalidOverlandFlowUpdateThreshold;

    const apply_overland = @abs(net_runoff_volume_m3_per_step) > negligible_volume_m3;
    return .{
        .apply_overland_updates = apply_overland,
        .apply_salt_updates = apply_overland and salt_simulation == .enabled,
    };
}

test "REDIST overland update gate uses strict absolute runoff threshold" {
    const positive = try decide(2, 1, .disabled);
    const negative = try decide(-2, 1, .disabled);
    try std.testing.expect(positive.apply_overland_updates);
    try std.testing.expect(negative.apply_overland_updates);
    try std.testing.expect(!positive.apply_salt_updates);
}

test "runoff equal to threshold does not enter REDIST update block" {
    const positive = try decide(1, 1, .enabled);
    const negative = try decide(-1, 1, .enabled);
    try std.testing.expect(!positive.apply_overland_updates);
    try std.testing.expect(!negative.apply_overland_updates);
    try std.testing.expect(!positive.apply_salt_updates);
}

test "nested salt update requires both outer runoff gate and ISALTG" {
    const enabled = try decide(2, 1, .enabled);
    const disabled = try decide(2, 1, .disabled);
    const no_runoff = try decide(0, 1, .enabled);
    try std.testing.expect(enabled.apply_salt_updates);
    try std.testing.expect(!disabled.apply_salt_updates);
    try std.testing.expect(!no_runoff.apply_salt_updates);
}

test "zero threshold retains strict nonzero runoff behavior" {
    try std.testing.expect(!(try decide(0, 0, .enabled)).apply_overland_updates);
    try std.testing.expect((try decide(-1.0e-20, 0, .enabled)).apply_overland_updates);
}

test "invalid gate inputs fail immediately" {
    try std.testing.expectError(error.NonFiniteOverlandFlowUpdateGateInput, decide(std.math.nan(f64), 0, .enabled));
    try std.testing.expectError(error.NonFiniteOverlandFlowUpdateGateInput, decide(1, std.math.inf(f64), .enabled));
    try std.testing.expectError(error.InvalidOverlandFlowUpdateThreshold, decide(1, -1, .enabled));
}
