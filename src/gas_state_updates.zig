const std = @import("std");
const gas = @import("gas_transport.zig");

/// Extensive changes for one nonlinear gas solve. Every slice is cell-major
/// and has `cell_count * 7` entries. Positive phase exchange is gas -> water;
/// bubbling is zero or negative and represents aqueous mass vented externally.
pub const Changes = struct {
    spatial_net_g: []const f64,
    phase_to_dissolved_g: []const f64,
    phase_to_band_dissolved_g: []const f64,
    dissolved_external_net_g: []const f64,
    band_dissolved_external_net_g: []const f64,
    bubbling_from_dissolved_g: []const f64,
    bubbling_from_band_dissolved_g: []const f64,
};

pub const CommitSummary = struct {
    atmospheric_bubble_loss_g: [gas.species_count]f64,
};

/// Applies the TRNSFRS state equations atomically. Validation precedes every
/// write, so a negative or non-finite candidate cannot partially corrupt state.
pub fn commit(state: *gas.State, changes: Changes) !CommitSummary {
    const n = state.gaseous_mass_g.len;
    inline for (@typeInfo(Changes).@"struct".fields) |field| {
        if (@field(changes, field.name).len != n) return error.GasChangeSizeMismatch;
    }
    var summary = CommitSummary{ .atmospheric_bubble_loss_g = [_]f64{0} ** gas.species_count };
    for (0..n) |index| {
        const spatial = changes.spatial_net_g[index];
        const phase = changes.phase_to_dissolved_g[index];
        const band_phase = changes.phase_to_band_dissolved_g[index];
        const dissolved_external = changes.dissolved_external_net_g[index];
        const band_external = changes.band_dissolved_external_net_g[index];
        const bubble = changes.bubbling_from_dissolved_g[index];
        const band_bubble = changes.bubbling_from_band_dissolved_g[index];
        const values = [_]f64{ spatial, phase, band_phase, dissolved_external, band_external, bubble, band_bubble };
        for (values) |value| if (!std.math.isFinite(value)) return error.NonFiniteGasStateChange;
        if (bubble > 0 or band_bubble > 0) return error.InvalidGasBubblingSign;
        const gaseous_candidate = state.gaseous_mass_g[index] + spatial - phase - band_phase;
        const dissolved_candidate = state.dissolved_mass_g[index] + dissolved_external + phase + bubble;
        const band_candidate = state.band_dissolved_mass_g[index] + band_external + band_phase + band_bubble;
        if (!std.math.isFinite(gaseous_candidate) or !std.math.isFinite(dissolved_candidate) or !std.math.isFinite(band_candidate)) return error.NonFiniteGasStateCandidate;
        if (gaseous_candidate < -1e-12 or dissolved_candidate < -1e-12 or band_candidate < -1e-12) return error.NegativeGasStateCandidate;
    }
    for (0..n) |index| {
        state.gaseous_mass_g[index] = @max(0, state.gaseous_mass_g[index] + changes.spatial_net_g[index] - changes.phase_to_dissolved_g[index] - changes.phase_to_band_dissolved_g[index]);
        state.dissolved_mass_g[index] = @max(0, state.dissolved_mass_g[index] + changes.dissolved_external_net_g[index] + changes.phase_to_dissolved_g[index] + changes.bubbling_from_dissolved_g[index]);
        state.band_dissolved_mass_g[index] = @max(0, state.band_dissolved_mass_g[index] + changes.band_dissolved_external_net_g[index] + changes.phase_to_band_dissolved_g[index] + changes.bubbling_from_band_dissolved_g[index]);
        summary.atmospheric_bubble_loss_g[index % gas.species_count] -= changes.bubbling_from_dissolved_g[index] + changes.bubbling_from_band_dissolved_g[index];
    }
    return summary;
}

test "phase exchange conserves gas and both aqueous pools" {
    var state = try gas.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.gaseous_mass_g[5] = 10;
    state.dissolved_mass_g[5] = 2;
    state.band_dissolved_mass_g[5] = 1;
    const zero = [_]f64{0} ** gas.species_count;
    var nonband_phase = zero;
    var band_phase = zero;
    nonband_phase[5] = 2;
    band_phase[5] = 1;
    _ = try commit(&state, .{ .spatial_net_g = &zero, .phase_to_dissolved_g = &nonband_phase, .phase_to_band_dissolved_g = &band_phase, .dissolved_external_net_g = &zero, .band_dissolved_external_net_g = &zero, .bubbling_from_dissolved_g = &zero, .bubbling_from_band_dissolved_g = &zero });
    try std.testing.expectEqual(@as(f64, 7), state.gaseous_mass_g[5]);
    try std.testing.expectEqual(@as(f64, 4), state.dissolved_mass_g[5]);
    try std.testing.expectEqual(@as(f64, 2), state.band_dissolved_mass_g[5]);
}

test "bubbling is reported as external loss and failure is atomic" {
    var state = try gas.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.dissolved_mass_g[0] = 3;
    const zero = [_]f64{0} ** gas.species_count;
    var bubble = zero;
    bubble[0] = -1;
    const summary = try commit(&state, .{ .spatial_net_g = &zero, .phase_to_dissolved_g = &zero, .phase_to_band_dissolved_g = &zero, .dissolved_external_net_g = &zero, .band_dissolved_external_net_g = &zero, .bubbling_from_dissolved_g = &bubble, .bubbling_from_band_dissolved_g = &zero });
    try std.testing.expectEqual(@as(f64, 2), state.dissolved_mass_g[0]);
    try std.testing.expectEqual(@as(f64, 1), summary.atmospheric_bubble_loss_g[0]);
    bubble[0] = -4;
    try std.testing.expectError(error.NegativeGasStateCandidate, commit(&state, .{ .spatial_net_g = &zero, .phase_to_dissolved_g = &zero, .phase_to_band_dissolved_g = &zero, .dissolved_external_net_g = &zero, .band_dissolved_external_net_g = &zero, .bubbling_from_dissolved_g = &bubble, .bubbling_from_band_dissolved_g = &zero }));
    try std.testing.expectEqual(@as(f64, 2), state.dissolved_mass_g[0]);
}
