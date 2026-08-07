const std = @import("std");
const phosphate = @import("dihydrogen_phosphate_non_band_uptake.zig");

pub const State = struct {
    dihydrogen_phosphate_non_band: phosphate.Result,
    dihydrogen_phosphate_band: phosphate.Result,
    hydrogen_phosphate_non_band: phosphate.Result,
    hydrogen_phosphate_band: phosphate.Result,
};

/// UPTAKE.F 3695--3711. When the enclosing FPUP gate is inactive, reset
/// non-band H2PO4, band H2PO4, non-band HPO4, then band HPO4 in source order.
pub fn applyRuntimeAxes(
    phosphorus_uptake_active: []const bool,
    initial: []const State,
    scratch: []State,
    destination: []State,
) !void {
    if (phosphorus_uptake_active.len != initial.len or
        initial.len != scratch.len or
        initial.len != destination.len)
        return error.RootPhosphorusResetDimensionMismatch;
    for (phosphorus_uptake_active, initial, scratch) |active, before, *candidate| {
        candidate.* = try apply(active, before);
    }
    @memcpy(destination, scratch);
}

pub fn apply(phosphorus_uptake_active: bool, before: State) !State {
    if (phosphorus_uptake_active) {
        try validateState(before);
        return before;
    }
    var result = before;
    result.dihydrogen_phosphate_non_band = std.mem.zeroes(phosphate.Result);
    result.dihydrogen_phosphate_band = std.mem.zeroes(phosphate.Result);
    result.hydrogen_phosphate_non_band = std.mem.zeroes(phosphate.Result);
    result.hydrogen_phosphate_band = std.mem.zeroes(phosphate.Result);
    return result;
}

fn validateState(state: State) !void {
    try validateZone(state.dihydrogen_phosphate_non_band);
    try validateZone(state.dihydrogen_phosphate_band);
    try validateZone(state.hydrogen_phosphate_non_band);
    try validateZone(state.hydrogen_phosphate_band);
}

fn validateZone(zone: phosphate.Result) !void {
    inline for (@typeInfo(phosphate.Result).@"struct".fields) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(zone, field.name)))
            return error.NonFiniteRootPhosphorusUptakeState;
    }
}

fn testState(value_g_p_per_step: f64) State {
    var state = std.mem.zeroes(State);
    state.dihydrogen_phosphate_non_band.active = true;
    state.dihydrogen_phosphate_non_band.population_uptake_g_p_per_step =
        value_g_p_per_step;
    state.dihydrogen_phosphate_band.active = true;
    state.dihydrogen_phosphate_band.population_uptake_g_p_per_step =
        value_g_p_per_step + 1;
    state.hydrogen_phosphate_non_band.active = true;
    state.hydrogen_phosphate_non_band.population_uptake_g_p_per_step =
        value_g_p_per_step + 2;
    state.hydrogen_phosphate_band.active = true;
    state.hydrogen_phosphate_band.population_uptake_g_p_per_step =
        value_g_p_per_step + 3;
    return state;
}

test "inactive phosphorus gate resets four zone groups in source order" {
    const result = try apply(false, testState(1));
    try std.testing.expectEqualDeep(std.mem.zeroes(State), result);
}

test "active phosphorus gate preserves all four zone groups" {
    const before = testState(1);
    try std.testing.expectEqualDeep(before, try apply(true, before));
}

test "runtime axes preserve independent active and inactive states" {
    const active = [_]bool{ false, true };
    const initial = [_]State{ testState(1), testState(10) };
    var scratch: [2]State = undefined;
    var destination: [2]State = undefined;
    try applyRuntimeAxes(&active, &initial, &scratch, &destination);
    try std.testing.expectEqualDeep(std.mem.zeroes(State), destination[0]);
    try std.testing.expectEqualDeep(initial[1], destination[1]);
}

test "runtime axes reject mismatched dimensions before mutation" {
    const active = [_]bool{false};
    const initial = [_]State{ testState(1), testState(2) };
    var scratch: [2]State = undefined;
    var destination = [_]State{ testState(41), testState(42) };
    try std.testing.expectError(
        error.RootPhosphorusResetDimensionMismatch,
        applyRuntimeAxes(&active, &initial, &scratch, &destination),
    );
    try std.testing.expectEqualDeep(testState(41), destination[0]);
    try std.testing.expectEqualDeep(testState(42), destination[1]);
}

test "later active nonfinite state leaves destination unchanged" {
    const active = [_]bool{ false, true };
    var initial = [_]State{ testState(1), testState(2) };
    initial[1].hydrogen_phosphate_band.population_uptake_g_p_per_step =
        std.math.nan(f64);
    var scratch: [2]State = undefined;
    var destination = [_]State{ testState(41), testState(42) };
    try std.testing.expectError(
        error.NonFiniteRootPhosphorusUptakeState,
        applyRuntimeAxes(&active, &initial, &scratch, &destination),
    );
    try std.testing.expectEqual(
        @as(f64, 41),
        destination[0].dihydrogen_phosphate_non_band.population_uptake_g_p_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 42),
        destination[1].dihydrogen_phosphate_non_band.population_uptake_g_p_per_step,
    );
}
