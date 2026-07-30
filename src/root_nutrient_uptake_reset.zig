const std = @import("std");
const nitrogen_reset = @import("root_nitrogen_uptake_reset.zig");
const phosphorus_reset = @import("root_phosphorus_uptake_reset.zig");

pub const State = struct {
    nitrogen: nitrogen_reset.State,
    phosphorus: phosphorus_reset.State,
};

/// UPTAKE.F 3713--3745. When the enclosing root nutrient-admission gate is
/// inactive, reset all nitrogen groups and then all phosphorus groups in the
/// exact order established by their source-order reset kernels.
pub fn applyRuntimeAxes(
    nutrient_uptake_active: []const bool,
    initial: []const State,
    scratch: []State,
    destination: []State,
) !void {
    if (nutrient_uptake_active.len != initial.len or
        initial.len != scratch.len or
        initial.len != destination.len)
        return error.RootNutrientResetDimensionMismatch;
    for (nutrient_uptake_active, initial, scratch) |active, before, *candidate| {
        candidate.* = try apply(active, before);
    }
    @memcpy(destination, scratch);
}

pub fn apply(nutrient_uptake_active: bool, before: State) !State {
    if (nutrient_uptake_active) {
        return .{
            .nitrogen = try nitrogen_reset.apply(true, before.nitrogen),
            .phosphorus = try phosphorus_reset.apply(true, before.phosphorus),
        };
    }

    // Preserve UPTAKE.F 3713--3744 ordering: all NH4/NO3 groups precede
    // every H2PO4/HPO4 group.
    const cleared_nitrogen = try nitrogen_reset.apply(false, before.nitrogen);
    const cleared_phosphorus = try phosphorus_reset.apply(false, before.phosphorus);
    return .{
        .nitrogen = cleared_nitrogen,
        .phosphorus = cleared_phosphorus,
    };
}

fn testState(value: f64) State {
    var state = std.mem.zeroes(State);
    state.nitrogen.ammonium_non_band.active = true;
    state.nitrogen.ammonium_non_band.population_uptake_g_n_per_step = value;
    state.nitrogen.ammonium_band.active = true;
    state.nitrogen.ammonium_band.population_uptake_g_n_per_step = value + 1;
    state.nitrogen.nitrate_non_band.active = true;
    state.nitrogen.nitrate_non_band.population_uptake_g_n_per_step = value + 2;
    state.nitrogen.nitrate_band.active = true;
    state.nitrogen.nitrate_band.population_uptake_g_n_per_step = value + 3;
    state.phosphorus.dihydrogen_phosphate_non_band.active = true;
    state.phosphorus.dihydrogen_phosphate_non_band.population_uptake_g_p_per_step =
        value + 4;
    state.phosphorus.dihydrogen_phosphate_band.active = true;
    state.phosphorus.dihydrogen_phosphate_band.population_uptake_g_p_per_step =
        value + 5;
    state.phosphorus.hydrogen_phosphate_non_band.active = true;
    state.phosphorus.hydrogen_phosphate_non_band.population_uptake_g_p_per_step =
        value + 6;
    state.phosphorus.hydrogen_phosphate_band.active = true;
    state.phosphorus.hydrogen_phosphate_band.population_uptake_g_p_per_step =
        value + 7;
    return state;
}

test "inactive nutrient gate resets nitrogen then phosphorus groups" {
    try std.testing.expectEqualDeep(std.mem.zeroes(State), try apply(false, testState(1)));
}

test "active nutrient gate preserves every nutrient group" {
    const before = testState(1);
    try std.testing.expectEqualDeep(before, try apply(true, before));
}

test "runtime axes preserve independent gate states" {
    const active = [_]bool{ false, true };
    const initial = [_]State{ testState(1), testState(10) };
    var scratch: [2]State = undefined;
    var destination: [2]State = undefined;
    try applyRuntimeAxes(&active, &initial, &scratch, &destination);
    try std.testing.expectEqualDeep(std.mem.zeroes(State), destination[0]);
    try std.testing.expectEqualDeep(initial[1], destination[1]);
}

test "dimension mismatch leaves destination unchanged" {
    const active = [_]bool{false};
    const initial = [_]State{ testState(1), testState(2) };
    var scratch: [2]State = undefined;
    var destination = [_]State{ testState(41), testState(42) };
    try std.testing.expectError(
        error.RootNutrientResetDimensionMismatch,
        applyRuntimeAxes(&active, &initial, &scratch, &destination),
    );
    try std.testing.expectEqualDeep(testState(41), destination[0]);
    try std.testing.expectEqualDeep(testState(42), destination[1]);
}

test "later active phosphorus NaN prevents partial destination commit" {
    const active = [_]bool{ false, true };
    var initial = [_]State{ testState(1), testState(2) };
    initial[1].phosphorus.hydrogen_phosphate_band.population_uptake_g_p_per_step =
        std.math.nan(f64);
    var scratch: [2]State = undefined;
    var destination = [_]State{ testState(41), testState(42) };
    try std.testing.expectError(
        error.NonFiniteRootPhosphorusUptakeState,
        applyRuntimeAxes(&active, &initial, &scratch, &destination),
    );
    try std.testing.expectEqualDeep(testState(41), destination[0]);
    try std.testing.expectEqualDeep(testState(42), destination[1]);
}
