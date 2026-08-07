const std = @import("std");
const ammonium = @import("ammonium_non_band_uptake.zig");
const nitrate = @import("nitrate_non_band_uptake.zig");

pub const State = struct {
    ammonium_non_band: ammonium.Result,
    ammonium_band: ammonium.Result,
    nitrate_non_band: nitrate.Result,
    nitrate_band: nitrate.Result,
};

/// UPTAKE.F 3297--3313. When the enclosing FZUP gate is inactive, reset
/// non-band NH4, band NH4, non-band NO3, then band NO3 in source order.
pub fn applyRuntimeAxes(
    nitrogen_uptake_active: []const bool,
    initial: []const State,
    scratch: []State,
    destination: []State,
) !void {
    if (nitrogen_uptake_active.len != initial.len or
        initial.len != scratch.len or
        initial.len != destination.len)
        return error.RootNitrogenResetDimensionMismatch;
    for (nitrogen_uptake_active, initial, scratch) |active, before, *candidate| {
        candidate.* = try apply(active, before);
    }
    @memcpy(destination, scratch);
}

pub fn apply(nitrogen_uptake_active: bool, before: State) !State {
    if (nitrogen_uptake_active) {
        try validateState(before);
        return before;
    }
    var result = before;
    result.ammonium_non_band = std.mem.zeroes(ammonium.Result);
    result.ammonium_band = std.mem.zeroes(ammonium.Result);
    result.nitrate_non_band = std.mem.zeroes(nitrate.Result);
    result.nitrate_band = std.mem.zeroes(nitrate.Result);
    return result;
}

fn validateState(state: State) !void {
    try validateZone(ammonium.Result, state.ammonium_non_band);
    try validateZone(ammonium.Result, state.ammonium_band);
    try validateZone(nitrate.Result, state.nitrate_non_band);
    try validateZone(nitrate.Result, state.nitrate_band);
}

fn validateZone(comptime T: type, zone: T) !void {
    inline for (@typeInfo(T).@"struct".fields) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(zone, field.name)))
            return error.NonFiniteRootNitrogenUptakeState;
    }
}

fn testState(value: f64) State {
    var state = std.mem.zeroes(State);
    state.ammonium_non_band.active = true;
    state.ammonium_non_band.population_uptake_g_n_per_step = value;
    state.ammonium_band.active = true;
    state.ammonium_band.population_uptake_g_n_per_step = value + 1;
    state.nitrate_non_band.active = true;
    state.nitrate_non_band.population_uptake_g_n_per_step = value + 2;
    state.nitrate_band.active = true;
    state.nitrate_band.population_uptake_g_n_per_step = value + 3;
    return state;
}

test "inactive nitrogen gate resets four zone groups in source order" {
    const result = try apply(false, testState(1));
    try std.testing.expectEqualDeep(std.mem.zeroes(State), result);
}

test "active nitrogen gate preserves all four zone groups" {
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

test "later active nonfinite state leaves destination unchanged" {
    const active = [_]bool{ false, true };
    var initial = [_]State{ testState(1), testState(2) };
    initial[1].nitrate_band.population_uptake_g_n_per_step =
        std.math.nan(f64);
    var scratch: [2]State = undefined;
    var destination = [_]State{ testState(41), testState(42) };
    try std.testing.expectError(
        error.NonFiniteRootNitrogenUptakeState,
        applyRuntimeAxes(&active, &initial, &scratch, &destination),
    );
    try std.testing.expectEqual(
        @as(f64, 41),
        destination[0].ammonium_non_band.population_uptake_g_n_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 42),
        destination[1].ammonium_non_band.population_uptake_g_n_per_step,
    );
}
