const std = @import("std");
const nutrient_reset = @import("root_nutrient_uptake_reset.zig");
const organic_exchange = @import("root_organic_exchange.zig");

pub const GasFluxes = struct {
    carbon_dioxide_g_c_per_step: f64,
    oxygen_g_o_per_step: f64,
    methane_g_c_per_step: f64,
    nitrous_oxide_g_n_per_step: f64,
    ammonia_g_n_per_step: f64,
};

pub const SoilRootFluxes = struct {
    carbon_dioxide_g_c_per_step: f64,
    oxygen_g_o_per_step: f64,
    methane_g_c_per_step: f64,
    nitrous_oxide_g_n_per_step: f64,
    ammonia_g_n_per_step: f64,
};

pub const RootProcessFluxes = struct {
    carbon_dioxide_g_c_per_step: f64,
    oxygen_g_o_per_step: f64,
};

pub const State = struct {
    root_atmosphere_phase_flux: GasFluxes,
    root_atmosphere_diffusive_flux: GasFluxes,
    soil_root_flux: SoilRootFluxes,
    root_process_flux: RootProcessFluxes,
    organic_exchange: organic_exchange.Result,
    root_oxygen_constraint: f64,
    nutrient_uptake: nutrient_reset.State,
    primary_axis_nitrogen_uptake_g_n_per_step: f64,
};

/// UPTAKE.F 3747--3803. Reset inactive root gas, organic, and nutrient flux
/// state in source order. The final RUPNF reset applies only to root axis one.
pub fn applyRuntimeAxes(
    root_active: []const bool,
    primary_root_axis: []const bool,
    initial: []const State,
    scratch: []State,
    destination: []State,
) !void {
    if (root_active.len != primary_root_axis.len or
        root_active.len != initial.len or
        initial.len != scratch.len or
        initial.len != destination.len)
        return error.RootActivityResetDimensionMismatch;
    for (
        root_active,
        primary_root_axis,
        initial,
        scratch,
    ) |active, is_primary, before, *candidate| {
        candidate.* = try apply(active, is_primary, before);
    }
    @memcpy(destination, scratch);
}

pub fn apply(root_active: bool, primary_root_axis: bool, before: State) !State {
    if (root_active) {
        try validateFinite(State, before);
        return before;
    }

    var result = before;
    // UPTAKE.F 3747--3763: preserve scalar assignment order.
    clearGasFluxes(&result.root_atmosphere_phase_flux);
    clearGasFluxes(&result.root_atmosphere_diffusive_flux);
    clearSoilRootFluxes(&result.soil_root_flux);
    result.root_process_flux.carbon_dioxide_g_c_per_step = 0;
    result.root_process_flux.oxygen_g_o_per_step = 0;

    // UPTAKE.F 3764--3768: K-major, then C/N/P within each substrate.
    for (0..organic_exchange.substrate_count) |substrate| {
        result.organic_exchange.carbon_exchange_g_c_per_step[substrate] = 0;
        result.organic_exchange.nitrogen_exchange_g_n_per_step[substrate] = 0;
        result.organic_exchange.phosphorus_exchange_g_p_per_step[substrate] = 0;
    }
    result.root_oxygen_constraint = 1;

    // UPTAKE.F 3770--3801: four nitrogen groups precede four phosphorus
    // groups through the already traced source-order reset.
    result.nutrient_uptake =
        try nutrient_reset.apply(false, result.nutrient_uptake);
    if (primary_root_axis)
        result.primary_axis_nitrogen_uptake_g_n_per_step = 0;
    return result;
}

fn clearGasFluxes(fluxes: *GasFluxes) void {
    fluxes.carbon_dioxide_g_c_per_step = 0;
    fluxes.oxygen_g_o_per_step = 0;
    fluxes.methane_g_c_per_step = 0;
    fluxes.nitrous_oxide_g_n_per_step = 0;
    fluxes.ammonia_g_n_per_step = 0;
}

fn clearSoilRootFluxes(fluxes: *SoilRootFluxes) void {
    fluxes.carbon_dioxide_g_c_per_step = 0;
    fluxes.oxygen_g_o_per_step = 0;
    fluxes.methane_g_c_per_step = 0;
    fluxes.nitrous_oxide_g_n_per_step = 0;
    fluxes.ammonia_g_n_per_step = 0;
}

fn validateFinite(comptime T: type, value: T) !void {
    switch (@typeInfo(T)) {
        .float => if (!std.math.isFinite(value))
            return error.NonFiniteRootActivityFluxState,
        .array => |array| for (value) |element|
            try validateFinite(array.child, element),
        .@"struct" => |structure| inline for (structure.fields) |field|
            try validateFinite(field.type, @field(value, field.name)),
        else => {},
    }
}

fn testState(value: f64) State {
    var state = std.mem.zeroes(State);
    state.root_atmosphere_phase_flux = testGasFluxes(value);
    state.root_atmosphere_diffusive_flux = testGasFluxes(value + 1);
    state.soil_root_flux = .{
        .carbon_dioxide_g_c_per_step = value + 2,
        .oxygen_g_o_per_step = value + 2,
        .methane_g_c_per_step = value + 2,
        .nitrous_oxide_g_n_per_step = value + 2,
        .ammonia_g_n_per_step = value + 2,
    };
    state.root_process_flux = .{
        .carbon_dioxide_g_c_per_step = value + 3,
        .oxygen_g_o_per_step = value + 3,
    };
    state.organic_exchange.carbon_exchange_g_c_per_step = @splat(value + 4);
    state.organic_exchange.nitrogen_exchange_g_n_per_step = @splat(value + 5);
    state.organic_exchange.phosphorus_exchange_g_p_per_step = @splat(value + 6);
    state.root_oxygen_constraint = value + 7;
    state.nutrient_uptake.nitrogen.ammonium_non_band.population_uptake_g_n_per_step =
        value + 8;
    state.nutrient_uptake.phosphorus.dihydrogen_phosphate_non_band.population_uptake_g_p_per_step =
        value + 9;
    state.primary_axis_nitrogen_uptake_g_n_per_step = value + 10;
    return state;
}

fn testGasFluxes(value: f64) GasFluxes {
    return .{
        .carbon_dioxide_g_c_per_step = value,
        .oxygen_g_o_per_step = value,
        .methane_g_c_per_step = value,
        .nitrous_oxide_g_n_per_step = value,
        .ammonia_g_n_per_step = value,
    };
}

test "inactive primary root axis resets complete source state" {
    const result = try apply(false, true, testState(1));
    var expected = std.mem.zeroes(State);
    expected.root_oxygen_constraint = 1;
    try std.testing.expectEqualDeep(expected, result);
}

test "inactive nonprimary axis preserves primary-only nitrogen total" {
    const before = testState(1);
    const result = try apply(false, false, before);
    try std.testing.expectEqual(
        before.primary_axis_nitrogen_uptake_g_n_per_step,
        result.primary_axis_nitrogen_uptake_g_n_per_step,
    );
    try std.testing.expectEqual(@as(f64, 1), result.root_oxygen_constraint);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.organic_exchange.phosphorus_exchange_g_p_per_step[4],
    );
}

test "active root axis preserves finite state" {
    const before = testState(1);
    try std.testing.expectEqualDeep(before, try apply(true, true, before));
}

test "runtime axes preserve independent active and primary-axis controls" {
    const active = [_]bool{ false, false, true };
    const primary = [_]bool{ true, false, true };
    const initial = [_]State{ testState(1), testState(10), testState(20) };
    var scratch: [3]State = undefined;
    var destination: [3]State = undefined;
    try applyRuntimeAxes(&active, &primary, &initial, &scratch, &destination);
    try std.testing.expectEqual(@as(f64, 0), destination[0].primary_axis_nitrogen_uptake_g_n_per_step);
    try std.testing.expectEqual(
        initial[1].primary_axis_nitrogen_uptake_g_n_per_step,
        destination[1].primary_axis_nitrogen_uptake_g_n_per_step,
    );
    try std.testing.expectEqualDeep(initial[2], destination[2]);
}

test "later active NaN leaves destination unchanged" {
    const active = [_]bool{ false, true };
    const primary = [_]bool{ true, true };
    var initial = [_]State{ testState(1), testState(2) };
    initial[1].organic_exchange.nitrogen_exchange_g_n_per_step[4] =
        std.math.nan(f64);
    var scratch: [2]State = undefined;
    var destination = [_]State{ testState(41), testState(42) };
    try std.testing.expectError(
        error.NonFiniteRootActivityFluxState,
        applyRuntimeAxes(&active, &primary, &initial, &scratch, &destination),
    );
    try std.testing.expectEqualDeep(testState(41), destination[0]);
    try std.testing.expectEqualDeep(testState(42), destination[1]);
}
