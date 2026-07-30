const std = @import("std");

pub const Diagnostics = struct {
    net_radiation_mj_per_step: f64,
    latent_heat_mj_per_step: f64,
    sensible_heat_mj_per_step: f64,
    storage_heat_mj_per_step: f64,
    vapor_convective_heat_mj_per_step: f64,
    emitted_thermal_radiation_mj_per_step: f64,
    evaporation_m3_per_step: f64,
    transpiration_m3_per_step: f64,
    ground_vapor_flux_m3_per_step: f64,
    ground_sensible_heat_mj_per_step: f64,
};

pub const State = struct {
    intercepted_water_m3: f64,
    diagnostics: Diagnostics,
    canopy_air_temperature_k: f64,
    canopy_air_vapor_volume_fraction: f64,
    surface_temperature_k: f64,
    surface_temperature_c: f64,
};

pub const FallbackInputs = struct {
    precipitation_retention_rate_m3_per_h: f64,
    timestep_h: f64,
};

/// UPTAKE.F 4284--4295. An admitted energy balance preserves the supplied
/// state; a rejected balance adds full-step precipitation and clears ten
/// diagnostics while retaining temperature and vapor state.
pub fn apply(
    energy_balance_active: bool,
    before: State,
    fallback_inputs: FallbackInputs,
) !State {
    if (energy_balance_active) {
        try validateState(before);
        return before;
    }
    try validatePreservedState(before);
    try validateFallbackInputs(fallback_inputs);
    var result = before;
    result.intercepted_water_m3 =
        before.intercepted_water_m3 +
        fallback_inputs.precipitation_retention_rate_m3_per_h *
            fallback_inputs.timestep_h;
    result.diagnostics.net_radiation_mj_per_step = 0;
    result.diagnostics.latent_heat_mj_per_step = 0;
    result.diagnostics.sensible_heat_mj_per_step = 0;
    result.diagnostics.storage_heat_mj_per_step = 0;
    result.diagnostics.vapor_convective_heat_mj_per_step = 0;
    result.diagnostics.emitted_thermal_radiation_mj_per_step = 0;
    result.diagnostics.evaporation_m3_per_step = 0;
    result.diagnostics.transpiration_m3_per_step = 0;
    result.diagnostics.ground_vapor_flux_m3_per_step = 0;
    result.diagnostics.ground_sensible_heat_mj_per_step = 0;
    if (!std.math.isFinite(result.intercepted_water_m3) or
        result.intercepted_water_m3 < 0)
        return error.InvalidStandingDeadEnergyFallbackResult;
    return result;
}

pub fn applyRuntimeSpecies(
    energy_balance_active: []const bool,
    initial: []const State,
    fallback_inputs: []const FallbackInputs,
    scratch: []State,
    destination: []State,
) !void {
    if (energy_balance_active.len != initial.len or
        initial.len != fallback_inputs.len or initial.len != scratch.len or
        initial.len != destination.len)
        return error.StandingDeadEnergyFallbackDimensionMismatch;
    for (
        energy_balance_active,
        initial,
        fallback_inputs,
        scratch,
    ) |active, before, inputs, *candidate| {
        candidate.* = try apply(active, before, inputs);
    }
    @memcpy(destination, scratch);
}

fn validateFallbackInputs(inputs: FallbackInputs) !void {
    if (!std.math.isFinite(inputs.precipitation_retention_rate_m3_per_h) or
        !std.math.isFinite(inputs.timestep_h))
        return error.NonFiniteStandingDeadEnergyFallbackInput;
    if (inputs.precipitation_retention_rate_m3_per_h < 0 or inputs.timestep_h < 0)
        return error.InvalidStandingDeadEnergyFallbackInput;
}

fn validatePreservedState(state: State) !void {
    inline for (.{
        state.intercepted_water_m3,
        state.canopy_air_temperature_k,
        state.canopy_air_vapor_volume_fraction,
        state.surface_temperature_k,
        state.surface_temperature_c,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteStandingDeadEnergyFallbackState;
    if (state.intercepted_water_m3 < 0 or state.canopy_air_temperature_k <= 0 or
        state.canopy_air_vapor_volume_fraction < 0 or
        state.surface_temperature_k <= 0)
        return error.InvalidStandingDeadEnergyFallbackState;
}

fn validateState(state: State) !void {
    try validatePreservedState(state);
    inline for (@typeInfo(Diagnostics).@"struct".fields) |field|
        if (!std.math.isFinite(@field(state.diagnostics, field.name)))
            return error.NonFiniteStandingDeadEnergyFallbackState;
}

fn testState(value: f64) State {
    return .{
        .intercepted_water_m3 = value,
        .diagnostics = .{
            .net_radiation_mj_per_step = value,
            .latent_heat_mj_per_step = value,
            .sensible_heat_mj_per_step = value,
            .storage_heat_mj_per_step = value,
            .vapor_convective_heat_mj_per_step = value,
            .emitted_thermal_radiation_mj_per_step = value,
            .evaporation_m3_per_step = value,
            .transpiration_m3_per_step = value,
            .ground_vapor_flux_m3_per_step = value,
            .ground_sensible_heat_mj_per_step = value,
        },
        .canopy_air_temperature_k = 280 + value,
        .canopy_air_vapor_volume_fraction = 0.01,
        .surface_temperature_k = 275 + value,
        .surface_temperature_c = 1.85 + value,
    };
}

test "rejected energy gate preserves source fallback assignment order" {
    const before = testState(1);
    const result = try apply(false, before, .{
        .precipitation_retention_rate_m3_per_h = 0.4,
        .timestep_h = 0.25,
    });
    try std.testing.expectEqual(@as(f64, 1.1), result.intercepted_water_m3);
    try std.testing.expectEqualDeep(std.mem.zeroes(Diagnostics), result.diagnostics);
    try std.testing.expectEqual(before.canopy_air_temperature_k, result.canopy_air_temperature_k);
    try std.testing.expectEqual(before.surface_temperature_k, result.surface_temperature_k);
}

test "active energy gate does not evaluate fallback-only inputs" {
    const before = testState(1);
    try std.testing.expectEqualDeep(before, try apply(true, before, .{
        .precipitation_retention_rate_m3_per_h = std.math.nan(f64),
        .timestep_h = std.math.nan(f64),
    }));
}

test "inactive branch overwrites invalid diagnostics but validates preserved state" {
    var before = testState(1);
    before.diagnostics.storage_heat_mj_per_step = std.math.nan(f64);
    const result = try apply(false, before, .{
        .precipitation_retention_rate_m3_per_h = 0,
        .timestep_h = 1,
    });
    try std.testing.expectEqual(@as(f64, 0), result.diagnostics.storage_heat_mj_per_step);
}

test "later invalid fallback input leaves runtime destination unchanged" {
    const active = [_]bool{ false, false };
    const initial = [_]State{ testState(1), testState(2) };
    const inputs = [_]FallbackInputs{
        .{ .precipitation_retention_rate_m3_per_h = 0.1, .timestep_h = 1 },
        .{ .precipitation_retention_rate_m3_per_h = -1, .timestep_h = 1 },
    };
    var scratch: [2]State = undefined;
    var destination = [_]State{ testState(41), testState(42) };
    try std.testing.expectError(
        error.InvalidStandingDeadEnergyFallbackInput,
        applyRuntimeSpecies(&active, &initial, &inputs, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 41), destination[0].intercepted_water_m3);
    try std.testing.expectEqual(@as(f64, 42), destination[1].intercepted_water_m3);
}
