const std = @import("std");

pub const ComponentFluxes = struct {
    net_radiation_megajoules_per_step: f64,
    latent_heat_megajoules_per_step: f64,
    sensible_heat_megajoules_per_step: f64,
    storage_heat_megajoules_per_step: f64,
    vapor_convective_heat_megajoules_per_step: f64,
    emitted_thermal_radiation_megajoules_per_step: f64,
    evaporation_m3_per_step: f64,
    transpiration_m3_per_step: f64,
    precipitation_convective_heat_megajoules_per_step: f64,
};

pub const Inputs = struct {
    live_canopy: ComponentFluxes,
    standing_dead: ComponentFluxes,
};

pub const Result = struct {
    net_radiation_megajoules_per_step: f64,
    latent_heat_megajoules_per_step: f64,
    sensible_heat_megajoules_per_step: f64,
    storage_heat_megajoules_per_step: f64,
    vapor_convective_heat_megajoules_per_step: f64,
    emitted_thermal_radiation_megajoules_per_step: f64,
    evaporation_m3_per_step: f64,
    transpiration_m3_per_step: f64,
    energy_balance_residual_megajoules_per_step: f64,
};

/// UPTAKE.F 4329--4336 total live-canopy plus standing-dead energy and water
/// diagnostics. Each runtime entry is exclusively owned by one cell/species.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    var result: Result = undefined;
    result.net_radiation_megajoules_per_step =
        inputs.live_canopy.net_radiation_megajoules_per_step +
        inputs.standing_dead.net_radiation_megajoules_per_step;
    result.latent_heat_megajoules_per_step =
        inputs.live_canopy.latent_heat_megajoules_per_step +
        inputs.standing_dead.latent_heat_megajoules_per_step;
    result.sensible_heat_megajoules_per_step =
        inputs.live_canopy.sensible_heat_megajoules_per_step +
        inputs.standing_dead.sensible_heat_megajoules_per_step;
    result.storage_heat_megajoules_per_step =
        inputs.live_canopy.storage_heat_megajoules_per_step +
        inputs.standing_dead.storage_heat_megajoules_per_step;
    result.vapor_convective_heat_megajoules_per_step =
        inputs.live_canopy.vapor_convective_heat_megajoules_per_step +
        inputs.standing_dead.vapor_convective_heat_megajoules_per_step;
    result.emitted_thermal_radiation_megajoules_per_step =
        inputs.live_canopy.emitted_thermal_radiation_megajoules_per_step +
        inputs.standing_dead.emitted_thermal_radiation_megajoules_per_step;
    result.evaporation_m3_per_step =
        inputs.live_canopy.evaporation_m3_per_step +
        inputs.standing_dead.evaporation_m3_per_step;
    result.transpiration_m3_per_step =
        inputs.live_canopy.transpiration_m3_per_step +
        inputs.standing_dead.transpiration_m3_per_step;
    result.energy_balance_residual_megajoules_per_step =
        result.storage_heat_megajoules_per_step -
        (result.net_radiation_megajoules_per_step +
            result.latent_heat_megajoules_per_step +
            result.sensible_heat_megajoules_per_step +
            result.vapor_convective_heat_megajoules_per_step +
            inputs.live_canopy.precipitation_convective_heat_megajoules_per_step +
            inputs.standing_dead.precipitation_convective_heat_megajoules_per_step);
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteCanopyEnergyAggregationResult;
    return result;
}

pub fn computeRuntimeSpecies(
    inputs: []const Inputs,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.CanopyEnergyAggregationDimensionMismatch;
    for (inputs, scratch) |species_inputs, *candidate|
        candidate.* = try compute(species_inputs);
    @memcpy(destination, scratch);
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(ComponentFluxes).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs.live_canopy, field.name)) or
            !std.math.isFinite(@field(inputs.standing_dead, field.name)))
            return error.NonFiniteCanopyEnergyAggregationInput;
    }
}

fn testComponent(scale: f64) ComponentFluxes {
    return .{
        .net_radiation_megajoules_per_step = 5 * scale,
        .latent_heat_megajoules_per_step = -2 * scale,
        .sensible_heat_megajoules_per_step = -1 * scale,
        .storage_heat_megajoules_per_step = 1.5 * scale,
        .vapor_convective_heat_megajoules_per_step = -0.5 * scale,
        .emitted_thermal_radiation_megajoules_per_step = 0.3 * scale,
        .evaporation_m3_per_step = -0.01 * scale,
        .transpiration_m3_per_step = -0.02 * scale,
        .precipitation_convective_heat_megajoules_per_step = 0,
    };
}

test "live plus standing-dead aggregation preserves eight source assignments" {
    const result = try compute(.{
        .live_canopy = testComponent(1),
        .standing_dead = testComponent(2),
    });
    try std.testing.expectEqual(@as(f64, 15), result.net_radiation_megajoules_per_step);
    try std.testing.expectEqual(@as(f64, -6), result.latent_heat_megajoules_per_step);
    try std.testing.expectEqual(@as(f64, -3), result.sensible_heat_megajoules_per_step);
    try std.testing.expectEqual(@as(f64, 4.5), result.storage_heat_megajoules_per_step);
    try std.testing.expectEqual(@as(f64, -1.5), result.vapor_convective_heat_megajoules_per_step);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.9),
        result.emitted_thermal_radiation_megajoules_per_step,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(@as(f64, -0.03), result.evaporation_m3_per_step, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.06), result.transpiration_m3_per_step, 1e-15);
}

test "combined canopy energy closure is conserved" {
    const result = try compute(.{
        .live_canopy = testComponent(1),
        .standing_dead = testComponent(2),
    });
    try std.testing.expectEqual(@as(f64, 0), result.energy_balance_residual_megajoules_per_step);
}

test "runtime cell-species entries have independent ownership" {
    const inputs = [_]Inputs{
        .{ .live_canopy = testComponent(1), .standing_dead = testComponent(2) },
        .{ .live_canopy = testComponent(10), .standing_dead = testComponent(20) },
    };
    var scratch: [2]Result = undefined;
    var destination: [2]Result = undefined;
    try computeRuntimeSpecies(&inputs, &scratch, &destination);
    try std.testing.expectEqual(@as(f64, 15), destination[0].net_radiation_megajoules_per_step);
    try std.testing.expectEqual(@as(f64, 150), destination[1].net_radiation_megajoules_per_step);
}

test "later nonfinite species leaves every destination unchanged" {
    var second: Inputs = .{
        .live_canopy = testComponent(10),
        .standing_dead = testComponent(20),
    };
    second.standing_dead.evaporation_m3_per_step = std.math.nan(f64);
    const inputs = [_]Inputs{
        .{ .live_canopy = testComponent(1), .standing_dead = testComponent(2) },
        second,
    };
    var scratch: [2]Result = undefined;
    var destination: [2]Result = undefined;
    destination[0] = try compute(inputs[0]);
    destination[1] = destination[0];
    destination[0].net_radiation_megajoules_per_step = 41;
    destination[1].net_radiation_megajoules_per_step = 42;
    try std.testing.expectError(
        error.NonFiniteCanopyEnergyAggregationInput,
        computeRuntimeSpecies(&inputs, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 41), destination[0].net_radiation_megajoules_per_step);
    try std.testing.expectEqual(@as(f64, 42), destination[1].net_radiation_megajoules_per_step);
}
