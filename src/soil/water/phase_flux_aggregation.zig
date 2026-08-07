const std = @import("std");

/// Phase and vapor-change diagnostics in exact REDIST publication order.
pub const PhaseFlux = struct {
    micropore_phase_change_water_m3_per_step: f64 = 0,
    macropore_phase_change_water_m3_per_step: f64 = 0,
    phase_change_latent_heat_megajoules_per_step: f64 = 0,
    evaporation_condensation_water_m3_per_step: f64 = 0,
    evaporation_condensation_latent_heat_megajoules_per_step: f64 = 0,
};

pub const State = struct {
    net_flux: PhaseFlux,
};

/// Publishes one layer's phase-change and vapor-change flux increments.
///
/// Traceability: REDIST.F lines 3950--3954. This is an accounting kernel, not
/// the legacy freeze-thaw formulation: `increment` is supplied by ecosys-ng's
/// conservative Dall'Amico enthalpy/Newton-Picard phase solver. The operation
/// order remains micropore water, macropore water, phase latent heat,
/// evaporation/condensation water, then its latent heat. Candidate state is
/// committed only after all five additions are finite.
pub fn aggregate(increment: PhaseFlux, state: *State) !void {
    try validateFlux(increment);
    try validateFlux(state.net_flux);
    var candidate = state.net_flux;
    inline for (@typeInfo(PhaseFlux).@"struct".fields) |field| {
        const result = @field(candidate, field.name) +
            @field(increment, field.name);
        if (!std.math.isFinite(result))
            return error.NonFiniteSoilPhaseFluxResult;
        @field(candidate, field.name) = result;
    }
    state.net_flux = candidate;
}

fn validateFlux(flux: PhaseFlux) !void {
    inline for (@typeInfo(PhaseFlux).@"struct".fields) |field|
        if (!std.math.isFinite(@field(flux, field.name)))
            return error.NonFiniteSoilPhaseFluxInput;
}

fn filledFlux(value: f64) PhaseFlux {
    var result: PhaseFlux = undefined;
    inline for (@typeInfo(PhaseFlux).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn expectFlux(actual: PhaseFlux, expected: f64) !void {
    inline for (@typeInfo(PhaseFlux).@"struct".fields) |field|
        try std.testing.expectEqual(expected, @field(actual, field.name));
}

test "all five phase diagnostics accumulate in source order" {
    var state = State{ .net_flux = .{
        .micropore_phase_change_water_m3_per_step = 10,
        .macropore_phase_change_water_m3_per_step = 20,
        .phase_change_latent_heat_megajoules_per_step = 30,
        .evaporation_condensation_water_m3_per_step = 40,
        .evaporation_condensation_latent_heat_megajoules_per_step = 50,
    } };
    try aggregate(.{
        .micropore_phase_change_water_m3_per_step = 1,
        .macropore_phase_change_water_m3_per_step = 2,
        .phase_change_latent_heat_megajoules_per_step = 3,
        .evaporation_condensation_water_m3_per_step = 4,
        .evaporation_condensation_latent_heat_megajoules_per_step = 5,
    }, &state);
    try std.testing.expectEqual(
        @as(f64, 11),
        state.net_flux.micropore_phase_change_water_m3_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 22),
        state.net_flux.macropore_phase_change_water_m3_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 33),
        state.net_flux.phase_change_latent_heat_megajoules_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 44),
        state.net_flux.evaporation_condensation_water_m3_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 55),
        state.net_flux.evaporation_condensation_latent_heat_megajoules_per_step,
    );
}

test "signed increments retain exact accumulated component balances" {
    var state = State{ .net_flux = .{} };
    try aggregate(filledFlux(7), &state);
    try aggregate(filledFlux(-3), &state);
    try expectFlux(state.net_flux, 4);
}

test "zero phase publication leaves prior diagnostics unchanged" {
    var state = State{ .net_flux = filledFlux(9) };
    try aggregate(.{}, &state);
    try expectFlux(state.net_flux, 9);
}

test "nonfinite source and state fail before mutation" {
    var increment = filledFlux(3);
    increment.evaporation_condensation_latent_heat_megajoules_per_step =
        std.math.nan(f64);
    var state = State{ .net_flux = filledFlux(5) };
    try std.testing.expectError(
        error.NonFiniteSoilPhaseFluxInput,
        aggregate(increment, &state),
    );
    try expectFlux(state.net_flux, 5);

    increment = filledFlux(3);
    state.net_flux.phase_change_latent_heat_megajoules_per_step =
        std.math.inf(f64);
    try std.testing.expectError(
        error.NonFiniteSoilPhaseFluxInput,
        aggregate(increment, &state),
    );
    try std.testing.expect(std.math.isInf(
        state.net_flux.phase_change_latent_heat_megajoules_per_step,
    ));
}

test "late arithmetic overflow preserves every diagnostic atomically" {
    var increment = filledFlux(1);
    increment.evaporation_condensation_latent_heat_megajoules_per_step =
        std.math.floatMax(f64);
    var state = State{ .net_flux = filledFlux(10) };
    state.net_flux.evaporation_condensation_latent_heat_megajoules_per_step =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteSoilPhaseFluxResult,
        aggregate(increment, &state),
    );
    try std.testing.expectEqual(
        @as(f64, 10),
        state.net_flux.micropore_phase_change_water_m3_per_step,
    );
    try std.testing.expectEqual(
        std.math.floatMax(f64),
        state.net_flux.evaporation_condensation_latent_heat_megajoules_per_step,
    );
}
