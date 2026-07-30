const std = @import("std");

pub const Parameters = struct {
    solid_snow_heat_capacity_mj_per_m3_k: f64,
    liquid_water_heat_capacity_mj_per_m3_k: f64,
    ice_heat_capacity_mj_per_m3_k: f64,
    active_heat_capacity_threshold_mj_per_k: f64,
    atmospheric_fallback_temperature_k: f64,
    absolute_heat_capacity_tolerance_mj_per_k: f64,
    relative_heat_capacity_tolerance: f64,
};

pub const AcceptedDrift = struct {
    solid_snow_water_equivalent_m3: f64,
    liquid_water_m3: f64,
    ice_volume_m3: f64,
    convective_heat_mj: f64,
};

pub const TopLayerState = struct {
    active: bool,
    solid_snow_water_equivalent_m3: f64,
    liquid_water_m3: f64,
    water_vapor_equivalent_m3: f64,
    ice_volume_m3: f64,
    heat_capacity_mj_per_k: f64,
    temperature_k: f64,
};

pub const Report = struct {
    sensible_energy_before_mj: f64,
    accepted_convective_heat_mj: f64,
    sensible_energy_after_mj: f64,
    fallback_energy_adjustment_mj: f64,
    used_atmospheric_fallback: bool,
};

/// Publishes accepted snow-drift water and heat into one top snow layer.
///
/// Traceability: REDIST.F lines 4228--4240. Solid snow, liquid water, and ice
/// increments are applied in source order. Previous sensible energy is retained
/// as `heat_capacity * temperature`; runtime phase heat capacities reconstruct
/// the new capacity before accepted drift heat determines temperature.
///
/// The preceding horizontal aggregation owns routing. This scalar, cell-local
/// kernel only commits its signed accepted result, so separate cells can execute
/// in parallel. Candidate calculation completes before mutation. When the new
/// capacity is below the runtime activity threshold, the source atmospheric
/// fallback is retained and its energy adjustment is reported explicitly.
pub fn publishAcceptedDrift(
    state: *TopLayerState,
    accepted: AcceptedDrift,
    parameters: Parameters,
) !Report {
    try validateInputs(state.*, accepted, parameters);
    const expected_old_capacity = try phaseHeatCapacity(state.*, parameters);
    try requireConsistentHeatCapacity(
        state.heat_capacity_mj_per_k,
        expected_old_capacity,
        parameters,
    );

    var candidate = state.*;
    candidate.solid_snow_water_equivalent_m3 = try addPhaseVolume(
        candidate.solid_snow_water_equivalent_m3,
        accepted.solid_snow_water_equivalent_m3,
    );
    candidate.liquid_water_m3 = try addPhaseVolume(
        candidate.liquid_water_m3,
        accepted.liquid_water_m3,
    );
    candidate.ice_volume_m3 = try addPhaseVolume(
        candidate.ice_volume_m3,
        accepted.ice_volume_m3,
    );

    const energy_before_mj =
        state.heat_capacity_mj_per_k * state.temperature_k;
    if (!std.math.isFinite(energy_before_mj))
        return error.NonFiniteSnowDriftPublicationResult;
    candidate.heat_capacity_mj_per_k =
        try phaseHeatCapacity(candidate, parameters);
    const expected_energy_mj =
        energy_before_mj + accepted.convective_heat_mj;
    if (!std.math.isFinite(expected_energy_mj))
        return error.NonFiniteSnowDriftPublicationResult;

    const used_fallback = candidate.heat_capacity_mj_per_k <=
        parameters.active_heat_capacity_threshold_mj_per_k;
    candidate.active = !used_fallback;
    candidate.temperature_k = if (used_fallback)
        parameters.atmospheric_fallback_temperature_k
    else
        expected_energy_mj / candidate.heat_capacity_mj_per_k;
    if (!std.math.isFinite(candidate.temperature_k) or
        candidate.temperature_k <= 0)
    {
        return error.InvalidSnowDriftTemperature;
    }
    const energy_after_mj =
        candidate.heat_capacity_mj_per_k * candidate.temperature_k;
    const fallback_adjustment_mj = energy_after_mj - expected_energy_mj;
    inline for (.{ energy_after_mj, fallback_adjustment_mj }) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteSnowDriftPublicationResult;
    }

    state.* = candidate;
    return .{
        .sensible_energy_before_mj = energy_before_mj,
        .accepted_convective_heat_mj = accepted.convective_heat_mj,
        .sensible_energy_after_mj = energy_after_mj,
        .fallback_energy_adjustment_mj = fallback_adjustment_mj,
        .used_atmospheric_fallback = used_fallback,
    };
}

fn phaseHeatCapacity(
    state: TopLayerState,
    parameters: Parameters,
) !f64 {
    const capacity_mj_per_k =
        parameters.solid_snow_heat_capacity_mj_per_m3_k *
        state.solid_snow_water_equivalent_m3 +
        parameters.liquid_water_heat_capacity_mj_per_m3_k *
            (state.liquid_water_m3 +
                state.water_vapor_equivalent_m3) +
        parameters.ice_heat_capacity_mj_per_m3_k *
            state.ice_volume_m3;
    if (!std.math.isFinite(capacity_mj_per_k))
        return error.NonFiniteSnowDriftPublicationResult;
    if (capacity_mj_per_k < 0)
        return error.InvalidSnowDriftHeatCapacity;
    return capacity_mj_per_k;
}

fn addPhaseVolume(current_m3: f64, increment_m3: f64) !f64 {
    const result_m3 = current_m3 + increment_m3;
    if (!std.math.isFinite(result_m3))
        return error.NonFiniteSnowDriftPublicationResult;
    if (result_m3 < 0) return error.NegativeSnowDriftPhaseInventory;
    return result_m3;
}

fn requireConsistentHeatCapacity(
    committed_mj_per_k: f64,
    expected_mj_per_k: f64,
    parameters: Parameters,
) !void {
    const tolerance_mj_per_k =
        parameters.absolute_heat_capacity_tolerance_mj_per_k +
        parameters.relative_heat_capacity_tolerance *
            @max(@abs(committed_mj_per_k), @abs(expected_mj_per_k));
    if (!std.math.isFinite(tolerance_mj_per_k))
        return error.NonFiniteSnowDriftPublicationResult;
    if (@abs(committed_mj_per_k - expected_mj_per_k) >
        tolerance_mj_per_k)
    {
        return error.InconsistentSnowDriftHeatCapacity;
    }
}

fn validateInputs(
    state: TopLayerState,
    accepted: AcceptedDrift,
    parameters: Parameters,
) !void {
    inline for (@typeInfo(TopLayerState).@"struct".fields) |field| {
        if (field.type == f64 and
            !std.math.isFinite(@field(state, field.name)))
        {
            return error.NonFiniteSnowDriftPublicationInput;
        }
    }
    inline for (@typeInfo(AcceptedDrift).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(accepted, field.name)))
            return error.NonFiniteSnowDriftPublicationInput;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(parameters, field.name)))
            return error.NonFiniteSnowDriftPublicationInput;
    }
    if (state.solid_snow_water_equivalent_m3 < 0 or
        state.liquid_water_m3 < 0 or
        state.water_vapor_equivalent_m3 < 0 or
        state.ice_volume_m3 < 0 or
        state.heat_capacity_mj_per_k < 0 or
        state.temperature_k <= 0 or
        parameters.solid_snow_heat_capacity_mj_per_m3_k <= 0 or
        parameters.liquid_water_heat_capacity_mj_per_m3_k <= 0 or
        parameters.ice_heat_capacity_mj_per_m3_k <= 0 or
        parameters.active_heat_capacity_threshold_mj_per_k < 0 or
        parameters.atmospheric_fallback_temperature_k <= 0 or
        parameters.absolute_heat_capacity_tolerance_mj_per_k < 0 or
        parameters.relative_heat_capacity_tolerance < 0)
    {
        return error.InvalidSnowDriftPublicationInput;
    }
}

const test_parameters: Parameters = .{
    .solid_snow_heat_capacity_mj_per_m3_k = 2,
    .liquid_water_heat_capacity_mj_per_m3_k = 4,
    .ice_heat_capacity_mj_per_m3_k = 1.5,
    .active_heat_capacity_threshold_mj_per_k = 0.01,
    .atmospheric_fallback_temperature_k = 265,
    .absolute_heat_capacity_tolerance_mj_per_k = 1e-12,
    .relative_heat_capacity_tolerance = 1e-12,
};

fn activeTestState() TopLayerState {
    return .{
        .active = true,
        .solid_snow_water_equivalent_m3 = 1,
        .liquid_water_m3 = 1,
        .water_vapor_equivalent_m3 = 0.5,
        .ice_volume_m3 = 1,
        .heat_capacity_mj_per_k = 9.5,
        .temperature_k = 260,
    };
}

test "REDIST accepted drift preserves source mass and energy order" {
    var state = activeTestState();
    const report = try publishAcceptedDrift(
        &state,
        .{
            .solid_snow_water_equivalent_m3 = 0.2,
            .liquid_water_m3 = -0.1,
            .ice_volume_m3 = 0.3,
            .convective_heat_mj = 50,
        },
        test_parameters,
    );

    try std.testing.expectEqual(@as(f64, 1.2), state.solid_snow_water_equivalent_m3);
    try std.testing.expectEqual(@as(f64, 0.9), state.liquid_water_m3);
    try std.testing.expectEqual(@as(f64, 1.3), state.ice_volume_m3);
    try std.testing.expectApproxEqAbs(
        @as(f64, 9.95),
        state.heat_capacity_mj_per_k,
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        report.sensible_energy_before_mj +
            report.accepted_convective_heat_mj,
        report.sensible_energy_after_mj,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 2520.0 / 9.95),
        state.temperature_k,
        1e-12,
    );
    try std.testing.expectEqual(@as(f64, 0), report.fallback_energy_adjustment_mj);
    try std.testing.expect(!report.used_atmospheric_fallback);
    try std.testing.expect(state.active);
}

test "signed drift export removes phases without losing sensible energy" {
    var state = activeTestState();
    const report = try publishAcceptedDrift(
        &state,
        .{
            .solid_snow_water_equivalent_m3 = -0.25,
            .liquid_water_m3 = -0.5,
            .ice_volume_m3 = -0.25,
            .convective_heat_mj = -100,
        },
        test_parameters,
    );

    try std.testing.expectEqual(@as(f64, 0.75), state.solid_snow_water_equivalent_m3);
    try std.testing.expectEqual(@as(f64, 0.5), state.liquid_water_m3);
    try std.testing.expectEqual(@as(f64, 0.75), state.ice_volume_m3);
    try std.testing.expectApproxEqAbs(
        report.sensible_energy_before_mj - 100,
        report.sensible_energy_after_mj,
        1e-12,
    );
}

test "sub-threshold layer reports atmospheric fallback energy adjustment" {
    var state = TopLayerState{
        .active = true,
        .solid_snow_water_equivalent_m3 = 0.01,
        .liquid_water_m3 = 0,
        .water_vapor_equivalent_m3 = 0,
        .ice_volume_m3 = 0,
        .heat_capacity_mj_per_k = 0.02,
        .temperature_k = 260,
    };
    var parameters = test_parameters;
    parameters.active_heat_capacity_threshold_mj_per_k = 0.01;
    const report = try publishAcceptedDrift(
        &state,
        .{
            .solid_snow_water_equivalent_m3 = -0.006,
            .liquid_water_m3 = 0,
            .ice_volume_m3 = 0,
            .convective_heat_mj = -1,
        },
        parameters,
    );

    try std.testing.expect(!state.active);
    try std.testing.expect(report.used_atmospheric_fallback);
    try std.testing.expectEqual(
        parameters.atmospheric_fallback_temperature_k,
        state.temperature_k,
    );
    try std.testing.expectApproxEqAbs(
        report.sensible_energy_after_mj -
            (report.sensible_energy_before_mj +
                report.accepted_convective_heat_mj),
        report.fallback_energy_adjustment_mj,
        1e-14,
    );
}

test "negative phase candidate fails atomically" {
    var state = activeTestState();
    const before = state;
    try std.testing.expectError(
        error.NegativeSnowDriftPhaseInventory,
        publishAcceptedDrift(
            &state,
            .{
                .solid_snow_water_equivalent_m3 = 1,
                .liquid_water_m3 = 1,
                .ice_volume_m3 = -2,
                .convective_heat_mj = 0,
            },
            test_parameters,
        ),
    );
    try std.testing.expectEqualDeep(before, state);
}

test "stale capacity and non-finite heat fail before mutation" {
    var state = activeTestState();
    state.heat_capacity_mj_per_k += 1;
    const before = state;
    try std.testing.expectError(
        error.InconsistentSnowDriftHeatCapacity,
        publishAcceptedDrift(
            &state,
            .{
                .solid_snow_water_equivalent_m3 = 0,
                .liquid_water_m3 = 0,
                .ice_volume_m3 = 0,
                .convective_heat_mj = 0,
            },
            test_parameters,
        ),
    );
    try std.testing.expectEqualDeep(before, state);

    state = activeTestState();
    const finite_before = state;
    try std.testing.expectError(
        error.NonFiniteSnowDriftPublicationInput,
        publishAcceptedDrift(
            &state,
            .{
                .solid_snow_water_equivalent_m3 = 0,
                .liquid_water_m3 = 0,
                .ice_volume_m3 = 0,
                .convective_heat_mj = std.math.nan(f64),
            },
            test_parameters,
        ),
    );
    try std.testing.expectEqualDeep(finite_before, state);
}
