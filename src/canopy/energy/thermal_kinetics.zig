const std = @import("std");

pub const Inputs = struct {
    canopy_temperature_k: f64,
    thermal_adaptation_offset_k: f64,
    rubisco_co2_half_saturation_25c_umol_per_l: f64,
    rubisco_o2_half_saturation_25c_umol_per_l: f64,
    dissolved_o2_umol_per_l: f64,
};

pub const Result = struct {
    adapted_canopy_temperature_k: f64,
    gas_constant_temperature_j_per_mol: f64,
    entropy_temperature_j_per_mol: f64,
    high_low_temperature_inactivation: f64,
    carboxylation_temperature_factor: f64,
    oxygenation_temperature_factor: f64,
    electron_transport_temperature_factor: f64,
    co2_half_saturation_without_o2_umol_per_l: f64,
    o2_half_saturation_umol_per_l: f64,
    co2_half_saturation_with_o2_umol_per_l: f64,
};

/// STOMATE.F 88--106 thermal response and temperature-adjusted Michaelis
/// constants in exact source assignment order.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    var result: Result = undefined;
    result.adapted_canopy_temperature_k =
        inputs.canopy_temperature_k + inputs.thermal_adaptation_offset_k;
    result.gas_constant_temperature_j_per_mol =
        8.3143 * result.adapted_canopy_temperature_k;
    result.entropy_temperature_j_per_mol =
        710 * result.adapted_canopy_temperature_k;
    result.high_low_temperature_inactivation =
        1 +
        @exp(
            (197_500 - result.entropy_temperature_j_per_mol) /
                result.gas_constant_temperature_j_per_mol,
        ) +
        @exp(
            (result.entropy_temperature_j_per_mol - 222_500) /
                result.gas_constant_temperature_j_per_mol,
        );
    result.carboxylation_temperature_factor =
        @exp(
            26.237 -
                65_000 / result.gas_constant_temperature_j_per_mol,
        ) / result.high_low_temperature_inactivation;
    result.oxygenation_temperature_factor =
        @exp(
            24.220 -
                60_000 / result.gas_constant_temperature_j_per_mol,
        ) / result.high_low_temperature_inactivation;
    result.electron_transport_temperature_factor =
        @exp(
            17.362 -
                43_000 / result.gas_constant_temperature_j_per_mol,
        ) / result.high_low_temperature_inactivation;
    result.co2_half_saturation_without_o2_umol_per_l =
        inputs.rubisco_co2_half_saturation_25c_umol_per_l *
        @exp(
            16.136 -
                40_000 / result.gas_constant_temperature_j_per_mol,
        );
    result.o2_half_saturation_umol_per_l =
        inputs.rubisco_o2_half_saturation_25c_umol_per_l *
        @exp(
            8.067 -
                20_000 / result.gas_constant_temperature_j_per_mol,
        );
    result.co2_half_saturation_with_o2_umol_per_l =
        result.co2_half_saturation_without_o2_umol_per_l *
        (1 + inputs.dissolved_o2_umol_per_l /
            result.o2_half_saturation_umol_per_l);
    inline for (@typeInfo(Result).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteCanopyThermalKineticsResult;
    return result;
}

pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.CanopyThermalKineticsDimensionMismatch;
    for (inputs, scratch) |axis_inputs, *candidate|
        candidate.* = try compute(axis_inputs);
    @memcpy(destination, scratch);
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteCanopyThermalKineticsInput;
    if (inputs.canopy_temperature_k <= 0 or
        inputs.canopy_temperature_k + inputs.thermal_adaptation_offset_k <= 0 or
        inputs.rubisco_co2_half_saturation_25c_umol_per_l <= 0 or
        inputs.rubisco_o2_half_saturation_25c_umol_per_l <= 0 or
        inputs.dissolved_o2_umol_per_l < 0)
        return error.InvalidCanopyThermalKineticsInput;
}

fn testInputs() Inputs {
    return .{
        .canopy_temperature_k = 298.15,
        .thermal_adaptation_offset_k = 0,
        .rubisco_co2_half_saturation_25c_umol_per_l = 30,
        .rubisco_o2_half_saturation_25c_umol_per_l = 300,
        .dissolved_o2_umol_per_l = 210_000 * @exp(-6.175 - 0.0211 * 25),
    };
}

test "thermal kinetics preserves every STOMATE source intermediate" {
    const inputs = testInputs();
    const result = try compute(inputs);
    const rtk: f64 = 8.3143 * 298.15;
    const stk: f64 = 710 * 298.15;
    const inactivation =
        1 + @exp((197_500 - stk) / rtk) + @exp((stk - 222_500) / rtk);
    try std.testing.expectEqual(@as(f64, 298.15), result.adapted_canopy_temperature_k);
    try std.testing.expectEqual(rtk, result.gas_constant_temperature_j_per_mol);
    try std.testing.expectApproxEqAbs(
        stk,
        result.entropy_temperature_j_per_mol,
        1e-10,
    );
    try std.testing.expectApproxEqAbs(
        inactivation,
        result.high_low_temperature_inactivation,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 300) * @exp(8.067 - 20_000 / rtk),
        result.o2_half_saturation_umol_per_l,
        1e-12,
    );
}

test "shared outputs match existing gas environment owner" {
    const existing = @import("stomatal_resistance.zig");
    const gas = try existing.gasEnvironment(
        298.15,
        0,
        400,
        0.7,
        210_000,
        30,
        300,
    );
    const result = try compute(testInputs());
    try std.testing.expectEqual(
        gas.rubisco_carboxylation_temperature_factor,
        result.carboxylation_temperature_factor,
    );
    try std.testing.expectEqual(
        gas.rubisco_oxygenation_temperature_factor,
        result.oxygenation_temperature_factor,
    );
    try std.testing.expectEqual(
        gas.electron_transport_temperature_factor,
        result.electron_transport_temperature_factor,
    );
    try std.testing.expectEqual(
        gas.rubisco_co2_half_saturation_umol_per_l,
        result.co2_half_saturation_without_o2_umol_per_l,
    );
    try std.testing.expectEqual(
        gas.rubisco_co2_half_saturation_with_o2_umol_per_l,
        result.co2_half_saturation_with_o2_umol_per_l,
    );
}

test "runtime axes preserve independent thermal adaptation" {
    var second = testInputs();
    second.thermal_adaptation_offset_k = 5;
    const inputs = [_]Inputs{ testInputs(), second };
    var scratch: [2]Result = undefined;
    var destination: [2]Result = undefined;
    try computeRuntimeAxes(&inputs, &scratch, &destination);
    try std.testing.expectEqual(@as(f64, 298.15), destination[0].adapted_canopy_temperature_k);
    try std.testing.expectEqual(@as(f64, 303.15), destination[1].adapted_canopy_temperature_k);
}

test "later invalid axis leaves all destinations unchanged" {
    var second = testInputs();
    second.rubisco_o2_half_saturation_25c_umol_per_l = 0;
    const inputs = [_]Inputs{ testInputs(), second };
    var scratch: [2]Result = undefined;
    var destination: [2]Result = undefined;
    destination[0] = try compute(testInputs());
    destination[1] = destination[0];
    destination[0].adapted_canopy_temperature_k = 241;
    destination[1].adapted_canopy_temperature_k = 242;
    try std.testing.expectError(
        error.InvalidCanopyThermalKineticsInput,
        computeRuntimeAxes(&inputs, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 241), destination[0].adapted_canopy_temperature_k);
    try std.testing.expectEqual(@as(f64, 242), destination[1].adapted_canopy_temperature_k);
}
