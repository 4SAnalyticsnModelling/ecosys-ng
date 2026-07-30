const std = @import("std");

pub const Salt = enum(u8) {
    aluminum,
    iron,
    calcium,
    magnesium,
    sodium,
    potassium,
    sulfate,
    chloride,
};

pub const salt_count = @typeInfo(Salt).@"enum".fields.len;
pub const SaltValues = [salt_count]f64;

pub const Inputs = struct {
    dynamic_salts_enabled: bool,
    soil_allocated_water_m3: f64,
    root_water_m3: f64,
    combined_water_m3: f64,
    water_transport_factor: f64,
    root_area_path_m: f64,
    root_water_uptake_m3_per_step: f64,
    population: f64,
    shorter_time_fraction: f64,
    aqueous_diffusivity_m2_per_step: SaltValues,
    inhibition_concentration_mol_m3: SaltValues,
};

pub const State = struct {
    soil_mol: SaltValues,
    root_mol: SaltValues,
    cumulative_uptake_mol: SaltValues,
};

pub const Diagnostics = struct {
    soil_concentration_mol_m3: SaltValues,
    root_concentration_mol_m3: SaltValues,
    diffusivity_m3_per_step: SaltValues,
    convective_flow_mol_per_step: SaltValues,
    diffusive_convective_flow_mol_per_step: SaltValues,
    concentration_gradient_transfer_mol_per_step: SaltValues,
    uptake_mol_per_step: SaltValues,
};

pub const Result = struct {
    state: State,
    diagnostics: Diagnostics,
};

const source_order = [_]Salt{
    .aluminum,
    .iron,
    .calcium,
    .magnesium,
    .sodium,
    .potassium,
    .sulfate,
    .chloride,
};

/// UPTAKE.F 2647--2804 dynamic salt gate. Runtime axis entries are computed
/// into caller-owned scratch and committed only after all entries succeed.
pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    initial: []const State,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != initial.len or
        inputs.len != scratch.len or
        inputs.len != destination.len)
        return error.RootSaltUptakeDimensionMismatch;
    for (inputs, initial, scratch) |axis_inputs, axis_state, *candidate|
        candidate.* = try compute(axis_inputs, axis_state);
    @memcpy(destination, scratch);
}

pub fn compute(inputs: Inputs, initial: State) !Result {
    var result = Result{
        .state = initial,
        .diagnostics = std.mem.zeroes(Diagnostics),
    };
    if (!inputs.dynamic_salts_enabled) {
        try validateResult(result);
        return result;
    }
    try validateActiveInputs(inputs, initial);
    inline for (source_order) |salt| {
        const index = @intFromEnum(salt);
        const soil_concentration =
            @max(0, result.state.soil_mol[index] / inputs.soil_allocated_water_m3);
        const root_concentration =
            @max(0, result.state.root_mol[index] / inputs.root_water_m3);
        const diffusivity = inputs.water_transport_factor *
            inputs.aqueous_diffusivity_m2_per_step[index] *
            inputs.root_area_path_m;
        const convective_flow =
            inputs.root_water_uptake_m3_per_step * soil_concentration;
        const combined_flow = convective_flow +
            diffusivity * (soil_concentration - root_concentration);
        const gradient_transfer =
            (inputs.root_water_m3 * @max(0, result.state.soil_mol[index]) -
                inputs.soil_allocated_water_m3 *
                    @max(0, result.state.root_mol[index])) /
            inputs.combined_water_m3 * inputs.shorter_time_fraction;
        const population_flow = combined_flow * inputs.population;
        const limited = if (combined_flow > 0)
            @min(@max(0, gradient_transfer), population_flow)
        else
            @max(@min(0, gradient_transfer), population_flow);
        const uptake = limited /
            (1 + root_concentration /
                inputs.inhibition_concentration_mol_m3[index]);
        result.diagnostics.soil_concentration_mol_m3[index] = soil_concentration;
        result.diagnostics.root_concentration_mol_m3[index] = root_concentration;
        result.diagnostics.diffusivity_m3_per_step[index] = diffusivity;
        result.diagnostics.convective_flow_mol_per_step[index] = convective_flow;
        result.diagnostics.diffusive_convective_flow_mol_per_step[index] =
            combined_flow;
        result.diagnostics.concentration_gradient_transfer_mol_per_step[index] =
            gradient_transfer;
        result.diagnostics.uptake_mol_per_step[index] = uptake;
        result.state.soil_mol[index] -= uptake;
        result.state.root_mol[index] += uptake;
        result.state.cumulative_uptake_mol[index] += uptake;
    }
    try validateResult(result);
    return result;
}

fn validateActiveInputs(inputs: Inputs, initial: State) !void {
    inline for (.{
        inputs.soil_allocated_water_m3,
        inputs.root_water_m3,
        inputs.combined_water_m3,
        inputs.population,
        inputs.shorter_time_fraction,
    }) |value|
        if (!std.math.isFinite(value) or value <= 0)
            return error.InvalidRootSaltUptakeInput;
    inline for (.{
        inputs.water_transport_factor,
        inputs.root_area_path_m,
    }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootSaltUptakeInput;
    if (!std.math.isFinite(inputs.root_water_uptake_m3_per_step))
        return error.InvalidRootSaltUptakeInput;
    inline for (0..salt_count) |index| {
        if (!std.math.isFinite(inputs.aqueous_diffusivity_m2_per_step[index]) or
            inputs.aqueous_diffusivity_m2_per_step[index] < 0 or
            !std.math.isFinite(inputs.inhibition_concentration_mol_m3[index]) or
            inputs.inhibition_concentration_mol_m3[index] <= 0 or
            !std.math.isFinite(initial.soil_mol[index]) or
            initial.soil_mol[index] < 0 or
            !std.math.isFinite(initial.root_mol[index]) or
            initial.root_mol[index] < 0 or
            !std.math.isFinite(initial.cumulative_uptake_mol[index]))
            return error.InvalidRootSaltUptakeInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (0..salt_count) |index| {
        inline for (.{
            result.state.soil_mol[index],
            result.state.root_mol[index],
            result.state.cumulative_uptake_mol[index],
            result.diagnostics.soil_concentration_mol_m3[index],
            result.diagnostics.root_concentration_mol_m3[index],
            result.diagnostics.diffusivity_m3_per_step[index],
            result.diagnostics.convective_flow_mol_per_step[index],
            result.diagnostics.diffusive_convective_flow_mol_per_step[index],
            result.diagnostics.concentration_gradient_transfer_mol_per_step[index],
            result.diagnostics.uptake_mol_per_step[index],
        }) |value|
            if (!std.math.isFinite(value))
                return error.NonFiniteRootSaltUptakeResult;
    }
}

fn testInputs() Inputs {
    return .{
        .dynamic_salts_enabled = true,
        .soil_allocated_water_m3 = 2,
        .root_water_m3 = 1,
        .combined_water_m3 = 3,
        .water_transport_factor = 0.5,
        .root_area_path_m = 1,
        .root_water_uptake_m3_per_step = 0.1,
        .population = 2,
        .shorter_time_fraction = 0.25,
        .aqueous_diffusivity_m2_per_step = .{1} ** salt_count,
        .inhibition_concentration_mol_m3 = .{2} ** salt_count,
    };
}

fn testState() State {
    return .{
        .soil_mol = .{4} ** salt_count,
        .root_mol = .{1} ** salt_count,
        .cumulative_uptake_mol = .{0} ** salt_count,
    };
}

test "all eight salts preserve soil-root mole balance" {
    const initial = testState();
    const result = try compute(testInputs(), initial);
    inline for (0..salt_count) |index| {
        try std.testing.expectApproxEqAbs(
            initial.soil_mol[index] + initial.root_mol[index],
            result.state.soil_mol[index] + result.state.root_mol[index],
            1e-12,
        );
        try std.testing.expect(result.diagnostics.uptake_mol_per_step[index] > 0);
        try std.testing.expectApproxEqAbs(
            result.diagnostics.uptake_mol_per_step[index],
            result.state.cumulative_uptake_mol[index],
            1e-12,
        );
    }
}

test "negative combined flow preserves signed legacy limiter" {
    var inputs = testInputs();
    inputs.root_water_uptake_m3_per_step = 0;
    var initial = testState();
    initial.soil_mol = .{0.2} ** salt_count;
    initial.root_mol = .{3} ** salt_count;
    const result = try compute(inputs, initial);
    try std.testing.expect(
        result.diagnostics.diffusive_convective_flow_mol_per_step[0] < 0,
    );
    try std.testing.expect(result.diagnostics.uptake_mol_per_step[0] < 0);
}

test "disabled dynamic salts leave state unchanged" {
    var inputs = testInputs();
    inputs.dynamic_salts_enabled = false;
    const initial = testState();
    const result = try compute(inputs, initial);
    try std.testing.expectEqualDeep(initial, result.state);
    try std.testing.expectEqualDeep(std.mem.zeroes(Diagnostics), result.diagnostics);
}

test "runtime axes fail atomically on later invalid inhibition" {
    var inputs = [_]Inputs{ testInputs(), testInputs() };
    inputs[1].inhibition_concentration_mol_m3[3] = 0;
    const initial = [_]State{ testState(), testState() };
    var scratch: [2]Result = undefined;
    var destination = [_]Result{
        .{ .state = testState(), .diagnostics = std.mem.zeroes(Diagnostics) },
        .{ .state = testState(), .diagnostics = std.mem.zeroes(Diagnostics) },
    };
    destination[0].state.soil_mol[0] = 41;
    destination[1].state.soil_mol[0] = 42;
    try std.testing.expectError(
        error.InvalidRootSaltUptakeInput,
        computeRuntimeAxes(&inputs, &initial, &scratch, &destination),
    );
    try std.testing.expectEqual(@as(f64, 41), destination[0].state.soil_mol[0]);
    try std.testing.expectEqual(@as(f64, 42), destination[1].state.soil_mol[0]);
}
