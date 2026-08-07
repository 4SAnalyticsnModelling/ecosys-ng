const std = @import("std");

pub const substrate_count = 5;
pub const SubstrateValues = [substrate_count]f64;

pub const RootPools = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const SoilDissolvedPools = struct {
    carbon_g_c: SubstrateValues,
    nitrogen_g_n: SubstrateValues,
    phosphorus_g_p: SubstrateValues,
};

pub const Inputs = struct {
    soil_micropore_water_m3: f64,
    substrate_fraction: SubstrateValues,
    root_water_m3: f64,
    minimum_soil_water_m3: f64,
    minimum_root_water_m3: f64,
    minimum_soil_carbon_g_c: f64,
    minimum_root_carbon_g_c: f64,
    soil: SoilDissolvedPools,
    root: RootPools,
    exchange_rate_per_h: f64,
    solute_timestep_h: f64,
};

pub const Result = struct {
    carbon_exchange_g_c_per_step: SubstrateValues,
    nitrogen_exchange_g_n_per_step: SubstrateValues,
    phosphorus_exchange_g_p_per_step: SubstrateValues,
};

/// UPTAKE.F 2855--2904 fixed scientific substrate loop. All model axes
/// outside the five substrate classes remain caller-sized runtime slices.
pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.RootOrganicExchangeDimensionMismatch;
    for (inputs, scratch) |axis_inputs, *candidate|
        candidate.* = try compute(axis_inputs);
    @memcpy(destination, scratch);
}

pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    var result = std.mem.zeroes(Result);
    inline for (0..substrate_count) |substrate| {
        const substrate_water =
            inputs.soil_micropore_water_m3 *
            inputs.substrate_fraction[substrate];
        if (substrate_water > inputs.minimum_soil_water_m3 and
            inputs.root_water_m3 > inputs.minimum_root_water_m3)
        {
            const combined_water = substrate_water + inputs.root_water_m3;
            const soil_carbon = @max(0, inputs.soil.carbon_g_c[substrate]);
            const root_carbon =
                @min(1.0e3 * inputs.root_water_m3, inputs.root.carbon_g_c);
            const carbon_equilibrium =
                (soil_carbon * inputs.root_water_m3 -
                    root_carbon * substrate_water) /
                combined_water;
            result.carbon_exchange_g_c_per_step[substrate] =
                inputs.exchange_rate_per_h *
                carbon_equilibrium *
                inputs.solute_timestep_h;
            if (inputs.soil.carbon_g_c[substrate] >
                inputs.minimum_soil_carbon_g_c and
                inputs.root.carbon_g_c > inputs.minimum_root_carbon_g_c)
            {
                const combined_carbon =
                    inputs.soil.carbon_g_c[substrate] +
                    inputs.root.carbon_g_c;
                const soil_nitrogen =
                    @max(0, inputs.soil.nitrogen_g_n[substrate]);
                const soil_phosphorus =
                    @max(0, inputs.soil.phosphorus_g_p[substrate]);
                const root_nitrogen = @max(0, 0.1 * inputs.root.nitrogen_g_n);
                const root_phosphorus = @max(0, 0.1 * inputs.root.phosphorus_g_p);
                const nitrogen_equilibrium =
                    (soil_nitrogen * inputs.root.carbon_g_c -
                        root_nitrogen *
                            inputs.soil.carbon_g_c[substrate]) /
                    combined_carbon;
                const phosphorus_equilibrium =
                    (soil_phosphorus * inputs.root.carbon_g_c -
                        root_phosphorus *
                            inputs.soil.carbon_g_c[substrate]) /
                    combined_carbon;
                result.nitrogen_exchange_g_n_per_step[substrate] =
                    inputs.exchange_rate_per_h *
                    nitrogen_equilibrium *
                    inputs.solute_timestep_h;
                result.phosphorus_exchange_g_p_per_step[substrate] =
                    inputs.exchange_rate_per_h *
                    phosphorus_equilibrium *
                    inputs.solute_timestep_h;
            }
        }
    }
    try validateResult(result);
    return result;
}

fn validate(inputs: Inputs) !void {
    inline for (.{
        inputs.soil_micropore_water_m3,
        inputs.root_water_m3,
        inputs.minimum_soil_water_m3,
        inputs.minimum_root_water_m3,
        inputs.minimum_soil_carbon_g_c,
        inputs.minimum_root_carbon_g_c,
        inputs.root.carbon_g_c,
        inputs.root.nitrogen_g_n,
        inputs.root.phosphorus_g_p,
        inputs.exchange_rate_per_h,
        inputs.solute_timestep_h,
    }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootOrganicExchangeInput;
    inline for (0..substrate_count) |substrate| {
        if (!std.math.isFinite(inputs.substrate_fraction[substrate]) or
            inputs.substrate_fraction[substrate] < 0 or
            inputs.substrate_fraction[substrate] > 1 or
            !std.math.isFinite(inputs.soil.carbon_g_c[substrate]) or
            !std.math.isFinite(inputs.soil.nitrogen_g_n[substrate]) or
            !std.math.isFinite(inputs.soil.phosphorus_g_p[substrate]))
            return error.InvalidRootOrganicExchangeInput;
    }
}

fn validateResult(result: Result) !void {
    inline for (0..substrate_count) |substrate| {
        inline for (.{
            result.carbon_exchange_g_c_per_step[substrate],
            result.nitrogen_exchange_g_n_per_step[substrate],
            result.phosphorus_exchange_g_p_per_step[substrate],
        }) |value|
            if (!std.math.isFinite(value))
                return error.NonFiniteRootOrganicExchangeResult;
    }
}

fn testInputs() Inputs {
    return .{
        .soil_micropore_water_m3 = 2,
        .substrate_fraction = .{0.2} ** substrate_count,
        .root_water_m3 = 1,
        .minimum_soil_water_m3 = 1e-12,
        .minimum_root_water_m3 = 1e-12,
        .minimum_soil_carbon_g_c = 1e-12,
        .minimum_root_carbon_g_c = 1e-12,
        .soil = .{
            .carbon_g_c = .{2} ** substrate_count,
            .nitrogen_g_n = .{0.4} ** substrate_count,
            .phosphorus_g_p = .{0.2} ** substrate_count,
        },
        .root = .{
            .carbon_g_c = 1,
            .nitrogen_g_n = 0.5,
            .phosphorus_g_p = 0.2,
        },
        .exchange_rate_per_h = 0.1,
        .solute_timestep_h = 0.5,
    };
}

test "five substrate exchanges preserve source-order equations" {
    const inputs = testInputs();
    const result = try compute(inputs);
    const substrate_water = 0.4;
    const expected_carbon =
        inputs.exchange_rate_per_h *
        ((2 * 1 - 1 * substrate_water) / (substrate_water + 1)) *
        inputs.solute_timestep_h;
    try std.testing.expectApproxEqAbs(
        expected_carbon,
        result.carbon_exchange_g_c_per_step[0],
        1e-12,
    );
    const expected_nitrogen: f64 = 0.1 * ((0.4 * 1 - 0.05 * 2) / 3.0) * 0.5;
    const expected_phosphorus: f64 = 0.1 * ((0.2 * 1 - 0.02 * 2) / 3.0) * 0.5;
    try std.testing.expectApproxEqAbs(
        expected_nitrogen,
        result.nitrogen_exchange_g_n_per_step[0],
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        expected_phosphorus,
        result.phosphorus_exchange_g_p_per_step[0],
        1e-12,
    );
    inline for (1..substrate_count) |substrate|
        try std.testing.expectApproxEqAbs(
            result.carbon_exchange_g_c_per_step[0],
            result.carbon_exchange_g_c_per_step[substrate],
            1e-12,
        );
}

test "root-rich carbon produces signed exudation" {
    var inputs = testInputs();
    inputs.root.carbon_g_c = 20;
    inputs.soil.carbon_g_c = .{0.1} ** substrate_count;
    const result = try compute(inputs);
    try std.testing.expect(result.carbon_exchange_g_c_per_step[0] < 0);
}

test "carbon threshold zeros nitrogen and phosphorus only" {
    var inputs = testInputs();
    inputs.soil.carbon_g_c = .{0} ** substrate_count;
    const result = try compute(inputs);
    try std.testing.expect(result.carbon_exchange_g_c_per_step[0] < 0);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.nitrogen_exchange_g_n_per_step[0],
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        result.phosphorus_exchange_g_p_per_step[0],
    );
}

test "later runtime-axis error leaves destination unchanged" {
    var inputs = [_]Inputs{ testInputs(), testInputs() };
    inputs[1].substrate_fraction[3] = std.math.nan(f64);
    var scratch: [2]Result = undefined;
    var destination = [_]Result{std.mem.zeroes(Result)} ** 2;
    destination[0].carbon_exchange_g_c_per_step[0] = 41;
    destination[1].carbon_exchange_g_c_per_step[0] = 42;
    try std.testing.expectError(
        error.InvalidRootOrganicExchangeInput,
        computeRuntimeAxes(&inputs, &scratch, &destination),
    );
    try std.testing.expectEqual(
        @as(f64, 41),
        destination[0].carbon_exchange_g_c_per_step[0],
    );
    try std.testing.expectEqual(
        @as(f64, 42),
        destination[1].carbon_exchange_g_c_per_step[0],
    );
}
