const std = @import("std");
const numerics = @import("../../core/numerics.zig");

pub const Inputs = struct {
    non_band_fraction: f64,
    soil_ammonium_concentration_g_n_m3: f64,
    minimum_ammonium_concentration_g_n_m3: f64,
    root_water_uptake_m3_per_step: f64,
    radial_diffusivity_m3_per_step: f64,
    maximum_uptake_g_n_m2_h: f64,
    root_surface_area_m2_per_plant: f64,
    relative_protein_concentration: f64,
    temperature_response: f64,
    carbon_uptake_constraint: f64,
    nitrogen_uptake_constraint: f64,
    oxygen_constraint: f64,
    solute_timestep_h: f64,
    half_saturation_g_n_m3: f64,
    population: f64,
    total_soil_water_m3: f64,
    soil_ammonium_mass_g_n: f64,
    competition_fraction: f64,
};

pub const SolverPolicy = struct {
    maximum_iterations: u16,
    absolute_tolerance_g_n_per_step: f64,
    relative_tolerance: f64,
    picard_relaxation: f64,
};

pub const Result = struct {
    active: bool,
    convective_flow_g_n_per_plant_step: f64,
    diffusivity_m3_per_step: f64,
    unconstrained_demand_g_n_per_plant_step: f64,
    oxygen_limited_demand_g_n_per_plant_step: f64,
    oxygen_limited_root_uptake_g_n_per_plant_step: f64,
    oxygen_unlimited_root_uptake_g_n_per_plant_step: f64,
    available_ammonium_g_n_per_step: f64,
    population_uptake_g_n_per_step: f64,
    oxygen_unlimited_population_uptake_g_n_per_step: f64,
    carbon_unlimited_population_uptake_g_n_per_step: f64,
    iterations: u32,
};

const Quadratic = struct { b: f64, c: f64 };

/// UPTAKE.F 2948--3029. The two closed-form roots at 2996 and 2999 are
/// solved by bounded local Newton-Raphson/Picard iteration.
pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    policy: SolverPolicy,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.RootAmmoniumUptakeDimensionMismatch;
    for (inputs, scratch) |axis_inputs, *candidate|
        candidate.* = try compute(axis_inputs, policy);
    @memcpy(destination, scratch);
}

pub fn compute(inputs: Inputs, policy: SolverPolicy) !Result {
    try validate(inputs, policy);
    if (inputs.non_band_fraction <= 0 or
        inputs.soil_ammonium_concentration_g_n_m3 <=
            inputs.minimum_ammonium_concentration_g_n_m3)
        return std.mem.zeroes(Result);
    const convective_flow =
        inputs.root_water_uptake_m3_per_step *
        inputs.soil_ammonium_concentration_g_n_m3 *
        inputs.non_band_fraction;
    const diffusivity =
        inputs.radial_diffusivity_m3_per_step * inputs.non_band_fraction;
    const unconstrained_demand =
        inputs.maximum_uptake_g_n_m2_h *
        inputs.root_surface_area_m2_per_plant *
        inputs.relative_protein_concentration *
        inputs.temperature_response *
        inputs.non_band_fraction *
        @min(inputs.carbon_uptake_constraint, inputs.nitrogen_uptake_constraint) *
        inputs.solute_timestep_h;
    const oxygen_limited_demand = unconstrained_demand * inputs.oxygen_constraint;
    const x = (diffusivity + convective_flow) *
        inputs.soil_ammonium_concentration_g_n_m3;
    const y = diffusivity * inputs.minimum_ammonium_concentration_g_n_m3;
    const limited = try solveRoot(
        oxygen_limited_demand,
        diffusivity,
        x,
        y,
        inputs.half_saturation_g_n_m3,
        policy,
    );
    const unlimited = try solveRoot(
        unconstrained_demand,
        diffusivity,
        x,
        y,
        inputs.half_saturation_g_n_m3,
        policy,
    );
    const minimum_mass =
        inputs.minimum_ammonium_concentration_g_n_m3 *
        inputs.total_soil_water_m3 *
        inputs.non_band_fraction;
    const available = @max(
        0,
        inputs.competition_fraction *
            (inputs.soil_ammonium_mass_g_n - minimum_mass) *
            inputs.solute_timestep_h,
    );
    const population_limited = @max(0, limited.root * inputs.population);
    const population_uptake = @min(available, population_limited);
    const population_unlimited = @min(
        available,
        @max(0, unlimited.root * inputs.population),
    );
    const result = Result{
        .active = true,
        .convective_flow_g_n_per_plant_step = convective_flow,
        .diffusivity_m3_per_step = diffusivity,
        .unconstrained_demand_g_n_per_plant_step = unconstrained_demand,
        .oxygen_limited_demand_g_n_per_plant_step = oxygen_limited_demand,
        .oxygen_limited_root_uptake_g_n_per_plant_step = limited.root,
        .oxygen_unlimited_root_uptake_g_n_per_plant_step = unlimited.root,
        .available_ammonium_g_n_per_step = available,
        .population_uptake_g_n_per_step = population_uptake,
        .oxygen_unlimited_population_uptake_g_n_per_step = population_unlimited,
        .carbon_unlimited_population_uptake_g_n_per_step = population_uptake / inputs.carbon_uptake_constraint,
        .iterations = @as(u32, limited.iterations) + unlimited.iterations,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteRootAmmoniumUptakeResult;
    }
    return result;
}

fn solveRoot(
    demand: f64,
    diffusivity: f64,
    x: f64,
    y: f64,
    half_saturation: f64,
    policy: SolverPolicy,
) !numerics.SolveResult {
    if (demand <= 0 or x - y <= 0) return .{
        .root = 0,
        .residual = 0,
        .iterations = 0,
        .newton_raphson_steps = 0,
        .picard_steps = 0,
    };
    const quadratic = Quadratic{
        .b = -demand - diffusivity * half_saturation - x + y,
        .c = (x - y) * demand,
    };
    const Functions = struct {
        fn residual(context: Quadratic, uptake: f64) f64 {
            return uptake * uptake + context.b * uptake + context.c;
        }
        fn derivative(context: Quadratic, uptake: f64) f64 {
            return 2 * uptake + context.b;
        }
        fn picard(context: Quadratic, uptake: f64) f64 {
            const denominator = -context.b - uptake;
            return if (denominator > 0) context.c / denominator else 0;
        }
    };
    const upper = @max(demand, x - y);
    return numerics.newtonPicard(
        quadratic,
        Functions.residual,
        Functions.derivative,
        Functions.picard,
        0,
        upper,
        0,
        .{
            .max_iterations = policy.maximum_iterations,
            .absolute_tolerance = policy.absolute_tolerance_g_n_per_step,
            .relative_tolerance = policy.relative_tolerance,
            .picard_relaxation = policy.picard_relaxation,
            .residual_scale = @max(upper * upper, 1e-30),
            .safeguard_with_bracket = true,
        },
    );
}

fn validate(inputs: Inputs, policy: SolverPolicy) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootAmmoniumUptakeInput;
    }
    if (inputs.non_band_fraction > 1 or
        inputs.carbon_uptake_constraint <= 0 or
        inputs.carbon_uptake_constraint > 1 or
        inputs.nitrogen_uptake_constraint > 1 or
        inputs.oxygen_constraint > 1 or
        inputs.competition_fraction > 1 or
        policy.maximum_iterations == 0 or
        !std.math.isFinite(policy.absolute_tolerance_g_n_per_step) or
        policy.absolute_tolerance_g_n_per_step <= 0 or
        !std.math.isFinite(policy.relative_tolerance) or
        policy.relative_tolerance <= 0 or
        !std.math.isFinite(policy.picard_relaxation) or
        policy.picard_relaxation <= 0 or
        policy.picard_relaxation > 1)
        return error.InvalidRootAmmoniumUptakeInput;
}

fn testInputs() Inputs {
    return .{
        .non_band_fraction = 0.7,
        .soil_ammonium_concentration_g_n_m3 = 2,
        .minimum_ammonium_concentration_g_n_m3 = 0.1,
        .root_water_uptake_m3_per_step = 0.1,
        .radial_diffusivity_m3_per_step = 0.5,
        .maximum_uptake_g_n_m2_h = 0.8,
        .root_surface_area_m2_per_plant = 1,
        .relative_protein_concentration = 0.9,
        .temperature_response = 0.8,
        .carbon_uptake_constraint = 0.7,
        .nitrogen_uptake_constraint = 0.6,
        .oxygen_constraint = 0.5,
        .solute_timestep_h = 0.25,
        .half_saturation_g_n_m3 = 0.2,
        .population = 2,
        .total_soil_water_m3 = 3,
        .soil_ammonium_mass_g_n = 4,
        .competition_fraction = 0.4,
    };
}

fn testPolicy() SolverPolicy {
    return .{
        .maximum_iterations = 30,
        .absolute_tolerance_g_n_per_step = 1e-12,
        .relative_tolerance = 1e-10,
        .picard_relaxation = 0.5,
    };
}

test "non-band ammonium solver converges and respects availability" {
    const result = try compute(testInputs(), testPolicy());
    try std.testing.expect(result.active);
    try std.testing.expect(result.iterations < 60);
    try std.testing.expect(result.population_uptake_g_n_per_step <=
        result.available_ammonium_g_n_per_step);
    try std.testing.expect(result.population_uptake_g_n_per_step <=
        result.oxygen_unlimited_population_uptake_g_n_per_step);
}

test "inactive concentration branch zeros every diagnostic" {
    var inputs = testInputs();
    inputs.soil_ammonium_concentration_g_n_m3 =
        inputs.minimum_ammonium_concentration_g_n_m3;
    try std.testing.expectEqualDeep(
        std.mem.zeroes(Result),
        try compute(inputs, testPolicy()),
    );
}

test "carbon-unlimited diagnostic preserves legacy division" {
    const inputs = testInputs();
    const result = try compute(inputs, testPolicy());
    try std.testing.expectApproxEqAbs(
        result.population_uptake_g_n_per_step / inputs.carbon_uptake_constraint,
        result.carbon_unlimited_population_uptake_g_n_per_step,
        1e-12,
    );
}

test "runtime axes fail atomically on later invalid population" {
    var inputs = [_]Inputs{ testInputs(), testInputs() };
    inputs[1].population = std.math.nan(f64);
    var scratch: [2]Result = undefined;
    var destination = [_]Result{std.mem.zeroes(Result)} ** 2;
    destination[0].population_uptake_g_n_per_step = 41;
    destination[1].population_uptake_g_n_per_step = 42;
    try std.testing.expectError(
        error.InvalidRootAmmoniumUptakeInput,
        computeRuntimeAxes(&inputs, testPolicy(), &scratch, &destination),
    );
    try std.testing.expectEqual(
        @as(f64, 41),
        destination[0].population_uptake_g_n_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 42),
        destination[1].population_uptake_g_n_per_step,
    );
}
