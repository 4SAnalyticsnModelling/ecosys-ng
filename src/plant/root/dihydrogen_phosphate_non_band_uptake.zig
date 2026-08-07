const std = @import("std");
const zone_solver = @import("ammonium_non_band_uptake.zig");

pub const SolverPolicy = struct {
    maximum_iterations: u16,
    absolute_tolerance_g_p_per_step: f64,
    relative_tolerance: f64,
    picard_relaxation: f64,
};

pub const Inputs = struct {
    non_band_fraction: f64,
    soil_concentration_g_p_m3: f64,
    minimum_concentration_g_p_m3: f64,
    root_water_uptake_m3_per_step: f64,
    radial_diffusivity_m3_per_step: f64,
    maximum_uptake_g_p_m2_h: f64,
    root_surface_area_m2_per_plant: f64,
    relative_protein_concentration: f64,
    temperature_response: f64,
    carbon_uptake_constraint: f64,
    phosphorus_uptake_constraint: f64,
    oxygen_constraint: f64,
    solute_timestep_h: f64,
    half_saturation_g_p_m3: f64,
    population: f64,
    total_soil_water_m3: f64,
    soil_mass_g_p: f64,
    competition_fraction: f64,
};

pub const Result = struct {
    active: bool,
    convective_flow_g_p_per_plant_step: f64,
    diffusivity_m3_per_step: f64,
    unconstrained_demand_g_p_per_plant_step: f64,
    oxygen_limited_demand_g_p_per_plant_step: f64,
    oxygen_limited_root_uptake_g_p_per_plant_step: f64,
    oxygen_unlimited_root_uptake_g_p_per_plant_step: f64,
    available_phosphorus_g_p_per_step: f64,
    population_uptake_g_p_per_step: f64,
    oxygen_unlimited_population_uptake_g_p_per_step: f64,
    carbon_unlimited_population_uptake_g_p_per_step: f64,
    iterations: u32,
};

/// UPTAKE.F 3346--3421 non-band H2PO4 uptake.
pub fn compute(inputs: Inputs, policy: SolverPolicy) !Result {
    const solved = zone_solver.compute(asZoneInputs(inputs), asZonePolicy(policy)) catch |err|
        return switch (err) {
            error.InvalidRootAmmoniumUptakeInput => error.InvalidRootDihydrogenPhosphateUptakeInput,
            else => err,
        };
    return fromZoneResult(solved);
}

fn asZonePolicy(policy: SolverPolicy) zone_solver.SolverPolicy {
    return .{
        .maximum_iterations = policy.maximum_iterations,
        .absolute_tolerance_g_n_per_step = policy.absolute_tolerance_g_p_per_step,
        .relative_tolerance = policy.relative_tolerance,
        .picard_relaxation = policy.picard_relaxation,
    };
}

pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    policy: SolverPolicy,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.RootDihydrogenPhosphateUptakeDimensionMismatch;
    for (inputs, scratch) |axis_inputs, *candidate|
        candidate.* = try compute(axis_inputs, policy);
    @memcpy(destination, scratch);
}

fn asZoneInputs(inputs: Inputs) zone_solver.Inputs {
    return .{
        .non_band_fraction = inputs.non_band_fraction,
        .soil_ammonium_concentration_g_n_m3 = inputs.soil_concentration_g_p_m3,
        .minimum_ammonium_concentration_g_n_m3 = inputs.minimum_concentration_g_p_m3,
        .root_water_uptake_m3_per_step = inputs.root_water_uptake_m3_per_step,
        .radial_diffusivity_m3_per_step = inputs.radial_diffusivity_m3_per_step,
        .maximum_uptake_g_n_m2_h = inputs.maximum_uptake_g_p_m2_h,
        .root_surface_area_m2_per_plant = inputs.root_surface_area_m2_per_plant,
        .relative_protein_concentration = inputs.relative_protein_concentration,
        .temperature_response = inputs.temperature_response,
        .carbon_uptake_constraint = inputs.carbon_uptake_constraint,
        .nitrogen_uptake_constraint = inputs.phosphorus_uptake_constraint,
        .oxygen_constraint = inputs.oxygen_constraint,
        .solute_timestep_h = inputs.solute_timestep_h,
        .half_saturation_g_n_m3 = inputs.half_saturation_g_p_m3,
        .population = inputs.population,
        .total_soil_water_m3 = inputs.total_soil_water_m3,
        .soil_ammonium_mass_g_n = inputs.soil_mass_g_p,
        .competition_fraction = inputs.competition_fraction,
    };
}

fn fromZoneResult(result: zone_solver.Result) Result {
    return .{
        .active = result.active,
        .convective_flow_g_p_per_plant_step = result.convective_flow_g_n_per_plant_step,
        .diffusivity_m3_per_step = result.diffusivity_m3_per_step,
        .unconstrained_demand_g_p_per_plant_step = result.unconstrained_demand_g_n_per_plant_step,
        .oxygen_limited_demand_g_p_per_plant_step = result.oxygen_limited_demand_g_n_per_plant_step,
        .oxygen_limited_root_uptake_g_p_per_plant_step = result.oxygen_limited_root_uptake_g_n_per_plant_step,
        .oxygen_unlimited_root_uptake_g_p_per_plant_step = result.oxygen_unlimited_root_uptake_g_n_per_plant_step,
        .available_phosphorus_g_p_per_step = result.available_ammonium_g_n_per_step,
        .population_uptake_g_p_per_step = result.population_uptake_g_n_per_step,
        .oxygen_unlimited_population_uptake_g_p_per_step = result.oxygen_unlimited_population_uptake_g_n_per_step,
        .carbon_unlimited_population_uptake_g_p_per_step = result.carbon_unlimited_population_uptake_g_n_per_step,
        .iterations = result.iterations,
    };
}

fn testInputs() Inputs {
    return .{
        .non_band_fraction = 0.7,
        .soil_concentration_g_p_m3 = 1.5,
        .minimum_concentration_g_p_m3 = 0.05,
        .root_water_uptake_m3_per_step = 0.1,
        .radial_diffusivity_m3_per_step = 0.3,
        .maximum_uptake_g_p_m2_h = 0.4,
        .root_surface_area_m2_per_plant = 1.2,
        .relative_protein_concentration = 0.9,
        .temperature_response = 0.8,
        .carbon_uptake_constraint = 0.7,
        .phosphorus_uptake_constraint = 0.6,
        .oxygen_constraint = 0.5,
        .solute_timestep_h = 0.25,
        .half_saturation_g_p_m3 = 0.1,
        .population = 2,
        .total_soil_water_m3 = 3,
        .soil_mass_g_p = 2,
        .competition_fraction = 0.4,
    };
}

fn testPolicy() SolverPolicy {
    return .{
        .maximum_iterations = 30,
        .absolute_tolerance_g_p_per_step = 1e-12,
        .relative_tolerance = 1e-10,
        .picard_relaxation = 0.5,
    };
}

test "non-band H2PO4 solve converges and preserves uptake bounds" {
    const result = try compute(testInputs(), testPolicy());
    try std.testing.expect(result.active);
    try std.testing.expect(result.iterations < 60);
    try std.testing.expect(result.population_uptake_g_p_per_step <=
        result.available_phosphorus_g_p_per_step);
    try std.testing.expect(result.population_uptake_g_p_per_step <=
        result.oxygen_unlimited_population_uptake_g_p_per_step);
}

test "H2PO4 availability preserves retained phosphorus mass" {
    const inputs = testInputs();
    const result = try compute(inputs, testPolicy());
    const retained =
        inputs.minimum_concentration_g_p_m3 *
        inputs.total_soil_water_m3 *
        inputs.non_band_fraction;
    const expected = @max(
        0,
        inputs.competition_fraction *
            (inputs.soil_mass_g_p - retained) *
            inputs.solute_timestep_h,
    );
    try std.testing.expectApproxEqAbs(
        expected,
        result.available_phosphorus_g_p_per_step,
        1e-12,
    );
}

test "minimum H2PO4 concentration gate resets diagnostics" {
    var inputs = testInputs();
    inputs.soil_concentration_g_p_m3 = inputs.minimum_concentration_g_p_m3;
    try std.testing.expectEqualDeep(
        std.mem.zeroes(Result),
        try compute(inputs, testPolicy()),
    );
}

test "runtime H2PO4 axes fail atomically on later invalid fraction" {
    var inputs = [_]Inputs{ testInputs(), testInputs() };
    inputs[1].non_band_fraction = 2;
    var scratch: [2]Result = undefined;
    var destination = [_]Result{std.mem.zeroes(Result)} ** 2;
    destination[0].population_uptake_g_p_per_step = 41;
    destination[1].population_uptake_g_p_per_step = 42;
    try std.testing.expectError(
        error.InvalidRootDihydrogenPhosphateUptakeInput,
        computeRuntimeAxes(&inputs, testPolicy(), &scratch, &destination),
    );
    try std.testing.expectEqual(
        @as(f64, 41),
        destination[0].population_uptake_g_p_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 42),
        destination[1].population_uptake_g_p_per_step,
    );
}
