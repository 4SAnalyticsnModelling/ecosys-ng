const std = @import("std");
const nitrate = @import("root_nitrate_non_band_uptake.zig");

pub const SolverPolicy = nitrate.SolverPolicy;
pub const Result = nitrate.Result;

pub const Inputs = struct {
    band_fraction: f64,
    soil_nitrate_concentration_g_n_m3: f64,
    minimum_nitrate_concentration_g_n_m3: f64,
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
    band_nitrate_mass_g_n: f64,
    competition_fraction: f64,
};

/// UPTAKE.F 3232--3295 band NO3 branch.
pub fn compute(inputs: Inputs, policy: SolverPolicy) !Result {
    return nitrate.compute(asZoneInputs(inputs), policy);
}

pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    policy: SolverPolicy,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.RootNitrateBandUptakeDimensionMismatch;
    for (inputs, scratch) |axis_inputs, *candidate|
        candidate.* = try compute(axis_inputs, policy);
    @memcpy(destination, scratch);
}

fn asZoneInputs(inputs: Inputs) nitrate.Inputs {
    return .{
        .non_band_fraction = inputs.band_fraction,
        .soil_nitrate_concentration_g_n_m3 = inputs.soil_nitrate_concentration_g_n_m3,
        .minimum_nitrate_concentration_g_n_m3 = inputs.minimum_nitrate_concentration_g_n_m3,
        .root_water_uptake_m3_per_step = inputs.root_water_uptake_m3_per_step,
        .radial_diffusivity_m3_per_step = inputs.radial_diffusivity_m3_per_step,
        .maximum_uptake_g_n_m2_h = inputs.maximum_uptake_g_n_m2_h,
        .root_surface_area_m2_per_plant = inputs.root_surface_area_m2_per_plant,
        .relative_protein_concentration = inputs.relative_protein_concentration,
        .temperature_response = inputs.temperature_response,
        .carbon_uptake_constraint = inputs.carbon_uptake_constraint,
        .nitrogen_uptake_constraint = inputs.nitrogen_uptake_constraint,
        .oxygen_constraint = inputs.oxygen_constraint,
        .solute_timestep_h = inputs.solute_timestep_h,
        .half_saturation_g_n_m3 = inputs.half_saturation_g_n_m3,
        .population = inputs.population,
        .total_soil_water_m3 = inputs.total_soil_water_m3,
        .soil_nitrate_mass_g_n = inputs.band_nitrate_mass_g_n,
        .competition_fraction = inputs.competition_fraction,
    };
}

fn testInputs() Inputs {
    return .{
        .band_fraction = 0.35,
        .soil_nitrate_concentration_g_n_m3 = 2.5,
        .minimum_nitrate_concentration_g_n_m3 = 0.1,
        .root_water_uptake_m3_per_step = 0.1,
        .radial_diffusivity_m3_per_step = 0.4,
        .maximum_uptake_g_n_m2_h = 0.7,
        .root_surface_area_m2_per_plant = 1.2,
        .relative_protein_concentration = 0.9,
        .temperature_response = 0.8,
        .carbon_uptake_constraint = 0.7,
        .nitrogen_uptake_constraint = 0.6,
        .oxygen_constraint = 0.5,
        .solute_timestep_h = 0.25,
        .half_saturation_g_n_m3 = 0.2,
        .population = 2,
        .total_soil_water_m3 = 3,
        .band_nitrate_mass_g_n = 2,
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

test "band nitrate solver converges before both runtime ceilings" {
    const result = try compute(testInputs(), testPolicy());
    try std.testing.expect(result.active);
    try std.testing.expect(result.iterations < 60);
    try std.testing.expect(result.population_uptake_g_n_per_step <=
        result.available_nitrate_g_n_per_step);
    try std.testing.expect(result.population_uptake_g_n_per_step <=
        result.oxygen_unlimited_population_uptake_g_n_per_step);
}

test "band nitrate availability uses band fraction and pool" {
    const inputs = testInputs();
    const result = try compute(inputs, testPolicy());
    const minimum_mass =
        inputs.minimum_nitrate_concentration_g_n_m3 *
        inputs.total_soil_water_m3 *
        inputs.band_fraction;
    const expected = @max(
        0,
        inputs.competition_fraction *
            (inputs.band_nitrate_mass_g_n - minimum_mass) *
            inputs.solute_timestep_h,
    );
    try std.testing.expectApproxEqAbs(
        expected,
        result.available_nitrate_g_n_per_step,
        1e-12,
    );
}

test "band nitrate gate resets all diagnostics" {
    var inputs = testInputs();
    inputs.band_fraction = 0;
    try std.testing.expectEqualDeep(
        std.mem.zeroes(Result),
        try compute(inputs, testPolicy()),
    );
}

test "runtime band nitrate axes fail atomically" {
    var inputs = [_]Inputs{ testInputs(), testInputs() };
    inputs[1].band_fraction = 2;
    var scratch: [2]Result = undefined;
    var destination = [_]Result{std.mem.zeroes(Result)} ** 2;
    destination[0].population_uptake_g_n_per_step = 41;
    destination[1].population_uptake_g_n_per_step = 42;
    try std.testing.expectError(
        error.InvalidRootNitrateUptakeInput,
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
