const std = @import("std");
const phosphate = @import("dihydrogen_phosphate_non_band_uptake.zig");

pub const SolverPolicy = phosphate.SolverPolicy;
pub const Result = phosphate.Result;

pub const Inputs = struct {
    band_fraction: f64,
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
    band_mass_g_p: f64,
    competition_fraction: f64,
};

/// UPTAKE.F 3435--3509 band H2PO4 uptake.
pub fn compute(inputs: Inputs, policy: SolverPolicy) !Result {
    return phosphate.compute(asZoneInputs(inputs), policy);
}

pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    policy: SolverPolicy,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.RootDihydrogenPhosphateBandDimensionMismatch;
    for (inputs, scratch) |axis_inputs, *candidate|
        candidate.* = try compute(axis_inputs, policy);
    @memcpy(destination, scratch);
}

fn asZoneInputs(inputs: Inputs) phosphate.Inputs {
    return .{
        .non_band_fraction = inputs.band_fraction,
        .soil_concentration_g_p_m3 = inputs.soil_concentration_g_p_m3,
        .minimum_concentration_g_p_m3 = inputs.minimum_concentration_g_p_m3,
        .root_water_uptake_m3_per_step = inputs.root_water_uptake_m3_per_step,
        .radial_diffusivity_m3_per_step = inputs.radial_diffusivity_m3_per_step,
        .maximum_uptake_g_p_m2_h = inputs.maximum_uptake_g_p_m2_h,
        .root_surface_area_m2_per_plant = inputs.root_surface_area_m2_per_plant,
        .relative_protein_concentration = inputs.relative_protein_concentration,
        .temperature_response = inputs.temperature_response,
        .carbon_uptake_constraint = inputs.carbon_uptake_constraint,
        .phosphorus_uptake_constraint = inputs.phosphorus_uptake_constraint,
        .oxygen_constraint = inputs.oxygen_constraint,
        .solute_timestep_h = inputs.solute_timestep_h,
        .half_saturation_g_p_m3 = inputs.half_saturation_g_p_m3,
        .population = inputs.population,
        .total_soil_water_m3 = inputs.total_soil_water_m3,
        .soil_mass_g_p = inputs.band_mass_g_p,
        .competition_fraction = inputs.competition_fraction,
    };
}

fn testInputs() Inputs {
    return .{
        .band_fraction = 0.3,
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
        .band_mass_g_p = 1.5,
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

test "band H2PO4 solver converges before both runtime ceilings" {
    const result = try compute(testInputs(), testPolicy());
    try std.testing.expect(result.active);
    try std.testing.expect(result.iterations < 60);
    try std.testing.expect(result.population_uptake_g_p_per_step <=
        result.available_phosphorus_g_p_per_step);
    try std.testing.expect(result.population_uptake_g_p_per_step <=
        result.oxygen_unlimited_population_uptake_g_p_per_step);
}

test "band H2PO4 availability uses band fraction and mass" {
    const inputs = testInputs();
    const result = try compute(inputs, testPolicy());
    const retained =
        inputs.minimum_concentration_g_p_m3 *
        inputs.total_soil_water_m3 *
        inputs.band_fraction;
    const expected = @max(
        0,
        inputs.competition_fraction *
            (inputs.band_mass_g_p - retained) *
            inputs.solute_timestep_h,
    );
    try std.testing.expectApproxEqAbs(
        expected,
        result.available_phosphorus_g_p_per_step,
        1e-12,
    );
}

test "band H2PO4 gate resets diagnostics" {
    var inputs = testInputs();
    inputs.band_fraction = 0;
    try std.testing.expectEqualDeep(
        std.mem.zeroes(Result),
        try compute(inputs, testPolicy()),
    );
}

test "runtime band H2PO4 axes fail atomically" {
    var inputs = [_]Inputs{ testInputs(), testInputs() };
    inputs[1].band_fraction = 2;
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
