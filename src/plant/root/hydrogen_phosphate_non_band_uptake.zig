const std = @import("std");
const phosphate = @import("dihydrogen_phosphate_non_band_uptake.zig");

pub const SolverPolicy = phosphate.SolverPolicy;
pub const Result = phosphate.Result;

pub const AliasedKinetics = struct {
    /// UPTAKE 3524: HPO4 maximum uptake aliases the H2PO4 parameter.
    maximum_uptake_g_p_m2_h: f64,
    /// UPTAKE 3525: HPO4 minimum concentration aliases H2PO4.
    minimum_concentration_g_p_m3: f64,
    /// UPTAKE 3526: HPO4 half saturation aliases H2PO4.
    half_saturation_g_p_m3: f64,
};

pub const Inputs = struct {
    kinetics_from_dihydrogen_phosphate: AliasedKinetics,
    non_band_fraction: f64,
    soil_concentration_g_p_m3: f64,
    root_water_uptake_m3_per_step: f64,
    radial_diffusivity_m3_per_step: f64,
    root_surface_area_m2_per_plant: f64,
    relative_protein_concentration: f64,
    temperature_response: f64,
    carbon_uptake_constraint: f64,
    phosphorus_uptake_constraint: f64,
    oxygen_constraint: f64,
    solute_timestep_h: f64,
    population: f64,
    total_soil_water_m3: f64,
    soil_mass_g_p: f64,
    competition_fraction: f64,
};

/// UPTAKE.F 3524--3603 parameter aliases and non-band HPO4 uptake.
pub fn compute(inputs: Inputs, policy: SolverPolicy) !Result {
    return phosphate.compute(asZoneInputs(inputs), policy) catch |err|
        return switch (err) {
            error.InvalidRootDihydrogenPhosphateUptakeInput => error.InvalidRootHydrogenPhosphateUptakeInput,
            else => err,
        };
}

pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    policy: SolverPolicy,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.RootHydrogenPhosphateUptakeDimensionMismatch;
    for (inputs, scratch) |axis_inputs, *candidate|
        candidate.* = try compute(axis_inputs, policy);
    @memcpy(destination, scratch);
}

fn asZoneInputs(inputs: Inputs) phosphate.Inputs {
    const kinetics = inputs.kinetics_from_dihydrogen_phosphate;
    return .{
        .non_band_fraction = inputs.non_band_fraction,
        .soil_concentration_g_p_m3 = inputs.soil_concentration_g_p_m3,
        .minimum_concentration_g_p_m3 = kinetics.minimum_concentration_g_p_m3,
        .root_water_uptake_m3_per_step = inputs.root_water_uptake_m3_per_step,
        .radial_diffusivity_m3_per_step = inputs.radial_diffusivity_m3_per_step,
        .maximum_uptake_g_p_m2_h = kinetics.maximum_uptake_g_p_m2_h,
        .root_surface_area_m2_per_plant = inputs.root_surface_area_m2_per_plant,
        .relative_protein_concentration = inputs.relative_protein_concentration,
        .temperature_response = inputs.temperature_response,
        .carbon_uptake_constraint = inputs.carbon_uptake_constraint,
        .phosphorus_uptake_constraint = inputs.phosphorus_uptake_constraint,
        .oxygen_constraint = inputs.oxygen_constraint,
        .solute_timestep_h = inputs.solute_timestep_h,
        .half_saturation_g_p_m3 = kinetics.half_saturation_g_p_m3,
        .population = inputs.population,
        .total_soil_water_m3 = inputs.total_soil_water_m3,
        .soil_mass_g_p = inputs.soil_mass_g_p,
        .competition_fraction = inputs.competition_fraction,
    };
}

fn testInputs() Inputs {
    return .{
        .kinetics_from_dihydrogen_phosphate = .{
            .maximum_uptake_g_p_m2_h = 0.4,
            .minimum_concentration_g_p_m3 = 0.05,
            .half_saturation_g_p_m3 = 0.1,
        },
        .non_band_fraction = 0.7,
        .soil_concentration_g_p_m3 = 1.5,
        .root_water_uptake_m3_per_step = 0.1,
        .radial_diffusivity_m3_per_step = 0.3,
        .root_surface_area_m2_per_plant = 1.2,
        .relative_protein_concentration = 0.9,
        .temperature_response = 0.8,
        .carbon_uptake_constraint = 0.7,
        .phosphorus_uptake_constraint = 0.6,
        .oxygen_constraint = 0.5,
        .solute_timestep_h = 0.25,
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

test "HPO4 aliases reproduce the matching H2PO4 zone solution" {
    const inputs = testInputs();
    const hpo4 = try compute(inputs, testPolicy());
    const h2po4 = try phosphate.compute(asZoneInputs(inputs), testPolicy());
    try std.testing.expectEqualDeep(h2po4, hpo4);
    try std.testing.expect(hpo4.iterations < 60);
}

test "non-band HPO4 retained mass controls availability" {
    const inputs = testInputs();
    const result = try compute(inputs, testPolicy());
    const retained =
        inputs.kinetics_from_dihydrogen_phosphate.minimum_concentration_g_p_m3 *
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

test "HPO4 concentration gate resets diagnostics" {
    var inputs = testInputs();
    inputs.soil_concentration_g_p_m3 =
        inputs.kinetics_from_dihydrogen_phosphate.minimum_concentration_g_p_m3;
    try std.testing.expectEqualDeep(
        std.mem.zeroes(Result),
        try compute(inputs, testPolicy()),
    );
}

test "runtime HPO4 axes fail atomically on later invalid fraction" {
    var inputs = [_]Inputs{ testInputs(), testInputs() };
    inputs[1].non_band_fraction = 2;
    var scratch: [2]Result = undefined;
    var destination = [_]Result{std.mem.zeroes(Result)} ** 2;
    destination[0].population_uptake_g_p_per_step = 41;
    destination[1].population_uptake_g_p_per_step = 42;
    try std.testing.expectError(
        error.InvalidRootHydrogenPhosphateUptakeInput,
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
