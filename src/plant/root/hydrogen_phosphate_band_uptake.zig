const std = @import("std");
const hpo4 = @import("hydrogen_phosphate_non_band_uptake.zig");
const phosphate = @import("dihydrogen_phosphate_non_band_uptake.zig");

pub const SolverPolicy = phosphate.SolverPolicy;
pub const Result = phosphate.Result;

pub const Inputs = struct {
    kinetics_from_dihydrogen_phosphate: hpo4.AliasedKinetics,
    band_fraction: f64,
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
    band_mass_g_p: f64,
    competition_fraction: f64,
    /// UPTAKE line 3664 uses RMFH2B rather than the RMFH1B calculated at 3620.
    dihydrogen_phosphate_convective_flow_g_p_per_plant_step: f64,
};

/// UPTAKE.F 3618--3693 band HPO4 uptake, including the source's line-3664
/// dependency on band H2PO4 convective flow.
pub fn compute(inputs: Inputs, policy: SolverPolicy) !Result {
    const active = inputs.band_fraction > 0 and
        inputs.soil_concentration_g_p_m3 >
            inputs.kinetics_from_dihydrogen_phosphate.minimum_concentration_g_p_m3;
    var mapped = asZoneInputs(inputs);
    const own_convective_flow =
        inputs.root_water_uptake_m3_per_step *
        inputs.soil_concentration_g_p_m3 *
        inputs.band_fraction;
    if (active) {
        const concentration_fraction =
            inputs.soil_concentration_g_p_m3 * inputs.band_fraction;
        if (concentration_fraction <= 0)
            return error.InvalidRootHydrogenPhosphateBandUptakeInput;
        mapped.root_water_uptake_m3_per_step =
            inputs.dihydrogen_phosphate_convective_flow_g_p_per_plant_step /
            concentration_fraction;
    }
    var result = phosphate.compute(mapped, policy) catch |err|
        return switch (err) {
            error.InvalidRootDihydrogenPhosphateUptakeInput => error.InvalidRootHydrogenPhosphateBandUptakeInput,
            else => err,
        };
    // Report RMFH1B from line 3620 even though line 3664 uses RMFH2B.
    result.convective_flow_g_p_per_plant_step =
        if (active) own_convective_flow else 0;
    return result;
}

pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    policy: SolverPolicy,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.RootHydrogenPhosphateBandDimensionMismatch;
    for (inputs, scratch) |axis_inputs, *candidate|
        candidate.* = try compute(axis_inputs, policy);
    @memcpy(destination, scratch);
}

fn asZoneInputs(inputs: Inputs) phosphate.Inputs {
    const kinetics = inputs.kinetics_from_dihydrogen_phosphate;
    return .{
        .non_band_fraction = inputs.band_fraction,
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
        .soil_mass_g_p = inputs.band_mass_g_p,
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
        .band_fraction = 0.3,
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
        .band_mass_g_p = 1.5,
        .competition_fraction = 0.4,
        .dihydrogen_phosphate_convective_flow_g_p_per_plant_step = 0.2,
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

test "band HPO4 preserves distinct reported and nonlinear convective flows" {
    const inputs = testInputs();
    const result = try compute(inputs, testPolicy());
    const own_flow =
        inputs.root_water_uptake_m3_per_step *
        inputs.soil_concentration_g_p_m3 *
        inputs.band_fraction;
    try std.testing.expectApproxEqAbs(
        own_flow,
        result.convective_flow_g_p_per_plant_step,
        1e-12,
    );
    try std.testing.expect(result.iterations < 60);
}

test "line 3664 H2PO4 convective input changes HPO4 nonlinear uptake" {
    const baseline = try compute(testInputs(), testPolicy());
    var changed = testInputs();
    changed.dihydrogen_phosphate_convective_flow_g_p_per_plant_step = 0.8;
    const result = try compute(changed, testPolicy());
    try std.testing.expect(
        result.oxygen_limited_root_uptake_g_p_per_plant_step !=
            baseline.oxygen_limited_root_uptake_g_p_per_plant_step,
    );
}

test "band HPO4 gate resets diagnostics" {
    var inputs = testInputs();
    inputs.band_fraction = 0;
    try std.testing.expectEqualDeep(
        std.mem.zeroes(Result),
        try compute(inputs, testPolicy()),
    );
}

test "runtime band HPO4 axes fail atomically" {
    var inputs = [_]Inputs{ testInputs(), testInputs() };
    inputs[1].dihydrogen_phosphate_convective_flow_g_p_per_plant_step =
        std.math.nan(f64);
    var scratch: [2]Result = undefined;
    var destination = [_]Result{std.mem.zeroes(Result)} ** 2;
    destination[0].population_uptake_g_p_per_step = 41;
    destination[1].population_uptake_g_p_per_step = 42;
    try std.testing.expectError(
        error.InvalidRootHydrogenPhosphateBandUptakeInput,
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
