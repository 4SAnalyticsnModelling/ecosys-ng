const std = @import("std");

pub const Inputs = struct {
    oxygen_constraint: f64,
    carbon_uptake_constraint: f64,
    relative_protein_concentration: f64,
    root_length_m_per_plant: f64,
    negligible_root_length_m_per_plant: f64,
    nitrogen_uptake_constraint: f64,
    minimum_nitrogen_uptake_constraint: f64,
    ammonium_diffusivity_m2_per_step: f64,
    soil_tortuosity: f64,
    solute_timestep_h: f64,
    water_nutrient_path_m: f64,
    root_radius_m: f64,
    root_surface_area_per_radius_m: f64,
};

pub const Result = struct {
    nutrient_uptake_active: bool,
    nitrogen_uptake_active: bool,
    tortuosity_adjusted_diffusivity_m2_per_step: f64,
    limited_diffusion_path_m: f64,
    radial_ammonium_diffusivity_m3_per_step: f64,
};

/// UPTAKE.F 2913--2934 admission gates and NH4/NO3 radial transport setup.
pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.RootNutrientAdmissionDimensionMismatch;
    for (inputs, scratch) |axis_inputs, *candidate|
        candidate.* = try compute(axis_inputs);
    @memcpy(destination, scratch);
}

pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    const nutrient_active =
        inputs.oxygen_constraint > 0 and
        inputs.carbon_uptake_constraint > 0 and
        inputs.relative_protein_concentration > 0 and
        inputs.root_length_m_per_plant >
            inputs.negligible_root_length_m_per_plant;
    const nitrogen_active = nutrient_active and
        inputs.nitrogen_uptake_constraint >
            inputs.minimum_nitrogen_uptake_constraint;
    if (!nitrogen_active) return .{
        .nutrient_uptake_active = nutrient_active,
        .nitrogen_uptake_active = false,
        .tortuosity_adjusted_diffusivity_m2_per_step = 0,
        .limited_diffusion_path_m = 0,
        .radial_ammonium_diffusivity_m3_per_step = 0,
    };
    const adjusted_diffusivity =
        inputs.ammonium_diffusivity_m2_per_step *
        inputs.soil_tortuosity *
        inputs.solute_timestep_h;
    if (adjusted_diffusivity <= 0)
        return error.DegenerateRootNutrientDiffusivity;
    const limited_path = @min(
        inputs.water_nutrient_path_m,
        @sqrt(2 * adjusted_diffusivity),
    );
    const path_log = @log(
        (limited_path + inputs.root_radius_m) / inputs.root_radius_m,
    );
    if (!std.math.isFinite(path_log) or path_log <= 0)
        return error.DegenerateRootNutrientDiffusionPath;
    const radial_diffusivity =
        adjusted_diffusivity *
        inputs.root_surface_area_per_radius_m /
        path_log;
    if (!std.math.isFinite(radial_diffusivity))
        return error.NonFiniteRootNutrientAdmissionResult;
    return .{
        .nutrient_uptake_active = true,
        .nitrogen_uptake_active = true,
        .tortuosity_adjusted_diffusivity_m2_per_step = adjusted_diffusivity,
        .limited_diffusion_path_m = limited_path,
        .radial_ammonium_diffusivity_m3_per_step = radial_diffusivity,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootNutrientAdmissionInput;
    }
    if (inputs.oxygen_constraint > 1 or
        inputs.carbon_uptake_constraint > 1 or
        inputs.nitrogen_uptake_constraint > 1 or
        inputs.root_radius_m <= 0)
        return error.InvalidRootNutrientAdmissionInput;
}

fn testInputs() Inputs {
    return .{
        .oxygen_constraint = 0.8,
        .carbon_uptake_constraint = 0.7,
        .relative_protein_concentration = 0.9,
        .root_length_m_per_plant = 2,
        .negligible_root_length_m_per_plant = 1e-12,
        .nitrogen_uptake_constraint = 0.6,
        .minimum_nitrogen_uptake_constraint = 1e-12,
        .ammonium_diffusivity_m2_per_step = 0.5,
        .soil_tortuosity = 0.4,
        .solute_timestep_h = 0.5,
        .water_nutrient_path_m = 1,
        .root_radius_m = 0.01,
        .root_surface_area_per_radius_m = 2,
    };
}

test "active nitrogen gate preserves radial transport equations" {
    const inputs = testInputs();
    const result = try compute(inputs);
    const adjusted = 0.5 * 0.4 * 0.5;
    const path = @min(@as(f64, 1), @sqrt(2 * adjusted));
    const expected = adjusted * 2 / @log((path + 0.01) / 0.01);
    try std.testing.expect(result.nutrient_uptake_active);
    try std.testing.expect(result.nitrogen_uptake_active);
    try std.testing.expectApproxEqAbs(
        adjusted,
        result.tortuosity_adjusted_diffusivity_m2_per_step,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(path, result.limited_diffusion_path_m, 1e-12);
    try std.testing.expectApproxEqAbs(
        expected,
        result.radial_ammonium_diffusivity_m3_per_step,
        1e-12,
    );
}

test "outer nutrient gate failure zeros nitrogen transport" {
    var inputs = testInputs();
    inputs.oxygen_constraint = 0;
    const result = try compute(inputs);
    try std.testing.expect(!result.nutrient_uptake_active);
    try std.testing.expect(!result.nitrogen_uptake_active);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.radial_ammonium_diffusivity_m3_per_step,
    );
}

test "nitrogen gate can fail while outer nutrient gate remains active" {
    var inputs = testInputs();
    inputs.nitrogen_uptake_constraint = 0;
    const result = try compute(inputs);
    try std.testing.expect(result.nutrient_uptake_active);
    try std.testing.expect(!result.nitrogen_uptake_active);
}

test "active degenerate diffusion fails instead of producing NaN" {
    var inputs = testInputs();
    inputs.soil_tortuosity = 0;
    try std.testing.expectError(
        error.DegenerateRootNutrientDiffusivity,
        compute(inputs),
    );
}

test "runtime axes fail atomically on later invalid radius" {
    var inputs = [_]Inputs{ testInputs(), testInputs() };
    inputs[1].root_radius_m = 0;
    var scratch: [2]Result = undefined;
    var destination = [_]Result{std.mem.zeroes(Result)} ** 2;
    destination[0].radial_ammonium_diffusivity_m3_per_step = 41;
    destination[1].radial_ammonium_diffusivity_m3_per_step = 42;
    try std.testing.expectError(
        error.InvalidRootNutrientAdmissionInput,
        computeRuntimeAxes(&inputs, &scratch, &destination),
    );
    try std.testing.expectEqual(
        @as(f64, 41),
        destination[0].radial_ammonium_diffusivity_m3_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 42),
        destination[1].radial_ammonium_diffusivity_m3_per_step,
    );
}
