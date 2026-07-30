const std = @import("std");

pub const Inputs = struct {
    nitrate_diffusivity_m2_per_step: f64,
    soil_tortuosity: f64,
    solute_timestep_h: f64,
    water_nutrient_path_m: f64,
    root_radius_m: f64,
    root_surface_area_per_radius_m: f64,
};

pub const Result = struct {
    tortuosity_adjusted_diffusivity_m2_per_step: f64,
    limited_diffusion_path_m: f64,
    radial_nitrate_diffusivity_m3_per_step: f64,
};

/// UPTAKE.F 3126--3128 radial nitrate transport preparation.
pub fn computeRuntimeAxes(
    inputs: []const Inputs,
    scratch: []Result,
    destination: []Result,
) !void {
    if (inputs.len != scratch.len or inputs.len != destination.len)
        return error.RootNitrateTransportDimensionMismatch;
    for (inputs, scratch) |axis_inputs, *candidate|
        candidate.* = try compute(axis_inputs);
    @memcpy(destination, scratch);
}

pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    const adjusted_diffusivity =
        inputs.nitrate_diffusivity_m2_per_step *
        inputs.soil_tortuosity *
        inputs.solute_timestep_h;
    if (adjusted_diffusivity <= 0)
        return error.DegenerateRootNitrateDiffusivity;
    const limited_path = @min(
        inputs.water_nutrient_path_m,
        @sqrt(2 * adjusted_diffusivity),
    );
    const path_log = @log(
        (limited_path + inputs.root_radius_m) / inputs.root_radius_m,
    );
    if (!std.math.isFinite(path_log) or path_log <= 0)
        return error.DegenerateRootNitrateDiffusionPath;
    const radial_diffusivity =
        adjusted_diffusivity *
        inputs.root_surface_area_per_radius_m /
        path_log;
    if (!std.math.isFinite(radial_diffusivity))
        return error.NonFiniteRootNitrateTransportResult;
    return .{
        .tortuosity_adjusted_diffusivity_m2_per_step = adjusted_diffusivity,
        .limited_diffusion_path_m = limited_path,
        .radial_nitrate_diffusivity_m3_per_step = radial_diffusivity,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootNitrateTransportInput;
    }
    if (inputs.root_radius_m <= 0)
        return error.InvalidRootNitrateTransportInput;
}

fn testInputs() Inputs {
    return .{
        .nitrate_diffusivity_m2_per_step = 0.4,
        .soil_tortuosity = 0.5,
        .solute_timestep_h = 0.25,
        .water_nutrient_path_m = 1,
        .root_radius_m = 0.01,
        .root_surface_area_per_radius_m = 2,
    };
}

test "nitrate transport preserves adjusted path and radial equations" {
    const result = try compute(testInputs());
    const adjusted: f64 = 0.4 * 0.5 * 0.25;
    const path = @min(@as(f64, 1), @sqrt(2 * adjusted));
    const expected = adjusted * 2 / @log((path + 0.01) / 0.01);
    try std.testing.expectApproxEqAbs(
        adjusted,
        result.tortuosity_adjusted_diffusivity_m2_per_step,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(path, result.limited_diffusion_path_m, 1e-12);
    try std.testing.expectApproxEqAbs(
        expected,
        result.radial_nitrate_diffusivity_m3_per_step,
        1e-12,
    );
}

test "model path limits diffusion before the diffusive distance" {
    var inputs = testInputs();
    inputs.water_nutrient_path_m = 0.02;
    const result = try compute(inputs);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.02),
        result.limited_diffusion_path_m,
        1e-12,
    );
}

test "degenerate active nitrate diffusion fails explicitly" {
    var inputs = testInputs();
    inputs.soil_tortuosity = 0;
    try std.testing.expectError(
        error.DegenerateRootNitrateDiffusivity,
        compute(inputs),
    );
}

test "runtime axes fail atomically on a later invalid radius" {
    var inputs = [_]Inputs{ testInputs(), testInputs() };
    inputs[1].root_radius_m = 0;
    var scratch: [2]Result = undefined;
    var destination = [_]Result{std.mem.zeroes(Result)} ** 2;
    destination[0].radial_nitrate_diffusivity_m3_per_step = 41;
    destination[1].radial_nitrate_diffusivity_m3_per_step = 42;
    try std.testing.expectError(
        error.InvalidRootNitrateTransportInput,
        computeRuntimeAxes(&inputs, &scratch, &destination),
    );
    try std.testing.expectEqual(
        @as(f64, 41),
        destination[0].radial_nitrate_diffusivity_m3_per_step,
    );
    try std.testing.expectEqual(
        @as(f64, 42),
        destination[1].radial_nitrate_diffusivity_m3_per_step,
    );
}
