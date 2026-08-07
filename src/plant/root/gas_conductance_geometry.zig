const std = @import("std");

pub const Inputs = struct {
    total_root_and_mycorrhizal_carbon_g_c: f64,
    negligible_root_carbon_g_c: f64,
    primary_root_layer_fraction: f64,
    significant_primary_root_fraction: f64,
    plant_population: f64,
    primary_axis_count: f64,
    primary_root_radius_m: f64,
    primary_root_depth_m: f64,
    secondary_axis_count: f64,
    secondary_root_radius_m: f64,
    average_secondary_root_length_m: f64,
    circle_pi_compatibility: f64,
};

pub const Result = struct {
    geometry_active: bool,
    primary_cross_section_per_length_m: f64,
    secondary_cross_section_per_length_m: f64,
    effective_cross_section_per_length_m: f64,
};

/// UPTAKE.F 2015--2029. Computes root conductance geometry after the source's
/// strict root-mass and primary-layer-fraction admission gates.
pub fn calculate(inputs: Inputs) !Result {
    try validateCommon(inputs);
    if (!(inputs.total_root_and_mycorrhizal_carbon_g_c >
        inputs.negligible_root_carbon_g_c and
        inputs.primary_root_layer_fraction >
            inputs.significant_primary_root_fraction))
        return .{
            .geometry_active = false,
            .primary_cross_section_per_length_m = 0,
            .secondary_cross_section_per_length_m = 0,
            .effective_cross_section_per_length_m = 0,
        };
    if (inputs.primary_root_depth_m <= 0 or
        inputs.average_secondary_root_length_m <= 0 or
        inputs.primary_root_layer_fraction <= 0)
        return error.InvalidRootGasConductanceGeometryDivisor;
    const primary =
        @max(inputs.plant_population, inputs.primary_axis_count) *
        inputs.circle_pi_compatibility *
        inputs.primary_root_radius_m * inputs.primary_root_radius_m /
        inputs.primary_root_depth_m;
    const secondary =
        (inputs.secondary_axis_count *
            inputs.circle_pi_compatibility *
            inputs.secondary_root_radius_m *
            inputs.secondary_root_radius_m /
            inputs.average_secondary_root_length_m) /
        inputs.primary_root_layer_fraction;
    const effective =
        if (secondary > primary)
            primary * secondary / (primary + secondary)
        else
            primary;
    inline for (.{ primary, secondary, effective }) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteRootGasConductanceGeometryResult;
    return .{
        .geometry_active = true,
        .primary_cross_section_per_length_m = primary,
        .secondary_cross_section_per_length_m = secondary,
        .effective_cross_section_per_length_m = effective,
    };
}

fn validateCommon(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidRootGasConductanceGeometryInput;
    if (inputs.total_root_and_mycorrhizal_carbon_g_c < 0 or
        inputs.negligible_root_carbon_g_c < 0 or
        inputs.primary_root_layer_fraction < 0 or
        inputs.significant_primary_root_fraction < 0 or
        inputs.plant_population < 0 or
        inputs.primary_axis_count < 0 or
        inputs.primary_root_radius_m < 0 or
        inputs.secondary_axis_count < 0 or
        inputs.secondary_root_radius_m < 0 or
        inputs.circle_pi_compatibility <= 0)
        return error.InvalidRootGasConductanceGeometryInput;
}

fn sourceInputs() Inputs {
    return .{
        .total_root_and_mycorrhizal_carbon_g_c = 10,
        .negligible_root_carbon_g_c = 1e-12,
        .primary_root_layer_fraction = 0.5,
        .significant_primary_root_fraction = 1e-12,
        .plant_population = 2,
        .primary_axis_count = 3,
        .primary_root_radius_m = 0.01,
        .primary_root_depth_m = 1,
        .secondary_axis_count = 100,
        .secondary_root_radius_m = 0.005,
        .average_secondary_root_length_m = 0.2,
        .circle_pi_compatibility = 3.1416,
    };
}

test "UPTAKE secondary-larger branch uses source harmonic cross section" {
    const inputs = sourceInputs();
    const result = try calculate(inputs);
    const primary = 3.0 * 3.1416 * 0.01 * 0.01 / 1.0;
    const secondary = (100.0 * 3.1416 * 0.005 * 0.005 / 0.2) / 0.5;
    try std.testing.expect(result.geometry_active);
    try std.testing.expectEqual(primary, result.primary_cross_section_per_length_m);
    try std.testing.expectEqual(secondary, result.secondary_cross_section_per_length_m);
    try std.testing.expectApproxEqAbs(
        primary * secondary / (primary + secondary),
        result.effective_cross_section_per_length_m,
        2e-19,
    );
}

test "secondary-not-larger branch retains primary cross section" {
    var inputs = sourceInputs();
    inputs.secondary_axis_count = 0;
    const result = try calculate(inputs);
    try std.testing.expectEqual(
        result.primary_cross_section_per_length_m,
        result.effective_cross_section_per_length_m,
    );
}

test "inactive strict gate returns zero without evaluating geometry divisors" {
    var inputs = sourceInputs();
    inputs.total_root_and_mycorrhizal_carbon_g_c =
        inputs.negligible_root_carbon_g_c;
    inputs.primary_root_depth_m = 0;
    inputs.average_secondary_root_length_m = 0;
    const result = try calculate(inputs);
    try std.testing.expect(!result.geometry_active);
    try std.testing.expectEqual(@as(f64, 0), result.effective_cross_section_per_length_m);
}
