const std = @import("std");

pub const Parameters = struct {
    seed_volume_per_carbon_m3_g_c: f64,
    seed_length_multiplier: f64,
    seed_shape_coefficient: f64,
    pi_approximation: f64,
    seed_length_exponent: f64,
    seed_area_multiplier: f64,
};

pub const Inputs = struct {
    seed_carbon_mass_g_c: f64,
    initial_seeding_depth_m: f64,
    layer_boundary_depths_m: []const f64,
    first_active_layer_index: usize,
    root_nitrogen_to_carbon_g_n_g_c: f64,
    root_phosphorus_to_carbon_g_p_g_c: f64,
    root_dry_matter_growth_yield: f64,
};

pub const Result = struct {
    seed_volume_m3: f64,
    seed_length_m: f64,
    seed_surface_area_m2: f64,
    seeding_depth_m: f64,
    seeding_layer_index: usize,
    upper_rooting_layer_index: usize,
    root_nitrogen_growth_yield_g_n_g_c: f64,
    root_phosphorus_growth_yield_g_p_g_c: f64,
};

pub const InitializationError = error{
    NonFiniteInput,
    InvalidParameter,
    InvalidLayerBoundaries,
    EmptyRootOrderExtent,
    SeedingDepthOutsideActiveProfile,
    NonFiniteResult,
};

/// Translates `startq.f` lines 360--375 for one plant species.
///
/// Layer indexes are zero-based Zig indexes. `root_order_layer_indexes` has a
/// runtime extent and replaces the legacy fixed ten root orders.
pub fn initialize(
    inputs: Inputs,
    parameters: Parameters,
    root_order_layer_indexes: []usize,
) InitializationError!Result {
    try validate(inputs, parameters);
    if (root_order_layer_indexes.len == 0) return error.EmptyRootOrderExtent;

    const seed_volume_m3 =
        inputs.seed_carbon_mass_g_c * parameters.seed_volume_per_carbon_m3_g_c;
    const seed_length_m = parameters.seed_length_multiplier * std.math.pow(
        f64,
        parameters.seed_shape_coefficient * seed_volume_m3 / parameters.pi_approximation,
        parameters.seed_length_exponent,
    );
    const seed_surface_area_m2 = parameters.seed_area_multiplier *
        parameters.pi_approximation *
        std.math.pow(f64, seed_length_m / 2.0, 2.0);

    var seeding_layer_index: ?usize = null;
    const layer_count = inputs.layer_boundary_depths_m.len - 1;
    for (inputs.first_active_layer_index..layer_count) |layer_index| {
        if (inputs.initial_seeding_depth_m >= inputs.layer_boundary_depths_m[layer_index] and
            inputs.initial_seeding_depth_m < inputs.layer_boundary_depths_m[layer_index + 1])
        {
            seeding_layer_index = layer_index;
        }
    }
    const selected_layer = seeding_layer_index orelse
        return error.SeedingDepthOutsideActiveProfile;

    const nitrogen_growth_yield = inputs.root_nitrogen_to_carbon_g_n_g_c *
        inputs.root_dry_matter_growth_yield;
    const phosphorus_growth_yield = inputs.root_phosphorus_to_carbon_g_p_g_c *
        inputs.root_dry_matter_growth_yield;
    inline for (.{
        seed_volume_m3,
        seed_length_m,
        seed_surface_area_m2,
        nitrogen_growth_yield,
        phosphorus_growth_yield,
    }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteResult;
    }

    @memset(root_order_layer_indexes, selected_layer);
    return .{
        .seed_volume_m3 = seed_volume_m3,
        .seed_length_m = seed_length_m,
        .seed_surface_area_m2 = seed_surface_area_m2,
        .seeding_depth_m = inputs.initial_seeding_depth_m,
        .seeding_layer_index = selected_layer,
        .upper_rooting_layer_index = selected_layer,
        .root_nitrogen_growth_yield_g_n_g_c = nitrogen_growth_yield,
        .root_phosphorus_growth_yield_g_p_g_c = phosphorus_growth_yield,
    };
}

fn validate(inputs: Inputs, parameters: Parameters) InitializationError!void {
    inline for (std.meta.fields(Inputs)) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) {
            return error.NonFiniteInput;
        }
    }
    inline for (std.meta.fields(Parameters)) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (value <= 0.0) return error.InvalidParameter;
    }
    if (inputs.seed_carbon_mass_g_c < 0.0 or inputs.initial_seeding_depth_m < 0.0 or
        inputs.root_nitrogen_to_carbon_g_n_g_c < 0.0 or
        inputs.root_phosphorus_to_carbon_g_p_g_c < 0.0 or
        inputs.root_dry_matter_growth_yield < 0.0)
    {
        return error.InvalidParameter;
    }
    if (inputs.layer_boundary_depths_m.len < 2 or
        inputs.first_active_layer_index >= inputs.layer_boundary_depths_m.len - 1)
    {
        return error.InvalidLayerBoundaries;
    }
    var previous_depth_m = inputs.layer_boundary_depths_m[0];
    if (!std.math.isFinite(previous_depth_m) or previous_depth_m < 0.0) {
        return error.InvalidLayerBoundaries;
    }
    for (inputs.layer_boundary_depths_m[1..]) |depth_m| {
        if (!std.math.isFinite(depth_m) or depth_m <= previous_depth_m) {
            return error.InvalidLayerBoundaries;
        }
        previous_depth_m = depth_m;
    }
}

fn legacyParameters() Parameters {
    return .{
        .seed_volume_per_carbon_m3_g_c = 5.0e-6,
        .seed_length_multiplier = 2.0,
        .seed_shape_coefficient = 0.75,
        .pi_approximation = 3.1416,
        .seed_length_exponent = 0.33,
        .seed_area_multiplier = 4.0,
    };
}

test "seed geometry and runtime root orders preserve STARTQ equations" {
    var root_order_layers: [3]usize = undefined;
    const result = try initialize(.{
        .seed_carbon_mass_g_c = 0.2,
        .initial_seeding_depth_m = 0.12,
        .layer_boundary_depths_m = &.{ 0.0, 0.1, 0.2, 0.4 },
        .first_active_layer_index = 0,
        .root_nitrogen_to_carbon_g_n_g_c = 0.04,
        .root_phosphorus_to_carbon_g_p_g_c = 0.005,
        .root_dry_matter_growth_yield = 0.8,
    }, legacyParameters(), &root_order_layers);

    const expected_volume_m3 = 0.2 * 5.0e-6;
    const expected_length_m =
        2.0 * std.math.pow(f64, 0.75 * expected_volume_m3 / 3.1416, 0.33);
    try std.testing.expectApproxEqRel(expected_volume_m3, result.seed_volume_m3, 1.0e-14);
    try std.testing.expectApproxEqRel(expected_length_m, result.seed_length_m, 1.0e-14);
    try std.testing.expectEqual(@as(usize, 1), result.seeding_layer_index);
    try std.testing.expectEqualSlices(usize, &.{ 1, 1, 1 }, &root_order_layers);
    try std.testing.expectEqual(
        @as(f64, 0.032),
        result.root_nitrogen_growth_yield_g_n_g_c,
    );
}

test "seeding depth outside profile fails before root-order mutation" {
    var root_order_layers = [_]usize{ 8, 9 };
    try std.testing.expectError(error.SeedingDepthOutsideActiveProfile, initialize(.{
        .seed_carbon_mass_g_c = 0.2,
        .initial_seeding_depth_m = 0.5,
        .layer_boundary_depths_m = &.{ 0.0, 0.1, 0.2 },
        .first_active_layer_index = 0,
        .root_nitrogen_to_carbon_g_n_g_c = 0.04,
        .root_phosphorus_to_carbon_g_p_g_c = 0.005,
        .root_dry_matter_growth_yield = 0.8,
    }, legacyParameters(), &root_order_layers));
    try std.testing.expectEqualSlices(usize, &.{ 8, 9 }, &root_order_layers);
}
