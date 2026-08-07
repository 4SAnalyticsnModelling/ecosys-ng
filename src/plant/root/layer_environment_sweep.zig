const std = @import("std");
const metabolism = @import("plant_root_metabolism.zig");

pub const State = struct {
    /// Flattened `[domain][layer]`; inactive entries remain untouched.
    next_lower_layer_by_domain_layer: []usize,
    environment_by_domain_layer: []metabolism.RootEnvironmentResponses,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    soil_temperature_k_by_layer: []const f64,
    soil_ph_by_layer: []const f64,
    soil_penetration_resistance_mpa_by_layer: []const f64,
    root_total_water_potential_mpa_by_domain_layer: []const f64,
    root_turgor_water_potential_mpa_by_domain_layer: []const f64,
    secondary_root_radius_m_by_domain_layer: []const f64,
    thermal_adaptation_offset_k: f64,
    minimum_extension_water_potential_megapascal: f64,
    shallow_root_profile: bool,
    parameters: metabolism.SecondaryRootParameters,
};

/// grosub.f lines 5975--6034 runtime N,L traversal.
///
/// Temperatures are K, water potentials and penetration resistance are MPa,
/// root radius/layer thickness are m, and responses are dimensionless. The
/// lower-layer identity is a logical zero-based layer, never a storage index.
pub fn apply(state: State, inputs: Inputs) !void {
    const domain_layer_count = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch
        return error.RootEnvironmentSweepDimensionOverflow;
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count < 2 or
        inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.soil_temperature_k_by_layer.len != inputs.soil_layer_count or
        inputs.soil_ph_by_layer.len != inputs.soil_layer_count or
        inputs.soil_penetration_resistance_mpa_by_layer.len != inputs.soil_layer_count or
        inputs.root_total_water_potential_mpa_by_domain_layer.len != domain_layer_count or
        inputs.root_turgor_water_potential_mpa_by_domain_layer.len != domain_layer_count or
        inputs.secondary_root_radius_m_by_domain_layer.len != domain_layer_count or
        state.next_lower_layer_by_domain_layer.len != domain_layer_count or
        state.environment_by_domain_layer.len != domain_layer_count)
        return error.RootEnvironmentSweepDimensionMismatch;
    if (!std.math.isFinite(inputs.minimum_active_layer_thickness_m) or
        inputs.minimum_active_layer_thickness_m < 0)
        return error.InvalidRootLayerThickness;
    for (inputs.layer_thickness_m) |thickness_m|
        if (!std.math.isFinite(thickness_m) or thickness_m < 0)
            return error.InvalidRootLayerThickness;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count - 1)
        return error.InvalidRootEnvironmentLayerRange;

    // Preflight in exact source N outer, L inner order.
    for (0..inputs.biological_domain_count) |domain| {
        const deepest_layer = inputs.deepest_rooted_layer_index;
        for (inputs.planting_layer_index..deepest_layer + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const index = domain * inputs.soil_layer_count + layer;
            _ = try metabolism.nextLowerRootLayer(
                inputs.layer_thickness_m,
                layer,
                inputs.minimum_active_layer_thickness_m,
            );
            _ = try metabolism.rootEnvironmentResponses(
                inputs.parameters,
                inputs.soil_temperature_k_by_layer[layer],
                inputs.thermal_adaptation_offset_k,
                inputs.soil_ph_by_layer[layer],
                inputs.root_total_water_potential_mpa_by_domain_layer[index],
                inputs.root_turgor_water_potential_mpa_by_domain_layer[index],
                inputs.minimum_extension_water_potential_megapascal,
                inputs.soil_penetration_resistance_mpa_by_layer[layer],
                inputs.secondary_root_radius_m_by_domain_layer[index],
                inputs.shallow_root_profile,
            );
        }
    }

    for (0..inputs.biological_domain_count) |domain| {
        const deepest_layer = inputs.deepest_rooted_layer_index;
        for (inputs.planting_layer_index..deepest_layer + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const index = domain * inputs.soil_layer_count + layer;
            state.next_lower_layer_by_domain_layer[index] = try metabolism.nextLowerRootLayer(
                inputs.layer_thickness_m,
                layer,
                inputs.minimum_active_layer_thickness_m,
            );
            state.environment_by_domain_layer[index] = try metabolism.rootEnvironmentResponses(
                inputs.parameters,
                inputs.soil_temperature_k_by_layer[layer],
                inputs.thermal_adaptation_offset_k,
                inputs.soil_ph_by_layer[layer],
                inputs.root_total_water_potential_mpa_by_domain_layer[index],
                inputs.root_turgor_water_potential_mpa_by_domain_layer[index],
                inputs.minimum_extension_water_potential_megapascal,
                inputs.soil_penetration_resistance_mpa_by_layer[layer],
                inputs.secondary_root_radius_m_by_domain_layer[index],
                inputs.shallow_root_profile,
            );
        }
    }
}

fn validInputs(
    total_water: []const f64,
    turgor_water: []const f64,
    radii: []const f64,
) Inputs {
    return .{
        .biological_domain_count = 2,
        .soil_layer_count = 4,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 2,
        .layer_thickness_m = &.{ 0.2, 0.0005, 0.2, 0.2 },
        .minimum_active_layer_thickness_m = 0.001,
        .soil_temperature_k_by_layer = &.{ 290, 291, 292, 293 },
        .soil_ph_by_layer = &.{ 6, 6.5, 7, 7.5 },
        .soil_penetration_resistance_mpa_by_layer = &.{ 0.1, 0.2, 0.3, 0.4 },
        .root_total_water_potential_mpa_by_domain_layer = total_water,
        .root_turgor_water_potential_mpa_by_domain_layer = turgor_water,
        .secondary_root_radius_m_by_domain_layer = radii,
        .thermal_adaptation_offset_k = 0,
        .minimum_extension_water_potential_megapascal = -1.5,
        .shallow_root_profile = false,
        .parameters = metabolism.compatibilitySecondaryRootParameters(),
    };
}

test "GROSUB environment sweep preserves plant-wide rooted bound and thin-layer scan" {
    const total_water = [_]f64{-0.5} ** 8;
    const turgor_water = [_]f64{0.5} ** 8;
    const radii = [_]f64{0.001} ** 8;
    var next_lower = [_]usize{99} ** 8;
    var environments = [_]metabolism.RootEnvironmentResponses{undefined} ** 8;
    try apply(.{ .next_lower_layer_by_domain_layer = &next_lower, .environment_by_domain_layer = &environments }, validInputs(&total_water, &turgor_water, &radii));

    try std.testing.expectEqual(@as(usize, 2), next_lower[0]);
    try std.testing.expectEqual(@as(usize, 99), next_lower[1]);
    try std.testing.expectEqual(@as(usize, 2), next_lower[4]);
    try std.testing.expectEqual(@as(usize, 3), next_lower[6]);
    try std.testing.expect(environments[0].growth_water > 0);
    try std.testing.expect(environments[6].scaled_penetration_resistance_megapascal > environments[4].scaled_penetration_resistance_megapascal);
}

test "GROSUB environment sweep is atomic on invalid late domain-layer input" {
    var total_water = [_]f64{-0.5} ** 8;
    const turgor_water = [_]f64{0.5} ** 8;
    const radii = [_]f64{0.001} ** 8;
    total_water[6] = std.math.nan(f64);
    var next_lower = [_]usize{77} ** 8;
    var environments = [_]metabolism.RootEnvironmentResponses{.{
        .maintenance_temperature = 1,
        .acidity = 2,
        .growth_water = 3,
        .maintenance_water = 4,
        .extension_water = 5,
        .scaled_penetration_resistance_megapascal = 6,
    }} ** 8;
    const before_next = next_lower;
    const before_environment = environments;
    try std.testing.expectError(error.NonFiniteRootEnvironmentInput, apply(
        .{ .next_lower_layer_by_domain_layer = &next_lower, .environment_by_domain_layer = &environments },
        validInputs(&total_water, &turgor_water, &radii),
    ));
    try std.testing.expectEqualDeep(before_next, next_lower);
    try std.testing.expectEqualDeep(before_environment, environments);
}

test "GROSUB environment sweep rejects non-finite admission geometry atomically" {
    const total_water = [_]f64{-0.5} ** 8;
    const turgor_water = [_]f64{0.5} ** 8;
    const radii = [_]f64{0.001} ** 8;
    var inputs = validInputs(&total_water, &turgor_water, &radii);
    inputs.layer_thickness_m = &.{ 0.2, std.math.nan(f64), 0.2, 0.2 };
    var next_lower = [_]usize{55} ** 8;
    var environments = [_]metabolism.RootEnvironmentResponses{.{ .maintenance_temperature = 1, .acidity = 2, .growth_water = 3, .maintenance_water = 4, .extension_water = 5, .scaled_penetration_resistance_megapascal = 6 }} ** 8;
    const before_next = next_lower;
    const before_environment = environments;
    try std.testing.expectError(error.InvalidRootLayerThickness, apply(.{ .next_lower_layer_by_domain_layer = &next_lower, .environment_by_domain_layer = &environments }, inputs));
    try std.testing.expectEqualDeep(before_next, next_lower);
    try std.testing.expectEqualDeep(before_environment, environments);
}
