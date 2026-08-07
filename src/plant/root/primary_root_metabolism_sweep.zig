const std = @import("std");
const metabolism = @import("plant_root_metabolism.zig");

pub const State = struct {
    /// Flattened `[domain][layer][axis]`; inactive entries remain untouched.
    result_by_domain_layer_axis: []metabolism.SecondaryRootResult,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_index: usize,
    /// Source `NI(NZ,NY,NX)`: shared by all biological domains of the plant.
    deepest_rooted_layer_index: usize,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    tip_active_by_domain_layer_axis: []const bool,
    calculation_by_domain_layer_axis: []const metabolism.PrimaryRootInputs,
    parameters: metabolism.SecondaryRootParameters,
};

/// grosub.f lines 6498--6632 primary-root CNPG, respiration, growth, and N/P
/// demand staging in exact N,L,NR order. C/N/P rates are g C/N/P per
/// biological hour; response and sink operands are dimensionless.
pub fn apply(state: State, inputs: Inputs) !void {
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch
        return error.PrimaryRootMetabolismSweepDimensionOverflow;
    const values = std.math.mul(usize, domain_layers, inputs.root_axis_count) catch
        return error.PrimaryRootMetabolismSweepDimensionOverflow;
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or
        inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.tip_active_by_domain_layer_axis.len != values or
        inputs.calculation_by_domain_layer_axis.len != values or
        state.result_by_domain_layer_axis.len != values)
        return error.PrimaryRootMetabolismSweepDimensionMismatch;
    if (!std.math.isFinite(inputs.minimum_active_layer_thickness_m) or inputs.minimum_active_layer_thickness_m < 0)
        return error.InvalidPrimaryRootMetabolismSweepThreshold;
    for (inputs.layer_thickness_m) |thickness_m|
        if (!std.math.isFinite(thickness_m) or thickness_m < 0)
            return error.InvalidPrimaryRootMetabolismSweepLayerThickness;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count)
        return error.InvalidPrimaryRootMetabolismSweepLayerRange;
    for (1..inputs.biological_domain_count) |domain| for (0..inputs.soil_layer_count) |layer| for (0..inputs.root_axis_count) |axis| {
        const index = (domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
        if (inputs.tip_active_by_domain_layer_axis[index]) return error.PrimaryRootTipOutsideHostDomain;
    };

    // Preflight preserves exact source N outer, L middle, NR inner order.
    for (0..inputs.biological_domain_count) |domain| {
        const deepest = inputs.deepest_rooted_layer_index;
        for (inputs.planting_layer_index..deepest + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            for (0..inputs.root_axis_count) |axis| {
                const index = (domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
                if (!inputs.tip_active_by_domain_layer_axis[index]) continue;
                _ = try metabolism.primaryRootMetabolism(inputs.parameters, inputs.calculation_by_domain_layer_axis[index]);
            }
        }
    }
    for (0..inputs.biological_domain_count) |domain| {
        const deepest = inputs.deepest_rooted_layer_index;
        for (inputs.planting_layer_index..deepest + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            for (0..inputs.root_axis_count) |axis| {
                const index = (domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
                if (!inputs.tip_active_by_domain_layer_axis[index]) continue;
                state.result_by_domain_layer_axis[index] = try metabolism.primaryRootMetabolism(
                    inputs.parameters,
                    inputs.calculation_by_domain_layer_axis[index],
                );
            }
        }
    }
}

fn calculation(fraction: f64, bottom_tip: bool) metabolism.PrimaryRootInputs {
    return .{
        .shared = .{
            .mobile_carbon_g_c = 10,
            .nonstructural_nitrogen_g_n = 1,
            .nonstructural_phosphorus_g_p = 0.1,
            .root_carbon_g_c = 20,
            .root_nitrogen_g_n = 0.4,
            .root_nitrogen_to_carbon_ratio_g_n_per_g_c = 0.02,
            .root_phosphorus_to_carbon_ratio_g_p_per_g_c = 0.002,
            .root_growth_yield_g_c_per_g_c = 0.75,
            .active_root_fraction = fraction,
            .biological_timestep_h = 1,
            .substrate_temperature_response = 1,
            .maintenance_temperature_response = 1,
            .acidity_response = 1,
            .substrate_feedback = 1,
            .oxygen_limitation = 0.8,
            .substrate_water_response = 1,
            .maintenance_water_response = 1,
        },
        .primary_tip_at_or_below_profile_bottom = bottom_tip,
    };
}

test "GROSUB primary metabolism stages runtime tips beyond five axes" {
    const domains = 2;
    const layers = 2;
    const axes = 7;
    const count = domains * layers * axes;
    var calculations: [count]metabolism.PrimaryRootInputs = undefined;
    for (&calculations, 0..) |*item, index| item.* = calculation(@as(f64, @floatFromInt(index + 1)) / 100, false);
    calculations[13].primary_tip_at_or_below_profile_bottom = true;
    var active = [_]bool{false} ** count;
    active[0] = true;
    active[13] = true;
    var results = [_]metabolism.SecondaryRootResult{std.mem.zeroes(metabolism.SecondaryRootResult)} ** count;
    results[7].nutrient_feedback = -1;
    try apply(.{ .result_by_domain_layer_axis = &results }, .{
        .biological_domain_count = domains,
        .soil_layer_count = layers,
        .root_axis_count = axes,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 1,
        .layer_thickness_m = &.{ 0.2, 0.2 },
        .minimum_active_layer_thickness_m = 0.001,
        .tip_active_by_domain_layer_axis = &active,
        .calculation_by_domain_layer_axis = &calculations,
        .parameters = metabolism.compatibilitySecondaryRootParameters(),
    });
    try std.testing.expect(results[0].nutrient_feedback > 0);
    try std.testing.expectEqual(@as(f64, -1), results[7].nutrient_feedback);
    try std.testing.expect(results[13].substrate_respiration_oxygen_unlimited_g_c_per_h <= results[13].maintenance_respiration_g_c_per_h);
}

test "GROSUB primary metabolism closes carbon growth and CNP demand identities" {
    var result = [_]metabolism.SecondaryRootResult{std.mem.zeroes(metabolism.SecondaryRootResult)};
    const input = calculation(0.4, false);
    try apply(.{ .result_by_domain_layer_axis = &result }, .{
        .biological_domain_count = 1,
        .soil_layer_count = 1,
        .root_axis_count = 1,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 0,
        .layer_thickness_m = &.{0.2},
        .minimum_active_layer_thickness_m = 0.001,
        .tip_active_by_domain_layer_axis = &.{true},
        .calculation_by_domain_layer_axis = &.{input},
        .parameters = metabolism.compatibilitySecondaryRootParameters(),
    });
    const value = result[0];
    try std.testing.expectApproxEqAbs(value.growth_and_respiration_carbon_actual_g_c_per_h, value.root_growth_actual_g_c_per_h + value.growth_respiration_actual_g_c_per_h, 1e-15);
    try std.testing.expect(value.nitrogen_growth_actual_g_n_per_h <= input.shared.active_root_fraction * input.shared.nonstructural_nitrogen_g_n);
    try std.testing.expect(value.phosphorus_growth_actual_g_p_per_h <= input.shared.active_root_fraction * input.shared.nonstructural_phosphorus_g_p);
}

test "GROSUB primary metabolism sweep rolls back on invalid late tip" {
    var calculations = [_]metabolism.PrimaryRootInputs{calculation(0.2, false)} ** 6;
    calculations[5].shared.mobile_carbon_g_c = std.math.nan(f64);
    var results = [_]metabolism.SecondaryRootResult{std.mem.zeroes(metabolism.SecondaryRootResult)} ** 6;
    results[0].nutrient_feedback = 0.25;
    const before = results;
    try std.testing.expectError(error.InvalidSecondaryRootMetabolismInput, apply(.{ .result_by_domain_layer_axis = &results }, .{
        .biological_domain_count = 1,
        .soil_layer_count = 1,
        .root_axis_count = 6,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 0,
        .layer_thickness_m = &.{0.2},
        .minimum_active_layer_thickness_m = 0.001,
        .tip_active_by_domain_layer_axis = &.{ true, true, true, true, true, true },
        .calculation_by_domain_layer_axis = &calculations,
        .parameters = metabolism.compatibilitySecondaryRootParameters(),
    }));
    try std.testing.expectEqualDeep(before, results);
}
