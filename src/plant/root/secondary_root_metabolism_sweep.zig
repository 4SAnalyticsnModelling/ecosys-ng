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
    deepest_rooted_layer_index: usize,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    active_by_domain_layer_axis: []const bool,
    calculation_by_domain_layer_axis: []const metabolism.SecondaryRootInputs,
    parameters: metabolism.SecondaryRootParameters,
};

/// Exact grosub.f lines 6066--6197 secondary-root calculation traversal.
///
/// C, N, and P pools/rates use g C, g N, and g P (rates per biological
/// hour). All response, feedback, oxygen, and sink-fraction operands are
/// dimensionless. This stages calculations only; source pool publication at
/// lines 6356--6427 is a separate sequential transaction.
pub fn apply(state: State, inputs: Inputs) !void {
    const domain_layer_count = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch
        return error.SecondaryRootMetabolismSweepDimensionOverflow;
    const value_count = std.math.mul(usize, domain_layer_count, inputs.root_axis_count) catch
        return error.SecondaryRootMetabolismSweepDimensionOverflow;
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or
        inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.active_by_domain_layer_axis.len != value_count or
        inputs.calculation_by_domain_layer_axis.len != value_count or
        state.result_by_domain_layer_axis.len != value_count)
        return error.SecondaryRootMetabolismSweepDimensionMismatch;
    if (!std.math.isFinite(inputs.minimum_active_layer_thickness_m) or inputs.minimum_active_layer_thickness_m < 0)
        return error.InvalidSecondaryRootMetabolismSweepThreshold;
    for (inputs.layer_thickness_m) |thickness_m|
        if (!std.math.isFinite(thickness_m) or thickness_m < 0)
            return error.InvalidSecondaryRootMetabolismSweepLayerThickness;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count)
        return error.InvalidSecondaryRootMetabolismSweepLayerRange;

    // Full validation in exact source N outer, L middle, NR inner order.
    for (0..inputs.biological_domain_count) |domain| {
        const deepest_layer = inputs.deepest_rooted_layer_index;
        for (inputs.planting_layer_index..deepest_layer + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const domain_layer = domain * inputs.soil_layer_count + layer;
            for (0..inputs.root_axis_count) |axis| {
                const index = domain_layer * inputs.root_axis_count + axis;
                if (!inputs.active_by_domain_layer_axis[index]) continue;
                _ = try metabolism.secondaryRootMetabolism(inputs.parameters, inputs.calculation_by_domain_layer_axis[index]);
            }
        }
    }

    for (0..inputs.biological_domain_count) |domain| {
        const deepest_layer = inputs.deepest_rooted_layer_index;
        for (inputs.planting_layer_index..deepest_layer + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const domain_layer = domain * inputs.soil_layer_count + layer;
            for (0..inputs.root_axis_count) |axis| {
                const index = domain_layer * inputs.root_axis_count + axis;
                if (!inputs.active_by_domain_layer_axis[index]) continue;
                state.result_by_domain_layer_axis[index] = try metabolism.secondaryRootMetabolism(
                    inputs.parameters,
                    inputs.calculation_by_domain_layer_axis[index],
                );
            }
        }
    }
}

fn calculation(fraction: f64, mobile_c: f64, mobile_n: f64, mobile_p: f64) metabolism.SecondaryRootInputs {
    return .{
        .mobile_carbon_g_c = mobile_c,
        .nonstructural_nitrogen_g_n = mobile_n,
        .nonstructural_phosphorus_g_p = mobile_p,
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
    };
}

test "GROSUB secondary metabolism sweep scales beyond five axes in N L NR order" {
    const domains = 2;
    const layers = 2;
    const axes = 7;
    const count = domains * layers * axes;
    var calculations: [count]metabolism.SecondaryRootInputs = undefined;
    for (&calculations, 0..) |*item, index| item.* = calculation(
        @as(f64, @floatFromInt(index + 1)) / 100.0,
        10,
        1,
        0.1,
    );
    var active = [_]bool{true} ** count;
    active[3] = false;
    var results = [_]metabolism.SecondaryRootResult{std.mem.zeroes(metabolism.SecondaryRootResult)} ** count;
    results[3].nutrient_feedback = -1;
    try apply(.{ .result_by_domain_layer_axis = &results }, .{
        .biological_domain_count = domains,
        .soil_layer_count = layers,
        .root_axis_count = axes,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 1,
        .layer_thickness_m = &.{ 0.2, 0.2 },
        .minimum_active_layer_thickness_m = 0.001,
        .active_by_domain_layer_axis = &active,
        .calculation_by_domain_layer_axis = &calculations,
        .parameters = metabolism.compatibilitySecondaryRootParameters(),
    });

    try std.testing.expectEqual(@as(f64, -1), results[3].nutrient_feedback);
    try std.testing.expect(results[7].nutrient_feedback > 0); // Plant-wide NI includes layer 1 for every domain.
    try std.testing.expect(results[27].growth_respiration_actual_g_c_per_h > results[14].growth_respiration_actual_g_c_per_h);
}

test "GROSUB secondary metabolism closes C growth and CNP demand identities" {
    var result = [_]metabolism.SecondaryRootResult{std.mem.zeroes(metabolism.SecondaryRootResult)};
    const input = calculation(0.4, 12, 2, 0.2);
    try apply(.{ .result_by_domain_layer_axis = &result }, .{
        .biological_domain_count = 1,
        .soil_layer_count = 1,
        .root_axis_count = 1,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 0,
        .layer_thickness_m = &.{0.2},
        .minimum_active_layer_thickness_m = 0.001,
        .active_by_domain_layer_axis = &.{true},
        .calculation_by_domain_layer_axis = &.{input},
        .parameters = metabolism.compatibilitySecondaryRootParameters(),
    });
    const value = result[0];
    try std.testing.expectApproxEqAbs(
        value.growth_and_respiration_carbon_actual_g_c_per_h,
        value.root_growth_actual_g_c_per_h + value.growth_respiration_actual_g_c_per_h,
        1e-15,
    );
    try std.testing.expect(value.nitrogen_growth_actual_g_n_per_h <= input.active_root_fraction * input.nonstructural_nitrogen_g_n);
    try std.testing.expect(value.phosphorus_growth_actual_g_p_per_h <= input.active_root_fraction * input.nonstructural_phosphorus_g_p);
    try std.testing.expectApproxEqAbs(
        value.nitrogen_assimilation_respiration_actual_g_c_per_h,
        metabolism.compatibilitySecondaryRootParameters().nitrogen_assimilation_respiration_g_c_per_g_n * value.nitrogen_growth_actual_g_n_per_h,
        1e-15,
    );
}

test "GROSUB secondary metabolism sweep rolls back on invalid late NR" {
    var calculations = [_]metabolism.SecondaryRootInputs{calculation(0.2, 10, 1, 0.1)} ** 8;
    calculations[7].mobile_carbon_g_c = std.math.nan(f64);
    var results = [_]metabolism.SecondaryRootResult{std.mem.zeroes(metabolism.SecondaryRootResult)} ** 8;
    results[0].nutrient_feedback = 0.25;
    const before = results;
    try std.testing.expectError(error.InvalidSecondaryRootMetabolismInput, apply(
        .{ .result_by_domain_layer_axis = &results },
        .{
            .biological_domain_count = 2,
            .soil_layer_count = 1,
            .root_axis_count = 4,
            .planting_layer_index = 0,
            .deepest_rooted_layer_index = 0,
            .layer_thickness_m = &.{0.2},
            .minimum_active_layer_thickness_m = 0.001,
            .active_by_domain_layer_axis = &.{ true, true, true, true, true, true, true, true },
            .calculation_by_domain_layer_axis = &calculations,
            .parameters = metabolism.compatibilitySecondaryRootParameters(),
        },
    ));
    try std.testing.expectEqualDeep(before, results);
}
