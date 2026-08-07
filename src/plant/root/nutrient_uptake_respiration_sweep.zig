const std = @import("std");

pub const NutrientUptake = struct {
    /// NH4 non-band, NH4 band, NO3 non-band, NO3 band, H2PO4 non-band,
    /// H2PO4 band, HPO4 non-band, HPO4 band, in source summation order.
    actual_g_element: [8]f64,
    oxygen_unlimited_g_element: [8]f64,
    carbon_unlimited_g_element: [8]f64,
};

pub const RootLayerState = struct {
    mobile_carbon_g_c: []f64,
    actual_respiration_carbon_exchange_g_c: []f64,
    oxygen_unlimited_respiration_g_c: []f64,
    carbon_unlimited_respiration_g_c: []f64,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    nutrient_uptake_by_domain_layer: []const NutrientUptake,
    nutrient_uptake_respiration_g_c_per_g_element: f64,
};

/// Exact GROSUB root-loop publication from NIX=NG and lines 5742--5794:
/// domains N are outermost, rooted layers L are inner, and only layers with
/// `DLYR(3,L) > DLYRM` participate. RCO2A retains the source negative carbon
/// exchange sign while RCO2M/RCO2N accumulate positive potential respiration.
pub fn apply(state: RootLayerState, inputs: Inputs) !void {
    const item_count = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch
        return error.RootNutrientRespirationDimensionOverflow;
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or
        inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.nutrient_uptake_by_domain_layer.len != item_count or
        state.mobile_carbon_g_c.len != item_count or
        state.actual_respiration_carbon_exchange_g_c.len != item_count or
        state.oxygen_unlimited_respiration_g_c.len != item_count or
        state.carbon_unlimited_respiration_g_c.len != item_count)
        return error.RootNutrientRespirationDimensionMismatch;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or
        inputs.deepest_rooted_layer_index >= inputs.soil_layer_count)
        return error.InvalidRootedLayerRange;
    inline for (.{ inputs.minimum_active_layer_thickness_m, inputs.nutrient_uptake_respiration_g_c_per_g_element }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootNutrientRespirationInput;

    for (inputs.layer_thickness_m) |thickness_m|
        if (!std.math.isFinite(thickness_m) or thickness_m < 0)
            return error.InvalidRootLayerThickness;

    // Preflight in exact N-outer/L-inner source order.
    for (0..inputs.biological_domain_count) |domain| {
        for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const index = domain * inputs.soil_layer_count + layer;
            const respiration = try calculate(inputs.nutrient_uptake_by_domain_layer[index], inputs.nutrient_uptake_respiration_g_c_per_g_element);
            inline for (.{
                state.mobile_carbon_g_c[index] - respiration.actual_g_c,
                state.actual_respiration_carbon_exchange_g_c[index] - respiration.actual_g_c,
                state.oxygen_unlimited_respiration_g_c[index] + respiration.oxygen_unlimited_g_c,
                state.carbon_unlimited_respiration_g_c[index] + respiration.carbon_unlimited_g_c,
            }) |value| if (!std.math.isFinite(value))
                return error.NonFiniteRootNutrientRespirationPublication;
            if (state.mobile_carbon_g_c[index] < respiration.actual_g_c)
                return error.RootNutrientRespirationWouldOverdrawCarbon;
        }
    }

    for (0..inputs.biological_domain_count) |domain| {
        for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const index = domain * inputs.soil_layer_count + layer;
            const respiration = try calculate(inputs.nutrient_uptake_by_domain_layer[index], inputs.nutrient_uptake_respiration_g_c_per_g_element);
            state.oxygen_unlimited_respiration_g_c[index] += respiration.oxygen_unlimited_g_c;
            state.carbon_unlimited_respiration_g_c[index] += respiration.carbon_unlimited_g_c;
            state.actual_respiration_carbon_exchange_g_c[index] -= respiration.actual_g_c;
            state.mobile_carbon_g_c[index] -= respiration.actual_g_c;
        }
    }
}

const Respiration = struct {
    actual_g_c: f64,
    oxygen_unlimited_g_c: f64,
    carbon_unlimited_g_c: f64,
};

fn calculate(uptake: NutrientUptake, coefficient: f64) !Respiration {
    var actual: f64 = 0;
    var oxygen_unlimited: f64 = 0;
    var carbon_unlimited: f64 = 0;
    for (0..8) |pool| {
        inline for (.{ uptake.actual_g_element[pool], uptake.oxygen_unlimited_g_element[pool], uptake.carbon_unlimited_g_element[pool] }) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidRootNutrientUptakeFlux;
        actual += uptake.actual_g_element[pool];
        oxygen_unlimited += uptake.oxygen_unlimited_g_element[pool];
        carbon_unlimited += uptake.carbon_unlimited_g_element[pool];
    }
    const result = Respiration{
        .actual_g_c = coefficient * actual,
        .oxygen_unlimited_g_c = coefficient * oxygen_unlimited,
        .carbon_unlimited_g_c = coefficient * carbon_unlimited,
    };
    inline for (.{ result.actual_g_c, result.oxygen_unlimited_g_c, result.carbon_unlimited_g_c }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteRootNutrientRespiration;
    return result;
}

fn testUptake(actual: f64, oxygen: f64, carbon: f64) NutrientUptake {
    return .{
        .actual_g_element = [_]f64{actual} ** 8,
        .oxygen_unlimited_g_element = [_]f64{oxygen} ** 8,
        .carbon_unlimited_g_element = [_]f64{carbon} ** 8,
    };
}

test "GROSUB CUPRL CUPRO CUPRC retain eight-pool source order and signs" {
    var mobile = [_]f64{20};
    var actual = [_]f64{5};
    var oxygen = [_]f64{7};
    var carbon = [_]f64{9};
    try apply(.{ .mobile_carbon_g_c = &mobile, .actual_respiration_carbon_exchange_g_c = &actual, .oxygen_unlimited_respiration_g_c = &oxygen, .carbon_unlimited_respiration_g_c = &carbon }, .{
        .biological_domain_count = 1,
        .soil_layer_count = 1,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 0,
        .layer_thickness_m = &.{1},
        .minimum_active_layer_thickness_m = 0.001,
        .nutrient_uptake_by_domain_layer = &.{testUptake(0.5, 0.75, 1)},
        .nutrient_uptake_respiration_g_c_per_g_element = 0.86,
    });
    const respiration_g_c = 0.86 * 8 * 0.5;
    try std.testing.expectApproxEqAbs(20 - respiration_g_c, mobile[0], 1e-15);
    try std.testing.expectApproxEqAbs(5 - respiration_g_c, actual[0], 1e-15);
    try std.testing.expectApproxEqAbs(7 + 0.86 * 8 * 0.75, oxygen[0], 1e-15);
    try std.testing.expectApproxEqAbs(9 + 0.86 * 8, carbon[0], 2e-15);
    // RCO2A and CPOOLR receive the same negative CUPRL increment.
    try std.testing.expectApproxEqAbs(20 - mobile[0], 5 - actual[0], 2e-15);
}

test "GROSUB layer thickness gate is strict and does not read inactive uptake" {
    var mobile = [_]f64{std.math.nan(f64)};
    var actual = [_]f64{std.math.nan(f64)};
    var oxygen = [_]f64{std.math.nan(f64)};
    var carbon = [_]f64{std.math.nan(f64)};
    try apply(.{ .mobile_carbon_g_c = &mobile, .actual_respiration_carbon_exchange_g_c = &actual, .oxygen_unlimited_respiration_g_c = &oxygen, .carbon_unlimited_respiration_g_c = &carbon }, .{
        .biological_domain_count = 1,
        .soil_layer_count = 1,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 0,
        .layer_thickness_m = &.{0.001},
        .minimum_active_layer_thickness_m = 0.001,
        .nutrient_uptake_by_domain_layer = &.{testUptake(std.math.nan(f64), 0, 0)},
        .nutrient_uptake_respiration_g_c_per_g_element = 0.86,
    });
    try std.testing.expect(std.math.isNan(mobile[0]));
}

test "GROSUB runtime domain layer sweep is atomic on late carbon overdraw" {
    const domains = 5;
    const layers = 11;
    var mobile = [_]f64{100} ** (domains * layers);
    var actual = [_]f64{0} ** (domains * layers);
    var oxygen = [_]f64{0} ** (domains * layers);
    var carbon = [_]f64{0} ** (domains * layers);
    var uptake_flux = [_]NutrientUptake{testUptake(0.1, 0.2, 0.3)} ** (domains * layers);
    const thickness = [_]f64{0.1} ** layers;
    const state = RootLayerState{ .mobile_carbon_g_c = &mobile, .actual_respiration_carbon_exchange_g_c = &actual, .oxygen_unlimited_respiration_g_c = &oxygen, .carbon_unlimited_respiration_g_c = &carbon };
    const inputs = Inputs{ .biological_domain_count = domains, .soil_layer_count = layers, .planting_layer_index = 2, .deepest_rooted_layer_index = 9, .layer_thickness_m = &thickness, .minimum_active_layer_thickness_m = 0.001, .nutrient_uptake_by_domain_layer = &uptake_flux, .nutrient_uptake_respiration_g_c_per_g_element = 0.86 };
    try apply(state, inputs);
    const before_mobile = mobile;
    const before_actual = actual;
    uptake_flux[4 * layers + 9] = testUptake(100, 100, 100);
    try std.testing.expectError(error.RootNutrientRespirationWouldOverdrawCarbon, apply(state, inputs));
    try std.testing.expectEqualDeep(before_mobile, mobile);
    try std.testing.expectEqualDeep(before_actual, actual);
}
