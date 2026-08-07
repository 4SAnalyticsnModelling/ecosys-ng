const std = @import("std");
const metabolism = @import("plant_root_metabolism.zig");

pub const AxisInputs = struct {
    emerged: bool,
    mobile_carbon_concentration_g_c_per_g_c: f64,
    mobile_nitrogen_concentration_g_n_per_g_c: f64,
    mobile_phosphorus_concentration_g_p_per_g_c: f64,
    senescence: metabolism.SecondaryRootSenescenceInputs,
    woody_carbon_fraction: [2]f64,
    woody_nitrogen_fraction: [2]f64,
    woody_phosphorus_fraction: [2]f64,
    kinetics: metabolism.RootLitterFractions,
};

pub const State = struct {
    recycling_by_domain_layer_axis: []metabolism.RecyclingFractions,
    senescence_by_domain_layer_axis: []metabolism.SecondaryRootSenescence,
    litter_by_domain_layer_axis: []metabolism.RootLitter,
    litter_by_domain_layer: []metabolism.RootLitter,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_index: usize,
    /// Source `NI(NZ,NY,NX)`: one plant-wide inclusive rooted bound.
    deepest_rooted_layer_index: usize,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    tip_active_by_domain_layer_axis: []const bool,
    axis_inputs_by_domain_layer_axis: []const AxisInputs,
    parameters: metabolism.SecondaryRootParameters,
};

const Calculation = struct {
    recycling: metabolism.RecyclingFractions,
    senescence: metabolism.SecondaryRootSenescence,
    litter: metabolism.RootLitter,
};

fn calculate(parameters: metabolism.SecondaryRootParameters, input: AxisInputs) !Calculation {
    const recycling = try metabolism.secondaryRootRecyclingFractions(input.emerged, input.mobile_carbon_concentration_g_c_per_g_c, input.mobile_nitrogen_concentration_g_n_per_g_c, input.mobile_phosphorus_concentration_g_p_per_g_c, parameters);
    const senescence = try metabolism.primaryRootSenescence(input.senescence, recycling);
    return .{ .recycling = recycling, .senescence = senescence, .litter = try metabolism.secondaryRootLitter(senescence, input.senescence.root_carbon_g_c, input.senescence.root_nitrogen_g_n, input.senescence.root_phosphorus_g_p, input.woody_carbon_fraction, input.woody_nitrogen_fraction, input.woody_phosphorus_fraction, input.kinetics) };
}

fn addLitter(total: *metabolism.RootLitter, addition: metabolism.RootLitter) !void {
    inline for (@typeInfo(metabolism.RootLitter).@"struct".fields) |field| for (&@field(total, field.name), @field(addition, field.name)) |*target, value| {
        if (!std.math.isFinite(target.*) or target.* < 0) return error.InvalidPrimaryRootLitterAccumulator;
        target.* += value;
        if (!std.math.isFinite(target.*) or target.* < 0) return error.NonFinitePrimaryRootLitterAccumulator;
    };
}

/// Exact grosub.f lines 6645--6747 primary recycling, non-phenological
/// senescence, and four-kinetic-pool litter sweep. Masses are g C/N/P per
/// biological hour; fractions are dimensionless.
pub fn apply(state: State, inputs: Inputs) !void {
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.PrimaryRootSenescenceSweepDimensionOverflow;
    const values = std.math.mul(usize, domain_layers, inputs.root_axis_count) catch return error.PrimaryRootSenescenceSweepDimensionOverflow;
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or
        inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.tip_active_by_domain_layer_axis.len != values or inputs.axis_inputs_by_domain_layer_axis.len != values or
        state.recycling_by_domain_layer_axis.len != values or state.senescence_by_domain_layer_axis.len != values or state.litter_by_domain_layer_axis.len != values or state.litter_by_domain_layer.len != domain_layers)
        return error.PrimaryRootSenescenceSweepDimensionMismatch;
    if (!std.math.isFinite(inputs.minimum_active_layer_thickness_m) or inputs.minimum_active_layer_thickness_m < 0) return error.InvalidPrimaryRootSenescenceSweepThreshold;
    for (inputs.layer_thickness_m) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootSenescenceSweepLayerThickness;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count) return error.InvalidPrimaryRootSenescenceSweepLayerRange;
    for (1..inputs.biological_domain_count) |domain| for (0..inputs.soil_layer_count) |layer| for (0..inputs.root_axis_count) |axis| {
        if (inputs.tip_active_by_domain_layer_axis[(domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis]) return error.PrimaryRootTipOutsideHostDomain;
    };

    for (0..inputs.biological_domain_count) |domain| for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
        if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
        const domain_layer = domain * inputs.soil_layer_count + layer;
        var total = state.litter_by_domain_layer[domain_layer];
        for (0..inputs.root_axis_count) |axis| {
            const index = domain_layer * inputs.root_axis_count + axis;
            if (!inputs.tip_active_by_domain_layer_axis[index]) continue;
            const result = try calculate(inputs.parameters, inputs.axis_inputs_by_domain_layer_axis[index]);
            try addLitter(&total, result.litter);
        }
    };
    for (0..inputs.biological_domain_count) |domain| for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
        if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
        const domain_layer = domain * inputs.soil_layer_count + layer;
        for (0..inputs.root_axis_count) |axis| {
            const index = domain_layer * inputs.root_axis_count + axis;
            if (!inputs.tip_active_by_domain_layer_axis[index]) continue;
            const result = try calculate(inputs.parameters, inputs.axis_inputs_by_domain_layer_axis[index]);
            state.recycling_by_domain_layer_axis[index] = result.recycling;
            state.senescence_by_domain_layer_axis[index] = result.senescence;
            state.litter_by_domain_layer_axis[index] = result.litter;
            try addLitter(&state.litter_by_domain_layer[domain_layer], result.litter);
        }
    };
}

fn zeroLitter() metabolism.RootLitter {
    return std.mem.zeroes(metabolism.RootLitter);
}

fn axisInput(carbon_g_c: f64) AxisInputs {
    const kinetics = [4]f64{ 0.1, 0.2, 0.3, 0.4 };
    return .{
        .emerged = true,
        .mobile_carbon_concentration_g_c_per_g_c = 0.1,
        .mobile_nitrogen_concentration_g_n_per_g_c = 0.02,
        .mobile_phosphorus_concentration_g_p_per_g_c = 0.002,
        .senescence = .{ .oxygen_unlimited_substrate_minus_maintenance_g_c_per_h = -10, .actual_substrate_minus_maintenance_g_c_per_h = -10, .root_carbon_g_c = carbon_g_c, .root_nitrogen_g_n = carbon_g_c * 0.1, .root_phosphorus_g_p = carbon_g_c * 0.01, .oxygen_limitation = 1, .phenological_remobilization_enabled = true, .root_remobilization_enabled = true, .storage_exchange_fraction_per_h = 1, .remobilization_elapsed_h = 100, .full_senescence_h = 480, .biological_timestep_h = 1, .structural_presence_threshold_g_c = 1e-12 },
        .woody_carbon_fraction = .{ 0.25, 0.75 },
        .woody_nitrogen_fraction = .{ 0.25, 0.75 },
        .woody_phosphorus_fraction = .{ 0.25, 0.75 },
        .kinetics = .{ .woody_carbon = kinetics, .woody_nitrogen = kinetics, .woody_phosphorus = kinetics, .nonwoody_carbon = kinetics, .nonwoody_nitrogen = kinetics, .nonwoody_phosphorus = kinetics },
    };
}

test "GROSUB primary senescence suppresses phenology and accumulates tip litter" {
    var axis_inputs = [_]AxisInputs{axisInput(4)} ** 12;
    axis_inputs[4] = axisInput(12);
    var recycling = [_]metabolism.RecyclingFractions{std.mem.zeroes(metabolism.RecyclingFractions)} ** 12;
    var losses = [_]metabolism.SecondaryRootSenescence{std.mem.zeroes(metabolism.SecondaryRootSenescence)} ** 12;
    var axis_litter = [_]metabolism.RootLitter{zeroLitter()} ** 12;
    var layer_litter = [_]metabolism.RootLitter{zeroLitter()} ** 4;
    try apply(.{ .recycling_by_domain_layer_axis = &recycling, .senescence_by_domain_layer_axis = &losses, .litter_by_domain_layer_axis = &axis_litter, .litter_by_domain_layer = &layer_litter }, .{
        .biological_domain_count = 2,
        .soil_layer_count = 2,
        .root_axis_count = 3,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 1,
        .layer_thickness_m = &.{ 0.2, 0.2 },
        .minimum_active_layer_thickness_m = 0.001,
        .tip_active_by_domain_layer_axis = &.{ true, false, false, false, true, false, false, false, false, false, false, false },
        .axis_inputs_by_domain_layer_axis = &axis_inputs,
        .parameters = metabolism.compatibilitySecondaryRootParameters(),
    });
    try std.testing.expectEqual(@as(f64, 0), losses[0].phenological_senescence_g_c_per_h);
    try std.testing.expectEqual(@as(f64, 0), losses[4].phenological_senescence_g_c_per_h);
    try std.testing.expectApproxEqAbs(axis_litter[0].woody_carbon_g_c[0], layer_litter[0].woody_carbon_g_c[0], 1e-15);
    try std.testing.expectEqual(@as(f64, 0), layer_litter[2].woody_carbon_g_c[0]);
}

test "GROSUB primary litter and recyclable return partition senesced C N P" {
    const input = axisInput(8);
    const result = try calculate(metabolism.compatibilitySecondaryRootParameters(), input);
    inline for (.{
        .{ input.senescence.root_carbon_g_c, result.senescence.recyclable_carbon_g_c, result.litter.woody_carbon_g_c, result.litter.nonwoody_carbon_g_c, input.woody_carbon_fraction[1] },
        .{ input.senescence.root_nitrogen_g_n, result.senescence.recyclable_nitrogen_g_n, result.litter.woody_nitrogen_g_n, result.litter.nonwoody_nitrogen_g_n, input.woody_nitrogen_fraction[1] },
        .{ input.senescence.root_phosphorus_g_p, result.senescence.recyclable_phosphorus_g_p, result.litter.woody_phosphorus_g_p, result.litter.nonwoody_phosphorus_g_p, input.woody_phosphorus_fraction[1] },
    }) |balance| {
        var litter: f64 = 0;
        for (balance[2]) |value| litter += value;
        for (balance[3]) |value| litter += value;
        try std.testing.expectApproxEqAbs(result.senescence.senesced_fraction * balance[0], litter + result.senescence.senesced_fraction * balance[1] * balance[4], 1e-12);
    }
}

test "GROSUB primary senescence litter rolls back on invalid late tip" {
    var axis_inputs = [_]AxisInputs{axisInput(4)} ** 6;
    axis_inputs[5].senescence.root_carbon_g_c = std.math.nan(f64);
    var recycling = [_]metabolism.RecyclingFractions{.{ .carbon = 0.1, .nitrogen = 0.2, .phosphorus = 0.3 }} ** 6;
    var losses = [_]metabolism.SecondaryRootSenescence{std.mem.zeroes(metabolism.SecondaryRootSenescence)} ** 6;
    var axis_litter = [_]metabolism.RootLitter{zeroLitter()} ** 6;
    var layer_litter = [_]metabolism.RootLitter{zeroLitter()};
    const before = recycling;
    try std.testing.expectError(error.NonFiniteSecondaryRootSenescenceInput, apply(.{ .recycling_by_domain_layer_axis = &recycling, .senescence_by_domain_layer_axis = &losses, .litter_by_domain_layer_axis = &axis_litter, .litter_by_domain_layer = &layer_litter }, .{
        .biological_domain_count = 1,
        .soil_layer_count = 1,
        .root_axis_count = 6,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 0,
        .layer_thickness_m = &.{0.2},
        .minimum_active_layer_thickness_m = 0.001,
        .tip_active_by_domain_layer_axis = &.{ true, true, true, true, true, true },
        .axis_inputs_by_domain_layer_axis = &axis_inputs,
        .parameters = metabolism.compatibilitySecondaryRootParameters(),
    }));
    try std.testing.expectEqualDeep(before, recycling);
}
