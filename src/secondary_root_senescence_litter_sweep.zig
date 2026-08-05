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
    /// Existing source CSNC/ZSNC/PSNC totals, flattened `[domain][layer]`.
    litter_by_domain_layer: []metabolism.RootLitter,
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
    axis_inputs_by_domain_layer_axis: []const AxisInputs,
    parameters: metabolism.SecondaryRootParameters,
};

const Calculation = struct {
    recycling: metabolism.RecyclingFractions,
    senescence: metabolism.SecondaryRootSenescence,
    litter: metabolism.RootLitter,
};

fn calculate(parameters: metabolism.SecondaryRootParameters, input: AxisInputs) !Calculation {
    const recycling = try metabolism.secondaryRootRecyclingFractions(
        input.emerged,
        input.mobile_carbon_concentration_g_c_per_g_c,
        input.mobile_nitrogen_concentration_g_n_per_g_c,
        input.mobile_phosphorus_concentration_g_p_per_g_c,
        parameters,
    );
    const senescence = try metabolism.secondaryRootSenescence(input.senescence, recycling);
    return .{
        .recycling = recycling,
        .senescence = senescence,
        .litter = try metabolism.secondaryRootLitter(
            senescence,
            input.senescence.root_carbon_g_c,
            input.senescence.root_nitrogen_g_n,
            input.senescence.root_phosphorus_g_p,
            input.woody_carbon_fraction,
            input.woody_nitrogen_fraction,
            input.woody_phosphorus_fraction,
            input.kinetics,
        ),
    };
}

fn addLitter(total: *metabolism.RootLitter, addition: metabolism.RootLitter) !void {
    inline for (@typeInfo(metabolism.RootLitter).@"struct".fields) |field| {
        for (&@field(total, field.name), @field(addition, field.name)) |*target, value| {
            if (!std.math.isFinite(target.*) or target.* < 0) return error.InvalidSecondaryRootLitterAccumulator;
            const next = target.* + value;
            if (!std.math.isFinite(next) or next < 0) return error.NonFiniteSecondaryRootLitterAccumulator;
            target.* = next;
        }
    }
}

/// Exact GROSUB lines 6210--6337 N,L,NR recycling, senescence, and four-
/// kinetic-pool litter sweep. All masses use g C, g N, or g P per biological
/// hour; recycling and senesced fractions are dimensionless.
pub fn apply(state: State, inputs: Inputs) !void {
    const domain_layer_count = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch
        return error.SecondaryRootSenescenceSweepDimensionOverflow;
    const value_count = std.math.mul(usize, domain_layer_count, inputs.root_axis_count) catch
        return error.SecondaryRootSenescenceSweepDimensionOverflow;
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or
        inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.active_by_domain_layer_axis.len != value_count or
        inputs.axis_inputs_by_domain_layer_axis.len != value_count or
        state.recycling_by_domain_layer_axis.len != value_count or
        state.senescence_by_domain_layer_axis.len != value_count or
        state.litter_by_domain_layer_axis.len != value_count or
        state.litter_by_domain_layer.len != domain_layer_count)
        return error.SecondaryRootSenescenceSweepDimensionMismatch;
    if (!std.math.isFinite(inputs.minimum_active_layer_thickness_m) or inputs.minimum_active_layer_thickness_m < 0)
        return error.InvalidSecondaryRootSenescenceSweepThreshold;
    for (inputs.layer_thickness_m) |thickness_m|
        if (!std.math.isFinite(thickness_m) or thickness_m < 0)
            return error.InvalidSecondaryRootSenescenceSweepLayerThickness;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count)
        return error.InvalidSecondaryRootSenescenceSweepLayerRange;

    // Preflight calculations and source-order CSNC/ZSNC/PSNC accumulation.
    for (0..inputs.biological_domain_count) |domain| {
        const deepest_layer = inputs.deepest_rooted_layer_index;
        for (inputs.planting_layer_index..deepest_layer + 1) |layer| {
            if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
            const domain_layer = domain * inputs.soil_layer_count + layer;
            var next_litter = state.litter_by_domain_layer[domain_layer];
            for (0..inputs.root_axis_count) |axis| {
                const index = domain_layer * inputs.root_axis_count + axis;
                if (!inputs.active_by_domain_layer_axis[index]) continue;
                const result = try calculate(inputs.parameters, inputs.axis_inputs_by_domain_layer_axis[index]);
                try addLitter(&next_litter, result.litter);
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
                const result = try calculate(inputs.parameters, inputs.axis_inputs_by_domain_layer_axis[index]);
                state.recycling_by_domain_layer_axis[index] = result.recycling;
                state.senescence_by_domain_layer_axis[index] = result.senescence;
                state.litter_by_domain_layer_axis[index] = result.litter;
                try addLitter(&state.litter_by_domain_layer[domain_layer], result.litter);
            }
        }
    }
}

fn zeroLitter() metabolism.RootLitter {
    return std.mem.zeroes(metabolism.RootLitter);
}

fn axisInput(root_carbon_g_c: f64) AxisInputs {
    const unit = [4]f64{ 0.1, 0.2, 0.3, 0.4 };
    return .{
        .emerged = true,
        .mobile_carbon_concentration_g_c_per_g_c = 0.1,
        .mobile_nitrogen_concentration_g_n_per_g_c = 0.02,
        .mobile_phosphorus_concentration_g_p_per_g_c = 0.002,
        .senescence = .{
            .oxygen_unlimited_substrate_minus_maintenance_g_c_per_h = -10,
            .actual_substrate_minus_maintenance_g_c_per_h = -10,
            .root_carbon_g_c = root_carbon_g_c,
            .root_nitrogen_g_n = root_carbon_g_c * 0.1,
            .root_phosphorus_g_p = root_carbon_g_c * 0.01,
            .oxygen_limitation = 1,
            .phenological_remobilization_enabled = false,
            .root_remobilization_enabled = false,
            .storage_exchange_fraction_per_h = 0,
            .remobilization_elapsed_h = 0,
            .full_senescence_h = 480,
            .biological_timestep_h = 1,
            .structural_presence_threshold_g_c = 1e-12,
        },
        .woody_carbon_fraction = .{ 0.25, 0.75 },
        .woody_nitrogen_fraction = .{ 0.25, 0.75 },
        .woody_phosphorus_fraction = .{ 0.25, 0.75 },
        .kinetics = .{
            .woody_carbon = unit,
            .woody_nitrogen = unit,
            .woody_phosphorus = unit,
            .nonwoody_carbon = unit,
            .nonwoody_nitrogen = unit,
            .nonwoody_phosphorus = unit,
        },
    };
}

test "GROSUB senescence litter sweep accumulates runtime axes in N L NR order" {
    const count = 12;
    var axis_inputs = [_]AxisInputs{axisInput(2)} ** count;
    for (&axis_inputs, 0..) |*input, axis| input.senescence.root_carbon_g_c = @as(f64, @floatFromInt(axis + 1));
    var recycling = [_]metabolism.RecyclingFractions{std.mem.zeroes(metabolism.RecyclingFractions)} ** count;
    var senescence = [_]metabolism.SecondaryRootSenescence{std.mem.zeroes(metabolism.SecondaryRootSenescence)} ** count;
    var axis_litter = [_]metabolism.RootLitter{zeroLitter()} ** count;
    var layer_litter = [_]metabolism.RootLitter{zeroLitter()} ** 4;
    try apply(.{ .recycling_by_domain_layer_axis = &recycling, .senescence_by_domain_layer_axis = &senescence, .litter_by_domain_layer_axis = &axis_litter, .litter_by_domain_layer = &layer_litter }, .{
        .biological_domain_count = 2,
        .soil_layer_count = 2,
        .root_axis_count = 3,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 1,
        .layer_thickness_m = &.{ 0.2, 0.2 },
        .minimum_active_layer_thickness_m = 0.001,
        .active_by_domain_layer_axis = &.{ true, true, false, true, true, true, true, true, true, true, true, true },
        .axis_inputs_by_domain_layer_axis = &axis_inputs,
        .parameters = metabolism.compatibilitySecondaryRootParameters(),
    });
    try std.testing.expectApproxEqAbs(
        axis_litter[3].woody_carbon_g_c[0] + axis_litter[4].woody_carbon_g_c[0] + axis_litter[5].woody_carbon_g_c[0],
        layer_litter[1].woody_carbon_g_c[0],
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        axis_litter[0].woody_carbon_g_c[0] + axis_litter[1].woody_carbon_g_c[0],
        layer_litter[0].woody_carbon_g_c[0],
        1e-15,
    );
}

test "GROSUB secondary litter plus recycled mobile mass closes C N P senescence" {
    const input = axisInput(8);
    const result = try calculate(metabolism.compatibilitySecondaryRootParameters(), input);
    inline for (.{
        .{ input.senescence.root_carbon_g_c, result.senescence.recyclable_carbon_g_c, result.litter.woody_carbon_g_c, result.litter.nonwoody_carbon_g_c, input.woody_carbon_fraction[1] },
        .{ input.senescence.root_nitrogen_g_n, result.senescence.recyclable_nitrogen_g_n, result.litter.woody_nitrogen_g_n, result.litter.nonwoody_nitrogen_g_n, input.woody_nitrogen_fraction[1] },
        .{ input.senescence.root_phosphorus_g_p, result.senescence.recyclable_phosphorus_g_p, result.litter.woody_phosphorus_g_p, result.litter.nonwoody_phosphorus_g_p, input.woody_phosphorus_fraction[1] },
    }) |balance| {
        var litter_total: f64 = 0;
        for (balance[2]) |value| litter_total += value;
        for (balance[3]) |value| litter_total += value;
        const recycled_to_mobile = result.senescence.senesced_fraction * balance[1] * balance[4];
        try std.testing.expectApproxEqAbs(
            result.senescence.senesced_fraction * balance[0],
            litter_total + recycled_to_mobile,
            1e-12,
        );
    }
}

test "GROSUB senescence litter sweep is atomic on invalid late axis" {
    var axis_inputs = [_]AxisInputs{axisInput(2)} ** 6;
    axis_inputs[5].senescence.root_carbon_g_c = std.math.nan(f64);
    var recycling = [_]metabolism.RecyclingFractions{.{ .carbon = 0.1, .nitrogen = 0.2, .phosphorus = 0.3 }} ** 6;
    var senescence = [_]metabolism.SecondaryRootSenescence{std.mem.zeroes(metabolism.SecondaryRootSenescence)} ** 6;
    var axis_litter = [_]metabolism.RootLitter{zeroLitter()} ** 6;
    var layer_litter = [_]metabolism.RootLitter{zeroLitter()} ** 2;
    const before_recycling = recycling;
    const before_layer_litter = layer_litter;
    try std.testing.expectError(error.NonFiniteSecondaryRootSenescenceInput, apply(.{ .recycling_by_domain_layer_axis = &recycling, .senescence_by_domain_layer_axis = &senescence, .litter_by_domain_layer_axis = &axis_litter, .litter_by_domain_layer = &layer_litter }, .{
        .biological_domain_count = 2,
        .soil_layer_count = 1,
        .root_axis_count = 3,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 0,
        .layer_thickness_m = &.{0.2},
        .minimum_active_layer_thickness_m = 0.001,
        .active_by_domain_layer_axis = &.{ true, true, true, true, true, true },
        .axis_inputs_by_domain_layer_axis = &axis_inputs,
        .parameters = metabolism.compatibilitySecondaryRootParameters(),
    }));
    try std.testing.expectEqualDeep(before_recycling, recycling);
    try std.testing.expectEqualDeep(before_layer_litter, layer_litter);
}
