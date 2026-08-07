const std = @import("std");
const metabolism = @import("plant_root_metabolism.zig");

pub const State = struct {
    mobile_carbon_g_c_by_domain_layer: []f64,
    mobile_nitrogen_g_n_by_domain_layer: []f64,
    mobile_phosphorus_g_p_by_domain_layer: []f64,
    protein_carbon_g_c_by_domain_layer: []f64,
    actual_respiration_g_c_by_domain_layer: []f64,
    oxygen_unlimited_respiration_g_c_by_domain_layer: []f64,
    carbon_unlimited_respiration_g_c_by_domain_layer: []f64,
    secondary_carbon_g_c_by_domain_layer_axis: []f64,
    secondary_nitrogen_g_n_by_domain_layer_axis: []f64,
    secondary_phosphorus_g_p_by_domain_layer_axis: []f64,
    secondary_length_m_by_domain_layer_axis: []f64,
    secondary_axis_count_by_domain_layer_axis: []f64,
    total_secondary_axis_count_by_domain_layer: []f64,
};

/// Caller-owned runtime memory used for rollback-safe sequential simulation.
pub const Workspace = State;

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    active_by_domain_layer_axis: []const bool,
    metabolism_by_domain_layer_axis: []const metabolism.SecondaryRootResult,
    senescence_by_domain_layer_axis: []const metabolism.SecondaryRootSenescence,
    specific_root_length_m_per_g_c_by_domain: []const f64,
    extension_water_response_by_domain_layer: []const f64,
    nonwoody_carbon_fraction_by_domain: []const f64,
    nonwoody_nitrogen_fraction_by_domain: []const f64,
    nonwoody_phosphorus_fraction_by_domain: []const f64,
    protein_carbon_per_nitrogen_g_c_per_g_n: f64,
    protein_carbon_per_phosphorus_g_c_per_g_p: f64,
    secondary_branching_per_m: f64,
    primary_axis_count_multiplier: f64,
};

fn validateStateShape(state: State, domain_layers: usize, values: usize) !void {
    inline for (.{
        state.mobile_carbon_g_c_by_domain_layer,
        state.mobile_nitrogen_g_n_by_domain_layer,
        state.mobile_phosphorus_g_p_by_domain_layer,
        state.protein_carbon_g_c_by_domain_layer,
        state.actual_respiration_g_c_by_domain_layer,
        state.oxygen_unlimited_respiration_g_c_by_domain_layer,
        state.carbon_unlimited_respiration_g_c_by_domain_layer,
        state.total_secondary_axis_count_by_domain_layer,
    }) |slice| if (slice.len != domain_layers) return error.SecondaryRootStateUpdateDimensionMismatch;
    inline for (.{
        state.secondary_carbon_g_c_by_domain_layer_axis,
        state.secondary_nitrogen_g_n_by_domain_layer_axis,
        state.secondary_phosphorus_g_p_by_domain_layer_axis,
        state.secondary_length_m_by_domain_layer_axis,
        state.secondary_axis_count_by_domain_layer_axis,
    }) |slice| if (slice.len != values) return error.SecondaryRootStateUpdateDimensionMismatch;
}

fn copyState(target: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        @memcpy(@field(target, field.name), @field(source, field.name));
}

fn validateFiniteNonnegative(state: State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| for (@field(state, field.name)) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootState;
}

fn updateAxis(state: State, inputs: Inputs, domain: usize, layer: usize, axis: usize) !void {
    const domain_layer = domain * inputs.soil_layer_count + layer;
    const index = domain_layer * inputs.root_axis_count + axis;
    const process = inputs.metabolism_by_domain_layer_axis[index];
    const senescence = inputs.senescence_by_domain_layer_axis[index];
    inline for (@typeInfo(metabolism.SecondaryRootResult).@"struct".fields) |field| {
        const value = @field(process, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootMetabolismResult;
    }
    inline for (@typeInfo(metabolism.SecondaryRootSenescence).@"struct".fields) |field| {
        const value = @field(senescence, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootSenescenceResult;
    }
    if (senescence.senesced_fraction > 1) return error.InvalidSecondaryRootSenescenceResult;

    const old_c = state.secondary_carbon_g_c_by_domain_layer_axis[index];
    const old_n = state.secondary_nitrogen_g_n_by_domain_layer_axis[index];
    const old_p = state.secondary_phosphorus_g_p_by_domain_layer_axis[index];
    const old_length = state.secondary_length_m_by_domain_layer_axis[index];
    const fraction = senescence.senesced_fraction;
    const recovered_c = fraction * senescence.recyclable_carbon_g_c * inputs.nonwoody_carbon_fraction_by_domain[domain];
    const recovered_n = fraction * senescence.recyclable_nitrogen_g_n * inputs.nonwoody_nitrogen_fraction_by_domain[domain];
    const recovered_p = fraction * senescence.recyclable_phosphorus_g_p * inputs.nonwoody_phosphorus_fraction_by_domain[domain];
    const respiration = try metabolism.assemble(.{
        .maintenance_demand_g_c = process.maintenance_respiration_g_c_per_h,
        .substrate_respiration_actual_g_c = process.substrate_respiration_actual_g_c_per_h,
        .substrate_respiration_oxygen_unlimited_g_c = process.substrate_respiration_oxygen_unlimited_g_c_per_h,
        .growth_respiration_actual_g_c = process.growth_respiration_actual_g_c_per_h,
        .growth_respiration_oxygen_unlimited_g_c = process.growth_respiration_oxygen_unlimited_g_c_per_h,
        .senescence_respiration_actual_g_c = senescence.respiration_actual_g_c_per_h,
        .senescence_respiration_oxygen_unlimited_g_c = senescence.respiration_oxygen_unlimited_g_c_per_h,
        .nitrogen_assimilation_respiration_actual_g_c = process.nitrogen_assimilation_respiration_actual_g_c_per_h,
        .nitrogen_assimilation_respiration_oxygen_unlimited_g_c = process.nitrogen_assimilation_respiration_oxygen_unlimited_g_c_per_h,
    });

    state.mobile_carbon_g_c_by_domain_layer[domain_layer] += -@min(process.maintenance_respiration_g_c_per_h, process.substrate_respiration_actual_g_c_per_h) - process.growth_and_respiration_carbon_actual_g_c_per_h - process.nitrogen_assimilation_respiration_actual_g_c_per_h - senescence.respiration_actual_g_c_per_h + recovered_c;
    state.mobile_nitrogen_g_n_by_domain_layer[domain_layer] += -process.nitrogen_growth_actual_g_n_per_h + recovered_n;
    state.mobile_phosphorus_g_p_by_domain_layer[domain_layer] += -process.phosphorus_growth_actual_g_p_per_h + recovered_p;
    const next_c = old_c + process.root_growth_actual_g_c_per_h - fraction * old_c;
    const next_n = old_n + process.nitrogen_growth_actual_g_n_per_h - fraction * old_n;
    const next_p = old_p + process.phosphorus_growth_actual_g_p_per_h - fraction * old_p;
    const next_length = old_length + process.root_growth_actual_g_c_per_h * inputs.specific_root_length_m_per_g_c_by_domain[domain] * inputs.extension_water_response_by_domain_layer[domain_layer] - fraction * old_length;
    state.secondary_carbon_g_c_by_domain_layer_axis[index] = next_c;
    state.secondary_nitrogen_g_n_by_domain_layer_axis[index] = next_n;
    state.secondary_phosphorus_g_p_by_domain_layer_axis[index] = next_p;
    state.secondary_length_m_by_domain_layer_axis[index] = next_length;
    state.protein_carbon_g_c_by_domain_layer[domain_layer] += @min(inputs.protein_carbon_per_nitrogen_g_c_per_g_n * next_n, inputs.protein_carbon_per_phosphorus_g_c_per_g_p * next_p);
    state.actual_respiration_g_c_by_domain_layer[domain_layer] += respiration.actual_g_c;
    state.oxygen_unlimited_respiration_g_c_by_domain_layer[domain_layer] += respiration.oxygen_unlimited_g_c;
    state.carbon_unlimited_respiration_g_c_by_domain_layer[domain_layer] += respiration.carbon_unlimited_g_c;
    const first_order = inputs.secondary_branching_per_m * inputs.primary_axis_count_multiplier;
    const axis_count = (first_order + inputs.secondary_branching_per_m * first_order) * inputs.layer_thickness_m[layer];
    state.secondary_axis_count_by_domain_layer_axis[index] = axis_count;
    state.total_secondary_axis_count_by_domain_layer[domain_layer] += axis_count;
    try validateFiniteNonnegative(state);
}

/// Exact grosub.f lines 6356--6427 sequential N,L,NR state transaction.
/// Masses are g C/N/P, root length is m, and axis counts are dimensionless.
pub fn apply(state: State, workspace: Workspace, inputs: Inputs) !void {
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.SecondaryRootStateUpdateDimensionOverflow;
    const values = std.math.mul(usize, domain_layers, inputs.root_axis_count) catch return error.SecondaryRootStateUpdateDimensionOverflow;
    try validateStateShape(state, domain_layers, values);
    try validateStateShape(workspace, domain_layers, values);
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or
        inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.active_by_domain_layer_axis.len != values or inputs.metabolism_by_domain_layer_axis.len != values or inputs.senescence_by_domain_layer_axis.len != values or
        inputs.specific_root_length_m_per_g_c_by_domain.len != inputs.biological_domain_count or inputs.extension_water_response_by_domain_layer.len != domain_layers or
        inputs.nonwoody_carbon_fraction_by_domain.len != inputs.biological_domain_count or inputs.nonwoody_nitrogen_fraction_by_domain.len != inputs.biological_domain_count or inputs.nonwoody_phosphorus_fraction_by_domain.len != inputs.biological_domain_count)
        return error.SecondaryRootStateUpdateDimensionMismatch;
    inline for (.{ inputs.minimum_active_layer_thickness_m, inputs.protein_carbon_per_nitrogen_g_c_per_g_n, inputs.protein_carbon_per_phosphorus_g_c_per_g_p, inputs.secondary_branching_per_m, inputs.primary_axis_count_multiplier }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootStateUpdateInput;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count) return error.InvalidSecondaryRootStateUpdateLayerRange;
    for (inputs.layer_thickness_m) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootStateUpdateInput;
    for (inputs.specific_root_length_m_per_g_c_by_domain) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootStateUpdateInput;
    inline for (.{ inputs.extension_water_response_by_domain_layer, inputs.nonwoody_carbon_fraction_by_domain, inputs.nonwoody_nitrogen_fraction_by_domain, inputs.nonwoody_phosphorus_fraction_by_domain }) |slice| for (slice) |value|
        if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidSecondaryRootStateUpdateInput;
    try validateFiniteNonnegative(state);
    copyState(workspace, state);
    for (0..inputs.biological_domain_count) |domain| for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
        if (!(inputs.layer_thickness_m[layer] > inputs.minimum_active_layer_thickness_m)) continue;
        for (0..inputs.root_axis_count) |axis| {
            const index = (domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
            if (inputs.active_by_domain_layer_axis[index]) try updateAxis(workspace, inputs, domain, layer, axis);
        }
    };
    copyState(state, workspace);
}

fn testProcess() metabolism.SecondaryRootResult {
    return .{
        .nutrient_feedback = 1,
        .substrate_respiration_oxygen_unlimited_g_c_per_h = 0.6,
        .maintenance_respiration_g_c_per_h = 0.1,
        .substrate_respiration_actual_g_c_per_h = 0.5,
        .growth_respiration_oxygen_unlimited_g_c_per_h = 0.25,
        .growth_respiration_actual_g_c_per_h = 0.2,
        .growth_and_respiration_carbon_oxygen_unlimited_g_c_per_h = 1.0,
        .growth_and_respiration_carbon_actual_g_c_per_h = 0.8,
        .root_growth_oxygen_unlimited_g_c_per_h = 0.75,
        .root_growth_actual_g_c_per_h = 0.6,
        .nitrogen_growth_demand_g_n_per_h = 0.075,
        .nitrogen_growth_actual_g_n_per_h = 0.06,
        .phosphorus_growth_actual_g_p_per_h = 0.006,
        .nitrogen_assimilation_respiration_oxygen_unlimited_g_c_per_h = 0.1275,
        .nitrogen_assimilation_respiration_actual_g_c_per_h = 0.102,
    };
}

fn testSenescence() metabolism.SecondaryRootSenescence {
    return .{
        .respiration_oxygen_unlimited_g_c_per_h = 0.06,
        .respiration_actual_g_c_per_h = 0.05,
        .phenological_senescence_g_c_per_h = 0,
        .senesced_fraction = 0.1,
        .recyclable_carbon_g_c = 1,
        .recyclable_nitrogen_g_n = 0.1,
        .recyclable_phosphorus_g_p = 0.01,
    };
}

fn testInputs(processes: []const metabolism.SecondaryRootResult, losses: []const metabolism.SecondaryRootSenescence) Inputs {
    return .{
        .biological_domain_count = 1,
        .soil_layer_count = 1,
        .root_axis_count = 2,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 0,
        .layer_thickness_m = &.{0.2},
        .minimum_active_layer_thickness_m = 0.001,
        .active_by_domain_layer_axis = &.{ true, true },
        .metabolism_by_domain_layer_axis = processes,
        .senescence_by_domain_layer_axis = losses,
        .specific_root_length_m_per_g_c_by_domain = &.{2},
        .extension_water_response_by_domain_layer = &.{0.5},
        .nonwoody_carbon_fraction_by_domain = &.{0.5},
        .nonwoody_nitrogen_fraction_by_domain = &.{0.5},
        .nonwoody_phosphorus_fraction_by_domain = &.{0.5},
        .protein_carbon_per_nitrogen_g_c_per_g_n = 2,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 20,
        .secondary_branching_per_m = 0.5,
        .primary_axis_count_multiplier = 2,
    };
}

fn makeState(
    mobile_c: []f64,
    mobile_n: []f64,
    mobile_p: []f64,
    protein: []f64,
    actual: []f64,
    oxygen: []f64,
    carbon: []f64,
    structural_c: []f64,
    structural_n: []f64,
    structural_p: []f64,
    length: []f64,
    axis_count: []f64,
    total_axis_count: []f64,
) State {
    return .{
        .mobile_carbon_g_c_by_domain_layer = mobile_c,
        .mobile_nitrogen_g_n_by_domain_layer = mobile_n,
        .mobile_phosphorus_g_p_by_domain_layer = mobile_p,
        .protein_carbon_g_c_by_domain_layer = protein,
        .actual_respiration_g_c_by_domain_layer = actual,
        .oxygen_unlimited_respiration_g_c_by_domain_layer = oxygen,
        .carbon_unlimited_respiration_g_c_by_domain_layer = carbon,
        .secondary_carbon_g_c_by_domain_layer_axis = structural_c,
        .secondary_nitrogen_g_n_by_domain_layer_axis = structural_n,
        .secondary_phosphorus_g_p_by_domain_layer_axis = structural_p,
        .secondary_length_m_by_domain_layer_axis = length,
        .secondary_axis_count_by_domain_layer_axis = axis_count,
        .total_secondary_axis_count_by_domain_layer = total_axis_count,
    };
}

test "GROSUB sequential secondary update closes C N P and axis totals" {
    var mobile_c = [_]f64{5};
    var mobile_n = [_]f64{1};
    var mobile_p = [_]f64{0.1};
    var protein = [_]f64{0};
    var actual = [_]f64{0};
    var oxygen = [_]f64{0};
    var carbon = [_]f64{0};
    var structural_c = [_]f64{ 2, 2 };
    var structural_n = [_]f64{ 0.2, 0.2 };
    var structural_p = [_]f64{ 0.02, 0.02 };
    var length = [_]f64{ 1, 1 };
    var axis_count = [_]f64{ 0, 0 };
    var total_axis_count = [_]f64{0};
    const live = makeState(&mobile_c, &mobile_n, &mobile_p, &protein, &actual, &oxygen, &carbon, &structural_c, &structural_n, &structural_p, &length, &axis_count, &total_axis_count);
    var workspace_storage: [13][2]f64 = std.mem.zeroes([13][2]f64);
    const work = makeState(workspace_storage[0][0..1], workspace_storage[1][0..1], workspace_storage[2][0..1], workspace_storage[3][0..1], workspace_storage[4][0..1], workspace_storage[5][0..1], workspace_storage[6][0..1], &workspace_storage[7], &workspace_storage[8], &workspace_storage[9], &workspace_storage[10], &workspace_storage[11], workspace_storage[12][0..1]);
    const processes = [_]metabolism.SecondaryRootResult{ testProcess(), testProcess() };
    const losses = [_]metabolism.SecondaryRootSenescence{ testSenescence(), testSenescence() };
    const initial_c = mobile_c[0] + structural_c[0] + structural_c[1];
    const initial_n = mobile_n[0] + structural_n[0] + structural_n[1];
    const initial_p = mobile_p[0] + structural_p[0] + structural_p[1];
    try apply(live, work, testInputs(&processes, &losses));
    const litter_c = 2 * (0.1 * 2 - 0.1 * 1 * 0.5);
    const litter_n = 2 * (0.1 * 0.2 - 0.1 * 0.1 * 0.5);
    const litter_p = 2 * (0.1 * 0.02 - 0.1 * 0.01 * 0.5);
    try std.testing.expectApproxEqAbs(initial_c, mobile_c[0] + structural_c[0] + structural_c[1] + actual[0] + litter_c, 1e-12);
    try std.testing.expectApproxEqAbs(initial_n, mobile_n[0] + structural_n[0] + structural_n[1] + litter_n, 1e-12);
    try std.testing.expectApproxEqAbs(initial_p, mobile_p[0] + structural_p[0] + structural_p[1] + litter_p, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), axis_count[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), total_axis_count[0], 1e-15);
}

test "GROSUB sequential shared-pool overdraw rolls back every axis" {
    var arrays: [13][2]f64 = std.mem.zeroes([13][2]f64);
    arrays[0][0] = 1.5;
    arrays[1][0] = 1;
    arrays[2][0] = 0.1;
    arrays[7] = .{ 2, 2 };
    arrays[8] = .{ 0.2, 0.2 };
    arrays[9] = .{ 0.02, 0.02 };
    arrays[10] = .{ 1, 1 };
    const live = makeState(arrays[0][0..1], arrays[1][0..1], arrays[2][0..1], arrays[3][0..1], arrays[4][0..1], arrays[5][0..1], arrays[6][0..1], &arrays[7], &arrays[8], &arrays[9], &arrays[10], &arrays[11], arrays[12][0..1]);
    var workspace_arrays: [13][2]f64 = std.mem.zeroes([13][2]f64);
    const work = makeState(workspace_arrays[0][0..1], workspace_arrays[1][0..1], workspace_arrays[2][0..1], workspace_arrays[3][0..1], workspace_arrays[4][0..1], workspace_arrays[5][0..1], workspace_arrays[6][0..1], &workspace_arrays[7], &workspace_arrays[8], &workspace_arrays[9], &workspace_arrays[10], &workspace_arrays[11], workspace_arrays[12][0..1]);
    const before = arrays;
    const processes = [_]metabolism.SecondaryRootResult{ testProcess(), testProcess() };
    const losses = [_]metabolism.SecondaryRootSenescence{ testSenescence(), testSenescence() };
    try std.testing.expectError(error.InvalidSecondaryRootState, apply(live, work, testInputs(&processes, &losses)));
    try std.testing.expectEqualDeep(before, arrays);
}
