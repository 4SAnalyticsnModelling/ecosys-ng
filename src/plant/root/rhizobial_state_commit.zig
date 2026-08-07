const std = @import("std");

pub const State = struct {
    respiration_oxygen_unlimited_g_c_by_layer: []f64,
    respiration_substrate_limited_g_c_by_layer: []f64,
    respiration_actual_g_c_by_layer: []f64,
    litter_carbon_g_c_by_layer_pool: []f64,
    litter_nitrogen_g_n_by_layer_pool: []f64,
    litter_phosphorus_g_p_by_layer_pool: []f64,
    mobile_carbon_g_c_by_layer: []f64,
    mobile_nitrogen_g_n_by_layer: []f64,
    mobile_phosphorus_g_p_by_layer: []f64,
    structural_carbon_g_c_by_layer: []f64,
    structural_nitrogen_g_n_by_layer: []f64,
    structural_phosphorus_g_p_by_layer: []f64,
};

pub const LitterFractions = struct {
    carbon: []const f64,
    nitrogen: []const f64,
    phosphorus: []const f64,
};

pub const Inputs = struct {
    soil_layer_count: usize,
    kinetic_pool_count: usize,
    /// Source `NU`, the first active soil layer (not planting layer `NG`).
    planting_layer_index: usize,
    plant_deepest_rooted_layer_index: usize,
    fixation_type: i32,
    layer_active: []const bool,
    maintenance_respiration_g_c_by_layer: []const f64,
    substrate_respiration_oxygen_unlimited_g_c_by_layer: []const f64,
    substrate_respiration_g_c_by_layer: []const f64,
    growth_respiration_oxygen_unlimited_g_c_by_layer: []const f64,
    growth_respiration_g_c_by_layer: []const f64,
    senescence_recycled_carbon_g_c_by_layer: []const f64,
    decomposition_litter_carbon_g_c_by_layer: []const f64,
    decomposition_litter_nitrogen_g_n_by_layer: []const f64,
    decomposition_litter_phosphorus_g_p_by_layer: []const f64,
    senescence_litter_carbon_g_c_by_layer: []const f64,
    senescence_litter_nitrogen_g_n_by_layer: []const f64,
    senescence_litter_phosphorus_g_p_by_layer: []const f64,
    fixation_respiration_g_c_by_layer: []const f64,
    growth_substrate_carbon_g_c_by_layer: []const f64,
    decomposition_recycled_carbon_g_c_by_layer: []const f64,
    decomposition_recycled_nitrogen_g_n_by_layer: []const f64,
    decomposition_recycled_phosphorus_g_p_by_layer: []const f64,
    senescence_recycled_nitrogen_g_n_by_layer: []const f64,
    senescence_recycled_phosphorus_g_p_by_layer: []const f64,
    growth_nitrogen_g_n_by_layer: []const f64,
    growth_phosphorus_g_p_by_layer: []const f64,
    nitrogen_fixation_g_n_by_layer: []const f64,
    structural_growth_carbon_g_c_by_layer: []const f64,
    decomposition_loss_carbon_g_c_by_layer: []const f64,
    decomposition_loss_nitrogen_g_n_by_layer: []const f64,
    decomposition_loss_phosphorus_g_p_by_layer: []const f64,
    senescence_loss_carbon_g_c_by_layer: []const f64,
    senescence_loss_nitrogen_g_n_by_layer: []const f64,
    senescence_loss_phosphorus_g_p_by_layer: []const f64,
    litter_fractions: LitterFractions,
};

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateState(state: State, inputs: Inputs) !void {
    const litter_count = std.math.mul(usize, inputs.soil_layer_count, inputs.kinetic_pool_count) catch return error.RootRhizobialCommitDimensionOverflow;
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        const expected = if (std.mem.startsWith(u8, field.name, "litter_")) litter_count else inputs.soil_layer_count;
        if (values.len != expected) return error.RootRhizobialCommitDimensionMismatch;
        const signed = std.mem.eql(u8, field.name, "respiration_actual_g_c_by_layer");
        for (values) |value| if (!std.math.isFinite(value) or (!signed and value < 0)) return error.InvalidRootRhizobialCommitState;
    }
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.soil_layer_count == 0 or inputs.kinetic_pool_count == 0 or inputs.planting_layer_index > inputs.plant_deepest_rooted_layer_index or inputs.plant_deepest_rooted_layer_index >= inputs.soil_layer_count or inputs.layer_active.len != inputs.soil_layer_count) return error.RootRhizobialCommitDimensionMismatch;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != inputs.soil_layer_count) return error.RootRhizobialCommitDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value)) return error.InvalidRootRhizobialCommitInput;
    };
    inline for (@typeInfo(LitterFractions).@"struct".fields) |field| {
        const values = @field(inputs.litter_fractions, field.name);
        if (values.len != inputs.kinetic_pool_count) return error.RootRhizobialCommitDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootRhizobialCommitInput;
    }
}

/// Exact GROSUB 7848--7907 rhizobial respiration publication, runtime-pool
/// litter allocation, and mobile/structural C:N:P commit. Traversal is L then M.
/// Fluxes and pools use g C/N/P per biological timestep.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateState(state, inputs);
    try validateState(workspace, inputs);
    try validateInputs(inputs);
    copyState(workspace, state);
    if (inputs.fixation_type >= 1 and inputs.fixation_type <= 3) for (inputs.planting_layer_index..inputs.plant_deepest_rooted_layer_index + 1) |layer| {
        if (!inputs.layer_active[layer]) continue;
        const unlimited = @min(inputs.maintenance_respiration_g_c_by_layer[layer], inputs.substrate_respiration_oxygen_unlimited_g_c_by_layer[layer]) + inputs.growth_respiration_oxygen_unlimited_g_c_by_layer[layer] + inputs.senescence_recycled_carbon_g_c_by_layer[layer];
        const actual = @min(inputs.maintenance_respiration_g_c_by_layer[layer], inputs.substrate_respiration_g_c_by_layer[layer]) + inputs.growth_respiration_g_c_by_layer[layer] + inputs.senescence_recycled_carbon_g_c_by_layer[layer];
        workspace.respiration_oxygen_unlimited_g_c_by_layer[layer] += unlimited;
        workspace.respiration_substrate_limited_g_c_by_layer[layer] += actual;
        workspace.respiration_actual_g_c_by_layer[layer] -= actual;
        for (0..inputs.kinetic_pool_count) |pool| {
            const index = layer * inputs.kinetic_pool_count + pool;
            workspace.litter_carbon_g_c_by_layer_pool[index] += inputs.litter_fractions.carbon[pool] * (inputs.decomposition_litter_carbon_g_c_by_layer[layer] + inputs.senescence_litter_carbon_g_c_by_layer[layer]);
            workspace.litter_nitrogen_g_n_by_layer_pool[index] += inputs.litter_fractions.nitrogen[pool] * (inputs.decomposition_litter_nitrogen_g_n_by_layer[layer] + inputs.senescence_litter_nitrogen_g_n_by_layer[layer]);
            workspace.litter_phosphorus_g_p_by_layer_pool[index] += inputs.litter_fractions.phosphorus[pool] * (inputs.decomposition_litter_phosphorus_g_p_by_layer[layer] + inputs.senescence_litter_phosphorus_g_p_by_layer[layer]);
        }
        workspace.mobile_carbon_g_c_by_layer[layer] = workspace.mobile_carbon_g_c_by_layer[layer] - @min(inputs.maintenance_respiration_g_c_by_layer[layer], inputs.substrate_respiration_g_c_by_layer[layer]) - inputs.fixation_respiration_g_c_by_layer[layer] - inputs.growth_substrate_carbon_g_c_by_layer[layer] + inputs.decomposition_recycled_carbon_g_c_by_layer[layer];
        workspace.mobile_nitrogen_g_n_by_layer[layer] = workspace.mobile_nitrogen_g_n_by_layer[layer] - inputs.growth_nitrogen_g_n_by_layer[layer] + inputs.decomposition_recycled_nitrogen_g_n_by_layer[layer] + inputs.senescence_recycled_nitrogen_g_n_by_layer[layer] + inputs.nitrogen_fixation_g_n_by_layer[layer];
        workspace.mobile_phosphorus_g_p_by_layer[layer] = workspace.mobile_phosphorus_g_p_by_layer[layer] - inputs.growth_phosphorus_g_p_by_layer[layer] + inputs.decomposition_recycled_phosphorus_g_p_by_layer[layer] + inputs.senescence_recycled_phosphorus_g_p_by_layer[layer];
        workspace.structural_carbon_g_c_by_layer[layer] += inputs.structural_growth_carbon_g_c_by_layer[layer] - inputs.decomposition_loss_carbon_g_c_by_layer[layer] - inputs.senescence_loss_carbon_g_c_by_layer[layer];
        workspace.structural_nitrogen_g_n_by_layer[layer] += inputs.growth_nitrogen_g_n_by_layer[layer] - inputs.decomposition_loss_nitrogen_g_n_by_layer[layer] - inputs.senescence_loss_nitrogen_g_n_by_layer[layer];
        workspace.structural_phosphorus_g_p_by_layer[layer] += inputs.growth_phosphorus_g_p_by_layer[layer] - inputs.decomposition_loss_phosphorus_g_p_by_layer[layer] - inputs.senescence_loss_phosphorus_g_p_by_layer[layer];
    };
    try validateState(workspace, inputs);
    copyState(state, workspace);
}

fn makeState(layer_values: *[9][1]f64, litter_values: *[3][3]f64) State {
    return .{ .respiration_oxygen_unlimited_g_c_by_layer = &layer_values[0], .respiration_substrate_limited_g_c_by_layer = &layer_values[1], .respiration_actual_g_c_by_layer = &layer_values[2], .litter_carbon_g_c_by_layer_pool = &litter_values[0], .litter_nitrogen_g_n_by_layer_pool = &litter_values[1], .litter_phosphorus_g_p_by_layer_pool = &litter_values[2], .mobile_carbon_g_c_by_layer = &layer_values[3], .mobile_nitrogen_g_n_by_layer = &layer_values[4], .mobile_phosphorus_g_p_by_layer = &layer_values[5], .structural_carbon_g_c_by_layer = &layer_values[6], .structural_nitrogen_g_n_by_layer = &layer_values[7], .structural_phosphorus_g_p_by_layer = &layer_values[8] };
}

fn testInputs() Inputs {
    return .{ .soil_layer_count = 1, .kinetic_pool_count = 3, .planting_layer_index = 0, .plant_deepest_rooted_layer_index = 0, .fixation_type = 2, .layer_active = &.{true}, .maintenance_respiration_g_c_by_layer = &.{0.2}, .substrate_respiration_oxygen_unlimited_g_c_by_layer = &.{1}, .substrate_respiration_g_c_by_layer = &.{0.8}, .growth_respiration_oxygen_unlimited_g_c_by_layer = &.{0.5}, .growth_respiration_g_c_by_layer = &.{0.4}, .senescence_recycled_carbon_g_c_by_layer = &.{0.1}, .decomposition_litter_carbon_g_c_by_layer = &.{0.6}, .decomposition_litter_nitrogen_g_n_by_layer = &.{0.3}, .decomposition_litter_phosphorus_g_p_by_layer = &.{0.15}, .senescence_litter_carbon_g_c_by_layer = &.{0.3}, .senescence_litter_nitrogen_g_n_by_layer = &.{0.15}, .senescence_litter_phosphorus_g_p_by_layer = &.{0.075}, .fixation_respiration_g_c_by_layer = &.{0.1}, .growth_substrate_carbon_g_c_by_layer = &.{0.5}, .decomposition_recycled_carbon_g_c_by_layer = &.{0.2}, .decomposition_recycled_nitrogen_g_n_by_layer = &.{0.1}, .decomposition_recycled_phosphorus_g_p_by_layer = &.{0.05}, .senescence_recycled_nitrogen_g_n_by_layer = &.{0.05}, .senescence_recycled_phosphorus_g_p_by_layer = &.{0.025}, .growth_nitrogen_g_n_by_layer = &.{0.2}, .growth_phosphorus_g_p_by_layer = &.{0.1}, .nitrogen_fixation_g_n_by_layer = &.{0.04}, .structural_growth_carbon_g_c_by_layer = &.{0.5}, .decomposition_loss_carbon_g_c_by_layer = &.{1}, .decomposition_loss_nitrogen_g_n_by_layer = &.{0.5}, .decomposition_loss_phosphorus_g_p_by_layer = &.{0.25}, .senescence_loss_carbon_g_c_by_layer = &.{0.5}, .senescence_loss_nitrogen_g_n_by_layer = &.{0.25}, .senescence_loss_phosphorus_g_p_by_layer = &.{0.125}, .litter_fractions = .{ .carbon = &.{ 0.2, 0.3, 0.5 }, .nitrogen = &.{ 0.2, 0.3, 0.5 }, .phosphorus = &.{ 0.2, 0.3, 0.5 } } };
}

test "GROSUB rhizobial commit routes runtime litter pools and closes respiration" {
    var layers: [9][1]f64 = .{ .{0}, .{0}, .{0}, .{5}, .{2}, .{1}, .{10}, .{2}, .{1} };
    var litter: [3][3]f64 = std.mem.zeroes([3][3]f64);
    var work_layers: [9][1]f64 = std.mem.zeroes([9][1]f64);
    var work_litter: [3][3]f64 = std.mem.zeroes([3][3]f64);
    try apply(makeState(&layers, &litter), makeState(&work_layers, &work_litter), testInputs());
    try std.testing.expectApproxEqAbs(0.8, layers[0][0], 1e-12);
    try std.testing.expectApproxEqAbs(0.7, layers[1][0], 1e-12);
    try std.testing.expectApproxEqAbs(-0.7, layers[2][0], 1e-12);
    var litter_c: f64 = 0;
    for (litter[0]) |value| litter_c += value;
    try std.testing.expectApproxEqAbs(0.9, litter_c, 1e-12);
    try std.testing.expectApproxEqAbs(4.4, layers[3][0], 1e-12);
    try std.testing.expectApproxEqAbs(9, layers[6][0], 1e-12);
}

test "GROSUB rhizobial commit rolls back structural overdraw" {
    var layers: [9][1]f64 = .{ .{0}, .{0}, .{0}, .{5}, .{2}, .{1}, .{0.1}, .{0.1}, .{0.1} };
    const before = layers;
    var litter: [3][3]f64 = std.mem.zeroes([3][3]f64);
    var work_layers: [9][1]f64 = std.mem.zeroes([9][1]f64);
    var work_litter: [3][3]f64 = std.mem.zeroes([3][3]f64);
    try std.testing.expectError(error.InvalidRootRhizobialCommitState, apply(makeState(&layers, &litter), makeState(&work_layers, &work_litter), testInputs()));
    try std.testing.expectEqualDeep(before, layers);
}
