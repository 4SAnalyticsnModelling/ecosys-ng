const std = @import("std");

pub const State = struct {
    gaseous_mass_g_by_domain_layer_species: []f64,
    aqueous_mass_g_by_domain_layer_species: []f64,
    withdrawal_ledger_g_by_species: []f64,
    total_secondary_axis_count_by_domain_layer: []f64,
    secondary_axis_count_by_domain_layer_axis: []f64,
    primary_axis_count_by_domain_layer: []f64,
    primary_length_m_by_domain_layer_axis: []f64,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    gas_species_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    trigger_active_by_domain_layer_axis: []const bool,
    deepest_rooted_layer_index_by_axis: []const usize,
    retracted_depth_m_by_domain_axis: []const f64,
    seeding_depth_m: f64,
    depth_above_planting_layer_bottom_m: f64,
    layer_bottom_depth_m: []const f64,
    layer_thickness_m: []const f64,
    minimum_active_layer_thickness_m: f64,
    axis_primary_sink_strength_m_by_domain_layer_axis: []const f64,
    axis_secondary_sink_strength_m_by_domain_layer_axis: []const f64,
    total_root_sink_strength_m_by_domain_layer: []const f64,
    sink_presence_threshold_m: f64,
    primary_axis_withdrawal_count_by_domain_layer_axis: []const f64,
};

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateState(model_state: State, species_values: usize, domain_layers: usize, axis_layers: usize, species: usize) !void {
    inline for (.{ model_state.gaseous_mass_g_by_domain_layer_species, model_state.aqueous_mass_g_by_domain_layer_species }) |values| {
        if (values.len != species_values) return error.PrimaryRootWithdrawalDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootWithdrawalState;
    }
    if (model_state.withdrawal_ledger_g_by_species.len != species) return error.PrimaryRootWithdrawalDimensionMismatch;
    for (model_state.withdrawal_ledger_g_by_species) |value| if (!std.math.isFinite(value)) return error.InvalidPrimaryRootWithdrawalState;
    inline for (.{ model_state.total_secondary_axis_count_by_domain_layer, model_state.primary_axis_count_by_domain_layer }) |values| {
        if (values.len != domain_layers) return error.PrimaryRootWithdrawalDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootWithdrawalState;
    }
    inline for (.{ model_state.secondary_axis_count_by_domain_layer_axis, model_state.primary_length_m_by_domain_layer_axis }) |values| {
        if (values.len != axis_layers) return error.PrimaryRootWithdrawalDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootWithdrawalState;
    }
}

fn validateInputs(inputs: Inputs, domain_layers: usize, axis_layers: usize, domain_axes: usize) !void {
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or inputs.gas_species_count == 0 or inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count or
        inputs.trigger_active_by_domain_layer_axis.len != axis_layers or inputs.deepest_rooted_layer_index_by_axis.len != inputs.root_axis_count or inputs.retracted_depth_m_by_domain_axis.len != domain_axes or inputs.layer_bottom_depth_m.len != inputs.soil_layer_count or inputs.layer_thickness_m.len != inputs.soil_layer_count or
        inputs.axis_primary_sink_strength_m_by_domain_layer_axis.len != axis_layers or inputs.axis_secondary_sink_strength_m_by_domain_layer_axis.len != axis_layers or inputs.total_root_sink_strength_m_by_domain_layer.len != domain_layers or inputs.primary_axis_withdrawal_count_by_domain_layer_axis.len != axis_layers)
        return error.PrimaryRootWithdrawalDimensionMismatch;
    inline for (.{ inputs.seeding_depth_m, inputs.depth_above_planting_layer_bottom_m, inputs.minimum_active_layer_thickness_m, inputs.sink_presence_threshold_m }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootWithdrawalInput;
    for (inputs.deepest_rooted_layer_index_by_axis) |layer| if (layer < inputs.planting_layer_index or layer >= inputs.soil_layer_count) return error.InvalidPrimaryRootWithdrawalInput;
    inline for (.{ inputs.retracted_depth_m_by_domain_axis, inputs.layer_bottom_depth_m, inputs.layer_thickness_m, inputs.axis_primary_sink_strength_m_by_domain_layer_axis, inputs.axis_secondary_sink_strength_m_by_domain_layer_axis, inputs.total_root_sink_strength_m_by_domain_layer, inputs.primary_axis_withdrawal_count_by_domain_layer_axis }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootWithdrawalInput;
}

fn layerIndex(inputs: Inputs, domain: usize, layer: usize) usize {
    return domain * inputs.soil_layer_count + layer;
}
fn axisLayerIndex(inputs: Inputs, domain: usize, layer: usize, axis: usize) usize {
    return (domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
}
fn speciesIndex(inputs: Inputs, domain: usize, layer: usize, species: usize) usize {
    return (domain * inputs.soil_layer_count + layer) * inputs.gas_species_count + species;
}

fn applyWithdrawalLayer(model_state: State, inputs: Inputs, controlling_domain: usize, trigger_layer: usize, layer: usize, axis: usize) !void {
    const controller_axis_layer = axisLayerIndex(inputs, controlling_domain, layer, axis);
    const controller_layer = layerIndex(inputs, controlling_domain, layer);
    const total_sink = inputs.total_root_sink_strength_m_by_domain_layer[controller_layer];
    const fraction = if (total_sink > inputs.sink_presence_threshold_m)
        (inputs.axis_primary_sink_strength_m_by_domain_layer_axis[controller_axis_layer] + inputs.axis_secondary_sink_strength_m_by_domain_layer_axis[controller_axis_layer]) / total_sink
    else
        1;
    if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidPrimaryRootWithdrawalFraction;

    // Exact NN loop and negative RCO2Z/... withdrawal-ledger convention.
    for (0..inputs.biological_domain_count) |domain| {
        for (0..inputs.gas_species_count) |species| {
            const index = speciesIndex(inputs, domain, layer, species);
            const removed_g = fraction * (model_state.gaseous_mass_g_by_domain_layer_species[index] + model_state.aqueous_mass_g_by_domain_layer_species[index]);
            model_state.withdrawal_ledger_g_by_species[species] -= removed_g;
            model_state.gaseous_mass_g_by_domain_layer_species[index] *= 1 - fraction;
            model_state.aqueous_mass_g_by_domain_layer_species[index] *= 1 - fraction;
        }
    }

    const count = model_state.secondary_axis_count_by_domain_layer_axis[controller_axis_layer];
    const upper_controller_layer = layerIndex(inputs, controlling_domain, layer - 1);
    model_state.total_secondary_axis_count_by_domain_layer[controller_layer] -= count;
    model_state.total_secondary_axis_count_by_domain_layer[upper_controller_layer] += count;
    model_state.secondary_axis_count_by_domain_layer_axis[controller_axis_layer] = 0;
    model_state.primary_axis_count_by_domain_layer[controller_layer] -= inputs.primary_axis_withdrawal_count_by_domain_layer_axis[controller_axis_layer];
    const upper_axis_layer = axisLayerIndex(inputs, controlling_domain, layer - 1, axis);
    const depth_m = inputs.retracted_depth_m_by_domain_axis[controlling_domain * inputs.root_axis_count + axis];
    var reconstructed_length_m = inputs.layer_thickness_m[layer - 1] - (inputs.layer_bottom_depth_m[layer - 1] - depth_m);
    if (layer - 1 == inputs.planting_layer_index) reconstructed_length_m -= inputs.seeding_depth_m - inputs.depth_above_planting_layer_bottom_m;
    if (!std.math.isFinite(reconstructed_length_m) or reconstructed_length_m < 0 or model_state.total_secondary_axis_count_by_domain_layer[controller_layer] < 0 or model_state.primary_axis_count_by_domain_layer[controller_layer] < 0) return error.PrimaryRootWithdrawalWouldOverdrawState;
    model_state.primary_length_m_by_domain_layer_axis[upper_axis_layer] = reconstructed_length_m;
    _ = trigger_layer;
}

/// Exact GROSUB 7191--7251 gas withdrawal and root-count/length reset.
/// Traversal is N, shared-NI L, runtime NR, descending LL, runtime gas species,
/// then NN. Gas masses use each species' elemental gram basis; length is m.
pub fn apply(model_state: State, workspace: State, inputs: Inputs) !void {
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.PrimaryRootWithdrawalDimensionOverflow;
    const axis_layers = std.math.mul(usize, domain_layers, inputs.root_axis_count) catch return error.PrimaryRootWithdrawalDimensionOverflow;
    const domain_axes = std.math.mul(usize, inputs.biological_domain_count, inputs.root_axis_count) catch return error.PrimaryRootWithdrawalDimensionOverflow;
    const species_values = std.math.mul(usize, domain_layers, inputs.gas_species_count) catch return error.PrimaryRootWithdrawalDimensionOverflow;
    try validateState(model_state, species_values, domain_layers, axis_layers, inputs.gas_species_count);
    try validateState(workspace, species_values, domain_layers, axis_layers, inputs.gas_species_count);
    try validateInputs(inputs, domain_layers, axis_layers, domain_axes);
    copyState(workspace, model_state);
    for (0..inputs.biological_domain_count) |domain| for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |trigger_layer| for (0..inputs.root_axis_count) |axis| {
        const trigger = axisLayerIndex(inputs, domain, trigger_layer, axis);
        if (!inputs.trigger_active_by_domain_layer_axis[trigger] or trigger_layer != inputs.deepest_rooted_layer_index_by_axis[axis]) continue;
        var layer = trigger_layer;
        while (layer > inputs.planting_layer_index) : (layer -= 1) {
            const depth_m = inputs.retracted_depth_m_by_domain_axis[domain * inputs.root_axis_count + axis];
            if (!(inputs.layer_thickness_m[layer - 1] > inputs.minimum_active_layer_thickness_m and (depth_m < inputs.layer_bottom_depth_m[layer - 1] or depth_m < inputs.seeding_depth_m))) break;
            try applyWithdrawalLayer(workspace, inputs, domain, trigger_layer, layer, axis);
        }
    };
    try validateState(workspace, species_values, domain_layers, axis_layers, inputs.gas_species_count);
    copyState(model_state, workspace);
}

test "GROSUB root withdrawal conserves runtime gas species with signed ledger and resets topology" {
    const domains = 2;
    const layers = 2;
    const axes = 1;
    const gases = 3;
    var gaseous = [_]f64{ 0, 0, 0, 2, 4, 6, 0, 0, 0, 1, 3, 5 };
    var aqueous = [_]f64{ 0, 0, 0, 1, 2, 3, 0, 0, 0, 1, 1, 1 };
    var ledger = [_]f64{ 0, 0, 0 };
    var total_secondary = [_]f64{ 0, 4, 0, 2 };
    var secondary = [_]f64{ 0, 2, 0, 1 };
    var primary = [_]f64{ 0, 2, 0, 1 };
    var length = [_]f64{ 0, 1, 0, 1 };
    var wg = [_]f64{0} ** (domains * layers * gases);
    var wa = wg;
    var wl = [_]f64{0} ** gases;
    var wts = [_]f64{0} ** (domains * layers);
    var ws = wts;
    var wp = wts;
    var wlen = [_]f64{0} ** (domains * layers * axes);
    const state: State = .{ .gaseous_mass_g_by_domain_layer_species = &gaseous, .aqueous_mass_g_by_domain_layer_species = &aqueous, .withdrawal_ledger_g_by_species = &ledger, .total_secondary_axis_count_by_domain_layer = &total_secondary, .secondary_axis_count_by_domain_layer_axis = &secondary, .primary_axis_count_by_domain_layer = &primary, .primary_length_m_by_domain_layer_axis = &length };
    const work: State = .{ .gaseous_mass_g_by_domain_layer_species = &wg, .aqueous_mass_g_by_domain_layer_species = &wa, .withdrawal_ledger_g_by_species = &wl, .total_secondary_axis_count_by_domain_layer = &wts, .secondary_axis_count_by_domain_layer_axis = &ws, .primary_axis_count_by_domain_layer = &wp, .primary_length_m_by_domain_layer_axis = &wlen };
    try apply(state, work, .{ .biological_domain_count = domains, .soil_layer_count = layers, .root_axis_count = axes, .gas_species_count = gases, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .trigger_active_by_domain_layer_axis = &.{ false, true, false, false }, .deepest_rooted_layer_index_by_axis = &.{1}, .retracted_depth_m_by_domain_axis = &.{ 0.1, 0.1 }, .seeding_depth_m = 0.15, .depth_above_planting_layer_bottom_m = 0.1, .layer_bottom_depth_m = &.{ 0.2, 0.4 }, .layer_thickness_m = &.{ 0.2, 0.2 }, .minimum_active_layer_thickness_m = 0.01, .axis_primary_sink_strength_m_by_domain_layer_axis = &.{ 0, 1, 0, 0 }, .axis_secondary_sink_strength_m_by_domain_layer_axis = &.{ 0, 1, 0, 0 }, .total_root_sink_strength_m_by_domain_layer = &.{ 0, 4, 0, 0 }, .sink_presence_threshold_m = 0, .primary_axis_withdrawal_count_by_domain_layer_axis = &.{ 0, 1, 0, 0 } });
    const initial = [_]f64{ 5, 10, 15 };
    for (0..gases) |species| {
        var remaining: f64 = 0;
        for (0..domains) |domain| {
            const index = (domain * layers + 1) * gases + species;
            remaining += gaseous[index] + aqueous[index];
        }
        try std.testing.expectApproxEqAbs(initial[species], remaining - ledger[species], 1e-12);
    }
    try std.testing.expectEqual(@as(f64, 0), secondary[1]);
    try std.testing.expectApproxEqAbs(0.05, length[0], 1e-12);
}

test "GROSUB root gas and topology withdrawal is atomic on count overdraw" {
    var gas = [_]f64{ 0, 1 };
    var aqueous = [_]f64{ 0, 1 };
    var ledger = [_]f64{0};
    var totals = [_]f64{ 0, 0.5 };
    var secondary = [_]f64{ 0, 1 };
    var primary = [_]f64{ 0, 1 };
    var length = [_]f64{ 0, 1 };
    var wg = [_]f64{0} ** 2;
    var wa = wg;
    var wl = [_]f64{0};
    var wt = [_]f64{0} ** 2;
    var ws = wt;
    var wp = wt;
    var wlen = wt;
    const before = totals;
    try std.testing.expectError(error.PrimaryRootWithdrawalWouldOverdrawState, apply(.{ .gaseous_mass_g_by_domain_layer_species = &gas, .aqueous_mass_g_by_domain_layer_species = &aqueous, .withdrawal_ledger_g_by_species = &ledger, .total_secondary_axis_count_by_domain_layer = &totals, .secondary_axis_count_by_domain_layer_axis = &secondary, .primary_axis_count_by_domain_layer = &primary, .primary_length_m_by_domain_layer_axis = &length }, .{ .gaseous_mass_g_by_domain_layer_species = &wg, .aqueous_mass_g_by_domain_layer_species = &wa, .withdrawal_ledger_g_by_species = &wl, .total_secondary_axis_count_by_domain_layer = &wt, .secondary_axis_count_by_domain_layer_axis = &ws, .primary_axis_count_by_domain_layer = &wp, .primary_length_m_by_domain_layer_axis = &wlen }, .{ .biological_domain_count = 1, .soil_layer_count = 2, .root_axis_count = 1, .gas_species_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .trigger_active_by_domain_layer_axis = &.{ false, true }, .deepest_rooted_layer_index_by_axis = &.{1}, .retracted_depth_m_by_domain_axis = &.{0}, .seeding_depth_m = 0.1, .depth_above_planting_layer_bottom_m = 0, .layer_bottom_depth_m = &.{ 0.2, 0.4 }, .layer_thickness_m = &.{ 0.2, 0.2 }, .minimum_active_layer_thickness_m = 0, .axis_primary_sink_strength_m_by_domain_layer_axis = &.{ 0, 1 }, .axis_secondary_sink_strength_m_by_domain_layer_axis = &.{ 0, 0 }, .total_root_sink_strength_m_by_domain_layer = &.{ 0, 1 }, .sink_presence_threshold_m = 0, .primary_axis_withdrawal_count_by_domain_layer_axis = &.{ 0, 0 } }));
    try std.testing.expectEqualDeep(before, totals);
}
