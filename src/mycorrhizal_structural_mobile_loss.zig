const std = @import("std");

pub const State = struct {
    structural_carbon_g_c: []f64,
    structural_nitrogen_g_n: []f64,
    structural_phosphorus_g_p: []f64,
    length_m: []f64,
    mobile_carbon_g_c: []f64,
    mobile_nitrogen_g_n: []f64,
    mobile_phosphorus_g_p: []f64,
    litter_carbon_g_c: []f64,
    litter_nitrogen_g_n: []f64,
    litter_phosphorus_g_p: []f64,
};

pub const KineticFractions = struct {
    woody_carbon: []const f64,
    woody_nitrogen: []const f64,
    woody_phosphorus: []const f64,
    nonwoody_carbon: []const f64,
    nonwoody_nitrogen: []const f64,
    nonwoody_phosphorus: []const f64,
};

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    kinetic_pool_count: usize,
    host_domain_index: usize,
    mycorrhizal_domain_index: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    tip_active_by_layer_axis: []const bool,
    current_entering_deficit_g_c_by_layer_axis: []const f64,
    upper_entering_deficit_g_c_by_layer_axis: []const f64,
    current_host_secondary_carbon_g_c_by_layer_axis: []const f64,
    upper_host_secondary_carbon_g_c_by_layer_axis: []const f64,
    host_active_root_carbon_g_c_by_layer: []const f64,
    presence_threshold_g_c: f64,
    woody_fraction_by_element: [3][2]f64,
    kinetics: KineticFractions,
};

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        @memcpy(@field(destination, field.name), @field(source, field.name));
    }
}

fn validateState(state: State, axis_values: usize, layer_values: usize, litter_values: usize) !void {
    inline for (.{ state.structural_carbon_g_c, state.structural_nitrogen_g_n, state.structural_phosphorus_g_p, state.length_m }) |values| {
        if (values.len != axis_values) return error.MycorrhizalLossDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidMycorrhizalLossState;
    }
    inline for (.{ state.mobile_carbon_g_c, state.mobile_nitrogen_g_n, state.mobile_phosphorus_g_p }) |values| {
        if (values.len != layer_values) return error.MycorrhizalLossDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidMycorrhizalLossState;
    }
    inline for (.{ state.litter_carbon_g_c, state.litter_nitrogen_g_n, state.litter_phosphorus_g_p }) |values| {
        if (values.len != litter_values) return error.MycorrhizalLossDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidMycorrhizalLossState;
    }
}

fn validateInputs(inputs: Inputs, host_axis_values: usize) !void {
    if (inputs.biological_domain_count < 2 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or inputs.kinetic_pool_count == 0 or
        inputs.host_domain_index >= inputs.biological_domain_count or inputs.mycorrhizal_domain_index >= inputs.biological_domain_count or
        inputs.host_domain_index == inputs.mycorrhizal_domain_index or inputs.planting_layer_index > inputs.deepest_rooted_layer_index or
        inputs.deepest_rooted_layer_index >= inputs.soil_layer_count or inputs.tip_active_by_layer_axis.len != host_axis_values or
        inputs.current_entering_deficit_g_c_by_layer_axis.len != host_axis_values or inputs.upper_entering_deficit_g_c_by_layer_axis.len != host_axis_values or
        inputs.current_host_secondary_carbon_g_c_by_layer_axis.len != host_axis_values or inputs.upper_host_secondary_carbon_g_c_by_layer_axis.len != host_axis_values or
        inputs.host_active_root_carbon_g_c_by_layer.len != inputs.soil_layer_count)
        return error.MycorrhizalLossDimensionMismatch;
    if (!std.math.isFinite(inputs.presence_threshold_g_c) or inputs.presence_threshold_g_c < 0) return error.InvalidMycorrhizalLossInput;
    inline for (@typeInfo(KineticFractions).@"struct".fields) |field| {
        const values = @field(inputs.kinetics, field.name);
        if (values.len != inputs.kinetic_pool_count) return error.MycorrhizalLossDimensionMismatch;
        var total: f64 = 0;
        for (values) |value| {
            if (!std.math.isFinite(value) or value < 0) return error.InvalidMycorrhizalLossInput;
            total += value;
        }
        if (!std.math.isFinite(total) or @abs(total - 1) > 1e-12) return error.InvalidMycorrhizalKineticFractions;
    }
    for (inputs.woody_fraction_by_element) |fractions| {
        for (fractions) |value| if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidMycorrhizalLossInput;
        if (@abs(fractions[0] + fractions[1] - 1) > 1e-12) return error.InvalidMycorrhizalWoodFractions;
    }
    inline for (.{ inputs.current_entering_deficit_g_c_by_layer_axis, inputs.upper_entering_deficit_g_c_by_layer_axis, inputs.current_host_secondary_carbon_g_c_by_layer_axis, inputs.upper_host_secondary_carbon_g_c_by_layer_axis }) |values| {
        for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidMycorrhizalLossInput;
    }
    for (inputs.host_active_root_carbon_g_c_by_layer) |carbon| if (!std.math.isFinite(carbon) or carbon < 0) return error.InvalidMycorrhizalLossInput;
}

fn litterIndex(inputs: Inputs, domain: usize, layer: usize, category: usize, pool: usize) usize {
    return (((domain * inputs.soil_layer_count + layer) * 2 + category) * inputs.kinetic_pool_count) + pool;
}

fn applyLayerAxis(state: State, inputs: Inputs, layer: usize, axis: usize, deficit_g_c: f64, host_secondary_carbon_g_c: f64) !void {
    if (deficit_g_c <= 0) return;
    const target_layer = inputs.mycorrhizal_domain_index * inputs.soil_layer_count + layer;
    const target_axis = target_layer * inputs.root_axis_count + axis;
    const structural_fraction = if (host_secondary_carbon_g_c > inputs.presence_threshold_g_c)
        @min(1, deficit_g_c / host_secondary_carbon_g_c)
    else
        1;
    const mobile_fraction = if (inputs.host_active_root_carbon_g_c_by_layer[layer] > inputs.presence_threshold_g_c)
        @min(1, deficit_g_c / inputs.host_active_root_carbon_g_c_by_layer[layer])
    else
        1;

    const structure = [3]f64{ state.structural_carbon_g_c[target_axis], state.structural_nitrogen_g_n[target_axis], state.structural_phosphorus_g_p[target_axis] };
    const mobile = [3]f64{ state.mobile_carbon_g_c[target_layer], state.mobile_nitrogen_g_n[target_layer], state.mobile_phosphorus_g_p[target_layer] };
    const kinetic_sets = [3][2][]const f64{
        .{ inputs.kinetics.woody_carbon, inputs.kinetics.nonwoody_carbon },
        .{ inputs.kinetics.woody_nitrogen, inputs.kinetics.nonwoody_nitrogen },
        .{ inputs.kinetics.woody_phosphorus, inputs.kinetics.nonwoody_phosphorus },
    };
    const litter_sets = [3][]f64{ state.litter_carbon_g_c, state.litter_nitrogen_g_n, state.litter_phosphorus_g_p };
    for (0..inputs.kinetic_pool_count) |pool| {
        for (0..3) |element| {
            const woody = kinetic_sets[element][0][pool] * structural_fraction * structure[element] * inputs.woody_fraction_by_element[element][0];
            const nonwoody = kinetic_sets[element][1][pool] * (structural_fraction * structure[element] * inputs.woody_fraction_by_element[element][1] + mobile_fraction * mobile[element]);
            const woody_index = litterIndex(inputs, inputs.mycorrhizal_domain_index, layer, 0, pool);
            const nonwoody_index = litterIndex(inputs, inputs.mycorrhizal_domain_index, layer, 1, pool);
            litter_sets[element][woody_index] += woody;
            litter_sets[element][nonwoody_index] += nonwoody;
            if (!std.math.isFinite(litter_sets[element][woody_index]) or !std.math.isFinite(litter_sets[element][nonwoody_index])) return error.NonFiniteMycorrhizalLoss;
        }
    }
    const structural_retained = 1 - structural_fraction;
    const mobile_retained = 1 - mobile_fraction;
    state.structural_carbon_g_c[target_axis] *= structural_retained;
    state.structural_nitrogen_g_n[target_axis] *= structural_retained;
    state.structural_phosphorus_g_p[target_axis] *= structural_retained;
    state.length_m[target_axis] *= structural_retained;
    state.mobile_carbon_g_c[target_layer] *= mobile_retained;
    state.mobile_nitrogen_g_n[target_layer] *= mobile_retained;
    state.mobile_phosphorus_g_p[target_layer] *= mobile_retained;
}

/// Exact GROSUB 6907--6951 concurrent mycorrhizal structural/mobile loss.
/// Source order is layer (L), root axis (NR), current then upper layer (LL),
/// then kinetic pool (M); the caller invokes this once for each plant (N).
/// Masses are g C, g N, and g P; length is m.
/// `workspace` makes the complete publication atomic.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    const layer_values = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.MycorrhizalLossDimensionOverflow;
    const axis_values = std.math.mul(usize, layer_values, inputs.root_axis_count) catch return error.MycorrhizalLossDimensionOverflow;
    const litter_values = std.math.mul(usize, std.math.mul(usize, layer_values, 2) catch return error.MycorrhizalLossDimensionOverflow, inputs.kinetic_pool_count) catch return error.MycorrhizalLossDimensionOverflow;
    const host_axis_values = std.math.mul(usize, inputs.soil_layer_count, inputs.root_axis_count) catch return error.MycorrhizalLossDimensionOverflow;
    try validateState(state, axis_values, layer_values, litter_values);
    try validateState(workspace, axis_values, layer_values, litter_values);
    try validateInputs(inputs, host_axis_values);
    copyState(workspace, state);
    for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| {
        for (0..inputs.root_axis_count) |axis| {
            const host_axis = layer * inputs.root_axis_count + axis;
            if (!inputs.tip_active_by_layer_axis[host_axis]) continue;
            try applyLayerAxis(workspace, inputs, layer, axis, inputs.current_entering_deficit_g_c_by_layer_axis[host_axis], inputs.current_host_secondary_carbon_g_c_by_layer_axis[host_axis]);
            // Source `LX=MAX(1,L-1)` reaches the physical upper layer even
            // when the planting/rooted-loop lower bound is deeper than it.
            if (layer > 0) try applyLayerAxis(workspace, inputs, layer - 1, axis, inputs.upper_entering_deficit_g_c_by_layer_axis[host_axis], inputs.upper_host_secondary_carbon_g_c_by_layer_axis[host_axis]);
        }
    }
    copyState(state, workspace);
}

fn sum(values: []const f64) f64 {
    var total: f64 = 0;
    for (values) |value| total += value;
    return total;
}

test "GROSUB dynamic kinetic pools conserve concurrent mycorrhizal C N P" {
    const domains = 2;
    const layers = 1;
    const axes = 1;
    const pools = 3;
    var values = [4][2]f64{ .{ 0, 8 }, .{ 0, 0.8 }, .{ 0, 0.08 }, .{ 0, 4 } };
    var mobile = [3][2]f64{ .{ 0, 5 }, .{ 0, 0.5 }, .{ 0, 0.05 } };
    var litter = [3][domains * layers * 2 * pools]f64{ [_]f64{0} ** (domains * layers * 2 * pools), [_]f64{0} ** (domains * layers * 2 * pools), [_]f64{0} ** (domains * layers * 2 * pools) };
    var work_values: [4][2]f64 = std.mem.zeroes([4][2]f64);
    var work_mobile: [3][2]f64 = std.mem.zeroes([3][2]f64);
    var work_litter: [3][domains * layers * 2 * pools]f64 = std.mem.zeroes([3][domains * layers * 2 * pools]f64);
    const fractions = [_]f64{ 0.2, 0.3, 0.5 };
    const state = State{ .structural_carbon_g_c = &values[0], .structural_nitrogen_g_n = &values[1], .structural_phosphorus_g_p = &values[2], .length_m = &values[3], .mobile_carbon_g_c = &mobile[0], .mobile_nitrogen_g_n = &mobile[1], .mobile_phosphorus_g_p = &mobile[2], .litter_carbon_g_c = &litter[0], .litter_nitrogen_g_n = &litter[1], .litter_phosphorus_g_p = &litter[2] };
    const workspace = State{ .structural_carbon_g_c = &work_values[0], .structural_nitrogen_g_n = &work_values[1], .structural_phosphorus_g_p = &work_values[2], .length_m = &work_values[3], .mobile_carbon_g_c = &work_mobile[0], .mobile_nitrogen_g_n = &work_mobile[1], .mobile_phosphorus_g_p = &work_mobile[2], .litter_carbon_g_c = &work_litter[0], .litter_nitrogen_g_n = &work_litter[1], .litter_phosphorus_g_p = &work_litter[2] };
    try apply(state, workspace, .{ .biological_domain_count = domains, .soil_layer_count = layers, .root_axis_count = axes, .kinetic_pool_count = pools, .host_domain_index = 0, .mycorrhizal_domain_index = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 0, .tip_active_by_layer_axis = &.{true}, .current_entering_deficit_g_c_by_layer_axis = &.{2}, .upper_entering_deficit_g_c_by_layer_axis = &.{0}, .current_host_secondary_carbon_g_c_by_layer_axis = &.{8}, .upper_host_secondary_carbon_g_c_by_layer_axis = &.{0}, .host_active_root_carbon_g_c_by_layer = &.{10}, .presence_threshold_g_c = 1e-12, .woody_fraction_by_element = .{ .{ 0.25, 0.75 }, .{ 0.25, 0.75 }, .{ 0.25, 0.75 } }, .kinetics = .{ .woody_carbon = &fractions, .woody_nitrogen = &fractions, .woody_phosphorus = &fractions, .nonwoody_carbon = &fractions, .nonwoody_nitrogen = &fractions, .nonwoody_phosphorus = &fractions } });
    try std.testing.expectApproxEqAbs(13, values[0][1] + mobile[0][1] + sum(&litter[0]), 1e-12);
    try std.testing.expectApproxEqAbs(1.3, values[1][1] + mobile[1][1] + sum(&litter[1]), 1e-12);
    try std.testing.expectApproxEqAbs(0.13, values[2][1] + mobile[2][1] + sum(&litter[2]), 1e-12);
}

test "GROSUB dynamic mycorrhizal loss rejects late invalid state atomically" {
    var values = [4][2]f64{ .{ 0, 8 }, .{ 0, 0.8 }, .{ 0, 0.08 }, .{ 0, 4 } };
    var mobile = [3][2]f64{ .{ 0, 5 }, .{ 0, 0.5 }, .{ 0, 0.05 } };
    var litter: [3][4]f64 = std.mem.zeroes([3][4]f64);
    litter[2][3] = std.math.nan(f64);
    var work_values: [4][2]f64 = std.mem.zeroes([4][2]f64);
    var work_mobile: [3][2]f64 = std.mem.zeroes([3][2]f64);
    var work_litter: [3][4]f64 = std.mem.zeroes([3][4]f64);
    const before = values;
    const fractions = [_]f64{1};
    const state = State{ .structural_carbon_g_c = &values[0], .structural_nitrogen_g_n = &values[1], .structural_phosphorus_g_p = &values[2], .length_m = &values[3], .mobile_carbon_g_c = &mobile[0], .mobile_nitrogen_g_n = &mobile[1], .mobile_phosphorus_g_p = &mobile[2], .litter_carbon_g_c = &litter[0], .litter_nitrogen_g_n = &litter[1], .litter_phosphorus_g_p = &litter[2] };
    const workspace = State{ .structural_carbon_g_c = &work_values[0], .structural_nitrogen_g_n = &work_values[1], .structural_phosphorus_g_p = &work_values[2], .length_m = &work_values[3], .mobile_carbon_g_c = &work_mobile[0], .mobile_nitrogen_g_n = &work_mobile[1], .mobile_phosphorus_g_p = &work_mobile[2], .litter_carbon_g_c = &work_litter[0], .litter_nitrogen_g_n = &work_litter[1], .litter_phosphorus_g_p = &work_litter[2] };
    try std.testing.expectError(error.InvalidMycorrhizalLossState, apply(state, workspace, .{ .biological_domain_count = 2, .soil_layer_count = 1, .root_axis_count = 1, .kinetic_pool_count = 1, .host_domain_index = 0, .mycorrhizal_domain_index = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 0, .tip_active_by_layer_axis = &.{true}, .current_entering_deficit_g_c_by_layer_axis = &.{1}, .upper_entering_deficit_g_c_by_layer_axis = &.{0}, .current_host_secondary_carbon_g_c_by_layer_axis = &.{8}, .upper_host_secondary_carbon_g_c_by_layer_axis = &.{0}, .host_active_root_carbon_g_c_by_layer = &.{10}, .presence_threshold_g_c = 0, .woody_fraction_by_element = .{ .{ 0.5, 0.5 }, .{ 0.5, 0.5 }, .{ 0.5, 0.5 } }, .kinetics = .{ .woody_carbon = &fractions, .woody_nitrogen = &fractions, .woody_phosphorus = &fractions, .nonwoody_carbon = &fractions, .nonwoody_nitrogen = &fractions, .nonwoody_phosphorus = &fractions } }));
    try std.testing.expectEqualDeep(before, values);
}
