const std = @import("std");
const metabolism = @import("plant_root_metabolism.zig");

pub const State = struct {
    structural_carbon_g_c_by_domain_layer_axis: []f64,
    structural_nitrogen_g_n_by_domain_layer_axis: []f64,
    structural_phosphorus_g_p_by_domain_layer_axis: []f64,
    length_m_by_domain_layer_axis: []f64,
    mobile_carbon_g_c_by_domain_layer: []f64,
    mobile_nitrogen_g_n_by_domain_layer: []f64,
    mobile_phosphorus_g_p_by_domain_layer: []f64,
    litter_by_domain_layer: []metabolism.RootLitter,
};

pub const Workspace = State;

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    host_domain_index: usize,
    mycorrhizal_domain_index: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    tip_active_by_domain_layer_axis: []const bool,
    deficit_absorption_by_domain_layer_axis: []const metabolism.SecondaryRootDeficitAbsorption,
    host_active_root_carbon_g_c_by_domain_layer: []const f64,
    presence_threshold_g_c: f64,
    woody_fraction: [3][2]f64,
    kinetics: metabolism.RootLitterFractions,
};

fn copyState(target: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(target, field.name), @field(source, field.name));
}

fn validateShape(state: State, domain_layers: usize, values: usize) !void {
    inline for (.{ state.structural_carbon_g_c_by_domain_layer_axis, state.structural_nitrogen_g_n_by_domain_layer_axis, state.structural_phosphorus_g_p_by_domain_layer_axis, state.length_m_by_domain_layer_axis }) |slice|
        if (slice.len != values) return error.MycorrhizalLossSweepDimensionMismatch;
    inline for (.{ state.mobile_carbon_g_c_by_domain_layer, state.mobile_nitrogen_g_n_by_domain_layer, state.mobile_phosphorus_g_p_by_domain_layer }) |slice|
        if (slice.len != domain_layers) return error.MycorrhizalLossSweepDimensionMismatch;
    if (state.litter_by_domain_layer.len != domain_layers) return error.MycorrhizalLossSweepDimensionMismatch;
}

fn addLitter(total: *metabolism.RootLitter, addition: metabolism.RootLitter) !void {
    inline for (@typeInfo(metabolism.RootLitter).@"struct".fields) |field| for (&@field(total, field.name), @field(addition, field.name)) |*target, value| {
        target.* += value;
        if (!std.math.isFinite(target.*) or target.* < 0) return error.NonFiniteMycorrhizalLoss;
    };
}

fn validateDeficitPlan(plan: metabolism.SecondaryRootDeficitAbsorption) !void {
    inline for (.{ plan.current_entering_carbon_deficit_g_c, plan.upper_entering_carbon_deficit_g_c, plan.residual_carbon_deficit_g_c, plan.residual_nitrogen_deficit_g_n, plan.residual_phosphorus_deficit_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.NonFiniteMycorrhizalLossInput;
    inline for (.{ plan.current, plan.upper }) |layer| inline for (@typeInfo(metabolism.SecondaryRootDeficitLayer).@"struct".fields) |field| {
        const value = @field(layer, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.NonFiniteMycorrhizalLossInput;
    };
}

fn updateLayer(state: State, inputs: Inputs, axis: usize, layer: usize, entering_deficit_g_c: f64, host_remaining_g_c: f64) !void {
    const host_layer = inputs.host_domain_index * inputs.soil_layer_count + layer;
    const target_layer = inputs.mycorrhizal_domain_index * inputs.soil_layer_count + layer;
    const target = target_layer * inputs.root_axis_count + axis;
    const result = try metabolism.mycorrhizalLossWithSecondaryRoots(
        -entering_deficit_g_c,
        host_remaining_g_c,
        inputs.host_active_root_carbon_g_c_by_domain_layer[host_layer],
        inputs.presence_threshold_g_c,
        .{
            .structural_carbon_g_c = state.structural_carbon_g_c_by_domain_layer_axis[target],
            .structural_nitrogen_g_n = state.structural_nitrogen_g_n_by_domain_layer_axis[target],
            .structural_phosphorus_g_p = state.structural_phosphorus_g_p_by_domain_layer_axis[target],
            .length_m = state.length_m_by_domain_layer_axis[target],
            .mobile_carbon_g_c = state.mobile_carbon_g_c_by_domain_layer[target_layer],
            .mobile_nitrogen_g_n = state.mobile_nitrogen_g_n_by_domain_layer[target_layer],
            .mobile_phosphorus_g_p = state.mobile_phosphorus_g_p_by_domain_layer[target_layer],
        },
        inputs.woody_fraction,
        inputs.kinetics,
    );
    state.structural_carbon_g_c_by_domain_layer_axis[target] = result.remaining.structural_carbon_g_c;
    state.structural_nitrogen_g_n_by_domain_layer_axis[target] = result.remaining.structural_nitrogen_g_n;
    state.structural_phosphorus_g_p_by_domain_layer_axis[target] = result.remaining.structural_phosphorus_g_p;
    state.length_m_by_domain_layer_axis[target] = result.remaining.length_m;
    state.mobile_carbon_g_c_by_domain_layer[target_layer] = result.remaining.mobile_carbon_g_c;
    state.mobile_nitrogen_g_n_by_domain_layer[target_layer] = result.remaining.mobile_nitrogen_g_n;
    state.mobile_phosphorus_g_p_by_domain_layer[target_layer] = result.remaining.mobile_phosphorus_g_p;
    try addLitter(&state.litter_by_domain_layer[target_layer], result.litter);
}

/// Exact GROSUB lines 6907--6951 N,L,NR then current/upper LL concurrent
/// mycorrhizal loss. Masses are g C/N/P, length is m, and litter retains the
/// source four kinetic pools.
pub fn apply(state: State, workspace: Workspace, inputs: Inputs) !void {
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.MycorrhizalLossSweepDimensionOverflow;
    const values = std.math.mul(usize, domain_layers, inputs.root_axis_count) catch return error.MycorrhizalLossSweepDimensionOverflow;
    try validateShape(state, domain_layers, values);
    try validateShape(workspace, domain_layers, values);
    if (inputs.biological_domain_count < 2 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or inputs.host_domain_index >= inputs.biological_domain_count or inputs.mycorrhizal_domain_index >= inputs.biological_domain_count or inputs.host_domain_index == inputs.mycorrhizal_domain_index or
        inputs.tip_active_by_domain_layer_axis.len != values or inputs.deficit_absorption_by_domain_layer_axis.len != values or inputs.host_active_root_carbon_g_c_by_domain_layer.len != domain_layers)
        return error.MycorrhizalLossSweepDimensionMismatch;
    if (!std.math.isFinite(inputs.presence_threshold_g_c) or inputs.presence_threshold_g_c < 0) return error.InvalidMycorrhizalLossInput;
    for (inputs.deficit_absorption_by_domain_layer_axis) |plan| try validateDeficitPlan(plan);
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count) return error.InvalidMycorrhizalLossLayerRange;
    copyState(workspace, state);
    const host = inputs.host_domain_index;
    for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| for (0..inputs.root_axis_count) |axis| {
        const index = (host * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
        if (!inputs.tip_active_by_domain_layer_axis[index]) continue;
        const deficit_plan = inputs.deficit_absorption_by_domain_layer_axis[index];
        if (deficit_plan.current_entering_carbon_deficit_g_c > 0) try updateLayer(workspace, inputs, axis, layer, deficit_plan.current_entering_carbon_deficit_g_c, deficit_plan.current.carbon_g_c);
        if (layer > inputs.planting_layer_index and deficit_plan.upper_entering_carbon_deficit_g_c > 0) try updateLayer(workspace, inputs, axis, layer - 1, deficit_plan.upper_entering_carbon_deficit_g_c, deficit_plan.upper.carbon_g_c);
    };
    copyState(state, workspace);
}

fn zeroLitter() metabolism.RootLitter {
    return std.mem.zeroes(metabolism.RootLitter);
}

fn absorption(current_deficit: f64, upper_deficit: f64) metabolism.SecondaryRootDeficitAbsorption {
    return .{ .current = .{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01, .length_m = 1 }, .upper = .{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01, .length_m = 1 }, .current_entering_carbon_deficit_g_c = current_deficit, .upper_entering_carbon_deficit_g_c = upper_deficit, .residual_carbon_deficit_g_c = 0, .residual_nitrogen_deficit_g_n = 0, .residual_phosphorus_deficit_g_p = 0 };
}

test "GROSUB mycorrhizal current upper losses conserve C N P into litter" {
    const domains = 2;
    const layers = 2;
    const axes = 2;
    const values = domains * layers * axes;
    var c = [_]f64{2} ** values;
    var n = [_]f64{0.2} ** values;
    var p = [_]f64{0.02} ** values;
    var length = [_]f64{1} ** values;
    var mc = [_]f64{4} ** (domains * layers);
    var mn = [_]f64{0.4} ** (domains * layers);
    var mp = [_]f64{0.04} ** (domains * layers);
    var litter = [_]metabolism.RootLitter{zeroLitter()} ** (domains * layers);
    var wc = [_]f64{0} ** values;
    var wn = [_]f64{0} ** values;
    var wp = [_]f64{0} ** values;
    var wl = [_]f64{0} ** values;
    var wmc = [_]f64{0} ** (domains * layers);
    var wmn = [_]f64{0} ** (domains * layers);
    var wmp = [_]f64{0} ** (domains * layers);
    var wlit = [_]metabolism.RootLitter{zeroLitter()} ** (domains * layers);
    var plans = [_]metabolism.SecondaryRootDeficitAbsorption{absorption(0, 0)} ** values;
    plans[2] = absorption(1, 0.5);
    const target_before_c = c[4] + c[6] + mc[2] + mc[3];
    const kinetics = [4]f64{ 0.1, 0.2, 0.3, 0.4 };
    try apply(.{ .structural_carbon_g_c_by_domain_layer_axis = &c, .structural_nitrogen_g_n_by_domain_layer_axis = &n, .structural_phosphorus_g_p_by_domain_layer_axis = &p, .length_m_by_domain_layer_axis = &length, .mobile_carbon_g_c_by_domain_layer = &mc, .mobile_nitrogen_g_n_by_domain_layer = &mn, .mobile_phosphorus_g_p_by_domain_layer = &mp, .litter_by_domain_layer = &litter }, .{ .structural_carbon_g_c_by_domain_layer_axis = &wc, .structural_nitrogen_g_n_by_domain_layer_axis = &wn, .structural_phosphorus_g_p_by_domain_layer_axis = &wp, .length_m_by_domain_layer_axis = &wl, .mobile_carbon_g_c_by_domain_layer = &wmc, .mobile_nitrogen_g_n_by_domain_layer = &wmn, .mobile_phosphorus_g_p_by_domain_layer = &wmp, .litter_by_domain_layer = &wlit }, .{
        .biological_domain_count = domains,
        .soil_layer_count = layers,
        .root_axis_count = axes,
        .host_domain_index = 0,
        .mycorrhizal_domain_index = 1,
        .planting_layer_index = 0,
        .deepest_rooted_layer_index = 1,
        .tip_active_by_domain_layer_axis = &.{ false, false, true, false, false, false, false, false },
        .deficit_absorption_by_domain_layer_axis = &plans,
        .host_active_root_carbon_g_c_by_domain_layer = &.{ 3, 3, 0, 0 },
        .presence_threshold_g_c = 1e-12,
        .woody_fraction = .{ .{ 0.25, 0.75 }, .{ 0.25, 0.75 }, .{ 0.25, 0.75 } },
        .kinetics = .{ .woody_carbon = kinetics, .woody_nitrogen = kinetics, .woody_phosphorus = kinetics, .nonwoody_carbon = kinetics, .nonwoody_nitrogen = kinetics, .nonwoody_phosphorus = kinetics },
    });
    var litter_c: f64 = 0;
    for (litter[2].woody_carbon_g_c) |value| litter_c += value;
    for (litter[2].nonwoody_carbon_g_c) |value| litter_c += value;
    for (litter[3].woody_carbon_g_c) |value| litter_c += value;
    for (litter[3].nonwoody_carbon_g_c) |value| litter_c += value;
    try std.testing.expectApproxEqAbs(target_before_c, c[4] + c[6] + mc[2] + mc[3] + litter_c, 1e-12);
}

test "GROSUB invalid late mycorrhizal loss rolls back earlier axis" {
    const values = 4;
    var c = [_]f64{2} ** values;
    var n = [_]f64{0.2} ** values;
    var p = [_]f64{0.02} ** values;
    var length = [_]f64{1} ** values;
    var mc = [_]f64{4} ** 2;
    var mn = [_]f64{0.4} ** 2;
    var mp = [_]f64{0.04} ** 2;
    var litter = [_]metabolism.RootLitter{zeroLitter()} ** 2;
    var workspace_values: [4][4]f64 = std.mem.zeroes([4][4]f64);
    var workspace_layers: [3][2]f64 = std.mem.zeroes([3][2]f64);
    var workspace_litter = [_]metabolism.RootLitter{zeroLitter()} ** 2;
    var plans = [_]metabolism.SecondaryRootDeficitAbsorption{absorption(1, 0)} ** values;
    plans[1].current_entering_carbon_deficit_g_c = std.math.nan(f64);
    const before_c = c;
    const kinetics = [4]f64{ 0.1, 0.2, 0.3, 0.4 };
    try std.testing.expectError(error.NonFiniteMycorrhizalLossInput, apply(.{ .structural_carbon_g_c_by_domain_layer_axis = &c, .structural_nitrogen_g_n_by_domain_layer_axis = &n, .structural_phosphorus_g_p_by_domain_layer_axis = &p, .length_m_by_domain_layer_axis = &length, .mobile_carbon_g_c_by_domain_layer = &mc, .mobile_nitrogen_g_n_by_domain_layer = &mn, .mobile_phosphorus_g_p_by_domain_layer = &mp, .litter_by_domain_layer = &litter }, .{ .structural_carbon_g_c_by_domain_layer_axis = &workspace_values[0], .structural_nitrogen_g_n_by_domain_layer_axis = &workspace_values[1], .structural_phosphorus_g_p_by_domain_layer_axis = &workspace_values[2], .length_m_by_domain_layer_axis = &workspace_values[3], .mobile_carbon_g_c_by_domain_layer = &workspace_layers[0], .mobile_nitrogen_g_n_by_domain_layer = &workspace_layers[1], .mobile_phosphorus_g_p_by_domain_layer = &workspace_layers[2], .litter_by_domain_layer = &workspace_litter }, .{ .biological_domain_count = 2, .soil_layer_count = 1, .root_axis_count = 2, .host_domain_index = 0, .mycorrhizal_domain_index = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 0, .tip_active_by_domain_layer_axis = &.{ true, true, false, false }, .deficit_absorption_by_domain_layer_axis = &plans, .host_active_root_carbon_g_c_by_domain_layer = &.{ 3, 0 }, .presence_threshold_g_c = 1e-12, .woody_fraction = .{ .{ 0.25, 0.75 }, .{ 0.25, 0.75 }, .{ 0.25, 0.75 } }, .kinetics = .{ .woody_carbon = kinetics, .woody_nitrogen = kinetics, .woody_phosphorus = kinetics, .nonwoody_carbon = kinetics, .nonwoody_nitrogen = kinetics, .nonwoody_phosphorus = kinetics } }));
    try std.testing.expectEqualDeep(before_c, c);
}
