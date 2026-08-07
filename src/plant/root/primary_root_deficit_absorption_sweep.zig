const std = @import("std");
const metabolism = @import("plant_root_metabolism.zig");

pub const State = struct {
    secondary_carbon_g_c_by_domain_layer_axis: []f64,
    secondary_nitrogen_g_n_by_domain_layer_axis: []f64,
    secondary_phosphorus_g_p_by_domain_layer_axis: []f64,
    secondary_length_m_by_domain_layer_axis: []f64,
    /// Residual signed primary growth after host-secondary absorption.
    primary_growth_carbon_g_c_by_domain_layer_axis: []f64,
    primary_growth_nitrogen_g_n_by_domain_layer_axis: []f64,
    primary_growth_phosphorus_g_p_by_domain_layer_axis: []f64,
};

pub const Workspace = State;

pub const Inputs = struct {
    biological_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_index: usize,
    deepest_rooted_layer_index: usize,
    tip_active_by_domain_layer_axis: []const bool,
};

fn copyState(target: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(target, field.name), @field(source, field.name));
}

fn validateShape(state: State, values: usize) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (@field(state, field.name).len != values) return error.PrimaryRootDeficitSweepDimensionMismatch;
}

fn absorb(state: State, inputs: Inputs, domain: usize, layer: usize, axis: usize) !void {
    const index = (domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
    const initial_c = state.primary_growth_carbon_g_c_by_domain_layer_axis[index];
    const initial_n = state.primary_growth_nitrogen_g_n_by_domain_layer_axis[index];
    const initial_p = state.primary_growth_phosphorus_g_p_by_domain_layer_axis[index];
    inline for (.{ initial_c, initial_n, initial_p }) |value| if (!std.math.isFinite(value)) return error.NonFinitePrimaryRootGrowth;
    if (!(initial_c < 0)) return;
    // Source `LX=MAX(1,L-1)` is bounded by the first soil layer, not NU.
    const upper_layer = if (layer > 0) layer - 1 else layer;
    const current_index = index;
    const upper_index = (domain * inputs.soil_layer_count + upper_layer) * inputs.root_axis_count + axis;
    const result = try metabolism.absorbPrimaryDeficitFromSecondaryRoots(
        -initial_c,
        @max(0, -initial_n),
        @max(0, -initial_p),
        .{ .carbon_g_c = state.secondary_carbon_g_c_by_domain_layer_axis[current_index], .nitrogen_g_n = state.secondary_nitrogen_g_n_by_domain_layer_axis[current_index], .phosphorus_g_p = state.secondary_phosphorus_g_p_by_domain_layer_axis[current_index], .length_m = state.secondary_length_m_by_domain_layer_axis[current_index] },
        if (upper_index == current_index) .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0, .length_m = 0 } else .{ .carbon_g_c = state.secondary_carbon_g_c_by_domain_layer_axis[upper_index], .nitrogen_g_n = state.secondary_nitrogen_g_n_by_domain_layer_axis[upper_index], .phosphorus_g_p = state.secondary_phosphorus_g_p_by_domain_layer_axis[upper_index], .length_m = state.secondary_length_m_by_domain_layer_axis[upper_index] },
    );
    state.secondary_carbon_g_c_by_domain_layer_axis[current_index] = result.current.carbon_g_c;
    state.secondary_nitrogen_g_n_by_domain_layer_axis[current_index] = result.current.nitrogen_g_n;
    state.secondary_phosphorus_g_p_by_domain_layer_axis[current_index] = result.current.phosphorus_g_p;
    state.secondary_length_m_by_domain_layer_axis[current_index] = result.current.length_m;
    if (upper_index != current_index) {
        state.secondary_carbon_g_c_by_domain_layer_axis[upper_index] = result.upper.carbon_g_c;
        state.secondary_nitrogen_g_n_by_domain_layer_axis[upper_index] = result.upper.nitrogen_g_n;
        state.secondary_phosphorus_g_p_by_domain_layer_axis[upper_index] = result.upper.phosphorus_g_p;
        state.secondary_length_m_by_domain_layer_axis[upper_index] = result.upper.length_m;
    }
    state.primary_growth_carbon_g_c_by_domain_layer_axis[index] = -result.residual_carbon_deficit_g_c;
    if (initial_n < 0) state.primary_growth_nitrogen_g_n_by_domain_layer_axis[index] = -result.residual_nitrogen_deficit_g_n;
    if (initial_p < 0) state.primary_growth_phosphorus_g_p_by_domain_layer_axis[index] = -result.residual_phosphorus_deficit_g_p;
}

/// Exact grosub.f lines 6843--6879 N,L,NR then descending LL host-secondary
/// deficit absorption. C/N/P masses are grams and length is m.
pub fn apply(state: State, workspace: Workspace, inputs: Inputs) !void {
    const domain_layers = std.math.mul(usize, inputs.biological_domain_count, inputs.soil_layer_count) catch return error.PrimaryRootDeficitSweepDimensionOverflow;
    const values = std.math.mul(usize, domain_layers, inputs.root_axis_count) catch return error.PrimaryRootDeficitSweepDimensionOverflow;
    try validateShape(state, values);
    try validateShape(workspace, values);
    if (inputs.biological_domain_count == 0 or inputs.soil_layer_count == 0 or inputs.root_axis_count == 0 or inputs.tip_active_by_domain_layer_axis.len != values) return error.PrimaryRootDeficitSweepDimensionMismatch;
    if (inputs.planting_layer_index > inputs.deepest_rooted_layer_index or inputs.deepest_rooted_layer_index >= inputs.soil_layer_count) return error.InvalidPrimaryRootDeficitSweepLayerRange;
    for (1..inputs.biological_domain_count) |domain| for (0..inputs.soil_layer_count) |layer| for (0..inputs.root_axis_count) |axis| if (inputs.tip_active_by_domain_layer_axis[(domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis]) return error.PrimaryRootTipOutsideHostDomain;
    inline for (.{ state.secondary_carbon_g_c_by_domain_layer_axis, state.secondary_nitrogen_g_n_by_domain_layer_axis, state.secondary_phosphorus_g_p_by_domain_layer_axis, state.secondary_length_m_by_domain_layer_axis }) |slice| for (slice) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootDeficitPool;
    copyState(workspace, state);
    for (0..inputs.biological_domain_count) |domain| for (inputs.planting_layer_index..inputs.deepest_rooted_layer_index + 1) |layer| for (0..inputs.root_axis_count) |axis| {
        const index = (domain * inputs.soil_layer_count + layer) * inputs.root_axis_count + axis;
        if (inputs.tip_active_by_domain_layer_axis[index]) try absorb(workspace, inputs, domain, layer, axis);
    };
    copyState(state, workspace);
}

fn makeState(c: []f64, n: []f64, p: []f64, length: []f64, growth_c: []f64, growth_n: []f64, growth_p: []f64) State {
    return .{ .secondary_carbon_g_c_by_domain_layer_axis = c, .secondary_nitrogen_g_n_by_domain_layer_axis = n, .secondary_phosphorus_g_p_by_domain_layer_axis = p, .secondary_length_m_by_domain_layer_axis = length, .primary_growth_carbon_g_c_by_domain_layer_axis = growth_c, .primary_growth_nitrogen_g_n_by_domain_layer_axis = growth_n, .primary_growth_phosphorus_g_p_by_domain_layer_axis = growth_p };
}

test "GROSUB primary deficits consume current then upper secondary roots conservatively" {
    var c = [_]f64{ 2, 1 };
    var n = [_]f64{ 0.2, 0.1 };
    var p = [_]f64{ 0.02, 0.01 };
    var length = [_]f64{ 4, 2 };
    var growth_c = [_]f64{ 0, -2.5 };
    var growth_n = [_]f64{ 0, -0.25 };
    var growth_p = [_]f64{ 0, -0.025 };
    var work: [7][2]f64 = std.mem.zeroes([7][2]f64);
    const initial_c = c[0] + c[1] + growth_c[1];
    const initial_n = n[0] + n[1] + growth_n[1];
    const initial_p = p[0] + p[1] + growth_p[1];
    try apply(makeState(&c, &n, &p, &length, &growth_c, &growth_n, &growth_p), makeState(&work[0], &work[1], &work[2], &work[3], &work[4], &work[5], &work[6]), .{ .biological_domain_count = 1, .soil_layer_count = 2, .root_axis_count = 1, .planting_layer_index = 0, .deepest_rooted_layer_index = 1, .tip_active_by_domain_layer_axis = &.{ false, true } });
    try std.testing.expectApproxEqAbs(initial_c, c[0] + c[1] + growth_c[1], 1e-15);
    try std.testing.expectApproxEqAbs(initial_n, n[0] + n[1] + growth_n[1], 1e-15);
    try std.testing.expectApproxEqAbs(initial_p, p[0] + p[1] + growth_p[1], 1e-15);
    try std.testing.expectEqual(@as(f64, 0), c[1]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), length[0], 1e-15);
}

test "GROSUB invalid late primary deficit rolls back earlier axes" {
    var c = [_]f64{ 1, 1 };
    var n = [_]f64{ 0.1, 0.1 };
    var p = [_]f64{ 0.01, 0.01 };
    var length = [_]f64{ 1, 1 };
    var growth_c = [_]f64{ -0.5, std.math.nan(f64) };
    var growth_n = [_]f64{ 0, 0 };
    var growth_p = [_]f64{ 0, 0 };
    var work: [7][2]f64 = std.mem.zeroes([7][2]f64);
    const before_c = c;
    try std.testing.expectError(error.NonFinitePrimaryRootGrowth, apply(makeState(&c, &n, &p, &length, &growth_c, &growth_n, &growth_p), makeState(&work[0], &work[1], &work[2], &work[3], &work[4], &work[5], &work[6]), .{ .biological_domain_count = 1, .soil_layer_count = 1, .root_axis_count = 2, .planting_layer_index = 0, .deepest_rooted_layer_index = 0, .tip_active_by_domain_layer_axis = &.{ true, true } }));
    try std.testing.expectEqualDeep(before_c, c);
    try std.testing.expectEqual(@as(f64, -0.5), growth_c[0]);
    try std.testing.expect(std.math.isNan(growth_c[1]));
}
