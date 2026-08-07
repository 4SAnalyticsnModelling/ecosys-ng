const std = @import("std");

pub const EnableStatus = enum { enabled, disabled };
pub const Turnover = enum { fully_deciduous, other };
pub const Elements = struct { carbon: f64, nitrogen: f64, phosphorus: f64 };
pub const Fractions = struct { carbon: []const f64, nitrogen: []const f64, phosphorus: []const f64 };
pub const Kinetics = struct { woody: Fractions, leaf: Fractions, sheath: Fractions };
pub const Composition = struct {
    carbon: [2]f64,
    leaf_nitrogen: [2]f64,
    leaf_phosphorus: [2]f64,
    sheath_nitrogen: [2]f64,
    sheath_phosphorus: [2]f64,
};
pub const Litter = struct {
    woody_carbon_g_c: []f64,
    woody_nitrogen_g_n: []f64,
    woody_phosphorus_g_p: []f64,
    nonwoody_carbon_g_c: []f64,
    nonwoody_nitrogen_g_n: []f64,
    nonwoody_phosphorus_g_p: []f64,
};
pub const Nodes = struct {
    leaf_area_m2: []f64,
    sheath_height_m: []f64,
    leaf_carbon_g_c: []f64,
    leaf_protein_g: []f64,
    leaf_nitrogen_g_n: []f64,
    leaf_phosphorus_g_p: []f64,
    sheath_carbon_g_c: []f64,
    sheath_protein_g: []f64,
    sheath_nitrogen_g_n: []f64,
    sheath_phosphorus_g_p: []f64,
};
pub const State = struct {
    phenological_node: *f64,
    appeared_leaf_count: *f64,
    next_leaf_index: *usize,
    latest_leaf_index: *usize,
    no_grain_fill_h: *f64,
    branch_leaf_area_m2: *f64,
    branch_leaf: *Elements,
    branch_sheath: *Elements,
    nodes: Nodes,
    litter: Litter,
};
pub const Inputs = struct {
    leafout_status: EnableStatus,
    accumulated_leafout_h: f64,
    required_leafout_h: f64,
    perennial: bool,
    turnover: Turnover,
    initial_phenological_node: f64,
    kinetics: Kinetics,
    composition: Composition,
};

/// grosub.f lines 4249--4297. For completed perennial leafout with IBTYP=0,
/// resets leaf development, routes all branch leaf/sheath C/N/P to runtime
/// woody and nonwoody litter pools, clears branch totals, then clears every
/// runtime logical node. Full projection validation makes the reset atomic.
pub fn apply(state: State, inputs: Inputs) !bool {
    if (inputs.leafout_status != .enabled or !inputs.perennial or inputs.accumulated_leafout_h < inputs.required_leafout_h) return false;
    if (inputs.turnover != .fully_deciduous) return false;
    const count = try validate(state, inputs);
    const leaf = state.branch_leaf.*;
    const sheath = state.branch_sheath.*;
    for (0..count) |k| try validateAddition(state, inputs, leaf, sheath, k);

    state.phenological_node.* = inputs.initial_phenological_node;
    state.appeared_leaf_count.* = 0;
    state.next_leaf_index.* = 1;
    state.latest_leaf_index.* = 1;
    state.no_grain_fill_h.* = 0;
    for (0..count) |k| commitAddition(state, inputs, leaf, sheath, k);
    state.branch_leaf_area_m2.* = 0;
    state.branch_leaf.* = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 };
    state.branch_sheath.* = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 };
    inline for (@typeInfo(Nodes).@"struct".fields) |field| @memset(@field(state.nodes, field.name), 0);
    return true;
}

fn validate(state: State, inputs: Inputs) !usize {
    inline for (.{ inputs.accumulated_leafout_h, inputs.required_leafout_h, inputs.initial_phenological_node, state.branch_leaf_area_m2.* }) |v| if (!std.math.isFinite(v) or v < 0) return error.InvalidSpringLeafResetInput;
    inline for (.{ state.branch_leaf.*, state.branch_sheath.* }) |e| inline for (@typeInfo(Elements).@"struct".fields) |f| {
        const v = @field(e, f.name);
        if (!std.math.isFinite(v) or v < 0) return error.InvalidSpringLeafResetState;
    };
    const nodes = state.nodes.leaf_area_m2.len;
    inline for (@typeInfo(Nodes).@"struct".fields[1..]) |f| if (@field(state.nodes, f.name).len != nodes) return error.SpringLeafResetNodeDimensionMismatch;
    inline for (@typeInfo(Nodes).@"struct".fields) |f| for (@field(state.nodes, f.name)) |v| if (!std.math.isFinite(v) or v < 0) return error.InvalidSpringLeafResetState;
    const count = inputs.kinetics.woody.carbon.len;
    if (count == 0) return error.ZeroSpringLeafResetKineticPools;
    inline for (.{ inputs.kinetics.woody, inputs.kinetics.leaf, inputs.kinetics.sheath }) |kin| inline for (@typeInfo(Fractions).@"struct".fields) |f| {
        const values = @field(kin, f.name);
        if (values.len != count) return error.SpringLeafResetKineticDimensionMismatch;
        var total: f64 = 0;
        for (values) |v| {
            if (!std.math.isFinite(v) or v < 0 or v > 1) return error.InvalidSpringLeafResetKinetics;
            total += v;
        }
        if (@abs(total - 1) > 1e-12) return error.InvalidSpringLeafResetKinetics;
    };
    inline for (@typeInfo(Litter).@"struct".fields) |f| if (@field(state.litter, f.name).len != count) return error.SpringLeafResetKineticDimensionMismatch;
    inline for (.{ inputs.composition.carbon, inputs.composition.leaf_nitrogen, inputs.composition.leaf_phosphorus, inputs.composition.sheath_nitrogen, inputs.composition.sheath_phosphorus }) |pair| {
        inline for (pair) |v| if (!std.math.isFinite(v) or v < 0 or v > 1) return error.InvalidSpringLeafResetComposition;
        if (@abs(pair[0] + pair[1] - 1) > 1e-12) return error.InvalidSpringLeafResetComposition;
    }
    return count;
}
fn additions(inputs: Inputs, leaf: Elements, sheath: Elements, k: usize) [6]f64 {
    return .{
        inputs.kinetics.woody.carbon[k] * (leaf.carbon + sheath.carbon) * inputs.composition.carbon[0],
        inputs.kinetics.woody.nitrogen[k] * (leaf.nitrogen * inputs.composition.leaf_nitrogen[0] + sheath.nitrogen * inputs.composition.sheath_nitrogen[0]),
        inputs.kinetics.woody.phosphorus[k] * (leaf.phosphorus * inputs.composition.leaf_phosphorus[0] + sheath.phosphorus * inputs.composition.sheath_phosphorus[0]),
        inputs.kinetics.leaf.carbon[k] * leaf.carbon * inputs.composition.carbon[1] + inputs.kinetics.sheath.carbon[k] * sheath.carbon * inputs.composition.carbon[1],
        inputs.kinetics.leaf.nitrogen[k] * leaf.nitrogen * inputs.composition.leaf_nitrogen[1] + inputs.kinetics.sheath.nitrogen[k] * sheath.nitrogen * inputs.composition.sheath_nitrogen[1],
        inputs.kinetics.leaf.phosphorus[k] * leaf.phosphorus * inputs.composition.leaf_phosphorus[1] + inputs.kinetics.sheath.phosphorus[k] * sheath.phosphorus * inputs.composition.sheath_phosphorus[1],
    };
}
fn validateAddition(state: State, inputs: Inputs, leaf: Elements, sheath: Elements, k: usize) !void {
    const add = additions(inputs, leaf, sheath, k);
    inline for (@typeInfo(Litter).@"struct".fields, 0..) |f, i| {
        const old = @field(state.litter, f.name)[k];
        if (!std.math.isFinite(old) or old < 0 or !std.math.isFinite(add[i]) or add[i] < 0 or !std.math.isFinite(old + add[i])) return error.InvalidSpringLeafResetLitterResult;
    }
}
fn commitAddition(state: State, inputs: Inputs, leaf: Elements, sheath: Elements, k: usize) void {
    const add = additions(inputs, leaf, sheath, k);
    inline for (@typeInfo(Litter).@"struct".fields, 0..) |f, i| @field(state.litter, f.name)[k] += add[i];
}

fn fractions(v: []const f64) Fractions {
    return .{ .carbon = v, .nitrogen = v, .phosphorus = v };
}
fn testInputs() Inputs {
    const f = &[_]f64{ 0.25, 0.75 };
    return .{ .leafout_status = .enabled, .accumulated_leafout_h = 10, .required_leafout_h = 10, .perennial = true, .turnover = .fully_deciduous, .initial_phenological_node = 2, .kinetics = .{ .woody = fractions(f), .leaf = fractions(f), .sheath = fractions(f) }, .composition = .{ .carbon = .{ 0.2, 0.8 }, .leaf_nitrogen = .{ 0.3, 0.7 }, .leaf_phosphorus = .{ 0.4, 0.6 }, .sheath_nitrogen = .{ 0.5, 0.5 }, .sheath_phosphorus = .{ 0.6, 0.4 } } };
}
fn makeState(dev: *[3]f64, index: *[2]usize, leaf_area: *f64, leaf: *Elements, sheath: *Elements, node_store: *[10][31]f64, litter_store: *[6][2]f64, node_count: usize) State {
    return .{ .phenological_node = &dev[0], .appeared_leaf_count = &dev[1], .next_leaf_index = &index[0], .latest_leaf_index = &index[1], .no_grain_fill_h = &dev[2], .branch_leaf_area_m2 = leaf_area, .branch_leaf = leaf, .branch_sheath = sheath, .nodes = .{ .leaf_area_m2 = node_store[0][0..node_count], .sheath_height_m = node_store[1][0..node_count], .leaf_carbon_g_c = node_store[2][0..node_count], .leaf_protein_g = node_store[3][0..node_count], .leaf_nitrogen_g_n = node_store[4][0..node_count], .leaf_phosphorus_g_p = node_store[5][0..node_count], .sheath_carbon_g_c = node_store[6][0..node_count], .sheath_protein_g = node_store[7][0..node_count], .sheath_nitrogen_g_n = node_store[8][0..node_count], .sheath_phosphorus_g_p = node_store[9][0..node_count] }, .litter = .{ .woody_carbon_g_c = &litter_store[0], .woody_nitrogen_g_n = &litter_store[1], .woody_phosphorus_g_p = &litter_store[2], .nonwoody_carbon_g_c = &litter_store[3], .nonwoody_nitrogen_g_n = &litter_store[4], .nonwoody_phosphorus_g_p = &litter_store[5] } };
}

test "spring reset conserves branch leaf and sheath C N P in runtime litter pools" {
    var dev = [_]f64{ 9, 8, 7 };
    var index = [_]usize{ 9, 9 };
    var area: f64 = 3;
    var leaf: Elements = .{ .carbon = 10, .nitrogen = 2, .phosphorus = 0.5 };
    var sheath: Elements = .{ .carbon = 5, .nitrogen = 1, .phosphorus = 0.25 };
    var nodes: [10][31]f64 = @splat(@splat(1));
    var litter: [6][2]f64 = @splat(@splat(0));
    try std.testing.expect(try apply(makeState(&dev, &index, &area, &leaf, &sheath, &nodes, &litter, 31), testInputs()));
    try std.testing.expectEqual(@as(f64, 2), dev[0]);
    try std.testing.expectEqual(@as(usize, 1), index[0]);
    try std.testing.expectEqual(@as(f64, 0), nodes[9][30]);
    inline for (@typeInfo(Elements).@"struct".fields, 0..) |field, i| {
        const total = litter[i][0] + litter[i][1] + litter[i + 3][0] + litter[i + 3][1];
        const expected = @field(Elements{ .carbon = 15, .nitrogen = 3, .phosphorus = 0.75 }, field.name);
        try std.testing.expectApproxEqAbs(expected, total, 1e-15);
    }
}
test "disabled and nondeciduous gates are strict no-ops" {
    var dev = [_]f64{ 9, 8, 7 };
    var index = [_]usize{ 9, 9 };
    var area: f64 = 3;
    var leaf: Elements = .{ .carbon = 1, .nitrogen = 1, .phosphorus = 1 };
    var sheath = leaf;
    var nodes: [10][31]f64 = @splat(@splat(std.math.nan(f64)));
    var litter: [6][2]f64 = @splat(@splat(std.math.nan(f64)));
    var input = testInputs();
    input.leafout_status = .disabled;
    try std.testing.expect(!try apply(makeState(&dev, &index, &area, &leaf, &sheath, &nodes, &litter, 31), input));
    try std.testing.expectEqual(@as(f64, 9), dev[0]);
    input = testInputs();
    input.turnover = .other;
    try std.testing.expect(!try apply(makeState(&dev, &index, &area, &leaf, &sheath, &nodes, &litter, 31), input));
}
test "late invalid litter leaves development branch and nodes unchanged" {
    var dev = [_]f64{ 9, 8, 7 };
    var index = [_]usize{ 9, 9 };
    var area: f64 = 3;
    var leaf: Elements = .{ .carbon = 10, .nitrogen = 2, .phosphorus = 0.5 };
    var sheath: Elements = .{ .carbon = 5, .nitrogen = 1, .phosphorus = 0.25 };
    var nodes: [10][31]f64 = @splat(@splat(1));
    var litter: [6][2]f64 = @splat(@splat(0));
    litter[5][1] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSpringLeafResetLitterResult, apply(makeState(&dev, &index, &area, &leaf, &sheath, &nodes, &litter, 31), testInputs()));
    try std.testing.expectEqual(@as(f64, 9), dev[0]);
    try std.testing.expectEqual(@as(f64, 10), leaf.carbon);
    try std.testing.expectEqual(@as(f64, 1), nodes[0][30]);
}
