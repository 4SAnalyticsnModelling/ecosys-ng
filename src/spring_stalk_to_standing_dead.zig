const std = @import("std");

pub const EnableStatus = enum { enabled, disabled };
pub const Turnover = enum { fully_deciduous, other };
pub const RootProfile = enum { shallow, intermediate, deep, deeper };
pub const Elements = struct { carbon: f64, nitrogen: f64, phosphorus: f64 };
pub const Kinetics = struct { carbon: []const f64, nitrogen: []const f64, phosphorus: []const f64 };
pub const Nodes = struct {
    height_m: []f64,
    senescing_internode_height_m: []f64,
    internode_carbon_g_c: []f64,
    internode_nitrogen_g_n: []f64,
    internode_phosphorus_g_p: []f64,
};
pub const State = struct {
    branch_stalk: *Elements,
    residual_stalk: *Elements,
    branch_reserve: *Elements,
    plant_storage: *Elements,
    standing_dead_carbon_g_c: []f64,
    standing_dead_nitrogen_g_n: []f64,
    standing_dead_phosphorus_g_p: []f64,
    nodes: Nodes,
};
pub const Inputs = struct {
    leafout_status: EnableStatus,
    perennial: bool,
    accumulated_leafout_h: f64,
    required_leafout_h: f64,
    turnover: Turnover,
    root_profile: RootProfile,
    stalk_litter_kinetics: Kinetics,
};

/// GROSUB lines 4350--4378. At completed perennial leafout for fully
/// deciduous or intermediate-root plants, transfers branch stalk C/N/P to
/// standing-dead runtime kinetic pools, transfers reserve C/N/P to plant
/// storage, clears branch/residual/reserve pools, then clears all runtime nodes.
pub fn apply(state: State, inputs: Inputs) !bool {
    if (inputs.leafout_status != .enabled or !inputs.perennial or
        inputs.accumulated_leafout_h < inputs.required_leafout_h)
        return false;
    if (inputs.turnover != .fully_deciduous and inputs.root_profile != .intermediate) return false;
    const count = try validate(state, inputs);
    const stalk = state.branch_stalk.*;
    const reserve = state.branch_reserve.*;
    const storage_after: Elements = .{
        .carbon = state.plant_storage.carbon + reserve.carbon,
        .nitrogen = state.plant_storage.nitrogen + reserve.nitrogen,
        .phosphorus = state.plant_storage.phosphorus + reserve.phosphorus,
    };
    inline for (@typeInfo(Elements).@"struct".fields) |field| {
        const value = @field(storage_after, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSpringStalkTransferResult;
    }
    for (0..count) |kinetic| {
        inline for (.{
            state.standing_dead_carbon_g_c[kinetic] + inputs.stalk_litter_kinetics.carbon[kinetic] * stalk.carbon,
            state.standing_dead_nitrogen_g_n[kinetic] + inputs.stalk_litter_kinetics.nitrogen[kinetic] * stalk.nitrogen,
            state.standing_dead_phosphorus_g_p[kinetic] + inputs.stalk_litter_kinetics.phosphorus[kinetic] * stalk.phosphorus,
        }) |updated| if (!std.math.isFinite(updated) or updated < 0) return error.InvalidSpringStalkTransferResult;
    }

    for (0..count) |kinetic| {
        state.standing_dead_carbon_g_c[kinetic] += inputs.stalk_litter_kinetics.carbon[kinetic] * stalk.carbon;
        state.standing_dead_nitrogen_g_n[kinetic] += inputs.stalk_litter_kinetics.nitrogen[kinetic] * stalk.nitrogen;
        state.standing_dead_phosphorus_g_p[kinetic] += inputs.stalk_litter_kinetics.phosphorus[kinetic] * stalk.phosphorus;
    }
    state.plant_storage.* = storage_after;
    state.branch_stalk.* = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 };
    state.residual_stalk.* = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 };
    state.branch_reserve.* = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 };
    inline for (@typeInfo(Nodes).@"struct".fields) |field| @memset(@field(state.nodes, field.name), 0);
    return true;
}

fn validate(state: State, inputs: Inputs) !usize {
    inline for (.{ inputs.accumulated_leafout_h, inputs.required_leafout_h }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSpringStalkTransferInput;
    inline for (.{ state.branch_stalk.*, state.residual_stalk.*, state.branch_reserve.*, state.plant_storage.* }) |elements| inline for (@typeInfo(Elements).@"struct".fields) |field| {
        const value = @field(elements, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSpringStalkTransferState;
    };
    const node_count = state.nodes.height_m.len;
    inline for (@typeInfo(Nodes).@"struct".fields[1..]) |field| if (@field(state.nodes, field.name).len != node_count) return error.SpringStalkTransferNodeDimensionMismatch;
    inline for (@typeInfo(Nodes).@"struct".fields) |field| for (@field(state.nodes, field.name)) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSpringStalkTransferState;
    const count = inputs.stalk_litter_kinetics.carbon.len;
    if (count == 0) return error.ZeroSpringStalkKineticPools;
    inline for (.{ inputs.stalk_litter_kinetics.nitrogen, inputs.stalk_litter_kinetics.phosphorus, state.standing_dead_carbon_g_c, state.standing_dead_nitrogen_g_n, state.standing_dead_phosphorus_g_p }) |values| if (values.len != count) return error.SpringStalkKineticDimensionMismatch;
    inline for (@typeInfo(Kinetics).@"struct".fields) |field| {
        var total: f64 = 0;
        for (@field(inputs.stalk_litter_kinetics, field.name)) |fraction| {
            if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidSpringStalkKinetics;
            total += fraction;
        }
        if (!std.math.isFinite(total) or @abs(total - 1) > 1e-12) return error.InvalidSpringStalkKinetics;
    }
    inline for (.{ state.standing_dead_carbon_g_c, state.standing_dead_nitrogen_g_n, state.standing_dead_phosphorus_g_p }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSpringStalkTransferState;
    return count;
}

fn makeState(pools: *[4]Elements, standing: *[3][3]f64, node_store: *[5][31]f64, node_count: usize) State {
    return .{ .branch_stalk = &pools[0], .residual_stalk = &pools[1], .branch_reserve = &pools[2], .plant_storage = &pools[3], .standing_dead_carbon_g_c = &standing[0], .standing_dead_nitrogen_g_n = &standing[1], .standing_dead_phosphorus_g_p = &standing[2], .nodes = .{ .height_m = node_store[0][0..node_count], .senescing_internode_height_m = node_store[1][0..node_count], .internode_carbon_g_c = node_store[2][0..node_count], .internode_nitrogen_g_n = node_store[3][0..node_count], .internode_phosphorus_g_p = node_store[4][0..node_count] } };
}
fn testInputs() Inputs {
    const f = &[_]f64{ 0.2, 0.3, 0.5 };
    return .{ .leafout_status = .enabled, .perennial = true, .accumulated_leafout_h = 10, .required_leafout_h = 10, .turnover = .fully_deciduous, .root_profile = .deep, .stalk_litter_kinetics = .{ .carbon = f, .nitrogen = f, .phosphorus = f } };
}

test "spring transfer conserves stalk and reserve C N P and clears runtime nodes" {
    var pools = [_]Elements{ .{ .carbon = 10, .nitrogen = 2, .phosphorus = 0.5 }, .{ .carbon = 3, .nitrogen = 0.6, .phosphorus = 0.1 }, .{ .carbon = 4, .nitrogen = 0.8, .phosphorus = 0.2 }, .{ .carbon = 1, .nitrogen = 0.2, .phosphorus = 0.05 } };
    var standing: [3][3]f64 = @splat(@splat(0));
    var nodes: [5][31]f64 = @splat(@splat(1));
    try std.testing.expect(try apply(makeState(&pools, &standing, &nodes, 31), testInputs()));
    try std.testing.expectApproxEqAbs(@as(f64, 10), standing[0][0] + standing[0][1] + standing[0][2], 1e-15);
    try std.testing.expectEqual(@as(f64, 5), pools[3].carbon);
    try std.testing.expectEqual(@as(f64, 0), pools[0].carbon);
    try std.testing.expectEqual(@as(f64, 0), nodes[4][30]);
}
test "intermediate roots activate transfer independent of turnover" {
    var pools: [4]Elements = @splat(.{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 });
    var standing: [3][3]f64 = @splat(@splat(0));
    var nodes: [5][31]f64 = @splat(@splat(0));
    var input = testInputs();
    input.turnover = .other;
    input.root_profile = .intermediate;
    try std.testing.expect(try apply(makeState(&pools, &standing, &nodes, 1), input));
}
test "inactive gate is strict no-op before state reads" {
    var pools: [4]Elements = @splat(.{ .carbon = std.math.nan(f64), .nitrogen = 0, .phosphorus = 0 });
    var standing: [3][3]f64 = @splat(@splat(std.math.nan(f64)));
    var nodes: [5][31]f64 = @splat(@splat(std.math.nan(f64)));
    var input = testInputs();
    input.turnover = .other;
    input.root_profile = .deep;
    try std.testing.expect(!try apply(makeState(&pools, &standing, &nodes, 31), input));
}
test "late invalid standing dead destination is atomic" {
    var pools = [_]Elements{ .{ .carbon = 10, .nitrogen = 2, .phosphorus = 0.5 }, .{ .carbon = 3, .nitrogen = 0.6, .phosphorus = 0.1 }, .{ .carbon = 4, .nitrogen = 0.8, .phosphorus = 0.2 }, .{ .carbon = 1, .nitrogen = 0.2, .phosphorus = 0.05 } };
    var standing: [3][3]f64 = @splat(@splat(0));
    standing[2][2] = std.math.nan(f64);
    var nodes: [5][31]f64 = @splat(@splat(1));
    try std.testing.expectError(error.InvalidSpringStalkTransferState, apply(makeState(&pools, &standing, &nodes, 31), testInputs()));
    try std.testing.expectEqual(@as(f64, 10), pools[0].carbon);
    try std.testing.expectEqual(@as(f64, 1), nodes[0][30]);
    try std.testing.expectEqual(@as(f64, 0), standing[0][0]);
}
