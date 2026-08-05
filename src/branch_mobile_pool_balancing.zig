const std = @import("std");

pub const State = struct {
    mobile_carbon_g_c_by_branch: []f64,
    mobile_nitrogen_g_n_by_branch: []f64,
    mobile_phosphorus_g_p_by_branch: []f64,
};

pub const Inputs = struct {
    branch_count: usize,
    branch_is_alive: []const bool,
    hours_since_germination_by_branch: []const f64,
    leaf_petiole_carbon_g_c_by_branch: []const f64,
    storage_remobilization_end_h: f64,
    minimum_pool_g_c: f64,
};

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateState(state: State, branch_count: usize) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        if (values.len != branch_count) return error.BranchMobilePoolDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidBranchMobilePoolState;
    }
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.branch_count == 0 or
        inputs.branch_is_alive.len != inputs.branch_count or
        inputs.hours_since_germination_by_branch.len != inputs.branch_count or
        inputs.leaf_petiole_carbon_g_c_by_branch.len != inputs.branch_count)
        return error.BranchMobilePoolDimensionMismatch;
    if (!std.math.isFinite(inputs.storage_remobilization_end_h) or
        !std.math.isFinite(inputs.minimum_pool_g_c) or
        inputs.minimum_pool_g_c < 0.0)
        return error.InvalidBranchMobilePoolInput;
    for (inputs.hours_since_germination_by_branch) |value| if (!std.math.isFinite(value)) return error.InvalidBranchMobilePoolInput;
    for (inputs.leaf_petiole_carbon_g_c_by_branch) |value| if (!std.math.isFinite(value)) return error.InvalidBranchMobilePoolInput;
}

/// Exact GROSUB 7994--8031 balancing of non-structural C:N:P among living
/// branches after seasonal storage remobilization. Runtime arrays replace NBR;
/// masses are g C/N/P and time is h. The fixed legacy exchange fraction is 0.01.
pub fn apply(allocator: std.mem.Allocator, state: State, workspace: State, inputs: Inputs) !void {
    try validateState(state, inputs.branch_count);
    try validateState(workspace, inputs.branch_count);
    try validateInputs(inputs);
    copyState(workspace, state);
    if (inputs.branch_count <= 1) return;

    const leaf_carbon_g_c = try allocator.alloc(f64, inputs.branch_count);
    defer allocator.free(leaf_carbon_g_c);
    const mobile_carbon_g_c = try allocator.alloc(f64, inputs.branch_count);
    defer allocator.free(mobile_carbon_g_c);
    const mobile_nitrogen_g_n = try allocator.alloc(f64, inputs.branch_count);
    defer allocator.free(mobile_nitrogen_g_n);
    const mobile_phosphorus_g_p = try allocator.alloc(f64, inputs.branch_count);
    defer allocator.free(mobile_phosphorus_g_p);
    @memset(leaf_carbon_g_c, 0.0);
    @memset(mobile_carbon_g_c, 0.0);
    @memset(mobile_nitrogen_g_n, 0.0);
    @memset(mobile_phosphorus_g_p, 0.0);

    var total_leaf_carbon_g_c: f64 = 0.0;
    var total_mobile_carbon_g_c: f64 = 0.0;
    var total_mobile_nitrogen_g_n: f64 = 0.0;
    var total_mobile_phosphorus_g_p: f64 = 0.0;
    for (0..inputs.branch_count) |branch| {
        if (!inputs.branch_is_alive[branch] or inputs.hours_since_germination_by_branch[branch] <= inputs.storage_remobilization_end_h) continue;
        leaf_carbon_g_c[branch] = @max(0.0, inputs.leaf_petiole_carbon_g_c_by_branch[branch]);
        mobile_carbon_g_c[branch] = @max(0.0, state.mobile_carbon_g_c_by_branch[branch]);
        mobile_nitrogen_g_n[branch] = @max(0.0, state.mobile_nitrogen_g_n_by_branch[branch]);
        mobile_phosphorus_g_p[branch] = @max(0.0, state.mobile_phosphorus_g_p_by_branch[branch]);
        total_leaf_carbon_g_c += leaf_carbon_g_c[branch];
        total_mobile_carbon_g_c += mobile_carbon_g_c[branch];
        total_mobile_nitrogen_g_n += mobile_nitrogen_g_n[branch];
        total_mobile_phosphorus_g_p += mobile_phosphorus_g_p[branch];
    }
    for (0..inputs.branch_count) |branch| {
        if (!inputs.branch_is_alive[branch] or inputs.hours_since_germination_by_branch[branch] <= inputs.storage_remobilization_end_h) continue;
        if (total_leaf_carbon_g_c > inputs.minimum_pool_g_c and total_mobile_carbon_g_c > inputs.minimum_pool_g_c) {
            const carbon_difference_g_c2 = total_mobile_carbon_g_c * leaf_carbon_g_c[branch] - mobile_carbon_g_c[branch] * total_leaf_carbon_g_c;
            const nitrogen_difference_g_n_c = total_mobile_nitrogen_g_n * mobile_carbon_g_c[branch] - mobile_nitrogen_g_n[branch] * total_mobile_carbon_g_c;
            const phosphorus_difference_g_p_c = total_mobile_phosphorus_g_p * mobile_carbon_g_c[branch] - mobile_phosphorus_g_p[branch] * total_mobile_carbon_g_c;
            workspace.mobile_carbon_g_c_by_branch[branch] += 0.01 * carbon_difference_g_c2 / total_leaf_carbon_g_c;
            workspace.mobile_nitrogen_g_n_by_branch[branch] += 0.01 * nitrogen_difference_g_n_c / total_mobile_carbon_g_c;
            workspace.mobile_phosphorus_g_p_by_branch[branch] += 0.01 * phosphorus_difference_g_p_c / total_mobile_carbon_g_c;
        }
    }
    try validateState(workspace, inputs.branch_count);
    copyState(state, workspace);
}

test "GROSUB branch balancing conserves eligible mobile C N and P" {
    var carbon = [_]f64{ 8.0, 2.0, 4.0 };
    var nitrogen = [_]f64{ 0.4, 0.3, 0.2 };
    var phosphorus = [_]f64{ 0.08, 0.03, 0.02 };
    var work_carbon = [_]f64{0} ** 3;
    var work_nitrogen = [_]f64{0} ** 3;
    var work_phosphorus = [_]f64{0} ** 3;
    const state: State = .{ .mobile_carbon_g_c_by_branch = &carbon, .mobile_nitrogen_g_n_by_branch = &nitrogen, .mobile_phosphorus_g_p_by_branch = &phosphorus };
    const workspace: State = .{ .mobile_carbon_g_c_by_branch = &work_carbon, .mobile_nitrogen_g_n_by_branch = &work_nitrogen, .mobile_phosphorus_g_p_by_branch = &work_phosphorus };
    const inputs: Inputs = .{
        .branch_count = 3,
        .branch_is_alive = &.{ true, true, false },
        .hours_since_germination_by_branch = &.{ 100.0, 100.0, 100.0 },
        .leaf_petiole_carbon_g_c_by_branch = &.{ 2.0, 8.0, 5.0 },
        .storage_remobilization_end_h = 50.0,
        .minimum_pool_g_c = 1e-12,
    };
    const before = [_]f64{ carbon[0] + carbon[1], nitrogen[0] + nitrogen[1], phosphorus[0] + phosphorus[1] };
    try apply(std.testing.allocator, state, workspace, inputs);
    try std.testing.expectApproxEqAbs(before[0], carbon[0] + carbon[1], 1e-14);
    try std.testing.expectApproxEqAbs(before[1], nitrogen[0] + nitrogen[1], 1e-14);
    try std.testing.expectApproxEqAbs(before[2], phosphorus[0] + phosphorus[1], 1e-14);
    try std.testing.expectEqual(@as(f64, 4.0), carbon[2]);
    try std.testing.expect(carbon[0] < 8.0 and carbon[1] > 2.0);
}
