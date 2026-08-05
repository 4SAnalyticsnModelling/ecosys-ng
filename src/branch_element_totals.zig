const std = @import("std");

pub const ElementOrgans = struct {
    leaf: []const f64,
    petiole: []const f64,
    stalk: []const f64,
    reserve: []const f64,
    husk: []const f64,
    ear: []const f64,
    grain: []const f64,
    mobile: []const f64,
};

pub const State = struct {
    c4_mobile_carbon_g_c_by_branch: []f64,
    total_carbon_g_c_by_branch: []f64,
    total_nitrogen_g_n_by_branch: []f64,
    total_phosphorus_g_p_by_branch: []f64,
    canopy_c3_carbohydrate_g_c_by_pool: []f64,
    canopy_c4_carbohydrate_g_c_by_pool: []f64,
};

pub const Inputs = struct {
    branch_count: usize,
    biochemical_pool_count: usize,
    c3_mobile_carbon_g_c_by_pool_branch: []const f64,
    c4_mobile_carbon_g_c_by_pool_branch: []const f64,
    aqueous_co2_carbon_g_c_by_pool_branch: []const f64,
    bicarbonate_carbon_g_c_by_pool_branch: []const f64,
    carbon: ElementOrgans,
    nitrogen: ElementOrgans,
    phosphorus: ElementOrgans,
};

fn poolBranchIndex(inputs: Inputs, pool: usize, branch: usize) usize {
    return pool * inputs.branch_count + branch;
}

fn copyState(destination: State, source: State) void {
    inline for (@typeInfo(State).@"struct".fields) |field| @memcpy(@field(destination, field.name), @field(source, field.name));
}

fn validateElementOrgans(organs: ElementOrgans, branch_count: usize) !void {
    inline for (@typeInfo(ElementOrgans).@"struct".fields) |field| {
        const values = @field(organs, field.name);
        if (values.len != branch_count) return error.BranchTotalsDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidBranchTotalsInput;
    }
}

fn validateInputs(inputs: Inputs) !void {
    if (inputs.branch_count == 0 or inputs.biochemical_pool_count == 0) return error.BranchTotalsDimensionMismatch;
    const pool_branch_count = std.math.mul(usize, inputs.biochemical_pool_count, inputs.branch_count) catch return error.BranchTotalsDimensionOverflow;
    inline for (.{ inputs.c3_mobile_carbon_g_c_by_pool_branch, inputs.c4_mobile_carbon_g_c_by_pool_branch, inputs.aqueous_co2_carbon_g_c_by_pool_branch, inputs.bicarbonate_carbon_g_c_by_pool_branch }) |values| {
        if (values.len != pool_branch_count) return error.BranchTotalsDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidBranchTotalsInput;
    }
    try validateElementOrgans(inputs.carbon, inputs.branch_count);
    try validateElementOrgans(inputs.nitrogen, inputs.branch_count);
    try validateElementOrgans(inputs.phosphorus, inputs.branch_count);
}

fn validateState(state: State, inputs: Inputs) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        const values = @field(state, field.name);
        const expected = if (std.mem.indexOf(u8, field.name, "by_branch") != null) inputs.branch_count else inputs.biochemical_pool_count;
        if (values.len != expected) return error.BranchTotalsDimensionMismatch;
        for (values) |value| if (!std.math.isFinite(value) or value < 0.0) return error.InvalidBranchTotalsState;
    }
}

/// Exact GROSUB 8466--8501 branch C:N:P totals. Runtime biochemical pools
/// replace legacy K=1..25. Traversal remains branch then pool; masses are g C/N/P.
pub fn apply(state: State, workspace: State, inputs: Inputs) !void {
    try validateInputs(inputs);
    try validateState(state, inputs);
    try validateState(workspace, inputs);
    copyState(workspace, state);
    for (0..inputs.branch_count) |branch| {
        workspace.c4_mobile_carbon_g_c_by_branch[branch] = 0.0;
        for (0..inputs.biochemical_pool_count) |pool| {
            const index = poolBranchIndex(inputs, pool, branch);
            workspace.c4_mobile_carbon_g_c_by_branch[branch] = workspace.c4_mobile_carbon_g_c_by_branch[branch] +
                inputs.c3_mobile_carbon_g_c_by_pool_branch[index] + inputs.c4_mobile_carbon_g_c_by_pool_branch[index] +
                inputs.aqueous_co2_carbon_g_c_by_pool_branch[index] + inputs.bicarbonate_carbon_g_c_by_pool_branch[index];
            if (branch == 0) {
                workspace.canopy_c3_carbohydrate_g_c_by_pool[pool] = 0.0;
                workspace.canopy_c4_carbohydrate_g_c_by_pool[pool] = 0.0;
            }
        }
        workspace.total_carbon_g_c_by_branch[branch] = inputs.carbon.leaf[branch] + inputs.carbon.petiole[branch] + inputs.carbon.stalk[branch] + inputs.carbon.reserve[branch] + inputs.carbon.husk[branch] + inputs.carbon.ear[branch] + inputs.carbon.grain[branch] + inputs.carbon.mobile[branch] + workspace.c4_mobile_carbon_g_c_by_branch[branch];
        workspace.total_nitrogen_g_n_by_branch[branch] = inputs.nitrogen.leaf[branch] + inputs.nitrogen.petiole[branch] + inputs.nitrogen.stalk[branch] + inputs.nitrogen.reserve[branch] + inputs.nitrogen.husk[branch] + inputs.nitrogen.ear[branch] + inputs.nitrogen.grain[branch] + inputs.nitrogen.mobile[branch];
        workspace.total_phosphorus_g_p_by_branch[branch] = inputs.phosphorus.leaf[branch] + inputs.phosphorus.petiole[branch] + inputs.phosphorus.stalk[branch] + inputs.phosphorus.reserve[branch] + inputs.phosphorus.husk[branch] + inputs.phosphorus.ear[branch] + inputs.phosphorus.grain[branch] + inputs.phosphorus.mobile[branch];
    }
    try validateState(workspace, inputs);
    copyState(state, workspace);
}

fn testOrgans(values: []const f64) ElementOrgans {
    return .{ .leaf = values, .petiole = values, .stalk = values, .reserve = values, .husk = values, .ear = values, .grain = values, .mobile = values };
}

test "GROSUB branch totals preserve pool and organ summation definitions" {
    var c4_total = [_]f64{ 99, 99 };
    var total_c = [_]f64{0} ** 2;
    var total_n = [_]f64{0} ** 2;
    var total_p = [_]f64{0} ** 2;
    var canopy_c3 = [_]f64{ 4, 5, 6 };
    var canopy_c4 = [_]f64{ 7, 8, 9 };
    var work_c4 = [_]f64{0} ** 2;
    var work_c = [_]f64{0} ** 2;
    var work_n = [_]f64{0} ** 2;
    var work_p = [_]f64{0} ** 2;
    var work_canopy3 = [_]f64{0} ** 3;
    var work_canopy4 = [_]f64{0} ** 3;
    const state: State = .{ .c4_mobile_carbon_g_c_by_branch = &c4_total, .total_carbon_g_c_by_branch = &total_c, .total_nitrogen_g_n_by_branch = &total_n, .total_phosphorus_g_p_by_branch = &total_p, .canopy_c3_carbohydrate_g_c_by_pool = &canopy_c3, .canopy_c4_carbohydrate_g_c_by_pool = &canopy_c4 };
    const workspace: State = .{ .c4_mobile_carbon_g_c_by_branch = &work_c4, .total_carbon_g_c_by_branch = &work_c, .total_nitrogen_g_n_by_branch = &work_n, .total_phosphorus_g_p_by_branch = &work_p, .canopy_c3_carbohydrate_g_c_by_pool = &work_canopy3, .canopy_c4_carbohydrate_g_c_by_pool = &work_canopy4 };
    const organ_values = [_]f64{ 1, 2 };
    const inputs: Inputs = .{ .branch_count = 2, .biochemical_pool_count = 3, .c3_mobile_carbon_g_c_by_pool_branch = &.{ 1, 2, 3, 4, 5, 6 }, .c4_mobile_carbon_g_c_by_pool_branch = &.{ 1, 1, 1, 1, 1, 1 }, .aqueous_co2_carbon_g_c_by_pool_branch = &.{ 0, 0, 0, 0, 0, 0 }, .bicarbonate_carbon_g_c_by_pool_branch = &.{ 0, 0, 0, 0, 0, 0 }, .carbon = testOrgans(&organ_values), .nitrogen = testOrgans(&organ_values), .phosphorus = testOrgans(&organ_values) };
    try apply(state, workspace, inputs);
    try std.testing.expectEqual(@as(f64, 12), c4_total[0]);
    try std.testing.expectEqual(@as(f64, 15), c4_total[1]);
    try std.testing.expectEqual(@as(f64, 20), total_c[0]);
    try std.testing.expectEqual(@as(f64, 31), total_c[1]);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, &canopy_c3);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, &canopy_c4);
}

test "GROSUB branch totals roll back a non-finite derived total" {
    var c4_total = [_]f64{3};
    var total_c = [_]f64{4};
    var total_n = [_]f64{5};
    var total_p = [_]f64{6};
    var canopy_c3 = [_]f64{7};
    var canopy_c4 = [_]f64{8};
    var work_c4 = [_]f64{0};
    var work_c = [_]f64{0};
    var work_n = [_]f64{0};
    var work_p = [_]f64{0};
    var work_canopy3 = [_]f64{0};
    var work_canopy4 = [_]f64{0};
    const state: State = .{ .c4_mobile_carbon_g_c_by_branch = &c4_total, .total_carbon_g_c_by_branch = &total_c, .total_nitrogen_g_n_by_branch = &total_n, .total_phosphorus_g_p_by_branch = &total_p, .canopy_c3_carbohydrate_g_c_by_pool = &canopy_c3, .canopy_c4_carbohydrate_g_c_by_pool = &canopy_c4 };
    const workspace: State = .{ .c4_mobile_carbon_g_c_by_branch = &work_c4, .total_carbon_g_c_by_branch = &work_c, .total_nitrogen_g_n_by_branch = &work_n, .total_phosphorus_g_p_by_branch = &work_p, .canopy_c3_carbohydrate_g_c_by_pool = &work_canopy3, .canopy_c4_carbohydrate_g_c_by_pool = &work_canopy4 };
    const huge = [_]f64{std.math.floatMax(f64)};
    const zero = [_]f64{0};
    const inputs: Inputs = .{ .branch_count = 1, .biochemical_pool_count = 1, .c3_mobile_carbon_g_c_by_pool_branch = &zero, .c4_mobile_carbon_g_c_by_pool_branch = &zero, .aqueous_co2_carbon_g_c_by_pool_branch = &zero, .bicarbonate_carbon_g_c_by_pool_branch = &zero, .carbon = testOrgans(&huge), .nitrogen = testOrgans(&zero), .phosphorus = testOrgans(&zero) };
    try std.testing.expectError(error.InvalidBranchTotalsState, apply(state, workspace, inputs));
    try std.testing.expectEqual(@as(f64, 4), total_c[0]);
    try std.testing.expectEqual(@as(f64, 7), canopy_c3[0]);
}
