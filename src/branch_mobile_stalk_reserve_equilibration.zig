const std = @import("std");
const growth_stages = @import("plant_growth_stages.zig");

pub const ElementMass = struct {
    carbon_g_c: f64 = 0,
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
};

pub const State = struct {
    branch_mobile: []ElementMass,
    stalk_reserve: []ElementMass,
};

pub const Inputs = struct {
    growth_habit: growth_stages.GrowthHabit,
    final_seed_number_set_by_branch: []const bool,
    stem_elongation_started_by_branch: []const bool,
    leaf_and_petiole_carbon_g_c_by_branch: []const f64,
    stalk_sapwood_carbon_g_c_by_branch: []const f64,
    carbon_exchange_fraction_per_h: f64,
    nutrient_exchange_fraction_per_h: f64,
    biological_timestep_h: f64,
    maximum_nitrogen_per_carbon_g_n_per_g_c: f64,
    maximum_phosphorus_per_carbon_g_p_per_g_c: f64,
    carbon_presence_threshold_g_c: f64,
};

/// Exact GROSUB lines 4987--5020 signed leaf/mobile-to-stalk-reserve
/// equilibration. Runtime branches are swept in ascending NB order. Carbon is
/// solved first; its prospective pools enter the N/P gradients and excess
/// reserve corrections. Units are g C, g N, g P, and h.
pub fn apply(state: State, inputs: Inputs) !usize {
    try validateDimensionsAndParameters(state, inputs);
    var admitted_count: usize = 0;
    for (0..state.branch_mobile.len) |branch| {
        if (!lifecycleEnabled(inputs, branch)) continue;
        const leaf_carbon_g_c = inputs.leaf_and_petiole_carbon_g_c_by_branch[branch];
        if (!std.math.isFinite(leaf_carbon_g_c) or leaf_carbon_g_c < 0)
            return error.InvalidBranchReserveEquilibrationState;
        if (leaf_carbon_g_c <= inputs.carbon_presence_threshold_g_c) continue;
        admitted_count += 1;
        _ = try calculate(state.branch_mobile[branch], state.stalk_reserve[branch], leaf_carbon_g_c, inputs.stalk_sapwood_carbon_g_c_by_branch[branch], inputs);
    }
    if (admitted_count == 0) return 0;

    for (0..state.branch_mobile.len) |branch| {
        if (!lifecycleEnabled(inputs, branch) or
            inputs.leaf_and_petiole_carbon_g_c_by_branch[branch] <=
                inputs.carbon_presence_threshold_g_c)
            continue;
        const result = try calculate(
            state.branch_mobile[branch],
            state.stalk_reserve[branch],
            inputs.leaf_and_petiole_carbon_g_c_by_branch[branch],
            inputs.stalk_sapwood_carbon_g_c_by_branch[branch],
            inputs,
        );
        state.branch_mobile[branch].carbon_g_c = result.mobile_after.carbon_g_c;
        state.stalk_reserve[branch].carbon_g_c = result.reserve_after.carbon_g_c;
        if (result.mobile_carbon_total_g_c > inputs.carbon_presence_threshold_g_c) {
            state.branch_mobile[branch].nitrogen_g_n = result.mobile_after.nitrogen_g_n;
            state.stalk_reserve[branch].nitrogen_g_n = result.reserve_after.nitrogen_g_n;
            state.branch_mobile[branch].phosphorus_g_p = result.mobile_after.phosphorus_g_p;
            state.stalk_reserve[branch].phosphorus_g_p = result.reserve_after.phosphorus_g_p;
        }
    }
    return admitted_count;
}

const Calculation = struct {
    mobile_after: ElementMass,
    reserve_after: ElementMass,
    mobile_carbon_total_g_c: f64,
};

fn calculate(
    mobile: ElementMass,
    reserve: ElementMass,
    leaf_carbon_g_c: f64,
    sapwood_carbon_g_c: f64,
    inputs: Inputs,
) !Calculation {
    try validateMass(mobile);
    try validateMass(reserve);
    if (!std.math.isFinite(sapwood_carbon_g_c) or sapwood_carbon_g_c < 0)
        return error.InvalidBranchReserveEquilibrationState;
    const structural_carbon_g_c = leaf_carbon_g_c + sapwood_carbon_g_c;
    const mobile_carbon_total_g_c = mobile.carbon_g_c + reserve.carbon_g_c;
    const carbon_difference_g_c =
        (mobile.carbon_g_c * sapwood_carbon_g_c -
            reserve.carbon_g_c * leaf_carbon_g_c) /
        structural_carbon_g_c;
    const carbon_transfer_g_c = inputs.carbon_exchange_fraction_per_h *
        carbon_difference_g_c * inputs.biological_timestep_h;
    var mobile_after = mobile;
    var reserve_after = reserve;
    mobile_after.carbon_g_c = mobile.carbon_g_c - carbon_transfer_g_c;
    reserve_after.carbon_g_c = reserve.carbon_g_c + carbon_transfer_g_c;
    if (mobile_carbon_total_g_c > inputs.carbon_presence_threshold_g_c) {
        const nitrogen_difference_g_n =
            (mobile.nitrogen_g_n * reserve_after.carbon_g_c -
                reserve.nitrogen_g_n * mobile_after.carbon_g_c) /
            mobile_carbon_total_g_c;
        const phosphorus_difference_g_p =
            (mobile.phosphorus_g_p * reserve_after.carbon_g_c -
                reserve.phosphorus_g_p * mobile_after.carbon_g_c) /
            mobile_carbon_total_g_c;
        const nitrogen_transfer_g_n = inputs.nutrient_exchange_fraction_per_h *
            nitrogen_difference_g_n * inputs.biological_timestep_h -
            @max(0.0, reserve.nitrogen_g_n - reserve_after.carbon_g_c *
                inputs.maximum_nitrogen_per_carbon_g_n_per_g_c);
        const phosphorus_transfer_g_p = inputs.nutrient_exchange_fraction_per_h *
            phosphorus_difference_g_p * inputs.biological_timestep_h -
            @max(0.0, reserve.phosphorus_g_p - reserve_after.carbon_g_c *
                inputs.maximum_phosphorus_per_carbon_g_p_per_g_c);
        mobile_after.nitrogen_g_n = mobile.nitrogen_g_n - nitrogen_transfer_g_n;
        reserve_after.nitrogen_g_n = reserve.nitrogen_g_n + nitrogen_transfer_g_n;
        mobile_after.phosphorus_g_p = mobile.phosphorus_g_p - phosphorus_transfer_g_p;
        reserve_after.phosphorus_g_p = reserve.phosphorus_g_p + phosphorus_transfer_g_p;
    }
    try validateResultMass(mobile_after);
    try validateResultMass(reserve_after);
    return .{
        .mobile_after = mobile_after,
        .reserve_after = reserve_after,
        .mobile_carbon_total_g_c = mobile_carbon_total_g_c,
    };
}

fn lifecycleEnabled(inputs: Inputs, branch: usize) bool {
    return switch (inputs.growth_habit) {
        .annual => inputs.final_seed_number_set_by_branch[branch],
        .perennial => inputs.stem_elongation_started_by_branch[branch],
    };
}

fn validateDimensionsAndParameters(state: State, inputs: Inputs) !void {
    const count = state.branch_mobile.len;
    if (count == 0 or state.stalk_reserve.len != count or
        inputs.final_seed_number_set_by_branch.len != count or
        inputs.stem_elongation_started_by_branch.len != count or
        inputs.leaf_and_petiole_carbon_g_c_by_branch.len != count or
        inputs.stalk_sapwood_carbon_g_c_by_branch.len != count)
        return error.BranchReserveEquilibrationDimensionMismatch;
    inline for (.{
        inputs.carbon_exchange_fraction_per_h,
        inputs.nutrient_exchange_fraction_per_h,
        inputs.biological_timestep_h,
        inputs.maximum_nitrogen_per_carbon_g_n_per_g_c,
        inputs.maximum_phosphorus_per_carbon_g_p_per_g_c,
        inputs.carbon_presence_threshold_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidBranchReserveEquilibrationInput;
    if (inputs.biological_timestep_h == 0 or
        inputs.maximum_nitrogen_per_carbon_g_n_per_g_c == 0 or
        inputs.maximum_phosphorus_per_carbon_g_p_per_g_c == 0)
        return error.InvalidBranchReserveEquilibrationInput;
}

fn validateMass(mass: ElementMass) !void {
    inline for (.{ mass.carbon_g_c, mass.nitrogen_g_n, mass.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidBranchReserveEquilibrationState;
}

fn validateResultMass(mass: ElementMass) !void {
    inline for (.{ mass.carbon_g_c, mass.nitrogen_g_n, mass.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.BranchReserveEquilibrationOverdraw;
}

fn testInputs(count: usize) Inputs {
    return .{
        .growth_habit = .perennial,
        .final_seed_number_set_by_branch = if (count == 2) &.{ true, true } else &.{},
        .stem_elongation_started_by_branch = if (count == 2) &.{ true, true } else &.{},
        .leaf_and_petiole_carbon_g_c_by_branch = if (count == 2) &.{ 10, 20 } else &.{},
        .stalk_sapwood_carbon_g_c_by_branch = if (count == 2) &.{ 10, 20 } else &.{},
        .carbon_exchange_fraction_per_h = 0.1,
        .nutrient_exchange_fraction_per_h = 0.1,
        .biological_timestep_h = 1,
        .maximum_nitrogen_per_carbon_g_n_per_g_c = 0.2,
        .maximum_phosphorus_per_carbon_g_p_per_g_c = 0.02,
        .carbon_presence_threshold_g_c = 1e-12,
    };
}

fn total(mobile: []const ElementMass, reserve: []const ElementMass) ElementMass {
    var result: ElementMass = .{};
    for (mobile) |mass| {
        result.carbon_g_c += mass.carbon_g_c;
        result.nitrogen_g_n += mass.nitrogen_g_n;
        result.phosphorus_g_p += mass.phosphorus_g_p;
    }
    for (reserve) |mass| {
        result.carbon_g_c += mass.carbon_g_c;
        result.nitrogen_g_n += mass.nitrogen_g_n;
        result.phosphorus_g_p += mass.phosphorus_g_p;
    }
    return result;
}

test "GROSUB branch sweep conserves signed C N P equilibration" {
    var mobile = [_]ElementMass{
        .{ .carbon_g_c = 8, .nitrogen_g_n = 0.8, .phosphorus_g_p = 0.08 },
        .{ .carbon_g_c = 4, .nitrogen_g_n = 0.4, .phosphorus_g_p = 0.04 },
    };
    var reserve = [_]ElementMass{
        .{ .carbon_g_c = 2, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.02 },
        .{ .carbon_g_c = 6, .nitrogen_g_n = 0.6, .phosphorus_g_p = 0.06 },
    };
    const before = total(&mobile, &reserve);
    try std.testing.expectEqual(@as(usize, 2), try apply(.{
        .branch_mobile = &mobile,
        .stalk_reserve = &reserve,
    }, testInputs(2)));
    const after = total(&mobile, &reserve);
    try std.testing.expectApproxEqAbs(before.carbon_g_c, after.carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(before.nitrogen_g_n, after.nitrogen_g_n, 1e-15);
    try std.testing.expectApproxEqAbs(before.phosphorus_g_p, after.phosphorus_g_p, 1e-15);
    try std.testing.expect(reserve[0].carbon_g_c > 2);
    try std.testing.expect(reserve[1].carbon_g_c < 6);
}

test "nutrient gradient uses carbon pools after signed carbon transfer" {
    var mobile = [_]ElementMass{.{ .carbon_g_c = 8, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 }};
    var reserve = [_]ElementMass{.{ .carbon_g_c = 2, .nitrogen_g_n = 0, .phosphorus_g_p = 0 }};
    var inputs = testInputs(2);
    inputs.final_seed_number_set_by_branch = &.{true};
    inputs.stem_elongation_started_by_branch = &.{true};
    inputs.leaf_and_petiole_carbon_g_c_by_branch = &.{10};
    inputs.stalk_sapwood_carbon_g_c_by_branch = &.{10};
    try std.testing.expectEqual(@as(usize, 1), try apply(.{
        .branch_mobile = &mobile,
        .stalk_reserve = &reserve,
    }, inputs));
    // Carbon transfer is 0.3, so N gradient uses reserve C=2.3 and mobile C=7.7.
    try std.testing.expectApproxEqAbs(@as(f64, 0.023), reserve[0].nitrogen_g_n, 1e-15);
}

test "lifecycle and leaf-presence gates do not read dormant pools" {
    var mobile = [_]ElementMass{.{ .carbon_g_c = std.math.nan(f64) }};
    var reserve = [_]ElementMass{.{ .carbon_g_c = std.math.nan(f64) }};
    var inputs = testInputs(2);
    inputs.final_seed_number_set_by_branch = &.{false};
    inputs.stem_elongation_started_by_branch = &.{false};
    inputs.leaf_and_petiole_carbon_g_c_by_branch = &.{std.math.nan(f64)};
    inputs.stalk_sapwood_carbon_g_c_by_branch = &.{std.math.nan(f64)};
    try std.testing.expectEqual(@as(usize, 0), try apply(.{
        .branch_mobile = &mobile,
        .stalk_reserve = &reserve,
    }, inputs));
}

test "runtime branch count has no source ceiling" {
    const allocator = std.testing.allocator;
    const count = 42;
    const mobile = try allocator.alloc(ElementMass, count);
    defer allocator.free(mobile);
    const reserve = try allocator.alloc(ElementMass, count);
    defer allocator.free(reserve);
    const flags = try allocator.alloc(bool, count);
    defer allocator.free(flags);
    const leaf = try allocator.alloc(f64, count);
    defer allocator.free(leaf);
    const sapwood = try allocator.alloc(f64, count);
    defer allocator.free(sapwood);
    @memset(mobile, .{ .carbon_g_c = 2, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.02 });
    @memset(reserve, .{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01 });
    @memset(flags, true);
    @memset(leaf, 1);
    @memset(sapwood, 1);
    var inputs = testInputs(2);
    inputs.final_seed_number_set_by_branch = flags;
    inputs.stem_elongation_started_by_branch = flags;
    inputs.leaf_and_petiole_carbon_g_c_by_branch = leaf;
    inputs.stalk_sapwood_carbon_g_c_by_branch = sapwood;
    try std.testing.expectEqual(count, try apply(.{
        .branch_mobile = mobile,
        .stalk_reserve = reserve,
    }, inputs));
}

test "late invalid admitted branch and overdraw leave sweep atomic" {
    var mobile = [_]ElementMass{
        .{ .carbon_g_c = 8, .nitrogen_g_n = 0.8, .phosphorus_g_p = 0.08 },
        .{ .carbon_g_c = 4, .nitrogen_g_n = 0.4, .phosphorus_g_p = std.math.nan(f64) },
    };
    var reserve = [_]ElementMass{ .{}, .{} };
    const mobile_before = mobile;
    try std.testing.expectError(
        error.InvalidBranchReserveEquilibrationState,
        apply(.{ .branch_mobile = &mobile, .stalk_reserve = &reserve }, testInputs(2)),
    );
    try std.testing.expectEqualSlices(u8, std.mem.asBytes(&mobile_before), std.mem.asBytes(&mobile));

    mobile = .{ .{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01 }, .{} };
    reserve = .{ .{}, .{} };
    var inputs = testInputs(2);
    inputs.carbon_exchange_fraction_per_h = 100;
    try std.testing.expectError(
        error.BranchReserveEquilibrationOverdraw,
        apply(.{ .branch_mobile = &mobile, .stalk_reserve = &reserve }, inputs),
    );
}
