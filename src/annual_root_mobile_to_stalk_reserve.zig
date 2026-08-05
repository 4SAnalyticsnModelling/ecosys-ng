const std = @import("std");
const growth_stages = @import("plant_growth_stages.zig");

pub const ElementMass = struct {
    carbon_g_c: f64 = 0,
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
};

pub const State = struct {
    root_mobile_by_layer: []ElementMass,
    stalk_reserve_by_branch: []ElementMass,
};

pub const Inputs = struct {
    growth_habit: growth_stages.GrowthHabit,
    final_seed_number_set_by_branch: []const bool,
    stalk_sapwood_carbon_g_c_by_branch: []const f64,
    soil_layer_thickness_m: []const f64,
    active_root_carbon_g_c_by_layer: []const f64,
    root_woody_carbon_fraction: f64,
    minimum_soil_layer_thickness_m: f64,
    carbon_exchange_fraction_per_h: f64,
    nutrient_exchange_fraction_per_h: f64,
    biological_timestep_h: f64,
    carbon_presence_threshold_g_c: f64,
};

pub const Result = struct {
    transferred: ElementMass,
    branch_layer_transfer_count: usize,
};

/// Exact GROSUB lines 5068--5110. Runtime branches are outer and root layers
/// inner, both ascending. Carbon publication precedes the nutrient gradients.
/// XFRP deliberately retains the source ZPOOLR (root N) cap. Units are m, h,
/// g C, g N, and g P.
pub fn apply(
    allocator: std.mem.Allocator,
    state: State,
    inputs: Inputs,
) !Result {
    try validateDimensionsAndParameters(state, inputs);
    if (inputs.growth_habit != .annual) return .{
        .transferred = .{},
        .branch_layer_transfer_count = 0,
    };

    const root_scratch = try allocator.dupe(ElementMass, state.root_mobile_by_layer);
    defer allocator.free(root_scratch);
    const reserve_scratch = try allocator.dupe(ElementMass, state.stalk_reserve_by_branch);
    defer allocator.free(reserve_scratch);
    const result = try sweep(root_scratch, reserve_scratch, inputs);
    if (result.branch_layer_transfer_count == 0) return result;
    const committed = try sweep(
        state.root_mobile_by_layer,
        state.stalk_reserve_by_branch,
        inputs,
    );
    std.debug.assert(committed.branch_layer_transfer_count == result.branch_layer_transfer_count);
    return committed;
}

fn sweep(
    root_mobile_by_layer: []ElementMass,
    stalk_reserve_by_branch: []ElementMass,
    inputs: Inputs,
) !Result {
    var result: Result = .{ .transferred = .{}, .branch_layer_transfer_count = 0 };
    for (stalk_reserve_by_branch, 0..) |*reserve, branch| {
        if (!inputs.final_seed_number_set_by_branch[branch]) continue;
        const sapwood_carbon_g_c = inputs.stalk_sapwood_carbon_g_c_by_branch[branch];
        if (!std.math.isFinite(sapwood_carbon_g_c) or sapwood_carbon_g_c < 0)
            return error.InvalidAnnualRootReserveState;
        try validateMass(reserve.*);
        for (root_mobile_by_layer, 0..) |*root, layer| {
            const thickness_m = inputs.soil_layer_thickness_m[layer];
            const active_root_carbon_g_c = inputs.active_root_carbon_g_c_by_layer[layer];
            if (!std.math.isFinite(thickness_m) or thickness_m < 0 or
                !std.math.isFinite(active_root_carbon_g_c) or active_root_carbon_g_c < 0)
                return error.InvalidAnnualRootReserveState;
            if (thickness_m <= inputs.minimum_soil_layer_thickness_m or
                active_root_carbon_g_c <= inputs.carbon_presence_threshold_g_c)
                continue;
            try validateMass(root.*);
            result.branch_layer_transfer_count += 1;

            const effective_root_carbon_g_c = @max(
                inputs.carbon_presence_threshold_g_c,
                active_root_carbon_g_c * inputs.root_woody_carbon_fraction,
            );
            const structural_carbon_g_c = effective_root_carbon_g_c + sapwood_carbon_g_c;
            const carbon_difference_g_c =
                (root.carbon_g_c * sapwood_carbon_g_c -
                    reserve.carbon_g_c * effective_root_carbon_g_c) /
                structural_carbon_g_c;
            const carbon_transfer_g_c = @max(0.0, inputs.carbon_exchange_fraction_per_h * carbon_difference_g_c) *
                inputs.biological_timestep_h;
            root.carbon_g_c = root.carbon_g_c - carbon_transfer_g_c;
            reserve.carbon_g_c = reserve.carbon_g_c + carbon_transfer_g_c;
            result.transferred.carbon_g_c += carbon_transfer_g_c;
            const mobile_carbon_total_g_c = root.carbon_g_c + reserve.carbon_g_c;
            if (mobile_carbon_total_g_c > inputs.carbon_presence_threshold_g_c) {
                const nitrogen_difference_g_n =
                    (root.nitrogen_g_n * reserve.carbon_g_c -
                        reserve.nitrogen_g_n * root.carbon_g_c) /
                    mobile_carbon_total_g_c;
                const phosphorus_difference_g_p =
                    (root.phosphorus_g_p * reserve.carbon_g_c -
                        reserve.phosphorus_g_p * root.carbon_g_c) /
                    mobile_carbon_total_g_c;
                const nitrogen_transfer_g_n = @max(0.0, @min(root.nitrogen_g_n, inputs.nutrient_exchange_fraction_per_h * nitrogen_difference_g_n)) *
                    inputs.biological_timestep_h;
                const phosphorus_transfer_g_p = @max(0.0, @min(root.nitrogen_g_n, inputs.nutrient_exchange_fraction_per_h * phosphorus_difference_g_p)) *
                    inputs.biological_timestep_h;
                root.nitrogen_g_n = root.nitrogen_g_n - nitrogen_transfer_g_n;
                reserve.nitrogen_g_n = reserve.nitrogen_g_n + nitrogen_transfer_g_n;
                root.phosphorus_g_p = root.phosphorus_g_p - phosphorus_transfer_g_p;
                reserve.phosphorus_g_p = reserve.phosphorus_g_p + phosphorus_transfer_g_p;
                result.transferred.nitrogen_g_n += nitrogen_transfer_g_n;
                result.transferred.phosphorus_g_p += phosphorus_transfer_g_p;
            }
            try validateResultMass(root.*);
            try validateResultMass(reserve.*);
            try validateResultMass(result.transferred);
        }
    }
    return result;
}

fn validateDimensionsAndParameters(state: State, inputs: Inputs) !void {
    const branch_count = state.stalk_reserve_by_branch.len;
    const layer_count = state.root_mobile_by_layer.len;
    if (branch_count == 0 or layer_count == 0 or
        inputs.final_seed_number_set_by_branch.len != branch_count or
        inputs.stalk_sapwood_carbon_g_c_by_branch.len != branch_count or
        inputs.soil_layer_thickness_m.len != layer_count or
        inputs.active_root_carbon_g_c_by_layer.len != layer_count)
        return error.AnnualRootReserveDimensionMismatch;
    inline for (.{
        inputs.root_woody_carbon_fraction,
        inputs.minimum_soil_layer_thickness_m,
        inputs.carbon_exchange_fraction_per_h,
        inputs.nutrient_exchange_fraction_per_h,
        inputs.biological_timestep_h,
        inputs.carbon_presence_threshold_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidAnnualRootReserveInput;
    if (inputs.root_woody_carbon_fraction > 1 or inputs.biological_timestep_h == 0)
        return error.InvalidAnnualRootReserveInput;
}

fn validateMass(mass: ElementMass) !void {
    inline for (.{ mass.carbon_g_c, mass.nitrogen_g_n, mass.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidAnnualRootReserveState;
}

fn validateResultMass(mass: ElementMass) !void {
    inline for (.{ mass.carbon_g_c, mass.nitrogen_g_n, mass.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.AnnualRootReserveOverdraw;
}

fn testInputs(branch_count: usize, layer_count: usize) Inputs {
    return .{
        .growth_habit = .annual,
        .final_seed_number_set_by_branch = if (branch_count == 2) &.{ true, true } else &.{},
        .stalk_sapwood_carbon_g_c_by_branch = if (branch_count == 2) &.{ 5, 10 } else &.{},
        .soil_layer_thickness_m = if (layer_count == 2) &.{ 0.1, 0.2 } else &.{},
        .active_root_carbon_g_c_by_layer = if (layer_count == 2) &.{ 10, 20 } else &.{},
        .root_woody_carbon_fraction = 0.5,
        .minimum_soil_layer_thickness_m = 0.01,
        .carbon_exchange_fraction_per_h = 0.1,
        .nutrient_exchange_fraction_per_h = 0.1,
        .biological_timestep_h = 1,
        .carbon_presence_threshold_g_c = 1e-12,
    };
}

fn total(roots: []const ElementMass, reserves: []const ElementMass) ElementMass {
    var result: ElementMass = .{};
    for (roots) |mass| {
        result.carbon_g_c += mass.carbon_g_c;
        result.nitrogen_g_n += mass.nitrogen_g_n;
        result.phosphorus_g_p += mass.phosphorus_g_p;
    }
    for (reserves) |mass| {
        result.carbon_g_c += mass.carbon_g_c;
        result.nitrogen_g_n += mass.nitrogen_g_n;
        result.phosphorus_g_p += mass.phosphorus_g_p;
    }
    return result;
}

test "GROSUB branch-outer layer-inner sweep conserves C N P" {
    var roots = [_]ElementMass{
        .{ .carbon_g_c = 8, .nitrogen_g_n = 0.8, .phosphorus_g_p = 0.08 },
        .{ .carbon_g_c = 12, .nitrogen_g_n = 1.2, .phosphorus_g_p = 0.12 },
    };
    var reserves = [_]ElementMass{ .{}, .{} };
    const before = total(&roots, &reserves);
    const result = try apply(std.testing.allocator, .{
        .root_mobile_by_layer = &roots,
        .stalk_reserve_by_branch = &reserves,
    }, testInputs(2, 2));
    const after = total(&roots, &reserves);
    try std.testing.expectEqual(@as(usize, 4), result.branch_layer_transfer_count);
    try std.testing.expectApproxEqAbs(before.carbon_g_c, after.carbon_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(before.nitrogen_g_n, after.nitrogen_g_n, 1e-14);
    try std.testing.expectApproxEqAbs(before.phosphorus_g_p, after.phosphorus_g_p, 1e-14);
}

test "source phosphorus transfer is capped by root nitrogen" {
    var roots = [_]ElementMass{.{ .carbon_g_c = 10, .nitrogen_g_n = 0.01, .phosphorus_g_p = 10 }};
    var reserves = [_]ElementMass{.{}};
    var inputs = testInputs(2, 2);
    inputs.final_seed_number_set_by_branch = &.{true};
    inputs.stalk_sapwood_carbon_g_c_by_branch = &.{10};
    inputs.soil_layer_thickness_m = &.{0.1};
    inputs.active_root_carbon_g_c_by_layer = &.{10};
    inputs.nutrient_exchange_fraction_per_h = 100;
    const result = try apply(std.testing.allocator, .{
        .root_mobile_by_layer = &roots,
        .stalk_reserve_by_branch = &reserves,
    }, inputs);
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), result.transferred.phosphorus_g_p, 1e-15);
}

test "annual final-seed and layer gates avoid inactive pool reads" {
    var roots = [_]ElementMass{.{ .carbon_g_c = std.math.nan(f64) }};
    var reserves = [_]ElementMass{.{ .carbon_g_c = std.math.nan(f64) }};
    var inputs = testInputs(2, 2);
    inputs.final_seed_number_set_by_branch = &.{false};
    inputs.stalk_sapwood_carbon_g_c_by_branch = &.{std.math.nan(f64)};
    inputs.soil_layer_thickness_m = &.{0};
    inputs.active_root_carbon_g_c_by_layer = &.{std.math.nan(f64)};
    const result = try apply(std.testing.allocator, .{
        .root_mobile_by_layer = &roots,
        .stalk_reserve_by_branch = &reserves,
    }, inputs);
    try std.testing.expectEqual(@as(usize, 0), result.branch_layer_transfer_count);
}

test "runtime branch and layer dimensions have no source ceilings" {
    const allocator = std.testing.allocator;
    const branch_count = 37;
    const layer_count = 41;
    const roots = try allocator.alloc(ElementMass, layer_count);
    defer allocator.free(roots);
    const reserves = try allocator.alloc(ElementMass, branch_count);
    defer allocator.free(reserves);
    const flags = try allocator.alloc(bool, branch_count);
    defer allocator.free(flags);
    const sapwood = try allocator.alloc(f64, branch_count);
    defer allocator.free(sapwood);
    const thickness = try allocator.alloc(f64, layer_count);
    defer allocator.free(thickness);
    const active_root = try allocator.alloc(f64, layer_count);
    defer allocator.free(active_root);
    @memset(roots, .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 });
    @memset(reserves, .{});
    @memset(flags, true);
    @memset(sapwood, 1);
    @memset(thickness, 0.1);
    @memset(active_root, 1);
    var inputs = testInputs(2, 2);
    inputs.final_seed_number_set_by_branch = flags;
    inputs.stalk_sapwood_carbon_g_c_by_branch = sapwood;
    inputs.soil_layer_thickness_m = thickness;
    inputs.active_root_carbon_g_c_by_layer = active_root;
    const result = try apply(allocator, .{
        .root_mobile_by_layer = roots,
        .stalk_reserve_by_branch = reserves,
    }, inputs);
    try std.testing.expectEqual(branch_count * layer_count, result.branch_layer_transfer_count);
}

test "late invalid layer and phosphorus overdraw leave real state atomic" {
    var roots = [_]ElementMass{
        .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 },
        .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = std.math.nan(f64) },
    };
    var reserves = [_]ElementMass{ .{}, .{} };
    const roots_before = roots;
    try std.testing.expectError(
        error.InvalidAnnualRootReserveState,
        apply(std.testing.allocator, .{
            .root_mobile_by_layer = &roots,
            .stalk_reserve_by_branch = &reserves,
        }, testInputs(2, 2)),
    );
    try std.testing.expectEqualSlices(u8, std.mem.asBytes(&roots_before), std.mem.asBytes(&roots));

    roots = .{
        .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 1e-6 },
        .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 },
    };
    var inputs = testInputs(2, 2);
    inputs.nutrient_exchange_fraction_per_h = 100;
    try std.testing.expectError(
        error.AnnualRootReserveOverdraw,
        apply(std.testing.allocator, .{
            .root_mobile_by_layer = &roots,
            .stalk_reserve_by_branch = &reserves,
        }, inputs),
    );
}
