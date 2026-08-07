const std = @import("std");
const growth_stages = @import("../lifecycle/growth_stages.zig");

pub const ElementMass = struct {
    carbon_g_c: f64 = 0,
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
};

pub const State = struct {
    seasonal_storage: *ElementMass,
    branch_reserve: []ElementMass,
    branch_mobile: []ElementMass,
};

pub const Inputs = struct {
    growth_habit: growth_stages.GrowthHabit,
    shoot_remobilization_enabled_by_branch: []const bool,
    exchange_fraction_per_h: f64,
    biological_timestep_h: f64,
    minimum_nitrogen_per_carbon_g_n_per_g_c: f64,
    maximum_nitrogen_per_carbon_g_n_per_g_c: f64,
    minimum_phosphorus_per_carbon_g_p_per_g_c: f64,
    maximum_phosphorus_per_carbon_g_p_per_g_c: f64,
};

/// Exact grosub.f lines 4926--4964. Runtime branches are processed in ascending
/// NB order. Within each branch, stalk reserve precedes shoot mobile, and each
/// category commits C, N, then P to seasonal storage. Mass units are g C,
/// g N, and g P; rates are h^-1 and the timestep is h.
pub fn apply(state: State, inputs: Inputs) !usize {
    try validateDimensionsAndParameters(state, inputs);
    if (inputs.growth_habit != .perennial) return 0;

    var admitted_count: usize = 0;
    var prospective_storage = state.seasonal_storage.*;
    try validateMass(prospective_storage);
    for (inputs.shoot_remobilization_enabled_by_branch, state.branch_reserve, state.branch_mobile) |enabled, reserve, mobile| {
        if (!enabled) continue;
        admitted_count += 1;
        try validateMass(reserve);
        try validateMass(mobile);
        const reserve_transfer = calculateTransfer(reserve, inputs);
        try validateTransferAndRemainder(reserve, reserve_transfer);
        prospective_storage.carbon_g_c += reserve_transfer.carbon_g_c;
        prospective_storage.nitrogen_g_n += reserve_transfer.nitrogen_g_n;
        prospective_storage.phosphorus_g_p += reserve_transfer.phosphorus_g_p;
        try validateResultMass(prospective_storage);

        const mobile_transfer = calculateTransfer(mobile, inputs);
        try validateTransferAndRemainder(mobile, mobile_transfer);
        prospective_storage.carbon_g_c += mobile_transfer.carbon_g_c;
        prospective_storage.nitrogen_g_n += mobile_transfer.nitrogen_g_n;
        prospective_storage.phosphorus_g_p += mobile_transfer.phosphorus_g_p;
        try validateResultMass(prospective_storage);
    }
    if (admitted_count == 0) return 0;

    for (inputs.shoot_remobilization_enabled_by_branch, state.branch_reserve, state.branch_mobile) |enabled, *reserve, *mobile| {
        if (!enabled) continue;
        const reserve_transfer = calculateTransfer(reserve.*, inputs);
        reserve.carbon_g_c = reserve.carbon_g_c - reserve_transfer.carbon_g_c;
        state.seasonal_storage.carbon_g_c = state.seasonal_storage.carbon_g_c +
            reserve_transfer.carbon_g_c;
        reserve.nitrogen_g_n = reserve.nitrogen_g_n - reserve_transfer.nitrogen_g_n;
        state.seasonal_storage.nitrogen_g_n = state.seasonal_storage.nitrogen_g_n +
            reserve_transfer.nitrogen_g_n;
        reserve.phosphorus_g_p = reserve.phosphorus_g_p - reserve_transfer.phosphorus_g_p;
        state.seasonal_storage.phosphorus_g_p = state.seasonal_storage.phosphorus_g_p +
            reserve_transfer.phosphorus_g_p;

        const mobile_transfer = calculateTransfer(mobile.*, inputs);
        mobile.carbon_g_c = mobile.carbon_g_c - mobile_transfer.carbon_g_c;
        state.seasonal_storage.carbon_g_c = state.seasonal_storage.carbon_g_c +
            mobile_transfer.carbon_g_c;
        mobile.nitrogen_g_n = mobile.nitrogen_g_n - mobile_transfer.nitrogen_g_n;
        state.seasonal_storage.nitrogen_g_n = state.seasonal_storage.nitrogen_g_n +
            mobile_transfer.nitrogen_g_n;
        mobile.phosphorus_g_p = mobile.phosphorus_g_p - mobile_transfer.phosphorus_g_p;
        state.seasonal_storage.phosphorus_g_p = state.seasonal_storage.phosphorus_g_p +
            mobile_transfer.phosphorus_g_p;
    }
    return admitted_count;
}

fn calculateTransfer(pool: ElementMass, inputs: Inputs) ElementMass {
    const unconstrained_carbon_g_c = inputs.exchange_fraction_per_h *
        @max(0.0, pool.carbon_g_c) * inputs.biological_timestep_h;
    const unconstrained_nitrogen_g_n = inputs.exchange_fraction_per_h *
        @max(0.0, pool.nitrogen_g_n) * inputs.biological_timestep_h;
    const unconstrained_phosphorus_g_p = inputs.exchange_fraction_per_h *
        @max(0.0, pool.phosphorus_g_p) * inputs.biological_timestep_h;
    const carbon_g_c = @min(unconstrained_carbon_g_c, @min(unconstrained_nitrogen_g_n /
        inputs.minimum_nitrogen_per_carbon_g_n_per_g_c, unconstrained_phosphorus_g_p /
        inputs.minimum_phosphorus_per_carbon_g_p_per_g_c));
    const nitrogen_g_n = @min(unconstrained_nitrogen_g_n, @min(carbon_g_c * inputs.maximum_nitrogen_per_carbon_g_n_per_g_c, unconstrained_phosphorus_g_p *
        inputs.maximum_nitrogen_per_carbon_g_n_per_g_c /
        inputs.minimum_phosphorus_per_carbon_g_p_per_g_c));
    const phosphorus_g_p = @min(unconstrained_phosphorus_g_p, @min(carbon_g_c * inputs.maximum_phosphorus_per_carbon_g_p_per_g_c, unconstrained_nitrogen_g_n *
        inputs.maximum_phosphorus_per_carbon_g_p_per_g_c /
        inputs.minimum_nitrogen_per_carbon_g_n_per_g_c));
    return .{
        .carbon_g_c = carbon_g_c,
        .nitrogen_g_n = nitrogen_g_n,
        .phosphorus_g_p = phosphorus_g_p,
    };
}

fn validateDimensionsAndParameters(state: State, inputs: Inputs) !void {
    const branch_count = state.branch_reserve.len;
    if (branch_count == 0 or state.branch_mobile.len != branch_count or
        inputs.shoot_remobilization_enabled_by_branch.len != branch_count)
        return error.PerennialStorageRemobilizationDimensionMismatch;
    inline for (.{
        inputs.exchange_fraction_per_h,
        inputs.biological_timestep_h,
        inputs.minimum_nitrogen_per_carbon_g_n_per_g_c,
        inputs.maximum_nitrogen_per_carbon_g_n_per_g_c,
        inputs.minimum_phosphorus_per_carbon_g_p_per_g_c,
        inputs.maximum_phosphorus_per_carbon_g_p_per_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidPerennialStorageRemobilizationInput;
    if (inputs.biological_timestep_h == 0 or
        inputs.minimum_nitrogen_per_carbon_g_n_per_g_c == 0 or
        inputs.minimum_phosphorus_per_carbon_g_p_per_g_c == 0 or
        inputs.maximum_nitrogen_per_carbon_g_n_per_g_c <
            inputs.minimum_nitrogen_per_carbon_g_n_per_g_c or
        inputs.maximum_phosphorus_per_carbon_g_p_per_g_c <
            inputs.minimum_phosphorus_per_carbon_g_p_per_g_c)
        return error.InvalidPerennialStorageRemobilizationInput;
}

fn validateMass(mass: ElementMass) !void {
    inline for (.{ mass.carbon_g_c, mass.nitrogen_g_n, mass.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidPerennialStorageRemobilizationState;
}

fn validateTransferAndRemainder(pool: ElementMass, transfer: ElementMass) !void {
    inline for (.{
        transfer.carbon_g_c,
        transfer.nitrogen_g_n,
        transfer.phosphorus_g_p,
        pool.carbon_g_c - transfer.carbon_g_c,
        pool.nitrogen_g_n - transfer.nitrogen_g_n,
        pool.phosphorus_g_p - transfer.phosphorus_g_p,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.PerennialStorageRemobilizationOverdraw;
}

fn validateResultMass(mass: ElementMass) !void {
    inline for (.{ mass.carbon_g_c, mass.nitrogen_g_n, mass.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.NonFinitePerennialStorageRemobilizationResult;
}

fn testInputs(enabled: []const bool) Inputs {
    return .{
        .growth_habit = .perennial,
        .shoot_remobilization_enabled_by_branch = enabled,
        .exchange_fraction_per_h = 0.1,
        .biological_timestep_h = 1,
        .minimum_nitrogen_per_carbon_g_n_per_g_c = 0.05,
        .maximum_nitrogen_per_carbon_g_n_per_g_c = 0.2,
        .minimum_phosphorus_per_carbon_g_p_per_g_c = 0.005,
        .maximum_phosphorus_per_carbon_g_p_per_g_c = 0.02,
    };
}

fn total(storage: ElementMass, reserve: []const ElementMass, mobile: []const ElementMass) ElementMass {
    var result = storage;
    for (reserve) |mass| {
        result.carbon_g_c += mass.carbon_g_c;
        result.nitrogen_g_n += mass.nitrogen_g_n;
        result.phosphorus_g_p += mass.phosphorus_g_p;
    }
    for (mobile) |mass| {
        result.carbon_g_c += mass.carbon_g_c;
        result.nitrogen_g_n += mass.nitrogen_g_n;
        result.phosphorus_g_p += mass.phosphorus_g_p;
    }
    return result;
}

test "GROSUB ascending reserve then mobile transfers conserve C N P" {
    var storage: ElementMass = .{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01 };
    var reserve = [_]ElementMass{
        .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 },
        .{ .carbon_g_c = 20, .nitrogen_g_n = 2, .phosphorus_g_p = 0.2 },
    };
    var mobile = [_]ElementMass{
        .{ .carbon_g_c = 4, .nitrogen_g_n = 0.4, .phosphorus_g_p = 0.04 },
        .{ .carbon_g_c = 8, .nitrogen_g_n = 0.8, .phosphorus_g_p = 0.08 },
    };
    const before = total(storage, &reserve, &mobile);
    try std.testing.expectEqual(@as(usize, 2), try apply(.{
        .seasonal_storage = &storage,
        .branch_reserve = &reserve,
        .branch_mobile = &mobile,
    }, testInputs(&.{ true, true })));
    const after = total(storage, &reserve, &mobile);
    try std.testing.expectApproxEqAbs(before.carbon_g_c, after.carbon_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(before.nitrogen_g_n, after.nitrogen_g_n, 1e-14);
    try std.testing.expectApproxEqAbs(before.phosphorus_g_p, after.phosphorus_g_p, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 5.2), storage.carbon_g_c, 1e-15);
}

test "element constraints retain source three-way minima" {
    var storage: ElementMass = .{};
    var reserve = [_]ElementMass{.{ .carbon_g_c = 100, .nitrogen_g_n = 0.1, .phosphorus_g_p = 10 }};
    var mobile = [_]ElementMass{.{}};
    try std.testing.expectEqual(@as(usize, 1), try apply(.{
        .seasonal_storage = &storage,
        .branch_reserve = &reserve,
        .branch_mobile = &mobile,
    }, testInputs(&.{true})));
    // C is constrained by unconstrained N / CNMN = 0.01 / 0.05.
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), storage.carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), storage.nitrogen_g_n, 1e-15);
}

test "annual and disabled gates do not read branch pools" {
    var storage: ElementMass = .{};
    var reserve = [_]ElementMass{.{ .carbon_g_c = std.math.nan(f64) }};
    var mobile = [_]ElementMass{.{ .carbon_g_c = std.math.nan(f64) }};
    var inputs = testInputs(&.{false});
    try std.testing.expectEqual(@as(usize, 0), try apply(.{
        .seasonal_storage = &storage,
        .branch_reserve = &reserve,
        .branch_mobile = &mobile,
    }, inputs));
    inputs.growth_habit = .annual;
    inputs.shoot_remobilization_enabled_by_branch = &.{true};
    try std.testing.expectEqual(@as(usize, 0), try apply(.{
        .seasonal_storage = &storage,
        .branch_reserve = &reserve,
        .branch_mobile = &mobile,
    }, inputs));
}

test "runtime branch sweep has no five-species or branch ceiling" {
    const allocator = std.testing.allocator;
    const count = 43;
    const flags = try allocator.alloc(bool, count);
    defer allocator.free(flags);
    const reserve = try allocator.alloc(ElementMass, count);
    defer allocator.free(reserve);
    const mobile = try allocator.alloc(ElementMass, count);
    defer allocator.free(mobile);
    @memset(flags, true);
    @memset(reserve, .{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01 });
    @memset(mobile, .{});
    var storage: ElementMass = .{};
    try std.testing.expectEqual(count, try apply(.{
        .seasonal_storage = &storage,
        .branch_reserve = reserve,
        .branch_mobile = mobile,
    }, testInputs(flags)));
}

test "late invalid enabled branch and overdraw leave sweep atomic" {
    var storage: ElementMass = .{};
    var reserve = [_]ElementMass{
        .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 },
        .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = std.math.nan(f64) },
    };
    var mobile = [_]ElementMass{ .{}, .{} };
    const reserve_before = reserve;
    try std.testing.expectError(
        error.InvalidPerennialStorageRemobilizationState,
        apply(.{
            .seasonal_storage = &storage,
            .branch_reserve = &reserve,
            .branch_mobile = &mobile,
        }, testInputs(&.{ true, true })),
    );
    try std.testing.expectEqual(@as(ElementMass, .{}), storage);
    try std.testing.expectEqualSlices(u8, std.mem.asBytes(&reserve_before), std.mem.asBytes(&reserve));

    reserve = .{ .{ .carbon_g_c = 1, .nitrogen_g_n = 1, .phosphorus_g_p = 1 }, .{} };
    var inputs = testInputs(&.{ true, false });
    inputs.exchange_fraction_per_h = 2;
    try std.testing.expectError(
        error.PerennialStorageRemobilizationOverdraw,
        apply(.{
            .seasonal_storage = &storage,
            .branch_reserve = &reserve,
            .branch_mobile = &mobile,
        }, inputs),
    );
}
