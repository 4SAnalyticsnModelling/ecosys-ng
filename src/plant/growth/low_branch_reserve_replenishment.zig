const std = @import("std");

pub const State = struct {
    seasonal_storage_carbon_g_c: *f64,
    branch_reserve_carbon_g_c: []f64,
};

pub const Inputs = struct {
    branch_sapwood_carbon_g_c: []const f64,
    plant_total_sapwood_carbon_g_c: f64,
    plant_total_root_carbon_g_c: f64,
    low_reserve_threshold_g_c_per_g_sapwood_c: f64,
    exchange_fraction_per_h: f64,
    biological_timestep_h: f64,
    carbon_presence_threshold_g_c: f64,
};

/// Exact grosub.f lines 5034--5048. Branches are swept in ascending NB order;
/// each branch's WTRVCX therefore uses seasonal storage after all preceding
/// transfers. All extensive quantities are g C, rates are h^-1, and timestep
/// is h.
pub fn apply(state: State, inputs: Inputs) !usize {
    try validateDimensionsAndParameters(state, inputs);
    if (inputs.plant_total_sapwood_carbon_g_c <= inputs.carbon_presence_threshold_g_c or
        inputs.plant_total_root_carbon_g_c <= inputs.carbon_presence_threshold_g_c)
        return 0;
    if (!std.math.isFinite(state.seasonal_storage_carbon_g_c.*) or
        state.seasonal_storage_carbon_g_c.* < 0)
        return error.InvalidLowReserveReplenishmentState;

    var admitted_count: usize = 0;
    var prospective_storage_g_c = state.seasonal_storage_carbon_g_c.*;
    for (inputs.branch_sapwood_carbon_g_c, state.branch_reserve_carbon_g_c) |sapwood_g_c, reserve_g_c| {
        if (!std.math.isFinite(sapwood_g_c) or sapwood_g_c < 0)
            return error.InvalidLowReserveReplenishmentState;
        if (sapwood_g_c <= inputs.carbon_presence_threshold_g_c) continue;
        if (!std.math.isFinite(reserve_g_c) or reserve_g_c < 0)
            return error.InvalidLowReserveReplenishmentState;
        if (reserve_g_c > inputs.low_reserve_threshold_g_c_per_g_sapwood_c * sapwood_g_c)
            continue;
        admitted_count += 1;
        const transfer_g_c = calculateTransfer(
            prospective_storage_g_c,
            reserve_g_c,
            sapwood_g_c,
            inputs,
        );
        const next_reserve_g_c = reserve_g_c + transfer_g_c;
        prospective_storage_g_c -= transfer_g_c;
        if (!std.math.isFinite(next_reserve_g_c) or next_reserve_g_c < 0 or
            !std.math.isFinite(prospective_storage_g_c) or prospective_storage_g_c < 0)
            return error.LowReserveReplenishmentOverdraw;
    }
    if (admitted_count == 0) return 0;

    for (inputs.branch_sapwood_carbon_g_c, state.branch_reserve_carbon_g_c) |sapwood_g_c, *reserve_g_c| {
        if (sapwood_g_c <= inputs.carbon_presence_threshold_g_c or
            reserve_g_c.* > inputs.low_reserve_threshold_g_c_per_g_sapwood_c * sapwood_g_c)
            continue;
        const transfer_g_c = calculateTransfer(
            state.seasonal_storage_carbon_g_c.*,
            reserve_g_c.*,
            sapwood_g_c,
            inputs,
        );
        reserve_g_c.* = reserve_g_c.* + transfer_g_c;
        state.seasonal_storage_carbon_g_c.* =
            state.seasonal_storage_carbon_g_c.* - transfer_g_c;
    }
    return admitted_count;
}

fn calculateTransfer(
    seasonal_storage_g_c: f64,
    reserve_g_c: f64,
    branch_sapwood_g_c: f64,
    inputs: Inputs,
) f64 {
    const branch_sapwood_fraction = branch_sapwood_g_c /
        inputs.plant_total_sapwood_carbon_g_c;
    const allocated_root_carbon_g_c =
        inputs.plant_total_root_carbon_g_c * branch_sapwood_fraction;
    const structural_carbon_g_c = branch_sapwood_g_c + allocated_root_carbon_g_c;
    const nonnegative_reserve_g_c = @max(0.0, reserve_g_c);
    const allocated_storage_g_c = @max(0.0, seasonal_storage_g_c * branch_sapwood_fraction);
    const carbon_difference_g_c =
        (allocated_storage_g_c * branch_sapwood_g_c -
            nonnegative_reserve_g_c * allocated_root_carbon_g_c) /
        structural_carbon_g_c;
    return @max(0.0, inputs.exchange_fraction_per_h * carbon_difference_g_c *
        inputs.biological_timestep_h);
}

fn validateDimensionsAndParameters(state: State, inputs: Inputs) !void {
    if (state.branch_reserve_carbon_g_c.len == 0 or
        inputs.branch_sapwood_carbon_g_c.len != state.branch_reserve_carbon_g_c.len)
        return error.LowReserveReplenishmentDimensionMismatch;
    inline for (.{
        inputs.plant_total_sapwood_carbon_g_c,
        inputs.plant_total_root_carbon_g_c,
        inputs.low_reserve_threshold_g_c_per_g_sapwood_c,
        inputs.exchange_fraction_per_h,
        inputs.biological_timestep_h,
        inputs.carbon_presence_threshold_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidLowReserveReplenishmentInput;
    if (inputs.biological_timestep_h == 0)
        return error.InvalidLowReserveReplenishmentInput;
}

fn testInputs(sapwood: []const f64) Inputs {
    return .{
        .branch_sapwood_carbon_g_c = sapwood,
        .plant_total_sapwood_carbon_g_c = 20,
        .plant_total_root_carbon_g_c = 20,
        .low_reserve_threshold_g_c_per_g_sapwood_c = 1,
        .exchange_fraction_per_h = 0.1,
        .biological_timestep_h = 1,
        .carbon_presence_threshold_g_c = 1e-12,
    };
}

fn sum(values: []const f64) f64 {
    var result: f64 = 0;
    for (values) |value| result += value;
    return result;
}

test "GROSUB ascending branches use progressively depleted storage" {
    var storage_g_c: f64 = 10;
    var reserve = [_]f64{ 0, 0 };
    const before = storage_g_c + sum(&reserve);
    try std.testing.expectEqual(@as(usize, 2), try apply(.{
        .seasonal_storage_carbon_g_c = &storage_g_c,
        .branch_reserve_carbon_g_c = &reserve,
    }, testInputs(&.{ 10, 10 })));
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), reserve[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.24375), reserve[1], 1e-15);
    try std.testing.expectApproxEqAbs(before, storage_g_c + sum(&reserve), 1e-15);
}

test "reserve concentration equality remains admitted" {
    var storage_g_c: f64 = 10;
    var reserve = [_]f64{10};
    var inputs = testInputs(&.{10});
    inputs.plant_total_sapwood_carbon_g_c = 10;
    try std.testing.expectEqual(@as(usize, 1), try apply(.{
        .seasonal_storage_carbon_g_c = &storage_g_c,
        .branch_reserve_carbon_g_c = &reserve,
    }, inputs));
}

test "global and branch structural gates do not read inactive reserves" {
    var storage_g_c: f64 = std.math.nan(f64);
    var reserve = [_]f64{std.math.nan(f64)};
    var inputs = testInputs(&.{std.math.nan(f64)});
    inputs.plant_total_root_carbon_g_c = 0;
    try std.testing.expectEqual(@as(usize, 0), try apply(.{
        .seasonal_storage_carbon_g_c = &storage_g_c,
        .branch_reserve_carbon_g_c = &reserve,
    }, inputs));

    storage_g_c = 1;
    inputs = testInputs(&.{0});
    try std.testing.expectEqual(@as(usize, 0), try apply(.{
        .seasonal_storage_carbon_g_c = &storage_g_c,
        .branch_reserve_carbon_g_c = &reserve,
    }, inputs));
}

test "runtime branch count has no source ceiling" {
    const allocator = std.testing.allocator;
    const count = 44;
    const sapwood = try allocator.alloc(f64, count);
    defer allocator.free(sapwood);
    const reserve = try allocator.alloc(f64, count);
    defer allocator.free(reserve);
    @memset(sapwood, 1);
    @memset(reserve, 0);
    var storage_g_c: f64 = 100;
    var inputs = testInputs(sapwood);
    inputs.plant_total_sapwood_carbon_g_c = @floatFromInt(count);
    try std.testing.expectEqual(count, try apply(.{
        .seasonal_storage_carbon_g_c = &storage_g_c,
        .branch_reserve_carbon_g_c = reserve,
    }, inputs));
}

test "late invalid branch and storage overdraw leave sweep atomic" {
    var storage_g_c: f64 = 10;
    var reserve = [_]f64{ 0, std.math.nan(f64) };
    const before = reserve;
    try std.testing.expectError(
        error.InvalidLowReserveReplenishmentState,
        apply(.{
            .seasonal_storage_carbon_g_c = &storage_g_c,
            .branch_reserve_carbon_g_c = &reserve,
        }, testInputs(&.{ 10, 10 })),
    );
    try std.testing.expectEqual(@as(f64, 10), storage_g_c);
    try std.testing.expectEqualSlices(u8, std.mem.asBytes(&before), std.mem.asBytes(&reserve));

    reserve = .{ 0, 0 };
    var inputs = testInputs(&.{ 10, 10 });
    inputs.exchange_fraction_per_h = 100;
    try std.testing.expectError(
        error.LowReserveReplenishmentOverdraw,
        apply(.{
            .seasonal_storage_carbon_g_c = &storage_g_c,
            .branch_reserve_carbon_g_c = &reserve,
        }, inputs),
    );
}
