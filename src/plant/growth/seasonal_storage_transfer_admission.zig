const std = @import("std");
const execution_calendar = @import("../../driver/execution_calendar_date.zig");
const growth_stages = @import("../lifecycle/growth_stages.zig");

pub const SimulationDate = struct {
    day_of_year: u16,
    year: u16,
};

pub const State = struct {
    time_since_germination_h: *f64,
    remobilization_counter_initialized: *bool,
};

pub const Inputs = struct {
    growth_habit: growth_stages.GrowthHabit,
    lifecycle_initialized: bool,
    current_date: SimulationDate,
    planting_date: ?SimulationDate,
    current_branch_accumulated_leafoff_h: f64,
    current_branch_required_leafoff_h: f64,
    leafoff_remobilization_start_fraction: f64,
    main_branch_accumulated_leafout_h: f64,
    current_branch_required_leafout_h: f64,
    root_structural_carbon_g_c_by_layer: []const f64,
    root_mobile_carbon_g_c_by_layer: []const f64,
};

pub const RootCarbonTotals = struct {
    structural_carbon_g_c: f64,
    mobile_carbon_g_c: f64,
};

/// Exact grosub.f lines 4632--4654 seasonal-storage transfer admission and
/// preparation. Root C is g C and the counter is h. The main branch supplies
/// VRNS while the current branch supplies VRNL/VRNF/VRNX, matching the source
/// subscripts. A null planting date explicitly represents the source sentinel.
pub fn apply(state: State, inputs: Inputs) !?RootCarbonTotals {
    try validateGateInputs(inputs);
    const before_leafoff_remobilization =
        inputs.current_branch_accumulated_leafoff_h <
        inputs.leafoff_remobilization_start_fraction *
            inputs.current_branch_required_leafoff_h;
    const planting_gate = if (inputs.planting_date) |planting_date|
        inputs.current_date.year == planting_date.year and
            inputs.current_date.day_of_year >= planting_date.day_of_year and
            before_leafoff_remobilization
    else
        false;
    const admitted = (inputs.growth_habit == .annual and
        !inputs.lifecycle_initialized) or planting_gate or
        (inputs.main_branch_accumulated_leafout_h >=
            inputs.current_branch_required_leafout_h and
            before_leafoff_remobilization);
    if (!admitted) return null;

    if (inputs.root_structural_carbon_g_c_by_layer.len == 0 or
        inputs.root_structural_carbon_g_c_by_layer.len !=
            inputs.root_mobile_carbon_g_c_by_layer.len)
        return error.SeasonalStorageAdmissionDimensionMismatch;
    var structural_carbon_g_c: f64 = 0;
    var mobile_carbon_g_c: f64 = 0;
    for (inputs.root_structural_carbon_g_c_by_layer, inputs.root_mobile_carbon_g_c_by_layer) |structural, mobile| {
        inline for (.{ structural, mobile }) |value| {
            if (!std.math.isFinite(value))
                return error.NonFiniteSeasonalStorageAdmissionState;
            // The source AMAX1 protects totals from corrupt negative pools;
            // ecosys-ng reports the impossible pool at its origin instead.
            if (value < 0)
                return error.NegativeSeasonalStorageAdmissionState;
        }
        structural_carbon_g_c += @max(0.0, structural);
        mobile_carbon_g_c += @max(0.0, mobile);
        if (!std.math.isFinite(structural_carbon_g_c) or
            !std.math.isFinite(mobile_carbon_g_c))
            return error.NonFiniteSeasonalStorageAdmissionResult;
    }
    if (!std.math.isFinite(state.time_since_germination_h.*) or
        state.time_since_germination_h.* < 0)
        return error.InvalidSeasonalStorageAdmissionCounter;

    if (!state.remobilization_counter_initialized.*) {
        state.time_since_germination_h.* = 0;
        state.remobilization_counter_initialized.* = true;
    }
    return .{
        .structural_carbon_g_c = structural_carbon_g_c,
        .mobile_carbon_g_c = mobile_carbon_g_c,
    };
}

fn validateGateInputs(inputs: Inputs) !void {
    _ = try execution_calendar.fromDayOfYear(
        inputs.current_date.day_of_year,
        inputs.current_date.year,
    );
    if (inputs.planting_date) |date|
        _ = try execution_calendar.fromDayOfYear(date.day_of_year, date.year);
    inline for (.{
        inputs.current_branch_accumulated_leafoff_h,
        inputs.current_branch_required_leafoff_h,
        inputs.leafoff_remobilization_start_fraction,
        inputs.main_branch_accumulated_leafout_h,
        inputs.current_branch_required_leafout_h,
    }) |value| {
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSeasonalStorageAdmissionInput;
    }
    if (inputs.leafoff_remobilization_start_fraction > 1)
        return error.InvalidSeasonalStorageAdmissionInput;
}

fn testInputs() Inputs {
    return .{
        .growth_habit = .perennial,
        .lifecycle_initialized = true,
        .current_date = .{ .day_of_year = 150, .year = 2000 },
        .planting_date = null,
        .current_branch_accumulated_leafoff_h = 10,
        .current_branch_required_leafoff_h = 100,
        .leafoff_remobilization_start_fraction = 0.5,
        .main_branch_accumulated_leafout_h = 100,
        .current_branch_required_leafout_h = 100,
        .root_structural_carbon_g_c_by_layer = &.{ 1, 2, 3 },
        .root_mobile_carbon_g_c_by_layer = &.{ 0.5, 1, 1.5 },
    };
}

test "GROSUB main leafout gate totals runtime root carbon and resets latch" {
    var counter_h: f64 = 27;
    var initialized = false;
    const totals = (try apply(.{
        .time_since_germination_h = &counter_h,
        .remobilization_counter_initialized = &initialized,
    }, testInputs())).?;
    try std.testing.expectEqual(@as(f64, 6), totals.structural_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 3), totals.mobile_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), counter_h);
    try std.testing.expect(initialized);
}

test "annual initialization and planting-date gates are independent" {
    var counter_h: f64 = 5;
    var initialized = true;
    var inputs = testInputs();
    inputs.growth_habit = .annual;
    inputs.lifecycle_initialized = false;
    inputs.main_branch_accumulated_leafout_h = 0;
    try std.testing.expect((try apply(.{
        .time_since_germination_h = &counter_h,
        .remobilization_counter_initialized = &initialized,
    }, inputs)) != null);
    try std.testing.expectEqual(@as(f64, 5), counter_h);

    inputs.growth_habit = .perennial;
    inputs.lifecycle_initialized = true;
    inputs.planting_date = .{ .day_of_year = 150, .year = 2000 };
    try std.testing.expect((try apply(.{
        .time_since_germination_h = &counter_h,
        .remobilization_counter_initialized = &initialized,
    }, inputs)) != null);
}

test "leafoff threshold closes planting and leafout gates exactly at equality" {
    var counter_h: f64 = std.math.nan(f64);
    var initialized = false;
    var inputs = testInputs();
    inputs.planting_date = inputs.current_date;
    inputs.current_branch_accumulated_leafoff_h = 50;
    try std.testing.expectEqual(@as(?RootCarbonTotals, null), try apply(.{
        .time_since_germination_h = &counter_h,
        .remobilization_counter_initialized = &initialized,
    }, inputs));
    try std.testing.expect(!initialized);
}

test "runtime layer count has no source ceiling and preserves sum order" {
    const allocator = std.testing.allocator;
    const count = 47;
    const structural = try allocator.alloc(f64, count);
    defer allocator.free(structural);
    const mobile = try allocator.alloc(f64, count);
    defer allocator.free(mobile);
    for (structural, mobile, 0..) |*structural_g_c, *mobile_g_c, layer| {
        structural_g_c.* = @as(f64, @floatFromInt(layer + 1)) * 0.25;
        mobile_g_c.* = @as(f64, @floatFromInt(layer + 1)) * 0.125;
    }
    var inputs = testInputs();
    inputs.root_structural_carbon_g_c_by_layer = structural;
    inputs.root_mobile_carbon_g_c_by_layer = mobile;
    var counter_h: f64 = 1;
    var initialized = false;
    const totals = (try apply(.{
        .time_since_germination_h = &counter_h,
        .remobilization_counter_initialized = &initialized,
    }, inputs)).?;
    var expected_structural: f64 = 0;
    var expected_mobile: f64 = 0;
    for (structural, mobile) |structural_g_c, mobile_g_c| {
        expected_structural += structural_g_c;
        expected_mobile += mobile_g_c;
    }
    try std.testing.expectEqual(expected_structural, totals.structural_carbon_g_c);
    try std.testing.expectEqual(expected_mobile, totals.mobile_carbon_g_c);
}

test "invalid late root layer leaves initialization state atomic" {
    var counter_h: f64 = 9;
    var initialized = false;
    var inputs = testInputs();
    inputs.root_mobile_carbon_g_c_by_layer = &.{ 1, 2, std.math.nan(f64) };
    try std.testing.expectError(
        error.NonFiniteSeasonalStorageAdmissionState,
        apply(.{
            .time_since_germination_h = &counter_h,
            .remobilization_counter_initialized = &initialized,
        }, inputs),
    );
    try std.testing.expectEqual(@as(f64, 9), counter_h);
    try std.testing.expect(!initialized);
}
