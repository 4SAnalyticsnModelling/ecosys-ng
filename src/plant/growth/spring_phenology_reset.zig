const std = @import("std");

pub const GrowthHabit = enum { annual, perennial };
pub const Phenology = enum { evergreen, cold_deciduous, drought_deciduous, cold_and_drought_deciduous };
pub const EnableStatus = enum { enabled, disabled };

pub const State = struct {
    floral_initiation_group: []f64,
    current_phenological_node: []const f64,
    phenological_node_at_initiation: []f64,
    phenological_node_at_flowering: []f64,
    leaf_number_at_initiation: []f64,
    total_vegetative_stage_change: []f64,
    total_reproductive_stage_change: []f64,
    event_day: [][10]i32,
    plant_water_stress_accumulator: *f64,
};

pub const Inputs = struct {
    first_biological_substep: bool,
    first_branch: usize,
    end_branch: usize,
    main_branch: usize,
    growth_habit: GrowthHabit,
    phenology: Phenology,
    leafout_status: []const EnableStatus,
    leafoff_status: []const EnableStatus,
    accumulated_leafout_h: []const f64,
    required_leafout_h: []const f64,
    accumulated_leafoff_h: []const f64,
    required_leafoff_h: []const f64,
    initial_floral_initiation_group: f64,
    current_day: i32,
};

/// grosub.f lines 4182--4222. Applies outer first-substep and seasonal-completion
/// gates, then the perennial spring leafout reset. IDAY(1) receives current day,
/// IDAY(2..10) clear, and only the explicit main branch resets plant WSTR.
pub fn apply(state: State, inputs: Inputs) !usize {
    if (!inputs.first_biological_substep) return 0;
    const branch_count = state.floral_initiation_group.len;
    try validateDimensions(state, inputs, branch_count);
    if (!std.math.isFinite(inputs.initial_floral_initiation_group) or inputs.initial_floral_initiation_group < 0)
        return error.InvalidSpringPhenologyResetInput;
    if (inputs.first_branch > inputs.end_branch or inputs.end_branch > branch_count or inputs.main_branch >= branch_count)
        return error.SpringPhenologyBranchRangeOutOfBounds;

    const growth_outer_gate = inputs.growth_habit == .perennial or
        (inputs.growth_habit == .annual and
            (inputs.phenology == .drought_deciduous or inputs.phenology == .cold_and_drought_deciduous));
    if (!growth_outer_gate) return 0;

    // Validate the full selected range before any reset, so a late invalid
    // branch cannot leave earlier phenology partially reset.
    for (inputs.first_branch..inputs.end_branch) |branch| {
        inline for (.{ inputs.accumulated_leafout_h[branch], inputs.required_leafout_h[branch], inputs.accumulated_leafoff_h[branch], inputs.required_leafoff_h[branch], state.current_phenological_node[branch] }) |value|
            if (!std.math.isFinite(value) or value < 0) return error.InvalidSpringPhenologyResetState;
    }

    var reset_count: usize = 0;
    for (inputs.first_branch..inputs.end_branch) |branch| {
        const seasonal_trigger =
            (inputs.leafout_status[branch] == .enabled and inputs.accumulated_leafout_h[branch] >= inputs.required_leafout_h[branch]) or
            (inputs.leafoff_status[branch] == .enabled and inputs.accumulated_leafoff_h[branch] >= inputs.required_leafoff_h[branch]);
        if (!seasonal_trigger) continue;
        const spring_trigger = inputs.leafout_status[branch] == .enabled and
            inputs.growth_habit == .perennial and
            inputs.accumulated_leafout_h[branch] >= inputs.required_leafout_h[branch];
        if (!spring_trigger) continue;

        state.floral_initiation_group[branch] = inputs.initial_floral_initiation_group;
        state.phenological_node_at_initiation[branch] = state.current_phenological_node[branch];
        state.phenological_node_at_flowering[branch] = 0;
        state.leaf_number_at_initiation[branch] = 0;
        state.total_vegetative_stage_change[branch] = 0;
        state.total_reproductive_stage_change[branch] = 0;
        state.event_day[branch][0] = inputs.current_day;
        for (1..10) |event| state.event_day[branch][event] = 0;
        if (branch == inputs.main_branch) state.plant_water_stress_accumulator.* = 0;
        reset_count += 1;
    }
    return reset_count;
}

fn validateDimensions(state: State, inputs: Inputs, count: usize) !void {
    inline for (.{ state.current_phenological_node, state.phenological_node_at_initiation, state.phenological_node_at_flowering, state.leaf_number_at_initiation, state.total_vegetative_stage_change, state.total_reproductive_stage_change }) |values|
        if (values.len != count) return error.SpringPhenologyResetDimensionMismatch;
    if (state.event_day.len != count) return error.SpringPhenologyResetDimensionMismatch;
    inline for (.{ inputs.leafout_status, inputs.leafoff_status, inputs.accumulated_leafout_h, inputs.required_leafout_h, inputs.accumulated_leafoff_h, inputs.required_leafoff_h }) |values|
        if (values.len != count) return error.SpringPhenologyResetDimensionMismatch;
}

fn makeState(storage: *[7][3]f64, dates: *[3][10]i32, water: *f64) State {
    return .{ .floral_initiation_group = &storage[0], .current_phenological_node = &storage[1], .phenological_node_at_initiation = &storage[2], .phenological_node_at_flowering = &storage[3], .leaf_number_at_initiation = &storage[4], .total_vegetative_stage_change = &storage[5], .total_reproductive_stage_change = &storage[6], .event_day = dates, .plant_water_stress_accumulator = water };
}

test "perennial completed leafout resets runtime branches and explicit main water stress" {
    var s: [7][3]f64 = @splat(@splat(9));
    s[1] = .{ 2, 3, 4 };
    var dates: [3][10]i32 = @splat(@splat(7));
    var water: f64 = 5;
    const count = try apply(makeState(&s, &dates, &water), .{ .first_biological_substep = true, .first_branch = 0, .end_branch = 3, .main_branch = 2, .growth_habit = .perennial, .phenology = .cold_deciduous, .leafout_status = &.{ .enabled, .disabled, .enabled }, .leafoff_status = &.{ .disabled, .disabled, .disabled }, .accumulated_leafout_h = &.{ 10, 10, 10 }, .required_leafout_h = &.{ 10, 10, 10 }, .accumulated_leafoff_h = &.{ 0, 0, 0 }, .required_leafoff_h = &.{ 10, 10, 10 }, .initial_floral_initiation_group = 6, .current_day = 123 });
    try std.testing.expectEqual(@as(usize, 2), count);
    try std.testing.expectEqual(@as(f64, 6), s[0][2]);
    try std.testing.expectEqual(@as(f64, 4), s[2][2]);
    try std.testing.expectEqual(@as(i32, 123), dates[2][0]);
    try std.testing.expectEqual(@as(i32, 0), dates[2][9]);
    try std.testing.expectEqual(@as(f64, 0), water);
}
test "leafoff-only outer trigger does not execute spring reset" {
    var s: [7][3]f64 = @splat(@splat(9));
    var dates: [3][10]i32 = @splat(@splat(7));
    var water: f64 = 5;
    const count = try apply(makeState(&s, &dates, &water), .{ .first_biological_substep = true, .first_branch = 0, .end_branch = 1, .main_branch = 0, .growth_habit = .perennial, .phenology = .cold_deciduous, .leafout_status = &.{ .disabled, .disabled, .disabled }, .leafoff_status = &.{ .enabled, .disabled, .disabled }, .accumulated_leafout_h = &.{ 0, 0, 0 }, .required_leafout_h = &.{ 10, 10, 10 }, .accumulated_leafoff_h = &.{ 10, 0, 0 }, .required_leafoff_h = &.{ 10, 10, 10 }, .initial_floral_initiation_group = 6, .current_day = 123 });
    try std.testing.expectEqual(@as(usize, 0), count);
    try std.testing.expectEqual(@as(f64, 9), s[0][0]);
}
test "non-first biological substep is strict no-op before topology reads" {
    var x: [0]f64 = .{};
    var dates: [0][10]i32 = .{};
    var water: f64 = std.math.nan(f64);
    const state: State = .{ .floral_initiation_group = &x, .current_phenological_node = &x, .phenological_node_at_initiation = &x, .phenological_node_at_flowering = &x, .leaf_number_at_initiation = &x, .total_vegetative_stage_change = &x, .total_reproductive_stage_change = &x, .event_day = &dates, .plant_water_stress_accumulator = &water };
    try std.testing.expectEqual(@as(usize, 0), try apply(state, undefined));
}
test "late invalid branch leaves earlier state unchanged" {
    var s: [7][3]f64 = @splat(@splat(9));
    s[1][2] = std.math.nan(f64);
    var dates: [3][10]i32 = @splat(@splat(7));
    var water: f64 = 5;
    try std.testing.expectError(error.InvalidSpringPhenologyResetState, apply(makeState(&s, &dates, &water), .{ .first_biological_substep = true, .first_branch = 0, .end_branch = 3, .main_branch = 2, .growth_habit = .perennial, .phenology = .cold_deciduous, .leafout_status = &.{ .enabled, .enabled, .enabled }, .leafoff_status = &.{ .disabled, .disabled, .disabled }, .accumulated_leafout_h = &.{ 10, 10, 10 }, .required_leafout_h = &.{ 10, 10, 10 }, .accumulated_leafoff_h = &.{ 0, 0, 0 }, .required_leafoff_h = &.{ 10, 10, 10 }, .initial_floral_initiation_group = 6, .current_day = 123 }));
    try std.testing.expectEqual(@as(f64, 9), s[0][0]);
    try std.testing.expectEqual(@as(i32, 7), dates[0][0]);
    try std.testing.expectEqual(@as(f64, 5), water);
}
