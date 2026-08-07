const std = @import("std");
const growth_stages = @import("../lifecycle/growth_stages.zig");

pub const BranchState = struct {
    /// FLG4: consecutive biological hours without detectable grain fill (h).
    hours_without_grain_fill_h: f64,
    /// IDAY(10): first day of year on which physiological maturity was reached.
    physiological_maturity_day_of_year: ?u16,
    /// VRNF: hours accumulated toward leafoff (h).
    accumulated_leafoff_h: f64,
};

pub const Inputs = struct {
    /// IDAY(8) != 0: final seed-number setting has ended for this branch.
    final_seed_number_set_by_branch: []const bool,
    /// XLOCC, indexed by branch (g C biological-timestep^-1).
    grain_carbon_translocation_g_c_by_branch: []const f64,
    /// VRNX, indexed by branch (h).
    required_leafoff_h_by_branch: []const f64,
    grain_fill_detection_threshold_g_c: f64,
    biological_timestep_h: f64,
    physiological_maturity_no_fill_h: f64,
    annual_leafoff_delay_h: f64,
    current_day_of_year: u16,
    growth_habit: growth_stages.GrowthHabit,
    phenology: growth_stages.PhenologyType,
};

pub const Result = struct {
    active_branch_count: usize,
    maturity_dates_set: usize,
    annual_leafoff_forced: usize,
};

/// Exact grosub.f lines 5242--5276 no-grain-fill maturity and annual leafoff
/// transition in ascending NB order. Array dimensions are runtime branch
/// counts. A full active-branch preflight makes the sweep atomic.
pub fn apply(branches: []BranchState, inputs: Inputs) !Result {
    if (branches.len == 0 or
        inputs.final_seed_number_set_by_branch.len != branches.len or
        inputs.grain_carbon_translocation_g_c_by_branch.len != branches.len or
        inputs.required_leafoff_h_by_branch.len != branches.len)
        return error.GrainFillMaturityDimensionMismatch;

    inline for (.{
        inputs.grain_fill_detection_threshold_g_c,
        inputs.biological_timestep_h,
        inputs.physiological_maturity_no_fill_h,
        inputs.annual_leafoff_delay_h,
    }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteGrainFillMaturityInput;
        if (value < 0) return error.NegativeGrainFillMaturityInput;
    }
    if (inputs.biological_timestep_h == 0 or
        inputs.physiological_maturity_no_fill_h == 0 or
        inputs.current_day_of_year == 0 or inputs.current_day_of_year > 366)
        return error.InvalidGrainFillMaturityInput;

    for (branches, 0..) |branch, index| {
        if (!inputs.final_seed_number_set_by_branch[index]) continue;
        const translocation_g_c = inputs.grain_carbon_translocation_g_c_by_branch[index];
        if (!std.math.isFinite(translocation_g_c) or
            !std.math.isFinite(branch.hours_without_grain_fill_h))
            return error.NonFiniteGrainFillMaturityState;
        if (translocation_g_c < 0 or branch.hours_without_grain_fill_h < 0)
            return error.NegativeGrainFillMaturityState;
        if (branch.physiological_maturity_day_of_year) |day| {
            if (day == 0 or day > 366) return error.InvalidPhysiologicalMaturityDay;
        }

        const next_no_fill_h = if (translocation_g_c <= inputs.grain_fill_detection_threshold_g_c)
            branch.hours_without_grain_fill_h + inputs.biological_timestep_h
        else
            0;
        if (!std.math.isFinite(next_no_fill_h)) return error.NonFiniteGrainFillMaturityState;

        if (inputs.growth_habit == .annual and inputs.phenology != .evergreen and
            next_no_fill_h > inputs.physiological_maturity_no_fill_h + inputs.annual_leafoff_delay_h)
        {
            const required_h = inputs.required_leafoff_h_by_branch[index];
            if (!std.math.isFinite(required_h)) return error.NonFiniteGrainFillMaturityInput;
            if (required_h < 0 or !std.math.isFinite(required_h + 0.5))
                return error.InvalidGrainFillMaturityInput;
        }
    }

    var result = Result{ .active_branch_count = 0, .maturity_dates_set = 0, .annual_leafoff_forced = 0 };
    for (branches, 0..) |*branch, index| {
        if (!inputs.final_seed_number_set_by_branch[index]) continue;
        result.active_branch_count += 1;
        branch.hours_without_grain_fill_h = if (inputs.grain_carbon_translocation_g_c_by_branch[index] <= inputs.grain_fill_detection_threshold_g_c)
            branch.hours_without_grain_fill_h + inputs.biological_timestep_h
        else
            0;
        if (branch.hours_without_grain_fill_h >= inputs.physiological_maturity_no_fill_h and
            branch.physiological_maturity_day_of_year == null)
        {
            branch.physiological_maturity_day_of_year = inputs.current_day_of_year;
            result.maturity_dates_set += 1;
        }
        if (inputs.growth_habit == .annual and inputs.phenology != .evergreen and
            branch.hours_without_grain_fill_h > inputs.physiological_maturity_no_fill_h + inputs.annual_leafoff_delay_h)
        {
            branch.accumulated_leafoff_h = inputs.required_leafoff_h_by_branch[index] + 0.5;
            result.annual_leafoff_forced += 1;
        }
    }
    return result;
}

fn testInputs(final_set: []const bool, translocation: []const f64, required: []const f64) Inputs {
    return .{
        .final_seed_number_set_by_branch = final_set,
        .grain_carbon_translocation_g_c_by_branch = translocation,
        .required_leafoff_h_by_branch = required,
        .grain_fill_detection_threshold_g_c = 0.01,
        .biological_timestep_h = 0.5,
        .physiological_maturity_no_fill_h = 2,
        .annual_leafoff_delay_h = 1,
        .current_day_of_year = 240,
        .growth_habit = .annual,
        .phenology = .winter_deciduous,
    };
}

test "GROSUB threshold equality increments and maturity equality sets first date" {
    var branches = [_]BranchState{.{ .hours_without_grain_fill_h = 1.5, .physiological_maturity_day_of_year = null, .accumulated_leafoff_h = 0 }};
    const result = try apply(&branches, testInputs(&.{true}, &.{0.01}, &.{10}));
    try std.testing.expectEqual(@as(f64, 2), branches[0].hours_without_grain_fill_h);
    try std.testing.expectEqual(@as(?u16, 240), branches[0].physiological_maturity_day_of_year);
    try std.testing.expectEqual(@as(usize, 1), result.maturity_dates_set);
    try std.testing.expectEqual(@as(f64, 0), branches[0].accumulated_leafoff_h);
}

test "GROSUB detectable fill resets counter and preserves first maturity date" {
    var branches = [_]BranchState{.{ .hours_without_grain_fill_h = 9, .physiological_maturity_day_of_year = 190, .accumulated_leafoff_h = 4 }};
    _ = try apply(&branches, testInputs(&.{true}, &.{0.011}, &.{10}));
    try std.testing.expectEqual(@as(f64, 0), branches[0].hours_without_grain_fill_h);
    try std.testing.expectEqual(@as(?u16, 190), branches[0].physiological_maturity_day_of_year);
    try std.testing.expectEqual(@as(f64, 4), branches[0].accumulated_leafoff_h);
}

test "GROSUB annual leafoff uses strict greater than threshold plus delay" {
    var branches = [_]BranchState{
        .{ .hours_without_grain_fill_h = 2.5, .physiological_maturity_day_of_year = null, .accumulated_leafoff_h = 1 },
        .{ .hours_without_grain_fill_h = 3.0, .physiological_maturity_day_of_year = null, .accumulated_leafoff_h = 1 },
    };
    const result = try apply(&branches, testInputs(&.{ true, true }, &.{ 0, 0 }, &.{ 8, 9 }));
    try std.testing.expectEqual(@as(f64, 1), branches[0].accumulated_leafoff_h);
    try std.testing.expectEqual(@as(f64, 9.5), branches[1].accumulated_leafoff_h);
    try std.testing.expectEqual(@as(usize, 1), result.annual_leafoff_forced);
}

test "GROSUB inactive seed stage does not read branch scientific values" {
    var branches = [_]BranchState{.{ .hours_without_grain_fill_h = std.math.nan(f64), .physiological_maturity_day_of_year = null, .accumulated_leafoff_h = std.math.nan(f64) }};
    const result = try apply(&branches, testInputs(&.{false}, &.{std.math.nan(f64)}, &.{std.math.nan(f64)}));
    try std.testing.expectEqual(@as(usize, 0), result.active_branch_count);
    try std.testing.expect(std.math.isNan(branches[0].hours_without_grain_fill_h));
}

test "GROSUB sweep supports arbitrary runtime branch counts" {
    var branches: [41]BranchState = undefined;
    var final_set = [_]bool{true} ** 41;
    var translocation = [_]f64{0} ** 41;
    var required = [_]f64{12} ** 41;
    for (&branches) |*branch| branch.* = .{ .hours_without_grain_fill_h = 1.5, .physiological_maturity_day_of_year = null, .accumulated_leafoff_h = 0 };
    const result = try apply(&branches, testInputs(&final_set, &translocation, &required));
    try std.testing.expectEqual(@as(usize, 41), result.active_branch_count);
    try std.testing.expectEqual(@as(usize, 41), result.maturity_dates_set);
}

test "GROSUB invalid late active branch leaves the complete sweep unchanged" {
    var branches = [_]BranchState{
        .{ .hours_without_grain_fill_h = 1.5, .physiological_maturity_day_of_year = null, .accumulated_leafoff_h = 0 },
        .{ .hours_without_grain_fill_h = std.math.nan(f64), .physiological_maturity_day_of_year = null, .accumulated_leafoff_h = 0 },
    };
    try std.testing.expectError(error.NonFiniteGrainFillMaturityState, apply(&branches, testInputs(&.{ true, true }, &.{ 0, 0 }, &.{ 4, 4 })));
    try std.testing.expectEqual(@as(f64, 1.5), branches[0].hours_without_grain_fill_h);
    try std.testing.expectEqual(@as(?u16, null), branches[0].physiological_maturity_day_of_year);
    try std.testing.expectEqual(@as(f64, 0), branches[0].accumulated_leafoff_h);
    try std.testing.expect(std.math.isNan(branches[1].hours_without_grain_fill_h));
}
