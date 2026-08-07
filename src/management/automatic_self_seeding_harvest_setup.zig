const std = @import("std");
const execution_calendar = @import("../driver/execution_calendar_date.zig");
const growth_stages = @import("../plant/lifecycle/growth_stages.zig");

pub const HarvestKind = enum { none, grain, aboveground, pruning, animal_grazing, insect_grazing };
pub const Termination = enum { retain, terminate, terminate_and_reseed };

pub const RemovalFractions = struct {
    leaf: f64,
    nonfoliar: f64,
    woody: f64,
    standing_dead: f64,
};

pub const SimulationDate = struct {
    day_of_year: u16,
    year: u16,
};

pub const State = struct {
    harvest_date: SimulationDate,
    harvest_kind: HarvestKind,
    termination: Termination,
    cutting_height_m_or_leaf_area_fraction: f64,
    thinning_fraction_or_consumption_rate: f64,
    harvested_fraction: RemovalFractions,
    ecosystem_export_fraction: RemovalFractions,
    planting_date: ?SimulationDate,
    plant_initialization_requested: bool,
};

pub const Inputs = struct {
    branch_count: usize,
    current_branch: usize,
    main_branch: usize,
    current_day_of_year: u16,
    current_year: u16,
    growth_habit: growth_stages.GrowthHabit,
    phenology: growth_stages.PhenologyType,
};

/// Exact grosub.f lines 4572--4605 automatic harvest setup. The caller retains
/// the enclosing, already-admitted IFLGR transaction. Branch identity is
/// explicit: only NB == NB1 schedules the grain harvest. `planting_date=null`
/// is the typed ecosys-ng representation of the source -1E+06 sentinels.
pub fn apply(state: *State, inputs: Inputs) !bool {
    if (inputs.branch_count == 0 or inputs.current_branch >= inputs.branch_count or
        inputs.main_branch >= inputs.branch_count)
        return error.AutomaticSelfSeedingBranchIndexOutOfBounds;
    if (inputs.current_branch != inputs.main_branch) return false;
    if (inputs.growth_habit != .annual or inputs.phenology == .evergreen)
        return false;
    _ = try execution_calendar.fromDayOfYear(
        inputs.current_day_of_year,
        inputs.current_year,
    );

    state.harvest_date.day_of_year = inputs.current_day_of_year;
    state.harvest_date.year = inputs.current_year;
    state.harvest_kind = .grain;
    state.termination = .terminate_and_reseed;
    state.cutting_height_m_or_leaf_area_fraction = 0;
    state.thinning_fraction_or_consumption_rate = 0;
    state.harvested_fraction.leaf = 1;
    state.harvested_fraction.nonfoliar = 1;
    state.harvested_fraction.woody = 1;
    state.harvested_fraction.standing_dead = 1;
    state.ecosystem_export_fraction.leaf = 0;
    state.ecosystem_export_fraction.nonfoliar = 1;
    state.ecosystem_export_fraction.woody = 0;
    state.ecosystem_export_fraction.standing_dead = 0;
    state.planting_date = null;
    state.plant_initialization_requested = true;
    return true;
}

fn initialState() State {
    return .{
        .harvest_date = .{ .day_of_year = 200, .year = 1999 },
        .harvest_kind = .pruning,
        .termination = .retain,
        .cutting_height_m_or_leaf_area_fraction = 2,
        .thinning_fraction_or_consumption_rate = 0.5,
        .harvested_fraction = .{ .leaf = 0.1, .nonfoliar = 0.2, .woody = 0.3, .standing_dead = 0.4 },
        .ecosystem_export_fraction = .{ .leaf = 0.9, .nonfoliar = 0.8, .woody = 0.7, .standing_dead = 0.6 },
        .planting_date = .{ .day_of_year = 100, .year = 1998 },
        .plant_initialization_requested = false,
    };
}

fn request(branch_count: usize, current_branch: usize, main_branch: usize) Inputs {
    return .{
        .branch_count = branch_count,
        .current_branch = current_branch,
        .main_branch = main_branch,
        .current_day_of_year = 61,
        .current_year = 2000,
        .growth_habit = .annual,
        .phenology = .winter_deciduous,
    };
}

test "GROSUB main deciduous annual schedules exact grain reseed harvest" {
    var state = initialState();
    try std.testing.expect(try apply(&state, request(7, 4, 4)));
    try std.testing.expectEqual(SimulationDate{ .day_of_year = 61, .year = 2000 }, state.harvest_date);
    try std.testing.expectEqual(HarvestKind.grain, state.harvest_kind);
    try std.testing.expectEqual(Termination.terminate_and_reseed, state.termination);
    try std.testing.expectEqual(@as(f64, 0), state.cutting_height_m_or_leaf_area_fraction);
    try std.testing.expectEqual(@as(f64, 0), state.thinning_fraction_or_consumption_rate);
    try std.testing.expectEqual(RemovalFractions{ .leaf = 1, .nonfoliar = 1, .woody = 1, .standing_dead = 1 }, state.harvested_fraction);
    try std.testing.expectEqual(RemovalFractions{ .leaf = 0, .nonfoliar = 1, .woody = 0, .standing_dead = 0 }, state.ecosystem_export_fraction);
    try std.testing.expectEqual(@as(?SimulationDate, null), state.planting_date);
    try std.testing.expect(state.plant_initialization_requested);
}

test "explicit main identity selects one branch without a source ceiling" {
    const branch_count = 43;
    const main_branch = 37;
    var scheduled: usize = 0;
    var state = initialState();
    for (0..branch_count) |branch| {
        if (try apply(&state, request(branch_count, branch, main_branch))) {
            scheduled += 1;
        }
    }
    try std.testing.expectEqual(@as(usize, 1), scheduled);
    try std.testing.expectEqual(HarvestKind.grain, state.harvest_kind);
}

test "non-main perennial and evergreen gates do not read or mutate schedule" {
    var state = initialState();
    const before = state;
    var inputs = request(5, 2, 3);
    inputs.current_day_of_year = 0;
    try std.testing.expect(!try apply(&state, inputs));
    try std.testing.expectEqual(before, state);

    inputs = request(5, 3, 3);
    inputs.growth_habit = .perennial;
    inputs.current_day_of_year = 0;
    try std.testing.expect(!try apply(&state, inputs));
    try std.testing.expectEqual(before, state);

    inputs = request(5, 3, 3);
    inputs.phenology = .evergreen;
    inputs.current_day_of_year = 0;
    try std.testing.expect(!try apply(&state, inputs));
    try std.testing.expectEqual(before, state);
}

test "invalid admitted date leaves schedule unchanged" {
    var state = initialState();
    const before = state;
    var inputs = request(5, 3, 3);
    inputs.current_day_of_year = 366;
    inputs.current_year = 2001;
    try std.testing.expectError(
        error.InvalidExecutionCalendarDay,
        apply(&state, inputs),
    );
    try std.testing.expectEqual(before, state);
}

test "invalid explicit branch identity fails before schedule mutation" {
    var state = initialState();
    const before = state;
    try std.testing.expectError(
        error.AutomaticSelfSeedingBranchIndexOutOfBounds,
        apply(&state, request(4, 4, 2)),
    );
    try std.testing.expectError(
        error.AutomaticSelfSeedingBranchIndexOutOfBounds,
        apply(&state, request(4, 2, 4)),
    );
    try std.testing.expectEqual(before, state);
}
