const std = @import("std");
const growth_stages = @import("plant_growth_stages.zig");

pub const TimeIncrementWorkspace = struct {
    initialized: bool = false,
    increment_h: f64 = 0,
};

pub const State = struct {
    time_since_germination_h_by_branch: []f64,
    time_increment_workspace: *TimeIncrementWorkspace,
    seasonal_storage_carbon_g_c: *f64,
    branch_mobile_carbon_g_c: []f64,
    root_structural_carbon_g_c_by_layer: []const f64,
    root_mobile_carbon_g_c_by_layer: []f64,
};

pub const Inputs = struct {
    current_branch: usize,
    main_branch: usize,
    photoperiod_type: growth_stages.PhotoperiodType,
    phenology: growth_stages.PhenologyType,
    growth_habit: growth_stages.GrowthHabit,
    critical_photoperiod_h: f64,
    photoperiod_induction_difference_h: f64,
    daylength_h: f64,
    canopy_growth_temperature_response: f64,
    canopy_growth_water_fraction: f64,
    biological_timestep_h: f64,
    remobilization_duration_h: f64,
    storage_carbon_oxidation_fraction_per_h: f64,
    shoot_partition_fraction: f64,
    root_partition_fraction: f64,
    carbon_presence_threshold_g_c: f64,
    planting_layer: usize,
};

pub const Result = struct {
    remobilization_time_increment_h: f64,
    shoot_and_root_mobile_carbon_g_c: f64,
    oxidized_storage_carbon_g_c: f64,
    used_structural_root_distribution: bool,
};

/// Exact GROSUB lines 4676--4752. The main branch forms DATRP and advances
/// ATRP; subsequent branches consume that retained hourly workspace value.
/// Storage and mobile pools are g C, time is h, and root layers are runtime
/// sized. The caller owns the outer seasonal-storage admission gate.
pub fn apply(state: State, inputs: Inputs) !Result {
    try validateDimensions(state, inputs);
    const main_branch = inputs.current_branch == inputs.main_branch;
    var time_increment_h = state.time_increment_workspace.increment_h;
    var next_time_since_germination_h =
        state.time_since_germination_h_by_branch[inputs.current_branch];
    if (main_branch) {
        try validateMainTimeInputs(inputs);
        const cold_sensitive_long_day = inputs.photoperiod_type == .long_day and
            (inputs.phenology == .winter_deciduous or
                inputs.phenology == .winter_and_drought_deciduous);
        const photoperiod_response = if (cold_sensitive_long_day) response: {
            const photoperiod_deficit_h = @max(0.0, inputs.critical_photoperiod_h -
                inputs.photoperiod_induction_difference_h - inputs.daylength_h);
            break :response @exp(-0.0 * photoperiod_deficit_h);
        } else 1.0;
        time_increment_h = photoperiod_response *
            inputs.canopy_growth_temperature_response *
            inputs.canopy_growth_water_fraction * inputs.biological_timestep_h;
        next_time_since_germination_h += time_increment_h;
        if (!std.math.isFinite(time_increment_h) or time_increment_h < 0 or
            !std.math.isFinite(next_time_since_germination_h))
            return error.NonFiniteSeasonalStorageCarbonResult;
    } else if (!state.time_increment_workspace.initialized) {
        return error.UninitializedSeasonalStorageTimeIncrement;
    }
    if (!std.math.isFinite(time_increment_h) or time_increment_h < 0)
        return error.InvalidSeasonalStorageTimeIncrement;

    try validateCarbonState(state, inputs);
    var root_structural_total_g_c: f64 = 0;
    var root_mobile_total_g_c: f64 = 0;
    for (state.root_structural_carbon_g_c_by_layer, state.root_mobile_carbon_g_c_by_layer) |structural, mobile| {
        root_structural_total_g_c += structural;
        root_mobile_total_g_c += mobile;
        if (!std.math.isFinite(root_structural_total_g_c) or
            !std.math.isFinite(root_mobile_total_g_c))
            return error.NonFiniteSeasonalStorageCarbonResult;
    }
    const shoot_and_root_mobile_carbon_g_c = root_mobile_total_g_c +
        state.branch_mobile_carbon_g_c[inputs.current_branch];
    if (!std.math.isFinite(shoot_and_root_mobile_carbon_g_c))
        return error.NonFiniteSeasonalStorageCarbonResult;

    const annual_cold_continuation = inputs.growth_habit == .annual and
        (inputs.phenology == .evergreen or inputs.phenology == .winter_deciduous);
    const duration_admitted = next_time_since_germination_h <=
        inputs.remobilization_duration_h or annual_cold_continuation;
    var oxidized_storage_carbon_g_c: f64 = 0;
    if (duration_admitted and
        state.seasonal_storage_carbon_g_c.* > inputs.carbon_presence_threshold_g_c)
    {
        const oxidation_fraction =
            inputs.storage_carbon_oxidation_fraction_per_h * time_increment_h;
        oxidized_storage_carbon_g_c = @max(0.0, oxidation_fraction * state.seasonal_storage_carbon_g_c.*);
    }
    const next_storage_carbon_g_c =
        state.seasonal_storage_carbon_g_c.* - oxidized_storage_carbon_g_c;
    const shoot_transfer_g_c =
        oxidized_storage_carbon_g_c * inputs.shoot_partition_fraction;
    if (!std.math.isFinite(next_storage_carbon_g_c) or
        next_storage_carbon_g_c < 0)
        return error.SeasonalStorageCarbonOverdraw;
    const next_shoot_mobile_g_c =
        state.branch_mobile_carbon_g_c[inputs.current_branch] + shoot_transfer_g_c;
    if (!std.math.isFinite(next_shoot_mobile_g_c) or next_shoot_mobile_g_c < 0)
        return error.NonFiniteSeasonalStorageCarbonResult;

    const distributed = root_structural_total_g_c > inputs.carbon_presence_threshold_g_c and
        root_mobile_total_g_c > inputs.carbon_presence_threshold_g_c;
    for (state.root_mobile_carbon_g_c_by_layer, state.root_structural_carbon_g_c_by_layer, 0..) |mobile, structural, layer| {
        const transfer = if (distributed)
            structural / root_structural_total_g_c * oxidized_storage_carbon_g_c *
                inputs.root_partition_fraction
        else if (layer == inputs.planting_layer)
            oxidized_storage_carbon_g_c * inputs.root_partition_fraction
        else
            0;
        const next_mobile_g_c = mobile + transfer;
        if (!std.math.isFinite(next_mobile_g_c) or next_mobile_g_c < 0)
            return error.NonFiniteSeasonalStorageCarbonResult;
    }

    if (main_branch) {
        state.time_increment_workspace.increment_h = time_increment_h;
        state.time_increment_workspace.initialized = true;
        state.time_since_germination_h_by_branch[inputs.current_branch] =
            next_time_since_germination_h;
    }
    state.seasonal_storage_carbon_g_c.* = next_storage_carbon_g_c;
    state.branch_mobile_carbon_g_c[inputs.current_branch] = next_shoot_mobile_g_c;
    for (state.root_mobile_carbon_g_c_by_layer, state.root_structural_carbon_g_c_by_layer, 0..) |*mobile, structural, layer| {
        if (distributed) {
            mobile.* += structural / root_structural_total_g_c *
                oxidized_storage_carbon_g_c * inputs.root_partition_fraction;
        } else if (layer == inputs.planting_layer) {
            mobile.* += oxidized_storage_carbon_g_c * inputs.root_partition_fraction;
        }
    }
    return .{
        .remobilization_time_increment_h = time_increment_h,
        .shoot_and_root_mobile_carbon_g_c = shoot_and_root_mobile_carbon_g_c,
        .oxidized_storage_carbon_g_c = oxidized_storage_carbon_g_c,
        .used_structural_root_distribution = distributed,
    };
}

fn validateDimensions(state: State, inputs: Inputs) !void {
    const branch_count = state.time_since_germination_h_by_branch.len;
    if (branch_count == 0 or state.branch_mobile_carbon_g_c.len != branch_count or
        inputs.current_branch >= branch_count or inputs.main_branch >= branch_count or
        state.root_structural_carbon_g_c_by_layer.len == 0 or
        state.root_structural_carbon_g_c_by_layer.len !=
            state.root_mobile_carbon_g_c_by_layer.len or
        inputs.planting_layer >= state.root_mobile_carbon_g_c_by_layer.len)
        return error.SeasonalStorageCarbonDimensionMismatch;
}

fn validateMainTimeInputs(inputs: Inputs) !void {
    inline for (.{
        inputs.critical_photoperiod_h,
        inputs.photoperiod_induction_difference_h,
        inputs.daylength_h,
        inputs.canopy_growth_temperature_response,
        inputs.canopy_growth_water_fraction,
        inputs.biological_timestep_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSeasonalStorageTimeInput;
    if (inputs.daylength_h > 24 or inputs.biological_timestep_h == 0)
        return error.InvalidSeasonalStorageTimeInput;
}

fn validateCarbonState(state: State, inputs: Inputs) !void {
    inline for (.{
        inputs.remobilization_duration_h,
        inputs.storage_carbon_oxidation_fraction_per_h,
        inputs.shoot_partition_fraction,
        inputs.root_partition_fraction,
        inputs.carbon_presence_threshold_g_c,
        state.time_since_germination_h_by_branch[inputs.current_branch],
        state.seasonal_storage_carbon_g_c.*,
        state.branch_mobile_carbon_g_c[inputs.current_branch],
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSeasonalStorageCarbonState;
    if (inputs.remobilization_duration_h == 0 or
        @abs(inputs.shoot_partition_fraction + inputs.root_partition_fraction - 1) > 1e-12)
        return error.InvalidSeasonalStorageCarbonPartition;
    for (state.root_structural_carbon_g_c_by_layer, state.root_mobile_carbon_g_c_by_layer) |structural, mobile|
        inline for (.{ structural, mobile }) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidSeasonalStorageCarbonState;
}

const Fixture = struct {
    times_h: [3]f64 = .{ 1, 2, 3 },
    workspace: TimeIncrementWorkspace = .{},
    storage_g_c: f64 = 100,
    shoot_g_c: [3]f64 = .{ 10, 20, 30 },
    structural_g_c: [3]f64 = .{ 1, 2, 3 },
    root_mobile_g_c: [3]f64 = .{ 3, 2, 1 },

    fn state(self: *Fixture) State {
        return .{
            .time_since_germination_h_by_branch = &self.times_h,
            .time_increment_workspace = &self.workspace,
            .seasonal_storage_carbon_g_c = &self.storage_g_c,
            .branch_mobile_carbon_g_c = &self.shoot_g_c,
            .root_structural_carbon_g_c_by_layer = &self.structural_g_c,
            .root_mobile_carbon_g_c_by_layer = &self.root_mobile_g_c,
        };
    }
};

fn testInputs(branch: usize, main: usize) Inputs {
    return .{
        .current_branch = branch,
        .main_branch = main,
        .photoperiod_type = .long_day,
        .phenology = .winter_deciduous,
        .growth_habit = .perennial,
        .critical_photoperiod_h = 14,
        .photoperiod_induction_difference_h = 1,
        .daylength_h = 10,
        .canopy_growth_temperature_response = 0.5,
        .canopy_growth_water_fraction = 0.4,
        .biological_timestep_h = 1,
        .remobilization_duration_h = 100,
        .storage_carbon_oxidation_fraction_per_h = 0.1,
        .shoot_partition_fraction = 0.25,
        .root_partition_fraction = 0.75,
        .carbon_presence_threshold_g_c = 1e-12,
        .planting_layer = 1,
    };
}

fn sum(values: []const f64) f64 {
    var total: f64 = 0;
    for (values) |value| total += value;
    return total;
}

test "GROSUB main DATRP and distributed root transfer conserve carbon" {
    var fixture: Fixture = .{};
    const before = fixture.storage_g_c + fixture.shoot_g_c[1] + sum(&fixture.root_mobile_g_c);
    const result = try apply(fixture.state(), testInputs(1, 1));
    try std.testing.expectEqual(@as(f64, 0.2), result.remobilization_time_increment_h);
    try std.testing.expectEqual(@as(f64, 2.2), fixture.times_h[1]);
    try std.testing.expect(fixture.workspace.initialized);
    try std.testing.expect(result.used_structural_root_distribution);
    const after = fixture.storage_g_c + fixture.shoot_g_c[1] + sum(&fixture.root_mobile_g_c);
    try std.testing.expectApproxEqAbs(before, after, 1e-13);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.oxidized_storage_carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 98), fixture.storage_g_c, 1e-15);
}

test "source zero photoperiod coefficient leaves DATRP unmodified" {
    var fixture: Fixture = .{};
    var inputs = testInputs(1, 1);
    inputs.critical_photoperiod_h = 24;
    inputs.daylength_h = 1;
    const result = try apply(fixture.state(), inputs);
    try std.testing.expectEqual(@as(f64, 0.2), result.remobilization_time_increment_h);
}

test "later branch consumes retained DATRP without advancing its counter" {
    var fixture: Fixture = .{};
    _ = try apply(fixture.state(), testInputs(1, 1));
    const before_time_h = fixture.times_h[2];
    const result = try apply(fixture.state(), testInputs(2, 1));
    try std.testing.expectEqual(@as(f64, 0.2), result.remobilization_time_increment_h);
    try std.testing.expectEqual(before_time_h, fixture.times_h[2]);
}

test "root fallback sends carbon only to runtime planting layer" {
    var fixture: Fixture = .{};
    fixture.structural_g_c = .{ 0, 0, 0 };
    const before = fixture.root_mobile_g_c;
    const result = try apply(fixture.state(), testInputs(1, 1));
    try std.testing.expect(!result.used_structural_root_distribution);
    try std.testing.expectEqual(before[0], fixture.root_mobile_g_c[0]);
    try std.testing.expect(fixture.root_mobile_g_c[1] > before[1]);
    try std.testing.expectEqual(before[2], fixture.root_mobile_g_c[2]);
}

test "duration gate and annual cold continuation preserve source inequality" {
    var fixture: Fixture = .{};
    fixture.times_h[1] = 101;
    var inputs = testInputs(1, 1);
    const stopped = try apply(fixture.state(), inputs);
    try std.testing.expectEqual(@as(f64, 0), stopped.oxidized_storage_carbon_g_c);

    fixture = .{};
    fixture.times_h[1] = 101;
    inputs.growth_habit = .annual;
    const continued = try apply(fixture.state(), inputs);
    try std.testing.expect(continued.oxidized_storage_carbon_g_c > 0);
}

test "uninitialized retained DATRP and late invalid root fail atomically" {
    var fixture: Fixture = .{};
    try std.testing.expectError(
        error.UninitializedSeasonalStorageTimeIncrement,
        apply(fixture.state(), testInputs(2, 1)),
    );
    fixture.root_mobile_g_c[2] = std.math.nan(f64);
    const before = fixture;
    try std.testing.expectError(
        error.InvalidSeasonalStorageCarbonState,
        apply(fixture.state(), testInputs(1, 1)),
    );
    try std.testing.expectEqualSlices(u8, std.mem.asBytes(&before), std.mem.asBytes(&fixture));
}
