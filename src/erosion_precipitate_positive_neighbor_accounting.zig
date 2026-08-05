const std = @import("std");

pub const DisturbanceMode = enum {
    no_profile_effects,
    freeze_thaw,
    freeze_thaw_and_erosion,
    freeze_thaw_and_organic_matter,
    freeze_thaw_erosion_and_organic_matter,
};

pub const TransportAxis = enum {
    east_west,
    north_south,
    vertical,
};

pub const NeighborConnection = enum {
    blocked,
    connected,
};

/// Precipitate and silicate pools in exact REDIST subtraction order.
///
/// The legacy `F` suffix on the second six silicate pools is not defined by
/// local science comments, so these retain a neutral secondary-pool name.
pub const MineralPool = enum {
    precipitated_aluminum_hydroxide,
    precipitated_iron_hydroxide,
    precipitated_calcium_carbonate,
    precipitated_calcium_sulfate,
    silicate_aluminum,
    silicate_iron,
    silicate_calcium,
    silicate_magnesium,
    silicate_sodium,
    silicate_potassium,
    secondary_silicate_aluminum,
    secondary_silicate_iron,
    secondary_silicate_calcium,
    secondary_silicate_magnesium,
    secondary_silicate_sodium,
    secondary_silicate_potassium,
    nonband_precipitated_aluminum_phosphate,
    nonband_precipitated_iron_phosphate,
    nonband_precipitated_calcium_hydrogen_phosphate,
    nonband_precipitated_apatite,
    nonband_precipitated_calcium_dihydrogen_phosphate,
    band_precipitated_aluminum_phosphate,
    band_precipitated_iron_phosphate,
    band_precipitated_calcium_hydrogen_phosphate,
    band_precipitated_apatite,
    band_precipitated_calcium_dihydrogen_phosphate,
};

pub const mineral_pool_count: usize =
    @typeInfo(MineralPool).@"enum".fields.len;

pub const NeighborSide = struct {
    local_total_sediment_megagrams_per_step: f64,
    positive_neighbor_total_sediment_megagrams_per_step: f64,
    /// Positive-neighbor `P*ER`, `Q*ER`, and `P*EB`, mol/step by pool.
    positive_neighbor_mol_per_step_by_pool: []const f64,
    connection: NeighborConnection,
};

pub const Inputs = struct {
    disturbance_mode: DisturbanceMode,
    transport_axis: TransportAxis,
    sediment_activity_threshold_megagrams_per_step: f64,
    neighbor_by_boundary_side: []const NeighborSide,
};

pub const State = struct {
    net_erosion_mol_per_step_by_pool: []f64,
};

pub const Workspace = struct {
    net_erosion_mol_per_step_by_pool: []f64,
};

/// Subtracts connected positive-neighbor precipitate and silicate erosion.
///
/// Traceability: REDIST.F lines 3123--3148, within gates 2888, 2894--2895,
/// and 3051. Runtime slices preserve all 26 molar inventories and their source
/// order. The strict activity gate is evaluated before the connection gate,
/// and state commits only after every finite subtraction succeeds.
pub fn account(inputs: Inputs, state: *State, workspace: Workspace) !void {
    if (!erosionEnabled(inputs.disturbance_mode) or
        inputs.transport_axis == .vertical)
    {
        return;
    }
    try validateInputs(inputs, state.*, workspace);
    @memcpy(
        workspace.net_erosion_mol_per_step_by_pool,
        state.net_erosion_mol_per_step_by_pool,
    );

    for (inputs.neighbor_by_boundary_side) |side| {
        if (@abs(side.local_total_sediment_megagrams_per_step) <=
            inputs.sediment_activity_threshold_megagrams_per_step and
            @abs(side.positive_neighbor_total_sediment_megagrams_per_step) <=
                inputs.sediment_activity_threshold_megagrams_per_step)
        {
            continue;
        }
        if (side.connection == .blocked) continue;
        for (
            workspace.net_erosion_mol_per_step_by_pool,
            side.positive_neighbor_mol_per_step_by_pool,
        ) |*candidate, contribution| {
            const result = candidate.* - contribution;
            if (!std.math.isFinite(result))
                return error.NonFinitePrecipitateNeighborErosionResult;
            candidate.* = result;
        }
    }

    @memcpy(
        state.net_erosion_mol_per_step_by_pool,
        workspace.net_erosion_mol_per_step_by_pool,
    );
}

fn erosionEnabled(mode: DisturbanceMode) bool {
    return mode == .freeze_thaw_and_erosion or
        mode == .freeze_thaw_erosion_and_organic_matter;
}

fn validateInputs(inputs: Inputs, state: State, workspace: Workspace) !void {
    if (inputs.neighbor_by_boundary_side.len == 0)
        return error.InvalidPrecipitateNeighborErosionDimensions;
    if (state.net_erosion_mol_per_step_by_pool.len != mineral_pool_count or
        workspace.net_erosion_mol_per_step_by_pool.len != mineral_pool_count)
    {
        return error.PrecipitateNeighborErosionDimensionMismatch;
    }
    if (!std.math.isFinite(inputs.sediment_activity_threshold_megagrams_per_step))
        return error.NonFinitePrecipitateNeighborErosionInput;
    if (inputs.sediment_activity_threshold_megagrams_per_step < 0)
        return error.InvalidErosionActivityThreshold;

    try validateFinite(state.net_erosion_mol_per_step_by_pool);
    const scratch = workspace.net_erosion_mol_per_step_by_pool;
    if (overlap(scratch, state.net_erosion_mol_per_step_by_pool))
        return error.PrecipitateNeighborErosionWorkspaceOverlap;
    for (inputs.neighbor_by_boundary_side) |side| {
        if (side.positive_neighbor_mol_per_step_by_pool.len !=
            mineral_pool_count)
        {
            return error.PrecipitateNeighborErosionDimensionMismatch;
        }
        if (!std.math.isFinite(side.local_total_sediment_megagrams_per_step) or
            !std.math.isFinite(
                side.positive_neighbor_total_sediment_megagrams_per_step,
            ))
        {
            return error.NonFinitePrecipitateNeighborErosionInput;
        }
        try validateFinite(side.positive_neighbor_mol_per_step_by_pool);
        if (overlap(scratch, side.positive_neighbor_mol_per_step_by_pool))
            return error.PrecipitateNeighborErosionWorkspaceOverlap;
    }
}

fn validateFinite(values: []const f64) !void {
    for (values) |value|
        if (!std.math.isFinite(value))
            return error.NonFinitePrecipitateNeighborErosionInput;
}

fn overlap(left: []const f64, right: []const f64) bool {
    if (left.len == 0 or right.len == 0) return false;
    const item_size = @sizeOf(f64);
    const left_start = @intFromPtr(left.ptr);
    const right_start = @intFromPtr(right.ptr);
    const left_bytes = std.math.mul(usize, left.len, item_size) catch return true;
    const right_bytes = std.math.mul(usize, right.len, item_size) catch return true;
    const left_end = std.math.add(usize, left_start, left_bytes) catch return true;
    const right_end = std.math.add(usize, right_start, right_bytes) catch return true;
    return left_start < right_end and right_start < left_end;
}

fn expectAll(values: []const f64, expected: f64) !void {
    for (values) |value| try std.testing.expectEqual(expected, value);
}

test "pool order preserves all 26 precipitate and silicate positions" {
    try std.testing.expectEqual(@as(usize, 26), mineral_pool_count);
    try std.testing.expectEqual(
        @as(usize, 10),
        @intFromEnum(MineralPool.secondary_silicate_aluminum),
    );
    try std.testing.expectEqual(
        @as(usize, 16),
        @intFromEnum(MineralPool.nonband_precipitated_aluminum_phosphate),
    );
    try std.testing.expectEqual(
        @as(usize, 25),
        @intFromEnum(MineralPool.band_precipitated_calcium_dihydrogen_phosphate),
    );
}

test "connected active sides subtract all precipitate inventories" {
    const three = [_]f64{3} ** mineral_pool_count;
    const five = [_]f64{5} ** mineral_pool_count;
    const sides = [_]NeighborSide{
        .{
            .local_total_sediment_megagrams_per_step = 2,
            .positive_neighbor_total_sediment_megagrams_per_step = 0,
            .positive_neighbor_mol_per_step_by_pool = &three,
            .connection = .connected,
        },
        .{
            .local_total_sediment_megagrams_per_step = 0,
            .positive_neighbor_total_sediment_megagrams_per_step = -4,
            .positive_neighbor_mol_per_step_by_pool = &five,
            .connection = .connected,
        },
    };
    var totals = [_]f64{100} ** mineral_pool_count;
    var scratch: [mineral_pool_count]f64 = undefined;
    var state = State{ .net_erosion_mol_per_step_by_pool = &totals };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    }, &state, .{ .net_erosion_mol_per_step_by_pool = &scratch });
    try expectAll(&totals, 92);
}

test "activity and connection gates remain independent" {
    const seven = [_]f64{7} ** mineral_pool_count;
    const four = [_]f64{4} ** mineral_pool_count;
    const sides = [_]NeighborSide{
        .{
            .local_total_sediment_megagrams_per_step = 2,
            .positive_neighbor_total_sediment_megagrams_per_step = 2,
            .positive_neighbor_mol_per_step_by_pool = &seven,
            .connection = .blocked,
        },
        .{
            .local_total_sediment_megagrams_per_step = 1,
            .positive_neighbor_total_sediment_megagrams_per_step = -1,
            .positive_neighbor_mol_per_step_by_pool = &seven,
            .connection = .connected,
        },
        .{
            .local_total_sediment_megagrams_per_step = 2,
            .positive_neighbor_total_sediment_megagrams_per_step = 1,
            .positive_neighbor_mol_per_step_by_pool = &four,
            .connection = .connected,
        },
    };
    var totals = [_]f64{10} ** mineral_pool_count;
    var scratch: [mineral_pool_count]f64 = undefined;
    var state = State{ .net_erosion_mol_per_step_by_pool = &totals };
    try account(.{
        .disturbance_mode = .freeze_thaw_erosion_and_organic_matter,
        .transport_axis = .north_south,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    }, &state, .{ .net_erosion_mol_per_step_by_pool = &scratch });
    try expectAll(&totals, 6);
}

test "shared face transfer conserves all precipitate inventories exactly" {
    var shared: [mineral_pool_count]f64 = undefined;
    for (&shared, 0..) |*value, index| value.* = @floatFromInt(index + 1);
    var source = shared;
    const sides = [_]NeighborSide{.{
        .local_total_sediment_megagrams_per_step = 0,
        .positive_neighbor_total_sediment_megagrams_per_step = 2,
        .positive_neighbor_mol_per_step_by_pool = &shared,
        .connection = .connected,
    }};
    var scratch: [mineral_pool_count]f64 = undefined;
    var state = State{ .net_erosion_mol_per_step_by_pool = &source };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    }, &state, .{ .net_erosion_mol_per_step_by_pool = &scratch });
    try expectAll(&source, 0);
}

test "inactive process and vertical axes bypass unused storage" {
    var state = State{ .net_erosion_mol_per_step_by_pool = &.{} };
    try account(.{
        .disturbance_mode = .freeze_thaw,
        .transport_axis = .east_west,
        .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
        .neighbor_by_boundary_side = &.{},
    }, &state, .{ .net_erosion_mol_per_step_by_pool = &.{} });
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .vertical,
        .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
        .neighbor_by_boundary_side = &.{},
    }, &state, .{ .net_erosion_mol_per_step_by_pool = &.{} });
}

test "invalid dimension alias and overflow errors preserve state atomically" {
    const three = [_]f64{3} ** mineral_pool_count;
    var sides = [_]NeighborSide{.{
        .local_total_sediment_megagrams_per_step = 2,
        .positive_neighbor_total_sediment_megagrams_per_step = 2,
        .positive_neighbor_mol_per_step_by_pool = &three,
        .connection = .connected,
    }};
    var totals = [_]f64{5} ** mineral_pool_count;
    var scratch: [mineral_pool_count]f64 = undefined;
    var state = State{ .net_erosion_mol_per_step_by_pool = &totals };
    const inputs = Inputs{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    };

    sides[0].positive_neighbor_total_sediment_megagrams_per_step =
        std.math.nan(f64);
    try std.testing.expectError(
        error.NonFinitePrecipitateNeighborErosionInput,
        account(
            inputs,
            &state,
            .{ .net_erosion_mol_per_step_by_pool = &scratch },
        ),
    );
    try expectAll(&totals, 5);

    sides[0].positive_neighbor_total_sediment_megagrams_per_step = 2;
    sides[0].positive_neighbor_mol_per_step_by_pool = three[0..25];
    try std.testing.expectError(
        error.PrecipitateNeighborErosionDimensionMismatch,
        account(
            inputs,
            &state,
            .{ .net_erosion_mol_per_step_by_pool = &scratch },
        ),
    );
    try expectAll(&totals, 5);

    sides[0].positive_neighbor_mol_per_step_by_pool = &three;
    try std.testing.expectError(
        error.PrecipitateNeighborErosionWorkspaceOverlap,
        account(
            inputs,
            &state,
            .{ .net_erosion_mol_per_step_by_pool = &totals },
        ),
    );
    try expectAll(&totals, 5);

    const negative_max = [_]f64{-std.math.floatMax(f64)} **
        mineral_pool_count;
    sides[0].positive_neighbor_mol_per_step_by_pool = &negative_max;
    @memset(&totals, std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFinitePrecipitateNeighborErosionResult,
        account(
            inputs,
            &state,
            .{ .net_erosion_mol_per_step_by_pool = &scratch },
        ),
    );
    try expectAll(&totals, std.math.floatMax(f64));
}
