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

/// Exchangeable and adsorbed pools in exact REDIST subtraction order.
///
/// The local source does not identify the chemistry of `XAL2ER` and `XFE2ER`;
/// they remain secondary forms rather than receiving guessed identities.
pub const ExchangeablePool = enum {
    nonband_adsorbed_ammonium,
    band_adsorbed_ammonium,
    adsorbed_hydrogen,
    adsorbed_aluminum,
    adsorbed_iron,
    adsorbed_calcium,
    adsorbed_magnesium,
    adsorbed_sodium,
    adsorbed_potassium,
    adsorbed_bicarbonate,
    adsorbed_aluminum_secondary_form,
    adsorbed_iron_secondary_form,
    nonband_deprotonated_surface_site,
    nonband_neutral_surface_site,
    nonband_protonated_surface_site,
    nonband_adsorbed_hydrogen_phosphate,
    nonband_adsorbed_dihydrogen_phosphate,
    band_deprotonated_surface_site,
    band_neutral_surface_site,
    band_protonated_surface_site,
    band_adsorbed_hydrogen_phosphate,
    band_adsorbed_dihydrogen_phosphate,
};

pub const exchangeable_pool_count: usize =
    @typeInfo(ExchangeablePool).@"enum".fields.len;

pub const NeighborSide = struct {
    local_total_sediment_megagrams_per_step: f64,
    positive_neighbor_total_sediment_megagrams_per_step: f64,
    /// Positive-neighbor `X*ER` and `X*EB` [exchangeable_pool], mol/step.
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

/// Subtracts connected positive-neighbor exchangeable erosion.
///
/// Traceability: REDIST.F lines 3090--3111, within gates 2888, 2894--2895,
/// and 3051. Runtime sides replace fixed `NN=1,2`; caller-owned pool slices
/// retain all 22 source positions. The strict sediment-activity gate and
/// `IFLBH == 0` connection gate are applied independently. State is committed
/// only after all finite, source-ordered subtractions succeed.
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
                return error.NonFiniteExchangeableNeighborErosionResult;
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
        return error.InvalidExchangeableNeighborErosionDimensions;
    if (state.net_erosion_mol_per_step_by_pool.len != exchangeable_pool_count or
        workspace.net_erosion_mol_per_step_by_pool.len != exchangeable_pool_count)
    {
        return error.ExchangeableNeighborErosionDimensionMismatch;
    }
    if (!std.math.isFinite(inputs.sediment_activity_threshold_megagrams_per_step))
        return error.NonFiniteExchangeableNeighborErosionInput;
    if (inputs.sediment_activity_threshold_megagrams_per_step < 0)
        return error.InvalidErosionActivityThreshold;

    try validateFinite(state.net_erosion_mol_per_step_by_pool);
    const scratch = workspace.net_erosion_mol_per_step_by_pool;
    if (overlap(scratch, state.net_erosion_mol_per_step_by_pool))
        return error.ExchangeableNeighborErosionWorkspaceOverlap;
    for (inputs.neighbor_by_boundary_side) |side| {
        if (side.positive_neighbor_mol_per_step_by_pool.len !=
            exchangeable_pool_count)
        {
            return error.ExchangeableNeighborErosionDimensionMismatch;
        }
        if (!std.math.isFinite(side.local_total_sediment_megagrams_per_step) or
            !std.math.isFinite(
                side.positive_neighbor_total_sediment_megagrams_per_step,
            ))
        {
            return error.NonFiniteExchangeableNeighborErosionInput;
        }
        try validateFinite(side.positive_neighbor_mol_per_step_by_pool);
        if (overlap(scratch, side.positive_neighbor_mol_per_step_by_pool))
            return error.ExchangeableNeighborErosionWorkspaceOverlap;
    }
}

fn validateFinite(values: []const f64) !void {
    for (values) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteExchangeableNeighborErosionInput;
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

test "pool order preserves all 22 source positions" {
    try std.testing.expectEqual(@as(usize, 22), exchangeable_pool_count);
    try std.testing.expectEqual(
        @as(usize, 0),
        @intFromEnum(ExchangeablePool.nonband_adsorbed_ammonium),
    );
    try std.testing.expectEqual(
        @as(usize, 10),
        @intFromEnum(ExchangeablePool.adsorbed_aluminum_secondary_form),
    );
    try std.testing.expectEqual(
        @as(usize, 21),
        @intFromEnum(ExchangeablePool.band_adsorbed_dihydrogen_phosphate),
    );
}

test "connected active sides subtract all exchangeable pools" {
    const three = [_]f64{3} ** exchangeable_pool_count;
    const five = [_]f64{5} ** exchangeable_pool_count;
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
    var totals = [_]f64{100} ** exchangeable_pool_count;
    var scratch: [exchangeable_pool_count]f64 = undefined;
    var state = State{ .net_erosion_mol_per_step_by_pool = &totals };

    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    }, &state, .{ .net_erosion_mol_per_step_by_pool = &scratch });
    try expectAll(&totals, 92);
}

test "strict activity and connection gates remain independent" {
    const seven = [_]f64{7} ** exchangeable_pool_count;
    const four = [_]f64{4} ** exchangeable_pool_count;
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
    var totals = [_]f64{10} ** exchangeable_pool_count;
    var scratch: [exchangeable_pool_count]f64 = undefined;
    var state = State{ .net_erosion_mol_per_step_by_pool = &totals };

    try account(.{
        .disturbance_mode = .freeze_thaw_erosion_and_organic_matter,
        .transport_axis = .north_south,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    }, &state, .{ .net_erosion_mol_per_step_by_pool = &scratch });
    try expectAll(&totals, 6);
}

test "shared face transfer conserves every exchangeable pool exactly" {
    var shared: [exchangeable_pool_count]f64 = undefined;
    for (&shared, 0..) |*value, index| value.* = @floatFromInt(index + 1);
    var source = shared;
    const sides = [_]NeighborSide{.{
        .local_total_sediment_megagrams_per_step = 0,
        .positive_neighbor_total_sediment_megagrams_per_step = 2,
        .positive_neighbor_mol_per_step_by_pool = &shared,
        .connection = .connected,
    }};
    var scratch: [exchangeable_pool_count]f64 = undefined;
    var state = State{ .net_erosion_mol_per_step_by_pool = &source };

    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .neighbor_by_boundary_side = &sides,
    }, &state, .{ .net_erosion_mol_per_step_by_pool = &scratch });
    try expectAll(&source, 0);
}

test "disabled and vertical modes bypass unused neighbor storage" {
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

test "dimension invalid finite overflow and alias errors preserve state" {
    const three = [_]f64{3} ** exchangeable_pool_count;
    var sides = [_]NeighborSide{.{
        .local_total_sediment_megagrams_per_step = 2,
        .positive_neighbor_total_sediment_megagrams_per_step = 2,
        .positive_neighbor_mol_per_step_by_pool = &three,
        .connection = .connected,
    }};
    var totals = [_]f64{5} ** exchangeable_pool_count;
    var scratch: [exchangeable_pool_count]f64 = undefined;
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
        error.NonFiniteExchangeableNeighborErosionInput,
        account(
            inputs,
            &state,
            .{ .net_erosion_mol_per_step_by_pool = &scratch },
        ),
    );
    try expectAll(&totals, 5);

    sides[0].positive_neighbor_total_sediment_megagrams_per_step = 2;
    sides[0].positive_neighbor_mol_per_step_by_pool = three[0..21];
    try std.testing.expectError(
        error.ExchangeableNeighborErosionDimensionMismatch,
        account(
            inputs,
            &state,
            .{ .net_erosion_mol_per_step_by_pool = &scratch },
        ),
    );
    try expectAll(&totals, 5);

    sides[0].positive_neighbor_mol_per_step_by_pool = &three;
    try std.testing.expectError(
        error.ExchangeableNeighborErosionWorkspaceOverlap,
        account(
            inputs,
            &state,
            .{ .net_erosion_mol_per_step_by_pool = &totals },
        ),
    );
    try expectAll(&totals, 5);

    const negative_max = [_]f64{-std.math.floatMax(f64)} **
        exchangeable_pool_count;
    sides[0].positive_neighbor_mol_per_step_by_pool = &negative_max;
    @memset(&totals, std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteExchangeableNeighborErosionResult,
        account(
            inputs,
            &state,
            .{ .net_erosion_mol_per_step_by_pool = &scratch },
        ),
    );
    try expectAll(&totals, std.math.floatMax(f64));
}
