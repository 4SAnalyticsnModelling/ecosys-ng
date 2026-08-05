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

/// Precipitate and silicate pools in exact REDIST update order.
///
/// The legacy `F` suffix on the second six silicate pools is not defined by
/// the local science comments, so these are named secondary pools rather than
/// assigned an unsupported chemical interpretation.
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

pub const BoundaryFlux = struct {
    total_sediment_megagrams_per_step: f64,
    /// `P*ER`, `Q*ER`, and band `P*EB` [mineral_pool], mol/step.
    mol_per_step_by_pool: []const f64,
};

pub const Inputs = struct {
    disturbance_mode: DisturbanceMode,
    transport_axis: TransportAxis,
    sediment_activity_threshold_megagrams_per_step: f64,
    local_flux_by_boundary_side: []const BoundaryFlux,
    positive_neighbor_total_sediment_megagrams_per_step_by_boundary_side: []const f64,
};

pub const State = struct {
    net_incoming_mol_per_step_by_pool: []f64,
};

pub const Workspace = struct {
    net_incoming_mol_per_step_by_pool: []f64,
};

/// Accounts precipitate and silicate erosion inflow on one horizontal axis.
///
/// Traceability: REDIST.F lines 2975--3000 under gates 2888--2895. Runtime
/// boundary sides replace fixed `NN=1,2`. All 26 local molar pools are added in
/// source order when either local or positive-neighbor sediment exceeds the
/// strict cell threshold. State commits only after finite candidate evaluation.
pub fn account(inputs: Inputs, state: *State, workspace: Workspace) !void {
    if (!erosionEnabled(inputs.disturbance_mode) or
        inputs.transport_axis == .vertical)
    {
        return;
    }
    try validateInputs(inputs, state.*, workspace);
    @memcpy(
        workspace.net_incoming_mol_per_step_by_pool,
        state.net_incoming_mol_per_step_by_pool,
    );
    for (
        inputs.local_flux_by_boundary_side,
        inputs.positive_neighbor_total_sediment_megagrams_per_step_by_boundary_side,
    ) |local, positive_neighbor_total| {
        if (@abs(local.total_sediment_megagrams_per_step) <=
            inputs.sediment_activity_threshold_megagrams_per_step and
            @abs(positive_neighbor_total) <=
                inputs.sediment_activity_threshold_megagrams_per_step)
        {
            continue;
        }
        for (
            workspace.net_incoming_mol_per_step_by_pool,
            local.mol_per_step_by_pool,
        ) |*candidate, contribution| {
            const result = candidate.* + contribution;
            if (!std.math.isFinite(result))
                return error.NonFinitePrecipitateErosionResult;
            candidate.* = result;
        }
    }
    @memcpy(
        state.net_incoming_mol_per_step_by_pool,
        workspace.net_incoming_mol_per_step_by_pool,
    );
}

fn erosionEnabled(mode: DisturbanceMode) bool {
    return mode == .freeze_thaw_and_erosion or
        mode == .freeze_thaw_erosion_and_organic_matter;
}

fn validateInputs(inputs: Inputs, state: State, workspace: Workspace) !void {
    if (inputs.local_flux_by_boundary_side.len == 0)
        return error.InvalidPrecipitateErosionDimensions;
    if (inputs.positive_neighbor_total_sediment_megagrams_per_step_by_boundary_side.len !=
        inputs.local_flux_by_boundary_side.len or
        state.net_incoming_mol_per_step_by_pool.len != mineral_pool_count or
        workspace.net_incoming_mol_per_step_by_pool.len != mineral_pool_count)
    {
        return error.PrecipitateErosionDimensionMismatch;
    }
    if (!std.math.isFinite(inputs.sediment_activity_threshold_megagrams_per_step))
        return error.NonFinitePrecipitateErosionInput;
    if (inputs.sediment_activity_threshold_megagrams_per_step < 0)
        return error.InvalidErosionActivityThreshold;
    try validateFinite(state.net_incoming_mol_per_step_by_pool);
    const scratch = workspace.net_incoming_mol_per_step_by_pool;
    if (overlap(scratch, state.net_incoming_mol_per_step_by_pool))
        return error.PrecipitateErosionWorkspaceOverlap;
    for (inputs.local_flux_by_boundary_side) |flux| {
        if (flux.mol_per_step_by_pool.len != mineral_pool_count)
            return error.PrecipitateErosionDimensionMismatch;
        if (!std.math.isFinite(flux.total_sediment_megagrams_per_step))
            return error.NonFinitePrecipitateErosionInput;
        try validateFinite(flux.mol_per_step_by_pool);
        if (overlap(scratch, flux.mol_per_step_by_pool))
            return error.PrecipitateErosionWorkspaceOverlap;
    }
    for (inputs.positive_neighbor_total_sediment_megagrams_per_step_by_boundary_side) |value|
        if (!std.math.isFinite(value))
            return error.NonFinitePrecipitateErosionInput;
}

fn validateFinite(values: []const f64) !void {
    for (values) |value|
        if (!std.math.isFinite(value))
            return error.NonFinitePrecipitateErosionInput;
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

test "active sides add all 26 precipitate and silicate pools" {
    try std.testing.expectEqual(@as(usize, 26), mineral_pool_count);
    try std.testing.expectEqual(
        @as(usize, 10),
        @intFromEnum(MineralPool.secondary_silicate_aluminum),
    );
    try std.testing.expectEqual(
        @as(usize, 25),
        @intFromEnum(MineralPool.band_precipitated_calcium_dihydrogen_phosphate),
    );
    const two = [_]f64{2} ** mineral_pool_count;
    const three = [_]f64{3} ** mineral_pool_count;
    const local = [_]BoundaryFlux{
        .{ .total_sediment_megagrams_per_step = 2, .mol_per_step_by_pool = &two },
        .{ .total_sediment_megagrams_per_step = -3, .mol_per_step_by_pool = &three },
    };
    const positive = [_]f64{ 0, 0 };
    var totals = [_]f64{100} ** mineral_pool_count;
    var scratch: [mineral_pool_count]f64 = undefined;
    var state = State{ .net_incoming_mol_per_step_by_pool = &totals };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_megagrams_per_step_by_boundary_side = &positive,
    }, &state, .{ .net_incoming_mol_per_step_by_pool = &scratch });
    try expectAll(&totals, 105);
}

test "runtime accepted side sum conserves every mineral pool exactly" {
    const two = [_]f64{2} ** mineral_pool_count;
    const three = [_]f64{-3} ** mineral_pool_count;
    const four = [_]f64{4} ** mineral_pool_count;
    const five = [_]f64{5} ** mineral_pool_count;
    const local = [_]BoundaryFlux{
        .{ .total_sediment_megagrams_per_step = 1, .mol_per_step_by_pool = &two },
        .{ .total_sediment_megagrams_per_step = -3, .mol_per_step_by_pool = &three },
        .{ .total_sediment_megagrams_per_step = 4, .mol_per_step_by_pool = &four },
        .{ .total_sediment_megagrams_per_step = 5, .mol_per_step_by_pool = &five },
    };
    const positive = [_]f64{ 1, 0, 0, 0 };
    var totals = [_]f64{0} ** mineral_pool_count;
    var scratch: [mineral_pool_count]f64 = undefined;
    var state = State{ .net_incoming_mol_per_step_by_pool = &totals };
    try account(.{
        .disturbance_mode = .freeze_thaw_erosion_and_organic_matter,
        .transport_axis = .north_south,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_megagrams_per_step_by_boundary_side = &positive,
    }, &state, .{ .net_incoming_mol_per_step_by_pool = &scratch });
    try expectAll(&totals, 6);
}

test "positive-neighbor sediment alone activates local precipitate flux" {
    const seven = [_]f64{7} ** mineral_pool_count;
    const local = [_]BoundaryFlux{.{
        .total_sediment_megagrams_per_step = 1,
        .mol_per_step_by_pool = &seven,
    }};
    const positive = [_]f64{2};
    var totals = [_]f64{0} ** mineral_pool_count;
    var scratch: [mineral_pool_count]f64 = undefined;
    var state = State{ .net_incoming_mol_per_step_by_pool = &totals };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_megagrams_per_step_by_boundary_side = &positive,
    }, &state, .{ .net_incoming_mol_per_step_by_pool = &scratch });
    try expectAll(&totals, 7);
}

test "disabled and vertical modes bypass unused pool storage" {
    var state = State{ .net_incoming_mol_per_step_by_pool = &.{} };
    try account(.{
        .disturbance_mode = .freeze_thaw,
        .transport_axis = .east_west,
        .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
        .local_flux_by_boundary_side = &.{},
        .positive_neighbor_total_sediment_megagrams_per_step_by_boundary_side = &.{},
    }, &state, .{ .net_incoming_mol_per_step_by_pool = &.{} });
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .vertical,
        .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
        .local_flux_by_boundary_side = &.{},
        .positive_neighbor_total_sediment_megagrams_per_step_by_boundary_side = &.{},
    }, &state, .{ .net_incoming_mol_per_step_by_pool = &.{} });
}

test "invalid overflow dimension and alias failures preserve state" {
    var contribution = [_]f64{2} ** mineral_pool_count;
    contribution[mineral_pool_count - 1] = std.math.nan(f64);
    const local = [_]BoundaryFlux{.{
        .total_sediment_megagrams_per_step = 2,
        .mol_per_step_by_pool = &contribution,
    }};
    const positive = [_]f64{0};
    var totals = [_]f64{5} ** mineral_pool_count;
    var scratch: [mineral_pool_count]f64 = undefined;
    var state = State{ .net_incoming_mol_per_step_by_pool = &totals };
    const inputs = Inputs{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_total_sediment_megagrams_per_step_by_boundary_side = &positive,
    };
    try std.testing.expectError(
        error.NonFinitePrecipitateErosionInput,
        account(inputs, &state, .{ .net_incoming_mol_per_step_by_pool = &scratch }),
    );
    try expectAll(&totals, 5);

    @memset(&contribution, std.math.floatMax(f64));
    @memset(&totals, std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFinitePrecipitateErosionResult,
        account(inputs, &state, .{ .net_incoming_mol_per_step_by_pool = &scratch }),
    );
    try expectAll(&totals, std.math.floatMax(f64));

    try std.testing.expectError(
        error.PrecipitateErosionWorkspaceOverlap,
        account(inputs, &state, .{ .net_incoming_mol_per_step_by_pool = &totals }),
    );
    var short_scratch: [mineral_pool_count - 1]f64 = undefined;
    try std.testing.expectError(
        error.PrecipitateErosionDimensionMismatch,
        account(inputs, &state, .{
            .net_incoming_mol_per_step_by_pool = &short_scratch,
        }),
    );
    try expectAll(&totals, std.math.floatMax(f64));
}
