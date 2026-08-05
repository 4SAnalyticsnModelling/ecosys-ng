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

pub const BoundarySide = enum {
    first,
    second,
};

/// Exchangeable and adsorbed pools in exact REDIST subtraction order.
///
/// The local source does not identify `XAL2ER` and `XFE2ER` chemistry, so
/// these retain neutral secondary-form names.
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

pub const OppositeNeighborFlux = struct {
    total_sediment_megagrams_per_step: f64,
    /// `X*ER` and `X*EB`, mol/step by [exchangeable_pool].
    mol_per_step_by_pool: []const f64,
};

pub const Inputs = struct {
    disturbance_mode: DisturbanceMode,
    transport_axis: TransportAxis,
    boundary_side: BoundarySide,
    sediment_activity_threshold_megagrams_per_step: f64,
    /// Null when the geometry-derived opposite-neighbor coordinate is absent.
    opposite_neighbor_first_side_flux: ?OppositeNeighborFlux,
};

pub const State = struct {
    net_erosion_mol_per_step_by_pool: []f64,
};

pub const Workspace = struct {
    net_erosion_mol_per_step_by_pool: []f64,
};

/// Subtracts exchangeable erosion from the first side's opposite neighbor.
///
/// Traceability: REDIST.F lines 3238--3259 under gates 3192--3193. Optional
/// geometry replaces `N4B,N5B`; `.first` replaces `NN == 1`. All 22 molar
/// inventories retain source order. State commits only after finite candidate
/// evaluation, using caller-owned runtime workspace.
pub fn account(inputs: Inputs, state: *State, workspace: Workspace) !void {
    if (!erosionEnabled(inputs.disturbance_mode) or
        inputs.transport_axis == .vertical or
        inputs.boundary_side != .first)
    {
        return;
    }
    const flux = inputs.opposite_neighbor_first_side_flux orelse return;
    try validateInputs(
        inputs.sediment_activity_threshold_megagrams_per_step,
        flux,
        state.*,
        workspace,
    );
    if (@abs(flux.total_sediment_megagrams_per_step) <=
        inputs.sediment_activity_threshold_megagrams_per_step)
    {
        return;
    }

    @memcpy(
        workspace.net_erosion_mol_per_step_by_pool,
        state.net_erosion_mol_per_step_by_pool,
    );
    for (
        workspace.net_erosion_mol_per_step_by_pool,
        flux.mol_per_step_by_pool,
    ) |*candidate, contribution| {
        const result = candidate.* - contribution;
        if (!std.math.isFinite(result))
            return error.NonFiniteOppositeNeighborExchangeableErosionResult;
        candidate.* = result;
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

fn validateInputs(
    threshold: f64,
    flux: OppositeNeighborFlux,
    state: State,
    workspace: Workspace,
) !void {
    if (state.net_erosion_mol_per_step_by_pool.len != exchangeable_pool_count or
        workspace.net_erosion_mol_per_step_by_pool.len !=
            exchangeable_pool_count or
        flux.mol_per_step_by_pool.len != exchangeable_pool_count)
    {
        return error.OppositeNeighborExchangeableErosionDimensionMismatch;
    }
    if (!std.math.isFinite(threshold) or
        !std.math.isFinite(flux.total_sediment_megagrams_per_step))
    {
        return error.NonFiniteOppositeNeighborExchangeableErosionInput;
    }
    if (threshold < 0) return error.InvalidErosionActivityThreshold;
    try validateFinite(state.net_erosion_mol_per_step_by_pool);
    try validateFinite(flux.mol_per_step_by_pool);
    const scratch = workspace.net_erosion_mol_per_step_by_pool;
    if (overlap(scratch, state.net_erosion_mol_per_step_by_pool) or
        overlap(scratch, flux.mol_per_step_by_pool))
    {
        return error.OppositeNeighborExchangeableErosionWorkspaceOverlap;
    }
}

fn validateFinite(values: []const f64) !void {
    for (values) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteOppositeNeighborExchangeableErosionInput;
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

test "pool order preserves all 22 opposite-neighbor source positions" {
    try std.testing.expectEqual(@as(usize, 22), exchangeable_pool_count);
    try std.testing.expectEqual(
        @as(usize, 10),
        @intFromEnum(ExchangeablePool.adsorbed_aluminum_secondary_form),
    );
    try std.testing.expectEqual(
        @as(usize, 21),
        @intFromEnum(ExchangeablePool.band_adsorbed_dihydrogen_phosphate),
    );
}

test "active opposite neighbor subtracts all exchangeable pools" {
    const three = [_]f64{3} ** exchangeable_pool_count;
    var totals = [_]f64{10} ** exchangeable_pool_count;
    var scratch: [exchangeable_pool_count]f64 = undefined;
    var state = State{ .net_erosion_mol_per_step_by_pool = &totals };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .opposite_neighbor_first_side_flux = .{
            .total_sediment_megagrams_per_step = -2,
            .mol_per_step_by_pool = &three,
        },
    }, &state, .{ .net_erosion_mol_per_step_by_pool = &scratch });
    try expectAll(&totals, 7);
}

test "shared face exchangeable transfer conserves every pool exactly" {
    var shared: [exchangeable_pool_count]f64 = undefined;
    for (&shared, 0..) |*value, index| value.* = @floatFromInt(index + 1);
    var source = shared;
    var scratch: [exchangeable_pool_count]f64 = undefined;
    var state = State{ .net_erosion_mol_per_step_by_pool = &source };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .north_south,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .opposite_neighbor_first_side_flux = .{
            .total_sediment_megagrams_per_step = 2,
            .mol_per_step_by_pool = &shared,
        },
    }, &state, .{ .net_erosion_mol_per_step_by_pool = &scratch });
    try expectAll(&source, 0);
}

test "strict sediment outer side and geometry gates bypass pools" {
    const huge = [_]f64{100} ** exchangeable_pool_count;
    var totals = [_]f64{9} ** exchangeable_pool_count;
    var scratch: [exchangeable_pool_count]f64 = undefined;
    var state = State{ .net_erosion_mol_per_step_by_pool = &totals };
    try account(.{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .opposite_neighbor_first_side_flux = .{
            .total_sediment_megagrams_per_step = 1,
            .mol_per_step_by_pool = &huge,
        },
    }, &state, .{ .net_erosion_mol_per_step_by_pool = &scratch });
    try expectAll(&totals, 9);

    const bypass = [_]Inputs{
        .{
            .disturbance_mode = .freeze_thaw,
            .transport_axis = .east_west,
            .boundary_side = .first,
            .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
            .opposite_neighbor_first_side_flux = null,
        },
        .{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .vertical,
            .boundary_side = .first,
            .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
            .opposite_neighbor_first_side_flux = null,
        },
        .{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .east_west,
            .boundary_side = .second,
            .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
            .opposite_neighbor_first_side_flux = null,
        },
        .{
            .disturbance_mode = .freeze_thaw_and_erosion,
            .transport_axis = .east_west,
            .boundary_side = .first,
            .sediment_activity_threshold_megagrams_per_step = std.math.nan(f64),
            .opposite_neighbor_first_side_flux = null,
        },
    };
    for (bypass) |inputs|
        try account(
            inputs,
            &state,
            .{ .net_erosion_mol_per_step_by_pool = &.{} },
        );
    try expectAll(&totals, 9);
}

test "dimension alias invalid and overflow errors preserve state atomically" {
    const three = [_]f64{3} ** exchangeable_pool_count;
    var totals = [_]f64{5} ** exchangeable_pool_count;
    var scratch: [exchangeable_pool_count]f64 = undefined;
    var state = State{ .net_erosion_mol_per_step_by_pool = &totals };
    var inputs = Inputs{
        .disturbance_mode = .freeze_thaw_and_erosion,
        .transport_axis = .east_west,
        .boundary_side = .first,
        .sediment_activity_threshold_megagrams_per_step = 1,
        .opposite_neighbor_first_side_flux = .{
            .total_sediment_megagrams_per_step = 2,
            .mol_per_step_by_pool = three[0..21],
        },
    };
    try std.testing.expectError(
        error.OppositeNeighborExchangeableErosionDimensionMismatch,
        account(
            inputs,
            &state,
            .{ .net_erosion_mol_per_step_by_pool = &scratch },
        ),
    );
    try expectAll(&totals, 5);

    inputs.opposite_neighbor_first_side_flux.?.mol_per_step_by_pool = &three;
    try std.testing.expectError(
        error.OppositeNeighborExchangeableErosionWorkspaceOverlap,
        account(
            inputs,
            &state,
            .{ .net_erosion_mol_per_step_by_pool = &totals },
        ),
    );
    try expectAll(&totals, 5);

    const negative_max = [_]f64{-std.math.floatMax(f64)} **
        exchangeable_pool_count;
    inputs.opposite_neighbor_first_side_flux.?.mol_per_step_by_pool =
        &negative_max;
    @memset(&totals, std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteOppositeNeighborExchangeableErosionResult,
        account(
            inputs,
            &state,
            .{ .net_erosion_mol_per_step_by_pool = &scratch },
        ),
    );
    try expectAll(&totals, std.math.floatMax(f64));
}
