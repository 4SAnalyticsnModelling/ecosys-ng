const std = @import("std");

pub const Parameters = struct {
    solid_snow_heat_capacity_mj_per_m3_k: f64,
    liquid_water_heat_capacity_mj_per_m3_k: f64,
    ice_heat_capacity_mj_per_m3_k: f64,
    celsius_zero_temperature_k: f64,
    absolute_heat_capacity_tolerance_mj_per_k: f64,
    relative_heat_capacity_tolerance: f64,
};

pub const CellInputs = struct {
    horizontal_area_m2: f64,
    atmospheric_temperature_k: f64,
    initial_snow_density_Mg_per_m3: f64,
    inactive_phase_volume_threshold_m3: f64,
};

pub const LayerState = struct {
    active: []bool,
    solid_snow_water_equivalent_m3: []f64,
    liquid_water_m3: []f64,
    water_vapor_equivalent_m3: []f64,
    ice_volume_m3: []f64,
    snow_density_Mg_per_m3: []f64,
    temperature_k: []f64,
    temperature_c: []f64,
    heat_capacity_mj_per_k: []f64,
    total_layer_volume_m3: []f64,
    layer_thickness_m: []f64,
    cumulative_depth_m: []f64,
};

pub const ColumnTotals = struct {
    solid_snow_water_equivalent_m3: f64 = 0,
    liquid_water_m3: f64 = 0,
    ice_volume_m3: f64 = 0,
    snowpack_volume_m3: f64 = 0,
    snowpack_depth_m: f64 = 0,
};

pub const ClearedInactiveStorage = struct {
    solid_snow_water_equivalent_m3: f64 = 0,
    liquid_water_m3: f64 = 0,
    water_vapor_equivalent_m3: f64 = 0,
    ice_volume_m3: f64 = 0,
    sensible_energy_mj: f64 = 0,
};

pub const Report = struct {
    column_totals: ColumnTotals = .{},
    cleared_inactive_storage: ClearedInactiveStorage = .{},
    inactive_layer_count: usize = 0,
};

pub const LayerCandidate = struct {
    active: bool,
    solid_snow_water_equivalent_m3: f64,
    liquid_water_m3: f64,
    water_vapor_equivalent_m3: f64,
    ice_volume_m3: f64,
    snow_density_Mg_per_m3: f64,
    temperature_k: f64,
    temperature_c: f64,
    heat_capacity_mj_per_k: f64,
    total_layer_volume_m3: f64,
    layer_thickness_m: f64,
    cumulative_depth_m: f64,
};

/// Publishes enhanced snow-solver state into layer geometry and column totals.
///
/// Traceability: `REDIST` lines 4046--4133. Geometry, the active-layer gate,
/// column accumulation, inactive-layer reset, and Celsius conversion retain
/// source order. The legacy `THFLWW + XHFLF0 + XHFLV0` temperature update is
/// intentionally absent: ecosys-ng's heat, vapor, and Dall'Amico phase solvers
/// have already committed `temperature_k` and `heat_capacity_mj_per_k`.
/// Reapplying those energy increments here would count them twice.
///
/// One invocation owns one cell's runtime layer slices. Independent cells can
/// therefore run in parallel without shared mutation. The caller supplies
/// `workspace` so validation and calculation complete before atomic commit.
pub fn publishCell(
    state: LayerState,
    inputs: CellInputs,
    parameters: Parameters,
    workspace: []LayerCandidate,
) !Report {
    try validateDimensions(state, workspace);
    try validateParameters(inputs, parameters);
    try validateWorkspaceOwnership(state, workspace);

    const report = try prepareCandidates(state, inputs, parameters, workspace);
    commitCandidates(state, workspace);
    return report;
}

fn prepareCandidates(
    state: LayerState,
    inputs: CellInputs,
    parameters: Parameters,
    workspace: []LayerCandidate,
) !Report {
    var report: Report = .{};
    var cumulative_depth_m: f64 = 0;
    for (workspace, 0..) |*candidate, layer| {
        try validateLayerInputs(state, layer);
        const phase_volume_m3 =
            state.solid_snow_water_equivalent_m3[layer] +
            state.liquid_water_m3[layer] +
            state.ice_volume_m3[layer];
        if (!std.math.isFinite(phase_volume_m3))
            return error.NonFiniteSnowpackPublicationResult;

        if (phase_volume_m3 > inputs.inactive_phase_volume_threshold_m3) {
            candidate.* = try prepareActiveLayer(
                state,
                inputs,
                parameters,
                layer,
                cumulative_depth_m,
            );
            cumulative_depth_m = candidate.cumulative_depth_m;
            try addActiveLayerToTotals(&report.column_totals, candidate.*);
        } else {
            candidate.* = try prepareInactiveLayer(
                state,
                inputs,
                parameters,
                workspace,
                layer,
                cumulative_depth_m,
            );
            try recordInactiveCleanup(&report, state, layer);
        }
    }
    try validateReport(report);
    return report;
}

fn prepareActiveLayer(
    state: LayerState,
    inputs: CellInputs,
    parameters: Parameters,
    layer: usize,
    previous_depth_m: f64,
) !LayerCandidate {
    const density_Mg_per_m3 = state.snow_density_Mg_per_m3[layer];
    if (density_Mg_per_m3 <= 0) return error.InvalidActiveSnowDensity;
    const total_volume_m3 =
        state.solid_snow_water_equivalent_m3[layer] / density_Mg_per_m3 +
        state.liquid_water_m3[layer] +
        state.ice_volume_m3[layer];
    const thickness_m = @max(0, total_volume_m3) / inputs.horizontal_area_m2;
    const cumulative_depth_m = previous_depth_m + thickness_m;
    const calculated_capacity_mj_per_k =
        parameters.solid_snow_heat_capacity_mj_per_m3_k *
        state.solid_snow_water_equivalent_m3[layer] +
        parameters.liquid_water_heat_capacity_mj_per_m3_k *
            (state.liquid_water_m3[layer] +
                state.water_vapor_equivalent_m3[layer]) +
        parameters.ice_heat_capacity_mj_per_m3_k *
            state.ice_volume_m3[layer];
    inline for (.{
        total_volume_m3,
        thickness_m,
        cumulative_depth_m,
        calculated_capacity_mj_per_k,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteSnowpackPublicationResult;
    try requireConsistentHeatCapacity(
        state.heat_capacity_mj_per_k[layer],
        calculated_capacity_mj_per_k,
        parameters,
    );
    const temperature_c =
        state.temperature_k[layer] - parameters.celsius_zero_temperature_k;
    if (!std.math.isFinite(temperature_c))
        return error.NonFiniteSnowpackPublicationResult;

    return .{
        .active = true,
        .solid_snow_water_equivalent_m3 = state.solid_snow_water_equivalent_m3[layer],
        .liquid_water_m3 = state.liquid_water_m3[layer],
        .water_vapor_equivalent_m3 = state.water_vapor_equivalent_m3[layer],
        .ice_volume_m3 = state.ice_volume_m3[layer],
        .snow_density_Mg_per_m3 = density_Mg_per_m3,
        .temperature_k = state.temperature_k[layer],
        .temperature_c = temperature_c,
        .heat_capacity_mj_per_k = calculated_capacity_mj_per_k,
        .total_layer_volume_m3 = total_volume_m3,
        .layer_thickness_m = thickness_m,
        .cumulative_depth_m = cumulative_depth_m,
    };
}

fn prepareInactiveLayer(
    state: LayerState,
    inputs: CellInputs,
    parameters: Parameters,
    workspace: []const LayerCandidate,
    layer: usize,
    cumulative_depth_m: f64,
) !LayerCandidate {
    const temperature_k = if (layer == 0)
        inputs.atmospheric_temperature_k
    else
        workspace[layer - 1].temperature_k;
    const density_Mg_per_m3 = if (layer == 0)
        inputs.initial_snow_density_Mg_per_m3
    else
        workspace[layer - 1].snow_density_Mg_per_m3;
    const temperature_c =
        temperature_k - parameters.celsius_zero_temperature_k;
    inline for (.{ temperature_k, density_Mg_per_m3, temperature_c }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteSnowpackPublicationResult;
    _ = state;

    return .{
        .active = false,
        .solid_snow_water_equivalent_m3 = 0,
        .liquid_water_m3 = 0,
        .water_vapor_equivalent_m3 = 0,
        .ice_volume_m3 = 0,
        .snow_density_Mg_per_m3 = density_Mg_per_m3,
        .temperature_k = temperature_k,
        .temperature_c = temperature_c,
        .heat_capacity_mj_per_k = 0,
        .total_layer_volume_m3 = 0,
        .layer_thickness_m = 0,
        .cumulative_depth_m = cumulative_depth_m,
    };
}

fn addActiveLayerToTotals(
    totals: *ColumnTotals,
    layer: LayerCandidate,
) !void {
    totals.solid_snow_water_equivalent_m3 +=
        layer.solid_snow_water_equivalent_m3;
    totals.liquid_water_m3 += layer.liquid_water_m3;
    totals.ice_volume_m3 += layer.ice_volume_m3;
    totals.snowpack_volume_m3 += layer.total_layer_volume_m3;
    totals.snowpack_depth_m += layer.layer_thickness_m;
    inline for (@typeInfo(ColumnTotals).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(totals.*, field.name)))
            return error.NonFiniteSnowpackPublicationResult;
    }
}

fn recordInactiveCleanup(
    report: *Report,
    state: LayerState,
    layer: usize,
) !void {
    report.inactive_layer_count = std.math.add(
        usize,
        report.inactive_layer_count,
        1,
    ) catch return error.SnowpackPublicationCountOverflow;
    report.cleared_inactive_storage.solid_snow_water_equivalent_m3 +=
        state.solid_snow_water_equivalent_m3[layer];
    report.cleared_inactive_storage.liquid_water_m3 +=
        state.liquid_water_m3[layer];
    report.cleared_inactive_storage.water_vapor_equivalent_m3 +=
        state.water_vapor_equivalent_m3[layer];
    report.cleared_inactive_storage.ice_volume_m3 +=
        state.ice_volume_m3[layer];
    report.cleared_inactive_storage.sensible_energy_mj +=
        state.heat_capacity_mj_per_k[layer] * state.temperature_k[layer];
}

fn commitCandidates(state: LayerState, candidates: []const LayerCandidate) void {
    for (candidates, 0..) |candidate, layer| {
        state.active[layer] = candidate.active;
        inline for (.{
            .{ state.solid_snow_water_equivalent_m3, candidate.solid_snow_water_equivalent_m3 },
            .{ state.liquid_water_m3, candidate.liquid_water_m3 },
            .{ state.water_vapor_equivalent_m3, candidate.water_vapor_equivalent_m3 },
            .{ state.ice_volume_m3, candidate.ice_volume_m3 },
            .{ state.snow_density_Mg_per_m3, candidate.snow_density_Mg_per_m3 },
            .{ state.temperature_k, candidate.temperature_k },
            .{ state.temperature_c, candidate.temperature_c },
            .{ state.heat_capacity_mj_per_k, candidate.heat_capacity_mj_per_k },
            .{ state.total_layer_volume_m3, candidate.total_layer_volume_m3 },
            .{ state.layer_thickness_m, candidate.layer_thickness_m },
            .{ state.cumulative_depth_m, candidate.cumulative_depth_m },
        }) |destination_and_value| {
            destination_and_value[0][layer] = destination_and_value[1];
        }
    }
}

fn validateDimensions(state: LayerState, workspace: []LayerCandidate) !void {
    const layer_count = state.active.len;
    if (layer_count == 0) return error.InvalidSnowpackPublicationDimensions;
    inline for (@typeInfo(LayerState).@"struct".fields) |field| {
        if (@field(state, field.name).len != layer_count)
            return error.SnowpackPublicationDimensionMismatch;
    }
    if (workspace.len != layer_count)
        return error.SnowpackPublicationDimensionMismatch;
}

fn validateParameters(inputs: CellInputs, parameters: Parameters) !void {
    inline for (@typeInfo(CellInputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteSnowpackPublicationParameter;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(parameters, field.name)))
            return error.NonFiniteSnowpackPublicationParameter;
    }
    if (inputs.horizontal_area_m2 <= 0 or
        inputs.atmospheric_temperature_k <= 0 or
        inputs.initial_snow_density_Mg_per_m3 <= 0 or
        inputs.inactive_phase_volume_threshold_m3 < 0 or
        parameters.solid_snow_heat_capacity_mj_per_m3_k <= 0 or
        parameters.liquid_water_heat_capacity_mj_per_m3_k <= 0 or
        parameters.ice_heat_capacity_mj_per_m3_k <= 0 or
        parameters.celsius_zero_temperature_k <= 0 or
        parameters.absolute_heat_capacity_tolerance_mj_per_k < 0 or
        parameters.relative_heat_capacity_tolerance < 0)
    {
        return error.InvalidSnowpackPublicationParameter;
    }
}

fn validateLayerInputs(state: LayerState, layer: usize) !void {
    inline for (.{
        state.solid_snow_water_equivalent_m3[layer],
        state.liquid_water_m3[layer],
        state.water_vapor_equivalent_m3[layer],
        state.ice_volume_m3[layer],
        state.snow_density_Mg_per_m3[layer],
        state.temperature_k[layer],
        state.heat_capacity_mj_per_k[layer],
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteSnowpackPublicationState;
    if (state.solid_snow_water_equivalent_m3[layer] < 0 or
        state.liquid_water_m3[layer] < 0 or
        state.water_vapor_equivalent_m3[layer] < 0 or
        state.ice_volume_m3[layer] < 0 or
        state.snow_density_Mg_per_m3[layer] < 0 or
        state.temperature_k[layer] <= 0 or
        state.heat_capacity_mj_per_k[layer] < 0)
    {
        return error.InvalidSnowpackPublicationState;
    }
}

fn requireConsistentHeatCapacity(
    committed_mj_per_k: f64,
    calculated_mj_per_k: f64,
    parameters: Parameters,
) !void {
    const tolerance_mj_per_k =
        parameters.absolute_heat_capacity_tolerance_mj_per_k +
        parameters.relative_heat_capacity_tolerance *
            @max(@abs(committed_mj_per_k), @abs(calculated_mj_per_k));
    if (!std.math.isFinite(tolerance_mj_per_k))
        return error.NonFiniteSnowpackPublicationResult;
    if (@abs(committed_mj_per_k - calculated_mj_per_k) >
        tolerance_mj_per_k)
    {
        return error.InconsistentCommittedSnowHeatCapacity;
    }
}

fn validateReport(report: Report) !void {
    inline for (@typeInfo(ColumnTotals).@"struct".fields) |field| {
        const value = @field(report.column_totals, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSnowpackPublicationResult;
        if (value < 0) return error.InvalidSnowpackPublicationResult;
    }
    inline for (@typeInfo(ClearedInactiveStorage).@"struct".fields) |field| {
        const value = @field(report.cleared_inactive_storage, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteSnowpackPublicationResult;
        if (value < 0) return error.InvalidSnowpackPublicationResult;
    }
}

fn validateWorkspaceOwnership(
    state: LayerState,
    workspace: []LayerCandidate,
) !void {
    const workspace_start = @intFromPtr(workspace.ptr);
    const workspace_bytes = std.math.mul(
        usize,
        workspace.len,
        @sizeOf(LayerCandidate),
    ) catch return error.SnowpackPublicationDimensionOverflow;
    inline for (@typeInfo(LayerState).@"struct".fields) |field| {
        const values = @field(state, field.name);
        const value_type = @typeInfo(@TypeOf(values)).pointer.child;
        const state_bytes = std.math.mul(
            usize,
            values.len,
            @sizeOf(value_type),
        ) catch return error.SnowpackPublicationDimensionOverflow;
        if (rangesOverlap(
            workspace_start,
            workspace_bytes,
            @intFromPtr(values.ptr),
            state_bytes,
        )) return error.SnowpackPublicationWorkspaceOverlap;
    }
}

fn rangesOverlap(
    left_start: usize,
    left_bytes: usize,
    right_start: usize,
    right_bytes: usize,
) bool {
    const left_end = std.math.add(usize, left_start, left_bytes) catch
        return true;
    const right_end = std.math.add(usize, right_start, right_bytes) catch
        return true;
    return left_start < right_end and right_start < left_end;
}

const TestStorage = struct {
    allocator: std.mem.Allocator,
    state: LayerState,

    fn init(allocator: std.mem.Allocator, layer_count: usize) !TestStorage {
        var storage: TestStorage = undefined;
        storage.allocator = allocator;
        inline for (@typeInfo(LayerState).@"struct".fields) |field| {
            const slice_type = field.type;
            const value_type = @typeInfo(slice_type).pointer.child;
            @field(storage.state, field.name) =
                try allocator.alloc(value_type, layer_count);
        }
        @memset(storage.state.active, false);
        inline for (@typeInfo(LayerState).@"struct".fields[1..]) |field|
            @memset(@field(storage.state, field.name), 0);
        return storage;
    }

    fn deinit(self: *TestStorage) void {
        inline for (@typeInfo(LayerState).@"struct".fields) |field|
            self.allocator.free(@field(self.state, field.name));
        self.* = undefined;
    }
};

const test_parameters: Parameters = .{
    .solid_snow_heat_capacity_mj_per_m3_k = 2,
    .liquid_water_heat_capacity_mj_per_m3_k = 4,
    .ice_heat_capacity_mj_per_m3_k = 1.5,
    .celsius_zero_temperature_k = 273.15,
    .absolute_heat_capacity_tolerance_mj_per_k = 1e-12,
    .relative_heat_capacity_tolerance = 1e-12,
};

const test_inputs: CellInputs = .{
    .horizontal_area_m2 = 2,
    .atmospheric_temperature_k = 265,
    .initial_snow_density_Mg_per_m3 = 0.2,
    .inactive_phase_volume_threshold_m3 = 1e-9,
};

fn setLayer(
    state: LayerState,
    layer: usize,
    solid_m3: f64,
    liquid_m3: f64,
    vapor_m3: f64,
    ice_m3: f64,
    density_Mg_per_m3: f64,
    temperature_k: f64,
) void {
    state.solid_snow_water_equivalent_m3[layer] = solid_m3;
    state.liquid_water_m3[layer] = liquid_m3;
    state.water_vapor_equivalent_m3[layer] = vapor_m3;
    state.ice_volume_m3[layer] = ice_m3;
    state.snow_density_Mg_per_m3[layer] = density_Mg_per_m3;
    state.temperature_k[layer] = temperature_k;
    state.heat_capacity_mj_per_k[layer] =
        test_parameters.solid_snow_heat_capacity_mj_per_m3_k * solid_m3 +
        test_parameters.liquid_water_heat_capacity_mj_per_m3_k *
            (liquid_m3 + vapor_m3) +
        test_parameters.ice_heat_capacity_mj_per_m3_k * ice_m3;
}

test "REDIST active snow geometry preserves enhanced thermal publication" {
    var storage = try TestStorage.init(std.testing.allocator, 2);
    defer storage.deinit();
    setLayer(storage.state, 0, 0.2, 0.1, 0.01, 0.05, 0.4, 270);
    setLayer(storage.state, 1, 0.1, 0, 0.02, 0, 0.5, 268);
    const workspace =
        try std.testing.allocator.alloc(LayerCandidate, 2);
    defer std.testing.allocator.free(workspace);

    const report = try publishCell(
        storage.state,
        test_inputs,
        test_parameters,
        workspace,
    );

    try std.testing.expect(storage.state.active[0]);
    try std.testing.expect(storage.state.active[1]);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.65),
        storage.state.total_layer_volume_m3[0],
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.325),
        storage.state.layer_thickness_m[0],
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.425),
        storage.state.cumulative_depth_m[1],
        1e-14,
    );
    try std.testing.expectEqual(@as(f64, 270), storage.state.temperature_k[0]);
    try std.testing.expectApproxEqAbs(
        @as(f64, -3.15),
        storage.state.temperature_c[0],
        1e-13,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.3),
        report.column_totals.solid_snow_water_equivalent_m3,
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.85),
        report.column_totals.snowpack_volume_m3,
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        report.column_totals.snowpack_depth_m,
        storage.state.cumulative_depth_m[1],
        1e-14,
    );
}

test "inactive layers reset explicitly and inherit source-order boundary state" {
    var storage = try TestStorage.init(std.testing.allocator, 3);
    defer storage.deinit();
    setLayer(storage.state, 0, 0, 0, 4e-10, 0, 0, 275);
    setLayer(storage.state, 1, 0.12, 0, 0, 0, 0.3, 269);
    setLayer(storage.state, 2, 2e-10, 0, 3e-10, 0, 0.1, 271);
    const workspace =
        try std.testing.allocator.alloc(LayerCandidate, 3);
    defer std.testing.allocator.free(workspace);

    const report = try publishCell(
        storage.state,
        test_inputs,
        test_parameters,
        workspace,
    );

    try std.testing.expect(!storage.state.active[0]);
    try std.testing.expectEqual(
        test_inputs.atmospheric_temperature_k,
        storage.state.temperature_k[0],
    );
    try std.testing.expectEqual(
        test_inputs.initial_snow_density_Mg_per_m3,
        storage.state.snow_density_Mg_per_m3[0],
    );
    try std.testing.expect(storage.state.active[1]);
    try std.testing.expect(!storage.state.active[2]);
    try std.testing.expectEqual(
        storage.state.temperature_k[1],
        storage.state.temperature_k[2],
    );
    try std.testing.expectEqual(
        storage.state.snow_density_Mg_per_m3[1],
        storage.state.snow_density_Mg_per_m3[2],
    );
    try std.testing.expectEqual(@as(f64, 0), storage.state.temperature_c[2] -
        (storage.state.temperature_k[2] -
            test_parameters.celsius_zero_temperature_k));
    try std.testing.expectEqual(@as(usize, 2), report.inactive_layer_count);
    try std.testing.expectApproxEqAbs(
        @as(f64, 7e-10),
        report.cleared_inactive_storage.water_vapor_equivalent_m3,
        1e-20,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 2e-10),
        report.cleared_inactive_storage.solid_snow_water_equivalent_m3,
        1e-20,
    );
    try std.testing.expectEqual(
        storage.state.cumulative_depth_m[1],
        storage.state.cumulative_depth_m[2],
    );
}

test "runtime layer count has no fixed snowpack extent" {
    const layer_count: usize = 17;
    var storage = try TestStorage.init(std.testing.allocator, layer_count);
    defer storage.deinit();
    for (0..layer_count) |layer|
        setLayer(storage.state, layer, 0.01, 0, 0, 0, 0.25, 266);
    const workspace =
        try std.testing.allocator.alloc(LayerCandidate, layer_count);
    defer std.testing.allocator.free(workspace);

    const report = try publishCell(
        storage.state,
        test_inputs,
        test_parameters,
        workspace,
    );

    try std.testing.expectApproxEqAbs(
        @as(f64, 0.17),
        report.column_totals.solid_snow_water_equivalent_m3,
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.34),
        report.column_totals.snowpack_depth_m,
        1e-14,
    );
}

test "heat-capacity discrepancy fails without partial state commit" {
    var storage = try TestStorage.init(std.testing.allocator, 2);
    defer storage.deinit();
    setLayer(storage.state, 0, 0.1, 0, 0, 0, 0.25, 268);
    setLayer(storage.state, 1, 0.1, 0, 0, 0, 0.25, 268);
    storage.state.heat_capacity_mj_per_k[1] += 0.1;
    @memset(storage.state.total_layer_volume_m3, 99);
    @memset(storage.state.layer_thickness_m, 98);
    const workspace =
        try std.testing.allocator.alloc(LayerCandidate, 2);
    defer std.testing.allocator.free(workspace);

    try std.testing.expectError(
        error.InconsistentCommittedSnowHeatCapacity,
        publishCell(storage.state, test_inputs, test_parameters, workspace),
    );
    try std.testing.expectEqual(
        @as(f64, 99),
        storage.state.total_layer_volume_m3[0],
    );
    try std.testing.expectEqual(
        @as(f64, 98),
        storage.state.layer_thickness_m[0],
    );
    try std.testing.expect(!storage.state.active[0]);
}

test "dimension and non-finite failures are diagnostic" {
    var storage = try TestStorage.init(std.testing.allocator, 2);
    defer storage.deinit();
    setLayer(storage.state, 0, 0.1, 0, 0, 0, 0.25, 268);
    setLayer(storage.state, 1, 0.1, 0, 0, 0, 0.25, 268);
    const short_workspace =
        try std.testing.allocator.alloc(LayerCandidate, 1);
    defer std.testing.allocator.free(short_workspace);
    try std.testing.expectError(
        error.SnowpackPublicationDimensionMismatch,
        publishCell(
            storage.state,
            test_inputs,
            test_parameters,
            short_workspace,
        ),
    );

    const workspace =
        try std.testing.allocator.alloc(LayerCandidate, 2);
    defer std.testing.allocator.free(workspace);
    storage.state.temperature_k[1] = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteSnowpackPublicationState,
        publishCell(storage.state, test_inputs, test_parameters, workspace),
    );
}
