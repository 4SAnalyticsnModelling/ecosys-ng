const std = @import("std");

pub const NeighborFluxMode = enum {
    blocked,
    subtract,
};

pub const RunoffWaterHeatFlux = struct {
    water_m3_per_step: f64 = 0,
    convective_heat_megajoules_per_step: f64 = 0,
};

pub const SnowWaterHeatFlux = struct {
    snow_m3_per_step: f64 = 0,
    liquid_water_m3_per_step: f64 = 0,
    ice_m3_per_step: f64 = 0,
    convective_heat_megajoules_per_step: f64 = 0,
};

pub const OrganicFlux = struct {
    carbon_g_c_per_step: f64 = 0,
    nitrogen_g_n_per_step: f64 = 0,
    phosphorus_g_p_per_step: f64 = 0,
    acetate_carbon_g_c_per_step: f64 = 0,
};

pub const RunoffElementFlux = struct {
    carbon_dioxide_g_c_per_step: f64 = 0,
    methane_g_c_per_step: f64 = 0,
    oxygen_g_o_per_step: f64 = 0,
    dinitrogen_g_n_per_step: f64 = 0,
    nitrous_oxide_g_n_per_step: f64 = 0,
    hydrogen_g_h_per_step: f64 = 0,
    ammonium_g_n_per_step: f64 = 0,
    ammonia_g_n_per_step: f64 = 0,
    nitrate_g_n_per_step: f64 = 0,
    nitrite_g_n_per_step: f64 = 0,
    hpo4_g_p_per_step: f64 = 0,
    h2po4_g_p_per_step: f64 = 0,
};

pub const SnowElementFlux = struct {
    carbon_dioxide_g_c_per_step: f64 = 0,
    methane_g_c_per_step: f64 = 0,
    oxygen_g_o_per_step: f64 = 0,
    dinitrogen_g_n_per_step: f64 = 0,
    nitrous_oxide_g_n_per_step: f64 = 0,
    ammonium_g_n_per_step: f64 = 0,
    ammonia_g_n_per_step: f64 = 0,
    nitrate_g_n_per_step: f64 = 0,
    hpo4_g_p_per_step: f64 = 0,
    h2po4_g_p_per_step: f64 = 0,
};

pub const SurfaceFlux = struct {
    runoff_water_heat: RunoffWaterHeatFlux = .{},
    snow_water_heat: SnowWaterHeatFlux = .{},
    organic_by_class: []const OrganicFlux,
    runoff_element: RunoffElementFlux = .{},
    snow_element: SnowElementFlux = .{},
};

pub const PositiveNeighbor = struct {
    flux: ?SurfaceFlux,
    runoff_mode: NeighborFluxMode,
    snow_mode: NeighborFluxMode,
};

pub const Inputs = struct {
    organic_class_count: usize,
    local_flux_by_boundary_side: []const SurfaceFlux,
    positive_neighbor_by_boundary_side: []const PositiveNeighbor,
    /// REDIST subtracts this only during `NN=1`.
    opposite_neighbor_first_side: ?SurfaceFlux,
};

pub const State = struct {
    runoff_water_heat: RunoffWaterHeatFlux,
    snow_water_heat: SnowWaterHeatFlux,
    organic_by_class: []OrganicFlux,
    runoff_element: RunoffElementFlux,
    snow_element: SnowElementFlux,
};

/// Caller-owned scratch makes candidate evaluation atomic without allocating
/// inside a per-cell transport kernel.
pub const Workspace = struct {
    organic_by_class: []OrganicFlux,
};

const Operation = enum {
    add,
    subtract,
};

/// Aggregates one surface cell along one horizontal axis.
///
/// Traceability: REDIST.F lines 1957--2111 (`TQR`--`TPOQSS`). The caller
/// invokes this only for the surface soil layer and column/row axes, matching
/// the source gates. Each runtime boundary side is applied in source order:
/// local addition, gated positive-neighbor runoff subtraction, gated
/// positive-neighbor snow subtraction, then the first-side opposite-neighbor
/// subtraction. State is committed only after every candidate remains finite.
pub fn aggregate(inputs: Inputs, state: *State, workspace: Workspace) !void {
    try validateDimensions(inputs, state.*, workspace);
    try validateInputs(inputs, state.*, workspace);

    var runoff_water_heat = state.runoff_water_heat;
    var snow_water_heat = state.snow_water_heat;
    var runoff_element = state.runoff_element;
    var snow_element = state.snow_element;
    @memcpy(workspace.organic_by_class, state.organic_by_class);

    for (inputs.local_flux_by_boundary_side, 0..) |local, side| {
        try updateWhole(
            &runoff_water_heat,
            &snow_water_heat,
            workspace.organic_by_class,
            &runoff_element,
            &snow_element,
            local,
            .add,
        );
        const positive = inputs.positive_neighbor_by_boundary_side[side];
        if (positive.runoff_mode == .subtract) {
            try updateRunoff(
                &runoff_water_heat,
                workspace.organic_by_class,
                &runoff_element,
                positive.flux.?,
                .subtract,
            );
        }
        if (positive.snow_mode == .subtract) {
            try updateSnow(
                &snow_water_heat,
                &snow_element,
                positive.flux.?,
                .subtract,
            );
        }
        if (side == 0) {
            if (inputs.opposite_neighbor_first_side) |opposite| {
                try updateWhole(
                    &runoff_water_heat,
                    &snow_water_heat,
                    workspace.organic_by_class,
                    &runoff_element,
                    &snow_element,
                    opposite,
                    .subtract,
                );
            }
        }
    }

    state.runoff_water_heat = runoff_water_heat;
    state.snow_water_heat = snow_water_heat;
    @memcpy(state.organic_by_class, workspace.organic_by_class);
    state.runoff_element = runoff_element;
    state.snow_element = snow_element;
}

fn validateDimensions(inputs: Inputs, state: State, workspace: Workspace) !void {
    if (inputs.organic_class_count == 0 or
        inputs.local_flux_by_boundary_side.len == 0)
    {
        return error.InvalidSurfaceAggregationDimensions;
    }
    const side_count = inputs.local_flux_by_boundary_side.len;
    if (inputs.positive_neighbor_by_boundary_side.len != side_count or
        state.organic_by_class.len != inputs.organic_class_count or
        workspace.organic_by_class.len != inputs.organic_class_count)
    {
        return error.SurfaceAggregationDimensionMismatch;
    }
    for (inputs.local_flux_by_boundary_side) |flux|
        if (flux.organic_by_class.len != inputs.organic_class_count)
            return error.SurfaceAggregationDimensionMismatch;
    for (inputs.positive_neighbor_by_boundary_side) |neighbor| {
        if (neighbor.flux) |flux| {
            if (flux.organic_by_class.len != inputs.organic_class_count)
                return error.SurfaceAggregationDimensionMismatch;
        } else if (neighbor.runoff_mode == .subtract or neighbor.snow_mode == .subtract) {
            return error.MissingPositiveNeighborFlux;
        }
    }
    if (inputs.opposite_neighbor_first_side) |flux|
        if (flux.organic_by_class.len != inputs.organic_class_count)
            return error.SurfaceAggregationDimensionMismatch;
}

fn validateInputs(inputs: Inputs, state: State, workspace: Workspace) !void {
    try validateFinite(state.runoff_water_heat);
    try validateFinite(state.snow_water_heat);
    try validateFinite(state.runoff_element);
    try validateFinite(state.snow_element);
    try validateFiniteSlice(state.organic_by_class);
    if (overlap(workspace.organic_by_class, state.organic_by_class))
        return error.SurfaceAggregationWorkspaceOverlap;

    for (inputs.local_flux_by_boundary_side) |flux| {
        try validateFlux(flux);
        if (overlap(workspace.organic_by_class, flux.organic_by_class))
            return error.SurfaceAggregationWorkspaceOverlap;
    }
    for (inputs.positive_neighbor_by_boundary_side) |neighbor| {
        if (neighbor.flux) |flux| {
            try validateFlux(flux);
            if (overlap(workspace.organic_by_class, flux.organic_by_class))
                return error.SurfaceAggregationWorkspaceOverlap;
        }
    }
    if (inputs.opposite_neighbor_first_side) |flux| {
        try validateFlux(flux);
        if (overlap(workspace.organic_by_class, flux.organic_by_class))
            return error.SurfaceAggregationWorkspaceOverlap;
    }
}

fn validateFlux(flux: SurfaceFlux) !void {
    try validateFinite(flux.runoff_water_heat);
    try validateFinite(flux.snow_water_heat);
    try validateFiniteSlice(flux.organic_by_class);
    try validateFinite(flux.runoff_element);
    try validateFinite(flux.snow_element);
}

fn updateWhole(
    runoff_water_heat: *RunoffWaterHeatFlux,
    snow_water_heat: *SnowWaterHeatFlux,
    organic_by_class: []OrganicFlux,
    runoff_element: *RunoffElementFlux,
    snow_element: *SnowElementFlux,
    flux: SurfaceFlux,
    operation: Operation,
) !void {
    try updateStruct(runoff_water_heat, flux.runoff_water_heat, operation);
    try updateStruct(snow_water_heat, flux.snow_water_heat, operation);
    for (organic_by_class, flux.organic_by_class) |*candidate, contribution|
        try updateStruct(candidate, contribution, operation);
    try updateStruct(runoff_element, flux.runoff_element, operation);
    try updateStruct(snow_element, flux.snow_element, operation);
}

fn updateRunoff(
    runoff_water_heat: *RunoffWaterHeatFlux,
    organic_by_class: []OrganicFlux,
    runoff_element: *RunoffElementFlux,
    flux: SurfaceFlux,
    operation: Operation,
) !void {
    try updateStruct(runoff_water_heat, flux.runoff_water_heat, operation);
    for (organic_by_class, flux.organic_by_class) |*candidate, contribution|
        try updateStruct(candidate, contribution, operation);
    try updateStruct(runoff_element, flux.runoff_element, operation);
}

fn updateSnow(
    snow_water_heat: *SnowWaterHeatFlux,
    snow_element: *SnowElementFlux,
    flux: SurfaceFlux,
    operation: Operation,
) !void {
    try updateStruct(snow_water_heat, flux.snow_water_heat, operation);
    try updateStruct(snow_element, flux.snow_element, operation);
}

fn updateStruct(candidate: anytype, contribution: anytype, operation: Operation) !void {
    inline for (@typeInfo(@TypeOf(candidate.*)).@"struct".fields) |field| {
        const current = @field(candidate.*, field.name);
        const change = @field(contribution, field.name);
        const result = switch (operation) {
            .add => current + change,
            .subtract => current - change,
        };
        if (!std.math.isFinite(result)) return error.NonFiniteSurfaceAggregationResult;
        @field(candidate.*, field.name) = result;
    }
}

fn validateFiniteSlice(values: []const OrganicFlux) !void {
    for (values) |value| try validateFinite(value);
}

fn validateFinite(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name)))
            return error.NonFiniteSurfaceAggregationInput;
}

fn overlap(left: []const OrganicFlux, right: []const OrganicFlux) bool {
    if (left.len == 0 or right.len == 0) return false;
    const item_size = @sizeOf(OrganicFlux);
    const left_start = @intFromPtr(left.ptr);
    const right_start = @intFromPtr(right.ptr);
    const left_bytes = std.math.mul(usize, left.len, item_size) catch return true;
    const right_bytes = std.math.mul(usize, right.len, item_size) catch return true;
    const left_end = std.math.add(usize, left_start, left_bytes) catch return true;
    const right_end = std.math.add(usize, right_start, right_bytes) catch return true;
    return left_start < right_end and right_start < left_end;
}

fn filled(comptime T: type, value: f64) T {
    var result: T = undefined;
    inline for (@typeInfo(T).@"struct".fields) |field|
        @field(result, field.name) = value;
    return result;
}

fn filledFlux(organic: []const OrganicFlux, value: f64) SurfaceFlux {
    return .{
        .runoff_water_heat = filled(RunoffWaterHeatFlux, value),
        .snow_water_heat = filled(SnowWaterHeatFlux, value),
        .organic_by_class = organic,
        .runoff_element = filled(RunoffElementFlux, value),
        .snow_element = filled(SnowElementFlux, value),
    };
}

fn expectStructValue(actual: anytype, expected: f64) !void {
    inline for (@typeInfo(@TypeOf(actual)).@"struct".fields) |field|
        try std.testing.expectEqual(expected, @field(actual, field.name));
}

fn expectOppositeSumZero(first: anytype, second: @TypeOf(first)) !void {
    inline for (@typeInfo(@TypeOf(first)).@"struct".fields) |field|
        try std.testing.expectEqual(
            @as(f64, 0),
            @field(first, field.name) + @field(second, field.name),
        );
}

test "source side order and independent runoff snow gates are exact" {
    const organic_count = 3;
    const one = [_]OrganicFlux{filled(OrganicFlux, 1)} ** organic_count;
    const two = [_]OrganicFlux{filled(OrganicFlux, 2)} ** organic_count;
    const three = [_]OrganicFlux{filled(OrganicFlux, 3)} ** organic_count;
    const four = [_]OrganicFlux{filled(OrganicFlux, 4)} ** organic_count;
    const five = [_]OrganicFlux{filled(OrganicFlux, 5)} ** organic_count;
    const local = [_]SurfaceFlux{ filledFlux(&one, 1), filledFlux(&four, 4) };
    const positive = [_]PositiveNeighbor{
        .{ .flux = filledFlux(&two, 2), .runoff_mode = .subtract, .snow_mode = .blocked },
        .{ .flux = filledFlux(&five, 5), .runoff_mode = .blocked, .snow_mode = .subtract },
    };
    var organic = [_]OrganicFlux{filled(OrganicFlux, 100)} ** organic_count;
    var scratch: [organic_count]OrganicFlux = undefined;
    var state = State{
        .runoff_water_heat = filled(RunoffWaterHeatFlux, 100),
        .snow_water_heat = filled(SnowWaterHeatFlux, 100),
        .organic_by_class = &organic,
        .runoff_element = filled(RunoffElementFlux, 100),
        .snow_element = filled(SnowElementFlux, 100),
    };

    try aggregate(.{
        .organic_class_count = organic_count,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_by_boundary_side = &positive,
        .opposite_neighbor_first_side = filledFlux(&three, 3),
    }, &state, .{ .organic_by_class = &scratch });

    try expectStructValue(state.runoff_water_heat, 100);
    try expectStructValue(state.snow_water_heat, 97);
    try expectStructValue(state.runoff_element, 100);
    try expectStructValue(state.snow_element, 97);
    for (state.organic_by_class) |value| try expectStructValue(value, 100);
}

test "shared horizontal face conserves every transported inventory exactly" {
    const organic_count = 2;
    const seven = [_]OrganicFlux{filled(OrganicFlux, 7)} ** organic_count;
    const zero = [_]OrganicFlux{filled(OrganicFlux, 0)} ** organic_count;
    const first_local = [_]SurfaceFlux{filledFlux(&seven, 7)};
    const second_local = [_]SurfaceFlux{filledFlux(&zero, 0)};
    const blocked = [_]PositiveNeighbor{.{
        .flux = null,
        .runoff_mode = .blocked,
        .snow_mode = .blocked,
    }};
    var first_organic = zero;
    var second_organic = zero;
    var first_scratch: [organic_count]OrganicFlux = undefined;
    var second_scratch: [organic_count]OrganicFlux = undefined;
    var first = State{
        .runoff_water_heat = .{},
        .snow_water_heat = .{},
        .organic_by_class = &first_organic,
        .runoff_element = .{},
        .snow_element = .{},
    };
    var second = State{
        .runoff_water_heat = .{},
        .snow_water_heat = .{},
        .organic_by_class = &second_organic,
        .runoff_element = .{},
        .snow_element = .{},
    };

    try aggregate(.{
        .organic_class_count = organic_count,
        .local_flux_by_boundary_side = &first_local,
        .positive_neighbor_by_boundary_side = &blocked,
        .opposite_neighbor_first_side = null,
    }, &first, .{ .organic_by_class = &first_scratch });
    try aggregate(.{
        .organic_class_count = organic_count,
        .local_flux_by_boundary_side = &second_local,
        .positive_neighbor_by_boundary_side = &blocked,
        .opposite_neighbor_first_side = filledFlux(&seven, 7),
    }, &second, .{ .organic_by_class = &second_scratch });

    try expectOppositeSumZero(first.runoff_water_heat, second.runoff_water_heat);
    try expectOppositeSumZero(first.snow_water_heat, second.snow_water_heat);
    try expectOppositeSumZero(first.runoff_element, second.runoff_element);
    try expectOppositeSumZero(first.snow_element, second.snow_element);
    for (first.organic_by_class, second.organic_by_class) |a, b|
        try expectOppositeSumZero(a, b);
}

test "runtime organic extent and source scalar counts are explicit" {
    try std.testing.expectEqual(
        @as(usize, 2),
        @typeInfo(RunoffWaterHeatFlux).@"struct".fields.len,
    );
    try std.testing.expectEqual(
        @as(usize, 4),
        @typeInfo(SnowWaterHeatFlux).@"struct".fields.len,
    );
    try std.testing.expectEqual(
        @as(usize, 12),
        @typeInfo(RunoffElementFlux).@"struct".fields.len,
    );
    try std.testing.expectEqual(
        @as(usize, 10),
        @typeInfo(SnowElementFlux).@"struct".fields.len,
    );
    const organic = [_]OrganicFlux{filled(OrganicFlux, 1)} ** 7;
    const local = [_]SurfaceFlux{filledFlux(&organic, 1)};
    const blocked = [_]PositiveNeighbor{.{
        .flux = null,
        .runoff_mode = .blocked,
        .snow_mode = .blocked,
    }};
    var state_organic = [_]OrganicFlux{filled(OrganicFlux, 0)} ** 7;
    var scratch: [7]OrganicFlux = undefined;
    var state = State{
        .runoff_water_heat = .{},
        .snow_water_heat = .{},
        .organic_by_class = &state_organic,
        .runoff_element = .{},
        .snow_element = .{},
    };
    try aggregate(.{
        .organic_class_count = 7,
        .local_flux_by_boundary_side = &local,
        .positive_neighbor_by_boundary_side = &blocked,
        .opposite_neighbor_first_side = null,
    }, &state, .{ .organic_by_class = &scratch });
    for (state.organic_by_class) |value| try expectStructValue(value, 1);
}

test "late invalid input and arithmetic overflow preserve state atomically" {
    var local_organic = [_]OrganicFlux{ filled(OrganicFlux, 1), filled(OrganicFlux, 1) };
    local_organic[1].acetate_carbon_g_c_per_step = std.math.nan(f64);
    const local = [_]SurfaceFlux{filledFlux(&local_organic, 1)};
    const blocked = [_]PositiveNeighbor{.{
        .flux = null,
        .runoff_mode = .blocked,
        .snow_mode = .blocked,
    }};
    var organic = [_]OrganicFlux{filled(OrganicFlux, 3)} ** 2;
    var scratch = [_]OrganicFlux{filled(OrganicFlux, 9)} ** 2;
    var state = State{
        .runoff_water_heat = filled(RunoffWaterHeatFlux, 3),
        .snow_water_heat = filled(SnowWaterHeatFlux, 3),
        .organic_by_class = &organic,
        .runoff_element = filled(RunoffElementFlux, 3),
        .snow_element = filled(SnowElementFlux, 3),
    };
    try std.testing.expectError(
        error.NonFiniteSurfaceAggregationInput,
        aggregate(.{
            .organic_class_count = 2,
            .local_flux_by_boundary_side = &local,
            .positive_neighbor_by_boundary_side = &blocked,
            .opposite_neighbor_first_side = null,
        }, &state, .{ .organic_by_class = &scratch }),
    );
    try expectStructValue(state.runoff_water_heat, 3);
    for (state.organic_by_class) |value| try expectStructValue(value, 3);

    local_organic[1] = filled(OrganicFlux, std.math.floatMax(f64));
    local_organic[0] = filled(OrganicFlux, std.math.floatMax(f64));
    const overflowing_flux = filledFlux(&local_organic, std.math.floatMax(f64));
    const overflowing_local = [_]SurfaceFlux{overflowing_flux};
    state.runoff_water_heat = filled(RunoffWaterHeatFlux, std.math.floatMax(f64));
    try std.testing.expectError(
        error.NonFiniteSurfaceAggregationResult,
        aggregate(.{
            .organic_class_count = 2,
            .local_flux_by_boundary_side = &overflowing_local,
            .positive_neighbor_by_boundary_side = &blocked,
            .opposite_neighbor_first_side = null,
        }, &state, .{ .organic_by_class = &scratch }),
    );
    try expectStructValue(state.runoff_water_heat, std.math.floatMax(f64));
}

test "dimension and workspace alias errors precede mutation" {
    const organic = [_]OrganicFlux{filled(OrganicFlux, 1)} ** 2;
    const local = [_]SurfaceFlux{filledFlux(&organic, 1)};
    const blocked = [_]PositiveNeighbor{.{
        .flux = null,
        .runoff_mode = .blocked,
        .snow_mode = .blocked,
    }};
    var state_organic = [_]OrganicFlux{filled(OrganicFlux, 4)} ** 2;
    var state = State{
        .runoff_water_heat = filled(RunoffWaterHeatFlux, 4),
        .snow_water_heat = filled(SnowWaterHeatFlux, 4),
        .organic_by_class = &state_organic,
        .runoff_element = filled(RunoffElementFlux, 4),
        .snow_element = filled(SnowElementFlux, 4),
    };
    try std.testing.expectError(
        error.SurfaceAggregationWorkspaceOverlap,
        aggregate(.{
            .organic_class_count = 2,
            .local_flux_by_boundary_side = &local,
            .positive_neighbor_by_boundary_side = &blocked,
            .opposite_neighbor_first_side = null,
        }, &state, .{ .organic_by_class = &state_organic }),
    );
    try expectStructValue(state.runoff_water_heat, 4);

    var short_scratch = [_]OrganicFlux{filled(OrganicFlux, 0)};
    try std.testing.expectError(
        error.SurfaceAggregationDimensionMismatch,
        aggregate(.{
            .organic_class_count = 2,
            .local_flux_by_boundary_side = &local,
            .positive_neighbor_by_boundary_side = &blocked,
            .opposite_neighbor_first_side = null,
        }, &state, .{ .organic_by_class = &short_scratch }),
    );
    try expectStructValue(state.runoff_water_heat, 4);
}
