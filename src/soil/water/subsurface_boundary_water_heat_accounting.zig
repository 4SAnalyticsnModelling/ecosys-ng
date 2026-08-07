const std = @import("std");

pub const HorizontalSubsurfaceExchange = enum {
    enabled,
    standalone,
};

pub const ExternalBoundaryWindow = struct {
    first_column: usize,
    last_column_inclusive: usize,
    first_row: usize,
    last_row_inclusive: usize,
};

/// REDIST evaluates axes and sides in this order for every soil layer.
pub const Face = enum(u8) {
    east,
    west,
    south,
    north,
    lower,

    pub const count: usize = @typeInfo(Face).@"enum".fields.len;
};

pub const BoundaryFlux = struct {
    micropore_water_m3_per_step: f64 = 0,
    macropore_water_m3_per_step: f64 = 0,
    convective_heat_megajoules_per_step: f64 = 0,
    artificial_drainage_micropore_water_m3_per_step: f64 = 0,
    artificial_drainage_macropore_water_m3_per_step: f64 = 0,
};

pub const Inputs = struct {
    lon_count: usize,
    lat_count: usize,
    soil_layer_count: usize,
    external_boundary_window: ExternalBoundaryWindow,
    horizontal_exchange_by_cell: []const HorizontalSubsurfaceExchange,
    /// Each face slice is flattened `(cell, soil_layer)`.
    flux_by_face: [Face.count][]const BoundaryFlux,
};

pub const State = struct {
    cumulative_outward_water_m3: f64,
    cumulative_outward_heat_megajoules: f64,
    net_outward_water_m3_per_step_by_cell: []f64,
    cumulative_outward_water_m3_by_cell: []f64,
    artificial_drainage_outward_m3_per_step_by_cell: []f64,
};

/// Accounts for subsurface water, heat, and artificial-drainage boundaries.
///
/// Traceability: REDIST.F lines 1152--1160 (`HFLW`, `WO`, `VOLWOU`,
/// `HVOLO`, `UVOLO`, `WY`, and `UVOLY`). Traversal retains source order:
/// column, row, soil layer, east, west, south, north, then the lower face at
/// the last layer. A standalone cell suppresses horizontal exchange while its
/// lower boundary remains active, matching `NCN == 3` with `N == 3`.
///
/// Heat is accumulated before and independently of the exact `WO != 0` water
/// gate. Artificial drainage remains inside that gate, including when its own
/// flux is nonzero but micropore and macropore ordinary water cancel.
pub fn apply(inputs: Inputs, state: *State) !void {
    const dimensions = try validateDimensions(inputs, state.*);
    try validateInputsAndState(inputs, state.*, dimensions);
    try preflightUpdates(inputs, state.*);

    const window = inputs.external_boundary_window;
    for (window.first_column..window.last_column_inclusive + 1) |column| {
        for (window.first_row..window.last_row_inclusive + 1) |row| {
            const cell = row * inputs.lon_count + column;
            for (0..inputs.soil_layer_count) |layer| {
                const record = cell * inputs.soil_layer_count + layer;
                for (std.meta.tags(Face)) |face| {
                    if (!isActiveBoundary(inputs, face, column, row, cell, layer))
                        continue;
                    const flux = inputs.flux_by_face[@intFromEnum(face)][record];
                    const direction = inwardSign(face);
                    const signed_heat_megajoules =
                        checkedProduct(direction, flux.convective_heat_megajoules_per_step) catch
                            unreachable;
                    state.cumulative_outward_heat_megajoules =
                        checkedSum(
                            state.cumulative_outward_heat_megajoules,
                            -signed_heat_megajoules,
                        ) catch unreachable;
                    const signed_water_m3 =
                        calculateSignedOrdinaryWater(flux, direction) catch unreachable;
                    if (signed_water_m3 == 0) continue;
                    state.cumulative_outward_water_m3 =
                        checkedSum(
                            state.cumulative_outward_water_m3,
                            -signed_water_m3,
                        ) catch unreachable;
                    state.net_outward_water_m3_per_step_by_cell[cell] =
                        checkedSum(
                            state.net_outward_water_m3_per_step_by_cell[cell],
                            -signed_water_m3,
                        ) catch unreachable;
                    state.cumulative_outward_water_m3_by_cell[cell] =
                        checkedSum(
                            state.cumulative_outward_water_m3_by_cell[cell],
                            -signed_water_m3,
                        ) catch unreachable;
                    const signed_drainage_m3 =
                        calculateSignedArtificialDrainage(flux, direction) catch
                            unreachable;
                    state.artificial_drainage_outward_m3_per_step_by_cell[cell] =
                        checkedSum(
                            state.artificial_drainage_outward_m3_per_step_by_cell[cell],
                            -signed_drainage_m3,
                        ) catch unreachable;
                }
            }
        }
    }
}

const Dimensions = struct {
    cell_count: usize,
    record_count: usize,
};

fn preflightUpdates(inputs: Inputs, state: State) !void {
    var cumulative_outward_water_m3 = state.cumulative_outward_water_m3;
    var cumulative_outward_heat_megajoules = state.cumulative_outward_heat_megajoules;
    const window = inputs.external_boundary_window;
    for (window.first_column..window.last_column_inclusive + 1) |column| {
        for (window.first_row..window.last_row_inclusive + 1) |row| {
            const cell = row * inputs.lon_count + column;
            var net_outward_water_m3 =
                state.net_outward_water_m3_per_step_by_cell[cell];
            var cell_cumulative_outward_water_m3 =
                state.cumulative_outward_water_m3_by_cell[cell];
            var artificial_drainage_outward_m3 =
                state.artificial_drainage_outward_m3_per_step_by_cell[cell];
            for (0..inputs.soil_layer_count) |layer| {
                const record = cell * inputs.soil_layer_count + layer;
                for (std.meta.tags(Face)) |face| {
                    if (!isActiveBoundary(inputs, face, column, row, cell, layer))
                        continue;
                    const flux = inputs.flux_by_face[@intFromEnum(face)][record];
                    const direction = inwardSign(face);
                    const signed_heat_megajoules =
                        try checkedProduct(direction, flux.convective_heat_megajoules_per_step);
                    cumulative_outward_heat_megajoules = try checkedSum(
                        cumulative_outward_heat_megajoules,
                        -signed_heat_megajoules,
                    );
                    const signed_water_m3 =
                        try calculateSignedOrdinaryWater(flux, direction);
                    if (signed_water_m3 == 0) continue;
                    cumulative_outward_water_m3 = try checkedSum(
                        cumulative_outward_water_m3,
                        -signed_water_m3,
                    );
                    net_outward_water_m3 =
                        try checkedSum(net_outward_water_m3, -signed_water_m3);
                    cell_cumulative_outward_water_m3 = try checkedSum(
                        cell_cumulative_outward_water_m3,
                        -signed_water_m3,
                    );
                    const signed_drainage_m3 =
                        try calculateSignedArtificialDrainage(flux, direction);
                    artificial_drainage_outward_m3 = try checkedSum(
                        artificial_drainage_outward_m3,
                        -signed_drainage_m3,
                    );
                }
            }
        }
    }
}

fn calculateSignedOrdinaryWater(
    flux: BoundaryFlux,
    direction: f64,
) !f64 {
    return checkedProduct(
        direction,
        try checkedSum(
            flux.micropore_water_m3_per_step,
            flux.macropore_water_m3_per_step,
        ),
    );
}

fn calculateSignedArtificialDrainage(
    flux: BoundaryFlux,
    direction: f64,
) !f64 {
    return checkedProduct(
        direction,
        try checkedSum(
            flux.artificial_drainage_micropore_water_m3_per_step,
            flux.artificial_drainage_macropore_water_m3_per_step,
        ),
    );
}

fn validateDimensions(inputs: Inputs, state: State) !Dimensions {
    if (inputs.lon_count == 0 or
        inputs.lat_count == 0 or
        inputs.soil_layer_count == 0)
        return error.InvalidSubsurfaceBoundaryDimensions;
    const cell_count = std.math.mul(
        usize,
        inputs.lon_count,
        inputs.lat_count,
    ) catch return error.InvalidSubsurfaceBoundaryDimensions;
    const record_count = std.math.mul(
        usize,
        cell_count,
        inputs.soil_layer_count,
    ) catch return error.InvalidSubsurfaceBoundaryDimensions;
    inline for (inputs.flux_by_face) |values|
        if (values.len != record_count)
            return error.InvalidSubsurfaceBoundaryDimensions;
    if (inputs.horizontal_exchange_by_cell.len != cell_count or
        state.net_outward_water_m3_per_step_by_cell.len != cell_count or
        state.cumulative_outward_water_m3_by_cell.len != cell_count or
        state.artificial_drainage_outward_m3_per_step_by_cell.len != cell_count)
        return error.InvalidSubsurfaceBoundaryDimensions;
    const window = inputs.external_boundary_window;
    if (window.first_column > window.last_column_inclusive or
        window.first_row > window.last_row_inclusive or
        window.last_column_inclusive >= inputs.lon_count or
        window.last_row_inclusive >= inputs.lat_count)
        return error.InvalidSubsurfaceBoundaryWindow;
    return .{ .cell_count = cell_count, .record_count = record_count };
}

fn validateInputsAndState(
    inputs: Inputs,
    state: State,
    dimensions: Dimensions,
) !void {
    inline for (inputs.flux_by_face) |values|
        for (values) |flux|
            inline for (@typeInfo(BoundaryFlux).@"struct".fields) |field|
                if (!std.math.isFinite(@field(flux, field.name)))
                    return error.InvalidSubsurfaceBoundaryInput;
    if (!std.math.isFinite(state.cumulative_outward_water_m3) or
        !std.math.isFinite(state.cumulative_outward_heat_megajoules))
        return error.InvalidSubsurfaceBoundaryState;
    for (0..dimensions.cell_count) |cell|
        inline for (.{
            state.net_outward_water_m3_per_step_by_cell[cell],
            state.cumulative_outward_water_m3_by_cell[cell],
            state.artificial_drainage_outward_m3_per_step_by_cell[cell],
        }) |value|
            if (!std.math.isFinite(value))
                return error.InvalidSubsurfaceBoundaryState;
}

fn isActiveBoundary(
    inputs: Inputs,
    face: Face,
    column: usize,
    row: usize,
    cell: usize,
    layer: usize,
) bool {
    if (face == .lower)
        return layer + 1 == inputs.soil_layer_count;
    if (inputs.horizontal_exchange_by_cell[cell] == .standalone)
        return false;
    const window = inputs.external_boundary_window;
    return switch (face) {
        .east => column == window.last_column_inclusive,
        .west => column == window.first_column,
        .south => row == window.last_row_inclusive,
        .north => row == window.first_row,
        .lower => unreachable,
    };
}

fn inwardSign(face: Face) f64 {
    return switch (face) {
        .east, .south, .lower => -1,
        .west, .north => 1,
    };
}

fn checkedProduct(left: f64, right: f64) !f64 {
    const result = left * right;
    if (!std.math.isFinite(result))
        return error.NonFiniteSubsurfaceBoundaryResult;
    return result;
}

fn checkedSum(current: f64, increment: f64) !f64 {
    const result = current + increment;
    if (!std.math.isFinite(result))
        return error.NonFiniteSubsurfaceBoundaryResult;
    return result;
}

fn zeroState(
    net: []f64,
    cumulative: []f64,
    drainage: []f64,
) State {
    return .{
        .cumulative_outward_water_m3 = 0,
        .cumulative_outward_heat_megajoules = 0,
        .net_outward_water_m3_per_step_by_cell = net,
        .cumulative_outward_water_m3_by_cell = cumulative,
        .artificial_drainage_outward_m3_per_step_by_cell = drainage,
    };
}

test "REDIST subsurface boundary accounts layers and lower face in source order" {
    const empty = [_]BoundaryFlux{ .{}, .{} };
    const east = [_]BoundaryFlux{
        .{
            .micropore_water_m3_per_step = 1,
            .macropore_water_m3_per_step = 2,
            .convective_heat_megajoules_per_step = 3,
            .artificial_drainage_micropore_water_m3_per_step = 4,
            .artificial_drainage_macropore_water_m3_per_step = 5,
        },
        .{},
    };
    const lower = [_]BoundaryFlux{
        .{},
        .{
            .micropore_water_m3_per_step = 6,
            .macropore_water_m3_per_step = 7,
            .convective_heat_megajoules_per_step = 8,
            .artificial_drainage_micropore_water_m3_per_step = 9,
            .artificial_drainage_macropore_water_m3_per_step = 10,
        },
    };
    var net = [_]f64{0};
    var cumulative = [_]f64{0};
    var drainage = [_]f64{0};
    var state = zeroState(&net, &cumulative, &drainage);
    try apply(.{
        .lon_count = 1,
        .lat_count = 1,
        .soil_layer_count = 2,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .horizontal_exchange_by_cell = &.{.enabled},
        .flux_by_face = .{ &east, &empty, &empty, &empty, &lower },
    }, &state);

    try std.testing.expectEqual(@as(f64, 16), state.cumulative_outward_water_m3);
    try std.testing.expectEqual(@as(f64, 11), state.cumulative_outward_heat_megajoules);
    try std.testing.expectEqual(@as(f64, 16), net[0]);
    try std.testing.expectEqual(@as(f64, 16), cumulative[0]);
    try std.testing.expectEqual(@as(f64, 28), drainage[0]);
}

test "standalone cell suppresses lateral exchange but retains lower drainage" {
    const lateral = [_]BoundaryFlux{
        .{ .micropore_water_m3_per_step = 20, .convective_heat_megajoules_per_step = 20 },
    };
    const empty = [_]BoundaryFlux{.{}};
    const lower = [_]BoundaryFlux{
        .{ .micropore_water_m3_per_step = 2, .convective_heat_megajoules_per_step = 3 },
    };
    var net = [_]f64{0};
    var cumulative = [_]f64{0};
    var drainage = [_]f64{0};
    var state = zeroState(&net, &cumulative, &drainage);
    try apply(.{
        .lon_count = 1,
        .lat_count = 1,
        .soil_layer_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .horizontal_exchange_by_cell = &.{.standalone},
        .flux_by_face = .{ &lateral, &empty, &empty, &empty, &lower },
    }, &state);
    try std.testing.expectEqual(@as(f64, 2), state.cumulative_outward_water_m3);
    try std.testing.expectEqual(@as(f64, 3), state.cumulative_outward_heat_megajoules);
}

test "heat advances while exact ordinary-water cancellation gates drainage" {
    const cancelling = [_]BoundaryFlux{.{
        .micropore_water_m3_per_step = 1,
        .macropore_water_m3_per_step = -1,
        .convective_heat_megajoules_per_step = 2,
        .artificial_drainage_micropore_water_m3_per_step = 5,
        .artificial_drainage_macropore_water_m3_per_step = 6,
    }};
    const empty = [_]BoundaryFlux{.{}};
    var net = [_]f64{0};
    var cumulative = [_]f64{0};
    var drainage = [_]f64{0};
    var state = zeroState(&net, &cumulative, &drainage);
    try apply(.{
        .lon_count = 1,
        .lat_count = 1,
        .soil_layer_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .horizontal_exchange_by_cell = &.{.enabled},
        .flux_by_face = .{ &cancelling, &empty, &empty, &empty, &empty },
    }, &state);
    try std.testing.expectEqual(@as(f64, 0), state.cumulative_outward_water_m3);
    try std.testing.expectEqual(@as(f64, 2), state.cumulative_outward_heat_megajoules);
    try std.testing.expectEqual(@as(f64, 0), drainage[0]);
}

test "runtime grid water ledger closes against per-cell cumulative values" {
    const active: BoundaryFlux = .{
        .micropore_water_m3_per_step = 1,
        .convective_heat_megajoules_per_step = 2,
    };
    const east = [_]BoundaryFlux{
        .{}, .{}, .{},    active,
        .{}, .{}, active, .{},
    };
    const empty = [_]BoundaryFlux{.{}} ** 8;
    var net = [_]f64{ 0, 0, 0, 0 };
    var cumulative = [_]f64{ 0, 0, 0, 0 };
    var drainage = [_]f64{ 0, 0, 0, 0 };
    var state = zeroState(&net, &cumulative, &drainage);
    try apply(.{
        .lon_count = 2,
        .lat_count = 2,
        .soil_layer_count = 2,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 1,
            .first_row = 0,
            .last_row_inclusive = 1,
        },
        .horizontal_exchange_by_cell = &.{
            .enabled, .enabled, .enabled, .enabled,
        },
        .flux_by_face = .{ &east, &empty, &empty, &empty, &empty },
    }, &state);
    var cell_total_m3: f64 = 0;
    for (cumulative) |value| cell_total_m3 += value;
    try std.testing.expectEqual(state.cumulative_outward_water_m3, cell_total_m3);
    try std.testing.expectEqualSlices(f64, &.{ 0, 1, 0, 1 }, &cumulative);
    try std.testing.expectEqual(@as(f64, 4), state.cumulative_outward_heat_megajoules);
}

test "late non-finite boundary flux leaves every ledger unchanged" {
    var east = [_]BoundaryFlux{ .{}, .{} };
    east[1].convective_heat_megajoules_per_step = std.math.nan(f64);
    const empty = [_]BoundaryFlux{ .{}, .{} };
    var net = [_]f64{ 1, 2 };
    var cumulative = [_]f64{ 3, 4 };
    var drainage = [_]f64{ 5, 6 };
    var state: State = .{
        .cumulative_outward_water_m3 = 7,
        .cumulative_outward_heat_megajoules = 8,
        .net_outward_water_m3_per_step_by_cell = &net,
        .cumulative_outward_water_m3_by_cell = &cumulative,
        .artificial_drainage_outward_m3_per_step_by_cell = &drainage,
    };
    try std.testing.expectError(error.InvalidSubsurfaceBoundaryInput, apply(.{
        .lon_count = 2,
        .lat_count = 1,
        .soil_layer_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 1,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .horizontal_exchange_by_cell = &.{ .enabled, .enabled },
        .flux_by_face = .{ &east, &empty, &empty, &empty, &empty },
    }, &state));
    try std.testing.expectEqual(@as(f64, 7), state.cumulative_outward_water_m3);
    try std.testing.expectEqual(@as(f64, 8), state.cumulative_outward_heat_megajoules);
    try std.testing.expectEqualSlices(f64, &.{ 1, 2 }, &net);
    try std.testing.expectEqualSlices(f64, &.{ 3, 4 }, &cumulative);
    try std.testing.expectEqualSlices(f64, &.{ 5, 6 }, &drainage);
}

test "invalid runtime layer extent and boundary window fail explicitly" {
    const empty = [_]BoundaryFlux{.{}};
    var net = [_]f64{0};
    var cumulative = [_]f64{0};
    var drainage = [_]f64{0};
    var state = zeroState(&net, &cumulative, &drainage);
    const base: Inputs = .{
        .lon_count = 1,
        .lat_count = 1,
        .soil_layer_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .horizontal_exchange_by_cell = &.{.enabled},
        .flux_by_face = .{ &empty, &empty, &empty, &empty, &empty },
    };
    var invalid_layers = base;
    invalid_layers.soil_layer_count = 2;
    try std.testing.expectError(
        error.InvalidSubsurfaceBoundaryDimensions,
        apply(invalid_layers, &state),
    );
    var invalid_window = base;
    invalid_window.external_boundary_window.first_column = 1;
    try std.testing.expectError(
        error.InvalidSubsurfaceBoundaryWindow,
        apply(invalid_window, &state),
    );
}
