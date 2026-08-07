const std = @import("std");

pub const ExternalBoundaryWindow = struct {
    first_column: usize,
    last_column_inclusive: usize,
    first_row: usize,
    last_row_inclusive: usize,
};

/// REDIST evaluates these face records in this exact order.
pub const FaceFlux = enum(u8) {
    east,
    west,
    south,
    north,

    pub const count: usize = @typeInfo(FaceFlux).@"enum".fields.len;
};

pub const Inputs = struct {
    lon_count: usize,
    lat_count: usize,
    external_boundary_window: ExternalBoundaryWindow,
    water_m3_per_step_by_face: [FaceFlux.count][]const f64,
    heat_megajoules_per_step_by_face: [FaceFlux.count][]const f64,
    negligible_water_m3_by_cell: []const f64,
};

pub const State = struct {
    signed_runoff_water_change_m3_by_cell: []f64,
    cumulative_outward_runoff_m3_by_cell: []f64,
    cumulative_outward_runoff_m3: f64,
    cumulative_outward_runoff_heat_megajoules: f64,
};

/// Accounts for surface-runoff water and associated heat at external faces.
///
/// Traceability: REDIST.F lines 618--704 (`QR`, `HQR`, `WQRN`, `WQRH`,
/// `CRUN`, `URUN`, and `HEATOU`). The source traverses column, row, layer,
/// axis, then side. Because runoff is admitted only at the top layer and only
/// on horizontal axes, the equivalent retained order is column, row, east,
/// west, south, then north. Positive signed cell change is inward; cumulative
/// outward ledgers use the opposite sign. Heat is accounted only when the
/// corresponding water flux exceeds the cell's source threshold.
pub fn apply(inputs: Inputs, state: *State) !void {
    const cell_count = try validateDimensions(inputs, state.*);
    try validateInputsAndState(inputs, state.*, cell_count);
    try preflightUpdates(inputs, state.*);

    const window = inputs.external_boundary_window;
    for (window.first_column..window.last_column_inclusive + 1) |column| {
        for (window.first_row..window.last_row_inclusive + 1) |row| {
            const cell = row * inputs.lon_count + column;
            for (std.meta.tags(FaceFlux)) |face| {
                if (!isExternalFace(face, column, row, window)) continue;
                const face_index = @intFromEnum(face);
                const direction = inwardSign(face);
                const signed_water_change_m3 =
                    direction * inputs.water_m3_per_step_by_face[face_index][cell];
                state.signed_runoff_water_change_m3_by_cell[cell] +=
                    signed_water_change_m3;
                if (@abs(signed_water_change_m3) <=
                    inputs.negligible_water_m3_by_cell[cell]) continue;
                state.cumulative_outward_runoff_m3 -= signed_water_change_m3;
                state.cumulative_outward_runoff_m3_by_cell[cell] -=
                    signed_water_change_m3;
                const signed_heat_change_megajoules =
                    direction * inputs.heat_megajoules_per_step_by_face[face_index][cell];
                state.cumulative_outward_runoff_heat_megajoules -= signed_heat_change_megajoules;
            }
        }
    }
}

fn validateDimensions(inputs: Inputs, state: State) !usize {
    if (inputs.lon_count == 0 or inputs.lat_count == 0)
        return error.InvalidRunoffBoundaryDimensions;
    const cell_count = std.math.mul(
        usize,
        inputs.lon_count,
        inputs.lat_count,
    ) catch return error.InvalidRunoffBoundaryDimensions;
    inline for (inputs.water_m3_per_step_by_face ++ inputs.heat_megajoules_per_step_by_face) |values|
        if (values.len != cell_count)
            return error.InvalidRunoffBoundaryDimensions;
    if (inputs.negligible_water_m3_by_cell.len != cell_count or
        state.signed_runoff_water_change_m3_by_cell.len != cell_count or
        state.cumulative_outward_runoff_m3_by_cell.len != cell_count)
        return error.InvalidRunoffBoundaryDimensions;
    const window = inputs.external_boundary_window;
    if (window.first_column > window.last_column_inclusive or
        window.first_row > window.last_row_inclusive or
        window.last_column_inclusive >= inputs.lon_count or
        window.last_row_inclusive >= inputs.lat_count)
        return error.InvalidRunoffBoundaryWindow;
    return cell_count;
}

fn validateInputsAndState(inputs: Inputs, state: State, cell_count: usize) !void {
    inline for (inputs.water_m3_per_step_by_face) |values|
        for (values) |value|
            if (!std.math.isFinite(value))
                return error.InvalidRunoffBoundaryInput;
    inline for (inputs.heat_megajoules_per_step_by_face) |values|
        for (values) |value|
            if (!std.math.isFinite(value))
                return error.InvalidRunoffBoundaryInput;
    for (0..cell_count) |cell| {
        if (!nonnegativeFinite(inputs.negligible_water_m3_by_cell[cell]) or
            !std.math.isFinite(state.signed_runoff_water_change_m3_by_cell[cell]) or
            !std.math.isFinite(state.cumulative_outward_runoff_m3_by_cell[cell]))
            return error.InvalidRunoffBoundaryState;
    }
    if (!std.math.isFinite(state.cumulative_outward_runoff_m3) or
        !std.math.isFinite(state.cumulative_outward_runoff_heat_megajoules))
        return error.InvalidRunoffBoundaryState;
}

fn preflightUpdates(inputs: Inputs, state: State) !void {
    var cumulative_outward_runoff_m3 = state.cumulative_outward_runoff_m3;
    var cumulative_outward_runoff_heat_megajoules =
        state.cumulative_outward_runoff_heat_megajoules;
    const window = inputs.external_boundary_window;
    for (window.first_column..window.last_column_inclusive + 1) |column| {
        for (window.first_row..window.last_row_inclusive + 1) |row| {
            const cell = row * inputs.lon_count + column;
            var signed_runoff_water_change_m3 =
                state.signed_runoff_water_change_m3_by_cell[cell];
            var cumulative_outward_runoff_m3_for_cell =
                state.cumulative_outward_runoff_m3_by_cell[cell];
            for (std.meta.tags(FaceFlux)) |face| {
                if (!isExternalFace(face, column, row, window)) continue;
                const face_index = @intFromEnum(face);
                const direction = inwardSign(face);
                const signed_water_change_m3 =
                    direction * inputs.water_m3_per_step_by_face[face_index][cell];
                signed_runoff_water_change_m3 = try checkedSum(
                    signed_runoff_water_change_m3,
                    signed_water_change_m3,
                );
                if (@abs(signed_water_change_m3) <=
                    inputs.negligible_water_m3_by_cell[cell]) continue;
                cumulative_outward_runoff_m3 = try checkedSum(
                    cumulative_outward_runoff_m3,
                    -signed_water_change_m3,
                );
                cumulative_outward_runoff_m3_for_cell = try checkedSum(
                    cumulative_outward_runoff_m3_for_cell,
                    -signed_water_change_m3,
                );
                const signed_heat_change_megajoules =
                    direction * inputs.heat_megajoules_per_step_by_face[face_index][cell];
                cumulative_outward_runoff_heat_megajoules = try checkedSum(
                    cumulative_outward_runoff_heat_megajoules,
                    -signed_heat_change_megajoules,
                );
            }
        }
    }
}

fn isExternalFace(
    face: FaceFlux,
    column: usize,
    row: usize,
    window: ExternalBoundaryWindow,
) bool {
    return switch (face) {
        .east => column == window.last_column_inclusive,
        .west => column == window.first_column,
        .south => row == window.last_row_inclusive,
        .north => row == window.first_row,
    };
}

fn inwardSign(face: FaceFlux) f64 {
    return switch (face) {
        .east, .south => -1,
        .west, .north => 1,
    };
}

fn checkedSum(current: f64, increment: f64) !f64 {
    const result = current + increment;
    if (!std.math.isFinite(result))
        return error.NonFiniteRunoffBoundaryResult;
    return result;
}

fn nonnegativeFinite(value: f64) bool {
    return std.math.isFinite(value) and value >= 0;
}

test "REDIST runoff boundary water and heat preserve signed conservation" {
    const east_water = [_]f64{ 0, 1, 0, 1 };
    const west_water = [_]f64{ 2, 0, 2, 0 };
    const south_water = [_]f64{ 0, 0, 3, 3 };
    const north_water = [_]f64{ 4, 4, 0, 0 };
    const east_heat = [_]f64{ 0, 10, 0, 10 };
    const west_heat = [_]f64{ 20, 0, 20, 0 };
    const south_heat = [_]f64{ 0, 0, 30, 30 };
    const north_heat = [_]f64{ 40, 40, 0, 0 };
    var signed_water = [_]f64{0} ** 4;
    var cell_outward_water = [_]f64{0} ** 4;
    var state: State = .{
        .signed_runoff_water_change_m3_by_cell = &signed_water,
        .cumulative_outward_runoff_m3_by_cell = &cell_outward_water,
        .cumulative_outward_runoff_m3 = 10,
        .cumulative_outward_runoff_heat_megajoules = 100,
    };
    try apply(.{
        .lon_count = 2,
        .lat_count = 2,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 1,
            .first_row = 0,
            .last_row_inclusive = 1,
        },
        .water_m3_per_step_by_face = .{ &east_water, &west_water, &south_water, &north_water },
        .heat_megajoules_per_step_by_face = .{ &east_heat, &west_heat, &south_heat, &north_heat },
        .negligible_water_m3_by_cell = &.{ 0, 0, 0, 0 },
    }, &state);

    try std.testing.expectEqualSlices(f64, &.{ 6, 3, -1, -4 }, &signed_water);
    try std.testing.expectEqualSlices(
        f64,
        &.{ -6, -3, 1, 4 },
        &cell_outward_water,
    );
    try std.testing.expectEqual(@as(f64, 6), state.cumulative_outward_runoff_m3);
    try std.testing.expectEqual(
        @as(f64, 60),
        state.cumulative_outward_runoff_heat_megajoules,
    );
    var signed_change_m3: f64 = 0;
    var outward_change_m3: f64 = 0;
    for (signed_water, cell_outward_water) |change, outward| {
        signed_change_m3 += change;
        outward_change_m3 += outward;
    }
    try std.testing.expectEqual(@as(f64, 0), signed_change_m3 + outward_change_m3);
    try std.testing.expectEqual(
        outward_change_m3,
        state.cumulative_outward_runoff_m3 - 10,
    );
}

test "water threshold retains signed hourly flux but suppresses cumulative ledgers" {
    const water = [_]f64{0.5};
    const heat = [_]f64{9};
    var signed_water = [_]f64{0};
    var cell_outward_water = [_]f64{7};
    var state: State = .{
        .signed_runoff_water_change_m3_by_cell = &signed_water,
        .cumulative_outward_runoff_m3_by_cell = &cell_outward_water,
        .cumulative_outward_runoff_m3 = 8,
        .cumulative_outward_runoff_heat_megajoules = 10,
    };
    try apply(.{
        .lon_count = 1,
        .lat_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .water_m3_per_step_by_face = .{ &water, &.{0}, &.{0}, &.{0} },
        .heat_megajoules_per_step_by_face = .{ &heat, &.{0}, &.{0}, &.{0} },
        .negligible_water_m3_by_cell = &.{0.5},
    }, &state);
    try std.testing.expectEqual(@as(f64, -0.5), signed_water[0]);
    try std.testing.expectEqual(@as(f64, 7), cell_outward_water[0]);
    try std.testing.expectEqual(@as(f64, 8), state.cumulative_outward_runoff_m3);
    try std.testing.expectEqual(
        @as(f64, 10),
        state.cumulative_outward_runoff_heat_megajoules,
    );
}

test "signed source face flux retains independent REDIST direction multiplier" {
    const west_water = [_]f64{-2};
    const west_heat = [_]f64{-20};
    var signed_water = [_]f64{0};
    var cell_outward_water = [_]f64{0};
    var state: State = .{
        .signed_runoff_water_change_m3_by_cell = &signed_water,
        .cumulative_outward_runoff_m3_by_cell = &cell_outward_water,
        .cumulative_outward_runoff_m3 = 0,
        .cumulative_outward_runoff_heat_megajoules = 0,
    };
    try apply(.{
        .lon_count = 1,
        .lat_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .water_m3_per_step_by_face = .{ &.{0}, &west_water, &.{0}, &.{0} },
        .heat_megajoules_per_step_by_face = .{ &.{0}, &west_heat, &.{0}, &.{0} },
        .negligible_water_m3_by_cell = &.{0},
    }, &state);
    try std.testing.expectEqual(@as(f64, -2), signed_water[0]);
    try std.testing.expectEqual(@as(f64, 2), cell_outward_water[0]);
    try std.testing.expectEqual(@as(f64, 2), state.cumulative_outward_runoff_m3);
    try std.testing.expectEqual(
        @as(f64, 20),
        state.cumulative_outward_runoff_heat_megajoules,
    );
}

test "runtime boundary window does not treat outer storage as its edge" {
    const cell_count = 12;
    const east_water = [_]f64{1} ** cell_count;
    const west_water = [_]f64{2} ** cell_count;
    const south_water = [_]f64{3} ** cell_count;
    const north_water = [_]f64{4} ** cell_count;
    const zero_heat = [_]f64{0} ** cell_count;
    const negligible = [_]f64{0} ** cell_count;
    var signed_water = [_]f64{0} ** cell_count;
    var cell_outward_water = [_]f64{0} ** cell_count;
    var state: State = .{
        .signed_runoff_water_change_m3_by_cell = &signed_water,
        .cumulative_outward_runoff_m3_by_cell = &cell_outward_water,
        .cumulative_outward_runoff_m3 = 0,
        .cumulative_outward_runoff_heat_megajoules = 0,
    };
    try apply(.{
        .lon_count = 4,
        .lat_count = 3,
        .external_boundary_window = .{
            .first_column = 1,
            .last_column_inclusive = 2,
            .first_row = 1,
            .last_row_inclusive = 2,
        },
        .water_m3_per_step_by_face = .{ &east_water, &west_water, &south_water, &north_water },
        .heat_megajoules_per_step_by_face = .{ &zero_heat, &zero_heat, &zero_heat, &zero_heat },
        .negligible_water_m3_by_cell = &negligible,
    }, &state);
    try std.testing.expectEqual(@as(f64, 0), signed_water[0]);
    try std.testing.expectEqual(@as(f64, 6), signed_water[1 * 4 + 1]);
    try std.testing.expectEqual(@as(f64, 3), signed_water[1 * 4 + 2]);
    try std.testing.expectEqual(@as(f64, -1), signed_water[2 * 4 + 1]);
    try std.testing.expectEqual(@as(f64, -4), signed_water[2 * 4 + 2]);
}

test "late invalid boundary flux leaves all accounting unchanged" {
    const valid = [_]f64{ 1, 2 };
    const invalid = [_]f64{ 3, std.math.nan(f64) };
    var signed_water = [_]f64{ 4, 5 };
    var cell_outward_water = [_]f64{ 6, 7 };
    var state: State = .{
        .signed_runoff_water_change_m3_by_cell = &signed_water,
        .cumulative_outward_runoff_m3_by_cell = &cell_outward_water,
        .cumulative_outward_runoff_m3 = 8,
        .cumulative_outward_runoff_heat_megajoules = 9,
    };
    try std.testing.expectError(error.InvalidRunoffBoundaryInput, apply(.{
        .lon_count = 2,
        .lat_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 1,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .water_m3_per_step_by_face = .{ &valid, &valid, &valid, &invalid },
        .heat_megajoules_per_step_by_face = .{ &valid, &valid, &valid, &valid },
        .negligible_water_m3_by_cell = &.{ 0, 0 },
    }, &state));
    try std.testing.expectEqualSlices(f64, &.{ 4, 5 }, &signed_water);
    try std.testing.expectEqualSlices(f64, &.{ 6, 7 }, &cell_outward_water);
    try std.testing.expectEqual(@as(f64, 8), state.cumulative_outward_runoff_m3);
    try std.testing.expectEqual(
        @as(f64, 9),
        state.cumulative_outward_runoff_heat_megajoules,
    );
}

test "invalid boundary dimensions and window fail explicitly" {
    var signed_water = [_]f64{0};
    var cell_outward_water = [_]f64{0};
    var state: State = .{
        .signed_runoff_water_change_m3_by_cell = &signed_water,
        .cumulative_outward_runoff_m3_by_cell = &cell_outward_water,
        .cumulative_outward_runoff_m3 = 0,
        .cumulative_outward_runoff_heat_megajoules = 0,
    };
    try std.testing.expectError(error.InvalidRunoffBoundaryWindow, apply(.{
        .lon_count = 1,
        .lat_count = 1,
        .external_boundary_window = .{
            .first_column = 1,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .water_m3_per_step_by_face = .{ &.{0}, &.{0}, &.{0}, &.{0} },
        .heat_megajoules_per_step_by_face = .{ &.{0}, &.{0}, &.{0}, &.{0} },
        .negligible_water_m3_by_cell = &.{0},
    }, &state));
}
