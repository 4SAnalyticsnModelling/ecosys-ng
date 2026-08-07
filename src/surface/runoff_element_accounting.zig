const std = @import("std");

pub const dissolved_organic_fraction_count: usize = 5;

pub const ExternalBoundaryWindow = struct {
    first_column: usize,
    last_column_inclusive: usize,
    first_row: usize,
    last_row_inclusive: usize,
};

/// REDIST evaluates these face records in this exact order.
pub const Face = enum(u8) {
    east,
    west,
    south,
    north,

    pub const count: usize = @typeInfo(Face).@"enum".fields.len;
};

pub const InorganicFlux = struct {
    carbon_dioxide_g_c_per_step: f64 = 0,
    methane_g_c_per_step: f64 = 0,
    ammonium_g_n_per_step: f64 = 0,
    ammonia_g_n_per_step: f64 = 0,
    nitrate_g_n_per_step: f64 = 0,
    nitrite_g_n_per_step: f64 = 0,
    nitrous_oxide_g_n_per_step: f64 = 0,
    dinitrogen_g_n_per_step: f64 = 0,
    h2po4_g_p_per_step: f64 = 0,
    hpo4_g_p_per_step: f64 = 0,
    oxygen_g_o_per_step: f64 = 0,
    hydrogen_g_h_per_step: f64 = 0,
};

pub const OrganicFlux = struct {
    dissolved_organic_carbon_g_c_per_step: f64 = 0,
    dissolved_acetate_carbon_g_c_per_step: f64 = 0,
    dissolved_organic_nitrogen_g_n_per_step: f64 = 0,
    dissolved_organic_phosphorus_g_p_per_step: f64 = 0,
};

pub const BoundaryFlux = struct {
    water_m3_per_step: f64 = 0,
    inorganic: InorganicFlux = .{},
    organic_by_fraction: [dissolved_organic_fraction_count]OrganicFlux =
        [_]OrganicFlux{.{}} ** dissolved_organic_fraction_count,
};

pub const CellExport = struct {
    dissolved_organic_carbon_g_c: f64 = 0,
    dissolved_inorganic_carbon_g_c: f64 = 0,
    dissolved_organic_nitrogen_g_n: f64 = 0,
    dissolved_inorganic_nitrogen_g_n: f64 = 0,
    dissolved_organic_phosphorus_g_p: f64 = 0,
    dissolved_inorganic_phosphorus_g_p: f64 = 0,
};

pub const CumulativeBoundaryElements = struct {
    carbon_g_c: f64 = 0,
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
    oxygen_g_o: f64 = 0,
    /// Retains REDIST's `H2GOU = H2GOU + HGR` sign convention.
    hydrogen_g_h: f64 = 0,
};

pub const Inputs = struct {
    lon_count: usize,
    lat_count: usize,
    external_boundary_window: ExternalBoundaryWindow,
    flux_by_face: [Face.count][]const BoundaryFlux,
    negligible_water_m3_by_cell: []const f64,
};

pub const State = struct {
    cumulative: CumulativeBoundaryElements,
    export_by_cell: []CellExport,
};

/// Accounts for runoff C, N, P, O, and H at external horizontal faces.
///
/// Traceability: REDIST.F lines 726--751 (`CXR`, `ZXR`, `ZGR`, `PXR`,
/// `COR`, `ZOR`, `POR`, `TCOU`, `TZOU`, `TPOU`, `UDOCQ`--`UDIPQ`,
/// `OXYGOU`, and `H2GOU`). This block remains inside the source's strict
/// runoff-water threshold branch. Traversal retains column, row, face, then
/// the five dissolved-organic fractions. Flux records are signed before the
/// independent source direction multiplier is applied.
pub fn apply(inputs: Inputs, state: *State) !void {
    const cell_count = try validateDimensions(inputs, state.*);
    try validateInputsAndState(inputs, state.*, cell_count);
    try preflightUpdates(inputs, state.*);

    const window = inputs.external_boundary_window;
    for (window.first_column..window.last_column_inclusive + 1) |column| {
        for (window.first_row..window.last_row_inclusive + 1) |row| {
            const cell = row * inputs.lon_count + column;
            for (std.meta.tags(Face)) |face| {
                if (!isExternalFace(face, column, row, window)) continue;
                const flux = inputs.flux_by_face[@intFromEnum(face)][cell];
                const direction = inwardSign(face);
                if (@abs(direction * flux.water_m3_per_step) <=
                    inputs.negligible_water_m3_by_cell[cell]) continue;
                const groups = groupsForFlux(flux, direction) catch unreachable;
                commitGroups(&state.cumulative, &state.export_by_cell[cell], groups);
            }
        }
    }
}

const SignedGroups = struct {
    inorganic_carbon_g_c: f64,
    mineral_nitrogen_g_n: f64,
    gaseous_nitrogen_g_n: f64,
    inorganic_phosphorus_g_p: f64,
    organic_carbon_g_c: f64,
    organic_nitrogen_g_n: f64,
    organic_phosphorus_g_p: f64,
    oxygen_g_o: f64,
    hydrogen_g_h: f64,
};

fn groupsForFlux(flux: BoundaryFlux, direction: f64) !SignedGroups {
    const inorganic = flux.inorganic;
    const inorganic_carbon_g_c = try signedSum(
        direction,
        &.{ inorganic.carbon_dioxide_g_c_per_step, inorganic.methane_g_c_per_step },
    );
    const mineral_nitrogen_g_n = try signedSum(direction, &.{
        inorganic.ammonium_g_n_per_step,
        inorganic.ammonia_g_n_per_step,
        inorganic.nitrate_g_n_per_step,
        inorganic.nitrite_g_n_per_step,
    });
    const gaseous_nitrogen_g_n = try signedSum(
        direction,
        &.{ inorganic.nitrous_oxide_g_n_per_step, inorganic.dinitrogen_g_n_per_step },
    );
    const inorganic_phosphorus_g_p = try signedSum(
        direction,
        &.{ inorganic.h2po4_g_p_per_step, inorganic.hpo4_g_p_per_step },
    );
    var organic_carbon_g_c: f64 = 0;
    var organic_nitrogen_g_n: f64 = 0;
    var organic_phosphorus_g_p: f64 = 0;
    for (flux.organic_by_fraction) |organic| {
        organic_carbon_g_c = try checkedSum(
            organic_carbon_g_c,
            try signedSum(direction, &.{
                organic.dissolved_organic_carbon_g_c_per_step,
                organic.dissolved_acetate_carbon_g_c_per_step,
            }),
        );
        organic_nitrogen_g_n = try checkedSum(
            organic_nitrogen_g_n,
            try checkedProduct(
                direction,
                organic.dissolved_organic_nitrogen_g_n_per_step,
            ),
        );
        organic_phosphorus_g_p = try checkedSum(
            organic_phosphorus_g_p,
            try checkedProduct(
                direction,
                organic.dissolved_organic_phosphorus_g_p_per_step,
            ),
        );
    }
    return .{
        .inorganic_carbon_g_c = inorganic_carbon_g_c,
        .mineral_nitrogen_g_n = mineral_nitrogen_g_n,
        .gaseous_nitrogen_g_n = gaseous_nitrogen_g_n,
        .inorganic_phosphorus_g_p = inorganic_phosphorus_g_p,
        .organic_carbon_g_c = organic_carbon_g_c,
        .organic_nitrogen_g_n = organic_nitrogen_g_n,
        .organic_phosphorus_g_p = organic_phosphorus_g_p,
        .oxygen_g_o = try checkedProduct(direction, inorganic.oxygen_g_o_per_step),
        .hydrogen_g_h = try checkedProduct(direction, inorganic.hydrogen_g_h_per_step),
    };
}

fn preflightUpdates(inputs: Inputs, state: State) !void {
    var cumulative = state.cumulative;
    const window = inputs.external_boundary_window;
    for (window.first_column..window.last_column_inclusive + 1) |column| {
        for (window.first_row..window.last_row_inclusive + 1) |row| {
            const cell = row * inputs.lon_count + column;
            var cell_export = state.export_by_cell[cell];
            for (std.meta.tags(Face)) |face| {
                if (!isExternalFace(face, column, row, window)) continue;
                const flux = inputs.flux_by_face[@intFromEnum(face)][cell];
                const direction = inwardSign(face);
                if (@abs(direction * flux.water_m3_per_step) <=
                    inputs.negligible_water_m3_by_cell[cell]) continue;
                const groups = try groupsForFlux(flux, direction);
                try accumulateGroups(&cumulative, &cell_export, groups);
            }
        }
    }
}

fn accumulateGroups(
    cumulative: *CumulativeBoundaryElements,
    cell_export: *CellExport,
    groups: SignedGroups,
) !void {
    cumulative.carbon_g_c =
        try checkedSum(cumulative.carbon_g_c, -groups.inorganic_carbon_g_c);
    cumulative.carbon_g_c =
        try checkedSum(cumulative.carbon_g_c, -groups.organic_carbon_g_c);
    cumulative.nitrogen_g_n =
        try checkedSum(cumulative.nitrogen_g_n, -groups.mineral_nitrogen_g_n);
    cumulative.nitrogen_g_n =
        try checkedSum(cumulative.nitrogen_g_n, -groups.organic_nitrogen_g_n);
    cumulative.nitrogen_g_n =
        try checkedSum(cumulative.nitrogen_g_n, -groups.gaseous_nitrogen_g_n);
    cumulative.phosphorus_g_p = try checkedSum(
        cumulative.phosphorus_g_p,
        -groups.inorganic_phosphorus_g_p,
    );
    cumulative.phosphorus_g_p =
        try checkedSum(cumulative.phosphorus_g_p, -groups.organic_phosphorus_g_p);
    cell_export.dissolved_organic_carbon_g_c = try checkedSum(
        cell_export.dissolved_organic_carbon_g_c,
        -groups.organic_carbon_g_c,
    );
    cell_export.dissolved_inorganic_carbon_g_c = try checkedSum(
        cell_export.dissolved_inorganic_carbon_g_c,
        -groups.inorganic_carbon_g_c,
    );
    cell_export.dissolved_organic_nitrogen_g_n = try checkedSum(
        cell_export.dissolved_organic_nitrogen_g_n,
        -groups.organic_nitrogen_g_n,
    );
    cell_export.dissolved_inorganic_nitrogen_g_n = try checkedSum(
        cell_export.dissolved_inorganic_nitrogen_g_n,
        -groups.mineral_nitrogen_g_n,
    );
    cell_export.dissolved_organic_phosphorus_g_p = try checkedSum(
        cell_export.dissolved_organic_phosphorus_g_p,
        -groups.organic_phosphorus_g_p,
    );
    cell_export.dissolved_inorganic_phosphorus_g_p = try checkedSum(
        cell_export.dissolved_inorganic_phosphorus_g_p,
        -groups.inorganic_phosphorus_g_p,
    );
    cumulative.oxygen_g_o =
        try checkedSum(cumulative.oxygen_g_o, -groups.oxygen_g_o);
    cumulative.hydrogen_g_h =
        try checkedSum(cumulative.hydrogen_g_h, groups.hydrogen_g_h);
}

fn commitGroups(
    cumulative: *CumulativeBoundaryElements,
    cell_export: *CellExport,
    groups: SignedGroups,
) void {
    cumulative.carbon_g_c -= groups.inorganic_carbon_g_c;
    cumulative.carbon_g_c -= groups.organic_carbon_g_c;
    cumulative.nitrogen_g_n -= groups.mineral_nitrogen_g_n;
    cumulative.nitrogen_g_n -= groups.organic_nitrogen_g_n;
    cumulative.nitrogen_g_n -= groups.gaseous_nitrogen_g_n;
    cumulative.phosphorus_g_p -= groups.inorganic_phosphorus_g_p;
    cumulative.phosphorus_g_p -= groups.organic_phosphorus_g_p;
    cell_export.dissolved_organic_carbon_g_c -= groups.organic_carbon_g_c;
    cell_export.dissolved_inorganic_carbon_g_c -= groups.inorganic_carbon_g_c;
    cell_export.dissolved_organic_nitrogen_g_n -= groups.organic_nitrogen_g_n;
    cell_export.dissolved_inorganic_nitrogen_g_n -= groups.mineral_nitrogen_g_n;
    cell_export.dissolved_organic_phosphorus_g_p -= groups.organic_phosphorus_g_p;
    cell_export.dissolved_inorganic_phosphorus_g_p -=
        groups.inorganic_phosphorus_g_p;
    cumulative.oxygen_g_o -= groups.oxygen_g_o;
    cumulative.hydrogen_g_h += groups.hydrogen_g_h;
}

fn validateDimensions(inputs: Inputs, state: State) !usize {
    if (inputs.lon_count == 0 or inputs.lat_count == 0)
        return error.InvalidRunoffElementDimensions;
    const cell_count = std.math.mul(
        usize,
        inputs.lon_count,
        inputs.lat_count,
    ) catch return error.InvalidRunoffElementDimensions;
    inline for (inputs.flux_by_face) |values|
        if (values.len != cell_count)
            return error.InvalidRunoffElementDimensions;
    if (inputs.negligible_water_m3_by_cell.len != cell_count or
        state.export_by_cell.len != cell_count)
        return error.InvalidRunoffElementDimensions;
    const window = inputs.external_boundary_window;
    if (window.first_column > window.last_column_inclusive or
        window.first_row > window.last_row_inclusive or
        window.last_column_inclusive >= inputs.lon_count or
        window.last_row_inclusive >= inputs.lat_count)
        return error.InvalidRunoffElementBoundaryWindow;
    return cell_count;
}

fn validateInputsAndState(inputs: Inputs, state: State, cell_count: usize) !void {
    inline for (inputs.flux_by_face) |values| {
        for (values) |flux| {
            if (!std.math.isFinite(flux.water_m3_per_step))
                return error.InvalidRunoffElementInput;
            try validateNumericStruct(flux.inorganic, error.InvalidRunoffElementInput);
            for (flux.organic_by_fraction) |organic|
                try validateNumericStruct(organic, error.InvalidRunoffElementInput);
        }
    }
    for (0..cell_count) |cell| {
        if (!nonnegativeFinite(inputs.negligible_water_m3_by_cell[cell]))
            return error.InvalidRunoffElementInput;
        try validateNumericStruct(
            state.export_by_cell[cell],
            error.InvalidRunoffElementState,
        );
    }
    try validateNumericStruct(state.cumulative, error.InvalidRunoffElementState);
}

fn isExternalFace(
    face: Face,
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

fn inwardSign(face: Face) f64 {
    return switch (face) {
        .east, .south => -1,
        .west, .north => 1,
    };
}

fn signedSum(direction: f64, values: []const f64) !f64 {
    var sum: f64 = 0;
    for (values) |value| sum = try checkedSum(sum, value);
    return checkedProduct(direction, sum);
}

fn checkedProduct(left: f64, right: f64) !f64 {
    const result = left * right;
    if (!std.math.isFinite(result)) return error.NonFiniteRunoffElementResult;
    return result;
}

fn checkedSum(current: f64, increment: f64) !f64 {
    const result = current + increment;
    if (!std.math.isFinite(result)) return error.NonFiniteRunoffElementResult;
    return result;
}

fn validateNumericStruct(value: anytype, comptime failure: anyerror) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name))) return failure;
}

fn nonnegativeFinite(value: f64) bool {
    return std.math.isFinite(value) and value >= 0;
}

fn sourceFixtureFlux() BoundaryFlux {
    var result: BoundaryFlux = .{
        .water_m3_per_step = 1,
        .inorganic = .{
            .carbon_dioxide_g_c_per_step = 1,
            .methane_g_c_per_step = 2,
            .ammonium_g_n_per_step = 1,
            .ammonia_g_n_per_step = 2,
            .nitrate_g_n_per_step = 3,
            .nitrite_g_n_per_step = 4,
            .nitrous_oxide_g_n_per_step = 5,
            .dinitrogen_g_n_per_step = 6,
            .h2po4_g_p_per_step = 8,
            .hpo4_g_p_per_step = 9,
            .oxygen_g_o_per_step = 11,
            .hydrogen_g_h_per_step = 12,
        },
    };
    for (&result.organic_by_fraction) |*organic| organic.* = .{
        .dissolved_organic_carbon_g_c_per_step = 1,
        .dissolved_acetate_carbon_g_c_per_step = 2,
        .dissolved_organic_nitrogen_g_n_per_step = 7,
        .dissolved_organic_phosphorus_g_p_per_step = 10,
    };
    return result;
}

test "REDIST runoff element groups preserve source signs and closure" {
    const east = [_]BoundaryFlux{sourceFixtureFlux()};
    const zero = [_]BoundaryFlux{.{}};
    var exports = [_]CellExport{.{}};
    var state: State = .{ .cumulative = .{}, .export_by_cell = &exports };
    try apply(.{
        .lon_count = 1,
        .lat_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .flux_by_face = .{ &east, &zero, &zero, &zero },
        .negligible_water_m3_by_cell = &.{0},
    }, &state);

    try std.testing.expectEqual(@as(f64, 18), state.cumulative.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 56), state.cumulative.nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 67), state.cumulative.phosphorus_g_p);
    try std.testing.expectEqual(@as(f64, 11), state.cumulative.oxygen_g_o);
    try std.testing.expectEqual(@as(f64, -12), state.cumulative.hydrogen_g_h);
    try std.testing.expectEqual(
        @as(f64, 15),
        exports[0].dissolved_organic_carbon_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 3),
        exports[0].dissolved_inorganic_carbon_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 35),
        exports[0].dissolved_organic_nitrogen_g_n,
    );
    try std.testing.expectEqual(
        @as(f64, 10),
        exports[0].dissolved_inorganic_nitrogen_g_n,
    );
    try std.testing.expectEqual(
        @as(f64, 50),
        exports[0].dissolved_organic_phosphorus_g_p,
    );
    try std.testing.expectEqual(
        @as(f64, 17),
        exports[0].dissolved_inorganic_phosphorus_g_p,
    );
    try std.testing.expectEqual(
        state.cumulative.carbon_g_c,
        exports[0].dissolved_organic_carbon_g_c +
            exports[0].dissolved_inorganic_carbon_g_c,
    );
    try std.testing.expectEqual(
        state.cumulative.phosphorus_g_p,
        exports[0].dissolved_organic_phosphorus_g_p +
            exports[0].dissolved_inorganic_phosphorus_g_p,
    );
    try std.testing.expectEqual(
        state.cumulative.nitrogen_g_n,
        exports[0].dissolved_organic_nitrogen_g_n +
            exports[0].dissolved_inorganic_nitrogen_g_n + 11,
    );
}

test "west face reverses source group signs and adds to existing state" {
    const west = [_]BoundaryFlux{sourceFixtureFlux()};
    const zero = [_]BoundaryFlux{.{}};
    var exports = [_]CellExport{.{ .dissolved_organic_carbon_g_c = 20 }};
    var state: State = .{
        .cumulative = .{ .carbon_g_c = 30 },
        .export_by_cell = &exports,
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
        .flux_by_face = .{ &zero, &west, &zero, &zero },
        .negligible_water_m3_by_cell = &.{0},
    }, &state);
    try std.testing.expectEqual(@as(f64, 12), state.cumulative.carbon_g_c);
    try std.testing.expectEqual(
        @as(f64, 5),
        exports[0].dissolved_organic_carbon_g_c,
    );
    try std.testing.expectEqual(@as(f64, 12), state.cumulative.hydrogen_g_h);
}

test "runtime boundary window and water threshold gate element accounting" {
    const cell_count = 6;
    var east = [_]BoundaryFlux{.{}} ** cell_count;
    east[4] = sourceFixtureFlux();
    east[5] = sourceFixtureFlux();
    east[5].water_m3_per_step = 0.25;
    const zero = [_]BoundaryFlux{.{}} ** cell_count;
    var exports = [_]CellExport{.{}} ** cell_count;
    var state: State = .{ .cumulative = .{}, .export_by_cell = &exports };
    try apply(.{
        .lon_count = 3,
        .lat_count = 2,
        .external_boundary_window = .{
            .first_column = 1,
            .last_column_inclusive = 2,
            .first_row = 1,
            .last_row_inclusive = 1,
        },
        .flux_by_face = .{ &east, &zero, &zero, &zero },
        .negligible_water_m3_by_cell = &.{ 0, 0, 0, 0, 0, 0.25 },
    }, &state);
    try std.testing.expectEqual(CumulativeBoundaryElements{}, state.cumulative);
    for (exports) |value| try std.testing.expectEqual(CellExport{}, value);
}

test "late invalid runoff element flux leaves state unchanged" {
    var east = [_]BoundaryFlux{ .{}, .{} };
    east[0] = sourceFixtureFlux();
    east[1].organic_by_fraction[4].dissolved_organic_phosphorus_g_p_per_step =
        std.math.nan(f64);
    const zero = [_]BoundaryFlux{ .{}, .{} };
    var exports = [_]CellExport{
        .{ .dissolved_organic_carbon_g_c = 1 },
        .{ .dissolved_organic_carbon_g_c = 2 },
    };
    var state: State = .{
        .cumulative = .{ .carbon_g_c = 3 },
        .export_by_cell = &exports,
    };
    try std.testing.expectError(error.InvalidRunoffElementInput, apply(.{
        .lon_count = 2,
        .lat_count = 1,
        .external_boundary_window = .{
            .first_column = 0,
            .last_column_inclusive = 1,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .flux_by_face = .{ &east, &zero, &zero, &zero },
        .negligible_water_m3_by_cell = &.{ 0, 0 },
    }, &state));
    try std.testing.expectEqual(@as(f64, 3), state.cumulative.carbon_g_c);
    try std.testing.expectEqual(
        @as(f64, 1),
        exports[0].dissolved_organic_carbon_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 2),
        exports[1].dissolved_organic_carbon_g_c,
    );
}

test "invalid runoff element dimensions and window fail explicitly" {
    const zero = [_]BoundaryFlux{.{}};
    var exports = [_]CellExport{.{}};
    var state: State = .{ .cumulative = .{}, .export_by_cell = &exports };
    try std.testing.expectError(error.InvalidRunoffElementBoundaryWindow, apply(.{
        .lon_count = 1,
        .lat_count = 1,
        .external_boundary_window = .{
            .first_column = 1,
            .last_column_inclusive = 0,
            .first_row = 0,
            .last_row_inclusive = 0,
        },
        .flux_by_face = .{ &zero, &zero, &zero, &zero },
        .negligible_water_m3_by_cell = &.{0},
    }, &state));
}
