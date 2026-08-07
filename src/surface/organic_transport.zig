const std = @import("std");
const organic = @import("../soil/organic/initialization.zig");
const runoff_carrier = @import("runoff_carrier.zig");

pub const Directions = struct {
    east_m3: []const f64,
    west_m3: []const f64,
    south_m3: []const f64,
    north_m3: []const f64,
};

pub const Output = struct {
    dissolved_organic_carbon_export_g_c_by_cell: []f64,
    dissolved_organic_nitrogen_export_g_n_by_cell: []f64,
    dissolved_organic_phosphorus_export_g_p_by_cell: []f64,
};

const components_per_substrate: usize = 4;
const component_count: usize = organic.substrate_count * components_per_substrate;

/// TRNSFR surface transport of five DOC/DON/DOP/acetate complexes using the
/// converged hourly runoff. Every directional transfer is based on the same
/// pre-runoff donor inventory and the complete update commits atomically.
pub fn advance(
    allocator: std.mem.Allocator,
    state: *organic.State,
    columns: usize,
    rows: usize,
    post_runoff_water_m3: []const f64,
    runoff_water_change_m3: []const f64,
    directions: Directions,
    maximum_transport_fraction: f64,
    output: Output,
) !void {
    const cells = try std.math.mul(usize, columns, rows);
    try validate(state, cells, post_runoff_water_m3, runoff_water_change_m3, directions, maximum_transport_fraction, output);
    const original = try allocator.alloc(f64, cells * component_count);
    defer allocator.free(original);
    const candidate = try allocator.alloc(f64, original.len);
    defer allocator.free(candidate);
    const boundary_export = try allocator.alloc(f64, original.len);
    defer allocator.free(boundary_export);
    const pre_runoff_water_m3 = try allocator.alloc(f64, cells);
    defer allocator.free(pre_runoff_water_m3);
    @memset(output.dissolved_organic_carbon_export_g_c_by_cell, 0);
    @memset(output.dissolved_organic_nitrogen_export_g_n_by_cell, 0);
    @memset(output.dissolved_organic_phosphorus_export_g_p_by_cell, 0);

    for (0..cells) |cell| {
        const before = post_runoff_water_m3[cell] - runoff_water_change_m3[cell];
        if (!std.math.isFinite(before) or before < -1e-12 or !std.math.isFinite(post_runoff_water_m3[cell]) or post_runoff_water_m3[cell] < 0) return error.InvalidSurfaceOrganicWaterState;
        pre_runoff_water_m3[cell] = @max(0, before);
        for (0..organic.substrate_count) |substrate| {
            const pool = state.dissolved[cell * organic.substrate_count + substrate];
            const acetate = state.dissolved_acetate_carbon_g_c[cell * organic.substrate_count + substrate];
            inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p, acetate }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSurfaceOrganicPool;
            const base = cell * component_count + substrate * components_per_substrate;
            original[base] = pool.carbon_g_c;
            original[base + 1] = pool.nitrogen_g_n;
            original[base + 2] = pool.phosphorus_g_p;
            original[base + 3] = acetate;
        }
    }
    try runoff_carrier.calculateChanges(
        columns,
        rows,
        component_count,
        original,
        pre_runoff_water_m3,
        .{
            .east_m3 = directions.east_m3,
            .west_m3 = directions.west_m3,
            .south_m3 = directions.south_m3,
            .north_m3 = directions.north_m3,
        },
        maximum_transport_fraction,
        candidate,
        boundary_export,
    );
    for (candidate, original) |*change, amount| change.* += amount;
    for (0..cells) |cell| for (0..component_count) |component| {
        const exported =
            boundary_export[cell * component_count + component];
        switch (component % components_per_substrate) {
            0, 3 => output.dissolved_organic_carbon_export_g_c_by_cell[cell] +=
                exported,
            1 => output.dissolved_organic_nitrogen_export_g_n_by_cell[cell] +=
                exported,
            2 => output.dissolved_organic_phosphorus_export_g_p_by_cell[cell] +=
                exported,
            else => unreachable,
        }
    };
    for (candidate) |value| if (!std.math.isFinite(value) or value < -1e-12) return error.InvalidSurfaceOrganicTransportCandidate;
    for (0..cells) |cell| for (0..organic.substrate_count) |substrate| {
        const base = cell * component_count + substrate * components_per_substrate;
        state.dissolved[cell * organic.substrate_count + substrate] = .{ .carbon_g_c = @max(0, candidate[base]), .nitrogen_g_n = @max(0, candidate[base + 1]), .phosphorus_g_p = @max(0, candidate[base + 2]) };
        state.dissolved_acetate_carbon_g_c[cell * organic.substrate_count + substrate] = @max(0, candidate[base + 3]);
    };
}

fn validate(state: *const organic.State, cells: usize, water: []const f64, water_change: []const f64, directions: Directions, maximum_fraction: f64, output: Output) !void {
    if (state.layer_count != cells or water.len != cells or water_change.len != cells or output.dissolved_organic_carbon_export_g_c_by_cell.len != cells or output.dissolved_organic_nitrogen_export_g_n_by_cell.len != cells or output.dissolved_organic_phosphorus_export_g_p_by_cell.len != cells) return error.SurfaceOrganicTransportDimensionMismatch;
    inline for (.{ directions.east_m3, directions.west_m3, directions.south_m3, directions.north_m3 }) |values| if (values.len != cells) return error.SurfaceOrganicTransportDimensionMismatch;
    if (!std.math.isFinite(maximum_fraction) or maximum_fraction < 0 or maximum_fraction > 1) return error.InvalidSurfaceOrganicTransportParameter;
}

test "surface organic runoff conserves internal transfer and reports source-cell export" {
    var state = try organic.State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.dissolved[0] = .{ .carbon_g_c = 8, .nitrogen_g_n = 4, .phosphorus_g_p = 2 };
    state.dissolved_acetate_carbon_g_c[0] = 6;
    state.dissolved[organic.substrate_count] = .{ .carbon_g_c = 5, .nitrogen_g_n = 2, .phosphorus_g_p = 1 };
    state.dissolved_acetate_carbon_g_c[organic.substrate_count] = 2;
    var carbon = [_]f64{ 0, 0 };
    var nitrogen = [_]f64{ 0, 0 };
    var phosphorus = [_]f64{ 0, 0 };
    const zero = [_]f64{ 0, 0 };
    try advance(std.testing.allocator, &state, 2, 1, &.{ 0.5, 1 }, &.{ -0.5, 0 }, .{ .east_m3 = &.{ 0.5, 1 }, .west_m3 = &zero, .south_m3 = &zero, .north_m3 = &zero }, 1, .{
        .dissolved_organic_carbon_export_g_c_by_cell = &carbon,
        .dissolved_organic_nitrogen_export_g_n_by_cell = &nitrogen,
        .dissolved_organic_phosphorus_export_g_p_by_cell = &phosphorus,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0), carbon[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 7), carbon[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2), nitrogen[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), phosphorus[1], 1e-15);
}
