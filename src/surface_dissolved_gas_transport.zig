const std = @import("std");
const gas = @import("gas_transport.zig");
const runoff_carrier = @import("surface_runoff_carrier.zig");

pub const Directions = struct {
    east_m3: []const f64,
    west_m3: []const f64,
    south_m3: []const f64,
    north_m3: []const f64,
};

/// TRNSFR `RQRCOS/RQRCHS/...`: route every dissolved litter gas with the
/// converged surface runoff. Directional transfers share the pre-runoff
/// inventory, so routing is order independent and commits atomically.
pub fn advance(
    allocator: std.mem.Allocator,
    state: *gas.State,
    columns: usize,
    rows: usize,
    post_runoff_water_m3: []const f64,
    runoff_water_change_m3: []const f64,
    directions: Directions,
    maximum_transport_fraction: f64,
    inorganic_carbon_export_g_c_by_cell: []f64,
) !void {
    const cells = try std.math.mul(usize, columns, rows);
    if (state.cell_count != cells or post_runoff_water_m3.len != cells or runoff_water_change_m3.len != cells or inorganic_carbon_export_g_c_by_cell.len != cells) return error.SurfaceDissolvedGasDimensionMismatch;
    inline for (.{ directions.east_m3, directions.west_m3, directions.south_m3, directions.north_m3 }) |values| if (values.len != cells) return error.SurfaceDissolvedGasDimensionMismatch;
    if (!std.math.isFinite(maximum_transport_fraction) or maximum_transport_fraction < 0 or maximum_transport_fraction > 1) return error.InvalidSurfaceDissolvedGasTransportParameter;

    const original = try allocator.dupe(f64, state.dissolved_mass_g);
    defer allocator.free(original);
    const candidate = try allocator.dupe(f64, original);
    defer allocator.free(candidate);
    const boundary_export = try allocator.alloc(f64, original.len);
    defer allocator.free(boundary_export);
    const pre_runoff_water_m3 = try allocator.alloc(f64, cells);
    defer allocator.free(pre_runoff_water_m3);
    @memset(inorganic_carbon_export_g_c_by_cell, 0);

    for (0..cells) |cell| {
        const before = post_runoff_water_m3[cell] - runoff_water_change_m3[cell];
        if (!std.math.isFinite(before) or before < -1e-12 or !std.math.isFinite(post_runoff_water_m3[cell]) or post_runoff_water_m3[cell] < 0) return error.InvalidSurfaceDissolvedGasWaterState;
        pre_runoff_water_m3[cell] = @max(0, before);
        for (original[cell * gas.species_count ..][0..gas.species_count]) |mass| if (!std.math.isFinite(mass) or mass < 0) return error.InvalidSurfaceDissolvedGasPool;
    }

    try runoff_carrier.calculateChanges(
        columns,
        rows,
        gas.species_count,
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
    for (0..cells) |cell| {
        inline for (.{
            gas.Species.carbon_dioxide,
            gas.Species.methane,
        }) |species| {
            inorganic_carbon_export_g_c_by_cell[cell] += boundary_export[
                cell * gas.species_count + @intFromEnum(species)
            ];
        }
    }
    for (candidate) |mass| if (!std.math.isFinite(mass) or mass < -1e-12) return error.InvalidSurfaceDissolvedGasCandidate;
    for (state.dissolved_mass_g, candidate) |*mass, next| mass.* = @max(0, next);
}

test "surface dissolved gases conserve internal routing and report only boundary carbon" {
    var state = try gas.State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.dissolved_mass_g[@intFromEnum(gas.Species.carbon_dioxide)] = 8;
    state.dissolved_mass_g[@intFromEnum(gas.Species.methane)] = 4;
    state.dissolved_mass_g[@intFromEnum(gas.Species.oxygen)] = 2;
    const second = gas.species_count;
    state.dissolved_mass_g[second + @intFromEnum(gas.Species.carbon_dioxide)] = 2;
    var export_g_c = [_]f64{ 0, 0 };
    const zero = [_]f64{ 0, 0 };
    try advance(std.testing.allocator, &state, 2, 1, &.{ 0.5, 1 }, &.{ -0.5, 0 }, .{
        .east_m3 = &.{ 0.5, 1 },
        .west_m3 = &zero,
        .south_m3 = &zero,
        .north_m3 = &zero,
    }, 1, &export_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 0), export_g_c[0], 1e-15);
    // Simultaneous routing does not cascade cell 0's incoming mass through
    // cell 1 during the same accepted transport step.
    try std.testing.expectApproxEqAbs(@as(f64, 2), export_g_c[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), state.dissolved_mass_g[second + @intFromEnum(gas.Species.oxygen)], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 6), state.dissolved_mass_g[second + @intFromEnum(gas.Species.carbon_dioxide)] + state.dissolved_mass_g[second + @intFromEnum(gas.Species.methane)], 1e-15);
}
