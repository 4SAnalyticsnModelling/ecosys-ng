const std = @import("std");
const snow = @import("snow_solute_transport.zig");
const routing = @import("surface_solute_routing.zig");

/// Routes the ten tracked inventories in snow layer 1 using QST/QSTN and
/// VOLSL(1). Incoming drift activates the destination surface snow layer.
pub fn route(allocator: std.mem.Allocator, state: *snow.State, columns: usize, rows: usize, carrier_volume_m3: []const f64, total_transfer_m3: []const f64, directions: routing.Directions, boundaries: routing.BoundaryConditions, maximum_transport_fraction: f64, exported_g: []f64) !void {
    const cells = try std.math.mul(usize, columns, rows);
    if (cells != state.cell_count or carrier_volume_m3.len != cells or total_transfer_m3.len != cells or directions.east_m3.len != cells or directions.west_m3.len != cells or directions.south_m3.len != cells or directions.north_m3.len != cells or boundaries.east_open.len != cells or boundaries.west_open.len != cells or boundaries.south_open.len != cells or boundaries.north_open.len != cells or exported_g.len != snow.species_count) return error.SnowDriftDimensionMismatch;
    if (!std.math.isFinite(maximum_transport_fraction) or maximum_transport_fraction < 0 or maximum_transport_fraction > 1) return error.InvalidSnowDriftInput;
    const candidate = try allocator.dupe(f64, state.amount_g);
    defer allocator.free(candidate);
    const active_candidate = try allocator.dupe(bool, state.active);
    defer allocator.free(active_candidate);
    @memset(exported_g, 0);
    for (0..cells) |cell| {
        const carrier = carrier_volume_m3[cell];
        const total = total_transfer_m3[cell];
        if (!std.math.isFinite(carrier) or carrier < 0 or !std.math.isFinite(total) or total < 0) return error.InvalidSnowDriftInput;
        const directional = [_]f64{ directions.east_m3[cell], directions.west_m3[cell], directions.south_m3[cell], directions.north_m3[cell] };
        var directional_sum: f64 = 0;
        for (directional) |value| {
            if (!std.math.isFinite(value) or value < 0) return error.InvalidSnowDriftInput;
            directional_sum += value;
        }
        if (directional_sum > total + 64 * std.math.floatEps(f64) * @max(1, total)) return error.DirectionalSnowDriftExceedsTotal;
        if (total == 0) continue;
        const transported_fraction = if (carrier > 0) @min(maximum_transport_fraction, total / carrier) else maximum_transport_fraction;
        const column = cell % columns;
        const row = cell / columns;
        const neighbors = [_]?usize{ if (column + 1 < columns) cell + 1 else null, if (column > 0) cell - 1 else null, if (row + 1 < rows) cell + columns else null, if (row > 0) cell - columns else null };
        const open = [_]bool{ boundaries.east_open[cell], boundaries.west_open[cell], boundaries.south_open[cell], boundaries.north_open[cell] };
        const source_layer = try state.layerIndex(cell, 0);
        for (directional, neighbors, open) |direction_volume, neighbor, is_open| {
            if (direction_volume == 0 or (neighbor == null and !is_open)) continue;
            const fraction = transported_fraction * direction_volume / total;
            for (0..snow.species_count) |species| {
                const source = source_layer * snow.species_count + species;
                const flux = state.amount_g[source] * fraction;
                candidate[source] -= flux;
                if (neighbor) |destination_cell| {
                    const destination_layer = try state.layerIndex(destination_cell, 0);
                    candidate[destination_layer * snow.species_count + species] += flux;
                    if (flux > 0) active_candidate[destination_layer] = true;
                } else exported_g[species] += flux;
            }
        }
    }
    for (candidate) |value| if (!std.math.isFinite(value) or value < -1e-12) return error.InvalidSnowDriftCandidate;
    for (candidate, state.amount_g) |value, *amount| amount.* = @max(0, value);
    @memcpy(state.active, active_candidate);
}

test "snow drift conserves ten species and activates receiving snow" {
    var state = try snow.State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    state.active[0] = true;
    @memset(try state.amounts(0, 0), 4);
    const east = [_]f64{ 0.5, 0 };
    const zero = [_]f64{ 0, 0 };
    const no = [_]bool{ false, false };
    var exported: [snow.species_count]f64 = undefined;
    try route(std.testing.allocator, &state, 2, 1, &[_]f64{ 2, 0 }, &[_]f64{ 0.5, 0 }, .{ .east_m3 = &east, .west_m3 = &zero, .south_m3 = &zero, .north_m3 = &zero }, .{ .east_open = &no, .west_open = &no, .south_open = &no, .north_open = &no }, 0.5, &exported);
    try std.testing.expect(state.active[2]);
    for (0..snow.species_count) |species| try std.testing.expectApproxEqAbs(@as(f64, 4), state.amount_g[species] + state.amount_g[2 * snow.species_count + species], 1e-14);
}

test "open snow boundary accounts exported mass" {
    var state = try snow.State.init(std.testing.allocator, 1, 1);
    defer state.deinit();
    state.active[0] = true;
    state.amount_g[0] = 8;
    const one = [_]f64{1};
    const zero = [_]f64{0};
    const yes = [_]bool{true};
    const no = [_]bool{false};
    var exported: [snow.species_count]f64 = undefined;
    try route(std.testing.allocator, &state, 1, 1, &[_]f64{2}, &one, .{ .east_m3 = &one, .west_m3 = &zero, .south_m3 = &zero, .north_m3 = &zero }, .{ .east_open = &yes, .west_open = &no, .south_open = &no, .north_open = &no }, 1, &exported);
    try std.testing.expectEqual(@as(f64, 4), state.amount_g[0]);
    try std.testing.expectEqual(@as(f64, 4), exported[0]);
}
