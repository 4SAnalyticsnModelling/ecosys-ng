const std = @import("std");
const solute = @import("solute_transport.zig");

pub const Boundary = struct {
    cell_index: usize,
    /// Positive leaves the modeled domain; negative recharges from outside.
    outward_water_flux_m3_per_step: f64,
};

/// Calculates the TRNSFRS external micropore/macropore convection rule in a
/// direction-independent form. Positive output is a gain by the model.
pub fn calculateNetFluxMol(state: *const solute.State, boundary: Boundary, maximum_transport_fraction: f64, discharge_mobility_fraction: []const f64, recharge_concentration_mol_per_m3: []const f64, output_net_flux_mol: []f64) !void {
    if (boundary.cell_index >= state.cell_count) return error.SoluteBoundaryCellIndexOutOfBounds;
    if (discharge_mobility_fraction.len != state.species_count or recharge_concentration_mol_per_m3.len != state.species_count or output_net_flux_mol.len != state.species_count) return error.TransportSpeciesCountMismatch;
    if (!std.math.isFinite(boundary.outward_water_flux_m3_per_step) or !std.math.isFinite(maximum_transport_fraction) or maximum_transport_fraction < 0 or maximum_transport_fraction > 1) return error.InvalidSoluteBoundaryInput;
    const water = state.water_volume_m3[boundary.cell_index];
    if (!std.math.isFinite(water) or water < 0) return error.InvalidSoluteBoundaryInput;
    const amounts = try state.cellAmountsConst(boundary.cell_index);
    if (boundary.outward_water_flux_m3_per_step > 0) {
        const fraction = if (water > 0) @min(maximum_transport_fraction, boundary.outward_water_flux_m3_per_step / water) else 0;
        for (amounts, discharge_mobility_fraction, recharge_concentration_mol_per_m3, output_net_flux_mol) |amount, mobility, concentration, *net| {
            try validateSpeciesInput(amount, mobility, concentration);
            net.* = -amount * mobility * fraction;
        }
    } else if (boundary.outward_water_flux_m3_per_step < 0) {
        for (amounts, discharge_mobility_fraction, recharge_concentration_mol_per_m3, output_net_flux_mol) |amount, mobility, concentration, *net| {
            try validateSpeciesInput(amount, mobility, concentration);
            net.* = -boundary.outward_water_flux_m3_per_step * concentration * mobility;
            if (!std.math.isFinite(net.*)) return error.NonFiniteSoluteBoundaryFlux;
        }
    } else {
        @memset(output_net_flux_mol, 0);
    }
}

/// Sums all external faces before validating and committing, preventing one
/// failed face from leaving a partial boundary update.
pub fn commitNetFluxes(allocator: std.mem.Allocator, state: *solute.State, boundaries: []const Boundary, boundary_net_flux_mol: []const f64) !void {
    if (boundary_net_flux_mol.len != try std.math.mul(usize, boundaries.len, state.species_count)) return error.SoluteBoundaryFluxSizeMismatch;
    const net = try allocator.alloc(f64, state.amount_mol.len);
    defer allocator.free(net);
    @memset(net, 0);
    for (boundaries, 0..) |boundary, boundary_index| {
        if (boundary.cell_index >= state.cell_count) return error.SoluteBoundaryCellIndexOutOfBounds;
        for (0..state.species_count) |species| {
            const value = boundary_net_flux_mol[boundary_index * state.species_count + species];
            if (!std.math.isFinite(value)) return error.NonFiniteSoluteBoundaryFlux;
            net[boundary.cell_index * state.species_count + species] += value;
        }
    }
    for (state.amount_mol, net) |amount, change| if (!std.math.isFinite(amount + change) or amount + change < -1e-12) return error.NegativeSoluteBoundaryCandidate;
    for (state.amount_mol, net) |*amount, change| amount.* = @max(0, amount.* + change);
}

pub fn commitCellNetFlux(state: *solute.State, cell_index: usize, net_flux_mol: []const f64) !void {
    if (cell_index >= state.cell_count) return error.SoluteBoundaryCellIndexOutOfBounds;
    if (net_flux_mol.len != state.species_count) return error.TransportSpeciesCountMismatch;
    const amounts = try state.cellAmounts(cell_index);
    for (amounts, net_flux_mol) |amount, change| if (!std.math.isFinite(change) or !std.math.isFinite(amount + change) or amount + change < -1e-12) return error.NegativeSoluteBoundaryCandidate;
    for (amounts, net_flux_mol) |*amount, change| amount.* = @max(0, amount.* + change);
}

fn validateSpeciesInput(amount: f64, mobility: f64, concentration: f64) !void {
    if (!std.math.isFinite(amount) or amount < 0 or !std.math.isFinite(mobility) or mobility < 0 or mobility > 1 or !std.math.isFinite(concentration) or concentration < 0) return error.InvalidSoluteBoundaryInput;
}

test "external discharge uses donor inventory and VFLWX ceiling" {
    var state = try solute.State.init(std.testing.allocator, 1, 3);
    defer state.deinit();
    state.water_volume_m3[0] = 2;
    @memcpy(try state.cellAmounts(0), &[_]f64{ 4, 6, 8 });
    var net: [3]f64 = undefined;
    try calculateNetFluxMol(&state, .{ .cell_index = 0, .outward_water_flux_m3_per_step = 10 }, 0.25, &[_]f64{ 1, 0.5, 0 }, &[_]f64{ 0, 0, 0 }, &net);
    try std.testing.expectEqualSlices(f64, &[_]f64{ -1, -0.75, 0 }, &net);
}

test "external recharge uses prescribed concentration and mobility" {
    var state = try solute.State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    var net: [2]f64 = undefined;
    try calculateNetFluxMol(&state, .{ .cell_index = 0, .outward_water_flux_m3_per_step = -2 }, 0.5, &[_]f64{ 0.75, 0.25 }, &[_]f64{ 3, 4 }, &net);
    try std.testing.expectEqualSlices(f64, &[_]f64{ 4.5, 2 }, &net);
}

test "multiple external faces commit atomically" {
    var state = try solute.State.init(std.testing.allocator, 1, 1);
    defer state.deinit();
    state.amount_mol[0] = 2;
    const boundaries = [_]Boundary{ .{ .cell_index = 0, .outward_water_flux_m3_per_step = 1 }, .{ .cell_index = 0, .outward_water_flux_m3_per_step = 1 } };
    try commitNetFluxes(std.testing.allocator, &state, &boundaries, &[_]f64{ -0.5, -0.5 });
    try std.testing.expectEqual(@as(f64, 1), state.amount_mol[0]);
    try std.testing.expectError(error.NegativeSoluteBoundaryCandidate, commitNetFluxes(std.testing.allocator, &state, &boundaries, &[_]f64{ -1, -1 }));
    try std.testing.expectEqual(@as(f64, 1), state.amount_mol[0]);
}
