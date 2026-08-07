//! Tests for `litter_chemistry.zig`.
//!
//! Extracted verbatim so the module beside it contains only the model
//! code. Tests that use private declarations of that module stay there,
//! since a sibling file can only reach `pub` declarations.

const activity_coefficients = @import("../soil/solute/activity_coefficients.zig");
const ledger = @import("../surface/litter_reaction_ledger.zig");
const numerics = @import("../core/numerics.zig");
const std = @import("std");
const surface_litter_chemistry = @import("../surface/litter_chemistry.zig");
test "solid mineral inventory survives wet dry and rewet carrier changes" {
    var state = try surface_litter_chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].salt_minerals.calcite_mol_per_m3 = 2;
    state.cells[0].phosphate_minerals.hydroxyapatite_mol_per_m3 = 3;
    try state.bindMineralReferenceWater(&.{1});

    try state.renormalizeMinerals(0, 0);
    try std.testing.expectEqual(@as(f64, 2), state.cells[0].salt_minerals.calcite_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), state.mineral_reference_water_m3[0]);

    try state.renormalizeMinerals(0, 0.25);
    try std.testing.expectEqual(@as(f64, 8), state.cells[0].salt_minerals.calcite_mol_per_m3);
    try std.testing.expectEqual(
        @as(f64, 12),
        state.cells[0].phosphate_minerals.hydroxyapatite_mol_per_m3,
    );
    try std.testing.expectEqual(
        @as(f64, 2),
        state.cells[0].salt_minerals.calcite_mol_per_m3 * 0.25,
    );
}

fn simultaneousPhosphateSinkEvaluator(
    _: *const anyopaque,
    _: surface_litter_chemistry.Cell,
) !ledger.ReactionExtents {
    var extents = std.mem.zeroes(ledger.ReactionExtents);
    extents.phosphate_minerals.aluminum_phosphate_mol_per_m3 = 0.5;
    extents.phosphate_minerals.iron_phosphate_mol_per_m3 = 0.5;
    extents.phosphate_minerals.dicalcium_phosphate_mol_per_m3 = 0.5;
    return extents;
}

test "fixed-pH simultaneous sinks share one exact admissible substrate fraction" {
    var state = try surface_litter_chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].h2po4_mol_p_per_m3 = 1;
    const context: u8 = 0;
    const result = try surface_litter_chemistry.solveCell(
        &state,
        0,
        .{
            .litter_mass_per_water_volume_megagrams_per_m3 = 1,
            .dynamic_salts = false,
        },
        .{
            .context = &context,
            .evaluate = simultaneousPhosphateSinkEvaluator,
        },
        .{},
    );
    try std.testing.expectEqual(@as(u16, 1), result.iterations);
    try std.testing.expectEqual(@as(f64, 0), state.cells[0].h2po4_mol_p_per_m3);
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.0 / 3.0),
        state.cells[0].phosphate_minerals.aluminum_phosphate_mol_per_m3,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.0 / 3.0),
        state.cells[0].phosphate_minerals.iron_phosphate_mol_per_m3,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.0 / 3.0),
        state.cells[0].phosphate_minerals.dicalcium_phosphate_mol_per_m3,
        1e-15,
    );
}

/// Reproduces the Ottawa surface-litter phosphate tail: a capped, locally
/// constant precipitation rate whose directional derivative is zero, so the
/// solver's probe finds no usable direction and only a phosphate-aware
/// escape can descend. Supersaturation is reported against the aqueous
/// H2PO4 pool, exactly like the production evaluator.
const CappedPhosphatePlateauContext = struct {
    /// mol P m-3 per iteration; the recorded failure had rates in this range
    /// against inventories of order 10 mol m-3.
    capped_rate: f64,
    target_h2po4_mol_p_per_m3: f64,
};

fn cappedPhosphatePlateauEvaluator(
    raw: *const anyopaque,
    cell: surface_litter_chemistry.Cell,
) !ledger.ReactionExtents {
    const context: *const CappedPhosphatePlateauContext =
        @ptrCast(@alignCast(raw));
    var extents = std.mem.zeroes(ledger.ReactionExtents);
    // Clipped to a constant whenever supersaturated: no local gradient.
    if (cell.h2po4_mol_p_per_m3 > context.target_h2po4_mol_p_per_m3)
        extents.phosphate_minerals.aluminum_phosphate_mol_per_m3 =
            context.capped_rate;
    return extents;
}

fn cappedPhosphatePlateauResiduals(
    raw: *const anyopaque,
    cell: surface_litter_chemistry.Cell,
) !ledger.PhosphateMineralExtents {
    const context: *const CappedPhosphatePlateauContext =
        @ptrCast(@alignCast(raw));
    const supersaturation =
        cell.h2po4_mol_p_per_m3 - context.target_h2po4_mol_p_per_m3;
    return .{
        .aluminum_phosphate_mol_per_m3 =
            if (cell.phosphate_minerals.aluminum_phosphate_mol_per_m3 <= 0)
                @max(0, supersaturation)
            else
                supersaturation,
        .iron_phosphate_mol_per_m3 = 0,
        .dicalcium_phosphate_mol_per_m3 = 0,
        .hydroxyapatite_mol_per_m3 = 0,
        .monocalcium_phosphate_mol_per_m3 = 0,
    };
}

test "capped phosphate plateau escapes are reachable in the iterative branch" {
    // Regression for SOLUTE-PHOSPHATE-EQUIL. Every phosphate escape used to be
    // guarded by `!environment.dynamic_salts`, which the fixed-pH early return
    // already excludes from the loop, so the iterative branch had no phosphate
    // Newton direction and degenerated to pure Picard on a zero-gradient
    // plateau. That produced `newton_steps=0` with the ceiling exhausted.
    var state = try surface_litter_chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].h2po4_mol_p_per_m3 = 12;
    state.cells[0].hpo4_mol_p_per_m3 = 2;
    state.cells[0].phosphate_minerals.aluminum_phosphate_mol_per_m3 = 14;
    const phosphorus_before =
        state.cells[0].h2po4_mol_p_per_m3 +
        state.cells[0].hpo4_mol_p_per_m3 +
        state.cells[0].phosphate_minerals.aluminum_phosphate_mol_per_m3;
    const context = CappedPhosphatePlateauContext{
        .capped_rate = 2.5e-3,
        .target_h2po4_mol_p_per_m3 = 4,
    };
    const result = try surface_litter_chemistry.solveCell(
        &state,
        0,
        .{
            .litter_mass_per_water_volume_megagrams_per_m3 = 1,
            .dynamic_salts = true,
        },
        .{
            .context = &context,
            .evaluate = cappedPhosphatePlateauEvaluator,
            .phosphate_mineral_equilibrium_residuals =
                cappedPhosphatePlateauResiduals,
        },
        // Deliberately far below the production ceiling: the fix must converge
        // by finding a direction, never by spending more iterations.
        .{ .max_iterations = 24 },
    );
    try std.testing.expect(result.maximum_scaled_residual <= 1);
    try std.testing.expect(result.newton_raphson_steps > 0);
    // Phosphorus is conserved across the aqueous/mineral split.
    try std.testing.expectApproxEqAbs(
        phosphorus_before,
        state.cells[0].h2po4_mol_p_per_m3 +
            state.cells[0].hpo4_mol_p_per_m3 +
            state.cells[0].phosphate_minerals.aluminum_phosphate_mol_per_m3,
        1e-9,
    );
    // Every inventory stays in its physical domain.
    try std.testing.expect(state.cells[0].h2po4_mol_p_per_m3 >= 0);
    try std.testing.expect(state.cells[0].hpo4_mol_p_per_m3 >= 0);
    try std.testing.expect(
        state.cells[0].phosphate_minerals.aluminum_phosphate_mol_per_m3 >= 0,
    );
}

test "fixed-pH litter never iterates undivided hourly rates" {
    // SOLUTE.F 3996--5250 has no `DO M=` statement and uses the undivided
    // hourly TPD/TSL constants, so the fixed-pH litter branch must commit one
    // increment even when offered a large ceiling. Guards the MRXN
    // multiplication that iterating this branch would introduce.
    var state = try surface_litter_chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cells[0].h2po4_mol_p_per_m3 = 12;
    state.cells[0].phosphate_minerals.aluminum_phosphate_mol_per_m3 = 14;
    const context = CappedPhosphatePlateauContext{
        .capped_rate = 2.5e-3,
        .target_h2po4_mol_p_per_m3 = 4,
    };
    const result = try surface_litter_chemistry.solveCell(
        &state,
        0,
        .{
            .litter_mass_per_water_volume_megagrams_per_m3 = 1,
            .dynamic_salts = false,
        },
        .{
            .context = &context,
            .evaluate = cappedPhosphatePlateauEvaluator,
            .phosphate_mineral_equilibrium_residuals =
                cappedPhosphatePlateauResiduals,
        },
        .{ .max_iterations = 1000 },
    );
    try std.testing.expectEqual(@as(u16, 1), result.iterations);
    try std.testing.expectEqual(@as(u16, 0), result.newton_raphson_steps);
    try std.testing.expectEqual(@as(u16, 0), result.picard_steps);
    // Exactly one hourly capped increment moved from aqueous to solid.
    try std.testing.expectApproxEqAbs(
        @as(f64, 12 - 2.5e-3),
        state.cells[0].h2po4_mol_p_per_m3,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 14 + 2.5e-3),
        state.cells[0].phosphate_minerals.aluminum_phosphate_mol_per_m3,
        1e-15,
    );
}
