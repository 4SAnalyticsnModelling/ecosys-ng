const std = @import("std");

pub const Pool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const BranchState = struct {
    nonstructural: Pool,
    structural: Pool,
};

pub const BranchFluxes = struct {
    maintenance_respiration_g_c: f64,
    oxygen_unconstrained_respiration_g_c: f64,
    nitrogen_fixation_respiration_g_c: f64,
    total_growth_carbon_use_g_c: f64,
    recycled_decomposition: Pool,
    recycled_senescence: Pool,
    fixed_nitrogen_g_n: f64,
    structural_growth_g_c: f64,
    nitrogen_used_for_growth_g_n: f64,
    phosphorus_used_for_growth_g_p: f64,
    decomposed_structural: Pool,
    senesced_structural: Pool,
};

pub const Inputs = struct {
    fixation_type: u8,
};

/// Exact GROSUB lines 5609--5627 canopy symbiont nonstructural then
/// structural C:N:P update in ascending NB order. Every next pool is derived
/// from the preceding state. Publication is deferred until every runtime
/// branch has a finite, nonnegative next state.
pub fn applyAll(states: []BranchState, fluxes: []const BranchFluxes, inputs: Inputs) !void {
    if (states.len == 0 or states.len != fluxes.len)
        return error.CanopySymbioticPoolUpdateDimensionMismatch;
    if (inputs.fixation_type > 6) return error.InvalidSymbioticFixationType;
    if (inputs.fixation_type < 4) return;
    for (states, fluxes) |state, branch_fluxes| _ = try nextState(state, branch_fluxes);
    for (states, fluxes) |*state, branch_fluxes| state.* = try nextState(state.*, branch_fluxes);
}

fn nextState(state: BranchState, fluxes: BranchFluxes) !BranchState {
    try validatePool(state.nonstructural);
    try validatePool(state.structural);
    inline for (@typeInfo(BranchFluxes).@"struct".fields) |field| {
        if (field.type == Pool) {
            try validatePool(@field(fluxes, field.name));
        } else {
            const value = @field(fluxes, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidCanopySymbioticPoolUpdateFlux;
        }
    }

    const next_nonstructural = Pool{
        .carbon_g_c = state.nonstructural.carbon_g_c -
            @min(fluxes.maintenance_respiration_g_c, fluxes.oxygen_unconstrained_respiration_g_c) -
            fluxes.nitrogen_fixation_respiration_g_c - fluxes.total_growth_carbon_use_g_c +
            fluxes.recycled_decomposition.carbon_g_c,
        .nitrogen_g_n = state.nonstructural.nitrogen_g_n - fluxes.nitrogen_used_for_growth_g_n +
            fluxes.recycled_decomposition.nitrogen_g_n + fluxes.recycled_senescence.nitrogen_g_n +
            fluxes.fixed_nitrogen_g_n,
        .phosphorus_g_p = state.nonstructural.phosphorus_g_p - fluxes.phosphorus_used_for_growth_g_p +
            fluxes.recycled_decomposition.phosphorus_g_p + fluxes.recycled_senescence.phosphorus_g_p,
    };
    const next_structural = Pool{
        .carbon_g_c = state.structural.carbon_g_c + fluxes.structural_growth_g_c -
            fluxes.decomposed_structural.carbon_g_c - fluxes.senesced_structural.carbon_g_c,
        .nitrogen_g_n = state.structural.nitrogen_g_n + fluxes.nitrogen_used_for_growth_g_n -
            fluxes.decomposed_structural.nitrogen_g_n - fluxes.senesced_structural.nitrogen_g_n,
        .phosphorus_g_p = state.structural.phosphorus_g_p + fluxes.phosphorus_used_for_growth_g_p -
            fluxes.decomposed_structural.phosphorus_g_p - fluxes.senesced_structural.phosphorus_g_p,
    };
    validateNextPool(next_nonstructural) catch return error.CanopySymbioticNonstructuralPoolOverdraw;
    validateNextPool(next_structural) catch return error.CanopySymbioticStructuralPoolOverdraw;
    return .{ .nonstructural = next_nonstructural, .structural = next_structural };
}

fn validatePool(pool: Pool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopySymbioticPoolUpdateState;
}

fn validateNextPool(pool: Pool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidNextCanopySymbioticPool;
}

fn testState() BranchState {
    return .{
        .nonstructural = .{ .carbon_g_c = 5, .nitrogen_g_n = 0.4, .phosphorus_g_p = 0.08 },
        .structural = .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.2 },
    };
}

fn testFluxes() BranchFluxes {
    return .{
        .maintenance_respiration_g_c = 0.2,
        .oxygen_unconstrained_respiration_g_c = 0.3,
        .nitrogen_fixation_respiration_g_c = 0.1,
        .total_growth_carbon_use_g_c = 1,
        .recycled_decomposition = .{ .carbon_g_c = 0.2, .nitrogen_g_n = 0.02, .phosphorus_g_p = 0.004 },
        .recycled_senescence = .{ .carbon_g_c = 0.1, .nitrogen_g_n = 0.01, .phosphorus_g_p = 0.002 },
        .fixed_nitrogen_g_n = 0.025,
        .structural_growth_g_c = 0.4,
        .nitrogen_used_for_growth_g_n = 0.04,
        .phosphorus_used_for_growth_g_p = 0.008,
        .decomposed_structural = .{ .carbon_g_c = 0.5, .nitrogen_g_n = 0.05, .phosphorus_g_p = 0.01 },
        .senesced_structural = .{ .carbon_g_c = 0.25, .nitrogen_g_n = 0.025, .phosphorus_g_p = 0.005 },
    };
}

test "GROSUB publishes nonstructural then structural equations from preceding state" {
    var states = [_]BranchState{testState()};
    const fluxes = testFluxes();
    try applyAll(&states, &.{fluxes}, .{ .fixation_type = 4 });
    try std.testing.expectApproxEqAbs(@as(f64, 3.9), states[0].nonstructural.carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.415), states[0].nonstructural.nitrogen_g_n, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.078), states[0].nonstructural.phosphorus_g_p, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 9.65), states[0].structural.carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.965), states[0].structural.nitrogen_g_n, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.193), states[0].structural.phosphorus_g_p, 1e-15);
}

test "GROSUB pool update closes C N P with respiration litter and fixed atmospheric N" {
    var states = [_]BranchState{testState()};
    const before = states[0];
    const fluxes = testFluxes();
    try applyAll(&states, &.{fluxes}, .{ .fixation_type = 5 });
    const decomposition_litter = Pool{
        .carbon_g_c = fluxes.decomposed_structural.carbon_g_c - fluxes.recycled_decomposition.carbon_g_c,
        .nitrogen_g_n = fluxes.decomposed_structural.nitrogen_g_n - fluxes.recycled_decomposition.nitrogen_g_n,
        .phosphorus_g_p = fluxes.decomposed_structural.phosphorus_g_p - fluxes.recycled_decomposition.phosphorus_g_p,
    };
    const senescence_litter = Pool{
        .carbon_g_c = fluxes.senesced_structural.carbon_g_c - fluxes.recycled_senescence.carbon_g_c,
        .nitrogen_g_n = fluxes.senesced_structural.nitrogen_g_n - fluxes.recycled_senescence.nitrogen_g_n,
        .phosphorus_g_p = fluxes.senesced_structural.phosphorus_g_p - fluxes.recycled_senescence.phosphorus_g_p,
    };
    const growth_respiration_g_c = fluxes.total_growth_carbon_use_g_c - fluxes.structural_growth_g_c + fluxes.nitrogen_fixation_respiration_g_c;
    const total_respiration_g_c = @min(fluxes.maintenance_respiration_g_c, fluxes.oxygen_unconstrained_respiration_g_c) + growth_respiration_g_c + fluxes.recycled_senescence.carbon_g_c;
    try std.testing.expectApproxEqAbs(before.nonstructural.carbon_g_c + before.structural.carbon_g_c, states[0].nonstructural.carbon_g_c + states[0].structural.carbon_g_c + decomposition_litter.carbon_g_c + senescence_litter.carbon_g_c + total_respiration_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(before.nonstructural.nitrogen_g_n + before.structural.nitrogen_g_n + fluxes.fixed_nitrogen_g_n, states[0].nonstructural.nitrogen_g_n + states[0].structural.nitrogen_g_n + decomposition_litter.nitrogen_g_n + senescence_litter.nitrogen_g_n, 1e-14);
    try std.testing.expectApproxEqAbs(before.nonstructural.phosphorus_g_p + before.structural.phosphorus_g_p, states[0].nonstructural.phosphorus_g_p + states[0].structural.phosphorus_g_p + decomposition_litter.phosphorus_g_p + senescence_litter.phosphorus_g_p, 1e-14);
}

test "GROSUB root fixation gate leaves canopy pools unread and unchanged" {
    var states = [_]BranchState{testState()};
    states[0].structural.carbon_g_c = std.math.nan(f64);
    var fluxes = testFluxes();
    fluxes.fixed_nitrogen_g_n = std.math.nan(f64);
    try applyAll(&states, &.{fluxes}, .{ .fixation_type = 2 });
    try std.testing.expect(std.math.isNan(states[0].structural.carbon_g_c));
}

test "GROSUB pool update supports runtime branches and atomic late overdraw" {
    var states = [_]BranchState{testState()} ** 49;
    var fluxes = [_]BranchFluxes{testFluxes()} ** 49;
    try applyAll(&states, &fluxes, .{ .fixation_type = 6 });
    const before = states;
    fluxes[48].total_growth_carbon_use_g_c = 100;
    try std.testing.expectError(error.CanopySymbioticNonstructuralPoolOverdraw, applyAll(&states, &fluxes, .{ .fixation_type = 6 }));
    try std.testing.expectEqualDeep(before, states);
}
