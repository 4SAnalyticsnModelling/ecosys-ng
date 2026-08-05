const std = @import("std");

pub const Pool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const BranchInputs = struct {
    decomposition_litterfall: Pool,
    senescence_litterfall: Pool,
};

pub const KineticFractions = struct {
    carbon: []const f64,
    nitrogen: []const f64,
    phosphorus: []const f64,
};

pub const LitterState = struct {
    carbon_g_c_by_kinetic_class: []f64,
    nitrogen_g_n_by_kinetic_class: []f64,
    phosphorus_g_p_by_kinetic_class: []f64,
};

pub const Inputs = struct {
    fixation_type: u8,
    foliar_structural_fractions: KineticFractions,
    fraction_sum_tolerance: f64,
};

/// Exact GROSUB lines 5581--5592 canopy symbiont structural litter routing.
/// Runtime kinetic indexes correspond to source M=1..4; no array extent is
/// compiled into this kernel. Branches remain the outer loop and kinetic
/// classes the inner loop, preserving the source accumulation order.
pub fn routeAll(branches: []const BranchInputs, state: LitterState, inputs: Inputs) !void {
    const kinetic_count = state.carbon_g_c_by_kinetic_class.len;
    if (branches.len == 0 or kinetic_count == 0 or
        state.nitrogen_g_n_by_kinetic_class.len != kinetic_count or
        state.phosphorus_g_p_by_kinetic_class.len != kinetic_count or
        inputs.foliar_structural_fractions.carbon.len != kinetic_count or
        inputs.foliar_structural_fractions.nitrogen.len != kinetic_count or
        inputs.foliar_structural_fractions.phosphorus.len != kinetic_count)
        return error.CanopySymbioticLitterDimensionMismatch;
    if (inputs.fixation_type > 6) return error.InvalidSymbioticFixationType;
    if (inputs.fixation_type < 4) return;
    if (!std.math.isFinite(inputs.fraction_sum_tolerance) or inputs.fraction_sum_tolerance < 0)
        return error.InvalidCanopySymbioticLitterInput;

    try validateFractions(inputs.foliar_structural_fractions, inputs.fraction_sum_tolerance);
    inline for (.{
        state.carbon_g_c_by_kinetic_class,
        state.nitrogen_g_n_by_kinetic_class,
        state.phosphorus_g_p_by_kinetic_class,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopySymbioticLitterState;
    for (branches) |branch| {
        try validatePool(branch.decomposition_litterfall);
        try validatePool(branch.senescence_litterfall);
    }

    // Preflight each owned kinetic pool using the same ascending branch sum.
    for (0..kinetic_count) |kinetic| {
        var next_carbon = state.carbon_g_c_by_kinetic_class[kinetic];
        var next_nitrogen = state.nitrogen_g_n_by_kinetic_class[kinetic];
        var next_phosphorus = state.phosphorus_g_p_by_kinetic_class[kinetic];
        for (branches) |branch| {
            next_carbon += inputs.foliar_structural_fractions.carbon[kinetic] *
                (branch.decomposition_litterfall.carbon_g_c + branch.senescence_litterfall.carbon_g_c);
            next_nitrogen += inputs.foliar_structural_fractions.nitrogen[kinetic] *
                (branch.decomposition_litterfall.nitrogen_g_n + branch.senescence_litterfall.nitrogen_g_n);
            next_phosphorus += inputs.foliar_structural_fractions.phosphorus[kinetic] *
                (branch.decomposition_litterfall.phosphorus_g_p + branch.senescence_litterfall.phosphorus_g_p);
            inline for (.{ next_carbon, next_nitrogen, next_phosphorus }) |value|
                if (!std.math.isFinite(value)) return error.NonFiniteCanopySymbioticLitterState;
        }
    }

    for (branches) |branch| {
        for (0..kinetic_count) |kinetic| {
            state.carbon_g_c_by_kinetic_class[kinetic] +=
                inputs.foliar_structural_fractions.carbon[kinetic] *
                (branch.decomposition_litterfall.carbon_g_c + branch.senescence_litterfall.carbon_g_c);
            state.nitrogen_g_n_by_kinetic_class[kinetic] +=
                inputs.foliar_structural_fractions.nitrogen[kinetic] *
                (branch.decomposition_litterfall.nitrogen_g_n + branch.senescence_litterfall.nitrogen_g_n);
            state.phosphorus_g_p_by_kinetic_class[kinetic] +=
                inputs.foliar_structural_fractions.phosphorus[kinetic] *
                (branch.decomposition_litterfall.phosphorus_g_p + branch.senescence_litterfall.phosphorus_g_p);
        }
    }
}

fn validateFractions(fractions: KineticFractions, tolerance: f64) !void {
    inline for (.{ fractions.carbon, fractions.nitrogen, fractions.phosphorus }) |values| {
        var sum: f64 = 0;
        for (values) |value| {
            if (!std.math.isFinite(value) or value < 0 or value > 1)
                return error.InvalidCanopySymbioticLitterFraction;
            sum += value;
        }
        if (!std.math.isFinite(sum) or @abs(sum - 1) > tolerance)
            return error.InvalidCanopySymbioticLitterFractionSum;
    }
}

fn validatePool(pool: Pool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopySymbioticLitterPool;
}

fn testBranch(scale: f64) BranchInputs {
    return .{
        .decomposition_litterfall = .{ .carbon_g_c = 2 * scale, .nitrogen_g_n = 0.2 * scale, .phosphorus_g_p = 0.02 * scale },
        .senescence_litterfall = .{ .carbon_g_c = 1 * scale, .nitrogen_g_n = 0.1 * scale, .phosphorus_g_p = 0.01 * scale },
    };
}

test "GROSUB structural litter routes decomposition plus senescence by element" {
    var carbon = [_]f64{ 0, 0, 0, 0 };
    var nitrogen = [_]f64{ 0, 0, 0, 0 };
    var phosphorus = [_]f64{ 0, 0, 0, 0 };
    const fractions = [_]f64{ 0.1, 0.2, 0.3, 0.4 };
    try routeAll(&.{testBranch(1)}, .{
        .carbon_g_c_by_kinetic_class = &carbon,
        .nitrogen_g_n_by_kinetic_class = &nitrogen,
        .phosphorus_g_p_by_kinetic_class = &phosphorus,
    }, .{
        .fixation_type = 4,
        .foliar_structural_fractions = .{ .carbon = &fractions, .nitrogen = &fractions, .phosphorus = &fractions },
        .fraction_sum_tolerance = 1e-12,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), carbon[3], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.12), nitrogen[3], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.012), phosphorus[3], 1e-15);
}

test "GROSUB kinetic routing conserves structural litter C N P" {
    var carbon = [_]f64{ 0, 0, 0 };
    var nitrogen = [_]f64{ 0, 0, 0 };
    var phosphorus = [_]f64{ 0, 0, 0 };
    const fractions = [_]f64{ 0.2, 0.3, 0.5 };
    const branches = [_]BranchInputs{ testBranch(1), testBranch(2) };
    try routeAll(&branches, .{ .carbon_g_c_by_kinetic_class = &carbon, .nitrogen_g_n_by_kinetic_class = &nitrogen, .phosphorus_g_p_by_kinetic_class = &phosphorus }, .{ .fixation_type = 5, .foliar_structural_fractions = .{ .carbon = &fractions, .nitrogen = &fractions, .phosphorus = &fractions }, .fraction_sum_tolerance = 1e-12 });
    var routed = Pool{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
    for (carbon, nitrogen, phosphorus) |c, n, p| {
        routed.carbon_g_c += c;
        routed.nitrogen_g_n += n;
        routed.phosphorus_g_p += p;
    }
    try std.testing.expectApproxEqAbs(@as(f64, 9), routed.carbon_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.9), routed.nitrogen_g_n, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.09), routed.phosphorus_g_p, 1e-14);
}

test "GROSUB structural litter supports runtime branches and kinetic topology atomically" {
    var branches = [_]BranchInputs{testBranch(1)} ** 43;
    var carbon = [_]f64{ 1, 2, 3, 4, 5 };
    var nitrogen = [_]f64{ 1, 2, 3, 4, 5 };
    var phosphorus = [_]f64{ 1, 2, 3, 4, 5 };
    const fractions = [_]f64{ 0.1, 0.2, 0.25, 0.15, 0.3 };
    const state = LitterState{ .carbon_g_c_by_kinetic_class = &carbon, .nitrogen_g_n_by_kinetic_class = &nitrogen, .phosphorus_g_p_by_kinetic_class = &phosphorus };
    const inputs = Inputs{ .fixation_type = 6, .foliar_structural_fractions = .{ .carbon = &fractions, .nitrogen = &fractions, .phosphorus = &fractions }, .fraction_sum_tolerance = 1e-12 };
    try routeAll(&branches, state, inputs);
    const before_carbon = carbon;
    const before_nitrogen = nitrogen;
    const before_phosphorus = phosphorus;
    branches[42].senescence_litterfall.nitrogen_g_n = std.math.nan(f64);
    try std.testing.expectError(error.InvalidCanopySymbioticLitterPool, routeAll(&branches, state, inputs));
    try std.testing.expectEqualDeep(before_carbon, carbon);
    try std.testing.expectEqualDeep(before_nitrogen, nitrogen);
    try std.testing.expectEqualDeep(before_phosphorus, phosphorus);
}

test "GROSUB root fixation gate leaves structural litter state unread" {
    var carbon = [_]f64{std.math.nan(f64)};
    var nitrogen = [_]f64{std.math.nan(f64)};
    var phosphorus = [_]f64{std.math.nan(f64)};
    const nan_fraction = [_]f64{std.math.nan(f64)};
    try routeAll(&.{testBranch(1)}, .{ .carbon_g_c_by_kinetic_class = &carbon, .nitrogen_g_n_by_kinetic_class = &nitrogen, .phosphorus_g_p_by_kinetic_class = &phosphorus }, .{ .fixation_type = 3, .foliar_structural_fractions = .{ .carbon = &nan_fraction, .nitrogen = &nan_fraction, .phosphorus = &nan_fraction }, .fraction_sum_tolerance = std.math.nan(f64) });
    try std.testing.expect(std.math.isNan(carbon[0]));
}
