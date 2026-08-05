const std = @import("std");

pub const Pool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const Inputs = struct {
    /// INTYP: 0 disables fixation, 1--3 are root symbioses, 4--6 canopy.
    fixation_type: u8,
    /// NFZ == 1, the first biological subhour.
    first_biological_subhour: bool,
    /// ICHKF != 0 while restoring a checkpoint.
    restoring_checkpoint: bool,
    initial_bacterial_carbon_g_c_per_m2: f64,
    horizontal_cell_area_m2: f64,
    target_nitrogen_per_carbon_g_n_per_g_c: f64,
    target_phosphorus_per_carbon_g_p_per_g_c: f64,
};

pub const Result = struct {
    infected_branch_count: usize,
    introduced: Pool,
};

/// Exact GROSUB lines 5283--5301 canopy cyanobacteria infection sweep.
/// The introduced pool is external inoculum, not a transfer from host tissue.
/// All eligible branches are validated before any branch is mutated.
pub fn apply(structural_by_branch: []Pool, inputs: Inputs) !Result {
    if (structural_by_branch.len == 0) return error.CanopySymbioticInfectionDimensionMismatch;
    if (inputs.fixation_type > 6) return error.InvalidSymbioticFixationType;
    if (inputs.fixation_type < 4 or
        !inputs.first_biological_subhour or inputs.restoring_checkpoint)
        return .{ .infected_branch_count = 0, .introduced = zeroPool() };

    inline for (.{
        inputs.initial_bacterial_carbon_g_c_per_m2,
        inputs.horizontal_cell_area_m2,
        inputs.target_nitrogen_per_carbon_g_n_per_g_c,
        inputs.target_phosphorus_per_carbon_g_p_per_g_c,
    }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteCanopySymbioticInfectionInput;
        if (value < 0) return error.NegativeCanopySymbioticInfectionInput;
    }
    if (inputs.horizontal_cell_area_m2 <= 0)
        return error.InvalidCanopySymbioticInfectionInput;

    const carbon_g_c = inputs.initial_bacterial_carbon_g_c_per_m2 * inputs.horizontal_cell_area_m2;
    const inoculum = Pool{
        .carbon_g_c = carbon_g_c,
        .nitrogen_g_n = carbon_g_c * inputs.target_nitrogen_per_carbon_g_n_per_g_c,
        .phosphorus_g_p = carbon_g_c * inputs.target_phosphorus_per_carbon_g_p_per_g_c,
    };
    try validatePool(inoculum);

    for (structural_by_branch) |pool| {
        try validatePool(pool);
        if (pool.carbon_g_c > 0) continue;
        // Safe initialized state turns the source WTNDB <= 0 gate into zero.
        if (pool.carbon_g_c != 0) return error.NegativeCanopySymbioticStructuralCarbon;
        const next = add(pool, inoculum);
        try validatePool(next);
    }

    var result = Result{ .infected_branch_count = 0, .introduced = zeroPool() };
    for (structural_by_branch) |*pool| {
        if (pool.carbon_g_c > 0) continue;
        pool.* = add(pool.*, inoculum);
        result.infected_branch_count += 1;
        result.introduced = add(result.introduced, inoculum);
    }
    try validatePool(result.introduced);
    return result;
}

fn zeroPool() Pool {
    return .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
}

fn add(a: Pool, b: Pool) Pool {
    return .{
        .carbon_g_c = a.carbon_g_c + b.carbon_g_c,
        .nitrogen_g_n = a.nitrogen_g_n + b.nitrogen_g_n,
        .phosphorus_g_p = a.phosphorus_g_p + b.phosphorus_g_p,
    };
}

fn validatePool(pool: Pool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteCanopySymbioticPool;
        if (value < 0) return error.NegativeCanopySymbioticPool;
    }
}

fn testInputs() Inputs {
    return .{
        .fixation_type = 4,
        .first_biological_subhour = true,
        .restoring_checkpoint = false,
        .initial_bacterial_carbon_g_c_per_m2 = 0.02,
        .horizontal_cell_area_m2 = 100,
        .target_nitrogen_per_carbon_g_n_per_g_c = 0.1,
        .target_phosphorus_per_carbon_g_p_per_g_c = 0.02,
    };
}

test "GROSUB canopy infection adds area-scaled inoculum at exact C N P ratios" {
    var branches = [_]Pool{ zeroPool(), .{ .carbon_g_c = 3, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.04 }, zeroPool() };
    const before = branches;
    const result = try apply(&branches, testInputs());
    try std.testing.expectEqual(@as(usize, 2), result.infected_branch_count);
    try std.testing.expectEqual(Pool{ .carbon_g_c = 4, .nitrogen_g_n = 0.4, .phosphorus_g_p = 0.08 }, result.introduced);
    try std.testing.expectEqual(Pool{ .carbon_g_c = 2, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.04 }, branches[0]);
    try std.testing.expectEqual(before[1], branches[1]);
    try std.testing.expectEqual(branches[0], branches[2]);
    try std.testing.expectEqual(before[0].carbon_g_c + result.introduced.carbon_g_c / 2, branches[0].carbon_g_c);
}

test "GROSUB canopy infection gate excludes root types later subhours and checkpoints" {
    const nan_pool = Pool{ .carbon_g_c = std.math.nan(f64), .nitrogen_g_n = std.math.nan(f64), .phosphorus_g_p = std.math.nan(f64) };
    var root_type = [_]Pool{nan_pool};
    var root_inputs = testInputs();
    root_inputs.fixation_type = 3;
    try std.testing.expectEqual(@as(usize, 0), (try apply(&root_type, root_inputs)).infected_branch_count);
    try std.testing.expect(std.math.isNan(root_type[0].carbon_g_c));

    var later = [_]Pool{nan_pool};
    var later_inputs = testInputs();
    later_inputs.first_biological_subhour = false;
    try std.testing.expectEqual(@as(usize, 0), (try apply(&later, later_inputs)).infected_branch_count);

    var checkpoint = [_]Pool{nan_pool};
    var checkpoint_inputs = testInputs();
    checkpoint_inputs.restoring_checkpoint = true;
    try std.testing.expectEqual(@as(usize, 0), (try apply(&checkpoint, checkpoint_inputs)).infected_branch_count);
}

test "GROSUB canopy infection supports arbitrary runtime branch counts" {
    var branches = [_]Pool{zeroPool()} ** 47;
    const result = try apply(&branches, testInputs());
    try std.testing.expectEqual(@as(usize, 47), result.infected_branch_count);
    try std.testing.expectEqual(@as(f64, 94), result.introduced.carbon_g_c);
}

test "GROSUB invalid late canopy branch leaves the complete sweep unchanged" {
    var branches = [_]Pool{ zeroPool(), .{ .carbon_g_c = 0, .nitrogen_g_n = std.math.nan(f64), .phosphorus_g_p = 0 } };
    try std.testing.expectError(error.NonFiniteCanopySymbioticPool, apply(&branches, testInputs()));
    try std.testing.expectEqual(zeroPool(), branches[0]);
    try std.testing.expect(std.math.isNan(branches[1].nitrogen_g_n));
}

test "GROSUB infection reports introduced mass as exact state delta" {
    var branches = [_]Pool{ zeroPool(), .{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.02 } };
    const before = branches;
    const result = try apply(&branches, testInputs());
    const delta = Pool{
        .carbon_g_c = branches[0].carbon_g_c + branches[1].carbon_g_c - before[0].carbon_g_c - before[1].carbon_g_c,
        .nitrogen_g_n = branches[0].nitrogen_g_n + branches[1].nitrogen_g_n - before[0].nitrogen_g_n - before[1].nitrogen_g_n,
        .phosphorus_g_p = branches[0].phosphorus_g_p + branches[1].phosphorus_g_p - before[0].phosphorus_g_p - before[1].phosphorus_g_p,
    };
    try std.testing.expectApproxEqAbs(result.introduced.carbon_g_c, delta.carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(result.introduced.nitrogen_g_n, delta.nitrogen_g_n, 1e-15);
    try std.testing.expectApproxEqAbs(result.introduced.phosphorus_g_p, delta.phosphorus_g_p, 1e-15);
}
