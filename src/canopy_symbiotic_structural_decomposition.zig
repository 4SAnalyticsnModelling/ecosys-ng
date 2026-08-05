const std = @import("std");

pub const Pool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const BranchInputs = struct {
    structural: Pool,
    carbon_recovery_constraint_fraction: f64,
    nitrogen_recovery_constraint_fraction: f64,
    phosphorus_recovery_constraint_fraction: f64,
    bacteria_per_leaf_area_g_c_per_m2: f64,
    canopy_growth_temperature_response: f64,
    canopy_growth_water_response: f64,
};

pub const Fluxes = struct {
    carbon_recycling_fraction: f64,
    nitrogen_recycling_fraction: f64,
    phosphorus_recycling_fraction: f64,
    structural_decomposition_fraction: f64,
    decomposed: Pool,
    recycled: Pool,
    litterfall: Pool,
};

pub const Inputs = struct {
    fixation_type: u8,
    timestep_h: f64,
    specific_decomposition_rate_per_h: f64,
    bacteria_to_leaf_decomposition_control_g_c_per_g_c: f64,
    minimum_carbon_recycling_fraction: f64,
    carbon_recycling_range_fraction: f64,
    maximum_nitrogen_recycling_fraction: f64,
    maximum_phosphorus_recycling_fraction: f64,
};

/// Exact GROSUB lines 5469--5482 canopy bacterial structural decomposition
/// and recycled/litter C:N:P partitioning in ascending NB order. The source
/// caps only `SPNDL*CNDLB/CNDLI` before multiplying its environmental and
/// timestep terms. ecosys-ng fails if that final fraction would overdraw a
/// structural pool instead of silently propagating negative mass.
pub fn calculateAll(branches: []const BranchInputs, outputs: []Fluxes, inputs: Inputs) !void {
    if (branches.len == 0 or branches.len != outputs.len)
        return error.CanopySymbioticDecompositionDimensionMismatch;
    if (inputs.fixation_type > 6) return error.InvalidSymbioticFixationType;
    if (inputs.fixation_type < 4) return;
    try validateInputs(inputs);
    for (branches) |branch| _ = try calculateOne(branch, inputs);
    for (branches, outputs) |branch, *output| output.* = try calculateOne(branch, inputs);
}

fn calculateOne(branch: BranchInputs, inputs: Inputs) !Fluxes {
    try validatePool(branch.structural);
    inline for (.{
        branch.carbon_recovery_constraint_fraction,
        branch.nitrogen_recovery_constraint_fraction,
        branch.phosphorus_recovery_constraint_fraction,
        branch.bacteria_per_leaf_area_g_c_per_m2,
        branch.canopy_growth_temperature_response,
        branch.canopy_growth_water_response,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopySymbioticDecompositionState;
    inline for (.{
        branch.carbon_recovery_constraint_fraction,
        branch.nitrogen_recovery_constraint_fraction,
        branch.phosphorus_recovery_constraint_fraction,
    }) |fraction| if (fraction > 1)
        return error.InvalidCanopySymbioticDecompositionState;

    const carbon_recycling_fraction = inputs.minimum_carbon_recycling_fraction +
        branch.carbon_recovery_constraint_fraction * inputs.carbon_recycling_range_fraction;
    const nitrogen_recycling_fraction = branch.nitrogen_recovery_constraint_fraction *
        inputs.maximum_nitrogen_recycling_fraction;
    const phosphorus_recycling_fraction = branch.phosphorus_recovery_constraint_fraction *
        inputs.maximum_phosphorus_recycling_fraction;
    inline for (.{ carbon_recycling_fraction, nitrogen_recycling_fraction, phosphorus_recycling_fraction }) |fraction|
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.InvalidCanopySymbioticRecyclingFraction;

    const decomposition_fraction = @min(
        1,
        inputs.specific_decomposition_rate_per_h * branch.bacteria_per_leaf_area_g_c_per_m2 /
            inputs.bacteria_to_leaf_decomposition_control_g_c_per_g_c,
    ) * @sqrt(branch.canopy_growth_temperature_response * branch.canopy_growth_water_response) *
        inputs.timestep_h;
    if (!std.math.isFinite(decomposition_fraction) or decomposition_fraction < 0)
        return error.NonFiniteCanopySymbioticDecomposition;
    if (decomposition_fraction > 1)
        return error.CanopySymbioticDecompositionWouldOverdraw;

    const decomposed = scale(branch.structural, decomposition_fraction);
    const recycled = Pool{
        .carbon_g_c = decomposed.carbon_g_c * carbon_recycling_fraction,
        .nitrogen_g_n = decomposed.nitrogen_g_n * (nitrogen_recycling_fraction + (1 - nitrogen_recycling_fraction) * carbon_recycling_fraction),
        .phosphorus_g_p = decomposed.phosphorus_g_p * (phosphorus_recycling_fraction + (1 - phosphorus_recycling_fraction) * carbon_recycling_fraction),
    };
    const litterfall = subtract(decomposed, recycled);
    try validatePool(decomposed);
    try validatePool(recycled);
    try validatePool(litterfall);
    return .{
        .carbon_recycling_fraction = carbon_recycling_fraction,
        .nitrogen_recycling_fraction = nitrogen_recycling_fraction,
        .phosphorus_recycling_fraction = phosphorus_recycling_fraction,
        .structural_decomposition_fraction = decomposition_fraction,
        .decomposed = decomposed,
        .recycled = recycled,
        .litterfall = litterfall,
    };
}

fn validateInputs(inputs: Inputs) !void {
    inline for (.{
        inputs.timestep_h,
        inputs.specific_decomposition_rate_per_h,
        inputs.bacteria_to_leaf_decomposition_control_g_c_per_g_c,
        inputs.minimum_carbon_recycling_fraction,
        inputs.carbon_recycling_range_fraction,
        inputs.maximum_nitrogen_recycling_fraction,
        inputs.maximum_phosphorus_recycling_fraction,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopySymbioticDecompositionInput;
    if (inputs.timestep_h == 0 or inputs.bacteria_to_leaf_decomposition_control_g_c_per_g_c == 0 or
        inputs.minimum_carbon_recycling_fraction > 1 or
        inputs.minimum_carbon_recycling_fraction + inputs.carbon_recycling_range_fraction > 1 or
        inputs.maximum_nitrogen_recycling_fraction > 1 or
        inputs.maximum_phosphorus_recycling_fraction > 1)
        return error.InvalidCanopySymbioticDecompositionInput;
}

fn validatePool(pool: Pool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopySymbioticDecompositionPool;
}

fn scale(pool: Pool, fraction: f64) Pool {
    return .{
        .carbon_g_c = pool.carbon_g_c * fraction,
        .nitrogen_g_n = pool.nitrogen_g_n * fraction,
        .phosphorus_g_p = pool.phosphorus_g_p * fraction,
    };
}

fn subtract(a: Pool, b: Pool) Pool {
    return .{
        .carbon_g_c = a.carbon_g_c - b.carbon_g_c,
        .nitrogen_g_n = a.nitrogen_g_n - b.nitrogen_g_n,
        .phosphorus_g_p = a.phosphorus_g_p - b.phosphorus_g_p,
    };
}

fn testInputs() Inputs {
    return .{
        .fixation_type = 4,
        .timestep_h = 0.5,
        .specific_decomposition_rate_per_h = 0.01,
        .bacteria_to_leaf_decomposition_control_g_c_per_g_c = 0.05,
        .minimum_carbon_recycling_fraction = 0.167,
        .carbon_recycling_range_fraction = 0.333,
        .maximum_nitrogen_recycling_fraction = 0.333,
        .maximum_phosphorus_recycling_fraction = 0.333,
    };
}

fn testBranch() BranchInputs {
    return .{
        .structural = .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.2 },
        .carbon_recovery_constraint_fraction = 0.6,
        .nitrogen_recovery_constraint_fraction = 0.4,
        .phosphorus_recovery_constraint_fraction = 0.2,
        .bacteria_per_leaf_area_g_c_per_m2 = 0.5,
        .canopy_growth_temperature_response = 0.81,
        .canopy_growth_water_response = 1,
    };
}

test "GROSUB structural decomposition preserves cap then environment then timestep order" {
    var inputs = testInputs();
    inputs.specific_decomposition_rate_per_h = 1;
    var branch = testBranch();
    branch.canopy_growth_temperature_response = 0.25;
    var outputs: [1]Fluxes = undefined;
    try calculateAll(&.{branch}, &outputs, inputs);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), outputs[0].structural_decomposition_fraction, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.5), outputs[0].decomposed.carbon_g_c, 1e-15);
}

test "GROSUB decomposed C N P exactly partition between recycling and litterfall" {
    var outputs: [1]Fluxes = undefined;
    try calculateAll(&.{testBranch()}, &outputs, testInputs());
    inline for (.{ "carbon_g_c", "nitrogen_g_n", "phosphorus_g_p" }) |field_name|
        try std.testing.expectApproxEqAbs(
            @field(outputs[0].decomposed, field_name),
            @field(outputs[0].recycled, field_name) + @field(outputs[0].litterfall, field_name),
            1e-15,
        );
}

test "GROSUB canopy gate leaves root-type outputs unread and unchanged" {
    var outputs: [1]Fluxes = undefined;
    outputs[0].structural_decomposition_fraction = 7;
    var inputs = testInputs();
    inputs.fixation_type = 2;
    var branch = testBranch();
    branch.structural.carbon_g_c = std.math.nan(f64);
    try calculateAll(&.{branch}, &outputs, inputs);
    try std.testing.expectEqual(@as(f64, 7), outputs[0].structural_decomposition_fraction);
}

test "GROSUB decomposition runtime sweep is atomic on unsafe late overdraw" {
    var branches = [_]BranchInputs{testBranch()} ** 39;
    var outputs: [39]Fluxes = undefined;
    try calculateAll(&branches, &outputs, testInputs());
    const before = outputs;
    branches[38].canopy_growth_temperature_response = 25;
    branches[38].canopy_growth_water_response = 25;
    try std.testing.expectError(error.CanopySymbioticDecompositionWouldOverdraw, calculateAll(&branches, &outputs, testInputs()));
    try std.testing.expectEqualDeep(before, outputs);
}
