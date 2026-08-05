const std = @import("std");

pub const Pool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const BranchInputs = struct {
    nonstructural: Pool,
    maintenance_respiration_g_c: f64,
    oxygen_unconstrained_respiration_g_c: f64,
    recycled_decomposition_carbon_g_c: f64,
    growth_respiration_g_c: f64,
    nitrogen_fixation_respiration_g_c: f64,
    nonstructural_nitrogen_g_n_per_g_c: f64,
    nonstructural_phosphorus_g_p_per_g_c: f64,
};

pub const Allocation = struct {
    total_growth_carbon_use_g_c: f64,
    structural_growth_g_c: f64,
    growth_and_fixation_respiration_g_c: f64,
    nitrogen_used_for_growth_g_n: f64,
    phosphorus_used_for_growth_g_p: f64,
};

pub const Inputs = struct {
    fixation_type: u8,
    growth_yield_g_c_per_g_c: f64,
    target_nitrogen_per_carbon_g_n_per_g_c: f64,
    target_phosphorus_per_carbon_g_p_per_g_c: f64,
    nitrogen_uptake_half_saturation_g_n_per_g_c: f64,
    phosphorus_uptake_half_saturation_g_p_per_g_c: f64,
};

/// Exact GROSUB lines 5506--5513 canopy diazotroph growth allocation in
/// ascending NB order. Flux units are g C, g N, or g P per biological step.
/// A negative source `CGNDL` is rejected before any runtime branch output is
/// committed because it would imply an overdraw of nonstructural carbon.
pub fn calculateAll(branches: []const BranchInputs, outputs: []Allocation, inputs: Inputs) !void {
    if (branches.len == 0 or branches.len != outputs.len)
        return error.CanopySymbioticGrowthDimensionMismatch;
    if (inputs.fixation_type > 6) return error.InvalidSymbioticFixationType;
    if (inputs.fixation_type < 4) return;
    try validateInputs(inputs);
    for (branches) |branch| _ = try calculateOne(branch, inputs);
    for (branches, outputs) |branch, *output| output.* = try calculateOne(branch, inputs);
}

fn calculateOne(branch: BranchInputs, inputs: Inputs) !Allocation {
    try validatePool(branch.nonstructural);
    inline for (@typeInfo(BranchInputs).@"struct".fields) |field| if (field.type == f64) {
        const value = @field(branch, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopySymbioticGrowthState;
    };
    if (branch.nitrogen_fixation_respiration_g_c > branch.growth_respiration_g_c)
        return error.InvalidCanopySymbioticGrowthState;

    const carbon_available = branch.nonstructural.carbon_g_c -
        @min(branch.maintenance_respiration_g_c, branch.oxygen_unconstrained_respiration_g_c) -
        branch.nitrogen_fixation_respiration_g_c + branch.recycled_decomposition_carbon_g_c;
    const respiration_limited_carbon = (branch.growth_respiration_g_c -
        branch.nitrogen_fixation_respiration_g_c) / (1 - inputs.growth_yield_g_c_per_g_c);
    const total_growth_carbon_use = @min(carbon_available, respiration_limited_carbon);
    if (!std.math.isFinite(total_growth_carbon_use) or total_growth_carbon_use < 0)
        return error.CanopySymbioticGrowthWouldOverdrawCarbon;

    const structural_growth = total_growth_carbon_use * inputs.growth_yield_g_c_per_g_c;
    const growth_and_fixation_respiration = branch.nitrogen_fixation_respiration_g_c +
        total_growth_carbon_use * (1 - inputs.growth_yield_g_c_per_g_c);
    const nitrogen_used = @max(
        0,
        @min(branch.nonstructural.nitrogen_g_n, structural_growth * inputs.target_nitrogen_per_carbon_g_n_per_g_c),
    ) * branch.nonstructural_nitrogen_g_n_per_g_c /
        (branch.nonstructural_nitrogen_g_n_per_g_c + inputs.nitrogen_uptake_half_saturation_g_n_per_g_c);
    const phosphorus_used = @max(
        0,
        @min(branch.nonstructural.phosphorus_g_p, structural_growth * inputs.target_phosphorus_per_carbon_g_p_per_g_c),
    ) * branch.nonstructural_phosphorus_g_p_per_g_c /
        (branch.nonstructural_phosphorus_g_p_per_g_c + inputs.phosphorus_uptake_half_saturation_g_p_per_g_c);
    const result = Allocation{
        .total_growth_carbon_use_g_c = total_growth_carbon_use,
        .structural_growth_g_c = structural_growth,
        .growth_and_fixation_respiration_g_c = growth_and_fixation_respiration,
        .nitrogen_used_for_growth_g_n = nitrogen_used,
        .phosphorus_used_for_growth_g_p = phosphorus_used,
    };
    inline for (@typeInfo(Allocation).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)) or @field(result, field.name) < 0)
            return error.InvalidCanopySymbioticGrowthFlux;
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (.{
        inputs.growth_yield_g_c_per_g_c,
        inputs.target_nitrogen_per_carbon_g_n_per_g_c,
        inputs.target_phosphorus_per_carbon_g_p_per_g_c,
        inputs.nitrogen_uptake_half_saturation_g_n_per_g_c,
        inputs.phosphorus_uptake_half_saturation_g_p_per_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopySymbioticGrowthInput;
    if (inputs.growth_yield_g_c_per_g_c >= 1 or
        inputs.target_nitrogen_per_carbon_g_n_per_g_c == 0 or
        inputs.target_phosphorus_per_carbon_g_p_per_g_c == 0 or
        inputs.nitrogen_uptake_half_saturation_g_n_per_g_c == 0 or
        inputs.phosphorus_uptake_half_saturation_g_p_per_g_c == 0)
        return error.InvalidCanopySymbioticGrowthInput;
}

fn validatePool(pool: Pool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopySymbioticGrowthPool;
}

fn testInputs() Inputs {
    return .{
        .fixation_type = 4,
        .growth_yield_g_c_per_g_c = 0.4,
        .target_nitrogen_per_carbon_g_n_per_g_c = 0.1,
        .target_phosphorus_per_carbon_g_p_per_g_c = 0.02,
        .nitrogen_uptake_half_saturation_g_n_per_g_c = 1.0e-4,
        .phosphorus_uptake_half_saturation_g_p_per_g_c = 1.0e-5,
    };
}

fn testBranch() BranchInputs {
    return .{
        .nonstructural = .{ .carbon_g_c = 2, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.02 },
        .maintenance_respiration_g_c = 0.01,
        .oxygen_unconstrained_respiration_g_c = 0.5,
        .recycled_decomposition_carbon_g_c = 0.02,
        .growth_respiration_g_c = 0.49,
        .nitrogen_fixation_respiration_g_c = 0.09,
        .nonstructural_nitrogen_g_n_per_g_c = 0.05,
        .nonstructural_phosphorus_g_p_per_g_c = 0.01,
    };
}

test "GROSUB CGNDL allocation preserves availability and respiration constraints" {
    var outputs: [1]Allocation = undefined;
    const branch = testBranch();
    const inputs = testInputs();
    try calculateAll(&.{branch}, &outputs, inputs);
    const availability = branch.nonstructural.carbon_g_c - @min(branch.maintenance_respiration_g_c, branch.oxygen_unconstrained_respiration_g_c) - branch.nitrogen_fixation_respiration_g_c + branch.recycled_decomposition_carbon_g_c;
    const respiration_limit = (branch.growth_respiration_g_c - branch.nitrogen_fixation_respiration_g_c) / (1 - inputs.growth_yield_g_c_per_g_c);
    try std.testing.expectEqual(@min(availability, respiration_limit), outputs[0].total_growth_carbon_use_g_c);
}

test "GROSUB bacterial growth carbon partitions into structure and respiration" {
    var outputs: [1]Allocation = undefined;
    const branch = testBranch();
    try calculateAll(&.{branch}, &outputs, testInputs());
    try std.testing.expectApproxEqAbs(
        outputs[0].total_growth_carbon_use_g_c + branch.nitrogen_fixation_respiration_g_c,
        outputs[0].structural_growth_g_c + outputs[0].growth_and_fixation_respiration_g_c,
        1e-15,
    );
    try std.testing.expect(outputs[0].nitrogen_used_for_growth_g_n <= branch.nonstructural.nitrogen_g_n);
    try std.testing.expect(outputs[0].phosphorus_used_for_growth_g_p <= branch.nonstructural.phosphorus_g_p);
}

test "GROSUB root fixation gate leaves canopy allocation unread and unchanged" {
    var outputs: [1]Allocation = undefined;
    outputs[0].structural_growth_g_c = 23;
    var inputs = testInputs();
    inputs.fixation_type = 1;
    var branch = testBranch();
    branch.nonstructural.carbon_g_c = std.math.nan(f64);
    try calculateAll(&.{branch}, &outputs, inputs);
    try std.testing.expectEqual(@as(f64, 23), outputs[0].structural_growth_g_c);
}

test "GROSUB growth allocation runtime sweep is atomic on late carbon overdraw" {
    var branches = [_]BranchInputs{testBranch()} ** 41;
    var outputs: [41]Allocation = undefined;
    try calculateAll(&branches, &outputs, testInputs());
    const before = outputs;
    branches[40].nonstructural.carbon_g_c = 0;
    branches[40].recycled_decomposition_carbon_g_c = 0;
    try std.testing.expectError(error.CanopySymbioticGrowthWouldOverdrawCarbon, calculateAll(&branches, &outputs, testInputs()));
    try std.testing.expectEqualDeep(before, outputs);
}
