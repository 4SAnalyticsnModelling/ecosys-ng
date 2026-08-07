const std = @import("std");

pub const Pool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const BranchInputs = struct {
    structural: Pool,
    nonstructural: Pool,
    leaf_and_petiole_area_m2: f64,
};

pub const Constraints = struct {
    nonstructural_carbon_g_c_per_g_c: f64,
    nonstructural_nitrogen_g_n_per_g_c: f64,
    nonstructural_phosphorus_g_p_per_g_c: f64,
    nitrogen_per_nonstructural_carbon_g_n_per_g_c: f64,
    nitrogen_per_nonstructural_phosphorus_g_n_per_g_p: f64,
    carbon_recovery_fraction: f64,
    nitrogen_recovery_constraint_fraction: f64,
    phosphorus_recovery_constraint_fraction: f64,
    nutrient_activity_fraction: f64,
    bacteria_per_leaf_area_g_c_per_m2: f64,
};

pub const Inputs = struct {
    fixation_type: u8,
    structural_presence_threshold_g_c: f64,
    ratio_division_threshold: f64,
    leaf_area_presence_threshold_m2: f64,
    nitrogen_inhibition_g_n_per_g_c: f64,
    phosphorus_inhibition_g_p_per_g_c: f64,
    target_nitrogen_per_carbon_g_n_per_g_c: f64,
    target_phosphorus_per_carbon_g_p_per_g_c: f64,
};

/// Exact grosub.f lines 5312--5373 canopy bacterial activity and recycling
/// constraints in ascending NB order. `CNKI`/`CPKI` remain distinct from the
/// structural `CNND`/`CPND` targets. The read-only sweep commits outputs only
/// after every active branch has passed validation and calculation.
pub fn calculateAll(branches: []const BranchInputs, outputs: []Constraints, inputs: Inputs) !void {
    if (branches.len == 0 or outputs.len != branches.len)
        return error.CanopySymbioticConstraintDimensionMismatch;
    if (inputs.fixation_type > 6) return error.InvalidSymbioticFixationType;
    if (inputs.fixation_type < 4) return;
    try validateInputs(inputs);
    for (branches) |branch| _ = try calculateOne(branch, inputs);
    for (branches, outputs) |branch, *output| output.* = try calculateOne(branch, inputs);
}

fn calculateOne(branch: BranchInputs, inputs: Inputs) !Constraints {
    try validatePool(branch.structural);
    try validatePool(branch.nonstructural);
    if (!std.math.isFinite(branch.leaf_and_petiole_area_m2) or branch.leaf_and_petiole_area_m2 < 0)
        return error.InvalidCanopySymbioticConstraintState;

    var carbon_concentration: f64 = 1;
    var nitrogen_concentration: f64 = 1;
    var phosphorus_concentration: f64 = 1;
    if (branch.structural.carbon_g_c > inputs.structural_presence_threshold_g_c) {
        carbon_concentration = @max(0, branch.nonstructural.carbon_g_c / branch.structural.carbon_g_c);
        nitrogen_concentration = @max(0, branch.nonstructural.nitrogen_g_n / branch.structural.carbon_g_c);
        phosphorus_concentration = @max(0, branch.nonstructural.phosphorus_g_p / branch.structural.carbon_g_c);
    }
    const nitrogen_per_carbon = if (carbon_concentration > inputs.ratio_division_threshold)
        nitrogen_concentration / carbon_concentration
    else
        0;
    const nitrogen_per_phosphorus = if (phosphorus_concentration > inputs.ratio_division_threshold)
        nitrogen_concentration / phosphorus_concentration
    else
        0;

    var carbon_recovery: f64 = 1;
    var nitrogen_recovery: f64 = 0;
    var phosphorus_recovery: f64 = 0;
    if (carbon_concentration > inputs.ratio_division_threshold) {
        carbon_recovery = @max(0, @min(
            1,
            nitrogen_concentration / (nitrogen_concentration + carbon_concentration * inputs.nitrogen_inhibition_g_n_per_g_c),
            phosphorus_concentration / (phosphorus_concentration + carbon_concentration * inputs.phosphorus_inhibition_g_p_per_g_c),
        ));
        nitrogen_recovery = @max(0, @min(1, carbon_concentration / (carbon_concentration + nitrogen_concentration / inputs.nitrogen_inhibition_g_n_per_g_c)));
        phosphorus_recovery = @max(0, @min(1, carbon_concentration / (carbon_concentration + phosphorus_concentration / inputs.phosphorus_inhibition_g_p_per_g_c)));
    }

    const nutrient_activity = if (branch.structural.carbon_g_c > inputs.structural_presence_threshold_g_c)
        @min(
            1,
            @sqrt(branch.structural.nitrogen_g_n / (branch.structural.carbon_g_c * inputs.target_nitrogen_per_carbon_g_n_per_g_c)),
            @sqrt(branch.structural.phosphorus_g_p / (branch.structural.carbon_g_c * inputs.target_phosphorus_per_carbon_g_p_per_g_c)),
        )
    else
        1;
    const decomposition_density = if (branch.leaf_and_petiole_area_m2 > inputs.leaf_area_presence_threshold_m2)
        @max(0, branch.structural.carbon_g_c / branch.leaf_and_petiole_area_m2)
    else
        0;
    const result = Constraints{
        .nonstructural_carbon_g_c_per_g_c = carbon_concentration,
        .nonstructural_nitrogen_g_n_per_g_c = nitrogen_concentration,
        .nonstructural_phosphorus_g_p_per_g_c = phosphorus_concentration,
        .nitrogen_per_nonstructural_carbon_g_n_per_g_c = nitrogen_per_carbon,
        .nitrogen_per_nonstructural_phosphorus_g_n_per_g_p = nitrogen_per_phosphorus,
        .carbon_recovery_fraction = carbon_recovery,
        .nitrogen_recovery_constraint_fraction = nitrogen_recovery,
        .phosphorus_recovery_constraint_fraction = phosphorus_recovery,
        .nutrient_activity_fraction = nutrient_activity,
        .bacteria_per_leaf_area_g_c_per_m2 = decomposition_density,
    };
    inline for (@typeInfo(Constraints).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteCanopySymbioticConstraint;
    return result;
}

fn validateInputs(inputs: Inputs) !void {
    inline for (.{
        inputs.structural_presence_threshold_g_c,
        inputs.ratio_division_threshold,
        inputs.leaf_area_presence_threshold_m2,
        inputs.nitrogen_inhibition_g_n_per_g_c,
        inputs.phosphorus_inhibition_g_p_per_g_c,
        inputs.target_nitrogen_per_carbon_g_n_per_g_c,
        inputs.target_phosphorus_per_carbon_g_p_per_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopySymbioticConstraintInput;
    if (inputs.nitrogen_inhibition_g_n_per_g_c == 0 or
        inputs.phosphorus_inhibition_g_p_per_g_c == 0 or
        inputs.target_nitrogen_per_carbon_g_n_per_g_c == 0 or
        inputs.target_phosphorus_per_carbon_g_p_per_g_c == 0)
        return error.InvalidCanopySymbioticConstraintInput;
}

fn validatePool(pool: Pool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopySymbioticConstraintState;
}

fn testInputs() Inputs {
    return .{
        .fixation_type = 4,
        .structural_presence_threshold_g_c = 1.0e-12,
        .ratio_division_threshold = 1.0e-20,
        .leaf_area_presence_threshold_m2 = 1.0e-10,
        .nitrogen_inhibition_g_n_per_g_c = 0.1,
        .phosphorus_inhibition_g_p_per_g_c = 0.01,
        .target_nitrogen_per_carbon_g_n_per_g_c = 0.1,
        .target_phosphorus_per_carbon_g_p_per_g_c = 0.02,
    };
}

fn testBranch() BranchInputs {
    return .{
        .structural = .{ .carbon_g_c = 10, .nitrogen_g_n = 0.5, .phosphorus_g_p = 0.1 },
        .nonstructural = .{ .carbon_g_c = 2, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.02 },
        .leaf_and_petiole_area_m2 = 5,
    };
}

test "GROSUB bacterial concentrations ratios activity and density preserve equation order" {
    var outputs: [1]Constraints = undefined;
    try calculateAll(&.{testBranch()}, &outputs, testInputs());
    const result = outputs[0];
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), result.nonstructural_carbon_g_c_per_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.01), result.nonstructural_nitrogen_g_n_per_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.002), result.nonstructural_phosphorus_g_p_per_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), result.nitrogen_per_nonstructural_carbon_g_n_per_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 5), result.nitrogen_per_nonstructural_phosphorus_g_n_per_g_p, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.bacteria_per_leaf_area_g_c_per_m2, 1e-15);
    try std.testing.expectApproxEqAbs(@sqrt(@as(f64, 0.5)), result.nutrient_activity_fraction, 1e-15);
}

test "GROSUB CNKI CPKI are not structural bacterial target ratios" {
    var outputs: [1]Constraints = undefined;
    var exact = testInputs();
    try calculateAll(&.{testBranch()}, &outputs, exact);
    const exact_carbon_recovery = outputs[0].carbon_recovery_fraction;
    exact.nitrogen_inhibition_g_n_per_g_c = exact.target_nitrogen_per_carbon_g_n_per_g_c * 0.25;
    exact.phosphorus_inhibition_g_p_per_g_c = exact.target_phosphorus_per_carbon_g_p_per_g_c * 0.25;
    try calculateAll(&.{testBranch()}, &outputs, exact);
    try std.testing.expect(outputs[0].carbon_recovery_fraction > exact_carbon_recovery);
}

test "GROSUB below structural threshold uses unit concentrations and activity" {
    var low = testBranch();
    low.structural.carbon_g_c = 1.0e-13;
    var outputs: [1]Constraints = undefined;
    try calculateAll(&.{low}, &outputs, testInputs());
    try std.testing.expectEqual(@as(f64, 1), outputs[0].nonstructural_carbon_g_c_per_g_c);
    try std.testing.expectEqual(@as(f64, 1), outputs[0].nonstructural_nitrogen_g_n_per_g_c);
    try std.testing.expectEqual(@as(f64, 1), outputs[0].nonstructural_phosphorus_g_p_per_g_c);
    try std.testing.expectEqual(@as(f64, 1), outputs[0].nutrient_activity_fraction);
}

test "GROSUB root fixation gate leaves canopy constraint outputs unread and unchanged" {
    var outputs = [_]Constraints{undefined};
    outputs[0].carbon_recovery_fraction = 17;
    var inputs = testInputs();
    inputs.fixation_type = 3;
    const invalid = BranchInputs{ .structural = .{ .carbon_g_c = std.math.nan(f64), .nitrogen_g_n = 0, .phosphorus_g_p = 0 }, .nonstructural = .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 }, .leaf_and_petiole_area_m2 = 0 };
    try calculateAll(&.{invalid}, &outputs, inputs);
    try std.testing.expectEqual(@as(f64, 17), outputs[0].carbon_recovery_fraction);
}

test "GROSUB constraint sweep supports runtime branches and atomic late failure" {
    var branches = [_]BranchInputs{testBranch()} ** 37;
    var outputs: [37]Constraints = undefined;
    try calculateAll(&branches, &outputs, testInputs());
    try std.testing.expectApproxEqAbs(outputs[0].carbon_recovery_fraction, outputs[36].carbon_recovery_fraction, 0);

    const before = outputs;
    branches[36].nonstructural.nitrogen_g_n = std.math.nan(f64);
    try std.testing.expectError(error.InvalidCanopySymbioticConstraintState, calculateAll(&branches, &outputs, testInputs()));
    try std.testing.expectEqualDeep(before, outputs);
}
