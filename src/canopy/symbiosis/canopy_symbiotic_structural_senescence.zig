const std = @import("std");

pub const Pool = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
};

pub const BranchInputs = struct {
    structural: Pool,
    excess_maintenance_respiration_g_c: f64,
    carbon_recycling_fraction: f64,
    nitrogen_recycling_fraction: f64,
    phosphorus_recycling_fraction: f64,
};

pub const Fluxes = struct {
    senesced: Pool,
    recycled: Pool,
    litterfall: Pool,
};

pub const Inputs = struct {
    fixation_type: u8,
    structural_presence_threshold_g_c: f64,
};

/// Exact grosub.f lines 5527--5547 maintenance-driven canopy bacterial
/// senescence in ascending NB order. Fluxes are g C, g N, and g P per
/// biological timestep. A senescence C demand larger than structural C fails
/// atomically rather than producing negative structural mass downstream.
pub fn calculateAll(branches: []const BranchInputs, outputs: []Fluxes, inputs: Inputs) !void {
    if (branches.len == 0 or branches.len != outputs.len)
        return error.CanopySymbioticSenescenceDimensionMismatch;
    if (inputs.fixation_type > 6) return error.InvalidSymbioticFixationType;
    if (inputs.fixation_type < 4) return;
    if (!std.math.isFinite(inputs.structural_presence_threshold_g_c) or
        inputs.structural_presence_threshold_g_c < 0)
        return error.InvalidCanopySymbioticSenescenceInput;
    for (branches) |branch| _ = try calculateOne(branch, inputs);
    for (branches, outputs) |branch, *output| output.* = try calculateOne(branch, inputs);
}

fn calculateOne(branch: BranchInputs, inputs: Inputs) !Fluxes {
    try validatePool(branch.structural);
    inline for (.{
        branch.excess_maintenance_respiration_g_c,
        branch.carbon_recycling_fraction,
        branch.nitrogen_recycling_fraction,
        branch.phosphorus_recycling_fraction,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCanopySymbioticSenescenceState;
    inline for (.{
        branch.carbon_recycling_fraction,
        branch.nitrogen_recycling_fraction,
        branch.phosphorus_recycling_fraction,
    }) |fraction| if (fraction > 1)
        return error.InvalidCanopySymbioticSenescenceState;

    if (!(branch.excess_maintenance_respiration_g_c > 0 and
        branch.structural.carbon_g_c > inputs.structural_presence_threshold_g_c))
        return zeroFluxes();
    if (branch.excess_maintenance_respiration_g_c > branch.structural.carbon_g_c)
        return error.CanopySymbioticSenescenceWouldOverdraw;

    const senesced = Pool{
        .carbon_g_c = branch.excess_maintenance_respiration_g_c,
        .nitrogen_g_n = branch.excess_maintenance_respiration_g_c *
            branch.structural.nitrogen_g_n / branch.structural.carbon_g_c,
        .phosphorus_g_p = branch.excess_maintenance_respiration_g_c *
            branch.structural.phosphorus_g_p / branch.structural.carbon_g_c,
    };
    const recycled = Pool{
        .carbon_g_c = senesced.carbon_g_c * branch.carbon_recycling_fraction,
        .nitrogen_g_n = senesced.nitrogen_g_n * (branch.nitrogen_recycling_fraction +
            (1 - branch.nitrogen_recycling_fraction) * branch.carbon_recycling_fraction),
        .phosphorus_g_p = senesced.phosphorus_g_p * (branch.phosphorus_recycling_fraction +
            (1 - branch.phosphorus_recycling_fraction) * branch.carbon_recycling_fraction),
    };
    const litterfall = subtract(senesced, recycled);
    try validatePool(senesced);
    try validatePool(recycled);
    try validatePool(litterfall);
    return .{ .senesced = senesced, .recycled = recycled, .litterfall = litterfall };
}

fn zeroPool() Pool {
    return .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
}

fn zeroFluxes() Fluxes {
    return .{ .senesced = zeroPool(), .recycled = zeroPool(), .litterfall = zeroPool() };
}

fn subtract(a: Pool, b: Pool) Pool {
    return .{
        .carbon_g_c = a.carbon_g_c - b.carbon_g_c,
        .nitrogen_g_n = a.nitrogen_g_n - b.nitrogen_g_n,
        .phosphorus_g_p = a.phosphorus_g_p - b.phosphorus_g_p,
    };
}

fn validatePool(pool: Pool) !void {
    inline for (.{ pool.carbon_g_c, pool.nitrogen_g_n, pool.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopySymbioticSenescencePool;
}

fn testInputs() Inputs {
    return .{ .fixation_type = 4, .structural_presence_threshold_g_c = 1.0e-12 };
}

fn testBranch() BranchInputs {
    return .{
        .structural = .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.2 },
        .excess_maintenance_respiration_g_c = 2,
        .carbon_recycling_fraction = 0.4,
        .nitrogen_recycling_fraction = 0.3,
        .phosphorus_recycling_fraction = 0.2,
    };
}

test "GROSUB senescence removes structural N and P in structural C proportion" {
    var outputs: [1]Fluxes = undefined;
    try calculateAll(&.{testBranch()}, &outputs, testInputs());
    try std.testing.expectEqual(@as(f64, 2), outputs[0].senesced.carbon_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), outputs[0].senesced.nitrogen_g_n, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.04), outputs[0].senesced.phosphorus_g_p, 1e-15);
}

test "GROSUB senesced C N P exactly partition into recycling and litterfall" {
    var outputs: [1]Fluxes = undefined;
    try calculateAll(&.{testBranch()}, &outputs, testInputs());
    inline for (.{ "carbon_g_c", "nitrogen_g_n", "phosphorus_g_p" }) |field_name|
        try std.testing.expectApproxEqAbs(
            @field(outputs[0].senesced, field_name),
            @field(outputs[0].recycled, field_name) + @field(outputs[0].litterfall, field_name),
            1e-15,
        );
}

test "GROSUB zero excess or threshold structural carbon writes zero fluxes" {
    var branches = [_]BranchInputs{ testBranch(), testBranch() };
    branches[0].excess_maintenance_respiration_g_c = 0;
    branches[1].structural.carbon_g_c = 1.0e-12;
    var outputs: [2]Fluxes = undefined;
    try calculateAll(&branches, &outputs, testInputs());
    try std.testing.expectEqualDeep(zeroFluxes(), outputs[0]);
    try std.testing.expectEqualDeep(zeroFluxes(), outputs[1]);
}

test "GROSUB root fixation gate leaves canopy senescence output unchanged" {
    var outputs: [1]Fluxes = undefined;
    outputs[0].senesced.carbon_g_c = 29;
    var inputs = testInputs();
    inputs.fixation_type = 3;
    var branch = testBranch();
    branch.structural.carbon_g_c = std.math.nan(f64);
    try calculateAll(&.{branch}, &outputs, inputs);
    try std.testing.expectEqual(@as(f64, 29), outputs[0].senesced.carbon_g_c);
}

test "GROSUB senescence runtime sweep is atomic on late structural overdraw" {
    var branches = [_]BranchInputs{testBranch()} ** 45;
    var outputs: [45]Fluxes = undefined;
    try calculateAll(&branches, &outputs, testInputs());
    const before = outputs;
    branches[44].excess_maintenance_respiration_g_c = 11;
    try std.testing.expectError(error.CanopySymbioticSenescenceWouldOverdraw, calculateAll(&branches, &outputs, testInputs()));
    try std.testing.expectEqualDeep(before, outputs);
}
