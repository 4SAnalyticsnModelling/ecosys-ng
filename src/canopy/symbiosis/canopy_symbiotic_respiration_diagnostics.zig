const std = @import("std");

pub const BranchInputs = struct {
    maintenance_respiration_g_c: f64,
    oxygen_unconstrained_respiration_g_c: f64,
    growth_and_fixation_respiration_g_c: f64,
    recycled_senescence_carbon_g_c: f64,
};

pub const Diagnostics = struct {
    total_plant_carbon_exchange_g_c: f64,
    aboveground_plant_carbon_exchange_g_c: f64,
    net_plant_carbon_fixation_g_c: f64,
    hourly_net_plant_carbon_fixation_g_c: f64,
    ecosystem_carbon_exchange_g_c: f64,
    autotrophic_carbon_exchange_g_c: f64,
};

pub const Inputs = struct {
    fixation_type: u8,
};

/// Exact grosub.f lines 5561--5567 total canopy nodule respiration and signed
/// diagnostic accumulation in ascending NB order. Branch respiration is a
/// positive g C loss per biological timestep; source diagnostic ledgers use
/// the opposite sign and therefore each subtract the same RCO2T.
pub fn accumulate(
    branches: []const BranchInputs,
    respiration_g_c_by_branch: []f64,
    diagnostics: *Diagnostics,
    inputs: Inputs,
) !void {
    if (branches.len == 0 or branches.len != respiration_g_c_by_branch.len)
        return error.CanopySymbioticRespirationDiagnosticDimensionMismatch;
    if (inputs.fixation_type > 6) return error.InvalidSymbioticFixationType;
    if (inputs.fixation_type < 4) return;
    try validateDiagnostics(diagnostics.*);

    var next = diagnostics.*;
    for (branches) |branch| {
        const respiration_g_c = try calculateRespiration(branch);
        subtractRespiration(&next, respiration_g_c);
        try validateDiagnostics(next);
    }
    for (branches, respiration_g_c_by_branch) |branch, *respiration_g_c|
        respiration_g_c.* = try calculateRespiration(branch);
    diagnostics.* = next;
}

fn calculateRespiration(branch: BranchInputs) !f64 {
    inline for (@typeInfo(BranchInputs).@"struct".fields) |field| {
        const value = @field(branch, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCanopySymbioticRespirationDiagnosticState;
    }
    const result = @min(branch.maintenance_respiration_g_c, branch.oxygen_unconstrained_respiration_g_c) +
        branch.growth_and_fixation_respiration_g_c + branch.recycled_senescence_carbon_g_c;
    if (!std.math.isFinite(result)) return error.NonFiniteCanopySymbioticRespirationDiagnostic;
    return result;
}

fn subtractRespiration(diagnostics: *Diagnostics, respiration_g_c: f64) void {
    diagnostics.total_plant_carbon_exchange_g_c -= respiration_g_c;
    diagnostics.aboveground_plant_carbon_exchange_g_c -= respiration_g_c;
    diagnostics.net_plant_carbon_fixation_g_c -= respiration_g_c;
    diagnostics.hourly_net_plant_carbon_fixation_g_c -= respiration_g_c;
    diagnostics.ecosystem_carbon_exchange_g_c -= respiration_g_c;
    diagnostics.autotrophic_carbon_exchange_g_c -= respiration_g_c;
}

fn validateDiagnostics(diagnostics: Diagnostics) !void {
    inline for (@typeInfo(Diagnostics).@"struct".fields) |field|
        if (!std.math.isFinite(@field(diagnostics, field.name)))
            return error.NonFiniteCanopySymbioticRespirationDiagnostic;
}

fn testBranch(scale: f64) BranchInputs {
    return .{
        .maintenance_respiration_g_c = 2 * scale,
        .oxygen_unconstrained_respiration_g_c = 1.5 * scale,
        .growth_and_fixation_respiration_g_c = 0.5 * scale,
        .recycled_senescence_carbon_g_c = 0.25 * scale,
    };
}

fn testDiagnostics() Diagnostics {
    return .{
        .total_plant_carbon_exchange_g_c = 10,
        .aboveground_plant_carbon_exchange_g_c = 9,
        .net_plant_carbon_fixation_g_c = 8,
        .hourly_net_plant_carbon_fixation_g_c = 7,
        .ecosystem_carbon_exchange_g_c = 6,
        .autotrophic_carbon_exchange_g_c = 5,
    };
}

test "GROSUB RCO2T sums constrained maintenance growth and recycled senescence carbon" {
    var respiration: [1]f64 = undefined;
    var diagnostics = testDiagnostics();
    const before = diagnostics;
    try accumulate(&.{testBranch(1)}, &respiration, &diagnostics, .{ .fixation_type = 4 });
    try std.testing.expectEqual(@as(f64, 2.25), respiration[0]);
    inline for (@typeInfo(Diagnostics).@"struct".fields) |field|
        try std.testing.expectEqual(@field(before, field.name) - respiration[0], @field(diagnostics, field.name));
}

test "GROSUB all six respiration ledgers preserve the same negative sign" {
    var respiration: [2]f64 = undefined;
    var diagnostics = testDiagnostics();
    const before = diagnostics;
    try accumulate(&.{ testBranch(1), testBranch(2) }, &respiration, &diagnostics, .{ .fixation_type = 6 });
    const loss = respiration[0] + respiration[1];
    inline for (@typeInfo(Diagnostics).@"struct".fields) |field|
        try std.testing.expectApproxEqAbs(@field(before, field.name) - loss, @field(diagnostics, field.name), 1e-15);
}

test "GROSUB root fixation gate leaves canopy respiration outputs unread and unchanged" {
    var respiration = [_]f64{31};
    var diagnostics: Diagnostics = undefined;
    diagnostics.total_plant_carbon_exchange_g_c = std.math.nan(f64);
    var branch = testBranch(1);
    branch.maintenance_respiration_g_c = std.math.nan(f64);
    try accumulate(&.{branch}, &respiration, &diagnostics, .{ .fixation_type = 3 });
    try std.testing.expectEqual(@as(f64, 31), respiration[0]);
    try std.testing.expect(std.math.isNan(diagnostics.total_plant_carbon_exchange_g_c));
}

test "GROSUB respiration diagnostic sweep preserves runtime order and atomicity" {
    var branches = [_]BranchInputs{testBranch(1)} ** 47;
    var respiration: [47]f64 = undefined;
    var diagnostics = testDiagnostics();
    try accumulate(&branches, &respiration, &diagnostics, .{ .fixation_type = 4 });
    const before_respiration = respiration;
    const before_diagnostics = diagnostics;
    branches[46].growth_and_fixation_respiration_g_c = std.math.nan(f64);
    try std.testing.expectError(error.InvalidCanopySymbioticRespirationDiagnosticState, accumulate(&branches, &respiration, &diagnostics, .{ .fixation_type = 4 }));
    try std.testing.expectEqualDeep(before_respiration, respiration);
    try std.testing.expectEqualDeep(before_diagnostics, diagnostics);
}
