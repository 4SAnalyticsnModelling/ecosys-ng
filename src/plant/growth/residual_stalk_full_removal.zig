const std = @import("std");

pub const Elements = struct { carbon: f64, nitrogen: f64, phosphorus: f64 };
pub const Fractions = struct { carbon: []const f64, nitrogen: []const f64, phosphorus: []const f64 };
pub const FireStatus = enum { absent, active };
pub const Litter = struct {
    woody_carbon_g_c: []f64,
    woody_nitrogen_g_n: []f64,
    woody_phosphorus_g_p: []f64,
    stalk_carbon_g_c: []f64,
    stalk_nitrogen_g_n: []f64,
    stalk_phosphorus_g_p: []f64,
};
pub const State = struct {
    branch_stalk_carbon_g_c: *f64,
    branch_stalk_nitrogen_g_n: *f64,
    branch_stalk_phosphorus_g_p: *f64,
    residual_stalk_carbon_g_c: *f64,
    residual_stalk_nitrogen_g_n: *f64,
    residual_stalk_phosphorus_g_p: *f64,
    litter: Litter,
};
pub const Inputs = struct {
    /// True means the GROSUB 3285 perennial-stalk gate was entered, so its
    /// residual-remobilization branch owns state and this fallback does nothing.
    stalk_senescence_active: bool,
    fire_status: FireStatus,
    woody_fraction: Elements,
    nonwoody_fraction: Elements,
    woody_kinetics: Fractions,
    stalk_kinetics: Fractions,
};
pub const Result = struct {
    removed_from_residual_g: Elements,
    routed_to_litter_g: Elements,
    suppressed_by_fire_g: Elements,
};

/// grosub.f lines 3577--3608. When the earlier stalk-senescence gate is false,
/// publishes all residual stalk to woody then stalk litter C/N/P if no fire,
/// subtracts it from branch totals, and clears residual pools. Fire suppresses
/// litter only; state clearing still occurs. All destinations validate first.
pub fn apply(state: State, inputs: Inputs) !Result {
    if (inputs.stalk_senescence_active) return .{
        .removed_from_residual_g = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 },
        .routed_to_litter_g = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 },
        .suppressed_by_fire_g = .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 },
    };
    const count = try validateTopology(state, inputs);
    const residual: Elements = .{ .carbon = state.residual_stalk_carbon_g_c.*, .nitrogen = state.residual_stalk_nitrogen_g_n.*, .phosphorus = state.residual_stalk_phosphorus_g_p.* };
    const branch_after: Elements = .{
        .carbon = @max(0, state.branch_stalk_carbon_g_c.* - residual.carbon),
        .nitrogen = @max(0, state.branch_stalk_nitrogen_g_n.* - residual.nitrogen),
        .phosphorus = @max(0, state.branch_stalk_phosphorus_g_p.* - residual.phosphorus),
    };
    inline for (.{ residual, branch_after }) |values| inline for (@typeInfo(Elements).@"struct".fields) |field| {
        const value = @field(values, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidResidualStalkFullRemovalState;
    };
    if (inputs.fire_status == .absent) for (0..count) |kinetic| try validateLitter(state, inputs, residual, kinetic);

    if (inputs.fire_status == .absent) for (0..count) |kinetic| commitLitter(state, inputs, residual, kinetic);
    state.branch_stalk_carbon_g_c.* = branch_after.carbon;
    state.branch_stalk_nitrogen_g_n.* = branch_after.nitrogen;
    state.branch_stalk_phosphorus_g_p.* = branch_after.phosphorus;
    state.residual_stalk_carbon_g_c.* = 0;
    state.residual_stalk_nitrogen_g_n.* = 0;
    state.residual_stalk_phosphorus_g_p.* = 0;
    return .{
        .removed_from_residual_g = residual,
        .routed_to_litter_g = if (inputs.fire_status == .absent) residual else .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 },
        .suppressed_by_fire_g = if (inputs.fire_status == .active) residual else .{ .carbon = 0, .nitrogen = 0, .phosphorus = 0 },
    };
}

fn validateTopology(state: State, inputs: Inputs) !usize {
    inline for (.{ inputs.woody_fraction, inputs.nonwoody_fraction }) |values| inline for (@typeInfo(Elements).@"struct".fields) |field| {
        const value = @field(values, field.name);
        if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidResidualStalkFullRemovalInput;
    };
    inline for (@typeInfo(Elements).@"struct".fields) |field| if (@abs(@field(inputs.woody_fraction, field.name) + @field(inputs.nonwoody_fraction, field.name) - 1) > 1e-12) return error.InvalidResidualStalkFullRemovalInput;
    const count = inputs.woody_kinetics.carbon.len;
    if (count == 0) return error.ZeroResidualStalkKineticPools;
    inline for (.{ inputs.woody_kinetics.nitrogen, inputs.woody_kinetics.phosphorus, inputs.stalk_kinetics.carbon, inputs.stalk_kinetics.nitrogen, inputs.stalk_kinetics.phosphorus }) |values| if (values.len != count) return error.ResidualStalkKineticDimensionMismatch;
    inline for (@typeInfo(Litter).@"struct".fields) |field| if (@field(state.litter, field.name).len != count) return error.ResidualStalkKineticDimensionMismatch;
    inline for (.{ inputs.woody_kinetics, inputs.stalk_kinetics }) |kinetics| inline for (@typeInfo(Fractions).@"struct".fields) |field| {
        var total: f64 = 0;
        for (@field(kinetics, field.name)) |value| {
            if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidResidualStalkKineticFraction;
            total += value;
        }
        if (!std.math.isFinite(total) or @abs(total - 1) > 1e-12) return error.InvalidResidualStalkKineticFraction;
    };
    return count;
}
fn additions(inputs: Inputs, residual: Elements, k: usize) [6]f64 {
    return .{
        inputs.woody_kinetics.carbon[k] * residual.carbon * inputs.woody_fraction.carbon,
        inputs.woody_kinetics.nitrogen[k] * residual.nitrogen * inputs.woody_fraction.nitrogen,
        inputs.woody_kinetics.phosphorus[k] * residual.phosphorus * inputs.woody_fraction.phosphorus,
        inputs.stalk_kinetics.carbon[k] * residual.carbon * inputs.nonwoody_fraction.carbon,
        inputs.stalk_kinetics.nitrogen[k] * residual.nitrogen * inputs.nonwoody_fraction.nitrogen,
        inputs.stalk_kinetics.phosphorus[k] * residual.phosphorus * inputs.nonwoody_fraction.phosphorus,
    };
}
fn validateLitter(state: State, inputs: Inputs, residual: Elements, k: usize) !void {
    const add = additions(inputs, residual, k);
    inline for (@typeInfo(Litter).@"struct".fields, 0..) |field, i| {
        const current = @field(state.litter, field.name)[k];
        if (!std.math.isFinite(current) or current < 0 or !std.math.isFinite(add[i]) or add[i] < 0 or !std.math.isFinite(current + add[i])) return error.InvalidResidualStalkLitterResult;
    }
}
fn commitLitter(state: State, inputs: Inputs, residual: Elements, k: usize) void {
    const add = additions(inputs, residual, k);
    inline for (@typeInfo(Litter).@"struct".fields, 0..) |field, i| @field(state.litter, field.name)[k] += add[i];
}
fn fixture(s: *[6]f64, l: *[6][2]f64) State {
    return .{ .branch_stalk_carbon_g_c = &s[0], .branch_stalk_nitrogen_g_n = &s[1], .branch_stalk_phosphorus_g_p = &s[2], .residual_stalk_carbon_g_c = &s[3], .residual_stalk_nitrogen_g_n = &s[4], .residual_stalk_phosphorus_g_p = &s[5], .litter = .{ .woody_carbon_g_c = &l[0], .woody_nitrogen_g_n = &l[1], .woody_phosphorus_g_p = &l[2], .stalk_carbon_g_c = &l[3], .stalk_nitrogen_g_n = &l[4], .stalk_phosphorus_g_p = &l[5] } };
}
fn testInputs(fire: FireStatus) Inputs {
    const f = &[_]f64{ 0.25, 0.75 };
    return .{ .stalk_senescence_active = false, .fire_status = fire, .woody_fraction = .{ .carbon = 0.2, .nitrogen = 0.3, .phosphorus = 0.4 }, .nonwoody_fraction = .{ .carbon = 0.8, .nitrogen = 0.7, .phosphorus = 0.6 }, .woody_kinetics = .{ .carbon = f, .nitrogen = f, .phosphorus = f }, .stalk_kinetics = .{ .carbon = f, .nitrogen = f, .phosphorus = f } };
}

test "no-fire fallback conserves C N P through runtime litter pools" {
    var s = [_]f64{ 10, 2, 0.5, 8, 1, 0.2 };
    var l: [6][2]f64 = @splat(@splat(0));
    const result = try apply(fixture(&s, &l), testInputs(.absent));
    try std.testing.expectEqual(@as(f64, 2), s[0]);
    try std.testing.expectEqual(@as(f64, 0), s[3]);
    inline for (@typeInfo(Elements).@"struct".fields, 0..) |field, i| {
        var total: f64 = 0;
        total += l[i][0] + l[i][1] + l[i + 3][0] + l[i + 3][1];
        try std.testing.expectApproxEqAbs(@field(result.removed_from_residual_g, field.name), total, 1e-15);
    }
}
test "fire suppresses litter but still clears residual and branch totals" {
    var s = [_]f64{ 10, 2, 0.5, 8, 1, 0.2 };
    var l: [6][2]f64 = @splat(@splat(0));
    const result = try apply(fixture(&s, &l), testInputs(.active));
    try std.testing.expectEqual(@as(f64, 2), s[0]);
    try std.testing.expectEqual(@as(f64, 0), s[3]);
    try std.testing.expectEqual(@as(f64, 8), result.suppressed_by_fire_g.carbon);
    try std.testing.expectEqual(@as(f64, 0), l[0][0]);
}
test "active stalk branch is a strict no-op before topology validation" {
    var s = [_]f64{ 10, 2, 0.5, 8, 1, 0.2 };
    var l: [6][2]f64 = @splat(@splat(std.math.nan(f64)));
    var inputs = testInputs(.absent);
    inputs.stalk_senescence_active = true;
    _ = try apply(fixture(&s, &l), inputs);
    try std.testing.expectEqual(@as(f64, 8), s[3]);
}
test "late invalid litter is atomic" {
    var s = [_]f64{ 10, 2, 0.5, 8, 1, 0.2 };
    var l: [6][2]f64 = @splat(@splat(0));
    l[5][1] = std.math.nan(f64);
    try std.testing.expectError(error.InvalidResidualStalkLitterResult, apply(fixture(&s, &l), testInputs(.absent)));
    try std.testing.expectEqual(@as(f64, 8), s[3]);
    try std.testing.expectEqual(@as(f64, 0), l[0][0]);
}
