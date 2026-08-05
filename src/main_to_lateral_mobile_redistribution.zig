const std = @import("std");

pub const ElementMass = struct {
    carbon_g_c: f64 = 0,
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
};

pub const State = struct {
    time_since_germination_h_by_branch: []f64,
    mobile_mass_by_branch: []ElementMass,
};

pub const Inputs = struct {
    main_branch: usize,
    remobilization_duration_h: f64,
    redistribution_fraction_per_biological_step: f64,
    canopy_growth_temperature_response: f64,
    canopy_growth_water_fraction: f64,
    biological_timestep_h: f64,
};

/// Exact GROSUB lines 4881--4899 main-to-lateral mobile-pool redistribution.
/// Branches are swept in ascending runtime NB order. Each admitted lateral
/// branch observes the main pool after earlier transfers. Time is h and masses
/// are g C, g N, and g P. A runtime redistribution fraction of 0.05 per
/// biological invocation reproduces the source literal without making it a
/// compile-time parameter. GROSUB deliberately does not multiply this transfer
/// by XNFH; XNFH affects only the lateral ATRP increment here.
pub fn apply(state: State, inputs: Inputs) !usize {
    try validateDimensionsAndInputs(state, inputs);
    var admitted_count: usize = 0;
    for (state.time_since_germination_h_by_branch, 0..) |time_h, branch| {
        if (!std.math.isFinite(time_h) or time_h < 0)
            return error.InvalidMainToLateralRedistributionState;
        if (branch != inputs.main_branch and time_h <= inputs.remobilization_duration_h)
            admitted_count += 1;
    }
    if (admitted_count == 0) return 0;

    try validateMass(state.mobile_mass_by_branch[inputs.main_branch]);
    var prospective_main = state.mobile_mass_by_branch[inputs.main_branch];
    const time_increment_h = inputs.canopy_growth_temperature_response *
        inputs.canopy_growth_water_fraction * inputs.biological_timestep_h;
    if (!std.math.isFinite(time_increment_h) or time_increment_h < 0)
        return error.NonFiniteMainToLateralRedistributionResult;
    for (state.mobile_mass_by_branch, state.time_since_germination_h_by_branch, 0..) |lateral, time_h, branch| {
        if (branch == inputs.main_branch or time_h > inputs.remobilization_duration_h)
            continue;
        try validateMass(lateral);
        const next_time_h = time_h + time_increment_h;
        const transfer = transfers(
            prospective_main,
            lateral,
            inputs.redistribution_fraction_per_biological_step,
            inputs.canopy_growth_temperature_response,
        );
        const next_lateral = add(lateral, transfer);
        prospective_main = subtract(prospective_main, transfer);
        if (!std.math.isFinite(next_time_h))
            return error.NonFiniteMainToLateralRedistributionResult;
        try validateResultMass(next_lateral);
        try validateResultMass(prospective_main);
    }

    for (state.mobile_mass_by_branch, state.time_since_germination_h_by_branch, 0..) |*lateral, *time_h, branch| {
        if (branch == inputs.main_branch or time_h.* > inputs.remobilization_duration_h)
            continue;
        time_h.* = time_h.* + inputs.canopy_growth_temperature_response *
            inputs.canopy_growth_water_fraction * inputs.biological_timestep_h;
        const transfer = transfers(
            state.mobile_mass_by_branch[inputs.main_branch],
            lateral.*,
            inputs.redistribution_fraction_per_biological_step,
            inputs.canopy_growth_temperature_response,
        );
        lateral.carbon_g_c = lateral.carbon_g_c + transfer.carbon_g_c;
        lateral.nitrogen_g_n = lateral.nitrogen_g_n + transfer.nitrogen_g_n;
        lateral.phosphorus_g_p = lateral.phosphorus_g_p + transfer.phosphorus_g_p;
        state.mobile_mass_by_branch[inputs.main_branch].carbon_g_c =
            state.mobile_mass_by_branch[inputs.main_branch].carbon_g_c - transfer.carbon_g_c;
        state.mobile_mass_by_branch[inputs.main_branch].nitrogen_g_n =
            state.mobile_mass_by_branch[inputs.main_branch].nitrogen_g_n - transfer.nitrogen_g_n;
        state.mobile_mass_by_branch[inputs.main_branch].phosphorus_g_p =
            state.mobile_mass_by_branch[inputs.main_branch].phosphorus_g_p - transfer.phosphorus_g_p;
    }
    return admitted_count;
}

fn transfers(
    main: ElementMass,
    lateral: ElementMass,
    redistribution_fraction_per_biological_step: f64,
    temperature_response: f64,
) ElementMass {
    return .{
        .carbon_g_c = @max(0.0, redistribution_fraction_per_biological_step * temperature_response *
            (0.5 * (main.carbon_g_c + lateral.carbon_g_c) - lateral.carbon_g_c)),
        .nitrogen_g_n = @max(0.0, redistribution_fraction_per_biological_step * temperature_response *
            (0.5 * (main.nitrogen_g_n + lateral.nitrogen_g_n) - lateral.nitrogen_g_n)),
        .phosphorus_g_p = @max(0.0, redistribution_fraction_per_biological_step * temperature_response *
            (0.5 * (main.phosphorus_g_p + lateral.phosphorus_g_p) - lateral.phosphorus_g_p)),
    };
}

fn add(left: ElementMass, right: ElementMass) ElementMass {
    return .{
        .carbon_g_c = left.carbon_g_c + right.carbon_g_c,
        .nitrogen_g_n = left.nitrogen_g_n + right.nitrogen_g_n,
        .phosphorus_g_p = left.phosphorus_g_p + right.phosphorus_g_p,
    };
}

fn subtract(left: ElementMass, right: ElementMass) ElementMass {
    return .{
        .carbon_g_c = left.carbon_g_c - right.carbon_g_c,
        .nitrogen_g_n = left.nitrogen_g_n - right.nitrogen_g_n,
        .phosphorus_g_p = left.phosphorus_g_p - right.phosphorus_g_p,
    };
}

fn validateDimensionsAndInputs(state: State, inputs: Inputs) !void {
    if (state.mobile_mass_by_branch.len == 0 or
        state.time_since_germination_h_by_branch.len != state.mobile_mass_by_branch.len or
        inputs.main_branch >= state.mobile_mass_by_branch.len)
        return error.MainToLateralRedistributionDimensionMismatch;
    inline for (.{
        inputs.remobilization_duration_h,
        inputs.redistribution_fraction_per_biological_step,
        inputs.canopy_growth_temperature_response,
        inputs.canopy_growth_water_fraction,
        inputs.biological_timestep_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidMainToLateralRedistributionInput;
    if (inputs.remobilization_duration_h == 0 or inputs.biological_timestep_h == 0)
        return error.InvalidMainToLateralRedistributionInput;
}

fn validateMass(mass: ElementMass) !void {
    inline for (.{ mass.carbon_g_c, mass.nitrogen_g_n, mass.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidMainToLateralRedistributionState;
}

fn validateResultMass(mass: ElementMass) !void {
    inline for (.{ mass.carbon_g_c, mass.nitrogen_g_n, mass.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.MainToLateralRedistributionOverdraw;
}

fn testInputs(main_branch: usize) Inputs {
    return .{
        .main_branch = main_branch,
        .remobilization_duration_h = 10,
        .redistribution_fraction_per_biological_step = 0.05,
        .canopy_growth_temperature_response = 1,
        .canopy_growth_water_fraction = 0.5,
        .biological_timestep_h = 2,
    };
}

fn sum(masses: []const ElementMass) ElementMass {
    var result: ElementMass = .{};
    for (masses) |mass| result = add(result, mass);
    return result;
}

test "GROSUB ascending lateral sweep conserves mobile C N P" {
    var times = [_]f64{ 0, 0, 0, 0 };
    var masses = [_]ElementMass{
        .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 },
        .{ .carbon_g_c = 100, .nitrogen_g_n = 10, .phosphorus_g_p = 1 },
        .{ .carbon_g_c = 20, .nitrogen_g_n = 2, .phosphorus_g_p = 0.2 },
        .{ .carbon_g_c = 40, .nitrogen_g_n = 4, .phosphorus_g_p = 0.4 },
    };
    const before = sum(&masses);
    try std.testing.expectEqual(@as(usize, 3), try apply(.{
        .time_since_germination_h_by_branch = &times,
        .mobile_mass_by_branch = &masses,
    }, testInputs(1)));
    const after = sum(&masses);
    try std.testing.expectApproxEqAbs(before.carbon_g_c, after.carbon_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(before.nitrogen_g_n, after.nitrogen_g_n, 1e-14);
    try std.testing.expectApproxEqAbs(before.phosphorus_g_p, after.phosphorus_g_p, 1e-14);
    try std.testing.expectEqual(@as(f64, 0), times[1]);
    try std.testing.expectEqual(@as(f64, 1), times[0]);
    try std.testing.expectEqual(@as(f64, 1), times[2]);
    try std.testing.expectEqual(@as(f64, 1), times[3]);
}

test "later branch observes main depletion from earlier branch" {
    var times = [_]f64{ 0, 0, 0 };
    var masses = [_]ElementMass{
        .{ .carbon_g_c = 0 },
        .{ .carbon_g_c = 100 },
        .{ .carbon_g_c = 0 },
    };
    _ = try apply(.{
        .time_since_germination_h_by_branch = &times,
        .mobile_mass_by_branch = &masses,
    }, testInputs(1));
    // NB=0 takes 2.5 first; NB=2 then takes 0.025 * 97.5.
    try std.testing.expectApproxEqAbs(@as(f64, 2.5), masses[0].carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.4375), masses[2].carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 95.0625), masses[1].carbon_g_c, 1e-15);
}

test "runtime redistribution fraction replaces the source literal" {
    var times = [_]f64{ 0, 0 };
    var masses = [_]ElementMass{
        .{ .carbon_g_c = 100 },
        .{},
    };
    var inputs = testInputs(0);
    inputs.redistribution_fraction_per_biological_step = 0.10;
    _ = try apply(.{
        .time_since_germination_h_by_branch = &times,
        .mobile_mass_by_branch = &masses,
    }, inputs);
    try std.testing.expectApproxEqAbs(@as(f64, 5), masses[1].carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 95), masses[0].carbon_g_c, 1e-15);
}

test "duration and main gates are strict no-ops" {
    var times = [_]f64{ 11, 7 };
    var masses = [_]ElementMass{
        .{ .carbon_g_c = std.math.nan(f64) },
        .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 },
    };
    try std.testing.expectEqual(@as(usize, 0), try apply(.{
        .time_since_germination_h_by_branch = &times,
        .mobile_mass_by_branch = &masses,
    }, testInputs(1)));
}

test "runtime branch count and explicit nonzero main have no source ceiling" {
    const allocator = std.testing.allocator;
    const count = 45;
    const main_branch = 38;
    const times = try allocator.alloc(f64, count);
    defer allocator.free(times);
    const masses = try allocator.alloc(ElementMass, count);
    defer allocator.free(masses);
    @memset(times, 0);
    @memset(masses, .{});
    masses[main_branch] = .{ .carbon_g_c = 100, .nitrogen_g_n = 10, .phosphorus_g_p = 1 };
    try std.testing.expectEqual(count - 1, try apply(.{
        .time_since_germination_h_by_branch = times,
        .mobile_mass_by_branch = masses,
    }, testInputs(main_branch)));
}

test "invalid late admitted branch and overdraw leave sweep atomic" {
    var times = [_]f64{ 0, 0, 0 };
    var masses = [_]ElementMass{
        .{},
        .{ .carbon_g_c = 100, .nitrogen_g_n = 10, .phosphorus_g_p = 1 },
        .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = std.math.nan(f64) },
    };
    const times_before = times;
    const masses_before = masses;
    try std.testing.expectError(
        error.InvalidMainToLateralRedistributionState,
        apply(.{
            .time_since_germination_h_by_branch = &times,
            .mobile_mass_by_branch = &masses,
        }, testInputs(1)),
    );
    try std.testing.expectEqualSlices(f64, &times_before, &times);
    try std.testing.expectEqualSlices(u8, std.mem.asBytes(&masses_before), std.mem.asBytes(&masses));

    masses = .{ .{}, .{ .carbon_g_c = 1, .nitrogen_g_n = 1, .phosphorus_g_p = 1 }, .{} };
    var inputs = testInputs(1);
    inputs.canopy_growth_temperature_response = 100;
    try std.testing.expectError(
        error.MainToLateralRedistributionOverdraw,
        apply(.{
            .time_since_germination_h_by_branch = &times,
            .mobile_mass_by_branch = &masses,
        }, inputs),
    );
}
