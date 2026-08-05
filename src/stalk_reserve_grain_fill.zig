const std = @import("std");

pub const ElementMass = struct {
    carbon_g_c: f64 = 0,
    nitrogen_g_n: f64 = 0,
    phosphorus_g_p: f64 = 0,
};

pub const State = struct {
    stalk_reserve_by_branch: []ElementMass,
    grain_by_branch: []ElementMass,
};

pub const Inputs = struct {
    grain_fill_started_by_branch: []const bool,
    seed_count_by_branch: []const f64,
    maximum_individual_seed_carbon_g_c_by_branch: []const f64,
    grain_precursor_growth_by_branch: []const ElementMass,
    root_profile_type: u8,
    grain_fill_carbon_g_c_per_seed_h_at_25c: f64,
    canopy_growth_temperature_response: f64,
    planting_layer_root_growth_temperature_response: f64,
    biological_timestep_h: f64,
    minimum_grain_nutrient_fraction: f64,
    maximum_grain_nitrogen_per_carbon_g_n_per_g_c: f64,
    maximum_grain_phosphorus_per_carbon_g_p_per_g_c: f64,
    reserve_nitrogen_half_saturation_g_n_per_g_c: f64,
    reserve_phosphorus_half_saturation_g_p_per_g_c: f64,
    carbon_presence_threshold_g_c: f64,
};

pub const Result = struct {
    carbon_translocated_g_c: f64 = 0,
    nitrogen_translocated_g_n: f64 = 0,
    phosphorus_translocated_g_p: f64 = 0,
    active_branch_count: usize = 0,
};

/// Exact GROSUB lines 5126--5228 grain fill from stalk reserves. Runtime
/// branches are swept in ascending NB order. The source three-argument AMIN1
/// constraints and sequential N:P coupling are retained, including signed
/// nutrient transfers. Units are g C, g N, g P, and h.
pub fn apply(state: State, inputs: Inputs) !Result {
    try validateDimensionsAndParameters(state, inputs);
    var result: Result = .{};
    for (0..state.stalk_reserve_by_branch.len) |branch| {
        if (!inputs.grain_fill_started_by_branch[branch]) continue;
        const branch_result = try calculateBranch(
            state.stalk_reserve_by_branch[branch],
            state.grain_by_branch[branch],
            inputs.seed_count_by_branch[branch],
            inputs.maximum_individual_seed_carbon_g_c_by_branch[branch],
            inputs.grain_precursor_growth_by_branch[branch],
            inputs,
        );
        result.carbon_translocated_g_c += branch_result.carbon_translocated_g_c;
        result.nitrogen_translocated_g_n += branch_result.nitrogen_translocated_g_n;
        result.phosphorus_translocated_g_p += branch_result.phosphorus_translocated_g_p;
        result.active_branch_count += 1;
        inline for (.{
            result.carbon_translocated_g_c,
            result.nitrogen_translocated_g_n,
            result.phosphorus_translocated_g_p,
        }) |value| if (!std.math.isFinite(value))
            return error.NonFiniteGrainFillSweepResult;
    }

    for (0..state.stalk_reserve_by_branch.len) |branch| {
        if (!inputs.grain_fill_started_by_branch[branch]) continue;
        const branch_result = try calculateBranch(
            state.stalk_reserve_by_branch[branch],
            state.grain_by_branch[branch],
            inputs.seed_count_by_branch[branch],
            inputs.maximum_individual_seed_carbon_g_c_by_branch[branch],
            inputs.grain_precursor_growth_by_branch[branch],
            inputs,
        );
        state.stalk_reserve_by_branch[branch].carbon_g_c = branch_result.reserve_after.carbon_g_c;
        state.stalk_reserve_by_branch[branch].nitrogen_g_n = branch_result.reserve_after.nitrogen_g_n;
        state.stalk_reserve_by_branch[branch].phosphorus_g_p = branch_result.reserve_after.phosphorus_g_p;
        state.grain_by_branch[branch].carbon_g_c = branch_result.grain_after.carbon_g_c;
        state.grain_by_branch[branch].nitrogen_g_n = branch_result.grain_after.nitrogen_g_n;
        state.grain_by_branch[branch].phosphorus_g_p = branch_result.grain_after.phosphorus_g_p;
    }
    return result;
}

const BranchResult = struct {
    reserve_after: ElementMass,
    grain_after: ElementMass,
    carbon_translocated_g_c: f64,
    nitrogen_translocated_g_n: f64,
    phosphorus_translocated_g_p: f64,
};

fn calculateBranch(
    reserve: ElementMass,
    grain: ElementMass,
    seed_count: f64,
    maximum_individual_seed_carbon_g_c: f64,
    precursor_growth: ElementMass,
    inputs: Inputs,
) !BranchResult {
    try validateMass(reserve);
    try validateMass(grain);
    try validateMass(precursor_growth);
    inline for (.{ seed_count, maximum_individual_seed_carbon_g_c }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidGrainFillBranchState;

    const sink_carbon_g_c = maximum_individual_seed_carbon_g_c * seed_count;
    const temperature_response = if (inputs.root_profile_type == 0)
        inputs.canopy_growth_temperature_response
    else
        inputs.planting_layer_root_growth_temperature_response;
    const maximum_grain_fill_g_c = if (grain.carbon_g_c >= sink_carbon_g_c)
        0
    else
        @max(0.0, inputs.grain_fill_carbon_g_c_per_seed_h_at_25c * seed_count *
            @sqrt(temperature_response) * inputs.biological_timestep_h);
    const nutrient_deficient = grain.nitrogen_g_n <
        inputs.minimum_grain_nutrient_fraction *
            inputs.maximum_grain_nitrogen_per_carbon_g_n_per_g_c * grain.carbon_g_c or
        grain.phosphorus_g_p < inputs.minimum_grain_nutrient_fraction *
            inputs.maximum_grain_phosphorus_per_carbon_g_p_per_g_c * grain.carbon_g_c;
    const actual_grain_fill_g_c = if (nutrient_deficient) 0 else maximum_grain_fill_g_c;
    const maximum_carbon_translocation_g_c = @min(maximum_grain_fill_g_c, reserve.carbon_g_c);
    const carbon_translocation_g_c = @min(actual_grain_fill_g_c, reserve.carbon_g_c);
    const responsive_fraction = 1.0 - inputs.minimum_grain_nutrient_fraction;

    var nitrogen_translocation_g_n: f64 = 0;
    if (reserve.nitrogen_g_n > inputs.carbon_presence_threshold_g_c) {
        const limitation = reserve.nitrogen_g_n /
            (reserve.nitrogen_g_n + inputs.reserve_nitrogen_half_saturation_g_n_per_g_c *
                reserve.carbon_g_c);
        const grain_ratio = inputs.minimum_grain_nutrient_fraction + responsive_fraction *
            @max(0.0, @min(1.0, limitation));
        nitrogen_translocation_g_n = @min(
            maximum_carbon_translocation_g_c *
                inputs.maximum_grain_nitrogen_per_carbon_g_n_per_g_c,
            @min(
                @max(0.0, reserve.nitrogen_g_n * grain_ratio),
                (grain.carbon_g_c + carbon_translocation_g_c) *
                    inputs.maximum_grain_nitrogen_per_carbon_g_n_per_g_c -
                    grain.nitrogen_g_n,
            ),
        );
    }
    var phosphorus_translocation_g_p: f64 = 0;
    if (reserve.phosphorus_g_p > inputs.carbon_presence_threshold_g_c) {
        const limitation = reserve.phosphorus_g_p /
            (reserve.phosphorus_g_p +
                inputs.reserve_phosphorus_half_saturation_g_p_per_g_c * reserve.carbon_g_c);
        const grain_ratio = inputs.minimum_grain_nutrient_fraction + responsive_fraction *
            @max(0.0, @min(1.0, limitation));
        phosphorus_translocation_g_p = @min(
            maximum_carbon_translocation_g_c *
                inputs.maximum_grain_phosphorus_per_carbon_g_p_per_g_c,
            @min(
                @max(0.0, reserve.phosphorus_g_p * grain_ratio),
                (grain.carbon_g_c + carbon_translocation_g_c) *
                    inputs.maximum_grain_phosphorus_per_carbon_g_p_per_g_c -
                    grain.phosphorus_g_p,
            ),
        );
    }
    nitrogen_translocation_g_n = @min(
        nitrogen_translocation_g_n,
        phosphorus_translocation_g_p *
            inputs.maximum_grain_nitrogen_per_carbon_g_n_per_g_c /
            inputs.maximum_grain_phosphorus_per_carbon_g_p_per_g_c,
    );
    phosphorus_translocation_g_p = @min(
        phosphorus_translocation_g_p,
        nitrogen_translocation_g_n *
            inputs.maximum_grain_phosphorus_per_carbon_g_p_per_g_c /
            inputs.maximum_grain_nitrogen_per_carbon_g_n_per_g_c,
    );
    const reserve_after: ElementMass = .{
        .carbon_g_c = reserve.carbon_g_c + precursor_growth.carbon_g_c -
            carbon_translocation_g_c,
        .nitrogen_g_n = reserve.nitrogen_g_n + precursor_growth.nitrogen_g_n -
            nitrogen_translocation_g_n,
        .phosphorus_g_p = reserve.phosphorus_g_p + precursor_growth.phosphorus_g_p -
            phosphorus_translocation_g_p,
    };
    const grain_after: ElementMass = .{
        .carbon_g_c = grain.carbon_g_c + carbon_translocation_g_c,
        .nitrogen_g_n = grain.nitrogen_g_n + nitrogen_translocation_g_n,
        .phosphorus_g_p = grain.phosphorus_g_p + phosphorus_translocation_g_p,
    };
    try validateResultMass(reserve_after);
    try validateResultMass(grain_after);
    return .{
        .reserve_after = reserve_after,
        .grain_after = grain_after,
        .carbon_translocated_g_c = carbon_translocation_g_c,
        .nitrogen_translocated_g_n = nitrogen_translocation_g_n,
        .phosphorus_translocated_g_p = phosphorus_translocation_g_p,
    };
}

fn validateDimensionsAndParameters(state: State, inputs: Inputs) !void {
    const count = state.stalk_reserve_by_branch.len;
    if (count == 0 or state.grain_by_branch.len != count or
        inputs.grain_fill_started_by_branch.len != count or
        inputs.seed_count_by_branch.len != count or
        inputs.maximum_individual_seed_carbon_g_c_by_branch.len != count or
        inputs.grain_precursor_growth_by_branch.len != count)
        return error.GrainFillSweepDimensionMismatch;
    inline for (.{
        inputs.grain_fill_carbon_g_c_per_seed_h_at_25c,
        inputs.canopy_growth_temperature_response,
        inputs.planting_layer_root_growth_temperature_response,
        inputs.biological_timestep_h,
        inputs.minimum_grain_nutrient_fraction,
        inputs.maximum_grain_nitrogen_per_carbon_g_n_per_g_c,
        inputs.maximum_grain_phosphorus_per_carbon_g_p_per_g_c,
        inputs.reserve_nitrogen_half_saturation_g_n_per_g_c,
        inputs.reserve_phosphorus_half_saturation_g_p_per_g_c,
        inputs.carbon_presence_threshold_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidGrainFillSweepInput;
    if (inputs.root_profile_type > 3 or inputs.biological_timestep_h == 0 or
        inputs.minimum_grain_nutrient_fraction > 1 or
        inputs.maximum_grain_nitrogen_per_carbon_g_n_per_g_c == 0 or
        inputs.maximum_grain_phosphorus_per_carbon_g_p_per_g_c == 0)
        return error.InvalidGrainFillSweepInput;
}

fn validateMass(mass: ElementMass) !void {
    inline for (.{ mass.carbon_g_c, mass.nitrogen_g_n, mass.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidGrainFillBranchState;
}

fn validateResultMass(mass: ElementMass) !void {
    inline for (.{ mass.carbon_g_c, mass.nitrogen_g_n, mass.phosphorus_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.GrainFillTranslocationOverdraw;
}

fn testInputs(count: usize) Inputs {
    return .{
        .grain_fill_started_by_branch = if (count == 2) &.{ true, true } else &.{},
        .seed_count_by_branch = if (count == 2) &.{ 10, 20 } else &.{},
        .maximum_individual_seed_carbon_g_c_by_branch = if (count == 2) &.{ 5, 5 } else &.{},
        .grain_precursor_growth_by_branch = if (count == 2) &.{ .{}, .{} } else &.{},
        .root_profile_type = 0,
        .grain_fill_carbon_g_c_per_seed_h_at_25c = 0.1,
        .canopy_growth_temperature_response = 1,
        .planting_layer_root_growth_temperature_response = 0.25,
        .biological_timestep_h = 1,
        .minimum_grain_nutrient_fraction = 0.5,
        .maximum_grain_nitrogen_per_carbon_g_n_per_g_c = 0.1,
        .maximum_grain_phosphorus_per_carbon_g_p_per_g_c = 0.01,
        .reserve_nitrogen_half_saturation_g_n_per_g_c = 0.02,
        .reserve_phosphorus_half_saturation_g_p_per_g_c = 0.002,
        .carbon_presence_threshold_g_c = 1e-12,
    };
}

fn total(reserves: []const ElementMass, grains: []const ElementMass) ElementMass {
    var result: ElementMass = .{};
    for (reserves) |mass| {
        result.carbon_g_c += mass.carbon_g_c;
        result.nitrogen_g_n += mass.nitrogen_g_n;
        result.phosphorus_g_p += mass.phosphorus_g_p;
    }
    for (grains) |mass| {
        result.carbon_g_c += mass.carbon_g_c;
        result.nitrogen_g_n += mass.nitrogen_g_n;
        result.phosphorus_g_p += mass.phosphorus_g_p;
    }
    return result;
}

test "GROSUB grain fill conserves reserve grain and precursor C N P" {
    var reserves = [_]ElementMass{
        .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 },
        .{ .carbon_g_c = 20, .nitrogen_g_n = 2, .phosphorus_g_p = 0.2 },
    };
    var grains = [_]ElementMass{
        .{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01 },
        .{ .carbon_g_c = 2, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.02 },
    };
    var inputs = testInputs(2);
    inputs.grain_precursor_growth_by_branch = &.{
        .{ .carbon_g_c = 0.2, .nitrogen_g_n = 0.02, .phosphorus_g_p = 0.002 },
        .{ .carbon_g_c = 0.4, .nitrogen_g_n = 0.04, .phosphorus_g_p = 0.004 },
    };
    const before = total(&reserves, &grains);
    const result = try apply(.{
        .stalk_reserve_by_branch = &reserves,
        .grain_by_branch = &grains,
    }, inputs);
    const after = total(&reserves, &grains);
    try std.testing.expectEqual(@as(usize, 2), result.active_branch_count);
    try std.testing.expectApproxEqAbs(before.carbon_g_c + 0.6, after.carbon_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(before.nitrogen_g_n + 0.06, after.nitrogen_g_n, 1e-14);
    try std.testing.expectApproxEqAbs(before.phosphorus_g_p + 0.006, after.phosphorus_g_p, 1e-14);
}

test "source three-way nutrient minimum uses grain deficit not maximum" {
    var reserves = [_]ElementMass{.{ .carbon_g_c = 100, .nitrogen_g_n = 10, .phosphorus_g_p = 10 }};
    var grains = [_]ElementMass{.{ .carbon_g_c = 10, .nitrogen_g_n = 1.8, .phosphorus_g_p = 1.8 }};
    var inputs = testInputs(2);
    inputs.grain_fill_started_by_branch = &.{true};
    inputs.seed_count_by_branch = &.{100};
    inputs.maximum_individual_seed_carbon_g_c_by_branch = &.{5};
    inputs.grain_precursor_growth_by_branch = &.{.{}};
    inputs.grain_fill_carbon_g_c_per_seed_h_at_25c = 0.1;
    inputs.maximum_grain_phosphorus_per_carbon_g_p_per_g_c = 0.1;
    const result = try apply(.{
        .stalk_reserve_by_branch = &reserves,
        .grain_by_branch = &grains,
    }, inputs);
    // XLOCC=10, target nutrient=2, so the third AMIN1 argument is 0.2.
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), result.nitrogen_translocated_g_n, 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), result.phosphorus_translocated_g_p, 1e-14);
}

test "shallow and rooted profiles select different temperature responses" {
    var reserves = [_]ElementMass{.{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 }};
    var grains = [_]ElementMass{.{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01 }};
    var inputs = testInputs(2);
    inputs.grain_fill_started_by_branch = &.{true};
    inputs.seed_count_by_branch = &.{10};
    inputs.maximum_individual_seed_carbon_g_c_by_branch = &.{5};
    inputs.grain_precursor_growth_by_branch = &.{.{}};
    const shallow = try apply(.{ .stalk_reserve_by_branch = &reserves, .grain_by_branch = &grains }, inputs);
    reserves = .{.{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 }};
    grains = .{.{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01 }};
    inputs.root_profile_type = 1;
    const rooted = try apply(.{ .stalk_reserve_by_branch = &reserves, .grain_by_branch = &grains }, inputs);
    try std.testing.expectApproxEqAbs(shallow.carbon_translocated_g_c * 0.5, rooted.carbon_translocated_g_c, 1e-15);
}

test "inactive stage does not read branch pools" {
    var reserves = [_]ElementMass{.{ .carbon_g_c = std.math.nan(f64) }};
    var grains = [_]ElementMass{.{ .carbon_g_c = std.math.nan(f64) }};
    var inputs = testInputs(2);
    inputs.grain_fill_started_by_branch = &.{false};
    inputs.seed_count_by_branch = &.{std.math.nan(f64)};
    inputs.maximum_individual_seed_carbon_g_c_by_branch = &.{std.math.nan(f64)};
    inputs.grain_precursor_growth_by_branch = &.{.{ .carbon_g_c = std.math.nan(f64) }};
    const result = try apply(.{
        .stalk_reserve_by_branch = &reserves,
        .grain_by_branch = &grains,
    }, inputs);
    try std.testing.expectEqual(@as(usize, 0), result.active_branch_count);
}

test "runtime branch and seed counts have no source ceiling" {
    const allocator = std.testing.allocator;
    const count = 41;
    const reserves = try allocator.alloc(ElementMass, count);
    defer allocator.free(reserves);
    const grains = try allocator.alloc(ElementMass, count);
    defer allocator.free(grains);
    const flags = try allocator.alloc(bool, count);
    defer allocator.free(flags);
    const seed_counts = try allocator.alloc(f64, count);
    defer allocator.free(seed_counts);
    const seed_sizes = try allocator.alloc(f64, count);
    defer allocator.free(seed_sizes);
    const precursors = try allocator.alloc(ElementMass, count);
    defer allocator.free(precursors);
    @memset(reserves, .{ .carbon_g_c = 100, .nitrogen_g_n = 10, .phosphorus_g_p = 1 });
    @memset(grains, .{});
    @memset(flags, true);
    @memset(seed_counts, 1000);
    @memset(seed_sizes, 1);
    @memset(precursors, .{});
    var inputs = testInputs(2);
    inputs.grain_fill_started_by_branch = flags;
    inputs.seed_count_by_branch = seed_counts;
    inputs.maximum_individual_seed_carbon_g_c_by_branch = seed_sizes;
    inputs.grain_precursor_growth_by_branch = precursors;
    const result = try apply(.{ .stalk_reserve_by_branch = reserves, .grain_by_branch = grains }, inputs);
    try std.testing.expectEqual(count, result.active_branch_count);
}

test "late invalid branch leaves full sweep atomic" {
    var reserves = [_]ElementMass{
        .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.1 },
        .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = std.math.nan(f64) },
    };
    var grains = [_]ElementMass{ .{}, .{} };
    const reserve_before = reserves;
    try std.testing.expectError(
        error.InvalidGrainFillBranchState,
        apply(.{ .stalk_reserve_by_branch = &reserves, .grain_by_branch = &grains }, testInputs(2)),
    );
    try std.testing.expectEqualSlices(u8, std.mem.asBytes(&reserve_before), std.mem.asBytes(&reserves));
}
