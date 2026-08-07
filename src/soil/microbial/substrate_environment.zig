const std = @import("std");

pub const Inputs = struct {
    substrate_complex_count: usize,
    structural_fraction_count: usize,
    residue_fraction_count: usize,
    structural_carbon_g_c: []const f64,
    colonized_structural_carbon_g_c: []const f64,
    microbial_residue_carbon_g_c: []const f64,
    adsorbed_organic_carbon_g_c: []const f64,
    adsorbed_acetate_carbon_g_c: []const f64,
    dissolved_organic_carbon_g_c: []const f64,
    dissolved_acetate_carbon_g_c: []const f64,
    dissolved_organic_nitrogen_g_n: []const f64,
    dissolved_organic_phosphorus_g_p: []const f64,
    biologically_active_water_m3: f64,
    negligible_carbon_g_c: f64,
    negligible_water_m3: f64,
    negligible_complex_fraction: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    substrate_complex_count: usize,
    total_solid_carbon_g_c: []f64,
    total_colonized_structural_carbon_g_c: []f64,
    total_microbial_residue_carbon_g_c: []f64,
    total_carbon_with_residue_and_adsorbates_g_c: []f64,
    colonized_carbon_with_residue_and_adsorbates_g_c: []f64,
    substrate_complex_fraction: []f64,
    dissolved_organic_carbon_concentration_g_c_per_m3: []f64,
    dissolved_acetate_carbon_concentration_g_c_per_m3: []f64,
    dissolved_nitrogen_per_organic_carbon_g_n_per_g_c: []f64,
    dissolved_phosphorus_per_organic_carbon_g_p_per_g_c: []f64,
    dissolved_organic_carbon_fraction: []f64,
    dissolved_acetate_carbon_fraction: []f64,
    all_solid_carbon_g_c: f64,
    all_colonized_structural_carbon_g_c: f64,
    all_microbial_residue_carbon_g_c: f64,
    all_adsorbed_carbon_g_c: f64,
    all_colonized_carbon_with_residue_and_adsorbates_g_c: f64,

    pub fn init(allocator: std.mem.Allocator, substrate_complex_count: usize) !State {
        if (substrate_complex_count == 0) return error.InvalidSubstrateComplexCount;
        var state: State = undefined;
        state.allocator = allocator;
        state.substrate_complex_count = substrate_complex_count;
        state.all_solid_carbon_g_c = 0;
        state.all_colonized_structural_carbon_g_c = 0;
        state.all_microbial_residue_carbon_g_c = 0;
        state.all_adsorbed_carbon_g_c = 0;
        state.all_colonized_carbon_with_residue_and_adsorbates_g_c = 0;
        var allocated: usize = 0;
        errdefer {
            inline for (@typeInfo(State).@"struct".fields) |field| {
                if (field.type == []f64 and allocated > 0) {
                    allocated -= 1;
                    allocator.free(@field(state, field.name));
                }
            }
        }
        inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64) {
                @field(state, field.name) =
                    try allocator.alloc(f64, substrate_complex_count);
                @memset(@field(state, field.name), 0);
                allocated += 1;
            }
        }
        return state;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

/// Ports NITRO.F 349--410 and 481--539. The legacy dimensions were five
/// structural and two microbial-residue fractions; both are runtime axes
/// here. `biologically_active_water_m3` is the VOLWY result from the adjacent
/// mineral-layer microbial environment calculation.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.substrate_complex_count);
    defer staged.deinit();
    calculateValidated(&staged, inputs);
    try validateResult(&staged);
    inline for (@typeInfo(State).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .pointer => if (field.type == []f64)
            @memcpy(@field(state, field.name), @field(staged, field.name)),
        .float => @field(state, field.name) = @field(staged, field.name),
        else => {},
    };
}

fn calculateValidated(state: *State, inputs: Inputs) void {
    clear(state);

    for (0..inputs.substrate_complex_count) |complex| {
        const structural_start = complex * inputs.structural_fraction_count;
        for (inputs.structural_carbon_g_c[structural_start..][0..inputs.structural_fraction_count]) |carbon|
            state.total_solid_carbon_g_c[complex] += carbon;
        for (inputs.colonized_structural_carbon_g_c[structural_start..][0..inputs.structural_fraction_count]) |carbon|
            state.total_colonized_structural_carbon_g_c[complex] += carbon;
        state.all_solid_carbon_g_c += state.total_solid_carbon_g_c[complex];
        state.all_colonized_structural_carbon_g_c +=
            state.total_colonized_structural_carbon_g_c[complex];
    }

    for (0..inputs.substrate_complex_count) |complex| {
        const residue_start = complex * inputs.residue_fraction_count;
        for (inputs.microbial_residue_carbon_g_c[residue_start..][0..inputs.residue_fraction_count]) |carbon|
            state.total_microbial_residue_carbon_g_c[complex] += carbon;
        state.all_microbial_residue_carbon_g_c +=
            state.total_microbial_residue_carbon_g_c[complex];
        state.all_adsorbed_carbon_g_c +=
            inputs.adsorbed_organic_carbon_g_c[complex];
        state.all_adsorbed_carbon_g_c +=
            inputs.adsorbed_acetate_carbon_g_c[complex];
    }

    for (0..inputs.substrate_complex_count) |complex| {
        state.total_carbon_with_residue_and_adsorbates_g_c[complex] =
            state.total_solid_carbon_g_c[complex] +
            state.total_microbial_residue_carbon_g_c[complex] +
            inputs.adsorbed_organic_carbon_g_c[complex] +
            inputs.adsorbed_acetate_carbon_g_c[complex];
        state.colonized_carbon_with_residue_and_adsorbates_g_c[complex] =
            state.total_colonized_structural_carbon_g_c[complex] +
            state.total_microbial_residue_carbon_g_c[complex] +
            inputs.adsorbed_organic_carbon_g_c[complex] +
            inputs.adsorbed_acetate_carbon_g_c[complex];
    }
    state.all_colonized_carbon_with_residue_and_adsorbates_g_c =
        state.all_colonized_structural_carbon_g_c +
        state.all_microbial_residue_carbon_g_c +
        state.all_adsorbed_carbon_g_c;

    for (0..inputs.substrate_complex_count) |complex| {
        state.substrate_complex_fraction[complex] =
            if (state.all_colonized_carbon_with_residue_and_adsorbates_g_c >
            inputs.negligible_carbon_g_c)
                state.colonized_carbon_with_residue_and_adsorbates_g_c[complex] /
                    state.all_colonized_carbon_with_residue_and_adsorbates_g_c
            else
                1;

        if (inputs.biologically_active_water_m3 > inputs.negligible_water_m3) {
            const fraction = state.substrate_complex_fraction[complex];
            const concentration_water_m3 =
                if (fraction > inputs.negligible_complex_fraction)
                    inputs.biologically_active_water_m3 * fraction
                else
                    inputs.biologically_active_water_m3;
            state.dissolved_organic_carbon_concentration_g_c_per_m3[complex] =
                @max(
                    0,
                    inputs.dissolved_organic_carbon_g_c[complex] /
                        concentration_water_m3,
                );
            state.dissolved_acetate_carbon_concentration_g_c_per_m3[complex] =
                @max(
                    0,
                    inputs.dissolved_acetate_carbon_g_c[complex] /
                        concentration_water_m3,
                );
        }

        const doc = inputs.dissolved_organic_carbon_g_c[complex];
        const acetate = inputs.dissolved_acetate_carbon_g_c[complex];
        if (doc > inputs.negligible_carbon_g_c) {
            state.dissolved_nitrogen_per_organic_carbon_g_n_per_g_c[complex] =
                @max(
                    0,
                    inputs.dissolved_organic_nitrogen_g_n[complex] / doc,
                );
            state.dissolved_phosphorus_per_organic_carbon_g_p_per_g_c[complex] =
                @max(
                    0,
                    inputs.dissolved_organic_phosphorus_g_p[complex] / doc,
                );
        }
        if (doc > inputs.negligible_carbon_g_c and
            acetate > inputs.negligible_carbon_g_c)
        {
            state.dissolved_organic_carbon_fraction[complex] =
                doc / (doc + acetate);
            state.dissolved_acetate_carbon_fraction[complex] =
                1 - state.dissolved_organic_carbon_fraction[complex];
        } else if (doc > inputs.negligible_carbon_g_c) {
            state.dissolved_organic_carbon_fraction[complex] = 1;
        } else {
            // Exact source fallback, including the all-zero case.
            state.dissolved_acetate_carbon_fraction[complex] = 1;
        }
    }
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| switch (@typeInfo(field.type)) {
        .pointer => if (field.type == []f64) {
            for (@field(state, field.name)) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteSubstrateEnvironmentResult;
        },
        .float => {
            const value = @field(state, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteSubstrateEnvironmentResult;
        },
        else => {},
    };
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.substrate_complex_count == 0 or
        inputs.structural_fraction_count == 0 or
        inputs.residue_fraction_count == 0 or
        state.substrate_complex_count != inputs.substrate_complex_count)
        return error.InvalidSubstrateEnvironmentDimensions;
    const structural_count = try std.math.mul(
        usize,
        inputs.substrate_complex_count,
        inputs.structural_fraction_count,
    );
    const residue_count = try std.math.mul(
        usize,
        inputs.substrate_complex_count,
        inputs.residue_fraction_count,
    );
    if (inputs.structural_carbon_g_c.len != structural_count or
        inputs.colonized_structural_carbon_g_c.len != structural_count or
        inputs.microbial_residue_carbon_g_c.len != residue_count)
        return error.InvalidSubstrateEnvironmentDimensions;
    inline for (.{
        inputs.adsorbed_organic_carbon_g_c,
        inputs.adsorbed_acetate_carbon_g_c,
        inputs.dissolved_organic_carbon_g_c,
        inputs.dissolved_acetate_carbon_g_c,
        inputs.dissolved_organic_nitrogen_g_n,
        inputs.dissolved_organic_phosphorus_g_p,
    }) |values| if (values.len != inputs.substrate_complex_count)
        return error.InvalidSubstrateEnvironmentDimensions;
    inline for (.{
        inputs.structural_carbon_g_c,
        inputs.colonized_structural_carbon_g_c,
        inputs.microbial_residue_carbon_g_c,
        inputs.adsorbed_organic_carbon_g_c,
        inputs.adsorbed_acetate_carbon_g_c,
        inputs.dissolved_organic_carbon_g_c,
        inputs.dissolved_acetate_carbon_g_c,
        inputs.dissolved_organic_nitrogen_g_n,
        inputs.dissolved_organic_phosphorus_g_p,
    }) |values| for (values) |value| {
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSubstrateEnvironmentInput;
    };
    inline for (.{
        inputs.biologically_active_water_m3,
        inputs.negligible_carbon_g_c,
        inputs.negligible_water_m3,
        inputs.negligible_complex_fraction,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSubstrateEnvironmentInput;
}

fn clear(state: *State) void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) @memset(@field(state, field.name), 0);
    state.all_solid_carbon_g_c = 0;
    state.all_colonized_structural_carbon_g_c = 0;
    state.all_microbial_residue_carbon_g_c = 0;
    state.all_adsorbed_carbon_g_c = 0;
    state.all_colonized_carbon_with_residue_and_adsorbates_g_c = 0;
}

test "runtime substrate aggregation reproduces NITRO totals concentrations and ratios" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, .{
        .substrate_complex_count = 2,
        .structural_fraction_count = 3,
        .residue_fraction_count = 2,
        .structural_carbon_g_c = &.{ 1, 2, 3, 4, 5, 6 },
        .colonized_structural_carbon_g_c = &.{ 0.5, 1, 1.5, 2, 2.5, 3 },
        .microbial_residue_carbon_g_c = &.{ 1, 2, 3, 4 },
        .adsorbed_organic_carbon_g_c = &.{ 1, 2 },
        .adsorbed_acetate_carbon_g_c = &.{ 0.5, 1 },
        .dissolved_organic_carbon_g_c = &.{ 8, 0 },
        .dissolved_acetate_carbon_g_c = &.{ 2, 5 },
        .dissolved_organic_nitrogen_g_n = &.{ 0.8, 4 },
        .dissolved_organic_phosphorus_g_p = &.{ 0.4, 2 },
        .biologically_active_water_m3 = 10,
        .negligible_carbon_g_c = 1e-12,
        .negligible_water_m3 = 1e-12,
        .negligible_complex_fraction = 1e-12,
    });
    try std.testing.expectEqualSlices(f64, &.{ 6, 15 }, state.total_solid_carbon_g_c);
    try std.testing.expectEqualSlices(f64, &.{ 3, 7 }, state.total_microbial_residue_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 21), state.all_solid_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 10), state.all_microbial_residue_carbon_g_c);
    const first_carrier: f64 = 3 + 3 + 1.5;
    const second_carrier: f64 = 7.5 + 7 + 3;
    const total_carrier = first_carrier + second_carrier;
    try std.testing.expectApproxEqAbs(first_carrier / total_carrier, state.substrate_complex_fraction[0], 1e-15);
    try std.testing.expectApproxEqAbs(8 / (10 * first_carrier / total_carrier), state.dissolved_organic_carbon_concentration_g_c_per_m3[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), state.dissolved_nitrogen_per_organic_carbon_g_n_per_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), state.dissolved_phosphorus_per_organic_carbon_g_p_per_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), state.dissolved_organic_carbon_fraction[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), state.dissolved_acetate_carbon_fraction[0], 1e-15);
    try std.testing.expectEqual(@as(f64, 0), state.dissolved_organic_carbon_fraction[1]);
    try std.testing.expectEqual(@as(f64, 1), state.dissolved_acetate_carbon_fraction[1]);
}

test "zero colonized substrate retains exact source fallback fractions" {
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    try calculate(&state, .{
        .substrate_complex_count = 3,
        .structural_fraction_count = 1,
        .residue_fraction_count = 1,
        .structural_carbon_g_c = &.{ 0, 0, 0 },
        .colonized_structural_carbon_g_c = &.{ 0, 0, 0 },
        .microbial_residue_carbon_g_c = &.{ 0, 0, 0 },
        .adsorbed_organic_carbon_g_c = &.{ 0, 0, 0 },
        .adsorbed_acetate_carbon_g_c = &.{ 0, 0, 0 },
        .dissolved_organic_carbon_g_c = &.{ 0, 0, 0 },
        .dissolved_acetate_carbon_g_c = &.{ 0, 0, 0 },
        .dissolved_organic_nitrogen_g_n = &.{ 0, 0, 0 },
        .dissolved_organic_phosphorus_g_p = &.{ 0, 0, 0 },
        .biologically_active_water_m3 = 0,
        .negligible_carbon_g_c = 1e-12,
        .negligible_water_m3 = 1e-12,
        .negligible_complex_fraction = 1e-12,
    });
    try std.testing.expectEqualSlices(f64, &.{ 1, 1, 1 }, state.substrate_complex_fraction);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0 }, state.dissolved_organic_carbon_concentration_g_c_per_m3);
    try std.testing.expectEqualSlices(f64, &.{ 1, 1, 1 }, state.dissolved_acetate_carbon_fraction);
}

test "NITRO 481-539 uses strict independent carbon water and fraction thresholds" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, .{
        .substrate_complex_count = 1,
        .structural_fraction_count = 1,
        .residue_fraction_count = 1,
        .structural_carbon_g_c = &.{0},
        .colonized_structural_carbon_g_c = &.{0},
        .microbial_residue_carbon_g_c = &.{0.5},
        .adsorbed_organic_carbon_g_c = &.{0},
        .adsorbed_acetate_carbon_g_c = &.{0},
        .dissolved_organic_carbon_g_c = &.{0.5},
        .dissolved_acetate_carbon_g_c = &.{0.5},
        .dissolved_organic_nitrogen_g_n = &.{0.1},
        .dissolved_organic_phosphorus_g_p = &.{0.01},
        .biologically_active_water_m3 = 2,
        .negligible_carbon_g_c = 0.5,
        .negligible_water_m3 = 1,
        .negligible_complex_fraction = 1,
    });
    try std.testing.expectEqual(@as(f64, 1), state.substrate_complex_fraction[0]);
    try std.testing.expectEqual(
        @as(f64, 0.25),
        state.dissolved_organic_carbon_concentration_g_c_per_m3[0],
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        state.dissolved_nitrogen_per_organic_carbon_g_n_per_g_c[0],
    );
    try std.testing.expectEqual(
        @as(f64, 1),
        state.dissolved_acetate_carbon_fraction[0],
    );
}

test "invalid late complex input fails before replacing prior state" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.total_solid_carbon_g_c[0] = 9;
    try std.testing.expectError(error.InvalidSubstrateEnvironmentInput, calculate(&state, .{
        .substrate_complex_count = 2,
        .structural_fraction_count = 1,
        .residue_fraction_count = 1,
        .structural_carbon_g_c = &.{ 1, -1 },
        .colonized_structural_carbon_g_c = &.{ 0, 0 },
        .microbial_residue_carbon_g_c = &.{ 0, 0 },
        .adsorbed_organic_carbon_g_c = &.{ 0, 0 },
        .adsorbed_acetate_carbon_g_c = &.{ 0, 0 },
        .dissolved_organic_carbon_g_c = &.{ 0, 0 },
        .dissolved_acetate_carbon_g_c = &.{ 0, 0 },
        .dissolved_organic_nitrogen_g_n = &.{ 0, 0 },
        .dissolved_organic_phosphorus_g_p = &.{ 0, 0 },
        .biologically_active_water_m3 = 1,
        .negligible_carbon_g_c = 1e-12,
        .negligible_water_m3 = 1e-12,
        .negligible_complex_fraction = 1e-12,
    }));
    try std.testing.expectEqual(@as(f64, 9), state.total_solid_carbon_g_c[0]);
}

test "NITRO 349-410 aggregation rolls back on derived overflow" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.total_solid_carbon_g_c[0] = 9;
    try std.testing.expectError(error.NonFiniteSubstrateEnvironmentResult, calculate(&state, .{
        .substrate_complex_count = 1,
        .structural_fraction_count = 2,
        .residue_fraction_count = 1,
        .structural_carbon_g_c = &.{ std.math.floatMax(f64), std.math.floatMax(f64) },
        .colonized_structural_carbon_g_c = &.{ 0, 0 },
        .microbial_residue_carbon_g_c = &.{0},
        .adsorbed_organic_carbon_g_c = &.{0},
        .adsorbed_acetate_carbon_g_c = &.{0},
        .dissolved_organic_carbon_g_c = &.{0},
        .dissolved_acetate_carbon_g_c = &.{0},
        .dissolved_organic_nitrogen_g_n = &.{0},
        .dissolved_organic_phosphorus_g_p = &.{0},
        .biologically_active_water_m3 = 1,
        .negligible_carbon_g_c = 1e-12,
        .negligible_water_m3 = 1e-12,
        .negligible_complex_fraction = 1e-12,
    }));
    try std.testing.expectEqual(@as(f64, 9), state.total_solid_carbon_g_c[0]);
}
