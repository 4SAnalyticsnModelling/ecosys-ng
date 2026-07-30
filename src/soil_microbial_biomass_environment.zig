const std = @import("std");

pub const Inputs = struct {
    substrate_complex_count: usize,
    population_count: usize,
    enabled: []const bool,
    density_reference_population: []const bool,
    labile_carbon_g_c: []const f64,
    labile_nitrogen_g_n: []const f64,
    labile_phosphorus_g_p: []const f64,
    resistant_carbon_g_c: []const f64,
    resistant_nitrogen_g_n: []const f64,
    target_nitrogen_per_carbon_g_n_per_g_c: []const f64,
    target_phosphorus_per_carbon_g_p_per_g_c: []const f64,
    labile_biomass_fraction: f64,
    resistant_biomass_fraction: f64,
    negligible_carbon_g_c: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    substrate_complex_count: usize,
    population_count: usize,
    actual_nitrogen_per_carbon_g_n_per_g_c: []f64,
    actual_phosphorus_per_carbon_g_p_per_g_c: []f64,
    active_biomass_g_c: []f64,
    nitrogen_limitation_fraction: []f64,
    phosphorus_limitation_fraction: []f64,
    combined_nutrient_limitation_fraction: []f64,
    active_resistant_biomass_g_c: []f64,
    active_resistant_nitrogen_g_n: []f64,
    active_carbon_by_complex_g_c: []f64,
    active_nitrogen_by_complex_g_n: []f64,
    active_phosphorus_by_complex_g_p: []f64,
    maximum_nitrogen_by_complex_g_n: []f64,
    maximum_phosphorus_by_complex_g_p: []f64,
    total_active_biomass_g_c: f64,
    density_reference_biomass_g_c: f64,

    pub fn init(
        allocator: std.mem.Allocator,
        substrate_complex_count: usize,
        population_count: usize,
    ) !State {
        if (substrate_complex_count == 0 or population_count == 0)
            return error.InvalidMicrobialBiomassDimensions;
        const unit_count = try std.math.mul(
            usize,
            substrate_complex_count,
            population_count,
        );
        var state: State = undefined;
        state.allocator = allocator;
        state.substrate_complex_count = substrate_complex_count;
        state.population_count = population_count;
        state.total_active_biomass_g_c = 0;
        state.density_reference_biomass_g_c = 0;
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
                const count =
                    if (std.mem.endsWith(u8, field.name, "_by_complex_g_c") or
                    std.mem.endsWith(u8, field.name, "_by_complex_g_n") or
                    std.mem.endsWith(u8, field.name, "_by_complex_g_p"))
                        substrate_complex_count
                    else
                        unit_count;
                @field(state, field.name) = try allocator.alloc(f64, count);
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

/// NITRO.F 418--479. Admission is supplied as runtime masks so the source
/// L=0 complex exclusions and K=5 population exclusions remain exact without
/// retaining fixed K=0:5 and N=1:7 storage.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(
        state.allocator,
        state.substrate_complex_count,
        state.population_count,
    );
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
        for (0..inputs.population_count) |population| {
            const unit = complex * inputs.population_count + population;
            if (!inputs.enabled[unit]) continue;
            const labile_carbon = inputs.labile_carbon_g_c[unit];
            const target_n =
                inputs.target_nitrogen_per_carbon_g_n_per_g_c[unit];
            const target_p =
                inputs.target_phosphorus_per_carbon_g_p_per_g_c[unit];
            state.actual_nitrogen_per_carbon_g_n_per_g_c[unit] =
                if (labile_carbon > inputs.negligible_carbon_g_c)
                    inputs.labile_nitrogen_g_n[unit] / labile_carbon
                else
                    target_n;
            state.actual_phosphorus_per_carbon_g_p_per_g_c[unit] =
                if (labile_carbon > inputs.negligible_carbon_g_c)
                    inputs.labile_phosphorus_g_p[unit] / labile_carbon
                else
                    target_p;
            state.active_biomass_g_c[unit] =
                labile_carbon / inputs.labile_biomass_fraction;
            state.nitrogen_limitation_fraction[unit] = @min(
                1,
                @max(0.1, std.math.pow(
                    f64,
                    state.actual_nitrogen_per_carbon_g_n_per_g_c[unit] /
                        target_n,
                    0.25,
                )),
            );
            state.phosphorus_limitation_fraction[unit] = @min(
                1,
                @max(0.1, std.math.pow(
                    f64,
                    state.actual_phosphorus_per_carbon_g_p_per_g_c[unit] /
                        target_p,
                    0.25,
                )),
            );
            state.combined_nutrient_limitation_fraction[unit] = @min(
                state.nitrogen_limitation_fraction[unit],
                state.phosphorus_limitation_fraction[unit],
            );
            state.total_active_biomass_g_c += state.active_biomass_g_c[unit];
            if (inputs.density_reference_population[unit])
                state.density_reference_biomass_g_c +=
                    state.active_biomass_g_c[unit];
            state.active_resistant_biomass_g_c[unit] = @min(
                state.active_biomass_g_c[unit] *
                    inputs.resistant_biomass_fraction,
                inputs.resistant_carbon_g_c[unit],
            );
            state.active_resistant_nitrogen_g_n[unit] =
                if (inputs.resistant_carbon_g_c[unit] >
                inputs.negligible_carbon_g_c)
                    state.active_resistant_biomass_g_c[unit] /
                        inputs.resistant_carbon_g_c[unit] *
                        inputs.resistant_nitrogen_g_n[unit]
                else
                    0;
        }
    }

    for (0..inputs.substrate_complex_count) |complex| {
        for (0..inputs.population_count) |population| {
            const unit = complex * inputs.population_count + population;
            const active_carbon = state.active_biomass_g_c[unit];
            state.active_carbon_by_complex_g_c[complex] += active_carbon;
            state.active_nitrogen_by_complex_g_n[complex] +=
                active_carbon *
                state.actual_nitrogen_per_carbon_g_n_per_g_c[unit];
            state.active_phosphorus_by_complex_g_p[complex] +=
                active_carbon *
                state.actual_phosphorus_per_carbon_g_p_per_g_c[unit];
            state.maximum_nitrogen_by_complex_g_n[complex] +=
                active_carbon *
                inputs.target_nitrogen_per_carbon_g_n_per_g_c[unit];
            state.maximum_phosphorus_by_complex_g_p[complex] +=
                active_carbon *
                inputs.target_phosphorus_per_carbon_g_p_per_g_c[unit];
        }
    }
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (field.type == f64 and
            (!std.math.isFinite(@field(state.*, field.name)) or
                @field(state.*, field.name) < 0))
            return error.NonFiniteMicrobialBiomassResult;
        if (field.type == []f64) for (@field(state.*, field.name)) |value| {
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteMicrobialBiomassResult;
        };
    }
}

/// Exact NITRO source admission: litter excludes complexes K=3,4; the
/// autotrophic complex K=5 admits populations N=1,2,3,5 only.
pub fn sourceEnabled(
    is_surface_litter: bool,
    zero_based_complex: usize,
    zero_based_population: usize,
) bool {
    if (is_surface_litter and
        (zero_based_complex == 3 or zero_based_complex == 4)) return false;
    return zero_based_complex != 5 or
        zero_based_population <= 2 or zero_based_population == 4;
}

/// Exact TOMA companion TOMN selector at NITRO.F 451--453.
pub fn sourceDensityReference(
    zero_based_complex: usize,
    zero_based_population: usize,
) bool {
    return (zero_based_complex <= 4 and zero_based_population == 1) or
        (zero_based_complex == 5 and zero_based_population == 0);
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.substrate_complex_count == 0 or inputs.population_count == 0 or
        state.substrate_complex_count != inputs.substrate_complex_count or
        state.population_count != inputs.population_count)
        return error.InvalidMicrobialBiomassDimensions;
    const count = try std.math.mul(
        usize,
        inputs.substrate_complex_count,
        inputs.population_count,
    );
    inline for (.{
        inputs.enabled,
        inputs.density_reference_population,
    }) |values| if (values.len != count)
        return error.InvalidMicrobialBiomassDimensions;
    inline for (.{
        inputs.labile_carbon_g_c,
        inputs.labile_nitrogen_g_n,
        inputs.labile_phosphorus_g_p,
        inputs.resistant_carbon_g_c,
        inputs.resistant_nitrogen_g_n,
        inputs.target_nitrogen_per_carbon_g_n_per_g_c,
        inputs.target_phosphorus_per_carbon_g_p_per_g_c,
    }) |values| {
        if (values.len != count) return error.InvalidMicrobialBiomassDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidMicrobialBiomassInput;
    }
    for (0..count) |unit| if (inputs.enabled[unit] and
        (inputs.target_nitrogen_per_carbon_g_n_per_g_c[unit] <= 0 or
            inputs.target_phosphorus_per_carbon_g_p_per_g_c[unit] <= 0))
        return error.InvalidMicrobialBiomassInput;
    inline for (.{
        inputs.labile_biomass_fraction,
        inputs.resistant_biomass_fraction,
        inputs.negligible_carbon_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidMicrobialBiomassInput;
    if (inputs.labile_biomass_fraction <= 0 or
        inputs.labile_biomass_fraction > 1 or
        inputs.resistant_biomass_fraction > 1)
        return error.InvalidMicrobialBiomassInput;
}

fn clear(state: *State) void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) @memset(@field(state, field.name), 0);
    state.total_active_biomass_g_c = 0;
    state.density_reference_biomass_g_c = 0;
}

test "source masks reproduce surface and autotrophic NITRO admissions" {
    try std.testing.expect(sourceEnabled(true, 2, 6));
    try std.testing.expect(!sourceEnabled(true, 3, 0));
    try std.testing.expect(!sourceEnabled(true, 4, 0));
    try std.testing.expect(sourceEnabled(false, 4, 6));
    try std.testing.expect(sourceEnabled(false, 5, 0));
    try std.testing.expect(sourceEnabled(false, 5, 2));
    try std.testing.expect(sourceEnabled(false, 5, 4));
    try std.testing.expect(!sourceEnabled(false, 5, 3));
    try std.testing.expect(!sourceEnabled(false, 5, 6));
    try std.testing.expect(sourceDensityReference(0, 1));
    try std.testing.expect(sourceDensityReference(5, 0));
    try std.testing.expect(!sourceDensityReference(5, 1));
}

test "active biomass and nutrient factors reproduce NITRO equations on runtime axes" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    try calculate(&state, .{
        .substrate_complex_count = 2,
        .population_count = 2,
        .enabled = &.{ true, true, true, false },
        .density_reference_population = &.{ false, true, true, false },
        .labile_carbon_g_c = &.{ 2, 1, 0, 99 },
        .labile_nitrogen_g_n = &.{ 0.1, 0.2, 0, 99 },
        .labile_phosphorus_g_p = &.{ 0.01, 0.001, 0, 99 },
        .resistant_carbon_g_c = &.{ 5, 1, 2, 99 },
        .resistant_nitrogen_g_n = &.{ 0.5, 0.4, 0.2, 99 },
        .target_nitrogen_per_carbon_g_n_per_g_c = &.{ 0.1, 0.1, 0.2, 0 },
        .target_phosphorus_per_carbon_g_p_per_g_c = &.{ 0.01, 0.01, 0.02, 0 },
        .labile_biomass_fraction = 0.5,
        .resistant_biomass_fraction = 0.5,
        .negligible_carbon_g_c = 1e-12,
    });
    try std.testing.expectEqualSlices(f64, &.{ 4, 2, 0, 0 }, state.active_biomass_g_c);
    try std.testing.expectApproxEqAbs(std.math.pow(f64, 0.5, 0.25), state.nitrogen_limitation_fraction[0], 1e-15);
    try std.testing.expectEqual(@as(f64, 1), state.nitrogen_limitation_fraction[1]);
    try std.testing.expectApproxEqAbs(std.math.pow(f64, 0.1, 0.25), state.phosphorus_limitation_fraction[1], 1e-15);
    try std.testing.expectEqual(@as(f64, 1), state.combined_nutrient_limitation_fraction[2]);
    try std.testing.expectEqualSlices(f64, &.{ 2, 1, 0, 0 }, state.active_resistant_biomass_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), state.active_resistant_nitrogen_g_n[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), state.active_resistant_nitrogen_g_n[1], 1e-15);
    try std.testing.expectEqualSlices(f64, &.{ 6, 0 }, state.active_carbon_by_complex_g_c);
    try std.testing.expectEqual(@as(f64, 6), state.total_active_biomass_g_c);
    try std.testing.expectEqual(@as(f64, 2), state.density_reference_biomass_g_c);
}

test "invalid late population leaves prior biomass state unchanged" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.active_biomass_g_c[0] = 7;
    try std.testing.expectError(error.InvalidMicrobialBiomassInput, calculate(&state, .{
        .substrate_complex_count = 1,
        .population_count = 2,
        .enabled = &.{ true, true },
        .density_reference_population = &.{ false, false },
        .labile_carbon_g_c = &.{ 1, 1 },
        .labile_nitrogen_g_n = &.{ 0.1, 0.1 },
        .labile_phosphorus_g_p = &.{ 0.01, 0.01 },
        .resistant_carbon_g_c = &.{ 1, 1 },
        .resistant_nitrogen_g_n = &.{ 0.1, 0.1 },
        .target_nitrogen_per_carbon_g_n_per_g_c = &.{ 0.1, 0 },
        .target_phosphorus_per_carbon_g_p_per_g_c = &.{ 0.01, 0.01 },
        .labile_biomass_fraction = 0.5,
        .resistant_biomass_fraction = 0.5,
        .negligible_carbon_g_c = 1e-12,
    }));
    try std.testing.expectEqual(@as(f64, 7), state.active_biomass_g_c[0]);
}

test "NITRO 418-479 derived overflow leaves prior biomass state unchanged" {
    var state = try State.init(std.testing.allocator, 1, 1);
    defer state.deinit();
    state.active_biomass_g_c[0] = 7;
    try std.testing.expectError(error.NonFiniteMicrobialBiomassResult, calculate(&state, .{
        .substrate_complex_count = 1,
        .population_count = 1,
        .enabled = &.{true},
        .density_reference_population = &.{false},
        .labile_carbon_g_c = &.{std.math.floatMax(f64)},
        .labile_nitrogen_g_n = &.{0},
        .labile_phosphorus_g_p = &.{0},
        .resistant_carbon_g_c = &.{0},
        .resistant_nitrogen_g_n = &.{0},
        .target_nitrogen_per_carbon_g_n_per_g_c = &.{0.1},
        .target_phosphorus_per_carbon_g_p_per_g_c = &.{0.01},
        .labile_biomass_fraction = 0.5,
        .resistant_biomass_fraction = 0.5,
        .negligible_carbon_g_c = 1e-12,
    }));
    try std.testing.expectEqual(@as(f64, 7), state.active_biomass_g_c[0]);
}
