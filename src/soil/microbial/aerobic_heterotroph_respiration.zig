const std = @import("std");

pub const Inputs = struct {
    enabled: []const bool,
    dissolved_organic_carbon_concentration_g_c_per_m3: []const f64,
    dissolved_acetate_carbon_concentration_g_c_per_m3: []const f64,
    dissolved_organic_carbon_g_c: []const f64,
    dissolved_acetate_carbon_g_c: []const f64,
    dissolved_organic_carbon_fraction: []const f64,
    dissolved_acetate_carbon_fraction: []const f64,
    dissolved_organic_carbon_competition_fraction: []const f64,
    dissolved_acetate_carbon_competition_fraction: []const f64,
    combined_nutrient_limitation_fraction: []const f64,
    water_response: []const f64,
    active_biomass_g_c: []const f64,
    growth_temperature_response: []const f64,
    specific_oxidation_rate_g_c_per_g_c_h: []const f64,
    dissolved_organic_carbon_half_saturation_g_c_per_m3: []const f64,
    dissolved_acetate_carbon_half_saturation_g_c_per_m3: []const f64,
    dissolved_organic_carbon_respiration_requirement_g_c_per_g_c: []const f64,
    dissolved_acetate_carbon_respiration_requirement_g_c_per_g_c: []const f64,
    oxygen_per_respired_carbon_g_o_per_g_c: f64,
    timestep_h: f64,
    negligible_respiration_g_c: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    dissolved_organic_carbon_monod_fraction: []f64,
    dissolved_acetate_carbon_monod_fraction: []f64,
    combined_substrate_response: []f64,
    substrate_unlimited_respiration_g_c: []f64,
    dissolved_organic_carbon_microbial_limit_g_c: []f64,
    dissolved_acetate_carbon_microbial_limit_g_c: []f64,
    dissolved_organic_carbon_supply_limit_g_c: []f64,
    dissolved_acetate_carbon_supply_limit_g_c: []f64,
    dissolved_organic_carbon_respiration_g_c: []f64,
    dissolved_acetate_carbon_respiration_g_c: []f64,
    total_respiration_g_c: []f64,
    dissolved_organic_carbon_respiration_fraction: []f64,
    dissolved_acetate_carbon_respiration_fraction: []f64,
    weighted_respiration_requirement_g_c_per_g_c: []f64,
    oxygen_demand_g_o: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidAerobicRespirationDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
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
                @field(state, field.name) = try allocator.alloc(f64, unit_count);
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

/// NITRO.F 789--885 (source populations N=1,2,3,6 in K<=4). Population
/// admission and its EO2Q/EO2A requirements are runtime data, retaining the
/// exact equations without fixed population roles or complex dimensions.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.unit_count);
    defer staged.deinit();
    calculateValidated(&staged, inputs);
    try validateResult(&staged);
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64)
            @memcpy(@field(state, field.name), @field(staged, field.name));
}

fn calculateValidated(state: *State, inputs: Inputs) void {
    clear(state);
    for (0..state.unit_count) |unit| {
        if (!inputs.enabled[unit]) continue;
        const doc_concentration =
            inputs.dissolved_organic_carbon_concentration_g_c_per_m3[unit];
        const acetate_concentration =
            inputs.dissolved_acetate_carbon_concentration_g_c_per_m3[unit];
        const doc_monod = doc_concentration /
            (doc_concentration +
                inputs.dissolved_organic_carbon_half_saturation_g_c_per_m3[unit]);
        const acetate_monod = acetate_concentration /
            (acetate_concentration +
                inputs.dissolved_acetate_carbon_half_saturation_g_c_per_m3[unit]);
        const combined_substrate_response =
            inputs.dissolved_organic_carbon_fraction[unit] * doc_monod +
            inputs.dissolved_acetate_carbon_fraction[unit] * acetate_monod;
        const unlimited = @max(
            0,
            inputs.specific_oxidation_rate_g_c_per_g_c_h[unit] *
                inputs.combined_nutrient_limitation_fraction[unit] *
                inputs.water_response[unit] *
                inputs.active_biomass_g_c[unit],
        ) * inputs.timestep_h;
        const doc_microbial_limit =
            unlimited * doc_monod *
            inputs.dissolved_organic_carbon_fraction[unit] *
            inputs.growth_temperature_response[unit];
        const acetate_microbial_limit =
            unlimited * acetate_monod *
            inputs.dissolved_acetate_carbon_fraction[unit] *
            inputs.growth_temperature_response[unit];
        const doc_supply_limit = @max(
            0,
            inputs.dissolved_organic_carbon_g_c[unit] *
                inputs.dissolved_organic_carbon_competition_fraction[unit] *
                inputs.dissolved_organic_carbon_respiration_requirement_g_c_per_g_c[unit] *
                inputs.timestep_h,
        );
        const acetate_supply_limit = @max(
            0,
            inputs.dissolved_acetate_carbon_g_c[unit] *
                inputs.dissolved_acetate_carbon_competition_fraction[unit] *
                inputs.dissolved_acetate_carbon_respiration_requirement_g_c_per_g_c[unit] *
                inputs.timestep_h,
        );
        const doc_respiration = @min(doc_supply_limit, doc_microbial_limit);
        const acetate_respiration =
            @min(acetate_supply_limit, acetate_microbial_limit);
        const total_respiration = doc_respiration + acetate_respiration;
        const doc_fraction =
            if (total_respiration > inputs.negligible_respiration_g_c)
                doc_respiration / total_respiration
            else
                1;
        const acetate_fraction =
            if (total_respiration > inputs.negligible_respiration_g_c)
                acetate_respiration / total_respiration
            else
                0;

        state.dissolved_organic_carbon_monod_fraction[unit] = doc_monod;
        state.dissolved_acetate_carbon_monod_fraction[unit] = acetate_monod;
        state.combined_substrate_response[unit] =
            combined_substrate_response;
        state.substrate_unlimited_respiration_g_c[unit] = unlimited;
        state.dissolved_organic_carbon_microbial_limit_g_c[unit] =
            doc_microbial_limit;
        state.dissolved_acetate_carbon_microbial_limit_g_c[unit] =
            acetate_microbial_limit;
        state.dissolved_organic_carbon_supply_limit_g_c[unit] =
            doc_supply_limit;
        state.dissolved_acetate_carbon_supply_limit_g_c[unit] =
            acetate_supply_limit;
        state.dissolved_organic_carbon_respiration_g_c[unit] = doc_respiration;
        state.dissolved_acetate_carbon_respiration_g_c[unit] =
            acetate_respiration;
        state.total_respiration_g_c[unit] = total_respiration;
        state.dissolved_organic_carbon_respiration_fraction[unit] =
            doc_fraction;
        state.dissolved_acetate_carbon_respiration_fraction[unit] =
            acetate_fraction;
        state.weighted_respiration_requirement_g_c_per_g_c[unit] =
            inputs.dissolved_organic_carbon_respiration_requirement_g_c_per_g_c[unit] *
            doc_fraction +
            inputs.dissolved_acetate_carbon_respiration_requirement_g_c_per_g_c[unit] *
                acetate_fraction;
        state.oxygen_demand_g_o[unit] =
            inputs.oxygen_per_respired_carbon_g_o_per_g_c *
            total_respiration;
    }
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state.*, field.name)) |value| {
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteAerobicRespirationResult;
        };
}

/// Exact source aerobic heterotroph selector for zero-based K,N.
pub fn sourceEnabled(
    zero_based_complex: usize,
    zero_based_population: usize,
) bool {
    return zero_based_complex <= 4 and
        (zero_based_population <= 2 or zero_based_population == 5);
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.enabled.len != n)
        return error.InvalidAerobicRespirationDimensions;
    inline for (.{
        inputs.dissolved_organic_carbon_concentration_g_c_per_m3,
        inputs.dissolved_acetate_carbon_concentration_g_c_per_m3,
        inputs.dissolved_organic_carbon_g_c,
        inputs.dissolved_acetate_carbon_g_c,
        inputs.dissolved_organic_carbon_fraction,
        inputs.dissolved_acetate_carbon_fraction,
        inputs.dissolved_organic_carbon_competition_fraction,
        inputs.dissolved_acetate_carbon_competition_fraction,
        inputs.combined_nutrient_limitation_fraction,
        inputs.water_response,
        inputs.active_biomass_g_c,
        inputs.growth_temperature_response,
        inputs.specific_oxidation_rate_g_c_per_g_c_h,
        inputs.dissolved_organic_carbon_half_saturation_g_c_per_m3,
        inputs.dissolved_acetate_carbon_half_saturation_g_c_per_m3,
        inputs.dissolved_organic_carbon_respiration_requirement_g_c_per_g_c,
        inputs.dissolved_acetate_carbon_respiration_requirement_g_c_per_g_c,
    }) |values| {
        if (values.len != n) return error.InvalidAerobicRespirationDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidAerobicRespirationInput;
    }
    for (0..n) |unit| if (inputs.enabled[unit]) {
        if (inputs.dissolved_organic_carbon_fraction[unit] > 1 or
            inputs.dissolved_acetate_carbon_fraction[unit] > 1 or
            @abs(inputs.dissolved_organic_carbon_fraction[unit] +
                inputs.dissolved_acetate_carbon_fraction[unit] - 1) >
                64 * std.math.floatEps(f64) or
            inputs.combined_nutrient_limitation_fraction[unit] > 1 or
            inputs.dissolved_organic_carbon_half_saturation_g_c_per_m3[unit] <= 0 or
            inputs.dissolved_acetate_carbon_half_saturation_g_c_per_m3[unit] <= 0)
            return error.InvalidAerobicRespirationInput;
    };
    inline for (.{
        inputs.oxygen_per_respired_carbon_g_o_per_g_c,
        inputs.timestep_h,
        inputs.negligible_respiration_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidAerobicRespirationInput;
    if (inputs.oxygen_per_respired_carbon_g_o_per_g_c <= 0 or
        inputs.timestep_h <= 0)
        return error.InvalidAerobicRespirationInput;
}

fn clear(state: *State) void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) @memset(@field(state, field.name), 0);
}

test "source aerobic population roles reproduce NITRO selector" {
    try std.testing.expect(sourceEnabled(0, 0));
    try std.testing.expect(sourceEnabled(4, 2));
    try std.testing.expect(sourceEnabled(2, 5));
    try std.testing.expect(!sourceEnabled(2, 3));
    try std.testing.expect(!sourceEnabled(5, 0));
}

test "DOC acetate caps and oxygen demand reproduce NITRO aerobic equations" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, .{
        .enabled = &.{ true, false },
        .dissolved_organic_carbon_concentration_g_c_per_m3 = &.{ 3, 0 },
        .dissolved_acetate_carbon_concentration_g_c_per_m3 = &.{ 1, 0 },
        .dissolved_organic_carbon_g_c = &.{ 4, 0 },
        .dissolved_acetate_carbon_g_c = &.{ 2, 0 },
        .dissolved_organic_carbon_fraction = &.{ 0.75, 0 },
        .dissolved_acetate_carbon_fraction = &.{ 0.25, 0 },
        .dissolved_organic_carbon_competition_fraction = &.{ 0.5, 0 },
        .dissolved_acetate_carbon_competition_fraction = &.{ 0.5, 0 },
        .combined_nutrient_limitation_fraction = &.{ 0.8, 0 },
        .water_response = &.{ 0.5, 0 },
        .active_biomass_g_c = &.{ 10, 0 },
        .growth_temperature_response = &.{ 2, 0 },
        .specific_oxidation_rate_g_c_per_g_c_h = &.{ 0.2, 0 },
        .dissolved_organic_carbon_half_saturation_g_c_per_m3 = &.{ 1, 1 },
        .dissolved_acetate_carbon_half_saturation_g_c_per_m3 = &.{ 1, 1 },
        .dissolved_organic_carbon_respiration_requirement_g_c_per_g_c = &.{ 0.4, 0 },
        .dissolved_acetate_carbon_respiration_requirement_g_c_per_g_c = &.{ 0.2, 0 },
        .oxygen_per_respired_carbon_g_o_per_g_c = 2.667,
        .timestep_h = 1,
        .negligible_respiration_g_c = 1e-12,
    });
    // RGOCY=.8; RGOCZ=.8*(3/4)*.75*2=.9; supply=.8 => .8.
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), state.substrate_unlimited_respiration_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.9), state.dissolved_organic_carbon_microbial_limit_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), state.dissolved_organic_carbon_respiration_g_c[0], 1e-15);
    // RGOAZ=.8*(1/2)*.25*2=.2; supply=.2 => .2.
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), state.dissolved_acetate_carbon_respiration_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), state.total_respiration_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), state.dissolved_organic_carbon_respiration_fraction[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.36), state.weighted_respiration_requirement_g_c_per_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2.667), state.oxygen_demand_g_o[0], 1e-15);
}

test "zero respiration uses exact source DOC fraction fallback" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, .{
        .enabled = &.{true},
        .dissolved_organic_carbon_concentration_g_c_per_m3 = &.{0},
        .dissolved_acetate_carbon_concentration_g_c_per_m3 = &.{0},
        .dissolved_organic_carbon_g_c = &.{0},
        .dissolved_acetate_carbon_g_c = &.{0},
        .dissolved_organic_carbon_fraction = &.{0},
        .dissolved_acetate_carbon_fraction = &.{1},
        .dissolved_organic_carbon_competition_fraction = &.{0},
        .dissolved_acetate_carbon_competition_fraction = &.{0},
        .combined_nutrient_limitation_fraction = &.{1},
        .water_response = &.{1},
        .active_biomass_g_c = &.{0},
        .growth_temperature_response = &.{1},
        .specific_oxidation_rate_g_c_per_g_c_h = &.{1},
        .dissolved_organic_carbon_half_saturation_g_c_per_m3 = &.{1},
        .dissolved_acetate_carbon_half_saturation_g_c_per_m3 = &.{1},
        .dissolved_organic_carbon_respiration_requirement_g_c_per_g_c = &.{0.4},
        .dissolved_acetate_carbon_respiration_requirement_g_c_per_g_c = &.{0.2},
        .oxygen_per_respired_carbon_g_o_per_g_c = 2.667,
        .timestep_h = 1,
        .negligible_respiration_g_c = 1e-12,
    });
    try std.testing.expectEqual(@as(f64, 1), state.dissolved_organic_carbon_respiration_fraction[0]);
    try std.testing.expectEqual(@as(f64, 0), state.dissolved_acetate_carbon_respiration_fraction[0]);
}

test "invalid late unit leaves respiration state unchanged" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.total_respiration_g_c[0] = 7;
    var zeros = [_]f64{ 0, 0 };
    var ones = [_]f64{ 1, 1 };
    var fractions = [_]f64{ 1, 2 };
    try std.testing.expectError(error.InvalidAerobicRespirationInput, calculate(&state, .{
        .enabled = &.{ true, true },
        .dissolved_organic_carbon_concentration_g_c_per_m3 = &zeros,
        .dissolved_acetate_carbon_concentration_g_c_per_m3 = &zeros,
        .dissolved_organic_carbon_g_c = &zeros,
        .dissolved_acetate_carbon_g_c = &zeros,
        .dissolved_organic_carbon_fraction = &fractions,
        .dissolved_acetate_carbon_fraction = &zeros,
        .dissolved_organic_carbon_competition_fraction = &zeros,
        .dissolved_acetate_carbon_competition_fraction = &zeros,
        .combined_nutrient_limitation_fraction = &ones,
        .water_response = &ones,
        .active_biomass_g_c = &zeros,
        .growth_temperature_response = &ones,
        .specific_oxidation_rate_g_c_per_g_c_h = &ones,
        .dissolved_organic_carbon_half_saturation_g_c_per_m3 = &ones,
        .dissolved_acetate_carbon_half_saturation_g_c_per_m3 = &ones,
        .dissolved_organic_carbon_respiration_requirement_g_c_per_g_c = &ones,
        .dissolved_acetate_carbon_respiration_requirement_g_c_per_g_c = &ones,
        .oxygen_per_respired_carbon_g_o_per_g_c = 2.667,
        .timestep_h = 1,
        .negligible_respiration_g_c = 1e-12,
    }));
    try std.testing.expectEqual(@as(f64, 7), state.total_respiration_g_c[0]);
}

test "NITRO 789-885 derived overflow leaves aerobic state unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.total_respiration_g_c[0] = 7;
    const one = [_]f64{1};
    const zero = [_]f64{0};
    const maximum = [_]f64{std.math.floatMax(f64)};
    try std.testing.expectError(
        error.NonFiniteAerobicRespirationResult,
        calculate(&state, .{
            .enabled = &.{true},
            .dissolved_organic_carbon_concentration_g_c_per_m3 = &one,
            .dissolved_acetate_carbon_concentration_g_c_per_m3 = &zero,
            .dissolved_organic_carbon_g_c = &one,
            .dissolved_acetate_carbon_g_c = &zero,
            .dissolved_organic_carbon_fraction = &one,
            .dissolved_acetate_carbon_fraction = &zero,
            .dissolved_organic_carbon_competition_fraction = &one,
            .dissolved_acetate_carbon_competition_fraction = &zero,
            .combined_nutrient_limitation_fraction = &one,
            .water_response = &one,
            .active_biomass_g_c = &maximum,
            .growth_temperature_response = &one,
            .specific_oxidation_rate_g_c_per_g_c_h = &maximum,
            .dissolved_organic_carbon_half_saturation_g_c_per_m3 = &one,
            .dissolved_acetate_carbon_half_saturation_g_c_per_m3 = &one,
            .dissolved_organic_carbon_respiration_requirement_g_c_per_g_c = &one,
            .dissolved_acetate_carbon_respiration_requirement_g_c_per_g_c = &one,
            .oxygen_per_respired_carbon_g_o_per_g_c = 2.667,
            .timestep_h = 1,
            .negligible_respiration_g_c = 1e-12,
        }),
    );
    try std.testing.expectEqual(@as(f64, 7), state.total_respiration_g_c[0]);
}
