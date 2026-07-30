const std = @import("std");

pub const Inputs = struct {
    enabled: []const bool,
    previous_total_non_band_nitrite_demand_g_n: []const f64,
    previous_total_band_nitrite_demand_g_n: []const f64,
    previous_non_band_potential_oxidation_g_n: []const f64,
    previous_band_potential_oxidation_g_n: []const f64,
    fraction_of_density_reference_biomass: []const f64,
    non_band_nitrate_solution_fraction: []const f64,
    band_nitrate_solution_fraction: []const f64,
    non_band_nitrite_concentration_g_n_per_m3: []const f64,
    band_nitrite_concentration_g_n_per_m3: []const f64,
    non_band_nitrite_g_n: []const f64,
    band_nitrite_g_n: []const f64,
    temperature_water_response: []const f64,
    combined_nutrient_limitation_fraction: []const f64,
    aqueous_carbon_dioxide_limitation_fraction: []const f64,
    active_biomass_g_c: []const f64,
    specific_nitrite_oxidation_g_n_per_g_c_h: []const f64,
    nitrite_half_saturation_g_n_per_m3: []const f64,
    respiration_requirement_g_c_per_g_n: []const f64,
    carbon_conversion_efficiency_g_c_per_g_n: []const f64,
    oxygen_per_respired_carbon_g_o_per_g_c: f64,
    oxygen_per_oxidized_nitrogen_g_o_per_g_n: f64,
    minimum_competition_fraction: f64,
    negligible_demand_g_n: f64,
    timestep_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    non_band_competition_fraction: []f64,
    band_competition_fraction: []f64,
    nitrite_unlimited_oxidation_g_n: []f64,
    non_band_substrate_response: []f64,
    band_substrate_response: []f64,
    combined_substrate_response: []f64,
    non_band_potential_oxidation_g_n: []f64,
    band_potential_oxidation_g_n: []f64,
    non_band_nitrite_oxidation_g_n: []f64,
    band_nitrite_oxidation_g_n: []f64,
    total_nitrite_oxidation_g_n: []f64,
    carbon_respiration_g_c: []f64,
    respiration_oxygen_demand_g_o: []f64,
    total_oxygen_demand_g_o: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidNitriteOxidationDimensions;
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

/// NITRO.F 1195--1274. Historical K=5,N=2 admission is runtime data.
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
        const non_band_competition = @max(
            inputs.minimum_competition_fraction,
            if (inputs.previous_total_non_band_nitrite_demand_g_n[unit] >
                inputs.negligible_demand_g_n)
                inputs.previous_non_band_potential_oxidation_g_n[unit] /
                    inputs.previous_total_non_band_nitrite_demand_g_n[unit]
            else
                inputs.fraction_of_density_reference_biomass[unit] *
                    inputs.non_band_nitrate_solution_fraction[unit],
        );
        const band_competition = @max(
            inputs.minimum_competition_fraction,
            if (inputs.previous_total_band_nitrite_demand_g_n[unit] >
                inputs.negligible_demand_g_n)
                inputs.previous_band_potential_oxidation_g_n[unit] /
                    inputs.previous_total_band_nitrite_demand_g_n[unit]
            else
                inputs.fraction_of_density_reference_biomass[unit] *
                    inputs.band_nitrate_solution_fraction[unit],
        );
        const unlimited =
            inputs.specific_nitrite_oxidation_g_n_per_g_c_h[unit] *
            inputs.temperature_water_response[unit] *
            inputs.combined_nutrient_limitation_fraction[unit] *
            inputs.aqueous_carbon_dioxide_limitation_fraction[unit] *
            inputs.active_biomass_g_c[unit] *
            inputs.timestep_h;
        const non_band_response =
            inputs.non_band_nitrate_solution_fraction[unit] *
            inputs.non_band_nitrite_concentration_g_n_per_m3[unit] /
            (inputs.non_band_nitrite_concentration_g_n_per_m3[unit] +
                inputs.nitrite_half_saturation_g_n_per_m3[unit]);
        const band_response =
            inputs.band_nitrate_solution_fraction[unit] *
            inputs.band_nitrite_concentration_g_n_per_m3[unit] /
            (inputs.band_nitrite_concentration_g_n_per_m3[unit] +
                inputs.nitrite_half_saturation_g_n_per_m3[unit]);
        const combined_substrate_response =
            non_band_response + band_response;
        const non_band_potential = unlimited * non_band_response;
        const band_potential = unlimited * band_response;
        const non_band_oxidation = @max(
            0,
            @min(
                non_band_potential,
                non_band_competition *
                    inputs.non_band_nitrite_g_n[unit] *
                    inputs.timestep_h,
            ),
        );
        const band_oxidation = @max(
            0,
            @min(
                band_potential,
                band_competition *
                    inputs.band_nitrite_g_n[unit] *
                    inputs.timestep_h,
            ),
        );
        const total_oxidation = non_band_oxidation + band_oxidation;
        const carbon_respiration = @max(
            0,
            total_oxidation *
                inputs.carbon_conversion_efficiency_g_c_per_g_n[unit] *
                inputs.respiration_requirement_g_c_per_g_n[unit],
        );
        const respiration_oxygen =
            inputs.oxygen_per_respired_carbon_g_o_per_g_c *
            carbon_respiration;
        const total_oxygen =
            respiration_oxygen +
            inputs.oxygen_per_oxidized_nitrogen_g_o_per_g_n *
                total_oxidation;

        state.non_band_competition_fraction[unit] = non_band_competition;
        state.band_competition_fraction[unit] = band_competition;
        state.nitrite_unlimited_oxidation_g_n[unit] = unlimited;
        state.non_band_substrate_response[unit] = non_band_response;
        state.band_substrate_response[unit] = band_response;
        state.combined_substrate_response[unit] =
            combined_substrate_response;
        state.non_band_potential_oxidation_g_n[unit] = non_band_potential;
        state.band_potential_oxidation_g_n[unit] = band_potential;
        state.non_band_nitrite_oxidation_g_n[unit] = non_band_oxidation;
        state.band_nitrite_oxidation_g_n[unit] = band_oxidation;
        state.total_nitrite_oxidation_g_n[unit] = total_oxidation;
        state.carbon_respiration_g_c[unit] = carbon_respiration;
        state.respiration_oxygen_demand_g_o[unit] = respiration_oxygen;
        state.total_oxygen_demand_g_o[unit] = total_oxygen;
    }
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state.*, field.name)) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteNitriteOxidationResult;
}

pub fn sourceEnabled(
    zero_based_complex: usize,
    zero_based_population: usize,
) bool {
    return zero_based_complex == 5 and zero_based_population == 1;
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.enabled.len != n)
        return error.InvalidNitriteOxidationDimensions;
    inline for (.{
        inputs.previous_total_non_band_nitrite_demand_g_n,
        inputs.previous_total_band_nitrite_demand_g_n,
        inputs.previous_non_band_potential_oxidation_g_n,
        inputs.previous_band_potential_oxidation_g_n,
        inputs.fraction_of_density_reference_biomass,
        inputs.non_band_nitrate_solution_fraction,
        inputs.band_nitrate_solution_fraction,
        inputs.non_band_nitrite_concentration_g_n_per_m3,
        inputs.band_nitrite_concentration_g_n_per_m3,
        inputs.non_band_nitrite_g_n,
        inputs.band_nitrite_g_n,
        inputs.temperature_water_response,
        inputs.combined_nutrient_limitation_fraction,
        inputs.aqueous_carbon_dioxide_limitation_fraction,
        inputs.active_biomass_g_c,
        inputs.specific_nitrite_oxidation_g_n_per_g_c_h,
        inputs.nitrite_half_saturation_g_n_per_m3,
        inputs.respiration_requirement_g_c_per_g_n,
        inputs.carbon_conversion_efficiency_g_c_per_g_n,
    }) |values| {
        if (values.len != n) return error.InvalidNitriteOxidationDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidNitriteOxidationInput;
    }
    for (0..n) |unit| if (inputs.enabled[unit]) {
        inline for (.{
            inputs.fraction_of_density_reference_biomass[unit],
            inputs.non_band_nitrate_solution_fraction[unit],
            inputs.band_nitrate_solution_fraction[unit],
            inputs.combined_nutrient_limitation_fraction[unit],
            inputs.aqueous_carbon_dioxide_limitation_fraction[unit],
        }) |fraction| if (fraction > 1)
            return error.InvalidNitriteOxidationInput;
        if (inputs.nitrite_half_saturation_g_n_per_m3[unit] <= 0 or
            inputs.respiration_requirement_g_c_per_g_n[unit] <= 0)
            return error.InvalidNitriteOxidationInput;
    };
    inline for (.{
        inputs.oxygen_per_respired_carbon_g_o_per_g_c,
        inputs.oxygen_per_oxidized_nitrogen_g_o_per_g_n,
        inputs.timestep_h,
    }) |value| if (!std.math.isFinite(value) or value <= 0)
        return error.InvalidNitriteOxidationInput;
    inline for (.{
        inputs.minimum_competition_fraction,
        inputs.negligible_demand_g_n,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidNitriteOxidationInput;
    if (inputs.minimum_competition_fraction > 1)
        return error.InvalidNitriteOxidationInput;
}

fn clear(state: *State) void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) @memset(@field(state, field.name), 0);
}

test "source nitrite oxidizer role reproduces NITRO selector" {
    try std.testing.expect(sourceEnabled(5, 1));
    try std.testing.expect(!sourceEnabled(5, 0));
    try std.testing.expect(!sourceEnabled(4, 1));
}

test "competition two-zone kinetics and oxygen demand reproduce NITRO" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, .{
        .enabled = &.{true},
        .previous_total_non_band_nitrite_demand_g_n = &.{4},
        .previous_total_band_nitrite_demand_g_n = &.{0},
        .previous_non_band_potential_oxidation_g_n = &.{1},
        .previous_band_potential_oxidation_g_n = &.{0},
        .fraction_of_density_reference_biomass = &.{0.2},
        .non_band_nitrate_solution_fraction = &.{0.75},
        .band_nitrate_solution_fraction = &.{0.25},
        .non_band_nitrite_concentration_g_n_per_m3 = &.{3},
        .band_nitrite_concentration_g_n_per_m3 = &.{1},
        .non_band_nitrite_g_n = &.{10},
        .band_nitrite_g_n = &.{10},
        .temperature_water_response = &.{0.5},
        .combined_nutrient_limitation_fraction = &.{0.8},
        .aqueous_carbon_dioxide_limitation_fraction = &.{0.5},
        .active_biomass_g_c = &.{10},
        .specific_nitrite_oxidation_g_n_per_g_c_h = &.{0.2},
        .nitrite_half_saturation_g_n_per_m3 = &.{1},
        .respiration_requirement_g_c_per_g_n = &.{0.4},
        .carbon_conversion_efficiency_g_c_per_g_n = &.{0.2},
        .oxygen_per_respired_carbon_g_o_per_g_c = 2.667,
        .oxygen_per_oxidized_nitrogen_g_o_per_g_n = 1.143,
        .minimum_competition_fraction = 0.01,
        .negligible_demand_g_n = 1e-12,
        .timestep_h = 1,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), state.non_band_competition_fraction[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), state.band_competition_fraction[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), state.nitrite_unlimited_oxidation_g_n[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.225), state.non_band_potential_oxidation_g_n[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), state.band_potential_oxidation_g_n[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.275), state.total_nitrite_oxidation_g_n[0], 1e-15);
    try std.testing.expect(state.total_oxygen_demand_g_o[0] > state.respiration_oxygen_demand_g_o[0]);
}

test "invalid late nitrite unit preserves prior state" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.total_nitrite_oxidation_g_n[0] = 7;
    const zeros = [_]f64{ 0, 0 };
    const ones = [_]f64{ 1, 1 };
    try std.testing.expectError(error.InvalidNitriteOxidationInput, calculate(&state, .{
        .enabled = &.{ true, true },
        .previous_total_non_band_nitrite_demand_g_n = &zeros,
        .previous_total_band_nitrite_demand_g_n = &zeros,
        .previous_non_band_potential_oxidation_g_n = &zeros,
        .previous_band_potential_oxidation_g_n = &zeros,
        .fraction_of_density_reference_biomass = &.{ 0.5, 2 },
        .non_band_nitrate_solution_fraction = &ones,
        .band_nitrate_solution_fraction = &zeros,
        .non_band_nitrite_concentration_g_n_per_m3 = &ones,
        .band_nitrite_concentration_g_n_per_m3 = &ones,
        .non_band_nitrite_g_n = &ones,
        .band_nitrite_g_n = &ones,
        .temperature_water_response = &ones,
        .combined_nutrient_limitation_fraction = &ones,
        .aqueous_carbon_dioxide_limitation_fraction = &ones,
        .active_biomass_g_c = &ones,
        .specific_nitrite_oxidation_g_n_per_g_c_h = &ones,
        .nitrite_half_saturation_g_n_per_m3 = &ones,
        .respiration_requirement_g_c_per_g_n = &ones,
        .carbon_conversion_efficiency_g_c_per_g_n = &ones,
        .oxygen_per_respired_carbon_g_o_per_g_c = 2.667,
        .oxygen_per_oxidized_nitrogen_g_o_per_g_n = 1.143,
        .minimum_competition_fraction = 0.01,
        .negligible_demand_g_n = 1e-12,
        .timestep_h = 1,
    }));
    try std.testing.expectEqual(@as(f64, 7), state.total_nitrite_oxidation_g_n[0]);
}

test "NITRO 1195-1274 derived overflow leaves nitrite oxidation unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.total_nitrite_oxidation_g_n[0] = 7;
    const zero = [_]f64{0};
    const one = [_]f64{1};
    try std.testing.expectError(
        error.NonFiniteNitriteOxidationResult,
        calculate(&state, .{
            .enabled = &.{true},
            .previous_total_non_band_nitrite_demand_g_n = &zero,
            .previous_total_band_nitrite_demand_g_n = &zero,
            .previous_non_band_potential_oxidation_g_n = &zero,
            .previous_band_potential_oxidation_g_n = &zero,
            .fraction_of_density_reference_biomass = &one,
            .non_band_nitrate_solution_fraction = &one,
            .band_nitrate_solution_fraction = &zero,
            .non_band_nitrite_concentration_g_n_per_m3 = &one,
            .band_nitrite_concentration_g_n_per_m3 = &zero,
            .non_band_nitrite_g_n = &one,
            .band_nitrite_g_n = &zero,
            .temperature_water_response = &one,
            .combined_nutrient_limitation_fraction = &one,
            .aqueous_carbon_dioxide_limitation_fraction = &one,
            .active_biomass_g_c = &.{std.math.floatMax(f64)},
            .specific_nitrite_oxidation_g_n_per_g_c_h = &.{std.math.floatMax(f64)},
            .nitrite_half_saturation_g_n_per_m3 = &one,
            .respiration_requirement_g_c_per_g_n = &one,
            .carbon_conversion_efficiency_g_c_per_g_n = &one,
            .oxygen_per_respired_carbon_g_o_per_g_c = 2.667,
            .oxygen_per_oxidized_nitrogen_g_o_per_g_n = 1.143,
            .minimum_competition_fraction = 0.01,
            .negligible_demand_g_n = 1e-12,
            .timestep_h = 1,
        }),
    );
    try std.testing.expectEqual(@as(f64, 7), state.total_nitrite_oxidation_g_n[0]);
}
