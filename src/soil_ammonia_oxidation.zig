const std = @import("std");

pub const Inputs = struct {
    enabled: []const bool,
    previous_total_non_band_ammonium_demand_g_n: []const f64,
    previous_total_band_ammonium_demand_g_n: []const f64,
    previous_non_band_potential_oxidation_g_n: []const f64,
    previous_band_potential_oxidation_g_n: []const f64,
    fraction_of_total_active_biomass: []const f64,
    non_band_ammonium_solution_fraction: []const f64,
    band_ammonium_solution_fraction: []const f64,
    initial_nitrification_inhibitor_activity: []const f64,
    current_nitrification_inhibitor_activity: []const f64,
    non_band_ammonium_concentration_g_n_per_m3: []const f64,
    band_ammonium_concentration_g_n_per_m3: []const f64,
    non_band_ammonia_concentration_g_n_per_m3: []const f64,
    band_ammonia_concentration_g_n_per_m3: []const f64,
    non_band_ammonium_g_n: []const f64,
    band_ammonium_g_n: []const f64,
    temperature_water_response: []const f64,
    growth_temperature_response: []const f64,
    combined_nutrient_limitation_fraction: []const f64,
    aqueous_carbon_dioxide_limitation_fraction: []const f64,
    active_biomass_g_c: []const f64,
    specific_ammonia_oxidation_g_n_per_g_c_h: []const f64,
    inhibitor_decay_rate_per_h: []const f64,
    ammonium_inhibition_of_decay_g_n_per_m3: []const f64,
    ammonia_inhibition_of_oxidation_g_n_per_m3: []const f64,
    ammonium_half_saturation_g_n_per_m3: []const f64,
    respiration_requirement_g_c_per_g_n: []const f64,
    carbon_conversion_efficiency_g_c_per_g_n: []const f64,
    oxygen_per_respired_carbon_g_o_per_g_c: f64,
    oxygen_per_oxidized_nitrogen_g_o_per_g_n: f64,
    minimum_competition_fraction: f64,
    negligible_demand_g_n: f64,
    negligible_inhibitor_activity: f64,
    timestep_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    non_band_competition_fraction: []f64,
    band_competition_fraction: []f64,
    next_nitrification_inhibitor_activity: []f64,
    non_band_inhibitor_response: []f64,
    band_inhibitor_response: []f64,
    ammonium_unlimited_oxidation_g_n: []f64,
    non_band_ammonia_inhibited_oxidation_g_n: []f64,
    band_ammonia_inhibited_oxidation_g_n: []f64,
    non_band_ammonium_monod_fraction: []f64,
    band_ammonium_monod_fraction: []f64,
    non_band_potential_oxidation_g_n: []f64,
    band_potential_oxidation_g_n: []f64,
    non_band_ammonia_oxidation_g_n: []f64,
    band_ammonia_oxidation_g_n: []f64,
    total_ammonia_oxidation_g_n: []f64,
    carbon_respiration_g_c: []f64,
    respiration_oxygen_demand_g_o: []f64,
    total_oxygen_demand_g_o: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidAmmoniaOxidationDimensions;
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

/// NITRO.F 1074--1180. The historical K=5,N=1 admission is a runtime mask;
/// all prior-demand, inhibitor and zone coordinates are explicit.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(state.allocator, state.unit_count);
    defer staged.deinit();
    try calculateValidated(&staged, inputs);
    try validateResult(&staged);
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64)
            @memcpy(@field(state, field.name), @field(staged, field.name));
}

fn calculateValidated(state: *State, inputs: Inputs) !void {
    clear(state);
    for (0..state.unit_count) |unit| {
        if (!inputs.enabled[unit]) continue;
        const non_band_competition = @max(
            inputs.minimum_competition_fraction,
            if (inputs.previous_total_non_band_ammonium_demand_g_n[unit] >
                inputs.negligible_demand_g_n)
                inputs.previous_non_band_potential_oxidation_g_n[unit] /
                    inputs.previous_total_non_band_ammonium_demand_g_n[unit]
            else
                inputs.non_band_ammonium_solution_fraction[unit] *
                    inputs.fraction_of_total_active_biomass[unit],
        );
        const band_competition = @max(
            inputs.minimum_competition_fraction,
            if (inputs.previous_total_band_ammonium_demand_g_n[unit] >
                inputs.negligible_demand_g_n)
                inputs.previous_band_potential_oxidation_g_n[unit] /
                    inputs.previous_total_band_ammonium_demand_g_n[unit]
            else
                inputs.band_ammonium_solution_fraction[unit] *
                    inputs.fraction_of_total_active_biomass[unit],
        );
        var next_inhibitor =
            inputs.current_nitrification_inhibitor_activity[unit];
        var non_band_inhibitor_response: f64 = 1;
        var band_inhibitor_response: f64 = 1;
        if (inputs.initial_nitrification_inhibitor_activity[unit] >
            inputs.negligible_inhibitor_activity)
        {
            next_inhibitor *=
                1 -
                inputs.inhibitor_decay_rate_per_h[unit] *
                    inputs.growth_temperature_response[unit] *
                    inputs.timestep_h;
            non_band_inhibitor_response =
                inputs.initial_nitrification_inhibitor_activity[unit] -
                next_inhibitor /
                    (1 +
                        inputs.non_band_ammonium_concentration_g_n_per_m3[unit] /
                            inputs.ammonium_inhibition_of_decay_g_n_per_m3[unit]);
            band_inhibitor_response =
                inputs.initial_nitrification_inhibitor_activity[unit] -
                next_inhibitor /
                    (1 +
                        inputs.band_ammonium_concentration_g_n_per_m3[unit] /
                            inputs.ammonium_inhibition_of_decay_g_n_per_m3[unit]);
        }
        if (next_inhibitor < 0 or non_band_inhibitor_response < 0 or
            band_inhibitor_response < 0)
            return error.InvalidNitrificationInhibitorProgression;

        const unlimited =
            inputs.specific_ammonia_oxidation_g_n_per_g_c_h[unit] *
            inputs.temperature_water_response[unit] *
            inputs.combined_nutrient_limitation_fraction[unit] *
            inputs.aqueous_carbon_dioxide_limitation_fraction[unit] *
            inputs.active_biomass_g_c[unit] *
            inputs.timestep_h;
        const non_band_inhibited =
            inputs.non_band_ammonium_solution_fraction[unit] * unlimited /
            (1 +
                inputs.non_band_ammonia_concentration_g_n_per_m3[unit] /
                    inputs.ammonia_inhibition_of_oxidation_g_n_per_m3[unit]);
        const band_inhibited =
            inputs.band_ammonium_solution_fraction[unit] * unlimited /
            (1 +
                inputs.band_ammonia_concentration_g_n_per_m3[unit] /
                    inputs.ammonia_inhibition_of_oxidation_g_n_per_m3[unit]);
        const non_band_monod =
            inputs.non_band_ammonium_concentration_g_n_per_m3[unit] /
            (inputs.non_band_ammonium_concentration_g_n_per_m3[unit] +
                inputs.ammonium_half_saturation_g_n_per_m3[unit]);
        const band_monod =
            inputs.band_ammonium_concentration_g_n_per_m3[unit] /
            (inputs.band_ammonium_concentration_g_n_per_m3[unit] +
                inputs.ammonium_half_saturation_g_n_per_m3[unit]);
        const non_band_potential = non_band_inhibited * non_band_monod;
        const band_potential = band_inhibited * band_monod;
        const non_band_oxidation = @max(
            0,
            @min(
                non_band_potential,
                non_band_competition *
                    inputs.non_band_ammonium_g_n[unit] *
                    inputs.timestep_h,
            ),
        ) * non_band_inhibitor_response;
        const band_oxidation = @max(
            0,
            @min(
                band_potential,
                band_competition *
                    inputs.band_ammonium_g_n[unit] *
                    inputs.timestep_h,
            ),
        ) * band_inhibitor_response;
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
        state.next_nitrification_inhibitor_activity[unit] = next_inhibitor;
        state.non_band_inhibitor_response[unit] =
            non_band_inhibitor_response;
        state.band_inhibitor_response[unit] = band_inhibitor_response;
        state.ammonium_unlimited_oxidation_g_n[unit] = unlimited;
        state.non_band_ammonia_inhibited_oxidation_g_n[unit] =
            non_band_inhibited;
        state.band_ammonia_inhibited_oxidation_g_n[unit] = band_inhibited;
        state.non_band_ammonium_monod_fraction[unit] = non_band_monod;
        state.band_ammonium_monod_fraction[unit] = band_monod;
        state.non_band_potential_oxidation_g_n[unit] = non_band_potential;
        state.band_potential_oxidation_g_n[unit] = band_potential;
        state.non_band_ammonia_oxidation_g_n[unit] = non_band_oxidation;
        state.band_ammonia_oxidation_g_n[unit] = band_oxidation;
        state.total_ammonia_oxidation_g_n[unit] = total_oxidation;
        state.carbon_respiration_g_c[unit] = carbon_respiration;
        state.respiration_oxygen_demand_g_o[unit] = respiration_oxygen;
        state.total_oxygen_demand_g_o[unit] = total_oxygen;
    }
}

fn validateResult(state: *const State) !void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) for (@field(state.*, field.name)) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteAmmoniaOxidationResult;
}

pub fn sourceEnabled(
    zero_based_complex: usize,
    zero_based_population: usize,
) bool {
    return zero_based_complex == 5 and zero_based_population == 0;
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (n == 0 or inputs.enabled.len != n)
        return error.InvalidAmmoniaOxidationDimensions;
    inline for (.{
        inputs.previous_total_non_band_ammonium_demand_g_n,
        inputs.previous_total_band_ammonium_demand_g_n,
        inputs.previous_non_band_potential_oxidation_g_n,
        inputs.previous_band_potential_oxidation_g_n,
        inputs.fraction_of_total_active_biomass,
        inputs.non_band_ammonium_solution_fraction,
        inputs.band_ammonium_solution_fraction,
        inputs.initial_nitrification_inhibitor_activity,
        inputs.current_nitrification_inhibitor_activity,
        inputs.non_band_ammonium_concentration_g_n_per_m3,
        inputs.band_ammonium_concentration_g_n_per_m3,
        inputs.non_band_ammonia_concentration_g_n_per_m3,
        inputs.band_ammonia_concentration_g_n_per_m3,
        inputs.non_band_ammonium_g_n,
        inputs.band_ammonium_g_n,
        inputs.temperature_water_response,
        inputs.growth_temperature_response,
        inputs.combined_nutrient_limitation_fraction,
        inputs.aqueous_carbon_dioxide_limitation_fraction,
        inputs.active_biomass_g_c,
        inputs.specific_ammonia_oxidation_g_n_per_g_c_h,
        inputs.inhibitor_decay_rate_per_h,
        inputs.ammonium_inhibition_of_decay_g_n_per_m3,
        inputs.ammonia_inhibition_of_oxidation_g_n_per_m3,
        inputs.ammonium_half_saturation_g_n_per_m3,
        inputs.respiration_requirement_g_c_per_g_n,
        inputs.carbon_conversion_efficiency_g_c_per_g_n,
    }) |values| {
        if (values.len != n) return error.InvalidAmmoniaOxidationDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidAmmoniaOxidationInput;
    }
    for (0..n) |unit| if (inputs.enabled[unit]) {
        inline for (.{
            inputs.fraction_of_total_active_biomass[unit],
            inputs.non_band_ammonium_solution_fraction[unit],
            inputs.band_ammonium_solution_fraction[unit],
            inputs.combined_nutrient_limitation_fraction[unit],
            inputs.aqueous_carbon_dioxide_limitation_fraction[unit],
        }) |fraction| if (fraction > 1)
            return error.InvalidAmmoniaOxidationInput;
        if (inputs.ammonium_inhibition_of_decay_g_n_per_m3[unit] <= 0 or
            inputs.ammonia_inhibition_of_oxidation_g_n_per_m3[unit] <= 0 or
            inputs.ammonium_half_saturation_g_n_per_m3[unit] <= 0 or
            inputs.respiration_requirement_g_c_per_g_n[unit] <= 0)
            return error.InvalidAmmoniaOxidationInput;
        const decay =
            inputs.inhibitor_decay_rate_per_h[unit] *
            inputs.growth_temperature_response[unit] *
            inputs.timestep_h;
        if (decay > 1) return error.InvalidNitrificationInhibitorProgression;
    };
    inline for (.{
        inputs.oxygen_per_respired_carbon_g_o_per_g_c,
        inputs.oxygen_per_oxidized_nitrogen_g_o_per_g_n,
        inputs.timestep_h,
    }) |value| if (!std.math.isFinite(value) or value <= 0)
        return error.InvalidAmmoniaOxidationInput;
    inline for (.{
        inputs.minimum_competition_fraction,
        inputs.negligible_demand_g_n,
        inputs.negligible_inhibitor_activity,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidAmmoniaOxidationInput;
    if (inputs.minimum_competition_fraction > 1)
        return error.InvalidAmmoniaOxidationInput;
}

fn clear(state: *State) void {
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (field.type == []f64) @memset(@field(state, field.name), 0);
}

test "source ammonia oxidizer role reproduces NITRO selector" {
    try std.testing.expect(sourceEnabled(5, 0));
    try std.testing.expect(!sourceEnabled(5, 1));
    try std.testing.expect(!sourceEnabled(4, 0));
}

test "inhibition competition oxidation and oxygen demand reproduce NITRO" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, .{
        .enabled = &.{true},
        .previous_total_non_band_ammonium_demand_g_n = &.{4},
        .previous_total_band_ammonium_demand_g_n = &.{0},
        .previous_non_band_potential_oxidation_g_n = &.{1},
        .previous_band_potential_oxidation_g_n = &.{0},
        .fraction_of_total_active_biomass = &.{0.2},
        .non_band_ammonium_solution_fraction = &.{0.75},
        .band_ammonium_solution_fraction = &.{0.25},
        .initial_nitrification_inhibitor_activity = &.{1},
        .current_nitrification_inhibitor_activity = &.{0.8},
        .non_band_ammonium_concentration_g_n_per_m3 = &.{3},
        .band_ammonium_concentration_g_n_per_m3 = &.{1},
        .non_band_ammonia_concentration_g_n_per_m3 = &.{1},
        .band_ammonia_concentration_g_n_per_m3 = &.{0},
        .non_band_ammonium_g_n = &.{10},
        .band_ammonium_g_n = &.{10},
        .temperature_water_response = &.{0.5},
        .growth_temperature_response = &.{2},
        .combined_nutrient_limitation_fraction = &.{0.8},
        .aqueous_carbon_dioxide_limitation_fraction = &.{0.5},
        .active_biomass_g_c = &.{10},
        .specific_ammonia_oxidation_g_n_per_g_c_h = &.{0.2},
        .inhibitor_decay_rate_per_h = &.{0.1},
        .ammonium_inhibition_of_decay_g_n_per_m3 = &.{1},
        .ammonia_inhibition_of_oxidation_g_n_per_m3 = &.{1},
        .ammonium_half_saturation_g_n_per_m3 = &.{1},
        .respiration_requirement_g_c_per_g_n = &.{0.4},
        .carbon_conversion_efficiency_g_c_per_g_n = &.{0.2},
        .oxygen_per_respired_carbon_g_o_per_g_c = 2.667,
        .oxygen_per_oxidized_nitrogen_g_o_per_g_n = 3.429,
        .minimum_competition_fraction = 0.01,
        .negligible_demand_g_n = 1e-12,
        .negligible_inhibitor_activity = 1e-12,
        .timestep_h = 1,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), state.non_band_competition_fraction[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), state.band_competition_fraction[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.64), state.next_nitrification_inhibitor_activity[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.84), state.non_band_inhibitor_response[0], 1e-15);
    // VMXX=.4; VMXA=.75*.4/(1+1)=.15; FCN3S=.75; VMX4S=.1125.
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), state.ammonium_unlimited_oxidation_g_n[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1125), state.non_band_potential_oxidation_g_n[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0945), state.non_band_ammonia_oxidation_g_n[0], 1e-15);
    try std.testing.expect(state.total_oxygen_demand_g_o[0] > state.respiration_oxygen_demand_g_o[0]);
}

test "invalid late inhibitor step leaves prior state unchanged" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.total_ammonia_oxidation_g_n[0] = 7;
    const zeros = [_]f64{ 0, 0 };
    const ones = [_]f64{ 1, 1 };
    try std.testing.expectError(error.InvalidNitrificationInhibitorProgression, calculate(&state, .{
        .enabled = &.{ true, true },
        .previous_total_non_band_ammonium_demand_g_n = &zeros,
        .previous_total_band_ammonium_demand_g_n = &zeros,
        .previous_non_band_potential_oxidation_g_n = &zeros,
        .previous_band_potential_oxidation_g_n = &zeros,
        .fraction_of_total_active_biomass = &ones,
        .non_band_ammonium_solution_fraction = &ones,
        .band_ammonium_solution_fraction = &zeros,
        .initial_nitrification_inhibitor_activity = &ones,
        .current_nitrification_inhibitor_activity = &ones,
        .non_band_ammonium_concentration_g_n_per_m3 = &ones,
        .band_ammonium_concentration_g_n_per_m3 = &ones,
        .non_band_ammonia_concentration_g_n_per_m3 = &zeros,
        .band_ammonia_concentration_g_n_per_m3 = &zeros,
        .non_band_ammonium_g_n = &ones,
        .band_ammonium_g_n = &zeros,
        .temperature_water_response = &ones,
        .growth_temperature_response = &ones,
        .combined_nutrient_limitation_fraction = &ones,
        .aqueous_carbon_dioxide_limitation_fraction = &ones,
        .active_biomass_g_c = &ones,
        .specific_ammonia_oxidation_g_n_per_g_c_h = &ones,
        .inhibitor_decay_rate_per_h = &.{ 0.1, 2 },
        .ammonium_inhibition_of_decay_g_n_per_m3 = &ones,
        .ammonia_inhibition_of_oxidation_g_n_per_m3 = &ones,
        .ammonium_half_saturation_g_n_per_m3 = &ones,
        .respiration_requirement_g_c_per_g_n = &ones,
        .carbon_conversion_efficiency_g_c_per_g_n = &ones,
        .oxygen_per_respired_carbon_g_o_per_g_c = 2.667,
        .oxygen_per_oxidized_nitrogen_g_o_per_g_n = 3.429,
        .minimum_competition_fraction = 0.01,
        .negligible_demand_g_n = 1e-12,
        .negligible_inhibitor_activity = 1e-12,
        .timestep_h = 1,
    }));
    try std.testing.expectEqual(@as(f64, 7), state.total_ammonia_oxidation_g_n[0]);
}

test "NITRO 1072-1193 derived overflow leaves ammonia oxidation unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.total_ammonia_oxidation_g_n[0] = 7;
    const zero = [_]f64{0};
    const one = [_]f64{1};
    try std.testing.expectError(
        error.NonFiniteAmmoniaOxidationResult,
        calculate(&state, .{
            .enabled = &.{true},
            .previous_total_non_band_ammonium_demand_g_n = &zero,
            .previous_total_band_ammonium_demand_g_n = &zero,
            .previous_non_band_potential_oxidation_g_n = &zero,
            .previous_band_potential_oxidation_g_n = &zero,
            .fraction_of_total_active_biomass = &one,
            .non_band_ammonium_solution_fraction = &one,
            .band_ammonium_solution_fraction = &zero,
            .initial_nitrification_inhibitor_activity = &zero,
            .current_nitrification_inhibitor_activity = &zero,
            .non_band_ammonium_concentration_g_n_per_m3 = &one,
            .band_ammonium_concentration_g_n_per_m3 = &zero,
            .non_band_ammonia_concentration_g_n_per_m3 = &zero,
            .band_ammonia_concentration_g_n_per_m3 = &zero,
            .non_band_ammonium_g_n = &one,
            .band_ammonium_g_n = &zero,
            .temperature_water_response = &one,
            .growth_temperature_response = &one,
            .combined_nutrient_limitation_fraction = &one,
            .aqueous_carbon_dioxide_limitation_fraction = &one,
            .active_biomass_g_c = &.{std.math.floatMax(f64)},
            .specific_ammonia_oxidation_g_n_per_g_c_h = &.{std.math.floatMax(f64)},
            .inhibitor_decay_rate_per_h = &zero,
            .ammonium_inhibition_of_decay_g_n_per_m3 = &one,
            .ammonia_inhibition_of_oxidation_g_n_per_m3 = &one,
            .ammonium_half_saturation_g_n_per_m3 = &one,
            .respiration_requirement_g_c_per_g_n = &one,
            .carbon_conversion_efficiency_g_c_per_g_n = &one,
            .oxygen_per_respired_carbon_g_o_per_g_c = 2.667,
            .oxygen_per_oxidized_nitrogen_g_o_per_g_n = 3.429,
            .minimum_competition_fraction = 0.01,
            .negligible_demand_g_n = 1e-12,
            .negligible_inhibitor_activity = 1e-12,
            .timestep_h = 1,
        }),
    );
    try std.testing.expectEqual(@as(f64, 7), state.total_ammonia_oxidation_g_n[0]);
}
