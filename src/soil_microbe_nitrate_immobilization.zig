const std = @import("std");

pub const Inputs = struct {
    preceding_total_mineral_exchange_g_n: f64,
    ammonium_stoichiometric_demand_g_n: []const f64,
    non_band_ammonium_exchange_g_n: []const f64,
    band_ammonium_exchange_g_n: []const f64,
    microbial_surface_area_m2_per_g_c: []const f64,
    active_biomass_g_c: []const f64,
    biological_activity_fraction: []const f64,
    non_band_nitrate_fraction: []const f64,
    band_nitrate_fraction: []const f64,
    non_band_nitrate_competition: []const f64,
    band_nitrate_competition: []const f64,
    non_band_nitrate_concentration_g_n_m3: []const f64,
    band_nitrate_concentration_g_n_m3: []const f64,
    non_band_nitrate_g_n: []const f64,
    band_nitrate_g_n: []const f64,
    soil_water_m3: []const f64,
    maximum_uptake_g_n_m2_h: f64,
    protected_concentration_g_n_m3: f64,
    uptake_half_saturation_g_n_m3: f64,
    biochemical_time_fraction_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    total_net_mineral_exchange_g_n: f64,
    residual_demand_g_n: []f64,
    rate_limited_demand_g_n: []f64,
    non_band_concentration_excess_g_n_m3: []f64,
    band_concentration_excess_g_n_m3: []f64,
    non_band_unlimited_immobilization_g_n: []f64,
    band_unlimited_immobilization_g_n: []f64,
    non_band_protected_nitrate_g_n: []f64,
    band_protected_nitrate_g_n: []f64,
    non_band_available_nitrate_g_n: []f64,
    band_available_nitrate_g_n: []f64,
    non_band_immobilization_g_n: []f64,
    band_immobilization_g_n: []f64,
    total_immobilization_g_n: []f64,
    total_mineral_exchange_g_n: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidNitrateImmobilizationDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
        state.total_net_mineral_exchange_g_n = 0;
        var allocated: usize = 0;
        errdefer {
            inline for (@typeInfo(State).@"struct".fields) |field| {
                if (field.type == []f64 and allocated > 0) {
                    allocated -= 1;
                    allocator.free(@field(state, field.name));
                }
            }
        }
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            @field(state, field.name) = try allocator.alloc(f64, unit_count);
            @memset(@field(state, field.name), 0);
            allocated += 1;
        };
        return state;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field|
            if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

/// Exact NITRO.F 2128--2173 residual microbial NO3 immobilization.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([14]f64, state.unit_count);
    defer state.allocator.free(temporary);
    var total_net_exchange = inputs.preceding_total_mineral_exchange_g_n;
    for (0..state.unit_count) |unit| {
        const residual = @max(0, inputs.ammonium_stoichiometric_demand_g_n[unit] -
            inputs.non_band_ammonium_exchange_g_n[unit] -
            inputs.band_ammonium_exchange_g_n[unit]);
        var rate_limited: f64 = 0;
        var non_band_excess: f64 = 0;
        var band_excess: f64 = 0;
        var non_band_unlimited: f64 = 0;
        var band_unlimited: f64 = 0;
        var non_band_protected: f64 = 0;
        var band_protected: f64 = 0;
        var non_band_available: f64 = 0;
        var band_available: f64 = 0;
        var non_band_actual: f64 = 0;
        var band_actual: f64 = 0;
        if (residual > 0) {
            non_band_excess = @max(0, inputs.non_band_nitrate_concentration_g_n_m3[unit] -
                inputs.protected_concentration_g_n_m3);
            band_excess = @max(0, inputs.band_nitrate_concentration_g_n_m3[unit] -
                inputs.protected_concentration_g_n_m3);
            rate_limited = @min(residual, inputs.microbial_surface_area_m2_per_g_c[unit] *
                inputs.active_biomass_g_c[unit] *
                inputs.biological_activity_fraction[unit] *
                inputs.maximum_uptake_g_n_m2_h *
                inputs.biochemical_time_fraction_h);
            non_band_unlimited = inputs.non_band_nitrate_fraction[unit] * rate_limited *
                non_band_excess / (non_band_excess + inputs.uptake_half_saturation_g_n_m3);
            band_unlimited = inputs.band_nitrate_fraction[unit] * rate_limited *
                band_excess / (band_excess + inputs.uptake_half_saturation_g_n_m3);
            non_band_protected = inputs.protected_concentration_g_n_m3 *
                inputs.soil_water_m3[unit] * inputs.non_band_nitrate_fraction[unit];
            band_protected = inputs.protected_concentration_g_n_m3 *
                inputs.soil_water_m3[unit] * inputs.band_nitrate_fraction[unit];
            non_band_available = inputs.non_band_nitrate_competition[unit] *
                @max(0, (inputs.non_band_nitrate_g_n[unit] -
                    non_band_protected) *
                    inputs.biochemical_time_fraction_h);
            band_available = inputs.band_nitrate_competition[unit] *
                @max(0, (inputs.band_nitrate_g_n[unit] - band_protected) *
                    inputs.biochemical_time_fraction_h);
            non_band_actual = @min(non_band_available, non_band_unlimited);
            band_actual = @min(band_available, band_unlimited);
        } else {
            // Retained literally from the source; RINOP is zero in this branch.
            non_band_actual = residual * inputs.non_band_nitrate_fraction[unit];
            band_actual = residual * inputs.band_nitrate_fraction[unit];
        }
        const nitrate_total = non_band_actual + band_actual;
        temporary[unit] = .{
            residual,           rate_limited,   non_band_excess,    band_excess,
            non_band_unlimited, band_unlimited, non_band_protected, band_protected,
            non_band_available, band_available, non_band_actual,    band_actual,
            nitrate_total,
            inputs.non_band_ammonium_exchange_g_n[unit] +
                inputs.band_ammonium_exchange_g_n[unit] + nitrate_total,
        };
        for (temporary[unit], 0..) |value, index| {
            if (!std.math.isFinite(value))
                return error.NonFiniteNitrateImmobilizationResult;
            if (index != 13 and value < 0)
                return error.NonFiniteNitrateImmobilizationResult;
        }
        total_net_exchange += nitrate_total;
    }
    if (!std.math.isFinite(total_net_exchange))
        return error.NonFiniteNitrateImmobilizationResult;
    state.total_net_mineral_exchange_g_n = total_net_exchange;
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, unit| @field(state, field.name)[unit] = values[index];
    };
}

fn stateFieldIndex(comptime name: []const u8) usize {
    var index: usize = 0;
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        if (comptime std.mem.eql(u8, field.name, name)) return index;
        index += 1;
    };
    unreachable;
}

fn validate(state: *const State, inputs: Inputs) !void {
    const n = state.unit_count;
    if (!std.math.isFinite(inputs.preceding_total_mineral_exchange_g_n))
        return error.InvalidNitrateImmobilizationInput;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != n) return error.InvalidNitrateImmobilizationDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidNitrateImmobilizationInput;
        const signed = comptime std.mem.eql(u8, field.name, "ammonium_stoichiometric_demand_g_n") or
            std.mem.eql(u8, field.name, "non_band_ammonium_exchange_g_n") or
            std.mem.eql(u8, field.name, "band_ammonium_exchange_g_n");
        if (!signed) for (values) |value| if (value < 0)
            return error.InvalidNitrateImmobilizationInput;
    };
    inline for (.{
        inputs.maximum_uptake_g_n_m2_h,
        inputs.protected_concentration_g_n_m3,
        inputs.uptake_half_saturation_g_n_m3,
        inputs.biochemical_time_fraction_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidNitrateImmobilizationInput;
    if (inputs.uptake_half_saturation_g_n_m3 == 0)
        return error.InvalidNitrateImmobilizationInput;
}

fn fixture() Inputs {
    return .{
        .preceding_total_mineral_exchange_g_n = 0.6,
        .ammonium_stoichiometric_demand_g_n = &.{2},
        .non_band_ammonium_exchange_g_n = &.{0.4},
        .band_ammonium_exchange_g_n = &.{0.2},
        .microbial_surface_area_m2_per_g_c = &.{2},
        .active_biomass_g_c = &.{4},
        .biological_activity_fraction = &.{0.5},
        .non_band_nitrate_fraction = &.{0.6},
        .band_nitrate_fraction = &.{0.4},
        .non_band_nitrate_competition = &.{0.5},
        .band_nitrate_competition = &.{0.5},
        .non_band_nitrate_concentration_g_n_m3 = &.{4},
        .band_nitrate_concentration_g_n_m3 = &.{2},
        .non_band_nitrate_g_n = &.{20},
        .band_nitrate_g_n = &.{10},
        .soil_water_m3 = &.{2},
        .maximum_uptake_g_n_m2_h = 1,
        .protected_concentration_g_n_m3 = 1,
        .uptake_half_saturation_g_n_m3 = 2,
        .biochemical_time_fraction_h = 1,
    };
}

test "nitrate immobilization consumes demand left after ammonium" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(1.4, state.residual_demand_g_n[0], 1e-12);
    try std.testing.expect(state.total_immobilization_g_n[0] > 0);
    try std.testing.expectApproxEqAbs(
        0.6 + state.total_immobilization_g_n[0],
        state.total_mineral_exchange_g_n[0],
        1e-12,
    );
    var expected_source_total: f64 = 0.6;
    expected_source_total += state.non_band_immobilization_g_n[0] +
        state.band_immobilization_g_n[0];
    try std.testing.expectEqual(expected_source_total, state.total_net_mineral_exchange_g_n);
}

test "satisfied ammonium demand disables nitrate immobilization" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.non_band_ammonium_exchange_g_n = &.{1.5};
    inputs.band_ammonium_exchange_g_n = &.{0.5};
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.total_immobilization_g_n[0]);
}

test "protected nitrate cannot be immobilized" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.non_band_nitrate_g_n = &.{0.5};
    inputs.band_nitrate_g_n = &.{0.5};
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.total_immobilization_g_n[0]);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.total_immobilization_g_n[0] = 7;
    var inputs = fixture();
    inputs.band_nitrate_g_n = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidNitrateImmobilizationInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.total_immobilization_g_n[0]);
}

test "NITRO 2128-2173 derived overflow preserves nitrate exchange" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.total_immobilization_g_n[0] = 7;
    var inputs = fixture();
    inputs.ammonium_stoichiometric_demand_g_n =
        &.{std.math.floatMax(f64)};
    inputs.microbial_surface_area_m2_per_g_c =
        &.{std.math.floatMax(f64)};
    inputs.active_biomass_g_c = &.{std.math.floatMax(f64)};
    inputs.non_band_nitrate_concentration_g_n_m3 =
        &.{std.math.floatMax(f64)};
    try std.testing.expectError(
        error.NonFiniteNitrateImmobilizationResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(@as(f64, 7), state.total_immobilization_g_n[0]);
}
