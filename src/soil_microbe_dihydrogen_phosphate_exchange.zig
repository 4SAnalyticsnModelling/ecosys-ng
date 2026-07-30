const std = @import("std");

pub const Inputs = struct {
    nonstructural_carbon_g_c: []const f64,
    nonstructural_phosphorus_g_p: []const f64,
    maximum_phosphorus_to_carbon_ratio_g_p_per_g_c: []const f64,
    microbial_surface_area_m2_per_g_c: []const f64,
    active_biomass_g_c: []const f64,
    biological_activity_fraction: []const f64,
    non_band_dihydrogen_phosphate_fraction: []const f64,
    band_dihydrogen_phosphate_fraction: []const f64,
    non_band_competition: []const f64,
    band_competition: []const f64,
    non_band_concentration_g_p_m3: []const f64,
    band_concentration_g_p_m3: []const f64,
    non_band_dihydrogen_phosphate_g_p: []const f64,
    band_dihydrogen_phosphate_g_p: []const f64,
    soil_water_m3: []const f64,
    maximum_uptake_g_p_m2_h: f64,
    protected_concentration_g_p_m3: f64,
    uptake_half_saturation_g_p_m3: f64,
    biochemical_time_fraction_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    total_net_dihydrogen_phosphate_exchange_g_p: f64,
    stoichiometric_demand_g_p: []f64,
    rate_limited_demand_g_p: []f64,
    non_band_concentration_excess_g_p_m3: []f64,
    band_concentration_excess_g_p_m3: []f64,
    non_band_unlimited_exchange_g_p: []f64,
    band_unlimited_exchange_g_p: []f64,
    non_band_protected_phosphate_g_p: []f64,
    band_protected_phosphate_g_p: []f64,
    non_band_available_phosphate_g_p: []f64,
    band_available_phosphate_g_p: []f64,
    non_band_exchange_g_p: []f64,
    band_exchange_g_p: []f64,
    total_exchange_g_p: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidDihydrogenPhosphateExchangeDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
        state.total_net_dihydrogen_phosphate_exchange_g_p = 0;
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

/// Exact NITRO.F 2184--2232 microbial H2PO4 mineralization/immobilization.
/// Negative exchange mineralizes phosphorus; positive exchange immobilizes it.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([13]f64, state.unit_count);
    defer state.allocator.free(temporary);
    var total_net_exchange: f64 = 0;
    for (0..state.unit_count) |unit| {
        const carbon = @max(0, inputs.nonstructural_carbon_g_c[unit]);
        const phosphorus = @max(0, inputs.nonstructural_phosphorus_g_p[unit]);
        const demand = carbon *
            inputs.maximum_phosphorus_to_carbon_ratio_g_p_per_g_c[unit] - phosphorus;
        var rate_limited: f64 = 0;
        var non_band_excess: f64 = 0;
        var band_excess: f64 = 0;
        var non_band_unlimited: f64 = 0;
        var band_unlimited: f64 = 0;
        var non_band_protected: f64 = 0;
        var band_protected: f64 = 0;
        var non_band_available: f64 = 0;
        var band_available: f64 = 0;
        var non_band_exchange: f64 = 0;
        var band_exchange: f64 = 0;
        if (demand > 0) {
            non_band_excess = @max(0, inputs.non_band_concentration_g_p_m3[unit] -
                inputs.protected_concentration_g_p_m3);
            band_excess = @max(0, inputs.band_concentration_g_p_m3[unit] -
                inputs.protected_concentration_g_p_m3);
            rate_limited = @min(demand, inputs.microbial_surface_area_m2_per_g_c[unit] *
                inputs.active_biomass_g_c[unit] *
                inputs.biological_activity_fraction[unit] *
                inputs.maximum_uptake_g_p_m2_h *
                inputs.biochemical_time_fraction_h);
            non_band_unlimited =
                inputs.non_band_dihydrogen_phosphate_fraction[unit] * rate_limited *
                non_band_excess / (non_band_excess + inputs.uptake_half_saturation_g_p_m3);
            band_unlimited = inputs.band_dihydrogen_phosphate_fraction[unit] * rate_limited *
                band_excess / (band_excess + inputs.uptake_half_saturation_g_p_m3);
            non_band_protected = inputs.protected_concentration_g_p_m3 *
                inputs.soil_water_m3[unit] *
                inputs.non_band_dihydrogen_phosphate_fraction[unit];
            band_protected = inputs.protected_concentration_g_p_m3 *
                inputs.soil_water_m3[unit] *
                inputs.band_dihydrogen_phosphate_fraction[unit];
            non_band_available = inputs.non_band_competition[unit] *
                @max(0, (inputs.non_band_dihydrogen_phosphate_g_p[unit] -
                    non_band_protected) *
                    inputs.biochemical_time_fraction_h);
            band_available = inputs.band_competition[unit] *
                @max(0, (inputs.band_dihydrogen_phosphate_g_p[unit] -
                    band_protected) *
                    inputs.biochemical_time_fraction_h);
            non_band_exchange = @min(non_band_available, non_band_unlimited);
            band_exchange = @min(band_available, band_unlimited);
        } else {
            non_band_exchange =
                demand * inputs.non_band_dihydrogen_phosphate_fraction[unit];
            band_exchange = demand * inputs.band_dihydrogen_phosphate_fraction[unit];
        }
        const total_exchange = non_band_exchange + band_exchange;
        temporary[unit] = .{
            demand,             rate_limited,   non_band_excess,    band_excess,
            non_band_unlimited, band_unlimited, non_band_protected, band_protected,
            non_band_available, band_available, non_band_exchange,  band_exchange,
            total_exchange,
        };
        for (temporary[unit], 0..) |value, index| {
            if (!std.math.isFinite(value))
                return error.NonFiniteDihydrogenPhosphateExchangeResult;
            if (index != 0 and index != 10 and index != 11 and index != 12 and
                value < 0)
                return error.NonFiniteDihydrogenPhosphateExchangeResult;
        }
        total_net_exchange += total_exchange;
    }
    if (!std.math.isFinite(total_net_exchange))
        return error.NonFiniteDihydrogenPhosphateExchangeResult;
    state.total_net_dihydrogen_phosphate_exchange_g_p = total_net_exchange;
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
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != n) return error.InvalidDihydrogenPhosphateExchangeDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidDihydrogenPhosphateExchangeInput;
        const clamped = comptime std.mem.eql(u8, field.name, "nonstructural_carbon_g_c") or
            std.mem.eql(u8, field.name, "nonstructural_phosphorus_g_p");
        if (!clamped) for (values) |value| if (value < 0)
            return error.InvalidDihydrogenPhosphateExchangeInput;
    };
    inline for (.{
        inputs.maximum_uptake_g_p_m2_h,
        inputs.protected_concentration_g_p_m3,
        inputs.uptake_half_saturation_g_p_m3,
        inputs.biochemical_time_fraction_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidDihydrogenPhosphateExchangeInput;
    if (inputs.uptake_half_saturation_g_p_m3 == 0)
        return error.InvalidDihydrogenPhosphateExchangeInput;
}

fn fixture() Inputs {
    return .{
        .nonstructural_carbon_g_c = &.{10},
        .nonstructural_phosphorus_g_p = &.{1},
        .maximum_phosphorus_to_carbon_ratio_g_p_per_g_c = &.{0.2},
        .microbial_surface_area_m2_per_g_c = &.{2},
        .active_biomass_g_c = &.{4},
        .biological_activity_fraction = &.{0.5},
        .non_band_dihydrogen_phosphate_fraction = &.{0.6},
        .band_dihydrogen_phosphate_fraction = &.{0.4},
        .non_band_competition = &.{0.5},
        .band_competition = &.{0.5},
        .non_band_concentration_g_p_m3 = &.{4},
        .band_concentration_g_p_m3 = &.{2},
        .non_band_dihydrogen_phosphate_g_p = &.{20},
        .band_dihydrogen_phosphate_g_p = &.{10},
        .soil_water_m3 = &.{2},
        .maximum_uptake_g_p_m2_h = 1,
        .protected_concentration_g_p_m3 = 1,
        .uptake_half_saturation_g_p_m3 = 2,
        .biochemical_time_fraction_h = 1,
    };
}

test "phosphorus deficit immobilizes dihydrogen phosphate in two zones" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(1, state.stoichiometric_demand_g_p[0], 1e-12);
    try std.testing.expect(state.non_band_exchange_g_p[0] > 0);
    try std.testing.expect(state.band_exchange_g_p[0] > 0);
    try std.testing.expectEqual(
        state.total_exchange_g_p[0],
        state.total_net_dihydrogen_phosphate_exchange_g_p,
    );
}

test "phosphorus excess mineralizes by zone fractions" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.nonstructural_phosphorus_g_p = &.{3};
    try calculate(&state, inputs);
    try std.testing.expectApproxEqAbs(-0.6, state.non_band_exchange_g_p[0], 1e-12);
    try std.testing.expectApproxEqAbs(-0.4, state.band_exchange_g_p[0], 1e-12);
}

test "protected phosphate cannot be immobilized" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.non_band_dihydrogen_phosphate_g_p = &.{0.5};
    inputs.band_dihydrogen_phosphate_g_p = &.{0.5};
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.total_exchange_g_p[0]);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.total_exchange_g_p[0] = 7;
    var inputs = fixture();
    inputs.soil_water_m3 = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidDihydrogenPhosphateExchangeInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.total_exchange_g_p[0]);
}

test "NITRO 2184-2232 derived overflow preserves phosphate exchange" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.total_exchange_g_p[0] = 7;
    var inputs = fixture();
    inputs.nonstructural_carbon_g_c = &.{std.math.floatMax(f64)};
    inputs.maximum_phosphorus_to_carbon_ratio_g_p_per_g_c =
        &.{std.math.floatMax(f64)};
    try std.testing.expectError(
        error.NonFiniteDihydrogenPhosphateExchangeResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(@as(f64, 7), state.total_exchange_g_p[0]);
}
