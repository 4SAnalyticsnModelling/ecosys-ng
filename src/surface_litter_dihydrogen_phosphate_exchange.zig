const std = @import("std");

pub const Inputs = struct {
    preceding_total_phosphate_exchange_g_p: f64,
    surface_litter_unit: []const bool,
    stoichiometric_phosphorus_demand_g_p: []const f64,
    upstream_dihydrogen_phosphate_exchange_g_p: []const f64,
    microbial_surface_area_m2_per_g_c: []const f64,
    active_biomass_g_c: []const f64,
    biological_activity_fraction: []const f64,
    non_band_phosphate_fraction: []const f64,
    band_phosphate_fraction: []const f64,
    surface_non_band_concentration_g_p_m3: []const f64,
    surface_band_concentration_g_p_m3: []const f64,
    surface_total_dihydrogen_phosphate_g_p: []const f64,
    surface_soil_water_m3: []const f64,
    litter_phosphate_competition: []const f64,
    maximum_uptake_g_p_m2_h: f64,
    protected_concentration_g_p_m3: f64,
    uptake_half_saturation_g_p_m3: f64,
    biochemical_time_fraction_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    total_net_phosphate_exchange_g_p: f64,
    residual_litter_demand_g_p: []f64,
    non_band_concentration_excess_g_p_m3: []f64,
    band_concentration_excess_g_p_m3: []f64,
    rate_limited_demand_g_p: []f64,
    concentration_response: []f64,
    unlimited_exchange_g_p: []f64,
    protected_phosphate_g_p: []f64,
    available_phosphate_g_p: []f64,
    exchange_g_p: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidSurfaceLitterDihydrogenPhosphateDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
        state.total_net_phosphate_exchange_g_p = 0;
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

/// Exact NITRO.F 2390--2430 additional surface-litter H2PO4 exchange.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([9]f64, state.unit_count);
    defer state.allocator.free(temporary);
    var total_net_exchange = inputs.preceding_total_phosphate_exchange_g_p;
    for (0..state.unit_count) |unit| {
        if (!inputs.surface_litter_unit[unit]) {
            temporary[unit] = @splat(0);
            continue;
        }
        const residual = inputs.stoichiometric_phosphorus_demand_g_p[unit] -
            inputs.upstream_dihydrogen_phosphate_exchange_g_p[unit];
        var non_band_excess: f64 = 0;
        var band_excess: f64 = 0;
        var rate_limited: f64 = 0;
        var response: f64 = 0;
        var unlimited: f64 = 0;
        var protected: f64 = 0;
        var available: f64 = 0;
        var exchange = residual;
        if (residual > 0) {
            non_band_excess = @max(0, inputs.surface_non_band_concentration_g_p_m3[unit] -
                inputs.protected_concentration_g_p_m3);
            band_excess = @max(0, inputs.surface_band_concentration_g_p_m3[unit] -
                inputs.protected_concentration_g_p_m3);
            rate_limited = @min(residual, inputs.microbial_surface_area_m2_per_g_c[unit] *
                inputs.active_biomass_g_c[unit] *
                inputs.biological_activity_fraction[unit] *
                inputs.maximum_uptake_g_p_m2_h *
                inputs.biochemical_time_fraction_h);
            response = inputs.non_band_phosphate_fraction[unit] * non_band_excess /
                (non_band_excess + inputs.uptake_half_saturation_g_p_m3) +
                inputs.band_phosphate_fraction[unit] * band_excess /
                    (band_excess + inputs.uptake_half_saturation_g_p_m3);
            unlimited = rate_limited * response;
            protected = inputs.protected_concentration_g_p_m3 *
                inputs.surface_soil_water_m3[unit];
            available = inputs.litter_phosphate_competition[unit] *
                @max(0, inputs.surface_total_dihydrogen_phosphate_g_p[unit] - protected);
            exchange = @min(available, unlimited);
        }
        temporary[unit] = .{
            residual,  non_band_excess, band_excess, rate_limited, response,
            unlimited, protected,       available,   exchange,
        };
        for (temporary[unit], 0..) |value, index| {
            const signed = index == 0 or index == 8;
            if (!std.math.isFinite(value) or (!signed and value < 0))
                return error.NonFiniteSurfaceLitterDihydrogenPhosphateResult;
        }
        total_net_exchange += exchange;
    }
    if (!std.math.isFinite(total_net_exchange))
        return error.NonFiniteSurfaceLitterDihydrogenPhosphateResult;
    state.total_net_phosphate_exchange_g_p = total_net_exchange;
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
    if (!std.math.isFinite(inputs.preceding_total_phosphate_exchange_g_p))
        return error.InvalidSurfaceLitterDihydrogenPhosphateInput;
    if (inputs.surface_litter_unit.len != n)
        return error.InvalidSurfaceLitterDihydrogenPhosphateDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != n) return error.InvalidSurfaceLitterDihydrogenPhosphateDimensions;
        for (values) |value| if (!std.math.isFinite(value))
            return error.InvalidSurfaceLitterDihydrogenPhosphateInput;
        const signed = comptime std.mem.eql(u8, field.name, "stoichiometric_phosphorus_demand_g_p") or
            std.mem.eql(u8, field.name, "upstream_dihydrogen_phosphate_exchange_g_p");
        if (!signed) for (values) |value| if (value < 0)
            return error.InvalidSurfaceLitterDihydrogenPhosphateInput;
    };
    inline for (.{
        inputs.maximum_uptake_g_p_m2_h,
        inputs.protected_concentration_g_p_m3,
        inputs.uptake_half_saturation_g_p_m3,
        inputs.biochemical_time_fraction_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidSurfaceLitterDihydrogenPhosphateInput;
    if (inputs.uptake_half_saturation_g_p_m3 == 0)
        return error.InvalidSurfaceLitterDihydrogenPhosphateInput;
}

fn fixture() Inputs {
    return .{
        .preceding_total_phosphate_exchange_g_p = 0.6,
        .surface_litter_unit = &.{true},
        .stoichiometric_phosphorus_demand_g_p = &.{2},
        .upstream_dihydrogen_phosphate_exchange_g_p = &.{0.6},
        .microbial_surface_area_m2_per_g_c = &.{2},
        .active_biomass_g_c = &.{4},
        .biological_activity_fraction = &.{0.5},
        .non_band_phosphate_fraction = &.{0.6},
        .band_phosphate_fraction = &.{0.4},
        .surface_non_band_concentration_g_p_m3 = &.{4},
        .surface_band_concentration_g_p_m3 = &.{2},
        .surface_total_dihydrogen_phosphate_g_p = &.{20},
        .surface_soil_water_m3 = &.{2},
        .litter_phosphate_competition = &.{0.5},
        .maximum_uptake_g_p_m2_h = 1,
        .protected_concentration_g_p_m3 = 1,
        .uptake_half_saturation_g_p_m3 = 2,
        .biochemical_time_fraction_h = 1,
    };
}

test "surface litter immobilizes residual dihydrogen phosphate demand" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(1.4, state.residual_litter_demand_g_p[0], 1e-12);
    try std.testing.expect(state.exchange_g_p[0] > 0);
    try std.testing.expectApproxEqAbs(
        0.6 + state.exchange_g_p[0],
        state.total_net_phosphate_exchange_g_p,
        1e-12,
    );
}

test "negative litter phosphorus residual mineralizes unrestricted" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.stoichiometric_phosphorus_demand_g_p = &.{0};
    try calculate(&state, inputs);
    try std.testing.expectApproxEqAbs(-0.6, state.exchange_g_p[0], 1e-12);
}

test "non-surface units do not execute litter phosphate exchange" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.surface_litter_unit = &.{false};
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.exchange_g_p[0]);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.exchange_g_p[0] = 7;
    var inputs = fixture();
    inputs.surface_soil_water_m3 = &.{std.math.nan(f64)};
    try std.testing.expectError(
        error.InvalidSurfaceLitterDihydrogenPhosphateInput,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.exchange_g_p[0]);
}

test "NITRO 2390-2430 derived overflow preserves litter phosphate state" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.exchange_g_p[0] = 7;
    state.total_net_phosphate_exchange_g_p = 8;
    var inputs = fixture();
    inputs.stoichiometric_phosphorus_demand_g_p = &.{std.math.floatMax(f64)};
    inputs.non_band_phosphate_fraction = &.{std.math.floatMax(f64)};
    inputs.band_phosphate_fraction = &.{std.math.floatMax(f64)};
    try std.testing.expectError(
        error.NonFiniteSurfaceLitterDihydrogenPhosphateResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.exchange_g_p[0]);
    try std.testing.expectEqual(8, state.total_net_phosphate_exchange_g_p);
}
