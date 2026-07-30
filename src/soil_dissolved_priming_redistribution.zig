const std = @import("std");

pub const Inputs = struct {
    soil_organic_carbon_g_c: []const f64,
    microbial_activity_respiration_g_c: []const f64,
    dissolved_organic_carbon_g_c: []const f64,
    dissolved_organic_nitrogen_g_n: []const f64,
    dissolved_organic_phosphorus_g_p: []const f64,
    acetate_g_c: []const f64,
    temperature_response: f64,
    priming_transfer_rate_h: f64,
    biochemical_time_fraction_h: f64,
    negligible_soil_organic_carbon_g_c: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    complex_count: usize,
    activity_transfer_g_c: []f64,
    dissolved_carbon_transfer_g_c: []f64,
    dissolved_nitrogen_transfer_g_n: []f64,
    dissolved_phosphorus_transfer_g_p: []f64,
    acetate_transfer_g_c: []f64,
    redistributed_activity_respiration_g_c: []f64,
    redistributed_dissolved_organic_carbon_g_c: []f64,
    redistributed_dissolved_organic_nitrogen_g_n: []f64,
    redistributed_dissolved_organic_phosphorus_g_p: []f64,
    redistributed_acetate_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, complex_count: usize) !State {
        if (complex_count == 0) return error.InvalidPrimingRedistributionDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.complex_count = complex_count;
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
            @field(state, field.name) = try allocator.alloc(f64, complex_count);
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

/// Exact dissolved-pool portion of NITRO.F 3003--3080 priming redistribution.
/// Pair order follows the source K ascending, then KK ascending.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const deltas = try state.allocator.alloc([5]f64, state.complex_count);
    defer state.allocator.free(deltas);
    @memset(deltas, @splat(0));
    for (0..state.complex_count) |donor| {
        for (donor + 1..state.complex_count) |receiver| {
            const donor_soc = inputs.soil_organic_carbon_g_c[donor];
            const receiver_soc = inputs.soil_organic_carbon_g_c[receiver];
            if (donor_soc <= inputs.negligible_soil_organic_carbon_g_c or
                receiver_soc <= inputs.negligible_soil_organic_carbon_g_c)
                continue;
            const combined_soc = donor_soc + receiver_soc;
            const fluxes: [5]f64 = .{
                try pairFlux(inputs, inputs.microbial_activity_respiration_g_c, donor_soc, receiver_soc, donor, receiver, combined_soc),
                try pairFlux(inputs, inputs.dissolved_organic_carbon_g_c, donor_soc, receiver_soc, donor, receiver, combined_soc),
                try pairFlux(inputs, inputs.dissolved_organic_nitrogen_g_n, donor_soc, receiver_soc, donor, receiver, combined_soc),
                try pairFlux(inputs, inputs.dissolved_organic_phosphorus_g_p, donor_soc, receiver_soc, donor, receiver, combined_soc),
                try pairFlux(inputs, inputs.acetate_g_c, donor_soc, receiver_soc, donor, receiver, combined_soc),
            };
            inline for (0..5) |pool| {
                const values = poolValues(inputs, pool);
                const flux = fluxes[pool];
                if (values[donor] + deltas[donor][pool] - flux > 0 and
                    values[receiver] + deltas[receiver][pool] + flux > 0)
                {
                    const donor_delta = deltas[donor][pool] - flux;
                    const receiver_delta = deltas[receiver][pool] + flux;
                    if (!std.math.isFinite(donor_delta) or
                        !std.math.isFinite(receiver_delta))
                        return error.NonFinitePrimingRedistributionResult;
                    deltas[donor][pool] = donor_delta;
                    deltas[receiver][pool] = receiver_delta;
                }
            }
        }
    }
    const temporary = try state.allocator.alloc([10]f64, state.complex_count);
    defer state.allocator.free(temporary);
    for (0..state.complex_count) |complex| {
        temporary[complex] = .{
            deltas[complex][0],
            deltas[complex][1],
            deltas[complex][2],
            deltas[complex][3],
            deltas[complex][4],
            inputs.microbial_activity_respiration_g_c[complex] + deltas[complex][0],
            inputs.dissolved_organic_carbon_g_c[complex] + deltas[complex][1],
            inputs.dissolved_organic_nitrogen_g_n[complex] + deltas[complex][2],
            inputs.dissolved_organic_phosphorus_g_p[complex] + deltas[complex][3],
            inputs.acetate_g_c[complex] + deltas[complex][4],
        };
        for (temporary[complex]) |value|
            if (!std.math.isFinite(value))
                return error.NonFinitePrimingRedistributionResult;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, complex|
            @field(state, field.name)[complex] = values[index];
    };
}

fn pairFlux(
    inputs: Inputs,
    values: []const f64,
    donor_soc: f64,
    receiver_soc: f64,
    donor: usize,
    receiver: usize,
    combined_soc: f64,
) !f64 {
    const left = values[donor] * receiver_soc;
    const right = values[receiver] * donor_soc;
    const flux = inputs.priming_transfer_rate_h * inputs.temperature_response *
        (left - right) / combined_soc * inputs.biochemical_time_fraction_h;
    inline for (.{ left, right, flux }) |value| if (!std.math.isFinite(value))
        return error.NonFinitePrimingRedistributionResult;
    return flux;
}

fn stateFieldIndex(comptime name: []const u8) usize {
    var index: usize = 0;
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        if (comptime std.mem.eql(u8, field.name, name)) return index;
        index += 1;
    };
    unreachable;
}

fn poolValues(inputs: Inputs, comptime pool: usize) []const f64 {
    return switch (pool) {
        0 => inputs.microbial_activity_respiration_g_c,
        1 => inputs.dissolved_organic_carbon_g_c,
        2 => inputs.dissolved_organic_nitrogen_g_n,
        3 => inputs.dissolved_organic_phosphorus_g_p,
        4 => inputs.acetate_g_c,
        else => unreachable,
    };
}

fn validate(state: *const State, inputs: Inputs) !void {
    inline for (.{
        inputs.soil_organic_carbon_g_c,          inputs.microbial_activity_respiration_g_c,
        inputs.dissolved_organic_carbon_g_c,     inputs.dissolved_organic_nitrogen_g_n,
        inputs.dissolved_organic_phosphorus_g_p, inputs.acetate_g_c,
    }) |values| {
        if (values.len != state.complex_count) return error.InvalidPrimingRedistributionDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidPrimingRedistributionInput;
    }
    inline for (.{
        inputs.temperature_response,        inputs.priming_transfer_rate_h,
        inputs.biochemical_time_fraction_h, inputs.negligible_soil_organic_carbon_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidPrimingRedistributionInput;
}

fn fixture() Inputs {
    return .{
        .soil_organic_carbon_g_c = &.{ 10, 20, 30 },
        .microbial_activity_respiration_g_c = &.{ 4, 2, 1 },
        .dissolved_organic_carbon_g_c = &.{ 8, 4, 2 },
        .dissolved_organic_nitrogen_g_n = &.{ 2, 1, 0.5 },
        .dissolved_organic_phosphorus_g_p = &.{ 1, 0.5, 0.25 },
        .acetate_g_c = &.{ 3, 1.5, 0.75 },
        .temperature_response = 0.5,
        .priming_transfer_rate_h = 0.1,
        .biochemical_time_fraction_h = 1,
        .negligible_soil_organic_carbon_g_c = 1e-12,
    };
}

test "priming redistribution conserves every dissolved pool" {
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    try calculate(&state, fixture());
    inline for (.{
        state.activity_transfer_g_c,           state.dissolved_carbon_transfer_g_c,
        state.dissolved_nitrogen_transfer_g_n, state.dissolved_phosphorus_transfer_g_p,
        state.acetate_transfer_g_c,
    }) |deltas| {
        var total: f64 = 0;
        for (deltas) |value| total += value;
        try std.testing.expectApproxEqAbs(0, total, 1e-12);
    }
}

test "priming moves high pool to SOC ratio toward lower ratio complexes" {
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expect(state.dissolved_carbon_transfer_g_c[0] < 0);
    try std.testing.expect(state.dissolved_carbon_transfer_g_c[2] > 0);
}

test "zero SOC complex is excluded from pair transfers" {
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    var inputs = fixture();
    inputs.soil_organic_carbon_g_c = &.{ 10, 0, 30 };
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.dissolved_carbon_transfer_g_c[1]);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    state.dissolved_carbon_transfer_g_c[0] = 7;
    var inputs = fixture();
    inputs.temperature_response = std.math.nan(f64);
    try std.testing.expectError(error.InvalidPrimingRedistributionInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.dissolved_carbon_transfer_g_c[0]);
}

test "NITRO 3003-3080 pre-cap overflow preserves redistribution state" {
    var state = try State.init(std.testing.allocator, 3);
    defer state.deinit();
    state.dissolved_carbon_transfer_g_c[0] = 7;
    var inputs = fixture();
    inputs.dissolved_organic_carbon_g_c =
        &.{ std.math.floatMax(f64), 4, 2 };
    inputs.soil_organic_carbon_g_c =
        &.{ 10, std.math.floatMax(f64), 30 };
    try std.testing.expectError(
        error.NonFinitePrimingRedistributionResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.dissolved_carbon_transfer_g_c[0]);
}
