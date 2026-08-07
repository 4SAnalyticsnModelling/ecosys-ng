const std = @import("std");

pub const Inputs = struct {
    complex_count: usize,
    population_count: usize,
    biomass_fraction_count: usize,
    soil_organic_carbon_g_c: []const f64,
    biological_activity_fraction: []const f64,
    microbial_carbon_g_c: []const f64,
    microbial_nitrogen_g_n: []const f64,
    microbial_phosphorus_g_p: []const f64,
    priming_transfer_rate_h: f64,
    biochemical_time_fraction_h: f64,
    negligible_soil_organic_carbon_g_c: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    item_count: usize,
    carbon_transfer_g_c: []f64,
    nitrogen_transfer_g_n: []f64,
    phosphorus_transfer_g_p: []f64,
    redistributed_microbial_carbon_g_c: []f64,
    redistributed_microbial_nitrogen_g_n: []f64,
    redistributed_microbial_phosphorus_g_p: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        complex_count: usize,
        population_count: usize,
        biomass_fraction_count: usize,
    ) !State {
        if (complex_count == 0 or population_count == 0 or biomass_fraction_count == 0)
            return error.InvalidMicrobialPrimingDimensions;
        const cp = std.math.mul(usize, complex_count, population_count) catch
            return error.InvalidMicrobialPrimingDimensions;
        const item_count = std.math.mul(usize, cp, biomass_fraction_count) catch
            return error.InvalidMicrobialPrimingDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.item_count = item_count;
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
            @field(state, field.name) = try allocator.alloc(f64, item_count);
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

/// Exact NITRO.F 3082--3138 microbial C/N/P priming redistribution.
/// Flattening order is complex, population, biomass fraction.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const deltas = try state.allocator.alloc([3]f64, state.item_count);
    defer state.allocator.free(deltas);
    @memset(deltas, @splat(0));
    for (0..inputs.complex_count) |donor_complex| {
        for (donor_complex + 1..inputs.complex_count) |receiver_complex| {
            const donor_soc = inputs.soil_organic_carbon_g_c[donor_complex];
            const receiver_soc = inputs.soil_organic_carbon_g_c[receiver_complex];
            if (donor_soc <= inputs.negligible_soil_organic_carbon_g_c or
                receiver_soc <= inputs.negligible_soil_organic_carbon_g_c)
                continue;
            const combined_soc = donor_soc + receiver_soc;
            for (0..inputs.population_count) |population| {
                const activity_index = donor_complex * inputs.population_count + population;
                for (0..inputs.biomass_fraction_count) |fraction| {
                    const donor = itemIndex(inputs, donor_complex, population, fraction);
                    const receiver = itemIndex(inputs, receiver_complex, population, fraction);
                    const fluxes: [3]f64 = .{
                        try pairFlux(inputs, activity_index, inputs.microbial_carbon_g_c, donor_soc, receiver_soc, donor, receiver, combined_soc),
                        try pairFlux(inputs, activity_index, inputs.microbial_nitrogen_g_n, donor_soc, receiver_soc, donor, receiver, combined_soc),
                        try pairFlux(inputs, activity_index, inputs.microbial_phosphorus_g_p, donor_soc, receiver_soc, donor, receiver, combined_soc),
                    };
                    inline for (0..3) |element| {
                        const values = poolValues(inputs, element);
                        const flux = fluxes[element];
                        if (values[donor] + deltas[donor][element] - flux > 0 and
                            values[receiver] + deltas[receiver][element] + flux > 0)
                        {
                            const donor_delta = deltas[donor][element] - flux;
                            const receiver_delta = deltas[receiver][element] + flux;
                            if (!std.math.isFinite(donor_delta) or
                                !std.math.isFinite(receiver_delta))
                                return error.NonFiniteMicrobialPrimingResult;
                            deltas[donor][element] = donor_delta;
                            deltas[receiver][element] = receiver_delta;
                        }
                    }
                }
            }
        }
    }
    const temporary = try state.allocator.alloc([6]f64, state.item_count);
    defer state.allocator.free(temporary);
    for (0..state.item_count) |item| {
        temporary[item] = .{
            deltas[item][0],
            deltas[item][1],
            deltas[item][2],
            inputs.microbial_carbon_g_c[item] + deltas[item][0],
            inputs.microbial_nitrogen_g_n[item] + deltas[item][1],
            inputs.microbial_phosphorus_g_p[item] + deltas[item][2],
        };
        for (temporary[item]) |value|
            if (!std.math.isFinite(value))
                return error.NonFiniteMicrobialPrimingResult;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const field_index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, item|
            @field(state, field.name)[item] = values[field_index];
    };
}

fn itemIndex(inputs: Inputs, complex: usize, population: usize, fraction: usize) usize {
    return (complex * inputs.population_count + population) *
        inputs.biomass_fraction_count + fraction;
}

fn pairFlux(
    inputs: Inputs,
    activity_index: usize,
    values: []const f64,
    donor_soc: f64,
    receiver_soc: f64,
    donor: usize,
    receiver: usize,
    combined_soc: f64,
) !f64 {
    const left = values[donor] * receiver_soc;
    const right = values[receiver] * donor_soc;
    const flux = inputs.priming_transfer_rate_h *
        inputs.biological_activity_fraction[activity_index] *
        (left - right) / combined_soc * inputs.biochemical_time_fraction_h;
    inline for (.{ left, right, flux }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteMicrobialPrimingResult;
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

fn poolValues(inputs: Inputs, comptime element: usize) []const f64 {
    return switch (element) {
        0 => inputs.microbial_carbon_g_c,
        1 => inputs.microbial_nitrogen_g_n,
        2 => inputs.microbial_phosphorus_g_p,
        else => unreachable,
    };
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.complex_count == 0 or inputs.population_count == 0 or
        inputs.biomass_fraction_count == 0)
        return error.InvalidMicrobialPrimingDimensions;
    const cp = std.math.mul(usize, inputs.complex_count, inputs.population_count) catch
        return error.InvalidMicrobialPrimingDimensions;
    const expected = std.math.mul(usize, cp, inputs.biomass_fraction_count) catch
        return error.InvalidMicrobialPrimingDimensions;
    if (expected != state.item_count or inputs.soil_organic_carbon_g_c.len != inputs.complex_count or
        inputs.biological_activity_fraction.len != cp)
        return error.InvalidMicrobialPrimingDimensions;
    inline for (.{
        inputs.soil_organic_carbon_g_c,  inputs.biological_activity_fraction,
        inputs.microbial_carbon_g_c,     inputs.microbial_nitrogen_g_n,
        inputs.microbial_phosphorus_g_p,
    }) |values| for (values) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidMicrobialPrimingInput;
    inline for (.{
        inputs.microbial_carbon_g_c,     inputs.microbial_nitrogen_g_n,
        inputs.microbial_phosphorus_g_p,
    }) |values| if (values.len != expected) return error.InvalidMicrobialPrimingDimensions;
    inline for (.{
        inputs.priming_transfer_rate_h,            inputs.biochemical_time_fraction_h,
        inputs.negligible_soil_organic_carbon_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidMicrobialPrimingInput;
}

fn fixture() Inputs {
    return .{
        .complex_count = 3,
        .population_count = 1,
        .biomass_fraction_count = 2,
        .soil_organic_carbon_g_c = &.{ 10, 20, 30 },
        .biological_activity_fraction = &.{ 0.5, 0.4, 0.3 },
        .microbial_carbon_g_c = &.{ 4, 2, 2, 1, 1, 0.5 },
        .microbial_nitrogen_g_n = &.{ 2, 1, 1, 0.5, 0.5, 0.25 },
        .microbial_phosphorus_g_p = &.{ 1, 0.5, 0.5, 0.25, 0.25, 0.125 },
        .priming_transfer_rate_h = 0.1,
        .biochemical_time_fraction_h = 1,
        .negligible_soil_organic_carbon_g_c = 1e-12,
    };
}

test "microbial priming conserves C N and P across runtime complexes" {
    var state = try State.init(std.testing.allocator, 3, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    inline for (.{
        state.carbon_transfer_g_c,     state.nitrogen_transfer_g_n,
        state.phosphorus_transfer_g_p,
    }) |values| {
        var total: f64 = 0;
        for (values) |value| total += value;
        try std.testing.expectApproxEqAbs(0, total, 1e-12);
    }
}

test "high biomass to SOC ratio transfers toward lower ratio complexes" {
    var state = try State.init(std.testing.allocator, 3, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expect(state.carbon_transfer_g_c[0] < 0);
    try std.testing.expect(state.carbon_transfer_g_c[4] > 0);
}

test "zero SOC complex is excluded from microbial transfers" {
    var state = try State.init(std.testing.allocator, 3, 1, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.soil_organic_carbon_g_c = &.{ 10, 0, 30 };
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.carbon_transfer_g_c[2]);
    try std.testing.expectEqual(0, state.carbon_transfer_g_c[3]);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 3, 1, 2);
    defer state.deinit();
    state.carbon_transfer_g_c[0] = 7;
    var inputs = fixture();
    inputs.priming_transfer_rate_h = std.math.nan(f64);
    try std.testing.expectError(error.InvalidMicrobialPrimingInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.carbon_transfer_g_c[0]);
}

test "NITRO 3082-3138 derived overflow preserves microbial priming state" {
    var state = try State.init(std.testing.allocator, 3, 1, 2);
    defer state.deinit();
    state.carbon_transfer_g_c[0] = 7;
    var inputs = fixture();
    inputs.microbial_carbon_g_c =
        &.{ std.math.floatMax(f64), 2, 2, 1, 1, 0.5 };
    inputs.soil_organic_carbon_g_c =
        &.{ 10, std.math.floatMax(f64), 30 };
    try std.testing.expectError(
        error.NonFiniteMicrobialPrimingResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.carbon_transfer_g_c[0]);
}
