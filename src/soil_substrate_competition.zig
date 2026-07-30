const std = @import("std");

pub const SourceChannel = enum(u8) {
    oxygen,
    ammonium_non_band,
    ammonium_band,
    nitrate_non_band,
    nitrate_band,
    dihydrogen_phosphate_non_band,
    dihydrogen_phosphate_band,
    hydrogen_phosphate_non_band,
    hydrogen_phosphate_band,
    dissolved_organic_carbon,
    dissolved_acetate_carbon,
};

pub const Inputs = struct {
    channel_count: usize,
    microbial_unit_count: usize,
    enabled: []const bool,
    previous_unit_demand_g: []const f64,
    previous_total_demand_g: []const f64,
    no_demand_fallback_fraction: []const f64,
    minimum_competition_fraction: f64,
    negligible_total_demand_g: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    channel_count: usize,
    microbial_unit_count: usize,
    competition_fraction: []f64,
    total_competition_fraction: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        channel_count: usize,
        microbial_unit_count: usize,
    ) !State {
        if (channel_count == 0 or microbial_unit_count == 0)
            return error.InvalidSubstrateCompetitionDimensions;
        const count = try std.math.mul(
            usize,
            channel_count,
            microbial_unit_count,
        );
        var state: State = undefined;
        state.allocator = allocator;
        state.channel_count = channel_count;
        state.microbial_unit_count = microbial_unit_count;
        state.competition_fraction = try allocator.alloc(f64, count);
        errdefer allocator.free(state.competition_fraction);
        state.total_competition_fraction =
            try allocator.alloc(f64, channel_count);
        @memset(state.competition_fraction, 0);
        @memset(state.total_competition_fraction, 0);
        return state;
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.competition_fraction);
        self.allocator.free(self.total_competition_fraction);
        self.* = undefined;
    }

    pub fn fraction(
        self: *const State,
        channel: usize,
        microbial_unit: usize,
    ) !f64 {
        if (channel >= self.channel_count or
            microbial_unit >= self.microbial_unit_count)
            return error.SubstrateCompetitionIndexOutOfRange;
        return self.competition_fraction[
            channel * self.microbial_unit_count + microbial_unit
        ];
    }
};

/// NITRO.F 666--786. Each historical substrate code is one runtime channel.
/// The caller supplies the source-specific fallback: FOMA, FOMA multiplied
/// by a solution-zone fraction, or FOMK. The source does not renormalize after
/// applying FMN, so `total_competition_fraction` intentionally may exceed one.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    var staged = try State.init(
        state.allocator,
        state.channel_count,
        state.microbial_unit_count,
    );
    defer staged.deinit();
    calculateValidated(&staged, inputs);
    try validateResult(&staged);
    @memcpy(state.competition_fraction, staged.competition_fraction);
    @memcpy(
        state.total_competition_fraction,
        staged.total_competition_fraction,
    );
}

fn calculateValidated(state: *State, inputs: Inputs) void {
    @memset(state.competition_fraction, 0);
    @memset(state.total_competition_fraction, 0);
    for (0..inputs.microbial_unit_count) |unit| {
        for (0..inputs.channel_count) |channel| {
            const total_demand = inputs.previous_total_demand_g[channel];
            const index = channel * inputs.microbial_unit_count + unit;
            if (!inputs.enabled[index]) continue;
            const raw_fraction =
                if (total_demand > inputs.negligible_total_demand_g)
                    inputs.previous_unit_demand_g[index] / total_demand
                else
                    inputs.no_demand_fallback_fraction[index];
            const fraction = @max(
                inputs.minimum_competition_fraction,
                raw_fraction,
            );
            state.competition_fraction[index] = fraction;
            state.total_competition_fraction[channel] += fraction;
        }
    }
}

fn validateResult(state: *const State) !void {
    for (state.competition_fraction) |fraction| {
        if (!std.math.isFinite(fraction) or fraction < 0)
            return error.NonFiniteSubstrateCompetitionResult;
    }
    for (state.total_competition_fraction) |fraction| {
        if (!std.math.isFinite(fraction) or fraction < 0)
            return error.NonFiniteSubstrateCompetitionResult;
    }
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.channel_count == 0 or inputs.microbial_unit_count == 0 or
        state.channel_count != inputs.channel_count or
        state.microbial_unit_count != inputs.microbial_unit_count)
        return error.InvalidSubstrateCompetitionDimensions;
    const count = try std.math.mul(
        usize,
        inputs.channel_count,
        inputs.microbial_unit_count,
    );
    if (inputs.enabled.len != count or
        inputs.previous_unit_demand_g.len != count or
        inputs.no_demand_fallback_fraction.len != count or
        inputs.previous_total_demand_g.len != inputs.channel_count)
        return error.InvalidSubstrateCompetitionDimensions;
    for (inputs.previous_total_demand_g) |value| {
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSubstrateCompetitionInput;
    }
    for (inputs.previous_unit_demand_g, 0..) |value, index| {
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSubstrateCompetitionInput;
        const channel = index / inputs.microbial_unit_count;
        const total = inputs.previous_total_demand_g[channel];
        const tolerance =
            64 * std.math.floatEps(f64) * @max(1, total);
        if (inputs.enabled[index] and
            total > inputs.negligible_total_demand_g and
            value > total + tolerance)
            return error.UnitDemandExceedsTotalDemand;
    }
    for (inputs.no_demand_fallback_fraction) |value| {
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return error.InvalidSubstrateCompetitionInput;
    }
    if (!std.math.isFinite(inputs.minimum_competition_fraction) or
        inputs.minimum_competition_fraction < 0 or
        inputs.minimum_competition_fraction > 1 or
        !std.math.isFinite(inputs.negligible_total_demand_g) or
        inputs.negligible_total_demand_g < 0)
        return error.InvalidSubstrateCompetitionInput;
}

test "previous demand and biomass fallback branches reproduce NITRO competition" {
    var state = try State.init(std.testing.allocator, 3, 3);
    defer state.deinit();
    try calculate(&state, .{
        .channel_count = 3,
        .microbial_unit_count = 3,
        .enabled = &.{
            true, true,  true,
            true, true,  true,
            true, false, true,
        },
        .previous_unit_demand_g = &.{
            2,    1,  0,
            0,    0,  0,
            0.01, 99, 0,
        },
        .previous_total_demand_g = &.{ 4, 0, 1 },
        .no_demand_fallback_fraction = &.{
            0.2, 0.3, 0.5,
            0.2, 0.3, 0.5,
            0.2, 0.3, 0.5,
        },
        .minimum_competition_fraction = 0.05,
        .negligible_total_demand_g = 1e-12,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), try state.fraction(0, 0), 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), try state.fraction(0, 1), 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), try state.fraction(0, 2), 1e-15);
    try std.testing.expectEqualSlices(f64, &.{ 0.2, 0.3, 0.5 }, state.competition_fraction[3..6]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), try state.fraction(2, 0), 1e-15);
    try std.testing.expectEqual(@as(f64, 0), try state.fraction(2, 1));
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), try state.fraction(2, 2), 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), state.total_competition_fraction[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), state.total_competition_fraction[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), state.total_competition_fraction[2], 1e-15);
}

test "minimum shares remain unnormalized exactly as source" {
    var state = try State.init(std.testing.allocator, 1, 4);
    defer state.deinit();
    try calculate(&state, .{
        .channel_count = 1,
        .microbial_unit_count = 4,
        .enabled = &.{ true, true, true, true },
        .previous_unit_demand_g = &.{ 0, 0, 0, 0 },
        .previous_total_demand_g = &.{1},
        .no_demand_fallback_fraction = &.{ 0, 0, 0, 0 },
        .minimum_competition_fraction = 0.3,
        .negligible_total_demand_g = 1e-12,
    });
    try std.testing.expectEqualSlices(f64, &.{ 0.3, 0.3, 0.3, 0.3 }, state.competition_fraction);
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), state.total_competition_fraction[0], 1e-15);
}

test "invalid late channel demand preserves prior competition state" {
    var state = try State.init(std.testing.allocator, 2, 1);
    defer state.deinit();
    state.competition_fraction[0] = 7;
    try std.testing.expectError(error.UnitDemandExceedsTotalDemand, calculate(&state, .{
        .channel_count = 2,
        .microbial_unit_count = 1,
        .enabled = &.{ true, true },
        .previous_unit_demand_g = &.{ 0.5, 2 },
        .previous_total_demand_g = &.{ 1, 1 },
        .no_demand_fallback_fraction = &.{ 1, 1 },
        .minimum_competition_fraction = 0.01,
        .negligible_total_demand_g = 1e-12,
    }));
    try std.testing.expectEqual(@as(f64, 7), state.competition_fraction[0]);
}

test "NITRO 678-745 competition publishes all channels atomically" {
    var state = try State.init(std.testing.allocator, 2, 2);
    defer state.deinit();
    @memset(state.competition_fraction, 7);
    @memset(state.total_competition_fraction, 8);
    try std.testing.expectError(
        error.UnitDemandExceedsTotalDemand,
        calculate(&state, .{
            .channel_count = 2,
            .microbial_unit_count = 2,
            .enabled = &.{ true, true, true, true },
            .previous_unit_demand_g = &.{ 0.5, 0.5, 0.5, 2 },
            .previous_total_demand_g = &.{ 1, 1 },
            .no_demand_fallback_fraction = &.{ 0.5, 0.5, 0.5, 0.5 },
            .minimum_competition_fraction = 0.01,
            .negligible_total_demand_g = 1e-12,
        }),
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 7, 7, 7, 7 },
        state.competition_fraction,
    );
    try std.testing.expectEqualSlices(
        f64,
        &.{ 8, 8 },
        state.total_competition_fraction,
    );
}
