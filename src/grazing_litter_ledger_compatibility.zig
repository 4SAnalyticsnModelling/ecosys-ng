const std = @import("std");

pub const State = struct {
    allocator: std.mem.Allocator,
    plant_count: usize,
    current_carbon_g_c: []f64,
    current_nitrogen_g_n: []f64,
    current_phosphorus_g_p: []f64,
    cumulative_carbon_g_c: []f64,
    cumulative_nitrogen_g_n: []f64,
    cumulative_phosphorus_g_p: []f64,
    cumulative_aboveground_carbon_g_c: []f64,
    cumulative_aboveground_nitrogen_g_n: []f64,
    cumulative_aboveground_phosphorus_g_p: []f64,
    current_surface_carbon_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, plant_count: usize) !State {
        if (plant_count == 0) return error.InvalidGrazingLitterLedgerDimensions;
        const values = try allocator.alloc(f64, try std.math.mul(usize, plant_count, 10));
        @memset(values, 0);
        var result: State = .{
            .allocator = allocator,
            .plant_count = plant_count,
            .current_carbon_g_c = undefined,
            .current_nitrogen_g_n = undefined,
            .current_phosphorus_g_p = undefined,
            .cumulative_carbon_g_c = undefined,
            .cumulative_nitrogen_g_n = undefined,
            .cumulative_phosphorus_g_p = undefined,
            .cumulative_aboveground_carbon_g_c = undefined,
            .cumulative_aboveground_nitrogen_g_n = undefined,
            .cumulative_aboveground_phosphorus_g_p = undefined,
            .current_surface_carbon_g_c = undefined,
        };
        inline for (@typeInfo(State).@"struct".fields, 0..) |field, field_index| if (field.type == []f64) {
            const index = field_index - 2;
            @field(result, field.name) = values[index * plant_count .. (index + 1) * plant_count];
        };
        return result;
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.current_carbon_g_c.ptr[0 .. self.plant_count * 10]);
        self.* = undefined;
    }
};

pub const Inputs = struct {
    returned_carbon_g_c: []const f64,
    returned_nitrogen_g_n: []const f64,
    returned_phosphorus_g_p: []const f64,
    harvested_litter_carbon_g_c: []const f64,
    harvested_litter_nitrogen_g_n: []const f64,
    harvested_litter_phosphorus_g_p: []const f64,
};

/// Exact GROSUB 10840–10854 compatibility publication. The source assigns
/// TZSN0/TPSN0 from the just-incremented total ledgers, unlike carbon.
pub fn accumulate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    for (0..state.plant_count) |plant| _ = try calculate(state, inputs, plant);
    for (0..state.plant_count) |plant| {
        const next = calculate(state, inputs, plant) catch unreachable;
        inline for (@typeInfo(Next).@"struct".fields) |field|
            @field(state, field.name)[plant] = @field(next, field.name);
    }
}

const Next = struct {
    current_carbon_g_c: f64,
    current_nitrogen_g_n: f64,
    current_phosphorus_g_p: f64,
    cumulative_carbon_g_c: f64,
    cumulative_nitrogen_g_n: f64,
    cumulative_phosphorus_g_p: f64,
    cumulative_aboveground_carbon_g_c: f64,
    cumulative_aboveground_nitrogen_g_n: f64,
    cumulative_aboveground_phosphorus_g_p: f64,
    current_surface_carbon_g_c: f64,
};

fn calculate(state: *const State, inputs: Inputs, plant: usize) !Next {
    const carbon = inputs.returned_carbon_g_c[plant] + inputs.harvested_litter_carbon_g_c[plant];
    const nitrogen = inputs.returned_nitrogen_g_n[plant] + inputs.harvested_litter_nitrogen_g_n[plant];
    const phosphorus = inputs.returned_phosphorus_g_p[plant] + inputs.harvested_litter_phosphorus_g_p[plant];
    const cumulative_carbon = state.cumulative_carbon_g_c[plant] + carbon;
    const cumulative_nitrogen = state.cumulative_nitrogen_g_n[plant] + nitrogen;
    const cumulative_phosphorus = state.cumulative_phosphorus_g_p[plant] + phosphorus;
    const next: Next = .{
        .current_carbon_g_c = state.current_carbon_g_c[plant] + carbon,
        .current_nitrogen_g_n = state.current_nitrogen_g_n[plant] + nitrogen,
        .current_phosphorus_g_p = state.current_phosphorus_g_p[plant] + phosphorus,
        .cumulative_carbon_g_c = cumulative_carbon,
        .cumulative_nitrogen_g_n = cumulative_nitrogen,
        .cumulative_phosphorus_g_p = cumulative_phosphorus,
        .cumulative_aboveground_carbon_g_c = state.cumulative_aboveground_carbon_g_c[plant] + carbon,
        .cumulative_aboveground_nitrogen_g_n = cumulative_nitrogen + nitrogen,
        .cumulative_aboveground_phosphorus_g_p = cumulative_phosphorus + phosphorus,
        .current_surface_carbon_g_c = state.current_surface_carbon_g_c[plant] + carbon,
    };
    inline for (@typeInfo(Next).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.NonFiniteGrazingLitterLedger;
    return next;
}

fn validate(state: *const State, inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const values = @field(inputs, field.name);
        if (values.len != state.plant_count) return error.InvalidGrazingLitterLedgerDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidGrazingLitterLedgerInput;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        if (@field(state, field.name).len != state.plant_count)
            return error.InvalidGrazingLitterLedgerDimensions;
        for (@field(state, field.name)) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidGrazingLitterLedgerState;
    };
}

test "compatibility ledger preserves source carbon versus nitrogen phosphorus asymmetry" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.cumulative_carbon_g_c[0] = 10;
    state.cumulative_nitrogen_g_n[0] = 20;
    state.cumulative_phosphorus_g_p[0] = 30;
    state.cumulative_aboveground_carbon_g_c[0] = 40;
    state.cumulative_aboveground_nitrogen_g_n[0] = 50;
    state.cumulative_aboveground_phosphorus_g_p[0] = 60;
    try accumulate(&state, .{
        .returned_carbon_g_c = &.{1},
        .returned_nitrogen_g_n = &.{2},
        .returned_phosphorus_g_p = &.{3},
        .harvested_litter_carbon_g_c = &.{4},
        .harvested_litter_nitrogen_g_n = &.{5},
        .harvested_litter_phosphorus_g_p = &.{6},
    });
    try std.testing.expectEqual(@as(f64, 45), state.cumulative_aboveground_carbon_g_c[0]);
    try std.testing.expectEqual(@as(f64, 34), state.cumulative_aboveground_nitrogen_g_n[0]);
    try std.testing.expectEqual(@as(f64, 48), state.cumulative_aboveground_phosphorus_g_p[0]);
}

test "current and cumulative totals preserve two-term source addition order" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try accumulate(&state, .{
        .returned_carbon_g_c = &.{ 1, 2 },
        .returned_nitrogen_g_n = &.{ 0.1, 0.2 },
        .returned_phosphorus_g_p = &.{ 0.01, 0.02 },
        .harvested_litter_carbon_g_c = &.{ 3, 4 },
        .harvested_litter_nitrogen_g_n = &.{ 0.3, 0.4 },
        .harvested_litter_phosphorus_g_p = &.{ 0.03, 0.04 },
    });
    try std.testing.expectEqualSlices(f64, &.{ 4, 6 }, state.current_carbon_g_c);
    try std.testing.expectEqualSlices(f64, &.{ 4, 6 }, state.cumulative_carbon_g_c);
    try std.testing.expectEqualSlices(f64, &.{ 4, 6 }, state.current_surface_carbon_g_c);
}

test "late overflow preserves all current and cumulative ledgers" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) @memset(@field(state, field.name), 7);
    const huge = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteGrazingLitterLedger, accumulate(&state, .{
        .returned_carbon_g_c = &.{ 1, huge },
        .returned_nitrogen_g_n = &.{ 1, 1 },
        .returned_phosphorus_g_p = &.{ 1, 1 },
        .harvested_litter_carbon_g_c = &.{ 1, huge },
        .harvested_litter_nitrogen_g_n = &.{ 1, 1 },
        .harvested_litter_phosphorus_g_p = &.{ 1, 1 },
    }));
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) for (@field(state, field.name)) |value| try std.testing.expectEqual(@as(f64, 7), value);
}
