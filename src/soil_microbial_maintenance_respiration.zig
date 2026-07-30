const std = @import("std");

pub const Inputs = struct {
    labile_microbial_nitrogen_g_n: []const f64,
    resistant_microbial_nitrogen_g_n: []const f64,
    temperature_response: []const f64,
    ph_response: []const f64,
    low_carbon_response: []const f64,
    aerobic_respiration_g_c: []const f64,
    specific_maintenance_respiration_g_c_per_g_n_h: f64,
    biochemical_time_fraction_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    effective_specific_maintenance_g_c_per_g_n: []f64,
    labile_maintenance_respiration_g_c: []f64,
    resistant_maintenance_respiration_g_c: []f64,
    total_maintenance_respiration_g_c: []f64,
    growth_respiration_g_c: []f64,
    senescence_respiration_g_c: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidMaintenanceRespirationDimensions;
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

/// Exact NITRO.F 2483--2509 microbial maintenance/growth/senescence split.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([6]f64, state.unit_count);
    defer state.allocator.free(temporary);
    for (0..state.unit_count) |unit| {
        const effective_specific = inputs.specific_maintenance_respiration_g_c_per_g_n_h *
            inputs.temperature_response[unit] * inputs.ph_response[unit] *
            inputs.biochemical_time_fraction_h * inputs.low_carbon_response[unit];
        const labile = inputs.labile_microbial_nitrogen_g_n[unit] * effective_specific;
        const resistant = inputs.resistant_microbial_nitrogen_g_n[unit] * effective_specific;
        const maintenance = labile + resistant;
        const aerobic = inputs.aerobic_respiration_g_c[unit];
        temporary[unit] = .{
            effective_specific,
            labile,
            resistant,
            maintenance,
            @max(0, aerobic - maintenance),
            @max(0, maintenance - aerobic),
        };
        for (temporary[unit]) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteMaintenanceRespirationResult;
    }
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
        if (values.len != n) return error.InvalidMaintenanceRespirationDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidMaintenanceRespirationInput;
    };
    inline for (.{
        inputs.specific_maintenance_respiration_g_c_per_g_n_h,
        inputs.biochemical_time_fraction_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidMaintenanceRespirationInput;
}

fn fixture() Inputs {
    return .{
        .labile_microbial_nitrogen_g_n = &.{2},
        .resistant_microbial_nitrogen_g_n = &.{3},
        .temperature_response = &.{0.5},
        .ph_response = &.{0.8},
        .low_carbon_response = &.{0.5},
        .aerobic_respiration_g_c = &.{4},
        .specific_maintenance_respiration_g_c_per_g_n_h = 1,
        .biochemical_time_fraction_h = 1,
    };
}

test "maintenance below aerobic respiration leaves growth respiration" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(0.2, state.effective_specific_maintenance_g_c_per_g_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(1, state.total_maintenance_respiration_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(3, state.growth_respiration_g_c[0], 1e-12);
    try std.testing.expectEqual(0, state.senescence_respiration_g_c[0]);
}

test "maintenance above aerobic respiration becomes senescence respiration" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.aerobic_respiration_g_c = &.{0.25};
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.growth_respiration_g_c[0]);
    try std.testing.expectApproxEqAbs(0.75, state.senescence_respiration_g_c[0], 1e-12);
}

test "runtime units are independent" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.labile_microbial_nitrogen_g_n = &.{ 2, 4 };
    inputs.resistant_microbial_nitrogen_g_n = &.{ 3, 6 };
    inputs.temperature_response = &.{ 0.5, 0.5 };
    inputs.ph_response = &.{ 0.8, 0.8 };
    inputs.low_carbon_response = &.{ 0.5, 0.5 };
    inputs.aerobic_respiration_g_c = &.{ 4, 4 };
    try calculate(&state, inputs);
    try std.testing.expectApproxEqAbs(2, state.total_maintenance_respiration_g_c[1], 1e-12);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.total_maintenance_respiration_g_c[0] = 7;
    var inputs = fixture();
    inputs.ph_response = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidMaintenanceRespirationInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.total_maintenance_respiration_g_c[0]);
}

test "NITRO 2483-2509 derived overflow preserves maintenance state" {
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.total_maintenance_respiration_g_c[0] = 7;
    var inputs = fixture();
    inputs.specific_maintenance_respiration_g_c_per_g_n_h = std.math.floatMax(f64);
    inputs.temperature_response = &.{2};
    try std.testing.expectError(
        error.NonFiniteMaintenanceRespirationResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.total_maintenance_respiration_g_c[0]);
}
