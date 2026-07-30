const std = @import("std");

pub const Inputs = struct {
    litterfall_carbon_g_c: []const f64,
    litterfall_nitrogen_g_n: []const f64,
    litterfall_phosphorus_g_p: []const f64,
    humus_carbon_g_c: []const f64,
    humus_nitrogen_g_n: []const f64,
    humus_phosphorus_g_p: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    component_count: usize,
    residue_carbon_g_c: []f64,
    residue_nitrogen_g_n: []f64,
    residue_phosphorus_g_p: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        unit_count: usize,
        component_count: usize,
    ) !State {
        if (unit_count == 0 or component_count == 0)
            return error.InvalidDecompositionResidueDimensions;
        const item_count = std.math.mul(usize, unit_count, component_count) catch
            return error.InvalidDecompositionResidueDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.unit_count = unit_count;
        state.component_count = component_count;
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

/// Exact NITRO.F 2737--2749 non-humified microbial-residue allocation.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const item_count = state.unit_count * state.component_count;
    const temporary = try state.allocator.alloc([3]f64, item_count);
    defer state.allocator.free(temporary);

    for (0..state.unit_count) |unit| {
        for (0..state.component_count) |component| {
            const item = unit * state.component_count + component;
            temporary[item] = .{
                inputs.litterfall_carbon_g_c[item] - inputs.humus_carbon_g_c[item],
                inputs.litterfall_nitrogen_g_n[item] - inputs.humus_nitrogen_g_n[item],
                inputs.litterfall_phosphorus_g_p[item] -
                    inputs.humus_phosphorus_g_p[item],
            };
            for (temporary[item]) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.InconsistentDecompositionResidueAllocation;
        }
    }

    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, item| @field(state, field.name)[item] = values[index];
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
    const item_count = std.math.mul(usize, state.unit_count, state.component_count) catch
        return error.InvalidDecompositionResidueDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const values = @field(inputs, field.name);
        if (values.len != item_count)
            return error.InvalidDecompositionResidueDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDecompositionResidueInput;
    }
}

fn fixture() Inputs {
    return .{
        .litterfall_carbon_g_c = &.{ 4, 2 },
        .litterfall_nitrogen_g_n = &.{ 1, 0.5 },
        .litterfall_phosphorus_g_p = &.{ 0.2, 0.1 },
        .humus_carbon_g_c = &.{ 2, 1 },
        .humus_nitrogen_g_n = &.{ 0.2, 0.1 },
        .humus_phosphorus_g_p = &.{ 0.04, 0.02 },
    };
}

test "NITRO 2737-2749 preserves component residue subtraction order" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(2, state.residue_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.8, state.residue_nitrogen_g_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.16, state.residue_phosphorus_g_p[0], 1e-12);
    try std.testing.expectApproxEqAbs(1, state.residue_carbon_g_c[1], 1e-12);
}

test "NITRO 2737-2749 humus plus residue closes each element" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    const inputs = fixture();
    try calculate(&state, inputs);
    for (0..2) |item| {
        try std.testing.expectApproxEqAbs(
            inputs.litterfall_carbon_g_c[item],
            inputs.humus_carbon_g_c[item] + state.residue_carbon_g_c[item],
            1e-12,
        );
        try std.testing.expectApproxEqAbs(
            inputs.litterfall_nitrogen_g_n[item],
            inputs.humus_nitrogen_g_n[item] + state.residue_nitrogen_g_n[item],
            1e-12,
        );
        try std.testing.expectApproxEqAbs(
            inputs.litterfall_phosphorus_g_p[item],
            inputs.humus_phosphorus_g_p[item] + state.residue_phosphorus_g_p[item],
            1e-12,
        );
    }
}

test "NITRO 2737-2749 inconsistent humification preserves residue state" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.residue_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.humus_nitrogen_g_n = &.{ 2, 0.1 };
    try std.testing.expectError(
        error.InconsistentDecompositionResidueAllocation,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.residue_carbon_g_c[0]);
}
