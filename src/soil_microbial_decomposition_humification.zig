const std = @import("std");

pub const Inputs = struct {
    litterfall_carbon_g_c: []const f64,
    litterfall_nitrogen_g_n: []const f64,
    litterfall_phosphorus_g_p: []const f64,
    humification_fraction: []const f64,
    humus_nitrogen_to_carbon_ratio_g_n_per_g_c: f64,
    humus_phosphorus_to_carbon_ratio_g_p_per_g_c: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    component_count: usize,
    humus_carbon_g_c: []f64,
    humus_nitrogen_g_n: []f64,
    humus_phosphorus_g_p: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        unit_count: usize,
        component_count: usize,
    ) !State {
        if (unit_count == 0 or component_count == 0)
            return error.InvalidDecompositionHumificationDimensions;
        const item_count = std.math.mul(usize, unit_count, component_count) catch
            return error.InvalidDecompositionHumificationDimensions;
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

/// Exact NITRO.F 2706--2722 humification of microbial decomposition products.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const item_count = state.unit_count * state.component_count;
    const temporary = try state.allocator.alloc([3]f64, item_count);
    defer state.allocator.free(temporary);

    for (0..state.unit_count) |unit| {
        for (0..state.component_count) |component| {
            const item = unit * state.component_count + component;
            const carbon_by_fraction = inputs.litterfall_carbon_g_c[item] *
                inputs.humification_fraction[unit];
            const humus_carbon = @max(0, carbon_by_fraction);
            const nitrogen_by_fraction = inputs.litterfall_nitrogen_g_n[item] *
                inputs.humification_fraction[unit];
            const nitrogen_by_ratio = humus_carbon *
                inputs.humus_nitrogen_to_carbon_ratio_g_n_per_g_c;
            const phosphorus_by_fraction = inputs.litterfall_phosphorus_g_p[item] *
                inputs.humification_fraction[unit];
            const phosphorus_by_ratio = humus_carbon *
                inputs.humus_phosphorus_to_carbon_ratio_g_p_per_g_c;
            inline for (.{
                carbon_by_fraction,
                nitrogen_by_fraction,
                nitrogen_by_ratio,
                phosphorus_by_fraction,
                phosphorus_by_ratio,
            }) |value| if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteDecompositionHumificationResult;
            const humus_nitrogen = @max(
                0,
                @min(nitrogen_by_fraction, nitrogen_by_ratio),
            );
            const humus_phosphorus = @max(
                0,
                @min(phosphorus_by_fraction, phosphorus_by_ratio),
            );
            temporary[item] = .{
                humus_carbon,
                humus_nitrogen,
                humus_phosphorus,
            };
            for (temporary[item]) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteDecompositionHumificationResult;
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
        return error.InvalidDecompositionHumificationDimensions;
    inline for (.{
        inputs.litterfall_carbon_g_c,
        inputs.litterfall_nitrogen_g_n,
        inputs.litterfall_phosphorus_g_p,
    }) |values| {
        if (values.len != item_count)
            return error.InvalidDecompositionHumificationDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDecompositionHumificationInput;
    }
    if (inputs.humification_fraction.len != state.unit_count)
        return error.InvalidDecompositionHumificationDimensions;
    for (inputs.humification_fraction) |value|
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return error.InvalidDecompositionHumificationInput;
    inline for (.{
        inputs.humus_nitrogen_to_carbon_ratio_g_n_per_g_c,
        inputs.humus_phosphorus_to_carbon_ratio_g_p_per_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidDecompositionHumificationInput;
}

fn fixture() Inputs {
    return .{
        .litterfall_carbon_g_c = &.{ 4, 2 },
        .litterfall_nitrogen_g_n = &.{ 1, 0.5 },
        .litterfall_phosphorus_g_p = &.{ 0.2, 0.1 },
        .humification_fraction = &.{0.5},
        .humus_nitrogen_to_carbon_ratio_g_n_per_g_c = 0.1,
        .humus_phosphorus_to_carbon_ratio_g_p_per_g_c = 0.02,
    };
}

test "NITRO 2706-2722 preserves component humification caps" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(2, state.humus_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.2, state.humus_nitrogen_g_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.04, state.humus_phosphorus_g_p[0], 1e-12);
    try std.testing.expectApproxEqAbs(1, state.humus_carbon_g_c[1], 1e-12);
}

test "NITRO 2706-2722 independently caps humus nutrients" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.litterfall_nitrogen_g_n = &.{ 0.1, 0.5 };
    try calculate(&state, inputs);
    try std.testing.expectApproxEqAbs(0.05, state.humus_nitrogen_g_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.04, state.humus_phosphorus_g_p[0], 1e-12);
}

test "NITRO 2706-2722 derived overflow preserves humification state" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.humus_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.litterfall_carbon_g_c = &.{ std.math.floatMax(f64), 2 };
    inputs.humus_nitrogen_to_carbon_ratio_g_n_per_g_c = std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteDecompositionHumificationResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.humus_carbon_g_c[0]);
}
