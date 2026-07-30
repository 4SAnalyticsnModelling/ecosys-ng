const std = @import("std");

pub const Inputs = struct {
    complex_count: usize,
    population_count: usize,
    structural_fraction_count: usize,
    surface_litter_layer: bool,
    humus_partition_fraction: []const f64,
    basal_humification_carbon_g_c: []const f64,
    basal_humification_nitrogen_g_n: []const f64,
    basal_humification_phosphorus_g_p: []const f64,
    senescence_humification_carbon_g_c: []const f64,
    senescence_humification_nitrogen_g_n: []const f64,
    senescence_humification_phosphorus_g_p: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    humus_fraction_count: usize,
    carbon_addition_g_c: []f64,
    colonized_carbon_addition_g_c: []f64,
    nitrogen_addition_g_n: []f64,
    phosphorus_addition_g_p: []f64,
    target_is_surface_soil_layer: bool,

    pub fn init(allocator: std.mem.Allocator, humus_fraction_count: usize) !State {
        if (humus_fraction_count != 2) return error.InvalidHumusTransferDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.humus_fraction_count = humus_fraction_count;
        state.target_is_surface_soil_layer = false;
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
            @field(state, field.name) = try allocator.alloc(f64, humus_fraction_count);
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

/// Exact NITRO.F 3741--3792 microbial humification transfer to humus K=4.
/// Traversal follows source K, then N, then structural fraction M.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([4]f64, state.humus_fraction_count);
    defer state.allocator.free(temporary);
    @memset(temporary, @splat(0));
    for (0..inputs.complex_count) |complex| {
        for (0..inputs.population_count) |population| {
            if (!sourceEnabled(inputs.surface_litter_layer, complex, population)) continue;
            for (0..inputs.structural_fraction_count) |fraction| {
                const item = (complex * inputs.population_count + population) *
                    inputs.structural_fraction_count + fraction;
                for (0..state.humus_fraction_count) |humus| {
                    temporary[humus][0] = temporary[humus][0] +
                        inputs.humus_partition_fraction[humus] *
                            (inputs.basal_humification_carbon_g_c[item] +
                                inputs.senescence_humification_carbon_g_c[item]);
                    temporary[humus][1] = temporary[humus][1] +
                        inputs.humus_partition_fraction[humus] *
                            (inputs.basal_humification_carbon_g_c[item] +
                                inputs.senescence_humification_carbon_g_c[item]);
                    temporary[humus][2] = temporary[humus][2] +
                        inputs.humus_partition_fraction[humus] *
                            (inputs.basal_humification_nitrogen_g_n[item] +
                                inputs.senescence_humification_nitrogen_g_n[item]);
                    temporary[humus][3] = temporary[humus][3] +
                        inputs.humus_partition_fraction[humus] *
                            (inputs.basal_humification_phosphorus_g_p[item] +
                                inputs.senescence_humification_phosphorus_g_p[item]);
                    inline for (temporary[humus]) |value|
                        if (!std.math.isFinite(value) or value < 0)
                            return error.NonFiniteHumusTransferResult;
                }
            }
        }
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const field_index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, humus|
            @field(state, field.name)[humus] = values[field_index];
    };
    state.target_is_surface_soil_layer = inputs.surface_litter_layer;
}

fn sourceEnabled(
    surface_litter_layer: bool,
    zero_based_complex: usize,
    zero_based_population: usize,
) bool {
    if (zero_based_complex > 5) return false;
    if (surface_litter_layer and
        (zero_based_complex == 3 or zero_based_complex == 4)) return false;
    return zero_based_complex != 5 or
        zero_based_population <= 2 or zero_based_population == 4;
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
    if (state.humus_fraction_count != 2 or
        inputs.humus_partition_fraction.len != state.humus_fraction_count or
        inputs.complex_count == 0 or inputs.population_count == 0 or
        inputs.structural_fraction_count == 0)
        return error.InvalidHumusTransferDimensions;
    const cp = std.math.mul(usize, inputs.complex_count, inputs.population_count) catch
        return error.InvalidHumusTransferDimensions;
    const items = std.math.mul(usize, cp, inputs.structural_fraction_count) catch
        return error.InvalidHumusTransferDimensions;
    inline for (.{
        inputs.basal_humification_carbon_g_c,
        inputs.basal_humification_nitrogen_g_n,
        inputs.basal_humification_phosphorus_g_p,
        inputs.senescence_humification_carbon_g_c,
        inputs.senescence_humification_nitrogen_g_n,
        inputs.senescence_humification_phosphorus_g_p,
    }) |values| if (values.len != items)
        return error.InvalidHumusTransferDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64)
        for (@field(inputs, field.name)) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidHumusTransferInput;
}

fn fixture() Inputs {
    return .{
        .complex_count = 1,
        .population_count = 2,
        .structural_fraction_count = 2,
        .surface_litter_layer = false,
        .humus_partition_fraction = &.{ 0.6, 0.4 },
        .basal_humification_carbon_g_c = &.{ 1, 2, 3, 4 },
        .basal_humification_nitrogen_g_n = &.{ 0.1, 0.2, 0.3, 0.4 },
        .basal_humification_phosphorus_g_p = &.{ 0.05, 0.1, 0.15, 0.2 },
        .senescence_humification_carbon_g_c = &.{ 0.5, 1, 1.5, 2 },
        .senescence_humification_nitrogen_g_n = &.{ 0.05, 0.1, 0.15, 0.2 },
        .senescence_humification_phosphorus_g_p = &.{ 0.025, 0.05, 0.075, 0.1 },
    };
}

test "humic and fulvic transfers preserve K N M source accumulation" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(9, state.carbon_addition_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(6, state.carbon_addition_g_c[1], 1e-12);
    try std.testing.expectApproxEqAbs(
        state.carbon_addition_g_c[0],
        state.colonized_carbon_addition_g_c[0],
        1e-12,
    );
}

test "partitioned humification closes C N and P" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    const inputs = fixture();
    try calculate(&state, inputs);
    inline for (.{
        .{ inputs.basal_humification_carbon_g_c, inputs.senescence_humification_carbon_g_c, state.carbon_addition_g_c },
        .{ inputs.basal_humification_nitrogen_g_n, inputs.senescence_humification_nitrogen_g_n, state.nitrogen_addition_g_n },
        .{ inputs.basal_humification_phosphorus_g_p, inputs.senescence_humification_phosphorus_g_p, state.phosphorus_addition_g_p },
    }) |element| {
        var source: f64 = 0;
        var destination: f64 = 0;
        for (element[0]) |value| source += value;
        for (element[1]) |value| source += value;
        for (element[2]) |value| destination += value;
        try std.testing.expectApproxEqAbs(source, destination, 1e-12);
    }
}

test "surface and autotrophic gates exclude source-disabled items" {
    try std.testing.expect(!sourceEnabled(true, 3, 0));
    try std.testing.expect(!sourceEnabled(true, 4, 0));
    try std.testing.expect(!sourceEnabled(false, 5, 3));
    try std.testing.expect(sourceEnabled(false, 5, 4));
    try std.testing.expect(!sourceEnabled(false, 6, 0));
}

test "surface litter records redirection to surface soil layer" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.surface_litter_layer = true;
    try calculate(&state, inputs);
    try std.testing.expect(state.target_is_surface_soil_layer);
}

test "invalid partition leaves prior humus additions unchanged" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.carbon_addition_g_c[0] = 7;
    var inputs = fixture();
    inputs.humus_partition_fraction = &.{ std.math.nan(f64), 0.4 };
    try std.testing.expectError(error.InvalidHumusTransferInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.carbon_addition_g_c[0]);
}

test "NITRO 3741-3792 late overflow preserves all humus additions" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.carbon_addition_g_c[0] = 7;
    state.phosphorus_addition_g_p[1] = 11;
    var inputs = fixture();
    inputs.senescence_humification_phosphorus_g_p =
        &.{ 0.025, 0.05, 0.075, std.math.floatMax(f64) };
    inputs.basal_humification_phosphorus_g_p =
        &.{ 0.05, 0.1, 0.15, std.math.floatMax(f64) };
    try std.testing.expectError(
        error.NonFiniteHumusTransferResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.carbon_addition_g_c[0]);
    try std.testing.expectEqual(11, state.phosphorus_addition_g_p[1]);
}
