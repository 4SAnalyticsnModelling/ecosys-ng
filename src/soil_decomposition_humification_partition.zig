const std = @import("std");

pub const SubstrateRole = enum {
    protein,
    carbohydrate,
    cellulose,
    lignin,
    non_humifying,
};

pub const Inputs = struct {
    substrate_class_count: usize,
    substrate_roles: []const SubstrateRole,
    decomposed_carbon_g_c: []const f64,
    decomposed_nitrogen_g_n: []const f64,
    decomposed_phosphorus_g_p: []const f64,
    particulate_nitrogen_to_carbon_ratio_g_n_per_g_c: f64,
    particulate_phosphorus_to_carbon_ratio_g_p_per_g_c: f64,
    companion_humification_fraction_of_lignin: f64 = 0.10,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    item_count: usize,
    humified_carbon_g_c: []f64,
    humified_nitrogen_g_n: []f64,
    humified_phosphorus_g_p: []f64,
    dissolved_carbon_g_c: []f64,
    dissolved_nitrogen_g_n: []f64,
    dissolved_phosphorus_g_p: []f64,

    pub fn init(allocator: std.mem.Allocator, complex_count: usize, substrate_class_count: usize) !State {
        if (complex_count == 0 or substrate_class_count == 0)
            return error.InvalidHumificationPartitionDimensions;
        const item_count = std.math.mul(usize, complex_count, substrate_class_count) catch
            return error.InvalidHumificationPartitionDimensions;
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

/// Exact NITRO.F 3300--3337 humification of residue lignin and companion
/// substrate products, with all remaining products sent to dissolved pools.
/// Traversal is complex-major and structural-fraction-minor.
pub fn calculate(state: *State, inputs: Inputs) !void {
    const roles = try validate(state, inputs);
    const complex_count = state.item_count / inputs.substrate_class_count;
    const temporary = try state.allocator.alloc([6]f64, state.item_count);
    defer state.allocator.free(temporary);
    @memset(temporary, @splat(0));
    for (0..complex_count) |complex| {
        const base = complex * inputs.substrate_class_count;
        if (complex <= 2) {
            const lignin = base + roles.lignin;
            const lignin_humified = constrainedCarbon(
                inputs.decomposed_carbon_g_c[lignin],
                inputs.decomposed_nitrogen_g_n[lignin],
                inputs.decomposed_phosphorus_g_p[lignin],
                inputs.particulate_nitrogen_to_carbon_ratio_g_n_per_g_c,
                inputs.particulate_phosphorus_to_carbon_ratio_g_p_per_g_c,
                std.math.inf(f64),
            );
            temporary[lignin][0] = lignin_humified;
            const companion_cap = inputs.companion_humification_fraction_of_lignin *
                lignin_humified;
            if (!std.math.isFinite(lignin_humified) or !std.math.isFinite(companion_cap))
                return error.NonFiniteHumificationPartitionResult;
            const protein = base + roles.protein;
            temporary[protein][0] = constrainedCarbon(
                inputs.decomposed_carbon_g_c[protein],
                inputs.decomposed_nitrogen_g_n[protein],
                inputs.decomposed_phosphorus_g_p[protein],
                inputs.particulate_nitrogen_to_carbon_ratio_g_n_per_g_c,
                inputs.particulate_phosphorus_to_carbon_ratio_g_p_per_g_c,
                companion_cap,
            );
            const carbohydrate = base + roles.carbohydrate;
            temporary[carbohydrate][0] = constrainedCarbon(
                inputs.decomposed_carbon_g_c[carbohydrate],
                inputs.decomposed_nitrogen_g_n[carbohydrate],
                inputs.decomposed_phosphorus_g_p[carbohydrate],
                inputs.particulate_nitrogen_to_carbon_ratio_g_n_per_g_c,
                inputs.particulate_phosphorus_to_carbon_ratio_g_p_per_g_c,
                companion_cap,
            );
            const cellulose = base + roles.cellulose;
            temporary[cellulose][0] = constrainedCarbon(
                inputs.decomposed_carbon_g_c[cellulose],
                inputs.decomposed_nitrogen_g_n[cellulose],
                inputs.decomposed_phosphorus_g_p[cellulose],
                inputs.particulate_nitrogen_to_carbon_ratio_g_n_per_g_c,
                inputs.particulate_phosphorus_to_carbon_ratio_g_p_per_g_c,
                @max(0, companion_cap - temporary[carbohydrate][0]),
            );
        }
        for (0..inputs.substrate_class_count) |class| {
            const item = base + class;
            const humified_carbon = temporary[item][0];
            const humified_nitrogen = @min(inputs.decomposed_nitrogen_g_n[item], humified_carbon * inputs.particulate_nitrogen_to_carbon_ratio_g_n_per_g_c);
            const humified_phosphorus = @min(inputs.decomposed_phosphorus_g_p[item], humified_carbon * inputs.particulate_phosphorus_to_carbon_ratio_g_p_per_g_c);
            temporary[item] = .{
                humified_carbon,                                      humified_nitrogen,                                        humified_phosphorus,
                inputs.decomposed_carbon_g_c[item] - humified_carbon, inputs.decomposed_nitrogen_g_n[item] - humified_nitrogen, inputs.decomposed_phosphorus_g_p[item] - humified_phosphorus,
            };
            for (temporary[item]) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteHumificationPartitionResult;
        }
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, item| @field(state, field.name)[item] = values[index];
    };
}

const RoleIndices = struct {
    protein: usize,
    carbohydrate: usize,
    cellulose: usize,
    lignin: usize,
};

fn constrainedCarbon(
    carbon: f64,
    nitrogen: f64,
    phosphorus: f64,
    nitrogen_ratio: f64,
    phosphorus_ratio: f64,
    ceiling: f64,
) f64 {
    return @max(0, @min(carbon, @min(nitrogen / nitrogen_ratio, @min(phosphorus / phosphorus_ratio, ceiling))));
}

fn stateFieldIndex(comptime name: []const u8) usize {
    var index: usize = 0;
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        if (comptime std.mem.eql(u8, field.name, name)) return index;
        index += 1;
    };
    unreachable;
}

fn validate(state: *const State, inputs: Inputs) !RoleIndices {
    if (inputs.substrate_class_count == 0 or state.item_count % inputs.substrate_class_count != 0 or
        inputs.substrate_roles.len != inputs.substrate_class_count)
        return error.InvalidHumificationPartitionDimensions;
    inline for (.{
        inputs.decomposed_carbon_g_c,     inputs.decomposed_nitrogen_g_n,
        inputs.decomposed_phosphorus_g_p,
    }) |values| {
        if (values.len != state.item_count) return error.InvalidHumificationPartitionDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidHumificationPartitionInput;
    }
    inline for (.{
        inputs.particulate_nitrogen_to_carbon_ratio_g_n_per_g_c,
        inputs.particulate_phosphorus_to_carbon_ratio_g_p_per_g_c,
        inputs.companion_humification_fraction_of_lignin,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidHumificationPartitionInput;
    if (inputs.particulate_nitrogen_to_carbon_ratio_g_n_per_g_c == 0 or
        inputs.particulate_phosphorus_to_carbon_ratio_g_p_per_g_c == 0)
        return error.InvalidHumificationPartitionInput;
    var protein: ?usize = null;
    var carbohydrate: ?usize = null;
    var cellulose: ?usize = null;
    var lignin: ?usize = null;
    for (inputs.substrate_roles, 0..) |role, index| switch (role) {
        .protein => if (protein != null) return error.DuplicateHumificationSubstrateRole else {
            protein = index;
        },
        .carbohydrate => if (carbohydrate != null) return error.DuplicateHumificationSubstrateRole else {
            carbohydrate = index;
        },
        .cellulose => if (cellulose != null) return error.DuplicateHumificationSubstrateRole else {
            cellulose = index;
        },
        .lignin => if (lignin != null) return error.DuplicateHumificationSubstrateRole else {
            lignin = index;
        },
        .non_humifying => {},
    };
    return .{
        .protein = protein orelse return error.MissingHumificationSubstrateRole,
        .carbohydrate = carbohydrate orelse return error.MissingHumificationSubstrateRole,
        .cellulose = cellulose orelse return error.MissingHumificationSubstrateRole,
        .lignin = lignin orelse return error.MissingHumificationSubstrateRole,
    };
}

fn fixture() Inputs {
    return .{
        .substrate_class_count = 5,
        .substrate_roles = &.{ .protein, .carbohydrate, .cellulose, .lignin, .non_humifying },
        .decomposed_carbon_g_c = &.{ 2, 2, 2, 10, 2 },
        .decomposed_nitrogen_g_n = &.{ 1, 1, 1, 2, 1 },
        .decomposed_phosphorus_g_p = &.{ 0.5, 0.5, 0.5, 1, 0.5 },
        .particulate_nitrogen_to_carbon_ratio_g_n_per_g_c = 0.1,
        .particulate_phosphorus_to_carbon_ratio_g_p_per_g_c = 0.05,
    };
}

test "lignin controls companion substrate humification" {
    var state = try State.init(std.testing.allocator, 1, 5);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(10, state.humified_carbon_g_c[3], 1e-12);
    try std.testing.expectApproxEqAbs(1, state.humified_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(1, state.humified_carbon_g_c[1], 1e-12);
    try std.testing.expectEqual(0, state.humified_carbon_g_c[2]);
    try std.testing.expectEqual(0, state.humified_carbon_g_c[4]);
}

test "humified plus dissolved products close C N and P" {
    var state = try State.init(std.testing.allocator, 1, 5);
    defer state.deinit();
    const inputs = fixture();
    try calculate(&state, inputs);
    for (0..5) |item| {
        try std.testing.expectApproxEqAbs(inputs.decomposed_carbon_g_c[item], state.humified_carbon_g_c[item] + state.dissolved_carbon_g_c[item], 1e-12);
        try std.testing.expectApproxEqAbs(inputs.decomposed_nitrogen_g_n[item], state.humified_nitrogen_g_n[item] + state.dissolved_nitrogen_g_n[item], 1e-12);
        try std.testing.expectApproxEqAbs(inputs.decomposed_phosphorus_g_p[item], state.humified_phosphorus_g_p[item] + state.dissolved_phosphorus_g_p[item], 1e-12);
    }
}

test "complexes above source residue range send all products dissolved" {
    var state = try State.init(std.testing.allocator, 4, 5);
    defer state.deinit();
    var carbon: [20]f64 = @splat(2);
    var nitrogen: [20]f64 = @splat(1);
    var phosphorus: [20]f64 = @splat(0.5);
    var inputs = fixture();
    inputs.decomposed_carbon_g_c = &carbon;
    inputs.decomposed_nitrogen_g_n = &nitrogen;
    inputs.decomposed_phosphorus_g_p = &phosphorus;
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.humified_carbon_g_c[15]);
    try std.testing.expectEqual(2, state.dissolved_carbon_g_c[15]);
}

test "duplicate role fails before publication" {
    var state = try State.init(std.testing.allocator, 1, 5);
    defer state.deinit();
    state.humified_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.substrate_roles = &.{ .protein, .protein, .cellulose, .lignin, .non_humifying };
    try std.testing.expectError(error.DuplicateHumificationSubstrateRole, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.humified_carbon_g_c[0]);
}

test "NITRO K greater than two sends every C N P product dissolved" {
    var state = try State.init(std.testing.allocator, 4, 5);
    defer state.deinit();
    var carbon: [20]f64 = @splat(2);
    var nitrogen: [20]f64 = @splat(1);
    var phosphorus: [20]f64 = @splat(0.5);
    var inputs = fixture();
    inputs.decomposed_carbon_g_c = &carbon;
    inputs.decomposed_nitrogen_g_n = &nitrogen;
    inputs.decomposed_phosphorus_g_p = &phosphorus;
    try calculate(&state, inputs);
    for (15..20) |item| {
        try std.testing.expectEqual(0, state.humified_carbon_g_c[item]);
        try std.testing.expectEqual(0, state.humified_nitrogen_g_n[item]);
        try std.testing.expectEqual(0, state.humified_phosphorus_g_p[item]);
        try std.testing.expectEqual(2, state.dissolved_carbon_g_c[item]);
        try std.testing.expectEqual(1, state.dissolved_nitrogen_g_n[item]);
        try std.testing.expectEqual(0.5, state.dissolved_phosphorus_g_p[item]);
    }
}

test "NITRO 3300-3337 late derived overflow preserves all products" {
    var state = try State.init(std.testing.allocator, 1, 5);
    defer state.deinit();
    state.humified_carbon_g_c[0] = 7;
    state.dissolved_phosphorus_g_p[4] = 11;
    var inputs = fixture();
    inputs.companion_humification_fraction_of_lignin =
        std.math.floatMax(f64);
    try std.testing.expectError(
        error.NonFiniteHumificationPartitionResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.humified_carbon_g_c[0]);
    try std.testing.expectEqual(11, state.dissolved_phosphorus_g_p[4]);
}
