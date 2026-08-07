const std = @import("std");

pub const Inputs = struct {
    complex_count: usize,
    population_count: usize,
    structural_fraction_count: usize,
    humus_fraction_count: usize,
    surface_litter_layer: bool,
    humus_allocation_fraction: []const f64,
    basal_humified_carbon_g_c: []const f64,
    basal_humified_nitrogen_g_n: []const f64,
    basal_humified_phosphorus_g_p: []const f64,
    senescence_humified_carbon_g_c: []const f64,
    senescence_humified_nitrogen_g_n: []const f64,
    senescence_humified_phosphorus_g_p: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    humus_fraction_count: usize,
    humus_carbon_addition_g_c: []f64,
    colonized_humus_carbon_addition_g_c: []f64,
    humus_nitrogen_addition_g_n: []f64,
    humus_phosphorus_addition_g_p: []f64,
    target_is_surface_soil_layer: bool,

    pub fn init(allocator: std.mem.Allocator, humus_fraction_count: usize) !State {
        if (humus_fraction_count == 0) return error.InvalidHumusTransferDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.humus_fraction_count = humus_fraction_count;
        state.target_is_surface_soil_layer = false;
        var allocated: usize = 0;
        errdefer inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64 and allocated > 0) {
                allocated -= 1;
                allocator.free(@field(state, field.name));
            }
        };
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

/// Exact NITRO.F 3734--3797 microbial humification-product transfer.
pub fn calculate(state: *State, inputs: Inputs) !void {
    const source_items = try validate(state, inputs);
    const temporary = try state.allocator.alloc([4]f64, state.humus_fraction_count);
    defer state.allocator.free(temporary);
    @memset(temporary, @splat(0));
    for (0..inputs.complex_count) |complex| {
        for (0..inputs.population_count) |population| {
            if (!sourceEnabled(inputs.surface_litter_layer, complex, population)) continue;
            for (0..inputs.structural_fraction_count) |fraction| {
                const item = (complex * inputs.population_count + population) *
                    inputs.structural_fraction_count + fraction;
                const carbon = inputs.basal_humified_carbon_g_c[item] +
                    inputs.senescence_humified_carbon_g_c[item];
                const nitrogen = inputs.basal_humified_nitrogen_g_n[item] +
                    inputs.senescence_humified_nitrogen_g_n[item];
                const phosphorus = inputs.basal_humified_phosphorus_g_p[item] +
                    inputs.senescence_humified_phosphorus_g_p[item];
                for (0..inputs.humus_fraction_count) |humus_fraction| {
                    const allocation = inputs.humus_allocation_fraction[humus_fraction];
                    temporary[humus_fraction][0] += allocation * carbon;
                    temporary[humus_fraction][1] += allocation * carbon;
                    temporary[humus_fraction][2] += allocation * nitrogen;
                    temporary[humus_fraction][3] += allocation * phosphorus;
                }
            }
        }
    }
    _ = source_items;
    for (temporary, 0..) |values, fraction| {
        inline for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidHumusTransferResult;
        state.humus_carbon_addition_g_c[fraction] = values[0];
        state.colonized_humus_carbon_addition_g_c[fraction] = values[1];
        state.humus_nitrogen_addition_g_n[fraction] = values[2];
        state.humus_phosphorus_addition_g_p[fraction] = values[3];
    }
    state.target_is_surface_soil_layer = inputs.surface_litter_layer;
}

pub fn sourceEnabled(
    surface_litter_layer: bool,
    zero_based_complex: usize,
    zero_based_population: usize,
) bool {
    return (!surface_litter_layer or
        (zero_based_complex != 3 and zero_based_complex != 4)) and
        (zero_based_complex != 5 or
            zero_based_population <= 2 or zero_based_population == 4);
}

fn validate(state: *const State, inputs: Inputs) !usize {
    if (inputs.complex_count == 0 or inputs.population_count == 0 or
        inputs.structural_fraction_count == 0 or
        inputs.humus_fraction_count != state.humus_fraction_count or
        inputs.humus_allocation_fraction.len != state.humus_fraction_count)
        return error.InvalidHumusTransferDimensions;
    const cp = std.math.mul(usize, inputs.complex_count, inputs.population_count) catch
        return error.InvalidHumusTransferDimensions;
    const source_items = std.math.mul(usize, cp, inputs.structural_fraction_count) catch
        return error.InvalidHumusTransferDimensions;
    inline for (.{
        inputs.basal_humified_carbon_g_c,        inputs.basal_humified_nitrogen_g_n,
        inputs.basal_humified_phosphorus_g_p,    inputs.senescence_humified_carbon_g_c,
        inputs.senescence_humified_nitrogen_g_n, inputs.senescence_humified_phosphorus_g_p,
    }) |values| {
        if (values.len != source_items) return error.InvalidHumusTransferDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidHumusTransferInput;
    }
    for (inputs.humus_allocation_fraction) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidHumusTransferInput;
    return source_items;
}

fn fixture() Inputs {
    return .{
        .complex_count = 1,
        .population_count = 2,
        .structural_fraction_count = 2,
        .humus_fraction_count = 2,
        .surface_litter_layer = false,
        .humus_allocation_fraction = &.{ 0.6, 0.4 },
        .basal_humified_carbon_g_c = &.{ 1, 2, 3, 4 },
        .basal_humified_nitrogen_g_n = &.{ 0.1, 0.2, 0.3, 0.4 },
        .basal_humified_phosphorus_g_p = &.{ 0.05, 0.1, 0.15, 0.2 },
        .senescence_humified_carbon_g_c = &.{ 0.5, 1, 1.5, 2 },
        .senescence_humified_nitrogen_g_n = &.{ 0.05, 0.1, 0.15, 0.2 },
        .senescence_humified_phosphorus_g_p = &.{ 0.025, 0.05, 0.075, 0.1 },
    };
}

test "humus products partition across runtime humic fractions" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(9, state.humus_carbon_addition_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(6, state.humus_carbon_addition_g_c[1], 1e-12);
    try std.testing.expectApproxEqAbs(
        state.humus_carbon_addition_g_c[0],
        state.colonized_humus_carbon_addition_g_c[0],
        1e-12,
    );
}

test "surface products are redirected to surface soil target" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.surface_litter_layer = true;
    try calculate(&state, inputs);
    try std.testing.expect(state.target_is_surface_soil_layer);
}

test "source gate excludes surface POC and humus complexes" {
    try std.testing.expect(!sourceEnabled(true, 3, 0));
    try std.testing.expect(!sourceEnabled(true, 4, 0));
    try std.testing.expect(sourceEnabled(false, 4, 0));
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.humus_carbon_addition_g_c[0] = 7;
    var inputs = fixture();
    inputs.humus_allocation_fraction = &.{ std.math.nan(f64), 0.4 };
    try std.testing.expectError(error.InvalidHumusTransferInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.humus_carbon_addition_g_c[0]);
}
