const std = @import("std");

pub const Inputs = struct {
    complex_count: usize,
    population_count: usize,
    structural_fraction_count: usize,
    surface_litter_layer: bool,
    structural_carbon_g_c: []const f64,
    structural_nitrogen_g_n: []const f64,
    structural_phosphorus_g_p: []const f64,
    structural_carbon_assimilation_g_c: []const f64,
    structural_nitrogen_assimilation_g_n: []const f64,
    structural_phosphorus_assimilation_g_p: []const f64,
    basal_decomposition_carbon_g_c: []const f64,
    basal_decomposition_nitrogen_g_n: []const f64,
    basal_decomposition_phosphorus_g_p: []const f64,
    maintenance_senescence_carbon_g_c: []const f64,
    maintenance_senescence_nitrogen_g_n: []const f64,
    maintenance_senescence_phosphorus_g_p: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    item_count: usize,
    structural_carbon_g_c: []f64,
    structural_nitrogen_g_n: []f64,
    structural_phosphorus_g_p: []f64,

    pub fn init(
        allocator: std.mem.Allocator,
        complex_count: usize,
        population_count: usize,
        fraction_count: usize,
    ) !State {
        if (complex_count == 0 or population_count == 0 or fraction_count == 0)
            return error.InvalidStructuralGrowthCommitDimensions;
        const cp = std.math.mul(usize, complex_count, population_count) catch
            return error.InvalidStructuralGrowthCommitDimensions;
        const item_count = std.math.mul(usize, cp, fraction_count) catch
            return error.InvalidStructuralGrowthCommitDimensions;
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

/// Exact NITRO.F 3704--3724 structural microbial growth, decay, and
/// maintenance-senescence state commit.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([3]f64, state.item_count);
    defer state.allocator.free(temporary);
    for (0..inputs.complex_count) |complex| {
        for (0..inputs.population_count) |population| {
            const enabled = sourceEnabled(
                inputs.surface_litter_layer,
                complex,
                population,
            );
            for (0..inputs.structural_fraction_count) |fraction| {
                const item = (complex * inputs.population_count + population) *
                    inputs.structural_fraction_count + fraction;
                if (!enabled) {
                    temporary[item] = .{
                        inputs.structural_carbon_g_c[item],
                        inputs.structural_nitrogen_g_n[item],
                        inputs.structural_phosphorus_g_p[item],
                    };
                    continue;
                }
                temporary[item] = .{
                    inputs.structural_carbon_g_c[item] +
                        inputs.structural_carbon_assimilation_g_c[item] -
                        inputs.basal_decomposition_carbon_g_c[item] -
                        inputs.maintenance_senescence_carbon_g_c[item],
                    inputs.structural_nitrogen_g_n[item] +
                        inputs.structural_nitrogen_assimilation_g_n[item] -
                        inputs.basal_decomposition_nitrogen_g_n[item] -
                        inputs.maintenance_senescence_nitrogen_g_n[item],
                    inputs.structural_phosphorus_g_p[item] +
                        inputs.structural_phosphorus_assimilation_g_p[item] -
                        inputs.basal_decomposition_phosphorus_g_p[item] -
                        inputs.maintenance_senescence_phosphorus_g_p[item],
                };
            }
        }
    }
    for (temporary) |values| inline for (values) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidStructuralGrowthCommitResult;
    for (temporary, 0..) |values, item| {
        state.structural_carbon_g_c[item] = values[0];
        state.structural_nitrogen_g_n[item] = values[1];
        state.structural_phosphorus_g_p[item] = values[2];
    }
}

pub fn sourceEnabled(
    surface_litter_layer: bool,
    zero_based_complex: usize,
    zero_based_population: usize,
) bool {
    if (zero_based_complex > 5) return false;
    const surface_complex_enabled = !surface_litter_layer or
        (zero_based_complex != 3 and zero_based_complex != 4);
    const autotroph_population_enabled = zero_based_complex != 5 or
        zero_based_population <= 2 or zero_based_population == 4;
    return surface_complex_enabled and autotroph_population_enabled;
}

fn validate(state: *const State, inputs: Inputs) !void {
    if (inputs.complex_count == 0 or inputs.population_count == 0 or
        inputs.structural_fraction_count == 0)
        return error.InvalidStructuralGrowthCommitDimensions;
    const cp = std.math.mul(usize, inputs.complex_count, inputs.population_count) catch
        return error.InvalidStructuralGrowthCommitDimensions;
    const items = std.math.mul(usize, cp, inputs.structural_fraction_count) catch
        return error.InvalidStructuralGrowthCommitDimensions;
    if (items != state.item_count) return error.InvalidStructuralGrowthCommitDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != items) return error.InvalidStructuralGrowthCommitDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidStructuralGrowthCommitInput;
    };
}

fn fixture() Inputs {
    return .{
        .complex_count = 1,
        .population_count = 2,
        .structural_fraction_count = 2,
        .surface_litter_layer = false,
        .structural_carbon_g_c = &.{ 10, 20, 30, 40 },
        .structural_nitrogen_g_n = &.{ 2, 4, 6, 8 },
        .structural_phosphorus_g_p = &.{ 1, 2, 3, 4 },
        .structural_carbon_assimilation_g_c = &.{ 2, 2, 2, 2 },
        .structural_nitrogen_assimilation_g_n = &.{ 0.4, 0.4, 0.4, 0.4 },
        .structural_phosphorus_assimilation_g_p = &.{ 0.2, 0.2, 0.2, 0.2 },
        .basal_decomposition_carbon_g_c = &.{ 1, 1, 1, 1 },
        .basal_decomposition_nitrogen_g_n = &.{ 0.2, 0.2, 0.2, 0.2 },
        .basal_decomposition_phosphorus_g_p = &.{ 0.1, 0.1, 0.1, 0.1 },
        .maintenance_senescence_carbon_g_c = &.{ 0.5, 0.5, 0.5, 0.5 },
        .maintenance_senescence_nitrogen_g_n = &.{ 0.1, 0.1, 0.1, 0.1 },
        .maintenance_senescence_phosphorus_g_p = &.{ 0.05, 0.05, 0.05, 0.05 },
    };
}

test "structural commit applies growth minus both turnover pathways" {
    var state = try State.init(std.testing.allocator, 1, 2, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(10.5, state.structural_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(2.1, state.structural_nitrogen_g_n[0], 1e-12);
}

test "source execution gates reproduce surface and autotroph exclusions" {
    try std.testing.expect(!sourceEnabled(true, 3, 0));
    try std.testing.expect(!sourceEnabled(true, 4, 0));
    try std.testing.expect(sourceEnabled(false, 3, 0));
    try std.testing.expect(!sourceEnabled(false, 5, 3));
    try std.testing.expect(sourceEnabled(false, 5, 4));
    try std.testing.expect(!sourceEnabled(false, 6, 0));
}

test "disabled surface complexes retain original pools" {
    var state = try State.init(std.testing.allocator, 4, 1, 1);
    defer state.deinit();
    var inputs = fixture();
    inputs.complex_count = 4;
    inputs.population_count = 1;
    inputs.structural_fraction_count = 1;
    inputs.surface_litter_layer = true;
    inputs.structural_carbon_g_c = &.{ 10, 20, 30, 40 };
    inputs.structural_nitrogen_g_n = &.{ 2, 4, 6, 8 };
    inputs.structural_phosphorus_g_p = &.{ 1, 2, 3, 4 };
    inputs.structural_carbon_assimilation_g_c = &.{ 2, 2, 2, 2 };
    inputs.structural_nitrogen_assimilation_g_n = &.{ 0.4, 0.4, 0.4, 0.4 };
    inputs.structural_phosphorus_assimilation_g_p = &.{ 0.2, 0.2, 0.2, 0.2 };
    inputs.basal_decomposition_carbon_g_c = &.{ 1, 1, 1, 1 };
    inputs.basal_decomposition_nitrogen_g_n = &.{ 0.2, 0.2, 0.2, 0.2 };
    inputs.basal_decomposition_phosphorus_g_p = &.{ 0.1, 0.1, 0.1, 0.1 };
    inputs.maintenance_senescence_carbon_g_c = &.{ 0.5, 0.5, 0.5, 0.5 };
    inputs.maintenance_senescence_nitrogen_g_n = &.{ 0.1, 0.1, 0.1, 0.1 };
    inputs.maintenance_senescence_phosphorus_g_p = &.{ 0.05, 0.05, 0.05, 0.05 };
    try calculate(&state, inputs);
    try std.testing.expectEqual(40, state.structural_carbon_g_c[3]);
}

test "turnover overdraw fails atomically" {
    var state = try State.init(std.testing.allocator, 1, 2, 2);
    defer state.deinit();
    state.structural_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.basal_decomposition_carbon_g_c = &.{ 20, 1, 1, 1 };
    try std.testing.expectError(error.InvalidStructuralGrowthCommitResult, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.structural_carbon_g_c[0]);
}

test "structural commit preserves exact C N P source balance" {
    var state = try State.init(std.testing.allocator, 1, 2, 2);
    defer state.deinit();
    const inputs = fixture();
    try calculate(&state, inputs);
    inline for (.{
        .{
            inputs.structural_carbon_g_c,
            inputs.structural_carbon_assimilation_g_c,
            inputs.basal_decomposition_carbon_g_c,
            inputs.maintenance_senescence_carbon_g_c,
            state.structural_carbon_g_c,
        },
        .{
            inputs.structural_nitrogen_g_n,
            inputs.structural_nitrogen_assimilation_g_n,
            inputs.basal_decomposition_nitrogen_g_n,
            inputs.maintenance_senescence_nitrogen_g_n,
            state.structural_nitrogen_g_n,
        },
        .{
            inputs.structural_phosphorus_g_p,
            inputs.structural_phosphorus_assimilation_g_p,
            inputs.basal_decomposition_phosphorus_g_p,
            inputs.maintenance_senescence_phosphorus_g_p,
            state.structural_phosphorus_g_p,
        },
    }) |element| for (0..state.item_count) |item|
        try std.testing.expectApproxEqAbs(
            element[0][item] + element[1][item] -
                element[2][item] - element[3][item],
            element[4][item],
            1e-12,
        );
}

test "NITRO 3704-3724 late derived overflow fails atomically" {
    var state = try State.init(std.testing.allocator, 1, 2, 2);
    defer state.deinit();
    state.structural_carbon_g_c[0] = 7;
    state.structural_phosphorus_g_p[3] = 11;
    var inputs = fixture();
    inputs.structural_phosphorus_g_p =
        &.{ 1, 2, 3, std.math.floatMax(f64) };
    inputs.structural_phosphorus_assimilation_g_p =
        &.{ 0.2, 0.2, 0.2, std.math.floatMax(f64) };
    try std.testing.expectError(
        error.InvalidStructuralGrowthCommitResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.structural_carbon_g_c[0]);
    try std.testing.expectEqual(11, state.structural_phosphorus_g_p[3]);
}
