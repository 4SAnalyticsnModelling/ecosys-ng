const std = @import("std");

pub const Inputs = struct {
    fraction_count: usize,
    nonstructural_carbon_g_c: []const f64,
    nonstructural_nitrogen_g_n: []const f64,
    nonstructural_phosphorus_g_p: []const f64,
    biological_activity_fraction: []const f64,
    temperature_response: []const f64,
    water_response: []const f64,
    decay_carbon_response: []const f64,
    carbon_recycling_fraction: []const f64,
    nitrogen_recycling_fraction: []const f64,
    phosphorus_recycling_fraction: []const f64,
    humification_fraction: []const f64,
    structural_partition_fraction: []const f64,
    maximum_nitrogen_to_carbon_ratio_g_n_per_g_c: []const f64,
    maximum_phosphorus_to_carbon_ratio_g_p_per_g_c: []const f64,
    basal_decomposition_rate_h: []const f64,
    structural_carbon_g_c: []const f64,
    structural_nitrogen_g_n: []const f64,
    structural_phosphorus_g_p: []const f64,
    nonstructural_to_structural_rate_h: f64,
    humus_nitrogen_to_carbon_ratio_g_n_per_g_c: f64,
    humus_phosphorus_to_carbon_ratio_g_p_per_g_c: f64,
    negligible_nonstructural_carbon_g_c: f64,
    biochemical_time_fraction_h: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    item_count: usize,
    structural_carbon_assimilation_g_c: []f64,
    structural_nitrogen_assimilation_g_n: []f64,
    structural_phosphorus_assimilation_g_p: []f64,
    decomposition_fraction: []f64,
    decomposed_carbon_g_c: []f64,
    decomposed_nitrogen_g_n: []f64,
    decomposed_phosphorus_g_p: []f64,
    recycled_carbon_g_c: []f64,
    recycled_nitrogen_g_n: []f64,
    recycled_phosphorus_g_p: []f64,
    litterfall_carbon_g_c: []f64,
    litterfall_nitrogen_g_n: []f64,
    litterfall_phosphorus_g_p: []f64,
    humus_carbon_g_c: []f64,
    humus_nitrogen_g_n: []f64,
    humus_phosphorus_g_p: []f64,
    residue_carbon_g_c: []f64,
    residue_nitrogen_g_n: []f64,
    residue_phosphorus_g_p: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize, fraction_count: usize) !State {
        if (unit_count == 0 or fraction_count == 0)
            return error.InvalidStructuralTurnoverDimensions;
        const item_count = std.math.mul(usize, unit_count, fraction_count) catch
            return error.InvalidStructuralTurnoverDimensions;
        var state: State = undefined;
        state.allocator = allocator;
        state.item_count = item_count;
        var allocated: usize = 0;
        errdefer inline for (@typeInfo(State).@"struct".fields) |field| {
            if (field.type == []f64 and allocated > 0) {
                allocated -= 1;
                allocator.free(@field(state, field.name));
            }
        };
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

/// Exact NITRO.F 2666--2747 runtime-fraction assimilation, turnover,
/// recycling, humification, and microbial-residue allocation.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const unit_count = state.item_count / inputs.fraction_count;
    const temporary = try state.allocator.alloc([19]f64, state.item_count);
    defer state.allocator.free(temporary);
    for (0..unit_count) |unit| {
        const nonstructural_carbon = inputs.nonstructural_carbon_g_c[unit];
        const total_carbon_assimilation = inputs.biological_activity_fraction[unit] *
            inputs.nonstructural_to_structural_rate_h * @max(0, nonstructural_carbon) *
            inputs.biochemical_time_fraction_h;
        for (0..inputs.fraction_count) |fraction| {
            const item = unit * inputs.fraction_count + fraction;
            const partition = inputs.structural_partition_fraction[item];
            const carbon_assimilation = partition * total_carbon_assimilation;
            var nitrogen_assimilation: f64 = 0;
            var phosphorus_assimilation: f64 = 0;
            if (nonstructural_carbon > inputs.negligible_nonstructural_carbon_g_c) {
                nitrogen_assimilation = @min(
                    partition * @max(0, inputs.nonstructural_nitrogen_g_n[unit]),
                    carbon_assimilation * @min(
                        inputs.maximum_nitrogen_to_carbon_ratio_g_n_per_g_c[item],
                        inputs.nonstructural_nitrogen_g_n[unit] / nonstructural_carbon,
                    ),
                );
                phosphorus_assimilation = @min(
                    partition * @max(0, inputs.nonstructural_phosphorus_g_p[unit]),
                    carbon_assimilation * @min(
                        inputs.maximum_phosphorus_to_carbon_ratio_g_p_per_g_c[item],
                        inputs.nonstructural_phosphorus_g_p[unit] / nonstructural_carbon,
                    ),
                );
            }
            const decomposition_fraction = @sqrt(inputs.temperature_response[unit]) *
                inputs.water_response[unit] * inputs.basal_decomposition_rate_h[item] *
                inputs.decay_carbon_response[unit] * inputs.biochemical_time_fraction_h;
            const decomposed_carbon = @max(0, inputs.structural_carbon_g_c[item] * decomposition_fraction);
            const decomposed_nitrogen = @max(0, inputs.structural_nitrogen_g_n[item] * decomposition_fraction);
            const decomposed_phosphorus = @max(0, inputs.structural_phosphorus_g_p[item] * decomposition_fraction);
            const recycled_carbon = decomposed_carbon *
                inputs.carbon_recycling_fraction[unit];
            const recycled_nitrogen = decomposed_nitrogen *
                (inputs.nitrogen_recycling_fraction[unit] +
                    (1 - inputs.nitrogen_recycling_fraction[unit]) *
                        inputs.carbon_recycling_fraction[unit]);
            const recycled_phosphorus = decomposed_phosphorus *
                (inputs.phosphorus_recycling_fraction[unit] +
                    (1 - inputs.phosphorus_recycling_fraction[unit]) *
                        inputs.carbon_recycling_fraction[unit]);
            const litter_carbon = decomposed_carbon - recycled_carbon;
            const litter_nitrogen = decomposed_nitrogen - recycled_nitrogen;
            const litter_phosphorus = decomposed_phosphorus - recycled_phosphorus;
            const humus_carbon = @max(0, litter_carbon * inputs.humification_fraction[unit]);
            const humus_nitrogen = @max(0, @min(
                litter_nitrogen * inputs.humification_fraction[unit],
                humus_carbon * inputs.humus_nitrogen_to_carbon_ratio_g_n_per_g_c,
            ));
            const humus_phosphorus = @max(0, @min(
                litter_phosphorus * inputs.humification_fraction[unit],
                humus_carbon * inputs.humus_phosphorus_to_carbon_ratio_g_p_per_g_c,
            ));
            temporary[item] = .{
                carbon_assimilation,                  nitrogen_assimilation,        phosphorus_assimilation,
                decomposition_fraction,               decomposed_carbon,            decomposed_nitrogen,
                decomposed_phosphorus,                recycled_carbon,              recycled_nitrogen,
                recycled_phosphorus,                  litter_carbon,                litter_nitrogen,
                litter_phosphorus,                    humus_carbon,                 humus_nitrogen,
                humus_phosphorus,                     litter_carbon - humus_carbon, litter_nitrogen - humus_nitrogen,
                litter_phosphorus - humus_phosphorus,
            };
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
    if (inputs.fraction_count == 0 or state.item_count % inputs.fraction_count != 0)
        return error.InvalidStructuralTurnoverDimensions;
    const units = state.item_count / inputs.fraction_count;
    inline for (.{
        inputs.nonstructural_carbon_g_c,     inputs.nonstructural_nitrogen_g_n,
        inputs.nonstructural_phosphorus_g_p, inputs.biological_activity_fraction,
        inputs.temperature_response,         inputs.water_response,
        inputs.decay_carbon_response,        inputs.carbon_recycling_fraction,
        inputs.nitrogen_recycling_fraction,  inputs.phosphorus_recycling_fraction,
        inputs.humification_fraction,
    }) |values| {
        if (values.len != units) return error.InvalidStructuralTurnoverDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidStructuralTurnoverInput;
    }
    inline for (.{
        inputs.structural_partition_fraction,
        inputs.maximum_nitrogen_to_carbon_ratio_g_n_per_g_c,
        inputs.maximum_phosphorus_to_carbon_ratio_g_p_per_g_c,
        inputs.basal_decomposition_rate_h,
        inputs.structural_carbon_g_c,
        inputs.structural_nitrogen_g_n,
        inputs.structural_phosphorus_g_p,
    }) |values| {
        if (values.len != state.item_count) return error.InvalidStructuralTurnoverDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidStructuralTurnoverInput;
    }
    inline for (.{
        inputs.carbon_recycling_fraction,     inputs.nitrogen_recycling_fraction,
        inputs.phosphorus_recycling_fraction, inputs.humification_fraction,
    }) |values| for (values) |value| if (value > 1)
        return error.InvalidStructuralTurnoverInput;
    inline for (.{
        inputs.nonstructural_to_structural_rate_h,
        inputs.humus_nitrogen_to_carbon_ratio_g_n_per_g_c,
        inputs.humus_phosphorus_to_carbon_ratio_g_p_per_g_c,
        inputs.negligible_nonstructural_carbon_g_c,
        inputs.biochemical_time_fraction_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidStructuralTurnoverInput;
}

fn fixture() Inputs {
    return .{
        .fraction_count = 2,
        .nonstructural_carbon_g_c = &.{10},
        .nonstructural_nitrogen_g_n = &.{2},
        .nonstructural_phosphorus_g_p = &.{1},
        .biological_activity_fraction = &.{0.5},
        .temperature_response = &.{0.64},
        .water_response = &.{0.5},
        .decay_carbon_response = &.{1},
        .carbon_recycling_fraction = &.{0.2},
        .nitrogen_recycling_fraction = &.{0.5},
        .phosphorus_recycling_fraction = &.{0.25},
        .humification_fraction = &.{0.5},
        .structural_partition_fraction = &.{ 0.6, 0.4 },
        .maximum_nitrogen_to_carbon_ratio_g_n_per_g_c = &.{ 0.2, 0.1 },
        .maximum_phosphorus_to_carbon_ratio_g_p_per_g_c = &.{ 0.1, 0.05 },
        .basal_decomposition_rate_h = &.{ 0.1, 0.05 },
        .structural_carbon_g_c = &.{ 10, 20 },
        .structural_nitrogen_g_n = &.{ 2, 2 },
        .structural_phosphorus_g_p = &.{ 1, 1 },
        .nonstructural_to_structural_rate_h = 0.2,
        .humus_nitrogen_to_carbon_ratio_g_n_per_g_c = 0.1,
        .humus_phosphorus_to_carbon_ratio_g_p_per_g_c = 0.05,
        .negligible_nonstructural_carbon_g_c = 1e-12,
        .biochemical_time_fraction_h = 1,
    };
}

test "runtime structural fractions assimilate and turn over independently" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(0.6, state.structural_carbon_assimilation_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.4, state.structural_carbon_assimilation_g_c[1], 1e-12);
    try std.testing.expect(state.decomposed_carbon_g_c[0] > 0);
    try std.testing.expect(state.humus_carbon_g_c[0] > 0);
}

test "turnover carbon partition closes exactly" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    for (0..2) |item| try std.testing.expectApproxEqAbs(
        state.decomposed_carbon_g_c[item],
        state.recycled_carbon_g_c[item] + state.humus_carbon_g_c[item] +
            state.residue_carbon_g_c[item],
        1e-12,
    );
}

test "negligible nonstructural carbon disables nutrient assimilation" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.nonstructural_carbon_g_c = &.{0};
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.structural_nitrogen_assimilation_g_n[0]);
    try std.testing.expectEqual(0, state.structural_phosphorus_assimilation_g_p[0]);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.decomposed_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.water_response = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidStructuralTurnoverInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.decomposed_carbon_g_c[0]);
}
