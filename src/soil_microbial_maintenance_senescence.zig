const std = @import("std");

pub const Inputs = struct {
    fraction_count: usize,
    senescence_respiration_g_c: []const f64,
    total_maintenance_respiration_g_c: []const f64,
    carbon_recycling_fraction: []const f64,
    nitrogen_recycling_fraction: []const f64,
    phosphorus_recycling_fraction: []const f64,
    humification_fraction: []const f64,
    active_nitrogen_to_carbon_ratio_g_n_per_g_c: []const f64,
    active_phosphorus_to_carbon_ratio_g_p_per_g_c: []const f64,
    fraction_maintenance_respiration_g_c: []const f64,
    structural_carbon_g_c: []const f64,
    structural_nitrogen_g_n: []const f64,
    structural_phosphorus_g_p: []const f64,
    humus_nitrogen_to_carbon_ratio_g_n_per_g_c: f64,
    humus_phosphorus_to_carbon_ratio_g_p_per_g_c: f64,
    negligible_respiration_g_c: f64,
    negligible_recycling_fraction: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    item_count: usize,
    maintenance_deficit_fraction: []f64,
    senesced_carbon_g_c: []f64,
    senesced_nitrogen_g_n: []f64,
    senesced_phosphorus_g_p: []f64,
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
            return error.InvalidMaintenanceSenescenceDimensions;
        const item_count = std.math.mul(usize, unit_count, fraction_count) catch
            return error.InvalidMaintenanceSenescenceDimensions;
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

/// Exact NITRO.F 2751--2832 maintenance-deficit microbial senescence.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const units = state.item_count / inputs.fraction_count;
    const temporary = try state.allocator.alloc([16]f64, state.item_count);
    defer state.allocator.free(temporary);
    for (0..units) |unit| {
        const active = inputs.senescence_respiration_g_c[unit] >
            inputs.negligible_respiration_g_c and
            inputs.total_maintenance_respiration_g_c[unit] >
                inputs.negligible_respiration_g_c and
            inputs.carbon_recycling_fraction[unit] >
                inputs.negligible_recycling_fraction;
        const deficit_fraction = if (active)
            inputs.senescence_respiration_g_c[unit] /
                inputs.total_maintenance_respiration_g_c[unit]
        else
            0;
        for (0..inputs.fraction_count) |fraction| {
            const item = unit * inputs.fraction_count + fraction;
            if (!active) {
                temporary[item] = @splat(0);
                continue;
            }
            const carbon_demand = deficit_fraction *
                inputs.fraction_maintenance_respiration_g_c[item] /
                inputs.carbon_recycling_fraction[unit];
            const senesced_carbon = @min(
                inputs.structural_carbon_g_c[item],
                @max(0, carbon_demand),
            );
            const nitrogen_demand = senesced_carbon *
                inputs.active_nitrogen_to_carbon_ratio_g_n_per_g_c[unit];
            const senesced_nitrogen = @min(
                inputs.structural_nitrogen_g_n[item],
                @max(0, nitrogen_demand),
            );
            const phosphorus_demand = senesced_carbon *
                inputs.active_phosphorus_to_carbon_ratio_g_p_per_g_c[unit];
            const senesced_phosphorus = @min(
                inputs.structural_phosphorus_g_p[item],
                @max(0, phosphorus_demand),
            );
            const recycled_carbon = senesced_carbon *
                inputs.carbon_recycling_fraction[unit];
            const recycled_nitrogen = senesced_nitrogen *
                (inputs.nitrogen_recycling_fraction[unit] +
                    (1 - inputs.nitrogen_recycling_fraction[unit]) *
                        inputs.carbon_recycling_fraction[unit]);
            const recycled_phosphorus = senesced_phosphorus *
                (inputs.phosphorus_recycling_fraction[unit] +
                    (1 - inputs.phosphorus_recycling_fraction[unit]) *
                        inputs.carbon_recycling_fraction[unit]);
            const litter_carbon = senesced_carbon - recycled_carbon;
            const litter_nitrogen = senesced_nitrogen - recycled_nitrogen;
            const litter_phosphorus = senesced_phosphorus - recycled_phosphorus;
            const humus_carbon_candidate =
                litter_carbon * inputs.humification_fraction[unit];
            const humus_carbon = @max(0, humus_carbon_candidate);
            // Source caps humus nutrients against total litterfall C, not humified C.
            const humus_nitrogen_by_fraction =
                litter_nitrogen * inputs.humification_fraction[unit];
            const humus_nitrogen_by_ratio = litter_carbon *
                inputs.humus_nitrogen_to_carbon_ratio_g_n_per_g_c;
            const humus_nitrogen = @max(
                0,
                @min(humus_nitrogen_by_fraction, humus_nitrogen_by_ratio),
            );
            const humus_phosphorus_by_fraction =
                litter_phosphorus * inputs.humification_fraction[unit];
            const humus_phosphorus_by_ratio = litter_carbon *
                inputs.humus_phosphorus_to_carbon_ratio_g_p_per_g_c;
            const humus_phosphorus = @max(
                0,
                @min(humus_phosphorus_by_fraction, humus_phosphorus_by_ratio),
            );
            inline for (.{
                carbon_demand,
                nitrogen_demand,
                phosphorus_demand,
                humus_carbon_candidate,
                humus_nitrogen_by_fraction,
                humus_nitrogen_by_ratio,
                humus_phosphorus_by_fraction,
                humus_phosphorus_by_ratio,
            }) |value| if (!std.math.isFinite(value))
                return error.NonFiniteMaintenanceSenescenceResult;
            temporary[item] = .{
                deficit_fraction, senesced_carbon,              senesced_nitrogen,                senesced_phosphorus,
                recycled_carbon,  recycled_nitrogen,            recycled_phosphorus,              litter_carbon,
                litter_nitrogen,  litter_phosphorus,            humus_carbon,                     humus_nitrogen,
                humus_phosphorus, litter_carbon - humus_carbon, litter_nitrogen - humus_nitrogen, litter_phosphorus - humus_phosphorus,
            };
            for (temporary[item]) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteMaintenanceSenescenceResult;
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
        return error.InvalidMaintenanceSenescenceDimensions;
    const units = state.item_count / inputs.fraction_count;
    inline for (.{
        inputs.senescence_respiration_g_c,                  inputs.total_maintenance_respiration_g_c,
        inputs.carbon_recycling_fraction,                   inputs.nitrogen_recycling_fraction,
        inputs.phosphorus_recycling_fraction,               inputs.humification_fraction,
        inputs.active_nitrogen_to_carbon_ratio_g_n_per_g_c, inputs.active_phosphorus_to_carbon_ratio_g_p_per_g_c,
    }) |values| {
        if (values.len != units) return error.InvalidMaintenanceSenescenceDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidMaintenanceSenescenceInput;
    }
    inline for (.{
        inputs.fraction_maintenance_respiration_g_c, inputs.structural_carbon_g_c,
        inputs.structural_nitrogen_g_n,              inputs.structural_phosphorus_g_p,
    }) |values| {
        if (values.len != state.item_count) return error.InvalidMaintenanceSenescenceDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidMaintenanceSenescenceInput;
    }
    inline for (.{
        inputs.carbon_recycling_fraction,     inputs.nitrogen_recycling_fraction,
        inputs.phosphorus_recycling_fraction, inputs.humification_fraction,
    }) |values| for (values) |value| if (value > 1)
        return error.InvalidMaintenanceSenescenceInput;
    inline for (.{
        inputs.humus_nitrogen_to_carbon_ratio_g_n_per_g_c,
        inputs.humus_phosphorus_to_carbon_ratio_g_p_per_g_c,
        inputs.negligible_respiration_g_c,
        inputs.negligible_recycling_fraction,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidMaintenanceSenescenceInput;
}

fn fixture() Inputs {
    return .{
        .fraction_count = 2,
        .senescence_respiration_g_c = &.{1},
        .total_maintenance_respiration_g_c = &.{4},
        .carbon_recycling_fraction = &.{0.5},
        .nitrogen_recycling_fraction = &.{0.25},
        .phosphorus_recycling_fraction = &.{0.2},
        .humification_fraction = &.{0.5},
        .active_nitrogen_to_carbon_ratio_g_n_per_g_c = &.{0.1},
        .active_phosphorus_to_carbon_ratio_g_p_per_g_c = &.{0.05},
        .fraction_maintenance_respiration_g_c = &.{ 2, 2 },
        .structural_carbon_g_c = &.{ 10, 10 },
        .structural_nitrogen_g_n = &.{ 2, 2 },
        .structural_phosphorus_g_p = &.{ 1, 1 },
        .humus_nitrogen_to_carbon_ratio_g_n_per_g_c = 0.1,
        .humus_phosphorus_to_carbon_ratio_g_p_per_g_c = 0.05,
        .negligible_respiration_g_c = 1e-12,
        .negligible_recycling_fraction = 1e-12,
    };
}

test "maintenance deficit drives runtime-fraction senescence" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(0.25, state.maintenance_deficit_fraction[0], 1e-12);
    try std.testing.expectApproxEqAbs(1, state.senesced_carbon_g_c[0], 1e-12);
    try std.testing.expect(state.recycled_carbon_g_c[0] > 0);
}

test "senescence carbon partition closes exactly" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    for (0..2) |item| try std.testing.expectApproxEqAbs(
        state.senesced_carbon_g_c[item],
        state.recycled_carbon_g_c[item] + state.humus_carbon_g_c[item] +
            state.residue_carbon_g_c[item],
        1e-12,
    );
}

test "no maintenance deficit zeros all senescence fluxes" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.senescence_respiration_g_c = &.{0};
    try calculate(&state, inputs);
    try std.testing.expectEqual(0, state.senesced_carbon_g_c[0]);
}

test "NITRO ZERO threshold explicitly zeros all senescence outputs" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    var inputs = fixture();
    inputs.negligible_recycling_fraction = 0.5;
    try calculate(&state, inputs);
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64)
        for (@field(state, field.name)) |value| try std.testing.expectEqual(0, value);
}

test "invalid input leaves state unchanged" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.senesced_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.humification_fraction = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidMaintenanceSenescenceInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.senesced_carbon_g_c[0]);
}

test "NITRO 2751-2832 pre-cap overflow preserves senescence state" {
    var state = try State.init(std.testing.allocator, 1, 2);
    defer state.deinit();
    state.senesced_carbon_g_c[0] = 7;
    var inputs = fixture();
    inputs.active_nitrogen_to_carbon_ratio_g_n_per_g_c = &.{std.math.floatMax(f64)};
    inputs.structural_carbon_g_c = &.{ std.math.floatMax(f64), 10 };
    inputs.fraction_maintenance_respiration_g_c = &.{ std.math.floatMax(f64), 2 };
    try std.testing.expectError(
        error.NonFiniteMaintenanceSenescenceResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.senesced_carbon_g_c[0]);
}
