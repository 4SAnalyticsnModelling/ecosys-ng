const std = @import("std");

pub const Inputs = struct {
    heterotrophic_complex: []const bool,
    maintenance_respiration_g_c: []const f64,
    aerobic_respiration_g_c: []const f64,
    growth_respiration_g_c: []const f64,
    denitrification_respiration_g_c: []const f64,
    nitrogen_fixation_respiration_g_c: []const f64,
    aerobic_growth_efficiency: []const f64,
    denitrification_growth_efficiency: []const f64,
    dissolved_organic_carbon_fraction: []const f64,
    acetate_fraction: []const f64,
    dissolved_organic_nitrogen_g_n: []const f64,
    dissolved_organic_phosphorus_g_p: []const f64,
    active_biomass_fraction: []const f64,
    dissolved_nitrogen_to_carbon_ratio_g_n_per_g_c: []const f64,
    dissolved_phosphorus_to_carbon_ratio_g_p_per_g_c: []const f64,
    nitrogen_limitation_fraction: []const f64,
    phosphorus_limitation_fraction: []const f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    aerobic_carbon_uptake_g_c: []f64,
    denitrification_carbon_uptake_g_c: []f64,
    total_carbon_uptake_g_c: []f64,
    primary_carbon_substrate_uptake_g_c: []f64,
    acetate_uptake_g_c: []f64,
    dissolved_organic_nitrogen_uptake_g_n: []f64,
    dissolved_organic_phosphorus_uptake_g_p: []f64,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidSubstrateUptakeDimensions;
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

/// Exact NITRO.F 2555--2614 respiration-driven C, N, and P uptake.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(state, inputs);
    const temporary = try state.allocator.alloc([7]f64, state.unit_count);
    defer state.allocator.free(temporary);
    for (0..state.unit_count) |unit| {
        if (inputs.nitrogen_fixation_respiration_g_c[unit] >
            inputs.growth_respiration_g_c[unit])
            return error.InconsistentSubstrateUptakeRespirationLedger;
        const aerobic_carbon =
            @min(inputs.maintenance_respiration_g_c[unit], inputs.aerobic_respiration_g_c[unit]) +
            inputs.nitrogen_fixation_respiration_g_c[unit] +
            (inputs.growth_respiration_g_c[unit] -
                inputs.nitrogen_fixation_respiration_g_c[unit]) /
                inputs.aerobic_growth_efficiency[unit];
        const denitrification_carbon =
            inputs.denitrification_respiration_g_c[unit] /
            inputs.denitrification_growth_efficiency[unit];
        const total = aerobic_carbon + denitrification_carbon;
        var dissolved_carbon = total;
        var acetate: f64 = 0;
        var dissolved_nitrogen: f64 = 0;
        var dissolved_phosphorus: f64 = 0;
        if (inputs.heterotrophic_complex[unit]) {
            dissolved_carbon = aerobic_carbon *
                inputs.dissolved_organic_carbon_fraction[unit] +
                denitrification_carbon;
            acetate = aerobic_carbon * inputs.acetate_fraction[unit];
            const organic_carbon = dissolved_carbon + acetate;
            dissolved_nitrogen = @max(0, @min(
                inputs.dissolved_organic_nitrogen_g_n[unit] *
                    inputs.active_biomass_fraction[unit],
                organic_carbon *
                    inputs.dissolved_nitrogen_to_carbon_ratio_g_n_per_g_c[unit] /
                    inputs.nitrogen_limitation_fraction[unit],
            ));
            dissolved_phosphorus = @max(0, @min(
                inputs.dissolved_organic_phosphorus_g_p[unit] *
                    inputs.active_biomass_fraction[unit],
                organic_carbon *
                    inputs.dissolved_phosphorus_to_carbon_ratio_g_p_per_g_c[unit] /
                    inputs.phosphorus_limitation_fraction[unit],
            ));
        }
        temporary[unit] = .{
            aerobic_carbon, denitrification_carbon, total,                dissolved_carbon,
            acetate,        dissolved_nitrogen,     dissolved_phosphorus,
        };
        for (temporary[unit]) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.NonFiniteSubstrateUptakeResult;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        const index = comptime stateFieldIndex(field.name);
        for (temporary, 0..) |values, unit| @field(state, field.name)[unit] = values[index];
    };
}

pub fn sourceHeterotrophic(zero_based_complex: usize) bool {
    return zero_based_complex <= 3;
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
    if (inputs.heterotrophic_complex.len != n) return error.InvalidSubstrateUptakeDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| if (field.type == []const f64) {
        const values = @field(inputs, field.name);
        if (values.len != n) return error.InvalidSubstrateUptakeDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidSubstrateUptakeInput;
    };
    inline for (.{
        inputs.aerobic_growth_efficiency,
        inputs.denitrification_growth_efficiency,
        inputs.nitrogen_limitation_fraction,
        inputs.phosphorus_limitation_fraction,
    }) |values| for (values) |value| if (value == 0)
        return error.InvalidSubstrateUptakeInput;
}

fn fixture() Inputs {
    return .{
        .heterotrophic_complex = &.{ true, false },
        .maintenance_respiration_g_c = &.{ 1, 1 },
        .aerobic_respiration_g_c = &.{ 4, 4 },
        .growth_respiration_g_c = &.{ 3, 3 },
        .denitrification_respiration_g_c = &.{ 2, 2 },
        .nitrogen_fixation_respiration_g_c = &.{ 1, 1 },
        .aerobic_growth_efficiency = &.{ 0.5, 0.5 },
        .denitrification_growth_efficiency = &.{ 0.5, 0.5 },
        .dissolved_organic_carbon_fraction = &.{ 0.75, 0.75 },
        .acetate_fraction = &.{ 0.25, 0.25 },
        .dissolved_organic_nitrogen_g_n = &.{ 10, 10 },
        .dissolved_organic_phosphorus_g_p = &.{ 10, 10 },
        .active_biomass_fraction = &.{ 0.5, 0.5 },
        .dissolved_nitrogen_to_carbon_ratio_g_n_per_g_c = &.{ 0.1, 0.1 },
        .dissolved_phosphorus_to_carbon_ratio_g_p_per_g_c = &.{ 0.02, 0.02 },
        .nitrogen_limitation_fraction = &.{ 0.5, 0.5 },
        .phosphorus_limitation_fraction = &.{ 0.5, 0.5 },
    };
}

test "source complex split matches NITRO" {
    try std.testing.expect(sourceHeterotrophic(3));
    try std.testing.expect(!sourceHeterotrophic(4));
}

test "heterotroph partitions carbon and takes up dissolved nutrients" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(6, state.aerobic_carbon_uptake_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(4, state.denitrification_carbon_uptake_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(8.5, state.primary_carbon_substrate_uptake_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(1.5, state.acetate_uptake_g_c[0], 1e-12);
    try std.testing.expect(state.dissolved_organic_nitrogen_uptake_g_n[0] > 0);
}

test "autotroph assigns all carbon uptake to inorganic carbon ledger" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(10, state.primary_carbon_substrate_uptake_g_c[1], 1e-12);
    try std.testing.expectEqual(0, state.acetate_uptake_g_c[1]);
    try std.testing.expectEqual(0, state.dissolved_organic_nitrogen_uptake_g_n[1]);
}

test "inconsistent fixation ledger fails atomically" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.total_carbon_uptake_g_c[0] = 7;
    var inputs = fixture();
    inputs.nitrogen_fixation_respiration_g_c = &.{ 4, 1 };
    try std.testing.expectError(
        error.InconsistentSubstrateUptakeRespirationLedger,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.total_carbon_uptake_g_c[0]);
}

test "NITRO 2555-2614 derived overflow preserves uptake state" {
    var state = try State.init(std.testing.allocator, 2);
    defer state.deinit();
    state.total_carbon_uptake_g_c[0] = 7;
    var inputs = fixture();
    inputs.denitrification_respiration_g_c = &.{ std.math.floatMax(f64), 2 };
    try std.testing.expectError(
        error.NonFiniteSubstrateUptakeResult,
        calculate(&state, inputs),
    );
    try std.testing.expectEqual(7, state.total_carbon_uptake_g_c[0]);
}
