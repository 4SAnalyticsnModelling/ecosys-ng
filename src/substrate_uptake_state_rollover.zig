const std = @import("std");

/// Total substrate uptake for one layer and timestep, in g of the named
/// element or compound per timestep.
pub const ScalarUptake = struct {
    oxygen: f64, // ROX*
    ammonium_non_band: f64, // RNH4*
    nitrate_non_band: f64, // RNO3*
    nitrite_non_band: f64, // RNO2*
    nitrous_oxide_non_band: f64, // RN2O*
    dihydrogen_phosphate_non_band: f64, // RP14*
    hydrogen_phosphate_non_band: f64, // RPO4*
    ammonium_band: f64, // RNHB*
    ammonia_band: f64, // RN3B*
    nitrite_band: f64, // RN2B*
    dihydrogen_phosphate_band: f64, // RP1B*
    hydrogen_phosphate_band: f64, // RPOB*
};

pub const OrganicUptake = struct {
    dissolved_organic_carbon_g_c_timestep: []f64, // ROQC*
    acetate_g_c_timestep: []f64, // ROQA*
};

pub const UptakeState = struct {
    previous: ScalarUptake,
    current: ScalarUptake,
    previous_organic: OrganicUptake,
    current_organic: OrganicUptake,
};

pub const RolloverError = error{
    OrganicPoolLengthMismatch,
    NonFiniteCurrentUptake,
};

/// Translates HOUR1 lines 3903-3932 in source assignment order.
pub fn rollover(
    state: *UptakeState,
    organic_pool_count: usize,
) RolloverError!void {
    if (state.previous_organic.dissolved_organic_carbon_g_c_timestep.len != organic_pool_count or
        state.previous_organic.acetate_g_c_timestep.len != organic_pool_count or
        state.current_organic.dissolved_organic_carbon_g_c_timestep.len != organic_pool_count or
        state.current_organic.acetate_g_c_timestep.len != organic_pool_count)
    {
        return error.OrganicPoolLengthMismatch;
    }
    inline for (std.meta.fields(ScalarUptake)) |field| {
        if (!std.math.isFinite(@field(state.current, field.name))) {
            return error.NonFiniteCurrentUptake;
        }
    }
    for (0..organic_pool_count) |pool| {
        if (!std.math.isFinite(
            state.current_organic.dissolved_organic_carbon_g_c_timestep[pool],
        ) or !std.math.isFinite(state.current_organic.acetate_g_c_timestep[pool])) {
            return error.NonFiniteCurrentUptake;
        }
    }

    inline for (std.meta.fields(ScalarUptake)) |field| {
        @field(state.previous, field.name) = @field(state.current, field.name);
    }
    inline for (std.meta.fields(ScalarUptake)) |field| {
        @field(state.current, field.name) = 0.0;
    }
    for (0..organic_pool_count) |pool| {
        state.previous_organic.dissolved_organic_carbon_g_c_timestep[pool] =
            state.current_organic.dissolved_organic_carbon_g_c_timestep[pool];
        state.previous_organic.acetate_g_c_timestep[pool] =
            state.current_organic.acetate_g_c_timestep[pool];
        state.current_organic.dissolved_organic_carbon_g_c_timestep[pool] = 0.0;
        state.current_organic.acetate_g_c_timestep[pool] = 0.0;
    }
}

test "scalar and runtime organic uptake roll current into previous" {
    var previous_carbon = [_]f64{ 0.0, 0.0, 0.0 };
    var previous_acetate = [_]f64{ 0.0, 0.0, 0.0 };
    var current_carbon = [_]f64{ 1.0, 2.0, 3.0 };
    var current_acetate = [_]f64{ 4.0, 5.0, 6.0 };
    var state = UptakeState{
        .previous = std.mem.zeroes(ScalarUptake),
        .current = undefined,
        .previous_organic = .{
            .dissolved_organic_carbon_g_c_timestep = &previous_carbon,
            .acetate_g_c_timestep = &previous_acetate,
        },
        .current_organic = .{
            .dissolved_organic_carbon_g_c_timestep = &current_carbon,
            .acetate_g_c_timestep = &current_acetate,
        },
    };
    inline for (std.meta.fields(ScalarUptake), 0..) |field, index| {
        @field(state.current, field.name) = @floatFromInt(index + 1);
    }

    try rollover(&state, 3);

    inline for (std.meta.fields(ScalarUptake), 0..) |field, index| {
        try std.testing.expectEqual(
            @as(f64, @floatFromInt(index + 1)),
            @field(state.previous, field.name),
        );
        try std.testing.expectEqual(@as(f64, 0.0), @field(state.current, field.name));
    }
    try std.testing.expectEqual(@as(f64, 3.0), previous_carbon[2]);
    try std.testing.expectEqual(@as(f64, 6.0), previous_acetate[2]);
    try std.testing.expectEqual(@as(f64, 0.0), current_carbon[2]);
    try std.testing.expectEqual(@as(f64, 0.0), current_acetate[2]);
}

test "invalid current uptake fails before mutation" {
    var empty = [_]f64{};
    var state: UptakeState = undefined;
    state.previous.oxygen = 7.0;
    state.current = std.mem.zeroes(ScalarUptake);
    state.current.oxygen = std.math.nan(f64);
    state.previous_organic = .{
        .dissolved_organic_carbon_g_c_timestep = &empty,
        .acetate_g_c_timestep = &empty,
    };
    state.current_organic = state.previous_organic;
    try std.testing.expectError(error.NonFiniteCurrentUptake, rollover(&state, 0));
    try std.testing.expectEqual(@as(f64, 7.0), state.previous.oxygen);
}
