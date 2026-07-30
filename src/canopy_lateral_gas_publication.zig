const std = @import("std");

pub const Inputs = struct {
    lateral_carbon_dioxide_g_c_per_step: []const f64,
    lateral_methane_g_c_per_step: []const f64,
    lateral_oxygen_g_o_per_step: []const f64,
    biological_timestep_h: f64,
};

pub const State = struct {
    cumulative_carbon_dioxide_exchange_g_c: []f64,
    cumulative_methane_exchange_g_c: []f64,
    cumulative_oxygen_exchange_g_o: []f64,
};

/// UPTAKE.F 328--332. Publishes signed lateral canopy-gas exchange into the
/// existing cell totals in source order: CO2, CH4, then O2.
pub fn publish(state: State, inputs: Inputs) !void {
    const cell_count = state.cumulative_carbon_dioxide_exchange_g_c.len;
    if (cell_count == 0 or
        state.cumulative_methane_exchange_g_c.len != cell_count or
        state.cumulative_oxygen_exchange_g_o.len != cell_count or
        inputs.lateral_carbon_dioxide_g_c_per_step.len != cell_count or
        inputs.lateral_methane_g_c_per_step.len != cell_count or
        inputs.lateral_oxygen_g_o_per_step.len != cell_count)
        return error.InvalidCanopyLateralGasPublicationDimensions;
    if (!std.math.isFinite(inputs.biological_timestep_h) or
        inputs.biological_timestep_h < 0)
        return error.InvalidCanopyLateralGasPublicationInput;

    for (0..cell_count) |cell| {
        const carbon_dioxide = state.cumulative_carbon_dioxide_exchange_g_c[cell] +
            inputs.lateral_carbon_dioxide_g_c_per_step[cell] * inputs.biological_timestep_h;
        const methane = state.cumulative_methane_exchange_g_c[cell] +
            inputs.lateral_methane_g_c_per_step[cell] * inputs.biological_timestep_h;
        const oxygen = state.cumulative_oxygen_exchange_g_o[cell] +
            inputs.lateral_oxygen_g_o_per_step[cell] * inputs.biological_timestep_h;
        inline for (.{ carbon_dioxide, methane, oxygen }) |value|
            if (!std.math.isFinite(value))
                return error.NonFiniteCanopyLateralGasPublication;
    }

    for (0..cell_count) |cell| {
        state.cumulative_carbon_dioxide_exchange_g_c[cell] +=
            inputs.lateral_carbon_dioxide_g_c_per_step[cell] * inputs.biological_timestep_h;
        state.cumulative_methane_exchange_g_c[cell] +=
            inputs.lateral_methane_g_c_per_step[cell] * inputs.biological_timestep_h;
        state.cumulative_oxygen_exchange_g_o[cell] +=
            inputs.lateral_oxygen_g_o_per_step[cell] * inputs.biological_timestep_h;
    }
}

test "UPTAKE lateral gas publication preserves signs timestep and source order" {
    var carbon_dioxide = [_]f64{ 10, -5, 1 };
    var methane = [_]f64{ 20, 3, -2 };
    var oxygen = [_]f64{ -4, 7, 9 };
    try publish(.{
        .cumulative_carbon_dioxide_exchange_g_c = &carbon_dioxide,
        .cumulative_methane_exchange_g_c = &methane,
        .cumulative_oxygen_exchange_g_o = &oxygen,
    }, .{
        .lateral_carbon_dioxide_g_c_per_step = &.{ 2, -4, 0.5 },
        .lateral_methane_g_c_per_step = &.{ -1, 6, -0.5 },
        .lateral_oxygen_g_o_per_step = &.{ 8, -2, 1 },
        .biological_timestep_h = 0.25,
    });
    try std.testing.expectEqual([3]f64{ 10.5, -6, 1.125 }, carbon_dioxide);
    try std.testing.expectEqual([3]f64{ 19.75, 4.5, -2.125 }, methane);
    try std.testing.expectEqual([3]f64{ -2, 6.5, 9.25 }, oxygen);
}

test "zero lateral flux preserves cumulative gas ledgers" {
    var carbon_dioxide = [_]f64{3};
    var methane = [_]f64{-2};
    var oxygen = [_]f64{7};
    try publish(.{
        .cumulative_carbon_dioxide_exchange_g_c = &carbon_dioxide,
        .cumulative_methane_exchange_g_c = &methane,
        .cumulative_oxygen_exchange_g_o = &oxygen,
    }, .{
        .lateral_carbon_dioxide_g_c_per_step = &.{0},
        .lateral_methane_g_c_per_step = &.{0},
        .lateral_oxygen_g_o_per_step = &.{0},
        .biological_timestep_h = 1,
    });
    try std.testing.expectEqual(@as(f64, 3), carbon_dioxide[0]);
    try std.testing.expectEqual(@as(f64, -2), methane[0]);
    try std.testing.expectEqual(@as(f64, 7), oxygen[0]);
}

test "late oxygen overflow leaves all gas ledgers unchanged" {
    var carbon_dioxide = [_]f64{3};
    var methane = [_]f64{4};
    var oxygen = [_]f64{std.math.floatMax(f64)};
    try std.testing.expectError(
        error.NonFiniteCanopyLateralGasPublication,
        publish(.{
            .cumulative_carbon_dioxide_exchange_g_c = &carbon_dioxide,
            .cumulative_methane_exchange_g_c = &methane,
            .cumulative_oxygen_exchange_g_o = &oxygen,
        }, .{
            .lateral_carbon_dioxide_g_c_per_step = &.{1},
            .lateral_methane_g_c_per_step = &.{2},
            .lateral_oxygen_g_o_per_step = &.{std.math.floatMax(f64)},
            .biological_timestep_h = 2,
        }),
    );
    try std.testing.expectEqual(@as(f64, 3), carbon_dioxide[0]);
    try std.testing.expectEqual(@as(f64, 4), methane[0]);
    try std.testing.expectEqual(std.math.floatMax(f64), oxygen[0]);
}

test "lateral gas publication rejects mismatched runtime dimensions" {
    var carbon_dioxide = [_]f64{ 0, 0 };
    var methane = [_]f64{ 0, 0 };
    var oxygen = [_]f64{ 0, 0 };
    try std.testing.expectError(
        error.InvalidCanopyLateralGasPublicationDimensions,
        publish(.{
            .cumulative_carbon_dioxide_exchange_g_c = &carbon_dioxide,
            .cumulative_methane_exchange_g_c = &methane,
            .cumulative_oxygen_exchange_g_o = &oxygen,
        }, .{
            .lateral_carbon_dioxide_g_c_per_step = &.{0},
            .lateral_methane_g_c_per_step = &.{ 0, 0 },
            .lateral_oxygen_g_o_per_step = &.{ 0, 0 },
            .biological_timestep_h = 1,
        }),
    );
}
