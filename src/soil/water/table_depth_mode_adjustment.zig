const std = @import("std");

pub const Mode = enum(u8) {
    none = 0,
    absolute = 1,
    relative = 2,
    absolute_dynamic = 3,
    relative_dynamic = 4,
};

pub const Inputs = struct {
    mode: Mode,
    configured_water_table_depth_m: f64,
    configured_dynamic_depth_m: f64,
    cumulative_depth_above_uppermost_soil_layer_m: f64,
    initial_cumulative_depth_m: f64,
};

pub const Result = struct {
    adjusted_water_table_depth_m: f64,
    /// Null preserves the caller's DTBLY for source modes 0, 1, and 2.
    dynamic_water_table_depth_m: ?f64,
};

/// `hour1.f` lines 2348--2356. Preserves the absolute/relative adjustment
/// followed by the independent dynamic-mode assignment.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    const adjusted_water_table_depth_m = switch (inputs.mode) {
        .none, .absolute, .absolute_dynamic => inputs.configured_water_table_depth_m,
        .relative, .relative_dynamic => inputs.configured_water_table_depth_m +
            inputs.cumulative_depth_above_uppermost_soil_layer_m -
            inputs.initial_cumulative_depth_m,
    };
    if (!std.math.isFinite(adjusted_water_table_depth_m))
        return error.NonFiniteAdjustedWaterTableDepth;
    return .{
        .adjusted_water_table_depth_m = adjusted_water_table_depth_m,
        .dynamic_water_table_depth_m = switch (inputs.mode) {
            .absolute_dynamic, .relative_dynamic => inputs.configured_dynamic_depth_m,
            else => null,
        },
    };
}

fn validate(inputs: Inputs) !void {
    inline for (.{
        inputs.configured_water_table_depth_m,
        inputs.configured_dynamic_depth_m,
        inputs.cumulative_depth_above_uppermost_soil_layer_m,
        inputs.initial_cumulative_depth_m,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteWaterTableAdjustmentInput;
}

test "relative dynamic mode applies depth offset then dynamic target" {
    const result = try compute(.{
        .mode = .relative_dynamic,
        .configured_water_table_depth_m = 2,
        .configured_dynamic_depth_m = 3,
        .cumulative_depth_above_uppermost_soil_layer_m = 0.5,
        .initial_cumulative_depth_m = 0.2,
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, 2.3),
        result.adjusted_water_table_depth_m,
        1e-15,
    );
    try std.testing.expectEqual(
        @as(f64, 3),
        result.dynamic_water_table_depth_m.?,
    );
}

test "absolute nondynamic mode preserves dynamic caller state" {
    const result = try compute(.{
        .mode = .absolute,
        .configured_water_table_depth_m = 2,
        .configured_dynamic_depth_m = 9,
        .cumulative_depth_above_uppermost_soil_layer_m = 0.5,
        .initial_cumulative_depth_m = 0.2,
    });
    try std.testing.expectEqual(
        @as(f64, 2),
        result.adjusted_water_table_depth_m,
    );
    try std.testing.expect(result.dynamic_water_table_depth_m == null);
}
