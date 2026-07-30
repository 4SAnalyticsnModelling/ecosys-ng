const std = @import("std");

pub const WaterInputs = struct {
    rainfall_m_per_h: f64,
    snowfall_water_equivalent_m_per_h: f64,
    surface_irrigation_m_per_h: f64,
    subsurface_irrigation_m_per_h: f64,
};

pub const Multipliers = struct {
    precipitation_fraction: f64,
    irrigation_fraction: f64,
};

/// Exact WTHR water-input climate scaling from wthr.f:477-480.
///
/// Rain and snow use TDPRC. Surface and subsurface irrigation use TDIRI.
/// The update is transactional: every source, multiplier, and result validates
/// before any caller-owned carrier is replaced.
pub fn apply(inputs: *WaterInputs, multipliers: Multipliers) !void {
    inline for (@typeInfo(WaterInputs).@"struct".fields) |field| {
        const value = @field(inputs.*, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteWaterClimateInput;
        if (value < 0)
            return error.InvalidWaterClimateInput;
    }
    inline for (@typeInfo(Multipliers).@"struct".fields) |field| {
        const value = @field(multipliers, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteWaterClimateMultiplier;
        if (value < 0)
            return error.InvalidWaterClimateMultiplier;
    }

    const next: WaterInputs = .{
        .rainfall_m_per_h = inputs.rainfall_m_per_h * multipliers.precipitation_fraction,
        .snowfall_water_equivalent_m_per_h = inputs.snowfall_water_equivalent_m_per_h *
            multipliers.precipitation_fraction,
        .surface_irrigation_m_per_h = inputs.surface_irrigation_m_per_h *
            multipliers.irrigation_fraction,
        .subsurface_irrigation_m_per_h = inputs.subsurface_irrigation_m_per_h *
            multipliers.irrigation_fraction,
    };
    inline for (@typeInfo(WaterInputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.WaterClimateForcingOverflow;
    inputs.* = next;
}

fn exampleInputs() WaterInputs {
    return .{
        .rainfall_m_per_h = 0.001,
        .snowfall_water_equivalent_m_per_h = 0.002,
        .surface_irrigation_m_per_h = 0.003,
        .subsurface_irrigation_m_per_h = 0.004,
    };
}

test "precipitation and irrigation use distinct exact WTHR multipliers" {
    var inputs = exampleInputs();
    try apply(&inputs, .{
        .precipitation_fraction = 2,
        .irrigation_fraction = 3,
    });
    try std.testing.expectEqual(@as(f64, 0.002), inputs.rainfall_m_per_h);
    try std.testing.expectEqual(
        @as(f64, 0.004),
        inputs.snowfall_water_equivalent_m_per_h,
    );
    try std.testing.expectEqual(
        @as(f64, 0.009000000000000001),
        inputs.surface_irrigation_m_per_h,
    );
    try std.testing.expectEqual(
        @as(f64, 0.012),
        inputs.subsurface_irrigation_m_per_h,
    );
}

test "unity multipliers preserve all unmodified forcing" {
    var inputs = exampleInputs();
    const before = inputs;
    try apply(&inputs, .{
        .precipitation_fraction = 1,
        .irrigation_fraction = 1,
    });
    try std.testing.expectEqualDeep(before, inputs);
}

test "zero precipitation multiplier does not suppress irrigation" {
    var inputs = exampleInputs();
    try apply(&inputs, .{
        .precipitation_fraction = 0,
        .irrigation_fraction = 1,
    });
    try std.testing.expectEqual(@as(f64, 0), inputs.rainfall_m_per_h);
    try std.testing.expectEqual(
        @as(f64, 0),
        inputs.snowfall_water_equivalent_m_per_h,
    );
    try std.testing.expectEqual(
        @as(f64, 0.003),
        inputs.surface_irrigation_m_per_h,
    );
}

test "invalid late multiplier leaves every water carrier unchanged" {
    var inputs = exampleInputs();
    const before = inputs;
    try std.testing.expectError(
        error.NonFiniteWaterClimateMultiplier,
        apply(&inputs, .{
            .precipitation_fraction = 2,
            .irrigation_fraction = std.math.nan(f64),
        }),
    );
    try std.testing.expectEqualDeep(before, inputs);
}

test "overflow leaves every water carrier unchanged" {
    var inputs = exampleInputs();
    inputs.subsurface_irrigation_m_per_h = std.math.floatMax(f64);
    const before = inputs;
    try std.testing.expectError(
        error.WaterClimateForcingOverflow,
        apply(&inputs, .{
            .precipitation_fraction = 1,
            .irrigation_fraction = 2,
        }),
    );
    try std.testing.expectEqualDeep(before, inputs);
}
