const std = @import("std");

pub const Trigger = union(enum) {
    soil_water_fraction_above_wilting: f64,
    canopy_water_potential_mpa: f64,
};

pub const Inputs = struct {
    layer_top_depth_m: []const f64,
    layer_bottom_depth_m: []const f64,
    layer_total_volume_m3: []const f64,
    porosity_fraction: []const f64,
    field_capacity_fraction: []const f64,
    wilting_point_fraction: []const f64,
    liquid_water_m3: []const f64,
    ice_water_equivalent_m3: []const f64,
    evaluated_depth_m: f64,
    refill_fraction_of_field_capacity: f64,
    minimum_canopy_water_potential_mpa: f64,
    trigger: Trigger,
    first_source_hour: u8,
    last_source_hour: u8,
    application_depth_m: f64,
};

pub const Result = struct {
    target_water_m3: f64,
    wilting_water_m3: f64,
    current_water_m3: f64,
    required_water_m3: f64,
    hourly_water_m3: f64,
    application_depth_m: f64,
    triggered: bool,
};

/// Exact DAY automatic-irrigation criterion using authoritative VOLX-like
/// layer volumes rather than reconstructing volume from nominal geometry.
pub fn calculate(inputs: Inputs) !Result {
    const count = inputs.layer_total_volume_m3.len;
    if (count == 0 or inputs.layer_top_depth_m.len != count or
        inputs.layer_bottom_depth_m.len != count or
        inputs.porosity_fraction.len != count or
        inputs.field_capacity_fraction.len != count or
        inputs.wilting_point_fraction.len != count or
        inputs.liquid_water_m3.len != count or
        inputs.ice_water_equivalent_m3.len != count)
        return error.AutomaticIrrigationDimensionMismatch;
    inline for (.{
        inputs.evaluated_depth_m,
        inputs.refill_fraction_of_field_capacity,
        inputs.minimum_canopy_water_potential_mpa,
        inputs.application_depth_m,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteAutomaticIrrigationInput;
    if (inputs.evaluated_depth_m < 0 or
        inputs.refill_fraction_of_field_capacity < 0 or
        inputs.refill_fraction_of_field_capacity > 1 or
        inputs.application_depth_m < 0 or
        inputs.first_source_hour < 1 or inputs.first_source_hour > 24 or
        inputs.last_source_hour < inputs.first_source_hour or
        inputs.last_source_hour > 24)
        return error.InvalidAutomaticIrrigationInput;
    const trigger_value = switch (inputs.trigger) {
        .soil_water_fraction_above_wilting => |value| value,
        .canopy_water_potential_mpa => |value| value,
    };
    if (!std.math.isFinite(trigger_value))
        return error.NonFiniteAutomaticIrrigationInput;

    var target_m3: f64 = 0;
    var wilting_m3: f64 = 0;
    var current_m3: f64 = 0;
    for (0..count) |layer| {
        inline for (.{
            inputs.layer_top_depth_m[layer],
            inputs.layer_bottom_depth_m[layer],
            inputs.layer_total_volume_m3[layer],
            inputs.porosity_fraction[layer],
            inputs.field_capacity_fraction[layer],
            inputs.wilting_point_fraction[layer],
            inputs.liquid_water_m3[layer],
            inputs.ice_water_equivalent_m3[layer],
        }) |value| if (!std.math.isFinite(value))
            return error.NonFiniteAutomaticIrrigationLayer;
        const top_m = inputs.layer_top_depth_m[layer];
        const bottom_m = inputs.layer_bottom_depth_m[layer];
        const volume_m3 = inputs.layer_total_volume_m3[layer];
        const porosity = inputs.porosity_fraction[layer];
        const field_capacity = inputs.field_capacity_fraction[layer];
        const wilting = inputs.wilting_point_fraction[layer];
        if (bottom_m <= top_m or volume_m3 <= 0 or porosity < 0 or
            porosity > 1 or wilting < 0 or field_capacity < wilting or
            field_capacity > porosity or inputs.liquid_water_m3[layer] < 0 or
            inputs.ice_water_equivalent_m3[layer] < 0)
            return error.InvalidAutomaticIrrigationLayer;
        if (top_m >= inputs.evaluated_depth_m) continue;
        const fraction = @min(
            @as(f64, 1),
            (inputs.evaluated_depth_m - top_m) / (bottom_m - top_m),
        );
        const target_fraction = @min(
            porosity,
            wilting + inputs.refill_fraction_of_field_capacity *
                (field_capacity - wilting),
        );
        target_m3 += fraction * target_fraction * volume_m3;
        wilting_m3 += fraction * wilting * volume_m3;
        current_m3 += fraction *
            (inputs.liquid_water_m3[layer] +
                inputs.ice_water_equivalent_m3[layer]);
        inline for (.{ target_m3, wilting_m3, current_m3 }) |value|
            if (!std.math.isFinite(value))
                return error.AutomaticIrrigationCalculationOverflow;
    }
    const triggered = switch (inputs.trigger) {
        .soil_water_fraction_above_wilting => |fraction| current_m3 <
            wilting_m3 + fraction * (target_m3 - wilting_m3),
        .canopy_water_potential_mpa => |threshold_mpa| inputs.minimum_canopy_water_potential_mpa < threshold_mpa,
    };
    const required_m3 = if (triggered) @max(0, target_m3 - current_m3) else 0;
    const hour_count: u8 =
        inputs.last_source_hour - inputs.first_source_hour + 1;
    const hourly_m3 =
        required_m3 / @as(f64, @floatFromInt(hour_count));
    if (!std.math.isFinite(required_m3) or !std.math.isFinite(hourly_m3))
        return error.AutomaticIrrigationCalculationOverflow;
    return .{
        .target_water_m3 = target_m3,
        .wilting_water_m3 = wilting_m3,
        .current_water_m3 = current_m3,
        .required_water_m3 = required_m3,
        .hourly_water_m3 = hourly_m3,
        .application_depth_m = inputs.application_depth_m,
        .triggered = triggered,
    };
}

test "partial evaluated layer uses authoritative volume and liquid plus ice" {
    const result = try calculate(.{
        .layer_top_depth_m = &.{ 0, 0.1 },
        .layer_bottom_depth_m = &.{ 0.1, 0.3 },
        .layer_total_volume_m3 = &.{ 10, 40 },
        .porosity_fraction = &.{ 0.5, 0.5 },
        .field_capacity_fraction = &.{ 0.4, 0.4 },
        .wilting_point_fraction = &.{ 0.1, 0.1 },
        .liquid_water_m3 = &.{ 1, 2 },
        .ice_water_equivalent_m3 = &.{ 0.5, 1 },
        .evaluated_depth_m = 0.2,
        .refill_fraction_of_field_capacity = 1,
        .minimum_canopy_water_potential_mpa = -0.5,
        .trigger = .{ .soil_water_fraction_above_wilting = 1 },
        .first_source_hour = 3,
        .last_source_hour = 6,
        .application_depth_m = 0.15,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 12), result.target_water_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), result.wilting_water_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), result.current_water_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 9), result.required_water_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2.25), result.hourly_water_m3, 1e-12);
}

test "canopy criterion uses source PFT-one potential supplied by caller" {
    var inputs: Inputs = .{
        .layer_top_depth_m = &.{0},
        .layer_bottom_depth_m = &.{0.1},
        .layer_total_volume_m3 = &.{10},
        .porosity_fraction = &.{0.5},
        .field_capacity_fraction = &.{0.4},
        .wilting_point_fraction = &.{0.1},
        .liquid_water_m3 = &.{2},
        .ice_water_equivalent_m3 = &.{0},
        .evaluated_depth_m = 0.1,
        .refill_fraction_of_field_capacity = 1,
        .minimum_canopy_water_potential_mpa = -1.5,
        .trigger = .{ .canopy_water_potential_mpa = -1 },
        .first_source_hour = 1,
        .last_source_hour = 1,
        .application_depth_m = 0,
    };
    var result = try calculate(inputs);
    try std.testing.expect(result.triggered);
    try std.testing.expectEqual(@as(f64, 2), result.required_water_m3);
    inputs.minimum_canopy_water_potential_mpa = -0.5;
    result = try calculate(inputs);
    try std.testing.expect(!result.triggered);
    try std.testing.expectEqual(@as(f64, 0), result.required_water_m3);
}

test "strict soil threshold equality does not trigger" {
    const result = try calculate(.{
        .layer_top_depth_m = &.{0},
        .layer_bottom_depth_m = &.{0.1},
        .layer_total_volume_m3 = &.{10},
        .porosity_fraction = &.{0.5},
        .field_capacity_fraction = &.{0.4},
        .wilting_point_fraction = &.{0.1},
        .liquid_water_m3 = &.{2.5},
        .ice_water_equivalent_m3 = &.{0},
        .evaluated_depth_m = 0.1,
        .refill_fraction_of_field_capacity = 1,
        .minimum_canopy_water_potential_mpa = 0,
        .trigger = .{ .soil_water_fraction_above_wilting = 0.5 },
        .first_source_hour = 1,
        .last_source_hour = 2,
        .application_depth_m = 0,
    });
    try std.testing.expect(!result.triggered);
}

test "invalid late layer fails before returning a requirement" {
    try std.testing.expectError(
        error.NonFiniteAutomaticIrrigationLayer,
        calculate(.{
            .layer_top_depth_m = &.{ 0, 0.1 },
            .layer_bottom_depth_m = &.{ 0.1, 0.2 },
            .layer_total_volume_m3 = &.{ 10, 10 },
            .porosity_fraction = &.{ 0.5, 0.5 },
            .field_capacity_fraction = &.{ 0.4, 0.4 },
            .wilting_point_fraction = &.{ 0.1, 0.1 },
            .liquid_water_m3 = &.{ 1, std.math.nan(f64) },
            .ice_water_equivalent_m3 = &.{ 0, 0 },
            .evaluated_depth_m = 0.2,
            .refill_fraction_of_field_capacity = 1,
            .minimum_canopy_water_potential_mpa = 0,
            .trigger = .{ .soil_water_fraction_above_wilting = 1 },
            .first_source_hour = 1,
            .last_source_hour = 1,
            .application_depth_m = 0,
        }),
    );
}
