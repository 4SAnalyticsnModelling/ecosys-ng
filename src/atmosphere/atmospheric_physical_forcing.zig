const std = @import("std");

pub const Carriers = struct {
    direct_shortwave_megajoules_per_m2_h: f64,
    diffuse_shortwave_megajoules_per_m2_h: f64,
    direct_par_megajoules_per_m2_h: f64,
    diffuse_par_megajoules_per_m2_h: f64,
    wind_travel_m_per_h: f64,
    ambient_vapor_pressure_kpa: f64,
    saturated_vapor_pressure_kpa: f64,
};

pub const Multipliers = struct {
    radiation_fraction: f64,
    wind_fraction: f64,
    humidity_fraction: f64,
};

/// Exact WTHR atmospheric physical forcing from wthr.f:472-476.
pub fn apply(carriers: *Carriers, multipliers: Multipliers) !void {
    inline for (@typeInfo(Carriers).@"struct".fields) |field| {
        const value = @field(carriers.*, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteAtmosphericPhysicalCarrier;
        if (value < 0)
            return error.InvalidAtmosphericPhysicalCarrier;
    }
    inline for (@typeInfo(Multipliers).@"struct".fields) |field| {
        const value = @field(multipliers, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteAtmosphericPhysicalMultiplier;
        if (value < 0)
            return error.InvalidAtmosphericPhysicalMultiplier;
    }

    const next: Carriers = .{
        .direct_shortwave_megajoules_per_m2_h = carriers.direct_shortwave_megajoules_per_m2_h *
            multipliers.radiation_fraction,
        .diffuse_shortwave_megajoules_per_m2_h = carriers.diffuse_shortwave_megajoules_per_m2_h *
            multipliers.radiation_fraction,
        .direct_par_megajoules_per_m2_h = carriers.direct_par_megajoules_per_m2_h *
            multipliers.radiation_fraction,
        .diffuse_par_megajoules_per_m2_h = carriers.diffuse_par_megajoules_per_m2_h *
            multipliers.radiation_fraction,
        .wind_travel_m_per_h = carriers.wind_travel_m_per_h * multipliers.wind_fraction,
        .ambient_vapor_pressure_kpa = @min(
            carriers.saturated_vapor_pressure_kpa,
            carriers.ambient_vapor_pressure_kpa *
                multipliers.humidity_fraction,
        ),
        .saturated_vapor_pressure_kpa = carriers.saturated_vapor_pressure_kpa,
    };
    inline for (@typeInfo(Carriers).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.AtmosphericPhysicalForcingOverflow;
    carriers.* = next;
}

fn exampleCarriers() Carriers {
    return .{
        .direct_shortwave_megajoules_per_m2_h = 2,
        .diffuse_shortwave_megajoules_per_m2_h = 1,
        .direct_par_megajoules_per_m2_h = 0.8,
        .diffuse_par_megajoules_per_m2_h = 0.4,
        .wind_travel_m_per_h = 3600,
        .ambient_vapor_pressure_kpa = 1,
        .saturated_vapor_pressure_kpa = 2,
    };
}

test "radiation multiplier scales all four beams and no other carrier" {
    var carriers = exampleCarriers();
    try apply(&carriers, .{
        .radiation_fraction = 2,
        .wind_fraction = 1,
        .humidity_fraction = 1,
    });
    try std.testing.expectEqual(
        @as(f64, 4),
        carriers.direct_shortwave_megajoules_per_m2_h,
    );
    try std.testing.expectEqual(
        @as(f64, 2),
        carriers.diffuse_shortwave_megajoules_per_m2_h,
    );
    try std.testing.expectEqual(
        @as(f64, 1.6),
        carriers.direct_par_megajoules_per_m2_h,
    );
    try std.testing.expectEqual(
        @as(f64, 0.8),
        carriers.diffuse_par_megajoules_per_m2_h,
    );
    try std.testing.expectEqual(@as(f64, 3600), carriers.wind_travel_m_per_h);
}

test "humidity scaling is capped by current saturated vapor pressure" {
    var carriers = exampleCarriers();
    try apply(&carriers, .{
        .radiation_fraction = 1,
        .wind_fraction = 1,
        .humidity_fraction = 3,
    });
    try std.testing.expectEqual(
        @as(f64, 2),
        carriers.ambient_vapor_pressure_kpa,
    );
    try std.testing.expectEqual(
        @as(f64, 2),
        carriers.saturated_vapor_pressure_kpa,
    );
}

test "wind multiplier remains independent of radiation and humidity" {
    var carriers = exampleCarriers();
    try apply(&carriers, .{
        .radiation_fraction = 1,
        .wind_fraction = 0.5,
        .humidity_fraction = 1,
    });
    try std.testing.expectEqual(@as(f64, 1800), carriers.wind_travel_m_per_h);
    try std.testing.expectEqual(
        @as(f64, 2),
        carriers.direct_shortwave_megajoules_per_m2_h,
    );
    try std.testing.expectEqual(
        @as(f64, 1),
        carriers.ambient_vapor_pressure_kpa,
    );
}

test "invalid late multiplier leaves every physical carrier unchanged" {
    var carriers = exampleCarriers();
    const before = carriers;
    try std.testing.expectError(
        error.NonFiniteAtmosphericPhysicalMultiplier,
        apply(&carriers, .{
            .radiation_fraction = 2,
            .wind_fraction = 2,
            .humidity_fraction = std.math.nan(f64),
        }),
    );
    try std.testing.expectEqualDeep(before, carriers);
}

test "overflow leaves every physical carrier unchanged" {
    var carriers = exampleCarriers();
    carriers.diffuse_par_megajoules_per_m2_h = std.math.floatMax(f64);
    const before = carriers;
    try std.testing.expectError(
        error.AtmosphericPhysicalForcingOverflow,
        apply(&carriers, .{
            .radiation_fraction = 2,
            .wind_fraction = 1,
            .humidity_fraction = 1,
        }),
    );
    try std.testing.expectEqualDeep(before, carriers);
}
