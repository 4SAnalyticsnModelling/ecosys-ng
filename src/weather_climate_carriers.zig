const std = @import("std");

pub const Carriers = struct {
    direct_shortwave_megajoules_per_m2: f64,
    diffuse_shortwave_megajoules_per_m2: f64,
    direct_par_megajoules_per_m2: f64,
    diffuse_par_megajoules_per_m2: f64,
    wind_speed_m_per_h: f64,
    vapor_pressure_kpa: f64,
    saturated_vapor_pressure_kpa: f64,
    rainfall_m: f64,
    snowfall_water_equivalent_m: f64,
    atmospheric_co2_umol_per_mol: f64,
    precipitation_ammonium_g_per_m3: f64,
    precipitation_nitrate_g_per_m3: f64,
};

pub const Modifiers = struct {
    radiation_fraction: f64,
    wind_speed_fraction: f64,
    humidity_fraction: f64,
    precipitation_fraction: f64,
    atmospheric_co2_fraction: f64,
    precipitation_ammonium_fraction: f64,
    precipitation_nitrate_fraction: f64,
};

pub const Baselines = struct {
    atmospheric_co2_umol_per_mol: f64,
    precipitation_ammonium_g_per_m3: f64,
    precipitation_nitrate_g_per_m3: f64,
};

/// Atomic WTHR climate scaling for atmospheric carriers other than irrigation
/// (which is routed by the separate irrigation climate transaction).
pub fn apply(
    carriers: *Carriers,
    modifiers: Modifiers,
    baselines: Baselines,
) !void {
    inline for (@typeInfo(Carriers).@"struct".fields) |field| {
        const value = @field(carriers.*, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteWeatherClimateCarrier;
        if (value < 0) return error.InvalidWeatherClimateCarrier;
    }
    inline for (@typeInfo(Modifiers).@"struct".fields) |field| {
        const value = @field(modifiers, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteWeatherClimateModifier;
        if (value < 0) return error.InvalidWeatherClimateModifier;
    }
    inline for (@typeInfo(Baselines).@"struct".fields) |field| {
        const value = @field(baselines, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteWeatherClimateBaseline;
        if (value < 0) return error.InvalidWeatherClimateBaseline;
    }

    const next: Carriers = .{
        .direct_shortwave_megajoules_per_m2 = carriers.direct_shortwave_megajoules_per_m2 *
            modifiers.radiation_fraction,
        .diffuse_shortwave_megajoules_per_m2 = carriers.diffuse_shortwave_megajoules_per_m2 *
            modifiers.radiation_fraction,
        .direct_par_megajoules_per_m2 = carriers.direct_par_megajoules_per_m2 *
            modifiers.radiation_fraction,
        .diffuse_par_megajoules_per_m2 = carriers.diffuse_par_megajoules_per_m2 *
            modifiers.radiation_fraction,
        .wind_speed_m_per_h = carriers.wind_speed_m_per_h *
            modifiers.wind_speed_fraction,
        .vapor_pressure_kpa = @min(
            carriers.saturated_vapor_pressure_kpa,
            carriers.vapor_pressure_kpa * modifiers.humidity_fraction,
        ),
        .saturated_vapor_pressure_kpa = carriers.saturated_vapor_pressure_kpa,
        .rainfall_m = carriers.rainfall_m * modifiers.precipitation_fraction,
        .snowfall_water_equivalent_m = carriers.snowfall_water_equivalent_m *
            modifiers.precipitation_fraction,
        // Source resets these three from their unmodified hourly baselines,
        // rather than compounding the previously modified carrier.
        .atmospheric_co2_umol_per_mol = baselines.atmospheric_co2_umol_per_mol *
            modifiers.atmospheric_co2_fraction,
        .precipitation_ammonium_g_per_m3 = baselines.precipitation_ammonium_g_per_m3 *
            modifiers.precipitation_ammonium_fraction,
        .precipitation_nitrate_g_per_m3 = baselines.precipitation_nitrate_g_per_m3 *
            modifiers.precipitation_nitrate_fraction,
    };
    inline for (@typeInfo(Carriers).@"struct".fields) |field|
        if (!std.math.isFinite(@field(next, field.name)))
            return error.WeatherClimateCarrierOverflow;
    carriers.* = next;
}

fn exampleCarriers() Carriers {
    return .{
        .direct_shortwave_megajoules_per_m2 = 2,
        .diffuse_shortwave_megajoules_per_m2 = 1,
        .direct_par_megajoules_per_m2 = 0.8,
        .diffuse_par_megajoules_per_m2 = 0.4,
        .wind_speed_m_per_h = 100,
        .vapor_pressure_kpa = 1.5,
        .saturated_vapor_pressure_kpa = 2,
        .rainfall_m = 0.001,
        .snowfall_water_equivalent_m = 0.002,
        .atmospheric_co2_umol_per_mol = 999,
        .precipitation_ammonium_g_per_m3 = 999,
        .precipitation_nitrate_g_per_m3 = 999,
    };
}

fn exampleModifiers() Modifiers {
    return .{
        .radiation_fraction = 1.5,
        .wind_speed_fraction = 0.5,
        .humidity_fraction = 2,
        .precipitation_fraction = 1.2,
        .atmospheric_co2_fraction = 1.1,
        .precipitation_ammonium_fraction = 2,
        .precipitation_nitrate_fraction = 3,
    };
}

test "WTHR scales four radiation beams and clamps humidity" {
    var carriers = exampleCarriers();
    try apply(&carriers, exampleModifiers(), .{
        .atmospheric_co2_umol_per_mol = 400,
        .precipitation_ammonium_g_per_m3 = 0.1,
        .precipitation_nitrate_g_per_m3 = 0.2,
    });
    try std.testing.expectEqual(@as(f64, 3), carriers.direct_shortwave_megajoules_per_m2);
    try std.testing.expectEqual(@as(f64, 1.5), carriers.diffuse_shortwave_megajoules_per_m2);
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), carriers.direct_par_megajoules_per_m2, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), carriers.diffuse_par_megajoules_per_m2, 1e-15);
    try std.testing.expectEqual(@as(f64, 50), carriers.wind_speed_m_per_h);
    try std.testing.expectEqual(@as(f64, 2), carriers.vapor_pressure_kpa);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0012), carriers.rainfall_m, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 440), carriers.atmospheric_co2_umol_per_mol, 1e-12);
    try std.testing.expectEqual(@as(f64, 0.2), carriers.precipitation_ammonium_g_per_m3);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), carriers.precipitation_nitrate_g_per_m3, 1e-15);
}

test "baseline atmospheric carriers do not compound" {
    var carriers = exampleCarriers();
    const baseline: Baselines = .{
        .atmospheric_co2_umol_per_mol = 400,
        .precipitation_ammonium_g_per_m3 = 0.1,
        .precipitation_nitrate_g_per_m3 = 0.2,
    };
    try apply(&carriers, exampleModifiers(), baseline);
    try apply(&carriers, exampleModifiers(), baseline);
    try std.testing.expectApproxEqAbs(@as(f64, 440), carriers.atmospheric_co2_umol_per_mol, 1e-12);
    try std.testing.expectEqual(@as(f64, 0.2), carriers.precipitation_ammonium_g_per_m3);
}

test "invalid late modifier rolls back every weather carrier" {
    var carriers = exampleCarriers();
    const before = carriers;
    var modifiers = exampleModifiers();
    modifiers.precipitation_nitrate_fraction = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteWeatherClimateModifier,
        apply(&carriers, modifiers, .{
            .atmospheric_co2_umol_per_mol = 400,
            .precipitation_ammonium_g_per_m3 = 0.1,
            .precipitation_nitrate_g_per_m3 = 0.2,
        }),
    );
    try std.testing.expectEqualDeep(before, carriers);
}

test "overflowing beam scaling rolls back state" {
    var carriers = exampleCarriers();
    carriers.diffuse_par_megajoules_per_m2 = std.math.floatMax(f64);
    const before = carriers;
    var modifiers = exampleModifiers();
    modifiers.radiation_fraction = 2;
    try std.testing.expectError(
        error.WeatherClimateCarrierOverflow,
        apply(&carriers, modifiers, .{
            .atmospheric_co2_umol_per_mol = 400,
            .precipitation_ammonium_g_per_m3 = 0.1,
            .precipitation_nitrate_g_per_m3 = 0.2,
        }),
    );
    try std.testing.expectEqualDeep(before, carriers);
}
