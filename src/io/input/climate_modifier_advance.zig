const std = @import("std");

pub const Mode = enum {
    step,
    incremental,
};

pub const Modifier = struct {
    maximum_temperature_change_c: f64,
    minimum_temperature_change_c: f64,
    radiation_fraction: f64,
    wind_speed_fraction: f64,
    humidity_fraction: f64,
    precipitation_fraction: f64,
    irrigation_fraction: f64,
    atmospheric_co2_fraction: f64,
    precipitation_ammonium_fraction: f64,
    precipitation_nitrate_fraction: f64,
};

/// Atomic DAY climate-state update over a runtime number of forcing periods.
/// Four periods reproduce active seasonal inputs; the same kernel supports
/// future compulsory monthly records without a compile-time extent.
pub fn advance(
    active: []Modifier,
    targets: []const Modifier,
    mode: Mode,
    days_in_year: u16,
) !void {
    if (active.len == 0 or active.len != targets.len)
        return error.ClimateModifierDimensionMismatch;
    if (days_in_year != 365 and days_in_year != 366)
        return error.InvalidClimateModifierYearLength;
    for (active) |modifier| try validateActive(modifier);
    for (targets) |target| try validateTarget(target);

    // Derive every candidate before any write. Recalculation in the commit
    // pass is allocation-free and cannot fail after full validation.
    for (active, targets) |current, target|
        try validateActive(derive(current, target, mode, days_in_year));
    for (active, targets) |*current, target|
        current.* = derive(current.*, target, mode, days_in_year);
}

fn derive(
    current: Modifier,
    target: Modifier,
    mode: Mode,
    days_in_year: u16,
) Modifier {
    if (mode == .step) return target;
    const days: f64 = @floatFromInt(days_in_year);
    return .{
        .maximum_temperature_change_c = current.maximum_temperature_change_c +
            target.maximum_temperature_change_c / days,
        .minimum_temperature_change_c = current.minimum_temperature_change_c +
            target.minimum_temperature_change_c / days,
        .radiation_fraction = current.radiation_fraction +
            (target.radiation_fraction - 1) / days,
        .wind_speed_fraction = current.wind_speed_fraction +
            (target.wind_speed_fraction - 1) / days,
        .humidity_fraction = current.humidity_fraction +
            (target.humidity_fraction - 1) / days,
        .precipitation_fraction = current.precipitation_fraction +
            (target.precipitation_fraction - 1) / days,
        .irrigation_fraction = current.irrigation_fraction +
            (target.irrigation_fraction - 1) / days,
        .atmospheric_co2_fraction = current.atmospheric_co2_fraction *
            @exp(@log(target.atmospheric_co2_fraction) / days),
        .precipitation_ammonium_fraction = current.precipitation_ammonium_fraction +
            (target.precipitation_ammonium_fraction - 1) / days,
        .precipitation_nitrate_fraction = current.precipitation_nitrate_fraction +
            (target.precipitation_nitrate_fraction - 1) / days,
    };
}

fn validateActive(modifier: Modifier) !void {
    inline for (@typeInfo(Modifier).@"struct".fields) |field| {
        const value = @field(modifier, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteClimateModifier;
        if (comptime !std.mem.endsWith(u8, field.name, "_change_c"))
            if (value < 0) return error.InvalidClimateModifier;
    }
    if (modifier.atmospheric_co2_fraction <= 0)
        return error.InvalidClimateModifier;
}

fn validateTarget(target: Modifier) !void {
    try validateActive(target);
}

fn neutral() Modifier {
    return .{
        .maximum_temperature_change_c = 0,
        .minimum_temperature_change_c = 0,
        .radiation_fraction = 1,
        .wind_speed_fraction = 1,
        .humidity_fraction = 1,
        .precipitation_fraction = 1,
        .irrigation_fraction = 1,
        .atmospheric_co2_fraction = 1,
        .precipitation_ammonium_fraction = 1,
        .precipitation_nitrate_fraction = 1,
    };
}

fn exampleTarget() Modifier {
    return .{
        .maximum_temperature_change_c = 4,
        .minimum_temperature_change_c = 2,
        .radiation_fraction = 1.2,
        .wind_speed_fraction = 0.8,
        .humidity_fraction = 1.1,
        .precipitation_fraction = 1.3,
        .irrigation_fraction = 1.4,
        .atmospheric_co2_fraction = 2,
        .precipitation_ammonium_fraction = 1.5,
        .precipitation_nitrate_fraction = 0.5,
    };
}

test "step mode copies every runtime forcing period exactly" {
    var active = [_]Modifier{ neutral(), neutral(), neutral(), neutral() };
    const targets =
        [_]Modifier{ exampleTarget(), neutral(), exampleTarget(), neutral() };
    try advance(&active, &targets, .step, 365);
    try std.testing.expectEqualDeep(targets, active);
}

test "incremental daily changes reach annual additive and geometric targets" {
    var active = [_]Modifier{neutral()};
    const targets = [_]Modifier{exampleTarget()};
    for (0..365) |_|
        try advance(&active, &targets, .incremental, 365);
    try std.testing.expectApproxEqAbs(
        targets[0].maximum_temperature_change_c,
        active[0].maximum_temperature_change_c,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        targets[0].irrigation_fraction,
        active[0].irrigation_fraction,
        1e-12,
    );
    try std.testing.expectApproxEqAbs(
        targets[0].atmospheric_co2_fraction,
        active[0].atmospheric_co2_fraction,
        1e-12,
    );
}

test "runtime period count is not fixed to legacy storage" {
    var active = [_]Modifier{neutral()} ** 6;
    const targets = [_]Modifier{exampleTarget()} ** 6;
    try advance(&active, &targets, .incremental, 366);
    try std.testing.expect(active[5].radiation_fraction > 1);
}

test "invalid late target rolls back every earlier period" {
    var active = [_]Modifier{ neutral(), neutral(), neutral(), neutral() };
    const before = active;
    var targets = [_]Modifier{
        exampleTarget(),
        exampleTarget(),
        exampleTarget(),
        exampleTarget(),
    };
    targets[3].irrigation_fraction = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteClimateModifier,
        advance(&active, &targets, .incremental, 365),
    );
    try std.testing.expectEqualDeep(before, active);
}
