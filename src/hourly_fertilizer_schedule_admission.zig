const std = @import("std");

pub const OrganicAmendment = struct {
    carbon_g_c_per_m2: f64,
    nitrogen_g_n_per_m2: f64,
    phosphorus_g_p_per_m2: f64,
};

pub const Application = struct {
    broadcast_ammonium_g_n_per_m2: f64,
    broadcast_ammonia_g_n_per_m2: f64,
    broadcast_urea_g_n_per_m2: f64,
    broadcast_nitrate_g_n_per_m2: f64,
    banded_ammonium_g_n_per_m2: f64,
    banded_ammonia_g_n_per_m2: f64,
    banded_urea_g_n_per_m2: f64,
    banded_nitrate_g_n_per_m2: f64,
    broadcast_monocalcium_phosphate_g_p_per_m2: f64,
    banded_monocalcium_phosphate_g_p_per_m2: f64,
    broadcast_hydroxyapatite_g_p_per_m2: f64,
    calcium_carbonate_g_ca_per_m2: f64,
    calcium_sulfate_g_ca_or_rock_g_si_per_m2: f64,
    plant_residue: OrganicAmendment,
    animal_manure: OrganicAmendment,
};

/// HOUR1 lines 270--271 exclude organic amendments from the outer
/// application-presence gate and preserve the thirteen-term sum order.
pub fn hasMineralApplication(application: Application) !bool {
    const total =
        application.broadcast_ammonium_g_n_per_m2 +
        application.broadcast_ammonia_g_n_per_m2 +
        application.broadcast_urea_g_n_per_m2 +
        application.broadcast_nitrate_g_n_per_m2 +
        application.banded_ammonium_g_n_per_m2 +
        application.banded_ammonia_g_n_per_m2 +
        application.banded_urea_g_n_per_m2 +
        application.banded_nitrate_g_n_per_m2 +
        application.broadcast_monocalcium_phosphate_g_p_per_m2 +
        application.banded_monocalcium_phosphate_g_p_per_m2 +
        application.broadcast_hydroxyapatite_g_p_per_m2 +
        application.calcium_carbonate_g_ca_per_m2 +
        application.calcium_sulfate_g_ca_or_rock_g_si_per_m2;
    if (!std.math.isFinite(total)) return error.HourlyFertilizerTotalOverflow;
    return total > 0;
}

/// HOUR1 lines 218--261. Returns no application outside the source's
/// `J == INT(ZNOON)` gate; the 19 quantities retain exact record order.
pub fn admit(
    source_hour_one_through_twenty_four: u8,
    solar_noon_hour: f64,
    daily_fertilizer_values: []const f64,
) !?Application {
    if (source_hour_one_through_twenty_four == 0 or
        source_hour_one_through_twenty_four > 24 or
        !std.math.isFinite(solar_noon_hour) or
        solar_noon_hour < 0 or solar_noon_hour >= 25)
        return error.InvalidHourlyFertilizerSchedule;
    if (daily_fertilizer_values.len != 19)
        return error.HourlyFertilizerRecordArityMismatch;
    for (daily_fertilizer_values) |value| {
        if (!std.math.isFinite(value))
            return error.NonFiniteHourlyFertilizerValue;
        if (value < 0) return error.InvalidHourlyFertilizerValue;
    }
    const source_solar_noon: u8 = @intFromFloat(@trunc(solar_noon_hour));
    if (source_hour_one_through_twenty_four != source_solar_noon) return null;

    return .{
        .broadcast_ammonium_g_n_per_m2 = daily_fertilizer_values[0],
        .broadcast_ammonia_g_n_per_m2 = daily_fertilizer_values[1],
        .broadcast_urea_g_n_per_m2 = daily_fertilizer_values[2],
        .broadcast_nitrate_g_n_per_m2 = daily_fertilizer_values[3],
        .banded_ammonium_g_n_per_m2 = daily_fertilizer_values[4],
        .banded_ammonia_g_n_per_m2 = daily_fertilizer_values[5],
        .banded_urea_g_n_per_m2 = daily_fertilizer_values[6],
        .banded_nitrate_g_n_per_m2 = daily_fertilizer_values[7],
        .broadcast_monocalcium_phosphate_g_p_per_m2 = daily_fertilizer_values[8],
        .banded_monocalcium_phosphate_g_p_per_m2 = daily_fertilizer_values[9],
        .broadcast_hydroxyapatite_g_p_per_m2 = daily_fertilizer_values[10],
        .calcium_carbonate_g_ca_per_m2 = daily_fertilizer_values[11],
        .calcium_sulfate_g_ca_or_rock_g_si_per_m2 = daily_fertilizer_values[12],
        .plant_residue = .{
            .carbon_g_c_per_m2 = daily_fertilizer_values[13],
            .nitrogen_g_n_per_m2 = daily_fertilizer_values[14],
            .phosphorus_g_p_per_m2 = daily_fertilizer_values[15],
        },
        .animal_manure = .{
            .carbon_g_c_per_m2 = daily_fertilizer_values[16],
            .nitrogen_g_n_per_m2 = daily_fertilizer_values[17],
            .phosphorus_g_p_per_m2 = daily_fertilizer_values[18],
        },
    };
}

test "solar noon admission preserves all nineteen source positions" {
    const values = [_]f64{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
    const application = (try admit(12, 12.75, &values)).?;
    try std.testing.expectEqual(@as(f64, 1), application.broadcast_ammonium_g_n_per_m2);
    try std.testing.expectEqual(@as(f64, 8), application.banded_nitrate_g_n_per_m2);
    try std.testing.expectEqual(@as(f64, 13), application.calcium_sulfate_g_ca_or_rock_g_si_per_m2);
    try std.testing.expectEqual(@as(f64, 16), application.plant_residue.phosphorus_g_p_per_m2);
    try std.testing.expectEqual(@as(f64, 19), application.animal_manure.phosphorus_g_p_per_m2);
}

test "source INT solar noon gate truncates toward zero" {
    const values = [_]f64{0} ** 19;
    try std.testing.expect((try admit(12, 12.99, &values)) != null);
    try std.testing.expect((try admit(13, 12.99, &values)) == null);
}

test "mineral presence gate excludes organic-only amendments" {
    const values = [_]f64{
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 2, 3, 4, 5, 6,
    };
    const organic_only = (try admit(12, 12, &values)).?;
    try std.testing.expect(!try hasMineralApplication(organic_only));
    var mineral_values = values;
    mineral_values[12] = 1;
    const mineral = (try admit(12, 12, &mineral_values)).?;
    try std.testing.expect(try hasMineralApplication(mineral));
}

test "strict fertilizer arity and values fail before admission" {
    try std.testing.expectError(
        error.HourlyFertilizerRecordArityMismatch,
        admit(12, 12, &.{ 1, 2 }),
    );
    var values = [_]f64{0} ** 19;
    values[18] = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteHourlyFertilizerValue,
        admit(1, 12, &values),
    );
}
