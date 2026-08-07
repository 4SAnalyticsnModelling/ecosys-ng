const std = @import("std");

pub const PhosphorusMasses = struct {
    phosphate_0_mol: f64,
    phosphate_3_mol: f64,
    iron_phosphate_1_mol: f64,
    iron_phosphate_2_mol: f64,
    calcium_phosphate_0_mol: f64,
    calcium_phosphate_1_mol: f64,
    calcium_phosphate_2_mol: f64,
    magnesium_phosphate_1_mol: f64,
    dihydrogen_phosphate_g_p: f64,
    hydrogen_phosphate_g_p: f64,
};

pub const ZoneMasses = struct {
    ammonium_g_n: f64,
    ammonia_g_n: f64,
    nitrate_g_n: f64,
    nitrite_g_n: f64,
    phosphorus: PhosphorusMasses,
};

pub const ZoneFractions = struct {
    ammonium: f64,
    nitrate: f64,
    phosphate: f64,
};

pub const ZoneConcentrations = struct {
    ammonium_g_n_m3: f64,
    ammonia_g_n_m3: f64,
    nitrate_g_n_m3: f64,
    nitrite_g_n_m3: f64,
    dihydrogen_phosphate_g_p_m3: f64,
    hydrogen_phosphate_g_p_m3: f64,
    total_mineral_phosphorus_g_p_m3: f64,
};

pub const Result = struct {
    non_band: ZoneConcentrations,
    band: ZoneConcentrations,
};

pub const ConcentrationError = error{
    NonFiniteInput,
    NegativeWaterVolume,
    InvalidZoneFraction,
    NonFiniteResult,
};

/// Translates `hour1.f` lines 3804--3892 in non-band then band source order.
pub fn calculate(
    water_volume_m3: f64,
    water_volume_threshold_m3: f64,
    positive_fraction_threshold: f64,
    non_band_masses: ZoneMasses,
    non_band_fractions: ZoneFractions,
    band_masses: ZoneMasses,
    band_fractions: ZoneFractions,
) ConcentrationError!Result {
    const scalars = [_]f64{
        water_volume_m3,
        water_volume_threshold_m3,
        positive_fraction_threshold,
    };
    for (scalars) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (value < 0.0) return error.NegativeWaterVolume;
    }
    inline for ([_]ZoneMasses{ non_band_masses, band_masses }) |masses| {
        inline for (std.meta.fields(ZoneMasses)) |field| {
            if (field.type == f64) {
                if (!std.math.isFinite(@field(masses, field.name))) {
                    return error.NonFiniteInput;
                }
            } else {
                inline for (std.meta.fields(PhosphorusMasses)) |phosphorus_field| {
                    if (!std.math.isFinite(
                        @field(masses.phosphorus, phosphorus_field.name),
                    )) return error.NonFiniteInput;
                }
            }
        }
    }
    inline for ([_]ZoneFractions{ non_band_fractions, band_fractions }) |fractions| {
        inline for (std.meta.fields(ZoneFractions)) |field| {
            const fraction = @field(fractions, field.name);
            if (!std.math.isFinite(fraction)) return error.NonFiniteInput;
            if (fraction < 0.0 or fraction > 1.0) return error.InvalidZoneFraction;
        }
    }

    if (water_volume_m3 <= water_volume_threshold_m3) {
        return .{
            .non_band = std.mem.zeroes(ZoneConcentrations),
            .band = std.mem.zeroes(ZoneConcentrations),
        };
    }
    const result = Result{
        .non_band = zoneConcentrations(
            water_volume_m3,
            positive_fraction_threshold,
            non_band_masses,
            non_band_fractions,
        ),
        .band = zoneConcentrations(
            water_volume_m3,
            positive_fraction_threshold,
            band_masses,
            band_fractions,
        ),
    };
    inline for (std.meta.fields(Result)) |zone_field| {
        inline for (std.meta.fields(ZoneConcentrations)) |field| {
            if (!std.math.isFinite(@field(@field(result, zone_field.name), field.name))) {
                return error.NonFiniteResult;
            }
        }
    }
    return result;
}

fn zoneConcentrations(
    water_volume_m3: f64,
    positive_fraction_threshold: f64,
    masses: ZoneMasses,
    fractions: ZoneFractions,
) ZoneConcentrations {
    var result = std.mem.zeroes(ZoneConcentrations);
    if (fractions.ammonium > positive_fraction_threshold) {
        result.ammonium_g_n_m3 =
            @max(0.0, masses.ammonium_g_n / (water_volume_m3 * fractions.ammonium));
        result.ammonia_g_n_m3 =
            @max(0.0, masses.ammonia_g_n / (water_volume_m3 * fractions.ammonium));
    }
    if (fractions.nitrate > positive_fraction_threshold) {
        result.nitrate_g_n_m3 =
            @max(0.0, masses.nitrate_g_n / (water_volume_m3 * fractions.nitrate));
        result.nitrite_g_n_m3 =
            @max(0.0, masses.nitrite_g_n / (water_volume_m3 * fractions.nitrate));
    }
    if (fractions.phosphate > positive_fraction_threshold) {
        result.dihydrogen_phosphate_g_p_m3 = @max(
            0.0,
            masses.phosphorus.dihydrogen_phosphate_g_p /
                (water_volume_m3 * fractions.phosphate),
        );
        result.hydrogen_phosphate_g_p_m3 = @max(
            0.0,
            masses.phosphorus.hydrogen_phosphate_g_p /
                (water_volume_m3 * fractions.phosphate),
        );
        result.total_mineral_phosphorus_g_p_m3 = @max(
            0.0,
            ((masses.phosphorus.phosphate_0_mol +
                masses.phosphorus.phosphate_3_mol +
                masses.phosphorus.iron_phosphate_1_mol +
                masses.phosphorus.iron_phosphate_2_mol +
                masses.phosphorus.calcium_phosphate_0_mol +
                masses.phosphorus.calcium_phosphate_1_mol +
                masses.phosphorus.calcium_phosphate_2_mol +
                masses.phosphorus.magnesium_phosphate_1_mol) * 31.0 +
                masses.phosphorus.dihydrogen_phosphate_g_p +
                masses.phosphorus.hydrogen_phosphate_g_p) /
                (water_volume_m3 * fractions.phosphate),
        );
    }
    return result;
}

test "band and non-band concentrations preserve zone-volume divisions" {
    const phosphorus = PhosphorusMasses{
        .phosphate_0_mol = 1.0,
        .phosphate_3_mol = 0.0,
        .iron_phosphate_1_mol = 0.0,
        .iron_phosphate_2_mol = 0.0,
        .calcium_phosphate_0_mol = 0.0,
        .calcium_phosphate_1_mol = 0.0,
        .calcium_phosphate_2_mol = 0.0,
        .magnesium_phosphate_1_mol = 0.0,
        .dihydrogen_phosphate_g_p = 2.0,
        .hydrogen_phosphate_g_p = 3.0,
    };
    const masses = ZoneMasses{
        .ammonium_g_n = 4.0,
        .ammonia_g_n = 2.0,
        .nitrate_g_n = 6.0,
        .nitrite_g_n = 1.0,
        .phosphorus = phosphorus,
    };
    const fractions = ZoneFractions{ .ammonium = 0.5, .nitrate = 0.25, .phosphate = 0.5 };
    const result = try calculate(
        10.0,
        0.0,
        0.0,
        masses,
        fractions,
        masses,
        fractions,
    );
    try std.testing.expectEqual(@as(f64, 0.8), result.non_band.ammonium_g_n_m3);
    try std.testing.expectEqual(
        @as(f64, 7.2),
        result.non_band.total_mineral_phosphorus_g_p_m3,
    );
    try std.testing.expectEqual(result.non_band, result.band);
}

test "water at threshold zeros every zone concentration" {
    const masses = std.mem.zeroes(ZoneMasses);
    const fractions = std.mem.zeroes(ZoneFractions);
    const result = try calculate(0.0, 0.0, 0.0, masses, fractions, masses, fractions);
    try std.testing.expectEqual(std.mem.zeroes(Result), result);
}
