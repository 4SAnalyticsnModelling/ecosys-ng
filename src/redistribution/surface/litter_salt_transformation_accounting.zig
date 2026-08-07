const std = @import("std");

pub const Transformations = struct {
    co2_c_g: f64,
    water_nonband_mol: f64,
    water_band_mol: f64,
    ammonium_n_g: f64,
    ammonia_n_g: f64,
    nitrate_n_g: f64,
    hpo4_p_g: f64,
    h2po4_p_g: f64,
    silicic_acid_mol: f64,
};

pub const SenescenceInputs = struct {
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    sulfate_mol: f64,
    chloride_mol: f64,
};

pub const MolarMasses = struct {
    carbon_g_mol: f64,
    nitrogen_g_mol: f64,
    phosphorus_g_mol: f64,
};

pub const Accounting = struct {
    ion_output_mol: f64,
    ion_input_mol: f64,
};

pub const Result = struct {
    accounting: Accounting,
    transformation_output_mol: f64,
    senescence_input_mol: f64,
};

fn validate(value: anytype) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name)))
            return error.InvalidLitterSaltAccountingInput;
}

/// Exact standalone translation of redist.f lines 5842--5852.
/// The executable source multiplies both TRH2O and TBH2O by 2.0.
pub fn accumulate(
    initial: Accounting,
    transformations: Transformations,
    senescence: SenescenceInputs,
    molar_masses: MolarMasses,
) !Result {
    try validate(initial);
    try validate(transformations);
    try validate(senescence);
    try validate(molar_masses);
    if (molar_masses.carbon_g_mol <= 0 or molar_masses.nitrogen_g_mol <= 0 or
        molar_masses.phosphorus_g_mol <= 0)
        return error.InvalidLitterSaltAccountingInput;

    const output = transformations.co2_c_g / molar_masses.carbon_g_mol +
        2.0 * transformations.water_nonband_mol +
        2.0 * transformations.water_band_mol +
        (2.0 * transformations.ammonium_n_g + transformations.ammonia_n_g) /
            molar_masses.nitrogen_g_mol +
        transformations.nitrate_n_g / molar_masses.nitrogen_g_mol +
        (2.0 * transformations.hpo4_p_g + 3.0 * transformations.h2po4_p_g) /
            molar_masses.phosphorus_g_mol -
        transformations.silicic_acid_mol;
    const input = senescence.aluminum_mol + senescence.iron_mol +
        senescence.calcium_mol + senescence.magnesium_mol +
        senescence.sodium_mol + senescence.potassium_mol +
        senescence.sulfate_mol + senescence.chloride_mol;
    const accounting = Accounting{
        .ion_output_mol = initial.ion_output_mol + output,
        .ion_input_mol = initial.ion_input_mol + input,
    };
    if (!std.math.isFinite(output) or !std.math.isFinite(input) or
        !std.math.isFinite(accounting.ion_output_mol) or
        !std.math.isFinite(accounting.ion_input_mol))
        return error.NonFiniteLitterSaltAccounting;
    return .{ .accounting = accounting, .transformation_output_mol = output, .senescence_input_mol = input };
}

test "REDIST SSB preserves executable coefficients and source signs" {
    const result = try accumulate(
        .{ .ion_output_mol = 10, .ion_input_mol = 20 },
        .{
            .co2_c_g = 12,
            .water_nonband_mol = 1,
            .water_band_mol = 2,
            .ammonium_n_g = 14,
            .ammonia_n_g = 14,
            .nitrate_n_g = 14,
            .hpo4_p_g = 31,
            .h2po4_p_g = 31,
            .silicic_acid_mol = 1,
        },
        .{ .aluminum_mol = 1, .iron_mol = 2, .calcium_mol = 3, .magnesium_mol = 4, .sodium_mol = 5, .potassium_mol = 6, .sulfate_mol = 7, .chloride_mol = 8 },
        .{ .carbon_g_mol = 12, .nitrogen_g_mol = 14, .phosphorus_g_mol = 31 },
    );
    try std.testing.expectEqual(@as(f64, 15), result.transformation_output_mol);
    try std.testing.expectEqual(@as(f64, 36), result.senescence_input_mol);
    try std.testing.expectEqual(@as(f64, 25), result.accounting.ion_output_mol);
    try std.testing.expectEqual(@as(f64, 56), result.accounting.ion_input_mol);
}

test "REDIST SSB runtime molar masses and failures" {
    const zero_t = std.mem.zeroes(Transformations);
    const zero_s = std.mem.zeroes(SenescenceInputs);
    try std.testing.expectError(error.InvalidLitterSaltAccountingInput, accumulate(
        .{ .ion_output_mol = 0, .ion_input_mol = 0 },
        zero_t,
        zero_s,
        .{ .carbon_g_mol = 0, .nitrogen_g_mol = 14, .phosphorus_g_mol = 31 },
    ));
    const initial = Accounting{ .ion_output_mol = std.math.floatMax(f64), .ion_input_mol = 0 };
    var t = zero_t;
    t.water_nonband_mol = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteLitterSaltAccounting, accumulate(initial, t, zero_s, .{ .carbon_g_mol = 12, .nitrogen_g_mol = 14, .phosphorus_g_mol = 31 }));
}
