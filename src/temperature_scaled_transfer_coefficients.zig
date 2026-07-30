const std = @import("std");

pub const SaltMode = enum { static_equilibrium, dynamic_transport };

pub const GaseousCoefficients = struct {
    water_vapor_m2_h: f64,
    carbon_dioxide_m2_h: f64,
    methane_m2_h: f64,
    oxygen_m2_h: f64,
    nitrogen_m2_h: f64,
    nitrous_oxide_m2_h: f64,
    ammonia_m2_h: f64,
    hydrogen_m2_h: f64,
};

pub const AqueousCoefficients = struct {
    carbon_dioxide_m2_h: f64,
    methane_m2_h: f64,
    oxygen_m2_h: f64,
    nitrogen_m2_h: f64,
    ammonia_m2_h: f64,
    nitrous_oxide_m2_h: f64,
    nitrate_m2_h: f64,
    phosphate_m2_h: f64,
    dissolved_organic_carbon_m2_h: f64,
    dissolved_organic_nitrogen_m2_h: f64,
    dissolved_organic_phosphorus_m2_h: f64,
    acetate_m2_h: f64,
    hydrogen_m2_h: f64,
};

pub const SaltCoefficients = struct {
    aluminum_m2_h: f64,
    iron_m2_h: f64,
    hydrogen_m2_h: f64,
    calcium_m2_h: f64,
    magnesium_m2_h: f64,
    sodium_m2_h: f64,
    potassium_m2_h: f64,
    hydroxide_m2_h: f64,
    carbonate_m2_h: f64,
    bicarbonate_m2_h: f64,
    sulfate_m2_h: f64,
    chloride_m2_h: f64,
    hydrogen_sulfate_m2_h: f64,
};

pub const SaltTransferAccumulators = struct {
    aluminum: f64,
    iron: f64,
    calcium: f64,
    magnesium: f64,
    sodium: f64,
    potassium: f64,
    aluminum_fertilizer: f64,
    iron_fertilizer: f64,
    calcium_fertilizer: f64,
    magnesium_fertilizer: f64,
    sodium_fertilizer: f64,
    potassium_fertilizer: f64,
};

pub const DynamicSaltResult = struct {
    coefficients: SaltCoefficients,
    reset_transfers: SaltTransferAccumulators,
};

pub const Result = struct {
    gaseous_temperature_factor: f64,
    aqueous_temperature_factor: f64,
    nitrogen_diffusivity_factor: f64, // TFND
    gaseous: GaseousCoefficients,
    aqueous: AqueousCoefficients,
    dynamic_salt: ?DynamicSaltResult,
};

pub const ScalingError = error{
    NonFiniteInput,
    InvalidSoilTemperature,
    NegativeCoefficient,
    NonFiniteResult,
};

/// Translates HOUR1 lines 3947-3999. Base and scaled diffusivities are m2 h-1.
pub fn scale(
    salt_mode: SaltMode,
    soil_temperature_k: f64,
    gaseous_base: GaseousCoefficients,
    aqueous_base: AqueousCoefficients,
    salt_base: SaltCoefficients,
) ScalingError!Result {
    if (!std.math.isFinite(soil_temperature_k)) return error.NonFiniteInput;
    if (soil_temperature_k <= 0.0) return error.InvalidSoilTemperature;
    try validateCoefficients(gaseous_base);
    try validateCoefficients(aqueous_base);
    if (salt_mode == .dynamic_transport) try validateCoefficients(salt_base);

    const gaseous_factor =
        std.math.pow(f64, soil_temperature_k / 298.15, 1.75);
    const aqueous_factor =
        std.math.pow(f64, soil_temperature_k / 298.15, 6.0);
    const gaseous = multiplyCoefficients(gaseous_base, gaseous_factor);
    const aqueous = multiplyCoefficients(aqueous_base, aqueous_factor);
    const dynamic_salt = if (salt_mode == .dynamic_transport)
        DynamicSaltResult{
            .coefficients = multiplyCoefficients(salt_base, aqueous_factor),
            .reset_transfers = std.mem.zeroes(SaltTransferAccumulators),
        }
    else
        null;
    if (!std.math.isFinite(gaseous_factor) or !std.math.isFinite(aqueous_factor)) {
        return error.NonFiniteResult;
    }
    try validateCoefficients(gaseous);
    try validateCoefficients(aqueous);
    if (dynamic_salt) |salt| try validateCoefficients(salt.coefficients);
    return .{
        .gaseous_temperature_factor = gaseous_factor,
        .aqueous_temperature_factor = aqueous_factor,
        .nitrogen_diffusivity_factor = aqueous_factor,
        .gaseous = gaseous,
        .aqueous = aqueous,
        .dynamic_salt = dynamic_salt,
    };
}

fn validateCoefficients(coefficients: anytype) ScalingError!void {
    inline for (std.meta.fields(@TypeOf(coefficients))) |field| {
        const value = @field(coefficients, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        if (value < 0.0) return error.NegativeCoefficient;
    }
}

fn multiplyCoefficients(coefficients: anytype, factor: f64) @TypeOf(coefficients) {
    var scaled: @TypeOf(coefficients) = undefined;
    inline for (std.meta.fields(@TypeOf(coefficients))) |field| {
        @field(scaled, field.name) = @field(coefficients, field.name) * factor;
    }
    return scaled;
}

test "reference temperature preserves all base coefficients" {
    const gaseous: GaseousCoefficients = @import("std").mem.zeroInit(
        GaseousCoefficients,
        .{ .water_vapor_m2_h = 0.0896, .carbon_dioxide_m2_h = 0.0468 },
    );
    const aqueous: AqueousCoefficients = @import("std").mem.zeroInit(
        AqueousCoefficients,
        .{ .carbon_dioxide_m2_h = 4.25e-6, .acetate_m2_h = 3.64e-6 },
    );
    const salt: SaltCoefficients = @import("std").mem.zeroInit(
        SaltCoefficients,
        .{ .aluminum_m2_h = 5.0e-6, .chloride_m2_h = 5.0e-6 },
    );
    const result = try scale(.dynamic_transport, 298.15, gaseous, aqueous, salt);
    try std.testing.expectEqual(@as(f64, 1.0), result.gaseous_temperature_factor);
    try std.testing.expectEqual(gaseous, result.gaseous);
    try std.testing.expectEqual(aqueous, result.aqueous);
    try std.testing.expectEqual(salt, result.dynamic_salt.?.coefficients);
}

test "static salt mode does not produce dynamic salt updates" {
    const gaseous = std.mem.zeroes(GaseousCoefficients);
    const aqueous = std.mem.zeroes(AqueousCoefficients);
    const result = try scale(
        .static_equilibrium,
        280.0,
        gaseous,
        aqueous,
        undefined,
    );
    try std.testing.expect(result.dynamic_salt == null);
}
