const std = @import("std");

pub const Inputs = struct {
    vertical_was_missing: bool,
    lateral_was_missing: bool,
    existing_vertical_m2_h_megapascal: f64,
    existing_lateral_m2_h_megapascal: f64,
    organic_carbon_g_per_megagram: f64,
    organic_soil_threshold_g_per_megagram: f64,
    porosity_m3_m3: f64,
    log_saturation_potential: f64,
    log_porosity: f64,
    log_field_capacity: f64,
    saturation_to_field_potential_interval: f64,
    bulk_density_megagrams_m3: f64,
    macropore_factor: f64,
};

pub const Result = struct {
    vertical_m2_h_megapascal: f64,
    lateral_m2_h_megapascal: f64,
};

/// `hour1.f` lines 2197--2221. Vertical and lateral missing-value gates remain
/// independent and repeat the source default calculation in their order.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    var vertical = inputs.existing_vertical_m2_h_megapascal;
    var lateral = inputs.existing_lateral_m2_h_megapascal;
    if (inputs.vertical_was_missing)
        vertical = try defaultConductivity(inputs);
    if (inputs.lateral_was_missing)
        lateral = try defaultConductivity(inputs);
    return .{ .vertical_m2_h_megapascal = vertical, .lateral_m2_h_megapascal = lateral };
}

fn defaultConductivity(inputs: Inputs) !f64 {
    if (inputs.organic_carbon_g_per_megagram <
        inputs.organic_soil_threshold_g_per_megagram)
    {
        const field_equivalent_water_content = @min(
            inputs.porosity_m3_m3,
            @exp(
                (inputs.log_saturation_potential - @log(0.033)) *
                    (inputs.log_porosity - inputs.log_field_capacity) /
                    inputs.saturation_to_field_potential_interval +
                    inputs.log_porosity,
            ),
        );
        if (field_equivalent_water_content <= 0)
            return error.InvalidDefaultHydraulicConductivityDenominator;
        return 1.54 * std.math.pow(
            f64,
            (inputs.porosity_m3_m3 - field_equivalent_water_content) /
                field_equivalent_water_content,
            2.0,
        );
    }
    var conductivity = 0.10 + 75.0 *
        std.math.pow(f64, 1.0e-15, inputs.bulk_density_megagrams_m3);
    conductivity = conductivity * inputs.macropore_factor;
    return conductivity;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteDefaultHydraulicConductivityInput;
    inline for (.{
        inputs.existing_vertical_m2_h_megapascal,
        inputs.existing_lateral_m2_h_megapascal,
        inputs.organic_carbon_g_per_megagram,
        inputs.organic_soil_threshold_g_per_megagram,
        inputs.porosity_m3_m3,
        inputs.bulk_density_megagrams_m3,
        inputs.macropore_factor,
    }) |value| if (value < 0)
        return error.InvalidDefaultHydraulicConductivityInput;
    if (inputs.saturation_to_field_potential_interval == 0)
        return error.InvalidDefaultHydraulicConductivityInput;
}

test "mineral defaults independently replace missing directions" {
    const result = try compute(.{
        .vertical_was_missing = true,
        .lateral_was_missing = false,
        .existing_vertical_m2_h_megapascal = 0,
        .existing_lateral_m2_h_megapascal = 9,
        .organic_carbon_g_per_megagram = 1000,
        .organic_soil_threshold_g_per_megagram = 100_000,
        .porosity_m3_m3 = 0.5,
        .log_saturation_potential = @log(0.001),
        .log_porosity = @log(0.5),
        .log_field_capacity = @log(0.3),
        .saturation_to_field_potential_interval = 3,
        .bulk_density_megagrams_m3 = 1.2,
        .macropore_factor = 2,
    });
    try std.testing.expect(result.vertical_m2_h_megapascal >= 0);
    try std.testing.expectEqual(@as(f64, 9), result.lateral_m2_h_megapascal);
}

test "organic default applies bulk density power then macropore factor" {
    const result = try compute(.{
        .vertical_was_missing = true,
        .lateral_was_missing = true,
        .existing_vertical_m2_h_megapascal = 0,
        .existing_lateral_m2_h_megapascal = 0,
        .organic_carbon_g_per_megagram = 200_000,
        .organic_soil_threshold_g_per_megagram = 100_000,
        .porosity_m3_m3 = 0.8,
        .log_saturation_potential = 0,
        .log_porosity = 0,
        .log_field_capacity = 0,
        .saturation_to_field_potential_interval = 1,
        .bulk_density_megagrams_m3 = 1,
        .macropore_factor = 2,
    });
    const expected = (0.10 + 75.0e-15) * 2;
    try std.testing.expectApproxEqAbs(
        expected,
        result.vertical_m2_h_megapascal,
        1e-15,
    );
    try std.testing.expectEqual(
        result.vertical_m2_h_megapascal,
        result.lateral_m2_h_megapascal,
    );
}
