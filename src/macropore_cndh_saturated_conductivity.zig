const std = @import("std");

pub const Inputs = struct {
    macropore_volume_m3: f64,
    total_layer_volume_m3: f64,
    reference_water_viscosity_megagrams_m_s: f64,
    soil_temperature_c: f64,
};

pub const Result = struct {
    macropore_radius_m: f64,
    macropore_count: u64,
    macropore_spacing_m: f64,
    temperature_adjusted_viscosity_megagrams_m_s: f64,
    saturated_conductivity_m2_h_mpa: f64,
};

const source_pi = 3.1416;
const macropore_radius_m = 0.5e-3;

/// HOUR1 lines 2318--2328. This legacy CNDH calculation is retained as the
/// saturated conductivity supplied to ecosys-ng macropore Richards flow.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    const macropore_count_real = inputs.macropore_volume_m3 /
        (source_pi * macropore_radius_m * macropore_radius_m *
            inputs.total_layer_volume_m3);
    if (!std.math.isFinite(macropore_count_real) or macropore_count_real < 0 or
        macropore_count_real > @as(f64, @floatFromInt(std.math.maxInt(u64))))
        return error.InvalidMacroporeCount;
    const macropore_count: u64 = @intFromFloat(macropore_count_real);
    const macropore_spacing_m = if (macropore_count > 0)
        1.0 / @sqrt(source_pi * @as(f64, @floatFromInt(macropore_count)))
    else
        1.0;
    const temperature_adjusted_viscosity_megagrams_m_s =
        inputs.reference_water_viscosity_megagrams_m_s *
        @exp(0.533 - 0.0267 * inputs.soil_temperature_c);
    if (!std.math.isFinite(temperature_adjusted_viscosity_megagrams_m_s) or
        temperature_adjusted_viscosity_megagrams_m_s <= 0)
        return error.InvalidTemperatureAdjustedWaterViscosity;
    const saturated_conductivity_m2_h_mpa =
        3.6e3 * source_pi *
        @as(f64, @floatFromInt(macropore_count)) *
        std.math.pow(f64, macropore_radius_m, 4.0) /
        (8.0 * temperature_adjusted_viscosity_megagrams_m_s);
    if (!std.math.isFinite(saturated_conductivity_m2_h_mpa))
        return error.NonFiniteMacroporeSaturatedConductivity;
    return .{
        .macropore_radius_m = macropore_radius_m,
        .macropore_count = macropore_count,
        .macropore_spacing_m = macropore_spacing_m,
        .temperature_adjusted_viscosity_megagrams_m_s = temperature_adjusted_viscosity_megagrams_m_s,
        .saturated_conductivity_m2_h_mpa = saturated_conductivity_m2_h_mpa,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteMacroporeCndhInput;
    if (inputs.macropore_volume_m3 < 0 or inputs.total_layer_volume_m3 <= 0 or
        inputs.reference_water_viscosity_megagrams_m_s <= 0)
        return error.InvalidMacroporeCndhInput;
}

test "legacy CNDH preserves count spacing viscosity and conductivity order" {
    const result = try compute(.{
        .macropore_volume_m3 = 0.001,
        .total_layer_volume_m3 = 1,
        .reference_water_viscosity_megagrams_m_s = 1,
        .soil_temperature_c = 20,
    });
    try std.testing.expectEqual(@as(u64, 1273), result.macropore_count);
    try std.testing.expectApproxEqAbs(
        1.0 / @sqrt(source_pi * 1273.0),
        result.macropore_spacing_m,
        1e-15,
    );
    const expected_viscosity = @exp(0.533 - 0.0267 * 20.0);
    try std.testing.expectApproxEqAbs(
        expected_viscosity,
        result.temperature_adjusted_viscosity_megagrams_m_s,
        1e-15,
    );
    const expected_conductivity =
        3.6e3 * source_pi * 1273.0 *
        std.math.pow(f64, macropore_radius_m, 4.0) /
        (8.0 * expected_viscosity);
    try std.testing.expectApproxEqAbs(
        expected_conductivity,
        result.saturated_conductivity_m2_h_mpa,
        1e-24,
    );
}

test "zero macropore volume uses unit spacing and zero conductivity" {
    const result = try compute(.{
        .macropore_volume_m3 = 0,
        .total_layer_volume_m3 = 1,
        .reference_water_viscosity_megagrams_m_s = 1,
        .soil_temperature_c = 10,
    });
    try std.testing.expectEqual(@as(u64, 0), result.macropore_count);
    try std.testing.expectEqual(@as(f64, 1), result.macropore_spacing_m);
    try std.testing.expectEqual(
        @as(f64, 0),
        result.saturated_conductivity_m2_h_mpa,
    );
}
