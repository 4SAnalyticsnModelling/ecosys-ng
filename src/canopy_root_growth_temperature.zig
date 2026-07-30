const std = @import("std");

pub const Inputs = struct {
    phenology_type: u8,
    leafout_enabled: bool,
    accumulated_leafout_h: f64,
    required_leafout_h: f64,
    surface_soil_temperature_k: f64,
    canopy_surface_temperature_k: f64,
    thermal_adaptation_offset_k: f64,
    soil_temperature_k_by_layer: []const f64,
    first_soil_layer: usize,
    end_soil_layer: usize,
    previous_daily_minimum_canopy_water_potential_mpa: f64,
    canopy_total_water_potential_mpa: f64,
};

pub const Result = struct {
    canopy_growth_temperature_k: f64,
    canopy_growth_temperature_c: f64,
    canopy_growth_temperature_response: f64,
    daily_minimum_canopy_water_potential_mpa: f64,
};

/// UPTAKE.F 1546--1589. Selects growth temperature, evaluates canopy/root
/// Arrhenius responses, and updates daily minimum canopy water potential.
pub fn calculate(
    inputs: Inputs,
    root_growth_temperature_response_by_layer: []f64,
) !Result {
    try validate(inputs, root_growth_temperature_response_by_layer.len);
    const growth_temperature =
        if (inputs.phenology_type != 0 and
        inputs.leafout_enabled and
        inputs.accumulated_leafout_h <= inputs.required_leafout_h)
            inputs.surface_soil_temperature_k
        else
            inputs.canopy_surface_temperature_k;
    const growth_temperature_c = growth_temperature - 273.15;

    // Source TKCO uses TKC, not the selected TKG.
    const adapted_canopy_temperature =
        inputs.canopy_surface_temperature_k +
        inputs.thermal_adaptation_offset_k;
    const canopy_response = try arrheniusResponse(adapted_canopy_temperature);
    for (inputs.first_soil_layer..inputs.end_soil_layer) |layer| {
        root_growth_temperature_response_by_layer[layer] =
            try arrheniusResponse(
                inputs.soil_temperature_k_by_layer[layer] +
                    inputs.thermal_adaptation_offset_k,
            );
    }
    return .{
        .canopy_growth_temperature_k = growth_temperature,
        .canopy_growth_temperature_c = growth_temperature_c,
        .canopy_growth_temperature_response = canopy_response,
        .daily_minimum_canopy_water_potential_mpa = @min(
            inputs.previous_daily_minimum_canopy_water_potential_mpa,
            inputs.canopy_total_water_potential_mpa,
        ),
    };
}

fn arrheniusResponse(temperature_k: f64) !f64 {
    if (!std.math.isFinite(temperature_k) or temperature_k <= 0)
        return error.InvalidGrowthArrheniusTemperature;
    const gas_constant_temperature = 8.3143 * temperature_k;
    const entropy_temperature = 710.0 * temperature_k;
    const activity =
        1.0 +
        @exp((197_500.0 - entropy_temperature) / gas_constant_temperature) +
        @exp((entropy_temperature - 222_500.0) / gas_constant_temperature);
    const response =
        @exp(25.229 - 62_500.0 / gas_constant_temperature) /
        activity;
    if (!std.math.isFinite(response))
        return error.NonFiniteGrowthArrheniusResponse;
    return response;
}

fn validate(inputs: Inputs, output_count: usize) !void {
    inline for (.{
        inputs.accumulated_leafout_h,
        inputs.required_leafout_h,
        inputs.surface_soil_temperature_k,
        inputs.canopy_surface_temperature_k,
        inputs.thermal_adaptation_offset_k,
        inputs.previous_daily_minimum_canopy_water_potential_mpa,
        inputs.canopy_total_water_potential_mpa,
    }) |value|
        if (!std.math.isFinite(value))
            return error.InvalidGrowthTemperatureInput;
    if (inputs.accumulated_leafout_h < 0 or
        inputs.required_leafout_h < 0 or
        inputs.surface_soil_temperature_k <= 0 or
        inputs.canopy_surface_temperature_k <= 0 or
        inputs.first_soil_layer > inputs.end_soil_layer or
        inputs.end_soil_layer > inputs.soil_temperature_k_by_layer.len or
        output_count != inputs.soil_temperature_k_by_layer.len)
        return error.InvalidGrowthTemperatureInput;
    for (inputs.soil_temperature_k_by_layer) |temperature|
        if (!std.math.isFinite(temperature) or temperature <= 0)
            return error.InvalidGrowthTemperatureInput;
}

test "pre-leafout deciduous canopy selects soil but Arrhenius retains TKC" {
    var root_response = [_]f64{ 9, 9 };
    const result = try calculate(.{
        .phenology_type = 1,
        .leafout_enabled = true,
        .accumulated_leafout_h = 20,
        .required_leafout_h = 100,
        .surface_soil_temperature_k = 280,
        .canopy_surface_temperature_k = 300,
        .thermal_adaptation_offset_k = 2,
        .soil_temperature_k_by_layer = &.{ 280, 290 },
        .first_soil_layer = 0,
        .end_soil_layer = 2,
        .previous_daily_minimum_canopy_water_potential_mpa = -1,
        .canopy_total_water_potential_mpa = -2,
    }, &root_response);
    try std.testing.expectEqual(@as(f64, 280), result.canopy_growth_temperature_k);
    try std.testing.expectEqual(
        try arrheniusResponse(302),
        result.canopy_growth_temperature_response,
    );
    try std.testing.expectEqual(try arrheniusResponse(282), root_response[0]);
    try std.testing.expectEqual(try arrheniusResponse(292), root_response[1]);
    try std.testing.expectEqual(@as(f64, -2), result.daily_minimum_canopy_water_potential_mpa);
}

test "evergreen canopy selects canopy temperature and partial root range" {
    var root_response = [_]f64{ 9, 9 };
    const result = try calculate(.{
        .phenology_type = 0,
        .leafout_enabled = true,
        .accumulated_leafout_h = 0,
        .required_leafout_h = 100,
        .surface_soil_temperature_k = 280,
        .canopy_surface_temperature_k = 300,
        .thermal_adaptation_offset_k = 0,
        .soil_temperature_k_by_layer = &.{ 280, 290 },
        .first_soil_layer = 1,
        .end_soil_layer = 2,
        .previous_daily_minimum_canopy_water_potential_mpa = -3,
        .canopy_total_water_potential_mpa = -2,
    }, &root_response);
    try std.testing.expectEqual(@as(f64, 300), result.canopy_growth_temperature_k);
    try std.testing.expectEqual(@as(f64, 9), root_response[0]);
    try std.testing.expectEqual(@as(f64, -3), result.daily_minimum_canopy_water_potential_mpa);
}

test "nonpositive adapted temperature fails explicitly" {
    var root_response = [_]f64{0};
    try std.testing.expectError(
        error.InvalidGrowthArrheniusTemperature,
        calculate(.{
            .phenology_type = 0,
            .leafout_enabled = true,
            .accumulated_leafout_h = 0,
            .required_leafout_h = 0,
            .surface_soil_temperature_k = 280,
            .canopy_surface_temperature_k = 300,
            .thermal_adaptation_offset_k = -300,
            .soil_temperature_k_by_layer = &.{280},
            .first_soil_layer = 0,
            .end_soil_layer = 1,
            .previous_daily_minimum_canopy_water_potential_mpa = -1,
            .canopy_total_water_potential_mpa = -1,
        }, &root_response),
    );
}
