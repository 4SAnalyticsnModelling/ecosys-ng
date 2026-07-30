const std = @import("std");

pub const Inputs = struct {
    micropore_capacity_m3: f64,
    micropore_water_m3: f64,
    micropore_ice_m3: f64,
    macropore_capacity_m3: f64,
    macropore_water_m3: f64,
    macropore_ice_m3: f64,
    bulk_soil_volume_m3: f64,
    negligible_bulk_volume_m3: f64,
    bulk_density_megagrams_m3: f64,
    positive_density_threshold_megagrams_m3: f64,
    log_field_capacity_potential: f64,
    log_field_capacity_to_wilting_point: f64,
    field_to_wilting_potential_interval: f64,
    log_field_capacity: f64,
    hygroscopic_potential_mpa: f64,
    minimum_potential_mpa: f64,
    zero_water_content_m3_m3: f64,
};

pub const Result = struct {
    air_filled_pore_volume_m3: f64,
    air_filled_porosity_m3_m3: f64,
    hygroscopic_water_content_m3_m3: f64,
    minimum_water_content_m3_m3: f64,
};

/// HOUR1 lines 2166--2182. Preserves micropore then macropore air volume,
/// bulk-volume gate, and bulk-density gate for dry-end water contents.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    const air_filled_pore_volume_m3 =
        @max(
            0.0,
            inputs.micropore_capacity_m3 - inputs.micropore_water_m3 -
                inputs.micropore_ice_m3,
        ) +
        @max(
            0.0,
            inputs.macropore_capacity_m3 - inputs.macropore_water_m3 -
                inputs.macropore_ice_m3,
        );
    const air_filled_porosity_m3_m3 =
        if (inputs.bulk_soil_volume_m3 > inputs.negligible_bulk_volume_m3)
            air_filled_pore_volume_m3 / inputs.bulk_soil_volume_m3
        else
            0.0;
    var hygroscopic_water_content_m3_m3: f64 = undefined;
    var minimum_water_content_m3_m3: f64 = undefined;
    if (inputs.bulk_density_megagrams_m3 >
        inputs.positive_density_threshold_megagrams_m3)
    {
        hygroscopic_water_content_m3_m3 = @exp(
            (inputs.log_field_capacity_potential -
                @log(-inputs.hygroscopic_potential_mpa)) *
                inputs.log_field_capacity_to_wilting_point /
                inputs.field_to_wilting_potential_interval +
                inputs.log_field_capacity,
        );
        minimum_water_content_m3_m3 = @exp(
            (inputs.log_field_capacity_potential -
                @log(-inputs.minimum_potential_mpa)) *
                inputs.log_field_capacity_to_wilting_point /
                inputs.field_to_wilting_potential_interval +
                inputs.log_field_capacity,
        );
    } else {
        hygroscopic_water_content_m3_m3 =
            inputs.zero_water_content_m3_m3;
        minimum_water_content_m3_m3 = inputs.zero_water_content_m3_m3;
    }
    return .{
        .air_filled_pore_volume_m3 = air_filled_pore_volume_m3,
        .air_filled_porosity_m3_m3 = air_filled_porosity_m3_m3,
        .hygroscopic_water_content_m3_m3 = hygroscopic_water_content_m3_m3,
        .minimum_water_content_m3_m3 = minimum_water_content_m3_m3,
    };
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.NonFiniteSoilAirPorosityInput;
    inline for (.{
        inputs.micropore_capacity_m3,
        inputs.micropore_water_m3,
        inputs.micropore_ice_m3,
        inputs.macropore_capacity_m3,
        inputs.macropore_water_m3,
        inputs.macropore_ice_m3,
        inputs.bulk_soil_volume_m3,
        inputs.negligible_bulk_volume_m3,
        inputs.bulk_density_megagrams_m3,
        inputs.positive_density_threshold_megagrams_m3,
        inputs.zero_water_content_m3_m3,
    }) |value| if (value < 0) return error.InvalidSoilAirPorosityInput;
    if (inputs.hygroscopic_potential_mpa >= 0 or
        inputs.minimum_potential_mpa >= 0 or
        inputs.field_to_wilting_potential_interval == 0)
        return error.InvalidSoilAirPorosityInput;
}

test "air porosity and dry-end contents preserve source expressions" {
    const result = try compute(.{
        .micropore_capacity_m3 = 2,
        .micropore_water_m3 = 0.5,
        .micropore_ice_m3 = 0.25,
        .macropore_capacity_m3 = 1,
        .macropore_water_m3 = 0.2,
        .macropore_ice_m3 = 0.1,
        .bulk_soil_volume_m3 = 5,
        .negligible_bulk_volume_m3 = 1e-12,
        .bulk_density_megagrams_m3 = 1,
        .positive_density_threshold_megagrams_m3 = 0,
        .log_field_capacity_potential = 1,
        .log_field_capacity_to_wilting_point = 0.7,
        .field_to_wilting_potential_interval = 2,
        .log_field_capacity = @log(0.3),
        .hygroscopic_potential_mpa = -10,
        .minimum_potential_mpa = -100,
        .zero_water_content_m3_m3 = 1e-12,
    });
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.95),
        result.air_filled_pore_volume_m3,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.39),
        result.air_filled_porosity_m3_m3,
        1e-15,
    );
    try std.testing.expect(std.math.isFinite(result.hygroscopic_water_content_m3_m3));
}

test "zero-density branch publishes explicit zero-water content" {
    const result = try compute(.{
        .micropore_capacity_m3 = 0,
        .micropore_water_m3 = 0,
        .micropore_ice_m3 = 0,
        .macropore_capacity_m3 = 0,
        .macropore_water_m3 = 0,
        .macropore_ice_m3 = 0,
        .bulk_soil_volume_m3 = 0,
        .negligible_bulk_volume_m3 = 1e-12,
        .bulk_density_megagrams_m3 = 0,
        .positive_density_threshold_megagrams_m3 = 0,
        .log_field_capacity_potential = 1,
        .log_field_capacity_to_wilting_point = 1,
        .field_to_wilting_potential_interval = 1,
        .log_field_capacity = 1,
        .hygroscopic_potential_mpa = -1,
        .minimum_potential_mpa = -2,
        .zero_water_content_m3_m3 = 1e-10,
    });
    try std.testing.expectEqual(@as(f64, 1e-10), result.hygroscopic_water_content_m3_m3);
    try std.testing.expectEqual(@as(f64, 1e-10), result.minimum_water_content_m3_m3);
}
