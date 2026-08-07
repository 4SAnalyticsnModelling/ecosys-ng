const std = @import("std");

pub const Inputs = struct {
    turgor_water_potential_megapascal: f64,
    turgor_response_shape_per_megapascal: f64,
    minimum_stomatal_resistance_h_per_m: f64,
    cuticular_resistance_h_per_m: f64,
};

pub const Result = struct {
    turgor_response_fraction: f64,
    stomatal_resistance_h_per_m: f64,
};

/// UPTAKE.F 972--973. Applies the source exponential turgor response without
/// adding a clamp or reordering its affine resistance calculation.
pub fn calculate(inputs: Inputs) !Result {
    inline for (@typeInfo(Inputs).@"struct".fields) |field|
        if (!std.math.isFinite(@field(inputs, field.name)))
            return error.InvalidCanopyStomatalResistanceInput;
    if (inputs.turgor_water_potential_megapascal < 0 or
        inputs.minimum_stomatal_resistance_h_per_m < 0 or
        inputs.cuticular_resistance_h_per_m < 0)
        return error.InvalidCanopyStomatalResistanceInput;

    const response = @exp(
        inputs.turgor_response_shape_per_megapascal *
            inputs.turgor_water_potential_megapascal,
    );
    const resistance =
        inputs.minimum_stomatal_resistance_h_per_m +
        (inputs.cuticular_resistance_h_per_m -
            inputs.minimum_stomatal_resistance_h_per_m) *
            response;
    if (!std.math.isFinite(response) or !std.math.isFinite(resistance))
        return error.NonFiniteCanopyStomatalResistanceResult;
    return .{
        .turgor_response_fraction = response,
        .stomatal_resistance_h_per_m = resistance,
    };
}

test "UPTAKE canopy stomatal resistance preserves source arithmetic" {
    const inputs = Inputs{
        .turgor_water_potential_megapascal = 0.8,
        .turgor_response_shape_per_megapascal = -2,
        .minimum_stomatal_resistance_h_per_m = 0.01,
        .cuticular_resistance_h_per_m = 1.5,
    };
    const result = try calculate(inputs);
    const response = @exp(
        inputs.turgor_response_shape_per_megapascal *
            inputs.turgor_water_potential_megapascal,
    );
    try std.testing.expectEqual(response, result.turgor_response_fraction);
    try std.testing.expectEqual(
        inputs.minimum_stomatal_resistance_h_per_m +
            (inputs.cuticular_resistance_h_per_m -
                inputs.minimum_stomatal_resistance_h_per_m) *
                response,
        result.stomatal_resistance_h_per_m,
    );
}

test "zero turgor selects the source cuticular endpoint" {
    const result = try calculate(.{
        .turgor_water_potential_megapascal = 0,
        .turgor_response_shape_per_megapascal = -2,
        .minimum_stomatal_resistance_h_per_m = 0.01,
        .cuticular_resistance_h_per_m = 1.5,
    });
    try std.testing.expectEqual(@as(f64, 1), result.turgor_response_fraction);
    try std.testing.expectEqual(@as(f64, 1.5), result.stomatal_resistance_h_per_m);
}

test "non-finite exponential result fails explicitly" {
    try std.testing.expectError(
        error.NonFiniteCanopyStomatalResistanceResult,
        calculate(.{
            .turgor_water_potential_megapascal = 1.0e308,
            .turgor_response_shape_per_megapascal = 1,
            .minimum_stomatal_resistance_h_per_m = 0.01,
            .cuticular_resistance_h_per_m = 1.5,
        }),
    );
}
