const std = @import("std");

pub const RootProfile = enum {
    shallow,
    non_shallow,
};

pub const Inputs = struct {
    preliminary_carboxylation_umol_per_m2_s: f64,
    negligible_carboxylation_umol_per_m2_s: f64,
    atmospheric_to_intercellular_co2_umol_per_m3: f64,
    minimum_stomatal_resistance_s_per_m: f64,
    cuticular_resistance_s_per_m: f64,
    stomatal_turgor_response: f64,
    air_amount_mol_per_m3: f64,
    root_profile: RootProfile,
    shallow_root_growth_water_response: f64,
};

pub const Response = struct {
    admitted: bool,
    zero_water_stress_resistance_s_per_m: f64,
    current_resistance_s_per_m: f64,
    stomatal_conductance_mol_per_m2_s: f64,
    c3_water_response: f64,
    c4_water_response: f64,
};

/// Exact grosub.f lines 1026--1053, repeated at 1206--1234, 1412--1435, and
/// 1546--1569, for one illuminated leaf sample. The same operation applies to
/// direct and diffuse radiation and to both C4 mesophyll and C3 mesophyll
/// paths; C4 additionally copies the C3 response into WFN4.
pub fn calculate(inputs: Inputs) !Response {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == f64) {
            const value = @field(inputs, field.name);
            if (!std.math.isFinite(value))
                return error.NonFiniteCanopyCarboxylationWaterInput;
        }
    }
    if (inputs.preliminary_carboxylation_umol_per_m2_s < 0 or
        inputs.negligible_carboxylation_umol_per_m2_s < 0 or
        inputs.minimum_stomatal_resistance_s_per_m <= 0 or
        inputs.cuticular_resistance_s_per_m <
            inputs.minimum_stomatal_resistance_s_per_m or
        inputs.stomatal_turgor_response < 0 or
        inputs.air_amount_mol_per_m3 <= 0 or
        inputs.shallow_root_growth_water_response < 0)
        return error.InvalidCanopyCarboxylationWaterInput;
    if (inputs.preliminary_carboxylation_umol_per_m2_s <=
        inputs.negligible_carboxylation_umol_per_m2_s)
        return .{
            .admitted = false,
            .zero_water_stress_resistance_s_per_m = 0,
            .current_resistance_s_per_m = 0,
            .stomatal_conductance_mol_per_m2_s = 0,
            .c3_water_response = 0,
            .c4_water_response = 0,
        };

    const zero_water_stress_resistance_s_per_m = @min(
        inputs.cuticular_resistance_s_per_m,
        @max(
            inputs.minimum_stomatal_resistance_s_per_m,
            inputs.atmospheric_to_intercellular_co2_umol_per_m3 /
                inputs.preliminary_carboxylation_umol_per_m2_s,
        ),
    );
    const current_resistance_s_per_m =
        zero_water_stress_resistance_s_per_m +
        (inputs.cuticular_resistance_s_per_m -
            zero_water_stress_resistance_s_per_m) *
            inputs.stomatal_turgor_response;
    const stomatal_conductance_mol_per_m2_s =
        inputs.air_amount_mol_per_m3 / current_resistance_s_per_m;
    const water_response = switch (inputs.root_profile) {
        .non_shallow => std.math.pow(
            f64,
            zero_water_stress_resistance_s_per_m /
                current_resistance_s_per_m,
            0.667,
        ),
        .shallow => inputs.shallow_root_growth_water_response,
    };
    inline for (.{ current_resistance_s_per_m, stomatal_conductance_mol_per_m2_s, water_response }) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteCanopyCarboxylationWaterResponse;
    return .{
        .admitted = true,
        .zero_water_stress_resistance_s_per_m = zero_water_stress_resistance_s_per_m,
        .current_resistance_s_per_m = current_resistance_s_per_m,
        .stomatal_conductance_mol_per_m2_s = stomatal_conductance_mol_per_m2_s,
        .c3_water_response = water_response,
        .c4_water_response = water_response,
    };
}

fn sampleInputs() Inputs {
    return .{
        .preliminary_carboxylation_umol_per_m2_s = 20,
        .negligible_carboxylation_umol_per_m2_s = 1.0e-12,
        .atmospheric_to_intercellular_co2_umol_per_m3 = 2000,
        .minimum_stomatal_resistance_s_per_m = 50,
        .cuticular_resistance_s_per_m = 5000,
        .stomatal_turgor_response = 0.25,
        .air_amount_mol_per_m3 = 40,
        .root_profile = .non_shallow,
        .shallow_root_growth_water_response = 0.7,
    };
}

test "GROSUB deeper-root water response preserves source operation order" {
    const response = try calculate(sampleInputs());
    const source_resistance: f64 = 100;
    const current_resistance = source_resistance + (5000 - source_resistance) * 0.25;
    try std.testing.expect(response.admitted);
    try std.testing.expectEqual(source_resistance, response.zero_water_stress_resistance_s_per_m);
    try std.testing.expectEqual(current_resistance, response.current_resistance_s_per_m);
    try std.testing.expectEqual(40 / current_resistance, response.stomatal_conductance_mol_per_m2_s);
    try std.testing.expectEqual(
        std.math.pow(f64, source_resistance / current_resistance, 0.667),
        response.c4_water_response,
    );
    try std.testing.expectEqual(response.c3_water_response, response.c4_water_response);
}

test "GROSUB shallow roots use canopy growth water response" {
    var inputs = sampleInputs();
    inputs.root_profile = .shallow;
    const response = try calculate(inputs);
    try std.testing.expectEqual(@as(f64, 0.7), response.c3_water_response);
    try std.testing.expectEqual(@as(f64, 0.7), response.c4_water_response);
}

test "source permits responses above one and negative CO2 gradient before resistance bounding" {
    var inputs = sampleInputs();
    inputs.atmospheric_to_intercellular_co2_umol_per_m3 = -1;
    inputs.stomatal_turgor_response = 2;
    const non_shallow = try calculate(inputs);
    try std.testing.expectEqual(
        inputs.minimum_stomatal_resistance_s_per_m,
        non_shallow.zero_water_stress_resistance_s_per_m,
    );
    const expected_current = inputs.minimum_stomatal_resistance_s_per_m +
        (inputs.cuticular_resistance_s_per_m -
            inputs.minimum_stomatal_resistance_s_per_m) * 2;
    try std.testing.expectEqual(expected_current, non_shallow.current_resistance_s_per_m);

    inputs.root_profile = .shallow;
    inputs.shallow_root_growth_water_response = 1.2;
    const shallow = try calculate(inputs);
    try std.testing.expectEqual(@as(f64, 1.2), shallow.c3_water_response);
    try std.testing.expectEqual(@as(f64, 1.2), shallow.c4_water_response);
}

test "preliminary rate comparison remains strict at source threshold" {
    var inputs = sampleInputs();
    inputs.preliminary_carboxylation_umol_per_m2_s =
        inputs.negligible_carboxylation_umol_per_m2_s;
    const response = try calculate(inputs);
    try std.testing.expect(!response.admitted);
}

test "invalid canopy water-response input fails explicitly" {
    var inputs = sampleInputs();
    inputs.cuticular_resistance_s_per_m = 49;
    try std.testing.expectError(
        error.InvalidCanopyCarboxylationWaterInput,
        calculate(inputs),
    );
    inputs = sampleInputs();
    inputs.air_amount_mol_per_m3 = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteCanopyCarboxylationWaterInput,
        calculate(inputs),
    );
}
