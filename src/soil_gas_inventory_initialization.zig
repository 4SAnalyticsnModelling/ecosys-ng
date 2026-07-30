const std = @import("std");

pub const GasConcentrations = struct {
    carbon_dioxide_g_per_m3: f64,
    methane_g_per_m3: f64,
    oxygen_g_per_m3: f64,
    nitrogen_g_per_m3: f64,
    nitrous_oxide_g_per_m3: f64,
    ammonia_g_n_per_m3: f64,
    hydrogen_g_per_m3: f64,
};

pub const SolubilityParameters = struct {
    ratio_at_25_c: [6]f64,
    ionic_strength_exponent: [6]f64,
    temperature_intercept: [6]f64,
    temperature_slope_per_c: [6]f64,
};

pub const Inputs = struct {
    atmospheric: GasConcentrations,
    air_volume_m3: f64,
    field_capacity_water_m3: f64,
    layer_top_depth_m: f64,
    water_table_depth_m: f64,
    ionic_strength_mol_per_l: f64,
    mean_annual_air_temperature_c: f64,
};

pub const GasContents = struct {
    carbon_dioxide_g: f64,
    methane_g: f64,
    oxygen_g: f64,
    nitrogen_g: f64,
    nitrous_oxide_g: f64,
    ammonia_g_n: f64,
    hydrogen_g: f64,
};

pub const AqueousContents = struct {
    carbon_dioxide_g: f64,
    methane_g: f64,
    oxygen_g: f64,
    nitrogen_g: f64,
    nitrous_oxide_g: f64,
    hydrogen_g: f64,
};

pub const Result = struct {
    gaseous: GasContents,
    aqueous: AqueousContents,
};

/// Direct translation of STARTE lines 1410--1433. Callers apply the cold-start
/// gate (`DATA(20) == NO` and `IGO == 0`) before invoking this pure kernel.
pub fn calculate(inputs: Inputs, parameters: SolubilityParameters) !Result {
    inline for (@typeInfo(GasConcentrations).@"struct".fields) |field| {
        const value = @field(inputs.atmospheric, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteSoilGasInventoryInput;
        if (value < 0) return error.InvalidSoilGasInventoryInput;
    }
    inline for (.{
        inputs.air_volume_m3,
        inputs.field_capacity_water_m3,
        inputs.layer_top_depth_m,
        inputs.water_table_depth_m,
        inputs.ionic_strength_mol_per_l,
        inputs.mean_annual_air_temperature_c,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteSoilGasInventoryInput;
    if (inputs.air_volume_m3 < 0 or
        inputs.field_capacity_water_m3 < 0 or
        inputs.layer_top_depth_m < 0 or
        inputs.water_table_depth_m < 0 or
        inputs.ionic_strength_mol_per_l < 0 or
        inputs.mean_annual_air_temperature_c <= -273.15)
        return error.InvalidSoilGasInventoryInput;
    inline for (.{
        parameters.ratio_at_25_c,
        parameters.ionic_strength_exponent,
        parameters.temperature_intercept,
        parameters.temperature_slope_per_c,
    }) |values| for (values) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteSoilGasInventoryParameter;
    for (parameters.ratio_at_25_c) |value|
        if (value < 0) return error.InvalidSoilGasInventoryParameter;

    const atmospheric_aqueous = [6]f64{
        inputs.atmospheric.carbon_dioxide_g_per_m3,
        inputs.atmospheric.methane_g_per_m3,
        inputs.atmospheric.oxygen_g_per_m3,
        inputs.atmospheric.nitrogen_g_per_m3,
        inputs.atmospheric.nitrous_oxide_g_per_m3,
        inputs.atmospheric.hydrogen_g_per_m3,
    };
    var aqueous: [6]f64 = undefined;
    for (0..aqueous.len) |gas| {
        aqueous[gas] = atmospheric_aqueous[gas] *
            parameters.ratio_at_25_c[gas] /
            @exp(parameters.ionic_strength_exponent[gas] *
                inputs.ionic_strength_mol_per_l) *
            @exp(parameters.temperature_intercept[gas] -
                parameters.temperature_slope_per_c[gas] *
                    inputs.mean_annual_air_temperature_c) *
            inputs.field_capacity_water_m3;
    }
    if (!(inputs.layer_top_depth_m < inputs.water_table_depth_m))
        aqueous[2] = 0;

    const result: Result = .{
        .gaseous = .{
            .carbon_dioxide_g = inputs.atmospheric.carbon_dioxide_g_per_m3 *
                inputs.air_volume_m3,
            .methane_g = inputs.atmospheric.methane_g_per_m3 *
                inputs.air_volume_m3,
            .oxygen_g = inputs.atmospheric.oxygen_g_per_m3 *
                inputs.air_volume_m3,
            .nitrogen_g = inputs.atmospheric.nitrogen_g_per_m3 *
                inputs.air_volume_m3,
            .nitrous_oxide_g = inputs.atmospheric.nitrous_oxide_g_per_m3 *
                inputs.air_volume_m3,
            .ammonia_g_n = inputs.atmospheric.ammonia_g_n_per_m3 *
                inputs.air_volume_m3,
            .hydrogen_g = inputs.atmospheric.hydrogen_g_per_m3 *
                inputs.air_volume_m3,
        },
        .aqueous = .{
            .carbon_dioxide_g = aqueous[0],
            .methane_g = aqueous[1],
            .oxygen_g = aqueous[2],
            .nitrogen_g = aqueous[3],
            .nitrous_oxide_g = aqueous[4],
            .hydrogen_g = aqueous[5],
        },
    };
    inline for (@typeInfo(GasContents).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.gaseous, field.name)))
            return error.NonFiniteSoilGasInventoryResult;
    inline for (@typeInfo(AqueousContents).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.aqueous, field.name)))
            return error.NonFiniteSoilGasInventoryResult;
    return result;
}

fn sourceParameters() SolubilityParameters {
    return .{
        .ratio_at_25_c = .{ 0.7391, 0.03156, 0.02925, 0.01510, 0.5241, 0.03156 },
        .ionic_strength_exponent = .{ 0.14, 0.14, 0.31, 0.23, 0.23, 0.14 },
        .temperature_intercept = .{ 0.843, 0.597, 0.516, 0.456, 0.897, 0.597 },
        .temperature_slope_per_c = .{ 0.0281, 0.0199, 0.0172, 0.0152, 0.0299, 0.0199 },
    };
}

fn testInputs() Inputs {
    return .{
        .atmospheric = .{
            .carbon_dioxide_g_per_m3 = 1,
            .methane_g_per_m3 = 2,
            .oxygen_g_per_m3 = 3,
            .nitrogen_g_per_m3 = 4,
            .nitrous_oxide_g_per_m3 = 5,
            .ammonia_g_n_per_m3 = 6,
            .hydrogen_g_per_m3 = 7,
        },
        .air_volume_m3 = 2,
        .field_capacity_water_m3 = 3,
        .layer_top_depth_m = 0.5,
        .water_table_depth_m = 1,
        .ionic_strength_mol_per_l = 0.1,
        .mean_annual_air_temperature_c = 10,
    };
}

test "STARTE cold-start gas inventories preserve source equations" {
    const inputs = testInputs();
    const parameters = sourceParameters();
    const result = try calculate(inputs, parameters);
    try std.testing.expectEqual(@as(f64, 2), result.gaseous.carbon_dioxide_g);
    try std.testing.expectEqual(@as(f64, 12), result.gaseous.ammonia_g_n);
    try std.testing.expectEqual(@as(f64, 14), result.gaseous.hydrogen_g);
    try std.testing.expectApproxEqAbs(
        3 * 0.02925 / @exp(0.31 * 0.1) *
            @exp(0.516 - 0.0172 * 10) * 3,
        result.aqueous.oxygen_g,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        4 * 0.01510 / @exp(0.23 * 0.1) *
            @exp(0.456 - 0.0152 * 10) * 3,
        result.aqueous.nitrogen_g,
        1.0e-15,
    );
}

test "STARTE water-table gate zeros only aqueous oxygen" {
    var inputs = testInputs();
    inputs.layer_top_depth_m = inputs.water_table_depth_m;
    const result = try calculate(inputs, sourceParameters());
    try std.testing.expectEqual(@as(f64, 0), result.aqueous.oxygen_g);
    try std.testing.expect(result.aqueous.carbon_dioxide_g > 0);
    try std.testing.expect(result.gaseous.oxygen_g > 0);
}

test "STARTE gas initialization rejects non-finite state" {
    var inputs = testInputs();
    inputs.ionic_strength_mol_per_l = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteSoilGasInventoryInput,
        calculate(inputs, sourceParameters()),
    );
}
