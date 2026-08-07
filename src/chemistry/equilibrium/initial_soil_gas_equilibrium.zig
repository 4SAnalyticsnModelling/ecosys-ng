const std = @import("std");

pub const GasValues = struct {
    carbon_dioxide: f64,
    methane: f64,
    oxygen: f64,
    nitrogen: f64,
    nitrous_oxide: f64,
};

pub const TemperatureResponse = struct {
    intercept: GasValues,
    slope_per_c: GasValues,
};

pub const Parameters = struct {
    water_to_air_ratio_at_25_c: GasValues,
    activity_exponent: GasValues,
    temperature_response: TemperatureResponse,
    carbon_dioxide_to_bicarbonate_mol_per_m3: f64,
    bicarbonate_to_carbonate_mol_per_m3: f64,
};

pub const Inputs = struct {
    atmospheric_mol_per_m3: GasValues,
    mean_annual_air_temperature_c: f64,
    hydrogen_mol_per_m3: f64,
};

pub const Result = struct {
    aqueous_carbon_dioxide_mol_per_m3: f64,
    aqueous_carbonate_mol_per_m3: f64,
    aqueous_bicarbonate_mol_per_m3: f64,
    aqueous_methane_mol_per_m3: f64,
    aqueous_oxygen_mol_per_m3: f64,
    aqueous_nitrogen_mol_per_m3: f64,
    aqueous_nitrous_oxide_mol_per_m3: f64,
};

/// Direct translation of `starte.f` lines 234--250 and carbonate lines 241--242.
/// The source's dynamic-salt branch multiplies CO2 by exactly one, so both
/// branches intentionally share this calculation.
pub fn calculate(inputs: Inputs, parameters: Parameters) !Result {
    try validateGasValues(inputs.atmospheric_mol_per_m3);
    try validateGasValues(parameters.water_to_air_ratio_at_25_c);
    try validateFiniteGasValues(parameters.activity_exponent);
    try validateFiniteGasValues(parameters.temperature_response.intercept);
    try validateFiniteGasValues(parameters.temperature_response.slope_per_c);
    inline for (.{
        inputs.mean_annual_air_temperature_c,
        inputs.hydrogen_mol_per_m3,
        parameters.carbon_dioxide_to_bicarbonate_mol_per_m3,
        parameters.bicarbonate_to_carbonate_mol_per_m3,
    }) |value| if (!std.math.isFinite(value))
        return error.NonFiniteInitialSoilGasInput;
    if (inputs.mean_annual_air_temperature_c <= -273.15 or
        inputs.hydrogen_mol_per_m3 <= 0 or
        parameters.carbon_dioxide_to_bicarbonate_mol_per_m3 < 0 or
        parameters.bicarbonate_to_carbonate_mol_per_m3 < 0)
        return error.InvalidInitialSoilGasInput;

    const aqueous_carbon_dioxide = dissolved(
        inputs.atmospheric_mol_per_m3.carbon_dioxide,
        parameters.water_to_air_ratio_at_25_c.carbon_dioxide,
        parameters.activity_exponent.carbon_dioxide,
        parameters.temperature_response.intercept.carbon_dioxide,
        parameters.temperature_response.slope_per_c.carbon_dioxide,
        inputs.mean_annual_air_temperature_c,
    );
    const result: Result = .{
        .aqueous_carbon_dioxide_mol_per_m3 = aqueous_carbon_dioxide,
        .aqueous_carbonate_mol_per_m3 = aqueous_carbon_dioxide *
            parameters.carbon_dioxide_to_bicarbonate_mol_per_m3 *
            parameters.bicarbonate_to_carbonate_mol_per_m3 /
            (inputs.hydrogen_mol_per_m3 * inputs.hydrogen_mol_per_m3),
        .aqueous_bicarbonate_mol_per_m3 = aqueous_carbon_dioxide *
            parameters.carbon_dioxide_to_bicarbonate_mol_per_m3 /
            inputs.hydrogen_mol_per_m3,
        .aqueous_methane_mol_per_m3 = dissolvedGas(.methane, inputs, parameters),
        .aqueous_oxygen_mol_per_m3 = dissolvedGas(.oxygen, inputs, parameters),
        .aqueous_nitrogen_mol_per_m3 = dissolvedGas(.nitrogen, inputs, parameters),
        .aqueous_nitrous_oxide_mol_per_m3 = dissolvedGas(.nitrous_oxide, inputs, parameters),
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        const value = @field(result, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteInitialSoilGasResult;
        if (value < 0) return error.InvalidInitialSoilGasResult;
    }
    return result;
}

const Gas = enum { methane, oxygen, nitrogen, nitrous_oxide };

fn dissolvedGas(comptime gas: Gas, inputs: Inputs, parameters: Parameters) f64 {
    const name = @tagName(gas);
    return dissolved(
        @field(inputs.atmospheric_mol_per_m3, name),
        @field(parameters.water_to_air_ratio_at_25_c, name),
        @field(parameters.activity_exponent, name),
        @field(parameters.temperature_response.intercept, name),
        @field(parameters.temperature_response.slope_per_c, name),
        inputs.mean_annual_air_temperature_c,
    );
}

fn dissolved(
    atmospheric_mol_per_m3: f64,
    water_to_air_ratio_at_25_c: f64,
    activity_exponent: f64,
    temperature_intercept: f64,
    temperature_slope_per_c: f64,
    temperature_c: f64,
) f64 {
    return atmospheric_mol_per_m3 * water_to_air_ratio_at_25_c /
        @exp(activity_exponent) *
        @exp(temperature_intercept - temperature_slope_per_c * temperature_c);
}

fn validateGasValues(values: GasValues) !void {
    try validateFiniteGasValues(values);
    inline for (@typeInfo(GasValues).@"struct".fields) |field|
        if (@field(values, field.name) < 0) return error.InvalidInitialSoilGasInput;
}

fn validateFiniteGasValues(values: GasValues) !void {
    inline for (@typeInfo(GasValues).@"struct".fields) |field|
        if (!std.math.isFinite(@field(values, field.name)))
            return error.NonFiniteInitialSoilGasInput;
}

fn sourceParameters() Parameters {
    return .{
        .water_to_air_ratio_at_25_c = .{
            .carbon_dioxide = 0.7391,
            .methane = 0.03156,
            .oxygen = 0.02925,
            .nitrogen = 0.01510,
            .nitrous_oxide = 0.5241,
        },
        .activity_exponent = .{
            .carbon_dioxide = 0.14,
            .methane = 0.14,
            .oxygen = 0.31,
            .nitrogen = 0.23,
            .nitrous_oxide = 0.23,
        },
        .temperature_response = .{
            .intercept = .{
                .carbon_dioxide = 0.843,
                .methane = 0.597,
                .oxygen = 0.516,
                .nitrogen = 0.456,
                .nitrous_oxide = 0.897,
            },
            .slope_per_c = .{
                .carbon_dioxide = 0.0281,
                .methane = 0.0199,
                .oxygen = 0.0172,
                .nitrogen = 0.0152,
                .nitrous_oxide = 0.0299,
            },
        },
        .carbon_dioxide_to_bicarbonate_mol_per_m3 = 4.2e-4,
        .bicarbonate_to_carbonate_mol_per_m3 = 5.6e-8,
    };
}

test "STARTE initial gases preserve source operation order" {
    const inputs: Inputs = .{
        .atmospheric_mol_per_m3 = .{
            .carbon_dioxide = 1,
            .methane = 2,
            .oxygen = 3,
            .nitrogen = 4,
            .nitrous_oxide = 5,
        },
        .mean_annual_air_temperature_c = 10,
        .hydrogen_mol_per_m3 = 1.0e-4,
    };
    const parameters = sourceParameters();
    const result = try calculate(inputs, parameters);
    const carbon_dioxide = 0.7391 / @exp(0.14) *
        @exp(0.843 - 0.0281 * 10);
    try std.testing.expectApproxEqAbs(
        carbon_dioxide,
        result.aqueous_carbon_dioxide_mol_per_m3,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        carbon_dioxide * 4.2e-4 * 5.6e-8 / 1.0e-8,
        result.aqueous_carbonate_mol_per_m3,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        carbon_dioxide * 4.2e-4 / 1.0e-4,
        result.aqueous_bicarbonate_mol_per_m3,
        1.0e-15,
    );
    try std.testing.expectApproxEqAbs(
        4 * 0.01510 / @exp(0.23) * @exp(0.456 - 0.0152 * 10),
        result.aqueous_nitrogen_mol_per_m3,
        1.0e-15,
    );
}

test "STARTE initial gases accept runtime parameter values" {
    const unit: GasValues = .{
        .carbon_dioxide = 1,
        .methane = 1,
        .oxygen = 1,
        .nitrogen = 1,
        .nitrous_oxide = 1,
    };
    const zero: GasValues = .{
        .carbon_dioxide = 0,
        .methane = 0,
        .oxygen = 0,
        .nitrogen = 0,
        .nitrous_oxide = 0,
    };
    const result = try calculate(.{
        .atmospheric_mol_per_m3 = unit,
        .mean_annual_air_temperature_c = 20,
        .hydrogen_mol_per_m3 = 2,
    }, .{
        .water_to_air_ratio_at_25_c = unit,
        .activity_exponent = zero,
        .temperature_response = .{ .intercept = zero, .slope_per_c = zero },
        .carbon_dioxide_to_bicarbonate_mol_per_m3 = 2,
        .bicarbonate_to_carbonate_mol_per_m3 = 3,
    });
    try std.testing.expectEqual(@as(f64, 1), result.aqueous_carbon_dioxide_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), result.aqueous_bicarbonate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1.5), result.aqueous_carbonate_mol_per_m3);
    try std.testing.expectEqual(@as(f64, 1), result.aqueous_nitrous_oxide_mol_per_m3);
}

test "STARTE initial gases reject invalid hydrogen" {
    const unit: GasValues = .{
        .carbon_dioxide = 1,
        .methane = 1,
        .oxygen = 1,
        .nitrogen = 1,
        .nitrous_oxide = 1,
    };
    var inputs: Inputs = .{
        .atmospheric_mol_per_m3 = unit,
        .mean_annual_air_temperature_c = 20,
        .hydrogen_mol_per_m3 = 0,
    };
    try std.testing.expectError(
        error.InvalidInitialSoilGasInput,
        calculate(inputs, sourceParameters()),
    );
    inputs.hydrogen_mol_per_m3 = std.math.nan(f64);
    try std.testing.expectError(
        error.NonFiniteInitialSoilGasInput,
        calculate(inputs, sourceParameters()),
    );
}
