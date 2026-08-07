const std = @import("std");

pub const Inputs = struct {
    mean_annual_air_temperature_k: []const f64,
    initial_co2_umol_per_mol: []const f64,
    current_co2_umol_per_mol: []const f64,
    methane_umol_per_mol: []const f64,
    oxygen_umol_per_mol: []const f64,
    dinitrogen_umol_per_mol: []const f64,
    nitrous_oxide_umol_per_mol: []const f64,
    ammonia_umol_per_mol: []const f64,
    hydrogen_umol_per_mol: []const f64,
};

pub const Parameters = struct {
    reference_temperature_k: f64,
    carbon_g_per_m3_per_umol_per_mol: f64,
    oxygen_g_per_m3_per_umol_per_mol: f64,
    nitrogen_g_per_m3_per_umol_per_mol: f64,
    ammonia_nitrogen_g_per_m3_per_umol_per_mol: f64,
    hydrogen_g_per_m3_per_umol_per_mol: f64,
};

pub const State = struct {
    initial_co2_g_c_per_m3: []f64,
    current_co2_g_c_per_m3: []f64,
    methane_g_c_per_m3: []f64,
    oxygen_g_o_per_m3: []f64,
    dinitrogen_g_n_per_m3: []f64,
    nitrous_oxide_g_n_per_m3: []f64,
    ammonia_g_n_per_m3: []f64,
    hydrogen_g_h_per_m3: []f64,
};

/// Exact source-order translation of legacy `STARTS` lines 495--502.
pub fn convert(
    state: State,
    inputs: Inputs,
    parameters: Parameters,
) !void {
    const cell_count = inputs.mean_annual_air_temperature_k.len;
    if (cell_count == 0)
        return error.InvalidAtmosphericConcentrationDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (@field(inputs, field.name).len != cell_count)
            return error.AtmosphericConcentrationDimensionMismatch;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (@field(state, field.name).len != cell_count)
            return error.AtmosphericConcentrationDimensionMismatch;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteAtmosphericConversionParameter;
        if (value <= 0) return error.InvalidAtmosphericConversionParameter;
    }
    for (0..cell_count) |cell| {
        const temperature_k = inputs.mean_annual_air_temperature_k[cell];
        if (!std.math.isFinite(temperature_k))
            return error.NonFiniteAtmosphericConcentrationInput;
        if (temperature_k <= 0)
            return error.InvalidAtmosphericTemperature;
        inline for (@typeInfo(Inputs).@"struct".fields) |field| {
            if (!std.mem.eql(
                u8,
                field.name,
                "mean_annual_air_temperature_k",
            )) {
                const value = @field(inputs, field.name)[cell];
                if (!std.math.isFinite(value))
                    return error.NonFiniteAtmosphericConcentrationInput;
                if (value < 0)
                    return error.InvalidAtmosphericMixingRatio;
            }
        }
        const temperature_ratio =
            parameters.reference_temperature_k / temperature_k;
        inline for (.{
            inputs.initial_co2_umol_per_mol[cell] *
                parameters.carbon_g_per_m3_per_umol_per_mol *
                temperature_ratio,
            inputs.current_co2_umol_per_mol[cell] *
                parameters.carbon_g_per_m3_per_umol_per_mol *
                temperature_ratio,
            inputs.methane_umol_per_mol[cell] *
                parameters.carbon_g_per_m3_per_umol_per_mol *
                temperature_ratio,
            inputs.oxygen_umol_per_mol[cell] *
                parameters.oxygen_g_per_m3_per_umol_per_mol *
                temperature_ratio,
            inputs.dinitrogen_umol_per_mol[cell] *
                parameters.nitrogen_g_per_m3_per_umol_per_mol *
                temperature_ratio,
            inputs.nitrous_oxide_umol_per_mol[cell] *
                parameters.nitrogen_g_per_m3_per_umol_per_mol *
                temperature_ratio,
            inputs.ammonia_umol_per_mol[cell] *
                parameters.ammonia_nitrogen_g_per_m3_per_umol_per_mol *
                temperature_ratio,
            inputs.hydrogen_umol_per_mol[cell] *
                parameters.hydrogen_g_per_m3_per_umol_per_mol *
                temperature_ratio,
        }) |candidate| {
            if (!std.math.isFinite(candidate))
                return error.AtmosphericConcentrationOverflow;
        }
    }

    for (0..cell_count) |cell| {
        state.initial_co2_g_c_per_m3[cell] =
            inputs.initial_co2_umol_per_mol[cell] *
            parameters.carbon_g_per_m3_per_umol_per_mol *
            parameters.reference_temperature_k /
            inputs.mean_annual_air_temperature_k[cell];
        state.current_co2_g_c_per_m3[cell] =
            inputs.current_co2_umol_per_mol[cell] *
            parameters.carbon_g_per_m3_per_umol_per_mol *
            parameters.reference_temperature_k /
            inputs.mean_annual_air_temperature_k[cell];
        state.methane_g_c_per_m3[cell] =
            inputs.methane_umol_per_mol[cell] *
            parameters.carbon_g_per_m3_per_umol_per_mol *
            parameters.reference_temperature_k /
            inputs.mean_annual_air_temperature_k[cell];
        state.oxygen_g_o_per_m3[cell] =
            inputs.oxygen_umol_per_mol[cell] *
            parameters.oxygen_g_per_m3_per_umol_per_mol *
            parameters.reference_temperature_k /
            inputs.mean_annual_air_temperature_k[cell];
        state.dinitrogen_g_n_per_m3[cell] =
            inputs.dinitrogen_umol_per_mol[cell] *
            parameters.nitrogen_g_per_m3_per_umol_per_mol *
            parameters.reference_temperature_k /
            inputs.mean_annual_air_temperature_k[cell];
        state.nitrous_oxide_g_n_per_m3[cell] =
            inputs.nitrous_oxide_umol_per_mol[cell] *
            parameters.nitrogen_g_per_m3_per_umol_per_mol *
            parameters.reference_temperature_k /
            inputs.mean_annual_air_temperature_k[cell];
        state.ammonia_g_n_per_m3[cell] =
            inputs.ammonia_umol_per_mol[cell] *
            parameters.ammonia_nitrogen_g_per_m3_per_umol_per_mol *
            parameters.reference_temperature_k /
            inputs.mean_annual_air_temperature_k[cell];
        state.hydrogen_g_h_per_m3[cell] =
            inputs.hydrogen_umol_per_mol[cell] *
            parameters.hydrogen_g_per_m3_per_umol_per_mol *
            parameters.reference_temperature_k /
            inputs.mean_annual_air_temperature_k[cell];
    }
}

fn sourceParameters() Parameters {
    return .{
        .reference_temperature_k = 273.15,
        .carbon_g_per_m3_per_umol_per_mol = 5.36e-4,
        .oxygen_g_per_m3_per_umol_per_mol = 1.43e-3,
        .nitrogen_g_per_m3_per_umol_per_mol = 1.25e-3,
        .ammonia_nitrogen_g_per_m3_per_umol_per_mol = 6.25e-4,
        .hydrogen_g_per_m3_per_umol_per_mol = 8.92e-5,
    };
}

test "STARTS converts all atmospheric gases including both CO2 states" {
    var output = [_]f64{0.0} ** 8;
    try convert(.{
        .initial_co2_g_c_per_m3 = output[0..1],
        .current_co2_g_c_per_m3 = output[1..2],
        .methane_g_c_per_m3 = output[2..3],
        .oxygen_g_o_per_m3 = output[3..4],
        .dinitrogen_g_n_per_m3 = output[4..5],
        .nitrous_oxide_g_n_per_m3 = output[5..6],
        .ammonia_g_n_per_m3 = output[6..7],
        .hydrogen_g_h_per_m3 = output[7..8],
    }, .{
        .mean_annual_air_temperature_k = &.{273.15},
        .initial_co2_umol_per_mol = &.{350},
        .current_co2_umol_per_mol = &.{400},
        .methane_umol_per_mol = &.{2},
        .oxygen_umol_per_mol = &.{210_000},
        .dinitrogen_umol_per_mol = &.{780_000},
        .nitrous_oxide_umol_per_mol = &.{0.33},
        .ammonia_umol_per_mol = &.{0.02},
        .hydrogen_umol_per_mol = &.{0.5},
    }, sourceParameters());

    try std.testing.expectApproxEqAbs(@as(f64, 0.1876), output[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2144), output[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 300.3), output[3], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 975), output[4], 1e-12);
}

test "late invalid cell fails without partial conversion" {
    var output = [_]f64{7.0} ** 16;
    const before = output;
    try std.testing.expectError(
        error.InvalidAtmosphericTemperature,
        convert(.{
            .initial_co2_g_c_per_m3 = output[0..2],
            .current_co2_g_c_per_m3 = output[2..4],
            .methane_g_c_per_m3 = output[4..6],
            .oxygen_g_o_per_m3 = output[6..8],
            .dinitrogen_g_n_per_m3 = output[8..10],
            .nitrous_oxide_g_n_per_m3 = output[10..12],
            .ammonia_g_n_per_m3 = output[12..14],
            .hydrogen_g_h_per_m3 = output[14..16],
        }, .{
            .mean_annual_air_temperature_k = &.{ 273.15, 0 },
            .initial_co2_umol_per_mol = &.{ 350, 350 },
            .current_co2_umol_per_mol = &.{ 400, 400 },
            .methane_umol_per_mol = &.{ 2, 2 },
            .oxygen_umol_per_mol = &.{ 210_000, 210_000 },
            .dinitrogen_umol_per_mol = &.{ 780_000, 780_000 },
            .nitrous_oxide_umol_per_mol = &.{ 0.33, 0.33 },
            .ammonia_umol_per_mol = &.{ 0.02, 0.02 },
            .hydrogen_umol_per_mol = &.{ 0.5, 0.5 },
        }, sourceParameters()),
    );
    try std.testing.expectEqualSlices(f64, &before, &output);
}
