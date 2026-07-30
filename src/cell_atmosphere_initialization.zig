const std = @import("std");

pub const Inputs = struct {
    mean_annual_air_temperature_c: []const f64,
    atmospheric_co2_umol_per_mol: []const f64,
    atmospheric_ch4_umol_per_mol: []const f64,
    atmospheric_o2_umol_per_mol: []const f64,
    kelvin_offset: f64,
};

pub const State = struct {
    mean_annual_air_temperature_c: []f64,
    mean_annual_deep_soil_temperature_c: []f64,
    mean_annual_air_temperature_k: []f64,
    mean_annual_deep_soil_temperature_k: []f64,
    canopy_air_temperature_k: []f64,
    ground_reference_air_temperature_k: []f64,
    canopy_air_temperature_c: []f64,
    canopy_vapor_pressure_kpa: []f64,
    ground_reference_vapor_pressure_kpa: []f64,
    canopy_co2_umol_per_mol: []f64,
    canopy_ch4_umol_per_mol: []f64,
    canopy_o2_umol_per_mol: []f64,
    canopy_oxygen_exchange_g_o: []f64,
};

/// Exact source-order translation of legacy `STARTS` lines 397--409.
pub fn initialize(state: State, inputs: Inputs) !void {
    const cell_count = inputs.mean_annual_air_temperature_c.len;
    if (cell_count == 0) return error.InvalidCellAtmosphereDimensions;
    inline for (.{
        inputs.atmospheric_co2_umol_per_mol,
        inputs.atmospheric_ch4_umol_per_mol,
        inputs.atmospheric_o2_umol_per_mol,
    }) |values| {
        if (values.len != cell_count)
            return error.CellAtmosphereDimensionMismatch;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (@field(state, field.name).len != cell_count)
            return error.CellAtmosphereDimensionMismatch;
    }
    if (!std.math.isFinite(inputs.kelvin_offset) or inputs.kelvin_offset <= 0)
        return error.InvalidKelvinOffset;
    for (0..cell_count) |cell| {
        const temperature_c = inputs.mean_annual_air_temperature_c[cell];
        if (!std.math.isFinite(temperature_c) or
            temperature_c <= -inputs.kelvin_offset)
            return error.InvalidMeanAnnualTemperature;
        inline for (.{
            inputs.atmospheric_co2_umol_per_mol[cell],
            inputs.atmospheric_ch4_umol_per_mol[cell],
            inputs.atmospheric_o2_umol_per_mol[cell],
        }) |mixing_ratio| {
            if (!std.math.isFinite(mixing_ratio))
                return error.NonFiniteAtmosphericMixingRatio;
            if (mixing_ratio < 0)
                return error.InvalidAtmosphericMixingRatio;
        }
        if (!std.math.isFinite(temperature_c + inputs.kelvin_offset))
            return error.CellAtmosphereTemperatureOverflow;
    }

    for (0..cell_count) |cell| {
        state.mean_annual_air_temperature_c[cell] =
            inputs.mean_annual_air_temperature_c[cell];
        state.mean_annual_deep_soil_temperature_c[cell] =
            state.mean_annual_air_temperature_c[cell];
        state.mean_annual_air_temperature_k[cell] =
            state.mean_annual_air_temperature_c[cell] + inputs.kelvin_offset;
        state.mean_annual_deep_soil_temperature_k[cell] =
            state.mean_annual_deep_soil_temperature_c[cell] +
            inputs.kelvin_offset;
        state.canopy_air_temperature_k[cell] =
            state.mean_annual_air_temperature_k[cell];
        state.ground_reference_air_temperature_k[cell] =
            state.mean_annual_air_temperature_k[cell];
        state.canopy_air_temperature_c[cell] =
            state.mean_annual_air_temperature_c[cell];
        state.canopy_vapor_pressure_kpa[cell] = 0.0;
        state.ground_reference_vapor_pressure_kpa[cell] = 0.0;
        state.canopy_co2_umol_per_mol[cell] =
            inputs.atmospheric_co2_umol_per_mol[cell];
        state.canopy_ch4_umol_per_mol[cell] =
            inputs.atmospheric_ch4_umol_per_mol[cell];
        state.canopy_o2_umol_per_mol[cell] =
            inputs.atmospheric_o2_umol_per_mol[cell];
        state.canopy_oxygen_exchange_g_o[cell] = 0.0;
    }
}

test "STARTS initializes per-cell canopy atmosphere in source order" {
    var fields = [_]f64{99.0} ** (13 * 2);
    const state: State = .{
        .mean_annual_air_temperature_c = fields[0..2],
        .mean_annual_deep_soil_temperature_c = fields[2..4],
        .mean_annual_air_temperature_k = fields[4..6],
        .mean_annual_deep_soil_temperature_k = fields[6..8],
        .canopy_air_temperature_k = fields[8..10],
        .ground_reference_air_temperature_k = fields[10..12],
        .canopy_air_temperature_c = fields[12..14],
        .canopy_vapor_pressure_kpa = fields[14..16],
        .ground_reference_vapor_pressure_kpa = fields[16..18],
        .canopy_co2_umol_per_mol = fields[18..20],
        .canopy_ch4_umol_per_mol = fields[20..22],
        .canopy_o2_umol_per_mol = fields[22..24],
        .canopy_oxygen_exchange_g_o = fields[24..26],
    };
    try initialize(state, .{
        .mean_annual_air_temperature_c = &.{ 10.0, -5.0 },
        .atmospheric_co2_umol_per_mol = &.{ 400.0, 420.0 },
        .atmospheric_ch4_umol_per_mol = &.{ 1.8, 1.9 },
        .atmospheric_o2_umol_per_mol = &.{ 209_500.0, 209_400.0 },
        .kelvin_offset = 273.15,
    });

    try std.testing.expectEqualSlices(f64, &.{ 10, -5 }, state.mean_annual_air_temperature_c);
    try std.testing.expectEqualSlices(f64, &.{ 283.15, 268.15 }, state.mean_annual_air_temperature_k);
    try std.testing.expectEqualSlices(f64, &.{ 283.15, 268.15 }, state.canopy_air_temperature_k);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0 }, state.canopy_vapor_pressure_kpa);
    try std.testing.expectEqualSlices(f64, &.{ 400, 420 }, state.canopy_co2_umol_per_mol);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0 }, state.canopy_oxygen_exchange_g_o);
}

test "late invalid gas leaves every atmospheric field unchanged" {
    var fields = [_]f64{7.0} ** (13 * 2);
    const before = fields;
    const state: State = .{
        .mean_annual_air_temperature_c = fields[0..2],
        .mean_annual_deep_soil_temperature_c = fields[2..4],
        .mean_annual_air_temperature_k = fields[4..6],
        .mean_annual_deep_soil_temperature_k = fields[6..8],
        .canopy_air_temperature_k = fields[8..10],
        .ground_reference_air_temperature_k = fields[10..12],
        .canopy_air_temperature_c = fields[12..14],
        .canopy_vapor_pressure_kpa = fields[14..16],
        .ground_reference_vapor_pressure_kpa = fields[16..18],
        .canopy_co2_umol_per_mol = fields[18..20],
        .canopy_ch4_umol_per_mol = fields[20..22],
        .canopy_o2_umol_per_mol = fields[22..24],
        .canopy_oxygen_exchange_g_o = fields[24..26],
    };
    try std.testing.expectError(
        error.NonFiniteAtmosphericMixingRatio,
        initialize(state, .{
            .mean_annual_air_temperature_c = &.{ 10, 10 },
            .atmospheric_co2_umol_per_mol = &.{ 400, 400 },
            .atmospheric_ch4_umol_per_mol = &.{ 1.8, 1.8 },
            .atmospheric_o2_umol_per_mol = &.{ 209_500, std.math.nan(f64) },
            .kelvin_offset = 273.15,
        }),
    );
    try std.testing.expectEqualSlices(f64, &before, &fields);
}
