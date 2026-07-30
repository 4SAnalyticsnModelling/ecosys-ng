const std = @import("std");

pub const Inputs = struct {
    top_soil_horizontal_area_m2: []const f64,
    deepest_soil_bottom_depth_m: []const f64,
    mean_annual_soil_temperature_c: []const f64,
    mean_annual_soil_temperature_k: []const f64,
};

pub const Parameters = struct {
    minimum_snow_water_heat_capacity_mj_per_m2_k: f64,
    minimum_surface_litter_heat_capacity_mj_per_m2_k: f64,
    minimum_soil_water_heat_capacity_mj_per_m2_k: f64,
    minimum_source_depth_m: f64,
    source_depth_below_profile_m: f64,
    lower_boundary_conductivity_m_mj_per_h_k: f64,
    geothermal_flux_mj_per_m2_h: f64,
};

pub const State = struct {
    minimum_snow_water_heat_capacity_mj_per_k: []f64,
    minimum_surface_litter_heat_capacity_mj_per_k: []f64,
    minimum_soil_water_heat_capacity_mj_per_k: []f64,
    lower_heat_source_depth_m: []f64,
    surface_litter_temperature_c: []f64,
    surface_litter_temperature_k: []f64,
    deep_source_temperature_k: []f64,
};

/// Exact source-order translation of legacy `STARTS` lines 654--661.
pub fn initialize(
    state: State,
    inputs: Inputs,
    parameters: Parameters,
) !void {
    const cell_count = inputs.top_soil_horizontal_area_m2.len;
    if (cell_count == 0) return error.InvalidLowerHeatBoundaryDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (@field(inputs, field.name).len != cell_count)
            return error.LowerHeatBoundaryDimensionMismatch;
    }
    inline for (@typeInfo(State).@"struct".fields) |field| {
        if (@field(state, field.name).len != cell_count)
            return error.LowerHeatBoundaryDimensionMismatch;
    }
    inline for (@typeInfo(Parameters).@"struct".fields) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value))
            return error.NonFiniteLowerHeatBoundaryParameter;
    }
    if (parameters.minimum_snow_water_heat_capacity_mj_per_m2_k < 0 or
        parameters.minimum_surface_litter_heat_capacity_mj_per_m2_k < 0 or
        parameters.minimum_soil_water_heat_capacity_mj_per_m2_k < 0 or
        parameters.minimum_source_depth_m <= 0 or
        parameters.source_depth_below_profile_m < 0 or
        parameters.lower_boundary_conductivity_m_mj_per_h_k <= 0)
        return error.InvalidLowerHeatBoundaryParameter;

    for (0..cell_count) |cell| {
        inline for (.{
            inputs.top_soil_horizontal_area_m2[cell],
            inputs.deepest_soil_bottom_depth_m[cell],
            inputs.mean_annual_soil_temperature_c[cell],
            inputs.mean_annual_soil_temperature_k[cell],
        }) |value| if (!std.math.isFinite(value))
            return error.NonFiniteLowerHeatBoundaryInput;
        if (inputs.top_soil_horizontal_area_m2[cell] <= 0 or
            inputs.deepest_soil_bottom_depth_m[cell] < 0 or
            inputs.mean_annual_soil_temperature_k[cell] <= 0)
            return error.InvalidLowerHeatBoundaryInput;
        const source_depth_m = @max(
            parameters.minimum_source_depth_m,
            inputs.deepest_soil_bottom_depth_m[cell] +
                parameters.source_depth_below_profile_m,
        );
        inline for (.{
            parameters.minimum_snow_water_heat_capacity_mj_per_m2_k *
                inputs.top_soil_horizontal_area_m2[cell],
            parameters.minimum_surface_litter_heat_capacity_mj_per_m2_k *
                inputs.top_soil_horizontal_area_m2[cell],
            parameters.minimum_soil_water_heat_capacity_mj_per_m2_k *
                inputs.top_soil_horizontal_area_m2[cell],
            source_depth_m,
            inputs.mean_annual_soil_temperature_k[cell] +
                parameters.geothermal_flux_mj_per_m2_h * source_depth_m /
                    parameters.lower_boundary_conductivity_m_mj_per_h_k,
        }) |candidate| if (!std.math.isFinite(candidate))
            return error.LowerHeatBoundaryOverflow;
    }

    for (0..cell_count) |cell| {
        state.minimum_snow_water_heat_capacity_mj_per_k[cell] =
            parameters.minimum_snow_water_heat_capacity_mj_per_m2_k *
            inputs.top_soil_horizontal_area_m2[cell];
        state.minimum_surface_litter_heat_capacity_mj_per_k[cell] =
            parameters.minimum_surface_litter_heat_capacity_mj_per_m2_k *
            inputs.top_soil_horizontal_area_m2[cell];
        state.minimum_soil_water_heat_capacity_mj_per_k[cell] =
            parameters.minimum_soil_water_heat_capacity_mj_per_m2_k *
            inputs.top_soil_horizontal_area_m2[cell];
        state.lower_heat_source_depth_m[cell] = @max(
            parameters.minimum_source_depth_m,
            inputs.deepest_soil_bottom_depth_m[cell] +
                parameters.source_depth_below_profile_m,
        );
        state.surface_litter_temperature_c[cell] =
            inputs.mean_annual_soil_temperature_c[cell];
        state.surface_litter_temperature_k[cell] =
            inputs.mean_annual_soil_temperature_k[cell];
        state.deep_source_temperature_k[cell] =
            inputs.mean_annual_soil_temperature_k[cell] +
            parameters.geothermal_flux_mj_per_m2_h *
                state.lower_heat_source_depth_m[cell] /
                parameters.lower_boundary_conductivity_m_mj_per_h_k;
    }
}

fn sourceParameters() Parameters {
    return .{
        .minimum_snow_water_heat_capacity_mj_per_m2_k = 8.380e-4,
        .minimum_surface_litter_heat_capacity_mj_per_m2_k = 8.380e-5,
        .minimum_soil_water_heat_capacity_mj_per_m2_k = 8.380e-4,
        .minimum_source_depth_m = 10,
        .source_depth_below_profile_m = 1,
        .lower_boundary_conductivity_m_mj_per_h_k = 8.1e-3,
        .geothermal_flux_mj_per_m2_h = 2.052e-4,
    };
}

test "STARTS initializes minimum capacities and geothermal source" {
    var fields = [_]f64{0.0} ** 7;
    const state: State = .{
        .minimum_snow_water_heat_capacity_mj_per_k = fields[0..1],
        .minimum_surface_litter_heat_capacity_mj_per_k = fields[1..2],
        .minimum_soil_water_heat_capacity_mj_per_k = fields[2..3],
        .lower_heat_source_depth_m = fields[3..4],
        .surface_litter_temperature_c = fields[4..5],
        .surface_litter_temperature_k = fields[5..6],
        .deep_source_temperature_k = fields[6..7],
    };
    try initialize(state, .{
        .top_soil_horizontal_area_m2 = &.{200},
        .deepest_soil_bottom_depth_m = &.{3},
        .mean_annual_soil_temperature_c = &.{10},
        .mean_annual_soil_temperature_k = &.{283.15},
    }, sourceParameters());

    try std.testing.expectApproxEqAbs(@as(f64, 0.1676), fields[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.01676), fields[1], 1e-15);
    try std.testing.expectEqual(@as(f64, 10), fields[3]);
    try std.testing.expectApproxEqAbs(
        @as(f64, 283.15 + 2.052e-4 * 10 / 8.1e-3),
        fields[6],
        1e-13,
    );
}

test "runtime source depth follows deeper profile and late failure is atomic" {
    var fields = [_]f64{7.0} ** 14;
    const before = fields;
    const state: State = .{
        .minimum_snow_water_heat_capacity_mj_per_k = fields[0..2],
        .minimum_surface_litter_heat_capacity_mj_per_k = fields[2..4],
        .minimum_soil_water_heat_capacity_mj_per_k = fields[4..6],
        .lower_heat_source_depth_m = fields[6..8],
        .surface_litter_temperature_c = fields[8..10],
        .surface_litter_temperature_k = fields[10..12],
        .deep_source_temperature_k = fields[12..14],
    };
    try std.testing.expectError(
        error.InvalidLowerHeatBoundaryInput,
        initialize(state, .{
            .top_soil_horizontal_area_m2 = &.{ 100, 100 },
            .deepest_soil_bottom_depth_m = &.{ 15, -1 },
            .mean_annual_soil_temperature_c = &.{ 10, 10 },
            .mean_annual_soil_temperature_k = &.{ 283.15, 283.15 },
        }, sourceParameters()),
    );
    try std.testing.expectEqualSlices(f64, &before, &fields);
}
