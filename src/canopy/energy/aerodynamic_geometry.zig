const std = @import("std");

pub const WindReferenceMode = enum {
    measurement_height_above_displacement,
    at_least_two_m_above_displacement,
};

pub const Parameters = struct {
    canopy_area_threshold_m2: f64,
    height_tolerance_m: f64,
    water_surface_roughness_m: f64,
    snow_heat_capacity_threshold_megajoules_k: f64,
    minimum_air_column_height_m: f64,
    air_heat_capacity_megajoules_m3_k: f64,
    richardson_scale_h2_m: f64,
    resistance_denominator_factor: f64,
    inactive_boundary_layer_resistance_h_m: f64,
};

pub const Inputs = struct {
    ecosystem_type: i32,
    combined_leaf_area_m2: f64,
    combined_stem_area_m2: f64,
    combined_standing_dead_area_m2: f64,
    canopy_geometry_area_m2: f64,
    canopy_air_volume_area_m2: f64,
    canopy_height_m: f64,
    snow_depth_m: f64,
    surface_water_depth_m: f64,
    previous_roughness_height_m: f64,
    ground_surface_roughness_m: f64,
    wind_measurement_height_m: f64,
    wind_reference_mode: WindReferenceMode,
    first_snow_layer_heat_capacity_megajoules_k: f64,
    wind_speed_m_h: f64,
};

pub const Result = struct {
    zero_plane_displacement_m: f64,
    roughness_height_m: f64,
    wind_reference_height_m: f64,
    isothermal_richardson_number: f64,
    isothermal_boundary_layer_resistance_h_m: f64,
    canopy_air_heat_capacity_megajoules_k: f64,
    canopy_air_volume_m3: f64,
};

pub const CalculationError = error{
    NonFiniteInput,
    NegativeArea,
    NonPositiveGeometryArea,
    NonPositiveVolumeArea,
    NegativeHeight,
    InvalidParameter,
    NonPositiveWindSpeed,
    InvalidAerodynamicGeometry,
    NonFiniteResult,
};

/// Translates `hour1.f` lines 4829--4876 for one grid cell.
pub fn calculate(inputs: Inputs, parameters: Parameters) CalculationError!Result {
    try validate(inputs, parameters);

    var zero_plane_displacement_m: f64 = 0.0;
    var updated_roughness_height_m = inputs.previous_roughness_height_m;
    var wind_reference_height_m: f64 = 0.0;
    var richardson_number: f64 = 0.0;
    var boundary_layer_resistance_h_m =
        parameters.inactive_boundary_layer_resistance_h_m;

    if (inputs.ecosystem_type >= 0) {
        const combined_canopy_area_m2 = inputs.combined_leaf_area_m2 +
            inputs.combined_stem_area_m2 + inputs.combined_standing_dead_area_m2;
        var candidate_roughness_height_m: f64 = 0.0;
        if (combined_canopy_area_m2 > parameters.canopy_area_threshold_m2 and
            inputs.canopy_height_m >= inputs.snow_depth_m - parameters.height_tolerance_m and
            inputs.canopy_height_m >=
                inputs.surface_water_depth_m - parameters.height_tolerance_m)
        {
            const canopy_area_index = combined_canopy_area_m2 /
                inputs.canopy_geometry_area_m2;
            const attenuation = @exp(-0.5 * canopy_area_index);
            const intercepted_fraction = 1.0 - attenuation;
            zero_plane_displacement_m = inputs.previous_roughness_height_m +
                inputs.canopy_height_m *
                    @max(0.0, 1.0 - 2.0 / canopy_area_index * intercepted_fraction);

            // HOUR1 computes ZE (intermediate roughness) as:
            // ZE = ZT * max(ZS, exp(-0.5*ARLSG) * (1 - exp(-0.5*ARLSG))),
            // where ZT is canopy height and ZS is surface roughness.
            candidate_roughness_height_m = inputs.canopy_height_m *
                @max(inputs.ground_surface_roughness_m, attenuation * intercepted_fraction);
        } else {
            zero_plane_displacement_m = inputs.previous_roughness_height_m;
        }

        wind_reference_height_m = switch (inputs.wind_reference_mode) {
            .measurement_height_above_displacement => inputs.wind_measurement_height_m + zero_plane_displacement_m,
            .at_least_two_m_above_displacement => @max(inputs.wind_measurement_height_m, zero_plane_displacement_m + 2.0),
        };

        updated_roughness_height_m = if (inputs.first_snow_layer_heat_capacity_megajoules_k >
            parameters.snow_heat_capacity_threshold_megajoules_k)
            @max(candidate_roughness_height_m, parameters.water_surface_roughness_m)
        else
            @max(candidate_roughness_height_m, inputs.ground_surface_roughness_m);

        const aerodynamic_height_m = wind_reference_height_m - zero_plane_displacement_m;
        if (inputs.wind_speed_m_h <= 0.0) return error.NonPositiveWindSpeed;
        if (updated_roughness_height_m <= 0.0 or
            aerodynamic_height_m <= 0.0 or
            aerodynamic_height_m / updated_roughness_height_m <= 0.0)
        {
            return error.InvalidAerodynamicGeometry;
        }
        richardson_number = parameters.richardson_scale_h2_m *
            aerodynamic_height_m / (inputs.wind_speed_m_h * inputs.wind_speed_m_h);
        boundary_layer_resistance_h_m =
            std.math.pow(f64, @log(aerodynamic_height_m / updated_roughness_height_m), 2.0) /
            (parameters.resistance_denominator_factor * inputs.wind_speed_m_h);
    }

    const effective_air_column_height_m =
        @max(parameters.minimum_air_column_height_m, wind_reference_height_m);
    const canopy_air_volume_m3 =
        effective_air_column_height_m * inputs.canopy_air_volume_area_m2;
    const canopy_air_heat_capacity_megajoules_k =
        canopy_air_volume_m3 * parameters.air_heat_capacity_megajoules_m3_k;

    const result = Result{
        .zero_plane_displacement_m = zero_plane_displacement_m,
        .roughness_height_m = updated_roughness_height_m,
        .wind_reference_height_m = wind_reference_height_m,
        .isothermal_richardson_number = richardson_number,
        .isothermal_boundary_layer_resistance_h_m = boundary_layer_resistance_h_m,
        .canopy_air_heat_capacity_megajoules_k = canopy_air_heat_capacity_megajoules_k,
        .canopy_air_volume_m3 = canopy_air_volume_m3,
    };
    inline for (std.meta.fields(Result)) |field| {
        if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteResult;
    }
    return result;
}

fn validate(inputs: Inputs, parameters: Parameters) CalculationError!void {
    inline for (std.meta.fields(Inputs)) |field| {
        if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) {
            return error.NonFiniteInput;
        }
    }
    inline for (std.meta.fields(Parameters)) |field| {
        if (!std.math.isFinite(@field(parameters, field.name))) return error.NonFiniteInput;
    }
    if (inputs.combined_leaf_area_m2 < 0.0 or
        inputs.combined_stem_area_m2 < 0.0 or
        inputs.combined_standing_dead_area_m2 < 0.0)
    {
        return error.NegativeArea;
    }
    if (inputs.canopy_geometry_area_m2 <= 0.0) return error.NonPositiveGeometryArea;
    if (inputs.canopy_air_volume_area_m2 <= 0.0) return error.NonPositiveVolumeArea;
    if (inputs.canopy_height_m < 0.0 or
        inputs.snow_depth_m < 0.0 or
        inputs.surface_water_depth_m < 0.0 or
        inputs.previous_roughness_height_m < 0.0 or
        inputs.ground_surface_roughness_m < 0.0 or
        inputs.wind_measurement_height_m < 0.0 or
        inputs.first_snow_layer_heat_capacity_megajoules_k < 0.0)
    {
        return error.NegativeHeight;
    }
    if (parameters.canopy_area_threshold_m2 < 0.0 or
        parameters.height_tolerance_m < 0.0 or
        parameters.water_surface_roughness_m <= 0.0 or
        parameters.snow_heat_capacity_threshold_megajoules_k < 0.0 or
        parameters.minimum_air_column_height_m <= 0.0 or
        parameters.air_heat_capacity_megajoules_m3_k <= 0.0 or
        parameters.richardson_scale_h2_m <= 0.0 or
        parameters.resistance_denominator_factor <= 0.0 or
        parameters.inactive_boundary_layer_resistance_h_m < 0.0)
    {
        return error.InvalidParameter;
    }
}

fn testParameters() Parameters {
    return .{
        .canopy_area_threshold_m2 = 1.0e-12,
        .height_tolerance_m = 1.0e-12,
        .water_surface_roughness_m = 0.005,
        .snow_heat_capacity_threshold_megajoules_k = 0.1,
        .minimum_air_column_height_m = 5.0,
        .air_heat_capacity_megajoules_m3_k = 1.25e-3,
        .richardson_scale_h2_m = 1.27e8,
        .resistance_denominator_factor = 0.168,
        .inactive_boundary_layer_resistance_h_m = 100.0,
    };
}

test "active canopy preserves displacement roughness and resistance order" {
    const inputs = Inputs{
        .ecosystem_type = 1,
        .combined_leaf_area_m2 = 10.0,
        .combined_stem_area_m2 = 5.0,
        .combined_standing_dead_area_m2 = 5.0,
        .canopy_geometry_area_m2 = 10.0,
        .canopy_air_volume_area_m2 = 12.0,
        .canopy_height_m = 4.0,
        .snow_depth_m = 0.0,
        .surface_water_depth_m = 0.0,
        .previous_roughness_height_m = 0.1,
        .ground_surface_roughness_m = 0.025,
        .wind_measurement_height_m = 10.0,
        .wind_reference_mode = .measurement_height_above_displacement,
        .first_snow_layer_heat_capacity_megajoules_k = 0.0,
        .wind_speed_m_h = 3600.0,
    };
    const result = try calculate(inputs, testParameters());
    const area_index = 2.0;
    const attenuation = @exp(-0.5 * area_index);
    const intercepted = 1.0 - attenuation;
    const expected_displacement =
        0.1 + 4.0 * @max(0.0, 1.0 - 2.0 / area_index * intercepted);
    const expected_candidate_roughness = 4.0 * @max(0.025, attenuation * intercepted);

    try std.testing.expectApproxEqRel(
        expected_displacement,
        result.zero_plane_displacement_m,
        1.0e-14,
    );
    try std.testing.expectApproxEqRel(
        expected_candidate_roughness,
        result.roughness_height_m,
        1.0e-14,
    );
    try std.testing.expectApproxEqRel(
        @as(f64, 1.27e8 * 10.0 / (3600.0 * 3600.0)),
        result.isothermal_richardson_number,
        1.0e-14,
    );
    try std.testing.expectApproxEqRel(
        (10.0 + expected_displacement) * 12.0,
        result.canopy_air_volume_m3,
        1.0e-14,
    );
}

test "inactive ecosystem uses fallback resistance and minimum air column" {
    var inputs = std.mem.zeroes(Inputs);
    inputs.ecosystem_type = -1;
    inputs.canopy_geometry_area_m2 = 1.0;
    inputs.canopy_air_volume_area_m2 = 8.0;
    inputs.ground_surface_roughness_m = 0.025;
    inputs.wind_reference_mode = .at_least_two_m_above_displacement;
    const result = try calculate(inputs, testParameters());

    try std.testing.expectEqual(@as(f64, 0.0), result.isothermal_richardson_number);
    try std.testing.expectEqual(
        @as(f64, 100.0),
        result.isothermal_boundary_layer_resistance_h_m,
    );
    try std.testing.expectEqual(@as(f64, 40.0), result.canopy_air_volume_m3);
    try std.testing.expectEqual(@as(f64, 0.05), result.canopy_air_heat_capacity_megajoules_k);
}
