const std = @import("std");

pub const Inputs = struct {
    snowpack_heat_capacity_megajoules_k: f64,
    minimum_snowpack_heat_capacity_megajoules_k: f64,
    snow_volume_m3: f64,
    ice_volume_m3: f64,
    water_volume_m3: f64,
    snow_depth_m: f64,
    soil_total_volume_m3: f64,
    soil_water_volume_m3: f64,
    soil_porosity_m3_m3: f64,
    negligible_soil_volume_m3: f64,
    dry_soil_albedo: f64,
    maximum_soil_albedo: f64,
    ground_shortwave_megajoules_per_m2_h: f64,
    ground_par_umol_per_m2_s: f64,
    cell_area_m2: f64,
};

pub const Result = struct {
    ground_albedo: f64,
    backscattered_shortwave_megajoules_per_m2_h_per_sector: f64,
    backscattered_par_umol_per_m2_s_per_sector: f64,
    absorbed_ground_shortwave_megajoules_h: f64,
    absorbed_ground_par_umol_s: f64,
};

/// `hour1.f` lines 1632--1651. Preserves the snow/soil branch and reflected then
/// absorbed radiation assignment order.
pub fn compute(inputs: Inputs) !Result {
    try validate(inputs);
    const ground_albedo = if (inputs.snowpack_heat_capacity_megajoules_k >
        inputs.minimum_snowpack_heat_capacity_megajoules_k)
        try snowAlbedo(inputs)
    else
        soilAlbedo(inputs);
    const backscattered_shortwave_megajoules_per_m2_h_per_sector =
        inputs.ground_shortwave_megajoules_per_m2_h * ground_albedo * 0.25;
    const backscattered_par_umol_per_m2_s_per_sector =
        inputs.ground_par_umol_per_m2_s * ground_albedo * 0.25;
    const absorbed_ground_shortwave_megajoules_h =
        (1.0 - ground_albedo) * inputs.ground_shortwave_megajoules_per_m2_h *
        inputs.cell_area_m2;
    const absorbed_ground_par_umol_s =
        (1.0 - ground_albedo) * inputs.ground_par_umol_per_m2_s *
        inputs.cell_area_m2;
    return .{
        .ground_albedo = ground_albedo,
        .backscattered_shortwave_megajoules_per_m2_h_per_sector = backscattered_shortwave_megajoules_per_m2_h_per_sector,
        .backscattered_par_umol_per_m2_s_per_sector = backscattered_par_umol_per_m2_s_per_sector,
        .absorbed_ground_shortwave_megajoules_h = absorbed_ground_shortwave_megajoules_h,
        .absorbed_ground_par_umol_s = absorbed_ground_par_umol_s,
    };
}

fn snowAlbedo(inputs: Inputs) !f64 {
    const surface_volume_m3 =
        inputs.snow_volume_m3 + inputs.ice_volume_m3 + inputs.water_volume_m3;
    if (surface_volume_m3 <= 0) return error.ZeroSnowSurfaceVolume;
    const snow_surface_albedo =
        (0.90 * inputs.snow_volume_m3 + 0.30 * inputs.ice_volume_m3 +
            0.06 * inputs.water_volume_m3) /
        surface_volume_m3;
    const snow_depth_ratio = inputs.snow_depth_m / 0.07;
    const snow_cover_fraction = @min(
        std.math.pow(f64, snow_depth_ratio, 2.0),
        1.0,
    );
    return snow_cover_fraction * snow_surface_albedo +
        (1.0 - snow_cover_fraction) * inputs.dry_soil_albedo;
}

fn soilAlbedo(inputs: Inputs) f64 {
    const surface_water_content_m3_m3 =
        if (inputs.soil_total_volume_m3 > inputs.negligible_soil_volume_m3)
            @min(
                inputs.soil_porosity_m3_m3,
                inputs.soil_water_volume_m3 / inputs.soil_total_volume_m3,
            )
        else
            0.0;
    return @min(
        inputs.maximum_soil_albedo,
        inputs.dry_soil_albedo +
            @max(0.0, inputs.maximum_soil_albedo -
                surface_water_content_m3_m3),
    );
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidGroundAlbedoInput;
    }
    if (inputs.cell_area_m2 == 0 or inputs.soil_porosity_m3_m3 > 1 or
        inputs.dry_soil_albedo > 1 or inputs.maximum_soil_albedo > 1)
        return error.InvalidGroundAlbedoInput;
}

test "snow surface albedo weights snow ice and water before cover mixing" {
    const result = try compute(.{
        .snowpack_heat_capacity_megajoules_k = 2,
        .minimum_snowpack_heat_capacity_megajoules_k = 1,
        .snow_volume_m3 = 1,
        .ice_volume_m3 = 1,
        .water_volume_m3 = 0,
        .snow_depth_m = 0.07,
        .soil_total_volume_m3 = 1,
        .soil_water_volume_m3 = 0,
        .soil_porosity_m3_m3 = 0.5,
        .negligible_soil_volume_m3 = 1e-12,
        .dry_soil_albedo = 0.1,
        .maximum_soil_albedo = 0.7,
        .ground_shortwave_megajoules_per_m2_h = 2,
        .ground_par_umol_per_m2_s = 100,
        .cell_area_m2 = 10,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), result.ground_albedo, 1e-15);
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.3),
        result.backscattered_shortwave_megajoules_per_m2_h_per_sector,
        1e-15,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 8),
        result.absorbed_ground_shortwave_megajoules_h,
        1e-14,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 400),
        result.absorbed_ground_par_umol_s,
        1e-12,
    );
}

test "soil branch derives albedo from bounded surface water content" {
    const result = try compute(.{
        .snowpack_heat_capacity_megajoules_k = 0,
        .minimum_snowpack_heat_capacity_megajoules_k = 1,
        .snow_volume_m3 = 0,
        .ice_volume_m3 = 0,
        .water_volume_m3 = 0,
        .snow_depth_m = 0,
        .soil_total_volume_m3 = 1,
        .soil_water_volume_m3 = 0.2,
        .soil_porosity_m3_m3 = 0.5,
        .negligible_soil_volume_m3 = 1e-12,
        .dry_soil_albedo = 0.1,
        .maximum_soil_albedo = 0.7,
        .ground_shortwave_megajoules_per_m2_h = 0,
        .ground_par_umol_per_m2_s = 0,
        .cell_area_m2 = 1,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), result.ground_albedo, 1e-15);
}
