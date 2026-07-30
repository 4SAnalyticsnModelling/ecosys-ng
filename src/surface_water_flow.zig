const std = @import("std");
const water_flux = @import("soil_water_flux.zig");

pub const LitterSoilInputs = struct {
    litter_water_m3: f64,
    soil_matrix_water_m3: f64,
    litter_air_m3: f64,
    soil_matrix_air_m3: f64,
    litter_volume_m3: f64,
    soil_matrix_bulk_volume_m3: f64,
    litter_water_fraction: f64,
    soil_water_fraction: f64,
    litter_air_entry_water_fraction: f64,
    soil_air_entry_water_fraction: f64,
    litter_total_water_potential_mpa: f64,
    soil_total_water_potential_mpa: f64,
    litter_hydraulic_conductivity_m2_per_h_mpa: f64,
    soil_hydraulic_conductivity_m2_per_h_mpa: f64,
    litter_thickness_m: f64,
    soil_thickness_m: f64,
    soil_face_area_m2: f64,
    litter_cover_fraction: f64,
    wet_litter_cover_fraction: f64,
    time_fraction: f64,
    soil_excess_pore_volume_m3: f64,
    litter_temperature_k: f64,
    soil_temperature_k: f64,
};

pub const LitterSoilFlux = struct { water_m3: f64, unenhanced_water_m3: f64, convective_heat_mj: f64 };

/// NPR is now only the caller's convergence ceiling; this evaluates one
/// nonlinear litter-soil face for a whole-step Newton/Picard residual.
pub fn litterSoilFlux(inputs: LitterSoilInputs) !LitterSoilFlux {
    inline for (@typeInfo(LitterSoilInputs).@"struct".fields) |field| if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSurfaceWaterInput;
    if (inputs.litter_cover_fraction < 0 or inputs.litter_cover_fraction > 1 or inputs.wet_litter_cover_fraction < 0 or inputs.wet_litter_cover_fraction > 1 or inputs.litter_temperature_k <= 0 or inputs.soil_temperature_k <= 0) return error.InvalidSurfaceWaterInput;
    // Scaling both endpoint conductivities scales the harmonic conductance by
    // CVRD exactly; CVRDW scales the face area in FLQX.
    const flux = try water_flux.calculateMatrixFaceFlux(.{ .direction = .vertical, .source_water_m3 = inputs.litter_water_m3, .destination_water_m3 = inputs.soil_matrix_water_m3, .source_air_m3 = inputs.litter_air_m3, .destination_air_m3 = inputs.soil_matrix_air_m3, .source_micropore_volume_m3 = inputs.litter_volume_m3, .destination_micropore_volume_m3 = inputs.soil_matrix_bulk_volume_m3, .source_water_fraction = inputs.litter_water_fraction, .destination_water_fraction = inputs.soil_water_fraction, .source_air_entry_water_fraction = inputs.litter_air_entry_water_fraction, .destination_air_entry_water_fraction = inputs.soil_air_entry_water_fraction, .source_total_water_potential_mpa = inputs.litter_total_water_potential_mpa, .destination_total_water_potential_mpa = inputs.soil_total_water_potential_mpa, .source_hydraulic_conductivity_m2_per_h_mpa = inputs.litter_hydraulic_conductivity_m2_per_h_mpa * inputs.litter_cover_fraction, .destination_hydraulic_conductivity_m2_per_h_mpa = inputs.soil_hydraulic_conductivity_m2_per_h_mpa * inputs.litter_cover_fraction, .source_path_length_m = inputs.litter_thickness_m, .destination_path_length_m = inputs.soil_thickness_m, .face_area_m2 = inputs.soil_face_area_m2 * inputs.wet_litter_cover_fraction, .wetting_depth_factor = inputs.time_fraction, .time_fraction = inputs.time_fraction, .destination_excess_pore_volume_m3 = inputs.soil_excess_pore_volume_m3 });
    const donor_temperature = if (flux.limited_water_m3 > 0) inputs.litter_temperature_k else inputs.soil_temperature_k;
    return .{ .water_m3 = flux.limited_water_m3, .unenhanced_water_m3 = flux.transport_water_m3, .convective_heat_mj = 4.19 * donor_temperature * flux.limited_water_m3 };
}

pub fn pondToSoilWaterM3(pond_water_m3: f64, surface_area_m2: f64, time_fraction: f64) !f64 {
    inline for (.{ pond_water_m3, surface_area_m2, time_fraction }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSurfaceWaterInput;
    if (pond_water_m3 < 0 or surface_area_m2 <= 0 or time_fraction <= 0 or time_fraction > 1) return error.InvalidSurfaceWaterInput;
    return @max(0.0, pond_water_m3 - 0.01 / surface_area_m2) * time_fraction;
}

pub fn litterOverflowToMacroporeM3(excess_litter_water_m3: f64, macropore_air_m3: f64, time_fraction: f64) !f64 {
    inline for (.{ excess_litter_water_m3, macropore_air_m3, time_fraction }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSurfaceWaterInput;
    if (excess_litter_water_m3 < 0 or macropore_air_m3 < 0 or time_fraction <= 0 or time_fraction > 1) return error.InvalidSurfaceWaterInput;
    return @min(excess_litter_water_m3 * time_fraction, macropore_air_m3);
}

pub const StoragePartition = struct {
    retained_liquid_water_m3: f64,
    retained_ice_m3: f64,
    excess_liquid_water_m3: f64,
    excess_ice_m3: f64,
    total_excess_water_and_ice_m3: f64,
};

pub fn partitionSurfaceStorage(liquid_water_m3: f64, ice_m3: f64, litter_retention_capacity_m3: f64) !StoragePartition {
    inline for (.{ liquid_water_m3, ice_m3, litter_retention_capacity_m3 }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSurfaceWaterInput;
    if (liquid_water_m3 < 0 or ice_m3 < 0 or litter_retention_capacity_m3 < 0) return error.InvalidSurfaceWaterInput;
    const total = liquid_water_m3 + ice_m3;
    const retained_liquid = if (total > 0) liquid_water_m3 / total * litter_retention_capacity_m3 else 0;
    const retained_ice = if (total > 0) ice_m3 / total * litter_retention_capacity_m3 else 0;
    return .{ .retained_liquid_water_m3 = retained_liquid, .retained_ice_m3 = retained_ice, .excess_liquid_water_m3 = @max(0.0, liquid_water_m3 - retained_liquid), .excess_ice_m3 = @max(0.0, ice_m3 - retained_ice), .total_excess_water_and_ice_m3 = @max(0.0, total - litter_retention_capacity_m3) };
}

pub const RunoffInputs = struct {
    soil_surface_present: bool,
    total_excess_water_and_ice_m3: f64,
    excess_liquid_water_m3: f64,
    ground_surface_retention_capacity_m3: f64,
    soil_surface_depth_m: f64,
    natural_water_table_depth_m: f64,
    surface_area_m2: f64,
    surface_slope: f64,
    roughness_height_m: f64,
    flow_width_m: f64,
};

pub const Runoff = struct { water_m3_per_step: f64, velocity_m_per_s: f64, available_ponded_water_m3: f64 };

/// Exact WATSUB Manning runoff calculation and its 1e-3 per-step availability
/// cap. Directional partitioning remains in surface_solute_routing.zig.
pub fn runoff(inputs: RunoffInputs) !Runoff {
    inline for (@typeInfo(RunoffInputs).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSurfaceWaterInput;
    if (inputs.total_excess_water_and_ice_m3 < 0 or inputs.excess_liquid_water_m3 < 0 or inputs.ground_surface_retention_capacity_m3 < 0 or inputs.surface_area_m2 <= 0 or inputs.surface_slope < 0 or inputs.roughness_height_m <= 0 or inputs.flow_width_m < 0) return error.InvalidSurfaceWaterInput;
    const available = if (inputs.soil_surface_present and inputs.total_excess_water_and_ice_m3 > inputs.ground_surface_retention_capacity_m3)
        @min(1.0e-3, (inputs.total_excess_water_and_ice_m3 - inputs.ground_surface_retention_capacity_m3) * if (inputs.total_excess_water_and_ice_m3 > 0) inputs.excess_liquid_water_m3 / inputs.total_excess_water_and_ice_m3 else 0)
    else if (!inputs.soil_surface_present and inputs.soil_surface_depth_m <= inputs.natural_water_table_depth_m)
        @min(1.0e-3, inputs.natural_water_table_depth_m - inputs.soil_surface_depth_m)
    else
        0;
    if (available <= 0) return .{ .water_m3_per_step = 0, .velocity_m_per_s = 0, .available_ponded_water_m3 = 0 };
    const hydraulic_radius_m = available / inputs.surface_area_m2;
    const velocity = std.math.pow(f64, hydraulic_radius_m, 0.67) * @sqrt(inputs.surface_slope) / inputs.roughness_height_m;
    const requested = velocity * hydraulic_radius_m * inputs.flow_width_m * 3.6e3;
    return .{ .water_m3_per_step = @min(requested, available), .velocity_m_per_s = velocity, .available_ponded_water_m3 = available };
}

pub fn surfaceWaterFilmThicknessM(matric_potential_mpa: f64, heat_capacity_active: bool) !f64 {
    if (!std.math.isFinite(matric_potential_mpa) or matric_potential_mpa >= 0) return error.InvalidSurfaceMatricPotential;
    if (!heat_capacity_active) return 1e-6;
    return @max(1e-6, 0.5 * @exp(-13.650 - 0.857 * @log(-matric_potential_mpa)));
}

test "surface storage partitions liquid and ice proportionally" {
    const result = try partitionSurfaceStorage(0.8, 0.2, 0.5);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), result.retained_liquid_water_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), result.retained_ice_m3, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), result.total_excess_water_and_ice_m3, 1e-12);
}

test "Manning runoff is capped by available ponded water" {
    const result = try runoff(.{ .soil_surface_present = true, .total_excess_water_and_ice_m3 = 1, .excess_liquid_water_m3 = 1, .ground_surface_retention_capacity_m3 = 0, .soil_surface_depth_m = 0, .natural_water_table_depth_m = 1, .surface_area_m2 = 1, .surface_slope = 1, .roughness_height_m = 0.01, .flow_width_m = 1 });
    try std.testing.expectEqual(@as(f64, 1e-3), result.available_ponded_water_m3);
    try std.testing.expect(result.water_m3_per_step <= 1e-3);
    try std.testing.expect(result.velocity_m_per_s > 0);
}

test "pond and litter overflow preserve source limits" {
    try std.testing.expectApproxEqAbs(@as(f64, 0.09), try pondToSoilWaterM3(0.1, 1, 1), 1e-12);
    try std.testing.expectEqual(@as(f64, 0.02), try litterOverflowToMacroporeM3(0.1, 0.02, 1));
}
