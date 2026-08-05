const std = @import("std");

pub const Parameters = struct {
    minimum_richardson_number: f64,
    maximum_richardson_number: f64,
    richardson_resistance_multiplier: f64,
    minimum_aerodynamic_resistance_h_per_m: f64,
    maximum_aerodynamic_resistance_h_per_m: f64,
    canopy_drag_length_m: f64,
    minimum_air_fraction: f64,
};

pub const Inputs = struct {
    cell_area_m2: f64,
    air_temperature_k: f64,
    ground_air_temperature_k: f64,
    surface_temperature_k: f64,
    bulk_richardson_coefficient_k: f64,
    isothermal_atmospheric_resistance_h_per_m: f64,
    total_canopy_area_m2: f64,
    canopy_height_m: f64,
    atmospheric_vapor_diffusivity_m2_per_h: f64,
    isothermal_ground_surface_resistance_h_per_m: f64,
    bare_surface_fraction: f64,
    litter_surface_fraction: f64,
    litter_porous_resistance_h_per_m: f64,
    snow_layer_thickness_m: []const f64,
    snow_layer_total_volume_m3: []const f64,
    snow_layer_air_volume_m3: []const f64,
    snow_layer_vapor_diffusivity_m2_per_h: []const f64,
};

pub const Result = struct {
    atmosphere_to_ground_resistance_h_per_m: f64,
    ground_surface_resistance_h_per_m: f64,
    snowpack_resistance_h_per_m: f64,
    atmospheric_litter_gas_conductance_m3_per_h: f64,
    atmospheric_gas_conductance_m3_per_h: f64,
};

pub const SurfaceExchangeInputs = struct {
    cell_area_m2: f64,
    flux_timestep_h: f64,
    snow_flux_timestep_h: f64,
    litter_flux_timestep_h: f64,
    snow_depth_m: f64,
    full_snow_cover_depth_m: f64,
    nominal_bare_soil_fraction: f64,
    surface_excess_water_m3: f64,
    surface_water_capacity_m3: f64,
    current_ground_surface_resistance_h_per_m: f64,
    current_snow_surface_resistance_h_per_m: f64,
    litter_porous_resistance_h_per_m: f64,
    litter_air_fraction: f64,
    litter_porosity_m3_per_m3: f64,
    litter_tortuosity: f64,
    minimum_air_transport_factor: f64,
    soil_evaporation_pore_resistance_h_per_m: f64,
    soil_surface_air_fraction: f64,
    evaporation_surface_resistance_h_per_m: f64,
};

pub const SurfaceExchangeResult = struct {
    snow_cover_fraction: f64,
    snow_free_fraction: f64,
    water_adjusted_bare_soil_fraction: f64,
    water_adjusted_litter_fraction: f64,
    effective_ground_surface_resistance_h_per_m: f64,
    litter_evaporation_resistance_h_per_m: f64,
    soil_evaporation_resistance_h_per_m: f64,
    snow_latent_conductance_m3_per_step: f64,
    snow_sensible_conductance_megajoules_per_k_step: f64,
    litter_latent_conductance_m3_per_step: f64,
    litter_sensible_conductance_megajoules_per_k_step: f64,
    soil_latent_conductance_m3_per_step: f64,
    soil_sensible_conductance_megajoules_per_k_step: f64,
};

/// WATSUB surface fractions and `PAREWM/PARSWM`, `PARERM/PARSRM`, and
/// `PAREGM/PARSGM`. Conductances are whole-cell values for one nonlinear
/// flux step; callers divide by cell area only when an areal ledger is used.
pub fn calculateSurfaceExchange(inputs: SurfaceExchangeInputs) !SurfaceExchangeResult {
    inline for (@typeInfo(SurfaceExchangeInputs).@"struct".fields) |field| if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSurfaceExchangeInput;
    if (inputs.cell_area_m2 <= 0 or inputs.flux_timestep_h <= 0 or inputs.snow_flux_timestep_h <= 0 or inputs.litter_flux_timestep_h <= 0 or inputs.snow_depth_m < 0 or inputs.full_snow_cover_depth_m <= 0 or inputs.nominal_bare_soil_fraction < 0 or inputs.nominal_bare_soil_fraction > 1 or inputs.surface_excess_water_m3 < 0 or inputs.surface_water_capacity_m3 < 0 or inputs.current_ground_surface_resistance_h_per_m <= 0 or inputs.current_snow_surface_resistance_h_per_m <= 0 or inputs.litter_porous_resistance_h_per_m < 0 or inputs.litter_air_fraction < 0 or inputs.litter_porosity_m3_per_m3 <= 0 or inputs.litter_tortuosity <= 0 or inputs.minimum_air_transport_factor <= 0 or inputs.soil_evaporation_pore_resistance_h_per_m < 0 or inputs.soil_surface_air_fraction < 0 or inputs.evaporation_surface_resistance_h_per_m < 0) return error.InvalidSurfaceExchangeInput;

    var snow_cover = @min(1, @sqrt(inputs.snow_depth_m / inputs.full_snow_cover_depth_m));
    var snow_free: f64 = 0;
    if (snow_cover < 1) {
        snow_free = @max(1e-3, 1 - snow_cover);
        snow_cover = 1 - snow_free;
    }
    const water_capacity_fraction: f64 = if (inputs.surface_water_capacity_m3 > 0) std.math.clamp(inputs.surface_excess_water_m3 / inputs.surface_water_capacity_m3, 0, 1) else if (inputs.surface_excess_water_m3 > 0) 1.0 else 0.0;
    const bare = @max(0, inputs.nominal_bare_soil_fraction - water_capacity_fraction);
    const litter = 1 - bare;
    const air_transport_factor = @max(inputs.minimum_air_transport_factor, inputs.litter_tortuosity * inputs.litter_air_fraction * inputs.litter_air_fraction / inputs.litter_porosity_m3_per_m3);
    const porous_resistance = inputs.litter_porous_resistance_h_per_m / air_transport_factor;
    const litter_evaporation_resistance = inputs.litter_porous_resistance_h_per_m * inputs.litter_air_fraction;
    const soil_evaporation_resistance = inputs.soil_evaporation_pore_resistance_h_per_m * inputs.soil_surface_air_fraction;
    const combined_ground_resistance = 1 / (bare / inputs.current_ground_surface_resistance_h_per_m + litter / (inputs.current_ground_surface_resistance_h_per_m + porous_resistance));
    const latent_area_time_m2_h = inputs.cell_area_m2 * inputs.flux_timestep_h;
    const sensible_area_time_megajoules_m_per_k = inputs.cell_area_m2 * 1.25e-3 * inputs.flux_timestep_h;
    const result: SurfaceExchangeResult = .{
        .snow_cover_fraction = snow_cover,
        .snow_free_fraction = snow_free,
        .water_adjusted_bare_soil_fraction = bare,
        .water_adjusted_litter_fraction = litter,
        .effective_ground_surface_resistance_h_per_m = combined_ground_resistance,
        .litter_evaporation_resistance_h_per_m = litter_evaporation_resistance,
        .soil_evaporation_resistance_h_per_m = soil_evaporation_resistance,
        .snow_latent_conductance_m3_per_step = latent_area_time_m2_h * snow_cover * inputs.snow_flux_timestep_h / inputs.flux_timestep_h / (inputs.current_snow_surface_resistance_h_per_m + inputs.evaporation_surface_resistance_h_per_m),
        .snow_sensible_conductance_megajoules_per_k_step = sensible_area_time_megajoules_m_per_k * snow_cover * inputs.snow_flux_timestep_h / inputs.flux_timestep_h / inputs.current_snow_surface_resistance_h_per_m,
        .litter_latent_conductance_m3_per_step = latent_area_time_m2_h * snow_free * litter * inputs.litter_flux_timestep_h / inputs.flux_timestep_h / (inputs.current_ground_surface_resistance_h_per_m + inputs.evaporation_surface_resistance_h_per_m + 0.5 * litter_evaporation_resistance),
        .litter_sensible_conductance_megajoules_per_k_step = sensible_area_time_megajoules_m_per_k * snow_free * litter * inputs.litter_flux_timestep_h / inputs.flux_timestep_h / inputs.current_ground_surface_resistance_h_per_m,
        .soil_latent_conductance_m3_per_step = latent_area_time_m2_h * snow_free * bare / (combined_ground_resistance + inputs.evaporation_surface_resistance_h_per_m + soil_evaporation_resistance),
        .soil_sensible_conductance_megajoules_per_k_step = sensible_area_time_megajoules_m_per_k * snow_free * bare / combined_ground_resistance,
    };
    inline for (@typeInfo(SurfaceExchangeResult).@"struct".fields) |field| if (!std.math.isFinite(@field(result, field.name)) or @field(result, field.name) < 0) return error.InvalidSurfaceExchangeResult;
    return result;
}

/// Exact WATSUB `PARG=PAREX/(RATG+RAGS+RAS)` resistance chain for one
/// runtime cell. No constant boundary conductance is introduced.
pub fn calculate(inputs: Inputs, parameters: Parameters) !Result {
    inline for (.{ inputs.cell_area_m2, inputs.air_temperature_k, inputs.ground_air_temperature_k, inputs.surface_temperature_k, inputs.bulk_richardson_coefficient_k, inputs.isothermal_atmospheric_resistance_h_per_m, inputs.total_canopy_area_m2, inputs.canopy_height_m, inputs.atmospheric_vapor_diffusivity_m2_per_h, inputs.isothermal_ground_surface_resistance_h_per_m, inputs.bare_surface_fraction, inputs.litter_surface_fraction, inputs.litter_porous_resistance_h_per_m, parameters.minimum_richardson_number, parameters.maximum_richardson_number, parameters.richardson_resistance_multiplier, parameters.minimum_aerodynamic_resistance_h_per_m, parameters.maximum_aerodynamic_resistance_h_per_m, parameters.canopy_drag_length_m, parameters.minimum_air_fraction }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSurfaceGasResistanceInput;
    if (inputs.cell_area_m2 <= 0 or inputs.air_temperature_k <= 0 or inputs.ground_air_temperature_k <= 0 or inputs.atmospheric_vapor_diffusivity_m2_per_h <= 0 or inputs.isothermal_atmospheric_resistance_h_per_m < 0 or inputs.total_canopy_area_m2 < 0 or inputs.canopy_height_m < 0 or inputs.isothermal_ground_surface_resistance_h_per_m <= 0 or inputs.litter_porous_resistance_h_per_m < 0 or inputs.bare_surface_fraction < 0 or inputs.litter_surface_fraction < 0 or @abs(inputs.bare_surface_fraction + inputs.litter_surface_fraction - 1) > 1e-10 or parameters.maximum_richardson_number < parameters.minimum_richardson_number or parameters.richardson_resistance_multiplier <= 0 or parameters.minimum_aerodynamic_resistance_h_per_m <= 0 or parameters.maximum_aerodynamic_resistance_h_per_m < parameters.minimum_aerodynamic_resistance_h_per_m or parameters.canopy_drag_length_m < 0 or parameters.minimum_air_fraction < 0 or parameters.minimum_air_fraction > 1) return error.InvalidSurfaceGasResistanceInput;
    const layer_count = inputs.snow_layer_thickness_m.len;
    if (inputs.snow_layer_total_volume_m3.len != layer_count or inputs.snow_layer_air_volume_m3.len != layer_count or inputs.snow_layer_vapor_diffusivity_m2_per_h.len != layer_count) return error.SurfaceGasSnowDimensionMismatch;

    const atmosphere_temperature_difference_k = inputs.air_temperature_k - inputs.ground_air_temperature_k;
    const atmosphere_richardson = std.math.clamp(inputs.bulk_richardson_coefficient_k / inputs.air_temperature_k * atmosphere_temperature_difference_k, parameters.minimum_richardson_number, parameters.maximum_richardson_number);
    const atmosphere_stability = 1 - parameters.richardson_resistance_multiplier * atmosphere_richardson;
    if (atmosphere_stability <= 0) return error.InvalidAtmosphericStabilityCorrection;
    const canopy_area_density_m2_per_m3 = if (inputs.canopy_height_m > 0) inputs.total_canopy_area_m2 / (inputs.canopy_height_m * inputs.cell_area_m2) else 0;
    const canopy_resistance_h_per_m = if (inputs.canopy_height_m > 0) inputs.canopy_height_m * canopy_area_density_m2_per_m3 * parameters.canopy_drag_length_m / inputs.atmospheric_vapor_diffusivity_m2_per_h else 0;
    const atmosphere_to_ground = @min(parameters.maximum_aerodynamic_resistance_h_per_m, @max(parameters.minimum_aerodynamic_resistance_h_per_m, inputs.isothermal_atmospheric_resistance_h_per_m + canopy_resistance_h_per_m) / atmosphere_stability);

    const ground_temperature_difference_k = inputs.ground_air_temperature_k - inputs.surface_temperature_k;
    const ground_richardson = std.math.clamp(inputs.bulk_richardson_coefficient_k / inputs.ground_air_temperature_k * ground_temperature_difference_k, parameters.minimum_richardson_number, parameters.maximum_richardson_number);
    const ground_stability = 1 - parameters.richardson_resistance_multiplier * ground_richardson;
    if (ground_stability <= 0) return error.InvalidGroundStabilityCorrection;
    const current_surface_resistance = @min(parameters.maximum_aerodynamic_resistance_h_per_m, @max(parameters.minimum_aerodynamic_resistance_h_per_m, inputs.isothermal_ground_surface_resistance_h_per_m / ground_stability));
    const ground_conductance_m_per_h = inputs.bare_surface_fraction / current_surface_resistance + inputs.litter_surface_fraction / (current_surface_resistance + inputs.litter_porous_resistance_h_per_m);
    if (ground_conductance_m_per_h <= 0) return error.ZeroGroundSurfaceConductance;
    const ground_surface_resistance = 1 / ground_conductance_m_per_h;

    var snowpack_resistance: f64 = 0;
    for (inputs.snow_layer_thickness_m, inputs.snow_layer_total_volume_m3, inputs.snow_layer_air_volume_m3, inputs.snow_layer_vapor_diffusivity_m2_per_h) |thickness_m, total_volume_m3, air_volume_m3, diffusivity_m2_per_h| {
        inline for (.{ thickness_m, total_volume_m3, air_volume_m3, diffusivity_m2_per_h }) |value| if (!std.math.isFinite(value)) return error.NonFiniteSurfaceGasSnowState;
        if (thickness_m < 0 or total_volume_m3 < 0 or air_volume_m3 < 0 or air_volume_m3 > total_volume_m3 or diffusivity_m2_per_h <= 0) return error.InvalidSurfaceGasSnowState;
        if (total_volume_m3 <= 0 or thickness_m <= 0) continue;
        const air_fraction = @max(parameters.minimum_air_fraction, air_volume_m3 / total_volume_m3);
        if (air_fraction <= 0) return error.ZeroSnowAirFraction;
        snowpack_resistance += thickness_m / diffusivity_m2_per_h / (air_fraction * air_fraction);
    }
    const total_resistance = atmosphere_to_ground + ground_surface_resistance + snowpack_resistance;
    const result: Result = .{
        .atmosphere_to_ground_resistance_h_per_m = atmosphere_to_ground,
        .ground_surface_resistance_h_per_m = ground_surface_resistance,
        .snowpack_resistance_h_per_m = snowpack_resistance,
        .atmospheric_litter_gas_conductance_m3_per_h = inputs.cell_area_m2 / (atmosphere_to_ground + current_surface_resistance + snowpack_resistance),
        .atmospheric_gas_conductance_m3_per_h = inputs.cell_area_m2 / total_resistance,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| if (!std.math.isFinite(@field(result, field.name)) or @field(result, field.name) < 0) return error.InvalidSurfaceGasResistanceResult;
    return result;
}

test "PARG is area divided by exact three-part resistance chain" {
    const result = try calculate(.{ .cell_area_m2 = 10, .air_temperature_k = 280, .ground_air_temperature_k = 279, .surface_temperature_k = 278, .bulk_richardson_coefficient_k = 0, .isothermal_atmospheric_resistance_h_per_m = 2, .total_canopy_area_m2 = 20, .canopy_height_m = 2, .atmospheric_vapor_diffusivity_m2_per_h = 0.1, .isothermal_ground_surface_resistance_h_per_m = 1, .bare_surface_fraction = 0.25, .litter_surface_fraction = 0.75, .litter_porous_resistance_h_per_m = 3, .snow_layer_thickness_m = &.{ 0.05, 0.10 }, .snow_layer_total_volume_m3 = &.{ 1, 1 }, .snow_layer_air_volume_m3 = &.{ 0.5, 0.25 }, .snow_layer_vapor_diffusivity_m2_per_h = &.{ 0.1, 0.1 } }, .{ .minimum_richardson_number = -0.1, .maximum_richardson_number = 0.05, .richardson_resistance_multiplier = 10, .minimum_aerodynamic_resistance_h_per_m = 0.01, .maximum_aerodynamic_resistance_h_per_m = 100, .canopy_drag_length_m = 2e-4, .minimum_air_fraction = 1e-6 });
    const ratg = 2 + 2 * (20.0 / (2 * 10)) * 2e-4 / 0.1;
    const rags = 1.0 / (0.25 / 1.0 + 0.75 / 4.0);
    const ras = 0.05 / 0.1 / 0.25 + 0.10 / 0.1 / 0.0625;
    try std.testing.expectApproxEqAbs(ratg, result.atmosphere_to_ground_resistance_h_per_m, 1e-12);
    try std.testing.expectApproxEqAbs(rags, result.ground_surface_resistance_h_per_m, 1e-12);
    try std.testing.expectApproxEqAbs(ras, result.snowpack_resistance_h_per_m, 1e-12);
    try std.testing.expectApproxEqAbs(10.0 / (ratg + 1.0 + ras), result.atmospheric_litter_gas_conductance_m3_per_h, 1e-12);
    try std.testing.expectApproxEqAbs(10.0 / (ratg + rags + ras), result.atmospheric_gas_conductance_m3_per_h, 1e-12);
}

test "surface heat and vapor conductances reproduce WATSUB equations" {
    const value = try calculateSurfaceExchange(.{ .cell_area_m2 = 10, .flux_timestep_h = 0.25, .snow_flux_timestep_h = 0.5, .litter_flux_timestep_h = 0.125, .snow_depth_m = 0.01, .full_snow_cover_depth_m = 0.04, .nominal_bare_soil_fraction = 0.7, .surface_excess_water_m3 = 0.2, .surface_water_capacity_m3 = 1, .current_ground_surface_resistance_h_per_m = 0.01, .current_snow_surface_resistance_h_per_m = 0.02, .litter_porous_resistance_h_per_m = 0.03, .litter_air_fraction = 0.4, .litter_porosity_m3_per_m3 = 0.8, .litter_tortuosity = 0.5, .minimum_air_transport_factor = 1e-6, .soil_evaporation_pore_resistance_h_per_m = 0.04, .soil_surface_air_fraction = 0.25, .evaporation_surface_resistance_h_per_m = 0.0278 });
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), value.snow_cover_fraction, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), value.water_adjusted_bare_soil_fraction, 1e-15);
    const porous = 0.03 / (0.5 * 0.4 * 0.4 / 0.8);
    const ground = 1.0 / (0.5 / 0.01 + 0.5 / (0.01 + porous));
    try std.testing.expectApproxEqAbs(ground, value.effective_ground_surface_resistance_h_per_m, 1e-15);
    try std.testing.expectApproxEqAbs(10.0 * 0.5 * 0.5 / (0.02 + 0.0278), value.snow_latent_conductance_m3_per_step, 1e-12);
    try std.testing.expectApproxEqAbs(10.0 * 0.125 * 0.5 * 0.5 / (0.01 + 0.0278 + 0.5 * 0.03 * 0.4), value.litter_latent_conductance_m3_per_step, 1e-12);
}
