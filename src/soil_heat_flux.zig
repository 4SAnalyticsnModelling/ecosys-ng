const std = @import("std");

pub const TurbulenceParameters = struct {
    water_fraction_threshold: f64,
    air_fraction_threshold: f64,
    water_rayleigh_coefficient: f64,
    air_rayleigh_coefficient: f64,
    water_nusselt_denominator: f64,
    air_nusselt_denominator: f64,
    maximum_rayleigh_number: f64 = 1.0e4,
};

pub const CellConductivityInputs = struct {
    bulk_density_megagrams_per_m3: f64,
    liquid_water_fraction: f64,
    ice_fraction: f64,
    air_fraction: f64,
    fraction_of_pore_volume_air_filled: f64,
    solid_conductivity_numerator_m_megajoules_per_h_k: f64,
    solid_conductivity_denominator: f64,
    temperature_difference_k: f64,
};

/// Exact WATSUB TCND calculation, including Rayleigh/Nusselt enhancement of
/// liquid-water and air conductivity across the current face.
pub fn calculateCellConductivity(inputs: CellConductivityInputs, parameters: TurbulenceParameters) !f64 {
    inline for (@typeInfo(CellConductivityInputs).@"struct".fields) |field| if (!std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSoilHeatInput;
    inline for (@typeInfo(TurbulenceParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters, field.name))) return error.NonFiniteSoilHeatParameter;
    if (inputs.bulk_density_megagrams_per_m3 < 0 or inputs.liquid_water_fraction < 0 or inputs.ice_fraction < 0 or inputs.air_fraction < 0 or inputs.fraction_of_pore_volume_air_filled < 0 or inputs.solid_conductivity_numerator_m_megajoules_per_h_k < 0 or inputs.solid_conductivity_denominator < 0 or parameters.water_fraction_threshold < 0 or parameters.air_fraction_threshold < 0 or parameters.water_rayleigh_coefficient < 0 or parameters.air_rayleigh_coefficient < 0 or parameters.water_nusselt_denominator <= 0 or parameters.air_nusselt_denominator <= 0 or parameters.maximum_rayleigh_number <= 0) return error.InvalidSoilHeatInput;
    if (inputs.bulk_density_megagrams_per_m3 == 0 and inputs.liquid_water_fraction + inputs.ice_fraction == 0) return 0;
    const scaled_temperature_difference = @abs(inputs.temperature_difference_k) * 1.0e-6;
    const water_turbulent_fraction = std.math.pow(f64, @max(0.0, inputs.liquid_water_fraction - parameters.water_fraction_threshold), 3);
    const air_turbulent_fraction = std.math.pow(f64, @max(0.0, inputs.air_fraction - parameters.air_fraction_threshold), 3);
    const water_rayleigh = @min(parameters.maximum_rayleigh_number, parameters.water_rayleigh_coefficient * scaled_temperature_difference * water_turbulent_fraction);
    const air_rayleigh = @min(parameters.maximum_rayleigh_number, parameters.air_rayleigh_coefficient * scaled_temperature_difference * air_turbulent_fraction);
    const water_nusselt = @max(1.0, 0.68 + 0.67 * std.math.pow(f64, water_rayleigh, 0.25) / parameters.water_nusselt_denominator);
    const air_nusselt = @max(1.0, 0.68 + 0.67 * std.math.pow(f64, air_rayleigh, 0.25) / parameters.air_nusselt_denominator);
    const water_conductivity = 2.067e-3 * water_nusselt;
    const air_conductivity = 9.050e-5 * air_nusselt;
    const air_weight = 1.467 - 0.467 * inputs.fraction_of_pore_volume_air_filled;
    const numerator = inputs.solid_conductivity_numerator_m_megajoules_per_h_k + inputs.liquid_water_fraction * water_conductivity + 0.611 * inputs.ice_fraction * 7.844e-3 + air_weight * inputs.air_fraction * air_conductivity;
    const denominator = inputs.solid_conductivity_denominator + inputs.liquid_water_fraction + 0.611 * inputs.ice_fraction + air_weight * inputs.air_fraction;
    if (!std.math.isFinite(numerator) or !std.math.isFinite(denominator) or denominator <= 0) return error.InvalidSoilHeatConductivity;
    return numerator / denominator;
}

pub const FaceInputs = struct {
    source_temperature_k: f64,
    destination_temperature_k: f64,
    source_heat_capacity_megajoules_per_k: f64,
    destination_heat_capacity_megajoules_per_k: f64,
    source_minimum_heat_capacity_megajoules_per_k: f64,
    destination_minimum_heat_capacity_megajoules_per_k: f64,
    source_is_top_soil_layer: bool,
    top_snow_heat_capacity_megajoules_per_k: f64,
    maximum_negligible_snow_heat_capacity_megajoules_per_k: f64,
    snow_storage_heat_flux_megajoules: f64,
    liquid_water_flux_m3: f64,
    vapor_flux_m3: f64,
    macropore_water_flux_m3: f64,
    liquid_water_heat_capacity_megajoules_per_m3_k: f64,
    source_thermal_conductivity_m_megajoules_per_h_k: f64,
    destination_thermal_conductivity_m_megajoules_per_h_k: f64,
    source_path_length_m: f64,
    destination_path_length_m: f64,
    face_area_m2: f64,
    time_fraction: f64,
};

pub const FaceFlux = struct {
    conductive_unlimited_megajoules: f64,
    conductive_limited_megajoules: f64,
    convective_megajoules: f64,
    total_megajoules: f64,
    equilibrium_temperature_k: f64,
};

/// Exact HFLWL face calculation after liquid, vapor and macropore water fluxes
/// have been determined. Positive heat moves source to destination.
pub fn calculateFaceFlux(inputs: FaceInputs) !FaceFlux {
    inline for (@typeInfo(FaceInputs).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSoilHeatInput;
    if (inputs.source_temperature_k <= 0 or inputs.destination_temperature_k <= 0 or inputs.source_heat_capacity_megajoules_per_k <= 0 or inputs.destination_heat_capacity_megajoules_per_k <= 0 or inputs.source_minimum_heat_capacity_megajoules_per_k < 0 or inputs.destination_minimum_heat_capacity_megajoules_per_k < 0 or inputs.top_snow_heat_capacity_megajoules_per_k < 0 or inputs.maximum_negligible_snow_heat_capacity_megajoules_per_k < 0 or inputs.liquid_water_heat_capacity_megajoules_per_m3_k <= 0 or inputs.source_thermal_conductivity_m_megajoules_per_h_k < 0 or inputs.destination_thermal_conductivity_m_megajoules_per_h_k < 0 or inputs.source_path_length_m <= 0 or inputs.destination_path_length_m <= 0 or inputs.face_area_m2 < 0 or inputs.time_fraction <= 0 or inputs.time_fraction > 1) return error.InvalidSoilHeatInput;
    const vapor_donor_temperature = if (inputs.vapor_flux_m3 >= 0) inputs.source_temperature_k else inputs.destination_temperature_k;
    const liquid_donor_temperature = if (inputs.liquid_water_flux_m3 >= 0) inputs.source_temperature_k else inputs.destination_temperature_k;
    const macro_donor_temperature = if (inputs.macropore_water_flux_m3 >= 0) inputs.source_temperature_k else inputs.destination_temperature_k;
    const vapor_convective = inputs.liquid_water_heat_capacity_megajoules_per_m3_k * vapor_donor_temperature * inputs.vapor_flux_m3;
    const liquid_convective = inputs.liquid_water_heat_capacity_megajoules_per_m3_k * liquid_donor_temperature * inputs.liquid_water_flux_m3;
    const macro_convective = inputs.liquid_water_heat_capacity_megajoules_per_m3_k * macro_donor_temperature * inputs.macropore_water_flux_m3;
    var source_interim = inputs.source_temperature_k;
    if (inputs.source_heat_capacity_megajoules_per_k > inputs.source_minimum_heat_capacity_megajoules_per_k) {
        source_interim -= if (inputs.source_is_top_soil_layer and inputs.top_snow_heat_capacity_megajoules_per_k <= inputs.maximum_negligible_snow_heat_capacity_megajoules_per_k)
            (vapor_convective - inputs.snow_storage_heat_flux_megajoules) / inputs.source_heat_capacity_megajoules_per_k
        else
            vapor_convective / inputs.source_heat_capacity_megajoules_per_k;
    }
    var destination_interim = inputs.destination_temperature_k;
    if (inputs.destination_heat_capacity_megajoules_per_k > inputs.destination_minimum_heat_capacity_megajoules_per_k) destination_interim += vapor_convective / inputs.destination_heat_capacity_megajoules_per_k;
    const equilibrium = (inputs.source_heat_capacity_megajoules_per_k * source_interim + inputs.destination_heat_capacity_megajoules_per_k * destination_interim) / (inputs.source_heat_capacity_megajoules_per_k + inputs.destination_heat_capacity_megajoules_per_k);
    // `unlimited` is already integrated over `time_fraction`. The physical
    // pair-equilibration ceiling is the energy required to reach the shared
    // endpoint and must not be multiplied by the step fraction a second
    // time. At the ordinary whole-hour value of one this is exactly the
    // translated WATSUB expression; the correction matters for runtime
    // fine-step validation and other non-hourly kernels.
    const equilibration_limit =
        (source_interim - equilibrium) *
        inputs.source_heat_capacity_megajoules_per_k;
    const denominator = inputs.source_thermal_conductivity_m_megajoules_per_h_k * inputs.destination_path_length_m + inputs.destination_thermal_conductivity_m_megajoules_per_h_k * inputs.source_path_length_m;
    const conductance = if (denominator > 0) 2.0 * inputs.source_thermal_conductivity_m_megajoules_per_h_k * inputs.destination_thermal_conductivity_m_megajoules_per_h_k / denominator else 0.0;
    const unlimited = conductance * (source_interim - destination_interim) * inputs.face_area_m2 * inputs.time_fraction;
    const limited = if (unlimited >= 0) @max(0.0, @min(equilibration_limit, unlimited)) else @min(0.0, @max(equilibration_limit, unlimited));
    const convective = liquid_convective + vapor_convective + macro_convective;
    const total = convective + limited;
    if (!std.math.isFinite(total)) return error.NonFiniteSoilHeatFlux;
    return .{ .conductive_unlimited_megajoules = unlimited, .conductive_limited_megajoules = limited, .convective_megajoules = convective, .total_megajoules = total, .equilibrium_temperature_k = equilibrium };
}

pub fn waterFilmThicknessM(matric_potential_mpa: f64) !f64 {
    if (!std.math.isFinite(matric_potential_mpa) or matric_potential_mpa >= 0) return error.InvalidMatricPotentialForWaterFilm;
    return @max(1.0e-6, 0.5 * @exp(-13.833 - 0.857 * @log(-matric_potential_mpa)));
}

test "WATSUB conductive face is limited to pair equilibration" {
    const flux = try calculateFaceFlux(.{ .source_temperature_k = 300, .destination_temperature_k = 280, .source_heat_capacity_megajoules_per_k = 2, .destination_heat_capacity_megajoules_per_k = 2, .source_minimum_heat_capacity_megajoules_per_k = 0, .destination_minimum_heat_capacity_megajoules_per_k = 0, .source_is_top_soil_layer = false, .top_snow_heat_capacity_megajoules_per_k = 0, .maximum_negligible_snow_heat_capacity_megajoules_per_k = 0, .snow_storage_heat_flux_megajoules = 0, .liquid_water_flux_m3 = 0, .vapor_flux_m3 = 0, .macropore_water_flux_m3 = 0, .liquid_water_heat_capacity_megajoules_per_m3_k = 4.19, .source_thermal_conductivity_m_megajoules_per_h_k = 100, .destination_thermal_conductivity_m_megajoules_per_h_k = 100, .source_path_length_m = 1, .destination_path_length_m = 1, .face_area_m2 = 1, .time_fraction = 1 });
    try std.testing.expectApproxEqAbs(@as(f64, 290), flux.equilibrium_temperature_k, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 20), flux.conductive_limited_megajoules, 1e-12);
}

test "convective heat uses the upstream temperature for each phase" {
    const flux = try calculateFaceFlux(.{ .source_temperature_k = 300, .destination_temperature_k = 280, .source_heat_capacity_megajoules_per_k = 2, .destination_heat_capacity_megajoules_per_k = 2, .source_minimum_heat_capacity_megajoules_per_k = 0, .destination_minimum_heat_capacity_megajoules_per_k = 0, .source_is_top_soil_layer = false, .top_snow_heat_capacity_megajoules_per_k = 0, .maximum_negligible_snow_heat_capacity_megajoules_per_k = 0, .snow_storage_heat_flux_megajoules = 0, .liquid_water_flux_m3 = 0.01, .vapor_flux_m3 = -0.02, .macropore_water_flux_m3 = 0.03, .liquid_water_heat_capacity_megajoules_per_m3_k = 4.2, .source_thermal_conductivity_m_megajoules_per_h_k = 0, .destination_thermal_conductivity_m_megajoules_per_h_k = 0, .source_path_length_m = 1, .destination_path_length_m = 1, .face_area_m2 = 1, .time_fraction = 1 });
    const expected = 4.2 * (300 * 0.01 + 280 * -0.02 + 300 * 0.03);
    try std.testing.expectApproxEqAbs(expected, flux.convective_megajoules, 1e-12);
}

test "water film follows WATSUB lower bound" {
    try std.testing.expect(try waterFilmThicknessM(-0.01) >= 1e-6);
    try std.testing.expectEqual(@as(f64, 1e-6), try waterFilmThicknessM(-1e12));
}
