const std = @import("std");

pub const Zone = struct {
    dissolved_fraction: f64,
    concentration_g_per_m3: f64,
    dissolved_amount_g: f64,
    competition_fraction: f64,
};

pub const UptakeParameters = struct {
    minimum_concentration_g_per_m3: f64,
    half_saturation_g_per_m3: f64,
    maximum_uptake_g_per_m2_h: f64,
};

pub const CommonInputs = struct {
    microbial_nonstructural_c_g: f64,
    microbial_nonstructural_nutrient_g: f64,
    maximum_nutrient_per_c_g_per_g_c: f64,
    microbial_surface_area_m2_per_g_c: f64,
    active_biomass_g_c: f64,
    temperature_water_response: f64,
    water_volume_m3: f64,
    timestep_h: f64,
};

pub const TwoZoneFlux = struct {
    non_band_g: f64,
    band_g: f64,
    supply_unlimited_non_band_g: f64,
    supply_unlimited_band_g: f64,

    pub fn total(self: TwoZoneFlux) f64 {
        return self.non_band_g + self.band_g;
    }
};

pub const NitrogenResult = struct {
    ammonium: TwoZoneFlux,
    nitrate: TwoZoneFlux,
    initial_nitrogen_demand_g_n: f64,
    unmet_nitrogen_demand_g_n: f64,
};

/// Ports NITRO's signed NH4 mineralization/immobilization followed by NO3
/// immobilization of any remaining positive microbial N demand.
pub fn nitrogen(common: CommonInputs, ammonium_zones: [2]Zone, nitrate_zones: [2]Zone, ammonium_parameters: UptakeParameters, nitrate_parameters: UptakeParameters) !NitrogenResult {
    try validateCommon(common);
    try validateZones(ammonium_zones);
    try validateZones(nitrate_zones);
    try validateParameters(ammonium_parameters);
    try validateParameters(nitrate_parameters);
    const demand = @max(0, common.microbial_nonstructural_c_g) * common.maximum_nutrient_per_c_g_per_g_c - @max(0, common.microbial_nonstructural_nutrient_g);
    const ammonium = calculate(demand, common, ammonium_zones, ammonium_parameters);
    const remaining = @max(0, demand - ammonium.total());
    const nitrate = calculate(remaining, common, nitrate_zones, nitrate_parameters);
    return .{ .ammonium = ammonium, .nitrate = nitrate, .initial_nitrogen_demand_g_n = demand, .unmet_nitrogen_demand_g_n = @max(0, remaining - nitrate.total()) };
}

pub const PhosphorusResult = struct {
    dihydrogen_phosphate: TwoZoneFlux,
    hydrogen_phosphate: TwoZoneFlux,
    initial_phosphorus_demand_g_p: f64,
    unmet_phosphorus_demand_g_p: f64,
};

/// Ports the corresponding H2PO4 then HPO4 sequence. Negative H2PO4 flux is
/// mineralization; HPO4 is considered only for remaining immobilization.
pub fn phosphorus(common: CommonInputs, dihydrogen_phosphate_zones: [2]Zone, hydrogen_phosphate_zones: [2]Zone, dihydrogen_parameters: UptakeParameters, hydrogen_parameters: UptakeParameters) !PhosphorusResult {
    try validateCommon(common);
    try validateZones(dihydrogen_phosphate_zones);
    try validateZones(hydrogen_phosphate_zones);
    try validateParameters(dihydrogen_parameters);
    try validateParameters(hydrogen_parameters);
    const demand = @max(0, common.microbial_nonstructural_c_g) * common.maximum_nutrient_per_c_g_per_g_c - @max(0, common.microbial_nonstructural_nutrient_g);
    const dihydrogen = calculate(demand, common, dihydrogen_phosphate_zones, dihydrogen_parameters);
    const remaining = @max(0, demand - dihydrogen.total());
    const hydrogen = calculate(remaining, common, hydrogen_phosphate_zones, hydrogen_parameters);
    return .{ .dihydrogen_phosphate = dihydrogen, .hydrogen_phosphate = hydrogen, .initial_phosphorus_demand_g_p = demand, .unmet_phosphorus_demand_g_p = @max(0, remaining - hydrogen.total()) };
}

fn calculate(demand_g: f64, common: CommonInputs, zones: [2]Zone, parameters: UptakeParameters) TwoZoneFlux {
    if (demand_g <= 0) return .{ .non_band_g = demand_g * zones[0].dissolved_fraction, .band_g = demand_g * zones[1].dissolved_fraction, .supply_unlimited_non_band_g = 0, .supply_unlimited_band_g = 0 };
    const biological_limit = @min(demand_g, common.microbial_surface_area_m2_per_g_c * common.active_biomass_g_c * common.temperature_water_response * parameters.maximum_uptake_g_per_m2_h * common.timestep_h);
    var unlimited: [2]f64 = .{ 0, 0 };
    var actual: [2]f64 = .{ 0, 0 };
    for (zones, 0..) |zone, index| {
        const available_concentration = @max(0, zone.concentration_g_per_m3 - parameters.minimum_concentration_g_per_m3);
        unlimited[index] = zone.dissolved_fraction * biological_limit * available_concentration / (available_concentration + parameters.half_saturation_g_per_m3);
        const protected_amount_g = parameters.minimum_concentration_g_per_m3 * common.water_volume_m3 * zone.dissolved_fraction;
        const supply_g = zone.competition_fraction * @max(0, zone.dissolved_amount_g - protected_amount_g) * common.timestep_h;
        actual[index] = @min(supply_g, unlimited[index]);
    }
    return .{ .non_band_g = actual[0], .band_g = actual[1], .supply_unlimited_non_band_g = unlimited[0], .supply_unlimited_band_g = unlimited[1] };
}

pub const NitrogenState = struct { microbial_n_g: f64, non_band_ammonium_g_n: f64, band_ammonium_g_n: f64, non_band_nitrate_g_n: f64, band_nitrate_g_n: f64 };

pub fn commitNitrogen(state: *NitrogenState, result: NitrogenResult) !void {
    inline for (@typeInfo(NitrogenState).@"struct".fields) |field| if (!std.math.isFinite(@field(state.*, field.name)) or @field(state.*, field.name) < 0) return error.InvalidNitrogenExchangeState;
    const fluxes = .{ result.ammonium.non_band_g, result.ammonium.band_g, result.nitrate.non_band_g, result.nitrate.band_g };
    inline for (fluxes) |flux| if (!std.math.isFinite(flux)) return error.NonFiniteNitrogenExchangeFlux;
    if (state.non_band_ammonium_g_n < result.ammonium.non_band_g or state.band_ammonium_g_n < result.ammonium.band_g or state.non_band_nitrate_g_n < result.nitrate.non_band_g or state.band_nitrate_g_n < result.nitrate.band_g) return error.InsufficientMineralNitrogen;
    const microbial_after = state.microbial_n_g + result.ammonium.total() + result.nitrate.total();
    if (microbial_after < 0) return error.NegativeMicrobialNitrogen;
    state.microbial_n_g = microbial_after;
    state.non_band_ammonium_g_n -= result.ammonium.non_band_g;
    state.band_ammonium_g_n -= result.ammonium.band_g;
    state.non_band_nitrate_g_n -= result.nitrate.non_band_g;
    state.band_nitrate_g_n -= result.nitrate.band_g;
}

pub const PhosphorusState = struct { microbial_p_g: f64, non_band_dihydrogen_phosphate_g_p: f64, band_dihydrogen_phosphate_g_p: f64, non_band_hydrogen_phosphate_g_p: f64, band_hydrogen_phosphate_g_p: f64 };

pub fn commitPhosphorus(state: *PhosphorusState, result: PhosphorusResult) !void {
    inline for (@typeInfo(PhosphorusState).@"struct".fields) |field| if (!std.math.isFinite(@field(state.*, field.name)) or @field(state.*, field.name) < 0) return error.InvalidPhosphorusExchangeState;
    const fluxes = .{ result.dihydrogen_phosphate.non_band_g, result.dihydrogen_phosphate.band_g, result.hydrogen_phosphate.non_band_g, result.hydrogen_phosphate.band_g };
    inline for (fluxes) |flux| if (!std.math.isFinite(flux)) return error.NonFinitePhosphorusExchangeFlux;
    if (state.non_band_dihydrogen_phosphate_g_p < result.dihydrogen_phosphate.non_band_g or state.band_dihydrogen_phosphate_g_p < result.dihydrogen_phosphate.band_g or state.non_band_hydrogen_phosphate_g_p < result.hydrogen_phosphate.non_band_g or state.band_hydrogen_phosphate_g_p < result.hydrogen_phosphate.band_g) return error.InsufficientMineralPhosphorus;
    const microbial_after = state.microbial_p_g + result.dihydrogen_phosphate.total() + result.hydrogen_phosphate.total();
    if (microbial_after < 0) return error.NegativeMicrobialPhosphorus;
    state.microbial_p_g = microbial_after;
    state.non_band_dihydrogen_phosphate_g_p -= result.dihydrogen_phosphate.non_band_g;
    state.band_dihydrogen_phosphate_g_p -= result.dihydrogen_phosphate.band_g;
    state.non_band_hydrogen_phosphate_g_p -= result.hydrogen_phosphate.non_band_g;
    state.band_hydrogen_phosphate_g_p -= result.hydrogen_phosphate.band_g;
}

fn validateCommon(common: CommonInputs) !void {
    inline for (@typeInfo(CommonInputs).@"struct".fields) |field| if (!std.math.isFinite(@field(common, field.name))) return error.NonFiniteNutrientExchangeInput;
    if (common.microbial_nonstructural_c_g < 0 or common.microbial_nonstructural_nutrient_g < 0 or common.maximum_nutrient_per_c_g_per_g_c < 0 or common.microbial_surface_area_m2_per_g_c < 0 or common.active_biomass_g_c < 0 or common.temperature_water_response < 0 or common.water_volume_m3 < 0 or common.timestep_h <= 0) return error.InvalidNutrientExchangeInput;
}

fn validateZones(zones: [2]Zone) !void {
    var fraction_sum: f64 = 0;
    for (zones) |zone| {
        inline for (@typeInfo(Zone).@"struct".fields) |field| if (!std.math.isFinite(@field(zone, field.name)) or @field(zone, field.name) < 0) return error.InvalidNutrientExchangeZone;
        fraction_sum += zone.dissolved_fraction;
    }
    if (@abs(fraction_sum - 1) > 1e-12) return error.NutrientZoneFractionsDoNotSumToOne;
}

fn validateParameters(parameters: UptakeParameters) !void {
    inline for (@typeInfo(UptakeParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(parameters, field.name)) or @field(parameters, field.name) < 0) return error.InvalidNutrientUptakeParameter;
    if (parameters.half_saturation_g_per_m3 <= 0) return error.InvalidNutrientUptakeParameter;
}

fn testCommon(nutrient_g: f64) CommonInputs {
    return .{ .microbial_nonstructural_c_g = 10, .microbial_nonstructural_nutrient_g = nutrient_g, .maximum_nutrient_per_c_g_per_g_c = 0.1, .microbial_surface_area_m2_per_g_c = 1, .active_biomass_g_c = 1, .temperature_water_response = 1, .water_volume_m3 = 1, .timestep_h = 1 };
}

fn testZones() [2]Zone {
    return .{ .{ .dissolved_fraction = 0.7, .concentration_g_per_m3 = 2, .dissolved_amount_g = 2, .competition_fraction = 1 }, .{ .dissolved_fraction = 0.3, .concentration_g_per_m3 = 2, .dissolved_amount_g = 1, .competition_fraction = 1 } };
}

fn testParameters() UptakeParameters {
    return .{ .minimum_concentration_g_per_m3 = 0.01, .half_saturation_g_per_m3 = 0.1, .maximum_uptake_g_per_m2_h = 1 };
}

test "sequential ammonium nitrate immobilization conserves nitrogen" {
    const result = try nitrogen(testCommon(0.2), testZones(), testZones(), testParameters(), testParameters());
    var state: NitrogenState = .{ .microbial_n_g = 0.2, .non_band_ammonium_g_n = 2, .band_ammonium_g_n = 1, .non_band_nitrate_g_n = 2, .band_nitrate_g_n = 1 };
    const before = state.microbial_n_g + state.non_band_ammonium_g_n + state.band_ammonium_g_n + state.non_band_nitrate_g_n + state.band_nitrate_g_n;
    try commitNitrogen(&state, result);
    try std.testing.expectApproxEqAbs(before, state.microbial_n_g + state.non_band_ammonium_g_n + state.band_ammonium_g_n + state.non_band_nitrate_g_n + state.band_nitrate_g_n, 1e-14);
    try std.testing.expect(result.ammonium.total() > 0);
}

test "excess microbial nitrogen mineralizes ammonium by zone" {
    const result = try nitrogen(testCommon(2), testZones(), testZones(), testParameters(), testParameters());
    try std.testing.expect(result.ammonium.non_band_g < 0);
    try std.testing.expectApproxEqAbs(@as(f64, 0), result.nitrate.total(), 1e-15);
}

test "sequential phosphate immobilization conserves phosphorus" {
    const result = try phosphorus(testCommon(0.2), testZones(), testZones(), testParameters(), testParameters());
    var state: PhosphorusState = .{ .microbial_p_g = 0.2, .non_band_dihydrogen_phosphate_g_p = 2, .band_dihydrogen_phosphate_g_p = 1, .non_band_hydrogen_phosphate_g_p = 2, .band_hydrogen_phosphate_g_p = 1 };
    const before = state.microbial_p_g + state.non_band_dihydrogen_phosphate_g_p + state.band_dihydrogen_phosphate_g_p + state.non_band_hydrogen_phosphate_g_p + state.band_hydrogen_phosphate_g_p;
    try commitPhosphorus(&state, result);
    try std.testing.expectApproxEqAbs(before, state.microbial_p_g + state.non_band_dihydrogen_phosphate_g_p + state.band_dihydrogen_phosphate_g_p + state.non_band_hydrogen_phosphate_g_p + state.band_hydrogen_phosphate_g_p, 1e-14);
}
