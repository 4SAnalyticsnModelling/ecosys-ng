const std = @import("std");
const nutrient_exchange = @import("nutrient_exchange.zig");

pub const SurfacePool = struct {
    non_band_fraction: f64,
    band_fraction: f64,
    non_band_concentration_g_per_m3: f64,
    band_concentration_g_per_m3: f64,
    total_dissolved_g: f64,
    competition_fraction: f64,
};

pub const BiologicalCapacity = struct {
    microbial_surface_area_m2_per_g_c: f64,
    active_biomass_g_c: f64,
    temperature_water_response: f64,
    timestep_h: f64,
};

pub const Flux = struct { actual_g: f64, supply_unlimited_g: f64, remaining_demand_before_g: f64 };

pub const NitrogenInputs = struct {
    initial_demand_g_n: f64,
    soil_ammonium_exchange_g_n: f64,
    soil_nitrate_exchange_g_n: f64,
    surface_water_volume_m3: f64,
    ammonium: SurfacePool,
    nitrate: SurfacePool,
    biological_capacity: BiologicalCapacity,
};

pub const NitrogenResult = struct { ammonium: Flux, nitrate: Flux, unmet_demand_g_n: f64 };

/// Ports NITRO's L=0 litter exchange. The NO3 capacity deliberately retains
/// the source AMAX1(demand, biological capacity), unlike NH4 and phosphate.
pub fn nitrogen(inputs: NitrogenInputs, ammonium_parameters: nutrient_exchange.UptakeParameters, nitrate_parameters: nutrient_exchange.UptakeParameters) !NitrogenResult {
    try validateNitrogen(inputs);
    try validateParameters(ammonium_parameters);
    try validateParameters(nitrate_parameters);
    const ammonium_demand = inputs.initial_demand_g_n - inputs.soil_ammonium_exchange_g_n - inputs.soil_nitrate_exchange_g_n;
    const ammonium = surfaceFlux(ammonium_demand, inputs.surface_water_volume_m3, inputs.ammonium, inputs.biological_capacity, ammonium_parameters, false);
    const nitrate_demand = @max(0, ammonium_demand - ammonium.actual_g);
    const nitrate = surfaceFlux(nitrate_demand, inputs.surface_water_volume_m3, inputs.nitrate, inputs.biological_capacity, nitrate_parameters, true);
    return .{ .ammonium = ammonium, .nitrate = nitrate, .unmet_demand_g_n = @max(0, nitrate_demand - nitrate.actual_g) };
}

pub const PhosphorusInputs = struct {
    initial_demand_g_p: f64,
    soil_dihydrogen_phosphate_exchange_g_p: f64,
    surface_water_volume_m3: f64,
    dihydrogen_phosphate: SurfacePool,
    hydrogen_phosphate: SurfacePool,
    biological_capacity: BiologicalCapacity,
};

pub const PhosphorusResult = struct { dihydrogen_phosphate: Flux, hydrogen_phosphate: Flux, unmet_demand_g_p: f64 };

pub fn phosphorus(inputs: PhosphorusInputs, dihydrogen_parameters: nutrient_exchange.UptakeParameters, hydrogen_parameters: nutrient_exchange.UptakeParameters) !PhosphorusResult {
    try validatePhosphorus(inputs);
    try validateParameters(dihydrogen_parameters);
    try validateParameters(hydrogen_parameters);
    const dihydrogen_demand = inputs.initial_demand_g_p - inputs.soil_dihydrogen_phosphate_exchange_g_p;
    const dihydrogen = surfaceFlux(dihydrogen_demand, inputs.surface_water_volume_m3, inputs.dihydrogen_phosphate, inputs.biological_capacity, dihydrogen_parameters, false);
    const hydrogen_demand = @max(0, dihydrogen_demand - dihydrogen.actual_g);
    const hydrogen = surfaceFlux(hydrogen_demand, inputs.surface_water_volume_m3, inputs.hydrogen_phosphate, inputs.biological_capacity, hydrogen_parameters, false);
    return .{ .dihydrogen_phosphate = dihydrogen, .hydrogen_phosphate = hydrogen, .unmet_demand_g_p = @max(0, hydrogen_demand - hydrogen.actual_g) };
}

fn surfaceFlux(demand_g: f64, water_volume_m3: f64, pool: SurfacePool, biological: BiologicalCapacity, parameters: nutrient_exchange.UptakeParameters, use_fortran_no3_maximum: bool) Flux {
    if (demand_g <= 0) return .{ .actual_g = demand_g, .supply_unlimited_g = 0, .remaining_demand_before_g = demand_g };
    const non_band_available = @max(0, pool.non_band_concentration_g_per_m3 - parameters.minimum_concentration_g_per_m3);
    const band_available = @max(0, pool.band_concentration_g_per_m3 - parameters.minimum_concentration_g_per_m3);
    const biological_limit = biological.microbial_surface_area_m2_per_g_c * biological.active_biomass_g_c * biological.temperature_water_response * parameters.maximum_uptake_g_per_m2_h * biological.timestep_h;
    const capacity_basis = if (use_fortran_no3_maximum) @max(demand_g, biological_limit) else @min(demand_g, biological_limit);
    const weighted_monod = pool.non_band_fraction * non_band_available / (non_band_available + parameters.half_saturation_g_per_m3) + pool.band_fraction * band_available / (band_available + parameters.half_saturation_g_per_m3);
    const unlimited = capacity_basis * weighted_monod;
    const protected_g = parameters.minimum_concentration_g_per_m3 * water_volume_m3;
    const actual = @min(pool.competition_fraction * @max(0, pool.total_dissolved_g - protected_g), unlimited);
    return .{ .actual_g = actual, .supply_unlimited_g = unlimited, .remaining_demand_before_g = demand_g };
}

pub const State = struct {
    litter_microbial_n_g: f64,
    litter_microbial_p_g: f64,
    surface_ammonium_g_n: f64,
    surface_nitrate_g_n: f64,
    surface_dihydrogen_phosphate_g_p: f64,
    surface_hydrogen_phosphate_g_p: f64,
};

pub fn commit(state: *State, nitrogen_result: NitrogenResult, phosphorus_result: PhosphorusResult) !void {
    inline for (@typeInfo(State).@"struct".fields) |field| if (!std.math.isFinite(@field(state.*, field.name)) or @field(state.*, field.name) < 0) return error.InvalidSurfaceNutrientState;
    const n_flux = nitrogen_result.ammonium.actual_g + nitrogen_result.nitrate.actual_g;
    const p_flux = phosphorus_result.dihydrogen_phosphate.actual_g + phosphorus_result.hydrogen_phosphate.actual_g;
    inline for (.{ nitrogen_result.ammonium.actual_g, nitrogen_result.nitrate.actual_g, phosphorus_result.dihydrogen_phosphate.actual_g, phosphorus_result.hydrogen_phosphate.actual_g }) |flux| if (!std.math.isFinite(flux)) return error.NonFiniteSurfaceNutrientFlux;
    if (state.surface_ammonium_g_n < nitrogen_result.ammonium.actual_g or state.surface_nitrate_g_n < nitrogen_result.nitrate.actual_g or state.surface_dihydrogen_phosphate_g_p < phosphorus_result.dihydrogen_phosphate.actual_g or state.surface_hydrogen_phosphate_g_p < phosphorus_result.hydrogen_phosphate.actual_g) return error.InsufficientSurfaceNutrient;
    if (state.litter_microbial_n_g + n_flux < 0 or state.litter_microbial_p_g + p_flux < 0) return error.NegativeLitterMicrobialNutrient;
    state.litter_microbial_n_g += n_flux;
    state.litter_microbial_p_g += p_flux;
    state.surface_ammonium_g_n -= nitrogen_result.ammonium.actual_g;
    state.surface_nitrate_g_n -= nitrogen_result.nitrate.actual_g;
    state.surface_dihydrogen_phosphate_g_p -= phosphorus_result.dihydrogen_phosphate.actual_g;
    state.surface_hydrogen_phosphate_g_p -= phosphorus_result.hydrogen_phosphate.actual_g;
}

fn validatePool(pool: SurfacePool) !void {
    inline for (@typeInfo(SurfacePool).@"struct".fields) |field| if (!std.math.isFinite(@field(pool, field.name)) or @field(pool, field.name) < 0) return error.InvalidSurfaceNutrientPool;
    if (@abs(pool.non_band_fraction + pool.band_fraction - 1) > 1e-12) return error.SurfaceNutrientFractionsDoNotSumToOne;
}

fn validateBiological(value: BiologicalCapacity) !void {
    inline for (@typeInfo(BiologicalCapacity).@"struct".fields) |field| if (!std.math.isFinite(@field(value, field.name)) or @field(value, field.name) < 0) return error.InvalidSurfaceNutrientBiology;
    if (value.timestep_h <= 0) return error.InvalidSurfaceNutrientBiology;
}

fn validateParameters(value: nutrient_exchange.UptakeParameters) !void {
    inline for (@typeInfo(nutrient_exchange.UptakeParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(value, field.name)) or @field(value, field.name) < 0) return error.InvalidSurfaceNutrientParameter;
    if (value.half_saturation_g_per_m3 <= 0) return error.InvalidSurfaceNutrientParameter;
}

fn validateNitrogen(inputs: NitrogenInputs) !void {
    inline for (@typeInfo(NitrogenInputs).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSurfaceNitrogenInput;
    if (inputs.surface_water_volume_m3 < 0) return error.InvalidSurfaceNitrogenInput;
    try validatePool(inputs.ammonium);
    try validatePool(inputs.nitrate);
    try validateBiological(inputs.biological_capacity);
}

fn validatePhosphorus(inputs: PhosphorusInputs) !void {
    inline for (@typeInfo(PhosphorusInputs).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(inputs, field.name))) return error.NonFiniteSurfacePhosphorusInput;
    if (inputs.surface_water_volume_m3 < 0) return error.InvalidSurfacePhosphorusInput;
    try validatePool(inputs.dihydrogen_phosphate);
    try validatePool(inputs.hydrogen_phosphate);
    try validateBiological(inputs.biological_capacity);
}

fn testPool(amount_g: f64) SurfacePool {
    return .{ .non_band_fraction = 0.7, .band_fraction = 0.3, .non_band_concentration_g_per_m3 = 2, .band_concentration_g_per_m3 = 1, .total_dissolved_g = amount_g, .competition_fraction = 0.5 };
}

fn testBiology() BiologicalCapacity {
    return .{ .microbial_surface_area_m2_per_g_c = 1, .active_biomass_g_c = 1, .temperature_water_response = 1, .timestep_h = 1 };
}

fn testParameters() nutrient_exchange.UptakeParameters {
    return .{ .minimum_concentration_g_per_m3 = 0.01, .half_saturation_g_per_m3 = 0.1, .maximum_uptake_g_per_m2_h = 0.2 };
}

test "surface litter exchange conserves nitrogen and phosphorus" {
    const n = try nitrogen(.{ .initial_demand_g_n = 1, .soil_ammonium_exchange_g_n = 0.1, .soil_nitrate_exchange_g_n = 0.1, .surface_water_volume_m3 = 1, .ammonium = testPool(2), .nitrate = testPool(2), .biological_capacity = testBiology() }, testParameters(), testParameters());
    const p = try phosphorus(.{ .initial_demand_g_p = 0.5, .soil_dihydrogen_phosphate_exchange_g_p = 0.1, .surface_water_volume_m3 = 1, .dihydrogen_phosphate = testPool(2), .hydrogen_phosphate = testPool(2), .biological_capacity = testBiology() }, testParameters(), testParameters());
    var state: State = .{ .litter_microbial_n_g = 1, .litter_microbial_p_g = 1, .surface_ammonium_g_n = 2, .surface_nitrate_g_n = 2, .surface_dihydrogen_phosphate_g_p = 2, .surface_hydrogen_phosphate_g_p = 2 };
    const before_n = state.litter_microbial_n_g + state.surface_ammonium_g_n + state.surface_nitrate_g_n;
    const before_p = state.litter_microbial_p_g + state.surface_dihydrogen_phosphate_g_p + state.surface_hydrogen_phosphate_g_p;
    try commit(&state, n, p);
    try std.testing.expectApproxEqAbs(before_n, state.litter_microbial_n_g + state.surface_ammonium_g_n + state.surface_nitrate_g_n, 1e-14);
    try std.testing.expectApproxEqAbs(before_p, state.litter_microbial_p_g + state.surface_dihydrogen_phosphate_g_p + state.surface_hydrogen_phosphate_g_p, 1e-14);
}

test "surface nitrate preserves source AMAX1 capacity branch" {
    const result = try nitrogen(.{ .initial_demand_g_n = 1, .soil_ammonium_exchange_g_n = 0, .soil_nitrate_exchange_g_n = 0, .surface_water_volume_m3 = 1, .ammonium = testPool(0.01), .nitrate = testPool(10), .biological_capacity = testBiology() }, testParameters(), testParameters());
    try std.testing.expect(result.nitrate.supply_unlimited_g > 0.2);
}
