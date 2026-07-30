const std = @import("std");

pub const AcetotrophicInputs = struct {
    soil_temperature_k: f64,
    aqueous_acetate_concentration_g_c_per_m3: f64,
    aqueous_acetate_g_c: f64,
    acetate_competition_fraction: f64,
    nutrient_limitation_fraction: f64,
    water_stress_fraction: f64,
    temperature_response: f64,
    active_biomass_g_c: f64,
    timestep_h: f64,
};

pub const AcetotrophicParameters = struct {
    acetate_inhibition_concentration_g_c_per_m3: f64,
    acetate_half_saturation_g_c_per_m3: f64,
    specific_respiration_per_h: f64,
    reference_energy_yield_kj_per_g_c: f64,
    growth_energy_requirement_kj_per_g_c: f64,
    minimum_growth_respiration_fraction: f64,
};

pub const AcetotrophicResult = struct {
    acetate_oxidation_g_c: f64,
    methane_production_g_c: f64,
    growth_respiration_fraction: f64,
    substrate_unlimited_oxidation_g_c: f64,
};

/// Acetotrophic NITRO branch: acetate feedback changes energy yield, while
/// biomass kinetics and the finite acetate pool independently limit oxidation.
pub fn acetotrophic(inputs: AcetotrophicInputs, parameters: AcetotrophicParameters) !AcetotrophicResult {
    try validateStruct(inputs, error.InvalidAcetotrophicMethanogenesisInput);
    try validateStruct(parameters, error.InvalidAcetotrophicMethanogenesisParameter);
    if (inputs.soil_temperature_k <= 0 or inputs.aqueous_acetate_concentration_g_c_per_m3 < 0 or inputs.aqueous_acetate_g_c < 0 or inputs.acetate_competition_fraction < 0 or inputs.nutrient_limitation_fraction < 0 or inputs.water_stress_fraction < 0 or inputs.temperature_response < 0 or inputs.active_biomass_g_c < 0 or inputs.timestep_h <= 0 or parameters.acetate_inhibition_concentration_g_c_per_m3 <= 0 or parameters.acetate_half_saturation_g_c_per_m3 <= 0 or parameters.specific_respiration_per_h < 0 or parameters.growth_energy_requirement_kj_per_g_c <= 0 or parameters.minimum_growth_respiration_fraction < 0 or parameters.minimum_growth_respiration_fraction > 1) return error.InvalidAcetotrophicMethanogenesis;
    const feedback_kj_per_g_c = 8.3143e-3 * inputs.soil_temperature_k * @log(@max(std.math.floatMin(f64), inputs.aqueous_acetate_concentration_g_c_per_m3 / parameters.acetate_inhibition_concentration_g_c_per_m3)) / 24;
    const respiration_fraction = @max(parameters.minimum_growth_respiration_fraction, @min(1, 1 / (1 + @max(0, parameters.reference_energy_yield_kj_per_g_c + feedback_kj_per_g_c) / parameters.growth_energy_requirement_kj_per_g_c)));
    const monod = inputs.aqueous_acetate_concentration_g_c_per_m3 / (inputs.aqueous_acetate_concentration_g_c_per_m3 + parameters.acetate_half_saturation_g_c_per_m3);
    const unlimited = @max(0, parameters.specific_respiration_per_h * inputs.nutrient_limitation_fraction * inputs.water_stress_fraction * inputs.active_biomass_g_c * inputs.timestep_h);
    const kinetic_limit = unlimited * monod * inputs.temperature_response;
    const supply_limit = @max(0, inputs.aqueous_acetate_g_c * inputs.acetate_competition_fraction * respiration_fraction * inputs.timestep_h);
    const oxidation = @min(kinetic_limit, supply_limit);
    return .{ .acetate_oxidation_g_c = oxidation, .methane_production_g_c = 0.5 * oxidation, .growth_respiration_fraction = respiration_fraction, .substrate_unlimited_oxidation_g_c = unlimited };
}

pub const HydrogenotrophicInputs = struct {
    aqueous_hydrogen_concentration_g_h_per_m3: f64,
    aqueous_hydrogen_g_h: f64,
    fermentation_hydrogen_production_g_h: f64,
    temperature_water_response: f64,
    nutrient_limitation_fraction: f64,
    aqueous_co2_limitation_fraction: f64,
    active_biomass_g_c: f64,
    timestep_h: f64,
    hydrogen_feedback_energy_kj_per_mol: f64,
};

pub const HydrogenotrophicParameters = struct {
    hydrogen_half_saturation_g_h_per_m3: f64,
    specific_co2_reduction_g_c_per_g_c_h: f64,
    reference_energy_yield_kj_per_g_c: f64,
    growth_energy_requirement_kj_per_g_c: f64,
    minimum_growth_respiration_fraction: f64,
    hydrogen_supply_conversion_g_c_per_g_h: f64,
    fermentation_hydrogen_to_pool_fraction: f64,
};

pub const HydrogenotrophicResult = struct { co2_reduction_g_c: f64, methane_production_g_c: f64, growth_respiration_fraction: f64 };

pub fn hydrogenotrophic(inputs: HydrogenotrophicInputs, parameters: HydrogenotrophicParameters) !HydrogenotrophicResult {
    try validateStruct(inputs, error.InvalidHydrogenotrophicMethanogenesisInput);
    try validateStruct(parameters, error.InvalidHydrogenotrophicMethanogenesisParameter);
    if (inputs.aqueous_hydrogen_concentration_g_h_per_m3 < 0 or inputs.aqueous_hydrogen_g_h < 0 or inputs.fermentation_hydrogen_production_g_h < 0 or inputs.temperature_water_response < 0 or inputs.nutrient_limitation_fraction < 0 or inputs.aqueous_co2_limitation_fraction < 0 or inputs.active_biomass_g_c < 0 or inputs.timestep_h <= 0 or parameters.hydrogen_half_saturation_g_h_per_m3 <= 0 or parameters.specific_co2_reduction_g_c_per_g_c_h < 0 or parameters.growth_energy_requirement_kj_per_g_c <= 0 or parameters.minimum_growth_respiration_fraction < 0 or parameters.minimum_growth_respiration_fraction > 1 or parameters.hydrogen_supply_conversion_g_c_per_g_h <= 0 or parameters.fermentation_hydrogen_to_pool_fraction < 0) return error.InvalidHydrogenotrophicMethanogenesis;
    const respiration_fraction = @max(parameters.minimum_growth_respiration_fraction, @min(1, 1 / (1 + @max(0, parameters.reference_energy_yield_kj_per_g_c + inputs.hydrogen_feedback_energy_kj_per_mol / 12) / parameters.growth_energy_requirement_kj_per_g_c)));
    const unlimited = parameters.specific_co2_reduction_g_c_per_g_c_h * inputs.temperature_water_response * inputs.nutrient_limitation_fraction * inputs.aqueous_co2_limitation_fraction * inputs.active_biomass_g_c * inputs.timestep_h;
    const monod = inputs.aqueous_hydrogen_concentration_g_h_per_m3 / (inputs.aqueous_hydrogen_concentration_g_h_per_m3 + parameters.hydrogen_half_saturation_g_h_per_m3);
    const hydrogen_supply_g_h = inputs.aqueous_hydrogen_g_h + parameters.fermentation_hydrogen_to_pool_fraction * inputs.fermentation_hydrogen_production_g_h;
    const reduction = @max(0, @min(unlimited * monod, parameters.hydrogen_supply_conversion_g_c_per_g_h * hydrogen_supply_g_h * inputs.timestep_h));
    return .{ .co2_reduction_g_c = reduction, .methane_production_g_c = reduction, .growth_respiration_fraction = respiration_fraction };
}

fn validateStruct(value: anytype, comptime failure: anyerror) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| if (field.type == f64 and !std.math.isFinite(@field(value, field.name))) return failure;
}

test "acetotrophic methane is half of acetate oxidation" {
    const result = try acetotrophic(.{ .soil_temperature_k = 293, .aqueous_acetate_concentration_g_c_per_m3 = 2, .aqueous_acetate_g_c = 5, .acetate_competition_fraction = 0.8, .nutrient_limitation_fraction = 0.9, .water_stress_fraction = 1, .temperature_response = 0.8, .active_biomass_g_c = 1, .timestep_h = 1 }, .{ .acetate_inhibition_concentration_g_c_per_m3 = 1, .acetate_half_saturation_g_c_per_m3 = 0.5, .specific_respiration_per_h = 0.2, .reference_energy_yield_kj_per_g_c = 1, .growth_energy_requirement_kj_per_g_c = 10, .minimum_growth_respiration_fraction = 0.05 });
    try std.testing.expect(result.acetate_oxidation_g_c > 0);
    try std.testing.expectApproxEqAbs(0.5 * result.acetate_oxidation_g_c, result.methane_production_g_c, 1e-15);
}

test "hydrogenotrophic methane respects finite hydrogen supply" {
    const result = try hydrogenotrophic(.{ .aqueous_hydrogen_concentration_g_h_per_m3 = 1, .aqueous_hydrogen_g_h = 0.01, .fermentation_hydrogen_production_g_h = 0, .temperature_water_response = 1, .nutrient_limitation_fraction = 1, .aqueous_co2_limitation_fraction = 1, .active_biomass_g_c = 10, .timestep_h = 1, .hydrogen_feedback_energy_kj_per_mol = 0 }, .{ .hydrogen_half_saturation_g_h_per_m3 = 0.1, .specific_co2_reduction_g_c_per_g_c_h = 1, .reference_energy_yield_kj_per_g_c = 1, .growth_energy_requirement_kj_per_g_c = 10, .minimum_growth_respiration_fraction = 0.05, .hydrogen_supply_conversion_g_c_per_g_h = 1.5, .fermentation_hydrogen_to_pool_fraction = 0.111 });
    try std.testing.expectApproxEqAbs(@as(f64, 0.015), result.methane_production_g_c, 1e-15);
}
