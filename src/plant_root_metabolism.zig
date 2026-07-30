const std = @import("std");
const RootState = @import("plant_root_system.zig").State;
const NutrientResult = @import("plant_root_nutrient_uptake.zig").Result;
const root_domain_count = @import("plant_root_system.zig").biological_domain_count;

pub const PrimaryRootAxisScaling = struct {
    retained_root_carbon_g_c_per_plant: f64,
    primary_axis_count_multiplier: f64,
};

/// GROSUB WTRTA/XRTN1 retained root-mass state and primary-axis scaling.
/// The multiplication by the biological timestep preserves the source
/// operation exactly; it is not rewritten as an exponential decay.
pub fn primaryRootAxisScaling(
    previous_retained_root_carbon_g_c_per_plant: f64,
    total_root_carbon_g_c: f64,
    plant_population_count: f64,
    biological_timestep_h: f64,
) !PrimaryRootAxisScaling {
    inline for (.{ previous_retained_root_carbon_g_c_per_plant, total_root_carbon_g_c, plant_population_count, biological_timestep_h }) |value| {
        if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootAxisScalingInput;
    }
    if (biological_timestep_h <= 0) return error.InvalidPrimaryRootAxisScalingInput;

    const retained = if (plant_population_count > 0)
        @max(
            0.999992087 * previous_retained_root_carbon_g_c_per_plant * biological_timestep_h,
            total_root_carbon_g_c / plant_population_count,
        )
    else
        0.0;
    const axis_multiplier = @max(1.0, std.math.pow(f64, retained, 0.667)) * plant_population_count;
    if (!std.math.isFinite(retained) or !std.math.isFinite(axis_multiplier)) return error.NonFinitePrimaryRootAxisScaling;
    return .{
        .retained_root_carbon_g_c_per_plant = retained,
        .primary_axis_count_multiplier = axis_multiplier,
    };
}

pub const Respiration = struct {
    actual_g_c: f64,
    oxygen_unlimited_g_c: f64,
    carbon_unlimited_g_c: f64,
};

pub const Components = struct {
    maintenance_demand_g_c: f64,
    substrate_respiration_actual_g_c: f64,
    substrate_respiration_oxygen_unlimited_g_c: f64,
    growth_respiration_actual_g_c: f64,
    growth_respiration_oxygen_unlimited_g_c: f64,
    senescence_respiration_actual_g_c: f64,
    senescence_respiration_oxygen_unlimited_g_c: f64,
    nitrogen_assimilation_respiration_actual_g_c: f64,
    nitrogen_assimilation_respiration_oxygen_unlimited_g_c: f64,
};

pub const SecondaryRootParameters = struct {
    maximum_substrate_respiration_fraction_per_h: f64,
    substrate_respiration_half_saturation_g_c_per_g_c: f64,
    nitrogen_feedback_half_saturation_g_n_per_g_c: f64,
    phosphorus_feedback_half_saturation_g_p_per_g_c: f64,
    maintenance_respiration_g_c_per_g_n_h: f64,
    nitrogen_assimilation_respiration_g_c_per_g_n: f64,
    minimum_carbon_recycling_fraction: f64,
    responsive_carbon_recycling_fraction: f64,
    maximum_nitrogen_recycling_fraction: f64,
    maximum_phosphorus_recycling_fraction: f64,
    storage_exchange_fraction_per_h: f64,
    nonwoody_root_fraction_exponent: f64,
    maintenance_gas_constant_j_per_mol_k: f64,
    maintenance_enthalpy_j_per_mol_k: f64,
    maintenance_activation_energy_j_per_mol: f64,
    maintenance_low_temperature_inactivation_energy_j_per_mol: f64,
    maintenance_normalization_log_intercept: f64,
    maximum_maintenance_temperature_response: f64,
    shallow_root_water_response_per_mpa: f64,
    deep_root_water_response_per_mpa: f64,
    maintenance_water_response_exponent: f64,
    root_penetration_reference_radius_m: f64,
    acidity_half_effect_hydrogen_activity_mol_per_m3: f64,
    maximum_acidity_enhancement: f64,
    shallow_primary_root_sink_multiplier: f64,
    intermediate_primary_root_sink_multiplier: f64,
    deep_primary_root_sink_multiplier: f64,
    deeper_primary_root_sink_multiplier: f64,
    annual_termination_hours_without_grain_fill: f64,
    root_protein_carbon_per_nitrogen_g_c_per_g_n: f64,
    root_protein_carbon_per_phosphorus_g_c_per_g_p: f64,
    nutrient_uptake_respiration_g_c_per_g_element: f64,
    evergreen_leafoff_remobilization_start_fraction: f64,
    deciduous_leafoff_remobilization_start_fraction: f64,
    full_senescence_duration_h: f64,

    pub fn validate(self: SecondaryRootParameters) !void {
        inline for (@typeInfo(SecondaryRootParameters).@"struct".fields) |field| {
            const value = @field(self, field.name);
            if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootMetabolismParameter;
        }
        if (self.substrate_respiration_half_saturation_g_c_per_g_c == 0) return error.InvalidSecondaryRootMetabolismParameter;
        if (self.nitrogen_feedback_half_saturation_g_n_per_g_c == 0) return error.InvalidSecondaryRootMetabolismParameter;
        if (self.phosphorus_feedback_half_saturation_g_p_per_g_c == 0) return error.InvalidSecondaryRootMetabolismParameter;
        if (self.minimum_carbon_recycling_fraction > 1 or self.responsive_carbon_recycling_fraction > 1 or
            self.minimum_carbon_recycling_fraction + self.responsive_carbon_recycling_fraction > 1 or
            self.maximum_nitrogen_recycling_fraction > 1 or self.maximum_phosphorus_recycling_fraction > 1 or
            self.evergreen_leafoff_remobilization_start_fraction > 1 or self.deciduous_leafoff_remobilization_start_fraction > 1)
            return error.InvalidSecondaryRootMetabolismParameter;
        if (self.maintenance_gas_constant_j_per_mol_k == 0 or self.root_penetration_reference_radius_m == 0 or
            self.acidity_half_effect_hydrogen_activity_mol_per_m3 == 0)
            return error.InvalidSecondaryRootMetabolismParameter;
        if (self.full_senescence_duration_h == 0) return error.InvalidSecondaryRootMetabolismParameter;
    }
};

/// Historical GROSUB values, selected at runtime only when an older runscript
/// has not yet supplied the root_metabolism record.
pub fn compatibilitySecondaryRootParameters() SecondaryRootParameters {
    return .{
        .maximum_substrate_respiration_fraction_per_h = 0.015,
        .substrate_respiration_half_saturation_g_c_per_g_c = 0.025,
        .nitrogen_feedback_half_saturation_g_n_per_g_c = 0.1,
        .phosphorus_feedback_half_saturation_g_p_per_g_c = 0.01,
        .maintenance_respiration_g_c_per_g_n_h = 0.010,
        .nitrogen_assimilation_respiration_g_c_per_g_n = 1.70,
        .minimum_carbon_recycling_fraction = 0.167,
        .responsive_carbon_recycling_fraction = 0.333,
        .maximum_nitrogen_recycling_fraction = 0.667,
        .maximum_phosphorus_recycling_fraction = 0.667,
        .storage_exchange_fraction_per_h = 2.5e-5,
        .nonwoody_root_fraction_exponent = 0.167,
        .maintenance_gas_constant_j_per_mol_k = 8.3143,
        .maintenance_enthalpy_j_per_mol_k = 710,
        .maintenance_activation_energy_j_per_mol = 62500,
        .maintenance_low_temperature_inactivation_energy_j_per_mol = 197500,
        .maintenance_normalization_log_intercept = 25.216,
        .maximum_maintenance_temperature_response = 1.0e3,
        .shallow_root_water_response_per_mpa = 0.05,
        .deep_root_water_response_per_mpa = 0.10,
        .maintenance_water_response_exponent = 0.25,
        .root_penetration_reference_radius_m = 1.0e-3,
        .acidity_half_effect_hydrogen_activity_mol_per_m3 = 1,
        .maximum_acidity_enhancement = 4,
        .shallow_primary_root_sink_multiplier = 0.25,
        .intermediate_primary_root_sink_multiplier = 1,
        .deep_primary_root_sink_multiplier = 2,
        .deeper_primary_root_sink_multiplier = 4,
        .annual_termination_hours_without_grain_fill = 336,
        .root_protein_carbon_per_nitrogen_g_c_per_g_n = 2.5,
        .root_protein_carbon_per_phosphorus_g_c_per_g_p = 25,
        .nutrient_uptake_respiration_g_c_per_g_element = 0.86,
        .evergreen_leafoff_remobilization_start_fraction = 0.75,
        .deciduous_leafoff_remobilization_start_fraction = 0.5,
        .full_senescence_duration_h = 480,
    };
}

pub const SecondaryRootInputs = struct {
    mobile_carbon_g_c: f64,
    nonstructural_nitrogen_g_n: f64,
    nonstructural_phosphorus_g_p: f64,
    root_carbon_g_c: f64,
    root_nitrogen_g_n: f64,
    root_nitrogen_to_carbon_ratio_g_n_per_g_c: f64,
    root_phosphorus_to_carbon_ratio_g_p_per_g_c: f64,
    root_growth_yield_g_c_per_g_c: f64,
    active_root_fraction: f64,
    biological_timestep_h: f64,
    substrate_temperature_response: f64,
    maintenance_temperature_response: f64,
    acidity_response: f64,
    substrate_feedback: f64,
    oxygen_limitation: f64,
    substrate_water_response: f64,
    maintenance_water_response: f64,
};

/// Per-PFT values needed by the live GROSUB root kernel. These are copied
/// from the runtime plant catalog once, so hourly tiles never parse or retain
/// dependencies on input files.
pub const RuntimePlantParameters = struct {
    root_profile_type: u8,
    mycorrhizal_type: u8,
    growth_habit: u8,
    leaf_phenology_type: u8,
    root_growth_yield_g_c_per_g_c: f64,
    root_nitrogen_to_carbon_g_n_per_g_c: f64,
    root_phosphorus_to_carbon_g_p_per_g_c: f64,
    stalk_nitrogen_to_carbon_g_n_per_g_c: f64,
    stalk_phosphorus_to_carbon_g_p_per_g_c: f64,
    primary_root_radius_m: f64,
    secondary_root_radius_m: f64,
    primary_specific_length_m_per_g_c: f64,
    secondary_specific_length_m_per_g_c: f64,
    secondary_root_branching_per_m: f64,
    shoot_root_equilibration_fraction_per_h: f64,

    pub fn validate(self: RuntimePlantParameters) !void {
        if (self.root_profile_type > 3 or self.mycorrhizal_type > 2 or self.growth_habit > 1 or self.leaf_phenology_type > 5) return error.InvalidRootMetabolismPlantCode;
        inline for (@typeInfo(RuntimePlantParameters).@"struct".fields) |field| {
            if (field.type == u8) continue;
            const value = @field(self, field.name);
            if (!std.math.isFinite(value) or value < 0 or
                (value == 0 and !std.mem.eql(u8, field.name, "shoot_root_equilibration_fraction_per_h")))
                return error.InvalidRootMetabolismPlantParameter;
        }
        if (self.root_growth_yield_g_c_per_g_c >= 1) return error.InvalidRootMetabolismPlantParameter;
    }

    /// READQ MY is the source upper bound for every `DO N=1,MY` root loop.
    pub fn biologicalDomainCount(self: RuntimePlantParameters) usize {
        return if (self.mycorrhizal_type == 2) 2 else 1;
    }
};

pub const SecondaryRootResult = struct {
    nutrient_feedback: f64,
    substrate_respiration_oxygen_unlimited_g_c_per_h: f64,
    maintenance_respiration_g_c_per_h: f64,
    substrate_respiration_actual_g_c_per_h: f64,
    growth_respiration_oxygen_unlimited_g_c_per_h: f64,
    growth_respiration_actual_g_c_per_h: f64,
    growth_and_respiration_carbon_oxygen_unlimited_g_c_per_h: f64,
    growth_and_respiration_carbon_actual_g_c_per_h: f64,
    root_growth_oxygen_unlimited_g_c_per_h: f64,
    root_growth_actual_g_c_per_h: f64,
    nitrogen_growth_demand_g_n_per_h: f64,
    nitrogen_growth_actual_g_n_per_h: f64,
    phosphorus_growth_actual_g_p_per_h: f64,
    nitrogen_assimilation_respiration_oxygen_unlimited_g_c_per_h: f64,
    nitrogen_assimilation_respiration_actual_g_c_per_h: f64,
};

pub const PrimaryRootInputs = struct {
    shared: SecondaryRootInputs,
    primary_tip_at_or_below_profile_bottom: bool,
};

/// GROSUB `RGFNP` with STARTQ `CNRTS=CNRT*DMRT` and
/// `CPRTS=CPRT*DMRT`. The returned quantity is growth respiration, not
/// structural root growth.
pub fn nutrientLimitedRootGrowthRespiration(
    nonstructural_nitrogen_g_n: f64,
    nonstructural_phosphorus_g_p: f64,
    active_root_fraction: f64,
    respiration_fraction: f64,
    growth_yield_g_c_per_g_c: f64,
    nitrogen_to_carbon_g_n_per_g_c: f64,
    phosphorus_to_carbon_g_p_per_g_c: f64,
) !f64 {
    inline for (.{
        nonstructural_nitrogen_g_n,
        nonstructural_phosphorus_g_p,
        active_root_fraction,
        respiration_fraction,
        growth_yield_g_c_per_g_c,
        nitrogen_to_carbon_g_n_per_g_c,
        phosphorus_to_carbon_g_p_per_g_c,
    }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootGrowthNutrientLimitInput;
    if (nonstructural_nitrogen_g_n < 0 or nonstructural_phosphorus_g_p < 0 or active_root_fraction < 0 or active_root_fraction > 1 or respiration_fraction <= 0 or growth_yield_g_c_per_g_c <= 0 or nitrogen_to_carbon_g_n_per_g_c <= 0 or phosphorus_to_carbon_g_p_per_g_c <= 0) return error.InvalidRootGrowthNutrientLimitInput;
    const result = @min(
        nonstructural_nitrogen_g_n * respiration_fraction * active_root_fraction /
            (nitrogen_to_carbon_g_n_per_g_c * growth_yield_g_c_per_g_c),
        nonstructural_phosphorus_g_p * respiration_fraction * active_root_fraction /
            (phosphorus_to_carbon_g_p_per_g_c * growth_yield_g_c_per_g_c),
    );
    if (!std.math.isFinite(result)) return error.NonFiniteRootGrowthNutrientLimitResult;
    return result;
}

pub const RootEnvironmentResponses = struct {
    maintenance_temperature: f64,
    acidity: f64,
    growth_water: f64,
    maintenance_water: f64,
    extension_water: f64,
    scaled_penetration_resistance_mpa: f64,
};

/// GROSUB 5985--5995 lower-layer scan. Layers at or below the minimum
/// thickness are skipped, except that the bottom layer is always selected.
pub fn nextLowerRootLayer(
    layer_thickness_m: []const f64,
    current_layer: usize,
    minimum_active_layer_thickness_m: f64,
) !usize {
    if (layer_thickness_m.len == 0 or current_layer >= layer_thickness_m.len)
        return error.RootLayerIndexOutOfBounds;
    if (!std.math.isFinite(minimum_active_layer_thickness_m) or minimum_active_layer_thickness_m < 0)
        return error.InvalidRootLayerThickness;
    for (layer_thickness_m) |thickness_m|
        if (!std.math.isFinite(thickness_m) or thickness_m < 0)
            return error.InvalidRootLayerThickness;
    if (current_layer + 1 >= layer_thickness_m.len) return error.NoLowerRootLayer;
    const bottom_layer = layer_thickness_m.len - 1;
    for (current_layer + 1..layer_thickness_m.len) |candidate| {
        if (layer_thickness_m[candidate] > minimum_active_layer_thickness_m or candidate == bottom_layer)
            return candidate;
    }
    unreachable;
}

/// GROSUB TFN6, WFNRT, WFNRG/WFNRR and HOUR1 FPH derived directly from
/// live dimensional soil/root state.
pub fn rootEnvironmentResponses(
    parameters: SecondaryRootParameters,
    soil_temperature_k: f64,
    thermal_adaptation_offset_k: f64,
    soil_ph: f64,
    root_total_water_potential_mpa: f64,
    root_turgor_water_potential_mpa: f64,
    minimum_extension_water_potential_mpa: f64,
    soil_penetration_resistance_mpa: f64,
    secondary_root_radius_m: f64,
    shallow_root_profile: bool,
) !RootEnvironmentResponses {
    try parameters.validate();
    inline for (.{ soil_temperature_k, thermal_adaptation_offset_k, soil_ph, root_total_water_potential_mpa, root_turgor_water_potential_mpa, minimum_extension_water_potential_mpa, soil_penetration_resistance_mpa, secondary_root_radius_m }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootEnvironmentInput;
    if (soil_temperature_k <= 0 or soil_penetration_resistance_mpa < 0 or secondary_root_radius_m < 0) return error.InvalidRootEnvironmentInput;
    const adjusted_temperature_k = soil_temperature_k + thermal_adaptation_offset_k;
    if (adjusted_temperature_k <= 0) return error.InvalidRootEnvironmentInput;
    const gas_temperature_j_per_mol = parameters.maintenance_gas_constant_j_per_mol_k * adjusted_temperature_k;
    const enthalpy_temperature_j_per_mol = parameters.maintenance_enthalpy_j_per_mol_k * adjusted_temperature_k;
    const inactivation = 1 + std.math.exp((parameters.maintenance_low_temperature_inactivation_energy_j_per_mol - enthalpy_temperature_j_per_mol) / gas_temperature_j_per_mol);
    const maintenance_temperature = @min(
        parameters.maximum_maintenance_temperature_response,
        std.math.exp(parameters.maintenance_normalization_log_intercept - parameters.maintenance_activation_energy_j_per_mol / gas_temperature_j_per_mol) / inactivation,
    );
    const hydrogen_activity_mol_per_m3 = 1.0e3 * std.math.pow(f64, 10, -soil_ph);
    const acidity = 1 + @min(parameters.maximum_acidity_enhancement, hydrogen_activity_mol_per_m3 / parameters.acidity_half_effect_hydrogen_activity_mol_per_m3);
    const water_coefficient = if (shallow_root_profile) parameters.shallow_root_water_response_per_mpa else parameters.deep_root_water_response_per_mpa;
    const growth_water = std.math.exp(water_coefficient * root_total_water_potential_mpa);
    const maintenance_water = std.math.pow(f64, growth_water, parameters.maintenance_water_response_exponent);
    const radius_ratio = secondary_root_radius_m / parameters.root_penetration_reference_radius_m;
    const scaled_resistance = soil_penetration_resistance_mpa * radius_ratio * radius_ratio;
    const extension_water = std.math.clamp(root_turgor_water_potential_mpa - minimum_extension_water_potential_mpa - scaled_resistance, 0, 1);
    const result: RootEnvironmentResponses = .{
        .maintenance_temperature = maintenance_temperature,
        .acidity = acidity,
        .growth_water = growth_water,
        .maintenance_water = maintenance_water,
        .extension_water = extension_water,
        .scaled_penetration_resistance_mpa = scaled_resistance,
    };
    inline for (@typeInfo(RootEnvironmentResponses).@"struct".fields) |field| if (!std.math.isFinite(@field(result, field.name)) or @field(result, field.name) < 0) return error.NonFiniteRootEnvironmentResponse;
    return result;
}

pub const RootAxisSinkInputs = struct {
    root_profile_type: u8,
    primary_axis_count_multiplier: f64,
    primary_root_radius_m: f64,
    primary_root_depth_from_canopy_m: f64,
    secondary_root_depth_from_canopy_m: f64,
    secondary_axis_count: f64,
    secondary_root_radius_m: f64,
    average_secondary_root_length_m: f64,
    primary_biological_domain: bool,
};

pub const SourceOrderRootAxisSinkInputs = struct {
    root_profile_type: u8,
    primary_axis_count_multiplier: f64,
    primary_root_radius_m: f64,
    primary_root_depth_from_surface_m: f64,
    layer_top_depth_m: f64,
    layer_thickness_m: f64,
    secondary_root_origin_offset_m: f64,
    seeding_depth_m: f64,
    hypocotyledon_height_m: f64,
    canopy_height_m: f64,
    secondary_axis_count: f64,
    secondary_root_radius_m: f64,
    average_secondary_root_length_m: f64,
    negligible_sink_m: f64,
    primary_biological_domain: bool,
};

/// STOMATE FDBKX consumed by both shoot carboxylation and GROSUB root
/// substrate respiration.
pub fn annualTerminationFeedback(growth_habit: u8, hours_without_grain_fill: f64, termination_hours: f64) !f64 {
    inline for (.{ hours_without_grain_fill, termination_hours }) |value| if (!std.math.isFinite(value)) return error.NonFiniteAnnualTerminationInput;
    if (hours_without_grain_fill < 0 or termination_hours <= 0) return error.InvalidAnnualTerminationInput;
    return if (growth_habit == 0 and hours_without_grain_fill > 0)
        @max(0, 1 - hours_without_grain_fill / termination_hours)
    else
        1;
}

pub const RootAxisSinkStrength = struct {
    primary_m: f64,
    secondary_m: f64,
};

/// GROSUB 6042 secondary-root axis gate. Layer indexes are zero-based in
/// Zig; the source comparison remains inclusive.
pub fn secondaryRootAxisActive(
    current_layer: usize,
    deepest_secondary_root_layer: usize,
    axis_inactive: bool,
) bool {
    return current_layer <= deepest_secondary_root_layer and !axis_inactive;
}

/// GROSUB 6112--6129 annual physiological-maturity gate.
pub fn rootRespirationActive(
    perennial_growth_habit: bool,
    physiological_maturity_date_is_set: bool,
) bool {
    return perennial_growth_habit or !physiological_maturity_date_is_set;
}

pub const RootRespirationWaterResponses = struct {
    substrate: f64,
    maintenance: f64,
};

/// GROSUB 6115--6125 uses WFNRG for substrate respiration in every PFT.
/// Only maintenance switches to WFNRR outside shallow and drought-deciduous
/// plants.
pub fn sourceRootRespirationWaterResponses(
    root_profile_type: u8,
    leaf_phenology_type: u8,
    growth_water_response: f64,
    maintenance_water_response: f64,
) !RootRespirationWaterResponses {
    if (root_profile_type > 3 or leaf_phenology_type > 5)
        return error.InvalidRootMetabolismPlantCode;
    inline for (.{ growth_water_response, maintenance_water_response }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootEnvironmentResponse;
    return .{
        .substrate = growth_water_response,
        .maintenance = if (root_profile_type == 0 or leaf_phenology_type == 2)
            growth_water_response
        else
            maintenance_water_response,
    };
}

pub const AxisWorkspace = struct {
    allocator: std.mem.Allocator,
    axis_capacity: usize,
    sink_strengths: []RootAxisSinkStrength,
    primary_sink_fractions: []f64,
    secondary_sink_fractions: []f64,
    primary_metabolism: []SecondaryRootResult,
    secondary_metabolism: []SecondaryRootResult,
    primary_senescence: []SecondaryRootSenescence,
    secondary_senescence: []SecondaryRootSenescence,
    primary_deficit_absorption: []SecondaryRootDeficitAbsorption,
    primary_deficit_active: []bool,
    primary_active: []bool,
    secondary_active: []bool,
    primary_processed_by_domain: []bool,
    primary_respiration_allocation_fractions: []f64,
    primary_respiration_root_layer_indices: []usize,
    primary_length_m_by_layer: []f64,

    pub fn init(allocator: std.mem.Allocator, axis_capacity: usize, soil_layer_count: usize) !AxisWorkspace {
        if (axis_capacity == 0 or soil_layer_count == 0) return error.InvalidRootMetabolismWorkspaceDimensions;
        const strengths = try allocator.alloc(RootAxisSinkStrength, axis_capacity);
        errdefer allocator.free(strengths);
        const primary_fractions = try allocator.alloc(f64, axis_capacity);
        errdefer allocator.free(primary_fractions);
        const secondary_fractions = try allocator.alloc(f64, axis_capacity);
        errdefer allocator.free(secondary_fractions);
        const primary_metabolism = try allocator.alloc(SecondaryRootResult, axis_capacity);
        errdefer allocator.free(primary_metabolism);
        const secondary_metabolism = try allocator.alloc(SecondaryRootResult, axis_capacity);
        errdefer allocator.free(secondary_metabolism);
        const primary_senescence = try allocator.alloc(SecondaryRootSenescence, axis_capacity);
        errdefer allocator.free(primary_senescence);
        const secondary_senescence = try allocator.alloc(SecondaryRootSenescence, axis_capacity);
        errdefer allocator.free(secondary_senescence);
        const deficit_absorption = try allocator.alloc(SecondaryRootDeficitAbsorption, axis_capacity);
        errdefer allocator.free(deficit_absorption);
        const deficit_active = try allocator.alloc(bool, axis_capacity);
        errdefer allocator.free(deficit_active);
        const primary_active = try allocator.alloc(bool, axis_capacity);
        errdefer allocator.free(primary_active);
        const secondary_active = try allocator.alloc(bool, axis_capacity);
        errdefer allocator.free(secondary_active);
        const primary_processed = try allocator.alloc(bool, try std.math.mul(usize, axis_capacity, root_domain_count));
        errdefer allocator.free(primary_processed);
        const respiration_fractions = try allocator.alloc(f64, soil_layer_count);
        errdefer allocator.free(respiration_fractions);
        const respiration_indices = try allocator.alloc(usize, soil_layer_count);
        errdefer allocator.free(respiration_indices);
        const primary_lengths = try allocator.alloc(f64, soil_layer_count);
        @memset(strengths, .{ .primary_m = 0, .secondary_m = 0 });
        @memset(primary_fractions, 0);
        @memset(secondary_fractions, 0);
        @memset(primary_metabolism, std.mem.zeroes(SecondaryRootResult));
        @memset(secondary_metabolism, std.mem.zeroes(SecondaryRootResult));
        @memset(primary_senescence, std.mem.zeroes(SecondaryRootSenescence));
        @memset(secondary_senescence, std.mem.zeroes(SecondaryRootSenescence));
        @memset(deficit_absorption, std.mem.zeroes(SecondaryRootDeficitAbsorption));
        @memset(deficit_active, false);
        @memset(primary_active, false);
        @memset(secondary_active, false);
        @memset(primary_processed, false);
        @memset(respiration_fractions, 0);
        @memset(respiration_indices, 0);
        @memset(primary_lengths, 0);
        return .{
            .allocator = allocator,
            .axis_capacity = axis_capacity,
            .sink_strengths = strengths,
            .primary_sink_fractions = primary_fractions,
            .secondary_sink_fractions = secondary_fractions,
            .primary_metabolism = primary_metabolism,
            .secondary_metabolism = secondary_metabolism,
            .primary_senescence = primary_senescence,
            .secondary_senescence = secondary_senescence,
            .primary_deficit_absorption = deficit_absorption,
            .primary_deficit_active = deficit_active,
            .primary_active = primary_active,
            .secondary_active = secondary_active,
            .primary_processed_by_domain = primary_processed,
            .primary_respiration_allocation_fractions = respiration_fractions,
            .primary_respiration_root_layer_indices = respiration_indices,
            .primary_length_m_by_layer = primary_lengths,
        };
    }

    pub fn deinit(self: *AxisWorkspace) void {
        self.allocator.free(self.primary_length_m_by_layer);
        self.allocator.free(self.primary_respiration_root_layer_indices);
        self.allocator.free(self.primary_respiration_allocation_fractions);
        self.allocator.free(self.primary_processed_by_domain);
        self.allocator.free(self.secondary_active);
        self.allocator.free(self.primary_active);
        self.allocator.free(self.primary_deficit_active);
        self.allocator.free(self.primary_deficit_absorption);
        self.allocator.free(self.secondary_senescence);
        self.allocator.free(self.primary_senescence);
        self.allocator.free(self.secondary_metabolism);
        self.allocator.free(self.primary_metabolism);
        self.allocator.free(self.secondary_sink_fractions);
        self.allocator.free(self.primary_sink_fractions);
        self.allocator.free(self.sink_strengths);
        self.* = undefined;
    }

    pub fn resetAxes(self: *AxisWorkspace, active_axis_count: usize) !void {
        if (active_axis_count > self.axis_capacity) return error.RootMetabolismWorkspaceCapacityExceeded;
        @memset(self.sink_strengths[0..active_axis_count], .{ .primary_m = 0, .secondary_m = 0 });
        @memset(self.primary_sink_fractions[0..active_axis_count], 0);
        @memset(self.secondary_sink_fractions[0..active_axis_count], 0);
        @memset(self.primary_metabolism[0..active_axis_count], std.mem.zeroes(SecondaryRootResult));
        @memset(self.secondary_metabolism[0..active_axis_count], std.mem.zeroes(SecondaryRootResult));
        @memset(self.primary_senescence[0..active_axis_count], std.mem.zeroes(SecondaryRootSenescence));
        @memset(self.secondary_senescence[0..active_axis_count], std.mem.zeroes(SecondaryRootSenescence));
        @memset(self.primary_deficit_absorption[0..active_axis_count], std.mem.zeroes(SecondaryRootDeficitAbsorption));
        @memset(self.primary_deficit_active[0..active_axis_count], false);
        @memset(self.primary_active[0..active_axis_count], false);
        @memset(self.secondary_active[0..active_axis_count], false);
    }

    pub fn beginPlantHour(self: *AxisWorkspace, active_axis_count: usize) !void {
        if (active_axis_count > self.axis_capacity) return error.RootMetabolismWorkspaceCapacityExceeded;
        @memset(self.primary_processed_by_domain[0 .. active_axis_count * root_domain_count], false);
    }

    pub fn primaryWasProcessed(self: AxisWorkspace, domain: usize, axis: usize) !bool {
        if (domain >= root_domain_count or axis >= self.axis_capacity) return error.RootMetabolismWorkspaceCapacityExceeded;
        return self.primary_processed_by_domain[axis * root_domain_count + domain];
    }

    pub fn markPrimaryProcessed(self: *AxisWorkspace, domain: usize, axis: usize) !void {
        if (domain >= root_domain_count or axis >= self.axis_capacity) return error.RootMetabolismWorkspaceCapacityExceeded;
        self.primary_processed_by_domain[axis * root_domain_count + domain] = true;
    }
};

pub const GridWorkspace = struct {
    allocator: std.mem.Allocator,
    per_cell: []AxisWorkspace,

    pub fn init(allocator: std.mem.Allocator, cell_count: usize, axis_capacity: usize, soil_layer_count: usize) !GridWorkspace {
        if (cell_count == 0) return error.InvalidRootMetabolismWorkspaceDimensions;
        const cells = try allocator.alloc(AxisWorkspace, cell_count);
        errdefer allocator.free(cells);
        var initialized: usize = 0;
        errdefer for (cells[0..initialized]) |*cell| cell.deinit();
        for (cells) |*cell| {
            cell.* = try AxisWorkspace.init(allocator, axis_capacity, soil_layer_count);
            initialized += 1;
        }
        return .{ .allocator = allocator, .per_cell = cells };
    }

    pub fn deinit(self: *GridWorkspace) void {
        for (self.per_cell) |*cell| cell.deinit();
        self.allocator.free(self.per_cell);
        self.* = undefined;
    }
};

/// GROSUB RTSK1 and RTSK2 for one runtime root axis. Primary-domain
/// secondary roots retain the source harmonic series resistance.
pub fn rootAxisSinkStrength(parameters: SecondaryRootParameters, inputs: RootAxisSinkInputs) !RootAxisSinkStrength {
    try parameters.validate();
    if (inputs.root_profile_type > 3) return error.InvalidRootProfileType;
    inline for (@typeInfo(RootAxisSinkInputs).@"struct".fields) |field| {
        if (field.type == bool or field.type == u8) continue;
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidRootAxisSinkInput;
    }
    const profile_multiplier = switch (inputs.root_profile_type) {
        0 => parameters.shallow_primary_root_sink_multiplier,
        1 => parameters.intermediate_primary_root_sink_multiplier,
        2 => parameters.deep_primary_root_sink_multiplier,
        3 => parameters.deeper_primary_root_sink_multiplier,
        else => unreachable,
    };
    const primary = if (inputs.primary_root_depth_from_canopy_m > 0)
        profile_multiplier * inputs.primary_axis_count_multiplier * inputs.primary_root_radius_m * inputs.primary_root_radius_m / inputs.primary_root_depth_from_canopy_m
    else
        0;
    const secondary_parallel = if (inputs.average_secondary_root_length_m > 0)
        inputs.secondary_axis_count * inputs.secondary_root_radius_m * inputs.secondary_root_radius_m / inputs.average_secondary_root_length_m
    else
        0;
    const secondary = if (inputs.primary_biological_domain) blk: {
        if (inputs.secondary_root_depth_from_canopy_m <= 0) break :blk 0;
        const primary_series = inputs.primary_axis_count_multiplier * inputs.primary_root_radius_m * inputs.primary_root_radius_m / inputs.secondary_root_depth_from_canopy_m;
        break :blk if (primary_series + secondary_parallel > 0)
            primary_series * secondary_parallel / (primary_series + secondary_parallel)
        else
            0;
    } else secondary_parallel;
    inline for (.{ primary, secondary }) |value| if (!std.math.isFinite(value) or value < 0) return error.NonFiniteRootAxisSinkStrength;
    return .{ .primary_m = primary, .secondary_m = secondary };
}

/// Direct GROSUB 5888--5943 geometry comparator. It retains the primary-tip
/// layer gate and the clipped secondary-root midpoint depth that are not
/// represented by the simplified production operand contract.
pub fn sourceOrderRootAxisSinkStrength(parameters: SecondaryRootParameters, inputs: SourceOrderRootAxisSinkInputs) !RootAxisSinkStrength {
    try parameters.validate();
    if (inputs.root_profile_type > 3) return error.InvalidRootProfileType;
    inline for (@typeInfo(SourceOrderRootAxisSinkInputs).@"struct".fields) |field| {
        if (field.type == bool or field.type == u8) continue;
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidRootAxisSinkInput;
    }
    if (inputs.average_secondary_root_length_m <= 0) return error.InvalidRootAxisSinkInput;
    const profile_multiplier = switch (inputs.root_profile_type) {
        0 => parameters.shallow_primary_root_sink_multiplier,
        1 => parameters.intermediate_primary_root_sink_multiplier,
        2 => parameters.deep_primary_root_sink_multiplier,
        3 => parameters.deeper_primary_root_sink_multiplier,
        else => unreachable,
    };
    const layer_bottom_depth_m = inputs.layer_top_depth_m + inputs.layer_thickness_m;
    const primary_depth_from_canopy_m = inputs.primary_root_depth_from_surface_m + inputs.canopy_height_m;
    const tip_in_layer = inputs.primary_root_depth_from_surface_m > inputs.layer_top_depth_m and
        inputs.primary_root_depth_from_surface_m <= layer_bottom_depth_m;
    const primary = if (inputs.primary_biological_domain and tip_in_layer and primary_depth_from_canopy_m > 0)
        profile_multiplier * inputs.primary_axis_count_multiplier *
            inputs.primary_root_radius_m * inputs.primary_root_radius_m / primary_depth_from_canopy_m
    else
        0;
    var rooted_length_m = @max(0, inputs.primary_root_depth_from_surface_m -
        inputs.layer_top_depth_m - inputs.secondary_root_origin_offset_m);
    rooted_length_m = @max(0, @min(inputs.layer_thickness_m, rooted_length_m) -
        @max(0, inputs.seeding_depth_m - inputs.layer_top_depth_m - inputs.hypocotyledon_height_m));
    const secondary_depth_from_canopy_m = @max(inputs.seeding_depth_m, inputs.layer_top_depth_m) +
        0.5 * rooted_length_m + inputs.canopy_height_m;
    const secondary_parallel = inputs.secondary_axis_count *
        inputs.secondary_root_radius_m * inputs.secondary_root_radius_m /
        inputs.average_secondary_root_length_m;
    const secondary = if (inputs.primary_biological_domain) blk: {
        if (secondary_depth_from_canopy_m <= 0) break :blk 0;
        const primary_series = inputs.primary_axis_count_multiplier *
            inputs.primary_root_radius_m * inputs.primary_root_radius_m /
            secondary_depth_from_canopy_m;
        break :blk if (primary_series + secondary_parallel > inputs.negligible_sink_m)
            primary_series * secondary_parallel / (primary_series + secondary_parallel)
        else
            0;
    } else secondary_parallel;
    inline for (.{ primary, secondary }) |value| if (!std.math.isFinite(value) or value < 0) return error.NonFiniteRootAxisSinkStrength;
    return .{ .primary_m = primary, .secondary_m = secondary };
}

/// GROSUB FRTN normalization against RLNT. Caller-owned output slices avoid
/// hourly allocation and scale to any runtime axis count.
pub fn normalizeRootAxisSinkFractions(
    strengths: []const RootAxisSinkStrength,
    primary_fractions: []f64,
    secondary_fractions: []f64,
    negligible_sink_m: f64,
) !f64 {
    if (strengths.len != primary_fractions.len or strengths.len != secondary_fractions.len) return error.RootAxisSinkDimensionMismatch;
    if (!std.math.isFinite(negligible_sink_m) or negligible_sink_m < 0) return error.InvalidRootAxisSinkThreshold;
    var total: f64 = 0;
    for (strengths) |strength| {
        inline for (@typeInfo(RootAxisSinkStrength).@"struct".fields) |field| {
            const value = @field(strength, field.name);
            if (!std.math.isFinite(value) or value < 0) return error.InvalidRootAxisSinkStrength;
            total += value;
        }
    }
    if (!std.math.isFinite(total)) return error.NonFiniteRootAxisSinkTotal;
    for (strengths, 0..) |strength, axis| {
        // This per-branch fallback of 1.0 is intentionally source-exact.
        primary_fractions[axis] = if (total > negligible_sink_m) strength.primary_m / total else 1;
        secondary_fractions[axis] = if (total > negligible_sink_m) strength.secondary_m / total else 1;
    }
    return total;
}

/// Exact GROSUB secondary-root CNPG through CNRDM/CNRDA equations. The
/// caller supplies the source model's dimensionless temperature, water,
/// acidity, oxygen, and feedback responses so this kernel remains reusable.
fn rootMetabolism(parameters: SecondaryRootParameters, inputs: SecondaryRootInputs, cap_substrate_to_current_maintenance: bool) !SecondaryRootResult {
    try parameters.validate();
    inline for (@typeInfo(SecondaryRootInputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootMetabolismInput;
    }
    if (inputs.root_growth_yield_g_c_per_g_c <= 0 or inputs.root_growth_yield_g_c_per_g_c >= 1) return error.InvalidSecondaryRootMetabolismInput;
    if (inputs.root_nitrogen_to_carbon_ratio_g_n_per_g_c == 0 or inputs.root_phosphorus_to_carbon_ratio_g_p_per_g_c == 0) return error.InvalidSecondaryRootMetabolismInput;

    const nutrient_feedback = if (inputs.mobile_carbon_g_c > 0)
        @min(
            inputs.nonstructural_nitrogen_g_n /
                (inputs.nonstructural_nitrogen_g_n + inputs.mobile_carbon_g_c * parameters.nitrogen_feedback_half_saturation_g_n_per_g_c),
            inputs.nonstructural_phosphorus_g_p /
                (inputs.nonstructural_phosphorus_g_p + inputs.mobile_carbon_g_c * parameters.phosphorus_feedback_half_saturation_g_p_per_g_c),
        )
    else
        1.0;
    const mobile_carbon_concentration_g_c_per_g_c = if (inputs.root_carbon_g_c > 0)
        inputs.mobile_carbon_g_c / inputs.root_carbon_g_c
    else
        1;
    var substrate_respiration_oxygen_unlimited = @max(0, parameters.maximum_substrate_respiration_fraction_per_h *
        inputs.active_root_fraction * inputs.mobile_carbon_g_c * inputs.substrate_temperature_response *
        nutrient_feedback * inputs.substrate_feedback * inputs.substrate_water_response * inputs.biological_timestep_h) *
        mobile_carbon_concentration_g_c_per_g_c /
        (mobile_carbon_concentration_g_c_per_g_c + parameters.substrate_respiration_half_saturation_g_c_per_g_c);
    const maintenance_respiration = @max(0, parameters.maintenance_respiration_g_c_per_g_n_h *
        inputs.root_nitrogen_g_n * inputs.maintenance_temperature_response * inputs.acidity_response *
        inputs.biological_timestep_h * inputs.maintenance_water_response);
    // GROSUB places the bottom-tip cap before assigning RMNCR, which reuses a
    // stale loop value. Use the current axis maintenance demand: this is the
    // intended AMIN1 relation and removes traversal-order dependence.
    if (cap_substrate_to_current_maintenance) substrate_respiration_oxygen_unlimited = @min(substrate_respiration_oxygen_unlimited, maintenance_respiration);
    const substrate_respiration_actual = substrate_respiration_oxygen_unlimited * inputs.oxygen_limitation;
    const growth_energy_oxygen_unlimited = @max(0, substrate_respiration_oxygen_unlimited - maintenance_respiration);
    const growth_energy_actual = @max(0, substrate_respiration_actual - maintenance_respiration);
    const root_respiration_fraction = 1.0 - inputs.root_growth_yield_g_c_per_g_c;
    const nutrient_limited_growth_respiration = try nutrientLimitedRootGrowthRespiration(
        inputs.nonstructural_nitrogen_g_n,
        inputs.nonstructural_phosphorus_g_p,
        inputs.active_root_fraction,
        root_respiration_fraction,
        inputs.root_growth_yield_g_c_per_g_c,
        inputs.root_nitrogen_to_carbon_ratio_g_n_per_g_c,
        inputs.root_phosphorus_to_carbon_ratio_g_p_per_g_c,
    );
    const growth_respiration_oxygen_unlimited = @min(growth_energy_oxygen_unlimited, nutrient_limited_growth_respiration);
    const growth_respiration_actual = @min(growth_energy_actual, nutrient_limited_growth_respiration * inputs.oxygen_limitation);
    const root_growth_oxygen_unlimited = growth_respiration_oxygen_unlimited / root_respiration_fraction * inputs.root_growth_yield_g_c_per_g_c;
    const root_growth_actual = growth_respiration_actual / root_respiration_fraction * inputs.root_growth_yield_g_c_per_g_c;
    const nitrogen_growth_demand = @max(0, root_growth_oxygen_unlimited * inputs.root_nitrogen_to_carbon_ratio_g_n_per_g_c);
    const nitrogen_growth_actual = @max(0, @min(inputs.active_root_fraction * inputs.nonstructural_nitrogen_g_n, root_growth_actual * inputs.root_nitrogen_to_carbon_ratio_g_n_per_g_c));
    const phosphorus_growth_actual = @max(0, @min(inputs.active_root_fraction * inputs.nonstructural_phosphorus_g_p, root_growth_actual * inputs.root_phosphorus_to_carbon_ratio_g_p_per_g_c));
    const result: SecondaryRootResult = .{
        .nutrient_feedback = nutrient_feedback,
        .substrate_respiration_oxygen_unlimited_g_c_per_h = substrate_respiration_oxygen_unlimited,
        .maintenance_respiration_g_c_per_h = maintenance_respiration,
        .substrate_respiration_actual_g_c_per_h = substrate_respiration_actual,
        .growth_respiration_oxygen_unlimited_g_c_per_h = growth_respiration_oxygen_unlimited,
        .growth_respiration_actual_g_c_per_h = growth_respiration_actual,
        .growth_and_respiration_carbon_oxygen_unlimited_g_c_per_h = growth_respiration_oxygen_unlimited / root_respiration_fraction,
        .growth_and_respiration_carbon_actual_g_c_per_h = growth_respiration_actual / root_respiration_fraction,
        .root_growth_oxygen_unlimited_g_c_per_h = root_growth_oxygen_unlimited,
        .root_growth_actual_g_c_per_h = root_growth_actual,
        .nitrogen_growth_demand_g_n_per_h = nitrogen_growth_demand,
        .nitrogen_growth_actual_g_n_per_h = nitrogen_growth_actual,
        .phosphorus_growth_actual_g_p_per_h = phosphorus_growth_actual,
        .nitrogen_assimilation_respiration_oxygen_unlimited_g_c_per_h = parameters.nitrogen_assimilation_respiration_g_c_per_g_n * nitrogen_growth_demand,
        .nitrogen_assimilation_respiration_actual_g_c_per_h = parameters.nitrogen_assimilation_respiration_g_c_per_g_n * nitrogen_growth_actual,
    };
    inline for (@typeInfo(SecondaryRootResult).@"struct".fields) |field| if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteSecondaryRootMetabolism;
    return result;
}

pub fn secondaryRootMetabolism(parameters: SecondaryRootParameters, inputs: SecondaryRootInputs) !SecondaryRootResult {
    return rootMetabolism(parameters, inputs, false);
}

/// GROSUB primary-root RCO2RM through CNRDM/CNRDA. Primary and secondary
/// axes share the biochemical equations; only a primary tip at/below the
/// profile bottom receives the source respiration cap.
pub fn primaryRootMetabolism(parameters: SecondaryRootParameters, inputs: PrimaryRootInputs) !SecondaryRootResult {
    return rootMetabolism(parameters, inputs.shared, inputs.primary_tip_at_or_below_profile_bottom);
}

pub const RecyclingFractions = struct {
    carbon: f64,
    nitrogen: f64,
    phosphorus: f64,
};

pub const RootWoodComposition = struct {
    carbon_fraction: [2]f64,
    nitrogen_fraction: [2]f64,
    phosphorus_fraction: [2]f64,
    growth_nitrogen_to_carbon_g_n_per_g_c: f64,
    growth_phosphorus_to_carbon_g_p_per_g_c: f64,
};

/// GROSUB FWODR/FWODRN/FWODRP and CNRTW/CPRTW. Array index 0 is woody
/// and index 1 is nonwoody, matching the source equations.
pub fn rootWoodComposition(
    woody_growth_enabled: bool,
    deep_root_profile: bool,
    stalk_carbon_g_c: f64,
    sapwood_carbon_g_c: f64,
    stalk_nitrogen_to_carbon_g_n_per_g_c: f64,
    root_nitrogen_to_carbon_g_n_per_g_c: f64,
    stalk_phosphorus_to_carbon_g_p_per_g_c: f64,
    root_phosphorus_to_carbon_g_p_per_g_c: f64,
    structural_presence_threshold_g_c: f64,
    nonwoody_root_fraction_exponent: f64,
) !RootWoodComposition {
    inline for (.{ stalk_carbon_g_c, sapwood_carbon_g_c, stalk_nitrogen_to_carbon_g_n_per_g_c, root_nitrogen_to_carbon_g_n_per_g_c, stalk_phosphorus_to_carbon_g_p_per_g_c, root_phosphorus_to_carbon_g_p_per_g_c, structural_presence_threshold_g_c, nonwoody_root_fraction_exponent }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootWoodCompositionInput;
    if (sapwood_carbon_g_c > stalk_carbon_g_c + 1.0e-12) return error.RootSapwoodExceedsStalkCarbon;
    const nonwoody = if (!woody_growth_enabled or !deep_root_profile or stalk_carbon_g_c <= structural_presence_threshold_g_c)
        1.0
    else
        std.math.pow(f64, sapwood_carbon_g_c / stalk_carbon_g_c, nonwoody_root_fraction_exponent);
    if (!std.math.isFinite(nonwoody) or nonwoody < 0 or nonwoody > 1) return error.NonFiniteRootWoodComposition;
    const fractions = [2]f64{ 1 - nonwoody, nonwoody };
    return .{
        .carbon_fraction = fractions,
        .nitrogen_fraction = fractions,
        .phosphorus_fraction = fractions,
        .growth_nitrogen_to_carbon_g_n_per_g_c = fractions[0] * stalk_nitrogen_to_carbon_g_n_per_g_c + fractions[1] * root_nitrogen_to_carbon_g_n_per_g_c,
        .growth_phosphorus_to_carbon_g_p_per_g_c = fractions[0] * stalk_phosphorus_to_carbon_g_p_per_g_c + fractions[1] * root_phosphorus_to_carbon_g_p_per_g_c,
    };
}

pub fn secondaryRootRecyclingFractions(
    emerged: bool,
    mobile_carbon_concentration_g_c_per_g_c: f64,
    mobile_nitrogen_concentration_g_n_per_g_c: f64,
    mobile_phosphorus_concentration_g_p_per_g_c: f64,
    parameters: SecondaryRootParameters,
) !RecyclingFractions {
    try parameters.validate();
    inline for (.{ mobile_carbon_concentration_g_c_per_g_c, mobile_nitrogen_concentration_g_n_per_g_c, mobile_phosphorus_concentration_g_p_per_g_c }) |value| {
        if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootRecyclingInput;
    }
    var carbon_constraint: f64 = 1;
    var nitrogen_constraint: f64 = 0;
    var phosphorus_constraint: f64 = 0;
    if (emerged and mobile_carbon_concentration_g_c_per_g_c > 0) {
        carbon_constraint = std.math.clamp(@min(
            mobile_nitrogen_concentration_g_n_per_g_c / (mobile_nitrogen_concentration_g_n_per_g_c + mobile_carbon_concentration_g_c_per_g_c * parameters.nitrogen_feedback_half_saturation_g_n_per_g_c),
            mobile_phosphorus_concentration_g_p_per_g_c / (mobile_phosphorus_concentration_g_p_per_g_c + mobile_carbon_concentration_g_c_per_g_c * parameters.phosphorus_feedback_half_saturation_g_p_per_g_c),
        ), 0, 1);
        nitrogen_constraint = std.math.clamp(
            mobile_carbon_concentration_g_c_per_g_c / (mobile_carbon_concentration_g_c_per_g_c + mobile_nitrogen_concentration_g_n_per_g_c / parameters.nitrogen_feedback_half_saturation_g_n_per_g_c),
            0,
            1,
        );
        phosphorus_constraint = std.math.clamp(
            mobile_carbon_concentration_g_c_per_g_c / (mobile_carbon_concentration_g_c_per_g_c + mobile_phosphorus_concentration_g_p_per_g_c / parameters.phosphorus_feedback_half_saturation_g_p_per_g_c),
            0,
            1,
        );
    }
    const result: RecyclingFractions = .{
        .carbon = parameters.minimum_carbon_recycling_fraction + carbon_constraint * parameters.responsive_carbon_recycling_fraction,
        .nitrogen = nitrogen_constraint * parameters.maximum_nitrogen_recycling_fraction,
        .phosphorus = phosphorus_constraint * parameters.maximum_phosphorus_recycling_fraction,
    };
    inline for (@typeInfo(RecyclingFractions).@"struct".fields) |field| {
        const value = @field(result, field.name);
        if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidSecondaryRootRecyclingFraction;
    }
    return result;
}

pub const SecondaryRootSenescenceInputs = struct {
    oxygen_unlimited_substrate_minus_maintenance_g_c_per_h: f64,
    actual_substrate_minus_maintenance_g_c_per_h: f64,
    root_carbon_g_c: f64,
    root_nitrogen_g_n: f64,
    root_phosphorus_g_p: f64,
    oxygen_limitation: f64,
    phenological_remobilization_enabled: bool,
    root_remobilization_enabled: bool,
    storage_exchange_fraction_per_h: f64,
    remobilization_elapsed_h: f64,
    full_senescence_h: f64,
    biological_timestep_h: f64,
    structural_presence_threshold_g_c: f64,
};

pub const SecondaryRootSenescence = struct {
    respiration_oxygen_unlimited_g_c_per_h: f64,
    respiration_actual_g_c_per_h: f64,
    phenological_senescence_g_c_per_h: f64,
    senesced_fraction: f64,
    recyclable_carbon_g_c: f64,
    recyclable_nitrogen_g_n: f64,
    recyclable_phosphorus_g_p: f64,
};

/// Exact GROSUB SNCRM/SNCR/SNCZ and RCCR/RCZR/RCPR/FSNC2 block.
pub fn secondaryRootSenescence(inputs: SecondaryRootSenescenceInputs, recycling: RecyclingFractions) !SecondaryRootSenescence {
    inline for (@typeInfo(SecondaryRootSenescenceInputs).@"struct".fields) |field| {
        if (field.type == bool) continue;
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteSecondaryRootSenescenceInput;
    }
    inline for (.{ inputs.root_carbon_g_c, inputs.root_nitrogen_g_n, inputs.root_phosphorus_g_p, inputs.oxygen_limitation, inputs.storage_exchange_fraction_per_h, inputs.remobilization_elapsed_h, inputs.full_senescence_h, inputs.biological_timestep_h, inputs.structural_presence_threshold_g_c }) |value| if (value < 0) return error.InvalidSecondaryRootSenescenceInput;
    if (inputs.full_senescence_h == 0 or inputs.biological_timestep_h == 0 or inputs.oxygen_limitation > 1) return error.InvalidSecondaryRootSenescenceInput;
    inline for (@typeInfo(RecyclingFractions).@"struct".fields) |field| {
        const value = @field(recycling, field.name);
        if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidSecondaryRootRecyclingFraction;
    }

    const recyclable_carbon = inputs.root_carbon_g_c * recycling.carbon;
    const recyclable_nitrogen = inputs.root_nitrogen_g_n * (recycling.nitrogen + (1 - recycling.nitrogen) * recycling.carbon);
    const recyclable_phosphorus = inputs.root_phosphorus_g_p * (recycling.phosphorus + (1 - recycling.phosphorus) * recycling.carbon);
    const oxygen_unlimited_deficit = @max(0, -inputs.oxygen_unlimited_substrate_minus_maintenance_g_c_per_h);
    const actual_deficit = @max(0, -inputs.actual_substrate_minus_maintenance_g_c_per_h);
    const oxygen_unlimited_respiration = @min(oxygen_unlimited_deficit, recyclable_carbon);
    var actual_respiration = if (actual_deficit < recyclable_carbon)
        actual_deficit
    else
        recyclable_carbon * inputs.oxygen_limitation;
    const phenological = if (inputs.phenological_remobilization_enabled and inputs.root_remobilization_enabled)
        inputs.storage_exchange_fraction_per_h * inputs.root_carbon_g_c *
            @min(1, inputs.remobilization_elapsed_h / inputs.full_senescence_h) * inputs.biological_timestep_h
    else
        0;
    actual_respiration += phenological;
    const senesced_fraction = if (actual_respiration > 0 and inputs.root_carbon_g_c > inputs.structural_presence_threshold_g_c)
        (if (recyclable_carbon > inputs.structural_presence_threshold_g_c)
            std.math.clamp(actual_respiration / recyclable_carbon, 0, 1)
        else
            1)
    else
        0;
    return .{
        .respiration_oxygen_unlimited_g_c_per_h = oxygen_unlimited_respiration,
        .respiration_actual_g_c_per_h = actual_respiration,
        .phenological_senescence_g_c_per_h = phenological,
        .senesced_fraction = senesced_fraction,
        .recyclable_carbon_g_c = if (senesced_fraction > 0) recyclable_carbon else 0,
        .recyclable_nitrogen_g_n = if (senesced_fraction > 0) recyclable_nitrogen else 0,
        .recyclable_phosphorus_g_p = if (senesced_fraction > 0) recyclable_phosphorus else 0,
    };
}

/// GROSUB 6684--6718 primary-root senescence has no SNCZ phenological term.
/// Keep the shared deficit equations while explicitly disabling the
/// secondary-root-only remobilization branch.
pub fn primaryRootSenescence(inputs: SecondaryRootSenescenceInputs, recycling: RecyclingFractions) !SecondaryRootSenescence {
    var primary_inputs = inputs;
    primary_inputs.phenological_remobilization_enabled = false;
    primary_inputs.root_remobilization_enabled = false;
    return secondaryRootSenescence(primary_inputs, recycling);
}

pub const RootLitterFractions = struct {
    woody_carbon: [4]f64,
    woody_nitrogen: [4]f64,
    woody_phosphorus: [4]f64,
    nonwoody_carbon: [4]f64,
    nonwoody_nitrogen: [4]f64,
    nonwoody_phosphorus: [4]f64,
};

pub const RootLitter = struct {
    woody_carbon_g_c: [4]f64,
    woody_nitrogen_g_n: [4]f64,
    woody_phosphorus_g_p: [4]f64,
    nonwoody_carbon_g_c: [4]f64,
    nonwoody_nitrogen_g_n: [4]f64,
    nonwoody_phosphorus_g_p: [4]f64,
};

pub const MycorrhizalLossState = struct {
    structural_carbon_g_c: f64,
    structural_nitrogen_g_n: f64,
    structural_phosphorus_g_p: f64,
    length_m: f64,
    mobile_carbon_g_c: f64,
    mobile_nitrogen_g_n: f64,
    mobile_phosphorus_g_p: f64,
};

pub const MycorrhizalLossResult = struct {
    remaining: MycorrhizalLossState,
    litter: RootLitter,
    structural_loss_fraction: f64,
    mobile_loss_fraction: f64,
};

/// GROSUB concurrent mycorrhizal loss when negative primary-root growth is
/// absorbed by host secondary roots. Structural and mobile loss fractions use
/// the source model's distinct secondary-root and total-active-root C bases.
pub fn mycorrhizalLossWithSecondaryRoots(
    negative_primary_growth_g_c: f64,
    host_secondary_carbon_g_c: f64,
    host_active_root_carbon_g_c: f64,
    presence_threshold_g_c: f64,
    state: MycorrhizalLossState,
    woody_fraction: [3][2]f64,
    kinetics: RootLitterFractions,
) !MycorrhizalLossResult {
    inline for (.{ negative_primary_growth_g_c, host_secondary_carbon_g_c, host_active_root_carbon_g_c, presence_threshold_g_c }) |value|
        if (!std.math.isFinite(value)) return error.NonFiniteMycorrhizalLossInput;
    if (host_secondary_carbon_g_c < 0 or host_active_root_carbon_g_c < 0 or presence_threshold_g_c < 0)
        return error.InvalidMycorrhizalLossInput;
    inline for (@typeInfo(MycorrhizalLossState).@"struct".fields) |field| {
        const value = @field(state, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidMycorrhizalLossInput;
    }
    for (woody_fraction) |fractions| for (fractions) |value|
        if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidMycorrhizalLossInput;
    inline for (@typeInfo(RootLitterFractions).@"struct".fields) |field| for (@field(kinetics, field.name)) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidMycorrhizalLossInput;

    const deficit_g_c = @max(0, -negative_primary_growth_g_c);
    const structural_fraction = if (deficit_g_c == 0) 0 else if (host_secondary_carbon_g_c > presence_threshold_g_c)
        @min(1, deficit_g_c / host_secondary_carbon_g_c)
    else
        1;
    const mobile_fraction = if (deficit_g_c == 0) 0 else if (host_active_root_carbon_g_c > presence_threshold_g_c)
        @min(1, deficit_g_c / host_active_root_carbon_g_c)
    else
        1;

    var litter = std.mem.zeroes(RootLitter);
    for (0..4) |kinetic| {
        litter.woody_carbon_g_c[kinetic] = kinetics.woody_carbon[kinetic] * structural_fraction * state.structural_carbon_g_c * woody_fraction[0][0];
        litter.woody_nitrogen_g_n[kinetic] = kinetics.woody_nitrogen[kinetic] * structural_fraction * state.structural_nitrogen_g_n * woody_fraction[1][0];
        litter.woody_phosphorus_g_p[kinetic] = kinetics.woody_phosphorus[kinetic] * structural_fraction * state.structural_phosphorus_g_p * woody_fraction[2][0];
        litter.nonwoody_carbon_g_c[kinetic] = kinetics.nonwoody_carbon[kinetic] *
            (structural_fraction * state.structural_carbon_g_c * woody_fraction[0][1] + mobile_fraction * state.mobile_carbon_g_c);
        litter.nonwoody_nitrogen_g_n[kinetic] = kinetics.nonwoody_nitrogen[kinetic] *
            (structural_fraction * state.structural_nitrogen_g_n * woody_fraction[1][1] + mobile_fraction * state.mobile_nitrogen_g_n);
        litter.nonwoody_phosphorus_g_p[kinetic] = kinetics.nonwoody_phosphorus[kinetic] *
            (structural_fraction * state.structural_phosphorus_g_p * woody_fraction[2][1] + mobile_fraction * state.mobile_phosphorus_g_p);
    }
    const structural_retained = 1 - structural_fraction;
    const mobile_retained = 1 - mobile_fraction;
    const remaining: MycorrhizalLossState = .{
        .structural_carbon_g_c = state.structural_carbon_g_c * structural_retained,
        .structural_nitrogen_g_n = state.structural_nitrogen_g_n * structural_retained,
        .structural_phosphorus_g_p = state.structural_phosphorus_g_p * structural_retained,
        .length_m = state.length_m * structural_retained,
        .mobile_carbon_g_c = state.mobile_carbon_g_c * mobile_retained,
        .mobile_nitrogen_g_n = state.mobile_nitrogen_g_n * mobile_retained,
        .mobile_phosphorus_g_p = state.mobile_phosphorus_g_p * mobile_retained,
    };
    return .{
        .remaining = remaining,
        .litter = litter,
        .structural_loss_fraction = structural_fraction,
        .mobile_loss_fraction = mobile_fraction,
    };
}

pub const LayerPairRootLitter = struct {
    current: RootLitter = std.mem.zeroes(RootLitter),
    upper: RootLitter = std.mem.zeroes(RootLitter),
};

fn addRootLitter(total: *RootLitter, addition: RootLitter) !void {
    inline for (@typeInfo(RootLitter).@"struct".fields) |field| {
        for (&@field(total.*, field.name), @field(addition, field.name)) |*destination, value| {
            destination.* += value;
            if (!std.math.isFinite(destination.*) or destination.* < 0) return error.NonFiniteMycorrhizalLoss;
        }
    }
}

/// Applies all source-order GROSUB mycorrhizal losses associated with the
/// staged host-root deficit plans for one primary-root tip layer.
pub fn commitMycorrhizalLossWithSecondaryRoots(
    roots: *RootState,
    plant: usize,
    layer: usize,
    workspace: AxisWorkspace,
    active_axis_count: usize,
    host_active_root_carbon_g_c_by_layer: [2]f64,
    presence_threshold_g_c: f64,
    woody_fraction: [3][2]f64,
    kinetics: RootLitterFractions,
) !LayerPairRootLitter {
    if (active_axis_count > workspace.axis_capacity or active_axis_count > roots.root_axis_count)
        return error.RootMetabolismWorkspaceCapacityExceeded;
    if (layer >= roots.soil_layer_count) return error.RootMetabolismLayerOutOfBounds;
    var litter: LayerPairRootLitter = .{};
    const current_root = try roots.layerIndex(plant, 1, layer);
    const upper_root = if (layer > 0) try roots.layerIndex(plant, 1, layer - 1) else current_root;
    var virtual_mobile = [2][3]f64{
        .{ roots.mobile_carbon_g[current_root], roots.mobile_nitrogen_g[current_root], roots.mobile_phosphorus_g[current_root] },
        .{ roots.mobile_carbon_g[upper_root], roots.mobile_nitrogen_g[upper_root], roots.mobile_phosphorus_g[upper_root] },
    };

    // Validate and accumulate the complete two-layer transaction before any
    // state is changed. Shared mobile pools are advanced virtually in axis order.
    for (0..active_axis_count) |axis| {
        if (!workspace.primary_deficit_active[axis]) continue;
        const absorption = workspace.primary_deficit_absorption[axis];
        for (0..2) |layer_offset| {
            if (layer_offset == 1 and layer == 0) continue;
            const affected_layer = layer - layer_offset;
            const axis_layer = try roots.layerAxisIndex(plant, 1, affected_layer, axis);
            _ = try roots.layerIndex(plant, 1, affected_layer);
            const entering_deficit = if (layer_offset == 0)
                absorption.current_entering_carbon_deficit_g_c
            else
                absorption.upper_entering_carbon_deficit_g_c;
            const host_remaining = if (layer_offset == 0) absorption.current.carbon_g_c else absorption.upper.carbon_g_c;
            const result = try mycorrhizalLossWithSecondaryRoots(
                -entering_deficit,
                host_remaining,
                host_active_root_carbon_g_c_by_layer[layer_offset],
                presence_threshold_g_c,
                .{
                    .structural_carbon_g_c = roots.axis_secondary_carbon_g[axis_layer],
                    .structural_nitrogen_g_n = roots.axis_secondary_nitrogen_g[axis_layer],
                    .structural_phosphorus_g_p = roots.axis_secondary_phosphorus_g[axis_layer],
                    .length_m = roots.axis_secondary_length_m[axis_layer],
                    .mobile_carbon_g_c = virtual_mobile[layer_offset][0],
                    .mobile_nitrogen_g_n = virtual_mobile[layer_offset][1],
                    .mobile_phosphorus_g_p = virtual_mobile[layer_offset][2],
                },
                woody_fraction,
                kinetics,
            );
            virtual_mobile[layer_offset] = .{
                result.remaining.mobile_carbon_g_c,
                result.remaining.mobile_nitrogen_g_n,
                result.remaining.mobile_phosphorus_g_p,
            };
            try addRootLitter(if (layer_offset == 0) &litter.current else &litter.upper, result.litter);
        }
    }

    for (0..active_axis_count) |axis| {
        if (!workspace.primary_deficit_active[axis]) continue;
        const absorption = workspace.primary_deficit_absorption[axis];
        for (0..2) |layer_offset| {
            if (layer_offset == 1 and layer == 0) continue;
            const affected_layer = layer - layer_offset;
            const axis_layer = roots.layerAxisIndex(plant, 1, affected_layer, axis) catch unreachable;
            const root_layer = roots.layerIndex(plant, 1, affected_layer) catch unreachable;
            const entering_deficit = if (layer_offset == 0)
                absorption.current_entering_carbon_deficit_g_c
            else
                absorption.upper_entering_carbon_deficit_g_c;
            const host_remaining = if (layer_offset == 0) absorption.current.carbon_g_c else absorption.upper.carbon_g_c;
            const result = mycorrhizalLossWithSecondaryRoots(
                -entering_deficit,
                host_remaining,
                host_active_root_carbon_g_c_by_layer[layer_offset],
                presence_threshold_g_c,
                .{
                    .structural_carbon_g_c = roots.axis_secondary_carbon_g[axis_layer],
                    .structural_nitrogen_g_n = roots.axis_secondary_nitrogen_g[axis_layer],
                    .structural_phosphorus_g_p = roots.axis_secondary_phosphorus_g[axis_layer],
                    .length_m = roots.axis_secondary_length_m[axis_layer],
                    .mobile_carbon_g_c = roots.mobile_carbon_g[root_layer],
                    .mobile_nitrogen_g_n = roots.mobile_nitrogen_g[root_layer],
                    .mobile_phosphorus_g_p = roots.mobile_phosphorus_g[root_layer],
                },
                woody_fraction,
                kinetics,
            ) catch unreachable;
            roots.axis_secondary_carbon_g[axis_layer] = result.remaining.structural_carbon_g_c;
            roots.axis_secondary_nitrogen_g[axis_layer] = result.remaining.structural_nitrogen_g_n;
            roots.axis_secondary_phosphorus_g[axis_layer] = result.remaining.structural_phosphorus_g_p;
            roots.axis_secondary_length_m[axis_layer] = result.remaining.length_m;
            roots.mobile_carbon_g[root_layer] = result.remaining.mobile_carbon_g_c;
            roots.mobile_nitrogen_g[root_layer] = result.remaining.mobile_nitrogen_g_n;
            roots.mobile_phosphorus_g[root_layer] = result.remaining.mobile_phosphorus_g_p;
        }
    }
    return litter;
}

test "negative primary growth removes concurrent mycorrhizal structure and mobile pools" {
    const kinetics: RootLitterFractions = .{
        .woody_carbon = .{ 0.1, 0.2, 0.3, 0.4 },
        .woody_nitrogen = .{ 0.1, 0.2, 0.3, 0.4 },
        .woody_phosphorus = .{ 0.1, 0.2, 0.3, 0.4 },
        .nonwoody_carbon = .{ 0.4, 0.3, 0.2, 0.1 },
        .nonwoody_nitrogen = .{ 0.4, 0.3, 0.2, 0.1 },
        .nonwoody_phosphorus = .{ 0.4, 0.3, 0.2, 0.1 },
    };
    const result = try mycorrhizalLossWithSecondaryRoots(
        -2,
        8,
        20,
        1.0e-12,
        .{
            .structural_carbon_g_c = 4,
            .structural_nitrogen_g_n = 2,
            .structural_phosphorus_g_p = 1,
            .length_m = 12,
            .mobile_carbon_g_c = 5,
            .mobile_nitrogen_g_n = 3,
            .mobile_phosphorus_g_p = 2,
        },
        .{ .{ 0.6, 0.4 }, .{ 0.5, 0.5 }, .{ 0.25, 0.75 } },
        kinetics,
    );
    try std.testing.expectApproxEqAbs(0.25, result.structural_loss_fraction, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.1, result.mobile_loss_fraction, 1.0e-12);
    try std.testing.expectApproxEqAbs(3, result.remaining.structural_carbon_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(9, result.remaining.length_m, 1.0e-12);
    try std.testing.expectApproxEqAbs(4.5, result.remaining.mobile_carbon_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.06, result.litter.woody_carbon_g_c[0], 1.0e-12);
    try std.testing.expectApproxEqAbs(0.36, result.litter.nonwoody_carbon_g_c[0], 1.0e-12);

    const unchanged = try mycorrhizalLossWithSecondaryRoots(
        0.2,
        0,
        0,
        1.0e-12,
        result.remaining,
        .{ .{ 0.6, 0.4 }, .{ 0.5, 0.5 }, .{ 0.25, 0.75 } },
        kinetics,
    );
    try std.testing.expectEqual(@as(f64, 0), unchanged.structural_loss_fraction);
    try std.testing.expectEqual(result.remaining, unchanged.remaining);
}

test "live GROSUB mycorrhizal loss follows current then upper host-root deficits" {
    var roots = try RootState.init(std.testing.allocator, 1, 2, 1);
    defer roots.deinit();
    var workspace = try AxisWorkspace.init(std.testing.allocator, 1, 2);
    defer workspace.deinit();
    try workspace.resetAxes(1);
    workspace.primary_deficit_active[0] = true;
    workspace.primary_deficit_absorption[0] = try absorbPrimaryDeficitFromSecondaryRoots(
        3,
        0,
        0,
        .{ .carbon_g_c = 2, .nitrogen_g_n = 0, .phosphorus_g_p = 0, .length_m = 2 },
        .{ .carbon_g_c = 4, .nitrogen_g_n = 0, .phosphorus_g_p = 0, .length_m = 4 },
    );
    const upper_axis = try roots.layerAxisIndex(0, 1, 0, 0);
    const current_axis = try roots.layerAxisIndex(0, 1, 1, 0);
    const upper_root = try roots.layerIndex(0, 1, 0);
    const current_root = try roots.layerIndex(0, 1, 1);
    roots.axis_secondary_carbon_g[upper_axis] = 8;
    roots.axis_secondary_carbon_g[current_axis] = 4;
    roots.axis_secondary_length_m[upper_axis] = 16;
    roots.axis_secondary_length_m[current_axis] = 8;
    roots.mobile_carbon_g[upper_root] = 8;
    roots.mobile_carbon_g[current_root] = 10;
    const unit = [_]f64{1} ** 4;
    const zero = [_]f64{0} ** 4;
    const litter = try commitMycorrhizalLossWithSecondaryRoots(
        &roots,
        0,
        1,
        workspace,
        1,
        .{ 10, 8 },
        1.0e-12,
        .{ .{ 0, 1 }, .{ 0, 1 }, .{ 0, 1 } },
        .{
            .woody_carbon = zero,
            .woody_nitrogen = zero,
            .woody_phosphorus = zero,
            .nonwoody_carbon = unit,
            .nonwoody_nitrogen = unit,
            .nonwoody_phosphorus = unit,
        },
    );
    try std.testing.expectEqual(@as(f64, 0), roots.axis_secondary_carbon_g[current_axis]);
    try std.testing.expectEqual(@as(f64, 0), roots.axis_secondary_length_m[current_axis]);
    try std.testing.expectApproxEqAbs(@as(f64, 16.0 / 3.0), roots.axis_secondary_carbon_g[upper_axis], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 32.0 / 3.0), roots.axis_secondary_length_m[upper_axis], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 7), roots.mobile_carbon_g[current_root], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 7), roots.mobile_carbon_g[upper_root], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 7), litter.current.nonwoody_carbon_g_c[0], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 11.0 / 3.0), litter.upper.nonwoody_carbon_g_c[0], 1.0e-12);
}

/// GROSUB CSNC/ZSNC/PSNC secondary-root allocation across the four kinetic
/// litter fractions. Mg and mg do not appear here: all masses are grams.
pub fn secondaryRootLitter(
    senescence: SecondaryRootSenescence,
    root_carbon_g_c: f64,
    root_nitrogen_g_n: f64,
    root_phosphorus_g_p: f64,
    woody_carbon_fraction: [2]f64,
    woody_nitrogen_fraction: [2]f64,
    woody_phosphorus_fraction: [2]f64,
    kinetics: RootLitterFractions,
) !RootLitter {
    inline for (.{ root_carbon_g_c, root_nitrogen_g_n, root_phosphorus_g_p }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootLitterInput;
    inline for (.{ woody_carbon_fraction, woody_nitrogen_fraction, woody_phosphorus_fraction }) |fractions| for (fractions) |value| if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidSecondaryRootLitterInput;
    inline for (@typeInfo(RootLitterFractions).@"struct".fields) |field| for (@field(kinetics, field.name)) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootLitterInput;
    inline for (@typeInfo(SecondaryRootSenescence).@"struct".fields) |field| if (!std.math.isFinite(@field(senescence, field.name)) or @field(senescence, field.name) < 0) return error.InvalidSecondaryRootLitterInput;

    var result: RootLitter = undefined;
    for (0..4) |kinetic| {
        result.woody_carbon_g_c[kinetic] = kinetics.woody_carbon[kinetic] * senescence.senesced_fraction * root_carbon_g_c * woody_carbon_fraction[0];
        result.woody_nitrogen_g_n[kinetic] = kinetics.woody_nitrogen[kinetic] * senescence.senesced_fraction * root_nitrogen_g_n * woody_nitrogen_fraction[0];
        result.woody_phosphorus_g_p[kinetic] = kinetics.woody_phosphorus[kinetic] * senescence.senesced_fraction * root_phosphorus_g_p * woody_phosphorus_fraction[0];
        result.nonwoody_carbon_g_c[kinetic] = kinetics.nonwoody_carbon[kinetic] * senescence.senesced_fraction * (root_carbon_g_c - senescence.recyclable_carbon_g_c) * woody_carbon_fraction[1];
        result.nonwoody_nitrogen_g_n[kinetic] = kinetics.nonwoody_nitrogen[kinetic] * senescence.senesced_fraction * (root_nitrogen_g_n - senescence.recyclable_nitrogen_g_n) * woody_nitrogen_fraction[1];
        result.nonwoody_phosphorus_g_p[kinetic] = kinetics.nonwoody_phosphorus[kinetic] * senescence.senesced_fraction * (root_phosphorus_g_p - senescence.recyclable_phosphorus_g_p) * woody_phosphorus_fraction[1];
    }
    inline for (@typeInfo(RootLitter).@"struct".fields) |field| for (@field(result, field.name)) |value| if (!std.math.isFinite(value) or value < -1.0e-12) return error.NonFiniteSecondaryRootLitter;
    return result;
}

pub const SecondaryRootCommitInputs = struct {
    metabolism: SecondaryRootResult,
    senescence: SecondaryRootSenescence,
    root_specific_length_m_per_g_c: f64,
    root_extension_water_response: f64,
    nonwoody_carbon_fraction: f64,
    nonwoody_nitrogen_fraction: f64,
    nonwoody_phosphorus_fraction: f64,
    protein_carbon_per_nitrogen_g_c_per_g_n: f64,
    protein_carbon_per_phosphorus_g_c_per_g_p: f64,
};

pub const PrimaryRootCommitInputs = struct {
    metabolism: SecondaryRootResult,
    senescence: SecondaryRootSenescence,
    primary_specific_length_m_per_g_c: f64,
    root_extension_water_response: f64,
    nonwoody_carbon_fraction: f64,
    nonwoody_nitrogen_fraction: f64,
    protein_carbon_per_nitrogen_g_c_per_g_n: f64,
    protein_carbon_per_phosphorus_g_c_per_g_p: f64,
};

pub const StagedLayerCommitParameters = struct {
    primary_specific_length_m_per_g_c: f64,
    secondary_specific_length_m_per_g_c: f64,
    plant_population_count: f64,
    seeding_depth_m: f64,
    current_layer_bottom_depth_m: f64,
    next_layer_thickness_m: f64,
    extension_presence_threshold_m: f64,
    root_extension_water_response: f64,
    nonwoody_carbon_fraction: f64,
    nonwoody_nitrogen_fraction: f64,
    nonwoody_phosphorus_fraction: f64,
    protein_carbon_per_nitrogen_g_c_per_g_n: f64,
    protein_carbon_per_phosphorus_g_c_per_g_p: f64,
};

pub const PrimaryRootExtensionPlacement = struct {
    extension_m: f64,
    crosses_into_next_layer: bool,
};

pub const SecondaryRootDeficitLayer = struct {
    carbon_g_c: f64,
    nitrogen_g_n: f64,
    phosphorus_g_p: f64,
    length_m: f64,
};

pub const SecondaryRootDeficitAbsorption = struct {
    current: SecondaryRootDeficitLayer,
    upper: SecondaryRootDeficitLayer,
    current_entering_carbon_deficit_g_c: f64,
    upper_entering_carbon_deficit_g_c: f64,
    residual_carbon_deficit_g_c: f64,
    residual_nitrogen_deficit_g_n: f64,
    residual_phosphorus_deficit_g_p: f64,
};

/// GROSUB 5105 precursor to primary-tip withdrawal. Negative primary C/N/P
/// growth is absorbed by secondary roots in the tip layer, then the adjacent
/// upper layer. Carbon removal shortens secondary roots proportionally.
pub fn absorbPrimaryDeficitFromSecondaryRoots(
    carbon_deficit_g_c: f64,
    nitrogen_deficit_g_n: f64,
    phosphorus_deficit_g_p: f64,
    current: SecondaryRootDeficitLayer,
    upper: SecondaryRootDeficitLayer,
) !SecondaryRootDeficitAbsorption {
    inline for (.{ carbon_deficit_g_c, nitrogen_deficit_g_n, phosphorus_deficit_g_p }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootDeficit;
    var result: SecondaryRootDeficitAbsorption = .{
        .current = current,
        .upper = upper,
        .current_entering_carbon_deficit_g_c = carbon_deficit_g_c,
        .upper_entering_carbon_deficit_g_c = 0,
        .residual_carbon_deficit_g_c = carbon_deficit_g_c,
        .residual_nitrogen_deficit_g_n = nitrogen_deficit_g_n,
        .residual_phosphorus_deficit_g_p = phosphorus_deficit_g_p,
    };
    inline for (.{ "current", "upper" }) |field_name| {
        if (std.mem.eql(u8, field_name, "upper"))
            result.upper_entering_carbon_deficit_g_c = result.residual_carbon_deficit_g_c;
        var pool = &@field(result, field_name);
        inline for (@typeInfo(SecondaryRootDeficitLayer).@"struct".fields) |field| {
            const value = @field(pool.*, field.name);
            if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootDeficitPool;
        }
        const carbon_removed = @min(result.residual_carbon_deficit_g_c, pool.carbon_g_c);
        if (pool.carbon_g_c > 0) pool.length_m *= 1 - carbon_removed / pool.carbon_g_c;
        pool.carbon_g_c -= carbon_removed;
        result.residual_carbon_deficit_g_c -= carbon_removed;
        const nitrogen_removed = @min(result.residual_nitrogen_deficit_g_n, pool.nitrogen_g_n);
        pool.nitrogen_g_n -= nitrogen_removed;
        result.residual_nitrogen_deficit_g_n -= nitrogen_removed;
        const phosphorus_removed = @min(result.residual_phosphorus_deficit_g_p, pool.phosphorus_g_p);
        pool.phosphorus_g_p -= phosphorus_removed;
        result.residual_phosphorus_deficit_g_p -= phosphorus_removed;
    }
    inline for (@typeInfo(SecondaryRootDeficitAbsorption).@"struct".fields) |field| switch (field.type) {
        f64 => if (!std.math.isFinite(@field(result, field.name)) or @field(result, field.name) < 0) return error.NonFinitePrimaryRootDeficitAbsorption,
        else => {},
    };
    return result;
}

fn stagedSecondaryRootLayer(
    roots: RootState,
    axis_layer: usize,
    active: bool,
    metabolism: SecondaryRootResult,
    senescence: SecondaryRootSenescence,
    specific_length_m_per_g_c: f64,
    extension_water_response: f64,
) !SecondaryRootDeficitLayer {
    const fraction = if (active) senescence.senesced_fraction else 0;
    const result: SecondaryRootDeficitLayer = .{
        .carbon_g_c = roots.axis_secondary_carbon_g[axis_layer] + (if (active) metabolism.root_growth_actual_g_c_per_h else 0) - fraction * roots.axis_secondary_carbon_g[axis_layer],
        .nitrogen_g_n = roots.axis_secondary_nitrogen_g[axis_layer] + (if (active) metabolism.nitrogen_growth_actual_g_n_per_h else 0) - fraction * roots.axis_secondary_nitrogen_g[axis_layer],
        .phosphorus_g_p = roots.axis_secondary_phosphorus_g[axis_layer] + (if (active) metabolism.phosphorus_growth_actual_g_p_per_h else 0) - fraction * roots.axis_secondary_phosphorus_g[axis_layer],
        .length_m = roots.axis_secondary_length_m[axis_layer] +
            (if (active) metabolism.root_growth_actual_g_c_per_h * specific_length_m_per_g_c * extension_water_response else 0) -
            fraction * roots.axis_secondary_length_m[axis_layer],
    };
    inline for (@typeInfo(SecondaryRootDeficitLayer).@"struct".fields) |field| {
        const value = @field(result, field.name);
        if (!std.math.isFinite(value) or value < -1e-12) return error.StagedRootCommitWouldOverdrawPool;
    }
    return result;
}

/// GROSUB GRTLGL cap and binary FGROL/FGROZ placement. If any part of a
/// positive extension crosses the lower boundary, the source assigns all of
/// that hour's primary-root growth to the next layer.
pub fn primaryRootExtensionPlacement(
    requested_extension_m: f64,
    current_depth_m: f64,
    current_layer_bottom_depth_m: f64,
    next_layer_thickness_m: f64,
) !PrimaryRootExtensionPlacement {
    inline for (.{ requested_extension_m, current_depth_m, current_layer_bottom_depth_m, next_layer_thickness_m }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePrimaryRootExtensionPlacement;
    if (current_depth_m < 0 or current_layer_bottom_depth_m < 0 or next_layer_thickness_m < 0) return error.InvalidPrimaryRootExtensionPlacement;
    const extension_m = if (next_layer_thickness_m > 0) @min(next_layer_thickness_m, requested_extension_m) else requested_extension_m;
    const crosses = extension_m > 0 and next_layer_thickness_m > 0 and current_depth_m + extension_m > current_layer_bottom_depth_m;
    if (!std.math.isFinite(extension_m)) return error.NonFinitePrimaryRootExtensionPlacement;
    return .{ .extension_m = extension_m, .crosses_into_next_layer = crosses };
}

/// GROSUB 6986--7007 compatibility placement including the strict
/// population-scaled ZEROP gate before binary current/next-layer routing.
pub fn sourceOrderPrimaryRootExtensionPlacement(
    requested_extension_m: f64,
    current_depth_m: f64,
    current_layer_bottom_depth_m: f64,
    next_layer_thickness_m: f64,
    extension_presence_threshold_m: f64,
) !PrimaryRootExtensionPlacement {
    if (!std.math.isFinite(extension_presence_threshold_m) or extension_presence_threshold_m < 0)
        return error.InvalidPrimaryRootExtensionPlacement;
    const placement = try primaryRootExtensionPlacement(
        requested_extension_m,
        current_depth_m,
        current_layer_bottom_depth_m,
        next_layer_thickness_m,
    );
    if (placement.extension_m <= extension_presence_threshold_m)
        return .{ .extension_m = placement.extension_m, .crosses_into_next_layer = false };
    return placement;
}

/// GROSUB GRTLGL for a primary axis. Negative net C growth retracts the
/// existing rooted depth in proportion to total primary-axis C; gross growth
/// can offset that loss in the same hour.
pub fn primaryRootLengthChange(
    gross_growth_g_c: f64,
    net_growth_g_c: f64,
    total_primary_carbon_g_c: f64,
    current_depth_m: f64,
    seeding_depth_m: f64,
    specific_length_m_per_g_c: f64,
    plant_population_count: f64,
    root_extension_water_response: f64,
) !f64 {
    inline for (.{ gross_growth_g_c, net_growth_g_c, total_primary_carbon_g_c, current_depth_m, seeding_depth_m, specific_length_m_per_g_c, plant_population_count, root_extension_water_response }) |value|
        if (!std.math.isFinite(value)) return error.NonFinitePrimaryRootLengthChange;
    if (gross_growth_g_c < 0 or total_primary_carbon_g_c < 0 or current_depth_m < 0 or seeding_depth_m < 0 or specific_length_m_per_g_c < 0 or plant_population_count <= 0 or root_extension_water_response < 0 or root_extension_water_response > 1) return error.InvalidPrimaryRootLengthChange;
    var change_m = gross_growth_g_c * specific_length_m_per_g_c / plant_population_count * root_extension_water_response;
    if (net_growth_g_c < 0 and total_primary_carbon_g_c > 0)
        change_m += net_growth_g_c * (current_depth_m - seeding_depth_m) / total_primary_carbon_g_c;
    if (!std.math.isFinite(change_m)) return error.NonFinitePrimaryRootLengthChange;
    return change_m;
}

fn mobileChanges(metabolism: SecondaryRootResult, senescence: SecondaryRootSenescence, recovered_c: f64, recovered_n: f64, recovered_p: f64) [3]f64 {
    return .{
        -@min(metabolism.maintenance_respiration_g_c_per_h, metabolism.substrate_respiration_actual_g_c_per_h) -
            metabolism.growth_and_respiration_carbon_actual_g_c_per_h -
            metabolism.nitrogen_assimilation_respiration_actual_g_c_per_h -
            senescence.respiration_actual_g_c_per_h + recovered_c,
        -metabolism.nitrogen_growth_actual_g_n_per_h + recovered_n,
        -metabolism.phosphorus_growth_actual_g_p_per_h + recovered_p,
    };
}

/// Publishes all primary and secondary axes sharing one root-layer mobile
/// pool as one rollback-safe transaction. Every flux must have been staged
/// from the same pre-commit snapshot in `AxisWorkspace`.
pub fn commitStagedLayerAxes(
    roots: *RootState,
    plant: usize,
    domain: usize,
    layer: usize,
    workspace: *AxisWorkspace,
    active_axis_count: usize,
    parameters: StagedLayerCommitParameters,
) !void {
    if (active_axis_count > workspace.axis_capacity or active_axis_count > roots.root_axis_count) return error.RootMetabolismWorkspaceCapacityExceeded;
    inline for (@typeInfo(StagedLayerCommitParameters).@"struct".fields) |field| {
        const value = @field(parameters, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidStagedRootCommitParameter;
    }
    inline for (.{ parameters.root_extension_water_response, parameters.nonwoody_carbon_fraction, parameters.nonwoody_nitrogen_fraction, parameters.nonwoody_phosphorus_fraction }) |value| if (value > 1) return error.InvalidStagedRootCommitParameter;
    const root_layer = try roots.layerIndex(plant, domain, layer);
    var mobile_change = [_]f64{0} ** 3;
    var protein_change_current: f64 = 0;
    var protein_change_next: f64 = 0;
    var crossing_transfer_retained_fraction: f64 = 1;
    var any_primary_crossing = false;
    var secondary_respiration = Respiration{ .actual_g_c = 0, .oxygen_unlimited_g_c = 0, .carbon_unlimited_g_c = 0 };

    for (0..active_axis_count) |axis| {
        const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
        if (workspace.primary_active[axis]) {
            const metabolism = workspace.primary_metabolism[axis];
            const senescence = workspace.primary_senescence[axis];
            const fraction = senescence.senesced_fraction;
            const recovered = [3]f64{
                fraction * senescence.recyclable_carbon_g_c * parameters.nonwoody_carbon_fraction,
                fraction * senescence.recyclable_nitrogen_g_n * parameters.nonwoody_nitrogen_fraction,
                fraction * senescence.recyclable_phosphorus_g_p,
            };
            const changes = mobileChanges(metabolism, senescence, recovered[0], recovered[1], recovered[2]);
            for (&mobile_change, changes) |*total, change| total.* += change;
            var total_primary_carbon_g_c: f64 = 0;
            var total_primary_nitrogen_g_n: f64 = 0;
            var total_primary_phosphorus_g_p: f64 = 0;
            for (0..roots.soil_layer_count) |axis_carbon_layer| {
                const total_axis_layer = try roots.layerAxisIndex(plant, domain, axis_carbon_layer, axis);
                total_primary_carbon_g_c += roots.axis_primary_carbon_g[total_axis_layer];
                total_primary_nitrogen_g_n += roots.axis_primary_nitrogen_g[total_axis_layer];
                total_primary_phosphorus_g_p += roots.axis_primary_phosphorus_g[total_axis_layer];
            }
            var net_growth_g_c = metabolism.root_growth_actual_g_c_per_h - fraction * total_primary_carbon_g_c;
            var net_growth_g_n = metabolism.nitrogen_growth_actual_g_n_per_h - fraction * total_primary_nitrogen_g_n;
            var net_growth_g_p = metabolism.phosphorus_growth_actual_g_p_per_h - fraction * total_primary_phosphorus_g_p;
            if (net_growth_g_c < 0 or net_growth_g_n < 0 or net_growth_g_p < 0) {
                const current_secondary = try stagedSecondaryRootLayer(
                    roots.*,
                    axis_layer,
                    workspace.secondary_active[axis],
                    workspace.secondary_metabolism[axis],
                    workspace.secondary_senescence[axis],
                    parameters.secondary_specific_length_m_per_g_c,
                    parameters.root_extension_water_response,
                );
                const upper_secondary: SecondaryRootDeficitLayer = if (layer > 0) blk: {
                    const upper_axis_layer = try roots.layerAxisIndex(plant, domain, layer - 1, axis);
                    break :blk .{
                        .carbon_g_c = roots.axis_secondary_carbon_g[upper_axis_layer],
                        .nitrogen_g_n = roots.axis_secondary_nitrogen_g[upper_axis_layer],
                        .phosphorus_g_p = roots.axis_secondary_phosphorus_g[upper_axis_layer],
                        .length_m = roots.axis_secondary_length_m[upper_axis_layer],
                    };
                } else .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0, .length_m = 0 };
                const absorption = try absorbPrimaryDeficitFromSecondaryRoots(
                    @max(0, -net_growth_g_c),
                    @max(0, -net_growth_g_n),
                    @max(0, -net_growth_g_p),
                    current_secondary,
                    upper_secondary,
                );
                workspace.primary_deficit_absorption[axis] = absorption;
                workspace.primary_deficit_active[axis] = true;
                if (net_growth_g_c < 0) net_growth_g_c = -absorption.residual_carbon_deficit_g_c;
                if (net_growth_g_n < 0) net_growth_g_n = -absorption.residual_nitrogen_deficit_g_n;
                if (net_growth_g_p < 0) net_growth_g_p = -absorption.residual_phosphorus_deficit_g_p;
            }
            const axis_index = try roots.axisIndex(plant, domain, axis);
            const requested_extension_m = try primaryRootLengthChange(metabolism.root_growth_actual_g_c_per_h, net_growth_g_c, total_primary_carbon_g_c, roots.axis_depth_m[axis_index], parameters.seeding_depth_m, parameters.primary_specific_length_m_per_g_c, parameters.plant_population_count, parameters.root_extension_water_response);
            const placement = try sourceOrderPrimaryRootExtensionPlacement(requested_extension_m, roots.axis_depth_m[axis_index], parameters.current_layer_bottom_depth_m, parameters.next_layer_thickness_m, parameters.extension_presence_threshold_m);
            const target_layer = if (placement.crosses_into_next_layer) layer + 1 else layer;
            if (target_layer >= roots.soil_layer_count) return error.StagedRootCommitLayerOutOfBounds;
            const target_axis_layer = try roots.layerAxisIndex(plant, domain, target_layer, axis);
            const next_c = roots.axis_primary_carbon_g[target_axis_layer] + net_growth_g_c;
            const next_n = roots.axis_primary_nitrogen_g[target_axis_layer] + net_growth_g_n;
            const next_p = roots.axis_primary_phosphorus_g[target_axis_layer] + net_growth_g_p;
            const next_length = roots.axis_primary_length_m[target_axis_layer] + placement.extension_m;
            const next_depth = roots.axis_depth_m[axis_index] + placement.extension_m;
            inline for (.{ next_c, next_n, next_p, next_length, next_depth }) |value| if (!std.math.isFinite(value) or value < -1e-12) return error.StagedRootCommitWouldOverdrawPool;
            const protein = @min(parameters.protein_carbon_per_nitrogen_g_c_per_g_n * next_n, parameters.protein_carbon_per_phosphorus_g_c_per_g_p * next_p);
            if (placement.crosses_into_next_layer) {
                any_primary_crossing = true;
                const transfer_fraction = workspace.primary_sink_fractions[axis];
                if (!std.math.isFinite(transfer_fraction) or transfer_fraction < 0 or transfer_fraction > 1) return error.InvalidStagedRootCommitParameter;
                crossing_transfer_retained_fraction *= 1 - transfer_fraction;
                protein_change_current += @min(
                    parameters.protein_carbon_per_nitrogen_g_c_per_g_n * roots.axis_primary_nitrogen_g[axis_layer],
                    parameters.protein_carbon_per_phosphorus_g_c_per_g_p * roots.axis_primary_phosphorus_g[axis_layer],
                );
                protein_change_next += protein;
            } else {
                protein_change_current += protein;
            }
        }
        if (workspace.secondary_active[axis]) {
            const metabolism = workspace.secondary_metabolism[axis];
            const senescence = workspace.secondary_senescence[axis];
            const fraction = senescence.senesced_fraction;
            const recovered = [3]f64{
                fraction * senescence.recyclable_carbon_g_c * parameters.nonwoody_carbon_fraction,
                fraction * senescence.recyclable_nitrogen_g_n * parameters.nonwoody_nitrogen_fraction,
                fraction * senescence.recyclable_phosphorus_g_p * parameters.nonwoody_phosphorus_fraction,
            };
            const changes = mobileChanges(metabolism, senescence, recovered[0], recovered[1], recovered[2]);
            for (&mobile_change, changes) |*total, change| total.* += change;
            const next_c = roots.axis_secondary_carbon_g[axis_layer] + metabolism.root_growth_actual_g_c_per_h - fraction * roots.axis_secondary_carbon_g[axis_layer];
            const next_n = roots.axis_secondary_nitrogen_g[axis_layer] + metabolism.nitrogen_growth_actual_g_n_per_h - fraction * roots.axis_secondary_nitrogen_g[axis_layer];
            const next_p = roots.axis_secondary_phosphorus_g[axis_layer] + metabolism.phosphorus_growth_actual_g_p_per_h - fraction * roots.axis_secondary_phosphorus_g[axis_layer];
            const next_length = roots.axis_secondary_length_m[axis_layer] + metabolism.root_growth_actual_g_c_per_h * parameters.secondary_specific_length_m_per_g_c * parameters.root_extension_water_response - fraction * roots.axis_secondary_length_m[axis_layer];
            inline for (.{ next_c, next_n, next_p, next_length }) |value| if (!std.math.isFinite(value) or value < -1e-12) return error.StagedRootCommitWouldOverdrawPool;
            protein_change_current += @min(parameters.protein_carbon_per_nitrogen_g_c_per_g_n * next_n, parameters.protein_carbon_per_phosphorus_g_c_per_g_p * next_p);
            const respiration = try assemble(.{
                .maintenance_demand_g_c = metabolism.maintenance_respiration_g_c_per_h,
                .substrate_respiration_actual_g_c = metabolism.substrate_respiration_actual_g_c_per_h,
                .substrate_respiration_oxygen_unlimited_g_c = metabolism.substrate_respiration_oxygen_unlimited_g_c_per_h,
                .growth_respiration_actual_g_c = metabolism.growth_respiration_actual_g_c_per_h,
                .growth_respiration_oxygen_unlimited_g_c = metabolism.growth_respiration_oxygen_unlimited_g_c_per_h,
                .senescence_respiration_actual_g_c = senescence.respiration_actual_g_c_per_h,
                .senescence_respiration_oxygen_unlimited_g_c = senescence.respiration_oxygen_unlimited_g_c_per_h,
                .nitrogen_assimilation_respiration_actual_g_c = metabolism.nitrogen_assimilation_respiration_actual_g_c_per_h,
                .nitrogen_assimilation_respiration_oxygen_unlimited_g_c = metabolism.nitrogen_assimilation_respiration_oxygen_unlimited_g_c_per_h,
            });
            secondary_respiration.actual_g_c += respiration.actual_g_c;
            secondary_respiration.oxygen_unlimited_g_c += respiration.oxygen_unlimited_g_c;
            secondary_respiration.carbon_unlimited_g_c += respiration.carbon_unlimited_g_c;
        }
    }
    const next_mobile = [3]f64{
        roots.mobile_carbon_g[root_layer] + mobile_change[0],
        roots.mobile_nitrogen_g[root_layer] + mobile_change[1],
        roots.mobile_phosphorus_g[root_layer] + mobile_change[2],
    };
    const crosses_layer = any_primary_crossing;
    const next_root_layer = if (crosses_layer) try roots.layerIndex(plant, domain, layer + 1) else root_layer;
    const transfer_fraction = 1 - crossing_transfer_retained_fraction;
    const current_mobile = [3]f64{
        next_mobile[0] * crossing_transfer_retained_fraction,
        next_mobile[1] * crossing_transfer_retained_fraction,
        next_mobile[2] * crossing_transfer_retained_fraction,
    };
    const receiving_mobile = [3]f64{
        roots.mobile_carbon_g[next_root_layer] + next_mobile[0] * transfer_fraction,
        roots.mobile_nitrogen_g[next_root_layer] + next_mobile[1] * transfer_fraction,
        roots.mobile_phosphorus_g[next_root_layer] + next_mobile[2] * transfer_fraction,
    };
    const next_protein = roots.protein_carbon_g[root_layer] + protein_change_current;
    const receiving_protein = roots.protein_carbon_g[next_root_layer] + protein_change_next;
    const next_actual = roots.actual_respiration_g_c_per_h[root_layer] + secondary_respiration.actual_g_c;
    const next_oxygen_unlimited = roots.respiration_unlimited_by_oxygen_g_c_per_h[root_layer] + secondary_respiration.oxygen_unlimited_g_c;
    const next_carbon_unlimited = roots.respiration_unlimited_by_carbon_g_c_per_h[root_layer] + secondary_respiration.carbon_unlimited_g_c;
    inline for (.{
        current_mobile[0],
        current_mobile[1],
        current_mobile[2],
        receiving_mobile[0],
        receiving_mobile[1],
        receiving_mobile[2],
        next_protein,
        receiving_protein,
        next_actual,
        next_oxygen_unlimited,
        next_carbon_unlimited,
    }) |value| if (!std.math.isFinite(value) or value < -1e-12) return error.StagedRootCommitWouldOverdrawPool;

    roots.mobile_carbon_g[root_layer] = @max(0, current_mobile[0]);
    roots.mobile_nitrogen_g[root_layer] = @max(0, current_mobile[1]);
    roots.mobile_phosphorus_g[root_layer] = @max(0, current_mobile[2]);
    roots.protein_carbon_g[root_layer] = @max(0, next_protein);
    if (crosses_layer) {
        roots.mobile_carbon_g[next_root_layer] = @max(0, receiving_mobile[0]);
        roots.mobile_nitrogen_g[next_root_layer] = @max(0, receiving_mobile[1]);
        roots.mobile_phosphorus_g[next_root_layer] = @max(0, receiving_mobile[2]);
        roots.protein_carbon_g[next_root_layer] = @max(0, receiving_protein);
        roots.total_water_potential_mpa[next_root_layer] = roots.total_water_potential_mpa[root_layer];
        roots.osmotic_water_potential_mpa[next_root_layer] = roots.osmotic_water_potential_mpa[root_layer];
        roots.turgor_water_potential_mpa[next_root_layer] = roots.turgor_water_potential_mpa[root_layer];
        roots.primary_radius_m[next_root_layer] = roots.primary_radius_m[root_layer];
    }
    roots.actual_respiration_g_c_per_h[root_layer] = next_actual;
    roots.respiration_unlimited_by_oxygen_g_c_per_h[root_layer] = next_oxygen_unlimited;
    roots.respiration_unlimited_by_carbon_g_c_per_h[root_layer] = next_carbon_unlimited;
    for (0..active_axis_count) |axis| {
        const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
        if (workspace.primary_active[axis]) {
            const metabolism = workspace.primary_metabolism[axis];
            const fraction = workspace.primary_senescence[axis].senesced_fraction;
            var total_primary_carbon_g_c: f64 = 0;
            var total_primary_nitrogen_g_n: f64 = 0;
            var total_primary_phosphorus_g_p: f64 = 0;
            for (0..roots.soil_layer_count) |axis_carbon_layer| {
                const total_axis_layer = try roots.layerAxisIndex(plant, domain, axis_carbon_layer, axis);
                total_primary_carbon_g_c += roots.axis_primary_carbon_g[total_axis_layer];
                total_primary_nitrogen_g_n += roots.axis_primary_nitrogen_g[total_axis_layer];
                total_primary_phosphorus_g_p += roots.axis_primary_phosphorus_g[total_axis_layer];
            }
            var net_growth_g_c = metabolism.root_growth_actual_g_c_per_h - fraction * total_primary_carbon_g_c;
            var net_growth_g_n = metabolism.nitrogen_growth_actual_g_n_per_h - fraction * total_primary_nitrogen_g_n;
            var net_growth_g_p = metabolism.phosphorus_growth_actual_g_p_per_h - fraction * total_primary_phosphorus_g_p;
            const absorption = workspace.primary_deficit_absorption[axis];
            if (net_growth_g_c < 0) net_growth_g_c = -absorption.residual_carbon_deficit_g_c;
            if (net_growth_g_n < 0) net_growth_g_n = -absorption.residual_nitrogen_deficit_g_n;
            if (net_growth_g_p < 0) net_growth_g_p = -absorption.residual_phosphorus_deficit_g_p;
            const axis_index = try roots.axisIndex(plant, domain, axis);
            const requested_extension_m = try primaryRootLengthChange(metabolism.root_growth_actual_g_c_per_h, net_growth_g_c, total_primary_carbon_g_c, roots.axis_depth_m[axis_index], parameters.seeding_depth_m, parameters.primary_specific_length_m_per_g_c, parameters.plant_population_count, parameters.root_extension_water_response);
            const placement = try sourceOrderPrimaryRootExtensionPlacement(requested_extension_m, roots.axis_depth_m[axis_index], parameters.current_layer_bottom_depth_m, parameters.next_layer_thickness_m, parameters.extension_presence_threshold_m);
            const target_layer = if (placement.crosses_into_next_layer) layer + 1 else layer;
            const target_axis_layer = try roots.layerAxisIndex(plant, domain, target_layer, axis);
            roots.axis_primary_carbon_g[target_axis_layer] = @max(0, roots.axis_primary_carbon_g[target_axis_layer] + net_growth_g_c);
            roots.axis_primary_nitrogen_g[target_axis_layer] = @max(0, roots.axis_primary_nitrogen_g[target_axis_layer] + net_growth_g_n);
            roots.axis_primary_phosphorus_g[target_axis_layer] = @max(0, roots.axis_primary_phosphorus_g[target_axis_layer] + net_growth_g_p);
            roots.axis_primary_length_m[target_axis_layer] = @max(0, roots.axis_primary_length_m[target_axis_layer] + placement.extension_m);
            roots.axis_depth_m[axis_index] = @max(0, roots.axis_depth_m[axis_index] + placement.extension_m);
        }
        if (workspace.secondary_active[axis]) {
            const metabolism = workspace.secondary_metabolism[axis];
            const fraction = workspace.secondary_senescence[axis].senesced_fraction;
            const extension = metabolism.root_growth_actual_g_c_per_h * parameters.secondary_specific_length_m_per_g_c * parameters.root_extension_water_response;
            roots.axis_secondary_carbon_g[axis_layer] = @max(0, roots.axis_secondary_carbon_g[axis_layer] + metabolism.root_growth_actual_g_c_per_h - fraction * roots.axis_secondary_carbon_g[axis_layer]);
            roots.axis_secondary_nitrogen_g[axis_layer] = @max(0, roots.axis_secondary_nitrogen_g[axis_layer] + metabolism.nitrogen_growth_actual_g_n_per_h - fraction * roots.axis_secondary_nitrogen_g[axis_layer]);
            roots.axis_secondary_phosphorus_g[axis_layer] = @max(0, roots.axis_secondary_phosphorus_g[axis_layer] + metabolism.phosphorus_growth_actual_g_p_per_h - fraction * roots.axis_secondary_phosphorus_g[axis_layer]);
            roots.axis_secondary_length_m[axis_layer] = @max(0, roots.axis_secondary_length_m[axis_layer] + extension - fraction * roots.axis_secondary_length_m[axis_layer]);
        }
    }
    for (0..active_axis_count) |axis| {
        if (!workspace.primary_deficit_active[axis]) continue;
        const absorption = workspace.primary_deficit_absorption[axis];
        const current_axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
        roots.axis_secondary_carbon_g[current_axis_layer] = absorption.current.carbon_g_c;
        roots.axis_secondary_nitrogen_g[current_axis_layer] = absorption.current.nitrogen_g_n;
        roots.axis_secondary_phosphorus_g[current_axis_layer] = absorption.current.phosphorus_g_p;
        roots.axis_secondary_length_m[current_axis_layer] = absorption.current.length_m;
        if (layer > 0) {
            const upper_axis_layer = try roots.layerAxisIndex(plant, domain, layer - 1, axis);
            roots.axis_secondary_carbon_g[upper_axis_layer] = absorption.upper.carbon_g_c;
            roots.axis_secondary_nitrogen_g[upper_axis_layer] = absorption.upper.nitrogen_g_n;
            roots.axis_secondary_phosphorus_g[upper_axis_layer] = absorption.upper.phosphorus_g_p;
            roots.axis_secondary_length_m[upper_axis_layer] = absorption.upper.length_m;
        }
    }
}

/// Atomic GROSUB primary-axis C/N/P, length/depth, mobile pools, protein, and
/// RCO2 transaction. The caller stages all axis calculations before invoking
/// commits when several axes share one mobile pool.
pub fn commitPrimaryRoot(
    roots: *RootState,
    root_layer: usize,
    root_axis_layer: usize,
    root_axis: usize,
    inputs: PrimaryRootCommitInputs,
) !void {
    if (root_layer >= roots.mobile_carbon_g.len or root_axis_layer >= roots.axis_primary_carbon_g.len or root_axis >= roots.axis_depth_m.len) return error.PlantRootIndexOutOfBounds;
    inline for (.{ inputs.primary_specific_length_m_per_g_c, inputs.root_extension_water_response, inputs.nonwoody_carbon_fraction, inputs.nonwoody_nitrogen_fraction, inputs.protein_carbon_per_nitrogen_g_c_per_g_n, inputs.protein_carbon_per_phosphorus_g_c_per_g_p }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootCommitInput;
    inline for (.{ inputs.root_extension_water_response, inputs.nonwoody_carbon_fraction, inputs.nonwoody_nitrogen_fraction }) |value| if (value > 1) return error.InvalidPrimaryRootCommitInput;

    const old_c = roots.axis_primary_carbon_g[root_axis_layer];
    const old_n = roots.axis_primary_nitrogen_g[root_axis_layer];
    const old_p = roots.axis_primary_phosphorus_g[root_axis_layer];
    const old_length = roots.axis_primary_length_m[root_axis_layer];
    const old_depth = roots.axis_depth_m[root_axis];
    const fraction = inputs.senescence.senesced_fraction;
    const recovered_c = fraction * inputs.senescence.recyclable_carbon_g_c * inputs.nonwoody_carbon_fraction;
    const recovered_n = fraction * inputs.senescence.recyclable_nitrogen_g_n * inputs.nonwoody_nitrogen_fraction;
    // GROSUB primary P recovery has no FWODRP(1) multiplier.
    const recovered_p = fraction * inputs.senescence.recyclable_phosphorus_g_p;
    const next_mobile_c = roots.mobile_carbon_g[root_layer] -
        @min(inputs.metabolism.maintenance_respiration_g_c_per_h, inputs.metabolism.substrate_respiration_actual_g_c_per_h) -
        inputs.metabolism.growth_and_respiration_carbon_actual_g_c_per_h -
        inputs.metabolism.nitrogen_assimilation_respiration_actual_g_c_per_h -
        inputs.senescence.respiration_actual_g_c_per_h + recovered_c;
    const next_mobile_n = roots.mobile_nitrogen_g[root_layer] - inputs.metabolism.nitrogen_growth_actual_g_n_per_h + recovered_n;
    const next_mobile_p = roots.mobile_phosphorus_g[root_layer] - inputs.metabolism.phosphorus_growth_actual_g_p_per_h + recovered_p;
    const next_c = old_c + inputs.metabolism.root_growth_actual_g_c_per_h - fraction * old_c;
    const next_n = old_n + inputs.metabolism.nitrogen_growth_actual_g_n_per_h - fraction * old_n;
    const next_p = old_p + inputs.metabolism.phosphorus_growth_actual_g_p_per_h - fraction * old_p;
    const extension_m = inputs.metabolism.root_growth_actual_g_c_per_h * inputs.primary_specific_length_m_per_g_c * inputs.root_extension_water_response;
    const next_length = old_length + extension_m - fraction * old_length;
    const next_depth = old_depth + extension_m;
    const next_protein = roots.protein_carbon_g[root_layer] + @min(
        inputs.protein_carbon_per_nitrogen_g_c_per_g_n * next_n,
        inputs.protein_carbon_per_phosphorus_g_c_per_g_p * next_p,
    );
    inline for (.{ next_mobile_c, next_mobile_n, next_mobile_p, next_c, next_n, next_p, next_length, next_depth, next_protein }) |value| {
        if (!std.math.isFinite(value)) return error.NonFinitePrimaryRootCommit;
        if (value < -1.0e-12) return error.PrimaryRootCommitWouldOverdrawPool;
    }
    roots.mobile_carbon_g[root_layer] = @max(0, next_mobile_c);
    roots.mobile_nitrogen_g[root_layer] = @max(0, next_mobile_n);
    roots.mobile_phosphorus_g[root_layer] = @max(0, next_mobile_p);
    roots.axis_primary_carbon_g[root_axis_layer] = @max(0, next_c);
    roots.axis_primary_nitrogen_g[root_axis_layer] = @max(0, next_n);
    roots.axis_primary_phosphorus_g[root_axis_layer] = @max(0, next_p);
    roots.axis_primary_length_m[root_axis_layer] = @max(0, next_length);
    roots.axis_depth_m[root_axis] = @max(0, next_depth);
    roots.protein_carbon_g[root_layer] = @max(0, next_protein);
}

/// GROSUB FRCO2 allocation of one primary-axis respiration total through all
/// traversed soil layers. Caller-owned fractions avoid per-hour allocation.
pub fn allocatePrimaryRootRespiration(
    roots: *RootState,
    root_layer_indices: []const usize,
    primary_length_m_by_layer: []const f64,
    allocation_fractions: []f64,
    primary_depth_m: f64,
    seeding_depth_m: f64,
    traverses_multiple_layers: bool,
    respiration: Respiration,
) !void {
    if (root_layer_indices.len == 0 or root_layer_indices.len != primary_length_m_by_layer.len or root_layer_indices.len != allocation_fractions.len) return error.PrimaryRootRespirationDimensionMismatch;
    inline for (.{ primary_depth_m, seeding_depth_m, respiration.actual_g_c, respiration.oxygen_unlimited_g_c, respiration.carbon_unlimited_g_c }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPrimaryRootRespirationAllocation;
    const denominator = primary_depth_m - seeding_depth_m;
    if (traverses_multiple_layers and denominator <= 0) return error.InvalidPrimaryRootRespirationAllocation;
    var allocated: f64 = 0;
    for (primary_length_m_by_layer, 0..) |length, layer| {
        if (!std.math.isFinite(length) or length < 0 or root_layer_indices[layer] >= roots.actual_respiration_g_c_per_h.len) return error.InvalidPrimaryRootRespirationAllocation;
        allocation_fractions[layer] = if (!traverses_multiple_layers)
            @floatFromInt(@intFromBool(layer + 1 == primary_length_m_by_layer.len))
        else if (layer + 1 < primary_length_m_by_layer.len)
            @min(1, length / denominator)
        else
            @max(0, 1 - allocated);
        allocated += allocation_fractions[layer];
    }
    if (@abs(allocated - 1) > 1e-10) return error.InvalidPrimaryRootRespirationAllocation;
    for (root_layer_indices, allocation_fractions) |root, fraction| inline for (.{
        roots.actual_respiration_g_c_per_h[root] + respiration.actual_g_c * fraction,
        roots.respiration_unlimited_by_oxygen_g_c_per_h[root] + respiration.oxygen_unlimited_g_c * fraction,
        roots.respiration_unlimited_by_carbon_g_c_per_h[root] + respiration.carbon_unlimited_g_c * fraction,
    }) |value| if (!std.math.isFinite(value)) return error.NonFinitePrimaryRootRespirationAllocation;
    for (root_layer_indices, allocation_fractions) |root, fraction| {
        roots.actual_respiration_g_c_per_h[root] += respiration.actual_g_c * fraction;
        roots.respiration_unlimited_by_oxygen_g_c_per_h[root] += respiration.oxygen_unlimited_g_c * fraction;
        roots.respiration_unlimited_by_carbon_g_c_per_h[root] += respiration.carbon_unlimited_g_c * fraction;
    }
}

/// Atomic GROSUB CPOOLR/ZPOOLR/PPOOLR, RTLG2, WTRT2/N/P, WSRTL and
/// RCO2M/RCO2N/RCO2A transaction for one runtime root axis and soil layer.
pub fn commitSecondaryRoot(
    roots: *RootState,
    root_layer: usize,
    root_axis_layer: usize,
    inputs: SecondaryRootCommitInputs,
) !void {
    if (root_layer >= roots.mobile_carbon_g.len or root_axis_layer >= roots.axis_secondary_carbon_g.len) return error.PlantRootIndexOutOfBounds;
    inline for (.{ inputs.root_specific_length_m_per_g_c, inputs.root_extension_water_response, inputs.nonwoody_carbon_fraction, inputs.nonwoody_nitrogen_fraction, inputs.nonwoody_phosphorus_fraction, inputs.protein_carbon_per_nitrogen_g_c_per_g_n, inputs.protein_carbon_per_phosphorus_g_c_per_g_p }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSecondaryRootCommitInput;
    inline for (.{ inputs.nonwoody_carbon_fraction, inputs.nonwoody_nitrogen_fraction, inputs.nonwoody_phosphorus_fraction, inputs.root_extension_water_response }) |value| if (value > 1) return error.InvalidSecondaryRootCommitInput;

    const old_c = roots.axis_secondary_carbon_g[root_axis_layer];
    const old_n = roots.axis_secondary_nitrogen_g[root_axis_layer];
    const old_p = roots.axis_secondary_phosphorus_g[root_axis_layer];
    const old_length = roots.axis_secondary_length_m[root_axis_layer];
    const fraction = inputs.senescence.senesced_fraction;
    const recovered_c = fraction * inputs.senescence.recyclable_carbon_g_c * inputs.nonwoody_carbon_fraction;
    const recovered_n = fraction * inputs.senescence.recyclable_nitrogen_g_n * inputs.nonwoody_nitrogen_fraction;
    const recovered_p = fraction * inputs.senescence.recyclable_phosphorus_g_p * inputs.nonwoody_phosphorus_fraction;
    const next_mobile_c = roots.mobile_carbon_g[root_layer] -
        @min(inputs.metabolism.maintenance_respiration_g_c_per_h, inputs.metabolism.substrate_respiration_actual_g_c_per_h) -
        inputs.metabolism.growth_and_respiration_carbon_actual_g_c_per_h -
        inputs.metabolism.nitrogen_assimilation_respiration_actual_g_c_per_h -
        inputs.senescence.respiration_actual_g_c_per_h + recovered_c;
    const next_mobile_n = roots.mobile_nitrogen_g[root_layer] - inputs.metabolism.nitrogen_growth_actual_g_n_per_h + recovered_n;
    const next_mobile_p = roots.mobile_phosphorus_g[root_layer] - inputs.metabolism.phosphorus_growth_actual_g_p_per_h + recovered_p;
    const next_c = old_c + inputs.metabolism.root_growth_actual_g_c_per_h - fraction * old_c;
    const next_n = old_n + inputs.metabolism.nitrogen_growth_actual_g_n_per_h - fraction * old_n;
    const next_p = old_p + inputs.metabolism.phosphorus_growth_actual_g_p_per_h - fraction * old_p;
    const next_length = old_length +
        inputs.metabolism.root_growth_actual_g_c_per_h * inputs.root_specific_length_m_per_g_c * inputs.root_extension_water_response -
        fraction * old_length;
    const next_protein = roots.protein_carbon_g[root_layer] + @min(
        inputs.protein_carbon_per_nitrogen_g_c_per_g_n * next_n,
        inputs.protein_carbon_per_phosphorus_g_c_per_g_p * next_p,
    );
    const respiration = try assemble(.{
        .maintenance_demand_g_c = inputs.metabolism.maintenance_respiration_g_c_per_h,
        .substrate_respiration_actual_g_c = inputs.metabolism.substrate_respiration_actual_g_c_per_h,
        .substrate_respiration_oxygen_unlimited_g_c = inputs.metabolism.substrate_respiration_oxygen_unlimited_g_c_per_h,
        .growth_respiration_actual_g_c = inputs.metabolism.growth_respiration_actual_g_c_per_h,
        .growth_respiration_oxygen_unlimited_g_c = inputs.metabolism.growth_respiration_oxygen_unlimited_g_c_per_h,
        .senescence_respiration_actual_g_c = inputs.senescence.respiration_actual_g_c_per_h,
        .senescence_respiration_oxygen_unlimited_g_c = inputs.senescence.respiration_oxygen_unlimited_g_c_per_h,
        .nitrogen_assimilation_respiration_actual_g_c = inputs.metabolism.nitrogen_assimilation_respiration_actual_g_c_per_h,
        .nitrogen_assimilation_respiration_oxygen_unlimited_g_c = inputs.metabolism.nitrogen_assimilation_respiration_oxygen_unlimited_g_c_per_h,
    });
    const next_actual = roots.actual_respiration_g_c_per_h[root_layer] + respiration.actual_g_c;
    const next_oxygen_unlimited = roots.respiration_unlimited_by_oxygen_g_c_per_h[root_layer] + respiration.oxygen_unlimited_g_c;
    const next_carbon_unlimited = roots.respiration_unlimited_by_carbon_g_c_per_h[root_layer] + respiration.carbon_unlimited_g_c;
    inline for (.{ next_mobile_c, next_mobile_n, next_mobile_p, next_c, next_n, next_p, next_length, next_protein, next_actual, next_oxygen_unlimited, next_carbon_unlimited }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteSecondaryRootCommit;
        if (value < -1.0e-12) return error.SecondaryRootCommitWouldOverdrawPool;
    }
    roots.mobile_carbon_g[root_layer] = @max(0, next_mobile_c);
    roots.mobile_nitrogen_g[root_layer] = @max(0, next_mobile_n);
    roots.mobile_phosphorus_g[root_layer] = @max(0, next_mobile_p);
    roots.axis_secondary_carbon_g[root_axis_layer] = @max(0, next_c);
    roots.axis_secondary_nitrogen_g[root_axis_layer] = @max(0, next_n);
    roots.axis_secondary_phosphorus_g[root_axis_layer] = @max(0, next_p);
    roots.axis_secondary_length_m[root_axis_layer] = @max(0, next_length);
    roots.protein_carbon_g[root_layer] = @max(0, next_protein);
    roots.actual_respiration_g_c_per_h[root_layer] = next_actual;
    roots.respiration_unlimited_by_oxygen_g_c_per_h[root_layer] = next_oxygen_unlimited;
    roots.respiration_unlimited_by_carbon_g_c_per_h[root_layer] = next_carbon_unlimited;
}

/// GROSUB RCO2T/RCO2TM assembly for primary and secondary roots.
pub fn assemble(components: Components) !Respiration {
    inline for (@typeInfo(Components).@"struct".fields) |field| {
        const value = @field(components, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidRootRespirationComponent;
    }
    const actual =
        @min(components.maintenance_demand_g_c, components.substrate_respiration_actual_g_c) +
        components.growth_respiration_actual_g_c +
        components.senescence_respiration_actual_g_c +
        components.nitrogen_assimilation_respiration_actual_g_c;
    const oxygen_unlimited =
        @min(components.maintenance_demand_g_c, components.substrate_respiration_oxygen_unlimited_g_c) +
        components.growth_respiration_oxygen_unlimited_g_c +
        components.senescence_respiration_oxygen_unlimited_g_c +
        components.nitrogen_assimilation_respiration_oxygen_unlimited_g_c;
    const result: Respiration = .{
        .actual_g_c = actual,
        .oxygen_unlimited_g_c = oxygen_unlimited,
        // GROSUB publishes RCO2T into RCO2N. This is the demand subsequently
        // compared with the mobile-C pool by UPTAKE's FCUP.
        .carbon_unlimited_g_c = actual,
    };
    inline for (@typeInfo(Respiration).@"struct".fields) |field| if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteRootRespiration;
    return result;
}

/// GROSUB CUPRL/CUPRO/CUPRC. Results are the eight NH4, NO3, H2PO4,
/// and HPO4 band/non-band uptake pools for one root domain and layer.
pub fn nutrientUptakeRespiration(results: []const NutrientResult, respiration_g_c_per_g_element: f64) !Respiration {
    if (results.len != @import("plant_root_nutrient_uptake.zig").nutrient_pool_count) return error.RootNutrientPoolCountMismatch;
    if (!std.math.isFinite(respiration_g_c_per_g_element) or respiration_g_c_per_g_element < 0) return error.InvalidRootNutrientRespirationCoefficient;
    var actual: f64 = 0;
    var oxygen_unlimited: f64 = 0;
    var carbon_unlimited: f64 = 0;
    for (results) |result| {
        inline for (@typeInfo(NutrientResult).@"struct".fields) |field| if (!std.math.isFinite(@field(result, field.name)) or @field(result, field.name) < 0) return error.InvalidRootNutrientRespirationInput;
        actual += result.uptake_g_element;
        oxygen_unlimited += result.oxygen_unlimited_uptake_g_element;
        carbon_unlimited += result.carbon_unlimited_uptake_g_element;
    }
    return .{
        .actual_g_c = respiration_g_c_per_g_element * actual,
        .oxygen_unlimited_g_c = respiration_g_c_per_g_element * oxygen_unlimited,
        .carbon_unlimited_g_c = respiration_g_c_per_g_element * carbon_unlimited,
    };
}

/// Atomic root respiration publication. Actual respiration consumes mobile
/// carbon and is recorded as a positive CO2 production rate.
pub fn commit(roots: *RootState, root_layer: usize, respiration: Respiration) !void {
    if (root_layer >= roots.mobile_carbon_g.len) return error.PlantRootIndexOutOfBounds;
    inline for (@typeInfo(Respiration).@"struct".fields) |field| if (!std.math.isFinite(@field(respiration, field.name)) or @field(respiration, field.name) < 0) return error.InvalidRootRespirationCommit;
    const next_mobile = roots.mobile_carbon_g[root_layer] - respiration.actual_g_c;
    const next_actual = roots.actual_respiration_g_c_per_h[root_layer] + respiration.actual_g_c;
    const next_oxygen_unlimited = roots.respiration_unlimited_by_oxygen_g_c_per_h[root_layer] + respiration.oxygen_unlimited_g_c;
    const next_carbon_unlimited = roots.respiration_unlimited_by_carbon_g_c_per_h[root_layer] + respiration.carbon_unlimited_g_c;
    inline for (.{ next_mobile, next_actual, next_oxygen_unlimited, next_carbon_unlimited }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootRespirationCommit;
    if (next_mobile < -1.0e-12) return error.InsufficientRootMobileCarbonForRespiration;
    roots.mobile_carbon_g[root_layer] = @max(0, next_mobile);
    roots.actual_respiration_g_c_per_h[root_layer] = next_actual;
    roots.respiration_unlimited_by_oxygen_g_c_per_h[root_layer] = next_oxygen_unlimited;
    roots.respiration_unlimited_by_carbon_g_c_per_h[root_layer] = next_carbon_unlimited;
}

test "GROSUB root respiration assembly preserves RCO2T and RCO2TM equations" {
    const components: Components = .{
        .maintenance_demand_g_c = 3,
        .substrate_respiration_actual_g_c = 2,
        .substrate_respiration_oxygen_unlimited_g_c = 4,
        .growth_respiration_actual_g_c = 0.5,
        .growth_respiration_oxygen_unlimited_g_c = 0.8,
        .senescence_respiration_actual_g_c = 0.2,
        .senescence_respiration_oxygen_unlimited_g_c = 0.3,
        .nitrogen_assimilation_respiration_actual_g_c = 0.1,
        .nitrogen_assimilation_respiration_oxygen_unlimited_g_c = 0.15,
    };
    const result = try assemble(components);
    try std.testing.expectApproxEqAbs(@as(f64, 2.8), result.actual_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 4.25), result.oxygen_unlimited_g_c, 1.0e-12);
    try std.testing.expectEqual(result.actual_g_c, result.carbon_unlimited_g_c);
}

test "GROSUB secondary-root metabolism preserves source equations" {
    const parameters: SecondaryRootParameters = .{
        .maximum_substrate_respiration_fraction_per_h = 0.015,
        .substrate_respiration_half_saturation_g_c_per_g_c = 0.025,
        .nitrogen_feedback_half_saturation_g_n_per_g_c = 0.1,
        .phosphorus_feedback_half_saturation_g_p_per_g_c = 0.01,
        .maintenance_respiration_g_c_per_g_n_h = 0.010,
        .nitrogen_assimilation_respiration_g_c_per_g_n = 1.70,
        .minimum_carbon_recycling_fraction = 0.167,
        .responsive_carbon_recycling_fraction = 0.333,
        .maximum_nitrogen_recycling_fraction = 0.667,
        .maximum_phosphorus_recycling_fraction = 0.667,
        .storage_exchange_fraction_per_h = 2.5e-5,
        .nonwoody_root_fraction_exponent = 0.167,
        .maintenance_gas_constant_j_per_mol_k = 8.3143,
        .maintenance_enthalpy_j_per_mol_k = 710,
        .maintenance_activation_energy_j_per_mol = 62500,
        .maintenance_low_temperature_inactivation_energy_j_per_mol = 197500,
        .maintenance_normalization_log_intercept = 25.216,
        .maximum_maintenance_temperature_response = 1.0e3,
        .shallow_root_water_response_per_mpa = 0.05,
        .deep_root_water_response_per_mpa = 0.10,
        .maintenance_water_response_exponent = 0.25,
        .root_penetration_reference_radius_m = 1.0e-3,
        .acidity_half_effect_hydrogen_activity_mol_per_m3 = 1,
        .maximum_acidity_enhancement = 4,
        .shallow_primary_root_sink_multiplier = 0.25,
        .intermediate_primary_root_sink_multiplier = 1,
        .deep_primary_root_sink_multiplier = 2,
        .deeper_primary_root_sink_multiplier = 4,
        .annual_termination_hours_without_grain_fill = 336,
        .root_protein_carbon_per_nitrogen_g_c_per_g_n = 2.5,
        .root_protein_carbon_per_phosphorus_g_c_per_g_p = 25,
        .nutrient_uptake_respiration_g_c_per_g_element = 0.86,
        .evergreen_leafoff_remobilization_start_fraction = 0.75,
        .deciduous_leafoff_remobilization_start_fraction = 0.5,
        .full_senescence_duration_h = 480,
    };
    const inputs: SecondaryRootInputs = .{
        .mobile_carbon_g_c = 0.2,
        .nonstructural_nitrogen_g_n = 0.04,
        .nonstructural_phosphorus_g_p = 0.004,
        .root_carbon_g_c = 2,
        .root_nitrogen_g_n = 0.04,
        .root_nitrogen_to_carbon_ratio_g_n_per_g_c = 0.02,
        .root_phosphorus_to_carbon_ratio_g_p_per_g_c = 0.002,
        .root_growth_yield_g_c_per_g_c = 0.8,
        .active_root_fraction = 0.5,
        .biological_timestep_h = 1,
        .substrate_temperature_response = 0.9,
        .maintenance_temperature_response = 0.8,
        .acidity_response = 0.75,
        .substrate_feedback = 0.6,
        .oxygen_limitation = 0.7,
        .substrate_water_response = 0.5,
        .maintenance_water_response = 0.5,
    };
    const result = try secondaryRootMetabolism(parameters, inputs);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0 / 3.0), result.nutrient_feedback, 1.0e-12);
    const rco2rm = 0.015 * 0.5 * 0.2 * 0.9 * (2.0 / 3.0) * 0.6 * 0.5 * 0.1 / 0.125;
    const rmncr = 0.010 * 0.04 * 0.8 * 0.75 * 0.5;
    try std.testing.expectApproxEqAbs(rco2rm, result.substrate_respiration_oxygen_unlimited_g_c_per_h, 1.0e-12);
    try std.testing.expectApproxEqAbs(rmncr, result.maintenance_respiration_g_c_per_h, 1.0e-12);
    try std.testing.expectApproxEqAbs(1.70 * result.nitrogen_growth_actual_g_n_per_h, result.nitrogen_assimilation_respiration_actual_g_c_per_h, 1.0e-12);
    var singular = inputs;
    singular.root_growth_yield_g_c_per_g_c = 1;
    try std.testing.expectError(error.InvalidSecondaryRootMetabolismInput, secondaryRootMetabolism(parameters, singular));
}

test "GROSUB secondary-root entry and respiration water selectors preserve source gates" {
    try std.testing.expect(secondaryRootAxisActive(2, 2, false));
    try std.testing.expect(!secondaryRootAxisActive(3, 2, false));
    try std.testing.expect(!secondaryRootAxisActive(2, 2, true));
    try std.testing.expect(rootRespirationActive(true, true));
    try std.testing.expect(rootRespirationActive(false, false));
    try std.testing.expect(!rootRespirationActive(false, true));

    const evergreen_deep = try sourceRootRespirationWaterResponses(2, 0, 0.25, 0.70);
    try std.testing.expectEqual(@as(f64, 0.25), evergreen_deep.substrate);
    try std.testing.expectEqual(@as(f64, 0.70), evergreen_deep.maintenance);
    const drought_deciduous = try sourceRootRespirationWaterResponses(2, 2, 0.25, 0.70);
    try std.testing.expectEqual(@as(f64, 0.25), drought_deciduous.substrate);
    try std.testing.expectEqual(@as(f64, 0.25), drought_deciduous.maintenance);
}

test "runtime root metabolism plant parameters reject invalid dimensions and codes" {
    const valid: RuntimePlantParameters = .{
        .root_profile_type = 2,
        .mycorrhizal_type = 2,
        .growth_habit = 1,
        .leaf_phenology_type = 0,
        .root_growth_yield_g_c_per_g_c = 0.7,
        .root_nitrogen_to_carbon_g_n_per_g_c = 0.03,
        .root_phosphorus_to_carbon_g_p_per_g_c = 0.004,
        .stalk_nitrogen_to_carbon_g_n_per_g_c = 0.01,
        .stalk_phosphorus_to_carbon_g_p_per_g_c = 0.001,
        .primary_root_radius_m = 0.001,
        .secondary_root_radius_m = 0.0002,
        .primary_specific_length_m_per_g_c = 10,
        .secondary_specific_length_m_per_g_c = 100,
        .secondary_root_branching_per_m = 20,
        .shoot_root_equilibration_fraction_per_h = 0.1,
    };
    try valid.validate();
    var invalid = valid;
    invalid.root_profile_type = 4;
    try std.testing.expectError(error.InvalidRootMetabolismPlantCode, invalid.validate());
    invalid = valid;
    invalid.secondary_specific_length_m_per_g_c = 0;
    try std.testing.expectError(error.InvalidRootMetabolismPlantParameter, invalid.validate());
}

test "STARTQ CNRTS and CPRTS yield-scaled ratios bound GROSUB root growth respiration" {
    const respiration_fraction = 0.3;
    const growth_yield = 0.8;
    const active_fraction = 0.5;
    const nitrogen_ratio = 0.02;
    const phosphorus_ratio = 0.002;
    const nitrogen = 0.04;
    const phosphorus = 0.003;
    const result = try nutrientLimitedRootGrowthRespiration(
        nitrogen,
        phosphorus,
        active_fraction,
        respiration_fraction,
        growth_yield,
        nitrogen_ratio,
        phosphorus_ratio,
    );
    const source_n_limit = nitrogen * respiration_fraction * active_fraction / (nitrogen_ratio * growth_yield);
    const source_p_limit = phosphorus * respiration_fraction * active_fraction / (phosphorus_ratio * growth_yield);
    try std.testing.expect(source_p_limit < source_n_limit);
    try std.testing.expectApproxEqAbs(source_p_limit, result, 1.0e-15);
    // Converting source growth respiration back through DMRTD and DMRT
    // consumes exactly the limiting active phosphorus inventory.
    const structural_growth_g_c = result / respiration_fraction * growth_yield;
    try std.testing.expectApproxEqAbs(phosphorus * active_fraction, structural_growth_g_c * phosphorus_ratio, 1.0e-15);
}

test "GROSUB primary-root bottom cap uses current axis maintenance demand" {
    const result = try primaryRootMetabolism(compatibilitySecondaryRootParameters(), .{
        .shared = .{
            .mobile_carbon_g_c = 1,
            .nonstructural_nitrogen_g_n = 0.1,
            .nonstructural_phosphorus_g_p = 0.01,
            .root_carbon_g_c = 2,
            .root_nitrogen_g_n = 0.02,
            .root_nitrogen_to_carbon_ratio_g_n_per_g_c = 0.02,
            .root_phosphorus_to_carbon_ratio_g_p_per_g_c = 0.002,
            .root_growth_yield_g_c_per_g_c = 0.8,
            .active_root_fraction = 1,
            .biological_timestep_h = 1,
            .substrate_temperature_response = 1,
            .maintenance_temperature_response = 1,
            .acidity_response = 1,
            .substrate_feedback = 1,
            .oxygen_limitation = 1,
            .substrate_water_response = 1,
            .maintenance_water_response = 1,
        },
        .primary_tip_at_or_below_profile_bottom = true,
    });
    try std.testing.expectApproxEqAbs(result.maintenance_respiration_g_c_per_h, result.substrate_respiration_oxygen_unlimited_g_c_per_h, 1e-15);
    try std.testing.expectEqual(@as(f64, 0), result.growth_respiration_actual_g_c_per_h);
}

test "GROSUB primary-root commit updates axis and shared pools atomically" {
    var roots = try RootState.init(std.testing.allocator, 1, 1, 1);
    defer roots.deinit();
    roots.mobile_carbon_g[0] = 10;
    roots.mobile_nitrogen_g[0] = 1;
    roots.mobile_phosphorus_g[0] = 0.1;
    roots.axis_primary_carbon_g[0] = 2;
    roots.axis_primary_nitrogen_g[0] = 0.2;
    roots.axis_primary_phosphorus_g[0] = 0.02;
    roots.axis_primary_length_m[0] = 4;
    roots.axis_depth_m[0] = 1;
    const senescence: SecondaryRootSenescence = .{
        .respiration_oxygen_unlimited_g_c_per_h = 0.1,
        .respiration_actual_g_c_per_h = 0.1,
        .phenological_senescence_g_c_per_h = 0,
        .senesced_fraction = 0.25,
        .recyclable_carbon_g_c = 0.8,
        .recyclable_nitrogen_g_n = 0.08,
        .recyclable_phosphorus_g_p = 0.008,
    };
    try commitPrimaryRoot(&roots, 0, 0, 0, .{
        .metabolism = std.mem.zeroes(SecondaryRootResult),
        .senescence = senescence,
        .primary_specific_length_m_per_g_c = 10,
        .root_extension_water_response = 1,
        .nonwoody_carbon_fraction = 0.6,
        .nonwoody_nitrogen_fraction = 0.7,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 0,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 0,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 10.02), roots.mobile_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.014), roots.mobile_nitrogen_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.102), roots.mobile_phosphorus_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), roots.axis_primary_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), roots.axis_primary_length_m[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), roots.axis_depth_m[0], 1e-12);
    const mobile_before = roots.mobile_carbon_g[0];
    try std.testing.expectError(error.PlantRootIndexOutOfBounds, commitPrimaryRoot(&roots, 0, 0, 2, .{
        .metabolism = std.mem.zeroes(SecondaryRootResult),
        .senescence = senescence,
        .primary_specific_length_m_per_g_c = 10,
        .root_extension_water_response = 1,
        .nonwoody_carbon_fraction = 0.6,
        .nonwoody_nitrogen_fraction = 0.7,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 0,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 0,
    }));
    try std.testing.expectEqual(mobile_before, roots.mobile_carbon_g[0]);
}

test "GROSUB primary-root respiration is allocated across traversed layers" {
    var roots = try RootState.init(std.testing.allocator, 1, 3, 1);
    defer roots.deinit();
    const indices = [_]usize{
        try roots.layerIndex(0, 0, 0),
        try roots.layerIndex(0, 0, 1),
        try roots.layerIndex(0, 0, 2),
    };
    var fractions: [3]f64 = undefined;
    try allocatePrimaryRootRespiration(&roots, &indices, &.{ 0.2, 0.3, 0.1 }, &fractions, 1.1, 0.1, true, .{
        .actual_g_c = 2,
        .oxygen_unlimited_g_c = 3,
        .carbon_unlimited_g_c = 4,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), fractions[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), fractions[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), fractions[2], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), roots.actual_respiration_g_c_per_h[indices[0]], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), roots.actual_respiration_g_c_per_h[indices[2]], 1e-15);
    var actual_total: f64 = 0;
    for (indices) |root| actual_total += roots.actual_respiration_g_c_per_h[root];
    try std.testing.expectApproxEqAbs(@as(f64, 2), actual_total, 1e-15);
}

test "GROSUB primary root depth retracts under negative net growth" {
    const positive = try primaryRootLengthChange(0.2, 0.2, 4, 1.1, 0.1, 10, 2, 0.5);
    try std.testing.expectApproxEqAbs(0.5, positive, 1e-15);
    const retraction = try primaryRootLengthChange(0.02, -0.8, 4, 1.1, 0.1, 10, 2, 0.5);
    // Gross extension is 0.05 m; proportional withdrawal is -0.20 m.
    try std.testing.expectApproxEqAbs(-0.15, retraction, 1e-15);
}

test "GROSUB secondary-root recycling and senescence preserve source branches" {
    const parameters = compatibilitySecondaryRootParameters();
    const recycling = try secondaryRootRecyclingFractions(true, 0.2, 0.04, 0.004, parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 0.167 + 0.333 * (2.0 / 3.0)), recycling.carbon, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.667 / 3.0), recycling.nitrogen, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.667 / 3.0), recycling.phosphorus, 1.0e-12);
    const senescence = try secondaryRootSenescence(.{
        .oxygen_unlimited_substrate_minus_maintenance_g_c_per_h = -0.3,
        .actual_substrate_minus_maintenance_g_c_per_h = -0.4,
        .root_carbon_g_c = 1,
        .root_nitrogen_g_n = 0.02,
        .root_phosphorus_g_p = 0.002,
        .oxygen_limitation = 0.5,
        .phenological_remobilization_enabled = true,
        .root_remobilization_enabled = true,
        .storage_exchange_fraction_per_h = 0.01,
        .remobilization_elapsed_h = 50,
        .full_senescence_h = 100,
        .biological_timestep_h = 1,
        .structural_presence_threshold_g_c = 1.0e-12,
    }, recycling);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), senescence.respiration_oxygen_unlimited_g_c_per_h, 1.0e-12);
    try std.testing.expectApproxEqAbs(senescence.recyclable_carbon_g_c * 0.5 + 0.005, senescence.respiration_actual_g_c_per_h, 1.0e-12);
    try std.testing.expectApproxEqAbs(senescence.respiration_actual_g_c_per_h / senescence.recyclable_carbon_g_c, senescence.senesced_fraction, 1.0e-12);
}

test "GROSUB primary senescence excludes secondary phenological remobilization" {
    const recycling: RecyclingFractions = .{ .carbon = 0.5, .nitrogen = 0.5, .phosphorus = 0.5 };
    const inputs: SecondaryRootSenescenceInputs = .{
        .oxygen_unlimited_substrate_minus_maintenance_g_c_per_h = 0,
        .actual_substrate_minus_maintenance_g_c_per_h = 0,
        .root_carbon_g_c = 10,
        .root_nitrogen_g_n = 1,
        .root_phosphorus_g_p = 0.1,
        .oxygen_limitation = 1,
        .phenological_remobilization_enabled = true,
        .root_remobilization_enabled = true,
        .storage_exchange_fraction_per_h = 0.1,
        .remobilization_elapsed_h = 10,
        .full_senescence_h = 10,
        .biological_timestep_h = 1,
        .structural_presence_threshold_g_c = 1e-12,
    };
    const secondary = try secondaryRootSenescence(inputs, recycling);
    const primary = try primaryRootSenescence(inputs, recycling);
    try std.testing.expectEqual(@as(f64, 1), secondary.phenological_senescence_g_c_per_h);
    try std.testing.expectEqual(@as(f64, 0), primary.phenological_senescence_g_c_per_h);
    try std.testing.expectEqual(@as(f64, 0), primary.respiration_actual_g_c_per_h);
}

test "GROSUB root wood composition preserves FWODR and weighted growth ratios" {
    const composition = try rootWoodComposition(true, true, 8, 2, 0.01, 0.03, 0.001, 0.003, 1.0e-12, 0.167);
    const nonwoody = std.math.pow(f64, 0.25, 0.167);
    try std.testing.expectApproxEqAbs(nonwoody, composition.carbon_fraction[1], 1.0e-12);
    try std.testing.expectEqual(composition.carbon_fraction, composition.nitrogen_fraction);
    try std.testing.expectApproxEqAbs((1 - nonwoody) * 0.01 + nonwoody * 0.03, composition.growth_nitrogen_to_carbon_g_n_per_g_c, 1.0e-12);
    const herbaceous = try rootWoodComposition(false, false, 0, 0, 0.01, 0.03, 0.001, 0.003, 1.0e-12, 0.167);
    try std.testing.expectEqual([2]f64{ 0, 1 }, herbaceous.carbon_fraction);
}

test "GROSUB root environment preserves TFN6 FPH and water responses" {
    const parameters = compatibilitySecondaryRootParameters();
    const response = try rootEnvironmentResponses(parameters, 298.15, 0, 7, -1, 0.5, 0, 0.2, 0.5e-3, true);
    const adjusted_temperature_k = 298.15;
    const rtk = 8.3143 * adjusted_temperature_k;
    const expected_temperature = @min(@as(f64, 1.0e3), std.math.exp(25.216 - @as(f64, 62500) / rtk) / (1 + std.math.exp((197500 - 710 * adjusted_temperature_k) / rtk)));
    try std.testing.expectApproxEqAbs(expected_temperature, response.maintenance_temperature, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0001), response.acidity, 1.0e-12);
    try std.testing.expectApproxEqAbs(std.math.exp(@as(f64, -0.05)), response.growth_water, 1.0e-12);
    try std.testing.expectApproxEqAbs(std.math.pow(f64, response.growth_water, 0.25), response.maintenance_water, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), response.scaled_penetration_resistance_mpa, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.45), response.extension_water, 1.0e-12);
}

test "GROSUB next lower root layer skips thin layers but retains the bottom" {
    const thickness_m = [_]f64{ 0.1, 1e-8, 0.2, 0 };
    try std.testing.expectEqual(@as(usize, 2), try nextLowerRootLayer(&thickness_m, 0, 1e-6));
    try std.testing.expectEqual(@as(usize, 3), try nextLowerRootLayer(&thickness_m, 2, 1e-6));
    try std.testing.expectError(error.NoLowerRootLayer, nextLowerRootLayer(&thickness_m, 3, 1e-6));
    try std.testing.expectError(error.InvalidRootLayerThickness, nextLowerRootLayer(&.{ 0.1, std.math.nan(f64) }, 0, 1e-6));
}

test "GROSUB root axis sink strengths retain series equation and runtime axes" {
    const parameters = compatibilitySecondaryRootParameters();
    const first = try rootAxisSinkStrength(parameters, .{
        .root_profile_type = 2,
        .primary_axis_count_multiplier = 3,
        .primary_root_radius_m = 2e-3,
        .primary_root_depth_from_canopy_m = 0.4,
        .secondary_root_depth_from_canopy_m = 0.3,
        .secondary_axis_count = 8,
        .secondary_root_radius_m = 1e-3,
        .average_secondary_root_length_m = 0.2,
        .primary_biological_domain = true,
    });
    const primary_series = 3 * std.math.pow(f64, 2e-3, 2) / 0.3;
    const secondary_parallel = 8 * std.math.pow(f64, 1e-3, 2) / 0.2;
    try std.testing.expectApproxEqAbs(2 * 3 * std.math.pow(f64, 2e-3, 2) / 0.4, first.primary_m, 1.0e-15);
    try std.testing.expectApproxEqAbs(primary_series * secondary_parallel / (primary_series + secondary_parallel), first.secondary_m, 1.0e-15);
    var strengths = [_]RootAxisSinkStrength{first} ** 7;
    strengths[6] = .{ .primary_m = 0, .secondary_m = first.secondary_m };
    var primary_fractions: [7]f64 = undefined;
    var secondary_fractions: [7]f64 = undefined;
    const total = try normalizeRootAxisSinkFractions(&strengths, &primary_fractions, &secondary_fractions, 1e-20);
    try std.testing.expect(total > 0);
    var fraction_sum: f64 = 0;
    for (primary_fractions, secondary_fractions) |primary, secondary| fraction_sum += primary + secondary;
    try std.testing.expectApproxEqAbs(@as(f64, 1), fraction_sum, 1.0e-12);
}

test "GROSUB source sink comparator gates primary tips and uses rooted midpoint depth" {
    const parameters = compatibilitySecondaryRootParameters();
    const shared: SourceOrderRootAxisSinkInputs = .{
        .root_profile_type = 2,
        .primary_axis_count_multiplier = 3,
        .primary_root_radius_m = 2e-3,
        .primary_root_depth_from_surface_m = 0.45,
        .layer_top_depth_m = 0.2,
        .layer_thickness_m = 0.2,
        .secondary_root_origin_offset_m = 0.05,
        .seeding_depth_m = 0.1,
        .hypocotyledon_height_m = 0.02,
        .canopy_height_m = 0.3,
        .secondary_axis_count = 8,
        .secondary_root_radius_m = 1e-3,
        .average_secondary_root_length_m = 0.2,
        .negligible_sink_m = 1e-20,
        .primary_biological_domain = true,
    };
    const outside_tip_layer = try sourceOrderRootAxisSinkStrength(parameters, shared);
    try std.testing.expectEqual(@as(f64, 0), outside_tip_layer.primary_m);
    const rooted_length_m = 0.2;
    const secondary_depth_m = 0.2 + 0.5 * rooted_length_m + 0.3;
    const primary_series = 3 * std.math.pow(f64, 2e-3, 2) / secondary_depth_m;
    const secondary_parallel = 8 * std.math.pow(f64, 1e-3, 2) / 0.2;
    try std.testing.expectApproxEqAbs(primary_series * secondary_parallel / (primary_series + secondary_parallel), outside_tip_layer.secondary_m, 1e-15);

    var mycorrhiza = shared;
    mycorrhiza.primary_biological_domain = false;
    const mycorrhizal = try sourceOrderRootAxisSinkStrength(parameters, mycorrhiza);
    try std.testing.expectEqual(@as(f64, 0), mycorrhizal.primary_m);
    try std.testing.expectApproxEqAbs(secondary_parallel, mycorrhizal.secondary_m, 1e-15);
}

test "root metabolism grid workspace is runtime sized and cell independent" {
    var workspace = try GridWorkspace.init(std.testing.allocator, 4, 19, 7);
    defer workspace.deinit();
    try std.testing.expectEqual(@as(usize, 4), workspace.per_cell.len);
    for (workspace.per_cell) |cell| {
        try std.testing.expectEqual(@as(usize, 19), cell.sink_strengths.len);
        try std.testing.expectEqual(@as(usize, 7), cell.primary_respiration_allocation_fractions.len);
    }
    workspace.per_cell[0].primary_active[18] = true;
    try workspace.per_cell[0].beginPlantHour(19);
    try workspace.per_cell[0].markPrimaryProcessed(0, 18);
    try workspace.per_cell[0].resetAxes(19);
    try std.testing.expect(!workspace.per_cell[0].primary_active[18]);
    try std.testing.expect(try workspace.per_cell[0].primaryWasProcessed(0, 18));
    try workspace.per_cell[0].beginPlantHour(19);
    try std.testing.expect(!try workspace.per_cell[0].primaryWasProcessed(0, 18));
    try std.testing.expect(workspace.per_cell[0].sink_strengths.ptr != workspace.per_cell[1].sink_strengths.ptr);
    try std.testing.expectError(error.RootMetabolismWorkspaceCapacityExceeded, workspace.per_cell[0].resetAxes(20));
}

test "staged primary and secondary axes commit one shared mobile-pool transaction" {
    var roots = try RootState.init(std.testing.allocator, 1, 1, 2);
    defer roots.deinit();
    var workspace = try AxisWorkspace.init(std.testing.allocator, 2, 1);
    defer workspace.deinit();
    try workspace.resetAxes(2);
    roots.mobile_carbon_g[0] = 10;
    roots.mobile_nitrogen_g[0] = 1;
    roots.mobile_phosphorus_g[0] = 0.1;
    roots.axis_primary_carbon_g[0] = 1;
    roots.axis_secondary_carbon_g[1] = 1;
    workspace.primary_active[0] = true;
    workspace.secondary_active[1] = true;
    workspace.primary_metabolism[0].root_growth_actual_g_c_per_h = 0.2;
    workspace.primary_metabolism[0].growth_and_respiration_carbon_actual_g_c_per_h = 0.25;
    workspace.primary_metabolism[0].nitrogen_growth_actual_g_n_per_h = 0.02;
    workspace.primary_metabolism[0].phosphorus_growth_actual_g_p_per_h = 0.002;
    workspace.secondary_metabolism[1].root_growth_actual_g_c_per_h = 0.3;
    workspace.secondary_metabolism[1].growth_and_respiration_carbon_actual_g_c_per_h = 0.375;
    workspace.secondary_metabolism[1].nitrogen_growth_actual_g_n_per_h = 0.03;
    workspace.secondary_metabolism[1].phosphorus_growth_actual_g_p_per_h = 0.003;
    const parameters: StagedLayerCommitParameters = .{
        .primary_specific_length_m_per_g_c = 10,
        .secondary_specific_length_m_per_g_c = 20,
        .plant_population_count = 1,
        .seeding_depth_m = 0,
        .current_layer_bottom_depth_m = 2,
        .next_layer_thickness_m = 0,
        .extension_presence_threshold_m = 0,
        .root_extension_water_response = 0.5,
        .nonwoody_carbon_fraction = 1,
        .nonwoody_nitrogen_fraction = 1,
        .nonwoody_phosphorus_fraction = 1,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 0,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 0,
    };
    try commitStagedLayerAxes(&roots, 0, 0, 0, &workspace, 2, parameters);
    try std.testing.expectApproxEqAbs(@as(f64, 9.375), roots.mobile_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.95), roots.mobile_nitrogen_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.095), roots.mobile_phosphorus_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), roots.axis_primary_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.3), roots.axis_secondary_carbon_g[1], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), roots.axis_primary_length_m[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 3), roots.axis_secondary_length_m[1], 1e-15);

    const mobile_before = roots.mobile_carbon_g[0];
    const primary_before = roots.axis_primary_carbon_g[0];
    try workspace.resetAxes(2);
    workspace.primary_active[0] = true;
    workspace.primary_metabolism[0].growth_and_respiration_carbon_actual_g_c_per_h = 20;
    try std.testing.expectError(error.StagedRootCommitWouldOverdrawPool, commitStagedLayerAxes(&roots, 0, 0, 0, &workspace, 2, parameters));
    try std.testing.expectEqual(mobile_before, roots.mobile_carbon_g[0]);
    try std.testing.expectEqual(primary_before, roots.axis_primary_carbon_g[0]);
}

test "GROSUB primary crossing routes growth and mobile pools to the next runtime layer" {
    const placement = try primaryRootExtensionPlacement(2, 0.9, 1.0, 0.5);
    try std.testing.expect(placement.crosses_into_next_layer);
    try std.testing.expectEqual(@as(f64, 0.5), placement.extension_m);

    var roots = try RootState.init(std.testing.allocator, 1, 2, 1);
    defer roots.deinit();
    var workspace = try AxisWorkspace.init(std.testing.allocator, 1, 2);
    defer workspace.deinit();
    try workspace.resetAxes(1);
    const current_root = try roots.layerIndex(0, 0, 0);
    const next_root = try roots.layerIndex(0, 0, 1);
    const current_axis = try roots.layerAxisIndex(0, 0, 0, 0);
    const next_axis = try roots.layerAxisIndex(0, 0, 1, 0);
    roots.mobile_carbon_g[current_root] = 8;
    roots.mobile_nitrogen_g[current_root] = 1;
    roots.mobile_phosphorus_g[current_root] = 0.1;
    roots.total_water_potential_mpa[current_root] = -0.4;
    roots.osmotic_water_potential_mpa[current_root] = -0.8;
    roots.turgor_water_potential_mpa[current_root] = 0.4;
    roots.primary_radius_m[current_root] = 0.001;
    roots.axis_primary_carbon_g[current_axis] = 1;
    roots.axis_primary_nitrogen_g[current_axis] = 0.1;
    roots.axis_primary_phosphorus_g[current_axis] = 0.01;
    roots.axis_depth_m[try roots.axisIndex(0, 0, 0)] = 0.9;
    workspace.primary_active[0] = true;
    workspace.primary_sink_fractions[0] = 0.25;
    workspace.primary_metabolism[0].root_growth_actual_g_c_per_h = 0.2;
    workspace.primary_metabolism[0].nitrogen_growth_actual_g_n_per_h = 0.02;
    workspace.primary_metabolism[0].phosphorus_growth_actual_g_p_per_h = 0.002;
    try commitStagedLayerAxes(&roots, 0, 0, 0, &workspace, 1, .{
        .primary_specific_length_m_per_g_c = 10,
        .secondary_specific_length_m_per_g_c = 20,
        .plant_population_count = 1,
        .seeding_depth_m = 0,
        .current_layer_bottom_depth_m = 1,
        .next_layer_thickness_m = 0.5,
        .extension_presence_threshold_m = 0,
        .root_extension_water_response = 1,
        .nonwoody_carbon_fraction = 1,
        .nonwoody_nitrogen_fraction = 1,
        .nonwoody_phosphorus_fraction = 1,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 0,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 0,
    });
    try std.testing.expectEqual(@as(f64, 1), roots.axis_primary_carbon_g[current_axis]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), roots.axis_primary_carbon_g[next_axis], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), roots.axis_primary_length_m[next_axis], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1.4), roots.axis_depth_m[try roots.axisIndex(0, 0, 0)], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 6), roots.mobile_carbon_g[current_root], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 2), roots.mobile_carbon_g[next_root], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.735), roots.mobile_nitrogen_g[current_root], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.245), roots.mobile_nitrogen_g[next_root], 1e-15);
    try std.testing.expectEqual(roots.total_water_potential_mpa[current_root], roots.total_water_potential_mpa[next_root]);
    try std.testing.expectEqual(roots.primary_radius_m[current_root], roots.primary_radius_m[next_root]);
}

test "GROSUB primary crossing requires extension above ZEROP" {
    const below = try sourceOrderPrimaryRootExtensionPlacement(5e-7, 0.9999998, 1, 0.5, 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 5e-7), below.extension_m, 1e-18);
    try std.testing.expect(!below.crosses_into_next_layer);
    const above = try sourceOrderPrimaryRootExtensionPlacement(2e-6, 0.9999998, 1, 0.5, 1e-6);
    try std.testing.expect(above.crosses_into_next_layer);
}

test "GROSUB negative primary growth consumes secondary roots in tip then upper layer" {
    const result = try absorbPrimaryDeficitFromSecondaryRoots(
        3.5,
        0.35,
        0.035,
        .{ .carbon_g_c = 2, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.02, .length_m = 4 },
        .{ .carbon_g_c = 3, .nitrogen_g_n = 0.3, .phosphorus_g_p = 0.03, .length_m = 6 },
    );
    try std.testing.expectEqual(@as(f64, 0), result.current.carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), result.current.length_m);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), result.upper.carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 3), result.upper.length_m, 1e-15);
    try std.testing.expectEqual(@as(f64, 0), result.residual_carbon_deficit_g_c);
    try std.testing.expectEqual(@as(f64, 0), result.residual_nitrogen_deficit_g_n);
    try std.testing.expectEqual(@as(f64, 0), result.residual_phosphorus_deficit_g_p);

    const exhausted = try absorbPrimaryDeficitFromSecondaryRoots(
        7,
        0.7,
        0.07,
        .{ .carbon_g_c = 2, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.02, .length_m = 4 },
        .{ .carbon_g_c = 3, .nitrogen_g_n = 0.3, .phosphorus_g_p = 0.03, .length_m = 6 },
    );
    try std.testing.expectEqual(@as(f64, 2), exhausted.residual_carbon_deficit_g_c);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), exhausted.residual_nitrogen_deficit_g_n, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.02), exhausted.residual_phosphorus_deficit_g_p, 1e-15);
}

test "live staged GROSUB commit absorbs primary senescence from secondary layers first" {
    var roots = try RootState.init(std.testing.allocator, 1, 2, 1);
    defer roots.deinit();
    var workspace = try AxisWorkspace.init(std.testing.allocator, 1, 2);
    defer workspace.deinit();
    try workspace.resetAxes(1);
    const upper_axis = try roots.layerAxisIndex(0, 0, 0, 0);
    const tip_axis = try roots.layerAxisIndex(0, 0, 1, 0);
    roots.axis_primary_carbon_g[tip_axis] = 2;
    roots.axis_primary_nitrogen_g[tip_axis] = 0.2;
    roots.axis_primary_phosphorus_g[tip_axis] = 0.02;
    roots.axis_secondary_carbon_g[tip_axis] = 0.6;
    roots.axis_secondary_nitrogen_g[tip_axis] = 0.06;
    roots.axis_secondary_phosphorus_g[tip_axis] = 0.006;
    roots.axis_secondary_length_m[tip_axis] = 3;
    roots.axis_secondary_carbon_g[upper_axis] = 0.6;
    roots.axis_secondary_nitrogen_g[upper_axis] = 0.06;
    roots.axis_secondary_phosphorus_g[upper_axis] = 0.006;
    roots.axis_secondary_length_m[upper_axis] = 3;
    roots.axis_depth_m[try roots.axisIndex(0, 0, 0)] = 1.5;
    workspace.primary_active[0] = true;
    workspace.primary_senescence[0].senesced_fraction = 0.5;
    try commitStagedLayerAxes(&roots, 0, 0, 1, &workspace, 1, .{
        .primary_specific_length_m_per_g_c = 10,
        .secondary_specific_length_m_per_g_c = 20,
        .plant_population_count = 1,
        .seeding_depth_m = 0,
        .current_layer_bottom_depth_m = 2,
        .next_layer_thickness_m = 0,
        .extension_presence_threshold_m = 0,
        .root_extension_water_response = 1,
        .nonwoody_carbon_fraction = 1,
        .nonwoody_nitrogen_fraction = 1,
        .nonwoody_phosphorus_fraction = 1,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 0,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 0,
    });
    try std.testing.expectEqual(@as(f64, 2), roots.axis_primary_carbon_g[tip_axis]);
    try std.testing.expectEqual(@as(f64, 0), roots.axis_secondary_carbon_g[tip_axis]);
    try std.testing.expectEqual(@as(f64, 0), roots.axis_secondary_length_m[tip_axis]);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), roots.axis_secondary_carbon_g[upper_axis], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), roots.axis_secondary_length_m[upper_axis], 1e-15);
    try std.testing.expect(workspace.primary_deficit_active[0]);
}

test "STOMATE annual termination feedback is shared with root metabolism" {
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), try annualTerminationFeedback(0, 168, 336), 1.0e-12);
    try std.testing.expectEqual(@as(f64, 1), try annualTerminationFeedback(1, 168, 336));
    try std.testing.expectEqual(@as(f64, 0), try annualTerminationFeedback(0, 400, 336));
}

test "GROSUB secondary-root litter allocation and commit are atomic" {
    const senescence: SecondaryRootSenescence = .{
        .respiration_oxygen_unlimited_g_c_per_h = 0.1,
        .respiration_actual_g_c_per_h = 0.1,
        .phenological_senescence_g_c_per_h = 0,
        .senesced_fraction = 0.25,
        .recyclable_carbon_g_c = 0.8,
        .recyclable_nitrogen_g_n = 0.08,
        .recyclable_phosphorus_g_p = 0.008,
    };
    const quarter = [_]f64{0.25} ** 4;
    const litter = try secondaryRootLitter(senescence, 2, 0.2, 0.02, .{ 0.4, 0.6 }, .{ 0.3, 0.7 }, .{ 0.2, 0.8 }, .{
        .woody_carbon = quarter,
        .woody_nitrogen = quarter,
        .woody_phosphorus = quarter,
        .nonwoody_carbon = quarter,
        .nonwoody_nitrogen = quarter,
        .nonwoody_phosphorus = quarter,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), litter.woody_carbon_g_c[0], 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.045), litter.nonwoody_carbon_g_c[0], 1.0e-12);

    var roots = try RootState.init(std.testing.allocator, 1, 1, 1);
    defer roots.deinit();
    roots.mobile_carbon_g[0] = 10;
    roots.mobile_nitrogen_g[0] = 1;
    roots.mobile_phosphorus_g[0] = 0.1;
    roots.axis_secondary_carbon_g[0] = 2;
    roots.axis_secondary_nitrogen_g[0] = 0.2;
    roots.axis_secondary_phosphorus_g[0] = 0.02;
    roots.axis_secondary_length_m[0] = 4;
    const metabolism = try secondaryRootMetabolism(compatibilitySecondaryRootParameters(), .{
        .mobile_carbon_g_c = 0.2,
        .nonstructural_nitrogen_g_n = 0.04,
        .nonstructural_phosphorus_g_p = 0.004,
        .root_carbon_g_c = 2,
        .root_nitrogen_g_n = 0.04,
        .root_nitrogen_to_carbon_ratio_g_n_per_g_c = 0.02,
        .root_phosphorus_to_carbon_ratio_g_p_per_g_c = 0.002,
        .root_growth_yield_g_c_per_g_c = 0.8,
        .active_root_fraction = 0.5,
        .biological_timestep_h = 1,
        .substrate_temperature_response = 0.9,
        .maintenance_temperature_response = 0.8,
        .acidity_response = 0.75,
        .substrate_feedback = 0.6,
        .oxygen_limitation = 0.7,
        .substrate_water_response = 0.5,
        .maintenance_water_response = 0.5,
    });
    const commit_inputs: SecondaryRootCommitInputs = .{
        .metabolism = metabolism,
        .senescence = senescence,
        .root_specific_length_m_per_g_c = 10,
        .root_extension_water_response = 0.8,
        .nonwoody_carbon_fraction = 0.6,
        .nonwoody_nitrogen_fraction = 0.7,
        .nonwoody_phosphorus_fraction = 0.8,
        .protein_carbon_per_nitrogen_g_c_per_g_n = 2,
        .protein_carbon_per_phosphorus_g_c_per_g_p = 20,
    };
    try commitSecondaryRoot(&roots, 0, 0, commit_inputs);
    try std.testing.expect(roots.axis_secondary_carbon_g[0] < 2);
    try std.testing.expect(roots.actual_respiration_g_c_per_h[0] > 0);
    const mobile_before_failure = roots.mobile_carbon_g[0];
    roots.mobile_carbon_g[0] = 0;
    var failing_inputs = commit_inputs;
    failing_inputs.senescence.recyclable_carbon_g_c = 0;
    failing_inputs.senescence.recyclable_nitrogen_g_n = 0;
    failing_inputs.senescence.recyclable_phosphorus_g_p = 0;
    try std.testing.expectError(error.SecondaryRootCommitWouldOverdrawPool, commitSecondaryRoot(&roots, 0, 0, failing_inputs));
    try std.testing.expectEqual(@as(f64, 0), roots.mobile_carbon_g[0]);
    try std.testing.expect(mobile_before_failure > 0);
}

test "GROSUB nutrient uptake respiration retains 0.86 coefficient and three limits" {
    const result = NutrientResult{ .demand_g_element = 1, .uptake_g_element = 0.5, .oxygen_unlimited_uptake_g_element = 0.75, .carbon_unlimited_uptake_g_element = 1, .available_g_element = 2 };
    const respiration = try nutrientUptakeRespiration(&([_]NutrientResult{result} ** 8), 0.86);
    try std.testing.expectApproxEqAbs(@as(f64, 3.44), respiration.actual_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5.16), respiration.oxygen_unlimited_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 6.88), respiration.carbon_unlimited_g_c, 1.0e-12);
}

test "root respiration commit is conservative and rollback safe" {
    var roots = try RootState.init(std.testing.allocator, 1, 1, 1);
    defer roots.deinit();
    roots.mobile_carbon_g[0] = 5;
    try commit(&roots, 0, .{ .actual_g_c = 2, .oxygen_unlimited_g_c = 3, .carbon_unlimited_g_c = 4 });
    try std.testing.expectEqual(@as(f64, 3), roots.mobile_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 2), roots.actual_respiration_g_c_per_h[0]);
    try std.testing.expectError(error.InsufficientRootMobileCarbonForRespiration, commit(&roots, 0, .{ .actual_g_c = 4, .oxygen_unlimited_g_c = 4, .carbon_unlimited_g_c = 4 }));
    try std.testing.expectEqual(@as(f64, 3), roots.mobile_carbon_g[0]);
}
test "GROSUB primary-root axis scaling retains the larger carbon basis" {
    const decay_limited = try primaryRootAxisScaling(4, 2, 1, 1);
    const expected_retained = 0.999992087 * 4.0;
    try std.testing.expectApproxEqAbs(expected_retained, decay_limited.retained_root_carbon_g_c_per_plant, 1.0e-15);
    try std.testing.expectApproxEqAbs(std.math.pow(f64, expected_retained, 0.667), decay_limited.primary_axis_count_multiplier, 1.0e-15);

    const biomass_limited = try primaryRootAxisScaling(1, 18, 3, 1);
    try std.testing.expectEqual(@as(f64, 6), biomass_limited.retained_root_carbon_g_c_per_plant);
    try std.testing.expectApproxEqAbs(std.math.pow(f64, 6, 0.667) * 3, biomass_limited.primary_axis_count_multiplier, 1.0e-15);
}

test "GROSUB primary-root axis scaling clears state without population" {
    const result = try primaryRootAxisScaling(12, 30, 0, 1);
    try std.testing.expectEqual(@as(f64, 0), result.retained_root_carbon_g_c_per_plant);
    try std.testing.expectEqual(@as(f64, 0), result.primary_axis_count_multiplier);
}

test "GROSUB primary-root axis scaling rejects invalid state" {
    try std.testing.expectError(error.InvalidPrimaryRootAxisScalingInput, primaryRootAxisScaling(-1, 1, 1, 1));
    try std.testing.expectError(error.InvalidPrimaryRootAxisScalingInput, primaryRootAxisScaling(1, 1, 1, 0));
    try std.testing.expectError(error.InvalidPrimaryRootAxisScalingInput, primaryRootAxisScaling(1, std.math.nan(f64), 1, 1));
}
