const std = @import("std");
const PlantTraits = @import("plant_traits.zig").PlantTraits;

pub const biological_domain_count: usize = 2;
pub const salt_species_count: usize = 8;
pub const organic_substrate_count: usize = 5;
pub const transported_root_gas_count: usize = 6;

pub const InitializationParameters = struct {
    root_nitrogen_to_maximum_protein_multiplier: f64,
    root_phosphorus_to_maximum_protein_multiplier: f64,
    mycorrhizal_radius_m: f64,
    initial_total_water_potential_mpa: f64,
    osmotic_water_potential_decrement_mpa: f64,
    initial_active_length_m: f64,
    initial_water_fraction: f64,

    pub fn validate(self: InitializationParameters) !void {
        inline for (@typeInfo(InitializationParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(self, field.name))) return error.NonFiniteRootInitializationParameter;
        if (self.root_nitrogen_to_maximum_protein_multiplier <= 0 or self.root_phosphorus_to_maximum_protein_multiplier <= 0 or self.mycorrhizal_radius_m <= 0 or self.initial_total_water_potential_mpa > 0 or self.osmotic_water_potential_decrement_mpa > 0 or self.initial_active_length_m < 0 or self.initial_water_fraction < 0 or self.initial_water_fraction > 1) return error.InvalidRootInitializationParameter;
    }
};

pub fn compatibilityInitializationParameters() InitializationParameters {
    return .{
        .root_nitrogen_to_maximum_protein_multiplier = 2.5,
        .root_phosphorus_to_maximum_protein_multiplier = 25.0,
        .mycorrhizal_radius_m = 2.5e-6,
        .initial_total_water_potential_mpa = -0.01,
        .osmotic_water_potential_decrement_mpa = -0.01,
        .initial_active_length_m = 1.0e-3,
        .initial_water_fraction = 1.0,
    };
}

pub const MorphologyParameters = struct {
    minimum_average_secondary_length_m: f64,
    root_elastic_modulus_mpa: f64,

    pub fn validate(self: MorphologyParameters) !void {
        inline for (@typeInfo(MorphologyParameters).@"struct".fields) |field| if (!std.math.isFinite(@field(self, field.name))) return error.NonFiniteRootMorphologyParameter;
        if (self.minimum_average_secondary_length_m <= 0 or self.root_elastic_modulus_mpa <= 0) return error.InvalidRootMorphologyParameter;
    }
};

pub fn compatibilityMorphologyParameters() MorphologyParameters {
    return .{
        .minimum_average_secondary_length_m = 1.0e-2,
        .root_elastic_modulus_mpa = 5.0,
    };
}
pub const SaltSpecies = enum(u8) { aluminum, iron, calcium, magnesium, sodium, potassium, sulfate, chloride };

pub const UptakeCompetition = struct {
    oxygen: f64,
    ammonium_nonband: f64,
    ammonium_band: f64,
    nitrate_nonband: f64,
    nitrate_band: f64,
    phosphate_h2_nonband: f64,
    phosphate_h2_band: f64,
    phosphate_h_nonband: f64,
    phosphate_h_band: f64,
};

/// UPTAKE competition fractions for one root or mycorrhizal population.
/// Previous-hour population demand is used when the corresponding combined
/// microbial/root demand is significant; otherwise the biome root fraction is
/// retained. The source minimum population share is applied independently to
/// every nutrient and soil-zone pool.
pub fn uptakeCompetition(
    minimum_population_fraction: f64,
    root_biome_fraction: f64,
    significance_threshold_g: f64,
    population_demand_g_per_h: UptakeCompetition,
    combined_demand_g_per_h: UptakeCompetition,
) !UptakeCompetition {
    inline for (.{ minimum_population_fraction, root_biome_fraction, significance_threshold_g }) |value| {
        if (!std.math.isFinite(value)) return error.NonFiniteUptakeCompetitionInput;
    }
    if (minimum_population_fraction < 0 or root_biome_fraction < 0 or significance_threshold_g < 0) return error.InvalidUptakeCompetitionInput;

    var result: UptakeCompetition = undefined;
    inline for (@typeInfo(UptakeCompetition).@"struct".fields) |field| {
        const population_demand = @field(population_demand_g_per_h, field.name);
        const combined_demand = @field(combined_demand_g_per_h, field.name);
        if (!std.math.isFinite(population_demand) or population_demand < 0 or !std.math.isFinite(combined_demand) or combined_demand < 0) return error.InvalidUptakeCompetitionDemand;
        @field(result, field.name) = if (combined_demand > significance_threshold_g)
            @max(minimum_population_fraction, population_demand / combined_demand)
        else
            root_biome_fraction;
        if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteUptakeCompetitionResult;
    }
    return result;
}

/// Heap-owned STARTQ root/mycorrhizal state. Every extent is derived from the
/// runscript and mapped soil profiles; the historical NR=10 ceiling is absent.
pub const State = struct {
    allocator: std.mem.Allocator,
    plant_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    planting_layer_by_plant: []usize,
    active_root_axis_count: []usize,
    roots_dead: []bool,
    seed_volume_m3_per_plant: []f64,
    seed_length_m_per_plant: []f64,
    seed_surface_area_m2_per_plant: []f64,
    ammonium_uptake_g_n_per_h: []f64,
    nitrate_uptake_g_n_per_h: []f64,
    phosphate_uptake_g_p_per_h: []f64,
    fixation_uptake_g_n_per_h: []f64,
    /// Layer-resolved RUPNF owner; plant totals remain above for outputs.
    fixation_uptake_g_n_per_h_by_layer: []f64,
    current_porosity_fraction: []f64,
    initial_porosity_fraction: []f64,
    water_uptake_m3_per_h: []f64,
    total_water_potential_mpa: []f64,
    osmotic_water_potential_mpa: []f64,
    turgor_water_potential_mpa: []f64,
    mobile_carbon_g: []f64,
    mobile_nitrogen_g: []f64,
    mobile_phosphorus_g: []f64,
    symbiont_structural_carbon_g_c: []f64,
    symbiont_structural_nitrogen_g_n: []f64,
    symbiont_structural_phosphorus_g_p: []f64,
    symbiont_mobile_carbon_g_c: []f64,
    symbiont_mobile_nitrogen_g_n: []f64,
    symbiont_mobile_phosphorus_g_p: []f64,
    symbiotic_respiration_actual_g_c_per_h: []f64,
    symbiotic_respiration_oxygen_unlimited_g_c_per_h: []f64,
    mobile_carbon_concentration_g_per_g: []f64,
    mobile_nitrogen_concentration_g_per_g: []f64,
    mobile_phosphorus_concentration_g_per_g: []f64,
    salt_concentration_mol_per_g_c: []f64,
    maximum_protein_carbon_g_per_g_c: []f64,
    total_carbon_g: []f64,
    primary_root_carbon_g: []f64,
    protein_carbon_g: []f64,
    primary_radius_m: []f64,
    secondary_radius_m: []f64,
    reference_primary_radius_m: []f64,
    reference_secondary_radius_m: []f64,
    projected_area_m2: []f64,
    active_length_m: []f64,
    water_fraction: []f64,
    aqueous_volume_m3: []f64,
    gaseous_volume_m3: []f64,
    root_length_m_per_plant: []f64,
    root_length_density_m_per_m3: []f64,
    root_surface_area_m2_per_plant: []f64,
    average_secondary_length_m: []f64,
    sink_strength_m: []f64,
    ammonium_uptake_nonband_g_n_per_h: []f64,
    nitrate_uptake_nonband_g_n_per_h: []f64,
    phosphate_h2_uptake_nonband_g_p_per_h: []f64,
    phosphate_h_uptake_nonband_g_p_per_h: []f64,
    ammonium_uptake_band_g_n_per_h: []f64,
    nitrate_uptake_band_g_n_per_h: []f64,
    phosphate_h2_uptake_band_g_p_per_h: []f64,
    phosphate_h_uptake_band_g_p_per_h: []f64,
    previous_ammonium_uptake_nonband_g_n_per_h: []f64,
    previous_nitrate_uptake_nonband_g_n_per_h: []f64,
    previous_phosphate_h2_uptake_nonband_g_p_per_h: []f64,
    previous_phosphate_h_uptake_nonband_g_p_per_h: []f64,
    previous_ammonium_uptake_band_g_n_per_h: []f64,
    previous_nitrate_uptake_band_g_n_per_h: []f64,
    previous_phosphate_h2_uptake_band_g_p_per_h: []f64,
    previous_phosphate_h_uptake_band_g_p_per_h: []f64,
    ammonium_demand_nonband_g_n_per_h: []f64,
    nitrate_demand_nonband_g_n_per_h: []f64,
    phosphate_h2_demand_nonband_g_p_per_h: []f64,
    phosphate_h_demand_nonband_g_p_per_h: []f64,
    ammonium_demand_band_g_n_per_h: []f64,
    nitrate_demand_band_g_n_per_h: []f64,
    phosphate_h2_demand_band_g_p_per_h: []f64,
    phosphate_h_demand_band_g_p_per_h: []f64,
    previous_ammonium_demand_nonband_g_n_per_h: []f64,
    previous_nitrate_demand_nonband_g_n_per_h: []f64,
    previous_phosphate_h2_demand_nonband_g_p_per_h: []f64,
    previous_phosphate_h_demand_nonband_g_p_per_h: []f64,
    previous_ammonium_demand_band_g_n_per_h: []f64,
    previous_nitrate_demand_band_g_n_per_h: []f64,
    previous_phosphate_h2_demand_band_g_p_per_h: []f64,
    previous_phosphate_h_demand_band_g_p_per_h: []f64,
    oxygen_uptake_g_o_per_h: []f64,
    oxygen_uptake_from_soil_g_o_per_h: []f64,
    oxygen_uptake_from_root_pool_g_o_per_h: []f64,
    oxygen_demand_g_o_per_h: []f64,
    respiration_unlimited_by_oxygen_g_c_per_h: []f64,
    respiration_unlimited_by_carbon_g_c_per_h: []f64,
    actual_respiration_g_c_per_h: []f64,
    oxygen_process_constraint_fraction: []f64,
    ammonium_assimilation_g_n_per_h: []f64,
    band_ammonium_assimilation_g_n_per_h: []f64,
    nitrate_assimilation_g_n_per_h: []f64,
    band_nitrate_assimilation_g_n_per_h: []f64,
    phosphate_h2_assimilation_g_p_per_h: []f64,
    phosphate_h_assimilation_g_p_per_h: []f64,
    band_phosphate_h2_assimilation_g_p_per_h: []f64,
    band_phosphate_h_assimilation_g_p_per_h: []f64,
    gaseous_carbon_dioxide_g_c: []f64,
    aqueous_carbon_dioxide_g_c: []f64,
    gaseous_oxygen_g_o: []f64,
    aqueous_oxygen_g_o: []f64,
    gaseous_methane_g_c: []f64,
    aqueous_methane_g_c: []f64,
    gaseous_nitrous_oxide_g_n: []f64,
    aqueous_nitrous_oxide_g_n: []f64,
    gaseous_ammonia_g_n: []f64,
    aqueous_ammonia_g_n: []f64,
    gaseous_hydrogen_g_h: []f64,
    aqueous_hydrogen_g_h: []f64,
    /// Accepted extensive transactions, indexed by rootIndex * 6 + gas.
    /// Gas order is CO2-C, CH4-C, N2O-N, NH3-N, H2-H, and O2-O.
    soil_to_root_gas_exchange_g_per_h: []f64,
    aqueous_to_gaseous_root_exchange_g_per_h: []f64,
    atmosphere_to_root_gas_exchange_g_per_h: []f64,
    withdrawal_carbon_dioxide_loss_g_c_per_h: []f64,
    withdrawal_oxygen_loss_g_o_per_h: []f64,
    withdrawal_methane_loss_g_c_per_h: []f64,
    withdrawal_nitrous_oxide_loss_g_n_per_h: []f64,
    withdrawal_ammonia_loss_g_n_per_h: []f64,
    withdrawal_hydrogen_loss_g_h_per_h: []f64,
    combustion_carbon_loss_g_c_per_h: []f64,
    combustion_nitrogen_loss_g_n_per_h: []f64,
    combustion_phosphorus_loss_g_p_per_h: []f64,
    symbiont_combustion_g_c_per_h: []f64,
    root_combustion_g_c_per_h: []f64,
    combustion_salt_loss_mol_per_h: []f64,
    salt_content_mol: []f64,
    salt_uptake_mol_per_h: []f64,
    exudate_carbon_exchange_g_c_per_h: []f64,
    exudate_nitrogen_exchange_g_n_per_h: []f64,
    exudate_phosphorus_exchange_g_p_per_h: []f64,
    carbon_dioxide_advection_g_c_per_h: []f64,
    carbon_dioxide_diffusion_g_c_per_h: []f64,
    carbon_dioxide_solubilization_g_c_per_h: []f64,
    aqueous_carbon_dioxide_reaction_g_c_per_h: []f64,
    axis_depth_m: []f64,
    axis_primary_length_m: []f64,
    axis_primary_count: []f64,
    axis_primary_carbon_g: []f64,
    axis_primary_nitrogen_g: []f64,
    axis_primary_phosphorus_g: []f64,
    axis_secondary_length_m: []f64,
    axis_secondary_count: []f64,
    axis_secondary_carbon_g: []f64,
    axis_secondary_nitrogen_g: []f64,
    axis_secondary_phosphorus_g: []f64,

    pub fn init(allocator: std.mem.Allocator, plant_count: usize, soil_layer_count: usize, root_axis_count: usize) !State {
        if (plant_count == 0 or soil_layer_count == 0 or root_axis_count == 0) return error.InvalidPlantRootDimensions;
        const domain_layer_count = try std.math.mul(usize, try std.math.mul(usize, plant_count, biological_domain_count), soil_layer_count);
        const domain_axis_count = try std.math.mul(usize, try std.math.mul(usize, plant_count, biological_domain_count), root_axis_count);
        const domain_layer_axis_count = try std.math.mul(usize, domain_layer_count, root_axis_count);
        const salt_count = try std.math.mul(usize, domain_layer_count, salt_species_count);
        const exudate_count = try std.math.mul(usize, domain_layer_count, organic_substrate_count);
        var result: State = undefined;
        result.allocator = allocator;
        result.plant_count = plant_count;
        result.soil_layer_count = soil_layer_count;
        result.root_axis_count = root_axis_count;
        result.planting_layer_by_plant = try allocator.alloc(usize, plant_count);
        errdefer allocator.free(result.planting_layer_by_plant);
        @memset(result.planting_layer_by_plant, 0);
        result.active_root_axis_count = try allocator.alloc(usize, plant_count);
        errdefer allocator.free(result.active_root_axis_count);
        @memset(result.active_root_axis_count, 0);
        result.roots_dead = try allocator.alloc(bool, plant_count);
        errdefer allocator.free(result.roots_dead);
        @memset(result.roots_dead, true);
        var allocated: usize = 0;
        errdefer freeAllocated(&result, allocated);
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            const count = fieldCount(field.name, plant_count, domain_layer_count, domain_axis_count, domain_layer_axis_count, salt_count, exudate_count);
            @field(result, field.name) = try allocator.alloc(f64, count);
            @memset(@field(result, field.name), 0);
            allocated += 1;
        };
        return result;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.allocator.free(self.roots_dead);
        self.allocator.free(self.active_root_axis_count);
        self.allocator.free(self.planting_layer_by_plant);
        self.* = undefined;
    }

    pub fn layerIndex(self: State, plant: usize, domain: usize, layer: usize) !usize {
        if (plant >= self.plant_count or domain >= biological_domain_count or layer >= self.soil_layer_count) return error.PlantRootIndexOutOfBounds;
        return (plant * biological_domain_count + domain) * self.soil_layer_count + layer;
    }

    pub fn axisIndex(self: State, plant: usize, domain: usize, axis: usize) !usize {
        if (plant >= self.plant_count or domain >= biological_domain_count or axis >= self.root_axis_count) return error.PlantRootIndexOutOfBounds;
        return (plant * biological_domain_count + domain) * self.root_axis_count + axis;
    }

    pub fn layerAxisIndex(self: State, plant: usize, domain: usize, layer: usize, axis: usize) !usize {
        return (try self.layerIndex(plant, domain, layer)) * self.root_axis_count + axis;
    }

    pub fn saltIndex(self: State, plant: usize, domain: usize, layer: usize, species: SaltSpecies) !usize {
        return try std.math.add(usize, try std.math.mul(usize, try self.layerIndex(plant, domain, layer), salt_species_count), @intFromEnum(species));
    }

    pub fn substrateIndex(self: State, plant: usize, domain: usize, layer: usize, substrate: usize) !usize {
        if (substrate >= organic_substrate_count) return error.PlantRootSubstrateIndexOutOfBounds;
        return try std.math.add(usize, try std.math.mul(usize, try self.layerIndex(plant, domain, layer), organic_substrate_count), substrate);
    }

    pub fn initializePlant(self: *State, plant: usize, traits: PlantTraits, planting_layer: usize, seeding_depth_m: f64, parameters: InitializationParameters) !void {
        try parameters.validate();
        if (plant >= self.plant_count or planting_layer >= self.soil_layer_count) return error.PlantRootIndexOutOfBounds;
        if (!std.math.isFinite(seeding_depth_m) or seeding_depth_m < 0) return error.InvalidPlantingDepth;
        self.planting_layer_by_plant[plant] = planting_layer;
        self.active_root_axis_count[plant] = 0;
        self.roots_dead[plant] = true;
        self.current_porosity_fraction[plant] = traits.roots.root_porosity_fraction;
        self.initial_porosity_fraction[plant] = traits.roots.root_porosity_fraction;
        const maximum_protein = @min(
            traits.organ_nitrogen_to_carbon_ratio.root * parameters.root_nitrogen_to_maximum_protein_multiplier,
            traits.organ_phosphorus_to_carbon_ratio.root * parameters.root_phosphorus_to_maximum_protein_multiplier,
        );
        if (!std.math.isFinite(maximum_protein) or maximum_protein < 0) return error.InvalidRootProteinConcentration;
        for (0..biological_domain_count) |domain| {
            const primary_radius_m = if (domain == 0) traits.roots.primary_root_radius_m else parameters.mycorrhizal_radius_m;
            const secondary_radius_m = if (domain == 0) traits.roots.secondary_root_radius_m else parameters.mycorrhizal_radius_m;
            if (!std.math.isFinite(primary_radius_m) or primary_radius_m <= 0 or !std.math.isFinite(secondary_radius_m) or secondary_radius_m <= 0) return error.InvalidRootRadius;
            for (0..self.soil_layer_count) |layer| {
                const index = try self.layerIndex(plant, domain, layer);
                self.total_water_potential_mpa[index] = parameters.initial_total_water_potential_mpa;
                self.osmotic_water_potential_mpa[index] = traits.water_relations.osmotic_potential_mpa + parameters.osmotic_water_potential_decrement_mpa;
                self.turgor_water_potential_mpa[index] = @max(0.0, parameters.initial_total_water_potential_mpa - self.osmotic_water_potential_mpa[index]);
                self.maximum_protein_carbon_g_per_g_c[index] = maximum_protein;
                self.primary_radius_m[index] = primary_radius_m;
                self.secondary_radius_m[index] = secondary_radius_m;
                self.reference_primary_radius_m[index] = primary_radius_m;
                self.reference_secondary_radius_m[index] = secondary_radius_m;
                self.active_length_m[index] = parameters.initial_active_length_m;
                self.water_fraction[index] = parameters.initial_water_fraction;
            }
            for (0..self.root_axis_count) |axis| self.axis_depth_m[try self.axisIndex(plant, domain, axis)] = seeding_depth_m;
        }
    }

    /// Stores STARTQ SDVL/SDLG/SDAR seed geometry. GROSUB applies it only to
    /// the root domain in the planting layer.
    pub fn setSeedGeometry(self: *State, plant: usize, volume_m3_per_plant: f64, length_m_per_plant: f64, surface_area_m2_per_plant: f64) !void {
        if (plant >= self.plant_count) return error.PlantRootIndexOutOfBounds;
        inline for (.{ volume_m3_per_plant, length_m_per_plant, surface_area_m2_per_plant }) |value| {
            if (!std.math.isFinite(value) or value < 0) return error.InvalidSeedGeometry;
        }
        self.seed_volume_m3_per_plant[plant] = volume_m3_per_plant;
        self.seed_length_m_per_plant[plant] = length_m_per_plant;
        self.seed_surface_area_m2_per_plant[plant] = surface_area_m2_per_plant;
    }

    /// Reconstructs one runtime plant after a terminating harvest. Every
    /// plant-major persistent pool and diagnostic is cleared before STARTQ
    /// geometry and water status are applied, so a new crop cannot inherit
    /// C/N/P, gas, salt, exudate, demand, or root-axis history.
    pub fn reconstructPlant(self: *State, plant: usize, traits: PlantTraits, planting_layer: usize, seeding_depth_m: f64, parameters: InitializationParameters) !void {
        if (plant >= self.plant_count) return error.PlantRootIndexOutOfBounds;
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            const values = @field(self, field.name);
            if (values.len % self.plant_count != 0) return error.InvalidPlantRootFieldExtent;
            const per_plant = values.len / self.plant_count;
            @memset(values[plant * per_plant .. (plant + 1) * per_plant], 0);
        };
        self.planting_layer_by_plant[plant] = 0;
        self.active_root_axis_count[plant] = 0;
        self.roots_dead[plant] = true;
        try self.initializePlant(plant, traits, planting_layer, seeding_depth_m, parameters);
    }

    pub const LayerTopology = struct {
        primary_axis_count: f64,
        secondary_axis_count: f64,
        primary_length_m_per_plant: f64,
        secondary_length_m_per_cell: f64,
        total_root_length_m: f64,
        root_length_m_per_plant: f64,
        root_length_density_m_per_m3: f64,
        average_secondary_length_m: f64,
    };

    pub const SourceOrderLayerCarbonTotals = struct {
        active_secondary_carbon_g_c: f64,
        actual_primary_and_secondary_carbon_g_c: f64,
    };

    /// Exact GROSUB 8210-8218 WTRTL/WTRTD operands for one root domain/layer.
    pub fn sourceOrderLayerCarbonTotals(
        primary_carbon_g_c_by_axis: []const f64,
        secondary_carbon_g_c_by_axis: []const f64,
    ) !SourceOrderLayerCarbonTotals {
        if (primary_carbon_g_c_by_axis.len != secondary_carbon_g_c_by_axis.len)
            return error.RootCarbonTotalDimensionMismatch;
        var result: SourceOrderLayerCarbonTotals = .{
            .active_secondary_carbon_g_c = 0,
            .actual_primary_and_secondary_carbon_g_c = 0,
        };
        for (primary_carbon_g_c_by_axis, secondary_carbon_g_c_by_axis) |primary, secondary| {
            if (!std.math.isFinite(primary) or !std.math.isFinite(secondary) or primary < 0 or secondary < 0)
                return error.InvalidRootCarbonTotalInput;
            result.active_secondary_carbon_g_c += secondary;
            result.actual_primary_and_secondary_carbon_g_c += secondary + primary;
        }
        if (!std.math.isFinite(result.active_secondary_carbon_g_c) or
            !std.math.isFinite(result.actual_primary_and_secondary_carbon_g_c))
            return error.NonFiniteRootCarbonTotal;
        return result;
    }

    /// GROSUB WTRTL active-root C reconstruction. Secondary-root C remains in
    /// its physical layer; the complete primary C of each axis is assigned to
    /// that axis's deepest occupied tip layer (NINR).
    pub fn refreshActiveCarbonByLayer(self: *State, plant: usize, domain: usize) !void {
        if (plant >= self.plant_count or domain >= biological_domain_count) return error.PlantRootIndexOutOfBounds;
        for (0..self.soil_layer_count) |layer| {
            const root = try self.layerIndex(plant, domain, layer);
            self.total_carbon_g[root] = 0;
        }
        const active_axes = self.active_root_axis_count[plant];
        if (active_axes > self.root_axis_count) return error.InvalidActiveRootAxisCount;
        for (0..active_axes) |axis| {
            var total_primary_carbon_g_c: f64 = 0;
            var tip_layer: usize = self.planting_layer_by_plant[plant];
            if (tip_layer >= self.soil_layer_count) return error.PlantRootLayerOutOfBounds;
            for (0..self.soil_layer_count) |layer| {
                const axis_layer = try self.layerAxisIndex(plant, domain, layer, axis);
                const primary = self.axis_primary_carbon_g[axis_layer];
                const secondary = self.axis_secondary_carbon_g[axis_layer];
                const primary_length = self.axis_primary_length_m[axis_layer];
                inline for (.{ primary, secondary, primary_length }) |value|
                    if (!std.math.isFinite(value) or value < 0) return error.InvalidRootActiveCarbonState;
                total_primary_carbon_g_c += primary;
                self.total_carbon_g[try self.layerIndex(plant, domain, layer)] += secondary;
                if (primary > 0 or primary_length > 0) tip_layer = layer;
            }
            const tip = try self.layerIndex(plant, domain, tip_layer);
            self.total_carbon_g[tip] += total_primary_carbon_g_c;
            if (!std.math.isFinite(self.total_carbon_g[tip])) return error.NonFiniteRootActiveCarbon;
        }
        for (0..self.soil_layer_count) |layer| {
            const value = self.total_carbon_g[try self.layerIndex(plant, domain, layer)];
            if (!std.math.isFinite(value) or value < 0) return error.NonFiniteRootActiveCarbon;
        }
    }

    /// GROSUB RTLGT/RTLGP/RTDNP/RTLGA aggregation over runtime root axes.
    pub fn layerTopology(self: State, plant: usize, domain: usize, layer: usize, plant_population_count: f64, layer_thickness_m: f64, woody_root_fraction: f64, minimum_average_secondary_length_m: f64) !LayerTopology {
        inline for (.{ plant_population_count, layer_thickness_m, woody_root_fraction, minimum_average_secondary_length_m }) |value| if (!std.math.isFinite(value)) return error.NonFiniteRootTopologyInput;
        if (plant_population_count <= 0 or layer_thickness_m <= 0 or woody_root_fraction < 0 or woody_root_fraction > 1 or minimum_average_secondary_length_m < 0) return error.InvalidRootTopologyInput;
        var primary_count: f64 = 0;
        var secondary_count: f64 = 0;
        var primary_length: f64 = 0;
        var secondary_length: f64 = 0;
        for (0..self.root_axis_count) |axis| {
            const index = try self.layerAxisIndex(plant, domain, layer, axis);
            inline for (.{ self.axis_primary_count[index], self.axis_secondary_count[index], self.axis_primary_length_m[index], self.axis_secondary_length_m[index] }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootTopologyState;
            primary_count += self.axis_primary_count[index];
            secondary_count += self.axis_secondary_count[index];
            primary_length += self.axis_primary_length_m[index];
            secondary_length += self.axis_secondary_length_m[index];
        }
        const total_length = primary_length * plant_population_count + secondary_length;
        const seed_length_m = if (domain == 0 and layer == self.planting_layer_by_plant[plant]) self.seed_length_m_per_plant[plant] else 0;
        if (!std.math.isFinite(seed_length_m) or seed_length_m < 0) return error.InvalidSeedGeometry;
        const length_per_plant = total_length / plant_population_count * woody_root_fraction + seed_length_m;
        return .{
            .primary_axis_count = primary_count,
            .secondary_axis_count = secondary_count,
            .primary_length_m_per_plant = primary_length,
            .secondary_length_m_per_cell = secondary_length,
            .total_root_length_m = total_length,
            .root_length_m_per_plant = length_per_plant,
            .root_length_density_m_per_m3 = length_per_plant / layer_thickness_m,
            .average_secondary_length_m = if (secondary_count > 0) @max(minimum_average_secondary_length_m, secondary_length / secondary_count) else minimum_average_secondary_length_m,
        };
    }

    pub fn refreshLayerMorphology(self: *State, plant: usize, domain: usize, layer: usize, plant_population_count: f64, layer_thickness_m: f64, woody_root_fraction: f64, root_porosity_fraction: f64, volume_m3_per_g_c: f64, root_geometry_pi: f64, vascular_growth_habit: bool, parameters: MorphologyParameters) !LayerTopology {
        try parameters.validate();
        if (!std.math.isFinite(root_porosity_fraction) or root_porosity_fraction < 0 or root_porosity_fraction >= 1 or !std.math.isFinite(volume_m3_per_g_c) or volume_m3_per_g_c <= 0 or !std.math.isFinite(root_geometry_pi) or root_geometry_pi <= 0) return error.InvalidRootMorphologyInput;
        const topology = try self.layerTopology(plant, domain, layer, plant_population_count, layer_thickness_m, woody_root_fraction, parameters.minimum_average_secondary_length_m);
        const index = try self.layerIndex(plant, domain, layer);
        const total_carbon = self.total_carbon_g[index];
        const primary_carbon = self.primary_root_carbon_g[index];
        const turgor = self.turgor_water_potential_mpa[index];
        if (!std.math.isFinite(total_carbon) or total_carbon < 0 or !std.math.isFinite(primary_carbon) or primary_carbon < 0 or primary_carbon > total_carbon or !std.math.isFinite(turgor) or turgor < 0) return error.InvalidRootMorphologyState;
        const reference_primary_radius_m = self.reference_primary_radius_m[index];
        const reference_secondary_radius_m = self.reference_secondary_radius_m[index];
        if (!std.math.isFinite(reference_primary_radius_m) or reference_primary_radius_m <= 0 or !std.math.isFinite(reference_secondary_radius_m) or reference_secondary_radius_m <= 0) return error.InvalidRootMorphologyState;
        const has_live_roots = topology.total_root_length_m > 0 and total_carbon > 0;
        var root_volume_m3: f64 = 0;
        var surface_area_m2: f64 = 0;
        if (has_live_roots) {
            self.primary_radius_m[index] = @max(reference_primary_radius_m, (1.0 + turgor / parameters.root_elastic_modulus_mpa) * reference_primary_radius_m);
            self.secondary_radius_m[index] = @max(reference_secondary_radius_m, (1.0 + turgor / parameters.root_elastic_modulus_mpa) * reference_secondary_radius_m);
            const secondary_carbon_g = total_carbon - primary_carbon;
            root_volume_m3 = @max(
                root_geometry_pi * reference_secondary_radius_m * reference_secondary_radius_m * topology.secondary_length_m_per_cell,
                secondary_carbon_g * volume_m3_per_g_c * turgor,
            );
            const circumference_factor = 2.0 * root_geometry_pi;
            surface_area_m2 = circumference_factor * self.primary_radius_m[index] * topology.primary_length_m_per_plant * plant_population_count + circumference_factor * self.secondary_radius_m[index] * topology.secondary_length_m_per_cell;
            if (vascular_growth_habit) surface_area_m2 *= parameters.minimum_average_secondary_length_m / topology.average_secondary_length_m;
        } else {
            self.primary_radius_m[index] = reference_primary_radius_m;
            self.secondary_radius_m[index] = reference_secondary_radius_m;
            inline for (.{
                .{ "gaseous_carbon_dioxide_g_c", "aqueous_carbon_dioxide_g_c", "withdrawal_carbon_dioxide_loss_g_c_per_h" },
                .{ "gaseous_oxygen_g_o", "aqueous_oxygen_g_o", "withdrawal_oxygen_loss_g_o_per_h" },
                .{ "gaseous_methane_g_c", "aqueous_methane_g_c", "withdrawal_methane_loss_g_c_per_h" },
                .{ "gaseous_nitrous_oxide_g_n", "aqueous_nitrous_oxide_g_n", "withdrawal_nitrous_oxide_loss_g_n_per_h" },
                .{ "gaseous_ammonia_g_n", "aqueous_ammonia_g_n", "withdrawal_ammonia_loss_g_n_per_h" },
                .{ "gaseous_hydrogen_g_h", "aqueous_hydrogen_g_h", "withdrawal_hydrogen_loss_g_h_per_h" },
            }) |fields| {
                const released = @field(self, fields[0])[index] + @field(self, fields[1])[index];
                if (!std.math.isFinite(released)) return error.InvalidRootMorphologyState;
                @field(self, fields[2])[plant] -= released;
                @field(self, fields[0])[index] = 0;
                @field(self, fields[1])[index] = 0;
            }
        }
        const seed_volume_m3 = if (domain == 0 and layer == self.planting_layer_by_plant[plant]) self.seed_volume_m3_per_plant[plant] * plant_population_count else 0;
        if (!std.math.isFinite(seed_volume_m3) or seed_volume_m3 < 0) return error.InvalidSeedGeometry;
        const total_volume_m3 = root_volume_m3 + seed_volume_m3;
        self.gaseous_volume_m3[index] = root_porosity_fraction * total_volume_m3;
        self.aqueous_volume_m3[index] = (1.0 - root_porosity_fraction) * total_volume_m3;
        self.root_length_m_per_plant[index] = topology.root_length_m_per_plant;
        self.root_length_density_m_per_m3[index] = topology.root_length_density_m_per_m3;
        self.average_secondary_length_m[index] = topology.average_secondary_length_m;
        const seed_surface_area_m2 = if (domain == 0 and layer == self.planting_layer_by_plant[plant]) self.seed_surface_area_m2_per_plant[plant] else 0;
        if (!std.math.isFinite(seed_surface_area_m2) or seed_surface_area_m2 < 0) return error.InvalidSeedGeometry;
        self.root_surface_area_m2_per_plant[index] = surface_area_m2 / plant_population_count * woody_root_fraction + seed_surface_area_m2;
        return topology;
    }

    pub fn validateFinite(self: State) !void {
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) for (@field(self, field.name), 0..) |value, index| if (!std.math.isFinite(value)) {
            std.log.err("non-finite plant root state: field={s} index={d} value={e}", .{ field.name, index, value });
            return error.NonFinitePlantRootState;
        };
    }

    /// UPTAKE hourly flux initialization. Persistent C/N/P, gas inventories,
    /// morphology, potentials, and axis state are deliberately not cleared.
    pub fn resetHourlyFluxes(self: *State) void {
        @memcpy(self.previous_ammonium_uptake_nonband_g_n_per_h, self.ammonium_uptake_nonband_g_n_per_h);
        @memcpy(self.previous_nitrate_uptake_nonband_g_n_per_h, self.nitrate_uptake_nonband_g_n_per_h);
        @memcpy(self.previous_phosphate_h2_uptake_nonband_g_p_per_h, self.phosphate_h2_uptake_nonband_g_p_per_h);
        @memcpy(self.previous_phosphate_h_uptake_nonband_g_p_per_h, self.phosphate_h_uptake_nonband_g_p_per_h);
        @memcpy(self.previous_ammonium_uptake_band_g_n_per_h, self.ammonium_uptake_band_g_n_per_h);
        @memcpy(self.previous_nitrate_uptake_band_g_n_per_h, self.nitrate_uptake_band_g_n_per_h);
        @memcpy(self.previous_phosphate_h2_uptake_band_g_p_per_h, self.phosphate_h2_uptake_band_g_p_per_h);
        @memcpy(self.previous_phosphate_h_uptake_band_g_p_per_h, self.phosphate_h_uptake_band_g_p_per_h);
        @memcpy(self.previous_ammonium_demand_nonband_g_n_per_h, self.ammonium_demand_nonband_g_n_per_h);
        @memcpy(self.previous_nitrate_demand_nonband_g_n_per_h, self.nitrate_demand_nonband_g_n_per_h);
        @memcpy(self.previous_phosphate_h2_demand_nonband_g_p_per_h, self.phosphate_h2_demand_nonband_g_p_per_h);
        @memcpy(self.previous_phosphate_h_demand_nonband_g_p_per_h, self.phosphate_h_demand_nonband_g_p_per_h);
        @memcpy(self.previous_ammonium_demand_band_g_n_per_h, self.ammonium_demand_band_g_n_per_h);
        @memcpy(self.previous_nitrate_demand_band_g_n_per_h, self.nitrate_demand_band_g_n_per_h);
        @memcpy(self.previous_phosphate_h2_demand_band_g_p_per_h, self.phosphate_h2_demand_band_g_p_per_h);
        @memcpy(self.previous_phosphate_h_demand_band_g_p_per_h, self.phosphate_h_demand_band_g_p_per_h);
        inline for (.{
            self.ammonium_uptake_g_n_per_h,
            self.nitrate_uptake_g_n_per_h,
            self.phosphate_uptake_g_p_per_h,
            self.fixation_uptake_g_n_per_h,
            self.fixation_uptake_g_n_per_h_by_layer,
            self.water_uptake_m3_per_h,
            self.ammonium_uptake_nonband_g_n_per_h,
            self.nitrate_uptake_nonband_g_n_per_h,
            self.phosphate_h2_uptake_nonband_g_p_per_h,
            self.phosphate_h_uptake_nonband_g_p_per_h,
            self.ammonium_uptake_band_g_n_per_h,
            self.nitrate_uptake_band_g_n_per_h,
            self.phosphate_h2_uptake_band_g_p_per_h,
            self.phosphate_h_uptake_band_g_p_per_h,
            self.ammonium_demand_nonband_g_n_per_h,
            self.nitrate_demand_nonband_g_n_per_h,
            self.phosphate_h2_demand_nonband_g_p_per_h,
            self.phosphate_h_demand_nonband_g_p_per_h,
            self.ammonium_demand_band_g_n_per_h,
            self.nitrate_demand_band_g_n_per_h,
            self.phosphate_h2_demand_band_g_p_per_h,
            self.phosphate_h_demand_band_g_p_per_h,
            self.oxygen_uptake_g_o_per_h,
            self.oxygen_uptake_from_soil_g_o_per_h,
            self.oxygen_uptake_from_root_pool_g_o_per_h,
            self.oxygen_demand_g_o_per_h,
            self.respiration_unlimited_by_oxygen_g_c_per_h,
            self.respiration_unlimited_by_carbon_g_c_per_h,
            self.actual_respiration_g_c_per_h,
            self.symbiotic_respiration_actual_g_c_per_h,
            self.symbiotic_respiration_oxygen_unlimited_g_c_per_h,
            self.ammonium_assimilation_g_n_per_h,
            self.band_ammonium_assimilation_g_n_per_h,
            self.nitrate_assimilation_g_n_per_h,
            self.band_nitrate_assimilation_g_n_per_h,
            self.phosphate_h2_assimilation_g_p_per_h,
            self.phosphate_h_assimilation_g_p_per_h,
            self.band_phosphate_h2_assimilation_g_p_per_h,
            self.band_phosphate_h_assimilation_g_p_per_h,
            self.carbon_dioxide_advection_g_c_per_h,
            self.carbon_dioxide_diffusion_g_c_per_h,
            self.carbon_dioxide_solubilization_g_c_per_h,
            self.aqueous_carbon_dioxide_reaction_g_c_per_h,
            self.soil_to_root_gas_exchange_g_per_h,
            self.aqueous_to_gaseous_root_exchange_g_per_h,
            self.atmosphere_to_root_gas_exchange_g_per_h,
            self.withdrawal_carbon_dioxide_loss_g_c_per_h,
            self.withdrawal_oxygen_loss_g_o_per_h,
            self.withdrawal_methane_loss_g_c_per_h,
            self.withdrawal_nitrous_oxide_loss_g_n_per_h,
            self.withdrawal_ammonia_loss_g_n_per_h,
            self.withdrawal_hydrogen_loss_g_h_per_h,
            self.combustion_carbon_loss_g_c_per_h,
            self.combustion_nitrogen_loss_g_n_per_h,
            self.combustion_phosphorus_loss_g_p_per_h,
            self.symbiont_combustion_g_c_per_h,
            self.root_combustion_g_c_per_h,
            self.sink_strength_m,
            self.combustion_salt_loss_mol_per_h,
            self.salt_uptake_mol_per_h,
            self.exudate_carbon_exchange_g_c_per_h,
            self.exudate_nitrogen_exchange_g_n_per_h,
            self.exudate_phosphorus_exchange_g_p_per_h,
        }) |values| @memset(values, 0);
    }
};

pub const NegativeStructuralCarbonCleanup = struct {
    structural_carbon_g_c: f64,
    mobile_carbon_g_c: f64,
    removed_deficit_g_c: f64,
};

/// GROSUB lines 7301--7307 source-order characterization. The legacy model
/// zeroes a negative structural pool and charges the deficit to shared mobile
/// carbon. Production state validation intentionally rejects such a state
/// before it reaches this repair path.
pub fn sourceOrderNegativeStructuralCarbonCleanup(
    structural_carbon_g_c: f64,
    mobile_carbon_g_c: f64,
) !NegativeStructuralCarbonCleanup {
    if (!std.math.isFinite(structural_carbon_g_c) or
        !std.math.isFinite(mobile_carbon_g_c))
        return error.NonFiniteNegativeStructuralCarbonCleanup;
    if (structural_carbon_g_c >= 0) return .{
        .structural_carbon_g_c = structural_carbon_g_c,
        .mobile_carbon_g_c = mobile_carbon_g_c,
        .removed_deficit_g_c = 0,
    };
    return .{
        .structural_carbon_g_c = 0,
        .mobile_carbon_g_c = mobile_carbon_g_c + structural_carbon_g_c,
        .removed_deficit_g_c = -structural_carbon_g_c,
    };
}

/// GROSUB lines 7428--7429 gate for publishing live root morphology.
pub fn sourceOrderRootMorphologyIsActive(
    total_root_length_m: f64,
    total_root_carbon_g_c: f64,
    plant_population_count: f64,
    presence_threshold: f64,
) !bool {
    inline for (.{
        total_root_length_m,
        total_root_carbon_g_c,
        plant_population_count,
        presence_threshold,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidRootMorphologyGateInput;
    return total_root_length_m > presence_threshold and
        total_root_carbon_g_c > presence_threshold and
        plant_population_count > presence_threshold;
}

/// GROSUB lines 7431--7435 and 7512--7517. Numerically negligible layers
/// publish zero length density instead of dividing by their thickness.
pub fn sourceOrderRootLengthDensity(
    root_length_m_per_plant: f64,
    layer_thickness_m: f64,
    minimum_layer_thickness_m: f64,
) !f64 {
    inline for (.{ root_length_m_per_plant, layer_thickness_m, minimum_layer_thickness_m }) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootLengthDensityInput;
    return if (layer_thickness_m > minimum_layer_thickness_m)
        root_length_m_per_plant / layer_thickness_m
    else
        0;
}

fn fieldCount(comptime name: []const u8, plant_count: usize, domain_layer_count: usize, domain_axis_count: usize, domain_layer_axis_count: usize, salt_count: usize, exudate_count: usize) usize {
    if (comptime std.mem.startsWith(u8, name, "axis_depth_")) return domain_axis_count;
    if (comptime std.mem.startsWith(u8, name, "axis_")) return domain_layer_axis_count;
    if (comptime std.mem.eql(u8, name, "salt_content_mol") or std.mem.eql(u8, name, "salt_uptake_mol_per_h") or std.mem.eql(u8, name, "combustion_salt_loss_mol_per_h")) return salt_count;
    if (comptime std.mem.startsWith(u8, name, "exudate_")) return exudate_count;
    if (comptime std.mem.eql(u8, name, "soil_to_root_gas_exchange_g_per_h") or
        std.mem.eql(u8, name, "aqueous_to_gaseous_root_exchange_g_per_h") or
        std.mem.eql(u8, name, "atmosphere_to_root_gas_exchange_g_per_h"))
        return domain_layer_count * transported_root_gas_count;
    if (comptime std.mem.eql(u8, name, "ammonium_uptake_g_n_per_h") or
        std.mem.eql(u8, name, "nitrate_uptake_g_n_per_h") or
        std.mem.eql(u8, name, "phosphate_uptake_g_p_per_h") or
        std.mem.eql(u8, name, "fixation_uptake_g_n_per_h") or
        std.mem.startsWith(u8, name, "withdrawal_") or
        std.mem.startsWith(u8, name, "combustion_") or
        std.mem.eql(u8, name, "current_porosity_fraction") or
        std.mem.eql(u8, name, "initial_porosity_fraction") or
        std.mem.startsWith(u8, name, "seed_")) return plant_count;
    return domain_layer_count;
}

fn freeAllocated(state: *State, allocated_count: usize) void {
    var visited: usize = 0;
    inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
        if (visited < allocated_count) state.allocator.free(@field(state, field.name));
        visited += 1;
    };
    state.allocator.free(state.planting_layer_by_plant);
}

test "STARTQ root and mycorrhizal state has runtime axes and source initial values" {
    const traits = try @import("plant_traits.zig").parse(@import("test_fixtures.zig").plant_traits_source);
    var state = try State.init(std.testing.allocator, 2, 3, 17);
    defer state.deinit();
    try state.initializePlant(1, traits, 2, 0.08, compatibilityInitializationParameters());
    try std.testing.expectEqual(@as(usize, 17), state.root_axis_count);
    try std.testing.expectEqual(@as(usize, 2), state.planting_layer_by_plant[1]);
    const root_layer = try state.layerIndex(1, 0, 2);
    const mycorrhizal_layer = try state.layerIndex(1, 1, 2);
    try std.testing.expectEqual(@as(f64, -0.01), state.total_water_potential_mpa[root_layer]);
    try std.testing.expectEqual(@as(f64, 1.0e-3), state.active_length_m[root_layer]);
    try std.testing.expectEqual(traits.roots.primary_root_radius_m, state.primary_radius_m[root_layer]);
    try std.testing.expectEqual(@as(f64, 2.5e-6), state.primary_radius_m[mycorrhizal_layer]);
    try std.testing.expectEqual(@as(f64, 0.08), state.axis_depth_m[try state.axisIndex(1, 1, 16)]);
    try std.testing.expectEqual(traits.roots.root_porosity_fraction, state.current_porosity_fraction[1]);
    try std.testing.expectEqual(traits.roots.root_porosity_fraction, state.initial_porosity_fraction[1]);
    try state.validateFinite();
}

test "STARTQ runtime root initialization controls all domain seed values" {
    const traits = try @import("plant_traits.zig").parse(@import("test_fixtures.zig").plant_traits_source);
    var parameters = compatibilityInitializationParameters();
    parameters.mycorrhizal_radius_m = 4e-6;
    parameters.initial_total_water_potential_mpa = -0.02;
    parameters.osmotic_water_potential_decrement_mpa = -0.03;
    parameters.initial_active_length_m = 0.002;
    parameters.initial_water_fraction = 0.9;
    var state = try State.init(std.testing.allocator, 1, 2, 3);
    defer state.deinit();
    try state.initializePlant(0, traits, 1, 0.05, parameters);
    const root = try state.layerIndex(0, 0, 1);
    const mycorrhiza = try state.layerIndex(0, 1, 1);
    try std.testing.expectEqual(@as(f64, -0.02), state.total_water_potential_mpa[root]);
    try std.testing.expectEqual(@as(f64, 0.002), state.active_length_m[root]);
    try std.testing.expectEqual(@as(f64, 0.9), state.water_fraction[mycorrhiza]);
    try std.testing.expectEqual(@as(f64, 4e-6), state.primary_radius_m[mycorrhiza]);
    try std.testing.expectApproxEqAbs(traits.water_relations.osmotic_potential_mpa - 0.03, state.osmotic_water_potential_mpa[root], 1e-15);
}

test "replant reconstruction clears every root history without changing neighboring plants" {
    const traits = try @import("plant_traits.zig").parse(@import("test_fixtures.zig").plant_traits_source);
    var state = try State.init(std.testing.allocator, 2, 3, 4);
    defer state.deinit();
    @memset(state.mobile_carbon_g, 11);
    @memset(state.gaseous_oxygen_g_o, 12);
    @memset(state.salt_content_mol, 13);
    @memset(state.exudate_carbon_exchange_g_c_per_h, 14);
    @memset(state.axis_primary_carbon_g, 15);
    state.active_root_axis_count[0] = 3;
    state.roots_dead[0] = false;

    try state.reconstructPlant(0, traits, 2, 0.08, compatibilityInitializationParameters());

    const layer_extent = biological_domain_count * state.soil_layer_count;
    const salt_extent = layer_extent * salt_species_count;
    const substrate_extent = layer_extent * organic_substrate_count;
    const axis_extent = biological_domain_count * state.root_axis_count;
    for (state.mobile_carbon_g[0..layer_extent]) |value| try std.testing.expectEqual(@as(f64, 0), value);
    for (state.gaseous_oxygen_g_o[0..layer_extent]) |value| try std.testing.expectEqual(@as(f64, 0), value);
    for (state.salt_content_mol[0..salt_extent]) |value| try std.testing.expectEqual(@as(f64, 0), value);
    for (state.exudate_carbon_exchange_g_c_per_h[0..substrate_extent]) |value| try std.testing.expectEqual(@as(f64, 0), value);
    for (state.axis_primary_carbon_g[0..axis_extent]) |value| try std.testing.expectEqual(@as(f64, 0), value);
    try std.testing.expectEqual(@as(f64, 11), state.mobile_carbon_g[layer_extent]);
    try std.testing.expectEqual(@as(f64, 13), state.salt_content_mol[salt_extent]);
    try std.testing.expectEqual(@as(usize, 0), state.active_root_axis_count[0]);
    try std.testing.expect(state.roots_dead[0]);
    try std.testing.expectEqual(@as(usize, 2), state.planting_layer_by_plant[0]);
}

test "UPTAKE competition uses demand shares, source minima, and biome fallback" {
    const population = UptakeCompetition{
        .oxygen = 2,
        .ammonium_nonband = 0.01,
        .ammonium_band = 3,
        .nitrate_nonband = 4,
        .nitrate_band = 5,
        .phosphate_h2_nonband = 6,
        .phosphate_h2_band = 7,
        .phosphate_h_nonband = 8,
        .phosphate_h_band = 9,
    };
    const combined = UptakeCompetition{
        .oxygen = 10,
        .ammonium_nonband = 10,
        .ammonium_band = 0,
        .nitrate_nonband = 20,
        .nitrate_band = 20,
        .phosphate_h2_nonband = 20,
        .phosphate_h2_band = 20,
        .phosphate_h_nonband = 20,
        .phosphate_h_band = 20,
    };
    const fractions = try uptakeCompetition(0.05, 0.125, 1.0e-12, population, combined);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), fractions.oxygen, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), fractions.ammonium_nonband, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.125), fractions.ammonium_band, 1.0e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.45), fractions.phosphate_h_band, 1.0e-12);
}

test "root salt state uses runtime plant domain layer and typed species extents" {
    var state = try State.init(std.testing.allocator, 7, 3, 2);
    defer state.deinit();
    try std.testing.expectEqual(@as(usize, 7 * biological_domain_count * 3 * salt_species_count), state.salt_content_mol.len);
    const chloride = try state.saltIndex(6, 1, 2, .chloride);
    state.salt_content_mol[chloride] = 0.25;
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), state.salt_content_mol[state.salt_content_mol.len - 1], 1.0e-12);
}

test "hourly root reset clears fluxes but preserves pools and morphology" {
    var state = try State.init(std.testing.allocator, 1, 1, 1);
    defer state.deinit();
    state.mobile_carbon_g[0] = 3;
    state.root_length_m_per_plant[0] = 2;
    state.water_uptake_m3_per_h[0] = -0.1;
    state.ammonium_uptake_nonband_g_n_per_h[0] = 0.03;
    state.ammonium_demand_nonband_g_n_per_h[0] = 0.05;
    state.salt_uptake_mol_per_h[0] = 0.2;
    state.exudate_carbon_exchange_g_c_per_h[0] = -0.3;
    state.symbiont_structural_carbon_g_c[0] = 0.4;
    state.symbiont_mobile_nitrogen_g_n[0] = 0.02;
    state.symbiotic_respiration_actual_g_c_per_h[0] = 0.03;
    state.symbiotic_respiration_oxygen_unlimited_g_c_per_h[0] = 0.05;
    state.withdrawal_carbon_dioxide_loss_g_c_per_h[0] = -0.2;
    state.resetHourlyFluxes();
    try std.testing.expectEqual(@as(f64, 0), state.water_uptake_m3_per_h[0]);
    try std.testing.expectEqual(@as(f64, 0.03), state.previous_ammonium_uptake_nonband_g_n_per_h[0]);
    try std.testing.expectEqual(@as(f64, 0), state.ammonium_uptake_nonband_g_n_per_h[0]);
    try std.testing.expectEqual(@as(f64, 0.05), state.previous_ammonium_demand_nonband_g_n_per_h[0]);
    try std.testing.expectEqual(@as(f64, 0), state.ammonium_demand_nonband_g_n_per_h[0]);
    try std.testing.expectEqual(@as(f64, 0), state.salt_uptake_mol_per_h[0]);
    try std.testing.expectEqual(@as(f64, 0), state.exudate_carbon_exchange_g_c_per_h[0]);
    try std.testing.expectEqual(@as(f64, 3), state.mobile_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0.4), state.symbiont_structural_carbon_g_c[0]);
    try std.testing.expectEqual(@as(f64, 0.02), state.symbiont_mobile_nitrogen_g_n[0]);
    try std.testing.expectEqual(@as(f64, 0), state.symbiotic_respiration_actual_g_c_per_h[0]);
    try std.testing.expectEqual(@as(f64, 0), state.symbiotic_respiration_oxygen_unlimited_g_c_per_h[0]);
    try std.testing.expectEqual(@as(f64, 0), state.withdrawal_carbon_dioxide_loss_g_c_per_h[0]);
    try std.testing.expectEqual(@as(f64, 2), state.root_length_m_per_plant[0]);
}

test "GROSUB root topology aggregation has no axis ceiling" {
    var state = try State.init(std.testing.allocator, 1, 2, 12);
    defer state.deinit();
    for (0..12) |axis| {
        const index = try state.layerAxisIndex(0, 0, 1, axis);
        state.axis_primary_count[index] = 1;
        state.axis_secondary_count[index] = 2;
        state.axis_primary_length_m[index] = 0.1;
        state.axis_secondary_length_m[index] = 0.2;
    }
    const topology = try state.layerTopology(0, 0, 1, 6, 0.5, 0.75, 0.01);
    try std.testing.expectApproxEqAbs(12, topology.primary_axis_count, 1e-14);
    try std.testing.expectApproxEqAbs(24, topology.secondary_axis_count, 1e-14);
    try std.testing.expectApproxEqAbs(9.6, topology.total_root_length_m, 1e-14);
    try std.testing.expectApproxEqAbs(1.2, topology.root_length_m_per_plant, 1e-14);
    try std.testing.expectApproxEqAbs(2.4, topology.root_length_density_m_per_m3, 1e-14);
    try std.testing.expectApproxEqAbs(0.1, topology.average_secondary_length_m, 1e-14);
    const layer = try state.layerIndex(0, 0, 1);
    state.total_carbon_g[layer] = 2;
    state.turgor_water_potential_mpa[layer] = 1;
    state.reference_primary_radius_m[layer] = 4e-4;
    state.reference_secondary_radius_m[layer] = 2e-4;
    const parameters = compatibilityMorphologyParameters();
    const refreshed = try state.refreshLayerMorphology(0, 0, 1, 6, 0.5, 0.75, 0.2, 2e-5, 3.1416, true, parameters);
    try std.testing.expectEqual(topology.total_root_length_m, refreshed.total_root_length_m);
    try std.testing.expectApproxEqAbs(@as(f64, 4.8e-4), state.primary_radius_m[layer], 1e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 2.4e-4), state.secondary_radius_m[layer], 1e-14);
    try std.testing.expectApproxEqAbs(0.8 * @max(3.1416 * 2e-4 * 2e-4 * 2.4, 4e-5), state.aqueous_volume_m3[layer], 1e-14);
    try std.testing.expect(state.root_surface_area_m2_per_plant[layer] > 0);
}

test "STARTQ seed geometry augments only the root planting layer" {
    var state = try State.init(std.testing.allocator, 1, 2, 1);
    defer state.deinit();
    state.planting_layer_by_plant[0] = 1;
    try state.setSeedGeometry(0, 1.0e-6, 0.02, 0.003);

    const root_layer = try state.layerIndex(0, 0, 1);
    state.reference_primary_radius_m[root_layer] = 4.0e-4;
    state.reference_secondary_radius_m[root_layer] = 2.0e-4;
    state.turgor_water_potential_mpa[root_layer] = 1;
    const topology = try state.refreshLayerMorphology(0, 0, 1, 6, 0.5, 0.75, 0.2, 2.0e-5, 3.1416, true, compatibilityMorphologyParameters());
    try std.testing.expectApproxEqAbs(@as(f64, 0.02), topology.root_length_m_per_plant, 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.04), topology.root_length_density_m_per_m3, 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 6.0e-6), state.aqueous_volume_m3[root_layer] + state.gaseous_volume_m3[root_layer], 1.0e-14);
    try std.testing.expectApproxEqAbs(@as(f64, 0.003), state.root_surface_area_m2_per_plant[root_layer], 1.0e-14);

    const symbiont = try state.layerTopology(0, 1, 1, 6, 0.5, 0.75, 0.01);
    const other_layer = try state.layerTopology(0, 0, 0, 6, 0.5, 0.75, 0.01);
    try std.testing.expectEqual(@as(f64, 0), symbiont.root_length_m_per_plant);
    try std.testing.expectEqual(@as(f64, 0), other_layer.root_length_m_per_plant);
}

test "GROSUB empty root layer releases both gas phases and restores reference radii" {
    var state = try State.init(std.testing.allocator, 1, 1, 1);
    defer state.deinit();
    const root = try state.layerIndex(0, 0, 0);
    state.reference_primary_radius_m[root] = 4.0e-4;
    state.reference_secondary_radius_m[root] = 2.0e-4;
    state.primary_radius_m[root] = 8.0e-4;
    state.secondary_radius_m[root] = 6.0e-4;
    state.gaseous_carbon_dioxide_g_c[root] = 2;
    state.aqueous_carbon_dioxide_g_c[root] = 3;
    state.gaseous_oxygen_g_o[root] = 5;
    state.aqueous_oxygen_g_o[root] = 7;

    _ = try state.refreshLayerMorphology(0, 0, 0, 4, 0.25, 1, 0.2, 2.0e-5, 3.142, false, compatibilityMorphologyParameters());

    try std.testing.expectEqual(@as(f64, 4.0e-4), state.primary_radius_m[root]);
    try std.testing.expectEqual(@as(f64, 2.0e-4), state.secondary_radius_m[root]);
    try std.testing.expectEqual(@as(f64, -5), state.withdrawal_carbon_dioxide_loss_g_c_per_h[0]);
    try std.testing.expectEqual(@as(f64, -12), state.withdrawal_oxygen_loss_g_o_per_h[0]);
    try std.testing.expectEqual(@as(f64, 0), state.gaseous_carbon_dioxide_g_c[root]);
    try std.testing.expectEqual(@as(f64, 0), state.aqueous_oxygen_g_o[root]);
}

test "GROSUB active carbon assigns whole primary axis to its runtime tip layer" {
    var state = try State.init(std.testing.allocator, 1, 3, 2);
    defer state.deinit();
    state.active_root_axis_count[0] = 2;
    state.planting_layer_by_plant[0] = 0;
    const axis_0_layer_0 = try state.layerAxisIndex(0, 0, 0, 0);
    const axis_0_layer_1 = try state.layerAxisIndex(0, 0, 1, 0);
    const axis_1_layer_0 = try state.layerAxisIndex(0, 0, 0, 1);
    const axis_1_layer_2 = try state.layerAxisIndex(0, 0, 2, 1);
    state.axis_primary_carbon_g[axis_0_layer_0] = 1;
    state.axis_primary_carbon_g[axis_0_layer_1] = 2;
    state.axis_primary_length_m[axis_0_layer_1] = 0.1;
    state.axis_primary_carbon_g[axis_1_layer_0] = 4;
    state.axis_primary_carbon_g[axis_1_layer_2] = 5;
    state.axis_primary_length_m[axis_1_layer_2] = 0.1;
    state.axis_secondary_carbon_g[axis_0_layer_0] = 0.5;
    state.axis_secondary_carbon_g[axis_1_layer_2] = 0.25;
    try state.refreshActiveCarbonByLayer(0, 0);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), state.total_carbon_g[try state.layerIndex(0, 0, 0)], 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 3), state.total_carbon_g[try state.layerIndex(0, 0, 1)], 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 9.25), state.total_carbon_g[try state.layerIndex(0, 0, 2)], 1.0e-15);
}

test "GROSUB root totals distinguish active and actual layer carbon" {
    const totals = try State.sourceOrderLayerCarbonTotals(
        &.{ 1.0, 3.0, 5.0 },
        &.{ 0.5, 0.25, 0.125 },
    );
    try std.testing.expectApproxEqAbs(@as(f64, 0.875), totals.active_secondary_carbon_g_c, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 9.875), totals.actual_primary_and_secondary_carbon_g_c, 1.0e-15);
    try std.testing.expectError(error.RootCarbonTotalDimensionMismatch, State.sourceOrderLayerCarbonTotals(&.{1}, &.{ 1, 2 }));
}

test "GROSUB negative structural carbon cleanup charges the mobile pool" {
    const primary = try sourceOrderNegativeStructuralCarbonCleanup(-0.25, 2);
    try std.testing.expectEqual(@as(f64, 0), primary.structural_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 1.75), primary.mobile_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0.25), primary.removed_deficit_g_c);

    const secondary = try sourceOrderNegativeStructuralCarbonCleanup(0.5, 2);
    try std.testing.expectEqual(@as(f64, 0.5), secondary.structural_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 2), secondary.mobile_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), secondary.removed_deficit_g_c);
}

test "GROSUB negative cleanup exposes rather than hides a mobile overdraw" {
    const result = try sourceOrderNegativeStructuralCarbonCleanup(-2, 0.5);
    try std.testing.expectEqual(@as(f64, -1.5), result.mobile_carbon_g_c);
    try std.testing.expectError(
        error.NonFiniteNegativeStructuralCarbonCleanup,
        sourceOrderNegativeStructuralCarbonCleanup(std.math.nan(f64), 1),
    );
}

test "GROSUB morphology requires length carbon and population above ZEROP" {
    try std.testing.expect(try sourceOrderRootMorphologyIsActive(1, 2, 3, 1.0e-6));
    try std.testing.expect(!try sourceOrderRootMorphologyIsActive(1.0e-6, 2, 3, 1.0e-6));
    try std.testing.expect(!try sourceOrderRootMorphologyIsActive(1, 2, 0, 1.0e-6));
}

test "GROSUB thin layer publishes zero root length density" {
    try std.testing.expectEqual(@as(f64, 4), try sourceOrderRootLengthDensity(1, 0.25, 0.01));
    try std.testing.expectEqual(@as(f64, 0), try sourceOrderRootLengthDensity(1, 0.01, 0.01));
    try std.testing.expectEqual(@as(f64, 0), try sourceOrderRootLengthDensity(1, 0, 0.01));
}
