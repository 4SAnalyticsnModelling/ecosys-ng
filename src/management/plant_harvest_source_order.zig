//! Split out of plant_harvest_runtime.zig by tools/split_decl_group.py.
//! Pure code motion: every decl below is an exact line slice.

const std = @import("std");
const management = @import("plant_management.zig");
const canopy = @import("../canopy/photosynthesis/photosynthesis.zig");
const phenology = @import("../plant/lifecycle/phenology.zig");
const litter_partition = @import("../plant/partition/litter.zig");
const grazing_manure = @import("grazing_manure.zig");
const __parent = @import("plant_harvest_runtime.zig");
const HourlyDisturbanceReset = __parent.HourlyDisturbanceReset;
const PopulationAfterDisturbance = __parent.PopulationAfterDisturbance;
const PopulationScaledNumericalThresholds = __parent.PopulationScaledNumericalThresholds;
const TillageElementComposition = __parent.TillageElementComposition;
const boundedCombustionFraction = __parent.boundedCombustionFraction;
const harvest_product_component_count = __parent.harvest_product_component_count;
const remainingShootPools = __parent.remainingShootPools;
const scaleElementalMass = __parent.scaleElementalMass;
const scaleSaltInventory = __parent.scaleSaltInventory;
const source_order_standing_dead_component_count = __parent.source_order_standing_dead_component_count;
const subtractElementalMass = __parent.subtractElementalMass;
const subtractSaltInventory = __parent.subtractSaltInventory;
const sumHarvestProductComponents = __parent.sumHarvestProductComponents;
const uncombustedShootTotal = __parent.uncombustedShootTotal;
const validateCombustionNode = __parent.validateCombustionNode;
const validateCompleteDeathMass = __parent.validateCompleteDeathMass;
const validateCompleteDeathResultMass = __parent.validateCompleteDeathResultMass;
const validateDeadNoduleMass = __parent.validateDeadNoduleMass;
const validateDeadRootDepthMass = __parent.validateDeadRootDepthMass;
const validateDeadRootResetMass = __parent.validateDeadRootResetMass;
const validateSourceOrderRootMass = __parent.validateSourceOrderRootMass;
const validateTillageComposition = __parent.validateTillageComposition;
const validateTillageSlice = __parent.validateTillageSlice;

/// Exact GROSUB 8516-8524 first-substep disturbance reset.
pub fn sourceOrderHourlyDisturbanceReset(
    first_biological_substep: bool,
    cumulative_harvest_carbon_g_c: f64,
    current: HourlyDisturbanceReset,
) !HourlyDisturbanceReset {
    if (!std.math.isFinite(cumulative_harvest_carbon_g_c) or cumulative_harvest_carbon_g_c < 0)
        return error.InvalidHourlyDisturbanceResetInput;
    inline for (@typeInfo(HourlyDisturbanceReset).@"struct".fields) |field| {
        const value = @field(current, field.name);
        if (field.type == f64) {
            if (!std.math.isFinite(value) or value < 0) return error.InvalidHourlyDisturbanceResetInput;
        } else for (value) |item| {
            if (!std.math.isFinite(item) or item < 0) return error.InvalidHourlyDisturbanceResetInput;
        }
    }
    if (!first_biological_substep) return current;
    return .{
        .previous_cumulative_harvest_carbon_g_c = cumulative_harvest_carbon_g_c,
        .manure_organic_carbon_g_c = @splat(0),
        .manure_organic_nitrogen_g_n = @splat(0),
        .manure_organic_phosphorus_g_p = @splat(0),
        .manure_inorganic_nitrogen_g_n = 0,
        .manure_inorganic_phosphorus_g_p = 0,
    };
}

/// Exact GROSUB 8528-8531 calendar, trait, and management selector.
pub fn sourceOrderForestSelfThinningIsEnabled(
    day_of_year: u16,
    hour_of_day: u8,
    local_solar_noon_h: f64,
    biomass_turnover_type: u8,
    root_profile_type: u8,
    harvest_code: i8,
) !bool {
    if (day_of_year == 0 or day_of_year > 366 or hour_of_day > 23 or
        !std.math.isFinite(local_solar_noon_h) or local_solar_noon_h < 0 or local_solar_noon_h >= 24)
        return error.InvalidForestSelfThinningSelector;
    return day_of_year % 30 == 0 and
        hour_of_day == @as(u8, @intFromFloat(@floor(local_solar_noon_h))) and
        biomass_turnover_type != 0 and
        root_profile_type > 1 and
        (harvest_code < 0 or harvest_code == 4 or harvest_code == 6);
}

pub const SourceOrderTillagePopulationState = struct {
    living_population_per_m2: f64,
    living_population_count: f64,
    standing_dead_population_count: f64,
    canopy_radiation_fraction: f64,
};

pub const SourceOrderTillagePopulationInput = struct {
    hour_of_day: u8,
    local_solar_noon_h: f64,
    biomass_turnover_type: u8,
    root_profile_type: u8,
    current_day_of_year: u16,
    current_year: u32,
    planting_day_of_year: u16,
    planting_year: u32,
    tillage_code: u8,
    is_first_plant_population: bool,
    remaining_fraction: f64,
    zero_population_threshold: f64,
    state: SourceOrderTillagePopulationState,
};

pub const SourceOrderTillagePopulationResult = struct {
    state: SourceOrderTillagePopulationState,
    applied: bool,
    clear_leaf_sheath_and_sapwood_totals: bool,
    terminate_living_branches: bool,
};

/// Exact GROSUB 10031-10052 tillage population selector and scaling.
pub fn sourceOrderTillagePopulationReduction(
    input: SourceOrderTillagePopulationInput,
) !SourceOrderTillagePopulationResult {
    if (input.hour_of_day > 23 or
        !std.math.isFinite(input.local_solar_noon_h) or
        input.local_solar_noon_h < 0 or input.local_solar_noon_h >= 24 or
        input.current_day_of_year == 0 or input.current_day_of_year > 366 or
        input.planting_day_of_year == 0 or input.planting_day_of_year > 366 or
        input.current_year == 0 or input.planting_year == 0 or
        !std.math.isFinite(input.remaining_fraction) or
        input.remaining_fraction < 0 or input.remaining_fraction > 1 or
        !std.math.isFinite(input.zero_population_threshold) or
        input.zero_population_threshold < 0)
        return error.InvalidTillagePopulationInput;
    inline for (@typeInfo(SourceOrderTillagePopulationState).@"struct".fields) |field| {
        const value = @field(input.state, field.name);
        if (!std.math.isFinite(value) or value < 0) return error.InvalidTillagePopulationInput;
    }

    const at_solar_noon = input.hour_of_day ==
        @as(u8, @intFromFloat(@floor(input.local_solar_noon_h)));
    const herbaceous = input.biomass_turnover_type == 0 or input.root_profile_type <= 1;
    const not_planting_date = input.current_day_of_year != input.planting_day_of_year or
        input.current_year != input.planting_year;
    const valid_tillage_code = input.tillage_code > 0 and input.tillage_code <= 20;
    const includes_population = input.tillage_code <= 10 or !input.is_first_plant_population;
    // Preserve the source's disjunctive comparison, including its
    // non-chronological behavior when schedule years are out of sequence.
    const after_planting = input.current_day_of_year > input.planting_day_of_year or
        input.current_year > input.planting_year;
    const applied = at_solar_noon and herbaceous and not_planting_date and
        valid_tillage_code and includes_population and after_planting;
    if (!applied) return .{
        .state = input.state,
        .applied = false,
        .clear_leaf_sheath_and_sapwood_totals = false,
        .terminate_living_branches = false,
    };

    var state = input.state;
    inline for (@typeInfo(SourceOrderTillagePopulationState).@"struct".fields) |field|
        @field(state, field.name) *= input.remaining_fraction;
    return .{
        .state = state,
        .applied = true,
        .clear_leaf_sheath_and_sapwood_totals = true,
        .terminate_living_branches = state.living_population_count <= input.zero_population_threshold,
    };
}

pub const SourceOrderTillageBranchPools = struct {
    host_mobile: canopy.ElementalMass,
    symbiont_mobile: canopy.ElementalMass,
    c4_mobile_carbon_g_c: f64,
    stalk_reserve: canopy.ElementalMass,
    leaf: canopy.ElementalMass,
    symbiont_structural: canopy.ElementalMass,
    sheath: canopy.ElementalMass,
    husk: canopy.ElementalMass,
    ear: canopy.ElementalMass,
    grain: canopy.ElementalMass,
    stalk: canopy.ElementalMass,
};

pub const SourceOrderTillageBranchLitterInput = struct {
    remaining_fraction: f64,
    winter_annual: bool,
    pools: SourceOrderTillageBranchPools,
    leaf_composition: TillageElementComposition,
    sheath_composition: TillageElementComposition,
    stalk_composition: TillageElementComposition,
    nonstructural_kinetics: litter_partition.ElementFractions,
    foliar_kinetics: litter_partition.ElementFractions,
    nonfoliar_kinetics: litter_partition.ElementFractions,
    stalk_kinetics: litter_partition.ElementFractions,
    coarse_wood_kinetics: litter_partition.ElementFractions,
};

pub const SourceOrderTillageBranchLitterResult = struct {
    litter: canopy.SenescenceProducts,
    seasonal_storage: canopy.ElementalMass,
};

/// Exact GROSUB 10096-10157 branch litter allocation during tillage.
pub fn sourceOrderTillageBranchLitter(
    input: SourceOrderTillageBranchLitterInput,
) !SourceOrderTillageBranchLitterResult {
    if (!std.math.isFinite(input.remaining_fraction) or
        input.remaining_fraction < 0 or input.remaining_fraction > 1 or
        !std.math.isFinite(input.pools.c4_mobile_carbon_g_c) or
        input.pools.c4_mobile_carbon_g_c < 0)
        return error.InvalidTillageBranchLitterInput;
    inline for (@typeInfo(SourceOrderTillageBranchPools).@"struct".fields) |field| {
        if (field.type == canopy.ElementalMass) {
            const mass = @field(input.pools, field.name);
            inline for (.{ mass.carbon_g, mass.nitrogen_g, mass.phosphorus_g }) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.InvalidTillageBranchLitterInput;
        }
    }
    try validateTillageComposition(input.leaf_composition);
    try validateTillageComposition(input.sheath_composition);
    try validateTillageComposition(input.stalk_composition);
    inline for (.{
        input.nonstructural_kinetics,
        input.foliar_kinetics,
        input.nonfoliar_kinetics,
        input.stalk_kinetics,
        input.coarse_wood_kinetics,
    }) |kinetics| kinetics.validate() catch return error.InvalidTillageBranchLitterInput;

    const removed = 1 - input.remaining_fraction;
    const pools = input.pools;
    var result: SourceOrderTillageBranchLitterResult = .{
        .litter = .{},
        .seasonal_storage = .{},
    };
    for (0..litter_partition.kinetic_component_count) |kinetic| {
        const nonwoody_carbon_g_c = input.nonstructural_kinetics.carbon[kinetic] *
            (pools.host_mobile.carbon_g + pools.symbiont_mobile.carbon_g +
                pools.c4_mobile_carbon_g_c + pools.stalk_reserve.carbon_g) +
            input.foliar_kinetics.carbon[kinetic] *
                (pools.leaf.carbon_g * input.leaf_composition.carbon[1] +
                    pools.symbiont_structural.carbon_g) +
            input.nonfoliar_kinetics.carbon[kinetic] *
                (pools.sheath.carbon_g * input.sheath_composition.carbon[1] +
                    pools.husk.carbon_g + pools.ear.carbon_g) +
            input.stalk_kinetics.carbon[kinetic] *
                pools.stalk.carbon_g * input.stalk_composition.carbon[1];
        const nonwoody_nitrogen_g_n = input.nonstructural_kinetics.nitrogen[kinetic] *
            (pools.host_mobile.nitrogen_g + pools.symbiont_mobile.nitrogen_g +
                pools.stalk_reserve.nitrogen_g) +
            input.foliar_kinetics.nitrogen[kinetic] *
                (pools.leaf.nitrogen_g * input.leaf_composition.nitrogen[1] +
                    pools.symbiont_structural.nitrogen_g) +
            input.nonfoliar_kinetics.nitrogen[kinetic] *
                (pools.sheath.nitrogen_g * input.sheath_composition.nitrogen[1] +
                    pools.husk.nitrogen_g + pools.ear.nitrogen_g) +
            input.stalk_kinetics.nitrogen[kinetic] *
                pools.stalk.nitrogen_g * input.stalk_composition.nitrogen[1];
        const nonwoody_phosphorus_g_p = input.nonstructural_kinetics.phosphorus[kinetic] *
            (pools.host_mobile.phosphorus_g + pools.symbiont_mobile.phosphorus_g +
                pools.stalk_reserve.phosphorus_g) +
            input.foliar_kinetics.phosphorus[kinetic] *
                (pools.leaf.phosphorus_g * input.leaf_composition.phosphorus[1] +
                    pools.symbiont_structural.phosphorus_g) +
            input.nonfoliar_kinetics.phosphorus[kinetic] *
                (pools.sheath.phosphorus_g * input.sheath_composition.phosphorus[1] +
                    pools.husk.phosphorus_g + pools.ear.phosphorus_g) +
            input.stalk_kinetics.phosphorus[kinetic] *
                pools.stalk.phosphorus_g * input.stalk_composition.phosphorus[1];
        result.litter.nonwoody_carbon_g[kinetic] = removed * nonwoody_carbon_g_c;
        result.litter.nonwoody_nitrogen_g[kinetic] = removed * nonwoody_nitrogen_g_n;
        result.litter.nonwoody_phosphorus_g[kinetic] = removed * nonwoody_phosphorus_g_p;
        result.litter.woody_carbon_g[kinetic] = removed * input.coarse_wood_kinetics.carbon[kinetic] *
            (pools.leaf.carbon_g * input.leaf_composition.carbon[0] +
                pools.sheath.carbon_g * input.sheath_composition.carbon[0] +
                pools.stalk.carbon_g * input.stalk_composition.carbon[0]);
        result.litter.woody_nitrogen_g[kinetic] = removed * input.coarse_wood_kinetics.nitrogen[kinetic] *
            (pools.leaf.nitrogen_g * input.leaf_composition.nitrogen[0] +
                pools.sheath.nitrogen_g * input.sheath_composition.nitrogen[0] +
                pools.stalk.nitrogen_g * input.stalk_composition.nitrogen[0]);
        result.litter.woody_phosphorus_g[kinetic] = removed * input.coarse_wood_kinetics.phosphorus[kinetic] *
            (pools.leaf.phosphorus_g * input.leaf_composition.phosphorus[0] +
                pools.sheath.phosphorus_g * input.sheath_composition.phosphorus[0] +
                pools.stalk.phosphorus_g * input.stalk_composition.phosphorus[0]);
        if (input.winter_annual) {
            result.seasonal_storage.carbon_g += removed * input.nonfoliar_kinetics.carbon[kinetic] * pools.grain.carbon_g;
            result.seasonal_storage.nitrogen_g += removed * input.nonfoliar_kinetics.nitrogen[kinetic] * pools.grain.nitrogen_g;
            result.seasonal_storage.phosphorus_g += removed * input.nonfoliar_kinetics.phosphorus[kinetic] * pools.grain.phosphorus_g;
        } else {
            result.litter.nonwoody_carbon_g[kinetic] += removed * input.nonfoliar_kinetics.carbon[kinetic] * pools.grain.carbon_g;
            result.litter.nonwoody_nitrogen_g[kinetic] += removed * input.nonfoliar_kinetics.nitrogen[kinetic] * pools.grain.nitrogen_g;
            result.litter.nonwoody_phosphorus_g[kinetic] += removed * input.nonfoliar_kinetics.phosphorus[kinetic] * pools.grain.phosphorus_g;
        }
    }
    return result;
}

pub const SourceOrderTillageBranchScalarState = struct {
    host_mobile: canopy.ElementalMass,
    c4_mobile_carbon_g_c: f64,
    symbiont_mobile: canopy.ElementalMass,
    total_shoot: canopy.ElementalMass,
    leaf: canopy.ElementalMass,
    symbiont_structural: canopy.ElementalMass,
    sheath: canopy.ElementalMass,
    stalk: canopy.ElementalMass,
    sapwood_carbon_g_c: f64,
    stalk_reserve: canopy.ElementalMass,
    husk: canopy.ElementalMass,
    ear: canopy.ElementalMass,
    grain: canopy.ElementalMass,
    potential_seed_site_count: f64,
    seed_count: f64,
    individual_seed_carbon_g_c: f64,
    leaf_area_m2: f64,
    stalk_total: canopy.ElementalMass,
};

pub const SourceOrderTillageNodeState = struct {
    c3_mobile_carbon_g_c: []f64,
    c4_mobile_carbon_g_c: []f64,
    carbon_dioxide_g_c: []f64,
    bicarbonate_g_c: []f64,
    leaf_area_m2: []f64,
    growing_leaf_carbon_g_c: []f64,
    senescing_leaf_carbon_g_c: []f64,
    growing_sheath_carbon_g_c: []f64,
    senescing_sheath_carbon_g_c: []f64,
    growing_node_carbon_g_c: []f64,
    growing_leaf_nitrogen_g_n: []f64,
    growing_sheath_nitrogen_g_n: []f64,
    growing_node_nitrogen_g_n: []f64,
    growing_leaf_phosphorus_g_p: []f64,
    growing_sheath_phosphorus_g_p: []f64,
    growing_node_phosphorus_g_p: []f64,
};

pub const SourceOrderTillageLayerSampleState = struct {
    leaf_area_m2: []f64,
    growing_leaf_carbon_g_c: []f64,
    growing_leaf_nitrogen_g_n: []f64,
    growing_leaf_phosphorus_g_p: []f64,
};

pub const SourceOrderTillageBranchRetentionResult = struct {
    leaf_sheath_carbon_g_c: f64,
    sapwood_carbon_g_c: f64,
};

/// Exact GROSUB 10161-10235 branch state remaining after tillage.
pub fn sourceOrderRetainTillageBranchState(
    scalar: *SourceOrderTillageBranchScalarState,
    nodes: SourceOrderTillageNodeState,
    layer_samples: SourceOrderTillageLayerSampleState,
    remaining_fraction: f64,
) !SourceOrderTillageBranchRetentionResult {
    if (!std.math.isFinite(remaining_fraction) or remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidTillageBranchState;
    inline for (@typeInfo(SourceOrderTillageBranchScalarState).@"struct".fields) |field| {
        if (field.type == canopy.ElementalMass) {
            const mass = @field(scalar.*, field.name);
            inline for (.{ mass.carbon_g, mass.nitrogen_g, mass.phosphorus_g }) |value|
                if (!std.math.isFinite(value) or value < 0) return error.InvalidTillageBranchState;
        } else {
            const value = @field(scalar.*, field.name);
            if (!std.math.isFinite(value) or value < 0) return error.InvalidTillageBranchState;
        }
    }
    const node_count = nodes.leaf_area_m2.len;
    inline for (@typeInfo(SourceOrderTillageNodeState).@"struct".fields) |field|
        if (@field(nodes, field.name).len != node_count) return error.TillageBranchNodeDimensionMismatch;
    const sample_count = layer_samples.leaf_area_m2.len;
    inline for (@typeInfo(SourceOrderTillageLayerSampleState).@"struct".fields) |field|
        if (@field(layer_samples, field.name).len != sample_count) return error.TillageBranchSampleDimensionMismatch;
    inline for (@typeInfo(SourceOrderTillageNodeState).@"struct".fields) |field|
        try validateTillageSlice(@field(nodes, field.name));
    inline for (@typeInfo(SourceOrderTillageLayerSampleState).@"struct".fields) |field|
        try validateTillageSlice(@field(layer_samples, field.name));

    const seed_carbon_g_c = scalar.individual_seed_carbon_g_c;
    inline for (@typeInfo(SourceOrderTillageBranchScalarState).@"struct".fields) |field| {
        if (!std.mem.eql(u8, field.name, "individual_seed_carbon_g_c")) {
            if (field.type == canopy.ElementalMass) {
                inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element|
                    @field(@field(scalar.*, field.name), element.name) *= remaining_fraction;
            } else {
                @field(scalar.*, field.name) *= remaining_fraction;
            }
        }
    }
    scalar.individual_seed_carbon_g_c = seed_carbon_g_c;

    inline for (@typeInfo(SourceOrderTillageNodeState).@"struct".fields, 0..) |field, index| {
        const first: usize = if (index < 4) 1 else 0;
        for (@field(nodes, field.name)[first..]) |*value| value.* *= remaining_fraction;
    }
    inline for (@typeInfo(SourceOrderTillageLayerSampleState).@"struct".fields) |field| {
        for (@field(layer_samples, field.name)) |*value| value.* *= remaining_fraction;
    }

    return .{
        .leaf_sheath_carbon_g_c = @max(0, scalar.leaf.carbon_g + scalar.sheath.carbon_g),
        .sapwood_carbon_g_c = scalar.sapwood_carbon_g_c,
    };
}

pub const SourceOrderTillageStandingDeadInput = struct {
    remaining_fraction: f64,
    standing_dead_by_source_component: [litter_partition.kinetic_component_count]canopy.ElementalMass,
    composition: TillageElementComposition,
    stalk_kinetics: litter_partition.ElementFractions,
    coarse_wood_kinetics: litter_partition.ElementFractions,
};

pub const SourceOrderTillageStandingDeadResult = struct {
    litter: canopy.SenescenceProducts,
    remaining_by_source_component: [litter_partition.kinetic_component_count]canopy.ElementalMass,
};

/// Exact GROSUB 10250-10272 standing-dead litter and retention during tillage.
pub fn sourceOrderTillageStandingDead(
    input: SourceOrderTillageStandingDeadInput,
) !SourceOrderTillageStandingDeadResult {
    if (!std.math.isFinite(input.remaining_fraction) or
        input.remaining_fraction < 0 or input.remaining_fraction > 1)
        return error.InvalidTillageStandingDeadInput;
    try validateTillageComposition(input.composition);
    input.stalk_kinetics.validate() catch return error.InvalidTillageStandingDeadInput;
    input.coarse_wood_kinetics.validate() catch return error.InvalidTillageStandingDeadInput;

    var total: canopy.ElementalMass = .{};
    for (input.standing_dead_by_source_component) |mass| {
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
            const value = @field(mass, field.name);
            if (!std.math.isFinite(value) or value < 0) return error.InvalidTillageStandingDeadInput;
            @field(total, field.name) += value;
        }
    }
    const removed = 1 - input.remaining_fraction;
    var result: SourceOrderTillageStandingDeadResult = .{
        .litter = .{},
        .remaining_by_source_component = input.standing_dead_by_source_component,
    };
    for (0..litter_partition.kinetic_component_count) |kinetic| {
        result.litter.woody_carbon_g[kinetic] = removed *
            input.coarse_wood_kinetics.carbon[kinetic] * total.carbon_g * input.composition.carbon[0];
        result.litter.woody_nitrogen_g[kinetic] = removed *
            input.coarse_wood_kinetics.nitrogen[kinetic] * total.nitrogen_g * input.composition.nitrogen[0];
        result.litter.woody_phosphorus_g[kinetic] = removed *
            input.coarse_wood_kinetics.phosphorus[kinetic] * total.phosphorus_g * input.composition.phosphorus[0];
        result.litter.nonwoody_carbon_g[kinetic] = removed *
            input.stalk_kinetics.carbon[kinetic] * total.carbon_g * input.composition.carbon[1];
        result.litter.nonwoody_nitrogen_g[kinetic] = removed *
            input.stalk_kinetics.nitrogen[kinetic] * total.nitrogen_g * input.composition.nitrogen[1];
        result.litter.nonwoody_phosphorus_g[kinetic] = removed *
            input.stalk_kinetics.phosphorus[kinetic] * total.phosphorus_g * input.composition.phosphorus[1];
    }
    for (&result.remaining_by_source_component) |*mass| {
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
            @field(mass, field.name) *= input.remaining_fraction;
    }
    return result;
}

pub const SourceOrderTillageTerminationState = struct {
    roots_dead: bool,
    shoots_dead: bool,
    plant_dead: bool,
    harvest_termination_code: u8,
    harvest_day_of_year: u16,
    harvest_year: u32,
};

pub const SourceOrderTillageTerminationInput = struct {
    living_population_count: f64,
    zero_population_threshold: f64,
    current_day_of_year: u16,
    current_year: u32,
    state: SourceOrderTillageTerminationState,
};

pub const SourceOrderTillageTerminationResult = struct {
    state: SourceOrderTillageTerminationState,
    terminated: bool,
};

/// Exact GROSUB 10283-10290 zero-population tillage termination.
pub fn sourceOrderTillageTermination(
    input: SourceOrderTillageTerminationInput,
) !SourceOrderTillageTerminationResult {
    if (!std.math.isFinite(input.living_population_count) or
        input.living_population_count < 0 or
        !std.math.isFinite(input.zero_population_threshold) or
        input.zero_population_threshold < 0 or
        input.current_day_of_year == 0 or input.current_day_of_year > 366 or
        input.current_year == 0 or input.state.harvest_termination_code > 2 or
        input.state.harvest_day_of_year > 366)
        return error.InvalidTillageTerminationInput;
    if (input.living_population_count > input.zero_population_threshold)
        return .{ .state = input.state, .terminated = false };
    return .{
        .state = .{
            .roots_dead = true,
            .shoots_dead = true,
            .plant_dead = true,
            .harvest_termination_code = 1,
            .harvest_day_of_year = input.current_day_of_year,
            .harvest_year = input.current_year,
        },
        .terminated = true,
    };
}

pub const SourceOrderStandingDeadHarvestInput = struct {
    harvest_code: i8,
    hour_of_day: u8,
    local_solar_noon_h: f64,
    thinning_fraction_or_specific_consumption_rate: f64,
    standing_dead_removal_fraction: f64,
    grazer_live_mass_g_per_m2: f64,
    animal_accessible_area_m2: f64,
    insect_accessible_area_m2: f64,
    standing_dead_presence_threshold_g_c: f64,
    standing_dead_by_component: [source_order_standing_dead_component_count]canopy.ElementalMass,
};

pub const SourceOrderStandingDeadHarvestResult = struct {
    retained_fraction: f64,
    harvested_fraction: f64,
    harvested: canopy.ElementalMass,
    returned_to_litter: canopy.ElementalMass,
    remaining_by_component: [source_order_standing_dead_component_count]canopy.ElementalMass,
};

/// Exact GROSUB 10508-10551 standing-dead harvest selection and transaction.
pub fn sourceOrderStandingDeadHarvest(
    input: SourceOrderStandingDeadHarvestInput,
) !SourceOrderStandingDeadHarvestResult {
    if (input.hour_of_day > 23 or
        !std.math.isFinite(input.local_solar_noon_h) or
        input.local_solar_noon_h < 0 or input.local_solar_noon_h >= 24)
        return error.InvalidStandingDeadHarvestInput;
    inline for (.{
        input.thinning_fraction_or_specific_consumption_rate,
        input.standing_dead_removal_fraction,
        input.grazer_live_mass_g_per_m2,
        input.animal_accessible_area_m2,
        input.insect_accessible_area_m2,
        input.standing_dead_presence_threshold_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidStandingDeadHarvestInput;
    if (input.standing_dead_removal_fraction > 1) return error.InvalidStandingDeadHarvestInput;
    var total_standing_dead_carbon_g_c: f64 = 0;
    for (input.standing_dead_by_component) |mass| {
        inline for (.{ mass.carbon_g, mass.nitrogen_g, mass.phosphorus_g }) |value|
            if (!std.math.isFinite(value) or value < 0) return error.InvalidStandingDeadHarvestInput;
        total_standing_dead_carbon_g_c += mass.carbon_g;
    }

    var retained_fraction: f64 = 1;
    var harvested_fraction: f64 = 1;
    const at_solar_noon = input.hour_of_day ==
        @as(u8, @intFromFloat(@floor(input.local_solar_noon_h)));
    const grazing = input.harvest_code == 4 or input.harvest_code == 6;
    if (input.harvest_code >= 0 and at_solar_noon and !grazing) {
        if (input.thinning_fraction_or_specific_consumption_rate == 0) {
            retained_fraction = @max(0, 1 - input.standing_dead_removal_fraction);
            harvested_fraction = retained_fraction;
        } else {
            retained_fraction = @max(0, 1 - input.thinning_fraction_or_specific_consumption_rate);
            harvested_fraction = if (input.harvest_code == 0)
                @max(0, 1 - input.standing_dead_removal_fraction *
                    input.thinning_fraction_or_specific_consumption_rate)
            else
                retained_fraction;
        }
    } else if (grazing and total_standing_dead_carbon_g_c > input.standing_dead_presence_threshold_g_c) {
        const accessible_area_m2 = if (input.harvest_code == 4)
            input.animal_accessible_area_m2
        else
            input.insect_accessible_area_m2;
        const demand_g_c_per_h = input.grazer_live_mass_g_per_m2 *
            input.thinning_fraction_or_specific_consumption_rate * 0.5 / 24 *
            accessible_area_m2 * input.standing_dead_removal_fraction;
        retained_fraction = @max(0, 1 - demand_g_c_per_h / total_standing_dead_carbon_g_c);
        harvested_fraction = retained_fraction;
    }

    var result: SourceOrderStandingDeadHarvestResult = .{
        .retained_fraction = retained_fraction,
        .harvested_fraction = harvested_fraction,
        .harvested = .{},
        .returned_to_litter = .{},
        .remaining_by_component = input.standing_dead_by_component,
    };
    for (&result.remaining_by_component) |*mass| {
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
            const initial = @field(mass, field.name);
            @field(result.harvested, field.name) += (1 - harvested_fraction) * initial;
            @field(result.returned_to_litter, field.name) +=
                (harvested_fraction - retained_fraction) * initial;
            @field(mass, field.name) = retained_fraction * initial;
        }
    }
    return result;
}

pub const SourceOrderHarvestResidueInput = struct {
    harvest_code: i8,
    harvested_by_component: [harvest_product_component_count]canopy.ElementalMass,
    harvested_grain: canopy.ElementalMass,
    ecosystem_export_fraction: [4]f64,
};

/// Exact GROSUB 10585-10677 routing of five harvested-product components to
/// residue before ecosystem export totals are assembled.
pub fn sourceOrderHarvestResidueRouting(
    input: SourceOrderHarvestResidueInput,
) ![harvest_product_component_count]canopy.ElementalMass {
    if (input.harvest_code < 0 or input.harvest_code > 6 or input.harvest_code == 5)
        return error.InvalidHarvestResidueCode;
    for (input.ecosystem_export_fraction) |fraction|
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.InvalidHarvestResidueInput;
    inline for (.{
        input.harvested_grain.carbon_g,
        input.harvested_grain.nitrogen_g,
        input.harvested_grain.phosphorus_g,
    }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidHarvestResidueInput;
    for (input.harvested_by_component) |mass| {
        inline for (.{ mass.carbon_g, mass.nitrogen_g, mass.phosphorus_g }) |value|
            if (!std.math.isFinite(value) or value < 0) return error.InvalidHarvestResidueInput;
    }

    var residue = input.harvested_by_component;
    if (input.harvest_code == 1) {
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
            @field(residue[2], field.name) -=
                @field(input.harvested_grain, field.name) * input.ecosystem_export_fraction[1];
            if (!std.math.isFinite(@field(residue[2], field.name)) or
                @field(residue[2], field.name) < -1.0e-12)
                return error.HarvestResidueOverdraw;
            @field(residue[2], field.name) = @max(0, @field(residue[2], field.name));
        }
        return residue;
    }

    const fraction_by_component = [harvest_product_component_count]f64{
        input.ecosystem_export_fraction[0],
        input.ecosystem_export_fraction[0],
        input.ecosystem_export_fraction[1],
        input.ecosystem_export_fraction[2],
        input.ecosystem_export_fraction[3],
    };
    for (&residue, fraction_by_component) |*mass, export_fraction| {
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
            @field(mass, field.name) *= 1 - export_fraction;
    }
    return residue;
}

pub const SourceOrderDisturbanceRemovalInput = struct {
    harvest_code: i8,
    terminate_and_reseed: bool,
    grazer_growth_yield: f64,
    grazer_respiration_fraction: f64,
    harvested_by_component: [harvest_product_component_count]canopy.ElementalMass,
    residue_by_component: [harvest_product_component_count]canopy.ElementalMass,
    direct_litter_by_component: [harvest_product_component_count]canopy.ElementalMass,
};

pub const SourceOrderDisturbanceRemovalResult = struct {
    harvested_total: canopy.ElementalMass,
    residue_total: canopy.ElementalMass,
    direct_litter_total: canopy.ElementalMass,
    plant_ecosystem_removal: canopy.ElementalMass,
    grid_ecosystem_removal: canopy.ElementalMass,
    reseed_storage_addition: canopy.ElementalMass,
    net_biome_production_carbon_change_g_c_per_h: f64,
    plant_total_respiration_change_g_c_per_h: f64,
    plant_actual_respiration_change_g_c_per_h: f64,
    ecosystem_respiration_change_g_c_per_h: f64,
    autotrophic_respiration_change_g_c_per_h: f64,
};

/// Exact GROSUB 10699-10750 five-component disturbance totals and accounting.
pub fn sourceOrderTotalDisturbanceRemoval(
    input: SourceOrderDisturbanceRemovalInput,
) !SourceOrderDisturbanceRemovalResult {
    if (input.harvest_code < 0 or input.harvest_code > 6 or input.harvest_code == 5 or
        !std.math.isFinite(input.grazer_growth_yield) or input.grazer_growth_yield < 0 or
        !std.math.isFinite(input.grazer_respiration_fraction) or input.grazer_respiration_fraction < 0)
        return error.InvalidDisturbanceRemovalInput;
    const harvested = try sumHarvestProductComponents(input.harvested_by_component);
    const residue = try sumHarvestProductComponents(input.residue_by_component);
    const direct_litter = try sumHarvestProductComponents(input.direct_litter_by_component);
    var exported: canopy.ElementalMass = .{
        .carbon_g = harvested.carbon_g - residue.carbon_g,
        .nitrogen_g = harvested.nitrogen_g - residue.nitrogen_g,
        .phosphorus_g = harvested.phosphorus_g - residue.phosphorus_g,
    };
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(exported, field.name);
        if (!std.math.isFinite(value) or value < -1.0e-12)
            return error.DisturbanceRemovalResidueOverdraw;
        @field(exported, field.name) = @max(0, value);
    }
    var result: SourceOrderDisturbanceRemovalResult = .{
        .harvested_total = harvested,
        .residue_total = residue,
        .direct_litter_total = direct_litter,
        .plant_ecosystem_removal = .{},
        .grid_ecosystem_removal = .{},
        .reseed_storage_addition = .{},
        .net_biome_production_carbon_change_g_c_per_h = 0,
        .plant_total_respiration_change_g_c_per_h = 0,
        .plant_actual_respiration_change_g_c_per_h = 0,
        .ecosystem_respiration_change_g_c_per_h = 0,
        .autotrophic_respiration_change_g_c_per_h = 0,
    };
    const grazing = input.harvest_code == 4 or input.harvest_code == 6;
    if (!grazing) {
        if (input.terminate_and_reseed) {
            result.reseed_storage_addition = exported;
        } else {
            result.plant_ecosystem_removal = exported;
            result.grid_ecosystem_removal = exported;
            result.net_biome_production_carbon_change_g_c_per_h = -exported.carbon_g;
        }
        return result;
    }

    result.plant_ecosystem_removal = .{
        .carbon_g = input.grazer_growth_yield * exported.carbon_g,
        .nitrogen_g = exported.nitrogen_g,
        .phosphorus_g = exported.phosphorus_g,
    };
    result.grid_ecosystem_removal = result.plant_ecosystem_removal;
    const respired_carbon_g_c = input.grazer_respiration_fraction * exported.carbon_g;
    result.plant_total_respiration_change_g_c_per_h = -respired_carbon_g_c;
    result.plant_actual_respiration_change_g_c_per_h = -respired_carbon_g_c;
    result.ecosystem_respiration_change_g_c_per_h = -respired_carbon_g_c;
    result.autotrophic_respiration_change_g_c_per_h = -respired_carbon_g_c;
    return result;
}

pub const SourceOrderAbovegroundHarvestLitterInput = struct {
    harvest_code: i8,
    biomass_turnover_type: u8,
    root_profile_type: u8,
    residue_by_component: [harvest_product_component_count]canopy.ElementalMass,
    direct_litter_by_component: [harvest_product_component_count]canopy.ElementalMass,
    woody_composition: TillageElementComposition,
    nonstructural_kinetics: litter_partition.ElementFractions,
    foliar_kinetics: litter_partition.ElementFractions,
    nonfoliar_kinetics: litter_partition.ElementFractions,
    stalk_kinetics: litter_partition.ElementFractions,
    coarse_wood_kinetics: litter_partition.ElementFractions,
};

pub const SourceOrderAbovegroundHarvestLitterResult = struct {
    litter: canopy.SenescenceProducts,
    standing_dead_addition: canopy.SenescenceProducts,
};

/// Exact GROSUB 10796-10837 non-grazing above-ground harvest litterfall.
pub fn sourceOrderAbovegroundHarvestLitter(
    input: SourceOrderAbovegroundHarvestLitterInput,
) !SourceOrderAbovegroundHarvestLitterResult {
    if (input.harvest_code < 0 or input.harvest_code > 6 or
        input.harvest_code == 4 or input.harvest_code == 5 or input.harvest_code == 6)
        return error.InvalidAbovegroundHarvestLitterCode;
    _ = try sumHarvestProductComponents(input.residue_by_component);
    _ = try sumHarvestProductComponents(input.direct_litter_by_component);
    try validateTillageComposition(input.woody_composition);
    inline for (.{
        input.nonstructural_kinetics,
        input.foliar_kinetics,
        input.nonfoliar_kinetics,
        input.stalk_kinetics,
        input.coarse_wood_kinetics,
    }) |kinetics| kinetics.validate() catch return error.InvalidAbovegroundHarvestLitterInput;

    var result: SourceOrderAbovegroundHarvestLitterResult = .{
        .litter = .{},
        .standing_dead_addition = .{},
    };
    const residue = input.residue_by_component;
    const direct = input.direct_litter_by_component;
    const herbaceous = input.biomass_turnover_type == 0 or input.root_profile_type <= 1;
    for (0..litter_partition.kinetic_component_count) |kinetic| {
        result.litter.nonwoody_carbon_g[kinetic] =
            input.nonstructural_kinetics.carbon[kinetic] * (residue[0].carbon_g + direct[0].carbon_g) +
            input.foliar_kinetics.carbon[kinetic] * (residue[1].carbon_g + direct[1].carbon_g) +
            input.nonfoliar_kinetics.carbon[kinetic] * (residue[2].carbon_g + direct[2].carbon_g);
        result.litter.nonwoody_nitrogen_g[kinetic] =
            input.nonstructural_kinetics.nitrogen[kinetic] * (residue[0].nitrogen_g + direct[0].nitrogen_g) +
            input.foliar_kinetics.nitrogen[kinetic] * (residue[1].nitrogen_g + direct[1].nitrogen_g) +
            input.nonfoliar_kinetics.nitrogen[kinetic] * (residue[2].nitrogen_g + direct[2].nitrogen_g);
        result.litter.nonwoody_phosphorus_g[kinetic] =
            input.nonstructural_kinetics.phosphorus[kinetic] * (residue[0].phosphorus_g + direct[0].phosphorus_g) +
            input.foliar_kinetics.phosphorus[kinetic] * (residue[1].phosphorus_g + direct[1].phosphorus_g) +
            input.nonfoliar_kinetics.phosphorus[kinetic] * (residue[2].phosphorus_g + direct[2].phosphorus_g);
        if (herbaceous) {
            result.litter.nonwoody_carbon_g[kinetic] += input.stalk_kinetics.carbon[kinetic] *
                (residue[3].carbon_g + direct[3].carbon_g + residue[4].carbon_g + direct[4].carbon_g);
            result.litter.nonwoody_nitrogen_g[kinetic] += input.stalk_kinetics.nitrogen[kinetic] *
                (residue[3].nitrogen_g + direct[3].nitrogen_g + residue[4].nitrogen_g + direct[4].nitrogen_g);
            result.litter.nonwoody_phosphorus_g[kinetic] += input.stalk_kinetics.phosphorus[kinetic] *
                (residue[3].phosphorus_g + direct[3].phosphorus_g + residue[4].phosphorus_g + direct[4].phosphorus_g);
        } else {
            result.standing_dead_addition.woody_carbon_g[kinetic] =
                input.coarse_wood_kinetics.carbon[kinetic] * (direct[3].carbon_g + direct[4].carbon_g);
            result.standing_dead_addition.woody_nitrogen_g[kinetic] =
                input.coarse_wood_kinetics.nitrogen[kinetic] * (direct[3].nitrogen_g + direct[4].nitrogen_g);
            result.standing_dead_addition.woody_phosphorus_g[kinetic] =
                input.coarse_wood_kinetics.phosphorus[kinetic] * (direct[3].phosphorus_g + direct[4].phosphorus_g);
            result.litter.woody_carbon_g[kinetic] = input.coarse_wood_kinetics.carbon[kinetic] *
                (residue[3].carbon_g + residue[4].carbon_g) * input.woody_composition.carbon[0];
            result.litter.woody_nitrogen_g[kinetic] = input.coarse_wood_kinetics.nitrogen[kinetic] *
                (residue[3].nitrogen_g + residue[4].nitrogen_g) * input.woody_composition.nitrogen[0];
            result.litter.woody_phosphorus_g[kinetic] = input.coarse_wood_kinetics.phosphorus[kinetic] *
                (residue[3].phosphorus_g + residue[4].phosphorus_g) * input.woody_composition.phosphorus[0];
            result.litter.nonwoody_carbon_g[kinetic] += input.nonfoliar_kinetics.carbon[kinetic] *
                (residue[3].carbon_g + residue[4].carbon_g) * input.woody_composition.carbon[1];
            result.litter.nonwoody_nitrogen_g[kinetic] += input.nonfoliar_kinetics.nitrogen[kinetic] *
                (residue[3].nitrogen_g + residue[4].nitrogen_g) * input.woody_composition.nitrogen[1];
            result.litter.nonwoody_phosphorus_g[kinetic] += input.nonfoliar_kinetics.phosphorus[kinetic] *
                (residue[3].phosphorus_g + residue[4].phosphorus_g) * input.woody_composition.phosphorus[1];
        }
    }
    return result;
}

pub const SourceOrderGrazingLitterLedgerState = struct {
    hourly_litter: canopy.ElementalMass,
    cumulative_litter: canopy.ElementalMass,
    cumulative_aboveground_litter: canopy.ElementalMass,
    surface_litter_carbon_g_c: f64,
    accumulated_application: canopy.ElementalMass,
};

pub const SourceOrderGrazingLitterResult = struct {
    returned_mass: canopy.ElementalMass,
    state: SourceOrderGrazingLitterLedgerState,
    manure: grazing_manure.Products,
};

/// Exact GROSUB 10847-10906 grazing litter ledgers and manure deposition.
pub fn sourceOrderGrazingLitterLedgers(
    harvest_code: i8,
    residue_by_component: [harvest_product_component_count]canopy.ElementalMass,
    direct_litter_by_component: [harvest_product_component_count]canopy.ElementalMass,
    current: SourceOrderGrazingLitterLedgerState,
) !SourceOrderGrazingLitterResult {
    const kind: management.HarvestKind = switch (harvest_code) {
        4 => .animal_grazing,
        6 => .insect_grazing,
        else => return error.InvalidGrazingLitterCode,
    };
    const residue = try sumHarvestProductComponents(residue_by_component);
    const direct = try sumHarvestProductComponents(direct_litter_by_component);
    inline for (@typeInfo(SourceOrderGrazingLitterLedgerState).@"struct".fields) |field| {
        if (field.type == canopy.ElementalMass) {
            const mass = @field(current, field.name);
            inline for (.{ mass.carbon_g, mass.nitrogen_g, mass.phosphorus_g }) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.InvalidGrazingLitterLedger;
        } else {
            const value = @field(current, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidGrazingLitterLedger;
        }
    }
    const returned: canopy.ElementalMass = .{
        .carbon_g = residue.carbon_g + direct.carbon_g,
        .nitrogen_g = residue.nitrogen_g + direct.nitrogen_g,
        .phosphorus_g = residue.phosphorus_g + direct.phosphorus_g,
    };
    var next = current;
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        @field(next.hourly_litter, field.name) += @field(returned, field.name);
        @field(next.cumulative_litter, field.name) += @field(returned, field.name);
        @field(next.accumulated_application, field.name) += @field(returned, field.name);
    }
    next.cumulative_aboveground_litter.carbon_g += returned.carbon_g;
    // Preserve source lines 10854-10855: N/P use the just-updated total
    // cumulative ledger rather than their prior above-ground ledger.
    next.cumulative_aboveground_litter.nitrogen_g =
        next.cumulative_litter.nitrogen_g + returned.nitrogen_g;
    next.cumulative_aboveground_litter.phosphorus_g =
        next.cumulative_litter.phosphorus_g + returned.phosphorus_g;
    next.surface_litter_carbon_g_c += returned.carbon_g;
    inline for (@typeInfo(SourceOrderGrazingLitterLedgerState).@"struct".fields) |field| {
        if (field.type == canopy.ElementalMass) {
            const mass = @field(next, field.name);
            inline for (.{ mass.carbon_g, mass.nitrogen_g, mass.phosphorus_g }) |value|
                if (!std.math.isFinite(value) or value < 0)
                    return error.NonFiniteGrazingLitterLedger;
        } else if (!std.math.isFinite(@field(next, field.name)) or @field(next, field.name) < 0)
            return error.NonFiniteGrazingLitterLedger;
    }
    return .{
        .returned_mass = returned,
        .state = next,
        .manure = try grazing_manure.partition(kind, returned),
    };
}

/// Exact GROSUB 10914-10916 refresh after disturbance changes plant population.
pub fn sourceOrderPopulationScaledNumericalThresholds(
    plant_population_count: f64,
    horizontal_cell_area_m2: f64,
    mass_threshold_g_per_plant: f64,
    flux_threshold_g_per_plant_per_step: f64,
) !PopulationScaledNumericalThresholds {
    inline for (.{
        plant_population_count,
        horizontal_cell_area_m2,
        mass_threshold_g_per_plant,
        flux_threshold_g_per_plant_per_step,
    }) |value| {
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidPopulationScaledNumericalThresholdInput;
    }
    if (horizontal_cell_area_m2 == 0)
        return error.InvalidPopulationScaledNumericalThresholdInput;

    // Preserve the source evaluation order: multiply by population first,
    // then divide the mass threshold by horizontal area.
    const plant_mass_presence_g = mass_threshold_g_per_plant * plant_population_count;
    const result: PopulationScaledNumericalThresholds = .{
        .plant_mass_presence_g = plant_mass_presence_g,
        .plant_mass_density_g_m2 = plant_mass_presence_g / horizontal_cell_area_m2,
        .plant_flux_presence_g_per_step = flux_threshold_g_per_plant_per_step * plant_population_count,
    };
    inline for (@typeInfo(PopulationScaledNumericalThresholds).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFinitePopulationScaledNumericalThreshold;
    }
    return result;
}

pub const SourceOrderDeadBranchPhenologyState = struct {
    dead: bool,
    maturity_group_node_count: f64,
    initiated_node_count: f64,
    nodes_at_floral_initiation: f64,
    nodes_at_anthesis: f64,
    appeared_leaf_count: f64,
    leaves_at_floral_initiation: f64,
    current_leaf_ordinal: usize,
    current_growing_leaf_ordinal: usize,
    normalized_vegetative_node_change: f64,
    normalized_reproductive_node_change: f64,
    accumulated_leafout_h: f64,
    accumulated_leafoff_h: f64,
    lengthening_photoperiod_h: f64,
    shortening_photoperiod_h: f64,
    time_since_germination_h: f64,
    hours_without_grain_fill: f64,
    carbon_fixation_feedback: f64,
    carbon_fixation_feedback_previous: f64,
    leafout_initialization_enabled: bool,
    emergence_initialization_disabled: bool,
    leafoff_enabled: bool,
    remobilization_enabled: bool,
    hours_after_maturity_h: f64,
    new_branch_count: usize,
    stage_day_of_year: [10]u16,
};

pub const SourceOrderDeadBranchResetInput = struct {
    first_living_branch_emergence_day_of_year: u16,
    perennial_growth_habit: bool,
    current_day_of_year: u16,
    current_year: u32,
    harvest_day_of_year: u16,
    harvest_year: u32,
    hour_of_day: u8,
    local_solar_noon_h: f64,
    initial_maturity_group_node_count: f64,
    initial_node_count: f64,
};

/// Exact GROSUB 10955-10992 dead-branch phenology reset selector and loop.
pub fn sourceOrderResetDeadBranchPhenology(
    branches: []SourceOrderDeadBranchPhenologyState,
    input: SourceOrderDeadBranchResetInput,
) !bool {
    if (input.current_day_of_year == 0 or input.current_day_of_year > 366 or
        input.harvest_day_of_year == 0 or input.harvest_day_of_year > 366 or
        input.current_year == 0 or input.harvest_year == 0 or input.hour_of_day > 23 or
        !std.math.isFinite(input.local_solar_noon_h) or input.local_solar_noon_h < 0 or
        input.local_solar_noon_h >= 24 or
        !std.math.isFinite(input.initial_maturity_group_node_count) or
        !std.math.isFinite(input.initial_node_count))
        return error.InvalidDeadBranchPhenologyResetInput;
    for (branches) |branch| {
        inline for (@typeInfo(SourceOrderDeadBranchPhenologyState).@"struct".fields) |field| {
            if (field.type == f64 and !std.math.isFinite(@field(branch, field.name)))
                return error.InvalidDeadBranchPhenologyResetInput;
        }
    }

    const at_solar_noon = input.hour_of_day ==
        @as(u8, @intFromFloat(@floor(input.local_solar_noon_h)));
    const annual_harvest_reached =
        input.current_day_of_year >= input.harvest_day_of_year and
        input.current_year >= input.harvest_year and at_solar_noon;
    const selected = input.first_living_branch_emergence_day_of_year != 0 and
        (input.perennial_growth_habit or annual_harvest_reached);
    if (!selected) return false;

    for (branches) |*branch| {
        if (!branch.dead) continue;
        branch.maturity_group_node_count = input.initial_maturity_group_node_count;
        branch.initiated_node_count = input.initial_node_count;
        branch.nodes_at_floral_initiation = input.initial_node_count;
        branch.nodes_at_anthesis = 0;
        branch.appeared_leaf_count = 0;
        branch.leaves_at_floral_initiation = 0;
        branch.current_leaf_ordinal = 1;
        branch.current_growing_leaf_ordinal = 1;
        branch.normalized_vegetative_node_change = 0;
        branch.normalized_reproductive_node_change = 0;
        branch.accumulated_leafout_h = 0;
        branch.accumulated_leafoff_h = 0;
        branch.lengthening_photoperiod_h = 0;
        branch.shortening_photoperiod_h = 0;
        branch.time_since_germination_h = 0;
        branch.hours_without_grain_fill = 0;
        branch.carbon_fixation_feedback = 1;
        branch.carbon_fixation_feedback_previous = 1;
        branch.leafout_initialization_enabled = true;
        branch.emergence_initialization_disabled = true;
        branch.leafoff_enabled = true;
        branch.remobilization_enabled = true;
        branch.hours_after_maturity_h = 0;
        branch.new_branch_count = 0;
        branch.stage_day_of_year = @splat(0);
    }
    return true;
}

pub const SourceOrderDeadBranchLitterPools = struct {
    bacterial_nonstructural: canopy.ElementalMass,
    bacterial_structural: canopy.ElementalMass,
    leaf: canopy.ElementalMass,
    sheath: canopy.ElementalMass,
    husk: canopy.ElementalMass,
    ear: canopy.ElementalMass,
    grain: canopy.ElementalMass,
    stalk: canopy.ElementalMass,
};

pub const SourceOrderDeadBranchLitterInput = struct {
    annual_growth_habit: bool,
    deciduous_phenology: bool,
    pools: SourceOrderDeadBranchLitterPools,
    leaf_woody_fraction: TillageElementComposition,
    sheath_woody_fraction: TillageElementComposition,
    nonstructural_kinetics: litter_partition.ElementFractions,
    foliar_kinetics: litter_partition.ElementFractions,
    nonfoliar_kinetics: litter_partition.ElementFractions,
    stalk_kinetics: litter_partition.ElementFractions,
    coarse_wood_kinetics: litter_partition.ElementFractions,
};

pub const SourceOrderDeadBranchLitterResult = struct {
    nonwoody_litter: [4]canopy.ElementalMass,
    woody_litter: [4]canopy.ElementalMass,
    standing_dead_stalk: [4]canopy.ElementalMass,
    seasonal_storage_addition: canopy.ElementalMass,
};

/// Exact GROSUB 11031-11083 four-component dead-branch litter partition.
pub fn sourceOrderDeadBranchLitterfall(
    input: SourceOrderDeadBranchLitterInput,
) !SourceOrderDeadBranchLitterResult {
    inline for (@typeInfo(SourceOrderDeadBranchLitterPools).@"struct".fields) |field| {
        const mass = @field(input.pools, field.name);
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element| {
            const value = @field(mass, element.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidDeadBranchLitterfallInput;
        }
    }
    inline for (.{ input.leaf_woody_fraction, input.sheath_woody_fraction }) |fractions| {
        inline for (@typeInfo(TillageElementComposition).@"struct".fields) |field| {
            for (@field(fractions, field.name)) |value| {
                if (!std.math.isFinite(value) or value < 0 or value > 1)
                    return error.InvalidDeadBranchLitterfallInput;
            }
        }
    }
    inline for (.{
        input.nonstructural_kinetics,
        input.foliar_kinetics,
        input.nonfoliar_kinetics,
        input.stalk_kinetics,
        input.coarse_wood_kinetics,
    }) |kinetics| {
        inline for (@typeInfo(litter_partition.ElementFractions).@"struct".fields) |field| {
            for (@field(kinetics, field.name)) |value| {
                if (!std.math.isFinite(value) or value < 0 or value > 1)
                    return error.InvalidDeadBranchLitterfallInput;
            }
        }
    }

    var result: SourceOrderDeadBranchLitterResult = .{
        .nonwoody_litter = @splat(.{}),
        .woody_litter = @splat(.{}),
        .standing_dead_stalk = @splat(.{}),
        .seasonal_storage_addition = .{},
    };
    const winter_annual = input.annual_growth_habit and input.deciduous_phenology;
    inline for (0..4) |component| {
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields, 0..) |element, element_index| {
            const bacterial_nonstructural = @field(input.pools.bacterial_nonstructural, element.name);
            const bacterial_structural = @field(input.pools.bacterial_structural, element.name);
            const leaf = @field(input.pools.leaf, element.name);
            const sheath = @field(input.pools.sheath, element.name);
            const husk = @field(input.pools.husk, element.name);
            const ear = @field(input.pools.ear, element.name);
            const grain = @field(input.pools.grain, element.name);
            const stalk = @field(input.pools.stalk, element.name);
            const fraction_field =
                @typeInfo(litter_partition.ElementFractions).@"struct".fields[element_index].name;
            const nonstructural_fraction = @field(input.nonstructural_kinetics, fraction_field)[component];
            const foliar_fraction = @field(input.foliar_kinetics, fraction_field)[component];
            const nonfoliar_fraction = @field(input.nonfoliar_kinetics, fraction_field)[component];
            const stalk_fraction = @field(input.stalk_kinetics, fraction_field)[component];
            const coarse_wood_fraction = @field(input.coarse_wood_kinetics, fraction_field)[component];
            const leaf_wood = @field(input.leaf_woody_fraction, fraction_field);
            const sheath_wood = @field(input.sheath_woody_fraction, fraction_field);

            @field(result.nonwoody_litter[component], element.name) +=
                nonstructural_fraction * bacterial_nonstructural;
            @field(result.nonwoody_litter[component], element.name) +=
                foliar_fraction * (leaf * leaf_wood[1] + bacterial_structural);
            @field(result.nonwoody_litter[component], element.name) +=
                nonfoliar_fraction * (sheath * sheath_wood[1] + husk + ear);
            @field(result.woody_litter[component], element.name) +=
                coarse_wood_fraction * (leaf * leaf_wood[0] + sheath * sheath_wood[0]);
            if (winter_annual) {
                @field(result.seasonal_storage_addition, element.name) +=
                    nonfoliar_fraction * grain;
            } else {
                @field(result.nonwoody_litter[component], element.name) +=
                    nonfoliar_fraction * grain;
            }
            @field(result.standing_dead_stalk[component], element.name) +=
                stalk_fraction * stalk;
        }
    }
    inline for (@typeInfo(SourceOrderDeadBranchLitterResult).@"struct".fields) |field| {
        const value = @field(result, field.name);
        if (field.type == canopy.ElementalMass) {
            inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element|
                if (!std.math.isFinite(@field(value, element.name)))
                    return error.NonFiniteDeadBranchLitterfall;
        } else for (value) |mass| {
            inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element|
                if (!std.math.isFinite(@field(mass, element.name)))
                    return error.NonFiniteDeadBranchLitterfall;
        }
    }
    return result;
}

pub const SourceOrderDeadBranchStorageRecoveryInput = struct {
    current_seasonal_storage: canopy.ElementalMass,
    branch_mobile: canopy.ElementalMass,
    c4_intermediate_carbon_g_c: f64,
    stalk_reserve: canopy.ElementalMass,
};

pub const SourceOrderDeadBranchStorageRecoveryResult = struct {
    seasonal_storage: canopy.ElementalMass,
    recovered: canopy.ElementalMass,
};

/// Exact GROSUB 11100-11106 dead-branch mobile and reserve recovery.
pub fn sourceOrderDeadBranchStorageRecovery(
    input: SourceOrderDeadBranchStorageRecoveryInput,
) !SourceOrderDeadBranchStorageRecoveryResult {
    inline for (.{ input.current_seasonal_storage, input.branch_mobile, input.stalk_reserve }) |mass| {
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
            const value = @field(mass, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidDeadBranchStorageRecoveryInput;
        }
    }
    if (!std.math.isFinite(input.c4_intermediate_carbon_g_c) or
        input.c4_intermediate_carbon_g_c < 0)
        return error.InvalidDeadBranchStorageRecoveryInput;

    // Preserve the six source assignments, including the two sequential
    // additions to storage carbon.
    var seasonal_storage = input.current_seasonal_storage;
    seasonal_storage.carbon_g =
        seasonal_storage.carbon_g +
        input.branch_mobile.carbon_g +
        input.c4_intermediate_carbon_g_c;
    seasonal_storage.nitrogen_g =
        seasonal_storage.nitrogen_g + input.branch_mobile.nitrogen_g;
    seasonal_storage.phosphorus_g =
        seasonal_storage.phosphorus_g + input.branch_mobile.phosphorus_g;
    seasonal_storage.carbon_g =
        seasonal_storage.carbon_g + input.stalk_reserve.carbon_g;
    seasonal_storage.nitrogen_g =
        seasonal_storage.nitrogen_g + input.stalk_reserve.nitrogen_g;
    seasonal_storage.phosphorus_g =
        seasonal_storage.phosphorus_g + input.stalk_reserve.phosphorus_g;

    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(seasonal_storage, field.name)))
            return error.NonFiniteDeadBranchStorageRecovery;
    }
    return .{
        .seasonal_storage = seasonal_storage,
        .recovered = .{
            .carbon_g = input.branch_mobile.carbon_g +
                input.c4_intermediate_carbon_g_c +
                input.stalk_reserve.carbon_g,
            .nitrogen_g = input.branch_mobile.nitrogen_g +
                input.stalk_reserve.nitrogen_g,
            .phosphorus_g = input.branch_mobile.phosphorus_g +
                input.stalk_reserve.phosphorus_g,
        },
    };
}

pub const SourceOrderDeadBranchScalarResetState = struct {
    host_mobile: canopy.ElementalMass,
    c4_intermediate_carbon_g_c: f64,
    symbiont_mobile: canopy.ElementalMass,
    shoot_total: canopy.ElementalMass,
    leaf: canopy.ElementalMass,
    symbiont_structural: canopy.ElementalMass,
    sheath: canopy.ElementalMass,
    stalk: canopy.ElementalMass,
    vascular_stalk_carbon_g_c: f64,
    stalk_reserve: canopy.ElementalMass,
    husk: canopy.ElementalMass,
    ear: canopy.ElementalMass,
    grain: canopy.ElementalMass,
    live_symbiont_carbon_g_c: f64,
    potential_seed_sites: f64,
    seed_count: f64,
    individual_seed_carbon_g_c: f64,
    leaf_area_m2: f64,
    stale_stalk_total: canopy.ElementalMass,
};

pub const SourceOrderDeadBranchNodeResetState = struct {
    bundle_sheath_mobile_carbon_g_c: f64,
    mesophyll_mobile_carbon_g_c: f64,
    bundle_sheath_co2_carbon_g_c: f64,
    bundle_sheath_bicarbonate_carbon_g_c: f64,
    leaf_area_m2: f64,
    node_height_m: f64,
    node_height_previous_m: f64,
    sheath_height_m: f64,
    leaf: canopy.ElementalMass,
    leaf_protein_g: f64,
    sheath: canopy.ElementalMass,
    sheath_protein_g: f64,
    stalk: canopy.ElementalMass,
};

pub const SourceOrderDeadBranchLayerResetState = struct {
    leaf_area_m2: f64,
    leaf: canopy.ElementalMass,
    projected_leaf_surface_m2: [4]f64,
};

pub const SourceOrderDeadBranchCanopyReset = struct {
    scalar: *SourceOrderDeadBranchScalarResetState,
    nodes: []SourceOrderDeadBranchNodeResetState,
    node_layers: []SourceOrderDeadBranchLayerResetState,
    canopy_leaf_area_m2_by_layer: []f64,
    canopy_leaf_carbon_g_c_by_layer: []f64,
    branch_stalk_area_m2_by_layer: []f64,
    branch_projected_stalk_surface_m2: [][4]f64,
};

/// Exact GROSUB 11137-11220 dead-branch scalar, node, and canopy reset.
pub fn sourceOrderResetDeadBranchCanopy(state: SourceOrderDeadBranchCanopyReset) !void {
    const layer_count = state.canopy_leaf_area_m2_by_layer.len;
    const expected_node_layers = std.math.mul(usize, state.nodes.len, layer_count) catch
        return error.InvalidDeadBranchCanopyResetDimensions;
    if (layer_count == 0 or
        state.canopy_leaf_carbon_g_c_by_layer.len != layer_count or
        state.branch_stalk_area_m2_by_layer.len != layer_count or
        state.branch_projected_stalk_surface_m2.len != layer_count or
        state.node_layers.len != expected_node_layers)
        return error.InvalidDeadBranchCanopyResetDimensions;
    inline for (@typeInfo(SourceOrderDeadBranchScalarResetState).@"struct".fields) |field| {
        const value = @field(state.scalar.*, field.name);
        if (field.type == f64) {
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidDeadBranchCanopyResetInput;
        } else inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element| {
            const mass = @field(value, element.name);
            if (!std.math.isFinite(mass) or mass < 0)
                return error.InvalidDeadBranchCanopyResetInput;
        }
    }
    for (state.nodes) |node| inline for (@typeInfo(SourceOrderDeadBranchNodeResetState).@"struct".fields) |field| {
        const value = @field(node, field.name);
        if (field.type == f64) {
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidDeadBranchCanopyResetInput;
        } else inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element| {
            const mass = @field(value, element.name);
            if (!std.math.isFinite(mass) or mass < 0)
                return error.InvalidDeadBranchCanopyResetInput;
        }
    };
    for (0..layer_count) |layer| {
        var branch_area_m2: f64 = 0;
        var branch_carbon_g_c: f64 = 0;
        for (0..state.nodes.len) |node| {
            const contribution = state.node_layers[node * layer_count + layer];
            if (!std.math.isFinite(contribution.leaf_area_m2) or contribution.leaf_area_m2 < 0)
                return error.InvalidDeadBranchCanopyResetInput;
            inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element| {
                const mass = @field(contribution.leaf, element.name);
                if (!std.math.isFinite(mass) or mass < 0)
                    return error.InvalidDeadBranchCanopyResetInput;
            }
            for (contribution.projected_leaf_surface_m2) |surface|
                if (!std.math.isFinite(surface) or surface < 0)
                    return error.InvalidDeadBranchCanopyResetInput;
            branch_area_m2 += contribution.leaf_area_m2;
            branch_carbon_g_c += contribution.leaf.carbon_g;
        }
        if (!std.math.isFinite(state.canopy_leaf_area_m2_by_layer[layer]) or
            !std.math.isFinite(state.canopy_leaf_carbon_g_c_by_layer[layer]) or
            state.canopy_leaf_area_m2_by_layer[layer] < branch_area_m2 or
            state.canopy_leaf_carbon_g_c_by_layer[layer] < branch_carbon_g_c or
            !std.math.isFinite(state.branch_stalk_area_m2_by_layer[layer]) or
            state.branch_stalk_area_m2_by_layer[layer] < 0)
            return error.InvalidDeadBranchCanopyResetInput;
        for (state.branch_projected_stalk_surface_m2[layer]) |surface|
            if (!std.math.isFinite(surface) or surface < 0)
                return error.InvalidDeadBranchCanopyResetInput;
    }

    state.scalar.host_mobile = .{};
    state.scalar.c4_intermediate_carbon_g_c = 0;
    state.scalar.symbiont_mobile = .{};
    state.scalar.shoot_total = .{};
    state.scalar.leaf = .{};
    state.scalar.symbiont_structural = .{};
    state.scalar.sheath = .{};
    state.scalar.stalk = .{};
    state.scalar.vascular_stalk_carbon_g_c = 0;
    state.scalar.stalk_reserve = .{};
    state.scalar.husk = .{};
    state.scalar.ear = .{};
    state.scalar.grain = .{};
    state.scalar.live_symbiont_carbon_g_c = 0;
    state.scalar.potential_seed_sites = 0;
    state.scalar.seed_count = 0;
    state.scalar.individual_seed_carbon_g_c = 0;
    state.scalar.leaf_area_m2 = 0;
    state.scalar.stale_stalk_total = .{};
    for (state.nodes, 0..) |*node, node_index| {
        if (node_index != 0) {
            node.bundle_sheath_mobile_carbon_g_c = 0;
            node.mesophyll_mobile_carbon_g_c = 0;
            node.bundle_sheath_co2_carbon_g_c = 0;
            node.bundle_sheath_bicarbonate_carbon_g_c = 0;
        }
        node.leaf_area_m2 = 0;
        node.node_height_m = 0;
        node.node_height_previous_m = 0;
        node.sheath_height_m = 0;
        node.leaf = .{};
        node.leaf_protein_g = 0;
        node.sheath = .{};
        node.sheath_protein_g = 0;
        node.stalk = .{};
        for (0..layer_count) |layer| {
            const offset = node_index * layer_count + layer;
            const contribution = state.node_layers[offset];
            state.canopy_leaf_area_m2_by_layer[layer] -= contribution.leaf_area_m2;
            state.canopy_leaf_carbon_g_c_by_layer[layer] -= contribution.leaf.carbon_g;
            state.node_layers[offset].leaf_area_m2 = 0;
            state.node_layers[offset].leaf = .{};
            if (node_index != 0)
                state.node_layers[offset].projected_leaf_surface_m2 = @splat(0);
        }
    }
    @memset(state.branch_stalk_area_m2_by_layer, 0);
    @memset(state.branch_projected_stalk_surface_m2, @splat(0));
}

pub const SourceOrderWholePlantTerminationState = struct {
    shoot_alive: bool,
    root_alive: bool,
    total_node_count: usize,
    hours_below_leaf_turgor_threshold_h: f64,
    main_stalk_diameter_m: f64,
    branch_count: usize,
    living_population_per_m2: f64,
    living_population_count: f64,
    hypocotyl_height_m: f64,
};

pub const SourceOrderWholePlantTerminationResult = struct {
    state: SourceOrderWholePlantTerminationState,
    dead_branch_count: usize,
    all_branches_dead: bool,
};

/// Exact GROSUB 11221-11240 all-dead branch count and plant transition.
pub fn sourceOrderWholePlantTermination(
    branch_is_dead: []const bool,
    winter_annual: bool,
    current: SourceOrderWholePlantTerminationState,
) !SourceOrderWholePlantTerminationResult {
    if (current.branch_count != branch_is_dead.len or
        !std.math.isFinite(current.hours_below_leaf_turgor_threshold_h) or
        current.hours_below_leaf_turgor_threshold_h < 0 or
        !std.math.isFinite(current.main_stalk_diameter_m) or
        current.main_stalk_diameter_m < 0 or
        !std.math.isFinite(current.living_population_per_m2) or
        current.living_population_per_m2 < 0 or
        !std.math.isFinite(current.living_population_count) or
        current.living_population_count < 0 or
        !std.math.isFinite(current.hypocotyl_height_m) or
        current.hypocotyl_height_m < 0)
        return error.InvalidWholePlantTerminationInput;

    var dead_branch_count: usize = 0;
    for (branch_is_dead) |is_dead| {
        if (is_dead) dead_branch_count += 1;
    }
    if (dead_branch_count != current.branch_count) return .{
        .state = current,
        .dead_branch_count = dead_branch_count,
        .all_branches_dead = false,
    };

    var state = current;
    state.shoot_alive = false;
    state.root_alive = false;
    state.total_node_count = 0;
    state.hours_below_leaf_turgor_threshold_h = 0;
    state.main_stalk_diameter_m = 0;
    if (winter_annual) {
        state.branch_count = 1;
    } else {
        state.branch_count = 0;
        state.living_population_per_m2 = 0;
        state.living_population_count = 0;
    }
    state.hypocotyl_height_m = 0;
    return .{
        .state = state,
        .dead_branch_count = dead_branch_count,
        .all_branches_dead = true,
    };
}

pub const SourceOrderDeadRootAxisPools = struct {
    primary: canopy.ElementalMass,
    secondary: canopy.ElementalMass,
};

pub const SourceOrderDeadRootLitterInput = struct {
    roots_dead: bool,
    root_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    mobile_by_domain_layer: []const canopy.ElementalMass,
    structural_by_domain_layer_axis: []const SourceOrderDeadRootAxisPools,
    root_woody_fraction: TillageElementComposition,
    mobile_kinetics: litter_partition.ElementFractions,
    fine_root_kinetics: litter_partition.ElementFractions,
    coarse_root_kinetics: litter_partition.ElementFractions,
};

pub const SourceOrderDeadRootLayerLitter = struct {
    woody: [4]canopy.ElementalMass,
    nonwoody: [4]canopy.ElementalMass,
};

/// Exact GROSUB 11265-11288 dead-root C/N/P litterfall.
/// Caller owns the returned runtime-layer slice.
pub fn sourceOrderDeadRootLitterfall(
    allocator: std.mem.Allocator,
    input: SourceOrderDeadRootLitterInput,
) ![]SourceOrderDeadRootLayerLitter {
    if (input.root_domain_count == 0 or input.soil_layer_count == 0 or
        input.root_axis_count == 0)
        return error.InvalidDeadRootLitterfallDimensions;
    const domain_layers = std.math.mul(
        usize,
        input.root_domain_count,
        input.soil_layer_count,
    ) catch return error.InvalidDeadRootLitterfallDimensions;
    const structural_count = std.math.mul(
        usize,
        domain_layers,
        input.root_axis_count,
    ) catch return error.InvalidDeadRootLitterfallDimensions;
    if (input.mobile_by_domain_layer.len != domain_layers or
        input.structural_by_domain_layer_axis.len != structural_count)
        return error.InvalidDeadRootLitterfallDimensions;
    inline for (.{ input.mobile_kinetics, input.fine_root_kinetics, input.coarse_root_kinetics }) |kinetics| {
        inline for (@typeInfo(litter_partition.ElementFractions).@"struct".fields) |field| {
            for (@field(kinetics, field.name)) |fraction|
                if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
                    return error.InvalidDeadRootLitterfallInput;
        }
    }
    inline for (@typeInfo(TillageElementComposition).@"struct".fields) |field| {
        for (@field(input.root_woody_fraction, field.name)) |fraction|
            if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
                return error.InvalidDeadRootLitterfallInput;
    }
    for (input.mobile_by_domain_layer) |mass| try validateSourceOrderRootMass(mass);
    for (input.structural_by_domain_layer_axis) |axis| {
        try validateSourceOrderRootMass(axis.primary);
        try validateSourceOrderRootMass(axis.secondary);
    }

    const result = try allocator.alloc(SourceOrderDeadRootLayerLitter, input.soil_layer_count);
    errdefer allocator.free(result);
    @memset(result, .{ .woody = @splat(.{}), .nonwoody = @splat(.{}) });
    if (!input.roots_dead) return result;

    for (0..input.root_domain_count) |domain| {
        for (0..input.soil_layer_count) |layer| {
            const domain_layer = domain * input.soil_layer_count + layer;
            inline for (0..4) |component| {
                inline for (@typeInfo(canopy.ElementalMass).@"struct".fields, 0..) |element, element_index| {
                    const fraction_field =
                        @typeInfo(litter_partition.ElementFractions).@"struct".fields[element_index].name;
                    @field(result[layer].nonwoody[component], element.name) +=
                        @field(input.mobile_kinetics, fraction_field)[component] *
                        @field(input.mobile_by_domain_layer[domain_layer], element.name);
                }
                for (0..input.root_axis_count) |axis| {
                    const pools = input.structural_by_domain_layer_axis[
                        domain_layer * input.root_axis_count + axis
                    ];
                    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields, 0..) |element, element_index| {
                        const fraction_field =
                            @typeInfo(litter_partition.ElementFractions).@"struct".fields[element_index].name;
                        const structural_mass =
                            @field(pools.primary, element.name) +
                            @field(pools.secondary, element.name);
                        const wood = @field(input.root_woody_fraction, fraction_field);
                        @field(result[layer].woody[component], element.name) +=
                            @field(input.coarse_root_kinetics, fraction_field)[component] *
                            structural_mass * wood[0];
                        @field(result[layer].nonwoody[component], element.name) +=
                            @field(input.fine_root_kinetics, fraction_field)[component] *
                            structural_mass * wood[1];
                    }
                }
            }
        }
    }
    for (result) |layer| inline for (.{ layer.woody, layer.nonwoody }) |position| {
        for (position) |mass| inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element|
            if (!std.math.isFinite(@field(mass, element.name)))
                return error.NonFiniteDeadRootLitterfall;
    };
    return result;
}

pub const SourceOrderRootGasInventory = struct {
    carbon_dioxide_carbon_g_c: f64,
    oxygen_g_o: f64,
    methane_carbon_g_c: f64,
    nitrous_oxide_nitrogen_g_n: f64,
    ammonia_nitrogen_g_n: f64,
    hydrogen_g_h: f64,
};

pub const SourceOrderRootGasPhases = struct {
    gaseous: SourceOrderRootGasInventory,
    aqueous: SourceOrderRootGasInventory,
};

/// Exact GROSUB 11292-11315 dead-root gas release and phase clearing.
pub fn sourceOrderReleaseDeadRootGases(
    roots_dead: bool,
    domain_layer_phases: []SourceOrderRootGasPhases,
    current_disturbance_loss: SourceOrderRootGasInventory,
) !SourceOrderRootGasInventory {
    inline for (@typeInfo(SourceOrderRootGasInventory).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(current_disturbance_loss, field.name)))
            return error.InvalidDeadRootGasReleaseInput;
    }
    for (domain_layer_phases) |phases| {
        inline for (.{ phases.gaseous, phases.aqueous }) |phase| {
            inline for (@typeInfo(SourceOrderRootGasInventory).@"struct".fields) |field| {
                const value = @field(phase, field.name);
                if (!std.math.isFinite(value) or value < 0)
                    return error.InvalidDeadRootGasReleaseInput;
            }
        }
    }
    if (!roots_dead) return current_disturbance_loss;

    var disturbance_loss = current_disturbance_loss;
    // The caller supplies domain-major/layer-minor storage, matching the
    // enclosing source loops. Preserve gaseous-then-aqueous subtraction.
    for (domain_layer_phases) |phases| {
        inline for (@typeInfo(SourceOrderRootGasInventory).@"struct".fields) |field| {
            @field(disturbance_loss, field.name) =
                @field(disturbance_loss, field.name) -
                @field(phases.gaseous, field.name);
            @field(disturbance_loss, field.name) =
                @field(disturbance_loss, field.name) -
                @field(phases.aqueous, field.name);
            if (!std.math.isFinite(@field(disturbance_loss, field.name)))
                return error.NonFiniteDeadRootGasRelease;
        }
    }
    for (domain_layer_phases) |*phases| {
        phases.gaseous = std.mem.zeroes(SourceOrderRootGasInventory);
        phases.aqueous = std.mem.zeroes(SourceOrderRootGasInventory);
    }
    return disturbance_loss;
}

pub const SourceOrderDeadRootAxisLayerState = struct {
    primary: canopy.ElementalMass,
    secondary: canopy.ElementalMass,
    primary_length_m: f64,
    secondary_length_m: f64,
    secondary_axis_count: f64,
};

pub const SourceOrderDeadRootDomainAxisState = struct {
    primary_total: canopy.ElementalMass,
};

pub const SourceOrderDeadRootDomainLayerState = struct {
    mobile: canopy.ElementalMass,
    active_root_carbon_g_c: f64,
    actual_root_carbon_g_c: f64,
    root_protein_g: f64,
    primary_axis_count: f64,
    total_axis_count: f64,
    root_length_per_plant_m: f64,
    root_length_density_m_m3: f64,
    gaseous_volume_m3: f64,
    aqueous_volume_m3: f64,
    root_surface_area_per_plant_m2: f64,
    primary_radius_m: f64,
    secondary_radius_m: f64,
    average_secondary_root_length_m: f64,
};

pub const SourceOrderDeadRootResetState = struct {
    root_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    axis_layer: []SourceOrderDeadRootAxisLayerState,
    domain_axis: []SourceOrderDeadRootDomainAxisState,
    domain_layer: []SourceOrderDeadRootDomainLayerState,
    initial_primary_radius_m_by_domain: []const f64,
    initial_secondary_radius_m_by_domain: []const f64,
    initial_average_secondary_root_length_m: f64,
};

/// Exact GROSUB 11336-11365 root-axis and domain-layer reset.
pub fn sourceOrderResetDeadRootState(
    roots_dead: bool,
    state: SourceOrderDeadRootResetState,
) !void {
    if (state.root_domain_count == 0 or state.soil_layer_count == 0 or
        state.root_axis_count == 0)
        return error.InvalidDeadRootResetDimensions;
    const domain_layers = std.math.mul(
        usize,
        state.root_domain_count,
        state.soil_layer_count,
    ) catch return error.InvalidDeadRootResetDimensions;
    const domain_axes = std.math.mul(
        usize,
        state.root_domain_count,
        state.root_axis_count,
    ) catch return error.InvalidDeadRootResetDimensions;
    const axis_layers = std.math.mul(
        usize,
        domain_layers,
        state.root_axis_count,
    ) catch return error.InvalidDeadRootResetDimensions;
    if (state.axis_layer.len != axis_layers or
        state.domain_axis.len != domain_axes or
        state.domain_layer.len != domain_layers or
        state.initial_primary_radius_m_by_domain.len != state.root_domain_count or
        state.initial_secondary_radius_m_by_domain.len != state.root_domain_count or
        !std.math.isFinite(state.initial_average_secondary_root_length_m) or
        state.initial_average_secondary_root_length_m < 0)
        return error.InvalidDeadRootResetDimensions;

    for (state.initial_primary_radius_m_by_domain, state.initial_secondary_radius_m_by_domain) |primary, secondary| {
        if (!std.math.isFinite(primary) or primary < 0 or
            !std.math.isFinite(secondary) or secondary < 0)
            return error.InvalidDeadRootResetInput;
    }
    for (state.axis_layer) |axis| {
        try validateDeadRootResetMass(axis.primary);
        try validateDeadRootResetMass(axis.secondary);
        inline for (.{ axis.primary_length_m, axis.secondary_length_m, axis.secondary_axis_count }) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidDeadRootResetInput;
    }
    for (state.domain_axis) |axis| try validateDeadRootResetMass(axis.primary_total);
    for (state.domain_layer) |layer| {
        try validateDeadRootResetMass(layer.mobile);
        inline for (@typeInfo(SourceOrderDeadRootDomainLayerState).@"struct".fields) |field| {
            if (field.type != f64) continue;
            const value = @field(layer, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidDeadRootResetInput;
        }
    }
    if (!roots_dead) return;

    for (0..state.root_domain_count) |domain| {
        for (0..state.soil_layer_count) |layer| {
            const domain_layer = domain * state.soil_layer_count + layer;
            for (0..state.root_axis_count) |axis| {
                const axis_layer = domain_layer * state.root_axis_count + axis;
                state.axis_layer[axis_layer].primary = .{};
                state.axis_layer[axis_layer].secondary = .{};
                state.domain_axis[domain * state.root_axis_count + axis].primary_total = .{};
                state.axis_layer[axis_layer].primary_length_m = 0;
                state.axis_layer[axis_layer].secondary_length_m = 0;
                state.axis_layer[axis_layer].secondary_axis_count = 0;
            }
            state.domain_layer[domain_layer].mobile = .{};
            state.domain_layer[domain_layer].active_root_carbon_g_c = 0;
            state.domain_layer[domain_layer].actual_root_carbon_g_c = 0;
            state.domain_layer[domain_layer].root_protein_g = 0;
            state.domain_layer[domain_layer].primary_axis_count = 0;
            state.domain_layer[domain_layer].total_axis_count = 0;
            state.domain_layer[domain_layer].root_length_per_plant_m = 0;
            state.domain_layer[domain_layer].root_length_density_m_m3 = 0;
            state.domain_layer[domain_layer].gaseous_volume_m3 = 0;
            state.domain_layer[domain_layer].aqueous_volume_m3 = 0;
            state.domain_layer[domain_layer].primary_radius_m =
                state.initial_primary_radius_m_by_domain[domain];
            state.domain_layer[domain_layer].secondary_radius_m =
                state.initial_secondary_radius_m_by_domain[domain];
            state.domain_layer[domain_layer].root_surface_area_per_plant_m2 = 0;
            state.domain_layer[domain_layer].average_secondary_root_length_m =
                state.initial_average_secondary_root_length_m;
        }
    }
}

pub const SourceOrderDeadNoduleLayerPools = struct {
    structural: canopy.ElementalMass,
    mobile: canopy.ElementalMass,
};

pub const SourceOrderDeadNoduleLitterInput = struct {
    roots_dead: bool,
    nitrogen_fixation_enabled: bool,
    root_domain_count: usize,
    layer_pools: []SourceOrderDeadNoduleLayerPools,
    structural_kinetics: litter_partition.ElementFractions,
    mobile_kinetics: litter_partition.ElementFractions,
};

/// Exact GROSUB 11378-11393 dead-nodule litterfall and clearing.
/// Caller owns the returned runtime-layer slice.
pub fn sourceOrderDeadNoduleLitterfall(
    allocator: std.mem.Allocator,
    input: SourceOrderDeadNoduleLitterInput,
) ![][4]canopy.ElementalMass {
    if (input.root_domain_count == 0 or input.layer_pools.len == 0)
        return error.InvalidDeadNoduleLitterfallDimensions;
    inline for (.{ input.structural_kinetics, input.mobile_kinetics }) |kinetics| {
        inline for (@typeInfo(litter_partition.ElementFractions).@"struct".fields) |field| {
            for (@field(kinetics, field.name)) |fraction|
                if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
                    return error.InvalidDeadNoduleLitterfallInput;
        }
    }
    for (input.layer_pools) |pools| {
        try validateDeadNoduleMass(pools.structural);
        try validateDeadNoduleMass(pools.mobile);
    }

    const litter = try allocator.alloc([4]canopy.ElementalMass, input.layer_pools.len);
    errdefer allocator.free(litter);
    @memset(litter, @splat(.{}));
    if (!input.roots_dead or !input.nitrogen_fixation_enabled) return litter;

    // Preserve the enclosing domain/layer traversal and the source N == 1
    // condition. Nodule pools have layer ownership and are processed once.
    for (0..input.root_domain_count) |domain| {
        if (domain != 0) continue;
        for (input.layer_pools, 0..) |pools, layer| {
            inline for (0..4) |component| {
                inline for (@typeInfo(canopy.ElementalMass).@"struct".fields, 0..) |element, element_index| {
                    const fraction_field =
                        @typeInfo(litter_partition.ElementFractions).@"struct".fields[element_index].name;
                    @field(litter[layer][component], element.name) +=
                        @field(input.structural_kinetics, fraction_field)[component] *
                        @field(pools.structural, element.name) +
                        @field(input.mobile_kinetics, fraction_field)[component] *
                            @field(pools.mobile, element.name);
                    if (!std.math.isFinite(@field(litter[layer][component], element.name)))
                        return error.NonFiniteDeadNoduleLitterfall;
                }
            }
        }
    }
    for (input.layer_pools) |*pools| {
        pools.structural = .{};
        pools.mobile = .{};
    }
    return litter;
}

pub const SourceOrderDeadRootDepthResetState = struct {
    root_domain_count: usize,
    root_axis_count: usize,
    deepest_layer_by_axis: []usize,
    primary_depth_from_surface_m_by_axis_domain: []f64,
    primary_total_by_axis_domain: []canopy.ElementalMass,
    deepest_active_root_layer: *usize,
    active_root_axis_count: *usize,
};

/// Exact GROSUB 11403-11414 dead-root depth and axis-count reset.
pub fn sourceOrderResetDeadRootDepth(
    roots_dead: bool,
    planting_layer_index: usize,
    seed_depth_m: f64,
    state: SourceOrderDeadRootDepthResetState,
) !void {
    if (state.root_domain_count == 0 or
        state.root_axis_count != state.deepest_layer_by_axis.len or
        state.active_root_axis_count.* != state.root_axis_count or
        !std.math.isFinite(seed_depth_m) or seed_depth_m < 0)
        return error.InvalidDeadRootDepthResetDimensions;
    const axis_domains = std.math.mul(
        usize,
        state.root_axis_count,
        state.root_domain_count,
    ) catch return error.InvalidDeadRootDepthResetDimensions;
    if (state.primary_depth_from_surface_m_by_axis_domain.len != axis_domains or
        state.primary_total_by_axis_domain.len != axis_domains)
        return error.InvalidDeadRootDepthResetDimensions;
    for (state.primary_depth_from_surface_m_by_axis_domain) |depth|
        if (!std.math.isFinite(depth) or depth < 0)
            return error.InvalidDeadRootDepthResetInput;
    for (state.primary_total_by_axis_domain) |mass|
        try validateDeadRootDepthMass(mass);
    if (!roots_dead) return;

    for (0..state.root_axis_count) |axis| {
        state.deepest_layer_by_axis[axis] = planting_layer_index;
        for (0..state.root_domain_count) |domain| {
            const axis_domain = axis * state.root_domain_count + domain;
            state.primary_depth_from_surface_m_by_axis_domain[axis_domain] =
                seed_depth_m;
            state.primary_total_by_axis_domain[axis_domain] = .{};
        }
    }
    state.deepest_active_root_layer.* = planting_layer_index;
    state.active_root_axis_count.* = 0;
}

pub const SourceOrderCompleteDeathBranchPools = struct {
    host_mobile: canopy.ElementalMass,
    symbiont_mobile: canopy.ElementalMass,
    c4_intermediate_carbon_g_c: f64,
    leaf: canopy.ElementalMass,
    symbiont_structural: canopy.ElementalMass,
    sheath: canopy.ElementalMass,
    husk: canopy.ElementalMass,
    ear: canopy.ElementalMass,
    grain: canopy.ElementalMass,
    stalk: canopy.ElementalMass,
    stalk_reserve: canopy.ElementalMass,
};

pub const SourceOrderCompleteDeathShootInput = struct {
    shoot_dead: bool,
    roots_dead: bool,
    perennial_growth_habit: bool,
    deciduous_phenology: bool,
    seasonal_storage: canopy.ElementalMass,
    branches: []const SourceOrderCompleteDeathBranchPools,
    root_woody_fraction: TillageElementComposition,
    leaf_woody_fraction: TillageElementComposition,
    sheath_woody_fraction: TillageElementComposition,
    nonstructural_kinetics: litter_partition.ElementFractions,
    foliar_kinetics: litter_partition.ElementFractions,
    nonfoliar_kinetics: litter_partition.ElementFractions,
    stalk_kinetics: litter_partition.ElementFractions,
    coarse_wood_kinetics: litter_partition.ElementFractions,
};

pub const SourceOrderCompleteDeathShootResult = struct {
    planting_layer_woody_litter: [4]canopy.ElementalMass,
    planting_layer_nonwoody_litter: [4]canopy.ElementalMass,
    surface_woody_litter: [4]canopy.ElementalMass,
    surface_nonwoody_litter: [4]canopy.ElementalMass,
    standing_dead_stalk: [4]canopy.ElementalMass,
    seasonal_storage: canopy.ElementalMass,
    plant_death_initialized: bool,
};

/// Exact GROSUB 11443-11518 complete-death storage and shoot litterfall.
pub fn sourceOrderCompleteDeathShootLitterfall(
    input: SourceOrderCompleteDeathShootInput,
) !SourceOrderCompleteDeathShootResult {
    try validateCompleteDeathMass(input.seasonal_storage);
    for (input.branches) |branch| {
        inline for (@typeInfo(SourceOrderCompleteDeathBranchPools).@"struct".fields) |field| {
            const value = @field(branch, field.name);
            if (field.type == f64) {
                if (!std.math.isFinite(value) or value < 0)
                    return error.InvalidCompleteDeathShootLitterfallInput;
            } else try validateCompleteDeathMass(value);
        }
    }
    inline for (.{
        input.root_woody_fraction,
        input.leaf_woody_fraction,
        input.sheath_woody_fraction,
    }) |composition| inline for (@typeInfo(TillageElementComposition).@"struct".fields) |field|
        for (@field(composition, field.name)) |fraction|
            if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
                return error.InvalidCompleteDeathShootLitterfallInput;
    inline for (.{
        input.nonstructural_kinetics,
        input.foliar_kinetics,
        input.nonfoliar_kinetics,
        input.stalk_kinetics,
        input.coarse_wood_kinetics,
    }) |kinetics| inline for (@typeInfo(litter_partition.ElementFractions).@"struct".fields) |field|
        for (@field(kinetics, field.name)) |fraction|
            if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
                return error.InvalidCompleteDeathShootLitterfallInput;

    var result: SourceOrderCompleteDeathShootResult = .{
        .planting_layer_woody_litter = @splat(.{}),
        .planting_layer_nonwoody_litter = @splat(.{}),
        .surface_woody_litter = @splat(.{}),
        .surface_nonwoody_litter = @splat(.{}),
        .standing_dead_stalk = @splat(.{}),
        .seasonal_storage = input.seasonal_storage,
        .plant_death_initialized = false,
    };
    if (!input.shoot_dead or !input.roots_dead) return result;

    const winter_annual = !input.perennial_growth_habit and input.deciduous_phenology;
    if (input.perennial_growth_habit or !input.deciduous_phenology) {
        result.plant_death_initialized = true;
        inline for (0..4) |component| {
            inline for (@typeInfo(canopy.ElementalMass).@"struct".fields, 0..) |element, index| {
                const name =
                    @typeInfo(litter_partition.ElementFractions).@"struct".fields[index].name;
                const storage = @field(input.seasonal_storage, element.name);
                const wood = @field(input.root_woody_fraction, name);
                const kinetic = @field(input.nonstructural_kinetics, name)[component];
                @field(result.planting_layer_woody_litter[component], element.name) +=
                    kinetic * storage * wood[0];
                @field(result.planting_layer_nonwoody_litter[component], element.name) +=
                    kinetic * storage * wood[1];
            }
        }
        result.seasonal_storage = .{};
    }

    inline for (0..4) |component| {
        for (input.branches) |branch| {
            inline for (@typeInfo(canopy.ElementalMass).@"struct".fields, 0..) |element, index| {
                const name =
                    @typeInfo(litter_partition.ElementFractions).@"struct".fields[index].name;
                const leaf_wood = @field(input.leaf_woody_fraction, name);
                const sheath_wood = @field(input.sheath_woody_fraction, name);
                var mobile = @field(branch.host_mobile, element.name) +
                    @field(branch.symbiont_mobile, element.name);
                if (index == 0) mobile += branch.c4_intermediate_carbon_g_c;
                @field(result.surface_nonwoody_litter[component], element.name) +=
                    @field(input.nonstructural_kinetics, name)[component] * mobile;
                @field(result.surface_nonwoody_litter[component], element.name) +=
                    @field(input.foliar_kinetics, name)[component] *
                    (@field(branch.leaf, element.name) * leaf_wood[1] +
                        @field(branch.symbiont_structural, element.name));
                @field(result.surface_nonwoody_litter[component], element.name) +=
                    @field(input.nonfoliar_kinetics, name)[component] *
                    (@field(branch.sheath, element.name) * sheath_wood[1] +
                        @field(branch.husk, element.name) +
                        @field(branch.ear, element.name));
                @field(result.surface_woody_litter[component], element.name) +=
                    @field(input.coarse_wood_kinetics, name)[component] *
                    (@field(branch.leaf, element.name) * leaf_wood[0] +
                        @field(branch.sheath, element.name) * sheath_wood[0]);
                const grain = @field(input.nonfoliar_kinetics, name)[component] *
                    @field(branch.grain, element.name);
                if (winter_annual) {
                    @field(result.seasonal_storage, element.name) += grain;
                } else {
                    @field(result.surface_nonwoody_litter[component], element.name) += grain;
                }
                @field(result.standing_dead_stalk[component], element.name) +=
                    @field(input.stalk_kinetics, name)[component] *
                    (@field(branch.stalk, element.name) +
                        @field(branch.stalk_reserve, element.name));
            }
        }
    }
    inline for (@typeInfo(SourceOrderCompleteDeathShootResult).@"struct".fields) |field| {
        if (field.type == bool) continue;
        const value = @field(result, field.name);
        if (field.type == canopy.ElementalMass) {
            try validateCompleteDeathResultMass(value);
        } else for (value) |mass| try validateCompleteDeathResultMass(mass);
    }
    return result;
}

pub const SourceOrderCompleteDeathRootInput = struct {
    shoot_dead: bool,
    roots_dead: bool,
    root_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    mobile_by_domain_layer: []const canopy.ElementalMass,
    structural_by_domain_layer_axis: []const SourceOrderDeadRootAxisPools,
    root_woody_fraction: TillageElementComposition,
    nonstructural_kinetics: litter_partition.ElementFractions,
    fine_root_kinetics: litter_partition.ElementFractions,
    coarse_root_kinetics: litter_partition.ElementFractions,
};

/// Exact GROSUB 11535-11557 complete-death root litterfall.
/// The caller supplies only the active NU:NJ layer range and owns the result.
pub fn sourceOrderCompleteDeathRootLitterfall(
    allocator: std.mem.Allocator,
    input: SourceOrderCompleteDeathRootInput,
) ![]SourceOrderDeadRootLayerLitter {
    if (input.root_domain_count == 0 or input.soil_layer_count == 0 or
        input.root_axis_count == 0)
        return error.InvalidCompleteDeathRootLitterfallDimensions;
    const domain_layers = std.math.mul(
        usize,
        input.root_domain_count,
        input.soil_layer_count,
    ) catch return error.InvalidCompleteDeathRootLitterfallDimensions;
    const structural_count = std.math.mul(
        usize,
        domain_layers,
        input.root_axis_count,
    ) catch return error.InvalidCompleteDeathRootLitterfallDimensions;
    if (input.mobile_by_domain_layer.len != domain_layers or
        input.structural_by_domain_layer_axis.len != structural_count)
        return error.InvalidCompleteDeathRootLitterfallDimensions;
    inline for (.{
        input.nonstructural_kinetics,
        input.fine_root_kinetics,
        input.coarse_root_kinetics,
    }) |kinetics| inline for (@typeInfo(litter_partition.ElementFractions).@"struct".fields) |field|
        for (@field(kinetics, field.name)) |fraction|
            if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
                return error.InvalidCompleteDeathRootLitterfallInput;
    inline for (@typeInfo(TillageElementComposition).@"struct".fields) |field|
        for (@field(input.root_woody_fraction, field.name)) |fraction|
            if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
                return error.InvalidCompleteDeathRootLitterfallInput;
    for (input.mobile_by_domain_layer) |mass|
        validateSourceOrderRootMass(mass) catch
            return error.InvalidCompleteDeathRootLitterfallInput;
    for (input.structural_by_domain_layer_axis) |pools| {
        validateSourceOrderRootMass(pools.primary) catch
            return error.InvalidCompleteDeathRootLitterfallInput;
        validateSourceOrderRootMass(pools.secondary) catch
            return error.InvalidCompleteDeathRootLitterfallInput;
    }

    const result = try allocator.alloc(
        SourceOrderDeadRootLayerLitter,
        input.soil_layer_count,
    );
    errdefer allocator.free(result);
    @memset(result, .{ .woody = @splat(.{}), .nonwoody = @splat(.{}) });
    if (!input.shoot_dead or !input.roots_dead) return result;

    inline for (0..4) |component| {
        for (0..input.soil_layer_count) |layer| {
            for (0..input.root_domain_count) |domain| {
                const domain_layer = domain * input.soil_layer_count + layer;
                inline for (@typeInfo(canopy.ElementalMass).@"struct".fields, 0..) |element, index| {
                    const name =
                        @typeInfo(litter_partition.ElementFractions).@"struct".fields[index].name;
                    @field(result[layer].nonwoody[component], element.name) +=
                        @field(input.nonstructural_kinetics, name)[component] *
                        @field(input.mobile_by_domain_layer[domain_layer], element.name);
                }
                for (0..input.root_axis_count) |axis| {
                    const pools = input.structural_by_domain_layer_axis[
                        domain_layer * input.root_axis_count + axis
                    ];
                    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields, 0..) |element, index| {
                        const name =
                            @typeInfo(litter_partition.ElementFractions).@"struct".fields[index].name;
                        const mass = @field(pools.primary, element.name) +
                            @field(pools.secondary, element.name);
                        const woody_fraction = @field(input.root_woody_fraction, name);
                        @field(result[layer].woody[component], element.name) +=
                            @field(input.coarse_root_kinetics, name)[component] *
                            mass * woody_fraction[0];
                        @field(result[layer].nonwoody[component], element.name) +=
                            @field(input.fine_root_kinetics, name)[component] *
                            mass * woody_fraction[1];
                    }
                }
            }
        }
    }
    for (result) |layer| inline for (.{ layer.woody, layer.nonwoody }) |position|
        for (position) |mass| inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
            if (!std.math.isFinite(@field(mass, field.name)))
                return error.NonFiniteCompleteDeathRootLitterfall;
    return result;
}

pub const SourceOrderCompleteDeathBranchState = struct {
    host_mobile: canopy.ElementalMass,
    c4_intermediate_carbon_g_c: f64,
    symbiont_mobile: canopy.ElementalMass,
    shoot: canopy.ElementalMass,
    leaf: canopy.ElementalMass,
    nodule: canopy.ElementalMass,
    sheath: canopy.ElementalMass,
    stalk: canopy.ElementalMass,
    stalk_volume_m3: f64,
    reserve: canopy.ElementalMass,
    husk: canopy.ElementalMass,
    ear: canopy.ElementalMass,
    grain: canopy.ElementalMass,
    leaf_starch_carbon_g_c: f64,
    stalk_extra: canopy.ElementalMass,
};

/// Exact GROSUB 11561-11601 complete-death branch-state reset.
pub fn sourceOrderResetCompleteDeathBranches(
    shoot_dead: bool,
    roots_dead: bool,
    branches: []SourceOrderCompleteDeathBranchState,
) !void {
    for (branches) |branch| {
        inline for (@typeInfo(SourceOrderCompleteDeathBranchState).@"struct".fields) |field| {
            const value = @field(branch, field.name);
            if (field.type == f64) {
                if (!std.math.isFinite(value) or value < 0)
                    return error.InvalidCompleteDeathBranchState;
            } else inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element| {
                const mass = @field(value, element.name);
                if (!std.math.isFinite(mass) or mass < 0)
                    return error.InvalidCompleteDeathBranchState;
            }
        }
    }
    if (!shoot_dead or !roots_dead) return;
    for (branches) |*branch| branch.* = std.mem.zeroes(SourceOrderCompleteDeathBranchState);
}

pub const SourceOrderCompleteDeathRootResetState = struct {
    root_domain_count: usize,
    soil_layer_count: usize,
    root_axis_count: usize,
    mobile_by_domain_layer: []canopy.ElementalMass,
    structural_by_domain_layer_axis: []SourceOrderDeadRootAxisLayerState,
    primary_total_by_domain_axis: []SourceOrderDeadRootDomainAxisState,
};

/// Exact GROSUB 11605-11623 complete-death root-state reset.
pub fn sourceOrderResetCompleteDeathRoots(
    shoot_dead: bool,
    roots_dead: bool,
    state: SourceOrderCompleteDeathRootResetState,
) !void {
    if (state.root_domain_count == 0 or state.soil_layer_count == 0 or
        state.root_axis_count == 0)
        return error.InvalidCompleteDeathRootResetDimensions;
    const domain_layers = std.math.mul(
        usize,
        state.root_domain_count,
        state.soil_layer_count,
    ) catch return error.InvalidCompleteDeathRootResetDimensions;
    const axis_layers = std.math.mul(
        usize,
        domain_layers,
        state.root_axis_count,
    ) catch return error.InvalidCompleteDeathRootResetDimensions;
    const domain_axes = std.math.mul(
        usize,
        state.root_domain_count,
        state.root_axis_count,
    ) catch return error.InvalidCompleteDeathRootResetDimensions;
    if (state.mobile_by_domain_layer.len != domain_layers or
        state.structural_by_domain_layer_axis.len != axis_layers or
        state.primary_total_by_domain_axis.len != domain_axes)
        return error.InvalidCompleteDeathRootResetDimensions;

    for (state.mobile_by_domain_layer) |mass|
        validateDeadRootResetMass(mass) catch
            return error.InvalidCompleteDeathRootResetInput;
    for (state.structural_by_domain_layer_axis) |axis| {
        validateDeadRootResetMass(axis.primary) catch
            return error.InvalidCompleteDeathRootResetInput;
        validateDeadRootResetMass(axis.secondary) catch
            return error.InvalidCompleteDeathRootResetInput;
        inline for (.{ axis.primary_length_m, axis.secondary_length_m, axis.secondary_axis_count }) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidCompleteDeathRootResetInput;
    }
    for (state.primary_total_by_domain_axis) |axis|
        validateDeadRootResetMass(axis.primary_total) catch
            return error.InvalidCompleteDeathRootResetInput;
    if (!shoot_dead or !roots_dead) return;

    for (0..state.soil_layer_count) |layer| {
        for (0..state.root_domain_count) |domain| {
            const domain_layer = domain * state.soil_layer_count + layer;
            state.mobile_by_domain_layer[domain_layer] = .{};
            for (0..state.root_axis_count) |axis| {
                const axis_layer =
                    domain_layer * state.root_axis_count + axis;
                state.structural_by_domain_layer_axis[axis_layer].primary = .{};
                state.structural_by_domain_layer_axis[axis_layer].secondary = .{};
                state.primary_total_by_domain_axis[
                    domain * state.root_axis_count + axis
                ].primary_total = .{};
                state.structural_by_domain_layer_axis[axis_layer].primary_length_m = 0;
                state.structural_by_domain_layer_axis[axis_layer].secondary_length_m = 0;
                state.structural_by_domain_layer_axis[axis_layer].secondary_axis_count = 0;
            }
        }
    }
}

pub const SourceOrderReseedDate = struct {
    day_of_year: u16,
    year: u32,
};

pub const SourceOrderDeadPerennialReseedResult = struct {
    reseed_date: ?SourceOrderReseedDate,
    plant_death_flag: bool,
};

/// Exact GROSUB 11632-11644 dead-perennial reseeding decision.
pub fn sourceOrderScheduleDeadPerennialReseed(
    perennial_growth_habit: bool,
    terminate_on_death: bool,
    current_day_of_year: u16,
    days_in_current_year: u16,
    current_year: u32,
) !SourceOrderDeadPerennialReseedResult {
    if (current_year == 0 or days_in_current_year == 0 or
        current_day_of_year == 0 or current_day_of_year > days_in_current_year)
        return error.InvalidDeadPerennialReseedDate;
    if (!perennial_growth_habit or terminate_on_death) return .{
        .reseed_date = null,
        .plant_death_flag = false,
    };
    if (current_day_of_year < days_in_current_year) return .{
        .reseed_date = .{
            .day_of_year = current_day_of_year + 1,
            .year = current_year,
        },
        .plant_death_flag = true,
    };
    if (current_year == std.math.maxInt(u32))
        return error.DeadPerennialReseedYearOverflow;
    return .{
        .reseed_date = .{ .day_of_year = 1, .year = current_year + 1 },
        .plant_death_flag = true,
    };
}

pub const SourceOrderSoilPlantExchangeInput = struct {
    organic_carbon_exchange_g_c_step: f64,
    organic_nitrogen_exchange_g_n_step: f64,
    ammonium_uptake_g_n_step: f64,
    nitrate_uptake_g_n_step: f64,
    root_fixation_g_n_step: f64,
    canopy_fixation_g_n_step: f64,
    organic_phosphorus_exchange_g_p_step: f64,
    dihydrogen_phosphate_uptake_g_p_step: f64,
    hydrogen_phosphate_uptake_g_p_step: f64,
    cumulative_soil_exchange: canopy.ElementalMass,
    cumulative_fixation_g_n: f64,
    cumulative_plant_carbon_g_c: f64,
    cumulative_respired_carbon_g_c: f64,
};

pub const SourceOrderSoilPlantExchangeResult = struct {
    hourly_net_exchange: canopy.ElementalMass,
    cumulative_soil_exchange: canopy.ElementalMass,
    cumulative_fixation_g_n: f64,
    cumulative_net_primary_productivity_g_c: f64,
};

/// Exact GROSUB 11663-11675 hourly and cumulative soil-plant accounting.
pub fn sourceOrderAccumulateSoilPlantExchange(
    input: SourceOrderSoilPlantExchangeInput,
) !SourceOrderSoilPlantExchangeResult {
    inline for (@typeInfo(SourceOrderSoilPlantExchangeInput).@"struct".fields) |field| {
        const value = @field(input, field.name);
        if (field.type == f64) {
            if (!std.math.isFinite(value))
                return error.InvalidSoilPlantExchange;
        } else inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element|
            if (!std.math.isFinite(@field(value, element.name)))
                return error.InvalidSoilPlantExchange;
    }
    const result: SourceOrderSoilPlantExchangeResult = .{
        .hourly_net_exchange = .{
            .carbon_g = input.organic_carbon_exchange_g_c_step,
            .nitrogen_g = input.organic_nitrogen_exchange_g_n_step +
                input.ammonium_uptake_g_n_step +
                input.nitrate_uptake_g_n_step +
                input.root_fixation_g_n_step,
            .phosphorus_g = input.organic_phosphorus_exchange_g_p_step +
                input.dihydrogen_phosphate_uptake_g_p_step +
                input.hydrogen_phosphate_uptake_g_p_step,
        },
        .cumulative_soil_exchange = .{
            .carbon_g = input.cumulative_soil_exchange.carbon_g +
                input.organic_carbon_exchange_g_c_step,
            .nitrogen_g = input.cumulative_soil_exchange.nitrogen_g +
                input.organic_nitrogen_exchange_g_n_step +
                input.ammonium_uptake_g_n_step +
                input.nitrate_uptake_g_n_step,
            .phosphorus_g = input.cumulative_soil_exchange.phosphorus_g +
                input.organic_phosphorus_exchange_g_p_step +
                input.dihydrogen_phosphate_uptake_g_p_step +
                input.hydrogen_phosphate_uptake_g_p_step,
        },
        .cumulative_fixation_g_n = input.cumulative_fixation_g_n +
            input.root_fixation_g_n_step +
            input.canopy_fixation_g_n_step,
        .cumulative_net_primary_productivity_g_c = input.cumulative_plant_carbon_g_c +
            input.cumulative_respired_carbon_g_c,
    };
    inline for (@typeInfo(SourceOrderSoilPlantExchangeResult).@"struct".fields) |field| {
        const value = @field(result, field.name);
        if (field.type == f64) {
            if (!std.math.isFinite(value))
                return error.NonFiniteSoilPlantExchange;
        } else inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element|
            if (!std.math.isFinite(@field(value, element.name)))
                return error.NonFiniteSoilPlantExchange;
    }
    return result;
}

pub const SourceOrderStandingDeadGeometryInput = struct {
    components: []const canopy.ElementalMass,
    negligible_mass_g_c: f64,
    standing_dead_population_count: f64,
    previous_height_m: f64,
    canopy_height_m: f64,
    stalk_volume_per_carbon_m3_g_c: f64,
    canopy_layer_edges_m: []const f64,
};

pub const SourceOrderStandingDeadGeometryResult = struct {
    total_mass: canopy.ElementalMass,
    height_m: f64,
    total_surface_area_m2: f64,
    layer_surface_area_m2: []f64,
    projected_surface_area_m2: []f64,

    pub fn deinit(self: SourceOrderStandingDeadGeometryResult, allocator: std.mem.Allocator) void {
        allocator.free(self.layer_surface_area_m2);
        allocator.free(self.projected_surface_area_m2);
    }
};

/// Exact GROSUB 11687-11728 standing-dead totals and layer geometry.
pub fn sourceOrderStandingDeadGeometry(
    allocator: std.mem.Allocator,
    input: SourceOrderStandingDeadGeometryInput,
) !SourceOrderStandingDeadGeometryResult {
    if (input.components.len == 0 or input.canopy_layer_edges_m.len < 2)
        return error.InvalidStandingDeadGeometryDimensions;
    inline for (.{
        input.negligible_mass_g_c,
        input.standing_dead_population_count,
        input.previous_height_m,
        input.canopy_height_m,
        input.stalk_volume_per_carbon_m3_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidStandingDeadGeometry;
    for (input.components) |mass| inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidStandingDeadGeometry;
    };
    for (input.canopy_layer_edges_m, 0..) |edge, index| {
        if (!std.math.isFinite(edge) or edge < 0 or
            (index != 0 and edge <= input.canopy_layer_edges_m[index - 1]))
            return error.InvalidStandingDeadGeometry;
    }
    const layer_count = input.canopy_layer_edges_m.len - 1;
    const layer_area = try allocator.alloc(f64, layer_count);
    errdefer allocator.free(layer_area);
    const projected_area = try allocator.alloc(f64, layer_count);
    errdefer allocator.free(projected_area);
    @memset(layer_area, 0);
    @memset(projected_area, 0);

    var total_mass: canopy.ElementalMass = .{};
    for (input.components) |component| {
        total_mass.carbon_g += component.carbon_g;
        total_mass.nitrogen_g += component.nitrogen_g;
        total_mass.phosphorus_g += component.phosphorus_g;
    }
    if (total_mass.carbon_g <= input.negligible_mass_g_c or
        input.standing_dead_population_count <= input.negligible_mass_g_c)
        return .{
            .total_mass = total_mass,
            .height_m = 0,
            .total_surface_area_m2 = 0,
            .layer_surface_area_m2 = layer_area,
            .projected_surface_area_m2 = projected_area,
        };

    const height_m = @max(
        @as(f64, 1.0e-2),
        @max(input.previous_height_m, input.canopy_height_m),
    );
    const radius_m = @sqrt(input.stalk_volume_per_carbon_m3_g_c *
        (@max(@as(f64, 0), total_mass.carbon_g) /
            input.standing_dead_population_count) /
        (3.1416 * height_m));
    const total_surface_area_m2 = 6.2832 * radius_m * height_m *
        input.standing_dead_population_count;
    const domain_top_m = input.canopy_layer_edges_m[layer_count];
    const occupied_height_m = @min(height_m, domain_top_m);
    for (0..layer_count) |layer| {
        const lower_m = input.canopy_layer_edges_m[layer];
        const upper_m = input.canopy_layer_edges_m[layer + 1];
        if (height_m > 0 and lower_m < height_m and upper_m > lower_m) {
            const occupied_fraction = @min(
                @as(f64, 1),
                (height_m - lower_m) / (upper_m - lower_m),
            );
            layer_area[layer] = occupied_fraction * total_surface_area_m2 *
                (upper_m - lower_m) / occupied_height_m;
        }
        projected_area[layer] = 0.25 * layer_area[layer];
    }
    inline for (.{ height_m, total_surface_area_m2 }) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteStandingDeadGeometry;
    for (layer_area, projected_area) |area, projected|
        if (!std.math.isFinite(area) or !std.math.isFinite(projected))
            return error.NonFiniteStandingDeadGeometry;
    return .{
        .total_mass = total_mass,
        .height_m = height_m,
        .total_surface_area_m2 = total_surface_area_m2,
        .layer_surface_area_m2 = layer_area,
        .projected_surface_area_m2 = projected_area,
    };
}

pub const SourceOrderFireShootCarbon = struct {
    canopy_nonstructural_g_c: f64,
    leaf_g_c: f64,
    sheath_g_c: f64,
    stalk_g_c: f64,
    reserve_g_c: f64,
    husk_g_c: f64,
    ear_g_c: f64,
    grain_g_c: f64,
    symbiont_nonstructural_g_c: f64,
    symbiont_biomass_g_c: f64,
    seasonal_storage_g_c: f64,
    standing_dead_g_c: f64,
};

pub const SourceOrderFirePlantInventory = struct {
    shoot: SourceOrderFireShootCarbon,
    canopy_symbiont_included: bool,
    root_domain_count: usize,
    root_axis_count: usize,
    nodule_nonstructural_by_layer_g_c: []const f64,
    nodule_biomass_by_layer_g_c: []const f64,
    root_nonstructural_by_layer_domain_g_c: []const f64,
    root_structural_by_layer_domain_axis: []const SourceOrderDeadRootAxisPools,
};

pub const SourceOrderFireLayerCarbon = struct {
    root_nonstructural_g_c: f64,
    root_structural_g_c: f64,
    nodule_nonstructural_g_c: f64,
    nodule_biomass_g_c: f64,
};

pub const SourceOrderFireInventoryResult = struct {
    shoot: SourceOrderFireShootCarbon,
    layers: []SourceOrderFireLayerCarbon,

    pub fn deinit(self: SourceOrderFireInventoryResult, allocator: std.mem.Allocator) void {
        allocator.free(self.layers);
    }
};

/// Exact GROSUB 11735-11795 fire-event canopy and root C inventory.
pub fn sourceOrderAggregateFireCarbonInventory(
    allocator: std.mem.Allocator,
    fire_in_progress: bool,
    active_layer_count: usize,
    plants: []const SourceOrderFirePlantInventory,
) !?SourceOrderFireInventoryResult {
    if (!fire_in_progress) return null;
    if (active_layer_count == 0) return error.InvalidFireInventoryDimensions;
    for (plants) |plant| {
        if (plant.root_domain_count == 0 or plant.root_axis_count == 0 or
            plant.nodule_nonstructural_by_layer_g_c.len != active_layer_count or
            plant.nodule_biomass_by_layer_g_c.len != active_layer_count)
            return error.InvalidFireInventoryDimensions;
        const domain_layers = std.math.mul(
            usize,
            active_layer_count,
            plant.root_domain_count,
        ) catch return error.InvalidFireInventoryDimensions;
        const axis_layers = std.math.mul(
            usize,
            domain_layers,
            plant.root_axis_count,
        ) catch return error.InvalidFireInventoryDimensions;
        if (plant.root_nonstructural_by_layer_domain_g_c.len != domain_layers or
            plant.root_structural_by_layer_domain_axis.len != axis_layers)
            return error.InvalidFireInventoryDimensions;
        inline for (@typeInfo(SourceOrderFireShootCarbon).@"struct".fields) |field| {
            const value = @field(plant.shoot, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidFireInventory;
        }
        for (plant.nodule_nonstructural_by_layer_g_c) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidFireInventory;
        for (plant.nodule_biomass_by_layer_g_c) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidFireInventory;
        for (plant.root_nonstructural_by_layer_domain_g_c) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidFireInventory;
        for (plant.root_structural_by_layer_domain_axis) |axis| {
            try validateCompleteDeathMass(axis.primary);
            try validateCompleteDeathMass(axis.secondary);
        }
    }

    const layers = try allocator.alloc(SourceOrderFireLayerCarbon, active_layer_count);
    errdefer allocator.free(layers);
    @memset(layers, std.mem.zeroes(SourceOrderFireLayerCarbon));
    var shoot = std.mem.zeroes(SourceOrderFireShootCarbon);
    for (plants) |plant| {
        shoot.canopy_nonstructural_g_c += plant.shoot.canopy_nonstructural_g_c;
        shoot.leaf_g_c += plant.shoot.leaf_g_c;
        shoot.sheath_g_c += plant.shoot.sheath_g_c;
        shoot.stalk_g_c += plant.shoot.stalk_g_c;
        shoot.reserve_g_c += plant.shoot.reserve_g_c;
        shoot.husk_g_c += plant.shoot.husk_g_c;
        shoot.ear_g_c += plant.shoot.ear_g_c;
        shoot.grain_g_c += plant.shoot.grain_g_c;
        if (plant.canopy_symbiont_included) {
            shoot.symbiont_nonstructural_g_c += plant.shoot.symbiont_nonstructural_g_c;
            shoot.symbiont_biomass_g_c += plant.shoot.symbiont_biomass_g_c;
        }
        shoot.standing_dead_g_c += plant.shoot.standing_dead_g_c;
        shoot.seasonal_storage_g_c += plant.shoot.seasonal_storage_g_c;
        for (0..active_layer_count) |layer| {
            layers[layer].nodule_nonstructural_g_c +=
                plant.nodule_nonstructural_by_layer_g_c[layer];
            layers[layer].nodule_biomass_g_c +=
                plant.nodule_biomass_by_layer_g_c[layer];
            for (0..plant.root_domain_count) |domain| {
                const domain_layer = layer * plant.root_domain_count + domain;
                layers[layer].root_nonstructural_g_c +=
                    plant.root_nonstructural_by_layer_domain_g_c[domain_layer];
                for (0..plant.root_axis_count) |axis| {
                    const pools = plant.root_structural_by_layer_domain_axis[
                        domain_layer * plant.root_axis_count + axis
                    ];
                    layers[layer].root_structural_g_c +=
                        pools.primary.carbon_g + pools.secondary.carbon_g;
                }
            }
        }
    }
    inline for (@typeInfo(SourceOrderFireShootCarbon).@"struct".fields) |field|
        if (!std.math.isFinite(@field(shoot, field.name)))
            return error.NonFiniteFireInventory;
    for (layers) |layer| inline for (@typeInfo(SourceOrderFireLayerCarbon).@"struct".fields) |field|
        if (!std.math.isFinite(@field(layer, field.name)))
            return error.NonFiniteFireInventory;
    return .{ .shoot = shoot, .layers = layers };
}

pub const SourceOrderCombustionSpecificRates = struct {
    living_nonstructural_and_leaf_g_c_m2_h: f64,
    living_sheath_g_c_m2_h: f64,
    living_stalk_g_c_m2_h: f64,
    living_reproductive_g_c_m2_h: f64,
    standing_dead_g_c_m2_h: f64,
};

pub const SourceOrderCombustionRates = struct {
    living_temperature_fraction: f64,
    standing_dead_temperature_fraction: f64,
    living_nonstructural_and_leaf_g_c_step: f64,
    living_sheath_g_c_step: f64,
    living_stalk_g_c_step: f64,
    living_reproductive_g_c_step: f64,
    standing_dead_g_c_step: f64,
};

/// Exact GROSUB 11812-11839 living and standing-dead combustion rates.
pub fn sourceOrderCombustionRates(
    canopy_temperature_k: f64,
    standing_dead_temperature_k: f64,
    minimum_combustion_temperature_k: f64,
    maximum_temperature_response: f64,
    surface_area_m2: f64,
    biological_timestep_h: f64,
    specific: SourceOrderCombustionSpecificRates,
) !SourceOrderCombustionRates {
    inline for (.{
        canopy_temperature_k,
        standing_dead_temperature_k,
        minimum_combustion_temperature_k,
        maximum_temperature_response,
        surface_area_m2,
        biological_timestep_h,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCombustionRateInput;
    inline for (@typeInfo(SourceOrderCombustionSpecificRates).@"struct".fields) |field| {
        const value = @field(specific, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCombustionRateInput;
    }
    var result = std.mem.zeroes(SourceOrderCombustionRates);
    if (canopy_temperature_k <= minimum_combustion_temperature_k and
        standing_dead_temperature_k <= minimum_combustion_temperature_k)
        return result;
    if (canopy_temperature_k > minimum_combustion_temperature_k) {
        const gas_constant_temperature = 8.3143 * canopy_temperature_k;
        const response = @min(
            maximum_temperature_response,
            @exp(12.028 - 60000 / gas_constant_temperature),
        );
        result.living_temperature_fraction = @min(@as(f64, 1), response);
        const base_rate = response * surface_area_m2 * biological_timestep_h;
        result.living_nonstructural_and_leaf_g_c_step =
            specific.living_nonstructural_and_leaf_g_c_m2_h * base_rate;
        result.living_sheath_g_c_step =
            specific.living_sheath_g_c_m2_h * base_rate;
        result.living_stalk_g_c_step =
            specific.living_stalk_g_c_m2_h * base_rate;
        result.living_reproductive_g_c_step =
            specific.living_reproductive_g_c_m2_h * base_rate;
    }
    if (standing_dead_temperature_k > minimum_combustion_temperature_k) {
        const gas_constant_temperature = 8.3143 * standing_dead_temperature_k;
        const response = @min(
            maximum_temperature_response,
            @exp(12.028 - 60000 / gas_constant_temperature),
        );
        result.standing_dead_temperature_fraction = @min(@as(f64, 1), response);
        const base_rate = response * surface_area_m2 * biological_timestep_h;
        result.standing_dead_g_c_step =
            specific.standing_dead_g_c_m2_h * base_rate;
    }
    inline for (@typeInfo(SourceOrderCombustionRates).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.NonFiniteCombustionRate;
    return result;
}

pub const SourceOrderShootCombustionTotals = struct {
    canopy_nonstructural_g_c: f64,
    leaf_g_c: f64,
    sheath_g_c: f64,
    stalk_g_c: f64,
    husk_g_c: f64,
    ear_g_c: f64,
    grain_g_c: f64,
    symbiont_nonstructural_g_c: f64,
    symbiont_biomass_g_c: f64,
    standing_dead_g_c: f64,
};

pub const SourceOrderShootCombustionFractions = struct {
    canopy_nonstructural: f64,
    leaf: f64,
    sheath: f64,
    stalk: f64,
    reserve: f64,
    husk: f64,
    ear: f64,
    grain: f64,
    symbiont_nonstructural: f64,
    symbiont_biomass: f64,
    standing_dead: f64,
};

/// Exact GROSUB 11857-11915 shoot-pool combustion fractions.
pub fn sourceOrderShootCombustionFractions(
    totals: SourceOrderShootCombustionTotals,
    rates: SourceOrderCombustionRates,
    negligible_carbon_g_c: f64,
    biomass_type_code: i32,
    growth_type_code: i32,
) !SourceOrderShootCombustionFractions {
    if (!std.math.isFinite(negligible_carbon_g_c) or negligible_carbon_g_c < 0)
        return error.InvalidShootCombustionFractionInput;
    inline for (@typeInfo(SourceOrderShootCombustionTotals).@"struct".fields) |field| {
        const value = @field(totals, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidShootCombustionFractionInput;
    }
    inline for (@typeInfo(SourceOrderCombustionRates).@"struct".fields) |field| {
        const value = @field(rates, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidShootCombustionFractionInput;
    }
    const low_order_stem_routing =
        biomass_type_code == 0 or growth_type_code <= 1;
    const sheath_rate = if (low_order_stem_routing)
        rates.living_sheath_g_c_step
    else
        rates.living_stalk_g_c_step;
    const stalk_rate = if (low_order_stem_routing)
        rates.living_sheath_g_c_step
    else
        rates.living_reproductive_g_c_step;
    const fraction = struct {
        fn bounded(total_g_c: f64, rate_g_c_step: f64, threshold_g_c: f64) f64 {
            return if (total_g_c > threshold_g_c)
                @min(@as(f64, 1), rate_g_c_step / total_g_c)
            else
                0;
        }
    }.bounded;
    return .{
        .canopy_nonstructural = fraction(
            totals.canopy_nonstructural_g_c,
            rates.living_nonstructural_and_leaf_g_c_step,
            negligible_carbon_g_c,
        ),
        .leaf = fraction(
            totals.leaf_g_c,
            rates.living_nonstructural_and_leaf_g_c_step,
            negligible_carbon_g_c,
        ),
        .sheath = fraction(totals.sheath_g_c, sheath_rate, negligible_carbon_g_c),
        .stalk = fraction(totals.stalk_g_c, stalk_rate, negligible_carbon_g_c),
        .reserve = fraction(totals.stalk_g_c, stalk_rate, negligible_carbon_g_c),
        .husk = fraction(
            totals.husk_g_c,
            rates.living_nonstructural_and_leaf_g_c_step,
            negligible_carbon_g_c,
        ),
        .ear = fraction(
            totals.ear_g_c,
            rates.living_sheath_g_c_step,
            negligible_carbon_g_c,
        ),
        .grain = fraction(
            totals.grain_g_c,
            rates.living_sheath_g_c_step,
            negligible_carbon_g_c,
        ),
        .symbiont_nonstructural = fraction(
            totals.symbiont_nonstructural_g_c,
            rates.living_nonstructural_and_leaf_g_c_step,
            negligible_carbon_g_c,
        ),
        .symbiont_biomass = fraction(
            totals.symbiont_biomass_g_c,
            rates.living_sheath_g_c_step,
            negligible_carbon_g_c,
        ),
        .standing_dead = fraction(
            totals.standing_dead_g_c,
            rates.standing_dead_g_c_step,
            negligible_carbon_g_c,
        ),
    };
}

pub const SourceOrderShootCombustionBranchPools = struct {
    canopy_nonstructural: canopy.ElementalMass,
    leaf: canopy.ElementalMass,
    sheath: canopy.ElementalMass,
    stalk: canopy.ElementalMass,
    reserve: canopy.ElementalMass,
    husk: canopy.ElementalMass,
    ear: canopy.ElementalMass,
    grain: canopy.ElementalMass,
    symbiont_nonstructural: canopy.ElementalMass,
    symbiont_biomass: canopy.ElementalMass,
};

pub const SourceOrderShootCombustionBranchResult = struct {
    combusted: SourceOrderShootCombustionBranchPools,
    total_combusted: canopy.ElementalMass,
};

pub const SourceOrderShootCombustionResult = struct {
    branches: []SourceOrderShootCombustionBranchResult,
    cumulative_canopy_combustion_g_c: f64,
    disturbance_emission_ledger: canopy.ElementalMass,

    pub fn deinit(self: SourceOrderShootCombustionResult, allocator: std.mem.Allocator) void {
        allocator.free(self.branches);
    }
};

/// Exact GROSUB 11934-11988 per-branch shoot combustion and loss ledgers.
pub fn sourceOrderShootCombustionLosses(
    allocator: std.mem.Allocator,
    pools: []const SourceOrderShootCombustionBranchPools,
    fractions: SourceOrderShootCombustionFractions,
    preceding_cumulative_canopy_combustion_g_c: f64,
    preceding_disturbance_emission_ledger: canopy.ElementalMass,
) !SourceOrderShootCombustionResult {
    if (!std.math.isFinite(preceding_cumulative_canopy_combustion_g_c))
        return error.InvalidShootCombustionLossInput;
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(preceding_disturbance_emission_ledger, field.name)))
            return error.InvalidShootCombustionLossInput;
    inline for (@typeInfo(SourceOrderShootCombustionFractions).@"struct".fields) |field| {
        const value = @field(fractions, field.name);
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return error.InvalidShootCombustionLossInput;
    }
    for (pools) |branch| inline for (@typeInfo(SourceOrderShootCombustionBranchPools).@"struct".fields) |field|
        try validateCompleteDeathMass(@field(branch, field.name));

    const branches = try allocator.alloc(SourceOrderShootCombustionBranchResult, pools.len);
    errdefer allocator.free(branches);
    var cumulative = preceding_cumulative_canopy_combustion_g_c;
    var emissions = preceding_disturbance_emission_ledger;
    for (pools, branches) |branch, *result| {
        result.combusted = .{
            .canopy_nonstructural = scaleElementalMass(
                branch.canopy_nonstructural,
                fractions.canopy_nonstructural,
            ),
            .leaf = scaleElementalMass(branch.leaf, fractions.leaf),
            .sheath = scaleElementalMass(branch.sheath, fractions.sheath),
            .stalk = scaleElementalMass(branch.stalk, fractions.stalk),
            .reserve = scaleElementalMass(branch.reserve, fractions.reserve),
            .husk = scaleElementalMass(branch.husk, fractions.husk),
            .ear = scaleElementalMass(branch.ear, fractions.ear),
            .grain = scaleElementalMass(branch.grain, fractions.grain),
            .symbiont_nonstructural = scaleElementalMass(
                branch.symbiont_nonstructural,
                fractions.symbiont_nonstructural,
            ),
            .symbiont_biomass = scaleElementalMass(
                branch.symbiont_biomass,
                fractions.symbiont_biomass,
            ),
        };
        result.total_combusted = .{};
        inline for (@typeInfo(SourceOrderShootCombustionBranchPools).@"struct".fields) |field| {
            const mass = @field(result.combusted, field.name);
            result.total_combusted.carbon_g += mass.carbon_g;
            result.total_combusted.nitrogen_g += mass.nitrogen_g;
            result.total_combusted.phosphorus_g += mass.phosphorus_g;
        }
        cumulative += result.total_combusted.carbon_g;
        emissions.carbon_g -= result.total_combusted.carbon_g;
        emissions.nitrogen_g -= result.total_combusted.nitrogen_g;
        emissions.phosphorus_g -= result.total_combusted.phosphorus_g;
    }
    if (!std.math.isFinite(cumulative))
        return error.NonFiniteShootCombustionLoss;
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(emissions, field.name)))
            return error.NonFiniteShootCombustionLoss;
    return .{
        .branches = branches,
        .cumulative_canopy_combustion_g_c = cumulative,
        .disturbance_emission_ledger = emissions,
    };
}

pub const SourceOrderShootSaltInventory = struct {
    aluminum_mol: f64,
    iron_mol: f64,
    calcium_mol: f64,
    magnesium_mol: f64,
    sodium_mol: f64,
    potassium_mol: f64,
    sulfate_mol: f64,
    chloride_mol: f64,
};

pub const SourceOrderShootSaltCombustionBranchResult = struct {
    combusted: SourceOrderShootSaltInventory,
    remaining: SourceOrderShootSaltInventory,
};

/// Exact GROSUB 12011-12028 dynamic-salt shoot combustion.
pub fn sourceOrderShootSaltCombustion(
    allocator: std.mem.Allocator,
    dynamic_salt_enabled: bool,
    canopy_nonstructural_combustion_fraction: f64,
    branches: []const SourceOrderShootSaltInventory,
) !?[]SourceOrderShootSaltCombustionBranchResult {
    if (!dynamic_salt_enabled) return null;
    if (!std.math.isFinite(canopy_nonstructural_combustion_fraction) or
        canopy_nonstructural_combustion_fraction < 0 or
        canopy_nonstructural_combustion_fraction > 1)
        return error.InvalidShootSaltCombustionInput;
    for (branches) |branch| inline for (@typeInfo(SourceOrderShootSaltInventory).@"struct".fields) |field| {
        const value = @field(branch, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidShootSaltCombustionInput;
    };
    const result = try allocator.alloc(
        SourceOrderShootSaltCombustionBranchResult,
        branches.len,
    );
    errdefer allocator.free(result);
    for (branches, result) |branch, *branch_result| {
        branch_result.combusted = .{
            .aluminum_mol = branch.aluminum_mol * canopy_nonstructural_combustion_fraction,
            .iron_mol = branch.iron_mol * canopy_nonstructural_combustion_fraction,
            .calcium_mol = branch.calcium_mol * canopy_nonstructural_combustion_fraction,
            .magnesium_mol = branch.magnesium_mol * canopy_nonstructural_combustion_fraction,
            .sodium_mol = branch.sodium_mol * canopy_nonstructural_combustion_fraction,
            .potassium_mol = branch.potassium_mol * canopy_nonstructural_combustion_fraction,
            .sulfate_mol = branch.sulfate_mol * canopy_nonstructural_combustion_fraction,
            .chloride_mol = branch.chloride_mol * canopy_nonstructural_combustion_fraction,
        };
        branch_result.remaining = .{
            .aluminum_mol = branch.aluminum_mol - branch_result.combusted.aluminum_mol,
            .iron_mol = branch.iron_mol - branch_result.combusted.iron_mol,
            .calcium_mol = branch.calcium_mol - branch_result.combusted.calcium_mol,
            .magnesium_mol = branch.magnesium_mol - branch_result.combusted.magnesium_mol,
            .sodium_mol = branch.sodium_mol - branch_result.combusted.sodium_mol,
            .potassium_mol = branch.potassium_mol - branch_result.combusted.potassium_mol,
            .sulfate_mol = branch.sulfate_mol - branch_result.combusted.sulfate_mol,
            .chloride_mol = branch.chloride_mol - branch_result.combusted.chloride_mol,
        };
        inline for (.{ branch_result.combusted, branch_result.remaining }) |inventory|
            inline for (@typeInfo(SourceOrderShootSaltInventory).@"struct".fields) |field|
                if (!std.math.isFinite(@field(inventory, field.name)))
                    return error.NonFiniteShootSaltCombustion;
    }
    return result;
}

pub const SourceOrderShootCombustionNodeState = struct {
    leaf_area_m2: f64,
    sheath_height_m: f64,
    green_leaf: canopy.ElementalMass,
    senescent_leaf_carbon_g_c: f64,
    green_sheath: canopy.ElementalMass,
    senescent_sheath_carbon_g_c: f64,
    node: canopy.ElementalMass,
};

pub const SourceOrderShootCombustionNodeLayerState = struct {
    leaf_area_m2: f64,
    green_leaf: canopy.ElementalMass,
};

pub const SourceOrderUncombustedBranchState = struct {
    pools: SourceOrderShootCombustionBranchPools,
    c4_intermediate_carbon_g_c: f64,
    total_shoot: canopy.ElementalMass,
    leaf_area_m2: f64,
    nodes: []SourceOrderShootCombustionNodeState,
    node_layers: []SourceOrderShootCombustionNodeLayerState,
    canopy_layer_count: usize,
};

/// Exact GROSUB 12034-12102 remaining shoot pools and node attributes.
pub fn sourceOrderApplyUncombustedShootState(
    branches: []SourceOrderUncombustedBranchState,
    combusted: []const SourceOrderShootCombustionBranchPools,
    fractions: SourceOrderShootCombustionFractions,
) !void {
    if (branches.len != combusted.len)
        return error.InvalidUncombustedShootDimensions;
    inline for (.{ fractions.leaf, fractions.sheath, fractions.stalk }) |value|
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return error.InvalidUncombustedShootInput;
    for (branches, combusted) |branch, burned| {
        if (branch.canopy_layer_count == 0)
            return error.InvalidUncombustedShootDimensions;
        const node_layer_count = std.math.mul(
            usize,
            branch.nodes.len,
            branch.canopy_layer_count,
        ) catch return error.InvalidUncombustedShootDimensions;
        if (branch.node_layers.len != node_layer_count)
            return error.InvalidUncombustedShootDimensions;
        if (!std.math.isFinite(branch.c4_intermediate_carbon_g_c) or
            branch.c4_intermediate_carbon_g_c < 0 or
            !std.math.isFinite(branch.leaf_area_m2) or branch.leaf_area_m2 < 0)
            return error.InvalidUncombustedShootInput;
        inline for (@typeInfo(SourceOrderShootCombustionBranchPools).@"struct".fields) |field| {
            const pool = @field(branch.pools, field.name);
            const loss = @field(burned, field.name);
            try validateCompleteDeathMass(pool);
            try validateCompleteDeathMass(loss);
            inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element| {
                const remaining = @field(pool, element.name) - @field(loss, element.name);
                if (!std.math.isFinite(remaining) or remaining < 0)
                    return error.InvalidUncombustedShootLoss;
            }
        }
        for (branch.nodes) |node| try validateCombustionNode(node);
        for (branch.node_layers) |layer| {
            if (!std.math.isFinite(layer.leaf_area_m2) or layer.leaf_area_m2 < 0)
                return error.InvalidUncombustedShootInput;
            try validateCompleteDeathMass(layer.green_leaf);
        }
        const remaining = remainingShootPools(branch.pools, burned);
        const total = uncombustedShootTotal(
            remaining,
            branch.c4_intermediate_carbon_g_c,
        );
        try validateCompleteDeathResultMass(total);
    }

    for (branches, combusted) |*branch, burned| {
        branch.pools = remainingShootPools(branch.pools, burned);
        branch.total_shoot = uncombustedShootTotal(
            branch.pools,
            branch.c4_intermediate_carbon_g_c,
        );
        branch.leaf_area_m2 *= 1 - fractions.leaf;
        for (branch.nodes, 0..) |*node, node_index| {
            node.leaf_area_m2 *= 1 - fractions.leaf;
            node.sheath_height_m *= 1 - fractions.sheath;
            node.green_leaf = scaleElementalMass(node.green_leaf, 1 - fractions.leaf);
            node.senescent_leaf_carbon_g_c *= 1 - fractions.leaf;
            node.green_sheath = scaleElementalMass(node.green_sheath, 1 - fractions.sheath);
            node.senescent_sheath_carbon_g_c *= 1 - fractions.sheath;
            node.node = scaleElementalMass(node.node, 1 - fractions.stalk);
            for (0..branch.canopy_layer_count) |layer| {
                const state = &branch.node_layers[
                    node_index * branch.canopy_layer_count + layer
                ];
                state.leaf_area_m2 *= 1 - fractions.leaf;
                state.green_leaf = scaleElementalMass(
                    state.green_leaf,
                    1 - fractions.leaf,
                );
            }
        }
    }
}

pub const SourceOrderStandingDeadCombustionComponent = struct {
    combusted: canopy.ElementalMass,
    remaining: canopy.ElementalMass,
};

pub const SourceOrderStandingDeadCombustionResult = struct {
    components: []SourceOrderStandingDeadCombustionComponent,
    cumulative_standing_dead_combustion_g_c: f64,
    disturbance_emission_ledger: canopy.ElementalMass,

    pub fn deinit(
        self: SourceOrderStandingDeadCombustionResult,
        allocator: std.mem.Allocator,
    ) void {
        allocator.free(self.components);
    }
};

/// Exact GROSUB 12112-12138 non-charcoal standing-dead combustion.
pub fn sourceOrderStandingDeadCombustion(
    allocator: std.mem.Allocator,
    components: []const canopy.ElementalMass,
    combustion_fraction: f64,
    preceding_cumulative_standing_dead_combustion_g_c: f64,
    preceding_disturbance_emission_ledger: canopy.ElementalMass,
) !SourceOrderStandingDeadCombustionResult {
    if (components.len == 0 or !std.math.isFinite(combustion_fraction) or
        combustion_fraction < 0 or combustion_fraction > 1 or
        !std.math.isFinite(preceding_cumulative_standing_dead_combustion_g_c))
        return error.InvalidStandingDeadCombustionInput;
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(preceding_disturbance_emission_ledger, field.name)))
            return error.InvalidStandingDeadCombustionInput;
    for (components) |component| try validateCompleteDeathMass(component);

    const result_components = try allocator.alloc(
        SourceOrderStandingDeadCombustionComponent,
        components.len,
    );
    errdefer allocator.free(result_components);
    var cumulative = preceding_cumulative_standing_dead_combustion_g_c;
    var emissions = preceding_disturbance_emission_ledger;
    for (components, result_components) |component, *result| {
        result.combusted = scaleElementalMass(component, combustion_fraction);
        result.remaining = .{
            .carbon_g = component.carbon_g - result.combusted.carbon_g,
            .nitrogen_g = component.nitrogen_g - result.combusted.nitrogen_g,
            .phosphorus_g = component.phosphorus_g - result.combusted.phosphorus_g,
        };
        cumulative += result.combusted.carbon_g;
        emissions.carbon_g -= result.combusted.carbon_g;
        emissions.nitrogen_g -= result.combusted.nitrogen_g;
        emissions.phosphorus_g -= result.combusted.phosphorus_g;
        try validateCompleteDeathResultMass(result.combusted);
        try validateCompleteDeathResultMass(result.remaining);
    }
    if (!std.math.isFinite(cumulative))
        return error.NonFiniteStandingDeadCombustion;
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(emissions, field.name)))
            return error.NonFiniteStandingDeadCombustion;
    return .{
        .components = result_components,
        .cumulative_standing_dead_combustion_g_c = cumulative,
        .disturbance_emission_ledger = emissions,
    };
}

pub const SourceOrderCharcoalCombustionResult = struct {
    temperature_response: f64,
    potential_combustion_g_c_step: f64,
    combustion_fraction: f64,
    combusted: canopy.ElementalMass,
    remaining: canopy.ElementalMass,
    cumulative_standing_dead_combustion_g_c: f64,
    disturbance_emission_ledger: canopy.ElementalMass,
    grid_total_combustion_g_c: f64,
};

/// Exact GROSUB 12152-12178 charcoal combustion and loss ledgers.
pub fn sourceOrderCharcoalCombustion(
    standing_dead_temperature_k: f64,
    maximum_temperature_response: f64,
    surface_area_m2: f64,
    biological_timestep_h: f64,
    specific_charcoal_combustion_g_c_m2_h: f64,
    grid_standing_dead_carbon_g_c: f64,
    negligible_carbon_g_c: f64,
    charcoal: canopy.ElementalMass,
    preceding_cumulative_standing_dead_combustion_g_c: f64,
    preceding_disturbance_emission_ledger: canopy.ElementalMass,
    preceding_grid_total_combustion_g_c: f64,
    plant_canopy_combustion_g_c: f64,
) !SourceOrderCharcoalCombustionResult {
    inline for (.{
        standing_dead_temperature_k,
        maximum_temperature_response,
        surface_area_m2,
        biological_timestep_h,
        specific_charcoal_combustion_g_c_m2_h,
        grid_standing_dead_carbon_g_c,
        negligible_carbon_g_c,
        preceding_cumulative_standing_dead_combustion_g_c,
        preceding_grid_total_combustion_g_c,
        plant_canopy_combustion_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidCharcoalCombustionInput;
    if (standing_dead_temperature_k == 0)
        return error.InvalidCharcoalCombustionTemperature;
    try validateCompleteDeathMass(charcoal);
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(preceding_disturbance_emission_ledger, field.name)))
            return error.InvalidCharcoalCombustionInput;

    const gas_constant_temperature = 8.3143 * standing_dead_temperature_k;
    const response = @min(
        maximum_temperature_response,
        @exp(20.620 - 120000 / gas_constant_temperature),
    );
    const base_rate = response * surface_area_m2 * biological_timestep_h;
    const potential = specific_charcoal_combustion_g_c_m2_h * base_rate;
    const fraction = if (grid_standing_dead_carbon_g_c > negligible_carbon_g_c)
        @min(@as(f64, 1), potential / grid_standing_dead_carbon_g_c)
    else
        0;
    const combusted = scaleElementalMass(charcoal, fraction);
    const remaining: canopy.ElementalMass = .{
        .carbon_g = charcoal.carbon_g - combusted.carbon_g,
        .nitrogen_g = charcoal.nitrogen_g - combusted.nitrogen_g,
        .phosphorus_g = charcoal.phosphorus_g - combusted.phosphorus_g,
    };
    const cumulative = preceding_cumulative_standing_dead_combustion_g_c +
        combusted.carbon_g;
    const emissions: canopy.ElementalMass = .{
        .carbon_g = preceding_disturbance_emission_ledger.carbon_g -
            combusted.carbon_g,
        .nitrogen_g = preceding_disturbance_emission_ledger.nitrogen_g -
            combusted.nitrogen_g,
        .phosphorus_g = preceding_disturbance_emission_ledger.phosphorus_g -
            combusted.phosphorus_g,
    };
    const grid_total = preceding_grid_total_combustion_g_c +
        plant_canopy_combustion_g_c + cumulative;
    inline for (.{ response, potential, fraction, cumulative, grid_total }) |value|
        if (!std.math.isFinite(value))
            return error.NonFiniteCharcoalCombustion;
    try validateCompleteDeathResultMass(combusted);
    try validateCompleteDeathResultMass(remaining);
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(emissions, field.name)))
            return error.NonFiniteCharcoalCombustion;
    return .{
        .temperature_response = response,
        .potential_combustion_g_c_step = potential,
        .combustion_fraction = fraction,
        .combusted = combusted,
        .remaining = remaining,
        .cumulative_standing_dead_combustion_g_c = cumulative,
        .disturbance_emission_ledger = emissions,
        .grid_total_combustion_g_c = grid_total,
    };
}

pub const SourceOrderNoCombustionReset = struct {
    canopy_combustion_g_c_step: f64 = 0,
    standing_dead_combustion_g_c_step: f64 = 0,
    canopy_temperature_response: f64 = 0,
    standing_dead_temperature_response: f64 = 0,
    branch_combustion: []SourceOrderShootCombustionBranchPools,
    branch_salt_combustion: ?[]SourceOrderShootSaltInventory,
    standing_dead_combustion: []canopy.ElementalMass,

    pub fn deinit(self: SourceOrderNoCombustionReset, allocator: std.mem.Allocator) void {
        allocator.free(self.branch_combustion);
        if (self.branch_salt_combustion) |salt| allocator.free(salt);
        allocator.free(self.standing_dead_combustion);
    }
};

/// Exact GROSUB 12183-12234 cold-canopy combustion-rate reset.
pub fn sourceOrderResetNoCombustion(
    allocator: std.mem.Allocator,
    branch_count: usize,
    dynamic_salt_enabled: bool,
) !SourceOrderNoCombustionReset {
    const branch_combustion = try allocator.alloc(
        SourceOrderShootCombustionBranchPools,
        branch_count,
    );
    errdefer allocator.free(branch_combustion);
    @memset(branch_combustion, std.mem.zeroes(SourceOrderShootCombustionBranchPools));

    const branch_salt_combustion = if (dynamic_salt_enabled)
        try allocator.alloc(SourceOrderShootSaltInventory, branch_count)
    else
        null;
    errdefer if (branch_salt_combustion) |salt| allocator.free(salt);
    if (branch_salt_combustion) |salt|
        @memset(salt, std.mem.zeroes(SourceOrderShootSaltInventory));

    const standing_dead_combustion = try allocator.alloc(canopy.ElementalMass, 5);
    errdefer allocator.free(standing_dead_combustion);
    @memset(standing_dead_combustion, std.mem.zeroes(canopy.ElementalMass));
    return .{
        .branch_combustion = branch_combustion,
        .branch_salt_combustion = branch_salt_combustion,
        .standing_dead_combustion = standing_dead_combustion,
    };
}

pub const SourceOrderRootCombustionLayerTotals = struct {
    root_nonstructural_g_c: f64,
    active_root_g_c: f64,
    nodule_nonstructural_g_c: f64,
    nodule_biomass_g_c: f64,
};

pub const SourceOrderRootCombustionFractions = struct {
    root_nonstructural: f64,
    active_root: f64,
    nodule_nonstructural: f64,
    nodule_biomass: f64,
    storage: f64,
};

pub const SourceOrderRootCombustionPotentialRates = struct {
    root_nonstructural_g_c_step: f64,
    nodule_biomass_g_c_step: f64,
    active_root_g_c_step: f64,
    reproductive_g_c_step: f64,
    standing_dead_g_c_step: f64,
};

pub const SourceOrderRootStorageCombustionInput = struct {
    soil_temperature_k: f64,
    minimum_combustion_temperature_k: f64,
    maximum_temperature_response: f64,
    surface_area_m2: f64,
    biological_timestep_h: f64,
    specific_rates: SourceOrderCombustionSpecificRates,
    totals: SourceOrderRootCombustionLayerTotals,
    negligible_carbon_g_c: f64,
    is_surface_layer: bool,
    storage: canopy.ElementalMass,
    preceding_layer_combustion_g_c: f64,
    preceding_disturbance_emission_ledger: canopy.ElementalMass,
};

pub const SourceOrderRootStorageCombustionResult = struct {
    temperature_response: f64,
    potential_rates: SourceOrderRootCombustionPotentialRates,
    fractions: SourceOrderRootCombustionFractions,
    storage_combusted: canopy.ElementalMass,
    storage_remaining: canopy.ElementalMass,
    layer_combustion_g_c: f64,
    disturbance_emission_ledger: canopy.ElementalMass,
};

/// Exact GROSUB 12253-12332 soil-temperature rates, fractions, and storage loss.
pub fn sourceOrderRootStorageCombustion(
    input: SourceOrderRootStorageCombustionInput,
) !SourceOrderRootStorageCombustionResult {
    inline for (.{
        input.soil_temperature_k,
        input.minimum_combustion_temperature_k,
        input.maximum_temperature_response,
        input.surface_area_m2,
        input.biological_timestep_h,
        input.negligible_carbon_g_c,
        input.preceding_layer_combustion_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidRootStorageCombustionInput;
    inline for (@typeInfo(SourceOrderCombustionSpecificRates).@"struct".fields) |field| {
        const value = @field(input.specific_rates, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootStorageCombustionInput;
    }
    inline for (@typeInfo(SourceOrderRootCombustionLayerTotals).@"struct".fields) |field| {
        const value = @field(input.totals, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootStorageCombustionInput;
    }
    try validateCompleteDeathMass(input.storage);
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(input.preceding_disturbance_emission_ledger, field.name)))
            return error.InvalidRootStorageCombustionInput;

    var result = std.mem.zeroes(SourceOrderRootStorageCombustionResult);
    result.storage_remaining = input.storage;
    result.layer_combustion_g_c = input.preceding_layer_combustion_g_c;
    result.disturbance_emission_ledger = input.preceding_disturbance_emission_ledger;
    if (input.soil_temperature_k <= input.minimum_combustion_temperature_k)
        return result;

    result.temperature_response = @min(
        input.maximum_temperature_response,
        @exp(12.028 - 60000 / (8.3143 * input.soil_temperature_k)),
    );
    const rate_scale = result.temperature_response *
        input.surface_area_m2 * input.biological_timestep_h;
    result.potential_rates = .{
        .root_nonstructural_g_c_step = input.specific_rates.living_nonstructural_and_leaf_g_c_m2_h * rate_scale,
        .nodule_biomass_g_c_step = input.specific_rates.living_sheath_g_c_m2_h * rate_scale,
        .active_root_g_c_step = input.specific_rates.living_stalk_g_c_m2_h * rate_scale,
        .reproductive_g_c_step = input.specific_rates.living_reproductive_g_c_m2_h * rate_scale,
        .standing_dead_g_c_step = input.specific_rates.standing_dead_g_c_m2_h * rate_scale,
    };
    result.fractions = .{
        .root_nonstructural = boundedCombustionFraction(
            input.totals.root_nonstructural_g_c,
            result.potential_rates.root_nonstructural_g_c_step,
            input.negligible_carbon_g_c,
        ),
        .active_root = boundedCombustionFraction(
            input.totals.active_root_g_c,
            result.potential_rates.active_root_g_c_step,
            input.negligible_carbon_g_c,
        ),
        .nodule_nonstructural = boundedCombustionFraction(
            input.totals.nodule_nonstructural_g_c,
            result.potential_rates.root_nonstructural_g_c_step,
            input.negligible_carbon_g_c,
        ),
        .nodule_biomass = boundedCombustionFraction(
            input.totals.nodule_biomass_g_c,
            result.potential_rates.nodule_biomass_g_c_step,
            input.negligible_carbon_g_c,
        ),
        .storage = 0,
    };
    if (input.is_surface_layer) {
        result.fractions.storage = result.fractions.active_root;
        result.storage_combusted = scaleElementalMass(
            input.storage,
            result.fractions.storage,
        );
        result.storage_remaining = subtractElementalMass(
            input.storage,
            result.storage_combusted,
        );
        result.layer_combustion_g_c += result.storage_combusted.carbon_g;
        result.disturbance_emission_ledger.carbon_g -=
            result.storage_combusted.carbon_g;
        result.disturbance_emission_ledger.nitrogen_g -=
            result.storage_combusted.nitrogen_g;
        result.disturbance_emission_ledger.phosphorus_g -=
            result.storage_combusted.phosphorus_g;
    }
    inline for (@typeInfo(SourceOrderRootCombustionPotentialRates).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.potential_rates, field.name)))
            return error.NonFiniteRootStorageCombustion;
    inline for (@typeInfo(SourceOrderRootCombustionFractions).@"struct".fields) |field| {
        const value = @field(result.fractions, field.name);
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return error.NonFiniteRootStorageCombustion;
    }
    try validateCompleteDeathResultMass(result.storage_combusted);
    try validateCompleteDeathResultMass(result.storage_remaining);
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result.disturbance_emission_ledger, field.name)))
            return error.NonFiniteRootStorageCombustion;
    if (!std.math.isFinite(result.temperature_response) or
        !std.math.isFinite(result.layer_combustion_g_c))
        return error.NonFiniteRootStorageCombustion;
    return result;
}

pub const SourceOrderRootCombustionAxisState = struct {
    primary: canopy.ElementalMass,
    secondary: canopy.ElementalMass,
    whole_primary: canopy.ElementalMass,
    primary_length_m: f64,
    secondary_length_m: f64,
    secondary_root_number: f64,
};

pub const SourceOrderRootCombustionDomainState = struct {
    nonstructural: canopy.ElementalMass,
    salts: SourceOrderShootSaltInventory,
    active_root_carbon_g_c: f64,
    root_density_g_c_m3: f64,
    root_surface_area_m2: f64,
    primary_root_number: f64,
    root_length_m: f64,
    root_length_growth_m_step: f64,
    root_depth_growth_m_step: f64,
    root_volume_growth_m3_step: f64,
    root_volume_m3: f64,
    root_area_m2: f64,
    axes: []SourceOrderRootCombustionAxisState,
};

pub const SourceOrderRootCombustionDomainLoss = struct {
    nonstructural: canopy.ElementalMass,
    salts: ?SourceOrderShootSaltInventory,
    structural: canopy.ElementalMass,
};

pub const SourceOrderRootDomainCombustionResult = struct {
    losses: []SourceOrderRootCombustionDomainLoss,
    layer_combustion_g_c: f64,
    disturbance_emission_ledger: canopy.ElementalMass,

    pub fn deinit(self: SourceOrderRootDomainCombustionResult, allocator: std.mem.Allocator) void {
        allocator.free(self.losses);
    }
};

/// Exact GROSUB 12339-12454 domain/axis root combustion and topology scaling.
pub fn sourceOrderApplyRootDomainCombustion(
    allocator: std.mem.Allocator,
    domains: []SourceOrderRootCombustionDomainState,
    root_nonstructural_fraction: f64,
    active_root_fraction: f64,
    dynamic_salt_enabled: bool,
    preceding_layer_combustion_g_c: f64,
    preceding_disturbance_emission_ledger: canopy.ElementalMass,
) !SourceOrderRootDomainCombustionResult {
    inline for (.{ root_nonstructural_fraction, active_root_fraction }) |value|
        if (!std.math.isFinite(value) or value < 0 or value > 1)
            return error.InvalidRootDomainCombustionInput;
    if (!std.math.isFinite(preceding_layer_combustion_g_c))
        return error.InvalidRootDomainCombustionInput;
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(preceding_disturbance_emission_ledger, field.name)))
            return error.InvalidRootDomainCombustionInput;
    for (domains) |domain| {
        try validateCompleteDeathMass(domain.nonstructural);
        inline for (@typeInfo(SourceOrderShootSaltInventory).@"struct".fields) |field| {
            const value = @field(domain.salts, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidRootDomainCombustionInput;
        }
        inline for (.{
            domain.active_root_carbon_g_c,
            domain.root_density_g_c_m3,
            domain.root_surface_area_m2,
            domain.primary_root_number,
            domain.root_length_m,
            domain.root_length_growth_m_step,
            domain.root_depth_growth_m_step,
            domain.root_volume_growth_m3_step,
            domain.root_volume_m3,
            domain.root_area_m2,
        }) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidRootDomainCombustionInput;
        for (domain.axes) |axis| {
            try validateCompleteDeathMass(axis.primary);
            try validateCompleteDeathMass(axis.secondary);
            try validateCompleteDeathMass(axis.whole_primary);
            inline for (.{
                axis.primary_length_m,
                axis.secondary_length_m,
                axis.secondary_root_number,
            }) |value| if (!std.math.isFinite(value) or value < 0)
                return error.InvalidRootDomainCombustionInput;
        }
    }

    const losses = try allocator.alloc(SourceOrderRootCombustionDomainLoss, domains.len);
    errdefer allocator.free(losses);
    var layer_combustion = preceding_layer_combustion_g_c;
    var emissions = preceding_disturbance_emission_ledger;
    for (domains, losses) |*domain, *loss| {
        loss.nonstructural = scaleElementalMass(
            domain.nonstructural,
            root_nonstructural_fraction,
        );
        domain.nonstructural = subtractElementalMass(
            domain.nonstructural,
            loss.nonstructural,
        );
        layer_combustion += loss.nonstructural.carbon_g;
        emissions.carbon_g -= loss.nonstructural.carbon_g;
        emissions.nitrogen_g -= loss.nonstructural.nitrogen_g;
        emissions.phosphorus_g -= loss.nonstructural.phosphorus_g;
        loss.salts = if (dynamic_salt_enabled)
            scaleSaltInventory(domain.salts, root_nonstructural_fraction)
        else
            null;
        if (loss.salts) |salt_loss|
            domain.salts = subtractSaltInventory(domain.salts, salt_loss);

        const remaining_active_fraction = 1 - active_root_fraction;
        inline for (.{
            &domain.active_root_carbon_g_c,
            &domain.root_density_g_c_m3,
            &domain.root_surface_area_m2,
            &domain.primary_root_number,
            &domain.root_length_m,
            &domain.root_length_growth_m_step,
            &domain.root_depth_growth_m_step,
            &domain.root_volume_growth_m3_step,
            &domain.root_volume_m3,
            &domain.root_area_m2,
        }) |attribute| attribute.* *= remaining_active_fraction;
        loss.structural = .{};
        for (domain.axes) |*axis| {
            const primary_loss = scaleElementalMass(axis.primary, active_root_fraction);
            const secondary_loss = scaleElementalMass(axis.secondary, active_root_fraction);
            axis.primary = subtractElementalMass(axis.primary, primary_loss);
            axis.secondary = subtractElementalMass(axis.secondary, secondary_loss);
            inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
                @field(loss.structural, field.name) +=
                    @field(primary_loss, field.name) + @field(secondary_loss, field.name);
                @field(axis.whole_primary, field.name) *= remaining_active_fraction;
            }
            axis.primary_length_m *= remaining_active_fraction;
            axis.secondary_length_m *= remaining_active_fraction;
            axis.secondary_root_number *= remaining_active_fraction;
        }
        layer_combustion += loss.structural.carbon_g;
        emissions.carbon_g -= loss.structural.carbon_g;
        emissions.nitrogen_g -= loss.structural.nitrogen_g;
        emissions.phosphorus_g -= loss.structural.phosphorus_g;
    }
    if (!std.math.isFinite(layer_combustion))
        return error.NonFiniteRootDomainCombustion;
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(emissions, field.name)))
            return error.NonFiniteRootDomainCombustion;
    return .{
        .losses = losses,
        .layer_combustion_g_c = layer_combustion,
        .disturbance_emission_ledger = emissions,
    };
}

pub const SourceOrderRootNoduleCombustionResult = struct {
    nonstructural_combusted: canopy.ElementalMass,
    nonstructural_remaining: canopy.ElementalMass,
    biomass_combusted: canopy.ElementalMass,
    biomass_remaining: canopy.ElementalMass,
    layer_plant_combustion_g_c: f64,
    grid_layer_combustion_g_c: f64,
    disturbance_emission_ledger: canopy.ElementalMass,
};

/// Exact GROSUB 12461-12492 root-nodule combustion and layer/grid ledgers.
pub fn sourceOrderRootNoduleCombustion(
    nonstructural: canopy.ElementalMass,
    biomass: canopy.ElementalMass,
    nonstructural_fraction: f64,
    biomass_fraction: f64,
    preceding_layer_plant_combustion_g_c: f64,
    preceding_grid_layer_combustion_g_c: f64,
    preceding_disturbance_emission_ledger: canopy.ElementalMass,
) !SourceOrderRootNoduleCombustionResult {
    try validateCompleteDeathMass(nonstructural);
    try validateCompleteDeathMass(biomass);
    inline for (.{ nonstructural_fraction, biomass_fraction }) |fraction|
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
            return error.InvalidRootNoduleCombustionInput;
    inline for (.{
        preceding_layer_plant_combustion_g_c,
        preceding_grid_layer_combustion_g_c,
    }) |value| if (!std.math.isFinite(value))
        return error.InvalidRootNoduleCombustionInput;
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(preceding_disturbance_emission_ledger, field.name)))
            return error.InvalidRootNoduleCombustionInput;

    const nonstructural_loss = scaleElementalMass(nonstructural, nonstructural_fraction);
    const biomass_loss = scaleElementalMass(biomass, biomass_fraction);
    const combined_loss: canopy.ElementalMass = .{
        .carbon_g = nonstructural_loss.carbon_g + biomass_loss.carbon_g,
        .nitrogen_g = nonstructural_loss.nitrogen_g + biomass_loss.nitrogen_g,
        .phosphorus_g = nonstructural_loss.phosphorus_g + biomass_loss.phosphorus_g,
    };
    const plant_layer = preceding_layer_plant_combustion_g_c + combined_loss.carbon_g;
    const grid_layer = preceding_grid_layer_combustion_g_c + plant_layer;
    const emissions: canopy.ElementalMass = .{
        .carbon_g = preceding_disturbance_emission_ledger.carbon_g -
            combined_loss.carbon_g,
        .nitrogen_g = preceding_disturbance_emission_ledger.nitrogen_g -
            combined_loss.nitrogen_g,
        .phosphorus_g = preceding_disturbance_emission_ledger.phosphorus_g -
            combined_loss.phosphorus_g,
    };
    const result: SourceOrderRootNoduleCombustionResult = .{
        .nonstructural_combusted = nonstructural_loss,
        .nonstructural_remaining = subtractElementalMass(nonstructural, nonstructural_loss),
        .biomass_combusted = biomass_loss,
        .biomass_remaining = subtractElementalMass(biomass, biomass_loss),
        .layer_plant_combustion_g_c = plant_layer,
        .grid_layer_combustion_g_c = grid_layer,
        .disturbance_emission_ledger = emissions,
    };
    try validateCompleteDeathResultMass(result.nonstructural_remaining);
    try validateCompleteDeathResultMass(result.biomass_remaining);
    if (!std.math.isFinite(plant_layer) or !std.math.isFinite(grid_layer))
        return error.NonFiniteRootNoduleCombustion;
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(emissions, field.name)))
            return error.NonFiniteRootNoduleCombustion;
    return result;
}

pub const SourceOrderRootAxisCombustionReset = struct {
    primary: canopy.ElementalMass,
    secondary: canopy.ElementalMass,
};

pub const SourceOrderColdSoilCombustionReset = struct {
    layer_plant_combustion_g_c: f64 = 0,
    storage_combustion: ?canopy.ElementalMass,
    domains: []SourceOrderRootCombustionDomainLoss,
    axis_offsets: []usize,
    axes: []SourceOrderRootAxisCombustionReset,
    nodule_nonstructural_combustion: canopy.ElementalMass = .{},
    nodule_biomass_combustion: canopy.ElementalMass = .{},

    pub fn deinit(self: SourceOrderColdSoilCombustionReset, allocator: std.mem.Allocator) void {
        allocator.free(self.domains);
        allocator.free(self.axis_offsets);
        allocator.free(self.axes);
    }
};

/// Exact GROSUB 12506-12541 cold-soil root/storage/nodule rate reset.
pub fn sourceOrderResetColdSoilCombustion(
    allocator: std.mem.Allocator,
    root_axis_counts_by_domain: []const usize,
    is_surface_layer: bool,
    dynamic_salt_enabled: bool,
) !SourceOrderColdSoilCombustionReset {
    const offset_count = std.math.add(
        usize,
        root_axis_counts_by_domain.len,
        1,
    ) catch return error.InvalidColdSoilCombustionDimensions;
    const axis_offsets = try allocator.alloc(usize, offset_count);
    errdefer allocator.free(axis_offsets);
    axis_offsets[0] = 0;
    for (root_axis_counts_by_domain, 0..) |axis_count, domain|
        axis_offsets[domain + 1] = std.math.add(
            usize,
            axis_offsets[domain],
            axis_count,
        ) catch return error.InvalidColdSoilCombustionDimensions;

    const domains = try allocator.alloc(
        SourceOrderRootCombustionDomainLoss,
        root_axis_counts_by_domain.len,
    );
    errdefer allocator.free(domains);
    for (domains) |*domain| domain.* = .{
        .nonstructural = .{},
        .salts = if (dynamic_salt_enabled)
            std.mem.zeroes(SourceOrderShootSaltInventory)
        else
            null,
        .structural = .{},
    };

    const axes = try allocator.alloc(
        SourceOrderRootAxisCombustionReset,
        axis_offsets[axis_offsets.len - 1],
    );
    errdefer allocator.free(axes);
    @memset(axes, std.mem.zeroes(SourceOrderRootAxisCombustionReset));
    return .{
        .storage_combustion = if (is_surface_layer) .{} else null,
        .domains = domains,
        .axis_offsets = axis_offsets,
        .axes = axes,
    };
}

pub const SourceOrderDormantSeedBranch = struct {
    leafout_disabled: bool,
    accumulated_leafout_h: f64,
    required_leafout_h: f64,
};

pub const SourceOrderDormantSeedActivation = struct {
    qualifying_branch_index: usize,
    planting_day_of_year: u16,
    planting_year: i32,
    seeding_depth_m: f64,
    initialization_pending: bool = false,
};

/// Exact pure decision/state transition in GROSUB 12559-12572.
///
/// `initialization_pending` deliberately models legacy `IFLGI == 1`; callers
/// must not substitute the oppositely ordered production lifecycle flag.
pub fn sourceOrderDormantSeedActivation(
    first_subhour_iteration: bool,
    initialization_pending: bool,
    branches: []const SourceOrderDormantSeedBranch,
    current_day_of_year: u16,
    current_year: i32,
    soil_surface_boundary_depth_m: f64,
) !?SourceOrderDormantSeedActivation {
    if (current_day_of_year == 0 or current_day_of_year > 366 or current_year <= 0)
        return error.InvalidDormantSeedActivationDate;
    if (!std.math.isFinite(soil_surface_boundary_depth_m))
        return error.InvalidDormantSeedActivationDepth;
    const seeding_depth_m = 0.005 + soil_surface_boundary_depth_m;
    if (!std.math.isFinite(seeding_depth_m) or seeding_depth_m < 0)
        return error.InvalidDormantSeedActivationDepth;
    for (branches) |branch| {
        if (!std.math.isFinite(branch.accumulated_leafout_h) or
            !std.math.isFinite(branch.required_leafout_h) or
            branch.accumulated_leafout_h < 0 or branch.required_leafout_h < 0)
            return error.InvalidDormantSeedLeafoutHours;
    }
    if (!first_subhour_iteration or !initialization_pending) return null;
    for (branches, 0..) |branch, branch_index| {
        if (!branch.leafout_disabled and
            branch.accumulated_leafout_h >= branch.required_leafout_h)
            return .{
                .qualifying_branch_index = branch_index,
                .planting_day_of_year = current_day_of_year,
                .planting_year = current_year,
                .seeding_depth_m = seeding_depth_m,
            };
    }
    return null;
}

pub const SourceOrderLitterfallAccumulationInput = struct {
    carbon_g_c_by_position_fraction_layer: []const f64,
    nitrogen_g_n_by_position_fraction_layer: []const f64,
    phosphorus_g_p_by_position_fraction_layer: []const f64,
    layer_count_including_surface: usize,
    preceding_cumulative_surface_litter: canopy.ElementalMass,
    preceding_hourly_litter: canopy.ElementalMass,
    preceding_cumulative_litter: canopy.ElementalMass,
    preceding_layer_carbon_g_c: []const f64,
};

pub const SourceOrderLitterfallAccumulationResult = struct {
    cumulative_surface_litter: canopy.ElementalMass,
    hourly_litter: canopy.ElementalMass,
    cumulative_litter: canopy.ElementalMass,
    layer_carbon_g_c: []f64,

    pub fn deinit(self: SourceOrderLitterfallAccumulationResult, allocator: std.mem.Allocator) void {
        allocator.free(self.layer_carbon_g_c);
    }
};

/// Exact GROSUB 12636-12653 position/fraction/layer litter accumulation.
pub fn sourceOrderAccumulateLitterfall(
    allocator: std.mem.Allocator,
    input: SourceOrderLitterfallAccumulationInput,
) !SourceOrderLitterfallAccumulationResult {
    if (input.layer_count_including_surface == 0 or
        input.preceding_layer_carbon_g_c.len != input.layer_count_including_surface)
        return error.InvalidLitterfallAccumulationDimensions;
    const values_per_element = std.math.mul(
        usize,
        2 * 5,
        input.layer_count_including_surface,
    ) catch return error.InvalidLitterfallAccumulationDimensions;
    inline for (.{
        input.carbon_g_c_by_position_fraction_layer,
        input.nitrogen_g_n_by_position_fraction_layer,
        input.phosphorus_g_p_by_position_fraction_layer,
    }) |values| {
        if (values.len != values_per_element)
            return error.InvalidLitterfallAccumulationDimensions;
        for (values) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLitterfallAccumulationInput;
    }
    inline for (.{
        input.preceding_cumulative_surface_litter,
        input.preceding_hourly_litter,
        input.preceding_cumulative_litter,
    }) |ledger| inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(ledger, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLitterfallAccumulationInput;
    };
    for (input.preceding_layer_carbon_g_c) |value|
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidLitterfallAccumulationInput;

    const layer_carbon = try allocator.dupe(f64, input.preceding_layer_carbon_g_c);
    errdefer allocator.free(layer_carbon);
    var surface = input.preceding_cumulative_surface_litter;
    var hourly = input.preceding_hourly_litter;
    var cumulative = input.preceding_cumulative_litter;
    for (0..2) |position| {
        for (0..5) |fraction| {
            const first = (position * 5 + fraction) *
                input.layer_count_including_surface;
            surface.carbon_g += input.carbon_g_c_by_position_fraction_layer[first];
            surface.nitrogen_g += input.nitrogen_g_n_by_position_fraction_layer[first];
            surface.phosphorus_g += input.phosphorus_g_p_by_position_fraction_layer[first];
            for (0..input.layer_count_including_surface) |layer| {
                const index = first + layer;
                const carbon = input.carbon_g_c_by_position_fraction_layer[index];
                const nitrogen = input.nitrogen_g_n_by_position_fraction_layer[index];
                const phosphorus = input.phosphorus_g_p_by_position_fraction_layer[index];
                hourly.carbon_g += carbon;
                hourly.nitrogen_g += nitrogen;
                hourly.phosphorus_g += phosphorus;
                cumulative.carbon_g += carbon;
                cumulative.nitrogen_g += nitrogen;
                cumulative.phosphorus_g += phosphorus;
                layer_carbon[layer] += carbon;
            }
        }
    }
    inline for (.{ surface, hourly, cumulative }) |ledger|
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
            if (!std.math.isFinite(@field(ledger, field.name)))
                return error.NonFiniteLitterfallAccumulation;
    for (layer_carbon) |value| if (!std.math.isFinite(value))
        return error.NonFiniteLitterfallAccumulation;
    return .{
        .cumulative_surface_litter = surface,
        .hourly_litter = hourly,
        .cumulative_litter = cumulative,
        .layer_carbon_g_c = layer_carbon,
    };
}

/// Exact GROSUB 8603-8625 conversion using authoritative combined-canopy
/// ARLFC rather than recomputing it from ARLFT.
pub fn sourceOrderCuttingHeightFromLeafAreaRemoval(
    requested_fraction: f64,
    total_combined_leaf_area_m2: f64,
    boundary_height_m: []const f64,
    combined_leaf_area_m2: []const f64,
    presence_tolerance_m2: f64,
) !f64 {
    if (!std.math.isFinite(requested_fraction) or requested_fraction < 0 or requested_fraction > 1 or
        !std.math.isFinite(total_combined_leaf_area_m2) or total_combined_leaf_area_m2 < 0 or
        !std.math.isFinite(presence_tolerance_m2) or presence_tolerance_m2 < 0 or
        boundary_height_m.len != combined_leaf_area_m2.len + 1 or combined_leaf_area_m2.len == 0)
        return error.InvalidLeafAreaHarvestGeometry;
    for (combined_leaf_area_m2) |area|
        if (!std.math.isFinite(area) or area < 0) return error.InvalidLeafAreaHarvestGeometry;
    for (boundary_height_m, 0..) |height, index| {
        if (!std.math.isFinite(height) or height < 0 or (index > 0 and height < boundary_height_m[index - 1]))
            return error.InvalidLeafAreaHarvestGeometry;
    }
    const target_remaining_leaf_area_m2 = (1 - requested_fraction) * total_combined_leaf_area_m2;
    if (target_remaining_leaf_area_m2 <= 0) return 0;
    var accumulated_leaf_area_m2: f64 = 0;
    var cutting_height_m: f64 = 0;
    for (combined_leaf_area_m2, 0..) |layer_leaf_area_m2, layer| {
        if (boundary_height_m[layer + 1] > boundary_height_m[layer] and
            layer_leaf_area_m2 > presence_tolerance_m2 and
            accumulated_leaf_area_m2 < target_remaining_leaf_area_m2)
        {
            cutting_height_m = if (accumulated_leaf_area_m2 + layer_leaf_area_m2 > target_remaining_leaf_area_m2)
                boundary_height_m[layer] +
                    (target_remaining_leaf_area_m2 - accumulated_leaf_area_m2) / layer_leaf_area_m2 *
                        (boundary_height_m[layer + 1] - boundary_height_m[layer])
            else
                0;
            accumulated_leaf_area_m2 += layer_leaf_area_m2;
        }
    }
    return cutting_height_m;
}

/// Exact GROSUB 8566-8568 disturbance timing selector.
pub fn sourceOrderAbovegroundDisturbanceIsEnabled(
    harvest_code: i8,
    hour_of_day: u8,
    local_solar_noon_h: f64,
) !bool {
    if (hour_of_day > 23 or !std.math.isFinite(local_solar_noon_h) or local_solar_noon_h < 0 or local_solar_noon_h >= 24)
        return error.InvalidAbovegroundDisturbanceSelector;
    const grazing = harvest_code == 4 or harvest_code == 6;
    return grazing or (harvest_code >= 0 and hour_of_day == @as(u8, @intFromFloat(@floor(local_solar_noon_h))));
}

/// Exact GROSUB 8587-8596 population update for non-grazing disturbance.
pub fn sourceOrderPopulationAfterDisturbance(
    terminate_and_reseed: bool,
    thinning_fraction: f64,
    current: PopulationAfterDisturbance,
    reseed_population_per_m2: f64,
    cell_area_m2: f64,
) !PopulationAfterDisturbance {
    inline for (.{
        thinning_fraction,
        current.living_population_per_m2,
        current.living_population_count,
        current.standing_dead_population_count,
        reseed_population_per_m2,
        cell_area_m2,
    }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidDisturbancePopulationInput;
    if (thinning_fraction > 1 or cell_area_m2 <= 0) return error.InvalidDisturbancePopulationInput;
    if (terminate_and_reseed) {
        const population_count = reseed_population_per_m2 * cell_area_m2;
        if (!std.math.isFinite(population_count)) return error.NonFiniteDisturbancePopulation;
        return .{
            .living_population_per_m2 = reseed_population_per_m2,
            .living_population_count = population_count,
            .standing_dead_population_count = population_count,
        };
    }
    const retention = 1 - thinning_fraction;
    return .{
        .living_population_per_m2 = current.living_population_per_m2 * retention,
        .living_population_count = current.living_population_count * retention,
        .standing_dead_population_count = current.standing_dead_population_count * retention,
    };
}
