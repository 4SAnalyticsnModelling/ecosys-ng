const std = @import("std");
const builtin = @import("builtin");
const management = @import("plant_management.zig");
const canopy = @import("canopy_photosynthesis.zig");
const phenology = @import("plant_phenology.zig");
const growth_stages = @import("plant_growth_stages.zig");
const root_system = @import("plant_root_system.zig");
const root_disturbance = @import("plant_root_disturbance.zig");
const symbiotic_fixation = @import("plant_symbiotic_fixation.zig");
const root_litterfall = @import("plant_root_litterfall.zig");
const root_litter_ledger = @import("plant_root_litter_ledger.zig");
const litter_partition = @import("plant_litter_partition.zig");
const soil_organic = @import("soil_organic_initialization.zig");
const grid_module = @import("grid.zig");
const carbon_exchange = @import("canopy_carbon_exchange.zig");
const shoot_litter_bridge = @import("shoot_litter_bridge.zig");
const canopy_structure = @import("canopy_structure.zig");
const canopy_layers = @import("canopy_layer_distribution.zig");
const canopy_biochemistry = @import("canopy_biochemistry.zig");
const dormancy = @import("plant_dormancy.zig");
const grazing_manure = @import("grazing_manure.zig");
const surface_nutrients = @import("organic_matter_fire_exchange.zig");

pub const ScienceParameters = struct {
    nitrogen_fixation_type: u8 = 0,
    carbon_woody_fraction: [2]f64,
    leaf_nitrogen_woody_fraction: [2]f64,
    sheath_nitrogen_woody_fraction: [2]f64,
    leaf_phosphorus_woody_fraction: [2]f64,
    sheath_phosphorus_woody_fraction: [2]f64,
};

pub const ProductLedger = struct {
    direct_litter: canopy.SenescenceProducts = .{},
    nonstructural: canopy.HarvestProducts = .{},
    foliar: canopy.HarvestProducts = .{},
    nonfoliar: canopy.HarvestProducts = .{},
    woody: canopy.HarvestProducts = .{},
    harvested_grain: canopy.ElementalMass = .{},
    standing_dead_export: canopy.ElementalMass = .{},
    standing_dead_charcoal_litter: canopy.ElementalMass = .{},
    manure: grazing_manure.Products = .{},
};

pub const HourlyDisturbanceReset = struct {
    previous_cumulative_harvest_carbon_g_c: f64,
    manure_organic_carbon_g_c: [4]f64,
    manure_organic_nitrogen_g_n: [4]f64,
    manure_organic_phosphorus_g_p: [4]f64,
    manure_inorganic_nitrogen_g_n: f64,
    manure_inorganic_phosphorus_g_p: f64,
};

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

/// Exact GROSUB PPQ/PCUT monthly forest self-thinning equations.
pub fn forestSelfThinningFraction(stem_diameter_m: f64, living_population_per_m2: f64) !f64 {
    inline for (.{ stem_diameter_m, living_population_per_m2 }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidForestSelfThinningInput;
    if (stem_diameter_m == 0 or living_population_per_m2 == 0) return 0;
    const equilibrium_population_per_m2 = 0.1 * std.math.pow(f64, stem_diameter_m / 0.25, -1.6);
    const fraction = @max(0, 0.1 * (living_population_per_m2 - equilibrium_population_per_m2) / living_population_per_m2);
    if (!std.math.isFinite(fraction) or fraction > 1) return error.NonFiniteForestSelfThinning;
    return fraction;
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

pub const TillageElementComposition = struct {
    carbon: [2]f64,
    nitrogen: [2]f64,
    phosphorus: [2]f64,
};

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

fn validateTillageComposition(composition: TillageElementComposition) !void {
    inline for (.{ composition.carbon, composition.nitrogen, composition.phosphorus }) |fractions| {
        var sum: f64 = 0;
        for (fractions) |fraction| {
            if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1)
                return error.InvalidTillageBranchLitterInput;
            sum += fraction;
        }
        if (@abs(sum - 1) > 1.0e-12) return error.InvalidTillageBranchLitterInput;
    }
}

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

fn validateTillageSlice(values: []const f64) !void {
    for (values) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidTillageBranchState;
}

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

pub const source_order_standing_dead_component_count: usize = 5;

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

pub const harvest_product_component_count: usize = 5;

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

fn sumHarvestProductComponents(
    components: [harvest_product_component_count]canopy.ElementalMass,
) !canopy.ElementalMass {
    var total: canopy.ElementalMass = .{};
    for (components) |mass| {
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
            const value = @field(mass, field.name);
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidDisturbanceRemovalInput;
            @field(total, field.name) += value;
        }
    }
    return total;
}

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

pub fn prunedClumpingFactor(current_clumping_factor: f64, pruning_fraction: f64) !f64 {
    inline for (.{ current_clumping_factor, pruning_fraction }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidPruningClumpingFraction;
    const result = current_clumping_factor * pruning_fraction;
    if (!std.math.isFinite(result)) return error.NonFinitePruningClumpingFactor;
    return result;
}

pub const PopulationScaledNumericalThresholds = struct {
    plant_mass_presence_g: f64,
    plant_mass_density_g_m2: f64,
    plant_flux_presence_g_per_step: f64,
};

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

fn validateSourceOrderRootMass(mass: canopy.ElementalMass) !void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDeadRootLitterfallInput;
    }
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

fn validateDeadRootResetMass(mass: canopy.ElementalMass) !void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDeadRootResetInput;
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

fn validateDeadNoduleMass(mass: canopy.ElementalMass) !void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDeadNoduleLitterfallInput;
    }
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

fn validateDeadRootDepthMass(mass: canopy.ElementalMass) !void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDeadRootDepthResetInput;
    }
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

fn validateCompleteDeathMass(mass: canopy.ElementalMass) !void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCompleteDeathShootLitterfallInput;
    }
}

fn validateCompleteDeathResultMass(mass: canopy.ElementalMass) !void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(mass, field.name)))
            return error.NonFiniteCompleteDeathShootLitterfall;
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

fn scaleElementalMass(mass: canopy.ElementalMass, fraction: f64) canopy.ElementalMass {
    return .{
        .carbon_g = mass.carbon_g * fraction,
        .nitrogen_g = mass.nitrogen_g * fraction,
        .phosphorus_g = mass.phosphorus_g * fraction,
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

fn validateCombustionNode(node: SourceOrderShootCombustionNodeState) !void {
    inline for (.{
        node.leaf_area_m2,
        node.sheath_height_m,
        node.senescent_leaf_carbon_g_c,
        node.senescent_sheath_carbon_g_c,
    }) |value| if (!std.math.isFinite(value) or value < 0)
        return error.InvalidUncombustedShootInput;
    try validateCompleteDeathMass(node.green_leaf);
    try validateCompleteDeathMass(node.green_sheath);
    try validateCompleteDeathMass(node.node);
}

fn remainingShootPools(
    pools: SourceOrderShootCombustionBranchPools,
    burned: SourceOrderShootCombustionBranchPools,
) SourceOrderShootCombustionBranchPools {
    var result: SourceOrderShootCombustionBranchPools = undefined;
    inline for (@typeInfo(SourceOrderShootCombustionBranchPools).@"struct".fields) |field| {
        const pool = @field(pools, field.name);
        const loss = @field(burned, field.name);
        @field(result, field.name) = .{
            .carbon_g = pool.carbon_g - loss.carbon_g,
            .nitrogen_g = pool.nitrogen_g - loss.nitrogen_g,
            .phosphorus_g = pool.phosphorus_g - loss.phosphorus_g,
        };
    }
    return result;
}

fn uncombustedShootTotal(
    pools: SourceOrderShootCombustionBranchPools,
    c4_intermediate_carbon_g_c: f64,
) canopy.ElementalMass {
    var total = pools.leaf;
    inline for (.{
        pools.sheath,
        pools.stalk,
        pools.reserve,
        pools.husk,
        pools.ear,
        pools.grain,
        pools.canopy_nonstructural,
    }) |pool| {
        total.carbon_g += pool.carbon_g;
        total.nitrogen_g += pool.nitrogen_g;
        total.phosphorus_g += pool.phosphorus_g;
    }
    total.carbon_g += c4_intermediate_carbon_g_c;
    return total;
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

fn boundedCombustionFraction(total_g_c: f64, rate_g_c_step: f64, threshold_g_c: f64) f64 {
    return if (total_g_c > threshold_g_c)
        @min(@as(f64, 1), rate_g_c_step / total_g_c)
    else
        0;
}

fn subtractElementalMass(
    inventory: canopy.ElementalMass,
    loss: canopy.ElementalMass,
) canopy.ElementalMass {
    return .{
        .carbon_g = inventory.carbon_g - loss.carbon_g,
        .nitrogen_g = inventory.nitrogen_g - loss.nitrogen_g,
        .phosphorus_g = inventory.phosphorus_g - loss.phosphorus_g,
    };
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

fn scaleSaltInventory(
    inventory: SourceOrderShootSaltInventory,
    fraction: f64,
) SourceOrderShootSaltInventory {
    var result: SourceOrderShootSaltInventory = undefined;
    inline for (@typeInfo(SourceOrderShootSaltInventory).@"struct".fields) |field|
        @field(result, field.name) = @field(inventory, field.name) * fraction;
    return result;
}

fn subtractSaltInventory(
    inventory: SourceOrderShootSaltInventory,
    loss: SourceOrderShootSaltInventory,
) SourceOrderShootSaltInventory {
    var result: SourceOrderShootSaltInventory = undefined;
    inline for (@typeInfo(SourceOrderShootSaltInventory).@"struct".fields) |field|
        @field(result, field.name) = @field(inventory, field.name) - @field(loss, field.name);
    return result;
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

pub fn applyForestSelfThinning(context: *Context, plant: usize) !f64 {
    if (plant >= context.canopy_state.plant_population_per_m2.len or plant >= context.canopy_state.plant_stem_diameter_m.len)
        return error.PlantHarvestIndexOutOfBounds;
    const fraction = try forestSelfThinningFraction(
        context.canopy_state.plant_stem_diameter_m[plant],
        context.canopy_state.plant_population_per_m2[plant],
    );
    if (fraction == 0) return 0;
    try applyEvent(context, plant, .{
        .date = .{ .day = 1, .month = 1, .year = 9999 },
        .kind = .none,
        .termination = .retain,
        .cutting_height_m_or_lai_fraction = 1000,
        .thinning_fraction_or_consumption_rate = fraction,
        .harvested_fraction = .{ .leaf = 1, .nonfoliar = 1, .woody = 1, .standing_dead = 0 },
        .ecosystem_export_fraction = .{ .leaf = 0, .nonfoliar = 0, .woody = 0, .standing_dead = 0 },
    });
    return fraction;
}

pub const Context = struct {
    canopy_state: *canopy.State,
    canopy_structure_state: ?*canopy_structure.State = null,
    canopy_layer_state: ?*canopy_layers.State = null,
    branch_development: *phenology.BranchDevelopmentState,
    science_by_plant: []const ScienceParameters,
    products_by_plant: []ProductLedger,
    leaf_area_presence_tolerance_m2: f64,
    plant_structural_presence_threshold_g_per_plant: f64 = 0,
    plant_tissue_presence_threshold_g_per_plant: f64 = 0,
    canopy_biochemistry_parameters_by_plant: ?[]const canopy_biochemistry.Parameters = null,
    plant_phenology: ?*phenology.State = null,
    growth_stages: ?*growth_stages.State = null,
    emerged_by_plant: ?[]bool = null,
    root_state: ?*root_system.State = null,
    root_litter_partition: ?*const litter_partition.State = null,
    root_litter_carbon_ledger: ?*root_litter_ledger.State = null,
    shoot_litter_carbon_g_c_by_plant: ?[]f64 = null,
    shoot_litter_nitrogen_g_n_by_plant: ?[]f64 = null,
    shoot_litter_phosphorus_g_p_by_plant: ?[]f64 = null,
    soil_organic_state: ?*soil_organic.State = null,
    surface_organic_state: ?*soil_organic.State = null,
    surface_nutrient_state: ?*surface_nutrients.State = null,
    daily_manure_carbon_input_g_c: ?[]f64 = null,
    daily_manure_nitrogen_input_g_n: ?[]f64 = null,
    daily_manure_phosphorus_input_g_p: ?[]f64 = null,
    hourly_manure_products_by_plant: ?[]grazing_manure.Products = null,
    grid: ?*const grid_module.GridState = null,
    root_woody_fraction_by_plant: ?[]const f64 = null,
    carbon_exchange_state: ?*carbon_exchange.State = null,
    reseed_population_per_m2_by_plant: ?[]const f64 = null,
    cell_area_m2_by_cell: ?[]const f64 = null,
};

/// Exact GROSUB ARLFY/ARLFR conversion from a negative fractional combined-
/// canopy leaf-area request to a physical cutting height.
pub fn cuttingHeightFromLeafAreaRemoval(
    requested_fraction: f64,
    boundary_height_m: []const f64,
    combined_leaf_area_m2: []const f64,
    presence_tolerance_m2: f64,
) !f64 {
    if (!std.math.isFinite(requested_fraction) or requested_fraction < 0 or requested_fraction > 1 or
        !std.math.isFinite(presence_tolerance_m2) or presence_tolerance_m2 < 0 or
        boundary_height_m.len != combined_leaf_area_m2.len + 1 or combined_leaf_area_m2.len == 0)
        return error.InvalidLeafAreaHarvestGeometry;
    var total_leaf_area_m2: f64 = 0;
    for (combined_leaf_area_m2) |area| {
        if (!std.math.isFinite(area) or area < 0) return error.InvalidLeafAreaHarvestGeometry;
        total_leaf_area_m2 += area;
    }
    return sourceOrderCuttingHeightFromLeafAreaRemoval(
        requested_fraction,
        total_leaf_area_m2,
        boundary_height_m,
        combined_leaf_area_m2,
        presence_tolerance_m2,
    );
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

pub const PopulationAfterDisturbance = struct {
    living_population_per_m2: f64,
    living_population_count: f64,
    standing_dead_population_count: f64,
};

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

fn addHarvestLitterKinetics(
    destination: *canopy.SenescenceProducts,
    mass: canopy.ElementalMass,
    fractions: litter_partition.ElementFractions,
    woody: bool,
) !void {
    inline for (.{ mass.carbon_g, mass.nitrogen_g, mass.phosphorus_g }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidPlantHarvestProduct;
    for (0..4) |kinetic| {
        if (woody) {
            destination.woody_carbon_g[kinetic] += mass.carbon_g * fractions.carbon[kinetic];
            destination.woody_nitrogen_g[kinetic] += mass.nitrogen_g * fractions.nitrogen[kinetic];
            destination.woody_phosphorus_g[kinetic] += mass.phosphorus_g * fractions.phosphorus[kinetic];
        } else {
            destination.nonwoody_carbon_g[kinetic] += mass.carbon_g * fractions.carbon[kinetic];
            destination.nonwoody_nitrogen_g[kinetic] += mass.nitrogen_g * fractions.nitrogen[kinetic];
            destination.nonwoody_phosphorus_g[kinetic] += mass.phosphorus_g * fractions.phosphorus[kinetic];
        }
    }
}

pub fn harvestLitterToKinetics(
    products: ProductLedger,
    nonstructural: litter_partition.ElementFractions,
    foliar: litter_partition.ElementFractions,
    nonfoliar: litter_partition.ElementFractions,
    woody: litter_partition.ElementFractions,
) !canopy.SenescenceProducts {
    try nonstructural.validate();
    try foliar.validate();
    try nonfoliar.validate();
    try woody.validate();
    var result: canopy.SenescenceProducts = products.direct_litter;
    try addHarvestLitterKinetics(&result, products.nonstructural.litter, nonstructural, false);
    try addHarvestLitterKinetics(&result, products.foliar.litter, foliar, false);
    try addHarvestLitterKinetics(&result, products.nonfoliar.litter, nonfoliar, false);
    try addHarvestLitterKinetics(&result, products.woody.litter, woody, true);
    return result;
}

/// Publishes all above-ground harvest litter and returns ecosystem exports.
/// The product ledger is cleared only after the surface transaction succeeds.
pub fn publishPlantProducts(context: *Context, plant: usize) !canopy.ElementalMass {
    if (plant >= context.products_by_plant.len) return error.PlantHarvestIndexOutOfBounds;
    const partitions = context.root_litter_partition orelse return error.IncompletePlantHarvestLitterContext;
    const surface = context.surface_organic_state orelse return error.IncompletePlantHarvestLitterContext;
    const grid = context.grid orelse return error.IncompletePlantHarvestLitterContext;
    if (plant >= partitions.plant_count or context.canopy_state.species_count == 0) return error.PlantHarvestIndexOutOfBounds;
    const cell = plant / context.canopy_state.species_count;
    if (cell >= grid.cell_count) return error.PlantHarvestIndexOutOfBounds;
    const products = context.products_by_plant[plant];
    const litter = try harvestLitterToKinetics(
        products,
        try partitions.get(plant, .nonstructural),
        try partitions.get(plant, .foliar),
        try partitions.get(plant, .non_foliar),
        try partitions.get(plant, .coarse_wood),
    );
    var shoot_litter: canopy.ElementalMass = .{};
    for (0..litter.woody_carbon_g.len) |kinetic| {
        shoot_litter.carbon_g +=
            litter.woody_carbon_g[kinetic] +
            litter.nonwoody_carbon_g[kinetic];
        shoot_litter.nitrogen_g +=
            litter.woody_nitrogen_g[kinetic] +
            litter.nonwoody_nitrogen_g[kinetic];
        shoot_litter.phosphorus_g +=
            litter.woody_phosphorus_g[kinetic] +
            litter.nonwoody_phosphorus_g[kinetic];
        inline for (.{
            shoot_litter.carbon_g,
            shoot_litter.nitrogen_g,
            shoot_litter.phosphorus_g,
        }) |value| if (!std.math.isFinite(value) or value < 0)
            return error.InvalidPlantHarvestProduct;
    }
    addMass(&shoot_litter, products.standing_dead_charcoal_litter);
    const has_litter_publication =
        context.shoot_litter_carbon_g_c_by_plant != null or
        context.shoot_litter_nitrogen_g_n_by_plant != null or
        context.shoot_litter_phosphorus_g_p_by_plant != null;
    const next_shoot_litter: ?canopy.ElementalMass = if (has_litter_publication) blk: {
        const carbon = context.shoot_litter_carbon_g_c_by_plant orelse
            return error.IncompletePlantHarvestLitterPublication;
        const nitrogen = context.shoot_litter_nitrogen_g_n_by_plant orelse
            return error.IncompletePlantHarvestLitterPublication;
        const phosphorus = context.shoot_litter_phosphorus_g_p_by_plant orelse
            return error.IncompletePlantHarvestLitterPublication;
        if (carbon.len != context.products_by_plant.len or
            nitrogen.len != carbon.len or phosphorus.len != carbon.len)
            return error.PlantHarvestLitterPublicationDimensionMismatch;
        const next: canopy.ElementalMass = .{
            .carbon_g = carbon[plant] + shoot_litter.carbon_g,
            .nitrogen_g = nitrogen[plant] + shoot_litter.nitrogen_g,
            .phosphorus_g = phosphorus[plant] + shoot_litter.phosphorus_g,
        };
        inline for (.{ next.carbon_g, next.nitrogen_g, next.phosphorus_g }) |value|
            if (!std.math.isFinite(value) or value < 0)
                return error.InvalidPlantHarvestProduct;
        break :blk next;
    } else null;
    var exported: canopy.ElementalMass = .{};
    addMass(&exported, products.nonstructural.ecosystem_export);
    addMass(&exported, products.foliar.ecosystem_export);
    addMass(&exported, products.nonfoliar.ecosystem_export);
    addMass(&exported, products.woody.ecosystem_export);
    addMass(&exported, products.standing_dead_export);
    inline for (.{ exported.carbon_g, exported.nitrogen_g, exported.phosphorus_g }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidPlantHarvestProduct;
    var manure_carbon_g_c: f64 = 0;
    for (products.manure.organic_by_biochemical_fraction) |mass| manure_carbon_g_c += mass.carbon_g;
    const has_manure = manure_carbon_g_c > 0 or products.manure.inorganic_nitrogen_g_n > 0 or products.manure.inorganic_phosphorus_g_p > 0;
    const nutrients = if (has_manure)
        context.surface_nutrient_state orelse return error.IncompleteGrazingManureContext
    else
        null;
    var manure_nitrogen_g_n = products.manure.inorganic_nitrogen_g_n;
    var manure_phosphorus_g_p = products.manure.inorganic_phosphorus_g_p;
    for (products.manure.organic_by_biochemical_fraction) |mass| {
        manure_nitrogen_g_n += mass.nitrogen_g;
        manure_phosphorus_g_p += mass.phosphorus_g;
    }
    const daily_manure_next: ?canopy.ElementalMass = if (has_manure) blk: {
        const daily_c = context.daily_manure_carbon_input_g_c orelse return error.IncompleteGrazingManureContext;
        const daily_n = context.daily_manure_nitrogen_input_g_n orelse return error.IncompleteGrazingManureContext;
        const daily_p = context.daily_manure_phosphorus_input_g_p orelse return error.IncompleteGrazingManureContext;
        if (cell >= daily_c.len or cell >= daily_n.len or cell >= daily_p.len) return error.GrazingManureDailyLedgerDimensionMismatch;
        const next: canopy.ElementalMass = .{
            .carbon_g = daily_c[cell] + manure_carbon_g_c,
            .nitrogen_g = daily_n[cell] + manure_nitrogen_g_n,
            .phosphorus_g = daily_p[cell] + manure_phosphorus_g_p,
        };
        inline for (.{ next.carbon_g, next.nitrogen_g, next.phosphorus_g }) |value|
            if (!std.math.isFinite(value) or value < 0) return error.GrazingManureDailyLedgerOverflow;
        break :blk next;
    } else null;
    const hourly_manure_next: ?grazing_manure.Products = if (context.hourly_manure_products_by_plant) |hourly| blk: {
        if (hourly.len != context.products_by_plant.len)
            return error.GrazingManureHourlyLedgerDimensionMismatch;
        var next = hourly[plant];
        try grazing_manure.add(&next, products.manure);
        break :blk next;
    } else null;
    if (has_manure) {
        try grazing_manure.validateOrganicCommit(surface, cell, products.manure);
        try nutrients.?.validateSurfaceNutrients(
            cell,
            products.manure.inorganic_nitrogen_g_n / 14.0,
            products.manure.inorganic_phosphorus_g_p / 31.0,
        );
    }
    try shoot_litter_bridge.validateCharcoalCommit(
        surface,
        cell,
        products.standing_dead_charcoal_litter,
    );
    try shoot_litter_bridge.commitCell(surface, cell, litter);
    try shoot_litter_bridge.commitCharcoalCell(
        surface,
        cell,
        products.standing_dead_charcoal_litter,
    );
    if (has_manure) {
        try grazing_manure.commitOrganic(surface, cell, products.manure);
        try nutrients.?.addSurfaceNutrients(
            cell,
            products.manure.inorganic_nitrogen_g_n / 14.0,
            products.manure.inorganic_phosphorus_g_p / 31.0,
        );
        context.daily_manure_carbon_input_g_c.?[cell] = daily_manure_next.?.carbon_g;
        context.daily_manure_nitrogen_input_g_n.?[cell] = daily_manure_next.?.nitrogen_g;
        context.daily_manure_phosphorus_input_g_p.?[cell] = daily_manure_next.?.phosphorus_g;
    }
    if (next_shoot_litter) |next| {
        context.shoot_litter_carbon_g_c_by_plant.?[plant] = next.carbon_g;
        context.shoot_litter_nitrogen_g_n_by_plant.?[plant] = next.nitrogen_g;
        context.shoot_litter_phosphorus_g_p_by_plant.?[plant] = next.phosphorus_g;
    }
    if (hourly_manure_next) |next|
        context.hourly_manure_products_by_plant.?[plant] = next;
    context.products_by_plant[plant] = .{};
    return exported;
}

/// Management-dispatch callback for deterministic cutting, thinning, pruning,
/// and grain harvest. Grazing is rejected here and routed to its demand-driven
/// hourly kernel instead of being approximated as a fractional cut.
pub fn applyEvent(context: *Context, plant: usize, source_event: management.HarvestEvent) !void {
    return applyEventInternal(context, plant, source_event, true, true);
}

/// Publishes every remaining host-root, mycorrhizal, nodule, and root-gas
/// inventory after whole-plant mortality, before later reconstruction can
/// clear the runtime root topology.
pub fn releaseDeadRootsToLitter(context: *Context, plant: usize) !void {
    try applyRootSymbiontHarvest(context, plant, 0);
}

/// Automatic GROSUB winter-annual harvest generated when the main branch
/// enters end-of-season reproductive turnover.
pub fn applyAutomaticSelfSeedingHarvests(
    context: *Context,
    dormancy_state: *const dormancy.RuntimeState,
    growth_habit_by_plant: []const u8,
    leaf_phenology_type_by_plant: []const u8,
) !usize {
    const plant_state = context.plant_phenology orelse return error.IncompleteAutomaticSelfSeedingContext;
    const growth = context.growth_stages orelse return error.IncompleteAutomaticSelfSeedingContext;
    const plant_count = context.canopy_state.plant_branch_offsets.len - 1;
    if (growth_habit_by_plant.len != plant_count or leaf_phenology_type_by_plant.len != plant_count or
        plant_state.active.len != plant_count or plant_state.reseed_pending.len != plant_count or
        dormancy_state.branches.len != growth.branches.len)
        return error.AutomaticSelfSeedingDimensionMismatch;
    var applied: usize = 0;
    for (0..plant_count) |plant| {
        if (!plant_state.active[plant] or plant_state.reseed_pending[plant] or
            growth_habit_by_plant[plant] != 0 or leaf_phenology_type_by_plant[plant] == 0)
            continue;
        const branch = (try growth.mainLivingBranch(plant)) orelse continue;
        if (!dormancy_state.branches[branch].phenological_remobilization_enabled) continue;
        try applyEvent(context, plant, .{
            .date = .{ .day = 1, .month = 1, .year = 9999 },
            .kind = .grain,
            .termination = .terminate_and_reseed,
            .cutting_height_m_or_lai_fraction = 0,
            .thinning_fraction_or_consumption_rate = 0,
            .harvested_fraction = .{ .leaf = 1, .nonfoliar = 1, .woody = 1, .standing_dead = 1 },
            .ecosystem_export_fraction = .{ .leaf = 0, .nonfoliar = 1, .woody = 0, .standing_dead = 0 },
        });
        applied += 1;
    }
    return applied;
}

/// GROSUB spring perennial transition immediately before shoot topology is
/// reconstructed. Old deciduous foliage and all reproductive organs become
/// litter; herbaceous/shrub stalk enters standing dead and stalk reserve is
/// retained in seasonal storage.
pub fn applyStartOfSeasonResidue(
    context: *Context,
    plant: usize,
    biomass_turnover_type: u8,
    root_profile_type: u8,
) !void {
    const state = context.canopy_state;
    if (plant >= context.science_by_plant.len or plant >= context.products_by_plant.len)
        return error.PlantHarvestIndexOutOfBounds;
    const science = context.science_by_plant[plant];
    try validateScience(science);
    const partitions = context.root_litter_partition orelse return error.IncompleteStartOfSeasonResidueContext;
    if (plant >= partitions.plant_count) return error.PlantHarvestIndexOutOfBounds;
    const stalk_kinetics = try partitions.get(plant, .stalk);
    const branches = try state.branchRange(plant);
    for (branches.first..branches.end) |branch| {
        if (biomass_turnover_type == 0) {
            const nodes = try state.nodeRange(branch);
            for (nodes.first..nodes.end) |node| {
                const node_within_branch = node - nodes.first;
                const samples = try state.sampleRange(node);
                for (samples.first..samples.end) |sample| {
                    const products = try canopy.harvestLeafLayerSample(
                        state,
                        branch,
                        node_within_branch,
                        sample - samples.first,
                        .{ .remaining_fraction = 0, .unexported_fraction = 1, .height_below_cut_fraction = 0 },
                        science.carbon_woody_fraction,
                        science.leaf_nitrogen_woody_fraction,
                        science.leaf_phosphorus_woody_fraction,
                        node_within_branch == 1,
                    );
                    addProducts(&context.products_by_plant[plant].foliar, products.foliar);
                    addProducts(&context.products_by_plant[plant].woody, products.woody);
                }
                const sheath = try canopy.harvestNodeSheath(
                    state,
                    branch,
                    node_within_branch,
                    0,
                    1,
                    science.carbon_woody_fraction,
                    science.sheath_nitrogen_woody_fraction,
                    science.sheath_phosphorus_woody_fraction,
                    false,
                    0,
                );
                addProducts(&context.products_by_plant[plant].nonfoliar, sheath.nonwoody);
                addProducts(&context.products_by_plant[plant].woody, sheath.woody);
            }
        }
        const reproductive = try canopy.harvestReproductiveOrgans(state, branch, .{
            .husk_remaining = 0,
            .husk_unexported = 1,
            .ear_remaining = 0,
            .ear_unexported = 1,
            .grain_remaining = 0,
            .grain_unexported = 1,
        });
        addProducts(&context.products_by_plant[plant].nonfoliar, reproductive.products);

        if (biomass_turnover_type == 0 or root_profile_type == 1) {
            const stalk: canopy.ElementalMass = .{
                .carbon_g = state.branch_stalk_carbon_g[branch],
                .nitrogen_g = state.branch_stalk_nitrogen_g[branch],
                .phosphorus_g = state.branch_stalk_phosphorus_g[branch],
            };
            state.plant_standing_dead_carbon_g[plant] += stalk.carbon_g;
            state.plant_standing_dead_nitrogen_g[plant] += stalk.nitrogen_g;
            state.plant_standing_dead_phosphorus_g[plant] += stalk.phosphorus_g;
            for (0..4) |kinetic| {
                const index = plant * 4 + kinetic;
                state.plant_standing_dead_carbon_by_kinetic_g[index] += stalk.carbon_g * stalk_kinetics.carbon[kinetic];
                state.plant_standing_dead_nitrogen_by_kinetic_g[index] += stalk.nitrogen_g * stalk_kinetics.nitrogen[kinetic];
                state.plant_standing_dead_phosphorus_by_kinetic_g[index] += stalk.phosphorus_g * stalk_kinetics.phosphorus[kinetic];
            }
            state.plant_seed_storage_carbon_g[plant] += state.branch_reserve_carbon_g[branch];
            state.plant_seed_storage_nitrogen_g[plant] += state.branch_reserve_nitrogen_g[branch];
            state.plant_seed_storage_phosphorus_g[plant] += state.branch_reserve_phosphorus_g[branch];
            state.branch_stalk_carbon_g[branch] = 0;
            state.branch_stalk_nitrogen_g[branch] = 0;
            state.branch_stalk_phosphorus_g[branch] = 0;
            state.branch_senescing_stalk_carbon_g[branch] = 0;
            state.branch_senescing_stalk_nitrogen_g[branch] = 0;
            state.branch_senescing_stalk_phosphorus_g[branch] = 0;
            state.branch_reserve_carbon_g[branch] = 0;
            state.branch_reserve_nitrogen_g[branch] = 0;
            state.branch_reserve_phosphorus_g[branch] = 0;
            const nodes = try state.nodeRange(branch);
            for (nodes.first..nodes.end) |node| {
                state.node_internode_carbon_g[node] = 0;
                state.node_internode_nitrogen_g[node] = 0;
                state.node_internode_phosphorus_g[node] = 0;
                state.node_internode_length_m[node] = 0;
            }
        }
    }
    try state.validateFinite();
}

/// GROSUB whole-plant death transaction. All living shoot inventories are
/// removed before the dead perennial is scheduled for next-day reconstruction:
/// foliage and reproductive material enter surface litter, stalk plus reserve
/// enter standing dead, and seasonal storage enters nonstructural/woody litter.
pub fn applyWholePlantMortalityResidue(context: *Context, plant: usize) !void {
    const state = context.canopy_state;
    if (plant >= context.science_by_plant.len or plant >= context.products_by_plant.len)
        return error.PlantHarvestIndexOutOfBounds;
    const partitions = context.root_litter_partition orelse return error.IncompletePlantMortalityResidueContext;
    if (plant >= partitions.plant_count) return error.PlantHarvestIndexOutOfBounds;
    _ = try partitions.get(plant, .stalk);
    const science = context.science_by_plant[plant];
    try validateScience(science);

    const storage: canopy.ElementalMass = .{
        .carbon_g = state.plant_seed_storage_carbon_g[plant],
        .nitrogen_g = state.plant_seed_storage_nitrogen_g[plant],
        .phosphorus_g = state.plant_seed_storage_phosphorus_g[plant],
    };
    state.plant_seed_storage_carbon_g[plant] = 0;
    state.plant_seed_storage_nitrogen_g[plant] = 0;
    state.plant_seed_storage_phosphorus_g[plant] = 0;
    const woody_storage: canopy.ElementalMass = .{
        .carbon_g = storage.carbon_g * science.carbon_woody_fraction[0],
        .nitrogen_g = storage.nitrogen_g * science.leaf_nitrogen_woody_fraction[0],
        .phosphorus_g = storage.phosphorus_g * science.leaf_phosphorus_woody_fraction[0],
    };
    const nonwoody_storage: canopy.ElementalMass = .{
        .carbon_g = storage.carbon_g * science.carbon_woody_fraction[1],
        .nitrogen_g = storage.nitrogen_g * science.leaf_nitrogen_woody_fraction[1],
        .phosphorus_g = storage.phosphorus_g * science.leaf_phosphorus_woody_fraction[1],
    };
    const storage_kinetics = try partitions.get(plant, .nonstructural);
    try addHarvestLitterKinetics(&context.products_by_plant[plant].direct_litter, woody_storage, storage_kinetics, true);
    try addHarvestLitterKinetics(&context.products_by_plant[plant].direct_litter, nonwoody_storage, storage_kinetics, false);

    const branches = try state.branchRange(plant);
    for (branches.first..branches.end) |branch| {
        const reserve: canopy.ElementalMass = .{
            .carbon_g = state.branch_reserve_carbon_g[branch],
            .nitrogen_g = state.branch_reserve_nitrogen_g[branch],
            .phosphorus_g = state.branch_reserve_phosphorus_g[branch],
        };
        state.branch_reserve_carbon_g[branch] = 0;
        state.branch_reserve_nitrogen_g[branch] = 0;
        state.branch_reserve_phosphorus_g[branch] = 0;
        state.branch_stalk_carbon_g[branch] += reserve.carbon_g;
        state.branch_stalk_nitrogen_g[branch] += reserve.nitrogen_g;
        state.branch_stalk_phosphorus_g[branch] += reserve.phosphorus_g;

        const mobile = try canopy.harvestBranchMobilePools(state, branch, 0);
        addMass(&context.products_by_plant[plant].nonstructural.litter, mobile);
        const symbiont_mobile: canopy.ElementalMass = .{
            .carbon_g = state.branch_symbiont_mobile_carbon_g[branch],
            .nitrogen_g = state.branch_symbiont_mobile_nitrogen_g[branch],
            .phosphorus_g = state.branch_symbiont_mobile_phosphorus_g[branch],
        };
        const symbiont_structural: canopy.ElementalMass = .{
            .carbon_g = state.branch_symbiont_structural_carbon_g[branch],
            .nitrogen_g = state.branch_symbiont_structural_nitrogen_g[branch],
            .phosphorus_g = state.branch_symbiont_structural_phosphorus_g[branch],
        };
        addMass(&context.products_by_plant[plant].nonstructural.litter, symbiont_mobile);
        addMass(&context.products_by_plant[plant].foliar.litter, symbiont_structural);
        state.branch_symbiont_mobile_carbon_g[branch] = 0;
        state.branch_symbiont_mobile_nitrogen_g[branch] = 0;
        state.branch_symbiont_mobile_phosphorus_g[branch] = 0;
        state.branch_symbiont_structural_carbon_g[branch] = 0;
        state.branch_symbiont_structural_nitrogen_g[branch] = 0;
        state.branch_symbiont_structural_phosphorus_g[branch] = 0;
    }

    // Death removes foliage for every turnover type and routes every stalk.
    try applyStartOfSeasonResidue(context, plant, 0, 1);
    for (branches.first..branches.end) |branch| {
        state.branch_sapwood_carbon_g[branch] = 0;
        state.branch_fixed_carbon_g_c_per_h[branch] = 0;
        state.branch_shoot_carbohydrate_g_c_per_h[branch] = 0;
        state.branch_carboxylation_umol_per_s[branch] = 0;
        state.branch_potential_seed_site_count[branch] = 0;
        state.branch_seed_count[branch] = 0;
        const nodes = try state.nodeRange(branch);
        for (nodes.first..nodes.end) |node| {
            state.node_height_m[node] = 0;
            state.node_sheath_height_m[node] = 0;
            const samples = try state.sampleRange(node);
            @memset(state.sample_stalk_area_m2[samples.first..samples.end], 0);
        }
    }
    state.plant_carboxylation_umol_per_s[plant] = 0;
    state.plant_gross_primary_productivity_g_c_per_h[plant] = 0;
    state.plant_mobile_carbon_g[plant] = 0;
    state.plant_mobile_nitrogen_g[plant] = 0;
    state.plant_mobile_phosphorus_g[plant] = 0;
    state.plant_symbiont_mobile_carbon_g[plant] = 0;
    state.plant_symbiont_mobile_nitrogen_g[plant] = 0;
    state.plant_symbiont_mobile_phosphorus_g[plant] = 0;
    state.plant_total_shoot_carbon_g[plant] = 0;
    state.plant_shoot_growth_g_c_per_step[plant] = 0;
    try state.validateFinite();
}

/// GROSUB dead-branch transaction preceding whole-PFT death. Host mobile C/N/P,
/// C3/C4 intermediates, and stalk reserve return to seasonal storage; remaining
/// foliage and reproductive organs enter litter, stalk enters standing dead,
/// and canopy symbionts enter their source litter classes.
pub fn applyNaturalDeadBranchResidue(
    context: *Context,
    plant: usize,
    branch: usize,
    winter_annual: bool,
) !void {
    const state = context.canopy_state;
    if (plant >= context.science_by_plant.len or plant >= context.products_by_plant.len)
        return error.PlantHarvestIndexOutOfBounds;
    const branches = try state.branchRange(plant);
    if (branch < branches.first or branch >= branches.end) return error.CanopyBranchIndexOutOfBounds;
    const partitions = context.root_litter_partition orelse return error.IncompleteDeadBranchResidueContext;
    const stalk_kinetics = try partitions.get(plant, .stalk);
    const science = context.science_by_plant[plant];
    try validateScience(science);
    if (context.canopy_layer_state) |layers| try layers.clearDeadBranch(state, plant, branch);

    const recovered_mobile = try canopy.harvestBranchMobilePools(state, branch, 0);
    state.plant_seed_storage_carbon_g[plant] += recovered_mobile.carbon_g + state.branch_reserve_carbon_g[branch];
    state.plant_seed_storage_nitrogen_g[plant] += recovered_mobile.nitrogen_g + state.branch_reserve_nitrogen_g[branch];
    state.plant_seed_storage_phosphorus_g[plant] += recovered_mobile.phosphorus_g + state.branch_reserve_phosphorus_g[branch];
    state.branch_reserve_carbon_g[branch] = 0;
    state.branch_reserve_nitrogen_g[branch] = 0;
    state.branch_reserve_phosphorus_g[branch] = 0;

    const symbiont_mobile: canopy.ElementalMass = .{
        .carbon_g = state.branch_symbiont_mobile_carbon_g[branch],
        .nitrogen_g = state.branch_symbiont_mobile_nitrogen_g[branch],
        .phosphorus_g = state.branch_symbiont_mobile_phosphorus_g[branch],
    };
    const symbiont_structural: canopy.ElementalMass = .{
        .carbon_g = state.branch_symbiont_structural_carbon_g[branch],
        .nitrogen_g = state.branch_symbiont_structural_nitrogen_g[branch],
        .phosphorus_g = state.branch_symbiont_structural_phosphorus_g[branch],
    };
    addMass(&context.products_by_plant[plant].nonstructural.litter, symbiont_mobile);
    addMass(&context.products_by_plant[plant].foliar.litter, symbiont_structural);
    state.branch_symbiont_mobile_carbon_g[branch] = 0;
    state.branch_symbiont_mobile_nitrogen_g[branch] = 0;
    state.branch_symbiont_mobile_phosphorus_g[branch] = 0;
    state.branch_symbiont_structural_carbon_g[branch] = 0;
    state.branch_symbiont_structural_nitrogen_g[branch] = 0;
    state.branch_symbiont_structural_phosphorus_g[branch] = 0;

    const nodes = try state.nodeRange(branch);
    for (nodes.first..nodes.end) |node| {
        const node_within_branch = node - nodes.first;
        const samples = try state.sampleRange(node);
        for (samples.first..samples.end) |sample| {
            const leaf = try canopy.harvestLeafLayerSample(
                state,
                branch,
                node_within_branch,
                sample - samples.first,
                .{ .remaining_fraction = 0, .unexported_fraction = 1, .height_below_cut_fraction = 0 },
                science.carbon_woody_fraction,
                science.leaf_nitrogen_woody_fraction,
                science.leaf_phosphorus_woody_fraction,
                node_within_branch == 1,
            );
            addProducts(&context.products_by_plant[plant].foliar, leaf.foliar);
            addProducts(&context.products_by_plant[plant].woody, leaf.woody);
        }
        const sheath = try canopy.harvestNodeSheath(state, branch, node_within_branch, 0, 1, science.carbon_woody_fraction, science.sheath_nitrogen_woody_fraction, science.sheath_phosphorus_woody_fraction, false, 0);
        addProducts(&context.products_by_plant[plant].nonfoliar, sheath.nonwoody);
        addProducts(&context.products_by_plant[plant].woody, sheath.woody);
    }

    const grain: canopy.ElementalMass = .{
        .carbon_g = state.branch_grain_carbon_g[branch],
        .nitrogen_g = state.branch_grain_nitrogen_g[branch],
        .phosphorus_g = state.branch_grain_phosphorus_g[branch],
    };
    const reproductive = try canopy.harvestReproductiveOrgans(state, branch, .{
        .husk_remaining = 0,
        .husk_unexported = 1,
        .ear_remaining = 0,
        .ear_unexported = 1,
        .grain_remaining = 0,
        .grain_unexported = 1,
    });
    addProducts(&context.products_by_plant[plant].nonfoliar, reproductive.products);
    if (winter_annual) {
        subtractMass(&context.products_by_plant[plant].nonfoliar.litter, grain);
        addMassToStorage(state, plant, grain);
    }

    const stalk: canopy.ElementalMass = .{
        .carbon_g = state.branch_stalk_carbon_g[branch],
        .nitrogen_g = state.branch_stalk_nitrogen_g[branch],
        .phosphorus_g = state.branch_stalk_phosphorus_g[branch],
    };
    addMassToStandingDead(state, plant, stalk, stalk_kinetics);
    state.branch_stalk_carbon_g[branch] = 0;
    state.branch_stalk_nitrogen_g[branch] = 0;
    state.branch_stalk_phosphorus_g[branch] = 0;
    state.branch_sapwood_carbon_g[branch] = 0;
    state.branch_senescing_stalk_carbon_g[branch] = 0;
    state.branch_senescing_stalk_nitrogen_g[branch] = 0;
    state.branch_senescing_stalk_phosphorus_g[branch] = 0;
    state.branch_leaf_area_m2[branch] = 0;
    state.branch_potential_seed_site_count[branch] = 0;
    state.branch_seed_count[branch] = 0;
    state.branch_individual_seed_carbon_g[branch] = 0;
    state.branch_c3_feedback_fraction[branch] = 1;
    state.branch_c4_feedback_fraction[branch] = 1;
    try state.validateFinite();
}

/// Late-GROSUB tillage of an eligible herbaceous canopy. The common canopy
/// transaction is reused, but roots are excluded because the disturbance
/// dispatcher commits their layer-resolved transaction separately.
pub fn applyAbovegroundTillage(context: *Context, plant: usize, remaining_fraction: f64, winter_annual: bool) !void {
    if (!std.math.isFinite(remaining_fraction) or remaining_fraction < 0 or remaining_fraction > 1)
        return error.InvalidPlantTillageRetention;
    const state = context.canopy_state;
    if (plant >= state.plant_population_count.len or plant >= context.products_by_plant.len)
        return error.PlantHarvestIndexOutOfBounds;
    const root_fractions = context.root_woody_fraction_by_plant orelse return error.IncompletePlantTillageRootComposition;
    if (plant >= root_fractions.len) return error.PlantHarvestIndexOutOfBounds;
    const root_nonwoody = root_fractions[plant];
    if (!std.math.isFinite(root_nonwoody) or root_nonwoody < 0 or root_nonwoody > 1)
        return error.InvalidPlantTillageRootComposition;
    const standing: canopy.ElementalMass = .{
        .carbon_g = state.plant_standing_dead_carbon_g[plant],
        .nitrogen_g = state.plant_standing_dead_nitrogen_g[plant],
        .phosphorus_g = state.plant_standing_dead_phosphorus_g[plant],
    };
    inline for (.{ standing.carbon_g, standing.nitrogen_g, standing.phosphorus_g }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidPlantTillageStandingDead;
    const branches = try state.branchRange(plant);
    var grain_before: canopy.ElementalMass = .{};
    for (branches.first..branches.end) |branch| {
        grain_before.carbon_g += state.branch_grain_carbon_g[branch];
        grain_before.nitrogen_g += state.branch_grain_nitrogen_g[branch];
        grain_before.phosphorus_g += state.branch_grain_phosphorus_g[branch];
    }
    const storage_before: canopy.ElementalMass = .{
        .carbon_g = state.plant_seed_storage_carbon_g[plant],
        .nitrogen_g = state.plant_seed_storage_nitrogen_g[plant],
        .phosphorus_g = state.plant_seed_storage_phosphorus_g[plant],
    };
    inline for (.{
        grain_before.carbon_g,
        grain_before.nitrogen_g,
        grain_before.phosphorus_g,
        storage_before.carbon_g,
        storage_before.nitrogen_g,
        storage_before.phosphorus_g,
    }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidPlantTillageStorage;
    for (branches.first..branches.end) |branch| inline for (.{
        "branch_symbiont_mobile_carbon_g",
        "branch_symbiont_mobile_nitrogen_g",
        "branch_symbiont_mobile_phosphorus_g",
        "branch_symbiont_structural_carbon_g",
        "branch_symbiont_structural_nitrogen_g",
        "branch_symbiont_structural_phosphorus_g",
    }) |field_name| {
        const value = @field(state, field_name)[branch];
        if (!std.math.isFinite(value) or value < 0) return error.InvalidPlantTillageSymbiont;
    };
    const layers = context.canopy_layer_state;
    if (layers == null and (standing.carbon_g > 0 or standing.nitrogen_g > 0 or standing.phosphorus_g > 0))
        return error.IncompletePlantTillageStandingDeadContext;
    try applyEventInternal(context, plant, .{
        .date = .{ .day = 1, .month = 1, .year = 9999 },
        .kind = .none,
        .termination = if (remaining_fraction == 0) .terminate else .retain,
        .cutting_height_m_or_lai_fraction = 0,
        .thinning_fraction_or_consumption_rate = 1 - remaining_fraction,
        .harvested_fraction = .{ .leaf = 1, .nonfoliar = 1, .woody = 1, .standing_dead = 1 },
        .ecosystem_export_fraction = .{ .leaf = 0, .nonfoliar = 0, .woody = 0, .standing_dead = 0 },
    }, false, false);

    const removed_fraction = 1 - remaining_fraction;
    const grain_to_storage: canopy.ElementalMass = if (winter_annual) .{
        .carbon_g = grain_before.carbon_g * removed_fraction,
        .nitrogen_g = grain_before.nitrogen_g * removed_fraction,
        .phosphorus_g = grain_before.phosphorus_g * removed_fraction,
    } else .{};
    if (winter_annual) {
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
            const next = @field(context.products_by_plant[plant].nonfoliar.litter, field.name) - @field(grain_to_storage, field.name);
            if (!std.math.isFinite(next) or next < -1e-10) return error.InvalidPlantTillageStorage;
            @field(context.products_by_plant[plant].nonfoliar.litter, field.name) = @max(0, next);
        }
    }
    const combined_storage: canopy.ElementalMass = .{
        .carbon_g = storage_before.carbon_g + grain_to_storage.carbon_g,
        .nitrogen_g = storage_before.nitrogen_g + grain_to_storage.nitrogen_g,
        .phosphorus_g = storage_before.phosphorus_g + grain_to_storage.phosphorus_g,
    };
    const removed_storage: canopy.ElementalMass = .{
        .carbon_g = combined_storage.carbon_g * removed_fraction,
        .nitrogen_g = combined_storage.nitrogen_g * removed_fraction,
        .phosphorus_g = combined_storage.phosphorus_g * removed_fraction,
    };
    addScaledMass(&context.products_by_plant[plant].woody.litter, removed_storage, 1 - root_nonwoody);
    addScaledMass(&context.products_by_plant[plant].nonstructural.litter, removed_storage, root_nonwoody);
    state.plant_seed_storage_carbon_g[plant] = combined_storage.carbon_g * remaining_fraction;
    state.plant_seed_storage_nitrogen_g[plant] = combined_storage.nitrogen_g * remaining_fraction;
    state.plant_seed_storage_phosphorus_g[plant] = combined_storage.phosphorus_g * remaining_fraction;
    for (branches.first..branches.end) |branch| {
        const mobile_symbiont: canopy.ElementalMass = .{
            .carbon_g = state.branch_symbiont_mobile_carbon_g[branch] * removed_fraction,
            .nitrogen_g = state.branch_symbiont_mobile_nitrogen_g[branch] * removed_fraction,
            .phosphorus_g = state.branch_symbiont_mobile_phosphorus_g[branch] * removed_fraction,
        };
        const structural_symbiont: canopy.ElementalMass = .{
            .carbon_g = state.branch_symbiont_structural_carbon_g[branch] * removed_fraction,
            .nitrogen_g = state.branch_symbiont_structural_nitrogen_g[branch] * removed_fraction,
            .phosphorus_g = state.branch_symbiont_structural_phosphorus_g[branch] * removed_fraction,
        };
        addMass(&context.products_by_plant[plant].nonstructural.litter, mobile_symbiont);
        addMass(&context.products_by_plant[plant].foliar.litter, structural_symbiont);
        inline for (.{
            "branch_symbiont_mobile_carbon_g",
            "branch_symbiont_mobile_nitrogen_g",
            "branch_symbiont_mobile_phosphorus_g",
            "branch_symbiont_structural_carbon_g",
            "branch_symbiont_structural_nitrogen_g",
            "branch_symbiont_structural_phosphorus_g",
        }) |field_name| @field(state, field_name)[branch] *= remaining_fraction;
    }
    const removed_standing: canopy.ElementalMass = .{
        .carbon_g = standing.carbon_g * removed_fraction,
        .nitrogen_g = standing.nitrogen_g * removed_fraction,
        .phosphorus_g = standing.phosphorus_g * removed_fraction,
    };
    addScaledMass(&context.products_by_plant[plant].woody.litter, removed_standing, 1 - root_nonwoody);
    addScaledMass(&context.products_by_plant[plant].nonfoliar.litter, removed_standing, root_nonwoody);
    inline for (.{
        "plant_standing_dead_carbon_g",
        "plant_standing_dead_nitrogen_g",
        "plant_standing_dead_phosphorus_g",
    }) |field_name| @field(state, field_name)[plant] *= remaining_fraction;
    const kinetic_first = plant * 4;
    inline for (.{
        "plant_standing_dead_carbon_by_kinetic_g",
        "plant_standing_dead_nitrogen_by_kinetic_g",
        "plant_standing_dead_phosphorus_by_kinetic_g",
    }) |field_name| {
        for (@field(state, field_name)[kinetic_first..][0..4]) |*value|
            value.* *= remaining_fraction;
    }
    if (layers) |layer_state| {
        const cell = plant / state.species_count;
        for (0..layer_state.layer_count) |layer| {
            const plant_layer = plant * layer_state.layer_count + layer;
            const removed_area_m2 = layer_state.plant_standing_dead_area_m2[plant_layer] * removed_fraction;
            layer_state.plant_standing_dead_area_m2[plant_layer] *= remaining_fraction;
            layer_state.cell_standing_dead_area_m2[cell * layer_state.layer_count + layer] =
                @max(0, layer_state.cell_standing_dead_area_m2[cell * layer_state.layer_count + layer] - removed_area_m2);
            const projected_first = plant_layer * layer_state.inclination_count;
            for (layer_state.plant_standing_dead_projected_surface_m2[projected_first..][0..layer_state.inclination_count]) |*area_m2|
                area_m2.* *= remaining_fraction;
        }
    }
}

fn applyEventInternal(context: *Context, plant: usize, source_event: management.HarvestEvent, include_roots: bool, include_standing_dead: bool) !void {
    var event = source_event;
    if (event.kind == .animal_grazing or event.kind == .insect_grazing) return error.GrazingRequiresDemandDrivenKernel;
    const state = context.canopy_state;
    if (plant >= context.science_by_plant.len or plant >= context.products_by_plant.len) return error.PlantHarvestIndexOutOfBounds;
    if (@intFromEnum(event.kind) <= @intFromEnum(management.HarvestKind.above_ground) and event.cutting_height_m_or_lai_fraction < 0) {
        const layers = context.canopy_layer_state orelse return error.IncompleteLeafAreaHarvestContext;
        const cell = plant / state.species_count;
        if (cell >= layers.cell_count or layers.species_count != state.species_count) return error.PlantHarvestIndexOutOfBounds;
        const first = cell * layers.layer_count;
        event.cutting_height_m_or_lai_fraction = try cuttingHeightFromLeafAreaRemoval(
            @abs(event.cutting_height_m_or_lai_fraction),
            try layers.cellBoundaries(cell),
            layers.cell_leaf_area_m2[first..][0..layers.layer_count],
            context.leaf_area_presence_tolerance_m2,
        );
    }
    const branches = try state.branchRange(plant);
    if (context.branch_development.branch_count != state.branch_stalk_carbon_g.len) return error.BranchDevelopmentDimensionMismatch;
    const science = context.science_by_plant[plant];
    try validateScience(science);
    const thinning_retention = 1.0 - event.thinning_fraction_or_consumption_rate;
    if (!std.math.isFinite(thinning_retention) or thinning_retention < 0 or thinning_retention > 1)
        return error.InvalidPlantHarvestThinningFraction;
    const pruning = event.kind == .pruning;
    const grain_only = event.kind == .grain;
    const reseed = event.termination == .terminate_and_reseed;
    const reseed_population_per_m2: f64 = if (reseed) blk: {
        const targets = context.reseed_population_per_m2_by_plant orelse return error.IncompletePlantReseedContext;
        const areas = context.cell_area_m2_by_cell orelse return error.IncompletePlantReseedContext;
        const plant_state = context.plant_phenology orelse return error.IncompletePlantReseedContext;
        const cell = plant / state.species_count;
        if (plant >= targets.len or plant >= plant_state.reseed_pending.len or cell >= areas.len) return error.PlantHarvestIndexOutOfBounds;
        if (!std.math.isFinite(targets[plant]) or targets[plant] < 0 or
            !std.math.isFinite(areas[cell]) or areas[cell] <= 0)
            return error.InvalidPlantReseedPopulation;
        break :blk targets[plant];
    } else 0;
    const PruningCommit = struct { structure: *canopy_structure.State, initial: f64, effective: f64 };
    const pruning_commit: ?PruningCommit = if (pruning) blk: {
        const structure = context.canopy_structure_state orelse return error.IncompletePruningContext;
        if (plant >= structure.initial_clumping_factor.len or plant >= structure.effective_clumping_factor.len)
            return error.PlantHarvestIndexOutOfBounds;
        if (!std.math.isFinite(event.cutting_height_m_or_lai_fraction) or event.cutting_height_m_or_lai_fraction < 0)
            return error.InvalidPruningClumpingFraction;
        break :blk .{
            .structure = structure,
            .initial = try prunedClumpingFactor(structure.initial_clumping_factor[plant], event.cutting_height_m_or_lai_fraction),
            .effective = try prunedClumpingFactor(structure.effective_clumping_factor[plant], event.cutting_height_m_or_lai_fraction),
        };
    } else null;
    var maximum_internode_height_m: f64 = 0;
    for (branches.first..branches.end) |branch| {
        const nodes = try state.nodeRange(branch);
        for (state.node_height_m[nodes.first..nodes.end]) |height_m|
            maximum_internode_height_m = @max(maximum_internode_height_m, height_m);
    }
    const reproductive_organs_reached_by_cut =
        event.cutting_height_m_or_lai_fraction < maximum_internode_height_m;
    for (branches.first..branches.end) |branch| {
        if (!grain_only) try harvestVegetativeBranch(context, plant, branch, event, science, pruning);
        var retention = try canopy.reproductiveRetention(false, reproductive_organs_reached_by_cut, grain_only or pruning, event.thinning_fraction_or_consumption_rate, event.harvested_fraction.nonfoliar, state.branch_husk_carbon_g[branch], state.branch_ear_carbon_g[branch], state.branch_grain_carbon_g[branch], 0, 0, 0);
        retention.husk_unexported = unexportedFraction(retention.husk_remaining, event.ecosystem_export_fraction.nonfoliar);
        retention.ear_unexported = unexportedFraction(retention.ear_remaining, event.ecosystem_export_fraction.nonfoliar);
        retention.grain_unexported = unexportedFraction(retention.grain_remaining, event.ecosystem_export_fraction.nonfoliar);
        const reproductive = try canopy.harvestReproductiveOrgans(state, branch, retention);
        addProducts(&context.products_by_plant[plant].nonfoliar, reproductive.products);
        addMass(&context.products_by_plant[plant].harvested_grain, reproductive.harvested_grain);
    }
    if (include_standing_dead) try applyScheduledStandingDeadHarvest(context, plant, event);
    if (include_roots) try applyRootSymbiontHarvest(
        context,
        plant,
        if (event.termination == .retain) thinning_retention else 0,
    );
    if (!reseed) {
        state.plant_population_per_m2[plant] *= thinning_retention;
        state.plant_population_count[plant] *= thinning_retention;
        state.plant_population_change_count[plant] *= thinning_retention;
        state.plant_standing_dead_population_count[plant] *= thinning_retention;
    } else {
        const cell = plant / state.species_count;
        const population_count = reseed_population_per_m2 * context.cell_area_m2_by_cell.?[cell];
        state.plant_population_per_m2[plant] = reseed_population_per_m2;
        state.plant_population_count[plant] = population_count;
        state.plant_population_change_count[plant] = population_count;
        state.plant_standing_dead_population_count[plant] = population_count;
        try retainReseedProductsInSeedStorage(context, plant);
    }
    if (pruning) {
        const commit = pruning_commit.?;
        commit.structure.initial_clumping_factor[plant] = commit.initial;
        commit.structure.effective_clumping_factor[plant] = commit.effective;
    }
    if (event.termination != .retain) {
        if (context.root_state) |roots| {
            if (plant >= roots.roots_dead.len) return error.PlantHarvestIndexOutOfBounds;
            roots.roots_dead[plant] = true;
        }
        try phenology.terminatePlantBranches(context.branch_development, branches.first, branches.end, true);
        if (context.plant_phenology) |plant_state| {
            if (plant >= plant_state.active.len) return error.PlantHarvestIndexOutOfBounds;
            if (!reseed) plant_state.active[plant] = false;
        }
        if (context.emerged_by_plant) |emerged| {
            if (plant >= emerged.len) return error.PlantHarvestIndexOutOfBounds;
            emerged[plant] = false;
        }
        if (context.growth_stages) |growth| {
            const growth_branches = try growth.branchRange(plant);
            for (growth.branches[growth_branches.first..growth_branches.end]) |*branch| branch.dead = true;
        }
        if (reseed) {
            context.plant_phenology.?.reseed_pending[plant] = true;
        }
    }
    try state.validateFinite();
    try context.branch_development.validateFinite();
}

/// GROSUB JHVST=2 retains all physically exported harvested material as seed
/// storage for the next establishment instead of removing it from the ecosystem.
fn retainReseedProductsInSeedStorage(context: *Context, plant: usize) !void {
    const products = &context.products_by_plant[plant];
    var retained: canopy.ElementalMass = .{};
    addMass(&retained, products.nonstructural.ecosystem_export);
    addMass(&retained, products.foliar.ecosystem_export);
    addMass(&retained, products.nonfoliar.ecosystem_export);
    addMass(&retained, products.woody.ecosystem_export);
    addMass(&retained, products.standing_dead_export);
    const state = context.canopy_state;
    const next_carbon_g = state.plant_seed_storage_carbon_g[plant] + retained.carbon_g;
    const next_nitrogen_g = state.plant_seed_storage_nitrogen_g[plant] + retained.nitrogen_g;
    const next_phosphorus_g = state.plant_seed_storage_phosphorus_g[plant] + retained.phosphorus_g;
    inline for (.{ next_carbon_g, next_nitrogen_g, next_phosphorus_g }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidPlantReseedStorage;
    state.plant_seed_storage_carbon_g[plant] = next_carbon_g;
    state.plant_seed_storage_nitrogen_g[plant] = next_nitrogen_g;
    state.plant_seed_storage_phosphorus_g[plant] = next_phosphorus_g;
    products.nonstructural.ecosystem_export = .{};
    products.foliar.ecosystem_export = .{};
    products.nonfoliar.ecosystem_export = .{};
    products.woody.ecosystem_export = .{};
    products.standing_dead_export = .{};
}

fn applyScheduledStandingDeadHarvest(context: *Context, plant: usize, event: management.HarvestEvent) !void {
    const state = context.canopy_state;
    const total_removed_fraction = if (event.thinning_fraction_or_consumption_rate == 0)
        event.harvested_fraction.standing_dead
    else
        event.thinning_fraction_or_consumption_rate;
    const harvested_fraction = if (event.thinning_fraction_or_consumption_rate == 0)
        event.harvested_fraction.standing_dead
    else if (event.kind == .none)
        event.harvested_fraction.standing_dead * event.thinning_fraction_or_consumption_rate
    else
        event.thinning_fraction_or_consumption_rate;
    inline for (.{ total_removed_fraction, harvested_fraction, event.ecosystem_export_fraction.standing_dead }) |value|
        if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidStandingDeadHarvestFraction;
    if (total_removed_fraction == 0) return;
    const root_fractions = context.root_woody_fraction_by_plant orelse return error.IncompleteStandingDeadHarvestContext;
    const layers = context.canopy_layer_state orelse return error.IncompleteStandingDeadHarvestContext;
    if (plant >= root_fractions.len) return error.PlantHarvestIndexOutOfBounds;
    const nonwoody = root_fractions[plant];
    if (!std.math.isFinite(nonwoody) or nonwoody < 0 or nonwoody > 1) return error.InvalidStandingDeadHarvestFraction;
    const initial: canopy.ElementalMass = .{
        .carbon_g = state.plant_standing_dead_carbon_g[plant],
        .nitrogen_g = state.plant_standing_dead_nitrogen_g[plant],
        .phosphorus_g = state.plant_standing_dead_phosphorus_g[plant],
    };
    const initial_charcoal: canopy.ElementalMass = .{
        .carbon_g = state.plant_charcoal_carbon_g[plant],
        .nitrogen_g = state.plant_charcoal_nitrogen_g[plant],
        .phosphorus_g = state.plant_charcoal_phosphorus_g[plant],
    };
    inline for (.{ initial.carbon_g, initial.nitrogen_g, initial.phosphorus_g }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidStandingDeadHarvestState;
    inline for (.{ initial_charcoal.carbon_g, initial_charcoal.nitrogen_g, initial_charcoal.phosphorus_g }) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidStandingDeadHarvestState;
    // GROSUB IHVST=1 exports only grain; any simultaneous standing-dead
    // removal is returned completely to litter.
    const export_fraction = if (event.kind == .grain)
        0
    else
        harvested_fraction * event.ecosystem_export_fraction.standing_dead;
    const litter_fraction = total_removed_fraction - export_fraction;
    const exported: canopy.ElementalMass = .{
        .carbon_g = initial.carbon_g * export_fraction,
        .nitrogen_g = initial.nitrogen_g * export_fraction,
        .phosphorus_g = initial.phosphorus_g * export_fraction,
    };
    const litter: canopy.ElementalMass = .{
        .carbon_g = initial.carbon_g * litter_fraction,
        .nitrogen_g = initial.nitrogen_g * litter_fraction,
        .phosphorus_g = initial.phosphorus_g * litter_fraction,
    };
    const exported_charcoal: canopy.ElementalMass = .{
        .carbon_g = initial_charcoal.carbon_g * export_fraction,
        .nitrogen_g = initial_charcoal.nitrogen_g * export_fraction,
        .phosphorus_g = initial_charcoal.phosphorus_g * export_fraction,
    };
    const litter_charcoal: canopy.ElementalMass = .{
        .carbon_g = initial_charcoal.carbon_g * litter_fraction,
        .nitrogen_g = initial_charcoal.nitrogen_g * litter_fraction,
        .phosphorus_g = initial_charcoal.phosphorus_g * litter_fraction,
    };
    addMass(&context.products_by_plant[plant].standing_dead_export, exported);
    addMass(&context.products_by_plant[plant].standing_dead_export, exported_charcoal);
    addMass(&context.products_by_plant[plant].standing_dead_charcoal_litter, litter_charcoal);
    addScaledMass(&context.products_by_plant[plant].woody.litter, litter, 1 - nonwoody);
    addScaledMass(&context.products_by_plant[plant].nonfoliar.litter, litter, nonwoody);
    const retained = 1 - total_removed_fraction;
    inline for (.{
        "plant_standing_dead_carbon_g",
        "plant_standing_dead_nitrogen_g",
        "plant_standing_dead_phosphorus_g",
    }) |field_name| @field(state, field_name)[plant] *= retained;
    inline for (.{
        "plant_charcoal_carbon_g",
        "plant_charcoal_nitrogen_g",
        "plant_charcoal_phosphorus_g",
    }) |field_name| @field(state, field_name)[plant] *= retained;
    const kinetic_first = plant * 4;
    inline for (.{
        "plant_standing_dead_carbon_by_kinetic_g",
        "plant_standing_dead_nitrogen_by_kinetic_g",
        "plant_standing_dead_phosphorus_by_kinetic_g",
    }) |field_name| {
        for (@field(state, field_name)[kinetic_first..][0..4]) |*value| value.* *= retained;
    }
    const cell = plant / state.species_count;
    for (0..layers.layer_count) |layer| {
        const plant_layer = plant * layers.layer_count + layer;
        const removed_area_m2 = layers.plant_standing_dead_area_m2[plant_layer] * total_removed_fraction;
        layers.plant_standing_dead_area_m2[plant_layer] *= retained;
        layers.cell_standing_dead_area_m2[cell * layers.layer_count + layer] =
            @max(0, layers.cell_standing_dead_area_m2[cell * layers.layer_count + layer] - removed_area_m2);
        const projected_first = plant_layer * layers.inclination_count;
        for (layers.plant_standing_dead_projected_surface_m2[projected_first..][0..layers.inclination_count]) |*area_m2|
            area_m2.* *= retained;
    }
}

fn totalBranchPool(state: *const canopy.State, branches: canopy.Range, comptime field_name: []const u8) f64 {
    var total: f64 = 0;
    for (@field(state, field_name)[branches.first..branches.end]) |value| total += @max(0, value);
    return total;
}

fn removalFraction(target: f64, total: f64) f64 {
    return if (total > 0) std.math.clamp(target / total, 0, 1) else 0;
}

fn routeGrazedMass(products: *canopy.HarvestProducts, removed: canopy.ElementalMass, export_fraction: f64) void {
    addScaledMass(&products.ecosystem_export, removed, export_fraction);
    addScaledMass(&products.litter, removed, 1 - export_fraction);
}

fn productLedgerCarbonG(ledger: ProductLedger) f64 {
    // harvested_grain is a diagnostic subset of reproductive products, not
    // a second physical pool.
    var total_g_c = ledger.standing_dead_export.carbon_g +
        ledger.standing_dead_charcoal_litter.carbon_g;
    inline for (.{ "nonstructural", "foliar", "nonfoliar", "woody" }) |field_name| {
        total_g_c += @field(ledger, field_name).ecosystem_export.carbon_g;
        total_g_c += @field(ledger, field_name).litter.carbon_g;
    }
    for (ledger.direct_litter.woody_carbon_g, ledger.direct_litter.nonwoody_carbon_g) |woody, nonwoody|
        total_g_c += woody + nonwoody;
    for (ledger.manure.organic_by_biochemical_fraction) |mass| total_g_c += mass.carbon_g;
    return total_g_c;
}

fn validateGrazingFraction(value: f64) !void {
    if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidGrazingFraction;
}

fn validateNonnegativeFinite(comptime field_name: []const u8, values: []const f64, first: usize, end: usize) !void {
    if (first > end or end > values.len) return error.GrazingStateDimensionMismatch;
    for (values[first..end], first..) |value, index| {
        if (!std.math.isFinite(value) or value < 0) {
            if (!builtin.is_test) std.log.err("invalid grazing state: field={s} index={d} value={e}", .{ field_name, index, value });
            return error.InvalidGrazingState;
        }
    }
}

fn preflightGrazing(context: *const Context, plant: usize, event: management.HarvestEvent, layers: *const canopy_layers.State) !void {
    const state = context.canopy_state;
    if (!std.math.isFinite(event.cutting_height_m_or_lai_fraction) or event.cutting_height_m_or_lai_fraction < 0 or
        !std.math.isFinite(event.thinning_fraction_or_consumption_rate) or event.thinning_fraction_or_consumption_rate < 0)
        return error.InvalidGrazingEvent;
    inline for (@typeInfo(management.RemovalFractions).@"struct".fields) |field| {
        try validateGrazingFraction(@field(event.harvested_fraction, field.name));
        try validateGrazingFraction(@field(event.ecosystem_export_fraction, field.name));
    }
    try state.validateFinite();
    inline for (.{
        "plant_standing_dead_carbon_g",
        "plant_standing_dead_nitrogen_g",
        "plant_standing_dead_phosphorus_g",
    }) |field_name| try validateNonnegativeFinite(field_name, @field(state, field_name), plant, plant + 1);
    const kinetic_first = try std.math.mul(usize, plant, 4);
    inline for (.{
        "plant_standing_dead_carbon_by_kinetic_g",
        "plant_standing_dead_nitrogen_by_kinetic_g",
        "plant_standing_dead_phosphorus_by_kinetic_g",
    }) |field_name| try validateNonnegativeFinite(field_name, @field(state, field_name), kinetic_first, kinetic_first + 4);
    const standing_area_first = try std.math.mul(usize, plant, layers.layer_count);
    try validateNonnegativeFinite("plant_standing_dead_area_m2", layers.plant_standing_dead_area_m2, standing_area_first, standing_area_first + layers.layer_count);
    const branches = try state.branchRange(plant);
    inline for (.{
        "branch_leaf_carbon_g",                "branch_leaf_nitrogen_g",                "branch_leaf_phosphorus_g",
        "branch_sheath_carbon_g",              "branch_sheath_nitrogen_g",              "branch_sheath_phosphorus_g",
        "branch_husk_carbon_g",                "branch_husk_nitrogen_g",                "branch_husk_phosphorus_g",
        "branch_ear_carbon_g",                 "branch_ear_nitrogen_g",                 "branch_ear_phosphorus_g",
        "branch_grain_carbon_g",               "branch_grain_nitrogen_g",               "branch_grain_phosphorus_g",
        "branch_stalk_carbon_g",               "branch_stalk_nitrogen_g",               "branch_stalk_phosphorus_g",
        "branch_reserve_carbon_g",             "branch_reserve_nitrogen_g",             "branch_reserve_phosphorus_g",
        "branch_mobile_carbon_g",              "branch_mobile_nitrogen_g",              "branch_mobile_phosphorus_g",
        "branch_symbiont_mobile_carbon_g",     "branch_symbiont_mobile_nitrogen_g",     "branch_symbiont_mobile_phosphorus_g",
        "branch_symbiont_structural_carbon_g", "branch_symbiont_structural_nitrogen_g", "branch_symbiont_structural_phosphorus_g",
    }) |field_name| try validateNonnegativeFinite(field_name, @field(state, field_name), branches.first, branches.end);
    const expected_samples_per_node = try std.math.mul(usize, layers.layer_count, try std.math.mul(usize, layers.inclination_count, layers.azimuth_count));
    for (branches.first..branches.end) |branch| {
        const nodes = try state.nodeRange(branch);
        inline for (.{
            "node_leaf_area_m2",         "node_leaf_carbon_g",          "node_leaf_nitrogen_g",     "node_leaf_phosphorus_g",
            "node_sheath_carbon_g",      "node_sheath_nitrogen_g",      "node_sheath_phosphorus_g", "node_internode_carbon_g",
            "node_internode_nitrogen_g", "node_internode_phosphorus_g",
        }) |field_name| try validateNonnegativeFinite(field_name, @field(state, field_name), nodes.first, nodes.end);
        for (nodes.first..nodes.end) |node| {
            const samples = try state.sampleRange(node);
            if (samples.end - samples.first != expected_samples_per_node) return error.GrazingSampleTopologyMismatch;
            inline for (.{
                "sample_leaf_area_m2",    "sample_exposed_leaf_area_m2", "sample_leaf_carbon_g",
                "sample_leaf_nitrogen_g", "sample_leaf_phosphorus_g",    "sample_stalk_area_m2",
            }) |field_name| try validateNonnegativeFinite(field_name, @field(state, field_name), samples.first, samples.end);
            const layer_first = try std.math.mul(usize, node, layers.layer_count);
            inline for (.{ "node_leaf_area_m2", "node_leaf_carbon_g", "node_leaf_nitrogen_g", "node_leaf_phosphorus_g" }) |field_name|
                try validateNonnegativeFinite(field_name, @field(layers, field_name), layer_first, layer_first + layers.layer_count);
        }
    }
}

/// Live GROSUB animal/insect grazing transaction for one plant and hour.
pub fn applyGrazingEvent(
    context: *Context,
    plant: usize,
    event: management.HarvestEvent,
    landscape_average_shoot_carbon_g_c: f64,
    horizontal_cell_area_m2: f64,
) !f64 {
    if (event.kind != .animal_grazing and event.kind != .insect_grazing) return error.NotGrazingEvent;
    const state = context.canopy_state;
    if (plant >= state.plant_total_shoot_carbon_g.len or plant >= context.products_by_plant.len or
        !std.math.isFinite(landscape_average_shoot_carbon_g_c) or landscape_average_shoot_carbon_g_c < 0)
        return error.InvalidGrazingRuntimeInput;
    const layers = context.canopy_layer_state orelse return error.IncompleteGrazingContext;
    const branches = try state.branchRange(plant);
    const science = context.science_by_plant[plant];
    try validateScience(science);
    try preflightGrazing(context, plant, event, layers);
    const initial_product_carbon_g_c = productLedgerCarbonG(context.products_by_plant[plant]);
    var stalk_area_m2: f64 = 0;
    for (branches.first..branches.end) |branch| {
        const first = branch * layers.layer_count;
        for (layers.branch_stalk_area_m2[first..][0..layers.layer_count]) |area| stalk_area_m2 += area;
    }
    const leaf_area_m2 = try layers.plantLeafAreaM2(state, plant);
    const demand_g_c = try canopy.sourceOrderGrazingCarbonDemandGPerH(
        event.kind == .animal_grazing,
        event.cutting_height_m_or_lai_fraction,
        event.thinning_fraction_or_consumption_rate,
        horizontal_cell_area_m2,
        leaf_area_m2 + stalk_area_m2,
        state.plant_uptake_growth_temperature_response[plant],
        state.plant_total_shoot_carbon_g[plant],
        landscape_average_shoot_carbon_g_c,
        context.plant_structural_presence_threshold_g_per_plant *
            state.plant_population_count[plant],
    );
    const pools: canopy.GrazingPools = .{
        .leaf_carbon_g = totalBranchPool(state, branches, "branch_leaf_carbon_g"),
        .sheath_carbon_g = totalBranchPool(state, branches, "branch_sheath_carbon_g"),
        .husk_carbon_g = totalBranchPool(state, branches, "branch_husk_carbon_g"),
        .ear_carbon_g = totalBranchPool(state, branches, "branch_ear_carbon_g"),
        .grain_carbon_g = totalBranchPool(state, branches, "branch_grain_carbon_g"),
        .stalk_carbon_g = totalBranchPool(state, branches, "branch_stalk_carbon_g"),
        .reserve_carbon_g = totalBranchPool(state, branches, "branch_reserve_carbon_g"),
    };
    const allocation = try canopy.allocateGrazingDemand(
        demand_g_c,
        event.harvested_fraction.leaf,
        event.harvested_fraction.nonfoliar,
        event.harvested_fraction.woody,
        state.plant_mobile_carbon_concentration_g_per_g[plant],
        state.plant_symbiont_mobile_carbon_concentration_g_per_g[plant],
        pools,
    );

    // GROSUB 9865/9855/9845: canopy layers top-to-bottom, branches in
    // source order, and nodes newest-to-oldest. Each branch-layer receives
    // its share of the plant leaf demand from its pre-removal carbon.
    const angular_count = try std.math.mul(usize, layers.inclination_count, layers.azimuth_count);
    var layer_cursor = layers.layer_count;
    while (layer_cursor > 0) {
        layer_cursor -= 1;
        const layer = layer_cursor;
        for (branches.first..branches.end) |branch| {
            const nodes = try state.nodeRange(branch);
            var branch_layer_carbon_g_c: f64 = 0;
            for (nodes.first..nodes.end) |node|
                branch_layer_carbon_g_c += @max(0, layers.node_leaf_carbon_g[node * layers.layer_count + layer]);
            var branch_layer_demand_g_c = try canopy.sourceOrderBranchLayerLeafDemand(
                pools.leaf_carbon_g,
                allocation.structural_leaf_carbon_g,
                branch_layer_carbon_g_c,
                context.plant_structural_presence_threshold_g_per_plant *
                    state.plant_population_count[plant],
            );
            var node_cursor = nodes.end;
            while (node_cursor > nodes.first and branch_layer_demand_g_c > 0) {
                node_cursor -= 1;
                const node = node_cursor;
                const node_layer_carbon_g_c = @max(0, layers.node_leaf_carbon_g[node * layers.layer_count + layer]);
                if (node_layer_carbon_g_c <= 0) continue;
                const removed_carbon_g_c = @min(branch_layer_demand_g_c, node_layer_carbon_g_c);
                const leaf_remaining = std.math.clamp(1 - removed_carbon_g_c / node_layer_carbon_g_c, 0, 1);
                for (0..angular_count) |angular| {
                    const sample = layer * angular_count + angular;
                    const retention: canopy.LayerHarvestRetention = .{
                        .remaining_fraction = leaf_remaining,
                        .unexported_fraction = leaf_remaining + (1 - leaf_remaining) * (1 - event.ecosystem_export_fraction.leaf),
                        .height_below_cut_fraction = leaf_remaining,
                    };
                    const removed = try canopy.harvestLeafLayerSample(
                        state,
                        branch,
                        node - nodes.first,
                        sample,
                        retention,
                        science.carbon_woody_fraction,
                        science.leaf_nitrogen_woody_fraction,
                        science.leaf_phosphorus_woody_fraction,
                        node - nodes.first == 1,
                    );
                    addProducts(&context.products_by_plant[plant].foliar, removed.foliar);
                    addProducts(&context.products_by_plant[plant].woody, removed.woody);
                }
                branch_layer_demand_g_c = @max(0, branch_layer_demand_g_c - removed_carbon_g_c);
            }
        }
    }

    for (branches.first..branches.end) |branch| {
        const nodes = try state.nodeRange(branch);
        var initial_branch_leaf_carbon_g_c: f64 = 0;
        for (nodes.first..nodes.end) |node| {
            const first = node * layers.layer_count;
            for (layers.node_leaf_carbon_g[first..][0..layers.layer_count]) |carbon_g_c|
                initial_branch_leaf_carbon_g_c += @max(0, carbon_g_c);
        }
        const initial_branch_sheath_carbon_g_c = state.branch_sheath_carbon_g[branch];
        var branch_sheath_demand_g_c = if (pools.sheath_carbon_g > 0)
            allocation.structural_sheath_carbon_g * initial_branch_sheath_carbon_g_c / pools.sheath_carbon_g
        else
            0;
        var node_cursor = nodes.end;
        while (node_cursor > nodes.first and branch_sheath_demand_g_c > 0) {
            node_cursor -= 1;
            const node = node_cursor;
            const initial_sheath_carbon_g_c = state.node_sheath_carbon_g[node];
            if (initial_sheath_carbon_g_c <= 0) continue;
            const removed_carbon_g_c = @min(branch_sheath_demand_g_c, initial_sheath_carbon_g_c);
            const sheath_remaining = std.math.clamp(1 - removed_carbon_g_c / initial_sheath_carbon_g_c, 0, 1);
            const removed = try canopy.harvestNodeSheath(
                state,
                branch,
                node - nodes.first,
                sheath_remaining,
                sheath_remaining + (1 - sheath_remaining) * (1 - event.ecosystem_export_fraction.nonfoliar),
                science.carbon_woody_fraction,
                science.sheath_nitrogen_woody_fraction,
                science.sheath_phosphorus_woody_fraction,
                false,
                0,
            );
            addProducts(&context.products_by_plant[plant].nonfoliar, removed.nonwoody);
            addProducts(&context.products_by_plant[plant].woody, removed.woody);
            branch_sheath_demand_g_c = @max(0, branch_sheath_demand_g_c - removed_carbon_g_c);
        }

        const total_leaf_sheath_carbon_g_c = pools.leaf_carbon_g + pools.sheath_carbon_g;
        const branch_share = if (total_leaf_sheath_carbon_g_c >
            context.plant_tissue_presence_threshold_g_per_plant *
                state.plant_population_count[plant])
            @max(0, initial_branch_leaf_carbon_g_c + initial_branch_sheath_carbon_g_c) / total_leaf_sheath_carbon_g_c
        else
            0;
        const branch_mobile_target_g_c = allocation.mobile_carbon_g * branch_share;
        const branch_mobile_carbon_g_c = state.branch_mobile_carbon_g[branch];
        const host_mobile = try canopy.sourceOrderProportionalMobileRemoval(
            .{
                .carbon_g = branch_mobile_carbon_g_c,
                .nitrogen_g = state.branch_mobile_nitrogen_g[branch],
                .phosphorus_g = state.branch_mobile_phosphorus_g[branch],
            },
            branch_mobile_target_g_c,
            context.plant_structural_presence_threshold_g_per_plant *
                state.plant_population_count[plant],
        );
        const mobile_remaining = if (branch_mobile_carbon_g_c > 0)
            host_mobile.remaining.carbon_g / branch_mobile_carbon_g_c
        else
            0;
        const is_c4 = if (context.canopy_biochemistry_parameters_by_plant) |parameters|
            parameters[plant].pathway == .c4
        else
            true;
        const intermediate_remaining = try canopy.sourceOrderC4IntermediateRetention(
            is_c4,
            branch_mobile_carbon_g_c,
            host_mobile.remaining.carbon_g,
            context.plant_structural_presence_threshold_g_per_plant *
                state.plant_population_count[plant],
        );
        const mobile_removed = try canopy.harvestBranchMobilePoolsWithIntermediateRetention(
            state,
            branch,
            mobile_remaining,
            intermediate_remaining,
        );
        routeGrazedMass(&context.products_by_plant[plant].nonfoliar, mobile_removed, event.ecosystem_export_fraction.nonfoliar);

        const branch_symbiont_target_g_c = allocation.symbiont_mobile_carbon_g * branch_share;
        const branch_symbiont_mobile_carbon_g_c = state.branch_symbiont_mobile_carbon_g[branch];
        const symbiont_mobile = try canopy.sourceOrderProportionalMobileRemoval(
            .{
                .carbon_g = branch_symbiont_mobile_carbon_g_c,
                .nitrogen_g = state.branch_symbiont_mobile_nitrogen_g[branch],
                .phosphorus_g = state.branch_symbiont_mobile_phosphorus_g[branch],
            },
            branch_symbiont_target_g_c,
            context.plant_structural_presence_threshold_g_per_plant *
                state.plant_population_count[plant],
        );
        const symbiont_remaining = if (branch_symbiont_mobile_carbon_g_c > 0)
            std.math.clamp(
                symbiont_mobile.remaining.carbon_g / branch_symbiont_mobile_carbon_g_c,
                0,
                1,
            )
        else
            0;
        var symbiont_removed: canopy.ElementalMass = .{};
        inline for (
            .{
                "branch_symbiont_mobile_carbon_g",     "branch_symbiont_mobile_nitrogen_g",     "branch_symbiont_mobile_phosphorus_g",
                "branch_symbiont_structural_carbon_g", "branch_symbiont_structural_nitrogen_g", "branch_symbiont_structural_phosphorus_g",
            },
            .{ "carbon_g", "nitrogen_g", "phosphorus_g", "carbon_g", "nitrogen_g", "phosphorus_g" },
        ) |state_field, mass_field| {
            const initial = @field(state, state_field)[branch];
            @field(symbiont_removed, mass_field) += initial * (1 - symbiont_remaining);
            @field(state, state_field)[branch] = initial * symbiont_remaining;
        }
        routeGrazedMass(&context.products_by_plant[plant].nonstructural, symbiont_removed, event.ecosystem_export_fraction.nonfoliar);
    }

    const reproductive_retention = try canopy.sourceOrderReproductiveRetention(.{
        .grazing = true,
        .reproductive_organs_reached_by_cut = false,
        .grain_or_pruning = false,
        .thinning_fraction = 0,
        .harvested_nonfoliar_fraction = 0,
        .total_husk_carbon_g_c = pools.husk_carbon_g,
        .total_ear_carbon_g_c = pools.ear_carbon_g,
        .total_grain_carbon_g_c = pools.grain_carbon_g,
        .grazed_husk_carbon_g_c = allocation.husk_carbon_g,
        .grazed_ear_carbon_g_c = allocation.ear_carbon_g,
        .grazed_grain_carbon_g_c = allocation.grain_carbon_g,
        .plant_presence_threshold_g_c = context.plant_structural_presence_threshold_g_per_plant *
            state.plant_population_count[plant],
    });
    for (branches.first..branches.end) |branch| {
        const reproductive = try canopy.harvestReproductiveOrgans(
            state,
            branch,
            reproductive_retention,
        );
        // products already contains husk, ear, and grain removal;
        // harvested_grain is a diagnostic subset and must not be added twice.
        routeGrazedMass(&context.products_by_plant[plant].nonfoliar, reproductive.products.ecosystem_export, event.ecosystem_export_fraction.nonfoliar);
    }

    const plant_presence_threshold_g_c =
        context.plant_tissue_presence_threshold_g_per_plant *
        state.plant_population_count[plant];
    const stalk_remaining = if (pools.stalk_carbon_g > plant_presence_threshold_g_c)
        1 - removalFraction(allocation.stalk_carbon_g, pools.stalk_carbon_g)
    else
        1;
    for (branches.first..branches.end) |branch| {
        const nodes = try state.nodeRange(branch);
        for (nodes.first..nodes.end) |node|
            try canopy.commitInternodeHarvest(state, branch, node - nodes.first, stalk_remaining, false, 0);
        // GROSUB applies WHVRVH independently to each branch reserve rather
        // than distributing it by the plant-total reserve pool.
        const branch_reserve_carbon_g_c = state.branch_reserve_carbon_g[branch];
        const reserve_retention = try canopy.sourceOrderStalkReserveRetention(
            true,
            state.branch_stalk_carbon_g[branch] * stalk_remaining,
            .{
                .remaining_fraction = stalk_remaining,
                .unexported_fraction = stalk_remaining,
                .height_below_cut_fraction = 0,
            },
            branch_reserve_carbon_g_c,
            allocation.reserve_carbon_g,
            context.plant_structural_presence_threshold_g_per_plant *
                state.plant_population_count[plant],
        );
        const reserve_remaining = reserve_retention.remaining_fraction;
        const removed = try canopy.harvestBranchStalkAndReserve(
            state,
            branch,
            stalk_remaining,
            stalk_remaining + (1 - stalk_remaining) * (1 - event.ecosystem_export_fraction.woody),
            reserve_remaining,
            reserve_remaining + (1 - reserve_remaining) * (1 - event.ecosystem_export_fraction.woody),
        );
        addProducts(&context.products_by_plant[plant].woody, removed);
    }

    var returned_mass: canopy.ElementalMass = .{};
    inline for (.{ "foliar", "nonfoliar", "woody" }) |field_name| {
        addMass(&returned_mass, @field(context.products_by_plant[plant], field_name).litter);
        @field(context.products_by_plant[plant], field_name).litter = .{};
    }

    var standing_dead_area_m2: f64 = 0;
    const standing_area_first = plant * layers.layer_count;
    for (layers.plant_standing_dead_area_m2[standing_area_first..][0..layers.layer_count]) |area_m2|
        standing_dead_area_m2 += area_m2;
    const standing_dead_demand_g_c = try grazing_manure.standingDeadDemandGPerH(
        event.kind,
        event.cutting_height_m_or_lai_fraction,
        event.thinning_fraction_or_consumption_rate,
        horizontal_cell_area_m2,
        standing_dead_area_m2,
        event.harvested_fraction.standing_dead,
    );
    const standing_dead_carbon_g_c =
        state.plant_standing_dead_carbon_g[plant] +
        state.plant_charcoal_carbon_g[plant];
    const standing_dead_remaining = if (standing_dead_carbon_g_c > 0)
        std.math.clamp(1 - standing_dead_demand_g_c / standing_dead_carbon_g_c, 0, 1)
    else
        1;
    const removed_standing_dead: canopy.ElementalMass = .{
        .carbon_g = state.plant_standing_dead_carbon_g[plant] * (1 - standing_dead_remaining),
        .nitrogen_g = state.plant_standing_dead_nitrogen_g[plant] * (1 - standing_dead_remaining),
        .phosphorus_g = state.plant_standing_dead_phosphorus_g[plant] * (1 - standing_dead_remaining),
    };
    const removed_standing_dead_charcoal: canopy.ElementalMass = .{
        .carbon_g = state.plant_charcoal_carbon_g[plant] * (1 - standing_dead_remaining),
        .nitrogen_g = state.plant_charcoal_nitrogen_g[plant] * (1 - standing_dead_remaining),
        .phosphorus_g = state.plant_charcoal_phosphorus_g[plant] * (1 - standing_dead_remaining),
    };
    inline for (.{
        "plant_standing_dead_carbon_g",
        "plant_standing_dead_nitrogen_g",
        "plant_standing_dead_phosphorus_g",
    }) |field_name| @field(state, field_name)[plant] *= standing_dead_remaining;
    inline for (.{
        "plant_charcoal_carbon_g",
        "plant_charcoal_nitrogen_g",
        "plant_charcoal_phosphorus_g",
    }) |field_name| @field(state, field_name)[plant] *= standing_dead_remaining;
    const kinetic_first = plant * 4;
    inline for (.{
        "plant_standing_dead_carbon_by_kinetic_g",
        "plant_standing_dead_nitrogen_by_kinetic_g",
        "plant_standing_dead_phosphorus_by_kinetic_g",
    }) |field_name| {
        for (@field(state, field_name)[kinetic_first..][0..4]) |*value| value.* *= standing_dead_remaining;
    }
    const cell = plant / state.species_count;
    for (0..layers.layer_count) |layer| {
        const plant_layer = standing_area_first + layer;
        const removed_area_m2 = layers.plant_standing_dead_area_m2[plant_layer] * (1 - standing_dead_remaining);
        layers.plant_standing_dead_area_m2[plant_layer] *= standing_dead_remaining;
        layers.cell_standing_dead_area_m2[cell * layers.layer_count + layer] =
            @max(0, layers.cell_standing_dead_area_m2[cell * layers.layer_count + layer] - removed_area_m2);
        const projected_first = plant_layer * layers.inclination_count;
        for (layers.plant_standing_dead_projected_surface_m2[projected_first..][0..layers.inclination_count]) |*area_m2|
            area_m2.* *= standing_dead_remaining;
    }
    addScaledMass(&context.products_by_plant[plant].standing_dead_export, removed_standing_dead, event.ecosystem_export_fraction.standing_dead);
    addScaledMass(&context.products_by_plant[plant].standing_dead_export, removed_standing_dead_charcoal, event.ecosystem_export_fraction.standing_dead);
    addScaledMass(&returned_mass, removed_standing_dead, 1 - event.ecosystem_export_fraction.standing_dead);
    addScaledMass(&returned_mass, removed_standing_dead_charcoal, 1 - event.ecosystem_export_fraction.standing_dead);
    if (returned_mass.carbon_g > 0 or returned_mass.nitrogen_g > 0 or returned_mass.phosphorus_g > 0)
        try grazing_manure.add(&context.products_by_plant[plant].manure, try grazing_manure.partition(event.kind, returned_mass));

    try state.validateFinite();
    const removed_product_carbon_g_c = productLedgerCarbonG(context.products_by_plant[plant]) - initial_product_carbon_g_c;
    if (!std.math.isFinite(removed_product_carbon_g_c) or removed_product_carbon_g_c < -1e-10)
        return error.InvalidGrazingProductBalance;
    return @max(0, removed_product_carbon_g_c);
}

test "GROSUB forest self thinning retains exact density law and runtime units" {
    try std.testing.expectEqual(@as(f64, 0), try forestSelfThinningFraction(0, 10));
    try std.testing.expectEqual(@as(f64, 0), try forestSelfThinningFraction(0.25, 0.05));
    // At 0.25 m PPQ is exactly 0.1 plants m-2.
    try std.testing.expectApproxEqAbs(@as(f64, 0.09), try forestSelfThinningFraction(0.25, 1), 1.0e-15);
    try std.testing.expectError(error.InvalidForestSelfThinningInput, forestSelfThinningFraction(-0.1, 1));
}

test "GROSUB first substep resets only hourly disturbance products" {
    const current: HourlyDisturbanceReset = .{
        .previous_cumulative_harvest_carbon_g_c = 3,
        .manure_organic_carbon_g_c = .{ 1, 2, 3, 4 },
        .manure_organic_nitrogen_g_n = .{ 5, 6, 7, 8 },
        .manure_organic_phosphorus_g_p = .{ 9, 10, 11, 12 },
        .manure_inorganic_nitrogen_g_n = 13,
        .manure_inorganic_phosphorus_g_p = 14,
    };
    const unchanged = try sourceOrderHourlyDisturbanceReset(false, 20, current);
    try std.testing.expectEqualDeep(current, unchanged);
    const reset = try sourceOrderHourlyDisturbanceReset(true, 20, current);
    try std.testing.expectEqual(@as(f64, 20), reset.previous_cumulative_harvest_carbon_g_c);
    try std.testing.expectEqual([4]f64{ 0, 0, 0, 0 }, reset.manure_organic_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), reset.manure_inorganic_nitrogen_g_n);
}

test "GROSUB forest self thinning selector preserves monthly noon and event gates" {
    try std.testing.expect(try sourceOrderForestSelfThinningIsEnabled(30, 12, 12.75, 1, 2, -1));
    try std.testing.expect(try sourceOrderForestSelfThinningIsEnabled(360, 12, 12.75, 1, 2, 4));
    try std.testing.expect(try sourceOrderForestSelfThinningIsEnabled(30, 12, 12.75, 1, 2, 6));
    try std.testing.expect(!try sourceOrderForestSelfThinningIsEnabled(30, 12, 12.75, 1, 2, 0));
    try std.testing.expect(!try sourceOrderForestSelfThinningIsEnabled(29, 12, 12.75, 1, 2, -1));
    try std.testing.expect(!try sourceOrderForestSelfThinningIsEnabled(30, 13, 12.75, 1, 2, -1));
    try std.testing.expect(!try sourceOrderForestSelfThinningIsEnabled(30, 12, 12.75, 0, 2, -1));
    try std.testing.expect(!try sourceOrderForestSelfThinningIsEnabled(30, 12, 12.75, 1, 1, -1));
}

test "GROSUB pruning multiplies the persistent canopy clumping factor" {
    try std.testing.expectApproxEqAbs(@as(f64, 0.48), try prunedClumpingFactor(0.8, 0.6), 1.0e-15);
    try std.testing.expectError(error.InvalidPruningClumpingFraction, prunedClumpingFactor(0.8, -0.1));
}

test "GROSUB negative HVST interpolates combined canopy leaf area across runtime layers" {
    try std.testing.expectApproxEqAbs(
        @as(f64, 1.5),
        try cuttingHeightFromLeafAreaRemoval(0.5, &.{ 0, 1, 3 }, &.{ 2, 4 }, 1.0e-12),
        1.0e-15,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        try cuttingHeightFromLeafAreaRemoval(1, &.{ 0, 1, 3 }, &.{ 2, 4 }, 1.0e-12),
    );
    try std.testing.expectError(
        error.InvalidLeafAreaHarvestGeometry,
        cuttingHeightFromLeafAreaRemoval(1.1, &.{ 0, 1 }, &.{1}, 1.0e-12),
    );
}

test "GROSUB cutting height uses authoritative combined canopy leaf area" {
    const exact = try sourceOrderCuttingHeightFromLeafAreaRemoval(
        0.5,
        8,
        &.{ 0, 1, 3 },
        &.{ 2, 4 },
        1.0e-12,
    );
    try std.testing.expectEqual(@as(f64, 2), exact);
    const recomputed = try cuttingHeightFromLeafAreaRemoval(0.5, &.{ 0, 1, 3 }, &.{ 2, 4 }, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 1.5), recomputed);
}

test "GROSUB aboveground disturbance dispatch and population update are exact" {
    try std.testing.expect(try sourceOrderAbovegroundDisturbanceIsEnabled(4, 3, 12.5));
    try std.testing.expect(try sourceOrderAbovegroundDisturbanceIsEnabled(6, 3, 12.5));
    try std.testing.expect(try sourceOrderAbovegroundDisturbanceIsEnabled(2, 12, 12.5));
    try std.testing.expect(!try sourceOrderAbovegroundDisturbanceIsEnabled(2, 11, 12.5));
    try std.testing.expect(!try sourceOrderAbovegroundDisturbanceIsEnabled(-1, 12, 12.5));

    const current: PopulationAfterDisturbance = .{
        .living_population_per_m2 = 4,
        .living_population_count = 40,
        .standing_dead_population_count = 10,
    };
    const thinned = try sourceOrderPopulationAfterDisturbance(false, 0.25, current, 7, 10);
    try std.testing.expectEqual(@as(f64, 3), thinned.living_population_per_m2);
    try std.testing.expectEqual(@as(f64, 30), thinned.living_population_count);
    try std.testing.expectEqual(@as(f64, 7.5), thinned.standing_dead_population_count);
    const reseeded = try sourceOrderPopulationAfterDisturbance(true, 0.25, current, 7, 10);
    try std.testing.expectEqual(@as(f64, 7), reseeded.living_population_per_m2);
    try std.testing.expectEqual(@as(f64, 70), reseeded.living_population_count);
    try std.testing.expectEqual(@as(f64, 70), reseeded.standing_dead_population_count);
}

test "GROSUB harvest litter uses organ-specific runtime kinetics conservatively" {
    const foliar: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 0, 1, 0, 0 },
        .phosphorus = .{ 0, 0, 1, 0 },
    };
    const nonfoliar: litter_partition.ElementFractions = .{
        .carbon = .{ 0, 1, 0, 0 },
        .nitrogen = .{ 0, 0, 1, 0 },
        .phosphorus = .{ 0, 0, 0, 1 },
    };
    const woody: litter_partition.ElementFractions = .{
        .carbon = .{ 0, 0, 1, 0 },
        .nitrogen = .{ 0, 0, 0, 1 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const products: ProductLedger = .{
        .foliar = .{ .litter = .{ .carbon_g = 2, .nitrogen_g = 0.2, .phosphorus_g = 0.02 } },
        .nonfoliar = .{ .litter = .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 } },
        .woody = .{ .litter = .{ .carbon_g = 5, .nitrogen_g = 0.5, .phosphorus_g = 0.05 } },
    };
    const result = try harvestLitterToKinetics(products, foliar, foliar, nonfoliar, woody);
    try std.testing.expectEqual([4]f64{ 2, 3, 0, 0 }, result.nonwoody_carbon_g);
    try std.testing.expectEqual([4]f64{ 0, 0, 5, 0 }, result.woody_carbon_g);
    var carbon_g_c: f64 = 0;
    var nitrogen_g_n: f64 = 0;
    var phosphorus_g_p: f64 = 0;
    for (0..4) |kinetic| {
        carbon_g_c += result.nonwoody_carbon_g[kinetic] + result.woody_carbon_g[kinetic];
        nitrogen_g_n += result.nonwoody_nitrogen_g[kinetic] + result.woody_nitrogen_g[kinetic];
        phosphorus_g_p += result.nonwoody_phosphorus_g[kinetic] + result.woody_phosphorus_g[kinetic];
    }
    try std.testing.expectApproxEqAbs(@as(f64, 10), carbon_g_c, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 1), nitrogen_g_n, 1.0e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), phosphorus_g_p, 1.0e-15);
}

fn applyRootSymbiontHarvest(context: *Context, plant: usize, remaining_fraction: f64) !void {
    if (remaining_fraction == 1) return;
    const roots = context.root_state orelse return;
    const partitions = context.root_litter_partition orelse return error.IncompleteRootHarvestContext;
    const organic = context.soil_organic_state orelse return error.IncompleteRootHarvestContext;
    const grid = context.grid orelse return error.IncompleteRootHarvestContext;
    if (plant >= roots.plant_count or plant >= partitions.plant_count or plant >= context.science_by_plant.len) return error.PlantHarvestIndexOutOfBounds;
    if (!std.math.isFinite(remaining_fraction) or remaining_fraction < 0 or remaining_fraction > 1) return error.InvalidRootHarvestRetention;
    const cell = plant / context.canopy_state.species_count;
    if (cell >= grid.cell_count) return error.PlantHarvestIndexOutOfBounds;
    const fine = try partitions.get(plant, .fine_root);
    const coarse = try partitions.get(plant, .coarse_wood);
    const mobile = try partitions.get(plant, .nonstructural);
    const nitrogen_fixation_type = context.science_by_plant[plant].nitrogen_fixation_type;
    const woody_fraction = if (context.root_woody_fraction_by_plant) |values| blk: {
        if (plant >= values.len or !std.math.isFinite(values[plant]) or values[plant] < 0 or values[plant] > 1) return error.InvalidRootHarvestWoodyFraction;
        break :blk values[plant];
    } else 0;
    var removed_host_carbon_g_c: f64 = 0;

    // Validate every soil publication before changing any root or soil pool.
    for (0..grid.active_soil_layer_count[cell]) |layer| {
        const root = try roots.layerIndex(plant, 0, layer);
        const result = try noduleHarvestResult(
            nitrogen_fixation_type,
            .{
                .carbon_g_c = roots.symbiont_structural_carbon_g_c[root],
                .nitrogen_g_n = roots.symbiont_structural_nitrogen_g_n[root],
                .phosphorus_g_p = roots.symbiont_structural_phosphorus_g_p[root],
            },
            .{
                .carbon_g_c = roots.symbiont_mobile_carbon_g_c[root],
                .nitrogen_g_n = roots.symbiont_mobile_nitrogen_g_n[root],
                .phosphorus_g_p = roots.symbiont_mobile_phosphorus_g_p[root],
            },
            remaining_fraction,
            fine,
            mobile,
        );
        var publication: root_litterfall.LayerInput = .{};
        try publication.add(result.litterfall);
        const host = try hostLayerHarvest(roots, plant, layer, remaining_fraction, woody_fraction, fine, coarse, mobile);
        try publication.add(host.litterfall);
        removed_host_carbon_g_c += host.removed_carbon_g_c;
        if (context.root_litter_carbon_ledger) |ledger| {
            const host_domain_zero = try hostLayerHarvestDomain(roots, plant, 0, layer, remaining_fraction, woody_fraction, fine, coarse, mobile);
            const host_domain_one = try hostLayerHarvestDomain(roots, plant, 1, layer, remaining_fraction, woody_fraction, fine, coarse, mobile);
            try ledger.validateCarbonAdd(
                plant,
                0,
                layer,
                try root_litter_ledger.totalCarbon(host_domain_zero.litterfall) +
                    try root_litter_ledger.totalCarbon(result.litterfall),
            );
            try ledger.validateCarbonAdd(
                plant,
                1,
                layer,
                try root_litter_ledger.totalCarbon(host_domain_one.litterfall),
            );
        }
        try root_disturbance.validateRootGasRelease(roots, plant, layer, 1 - remaining_fraction);
        try root_litterfall.validatePublication(organic, try grid.layerIndex(cell, layer), publication);
    }
    if (!std.math.isFinite(removed_host_carbon_g_c)) return error.NonFiniteRootHarvest;
    if (context.carbon_exchange_state) |exchange| {
        const branches = try context.canopy_state.branchRange(plant);
        if (branches.first >= branches.end or exchange.branchCount() != context.canopy_state.branch_node_offsets.len - 1) return error.CanopyCarbonExchangeDimensionMismatch;
        const next = exchange.disturbance_carbon_g_c_per_h[branches.first] + removed_host_carbon_g_c;
        if (!std.math.isFinite(next)) return error.NonFiniteRootHarvest;
    }
    for (0..grid.active_soil_layer_count[cell]) |layer| {
        const root = try roots.layerIndex(plant, 0, layer);
        const soil = try grid.layerIndex(cell, layer);
        const result = try noduleHarvestResult(
            nitrogen_fixation_type,
            .{
                .carbon_g_c = roots.symbiont_structural_carbon_g_c[root],
                .nitrogen_g_n = roots.symbiont_structural_nitrogen_g_n[root],
                .phosphorus_g_p = roots.symbiont_structural_phosphorus_g_p[root],
            },
            .{
                .carbon_g_c = roots.symbiont_mobile_carbon_g_c[root],
                .nitrogen_g_n = roots.symbiont_mobile_nitrogen_g_n[root],
                .phosphorus_g_p = roots.symbiont_mobile_phosphorus_g_p[root],
            },
            remaining_fraction,
            fine,
            mobile,
        );
        var host_litter_by_domain = [_]@import("plant_root_metabolism.zig").RootLitter{
            std.mem.zeroes(@import("plant_root_metabolism.zig").RootLitter),
            std.mem.zeroes(@import("plant_root_metabolism.zig").RootLitter),
        };
        for (0..root_system.biological_domain_count) |domain|
            host_litter_by_domain[domain] = (hostLayerHarvestDomain(
                roots,
                plant,
                domain,
                layer,
                remaining_fraction,
                woody_fraction,
                fine,
                coarse,
                mobile,
            ) catch unreachable).litterfall;
        roots.symbiont_structural_carbon_g_c[root] = result.structural.carbon_g_c;
        roots.symbiont_structural_nitrogen_g_n[root] = result.structural.nitrogen_g_n;
        roots.symbiont_structural_phosphorus_g_p[root] = result.structural.phosphorus_g_p;
        roots.symbiont_mobile_carbon_g_c[root] = result.mobile.carbon_g_c;
        roots.symbiont_mobile_nitrogen_g_n[root] = result.mobile.nitrogen_g_n;
        roots.symbiont_mobile_phosphorus_g_p[root] = result.mobile.phosphorus_g_p;
        var publication: root_litterfall.LayerInput = .{};
        publication.add(result.litterfall) catch unreachable;
        const host = hostLayerHarvest(roots, plant, layer, remaining_fraction, woody_fraction, fine, coarse, mobile) catch unreachable;
        publication.add(host.litterfall) catch unreachable;
        commitHostLayerHarvest(roots, plant, layer, remaining_fraction);
        root_disturbance.releaseRootGasFraction(roots, plant, layer, 1 - remaining_fraction) catch unreachable;
        root_litterfall.publishValidated(organic, soil, publication);
        if (context.root_litter_carbon_ledger) |ledger| {
            for (0..root_system.biological_domain_count) |domain|
                ledger.addValidated(plant, domain, layer, host_litter_by_domain[domain]);
            ledger.addValidated(plant, 0, layer, result.litterfall);
        }
    }
    if (context.carbon_exchange_state) |exchange| {
        const branches = context.canopy_state.branchRange(plant) catch unreachable;
        exchange.disturbance_carbon_g_c_per_h[branches.first] += removed_host_carbon_g_c;
    }
}

fn noduleHarvestResult(
    nitrogen_fixation_type: u8,
    structural: symbiotic_fixation.Pool,
    mobile: symbiotic_fixation.Pool,
    remaining_fraction: f64,
    structural_partition: litter_partition.ElementFractions,
    mobile_partition: litter_partition.ElementFractions,
) !root_disturbance.Result {
    if (try root_disturbance.sourceOrderNoduleHarvestIsEnabled(
        nitrogen_fixation_type,
        0,
        root_system.biological_domain_count,
    )) {
        return root_disturbance.retainAndRelease(
            structural,
            mobile,
            root_disturbance.ElementRetention.uniform(remaining_fraction),
            structural_partition,
            mobile_partition,
        );
    }
    return .{
        .structural = structural,
        .mobile = mobile,
        .litterfall = std.mem.zeroes(@import("plant_root_metabolism.zig").RootLitter),
    };
}

const HostLayerHarvest = struct {
    litterfall: @import("plant_root_metabolism.zig").RootLitter,
    removed_carbon_g_c: f64,
};

fn hostLayerHarvest(
    roots: *const root_system.State,
    plant: usize,
    layer: usize,
    remaining_fraction: f64,
    woody_fraction: f64,
    fine: litter_partition.ElementFractions,
    coarse: litter_partition.ElementFractions,
    mobile: litter_partition.ElementFractions,
) !HostLayerHarvest {
    return hostLayerHarvestRange(
        roots,
        plant,
        0,
        root_system.biological_domain_count,
        layer,
        remaining_fraction,
        woody_fraction,
        fine,
        coarse,
        mobile,
    );
}

fn hostLayerHarvestDomain(
    roots: *const root_system.State,
    plant: usize,
    domain: usize,
    layer: usize,
    remaining_fraction: f64,
    woody_fraction: f64,
    fine: litter_partition.ElementFractions,
    coarse: litter_partition.ElementFractions,
    mobile: litter_partition.ElementFractions,
) !HostLayerHarvest {
    if (domain >= root_system.biological_domain_count)
        return error.PlantRootIndexOutOfBounds;
    return hostLayerHarvestRange(
        roots,
        plant,
        domain,
        domain + 1,
        layer,
        remaining_fraction,
        woody_fraction,
        fine,
        coarse,
        mobile,
    );
}

fn hostLayerHarvestRange(
    roots: *const root_system.State,
    plant: usize,
    first_domain: usize,
    end_domain: usize,
    layer: usize,
    remaining_fraction: f64,
    woody_fraction: f64,
    fine: litter_partition.ElementFractions,
    coarse: litter_partition.ElementFractions,
    mobile: litter_partition.ElementFractions,
) !HostLayerHarvest {
    var result: HostLayerHarvest = .{ .litterfall = std.mem.zeroes(@import("plant_root_metabolism.zig").RootLitter), .removed_carbon_g_c = 0 };
    const removed_fraction = 1 - remaining_fraction;
    for (first_domain..end_domain) |domain| {
        const root = try roots.layerIndex(plant, domain, layer);
        const removed_mobile = canopy.ElementalMass{
            .carbon_g = roots.mobile_carbon_g[root] * removed_fraction,
            .nitrogen_g = roots.mobile_nitrogen_g[root] * removed_fraction,
            .phosphorus_g = roots.mobile_phosphorus_g[root] * removed_fraction,
        };
        result.removed_carbon_g_c += removed_mobile.carbon_g;
        addElementPartition(&result.litterfall, removed_mobile, mobile, false, 1);
        for (0..roots.active_root_axis_count[plant]) |axis| {
            const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
            const removed_structural = canopy.ElementalMass{
                .carbon_g = (roots.axis_primary_carbon_g[axis_layer] + roots.axis_secondary_carbon_g[axis_layer]) * removed_fraction,
                .nitrogen_g = (roots.axis_primary_nitrogen_g[axis_layer] + roots.axis_secondary_nitrogen_g[axis_layer]) * removed_fraction,
                .phosphorus_g = (roots.axis_primary_phosphorus_g[axis_layer] + roots.axis_secondary_phosphorus_g[axis_layer]) * removed_fraction,
            };
            result.removed_carbon_g_c += removed_structural.carbon_g;
            addElementPartition(&result.litterfall, removed_structural, coarse, true, woody_fraction);
            addElementPartition(&result.litterfall, removed_structural, fine, false, 1 - woody_fraction);
        }
    }
    return result;
}

fn addElementPartition(litter: *@import("plant_root_metabolism.zig").RootLitter, mass: canopy.ElementalMass, fractions: litter_partition.ElementFractions, woody: bool, multiplier: f64) void {
    for (0..root_litterfall.kinetic_component_count) |component| {
        if (woody) {
            litter.woody_carbon_g_c[component] += mass.carbon_g * multiplier * fractions.carbon[component];
            litter.woody_nitrogen_g_n[component] += mass.nitrogen_g * multiplier * fractions.nitrogen[component];
            litter.woody_phosphorus_g_p[component] += mass.phosphorus_g * multiplier * fractions.phosphorus[component];
        } else {
            litter.nonwoody_carbon_g_c[component] += mass.carbon_g * multiplier * fractions.carbon[component];
            litter.nonwoody_nitrogen_g_n[component] += mass.nitrogen_g * multiplier * fractions.nitrogen[component];
            litter.nonwoody_phosphorus_g_p[component] += mass.phosphorus_g * multiplier * fractions.phosphorus[component];
        }
    }
}

fn commitHostLayerHarvest(roots: *root_system.State, plant: usize, layer: usize, remaining_fraction: f64) void {
    for (0..root_system.biological_domain_count) |domain| {
        const root = roots.layerIndex(plant, domain, layer) catch unreachable;
        inline for (.{ "mobile_carbon_g", "mobile_nitrogen_g", "mobile_phosphorus_g" }) |field_name|
            @field(roots, field_name)[root] *= remaining_fraction;
        for (0..roots.active_root_axis_count[plant]) |axis| {
            const axis_layer = roots.layerAxisIndex(plant, domain, layer, axis) catch unreachable;
            inline for (.{
                "axis_primary_carbon_g",   "axis_primary_nitrogen_g",   "axis_primary_phosphorus_g",
                "axis_secondary_carbon_g", "axis_secondary_nitrogen_g", "axis_secondary_phosphorus_g",
            }) |field_name| @field(roots, field_name)[axis_layer] *= remaining_fraction;
        }
    }
}

fn harvestVegetativeBranch(context: *Context, plant: usize, branch: usize, event: management.HarvestEvent, science: ScienceParameters, pruning: bool) !void {
    const state = context.canopy_state;
    const nodes = try state.nodeRange(branch);
    const initial_leaf_sheath_c = state.branch_leaf_carbon_g[branch] + state.branch_sheath_carbon_g[branch];
    var maximum_height_m: f64 = 0;
    for (state.node_height_m[nodes.first..nodes.end]) |height_m| maximum_height_m = @max(maximum_height_m, height_m);
    for (nodes.first..nodes.end) |node| {
        const node_within_branch = node - nodes.first;
        const initial_leaf_c = state.node_leaf_carbon_g[node];
        const samples = try state.sampleRange(node);
        for (samples.first..samples.end) |sample| {
            var retention = try canopy.layerHarvestRetention(state.sample_layer_lower_height_m[sample], state.sample_layer_upper_height_m[sample], event.cutting_height_m_or_lai_fraction, pruning, event.kind == .none, event.thinning_fraction_or_consumption_rate, event.harvested_fraction.leaf);
            retention.unexported_fraction = unexportedFraction(retention.remaining_fraction, event.ecosystem_export_fraction.leaf);
            const products = try canopy.harvestLeafLayerSample(state, branch, node_within_branch, sample - samples.first, retention, science.carbon_woody_fraction, science.leaf_nitrogen_woody_fraction, science.leaf_phosphorus_woody_fraction, node_within_branch == 1);
            addProducts(&context.products_by_plant[plant].foliar, products.foliar);
            addProducts(&context.products_by_plant[plant].woody, products.woody);
        }
        const retention = try canopy.sourceOrderNodeOrganRetention(
            false,
            event.kind == .none,
            initial_leaf_c,
            state.node_leaf_carbon_g[node],
            event.harvested_fraction.leaf,
            event.harvested_fraction.nonfoliar,
            event.thinning_fraction_or_consumption_rate,
            context.plant_structural_presence_threshold_g_per_plant *
                state.plant_population_count[plant],
        );
        const sheath_products = try canopy.harvestNodeSheath(state, branch, node_within_branch, retention.remaining_fraction, retention.unexported_fraction, science.carbon_woody_fraction, science.sheath_nitrogen_woody_fraction, science.sheath_phosphorus_woody_fraction, @intFromEnum(event.kind) <= @intFromEnum(management.HarvestKind.above_ground), event.cutting_height_m_or_lai_fraction);
        addProducts(&context.products_by_plant[plant].nonfoliar, sheath_products.nonwoody);
        addProducts(&context.products_by_plant[plant].woody, sheath_products.woody);
        const internode_remaining = try canopy.internodeHarvestRetention(state.node_height_m[node], state.node_internode_length_m[node], event.cutting_height_m_or_lai_fraction, pruning, event.thinning_fraction_or_consumption_rate, event.harvested_fraction.woody, false, 0, state.branch_stalk_carbon_g[branch]);
        try canopy.commitInternodeHarvest(state, branch, node_within_branch, internode_remaining, @intFromEnum(event.kind) <= @intFromEnum(management.HarvestKind.above_ground) and event.thinning_fraction_or_consumption_rate == 0, event.cutting_height_m_or_lai_fraction);
    }
    const stalk_retention = try canopy.sourceOrderBranchStalkRetention(
        false,
        event.kind == .none,
        pruning,
        maximum_height_m,
        event.cutting_height_m_or_lai_fraction,
        event.thinning_fraction_or_consumption_rate,
        event.harvested_fraction.woody,
        state.branch_stalk_carbon_g[branch],
        0,
        0,
        context.plant_tissue_presence_threshold_g_per_plant *
            state.plant_population_count[plant],
    );
    const reserve_retention = try canopy.sourceOrderStalkReserveRetention(
        false,
        state.branch_stalk_carbon_g[branch] * stalk_retention.remaining_fraction,
        stalk_retention,
        state.branch_reserve_carbon_g[branch],
        0,
        context.plant_structural_presence_threshold_g_per_plant *
            state.plant_population_count[plant],
    );
    const stalk_products = try canopy.harvestBranchStalkAndReserve(state, branch, stalk_retention.remaining_fraction, stalk_retention.unexported_fraction, reserve_retention.remaining_fraction, reserve_retention.unexported_fraction);
    addProducts(&context.products_by_plant[plant].woody, stalk_products);
    const mobile_remaining = try canopy.sourceOrderNonGrazingMobileRetention(
        initial_leaf_sheath_c,
        state.branch_leaf_carbon_g[branch] + state.branch_sheath_carbon_g[branch],
        context.plant_structural_presence_threshold_g_per_plant *
            state.plant_population_count[plant],
    );
    const initial_mobile_carbon_g_c = state.branch_mobile_carbon_g[branch];
    const is_c4 = if (context.canopy_biochemistry_parameters_by_plant) |parameters|
        parameters[plant].pathway == .c4
    else
        true;
    const intermediate_remaining = try canopy.sourceOrderC4IntermediateRetention(
        is_c4,
        initial_mobile_carbon_g_c,
        initial_mobile_carbon_g_c * mobile_remaining,
        context.plant_structural_presence_threshold_g_per_plant *
            state.plant_population_count[plant],
    );
    const mobile_removed = try canopy.harvestBranchMobilePoolsWithIntermediateRetention(
        state,
        branch,
        mobile_remaining,
        intermediate_remaining,
    );
    const export_fraction = event.ecosystem_export_fraction.nonfoliar;
    addScaledMass(&context.products_by_plant[plant].nonstructural.ecosystem_export, mobile_removed, export_fraction);
    addScaledMass(&context.products_by_plant[plant].nonstructural.litter, mobile_removed, 1.0 - export_fraction);
}

fn validateScience(science: ScienceParameters) !void {
    inline for (@typeInfo(ScienceParameters).@"struct".fields) |field| {
        if (field.type != [2]f64) continue;
        const fractions = @field(science, field.name);
        if (!std.math.isFinite(fractions[0]) or !std.math.isFinite(fractions[1]) or fractions[0] < 0 or fractions[1] < 0 or @abs(fractions[0] + fractions[1] - 1) > 1e-8) return error.InvalidPlantHarvestScience;
    }
}

fn unexportedFraction(remaining: f64, ecosystem_export_fraction: f64) f64 {
    return remaining + (1.0 - remaining) * (1.0 - ecosystem_export_fraction);
}

fn addProducts(target: *canopy.HarvestProducts, source: canopy.HarvestProducts) void {
    addMass(&target.ecosystem_export, source.ecosystem_export);
    addMass(&target.litter, source.litter);
}

fn addMass(target: *canopy.ElementalMass, source: canopy.ElementalMass) void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| @field(target, field.name) += @field(source, field.name);
}

fn subtractMass(target: *canopy.ElementalMass, source: canopy.ElementalMass) void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        @field(target, field.name) = @max(0, @field(target, field.name) - @field(source, field.name));
}

fn addMassToStorage(state: *canopy.State, plant: usize, mass: canopy.ElementalMass) void {
    state.plant_seed_storage_carbon_g[plant] += mass.carbon_g;
    state.plant_seed_storage_nitrogen_g[plant] += mass.nitrogen_g;
    state.plant_seed_storage_phosphorus_g[plant] += mass.phosphorus_g;
}

fn addMassToStandingDead(state: *canopy.State, plant: usize, mass: canopy.ElementalMass, kinetics: litter_partition.ElementFractions) void {
    state.plant_standing_dead_carbon_g[plant] += mass.carbon_g;
    state.plant_standing_dead_nitrogen_g[plant] += mass.nitrogen_g;
    state.plant_standing_dead_phosphorus_g[plant] += mass.phosphorus_g;
    for (0..4) |kinetic| {
        const index = plant * 4 + kinetic;
        state.plant_standing_dead_carbon_by_kinetic_g[index] += mass.carbon_g * kinetics.carbon[kinetic];
        state.plant_standing_dead_nitrogen_by_kinetic_g[index] += mass.nitrogen_g * kinetics.nitrogen[kinetic];
        state.plant_standing_dead_phosphorus_by_kinetic_g[index] += mass.phosphorus_g * kinetics.phosphorus[kinetic];
    }
}

fn addScaledMass(target: *canopy.ElementalMass, source: canopy.ElementalMass, fraction: f64) void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| @field(target, field.name) += fraction * @field(source, field.name);
}

test "GROSUB JHVST two resets population and retains exports for reseeding" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{0});
    defer state.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    state.plant_population_per_m2[0] = 3;
    state.plant_population_count[0] = 6;
    state.plant_population_change_count[0] = 6;
    state.plant_standing_dead_population_count[0] = 2;
    state.plant_seed_storage_carbon_g[0] = 1;
    state.plant_seed_storage_nitrogen_g[0] = 0.1;
    state.plant_seed_storage_phosphorus_g[0] = 0.01;
    const science = [_]ScienceParameters{.{ .carbon_woody_fraction = .{ 0, 1 }, .leaf_nitrogen_woody_fraction = .{ 0, 1 }, .sheath_nitrogen_woody_fraction = .{ 0, 1 }, .leaf_phosphorus_woody_fraction = .{ 0, 1 }, .sheath_phosphorus_woody_fraction = .{ 0, 1 } }};
    var ledgers = [_]ProductLedger{.{}};
    ledgers[0].foliar.ecosystem_export = .{ .carbon_g = 4, .nitrogen_g = 0.4, .phosphorus_g = 0.04 };
    ledgers[0].standing_dead_export = .{ .carbon_g = 2, .nitrogen_g = 0.2, .phosphorus_g = 0.02 };
    const target_population_per_m2 = [_]f64{7};
    const cell_area_m2 = [_]f64{2};
    var plant_state = try phenology.State.init(std.testing.allocator, 1, 1);
    defer plant_state.deinit();
    plant_state.active[0] = true;
    plant_state.lifecycle_initialized[0] = true;
    var context: Context = .{
        .canopy_state = &state,
        .branch_development = &development,
        .science_by_plant = &science,
        .products_by_plant = &ledgers,
        .leaf_area_presence_tolerance_m2 = 1e-12,
        .reseed_population_per_m2_by_plant = &target_population_per_m2,
        .cell_area_m2_by_cell = &cell_area_m2,
        .plant_phenology = &plant_state,
    };
    try applyEventInternal(&context, 0, .{
        .date = .{ .day = 1, .month = 1, .year = 9999 },
        .kind = .grain,
        .termination = .terminate_and_reseed,
        .cutting_height_m_or_lai_fraction = 0,
        .thinning_fraction_or_consumption_rate = 1,
        .harvested_fraction = .{ .leaf = 1, .nonfoliar = 1, .woody = 1, .standing_dead = 1 },
        .ecosystem_export_fraction = .{ .leaf = 1, .nonfoliar = 1, .woody = 1, .standing_dead = 1 },
    }, false, false);
    try std.testing.expectApproxEqAbs(@as(f64, 7), state.plant_population_per_m2[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 14), state.plant_population_count[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 14), state.plant_population_change_count[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 14), state.plant_standing_dead_population_count[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 7), state.plant_seed_storage_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.7), state.plant_seed_storage_nitrogen_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.07), state.plant_seed_storage_phosphorus_g[0], 1e-12);
    try std.testing.expectEqual(canopy.ElementalMass{}, ledgers[0].foliar.ecosystem_export);
    try std.testing.expectEqual(canopy.ElementalMass{}, ledgers[0].standing_dead_export);
    try std.testing.expect(plant_state.reseed_pending[0]);
    try std.testing.expect(plant_state.active[0]);
    try std.testing.expect(plant_state.lifecycle_initialized[0]);
}

test "GROSUB automatic deciduous annual harvest fires once at reproductive turnover" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{0}, &.{});
    defer state.deinit();
    var layers = try canopy_layers.State.init(std.testing.allocator, 1, 1, 1, 1, 1, &state);
    defer layers.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    var plant_state = try phenology.State.init(std.testing.allocator, 1, 1);
    defer plant_state.deinit();
    plant_state.active[0] = true;
    plant_state.lifecycle_initialized[0] = true;
    var growth = try growth_stages.State.init(std.testing.allocator, &.{1});
    defer growth.deinit();
    var dormant = try dormancy.RuntimeState.init(std.testing.allocator, 1);
    defer dormant.deinit();
    dormant.branches[0].phenological_remobilization_enabled = true;
    state.branch_grain_carbon_g[0] = 5;
    state.branch_grain_nitrogen_g[0] = 0.5;
    state.branch_grain_phosphorus_g[0] = 0.05;
    state.plant_population_per_m2[0] = 2;
    const science = [_]ScienceParameters{.{ .carbon_woody_fraction = .{ 0, 1 }, .leaf_nitrogen_woody_fraction = .{ 0, 1 }, .sheath_nitrogen_woody_fraction = .{ 0, 1 }, .leaf_phosphorus_woody_fraction = .{ 0, 1 }, .sheath_phosphorus_woody_fraction = .{ 0, 1 } }};
    var ledgers = [_]ProductLedger{.{}};
    const population = [_]f64{4};
    const area = [_]f64{2};
    const root_woody = [_]f64{0};
    var context: Context = .{
        .canopy_state = &state,
        .canopy_layer_state = &layers,
        .branch_development = &development,
        .science_by_plant = &science,
        .products_by_plant = &ledgers,
        .leaf_area_presence_tolerance_m2 = 1e-12,
        .plant_phenology = &plant_state,
        .growth_stages = &growth,
        .reseed_population_per_m2_by_plant = &population,
        .cell_area_m2_by_cell = &area,
        .root_woody_fraction_by_plant = &root_woody,
    };
    try std.testing.expectEqual(@as(usize, 1), try applyAutomaticSelfSeedingHarvests(&context, &dormant, &.{0}, &.{1}));
    try std.testing.expectApproxEqAbs(@as(f64, 5), state.plant_seed_storage_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 4), state.plant_population_per_m2[0], 1e-12);
    try std.testing.expect(plant_state.reseed_pending[0]);
    try std.testing.expectEqual(@as(usize, 0), try applyAutomaticSelfSeedingHarvests(&context, &dormant, &.{0}, &.{1}));
}

test "GROSUB perennial start-of-season residue is conserved before reconstruction" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{1});
    defer state.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    var partitions = try litter_partition.State.init(std.testing.allocator, 1);
    defer partitions.deinit();
    @memset(partitions.by_plant_and_organ, .{ .carbon = .{ 0.1, 0.2, 0.3, 0.4 }, .nitrogen = .{ 0.1, 0.2, 0.3, 0.4 }, .phosphorus = .{ 0.1, 0.2, 0.3, 0.4 } });
    state.sample_leaf_area_m2[0] = 1;
    state.sample_leaf_carbon_g[0] = 2;
    state.sample_leaf_nitrogen_g[0] = 0.2;
    state.sample_leaf_phosphorus_g[0] = 0.02;
    state.node_leaf_area_m2[0] = 1;
    state.node_leaf_carbon_g[0] = 2;
    state.node_leaf_nitrogen_g[0] = 0.2;
    state.node_leaf_phosphorus_g[0] = 0.02;
    state.branch_leaf_area_m2[0] = 1;
    state.branch_leaf_carbon_g[0] = 2;
    state.branch_leaf_nitrogen_g[0] = 0.2;
    state.branch_leaf_phosphorus_g[0] = 0.02;
    state.node_sheath_carbon_g[0] = 1;
    state.node_sheath_nitrogen_g[0] = 0.1;
    state.node_sheath_phosphorus_g[0] = 0.01;
    state.branch_sheath_carbon_g[0] = 1;
    state.branch_sheath_nitrogen_g[0] = 0.1;
    state.branch_sheath_phosphorus_g[0] = 0.01;
    state.branch_husk_carbon_g[0] = 1;
    state.branch_ear_carbon_g[0] = 1;
    state.branch_grain_carbon_g[0] = 1;
    state.branch_stalk_carbon_g[0] = 4;
    state.branch_stalk_nitrogen_g[0] = 0.4;
    state.branch_stalk_phosphorus_g[0] = 0.04;
    state.branch_reserve_carbon_g[0] = 2;
    state.branch_reserve_nitrogen_g[0] = 0.2;
    state.branch_reserve_phosphorus_g[0] = 0.02;
    const science = [_]ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.25, 0.75 }, .sheath_nitrogen_woody_fraction = .{ 0.25, 0.75 }, .leaf_phosphorus_woody_fraction = .{ 0.25, 0.75 }, .sheath_phosphorus_woody_fraction = .{ 0.25, 0.75 } }};
    var ledgers = [_]ProductLedger{.{}};
    var context: Context = .{ .canopy_state = &state, .branch_development = &development, .science_by_plant = &science, .products_by_plant = &ledgers, .leaf_area_presence_tolerance_m2 = 1e-12, .root_litter_partition = &partitions };
    try applyStartOfSeasonResidue(&context, 0, 0, 1);
    try std.testing.expectApproxEqAbs(@as(f64, 4), state.plant_standing_dead_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2), state.plant_seed_storage_carbon_g[0], 1e-12);
    const litter_carbon_g_c = ledgers[0].foliar.litter.carbon_g + ledgers[0].nonfoliar.litter.carbon_g + ledgers[0].woody.litter.carbon_g;
    try std.testing.expectApproxEqAbs(@as(f64, 6), litter_carbon_g_c, 1e-12);
    var standing_kinetic_carbon_g_c: f64 = 0;
    for (state.plant_standing_dead_carbon_by_kinetic_g[0..4]) |value| standing_kinetic_carbon_g_c += value;
    try std.testing.expectApproxEqAbs(@as(f64, 4), standing_kinetic_carbon_g_c, 1e-12);
    try std.testing.expectEqual(@as(f64, 0), state.branch_stalk_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_reserve_carbon_g[0]);
}

test "GROSUB whole-plant death conserves shoot storage symbiont and standing dead carbon" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{1});
    defer state.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    var partitions = try litter_partition.State.init(std.testing.allocator, 1);
    defer partitions.deinit();
    @memset(partitions.by_plant_and_organ, .{
        .carbon = .{ 0.1, 0.2, 0.3, 0.4 },
        .nitrogen = .{ 0.1, 0.2, 0.3, 0.4 },
        .phosphorus = .{ 0.1, 0.2, 0.3, 0.4 },
    });
    state.plant_seed_storage_carbon_g[0] = 4;
    state.branch_stalk_carbon_g[0] = 3;
    state.branch_reserve_carbon_g[0] = 2;
    state.branch_mobile_carbon_g[0] = 1;
    state.node_c4_mesophyll_nonstructural_carbon_g[0] = 1;
    state.branch_symbiont_mobile_carbon_g[0] = 0.5;
    state.branch_symbiont_structural_carbon_g[0] = 0.5;
    state.branch_husk_carbon_g[0] = 1;
    const science = [_]ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.25, 0.75 }, .sheath_nitrogen_woody_fraction = .{ 0.25, 0.75 }, .leaf_phosphorus_woody_fraction = .{ 0.25, 0.75 }, .sheath_phosphorus_woody_fraction = .{ 0.25, 0.75 } }};
    var ledgers = [_]ProductLedger{.{}};
    var context: Context = .{
        .canopy_state = &state,
        .branch_development = &development,
        .science_by_plant = &science,
        .products_by_plant = &ledgers,
        .leaf_area_presence_tolerance_m2 = 1e-12,
        .root_litter_partition = &partitions,
    };

    try applyWholePlantMortalityResidue(&context, 0);
    var direct_litter_carbon_g_c: f64 = 0;
    var direct_woody_carbon_g_c: f64 = 0;
    var direct_nonwoody_carbon_g_c: f64 = 0;
    for (ledgers[0].direct_litter.woody_carbon_g, ledgers[0].direct_litter.nonwoody_carbon_g) |woody, nonwoody| {
        direct_woody_carbon_g_c += woody;
        direct_nonwoody_carbon_g_c += nonwoody;
        direct_litter_carbon_g_c += woody + nonwoody;
    }
    const litter_carbon_g_c = direct_litter_carbon_g_c + ledgers[0].nonstructural.litter.carbon_g +
        ledgers[0].foliar.litter.carbon_g +
        ledgers[0].nonfoliar.litter.carbon_g +
        ledgers[0].woody.litter.carbon_g;
    try std.testing.expectApproxEqAbs(@as(f64, 5), state.plant_standing_dead_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 8), litter_carbon_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 13), state.plant_standing_dead_carbon_g[0] + litter_carbon_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), direct_woody_carbon_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), direct_nonwoody_carbon_g_c, 1e-12);
    try std.testing.expectEqual(@as(f64, 0), state.plant_seed_storage_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_stalk_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_reserve_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_mobile_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.node_c4_mesophyll_nonstructural_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_symbiont_mobile_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_symbiont_structural_carbon_g[0]);
}

test "GROSUB natural winter-annual branch death recovers mobile reserve and grain" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{0});
    defer state.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    var partitions = try litter_partition.State.init(std.testing.allocator, 1);
    defer partitions.deinit();
    @memset(partitions.by_plant_and_organ, .{
        .carbon = .{ 0.1, 0.2, 0.3, 0.4 },
        .nitrogen = .{ 0.1, 0.2, 0.3, 0.4 },
        .phosphorus = .{ 0.1, 0.2, 0.3, 0.4 },
    });
    state.branch_mobile_carbon_g[0] = 2;
    state.node_c4_mesophyll_nonstructural_carbon_g[0] = 1;
    state.branch_reserve_carbon_g[0] = 3;
    state.branch_grain_carbon_g[0] = 4;
    state.branch_stalk_carbon_g[0] = 5;
    state.branch_symbiont_structural_carbon_g[0] = 1;
    const science = [_]ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.25, 0.75 }, .sheath_nitrogen_woody_fraction = .{ 0.25, 0.75 }, .leaf_phosphorus_woody_fraction = .{ 0.25, 0.75 }, .sheath_phosphorus_woody_fraction = .{ 0.25, 0.75 } }};
    var ledgers = [_]ProductLedger{.{}};
    var context: Context = .{ .canopy_state = &state, .branch_development = &development, .science_by_plant = &science, .products_by_plant = &ledgers, .leaf_area_presence_tolerance_m2 = 1e-12, .root_litter_partition = &partitions };

    try applyNaturalDeadBranchResidue(&context, 0, 0, true);
    try std.testing.expectApproxEqAbs(@as(f64, 10), state.plant_seed_storage_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5), state.plant_standing_dead_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), ledgers[0].foliar.litter.carbon_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 16), state.plant_seed_storage_carbon_g[0] + state.plant_standing_dead_carbon_g[0] + ledgers[0].foliar.litter.carbon_g, 1e-12);
    try std.testing.expectEqual(@as(f64, 0), state.branch_mobile_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.node_c4_mesophyll_nonstructural_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), state.branch_grain_carbon_g[0]);
}

test "non-grazing runtime callback conserves branch carbon through cutting" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{1};
    const sample_counts = [_]usize{1};
    var canopy_state = try canopy.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer canopy_state.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    canopy_state.sample_leaf_area_m2[0] = 2;
    canopy_state.sample_layer_lower_height_m[0] = 0;
    canopy_state.sample_layer_upper_height_m[0] = 2;
    canopy_state.sample_leaf_carbon_g[0] = 4;
    canopy_state.sample_leaf_nitrogen_g[0] = 0.4;
    canopy_state.sample_leaf_phosphorus_g[0] = 0.04;
    canopy_state.node_leaf_area_m2[0] = 2;
    canopy_state.node_leaf_carbon_g[0] = 4;
    canopy_state.node_leaf_nitrogen_g[0] = 0.4;
    canopy_state.node_leaf_phosphorus_g[0] = 0.04;
    canopy_state.branch_leaf_area_m2[0] = 2;
    canopy_state.branch_leaf_carbon_g[0] = 4;
    canopy_state.branch_leaf_nitrogen_g[0] = 0.4;
    canopy_state.branch_leaf_phosphorus_g[0] = 0.04;
    canopy_state.node_height_m[0] = 2;
    canopy_state.node_sheath_height_m[0] = 1;
    canopy_state.node_sheath_carbon_g[0] = 2;
    canopy_state.node_sheath_nitrogen_g[0] = 0.2;
    canopy_state.node_sheath_phosphorus_g[0] = 0.02;
    canopy_state.branch_sheath_carbon_g[0] = 2;
    canopy_state.branch_sheath_nitrogen_g[0] = 0.2;
    canopy_state.branch_sheath_phosphorus_g[0] = 0.02;
    canopy_state.node_internode_length_m[0] = 2;
    canopy_state.node_internode_carbon_g[0] = 4;
    canopy_state.node_internode_nitrogen_g[0] = 0.4;
    canopy_state.node_internode_phosphorus_g[0] = 0.04;
    canopy_state.branch_stalk_carbon_g[0] = 4;
    canopy_state.branch_stalk_nitrogen_g[0] = 0.4;
    canopy_state.branch_stalk_phosphorus_g[0] = 0.04;
    canopy_state.branch_reserve_carbon_g[0] = 1;
    canopy_state.branch_reserve_nitrogen_g[0] = 0.1;
    canopy_state.branch_reserve_phosphorus_g[0] = 0.01;
    canopy_state.branch_mobile_carbon_g[0] = 2;
    canopy_state.branch_mobile_nitrogen_g[0] = 0.2;
    canopy_state.branch_mobile_phosphorus_g[0] = 0.02;
    const science = [_]ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .sheath_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .leaf_phosphorus_woody_fraction = .{ 0.1, 0.9 }, .sheath_phosphorus_woody_fraction = .{ 0.1, 0.9 } }};
    var ledgers = [_]ProductLedger{.{}};
    var context: Context = .{ .canopy_state = &canopy_state, .branch_development = &development, .science_by_plant = &science, .products_by_plant = &ledgers, .leaf_area_presence_tolerance_m2 = 1.0e-12 };
    const event: management.HarvestEvent = .{ .date = .{ .day = 1, .month = 1, .year = 9999 }, .kind = .above_ground, .termination = .retain, .cutting_height_m_or_lai_fraction = 1, .thinning_fraction_or_consumption_rate = 0, .harvested_fraction = .{ .leaf = 1, .nonfoliar = 1, .woody = 1, .standing_dead = 0 }, .ecosystem_export_fraction = .{ .leaf = 0.8, .nonfoliar = 0.8, .woody = 0.8, .standing_dead = 0 } };
    try applyEvent(&context, 0, event);
    const remaining_c = canopy_state.branch_leaf_carbon_g[0] + canopy_state.branch_sheath_carbon_g[0] + canopy_state.branch_stalk_carbon_g[0] + canopy_state.branch_reserve_carbon_g[0] + canopy_state.branch_mobile_carbon_g[0];
    const products_c = ledgers[0].nonstructural.ecosystem_export.carbon_g + ledgers[0].nonstructural.litter.carbon_g + ledgers[0].foliar.ecosystem_export.carbon_g + ledgers[0].foliar.litter.carbon_g + ledgers[0].nonfoliar.ecosystem_export.carbon_g + ledgers[0].nonfoliar.litter.carbon_g + ledgers[0].woody.ecosystem_export.carbon_g + ledgers[0].woody.litter.carbon_g;
    try std.testing.expectApproxEqAbs(13.0, remaining_c + products_c, 1e-12);
    try std.testing.expectApproxEqAbs(2.0, canopy_state.branch_leaf_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(1.0, canopy_state.branch_sheath_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(2.0, canopy_state.branch_stalk_carbon_g[0], 1e-15);
    try std.testing.expectApproxEqAbs(1.0, canopy_state.branch_mobile_carbon_g[0], 1e-15);

    canopy_state.branch_grain_carbon_g[0] = 3;
    canopy_state.branch_grain_nitrogen_g[0] = 0.3;
    canopy_state.branch_grain_phosphorus_g[0] = 0.03;
    var high_cut = event;
    high_cut.cutting_height_m_or_lai_fraction = 3;
    try applyEvent(&context, 0, high_cut);
    try std.testing.expectApproxEqAbs(
        3,
        canopy_state.branch_grain_carbon_g[0],
        1e-15,
    );
}

test "runtime callback rejects grazing approximation" {
    const event: management.HarvestEvent = .{ .date = .{ .day = 1, .month = 1, .year = 9999 }, .kind = .animal_grazing, .termination = .retain, .cutting_height_m_or_lai_fraction = 0, .thinning_fraction_or_consumption_rate = 0, .harvested_fraction = .{ .leaf = 0, .nonfoliar = 0, .woody = 0, .standing_dead = 0 }, .ecosystem_export_fraction = .{ .leaf = 0, .nonfoliar = 0, .woody = 0, .standing_dead = 0 } };
    var context: Context = undefined;
    try std.testing.expectError(error.GrazingRequiresDemandDrivenKernel, applyEvent(&context, 0, event));
}

test "GROSUB grazing removes top layers and newest nodes first and conserves C N P" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{2};
    const sample_counts = [_]usize{ 2, 2 };
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer state.deinit();
    var layers = try canopy_layers.State.init(std.testing.allocator, 1, 1, 2, 1, 1, &state);
    defer layers.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();

    // Old node: bottom=2, top=1. New node: bottom=1, top=2 g C.
    const carbon = [_]f64{ 2, 1, 1, 2 };
    for (0..2) |node| {
        var node_carbon_g_c: f64 = 0;
        for (0..2) |layer| {
            const node_layer = node * 2 + layer;
            const sample = node_layer;
            const carbon_g_c = carbon[node_layer];
            layers.node_leaf_area_m2[node_layer] = carbon_g_c;
            layers.node_leaf_carbon_g[node_layer] = carbon_g_c;
            layers.node_leaf_nitrogen_g[node_layer] = 0.1 * carbon_g_c;
            layers.node_leaf_phosphorus_g[node_layer] = 0.01 * carbon_g_c;
            state.sample_leaf_area_m2[sample] = carbon_g_c;
            state.sample_exposed_leaf_area_m2[sample] = carbon_g_c;
            state.sample_leaf_carbon_g[sample] = carbon_g_c;
            state.sample_leaf_nitrogen_g[sample] = 0.1 * carbon_g_c;
            state.sample_leaf_phosphorus_g[sample] = 0.01 * carbon_g_c;
            node_carbon_g_c += carbon_g_c;
        }
        state.node_leaf_area_m2[node] = node_carbon_g_c;
        state.node_leaf_carbon_g[node] = node_carbon_g_c;
        state.node_leaf_nitrogen_g[node] = 0.1 * node_carbon_g_c;
        state.node_leaf_phosphorus_g[node] = 0.01 * node_carbon_g_c;
    }
    state.node_sheath_carbon_g[0] = 3;
    state.node_sheath_nitrogen_g[0] = 0.3;
    state.node_sheath_phosphorus_g[0] = 0.03;
    state.node_sheath_carbon_g[1] = 1;
    state.node_sheath_nitrogen_g[1] = 0.1;
    state.node_sheath_phosphorus_g[1] = 0.01;
    state.branch_leaf_area_m2[0] = 6;
    state.branch_leaf_carbon_g[0] = 6;
    state.branch_leaf_nitrogen_g[0] = 0.6;
    state.branch_leaf_phosphorus_g[0] = 0.06;
    state.branch_sheath_carbon_g[0] = 4;
    state.branch_sheath_nitrogen_g[0] = 0.4;
    state.branch_sheath_phosphorus_g[0] = 0.04;
    state.plant_total_shoot_carbon_g[0] = 10;
    state.plant_uptake_growth_temperature_response[0] = 1;
    state.plant_mobile_carbon_concentration_g_per_g[0] = 0;
    state.plant_symbiont_mobile_carbon_concentration_g_per_g[0] = 0;

    const science = [_]ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .sheath_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .leaf_phosphorus_woody_fraction = .{ 0.1, 0.9 }, .sheath_phosphorus_woody_fraction = .{ 0.1, 0.9 } }};
    var ledgers = [_]ProductLedger{.{}};
    var context: Context = .{ .canopy_state = &state, .canopy_layer_state = &layers, .branch_development = &development, .science_by_plant = &science, .products_by_plant = &ledgers, .leaf_area_presence_tolerance_m2 = 1.0e-12 };
    const consumed_g_c = try applyGrazingEvent(&context, 0, .{
        .date = .{ .day = 1, .month = 1, .year = 9999 },
        .kind = .animal_grazing,
        .termination = .retain,
        .cutting_height_m_or_lai_fraction = 192,
        .thinning_fraction_or_consumption_rate = 1,
        .harvested_fraction = .{ .leaf = 0.5, .nonfoliar = 0.5, .woody = 0, .standing_dead = 0 },
        .ecosystem_export_fraction = .{ .leaf = 0.25, .nonfoliar = 0.25, .woody = 0, .standing_dead = 0 },
    }, 10, 1);
    try std.testing.expectApproxEqAbs(@as(f64, 4), consumed_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), state.node_leaf_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), state.node_leaf_carbon_g[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2), state.node_sheath_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), state.node_sheath_carbon_g[1], 1e-12);
    var products: canopy.ElementalMass = .{};
    inline for (.{ "foliar", "woody", "nonfoliar" }) |field_name| {
        addMass(&products, @field(ledgers[0], field_name).ecosystem_export);
        addMass(&products, @field(ledgers[0], field_name).litter);
    }
    addMass(&products, ledgers[0].standing_dead_export);
    for (ledgers[0].manure.organic_by_biochemical_fraction) |mass| addMass(&products, mass);
    products.nitrogen_g += ledgers[0].manure.inorganic_nitrogen_g_n;
    products.phosphorus_g += ledgers[0].manure.inorganic_phosphorus_g_p;
    try std.testing.expectApproxEqAbs(@as(f64, 10), state.branch_leaf_carbon_g[0] + state.branch_sheath_carbon_g[0] + products.carbon_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), state.branch_leaf_nitrogen_g[0] + state.branch_sheath_nitrogen_g[0] + products.nitrogen_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), state.branch_leaf_phosphorus_g[0] + state.branch_sheath_phosphorus_g[0] + products.phosphorus_g, 1e-12);
}

test "GROSUB grazing applies reserve demand independently to every runtime branch" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{2}, &.{ 1, 1 }, &.{ 1, 1 });
    defer state.deinit();
    var layers = try canopy_layers.State.init(std.testing.allocator, 1, 1, 1, 1, 1, &state);
    defer layers.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 2);
    defer development.deinit();
    for (0..2) |branch| {
        state.branch_stalk_carbon_g[branch] = 5;
        state.branch_stalk_nitrogen_g[branch] = 0.5;
        state.branch_stalk_phosphorus_g[branch] = 0.05;
        state.node_internode_carbon_g[branch] = 5;
        state.node_internode_nitrogen_g[branch] = 0.5;
        state.node_internode_phosphorus_g[branch] = 0.05;
    }
    state.branch_reserve_carbon_g[0] = 1;
    state.branch_reserve_nitrogen_g[0] = 0.1;
    state.branch_reserve_phosphorus_g[0] = 0.01;
    state.branch_reserve_carbon_g[1] = 3;
    state.branch_reserve_nitrogen_g[1] = 0.3;
    state.branch_reserve_phosphorus_g[1] = 0.03;
    state.plant_total_shoot_carbon_g[0] = 14;
    state.plant_uptake_growth_temperature_response[0] = 1;
    const science = [_]ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .sheath_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .leaf_phosphorus_woody_fraction = .{ 0.1, 0.9 }, .sheath_phosphorus_woody_fraction = .{ 0.1, 0.9 } }};
    var ledgers = [_]ProductLedger{.{}};
    var context: Context = .{ .canopy_state = &state, .canopy_layer_state = &layers, .branch_development = &development, .science_by_plant = &science, .products_by_plant = &ledgers, .leaf_area_presence_tolerance_m2 = 1e-12 };
    const event: management.HarvestEvent = .{
        .date = .{ .day = 1, .month = 1, .year = 9999 },
        .kind = .animal_grazing,
        .termination = .retain,
        .cutting_height_m_or_lai_fraction = 192,
        .thinning_fraction_or_consumption_rate = 1,
        .harvested_fraction = .{ .leaf = 0, .nonfoliar = 0, .woody = 1, .standing_dead = 0 },
        .ecosystem_export_fraction = .{ .leaf = 0, .nonfoliar = 0, .woody = 1, .standing_dead = 0 },
    };
    const removed_g_c = try applyGrazingEvent(&context, 0, event, 14, 1);
    const reserve_target_g_c = 4.0 * 4.0 / 14.0;
    try std.testing.expectApproxEqAbs(@as(f64, 0), state.branch_reserve_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(3.0 - reserve_target_g_c, state.branch_reserve_carbon_g[1], 1e-12);
    try std.testing.expectApproxEqAbs(4.0 * 10.0 / 14.0 + 1.0 + reserve_target_g_c, removed_g_c, 1e-12);
    const remaining_g_c = state.branch_stalk_carbon_g[0] + state.branch_stalk_carbon_g[1] +
        state.branch_reserve_carbon_g[0] + state.branch_reserve_carbon_g[1];
    try std.testing.expectApproxEqAbs(@as(f64, 14), remaining_g_c + ledgers[0].woody.ecosystem_export.carbon_g, 1e-12);

    // A late-organ invalid pool is rejected by preflight before any earlier
    // organ or product ledger can be changed.
    state.branch_reserve_carbon_g[0] = -1;
    const stalk_before = state.branch_stalk_carbon_g[0];
    const products_before = productLedgerCarbonG(ledgers[0]);
    try std.testing.expectError(error.InvalidGrazingState, applyGrazingEvent(&context, 0, event, 14, 1));
    try std.testing.expectEqual(stalk_before, state.branch_stalk_carbon_g[0]);
    try std.testing.expectEqual(products_before, productLedgerCarbonG(ledgers[0]));
}

test "GROSUB grazing counts reproductive grain removal exactly once" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{0}, &.{});
    defer state.deinit();
    var layers = try canopy_layers.State.init(std.testing.allocator, 1, 1, 1, 1, 1, &state);
    defer layers.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    state.branch_husk_carbon_g[0] = 2;
    state.branch_husk_nitrogen_g[0] = 0.2;
    state.branch_husk_phosphorus_g[0] = 0.02;
    state.branch_ear_carbon_g[0] = 2;
    state.branch_ear_nitrogen_g[0] = 0.2;
    state.branch_ear_phosphorus_g[0] = 0.02;
    state.branch_grain_carbon_g[0] = 2;
    state.branch_grain_nitrogen_g[0] = 0.2;
    state.branch_grain_phosphorus_g[0] = 0.02;
    state.plant_total_shoot_carbon_g[0] = 6;
    state.plant_uptake_growth_temperature_response[0] = 1;
    const science = [_]ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .sheath_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .leaf_phosphorus_woody_fraction = .{ 0.1, 0.9 }, .sheath_phosphorus_woody_fraction = .{ 0.1, 0.9 } }};
    var ledgers = [_]ProductLedger{.{}};
    var context: Context = .{ .canopy_state = &state, .canopy_layer_state = &layers, .branch_development = &development, .science_by_plant = &science, .products_by_plant = &ledgers, .leaf_area_presence_tolerance_m2 = 1e-12 };
    const event: management.HarvestEvent = .{
        .date = .{ .day = 1, .month = 1, .year = 9999 },
        .kind = .animal_grazing,
        .termination = .retain,
        .cutting_height_m_or_lai_fraction = 144,
        .thinning_fraction_or_consumption_rate = 1,
        .harvested_fraction = .{ .leaf = 0, .nonfoliar = 1, .woody = 0, .standing_dead = 0 },
        .ecosystem_export_fraction = .{ .leaf = 0, .nonfoliar = 0.5, .woody = 0, .standing_dead = 0 },
    };
    const removed_g_c = try applyGrazingEvent(&context, 0, event, 6, 1);
    const remaining_g_c = state.branch_husk_carbon_g[0] + state.branch_ear_carbon_g[0] + state.branch_grain_carbon_g[0];
    try std.testing.expectApproxEqAbs(@as(f64, 3), removed_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), remaining_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), ledgers[0].nonfoliar.ecosystem_export.carbon_g, 1e-12);
    var manure_carbon_g_c: f64 = 0;
    for (ledgers[0].manure.organic_by_biochemical_fraction) |mass| manure_carbon_g_c += mass.carbon_g;
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), manure_carbon_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 6), remaining_g_c + productLedgerCarbonG(ledgers[0]), 1e-12);

    state.branch_husk_carbon_g[0] = 2;
    state.branch_ear_carbon_g[0] = 2;
    state.branch_grain_carbon_g[0] = 2;
    state.plant_population_count[0] = 1;
    ledgers[0] = .{};
    context.plant_structural_presence_threshold_g_per_plant = 2;
    _ = try applyGrazingEvent(&context, 0, event, 6, 1);
    try std.testing.expectApproxEqAbs(
        6,
        state.branch_husk_carbon_g[0] +
            state.branch_ear_carbon_g[0] +
            state.branch_grain_carbon_g[0],
        1e-12,
    );
    try std.testing.expectApproxEqAbs(0, productLedgerCarbonG(ledgers[0]), 1e-12);
}

test "GROSUB grazing distributes host and symbiont pools by branch leaf sheath mass" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{2}, &.{ 1, 1 }, &.{ 1, 1 });
    defer state.deinit();
    var layers = try canopy_layers.State.init(std.testing.allocator, 1, 1, 1, 1, 1, &state);
    defer layers.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 2);
    defer development.deinit();
    for (0..2) |branch| {
        const leaf_g_c: f64 = if (branch == 0) 3 else 1;
        state.sample_leaf_area_m2[branch] = leaf_g_c;
        state.sample_exposed_leaf_area_m2[branch] = leaf_g_c;
        state.sample_leaf_carbon_g[branch] = leaf_g_c;
        state.node_leaf_area_m2[branch] = leaf_g_c;
        state.node_leaf_carbon_g[branch] = leaf_g_c;
        state.branch_leaf_area_m2[branch] = leaf_g_c;
        state.branch_leaf_carbon_g[branch] = leaf_g_c;
        layers.node_leaf_area_m2[branch] = leaf_g_c;
        layers.node_leaf_carbon_g[branch] = leaf_g_c;
        state.branch_mobile_carbon_g[branch] = if (branch == 0) 1 else 3;
        state.branch_symbiont_mobile_carbon_g[branch] = 1;
        state.branch_symbiont_structural_carbon_g[branch] = 2;
    }
    state.plant_total_shoot_carbon_g[0] = 14;
    state.plant_uptake_growth_temperature_response[0] = 1;
    state.plant_mobile_carbon_concentration_g_per_g[0] = 1;
    state.plant_symbiont_mobile_carbon_concentration_g_per_g[0] = 1;
    const science = [_]ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .sheath_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .leaf_phosphorus_woody_fraction = .{ 0.1, 0.9 }, .sheath_phosphorus_woody_fraction = .{ 0.1, 0.9 } }};
    var ledgers = [_]ProductLedger{.{}};
    var context: Context = .{ .canopy_state = &state, .canopy_layer_state = &layers, .branch_development = &development, .science_by_plant = &science, .products_by_plant = &ledgers, .leaf_area_presence_tolerance_m2 = 1e-12 };
    const removed_g_c = try applyGrazingEvent(&context, 0, .{
        .date = .{ .day = 1, .month = 1, .year = 9999 },
        .kind = .animal_grazing,
        .termination = .retain,
        .cutting_height_m_or_lai_fraction = 96,
        .thinning_fraction_or_consumption_rate = 1,
        .harvested_fraction = .{ .leaf = 1, .nonfoliar = 0, .woody = 0, .standing_dead = 0 },
        .ecosystem_export_fraction = .{ .leaf = 1, .nonfoliar = 1, .woody = 1, .standing_dead = 0 },
    }, 14, 1);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), state.branch_mobile_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2.75), state.branch_mobile_carbon_g[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), state.branch_symbiont_mobile_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.75), state.branch_symbiont_mobile_carbon_g[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), state.branch_symbiont_structural_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), state.branch_symbiont_structural_carbon_g[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5), removed_g_c, 1e-12);
}

test "GROSUB standing dead grazing updates mass area export and manure" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{0}, &.{});
    defer state.deinit();
    var layers = try canopy_layers.State.init(std.testing.allocator, 1, 1, 1, 1, 1, &state);
    defer layers.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    state.plant_standing_dead_carbon_g[0] = 8;
    state.plant_standing_dead_nitrogen_g[0] = 0.8;
    state.plant_standing_dead_phosphorus_g[0] = 0.08;
    state.plant_charcoal_carbon_g[0] = 4;
    state.plant_charcoal_nitrogen_g[0] = 0.4;
    state.plant_charcoal_phosphorus_g[0] = 0.04;
    for (0..4) |fraction| {
        state.plant_standing_dead_carbon_by_kinetic_g[fraction] = 2;
        state.plant_standing_dead_nitrogen_by_kinetic_g[fraction] = 0.2;
        state.plant_standing_dead_phosphorus_by_kinetic_g[fraction] = 0.02;
    }
    layers.plant_standing_dead_area_m2[0] = 5;
    layers.plant_standing_dead_projected_surface_m2[0] = 5;
    layers.cell_standing_dead_area_m2[0] = 5;
    state.plant_uptake_growth_temperature_response[0] = 1;
    const science = [_]ScienceParameters{.{ .carbon_woody_fraction = .{ 0.25, 0.75 }, .leaf_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .sheath_nitrogen_woody_fraction = .{ 0.2, 0.8 }, .leaf_phosphorus_woody_fraction = .{ 0.1, 0.9 }, .sheath_phosphorus_woody_fraction = .{ 0.1, 0.9 } }};
    var ledgers = [_]ProductLedger{.{}};
    var context: Context = .{ .canopy_state = &state, .canopy_layer_state = &layers, .branch_development = &development, .science_by_plant = &science, .products_by_plant = &ledgers, .leaf_area_presence_tolerance_m2 = 1e-12 };
    const removed_g_c = try applyGrazingEvent(&context, 0, .{
        .date = .{ .day = 1, .month = 1, .year = 9999 },
        .kind = .animal_grazing,
        .termination = .retain,
        .cutting_height_m_or_lai_fraction = 96,
        .thinning_fraction_or_consumption_rate = 1,
        .harvested_fraction = .{ .leaf = 0, .nonfoliar = 0, .woody = 0, .standing_dead = 0.5 },
        .ecosystem_export_fraction = .{ .leaf = 0, .nonfoliar = 0, .woody = 0, .standing_dead = 0.25 },
    }, 0, 2);
    try std.testing.expectApproxEqAbs(@as(f64, 2), removed_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 20.0 / 3.0), state.plant_standing_dead_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 10.0 / 3.0), state.plant_charcoal_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 25.0 / 6.0), layers.plant_standing_dead_area_m2[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), ledgers[0].standing_dead_export.carbon_g, 1e-12);
    var manure_carbon_g_c: f64 = 0;
    for (ledgers[0].manure.organic_by_biochemical_fraction) |mass| manure_carbon_g_c += mass.carbon_g;
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), manure_carbon_g_c, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 12), state.plant_standing_dead_carbon_g[0] +
        state.plant_charcoal_carbon_g[0] + productLedgerCarbonG(ledgers[0]), 1e-12);
}

test "GROSUB tillage removes shoot and standing dead but leaves root commit separate" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{1}, &.{0});
    defer state.deinit();
    var layers = try canopy_layers.State.init(std.testing.allocator, 1, 1, 1, 1, 1, &state);
    defer layers.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    var roots = try root_system.State.init(std.testing.allocator, 1, 1, 1);
    defer roots.deinit();
    const root = try roots.layerIndex(0, 0, 0);
    roots.mobile_carbon_g[root] = 7;
    state.branch_stalk_carbon_g[0] = 10;
    state.branch_stalk_nitrogen_g[0] = 1;
    state.branch_stalk_phosphorus_g[0] = 0.1;
    state.node_height_m[0] = 1;
    state.node_internode_length_m[0] = 1;
    state.branch_symbiont_mobile_carbon_g[0] = 2;
    state.branch_symbiont_structural_carbon_g[0] = 3;
    state.branch_grain_carbon_g[0] = 4;
    state.plant_seed_storage_carbon_g[0] = 2;
    state.plant_population_per_m2[0] = 20;
    state.plant_population_count[0] = 20;
    state.plant_standing_dead_population_count[0] = 4;
    state.plant_standing_dead_carbon_g[0] = 8;
    state.plant_standing_dead_nitrogen_g[0] = 0.8;
    state.plant_standing_dead_phosphorus_g[0] = 0.08;
    for (0..4) |kinetic| {
        state.plant_standing_dead_carbon_by_kinetic_g[kinetic] = 2;
        state.plant_standing_dead_nitrogen_by_kinetic_g[kinetic] = 0.2;
        state.plant_standing_dead_phosphorus_by_kinetic_g[kinetic] = 0.02;
    }
    layers.plant_standing_dead_area_m2[0] = 5;
    layers.cell_standing_dead_area_m2[0] = 5;
    layers.plant_standing_dead_projected_surface_m2[0] = 5;
    const science = [_]ScienceParameters{.{
        .carbon_woody_fraction = .{ 0.25, 0.75 },
        .leaf_nitrogen_woody_fraction = .{ 0.25, 0.75 },
        .sheath_nitrogen_woody_fraction = .{ 0.25, 0.75 },
        .leaf_phosphorus_woody_fraction = .{ 0.25, 0.75 },
        .sheath_phosphorus_woody_fraction = .{ 0.25, 0.75 },
    }};
    const root_nonwoody = [_]f64{0.4};
    var ledgers = [_]ProductLedger{.{}};
    var context: Context = .{
        .canopy_state = &state,
        .canopy_layer_state = &layers,
        .branch_development = &development,
        .science_by_plant = &science,
        .products_by_plant = &ledgers,
        .leaf_area_presence_tolerance_m2 = 1e-12,
        .root_state = &roots,
        .root_woody_fraction_by_plant = &root_nonwoody,
    };
    try applyAbovegroundTillage(&context, 0, 0.5, true);
    try std.testing.expectApproxEqAbs(@as(f64, 5), state.branch_stalk_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), state.branch_symbiont_mobile_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), state.branch_symbiont_structural_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2), state.branch_grain_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2), state.plant_seed_storage_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 4), state.plant_standing_dead_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2.5), layers.plant_standing_dead_area_m2[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 10), state.plant_population_count[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 7), roots.mobile_carbon_g[root], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 13.5), productLedgerCarbonG(ledgers[0]), 1e-12);
}

test "GROSUB kind zero standing dead thinning separates retained litter and export" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{0}, &.{});
    defer state.deinit();
    var layers = try canopy_layers.State.init(std.testing.allocator, 1, 1, 1, 1, 1, &state);
    defer layers.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    state.plant_standing_dead_carbon_g[0] = 10;
    state.plant_standing_dead_nitrogen_g[0] = 1;
    state.plant_standing_dead_phosphorus_g[0] = 0.1;
    state.plant_charcoal_carbon_g[0] = 5;
    state.plant_charcoal_nitrogen_g[0] = 0.5;
    state.plant_charcoal_phosphorus_g[0] = 0.05;
    for (0..4) |kinetic| state.plant_standing_dead_carbon_by_kinetic_g[kinetic] = 2.5;
    layers.plant_standing_dead_area_m2[0] = 5;
    layers.cell_standing_dead_area_m2[0] = 5;
    layers.plant_standing_dead_projected_surface_m2[0] = 5;
    const root_nonwoody = [_]f64{0.25};
    var ledgers = [_]ProductLedger{.{}};
    var context: Context = .{
        .canopy_state = &state,
        .canopy_layer_state = &layers,
        .branch_development = &development,
        .science_by_plant = &.{},
        .products_by_plant = &ledgers,
        .leaf_area_presence_tolerance_m2 = 1e-12,
        .root_woody_fraction_by_plant = &root_nonwoody,
    };
    try applyScheduledStandingDeadHarvest(&context, 0, .{
        .date = .{ .day = 1, .month = 1, .year = 9999 },
        .kind = .none,
        .termination = .retain,
        .cutting_height_m_or_lai_fraction = 0,
        .thinning_fraction_or_consumption_rate = 0.4,
        .harvested_fraction = .{ .leaf = 0, .nonfoliar = 0, .woody = 0, .standing_dead = 0.5 },
        .ecosystem_export_fraction = .{ .leaf = 0, .nonfoliar = 0, .woody = 0, .standing_dead = 0.5 },
    });
    try std.testing.expectApproxEqAbs(@as(f64, 6), state.plant_standing_dead_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), state.plant_charcoal_carbon_g[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), ledgers[0].standing_dead_export.carbon_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), ledgers[0].woody.litter.carbon_g + ledgers[0].nonfoliar.litter.carbon_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), ledgers[0].standing_dead_charcoal_litter.carbon_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 15), state.plant_standing_dead_carbon_g[0] +
        state.plant_charcoal_carbon_g[0] + productLedgerCarbonG(ledgers[0]), 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), layers.plant_standing_dead_area_m2[0], 1e-12);
}

test "grazing publication commits manure mineral nutrients and export without examples" {
    var state = try canopy.State.init(std.testing.allocator, 1, 1, &.{1}, &.{0}, &.{});
    defer state.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    var partitions = try litter_partition.State.init(std.testing.allocator, 1);
    defer partitions.deinit();
    const uniform: litter_partition.ElementFractions = .{
        .carbon = .{ 0.25, 0.25, 0.25, 0.25 },
        .nitrogen = .{ 0.25, 0.25, 0.25, 0.25 },
        .phosphorus = .{ 0.25, 0.25, 0.25, 0.25 },
    };
    @memset(partitions.by_plant_and_organ, uniform);
    const config = try @import("config.zig").SimulationConfig.init(
        .{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-12, .max_nonlinear_iterations = 10 },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var surface = try soil_organic.State.init(std.testing.allocator, 1);
    defer surface.deinit();
    var nutrients = try surface_nutrients.State.init(std.testing.allocator, 1, soil_organic.substrate_count);
    defer nutrients.deinit();
    var daily_manure_c = [_]f64{0};
    var daily_manure_n = [_]f64{0};
    var daily_manure_p = [_]f64{0};
    var shoot_litter_carbon_g_c = [_]f64{0};
    var shoot_litter_nitrogen_g_n = [_]f64{0};
    var shoot_litter_phosphorus_g_p = [_]f64{0};
    var hourly_manure_products = [_]grazing_manure.Products{.{}};
    const science = [_]ScienceParameters{.{ .carbon_woody_fraction = .{ 0, 1 }, .leaf_nitrogen_woody_fraction = .{ 0, 1 }, .sheath_nitrogen_woody_fraction = .{ 0, 1 }, .leaf_phosphorus_woody_fraction = .{ 0, 1 }, .sheath_phosphorus_woody_fraction = .{ 0, 1 } }};
    var ledgers = [_]ProductLedger{.{}};
    ledgers[0].direct_litter.nonwoody_carbon_g[0] = 2;
    ledgers[0].direct_litter.nonwoody_nitrogen_g[0] = 0.2;
    ledgers[0].direct_litter.nonwoody_phosphorus_g[0] = 0.02;
    ledgers[0].manure = try grazing_manure.partition(.animal_grazing, .{ .carbon_g = 10, .nitrogen_g = 1, .phosphorus_g = 0.2 });
    const expected_hourly_manure = ledgers[0].manure;
    ledgers[0].standing_dead_export = .{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 };
    var context: Context = .{
        .canopy_state = &state,
        .branch_development = &development,
        .science_by_plant = &science,
        .products_by_plant = &ledgers,
        .root_litter_partition = &partitions,
        .surface_organic_state = &surface,
        .surface_nutrient_state = &nutrients,
        .daily_manure_carbon_input_g_c = &daily_manure_c,
        .daily_manure_nitrogen_input_g_n = &daily_manure_n,
        .daily_manure_phosphorus_input_g_p = &daily_manure_p,
        .hourly_manure_products_by_plant = &hourly_manure_products,
        .shoot_litter_carbon_g_c_by_plant = &shoot_litter_carbon_g_c,
        .shoot_litter_nitrogen_g_n_by_plant = &shoot_litter_nitrogen_g_n,
        .shoot_litter_phosphorus_g_p_by_plant = &shoot_litter_phosphorus_g_p,
        .grid = &grid,
        .leaf_area_presence_tolerance_m2 = 1e-12,
    };
    const exported = try publishPlantProducts(&context, 0);
    try std.testing.expectApproxEqAbs(@as(f64, 1), exported.carbon_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2), shoot_litter_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), shoot_litter_nitrogen_g_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.02), shoot_litter_phosphorus_g_p[0], 1e-12);
    var manure: canopy.ElementalMass = .{};
    for (0..4) |fraction| {
        const index = (2 * soil_organic.structural_fraction_count) + fraction;
        manure.carbon_g += surface.structural[index].carbon_g_c;
        manure.nitrogen_g += surface.structural[index].nitrogen_g_n;
        manure.phosphorus_g += surface.structural[index].phosphorus_g_p;
    }
    try std.testing.expectApproxEqAbs(@as(f64, 10), manure.carbon_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), manure.nitrogen_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1), manure.phosphorus_g, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5 / 14.0), nutrients.pending_surface_ammonium_mol_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.1 / 31.0), nutrients.pending_surface_phosphate_mol_p[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 10), daily_manure_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), daily_manure_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.2), daily_manure_p[0], 1e-12);
    try std.testing.expectEqual(
        expected_hourly_manure,
        hourly_manure_products[0],
    );
    try std.testing.expectEqual(ProductLedger{}, ledgers[0]);

    const surface_before = try std.testing.allocator.dupe(soil_organic.ElementPool, surface.structural);
    defer std.testing.allocator.free(surface_before);
    shoot_litter_phosphorus_g_p[0] = std.math.floatMax(f64);
    ledgers[0].direct_litter.nonwoody_phosphorus_g[0] = std.math.floatMax(f64);
    try std.testing.expectError(
        error.InvalidPlantHarvestProduct,
        publishPlantProducts(&context, 0),
    );
    try std.testing.expectEqualSlices(soil_organic.ElementPool, surface_before, surface.structural);
    try std.testing.expectEqual(
        std.math.floatMax(f64),
        ledgers[0].direct_litter.nonwoody_phosphorus_g[0],
    );
}

test "thinning then complete mortality conserves host and nodule roots and publishes HCNET disturbance" {
    const branch_counts = [_]usize{1};
    const node_counts = [_]usize{1};
    const sample_counts = [_]usize{1};
    var canopy_state = try canopy.State.init(std.testing.allocator, 1, 1, &branch_counts, &node_counts, &sample_counts);
    defer canopy_state.deinit();
    var development = try phenology.BranchDevelopmentState.init(std.testing.allocator, 1);
    defer development.deinit();
    var roots = try root_system.State.init(std.testing.allocator, 1, 2, 1);
    defer roots.deinit();
    var partitions = try litter_partition.State.init(std.testing.allocator, 1);
    defer partitions.deinit();
    const uniform: litter_partition.ElementFractions = .{
        .carbon = .{ 0.25, 0.25, 0.25, 0.25 },
        .nitrogen = .{ 0.25, 0.25, 0.25, 0.25 },
        .phosphorus = .{ 0.25, 0.25, 0.25, 0.25 },
    };
    partitions.by_plant_and_organ[@intFromEnum(litter_partition.Organ.fine_root)] = uniform;
    partitions.by_plant_and_organ[@intFromEnum(litter_partition.Organ.nonstructural)] = uniform;
    const config = try @import("config.zig").SimulationConfig.init(
        .{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 2, .plant_populations = 1 },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-12, .max_nonlinear_iterations = 10 },
    );
    var grid = try grid_module.GridState.init(std.testing.allocator, config);
    defer grid.deinit();
    var organic = try soil_organic.State.init(std.testing.allocator, 2);
    defer organic.deinit();
    for (0..2) |layer| {
        const root = try roots.layerIndex(0, 0, layer);
        roots.symbiont_structural_carbon_g_c[root] = 2;
        roots.symbiont_structural_nitrogen_g_n[root] = 0.2;
        roots.symbiont_structural_phosphorus_g_p[root] = 0.02;
        roots.symbiont_mobile_carbon_g_c[root] = 1;
        roots.symbiont_mobile_nitrogen_g_n[root] = 0.1;
        roots.symbiont_mobile_phosphorus_g_p[root] = 0.01;
        roots.mobile_carbon_g[root] = 1;
        roots.mobile_nitrogen_g[root] = 0.1;
        roots.mobile_phosphorus_g[root] = 0.01;
        const axis_layer = try roots.layerAxisIndex(0, 0, layer, 0);
        roots.axis_primary_carbon_g[axis_layer] = 2;
        roots.axis_primary_nitrogen_g[axis_layer] = 0.2;
        roots.axis_primary_phosphorus_g[axis_layer] = 0.02;
        roots.axis_secondary_carbon_g[axis_layer] = 1;
        roots.axis_secondary_nitrogen_g[axis_layer] = 0.1;
        roots.axis_secondary_phosphorus_g[axis_layer] = 0.01;
        for (0..root_system.biological_domain_count) |domain| {
            const gas_root = try roots.layerIndex(0, domain, layer);
            roots.gaseous_carbon_dioxide_g_c[gas_root] = 1;
            roots.aqueous_carbon_dioxide_g_c[gas_root] = 3;
        }
    }
    roots.active_root_axis_count[0] = 1;
    var exchange = try carbon_exchange.State.init(std.testing.allocator, 1);
    defer exchange.deinit();
    var root_litter_carbon = try root_litter_ledger.State.init(
        std.testing.allocator,
        1,
        root_system.biological_domain_count,
        2,
    );
    defer root_litter_carbon.deinit();
    const woody_fraction = [_]f64{0};
    const science = [_]ScienceParameters{.{ .nitrogen_fixation_type = 1, .carbon_woody_fraction = .{ 0, 1 }, .leaf_nitrogen_woody_fraction = .{ 0, 1 }, .sheath_nitrogen_woody_fraction = .{ 0, 1 }, .leaf_phosphorus_woody_fraction = .{ 0, 1 }, .sheath_phosphorus_woody_fraction = .{ 0, 1 } }};
    var products = [_]ProductLedger{.{}};
    var context: Context = .{
        .canopy_state = &canopy_state,
        .branch_development = &development,
        .science_by_plant = &science,
        .products_by_plant = &products,
        .leaf_area_presence_tolerance_m2 = 1.0e-12,
        .root_state = &roots,
        .root_litter_partition = &partitions,
        .soil_organic_state = &organic,
        .grid = &grid,
        .root_woody_fraction_by_plant = &woody_fraction,
        .carbon_exchange_state = &exchange,
        .root_litter_carbon_ledger = &root_litter_carbon,
    };
    try applyRootSymbiontHarvest(&context, 0, 0.5);
    for (0..2) |layer| {
        try std.testing.expectApproxEqAbs(
            3.5,
            root_litter_carbon.carbon_g_c[try root_litter_carbon.index(0, 0, layer)],
            1e-14,
        );
        try std.testing.expectApproxEqAbs(
            0,
            root_litter_carbon.carbon_g_c[try root_litter_carbon.index(0, 1, layer)],
            1e-14,
        );
    }
    var remaining_carbon_g_c: f64 = 0;
    var litter_carbon_g_c: f64 = 0;
    for (0..2) |layer| {
        const root = try roots.layerIndex(0, 0, layer);
        const axis_layer = try roots.layerAxisIndex(0, 0, layer, 0);
        remaining_carbon_g_c += roots.symbiont_structural_carbon_g_c[root] + roots.symbiont_mobile_carbon_g_c[root] +
            roots.mobile_carbon_g[root] + roots.axis_primary_carbon_g[axis_layer] + roots.axis_secondary_carbon_g[axis_layer];
        for (0..root_litterfall.kinetic_component_count) |component|
            litter_carbon_g_c += organic.structural[(layer * soil_organic.substrate_count + 1) * soil_organic.structural_fraction_count + component].carbon_g_c;
    }
    try std.testing.expectApproxEqAbs(14, remaining_carbon_g_c + litter_carbon_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(7, remaining_carbon_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(4, exchange.disturbance_carbon_g_c_per_h[0], 1e-14);
    try std.testing.expectApproxEqAbs(-8, roots.withdrawal_carbon_dioxide_loss_g_c_per_h[0], 1e-14);
    var remaining_root_carbon_dioxide_g_c: f64 = 0;
    for (roots.gaseous_carbon_dioxide_g_c, roots.aqueous_carbon_dioxide_g_c) |gaseous, aqueous|
        remaining_root_carbon_dioxide_g_c += gaseous + aqueous;
    try std.testing.expectApproxEqAbs(8, remaining_root_carbon_dioxide_g_c, 1e-14);

    try releaseDeadRootsToLitter(&context, 0);
    remaining_carbon_g_c = 0;
    litter_carbon_g_c = 0;
    for (0..2) |layer| {
        const root = try roots.layerIndex(0, 0, layer);
        const axis_layer = try roots.layerAxisIndex(0, 0, layer, 0);
        remaining_carbon_g_c += roots.symbiont_structural_carbon_g_c[root] + roots.symbiont_mobile_carbon_g_c[root] +
            roots.mobile_carbon_g[root] + roots.axis_primary_carbon_g[axis_layer] + roots.axis_secondary_carbon_g[axis_layer];
        for (0..root_litterfall.kinetic_component_count) |component|
            litter_carbon_g_c += organic.structural[(layer * soil_organic.substrate_count + 1) * soil_organic.structural_fraction_count + component].carbon_g_c;
    }
    try std.testing.expectApproxEqAbs(0, remaining_carbon_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(14, litter_carbon_g_c, 1e-14);
    try std.testing.expectApproxEqAbs(8, exchange.disturbance_carbon_g_c_per_h[0], 1e-14);
    try std.testing.expectApproxEqAbs(-16, roots.withdrawal_carbon_dioxide_loss_g_c_per_h[0], 1e-14);
    remaining_root_carbon_dioxide_g_c = 0;
    for (roots.gaseous_carbon_dioxide_g_c, roots.aqueous_carbon_dioxide_g_c) |gaseous, aqueous|
        remaining_root_carbon_dioxide_g_c += gaseous + aqueous;
    try std.testing.expectApproxEqAbs(0, remaining_root_carbon_dioxide_g_c, 1e-14);
}

test "GROSUB harvest leaves nodule pools intact for non-fixing plants" {
    const pool: symbiotic_fixation.Pool = .{
        .carbon_g_c = 2,
        .nitrogen_g_n = 0.2,
        .phosphorus_g_p = 0.02,
    };
    const partition: litter_partition.ElementFractions = .{
        .carbon = .{ 0.25, 0.25, 0.25, 0.25 },
        .nitrogen = .{ 0.25, 0.25, 0.25, 0.25 },
        .phosphorus = .{ 0.25, 0.25, 0.25, 0.25 },
    };
    const result = try noduleHarvestResult(0, pool, pool, 0.5, partition, partition);
    try std.testing.expectEqual(pool, result.structural);
    try std.testing.expectEqual(pool, result.mobile);
    try std.testing.expectEqual(
        @as(f64, 0),
        try root_litter_ledger.totalCarbon(result.litterfall),
    );
}

test "source-order tillage population reduction scales all population fields" {
    const result = try sourceOrderTillagePopulationReduction(.{
        .hour_of_day = 12,
        .local_solar_noon_h = 12.75,
        .biomass_turnover_type = 0,
        .root_profile_type = 2,
        .current_day_of_year = 151,
        .current_year = 2001,
        .planting_day_of_year = 100,
        .planting_year = 2001,
        .tillage_code = 8,
        .is_first_plant_population = true,
        .remaining_fraction = 0.25,
        .zero_population_threshold = 1.0e-12,
        .state = .{
            .living_population_per_m2 = 8,
            .living_population_count = 40,
            .standing_dead_population_count = 12,
            .canopy_radiation_fraction = 0.8,
        },
    });
    try std.testing.expect(result.applied);
    try std.testing.expect(result.clear_leaf_sheath_and_sapwood_totals);
    try std.testing.expect(!result.terminate_living_branches);
    try std.testing.expectEqual(@as(f64, 2), result.state.living_population_per_m2);
    try std.testing.expectEqual(@as(f64, 10), result.state.living_population_count);
    try std.testing.expectEqual(@as(f64, 3), result.state.standing_dead_population_count);
    try std.testing.expectEqual(@as(f64, 0.2), result.state.canopy_radiation_fraction);
}

test "source-order tillage population selector preserves crop and date gates" {
    const base: SourceOrderTillagePopulationInput = .{
        .hour_of_day = 11,
        .local_solar_noon_h = 11.4,
        .biomass_turnover_type = 1,
        .root_profile_type = 1,
        .current_day_of_year = 200,
        .current_year = 2000,
        .planting_day_of_year = 100,
        .planting_year = 2001,
        .tillage_code = 15,
        .is_first_plant_population = true,
        .remaining_fraction = 0,
        .zero_population_threshold = 1.0e-9,
        .state = .{
            .living_population_per_m2 = 1,
            .living_population_count = 1,
            .standing_dead_population_count = 1,
            .canopy_radiation_fraction = 1,
        },
    };
    try std.testing.expect(!(try sourceOrderTillagePopulationReduction(base)).applied);

    var second_population = base;
    second_population.is_first_plant_population = false;
    const source_date_result = try sourceOrderTillagePopulationReduction(second_population);
    try std.testing.expect(source_date_result.applied);
    try std.testing.expect(source_date_result.terminate_living_branches);

    var planting_date = second_population;
    planting_date.current_day_of_year = planting_date.planting_day_of_year;
    planting_date.current_year = planting_date.planting_year;
    try std.testing.expect(!(try sourceOrderTillagePopulationReduction(planting_date)).applied);
}

test "source-order tillage branch litter conserves every element" {
    const kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 0.1, 0.2, 0.3, 0.4 },
        .nitrogen = .{ 0.4, 0.3, 0.2, 0.1 },
        .phosphorus = .{ 0.25, 0.25, 0.25, 0.25 },
    };
    const composition: TillageElementComposition = .{
        .carbon = .{ 0.25, 0.75 },
        .nitrogen = .{ 0.4, 0.6 },
        .phosphorus = .{ 0.5, 0.5 },
    };
    const pools: SourceOrderTillageBranchPools = .{
        .host_mobile = .{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 },
        .symbiont_mobile = .{ .carbon_g = 2, .nitrogen_g = 0.2, .phosphorus_g = 0.02 },
        .c4_mobile_carbon_g_c = 0.5,
        .stalk_reserve = .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 },
        .leaf = .{ .carbon_g = 4, .nitrogen_g = 0.4, .phosphorus_g = 0.04 },
        .symbiont_structural = .{ .carbon_g = 5, .nitrogen_g = 0.5, .phosphorus_g = 0.05 },
        .sheath = .{ .carbon_g = 6, .nitrogen_g = 0.6, .phosphorus_g = 0.06 },
        .husk = .{ .carbon_g = 7, .nitrogen_g = 0.7, .phosphorus_g = 0.07 },
        .ear = .{ .carbon_g = 8, .nitrogen_g = 0.8, .phosphorus_g = 0.08 },
        .grain = .{ .carbon_g = 9, .nitrogen_g = 0.9, .phosphorus_g = 0.09 },
        .stalk = .{ .carbon_g = 10, .nitrogen_g = 1, .phosphorus_g = 0.1 },
    };
    const result = try sourceOrderTillageBranchLitter(.{
        .remaining_fraction = 0.4,
        .winter_annual = true,
        .pools = pools,
        .leaf_composition = composition,
        .sheath_composition = composition,
        .stalk_composition = composition,
        .nonstructural_kinetics = kinetics,
        .foliar_kinetics = kinetics,
        .nonfoliar_kinetics = kinetics,
        .stalk_kinetics = kinetics,
        .coarse_wood_kinetics = kinetics,
    });
    var litter_carbon_g_c: f64 = 0;
    var litter_nitrogen_g_n: f64 = 0;
    var litter_phosphorus_g_p: f64 = 0;
    for (0..litter_partition.kinetic_component_count) |kinetic| {
        litter_carbon_g_c += result.litter.woody_carbon_g[kinetic] + result.litter.nonwoody_carbon_g[kinetic];
        litter_nitrogen_g_n += result.litter.woody_nitrogen_g[kinetic] + result.litter.nonwoody_nitrogen_g[kinetic];
        litter_phosphorus_g_p += result.litter.woody_phosphorus_g[kinetic] + result.litter.nonwoody_phosphorus_g[kinetic];
    }
    try std.testing.expectApproxEqAbs(0.6 * 55.5, litter_carbon_g_c + result.seasonal_storage.carbon_g, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.6 * 5.5, litter_nitrogen_g_n + result.seasonal_storage.nitrogen_g, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.6 * 0.55, litter_phosphorus_g_p + result.seasonal_storage.phosphorus_g, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.6 * pools.grain.carbon_g, result.seasonal_storage.carbon_g, 1.0e-12);
}

test "source-order tillage routes non-winter grain to nonwoody litter" {
    const one_hot: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const zero_mass: canopy.ElementalMass = .{};
    const result = try sourceOrderTillageBranchLitter(.{
        .remaining_fraction = 0.25,
        .winter_annual = false,
        .pools = .{
            .host_mobile = zero_mass,
            .symbiont_mobile = zero_mass,
            .c4_mobile_carbon_g_c = 0,
            .stalk_reserve = zero_mass,
            .leaf = zero_mass,
            .symbiont_structural = zero_mass,
            .sheath = zero_mass,
            .husk = zero_mass,
            .ear = zero_mass,
            .grain = .{ .carbon_g = 4, .nitrogen_g = 0.4, .phosphorus_g = 0.04 },
            .stalk = zero_mass,
        },
        .leaf_composition = .{ .carbon = .{ 0, 1 }, .nitrogen = .{ 0, 1 }, .phosphorus = .{ 0, 1 } },
        .sheath_composition = .{ .carbon = .{ 0, 1 }, .nitrogen = .{ 0, 1 }, .phosphorus = .{ 0, 1 } },
        .stalk_composition = .{ .carbon = .{ 0, 1 }, .nitrogen = .{ 0, 1 }, .phosphorus = .{ 0, 1 } },
        .nonstructural_kinetics = one_hot,
        .foliar_kinetics = one_hot,
        .nonfoliar_kinetics = one_hot,
        .stalk_kinetics = one_hot,
        .coarse_wood_kinetics = one_hot,
    });
    try std.testing.expectEqual(@as(f64, 3), result.litter.nonwoody_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), result.seasonal_storage.carbon_g);
}

test "source-order tillage retains complete branch state over runtime extents" {
    const mass: canopy.ElementalMass = .{ .carbon_g = 8, .nitrogen_g = 4, .phosphorus_g = 2 };
    var scalar: SourceOrderTillageBranchScalarState = .{
        .host_mobile = mass,
        .c4_mobile_carbon_g_c = 8,
        .symbiont_mobile = mass,
        .total_shoot = mass,
        .leaf = mass,
        .symbiont_structural = mass,
        .sheath = mass,
        .stalk = mass,
        .sapwood_carbon_g_c = 8,
        .stalk_reserve = mass,
        .husk = mass,
        .ear = mass,
        .grain = mass,
        .potential_seed_site_count = 8,
        .seed_count = 8,
        .individual_seed_carbon_g_c = 3,
        .leaf_area_m2 = 8,
        .stalk_total = mass,
    };
    var node_values: [16][3]f64 = @splat(.{ 2, 4, 6 });
    var sample_values: [4][2]f64 = @splat(.{ 10, 12 });
    const result = try sourceOrderRetainTillageBranchState(&scalar, .{
        .c3_mobile_carbon_g_c = &node_values[0],
        .c4_mobile_carbon_g_c = &node_values[1],
        .carbon_dioxide_g_c = &node_values[2],
        .bicarbonate_g_c = &node_values[3],
        .leaf_area_m2 = &node_values[4],
        .growing_leaf_carbon_g_c = &node_values[5],
        .senescing_leaf_carbon_g_c = &node_values[6],
        .growing_sheath_carbon_g_c = &node_values[7],
        .senescing_sheath_carbon_g_c = &node_values[8],
        .growing_node_carbon_g_c = &node_values[9],
        .growing_leaf_nitrogen_g_n = &node_values[10],
        .growing_sheath_nitrogen_g_n = &node_values[11],
        .growing_node_nitrogen_g_n = &node_values[12],
        .growing_leaf_phosphorus_g_p = &node_values[13],
        .growing_sheath_phosphorus_g_p = &node_values[14],
        .growing_node_phosphorus_g_p = &node_values[15],
    }, .{
        .leaf_area_m2 = &sample_values[0],
        .growing_leaf_carbon_g_c = &sample_values[1],
        .growing_leaf_nitrogen_g_n = &sample_values[2],
        .growing_leaf_phosphorus_g_p = &sample_values[3],
    }, 0.25);

    try std.testing.expectEqual(@as(f64, 2), scalar.host_mobile.carbon_g);
    try std.testing.expectEqual(@as(f64, 3), scalar.individual_seed_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 4), result.leaf_sheath_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 2), result.sapwood_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 2), node_values[0][0]);
    try std.testing.expectEqual(@as(f64, 1), node_values[0][1]);
    try std.testing.expectEqual(@as(f64, 0.5), node_values[4][0]);
    try std.testing.expectEqual(@as(f64, 2.5), sample_values[0][0]);
}

test "source-order tillage standing dead repartitions and conserves C N P" {
    const stalk_kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 0.1, 0.2, 0.3, 0.4 },
        .nitrogen = .{ 0.4, 0.3, 0.2, 0.1 },
        .phosphorus = .{ 0.25, 0.25, 0.25, 0.25 },
    };
    const coarse_kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 0.4, 0.3, 0.2, 0.1 },
        .nitrogen = .{ 0.1, 0.2, 0.3, 0.4 },
        .phosphorus = .{ 0.25, 0.25, 0.25, 0.25 },
    };
    const result = try sourceOrderTillageStandingDead(.{
        .remaining_fraction = 0.25,
        .standing_dead_by_source_component = .{
            .{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 },
            .{ .carbon_g = 2, .nitrogen_g = 0.2, .phosphorus_g = 0.02 },
            .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 },
            .{ .carbon_g = 4, .nitrogen_g = 0.4, .phosphorus_g = 0.04 },
        },
        .composition = .{
            .carbon = .{ 0.6, 0.4 },
            .nitrogen = .{ 0.7, 0.3 },
            .phosphorus = .{ 0.8, 0.2 },
        },
        .stalk_kinetics = stalk_kinetics,
        .coarse_wood_kinetics = coarse_kinetics,
    });
    var litter_carbon_g_c: f64 = 0;
    var litter_nitrogen_g_n: f64 = 0;
    var litter_phosphorus_g_p: f64 = 0;
    for (0..litter_partition.kinetic_component_count) |kinetic| {
        litter_carbon_g_c += result.litter.woody_carbon_g[kinetic] + result.litter.nonwoody_carbon_g[kinetic];
        litter_nitrogen_g_n += result.litter.woody_nitrogen_g[kinetic] + result.litter.nonwoody_nitrogen_g[kinetic];
        litter_phosphorus_g_p += result.litter.woody_phosphorus_g[kinetic] + result.litter.nonwoody_phosphorus_g[kinetic];
    }
    try std.testing.expectApproxEqAbs(7.5, litter_carbon_g_c, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.75, litter_nitrogen_g_n, 1.0e-12);
    try std.testing.expectApproxEqAbs(0.075, litter_phosphorus_g_p, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 1), result.remaining_by_source_component[3].carbon_g);
    try std.testing.expectApproxEqAbs(
        0.75 * 0.4 * 10 * 0.6,
        result.litter.woody_carbon_g[0],
        1.0e-12,
    );
}

test "source-order tillage termination sets all death flags and harvest date" {
    const prior: SourceOrderTillageTerminationState = .{
        .roots_dead = false,
        .shoots_dead = false,
        .plant_dead = false,
        .harvest_termination_code = 0,
        .harvest_day_of_year = 0,
        .harvest_year = 0,
    };
    const alive = try sourceOrderTillageTermination(.{
        .living_population_count = 1.1e-6,
        .zero_population_threshold = 1.0e-6,
        .current_day_of_year = 240,
        .current_year = 2004,
        .state = prior,
    });
    try std.testing.expect(!alive.terminated);
    try std.testing.expectEqualDeep(prior, alive.state);

    const terminated = try sourceOrderTillageTermination(.{
        .living_population_count = 1.0e-6,
        .zero_population_threshold = 1.0e-6,
        .current_day_of_year = 240,
        .current_year = 2004,
        .state = prior,
    });
    try std.testing.expect(terminated.terminated);
    try std.testing.expect(terminated.state.roots_dead);
    try std.testing.expect(terminated.state.shoots_dead);
    try std.testing.expect(terminated.state.plant_dead);
    try std.testing.expectEqual(@as(u8, 1), terminated.state.harvest_termination_code);
    try std.testing.expectEqual(@as(u16, 240), terminated.state.harvest_day_of_year);
    try std.testing.expectEqual(@as(u32, 2004), terminated.state.harvest_year);
}

test "source-order standing dead harvest separates harvest litter and retention" {
    const result = try sourceOrderStandingDeadHarvest(.{
        .harvest_code = 0,
        .hour_of_day = 12,
        .local_solar_noon_h = 12.8,
        .thinning_fraction_or_specific_consumption_rate = 0.5,
        .standing_dead_removal_fraction = 0.4,
        .grazer_live_mass_g_per_m2 = 0,
        .animal_accessible_area_m2 = 0,
        .insect_accessible_area_m2 = 0,
        .standing_dead_presence_threshold_g_c = 1.0e-12,
        .standing_dead_by_component = .{
            .{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 },
            .{ .carbon_g = 2, .nitrogen_g = 0.2, .phosphorus_g = 0.02 },
            .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 },
            .{ .carbon_g = 4, .nitrogen_g = 0.4, .phosphorus_g = 0.04 },
            .{ .carbon_g = 5, .nitrogen_g = 0.5, .phosphorus_g = 0.05 },
        },
    });
    try std.testing.expectEqual(@as(f64, 0.5), result.retained_fraction);
    try std.testing.expectEqual(@as(f64, 0.8), result.harvested_fraction);
    try std.testing.expectApproxEqAbs(3, result.harvested.carbon_g, 1.0e-14);
    try std.testing.expectApproxEqAbs(4.5, result.returned_to_litter.carbon_g, 1.0e-14);
    try std.testing.expectEqual(@as(f64, 2.5), result.remaining_by_component[4].carbon_g);
    try std.testing.expectApproxEqAbs(
        15,
        result.harvested.carbon_g + result.returned_to_litter.carbon_g +
            result.remaining_by_component[0].carbon_g + result.remaining_by_component[1].carbon_g +
            result.remaining_by_component[2].carbon_g + result.remaining_by_component[3].carbon_g +
            result.remaining_by_component[4].carbon_g,
        1.0e-14,
    );
}

test "source-order standing dead grazing selects animal and insect areas" {
    const base: SourceOrderStandingDeadHarvestInput = .{
        .harvest_code = 4,
        .hour_of_day = 3,
        .local_solar_noon_h = 12,
        .thinning_fraction_or_specific_consumption_rate = 1,
        .standing_dead_removal_fraction = 0.5,
        .grazer_live_mass_g_per_m2 = 96,
        .animal_accessible_area_m2 = 2,
        .insect_accessible_area_m2 = 5,
        .standing_dead_presence_threshold_g_c = 1,
        .standing_dead_by_component = .{
            .{ .carbon_g = 2 },
            .{ .carbon_g = 2 },
            .{ .carbon_g = 2 },
            .{ .carbon_g = 2 },
            .{ .carbon_g = 2 },
        },
    };
    const animal = try sourceOrderStandingDeadHarvest(base);
    try std.testing.expectEqual(@as(f64, 0.8), animal.retained_fraction);
    var insect_input = base;
    insect_input.harvest_code = 6;
    const insect = try sourceOrderStandingDeadHarvest(insect_input);
    try std.testing.expectEqual(@as(f64, 0.5), insect.retained_fraction);
    var absent = base;
    absent.standing_dead_presence_threshold_g_c = 10;
    try std.testing.expectEqual(
        @as(f64, 1),
        (try sourceOrderStandingDeadHarvest(absent)).retained_fraction,
    );
}

test "source-order harvest residue routes all five components" {
    const harvested = [harvest_product_component_count]canopy.ElementalMass{
        .{ .carbon_g = 10, .nitrogen_g = 1, .phosphorus_g = 0.1 },
        .{ .carbon_g = 20, .nitrogen_g = 2, .phosphorus_g = 0.2 },
        .{ .carbon_g = 30, .nitrogen_g = 3, .phosphorus_g = 0.3 },
        .{ .carbon_g = 40, .nitrogen_g = 4, .phosphorus_g = 0.4 },
        .{ .carbon_g = 50, .nitrogen_g = 5, .phosphorus_g = 0.5 },
    };
    const residue = try sourceOrderHarvestResidueRouting(.{
        .harvest_code = 0,
        .harvested_by_component = harvested,
        .harvested_grain = .{},
        .ecosystem_export_fraction = .{ 0.1, 0.2, 0.3, 0.4 },
    });
    try std.testing.expectEqual(@as(f64, 9), residue[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 18), residue[1].carbon_g);
    try std.testing.expectEqual(@as(f64, 24), residue[2].carbon_g);
    try std.testing.expectEqual(@as(f64, 28), residue[3].carbon_g);
    try std.testing.expectEqual(@as(f64, 30), residue[4].carbon_g);

    const grazing = try sourceOrderHarvestResidueRouting(.{
        .harvest_code = 6,
        .harvested_by_component = harvested,
        .harvested_grain = .{},
        .ecosystem_export_fraction = .{ 0.1, 0.2, 0.3, 0.4 },
    });
    try std.testing.expectEqualDeep(residue, grazing);
}

test "source-order grain harvest subtracts only exported grain from component two" {
    const harvested = [harvest_product_component_count]canopy.ElementalMass{
        .{ .carbon_g = 10, .nitrogen_g = 1, .phosphorus_g = 0.1 },
        .{ .carbon_g = 20, .nitrogen_g = 2, .phosphorus_g = 0.2 },
        .{ .carbon_g = 30, .nitrogen_g = 3, .phosphorus_g = 0.3 },
        .{ .carbon_g = 40, .nitrogen_g = 4, .phosphorus_g = 0.4 },
        .{ .carbon_g = 50, .nitrogen_g = 5, .phosphorus_g = 0.5 },
    };
    const residue = try sourceOrderHarvestResidueRouting(.{
        .harvest_code = 1,
        .harvested_by_component = harvested,
        .harvested_grain = .{ .carbon_g = 10, .nitrogen_g = 1, .phosphorus_g = 0.1 },
        .ecosystem_export_fraction = .{ 0.9, 0.5, 0.8, 0.7 },
    });
    try std.testing.expectEqualDeep(harvested[0], residue[0]);
    try std.testing.expectEqualDeep(harvested[1], residue[1]);
    try std.testing.expectEqual(@as(f64, 25), residue[2].carbon_g);
    try std.testing.expectEqualDeep(harvested[3], residue[3]);
    try std.testing.expectEqualDeep(harvested[4], residue[4]);
}

test "source-order disturbance totals route ordinary export and reseed storage" {
    const harvested: [harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 2, .nitrogen_g = 1, .phosphorus_g = 0.5 });
    const residue: [harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 1, .nitrogen_g = 0.5, .phosphorus_g = 0.25 });
    const direct_litter: [harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 0.2, .nitrogen_g = 0.1, .phosphorus_g = 0.05 });
    const ordinary = try sourceOrderTotalDisturbanceRemoval(.{
        .harvest_code = 2,
        .terminate_and_reseed = false,
        .grazer_growth_yield = 0,
        .grazer_respiration_fraction = 0,
        .harvested_by_component = harvested,
        .residue_by_component = residue,
        .direct_litter_by_component = direct_litter,
    });
    try std.testing.expectEqual(@as(f64, 10), ordinary.harvested_total.carbon_g);
    try std.testing.expectEqual(@as(f64, 5), ordinary.residue_total.carbon_g);
    try std.testing.expectEqual(@as(f64, 1), ordinary.direct_litter_total.carbon_g);
    try std.testing.expectEqual(@as(f64, 5), ordinary.plant_ecosystem_removal.carbon_g);
    try std.testing.expectEqual(@as(f64, -5), ordinary.net_biome_production_carbon_change_g_c_per_h);

    const reseed = try sourceOrderTotalDisturbanceRemoval(.{
        .harvest_code = 2,
        .terminate_and_reseed = true,
        .grazer_growth_yield = 0,
        .grazer_respiration_fraction = 0,
        .harvested_by_component = harvested,
        .residue_by_component = residue,
        .direct_litter_by_component = direct_litter,
    });
    try std.testing.expectEqual(@as(f64, 5), reseed.reseed_storage_addition.carbon_g);
    try std.testing.expectEqual(@as(f64, 0), reseed.plant_ecosystem_removal.carbon_g);
    try std.testing.expectEqual(@as(f64, 0), reseed.grid_ecosystem_removal.carbon_g);
}

test "source-order grazing totals split growth and respiration carbon" {
    const harvested: [harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 2, .nitrogen_g = 1, .phosphorus_g = 0.5 });
    const residue: [harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 1, .nitrogen_g = 0.5, .phosphorus_g = 0.25 });
    const result = try sourceOrderTotalDisturbanceRemoval(.{
        .harvest_code = 4,
        .terminate_and_reseed = false,
        .grazer_growth_yield = 0.4,
        .grazer_respiration_fraction = 0.6,
        .harvested_by_component = harvested,
        .residue_by_component = residue,
        .direct_litter_by_component = @splat(.{}),
    });
    try std.testing.expectEqual(@as(f64, 2), result.plant_ecosystem_removal.carbon_g);
    try std.testing.expectEqual(@as(f64, 2.5), result.plant_ecosystem_removal.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 1.25), result.plant_ecosystem_removal.phosphorus_g);
    try std.testing.expectEqual(@as(f64, -3), result.plant_total_respiration_change_g_c_per_h);
    try std.testing.expectEqual(@as(f64, -3), result.plant_actual_respiration_change_g_c_per_h);
    try std.testing.expectEqual(@as(f64, -3), result.ecosystem_respiration_change_g_c_per_h);
    try std.testing.expectEqual(@as(f64, -3), result.autotrophic_respiration_change_g_c_per_h);
}

test "source-order aboveground harvest litter preserves herbaceous and woody routing" {
    const residue = [harvest_product_component_count]canopy.ElementalMass{
        .{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 },
        .{ .carbon_g = 2, .nitrogen_g = 0.2, .phosphorus_g = 0.02 },
        .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 },
        .{ .carbon_g = 4, .nitrogen_g = 0.4, .phosphorus_g = 0.04 },
        .{ .carbon_g = 5, .nitrogen_g = 0.5, .phosphorus_g = 0.05 },
    };
    const direct: [harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 1, .nitrogen_g = 0.1, .phosphorus_g = 0.01 });
    const one_hot: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const common: SourceOrderAbovegroundHarvestLitterInput = .{
        .harvest_code = 2,
        .biomass_turnover_type = 0,
        .root_profile_type = 2,
        .residue_by_component = residue,
        .direct_litter_by_component = direct,
        .woody_composition = .{
            .carbon = .{ 0.25, 0.75 },
            .nitrogen = .{ 0.25, 0.75 },
            .phosphorus = .{ 0.25, 0.75 },
        },
        .nonstructural_kinetics = one_hot,
        .foliar_kinetics = one_hot,
        .nonfoliar_kinetics = one_hot,
        .stalk_kinetics = one_hot,
        .coarse_wood_kinetics = one_hot,
    };
    const herbaceous = try sourceOrderAbovegroundHarvestLitter(common);
    try std.testing.expectEqual(@as(f64, 20), herbaceous.litter.nonwoody_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), herbaceous.litter.woody_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 0), herbaceous.standing_dead_addition.woody_carbon_g[0]);

    var woody_input = common;
    woody_input.biomass_turnover_type = 2;
    const woody = try sourceOrderAbovegroundHarvestLitter(woody_input);
    try std.testing.expectEqual(@as(f64, 15.75), woody.litter.nonwoody_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 2.25), woody.litter.woody_carbon_g[0]);
    try std.testing.expectEqual(@as(f64, 2), woody.standing_dead_addition.woody_carbon_g[0]);
}

test "source-order grazing ledgers preserve five components and manure partition" {
    const residue: [harvest_product_component_count]canopy.ElementalMass =
        @splat(.{ .carbon_g = 1, .nitrogen_g = 0.2, .phosphorus_g = 0.04 });
    const direct = residue;
    const current: SourceOrderGrazingLitterLedgerState = .{
        .hourly_litter = .{},
        .cumulative_litter = .{ .carbon_g = 20, .nitrogen_g = 10, .phosphorus_g = 1 },
        .cumulative_aboveground_litter = .{ .carbon_g = 4, .nitrogen_g = 3, .phosphorus_g = 0.3 },
        .surface_litter_carbon_g_c = 7,
        .accumulated_application = .{},
    };
    const animal = try sourceOrderGrazingLitterLedgers(4, residue, direct, current);
    try std.testing.expectEqual(@as(f64, 10), animal.returned_mass.carbon_g);
    try std.testing.expectEqual(@as(f64, 2), animal.returned_mass.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 14), animal.state.cumulative_aboveground_litter.carbon_g);
    try std.testing.expectEqual(@as(f64, 14), animal.state.cumulative_aboveground_litter.nitrogen_g);
    try std.testing.expectApproxEqAbs(1.8, animal.state.cumulative_aboveground_litter.phosphorus_g, 1.0e-14);
    try std.testing.expectEqual(@as(f64, 17), animal.state.surface_litter_carbon_g_c);
    try std.testing.expectApproxEqAbs(0.36, animal.manure.organic_by_biochemical_fraction[0].carbon_g, 1.0e-14);
    try std.testing.expectEqual(@as(f64, 1), animal.manure.inorganic_nitrogen_g_n);
    try std.testing.expectEqual(@as(f64, 0.2), animal.manure.inorganic_phosphorus_g_p);

    const insect = try sourceOrderGrazingLitterLedgers(6, residue, direct, current);
    try std.testing.expectApproxEqAbs(1.38, insect.manure.organic_by_biochemical_fraction[0].carbon_g, 1.0e-14);
}

test "source-order population thresholds preserve runtime scaling and evaluation order" {
    const result = try sourceOrderPopulationScaledNumericalThresholds(
        250,
        50,
        1.0e-15,
        1.0e-6,
    );
    try std.testing.expectEqual(@as(f64, 2.5e-13), result.plant_mass_presence_g);
    try std.testing.expectEqual(@as(f64, 5.0e-15), result.plant_mass_density_g_m2);
    try std.testing.expectEqual(@as(f64, 2.5e-4), result.plant_flux_presence_g_per_step);

    const extinct = try sourceOrderPopulationScaledNumericalThresholds(0, 50, 1.0e-15, 1.0e-6);
    try std.testing.expectEqual(@as(f64, 0), extinct.plant_mass_presence_g);
    try std.testing.expectEqual(@as(f64, 0), extinct.plant_mass_density_g_m2);
    try std.testing.expectEqual(@as(f64, 0), extinct.plant_flux_presence_g_per_step);

    try std.testing.expectError(
        error.InvalidPopulationScaledNumericalThresholdInput,
        sourceOrderPopulationScaledNumericalThresholds(1, 0, 1.0e-15, 1.0e-6),
    );
    try std.testing.expectError(
        error.InvalidPopulationScaledNumericalThresholdInput,
        sourceOrderPopulationScaledNumericalThresholds(1, 1, -1, 1.0e-6),
    );
}

test "source-order dead branch reset preserves selector and runtime branch loop" {
    const stale: SourceOrderDeadBranchPhenologyState = .{
        .dead = true,
        .maturity_group_node_count = 9,
        .initiated_node_count = 8,
        .nodes_at_floral_initiation = 7,
        .nodes_at_anthesis = 6,
        .appeared_leaf_count = 5,
        .leaves_at_floral_initiation = 4,
        .current_leaf_ordinal = 3,
        .current_growing_leaf_ordinal = 2,
        .normalized_vegetative_node_change = 1,
        .normalized_reproductive_node_change = 1,
        .accumulated_leafout_h = 1,
        .accumulated_leafoff_h = 1,
        .lengthening_photoperiod_h = 1,
        .shortening_photoperiod_h = 1,
        .time_since_germination_h = 1,
        .hours_without_grain_fill = 1,
        .carbon_fixation_feedback = 0.5,
        .carbon_fixation_feedback_previous = 0.5,
        .leafout_initialization_enabled = false,
        .emergence_initialization_disabled = false,
        .leafoff_enabled = false,
        .remobilization_enabled = false,
        .hours_after_maturity_h = 1,
        .new_branch_count = 3,
        .stage_day_of_year = @splat(100),
    };
    var branches = [_]SourceOrderDeadBranchPhenologyState{ stale, stale };
    branches[1].dead = false;
    const live_before = branches[1];
    const applied = try sourceOrderResetDeadBranchPhenology(&branches, .{
        .first_living_branch_emergence_day_of_year = 20,
        .perennial_growth_habit = false,
        .current_day_of_year = 200,
        .current_year = 2025,
        .harvest_day_of_year = 200,
        .harvest_year = 2025,
        .hour_of_day = 12,
        .local_solar_noon_h = 12.8,
        .initial_maturity_group_node_count = 4,
        .initial_node_count = 0.25,
    });
    try std.testing.expect(applied);
    try std.testing.expectEqual(@as(f64, 4), branches[0].maturity_group_node_count);
    try std.testing.expectEqual(@as(f64, 0.25), branches[0].initiated_node_count);
    try std.testing.expectEqual(@as(f64, 0.25), branches[0].nodes_at_floral_initiation);
    try std.testing.expectEqual(@as(usize, 1), branches[0].current_leaf_ordinal);
    try std.testing.expectEqual(@as(f64, 1), branches[0].carbon_fixation_feedback);
    try std.testing.expect(branches[0].leafout_initialization_enabled);
    try std.testing.expect(branches[0].emergence_initialization_disabled);
    try std.testing.expectEqual([_]u16{0} ** 10, branches[0].stage_day_of_year);
    try std.testing.expectEqualDeep(live_before, branches[1]);

    var unselected = [_]SourceOrderDeadBranchPhenologyState{stale};
    const skipped = try sourceOrderResetDeadBranchPhenology(&unselected, .{
        .first_living_branch_emergence_day_of_year = 20,
        .perennial_growth_habit = false,
        .current_day_of_year = 199,
        .current_year = 2025,
        .harvest_day_of_year = 200,
        .harvest_year = 2025,
        .hour_of_day = 12,
        .local_solar_noon_h = 12,
        .initial_maturity_group_node_count = 4,
        .initial_node_count = 0.25,
    });
    try std.testing.expect(!skipped);
    try std.testing.expectEqualDeep(stale, unselected[0]);
}

test "source-order dead branch litterfall conserves runtime branch components" {
    const kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const composition: TillageElementComposition = .{
        .carbon = .{ 0.25, 0.75 },
        .nitrogen = .{ 0.25, 0.75 },
        .phosphorus = .{ 0.25, 0.75 },
    };
    const common: SourceOrderDeadBranchLitterInput = .{
        .annual_growth_habit = false,
        .deciduous_phenology = false,
        .pools = .{
            .bacterial_nonstructural = .{ .carbon_g = 1, .nitrogen_g = 1, .phosphorus_g = 1 },
            .bacterial_structural = .{ .carbon_g = 2, .nitrogen_g = 2, .phosphorus_g = 2 },
            .leaf = .{ .carbon_g = 4, .nitrogen_g = 4, .phosphorus_g = 4 },
            .sheath = .{ .carbon_g = 8, .nitrogen_g = 8, .phosphorus_g = 8 },
            .husk = .{ .carbon_g = 16, .nitrogen_g = 16, .phosphorus_g = 16 },
            .ear = .{ .carbon_g = 32, .nitrogen_g = 32, .phosphorus_g = 32 },
            .grain = .{ .carbon_g = 64, .nitrogen_g = 64, .phosphorus_g = 64 },
            .stalk = .{ .carbon_g = 128, .nitrogen_g = 128, .phosphorus_g = 128 },
        },
        .leaf_woody_fraction = composition,
        .sheath_woody_fraction = composition,
        .nonstructural_kinetics = kinetics,
        .foliar_kinetics = kinetics,
        .nonfoliar_kinetics = kinetics,
        .stalk_kinetics = kinetics,
        .coarse_wood_kinetics = kinetics,
    };
    const ordinary = try sourceOrderDeadBranchLitterfall(common);
    try std.testing.expectEqual(@as(f64, 124), ordinary.nonwoody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 3), ordinary.woody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 128), ordinary.standing_dead_stalk[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 0), ordinary.seasonal_storage_addition.carbon_g);

    var winter_input = common;
    winter_input.annual_growth_habit = true;
    winter_input.deciduous_phenology = true;
    const winter = try sourceOrderDeadBranchLitterfall(winter_input);
    try std.testing.expectEqual(@as(f64, 60), winter.nonwoody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 64), winter.seasonal_storage_addition.carbon_g);
    try std.testing.expectEqual(@as(f64, 255), winter.nonwoody_litter[0].carbon_g +
        winter.woody_litter[0].carbon_g +
        winter.standing_dead_stalk[0].carbon_g +
        winter.seasonal_storage_addition.carbon_g);
}

test "source-order dead branch storage recovery preserves six assignments" {
    const result = try sourceOrderDeadBranchStorageRecovery(.{
        .current_seasonal_storage = .{
            .carbon_g = 100,
            .nitrogen_g = 20,
            .phosphorus_g = 5,
        },
        .branch_mobile = .{
            .carbon_g = 7,
            .nitrogen_g = 3,
            .phosphorus_g = 0.4,
        },
        .c4_intermediate_carbon_g_c = 11,
        .stalk_reserve = .{
            .carbon_g = 13,
            .nitrogen_g = 2,
            .phosphorus_g = 0.6,
        },
    });
    try std.testing.expectEqual(@as(f64, 131), result.seasonal_storage.carbon_g);
    try std.testing.expectEqual(@as(f64, 25), result.seasonal_storage.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 6), result.seasonal_storage.phosphorus_g);
    try std.testing.expectEqual(@as(f64, 31), result.recovered.carbon_g);
    try std.testing.expectEqual(@as(f64, 5), result.recovered.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 1), result.recovered.phosphorus_g);

    try std.testing.expectError(
        error.InvalidDeadBranchStorageRecoveryInput,
        sourceOrderDeadBranchStorageRecovery(.{
            .current_seasonal_storage = .{},
            .branch_mobile = .{},
            .c4_intermediate_carbon_g_c = -1,
            .stalk_reserve = .{},
        }),
    );
}

test "source-order dead branch canopy reset preserves node zero exceptions" {
    var scalar: SourceOrderDeadBranchScalarResetState = undefined;
    inline for (@typeInfo(SourceOrderDeadBranchScalarResetState).@"struct".fields) |field| {
        if (field.type == f64) {
            @field(scalar, field.name) = 1;
        } else {
            @field(scalar, field.name) = .{
                .carbon_g = 1,
                .nitrogen_g = 1,
                .phosphorus_g = 1,
            };
        }
    }
    const node: SourceOrderDeadBranchNodeResetState = .{
        .bundle_sheath_mobile_carbon_g_c = 1,
        .mesophyll_mobile_carbon_g_c = 1,
        .bundle_sheath_co2_carbon_g_c = 1,
        .bundle_sheath_bicarbonate_carbon_g_c = 1,
        .leaf_area_m2 = 1,
        .node_height_m = 1,
        .node_height_previous_m = 1,
        .sheath_height_m = 1,
        .leaf = .{ .carbon_g = 1, .nitrogen_g = 1, .phosphorus_g = 1 },
        .leaf_protein_g = 1,
        .sheath = .{ .carbon_g = 1, .nitrogen_g = 1, .phosphorus_g = 1 },
        .sheath_protein_g = 1,
        .stalk = .{ .carbon_g = 1, .nitrogen_g = 1, .phosphorus_g = 1 },
    };
    var nodes = [_]SourceOrderDeadBranchNodeResetState{ node, node };
    var node_layers = [_]SourceOrderDeadBranchLayerResetState{
        .{
            .leaf_area_m2 = 2,
            .leaf = .{ .carbon_g = 3, .nitrogen_g = 0.3, .phosphorus_g = 0.03 },
            .projected_leaf_surface_m2 = @splat(4),
        },
        .{
            .leaf_area_m2 = 5,
            .leaf = .{ .carbon_g = 7, .nitrogen_g = 0.7, .phosphorus_g = 0.07 },
            .projected_leaf_surface_m2 = @splat(8),
        },
    };
    var canopy_area = [_]f64{10};
    var canopy_carbon = [_]f64{20};
    var stalk_area = [_]f64{6};
    var stalk_surface = [_][4]f64{@splat(9)};
    try sourceOrderResetDeadBranchCanopy(.{
        .scalar = &scalar,
        .nodes = &nodes,
        .node_layers = &node_layers,
        .canopy_leaf_area_m2_by_layer = &canopy_area,
        .canopy_leaf_carbon_g_c_by_layer = &canopy_carbon,
        .branch_stalk_area_m2_by_layer = &stalk_area,
        .branch_projected_stalk_surface_m2 = &stalk_surface,
    });
    try std.testing.expectEqual(@as(f64, 3), canopy_area[0]);
    try std.testing.expectEqual(@as(f64, 10), canopy_carbon[0]);
    try std.testing.expectEqual(@as(f64, 0), scalar.host_mobile.carbon_g);
    try std.testing.expectEqual(@as(f64, 1), nodes[0].bundle_sheath_mobile_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 0), nodes[1].bundle_sheath_mobile_carbon_g_c);
    try std.testing.expectEqual([_]f64{4} ** 4, node_layers[0].projected_leaf_surface_m2);
    try std.testing.expectEqual([_]f64{0} ** 4, node_layers[1].projected_leaf_surface_m2);
    try std.testing.expectEqual(@as(f64, 0), stalk_area[0]);
    try std.testing.expectEqual([_]f64{0} ** 4, stalk_surface[0]);
}

test "source-order whole plant termination preserves winter annual reseed branch" {
    const current: SourceOrderWholePlantTerminationState = .{
        .shoot_alive = true,
        .root_alive = true,
        .total_node_count = 12,
        .hours_below_leaf_turgor_threshold_h = 8,
        .main_stalk_diameter_m = 0.03,
        .branch_count = 3,
        .living_population_per_m2 = 20,
        .living_population_count = 200,
        .hypocotyl_height_m = 0.1,
    };
    const partial = try sourceOrderWholePlantTermination(&.{ true, false, true }, false, current);
    try std.testing.expect(!partial.all_branches_dead);
    try std.testing.expectEqual(@as(usize, 2), partial.dead_branch_count);
    try std.testing.expectEqualDeep(current, partial.state);

    const winter = try sourceOrderWholePlantTermination(&.{ true, true, true }, true, current);
    try std.testing.expect(winter.all_branches_dead);
    try std.testing.expectEqual(@as(usize, 1), winter.state.branch_count);
    try std.testing.expectEqual(@as(f64, 20), winter.state.living_population_per_m2);
    try std.testing.expectEqual(@as(f64, 200), winter.state.living_population_count);
    try std.testing.expect(!winter.state.shoot_alive);
    try std.testing.expect(!winter.state.root_alive);
    try std.testing.expectEqual(@as(f64, 0), winter.state.hypocotyl_height_m);

    const ordinary = try sourceOrderWholePlantTermination(&.{ true, true, true }, false, current);
    try std.testing.expectEqual(@as(usize, 0), ordinary.state.branch_count);
    try std.testing.expectEqual(@as(f64, 0), ordinary.state.living_population_per_m2);
    try std.testing.expectEqual(@as(f64, 0), ordinary.state.living_population_count);
}

test "source-order dead root litterfall conserves runtime domains and axes" {
    const kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const mobile = [_]canopy.ElementalMass{
        .{ .carbon_g = 1, .nitrogen_g = 1, .phosphorus_g = 1 },
        .{ .carbon_g = 2, .nitrogen_g = 2, .phosphorus_g = 2 },
    };
    const axis: SourceOrderDeadRootAxisPools = .{
        .primary = .{ .carbon_g = 1, .nitrogen_g = 1, .phosphorus_g = 1 },
        .secondary = .{ .carbon_g = 2, .nitrogen_g = 2, .phosphorus_g = 2 },
    };
    const structural = [_]SourceOrderDeadRootAxisPools{ axis, axis, axis, axis };
    const base: SourceOrderDeadRootLitterInput = .{
        .roots_dead = true,
        .root_domain_count = 2,
        .soil_layer_count = 1,
        .root_axis_count = 2,
        .mobile_by_domain_layer = &mobile,
        .structural_by_domain_layer_axis = &structural,
        .root_woody_fraction = .{
            .carbon = .{ 0.25, 0.75 },
            .nitrogen = .{ 0.25, 0.75 },
            .phosphorus = .{ 0.25, 0.75 },
        },
        .mobile_kinetics = kinetics,
        .fine_root_kinetics = kinetics,
        .coarse_root_kinetics = kinetics,
    };
    const dead = try sourceOrderDeadRootLitterfall(std.testing.allocator, base);
    defer std.testing.allocator.free(dead);
    try std.testing.expectEqual(@as(f64, 3), dead[0].woody[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 12), dead[0].nonwoody[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 15), dead[0].woody[0].nitrogen_g +
        dead[0].nonwoody[0].nitrogen_g);

    var live_input = base;
    live_input.roots_dead = false;
    const live = try sourceOrderDeadRootLitterfall(std.testing.allocator, live_input);
    defer std.testing.allocator.free(live);
    try std.testing.expectEqual(@as(f64, 0), live[0].woody[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 0), live[0].nonwoody[0].carbon_g);
}

test "source-order dead root gas release clears both phases conservatively" {
    const gaseous: SourceOrderRootGasInventory = .{
        .carbon_dioxide_carbon_g_c = 1,
        .oxygen_g_o = 2,
        .methane_carbon_g_c = 3,
        .nitrous_oxide_nitrogen_g_n = 4,
        .ammonia_nitrogen_g_n = 5,
        .hydrogen_g_h = 6,
    };
    const aqueous: SourceOrderRootGasInventory = .{
        .carbon_dioxide_carbon_g_c = 10,
        .oxygen_g_o = 20,
        .methane_carbon_g_c = 30,
        .nitrous_oxide_nitrogen_g_n = 40,
        .ammonia_nitrogen_g_n = 50,
        .hydrogen_g_h = 60,
    };
    var phases = [_]SourceOrderRootGasPhases{
        .{ .gaseous = gaseous, .aqueous = aqueous },
        .{ .gaseous = gaseous, .aqueous = aqueous },
    };
    const loss = try sourceOrderReleaseDeadRootGases(
        true,
        &phases,
        std.mem.zeroes(SourceOrderRootGasInventory),
    );
    try std.testing.expectEqual(@as(f64, -22), loss.carbon_dioxide_carbon_g_c);
    try std.testing.expectEqual(@as(f64, -44), loss.oxygen_g_o);
    try std.testing.expectEqual(@as(f64, -132), loss.hydrogen_g_h);
    try std.testing.expectEqualDeep(
        std.mem.zeroes(SourceOrderRootGasPhases),
        phases[0],
    );
    try std.testing.expectEqualDeep(
        std.mem.zeroes(SourceOrderRootGasPhases),
        phases[1],
    );

    var living = [_]SourceOrderRootGasPhases{.{ .gaseous = gaseous, .aqueous = aqueous }};
    const unchanged_loss = try sourceOrderReleaseDeadRootGases(
        false,
        &living,
        gaseous,
    );
    try std.testing.expectEqualDeep(gaseous, unchanged_loss);
    try std.testing.expectEqualDeep(gaseous, living[0].gaseous);
    try std.testing.expectEqualDeep(aqueous, living[0].aqueous);
}

test "source-order dead root reset preserves runtime extents and radius defaults" {
    const mass: canopy.ElementalMass = .{
        .carbon_g = 1,
        .nitrogen_g = 1,
        .phosphorus_g = 1,
    };
    const axis_layer_value: SourceOrderDeadRootAxisLayerState = .{
        .primary = mass,
        .secondary = mass,
        .primary_length_m = 1,
        .secondary_length_m = 1,
        .secondary_axis_count = 1,
    };
    var axis_layers = [_]SourceOrderDeadRootAxisLayerState{
        axis_layer_value,
        axis_layer_value,
        axis_layer_value,
        axis_layer_value,
    };
    var domain_axes = [_]SourceOrderDeadRootDomainAxisState{
        .{ .primary_total = mass },
        .{ .primary_total = mass },
    };
    const layer_value: SourceOrderDeadRootDomainLayerState = .{
        .mobile = mass,
        .active_root_carbon_g_c = 1,
        .actual_root_carbon_g_c = 1,
        .root_protein_g = 1,
        .primary_axis_count = 1,
        .total_axis_count = 1,
        .root_length_per_plant_m = 1,
        .root_length_density_m_m3 = 1,
        .gaseous_volume_m3 = 1,
        .aqueous_volume_m3 = 1,
        .root_surface_area_per_plant_m2 = 1,
        .primary_radius_m = 1,
        .secondary_radius_m = 1,
        .average_secondary_root_length_m = 1,
    };
    var domain_layers = [_]SourceOrderDeadRootDomainLayerState{
        layer_value,
        layer_value,
    };
    try sourceOrderResetDeadRootState(true, .{
        .root_domain_count = 1,
        .soil_layer_count = 2,
        .root_axis_count = 2,
        .axis_layer = &axis_layers,
        .domain_axis = &domain_axes,
        .domain_layer = &domain_layers,
        .initial_primary_radius_m_by_domain = &.{0.01},
        .initial_secondary_radius_m_by_domain = &.{0.002},
        .initial_average_secondary_root_length_m = 0.1,
    });
    for (axis_layers) |axis| {
        try std.testing.expectEqual(@as(f64, 0), axis.primary.carbon_g);
        try std.testing.expectEqual(@as(f64, 0), axis.secondary_length_m);
        try std.testing.expectEqual(@as(f64, 0), axis.secondary_axis_count);
    }
    for (domain_axes) |axis|
        try std.testing.expectEqual(@as(f64, 0), axis.primary_total.carbon_g);
    for (domain_layers) |layer| {
        try std.testing.expectEqual(@as(f64, 0), layer.mobile.carbon_g);
        try std.testing.expectEqual(@as(f64, 0.01), layer.primary_radius_m);
        try std.testing.expectEqual(@as(f64, 0.002), layer.secondary_radius_m);
        try std.testing.expectEqual(@as(f64, 0.1), layer.average_secondary_root_length_m);
        try std.testing.expectEqual(@as(f64, 0), layer.root_surface_area_per_plant_m2);
    }
}

test "source-order dead nodule litterfall uses only first root domain" {
    const kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const pools: SourceOrderDeadNoduleLayerPools = .{
        .structural = .{ .carbon_g = 10, .nitrogen_g = 5, .phosphorus_g = 1 },
        .mobile = .{ .carbon_g = 2, .nitrogen_g = 1, .phosphorus_g = 0.2 },
    };
    var layers = [_]SourceOrderDeadNoduleLayerPools{ pools, pools };
    const litter = try sourceOrderDeadNoduleLitterfall(std.testing.allocator, .{
        .roots_dead = true,
        .nitrogen_fixation_enabled = true,
        .root_domain_count = 3,
        .layer_pools = &layers,
        .structural_kinetics = kinetics,
        .mobile_kinetics = kinetics,
    });
    defer std.testing.allocator.free(litter);
    try std.testing.expectEqual(@as(f64, 12), litter[0][0].carbon_g);
    try std.testing.expectEqual(@as(f64, 6), litter[1][0].nitrogen_g);
    try std.testing.expectEqual(@as(f64, 0), layers[0].structural.carbon_g);
    try std.testing.expectEqual(@as(f64, 0), layers[1].mobile.carbon_g);

    var nonfixing = [_]SourceOrderDeadNoduleLayerPools{pools};
    const empty = try sourceOrderDeadNoduleLitterfall(std.testing.allocator, .{
        .roots_dead = true,
        .nitrogen_fixation_enabled = false,
        .root_domain_count = 2,
        .layer_pools = &nonfixing,
        .structural_kinetics = kinetics,
        .mobile_kinetics = kinetics,
    });
    defer std.testing.allocator.free(empty);
    try std.testing.expectEqual(@as(f64, 0), empty[0][0].carbon_g);
    try std.testing.expectEqualDeep(pools, nonfixing[0]);
}

test "source-order dead root depth reset preserves axis-major domain order" {
    var deepest_by_axis = [_]usize{ 8, 9 };
    var depths = [_]f64{ 1, 2, 3, 4, 5, 6 };
    const mass: canopy.ElementalMass = .{
        .carbon_g = 1,
        .nitrogen_g = 0.1,
        .phosphorus_g = 0.01,
    };
    var totals = [_]canopy.ElementalMass{mass} ** 6;
    var deepest_active: usize = 9;
    var active_axes: usize = 2;
    try sourceOrderResetDeadRootDepth(true, 3, 0.04, .{
        .root_domain_count = 3,
        .root_axis_count = 2,
        .deepest_layer_by_axis = &deepest_by_axis,
        .primary_depth_from_surface_m_by_axis_domain = &depths,
        .primary_total_by_axis_domain = &totals,
        .deepest_active_root_layer = &deepest_active,
        .active_root_axis_count = &active_axes,
    });
    try std.testing.expectEqual([_]usize{ 3, 3 }, deepest_by_axis);
    try std.testing.expectEqual([_]f64{0.04} ** 6, depths);
    for (totals) |total|
        try std.testing.expectEqual(@as(f64, 0), total.carbon_g);
    try std.testing.expectEqual(@as(usize, 3), deepest_active);
    try std.testing.expectEqual(@as(usize, 0), active_axes);
}

test "source-order complete death shoot litterfall conserves storage and branches" {
    const kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const composition: TillageElementComposition = .{
        .carbon = .{ 0.25, 0.75 },
        .nitrogen = .{ 0.25, 0.75 },
        .phosphorus = .{ 0.25, 0.75 },
    };
    const carbon = struct {
        fn mass(value: f64) canopy.ElementalMass {
            return .{ .carbon_g = value, .nitrogen_g = 0, .phosphorus_g = 0 };
        }
    }.mass;
    const branch: SourceOrderCompleteDeathBranchPools = .{
        .host_mobile = carbon(1),
        .symbiont_mobile = carbon(2),
        .c4_intermediate_carbon_g_c = 3,
        .leaf = carbon(4),
        .symbiont_structural = carbon(5),
        .sheath = carbon(6),
        .husk = carbon(7),
        .ear = carbon(8),
        .grain = carbon(9),
        .stalk = carbon(10),
        .stalk_reserve = carbon(11),
    };
    const base: SourceOrderCompleteDeathShootInput = .{
        .shoot_dead = true,
        .roots_dead = true,
        .perennial_growth_habit = true,
        .deciduous_phenology = true,
        .seasonal_storage = carbon(8),
        .branches = &.{branch},
        .root_woody_fraction = composition,
        .leaf_woody_fraction = composition,
        .sheath_woody_fraction = composition,
        .nonstructural_kinetics = kinetics,
        .foliar_kinetics = kinetics,
        .nonfoliar_kinetics = kinetics,
        .stalk_kinetics = kinetics,
        .coarse_wood_kinetics = kinetics,
    };
    const perennial = try sourceOrderCompleteDeathShootLitterfall(base);
    try std.testing.expect(perennial.plant_death_initialized);
    try std.testing.expectEqual(@as(f64, 2), perennial.planting_layer_woody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 6), perennial.planting_layer_nonwoody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 42.5), perennial.surface_nonwoody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 2.5), perennial.surface_woody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 21), perennial.standing_dead_stalk[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 0), perennial.seasonal_storage.carbon_g);

    var winter_input = base;
    winter_input.perennial_growth_habit = false;
    const winter = try sourceOrderCompleteDeathShootLitterfall(winter_input);
    try std.testing.expect(!winter.plant_death_initialized);
    try std.testing.expectEqual(@as(f64, 33.5), winter.surface_nonwoody_litter[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 17), winter.seasonal_storage.carbon_g);
    try std.testing.expectEqual(@as(f64, 74), winter.surface_nonwoody_litter[0].carbon_g +
        winter.surface_woody_litter[0].carbon_g +
        winter.standing_dead_stalk[0].carbon_g +
        winter.seasonal_storage.carbon_g);
}

test "source-order complete death root litterfall conserves domains axes and layers" {
    const kinetics: litter_partition.ElementFractions = .{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const mobile = [_]canopy.ElementalMass{
        .{ .carbon_g = 1, .nitrogen_g = 2, .phosphorus_g = 3 },
        .{ .carbon_g = 10, .nitrogen_g = 20, .phosphorus_g = 30 },
        .{ .carbon_g = 4, .nitrogen_g = 5, .phosphorus_g = 6 },
        .{ .carbon_g = 40, .nitrogen_g = 50, .phosphorus_g = 60 },
    };
    const axis: SourceOrderDeadRootAxisPools = .{
        .primary = .{ .carbon_g = 2, .nitrogen_g = 2, .phosphorus_g = 2 },
        .secondary = .{ .carbon_g = 2, .nitrogen_g = 2, .phosphorus_g = 2 },
    };
    const structural = [_]SourceOrderDeadRootAxisPools{
        axis, axis, axis, axis, axis, axis, axis, axis,
    };
    const base: SourceOrderCompleteDeathRootInput = .{
        .shoot_dead = true,
        .roots_dead = true,
        .root_domain_count = 2,
        .soil_layer_count = 2,
        .root_axis_count = 2,
        .mobile_by_domain_layer = &mobile,
        .structural_by_domain_layer_axis = &structural,
        .root_woody_fraction = .{
            .carbon = .{ 0.25, 0.75 },
            .nitrogen = .{ 0.25, 0.75 },
            .phosphorus = .{ 0.25, 0.75 },
        },
        .nonstructural_kinetics = kinetics,
        .fine_root_kinetics = kinetics,
        .coarse_root_kinetics = kinetics,
    };
    const dead = try sourceOrderCompleteDeathRootLitterfall(
        std.testing.allocator,
        base,
    );
    defer std.testing.allocator.free(dead);
    try std.testing.expectEqual(@as(f64, 4), dead[0].woody[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 17), dead[0].nonwoody[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 62), dead[1].nonwoody[0].carbon_g);
    try std.testing.expectEqual(@as(f64, 21), dead[0].woody[0].carbon_g +
        dead[0].nonwoody[0].carbon_g);

    var partial_death = base;
    partial_death.shoot_dead = false;
    const live = try sourceOrderCompleteDeathRootLitterfall(
        std.testing.allocator,
        partial_death,
    );
    defer std.testing.allocator.free(live);
    try std.testing.expectEqual(@as(f64, 0), live[0].nonwoody[0].carbon_g);
}

test "source-order complete death branch reset clears every runtime branch field" {
    const mass: canopy.ElementalMass = .{
        .carbon_g = 1,
        .nitrogen_g = 2,
        .phosphorus_g = 3,
    };
    const populated: SourceOrderCompleteDeathBranchState = .{
        .host_mobile = mass,
        .c4_intermediate_carbon_g_c = 4,
        .symbiont_mobile = mass,
        .shoot = mass,
        .leaf = mass,
        .nodule = mass,
        .sheath = mass,
        .stalk = mass,
        .stalk_volume_m3 = 5,
        .reserve = mass,
        .husk = mass,
        .ear = mass,
        .grain = mass,
        .leaf_starch_carbon_g_c = 6,
        .stalk_extra = mass,
    };
    var partial = [_]SourceOrderCompleteDeathBranchState{populated};
    try sourceOrderResetCompleteDeathBranches(true, false, &partial);
    try std.testing.expectEqual(@as(f64, 1), partial[0].grain.carbon_g);

    var dead = [_]SourceOrderCompleteDeathBranchState{ populated, populated };
    try sourceOrderResetCompleteDeathBranches(true, true, &dead);
    for (dead) |branch| {
        inline for (@typeInfo(SourceOrderCompleteDeathBranchState).@"struct".fields) |field| {
            const value = @field(branch, field.name);
            if (field.type == f64) {
                try std.testing.expectEqual(@as(f64, 0), value);
            } else inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |element|
                try std.testing.expectEqual(@as(f64, 0), @field(value, element.name));
        }
    }
}

test "source-order complete death root reset preserves layer domain axis extents" {
    const mass: canopy.ElementalMass = .{
        .carbon_g = 1,
        .nitrogen_g = 2,
        .phosphorus_g = 3,
    };
    const axis: SourceOrderDeadRootAxisLayerState = .{
        .primary = mass,
        .secondary = mass,
        .primary_length_m = 4,
        .secondary_length_m = 5,
        .secondary_axis_count = 6,
    };
    var mobile = [_]canopy.ElementalMass{ mass, mass, mass, mass };
    var structural = [_]SourceOrderDeadRootAxisLayerState{
        axis, axis, axis, axis, axis, axis, axis, axis,
    };
    var totals = [_]SourceOrderDeadRootDomainAxisState{
        .{ .primary_total = mass },
        .{ .primary_total = mass },
        .{ .primary_total = mass },
        .{ .primary_total = mass },
    };
    const state: SourceOrderCompleteDeathRootResetState = .{
        .root_domain_count = 2,
        .soil_layer_count = 2,
        .root_axis_count = 2,
        .mobile_by_domain_layer = &mobile,
        .structural_by_domain_layer_axis = &structural,
        .primary_total_by_domain_axis = &totals,
    };
    try sourceOrderResetCompleteDeathRoots(true, false, state);
    try std.testing.expectEqual(@as(f64, 1), mobile[0].carbon_g);

    try sourceOrderResetCompleteDeathRoots(true, true, state);
    for (mobile) |value|
        try std.testing.expectEqual(@as(f64, 0), value.carbon_g);
    for (structural) |value| {
        try std.testing.expectEqual(@as(f64, 0), value.primary.carbon_g);
        try std.testing.expectEqual(@as(f64, 0), value.secondary.nitrogen_g);
        try std.testing.expectEqual(@as(f64, 0), value.primary_length_m);
        try std.testing.expectEqual(@as(f64, 0), value.secondary_length_m);
        try std.testing.expectEqual(@as(f64, 0), value.secondary_axis_count);
    }
    for (totals) |value|
        try std.testing.expectEqual(@as(f64, 0), value.primary_total.phosphorus_g);
}

test "source-order dead perennial reseed preserves next-day year rollover and gates" {
    const same_year = try sourceOrderScheduleDeadPerennialReseed(
        true,
        false,
        200,
        365,
        2001,
    );
    try std.testing.expect(same_year.plant_death_flag);
    try std.testing.expectEqual(@as(u16, 201), same_year.reseed_date.?.day_of_year);
    try std.testing.expectEqual(@as(u32, 2001), same_year.reseed_date.?.year);

    const rollover = try sourceOrderScheduleDeadPerennialReseed(
        true,
        false,
        366,
        366,
        2004,
    );
    try std.testing.expectEqual(@as(u16, 1), rollover.reseed_date.?.day_of_year);
    try std.testing.expectEqual(@as(u32, 2005), rollover.reseed_date.?.year);

    const annual = try sourceOrderScheduleDeadPerennialReseed(
        false,
        false,
        100,
        365,
        2001,
    );
    try std.testing.expect(annual.reseed_date == null);
    try std.testing.expect(!annual.plant_death_flag);
    const terminated = try sourceOrderScheduleDeadPerennialReseed(
        true,
        true,
        100,
        365,
        2001,
    );
    try std.testing.expect(terminated.reseed_date == null);
}

test "source-order soil plant exchange separates uptake fixation and NPP ledgers" {
    const result = try sourceOrderAccumulateSoilPlantExchange(.{
        .organic_carbon_exchange_g_c_step = 1,
        .organic_nitrogen_exchange_g_n_step = 2,
        .ammonium_uptake_g_n_step = 3,
        .nitrate_uptake_g_n_step = 4,
        .root_fixation_g_n_step = 5,
        .canopy_fixation_g_n_step = 6,
        .organic_phosphorus_exchange_g_p_step = 7,
        .dihydrogen_phosphate_uptake_g_p_step = 8,
        .hydrogen_phosphate_uptake_g_p_step = 9,
        .cumulative_soil_exchange = .{
            .carbon_g = 10,
            .nitrogen_g = 20,
            .phosphorus_g = 30,
        },
        .cumulative_fixation_g_n = 40,
        .cumulative_plant_carbon_g_c = 50,
        .cumulative_respired_carbon_g_c = -12,
    });
    try std.testing.expectEqual(@as(f64, 1), result.hourly_net_exchange.carbon_g);
    try std.testing.expectEqual(@as(f64, 14), result.hourly_net_exchange.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 24), result.hourly_net_exchange.phosphorus_g);
    try std.testing.expectEqual(@as(f64, 11), result.cumulative_soil_exchange.carbon_g);
    try std.testing.expectEqual(@as(f64, 29), result.cumulative_soil_exchange.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 54), result.cumulative_soil_exchange.phosphorus_g);
    try std.testing.expectEqual(@as(f64, 51), result.cumulative_fixation_g_n);
    try std.testing.expectEqual(
        @as(f64, 38),
        result.cumulative_net_primary_productivity_g_c,
    );
    try std.testing.expectEqual(
        @as(f64, 5),
        result.hourly_net_exchange.nitrogen_g -
            (result.cumulative_soil_exchange.nitrogen_g - 20),
    );
}

test "source-order standing dead geometry aggregates runtime components and layers" {
    const components = [_]canopy.ElementalMass{
        .{ .carbon_g = 1.1416, .nitrogen_g = 2, .phosphorus_g = 3 },
        .{ .carbon_g = 2, .nitrogen_g = 4, .phosphorus_g = 5 },
    };
    const edges = [_]f64{ 0, 0.5, 1 };
    const result = try sourceOrderStandingDeadGeometry(std.testing.allocator, .{
        .components = &components,
        .negligible_mass_g_c = 1.0e-12,
        .standing_dead_population_count = 1,
        .previous_height_m = 0.75,
        .canopy_height_m = 1,
        .stalk_volume_per_carbon_m3_g_c = 1,
        .canopy_layer_edges_m = &edges,
    });
    defer result.deinit(std.testing.allocator);
    try std.testing.expectApproxEqAbs(@as(f64, 3.1416), result.total_mass.carbon_g, 1.0e-12);
    try std.testing.expectEqual(@as(f64, 6), result.total_mass.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 8), result.total_mass.phosphorus_g);
    try std.testing.expectApproxEqAbs(@as(f64, 1), result.height_m, 1.0e-12);
    try std.testing.expectApproxEqAbs(
        @as(f64, 6.2832),
        result.total_surface_area_m2,
        1.0e-12,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 3.1416),
        result.layer_surface_area_m2[0],
        1.0e-12,
    );
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.7854),
        result.projected_surface_area_m2[1],
        1.0e-12,
    );

    const empty = try sourceOrderStandingDeadGeometry(std.testing.allocator, .{
        .components = &.{.{}},
        .negligible_mass_g_c = 1.0e-12,
        .standing_dead_population_count = 1,
        .previous_height_m = 1,
        .canopy_height_m = 1,
        .stalk_volume_per_carbon_m3_g_c = 1,
        .canopy_layer_edges_m = &edges,
    });
    defer empty.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(f64, 0), empty.height_m);
    try std.testing.expectEqual(@as(f64, 0), empty.layer_surface_area_m2[1]);
}

test "source-order fire inventory preserves plant layer domain and axis sums" {
    const carbon = struct {
        fn pools(value: f64) SourceOrderDeadRootAxisPools {
            return .{
                .primary = .{ .carbon_g = value },
                .secondary = .{ .carbon_g = value },
            };
        }
    }.pools;
    const structural_a = [_]SourceOrderDeadRootAxisPools{
        carbon(1), carbon(1), carbon(1), carbon(1),
        carbon(1), carbon(1), carbon(1), carbon(1),
    };
    const structural_b = [_]SourceOrderDeadRootAxisPools{ carbon(2), carbon(3) };
    const shoot_a: SourceOrderFireShootCarbon = .{
        .canopy_nonstructural_g_c = 1,
        .leaf_g_c = 2,
        .sheath_g_c = 3,
        .stalk_g_c = 4,
        .reserve_g_c = 5,
        .husk_g_c = 6,
        .ear_g_c = 7,
        .grain_g_c = 8,
        .symbiont_nonstructural_g_c = 9,
        .symbiont_biomass_g_c = 10,
        .seasonal_storage_g_c = 11,
        .standing_dead_g_c = 12,
    };
    var shoot_b = shoot_a;
    shoot_b.symbiont_nonstructural_g_c = 90;
    shoot_b.symbiont_biomass_g_c = 100;
    const plants = [_]SourceOrderFirePlantInventory{
        .{
            .shoot = shoot_a,
            .canopy_symbiont_included = true,
            .root_domain_count = 2,
            .root_axis_count = 2,
            .nodule_nonstructural_by_layer_g_c = &.{ 1, 2 },
            .nodule_biomass_by_layer_g_c = &.{ 3, 4 },
            .root_nonstructural_by_layer_domain_g_c = &.{ 1, 2, 3, 4 },
            .root_structural_by_layer_domain_axis = &structural_a,
        },
        .{
            .shoot = shoot_b,
            .canopy_symbiont_included = false,
            .root_domain_count = 1,
            .root_axis_count = 1,
            .nodule_nonstructural_by_layer_g_c = &.{ 10, 20 },
            .nodule_biomass_by_layer_g_c = &.{ 30, 40 },
            .root_nonstructural_by_layer_domain_g_c = &.{ 10, 20 },
            .root_structural_by_layer_domain_axis = &structural_b,
        },
    };
    const optional = try sourceOrderAggregateFireCarbonInventory(
        std.testing.allocator,
        true,
        2,
        &plants,
    );
    const result = optional.?;
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(f64, 2), result.shoot.canopy_nonstructural_g_c);
    try std.testing.expectEqual(@as(f64, 9), result.shoot.symbiont_nonstructural_g_c);
    try std.testing.expectEqual(@as(f64, 22), result.shoot.seasonal_storage_g_c);
    try std.testing.expectEqual(@as(f64, 13), result.layers[0].root_nonstructural_g_c);
    try std.testing.expectEqual(@as(f64, 12), result.layers[0].root_structural_g_c);
    try std.testing.expectEqual(@as(f64, 11), result.layers[0].nodule_nonstructural_g_c);
    try std.testing.expectEqual(@as(f64, 44), result.layers[1].nodule_biomass_g_c);
    try std.testing.expect((try sourceOrderAggregateFireCarbonInventory(
        std.testing.allocator,
        false,
        0,
        &.{},
    )) == null);
}

test "source-order combustion rates preserve independent temperature gates" {
    const specific: SourceOrderCombustionSpecificRates = .{
        .living_nonstructural_and_leaf_g_c_m2_h = 1,
        .living_sheath_g_c_m2_h = 2,
        .living_stalk_g_c_m2_h = 3,
        .living_reproductive_g_c_m2_h = 4,
        .standing_dead_g_c_m2_h = 5,
    };
    const living = try sourceOrderCombustionRates(
        600,
        500,
        550,
        0.5,
        2,
        3,
        specific,
    );
    try std.testing.expectEqual(@as(f64, 0.5), living.living_temperature_fraction);
    try std.testing.expectEqual(@as(f64, 3), living.living_nonstructural_and_leaf_g_c_step);
    try std.testing.expectEqual(@as(f64, 6), living.living_sheath_g_c_step);
    try std.testing.expectEqual(@as(f64, 9), living.living_stalk_g_c_step);
    try std.testing.expectEqual(@as(f64, 12), living.living_reproductive_g_c_step);
    try std.testing.expectEqual(@as(f64, 0), living.standing_dead_g_c_step);

    const dead = try sourceOrderCombustionRates(
        500,
        600,
        550,
        0.5,
        2,
        3,
        specific,
    );
    try std.testing.expectEqual(@as(f64, 0), dead.living_temperature_fraction);
    try std.testing.expectEqual(@as(f64, 0.5), dead.standing_dead_temperature_fraction);
    try std.testing.expectEqual(@as(f64, 15), dead.standing_dead_g_c_step);

    const cold = try sourceOrderCombustionRates(
        550,
        550,
        550,
        0.5,
        2,
        3,
        specific,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        cold.living_nonstructural_and_leaf_g_c_step,
    );
    try std.testing.expectEqual(@as(f64, 0), cold.standing_dead_g_c_step);
}

test "source-order shoot combustion fractions preserve pool gates and routing" {
    const totals: SourceOrderShootCombustionTotals = .{
        .canopy_nonstructural_g_c = 10,
        .leaf_g_c = 10,
        .sheath_g_c = 10,
        .stalk_g_c = 10,
        .husk_g_c = 10,
        .ear_g_c = 10,
        .grain_g_c = 10,
        .symbiont_nonstructural_g_c = 10,
        .symbiont_biomass_g_c = 10,
        .standing_dead_g_c = 10,
    };
    const rates: SourceOrderCombustionRates = .{
        .living_temperature_fraction = 1,
        .standing_dead_temperature_fraction = 1,
        .living_nonstructural_and_leaf_g_c_step = 1,
        .living_sheath_g_c_step = 2,
        .living_stalk_g_c_step = 3,
        .living_reproductive_g_c_step = 4,
        .standing_dead_g_c_step = 5,
    };
    const woody = try sourceOrderShootCombustionFractions(
        totals,
        rates,
        0,
        1,
        2,
    );
    try std.testing.expectEqual(@as(f64, 0.1), woody.canopy_nonstructural);
    try std.testing.expectEqual(@as(f64, 0.3), woody.sheath);
    try std.testing.expectEqual(@as(f64, 0.4), woody.stalk);
    try std.testing.expectEqual(woody.stalk, woody.reserve);
    try std.testing.expectEqual(@as(f64, 0.2), woody.ear);
    try std.testing.expectEqual(@as(f64, 0.5), woody.standing_dead);

    const herbaceous = try sourceOrderShootCombustionFractions(
        totals,
        rates,
        0,
        0,
        2,
    );
    try std.testing.expectEqual(@as(f64, 0.2), herbaceous.sheath);
    try std.testing.expectEqual(@as(f64, 0.2), herbaceous.stalk);
    var sparse = totals;
    sparse.leaf_g_c = 1.0e-12;
    const gated = try sourceOrderShootCombustionFractions(
        sparse,
        rates,
        1.0e-12,
        0,
        2,
    );
    try std.testing.expectEqual(@as(f64, 0), gated.leaf);
    inline for (@typeInfo(SourceOrderShootCombustionFractions).@"struct".fields) |field| {
        const value = @field(woody, field.name);
        try std.testing.expect(value >= 0 and value <= 1);
    }
}

test "source-order branch combustion conserves C N P and signed ledgers" {
    const mass: canopy.ElementalMass = .{
        .carbon_g = 1,
        .nitrogen_g = 2,
        .phosphorus_g = 3,
    };
    const branch: SourceOrderShootCombustionBranchPools = .{
        .canopy_nonstructural = mass,
        .leaf = mass,
        .sheath = mass,
        .stalk = mass,
        .reserve = mass,
        .husk = mass,
        .ear = mass,
        .grain = mass,
        .symbiont_nonstructural = mass,
        .symbiont_biomass = mass,
    };
    const fractions: SourceOrderShootCombustionFractions = .{
        .canopy_nonstructural = 0.5,
        .leaf = 0.5,
        .sheath = 0.5,
        .stalk = 0.5,
        .reserve = 0.5,
        .husk = 0.5,
        .ear = 0.5,
        .grain = 0.5,
        .symbiont_nonstructural = 0.5,
        .symbiont_biomass = 0.5,
        .standing_dead = 0.5,
    };
    const result = try sourceOrderShootCombustionLosses(
        std.testing.allocator,
        &.{ branch, branch },
        fractions,
        1,
        .{ .carbon_g = 100, .nitrogen_g = 200, .phosphorus_g = 300 },
    );
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 2), result.branches.len);
    try std.testing.expectEqual(@as(f64, 5), result.branches[0].total_combusted.carbon_g);
    try std.testing.expectEqual(@as(f64, 10), result.branches[0].total_combusted.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 15), result.branches[0].total_combusted.phosphorus_g);
    try std.testing.expectEqual(@as(f64, 11), result.cumulative_canopy_combustion_g_c);
    try std.testing.expectEqual(@as(f64, 90), result.disturbance_emission_ledger.carbon_g);
    try std.testing.expectEqual(@as(f64, 180), result.disturbance_emission_ledger.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 270), result.disturbance_emission_ledger.phosphorus_g);
    try std.testing.expectEqual(
        @as(f64, 10),
        100 - result.disturbance_emission_ledger.carbon_g,
    );
}

test "source-order shoot salt combustion conserves all eight species" {
    const branch: SourceOrderShootSaltInventory = .{
        .aluminum_mol = 1,
        .iron_mol = 2,
        .calcium_mol = 3,
        .magnesium_mol = 4,
        .sodium_mol = 5,
        .potassium_mol = 6,
        .sulfate_mol = 7,
        .chloride_mol = 8,
    };
    const optional = try sourceOrderShootSaltCombustion(
        std.testing.allocator,
        true,
        0.25,
        &.{ branch, branch },
    );
    const result = optional.?;
    defer std.testing.allocator.free(result);
    try std.testing.expectEqual(@as(usize, 2), result.len);
    try std.testing.expectEqual(@as(f64, 0.25), result[0].combusted.aluminum_mol);
    try std.testing.expectEqual(@as(f64, 2), result[0].combusted.chloride_mol);
    try std.testing.expectEqual(@as(f64, 6), result[0].remaining.chloride_mol);
    inline for (@typeInfo(SourceOrderShootSaltInventory).@"struct".fields) |field| {
        try std.testing.expectEqual(
            @field(branch, field.name),
            @field(result[0].combusted, field.name) +
                @field(result[0].remaining, field.name),
        );
    }
    try std.testing.expect((try sourceOrderShootSaltCombustion(
        std.testing.allocator,
        false,
        2,
        &.{},
    )) == null);
}

test "source-order uncombusted shoot state conserves pools and runtime topology" {
    const pool: canopy.ElementalMass = .{
        .carbon_g = 10,
        .nitrogen_g = 20,
        .phosphorus_g = 30,
    };
    const burned_mass: canopy.ElementalMass = .{
        .carbon_g = 2,
        .nitrogen_g = 4,
        .phosphorus_g = 6,
    };
    const pools: SourceOrderShootCombustionBranchPools = .{
        .canopy_nonstructural = pool,
        .leaf = pool,
        .sheath = pool,
        .stalk = pool,
        .reserve = pool,
        .husk = pool,
        .ear = pool,
        .grain = pool,
        .symbiont_nonstructural = pool,
        .symbiont_biomass = pool,
    };
    const burned: SourceOrderShootCombustionBranchPools = .{
        .canopy_nonstructural = burned_mass,
        .leaf = burned_mass,
        .sheath = burned_mass,
        .stalk = burned_mass,
        .reserve = burned_mass,
        .husk = burned_mass,
        .ear = burned_mass,
        .grain = burned_mass,
        .symbiont_nonstructural = burned_mass,
        .symbiont_biomass = burned_mass,
    };
    var nodes = [_]SourceOrderShootCombustionNodeState{.{
        .leaf_area_m2 = 4,
        .sheath_height_m = 6,
        .green_leaf = pool,
        .senescent_leaf_carbon_g_c = 8,
        .green_sheath = pool,
        .senescent_sheath_carbon_g_c = 10,
        .node = pool,
    }};
    var layers = [_]SourceOrderShootCombustionNodeLayerState{
        .{ .leaf_area_m2 = 2, .green_leaf = pool },
        .{ .leaf_area_m2 = 4, .green_leaf = pool },
    };
    var branches = [_]SourceOrderUncombustedBranchState{.{
        .pools = pools,
        .c4_intermediate_carbon_g_c = 1,
        .total_shoot = .{},
        .leaf_area_m2 = 12,
        .nodes = &nodes,
        .node_layers = &layers,
        .canopy_layer_count = 2,
    }};
    var fractions = std.mem.zeroes(SourceOrderShootCombustionFractions);
    fractions.leaf = 0.25;
    fractions.sheath = 0.5;
    fractions.stalk = 0.75;
    try sourceOrderApplyUncombustedShootState(&branches, &.{burned}, fractions);
    try std.testing.expectEqual(@as(f64, 8), branches[0].pools.leaf.carbon_g);
    try std.testing.expectEqual(@as(f64, 65), branches[0].total_shoot.carbon_g);
    try std.testing.expectEqual(@as(f64, 128), branches[0].total_shoot.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 9), branches[0].leaf_area_m2);
    try std.testing.expectEqual(@as(f64, 3), nodes[0].leaf_area_m2);
    try std.testing.expectEqual(@as(f64, 3), nodes[0].sheath_height_m);
    try std.testing.expectEqual(@as(f64, 2.5), nodes[0].node.carbon_g);
    try std.testing.expectEqual(@as(f64, 1.5), layers[0].leaf_area_m2);
    try std.testing.expectEqual(@as(f64, 7.5), layers[1].green_leaf.carbon_g);
}

test "source-order standing dead combustion conserves components and ledgers" {
    const components = [_]canopy.ElementalMass{
        .{ .carbon_g = 1, .nitrogen_g = 2, .phosphorus_g = 3 },
        .{ .carbon_g = 2, .nitrogen_g = 4, .phosphorus_g = 6 },
        .{ .carbon_g = 3, .nitrogen_g = 6, .phosphorus_g = 9 },
        .{ .carbon_g = 4, .nitrogen_g = 8, .phosphorus_g = 12 },
    };
    const result = try sourceOrderStandingDeadCombustion(
        std.testing.allocator,
        &components,
        0.25,
        10,
        .{ .carbon_g = 100, .nitrogen_g = 200, .phosphorus_g = 300 },
    );
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 4), result.components.len);
    try std.testing.expectEqual(@as(f64, 0.25), result.components[0].combusted.carbon_g);
    try std.testing.expectEqual(@as(f64, 3), result.components[3].remaining.carbon_g);
    try std.testing.expectEqual(
        @as(f64, 12.5),
        result.cumulative_standing_dead_combustion_g_c,
    );
    try std.testing.expectEqual(@as(f64, 97.5), result.disturbance_emission_ledger.carbon_g);
    try std.testing.expectEqual(@as(f64, 195), result.disturbance_emission_ledger.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 292.5), result.disturbance_emission_ledger.phosphorus_g);
    for (components, result.components) |before, after| {
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
            try std.testing.expectEqual(
                @field(before, field.name),
                @field(after.combusted, field.name) +
                    @field(after.remaining, field.name),
            );
    }
}

test "source-order charcoal combustion preserves response fraction and ledgers" {
    const result = try sourceOrderCharcoalCombustion(
        700,
        0.5,
        2,
        3,
        4,
        24,
        1.0e-12,
        .{ .carbon_g = 10, .nitrogen_g = 2, .phosphorus_g = 1 },
        7,
        .{ .carbon_g = 100, .nitrogen_g = 20, .phosphorus_g = 10 },
        100,
        3,
    );
    try std.testing.expectEqual(@as(f64, 0.5), result.temperature_response);
    try std.testing.expectEqual(@as(f64, 12), result.potential_combustion_g_c_step);
    try std.testing.expectEqual(@as(f64, 0.5), result.combustion_fraction);
    try std.testing.expectEqual(@as(f64, 5), result.combusted.carbon_g);
    try std.testing.expectEqual(@as(f64, 5), result.remaining.carbon_g);
    try std.testing.expectEqual(@as(f64, 12), result.cumulative_standing_dead_combustion_g_c);
    try std.testing.expectEqual(@as(f64, 95), result.disturbance_emission_ledger.carbon_g);
    try std.testing.expectEqual(@as(f64, 19), result.disturbance_emission_ledger.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 9.5), result.disturbance_emission_ledger.phosphorus_g);
    try std.testing.expectEqual(@as(f64, 115), result.grid_total_combustion_g_c);
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        try std.testing.expectEqual(
            @field(result.combusted, field.name) +
                @field(result.remaining, field.name),
            @field(canopy.ElementalMass{
                .carbon_g = 10,
                .nitrogen_g = 2,
                .phosphorus_g = 1,
            }, field.name),
        );
}

test "source-order no-combustion branch resets every combustion rate" {
    const result = try sourceOrderResetNoCombustion(std.testing.allocator, 3, true);
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(f64, 0), result.canopy_combustion_g_c_step);
    try std.testing.expectEqual(@as(f64, 0), result.standing_dead_combustion_g_c_step);
    try std.testing.expectEqual(@as(f64, 0), result.canopy_temperature_response);
    try std.testing.expectEqual(@as(f64, 0), result.standing_dead_temperature_response);
    try std.testing.expectEqual(@as(usize, 3), result.branch_combustion.len);
    try std.testing.expectEqual(@as(usize, 5), result.standing_dead_combustion.len);
    for (result.branch_combustion) |branch|
        inline for (@typeInfo(SourceOrderShootCombustionBranchPools).@"struct".fields) |field|
            try std.testing.expectEqual(
                std.mem.zeroes(canopy.ElementalMass),
                @field(branch, field.name),
            );
    const salts = result.branch_salt_combustion.?;
    try std.testing.expectEqual(@as(usize, 3), salts.len);
    for (salts) |salt|
        inline for (@typeInfo(SourceOrderShootSaltInventory).@"struct".fields) |field|
            try std.testing.expectEqual(@as(f64, 0), @field(salt, field.name));
    for (result.standing_dead_combustion) |component|
        try std.testing.expectEqual(std.mem.zeroes(canopy.ElementalMass), component);

    const without_salt = try sourceOrderResetNoCombustion(
        std.testing.allocator,
        0,
        false,
    );
    defer without_salt.deinit(std.testing.allocator);
    try std.testing.expect(without_salt.branch_salt_combustion == null);
}

test "source-order root and surface storage combustion preserves source fractions" {
    const specific: SourceOrderCombustionSpecificRates = .{
        .living_nonstructural_and_leaf_g_c_m2_h = 1,
        .living_sheath_g_c_m2_h = 2,
        .living_stalk_g_c_m2_h = 3,
        .living_reproductive_g_c_m2_h = 4,
        .standing_dead_g_c_m2_h = 5,
    };
    const result = try sourceOrderRootStorageCombustion(.{
        .soil_temperature_k = 600,
        .minimum_combustion_temperature_k = 500,
        .maximum_temperature_response = 0.5,
        .surface_area_m2 = 2,
        .biological_timestep_h = 1,
        .specific_rates = specific,
        .totals = .{
            .root_nonstructural_g_c = 4,
            .active_root_g_c = 12,
            .nodule_nonstructural_g_c = 8,
            .nodule_biomass_g_c = 16,
        },
        .negligible_carbon_g_c = 1.0e-12,
        .is_surface_layer = true,
        .storage = .{ .carbon_g = 6, .nitrogen_g = 3, .phosphorus_g = 1.5 },
        .preceding_layer_combustion_g_c = 10,
        .preceding_disturbance_emission_ledger = .{ .carbon_g = 100, .nitrogen_g = 50, .phosphorus_g = 25 },
    });
    try std.testing.expectEqual(@as(f64, 0.5), result.temperature_response);
    try std.testing.expectEqual(@as(f64, 1), result.potential_rates.root_nonstructural_g_c_step);
    try std.testing.expectEqual(@as(f64, 2), result.potential_rates.nodule_biomass_g_c_step);
    try std.testing.expectEqual(@as(f64, 3), result.potential_rates.active_root_g_c_step);
    try std.testing.expectEqual(@as(f64, 0.25), result.fractions.root_nonstructural);
    try std.testing.expectEqual(@as(f64, 0.25), result.fractions.active_root);
    try std.testing.expectEqual(@as(f64, 0.125), result.fractions.nodule_nonstructural);
    try std.testing.expectEqual(@as(f64, 0.125), result.fractions.nodule_biomass);
    try std.testing.expectEqual(result.fractions.active_root, result.fractions.storage);
    try std.testing.expectEqual(@as(f64, 1.5), result.storage_combusted.carbon_g);
    try std.testing.expectEqual(@as(f64, 4.5), result.storage_remaining.carbon_g);
    try std.testing.expectEqual(@as(f64, 11.5), result.layer_combustion_g_c);
    try std.testing.expectEqual(@as(f64, 98.5), result.disturbance_emission_ledger.carbon_g);
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        try std.testing.expectEqual(
            @field(result.storage_combusted, field.name) +
                @field(result.storage_remaining, field.name),
            @field(canopy.ElementalMass{
                .carbon_g = 6,
                .nitrogen_g = 3,
                .phosphorus_g = 1.5,
            }, field.name),
        );

    const subsurface = try sourceOrderRootStorageCombustion(.{
        .soil_temperature_k = 600,
        .minimum_combustion_temperature_k = 500,
        .maximum_temperature_response = 0.5,
        .surface_area_m2 = 2,
        .biological_timestep_h = 1,
        .specific_rates = specific,
        .totals = .{
            .root_nonstructural_g_c = 4,
            .active_root_g_c = 12,
            .nodule_nonstructural_g_c = 8,
            .nodule_biomass_g_c = 16,
        },
        .negligible_carbon_g_c = 1.0e-12,
        .is_surface_layer = false,
        .storage = .{ .carbon_g = 6, .nitrogen_g = 3, .phosphorus_g = 1.5 },
        .preceding_layer_combustion_g_c = 10,
        .preceding_disturbance_emission_ledger = .{ .carbon_g = 100, .nitrogen_g = 50, .phosphorus_g = 25 },
    });
    try std.testing.expectEqual(@as(f64, 0), subsurface.fractions.storage);
    try std.testing.expectEqual(@as(f64, 6), subsurface.storage_remaining.carbon_g);
    try std.testing.expectEqual(@as(f64, 10), subsurface.layer_combustion_g_c);
}

test "source-order root-domain combustion conserves runtime domains and axes" {
    var axes = [_]SourceOrderRootCombustionAxisState{
        .{
            .primary = .{ .carbon_g = 8, .nitrogen_g = 4, .phosphorus_g = 2 },
            .secondary = .{ .carbon_g = 4, .nitrogen_g = 2, .phosphorus_g = 1 },
            .whole_primary = .{ .carbon_g = 16, .nitrogen_g = 8, .phosphorus_g = 4 },
            .primary_length_m = 10,
            .secondary_length_m = 6,
            .secondary_root_number = 2,
        },
        .{
            .primary = .{ .carbon_g = 12, .nitrogen_g = 6, .phosphorus_g = 3 },
            .secondary = .{ .carbon_g = 8, .nitrogen_g = 4, .phosphorus_g = 2 },
            .whole_primary = .{ .carbon_g = 20, .nitrogen_g = 10, .phosphorus_g = 5 },
            .primary_length_m = 14,
            .secondary_length_m = 8,
            .secondary_root_number = 4,
        },
    };
    var domains = [_]SourceOrderRootCombustionDomainState{.{
        .nonstructural = .{ .carbon_g = 10, .nitrogen_g = 5, .phosphorus_g = 2 },
        .salts = .{
            .aluminum_mol = 1,
            .iron_mol = 2,
            .calcium_mol = 3,
            .magnesium_mol = 4,
            .sodium_mol = 5,
            .potassium_mol = 6,
            .sulfate_mol = 7,
            .chloride_mol = 8,
        },
        .active_root_carbon_g_c = 100,
        .root_density_g_c_m3 = 20,
        .root_surface_area_m2 = 30,
        .primary_root_number = 4,
        .root_length_m = 5,
        .root_length_growth_m_step = 6,
        .root_depth_growth_m_step = 7,
        .root_volume_growth_m3_step = 8,
        .root_volume_m3 = 9,
        .root_area_m2 = 10,
        .axes = &axes,
    }};
    const result = try sourceOrderApplyRootDomainCombustion(
        std.testing.allocator,
        &domains,
        0.2,
        0.25,
        true,
        5,
        .{ .carbon_g = 100, .nitrogen_g = 50, .phosphorus_g = 25 },
    );
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 1), result.losses.len);
    try std.testing.expectEqual(@as(f64, 2), result.losses[0].nonstructural.carbon_g);
    try std.testing.expectEqual(@as(f64, 8), domains[0].nonstructural.carbon_g);
    try std.testing.expectEqual(@as(f64, 0.2), result.losses[0].salts.?.aluminum_mol);
    try std.testing.expectEqual(@as(f64, 0.8), domains[0].salts.aluminum_mol);
    try std.testing.expectEqual(@as(f64, 8), result.losses[0].structural.carbon_g);
    try std.testing.expectEqual(@as(f64, 6), axes[0].primary.carbon_g);
    try std.testing.expectEqual(@as(f64, 3), axes[0].secondary.carbon_g);
    try std.testing.expectEqual(@as(f64, 75), domains[0].active_root_carbon_g_c);
    try std.testing.expectEqual(@as(f64, 7.5), axes[0].primary_length_m);
    try std.testing.expectEqual(@as(f64, 15), result.layer_combustion_g_c);
    try std.testing.expectEqual(@as(f64, 90), result.disturbance_emission_ledger.carbon_g);
    try std.testing.expectEqual(@as(f64, 45), result.disturbance_emission_ledger.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 22.6), result.disturbance_emission_ledger.phosphorus_g);
}

test "source-order root-nodule combustion conserves pools and source ledgers" {
    const result = try sourceOrderRootNoduleCombustion(
        .{ .carbon_g = 8, .nitrogen_g = 4, .phosphorus_g = 2 },
        .{ .carbon_g = 12, .nitrogen_g = 6, .phosphorus_g = 3 },
        0.25,
        0.5,
        10,
        100,
        .{ .carbon_g = 50, .nitrogen_g = 25, .phosphorus_g = 12 },
    );
    try std.testing.expectEqual(@as(f64, 2), result.nonstructural_combusted.carbon_g);
    try std.testing.expectEqual(@as(f64, 6), result.nonstructural_remaining.carbon_g);
    try std.testing.expectEqual(@as(f64, 6), result.biomass_combusted.carbon_g);
    try std.testing.expectEqual(@as(f64, 6), result.biomass_remaining.carbon_g);
    try std.testing.expectEqual(@as(f64, 18), result.layer_plant_combustion_g_c);
    try std.testing.expectEqual(@as(f64, 118), result.grid_layer_combustion_g_c);
    try std.testing.expectEqual(@as(f64, 42), result.disturbance_emission_ledger.carbon_g);
    try std.testing.expectEqual(@as(f64, 21), result.disturbance_emission_ledger.nitrogen_g);
    try std.testing.expectEqual(@as(f64, 10), result.disturbance_emission_ledger.phosphorus_g);
    inline for (.{ .{
        result.nonstructural_combusted,
        result.nonstructural_remaining,
        canopy.ElementalMass{ .carbon_g = 8, .nitrogen_g = 4, .phosphorus_g = 2 },
    }, .{
        result.biomass_combusted,
        result.biomass_remaining,
        canopy.ElementalMass{ .carbon_g = 12, .nitrogen_g = 6, .phosphorus_g = 3 },
    } }) |pools|
        inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
            try std.testing.expectEqual(
                @field(pools[2], field.name),
                @field(pools[0], field.name) + @field(pools[1], field.name),
            );
}

test "source-order cold-soil reset preserves heterogeneous runtime topology" {
    const axis_counts = [_]usize{ 2, 0, 1 };
    const result = try sourceOrderResetColdSoilCombustion(
        std.testing.allocator,
        &axis_counts,
        true,
        true,
    );
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(f64, 0), result.layer_plant_combustion_g_c);
    try std.testing.expectEqual(std.mem.zeroes(canopy.ElementalMass), result.storage_combustion.?);
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 2, 3 }, result.axis_offsets);
    try std.testing.expectEqual(@as(usize, 3), result.domains.len);
    try std.testing.expectEqual(@as(usize, 3), result.axes.len);
    for (result.domains) |domain| {
        try std.testing.expectEqual(std.mem.zeroes(canopy.ElementalMass), domain.nonstructural);
        try std.testing.expectEqual(std.mem.zeroes(canopy.ElementalMass), domain.structural);
        const salts = domain.salts.?;
        inline for (@typeInfo(SourceOrderShootSaltInventory).@"struct".fields) |field|
            try std.testing.expectEqual(@as(f64, 0), @field(salts, field.name));
    }
    for (result.axes) |axis| {
        try std.testing.expectEqual(std.mem.zeroes(canopy.ElementalMass), axis.primary);
        try std.testing.expectEqual(std.mem.zeroes(canopy.ElementalMass), axis.secondary);
    }
    try std.testing.expectEqual(
        std.mem.zeroes(canopy.ElementalMass),
        result.nodule_nonstructural_combustion,
    );
    try std.testing.expectEqual(
        std.mem.zeroes(canopy.ElementalMass),
        result.nodule_biomass_combustion,
    );

    const subsurface = try sourceOrderResetColdSoilCombustion(
        std.testing.allocator,
        &.{},
        false,
        false,
    );
    defer subsurface.deinit(std.testing.allocator);
    try std.testing.expect(subsurface.storage_combustion == null);
    try std.testing.expectEqualSlices(usize, &.{0}, subsurface.axis_offsets);
}

test "source-order dormant-seed activation uses first qualifying branch" {
    const activation = (try sourceOrderDormantSeedActivation(
        true,
        true,
        &.{
            .{
                .leafout_disabled = true,
                .accumulated_leafout_h = 100,
                .required_leafout_h = 10,
            },
            .{
                .leafout_disabled = false,
                .accumulated_leafout_h = 9,
                .required_leafout_h = 10,
            },
            .{
                .leafout_disabled = false,
                .accumulated_leafout_h = 10,
                .required_leafout_h = 10,
            },
            .{
                .leafout_disabled = false,
                .accumulated_leafout_h = 20,
                .required_leafout_h = 10,
            },
        },
        150,
        2025,
        0.02,
    )).?;
    try std.testing.expectEqual(@as(usize, 2), activation.qualifying_branch_index);
    try std.testing.expectEqual(@as(u16, 150), activation.planting_day_of_year);
    try std.testing.expectEqual(@as(i32, 2025), activation.planting_year);
    try std.testing.expectEqual(@as(f64, 0.025), activation.seeding_depth_m);
    try std.testing.expect(!activation.initialization_pending);
}

test "source-order dormant-seed activation preserves iteration and pending gates" {
    const ready = [_]SourceOrderDormantSeedBranch{.{
        .leafout_disabled = false,
        .accumulated_leafout_h = 1,
        .required_leafout_h = 1,
    }};
    try std.testing.expect((try sourceOrderDormantSeedActivation(
        false,
        true,
        &ready,
        1,
        2025,
        0,
    )) == null);
    try std.testing.expect((try sourceOrderDormantSeedActivation(
        true,
        false,
        &ready,
        1,
        2025,
        0,
    )) == null);
    try std.testing.expectError(
        error.InvalidDormantSeedLeafoutHours,
        sourceOrderDormantSeedActivation(
            true,
            true,
            &.{.{
                .leafout_disabled = false,
                .accumulated_leafout_h = std.math.nan(f64),
                .required_leafout_h = 1,
            }},
            1,
            2025,
            0,
        ),
    );
}

test "source-order litterfall accumulation preserves position fraction layer order" {
    var carbon: [20]f64 = undefined;
    var nitrogen: [20]f64 = undefined;
    var phosphorus: [20]f64 = undefined;
    for (0..20) |index| {
        carbon[index] = @floatFromInt(index + 1);
        nitrogen[index] = carbon[index] * 0.1;
        phosphorus[index] = carbon[index] * 0.01;
    }
    const result = try sourceOrderAccumulateLitterfall(std.testing.allocator, .{
        .carbon_g_c_by_position_fraction_layer = &carbon,
        .nitrogen_g_n_by_position_fraction_layer = &nitrogen,
        .phosphorus_g_p_by_position_fraction_layer = &phosphorus,
        .layer_count_including_surface = 2,
        .preceding_cumulative_surface_litter = .{ .carbon_g = 5, .nitrogen_g = 0.5, .phosphorus_g = 0.05 },
        .preceding_hourly_litter = .{ .carbon_g = 7, .nitrogen_g = 0.7, .phosphorus_g = 0.07 },
        .preceding_cumulative_litter = .{ .carbon_g = 11, .nitrogen_g = 1.1, .phosphorus_g = 0.11 },
        .preceding_layer_carbon_g_c = &.{ 13, 17 },
    });
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(f64, 105), result.cumulative_surface_litter.carbon_g);
    try std.testing.expectApproxEqAbs(
        @as(f64, 10.5),
        result.cumulative_surface_litter.nitrogen_g,
        1.0e-14,
    );
    try std.testing.expectEqual(@as(f64, 217), result.hourly_litter.carbon_g);
    try std.testing.expectEqual(@as(f64, 221), result.cumulative_litter.carbon_g);
    try std.testing.expectEqualSlices(f64, &.{ 113, 127 }, result.layer_carbon_g_c);
}

test "source-order litterfall accumulation rejects a late invalid value atomically" {
    var values = [_]f64{0} ** 10;
    values[9] = std.math.nan(f64);
    try std.testing.expectError(
        error.InvalidLitterfallAccumulationInput,
        sourceOrderAccumulateLitterfall(std.testing.allocator, .{
            .carbon_g_c_by_position_fraction_layer = &values,
            .nitrogen_g_n_by_position_fraction_layer = &([_]f64{0} ** 10),
            .phosphorus_g_p_by_position_fraction_layer = &([_]f64{0} ** 10),
            .layer_count_including_surface = 1,
            .preceding_cumulative_surface_litter = .{},
            .preceding_hourly_litter = .{},
            .preceding_cumulative_litter = .{},
            .preceding_layer_carbon_g_c = &.{0},
        }),
    );
}
