const std = @import("std");

const builtin = @import("builtin");

const management = @import("plant_management.zig");

const canopy = @import("../canopy/photosynthesis/photosynthesis.zig");

const phenology = @import("../plant/lifecycle/phenology.zig");

const growth_stages = @import("../plant/lifecycle/growth_stages.zig");

const root_system = @import("../plant/root/plant_root_system.zig");

const root_disturbance = @import("../plant/root/plant_root_disturbance.zig");

const symbiotic_fixation = @import("../canopy/symbiosis/plant_symbiotic_fixation.zig");

const root_litterfall = @import("../plant/root/plant_root_litterfall.zig");

const root_litter_ledger = @import("../plant/root/plant_root_litter_ledger.zig");

const litter_partition = @import("../plant/partition/litter.zig");

const soil_organic = @import("../soil/organic/initialization.zig");

const grid_module = @import("../state/grid.zig");

const carbon_exchange = @import("../canopy/photosynthesis/carbon_exchange.zig");

const shoot_litter_bridge = @import("../plant/growth/shoot_litter_bridge.zig");

const canopy_structure = @import("../canopy/morphology/structure.zig");

const canopy_layers = @import("../canopy/radiation/layer_distribution.zig");

const canopy_biochemistry = @import("../canopy/photosynthesis/biochemistry.zig");

const dormancy = @import("../plant/lifecycle/dormancy.zig");

const grazing_manure = @import("grazing_manure.zig");

const surface_nutrients = @import("../soil/biogeochemistry/organic_matter_fire_exchange.zig");

const spring_reproductive_litterfall = @import("../plant/growth/spring_reproductive_litterfall.zig");

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





pub const TillageElementComposition = struct {
    carbon: [2]f64,
    nitrogen: [2]f64,
    phosphorus: [2]f64,
};




pub fn validateTillageComposition(composition: TillageElementComposition) !void {
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






pub fn validateTillageSlice(values: []const f64) !void {
    for (values) |value|
        if (!std.math.isFinite(value) or value < 0) return error.InvalidTillageBranchState;
}









pub const source_order_standing_dead_component_count: usize = 5;




pub const harvest_product_component_count: usize = 5;





pub fn sumHarvestProductComponents(
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
























pub fn validateSourceOrderRootMass(mass: canopy.ElementalMass) !void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDeadRootLitterfallInput;
    }
}









pub fn validateDeadRootResetMass(mass: canopy.ElementalMass) !void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDeadRootResetInput;
    }
}




pub fn validateDeadNoduleMass(mass: canopy.ElementalMass) !void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDeadNoduleLitterfallInput;
    }
}



pub fn validateDeadRootDepthMass(mass: canopy.ElementalMass) !void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidDeadRootDepthResetInput;
    }
}





pub fn validateCompleteDeathMass(mass: canopy.ElementalMass) !void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field| {
        const value = @field(mass, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidCompleteDeathShootLitterfallInput;
    }
}

pub fn validateCompleteDeathResultMass(mass: canopy.ElementalMass) !void {
    inline for (@typeInfo(canopy.ElementalMass).@"struct".fields) |field|
        if (!std.math.isFinite(@field(mass, field.name)))
            return error.NonFiniteCompleteDeathShootLitterfall;
}































pub fn scaleElementalMass(mass: canopy.ElementalMass, fraction: f64) canopy.ElementalMass {
    return .{
        .carbon_g = mass.carbon_g * fraction,
        .nitrogen_g = mass.nitrogen_g * fraction,
        .phosphorus_g = mass.phosphorus_g * fraction,
    };
}








pub fn validateCombustionNode(node: SourceOrderShootCombustionNodeState) !void {
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

pub fn remainingShootPools(
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

pub fn uncombustedShootTotal(
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














pub fn boundedCombustionFraction(total_g_c: f64, rate_g_c_step: f64, threshold_g_c: f64) f64 {
    return if (total_g_c > threshold_g_c)
        @min(@as(f64, 1), rate_g_c_step / total_g_c)
    else
        0;
}

pub fn subtractElementalMass(
    inventory: canopy.ElementalMass,
    loss: canopy.ElementalMass,
) canopy.ElementalMass {
    return .{
        .carbon_g = inventory.carbon_g - loss.carbon_g,
        .nitrogen_g = inventory.nitrogen_g - loss.nitrogen_g,
        .phosphorus_g = inventory.phosphorus_g - loss.phosphorus_g,
    };
}






pub fn scaleSaltInventory(
    inventory: SourceOrderShootSaltInventory,
    fraction: f64,
) SourceOrderShootSaltInventory {
    var result: SourceOrderShootSaltInventory = undefined;
    inline for (@typeInfo(SourceOrderShootSaltInventory).@"struct".fields) |field|
        @field(result, field.name) = @field(inventory, field.name) * fraction;
    return result;
}

pub fn subtractSaltInventory(
    inventory: SourceOrderShootSaltInventory,
    loss: SourceOrderShootSaltInventory,
) SourceOrderShootSaltInventory {
    var result: SourceOrderShootSaltInventory = undefined;
    inline for (@typeInfo(SourceOrderShootSaltInventory).@"struct".fields) |field|
        @field(result, field.name) = @field(inventory, field.name) - @field(loss, field.name);
    return result;
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
    /// Date assigned to source-generated automatic harvests. This mirrors
    /// IDAYH/IYRH and is distinct from user management schedules.
    automatic_harvest_date_by_plant: ?[]management.PackedDate = null,
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



pub const PopulationAfterDisturbance = struct {
    living_population_per_m2: f64,
    living_population_count: f64,
    standing_dead_population_count: f64,
};


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
    seasonal_turnover_event_by_plant: []const bool,
    growth_habit_by_plant: []const u8,
    leaf_phenology_type_by_plant: []const u8,
    current_date: management.PackedDate,
) !usize {
    const plant_state = context.plant_phenology orelse return error.IncompleteAutomaticSelfSeedingContext;
    const plant_count = context.canopy_state.plant_branch_offsets.len - 1;
    if (growth_habit_by_plant.len != plant_count or leaf_phenology_type_by_plant.len != plant_count or
        plant_state.active.len != plant_count or plant_state.reseed_pending.len != plant_count or
        seasonal_turnover_event_by_plant.len != plant_count or
        (context.automatic_harvest_date_by_plant != null and context.automatic_harvest_date_by_plant.?.len != plant_count))
        return error.AutomaticSelfSeedingDimensionMismatch;
    _ = try current_date.dayOfYear(current_date.year);
    var applied: usize = 0;
    for (0..plant_count) |plant| {
        if (!plant_state.active[plant] or plant_state.reseed_pending[plant] or
            growth_habit_by_plant[plant] != 0 or leaf_phenology_type_by_plant[plant] == 0)
            continue;
        if (!seasonal_turnover_event_by_plant[plant]) continue;
        try applyEvent(context, plant, .{
            .date = current_date,
            .kind = .grain,
            .termination = .terminate_and_reseed,
            .cutting_height_m_or_lai_fraction = 0,
            .thinning_fraction_or_consumption_rate = 0,
            .harvested_fraction = .{ .leaf = 1, .nonfoliar = 1, .woody = 1, .standing_dead = 1 },
            .ecosystem_export_fraction = .{ .leaf = 0, .nonfoliar = 1, .woody = 0, .standing_dead = 0 },
        });
        if (context.automatic_harvest_date_by_plant) |dates| dates[plant] = current_date;
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
    // This transition is rare (once per admitted leafout), so stage the full
    // canopy owner to guarantee that a late topology/arithmetic failure cannot
    // publish a partial C/N/P cleanup.
    var staged_canopy = try context.canopy_state.clone();
    defer staged_canopy.deinit();
    const staged_products = try context.canopy_state.allocator.dupe(ProductLedger, context.products_by_plant);
    defer context.canopy_state.allocator.free(staged_products);
    var staged_context = context.*;
    staged_context.canopy_state = &staged_canopy;
    staged_context.products_by_plant = staged_products;
    try applyStartOfSeasonResidueStaged(&staged_context, plant, biomass_turnover_type, root_profile_type);
    inline for (@typeInfo(canopy.State).@"struct".fields) |field| {
        if (field.type == []f64) @memcpy(@field(context.canopy_state, field.name), @field(staged_canopy, field.name));
    }
    @memcpy(context.products_by_plant, staged_products);
}

fn applyStartOfSeasonResidueStaged(
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
    const reproductive_kinetics = try partitions.get(plant, .non_foliar);
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
        var husk: spring_reproductive_litterfall.Elements = .{ .carbon = state.branch_husk_carbon_g[branch], .nitrogen = state.branch_husk_nitrogen_g[branch], .phosphorus = state.branch_husk_phosphorus_g[branch] };
        var ear: spring_reproductive_litterfall.Elements = .{ .carbon = state.branch_ear_carbon_g[branch], .nitrogen = state.branch_ear_nitrogen_g[branch], .phosphorus = state.branch_ear_phosphorus_g[branch] };
        var grain: spring_reproductive_litterfall.Elements = .{ .carbon = state.branch_grain_carbon_g[branch], .nitrogen = state.branch_grain_nitrogen_g[branch], .phosphorus = state.branch_grain_phosphorus_g[branch] };
        _ = try spring_reproductive_litterfall.apply(.{
            .husk = &husk,
            .ear = &ear,
            .grain = &grain,
            .potential_seed_site_count = &state.branch_potential_seed_site_count[branch],
            .grain_count = &state.branch_seed_count[branch],
            .individual_grain_carbon_g_c = &state.branch_individual_seed_carbon_g[branch],
            .litter_carbon_g_c = &context.products_by_plant[plant].direct_litter.nonwoody_carbon_g,
            .litter_nitrogen_g_n = &context.products_by_plant[plant].direct_litter.nonwoody_nitrogen_g,
            .litter_phosphorus_g_p = &context.products_by_plant[plant].direct_litter.nonwoody_phosphorus_g,
        }, .{
            .leafout_status = .enabled,
            .perennial = true,
            .accumulated_leafout_h = 1,
            .required_leafout_h = 1,
            .reproductive_litter_kinetics = .{ .carbon = &reproductive_kinetics.carbon, .nitrogen = &reproductive_kinetics.nitrogen, .phosphorus = &reproductive_kinetics.phosphorus },
        });
        state.branch_husk_carbon_g[branch] = husk.carbon;
        state.branch_husk_nitrogen_g[branch] = husk.nitrogen;
        state.branch_husk_phosphorus_g[branch] = husk.phosphorus;
        state.branch_ear_carbon_g[branch] = ear.carbon;
        state.branch_ear_nitrogen_g[branch] = ear.nitrogen;
        state.branch_ear_phosphorus_g[branch] = ear.phosphorus;
        state.branch_grain_carbon_g[branch] = grain.carbon;
        state.branch_grain_nitrogen_g[branch] = grain.nitrogen;
        state.branch_grain_phosphorus_g[branch] = grain.phosphorus;

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

    // Whole-plant death publishes reproductive organs through the nonfoliar
    // harvest ledger.  The spring-transition owner below deliberately routes
    // those organs directly to litter, so consume them here before sharing its
    // foliage and stalk cleanup.
    for (branches.first..branches.end) |branch| {
        const reproductive = try canopy.harvestReproductiveOrgans(state, branch, .{
            .husk_remaining = 0,
            .husk_unexported = 1,
            .ear_remaining = 0,
            .ear_unexported = 1,
            .grain_remaining = 0,
            .grain_unexported = 1,
        });
        addProducts(&context.products_by_plant[plant].nonfoliar, reproductive.products);
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
        var host_litter_by_domain = [_]@import("../plant/root/plant_root_metabolism.zig").RootLitter{
            std.mem.zeroes(@import("../plant/root/plant_root_metabolism.zig").RootLitter),
            std.mem.zeroes(@import("../plant/root/plant_root_metabolism.zig").RootLitter),
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
        .litterfall = std.mem.zeroes(@import("../plant/root/plant_root_metabolism.zig").RootLitter),
    };
}

const HostLayerHarvest = struct {
    litterfall: @import("../plant/root/plant_root_metabolism.zig").RootLitter,
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
    var result: HostLayerHarvest = .{ .litterfall = std.mem.zeroes(@import("../plant/root/plant_root_metabolism.zig").RootLitter), .removed_carbon_g_c = 0 };
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

fn addElementPartition(litter: *@import("../plant/root/plant_root_metabolism.zig").RootLitter, mass: canopy.ElementalMass, fractions: litter_partition.ElementFractions, woody: bool, multiplier: f64) void {
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
    const config = try @import("../core/config.zig").SimulationConfig.init(
        .{ .lon_count = 1, .lat_count = 1, .soil_layers = 2, .plant_populations = 1 },
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

test {
    // Keeps the extracted tests discoverable by `zig build test`,
    // which only reaches files reachable by import.
    _ = @import("../validation/plant_harvest_runtime_test.zig");
}

// Moved to plant_harvest_source_order.zig; re-exported so call sites are unchanged.
const __sourceOrder = @import("plant_harvest_source_order.zig");
pub const SourceOrderAbovegroundHarvestLitterInput = __sourceOrder.SourceOrderAbovegroundHarvestLitterInput;
pub const SourceOrderAbovegroundHarvestLitterResult = __sourceOrder.SourceOrderAbovegroundHarvestLitterResult;
pub const SourceOrderCharcoalCombustionResult = __sourceOrder.SourceOrderCharcoalCombustionResult;
pub const SourceOrderColdSoilCombustionReset = __sourceOrder.SourceOrderColdSoilCombustionReset;
pub const SourceOrderCombustionRates = __sourceOrder.SourceOrderCombustionRates;
pub const SourceOrderCombustionSpecificRates = __sourceOrder.SourceOrderCombustionSpecificRates;
pub const SourceOrderCompleteDeathBranchPools = __sourceOrder.SourceOrderCompleteDeathBranchPools;
pub const SourceOrderCompleteDeathBranchState = __sourceOrder.SourceOrderCompleteDeathBranchState;
pub const SourceOrderCompleteDeathRootInput = __sourceOrder.SourceOrderCompleteDeathRootInput;
pub const SourceOrderCompleteDeathRootResetState = __sourceOrder.SourceOrderCompleteDeathRootResetState;
pub const SourceOrderCompleteDeathShootInput = __sourceOrder.SourceOrderCompleteDeathShootInput;
pub const SourceOrderCompleteDeathShootResult = __sourceOrder.SourceOrderCompleteDeathShootResult;
pub const SourceOrderDeadBranchCanopyReset = __sourceOrder.SourceOrderDeadBranchCanopyReset;
pub const SourceOrderDeadBranchLayerResetState = __sourceOrder.SourceOrderDeadBranchLayerResetState;
pub const SourceOrderDeadBranchLitterInput = __sourceOrder.SourceOrderDeadBranchLitterInput;
pub const SourceOrderDeadBranchLitterPools = __sourceOrder.SourceOrderDeadBranchLitterPools;
pub const SourceOrderDeadBranchLitterResult = __sourceOrder.SourceOrderDeadBranchLitterResult;
pub const SourceOrderDeadBranchNodeResetState = __sourceOrder.SourceOrderDeadBranchNodeResetState;
pub const SourceOrderDeadBranchPhenologyState = __sourceOrder.SourceOrderDeadBranchPhenologyState;
pub const SourceOrderDeadBranchResetInput = __sourceOrder.SourceOrderDeadBranchResetInput;
pub const SourceOrderDeadBranchScalarResetState = __sourceOrder.SourceOrderDeadBranchScalarResetState;
pub const SourceOrderDeadBranchStorageRecoveryInput = __sourceOrder.SourceOrderDeadBranchStorageRecoveryInput;
pub const SourceOrderDeadBranchStorageRecoveryResult = __sourceOrder.SourceOrderDeadBranchStorageRecoveryResult;
pub const SourceOrderDeadNoduleLayerPools = __sourceOrder.SourceOrderDeadNoduleLayerPools;
pub const SourceOrderDeadNoduleLitterInput = __sourceOrder.SourceOrderDeadNoduleLitterInput;
pub const SourceOrderDeadPerennialReseedResult = __sourceOrder.SourceOrderDeadPerennialReseedResult;
pub const SourceOrderDeadRootAxisLayerState = __sourceOrder.SourceOrderDeadRootAxisLayerState;
pub const SourceOrderDeadRootAxisPools = __sourceOrder.SourceOrderDeadRootAxisPools;
pub const SourceOrderDeadRootDepthResetState = __sourceOrder.SourceOrderDeadRootDepthResetState;
pub const SourceOrderDeadRootDomainAxisState = __sourceOrder.SourceOrderDeadRootDomainAxisState;
pub const SourceOrderDeadRootDomainLayerState = __sourceOrder.SourceOrderDeadRootDomainLayerState;
pub const SourceOrderDeadRootLayerLitter = __sourceOrder.SourceOrderDeadRootLayerLitter;
pub const SourceOrderDeadRootLitterInput = __sourceOrder.SourceOrderDeadRootLitterInput;
pub const SourceOrderDeadRootResetState = __sourceOrder.SourceOrderDeadRootResetState;
pub const SourceOrderDisturbanceRemovalInput = __sourceOrder.SourceOrderDisturbanceRemovalInput;
pub const SourceOrderDisturbanceRemovalResult = __sourceOrder.SourceOrderDisturbanceRemovalResult;
pub const SourceOrderDormantSeedActivation = __sourceOrder.SourceOrderDormantSeedActivation;
pub const SourceOrderDormantSeedBranch = __sourceOrder.SourceOrderDormantSeedBranch;
pub const SourceOrderFireInventoryResult = __sourceOrder.SourceOrderFireInventoryResult;
pub const SourceOrderFireLayerCarbon = __sourceOrder.SourceOrderFireLayerCarbon;
pub const SourceOrderFirePlantInventory = __sourceOrder.SourceOrderFirePlantInventory;
pub const SourceOrderFireShootCarbon = __sourceOrder.SourceOrderFireShootCarbon;
pub const SourceOrderGrazingLitterLedgerState = __sourceOrder.SourceOrderGrazingLitterLedgerState;
pub const SourceOrderGrazingLitterResult = __sourceOrder.SourceOrderGrazingLitterResult;
pub const SourceOrderHarvestResidueInput = __sourceOrder.SourceOrderHarvestResidueInput;
pub const SourceOrderLitterfallAccumulationInput = __sourceOrder.SourceOrderLitterfallAccumulationInput;
pub const SourceOrderLitterfallAccumulationResult = __sourceOrder.SourceOrderLitterfallAccumulationResult;
pub const SourceOrderNoCombustionReset = __sourceOrder.SourceOrderNoCombustionReset;
pub const SourceOrderReseedDate = __sourceOrder.SourceOrderReseedDate;
pub const SourceOrderRootAxisCombustionReset = __sourceOrder.SourceOrderRootAxisCombustionReset;
pub const SourceOrderRootCombustionAxisState = __sourceOrder.SourceOrderRootCombustionAxisState;
pub const SourceOrderRootCombustionDomainLoss = __sourceOrder.SourceOrderRootCombustionDomainLoss;
pub const SourceOrderRootCombustionDomainState = __sourceOrder.SourceOrderRootCombustionDomainState;
pub const SourceOrderRootCombustionFractions = __sourceOrder.SourceOrderRootCombustionFractions;
pub const SourceOrderRootCombustionLayerTotals = __sourceOrder.SourceOrderRootCombustionLayerTotals;
pub const SourceOrderRootCombustionPotentialRates = __sourceOrder.SourceOrderRootCombustionPotentialRates;
pub const SourceOrderRootDomainCombustionResult = __sourceOrder.SourceOrderRootDomainCombustionResult;
pub const SourceOrderRootGasInventory = __sourceOrder.SourceOrderRootGasInventory;
pub const SourceOrderRootGasPhases = __sourceOrder.SourceOrderRootGasPhases;
pub const SourceOrderRootNoduleCombustionResult = __sourceOrder.SourceOrderRootNoduleCombustionResult;
pub const SourceOrderRootStorageCombustionInput = __sourceOrder.SourceOrderRootStorageCombustionInput;
pub const SourceOrderRootStorageCombustionResult = __sourceOrder.SourceOrderRootStorageCombustionResult;
pub const SourceOrderShootCombustionBranchPools = __sourceOrder.SourceOrderShootCombustionBranchPools;
pub const SourceOrderShootCombustionBranchResult = __sourceOrder.SourceOrderShootCombustionBranchResult;
pub const SourceOrderShootCombustionFractions = __sourceOrder.SourceOrderShootCombustionFractions;
pub const SourceOrderShootCombustionNodeLayerState = __sourceOrder.SourceOrderShootCombustionNodeLayerState;
pub const SourceOrderShootCombustionNodeState = __sourceOrder.SourceOrderShootCombustionNodeState;
pub const SourceOrderShootCombustionResult = __sourceOrder.SourceOrderShootCombustionResult;
pub const SourceOrderShootCombustionTotals = __sourceOrder.SourceOrderShootCombustionTotals;
pub const SourceOrderShootSaltCombustionBranchResult = __sourceOrder.SourceOrderShootSaltCombustionBranchResult;
pub const SourceOrderShootSaltInventory = __sourceOrder.SourceOrderShootSaltInventory;
pub const SourceOrderSoilPlantExchangeInput = __sourceOrder.SourceOrderSoilPlantExchangeInput;
pub const SourceOrderSoilPlantExchangeResult = __sourceOrder.SourceOrderSoilPlantExchangeResult;
pub const SourceOrderStandingDeadCombustionComponent = __sourceOrder.SourceOrderStandingDeadCombustionComponent;
pub const SourceOrderStandingDeadCombustionResult = __sourceOrder.SourceOrderStandingDeadCombustionResult;
pub const SourceOrderStandingDeadGeometryInput = __sourceOrder.SourceOrderStandingDeadGeometryInput;
pub const SourceOrderStandingDeadGeometryResult = __sourceOrder.SourceOrderStandingDeadGeometryResult;
pub const SourceOrderStandingDeadHarvestInput = __sourceOrder.SourceOrderStandingDeadHarvestInput;
pub const SourceOrderStandingDeadHarvestResult = __sourceOrder.SourceOrderStandingDeadHarvestResult;
pub const SourceOrderTillageBranchLitterInput = __sourceOrder.SourceOrderTillageBranchLitterInput;
pub const SourceOrderTillageBranchLitterResult = __sourceOrder.SourceOrderTillageBranchLitterResult;
pub const SourceOrderTillageBranchPools = __sourceOrder.SourceOrderTillageBranchPools;
pub const SourceOrderTillageBranchRetentionResult = __sourceOrder.SourceOrderTillageBranchRetentionResult;
pub const SourceOrderTillageBranchScalarState = __sourceOrder.SourceOrderTillageBranchScalarState;
pub const SourceOrderTillageLayerSampleState = __sourceOrder.SourceOrderTillageLayerSampleState;
pub const SourceOrderTillageNodeState = __sourceOrder.SourceOrderTillageNodeState;
pub const SourceOrderTillagePopulationInput = __sourceOrder.SourceOrderTillagePopulationInput;
pub const SourceOrderTillagePopulationResult = __sourceOrder.SourceOrderTillagePopulationResult;
pub const SourceOrderTillagePopulationState = __sourceOrder.SourceOrderTillagePopulationState;
pub const SourceOrderTillageStandingDeadInput = __sourceOrder.SourceOrderTillageStandingDeadInput;
pub const SourceOrderTillageStandingDeadResult = __sourceOrder.SourceOrderTillageStandingDeadResult;
pub const SourceOrderTillageTerminationInput = __sourceOrder.SourceOrderTillageTerminationInput;
pub const SourceOrderTillageTerminationResult = __sourceOrder.SourceOrderTillageTerminationResult;
pub const SourceOrderTillageTerminationState = __sourceOrder.SourceOrderTillageTerminationState;
pub const SourceOrderUncombustedBranchState = __sourceOrder.SourceOrderUncombustedBranchState;
pub const SourceOrderWholePlantTerminationResult = __sourceOrder.SourceOrderWholePlantTerminationResult;
pub const SourceOrderWholePlantTerminationState = __sourceOrder.SourceOrderWholePlantTerminationState;
pub const sourceOrderAbovegroundDisturbanceIsEnabled = __sourceOrder.sourceOrderAbovegroundDisturbanceIsEnabled;
pub const sourceOrderAbovegroundHarvestLitter = __sourceOrder.sourceOrderAbovegroundHarvestLitter;
pub const sourceOrderAccumulateLitterfall = __sourceOrder.sourceOrderAccumulateLitterfall;
pub const sourceOrderAccumulateSoilPlantExchange = __sourceOrder.sourceOrderAccumulateSoilPlantExchange;
pub const sourceOrderAggregateFireCarbonInventory = __sourceOrder.sourceOrderAggregateFireCarbonInventory;
pub const sourceOrderApplyRootDomainCombustion = __sourceOrder.sourceOrderApplyRootDomainCombustion;
pub const sourceOrderApplyUncombustedShootState = __sourceOrder.sourceOrderApplyUncombustedShootState;
pub const sourceOrderCharcoalCombustion = __sourceOrder.sourceOrderCharcoalCombustion;
pub const sourceOrderCombustionRates = __sourceOrder.sourceOrderCombustionRates;
pub const sourceOrderCompleteDeathRootLitterfall = __sourceOrder.sourceOrderCompleteDeathRootLitterfall;
pub const sourceOrderCompleteDeathShootLitterfall = __sourceOrder.sourceOrderCompleteDeathShootLitterfall;
pub const sourceOrderCuttingHeightFromLeafAreaRemoval = __sourceOrder.sourceOrderCuttingHeightFromLeafAreaRemoval;
pub const sourceOrderDeadBranchLitterfall = __sourceOrder.sourceOrderDeadBranchLitterfall;
pub const sourceOrderDeadBranchStorageRecovery = __sourceOrder.sourceOrderDeadBranchStorageRecovery;
pub const sourceOrderDeadNoduleLitterfall = __sourceOrder.sourceOrderDeadNoduleLitterfall;
pub const sourceOrderDeadRootLitterfall = __sourceOrder.sourceOrderDeadRootLitterfall;
pub const sourceOrderDormantSeedActivation = __sourceOrder.sourceOrderDormantSeedActivation;
pub const sourceOrderForestSelfThinningIsEnabled = __sourceOrder.sourceOrderForestSelfThinningIsEnabled;
pub const sourceOrderGrazingLitterLedgers = __sourceOrder.sourceOrderGrazingLitterLedgers;
pub const sourceOrderHarvestResidueRouting = __sourceOrder.sourceOrderHarvestResidueRouting;
pub const sourceOrderHourlyDisturbanceReset = __sourceOrder.sourceOrderHourlyDisturbanceReset;
pub const sourceOrderPopulationAfterDisturbance = __sourceOrder.sourceOrderPopulationAfterDisturbance;
pub const sourceOrderPopulationScaledNumericalThresholds = __sourceOrder.sourceOrderPopulationScaledNumericalThresholds;
pub const sourceOrderReleaseDeadRootGases = __sourceOrder.sourceOrderReleaseDeadRootGases;
pub const sourceOrderResetColdSoilCombustion = __sourceOrder.sourceOrderResetColdSoilCombustion;
pub const sourceOrderResetCompleteDeathBranches = __sourceOrder.sourceOrderResetCompleteDeathBranches;
pub const sourceOrderResetCompleteDeathRoots = __sourceOrder.sourceOrderResetCompleteDeathRoots;
pub const sourceOrderResetDeadBranchCanopy = __sourceOrder.sourceOrderResetDeadBranchCanopy;
pub const sourceOrderResetDeadBranchPhenology = __sourceOrder.sourceOrderResetDeadBranchPhenology;
pub const sourceOrderResetDeadRootDepth = __sourceOrder.sourceOrderResetDeadRootDepth;
pub const sourceOrderResetDeadRootState = __sourceOrder.sourceOrderResetDeadRootState;
pub const sourceOrderResetNoCombustion = __sourceOrder.sourceOrderResetNoCombustion;
pub const sourceOrderRetainTillageBranchState = __sourceOrder.sourceOrderRetainTillageBranchState;
pub const sourceOrderRootNoduleCombustion = __sourceOrder.sourceOrderRootNoduleCombustion;
pub const sourceOrderRootStorageCombustion = __sourceOrder.sourceOrderRootStorageCombustion;
pub const sourceOrderScheduleDeadPerennialReseed = __sourceOrder.sourceOrderScheduleDeadPerennialReseed;
pub const sourceOrderShootCombustionFractions = __sourceOrder.sourceOrderShootCombustionFractions;
pub const sourceOrderShootCombustionLosses = __sourceOrder.sourceOrderShootCombustionLosses;
pub const sourceOrderShootSaltCombustion = __sourceOrder.sourceOrderShootSaltCombustion;
pub const sourceOrderStandingDeadCombustion = __sourceOrder.sourceOrderStandingDeadCombustion;
pub const sourceOrderStandingDeadGeometry = __sourceOrder.sourceOrderStandingDeadGeometry;
pub const sourceOrderStandingDeadHarvest = __sourceOrder.sourceOrderStandingDeadHarvest;
pub const sourceOrderTillageBranchLitter = __sourceOrder.sourceOrderTillageBranchLitter;
pub const sourceOrderTillagePopulationReduction = __sourceOrder.sourceOrderTillagePopulationReduction;
pub const sourceOrderTillageStandingDead = __sourceOrder.sourceOrderTillageStandingDead;
pub const sourceOrderTillageTermination = __sourceOrder.sourceOrderTillageTermination;
pub const sourceOrderTotalDisturbanceRemoval = __sourceOrder.sourceOrderTotalDisturbanceRemoval;
pub const sourceOrderWholePlantTermination = __sourceOrder.sourceOrderWholePlantTermination;
