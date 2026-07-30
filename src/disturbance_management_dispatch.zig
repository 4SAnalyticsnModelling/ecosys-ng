const std = @import("std");
const delimited_input = @import("delimited_input.zig");
const disturbance = @import("disturbance_schedule.zig");
const land_management = @import("land_management.zig");
const Date = @import("options.zig").Date;
const RootSystem = @import("plant_root_system.zig");
const RootDisturbance = @import("plant_root_disturbance.zig");
const RootLitterfall = @import("plant_root_litterfall.zig");
const RootLitterLedger = @import("plant_root_litter_ledger.zig");
const RootMetabolism = @import("plant_root_metabolism.zig");
const LitterPartition = @import("plant_litter_partition.zig");
const SoilOrganic = @import("soil_organic_initialization.zig");
const Grid = @import("grid.zig").GridState;
const SurfaceEnergy = @import("surface_energy.zig").State;
const Canopy = @import("canopy_photosynthesis.zig");
const OrganicMatterFireExchange = @import("organic_matter_fire_exchange.zig");
const PlantHarvest = @import("plant_harvest_runtime.zig");
const execution_calendar_date = @import("execution_calendar_date.zig");

pub const ScheduleMap = struct {
    allocator: std.mem.Allocator,
    catalog_index_by_cell: []?usize,

    pub fn init(allocator: std.mem.Allocator, assignments: land_management.Assignments, unit_by_cell: []const usize, catalog: disturbance.Catalog) !ScheduleMap {
        const map = try allocator.alloc(?usize, unit_by_cell.len);
        errdefer allocator.free(map);
        for (unit_by_cell, 0..) |unit_index, cell| {
            if (unit_index >= assignments.units.len) return error.LandManagementUnitIndexOutOfBounds;
            const name = assignments.units[unit_index].tillage_file;
            map[cell] = if (delimited_input.isNo(name)) null else catalog.find(name) orelse return error.DisturbanceScheduleMissingFromCatalog;
        }
        return .{ .allocator = allocator, .catalog_index_by_cell = map };
    }

    pub fn deinit(self: *ScheduleMap) void {
        self.allocator.free(self.catalog_index_by_cell);
        self.* = undefined;
    }
};

pub const ApplyContext = struct {
    roots: *RootSystem.State,
    litter_partition: *const LitterPartition.State,
    soil_organic: *SoilOrganic.State,
    grid: *const Grid,
    species_count: usize,
    biological_domain_count_by_plant: []const u8,
    root_nonwoody_fraction_by_plant: []const f64,
    biomass_turnover_type_by_plant: []const u8,
    root_profile_type_by_plant: []const u8,
    growth_habit_by_plant: []const u8,
    leaf_phenology_type_by_plant: []const u8,
    planting_day_of_year_by_plant: []const u16,
    planting_year_by_plant: []const i32,
    current_day_of_year: u16,
    current_year: i32,
    plant_harvest: ?*PlantHarvest.Context,
    root_litter_carbon_ledger: ?*RootLitterLedger.State = null,
    surface_energy: *SurfaceEnergy,
    fire_active_this_hour: []bool,
};

pub fn dispatchDate(map: ScheduleMap, catalog: disturbance.Catalog, date: Date, context: *ApplyContext) !usize {
    const before = try dispatchDatePhase(map, catalog, date, context, .pre_science);
    return before + try dispatchDatePhase(map, catalog, date, context, .post_science);
}

pub const Phase = enum { pre_science, post_science };

pub fn dispatchDatePhase(map: ScheduleMap, catalog: disturbance.Catalog, date: Date, context: *ApplyContext, phase: Phase) !usize {
    try validateDispatchDate(date);
    var applied: usize = 0;
    for (map.catalog_index_by_cell, 0..) |maybe_schedule, cell| {
        const schedule = maybe_schedule orelse continue;
        if (schedule >= catalog.entries.items.len) return error.DisturbanceScheduleIndexOutOfBounds;
        for (catalog.entries.items[schedule].events) |event| {
            if (event.date.day != date.day or event.date.month != date.month or (!event.date.isRecurring() and event.date.year != date.year)) continue;
            const event_phase: Phase = switch (event.operation) {
                .fire => .pre_science,
                else => .post_science,
            };
            if (event_phase != phase) continue;
            try applyEvent(context, cell, event);
            applied += 1;
        }
    }
    return applied;
}

fn validateDispatchDate(date: Date) !void {
    if (date.year == 0) return error.InvalidDisturbanceDispatchDate;
    _ = execution_calendar_date.dayOfYear(.{ .day = date.day, .month = date.month, .year = date.year }) catch return error.InvalidDisturbanceDispatchDate;
}

pub fn applyEvent(context: *ApplyContext, cell: usize, event: disturbance.Event) !void {
    const tillage = switch (event.operation) {
        .tillage => |value| value,
        .fire => |fire| {
            if (cell >= context.surface_energy.cell_count or cell >= context.fire_active_this_hour.len) return error.DisturbanceCellOutOfBounds;
            const energy_mj_per_m2 = 3.6 * fire.energy_kw_per_m2;
            const next = context.surface_energy.fire_ignition_mj_per_m2[cell] + energy_mj_per_m2;
            if (!std.math.isFinite(next)) return error.NonFiniteFireIgnitionEnergy;
            context.surface_energy.fire_ignition_mj_per_m2[cell] = next;
            context.fire_active_this_hour[cell] = true;
            return;
        },
        else => return,
    };
    if (context.species_count == 0 or context.roots.plant_count != context.grid.cell_count * context.species_count or context.litter_partition.plant_count != context.roots.plant_count or
        context.biological_domain_count_by_plant.len != context.roots.plant_count or context.root_nonwoody_fraction_by_plant.len != context.roots.plant_count or
        context.biomass_turnover_type_by_plant.len != context.roots.plant_count or context.root_profile_type_by_plant.len != context.roots.plant_count or
        context.growth_habit_by_plant.len != context.roots.plant_count or context.leaf_phenology_type_by_plant.len != context.roots.plant_count or
        context.planting_day_of_year_by_plant.len != context.roots.plant_count or context.planting_year_by_plant.len != context.roots.plant_count)
        return error.DisturbancePlantDimensionMismatch;
    if (cell >= context.grid.cell_count) return error.DisturbanceCellOutOfBounds;
    const retention = RootDisturbance.ElementRetention.uniform(1 - tillage.mixing_fraction);
    for (0..context.species_count) |species| {
        if (!tillage.includes_crop and species == 0) continue;
        const plant = cell * context.species_count + species;
        const after_planting = context.current_year > context.planting_year_by_plant[plant] or
            (context.current_year == context.planting_year_by_plant[plant] and
                context.current_day_of_year > context.planting_day_of_year_by_plant[plant]);
        const tillage_eligible = (context.biomass_turnover_type_by_plant[plant] == 0 or
            context.root_profile_type_by_plant[plant] <= 1) and after_planting;
        if (!tillage_eligible) continue;
        const fine = try context.litter_partition.get(plant, .fine_root);
        const coarse = try context.litter_partition.get(plant, .coarse_wood);
        const mobile = try context.litter_partition.get(plant, .nonstructural);
        const nonwoody_fraction = context.root_nonwoody_fraction_by_plant[plant];
        if (!std.math.isFinite(nonwoody_fraction) or nonwoody_fraction < 0 or nonwoody_fraction > 1) return error.InvalidRootTillageInput;
        for (0..context.grid.active_soil_layer_count[cell]) |layer| {
            const root = try context.roots.layerIndex(plant, 0, layer);
            const result = try calculateLayer(context.roots, root, retention, fine, mobile);
            var publication: RootLitterfall.LayerInput = .{};
            try publication.add(result.litterfall);
            try publication.add(try calculateHostRootTillage(
                context.roots,
                plant,
                layer,
                context.biological_domain_count_by_plant[plant],
                tillage.mixing_fraction,
                nonwoody_fraction,
                coarse,
                fine,
                mobile,
            ));
            if (context.root_litter_carbon_ledger) |ledger| {
                const host_domain_zero = try calculateHostRootTillageDomain(context.roots, plant, 0, layer, tillage.mixing_fraction, nonwoody_fraction, coarse, fine, mobile);
                try ledger.validateAdd(plant, 0, layer, host_domain_zero);
                if (context.biological_domain_count_by_plant[plant] > 1) {
                    const host_domain_one = try calculateHostRootTillageDomain(context.roots, plant, 1, layer, tillage.mixing_fraction, nonwoody_fraction, coarse, fine, mobile);
                    try ledger.validateCarbonAdd(
                        plant,
                        1,
                        layer,
                        try RootLitterLedger.totalCarbon(host_domain_one) +
                            try RootLitterLedger.totalCarbon(result.litterfall),
                    );
                } else {
                    try ledger.validateAdd(plant, 1, layer, result.litterfall);
                }
            }
            try RootDisturbance.validateRootGasRelease(context.roots, plant, layer, tillage.mixing_fraction);
            try RootLitterfall.validatePublication(context.soil_organic, try context.grid.layerIndex(cell, layer), publication);
        }
        const harvest = context.plant_harvest orelse return error.IncompleteAbovegroundTillageContext;
        try PlantHarvest.applyAbovegroundTillage(
            harvest,
            plant,
            1 - tillage.mixing_fraction,
            context.growth_habit_by_plant[plant] == 0 and context.leaf_phenology_type_by_plant[plant] != 0,
        );
        for (0..context.grid.active_soil_layer_count[cell]) |layer| {
            const root = try context.roots.layerIndex(plant, 0, layer);
            const soil = try context.grid.layerIndex(cell, layer);
            const result = try calculateLayer(context.roots, root, retention, fine, mobile);
            var host_litter_by_domain = [_]RootMetabolism.RootLitter{
                std.mem.zeroes(RootMetabolism.RootLitter),
                std.mem.zeroes(RootMetabolism.RootLitter),
            };
            for (0..context.biological_domain_count_by_plant[plant]) |domain|
                host_litter_by_domain[domain] = calculateHostRootTillageDomain(
                    context.roots,
                    plant,
                    domain,
                    layer,
                    tillage.mixing_fraction,
                    nonwoody_fraction,
                    coarse,
                    fine,
                    mobile,
                ) catch unreachable;
            context.roots.symbiont_structural_carbon_g_c[root] = result.structural.carbon_g_c;
            context.roots.symbiont_structural_nitrogen_g_n[root] = result.structural.nitrogen_g_n;
            context.roots.symbiont_structural_phosphorus_g_p[root] = result.structural.phosphorus_g_p;
            context.roots.symbiont_mobile_carbon_g_c[root] = result.mobile.carbon_g_c;
            context.roots.symbiont_mobile_nitrogen_g_n[root] = result.mobile.nitrogen_g_n;
            context.roots.symbiont_mobile_phosphorus_g_p[root] = result.mobile.phosphorus_g_p;
            var publication: RootLitterfall.LayerInput = .{};
            publication.add(result.litterfall) catch unreachable;
            publication.add(calculateHostRootTillage(
                context.roots,
                plant,
                layer,
                context.biological_domain_count_by_plant[plant],
                tillage.mixing_fraction,
                nonwoody_fraction,
                coarse,
                fine,
                mobile,
            ) catch unreachable) catch unreachable;
            commitHostRootTillage(
                context.roots,
                plant,
                layer,
                context.biological_domain_count_by_plant[plant],
                tillage.mixing_fraction,
            ) catch unreachable;
            RootDisturbance.releaseRootGasFraction(context.roots, plant, layer, tillage.mixing_fraction) catch unreachable;
            RootLitterfall.publishValidated(context.soil_organic, soil, publication);
            if (context.root_litter_carbon_ledger) |ledger| {
                ledger.addValidated(plant, 1, layer, result.litterfall);
                for (0..context.biological_domain_count_by_plant[plant]) |domain|
                    ledger.addValidated(plant, domain, layer, host_litter_by_domain[domain]);
            }
        }
    }
}

fn calculateHostRootTillage(
    roots: *const RootSystem.State,
    plant: usize,
    layer: usize,
    domain_count: usize,
    removed_fraction: f64,
    nonwoody_fraction: f64,
    coarse: LitterPartition.ElementFractions,
    fine: LitterPartition.ElementFractions,
    mobile: LitterPartition.ElementFractions,
) !RootMetabolism.RootLitter {
    return calculateHostRootTillageRange(
        roots,
        plant,
        0,
        domain_count,
        layer,
        removed_fraction,
        nonwoody_fraction,
        coarse,
        fine,
        mobile,
    );
}

fn calculateHostRootTillageRange(
    roots: *const RootSystem.State,
    plant: usize,
    first_domain: usize,
    end_domain: usize,
    layer: usize,
    removed_fraction: f64,
    nonwoody_fraction: f64,
    coarse: LitterPartition.ElementFractions,
    fine: LitterPartition.ElementFractions,
    mobile: LitterPartition.ElementFractions,
) !RootMetabolism.RootLitter {
    if (first_domain >= end_domain or end_domain > RootSystem.biological_domain_count or
        !std.math.isFinite(removed_fraction) or removed_fraction < 0 or removed_fraction > 1 or
        !std.math.isFinite(nonwoody_fraction) or nonwoody_fraction < 0 or nonwoody_fraction > 1)
        return error.InvalidRootTillageInput;
    try coarse.validate();
    try fine.validate();
    try mobile.validate();
    var structural = [3]f64{ 0, 0, 0 };
    var mobile_pool = [3]f64{ 0, 0, 0 };
    for (first_domain..end_domain) |domain| {
        const root = try roots.layerIndex(plant, domain, layer);
        inline for (.{
            .{ "mobile_carbon_g", 0 },
            .{ "mobile_nitrogen_g", 1 },
            .{ "mobile_phosphorus_g", 2 },
        }) |entry| {
            const value = @field(roots, entry[0])[root];
            if (!std.math.isFinite(value) or value < 0) return error.InvalidRootTillageInput;
            mobile_pool[entry[1]] += value;
        }
        for (0..roots.active_root_axis_count[plant]) |axis| {
            const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
            inline for (.{
                .{ "axis_primary_carbon_g", "axis_secondary_carbon_g", 0 },
                .{ "axis_primary_nitrogen_g", "axis_secondary_nitrogen_g", 1 },
                .{ "axis_primary_phosphorus_g", "axis_secondary_phosphorus_g", 2 },
            }) |entry| {
                const primary = @field(roots, entry[0])[axis_layer];
                const secondary = @field(roots, entry[1])[axis_layer];
                if (!std.math.isFinite(primary) or primary < 0 or !std.math.isFinite(secondary) or secondary < 0)
                    return error.InvalidRootTillageInput;
                structural[entry[2]] += primary + secondary;
            }
        }
    }
    var litter = std.mem.zeroes(RootMetabolism.RootLitter);
    const woody_fraction = 1 - nonwoody_fraction;
    for (0..LitterPartition.kinetic_component_count) |component| {
        litter.woody_carbon_g_c[component] = removed_fraction * structural[0] * woody_fraction * coarse.carbon[component];
        litter.woody_nitrogen_g_n[component] = removed_fraction * structural[1] * woody_fraction * coarse.nitrogen[component];
        litter.woody_phosphorus_g_p[component] = removed_fraction * structural[2] * woody_fraction * coarse.phosphorus[component];
        litter.nonwoody_carbon_g_c[component] = removed_fraction * (structural[0] * nonwoody_fraction * fine.carbon[component] + mobile_pool[0] * mobile.carbon[component]);
        litter.nonwoody_nitrogen_g_n[component] = removed_fraction * (structural[1] * nonwoody_fraction * fine.nitrogen[component] + mobile_pool[1] * mobile.nitrogen[component]);
        litter.nonwoody_phosphorus_g_p[component] = removed_fraction * (structural[2] * nonwoody_fraction * fine.phosphorus[component] + mobile_pool[2] * mobile.phosphorus[component]);
    }
    return litter;
}

fn calculateHostRootTillageDomain(
    roots: *const RootSystem.State,
    plant: usize,
    domain: usize,
    layer: usize,
    removed_fraction: f64,
    nonwoody_fraction: f64,
    coarse: LitterPartition.ElementFractions,
    fine: LitterPartition.ElementFractions,
    mobile: LitterPartition.ElementFractions,
) !RootMetabolism.RootLitter {
    if (domain >= RootSystem.biological_domain_count)
        return error.InvalidRootTillageInput;
    return calculateHostRootTillageRange(
        roots,
        plant,
        domain,
        domain + 1,
        layer,
        removed_fraction,
        nonwoody_fraction,
        coarse,
        fine,
        mobile,
    );
}

fn commitHostRootTillage(
    roots: *RootSystem.State,
    plant: usize,
    layer: usize,
    domain_count: usize,
    removed_fraction: f64,
) !void {
    const retained = 1 - removed_fraction;
    for (0..domain_count) |domain| {
        const root = try roots.layerIndex(plant, domain, layer);
        inline for (.{
            "mobile_carbon_g",
            "mobile_nitrogen_g",
            "mobile_phosphorus_g",
            "protein_carbon_g",
            "total_carbon_g",
            "primary_root_carbon_g",
            "projected_area_m2",
            "active_length_m",
            "aqueous_volume_m3",
            "gaseous_volume_m3",
            "root_length_m_per_plant",
            "root_length_density_m_per_m3",
            "root_surface_area_m2_per_plant",
            "average_secondary_length_m",
            "symbiotic_respiration_actual_g_c_per_h",
            "symbiotic_respiration_oxygen_unlimited_g_c_per_h",
        }) |field_name| @field(roots, field_name)[root] *= retained;
        for (0..roots.active_root_axis_count[plant]) |axis| {
            const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
            inline for (.{
                "axis_primary_carbon_g",
                "axis_primary_nitrogen_g",
                "axis_primary_phosphorus_g",
                "axis_secondary_carbon_g",
                "axis_secondary_nitrogen_g",
                "axis_secondary_phosphorus_g",
                "axis_primary_length_m",
                "axis_secondary_length_m",
                "axis_primary_count",
                "axis_secondary_count",
            }) |field_name| @field(roots, field_name)[axis_layer] *= retained;
        }
    }
}

/// GROSUB nodule-fire transaction. Source FWPODL/FWTNDL fractions are formed
/// from the total mobile/structural nodule carbon across every runtime species
/// in a cell and layer, then applied uniformly to each species' C/N/P pools.
pub fn applyRootFireCombustion(
    roots: *RootSystem.State,
    canopy: ?*Canopy.State,
    grid: *const Grid,
    species_count: usize,
    biological_domain_count_by_plant: []const u8,
    cell_area_m2: []const f64,
    fire_active_this_hour: []const bool,
    timestep_h: f64,
    dynamic_salts: bool,
    parameters: RootDisturbance.CombustionParameters,
    fire_exchange: *OrganicMatterFireExchange.State,
) !void {
    try parameters.validate();
    if (species_count == 0 or roots.plant_count != grid.cell_count * species_count or biological_domain_count_by_plant.len != roots.plant_count) return error.DisturbancePlantDimensionMismatch;
    if (canopy) |state| if (state.cell_count != grid.cell_count or state.species_count != species_count) return error.DisturbancePlantDimensionMismatch;
    if (cell_area_m2.len != grid.cell_count or fire_active_this_hour.len != grid.cell_count or fire_exchange.layer_count != grid.layer_count) return error.DisturbanceCellDimensionMismatch;
    if (!std.math.isFinite(timestep_h) or timestep_h <= 0) return error.InvalidRootNoduleCombustionInput;
    for (0..grid.cell_count) |cell| {
        if (!fire_active_this_hour[cell]) continue;
        if (!std.math.isFinite(cell_area_m2[cell]) or cell_area_m2[cell] <= 0) return error.InvalidRootNoduleCombustionInput;
        for (0..grid.active_soil_layer_count[cell]) |layer| {
            const soil = try grid.layerIndex(cell, layer);
            var previous_carbon_loss_g_c: f64 = 0;
            var previous_nitrogen_loss_g_n: f64 = 0;
            var previous_phosphorus_loss_g_p: f64 = 0;
            var previous_salt_loss_mol: [OrganicMatterFireExchange.salt_species_count]f64 = @splat(0);
            for (0..species_count) |species| {
                const plant = cell * species_count + species;
                previous_carbon_loss_g_c += roots.combustion_carbon_loss_g_c_per_h[plant];
                previous_nitrogen_loss_g_n += roots.combustion_nitrogen_loss_g_n_per_h[plant];
                previous_phosphorus_loss_g_p += roots.combustion_phosphorus_loss_g_p_per_h[plant];
                const biological_domain_count = try domainCount(biological_domain_count_by_plant, plant);
                if (dynamic_salts) for (0..biological_domain_count) |domain| {
                    const root = try roots.layerIndex(plant, domain, layer);
                    for (0..OrganicMatterFireExchange.salt_species_count) |salt| previous_salt_loss_mol[salt] += roots.combustion_salt_loss_mol_per_h[root * OrganicMatterFireExchange.salt_species_count + salt];
                };
            }
            var total_symbiont_structural_carbon_g_c: f64 = 0;
            var total_symbiont_mobile_carbon_g_c: f64 = 0;
            var total_root_mobile_carbon_g_c: f64 = 0;
            var total_root_structural_carbon_g_c: f64 = 0;
            for (0..species_count) |species| {
                const plant = cell * species_count + species;
                const root = try roots.layerIndex(plant, 0, layer);
                const symbiont_structural = roots.symbiont_structural_carbon_g_c[root];
                const symbiont_mobile = roots.symbiont_mobile_carbon_g_c[root];
                if (!std.math.isFinite(symbiont_structural) or symbiont_structural < 0 or !std.math.isFinite(symbiont_mobile) or symbiont_mobile < 0) return error.InvalidRootSymbiontPool;
                total_symbiont_structural_carbon_g_c += symbiont_structural;
                total_symbiont_mobile_carbon_g_c += symbiont_mobile;
                const biological_domain_count = try domainCount(biological_domain_count_by_plant, plant);
                for (0..biological_domain_count) |domain| {
                    const domain_root = try roots.layerIndex(plant, domain, layer);
                    const mobile = roots.mobile_carbon_g[domain_root];
                    if (!std.math.isFinite(mobile) or mobile < 0) return error.InvalidRootCombustionPool;
                    total_root_mobile_carbon_g_c += mobile;
                    for (0..roots.root_axis_count) |axis| {
                        const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
                        const structural = roots.axis_primary_carbon_g[axis_layer] + roots.axis_secondary_carbon_g[axis_layer];
                        if (!std.math.isFinite(structural) or structural < 0) return error.InvalidRootCombustionPool;
                        total_root_structural_carbon_g_c += structural;
                    }
                }
            }
            const symbiont_structural_fraction = try RootDisturbance.combustionFraction(total_symbiont_structural_carbon_g_c, grid.soil_temperature_k[soil], cell_area_m2[cell], timestep_h, parameters.nonwoody_structural_specific_combustion_g_c_per_m2_h, parameters);
            const symbiont_mobile_fraction = try RootDisturbance.combustionFraction(total_symbiont_mobile_carbon_g_c, grid.soil_temperature_k[soil], cell_area_m2[cell], timestep_h, parameters.mobile_and_leaf_specific_combustion_g_c_per_m2_h, parameters);
            const root_mobile_fraction = try RootDisturbance.combustionFraction(total_root_mobile_carbon_g_c, grid.soil_temperature_k[soil], cell_area_m2[cell], timestep_h, parameters.mobile_and_leaf_specific_combustion_g_c_per_m2_h, parameters);
            const root_structural_fraction = try RootDisturbance.combustionFraction(total_root_structural_carbon_g_c, grid.soil_temperature_k[soil], cell_area_m2[cell], timestep_h, parameters.root_structural_specific_combustion_g_c_per_m2_h, parameters);

            // Validate the complete layer before publishing any mutation.
            for (0..species_count) |species| {
                const plant = cell * species_count + species;
                const root = try roots.layerIndex(plant, 0, layer);
                _ = try combustionResult(roots, root, symbiont_structural_fraction, symbiont_mobile_fraction);
                try validateRootCombustion(roots, plant, layer, try domainCount(biological_domain_count_by_plant, plant), root_mobile_fraction, root_structural_fraction, dynamic_salts);
                if (layer == 0) if (canopy) |state| try validateStorageCombustion(state, plant);
            }
            for (0..species_count) |species| {
                const plant = cell * species_count + species;
                const root = try roots.layerIndex(plant, 0, layer);
                const result = combustionResult(roots, root, symbiont_structural_fraction, symbiont_mobile_fraction) catch unreachable;
                roots.symbiont_structural_carbon_g_c[root] = result.structural.carbon_g_c;
                roots.symbiont_structural_nitrogen_g_n[root] = result.structural.nitrogen_g_n;
                roots.symbiont_structural_phosphorus_g_p[root] = result.structural.phosphorus_g_p;
                roots.symbiont_mobile_carbon_g_c[root] = result.mobile.carbon_g_c;
                roots.symbiont_mobile_nitrogen_g_n[root] = result.mobile.nitrogen_g_n;
                roots.symbiont_mobile_phosphorus_g_p[root] = result.mobile.phosphorus_g_p;
                roots.combustion_carbon_loss_g_c_per_h[plant] -= result.emitted.carbon_g_c;
                roots.combustion_nitrogen_loss_g_n_per_h[plant] -= result.emitted.nitrogen_g_n;
                roots.combustion_phosphorus_loss_g_p_per_h[plant] -= result.emitted.phosphorus_g_p;
                roots.symbiont_combustion_g_c_per_h[root] += result.emitted.carbon_g_c;
                commitRootCombustion(roots, plant, layer, biological_domain_count_by_plant[plant], root_mobile_fraction, root_structural_fraction, dynamic_salts);
                if (layer == 0) if (canopy) |state| commitStorageCombustion(state, roots, plant, root, root_structural_fraction);
            }
            var carbon_loss_g_c: f64 = 0;
            var nitrogen_loss_g_n: f64 = 0;
            var phosphorus_loss_g_p: f64 = 0;
            var salt_loss_mol: [OrganicMatterFireExchange.salt_species_count]f64 = @splat(0);
            for (0..species_count) |species| {
                const plant = cell * species_count + species;
                carbon_loss_g_c += roots.combustion_carbon_loss_g_c_per_h[plant];
                nitrogen_loss_g_n += roots.combustion_nitrogen_loss_g_n_per_h[plant];
                phosphorus_loss_g_p += roots.combustion_phosphorus_loss_g_p_per_h[plant];
                const biological_domain_count = try domainCount(biological_domain_count_by_plant, plant);
                if (dynamic_salts) for (0..biological_domain_count) |domain| {
                    const root = try roots.layerIndex(plant, domain, layer);
                    for (0..OrganicMatterFireExchange.salt_species_count) |salt| salt_loss_mol[salt] += roots.combustion_salt_loss_mol_per_h[root * OrganicMatterFireExchange.salt_species_count + salt];
                };
            }
            for (0..OrganicMatterFireExchange.salt_species_count) |salt| salt_loss_mol[salt] -= previous_salt_loss_mol[salt];
            try fire_exchange.addCombustedPoolsForSubstrate(
                soil,
                1,
                previous_carbon_loss_g_c - carbon_loss_g_c,
                previous_nitrogen_loss_g_n - nitrogen_loss_g_n,
                previous_phosphorus_loss_g_p - phosphorus_loss_g_p,
                &salt_loss_mol,
            );
        }
    }
}

fn validateStorageCombustion(canopy: *const Canopy.State, plant: usize) !void {
    inline for (.{
        canopy.plant_seed_storage_carbon_g[plant],
        canopy.plant_seed_storage_nitrogen_g[plant],
        canopy.plant_seed_storage_phosphorus_g[plant],
    }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootCombustionPool;
}

fn commitStorageCombustion(canopy: *Canopy.State, roots: *RootSystem.State, plant: usize, top_root: usize, fraction: f64) void {
    const carbon = canopy.plant_seed_storage_carbon_g[plant] * fraction;
    const nitrogen = canopy.plant_seed_storage_nitrogen_g[plant] * fraction;
    const phosphorus = canopy.plant_seed_storage_phosphorus_g[plant] * fraction;
    canopy.plant_seed_storage_carbon_g[plant] -= carbon;
    canopy.plant_seed_storage_nitrogen_g[plant] -= nitrogen;
    canopy.plant_seed_storage_phosphorus_g[plant] -= phosphorus;
    roots.combustion_carbon_loss_g_c_per_h[plant] -= carbon;
    roots.combustion_nitrogen_loss_g_n_per_h[plant] -= nitrogen;
    roots.combustion_phosphorus_loss_g_p_per_h[plant] -= phosphorus;
    roots.root_combustion_g_c_per_h[top_root] += carbon;
}

fn domainCount(counts: []const u8, plant: usize) !u8 {
    if (plant >= counts.len or counts[plant] < 1 or counts[plant] > RootSystem.biological_domain_count)
        return error.DisturbancePlantDimensionMismatch;
    return counts[plant];
}

fn validateRootCombustion(roots: *const RootSystem.State, plant: usize, layer: usize, biological_domain_count: u8, mobile_fraction: f64, structural_fraction: f64, dynamic_salts: bool) !void {
    for (0..biological_domain_count) |domain| {
        const root = try roots.layerIndex(plant, domain, layer);
        inline for (.{ roots.mobile_carbon_g[root], roots.mobile_nitrogen_g[root], roots.mobile_phosphorus_g[root] }) |value|
            if (!std.math.isFinite(value) or value < 0) return error.InvalidRootCombustionPool;
        for (0..roots.root_axis_count) |axis| {
            const axis_layer = try roots.layerAxisIndex(plant, domain, layer, axis);
            inline for (.{
                roots.axis_primary_carbon_g[axis_layer],   roots.axis_primary_nitrogen_g[axis_layer],   roots.axis_primary_phosphorus_g[axis_layer],
                roots.axis_secondary_carbon_g[axis_layer], roots.axis_secondary_nitrogen_g[axis_layer], roots.axis_secondary_phosphorus_g[axis_layer],
                roots.axis_primary_length_m[axis_layer],   roots.axis_secondary_length_m[axis_layer],   roots.axis_primary_count[axis_layer],
                roots.axis_secondary_count[axis_layer],
            }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidRootCombustionPool;
        }
        if (dynamic_salts) for (0..RootSystem.salt_species_count) |salt| {
            const index = root * RootSystem.salt_species_count + salt;
            if (!std.math.isFinite(roots.salt_content_mol[index]) or roots.salt_content_mol[index] < 0) return error.InvalidRootCombustionPool;
        };
    }
    inline for (.{ mobile_fraction, structural_fraction }) |fraction|
        if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidRootSymbiontCombustionFraction;
}

fn commitRootCombustion(roots: *RootSystem.State, plant: usize, layer: usize, biological_domain_count: u8, mobile_fraction: f64, structural_fraction: f64, dynamic_salts: bool) void {
    var emitted_c: f64 = 0;
    var emitted_n: f64 = 0;
    var emitted_p: f64 = 0;
    for (0..biological_domain_count) |domain| {
        const root = roots.layerIndex(plant, domain, layer) catch unreachable;
        var domain_emitted_c: f64 = 0;
        inline for (.{ "carbon", "nitrogen", "phosphorus" }) |element| {
            const field_name = "mobile_" ++ element ++ "_g";
            const burned = @field(roots, field_name)[root] * mobile_fraction;
            @field(roots, field_name)[root] -= burned;
            if (comptime std.mem.eql(u8, element, "carbon")) {
                emitted_c += burned;
                domain_emitted_c += burned;
            } else if (comptime std.mem.eql(u8, element, "nitrogen")) emitted_n += burned else emitted_p += burned;
        }
        for (0..roots.root_axis_count) |axis| {
            const axis_layer = roots.layerAxisIndex(plant, domain, layer, axis) catch unreachable;
            inline for (.{ "primary", "secondary" }) |order| {
                inline for (.{ "carbon", "nitrogen", "phosphorus" }) |element| {
                    const field_name = "axis_" ++ order ++ "_" ++ element ++ "_g";
                    const burned = @field(roots, field_name)[axis_layer] * structural_fraction;
                    @field(roots, field_name)[axis_layer] -= burned;
                    if (comptime std.mem.eql(u8, element, "carbon")) {
                        emitted_c += burned;
                        domain_emitted_c += burned;
                    } else if (comptime std.mem.eql(u8, element, "nitrogen")) emitted_n += burned else emitted_p += burned;
                }
                @field(roots, "axis_" ++ order ++ "_length_m")[axis_layer] *= 1 - structural_fraction;
                @field(roots, "axis_" ++ order ++ "_count")[axis_layer] *= 1 - structural_fraction;
            }
        }
        if (dynamic_salts) for (0..RootSystem.salt_species_count) |salt| {
            const index = root * RootSystem.salt_species_count + salt;
            const burned = roots.salt_content_mol[index] * mobile_fraction;
            roots.salt_content_mol[index] -= burned;
            roots.combustion_salt_loss_mol_per_h[index] += burned;
        };
        roots.root_combustion_g_c_per_h[root] += domain_emitted_c;
    }
    roots.combustion_carbon_loss_g_c_per_h[plant] -= emitted_c;
    roots.combustion_nitrogen_loss_g_n_per_h[plant] -= emitted_n;
    roots.combustion_phosphorus_loss_g_p_per_h[plant] -= emitted_p;
}

fn combustionResult(roots: *const RootSystem.State, root: usize, structural_fraction: f64, mobile_fraction: f64) !RootDisturbance.CombustionResult {
    return RootDisturbance.combustSymbiont(
        .{ .carbon_g_c = roots.symbiont_structural_carbon_g_c[root], .nitrogen_g_n = roots.symbiont_structural_nitrogen_g_n[root], .phosphorus_g_p = roots.symbiont_structural_phosphorus_g_p[root] },
        .{ .carbon_g_c = roots.symbiont_mobile_carbon_g_c[root], .nitrogen_g_n = roots.symbiont_mobile_nitrogen_g_n[root], .phosphorus_g_p = roots.symbiont_mobile_phosphorus_g_p[root] },
        structural_fraction,
        mobile_fraction,
    );
}

fn calculateLayer(roots: *const RootSystem.State, root: usize, retention: RootDisturbance.ElementRetention, fine: LitterPartition.ElementFractions, mobile: LitterPartition.ElementFractions) !RootDisturbance.Result {
    return RootDisturbance.retainAndRelease(
        .{ .carbon_g_c = roots.symbiont_structural_carbon_g_c[root], .nitrogen_g_n = roots.symbiont_structural_nitrogen_g_n[root], .phosphorus_g_p = roots.symbiont_structural_phosphorus_g_p[root] },
        .{ .carbon_g_c = roots.symbiont_mobile_carbon_g_c[root], .nitrogen_g_n = roots.symbiont_mobile_nitrogen_g_n[root], .phosphorus_g_p = roots.symbiont_mobile_phosphorus_g_p[root] },
        retention,
        fine,
        mobile,
    );
}

test "runtime disturbance map resolves schedules and case-insensitive no" {
    var assignments = try land_management.parse(std.testing.allocator, "1 1 1 1\nannual NO NO\n2 1 2 1\nnO NO NO\n");
    defer assignments.deinit();
    const unit_map = try assignments.buildCellUnitMap(std.testing.allocator, .{
        .west_column = 1,
        .north_row = 1,
        .east_column = 2,
        .south_row = 1,
    });
    defer std.testing.allocator.free(unit_map);
    var catalog = disturbance.Catalog.init(std.testing.allocator);
    defer catalog.deinit();
    _ = try catalog.appendFromSource("annual", "01010000,10,0.15\n");
    var map = try ScheduleMap.init(std.testing.allocator, assignments, unit_map, catalog);
    defer map.deinit();
    try std.testing.expectEqual(@as(?usize, 0), map.catalog_index_by_cell[0]);
    try std.testing.expectEqual(@as(?usize, null), map.catalog_index_by_cell[1]);
}

test "disturbance dispatch date preserves DAY modulo-four chronology" {
    try validateDispatchDate(.{ .day = 29, .month = 2, .year = 1900 });
    try validateDispatchDate(.{ .day = 1, .month = 3, .year = 1900 });
    try validateDispatchDate(.{ .day = 30, .month = 3, .year = 1900 });
    try std.testing.expectError(error.InvalidDisturbanceDispatchDate, validateDispatchDate(.{ .day = 0, .month = 1, .year = 1900 }));
    try std.testing.expectError(
        error.InvalidDisturbanceDispatchDate,
        validateDispatchDate(.{ .day = 1, .month = 1, .year = 0 }),
    );
}

test "GROSUB fire combustion conserves runtime species nodule C N P" {
    const SimulationConfig = @import("config.zig").SimulationConfig;
    const cfg = try SimulationConfig.init(
        .{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 7 },
        .{ .worker_threads = 1, .tile_cells = 1 },
        .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-12, .max_nonlinear_iterations = 10 },
    );
    var grid = try Grid.init(std.testing.allocator, cfg);
    defer grid.deinit();
    grid.soil_temperature_k[0] = 500;
    var roots = try RootSystem.State.init(std.testing.allocator, 7, 1, 1);
    defer roots.deinit();
    var canopy = try Canopy.State.init(std.testing.allocator, 1, 7, &([_]usize{1} ** 7), &([_]usize{1} ** 7), &([_]usize{1} ** 7));
    defer canopy.deinit();
    var initial_c: f64 = 0;
    var initial_n: f64 = 0;
    var initial_p: f64 = 0;
    var initial_salt_mol: f64 = 0;
    for (0..7) |plant| {
        const root = try roots.layerIndex(plant, 0, 0);
        const scale: f64 = @floatFromInt(plant + 1);
        roots.symbiont_structural_carbon_g_c[root] = 2 * scale;
        roots.symbiont_structural_nitrogen_g_n[root] = 0.2 * scale;
        roots.symbiont_structural_phosphorus_g_p[root] = 0.02 * scale;
        roots.symbiont_mobile_carbon_g_c[root] = scale;
        roots.symbiont_mobile_nitrogen_g_n[root] = 0.1 * scale;
        roots.symbiont_mobile_phosphorus_g_p[root] = 0.01 * scale;
        initial_c += 3 * scale;
        initial_n += 0.3 * scale;
        initial_p += 0.03 * scale;
        canopy.plant_seed_storage_carbon_g[plant] = 0.4 * scale;
        canopy.plant_seed_storage_nitrogen_g[plant] = 0.04 * scale;
        canopy.plant_seed_storage_phosphorus_g[plant] = 0.004 * scale;
        initial_c += 0.4 * scale;
        initial_n += 0.04 * scale;
        initial_p += 0.004 * scale;
        for (0..RootSystem.biological_domain_count) |domain| {
            const domain_root = try roots.layerIndex(plant, domain, 0);
            const axis_layer = try roots.layerAxisIndex(plant, domain, 0, 0);
            roots.mobile_carbon_g[domain_root] = 0.5 * scale;
            roots.mobile_nitrogen_g[domain_root] = 0.05 * scale;
            roots.mobile_phosphorus_g[domain_root] = 0.005 * scale;
            roots.axis_primary_carbon_g[axis_layer] = scale;
            roots.axis_primary_nitrogen_g[axis_layer] = 0.1 * scale;
            roots.axis_primary_phosphorus_g[axis_layer] = 0.01 * scale;
            roots.axis_secondary_carbon_g[axis_layer] = 0.25 * scale;
            roots.axis_secondary_nitrogen_g[axis_layer] = 0.025 * scale;
            roots.axis_secondary_phosphorus_g[axis_layer] = 0.0025 * scale;
            roots.axis_primary_length_m[axis_layer] = scale;
            roots.axis_secondary_length_m[axis_layer] = 2 * scale;
            roots.axis_primary_count[axis_layer] = scale;
            roots.axis_secondary_count[axis_layer] = 2 * scale;
            initial_c += 1.75 * scale;
            initial_n += 0.175 * scale;
            initial_p += 0.0175 * scale;
            for (0..RootSystem.salt_species_count) |salt| {
                const amount = 0.001 * scale * @as(f64, @floatFromInt(salt + 1));
                roots.salt_content_mol[domain_root * RootSystem.salt_species_count + salt] = amount;
                initial_salt_mol += amount;
            }
        }
    }
    var fire_exchange = try OrganicMatterFireExchange.State.init(std.testing.allocator, grid.layer_count, SoilOrganic.microbial_substrate_count);
    defer fire_exchange.deinit();
    var biological_domain_counts = [_]u8{2} ** 7;
    biological_domain_counts[0] = 1;
    try applyRootFireCombustion(&roots, &canopy, &grid, 7, &biological_domain_counts, &.{0.01}, &.{true}, 1, true, RootDisturbance.sourceCombustionParameters(), &fire_exchange);
    var remaining_c: f64 = 0;
    var remaining_n: f64 = 0;
    var remaining_p: f64 = 0;
    var emitted_c: f64 = 0;
    var emitted_n: f64 = 0;
    var emitted_p: f64 = 0;
    var remaining_salt_mol: f64 = 0;
    var emitted_salt_mol: f64 = 0;
    for (0..7) |plant| {
        const root = try roots.layerIndex(plant, 0, 0);
        remaining_c += roots.symbiont_structural_carbon_g_c[root] + roots.symbiont_mobile_carbon_g_c[root];
        remaining_n += roots.symbiont_structural_nitrogen_g_n[root] + roots.symbiont_mobile_nitrogen_g_n[root];
        remaining_p += roots.symbiont_structural_phosphorus_g_p[root] + roots.symbiont_mobile_phosphorus_g_p[root];
        remaining_c += canopy.plant_seed_storage_carbon_g[plant];
        remaining_n += canopy.plant_seed_storage_nitrogen_g[plant];
        remaining_p += canopy.plant_seed_storage_phosphorus_g[plant];
        for (0..RootSystem.biological_domain_count) |domain| {
            const domain_root = try roots.layerIndex(plant, domain, 0);
            const axis_layer = try roots.layerAxisIndex(plant, domain, 0, 0);
            remaining_c += roots.mobile_carbon_g[domain_root] + roots.axis_primary_carbon_g[axis_layer] + roots.axis_secondary_carbon_g[axis_layer];
            remaining_n += roots.mobile_nitrogen_g[domain_root] + roots.axis_primary_nitrogen_g[axis_layer] + roots.axis_secondary_nitrogen_g[axis_layer];
            remaining_p += roots.mobile_phosphorus_g[domain_root] + roots.axis_primary_phosphorus_g[axis_layer] + roots.axis_secondary_phosphorus_g[axis_layer];
            if (plant == 0 and domain == 1)
                try std.testing.expectEqual(@as(f64, 1), roots.axis_primary_length_m[axis_layer])
            else
                try std.testing.expect(roots.axis_primary_length_m[axis_layer] < @as(f64, @floatFromInt(plant + 1)));
            for (0..RootSystem.salt_species_count) |salt| {
                const salt_index = domain_root * RootSystem.salt_species_count + salt;
                remaining_salt_mol += roots.salt_content_mol[salt_index];
                emitted_salt_mol += roots.combustion_salt_loss_mol_per_h[salt_index];
            }
        }
        emitted_c -= roots.combustion_carbon_loss_g_c_per_h[plant];
        emitted_n -= roots.combustion_nitrogen_loss_g_n_per_h[plant];
        emitted_p -= roots.combustion_phosphorus_loss_g_p_per_h[plant];
        var root_emitted_c: f64 = 0;
        for (0..RootSystem.biological_domain_count) |domain| root_emitted_c += roots.root_combustion_g_c_per_h[try roots.layerIndex(plant, domain, 0)];
        try std.testing.expectApproxEqAbs(-roots.combustion_carbon_loss_g_c_per_h[plant], roots.symbiont_combustion_g_c_per_h[root] + root_emitted_c, 1e-12);
    }
    try std.testing.expectApproxEqAbs(initial_c, remaining_c + emitted_c, 1e-12);
    try std.testing.expectApproxEqAbs(initial_n, remaining_n + emitted_n, 1e-12);
    try std.testing.expectApproxEqAbs(initial_p, remaining_p + emitted_p, 1e-12);
    try std.testing.expectApproxEqAbs(initial_salt_mol, remaining_salt_mol + emitted_salt_mol, 1e-12);
    try std.testing.expectApproxEqAbs(emitted_c, fire_exchange.unlimited_combustion_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(emitted_n, fire_exchange.combusted_nitrogen_g_n[0], 1e-12);
    try std.testing.expectApproxEqAbs(emitted_p, fire_exchange.combusted_phosphorus_g_p[0], 1e-12);
    var ledger_salt_mol: f64 = 0;
    for (0..OrganicMatterFireExchange.salt_species_count) |salt| ledger_salt_mol += fire_exchange.released_salt_mol[salt];
    try std.testing.expectApproxEqAbs(emitted_salt_mol, ledger_salt_mol, 1e-12);
}

test "GROSUB tillage removes runtime host and mycorrhizal roots conservatively" {
    var roots = try RootSystem.State.init(std.testing.allocator, 1, 1, 2);
    defer roots.deinit();
    roots.active_root_axis_count[0] = 2;
    var structural = [3]f64{ 0, 0, 0 };
    var mobile_pool = [3]f64{ 0, 0, 0 };
    for (0..RootSystem.biological_domain_count) |domain| {
        const root = try roots.layerIndex(0, domain, 0);
        roots.mobile_carbon_g[root] = 4 + @as(f64, @floatFromInt(domain));
        roots.mobile_nitrogen_g[root] = 2 + @as(f64, @floatFromInt(domain));
        roots.mobile_phosphorus_g[root] = 1 + @as(f64, @floatFromInt(domain));
        roots.protein_carbon_g[root] = 3;
        mobile_pool[0] += roots.mobile_carbon_g[root];
        mobile_pool[1] += roots.mobile_nitrogen_g[root];
        mobile_pool[2] += roots.mobile_phosphorus_g[root];
        for (0..2) |axis| {
            const index = try roots.layerAxisIndex(0, domain, 0, axis);
            const scale = @as(f64, @floatFromInt(1 + domain + axis));
            roots.axis_primary_carbon_g[index] = 2 * scale;
            roots.axis_secondary_carbon_g[index] = scale;
            roots.axis_primary_nitrogen_g[index] = scale;
            roots.axis_secondary_nitrogen_g[index] = 0.5 * scale;
            roots.axis_primary_phosphorus_g[index] = 0.2 * scale;
            roots.axis_secondary_phosphorus_g[index] = 0.1 * scale;
            roots.axis_primary_length_m[index] = 10 * scale;
            roots.axis_secondary_length_m[index] = 5 * scale;
            roots.axis_primary_count[index] = scale;
            roots.axis_secondary_count[index] = 2 * scale;
            structural[0] += 3 * scale;
            structural[1] += 1.5 * scale;
            structural[2] += 0.3 * scale;
        }
    }
    const first = LitterPartition.ElementFractions{
        .carbon = .{ 1, 0, 0, 0 },
        .nitrogen = .{ 1, 0, 0, 0 },
        .phosphorus = .{ 1, 0, 0, 0 },
    };
    const second = LitterPartition.ElementFractions{
        .carbon = .{ 0, 1, 0, 0 },
        .nitrogen = .{ 0, 1, 0, 0 },
        .phosphorus = .{ 0, 1, 0, 0 },
    };
    const third = LitterPartition.ElementFractions{
        .carbon = .{ 0, 0, 1, 0 },
        .nitrogen = .{ 0, 0, 1, 0 },
        .phosphorus = .{ 0, 0, 1, 0 },
    };
    const removed_fraction = 0.25;
    const nonwoody_fraction = 0.4;
    const litter = try calculateHostRootTillage(
        &roots,
        0,
        0,
        RootSystem.biological_domain_count,
        removed_fraction,
        nonwoody_fraction,
        first,
        second,
        third,
    );
    const domain_zero_litter = try calculateHostRootTillageDomain(&roots, 0, 0, 0, removed_fraction, nonwoody_fraction, first, second, third);
    const domain_one_litter = try calculateHostRootTillageDomain(&roots, 0, 1, 0, removed_fraction, nonwoody_fraction, first, second, third);
    try std.testing.expectApproxEqAbs(
        try RootLitterLedger.totalCarbon(litter),
        try RootLitterLedger.totalCarbon(domain_zero_litter) +
            try RootLitterLedger.totalCarbon(domain_one_litter),
        1e-12,
    );
    try std.testing.expectApproxEqAbs(removed_fraction * structural[0] * (1 - nonwoody_fraction), litter.woody_carbon_g_c[0], 1e-12);
    try std.testing.expectApproxEqAbs(removed_fraction * structural[0] * nonwoody_fraction, litter.nonwoody_carbon_g_c[1], 1e-12);
    try std.testing.expectApproxEqAbs(removed_fraction * mobile_pool[0], litter.nonwoody_carbon_g_c[2], 1e-12);
    try std.testing.expectApproxEqAbs(removed_fraction * (structural[1] + mobile_pool[1]), sumElement(litter.woody_nitrogen_g_n) + sumElement(litter.nonwoody_nitrogen_g_n), 1e-12);
    try std.testing.expectApproxEqAbs(removed_fraction * (structural[2] + mobile_pool[2]), sumElement(litter.woody_phosphorus_g_p) + sumElement(litter.nonwoody_phosphorus_g_p), 1e-12);

    try commitHostRootTillage(&roots, 0, 0, RootSystem.biological_domain_count, removed_fraction);
    const retained = 1 - removed_fraction;
    const first_root = try roots.layerIndex(0, 0, 0);
    const first_axis = try roots.layerAxisIndex(0, 0, 0, 0);
    try std.testing.expectApproxEqAbs(@as(f64, 4) * retained, roots.mobile_carbon_g[first_root], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2) * retained, roots.axis_primary_carbon_g[first_axis], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 10) * retained, roots.axis_primary_length_m[first_axis], 1e-12);
    try std.testing.expectApproxEqAbs(retained, roots.axis_primary_count[first_axis], 1e-12);
}

fn sumElement(values: [LitterPartition.kinetic_component_count]f64) f64 {
    var total: f64 = 0;
    for (values) |value| total += value;
    return total;
}
