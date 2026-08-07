const std = @import("std");
const Canopy = @import("../../canopy/photosynthesis/photosynthesis.zig");
const CanopyControls = @import("../../canopy/radiation/layer_distribution.zig").Controls;
const Roots = @import("../root/plant_root_system.zig");
const FireScience = @import("../root/plant_root_disturbance.zig");
const PlantSoilExchange = @import("../exchange/soil.zig");
const SoilOrganic = @import("../../soil/organic/initialization.zig");

const Totals = struct {
    mobile_c: f64 = 0,
    leaf_c: f64 = 0,
    sheath_c: f64 = 0,
    stalk_c: f64 = 0,
    husk_c: f64 = 0,
    ear_c: f64 = 0,
    grain_c: f64 = 0,
    symbiont_mobile_c: f64 = 0,
    symbiont_structural_c: f64 = 0,
    standing_dead_c: f64 = 0,
};

const Fractions = struct {
    mobile: f64,
    leaf: f64,
    sheath: f64,
    stalk: f64,
    husk: f64,
    ear: f64,
    grain: f64,
    symbiont_mobile: f64,
    symbiont_structural: f64,
    standing_dead: f64,
};

/// GROSUB live-shoot and standing-dead combustion. Pool denominators are
/// cell-wide across all runtime species while temperature responses remain
/// plant-specific, exactly as in the source PFT loop.
pub fn apply(
    canopy: *Canopy.State,
    controls: *const CanopyControls,
    roots: *Roots.State,
    cell_area_m2: []const f64,
    fire_active_this_hour: []const bool,
    timestep_h: f64,
    dynamic_salts: bool,
    parameters: FireScience.CombustionParameters,
    canopy_air_temperature_k: []const f64,
    canopy_air_volume_m3: []const f64,
    oxygen_concentration_umol_per_mol: f64,
    methane_concentration_umol_per_mol: f64,
    oxygen_concentration_g_per_m3: f64,
    delayed_live_combustion_heat_megajoules: []f64,
    delayed_standing_dead_combustion_heat_megajoules: []f64,
    surface_organic: *SoilOrganic.State,
) !void {
    try parameters.validate();
    const plant_count = canopy.cell_count * canopy.species_count;
    if (controls.biomass_turnover_type.len != plant_count or controls.root_profile_type.len != plant_count or roots.plant_count != plant_count) return error.ShootFirePlantDimensionMismatch;
    if (cell_area_m2.len != canopy.cell_count or fire_active_this_hour.len != canopy.cell_count or canopy_air_temperature_k.len != canopy.cell_count or canopy_air_volume_m3.len != canopy.cell_count) return error.ShootFireCellDimensionMismatch;
    if (delayed_live_combustion_heat_megajoules.len != plant_count or delayed_standing_dead_combustion_heat_megajoules.len != plant_count) return error.ShootFirePlantDimensionMismatch;
    if (surface_organic.layer_count != canopy.cell_count) return error.ShootFireSurfaceOrganicDimensionMismatch;
    if (!std.math.isFinite(timestep_h) or timestep_h <= 0) return error.InvalidShootFireInput;
    inline for (.{ oxygen_concentration_umol_per_mol, methane_concentration_umol_per_mol, oxygen_concentration_g_per_m3 }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidShootFireInput;
    @memset(canopy.plant_combustion_carbon_loss_g_c_per_h, 0);
    @memset(canopy.plant_combustion_nitrogen_loss_g_n_per_h, 0);
    @memset(canopy.plant_combustion_phosphorus_loss_g_p_per_h, 0);
    @memset(canopy.plant_live_combustion_g_c_per_h, 0);
    @memset(canopy.plant_standing_dead_combustion_g_c_per_h, 0);
    @memset(canopy.branch_combustion_salt_loss_by_species_mol_per_h, 0);
    @memset(canopy.plant_fire_carbon_dioxide_emission_g_c_per_h, 0);
    @memset(canopy.plant_fire_methane_emission_g_c_per_h, 0);
    @memset(canopy.plant_fire_oxygen_consumption_g_o_per_h, 0);
    @memset(canopy.plant_fire_charcoal_production_g_c_per_h, 0);
    @memset(canopy.plant_fire_heat_release_megajoules_per_h, 0);
    for (0..canopy.cell_count) |cell| {
        if (!fire_active_this_hour[cell]) continue;
        const totals = try cellTotals(canopy, cell);
        for (0..canopy.species_count) |species| {
            const plant = try canopy.plantIndex(cell, species);
            const woody = controls.biomass_turnover_type[plant] != 0 and controls.root_profile_type[plant] > 1;
            const fractions = try calculateFractions(canopy, plant, totals, cell_area_m2[cell], timestep_h, woody, parameters);
            try validatePlant(canopy, plant, fractions, dynamic_salts);
        }
        for (0..canopy.species_count) |species| {
            const plant = canopy.plantIndex(cell, species) catch unreachable;
            const woody = controls.biomass_turnover_type[plant] != 0 and controls.root_profile_type[plant] > 1;
            const fractions = calculateFractions(canopy, plant, totals, cell_area_m2[cell], timestep_h, woody, parameters) catch unreachable;
            commitPlant(canopy, roots, plant, fractions, dynamic_salts);
        }
        try partitionCellProducts(
            canopy,
            cell,
            canopy_air_temperature_k[cell],
            canopy_air_volume_m3[cell] * oxygen_concentration_g_per_m3,
            oxygen_concentration_umol_per_mol,
            methane_concentration_umol_per_mol,
            parameters,
            delayed_live_combustion_heat_megajoules,
            delayed_standing_dead_combustion_heat_megajoules,
            surface_organic,
        );
    }
}

fn partitionCellProducts(canopy: *Canopy.State, cell: usize, temperature_k: f64, oxygen_content_g: f64, oxygen_umol_per_mol: f64, methane_umol_per_mol: f64, parameters: FireScience.CombustionParameters, delayed_live_heat_megajoules: []f64, delayed_dead_heat_megajoules: []f64, surface_organic: *SoilOrganic.State) !void {
    var total_combusted_c: f64 = 0;
    for (0..canopy.species_count) |species| total_combusted_c -= canopy.plant_combustion_carbon_loss_g_c_per_h[try canopy.plantIndex(cell, species)];
    if (total_combusted_c == 0) return;
    const products = try PlantSoilExchange.canopyFireCombustion(
        total_combusted_c,
        0,
        temperature_k,
        oxygen_umol_per_mol,
        oxygen_content_g,
        methane_umol_per_mol,
        .{
            .gas_constant_j_per_mol_k = parameters.gas_constant_j_per_mol_k,
            .combustion_activation_energy_j_per_mol = parameters.activation_energy_j_per_mol,
            .combustion_temperature_intercept = parameters.arrhenius_intercept,
            .oxygen_g_per_g_carbon = parameters.oxygen_g_per_g_combusted_carbon,
            .maximum_aerobic_charcoal_fraction = parameters.maximum_aerobic_charcoal_fraction,
            .maximum_anaerobic_charcoal_fraction = parameters.maximum_anaerobic_charcoal_fraction,
            .oxygen_half_saturation_umol_per_mol = parameters.oxygen_half_saturation_umol_per_mol,
            .methane_half_saturation_umol_per_mol = parameters.methane_half_saturation_umol_per_mol,
            .aerobic_combustion_energy_megajoules_per_g_carbon = parameters.aerobic_combustion_energy_megajoules_per_g_c,
            .anaerobic_combustion_energy_megajoules_per_g_carbon = parameters.anaerobic_combustion_energy_megajoules_per_g_c,
            .methane_combustion_energy_megajoules_per_g_carbon = parameters.methane_combustion_energy_megajoules_per_g_c,
        },
    );
    for (0..canopy.species_count) |species| {
        const plant = try canopy.plantIndex(cell, species);
        const share = -canopy.plant_combustion_carbon_loss_g_c_per_h[plant] / total_combusted_c;
        canopy.plant_fire_carbon_dioxide_emission_g_c_per_h[plant] = share * products.carbon_dioxide_emitted_g_carbon;
        canopy.plant_fire_methane_emission_g_c_per_h[plant] = share * products.methane_emitted_g_carbon;
        canopy.plant_fire_oxygen_consumption_g_o_per_h[plant] = share * products.oxygen_consumed_g;
        canopy.plant_fire_charcoal_production_g_c_per_h[plant] = share * products.charcoal_produced_g_carbon;
        canopy.plant_fire_heat_release_megajoules_per_h[plant] = share * products.heat_released_megajoules;
        const plant_combusted_c = canopy.plant_live_combustion_g_c_per_h[plant] + canopy.plant_standing_dead_combustion_g_c_per_h[plant];
        const live_fraction = if (plant_combusted_c > 0) canopy.plant_live_combustion_g_c_per_h[plant] / plant_combusted_c else 0;
        const live_heat = canopy.plant_fire_heat_release_megajoules_per_h[plant] * live_fraction;
        const dead_heat = canopy.plant_fire_heat_release_megajoules_per_h[plant] - live_heat;
        const next_live_heat = delayed_live_heat_megajoules[plant] + live_heat;
        const next_dead_heat = delayed_dead_heat_megajoules[plant] + dead_heat;
        if (!std.math.isFinite(next_live_heat) or !std.math.isFinite(next_dead_heat)) return error.NonFiniteShootFireHeat;
        delayed_live_heat_megajoules[plant] = next_live_heat;
        delayed_dead_heat_megajoules[plant] = next_dead_heat;
    }
    // EXTRACT routes canopy-fire charcoal to OSC(5,1,0): the most
    // resistant structural fraction of the first surface substrate.
    const charcoal_index = cell * SoilOrganic.substrate_count * SoilOrganic.structural_fraction_count +
        (SoilOrganic.structural_fraction_count - 1);
    const next_charcoal = surface_organic.structural[charcoal_index].carbon_g_c + products.charcoal_produced_g_carbon;
    if (!std.math.isFinite(next_charcoal) or next_charcoal < 0) return error.NonFiniteSurfaceCharcoal;
    surface_organic.structural[charcoal_index].carbon_g_c = next_charcoal;
}

fn cellTotals(canopy: *const Canopy.State, cell: usize) !Totals {
    var totals: Totals = .{};
    for (0..canopy.species_count) |species| {
        const plant = try canopy.plantIndex(cell, species);
        const branches = try canopy.branchRange(plant);
        for (branches.first..branches.end) |branch| {
            totals.mobile_c += try pool(canopy.branch_mobile_carbon_g[branch]);
            totals.leaf_c += try pool(canopy.branch_leaf_carbon_g[branch]);
            totals.sheath_c += try pool(canopy.branch_sheath_carbon_g[branch]);
            totals.stalk_c += try pool(canopy.branch_stalk_carbon_g[branch]);
            totals.husk_c += try pool(canopy.branch_husk_carbon_g[branch]);
            totals.ear_c += try pool(canopy.branch_ear_carbon_g[branch]);
            totals.grain_c += try pool(canopy.branch_grain_carbon_g[branch]);
            totals.symbiont_mobile_c += try pool(canopy.branch_symbiont_mobile_carbon_g[branch]);
            totals.symbiont_structural_c += try pool(canopy.branch_symbiont_structural_carbon_g[branch]);
        }
        totals.standing_dead_c += try pool(canopy.plant_standing_dead_carbon_g[plant]);
    }
    return totals;
}

fn calculateFractions(canopy: *const Canopy.State, plant: usize, totals: Totals, area_m2: f64, timestep_h: f64, woody: bool, parameters: FireScience.CombustionParameters) !Fractions {
    const live_temperature_k = canopy.plant_canopy_aerodynamic_temperature_k[plant];
    const dead_temperature_k = canopy.plant_standing_dead_surface_temperature_k[plant];
    const p1 = parameters.mobile_and_leaf_specific_combustion_g_c_per_m2_h;
    const p2 = parameters.nonwoody_structural_specific_combustion_g_c_per_m2_h;
    const sheath_rate = if (woody) parameters.root_structural_specific_combustion_g_c_per_m2_h else p2;
    const stalk_rate = if (woody) parameters.woody_structural_specific_combustion_g_c_per_m2_h else p2;
    return .{
        .mobile = try FireScience.combustionFraction(totals.mobile_c, live_temperature_k, area_m2, timestep_h, p1, parameters),
        .leaf = try FireScience.combustionFraction(totals.leaf_c, live_temperature_k, area_m2, timestep_h, p1, parameters),
        .sheath = try FireScience.combustionFraction(totals.sheath_c, live_temperature_k, area_m2, timestep_h, sheath_rate, parameters),
        .stalk = try FireScience.combustionFraction(totals.stalk_c, live_temperature_k, area_m2, timestep_h, stalk_rate, parameters),
        .husk = try FireScience.combustionFraction(totals.husk_c, live_temperature_k, area_m2, timestep_h, p1, parameters),
        .ear = try FireScience.combustionFraction(totals.ear_c, live_temperature_k, area_m2, timestep_h, p2, parameters),
        .grain = try FireScience.combustionFraction(totals.grain_c, live_temperature_k, area_m2, timestep_h, p2, parameters),
        .symbiont_mobile = try FireScience.combustionFraction(totals.symbiont_mobile_c, live_temperature_k, area_m2, timestep_h, p1, parameters),
        .symbiont_structural = try FireScience.combustionFraction(totals.symbiont_structural_c, live_temperature_k, area_m2, timestep_h, p2, parameters),
        .standing_dead = try FireScience.combustionFraction(totals.standing_dead_c, dead_temperature_k, area_m2, timestep_h, parameters.standing_dead_specific_combustion_g_c_per_m2_h, parameters),
    };
}

fn validatePlant(canopy: *const Canopy.State, plant: usize, fractions: Fractions, dynamic_salts: bool) !void {
    inline for (@typeInfo(Fractions).@"struct".fields) |field| {
        const value = @field(fractions, field.name);
        if (!std.math.isFinite(value) or value < 0 or value > 1) return error.InvalidShootFireFraction;
    }
    const branches = try canopy.branchRange(plant);
    for (branches.first..branches.end) |branch| {
        inline for (.{
            "branch_mobile_", "branch_leaf_", "branch_sheath_", "branch_stalk_",           "branch_reserve_",
            "branch_husk_",   "branch_ear_",  "branch_grain_",  "branch_symbiont_mobile_", "branch_symbiont_structural_",
        }) |prefix| inline for (.{ "carbon_g", "nitrogen_g", "phosphorus_g" }) |suffix| {
            const value = @field(canopy, prefix ++ suffix)[branch];
            if (!std.math.isFinite(value) or value < 0) return error.InvalidShootFirePool;
        };
        if (dynamic_salts) for (0..8) |salt| if (!std.math.isFinite(canopy.branch_salt_content_by_species_mol[branch * 8 + salt]) or canopy.branch_salt_content_by_species_mol[branch * 8 + salt] < 0) return error.InvalidShootFirePool;
    }
    inline for (.{
        canopy.plant_standing_dead_carbon_g[plant],
        canopy.plant_standing_dead_nitrogen_g[plant],
        canopy.plant_standing_dead_phosphorus_g[plant],
    }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidShootFirePool;
}

fn commitPlant(canopy: *Canopy.State, roots: *Roots.State, plant: usize, fractions: Fractions, dynamic_salts: bool) void {
    var emitted_c: f64 = 0;
    var emitted_n: f64 = 0;
    var emitted_p: f64 = 0;
    const branches = canopy.branchRange(plant) catch unreachable;
    for (branches.first..branches.end) |branch| {
        burnPool(canopy, branch, "branch_mobile_", fractions.mobile, &emitted_c, &emitted_n, &emitted_p);
        burnPool(canopy, branch, "branch_leaf_", fractions.leaf, &emitted_c, &emitted_n, &emitted_p);
        burnPool(canopy, branch, "branch_sheath_", fractions.sheath, &emitted_c, &emitted_n, &emitted_p);
        burnPool(canopy, branch, "branch_stalk_", fractions.stalk, &emitted_c, &emitted_n, &emitted_p);
        burnPool(canopy, branch, "branch_reserve_", fractions.stalk, &emitted_c, &emitted_n, &emitted_p);
        burnPool(canopy, branch, "branch_husk_", fractions.husk, &emitted_c, &emitted_n, &emitted_p);
        burnPool(canopy, branch, "branch_ear_", fractions.ear, &emitted_c, &emitted_n, &emitted_p);
        burnPool(canopy, branch, "branch_grain_", fractions.grain, &emitted_c, &emitted_n, &emitted_p);
        burnPool(canopy, branch, "branch_symbiont_mobile_", fractions.symbiont_mobile, &emitted_c, &emitted_n, &emitted_p);
        burnPool(canopy, branch, "branch_symbiont_structural_", fractions.symbiont_structural, &emitted_c, &emitted_n, &emitted_p);
        canopy.branch_leaf_area_m2[branch] *= 1 - fractions.leaf;
        canopy.branch_sapwood_carbon_g[branch] *= 1 - fractions.stalk;
        canopy.branch_senescing_stalk_carbon_g[branch] *= 1 - fractions.stalk;
        canopy.branch_senescing_stalk_nitrogen_g[branch] *= 1 - fractions.stalk;
        canopy.branch_senescing_stalk_phosphorus_g[branch] *= 1 - fractions.stalk;
        const nodes = canopy.nodeRange(branch) catch unreachable;
        for (nodes.first..nodes.end) |node| {
            scaleNode(canopy, node, "leaf", fractions.leaf);
            scaleNode(canopy, node, "sheath", fractions.sheath);
            scaleNode(canopy, node, "internode", fractions.stalk);
            canopy.node_leaf_protein_g[node] *= 1 - fractions.leaf;
            canopy.node_sheath_protein_g[node] *= 1 - fractions.sheath;
            canopy.node_c3_nonstructural_carbon_g[node] *= 1 - fractions.mobile;
            canopy.node_c4_mesophyll_nonstructural_carbon_g[node] *= 1 - fractions.mobile;
            canopy.node_bundle_sheath_co2_carbon_g[node] *= 1 - fractions.mobile;
            canopy.node_bundle_sheath_bicarbonate_carbon_g[node] *= 1 - fractions.mobile;
            canopy.node_leaf_area_m2[node] *= 1 - fractions.leaf;
            canopy.node_sheath_height_m[node] *= 1 - fractions.sheath;
            const samples = canopy.sampleRange(node) catch unreachable;
            for (samples.first..samples.end) |sample| {
                canopy.sample_leaf_area_m2[sample] *= 1 - fractions.leaf;
                canopy.sample_exposed_leaf_area_m2[sample] *= 1 - fractions.leaf;
                canopy.sample_leaf_carbon_g[sample] *= 1 - fractions.leaf;
                canopy.sample_leaf_nitrogen_g[sample] *= 1 - fractions.leaf;
                canopy.sample_leaf_phosphorus_g[sample] *= 1 - fractions.leaf;
                canopy.sample_stalk_area_m2[sample] *= 1 - fractions.stalk;
            }
        }
        if (dynamic_salts) {
            for (0..8) |salt| {
                const salt_index = branch * 8 + salt;
                const burned = canopy.branch_salt_content_by_species_mol[salt_index] * fractions.mobile;
                canopy.branch_salt_content_by_species_mol[salt_index] -= burned;
                canopy.branch_combustion_salt_loss_by_species_mol_per_h[salt_index] += burned;
            }
        }
    }
    const live_emitted_c = emitted_c;
    inline for (.{ "carbon", "nitrogen", "phosphorus" }, .{ &emitted_c, &emitted_n, &emitted_p }) |element, emitted| {
        const field_name = "plant_standing_dead_" ++ element ++ "_g";
        const burned = @field(canopy, field_name)[plant] * fractions.standing_dead;
        @field(canopy, field_name)[plant] -= burned;
        emitted.* += burned;
        for (0..4) |kinetic| @field(canopy, "plant_standing_dead_" ++ element ++ "_by_kinetic_g")[plant * 4 + kinetic] *= 1 - fractions.standing_dead;
    }
    canopy.plant_standing_dead_height_m[plant] *= 1 - fractions.standing_dead;
    canopy.plant_live_combustion_g_c_per_h[plant] = live_emitted_c;
    canopy.plant_standing_dead_combustion_g_c_per_h[plant] = emitted_c - live_emitted_c;
    inline for (.{ "carbon", "nitrogen", "phosphorus" }) |element| {
        @field(canopy, "plant_mobile_" ++ element ++ "_g")[plant] *= 1 - fractions.mobile;
        @field(canopy, "plant_symbiont_mobile_" ++ element ++ "_g")[plant] *= 1 - fractions.symbiont_mobile;
    }
    if (dynamic_salts) {
        canopy.plant_salt_content_mol[plant] *= 1 - fractions.mobile;
    }
    roots.combustion_carbon_loss_g_c_per_h[plant] -= emitted_c;
    roots.combustion_nitrogen_loss_g_n_per_h[plant] -= emitted_n;
    roots.combustion_phosphorus_loss_g_p_per_h[plant] -= emitted_p;
    canopy.plant_combustion_carbon_loss_g_c_per_h[plant] = -emitted_c;
    canopy.plant_combustion_nitrogen_loss_g_n_per_h[plant] = -emitted_n;
    canopy.plant_combustion_phosphorus_loss_g_p_per_h[plant] = -emitted_p;
}

fn burnPool(canopy: *Canopy.State, branch: usize, comptime prefix: []const u8, fraction: f64, emitted_c: *f64, emitted_n: *f64, emitted_p: *f64) void {
    inline for (.{ "carbon", "nitrogen", "phosphorus" }, .{ emitted_c, emitted_n, emitted_p }) |element, emitted| {
        const values = @field(canopy, prefix ++ element ++ "_g");
        const burned = values[branch] * fraction;
        values[branch] -= burned;
        emitted.* += burned;
    }
}

fn scaleNode(canopy: *Canopy.State, node: usize, comptime organ: []const u8, fraction: f64) void {
    inline for (.{ "carbon_g", "nitrogen_g", "phosphorus_g" }) |suffix| @field(canopy, "node_" ++ organ ++ "_" ++ suffix)[node] *= 1 - fraction;
}

fn pool(value: f64) !f64 {
    if (!std.math.isFinite(value) or value < 0) return error.InvalidShootFirePool;
    return value;
}

test "GROSUB shoot fire conserves arbitrary-species C N P and scales topology mirrors" {
    var canopy = try Canopy.State.init(std.testing.allocator, 1, 2, &.{ 1, 1 }, &.{ 1, 1 }, &.{ 1, 1 });
    defer canopy.deinit();
    var controls = try CanopyControls.init(std.testing.allocator, 2);
    defer controls.deinit();
    controls.biomass_turnover_type[1] = 2;
    controls.root_profile_type[1] = 2;
    var roots = try Roots.State.init(std.testing.allocator, 2, 1, 1);
    defer roots.deinit();
    var surface_organic = try SoilOrganic.State.init(std.testing.allocator, 1);
    defer surface_organic.deinit();
    var initial_c: f64 = 0;
    var initial_n: f64 = 0;
    var initial_p: f64 = 0;
    for (0..2) |plant| {
        canopy.plant_canopy_aerodynamic_temperature_k[plant] = 500;
        canopy.plant_standing_dead_surface_temperature_k[plant] = 500;
        const branch = (try canopy.branchRange(plant)).first;
        const scale: f64 = @floatFromInt(plant + 1);
        inline for (.{
            "branch_mobile_", "branch_leaf_", "branch_sheath_", "branch_stalk_",           "branch_reserve_",
            "branch_husk_",   "branch_ear_",  "branch_grain_",  "branch_symbiont_mobile_", "branch_symbiont_structural_",
        }) |prefix| {
            @field(&canopy, prefix ++ "carbon_g")[branch] = scale;
            @field(&canopy, prefix ++ "nitrogen_g")[branch] = 0.1 * scale;
            @field(&canopy, prefix ++ "phosphorus_g")[branch] = 0.01 * scale;
            initial_c += scale;
            initial_n += 0.1 * scale;
            initial_p += 0.01 * scale;
        }
        canopy.plant_standing_dead_carbon_g[plant] = 2 * scale;
        canopy.plant_standing_dead_nitrogen_g[plant] = 0.2 * scale;
        canopy.plant_standing_dead_phosphorus_g[plant] = 0.02 * scale;
        canopy.plant_charcoal_carbon_g[plant] = 0.5 * scale;
        canopy.plant_charcoal_nitrogen_g[plant] = 0.05 * scale;
        canopy.plant_charcoal_phosphorus_g[plant] = 0.005 * scale;
        initial_c += 2 * scale;
        initial_n += 0.2 * scale;
        initial_p += 0.02 * scale;
        initial_c += 0.5 * scale;
        initial_n += 0.05 * scale;
        initial_p += 0.005 * scale;
        canopy.branch_leaf_area_m2[branch] = scale;
        const node = (try canopy.nodeRange(branch)).first;
        canopy.node_leaf_area_m2[node] = scale;
        canopy.node_leaf_carbon_g[node] = scale;
        canopy.node_sheath_carbon_g[node] = scale;
        canopy.node_internode_carbon_g[node] = scale;
        const sample = (try canopy.sampleRange(node)).first;
        canopy.sample_leaf_area_m2[sample] = scale;
        canopy.sample_leaf_carbon_g[sample] = scale;
    }
    var delayed_live_heat = [_]f64{ 0, 0 };
    var delayed_dead_heat = [_]f64{ 0, 0 };
    var fire_parameters = FireScience.sourceCombustionParameters();
    fire_parameters.activation_energy_j_per_mol = 5000;
    fire_parameters.arrhenius_intercept = 0.3;
    // Oxygen-free combustion exercises the source anaerobic charcoal route.
    try apply(&canopy, &controls, &roots, &.{0.01}, &.{true}, 1, false, fire_parameters, &.{600}, &.{100}, 0, 2, 0, &delayed_live_heat, &delayed_dead_heat, &surface_organic);
    var remaining_c: f64 = 0;
    var remaining_n: f64 = 0;
    var remaining_p: f64 = 0;
    var emitted_c: f64 = 0;
    var emitted_n: f64 = 0;
    var emitted_p: f64 = 0;
    var gaseous_c: f64 = 0;
    for (0..2) |plant| {
        const branch = (try canopy.branchRange(plant)).first;
        inline for (.{
            "branch_mobile_", "branch_leaf_", "branch_sheath_", "branch_stalk_",           "branch_reserve_",
            "branch_husk_",   "branch_ear_",  "branch_grain_",  "branch_symbiont_mobile_", "branch_symbiont_structural_",
        }) |prefix| {
            remaining_c += @field(&canopy, prefix ++ "carbon_g")[branch];
            remaining_n += @field(&canopy, prefix ++ "nitrogen_g")[branch];
            remaining_p += @field(&canopy, prefix ++ "phosphorus_g")[branch];
        }
        remaining_c += canopy.plant_standing_dead_carbon_g[plant];
        remaining_n += canopy.plant_standing_dead_nitrogen_g[plant];
        remaining_p += canopy.plant_standing_dead_phosphorus_g[plant];
        remaining_c += canopy.plant_charcoal_carbon_g[plant];
        remaining_n += canopy.plant_charcoal_nitrogen_g[plant];
        remaining_p += canopy.plant_charcoal_phosphorus_g[plant];
        emitted_c -= roots.combustion_carbon_loss_g_c_per_h[plant];
        emitted_n -= roots.combustion_nitrogen_loss_g_n_per_h[plant];
        emitted_p -= roots.combustion_phosphorus_loss_g_p_per_h[plant];
        gaseous_c += canopy.plant_fire_carbon_dioxide_emission_g_c_per_h[plant] + canopy.plant_fire_methane_emission_g_c_per_h[plant];
        try std.testing.expectApproxEqAbs(roots.combustion_carbon_loss_g_c_per_h[plant], canopy.plant_combustion_carbon_loss_g_c_per_h[plant], 1e-12);
        const node = (try canopy.nodeRange(branch)).first;
        try std.testing.expect(canopy.branch_leaf_area_m2[branch] < @as(f64, @floatFromInt(plant + 1)));
        try std.testing.expect(canopy.node_leaf_area_m2[node] < @as(f64, @floatFromInt(plant + 1)));
    }
    remaining_c += try surface_organic.charcoalCarbon_g_c(0);
    try std.testing.expectApproxEqAbs(initial_c, remaining_c + gaseous_c, 1e-12);
    try std.testing.expectApproxEqAbs(initial_n, remaining_n + emitted_n, 1e-12);
    try std.testing.expectApproxEqAbs(initial_p, remaining_p + emitted_p, 1e-12);
    try std.testing.expect(delayed_live_heat[0] > 0);
    try std.testing.expect(delayed_dead_heat[0] > 0);
    try std.testing.expect(try surface_organic.charcoalCarbon_g_c(0) > 0);

    // EXTRACT-004 characterization. The three closures above are NOT symmetric,
    // and the asymmetry is the defect. Carbon closes against a real sink,
    // because line 410 adds back the charcoal that `apply` actually committed to
    // surface organic matter. Nitrogen and phosphorus close against
    // `roots.combustion_*_loss_g_*_per_h`, which is the loss ledger `apply`
    // wrote at :294--295 -- the same number, moved to the other side of the
    // equation. A budget closed against its own loss term is satisfied by
    // construction and cannot detect a missing sink.
    //
    // `extract.f:276--277` credits the mineral halves of the combusted N and P
    // back to the litter layer as NH4 and H2PO4. `combustPlantPools` already
    // computes that split with the exact source fractions from `extract.f:45`
    // and 256--259. `apply` never consumes it. The mass below is therefore
    // computed by production, discarded by production, and invisible to the
    // three assertions above.
    //
    // This test asserts the CURRENT behaviour so the gap has a number. When A1
    // binds the mineral return (`docs/binding_requests/A2_grosub_extract_batch4.md`
    // section 5), this block must be replaced by an assertion that the litter
    // ammonium and phosphate owners received exactly these quantities, and the
    // N and P closures above must be rewritten to balance against that sink
    // rather than against the loss ledger.
    {
        // Source-order fractions at the combustion temperature response the
        // fire above resolved to. Both plants share a cell, so one response
        // governs the whole cell partition.
        const fire_products = try PlantSoilExchange.canopyFireCombustion(
            emitted_c,
            0,
            600,
            0,
            2,
            0,
            .{
                .gas_constant_j_per_mol_k = fire_parameters.gas_constant_j_per_mol_k,
                .combustion_activation_energy_j_per_mol = fire_parameters.activation_energy_j_per_mol,
                .combustion_temperature_intercept = fire_parameters.arrhenius_intercept,
                .oxygen_g_per_g_carbon = fire_parameters.oxygen_g_per_g_combusted_carbon,
                .maximum_aerobic_charcoal_fraction = fire_parameters.maximum_aerobic_charcoal_fraction,
                .maximum_anaerobic_charcoal_fraction = fire_parameters.maximum_anaerobic_charcoal_fraction,
                .oxygen_half_saturation_umol_per_mol = fire_parameters.oxygen_half_saturation_umol_per_mol,
                .methane_half_saturation_umol_per_mol = fire_parameters.methane_half_saturation_umol_per_mol,
                .aerobic_combustion_energy_megajoules_per_g_carbon = fire_parameters.aerobic_combustion_energy_megajoules_per_g_c,
                .anaerobic_combustion_energy_megajoules_per_g_carbon = fire_parameters.anaerobic_combustion_energy_megajoules_per_g_c,
                .methane_combustion_energy_megajoules_per_g_carbon = fire_parameters.methane_combustion_energy_megajoules_per_g_c,
            },
        );
        const temperature_response = fire_products.combustion_temperature_response;
        // `extract.f:45` FCOMNY=0.1 FCOMNX=0.5 FCOMPY=0.7 FCOMPX=0.9, applied by
        // `extract.f:256` and `:258`.
        const ammonium_fraction = 0.1 + (0.5 - 0.1) * (1 - temperature_response);
        const phosphate_fraction = 0.7 + (0.9 - 0.7) * (1 - temperature_response);
        try std.testing.expect(ammonium_fraction >= 0.1 and ammonium_fraction <= 0.5);
        try std.testing.expect(phosphate_fraction >= 0.7 and phosphate_fraction <= 0.9);

        // What the reference returns to litter, and production does not.
        const discarded_ammonium_g_n = emitted_n * ammonium_fraction;
        const discarded_phosphate_g_p = emitted_p * phosphate_fraction;
        try std.testing.expect(discarded_ammonium_g_n > 0);
        try std.testing.expect(discarded_phosphate_g_p > 0);

        // The mineral share is a majority of the phosphorus and a large
        // minority of the nitrogen, so this is not a rounding term.
        try std.testing.expect(discarded_phosphate_g_p > 0.7 * emitted_p);
        try std.testing.expect(discarded_ammonium_g_n > 0.1 * emitted_n);

        // No litter mineral sink exists to receive it. Surface organic matter
        // gained only charcoal carbon: `apply` writes exactly one surface field
        // (`:150--154`) and it is carbon. There is no shoot-fire writer of any
        // litter ammonium or phosphate owner, which is what makes the mass above
        // lost rather than merely relocated.
        try std.testing.expectEqual(
            @as(f64, 0),
            surface_organic.structural[SoilOrganic.structural_fraction_count - 1].nitrogen_g_n,
        );
        try std.testing.expectEqual(
            @as(f64, 0),
            surface_organic.structural[SoilOrganic.structural_fraction_count - 1].phosphorus_g_p,
        );
    }
}
