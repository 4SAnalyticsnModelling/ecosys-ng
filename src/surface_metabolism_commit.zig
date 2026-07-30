const std = @import("std");
const compute = @import("compute.zig");
const gas = @import("gas_transport.zig");
const organic = @import("soil_organic_initialization.zig");
const chemistry = @import("surface_litter_chemistry.zig");
const respiration = @import("surface_microbial_respiration_step.zig");
const oxygen = @import("surface_microbial_oxygen_driver.zig");
const fixation = @import("surface_nonsymbiotic_nitrogen_fixation_step.zig");
const uptake = @import("surface_microbial_substrate_uptake_step.zig");
const denitrification = @import("surface_denitrification_step.zig");
const assimilation = @import("surface_microbial_assimilation_step.zig");
const mineral_exchange = @import("surface_microbial_mineral_exchange_step.zig");
const topsoil_exchange = @import("surface_topsoil_mineral_exchange_step.zig");
const turnover = @import("surface_microbial_turnover_step.zig");
const priming = @import("surface_organic_priming_step.zig");
const organic_decomposition = @import("surface_organic_decomposition_step.zig");
const organic_sorption = @import("surface_organic_sorption_step.zig");
const litter_colonization = @import("surface_litter_colonization_step.zig");
const metabolism = @import("soil_microbial_metabolism.zig");
const soil_chemistry = @import("solute_chemistry_state.zig");
const grid = @import("grid.zig");
const zone_classification = @import("solute_charge_classification.zig");

pub const ApplyContext = struct {
    surface_organic: *organic.State,
    litter_chemistry: *chemistry.State,
    litter_gas: *gas.State,
    litter_water_m3: []const f64,
    respiration: *const respiration.State,
    oxygen: *const oxygen.State,
    nitrogen_fixation: *const fixation.State,
    substrate_uptake: *const uptake.State,
    denitrification: *denitrification.State,
    assimilation: *const assimilation.State,
    mineral_exchange: *const mineral_exchange.State,
    topsoil_exchange: *const topsoil_exchange.State,
    turnover: *const turnover.State,
    priming: *const priming.State,
    organic_decomposition: *const organic_decomposition.State,
    organic_sorption: *const organic_sorption.State,
    litter_colonization: *const litter_colonization.State,
    topsoil_organic: *organic.State,
    topsoil_humus_partition: []const [2]f64,
    topsoil_chemistry: *soil_chemistry.State,
    model_grid: *const grid.GridState,
    zone_fractions: zone_classification.ZoneFractions,
    microbial_parameters: respiration.Parameters,
    nitrogen_molar_mass_g_per_mol: f64,
    tolerance: f64,
    phosphorus_molar_mass_g_per_mol: f64,
    hourly_signed_heterotrophic_respiration_g_c: ?[]f64 = null,
    hourly_carbon_dioxide_production_g_c: ?[]f64 = null,
};

/// Atomic per-cell NITRO redistribution for the translated surface metabolism
/// block. All source sufficiency and finite-result checks precede mutation.
pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    try validate(context.*, range);
    for (range.first..range.end) |cell| try commitCell(context, cell);
}

fn commitCell(context: *ApplyContext, cell: usize) !void {
    var dissolved_after: [respiration.litter_complex_count]organic.ElementPool = undefined;
    var acetate_after: [respiration.litter_complex_count]f64 = undefined;
    var nonstructural_after: [respiration.unit_count_per_cell]organic.ElementPool = undefined;
    var structural_after: [respiration.unit_count_per_cell * assimilation.structural_component_count]organic.ElementPool = undefined;
    var residue_after: [respiration.litter_complex_count * turnover.structural_component_count]organic.ElementPool = undefined;
    var substrate_structural_after: [respiration.litter_complex_count * organic.structural_fraction_count]organic.ElementPool = undefined;
    var colonized_structural_after: [respiration.litter_complex_count * organic.structural_fraction_count]f64 = undefined;
    var adsorbed_after: [respiration.litter_complex_count]organic.ElementPool = undefined;
    var adsorbed_acetate_after: [respiration.litter_complex_count]f64 = undefined;
    var total_co2_g_c: f64 = 0;
    var total_ch4_g_c: f64 = 0;
    var total_h2_g_h: f64 = 0;
    var total_fixed_n_g_n: f64 = 0;
    var total_ammonium_exchange_g_n: f64 = 0;
    var total_nitrate_exchange_g_n: f64 = 0;
    var total_h2po4_exchange_g_p: f64 = 0;
    var total_hpo4_exchange_g_p: f64 = 0;
    const top_organic_layer = try context.model_grid.layerIndex(cell, 0);
    var topsoil_humus_after = [2]organic.ElementPool{
        context.topsoil_organic.structural[(top_organic_layer * organic.substrate_count + 4) * organic.structural_fraction_count],
        context.topsoil_organic.structural[(top_organic_layer * organic.substrate_count + 4) * organic.structural_fraction_count + 1],
    };
    var topsoil_humus_colonized_after = [2]f64{
        context.topsoil_organic.colonized_structural_carbon_g_c[(top_organic_layer * organic.substrate_count + 4) * organic.structural_fraction_count],
        context.topsoil_organic.colonized_structural_carbon_g_c[(top_organic_layer * organic.substrate_count + 4) * organic.structural_fraction_count + 1],
    };
    const topsoil_particulate_index = (top_organic_layer * organic.substrate_count + 3) * organic.structural_fraction_count;
    var topsoil_particulate_after = context.topsoil_organic.structural[topsoil_particulate_index];
    var topsoil_particulate_colonized_after = context.topsoil_organic.colonized_structural_carbon_g_c[topsoil_particulate_index];

    for (0..respiration.litter_complex_count) |complex| {
        const compact = cell * respiration.litter_complex_count + complex;
        var decomposition_dissolved_input: organic.ElementPool = .{};
        for (0..organic.structural_fraction_count) |fraction| {
            const local = complex * organic.structural_fraction_count + fraction;
            const surface_index = (cell * organic.substrate_count + complex) * organic.structural_fraction_count + fraction;
            const decomposed = context.organic_decomposition.structural_decomposition[compact * organic.structural_fraction_count + fraction];
            substrate_structural_after[local] = subtractPool(context.surface_organic.structural[surface_index], decomposed);
            colonized_structural_after[local] = context.surface_organic.colonized_structural_carbon_g_c[surface_index] - decomposed.carbon_g_c + context.litter_colonization.colonized_carbon_increment_g_c[compact * organic.structural_fraction_count + fraction];
            decomposition_dissolved_input = addPool(decomposition_dissolved_input, context.organic_decomposition.dissolved_structural_products[compact * organic.structural_fraction_count + fraction]);
            const particulate = context.organic_decomposition.particulate_products[compact * organic.structural_fraction_count + fraction];
            topsoil_particulate_after = addPool(topsoil_particulate_after, particulate);
            topsoil_particulate_colonized_after += particulate.carbon_g_c;
        }
        for (0..turnover.structural_component_count) |component| residue_after[complex * turnover.structural_component_count + component] = context.surface_organic.residue[(cell * organic.substrate_count + complex) * organic.residue_fraction_count + component];
        for (0..organic.residue_fraction_count) |fraction| decomposition_dissolved_input = addPool(decomposition_dissolved_input, context.organic_decomposition.microbial_residue_decomposition[compact * organic.residue_fraction_count + fraction]);
        decomposition_dissolved_input = addPool(decomposition_dissolved_input, context.organic_decomposition.sorbed_organic_decomposition[compact]);
        const mobile = cell * organic.substrate_count + complex;
        const adsorbed_before_sorption = subtractPool(context.surface_organic.adsorbed[mobile], context.organic_decomposition.sorbed_organic_decomposition[compact]);
        const adsorbed_acetate_before_sorption = context.surface_organic.adsorbed_acetate_carbon_g_c[mobile] - context.organic_decomposition.sorbed_acetate_decomposition_g_c[compact];
        var sorption: organic.ElementPool = .{ .carbon_g_c = context.organic_sorption.doc_sorption_g_c[compact], .nitrogen_g_n = context.organic_sorption.don_sorption_g_n[compact], .phosphorus_g_p = context.organic_sorption.dop_sorption_g_p[compact] };
        var doc_c: f64 = 0;
        var don_n: f64 = 0;
        var dop_p: f64 = 0;
        var acetate_c: f64 = 0;
        var fermentation_c: f64 = 0;
        for (0..respiration.source_population_count) |population| {
            const unit_local = complex * respiration.source_population_count + population;
            const unit = cell * respiration.unit_count_per_cell + unit_local;
            doc_c += context.substrate_uptake.doc_uptake_g_c[unit];
            don_n += context.substrate_uptake.dissolved_organic_nitrogen_uptake_g_n[unit];
            dop_p += context.substrate_uptake.dissolved_organic_phosphorus_uptake_g_p[unit];
            acetate_c += context.substrate_uptake.acetate_uptake_g_c[unit];
            const oxygen_fraction = if (context.oxygen.populations[unit].is_aerobic) context.oxygen.allocation.demand_satisfaction_fraction[unit] else 1;
            const actual_respiration_g_c = context.respiration.substrate_limited_respiration_g_c[unit] * oxygen_fraction;
            switch (context.microbial_parameters.populations[population].metabolism) {
                .aerobic_heterotroph => total_co2_g_c += actual_respiration_g_c,
                .fermenting_heterotroph => {
                    total_co2_g_c += 0.333 * actual_respiration_g_c;
                    fermentation_c += 0.667 * actual_respiration_g_c;
                    total_h2_g_h += 0.111 * actual_respiration_g_c;
                },
                .acetotrophic_methanogen => {
                    total_co2_g_c += 0.5 * actual_respiration_g_c;
                    total_ch4_g_c += 0.5 * actual_respiration_g_c;
                },
            }
            total_co2_g_c += context.substrate_uptake.denitrification_respiration_g_c[unit] + context.nitrogen_fixation.fixation_respiration_g_c[unit];
            total_fixed_n_g_n += context.nitrogen_fixation.fixed_nitrogen_g_n[unit];
            const mineral_n = context.mineral_exchange.ammonium_exchange_g_n[unit] + context.mineral_exchange.nitrate_exchange_g_n[unit];
            const mineral_p = context.mineral_exchange.h2po4_exchange_g_p[unit] + context.mineral_exchange.hpo4_exchange_g_p[unit];
            const topsoil_n = context.topsoil_exchange.ammonium_exchange_g_n[unit] + context.topsoil_exchange.nitrate_exchange_g_n[unit];
            const topsoil_p = context.topsoil_exchange.h2po4_exchange_g_p[unit] + context.topsoil_exchange.hpo4_exchange_g_p[unit];
            total_ammonium_exchange_g_n += context.mineral_exchange.ammonium_exchange_g_n[unit];
            total_nitrate_exchange_g_n += context.mineral_exchange.nitrate_exchange_g_n[unit];
            total_h2po4_exchange_g_p += context.mineral_exchange.h2po4_exchange_g_p[unit];
            total_hpo4_exchange_g_p += context.mineral_exchange.hpo4_exchange_g_p[unit];

            const microbial = ((cell * organic.microbial_substrate_count + complex) * organic.microbial_population_count + population) * organic.kinetic_fraction_count + 2;
            const current = context.surface_organic.microbial[microbial];
            var structural_transfer: organic.ElementPool = .{};
            var turnover_recycled_to_nonstructural: organic.ElementPool = .{};
            for (0..assimilation.structural_component_count) |component| {
                const transfer = context.assimilation.transfer[(unit * assimilation.structural_component_count) + component];
                structural_transfer.carbon_g_c += transfer.carbon_g_c;
                structural_transfer.nitrogen_g_n += transfer.nitrogen_g_n;
                structural_transfer.phosphorus_g_p += transfer.phosphorus_g_p;
                const current_structural = context.surface_organic.microbial[microbial - 2 + component];
                const priming_change = context.priming.exchange.microbial_change[cell * respiration.unit_count_per_cell * organic.kinetic_fraction_count + unit_local * organic.kinetic_fraction_count + component];
                const turnover_index = unit * turnover.structural_component_count + component;
                const basal = context.turnover.basal[turnover_index];
                const senescence = context.turnover.senescence[turnover_index];
                structural_after[(unit_local * assimilation.structural_component_count) + component] = .{ .carbon_g_c = current_structural.carbon_g_c + transfer.carbon_g_c - basal.decomposed.carbon_g_c - senescence.decomposed.carbon_g_c + priming_change.carbon_g_c, .nitrogen_g_n = current_structural.nitrogen_g_n + transfer.nitrogen_g_n - basal.decomposed.nitrogen_g_n - senescence.decomposed.nitrogen_g_n + priming_change.nitrogen_g_n, .phosphorus_g_p = current_structural.phosphorus_g_p + transfer.phosphorus_g_p - basal.decomposed.phosphorus_g_p - senescence.decomposed.phosphorus_g_p + priming_change.phosphorus_g_p };
                try finitePool(structural_after[(unit_local * assimilation.structural_component_count) + component]);
                turnover_recycled_to_nonstructural = addPool(turnover_recycled_to_nonstructural, fromMetabolic(basal.recycled));
                turnover_recycled_to_nonstructural.nitrogen_g_n += senescence.recycled.nitrogen_g_n;
                turnover_recycled_to_nonstructural.phosphorus_g_p += senescence.recycled.phosphorus_g_p;
                total_co2_g_c += senescence.recycled.carbon_g_c;
                const residue_index = complex * turnover.structural_component_count + component;
                residue_after[residue_index] = addPool(residue_after[residue_index], fromMetabolic(basal.microbial_residue));
                residue_after[residue_index] = addPool(residue_after[residue_index], fromMetabolic(senescence.microbial_residue));
                const humified = addPool(fromMetabolic(basal.humified), fromMetabolic(senescence.humified));
                for (0..2) |humus_class| {
                    const humus_input = scalePool(humified, context.topsoil_humus_partition[cell][humus_class]);
                    topsoil_humus_after[humus_class] = addPool(topsoil_humus_after[humus_class], humus_input);
                    topsoil_humus_colonized_after[humus_class] += humus_input.carbon_g_c;
                }
            }
            if (structural_transfer.carbon_g_c > current.carbon_g_c + context.tolerance or structural_transfer.nitrogen_g_n > current.nitrogen_g_n + context.tolerance or structural_transfer.phosphorus_g_p > current.phosphorus_g_p + context.tolerance) return error.InsufficientNonstructuralSurfaceMicrobialPool;
            nonstructural_after[unit_local] = .{
                .carbon_g_c = current.carbon_g_c - structural_transfer.carbon_g_c + context.substrate_uptake.nonstructural_carbon_gain_g_c[unit],
                .nitrogen_g_n = current.nitrogen_g_n - structural_transfer.nitrogen_g_n + context.substrate_uptake.dissolved_organic_nitrogen_uptake_g_n[unit] + context.nitrogen_fixation.fixed_nitrogen_g_n[unit] + mineral_n,
                .phosphorus_g_p = current.phosphorus_g_p - structural_transfer.phosphorus_g_p + context.substrate_uptake.dissolved_organic_phosphorus_uptake_g_p[unit] + mineral_p + topsoil_p,
            };
            nonstructural_after[unit_local] = addPool(nonstructural_after[unit_local], turnover_recycled_to_nonstructural);
            nonstructural_after[unit_local] = addPool(nonstructural_after[unit_local], context.priming.exchange.microbial_change[cell * respiration.unit_count_per_cell * organic.kinetic_fraction_count + unit_local * organic.kinetic_fraction_count + 2]);
            nonstructural_after[unit_local].nitrogen_g_n += topsoil_n;
            try finitePool(nonstructural_after[unit_local]);
        }
        const current = context.surface_organic.dissolved[mobile];
        if (doc_c > current.carbon_g_c + context.tolerance or don_n > current.nitrogen_g_n + context.tolerance or dop_p > current.phosphorus_g_p + context.tolerance or acetate_c > context.surface_organic.dissolved_acetate_carbon_g_c[mobile] + context.tolerance) return error.InsufficientSurfaceMicrobialSubstrate;
        const priming_index = cell * priming.substrate_count + complex;
        dissolved_after[complex] = addPool(addPool(.{ .carbon_g_c = current.carbon_g_c - doc_c, .nitrogen_g_n = current.nitrogen_g_n - don_n, .phosphorus_g_p = current.phosphorus_g_p - dop_p }, context.priming.exchange.dissolved_change[priming_index]), decomposition_dissolved_input);
        sorption.carbon_g_c = boundedExchange(sorption.carbon_g_c, dissolved_after[complex].carbon_g_c, adsorbed_before_sorption.carbon_g_c);
        sorption.nitrogen_g_n = boundedExchange(sorption.nitrogen_g_n, dissolved_after[complex].nitrogen_g_n, adsorbed_before_sorption.nitrogen_g_n);
        sorption.phosphorus_g_p = boundedExchange(sorption.phosphorus_g_p, dissolved_after[complex].phosphorus_g_p, adsorbed_before_sorption.phosphorus_g_p);
        adsorbed_after[complex] = addPool(adsorbed_before_sorption, sorption);
        dissolved_after[complex] = subtractPool(dissolved_after[complex], sorption);
        acetate_after[complex] = context.surface_organic.dissolved_acetate_carbon_g_c[mobile] - acetate_c + fermentation_c + context.priming.exchange.acetate_change_g_c[priming_index] + context.organic_decomposition.sorbed_acetate_decomposition_g_c[compact];
        const acetate_sorption_g_c = boundedExchange(context.organic_sorption.acetate_sorption_g_c[compact], acetate_after[complex], adsorbed_acetate_before_sorption);
        adsorbed_acetate_after[complex] = adsorbed_acetate_before_sorption + acetate_sorption_g_c;
        acetate_after[complex] -= acetate_sorption_g_c;
        try finitePool(dissolved_after[complex]);
        if (!std.math.isFinite(acetate_after[complex]) or acetate_after[complex] < 0) return error.NonFiniteSurfaceMetabolismCommit;
    }

    const chemo_don_g_n = context.denitrification.chemodenitrification_dissolved_organic_nitrogen_production_g_n[cell];
    var total_residue_carbon_g_c: f64 = 0;
    for (residue_after) |pool| total_residue_carbon_g_c += pool.carbon_g_c;
    for (0..respiration.litter_complex_count) |complex| {
        var complex_residue_carbon_g_c: f64 = 0;
        for (0..turnover.structural_component_count) |component| complex_residue_carbon_g_c += residue_after[complex * turnover.structural_component_count + component].carbon_g_c;
        const residue_fraction: f64 = if (total_residue_carbon_g_c > context.tolerance) complex_residue_carbon_g_c / total_residue_carbon_g_c else if (complex == 0) 1.0 else 0.0;
        dissolved_after[complex].nitrogen_g_n += chemo_don_g_n * residue_fraction;
    }

    const water_m3 = context.litter_water_m3[cell];
    const nitrate_before_g_n = context.litter_chemistry.cells[cell].nitrate_mol_per_m3 * water_m3 * context.nitrogen_molar_mass_g_per_mol;
    const ammonium_before_g_n = context.litter_chemistry.cells[cell].ammonium_mol_per_m3 * water_m3 * context.nitrogen_molar_mass_g_per_mol;
    const h2po4_before_g_p = context.litter_chemistry.cells[cell].h2po4_mol_p_per_m3 * water_m3 * context.phosphorus_molar_mass_g_per_mol;
    const hpo4_before_g_p = context.litter_chemistry.cells[cell].hpo4_mol_p_per_m3 * water_m3 * context.phosphorus_molar_mass_g_per_mol;
    var nitrate_reduction_g_n: f64 = 0;
    var nitrite_reduction_g_n: f64 = 0;
    var n2o_reduction_g_n: f64 = 0;
    for (0..respiration.litter_complex_count) |complex| {
        const compact = cell * respiration.litter_complex_count + complex;
        nitrate_reduction_g_n += context.denitrification.nitrate_reduction_g_n[compact];
        nitrite_reduction_g_n += context.denitrification.nitrite_reduction_g_n[compact];
        n2o_reduction_g_n += context.denitrification.nitrous_oxide_reduction_g_n[compact];
    }
    const chemo_nitrite_reduction_g_n = context.denitrification.chemodenitrification_nitrite_reduction_g_n[cell];
    const chemo_n2o_production_g_n = context.denitrification.chemodenitrification_nitrous_oxide_production_g_n[cell];
    const nitrite_before_g_n = context.denitrification.nitrite_g_n[cell];
    const n2o_index = cell * gas.species_count + @intFromEnum(gas.Species.nitrous_oxide);
    const n2_index = cell * gas.species_count + @intFromEnum(gas.Species.nitrogen);
    const n2o_before_g_n = context.litter_gas.dissolved_mass_g[n2o_index];
    const n2_before_g_n = context.litter_gas.dissolved_mass_g[n2_index];
    const top = try context.model_grid.layerIndex(cell, 0);
    const top_water_m3 = context.model_grid.matrix_liquid_water_m3[top];
    var top_aqueous_after = context.topsoil_chemistry.aqueous[top];
    var top_non_band_phosphate_after = context.topsoil_chemistry.non_band_phosphate[top];
    var top_band_phosphate_after = context.topsoil_chemistry.band_phosphate[top];
    var top_ammonium_g_n: f64 = 0;
    var top_nitrate_g_n: f64 = 0;
    var top_h2po4_g_p: f64 = 0;
    var top_hpo4_g_p: f64 = 0;
    for (cell * respiration.unit_count_per_cell..(cell + 1) * respiration.unit_count_per_cell) |unit| {
        top_ammonium_g_n += context.topsoil_exchange.ammonium_exchange_g_n[unit];
        top_nitrate_g_n += context.topsoil_exchange.nitrate_exchange_g_n[unit];
        top_h2po4_g_p += context.topsoil_exchange.h2po4_exchange_g_p[unit];
        top_hpo4_g_p += context.topsoil_exchange.hpo4_exchange_g_p[unit];
    }
    const top_ammonium_available = top_water_m3 * context.nitrogen_molar_mass_g_per_mol * (context.zone_fractions.ammonium_non_band * top_aqueous_after.ammonium_non_band + context.zone_fractions.ammonium_band * top_aqueous_after.ammonium_band);
    const top_nitrate_available = top_water_m3 * context.nitrogen_molar_mass_g_per_mol * (context.zone_fractions.nitrate_non_band * top_aqueous_after.nitrate_non_band + context.zone_fractions.nitrate_band * top_aqueous_after.nitrate_band);
    const top_h2po4_available = top_water_m3 * context.phosphorus_molar_mass_g_per_mol * (context.zone_fractions.phosphate_non_band * top_non_band_phosphate_after.dissolved_h2po4_mol_p_per_m3 + context.zone_fractions.phosphate_band * top_band_phosphate_after.dissolved_h2po4_mol_p_per_m3);
    const top_hpo4_available = top_water_m3 * context.phosphorus_molar_mass_g_per_mol * (context.zone_fractions.phosphate_non_band * top_non_band_phosphate_after.dissolved_hpo4_mol_p_per_m3 + context.zone_fractions.phosphate_band * top_band_phosphate_after.dissolved_hpo4_mol_p_per_m3);
    if (top_ammonium_g_n > top_ammonium_available + context.tolerance or top_nitrate_g_n > top_nitrate_available + context.tolerance or top_h2po4_g_p > top_h2po4_available + context.tolerance or top_hpo4_g_p > top_hpo4_available + context.tolerance) return error.InsufficientTopsoilMineralNutrient;
    if (top_water_m3 <= 0 and top_ammonium_g_n + top_nitrate_g_n + top_h2po4_g_p + top_hpo4_g_p > context.tolerance) return error.InsufficientTopsoilMineralNutrient;
    if (top_water_m3 > 0) {
        const ammonium_delta = top_ammonium_g_n / (top_water_m3 * context.nitrogen_molar_mass_g_per_mol);
        const nitrate_delta = top_nitrate_g_n / (top_water_m3 * context.nitrogen_molar_mass_g_per_mol);
        const h2po4_delta = top_h2po4_g_p / (top_water_m3 * context.phosphorus_molar_mass_g_per_mol);
        const hpo4_delta = top_hpo4_g_p / (top_water_m3 * context.phosphorus_molar_mass_g_per_mol);
        if (context.zone_fractions.ammonium_non_band > 0) top_aqueous_after.ammonium_non_band -= ammonium_delta;
        if (context.zone_fractions.ammonium_band > 0) top_aqueous_after.ammonium_band -= ammonium_delta;
        if (context.zone_fractions.nitrate_non_band > 0) top_aqueous_after.nitrate_non_band -= nitrate_delta;
        if (context.zone_fractions.nitrate_band > 0) top_aqueous_after.nitrate_band -= nitrate_delta;
        if (context.zone_fractions.phosphate_non_band > 0) {
            top_non_band_phosphate_after.dissolved_h2po4_mol_p_per_m3 -= h2po4_delta;
            top_non_band_phosphate_after.dissolved_hpo4_mol_p_per_m3 -= hpo4_delta;
        }
        if (context.zone_fractions.phosphate_band > 0) {
            top_band_phosphate_after.dissolved_h2po4_mol_p_per_m3 -= h2po4_delta;
            top_band_phosphate_after.dissolved_hpo4_mol_p_per_m3 -= hpo4_delta;
        }
    }
    inline for (.{
        top_aqueous_after.ammonium_non_band,
        top_aqueous_after.ammonium_band,
        top_aqueous_after.nitrate_non_band,
        top_aqueous_after.nitrate_band,
        top_non_band_phosphate_after.dissolved_h2po4_mol_p_per_m3,
        top_non_band_phosphate_after.dissolved_hpo4_mol_p_per_m3,
        top_band_phosphate_after.dissolved_h2po4_mol_p_per_m3,
        top_band_phosphate_after.dissolved_hpo4_mol_p_per_m3,
    }) |value| if (!std.math.isFinite(value) or value < -context.tolerance) return error.InsufficientTopsoilMineralNutrient;
    if (total_ammonium_exchange_g_n > ammonium_before_g_n + context.tolerance or nitrate_reduction_g_n + total_nitrate_exchange_g_n > nitrate_before_g_n + context.tolerance or total_h2po4_exchange_g_p > h2po4_before_g_p + context.tolerance or total_hpo4_exchange_g_p > hpo4_before_g_p + context.tolerance or nitrite_reduction_g_n + chemo_nitrite_reduction_g_n > nitrite_before_g_n + nitrate_reduction_g_n + context.tolerance or n2o_reduction_g_n > n2o_before_g_n + nitrite_reduction_g_n + chemo_n2o_production_g_n + context.tolerance or total_fixed_n_g_n > n2_before_g_n + n2o_reduction_g_n + context.tolerance) return error.InsufficientSurfaceInorganicNitrogen;
    const ammonium_after_g_n = ammonium_before_g_n - total_ammonium_exchange_g_n;
    const nitrate_after_g_n = nitrate_before_g_n - nitrate_reduction_g_n - total_nitrate_exchange_g_n;
    const h2po4_after_g_p = h2po4_before_g_p - total_h2po4_exchange_g_p;
    const hpo4_after_g_p = hpo4_before_g_p - total_hpo4_exchange_g_p;
    const nitrite_after_g_n = @max(0, nitrite_before_g_n + nitrate_reduction_g_n - nitrite_reduction_g_n - chemo_nitrite_reduction_g_n);
    const n2o_after_g_n = @max(0, n2o_before_g_n + nitrite_reduction_g_n + chemo_n2o_production_g_n - n2o_reduction_g_n);
    const n2_after_g_n = @max(0, n2_before_g_n + n2o_reduction_g_n - total_fixed_n_g_n);
    inline for (.{ total_co2_g_c, total_ch4_g_c, total_h2_g_h, ammonium_after_g_n, nitrate_after_g_n, h2po4_after_g_p, hpo4_after_g_p, nitrite_after_g_n, n2o_after_g_n, n2_after_g_n }) |value| if (!std.math.isFinite(value) or value < -context.tolerance) return error.NonFiniteSurfaceMetabolismCommit;
    for (0..respiration.litter_complex_count) |complex| for (0..organic.residue_fraction_count) |component| {
        const residue_index = complex * organic.residue_fraction_count + component;
        residue_after[residue_index] = subtractPool(residue_after[residue_index], context.organic_decomposition.microbial_residue_decomposition[(cell * respiration.litter_complex_count + complex) * organic.residue_fraction_count + component]);
    };
    for (residue_after) |pool| try finitePool(pool);
    for (substrate_structural_after) |pool| try finitePool(pool);
    for (colonized_structural_after, substrate_structural_after) |value, substrate| if (!std.math.isFinite(value) or value < -context.tolerance or value > substrate.carbon_g_c + context.tolerance) return error.NonFiniteSurfaceMetabolismCommit;
    for (adsorbed_after) |pool| try finitePool(pool);
    for (adsorbed_acetate_after) |value| if (!std.math.isFinite(value) or value < -context.tolerance) return error.NonFiniteSurfaceMetabolismCommit;
    for (topsoil_humus_after) |pool| try finitePool(pool);
    for (topsoil_humus_colonized_after) |value| if (!std.math.isFinite(value) or value < 0) return error.NonFiniteSurfaceMetabolismCommit;
    try finitePool(topsoil_particulate_after);
    if (!std.math.isFinite(topsoil_particulate_colonized_after) or topsoil_particulate_colonized_after < 0) return error.NonFiniteSurfaceMetabolismCommit;

    for (0..respiration.litter_complex_count) |complex| {
        const mobile = cell * organic.substrate_count + complex;
        context.surface_organic.dissolved[mobile] = dissolved_after[complex];
        context.surface_organic.dissolved_acetate_carbon_g_c[mobile] = acetate_after[complex];
        context.surface_organic.adsorbed[mobile] = adsorbed_after[complex];
        context.surface_organic.adsorbed_acetate_carbon_g_c[mobile] = adsorbed_acetate_after[complex];
        for (0..organic.structural_fraction_count) |fraction| {
            const local = complex * organic.structural_fraction_count + fraction;
            const surface_index = (cell * organic.substrate_count + complex) * organic.structural_fraction_count + fraction;
            context.surface_organic.structural[surface_index] = substrate_structural_after[local];
            context.surface_organic.colonized_structural_carbon_g_c[surface_index] = @max(0, colonized_structural_after[local]);
        }
        for (0..turnover.structural_component_count) |component| context.surface_organic.residue[(cell * organic.substrate_count + complex) * organic.residue_fraction_count + component] = residue_after[complex * turnover.structural_component_count + component];
        for (0..respiration.source_population_count) |population| {
            const unit_local = complex * respiration.source_population_count + population;
            const microbial = ((cell * organic.microbial_substrate_count + complex) * organic.microbial_population_count + population) * organic.kinetic_fraction_count + 2;
            context.surface_organic.microbial[microbial] = nonstructural_after[unit_local];
            for (0..assimilation.structural_component_count) |component| context.surface_organic.microbial[microbial - 2 + component] = structural_after[(unit_local * assimilation.structural_component_count) + component];
        }
    }
    context.litter_chemistry.cells[cell].nitrate_mol_per_m3 = if (water_m3 > 0) nitrate_after_g_n / (water_m3 * context.nitrogen_molar_mass_g_per_mol) else 0;
    context.litter_chemistry.cells[cell].ammonium_mol_per_m3 = if (water_m3 > 0) ammonium_after_g_n / (water_m3 * context.nitrogen_molar_mass_g_per_mol) else 0;
    context.litter_chemistry.cells[cell].h2po4_mol_p_per_m3 = if (water_m3 > 0) h2po4_after_g_p / (water_m3 * context.phosphorus_molar_mass_g_per_mol) else 0;
    context.litter_chemistry.cells[cell].hpo4_mol_p_per_m3 = if (water_m3 > 0) hpo4_after_g_p / (water_m3 * context.phosphorus_molar_mass_g_per_mol) else 0;
    context.denitrification.nitrite_g_n[cell] = nitrite_after_g_n;
    context.litter_gas.dissolved_mass_g[n2o_index] = n2o_after_g_n;
    context.litter_gas.dissolved_mass_g[n2_index] = n2_after_g_n;
    context.topsoil_chemistry.aqueous[top] = top_aqueous_after;
    context.topsoil_chemistry.non_band_phosphate[top] = top_non_band_phosphate_after;
    context.topsoil_chemistry.band_phosphate[top] = top_band_phosphate_after;
    for (0..2) |humus_class| {
        const humus_index = (top_organic_layer * organic.substrate_count + 4) * organic.structural_fraction_count + humus_class;
        context.topsoil_organic.structural[humus_index] = topsoil_humus_after[humus_class];
        context.topsoil_organic.colonized_structural_carbon_g_c[humus_index] = topsoil_humus_colonized_after[humus_class];
    }
    context.topsoil_organic.structural[topsoil_particulate_index] = topsoil_particulate_after;
    context.topsoil_organic.colonized_structural_carbon_g_c[topsoil_particulate_index] = topsoil_particulate_colonized_after;
    context.litter_gas.dissolved_mass_g[cell * gas.species_count + @intFromEnum(gas.Species.carbon_dioxide)] += total_co2_g_c;
    context.litter_gas.dissolved_mass_g[cell * gas.species_count + @intFromEnum(gas.Species.methane)] += total_ch4_g_c;
    if (context.hourly_signed_heterotrophic_respiration_g_c) |ledger|
        ledger[cell] = -(total_co2_g_c + total_ch4_g_c);
    if (context.hourly_carbon_dioxide_production_g_c) |ledger|
        ledger[cell] = total_co2_g_c;
    context.litter_gas.dissolved_mass_g[cell * gas.species_count + @intFromEnum(gas.Species.hydrogen)] += total_h2_g_h;
}

fn finitePool(pool: organic.ElementPool) !void {
    inline for (@typeInfo(organic.ElementPool).@"struct".fields) |field| if (!std.math.isFinite(@field(pool, field.name)) or @field(pool, field.name) < 0) return error.NonFiniteSurfaceMetabolismCommit;
}

fn addPool(a: organic.ElementPool, b: organic.ElementPool) organic.ElementPool {
    return .{ .carbon_g_c = a.carbon_g_c + b.carbon_g_c, .nitrogen_g_n = a.nitrogen_g_n + b.nitrogen_g_n, .phosphorus_g_p = a.phosphorus_g_p + b.phosphorus_g_p };
}

fn scalePool(pool: organic.ElementPool, fraction: f64) organic.ElementPool {
    return .{ .carbon_g_c = pool.carbon_g_c * fraction, .nitrogen_g_n = pool.nitrogen_g_n * fraction, .phosphorus_g_p = pool.phosphorus_g_p * fraction };
}

fn subtractPool(a: organic.ElementPool, b: organic.ElementPool) organic.ElementPool {
    return .{ .carbon_g_c = a.carbon_g_c - b.carbon_g_c, .nitrogen_g_n = a.nitrogen_g_n - b.nitrogen_g_n, .phosphorus_g_p = a.phosphorus_g_p - b.phosphorus_g_p };
}

fn fromMetabolic(pool: metabolism.ElementalPool) organic.ElementPool {
    return .{ .carbon_g_c = pool.carbon_g_c, .nitrogen_g_n = pool.nitrogen_g_n, .phosphorus_g_p = pool.phosphorus_g_p };
}

/// Limits a simultaneously evaluated adsorption/desorption flux to the
/// inventory remaining after the other hourly source and sink terms.
fn boundedExchange(proposed: f64, dissolved_available: f64, sorbed_available: f64) f64 {
    if (proposed >= 0) return @min(proposed, @max(0, dissolved_available));
    return @max(proposed, -@max(0, sorbed_available));
}

fn validate(context: ApplyContext, range: compute.CellRange) !void {
    if (context.hourly_signed_heterotrophic_respiration_g_c) |ledger| if (ledger.len != context.model_grid.cell_count) return error.HeterotrophicRespirationLedgerDimensionMismatch;
    if (context.hourly_carbon_dioxide_production_g_c) |ledger| if (ledger.len != context.model_grid.cell_count) return error.CarbonDioxideProductionLedgerDimensionMismatch;
    const cells = context.surface_organic.layer_count;
    if (range.first > range.end or range.end > cells or context.litter_chemistry.cells.len != cells or context.litter_gas.cell_count != cells or context.litter_water_m3.len != cells or context.respiration.cell_count != cells or context.oxygen.cell_count != cells or context.nitrogen_fixation.cell_count != cells or context.substrate_uptake.cell_count != cells or context.denitrification.cell_count != cells or context.assimilation.cell_count != cells or context.mineral_exchange.cell_count != cells or context.topsoil_exchange.cell_count != cells or context.turnover.cell_count != cells or context.priming.cellCount() != cells or context.organic_decomposition.cell_count != cells or context.organic_sorption.cell_count != cells or context.litter_colonization.cell_count != cells or context.topsoil_humus_partition.len != cells or context.model_grid.cell_count != cells or context.topsoil_chemistry.cell_count != context.model_grid.layer_count or context.topsoil_organic.layer_count != context.model_grid.layer_count) return error.SurfaceMetabolismCommitDimensionMismatch;
    for (context.topsoil_humus_partition) |partition| if (!std.math.isFinite(partition[0]) or !std.math.isFinite(partition[1]) or partition[0] < 0 or partition[1] < 0 or @abs(partition[0] + partition[1] - 1) > context.tolerance) return error.InvalidSurfaceMetabolismCommitParameter;
    if (!std.math.isFinite(context.nitrogen_molar_mass_g_per_mol) or context.nitrogen_molar_mass_g_per_mol <= 0 or !std.math.isFinite(context.phosphorus_molar_mass_g_per_mol) or context.phosphorus_molar_mass_g_per_mol <= 0 or !std.math.isFinite(context.tolerance) or context.tolerance < 0) return error.InvalidSurfaceMetabolismCommitParameter;
}

test "surface sorption exchange cannot overdraw either shared pool" {
    try std.testing.expectEqual(@as(f64, 0.1), boundedExchange(0.4, 0.1, 2));
    try std.testing.expectEqual(@as(f64, -0.2), boundedExchange(-0.5, 2, 0.2));
}

test "surface metabolism commit conserves carbon and nitrogen" {
    var organic_state = try organic.State.init(std.testing.allocator, 1);
    defer organic_state.deinit();
    organic_state.dissolved[0] = .{ .carbon_g_c = 2, .nitrogen_g_n = 0.2, .phosphorus_g_p = 0.02 };
    organic_state.dissolved_acetate_carbon_g_c[0] = 1;
    var chemistry_state = try chemistry.State.init(std.testing.allocator, 1);
    defer chemistry_state.deinit();
    chemistry_state.cells[0].nitrate_mol_per_m3 = 0.1;
    var gas_state = try gas.State.init(std.testing.allocator, 1);
    defer gas_state.deinit();
    gas_state.dissolved_mass_g[@intFromEnum(gas.Species.nitrogen)] = 1;
    gas_state.dissolved_mass_g[@intFromEnum(gas.Species.nitrous_oxide)] = 0.2;
    var respiration_state = try respiration.State.init(std.testing.allocator, 1);
    defer respiration_state.deinit();
    respiration_state.substrate_limited_respiration_g_c[0] = 0.1;
    var oxygen_state = try oxygen.State.init(std.testing.allocator, 1);
    defer oxygen_state.deinit();
    oxygen_state.populations[0].is_aerobic = true;
    oxygen_state.allocation.demand_satisfaction_fraction[0] = 1;
    var fixation_state = try fixation.State.init(std.testing.allocator, 1);
    defer fixation_state.deinit();
    var uptake_state = try uptake.State.init(std.testing.allocator, 1);
    defer uptake_state.deinit();
    uptake_state.doc_uptake_g_c[0] = 0.3;
    uptake_state.dissolved_organic_nitrogen_uptake_g_n[0] = 0.03;
    uptake_state.dissolved_organic_phosphorus_uptake_g_p[0] = 0.003;
    uptake_state.nonstructural_carbon_gain_g_c[0] = 0.2;
    var denitrification_state = try denitrification.State.init(std.testing.allocator, 1);
    defer denitrification_state.deinit();
    var assimilation_state = try assimilation.State.init(std.testing.allocator, 1);
    defer assimilation_state.deinit();
    var mineral_exchange_state = try mineral_exchange.State.init(std.testing.allocator, 1);
    defer mineral_exchange_state.deinit();
    var topsoil_exchange_state = try topsoil_exchange.State.init(std.testing.allocator, 1);
    defer topsoil_exchange_state.deinit();
    var turnover_state = try turnover.State.init(std.testing.allocator, 1);
    defer turnover_state.deinit();
    var priming_state = try priming.State.init(std.testing.allocator, 1);
    defer priming_state.deinit();
    var organic_decomposition_state = try organic_decomposition.State.init(std.testing.allocator, 1);
    defer organic_decomposition_state.deinit();
    var organic_sorption_state = try organic_sorption.State.init(std.testing.allocator, 1);
    defer organic_sorption_state.deinit();
    var litter_colonization_state = try litter_colonization.State.init(std.testing.allocator, 1);
    defer litter_colonization_state.deinit();
    organic_sorption_state.doc_sorption_g_c[0] = 0.02;
    organic_sorption_state.acetate_sorption_g_c[0] = 0.01;
    organic_state.structural[3].carbon_g_c = 1;
    organic_state.colonized_structural_carbon_g_c[3] = 1;
    organic_state.residue[0].carbon_g_c = 0.2;
    organic_state.adsorbed[0].carbon_g_c = 0.1;
    organic_state.adsorbed_acetate_carbon_g_c[0] = 0.1;
    organic_decomposition_state.structural_decomposition[3].carbon_g_c = 0.1;
    organic_decomposition_state.particulate_products[3].carbon_g_c = 0.02;
    organic_decomposition_state.dissolved_structural_products[3].carbon_g_c = 0.08;
    organic_decomposition_state.microbial_residue_decomposition[0].carbon_g_c = 0.02;
    organic_decomposition_state.sorbed_organic_decomposition[0].carbon_g_c = 0.01;
    organic_decomposition_state.sorbed_acetate_decomposition_g_c[0] = 0.01;
    priming_state.exchange.dissolved_change[0].carbon_g_c = -0.01;
    priming_state.exchange.dissolved_change[1].carbon_g_c = 0.01;
    priming_state.exchange.microbial_change[0].carbon_g_c = -0.01;
    priming_state.exchange.microbial_change[respiration.source_population_count * organic.kinetic_fraction_count].carbon_g_c = 0.01;
    organic_state.microbial[0].carbon_g_c = 1;
    turnover_state.basal[0] = .{
        .decomposed = .{ .carbon_g_c = 0.1, .nitrogen_g_n = 0, .phosphorus_g_p = 0 },
        .recycled = .{ .carbon_g_c = 0.03, .nitrogen_g_n = 0, .phosphorus_g_p = 0 },
        .humified = .{ .carbon_g_c = 0.014, .nitrogen_g_n = 0, .phosphorus_g_p = 0 },
        .microbial_residue = .{ .carbon_g_c = 0.056, .nitrogen_g_n = 0, .phosphorus_g_p = 0 },
    };
    const runtime_config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var model_grid = try grid.GridState.init(std.testing.allocator, runtime_config);
    defer model_grid.deinit();
    model_grid.matrix_liquid_water_m3[0] = 1;
    var topsoil_chemistry = try soil_chemistry.State.init(std.testing.allocator, 1);
    defer topsoil_chemistry.deinit();
    var topsoil_organic = try organic.State.init(std.testing.allocator, 1);
    defer topsoil_organic.deinit();
    var parameters: respiration.Parameters = undefined;
    parameters.populations = [_]@import("soil_microbial_respiration_activity.zig").PopulationParameters{.{ .metabolism = .aerobic_heterotroph, .substrate_unlimited_respiration_per_h = 0 }} ** respiration.source_population_count;
    const carbon_before = try organic_state.totalCarbon_g_c(0) + try topsoil_organic.totalCarbon_g_c(0) + gas_state.dissolved_mass_g[@intFromEnum(gas.Species.carbon_dioxide)] + gas_state.dissolved_mass_g[@intFromEnum(gas.Species.methane)];
    const nitrogen_before = organic_state.dissolved[0].nitrogen_g_n + chemistry_state.cells[0].nitrate_mol_per_m3 * 14 + gas_state.dissolved_mass_g[@intFromEnum(gas.Species.nitrogen)] + gas_state.dissolved_mass_g[@intFromEnum(gas.Species.nitrous_oxide)];
    var context: ApplyContext = .{ .surface_organic = &organic_state, .litter_chemistry = &chemistry_state, .litter_gas = &gas_state, .litter_water_m3 = &.{1}, .respiration = &respiration_state, .oxygen = &oxygen_state, .nitrogen_fixation = &fixation_state, .substrate_uptake = &uptake_state, .denitrification = &denitrification_state, .assimilation = &assimilation_state, .mineral_exchange = &mineral_exchange_state, .topsoil_exchange = &topsoil_exchange_state, .turnover = &turnover_state, .priming = &priming_state, .organic_decomposition = &organic_decomposition_state, .organic_sorption = &organic_sorption_state, .litter_colonization = &litter_colonization_state, .topsoil_organic = &topsoil_organic, .topsoil_humus_partition = &.{.{ 0.5, 0.5 }}, .topsoil_chemistry = &topsoil_chemistry, .model_grid = &model_grid, .zone_fractions = .{ .ammonium_non_band = 1, .ammonium_band = 0, .nitrate_non_band = 1, .nitrate_band = 0, .phosphate_non_band = 1, .phosphate_band = 0 }, .microbial_parameters = parameters, .nitrogen_molar_mass_g_per_mol = 14, .phosphorus_molar_mass_g_per_mol = 31, .tolerance = 1e-12 };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    const carbon_after = try organic_state.totalCarbon_g_c(0) + try topsoil_organic.totalCarbon_g_c(0) + gas_state.dissolved_mass_g[@intFromEnum(gas.Species.carbon_dioxide)] + gas_state.dissolved_mass_g[@intFromEnum(gas.Species.methane)];
    const nitrogen_after = organic_state.dissolved[0].nitrogen_g_n + organic_state.microbial[2].nitrogen_g_n + chemistry_state.cells[0].nitrate_mol_per_m3 * 14 + gas_state.dissolved_mass_g[@intFromEnum(gas.Species.nitrogen)] + gas_state.dissolved_mass_g[@intFromEnum(gas.Species.nitrous_oxide)];
    try std.testing.expectApproxEqAbs(carbon_before, carbon_after, 1e-12);
    try std.testing.expectApproxEqAbs(nitrogen_before, nitrogen_after, 1e-12);
}

test "insufficient substrate leaves every surface metabolism pool unchanged" {
    var organic_state = try organic.State.init(std.testing.allocator, 1);
    defer organic_state.deinit();
    organic_state.dissolved[0].carbon_g_c = 0.1;
    var chemistry_state = try chemistry.State.init(std.testing.allocator, 1);
    defer chemistry_state.deinit();
    var gas_state = try gas.State.init(std.testing.allocator, 1);
    defer gas_state.deinit();
    var respiration_state = try respiration.State.init(std.testing.allocator, 1);
    defer respiration_state.deinit();
    var oxygen_state = try oxygen.State.init(std.testing.allocator, 1);
    defer oxygen_state.deinit();
    var fixation_state = try fixation.State.init(std.testing.allocator, 1);
    defer fixation_state.deinit();
    var uptake_state = try uptake.State.init(std.testing.allocator, 1);
    defer uptake_state.deinit();
    uptake_state.doc_uptake_g_c[0] = 0.2;
    var denitrification_state = try denitrification.State.init(std.testing.allocator, 1);
    defer denitrification_state.deinit();
    var assimilation_state = try assimilation.State.init(std.testing.allocator, 1);
    defer assimilation_state.deinit();
    var mineral_exchange_state = try mineral_exchange.State.init(std.testing.allocator, 1);
    defer mineral_exchange_state.deinit();
    var topsoil_exchange_state = try topsoil_exchange.State.init(std.testing.allocator, 1);
    defer topsoil_exchange_state.deinit();
    var turnover_state = try turnover.State.init(std.testing.allocator, 1);
    defer turnover_state.deinit();
    var priming_state = try priming.State.init(std.testing.allocator, 1);
    defer priming_state.deinit();
    var organic_decomposition_state = try organic_decomposition.State.init(std.testing.allocator, 1);
    defer organic_decomposition_state.deinit();
    var organic_sorption_state = try organic_sorption.State.init(std.testing.allocator, 1);
    defer organic_sorption_state.deinit();
    var litter_colonization_state = try litter_colonization.State.init(std.testing.allocator, 1);
    defer litter_colonization_state.deinit();
    const runtime_config = try @import("config.zig").SimulationConfig.init(.{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1 }, .{ .worker_threads = 1, .tile_cells = 1 }, .{ .relative_tolerance = 1e-8, .absolute_tolerance = 1e-11, .max_nonlinear_iterations = 20 });
    var model_grid = try grid.GridState.init(std.testing.allocator, runtime_config);
    defer model_grid.deinit();
    model_grid.matrix_liquid_water_m3[0] = 1;
    var topsoil_chemistry = try soil_chemistry.State.init(std.testing.allocator, 1);
    defer topsoil_chemistry.deinit();
    var topsoil_organic = try organic.State.init(std.testing.allocator, 1);
    defer topsoil_organic.deinit();
    var parameters: respiration.Parameters = undefined;
    parameters.populations = [_]@import("soil_microbial_respiration_activity.zig").PopulationParameters{.{ .metabolism = .aerobic_heterotroph, .substrate_unlimited_respiration_per_h = 0 }} ** respiration.source_population_count;
    const before_doc = organic_state.dissolved[0];
    const before_gas = gas_state.dissolved_mass_g;
    var gas_copy: [gas.species_count]f64 = undefined;
    @memcpy(&gas_copy, before_gas[0..gas.species_count]);
    var context: ApplyContext = .{ .surface_organic = &organic_state, .litter_chemistry = &chemistry_state, .litter_gas = &gas_state, .litter_water_m3 = &.{1}, .respiration = &respiration_state, .oxygen = &oxygen_state, .nitrogen_fixation = &fixation_state, .substrate_uptake = &uptake_state, .denitrification = &denitrification_state, .assimilation = &assimilation_state, .mineral_exchange = &mineral_exchange_state, .topsoil_exchange = &topsoil_exchange_state, .turnover = &turnover_state, .priming = &priming_state, .organic_decomposition = &organic_decomposition_state, .organic_sorption = &organic_sorption_state, .litter_colonization = &litter_colonization_state, .topsoil_organic = &topsoil_organic, .topsoil_humus_partition = &.{.{ 0.5, 0.5 }}, .topsoil_chemistry = &topsoil_chemistry, .model_grid = &model_grid, .zone_fractions = .{ .ammonium_non_band = 1, .ammonium_band = 0, .nitrate_non_band = 1, .nitrate_band = 0, .phosphate_non_band = 1, .phosphate_band = 0 }, .microbial_parameters = parameters, .nitrogen_molar_mass_g_per_mol = 14, .phosphorus_molar_mass_g_per_mol = 31, .tolerance = 1e-12 };
    try std.testing.expectError(error.InsufficientSurfaceMicrobialSubstrate, applyTile(&context, .{ .first = 0, .end = 1 }));
    try std.testing.expectEqual(before_doc, organic_state.dissolved[0]);
    try std.testing.expectEqualSlices(f64, &gas_copy, gas_state.dissolved_mass_g[0..gas.species_count]);
}
