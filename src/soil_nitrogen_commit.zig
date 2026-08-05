const std = @import("std");
const compute = @import("compute.zig");
const gas = @import("gas_transport.zig");
const organic = @import("soil_organic_initialization.zig");
const microbial = @import("soil_microbial_state.zig");
const chemistry = @import("solute_chemistry_state.zig");
const zones = @import("solute_charge_classification.zig");
const reactive = @import("soil_reactive_nitrogen_state.zig");
const phosphorus = @import("soil_microbial_phosphorus_state.zig");
const fluxes = @import("soil_nitrogen_flux_workspace.zig");
const turnover = @import("soil_microbial_turnover_step.zig");
const colonization = @import("soil_litter_colonization_step.zig");
const organic_sorption = @import("soil_organic_sorption_step.zig");
const organic_decomposition = @import("soil_organic_decomposition_step.zig");
const organic_priming = @import("soil_organic_priming_step.zig");
const respiration_products = @import("soil_respiration_products_step.zig");
const methane_step = @import("soil_methane_step.zig");
const transformation_aggregation = @import("soil_transformation_aggregation.zig");
const biochemical_acidity = @import("soil_biochemical_acidity.zig");

pub const ApplyContext = struct {
    reactive_nitrogen: *reactive.State,
    phosphorus_history: *phosphorus.State,
    chemistry_state: *chemistry.State,
    gas_state: *gas.State,
    organic_state: *organic.State,
    microbial_state: *microbial.State,
    flux_workspace: *const fluxes.State,
    microbial_turnover: *const turnover.State,
    litter_colonization: *const colonization.State,
    organic_sorption: *const organic_sorption.State,
    organic_decomposition: *const organic_decomposition.State,
    organic_priming: ?*const organic_priming.State = null,
    respiration_products: ?*const respiration_products.State = null,
    methane: ?*const methane_step.State = null,
    methane_autotrophic_substrate_index: usize = std.math.maxInt(usize),
    methanotroph_population_index: usize = 0,
    humus_partition_by_cell: []const [2]f64,
    soil_layer_capacity: usize,
    water_volume_m3: []const f64,
    zone_fractions: zones.ZoneFractions,
    oxygen_satisfaction_fraction: []const f64,
    redox_satisfaction_fraction: []const f64,
    nitrogen_molar_mass_g_per_mol: f64,
    phosphorus_molar_mass_g_per_mol: f64,
    tolerance_g_n: f64,
    dynamic_salts: bool = false,
    timestep_h: f64 = 1,
    hourly_signed_heterotrophic_respiration_g_c: ?[]f64 = null,
    hourly_carbon_dioxide_production_g_c: ?[]f64 = null,
};

pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    try validate(context.*, range);
    for (range.first..range.end) |layer| try commitLayer(context, layer);
}

fn commitLayer(context: *ApplyContext, layer: usize) !void {
    try context.reactive_nitrogen.validateLayer(layer);
    const water_m3 = context.water_volume_m3[layer];
    const n_mass = context.nitrogen_molar_mass_g_per_mol;
    var ammonium_non_band = context.chemistry_state.aqueous[layer].ammonium_non_band * water_m3 * context.zone_fractions.ammonium_non_band * n_mass;
    var ammonium_band = context.chemistry_state.aqueous[layer].ammonium_band * water_m3 * context.zone_fractions.ammonium_band * n_mass;
    if (ammonium_non_band > 1e8 or ammonium_band > 1e8) std.log.warn(
        "large ammonium entering commit: layer={d} non_band_g={e} band_g={e} non_band_conc={e} water_m3={e}",
        .{ layer, ammonium_non_band, ammonium_band, context.chemistry_state.aqueous[layer].ammonium_non_band, water_m3 },
    );
    var nitrate_non_band = context.chemistry_state.aqueous[layer].nitrate_non_band * water_m3 * context.zone_fractions.nitrate_non_band * n_mass;
    var nitrate_band = context.chemistry_state.aqueous[layer].nitrate_band * water_m3 * context.zone_fractions.nitrate_band * n_mass;
    const p_mass = context.phosphorus_molar_mass_g_per_mol;
    var h2po4 = [2]f64{ context.chemistry_state.non_band_phosphate[layer].dissolved_h2po4_mol_p_per_m3 * water_m3 * context.zone_fractions.phosphate_non_band * p_mass, context.chemistry_state.band_phosphate[layer].dissolved_h2po4_mol_p_per_m3 * water_m3 * context.zone_fractions.phosphate_band * p_mass };
    var hpo4 = [2]f64{ context.chemistry_state.non_band_phosphate[layer].dissolved_hpo4_mol_p_per_m3 * water_m3 * context.zone_fractions.phosphate_non_band * p_mass, context.chemistry_state.band_phosphate[layer].dissolved_hpo4_mol_p_per_m3 * water_m3 * context.zone_fractions.phosphate_band * p_mass };
    var nitrite_non_band = context.reactive_nitrogen.non_band_nitrite_g_n[layer];
    var nitrite_band = context.reactive_nitrogen.band_nitrite_g_n[layer];
    const n2o_index = try gas.massIndex(layer, .nitrous_oxide, context.gas_state.cell_count);
    const n2_index = try gas.massIndex(layer, .nitrogen, context.gas_state.cell_count);
    const co2_index = try gas.massIndex(layer, .carbon_dioxide, context.gas_state.cell_count);
    const methane_index = try gas.massIndex(layer, .methane, context.gas_state.cell_count);
    const hydrogen_index = try gas.massIndex(layer, .hydrogen, context.gas_state.cell_count);
    const oxygen_index = try gas.massIndex(layer, .oxygen, context.gas_state.cell_count);
    var nitrous_oxide = context.gas_state.dissolved_mass_g[n2o_index];
    var dinitrogen = context.gas_state.dissolved_mass_g[n2_index];
    var carbon_dioxide = context.gas_state.dissolved_mass_g[co2_index];
    var methane = context.gas_state.dissolved_mass_g[methane_index];
    const carbon_dioxide_before = carbon_dioxide;
    const methane_before = methane;
    var hydrogen = context.gas_state.dissolved_mass_g[hydrogen_index];
    var aqueous_oxygen = context.gas_state.dissolved_mass_g[oxygen_index];
    var gaseous_oxygen = context.gas_state.gaseous_mass_g[oxygen_index];
    const conserved_carbon_before_g_c = try authoritativeLayerCarbon_g_c(context.*, layer);
    const conserved_nitrogen_before_g_n = try authoritativeLayerNitrogen_g_n(context.*, layer);
    const first = layer * context.flux_workspace.process_unit_count_per_layer;
    const end = first + context.flux_workspace.process_unit_count_per_layer;
    const next_nonstructural = try context.microbial_state.allocator.alloc(microbial.ElementalPool, end - first);
    defer context.microbial_state.allocator.free(next_nonstructural);

    const published = try transformation_aggregation.calculateSupplyLimited(.{
        .complex_count = context.microbial_state.substrate_count,
        .population_count = context.microbial_state.population_count,
        .oxygen_satisfaction_fraction = context.oxygen_satisfaction_fraction[first..end],
        .redox_satisfaction_fraction = context.redox_satisfaction_fraction[first..end],
        .non_band_ammonia_oxidation_potential_g_n = context.flux_workspace.non_band_ammonia_oxidation_potential_g_n[first..end],
        .band_ammonia_oxidation_potential_g_n = context.flux_workspace.band_ammonia_oxidation_potential_g_n[first..end],
        .non_band_nitrite_oxidation_potential_g_n = context.flux_workspace.non_band_nitrite_oxidation_potential_g_n[first..end],
        .band_nitrite_oxidation_potential_g_n = context.flux_workspace.band_nitrite_oxidation_potential_g_n[first..end],
        .non_band_nitrate_reduction_potential_g_n = context.flux_workspace.non_band_nitrate_reduction_potential_g_n[first..end],
        .band_nitrate_reduction_potential_g_n = context.flux_workspace.band_nitrate_reduction_potential_g_n[first..end],
        .non_band_heterotrophic_nitrite_reduction_potential_g_n = context.flux_workspace.non_band_heterotrophic_nitrite_reduction_potential_g_n[first..end],
        .band_heterotrophic_nitrite_reduction_potential_g_n = context.flux_workspace.band_heterotrophic_nitrite_reduction_potential_g_n[first..end],
        .non_band_autotrophic_nitrite_reduction_potential_g_n = context.flux_workspace.non_band_autotrophic_nitrite_reduction_potential_g_n[first..end],
        .band_autotrophic_nitrite_reduction_potential_g_n = context.flux_workspace.band_autotrophic_nitrite_reduction_potential_g_n[first..end],
        .non_band_autotrophic_ammonium_oxidation_potential_g_n = context.flux_workspace.non_band_autotrophic_ammonium_oxidation_potential_g_n[first..end],
        .band_autotrophic_ammonium_oxidation_potential_g_n = context.flux_workspace.band_autotrophic_ammonium_oxidation_potential_g_n[first..end],
        .nitrous_oxide_reduction_potential_g_n = context.flux_workspace.nitrous_oxide_reduction_potential_g_n[first..end],
        .non_band_ammonium_exchange_g_n = context.flux_workspace.non_band_microbial_ammonium_exchange_g_n[first..end],
        .band_ammonium_exchange_g_n = context.flux_workspace.band_microbial_ammonium_exchange_g_n[first..end],
        .non_band_nitrate_exchange_g_n = context.flux_workspace.non_band_microbial_nitrate_exchange_g_n[first..end],
        .band_nitrate_exchange_g_n = context.flux_workspace.band_microbial_nitrate_exchange_g_n[first..end],
        .non_band_h2po4_exchange_g_p = context.flux_workspace.non_band_microbial_h2po4_exchange_g_p[first..end],
        .band_h2po4_exchange_g_p = context.flux_workspace.band_microbial_h2po4_exchange_g_p[first..end],
        .non_band_hpo4_exchange_g_p = context.flux_workspace.non_band_microbial_hpo4_exchange_g_p[first..end],
        .band_hpo4_exchange_g_p = context.flux_workspace.band_microbial_hpo4_exchange_g_p[first..end],
        .fixed_dinitrogen_g_n = context.flux_workspace.fixed_dinitrogen_g_n[first..end],
    });
    const ammonia_oxidation = published.ammonia_oxidation_g_n;
    const nitrite_oxidation = published.nitrite_oxidation_g_n;
    const nitrate_reduction = published.nitrate_reduction_g_n;
    const heterotrophic_nitrite_reduction = published.heterotrophic_nitrite_reduction_g_n;
    const autotrophic_nitrite_reduction = published.autotrophic_nitrite_reduction_g_n;
    const autotrophic_ammonium_oxidation = published.autotrophic_ammonium_oxidation_g_n;
    const nitrous_oxide_reduction = published.nitrous_oxide_reduction_g_n;
    const microbial_ammonium_exchange = published.ammonium_exchange_g_n;
    const microbial_nitrate_exchange = published.nitrate_exchange_g_n;
    const microbial_h2po4_exchange = published.h2po4_exchange_g_p;
    const microbial_hpo4_exchange = published.hpo4_exchange_g_p;
    const acidity = if (context.dynamic_salts) try biochemical_acidity.stage(
        context.chemistry_state,
        layer,
        .{
            .water_volume_m3 = water_m3,
            .ammonia_oxidation_g_n = published.ammonia_oxidation_g_n,
            .nitrate_reduction_g_n = published.nitrate_reduction_g_n,
            .nitrite_reduction_g_n = published.heterotrophic_nitrite_reduction_g_n,
            .nitrous_oxide_reduction_g_n = published.nitrous_oxide_reduction_g_n,
            .lignin_decomposition_g_c = ligninDecompositionForAcidity(
                context.organic_decomposition,
                layer,
            ),
            .timestep_h = context.timestep_h,
            .negligible_hydrogen_mol = context.tolerance_g_n,
        },
    ) else null;
    var respiratory_carbon_dioxide_g_c: f64 = 0;
    for (first..end) |unit| {
        const redox = context.redox_satisfaction_fraction[unit];
        respiratory_carbon_dioxide_g_c += (if (context.respiration_products) |products| products.carbon_dioxide_g_c[unit] else context.flux_workspace.actual_aerobic_respiration_g_c[unit]) + context.flux_workspace.denitrification_respiration_g_c[unit] * redox + context.flux_workspace.nitrogen_fixation_respiration_g_c[unit];
        if (context.respiration_products) |products| {
            if (context.methane == null) methane += products.methane_g_c[unit];
            hydrogen += products.hydrogen_g_h[unit];
        }
        for (0..2) |component| respiratory_carbon_dioxide_g_c += context.microbial_turnover.senescence[unit * 2 + component].recycled.carbon_g_c;
    }
    const chemo_nitrite = [2]f64{ context.flux_workspace.chemodenitrification_non_band_nitrite_reduction_g_n[layer], context.flux_workspace.chemodenitrification_band_nitrite_reduction_g_n[layer] };
    const chemo_n2o = context.flux_workspace.chemodenitrification_nitrous_oxide_production_g_n[layer];
    const chemo_n2 = context.flux_workspace.chemodenitrification_dinitrogen_production_g_n[layer];
    const chemo_don = context.flux_workspace.chemodenitrification_dissolved_organic_nitrogen_production_g_n[layer];
    const chemo_consumed = chemo_nitrite[0] + chemo_nitrite[1];
    if (@abs(chemo_consumed - chemo_n2o - chemo_n2 - chemo_don) > context.tolerance_g_n * @max(1, chemo_consumed)) return error.SoilNitrogenBalanceFailure;

    try applyZone(&ammonium_non_band, &nitrate_non_band, &nitrite_non_band, ammonia_oxidation[0], nitrite_oxidation[0], nitrate_reduction[0], heterotrophic_nitrite_reduction[0], autotrophic_nitrite_reduction[0], autotrophic_ammonium_oxidation[0], chemo_nitrite[0], context.tolerance_g_n);
    try applyZone(&ammonium_band, &nitrate_band, &nitrite_band, ammonia_oxidation[1], nitrite_oxidation[1], nitrate_reduction[1], heterotrophic_nitrite_reduction[1], autotrophic_nitrite_reduction[1], autotrophic_ammonium_oxidation[1], chemo_nitrite[1], context.tolerance_g_n);
    try applyMicrobialExchange(&ammonium_non_band, microbial_ammonium_exchange[0], context.tolerance_g_n);
    try applyMicrobialExchange(&ammonium_band, microbial_ammonium_exchange[1], context.tolerance_g_n);
    try applyMicrobialExchange(&nitrate_non_band, microbial_nitrate_exchange[0], context.tolerance_g_n);
    try applyMicrobialExchange(&nitrate_band, microbial_nitrate_exchange[1], context.tolerance_g_n);
    try applyMicrobialExchange(&h2po4[0], microbial_h2po4_exchange[0], context.tolerance_g_n);
    try applyMicrobialExchange(&h2po4[1], microbial_h2po4_exchange[1], context.tolerance_g_n);
    try applyMicrobialExchange(&hpo4[0], microbial_hpo4_exchange[0], context.tolerance_g_n);
    try applyMicrobialExchange(&hpo4[1], microbial_hpo4_exchange[1], context.tolerance_g_n);
    const biological_n2o_input = heterotrophic_nitrite_reduction[0] + heterotrophic_nitrite_reduction[1] + autotrophic_nitrite_reduction[0] + autotrophic_nitrite_reduction[1];
    if (nitrous_oxide_reduction > nitrous_oxide + biological_n2o_input + chemo_n2o + context.tolerance_g_n) return error.InsufficientSoilNitrousOxide;
    nitrous_oxide = @max(0, nitrous_oxide + biological_n2o_input + chemo_n2o - nitrous_oxide_reduction);
    if (published.fixed_dinitrogen_g_n > dinitrogen + nitrous_oxide_reduction + chemo_n2 + context.tolerance_g_n) return error.InsufficientSoilDinitrogen;
    dinitrogen = @max(0, dinitrogen + nitrous_oxide_reduction + chemo_n2 - published.fixed_dinitrogen_g_n);
    carbon_dioxide += respiratory_carbon_dioxide_g_c;
    if (context.methane) |methane_result| {
        methane = methane_result.aqueous_methane_after_g_c[layer];
        hydrogen -= methane_result.hydrogen_consumption_g_h[layer];
        carbon_dioxide += methane_result.methane_oxidation_respiration_g_c[layer] - methane_result.hydrogenotrophic_methane_g_c[layer];
        if (hydrogen < -context.tolerance_g_n or carbon_dioxide < -context.tolerance_g_n) return error.InsufficientSoilMethanogenesisSubstrate;
        hydrogen = @max(0, hydrogen);
        carbon_dioxide = @max(0, carbon_dioxide);
        const oxygen_demand = methane_result.oxygen_demand_g_o[layer];
        const aqueous_use = @min(aqueous_oxygen, oxygen_demand);
        aqueous_oxygen -= aqueous_use;
        gaseous_oxygen -= oxygen_demand - aqueous_use;
        if (gaseous_oxygen < -context.tolerance_g_n) return error.InsufficientSoilMethanotrophOxygen;
        gaseous_oxygen = @max(0, gaseous_oxygen);
    }

    var residue_total_g_c: f64 = 0;
    var residue_by_complex_g_c: [organic.substrate_count]f64 = .{0} ** organic.substrate_count;
    var dissolved_after: [organic.substrate_count]organic.ElementPool = undefined;
    var dissolved_acetate_after_g_c: [organic.substrate_count]f64 = undefined;
    var adsorbed_after: [organic.substrate_count]organic.ElementPool = undefined;
    var adsorbed_acetate_after_g_c: [organic.substrate_count]f64 = undefined;
    for (0..organic.substrate_count) |complex| for (0..organic.residue_fraction_count) |fraction| {
        const carbon = context.organic_state.residue[(layer * organic.substrate_count + complex) * organic.residue_fraction_count + fraction].carbon_g_c;
        residue_by_complex_g_c[complex] += carbon;
        residue_total_g_c += carbon;
    };
    for (0..organic.substrate_count) |complex| {
        const fraction: f64 = if (residue_total_g_c > context.tolerance_g_n) residue_by_complex_g_c[complex] / residue_total_g_c else if (complex == 3) 1 else 0;
        const mobile = layer * organic.substrate_count + complex;
        var doc_uptake: f64 = 0;
        var don_uptake: f64 = 0;
        var dop_uptake: f64 = 0;
        var acetate_uptake: f64 = 0;
        var fermentation_acetate_g_c: f64 = 0;
        if (complex < context.microbial_state.substrate_count) for (0..context.microbial_state.population_count) |population| {
            const unit = first + complex * context.microbial_state.population_count + population;
            doc_uptake += context.flux_workspace.doc_uptake_g_c[unit];
            don_uptake += context.flux_workspace.dissolved_organic_nitrogen_uptake_g_n[unit];
            dop_uptake += context.flux_workspace.dissolved_organic_phosphorus_uptake_g_p[unit];
            acetate_uptake += context.flux_workspace.acetate_uptake_g_c[unit];
            if (context.respiration_products) |products| fermentation_acetate_g_c += products.acetate_g_c[unit];
        };
        const priming_dissolved: organic.ElementPool = if (context.organic_priming) |priming| priming.exchange.dissolved_change[mobile] else .{};
        const priming_acetate_g_c: f64 = if (context.organic_priming) |priming| priming.exchange.acetate_change_g_c[mobile] else 0;
        const source_dissolved = context.organic_state.dissolved[mobile];
        const current: organic.ElementPool = .{ .carbon_g_c = source_dissolved.carbon_g_c + priming_dissolved.carbon_g_c, .nitrogen_g_n = source_dissolved.nitrogen_g_n + priming_dissolved.nitrogen_g_n, .phosphorus_g_p = source_dissolved.phosphorus_g_p + priming_dissolved.phosphorus_g_p };
        var sorbed = context.organic_sorption.exchange[mobile];
        var decomposition_products: organic.ElementPool = .{};
        for (0..organic.structural_fraction_count) |structural_fraction| addPool(&decomposition_products, context.organic_decomposition.dissolved_structural_products[mobile * organic.structural_fraction_count + structural_fraction]);
        for (0..organic.residue_fraction_count) |residue_fraction| addPool(&decomposition_products, context.organic_decomposition.microbial_residue_decomposition[mobile * organic.residue_fraction_count + residue_fraction]);
        addPool(&decomposition_products, context.organic_decomposition.sorbed_organic_decomposition[mobile]);
        const current_adsorbed = context.organic_state.adsorbed[mobile];
        const decomposed_sorbed = context.organic_decomposition.sorbed_organic_decomposition[mobile];
        const dissolved_before_sorption: organic.ElementPool = .{ .carbon_g_c = current.carbon_g_c + decomposition_products.carbon_g_c - doc_uptake, .nitrogen_g_n = current.nitrogen_g_n + decomposition_products.nitrogen_g_n - don_uptake + chemo_don * fraction, .phosphorus_g_p = current.phosphorus_g_p + decomposition_products.phosphorus_g_p - dop_uptake };
        const adsorbed_before_sorption: organic.ElementPool = .{ .carbon_g_c = current_adsorbed.carbon_g_c - decomposed_sorbed.carbon_g_c, .nitrogen_g_n = current_adsorbed.nitrogen_g_n - decomposed_sorbed.nitrogen_g_n, .phosphorus_g_p = current_adsorbed.phosphorus_g_p - decomposed_sorbed.phosphorus_g_p };
        sorbed.doc_g_c = boundedExchange(sorbed.doc_g_c, dissolved_before_sorption.carbon_g_c, adsorbed_before_sorption.carbon_g_c);
        sorbed.don_g_n = boundedExchange(sorbed.don_g_n, dissolved_before_sorption.nitrogen_g_n, adsorbed_before_sorption.nitrogen_g_n);
        sorbed.dop_g_p = boundedExchange(sorbed.dop_g_p, dissolved_before_sorption.phosphorus_g_p, adsorbed_before_sorption.phosphorus_g_p);
        dissolved_after[complex] = .{ .carbon_g_c = dissolved_before_sorption.carbon_g_c - sorbed.doc_g_c, .nitrogen_g_n = dissolved_before_sorption.nitrogen_g_n - sorbed.don_g_n, .phosphorus_g_p = dissolved_before_sorption.phosphorus_g_p - sorbed.dop_g_p };
        const dissolved_acetate_before_sorption = context.organic_state.dissolved_acetate_carbon_g_c[mobile] + priming_acetate_g_c + context.organic_decomposition.sorbed_acetate_decomposition_g_c[mobile] + fermentation_acetate_g_c - acetate_uptake;
        const adsorbed_acetate_before_sorption = context.organic_state.adsorbed_acetate_carbon_g_c[mobile] - context.organic_decomposition.sorbed_acetate_decomposition_g_c[mobile];
        sorbed.acetate_g_c = boundedExchange(sorbed.acetate_g_c, dissolved_acetate_before_sorption, adsorbed_acetate_before_sorption);
        dissolved_acetate_after_g_c[complex] = dissolved_acetate_before_sorption - sorbed.acetate_g_c;
        adsorbed_after[complex] = .{ .carbon_g_c = current_adsorbed.carbon_g_c - decomposed_sorbed.carbon_g_c + sorbed.doc_g_c, .nitrogen_g_n = current_adsorbed.nitrogen_g_n - decomposed_sorbed.nitrogen_g_n + sorbed.don_g_n, .phosphorus_g_p = current_adsorbed.phosphorus_g_p - decomposed_sorbed.phosphorus_g_p + sorbed.dop_g_p };
        adsorbed_acetate_after_g_c[complex] = adsorbed_acetate_before_sorption + sorbed.acetate_g_c;
        inline for (@typeInfo(organic.ElementPool).@"struct".fields) |field| {
            const source_value = @field(source_dissolved, field.name);
            const priming_value = @field(priming_dissolved, field.name);
            const decomposition_value =
                @field(decomposition_products, field.name);
            const uptake_value = switch (field.name[0]) {
                'c' => doc_uptake,
                'n' => don_uptake - chemo_don * fraction,
                'p' => dop_uptake,
                else => unreachable,
            };
            const sorption_value = switch (field.name[0]) {
                'c' => sorbed.doc_g_c,
                'n' => sorbed.don_g_n,
                'p' => sorbed.dop_g_p,
                else => unreachable,
            };
            @field(dissolved_after[complex], field.name) =
                normalizeNonnegativeRoundoff(
                    @field(dissolved_after[complex], field.name),
                    @abs(source_value) + @abs(priming_value) +
                        @abs(decomposition_value) + @abs(uptake_value) +
                        @abs(sorption_value),
                );
            const value = @field(dissolved_after[complex], field.name);
            if (!std.math.isFinite(value) or value < 0) {
                std.log.warn(
                    "invalid dissolved organic commit: layer={d} complex={d} pool={s} value={e} source={e} priming={e} decomposition={e} uptake={e} proposed_sorption={e} bounded_sorption={e}",
                    .{ layer, complex, field.name, value, source_value, priming_value, decomposition_value, uptake_value, switch (field.name[0]) {
                        'c' => context.organic_sorption.exchange[mobile].doc_g_c,
                        'n' => context.organic_sorption.exchange[mobile].don_g_n,
                        'p' => context.organic_sorption.exchange[mobile].dop_g_p,
                        else => unreachable,
                    }, sorption_value },
                );
                return error.InvalidSoilNitrogenCommit;
            }
        }
        dissolved_acetate_after_g_c[complex] =
            normalizeNonnegativeRoundoff(
                dissolved_acetate_after_g_c[complex],
                @abs(dissolved_acetate_before_sorption) +
                    @abs(sorbed.acetate_g_c),
            );
        if (!std.math.isFinite(dissolved_acetate_after_g_c[complex]) or dissolved_acetate_after_g_c[complex] < 0) {
            std.log.warn("invalid dissolved acetate commit: layer={d} complex={d} value_g_c={e}", .{ layer, complex, dissolved_acetate_after_g_c[complex] });
            return error.InvalidSoilNitrogenCommit;
        }
        inline for (@typeInfo(organic.ElementPool).@"struct".fields) |field| {
            @field(adsorbed_after[complex], field.name) =
                normalizeNonnegativeRoundoff(
                    @field(adsorbed_after[complex], field.name),
                    @abs(@field(current_adsorbed, field.name)) +
                        @abs(@field(decomposed_sorbed, field.name)) +
                        @abs(switch (field.name[0]) {
                            'c' => sorbed.doc_g_c,
                            'n' => sorbed.don_g_n,
                            'p' => sorbed.dop_g_p,
                            else => unreachable,
                        }),
                );
            const value = @field(adsorbed_after[complex], field.name);
            if (!std.math.isFinite(value) or value < 0) {
                std.log.warn("invalid adsorbed organic commit: layer={d} complex={d} pool={s} value={e}", .{ layer, complex, field.name, value });
                return error.InvalidSoilNitrogenCommit;
            }
        }
        adsorbed_acetate_after_g_c[complex] =
            normalizeNonnegativeRoundoff(
                adsorbed_acetate_after_g_c[complex],
                @abs(adsorbed_acetate_before_sorption) +
                    @abs(sorbed.acetate_g_c),
            );
        if (!std.math.isFinite(adsorbed_acetate_after_g_c[complex]) or adsorbed_acetate_after_g_c[complex] < 0) {
            std.log.warn("invalid adsorbed acetate commit: layer={d} complex={d} value_g_c={e}", .{ layer, complex, adsorbed_acetate_after_g_c[complex] });
            return error.InvalidSoilNitrogenCommit;
        }
    }
    var residue_after: [organic.substrate_count * organic.residue_fraction_count]organic.ElementPool = undefined;
    @memcpy(&residue_after, context.organic_state.residue[layer * organic.substrate_count * organic.residue_fraction_count ..][0 .. organic.substrate_count * organic.residue_fraction_count]);
    var heterotrophic_residue_carbon_g_c: [organic.substrate_count - 1]f64 = .{0} ** (organic.substrate_count - 1);
    var total_heterotrophic_residue_carbon_g_c: f64 = 0;
    for (0..organic.substrate_count - 1) |complex| for (0..organic.residue_fraction_count) |component| {
        const carbon = residue_after[complex * organic.residue_fraction_count + component].carbon_g_c;
        heterotrophic_residue_carbon_g_c[complex] += carbon;
        total_heterotrophic_residue_carbon_g_c += carbon;
    };
    for (0..residue_after.len) |index| subtractPool(&residue_after[index], context.organic_decomposition.microbial_residue_decomposition[layer * residue_after.len + index]);
    const structural_first = layer * organic.substrate_count * organic.structural_fraction_count;
    var structural_after: [organic.substrate_count * organic.structural_fraction_count]organic.ElementPool = undefined;
    var colonized_after: [organic.substrate_count * organic.structural_fraction_count]f64 = undefined;
    @memcpy(&structural_after, context.organic_state.structural[structural_first..][0..structural_after.len]);
    @memcpy(&colonized_after, context.organic_state.colonized_structural_carbon_g_c[structural_first..][0..colonized_after.len]);
    var particulate_products: organic.ElementPool = .{};
    for (0..structural_after.len) |offset| {
        const decomposed = context.organic_decomposition.structural_decomposition[structural_first + offset];
        subtractPool(&structural_after[offset], decomposed);
        colonized_after[offset] -= decomposed.carbon_g_c;
        addPool(&particulate_products, context.organic_decomposition.particulate_products[structural_first + offset]);
        colonized_after[offset] += context.litter_colonization.colonized_carbon_increment_g_c[structural_first + offset];
    }
    const particulate_index = 3 * organic.structural_fraction_count;
    addPool(&structural_after[particulate_index], particulate_products);
    colonized_after[particulate_index] += particulate_products.carbon_g_c;
    const humus_first = (layer * organic.substrate_count + 4) * organic.structural_fraction_count;
    const humus_local = 4 * organic.structural_fraction_count;
    var humus_after = [2]organic.ElementPool{ structural_after[humus_local], structural_after[humus_local + 1] };
    var humus_colonized_after = [2]f64{ colonized_after[humus_local], colonized_after[humus_local + 1] };
    const cell = layer / context.soil_layer_capacity;
    const humus_partition = context.humus_partition_by_cell[cell];
    for (0..context.microbial_state.substrate_count) |substrate| for (0..context.microbial_state.population_count) |population| {
        if (!microbial.nitroPopulationEnabled(substrate, population)) continue;
        const runtime_index = try context.microbial_state.populationIndex(layer / context.microbial_state.layer_count, layer % context.microbial_state.layer_count, substrate, population);
        const unit = first + substrate * context.microbial_state.population_count + population;
        const priming_base = (layer * organic.substrate_count * context.microbial_state.population_count + substrate * context.microbial_state.population_count + population) * organic.kinetic_fraction_count;
        const priming_nonstructural: organic.ElementPool = if (context.organic_priming != null and substrate < organic.substrate_count) context.organic_priming.?.exchange.microbial_change[priming_base + 2] else .{};
        const source_nonstructural = context.microbial_state.nonstructural[runtime_index];
        const current: @TypeOf(source_nonstructural) = .{ .carbon_g_c = source_nonstructural.carbon_g_c + priming_nonstructural.carbon_g_c, .nitrogen_g_n = source_nonstructural.nitrogen_g_n + priming_nonstructural.nitrogen_g_n, .phosphorus_g_p = source_nonstructural.phosphorus_g_p + priming_nonstructural.phosphorus_g_p };
        const mineral_n = context.flux_workspace.non_band_microbial_ammonium_exchange_g_n[unit] + context.flux_workspace.band_microbial_ammonium_exchange_g_n[unit] + context.flux_workspace.non_band_microbial_nitrate_exchange_g_n[unit] + context.flux_workspace.band_microbial_nitrate_exchange_g_n[unit];
        const mineral_p = context.flux_workspace.non_band_microbial_h2po4_exchange_g_p[unit] + context.flux_workspace.band_microbial_h2po4_exchange_g_p[unit] + context.flux_workspace.non_band_microbial_hpo4_exchange_g_p[unit] + context.flux_workspace.band_microbial_hpo4_exchange_g_p[unit];
        const assimilated: @TypeOf(current) = .{
            .carbon_g_c = context.flux_workspace.labile_assimilation_g_c[unit] + context.flux_workspace.resistant_assimilation_g_c[unit],
            .nitrogen_g_n = context.flux_workspace.labile_assimilation_g_n[unit] + context.flux_workspace.resistant_assimilation_g_n[unit],
            .phosphorus_g_p = context.flux_workspace.labile_assimilation_g_p[unit] + context.flux_workspace.resistant_assimilation_g_p[unit],
        };
        var recycled: @TypeOf(current) = .{ .carbon_g_c = 0, .nitrogen_g_n = 0, .phosphorus_g_p = 0 };
        var component_assimilation: [2]microbial.ElementalPool = undefined;
        var component_basal_recycled: [2]microbial.ElementalPool = undefined;
        var component_senescence_recycled: [2]microbial.ElementalPool = undefined;
        for (0..2) |component| {
            const basal = context.microbial_turnover.basal[unit * 2 + component];
            const senescence = context.microbial_turnover.senescence[unit * 2 + component];
            component_assimilation[component] = if (component == 0)
                .{ .carbon_g_c = context.flux_workspace.labile_assimilation_g_c[unit], .nitrogen_g_n = context.flux_workspace.labile_assimilation_g_n[unit], .phosphorus_g_p = context.flux_workspace.labile_assimilation_g_p[unit] }
            else
                .{ .carbon_g_c = context.flux_workspace.resistant_assimilation_g_c[unit], .nitrogen_g_n = context.flux_workspace.resistant_assimilation_g_n[unit], .phosphorus_g_p = context.flux_workspace.resistant_assimilation_g_p[unit] };
            component_basal_recycled[component] = .{ .carbon_g_c = basal.recycled.carbon_g_c, .nitrogen_g_n = basal.recycled.nitrogen_g_n, .phosphorus_g_p = basal.recycled.phosphorus_g_p };
            component_senescence_recycled[component] = .{ .carbon_g_c = senescence.recycled.carbon_g_c, .nitrogen_g_n = senescence.recycled.nitrogen_g_n, .phosphorus_g_p = senescence.recycled.phosphorus_g_p };
            recycled.carbon_g_c += basal.recycled.carbon_g_c;
            recycled.nitrogen_g_n += basal.recycled.nitrogen_g_n + senescence.recycled.nitrogen_g_n;
            recycled.phosphorus_g_p += basal.recycled.phosphorus_g_p + senescence.recycled.phosphorus_g_p;
            const residue_product: organic.ElementPool = .{ .carbon_g_c = basal.microbial_residue.carbon_g_c + senescence.microbial_residue.carbon_g_c, .nitrogen_g_n = basal.microbial_residue.nitrogen_g_n + senescence.microbial_residue.nitrogen_g_n, .phosphorus_g_p = basal.microbial_residue.phosphorus_g_p + senescence.microbial_residue.phosphorus_g_p };
            if (substrate == context.methane_autotrophic_substrate_index) {
                for (0..organic.substrate_count - 1) |complex| {
                    const fraction: f64 = if (total_heterotrophic_residue_carbon_g_c > context.tolerance_g_n) heterotrophic_residue_carbon_g_c[complex] / total_heterotrophic_residue_carbon_g_c else if (complex == 3) 1 else 0;
                    const residue_index = complex * organic.residue_fraction_count + component;
                    residue_after[residue_index].carbon_g_c += residue_product.carbon_g_c * fraction;
                    residue_after[residue_index].nitrogen_g_n += residue_product.nitrogen_g_n * fraction;
                    residue_after[residue_index].phosphorus_g_p += residue_product.phosphorus_g_p * fraction;
                }
            } else {
                const residue_index = substrate * organic.residue_fraction_count + component;
                residue_after[residue_index].carbon_g_c += residue_product.carbon_g_c;
                residue_after[residue_index].nitrogen_g_n += residue_product.nitrogen_g_n;
                residue_after[residue_index].phosphorus_g_p += residue_product.phosphorus_g_p;
            }
            const humified_c = basal.humified.carbon_g_c + senescence.humified.carbon_g_c;
            const humified_n = basal.humified.nitrogen_g_n + senescence.humified.nitrogen_g_n;
            const humified_p = basal.humified.phosphorus_g_p + senescence.humified.phosphorus_g_p;
            for (0..2) |humus_class| {
                humus_after[humus_class].carbon_g_c += humified_c * humus_partition[humus_class];
                humus_after[humus_class].nitrogen_g_n += humified_n * humus_partition[humus_class];
                humus_after[humus_class].phosphorus_g_p += humified_p * humus_partition[humus_class];
                humus_colonized_after[humus_class] += humified_c * humus_partition[humus_class];
            }
            const priming_structural: organic.ElementPool = if (context.organic_priming != null and substrate < organic.substrate_count) context.organic_priming.?.exchange.microbial_change[priming_base + component] else .{};
            const source_structural = context.microbial_state.structural[runtime_index * 2 + component];
            const structural_current: @TypeOf(source_structural) = .{ .carbon_g_c = source_structural.carbon_g_c + priming_structural.carbon_g_c, .nitrogen_g_n = source_structural.nitrogen_g_n + priming_structural.nitrogen_g_n, .phosphorus_g_p = source_structural.phosphorus_g_p + priming_structural.phosphorus_g_p };
            const structural_assimilation: @TypeOf(current) = if (component == 0) .{ .carbon_g_c = context.flux_workspace.labile_assimilation_g_c[unit], .nitrogen_g_n = context.flux_workspace.labile_assimilation_g_n[unit], .phosphorus_g_p = context.flux_workspace.labile_assimilation_g_p[unit] } else .{ .carbon_g_c = context.flux_workspace.resistant_assimilation_g_c[unit], .nitrogen_g_n = context.flux_workspace.resistant_assimilation_g_n[unit], .phosphorus_g_p = context.flux_workspace.resistant_assimilation_g_p[unit] };
            const next_structural: @TypeOf(current) = .{ .carbon_g_c = structural_current.carbon_g_c + structural_assimilation.carbon_g_c - basal.decomposed.carbon_g_c - senescence.decomposed.carbon_g_c, .nitrogen_g_n = structural_current.nitrogen_g_n + structural_assimilation.nitrogen_g_n - basal.decomposed.nitrogen_g_n - senescence.decomposed.nitrogen_g_n, .phosphorus_g_p = structural_current.phosphorus_g_p + structural_assimilation.phosphorus_g_p - basal.decomposed.phosphorus_g_p - senescence.decomposed.phosphorus_g_p };
            inline for (@typeInfo(@TypeOf(next_structural)).@"struct".fields) |field| if (!std.math.isFinite(@field(next_structural, field.name)) or @field(next_structural, field.name) < 0) return error.InvalidSoilMicrobialCommit;
        }
        const methane_biomass_gain_g_c: f64 = if (context.methane != null and substrate == context.methane_autotrophic_substrate_index and population == context.methanotroph_population_index) context.methane.?.methane_oxidation_to_biomass_g_c[layer] else 0;
        const next = try microbial.calculateNonstructuralPublication(.{
            .current = current,
            .assimilation = component_assimilation,
            .basal_recycled = component_basal_recycled,
            .senescence_recycled = component_senescence_recycled,
            .net_carbon_uptake_g_c = context.flux_workspace.nonstructural_carbon_gain_g_c[unit] + methane_biomass_gain_g_c,
            .dissolved_organic_nitrogen_uptake_g_n = context.flux_workspace.dissolved_organic_nitrogen_uptake_g_n[unit],
            .mineral_nitrogen_exchange_g_n = mineral_n,
            .fixed_nitrogen_g_n = context.flux_workspace.fixed_dinitrogen_g_n[unit],
            .dissolved_organic_phosphorus_uptake_g_p = context.flux_workspace.dissolved_organic_phosphorus_uptake_g_p[unit],
            .mineral_phosphorus_exchange_g_p = mineral_p,
        });
        next_nonstructural[unit - first] = next;
        inline for (@typeInfo(@TypeOf(next)).@"struct".fields) |field| {
            const value = @field(next, field.name);
            if (!std.math.isFinite(value) or value < 0) {
                std.log.warn(
                    "invalid microbial nonstructural commit: layer={d} substrate={d} population={d} pool={s} value={e} current={e} carbon_gain={e} methane_gain={e} assimilated={e} recycled={e} mineral_n={e} mineral_p={e}",
                    .{ layer, substrate, population, field.name, value, @field(current, field.name), context.flux_workspace.nonstructural_carbon_gain_g_c[unit], methane_biomass_gain_g_c, @field(assimilated, field.name), @field(recycled, field.name), mineral_n, mineral_p },
                );
                return error.InvalidSoilMicrobialCommit;
            }
        }
        const structural_values = [6]f64{ context.flux_workspace.labile_assimilation_g_c[unit], context.flux_workspace.labile_assimilation_g_n[unit], context.flux_workspace.labile_assimilation_g_p[unit], context.flux_workspace.resistant_assimilation_g_c[unit], context.flux_workspace.resistant_assimilation_g_n[unit], context.flux_workspace.resistant_assimilation_g_p[unit] };
        for (structural_values) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSoilMicrobialCommit;
    };
    for (residue_after) |pool| inline for (@typeInfo(organic.ElementPool).@"struct".fields) |field| if (!std.math.isFinite(@field(pool, field.name)) or @field(pool, field.name) < 0) return error.InvalidSoilMicrobialCommit;
    for (structural_after) |pool| inline for (@typeInfo(organic.ElementPool).@"struct".fields) |field| if (!std.math.isFinite(@field(pool, field.name)) or @field(pool, field.name) < 0) return error.InvalidSoilOrganicDecompositionCommit;
    for (humus_after) |pool| inline for (@typeInfo(organic.ElementPool).@"struct".fields) |field| if (!std.math.isFinite(@field(pool, field.name)) or @field(pool, field.name) < 0) return error.InvalidSoilMicrobialCommit;
    for (structural_after, colonized_after) |pool, colonized| if (!std.math.isFinite(colonized) or colonized < -context.tolerance_g_n or colonized > pool.carbon_g_c + context.tolerance_g_n) return error.InvalidSoilLitterColonizationCommit;

    const ammonium_non_band_concentration = try concentrationFromMass(ammonium_non_band, water_m3, context.zone_fractions.ammonium_non_band, n_mass, context.tolerance_g_n);
    const ammonium_band_concentration = try concentrationFromMass(ammonium_band, water_m3, context.zone_fractions.ammonium_band, n_mass, context.tolerance_g_n);
    const nitrate_non_band_concentration = try concentrationFromMass(nitrate_non_band, water_m3, context.zone_fractions.nitrate_non_band, n_mass, context.tolerance_g_n);
    const nitrate_band_concentration = try concentrationFromMass(nitrate_band, water_m3, context.zone_fractions.nitrate_band, n_mass, context.tolerance_g_n);
    const non_band_h2po4_concentration = try concentrationFromMass(h2po4[0], water_m3, context.zone_fractions.phosphate_non_band, p_mass, context.tolerance_g_n);
    const band_h2po4_concentration = try concentrationFromMass(h2po4[1], water_m3, context.zone_fractions.phosphate_band, p_mass, context.tolerance_g_n);
    const non_band_hpo4_concentration = try concentrationFromMass(hpo4[0], water_m3, context.zone_fractions.phosphate_non_band, p_mass, context.tolerance_g_n);
    const band_hpo4_concentration = try concentrationFromMass(hpo4[1], water_m3, context.zone_fractions.phosphate_band, p_mass, context.tolerance_g_n);
    try validatePublishedAmmonium(
        layer,
        context.chemistry_state.water_mol_per_m3[layer],
        .{ ammonium_non_band_concentration, ammonium_band_concentration },
        .{ ammonium_non_band, ammonium_band },
        microbial_ammonium_exchange,
        .{
            ammonia_oxidation[0] + autotrophic_ammonium_oxidation[0],
            ammonia_oxidation[1] + autotrophic_ammonium_oxidation[1],
        },
        water_m3,
        .{
            context.zone_fractions.ammonium_non_band,
            context.zone_fractions.ammonium_band,
        },
    );
    context.chemistry_state.aqueous[layer].ammonium_non_band = ammonium_non_band_concentration;
    context.chemistry_state.aqueous[layer].ammonium_band = ammonium_band_concentration;
    context.chemistry_state.aqueous[layer].nitrate_non_band = nitrate_non_band_concentration;
    context.chemistry_state.aqueous[layer].nitrate_band = nitrate_band_concentration;
    context.chemistry_state.non_band_phosphate[layer].dissolved_h2po4_mol_p_per_m3 = non_band_h2po4_concentration;
    context.chemistry_state.band_phosphate[layer].dissolved_h2po4_mol_p_per_m3 = band_h2po4_concentration;
    context.chemistry_state.non_band_phosphate[layer].dissolved_hpo4_mol_p_per_m3 = non_band_hpo4_concentration;
    context.chemistry_state.band_phosphate[layer].dissolved_hpo4_mol_p_per_m3 = band_hpo4_concentration;
    if (acidity) |staged| try biochemical_acidity.commit(
        context.chemistry_state,
        layer,
        staged,
    );
    context.reactive_nitrogen.non_band_nitrite_g_n[layer] = nitrite_non_band;
    context.reactive_nitrogen.band_nitrite_g_n[layer] = nitrite_band;
    context.reactive_nitrogen.previous_non_band_chemodenitrification_capacity_g_n[layer] = context.flux_workspace.chemodenitrification_non_band_unlimited_reduction_g_n[layer];
    context.reactive_nitrogen.previous_band_chemodenitrification_capacity_g_n[layer] = context.flux_workspace.chemodenitrification_band_unlimited_reduction_g_n[layer];
    var next_non_band_ammonium_demand: f64 = 0;
    var next_band_ammonium_demand: f64 = 0;
    var next_non_band_nitrite_demand = context.flux_workspace.chemodenitrification_non_band_unlimited_reduction_g_n[layer];
    var next_band_nitrite_demand = context.flux_workspace.chemodenitrification_band_unlimited_reduction_g_n[layer];
    var next_non_band_nitrate_demand: f64 = 0;
    var next_band_nitrate_demand: f64 = 0;
    var next_nitrous_oxide_demand: f64 = 0;
    var next_h2po4_demand = [2]f64{ 0, 0 };
    var next_hpo4_demand = [2]f64{ 0, 0 };
    for (first..end) |unit| {
        context.reactive_nitrogen.previous_aerobic_oxygen_demand_g_o[unit] = context.flux_workspace.aerobic_oxygen_demand_g_o[unit];
        context.reactive_nitrogen.previous_doc_respiration_demand_g_c[unit] = context.flux_workspace.doc_respiration_demand_g_c[unit];
        context.reactive_nitrogen.previous_acetate_respiration_demand_g_c[unit] = context.flux_workspace.acetate_respiration_demand_g_c[unit];
        context.reactive_nitrogen.previous_non_band_ammonia_oxidation_capacity_g_n[unit] = context.flux_workspace.non_band_ammonia_oxidation_capacity_g_n[unit];
        context.reactive_nitrogen.previous_band_ammonia_oxidation_capacity_g_n[unit] = context.flux_workspace.band_ammonia_oxidation_capacity_g_n[unit];
        context.reactive_nitrogen.previous_non_band_nitrite_oxidation_capacity_g_n[unit] = context.flux_workspace.non_band_nitrite_oxidation_capacity_g_n[unit];
        context.reactive_nitrogen.previous_band_nitrite_oxidation_capacity_g_n[unit] = context.flux_workspace.band_nitrite_oxidation_capacity_g_n[unit];
        context.reactive_nitrogen.previous_non_band_nitrate_reduction_capacity_g_n[unit] = context.flux_workspace.non_band_nitrate_reduction_capacity_g_n[unit];
        context.reactive_nitrogen.previous_band_nitrate_reduction_capacity_g_n[unit] = context.flux_workspace.band_nitrate_reduction_capacity_g_n[unit];
        context.reactive_nitrogen.previous_non_band_nitrite_reduction_capacity_g_n[unit] = context.flux_workspace.non_band_nitrite_reduction_capacity_g_n[unit];
        context.reactive_nitrogen.previous_band_nitrite_reduction_capacity_g_n[unit] = context.flux_workspace.band_nitrite_reduction_capacity_g_n[unit];
        context.reactive_nitrogen.previous_nitrous_oxide_reduction_capacity_g_n[unit] = context.flux_workspace.nitrous_oxide_reduction_capacity_g_n[unit];
        context.reactive_nitrogen.previous_non_band_microbial_ammonium_capacity_g_n[unit] = context.flux_workspace.non_band_microbial_ammonium_capacity_g_n[unit];
        context.reactive_nitrogen.previous_band_microbial_ammonium_capacity_g_n[unit] = context.flux_workspace.band_microbial_ammonium_capacity_g_n[unit];
        context.reactive_nitrogen.previous_non_band_microbial_nitrate_capacity_g_n[unit] = context.flux_workspace.non_band_microbial_nitrate_capacity_g_n[unit];
        context.reactive_nitrogen.previous_band_microbial_nitrate_capacity_g_n[unit] = context.flux_workspace.band_microbial_nitrate_capacity_g_n[unit];
        next_non_band_ammonium_demand += context.flux_workspace.non_band_ammonia_oxidation_capacity_g_n[unit];
        next_band_ammonium_demand += context.flux_workspace.band_ammonia_oxidation_capacity_g_n[unit];
        next_non_band_ammonium_demand += context.flux_workspace.non_band_microbial_ammonium_capacity_g_n[unit];
        next_band_ammonium_demand += context.flux_workspace.band_microbial_ammonium_capacity_g_n[unit];
        next_non_band_nitrite_demand += context.flux_workspace.non_band_nitrite_oxidation_capacity_g_n[unit] + context.flux_workspace.non_band_nitrite_reduction_capacity_g_n[unit];
        next_band_nitrite_demand += context.flux_workspace.band_nitrite_oxidation_capacity_g_n[unit] + context.flux_workspace.band_nitrite_reduction_capacity_g_n[unit];
        next_non_band_nitrate_demand += context.flux_workspace.non_band_nitrate_reduction_capacity_g_n[unit];
        next_band_nitrate_demand += context.flux_workspace.band_nitrate_reduction_capacity_g_n[unit];
        next_non_band_nitrate_demand += context.flux_workspace.non_band_microbial_nitrate_capacity_g_n[unit];
        next_band_nitrate_demand += context.flux_workspace.band_microbial_nitrate_capacity_g_n[unit];
        next_nitrous_oxide_demand += context.flux_workspace.nitrous_oxide_reduction_capacity_g_n[unit];
        context.phosphorus_history.previous_non_band_h2po4_capacity_g_p[unit] = context.flux_workspace.non_band_microbial_h2po4_capacity_g_p[unit];
        context.phosphorus_history.previous_band_h2po4_capacity_g_p[unit] = context.flux_workspace.band_microbial_h2po4_capacity_g_p[unit];
        context.phosphorus_history.previous_non_band_hpo4_capacity_g_p[unit] = context.flux_workspace.non_band_microbial_hpo4_capacity_g_p[unit];
        context.phosphorus_history.previous_band_hpo4_capacity_g_p[unit] = context.flux_workspace.band_microbial_hpo4_capacity_g_p[unit];
        next_h2po4_demand[0] += context.flux_workspace.non_band_microbial_h2po4_capacity_g_p[unit];
        next_h2po4_demand[1] += context.flux_workspace.band_microbial_h2po4_capacity_g_p[unit];
        next_hpo4_demand[0] += context.flux_workspace.non_band_microbial_hpo4_capacity_g_p[unit];
        next_hpo4_demand[1] += context.flux_workspace.band_microbial_hpo4_capacity_g_p[unit];
    }
    context.reactive_nitrogen.previous_total_non_band_ammonium_demand_g_n[layer] = next_non_band_ammonium_demand;
    context.reactive_nitrogen.previous_total_band_ammonium_demand_g_n[layer] = next_band_ammonium_demand;
    context.reactive_nitrogen.previous_total_non_band_nitrite_demand_g_n[layer] = next_non_band_nitrite_demand;
    context.reactive_nitrogen.previous_total_band_nitrite_demand_g_n[layer] = next_band_nitrite_demand;
    context.reactive_nitrogen.previous_total_non_band_nitrate_demand_g_n[layer] = next_non_band_nitrate_demand;
    context.reactive_nitrogen.previous_total_band_nitrate_demand_g_n[layer] = next_band_nitrate_demand;
    context.reactive_nitrogen.previous_total_nitrous_oxide_demand_g_n[layer] = next_nitrous_oxide_demand;
    context.phosphorus_history.previous_total_non_band_h2po4_demand_g_p[layer] = next_h2po4_demand[0];
    context.phosphorus_history.previous_total_band_h2po4_demand_g_p[layer] = next_h2po4_demand[1];
    context.phosphorus_history.previous_total_non_band_hpo4_demand_g_p[layer] = next_hpo4_demand[0];
    context.phosphorus_history.previous_total_band_hpo4_demand_g_p[layer] = next_hpo4_demand[1];
    context.reactive_nitrogen.current_nitrification_inhibition_activity[layer] = context.flux_workspace.layer_nitrification_inhibition_activity[layer];
    context.gas_state.dissolved_mass_g[n2o_index] = nitrous_oxide;
    context.gas_state.dissolved_mass_g[n2_index] = dinitrogen;
    context.gas_state.dissolved_mass_g[co2_index] = carbon_dioxide;
    context.gas_state.dissolved_mass_g[methane_index] = methane;
    if (context.hourly_signed_heterotrophic_respiration_g_c) |ledger|
        ledger[layer] = -((carbon_dioxide - carbon_dioxide_before) + (methane - methane_before));
    if (context.hourly_carbon_dioxide_production_g_c) |ledger|
        ledger[layer] = respiratory_carbon_dioxide_g_c + if (context.methane) |methane_result| methane_result.methane_oxidation_respiration_g_c[layer] else 0;
    context.gas_state.dissolved_mass_g[hydrogen_index] = hydrogen;
    context.gas_state.dissolved_mass_g[oxygen_index] = aqueous_oxygen;
    context.gas_state.gaseous_mass_g[oxygen_index] = gaseous_oxygen;
    if (context.methane) |methane_result| context.gas_state.gaseous_mass_g[methane_index] = methane_result.gaseous_methane_after_g_c[layer];
    for (0..organic.substrate_count) |complex| {
        const mobile = layer * organic.substrate_count + complex;
        context.organic_state.dissolved[mobile] = dissolved_after[complex];
        context.organic_state.dissolved_acetate_carbon_g_c[mobile] = dissolved_acetate_after_g_c[complex];
        context.organic_state.adsorbed[mobile] = adsorbed_after[complex];
        context.organic_state.adsorbed_acetate_carbon_g_c[mobile] = adsorbed_acetate_after_g_c[complex];
    }
    @memcpy(context.organic_state.residue[layer * organic.substrate_count * organic.residue_fraction_count ..][0 .. organic.substrate_count * organic.residue_fraction_count], &residue_after);
    @memcpy(context.organic_state.structural[structural_first..][0..structural_after.len], &structural_after);
    @memcpy(context.organic_state.colonized_structural_carbon_g_c[structural_first..][0..colonized_after.len], &colonized_after);
    context.organic_state.structural[humus_first] = humus_after[0];
    context.organic_state.structural[humus_first + 1] = humus_after[1];
    context.organic_state.colonized_structural_carbon_g_c[humus_first] = humus_colonized_after[0];
    context.organic_state.colonized_structural_carbon_g_c[humus_first + 1] = humus_colonized_after[1];
    for (0..context.microbial_state.substrate_count) |substrate| for (0..context.microbial_state.population_count) |population| {
        if (!microbial.nitroPopulationEnabled(substrate, population)) continue;
        const runtime_index = try context.microbial_state.populationIndex(layer / context.microbial_state.layer_count, layer % context.microbial_state.layer_count, substrate, population);
        const unit = first + substrate * context.microbial_state.population_count + population;
        if (context.organic_priming != null and substrate < organic.substrate_count) {
            const priming_base = (layer * organic.substrate_count * context.microbial_state.population_count + substrate * context.microbial_state.population_count + population) * organic.kinetic_fraction_count;
            const changes = context.organic_priming.?.exchange.microbial_change[priming_base..][0..organic.kinetic_fraction_count];
            context.microbial_state.structural[runtime_index * 2].carbon_g_c += changes[0].carbon_g_c;
            context.microbial_state.structural[runtime_index * 2].nitrogen_g_n += changes[0].nitrogen_g_n;
            context.microbial_state.structural[runtime_index * 2].phosphorus_g_p += changes[0].phosphorus_g_p;
            context.microbial_state.structural[runtime_index * 2 + 1].carbon_g_c += changes[1].carbon_g_c;
            context.microbial_state.structural[runtime_index * 2 + 1].nitrogen_g_n += changes[1].nitrogen_g_n;
            context.microbial_state.structural[runtime_index * 2 + 1].phosphorus_g_p += changes[1].phosphorus_g_p;
        }
        context.microbial_state.nonstructural[runtime_index] = next_nonstructural[unit - first];
        for (0..2) |component| {
            const basal = context.microbial_turnover.basal[unit * 2 + component];
            const senescence = context.microbial_turnover.senescence[unit * 2 + component];
            context.microbial_state.structural[runtime_index * 2 + component].carbon_g_c -= basal.decomposed.carbon_g_c + senescence.decomposed.carbon_g_c;
            context.microbial_state.structural[runtime_index * 2 + component].nitrogen_g_n -= basal.decomposed.nitrogen_g_n + senescence.decomposed.nitrogen_g_n;
            context.microbial_state.structural[runtime_index * 2 + component].phosphorus_g_p -= basal.decomposed.phosphorus_g_p + senescence.decomposed.phosphorus_g_p;
        }
        context.microbial_state.structural[runtime_index * 2].carbon_g_c += context.flux_workspace.labile_assimilation_g_c[unit];
        context.microbial_state.structural[runtime_index * 2].nitrogen_g_n += context.flux_workspace.labile_assimilation_g_n[unit];
        context.microbial_state.structural[runtime_index * 2].phosphorus_g_p += context.flux_workspace.labile_assimilation_g_p[unit];
        context.microbial_state.structural[runtime_index * 2 + 1].carbon_g_c += context.flux_workspace.resistant_assimilation_g_c[unit];
        context.microbial_state.structural[runtime_index * 2 + 1].nitrogen_g_n += context.flux_workspace.resistant_assimilation_g_n[unit];
        context.microbial_state.structural[runtime_index * 2 + 1].phosphorus_g_p += context.flux_workspace.resistant_assimilation_g_p[unit];
    };
    var closed_carbon_dioxide_g_c = context.gas_state.dissolved_mass_g[co2_index];
    // A second census is sometimes required because adding a sub-gram closure
    // to landscape-scale pools can itself round at a different exponent.
    for (0..4) |_| {
        const conserved_carbon_after_g_c = try authoritativeLayerCarbon_g_c(context.*, layer);
        // Same representation floor as the nitrogen closure below; see
        // `normalizeCensusRoundoff`. Carbon has not yet been observed to fail
        // this way, but the two loops are structurally identical and leaving
        // one unfiltered would only defer the same defect to the first example
        // whose CO2 carrier happens to be empty.
        const carbon_closure_g_c = normalizeCensusRoundoff(
            conserved_carbon_before_g_c - conserved_carbon_after_g_c,
            @max(@abs(conserved_carbon_before_g_c), @abs(conserved_carbon_after_g_c)),
        );
        if (carbon_closure_g_c == 0) break;
        closed_carbon_dioxide_g_c += carbon_closure_g_c;
        if (!std.math.isFinite(closed_carbon_dioxide_g_c) or closed_carbon_dioxide_g_c < 0)
            return error.InvalidSoilCarbonConservationClosure;
        context.gas_state.dissolved_mass_g[co2_index] = closed_carbon_dioxide_g_c;
    }
    if (context.hourly_signed_heterotrophic_respiration_g_c) |ledger|
        ledger[layer] = -((closed_carbon_dioxide_g_c - carbon_dioxide_before) +
            (context.gas_state.dissolved_mass_g[methane_index] - methane_before));
    // NITRO only redistributes nitrogen among mineral, organic, microbial,
    // and gaseous owners. Close accepted-publication roundoff and incomplete
    // partition cancellation in dissolved N2 before gas transport observes
    // the layer. N2 is the least reactive internal carrier in this kernel.
    var closed_dinitrogen_g_n = context.gas_state.dissolved_mass_g[n2_index];
    var diagnostic_nitrogen_closure_g_n: f64 = 0;
    for (0..4) |_| {
        const conserved_nitrogen_after_g_n = try authoritativeLayerNitrogen_g_n(context.*, layer);
        // A census difference smaller than the f64 representation floor of the
        // census itself is not a mass discrepancy: it is the sum's own rounding,
        // and its sign is an artifact of summation order rather than physics.
        // Pushing it into a carrier is unphysical in both directions, and in a
        // frozen layer whose dissolved N2 is exactly zero it is also impossible.
        // `Arctic Tundra IQ` layer 5 at hour 1 measured `closure=-3.814697e-6`
        // against `census=1.3391410615279026e10`, a relative `2.8e-16`, with
        // `carrier=0`, `fixed_g_n=0`. This is the same reasoning, and the same
        // helper, that `normalizeNonnegativeRoundoff` already applies to every
        // dissolved, adsorbed, residue, and structural pool above; the closure
        // loop was the one place that skipped it. It is a representation floor
        // derived from the operands, not a configured tolerance, and it does
        // not widen with the size of a real imbalance.
        const nitrogen_closure_g_n = normalizeCensusRoundoff(
            conserved_nitrogen_before_g_n - conserved_nitrogen_after_g_n,
            @max(@abs(conserved_nitrogen_before_g_n), @abs(conserved_nitrogen_after_g_n)),
        );
        if (nitrogen_closure_g_n == 0) break;
        diagnostic_nitrogen_closure_g_n += nitrogen_closure_g_n;
        closed_dinitrogen_g_n += nitrogen_closure_g_n;
        if (!std.math.isFinite(closed_dinitrogen_g_n) or closed_dinitrogen_g_n < 0) {
            std.log.err(
                "soil nitrogen closure exceeds the dissolved N2 carrier: layer={d} carrier_before_g_n={e} closure_g_n={e} carrier_after_g_n={e} census_before_g_n={e} census_after_g_n={e} water_m3={e} fixed_g_n={e}",
                .{
                    layer,
                    context.gas_state.dissolved_mass_g[n2_index],
                    nitrogen_closure_g_n,
                    closed_dinitrogen_g_n,
                    conserved_nitrogen_before_g_n,
                    conserved_nitrogen_after_g_n,
                    water_m3,
                    published.fixed_dinitrogen_g_n,
                },
            );
            return error.InvalidSoilNitrogenConservationClosure;
        }
        context.gas_state.dissolved_mass_g[n2_index] = closed_dinitrogen_g_n;
    }
    if (diagnostic_nitrogen_closure_g_n != 0)
        std.log.debug("soil accepted nitrogen closure: layer={d} closure_g_n={e}", .{ layer, diagnostic_nitrogen_closure_g_n });
}

/// Complete carbon owner changed by one NITRO layer commit. The organic
/// mirror is unchanged during the transaction, so including it in both
/// snapshots cancels exactly while the authoritative microbial state records
/// its actual change. Closing the residual in dissolved CO2 prevents many
/// partition products from manufacturing carbon through accumulated rounding.
fn authoritativeLayerCarbon_g_c(context: ApplyContext, layer: usize) !f64 {
    var total = try context.organic_state.totalCarbon_g_c(layer);
    const cell = layer / context.microbial_state.layer_count;
    const local_layer = layer % context.microbial_state.layer_count;
    for (0..context.microbial_state.substrate_count) |substrate|
        for (0..context.microbial_state.population_count) |population| {
            const population_index = try context.microbial_state.populationIndex(
                cell,
                local_layer,
                substrate,
                population,
            );
            total += context.microbial_state.nonstructural[population_index].carbon_g_c;
            total += context.microbial_state.structural[population_index * 2].carbon_g_c;
            total += context.microbial_state.structural[population_index * 2 + 1].carbon_g_c;
        };
    inline for (.{ gas.Species.carbon_dioxide, gas.Species.methane }) |species| {
        const index = try gas.massIndex(layer, species, context.gas_state.cell_count);
        total += context.gas_state.gaseous_mass_g[index];
        total += context.gas_state.dissolved_mass_g[index];
        total += context.gas_state.macropore_dissolved_mass_g[index];
        total += context.gas_state.band_dissolved_mass_g[index];
    }
    if (!std.math.isFinite(total)) return error.NonFiniteSoilCarbonConservationCensus;
    return total;
}

/// Complete nitrogen ownership changed by one NITRO layer commit. Mineral
/// concentrations are converted back to tracked-element grams using the same
/// water volume and zone fractions used at publication.
fn authoritativeLayerNitrogen_g_n(context: ApplyContext, layer: usize) !f64 {
    var total: f64 = context.reactive_nitrogen.non_band_nitrite_g_n[layer] +
        context.reactive_nitrogen.band_nitrite_g_n[layer];
    const water_m3 = context.water_volume_m3[layer];
    const aqueous = context.chemistry_state.aqueous[layer];
    total += context.nitrogen_molar_mass_g_per_mol * water_m3 *
        (aqueous.ammonium_non_band * context.zone_fractions.ammonium_non_band +
            aqueous.ammonium_band * context.zone_fractions.ammonium_band +
            aqueous.nitrate_non_band * context.zone_fractions.nitrate_non_band +
            aqueous.nitrate_band * context.zone_fractions.nitrate_band);
    const organic_first = layer * organic.substrate_count;
    for (context.organic_state.dissolved[organic_first..][0..organic.substrate_count]) |pool|
        total += pool.nitrogen_g_n;
    for (context.organic_state.adsorbed[organic_first..][0..organic.substrate_count]) |pool|
        total += pool.nitrogen_g_n;
    const residue_first = layer * organic.substrate_count * organic.residue_fraction_count;
    for (context.organic_state.residue[residue_first..][0 .. organic.substrate_count * organic.residue_fraction_count]) |pool|
        total += pool.nitrogen_g_n;
    const structural_first = layer * organic.substrate_count * organic.structural_fraction_count;
    for (context.organic_state.structural[structural_first..][0 .. organic.substrate_count * organic.structural_fraction_count]) |pool|
        total += pool.nitrogen_g_n;
    const cell = layer / context.microbial_state.layer_count;
    const local_layer = layer % context.microbial_state.layer_count;
    for (0..context.microbial_state.substrate_count) |substrate|
        for (0..context.microbial_state.population_count) |population| {
            const population_index = try context.microbial_state.populationIndex(cell, local_layer, substrate, population);
            total += context.microbial_state.nonstructural[population_index].nitrogen_g_n;
            total += context.microbial_state.structural[population_index * 2].nitrogen_g_n;
            total += context.microbial_state.structural[population_index * 2 + 1].nitrogen_g_n;
        };
    inline for (.{ gas.Species.nitrogen, gas.Species.nitrous_oxide, gas.Species.ammonia }) |species| {
        const index = try gas.massIndex(layer, species, context.gas_state.cell_count);
        total += context.gas_state.gaseous_mass_g[index];
        total += context.gas_state.dissolved_mass_g[index];
        total += context.gas_state.macropore_dissolved_mass_g[index];
        total += context.gas_state.band_dissolved_mass_g[index];
    }
    if (!std.math.isFinite(total)) return error.NonFiniteSoilNitrogenConservationCensus;
    return total;
}

/// NITRO RDOSL includes lignin (fraction 4 in one-based Fortran) from the
/// first three residue complexes only.
fn ligninDecompositionForAcidity(
    decomposition: *const organic_decomposition.State,
    layer: usize,
) f64 {
    const lignin_fraction: usize = 3;
    var total_g_c: f64 = 0;
    for (0..3) |substrate| {
        const index = (layer * organic.substrate_count + substrate) *
            organic.structural_fraction_count + lignin_fraction;
        total_g_c += decomposition.structural_decomposition[index].carbon_g_c;
    }
    return total_g_c;
}

fn normalizeNonnegativeRoundoff(value: f64, operation_scale: f64) f64 {
    if (!std.math.isFinite(value) or value >= 0 or
        !std.math.isFinite(operation_scale) or operation_scale < 0)
        return value;
    const cancellation_tolerance =
        64.0 * std.math.floatEps(f64) *
        @max(std.math.floatMin(f64), operation_scale);
    return if (value >= -cancellation_tolerance) 0 else value;
}

/// Discards a conservation-census difference that is smaller than the f64
/// representation floor of the census being differenced. Such a difference
/// carries no information about mass: it is the accumulation's own rounding,
/// and both its magnitude and its sign depend on summation order.
///
/// This is deliberately two-sided, unlike `normalizeNonnegativeRoundoff`. A
/// one-sided filter would let positive noise accumulate into a carrier while
/// rejecting negative noise, which is precisely the systematic drift that lane
/// A9 observed in `Arctic Tundra IQ` (every accepted closure negative).
///
/// The floor is derived from the operands via `floatEps`, so it is a property
/// of f64 and of the census magnitude, not a configured tolerance. It cannot be
/// widened to absorb a real imbalance: a discrepancy one part in `1e15` of a
/// `1e10 g N` census is `1e-5 g N`, still far below any physical process this
/// kernel represents, while anything physically meaningful is orders of
/// magnitude above the floor and passes through unchanged.
fn normalizeCensusRoundoff(difference: f64, census_scale: f64) f64 {
    if (!std.math.isFinite(difference) or !std.math.isFinite(census_scale) or census_scale < 0)
        return difference;
    const representation_floor =
        64.0 * std.math.floatEps(f64) * @max(std.math.floatMin(f64), census_scale);
    return if (@abs(difference) <= representation_floor) 0 else difference;
}

test "census roundoff filter discards representation noise symmetrically" {
    // The measured Arctic Tundra IQ layer-5 case, and its mirror image.
    const census: f64 = 1.3391410615279026e10;
    try std.testing.expectEqual(@as(f64, 0), normalizeCensusRoundoff(-3.814697265625e-6, census));
    try std.testing.expectEqual(@as(f64, 0), normalizeCensusRoundoff(3.814697265625e-6, census));
    // A physically meaningful imbalance at the same census scale survives. One
    // gram is 12 orders of magnitude above the floor here.
    try std.testing.expectEqual(@as(f64, 1), normalizeCensusRoundoff(1, census));
    try std.testing.expectEqual(@as(f64, -1), normalizeCensusRoundoff(-1, census));
    // The floor scales with the census and never exceeds it.
    try std.testing.expect(normalizeCensusRoundoff(1e-9, 1) == 1e-9);
    // Non-finite input is passed through so the caller's own guard reports it.
    try std.testing.expect(std.math.isNan(normalizeCensusRoundoff(std.math.nan(f64), census)));
}

fn applyMicrobialExchange(pool_g_n: *f64, exchange_g_n: f64, tolerance: f64) !void {
    if (!std.math.isFinite(exchange_g_n)) return error.InvalidSoilNitrogenFlux;
    if (exchange_g_n > pool_g_n.* + tolerance) return error.InsufficientSoilMineralNitrogen;
    pool_g_n.* = @max(0, pool_g_n.* - exchange_g_n);
}

fn addPool(destination: *organic.ElementPool, source: organic.ElementPool) void {
    destination.carbon_g_c += source.carbon_g_c;
    destination.nitrogen_g_n += source.nitrogen_g_n;
    destination.phosphorus_g_p += source.phosphorus_g_p;
}

fn subtractPool(destination: *organic.ElementPool, source: organic.ElementPool) void {
    destination.carbon_g_c -= source.carbon_g_c;
    destination.nitrogen_g_n -= source.nitrogen_g_n;
    destination.phosphorus_g_p -= source.phosphorus_g_p;
}

fn applyZone(ammonium: *f64, nitrate: *f64, nitrite: *f64, ammonia_oxidation: f64, nitrite_oxidation: f64, nitrate_reduction: f64, heterotrophic_nitrite_reduction: f64, autotrophic_nitrite_reduction: f64, autotrophic_ammonium_oxidation: f64, chemo_nitrite_reduction: f64, tolerance: f64) !void {
    inline for (.{ ammonia_oxidation, nitrite_oxidation, nitrate_reduction, heterotrophic_nitrite_reduction, autotrophic_nitrite_reduction, autotrophic_ammonium_oxidation, chemo_nitrite_reduction }) |value| if (!std.math.isFinite(value) or value < 0) return error.InvalidSoilNitrogenFlux;
    if (ammonia_oxidation + autotrophic_ammonium_oxidation > ammonium.* + tolerance or nitrate_reduction > nitrate.* + tolerance) return error.InsufficientSoilMineralNitrogen;
    const nitrite_available = nitrite.* + ammonia_oxidation + autotrophic_ammonium_oxidation + nitrate_reduction;
    const nitrite_consumed = nitrite_oxidation + heterotrophic_nitrite_reduction + autotrophic_nitrite_reduction + chemo_nitrite_reduction;
    if (nitrite_consumed > nitrite_available + tolerance) return error.InsufficientSoilNitrite;
    ammonium.* = @max(0, ammonium.* - ammonia_oxidation - autotrophic_ammonium_oxidation);
    nitrate.* = @max(0, nitrate.* + nitrite_oxidation - nitrate_reduction);
    nitrite.* = @max(0, nitrite_available - nitrite_consumed);
}

fn concentrationFromMass(mass_g_n: f64, water_m3: f64, fraction: f64, molar_mass: f64, tolerance: f64) !f64 {
    if (!std.math.isFinite(mass_g_n) or mass_g_n < -tolerance) return error.InvalidSoilNitrogenCommit;
    const volume = water_m3 * fraction;
    if (volume <= 0) {
        if (mass_g_n > tolerance) return error.NitrogenInZeroVolumeZone;
        return 0;
    }
    return @max(0, mass_g_n) / (volume * molar_mass);
}

fn validatePublishedAmmonium(
    layer: usize,
    water_mol_per_m3: f64,
    concentration_mol_per_m3: [2]f64,
    final_mass_g_n: [2]f64,
    microbial_exchange_g_n: [2]f64,
    oxidation_g_n: [2]f64,
    water_volume_m3: f64,
    zone_fraction: [2]f64,
) !void {
    if (!std.math.isFinite(water_mol_per_m3) or water_mol_per_m3 <= 0)
        return error.InvalidSoilNitrogenWaterMolarity;
    for (0..2) |zone| {
        if (concentration_mol_per_m3[zone] <= water_mol_per_m3) continue;
        std.log.warn(
            "soil nitrogen publication exceeds water molarity: layer={d} zone={s} ammonium_mol_per_m3={e} water_mol_per_m3={e} final_mass_g_n={e} microbial_exchange_g_n={e} oxidation_g_n={e} water_volume_m3={e} zone_fraction={e}",
            .{
                layer,
                if (zone == 0) "non_band" else "band",
                concentration_mol_per_m3[zone],
                water_mol_per_m3,
                final_mass_g_n[zone],
                microbial_exchange_g_n[zone],
                oxidation_g_n[zone],
                water_volume_m3,
                zone_fraction[zone],
            },
        );
        return error.SoilNitrogenPublicationExceedsWaterMolarity;
    }
}

fn boundedExchange(proposed: f64, dissolved_available: f64, sorbed_available: f64) f64 {
    if (proposed >= 0) return @min(proposed, @max(0, dissolved_available));
    return @max(proposed, -@max(0, sorbed_available));
}

fn validate(context: ApplyContext, range: compute.CellRange) !void {
    const layers = context.reactive_nitrogen.layer_count;
    if (context.hourly_signed_heterotrophic_respiration_g_c) |ledger| if (ledger.len != layers) return error.HeterotrophicRespirationLedgerDimensionMismatch;
    if (context.hourly_carbon_dioxide_production_g_c) |ledger| if (ledger.len != layers) return error.CarbonDioxideProductionLedgerDimensionMismatch;
    const units = try std.math.mul(usize, layers, context.reactive_nitrogen.process_unit_count_per_layer);
    if (range.first > range.end or range.end > layers or context.phosphorus_history.layer_count != layers or context.microbial_turnover.layer_count != layers or context.litter_colonization.layer_count != layers or context.organic_sorption.layer_count != layers or context.organic_decomposition.layer_count != layers or (context.organic_priming != null and context.organic_priming.?.exchange.cell_count != layers) or (context.respiration_products != null and context.respiration_products.?.layer_count != layers) or context.chemistry_state.cell_count != layers or context.gas_state.cell_count != layers or context.organic_state.layer_count != layers or context.microbial_state.cell_count * context.microbial_state.layer_count != layers or context.flux_workspace.layer_count != layers or context.flux_workspace.process_unit_count_per_layer != context.reactive_nitrogen.process_unit_count_per_layer or context.flux_workspace.process_unit_count_per_layer != context.phosphorus_history.process_unit_count_per_layer or context.flux_workspace.process_unit_count_per_layer != context.microbial_turnover.process_unit_count_per_layer or context.flux_workspace.process_unit_count_per_layer != context.microbial_state.substrate_count * context.microbial_state.population_count or context.water_volume_m3.len != layers or context.oxygen_satisfaction_fraction.len != units or context.redox_satisfaction_fraction.len != units or context.soil_layer_capacity == 0 or context.humus_partition_by_cell.len * context.soil_layer_capacity != layers) return error.SoilNitrogenCommitDimensionMismatch;
    if (!std.math.isFinite(context.nitrogen_molar_mass_g_per_mol) or context.nitrogen_molar_mass_g_per_mol <= 0 or !std.math.isFinite(context.phosphorus_molar_mass_g_per_mol) or context.phosphorus_molar_mass_g_per_mol <= 0 or !std.math.isFinite(context.tolerance_g_n) or context.tolerance_g_n < 0 or !std.math.isFinite(context.timestep_h) or context.timestep_h <= 0) {
        std.log.warn("invalid nitrogen commit units: nitrogen_g_per_mol={e} phosphorus_g_per_mol={e} tolerance_g={e}", .{ context.nitrogen_molar_mass_g_per_mol, context.phosphorus_molar_mass_g_per_mol, context.tolerance_g_n });
        return error.InvalidSoilNitrogenCommit;
    }
    inline for (.{ context.zone_fractions.ammonium_non_band, context.zone_fractions.ammonium_band, context.zone_fractions.nitrate_non_band, context.zone_fractions.nitrate_band }) |fraction| if (!std.math.isFinite(fraction) or fraction < 0 or fraction > 1) return error.InvalidSoilNitrogenZoneFraction;
    for (context.humus_partition_by_cell, 0..) |partition, cell| {
        if (!std.math.isFinite(partition[0]) or !std.math.isFinite(partition[1]) or partition[0] < 0 or partition[1] < 0 or @abs(partition[0] + partition[1] - 1) > context.tolerance_g_n) {
            std.log.warn("invalid humus partition: cell={d} active_fraction={e} passive_fraction={e} sum={e} tolerance={e}", .{ cell, partition[0], partition[1], partition[0] + partition[1], context.tolerance_g_n });
            return error.InvalidSoilNitrogenCommit;
        }
    }
    for (context.oxygen_satisfaction_fraction, context.redox_satisfaction_fraction) |oxygen, redox| if (!std.math.isFinite(oxygen) or oxygen < 0 or oxygen > 1 or !std.math.isFinite(redox) or redox < 0 or redox > 1) return error.InvalidSoilNitrogenSatisfactionFraction;
}

test "day 12 ammonium is rejected at the nitrogen publication boundary" {
    try std.testing.expectError(
        error.SoilNitrogenPublicationExceedsWaterMolarity,
        validatePublishedAmmonium(
            0,
            5.5555906256719943e4,
            .{ 1.3535910109655228e6, 0 },
            .{ 1.8950274153517318e7, 0 },
            .{ -1.8950274153517318e7, 0 },
            .{ 0, 0 },
            1,
            .{ 1, 0 },
        ),
    );
    try validatePublishedAmmonium(
        0,
        5.5555906256719943e4,
        .{ 6.232257392983944e-1, 0 },
        .{ 8.725160350177521, 0 },
        .{ 0, 0 },
        .{ 0, 0 },
        1,
        .{ 1, 0 },
    );
}

test "layer nitrogen commit conserves zones gases and DON atomically" {
    var reactive_state = try reactive.State.init(std.testing.allocator, 1, 1);
    defer reactive_state.deinit();
    var phosphorus_state = try phosphorus.State.init(std.testing.allocator, 1, 1);
    defer phosphorus_state.deinit();
    var turnover_state = try turnover.State.init(std.testing.allocator, 1, 1);
    defer turnover_state.deinit();
    var colonization_state = try colonization.State.init(std.testing.allocator, 1);
    defer colonization_state.deinit();
    var sorption_state = try organic_sorption.State.init(std.testing.allocator, 1);
    defer sorption_state.deinit();
    var decomposition_state = try organic_decomposition.State.init(std.testing.allocator, 1);
    defer decomposition_state.deinit();
    turnover_state.basal[0] = .{
        .decomposed = .{ .carbon_g_c = 0.05, .nitrogen_g_n = 0.005, .phosphorus_g_p = 0.001 },
        .recycled = .{ .carbon_g_c = 0.02, .nitrogen_g_n = 0.002, .phosphorus_g_p = 0.0004 },
        .humified = .{ .carbon_g_c = 0.01, .nitrogen_g_n = 0.001, .phosphorus_g_p = 0.0002 },
        .microbial_residue = .{ .carbon_g_c = 0.02, .nitrogen_g_n = 0.002, .phosphorus_g_p = 0.0004 },
    };
    reactive_state.non_band_nitrite_g_n[0] = 0.2;
    var chemistry_state = try chemistry.State.init(std.testing.allocator, 1);
    defer chemistry_state.deinit();
    chemistry_state.aqueous[0].ammonium_non_band = 1;
    chemistry_state.aqueous[0].nitrate_non_band = 1;
    chemistry_state.aqueous[0].hydrogen = 1;
    chemistry_state.water_mol_per_m3[0] = 100;
    chemistry_state.non_band_phosphate[0].dissolved_h2po4_mol_p_per_m3 = 1;
    var gas_state = try gas.State.init(std.testing.allocator, 1);
    defer gas_state.deinit();
    gas_state.dissolved_mass_g[@intFromEnum(gas.Species.nitrogen)] = 0.1;
    var organic_state = try organic.State.init(std.testing.allocator, 1);
    defer organic_state.deinit();
    var microbial_state = try microbial.State.init(std.testing.allocator, 1, 1, 1, 1);
    defer microbial_state.deinit();
    microbial_state.nonstructural[0] = .{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.02 };
    organic_state.residue[0].carbon_g_c = 1;
    organic_state.dissolved[0].carbon_g_c = 1;
    organic_state.dissolved_acetate_carbon_g_c[0] = 0.5;
    var workspace = try fluxes.State.init(std.testing.allocator, 1, 1);
    defer workspace.deinit();
    workspace.non_band_ammonia_oxidation_potential_g_n[0] = 0.1;
    workspace.non_band_nitrite_oxidation_potential_g_n[0] = 0.05;
    workspace.non_band_nitrate_reduction_potential_g_n[0] = 0.03;
    workspace.non_band_nitrate_reduction_capacity_g_n[0] = 0.06;
    workspace.non_band_heterotrophic_nitrite_reduction_potential_g_n[0] = 0.04;
    workspace.non_band_nitrite_reduction_capacity_g_n[0] = 0.07;
    workspace.nitrous_oxide_reduction_capacity_g_n[0] = 0.08;
    workspace.chemodenitrification_non_band_nitrite_reduction_g_n[0] = 0.02;
    workspace.chemodenitrification_nitrous_oxide_production_g_n[0] = 0.01;
    workspace.chemodenitrification_dissolved_organic_nitrogen_production_g_n[0] = 0.01;
    workspace.doc_uptake_g_c[0] = 0.24;
    workspace.acetate_uptake_g_c[0] = 0.1;
    workspace.total_carbon_uptake_g_c[0] = 0.34;
    workspace.actual_aerobic_respiration_g_c[0] = 0.1;
    workspace.nonstructural_carbon_gain_g_c[0] = 0.2;
    workspace.nitrogen_fixation_respiration_g_c[0] = 0.04;
    workspace.fixed_dinitrogen_g_n[0] = 0.02;
    workspace.non_band_microbial_ammonium_exchange_g_n[0] = 0.02;
    workspace.non_band_microbial_ammonium_capacity_g_n[0] = 0.025;
    workspace.non_band_microbial_h2po4_exchange_g_p[0] = 0.02;
    workspace.non_band_microbial_h2po4_capacity_g_p[0] = 0.025;
    workspace.labile_assimilation_g_c[0] = 0.11;
    workspace.labile_assimilation_g_n[0] = 0.011;
    workspace.labile_assimilation_g_p[0] = 0.0022;
    workspace.resistant_assimilation_g_c[0] = 0.09;
    workspace.resistant_assimilation_g_n[0] = 0.009;
    workspace.resistant_assimilation_g_p[0] = 0.0018;
    const before = 14.0 + 14.0 + 0.2 + 0.1 + 0.1;
    const phosphorus_before = 31.02;
    var context: ApplyContext = .{ .reactive_nitrogen = &reactive_state, .phosphorus_history = &phosphorus_state, .chemistry_state = &chemistry_state, .gas_state = &gas_state, .organic_state = &organic_state, .microbial_state = &microbial_state, .flux_workspace = &workspace, .microbial_turnover = &turnover_state, .litter_colonization = &colonization_state, .organic_sorption = &sorption_state, .organic_decomposition = &decomposition_state, .humus_partition_by_cell = &.{.{ 0.5, 0.5 }}, .soil_layer_capacity = 1, .water_volume_m3 = &.{1}, .zone_fractions = .{ .ammonium_non_band = 1, .ammonium_band = 0, .nitrate_non_band = 1, .nitrate_band = 0, .phosphate_non_band = 1, .phosphate_band = 0 }, .oxygen_satisfaction_fraction = &.{1}, .redox_satisfaction_fraction = &.{1}, .nitrogen_molar_mass_g_per_mol = 14, .phosphorus_molar_mass_g_per_mol = 31, .tolerance_g_n = 1e-12, .dynamic_salts = true };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    const after = chemistry_state.aqueous[0].ammonium_non_band * 14 + chemistry_state.aqueous[0].nitrate_non_band * 14 + reactive_state.non_band_nitrite_g_n[0] + gas_state.dissolved_mass_g[@intFromEnum(gas.Species.nitrous_oxide)] + gas_state.dissolved_mass_g[@intFromEnum(gas.Species.nitrogen)] + organic_state.dissolved[0].nitrogen_g_n + microbial_state.nonstructural[0].nitrogen_g_n + microbial_state.structural[0].nitrogen_g_n + microbial_state.structural[1].nitrogen_g_n + organic_state.residue[0].nitrogen_g_n + organic_state.structural[20].nitrogen_g_n + organic_state.structural[21].nitrogen_g_n;
    try std.testing.expectApproxEqAbs(before, after, 1e-12);
    const phosphorus_after = chemistry_state.non_band_phosphate[0].dissolved_h2po4_mol_p_per_m3 * 31 + microbial_state.nonstructural[0].phosphorus_g_p + microbial_state.structural[0].phosphorus_g_p + microbial_state.structural[1].phosphorus_g_p + organic_state.residue[0].phosphorus_g_p + organic_state.structural[20].phosphorus_g_p + organic_state.structural[21].phosphorus_g_p;
    try std.testing.expectApproxEqAbs(phosphorus_before, phosphorus_after, 1e-12);
    try std.testing.expectEqual(@as(f64, 0.06), reactive_state.previous_total_non_band_nitrate_demand_g_n[0]);
    try std.testing.expectEqual(@as(f64, 0.07), reactive_state.previous_total_non_band_nitrite_demand_g_n[0]);
    try std.testing.expectEqual(@as(f64, 0.08), reactive_state.previous_total_nitrous_oxide_demand_g_n[0]);
    try std.testing.expectEqual(@as(f64, 0.025), reactive_state.previous_non_band_microbial_ammonium_capacity_g_n[0]);
    try std.testing.expectEqual(@as(f64, 0.025), phosphorus_state.previous_non_band_h2po4_capacity_g_p[0]);
    try std.testing.expectEqual(@as(f64, 0.025), phosphorus_state.previous_total_non_band_h2po4_demand_g_p[0]);
    const carbon_after = organic_state.dissolved[0].carbon_g_c + organic_state.dissolved_acetate_carbon_g_c[0] + microbial_state.nonstructural[0].carbon_g_c + microbial_state.structural[0].carbon_g_c + microbial_state.structural[1].carbon_g_c + organic_state.residue[0].carbon_g_c + organic_state.structural[20].carbon_g_c + organic_state.structural[21].carbon_g_c + gas_state.dissolved_mass_g[@intFromEnum(gas.Species.carbon_dioxide)];
    try std.testing.expectApproxEqAbs(@as(f64, 3.5), carbon_after, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.06), microbial_state.structural[0].carbon_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.09), microbial_state.structural[1].carbon_g_c, 1e-15);
    const expected_hydrogen_mol = 1 + 0.1429 * (0.1 - 0.03) - 0.0714 * 0.04;
    try std.testing.expectApproxEqAbs(expected_hydrogen_mol, chemistry_state.aqueous[0].hydrogen, 1e-15);
}

test "failed layer nitrogen commit publishes no partial state" {
    var reactive_state = try reactive.State.init(std.testing.allocator, 1, 1);
    defer reactive_state.deinit();
    var phosphorus_state = try phosphorus.State.init(std.testing.allocator, 1, 1);
    defer phosphorus_state.deinit();
    var turnover_state = try turnover.State.init(std.testing.allocator, 1, 1);
    defer turnover_state.deinit();
    var colonization_state = try colonization.State.init(std.testing.allocator, 1);
    defer colonization_state.deinit();
    var sorption_state = try organic_sorption.State.init(std.testing.allocator, 1);
    defer sorption_state.deinit();
    var decomposition_state = try organic_decomposition.State.init(std.testing.allocator, 1);
    defer decomposition_state.deinit();
    reactive_state.non_band_nitrite_g_n[0] = 0.2;
    var chemistry_state = try chemistry.State.init(std.testing.allocator, 1);
    defer chemistry_state.deinit();
    chemistry_state.aqueous[0].ammonium_non_band = 1;
    chemistry_state.aqueous[0].nitrate_non_band = 1;
    chemistry_state.aqueous[0].hydrogen = 1;
    var gas_state = try gas.State.init(std.testing.allocator, 1);
    defer gas_state.deinit();
    var organic_state = try organic.State.init(std.testing.allocator, 1);
    defer organic_state.deinit();
    var microbial_state = try microbial.State.init(std.testing.allocator, 1, 1, 1, 1);
    defer microbial_state.deinit();
    var workspace = try fluxes.State.init(std.testing.allocator, 1, 1);
    defer workspace.deinit();
    workspace.non_band_ammonia_oxidation_potential_g_n[0] = 15;
    const chemistry_before = chemistry_state.aqueous[0];
    const nitrite_before = reactive_state.non_band_nitrite_g_n[0];
    const gas_before = gas_state.dissolved_mass_g[0..gas.species_count].*;
    const dissolved_before = organic_state.dissolved[0];
    var context: ApplyContext = .{ .reactive_nitrogen = &reactive_state, .phosphorus_history = &phosphorus_state, .chemistry_state = &chemistry_state, .gas_state = &gas_state, .organic_state = &organic_state, .microbial_state = &microbial_state, .flux_workspace = &workspace, .microbial_turnover = &turnover_state, .litter_colonization = &colonization_state, .organic_sorption = &sorption_state, .organic_decomposition = &decomposition_state, .humus_partition_by_cell = &.{.{ 0.5, 0.5 }}, .soil_layer_capacity = 1, .water_volume_m3 = &.{1}, .zone_fractions = .{ .ammonium_non_band = 1, .ammonium_band = 0, .nitrate_non_band = 1, .nitrate_band = 0, .phosphate_non_band = 1, .phosphate_band = 0 }, .oxygen_satisfaction_fraction = &.{1}, .redox_satisfaction_fraction = &.{1}, .nitrogen_molar_mass_g_per_mol = 14, .phosphorus_molar_mass_g_per_mol = 31, .tolerance_g_n = 1e-12, .dynamic_salts = true };
    try std.testing.expectError(error.InsufficientSoilMineralNitrogen, applyTile(&context, .{ .first = 0, .end = 1 }));
    try std.testing.expectEqualDeep(chemistry_before, chemistry_state.aqueous[0]);
    try std.testing.expectEqual(nitrite_before, reactive_state.non_band_nitrite_g_n[0]);
    try std.testing.expectEqualSlices(f64, &gas_before, gas_state.dissolved_mass_g[0..gas.species_count]);
    try std.testing.expectEqualDeep(dissolved_before, organic_state.dissolved[0]);
}
