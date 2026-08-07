const std = @import("std");
const compute = @import("../../core/compute.zig");
const organic = @import("../organic/initialization.zig");
const microbial = @import("state.zig");
const metabolism = @import("metabolism.zig");
const oxygen = @import("../gas/oxygen_allocation.zig");
const fluxes = @import("../nutrients/nitrogen_flux_workspace.zig");
const nitrogen_parameters = @import("../nutrients/nitrogen_parameters.zig");
const respiration = @import("respiration_activity.zig");

pub const ApplyContext = struct {
    result: *fluxes.State,
    organic_state: *const organic.State,
    microbial_state: *const microbial.State,
    oxygen_state: *const oxygen.State,
    parameters: nitrogen_parameters.Parameters,
    negligible_amount: f64,
};

/// Ports NITRO CGOMX/CGOMD/CGOMC and CGOQC/CGOAC/CGOMN/CGOMP.
pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    try validate(context.*, range);
    const populations = context.microbial_state.population_count;
    const substrates = @min(context.microbial_state.substrate_count, organic.substrate_count);
    const p = context.parameters.heterotrophic_respiration;
    for (range.first..range.end) |layer| for (0..substrates) |substrate| {
        var total_active_g_c: f64 = 0;
        for (0..populations) |population| {
            const index = try context.microbial_state.populationIndex(layer / context.microbial_state.layer_count, layer % context.microbial_state.layer_count, substrate, population);
            total_active_g_c += context.microbial_state.structural[index * 2].carbon_g_c / context.parameters.nitrifier_environment.labile_biomass_fraction;
        }
        const mobile = layer * organic.substrate_count + substrate;
        const dissolved = context.organic_state.dissolved[mobile];
        const dissolved_n_per_doc = if (dissolved.carbon_g_c > context.negligible_amount) dissolved.nitrogen_g_n / dissolved.carbon_g_c else 0;
        const dissolved_p_per_doc = if (dissolved.carbon_g_c > context.negligible_amount) dissolved.phosphorus_g_p / dissolved.carbon_g_c else 0;
        var aggregate_doc_uptake_g_c: f64 = 0;
        var aggregate_acetate_uptake_g_c: f64 = 0;
        for (0..populations) |population| {
            const runtime_index = try context.microbial_state.populationIndex(layer / context.microbial_state.layer_count, layer % context.microbial_state.layer_count, substrate, population);
            const unit = layer * context.result.process_unit_count_per_layer + substrate * populations + population;
            const active_g_c = context.microbial_state.structural[runtime_index * 2].carbon_g_c / context.parameters.nitrifier_environment.labile_biomass_fraction;
            const labile = context.microbial_state.structural[runtime_index * 2];
            const actual_n = if (labile.carbon_g_c > 0) labile.nitrogen_g_n / labile.carbon_g_c else p.target_nitrogen_per_carbon_g_n_per_g_c;
            const actual_p = if (labile.carbon_g_c > 0) labile.phosphorus_g_p / labile.carbon_g_c else p.target_phosphorus_per_carbon_g_p_per_g_c;
            const nitrogen_limitation = @min(1, @max(0.1, std.math.pow(f64, actual_n / p.target_nitrogen_per_carbon_g_n_per_g_c, 0.25)));
            const phosphorus_limitation = @min(1, @max(0.1, std.math.pow(f64, actual_p / p.target_phosphorus_per_carbon_g_p_per_g_c, 0.25)));
            const doc_respiration = context.result.doc_respiration_demand_g_c[unit];
            const acetate_respiration = context.result.acetate_respiration_demand_g_c[unit];
            const respiration_sum = doc_respiration + acetate_respiration;
            const doc_fraction: f64 = if (respiration_sum > context.negligible_amount) doc_respiration / respiration_sum else 1;
            const acetate_fraction = 1 - doc_fraction;
            const aerobic_requirement = p.doc_respiration_requirement_g_c_per_g_c * doc_fraction + p.acetate_respiration_requirement_g_c_per_g_c * acetate_fraction;
            const population_metabolism = respiration.sourceMetabolism(population);
            const actual_aerobic = context.result.substrate_limited_respiration_g_c[unit] * (if (population_metabolism == .aerobic_heterotroph) context.oxygen_state.demand_satisfaction_fraction[unit] else 1);
            const value = try metabolism.respirationDrivenSubstrateUptake(.{
                .is_heterotroph = true,
                .maintenance_respiration_g_c = context.result.total_maintenance_respiration_g_c[unit],
                .aerobic_respiration_g_c = actual_aerobic,
                .growth_respiration_g_c = context.result.growth_respiration_g_c[unit],
                .fixation_respiration_g_c = context.result.nitrogen_fixation_respiration_g_c[unit],
                .aerobic_growth_respiration_fraction_g_c_per_g_c = aerobic_requirement,
                .denitrification_respiration_g_c = context.result.denitrification_respiration_g_c[unit],
                .denitrification_growth_respiration_fraction_g_c_per_g_c = p.denitrification_growth_respiration_fraction_g_c_per_g_c,
                .doc_fraction_of_aerobic_carbon = doc_fraction,
                .acetate_fraction_of_aerobic_carbon = acetate_fraction,
                .dissolved_organic_nitrogen_g_n = dissolved.nitrogen_g_n,
                .dissolved_organic_phosphorus_g_p = dissolved.phosphorus_g_p,
                .population_biomass_fraction = if (total_active_g_c > context.negligible_amount) active_g_c / total_active_g_c else 0,
                .dissolved_nitrogen_per_doc_g_n_per_g_c = dissolved_n_per_doc,
                .dissolved_phosphorus_per_doc_g_p_per_g_c = dissolved_p_per_doc,
                .nitrogen_limitation_fraction = nitrogen_limitation,
                .phosphorus_limitation_fraction = phosphorus_limitation,
            });
            context.result.actual_aerobic_respiration_g_c[unit] = actual_aerobic;
            context.result.total_carbon_uptake_g_c[unit] = value.total_carbon_g_c;
            context.result.doc_uptake_g_c[unit] = value.doc_g_c;
            context.result.acetate_uptake_g_c[unit] = value.acetate_g_c;
            context.result.dissolved_organic_nitrogen_uptake_g_n[unit] = value.dissolved_organic_nitrogen_g_n;
            context.result.dissolved_organic_phosphorus_uptake_g_p[unit] = value.dissolved_organic_phosphorus_g_p;
            context.result.nonstructural_carbon_gain_g_c[unit] = value.total_carbon_g_c - actual_aerobic - context.result.denitrification_respiration_g_c[unit] - context.result.nitrogen_fixation_respiration_g_c[unit];
            aggregate_doc_uptake_g_c += value.doc_g_c;
            aggregate_acetate_uptake_g_c += value.acetate_g_c;
        }

        // NITRO limits each population's respiration against its competition
        // share before converting respiration back to CGOQC/CGOAC. Mixed
        // growth efficiencies and maintenance can make the converted sum
        // exceed a dissolved pool by roundoff-sized amounts. Apply the
        // corresponding simultaneous pool limit here so allocation remains
        // independent of population order and exactly conservative.
        const doc_scale = if (aggregate_doc_uptake_g_c > dissolved.carbon_g_c and aggregate_doc_uptake_g_c > 0) dissolved.carbon_g_c / aggregate_doc_uptake_g_c else 1;
        const available_acetate_g_c = context.organic_state.dissolved_acetate_carbon_g_c[mobile];
        const acetate_scale = if (aggregate_acetate_uptake_g_c > available_acetate_g_c and aggregate_acetate_uptake_g_c > 0) available_acetate_g_c / aggregate_acetate_uptake_g_c else 1;
        if (doc_scale < 1 or acetate_scale < 1) {
            for (0..populations) |population| {
                const unit = layer * context.result.process_unit_count_per_layer + substrate * populations + population;
                const previous_total_g_c = context.result.total_carbon_uptake_g_c[unit];
                context.result.doc_uptake_g_c[unit] *= doc_scale;
                context.result.acetate_uptake_g_c[unit] *= acetate_scale;
                const limited_total_g_c = context.result.doc_uptake_g_c[unit] + context.result.acetate_uptake_g_c[unit];
                const elemental_scale = if (previous_total_g_c > 0) limited_total_g_c / previous_total_g_c else 0;
                context.result.total_carbon_uptake_g_c[unit] = limited_total_g_c;
                context.result.dissolved_organic_nitrogen_uptake_g_n[unit] *= elemental_scale;
                context.result.dissolved_organic_phosphorus_uptake_g_p[unit] *= elemental_scale;
                const respired_g_c = context.result.actual_aerobic_respiration_g_c[unit] + context.result.denitrification_respiration_g_c[unit] + context.result.nitrogen_fixation_respiration_g_c[unit];
                // This is a signed pool change, not newly formed biomass.
                // Maintenance and, especially, nonsymbiotic N fixation may
                // draw on pre-existing nonstructural carbon when simultaneous
                // DOC/acetate competition limits current uptake. The fixation
                // kernel already bounds that draw by the available pool; the
                // atomic nitrogen commit performs the authoritative
                // nonnegative-pool check after all gains, recycling, and
                // assimilation are combined.
                context.result.nonstructural_carbon_gain_g_c[unit] =
                    limited_total_g_c - respired_g_c;
            }
        }
    };
}

fn validate(context: ApplyContext, range: compute.CellRange) !void {
    const layers = context.organic_state.layer_count;
    const units = try std.math.mul(usize, layers, context.result.process_unit_count_per_layer);
    if (range.first > range.end or range.end > layers or context.microbial_state.cell_count * context.microbial_state.layer_count != layers or context.result.process_unit_count_per_layer != context.microbial_state.substrate_count * context.microbial_state.population_count or context.oxygen_state.demand_satisfaction_fraction.len != units) return error.SoilMicrobialSubstrateUptakeDimensionMismatch;
    if (!std.math.isFinite(context.negligible_amount) or context.negligible_amount < 0) return error.InvalidSoilMicrobialSubstrateUptakeInput;
}

test "runtime soil uptake ports CGOQC CGOAC and microbial carbon gain" {
    var organic_state = try organic.State.init(std.testing.allocator, 1);
    defer organic_state.deinit();
    organic_state.dissolved[0] = .{ .carbon_g_c = 10, .nitrogen_g_n = 1, .phosphorus_g_p = 0.2 };
    var microbial_state = try microbial.State.init(std.testing.allocator, 1, 1, 6, 2);
    defer microbial_state.deinit();
    microbial_state.structural[0] = .{ .carbon_g_c = 0.55, .nitrogen_g_n = 0.055, .phosphorus_g_p = 0.0055 };
    var oxygen_state = try oxygen.State.init(std.testing.allocator, 1, 12, 1);
    defer oxygen_state.deinit();
    oxygen_state.demand_satisfaction_fraction[0] = 1;
    var result = try fluxes.State.init(std.testing.allocator, 1, 12);
    defer result.deinit();
    result.doc_respiration_demand_g_c[0] = 0.15;
    result.acetate_respiration_demand_g_c[0] = 0.05;
    result.substrate_limited_respiration_g_c[0] = 0.2;
    result.total_maintenance_respiration_g_c[0] = 0.02;
    result.growth_respiration_g_c[0] = 0.18;
    result.denitrification_respiration_g_c[0] = 0.1;
    const source =
        "soil_nitrification 0.001 0.0002 7000 14 1.4 1.4 0.125 0.125 0.3 0.1 0.5 2.667 3.429 1.143\n" ++
        "soil_denitrification 0.001 1.4 1.4 0.014 1 0.429 0.429 0.214 0.875\n" ++
        "soil_autotrophic_denitrification 0.5 0.333\nsoil_chemodenitrification 0.0005 0.001 1e-12 0.5 0 0.5\n" ++
        "nitrous_acid_dissociation_mol_per_m3 0.45\nsoil_microbial_thermal_adaptation_offset_k 0\n" ++
        "soil_nitrifier_indices 5 0 1 1\nsoil_nitrifier_environment 0.55 0.1 0.1 0.01 0.01 12 0.1\n" ++
        "soil_oxygen_uptake 1e-6 2.3866348449e11 0.064 -1.5e4 0.5 12 12 0.5 0.7 0.001 1e-12\n" ++
        "soil_heterotrophic_respiration 0.125 0.1 0.01 12 12 0.5 0.42016806722689076 0.1 2.667 0.01 0.01 1e-6 1 0.7142857142857143\n" ++
        "soil_microbial_mineral_exchange 0.014 0.0125 0.40 0.014 0.03 0.35 0.003 0.009 0.18 31\n" ++
        "soil_nonsymbiotic_nitrogen_fixation 5 6 0.25 0.02 0.14 0.25\n" ++
        "soil_microbial_turnover 0.01 0.001 0.167 0.333 0.333 0.333 0.150 0.300 0.333 0.25 2.0 5.0 1.0 0.5";
    const parsed = try nitrogen_parameters.parse(source);
    var context: ApplyContext = .{ .result = &result, .organic_state = &organic_state, .microbial_state = &microbial_state, .oxygen_state = &oxygen_state, .parameters = parsed, .negligible_amount = 1e-12 };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expect(result.total_carbon_uptake_g_c[0] > result.actual_aerobic_respiration_g_c[0] + result.denitrification_respiration_g_c[0]);
    try std.testing.expect(result.doc_uptake_g_c[0] > result.acetate_uptake_g_c[0]);
    try std.testing.expectApproxEqAbs(result.total_carbon_uptake_g_c[0] - result.actual_aerobic_respiration_g_c[0] - result.denitrification_respiration_g_c[0], result.nonstructural_carbon_gain_g_c[0], 1e-14);

    const unconstrained_doc_g_c = result.doc_uptake_g_c[0];
    organic_state.dissolved[0].carbon_g_c = unconstrained_doc_g_c * 0.999;
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expectApproxEqAbs(organic_state.dissolved[0].carbon_g_c, result.doc_uptake_g_c[0], 1e-14);
    try std.testing.expect(result.nonstructural_carbon_gain_g_c[0] >= 0);

    // Fixation is paid from the existing nonstructural pool. When current
    // dissolved supply is scarce, this hourly pool change must remain signed
    // so the atomic commit can apply and validate the drawdown.
    organic_state.dissolved[0].carbon_g_c = 1e-6;
    result.nitrogen_fixation_respiration_g_c[0] = 0.18;
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expect(result.nonstructural_carbon_gain_g_c[0] < 0);
}
