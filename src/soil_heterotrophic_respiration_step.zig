const std = @import("std");
const compute = @import("compute.zig");
const grid = @import("grid.zig");
const organic = @import("soil_organic_initialization.zig");
const microbial = @import("soil_microbial_state.zig");
const metabolism = @import("soil_microbial_metabolism.zig");
const respiration = @import("soil_microbial_respiration_activity.zig");
const reactive = @import("soil_reactive_nitrogen_state.zig");
const fluxes = @import("soil_nitrogen_flux_workspace.zig");
const nitrogen_parameters = @import("soil_nitrogen_parameters.zig");
const retention = @import("soil_water_retention.zig");

pub const ApplyContext = struct {
    result: *fluxes.State,
    reactive_nitrogen: *const reactive.State,
    organic_state: *const organic.State,
    microbial_state: *const microbial.State,
    model_grid: *const grid.GridState,
    retention_curve: []const retention.ResolvedCurve,
    matrix_bulk_volume_m3: []const f64,
    matric_plus_osmotic_potential_mpa: []const f64,
    parameters: nitrogen_parameters.Parameters,
    timestep_h: f64,
    negligible_amount: f64,
};

pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    try validate(context.*, range);
    const population_count = context.microbial_state.population_count;
    const substrate_count = @min(context.microbial_state.substrate_count, organic.substrate_count);
    const units_per_layer = context.result.process_unit_count_per_layer;
    for (range.first..range.end) |layer| {
        const active_water_m3 = try biologicallyActiveWaterM3(context.*, layer);
        context.result.layer_biologically_active_water_m3[layer] = active_water_m3;
        var total_colonized_g_c: f64 = 0;
        for (0..substrate_count) |substrate| total_colonized_g_c += colonizedAndSorbedCarbon(context.organic_state, layer, substrate);
        for (0..substrate_count) |substrate| {
            if (substrate == context.parameters.nitrifier_indices.autotrophic_substrate_index) continue;
            var total_active_g_c: f64 = 0;
            var total_previous_doc_g_c: f64 = 0;
            var total_previous_acetate_g_c: f64 = 0;
            for (0..population_count) |population| {
                const unit = layer * units_per_layer + substrate * population_count + population;
                const runtime_index = try context.microbial_state.populationIndex(layer / context.microbial_state.layer_count, layer % context.microbial_state.layer_count, substrate, population);
                total_active_g_c += context.microbial_state.structural[runtime_index * 2].carbon_g_c / context.parameters.nitrifier_environment.labile_biomass_fraction;
                total_previous_doc_g_c += context.reactive_nitrogen.previous_doc_respiration_demand_g_c[unit];
                total_previous_acetate_g_c += context.reactive_nitrogen.previous_acetate_respiration_demand_g_c[unit];
            }
            var total_doc_share: f64 = 0;
            var total_acetate_share: f64 = 0;
            for (0..population_count) |population| {
                const unit = layer * units_per_layer + substrate * population_count + population;
                const runtime_index = try context.microbial_state.populationIndex(layer / context.microbial_state.layer_count, layer % context.microbial_state.layer_count, substrate, population);
                const active_biomass_g_c = context.microbial_state.structural[runtime_index * 2].carbon_g_c / context.parameters.nitrifier_environment.labile_biomass_fraction;
                const active_fraction = if (total_active_g_c > context.negligible_amount) active_biomass_g_c / total_active_g_c else 0;
                total_doc_share += competition(context.reactive_nitrogen.previous_doc_respiration_demand_g_c[unit], total_previous_doc_g_c, active_fraction, context.parameters.denitrification.minimum_competition_fraction, context.negligible_amount);
                total_acetate_share += competition(context.reactive_nitrogen.previous_acetate_respiration_demand_g_c[unit], total_previous_acetate_g_c, active_fraction, context.parameters.denitrification.minimum_competition_fraction, context.negligible_amount);
            }
            for (0..population_count) |population| {
                const unit = layer * units_per_layer + substrate * population_count + population;
                const runtime_index = try context.microbial_state.populationIndex(layer / context.microbial_state.layer_count, layer % context.microbial_state.layer_count, substrate, population);
                const labile = context.microbial_state.structural[runtime_index * 2];
                const active_biomass_g_c = labile.carbon_g_c / context.parameters.nitrifier_environment.labile_biomass_fraction;
                const active_fraction = if (total_active_g_c > context.negligible_amount) active_biomass_g_c / total_active_g_c else 0;
                const doc_share = normalizedCompetition(competition(context.reactive_nitrogen.previous_doc_respiration_demand_g_c[unit], total_previous_doc_g_c, active_fraction, context.parameters.denitrification.minimum_competition_fraction, context.negligible_amount), total_doc_share);
                const acetate_share = normalizedCompetition(competition(context.reactive_nitrogen.previous_acetate_respiration_demand_g_c[unit], total_previous_acetate_g_c, active_fraction, context.parameters.denitrification.minimum_competition_fraction, context.negligible_amount), total_acetate_share);
                const parameters = context.parameters.heterotrophic_respiration;
                const actual_n = if (labile.carbon_g_c > 0) labile.nitrogen_g_n / labile.carbon_g_c else parameters.target_nitrogen_per_carbon_g_n_per_g_c;
                const actual_p = if (labile.carbon_g_c > 0) labile.phosphorus_g_p / labile.carbon_g_c else parameters.target_phosphorus_per_carbon_g_p_per_g_c;
                const nutrient = @min(@min(1, @max(0.1, std.math.pow(f64, actual_n / parameters.target_nitrogen_per_carbon_g_n_per_g_c, 0.25))), @min(1, @max(0.1, std.math.pow(f64, actual_p / parameters.target_phosphorus_per_carbon_g_p_per_g_c, 0.25))));
                const water = @exp(parameters.water_potential_sensitivity_per_mpa * context.matric_plus_osmotic_potential_mpa[layer]);
                const unlimited = parameters.substrate_unlimited_respiration_per_h * nutrient * water * active_biomass_g_c * context.timestep_h;
                const substrate_carbon_g_c = colonizedAndSorbedCarbon(context.organic_state, layer, substrate);
                const substrate_fraction = if (total_colonized_g_c > context.negligible_amount) substrate_carbon_g_c / total_colonized_g_c else 1;
                const mobile = layer * organic.substrate_count + substrate;
                const population_metabolism = respiration.sourceMetabolism(population);
                const temperature_response = try metabolism.growthTemperatureResponse(context.model_grid.soil_temperature_k[layer], context.parameters.microbial_thermal_adaptation_offset_k);
                const effective_water_m3 = active_water_m3 * substrate_fraction;
                const concentration_water_m3 = if (effective_water_m3 > 0) effective_water_m3 else active_water_m3;
                if (population_metabolism == .fermenting_heterotroph) {
                    const doc_concentration_g_c_per_m3 = if (concentration_water_m3 > 0) context.organic_state.dissolved[mobile].carbon_g_c / concentration_water_m3 else 0;
                    const doc_demand_g_c = unlimited * doc_concentration_g_c_per_m3 / (doc_concentration_g_c_per_m3 + parameters.doc_half_saturation_g_c_per_m3) * temperature_response;
                    const doc_supply_g_c = context.organic_state.dissolved[mobile].carbon_g_c * doc_share * parameters.doc_respiration_requirement_g_c_per_g_c * context.timestep_h;
                    const doc_respiration_g_c = @min(doc_supply_g_c, doc_demand_g_c);
                    context.result.substrate_unlimited_respiration_g_c[unit] = unlimited;
                    context.result.doc_respiration_demand_g_c[unit] = doc_respiration_g_c;
                    context.result.acetate_respiration_demand_g_c[unit] = 0;
                    context.result.substrate_limited_respiration_g_c[unit] = doc_respiration_g_c;
                    context.result.aerobic_oxygen_demand_g_o[unit] = 0;
                    context.result.aerobic_active_biomass_g_c[unit] = active_biomass_g_c;
                    context.result.aerobic_fallback_active_fraction[unit] = active_fraction;
                    context.result.doc_competition_fraction[unit] = doc_share;
                    context.result.substrate_complex_fraction[unit] = substrate_fraction;
                    continue;
                }
                if (population_metabolism == .acetotrophic_methanogen) {
                    const acetate_concentration_g_c_per_m3 = if (concentration_water_m3 > 0) context.organic_state.dissolved_acetate_carbon_g_c[mobile] / concentration_water_m3 else 0;
                    const acetate_demand_g_c = unlimited * acetate_concentration_g_c_per_m3 / (acetate_concentration_g_c_per_m3 + parameters.acetate_half_saturation_g_c_per_m3) * temperature_response;
                    const acetate_supply_g_c = context.organic_state.dissolved_acetate_carbon_g_c[mobile] * acetate_share * parameters.acetate_respiration_requirement_g_c_per_g_c * context.timestep_h;
                    const acetate_respiration_g_c = @min(acetate_supply_g_c, acetate_demand_g_c);
                    context.result.substrate_unlimited_respiration_g_c[unit] = 0;
                    context.result.doc_respiration_demand_g_c[unit] = 0;
                    context.result.acetate_respiration_demand_g_c[unit] = acetate_respiration_g_c;
                    context.result.substrate_limited_respiration_g_c[unit] = acetate_respiration_g_c;
                    context.result.aerobic_oxygen_demand_g_o[unit] = 0;
                    context.result.aerobic_active_biomass_g_c[unit] = active_biomass_g_c;
                    context.result.aerobic_fallback_active_fraction[unit] = active_fraction;
                    context.result.doc_competition_fraction[unit] = doc_share;
                    context.result.substrate_complex_fraction[unit] = substrate_fraction;
                    continue;
                }
                const aerobic_inputs: respiration.AerobicSubstrateInputs = .{
                    .unlimited_respiration_g_c = unlimited,
                    .dissolved_organic_carbon_g_c = context.organic_state.dissolved[mobile].carbon_g_c,
                    .dissolved_acetate_carbon_g_c = context.organic_state.dissolved_acetate_carbon_g_c[mobile],
                    .biologically_active_water_m3 = active_water_m3,
                    .substrate_complex_fraction = substrate_fraction,
                    .doc_half_saturation_g_c_per_m3 = parameters.doc_half_saturation_g_c_per_m3,
                    .acetate_half_saturation_g_c_per_m3 = parameters.acetate_half_saturation_g_c_per_m3,
                    .doc_biological_demand_fraction = doc_share,
                    .acetate_biological_demand_fraction = acetate_share,
                    .doc_respiration_requirement_g_c_per_g_c = parameters.doc_respiration_requirement_g_c_per_g_c,
                    .acetate_respiration_requirement_g_c_per_g_c = parameters.acetate_respiration_requirement_g_c_per_g_c,
                    .temperature_response = temperature_response,
                    .timestep_h = context.timestep_h,
                };
                const limited = respiration.aerobicSubstrateLimitedRespiration(aerobic_inputs) catch |err| {
                    std.log.err(
                        "invalid soil aerobic substrate inputs: layer={d} substrate={d} population={d} unlimited_g_c={e} doc_g_c={e} acetate_g_c={e} active_water_m3={e} substrate_fraction={e} doc_share={e} acetate_share={e} temperature_response={e}",
                        .{ layer, substrate, population, aerobic_inputs.unlimited_respiration_g_c, aerobic_inputs.dissolved_organic_carbon_g_c, aerobic_inputs.dissolved_acetate_carbon_g_c, aerobic_inputs.biologically_active_water_m3, aerobic_inputs.substrate_complex_fraction, aerobic_inputs.doc_biological_demand_fraction, aerobic_inputs.acetate_biological_demand_fraction, aerobic_inputs.temperature_response },
                    );
                    return err;
                };
                context.result.substrate_unlimited_respiration_g_c[unit] = unlimited;
                context.result.doc_respiration_demand_g_c[unit] = limited.doc_respiration_g_c;
                context.result.acetate_respiration_demand_g_c[unit] = limited.acetate_respiration_g_c;
                context.result.substrate_limited_respiration_g_c[unit] = limited.substrate_limited_respiration_g_c;
                context.result.aerobic_oxygen_demand_g_o[unit] = parameters.oxygen_per_respired_carbon_g_o_per_g_c * limited.substrate_limited_respiration_g_c;
                context.result.aerobic_active_biomass_g_c[unit] = active_biomass_g_c;
                context.result.aerobic_fallback_active_fraction[unit] = active_fraction;
                context.result.doc_competition_fraction[unit] = doc_share;
                context.result.substrate_complex_fraction[unit] = substrate_fraction;
            }
        }
    }
}

fn colonizedAndSorbedCarbon(state: *const organic.State, layer: usize, substrate: usize) f64 {
    var total = state.adsorbed[layer * organic.substrate_count + substrate].carbon_g_c + state.adsorbed_acetate_carbon_g_c[layer * organic.substrate_count + substrate];
    for (0..organic.structural_fraction_count) |fraction| total += state.colonized_structural_carbon_g_c[(layer * organic.substrate_count + substrate) * organic.structural_fraction_count + fraction];
    for (0..organic.residue_fraction_count) |fraction| total += state.residue[(layer * organic.substrate_count + substrate) * organic.residue_fraction_count + fraction].carbon_g_c;
    return total;
}

pub fn substrateComplexFractions(state: *const organic.State, layer: usize, output: []f64, negligible_carbon_g_c: f64) !void {
    if (layer >= state.layer_count or output.len != organic.substrate_count or !std.math.isFinite(negligible_carbon_g_c) or negligible_carbon_g_c < 0) return error.InvalidSubstrateComplexFractionRequest;
    var total: f64 = 0;
    for (0..organic.substrate_count) |substrate| {
        output[substrate] = colonizedAndSorbedCarbon(state, layer, substrate);
        if (!std.math.isFinite(output[substrate]) or output[substrate] < 0) return error.InvalidSubstrateComplexCarbon;
        total += output[substrate];
    }
    if (!std.math.isFinite(total)) return error.InvalidSubstrateComplexCarbon;
    for (output) |*fraction| fraction.* = if (total > negligible_carbon_g_c) fraction.* / total else 1;
}

fn biologicallyActiveWaterM3(context: ApplyContext, layer: usize) !f64 {
    const inactive = try context.retention_curve[layer].waterFractionAtPotentialMpa(context.parameters.oxygen_uptake.hygroscopic_water_potential_mpa);
    const water_fraction = context.model_grid.matrix_liquid_water_m3[layer] / context.matrix_bulk_volume_m3[layer];
    const active_fraction = @max(0, @min(@max(0.75 * context.retention_curve[layer].porosity_fraction, context.retention_curve[layer].curve.field_capacity_fraction), water_fraction) - inactive);
    return active_fraction / (1 + active_fraction) * context.matrix_bulk_volume_m3[layer];
}

fn competition(previous: f64, total: f64, fallback: f64, minimum: f64, negligible: f64) f64 {
    return @max(minimum, if (total > negligible) previous / total else fallback);
}

fn normalizedCompetition(raw_share: f64, total_raw_share: f64) f64 {
    return raw_share / @max(1, total_raw_share);
}

test "minimum competition floors remain a conservative shared allocation" {
    const raw = [_]f64{ 1, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001 };
    var total_raw: f64 = 0;
    for (raw) |share| total_raw += share;
    var normalized_total: f64 = 0;
    for (raw) |share| normalized_total += normalizedCompetition(share, total_raw);
    try std.testing.expectApproxEqAbs(@as(f64, 1), normalized_total, 1e-15);
}

fn validate(context: ApplyContext, range: compute.CellRange) !void {
    const layers = context.model_grid.layer_count;
    const units = try std.math.mul(usize, layers, context.result.process_unit_count_per_layer);
    if (range.first > range.end or range.end > layers or context.organic_state.layer_count != layers or context.reactive_nitrogen.layer_count != layers or context.retention_curve.len != layers or context.matrix_bulk_volume_m3.len != layers or context.matric_plus_osmotic_potential_mpa.len != layers or context.result.process_unit_count_per_layer != context.microbial_state.substrate_count * context.microbial_state.population_count or context.reactive_nitrogen.previous_doc_respiration_demand_g_c.len != units) return error.SoilHeterotrophicRespirationDimensionMismatch;
    if (context.parameters.nitrifier_indices.heterotrophic_denitrifier_population_index >= context.microbial_state.population_count or !std.math.isFinite(context.timestep_h) or context.timestep_h <= 0 or !std.math.isFinite(context.negligible_amount) or context.negligible_amount < 0) return error.InvalidSoilHeterotrophicRespirationInput;
}

test "runtime aerobic respiration evaluates every heterotrophic population" {
    const config = @import("config.zig").SimulationConfig{ .grid_columns = 1, .grid_rows = 1, .soil_layers = 1, .plant_populations = 1, .worker_threads = 1, .tile_cells = 1, .relative_tolerance = 1e-8, .absolute_tolerance = 1e-12, .max_nonlinear_iterations = 20, .picard_relaxation = 0.5 };
    var model_grid = try grid.GridState.init(std.testing.allocator, config);
    defer model_grid.deinit();
    model_grid.matrix_liquid_water_m3[0] = 0.3;
    model_grid.soil_temperature_k[0] = 293.15;
    var microbial_state = try microbial.State.init(std.testing.allocator, 1, 1, 2, 6);
    defer microbial_state.deinit();
    for (0..microbial_state.substrate_count * microbial_state.population_count) |population_index| {
        microbial_state.structural[population_index * 2] = .{ .carbon_g_c = 1, .nitrogen_g_n = 0.1, .phosphorus_g_p = 0.01 };
    }
    var organic_state = try organic.State.init(std.testing.allocator, 1);
    defer organic_state.deinit();
    organic_state.dissolved[0].carbon_g_c = 10;
    organic_state.dissolved_acetate_carbon_g_c[0] = 2;
    organic_state.adsorbed[0].carbon_g_c = 5;
    var reactive_state = try reactive.State.init(std.testing.allocator, 1, 12);
    defer reactive_state.deinit();
    var result = try fluxes.State.init(std.testing.allocator, 1, 12);
    defer result.deinit();
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
    const curve = retention.ResolvedCurve{ .porosity_fraction = 0.5, .curve = .{ .field_capacity_fraction = 0.25, .wilting_point_fraction = 0.1, .saturation_water_potential_mpa = -0.0005, .field_capacity_water_potential_mpa = -0.033, .wilting_point_water_potential_mpa = -1.5, .minimum_water_potential_mpa = -1.5e12, .saturation_to_field_shape = 0.5, .below_wilting_shape = 0.5 } };
    var context: ApplyContext = .{ .result = &result, .reactive_nitrogen = &reactive_state, .organic_state = &organic_state, .microbial_state = &microbial_state, .model_grid = &model_grid, .retention_curve = &.{curve}, .matrix_bulk_volume_m3 = &.{1}, .matric_plus_osmotic_potential_mpa = &.{0}, .parameters = parsed, .timestep_h = 1, .negligible_amount = 1e-12 };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expect(result.substrate_limited_respiration_g_c[0] > 0);
    try std.testing.expect(result.substrate_limited_respiration_g_c[1] > 0);
    try std.testing.expectApproxEqAbs(result.substrate_limited_respiration_g_c[0], result.substrate_limited_respiration_g_c[1], 1e-14);
    for (0..microbial_state.substrate_count) |substrate| {
        const first_unit = substrate * microbial_state.population_count;
        try std.testing.expectEqual(@as(f64, 0), result.acetate_respiration_demand_g_c[first_unit + 3]);
        try std.testing.expectEqual(@as(f64, 0), result.doc_respiration_demand_g_c[first_unit + 4]);
    }
}
