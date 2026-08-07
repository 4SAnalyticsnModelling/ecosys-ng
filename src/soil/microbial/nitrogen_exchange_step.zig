const std = @import("std");
const compute = @import("../../core/compute.zig");
const grid = @import("../../state/grid.zig");
const microbial = @import("state.zig");
const metabolism = @import("metabolism.zig");
const chemistry = @import("../solute/chemistry_state.zig");
const zones = @import("../solute/charge_classification.zig");
const reactive = @import("../nutrients/reactive_nitrogen_state.zig");
const fluxes = @import("../nutrients/nitrogen_flux_workspace.zig");
const nitrogen_parameters = @import("../nutrients/nitrogen_parameters.zig");

pub const ApplyContext = struct {
    result: *fluxes.State,
    reactive_nitrogen: *const reactive.State,
    microbial_state: *const microbial.State,
    chemistry_state: *const chemistry.State,
    model_grid: *const grid.GridState,
    matric_plus_osmotic_potential_megapascal: []const f64,
    zone_fractions: zones.ZoneFractions,
    parameters: nitrogen_parameters.Parameters,
    nitrogen_molar_mass_g_per_mol: f64,
    timestep_h: f64,
    negligible_amount: f64,
};

/// Ports NITRO RINHP/RINHO/RINHB/RINH4/RINB4 and the following
/// RINOP/RINOO/RINOB/RINO3/RINB3 sequence. Positive flux immobilizes
/// aqueous N; negative flux mineralizes microbial storage.
pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    try validate(context.*, range);
    const populations = context.microbial_state.population_count;
    const substrates = context.microbial_state.substrate_count;
    const p = context.parameters.microbial_mineral_exchange;
    const surface_area_m2_per_g_c = context.parameters.oxygen_uptake.microbial_count_per_g_c * 4 * std.math.pi * std.math.pow(f64, context.parameters.oxygen_uptake.microbial_radius_m, 2);
    for (range.first..range.end) |layer| {
        const water_m3 = context.model_grid.matrix_liquid_water_m3[layer];
        const aqueous = context.chemistry_state.aqueous[layer];
        const ammonium_concentration = [2]f64{ aqueous.ammonium_non_band * context.nitrogen_molar_mass_g_per_mol, aqueous.ammonium_band * context.nitrogen_molar_mass_g_per_mol };
        const nitrate_concentration = [2]f64{ aqueous.nitrate_non_band * context.nitrogen_molar_mass_g_per_mol, aqueous.nitrate_band * context.nitrogen_molar_mass_g_per_mol };
        const ammonium_fraction = [2]f64{ context.zone_fractions.ammonium_non_band, context.zone_fractions.ammonium_band };
        const nitrate_fraction = [2]f64{ context.zone_fractions.nitrate_non_band, context.zone_fractions.nitrate_band };
        var total_active_g_c: f64 = 0;
        const first = layer * context.result.process_unit_count_per_layer;
        for (0..substrates) |substrate| for (0..populations) |population| {
            const runtime_index = try context.microbial_state.populationIndex(layer / context.microbial_state.layer_count, layer % context.microbial_state.layer_count, substrate, population);
            total_active_g_c += context.microbial_state.structural[runtime_index * 2].carbon_g_c / context.parameters.nitrifier_environment.labile_biomass_fraction;
        };
        for (0..substrates) |substrate| for (0..populations) |population| {
            const runtime_index = try context.microbial_state.populationIndex(layer / context.microbial_state.layer_count, layer % context.microbial_state.layer_count, substrate, population);
            const unit = first + substrate * populations + population;
            const nonstructural = context.microbial_state.nonstructural[runtime_index];
            const active_g_c = context.microbial_state.structural[runtime_index * 2].carbon_g_c / context.parameters.nitrifier_environment.labile_biomass_fraction;
            const active_fraction = if (total_active_g_c > context.negligible_amount) active_g_c / total_active_g_c else 0;
            const activity = try metabolism.growthTemperatureResponse(context.model_grid.soil_temperature_k[layer], context.parameters.microbial_thermal_adaptation_offset_k) * @exp(context.parameters.heterotrophic_respiration.water_potential_sensitivity_per_megapascal * context.matric_plus_osmotic_potential_megapascal[layer]);
            const area_activity = surface_area_m2_per_g_c * active_g_c * activity * context.timestep_h;
            const nitrogen_demand = nonstructural.carbon_g_c * context.parameters.heterotrophic_respiration.target_nitrogen_per_carbon_g_n_per_g_c - nonstructural.nitrogen_g_n;
            var ammonium_exchange = [2]f64{ 0, 0 };
            var ammonium_capacity = [2]f64{ 0, 0 };
            if (nitrogen_demand > 0) {
                const limited_demand = @min(nitrogen_demand, p.ammonium_maximum_uptake_g_n_per_m2_h * area_activity);
                inline for (0..2) |zone| {
                    const available_concentration = @max(0, ammonium_concentration[zone] - p.ammonium_minimum_concentration_g_n_per_m3);
                    ammonium_capacity[zone] = ammonium_fraction[zone] * limited_demand * available_concentration / (available_concentration + p.ammonium_half_saturation_g_n_per_m3);
                    const total_previous = if (zone == 0) context.reactive_nitrogen.previous_total_non_band_ammonium_demand_g_n[layer] else context.reactive_nitrogen.previous_total_band_ammonium_demand_g_n[layer];
                    const previous = if (zone == 0) context.reactive_nitrogen.previous_non_band_microbial_ammonium_capacity_g_n[unit] else context.reactive_nitrogen.previous_band_microbial_ammonium_capacity_g_n[unit];
                    const share = competition(total_previous, previous, active_fraction * ammonium_fraction[zone], context.negligible_amount, context.parameters.nitrification.minimum_competition_fraction);
                    const available_g_n = @max(0, (ammonium_concentration[zone] - p.ammonium_minimum_concentration_g_n_per_m3) * water_m3 * ammonium_fraction[zone]);
                    ammonium_exchange[zone] = @min(ammonium_capacity[zone], share * available_g_n * context.timestep_h);
                }
            } else {
                ammonium_exchange = .{ nitrogen_demand * ammonium_fraction[0], nitrogen_demand * ammonium_fraction[1] };
            }
            const nitrate_demand = @max(0, nitrogen_demand - ammonium_exchange[0] - ammonium_exchange[1]);
            var nitrate_exchange = [2]f64{ 0, 0 };
            var nitrate_capacity = [2]f64{ 0, 0 };
            if (nitrate_demand > 0) {
                const limited_demand = @min(nitrate_demand, p.nitrate_maximum_uptake_g_n_per_m2_h * area_activity);
                inline for (0..2) |zone| {
                    const available_concentration = @max(0, nitrate_concentration[zone] - p.nitrate_minimum_concentration_g_n_per_m3);
                    nitrate_capacity[zone] = nitrate_fraction[zone] * limited_demand * available_concentration / (available_concentration + p.nitrate_half_saturation_g_n_per_m3);
                    const total_previous = if (zone == 0) context.reactive_nitrogen.previous_total_non_band_nitrate_demand_g_n[layer] else context.reactive_nitrogen.previous_total_band_nitrate_demand_g_n[layer];
                    const previous = if (zone == 0) context.reactive_nitrogen.previous_non_band_microbial_nitrate_capacity_g_n[unit] else context.reactive_nitrogen.previous_band_microbial_nitrate_capacity_g_n[unit];
                    const share = competition(total_previous, previous, active_fraction * nitrate_fraction[zone], context.negligible_amount, context.parameters.denitrification.minimum_competition_fraction);
                    const available_g_n = @max(0, (nitrate_concentration[zone] - p.nitrate_minimum_concentration_g_n_per_m3) * water_m3 * nitrate_fraction[zone]);
                    nitrate_exchange[zone] = @min(nitrate_capacity[zone], share * available_g_n * context.timestep_h);
                }
            }
            context.result.non_band_microbial_ammonium_exchange_g_n[unit] = ammonium_exchange[0];
            context.result.band_microbial_ammonium_exchange_g_n[unit] = ammonium_exchange[1];
            context.result.non_band_microbial_nitrate_exchange_g_n[unit] = nitrate_exchange[0];
            context.result.band_microbial_nitrate_exchange_g_n[unit] = nitrate_exchange[1];
            context.result.non_band_microbial_ammonium_capacity_g_n[unit] = ammonium_capacity[0];
            context.result.band_microbial_ammonium_capacity_g_n[unit] = ammonium_capacity[1];
            context.result.non_band_microbial_nitrate_capacity_g_n[unit] = nitrate_capacity[0];
            context.result.band_microbial_nitrate_capacity_g_n[unit] = nitrate_capacity[1];
        };
    }
}

fn competition(total: f64, previous: f64, fallback: f64, negligible: f64, minimum: f64) f64 {
    return @max(minimum, if (total > negligible) previous / total else fallback);
}

fn validate(context: ApplyContext, range: compute.CellRange) !void {
    const layers = context.model_grid.layer_count;
    if (range.first > range.end or range.end > layers or context.reactive_nitrogen.layer_count != layers or context.chemistry_state.cell_count != layers or context.microbial_state.cell_count * context.microbial_state.layer_count != layers or context.result.process_unit_count_per_layer != context.microbial_state.substrate_count * context.microbial_state.population_count or context.matric_plus_osmotic_potential_megapascal.len != layers) return error.SoilMicrobialNitrogenExchangeDimensionMismatch;
    inline for (.{ context.nitrogen_molar_mass_g_per_mol, context.timestep_h }) |value| if (!std.math.isFinite(value) or value <= 0) return error.InvalidSoilMicrobialNitrogenExchangeInput;
    if (!std.math.isFinite(context.negligible_amount) or context.negligible_amount < 0) return error.InvalidSoilMicrobialNitrogenExchangeInput;
}

test "runtime microbial nitrogen exchange immobilizes deficits and mineralizes surplus" {
    const config = @import("../../core/config.zig").SimulationConfig{ .lon_count = 1, .lat_count = 1, .soil_layers = 1, .plant_populations = 1, .worker_threads = 1, .tile_cells = 1, .relative_tolerance = 1e-8, .absolute_tolerance = 1e-12, .mass_balance_tolerance = 1e-12, .negligible_quantity_threshold = 1e-12, .max_nonlinear_iterations = 20, .picard_relaxation = 0.5 };
    var model_grid = try grid.GridState.init(std.testing.allocator, config);
    defer model_grid.deinit();
    model_grid.matrix_liquid_water_m3[0] = 1;
    model_grid.soil_temperature_k[0] = 293.15;
    var microbial_state = try microbial.State.init(std.testing.allocator, 1, 1, 1, 1);
    defer microbial_state.deinit();
    microbial_state.structural[0] = .{ .carbon_g_c = 0.55, .nitrogen_g_n = 0.055, .phosphorus_g_p = 0.0055 };
    microbial_state.nonstructural[0].carbon_g_c = 1;
    var chemistry_state = try chemistry.State.init(std.testing.allocator, 1);
    defer chemistry_state.deinit();
    chemistry_state.aqueous[0].ammonium_non_band = 1;
    chemistry_state.aqueous[0].nitrate_non_band = 1;
    var reactive_state = try reactive.State.init(std.testing.allocator, 1, 1);
    defer reactive_state.deinit();
    var result = try fluxes.State.init(std.testing.allocator, 1, 1);
    defer result.deinit();
    const source =
        "soil_nitrification 0.001 0.0002 7000 14 1.4 1.4 0.125 0.125 0.3 0.1 0.5 2.667 3.429 1.143\n" ++
        "soil_denitrification 0.001 1.4 1.4 0.014 1 0.429 0.429 0.214 0.875\nsoil_autotrophic_denitrification 0.5 0.333\n" ++
        "soil_chemodenitrification 0.0005 0.001 1e-12 0.5 0 0.5\nnitrous_acid_dissociation_mol_per_m3 0.45\nsoil_microbial_thermal_adaptation_offset_k 0\n" ++
        "soil_nitrifier_indices 5 0 1 1\nsoil_nitrifier_environment 0.55 0.1 0.1 0.01 0.01 12 0.1\n" ++
        "soil_oxygen_uptake 1e-6 2.3866348449e11 0.064 -1.5e4 0.5 12 12 0.5 0.7 0.001 1e-12\n" ++
        "soil_heterotrophic_respiration 0.125 0.1 0.01 12 12 0.5 0.42016806722689076 0.1 2.667 0.01 0.01 1e-6 1 0.7142857142857143\n" ++
        "soil_microbial_mineral_exchange 0.014 0.0125 0.40 0.014 0.03 0.35 0.003 0.009 0.18 31\n" ++
        "soil_nonsymbiotic_nitrogen_fixation 5 6 0.25 0.02 0.14 0.25\n" ++
        "soil_microbial_turnover 0.01 0.001 0.167 0.333 0.333 0.333 0.150 0.300 0.333 0.25 2.0 5.0 1.0 0.5";
    const parsed = try nitrogen_parameters.parse(source);
    var context: ApplyContext = .{ .result = &result, .reactive_nitrogen = &reactive_state, .microbial_state = &microbial_state, .chemistry_state = &chemistry_state, .model_grid = &model_grid, .matric_plus_osmotic_potential_megapascal = &.{0}, .zone_fractions = .{ .ammonium_non_band = 1, .ammonium_band = 0, .nitrate_non_band = 1, .nitrate_band = 0, .phosphate_non_band = 1, .phosphate_band = 0 }, .parameters = parsed, .nitrogen_molar_mass_g_per_mol = 14, .timestep_h = 1, .negligible_amount = 1e-12 };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expect(result.non_band_microbial_ammonium_exchange_g_n[0] > 0);
    try std.testing.expect(result.non_band_microbial_ammonium_capacity_g_n[0] > 0);
    microbial_state.nonstructural[0].nitrogen_g_n = 1;
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expect(result.non_band_microbial_ammonium_exchange_g_n[0] < 0);
    try std.testing.expectEqual(@as(f64, 0), result.non_band_microbial_nitrate_exchange_g_n[0]);
}
