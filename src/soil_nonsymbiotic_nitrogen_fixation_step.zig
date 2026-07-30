const std = @import("std");
const compute = @import("compute.zig");
const gas = @import("gas_transport.zig");
const microbial = @import("soil_microbial_state.zig");
const metabolism = @import("soil_microbial_metabolism.zig");
const fluxes = @import("soil_nitrogen_flux_workspace.zig");
const parameters_module = @import("soil_nitrogen_parameters.zig");

pub const ApplyContext = struct {
    result: *fluxes.State,
    microbial_state: *const microbial.State,
    gas_state: *const gas.State,
    water_volume_m3: []const f64,
    parameters: parameters_module.Parameters,
    timestep_h: f64,
};

/// NITRO RGN2P/RGN2F/RN2FX for user-selected aerobic and anaerobic
/// diazotroph populations in every runtime heterotrophic substrate.
pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    try validate(context.*, range);
    const population_count = context.microbial_state.population_count;
    const substrate_count = context.microbial_state.substrate_count;
    const p = context.parameters.nonsymbiotic_nitrogen_fixation;
    for (range.first..range.end) |layer| {
        const dinitrogen_index = try gas.massIndex(layer, .nitrogen, context.gas_state.cell_count);
        const dinitrogen_concentration = context.gas_state.dissolved_mass_g[dinitrogen_index] / context.water_volume_m3[layer];
        for (0..substrate_count) |substrate| for (0..population_count) |population| {
            const runtime_index = try context.microbial_state.populationIndex(layer / context.microbial_state.layer_count, layer % context.microbial_state.layer_count, substrate, population);
            const unit = layer * context.result.process_unit_count_per_layer + substrate * population_count + population;
            const is_aerobic = population == p.aerobic_diazotroph_population_index;
            const is_anaerobic = population == p.anaerobic_diazotroph_population_index;
            const value = try metabolism.nonsymbioticNitrogenFixation(.{
                .is_diazotroph = is_aerobic or is_anaerobic,
                .nonstructural_carbon_g_c = context.microbial_state.nonstructural[runtime_index].carbon_g_c,
                .nonstructural_nitrogen_g_n = context.microbial_state.nonstructural[runtime_index].nitrogen_g_n,
                .maximum_nitrogen_per_carbon_g_n_per_g_c = context.parameters.heterotrophic_respiration.target_nitrogen_per_carbon_g_n_per_g_c,
                .nitrogen_fixation_yield_g_n_per_g_c = if (is_aerobic) p.aerobic_yield_g_n_per_g_c else p.anaerobic_yield_g_n_per_g_c,
                .growth_respiration_g_c = context.result.growth_respiration_g_c[unit],
                .aqueous_dinitrogen_concentration_g_n_per_m3 = dinitrogen_concentration,
                .dinitrogen_half_saturation_g_n_per_m3 = p.dinitrogen_half_saturation_g_n_per_m3,
                .nonstructural_to_structural_rate_per_h = p.nonstructural_to_structural_rate_per_h,
                .timestep_h = context.timestep_h,
            });
            context.result.nitrogen_fixation_respiration_g_c[unit] = value.fixation_respiration_g_c;
            context.result.fixed_dinitrogen_g_n[unit] = value.fixed_nitrogen_g_n;
        };
    }
}

fn validate(context: ApplyContext, range: compute.CellRange) !void {
    const layers = context.gas_state.cell_count;
    const p = context.parameters.nonsymbiotic_nitrogen_fixation;
    if (range.first > range.end or range.end > layers or context.microbial_state.cell_count * context.microbial_state.layer_count != layers or context.result.layer_count != layers or context.result.process_unit_count_per_layer != context.microbial_state.substrate_count * context.microbial_state.population_count or context.water_volume_m3.len != layers) return error.SoilNitrogenFixationDimensionMismatch;
    if (p.aerobic_diazotroph_population_index >= context.microbial_state.population_count or p.anaerobic_diazotroph_population_index >= context.microbial_state.population_count) return error.SoilDiazotrophPopulationIndexOutOfBounds;
    if (!std.math.isFinite(context.timestep_h) or context.timestep_h <= 0) return error.InvalidSoilNitrogenFixationInput;
    for (range.first..range.end) |layer| if (!std.math.isFinite(context.water_volume_m3[layer]) or context.water_volume_m3[layer] <= 0) return error.InvalidSoilNitrogenFixationInput;
}

test "runtime soil diazotroph roles fix dissolved dinitrogen" {
    var microbial_state = try microbial.State.init(std.testing.allocator, 1, 1, 1, 7);
    defer microbial_state.deinit();
    for (0..7) |population| microbial_state.nonstructural[population].carbon_g_c = 1;
    var gas_state = try gas.State.init(std.testing.allocator, 1);
    defer gas_state.deinit();
    gas_state.dissolved_mass_g[try gas.massIndex(0, .nitrogen, 1)] = 1;
    var result = try fluxes.State.init(std.testing.allocator, 1, 7);
    defer result.deinit();
    result.growth_respiration_g_c[5] = 0.5;
    result.growth_respiration_g_c[6] = 0.5;
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
    var context: ApplyContext = .{ .result = &result, .microbial_state = &microbial_state, .gas_state = &gas_state, .water_volume_m3 = &.{1}, .parameters = try parameters_module.parse(source), .timestep_h = 1 };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expect(result.fixed_dinitrogen_g_n[5] > 0);
    try std.testing.expect(result.fixed_dinitrogen_g_n[6] > 0);
    try std.testing.expect(result.fixed_dinitrogen_g_n[5] > result.fixed_dinitrogen_g_n[6]);
    try std.testing.expectEqual(@as(f64, 0), result.fixed_dinitrogen_g_n[0]);
    try std.testing.expectApproxEqAbs(result.nitrogen_fixation_respiration_g_c[5] * 0.25, result.fixed_dinitrogen_g_n[5], 1e-15);
}
