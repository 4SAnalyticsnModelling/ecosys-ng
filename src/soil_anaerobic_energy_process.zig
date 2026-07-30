const std = @import("std");
const nitrogen_parameters = @import("soil_nitrogen_parameters.zig");
const energy_workspace = @import("soil_anaerobic_energy_workspace.zig");
const fermenter = @import("soil_fermenter_respiration.zig");
const acetotrophic = @import("soil_acetotrophic_methanogenesis.zig");

/// Read-only scientific state required by NITRO.F 900--1048. Every slice is
/// indexed by the same runtime microbial-unit axis as the parameter workspace.
pub const Inputs = struct {
    soil_temperature_k: []const f64,
    aqueous_acetate_concentration_g_c_per_m3: []const f64,
    aqueous_acetate_g_c: []const f64,
    acetate_competition_fraction: []const f64,
    aqueous_dissolved_organic_carbon_concentration_g_c_per_m3: []const f64,
    aqueous_dissolved_organic_carbon_g_c: []const f64,
    dissolved_organic_carbon_competition_fraction: []const f64,
    combined_nutrient_limitation_fraction: []const f64,
    water_response: []const f64,
    active_biomass_g_c: []const f64,
    growth_temperature_response: []const f64,
    fermentation_oxygen_inhibition_fraction: []const f64,
    hydrogen_feedback_energy_kj_per_mol: []const f64,
    minimum_acetate_concentration_g_c_per_m3: f64,
    timestep_h: f64,
};

/// Commitless paired results. Authoritative organic, gas, and competition
/// owners consume these values only after both kernels succeed.
pub const State = struct {
    allocator: std.mem.Allocator,
    unit_count: usize,
    fermenter: fermenter.State,
    acetotrophic: acetotrophic.State,

    pub fn init(allocator: std.mem.Allocator, unit_count: usize) !State {
        if (unit_count == 0) return error.InvalidAnaerobicEnergyProcessDimensions;
        var fermenter_state = try fermenter.State.init(allocator, unit_count);
        errdefer fermenter_state.deinit();
        const acetotrophic_state = try acetotrophic.State.init(allocator, unit_count);
        return .{
            .allocator = allocator,
            .unit_count = unit_count,
            .fermenter = fermenter_state,
            .acetotrophic = acetotrophic_state,
        };
    }

    pub fn deinit(self: *State) void {
        self.fermenter.deinit();
        self.acetotrophic.deinit();
        self.* = undefined;
    }
};

/// Runs both anaerobic branches atomically into a commitless process state.
pub fn calculate(
    state: *State,
    workspace: *energy_workspace.Workspace,
    parameters: nitrogen_parameters.AnaerobicEnergyParameters,
    inputs: Inputs,
) !void {
    try validateDimensions(state, workspace, inputs);
    try workspace.populate(parameters);
    var staged = try State.init(state.allocator, state.unit_count);
    errdefer staged.deinit();

    try fermenter.calculate(&staged.fermenter, .{
        .enabled = workspace.fermenter_enabled,
        .aqueous_acetate_concentration_g_c_per_m3 = inputs.aqueous_acetate_concentration_g_c_per_m3,
        .aqueous_dissolved_organic_carbon_concentration_g_c_per_m3 = inputs.aqueous_dissolved_organic_carbon_concentration_g_c_per_m3,
        .aqueous_dissolved_organic_carbon_g_c = inputs.aqueous_dissolved_organic_carbon_g_c,
        .dissolved_organic_carbon_competition_fraction = inputs.dissolved_organic_carbon_competition_fraction,
        .combined_nutrient_limitation_fraction = inputs.combined_nutrient_limitation_fraction,
        .water_response = inputs.water_response,
        .active_biomass_g_c = inputs.active_biomass_g_c,
        .growth_temperature_response = inputs.growth_temperature_response,
        .fermentation_oxygen_inhibition_fraction = inputs.fermentation_oxygen_inhibition_fraction,
        .soil_temperature_k = inputs.soil_temperature_k,
        .hydrogen_feedback_energy_kj_per_mol = inputs.hydrogen_feedback_energy_kj_per_mol,
        .reference_fermentation_energy_yield_kj_per_g_c = workspace.fermenter_reference_energy_yield_kj_per_g_c,
        .growth_energy_requirement_kj_per_g_c = workspace.fermenter_growth_energy_requirement_kj_per_g_c,
        .minimum_respiration_requirement_g_c_per_g_c = workspace.fermenter_minimum_respiration_requirement_g_c_per_g_c,
        .specific_oxidation_rate_g_c_per_g_c_h = workspace.fermenter_specific_oxidation_rate_g_c_per_g_c_h,
        .dissolved_organic_carbon_half_saturation_g_c_per_m3 = workspace.fermenter_doc_half_saturation_g_c_per_m3,
        .acetate_product_inhibition_g_c_per_m3 = workspace.fermenter_acetate_product_inhibition_g_c_per_m3,
        .minimum_acetate_concentration_g_c_per_m3 = inputs.minimum_acetate_concentration_g_c_per_m3,
        .gas_constant_kj_per_mol_k = parameters.fermenter.gas_constant_kj_per_mol_k,
        .acetate_feedback_stoichiometric_exponent = parameters.fermenter.acetate_feedback_stoichiometric_exponent,
        .feedback_carbon_conversion_g_c_per_mol = parameters.fermenter.feedback_carbon_conversion_g_c_per_mol,
        .timestep_h = inputs.timestep_h,
    });
    try acetotrophic.calculate(&staged.acetotrophic, .{
        .enabled = workspace.acetotrophic_enabled,
        .soil_temperature_k = inputs.soil_temperature_k,
        .aqueous_acetate_concentration_g_c_per_m3 = inputs.aqueous_acetate_concentration_g_c_per_m3,
        .aqueous_acetate_g_c = inputs.aqueous_acetate_g_c,
        .acetate_competition_fraction = inputs.acetate_competition_fraction,
        .combined_nutrient_limitation_fraction = inputs.combined_nutrient_limitation_fraction,
        .water_response = inputs.water_response,
        .active_biomass_g_c = inputs.active_biomass_g_c,
        .growth_temperature_response = inputs.growth_temperature_response,
        .acetate_product_inhibition_g_c_per_m3 = workspace.acetotrophic_acetate_product_inhibition_g_c_per_m3,
        .acetate_half_saturation_g_c_per_m3 = workspace.acetotrophic_acetate_half_saturation_g_c_per_m3,
        .specific_respiration_rate_g_c_per_g_c_h = workspace.acetotrophic_specific_respiration_rate_g_c_per_g_c_h,
        .reference_energy_yield_kj_per_g_c = workspace.acetotrophic_reference_energy_yield_kj_per_g_c,
        .growth_energy_requirement_kj_per_g_c = workspace.acetotrophic_growth_energy_requirement_kj_per_g_c,
        .minimum_growth_respiration_fraction = workspace.acetotrophic_minimum_growth_respiration_fraction,
        .minimum_acetate_concentration_g_c_per_m3 = inputs.minimum_acetate_concentration_g_c_per_m3,
        .gas_constant_kj_per_mol_k = parameters.acetotrophic_methanogenesis.gas_constant_kj_per_mol_k,
        .feedback_carbon_conversion_g_c_per_mol = parameters.acetotrophic_methanogenesis.feedback_carbon_conversion_g_c_per_mol,
        .methane_carbon_yield_g_c_per_g_c_oxidized = parameters.acetotrophic_methanogenesis.methane_carbon_yield_g_c_per_g_c_oxidized,
        .timestep_h = inputs.timestep_h,
    });

    state.deinit();
    state.* = staged;
}

fn validateDimensions(
    state: *const State,
    workspace: *const energy_workspace.Workspace,
    inputs: Inputs,
) !void {
    const count = state.unit_count;
    if (count == 0 or workspace.unit_count != count)
        return error.InvalidAnaerobicEnergyProcessDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (field.type == []const f64 and @field(inputs, field.name).len != count)
            return error.InvalidAnaerobicEnergyProcessDimensions;
    }
}

test "runtime contract drives each source role and leaves other units inactive" {
    const populations = [_]usize{ 3, 4, 6, 18 };
    var workspace = try energy_workspace.Workspace.init(std.testing.allocator, &populations);
    defer workspace.deinit();
    var state = try State.init(std.testing.allocator, populations.len);
    defer state.deinit();
    const ones = [_]f64{ 1, 1, 1, 1 };
    const temperatures = [_]f64{ 300, 300, 300, 300 };
    const acetate = [_]f64{ 12, 12, 12, 12 };
    const carbon = [_]f64{ 10, 10, 10, 10 };
    const hydrogen = [_]f64{ 0, 0, 0, 0 };
    try calculate(&state, &workspace, try nitrogen_parameters.sourceAnaerobicEnergyParameters(), .{
        .soil_temperature_k = &temperatures,
        .aqueous_acetate_concentration_g_c_per_m3 = &acetate,
        .aqueous_acetate_g_c = &carbon,
        .acetate_competition_fraction = &ones,
        .aqueous_dissolved_organic_carbon_concentration_g_c_per_m3 = &acetate,
        .aqueous_dissolved_organic_carbon_g_c = &carbon,
        .dissolved_organic_carbon_competition_fraction = &ones,
        .combined_nutrient_limitation_fraction = &ones,
        .water_response = &ones,
        .active_biomass_g_c = &ones,
        .growth_temperature_response = &ones,
        .fermentation_oxygen_inhibition_fraction = &ones,
        .hydrogen_feedback_energy_kj_per_mol = &hydrogen,
        .minimum_acetate_concentration_g_c_per_m3 = 1e-12,
        .timestep_h = 1,
    });
    try std.testing.expect(state.fermenter.actual_respiration_g_c[0] > 0);
    try std.testing.expect(state.acetotrophic.acetate_oxidation_g_c[1] > 0);
    try std.testing.expect(state.fermenter.actual_respiration_g_c[2] > 0);
    try std.testing.expectEqual(@as(f64, 0), state.fermenter.actual_respiration_g_c[3]);
    try std.testing.expectEqual(@as(f64, 0), state.acetotrophic.acetate_oxidation_g_c[3]);
}

test "late acetotrophic failure leaves paired process state unchanged" {
    const populations = [_]usize{ 3, 4 };
    var workspace = try energy_workspace.Workspace.init(std.testing.allocator, &populations);
    defer workspace.deinit();
    var state = try State.init(std.testing.allocator, populations.len);
    defer state.deinit();
    state.fermenter.actual_respiration_g_c[0] = 7;
    state.acetotrophic.acetate_oxidation_g_c[1] = 9;
    const ones = [_]f64{ 1, 1 };
    const temperatures = [_]f64{ 300, 300 };
    const acetate = [_]f64{ 12, 12 };
    const carbon = [_]f64{ 10, 10 };
    const hydrogen = [_]f64{ 0, 0 };
    const invalid_competition = [_]f64{ 1, -1 };
    try std.testing.expectError(
        error.InvalidAcetotrophicMethanogenesisInput,
        calculate(&state, &workspace, try nitrogen_parameters.sourceAnaerobicEnergyParameters(), .{
            .soil_temperature_k = &temperatures,
            .aqueous_acetate_concentration_g_c_per_m3 = &acetate,
            .aqueous_acetate_g_c = &carbon,
            .acetate_competition_fraction = &invalid_competition,
            .aqueous_dissolved_organic_carbon_concentration_g_c_per_m3 = &acetate,
            .aqueous_dissolved_organic_carbon_g_c = &carbon,
            .dissolved_organic_carbon_competition_fraction = &ones,
            .combined_nutrient_limitation_fraction = &ones,
            .water_response = &ones,
            .active_biomass_g_c = &ones,
            .growth_temperature_response = &ones,
            .fermentation_oxygen_inhibition_fraction = &ones,
            .hydrogen_feedback_energy_kj_per_mol = &hydrogen,
            .minimum_acetate_concentration_g_c_per_m3 = 1e-12,
            .timestep_h = 1,
        }),
    );
    try std.testing.expectEqual(@as(f64, 7), state.fermenter.actual_respiration_g_c[0]);
    try std.testing.expectEqual(@as(f64, 9), state.acetotrophic.acetate_oxidation_g_c[1]);
}

test "runtime contract rejects a short scientific input array" {
    const populations = [_]usize{ 3, 4 };
    var workspace = try energy_workspace.Workspace.init(std.testing.allocator, &populations);
    defer workspace.deinit();
    var state = try State.init(std.testing.allocator, populations.len);
    defer state.deinit();
    const full = [_]f64{ 1, 1 };
    const short = [_]f64{1};
    try std.testing.expectError(
        error.InvalidAnaerobicEnergyProcessDimensions,
        calculate(&state, &workspace, try nitrogen_parameters.sourceAnaerobicEnergyParameters(), .{
            .soil_temperature_k = &short,
            .aqueous_acetate_concentration_g_c_per_m3 = &full,
            .aqueous_acetate_g_c = &full,
            .acetate_competition_fraction = &full,
            .aqueous_dissolved_organic_carbon_concentration_g_c_per_m3 = &full,
            .aqueous_dissolved_organic_carbon_g_c = &full,
            .dissolved_organic_carbon_competition_fraction = &full,
            .combined_nutrient_limitation_fraction = &full,
            .water_response = &full,
            .active_biomass_g_c = &full,
            .growth_temperature_response = &full,
            .fermentation_oxygen_inhibition_fraction = &full,
            .hydrogen_feedback_energy_kj_per_mol = &full,
            .minimum_acetate_concentration_g_c_per_m3 = 1e-12,
            .timestep_h = 1,
        }),
    );
}
