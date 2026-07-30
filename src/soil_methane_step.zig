const std = @import("std");
const compute = @import("compute.zig");
const gas = @import("gas_transport.zig");
const methanogenesis = @import("soil_methanogenesis.zig");
const oxidation = @import("soil_methane_oxidation.zig");
const microbial = @import("soil_microbial_state.zig");
const products = @import("soil_respiration_products_step.zig");
const metabolism = @import("soil_microbial_metabolism.zig");

pub const Parameters = struct {
    hydrogenotrophic: methanogenesis.HydrogenotrophicParameters,
    methane_half_saturation_g_c_per_m3: f64,
    methane_solubility_water_to_air: f64,
    gas_exchange_rate_per_step: f64,
    biomass_conversion_efficiency_g_c_per_g_c: f64,
    methanotroph_growth_respiration_g_c_per_g_c: f64,
};

pub const State = struct {
    allocator: std.mem.Allocator,
    layer_count: usize,
    gaseous_methane_after_g_c: []f64,
    aqueous_methane_after_g_c: []f64,
    hydrogenotrophic_methane_g_c: []f64,
    hydrogen_consumption_g_h: []f64,
    methane_oxidation_to_biomass_g_c: []f64,
    methane_oxidation_respiration_g_c: []f64,
    oxygen_demand_g_o: []f64,
    solver_iterations: []u16,
    fermentation_hydrogen_production_g_h: []f64,
    acetotrophic_methane_production_g_c: []f64,
    hydrogenotroph_active_biomass_g_c: []f64,
    methanotroph_active_biomass_g_c: []f64,
    temperature_water_response: []f64,
    nutrient_limitation_fraction: []f64,
    aqueous_co2_limitation_fraction: []f64,
    hydrogen_feedback_energy_kj_per_mol: []f64,

    pub fn init(allocator: std.mem.Allocator, layer_count: usize) !State {
        if (layer_count == 0) return error.InvalidSoilMethaneDimensions;
        var result: State = undefined;
        result.allocator = allocator;
        result.layer_count = layer_count;
        var allocated_f64: usize = 0;
        errdefer inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64 and allocated_f64 > 0) {
            allocated_f64 -= 1;
            allocator.free(@field(result, field.name));
        };
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) {
            @field(result, field.name) = try allocator.alloc(f64, layer_count);
            allocated_f64 += 1;
            @memset(@field(result, field.name), 0);
        };
        result.solver_iterations = try allocator.alloc(u16, layer_count);
        @memset(result.solver_iterations, 0);
        return result;
    }

    pub fn deinit(self: *State) void {
        self.allocator.free(self.solver_iterations);
        inline for (@typeInfo(State).@"struct".fields) |field| if (field.type == []f64) self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

pub const PrepareContext = struct {
    result: *State,
    microbial_state: *const microbial.State,
    respiration_products: *const products.State,
    gas_state: *const gas.State,
    water_volume_m3: []const f64,
    soil_temperature_k: []const f64,
    matric_plus_osmotic_potential_mpa: []const f64,
    autotrophic_substrate_index: usize,
    hydrogenotroph_population_index: usize,
    methanotroph_population_index: usize,
    labile_biomass_fraction: f64,
    target_nitrogen_per_carbon_g_n_per_g_c: f64,
    target_phosphorus_per_carbon_g_p_per_g_c: f64,
    aqueous_co2_half_saturation_g_c_per_m3: f64,
    hydrogen_product_inhibition_g_h_per_m3: f64,
    gas_constant_kj_per_mol_k: f64,
    minimum_hydrogen_concentration_g_h_per_m3: f64,
    hydrogen_feedback_stoichiometric_exponent: f64,
    thermal_adaptation_offset_k: f64,
};

/// Exact NITRO.F 341--347 aqueous CO2 constraint on autotrophic uptake.
pub fn aqueousCo2LimitationFraction(
    aqueous_co2_concentration_g_c_per_m3: f64,
    half_saturation_g_c_per_m3: f64,
) !f64 {
    if (!std.math.isFinite(aqueous_co2_concentration_g_c_per_m3) or
        aqueous_co2_concentration_g_c_per_m3 < 0 or
        !std.math.isFinite(half_saturation_g_c_per_m3) or
        half_saturation_g_c_per_m3 <= 0)
        return error.InvalidAqueousCarbonDioxideLimitation;
    const result = aqueous_co2_concentration_g_c_per_m3 /
        (aqueous_co2_concentration_g_c_per_m3 + half_saturation_g_c_per_m3);
    if (!std.math.isFinite(result) or result < 0 or result > 1)
        return error.NonFiniteAqueousCarbonDioxideLimitation;
    return result;
}

/// Exact NITRO.F 556--557 GH2X operation order.
pub fn hydrogenFeedbackEnergy_kj_per_mol(
    soil_temperature_k: f64,
    aqueous_hydrogen_concentration_g_h_per_m3: f64,
    product_inhibition_g_h_per_m3: f64,
    gas_constant_kj_per_mol_k: f64,
    minimum_hydrogen_concentration_g_h_per_m3: f64,
    stoichiometric_exponent: f64,
) !f64 {
    inline for (.{
        soil_temperature_k,
        aqueous_hydrogen_concentration_g_h_per_m3,
        product_inhibition_g_h_per_m3,
        gas_constant_kj_per_mol_k,
        minimum_hydrogen_concentration_g_h_per_m3,
        stoichiometric_exponent,
    }) |value| if (!std.math.isFinite(value))
        return error.InvalidHydrogenFeedbackEnvironment;
    if (soil_temperature_k <= 0 or
        aqueous_hydrogen_concentration_g_h_per_m3 < 0 or
        product_inhibition_g_h_per_m3 <= 0 or
        gas_constant_kj_per_mol_k <= 0 or
        minimum_hydrogen_concentration_g_h_per_m3 <= 0 or
        stoichiometric_exponent <= 0)
        return error.InvalidHydrogenFeedbackEnvironment;
    const concentration_ratio =
        @max(
            minimum_hydrogen_concentration_g_h_per_m3,
            aqueous_hydrogen_concentration_g_h_per_m3,
        ) / product_inhibition_g_h_per_m3;
    const feedback = gas_constant_kj_per_mol_k * soil_temperature_k *
        @log(std.math.pow(f64, concentration_ratio, stoichiometric_exponent));
    if (!std.math.isFinite(feedback))
        return error.NonFiniteHydrogenFeedbackEnvironment;
    return feedback;
}

pub fn prepareTile(context: *PrepareContext, range: compute.CellRange) !void {
    const n = context.result.layer_count;
    if (range.first > range.end or range.end > n or context.gas_state.cell_count != n or context.respiration_products.layer_count != n or context.water_volume_m3.len != n or context.soil_temperature_k.len != n or context.matric_plus_osmotic_potential_mpa.len != n) return error.InvalidSoilMethaneDimensions;
    if (context.autotrophic_substrate_index >= context.microbial_state.substrate_count or context.hydrogenotroph_population_index >= context.microbial_state.population_count or context.methanotroph_population_index >= context.microbial_state.population_count or context.labile_biomass_fraction <= 0 or context.hydrogen_product_inhibition_g_h_per_m3 <= 0 or !std.math.isFinite(context.aqueous_co2_half_saturation_g_c_per_m3) or context.aqueous_co2_half_saturation_g_c_per_m3 <= 0 or !std.math.isFinite(context.gas_constant_kj_per_mol_k) or context.gas_constant_kj_per_mol_k <= 0 or !std.math.isFinite(context.minimum_hydrogen_concentration_g_h_per_m3) or context.minimum_hydrogen_concentration_g_h_per_m3 <= 0 or !std.math.isFinite(context.hydrogen_feedback_stoichiometric_exponent) or context.hydrogen_feedback_stoichiometric_exponent <= 0) return error.InvalidSoilMethaneInput;
    for (range.first..range.end) |layer| {
        var fermentation_hydrogen: f64 = 0;
        var acetotrophic_methane: f64 = 0;
        const first = layer * context.respiration_products.process_unit_count_per_layer;
        const end = first + context.respiration_products.process_unit_count_per_layer;
        for (first..end) |unit| {
            fermentation_hydrogen += context.respiration_products.hydrogen_g_h[unit];
            acetotrophic_methane += context.respiration_products.methane_g_c[unit];
        }
        const substrate = context.autotrophic_substrate_index;
        const hydrogen_index = try context.microbial_state.populationIndex(layer / context.microbial_state.layer_count, layer % context.microbial_state.layer_count, substrate, context.hydrogenotroph_population_index);
        const methane_index = try context.microbial_state.populationIndex(layer / context.microbial_state.layer_count, layer % context.microbial_state.layer_count, substrate, context.methanotroph_population_index);
        const hydrogen_labile = context.microbial_state.structural[hydrogen_index * 2];
        const methane_labile = context.microbial_state.structural[methane_index * 2];
        const active_hydrogen = hydrogen_labile.carbon_g_c / context.labile_biomass_fraction;
        const active_methane = methane_labile.carbon_g_c / context.labile_biomass_fraction;
        const n_ratio = if (hydrogen_labile.carbon_g_c > 0) hydrogen_labile.nitrogen_g_n / hydrogen_labile.carbon_g_c else context.target_nitrogen_per_carbon_g_n_per_g_c;
        const p_ratio = if (hydrogen_labile.carbon_g_c > 0) hydrogen_labile.phosphorus_g_p / hydrogen_labile.carbon_g_c else context.target_phosphorus_per_carbon_g_p_per_g_c;
        const nutrient = @min(@min(1.0, @max(0.1, std.math.pow(f64, n_ratio / context.target_nitrogen_per_carbon_g_n_per_g_c, 0.25))), @min(1.0, @max(0.1, std.math.pow(f64, p_ratio / context.target_phosphorus_per_carbon_g_p_per_g_c, 0.25))));
        const temperature_water = try metabolism.growthTemperatureResponse(context.soil_temperature_k[layer], context.thermal_adaptation_offset_k) * @exp(0.1 * context.matric_plus_osmotic_potential_mpa[layer]);
        const water_m3 = context.water_volume_m3[layer];
        if (water_m3 <= 0) return error.InvalidSoilMethaneInput;
        const co2 = context.gas_state.dissolved_mass_g[try gas.massIndex(layer, .carbon_dioxide, n)] / water_m3;
        const h2 = context.gas_state.dissolved_mass_g[try gas.massIndex(layer, .hydrogen, n)] / water_m3;
        const co2_limitation = try aqueousCo2LimitationFraction(
            co2,
            context.aqueous_co2_half_saturation_g_c_per_m3,
        );
        const hydrogen_feedback = try hydrogenFeedbackEnergy_kj_per_mol(
            context.soil_temperature_k[layer],
            h2,
            context.hydrogen_product_inhibition_g_h_per_m3,
            context.gas_constant_kj_per_mol_k,
            context.minimum_hydrogen_concentration_g_h_per_m3,
            context.hydrogen_feedback_stoichiometric_exponent,
        );
        inline for (.{ fermentation_hydrogen, acetotrophic_methane, active_hydrogen, active_methane, temperature_water, nutrient, co2_limitation, hydrogen_feedback }) |value|
            if (!std.math.isFinite(value)) return error.InvalidSoilMethaneInput;
        // All values for this runtime layer are valid before any owner field
        // is published, preventing a partial NITRO preparation record.
        context.result.fermentation_hydrogen_production_g_h[layer] = fermentation_hydrogen;
        context.result.acetotrophic_methane_production_g_c[layer] = acetotrophic_methane;
        context.result.hydrogenotroph_active_biomass_g_c[layer] = active_hydrogen;
        context.result.methanotroph_active_biomass_g_c[layer] = active_methane;
        context.result.temperature_water_response[layer] = temperature_water;
        context.result.nutrient_limitation_fraction[layer] = nutrient;
        context.result.aqueous_co2_limitation_fraction[layer] = co2_limitation;
        context.result.hydrogen_feedback_energy_kj_per_mol[layer] = hydrogen_feedback;
    }
}

test "NITRO 341-347 aqueous CO2 limitation preserves source quotient" {
    try std.testing.expectApproxEqAbs(
        @as(f64, 0.8),
        try aqueousCo2LimitationFraction(4, 1),
        1e-15,
    );
    try std.testing.expectEqual(
        @as(f64, 0),
        try aqueousCo2LimitationFraction(0, 1),
    );
}

test "NITRO aqueous CO2 limitation rejects invalid runtime half saturation" {
    try std.testing.expectError(
        error.InvalidAqueousCarbonDioxideLimitation,
        aqueousCo2LimitationFraction(1, 0),
    );
    try std.testing.expectError(
        error.InvalidAqueousCarbonDioxideLimitation,
        aqueousCo2LimitationFraction(std.math.nan(f64), 1),
    );
}

test "NITRO 556-557 hydrogen feedback preserves source operation order" {
    const expected = 8.3143e-3 * 300 *
        @log(std.math.pow(f64, @max(1.0e-3, 2.0) / 1.0, 4.0));
    try std.testing.expectEqual(
        expected,
        try hydrogenFeedbackEnergy_kj_per_mol(
            300,
            2,
            1,
            8.3143e-3,
            1.0e-3,
            4,
        ),
    );
    try std.testing.expectError(
        error.InvalidHydrogenFeedbackEnvironment,
        hydrogenFeedbackEnergy_kj_per_mol(300, -1, 1, 8.3143e-3, 1.0e-3, 4),
    );
}

pub const ApplyContext = struct {
    result: *State,
    gas_state: *const gas.State,
    water_volume_m3: []const f64,
    fermentation_hydrogen_production_g_h: []const f64,
    acetotrophic_methane_production_g_c: []const f64,
    hydrogenotroph_active_biomass_g_c: []const f64,
    methanotroph_active_biomass_g_c: []const f64,
    temperature_water_response: []const f64,
    nutrient_limitation_fraction: []const f64,
    aqueous_co2_limitation_fraction: []const f64,
    hydrogen_feedback_energy_kj_per_mol: []const f64,
    methanotroph_specific_oxidation_per_h: f64,
    parameters: Parameters,
    timestep_h: f64,
    solver_options: oxidation.Options,
};

/// Ports NITRO hydrogenotrophic `VMXA/H2GSX/RGOMP` and replaces its
/// NPH x NPT methane dissolution/oxidation loops with one bounded implicit
/// Newton-Raphson/Picard solve per runtime layer.
pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    try validate(context.*, range);
    for (range.first..range.end) |layer| {
        const water_m3 = context.water_volume_m3[layer];
        const hydrogen_index = try gas.massIndex(layer, .hydrogen, context.gas_state.cell_count);
        const methane_index = try gas.massIndex(layer, .methane, context.gas_state.cell_count);
        const aqueous_hydrogen_g_h = context.gas_state.dissolved_mass_g[hydrogen_index];
        const hydrogen_concentration = aqueous_hydrogen_g_h / water_m3;
        const hydrogen_result = try methanogenesis.hydrogenotrophic(.{
            .aqueous_hydrogen_concentration_g_h_per_m3 = hydrogen_concentration,
            .aqueous_hydrogen_g_h = aqueous_hydrogen_g_h,
            .fermentation_hydrogen_production_g_h = context.fermentation_hydrogen_production_g_h[layer],
            .temperature_water_response = context.temperature_water_response[layer],
            .nutrient_limitation_fraction = context.nutrient_limitation_fraction[layer],
            .aqueous_co2_limitation_fraction = context.aqueous_co2_limitation_fraction[layer],
            .active_biomass_g_c = context.hydrogenotroph_active_biomass_g_c[layer],
            .timestep_h = context.timestep_h,
            .hydrogen_feedback_energy_kj_per_mol = context.hydrogen_feedback_energy_kj_per_mol[layer],
        }, context.parameters.hydrogenotrophic);
        const hydrogen_consumption_g_h = hydrogen_result.co2_reduction_g_c / context.parameters.hydrogenotrophic.hydrogen_supply_conversion_g_c_per_g_h;
        const methane_production = context.acetotrophic_methane_production_g_c[layer] + hydrogen_result.methane_production_g_c;
        const kinetic_maximum_oxidation = context.methanotroph_specific_oxidation_per_h * context.temperature_water_response[layer] * context.nutrient_limitation_fraction[layer] * context.methanotroph_active_biomass_g_c[layer] * context.timestep_h;
        const oxygen_index = try gas.massIndex(layer, .oxygen, context.gas_state.cell_count);
        const oxygen_available_g_o = context.gas_state.gaseous_mass_g[oxygen_index] + context.gas_state.dissolved_mass_g[oxygen_index];
        const respiration_per_oxidation = context.parameters.biomass_conversion_efficiency_g_c_per_g_c * context.parameters.methanotroph_growth_respiration_g_c_per_g_c;
        const oxygen_per_oxidized_carbon_g_o_per_g_c = 5.333 + 2.667 * respiration_per_oxidation;
        const maximum_oxidation = @min(kinetic_maximum_oxidation, oxygen_available_g_o / oxygen_per_oxidized_carbon_g_o_per_g_c);
        const solved = try oxidation.solve(.{
            .gaseous_methane_g_c = context.gas_state.gaseous_mass_g[methane_index],
            .aqueous_methane_g_c = context.gas_state.dissolved_mass_g[methane_index],
            .gaseous_methane_flux_g_c = 0,
            .aqueous_methane_flux_g_c = 0,
            .methanogenesis_g_c = methane_production,
            .water_volume_m3 = water_m3,
            .air_volume_m3 = context.gas_state.air_volume_m3[layer],
            .methane_solubility_water_to_air = context.parameters.methane_solubility_water_to_air,
            .gas_exchange_rate_per_step = context.parameters.gas_exchange_rate_per_step,
            .gas_exchange_enabled = context.gas_state.air_volume_m3[layer] > 0,
            .minimum_gaseous_methane_g_c = context.solver_options.absolute_tolerance_g_c,
            .methane_half_saturation_g_c_per_m3 = context.parameters.methane_half_saturation_g_c_per_m3,
            .maximum_methane_oxidation_g_c = maximum_oxidation,
            .biomass_conversion_efficiency_g_c_per_g_c = context.parameters.biomass_conversion_efficiency_g_c_per_g_c,
            .growth_respiration_g_c_per_g_c = context.parameters.methanotroph_growth_respiration_g_c_per_g_c,
        }, context.solver_options);
        context.result.gaseous_methane_after_g_c[layer] = solved.gaseous_methane_g_c;
        context.result.aqueous_methane_after_g_c[layer] = solved.aqueous_methane_g_c;
        context.result.hydrogenotrophic_methane_g_c[layer] = hydrogen_result.methane_production_g_c;
        context.result.hydrogen_consumption_g_h[layer] = hydrogen_consumption_g_h;
        context.result.methane_oxidation_to_biomass_g_c[layer] = solved.methane_oxidation_g_c;
        context.result.methane_oxidation_respiration_g_c[layer] = solved.growth_respiration_g_c;
        context.result.oxygen_demand_g_o[layer] = solved.oxygen_demand_g_o;
        context.result.solver_iterations[layer] = solved.iterations;
        inline for (.{ hydrogen_consumption_g_h, methane_production, maximum_oxidation, oxygen_available_g_o }) |value| if (!std.math.isFinite(value) or value < 0) return error.NonFiniteSoilMethaneResult;
    }
}

fn validate(context: ApplyContext, range: compute.CellRange) !void {
    const n = context.result.layer_count;
    if (range.first > range.end or range.end > n or context.gas_state.cell_count != n) return error.InvalidSoilMethaneDimensions;
    inline for (.{ context.water_volume_m3, context.fermentation_hydrogen_production_g_h, context.acetotrophic_methane_production_g_c, context.hydrogenotroph_active_biomass_g_c, context.methanotroph_active_biomass_g_c, context.temperature_water_response, context.nutrient_limitation_fraction, context.aqueous_co2_limitation_fraction, context.hydrogen_feedback_energy_kj_per_mol }) |values| if (values.len != n) return error.InvalidSoilMethaneDimensions;
    if (!std.math.isFinite(context.timestep_h) or context.timestep_h <= 0 or !std.math.isFinite(context.methanotroph_specific_oxidation_per_h) or context.methanotroph_specific_oxidation_per_h < 0) return error.InvalidSoilMethaneInput;
}

test "tiled hydrogenotrophic methanogenesis and implicit methanotrophy conserve carbon" {
    var gas_state = try gas.State.init(std.testing.allocator, 1);
    defer gas_state.deinit();
    gas_state.air_volume_m3[0] = 1;
    gas_state.gaseous_mass_g[@intFromEnum(gas.Species.methane)] = 1;
    gas_state.dissolved_mass_g[@intFromEnum(gas.Species.methane)] = 0.2;
    gas_state.dissolved_mass_g[@intFromEnum(gas.Species.hydrogen)] = 0.1;
    gas_state.dissolved_mass_g[@intFromEnum(gas.Species.oxygen)] = 1;
    var state = try State.init(std.testing.allocator, 1);
    defer state.deinit();
    var context: ApplyContext = .{
        .result = &state,
        .gas_state = &gas_state,
        .water_volume_m3 = &.{1},
        .fermentation_hydrogen_production_g_h = &.{0.01},
        .acetotrophic_methane_production_g_c = &.{0.05},
        .hydrogenotroph_active_biomass_g_c = &.{1},
        .methanotroph_active_biomass_g_c = &.{1},
        .temperature_water_response = &.{1},
        .nutrient_limitation_fraction = &.{1},
        .aqueous_co2_limitation_fraction = &.{1},
        .hydrogen_feedback_energy_kj_per_mol = &.{0},
        .methanotroph_specific_oxidation_per_h = 0.1,
        .parameters = .{
            .hydrogenotrophic = .{ .hydrogen_half_saturation_g_h_per_m3 = 0.01, .specific_co2_reduction_g_c_per_g_c_h = 0.1, .reference_energy_yield_kj_per_g_c = 1, .growth_energy_requirement_kj_per_g_c = 10, .minimum_growth_respiration_fraction = 0.05, .hydrogen_supply_conversion_g_c_per_g_h = 1.5, .fermentation_hydrogen_to_pool_fraction = 0.111 },
            .methane_half_saturation_g_c_per_m3 = 0.01,
            .methane_solubility_water_to_air = 0.03,
            .gas_exchange_rate_per_step = 0.5,
            .biomass_conversion_efficiency_g_c_per_g_c = 0.4,
            .methanotroph_growth_respiration_g_c_per_g_c = 0.5,
        },
        .timestep_h = 1,
        .solver_options = .{ .absolute_tolerance_g_c = 1e-12, .relative_tolerance = 1e-10, .derivative_floor = 1e-14, .picard_relaxation = 0.5, .gas_max_iterations = 80 },
    };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    const methane_before = 1.2 + 0.05 + state.hydrogenotrophic_methane_g_c[0];
    const methane_after = state.gaseous_methane_after_g_c[0] + state.aqueous_methane_after_g_c[0] + state.methane_oxidation_to_biomass_g_c[0] + state.methane_oxidation_respiration_g_c[0];
    try std.testing.expectApproxEqAbs(methane_before, methane_after, 1e-9);
    try std.testing.expect(state.solver_iterations[0] < 80);
    try std.testing.expect(state.hydrogenotrophic_methane_g_c[0] > 0);
    try std.testing.expect(state.oxygen_demand_g_o[0] > 0);
}
