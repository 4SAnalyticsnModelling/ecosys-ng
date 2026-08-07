const std = @import("std");
const compute = @import("../../core/compute.zig");
const autotrophic_carbon = @import("../microbial/autotrophic_carbon_step.zig");
const methane = @import("methane_step.zig");
const nitrogen_fluxes = @import("../nutrients/nitrogen_flux_workspace.zig");
const oxygen = @import("oxygen_allocation.zig");
const respiration_products = @import("../microbial/respiration_products_step.zig");

/// Layer totals accumulated from the runtime-sized microbial populations and
/// substrate complexes before NITRO publishes its REDIST gas-flux arrays.
pub const Inputs = struct {
    autotrophic_carbon_dioxide_uptake_g_c: f64,
    aerobic_heterotrophic_carbon_dioxide_emission_g_c: f64,
    denitrification_carbon_dioxide_emission_g_c: f64,
    methane_oxidation_g_c: f64,
    methanotrophic_carbon_uptake_g_c: f64,
    methane_emission_g_c: f64,
    hydrogenotrophic_hydrogen_uptake_g_h: f64,
    fermentative_hydrogen_emission_g_h: f64,
    oxygen_uptake_g_o: f64,
    nitrous_oxide_reduction_g_n: f64,
    chemodenitrification_dinitrogen_production_g_n: f64,
    biological_nitrite_reduction_g_n: f64,
    chemodenitrification_nitrous_oxide_production_g_n: f64,
};

pub const Result = struct {
    /// Source `RCO2O`: positive is net biological CO2 uptake.
    source_signed_net_carbon_dioxide_uptake_g_c: f64,
    /// Source `RCH4O`: positive is net biological CH4 uptake.
    source_signed_net_methane_uptake_g_c: f64,
    /// Source `RH2GO`: positive is net biological H2 uptake.
    source_signed_net_hydrogen_uptake_g_h: f64,
    /// Source `RUPOXO`: positive oxygen removal from the gas/water domain.
    oxygen_uptake_g_o: f64,
    /// Source `RN2G`: production is negative under the transport flux sign.
    source_signed_dinitrogen_flux_g_n: f64,
    /// Source `RN2O`: production is negative and reduction is positive.
    source_signed_nitrous_oxide_flux_g_n: f64,
};

/// Runtime-sized, heap-backed transport fluxes. The structure-of-arrays layout
/// keeps each gas contiguous for tiled CPU kernels and a future GPU backend.
pub const State = struct {
    allocator: std.mem.Allocator,
    layer_count: usize,
    source_signed_net_carbon_dioxide_uptake_g_c: []f64,
    source_signed_net_methane_uptake_g_c: []f64,
    source_signed_net_hydrogen_uptake_g_h: []f64,
    oxygen_uptake_g_o: []f64,
    source_signed_dinitrogen_flux_g_n: []f64,
    source_signed_nitrous_oxide_flux_g_n: []f64,

    pub fn init(allocator: std.mem.Allocator, layer_count: usize) !State {
        if (layer_count == 0) return error.InvalidSoilLayerCount;
        var state: State = undefined;
        state.allocator = allocator;
        state.layer_count = layer_count;
        inline for (@typeInfo(State).@"struct".fields[2..]) |field| {
            @field(state, field.name) = try allocator.alloc(f64, layer_count);
            errdefer inline for (@typeInfo(State).@"struct".fields[2..]) |allocated_field| {
                if (@offsetOf(State, allocated_field.name) < @offsetOf(State, field.name))
                    allocator.free(@field(state, allocated_field.name));
            };
            @memset(@field(state, field.name), 0);
        }
        return state;
    }

    pub fn deinit(self: *State) void {
        inline for (@typeInfo(State).@"struct".fields[2..]) |field|
            self.allocator.free(@field(self, field.name));
        self.* = undefined;
    }
};

pub const ApplyContext = struct {
    layer_inputs: []const Inputs,
    result: *State,
};

/// Applies independent layers over a scheduler-provided tile. Inputs are
/// assembled from the prognostic process states by the owning model kernel.
pub fn applyTile(context: *ApplyContext, range: compute.CellRange) !void {
    if (context.layer_inputs.len != context.result.layer_count or range.end > context.result.layer_count or range.first > range.end)
        return error.InvalidSoilGasAggregationDimensions;
    for (range.first..range.end) |layer| {
        const result = try aggregate(context.layer_inputs[layer]);
        inline for (@typeInfo(Result).@"struct".fields) |field|
            @field(context.result, field.name)[layer] = @field(result, field.name);
    }
}

/// Runtime process-state binding for the NITRO label 640/645 aggregation.
/// It publishes diagnostics only; prognostic gas pools remain owned by the
/// nitrogen and autotrophic-carbon commit kernels.
pub const ProcessContext = struct {
    result: *State,
    autotrophic_carbon: *const autotrophic_carbon.State,
    respiration_products: *const respiration_products.State,
    methane: ?*const methane.State,
    oxygen: *const oxygen.State,
    nitrogen_fluxes: *const nitrogen_fluxes.State,
    redox_satisfaction_fraction: []const f64,
};

pub fn aggregateProcessTile(context: *ProcessContext, range: compute.CellRange) !void {
    const layers = context.result.layer_count;
    const units_per_layer = context.nitrogen_fluxes.process_unit_count_per_layer;
    if (range.first > range.end or range.end > layers or
        context.autotrophic_carbon.layer_count != layers or
        context.respiration_products.layer_count != layers or
        context.oxygen.cell_count * context.oxygen.layer_count != layers or
        context.nitrogen_fluxes.layer_count != layers or
        context.autotrophic_carbon.process_unit_count_per_layer != units_per_layer or
        context.respiration_products.process_unit_count_per_layer != units_per_layer or
        context.oxygen.population_count != units_per_layer or
        context.redox_satisfaction_fraction.len != layers * units_per_layer or
        (context.methane != null and context.methane.?.layer_count != layers))
        return error.InvalidSoilGasAggregationDimensions;

    for (range.first..range.end) |layer| {
        const first = layer * units_per_layer;
        const end = first + units_per_layer;
        var inputs: Inputs = std.mem.zeroes(Inputs);
        for (first..end) |unit| {
            const redox = context.redox_satisfaction_fraction[unit];
            if (!std.math.isFinite(redox) or redox < 0 or redox > 1)
                return error.InvalidSoilGasAggregationSatisfactionFraction;
            inputs.autotrophic_carbon_dioxide_uptake_g_c += context.autotrophic_carbon.carbon_dioxide_uptake_g_c[unit];
            inputs.aerobic_heterotrophic_carbon_dioxide_emission_g_c += context.respiration_products.carbon_dioxide_g_c[unit];
            inputs.denitrification_carbon_dioxide_emission_g_c += context.nitrogen_fluxes.denitrification_respiration_g_c[unit] * redox;
            inputs.oxygen_uptake_g_o += context.oxygen.oxygen_uptake_g_o[unit];
            inputs.nitrous_oxide_reduction_g_n += context.nitrogen_fluxes.nitrous_oxide_reduction_potential_g_n[unit] * redox;
            inputs.biological_nitrite_reduction_g_n +=
                (context.nitrogen_fluxes.non_band_heterotrophic_nitrite_reduction_potential_g_n[unit] +
                    context.nitrogen_fluxes.band_heterotrophic_nitrite_reduction_potential_g_n[unit] +
                    context.nitrogen_fluxes.non_band_autotrophic_nitrite_reduction_potential_g_n[unit] +
                    context.nitrogen_fluxes.band_autotrophic_nitrite_reduction_potential_g_n[unit]) * redox;
            inputs.methane_emission_g_c += context.respiration_products.methane_g_c[unit];
            inputs.fermentative_hydrogen_emission_g_h += context.respiration_products.hydrogen_g_h[unit];
        }
        inputs.chemodenitrification_dinitrogen_production_g_n =
            context.nitrogen_fluxes.chemodenitrification_dinitrogen_production_g_n[layer];
        inputs.chemodenitrification_nitrous_oxide_production_g_n =
            context.nitrogen_fluxes.chemodenitrification_nitrous_oxide_production_g_n[layer];
        if (context.methane) |methane_state| {
            inputs.methane_oxidation_g_c = methane_state.methane_oxidation_respiration_g_c[layer];
            inputs.methanotrophic_carbon_uptake_g_c = methane_state.methane_oxidation_to_biomass_g_c[layer];
            inputs.methane_emission_g_c += methane_state.hydrogenotrophic_methane_g_c[layer];
            inputs.hydrogenotrophic_hydrogen_uptake_g_h = methane_state.hydrogen_consumption_g_h[layer];
        }
        const result = try aggregate(inputs);
        inline for (@typeInfo(Result).@"struct".fields) |field|
            @field(context.result, field.name)[layer] = @field(result, field.name);
    }
}

/// Exact NITRO aggregation following labels 640/645:
/// RCO2O, RCH4O, RH2GO, RUPOXO, RN2G, and RN2O.
pub fn aggregate(inputs: Inputs) !Result {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value)) return error.NonFiniteSoilBiogeochemicalGasFlux;
        if (value < 0) return error.NegativeSoilBiogeochemicalGasComponent;
    }
    const result: Result = .{
        .source_signed_net_carbon_dioxide_uptake_g_c = inputs.autotrophic_carbon_dioxide_uptake_g_c -
            inputs.aerobic_heterotrophic_carbon_dioxide_emission_g_c -
            inputs.denitrification_carbon_dioxide_emission_g_c -
            inputs.methane_oxidation_g_c,
        .source_signed_net_methane_uptake_g_c = inputs.methane_oxidation_g_c +
            inputs.methanotrophic_carbon_uptake_g_c -
            inputs.methane_emission_g_c,
        .source_signed_net_hydrogen_uptake_g_h = inputs.hydrogenotrophic_hydrogen_uptake_g_h -
            inputs.fermentative_hydrogen_emission_g_h,
        .oxygen_uptake_g_o = inputs.oxygen_uptake_g_o,
        .source_signed_dinitrogen_flux_g_n = -inputs.nitrous_oxide_reduction_g_n -
            inputs.chemodenitrification_dinitrogen_production_g_n,
        .source_signed_nitrous_oxide_flux_g_n = -inputs.biological_nitrite_reduction_g_n -
            inputs.chemodenitrification_nitrous_oxide_production_g_n +
            inputs.nitrous_oxide_reduction_g_n,
    };
    inline for (@typeInfo(Result).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(result, field.name))) return error.NonFiniteSoilBiogeochemicalGasFlux;
    }
    return result;
}

test "NITRO layer gas aggregation preserves exact source signs and terms" {
    const result = try aggregate(.{
        .autotrophic_carbon_dioxide_uptake_g_c = 10,
        .aerobic_heterotrophic_carbon_dioxide_emission_g_c = 2,
        .denitrification_carbon_dioxide_emission_g_c = 1,
        .methane_oxidation_g_c = 0.5,
        .methanotrophic_carbon_uptake_g_c = 0.25,
        .methane_emission_g_c = 1.5,
        .hydrogenotrophic_hydrogen_uptake_g_h = 0.8,
        .fermentative_hydrogen_emission_g_h = 0.3,
        .oxygen_uptake_g_o = 4,
        .nitrous_oxide_reduction_g_n = 0.6,
        .chemodenitrification_dinitrogen_production_g_n = 0.1,
        .biological_nitrite_reduction_g_n = 1.2,
        .chemodenitrification_nitrous_oxide_production_g_n = 0.2,
    });
    try std.testing.expectApproxEqAbs(@as(f64, 6.5), result.source_signed_net_carbon_dioxide_uptake_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.75), result.source_signed_net_methane_uptake_g_c, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), result.source_signed_net_hydrogen_uptake_g_h, 1e-15);
    try std.testing.expectEqual(@as(f64, 4), result.oxygen_uptake_g_o);
    try std.testing.expectApproxEqAbs(@as(f64, -0.7), result.source_signed_dinitrogen_flux_g_n, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -0.8), result.source_signed_nitrous_oxide_flux_g_n, 1e-15);
}

test "gas aggregation rejects invalid components instead of propagating them" {
    const zero: Inputs = std.mem.zeroes(Inputs);
    var invalid = zero;
    invalid.methane_emission_g_c = std.math.nan(f64);
    try std.testing.expectError(error.NonFiniteSoilBiogeochemicalGasFlux, aggregate(invalid));
    invalid = zero;
    invalid.oxygen_uptake_g_o = -1;
    try std.testing.expectError(error.NegativeSoilBiogeochemicalGasComponent, aggregate(invalid));
}

test "runtime-sized tiled gas aggregation is heap backed and layer independent" {
    const allocator = std.testing.allocator;
    var state = try State.init(allocator, 2);
    defer state.deinit();
    const inputs = [_]Inputs{
        .{ .autotrophic_carbon_dioxide_uptake_g_c = 3, .aerobic_heterotrophic_carbon_dioxide_emission_g_c = 1, .denitrification_carbon_dioxide_emission_g_c = 0, .methane_oxidation_g_c = 0, .methanotrophic_carbon_uptake_g_c = 0, .methane_emission_g_c = 0, .hydrogenotrophic_hydrogen_uptake_g_h = 0, .fermentative_hydrogen_emission_g_h = 0, .oxygen_uptake_g_o = 2, .nitrous_oxide_reduction_g_n = 0, .chemodenitrification_dinitrogen_production_g_n = 0, .biological_nitrite_reduction_g_n = 0, .chemodenitrification_nitrous_oxide_production_g_n = 0 },
        .{ .autotrophic_carbon_dioxide_uptake_g_c = 0, .aerobic_heterotrophic_carbon_dioxide_emission_g_c = 0, .denitrification_carbon_dioxide_emission_g_c = 0, .methane_oxidation_g_c = 0, .methanotrophic_carbon_uptake_g_c = 0, .methane_emission_g_c = 0, .hydrogenotrophic_hydrogen_uptake_g_h = 0, .fermentative_hydrogen_emission_g_h = 0, .oxygen_uptake_g_o = 0, .nitrous_oxide_reduction_g_n = 0.25, .chemodenitrification_dinitrogen_production_g_n = 0.5, .biological_nitrite_reduction_g_n = 1, .chemodenitrification_nitrous_oxide_production_g_n = 0.5 },
    };
    var context: ApplyContext = .{ .layer_inputs = &inputs, .result = &state };
    try applyTile(&context, .{ .first = 0, .end = 1 });
    try std.testing.expectEqual(@as(f64, 2), state.source_signed_net_carbon_dioxide_uptake_g_c[0]);
    try std.testing.expectEqual(@as(f64, 0), state.source_signed_nitrous_oxide_flux_g_n[1]);
    try applyTile(&context, .{ .first = 1, .end = 2 });
    try std.testing.expectEqual(@as(f64, -0.75), state.source_signed_dinitrogen_flux_g_n[1]);
    try std.testing.expectEqual(@as(f64, -1.25), state.source_signed_nitrous_oxide_flux_g_n[1]);
}

test "runtime process binding preserves redox scaling and methane carbon roles" {
    const allocator = std.testing.allocator;
    var result = try State.init(allocator, 1);
    defer result.deinit();
    var autotroph = try autotrophic_carbon.State.init(allocator, 1, 1);
    defer autotroph.deinit();
    var products = try respiration_products.State.init(allocator, 1, 1);
    defer products.deinit();
    var methane_state = try methane.State.init(allocator, 1);
    defer methane_state.deinit();
    var oxygen_state = try oxygen.State.init(allocator, 1, 1, 1);
    defer oxygen_state.deinit();
    var fluxes = try nitrogen_fluxes.State.init(allocator, 1, 1);
    defer fluxes.deinit();

    autotroph.carbon_dioxide_uptake_g_c[0] = 10;
    products.carbon_dioxide_g_c[0] = 2;
    products.methane_g_c[0] = 0.25;
    products.hydrogen_g_h[0] = 0.1;
    methane_state.methane_oxidation_respiration_g_c[0] = 0.5;
    methane_state.methane_oxidation_to_biomass_g_c[0] = 0.75;
    methane_state.hydrogenotrophic_methane_g_c[0] = 0.2;
    methane_state.hydrogen_consumption_g_h[0] = 0.4;
    oxygen_state.oxygen_uptake_g_o[0] = 3;
    fluxes.denitrification_respiration_g_c[0] = 2;
    fluxes.nitrous_oxide_reduction_potential_g_n[0] = 0.8;
    fluxes.non_band_heterotrophic_nitrite_reduction_potential_g_n[0] = 0.4;
    fluxes.band_autotrophic_nitrite_reduction_potential_g_n[0] = 0.2;
    fluxes.chemodenitrification_dinitrogen_production_g_n[0] = 0.1;
    fluxes.chemodenitrification_nitrous_oxide_production_g_n[0] = 0.05;
    const redox = [_]f64{0.5};
    var context: ProcessContext = .{ .result = &result, .autotrophic_carbon = &autotroph, .respiration_products = &products, .methane = &methane_state, .oxygen = &oxygen_state, .nitrogen_fluxes = &fluxes, .redox_satisfaction_fraction = &redox };
    try aggregateProcessTile(&context, .{ .first = 0, .end = 1 });

    try std.testing.expectApproxEqAbs(@as(f64, 6.5), result.source_signed_net_carbon_dioxide_uptake_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), result.source_signed_net_methane_uptake_g_c[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), result.source_signed_net_hydrogen_uptake_g_h[0], 1e-15);
    try std.testing.expectEqual(@as(f64, 3), result.oxygen_uptake_g_o[0]);
    try std.testing.expectApproxEqAbs(@as(f64, -0.5), result.source_signed_dinitrogen_flux_g_n[0], 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.05), result.source_signed_nitrous_oxide_flux_g_n[0], 1e-15);
}
