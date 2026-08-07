const std = @import("std");

pub const Inputs = struct {
    autotrophic_carbon_uptake_g_c: f64,
    heterotrophic_carbon_dioxide_emission_g_c: f64,
    denitrification_carbon_dioxide_emission_g_c: f64,
    methane_oxidation_g_c: f64,
    methanotrophic_methane_uptake_g_c: f64,
    methane_emission_g_c: f64,
    hydrogen_uptake_g_h: f64,
    hydrogen_emission_g_h: f64,
    oxygen_uptake_g_o: f64,
    nitrous_oxide_reduction_g_n: f64,
    chemodenitrification_dinitrogen_production_g_n: f64,
    non_band_nitrite_reduction_g_n: f64,
    band_nitrite_reduction_g_n: f64,
    non_band_chemodenitrification_nitrous_oxide_g_n: f64,
    band_chemodenitrification_nitrous_oxide_g_n: f64,
};

pub const State = struct {
    net_carbon_dioxide_uptake_g_c: f64 = 0,
    net_methane_uptake_g_c: f64 = 0,
    net_hydrogen_uptake_g_h: f64 = 0,
    oxygen_uptake_g_o: f64 = 0,
    net_dinitrogen_uptake_g_n: f64 = 0,
    net_nitrous_oxide_uptake_g_n: f64 = 0,
};

/// Exact NITRO.F 4007--4022 layer gas-flux assembly.
/// Uptake is positive; production/emission therefore appears as negative flux.
pub fn calculate(state: *State, inputs: Inputs) !void {
    try validate(inputs);
    const result: State = .{
        .net_carbon_dioxide_uptake_g_c = inputs.autotrophic_carbon_uptake_g_c -
            inputs.heterotrophic_carbon_dioxide_emission_g_c -
            inputs.denitrification_carbon_dioxide_emission_g_c -
            inputs.methane_oxidation_g_c,
        .net_methane_uptake_g_c = inputs.methane_oxidation_g_c +
            inputs.methanotrophic_methane_uptake_g_c -
            inputs.methane_emission_g_c,
        .net_hydrogen_uptake_g_h = inputs.hydrogen_uptake_g_h - inputs.hydrogen_emission_g_h,
        .oxygen_uptake_g_o = inputs.oxygen_uptake_g_o,
        .net_dinitrogen_uptake_g_n = -inputs.nitrous_oxide_reduction_g_n -
            inputs.chemodenitrification_dinitrogen_production_g_n,
        .net_nitrous_oxide_uptake_g_n = -inputs.non_band_nitrite_reduction_g_n -
            inputs.band_nitrite_reduction_g_n -
            inputs.non_band_chemodenitrification_nitrous_oxide_g_n -
            inputs.band_chemodenitrification_nitrous_oxide_g_n +
            inputs.nitrous_oxide_reduction_g_n,
    };
    inline for (@typeInfo(State).@"struct".fields) |field|
        if (!std.math.isFinite(@field(result, field.name)))
            return error.InvalidNetGasFluxResult;
    state.* = result;
}

fn validate(inputs: Inputs) !void {
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        const value = @field(inputs, field.name);
        if (!std.math.isFinite(value) or value < 0)
            return error.InvalidNetGasFluxInput;
    }
}

fn fixture() Inputs {
    return .{
        .autotrophic_carbon_uptake_g_c = 10,
        .heterotrophic_carbon_dioxide_emission_g_c = 4,
        .denitrification_carbon_dioxide_emission_g_c = 1,
        .methane_oxidation_g_c = 2,
        .methanotrophic_methane_uptake_g_c = 3,
        .methane_emission_g_c = 1,
        .hydrogen_uptake_g_h = 4,
        .hydrogen_emission_g_h = 1,
        .oxygen_uptake_g_o = 5,
        .nitrous_oxide_reduction_g_n = 2,
        .chemodenitrification_dinitrogen_production_g_n = 0,
        .non_band_nitrite_reduction_g_n = 3,
        .band_nitrite_reduction_g_n = 1,
        .non_band_chemodenitrification_nitrous_oxide_g_n = 0.5,
        .band_chemodenitrification_nitrous_oxide_g_n = 0.25,
    };
}

test "net gas assembly preserves uptake-positive sign convention" {
    var state: State = .{};
    try calculate(&state, fixture());
    try std.testing.expectEqual(3, state.net_carbon_dioxide_uptake_g_c);
    try std.testing.expectEqual(4, state.net_methane_uptake_g_c);
    try std.testing.expectEqual(3, state.net_hydrogen_uptake_g_h);
    try std.testing.expectEqual(5, state.oxygen_uptake_g_o);
}

test "dinitrogen production is represented by negative uptake" {
    var state: State = .{};
    try calculate(&state, fixture());
    try std.testing.expectEqual(-2, state.net_dinitrogen_uptake_g_n);
}

test "N2O reduction offsets N2O production from nitrite reduction" {
    var state: State = .{};
    try calculate(&state, fixture());
    try std.testing.expectApproxEqAbs(-2.75, state.net_nitrous_oxide_uptake_g_n, 1e-12);
}

test "invalid input leaves state unchanged" {
    var state: State = .{ .net_carbon_dioxide_uptake_g_c = 7 };
    var inputs = fixture();
    inputs.oxygen_uptake_g_o = std.math.nan(f64);
    try std.testing.expectError(error.InvalidNetGasFluxInput, calculate(&state, inputs));
    try std.testing.expectEqual(7, state.net_carbon_dioxide_uptake_g_c);
}
