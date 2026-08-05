const std = @import("std");

pub const State = struct {
    net_primary_productivity_g_c: f64, // TNPP, cumulative
    net_carbon_exchange_g_c_h: f64, // TCNET
    ecosystem_respiration_g_c_h: f64, // RECO
    cumulative_canopy_exchange_g_c: f64, // TCAN
    net_biome_productivity_g_c: f64, // TNBP, cumulative
};

pub const Inputs = struct {
    gross_primary_productivity_g_c: f64, // TGPP
    signed_autotrophic_respiration_g_c: f64, // TRAU
    canopy_net_fixation_g_c_h: f64, // TCCAN
    ground_co2_exchange_g_c_h: f64, // HCO2G
    signed_heterotrophic_respiration_g_c: f64, // THRE
    dissolved_organic_runoff_g_c: f64, // UDOCQ
    dissolved_inorganic_runoff_g_c: f64, // UDICQ
    dissolved_organic_drainage_g_c: f64, // UDOCD
    dissolved_inorganic_drainage_g_c: f64, // UDICD
    harvested_carbon_g_c: f64, // XHVSTC
};

fn finiteStruct(value: anytype) bool {
    inline for (std.meta.fields(@TypeOf(value))) |field|
        if (!std.math.isFinite(@field(value, field.name))) return false;
    return true;
}

/// Direct translation of REDIST 10652--10668 under `NFZ == NFH`.
pub fn closeout(current_subcycle: usize, final_subcycle: usize, state: *State, inputs: Inputs) !void {
    if (!finiteStruct(state.*) or !finiteStruct(inputs)) return error.InvalidEcosystemCarbonCloseoutInput;
    if (current_subcycle != final_subcycle) return;

    var next = state.*;
    next.net_primary_productivity_g_c = inputs.gross_primary_productivity_g_c + inputs.signed_autotrophic_respiration_g_c;
    next.net_carbon_exchange_g_c_h = inputs.canopy_net_fixation_g_c_h + inputs.ground_co2_exchange_g_c_h;
    next.ecosystem_respiration_g_c_h = state.ecosystem_respiration_g_c_h + inputs.ground_co2_exchange_g_c_h;
    next.cumulative_canopy_exchange_g_c = state.cumulative_canopy_exchange_g_c + inputs.canopy_net_fixation_g_c_h;
    next.net_biome_productivity_g_c = next.net_primary_productivity_g_c + inputs.signed_heterotrophic_respiration_g_c -
        inputs.dissolved_organic_runoff_g_c - inputs.dissolved_inorganic_runoff_g_c -
        inputs.dissolved_organic_drainage_g_c - inputs.dissolved_inorganic_drainage_g_c -
        inputs.harvested_carbon_g_c;
    if (!finiteStruct(next)) return error.NonFiniteEcosystemCarbonCloseoutResult;
    state.* = next;
}

fn fixtureInputs() Inputs {
    return .{ .gross_primary_productivity_g_c = 100, .signed_autotrophic_respiration_g_c = -20, .canopy_net_fixation_g_c_h = 5, .ground_co2_exchange_g_c_h = -2, .signed_heterotrophic_respiration_g_c = -10, .dissolved_organic_runoff_g_c = 1, .dissolved_inorganic_runoff_g_c = 2, .dissolved_organic_drainage_g_c = 3, .dissolved_inorganic_drainage_g_c = 4, .harvested_carbon_g_c = 5 };
}

test "REDIST ecosystem carbon closeout preserves exact order and signs" {
    var state = State{ .net_primary_productivity_g_c = 0, .net_carbon_exchange_g_c_h = 0, .ecosystem_respiration_g_c_h = 7, .cumulative_canopy_exchange_g_c = 11, .net_biome_productivity_g_c = 0 };
    try closeout(4, 4, &state, fixtureInputs());
    try std.testing.expectEqual(@as(f64, 80), state.net_primary_productivity_g_c);
    try std.testing.expectEqual(@as(f64, 3), state.net_carbon_exchange_g_c_h);
    try std.testing.expectEqual(@as(f64, 5), state.ecosystem_respiration_g_c_h);
    try std.testing.expectEqual(@as(f64, 16), state.cumulative_canopy_exchange_g_c);
    try std.testing.expectEqual(@as(f64, 55), state.net_biome_productivity_g_c);
}

test "REDIST ecosystem carbon closeout gate and validation are atomic" {
    var state = State{ .net_primary_productivity_g_c = 1, .net_carbon_exchange_g_c_h = 2, .ecosystem_respiration_g_c_h = 3, .cumulative_canopy_exchange_g_c = 4, .net_biome_productivity_g_c = 5 };
    const before = state;
    try closeout(3, 4, &state, fixtureInputs());
    try std.testing.expectEqualDeep(before, state);
    var invalid = fixtureInputs();
    invalid.harvested_carbon_g_c = std.math.nan(f64);
    try std.testing.expectError(error.InvalidEcosystemCarbonCloseoutInput, closeout(4, 4, &state, invalid));
    try std.testing.expectEqualDeep(before, state);
}
