const std = @import("std");
const Organic = @import("soil_organic_initialization.zig");

pub fn captureHourStart(state: *const Organic.State, carbon_g_c_by_layer: []f64) !void {
    if (carbon_g_c_by_layer.len != state.layer_count) return error.SoilOrganicCarbonLedgerDimensionMismatch;
    for (0..state.layer_count) |layer| carbon_g_c_by_layer[layer] = try state.totalCarbon_g_c(layer);
}

/// Publishes REDIST `DORGC` only after validating every layer, so a late
/// non-finite pool cannot leave a partially advanced geometry driver.
pub fn publishAcceptedHourlyChange(state: *const Organic.State, hour_start_carbon_g_c: []const f64, change_g_c_per_h: []f64) !void {
    if (hour_start_carbon_g_c.len != state.layer_count or change_g_c_per_h.len != state.layer_count) return error.SoilOrganicCarbonLedgerDimensionMismatch;
    for (0..state.layer_count) |layer| {
        const final_carbon_g_c = try state.totalCarbon_g_c(layer);
        // Source sign: DORGC = ORGCX - ORGC - ORGCC. Positive is SOC loss,
        // which moves the upper boundary downward and contracts the layer.
        const change_g_c = hour_start_carbon_g_c[layer] - final_carbon_g_c;
        if (!std.math.isFinite(hour_start_carbon_g_c[layer]) or !std.math.isFinite(change_g_c)) return error.NonFiniteHourlySoilOrganicCarbonChange;
    }
    for (0..state.layer_count) |layer| change_g_c_per_h[layer] = hour_start_carbon_g_c[layer] - (try state.totalCarbon_g_c(layer));
}

test "REDIST DORGC includes every runtime organic carbon owner and publishes atomically" {
    var state = try Organic.State.init(std.testing.allocator, 2);
    defer state.deinit();
    var beginning = [_]f64{0} ** 2;
    var change = [_]f64{ 8, 9 };
    try captureHourStart(&state, &beginning);
    state.microbial[0].carbon_g_c = 1;
    state.residue[0].carbon_g_c = 2;
    state.dissolved[0].carbon_g_c = 3;
    state.adsorbed[0].carbon_g_c = 4;
    state.dissolved_acetate_carbon_g_c[0] = 5;
    state.adsorbed_acetate_carbon_g_c[0] = 6;
    state.structural[0].carbon_g_c = 7;
    try publishAcceptedHourlyChange(&state, &beginning, &change);
    try std.testing.expectEqual(@as(f64, -28), change[0]);
    try std.testing.expectEqual(@as(f64, 0), change[1]);

    change = .{ 8, 9 };
    const second_layer_structural = Organic.substrate_count * Organic.structural_fraction_count;
    state.structural[second_layer_structural].carbon_g_c = std.math.nan(f64);
    try std.testing.expectError(error.InvalidOrganicCarbonState, publishAcceptedHourlyChange(&state, &beginning, &change));
    try std.testing.expectEqualSlices(f64, &.{ 8, 9 }, &change);
}
