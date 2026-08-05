const std = @import("std");
const chemistry = @import("solute_chemistry_state.zig");

/// Preserves extensive aqueous and water-normalized mineral amounts after an
/// accepted pure-water change. Soil-mass-normalized exchange/site pools are
/// deliberately unchanged.
pub fn rebaseLayer(
    state: *chemistry.State,
    layer: usize,
    old_water_m3: f64,
    new_water_m3: f64,
) !void {
    if (layer >= state.cell_count) return error.SoilChemistryLayerOutOfBounds;
    if (!std.math.isFinite(old_water_m3) or !std.math.isFinite(new_water_m3) or
        old_water_m3 < 0 or new_water_m3 < 0)
        return error.InvalidSoilChemistryWaterCarrier;
    if (old_water_m3 > 0 and new_water_m3 == 0)
        return error.SoilChemistryMassWithoutWaterCarrier;
    const scale = if (old_water_m3 == 0) 0 else old_water_m3 / new_water_m3;
    try validateScalable(state.aqueous[layer], old_water_m3, new_water_m3);
    try validateScalable(state.non_band_phosphate[layer], old_water_m3, new_water_m3);
    try validateScalable(state.band_phosphate[layer], old_water_m3, new_water_m3);
    scaleScalable(&state.aqueous[layer], scale);
    scaleScalable(&state.non_band_phosphate[layer], scale);
    scaleScalable(&state.band_phosphate[layer], scale);
}

fn validateScalable(value: anytype, old_water_m3: f64, new_water_m3: f64) !void {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field| if (comptime scalable(field.name)) {
        const concentration = @field(value, field.name);
        if (!std.math.isFinite(concentration) or concentration < 0)
            return error.InvalidSoilChemistryWaterNormalizedPool;
        if (concentration > 0 and (old_water_m3 == 0 or new_water_m3 == 0))
            return error.SoilChemistryMassWithoutWaterCarrier;
    };
}

fn scaleScalable(value: anytype, scale: f64) void {
    inline for (@typeInfo(@TypeOf(value.*)).@"struct".fields) |field| {
        if (comptime scalable(field.name)) @field(value, field.name) *= scale;
    }
}

fn scalable(comptime name: []const u8) bool {
    return !std.mem.eql(u8, name, "water_mol_per_m3") and
        !std.mem.endsWith(u8, name, "_per_megagram");
}

test "pure-water evaporation preserves aqueous and phosphate extensive amounts" {
    var state = try chemistry.State.init(std.testing.allocator, 1);
    defer state.deinit();
    state.aqueous[0].calcium = 2;
    state.non_band_phosphate[0].dissolved_h2po4_mol_p_per_m3 = 3;
    state.non_band_phosphate[0].adsorbed_h2po4_mol_p_per_megagram = 5;
    try rebaseLayer(&state, 0, 4, 2);
    try std.testing.expectEqual(@as(f64, 4), state.aqueous[0].calcium);
    try std.testing.expectEqual(@as(f64, 6), state.non_band_phosphate[0].dissolved_h2po4_mol_p_per_m3);
    try std.testing.expectEqual(@as(f64, 5), state.non_band_phosphate[0].adsorbed_h2po4_mol_p_per_megagram);
}
