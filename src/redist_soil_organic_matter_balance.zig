const std = @import("std");

pub const ElementMass = struct { carbon_g: f64, nitrogen_g: f64, phosphorus_g: f64 };

pub const ClassSolubles = struct {
    doc_micropore_g_c: f64,
    doc_macropore_g_c: f64,
    adsorbed_carbon_g_c: f64,
    acetate_micropore_g_c: f64,
    acetate_macropore_g_c: f64,
    adsorbed_acetate_g_c: f64,
    don_micropore_g_n: f64,
    don_macropore_g_n: f64,
    adsorbed_nitrogen_g_n: f64,
    dop_micropore_g_p: f64,
    dop_macropore_g_p: f64,
    adsorbed_phosphorus_g_p: f64,
};

pub const LayerState = struct {
    microbial_carbon_g: f64, // OMCL
    microbial_nitrogen_g: f64, // OMNL
    previous_charcoal_carbon_g: f64, // ORGCCX
    organic_carbon_g: f64, // ORGC
    organic_nitrogen_g: f64, // ORGN
    charcoal_carbon_g: f64, // ORGCC
    charcoal_nitrogen_g: f64, // ORGNC
    residue_carbon_g: f64, // ORGR
    organic_carbon_change_g: f64, // DORGC
};

pub const Balance = struct {
    microbial_carbon_g: f64,
    microbial_nitrogen_g: f64,
    microbial_phosphorus_g: f64,
    landscape_residue_carbon_g: f64,
    grid_residue_carbon_g: f64,
    landscape_residue_nitrogen_g: f64,
    grid_residue_nitrogen_g: f64,
    landscape_residue_phosphorus_g: f64,
    grid_residue_phosphorus_g: f64,
    landscape_humus_carbon_g: f64,
    grid_humus_carbon_g: f64,
    landscape_humus_nitrogen_g: f64,
    grid_humus_nitrogen_g: f64,
    landscape_humus_phosphorus_g: f64,
    grid_humus_phosphorus_g: f64,
    landscape_sediment_megagrams: f64,
};

fn finiteStruct(value: anytype) bool {
    inline for (@typeInfo(@TypeOf(value)).@"struct".fields) |field|
        if (!std.math.isFinite(@field(value, field.name))) return false;
    return true;
}

/// Direct translation of REDIST 6778--6910. Runtime arrays are layer-major,
/// followed by legacy K, N, then M order (M varies fastest).
pub fn aggregateLayers(
    states: []LayerState,
    balance: *Balance,
    microbial_by_layer_k_n_m: []const ElementMass,
    residue_by_layer_k_m: []const ElementMass,
    solubles_by_layer_k: []const ClassSolubles,
    soc_by_layer_k_m: []const ElementMass,
    previous_organic_carbon_g: []const f64,
    erosion_mode: u8,
    first_soil_layer_index: usize,
    top_layer_erosion_carbon_g: f64,
) !void {
    const microbial_per_layer = 6 * 7 * 3;
    const residue_per_layer = 5 * 2;
    const soc_per_layer = 5 * 5;
    if (states.len == 0 or microbial_by_layer_k_n_m.len != states.len * microbial_per_layer or
        residue_by_layer_k_m.len != states.len * residue_per_layer or
        solubles_by_layer_k.len != states.len * 5 or soc_by_layer_k_m.len != states.len * soc_per_layer or
        previous_organic_carbon_g.len != states.len or first_soil_layer_index >= states.len)
        return error.SoilOrganicMatterDimensionMismatch;
    if (!finiteStruct(balance.*) or !std.math.isFinite(top_layer_erosion_carbon_g))
        return error.InvalidSoilOrganicMatter;
    for (states) |state| if (!finiteStruct(state)) return error.InvalidSoilOrganicMatter;
    for (microbial_by_layer_k_n_m) |value| if (!finiteStruct(value)) return error.InvalidSoilOrganicMatter;
    for (residue_by_layer_k_m) |value| if (!finiteStruct(value)) return error.InvalidSoilOrganicMatter;
    for (solubles_by_layer_k) |value| if (!finiteStruct(value)) return error.InvalidSoilOrganicMatter;
    for (soc_by_layer_k_m) |value| if (!finiteStruct(value)) return error.InvalidSoilOrganicMatter;
    for (previous_organic_carbon_g) |value| if (!std.math.isFinite(value)) return error.InvalidSoilOrganicMatter;

    var next_balance = balance.*;
    for (states, 0..) |*state, layer| {
        var residue = std.mem.zeroes(ElementMass);
        var humus = std.mem.zeroes(ElementMass);
        var residue_charcoal = std.mem.zeroes(ElementMass);
        var humus_charcoal = std.mem.zeroes(ElementMass);
        var microbial_c: f64 = 0.0;
        var microbial_n: f64 = 0.0;
        const microbial_start = layer * microbial_per_layer;
        for (0..6) |class| {
            for (0..7) |population| {
                for (0..3) |component| {
                    const value = microbial_by_layer_k_n_m[microbial_start + class * 21 + population * 3 + component];
                    const destination = if (class != 4) &residue else &humus;
                    destination.carbon_g = destination.carbon_g + value.carbon_g;
                    destination.nitrogen_g = destination.nitrogen_g + value.nitrogen_g;
                    destination.phosphorus_g = destination.phosphorus_g + value.phosphorus_g;
                    next_balance.microbial_carbon_g = next_balance.microbial_carbon_g + value.carbon_g;
                    next_balance.microbial_nitrogen_g = next_balance.microbial_nitrogen_g + value.nitrogen_g;
                    next_balance.microbial_phosphorus_g = next_balance.microbial_phosphorus_g + value.phosphorus_g;
                    microbial_c = microbial_c + value.carbon_g;
                    microbial_n = microbial_n + value.nitrogen_g;
                }
            }
        }

        const residue_start = layer * residue_per_layer;
        const soluble_start = layer * 5;
        const soc_start = layer * soc_per_layer;
        for (0..5) |class| {
            const destination = if (class != 4) &residue else &humus;
            for (0..2) |component| {
                const value = residue_by_layer_k_m[residue_start + class * 2 + component];
                destination.carbon_g = destination.carbon_g + value.carbon_g;
                destination.nitrogen_g = destination.nitrogen_g + value.nitrogen_g;
                destination.phosphorus_g = destination.phosphorus_g + value.phosphorus_g;
            }
            const soluble = solubles_by_layer_k[soluble_start + class];
            destination.carbon_g = destination.carbon_g + soluble.doc_micropore_g_c + soluble.doc_macropore_g_c + soluble.adsorbed_carbon_g_c +
                soluble.acetate_micropore_g_c + soluble.acetate_macropore_g_c + soluble.adsorbed_acetate_g_c;
            destination.nitrogen_g = destination.nitrogen_g + soluble.don_micropore_g_n + soluble.don_macropore_g_n + soluble.adsorbed_nitrogen_g_n;
            destination.phosphorus_g = destination.phosphorus_g + soluble.dop_micropore_g_p + soluble.dop_macropore_g_p + soluble.adsorbed_phosphorus_g_p;
            for (0..5) |fraction| {
                const value = soc_by_layer_k_m[soc_start + class * 5 + fraction];
                if (fraction <= 3) {
                    destination.carbon_g = destination.carbon_g + value.carbon_g;
                    destination.nitrogen_g = destination.nitrogen_g + value.nitrogen_g;
                    destination.phosphorus_g = destination.phosphorus_g + value.phosphorus_g;
                } else {
                    const charcoal = if (class != 4) &residue_charcoal else &humus_charcoal;
                    charcoal.carbon_g = charcoal.carbon_g + value.carbon_g;
                    charcoal.nitrogen_g = charcoal.nitrogen_g + value.nitrogen_g;
                    charcoal.phosphorus_g = charcoal.phosphorus_g + value.phosphorus_g;
                }
            }
        }

        const organic_carbon = residue.carbon_g + humus.carbon_g;
        const organic_nitrogen = residue.nitrogen_g + humus.nitrogen_g;
        const charcoal_carbon = residue_charcoal.carbon_g + humus_charcoal.carbon_g;
        const charcoal_nitrogen = residue_charcoal.nitrogen_g + humus_charcoal.nitrogen_g;
        var change: f64 = 0.0;
        if (erosion_mode == 2 or erosion_mode == 3) {
            change = previous_organic_carbon_g[layer] - organic_carbon - charcoal_carbon;
            if (layer == first_soil_layer_index) change = change + top_layer_erosion_carbon_g;
        }
        state.* = .{
            .microbial_carbon_g = microbial_c,
            .microbial_nitrogen_g = microbial_n,
            .previous_charcoal_carbon_g = state.charcoal_carbon_g,
            .organic_carbon_g = organic_carbon,
            .organic_nitrogen_g = organic_nitrogen,
            .charcoal_carbon_g = charcoal_carbon,
            .charcoal_nitrogen_g = charcoal_nitrogen,
            .residue_carbon_g = residue.carbon_g,
            .organic_carbon_change_g = change,
        };
        next_balance.landscape_residue_carbon_g = next_balance.landscape_residue_carbon_g + residue.carbon_g + residue_charcoal.carbon_g;
        next_balance.grid_residue_carbon_g = next_balance.grid_residue_carbon_g + residue.carbon_g + residue_charcoal.carbon_g;
        next_balance.landscape_residue_nitrogen_g = next_balance.landscape_residue_nitrogen_g + residue.nitrogen_g + residue_charcoal.nitrogen_g;
        next_balance.grid_residue_nitrogen_g = next_balance.grid_residue_nitrogen_g + residue.nitrogen_g + residue_charcoal.nitrogen_g;
        next_balance.landscape_residue_phosphorus_g = next_balance.landscape_residue_phosphorus_g + residue.phosphorus_g + residue_charcoal.phosphorus_g;
        next_balance.grid_residue_phosphorus_g = next_balance.grid_residue_phosphorus_g + residue.phosphorus_g + residue_charcoal.phosphorus_g;
        next_balance.landscape_humus_carbon_g = next_balance.landscape_humus_carbon_g + humus.carbon_g + humus_charcoal.carbon_g;
        next_balance.grid_humus_carbon_g = next_balance.grid_humus_carbon_g + humus.carbon_g + humus_charcoal.carbon_g;
        next_balance.landscape_humus_nitrogen_g = next_balance.landscape_humus_nitrogen_g + humus.nitrogen_g + humus_charcoal.nitrogen_g;
        next_balance.grid_humus_nitrogen_g = next_balance.grid_humus_nitrogen_g + humus.nitrogen_g + humus_charcoal.nitrogen_g;
        next_balance.landscape_humus_phosphorus_g = next_balance.landscape_humus_phosphorus_g + humus.phosphorus_g + humus_charcoal.phosphorus_g;
        next_balance.grid_humus_phosphorus_g = next_balance.grid_humus_phosphorus_g + humus.phosphorus_g + humus_charcoal.phosphorus_g;
        next_balance.landscape_sediment_megagrams = next_balance.landscape_sediment_megagrams + (residue.carbon_g + humus.carbon_g) * 1.0e-6;
        if (!finiteStruct(state.*) or !finiteStruct(next_balance)) return error.NonFiniteSoilOrganicMatter;
    }
    balance.* = next_balance;
}

test "REDIST organic matter routes class four and fraction five to charcoal owners" {
    var states = [_]LayerState{std.mem.zeroes(LayerState)};
    states[0].charcoal_carbon_g = 9.0;
    var balance = std.mem.zeroes(Balance);
    var microbial = [_]ElementMass{std.mem.zeroes(ElementMass)} ** 126;
    var residues = [_]ElementMass{std.mem.zeroes(ElementMass)} ** 10;
    var solubles = [_]ClassSolubles{std.mem.zeroes(ClassSolubles)} ** 5;
    var soc = [_]ElementMass{std.mem.zeroes(ElementMass)} ** 25;
    microbial[0] = .{ .carbon_g = 1, .nitrogen_g = 2, .phosphorus_g = 3 };
    microbial[4 * 21] = .{ .carbon_g = 4, .nitrogen_g = 5, .phosphorus_g = 6 };
    soc[4] = .{ .carbon_g = 7, .nitrogen_g = 8, .phosphorus_g = 9 };
    soc[4 * 5 + 4] = .{ .carbon_g = 10, .nitrogen_g = 11, .phosphorus_g = 12 };
    const previous = [_]f64{100};
    try aggregateLayers(&states, &balance, &microbial, &residues, &solubles, &soc, &previous, 2, 0, 3);
    try std.testing.expectEqual(@as(f64, 5), states[0].organic_carbon_g);
    try std.testing.expectEqual(@as(f64, 17), states[0].charcoal_carbon_g);
    try std.testing.expectEqual(@as(f64, 81), states[0].organic_carbon_change_g);
    try std.testing.expectEqual(@as(f64, 9), states[0].previous_charcoal_carbon_g);
    try std.testing.expectEqual(@as(f64, 5), balance.microbial_carbon_g);
}

test "REDIST organic matter preserves runtime layer and top erosion gate" {
    var states = [_]LayerState{ std.mem.zeroes(LayerState), std.mem.zeroes(LayerState) };
    var balance = std.mem.zeroes(Balance);
    const microbial = [_]ElementMass{std.mem.zeroes(ElementMass)} ** 252;
    const residues = [_]ElementMass{std.mem.zeroes(ElementMass)} ** 20;
    const solubles = [_]ClassSolubles{std.mem.zeroes(ClassSolubles)} ** 10;
    const soc = [_]ElementMass{std.mem.zeroes(ElementMass)} ** 50;
    const previous = [_]f64{ 10, 20 };
    try aggregateLayers(&states, &balance, &microbial, &residues, &solubles, &soc, &previous, 3, 1, 4);
    try std.testing.expectEqual(@as(f64, 10), states[0].organic_carbon_change_g);
    try std.testing.expectEqual(@as(f64, 24), states[1].organic_carbon_change_g);
}

test "REDIST organic matter non-erosion mode zeros change" {
    var states = [_]LayerState{std.mem.zeroes(LayerState)};
    var balance = std.mem.zeroes(Balance);
    const microbial = [_]ElementMass{std.mem.zeroes(ElementMass)} ** 126;
    const residues = [_]ElementMass{std.mem.zeroes(ElementMass)} ** 10;
    const solubles = [_]ClassSolubles{std.mem.zeroes(ClassSolubles)} ** 5;
    const soc = [_]ElementMass{std.mem.zeroes(ElementMass)} ** 25;
    const previous = [_]f64{10};
    try aggregateLayers(&states, &balance, &microbial, &residues, &solubles, &soc, &previous, 1, 0, 4);
    try std.testing.expectEqual(@as(f64, 0), states[0].organic_carbon_change_g);
}

test "REDIST organic matter rejects dimensions invalid input and overflow" {
    var states = [_]LayerState{std.mem.zeroes(LayerState)};
    var balance = std.mem.zeroes(Balance);
    var microbial = [_]ElementMass{std.mem.zeroes(ElementMass)} ** 126;
    const residues = [_]ElementMass{std.mem.zeroes(ElementMass)} ** 10;
    const solubles = [_]ClassSolubles{std.mem.zeroes(ClassSolubles)} ** 5;
    const soc = [_]ElementMass{std.mem.zeroes(ElementMass)} ** 25;
    const no_previous: [0]f64 = .{};
    try std.testing.expectError(error.SoilOrganicMatterDimensionMismatch, aggregateLayers(&states, &balance, &microbial, &residues, &solubles, &soc, &no_previous, 0, 0, 0));
    const previous = [_]f64{0};
    microbial[0].carbon_g = std.math.nan(f64);
    try std.testing.expectError(error.InvalidSoilOrganicMatter, aggregateLayers(&states, &balance, &microbial, &residues, &solubles, &soc, &previous, 0, 0, 0));
    microbial[0].carbon_g = std.math.floatMax(f64);
    microbial[1].carbon_g = std.math.floatMax(f64);
    try std.testing.expectError(error.NonFiniteSoilOrganicMatter, aggregateLayers(&states, &balance, &microbial, &residues, &solubles, &soc, &previous, 0, 0, 0));
}
