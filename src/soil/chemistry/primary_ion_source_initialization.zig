const std = @import("std");

pub const Dimensions = struct {
    column_count: usize,
    row_count: usize,
    layer_count_including_surface: usize,
};

pub const LoopPosition = struct {
    source_index: usize, // one-based K
    day_index: usize, // one-based I
    layer_index: usize, // L; zero is the surface layer
};

/// Cell-major, then layer storage for raw soil-file concentrations.
pub const Inputs = struct {
    aluminum_mol_per_megagram_by_cell_layer: []const f64, // CAL
    iron_mol_per_megagram_by_cell_layer: []const f64, // CFE
    calcium_mol_per_megagram_by_cell_layer: []const f64, // CCA
    magnesium_mol_per_megagram_by_cell_layer: []const f64, // CMG
    sodium_mol_per_megagram_by_cell_layer: []const f64, // CNA
    potassium_mol_per_megagram_by_cell_layer: []const f64, // CKA
    sulfate_mol_s_per_megagram_by_cell_layer: []const f64, // CSO4
    chloride_mol_per_megagram_by_cell_layer: []const f64, // CCL
    aluminum_phosphate_mol_per_megagram_by_cell_layer: []const f64, // CALPO
    iron_phosphate_mol_per_megagram_by_cell_layer: []const f64, // CFEPO
    calcium_hydrogen_phosphate_mol_per_megagram_by_cell_layer: []const f64, // CCAPD
    apatite_mol_per_megagram_by_cell_layer: []const f64, // CCAPH
    aluminum_hydroxide_mol_per_megagram_by_cell_layer: []const f64, // CALOH
    iron_hydroxide_mol_per_megagram_by_cell_layer: []const f64, // CFEOH
    calcium_carbonate_mol_per_megagram_by_cell_layer: []const f64, // CCACO
    calcium_sulfate_mol_per_megagram_by_cell_layer: []const f64, // CCASO
};

pub const State = struct {
    aluminum_mol_per_megagram: f64 = 0, // CALX
    iron_mol_per_megagram: f64 = 0, // CFEX
    calcium_mol_per_megagram: f64 = 0, // CCAX
    magnesium_mol_per_megagram: f64 = 0, // CMGX
    sodium_mol_per_megagram: f64 = 0, // CNAX
    potassium_mol_per_megagram: f64 = 0, // CKAX
    sulfate_mol_s_per_megagram: f64 = 0, // CSOX
    chloride_mol_per_megagram: f64 = 0, // CCLX
    aluminum_phosphate_mol_per_megagram: f64 = 0, // CALPOX
    iron_phosphate_mol_per_megagram: f64 = 0, // CFEPOX
    calcium_hydrogen_phosphate_mol_per_megagram: f64 = 0, // CCAPDX
    apatite_mol_per_megagram: f64 = 0, // CCAPHX
    aluminum_hydroxide_mol_per_megagram: f64 = 0, // CALOHX
    iron_hydroxide_mol_per_megagram: f64 = 0, // CFEOHX
    calcium_carbonate_mol_per_megagram: f64 = 0, // CCACOX
    calcium_sulfate_mol_per_megagram: f64 = 0, // CCASOX
};

/// Direct translation of `starte.f` lines 192--207 within the `K=3,I=1`
/// soil branch. Negative finite Al/Fe values are retained because subsequent
/// source branches interpret them; this unit rejects only non-finite data.
pub fn initialize(
    position: LoopPosition,
    dimensions: Dimensions,
    row: usize,
    column: usize,
    inputs: Inputs,
    state: *State,
) !bool {
    if (position.source_index != 3 or position.day_index != 1) return false;

    if (dimensions.column_count == 0 or dimensions.row_count == 0 or
        dimensions.layer_count_including_surface == 0 or
        row >= dimensions.row_count or column >= dimensions.column_count or
        position.layer_index >= dimensions.layer_count_including_surface)
        return error.InvalidSoilIonSourceDimensions;
    const cell_count = std.math.mul(usize, dimensions.column_count, dimensions.row_count) catch
        return error.InvalidSoilIonSourceDimensions;
    const value_count = std.math.mul(usize, cell_count, dimensions.layer_count_including_surface) catch
        return error.InvalidSoilIonSourceDimensions;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (@field(inputs, field.name).len != value_count)
            return error.InvalidSoilIonSourceDimensions;
    }

    const cell = row * dimensions.column_count + column;
    const value_index = cell * dimensions.layer_count_including_surface + position.layer_index;
    inline for (@typeInfo(Inputs).@"struct".fields) |field| {
        if (!std.math.isFinite(@field(inputs, field.name)[value_index]))
            return error.NonFiniteSoilIonSourceInput;
    }
    const next: State = .{
        .aluminum_mol_per_megagram = inputs.aluminum_mol_per_megagram_by_cell_layer[value_index],
        .iron_mol_per_megagram = inputs.iron_mol_per_megagram_by_cell_layer[value_index],
        .calcium_mol_per_megagram = inputs.calcium_mol_per_megagram_by_cell_layer[value_index],
        .magnesium_mol_per_megagram = inputs.magnesium_mol_per_megagram_by_cell_layer[value_index],
        .sodium_mol_per_megagram = inputs.sodium_mol_per_megagram_by_cell_layer[value_index],
        .potassium_mol_per_megagram = inputs.potassium_mol_per_megagram_by_cell_layer[value_index],
        .sulfate_mol_s_per_megagram = inputs.sulfate_mol_s_per_megagram_by_cell_layer[value_index],
        .chloride_mol_per_megagram = inputs.chloride_mol_per_megagram_by_cell_layer[value_index],
        .aluminum_phosphate_mol_per_megagram = inputs.aluminum_phosphate_mol_per_megagram_by_cell_layer[value_index],
        .iron_phosphate_mol_per_megagram = inputs.iron_phosphate_mol_per_megagram_by_cell_layer[value_index],
        .calcium_hydrogen_phosphate_mol_per_megagram = inputs.calcium_hydrogen_phosphate_mol_per_megagram_by_cell_layer[value_index],
        .apatite_mol_per_megagram = inputs.apatite_mol_per_megagram_by_cell_layer[value_index],
        .aluminum_hydroxide_mol_per_megagram = inputs.aluminum_hydroxide_mol_per_megagram_by_cell_layer[value_index],
        .iron_hydroxide_mol_per_megagram = inputs.iron_hydroxide_mol_per_megagram_by_cell_layer[value_index],
        .calcium_carbonate_mol_per_megagram = inputs.calcium_carbonate_mol_per_megagram_by_cell_layer[value_index],
        .calcium_sulfate_mol_per_megagram = inputs.calcium_sulfate_mol_per_megagram_by_cell_layer[value_index],
    };
    state.* = next;
    return true;
}

fn uniformInputs(values: []const f64) Inputs {
    return .{
        .aluminum_mol_per_megagram_by_cell_layer = values,
        .iron_mol_per_megagram_by_cell_layer = values,
        .calcium_mol_per_megagram_by_cell_layer = values,
        .magnesium_mol_per_megagram_by_cell_layer = values,
        .sodium_mol_per_megagram_by_cell_layer = values,
        .potassium_mol_per_megagram_by_cell_layer = values,
        .sulfate_mol_s_per_megagram_by_cell_layer = values,
        .chloride_mol_per_megagram_by_cell_layer = values,
        .aluminum_phosphate_mol_per_megagram_by_cell_layer = values,
        .iron_phosphate_mol_per_megagram_by_cell_layer = values,
        .calcium_hydrogen_phosphate_mol_per_megagram_by_cell_layer = values,
        .apatite_mol_per_megagram_by_cell_layer = values,
        .aluminum_hydroxide_mol_per_megagram_by_cell_layer = values,
        .iron_hydroxide_mol_per_megagram_by_cell_layer = values,
        .calcium_carbonate_mol_per_megagram_by_cell_layer = values,
        .calcium_sulfate_mol_per_megagram_by_cell_layer = values,
    };
}

test "STARTE soil source snapshot preserves runtime cell-layer topology and order" {
    const base = [_]f64{ 1, 2, 3, 4 };
    const aluminum = [_]f64{ 10, 20, 30, -40 };
    const sulfate = [_]f64{ 5, 6, 7, 8 };
    const apatite = [_]f64{ 9, 10, 11, 12 };
    const calcium_sulfate = [_]f64{ 13, 14, 15, 16 };
    var inputs = uniformInputs(&base);
    inputs.aluminum_mol_per_megagram_by_cell_layer = &aluminum;
    inputs.sulfate_mol_s_per_megagram_by_cell_layer = &sulfate;
    inputs.apatite_mol_per_megagram_by_cell_layer = &apatite;
    inputs.calcium_sulfate_mol_per_megagram_by_cell_layer = &calcium_sulfate;
    var state: State = .{};
    try std.testing.expect(try initialize(.{ .source_index = 3, .day_index = 1, .layer_index = 1 }, .{ .column_count = 2, .row_count = 1, .layer_count_including_surface = 2 }, 0, 1, inputs, &state));
    try std.testing.expectEqual(@as(f64, -40), state.aluminum_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 8), state.sulfate_mol_s_per_megagram);
    try std.testing.expectEqual(@as(f64, 12), state.apatite_mol_per_megagram);
    try std.testing.expectEqual(@as(f64, 16), state.calcium_sulfate_mol_per_megagram);
}

test "STARTE inactive soil-source guard leaves invalid dormant inputs untouched" {
    const empty = [_]f64{};
    var state: State = .{ .calcium_mol_per_megagram = 7 };
    const before = state;
    try std.testing.expect(!try initialize(.{ .source_index = 2, .day_index = 0, .layer_index = 99 }, .{ .column_count = 0, .row_count = 0, .layer_count_including_surface = 0 }, 99, 99, uniformInputs(&empty), &state));
    try std.testing.expectEqualDeep(before, state);
}

test "STARTE soil-source snapshot rejects late non-finite input atomically" {
    const valid = [_]f64{1};
    const invalid = [_]f64{std.math.nan(f64)};
    var inputs = uniformInputs(&valid);
    inputs.calcium_sulfate_mol_per_megagram_by_cell_layer = &invalid;
    var state: State = .{ .aluminum_mol_per_megagram = 2, .calcium_sulfate_mol_per_megagram = 3 };
    const before = state;
    try std.testing.expectError(error.NonFiniteSoilIonSourceInput, initialize(.{ .source_index = 3, .day_index = 1, .layer_index = 0 }, .{ .column_count = 1, .row_count = 1, .layer_count_including_surface = 1 }, 0, 0, inputs, &state));
    try std.testing.expectEqualDeep(before, state);
}
